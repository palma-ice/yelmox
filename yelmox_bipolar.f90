program yelmox

    ! == YelmoX-Bipolar, by Sergio Pérez-Montero == !
    ! Some notes:
    !   -> An additional type is built, yelmox. Contins everything that an individual yelmox domain needs to work
    !   -> The coding principle is to reproduce 1:1 the yelmox.f90 in those groups in order to be completely analogue to 1-domain yelmox executables
    !   -> This script is organized in "groups" in order to make the code more clear. 
    !      Each group is a subroutine that contains the specified lines of the original yelmox.f90 code
    !      
    !      Some exceptions are the groups that involve actions related with the coupling between domains and/or the ocean box model
    !      Groups description:
    !       1 - Parameters reading
    !       2 - Output definitions
    !       3 - Initialization:
    !           3.1 - Initialize ice sheet model
    !           3.2 - Initialize external models
    !           3.3 - Initialize boundary conditions  
    !           3.4 - Initialize output files
    !       4 - Loop:
    !           4.1 - Spinup procedure
    !           ... - Update ice sheet and external models (not grouped)
    !           4.2 - Update output

    ! Tools
    use nml 
    use ncio 
    use timer
    use timeout 
    use yelmo 
    use yelmo_tools, only : smooth_gauss_2D
    use ice_optimization
    use basal_dragging, only : calc_lambda_bed_lin, calc_lambda_bed_exp

    ! Yelmo libraries   (this could be removed in the near future! spm)
    use yelmo_ice
    use yelmo_grid
    use yelmo_timesteps
    use yelmo_topography
    use yelmo_dynamics
    use yelmo_material
    use yelmo_thermodynamics
    use yelmo_regions
    use yelmo_boundaries
    use yelmo_data 
    use mass_conservation

    ! External libraries
    use fastisostasy    ! also reexports barysealevel
    use snapclim
    use marine_shelf
    use smbpal
    use sediments
    use geothermal

    ! Ocean Box Model (OBM) libraries
    use obm_defs
    use stommel
    use nautilus
    use obm
    use ice2ocean
    use ocean2ice
    use atm2ocean

    implicit none

    ! Type definitions
    type reg_def_class 
        character(len=56)  :: name 
        character(len=512) :: fnm
        logical, allocatable :: mask(:,:) 
        logical :: write 
    end type

    type stats_class
        logical  :: defined
        real(wp) :: pd_rmse_H
        real(wp) :: pd_rmse_H_flt
        real(wp) :: V
        real(wp) :: A
        real(wp) :: A_grnd
        real(wp) :: A_flt
    end type

    type negis_params
        logical  :: use_negis_par
        real(wp) :: cf_0
        real(wp) :: cf_1
        real(wp) :: cf_centre
        real(wp) :: cf_north
        real(wp) :: cf_south

        real(wp) :: cf_x
        
    end type

    type yelmox_class   ! each yelmox domain contains all these objects

        ! Related with ´yelmo´
        type(yelmo_class)      :: yelmo1
        character(len=256)     :: domain

        ! External models
        type(snapclim_class)   :: snp1
        type(marshelf_class)   :: mshlf1
        type(smbpal_class)     :: smbpal1
        type(sediments_class)  :: sed1
        type(geothermal_class) :: gthrm1
        type(isos_class)       :: isos1

        ! Output files
        character(len=256) :: file1D, file2D, file2D_small, file_restart
        character(len=256) :: file_isos
        character(len=512) :: path_lgm

        character(len=512) :: path_par  ! domain specific nml
        character(len=512) :: path_hydro_mask
        real(wp), allocatable :: hydro_mask(:,:)

        ! Regions
        type(reg_def_class) :: reg1 
        type(reg_def_class) :: reg2 
        type(reg_def_class) :: reg3 
        character(len=512)    :: regions_mask_fnm
        real(wp), allocatable :: regions_mask(:,:) 

        ! Optimization
        type(ice_opt_params) :: opt  

        ! Statistics
        type(stats_class)    :: stats_pd_1, stats_lgm, stats_pd_2

        ! NEGIS
        type(negis_params)   :: ngs 

        ! Internal parameters
        logical  :: running_laurentide
        logical  :: laurentide_init_const_H 
        real(wp) :: laurentide_time_equil 

        logical  :: running_greenland
        logical  :: greenland_init_marine_H
        logical  :: scale_glacial_smb

    end type

    type ctrl_params
        character(len=56) :: run_step
        real(wp) :: time_init
        real(wp) :: time_end
        real(wp) :: time_equil      ! Only for spinup
        real(wp) :: time_const      ! Only for spinup
        real(wp) :: time_lgm_step 
        real(wp) :: time_pd_step
        real(wp) :: dtt
        real(wp) :: dt_restart
        real(wp) :: dt_clim

        logical  :: transient_clim
        logical  :: use_lgm_step
        logical  :: use_pd_step
        logical  :: with_ice_sheet 

        character(len=56) :: equil_method

        real(wp) :: isos_tau_1 
        real(wp) :: isos_tau_2 
        real(wp) :: isos_sigma 

        ! bipolar-specific control parameters
        logical :: active_north, active_south
        logical  :: obm2ism, ism2obm, atm2obm, use_yelmox, active_obm
        character(len=512) :: path_par_ocn, obm_name
        character(len=512) :: hydro_mask_north, hydro_mask_south
        
    end type 

    ! Run parameters/variables
    character(len=256) :: outfldr
    real(wp) :: time, time_bp, time_elapsed
    real(wp) :: dT_now 
    integer  :: n

    type(timeout_class) :: tm_1D, tm_2D, tm_2Dsm

    ! Model timing
    type(timer_class)  :: tmr
    type(timer_class)  :: tmrs
    character(len=512) :: tmr_file

    ! ctrl params
    type(ctrl_params)    :: ctl

    real(wp), parameter :: time_lgm = -19050.0_wp  ! [yr CE] == 21 kyr ago 
    real(wp), parameter :: time_pd  =   1950.0_wp  ! [yr CE] ==  0 kyr ago 

    ! OBM files
    character(len=256) :: obm_file, obm_file_restart

    ! Bipolar simulation objects
    type(yelmox_class) :: yelmox_north, yelmox_south  ! Northern and Southern hemisphere objects (each one is a yelmox domain)
    type(bsl_class)    :: bsl                   ! Barystatic sea level object is common to every domain
    character(len=256) :: file_bsl
    type(obm_class)    :: obm                   ! Ocean Box Model object
    character(len=512) :: bipolar_path_par      ! Generalistic namelist

    ! #### Start YelmoX-Bipolar simulation ##### !

    ! Group 1: Determine the general parameter file and read parameters
    call read_parameters(bipolar_path_par, ctl, yelmox_north, yelmox_south, tm_1D, tm_2D, tm_2Dsm)   ! yelmox.f90 lines 129-185

    ! Set initial time 
    time    = ctl%time_init 
    time_bp = time - 1950.0_wp 

    ! Start timing
    call timer_step(tmr,comp=-1) 
    
    ! Assume program is running from the output folder
    outfldr = "./"

    ! Group 2: Define input and output locations ! yelmox.f90 lines 196-200
    call define_input_and_output_locations(outfldr, "north", yelmox_north)
    call define_input_and_output_locations(outfldr, "south", yelmox_south)

    file_bsl = trim(outfldr)//"bsl.nc"
    tmr_file = trim(outfldr)//"timer_table.txt"

    ! Print summary of run settings 
    write(*,*)
    write(*,*) "transient_clim: ",  ctl%transient_clim
    write(*,*) "use_lgm_step:   ",  ctl%use_lgm_step
    write(*,*) "use_pd_step:    ",  ctl%use_pd_step
    write(*,*) "with_ice_sheet: ",  ctl%with_ice_sheet
    write(*,*) "equil_method:   ",  trim(ctl%equil_method)
    write(*,*)
    write(*,*) "time_init:  ",      ctl%time_init 
    write(*,*) "time_end:   ",      ctl%time_end 
    write(*,*) "dtt:        ",      ctl%dtt 
    write(*,*) "dt_restart: ",      ctl%dt_restart 
    write(*,*) 
    
    if (.not. ctl%transient_clim) then 
        write(*,*) "time_equil: ",    ctl%time_equil 
        write(*,*) "time_const: ",    ctl%time_const 

        if (ctl%use_lgm_step) then 
            write(*,*) "time_lgm_step: ", ctl%time_lgm_step 
        end if 

        if (ctl%use_pd_step) then 
            write(*,*) "time_pd_step: ", ctl%time_pd_step 
        end if 

        ! Set time before present == to constant climate time 
        time_bp = ctl%time_const - 1950.0_wp

    end if 

    write(*,*) "time    = ", time 
    write(*,*) "time_bp = ", time_bp 
    write(*,*) 

    ! === Initialize Yelmo X domains, Group 3
    ! Initialize barysealevel model (common to all the domains)
    call bsl_init(bsl, bipolar_path_par, time_bp)
    call bsl_update(bsl, year_bp=time_bp)
    call bsl_write_init(bsl, file_bsl, time)

    if (ctl%active_north) then
        call initialize_icesheet(bipolar_path_par, "north", time, ctl, yelmox_north)  ! yelmox.f90 lines 243-394
        
        ! BIPOLAR: allocate and load hydrographic mask
        allocate(yelmox_north%hydro_mask(yelmox_north%yelmo1%grd%nx,yelmox_north%yelmo1%grd%ny))
        call nc_read(yelmox_north%path_hydro_mask,"mask",yelmox_north%hydro_mask)
        
        call initialize_external_models(yelmox_north%path_par, yelmox_north)   ! yelmox.f90 lines 396-419
        call initialize_boundary_conditions(ctl, yelmox_north, bsl, time, time_bp)   ! yelmox.f90 lines 421-614

        call initialize_output_files(yelmox_north, time)

        write(*,*) "Northern Hemisphere domain initialized"
    end if

    if (ctl%active_south) then
        call initialize_icesheet(bipolar_path_par, "south", time, ctl, yelmox_south)  ! yelmox.f90 lines 243-394
        
        ! BIPOLAR: allocate and load hydrographic mask
        allocate(yelmox_south%hydro_mask(yelmox_south%yelmo1%grd%nx,yelmox_south%yelmo1%grd%ny))
        call nc_read(yelmox_south%path_hydro_mask,"mask",yelmox_south%hydro_mask)
        
        call initialize_external_models(yelmox_south%path_par, yelmox_south)   ! yelmox.f90 lines 396-419
        call initialize_boundary_conditions(ctl, yelmox_south, bsl, time, time_bp)   ! yelmox.f90 lines 421-614

        call initialize_output_files(yelmox_south, time)

        write(*,*) "Southern Hemisphere domain initialized"
    end if

    ! === Initialize ocean box model =====
    if (ctl%active_obm) then
        obm_file = trim(outfldr)//trim(ctl%obm_name)//".nc"
        obm_file_restart = trim(outfldr)//trim(ctl%obm_name)//"_restart.nc"

        call obm_init(obm, ctl%path_par_ocn, ctl%obm_name)
        call write_obm_init(obm_file, time, "years")
        call write_obm_update(obm, obm_file, ctl%obm_name, time)

        if (ctl%active_north) then
            if (ctl%atm2obm) then
                call update_theta_from_snapclim(obm%thetan, yelmox_north%snp1%now%tsl_ann, yelmox_north%hydro_mask, "north", ctl%obm_name)
            end if
        end if

        if (ctl%active_south) then
            if (ctl%atm2obm) then
                call update_theta_from_snapclim(obm%thetas, yelmox_south%snp1%now%tsl_ann, yelmox_south%hydro_mask, "south", ctl%obm_name)
            end if
        end if

    end if
    
    call timer_step(tmr,comp=1,label="initialization") 
    call timer_step(tmrs,comp=-1)

    ! ==== Begin main time loop ===== Group 4

    ! Advance timesteps
    do n = 0, ceiling((ctl%time_end-ctl%time_init)/ctl%dtt)

        ! Get current time 
        time = ctl%time_init + n*ctl%dtt

        if (ctl%transient_clim) then 
            time_bp = time - 1950.0_wp
        else

            if (ctl%use_lgm_step .and. time .ge. ctl%time_lgm_step) then 
                ctl%time_const = time_lgm 
            end if 

            if (ctl%use_pd_step .and. time .ge. ctl%time_pd_step) then 
                ctl%time_const = time_pd 
            end if 

            time_bp = ctl%time_const - 1950.0_wp 

        end if

        time_elapsed = time - ctl%time_init 

        ! Spin-up procedure - only relevant for time-time_init <= time_equil
        if (ctl%active_north) then
            call spinup_procedure(ctl, yelmox_north, time, time_elapsed)
        end if
        if (ctl%active_south) then
            call spinup_procedure(ctl, yelmox_south, time, time_elapsed)
        end if

        call timer_step(tmrs,comp=0) 

        ! == ISOSTASY and SEA LEVEL ======================================================
        call bsl_update(bsl, time_bp)
        if (ctl%active_north) then
            call isos_update(yelmox_north%isos1, yelmox_north%yelmo1%tpo%now%H_ice, time, bsl, dwdt_corr=yelmox_north%yelmo1%bnd%dzbdt_corr)
            yelmox_north%yelmo1%bnd%z_bed = yelmox_north%isos1%out%z_bed
            yelmox_north%yelmo1%bnd%z_sl  = yelmox_north%isos1%out%z_ss
        end if
        if (ctl%active_south) then
            call isos_update(yelmox_south%isos1, yelmox_south%yelmo1%tpo%now%H_ice, time, bsl, dwdt_corr=yelmox_south%yelmo1%bnd%dzbdt_corr)
            yelmox_south%yelmo1%bnd%z_bed = yelmox_south%isos1%out%z_bed
            yelmox_south%yelmo1%bnd%z_sl  = yelmox_south%isos1%out%z_ss
        end if

        call timer_step(tmrs,comp=1,time_mod=[time-ctl%dtt,time]*1e-3,label="isostasy") 

        ! == Ocean Box Model ===================================================
        if (ctl%active_obm) then
            if (ctl%ism2obm) then   ! apply fwf forcing to obm via ism domains
                if (ctl%active_north) then
                    obm%fn = calc_fwf(yelmox_north%yelmo1%bnd%c%rho_w,yelmox_north%yelmo1%bnd%c%rho_ice, &
                                        yelmox_north%yelmo1%bnd%c%sec_year, &
                                        yelmox_north%yelmo1%tpo%now%H_ice,yelmox_north%yelmo1%tpo%now%dHidt, &
                                        yelmox_north%yelmo1%tpo%now%f_grnd, &
                                        yelmox_north%yelmo1%tpo%par%dx,yelmox_north%yelmo1%tpo%par%dy, yelmox_north%hydro_mask, "north")
                else
                    obm%fn = 0.0
                end if
                if (ctl%active_south) then
                    obm%fs = calc_fwf(yelmox_south%yelmo1%bnd%c%rho_w,yelmox_south%yelmo1%bnd%c%rho_ice, &
                                        yelmox_south%yelmo1%bnd%c%sec_year, &
                                        yelmox_south%yelmo1%tpo%now%H_ice,yelmox_south%yelmo1%tpo%now%dHidt, &
                                        yelmox_south%yelmo1%tpo%now%f_grnd, &
                                        yelmox_south%yelmo1%tpo%par%dx,yelmox_south%yelmo1%tpo%par%dy, yelmox_south%hydro_mask, "south")
                else
                    obm%fs = 0.0
                end if
            else
                obm%fn = 0.0
                obm%fs = 0.0
            end if

            call obm_update(obm, ctl%dtt, ctl%obm_name)
        end if ! use_obm

        ! == ICE SHEET ===================================================
        if (ctl%active_north) then
            if (yelmox_north%running_greenland .and. yelmox_north%ngs%use_negis_par) then 
                ! Update cb_ref using negis parameters 
                call negis_update_cb_ref(yelmox_north%yelmo1,yelmox_north%ngs,time)
            end if
            
            ! Update Yelmo
            if (ctl%with_ice_sheet .and. (.not. (n .eq. 0 .and. yelmox_north%yelmo1%par%use_restart)) ) then
                call yelmo_update(yelmox_north%yelmo1,time)
            end if
        end if
        if (ctl%active_south) then
            if (yelmox_south%running_greenland .and. yelmox_south%ngs%use_negis_par) then 
                ! Update cb_ref using negis parameters 
                call negis_update_cb_ref(yelmox_south%yelmo1,yelmox_south%ngs,time)
            end if
            
            ! Update Yelmo
            if (ctl%with_ice_sheet .and. (.not. (n .eq. 0 .and. yelmox_south%yelmo1%par%use_restart)) ) then
                call yelmo_update(yelmox_south%yelmo1,time)
            end if
        end if

        call timer_step(tmrs,comp=2,time_mod=[time-ctl%dtt,time]*1e-3,label="yelmo")

        ! == CLIMATE (ATMOSPHERE AND OCEAN) ====================================
        if (mod(nint(time*100),nint(ctl%dt_clim*100))==0) then
                            ! Update snapclim

            if (ctl%active_north) then
                call snapclim_update(yelmox_north%snp1,z_srf=yelmox_north%yelmo1%tpo%now%z_srf,time=time_bp,domain=yelmox_north%domain,dx=yelmox_north%yelmo1%grd%dx,basins=yelmox_north%yelmo1%bnd%basins) 
            end if
            if (ctl%active_south) then
                call snapclim_update(yelmox_south%snp1,z_srf=yelmox_south%yelmo1%tpo%now%z_srf,time=time_bp,domain=yelmox_south%domain,dx=yelmox_south%yelmo1%grd%dx,basins=yelmox_south%yelmo1%bnd%basins) 
            end if

            ! atm2obm
            if (ctl%atm2obm) then
                if (ctl%active_north) then
                    call update_theta_from_snapclim(obm%thetan, yelmox_north%snp1%now%tsl_ann, yelmox_north%hydro_mask, "north", ctl%obm_name)
                end if
                if (ctl%active_south) then
                    call update_theta_from_snapclim(obm%thetas, yelmox_south%snp1%now%tsl_ann, yelmox_south%hydro_mask, "south", ctl%obm_name)
                end if
            end if

            ! obm2ism
            if (ctl%obm2ism) then
                if (ctl%active_north) then
                    call calc_ocean_temperature_field(yelmox_north%snp1%now%to_ann, obm%tn, ctl%obm_name)
                    call calc_ocean_salinity_field(yelmox_north%snp1%now%so_ann, obm%sn)
                end if
                if (ctl%active_south) then
                    call calc_ocean_temperature_field(yelmox_south%snp1%now%to_ann, obm%ts, ctl%obm_name)
                    call calc_ocean_salinity_field(yelmox_south%snp1%now%so_ann, obm%ss)
                end if
            end if

        end if

        ! == SURFACE MASS BALANCE ==============================================
        ! ajr: just testing...
        dT_now = 0.0 
        !if (time .gt. 7000.0) dT_now = 4.0 

        if (ctl%active_north) then
            call smbpal_update_monthly(yelmox_north%smbpal1,yelmox_north%snp1%now%tas,yelmox_north%snp1%now%pr, &
                                       yelmox_north%yelmo1%tpo%now%z_srf,yelmox_north%yelmo1%tpo%now%H_ice,time_bp) 
            yelmox_north%yelmo1%bnd%smb   = yelmox_north%smbpal1%ann%smb*yelmox_north%yelmo1%bnd%c%conv_we_ie*1e-3       ! [mm we/a] => [m ie/a]
            yelmox_north%yelmo1%bnd%T_srf = yelmox_north%smbpal1%ann%tsrf 

            if (trim(yelmox_north%yelmo1%par%domain) .eq. "Greenland" .and. yelmox_north%scale_glacial_smb) then 
                ! Modify glacial smb
                call calc_glacial_smb(yelmox_north%yelmo1%bnd%smb,yelmox_north%yelmo1%grd%lat,yelmox_north%snp1%now%ta_ann,yelmox_north%snp1%clim0%ta_ann)
            end if
        
            ! yelmox_north%yelmo1%bnd%smb   = yelmox_north%yelmo1%dta%pd%smb
            ! yelmox_north%yelmo1%bnd%T_srf = yelmox_north%yelmo1%dta%pd%t2m
            
            if (yelmox_north%running_laurentide .and. yelmox_north%laurentide_init_const_H  &
                    .and. (time-ctl%time_init) .lt. yelmox_north%laurentide_time_equil ) then 
                ! Additionally ensure smb is postive for land above 50degN in Laurentide region
                ! to make sure ice grows everywhere needed (Coridilleran ice sheet mainly)
                where (yelmox_north%yelmo1%bnd%regions .eq. 1.1 .and. yelmox_north%yelmo1%grd%lat .gt. 50.0 .and. &
                       yelmox_north%yelmo1%bnd%z_bed .gt. 0.0 .and. yelmox_north%yelmo1%bnd%smb .lt. 0.0 ) yelmox_north%yelmo1%bnd%smb = 0.5 
            end if
        end if
        if (ctl%active_south) then
            call smbpal_update_monthly(yelmox_south%smbpal1,yelmox_south%snp1%now%tas,yelmox_south%snp1%now%pr, &
                                       yelmox_south%yelmo1%tpo%now%z_srf,yelmox_south%yelmo1%tpo%now%H_ice,time_bp) 
            yelmox_south%yelmo1%bnd%smb   = yelmox_south%smbpal1%ann%smb*yelmox_south%yelmo1%bnd%c%conv_we_ie*1e-3       ! [mm we/a] => [m ie/a]
            yelmox_south%yelmo1%bnd%T_srf = yelmox_south%smbpal1%ann%tsrf 

            if (trim(yelmox_south%yelmo1%par%domain) .eq. "Greenland" .and. yelmox_south%scale_glacial_smb) then 
                ! Modify glacial smb
                call calc_glacial_smb(yelmox_south%yelmo1%bnd%smb,yelmox_south%yelmo1%grd%lat,yelmox_south%snp1%now%ta_ann,yelmox_south%snp1%clim0%ta_ann)
            end if
        
            ! yelmox_south%yelmo1%bnd%smb   = yelmox_south%yelmo1%dta%pd%smb
            ! yelmox_south%yelmo1%bnd%T_srf = yelmox_south%yelmo1%dta%pd%t2m
            
            if (yelmox_south%running_laurentide .and. yelmox_south%laurentide_init_const_H  &
                    .and. (time-ctl%time_init) .lt. yelmox_south%laurentide_time_equil ) then 
                ! Additionally ensure smb is postive for land above 50degN in Laurentide region
                ! to make sure ice grows everywhere needed (Coridilleran ice sheet mainly)
                where (yelmox_south%yelmo1%bnd%regions .eq. 1.1 .and. yelmox_south%yelmo1%grd%lat .gt. 50.0 .and. &
                       yelmox_south%yelmo1%bnd%z_bed .gt. 0.0 .and. yelmox_south%yelmo1%bnd%smb .lt. 0.0 ) yelmox_south%yelmo1%bnd%smb = 0.5 
            end if
        end if

        ! == MARINE AND TOTAL BASAL MASS BALANCE ===============================
        if (ctl%active_north) then
            call marshelf_update_shelf(yelmox_north%mshlf1,yelmox_north%yelmo1%tpo%now%H_ice,yelmox_north%yelmo1%bnd%z_bed,yelmox_north%yelmo1%tpo%now%f_grnd, &
                        yelmox_north%yelmo1%bnd%basins,yelmox_north%yelmo1%bnd%z_sl,yelmox_north%yelmo1%grd%dx,yelmox_north%snp1%now%depth, &
                        yelmox_north%snp1%now%to_ann,yelmox_north%snp1%now%so_ann,dto_ann=yelmox_north%snp1%now%to_ann-yelmox_north%snp1%clim0%to_ann)

            call marshelf_update(yelmox_north%mshlf1,yelmox_north%yelmo1%tpo%now%H_ice,yelmox_north%yelmo1%bnd%z_bed,yelmox_north%yelmo1%tpo%now%f_grnd, &
                                 yelmox_north%yelmo1%bnd%regions,yelmox_north%yelmo1%bnd%basins,yelmox_north%yelmo1%bnd%z_sl,dx=yelmox_north%yelmo1%grd%dx)

            yelmox_north%yelmo1%bnd%bmb_shlf = yelmox_north%mshlf1%now%bmb_shlf
            yelmox_north%yelmo1%bnd%T_shlf   = yelmox_north%mshlf1%now%T_shlf
        end if
        if (ctl%active_south) then
            call marshelf_update_shelf(yelmox_south%mshlf1,yelmox_south%yelmo1%tpo%now%H_ice,yelmox_south%yelmo1%bnd%z_bed,yelmox_south%yelmo1%tpo%now%f_grnd, &
                        yelmox_south%yelmo1%bnd%basins,yelmox_south%yelmo1%bnd%z_sl,yelmox_south%yelmo1%grd%dx,yelmox_south%snp1%now%depth, &
                        yelmox_south%snp1%now%to_ann,yelmox_south%snp1%now%so_ann,dto_ann=yelmox_south%snp1%now%to_ann-yelmox_south%snp1%clim0%to_ann)

            call marshelf_update(yelmox_south%mshlf1,yelmox_south%yelmo1%tpo%now%H_ice,yelmox_south%yelmo1%bnd%z_bed,yelmox_south%yelmo1%tpo%now%f_grnd, &
                                 yelmox_south%yelmo1%bnd%regions,yelmox_south%yelmo1%bnd%basins,yelmox_south%yelmo1%bnd%z_sl,dx=yelmox_south%yelmo1%grd%dx)

            yelmox_south%yelmo1%bnd%bmb_shlf = yelmox_south%mshlf1%now%bmb_shlf
            yelmox_south%yelmo1%bnd%T_shlf   = yelmox_south%mshlf1%now%T_shlf
        end if

        call timer_step(tmrs,comp=3,time_mod=[time-ctl%dtt,time]*1e-3,label="climate") 

        ! == MODEL OUTPUT =======================================================
        if (ctl%active_north) then
            call update_output(yelmox_north, bsl, time, tm_1D, tm_2D, tm_2Dsm,"north")
        end if
        if (ctl%active_south) then
            call update_output(yelmox_south, bsl, time, tm_1D, tm_2D, tm_2Dsm,"south")
        end if

        if (ctl%active_obm) then
            if (timeout_check(tm_1D,time)) then
                call write_obm_update(obm, obm_file, ctl%obm_name, time)
            end if

            if (mod(nint(time*100),nint(ctl%dt_restart*100))==0) then
                call write_obm_restart(obm, obm_file_restart, time, "years")
            end if
        end if

        call timer_step(tmrs,comp=4,time_mod=[time-ctl%dtt,time]*1e-3,label="io") 
        
        if (mod(time_elapsed,10.0)==0) then
            ! Print timestep timing info and write log table
            call timer_write_table(tmrs,[time,ctl%dtt]*1e-3,"m",tmr_file,init=time_elapsed .eq. 0.0)
        end if 

        if (mod(time,10.0)==0 .and. (.not. yelmo_log)) then
            write(*,"(a,f14.4)") "yelmo:: time = ", time
        end if 
        
    end do 

    ! Stop timing
    call timer_step(tmr,comp=2,time_mod=[ctl%time_init,time]*1e-3,label="timeloop") 

    ! Write the restart snapshot for the end of the simulation
    if (ctl%active_north) then
        call yelmox_bipolar_restart_write(bsl,yelmox_north%isos1,yelmox_north%yelmo1,time,hemisphere="north")
    end if
    if (ctl%active_south) then
        call yelmox_bipolar_restart_write(bsl,yelmox_south%isos1,yelmox_south%yelmo1,time,hemisphere="south")
    end if
    if (ctl%active_obm) then
        call write_obm_restart(obm, obm_file_restart, time, "years")
    end if

    ! Finalize program
    if (ctl%active_north) then
        call yelmo_end(yelmox_north%yelmo1,time=time)
    end if
    if (ctl%active_south) then
        call yelmo_end(yelmox_south%yelmo1,time=time)
    end if

    ! Print timing summary
    call timer_print_summary(tmr,units="m",units_mod="kyr",time_mod=time*1e-3)

contains

    ! Yelmo X original subroutines
    subroutine calc_glacial_smb(smb,lat2D,ta_ann,ta_ann_pd)

        implicit none

        real(wp), intent(INOUT) :: smb(:,:)
        real(wp), intent(IN)    :: lat2D(:,:)
        real(wp), intent(IN)    :: ta_ann(:,:)
        real(wp), intent(IN)    :: ta_ann_pd(:,:)

        ! Local variables 
        integer  :: i, j, nx, ny 
        real(wp) :: t0, tnow
        real(wp) :: at
        real(wp) :: fac

        real(wp), parameter :: dt_lgm  = -8.0 
        real(wp), parameter :: lat_lim = 55.0 
        real(wp), parameter :: fac_lim = 0.90

        nx = size(smb,1)
        ny = size(smb,2) 

        ! Determine a quasi glacial-interglacial index
        ! 0: interglacial
        ! 1: glacial 
        tnow = sum(ta_ann) / real(nx*ny,wp)
        t0   = sum(ta_ann_pd) / real(nx*ny,wp)

        at = (tnow-t0)/dt_lgm
        if (at .lt. 0.0) at = 0.0
        if (at .gt. 1.0) at = 1.0
        
        ! Now determine smb scaling as a function of glacial index and latitude
        ! fac==0: no change to smb
        ! 0<fac<=1: scale smb up to value of fac_lim

        do j = 1, ny 
        do i = 1, nx

            if (smb(i,j) .lt. 0.0 .and. lat2D(i,j) .gt. lat_lim) then
                ! Calculate scaling here

                smb(i,j) = smb(i,j) -  smb(i,j) * at * fac_lim

            end if

        end do
        end do

        return

    end subroutine calc_glacial_smb

    subroutine negis_update_cb_ref(ylmo,ngs,t)

        implicit none

        type(yelmo_class),  intent(INOUT) :: ylmo
        type(negis_params), intent(INOUT) :: ngs 
        real(wp), intent(IN) :: t 

        ! Local variables
        integer :: i, j, nx, ny 

        nx = ylmo%grd%nx 
        ny = ylmo%grd%ny 

        if (t .lt. -11e3) then 
            ngs%cf_x = ngs%cf_0
        else
            ! Linear interpolation from cf_0 to cf_1  
            ngs%cf_x = ngs%cf_0 + (t - (-11e3))/ ((0.0) - (-11e3)) * (ngs%cf_1 - ngs%cf_0)
        end if

        if (t .lt. -4e3) then
            ngs%cf_south = 1.0
        else
            ngs%cf_north = 1.0
        end if

        ! Update bed roughness coefficients cb_ref and c_bed (which are independent of velocity)
        ! like normal, using the default function defined in Yelmo:
        call calc_cb_ref(ylmo%dyn%now%cb_ref,ylmo%bnd%z_bed,ylmo%bnd%z_bed_sd,ylmo%bnd%z_sl, &
                            ylmo%bnd%H_sed,ylmo%dyn%par%till_f_sed,ylmo%dyn%par%till_sed_min, &
                            ylmo%dyn%par%till_sed_max,ylmo%dyn%par%till_cf_ref,ylmo%dyn%par%till_cf_min, &
                            ylmo%dyn%par%till_z0,ylmo%dyn%par%till_z1,ylmo%dyn%par%till_n_sd, &
                            ylmo%dyn%par%till_scale,ylmo%dyn%par%till_method)

        ! === Finally, apply NEGIS scaling =============================

        do j = 1, ny 
        do i = 1, nx 

            if(ylmo%bnd%basins(i,j) .eq. 9.1) ylmo%dyn%now%cb_ref(i,j) = ylmo%dyn%now%cb_ref(i,j) * ngs%cf_centre
            if(ylmo%bnd%basins(i,j) .eq. 9.2) ylmo%dyn%now%cb_ref(i,j) = ylmo%dyn%now%cb_ref(i,j) * ngs%cf_south
            if(ylmo%bnd%basins(i,j) .eq. 9.3) ylmo%dyn%now%cb_ref(i,j) = ylmo%dyn%now%cb_ref(i,j) * ngs%cf_north

        end do
        end do

        return

    end subroutine negis_update_cb_ref

    subroutine yelmox_write_step(ylmo,snp,mshlf,srf,filename,time)

        implicit none 
        
        type(yelmo_class),      intent(IN) :: ylmo
        type(snapclim_class),   intent(IN) :: snp 
        type(marshelf_class),   intent(IN) :: mshlf 
        type(smbpal_class),     intent(IN) :: srf 
        !type(sediments_class),  intent(IN) :: sed 
        !type(geothermal_class), intent(IN) :: gthrm
        !type(isos_class),       intent(IN) :: isos
        
        character(len=*),  intent(IN) :: filename
        real(wp), intent(IN) :: time

        ! Local variables
        integer  :: ncid, n
        real(wp) :: time_prev 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_time_index(filename,"time",time,ncid)

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Write model metrics (model speed, dt, eta)
        call yelmo_write_step_model_metrics(filename,ylmo,n,ncid)

        ! Write present-day data metrics (rmse[H],etc)
        call yelmo_write_step_pd_metrics(filename,ylmo,n,ncid)
                
        ! == yelmo_topography ==
        call yelmo_write_var(filename,"H_ice",ylmo,n,ncid)
        call yelmo_write_var(filename,"z_srf",ylmo,n,ncid)
        call yelmo_write_var(filename,"mask_bed",ylmo,n,ncid)
        call yelmo_write_var(filename,"mb_net",ylmo,n,ncid)
        call yelmo_write_var(filename,"mb_resid",ylmo,n,ncid)
        call yelmo_write_var(filename,"smb",ylmo,n,ncid)
        call yelmo_write_var(filename,"bmb",ylmo,n,ncid)
        call yelmo_write_var(filename,"fmb",ylmo,n,ncid)
        call yelmo_write_var(filename,"cmb",ylmo,n,ncid)
        call yelmo_write_var(filename,"H_grnd",ylmo,n,ncid)
        call yelmo_write_var(filename,"N_eff",ylmo,n,ncid)
        call yelmo_write_var(filename,"f_grnd",ylmo,n,ncid)
        call yelmo_write_var(filename,"f_ice",ylmo,n,ncid)
        call yelmo_write_var(filename,"dHidt",ylmo,n,ncid)
        
        call yelmo_write_var(filename,"cmb_flt",ylmo,n,ncid)
        call yelmo_write_var(filename,"cmb_grnd",ylmo,n,ncid)

        ! call yelmo_write_var(filename,"eps_eff",ylmo,n,ncid)
        ! call yelmo_write_var(filename,"tau_eff",ylmo,n,ncid)

        ! == yelmo_dynamics ==
        call yelmo_write_var(filename,"cb_ref",ylmo,n,ncid)
        call yelmo_write_var(filename,"c_bed",ylmo,n,ncid)
        call yelmo_write_var(filename,"beta",ylmo,n,ncid)
        call yelmo_write_var(filename,"visc_eff_int",ylmo,n,ncid)
        call yelmo_write_var(filename,"taud",ylmo,n,ncid)
        call yelmo_write_var(filename,"taub",ylmo,n,ncid)
        call yelmo_write_var(filename,"uxy_b",ylmo,n,ncid)
        call yelmo_write_var(filename,"uxy_bar",ylmo,n,ncid)
        call yelmo_write_var(filename,"uxy_s",ylmo,n,ncid)
        
        ! == yelmo_material ==
        call yelmo_write_var(filename,"enh_bar",ylmo,n,ncid)
        !call yelmo_write_var(filename,"ATT",ylmo,n,ncid)
        call yelmo_write_var(filename,"visc_int",ylmo,n,ncid)
        !call yelmo_write_var(filename,"strn_de",ylmo,n,ncid)
        !call yelmo_write_var(filename,"strn_te",ylmo,n,ncid)
        ! call yelmo_write_var(filename,"strn2D_de",ylmo,n,ncid)
        ! call yelmo_write_var(filename,"strn2D_div",ylmo,n,ncid)
        ! call yelmo_write_var(filename,"strs2D_te",ylmo,n,ncid)

        ! == yelmo_thermodynamics ==
        call yelmo_write_var(filename,"T_prime",ylmo,n,ncid)
        call yelmo_write_var(filename,"f_pmp",ylmo,n,ncid)
        call yelmo_write_var(filename,"Q_b",ylmo,n,ncid)
        call yelmo_write_var(filename,"bmb_grnd",ylmo,n,ncid)
        call yelmo_write_var(filename,"H_w",ylmo,n,ncid)
        
        ! == yelmo_boundaries ==
        call yelmo_write_var(filename,"z_bed",ylmo,n,ncid)
        call yelmo_write_var(filename,"z_sl",ylmo,n,ncid)
        call yelmo_write_var(filename,"smb_ref",ylmo,n,ncid)
        call yelmo_write_var(filename,"T_srf",ylmo,n,ncid)
        call yelmo_write_var(filename,"bmb_shlf",ylmo,n,ncid)
        call yelmo_write_var(filename,"Q_geo",ylmo,n,ncid)
        
        ! == yelmo_data (comparison with present-day) ==
        call yelmo_write_var(filename,"pd_err_H_ice",ylmo,n,ncid)
        call yelmo_write_var(filename,"pd_err_z_srf",ylmo,n,ncid)
        call yelmo_write_var(filename,"pd_err_uxy_s",ylmo,n,ncid)
        !call yelmo_write_var(filename,"pd_err_smb",ylmo,n,ncid)
        

        ! == yelmo extra fields ==

        call yelmo_write_var(filename,"ssa_mask_acx",ylmo,n,ncid)
        call yelmo_write_var(filename,"ssa_mask_acy",ylmo,n,ncid)
        call yelmo_write_var(filename,"dzsdx",ylmo,n,ncid)
        call yelmo_write_var(filename,"dzsdy",ylmo,n,ncid)
        call yelmo_write_var(filename,"f_grnd_acx",ylmo,n,ncid)
        call yelmo_write_var(filename,"f_grnd_acy",ylmo,n,ncid)
        call yelmo_write_var(filename,"taub_acx",ylmo,n,ncid)
        call yelmo_write_var(filename,"taub_acy",ylmo,n,ncid)
        call yelmo_write_var(filename,"taud_acx",ylmo,n,ncid)
        call yelmo_write_var(filename,"taud_acy",ylmo,n,ncid)
        call yelmo_write_var(filename,"ux_s",ylmo,n,ncid)
        call yelmo_write_var(filename,"uy_s",ylmo,n,ncid)
        call yelmo_write_var(filename,"ux_b",ylmo,n,ncid)
        call yelmo_write_var(filename,"uy_b",ylmo,n,ncid)
        ! call yelmo_write_var(filename,"ux_i_bar",ylmo,n,ncid)
        ! call yelmo_write_var(filename,"uy_i_bar",ylmo,n,ncid)
        call yelmo_write_var(filename,"ux_bar",ylmo,n,ncid)
        call yelmo_write_var(filename,"uy_bar",ylmo,n,ncid)
        call yelmo_write_var(filename,"beta_acx",ylmo,n,ncid)
        call yelmo_write_var(filename,"beta_acy",ylmo,n,ncid)

        call nc_write(filename,"err_pd_smb_ref",ylmo%bnd%smb-ylmo%dta%pd%smb,units="m/a ice equiv.",long_name="Surface mass balance error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"Q_strn_alt_units",ylmo%thrm%now%Q_strn/(ylmo%bnd%c%rho_ice*ylmo%thrm%now%cp),units="K a-1",long_name="Strain heating", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"Q_ice_b",ylmo%thrm%now%Q_ice_b,units="mW m-2",long_name="Basal ice heat flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="mW m-2",long_name="Basal frictional heating", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Static fields
        if (n .le. 1) then 
            call nc_write(filename,"H_sed",ylmo%bnd%H_sed,units="m",long_name="Sediment thickness", &
                        dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
        end if

        ! == snapclim ==

        call nc_write(filename,"Ta_ann",snp%now%ta_ann,units="K",long_name="Near-surface air temperature (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Ta_sum",snp%now%ta_sum,units="K",long_name="Near-surface air temperature (sum)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Pr_ann",snp%now%pr_ann*1e-3,units="m/a water equiv.",long_name="Precipitation (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        !call nc_write(filename,"pr",snp%now%pr*1e-3,units="m/a water equiv.",long_name="Precipitation (ann)", &
        !              dim1="xc",dim2="yc",dim3="month",dim4="time",start=[1,1,1,n],ncid=ncid)
              
        call nc_write(filename,"dTa_ann",snp%now%ta_ann-snp%clim0%ta_ann,units="K",long_name="Near-surface air temperature anomaly (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dTa_sum",snp%now%ta_sum-snp%clim0%ta_sum,units="K",long_name="Near-surface air temperature anomaly (sum)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dPr_ann",(snp%now%pr_ann-snp%clim0%pr_ann)*1e-3,units="m/a water equiv.",long_name="Precipitation anomaly (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! == smbpal ==

!         call nc_write(filename,"PDDs",srf%ann%PDDs,units="degC days",long_name="Positive degree days (annual total)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == marine_shelf ==

        call nc_write(filename,"T_shlf",mshlf%now%T_shlf,units="K",long_name="Shelf temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!        call nc_write(filename,"S_shlf",mshlf%now%S_shlf,units="PSU",long_name="Shelf salinity", &
!                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dT_shlf",mshlf%now%dT_shlf,units="K",long_name="Shelf temperature anomaly", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_ocn",mshlf%now%mask_ocn,units="", &
                     long_name="Ocean mask (0: land, 1: grline, 2: fltline, 3: open ocean, 4: deep ocean, 5: lakes)", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"tf_shlf",mshlf%now%tf_shlf,units="K",long_name="Shelf thermal forcing", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        if (trim(mshlf%par%bmb_method) .eq. "pico") then 
            call nc_write(filename,"d_shlf",mshlf%pico%now%d_shlf,units="km",long_name="Shelf distance to grounding line", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"d_if",mshlf%pico%now%d_if,units="km",long_name="Shelf distance to ice front", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"boxes",mshlf%pico%now%boxes,units="",long_name="Shelf boxes", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)      
            call nc_write(filename,"r_shlf",mshlf%pico%now%r_shlf,units="",long_name="Ratio of ice shelf", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"T_box",mshlf%pico%now%T_box,units="K?",long_name="Temperature of boxes", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"S_box",mshlf%pico%now%S_box,units="PSU",long_name="Salinity of boxes", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"A_box",mshlf%pico%now%A_box*1e-6,units="km2",long_name="Box area of ice shelf", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        else 
            call nc_write(filename,"tf_basin",mshlf%now%tf_basin,units="K",long_name="Mean basin thermal forcing", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"tf_shlf",mshlf%now%tf_shlf,units="K",long_name="Shelf thermal forcing", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"tf_corr",mshlf%now%tf_corr,units="K",long_name="Shelf thermal forcing applied correction factor", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

            call nc_write(filename,"tf_corr_basin",mshlf%now%tf_corr_basin,units="K",long_name="Shelf thermal forcing basin-wide correction factor", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"slope_base",mshlf%now%slope_base,units="",long_name="Shelf-base slope", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

            !call nc_write(filename,"slope_base",mshlf%now%slope_base,units="",long_name="Shelf-base slope", &
            !              dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        end if 

        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine yelmox_write_step

    ! Bipolar specific subroutines
    subroutine rewrite_common_parameters(filename, dom)
        
        implicit none

        character(len=512), intent(IN)    :: filename
        type(yelmo_class), intent(INOUT) :: dom

        ! == topography ==

        call ytopo_par_load(dom%tpo%par,filename,dom%grd%nx,dom%grd%ny,dom%grd%dx,init=.TRUE.)

        !call ytopo_alloc(dom%tpo%now,dom%tpo%par%nx,dom%tpo%par%ny)
        
        write(*,*) "rewrite_common_parameters:: topography rewritten."
        
        ! == dynamics == 

        call ydyn_par_load(dom%dyn%par,filename,dom%par%zeta_aa,dom%par%zeta_ac,dom%grd%nx,dom%grd%ny,dom%grd%dx,init=.TRUE.)

        !call ydyn_alloc(dom%dyn%now,dom%dyn%par%nx,dom%dyn%par%ny,dom%dyn%par%nz_aa,dom%dyn%par%nz_ac)
        
        write(*,*) "rewrite_common_parameters:: dynamics rewritten."
        
        ! == material == 

        call ymat_par_load(dom%mat%par,filename,dom%par%zeta_aa,dom%par%zeta_ac,dom%grd%nx,dom%grd%ny,dom%grd%dx,init=.TRUE.)

        !call ymat_alloc(dom%mat%now,dom%mat%par%nx,dom%mat%par%ny,dom%mat%par%nz_aa,dom%mat%par%nz_ac,dom%mat%par%n_iso)
        
        write(*,*) "rewrite_common_parameters:: material rewritten."
        
        ! == thermodynamics == 
        
        call ytherm_par_load(dom%thrm%par,filename,dom%par%zeta_aa,dom%par%zeta_ac,dom%grd%nx,dom%grd%ny,dom%grd%dx,init=.TRUE.)

        !call ytherm_alloc(dom%thrm%now,dom%thrm%par%nx,dom%thrm%par%ny,dom%thrm%par%nz_aa,dom%thrm%par%nz_ac,dom%thrm%par%nzr_aa)
        
        write(*,*) "rewrite_common_parameters:: thermodynamics rewritten."

        return

    end subroutine rewrite_common_parameters

    subroutine yelmox_bipolar_restart_write(bsl,isos,ylmo,time,fldr,hemisphere)

        implicit none

        type(bsl_class),    intent(IN) :: bsl
        type(isos_class),   intent(IN) :: isos
        type(yelmo_class),  intent(IN) :: ylmo
        real(wp),           intent(IN) :: time 
        character(len=*),   intent(IN), optional :: fldr
        character(len=5),   intent(IN), optional :: hemisphere

        ! Local variables
        real(wp) :: time_kyr
        character(len=32)   :: time_str
        character(len=1024) :: outfldr

        character(len=56) :: file_bsl   = "bsl_restart.nc"
        character(len=56) :: file_isos  = "isos_restart.nc"
        character(len=56) :: file_yelmo = "yelmo_restart.nc"

        file_bsl   = "bsl_"//trim(hemisphere)//"_restart.nc"
        file_isos  = "isos_"//trim(hemisphere)//"_restart.nc"
        file_yelmo = "yelmo_"//trim(hemisphere)//"_restart.nc"
        
        if (present(fldr)) then
            outfldr = trim(fldr)
        else
            time_kyr = time*1e-3
            write(time_str,"(f20.3)") time_kyr
            outfldr = "./"//"restart-"//trim(adjustl(time_str))//"-kyr"
        end if

        write(*,*) "yelmox_restart_write:: outfldr = ", trim(outfldr)

        ! Make directory (use -p to ignore if directory already exists)
        call execute_command_line('mkdir -p "' // trim(outfldr) // '"')
        
        call bsl_restart_write(bsl,trim(outfldr)//"/"//file_bsl,time)
        call isos_restart_write(isos,trim(outfldr)//"/"//file_isos,time)
        call yelmo_restart_write(ylmo,trim(outfldr)//"/"//file_yelmo,time) 

        return

    end subroutine yelmox_bipolar_restart_write


    ! Groups:
    subroutine read_parameters(path_par, ctl, yelmox1, yelmox2, tm_1D, tm_2D, tm_2Dsm)

        implicit none

        character(len=512),  intent(OUT) :: path_par     ! parameter namelist of the entire simulation
        type(ctrl_params),   intent(INOUT) :: ctl          ! run control parameters
        type(yelmox_class),  intent(INOUT) :: yelmox1, yelmox2
        type(timeout_class), intent(INOUT) :: tm_1D, tm_2D, tm_2Dsm
        

        ! Determine the parameter file from the command line 
        call yelmo_load_command_line_args(path_par)

        ! Timing and other parameters 
        call nml_read(path_par,"ctrl","time_init",      ctl%time_init)          ! [yr] Starting time
        call nml_read(path_par,"ctrl","time_end",       ctl%time_end)           ! [yr] Ending time
        call nml_read(path_par,"ctrl","time_equil",     ctl%time_equil)         ! [yr] Years to equilibrate first
        call nml_read(path_par,"ctrl","time_const",     ctl%time_const) 
        call nml_read(path_par,"ctrl","time_lgm_step",  ctl%time_lgm_step) 
        call nml_read(path_par,"ctrl","time_pd_step",   ctl%time_pd_step) 
        call nml_read(path_par,"ctrl","dtt",            ctl%dtt)                ! [yr] Main loop time step 
        call nml_read(path_par,"ctrl","dt_restart",     ctl%dt_restart)
        call nml_read(path_par,"ctrl","transient_clim", ctl%transient_clim)     ! Calculate transient climate? 
        call nml_read(path_par,"ctrl","use_lgm_step",   ctl%use_lgm_step)       ! Use lgm_step?
        call nml_read(path_par,"ctrl","use_pd_step",    ctl%use_pd_step)        ! Use pd_step?
        call nml_read(path_par,"ctrl","with_ice_sheet", ctl%with_ice_sheet)     ! Include an active ice sheet 
        call nml_read(path_par,"ctrl","equil_method",   ctl%equil_method)       ! What method should be used for spin-up?
        
        ! bipolar ctl parameters
        call nml_read(path_par, "ctrl",  "active_north",   ctl%active_north)       ! Is the Northern Hemisphere active?
        call nml_read(path_par, "ctrl",  "active_south",   ctl%active_south)       ! Is the Southern Hemisphere active?
        call nml_read(path_par, "ctrl", "active_obm",      ctl%active_obm)
        call nml_read(path_par, "ctrl", "ism2obm",         ctl%ism2obm)
        call nml_read(path_par, "ctrl", "obm2ism",         ctl%obm2ism)
        call nml_read(path_par, "ctrl", "atm2obm",         ctl%atm2obm)
        call nml_read(path_par, "ctrl", "path_par_ocn",    ctl%path_par_ocn)
        call nml_read(path_par, "ctrl", "obm_name",        ctl%obm_name)

        if (ctl%active_north) then
            call nml_read(path_par,"ctrl","path_par_north", yelmox1%path_par)
            call nml_read(path_par,"ctrl","hydro_mask_north", yelmox1%path_hydro_mask)
        end if
        if (ctl%active_south) then
            call nml_read(path_par,"ctrl","path_par_south", yelmox2%path_par)
            call nml_read(path_par,"ctrl","hydro_mask_south", yelmox2%path_hydro_mask)
        end if

        ! Get output times
        call timeout_init(tm_1D,  path_par,"tm_1D",  "small",  ctl%time_init,ctl%time_end)
        call timeout_init(tm_2D,  path_par,"tm_2D",  "heavy",  ctl%time_init,ctl%time_end)
        call timeout_init(tm_2Dsm,path_par,"tm_2Dsm","medium", ctl%time_init,ctl%time_end)

        
        ! Hard-coded for now:
        ctl%dt_clim = 10.0      ! [yrs] Frequency to update snapclim snapshot

        ! Consistency checks ===

        ! transient_clim overrides use_lgm_step 
        if (ctl%transient_clim) ctl%use_lgm_step = .FALSE. 

        ! lgm step should only come after time_equil is finished...
        if (ctl%time_lgm_step .lt. ctl%time_equil) then 
            write(5,*) ""
            write(5,*) "yelmox:: time_lgm_step must be greater than time_equil."
            write(5,*) "Try again."
            stop "Program stopped."
        end if 

        ! pd step should only come after lgm step 
        if (ctl%time_pd_step .lt. ctl%time_lgm_step) then 
            write(5,*) ""
            write(5,*) "yelmox:: time_pd_step must be greater than time_lgm_step."
            write(5,*) "Try again."
            stop "Program stopped."
        end if 

        write(*,*) "yelmox_bipolar :: Consistency checked"

        if (trim(ctl%equil_method) .eq. "opt") then 
            ! Load optimization parameters 

            ! Initially set to zero
            if (ctl%active_north) then
                yelmox1%opt%tf_basins = 0 
                call optimize_par_load(yelmox1%opt,yelmox1%path_par,"opt")
            end if
            if (ctl%active_south) then
                yelmox2%opt%tf_basins = 0 
                call optimize_par_load(yelmox2%opt,yelmox2%path_par,"opt")
            end if

        end if 

        write(*,*) "yelmox_bipolar :: opt parameters loaded"

        return

    end subroutine read_parameters

    subroutine define_input_and_output_locations(outf, hemisphere, yelmox)

        implicit none

        character(len=256), intent(IN)    :: outf
        character(len=5), intent(IN)    :: hemisphere
        type(yelmox_class), intent(INOUT) :: yelmox

        yelmox%file1D              = trim(outf)//trim(hemisphere)//"_yelmo1D.nc"
        yelmox%file2D              = trim(outf)//trim(hemisphere)//"_yelmo2D.nc"
        yelmox%file2D_small        = trim(outf)//trim(hemisphere)//"_yelmo2Dsm.nc"
        yelmox%file_isos           = trim(outf)//trim(hemisphere)//"_fastisostasy.nc"
        
        return

    end subroutine define_input_and_output_locations

    subroutine initialize_icesheet(path_par, hemisphere, time, ctl, yelmox)

        implicit none

        ! Arguments
        character(len=512),  intent(IN)    :: path_par     ! parameter namelist of the entire simulation
        character(len=5),    intent(IN)    :: hemisphere
        real(wp),            intent(IN)    :: time
        type(ctrl_params),   intent(INOUT) :: ctl          ! run control parameters
        type(yelmox_class),  intent(INOUT) :: yelmox

        ! === Initialize ice sheet model =====

        ! Initialize data objects and load initial topography
        call yelmo_init(yelmox%yelmo1,filename=yelmox%path_par,grid_def="file",time=time)
        call rewrite_common_parameters(filename=path_par, dom=yelmox%yelmo1)  ! Rewrite those parameters that are common between domains (remove in future version, spm)
        
        ! Store domain name as a shortcut 
        yelmox%domain = yelmox%yelmo1%par%domain  

        ! Ensure optimization fields are allocated and preassigned
        allocate(yelmox%opt%cf_min(yelmox%yelmo1%grd%nx,yelmox%yelmo1%grd%ny))
        allocate(yelmox%opt%cf_max(yelmox%yelmo1%grd%nx,yelmox%yelmo1%grd%ny))
        
        yelmox%opt%cf_min = yelmox%opt%cf_min_par 
        yelmox%opt%cf_max = yelmox%yelmo1%dyn%par%till_cf_ref

        ! Define specific regions of interest =====================

        ! Shortcut switches for later use
        yelmox%running_laurentide = .FALSE. 
        yelmox%running_greenland  = .FALSE. 

        select case(trim(yelmox%domain))

            case("Antarctica")

                ! Define base regions for whole domain first 
                yelmox%regions_mask_fnm = "ice_data/Antarctica/"//trim(yelmox%yelmo1%par%grid_name)//&
                                    "/"//trim(yelmox%yelmo1%par%grid_name)//"_BASINS-nasa.nc"
                allocate(yelmox%regions_mask(yelmox%yelmo1%grd%nx,yelmox%yelmo1%grd%ny))
                
                ! Load mask from file 
                call nc_read(yelmox%regions_mask_fnm,"mask_regions",yelmox%regions_mask)

                ! APIS region (region=3.0 in regions map)
                yelmox%reg1%write = .TRUE. 
                yelmox%reg1%name  = "APIS" 
                yelmox%reg1%fnm   = trim(outfldr)//trim(hemisphere)//"_yelmo1D_"//trim(yelmox%reg1%name)//".nc"

                allocate(yelmox%reg1%mask(yelmox%yelmo1%grd%nx,yelmox%yelmo1%grd%ny))
                yelmox%reg1%mask = .FALSE. 
                where(abs(yelmox%regions_mask - 3.0) .lt. 1e-3) yelmox%reg1%mask = .TRUE.

                ! WAIS region (region=4.0 in regions map)
                yelmox%reg2%write = .TRUE. 
                yelmox%reg2%name  = "WAIS" 
                yelmox%reg2%fnm   = trim(outfldr)//trim(hemisphere)//"_yelmo1D_"//trim(yelmox%reg2%name)//".nc"

                allocate(yelmox%reg2%mask(yelmox%yelmo1%grd%nx,yelmox%yelmo1%grd%ny))
                yelmox%reg2%mask = .FALSE. 
                where(abs(yelmox%regions_mask - 4.0) .lt. 1e-3) yelmox%reg2%mask = .TRUE.

                ! EAIS region (region=5.0 in regions map)
                yelmox%reg3%write = .TRUE. 
                yelmox%reg3%name  = "EAIS" 
                yelmox%reg3%fnm   = trim(outfldr)//trim(hemisphere)//"_yelmo1D_"//trim(yelmox%reg3%name)//".nc"

                allocate(yelmox%reg3%mask(yelmox%yelmo1%grd%nx,yelmox%yelmo1%grd%ny))
                yelmox%reg3%mask = .FALSE. 
                where(abs(yelmox%regions_mask - 5.0) .lt. 1e-3) yelmox%reg3%mask = .TRUE.

            case("Laurentide")

                yelmox%running_laurentide = .TRUE. 

                ctl%use_lgm_step        = .FALSE.
                ctl%use_pd_step         = .FALSE.

                yelmox%laurentide_init_const_H = .FALSE.

                if (yelmox%laurentide_init_const_H) then 
                    ! Make sure relaxation spinup is short, but transient spinup
                    ! with modified positive smb over North America is reasonably long.

                    ctl%time_equil        = 10.0 
                    yelmox%laurentide_time_equil = 5e3 

                else 
                    ! When starting from ice-6g, positive smb spinup is not necessary.
                    ! however ctl%time_equil should not be changed, as it may be useful
                    ! for spinning up thermodynamics.

                    yelmox%laurentide_time_equil = 0.0 

                end if 

                ! Make sure to set ice_allowed to prevent ice from growing in 
                ! Greenland (and on grid borders)

                where(abs(yelmox%yelmo1%bnd%regions - 1.30) .lt. 1e-3) yelmox%yelmo1%bnd%ice_allowed = .FALSE. 
                
                yelmox%yelmo1%bnd%ice_allowed(1,:)             = .FALSE. 
                yelmox%yelmo1%bnd%ice_allowed(yelmox%yelmo1%grd%nx,:) = .FALSE. 
                yelmox%yelmo1%bnd%ice_allowed(:,1)             = .FALSE. 
                yelmox%yelmo1%bnd%ice_allowed(:,yelmox%yelmo1%grd%ny) = .FALSE. 
                
                ! Hudson region (region=1.12 in regions map)
                yelmox%reg1%write = .TRUE. 
                yelmox%reg1%name  = "Hudson" 
                yelmox%reg1%fnm   = trim(outfldr)//trim(hemisphere)//"_yelmo1D_"//trim(yelmox%reg1%name)//".nc"

                allocate(yelmox%reg1%mask(yelmox%yelmo1%grd%nx,yelmox%yelmo1%grd%ny))
                yelmox%reg1%mask = .FALSE. 
                where(abs(yelmox%yelmo1%bnd%regions - 1.12) .lt. 1e-3) yelmox%reg1%mask = .TRUE.

                yelmox%reg2%write = .FALSE. 
                yelmox%reg3%write = .FALSE. 

            case("Greenland")

                yelmox%running_greenland = .TRUE.

                ! Should extra ice be imposed over continental shelf to mimic LGM state to start
                yelmox%greenland_init_marine_H = .TRUE. 
                
                ! Should glacial smb be modified to reduce negative smb values
                yelmox%scale_glacial_smb = .FALSE. 
                
                ! Should NEGIS parameter modifications be used
                yelmox%ngs%use_negis_par = .TRUE. 

                ! Load NEGIS parameters from file, if used
                if (yelmox%ngs%use_negis_par) then
                    
                    call nml_read(path_par,"negis","cf_0",       yelmox%ngs%cf_0)
                    call nml_read(path_par,"negis","cf_1",       yelmox%ngs%cf_1)
                    call nml_read(path_par,"negis","cf_centre",  yelmox%ngs%cf_centre)
                    call nml_read(path_par,"negis","cf_north",   yelmox%ngs%cf_north)
                    call nml_read(path_par,"negis","cf_north",   yelmox%ngs%cf_north)
                    
                end if 

                ! Make sure to set ice_allowed to prevent ice from growing in 
                ! Iceland and Svaalbard (on grid borders)

                where(abs(yelmox%yelmo1%bnd%regions - 1.20) .lt. 1e-3) yelmox%yelmo1%bnd%ice_allowed = .FALSE. 
                where(abs(yelmox%yelmo1%bnd%regions - 1.23) .lt. 1e-3) yelmox%yelmo1%bnd%ice_allowed = .FALSE. 
                where(abs(yelmox%yelmo1%bnd%regions - 1.31) .lt. 1e-3) yelmox%yelmo1%bnd%ice_allowed = .FALSE.            
                
                if (yelmox%yelmo1%dyn%par%till_method .eq. -1) then 
                    ! Initialize cb_ref to constant value to start
                    ! (will be updated in the time loop)

                    yelmox%yelmo1%dyn%now%cb_ref = yelmox%yelmo1%dyn%par%till_cf_ref 

                end if 

                yelmox%reg1%write = .FALSE. 
                yelmox%reg2%write = .FALSE. 
                yelmox%reg3%write = .FALSE. 
                
            case DEFAULT 

                yelmox%reg1%write = .FALSE. 
                yelmox%reg2%write = .FALSE. 
                yelmox%reg3%write = .FALSE. 

        end select

        return

    end subroutine initialize_icesheet

    subroutine initialize_external_models(path_par, yelmox)

        implicit none

        ! Arguments
        character(len=512),  intent(IN)    :: path_par     ! parameter namelist of the specific domain
        type(yelmox_class),  intent(INOUT) :: yelmox

        ! === Initialize external models (forcing for ice sheet) ======

        ! Initialize fastisosaty
        call isos_init(yelmox%isos1, path_par, "isos", yelmox%yelmo1%grd%nx, yelmox%yelmo1%grd%ny, &
            yelmox%yelmo1%grd%dx, yelmox%yelmo1%grd%dy)

        ! Initialize "climate" model (climate and ocean forcing)
        call snapclim_init(yelmox%snp1,path_par,yelmox%domain,yelmox%yelmo1%par%grid_name,yelmox%yelmo1%grd%nx,yelmox%yelmo1%grd%ny,yelmox%yelmo1%bnd%basins)
        
        ! Initialize surface mass balance model (bnd%smb, bnd%T_srf)
        call smbpal_init(yelmox%smbpal1,path_par,x=yelmox%yelmo1%grd%xc,y=yelmox%yelmo1%grd%yc,lats=yelmox%yelmo1%grd%lat)
        
        ! Initialize marine melt model (bnd%bmb_shlf)
        call marshelf_init(yelmox%mshlf1,path_par,"marine_shelf",yelmox%yelmo1%grd%nx,yelmox%yelmo1%grd%ny,yelmox%domain,yelmox%yelmo1%par%grid_name,yelmox%yelmo1%bnd%regions,yelmox%yelmo1%bnd%basins)
        
        ! Load other constant boundary variables (bnd%H_sed, bnd%Q_geo)
        call sediments_init(yelmox%sed1,path_par,yelmox%yelmo1%grd%nx,yelmox%yelmo1%grd%ny,yelmox%domain,yelmox%yelmo1%par%grid_name)
        yelmox%yelmo1%bnd%H_sed = yelmox%sed1%now%H 
        
        call geothermal_init(yelmox%gthrm1,path_par,yelmox%yelmo1%grd%nx,yelmox%yelmo1%grd%ny,yelmox%domain,yelmox%yelmo1%par%grid_name)
        yelmox%yelmo1%bnd%Q_geo    = yelmox%gthrm1%now%ghf 

        return

    end subroutine initialize_external_models

    subroutine initialize_boundary_conditions(ctl, yelmox, bsl, time, time_bp)

        implicit none

        ! Arguments
        type(ctrl_params),   intent(IN)    :: ctl          ! run control parameters
        type(yelmox_class),  intent(INOUT) :: yelmox
        type(bsl_class),     intent(INOUT) :: bsl
        real(wp),            intent(IN)    :: time, time_bp

        ! === Update initial boundary conditions for current time and yelmo state =====
        ! ybound: z_bed, z_sl, H_sed, smb, T_srf, bmb_shlf , Q_geo

        ! Initialize the isostasy reference state using reference topography fields
        call isos_init_ref(yelmox%isos1, yelmox%yelmo1%bnd%z_bed_ref, yelmox%yelmo1%bnd%H_ice_ref)
        call isos_init_state(yelmox%isos1, yelmox%yelmo1%bnd%z_bed, yelmox%yelmo1%tpo%now%H_ice, time, bsl)
        call isos_write_init_extended(yelmox%isos1, yelmox%file_isos, time)

        yelmox%yelmo1%bnd%z_bed = yelmox%isos1%out%z_bed
        yelmox%yelmo1%bnd%z_sl  = yelmox%isos1%out%z_ss

        ! Update snapclim
        call snapclim_update(yelmox%snp1,z_srf=yelmox%yelmo1%tpo%now%z_srf,time=time_bp,domain=yelmox%domain,dx=yelmox%yelmo1%grd%dx,basins=yelmox%yelmo1%bnd%basins)

        ! Equilibrate snowpack for itm
        if (trim(yelmox%smbpal1%par%abl_method) .eq. "itm") then 
            call smbpal_update_monthly_equil(yelmox%smbpal1,yelmox%snp1%now%tas,yelmox%snp1%now%pr, &
                                yelmox%yelmo1%tpo%now%z_srf,yelmox%yelmo1%tpo%now%H_ice,time_bp,time_equil=100.0)
        end if 

        ! Testing related to present-day surface mass balance
        !     yelmox%snp1%now%tas = yelmox%snp1%now%tas + 0.0 
        !     yelmox%snp1%now%pr  = yelmox%snp1%now%pr*exp(0.05*(0.0))

        !     call smbpal_update_monthly(yelmox%smbpal1,yelmox%snp1%now%tas,yelmox%snp1%now%pr, &
        !                                yelmox%yelmo1%tpo%now%z_srf,yelmox%yelmo1%tpo%now%H_ice,time, &
        !                                file_out=trim(outfldr)//"smbpal.nc",write_now=.TRUE.,write_init=.TRUE.) 

        !     stop 

        call smbpal_update_monthly(yelmox%smbpal1,yelmox%snp1%now%tas,yelmox%snp1%now%pr, &
                                yelmox%yelmo1%tpo%now%z_srf,yelmox%yelmo1%tpo%now%H_ice,time_bp) 
        yelmox%yelmo1%bnd%smb   = yelmox%smbpal1%ann%smb*yelmox%yelmo1%bnd%c%conv_we_ie*1e-3    ! [mm we/a] => [m ie/a]
        yelmox%yelmo1%bnd%T_srf = yelmox%smbpal1%ann%tsrf 

        if (trim(yelmox%yelmo1%par%domain) .eq. "Greenland" .and. yelmox%scale_glacial_smb) then 
            ! Modify glacial smb
            call calc_glacial_smb(yelmox%yelmo1%bnd%smb,yelmox%yelmo1%grd%lat,yelmox%snp1%now%ta_ann,yelmox%snp1%clim0%ta_ann)
        end if

        call marshelf_update_shelf(yelmox%mshlf1,yelmox%yelmo1%tpo%now%H_ice,yelmox%yelmo1%bnd%z_bed,yelmox%yelmo1%tpo%now%f_grnd, &
                            yelmox%yelmo1%bnd%basins,yelmox%yelmo1%bnd%z_sl,yelmox%yelmo1%grd%dx,yelmox%snp1%now%depth, &
                            yelmox%snp1%now%to_ann,yelmox%snp1%now%so_ann,dto_ann=yelmox%snp1%now%to_ann-yelmox%snp1%clim0%to_ann)

        call marshelf_update(yelmox%mshlf1,yelmox%yelmo1%tpo%now%H_ice,yelmox%yelmo1%bnd%z_bed,yelmox%yelmo1%tpo%now%f_grnd, &
                            yelmox%yelmo1%bnd%regions,yelmox%yelmo1%bnd%basins,yelmox%yelmo1%bnd%z_sl,dx=yelmox%yelmo1%grd%dx)

        yelmox%yelmo1%bnd%bmb_shlf = yelmox%mshlf1%now%bmb_shlf
        yelmox%yelmo1%bnd%T_shlf   = yelmox%mshlf1%now%T_shlf  

        call yelmo_print_bound(yelmox%yelmo1%bnd) 
        
        ! Initialize state variables (dyn,therm,mat)
        ! (initialize temps with robin method with a cold base)
        call yelmo_init_state(yelmox%yelmo1,time=time,thrm_method="robin-cold")

        ! ===== basal friction optimization ======
        if (trim(ctl%equil_method) .eq. "opt") then 
            
            ! Ensure that cb_ref will be optimized (till_method == set externally) 
            yelmox%yelmo1%dyn%par%till_method = -1  

            ! If not using restart, prescribe cb_ref to initial guess 
            if (.not. yelmox%yelmo1%par%use_restart) then
                yelmox%yelmo1%dyn%now%cb_ref = yelmox%opt%cf_init 
            end if 

        end if 
        ! ========================================

        if (.not. yelmox%yelmo1%par%use_restart) then 
            ! No restart file used, perform various initialization steps

            if (yelmox%running_laurentide) then 
                ! Start with some ice thickness for testing

                ! Load LGM reconstruction into reference ice thickness
                yelmox%path_lgm = "ice_data/Laurentide/"//trim(yelmox%yelmo1%par%grid_name)//&
                            "/"//trim(yelmox%yelmo1%par%grid_name)//"_TOPO-ICE-6G_C.nc"
                call nc_read(yelmox%path_lgm,"dz",yelmox%yelmo1%bnd%H_ice_ref,start=[1,1,1], &
                                    count=[yelmox%yelmo1%tpo%par%nx,yelmox%yelmo1%tpo%par%ny,1]) 

                ! Start with some ice cover to speed up initialization
                if (yelmox%laurentide_init_const_H) then
                
                    yelmox%yelmo1%tpo%now%H_ice = 0.0
                    where (yelmox%yelmo1%bnd%regions .eq. 1.1 .and. yelmox%yelmo1%bnd%z_bed .gt. 0.0) yelmox%yelmo1%tpo%now%H_ice = 1000.0 
                    where (yelmox%yelmo1%bnd%regions .eq. 1.12) yelmox%yelmo1%tpo%now%H_ice = 1000.0 

                else
                    ! Set LGM reconstruction as initial ice thickness over North America
                    where ( yelmox%yelmo1%bnd%z_bed .gt. -500.0 .and. &
                            (   yelmox%yelmo1%bnd%regions .eq. 1.1  .or. &
                                yelmox%yelmo1%bnd%regions .eq. 1.11 .or. &
                                yelmox%yelmo1%bnd%regions .eq. 1.12) )
                        yelmox%yelmo1%tpo%now%H_ice = yelmox%yelmo1%bnd%H_ice_ref
                    end where 

                    ! Apply Gaussian smoothing to keep things stable
                    call smooth_gauss_2D(yelmox%yelmo1%tpo%now%H_ice,dx=yelmox%yelmo1%grd%dx,f_sigma=2.0)
                
                end if 
                
                ! Load sediment mask 
                yelmox%path_lgm = "ice_data/Laurentide/"//trim(yelmox%yelmo1%par%grid_name)//&
                            "/"//trim(yelmox%yelmo1%par%grid_name)//"_SED-L97.nc"
                call nc_read(yelmox%path_lgm,"z_sed",yelmox%yelmo1%bnd%H_sed) 

                if (ctl%with_ice_sheet) then
                    ! Run Yelmo for briefly to update surface topography
                    call yelmo_update_equil(yelmox%yelmo1,time,time_tot=1.0_prec,dt=1.0,topo_fixed=.TRUE.)

                    ! Addtional cleanup - remove floating ice 
                    where( yelmox%yelmo1%tpo%now%mask_bed .eq. 5) yelmox%yelmo1%tpo%now%H_ice = 0.0 
                    call yelmo_update_equil(yelmox%yelmo1,time,time_tot=1.0_prec,dt=1.0,topo_fixed=.TRUE.)
                end if 

                ! Update snapclim to reflect new topography 
                call snapclim_update(yelmox%snp1,z_srf=yelmox%yelmo1%tpo%now%z_srf,time=time,domain=yelmox%domain,dx=yelmox%yelmo1%grd%dx,basins=yelmox%yelmo1%bnd%basins)

                ! Update smbpal
                call smbpal_update_monthly(yelmox%smbpal1,yelmox%snp1%now%tas,yelmox%snp1%now%pr, &
                                        yelmox%yelmo1%tpo%now%z_srf,yelmox%yelmo1%tpo%now%H_ice,time_bp) 
                yelmox%yelmo1%bnd%smb   = yelmox%smbpal1%ann%smb*yelmox%yelmo1%bnd%c%conv_we_ie*1e-3    ! [mm we/a] => [m ie/a]
                yelmox%yelmo1%bnd%T_srf = yelmox%smbpal1%ann%tsrf 

                if (yelmox%laurentide_init_const_H) then
                    ! Additionally ensure smb is postive for land above 50degN in Laurentide region
                    ! to make sure ice grows everywhere needed (Coridilleran ice sheet mainly)
                    where (yelmox%yelmo1%bnd%regions .eq. 1.1 .and. yelmox%yelmo1%grd%lat .gt. 50.0 .and. &
                            yelmox%yelmo1%bnd%z_bed .gt. 0.0 .and. yelmox%yelmo1%bnd%smb .lt. 0.0 ) yelmox%yelmo1%bnd%smb = 0.5 
                    
                    if (ctl%with_ice_sheet) then
                        ! Run yelmo for several years to ensure stable central ice dome
                        call yelmo_update_equil(yelmox%yelmo1,time,time_tot=5e3,dt=5.0,topo_fixed=.FALSE.)
                    end if 

                else 

                    if (ctl%with_ice_sheet) then
                        ! Run yelmo for several years with constant boundary conditions to stabilize fields
                        call yelmo_update_equil(yelmox%yelmo1,time,time_tot=1e2,dt=5.0,topo_fixed=.FALSE.)
                    end if 

                end if 

            else if (yelmox%running_greenland) then
                ! Special start-up steps for Greenland

                if (yelmox%ngs%use_negis_par) then 

                    ! Ensure till method is correct, since we are updating cb_ref externally
                    yelmox%yelmo1%dyn%par%till_method = -1 

                    ! Update cb_ref using negis parameters 
                    call negis_update_cb_ref(yelmox%yelmo1,yelmox%ngs,time)

                end if 

                if (yelmox%greenland_init_marine_H) then
                    ! Add extra ice-thickness over continental shelf to start with
                    ! an LGM-like state

                    where(yelmox%yelmo1%bnd%ice_allowed .and. yelmox%yelmo1%tpo%now%H_ice .lt. 600.0 &
                            .and. yelmox%yelmo1%bnd%z_bed .gt. -500.0)

                            yelmox%yelmo1%tpo%now%H_ice = 800.0 

                    end where

                    if (ctl%with_ice_sheet) then
                        ! Run yelmo for a few years with constant boundary conditions
                        ! to synchronize all model fields a bit
                        call yelmo_update_equil(yelmox%yelmo1,time,time_tot=10.0_prec, dt=1.0_prec,topo_fixed=.FALSE.)
                    end if 

                end if
                    
            else 
                ! Run simple startup equilibration step 
                
                if (ctl%with_ice_sheet) then
                    ! Run yelmo for a few years with constant boundary conditions
                    ! to synchronize all model fields a bit
                    call yelmo_update_equil(yelmox%yelmo1,time,time_tot=10.0_prec, dt=1.0_prec,topo_fixed=.FALSE.)
                end if 

            end if
            
        end if 

        return

    end subroutine initialize_boundary_conditions

    subroutine initialize_output_files(yelmox, time)

        implicit none

        ! Arguments
        type(yelmox_class),  intent(INOUT) :: yelmox
        real(wp),            intent(IN)    :: time

        call yelmo_write_init(yelmox%yelmo1,yelmox%file2D,time_init=time,units="years") 
        call yelmo_write_reg_init(yelmox%yelmo1,yelmox%file1D,time_init=time,units="years",mask=yelmox%yelmo1%bnd%ice_allowed)
        call yelmo_write_init(yelmox%yelmo1,yelmox%file2D_small,time_init=time,units="years") 

        if (yelmox%reg1%write) then 
            call yelmo_write_reg_init(yelmox%yelmo1,yelmox%reg1%fnm,time_init=time,units="years",mask=yelmox%reg1%mask)
        end if 

        if (yelmox%reg2%write) then
            call yelmo_write_reg_init(yelmox%yelmo1,yelmox%reg2%fnm,time_init=time,units="years",mask=yelmox%reg2%mask)
        end if

        if (yelmox%reg3%write) then
            call yelmo_write_reg_init(yelmox%yelmo1,yelmox%reg3%fnm,time_init=time,units="years",mask=yelmox%reg3%mask)
        end if

        ! Set stats
        yelmox%stats_pd_1%defined = .FALSE. 
        yelmox%stats_lgm%defined  = .FALSE. 
        yelmox%stats_pd_2%defined = .FALSE. 

        return

    end subroutine initialize_output_files

    subroutine spinup_procedure(ctl, yelmox, time, time_elapsed)

        implicit none

        ! Arguments
        type(ctrl_params),   intent(INOUT) :: ctl          ! run control parameters
        type(yelmox_class),  intent(INOUT) :: yelmox
        real(wp),            intent(IN)    :: time, time_elapsed

        select case(trim(ctl%equil_method))
            
            case("opt")

                if (time_elapsed .le. yelmox%opt%rel_time2) then 
                    ! Apply relaxation to the model 

                    ! Update model relaxation time scale and error scaling (in [m])
                    call optimize_set_transient_param(yelmox%opt%rel_tau,time_elapsed,time1=yelmox%opt%rel_time1,time2=yelmox%opt%rel_time2, &
                                                    p1=yelmox%opt%rel_tau1,p2=yelmox%opt%rel_tau2,m=yelmox%opt%rel_m)
                    
                    ! Set model tau, and set yelmo relaxation switch (4: gl line and grounding zone relaxing; 0: no relaxation)
                    yelmox%yelmo1%tpo%par%topo_rel_tau = yelmox%opt%rel_tau 
                    yelmox%yelmo1%tpo%par%topo_rel     = 4
                
                else 
                    ! Turn-off relaxation now

                    yelmox%yelmo1%tpo%par%topo_rel = 0 

                end if 

                ! === Optimization update step =========

                if (yelmox%opt%opt_cf .and. &
                        (time_elapsed .ge. yelmox%opt%cf_time_init .and. time_elapsed .le. yelmox%opt%cf_time_end) ) then
                    ! Perform cf_ref optimization
                
                    ! Update cb_ref based on error metric(s) 
                    call optimize_cb_ref(yelmox%yelmo1%dyn%now%cb_ref,yelmox%yelmo1%tpo%now%H_ice, &
                                            yelmox%yelmo1%tpo%now%dHidt,yelmox%yelmo1%bnd%z_bed,yelmox%yelmo1%bnd%z_sl,yelmox%yelmo1%dyn%now%ux_s,yelmox%yelmo1%dyn%now%uy_s, &
                                            yelmox%yelmo1%dta%pd%H_ice,yelmox%yelmo1%dta%pd%uxy_s,yelmox%yelmo1%dta%pd%H_grnd, &
                                            yelmox%opt%cf_min,yelmox%opt%cf_max,yelmox%yelmo1%tpo%par%dx,yelmox%opt%sigma_err,yelmox%opt%sigma_vel,yelmox%opt%tau_c,yelmox%opt%H0, &
                                            dt=ctl%dtt,fill_method=yelmox%opt%fill_method,fill_dist=yelmox%opt%sigma_err, &
                                            cb_tgt=yelmox%yelmo1%dyn%now%cb_tgt)

                end if

                if (yelmox%opt%opt_tf .and. &
                        (time_elapsed .ge. yelmox%opt%tf_time_init .and. time_elapsed .le. yelmox%opt%tf_time_end) ) then
                    ! Perform tf_corr optimization

                    call optimize_tf_corr(yelmox%mshlf1%now%tf_corr,yelmox%yelmo1%tpo%now%H_ice,yelmox%yelmo1%tpo%now%H_grnd,yelmox%yelmo1%tpo%now%dHidt, &
                                                yelmox%yelmo1%dta%pd%H_ice,yelmox%yelmo1%dta%pd%H_grnd,yelmox%opt%H_grnd_lim,yelmox%opt%tau_m,yelmox%opt%m_temp, &
                                                yelmox%opt%tf_min,yelmox%opt%tf_max,yelmox%yelmo1%tpo%par%dx,sigma=yelmox%opt%tf_sigma,dt=ctl%dtt)
                    ! call optimize_tf_corr(yelmox%mshlf1%now%tf_corr,yelmox%yelmo1%tpo%now%H_ice,yelmox%yelmo1%tpo%now%H_grnd,yelmox%yelmo1%tpo%now%dHidt, &
                    !                         yelmox%yelmo1%dta%pd%H_ice,yelmox%yelmo1%bnd%basins,opt%H_grnd_lim, &
                    !                         opt%tau_m,opt%m_temp,opt%tf_min,opt%tf_max,opt%tf_basins,dt=ctl%dtt)
                
                end if 

            case("relax")
                ! ===== relaxation spinup ==================

                if (time_elapsed .lt. ctl%time_equil) then 
                    ! Turn on relaxation for now, to let thermodynamics equilibrate
                    ! without changing the topography too much. Important when 
                    ! effective pressure = f(thermodynamics).

                    yelmox%yelmo1%tpo%par%topo_rel     = 3
                    yelmox%yelmo1%tpo%par%topo_rel_tau = 50.0 
                    write(*,*) "timelog, tau = ", yelmox%yelmo1%tpo%par%topo_rel_tau

                else if (time_elapsed .eq. ctl%time_equil) then 
                    ! Disable relaxation now... 

                    yelmox%yelmo1%tpo%par%topo_rel     = 0
                    write(*,*) "timelog, relation off..."

                end if 

            case DEFAULT   ! == "none", etc

                ! Pass - do nothing 

        end select 

        return

    end subroutine spinup_procedure

    subroutine update_output(yelmox, bsl, time, tm_1D, tm_2D, tm_2Dsm,hemisphere)

        implicit none

        type(yelmox_class),  intent(IN)    :: yelmox
        type(bsl_class),     intent(IN)    :: bsl
        real(wp),            intent(IN)    :: time
        type(timeout_class), intent(IN)    :: tm_1D, tm_2D, tm_2Dsm
        character(len=5),   intent(IN)     :: hemisphere

        ! == MODEL OUTPUT =======================================================

        if (timeout_check(tm_2D,time)) then
            call yelmox_write_step(yelmox%yelmo1,yelmox%snp1,yelmox%mshlf1,yelmox%smbpal1,yelmox%file2D,time=time)
        end if

        if (timeout_check(tm_2Dsm,time)) then 
            call yelmo_write_step(yelmox%yelmo1,yelmox%file2D_small,time,compare_pd=.FALSE.)
        end if

        if (timeout_check(tm_1D,time)) then 
            call yelmo_write_reg_step(yelmox%yelmo1,yelmox%file1D,time=time)

            if (yelmox%reg1%write) then 
                call yelmo_write_reg_step(yelmox%yelmo1,yelmox%reg1%fnm,time=time,mask=yelmox%reg1%mask)
            end if 

            if (yelmox%reg2%write) then
                call yelmo_write_reg_step(yelmox%yelmo1,yelmox%reg2%fnm,time=time,mask=yelmox%reg2%mask)
            end if

            if (yelmox%reg3%write) then
                call yelmo_write_reg_step(yelmox%yelmo1,yelmox%reg3%fnm,time=time,mask=yelmox%reg3%mask)
            end if

        end if 

        if (mod(nint(time*100),nint(ctl%dt_restart*100))==0) then
            call yelmox_bipolar_restart_write(bsl,yelmox%isos1,yelmox%yelmo1,time,hemisphere=hemisphere)
        end if 

        return

    end subroutine update_output

end program yelmox