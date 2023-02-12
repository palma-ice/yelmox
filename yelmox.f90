

program yelmox

    use ncio 
    use yelmo 
    use yelmo_tools, only : smooth_gauss_2D
    use ice_optimization

    ! External libraries
    use sealevel 
    use isostasy  
     
    use snapclim
    use marine_shelf 
    use smbpal   
    use sediments 
    use geothermal
    
    implicit none 

    type(yelmo_class)      :: yelmo1 
    
    type(sealevel_class)   :: sealev 
    type(snapclim_class)   :: snp1 
    type(marshelf_class)   :: mshlf1 
    type(smbpal_class)     :: smbpal1  
    type(sediments_class)  :: sed1 
    type(geothermal_class) :: gthrm1
    type(isos_class)       :: isos1
    
    character(len=256) :: outfldr, file1D, file2D, file_restart, domain
    character(len=256) :: file2D_small
    character(len=512) :: path_par, path_const
    character(len=512) :: path_lgm  
    real(prec) :: time, time_bp, time_elapsed
    real(wp)   :: dT_now 
    integer    :: n
    
    real(8) :: cpu_start_time, cpu_end_time, cpu_dtime  
    
    type reg_def_class 
        character(len=56)  :: name 
        character(len=512) :: fnm
        logical, allocatable :: mask(:,:) 
        logical :: write 
    end type
    
    type(reg_def_class) :: reg1 
    type(reg_def_class) :: reg2 
    type(reg_def_class) :: reg3 
    
    character(len=512)    :: regions_mask_fnm
    real(wp), allocatable :: regions_mask(:,:) 

    type ctrl_params
        character(len=56) :: run_step
        real(wp) :: time_init
        real(wp) :: time_end
        real(wp) :: time_equil      ! Only for spinup
        real(wp) :: time_const      ! Only for spinup
        real(wp) :: time_lgm_step 
        real(wp) :: time_pd_step
        real(wp) :: dtt
        real(wp) :: dt1D_out
        real(wp) :: dt2D_out
        real(wp) :: dt2D_small_out
        real(wp) :: dt_restart
        real(wp) :: dt_clim

        logical  :: transient_clim
        logical  :: use_lgm_step
        logical  :: use_pd_step
        logical  :: with_ice_sheet 

        character(len=56) :: equil_method
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


    type(ctrl_params)    :: ctl
    type(ice_opt_params) :: opt  
    type(stats_class)    :: stats_pd_1, stats_lgm, stats_pd_2

    ! Internal parameters
    logical  :: running_laurentide
    logical  :: laurentide_init_const_H 
    real(wp) :: laurentide_time_equil 

    logical  :: running_greenland
    logical  :: greenland_init_marine_H
    logical  :: scale_glacial_smb

    real(wp), parameter :: time_lgm = -19050.0_wp  ! [yr CE] == 21 kyr ago 
    real(wp), parameter :: time_pd  =   1950.0_wp  ! [yr CE] ==  0 kyr ago 
    
    ! Initially set to zero 
    opt%tf_basins = 0 

    ! Start timing 
    call yelmo_cpu_time(cpu_start_time)
    
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
    call nml_read(path_par,"ctrl","dt1D_out",       ctl%dt1D_out)           ! [yr] Frequency of 1D output 
    call nml_read(path_par,"ctrl","dt2D_out",       ctl%dt2D_out)           ! [yr] Frequency of 2D output 
    call nml_read(path_par,"ctrl","dt2D_small_out", ctl%dt2D_small_out)     ! [yr] Frequency of small 2D output 
    call nml_read(path_par,"ctrl","dt_restart",     ctl%dt_restart)
    call nml_read(path_par,"ctrl","transient_clim", ctl%transient_clim)     ! Calculate transient climate? 
    call nml_read(path_par,"ctrl","use_lgm_step",   ctl%use_lgm_step)       ! Use lgm_step?
    call nml_read(path_par,"ctrl","use_pd_step",    ctl%use_pd_step)        ! Use pd_step?
    call nml_read(path_par,"ctrl","with_ice_sheet", ctl%with_ice_sheet)     ! Include an active ice sheet 
    call nml_read(path_par,"ctrl","equil_method",   ctl%equil_method)       ! What method should be used for spin-up?

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

    if (trim(ctl%equil_method) .eq. "opt") then 
        ! Load optimization parameters 

        call optimize_par_load(opt,path_par,"opt")

    end if 

    ! Set initial time 
    time    = ctl%time_init 
    time_bp = time - 1950.0_wp 

    ! Assume program is running from the output folder
    outfldr = "./"

    ! Define input and output locations 
    path_const   = trim(outfldr)//"yelmo_const_Earth.nml"
    file1D       = trim(outfldr)//"yelmo1D.nc"
    file2D       = trim(outfldr)//"yelmo2D.nc"
    file_restart = trim(outfldr)//"yelmo_restart.nc"          

    file2D_small = trim(outfldr)//"yelmo2Dsm.nc"
    
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
    write(*,*) "dt1D_out:   ",      ctl%dt1D_out 
    write(*,*) "dt2D_out:   ",      ctl%dt2D_out 
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
    

    ! === Initialize ice sheet model =====

    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)

    ! Initialize data objects and load initial topography
    call yelmo_init(yelmo1,filename=path_par,grid_def="file",time=time)

    ! Store domain name as a shortcut 
    domain = yelmo1%par%domain 

    ! Ensure optimization fields are allocated and preassigned
    allocate(opt%cf_min(yelmo1%grd%nx,yelmo1%grd%ny))
    allocate(opt%cf_max(yelmo1%grd%nx,yelmo1%grd%ny))
    
    opt%cf_min = opt%cf_min_par 
    opt%cf_max = yelmo1%dyn%par%till_cf_ref

    ! Define specific regions of interest =====================

    ! Shortcut switches for later use
    running_laurentide = .FALSE. 
    running_greenland  = .FALSE. 

    select case(trim(domain))

        case("Antarctica")

            ! Define base regions for whole domain first 
            regions_mask_fnm = "ice_data/Antarctica/"//trim(yelmo1%par%grid_name)//&
                                "/"//trim(yelmo1%par%grid_name)//"_BASINS-nasa.nc"
            allocate(regions_mask(yelmo1%grd%nx,yelmo1%grd%ny))
            
            ! Load mask from file 
            call nc_read(regions_mask_fnm,"mask_regions",regions_mask)

            ! APIS region (region=3.0 in regions map)
            reg1%write = .TRUE. 
            reg1%name  = "APIS" 
            reg1%fnm   = trim(outfldr)//"yelmo1D_"//trim(reg1%name)//".nc"

            allocate(reg1%mask(yelmo1%grd%nx,yelmo1%grd%ny))
            reg1%mask = .FALSE. 
            where(abs(regions_mask - 3.0) .lt. 1e-3) reg1%mask = .TRUE.

            ! WAIS region (region=1.0 in regions map)
            reg2%write = .TRUE. 
            reg2%name  = "WAIS" 
            reg2%fnm   = trim(outfldr)//"yelmo1D_"//trim(reg2%name)//".nc"

            allocate(reg2%mask(yelmo1%grd%nx,yelmo1%grd%ny))
            reg2%mask = .FALSE. 
            where(abs(regions_mask - 1.0) .lt. 1e-3) reg2%mask = .TRUE.

            ! EAIS region (region=2.0 in regions map)
            reg3%write = .TRUE. 
            reg3%name  = "EAIS" 
            reg3%fnm   = trim(outfldr)//"yelmo1D_"//trim(reg3%name)//".nc"

            allocate(reg3%mask(yelmo1%grd%nx,yelmo1%grd%ny))
            reg3%mask = .FALSE. 
            where(abs(regions_mask - 2.0) .lt. 1e-3) reg3%mask = .TRUE.

        case("Laurentide")

            running_laurentide = .TRUE. 

            ctl%use_lgm_step        = .FALSE.
            ctl%use_pd_step         = .FALSE.

            laurentide_init_const_H = .FALSE.

            if (laurentide_init_const_H) then 
                ! Make sure relaxation spinup is short, but transient spinup
                ! with modified positive smb over North America is reasonably long.

                ctl%time_equil        = 10.0 
                laurentide_time_equil = 5e3 

            else 
                ! When starting from ice-6g, positive smb spinup is not necessary.
                ! however ctl%time_equil should not be changed, as it may be useful
                ! for spinning up thermodynamics.

                laurentide_time_equil = 0.0 

            end if 

            ! Make sure to set ice_allowed to prevent ice from growing in 
            ! Greenland (and on grid borders)

            where(abs(yelmo1%bnd%regions - 1.30) .lt. 1e-3) yelmo1%bnd%ice_allowed = .FALSE. 
            
            yelmo1%bnd%ice_allowed(1,:)             = .FALSE. 
            yelmo1%bnd%ice_allowed(yelmo1%grd%nx,:) = .FALSE. 
            yelmo1%bnd%ice_allowed(:,1)             = .FALSE. 
            yelmo1%bnd%ice_allowed(:,yelmo1%grd%ny) = .FALSE. 
            
            ! Hudson region (region=1.12 in regions map)
            reg1%write = .TRUE. 
            reg1%name  = "Hudson" 
            reg1%fnm   = trim(outfldr)//"yelmo1D_"//trim(reg1%name)//".nc"

            allocate(reg1%mask(yelmo1%grd%nx,yelmo1%grd%ny))
            reg1%mask = .FALSE. 
            where(abs(yelmo1%bnd%regions - 1.12) .lt. 1e-3) reg1%mask = .TRUE.

            reg2%write = .FALSE. 
            reg3%write = .FALSE. 

        case("Greenland")

            running_greenland = .TRUE.

            ! Should extra ice be imposed over continental shelf to mimic LGM state to start
            greenland_init_marine_H = .TRUE. 

            ! Should glacial smb be modified to reduce negative smb values
            scale_glacial_smb = .FALSE. 
            
            ! Make sure to set ice_allowed to prevent ice from growing in 
            ! Iceland and Svaalbard (on grid borders)

            where(abs(yelmo1%bnd%regions - 1.20) .lt. 1e-3) yelmo1%bnd%ice_allowed = .FALSE. 
            where(abs(yelmo1%bnd%regions - 1.23) .lt. 1e-3) yelmo1%bnd%ice_allowed = .FALSE. 
            where(abs(yelmo1%bnd%regions - 1.31) .lt. 1e-3) yelmo1%bnd%ice_allowed = .FALSE.            
            
        case DEFAULT 

            reg1%write = .FALSE. 
            reg2%write = .FALSE. 
            reg3%write = .FALSE. 

    end select

    ! === Initialize external models (forcing for ice sheet) ======

    ! Initialize global sea level model (bnd%z_sl)
    call sealevel_init(sealev,path_par)

    ! Initialize bedrock model 
    call isos_init(isos1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,yelmo1%grd%dx)

    ! Initialize "climate" model (climate and ocean forcing)
    call snapclim_init(snp1,path_par,domain,yelmo1%par%grid_name,yelmo1%grd%nx,yelmo1%grd%ny,yelmo1%bnd%basins)
    
    ! Initialize surface mass balance model (bnd%smb, bnd%T_srf)
    call smbpal_init(smbpal1,path_par,x=yelmo1%grd%xc,y=yelmo1%grd%yc,lats=yelmo1%grd%lat)
    
    ! Initialize marine melt model (bnd%bmb_shlf)
    call marshelf_init(mshlf1,path_par,"marine_shelf",yelmo1%grd%nx,yelmo1%grd%ny,domain,yelmo1%par%grid_name,yelmo1%bnd%regions,yelmo1%bnd%basins)
    
    ! Load other constant boundary variables (bnd%H_sed, bnd%Q_geo)
    call sediments_init(sed1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,yelmo1%par%grid_name)
    call geothermal_init(gthrm1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,yelmo1%par%grid_name)
    ! === Update initial boundary conditions for current time and yelmo state =====
    ! ybound: z_bed, z_sl, H_sed, smb, T_srf, bmb_shlf , Q_geo

    ! Initialize isostasy using present-day topography 
    ! values to calibrate the reference rebound
    call isos_init_state(isos1,z_bed=yelmo1%bnd%z_bed,H_ice=yelmo1%tpo%now%H_ice, &
                                    z_sl=yelmo1%bnd%z_sl,z_bed_ref=yelmo1%bnd%z_bed_ref, &
                                    H_ice_ref=yelmo1%bnd%H_ice_ref, &
                                    z_sl_ref=yelmo1%bnd%z_sl*0.0,time=time)
                                       

    call sealevel_update(sealev,year_bp=time_bp)
    yelmo1%bnd%z_sl  = sealev%z_sl 
    yelmo1%bnd%H_sed = sed1%now%H 

    ! Update snapclim
    call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=time_bp,domain=domain,dx=yelmo1%grd%dx,basins=yelmo1%bnd%basins)

    ! Equilibrate snowpack for itm
    if (trim(smbpal1%par%abl_method) .eq. "itm") then 
        call smbpal_update_monthly_equil(smbpal1,snp1%now%tas,snp1%now%pr, &
                               yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time_bp,time_equil=100.0)
    end if 

! Testing related to present-day surface mass balance
!     snp1%now%tas = snp1%now%tas + 0.0 
!     snp1%now%pr  = snp1%now%pr*exp(0.05*(0.0))

!     call smbpal_update_monthly(smbpal1,snp1%now%tas,snp1%now%pr, &
!                                yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time, &
!                                file_out=trim(outfldr)//"smbpal.nc",write_now=.TRUE.,write_init=.TRUE.) 

!     stop 

    call smbpal_update_monthly(smbpal1,snp1%now%tas,snp1%now%pr, &
                               yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time_bp) 
    yelmo1%bnd%smb   = smbpal1%ann%smb*conv_we_ie*1e-3    ! [mm we/a] => [m ie/a]
    yelmo1%bnd%T_srf = smbpal1%ann%tsrf 

    if (trim(yelmo1%par%domain) .eq. "Greenland" .and. scale_glacial_smb) then 
        ! Modify glacial smb
        call calc_glacial_smb(yelmo1%bnd%smb,yelmo1%grd%lat,snp1%now%ta_ann,snp1%clim0%ta_ann)
    end if

!     yelmo1%bnd%smb   = yelmo1%dta%pd%smb
!     yelmo1%bnd%T_srf = yelmo1%dta%pd%t2m
    
    call marshelf_update_shelf(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                        yelmo1%bnd%basins,yelmo1%bnd%z_sl,yelmo1%grd%dx,snp1%now%depth, &
                        snp1%now%to_ann,snp1%now%so_ann,dto_ann=snp1%now%to_ann-snp1%clim0%to_ann)

    call marshelf_update(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                         yelmo1%bnd%regions,yelmo1%bnd%basins,yelmo1%bnd%z_sl,dx=yelmo1%grd%dx)

    yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf
    yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  

    yelmo1%bnd%Q_geo    = gthrm1%now%ghf 
    
    call yelmo_print_bound(yelmo1%bnd) 
    
    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base)
    call yelmo_init_state(yelmo1,time=time,thrm_method="robin-cold")

    if (yelmo1%par%use_restart) then 
        ! If using restart file, set boundary module variables equal to restarted value 
         
        isos1%now%z_bed  = yelmo1%bnd%z_bed
      
    end if 

    ! ===== basal friction optimization ======
    if (trim(ctl%equil_method) .eq. "opt") then 
        
        ! Ensure that cb_ref will be optimized (till_method == set externally) 
        yelmo1%dyn%par%till_method = -1  

        ! If not using restart, prescribe cb_ref to initial guess 
        if (.not. yelmo1%par%use_restart) then
            yelmo1%dyn%now%cb_ref = opt%cf_init 
        end if 

    end if 
    ! ========================================

    if (ctl%with_ice_sheet .and. (.not. yelmo1%par%use_restart)) then 
        ! Ice sheet is active, and we have not loaded a restart file 

        if (running_laurentide) then 
            ! Start with some ice thickness for testing

            ! Load LGM reconstruction into reference ice thickness
            path_lgm = "ice_data/Laurentide/"//trim(yelmo1%par%grid_name)//&
                        "/"//trim(yelmo1%par%grid_name)//"_TOPO-ICE-6G_C.nc"
            call nc_read(path_lgm,"dz",yelmo1%bnd%H_ice_ref,start=[1,1,1], &
                                count=[yelmo1%tpo%par%nx,yelmo1%tpo%par%ny,1]) 

            ! Start with some ice cover to speed up initialization
            if (laurentide_init_const_H) then
            
                yelmo1%tpo%now%H_ice = 0.0
                where (yelmo1%bnd%regions .eq. 1.1 .and. yelmo1%bnd%z_bed .gt. 0.0) yelmo1%tpo%now%H_ice = 1000.0 
                where (yelmo1%bnd%regions .eq. 1.12) yelmo1%tpo%now%H_ice = 1000.0 

            else
                ! Set LGM reconstruction as initial ice thickness over North America
                where ( yelmo1%bnd%z_bed .gt. -500.0 .and. &
                        (   yelmo1%bnd%regions .eq. 1.1  .or. &
                            yelmo1%bnd%regions .eq. 1.11 .or. &
                            yelmo1%bnd%regions .eq. 1.12) )
                    yelmo1%tpo%now%H_ice = yelmo1%bnd%H_ice_ref
                end where 

                ! Apply Gaussian smoothing to keep things stable
                call smooth_gauss_2D(yelmo1%tpo%now%H_ice,dx=yelmo1%grd%dx,f_sigma=2.0)
            
            end if 
            
            ! Load sediment mask 
            path_lgm = "ice_data/Laurentide/"//trim(yelmo1%par%grid_name)//&
                        "/"//trim(yelmo1%par%grid_name)//"_SED-L97.nc"
            call nc_read(path_lgm,"z_sed",yelmo1%bnd%H_sed) 

            ! Run Yelmo for briefly to update surface topography
            call yelmo_update_equil(yelmo1,time,time_tot=1.0_prec,dt=1.0,topo_fixed=.TRUE.)

            ! Addtional cleanup - remove floating ice 
            where( yelmo1%tpo%now%mask_bed .eq. 5) yelmo1%tpo%now%H_ice = 0.0 
            call yelmo_update_equil(yelmo1,time,time_tot=1.0_prec,dt=1.0,topo_fixed=.TRUE.)

            ! Update snapclim to reflect new topography 
            call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=time,domain=domain,dx=yelmo1%grd%dx,basins=yelmo1%bnd%basins)

            ! Update smbpal
            call smbpal_update_monthly(smbpal1,snp1%now%tas,snp1%now%pr, &
                                       yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time_bp) 
            yelmo1%bnd%smb   = smbpal1%ann%smb*conv_we_ie*1e-3    ! [mm we/a] => [m ie/a]
            yelmo1%bnd%T_srf = smbpal1%ann%tsrf 

            if (laurentide_init_const_H) then
                ! Additionally ensure smb is postive for land above 50degN in Laurentide region
                ! to make sure ice grows everywhere needed (Coridilleran ice sheet mainly)
                where (yelmo1%bnd%regions .eq. 1.1 .and. yelmo1%grd%lat .gt. 50.0 .and. &
                        yelmo1%bnd%z_bed .gt. 0.0 .and. yelmo1%bnd%smb .lt. 0.0 ) yelmo1%bnd%smb = 0.5 
                
                ! Run yelmo for several years to ensure stable central ice dome
                call yelmo_update_equil(yelmo1,time,time_tot=5e3,dt=5.0,topo_fixed=.FALSE.)

            else 

                ! Run yelmo for several years with constant boundary conditions to stabilize fields
                call yelmo_update_equil(yelmo1,time,time_tot=1e2,dt=5.0,topo_fixed=.FALSE.)

            end if 

        else if (running_greenland) then
            ! Special start-up steps for Greenland

            if (greenland_init_marine_H) then
                ! Add extra ice-thickness over continental shelf to start with
                ! an LGM-like state

                where(yelmo1%bnd%ice_allowed .and. yelmo1%tpo%now%H_ice .lt. 600.0 &
                        .and. yelmo1%bnd%z_bed .gt. -500.0)

                        yelmo1%tpo%now%H_ice = 800.0 

                end where

                ! Run yelmo for a few years with constant boundary conditions
                ! to synchronize all model fields a bit
                call yelmo_update_equil(yelmo1,time,time_tot=10.0_prec, dt=1.0_prec,topo_fixed=.FALSE.)

            end if
                
        else 
            ! Run simple startup equilibration step 
            
            ! Run yelmo for a few years with constant boundary conditions
            ! to synchronize all model fields a bit
            call yelmo_update_equil(yelmo1,time,time_tot=10.0_prec, dt=1.0_prec,topo_fixed=.FALSE.)

        end if 

    else 
        ! Starting from restart file 

        call yelmo_update_equil(yelmo1,time,time_tot=10.0_prec, dt=1.0_prec,topo_fixed=.FALSE.,tpo_solver="none")

    end if 

    ! ===== Initialize output files ===== 
    
    call yelmo_write_init(yelmo1,file2D,time_init=time,units="years") 
    call yelmo_write_reg_init(yelmo1,file1D,time_init=time,units="years",mask=yelmo1%bnd%ice_allowed)
    
    call yelmo_write_init(yelmo1,file2D_small,time_init=time,units="years") 
    
    if (reg1%write) then 
        call yelmo_write_reg_init(yelmo1,reg1%fnm,time_init=time,units="years",mask=reg1%mask)
    end if 

    if (reg2%write) then
        call yelmo_write_reg_init(yelmo1,reg2%fnm,time_init=time,units="years",mask=reg2%mask)
    end if

    if (reg3%write) then
        call yelmo_write_reg_init(yelmo1,reg3%fnm,time_init=time,units="years",mask=reg3%mask)
    end if

    ! Set stats
    stats_pd_1%defined = .FALSE. 
    stats_lgm%defined  = .FALSE. 
    stats_pd_2%defined = .FALSE. 
    
    ! ==== Begin main time loop =====

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
        select case(trim(ctl%equil_method))
            
            case("opt")

                if (time_elapsed .le. opt%rel_time2) then 
                    ! Apply relaxation to the model 

                    ! Update model relaxation time scale and error scaling (in [m])
                    call optimize_set_transient_param(opt%rel_tau,time_elapsed,time1=opt%rel_time1,time2=opt%rel_time2, &
                                                    p1=opt%rel_tau1,p2=opt%rel_tau2,m=opt%rel_m)
                    
                    ! Set model tau, and set yelmo relaxation switch (4: gl line and grounding zone relaxing; 0: no relaxation)
                    yelmo1%tpo%par%topo_rel_tau = opt%rel_tau 
                    yelmo1%tpo%par%topo_rel     = 4
                
                else 
                    ! Turn-off relaxation now

                    yelmo1%tpo%par%topo_rel = 0 

                end if 

                ! === Optimization update step =========

                if (time_elapsed .le. opt%cf_time) then 
                    ! Perform cf_ref optimization
                
                    ! Update cb_ref based on error metric(s) 
                    call optimize_cb_ref(yelmo1%dyn%now%cb_ref,yelmo1%tpo%now%H_ice, &
                                            yelmo1%tpo%now%dHidt,yelmo1%bnd%z_bed,yelmo1%bnd%z_sl,yelmo1%dyn%now%ux_s,yelmo1%dyn%now%uy_s, &
                                            yelmo1%dta%pd%H_ice,yelmo1%dta%pd%uxy_s,yelmo1%dta%pd%H_grnd, &
                                            opt%cf_min,opt%cf_max,yelmo1%tpo%par%dx,opt%sigma_err,opt%sigma_vel,opt%tau_c,opt%H0, &
                                            dt=ctl%dtt,fill_method=opt%fill_method,fill_dist=opt%sigma_err, &
                                            cb_tgt=yelmo1%dyn%now%cb_tgt)

                end if

                if (opt%opt_tf .and. time_elapsed .le. opt%tf_time) then
                    ! Perform tf_corr optimization

                    call optimize_tf_corr(mshlf1%now%tf_corr,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%H_grnd,yelmo1%tpo%now%dHidt, &
                                                yelmo1%dta%pd%H_ice,yelmo1%dta%pd%H_grnd,opt%H_grnd_lim,opt%tau_m,opt%m_temp, &
                                                opt%tf_min,opt%tf_max,yelmo1%tpo%par%dx,sigma=opt%tf_sigma,dt=ctl%dtt)
                    ! call optimize_tf_corr(mshlf1%now%tf_corr,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%H_grnd,yelmo1%tpo%now%dHidt, &
                    !                         yelmo1%dta%pd%H_ice,yelmo1%bnd%basins,opt%H_grnd_lim, &
                    !                         opt%tau_m,opt%m_temp,opt%tf_min,opt%tf_max,opt%tf_basins,dt=ctl%dtt)
                
                end if 

            case("relax")
                ! ===== relaxation spinup ==================

                if (time_elapsed .lt. ctl%time_equil) then 
                    ! Turn on relaxation for now, to let thermodynamics equilibrate
                    ! without changing the topography too much. Important when 
                    ! effective pressure = f(thermodynamics).

                    yelmo1%tpo%par%topo_rel     = 3
                    yelmo1%tpo%par%topo_rel_tau = 50.0 
                    write(*,*) "timelog, tau = ", yelmo1%tpo%par%topo_rel_tau

                else if (time_elapsed .eq. ctl%time_equil) then 
                    ! Disable relaxation now... 

                    yelmo1%tpo%par%topo_rel     = 0
                    write(*,*) "timelog, relation off..."

                end if 

            case DEFAULT   ! == "none", etc

                ! Pass - do nothing 

        end select 

        ! == SEA LEVEL ==========================================================
        call sealevel_update(sealev,year_bp=time_bp)
        yelmo1%bnd%z_sl  = sealev%z_sl 

        ! == ISOSTASY ==========================================================
        call isos_update(isos1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_sl,time) 
        yelmo1%bnd%z_bed = isos1%now%z_bed
         
        ! == ICE SHEET ===================================================
        if (ctl%with_ice_sheet) call yelmo_update(yelmo1,time)

        
        ! == CLIMATE (ATMOSPHERE AND OCEAN) ====================================
        
        if (mod(nint(time*100),nint(ctl%dt_clim*100))==0) then
            ! Update snapclim
            call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=time_bp,domain=domain,dx=yelmo1%grd%dx,basins=yelmo1%bnd%basins) 
        end if 

        ! == SURFACE MASS BALANCE ==============================================

        ! ajr: just testing...
        dT_now = 0.0 
        !if (time .gt. 7000.0) dT_now = 4.0 

        call smbpal_update_monthly(smbpal1,snp1%now%tas+dT_now,snp1%now%pr, &
                                   yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time_bp) 
        yelmo1%bnd%smb   = smbpal1%ann%smb*conv_we_ie*1e-3       ! [mm we/a] => [m ie/a]
        yelmo1%bnd%T_srf = smbpal1%ann%tsrf 

        if (trim(yelmo1%par%domain) .eq. "Greenland" .and. scale_glacial_smb) then 
            ! Modify glacial smb
            call calc_glacial_smb(yelmo1%bnd%smb,yelmo1%grd%lat,snp1%now%ta_ann,snp1%clim0%ta_ann)
        end if
    
        ! yelmo1%bnd%smb   = yelmo1%dta%pd%smb
        ! yelmo1%bnd%T_srf = yelmo1%dta%pd%t2m
        
        if (running_laurentide .and. laurentide_init_const_H  &
                .and. (time-ctl%time_init) .lt. laurentide_time_equil ) then 
            ! Additionally ensure smb is postive for land above 50degN in Laurentide region
            ! to make sure ice grows everywhere needed (Coridilleran ice sheet mainly)
            where (yelmo1%bnd%regions .eq. 1.1 .and. yelmo1%grd%lat .gt. 50.0 .and. &
                        yelmo1%bnd%z_bed .gt. 0.0 .and. yelmo1%bnd%smb .lt. 0.0 ) yelmo1%bnd%smb = 0.5 
        end if 

        ! == MARINE AND TOTAL BASAL MASS BALANCE ===============================
        call marshelf_update_shelf(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                        yelmo1%bnd%basins,yelmo1%bnd%z_sl,yelmo1%grd%dx,snp1%now%depth, &
                        snp1%now%to_ann,snp1%now%so_ann,dto_ann=snp1%now%to_ann-snp1%clim0%to_ann)

        call marshelf_update(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                             yelmo1%bnd%regions,yelmo1%bnd%basins,yelmo1%bnd%z_sl,dx=yelmo1%grd%dx)

        yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf
        yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf
        
        ! == MODEL OUTPUT =======================================================

        if (mod(nint(time*100),nint(ctl%dt2D_out*100))==0) then
            call write_step_2D_combined(yelmo1,isos1,snp1,mshlf1,smbpal1,file2D,time=time)
        end if

        if (mod(nint(time*100),nint(ctl%dt2D_small_out*100))==0) then
               call yelmo_write_step(yelmo1,file2D_small,time,compare_pd=.FALSE.)
           end if

        if (mod(nint(time*100),nint(ctl%dt1D_out*100))==0) then
            call yelmo_write_reg_step(yelmo1,file1D,time=time)

            if (reg1%write) then 
                call yelmo_write_reg_step(yelmo1,reg1%fnm,time=time,mask=reg1%mask)
            end if 

            if (reg2%write) then
                call yelmo_write_reg_step(yelmo1,reg2%fnm,time=time,mask=reg2%mask)
            end if

            if (reg3%write) then
                call yelmo_write_reg_step(yelmo1,reg3%fnm,time=time,mask=reg3%mask)
            end if

        end if 

        if (mod(nint(time*100),nint(ctl%dt_restart*100))==0) then 
            call yelmo_restart_write(yelmo1,file_restart,time=time) 
        end if 

        if (mod(time,10.0)==0 .and. (.not. yelmo_log)) then
            write(*,"(a,f14.4)") "yelmo:: time = ", time
        end if 
        
    end do 

    ! Write the restart file for the end of the simulation
    call yelmo_restart_write(yelmo1,file_restart,time=time) 

    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)

    write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(ctl%time_end-ctl%time_init))/(cpu_dtime/3600.0), " kiloyears / hr"

contains
    
    subroutine write_step_2D_combined(ylmo,isos,snp,mshlf,srf,filename,time)

        implicit none 
        
        type(yelmo_class),      intent(IN) :: ylmo
        type(isos_class),       intent(IN) :: isos 
        type(snapclim_class),   intent(IN) :: snp 
        type(marshelf_class),   intent(IN) :: mshlf 
        type(smbpal_class),     intent(IN) :: srf 
        !type(sediments_class),  intent(IN) :: sed 
        !type(geothermal_class), intent(IN) :: gthrm
        !type(isos_class),       intent(IN) :: isos
        
        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: time

        ! Local variables
        integer    :: ncid, n
        real(prec) :: time_prev 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Write model metrics (model speed, dt, eta)
        call yelmo_write_step_model_metrics(filename,ylmo,n,ncid)

        ! Write present-day data metrics (rmse[H],etc)
        call yelmo_write_step_pd_metrics(filename,ylmo,n,ncid)
        
        if (n .eq. 1) then
            call nc_write(filename,"ice_allowed",ylmo%bnd%ice_allowed,units="",long_name="Ice allowed mask", &
                        dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
        end if 
        
        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"mb_applied",ylmo%tpo%now%mb_applied,units="m",long_name="Applied net mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mb_resid",ylmo%tpo%now%mb_resid,units="m",long_name="Residual mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"H_grnd",ylmo%tpo%now%H_grnd,units="m",long_name="Ice thickness overburden", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="bar",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"calv",ylmo%tpo%now%calv,units="m/a ice equiv.",long_name="Calving rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"calv_grnd",ylmo%tpo%now%calv_grnd,units="m/a ice equiv.",long_name="Calving rate (floating)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"calv_flt",ylmo%tpo%now%calv_flt,units="m/a ice equiv.",long_name="Calving rate (grounded)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"f_grnd",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",ylmo%tpo%now%f_ice,units="1",long_name="Ice fraction in grid cell", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"dist_grline",ylmo%tpo%now%dist_grline,units="km", &
!                       long_name="Distance to nearest grounding-line point", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"cb_ref",ylmo%dyn%now%cb_ref,units="--",long_name="Bed friction scalar", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"c_bed",ylmo%dyn%now%c_bed,units="Pa",long_name="Bed friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m^-1",long_name="Basal friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff_int",ylmo%dyn%now%visc_eff_int,units="Pa a m",long_name="Depth-integrated effective viscosity", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud",ylmo%dyn%now%taud,units="Pa",long_name="Driving stress", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)        
        call nc_write(filename,"taub",ylmo%dyn%now%taub,units="Pa",long_name="Basal stress", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Strain-rate and stress tensors 
        if (.FALSE.) then

            ! call nc_write(filename,"de",ylmo%mat%now%strn%de,units="a^-1",long_name="Effective strain rate", &
            !           dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
            ! call nc_write(filename,"te",ylmo%mat%now%strs%te,units="Pa",long_name="Effective stress", &
            !           dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
            ! call nc_write(filename,"visc_int",ylmo%mat%now%visc_int,units="Pa a m",long_name="Depth-integrated effective viscosity (SSA)", &
            !           dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
         
            call nc_write(filename,"de2D",ylmo%mat%now%strn2D%de,units="yr^-1",long_name="Effective strain rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"div2D",ylmo%mat%now%strn2D%div,units="yr^-1",long_name="Divergence strain rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"te2D",ylmo%mat%now%strs2D%te,units="Pa",long_name="Effective stress", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            
            call nc_write(filename,"eps_eig_1",ylmo%mat%now%strn2D%eps_eig_1,units="1/yr",long_name="Eigen strain 1", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"eps_eig_2",ylmo%mat%now%strn2D%eps_eig_2,units="1/yr",long_name="Eigen strain 2", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"eps_eff",ylmo%tpo%now%eps_eff,units="yr^-1",long_name="Effective calving strain", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            
            call nc_write(filename,"tau_eig_1",ylmo%mat%now%strs2D%tau_eig_1,units="Pa",long_name="Eigen stress 1", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"tau_eig_2",ylmo%mat%now%strs2D%tau_eig_2,units="Pa",long_name="Eigen stress 2", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"tau_eff",ylmo%tpo%now%tau_eff,units="Pa",long_name="Effective calving stress", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            
        end if 

!         call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically-averaged velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
                     dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
!        call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
!                     dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!        call nc_write(filename,"T_prime_b",ylmo%thrm%now%T_prime_b,units="deg C",long_name="Homologous basal ice temperature", &
!                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!        call nc_write(filename,"uz",ylmo%dyn%now%uz,units="m/a",long_name="Vertical velocity (z)", &
!                       dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",start=[1,1,1,n],ncid=ncid)
        
!        call nc_write(filename,"uz_star",ylmo%thrm%now%uz_star,units="m yr-1",long_name="Advection-adjusted vertical velocity", &
!                      dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",start=[1,1,1,n],ncid=ncid)
        
!        call nc_write(filename,"T_rock",ylmo%thrm%now%T_rock,units="K",long_name="Bedrock temperature", &
!                      dim1="xc",dim2="yc",dim3="zeta_rock",dim4="time",start=[1,1,1,n],ncid=ncid)
        
 !       call nc_write(filename,"Q_rock",ylmo%thrm%now%Q_rock,units="mW m-2",long_name="Bedrock surface heat flux", &
 !                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!        call nc_write(filename,"Q_ice_b",ylmo%thrm%now%Q_ice_b,units="mW m-2",long_name="Basal ice heat flux", &
!                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Q_strn",ylmo%thrm%now%Q_strn/(rho_ice*ylmo%thrm%now%cp),units="K a-1",long_name="Strain heating", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="mW m-2",long_name="Basal frictional heating", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"bmb_grnd",ylmo%thrm%now%bmb_grnd,units="m/a ice equiv.",long_name="Basal mass balance (grounded)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_w",ylmo%thrm%now%H_w,units="m",long_name="Basal water layer thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        !call nc_write(filename,"ATT",ylmo%mat%now%ATT,units="a^-1 Pa^-3",long_name="Rate factor", &
        !              dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"enh_bar",ylmo%mat%now%enh_bar,units="1",long_name="Vertically averaged enhancement factor", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_int",ylmo%mat%now%visc_int,units="Pa a m",long_name="Vertically integrated viscosity", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Boundaries
        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid) 
        call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"H_sed",ylmo%bnd%H_sed,units="m",long_name="Sediment thickness", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb",ylmo%bnd%smb,units="m/a ice equiv.",long_name="Surface mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb_errpd",ylmo%bnd%smb-ylmo%dta%pd%smb,units="m/a ice equiv.",long_name="Surface mass balance error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_srf",ylmo%bnd%T_srf,units="K",long_name="Surface temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/a ice equiv.",long_name="Basal mass balance (shelf)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        if (n .le. 1) then 
            call nc_write(filename,"H_sed",ylmo%bnd%Q_geo,units="m",long_name="Sediment thickness", &
                        dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
            call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW/m^2",long_name="Geothermal heat flux", &
                        dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
        end if 

        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a ice equiv.",long_name="Basal mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"fmb",ylmo%tpo%now%fmb,units="m/a ice equiv.",long_name="Margin-front mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! External data
        call nc_write(filename,"dzbdt",isos%now%dzbdt,units="m/a",long_name="Bedrock uplift rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"Ta_ann",snp%now%ta_ann,units="K",long_name="Near-surface air temperature (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Ta_sum",snp%now%ta_sum,units="K",long_name="Near-surface air temperature (sum)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Pr_ann",snp%now%pr_ann*1e-3,units="m/a water equiv.",long_name="Precipitation (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
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

        end if 

        !call nc_write(filename,"pr",snp%now%pr*1e-3,units="m/a water equiv.",long_name="Precipitation (ann)", &
        !              dim1="xc",dim2="yc",dim3="month",dim4="time",start=[1,1,1,n],ncid=ncid)
              
        call nc_write(filename,"dTa_ann",snp%now%ta_ann-snp%clim0%ta_ann,units="K",long_name="Near-surface air temperature anomaly (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dTa_sum",snp%now%ta_sum-snp%clim0%ta_sum,units="K",long_name="Near-surface air temperature anomaly (sum)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dPr_ann",(snp%now%pr_ann-snp%clim0%pr_ann)*1e-3,units="m/a water equiv.",long_name="Precipitation anomaly (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"PDDs",srf%ann%PDDs,units="degC days",long_name="Positive degree days (annual total)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Comparison with present-day 
        call nc_write(filename,"H_ice_pd_err",ylmo%dta%pd%err_H_ice,units="m",long_name="Ice thickness error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf_pd_err",ylmo%dta%pd%err_z_srf,units="m",long_name="Surface elevation error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s_pd_err",ylmo%dta%pd%err_uxy_s,units="m/a",long_name="Surface velocity error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"ssa_mask_acx",ylmo%dyn%now%ssa_mask_acx,units="1",long_name="SSA mask (acx)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"ssa_mask_acy",ylmo%dyn%now%ssa_mask_acy,units="1",long_name="SSA mask (acy)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"dzsdx",ylmo%tpo%now%dzsdx,units="m/m",long_name="Surface slope", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"dzsdy",ylmo%tpo%now%dzsdy,units="m/m",long_name="Surface slope", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"f_grnd_acx",ylmo%tpo%now%f_grnd_acx,units="1",long_name="Grounded fraction (acx)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"f_grnd_acy",ylmo%tpo%now%f_grnd_acy,units="1",long_name="Grounded fraction (acy)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"taub_acx",ylmo%dyn%now%taub_acx,units="Pa",long_name="Basal stress (x)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"taub_acy",ylmo%dyn%now%taub_acy,units="Pa",long_name="Basal stress (y)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"taud_acx",ylmo%dyn%now%taud_acx,units="Pa",long_name="Driving stress (x)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"taud_acy",ylmo%dyn%now%taud_acy,units="Pa",long_name="Driving stress (y)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
       
!        call nc_write(filename,"ux_s",ylmo%dyn%now%ux_s,units="m/a",long_name="Surface velocity (x)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!        call nc_write(filename,"uy_s",ylmo%dyn%now%uy_s,units="m/a",long_name="Surface velocity (y)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"ux_b",ylmo%dyn%now%ux_b,units="m/a",long_name="Basal sliding velocity (x)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_b",ylmo%dyn%now%uy_b,units="m/a",long_name="Basal sliding velocity (y)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"ux_i_bar",ylmo%dyn%now%ux_i_bar,units="m/a",long_name="Internal shear velocity (x)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_i_bar",ylmo%dyn%now%uy_i_bar,units="m/a",long_name="Internal shear velocity (y)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"ux_bar",ylmo%dyn%now%ux_bar,units="m/a",long_name="Vertically averaged velocity (x)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_bar",ylmo%dyn%now%uy_bar,units="m/a",long_name="Vertically averaged velocity (y)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"beta_acx",ylmo%dyn%now%beta_acx,units="Pa a m-1",long_name="Dragging coefficient (x)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"beta_acy",ylmo%dyn%now%beta_acy,units="Pa a m-1",long_name="Dragging coefficient (y)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D_combined



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

end program yelmox



