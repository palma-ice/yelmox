
program yelmox_rtip

    use nml
    use ncio
    use timer
    use timeout
    use yelmo
    use ice_optimization

    ! External libraries
    use geothermal
    use ismip6
    use fastisostasy
    use marine_shelf
    use sealevel
    use sediments
    use smbpal
    use snapclim

    use hyster

    implicit none

    character(len=256) :: outfldr, file1D, file2D, file2D_small, file2D_wais, file3D
    character(len=256) :: file1D_ismip6, file2D_ismip6
    character(len=256) :: file_restart
    character(len=256) :: file_isos_restart
    character(len=256) :: domain, grid_name
    character(len=512) :: path_par
    character(len=512) :: path_tf_corr
    character(len=512) :: ismip6_path_par
    integer  :: n, m, i1wais, i2wais, j1wais, j2wais, xmargin, ymargin
    real(wp) :: time, time_bp
    real(wp) :: time_elapsed
    real(wp) :: time_wt

    real(wp)            :: restart_every_df
    integer             :: counter_restart
    character(len=16)   :: counter_restart_str
    character(len=512)  :: filename_restarts

    real(sp) :: convert_km3_Gt

    type(yelmo_class)           :: yelmo1
    type(sealevel_class)        :: sealev
    type(snapclim_class)        :: snp1
    type(marshelf_class)        :: mshlf1
    type(smbpal_class)          :: smbpal1
    type(sediments_class)       :: sed1
    type(geothermal_class)      :: gthrm1
    type(isos_class)            :: isos1
    type(ismip6_forcing_class)  :: ismp1
    type(hyster_class)          :: hyst1

    type(timeout_class) :: tm_1D, tm_2D, tm_2Dsm, tm_2Dwais, tm_3D

    ! Model timing
    type(timer_class)  :: tmr
    type(timer_class)  :: tmrs
    character(len=512) :: tmr_file

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
        real(wp) :: dtt

        logical  :: kill_shelves
        logical  :: with_ice_sheet
        character(len=56) :: equil_method

        character(len=56) :: abumip_scenario
        real(wp)          :: abumip_bmb

        character(len=56) :: hyst_scenario
        real(wp)          :: hyst_f_to
        real(wp)          :: hyst_f_ta

        character(len=512) :: ismip6_par_file
        character(len=56)  :: ismip6_expname
        logical            :: ismip6_write_formatted
        real(wp)           :: ismip6_dt_formatted

        real(wp) :: isos_tau_1
        real(wp) :: isos_tau_2
        real(wp) :: isos_sigma

    end type

    type(ctrl_params)     :: ctl
    type(ice_opt_params)  :: opt
    type(ismip6_experiment_class) :: ismip6exp


    write(*,*) "Defining margins..."

    ! Read the margin params for cropped output
    ! call nml_read(path_par,"tm_2Dwais","xmargin", xmargin)
    ! call nml_read(path_par,"tm_2Dwais","ymargin", ymargin)
    xmargin = 0
    ymargin = -6

    write(*,*) "Loading command line args..."

    ! Determine the parameter file from the command line
    call yelmo_load_command_line_args(path_par)

    write(*,*) "Loading parameter file..."

    ! Control parameters
    call nml_read(path_par,"ctrl","run_step",   ctl%run_step)

    ! ISMIP6 parameters
    call nml_read(path_par,"ismip6","par_file",         ctl%ismip6_par_file)
    call nml_read(path_par,"ismip6","expname",          ctl%ismip6_expname)
    call nml_read(path_par,"ismip6","write_formatted",  ctl%ismip6_write_formatted)
    call nml_read(path_par,"ismip6","dt_formatted",     ctl%ismip6_dt_formatted)

    if (index(ctl%ismip6_par_file,"ant") .gt. 0) then
        ! Running Antarctica domain, load Antarctica specific parameters
        call ismip6_experiment_def(ismip6exp,ctl%ismip6_expname,ctl%ismip6_par_file,"UCM","YELMO")
    end if

    ! Read run_step specific control parameters
    call nml_read(path_par,trim(ctl%run_step),"time_init",  ctl%time_init)      ! [yr] Starting time
    call nml_read(path_par,trim(ctl%run_step),"time_end",   ctl%time_end)       ! [yr] Ending time
    call nml_read(path_par,trim(ctl%run_step),"dtt",        ctl%dtt)            ! [yr] Main loop time step
    call nml_read(path_par,trim(ctl%run_step),"time_equil", ctl%time_equil)     ! [yr] Years to equilibrate first
    call nml_read(path_par,trim(ctl%run_step),"time_const", ctl%time_const)

    call nml_read(path_par,trim(ctl%run_step),"kill_shelves",ctl%kill_shelves)      ! Kill shelves beyond pd?
    call nml_read(path_par,trim(ctl%run_step),"with_ice_sheet",ctl%with_ice_sheet)  ! Active ice sheet?
    call nml_read(path_par,trim(ctl%run_step),"equil_method",  ctl%equil_method)    ! What method should be used for spin-up?

    if (trim(ctl%equil_method) .eq. "opt") then
        ! Load optimization parameters

        call optimize_par_load(opt,path_par,"opt")

    end if

    write(*,*) "Loading special scenario parameters..."

    ! Check if special scenario is being run, make parameter adjustments
    select case(trim(ctl%run_step))

        case("abumip")

            call nml_read(path_par,"abumip","scenario",          ctl%abumip_scenario)
            call nml_read(path_par,"abumip","bmb",               ctl%abumip_bmb)

        case("hysteresis")

            call nml_read(path_par,"hysteresis","scenario",      ctl%hyst_scenario)
            call nml_read(path_par,"hysteresis","f_to",          ctl%hyst_f_to)
            call nml_read(path_par,"hysteresis","f_ta",          ctl%hyst_f_ta)

    end select

    write(*,*) "Loading output times..."

    ! Get output times
    call timeout_init(tm_1D, path_par, "tm_1D", "small", ctl%time_init, ctl%time_end)
    call timeout_init(tm_2Dsm,path_par,"tm_2Dsm", "medium", ctl%time_init,ctl%time_end)
    call timeout_init(tm_2Dwais, path_par, "tm_2Dwais", "medium", ctl%time_init, ctl%time_end)
    call timeout_init(tm_2D, path_par, "tm_2D", "heavy", ctl%time_init, ctl%time_end)
    call timeout_init(tm_3D, path_par, "tm_3D", "heavy", ctl%time_init, ctl%time_end)

    ! Assume program is running from the output folder
    outfldr = "./"

    write(*,*) "Setting output file names..."

    ! Define input and output locations
    file1D              = trim(outfldr)//"yelmo1D.nc"
    file2D              = trim(outfldr)//"yelmo2D.nc"
    file2D_small        = trim(outfldr)//"yelmo2Dsm.nc"
    file2D_wais         = trim(outfldr)//"yelmo2Dwais.nc"
    file3D              = trim(outfldr)//"yelmo3D.nc"
    file_restart        = trim(outfldr)//"yelmo_restart.nc"
    file_isos_restart   = trim(outfldr)//"isos_restart.nc"
    file1D_ismip6       = trim(outfldr)//"yelmo1D_ismip6.nc"
    file2D_ismip6       = trim(outfldr)//"yelmo2D_ismip6.nc"

    tmr_file            = trim(outfldr)//"timer_table.txt"

    ! Set initial model time
    time    = ctl%time_init
    time_bp = time - 1950.0_wp

    ! Init restart file tracking
    restart_every_df = 0.5
    counter_restart = 0

    !  =========================================================
    ! Print summary of run settings
    write(*,*)
    write(*,*) "run_step:  ", trim(ctl%run_step)
    write(*,*)
    write(*,*) "time_init: ",   ctl%time_init
    write(*,*) "time_end:  ",   ctl%time_end
    write(*,*) "dtt:       ",   ctl%dtt
    write(*,*)

    write(*,*) "ismip6_par_file:        ", trim(ctl%ismip6_par_file)
    write(*,*) "ismip6_expname:         ", trim(ctl%ismip6_expname)
    write(*,*) "ismip6_experiment:      ", trim(ismip6exp%experiment)

    select case(trim(ctl%run_step))

        case("spinup")

            write(*,*) "time_equil: ",    ctl%time_equil
            write(*,*) "time_const: ",    ctl%time_const

            time_bp = ctl%time_const - 1950.0_wp

        case("transient")

            write(*,*) "ismip6_shlf_collapse:   ", ismip6exp%shlf_collapse
            write(*,*) "ismip6_write_formatted: ", ctl%ismip6_write_formatted
            write(*,*) "ismip6_file_suffix:     ", trim(ismip6exp%file_suffix)

        case("abumip")

            write(*,*) "abumip_scenario: ", trim(ctl%abumip_scenario)

    end select

    write(*,*) "time    = ", time
    write(*,*) "time_bp = ", time_bp
    write(*,*)
    write(*,*) "with_ice_sheet: ",  ctl%with_ice_sheet
    write(*,*) "equil_method:   ",  trim(ctl%equil_method)

    ! Start timing
    call timer_step(tmr,comp=-1)

    ! === Initialize ice sheet model =====

    ! Initialize data objects and load initial topography
    call yelmo_init(yelmo1,filename=path_par,grid_def="file",time=time)
    ! yelmo1%par%restart

    ! Store domain and grid_name as shortcuts
    domain    = yelmo1%par%domain
    grid_name = yelmo1%par%grid_name

    ! Ensure optimization fields are allocated and preassigned
    allocate(opt%cf_min(yelmo1%grd%nx,yelmo1%grd%ny))
    allocate(opt%cf_max(yelmo1%grd%nx,yelmo1%grd%ny))

    opt%cf_min = opt%cf_min_par
    opt%cf_max = yelmo1%dyn%par%till_cf_ref

    ! Define specific regions of interest =====================

    select case(trim(domain))

        case("Antarctica")

            ! Define base regions for whole domain first
            regions_mask_fnm = "ice_data/Antarctica/"//trim(yelmo1%par%grid_name)//&
                                "/"//trim(yelmo1%par%grid_name)//"_BASINS-nasa.nc"
            allocate(regions_mask(yelmo1%grd%nx,yelmo1%grd%ny))

            ! Load mask from file
            call nc_read(regions_mask_fnm,"mask_regions",regions_mask)

            ! ajr (2023-03-13): files are now consistent, this fix should not be needed!
            ! ! ajr: fix mask inconsistency at 16km resolution
            ! ! Note: the files themselves should be fixed and made consistent!
            ! if (trim(yelmo1%par%grid_name) .eq. "ANT-16KM") then
            !     where(abs(regions_mask - 4.0) .lt. 1e-3) regions_mask = 1.0
            !     where(abs(regions_mask - 5.0) .lt. 1e-3) regions_mask = 2.0
            ! end if

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

            i1wais = 25
            j1wais = 50
            i2wais = i1wais + 63
            j2wais = j1wais + 63
            ! call get_cropwais_indices(regions_mask, 1.0, i1wais, i2wais, j1wais, j2wais, xmargin, ymargin)

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

            ! Adjust optimization cf_max field specifically for the WAIS and Amery
            !where(abs(yelmo1%bnd%basins - 14.0) .lt. 1e-3) opt%cf_max = min(opt%cf_ref_wais,opt%cf_max)
            !where(abs(yelmo1%bnd%basins -  6.0) .lt. 1e-3) opt%cf_max = min(opt%cf_ref_wais,opt%cf_max)
            ! Note ajr (2023-03-13): this is no longer needed as long as ytill.cf_ref is set to 0.1 instead of 1.0.
            ! Also, it was not nice having cf_ref_wais inside of the opt object in ice_optimization.f90. This
            ! will be removed. If it is needed in the future, cf_ref_wais should be included in a parameter
            ! section loaded in this program directly.

        case DEFAULT

            reg1%write = .FALSE.
            reg2%write = .FALSE.
            reg3%write = .FALSE.

    end select


    ! === Initialize external models (forcing for ice sheet) ======


    ! Initialize global sea level model (bnd%z_sl)
    call sealevel_init(sealev,path_par)

    ! Initialize bedrock model
    call isos_init(isos1, path_par, "isos", yelmo1%grd%nx, yelmo1%grd%ny, &
        yelmo1%grd%dx, yelmo1%grd%dy)

    ! ajr: for now, spatially variable tau is disabled, since it is not clear how to
    ! pass the information from an isos1%output field back to the correlary extended
    ! isos1%domain field.

    ! if (trim(domain) .eq. "Antarctica") then
    !     ! Redefine tau (asthenosphere relaxation constant) as spatially
    !     ! variable field using region mask loaded above (0=deepocean,1=wais,2=eais,3=apis)
    !     call nml_read(path_par,"isos_ant","tau",      ctl%isos_tau_1)
    !     call nml_read(path_par,"isos_ant","tau_eais", ctl%isos_tau_2)
    !     call nml_read(path_par,"isos_ant","sigma",    ctl%isos_sigma)

    !     call isos_set_field(isos1%output%tau, &
    !             [ctl%isos_tau_1,ctl%isos_tau_1,ctl%isos_tau_2,ctl%isos_tau_1], &
    !             [        0.0_wp,        1.0_wp,        2.0_wp,        3.0_wp], &
    !                                   regions_mask,yelmo1%grd%dx,ctl%isos_sigma)

    ! end if

    ! Initialize "climate" model (climate and ocean forcing)
    call snapclim_init(snp1,path_par,domain,yelmo1%par%grid_name,yelmo1%grd%nx,yelmo1%grd%ny,yelmo1%bnd%basins)

    ! Initialize surface mass balance model (bnd%smb, bnd%T_srf)
    call smbpal_init(smbpal1,path_par,x=yelmo1%grd%xc,y=yelmo1%grd%yc,lats=yelmo1%grd%lat)

    ! Initialize marine melt model (bnd%bmb_shlf)
    call marshelf_init(mshlf1,path_par,"marine_shelf",yelmo1%grd%nx,yelmo1%grd%ny,domain,grid_name,yelmo1%bnd%regions,yelmo1%bnd%basins)

    ! Initialize variables inside of ismip6 object
    ismip6_path_par = trim(outfldr)//"/"//trim(ctl%ismip6_par_file)
    call ismip6_forcing_init(ismp1,ismip6_path_par,domain,grid_name,experiment=ismip6exp%experiment, &
                                                                    shlf_collapse=ismip6exp%shlf_collapse)

    ! ===== tf_corr initialization ======

    ! Make sure that tf is prescribed externally
    mshlf1%par%tf_method = 0

    if (yelmo1%par%use_restart) then
        ! Load tf_corr field from file

        call load_tf_corr_from_restart(mshlf1%now%tf_corr,yelmo1%par%restart, &
                                                yelmo1%par%domain,yelmo1%par%grid_name)

    else
        ! Initialize tf_corr to be equal to tf_corr_basin, and
        ! set tf_corr_basin to zero (all corrections will be contained in one field)

        mshlf1%now%tf_corr       = mshlf1%now%tf_corr_basin
        mshlf1%now%tf_corr_basin = 0.0_wp

    end if

    ! === Update external modules and pass variables to yelmo boundaries =======

    ! Sediments
    call sediments_init(sed1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,grid_name)
    yelmo1%bnd%H_sed = sed1%now%H

    ! Geothermal heat flow
    call geothermal_init(gthrm1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,grid_name)
    yelmo1%bnd%Q_geo = gthrm1%now%ghf

    ! Barystatic sea level
    call sealevel_update(sealev,year_bp=time_bp)
    yelmo1%bnd%z_sl  = sealev%z_sl

    ! Initialize the isostasy reference state using reference topography fields
    call isos_init_ref(isos1,yelmo1%bnd%z_bed_ref, yelmo1%bnd%H_ice_ref)

    ! Initialize isostasy using current topography
    ! Optionally pass bsl (scalar) and dz_ss (2D sea-surface perturbation) too
    call isos_init_state(isos1, yelmo1%bnd%z_bed, yelmo1%tpo%now%H_ice, time, bsl=sealev%z_sl)
    
        ! (isos, z_bed, H_ice, time, w, we, &
        ! z_bed_ref, H_ice_ref, bsl, bsl_ref, z_ss_perturb, z_ss_perturb_ref, &
        ! w_ref, we_ref)

    yelmo1%bnd%z_bed = isos1%output%z_bed
    yelmo1%bnd%z_sl  = isos1%output%z_ss

    ! Update snapclim
    call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=time_bp,domain=domain,dx=yelmo1%grd%dx,basins=yelmo1%bnd%basins)

    ! Equilibrate snowpack for itm
    if (trim(smbpal1%par%abl_method) .eq. "itm") then
        call smbpal_update_monthly_equil(smbpal1,snp1%now%tas,snp1%now%pr, &
                               yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time_bp,time_equil=100.0)
    end if

    ! Update forcing to present-day reference using ISMIP6 forcing
    call calc_climate_ismip6(snp1,smbpal1,mshlf1,ismp1,yelmo1, &
                time=ctl%time_const,time_bp=ctl%time_const-1950.0_wp)

    yelmo1%bnd%smb      = smbpal1%ann%smb*yelmo1%bnd%c%conv_we_ie*1e-3   ! [mm we/a] => [m ie/a]
    yelmo1%bnd%T_srf    = smbpal1%ann%tsrf

    yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf
    yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf

    call yelmo_print_bound(yelmo1%bnd)

    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base)
    call yelmo_init_state(yelmo1,time=time,thrm_method="robin-cold")

! ================= RUN STEPS ===============================================


    select case(trim(ctl%run_step))

    case("spinup")
        ! Model can start from no spinup or equilibration (using restart file),
        ! here it is run under constant boundary conditions to spinup
        ! desired state.

        write(*,*)
        write(*,*) "Performing spinup."
        write(*,*)

        ! ===== basal friction optimization ======
        if (trim(ctl%equil_method) .eq. "opt") then

            ! Ensure that cb_ref will be optimized (till_method == set externally)
            yelmo1%dyn%par%till_method = -1

            ! If not using restart...
            if (.not. yelmo1%par%use_restart) then

                if (opt%cf_init .gt. 0.0) then
                    ! Prescribe cb_ref to initial guess
                    yelmo1%dyn%now%cb_ref = opt%cf_init
                else
                    ! Load cb_ref from calculated cb_tgt field
                    yelmo1%dyn%now%cb_ref = yelmo1%dyn%now%cb_tgt

                end if

            end if

        end if

        if (ctl%with_ice_sheet .and. .not. yelmo1%par%use_restart) then
            ! Run yelmo alone for one or a few years with constant boundary conditions
            ! to sort out inconsistencies from initialization.
            call yelmo_update_equil(yelmo1,time,time_tot=1.0_wp,dt=1.0_wp,topo_fixed=.FALSE.)
        end if

        if (trim(ctl%equil_method) .eq. "opt") then
            ! Additional initialization option when running 'opt' spinup...

            if (ctl%with_ice_sheet .and. ctl%time_equil .gt. 0.0) then
                ! Calculate thermodynamics with fixed ice sheet
                call yelmo_update_equil(yelmo1,time,time_tot=ctl%time_equil,dt=ctl%dtt,topo_fixed=.TRUE.)
            end if

        end if

        write(*,*) "Initialization complete."

        ! Initialize output files for checking progress
        call yelmo_write_reg_init(yelmo1, file1D, time_init=time, units="years", &
            mask=yelmo1%bnd%ice_allowed)

        call yelmo_write_init(yelmo1,file2D,time_init=time,units="years")
        call yelmo_write_init(yelmo1,file2D_small,time_init=time,units="years")
        call yelmo_initnc_cropped(yelmo1, file2D_wais, time, "years", &
                                  i1wais, i2wais, j1wais, j2wais)

        call yelmo_write_init3D(yelmo1, file3D, time_init=time, units="years")

        call timer_step(tmr,comp=1,label="initialization")
        call timer_step(tmrs,comp=-1)

        ! Next perform 'coupled' model simulations for desired time
        do n = 0, ceiling((ctl%time_end-ctl%time_init)/ctl%dtt)

            ! Get current time
            time         = ctl%time_init + n*ctl%dtt
            time_bp      = time - 1950.0_wp
            time_elapsed = time - ctl%time_init

            !!ajr: only update optimized fields if ice sheet is running
            if (ctl%with_ice_sheet) then

            select case(trim(ctl%equil_method))

                case("opt")

                    if (time_elapsed .le. opt%rel_time2) then
                        ! Apply relaxation to the model

                        ! Update model relaxation time scale and error scaling (in [m])
                        call optimize_set_transient_param(opt%rel_tau,time_elapsed,time1=opt%rel_time1, &
                                                time2=opt%rel_time2,p1=opt%rel_tau1,p2=opt%rel_tau2,m=opt%rel_m)

                        ! Set model tau, and set yelmo relaxation switch (4: gl line and grounding zone relaxing; 0: no relaxation)
                        yelmo1%tpo%par%topo_rel_tau = opt%rel_tau
                        yelmo1%tpo%par%topo_rel     = 3

                    else
                        ! Turn-off relaxation now

                        yelmo1%tpo%par%topo_rel = 0

                    end if

                    ! === Optimization update step =========

                    if (opt%opt_cf .and. &
                        (time_elapsed .ge. opt%cf_time_init .and. time_elapsed .le. opt%cf_time_end) ) then
                        ! Perform cf_ref optimization

                        ! Update cb_ref based on error metric(s)
                        call optimize_cb_ref(yelmo1%dyn%now%cb_ref,yelmo1%tpo%now%H_ice, &
                                            yelmo1%tpo%now%dHidt,yelmo1%bnd%z_bed,yelmo1%bnd%z_sl,yelmo1%dyn%now%ux_s,yelmo1%dyn%now%uy_s, &
                                            yelmo1%dta%pd%H_ice,yelmo1%dta%pd%uxy_s,yelmo1%dta%pd%H_grnd, &
                                            opt%cf_min,opt%cf_max,yelmo1%tpo%par%dx,opt%sigma_err,opt%sigma_vel,opt%tau_c,opt%H0, &
                                            dt=ctl%dtt,fill_method=opt%fill_method,fill_dist=opt%sigma_err, &
                                            cb_tgt=yelmo1%dyn%now%cb_tgt)

                    end if

                    if (opt%opt_tf .and. &
                        (time_elapsed .ge. opt%tf_time_init .and. time_elapsed .le. opt%tf_time_end) ) then
                        ! Perform tf_corr optimization

                        call optimize_tf_corr(mshlf1%now%tf_corr,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%H_grnd,yelmo1%tpo%now%dHidt, &
                                                yelmo1%dta%pd%H_ice,yelmo1%dta%pd%H_grnd,opt%H_grnd_lim,opt%tau_m,opt%m_temp, &
                                                opt%tf_min,opt%tf_max,yelmo1%tpo%par%dx,sigma=opt%tf_sigma,dt=ctl%dtt)

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

            end if

            ! ====================================================

            call timer_step(tmrs,comp=0)

            ! == SEA LEVEL (BARYSTATIC) ======================================================
            call sealevel_update(sealev,year_bp=time_bp)

            ! == ISOSTASY and SEA LEVEL (REGIONAL) ===========================================
            call isos_update(isos1, yelmo1%tpo%now%H_ice, time, &
                dwdt_corr=yelmo1%bnd%dzbdt_corr)
            yelmo1%bnd%z_bed = isos1%output%z_bed
            yelmo1%bnd%z_sl  = isos1%output%z_ss

            call timer_step(tmrs,comp=1,time_mod=[time-ctl%dtt,time]*1e-3,label="isostasy")

            ! == ICE SHEET ===================================================
            if (ctl%with_ice_sheet) call yelmo_update(yelmo1,time)
            if (ctl%kill_shelves) call maskkill_shelves(yelmo1%tpo%now%H_ice, yelmo1%dta%pd%mask)

            call timer_step(tmrs,comp=2,time_mod=[time-ctl%dtt,time]*1e-3,label="yelmo")

            ! == CLIMATE ===========================================================

            ! Update forcing to present-day reference, but
            ! adjusting to ice topography
            call calc_climate_ismip6(snp1,smbpal1,mshlf1,ismp1,yelmo1, &
                        time=ctl%time_const,time_bp=ctl%time_const-1950.0_wp)

            yelmo1%bnd%smb      = smbpal1%ann%smb*yelmo1%bnd%c%conv_we_ie*1e-3   ! [mm we/a] => [m ie/a]
            yelmo1%bnd%T_srf    = smbpal1%ann%tsrf

            yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf
            yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf

            call timer_step(tmrs,comp=3,time_mod=[time-ctl%dtt,time]*1e-3,label="climate")

            ! == MODEL OUTPUT ===================================

            if (timeout_check(tm_1D, time)) then
                call yelmo_write_reg_step(yelmo1, file1D, time=time)
            end if

            if (timeout_check(tm_2D, time)) then
                call write_step_2D_combined(yelmo1, isos1, snp1, mshlf1, smbpal1, file2D, time)
            end if

            if (timeout_check(tm_2Dsm, time)) then
                call yx_hyst_write_step_2D_combined_small(yelmo1, isos1, snp1, &
                    mshlf1,smbpal1,file2D_small,time)
            end if

            if ( timeout_check(tm_2Dwais,time) ) then
                call yx_hyst_write_step_2Dwais(yelmo1, mshlf1, file2D_wais, time, &
                    i1wais, i2wais, j1wais, j2wais)
            end if

            if (timeout_check(tm_3D, time)) then
                call write_step_3D(yelmo1, file3D, time)
            end if


            call timer_step(tmrs,comp=4,time_mod=[time-ctl%dtt,time]*1e-3,label="io")

            if (mod(time_elapsed,10.0)==0) then
                ! Print timestep timing info and write log table
                call timer_write_table(tmrs,[time,ctl%dtt]*1e-3,"m",tmr_file,init=time_elapsed .eq. 0.0)
            end if

            if (mod(time_elapsed,10.0)==0 .and. (.not. yelmo_log)) then
                write(*,"(a,f14.4)") "yelmo:: time = ", time
            end if

        end do

        write(*,*)
        write(*,*) "spinup complete."
        write(*,*)

        ! Write the restart file for the end of the simulation
        call yelmo_restart_write(yelmo1,file_restart,time=time_bp)
        call isos_restart_write(isos1,file_isos_restart,time_bp)

    case("transient")
        ! Here it is assumed that the model has gone through spinup
        ! and is ready for transient simulations

        write(*,*)
        write(*,*) "Performing transient."
        write(*,*)

        ! Get current time
        time    = ctl%time_init
        time_bp = time - 1950.0_wp

        ! Initialize output files
        call yelmo_write_init(yelmo1,file2D,time_init=time,units="years")
        call yelmo_write_reg_init(yelmo1,file1D,time_init=time,units="years",mask=yelmo1%bnd%ice_allowed)

        if (ctl%ismip6_write_formatted) then
            ! Initialize output files for ISMIP6
            call yelmo_write_init(yelmo1,file2D_ismip6,time_init=time,units="years")
            call yelmo_write_reg_init(yelmo1,file1D_ismip6,time_init=time,units="years",mask=yelmo1%bnd%ice_allowed)
        end if

        call timer_step(tmr,comp=1,label="initialization")
        call timer_step(tmrs,comp=-1)

        ! Perform 'coupled' model simulations for desired time
        do n = 0, ceiling((ctl%time_end-ctl%time_init)/ctl%dtt)

            ! Get current time

            time         = ctl%time_init + n*ctl%dtt
            time_bp      = time - 1950.0_wp
            time_elapsed = time - ctl%time_init

            if (ismip6exp%shlf_collapse) then
            ! Perform mask_shlf_collapse experiments
            ! Set H to zero where mask==1, then compute Yelmo.

                if(time .ge. 2015) then
                    !where((yelmo1%tpo%now%f_grnd .eq. 0.0) .and. (ismp1%mask_shlf%var(:,:,1,1) .eq. 1.0)) yelmo1%tpo%now%H_ice = 0.0
                    where((yelmo1%tpo%now%f_grnd .eq. 0.0) .and. (ismp1%mask_shlf%var(:,:,1,1) .eq. 1.0)) yelmo1%bnd%ice_allowed = .FALSE.
                end if
            end if

            call timer_step(tmrs,comp=0)

            ! == SEA LEVEL (BARYSTATIC) ======================================================
            call sealevel_update(sealev,year_bp=0.0_wp)

            ! == ISOSTASY and SEA LEVEL (REGIONAL) ===========================================
            call isos_update(isos1, yelmo1%tpo%now%H_ice, sealev%z_sl, time, &
                                                        dwdt_corr=yelmo1%bnd%dzbdt_corr)
            yelmo1%bnd%z_bed = isos1%output%z_bed
            yelmo1%bnd%z_sl  = isos1%output%z_ss

            call timer_step(tmrs,comp=1,time_mod=[time-ctl%dtt,time]*1e-3,label="isostasy")

            ! == ICE SHEET ===================================================
            if (ctl%with_ice_sheet) call yelmo_update(yelmo1,time)
            if (ctl%kill_shelves) call maskkill_shelves(yelmo1%tpo%now%H_ice, yelmo1%dta%pd%mask)

            if (ismip6exp%shlf_collapse) then
                        ! Clean up icebergs for mask_shlf_collapse experiments
                        call calc_iceberg_island(ismp1%iceberg_mask,yelmo1%tpo%now%f_grnd, &
                            yelmo1%tpo%now%H_ice)
                        where(ismp1%iceberg_mask .eq. 1.0) yelmo1%tpo%now%H_ice = 0.0
            end if

            call timer_step(tmrs,comp=2,time_mod=[time-ctl%dtt,time]*1e-3,label="yelmo")

            ! == CLIMATE and OCEAN ==========================================

            ! Get ISMIP6 climate and ocean forcing
            call calc_climate_ismip6(snp1,smbpal1,mshlf1,ismp1,yelmo1,time,time_bp)

            yelmo1%bnd%smb      = smbpal1%ann%smb*yelmo1%bnd%c%conv_we_ie*1e-3   ! [mm we/a] => [m ie/a]
            yelmo1%bnd%T_srf    = smbpal1%ann%tsrf

            yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf
            yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf

            call timer_step(tmrs,comp=3,time_mod=[time-ctl%dtt,time]*1e-3,label="climate")

            ! == MODEL OUTPUT ===================================

            if (timeout_check(tm_1D, time)) then
                call yelmo_write_reg_step(yelmo1, file1D, time=time)
            end if

            if (timeout_check(tm_2D, time)) then
                call write_step_2D_combined(yelmo1, isos1, snp1, mshlf1, smbpal1, file2D, time)
            end if

            if (timeout_check(tm_2Dsm,time)) then
                call yx_hyst_write_step_2D_combined_small(yelmo1, isos1, snp1, &
                    mshlf1,smbpal1,file2D_small,time)
            end if

            if ( timeout_check(tm_2Dwais,time) ) then
                call yx_hyst_write_step_2Dwais(yelmo1, mshlf1, file2D_wais, time, &
                    i1wais, i2wais, j1wais, j2wais)
            end if

            if (timeout_check(tm_3D, time)) then
                call write_step_3D(yelmo1, file3D, time)
            end if

            ! ISMIP6 output if desired:
            if (ctl%ismip6_write_formatted) then
                if (mod(nint(time_elapsed*100),nint(ctl%ismip6_dt_formatted*100))==0) then
                    call write_step_2D_ismip6(yelmo1,file2D_ismip6,time)
                    call write_1D_ismip6(yelmo1,file1D_ismip6,time)
                end if
            end if

            call timer_step(tmrs,comp=4,time_mod=[time-ctl%dtt,time]*1e-3,label="io")

            if (mod(time_elapsed,10.0)==0) then
                ! Print timestep timing info and write log table
                call timer_write_table(tmrs,[time,ctl%dtt]*1e-3,"m",tmr_file,init=time_elapsed .eq. 0.0)
            end if

            if (mod(time_elapsed,10.0)==0 .and. (.not. yelmo_log)) then
                write(*,"(a,f14.4)") "yelmo:: time = ", time
            end if

        end do

        write(*,*)
        write(*,*) "Transient complete."
        write(*,*)

        ! Write the restart file for the end of the transient simulation
        call yelmo_restart_write(yelmo1,file_restart,time=time)
        call isos_restart_write(isos1,file_isos_restart,time)

    case("abumip")
        ! Here it is assumed that the model has gone through spinup
        ! and is ready for transient simulations

        write(*,*)
        write(*,*) "Performing transient. [abumip]"
        write(*,*)

        ! Get current time
        time    = ctl%time_init
        time_bp = time - 1950.0_wp

        ! Initialize output files
        call yelmo_write_init(yelmo1,file2D,time_init=time,units="years")
        call yelmo_write_reg_init(yelmo1,file1D,time_init=time,units="years",mask=yelmo1%bnd%ice_allowed)

        call timer_step(tmr,comp=1,label="initialization")
        call timer_step(tmrs,comp=-1)

        ! Perform 'coupled' model simulations for desired time
        do n = 0, ceiling((ctl%time_end-ctl%time_init)/ctl%dtt)

            ! Get current time
            time         = ctl%time_init + n*ctl%dtt
            time_bp      = time - 1950.0_wp
            time_elapsed = time - ctl%time_init

            ! == ABUMIP =========================================================

            ! Make parameter changes relevant to abumip

            select case(trim(ctl%abumip_scenario))

                case("abuc")

                    ! Do nothing - control experiment

                case("abuk")
                    ! Ensure ice shelves are killed

                    yelmo1%tpo%par%calv_flt_method = "kill"

                case("abum")
                    ! Apply 400 m/yr melt rate on shelves

                    yelmo1%bnd%bmb_shlf = ctl%abumip_bmb     ! [m/yr]

                case DEFAULT

                    write(io_unit_err,*) ""
                    write(io_unit_err,*) "yelmox_ismip6:: error: abumip scenario not recognized."
                    write(io_unit_err,*) "abumip_scenario: ", trim(ctl%abumip_scenario)
                    stop 1

            end select

            call timer_step(tmrs,comp=0)

            ! == SEA LEVEL (BARYSTATIC) ======================================================
            call sealevel_update(sealev,year_bp=0.0_wp)

            ! == ISOSTASY and SEA LEVEL (REGIONAL) ===========================================
            call isos_update(isos1, yelmo1%tpo%now%H_ice, sealev%z_sl, time, &
                                                        dwdt_corr=yelmo1%bnd%dzbdt_corr)

            yelmo1%bnd%z_bed = isos1%output%z_bed
            yelmo1%bnd%z_sl  = isos1%output%z_ss

            call timer_step(tmrs,comp=1,time_mod=[time-ctl%dtt,time]*1e-3,label="isostasy")

            ! == ICE SHEET ===================================================
            if (ctl%with_ice_sheet) call yelmo_update(yelmo1,time)
            if (ctl%kill_shelves) call maskkill_shelves(yelmo1%tpo%now%H_ice, yelmo1%dta%pd%mask)

            call timer_step(tmrs,comp=2,time_mod=[time-ctl%dtt,time]*1e-3,label="yelmo")

            ! ISMIP6 forcing

            ! Update ismip6 forcing to current time
            call ismip6_forcing_update(ismp1,ctl%time_const)

            ! Set climate to present day
            snp1%now = snp1%clim0

            ! == SURFACE MASS BALANCE ==============================================

if (n .eq. 0) then
                ! Calculate smb for present day
                call smbpal_update_monthly(smbpal1,snp1%now%tas,snp1%now%pr, &
                                           yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,ctl%time_const)

                ! Apply ISMIP6 anomalies
                ! (apply to climate just for consistency)

                smbpal1%ann%smb  = smbpal1%ann%smb  + ismp1%smb%var(:,:,1,1)*1.0/(yelmo1%bnd%c%conv_we_ie*1e-3) ! [m ie/yr] => [mm we/a]
                smbpal1%ann%tsrf = smbpal1%ann%tsrf + ismp1%ts%var(:,:,1,1)

                do m = 1,12
                    snp1%now%tas(:,:,m) = snp1%now%tas(:,:,m) + ismp1%ts%var(:,:,1,1)
                    snp1%now%pr(:,:,m)  = snp1%now%pr(:,:,m)  + ismp1%pr%var(:,:,1,1)/365.0 ! [mm/yr] => [mm/d]
                end do

                snp1%now%ta_ann = sum(snp1%now%tas,dim=3) / 12.0_wp
                if (trim(domain) .eq. "Antarctica") then
                    snp1%now%ta_sum  = sum(snp1%now%tas(:,:,[12,1,2]),dim=3)/3.0  ! Antarctica summer
                else
                    snp1%now%ta_sum  = sum(snp1%now%tas(:,:,[6,7,8]),dim=3)/3.0  ! NH summer
                end if
                snp1%now%pr_ann = sum(snp1%now%pr,dim=3)  / 12.0 * 365.0     ! [mm/d] => [mm/a]

end if

            ! == MARINE AND TOTAL BASAL MASS BALANCE ===============================

            call marshelf_update_shelf(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                            yelmo1%bnd%basins,yelmo1%bnd%z_sl,yelmo1%grd%dx,-ismp1%to%z, &
                            ismp1%to%var(:,:,:,1),ismp1%so%var(:,:,:,1), &
                            dto_ann=ismp1%to%var(:,:,:,1)-ismp1%to_ref%var(:,:,:,1), &
                            tf_ann=ismp1%tf%var(:,:,:,1))

            ! Update temperature forcing field with tf_corr and tf_corr_basin
            mshlf1%now%tf_shlf = mshlf1%now%tf_shlf + mshlf1%now%tf_corr + mshlf1%now%tf_corr_basin

            call marshelf_update(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                                    yelmo1%bnd%regions,yelmo1%bnd%basins,yelmo1%bnd%z_sl,dx=yelmo1%grd%dx)

            yelmo1%bnd%smb      = smbpal1%ann%smb*yelmo1%bnd%c%conv_we_ie*1e-3
            yelmo1%bnd%T_srf    = smbpal1%ann%tsrf

            yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf
            yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf

            if (trim(ctl%abumip_scenario) .eq. "abum") then
                ! Ensure bmb_shlf output is consistent with what is applied

                yelmo1%bnd%bmb_shlf = ctl%abumip_bmb     ! [m/yr]

            end if

            call timer_step(tmrs,comp=3,time_mod=[time-ctl%dtt,time]*1e-3,label="climate")

            ! == MODEL OUTPUT ===================================

            if (timeout_check(tm_1D,time)) then
                call yx_hyst_write_step_1D_combined(yelmo1,hyst1,snp1,file1D,time=time)

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

            if (timeout_check(tm_2D, time)) then
                call write_step_2D_combined(yelmo1, isos1, snp1, mshlf1, smbpal1, file2D, time)
            end if

            if (timeout_check(tm_2Dsm,time)) then
                call yx_hyst_write_step_2D_combined_small(yelmo1, isos1, snp1, &
                    mshlf1,smbpal1,file2D_small,time)
            end if

            if ( timeout_check(tm_2Dwais,time) ) then
                call yx_hyst_write_step_2Dwais(yelmo1, mshlf1, file2D_wais, time, &
                    i1wais, i2wais, j1wais, j2wais)
            end if

            if (timeout_check(tm_3D, time)) then
                call write_step_3D(yelmo1, file3D, time)
            end if

            call timer_step(tmrs,comp=4,time_mod=[time-ctl%dtt,time]*1e-3,label="io")

            if (mod(time_elapsed,10.0)==0) then
                ! Print timestep timing info and write log table
                call timer_write_table(tmrs,[time,ctl%dtt]*1e-3,"m",tmr_file,init=time_elapsed .eq. 0.0)
            end if

            if (mod(time_elapsed,10.0)==0 .and. (.not. yelmo_log)) then
                write(*,"(a,f14.4)") "yelmo:: time = ", time
            end if

        end do

        write(*,*)
        write(*,*) "Transient complete."
        write(*,*)

        ! Write the restart file for the end of the transient simulation
        call yelmo_restart_write(yelmo1,file_restart,time=time)
        call isos_restart_write(isos1,file_isos_restart,time)

    case("hysteresis")
        ! Here it is assumed that the model has gone through spinup
        ! and is ready for transient simulations

        write(*,*)
        write(*,*) "Performing transient. [hysteresis]"
        write(*,*)

        ! Get current time
        time    = ctl%time_init
        time_bp = time - 1950.0_wp

        ! === HYST ============

        ! Initialize hysteresis module for transient forcing experiments
        call hyster_init(hyst1,path_par,time)
        convert_km3_Gt = yelmo1%bnd%c%rho_ice *1e-3

        ! =====================

        ! Update forcing to constant reference time with initial hyst forcing
        call calc_climate_ismip6(snp1,smbpal1,mshlf1,ismp1,yelmo1, &
                    time=ctl%time_const,time_bp=ctl%time_const-1950.0_wp, &
                    dTa=hyst1%f_now*ctl%hyst_f_ta,dTo=hyst1%f_now*ctl%hyst_f_to)

        yelmo1%bnd%smb      = smbpal1%ann%smb*yelmo1%bnd%c%conv_we_ie*1e-3   ! [mm we/a] => [m ie/a]
        yelmo1%bnd%T_srf    = smbpal1%ann%tsrf

        yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf
        yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf


        ! Initialize hysteresis output files
        call yx_hyst_write_yelmo_init_1D_combined(yelmo1, file1D, time, units="years", &
            mask=yelmo1%bnd%ice_allowed)

        call yelmo_write_init(yelmo1, file2D, time_init=time, units="years")
        call yelmo_initnc_cropped(yelmo1, file2D_wais, time, "years", &
                                  i1wais, i2wais, j1wais, j2wais)
        call yelmo_write_init(yelmo1, file2D_small, time_init=time, units="years")
        call yelmo_write_init3D(yelmo1, file3D, time_init=time, units="years")

        if (reg1%write) then
            call yelmo_write_reg_init(yelmo1,reg1%fnm,time_init=time,units="years",mask=reg1%mask)
        end if

        if (reg2%write) then
            call yelmo_write_reg_init(yelmo1,reg2%fnm,time_init=time,units="years",mask=reg2%mask)
        end if

        if (reg3%write) then
            call yelmo_write_reg_init(yelmo1,reg3%fnm,time_init=time,units="years",mask=reg3%mask)
        end if

        call timer_step(tmr,comp=1,label="initialization")
        call timer_step(tmrs,comp=-1)

        ! Perform 'coupled' model simulations for desired time
        do n = 0, ceiling((ctl%time_end-ctl%time_init)/ctl%dtt)

            ! Get current time
            time         = ctl%time_init + n*ctl%dtt
            time_bp      = time - 1950.0_wp
            time_elapsed = time - ctl%time_init

            ! == HYSTERESIS =========================================================

            ! Make parameter changes relevant to hysteresis runs

            ! snapclim should use the anomaly method
            snp1%par%atm_type = "anom"
            snp1%par%ocn_type = "anom"

            select case(trim(ctl%hyst_scenario))

                case("ctrl")

                    ! Do nothing - control experiment

                case("scenario1")
                    ! Possible scenario 1

                    ! To do

                case("scenario2")
                    ! Possible scenario 2

                    ! To do

                case DEFAULT

                    write(io_unit_err,*) ""
                    write(io_unit_err,*) "yelmox_ismip6:: error: hysteresis scenario not recognized."
                    write(io_unit_err,*) "hysteresis.scenario: ", trim(ctl%hyst_scenario)
                    stop 1

            end select

            call timer_step(tmrs,comp=0)

            ! == SEA LEVEL (BARYSTATIC) ======================================================
            call sealevel_update(sealev,year_bp=0.0_wp)

            ! == ISOSTASY and SEA LEVEL (REGIONAL) ===========================================
            call isos_update(isos1, yelmo1%tpo%now%H_ice, sealev%z_sl, time, &
                                                        dwdt_corr=yelmo1%bnd%dzbdt_corr)
            yelmo1%bnd%z_bed = isos1%output%z_bed
            yelmo1%bnd%z_sl  = isos1%output%z_ss

            call timer_step(tmrs,comp=1,time_mod=[time-ctl%dtt,time]*1e-3,label="isostasy")

            ! == ICE SHEET ===================================================
            if (ctl%with_ice_sheet) call yelmo_update(yelmo1,time)
            if (ctl%kill_shelves) call maskkill_shelves(yelmo1%tpo%now%H_ice, yelmo1%dta%pd%mask)

            call timer_step(tmrs,comp=2,time_mod=[time-ctl%dtt,time]*1e-3,label="yelmo")


            ! === HYST ============

            ! snapclim call using anomaly from the hyster package
            call hyster_calc_forcing(hyst1,time=time,var=yelmo1%reg%V_ice*convert_km3_Gt)

            ! =====================

            ! == CLIMATE (ATMOSPHERE, OCEAN and SMB) ====================================

            ! Update forcing to initial time with initial hyst forcing
            call calc_climate_ismip6(snp1,smbpal1,mshlf1,ismp1,yelmo1, &
                        time=ctl%time_const,time_bp=ctl%time_const-1950.0_wp, &
                        dTa=hyst1%f_now*ctl%hyst_f_ta,dTo=hyst1%f_now*ctl%hyst_f_to)

            yelmo1%bnd%smb      = smbpal1%ann%smb*yelmo1%bnd%c%conv_we_ie*1e-3   ! [mm we/a] => [m ie/a]
            yelmo1%bnd%T_srf    = smbpal1%ann%tsrf

            yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf
            yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf

            call timer_step(tmrs,comp=3,time_mod=[time-ctl%dtt,time]*1e-3,label="climate")

            ! == MODEL OUTPUT ===================================

            ! ** Using routines from yelmox_hysteresis_help **

            if (timeout_check(tm_2D, time)) then
                write(*,*) "writing 2D at t=", time
                call write_step_2D_combined(yelmo1, isos1, snp1, mshlf1, smbpal1, file2D, time)
            end if


            if (timeout_check(tm_2Dsm,time)) then
                write(*,*) "writing 2D small at t=", time
                call yx_hyst_write_step_2D_combined_small(yelmo1, isos1, snp1, &
                    mshlf1,smbpal1,file2D_small,time)
            end if

            if ( timeout_check(tm_2Dwais,time) ) then
                write(*,*) "writing 2D wais at t=", time
                call yx_hyst_write_step_2Dwais(yelmo1, mshlf1, file2D_wais, time, &
                    i1wais, i2wais, j1wais, j2wais)
            end if

            if (timeout_check(tm_3D, time)) then
                write(*,*) "writing 3D at t=", time
                call write_step_3D(yelmo1, file3D, time)
            end if

            if (timeout_check(tm_1D, time)) then
                write(*,*) "writing 1D at t=", time

                call yx_hyst_write_step_1D_combined(yelmo1,hyst1,snp1,file1D,time=time)

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


            call timer_step(tmrs,comp=4,time_mod=[time-ctl%dtt,time]*1e-3,label="io")

            if (mod(time_elapsed,10.0)==0) then
                ! Print timestep timing info and write log table
                call timer_write_table(tmrs,[time,ctl%dtt]*1e-3,"m",tmr_file,init=time_elapsed .eq. 0.0)
            end if

            if ((hyst1%f_now) .gt. ((counter_restart+1) * restart_every_df)) then
                write (*,*) "writing restart file at f=", hyst1%f_now
                counter_restart = counter_restart + 1
                if (counter_restart .lt. 10) then
                    write (counter_restart_str, '(i1)') counter_restart
                else if (counter_restart .lt. 100) then
                    write (counter_restart_str, '(i2)') counter_restart
                end if
                filename_restarts = trim(outfldr)//"yelmo_restart_"//trim(counter_restart_str)//".nc"
                call yelmo_restart_write(yelmo1, filename_restarts, time=time_bp)
            end if

            if (mod(time_elapsed,10.0)==0 .and. (.not. yelmo_log)) then
                write(*,"(a,f14.4)") "yelmo:: time = ", time
            end if

        end do

        ! Stop timing
        call timer_step(tmr,comp=2,time_mod=[ctl%time_init,time]*1e-3,label="timeloop")

        write(*,*)
        write(*,*) "Transient complete."
        write(*,*)

        ! Write the restart file for the end of the transient simulation
        call yelmo_restart_write(yelmo1,file_restart,time=time)
        call isos_restart_write(isos1,file_isos_restart,time)

    end select

    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Print timing summary
    call timer_print_summary(tmr,units="m",units_mod="kyr",time_mod=time*1e-3)
    
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

        character(len=*),       intent(IN) :: filename
        real(wp),               intent(IN) :: time

        ! Local variables
        integer  :: ncid, n
        real(wp) :: time_prev

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

        ! Write constant fields
        if (n .eq. 1) then
            call nc_write(filename,"isos_tau",isos%output%tau,units="yr",long_name="Asthenospheric relaxation timescale", &
                      dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
        end if

        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_grz",ylmo%tpo%now%mask_grz,units="",long_name="Grounding-zone mask", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ice_allowed",ylmo%bnd%ice_allowed,units="",long_name="Ice allowed mask", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        !call nc_write(filename,"mask_frnt",ylmo%tpo%now%mask_frnt,units="",long_name="Ice-front mask", &
         !             dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        !call nc_write(filename,"dist_grline",ylmo%tpo%now%dist_grline,units="km",long_name="Distance to grounding line", &
         !             dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dHidt",ylmo%tpo%now%dHidt,units="m/yr",long_name="Ice thickness rate of change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"mb_net",ylmo%tpo%now%mb_net,units="m",long_name="Applied net mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! call nc_write(filename,"taul_int_acx",ylmo%dyn%now%taul_int_acx,units="Pa m",long_name="Vertically integrated lateral stress (x)", &
        !             dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"taul_int_acy",ylmo%dyn%now%taul_int_acy,units="Pa m",long_name="Vertically integrated lateral stress (y)", &
        !             dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically-averaged velocity magnitude", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity magnitude", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"duxydt",ylmo%dyn%now%duxydt,units="m/yr^2",long_name="Velocity rate of change", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        !call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
        !              dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        !call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
        !              dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dist_grline",ylmo%tpo%now%dist_grline,units="km",long_name="Distance to grounding line", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dHidt",ylmo%tpo%now%dHidt,units="m/yr",long_name="Ice thickness rate of change", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"mb_net",ylmo%tpo%now%mb_net,units="m",long_name="Applied net mass balance", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"H_grnd",ylmo%tpo%now%H_grnd,units="m",long_name="Ice thickness overburden", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="bar",long_name="Effective pressure", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"cmb",ylmo%tpo%now%cmb,units="m/a ice equiv.",long_name="Calving mass balance rate", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"cmb_flt",ylmo%tpo%now%cmb_flt,units="m/a ice equiv.",long_name="Calving mass balance rate flt", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"f_grnd",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded fraction", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",ylmo%tpo%now%f_ice,units="1",long_name="Ice fraction in grid cell", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"f_grnd_bmb",ylmo%tpo%now%f_grnd_bmb,units="1",long_name="Grounded fraction (bmb)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        !call nc_write(filename,"dist_grline",ylmo%tpo%now%dist_grline,units="km", &
        !              long_name="Distance to nearest grounding-line point", &
        !              dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"cb_ref",ylmo%dyn%now%cb_ref,units="--",long_name="Bed friction scalar", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"c_bed",ylmo%dyn%now%c_bed,units="Pa",long_name="Bed friction coefficient", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m^-1",long_name="Basal friction coefficient", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff_int",ylmo%dyn%now%visc_eff_int,units="Pa a m",long_name="Depth-integrated effective viscosity (SSA)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud",ylmo%dyn%now%taud,units="Pa",long_name="Driving stress", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taub",ylmo%dyn%now%taub,units="Pa",long_name="Basal stress", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        !call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
        !              dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically-averaged velocity magnitude", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity magnitude", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"duxydt",ylmo%dyn%now%duxydt,units="m/yr^2",long_name="Velocity rate of change", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        !call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
        !              dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        !call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
        !              dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"T_prime_b",ylmo%thrm%now%T_prime_b,units="deg C",long_name="Homologous basal ice temperature", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        !call nc_write(filename,"uz",ylmo%dyn%now%uz,units="m/a",long_name="Vertical velocity (z)", &
        !               dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",start=[1,1,1,n],ncid=ncid)

        !call nc_write(filename,"Q_strn",ylmo%thrm%now%Q_strn/(rho_ice*ylmo%thrm%now%cp),units="K a-1",long_name="Strain heating", &
        !              dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="J a-1 m-2",long_name="Basal frictional heating", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"bmb_grnd",ylmo%thrm%now%bmb_grnd,units="m/a ice equiv.",long_name="Basal mass balance (grounded)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_w",ylmo%thrm%now%H_w,units="m",long_name="Basal water layer thickness", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        !call nc_write(filename,"ATT",ylmo%mat%now%ATT,units="a^-1 Pa^-3",long_name="Rate factor", &
        !              dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"f_shear_bar",ylmo%mat%now%f_shear_bar,units="1",long_name="Vertically averaged shearing fraction", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"enh_bar",ylmo%mat%now%enh_bar,units="1",long_name="Vertically averaged enhancement factor", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_int",ylmo%mat%now%visc_int,units="Pa a m",long_name="Vertically integrated viscosity", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Boundaries
        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        !call nc_write(filename,"H_sed",ylmo%bnd%H_sed,units="m",long_name="Sediment thickness", &
        !              dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
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

        call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW/m^2",long_name="Geothermal heat flux", &
                        dim1="xc",dim2="yc",start=[1,1],ncid=ncid)

        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a ice equiv.",long_name="Basal mass balance", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"fmb",ylmo%tpo%now%fmb,units="m/a ice equiv.",long_name="Margin-front mass balance", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! External data
        call nc_write(filename,"dzbdt",isos%output%dwdt,units="m/a",long_name="Bedrock uplift rate", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Comparison with present-day
        call nc_write(filename,"H_ice_pd_err",ylmo%dta%pd%err_H_ice,units="m",long_name="Ice thickness error wrt present day", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf_pd_err",ylmo%dta%pd%err_z_srf,units="m",long_name="Surface elevation error wrt present day", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s_pd_err",ylmo%dta%pd%err_uxy_s,units="m/a",long_name="Surface velocity error wrt present day", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ssa_mask_acx",ylmo%dyn%now%ssa_mask_acx,units="1",long_name="SSA mask (acx)", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_mask_acy",ylmo%dyn%now%ssa_mask_acy,units="1",long_name="SSA mask (acy)", &
                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Ta_ann",snp%now%ta_ann,units="K",long_name="Near-surface air temperature (ann)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Ta_sum",snp%now%ta_sum,units="K",long_name="Near-surface air temperature (sum)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Pr_ann",snp%now%pr_ann*1e-3,units="m/a water equiv.",long_name="Precipitation (ann)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"T_shlf",mshlf%now%T_shlf,units="K",long_name="Shelf temperature", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"S_shlf",mshlf%now%S_shlf,units="PSU",long_name="Shelf salinity", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dT_shlf",mshlf%now%dT_shlf,units="K",long_name="Shelf temperature anomaly", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_fp_shlf",mshlf%now%T_fp_shlf,units="K",long_name="Shelf freezing temperature", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"mask_ocn",mshlf%now%mask_ocn,units="", &
                        long_name="Ocean mask (0: land, 1: grline, 2: fltline, 3: open ocean, 4: deep ocean, 5: lakes)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"tf_basin",mshlf%now%tf_basin,units="K",long_name="Mean basin thermal forcing", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"tf_shlf",mshlf%now%tf_shlf,units="K",long_name="Shelf thermal forcing", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"tf_corr",mshlf%now%tf_corr,units="K",long_name="Shelf thermal forcing correction factor", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"tf_corr_basin",mshlf%now%tf_corr_basin,units="K",long_name="Shelf thermal forcing basin-wide correction factor", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"slope_base",mshlf%now%slope_base,units="",long_name="Shelf-base slope", &
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
        end if

        !call nc_write(filename,"pr",snp%now%pr*1e-3,units="m/a water equiv.",long_name="Precipitation (ann)", &
        !              dim1="xc",dim2="yc",dim3="month",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"dTa_ann",snp%now%ta_ann-snp%clim0%ta_ann,units="K",long_name="Near-surface air temperature anomaly (ann)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dTa_sum",snp%now%ta_sum-snp%clim0%ta_sum,units="K",long_name="Near-surface air temperature anomaly (sum)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dPr_ann",(snp%now%pr_ann-snp%clim0%pr_ann)*1e-3,units="m/a water equiv.",long_name="Precipitation anomaly (ann)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        !call nc_write(filename,"PDDs",srf%ann%PDDs,units="degC days",long_name="Positive degree days (annual total)", &
        !              dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Comparison with present-day
        call nc_write(filename,"H_ice_pd_err",ylmo%dta%pd%err_H_ice,units="m",long_name="Ice thickness error wrt present day", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf_pd_err",ylmo%dta%pd%err_z_srf,units="m",long_name="Surface elevation error wrt present day", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s_pd_err",ylmo%dta%pd%err_uxy_s,units="m/a",long_name="Surface velocity error wrt present day", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dzsdx",ylmo%tpo%now%dzsdx,units="m/m",long_name="Surface slope", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dzsdy",ylmo%tpo%now%dzsdy,units="m/m",long_name="Surface slope", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        !call nc_write(filename,"f_grnd_acx",ylmo%tpo%now%f_grnd_acx,units="1",long_name="Grounded fraction (acx)", &
        !              dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        !call nc_write(filename,"f_grnd_acy",ylmo%tpo%now%f_grnd_acy,units="1",long_name="Grounded fraction (acy)", &
        !              dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        !call nc_write(filename,"taub_acx",ylmo%dyn%now%taub_acx,units="Pa",long_name="Basal stress (x)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        !call nc_write(filename,"taub_acy",ylmo%dyn%now%taub_acy,units="Pa",long_name="Basal stress (y)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        !call nc_write(filename,"taud_acx",ylmo%dyn%now%taud_acx,units="Pa",long_name="Driving stress (x)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        !call nc_write(filename,"taud_acy",ylmo%dyn%now%taud_acy,units="Pa",long_name="Driving stress (y)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ux_s",ylmo%dyn%now%ux_s,units="m/a",long_name="Surface velocity (x)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_s",ylmo%dyn%now%uy_s,units="m/a",long_name="Surface velocity (y)", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Strain-rate and stress tensors
        if (.FALSE.) then

            !call nc_write(filename,"de",ylmo%mat%now%strn%de,units="a^-1",long_name="Effective strain rate", &
            !              dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
            !call nc_write(filename,"te",ylmo%mat%now%strs%te,units="Pa",long_name="Effective stress", &
            !              dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
            !call nc_write(filename,"visc_int",ylmo%mat%now%visc_int,units="Pa a m",long_name="Depth-integrated effective viscosity (SSA)", &
            !              dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

            !call nc_write(filename,"de2D",ylmo%mat%now%strn2D%de,units="yr^-1",long_name="Effective strain rate", &
            !              dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"div2D",ylmo%mat%now%strn2D%div,units="yr^-1",long_name="Divergence strain rate", &
                            dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            !call nc_write(filename,"te2D",ylmo%mat%now%strs2D%te,units="Pa",long_name="Effective stress", &
            !              dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

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

        ! Close the netcdf file
        call nc_close(ncid)

        return

    end subroutine write_step_2D_combined


    subroutine write_step_3D(ylmo, filename, time)

        implicit none

        type(yelmo_class), intent(IN) :: ylmo
        character(len=*), intent(IN) :: filename
        real(wp), intent(IN) :: time

        ! Local variables
        integer  :: ncid, n, nx, ny, nz
        real(wp) :: time_prev

        ! nx = size(ylmo%dyn%now%ux, 1)
        ! ny = size(ylmo%dyn%now%ux, 2)
        ! nz = size(ylmo%dyn%now%ux, 3)
        ! write (*, *) "nx, ny, nz = ", nx, ny, nz

        ! Open the file for writing
        call nc_open(filename, ncid, writable=.TRUE.)

        ! Determine current writing time step
        n = nc_size(filename, "time", ncid)
        call nc_read(filename, "time", time_prev, start=[n], count=[1], ncid=ncid)
        if (abs(time - time_prev) .gt. 1e-5) n = n + 1

        ! Update the time step
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)

        call nc_write(filename, "ux", ylmo%dyn%now%ux, units="m/yr", &
            long_name="Velocity in x", dim1="xc", dim2="yc", dim3="zc", dim4="time", &
            start=[1, 1, 1, n], ncid=ncid)
        call nc_write(filename, "uy", ylmo%dyn%now%uy, units="m/yr", &
            long_name="Velocity in y", dim1="xc", dim2="yc", dim3="zc", dim4="time", &
            start=[1, 1, 1, n], ncid=ncid)
        call nc_write(filename, "uz", ylmo%dyn%now%uz, units="m/yr", &
            long_name="Velocity in z", dim1="xc", dim2="yc", dim3="zc", dim4="time", &
            start=[1, 1, 1, n], ncid=ncid)

        ! For some reason this gives crappy output (missings or 0) --> use 2D output for now
        ! call nc_write(filename, "visc", ylmo%mat%now%visc, units="Pa a", &
        !     long_name="Ice viscosity", dim1="xc", dim2="yc", dim3="zc", dim4="time", &
        !     start=[1, 1, 1, n], ncid=ncid)

    end subroutine write_step_3D

    subroutine write_step_2D_small(ylmo,filename,time)

        implicit none

        type(yelmo_class), intent(IN) :: ylmo
!         type(snapclim_class),   intent(IN) :: snp
!         type(marshelf_class),   intent(IN) :: mshlf
!         type(smbpal_class),     intent(IN) :: srf
        !type(sediments_class),  intent(IN) :: sed
        !type(geothermal_class), intent(IN) :: gthrm
        !type(isos_class),       intent(IN) :: isos

        character(len=*), intent(IN) :: filename
        real(wp),         intent(IN) :: time

        ! Local variables
        integer    :: ncid, n
        real(wp) :: time_prev

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

        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_grz",ylmo%tpo%now%mask_grz,units="",long_name="Grounding-zone mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m-1",long_name="Dragging coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff",ylmo%dyn%now%visc_eff,units="Pa a m",long_name="Effective viscosity (SSA)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically integrated velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)

        return

    end subroutine write_step_2D_small


    ! === CLIMATE ====

    subroutine calc_climate_standard(snp,smbp,mshlf,ylmo,time,time_bp)

        implicit none

        type(snapclim_class),       intent(INOUT) :: snp
        type(smbpal_class),         intent(INOUT) :: smbp
        type(marshelf_class),       intent(INOUT) :: mshlf
        type(yelmo_class),          intent(IN)    :: ylmo
        real(wp),                   intent(IN)    :: time
        real(wp),                   intent(IN)    :: time_bp

        ! Step 1: udpate the climate to the present time

        call snapclim_update(snp,z_srf=ylmo%tpo%now%z_srf,time=time_bp,domain=ylmo%par%domain,dx=yelmo1%grd%dx,basins=yelmo1%bnd%basins)

        ! Step 2: update the smb fields

        call smbpal_update_monthly(smbp,snp%now%tas,snp%now%pr, &
                                   ylmo%tpo%now%z_srf,ylmo%tpo%now%H_ice,time)

        ! Step 3: update the marine_shelf fields

        call marshelf_update_shelf(mshlf,ylmo%tpo%now%H_ice,ylmo%bnd%z_bed,ylmo%tpo%now%f_grnd, &
                        ylmo%bnd%basins,ylmo%bnd%z_sl,ylmo%grd%dx,snp%now%depth, &
                        snp%now%to_ann,snp%now%so_ann,dto_ann=snp%now%to_ann-snp1%clim0%to_ann)

        call marshelf_update(mshlf,ylmo%tpo%now%H_ice,ylmo%bnd%z_bed,ylmo%tpo%now%f_grnd, &
                             ylmo%bnd%regions,ylmo%bnd%basins,ylmo%bnd%z_sl,dx=yelmo1%grd%dx)

        return

    end subroutine calc_climate_standard

    subroutine calc_climate_ismip6(snp,smbp,mshlf,ismp,ylmo,time,time_bp,dTa,dTo)

        implicit none

        type(snapclim_class),       intent(INOUT) :: snp
        type(smbpal_class),         intent(INOUT) :: smbp
        type(marshelf_class),       intent(INOUT) :: mshlf
        type(ismip6_forcing_class), intent(INOUT) :: ismp
        type(yelmo_class),          intent(IN)    :: ylmo
        real(wp),                   intent(IN)    :: time
        real(wp),                   intent(IN)    :: time_bp
        real(wp), intent(IN), optional :: dTa
        real(wp), intent(IN), optional :: dTo

        ! Local variables
        real(wp), allocatable :: dts_now(:,:)
        real(wp), allocatable :: dpr_now(:,:)
        real(wp), allocatable :: dsmb_now(:,:)

        allocate(dts_now(ylmo%grd%nx,ylmo%grd%ny))
        allocate(dpr_now(ylmo%grd%nx,ylmo%grd%ny))
        allocate(dsmb_now(ylmo%grd%nx,ylmo%grd%ny))

        ! Step 1: set climate to present day from input fields
        ! and calculate present-day smb based on this climate.

        ! Set present-day climate
        !snp%now = snp%clim0

        ! Set present-day climate with optional constant atmospheric anomaly
        call snapclim_update(snp,z_srf=ylmo%tpo%now%z_srf,time=0.0_wp, &
                                        dx=yelmo1%grd%dx,basins=yelmo1%bnd%basins, &
                                        domain=ylmo%par%domain,dTa=dTa,dTo=0.0_wp)

        ! Calculate smb for present day
        call smbpal_update_monthly(smbp,snp%now%tas,snp%now%pr, &
                                   ylmo%tpo%now%z_srf,ylmo%tpo%now%H_ice,time)


        ! Step 2: update the ISMIP6 forcing to the current year

        call ismip6_forcing_update(ismp,time)

        ! Calculate anomaly fields accounting for elevation difference with reference topo
        dts_now  = ismp%ts%var(:,:,1,1)  !+ ismp%dts_dz%var(:,:,1,1)*(ylmo%tpo%now%z_srf-ismp%z_srf%var(:,:,1,1))
        dpr_now  = ismp%pr%var(:,:,1,1)
        dsmb_now = ismp%smb%var(:,:,1,1) !+ ismp%dsmb_dz%var(:,:,1,1)*(ylmo%tpo%now%z_srf-ismp%z_srf%var(:,:,1,1))

        ! Step 3: apply ISMIP6 anomalies to climate and smb fields
        ! (apply to climate just for consistency)

        ! Update climatic fields
        do m = 1,12
            snp%now%tas(:,:,m) = snp%now%tas(:,:,m) + dts_now
            snp%now%pr(:,:,m)  = snp%now%pr(:,:,m)  + dpr_now/365.0 ! [mm/yr] => [mm/d]
        end do

        snp%now%ta_ann = sum(snp%now%tas,dim=3) / 12.0_wp
        if (trim(domain) .eq. "Antarctica") then
            snp%now%ta_sum  = sum(snp%now%tas(:,:,[12,1,2]),dim=3)/3.0  ! Antarctica summer
        else
            snp%now%ta_sum  = sum(snp%now%tas(:,:,[6,7,8]),dim=3)/3.0  ! NH summer
        end if
        snp%now%pr_ann = sum(snp%now%pr,dim=3)  / 12.0 * 365.0     ! [mm/d] => [mm/yr]

        ! Update smb fields
        smbp%ann%smb  = smbp%ann%smb  + dsmb_now*1.0/(yelmo1%bnd%c%conv_we_ie*1e-3) ! [m ie/yr] => [mm we/yr]
        smbp%ann%tsrf = smbp%ann%tsrf + dts_now

        ! Step 4: update marine_shelf based on ISMIP6 fields
        ! (no need to mix with present-day climate, since ismip6 includes the
        !  reference ocean fields with internal depth dimension)

        ! Update marine_shelf shelf fields
        ! jablasco: set anomaly to zero
        ! robinson: dto_ann=ismp%to%var(:,:,:,1)-ismp%to_ref%var(:,:,:,1)
        ! jablasco: volvamos al ppio! dto_ann=ismp%to%var(:,:,:,1)*0.0
        call marshelf_update_shelf(mshlf,ylmo%tpo%now%H_ice,ylmo%bnd%z_bed,ylmo%tpo%now%f_grnd, &
                        ylmo%bnd%basins,ylmo%bnd%z_sl,ylmo%grd%dx,-ismp%to%z, &
                        ismp%to%var(:,:,:,1),ismp%so%var(:,:,:,1), &
                        dto_ann=ismp%to%var(:,:,:,1)-ismp%to_ref%var(:,:,:,1), &
                        tf_ann=ismp%tf%var(:,:,:,1))

        ! Update temperature forcing field with tf_corr and tf_corr_basin
        mshlf%now%tf_shlf = mshlf%now%tf_shlf + mshlf%now%tf_corr + mshlf%now%tf_corr_basin

        if (present(dTo)) then
            ! Update temperature fields with hysteresis anomaly
            mshlf1%now%T_shlf  = mshlf1%now%T_shlf  + dTo
            mshlf1%now%dT_shlf = mshlf1%now%dT_shlf + dTo
            mshlf1%now%tf_shlf = mshlf1%now%tf_shlf + dTo
        end if

        ! Update bmb_shlf and mask_ocn
        call marshelf_update(mshlf,ylmo%tpo%now%H_ice,ylmo%bnd%z_bed,ylmo%tpo%now%f_grnd, &
                             ylmo%bnd%regions,ylmo%bnd%basins,ylmo%bnd%z_sl,dx=ylmo%grd%dx)

        return

    end subroutine calc_climate_ismip6

    subroutine load_tf_corr_from_restart(tf_corr,file_restart,domain,grid_name)

        use coordinates_mapping_scrip, only : map_scrip_class, map_scrip_init, map_scrip_field, &
                                            gen_map_filename, nc_read_interp

        implicit none

        real(wp),         intent(INOUT) :: tf_corr(:,:)
        character(len=*), intent(IN)    :: file_restart
        character(len=*), intent(IN)    :: domain
        character(len=*), intent(IN)    :: grid_name

        ! Local variables
        integer :: i0, i1, n, nx, ny
        character(len=1024) :: path_tf_corr

        character(len=56) :: restart_domain
        character(len=56) :: restart_grid_name
        type(map_scrip_class) :: mps
        logical :: restart_interpolated

        nx = size(tf_corr,1)
        ny = size(tf_corr,2)

        ! Load restart file grid attributes
        if (nc_exists_attr(file_restart,"domain")) then
            call nc_read_attr(file_restart,"domain",    restart_domain)
        else
            restart_domain = trim(domain)
        end if

        if (nc_exists_attr(file_restart,"grid_name")) then
            call nc_read_attr(file_restart,"grid_name", restart_grid_name)
        else
            restart_grid_name = trim(grid_name)
        end if

        if (trim(restart_grid_name) .ne. trim(grid_name)) then
            restart_interpolated = .TRUE.
        else
            restart_interpolated = .FALSE.
        end if

        if (restart_interpolated) then
            ! Load the scrip map from file (should already have been generated via cdo externally)
            call map_scrip_init(mps,restart_grid_name,grid_name, &
                                    method="con",fldr="maps",load=.TRUE.)
        end if

        ! Get folder holding the restart file and
        ! append filename holding the tf_corr field we want
        ! (should be in yelmo2D.nc)
        i0 = index(file_restart,"/",back=.TRUE.)
        path_tf_corr = file_restart(1:i0)//"yelmo2D.nc"

        ! Old alternative method that only works if restart file
        ! is named "yelmo_restart.nc"
        !path_tf_corr = file_restart
        !call nml_replace(path_tf_corr,"yelmo_restart.nc","yelmo2D.nc")


        ! Load the tf_corr field from the last timestep of the yelmo2D file
        n = nc_size(path_tf_corr,"time")

        if (restart_interpolated) then
            call nc_read_interp(path_tf_corr,"tf_corr",tf_corr,start=[1,1,n],count=[nx,ny,1],mps=mps)
        else
            call nc_read(path_tf_corr,"tf_corr",tf_corr,start=[1,1,n],count=[nx,ny,1])
        end if

        ! Write summary to screen
        write(*,*) "load_tf_corr_from_restart:: loaded tf_corr field."
        write(*,*) "path_tf_corr: ", trim(path_tf_corr)
        write(*,*) "tf_corr: ", minval(tf_corr), maxval(tf_corr)

        return

    end subroutine load_tf_corr_from_restart


! ======================================
!
! Hysteresis related output routines
!
! ======================================

! Get the indices of the cropped WAIS region
subroutine get_cropwais_indices(mask, val, i1, i2, j1, j2, xmargin, ymargin)
    implicit none
    real(wp), intent(IN)  :: mask(:,:), val
    integer,  intent(IN)  :: xmargin, ymargin
    integer, intent(OUT)  :: i1, i2, j1, j2

    ! Local variables
    integer :: i, j

    ! Initialize indices with values that will be overwritten
    i1 = 1000000
    i2 = 0
    j1 = 1000000
    j2 = 0

    do i = 1, size(mask,1)
    do j = 1, size(mask,2)
        if (mask(i,j) .eq. val) then
            i1 = min(i1,i)
            i2 = max(i2,i)
            j1 = min(j1,j)
            j2 = max(j2,j)
        end if
    end do
    end do

    ! Add margin
    i1 = max(i1-xmargin,1)
    i2 = min(i2+xmargin,size(mask,1))
    j1 = max(j1-ymargin,1)
    j2 = min(j2+ymargin,size(mask,2))

end subroutine get_cropwais_indices

subroutine yx_hyst_write_step_2Dwais(ylmo, mshlf, filename, time, i1, i2, j1, j2)

    implicit none

    type(yelmo_class),      intent(IN) :: ylmo
    type(marshelf_class),   intent(IN) :: mshlf
    character(len=*),       intent(IN) :: filename
    real(prec),             intent(IN) :: time
    integer,                intent(IN) :: i1, i2, j1, j2

    ! Local variables
    integer    :: ncid, n
    real(prec) :: time_prev
    character(len=56), allocatable :: names(:)

    ! Open the file for writing
    call nc_open(filename, ncid, writable=.TRUE.)

    ! Determine current writing time step
    n = nc_size(filename, "time", ncid)
    call nc_read(filename, "time", time_prev, start=[n], count=[1], ncid=ncid)
    if (abs(time-time_prev) .gt. 1e-5) n = n+1

    ! Update the time step
    call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

    call nc_write(filename, "H_ice", ylmo%tpo%now%H_ice(i1:i2, j1:j2), &
        units="m", long_name="Ice thickness", dim1="xc", dim2="yc", &
        dim3="time", start=[1,1,n], ncid=ncid)

    call nc_write(filename, "z_bed", ylmo%bnd%z_bed(i1:i2, j1:j2), &
        units="m", long_name="Bedrock elevation", dim1="xc", dim2="yc", &
        dim3="time", start=[1,1,n], ncid=ncid)

    call nc_write(filename, "f_grnd", ylmo%tpo%now%f_grnd(i1:i2, j1:j2), &
        units="1", long_name="Grounding line fraction", dim1="xc", dim2="yc", &
        dim3="time", start=[1,1,n], ncid=ncid)

    call nc_write(filename, "uxy_s", ylmo%dyn%now%uxy_s(i1:i2, j1:j2), &
        units="m/yr", long_name="Surface velocity magnitude", dim1="xc", dim2="yc", &
        dim3="time", start=[1,1,n], ncid=ncid)

    call nc_write(filename, "T_shlf", mshlf%now%T_shlf(i1:i2, j1:j2), &
        units="K", long_name="Ocean temperature at shelf interface", &
        dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)

    call nc_write(filename, "smb", ylmo%bnd%smb(i1:i2, j1:j2), &
        units="m/yr ice equiv.", long_name="Surface mass balance", dim1="xc", dim2="yc",&
        dim3="time", start=[1,1,n], ncid=ncid)

    call nc_write(filename, "bmb", ylmo%tpo%now%bmb(i1:i2, j1:j2), &
        units="m/yr ice equiv.", long_name="Basal mass balance", dim1="xc", dim2="yc", &
        dim3="time", start=[1,1,n], ncid=ncid)

    call nc_write(filename, "z_sl", ylmo%bnd%z_sl(i1:i2, j1:j2), &
        units="m", long_name="Sea level rel. to present", dim1="xc", dim2="yc", &
        dim3="time", start=[1,1,n], ncid=ncid)

    call nc_close(ncid)
    return
end subroutine yx_hyst_write_step_2Dwais

subroutine yx_hyst_write_step_2D_combined(ylmo,isos,snp,mshlf,srf,filename,time)

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
        real(wp), intent(IN) :: time

        ! Local variables
        integer  :: ncid, n
        real(wp) :: time_prev

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

        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"mb_net",ylmo%tpo%now%mb_net,units="m",long_name="Applied net mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"H_grnd",ylmo%tpo%now%H_grnd,units="m",long_name="Ice thickness overburden", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="bar",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"cmb",ylmo%tpo%now%cmb,units="m/a ice equiv.",long_name="Calving mass balance rate", &
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
        call nc_write(filename,"visc_eff_int",ylmo%dyn%now%visc_eff_int,units="Pa a m",long_name="Depth-integrated effective viscosity (SSA)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud",ylmo%dyn%now%taud,units="Pa",long_name="Driving stress", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taub",ylmo%dyn%now%taub,units="Pa",long_name="Basal stress", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically integrated velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        !call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
        !              dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"T_prime_b",ylmo%thrm%now%T_prime_b,units="K",long_name="Basal homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        !call nc_write(filename,"Q_strn",ylmo%thrm%now%Q_strn/(rho_ice*ylmo%thrm%now%cp),units="K a-1",long_name="Strain heating", &
        !              dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
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

        call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW m-2",long_name="Geothermal heat flux", &
                      dim1="xc",dim2="yc",start=[1,1],ncid=ncid)

        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a ice equiv.",long_name="Basal mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! External data
        call nc_write(filename,"dzbdt",isos%output%dwdt,units="m/a",long_name="Bedrock uplift rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dT_shlf",mshlf%now%dT_shlf,units="K",long_name="Shelf temperature anomaly", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_ocn",mshlf%now%mask_ocn,units="",long_name="Marine-shelf ocean mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

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

        call nc_write(filename,"dzsdx",ylmo%tpo%now%dzsdx,units="m/m",long_name="Surface slope", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dzsdy",ylmo%tpo%now%dzsdy,units="m/m",long_name="Surface slope", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

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

    end subroutine yx_hyst_write_step_2D_combined

    subroutine yx_hyst_write_step_2D_combined_small(ylmo,isos,snp,mshlf,srf,filename,time)

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
        real(wp), intent(IN) :: time

        ! Local variables
        integer  :: ncid, n
        real(wp) :: time_prev

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

        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m-1",long_name="Basal friction coefficient", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"T_prime_b",ylmo%thrm%now%T_prime_b,units="K",long_name="Basal homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Boundaries
        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb",ylmo%bnd%smb,units="m/a ice equiv.",long_name="Surface mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a ice equiv.",long_name="Basal mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! External data
!         call nc_write(filename,"dzbdt",isos%now%dzbdt,units="m/a",long_name="Bedrock uplift rate", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"Ta_ann",snp%now%ta_ann,units="K",long_name="Near-surface air temperature (ann)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"Ta_sum",snp%now%ta_sum,units="K",long_name="Near-surface air temperature (sum)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"Pr_ann",snp%now%pr_ann*1e-3,units="m/a water equiv.",long_name="Precipitation (ann)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dTa_ann",snp%now%ta_ann-snp%clim0%ta_ann,units="K",long_name="Near-surface air temperature anomaly (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dTa_sum",snp%now%ta_sum-snp%clim0%ta_sum,units="K",long_name="Near-surface air temperature anomaly (sum)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dPr_ann",(snp%now%pr_ann-snp%clim0%pr_ann)*1e-3,units="m/a water equiv.",long_name="Precipitation anomaly (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dT_shlf",mshlf%now%dT_shlf,units="K",long_name="Shelf temperature anomaly", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)

        return

    end subroutine yx_hyst_write_step_2D_combined_small

    subroutine yx_hyst_write_yelmo_init_1D_combined(dom,filename,time_init,units,mask)

        implicit none

        type(yelmo_class), intent(IN) :: dom
        character(len=*),  intent(IN) :: filename, units
        real(wp),          intent(IN) :: time_init
        logical,           intent(IN) :: mask(:,:)

        ! Initialize netcdf file and dimensions
        call nc_create(filename)
        call nc_write_dim(filename,"xc",        x=dom%grd%xc*1e-3,    units="kilometers")
        call nc_write_dim(filename,"yc",        x=dom%grd%yc*1e-3,    units="kilometers")
        call nc_write_dim(filename,"zeta",      x=dom%par%zeta_aa,    units="1")
        call nc_write_dim(filename,"time",      x=time_init,dx=1.0_wp,nx=1,units=trim(units),unlimited=.TRUE.)

        ! Static information
        call nc_write(filename,"mask", mask,  units="1",long_name="Region mask",dim1="xc",dim2="yc")

        return

    end subroutine yx_hyst_write_yelmo_init_1D_combined

    subroutine yx_hyst_write_step_1D_combined(ylm,hyst,snp,filename,time)

        implicit none

        type(yelmo_class),    intent(IN) :: ylm
        type(hyster_class),   intent(IN) :: hyst
        type(snapclim_class), intent(IN) :: snp
        character(len=*),     intent(IN) :: filename
        real(wp),             intent(IN) :: time

        ! Local variables
        integer  :: ncid, n, k
        real(wp) :: time_prev
        type(yregions_class) :: reg

        ! Assume region to write is the global region of yelmo
        reg = ylm%reg

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid)
        if (abs(time-time_prev).gt.1e-5) n = n+1

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)


        ! ===== Hyst / forcing variables =====

        call nc_write(filename,"hyst_f_now",hyst%f_now,units="K",long_name="hyst: forcing value", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"hyst_df_dt",hyst%df_dt*1e6,units="K/(1e6 a)",long_name="hyst: forcing rate of change", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"hyst_dv_dt",hyst%dv_dt,units="Gt/a",long_name="hyst: volume rate of change", &
                      dim1="time",start=[n],ncid=ncid)

        ! ===== Total ice variables =====

        call nc_write(filename,"H_ice",reg%H_ice,units="m",long_name="Mean ice thickness", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"z_srf",reg%z_srf,units="m",long_name="Mean surface elevation", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dHidt",reg%dHidt,units="m/a",long_name="Mean rate ice thickness change", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"H_ice_max",reg%H_ice_max,units="m/a",long_name="Max ice thickness", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dzsdt",reg%dzsdt,units="m/a",long_name="Mean rate surface elevation change", &
                      dim1="time",start=[n],ncid=ncid)

        call nc_write(filename,"V_ice",reg%V_ice*1e-6,units="1e6 km^3",long_name="Ice volume", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"A_ice",reg%A_ice*1e-6,units="1e6 km^2",long_name="Ice area", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dVidt",reg%dVidt,units="km^3/a",long_name="Rate volume change", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"fwf",reg%fwf,units="Sv",long_name="Rate volume change", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"V_sl",reg%V_sl*1e-6,units="1e6 km^3",long_name="Ice volume above flotation", &
                      dim1="time",start=[n],ncid=ncid)

        call nc_write(filename,"uxy_bar",reg%uxy_bar,units="m/a",long_name="Mean depth-ave velocity", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_s",reg%uxy_s,units="m/a",long_name="Mean surface velocity", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_b",reg%uxy_b,units="m/a",long_name="Mean basal velocity", &
                      dim1="time",start=[n],ncid=ncid)

        call nc_write(filename,"z_bed",reg%z_bed,units="m",long_name="Mean bedrock elevation", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"smb",reg%smb,units="m/a",long_name="Mean surface mass balance", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"T_srf",reg%T_srf,units="K",long_name="Mean surface temperature", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"bmb",reg%bmb,units="m/a",long_name="Mean total basal mass balance", &
                      dim1="time",start=[n],ncid=ncid)



        ! ===== Grounded ice variables =====

        call nc_write(filename,"H_ice_g",reg%H_ice_g,units="m",long_name="Mean ice thickness (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"z_srf_g",reg%z_srf_g,units="m",long_name="Mean surface elevation (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"V_ice_g",reg%V_ice_g*1e-6,units="1e6 km^3",long_name="Ice volume (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"A_ice_g",reg%A_ice_g*1e-6,units="1e6 km^2",long_name="Ice area (grounded)", &
                      dim1="time",start=[n],ncid=ncid)

        call nc_write(filename,"uxy_bar_g",reg%uxy_bar_g,units="m/a",long_name="Mean depth-ave velocity (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_s_g",reg%uxy_s_g,units="m/a",long_name="Mean surface velocity (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_b_g",reg%uxy_b_g,units="m/a",long_name="Mean basal velocity (grounded)", &
                      dim1="time",start=[n],ncid=ncid)

        call nc_write(filename,"f_pmp",reg%f_pmp,units="1",long_name="Temperate fraction (grounded)", &
                      dim1="time",start=[n],ncid=ncid)

        call nc_write(filename,"H_w",reg%H_w,units="m",long_name="Mean basal water thickness (grounded)", &
                      dim1="time",start=[n],ncid=ncid)

        call nc_write(filename,"bmb_g",reg%bmb_g,units="m/a",long_name="Mean basal mass balance (grounded)", &
                      dim1="time",start=[n],ncid=ncid)

        ! ===== Floating ice variables =====

        call nc_write(filename,"H_ice_f",reg%H_ice_f,units="m",long_name="Mean ice thickness (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"V_ice_f",reg%V_ice_f*1e-6,units="1e6 km^3",long_name="Ice volume (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"A_ice_f",reg%A_ice_f*1e-6,units="1e6 km^2",long_name="Ice area (floating)", &
                      dim1="time",start=[n],ncid=ncid)

        call nc_write(filename,"uxy_bar_f",reg%uxy_bar_f,units="m/a",long_name="Mean depth-ave velocity (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_s_f",reg%uxy_s_f,units="m/a",long_name="Mean surface velocity (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_b_f",reg%uxy_b_f,units="m/a",long_name="Mean basal velocity (floating)", &
                      dim1="time",start=[n],ncid=ncid)

        call nc_write(filename,"z_sl",reg%z_sl,units="m",long_name="Mean sea level (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"bmb_shlf",reg%bmb_shlf,units="m/a",long_name="Mean basal mass balance (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"T_shlf",reg%T_shlf,units="K",long_name="Mean marine shelf temperature (floating)", &
                      dim1="time",start=[n],ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)

        return

    end subroutine yx_hyst_write_step_1D_combined


    ! ===== ISMIP6 output routines =========

    subroutine write_step_2D_ismip6(ylmo,filename,time)

        implicit none

        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(wp),          intent(IN) :: time

        ! Local variables
        integer  :: ncid, n
        real(wp) :: time_prev
        real(wp) :: rho_ice
        real(wp) :: density_corr
        real(wp) :: m3yr_to_kgs
        real(wp) :: ismip6_correction
        real(wp) :: yr_to_sec

        ! ISMIP6 variables
        real(wp), allocatable :: bmb_grnd_masked(:,:)
        real(wp), allocatable :: bmb_shlf_masked(:,:)
        real(wp), allocatable :: z_base(:,:)
        real(wp), allocatable :: T_top_ice(:,:)
        real(wp), allocatable :: T_base_grnd(:,:)
        real(wp), allocatable :: T_base_flt(:,:)
        real(wp), allocatable :: flux_grl(:,:)

        ! Allocate and initialize local arrays
        allocate(bmb_grnd_masked(ylmo%grd%nx,ylmo%grd%ny))
        allocate(bmb_shlf_masked(ylmo%grd%nx,ylmo%grd%ny))
        allocate(z_base(ylmo%grd%nx,ylmo%grd%ny))
        allocate(T_top_ice(ylmo%grd%nx,ylmo%grd%ny))
        allocate(T_base_grnd(ylmo%grd%nx,ylmo%grd%ny))
        allocate(T_base_flt(ylmo%grd%nx,ylmo%grd%ny))
        allocate(flux_grl(ylmo%grd%nx,ylmo%grd%ny))

        bmb_shlf_masked   = 0.0_wp
        bmb_grnd_masked   = 0.0_wp
        z_base            = 0.0_wp
        T_top_ice         = 0.0_wp
        T_base_grnd       = 0.0_wp
        T_base_flt        = 0.0_wp
        flux_grl          = 0.0_wp

        ! === Data conversion factors ========================================

        rho_ice             = 917.0             ! ice density kg/m3
        m3yr_to_kgs         = 3.2e-5            ! m3/yr of pure water to kg/s
        density_corr        = rho_ice/1000.0    ! ice density correction with pure water
        ismip6_correction   = m3yr_to_kgs*density_corr
        yr_to_sec           = 31556952.0

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid)
        if (abs(time-time_prev).gt.1e-5) n = n+1

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)


        ! Write ismip6 variables

        ! == yelmo_topography ==
        call nc_write(filename,"lithk",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                        standard_name="land_ice_thickness", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"orog",ylmo%tpo%now%z_srf,units="m",long_name="Surface evelation", &
                        standard_name="surface_altitude", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! ajr, should we just be writing the Yelmo variable ylmo%tpo%now%z_base??
        where(ylmo%tpo%now%f_grnd_bmb .gt. 0.0) z_base = ylmo%bnd%z_bed
        call nc_write(filename,"base",z_base,units="m",long_name="Base elevation", &
                        standard_name="base_altitude", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"topg",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                        standard_name="bedrock_altitude", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_GHF ==
        call nc_write(filename,"hfgeoubed",ylmo%bnd%Q_geo*1e3,units="W/m^2",long_name="Geothermal heat flux", &
                        standard_name="upward_geothermal_heat_flux_at_ground_level", dim1="xc",dim2="yc",start=[1,1],ncid=ncid)

        ! == yelmo_mass_balances ==

        where(ylmo%tpo%now%f_grnd_bmb .gt. 0.0) bmb_grnd_masked = ylmo%thrm%now%bmb_grnd
        where(ylmo%tpo%now%H_ice .gt. 0.0 .and. ylmo%tpo%now%f_grnd_bmb .eq. 0.0) bmb_shlf_masked = ylmo%bnd%bmb_shlf
        call nc_write(filename,"acabf",ylmo%bnd%smb*ismip6_correction,units="kg m-2 s-1",long_name="Surface mass balance flux", &
                        standard_name="land_ice_surface_specific_mass_balance_flux", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"libmassbfgr",bmb_grnd_masked*ismip6_correction,units="kg m-2 s-1",long_name="Basal mass balance flux beneath grounded ice", &
                        standard_name="land_ice_basal_specific_mass_balance_flux", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"libmassbffl",bmb_shlf_masked*ismip6_correction,units="kg m-2 s-1",long_name="Basal mass balance flux beneath floating ice", &
                        standard_name="float_ice_basal_specific_mass_balance_flux", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_velocities ==
        call nc_write(filename,"dlithkdt",ylmo%tpo%now%dHidt/yr_to_sec,units="m s-1",long_name="Ice thickness imbalance", &
                        standard_name="tendency_of_land_ice_thickness", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! comparison dHdt m/yr
        !call nc_write(filename,"dHdt",ylmo%tpo%now%dHidt,units="m yr-1",long_name="Ice thickness imbalance", &
        !              standard_name="tendency_of_land_ice_thickness", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"xvelsurf",ylmo%dyn%now%ux_s/yr_to_sec,units="m s-1",long_name="Surface velocity in x", &
                        standard_name="land_ice_surface_x_velocity", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"yvelsurf",ylmo%dyn%now%uy_s/yr_to_sec,units="m s-1",long_name="Surface velocity in y", &
                        standard_name="land_ice_surface_y_velocity", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"zvelsurf",ylmo%dyn%now%uz_s/yr_to_sec,units="m s-1",long_name="Surface velocity in z", &
                        standard_name="land_ice_surface_upward_velocity", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"xvelbase",ylmo%dyn%now%ux_b/yr_to_sec,units="m s-1",long_name="Basal velocity in x", &
                        standard_name="land_ice_surface_x_velocity", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"yvelbase",ylmo%dyn%now%uy_b/yr_to_sec,units="m s-1",long_name="Basal velocity in y", &
                        standard_name="land_ice_surface_y_velocity", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"zvelbase",ylmo%dyn%now%uz_b/yr_to_sec,units="m s-1",long_name="Basal velocity in z", &
                        standard_name="land_ice_basal_upward_velocity", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"xvelmean",ylmo%dyn%now%ux_bar/yr_to_sec,units="m s-1",long_name="Mean velocity in x", &
                        standard_name="land_ice_vertical_mean_x_velocity", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"yvelmean",ylmo%dyn%now%uy_bar/yr_to_sec,units="m s-1",long_name="Mean velocity in y", &
                        standard_name="land_ice_vertical_mean_y_velocity", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! == yelmo_ice_temperature ==

        ! Local variables
        where(ylmo%tpo%now%f_grnd_bmb .gt. 0.0) T_base_grnd = ylmo%thrm%now%T_ice(:,:,1)
        where(ylmo%tpo%now%H_ice .gt. 0.0 .and. ylmo%tpo%now%f_grnd_bmb .eq. 0.0) T_base_flt = ylmo%thrm%now%T_ice(:,:,1)
        where(ylmo%tpo%now%H_ice .gt. 0.0) T_top_ice = ylmo%thrm%now%T_ice(:,:,ylmo%dyn%par%nz_aa)
        call nc_write(filename,"litemptop",T_top_ice,units="K",long_name="Surface temperature", &
                        standard_name="temperature_at_top_of_ice_sheet_model", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"litempbotgr",T_base_grnd,units="K",long_name="Basal temperature beneath grounded ice sheet", &
                        standard_name="land_ice_basal_temperature", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"litempbotfl",T_base_flt,units="K",long_name="Basal temperature beneath floating ice shelf", &
                        standard_name="floating_ice_basal_temperature", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_drag_and_fluxes ==
        ! ajr, flux at gl not calculated correctly?
        where(ylmo%tpo%now%f_grnd_bmb .gt. 0.0 .and. ylmo%tpo%now%f_grnd_bmb .lt. 1.0) flux_grl = ylmo%tpo%now%dHidt
        call nc_write(filename,"strbasemag",ylmo%dyn%now%taub,units="Pa",long_name="Basal drag", &
                        standard_name="Basal drag", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"licalvf",ylmo%tpo%now%cmb*ismip6_correction,units="kg m-2 s-1",long_name="Calving flux", &
                        standard_name="land_ice_specific_mass_flux_due_to_calving", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"lifmassbf",(ylmo%tpo%now%fmb+ylmo%tpo%now%cmb)*ismip6_correction,units="kg m-2 s-1",long_name="Ice front melt and calving flux", &
                        standard_name="land_ice_specific_mass_flux_due_to_calving_and_ice_front_melting", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ligroundf",flux_grl*ismip6_correction,units="kg m-2 s-1",long_name="Grounding line flux", &
                        standard_name="land_ice_specific_mass_flux_at_grounding_line", dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_masks ==
        call nc_write(filename,"sftgif",ylmo%tpo%now%f_ice,units="1",long_name="Land ice area fraction", &
                        standard_name="land_ice_area_fraction",dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"sftgrf",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded ice sheet area fraction", &
                        standard_name="grounded_ice_sheet_area_fraction",dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"sftflf",MAX(ylmo%tpo%now%f_ice-ylmo%tpo%now%f_grnd,0.0),units="1",long_name="Floating ice sheet area fraction", &
                        standard_name="floating_ice_sheet_area_fraction",dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_coll",ismp1%mask_shlf%var(:,:,1,1),units="1",long_name="Collapse mask", &
                        standard_name="",dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! ATM and OCN test files
        ! call nc_write(filename,"T_atm_snap",ismp1%ts%var(:,:,1,1),units="K",long_name="Surface temperature anomaly", &
        !                 dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"Tf_snap",ismp1%tf%var(:,:,1,1),units="K",long_name="TF map", &
        !                 dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"T_srf",yelmo1%bnd%T_srf,units="K",long_name="Surface temperature Yelmo", &
        !                 dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"Tf_shlf",mshlf%now%tf_shlf,units="K",long_name="Shelf temperature", &
        !                 dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"tf_corr",mshlf%now%tf_corr,units="K",long_name="Shelf thermal forcing correction factor", &
        !                 dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)

        return

    end subroutine write_step_2D_ismip6

    subroutine write_1D_ismip6(dom,filename,time)

        implicit none

        type(yelmo_class), intent(IN) :: dom
        character(len=*),  intent(IN) :: filename
        real(wp),          intent(IN) :: time

        ! Local variables
        type(yregions_class) :: reg

        integer  :: ncid, n
        real(wp) :: time_prev
        real(wp) :: rho_ice
        real(wp) :: density_corr
        real(wp) :: m3yr_to_kgs
        real(wp) :: ismip6_correction
        real(wp) :: yr_to_sec

        real(wp) :: m3_km3
        real(wp) :: m2_km2

        integer  :: npts_tot
        integer  :: npts_grnd
        integer  :: npts_flt
        integer  :: npts_grl
        integer  :: npts_frnt

        real(wp) :: dx
        real(wp) :: dy
        real(wp) :: smb_tot
        real(wp) :: bmb_tot
        real(wp) :: bmb_shlf_t

        real(wp) :: A_ice_grl
        real(wp) :: flux_grl
        real(wp) :: A_ice_frnt
        real(wp) :: calv_flt
        real(wp) :: flux_frnt

        logical, allocatable :: mask_tot(:,:)
        logical, allocatable :: mask_grnd(:,:)
        logical, allocatable :: mask_flt(:,:)
        logical, allocatable :: mask_grl(:,:)
        logical, allocatable :: mask_frnt(:,:)

        dx = dom%grd%dx
        dy = dom%grd%dy

        ! Allocate variables

        allocate(mask_tot(dom%grd%nx,dom%grd%ny))
        allocate(mask_grnd(dom%grd%nx,dom%grd%ny))
        allocate(mask_flt(dom%grd%nx,dom%grd%ny))
        allocate(mask_grl(dom%grd%nx,dom%grd%ny))
        allocate(mask_frnt(dom%grd%nx,dom%grd%ny))

        ! === Data conversion factors ========================================

        rho_ice             = 917.0             ! ice density kg/m3
        m3yr_to_kgs         = 3.2e-5            ! m3/yr of pure water to kg/s
        density_corr        = rho_ice/1000.0    ! ice density correction with pure water
        ismip6_correction   = m3yr_to_kgs*density_corr
        yr_to_sec           = 31556952.0

        m3_km3              = 1e-9
        m2_km2              = 1e-6

        ! 1. Determine regional values of variables

        ! Take the global regional data object that
        ! is calculated over the whole domain at each timestep
        reg = dom%reg

        ! Assign masks of interest

        mask_tot  = (dom%tpo%now%H_ice .gt. 0.0)
        mask_grnd = (dom%tpo%now%H_ice .gt. 0.0 .and. dom%tpo%now%f_grnd .gt. 0.0)
        mask_flt  = (dom%tpo%now%H_ice .gt. 0.0 .and. dom%tpo%now%f_grnd .eq. 0.0)
        mask_grl  = (dom%tpo%now%f_grnd_bmb .gt. 0.0 .and. dom%tpo%now%f_grnd_bmb .lt. 1.0)
        mask_frnt = (dom%tpo%now%H_ice .gt. 0.0 .and. dom%tpo%now%mask_bed .eq. 4)

        ! Determine number of points at grl and frnt
        npts_tot  = count(mask_tot)
        npts_grnd = count(mask_grnd)
        npts_flt  = count(mask_flt)
        npts_grl  = count(mask_grl)
        npts_frnt = count(mask_frnt)

        ! Calculate additional variables of interest for ISMIP6

        ! To do: these variables should be defined locally!
        ! ISMIP6 boundary: jablasco
        smb_tot     = sum(dom%bnd%smb,    mask=mask_tot)*(dx*dy)    ! m^3/yr: flux
        bmb_tot     = sum(dom%tpo%now%bmb,mask=mask_tot)*(dx*dy)    ! m^3/yr: flux
        bmb_shlf_t  = sum(dom%tpo%now%bmb,mask=mask_flt)*(dx*dy)    ! m^3/yr: flux

        if (npts_grl .gt. 0) then

            A_ice_grl = count(dom%tpo%now%H_ice .gt. 0.0 .and. mask_grl)*dx*dy*m2_km2 ! [km^2]
            !flux_grl  = sum(dom%tpo%now%H_ice,mask=mask_grl)*(dx*dy)                 ! m^3/yr: flux
            flux_grl  = sum(dom%dyn%now%uxy_bar*dom%tpo%now%H_ice,mask=mask_grl)*dx   ! m^3/yr: flux

            ! ajr, why only *dx above?

        else

            A_ice_grl = 0.0_wp
            flux_grl  = 0.0_wp

        end if

        ! ===== Frontal ice-shelves variables =====

        if (npts_frnt .gt. 0) then

            A_ice_frnt = count(dom%tpo%now%H_ice .gt. 0.0 .and. mask_frnt)*dx*dy*m2_km2         ! [km^2]
            calv_flt   = sum(dom%tpo%now%cmb_flt*dom%tpo%now%H_ice,mask=mask_frnt)*dx          ! m^3/yr: flux [m-1 yr-1]
            flux_frnt  = calv_flt+sum(dom%tpo%now%fmb*dom%tpo%now%H_ice,mask=mask_frnt)*dx      ! m^3/yr: flux [m-1 yr-1]

            ! ajr, why only *dx above?

        else

            A_ice_frnt = 0.0_wp
            calv_flt   = 0.0_wp
            flux_frnt  = 0.0_wp

        end if


        ! 2. Begin writing step

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid)
        if (abs(time-time_prev).gt.1e-5) n = n+1

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        call nc_write(filename,"V_sl",reg%V_sl*1e-6,units="1e6 km^3",long_name="Ice volume above flotation", &
                          dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"V_sle",reg%V_sle,units="m sle",long_name="Sea-level equivalent volume", &
                        dim1="time",start=[n],ncid=ncid)

        ! ===== Mass =====
        call nc_write(filename,"lim",reg%V_ice*rho_ice*1e9,units="kg",long_name="Total ice mass", &
                      standard_name="land_ice_mass",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"limnsw",reg%V_sl*rho_ice*1e9,units="kg",long_name="Mass above flotation", &
                      standard_name="land_ice_mass_not_displacing_sea_water",dim1="time",start=[n],ncid=ncid)

        ! ===== Area ====
        call nc_write(filename,"iareagr",reg%A_ice_g*1e6,units="m^2",long_name="Grounded ice area", &
                      standard_name="grounded_ice_sheet_area",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"iareafl",reg%A_ice_f*1e6,units="m^2",long_name="Floating ice area", &
                      standard_name="floating_ice_shelf_area",dim1="time",start=[n],ncid=ncid)

        ! ==== Fluxes ====
        call nc_write(filename,"tendacabf",smb_tot*ismip6_correction,units="kg s-1",long_name="Total SMB flux", &
                      standard_name="tendency_of_land_ice_mass_due_to_surface_mass_balance",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"tendlibmassbf ",bmb_tot*ismip6_correction,units="kg s-1",long_name="Total BMB flux", &
                      standard_name="tendency_of_land_ice_mass_due_to_basal_mass_balance",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"tendlibmassbffl",bmb_shlf_t*ismip6_correction,units="kg s-1",long_name="Total BMB flux beneath floating ice", &
                      standard_name="tendency_of_land_ice_mass_due_to_basal_mass_balance",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"tendlicalvf",calv_flt*ismip6_correction,units="kg s-1",long_name="Total calving flux", &
                      standard_name="tendency_of_land_ice_mass_due_to_calving",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"tendlifmassbf",flux_frnt*ismip6_correction,units="kg s-1",long_name="Total calving and ice front melting flux", &
                      standard_name="tendency_of_land_ice_mass_due_to_calving_and_ice_front_melting",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"tendligroundf",flux_grl*ismip6_correction,units="kg s-1",long_name="Total grounding line flux", &
                      standard_name="tendency_of_grounded_ice_mass",dim1="time",start=[n],ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)

        return

    end subroutine write_1D_ismip6

end program yelmox_rtip
