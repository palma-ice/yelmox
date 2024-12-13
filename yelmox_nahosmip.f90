
program yelmox_ismip6

    use nml
    use ncio
    use timer
    use timeout
    use yelmo
    use ice_optimization

    ! External libraries
    use fastisostasy    ! also reexports barysealevel
    use geothermal
    use ismip6
    use marine_shelf
    use sediments
    use smbpal
    use snapclim

    implicit none 

    character(len=256) :: outfldr, file1D, file2D, file2D_small
    character(len=256) :: file1D_ismip6, file2D_ismip6
    character(len=256) :: file_isos, file_bsl
    character(len=256) :: domain, grid_name 
    character(len=512) :: path_par  
    character(len=512) :: path_tf_corr 
    character(len=512) :: ismip6_path_par
    integer  :: n, m
    real(wp) :: time, time_bp
    real(wp) :: time_elapsed
    real(wp) :: time_wt 

    real(sp) :: convert_km3_Gt

    type(yelmo_class)           :: yelmo1
    type(bsl_class)             :: bsl
    type(snapclim_class)        :: snp1
    type(marshelf_class)        :: mshlf1
    type(smbpal_class)          :: smbpal1
    type(sediments_class)       :: sed1
    type(geothermal_class)      :: gthrm1
    type(isos_class)            :: isos1
    type(ismip6_forcing_class)  :: ismp1

    type(timeout_class) :: tm_1D, tm_2D, tm_2Dsm

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

        logical  :: with_ice_sheet 
        character(len=56) :: equil_method

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

    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)

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

    call nml_read(path_par,trim(ctl%run_step),"with_ice_sheet",ctl%with_ice_sheet)  ! Active ice sheet? 
    call nml_read(path_par,trim(ctl%run_step),"equil_method",  ctl%equil_method)    ! What method should be used for spin-up?

    if (trim(ctl%equil_method) .eq. "opt") then 
        ! Load optimization parameters 

        call optimize_par_load(opt,path_par,"opt")

    end if 

    ! Get output times
    call timeout_init(tm_1D,  path_par,"tm_1D",  "small",  ctl%time_init,ctl%time_end)
    call timeout_init(tm_2D,  path_par,"tm_2D",  "heavy",  ctl%time_init,ctl%time_end)
    call timeout_init(tm_2Dsm,path_par,"tm_2Dsm","medium", ctl%time_init,ctl%time_end)
    
    ! Assume program is running from the output folder
    outfldr = "./"

    ! Define input and output locations
    file1D              = trim(outfldr)//"yelmo1D.nc"
    file2D              = trim(outfldr)//"yelmo2D.nc"
    file2D_small        = trim(outfldr)//"yelmo2Dsm.nc"

    file_isos           = trim(outfldr)//"fastisostasy.nc"
    file_bsl            = trim(outfldr)//"bsl.nc"

    file1D_ismip6       = trim(outfldr)//"yelmo1D_ismip6.nc"
    file2D_ismip6       = trim(outfldr)//"yelmo2D_ismip6.nc"

    tmr_file            = trim(outfldr)//"timer_table.txt"

    ! Set initial model time 
    time    = ctl%time_init 
    time_bp = time - 1950.0_wp 

    !  =========================================================
    ! Print summary of run settings 
    write(*,*)
    write(*,*) "run_step:  ", trim(ctl%run_step) 
    write(*,*)
    write(*,*) "time_init: ",   ctl%time_init 
    write(*,*) "time_end:  ",   ctl%time_end 
    write(*,*) "dtt:       ",   ctl%dtt  
    write(*,*) 

! ajr: ismip6
    write(*,*) "ismip6_par_file:        ", trim(ctl%ismip6_par_file)
    write(*,*) "ismip6_expname:         ", trim(ctl%ismip6_expname)
    write(*,*) "ismip6_experiment:      ", trim(ismip6exp%experiment)
    
    select case(trim(ctl%run_step))

        case("spinup")

            write(*,*) "time_equil: ",    ctl%time_equil 
            write(*,*) "time_const: ",    ctl%time_const 

            time_bp = ctl%time_const - 1950.0_wp

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

        case DEFAULT 

            reg1%write = .FALSE.
            reg2%write = .FALSE.
            reg3%write = .FALSE.

    end select


    ! === Initialize external models (forcing for ice sheet) ======

    ! Initialize barysealevel model
    call bsl_init(bsl, path_par, time_bp)

    ! Initialize fastisosaty
    call isos_init(isos1, path_par, "isos", yelmo1%grd%nx, yelmo1%grd%ny, &
        yelmo1%grd%dx, yelmo1%grd%dy)

    ! Initialize "climate" model (climate and ocean forcing)
    call snapclim_init(snp1,path_par,domain,yelmo1%par%grid_name,yelmo1%grd%nx,yelmo1%grd%ny,yelmo1%bnd%basins)
    
    ! Initialize surface mass balance model (bnd%smb, bnd%T_srf)
    call smbpal_init(smbpal1,path_par,x=yelmo1%grd%xc,y=yelmo1%grd%yc,lats=yelmo1%grd%lat)
    
    ! Initialize marine melt model (bnd%bmb_shlf)
    call marshelf_init(mshlf1,path_par,"marine_shelf",yelmo1%grd%nx,yelmo1%grd%ny,domain,grid_name,yelmo1%bnd%regions,yelmo1%bnd%basins)
    
    ! Initialize variables inside of ismip6 object 
    ismip6_path_par = trim(outfldr)//"/"//trim(ctl%ismip6_par_file)
    call nahosmip_ant_forcing_init(ismp1,ismip6_path_par,domain,grid_name,experiment=ismip6exp%experiment)
    
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
    call bsl_update(bsl, year_bp=time_bp)
    call bsl_write_init(bsl, file_bsl, time)

    ! Initialize the isostasy reference state using reference topography fields
    call isos_init_ref(isos1, yelmo1%bnd%z_bed_ref, yelmo1%bnd%H_ice_ref)
    call isos_init_state(isos1, yelmo1%bnd%z_bed, yelmo1%tpo%now%H_ice, time, bsl)
    call isos_write_init_extended(isos1, file_isos, time)

    yelmo1%bnd%z_bed = isos1%out%z_bed
    yelmo1%bnd%z_sl  = isos1%out%z_ss

    ! Update snapclim
    call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=time_bp,domain=domain,dx=yelmo1%grd%dx,basins=yelmo1%bnd%basins)

    ! Equilibrate snowpack for itm
    if (trim(smbpal1%par%abl_method) .eq. "itm") then 
        call smbpal_update_monthly_equil(smbpal1,snp1%now%tas,snp1%now%pr, &
                               yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time_bp,time_equil=100.0)
    end if 
    
    ! Update forcing to present-day reference using ISMIP6 forcing
    call calc_climate_nahosmip(snp1,smbpal1,mshlf1,ismp1,yelmo1, &
                time=ctl%time_const,time_bp=ctl%time_const-1950.0_wp)
    
    yelmo1%bnd%smb      = smbpal1%ann%smb*yelmo1%bnd%c%conv_we_ie*1e-3   ! [mm we/a] => [m ie/a]
    yelmo1%bnd%T_srf    = smbpal1%ann%tsrf 

    yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
    yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  

    call yelmo_print_bound(yelmo1%bnd)

    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base)
    call yelmo_init_state(yelmo1,time=time,thrm_method="robin-cold")
    
    if (yelmo1%par%use_restart) then 
        ! Perform additional startup steps when using a restart

        ! Set boundary module variables equal to restarted value         
        isos1%now%z_bed  = yelmo1%bnd%z_bed
      
    end if


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
        call yelmo_write_init(yelmo1,file2D,time_init=time,units="years")  
        call yelmo_write_reg_init(yelmo1,file1D,time_init=time,units="years",mask=yelmo1%bnd%ice_allowed)
        
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
            
            ! == ISOSTASY and SEA LEVEL ===========================================
            call bsl_update(bsl, time_bp)
            call isos_update(isos1, yelmo1%tpo%now%H_ice, time, bsl, dwdt_corr=yelmo1%bnd%dzbdt_corr)
            yelmo1%bnd%z_bed = isos1%out%z_bed
            yelmo1%bnd%z_sl  = isos1%out%z_ss

            call timer_step(tmrs,comp=1,time_mod=[time-ctl%dtt,time]*1e-3,label="isostasy") 

            ! == ICE SHEET ===================================================
            if (ctl%with_ice_sheet) call yelmo_update(yelmo1,time)

            call timer_step(tmrs,comp=2,time_mod=[time-ctl%dtt,time]*1e-3,label="yelmo") 

            ! == CLIMATE ===========================================================

            ! Update forcing to present-day reference, but 
            ! adjusting to ice topography
            call calc_climate_nahosmip(snp1,smbpal1,mshlf1,ismp1,yelmo1, &
                        time=ctl%time_const,time_bp=ctl%time_const-1950.0_wp)

            yelmo1%bnd%smb      = smbpal1%ann%smb*yelmo1%bnd%c%conv_we_ie*1e-3   ! [mm we/a] => [m ie/a]
            yelmo1%bnd%T_srf    = smbpal1%ann%tsrf 

            yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
            yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  

            call timer_step(tmrs,comp=3,time_mod=[time-ctl%dtt,time]*1e-3,label="climate") 

            ! == MODEL OUTPUT ===================================

            if (timeout_check(tm_2D,time)) then
                call write_step_2D_combined(yelmo1,isos1,snp1,mshlf1,smbpal1,file2D,time)
            end if

            if (timeout_check(tm_1D,time)) then
                call yelmo_write_reg_step(yelmo1,file1D,time=time)
                 
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

        ! Write the restart snapshot for the end of the simulation
        call yelmox_restart_write(bsl,isos1,yelmo1,time_bp)

    case("transient")
        ! Here it is assumed that the model has gone through spinup 
        ! and is ready for transient simulations 

        write(*,*)
        write(*,*) "Performing transient."
        write(*,*) 

        ! Additionally make sure isostasy is updated every timestep 
        isos1%par%dt_prognostics = 1.0_wp 
        isos1%par%dt_diagnostics = 10.0_wp 
        
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

            call timer_step(tmrs,comp=0) 
            
            ! == ISOSTASY and SEA LEVEL ===========================================
            call bsl_update(bsl, time_bp)
            call isos_update(isos1, yelmo1%tpo%now%H_ice, time, bsl, dwdt_corr=yelmo1%bnd%dzbdt_corr)
            yelmo1%bnd%z_bed = isos1%out%z_bed
            yelmo1%bnd%z_sl  = isos1%out%z_ss

            call timer_step(tmrs,comp=1,time_mod=[time-ctl%dtt,time]*1e-3,label="isostasy") 

            ! == ICE SHEET ===================================================
            if (ctl%with_ice_sheet) call yelmo_update(yelmo1,time)

            call timer_step(tmrs,comp=2,time_mod=[time-ctl%dtt,time]*1e-3,label="yelmo") 

            ! == CLIMATE and OCEAN ==========================================

            ! Get ISMIP6 climate and ocean forcing
            call calc_climate_nahosmip(snp1,smbpal1,mshlf1,ismp1,yelmo1,time,time_bp)
            
            yelmo1%bnd%smb      = smbpal1%ann%smb*yelmo1%bnd%c%conv_we_ie*1e-3   ! [mm we/a] => [m ie/a]
            yelmo1%bnd%T_srf    = smbpal1%ann%tsrf 

            yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
            yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf   

            call timer_step(tmrs,comp=3,time_mod=[time-ctl%dtt,time]*1e-3,label="climate") 

            ! == MODEL OUTPUT ===================================

            if (timeout_check(tm_2D,time)) then
                call write_step_2D_combined(yelmo1,isos1,snp1,mshlf1,smbpal1,file2D,time)
            end if
           
             
            if (timeout_check(tm_1D,time)) then
                 call yelmo_write_reg_step(yelmo1,file1D,time=time)
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

        ! Write the restart snapshot for the end of the transient simulation
        call yelmox_restart_write(bsl,isos1,yelmo1,time)
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
        call nc_write(filename,"cmb",ylmo%tpo%now%cmb,units="m/a ice equiv.",long_name="Calving rate mb", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"cmb_flt",ylmo%tpo%now%cmb_flt,units="m/a ice equiv.",long_name="Calving rate mb flt", &
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
        call nc_write(filename,"smb",ylmo%tpo%now%smb,units="m/a ice equiv.",long_name="Net surface mass balance", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb_ref",ylmo%bnd%smb,units="m/a ice equiv.",long_name="Surface mass balance", &
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

        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a ice equiv.",long_name="Net basal mass balance", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"fmb",ylmo%tpo%now%fmb,units="m/a ice equiv.",long_name="Net margin-front mass balance", &
                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                        
        ! External data
        call nc_write(filename,"dzbdt",isos%out%dwdt,units="m/a",long_name="Bedrock uplift rate", &
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

    subroutine calc_climate_nahosmip(snp,smbp,mshlf,ismp,ylmo,time,time_bp,dTa,dTo)

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

        call nahosmip_ant_forcing_update(ismp,time)
        
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
            ! Update temperature fields with extra anomaly 
            mshlf1%now%T_shlf  = mshlf1%now%T_shlf  + dTo
            mshlf1%now%dT_shlf = mshlf1%now%dT_shlf + dTo
            mshlf1%now%tf_shlf = mshlf1%now%tf_shlf + dTo
        end if 

        ! Update bmb_shlf and mask_ocn
        call marshelf_update(mshlf,ylmo%tpo%now%H_ice,ylmo%bnd%z_bed,ylmo%tpo%now%f_grnd, &
                             ylmo%bnd%regions,ylmo%bnd%basins,ylmo%bnd%z_sl,dx=ylmo%grd%dx)

        return

    end subroutine calc_climate_nahosmip

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

    subroutine yelmox_restart_write(bsl,isos,ylmo,time,fldr)

        implicit none

        type(bsl_class),    intent(IN) :: bsl
        type(isos_class),   intent(IN) :: isos
        type(yelmo_class),  intent(IN) :: ylmo
        real(wp),           intent(IN) :: time 
        character(len=*),   intent(IN), optional :: fldr

        ! Local variables
        real(wp) :: time_kyr
        character(len=32)   :: time_str
        character(len=1024) :: outfldr

        character(len=56), parameter :: file_bsl   = "bsl_restart.nc"
        character(len=56), parameter :: file_isos  = "isos_restart.nc"
        character(len=56), parameter :: file_yelmo = "yelmo_restart.nc"
        
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

    end subroutine yelmox_restart_write

end program yelmox_ismip6
