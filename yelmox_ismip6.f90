

program yelmox_ismip6

    use ncio 
    use yelmo 
    use ice_optimization 

    ! External libraries
    use geothermal
    use ismip6
    use isostasy  
    use marine_shelf 
    use sealevel 
    use sediments 
    use smbpal   
    use snapclim
    
    implicit none 

    real(8) :: cpu_start_time, cpu_end_time, cpu_dtime  
    
    character(len=256) :: outfldr, file1D, file2D
    character(len=256) :: file_restart
    character(len=256) :: file_restart_hist
    character(len=256) :: file_restart_trans
    character(len=256) :: domain, grid_name 
    character(len=512) :: path_par, path_const  
    integer    :: n, m
    real(wp)   :: time, time_bp 
    real(wp)   :: time_wt 
    logical    :: init_output_files 

    type(yelmo_class)           :: yelmo1 
    type(sealevel_class)        :: sealev 
    type(snapclim_class)        :: snp1 
    type(marshelf_class)        :: mshlf1 
    type(smbpal_class)          :: smbpal1  
    type(sediments_class)       :: sed1 
    type(geothermal_class)      :: gthrm1
    type(isos_class)            :: isos1
    type(ismip6_forcing_class)  :: ismp1 

    type(snapclim_class)        :: snp2
    type(smbpal_class)          :: smbpal2
    type(marshelf_class)        :: mshlf2 

    type ctrl_params
        character(len=56) :: run_step
        real(wp) :: time_init
        real(wp) :: time_end
        real(wp) :: time_equil      ! Only for spinup
        real(wp) :: time_const      ! Only for spinup 
        real(wp) :: dtt
        real(wp) :: dt1D_out
        real(wp) :: dt2D_out

        logical  :: with_ice_sheet 
        logical  :: optimize
        real(wp) :: time0 
        real(wp) :: time1 

        character(len=256) :: scenario 

    end type 

    type opt_params
        real(wp) :: cf_init
        real(wp) :: cf_min
        real(wp) :: cf_max 
        real(wp) :: tau_c 
        real(wp) :: H0
        real(wp) :: sigma_err 
        real(wp) :: sigma_vel 

        real(wp) :: rel_tau 
        real(wp) :: rel_tau1 
        real(wp) :: rel_tau2
        real(wp) :: rel_time1
        real(wp) :: rel_time2
        real(wp) :: rel_m

        real(wp) :: H_grnd_lim
        real(wp) :: tau_m 
        real(wp) :: m_temp
        real(wp) :: tf_min 
        real(wp) :: tf_max
         
    end type 

    type(ctrl_params)   :: ctl
    type(opt_params)    :: opt 

    ! Set optimize to False by default, unless it is loaded later 
    ctl%optimize = .FALSE. 

    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)

    ! Control parameters 
    call nml_read(path_par,"ctrl","run_step",   ctl%run_step)
    call nml_read(path_par,"ctrl","time0",      ctl%time0)
    call nml_read(path_par,"ctrl","time1",      ctl%time1)
    
    ! Read run_step specific control parameters
    call nml_read(path_par,trim(ctl%run_step),"time_init",  ctl%time_init)      ! [yr] Starting time
    call nml_read(path_par,trim(ctl%run_step),"time_end",   ctl%time_end)       ! [yr] Ending time
    call nml_read(path_par,trim(ctl%run_step),"dtt",        ctl%dtt)            ! [yr] Main loop time step 
    call nml_read(path_par,trim(ctl%run_step),"time_equil", ctl%time_equil)     ! [yr] Years to equilibrate first
    call nml_read(path_par,trim(ctl%run_step),"time_const", ctl%time_const) 
    call nml_read(path_par,trim(ctl%run_step),"dt1D_out",ctl%dt1D_out)          ! [yr] Frequency of 1D output 
    call nml_read(path_par,trim(ctl%run_step),"dt2D_out",ctl%dt2D_out)          ! [yr] Frequency of 2D output 
    
    call nml_read(path_par,trim(ctl%run_step),"with_ice_sheet",ctl%with_ice_sheet)  ! Active ice sheet? 
    call nml_read(path_par,trim(ctl%run_step),"optimize",ctl%optimize)              ! Optimize basal friction?
    
    call nml_read(path_par,trim(ctl%run_step),"scenario",ctl%scenario)              ! Optimize basal friction?
    
    if (ctl%optimize) then 
        ! Load optimization parameters 

        call nml_read(path_par,"opt_L21","cf_init",     opt%cf_init)
        call nml_read(path_par,"opt_L21","cf_min",      opt%cf_min)
        call nml_read(path_par,"opt_L21","cf_max",      opt%cf_max)
        call nml_read(path_par,"opt_L21","tau_c",       opt%tau_c)
        call nml_read(path_par,"opt_L21","H0",          opt%H0)
        call nml_read(path_par,"opt_L21","sigma_err",   opt%sigma_err)   
        call nml_read(path_par,"opt_L21","sigma_vel",   opt%sigma_vel)   
        
        call nml_read(path_par,"opt_L21","rel_tau1",    opt%rel_tau1)   
        call nml_read(path_par,"opt_L21","rel_tau2",    opt%rel_tau2)  
        call nml_read(path_par,"opt_L21","rel_time1",   opt%rel_time1)    
        call nml_read(path_par,"opt_L21","rel_time2",   opt%rel_time2) 
        call nml_read(path_par,"opt_L21","rel_m",       opt%rel_m)

        call nml_read(path_par,"opt_L21","H_grnd_lim",  opt%H_grnd_lim)
        call nml_read(path_par,"opt_L21","tau_m",       opt%tau_m)
        call nml_read(path_par,"opt_L21","m_temp",      opt%m_temp)
        call nml_read(path_par,"opt_L21","tf_min",      opt%tf_min)
        call nml_read(path_par,"opt_L21","tf_max",      opt%tf_max)

    end if 

    ! Set initial time 
    time    = ctl%time_init 
    time_bp = time - 1950.0_wp 

    ! Assume program is running from the output folder
    outfldr = "./"

    ! Define input and output locations 
    path_const          = trim(outfldr)//"yelmo_const_Earth.nml"
    file1D              = trim(outfldr)//"yelmo1D.nc"
    file2D              = trim(outfldr)//"yelmo2D.nc"     
    file_restart        = trim(outfldr)//"yelmo_restart.nc" 
    file_restart_hist   = trim(outfldr)//"yelmo_restart_1950ce.nc"

    if (ctl%time0 .lt. 1e3) then 
        write(file_restart_trans,"(a,i3,a)") trim(outfldr)//"yelmo_restart_", int(ctl%time0), "ce.nc"
    else 
        write(file_restart_trans,"(a,i4,a)") trim(outfldr)//"yelmo_restart_", int(ctl%time0), "ce.nc"
    end if

    !  =========================================================


    ! Print summary of run settings 
    write(*,*)
    write(*,*) "run_step:  ", trim(ctl%run_step) 
    write(*,*)
    write(*,*) "time_init: ",   ctl%time_init 
    write(*,*) "time_end:  ",   ctl%time_end 
    write(*,*) "dtt:       ",   ctl%dtt 
    write(*,*) "dt1D_out:  ",   ctl%dt1D_out 
    write(*,*) "dt2D_out:  ",   ctl%dt2D_out 
    write(*,*) 
    
    if (trim(ctl%run_step) .eq. "spinup" .or. &
        trim(ctl%run_step) .eq. "spinup_ismip6") then 
        write(*,*) "time_equil: ",    ctl%time_equil 
        write(*,*) "time_const: ",    ctl%time_const 

        time_bp = ctl%time_const - 1950.0_wp

    end if 

    write(*,*) "time    = ", time 
    write(*,*) "time_bp = ", time_bp 
    write(*,*) 
    write(*,*) "with_ice_sheet: ",  ctl%with_ice_sheet
    write(*,*) "optimize:       ",  ctl%optimize
    

    ! === Initialize ice sheet model =====

    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)

    ! Initialize data objects and load initial topography
    call yelmo_init(yelmo1,filename=path_par,grid_def="file",time=time)

    ! === Initialize external models (forcing for ice sheet) ======

    ! Store domain and grid_name as shortcuts 
    domain    = yelmo1%par%domain 
    grid_name = yelmo1%par%grid_name 

    ! Initialize global sea level model (bnd%z_sl)
    call sealevel_init(sealev,path_par)

    ! Initialize bedrock model 
    call isos_init(isos1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,yelmo1%grd%dx)

    ! Initialize "climate" model (climate and ocean forcing)
    call snapclim_init(snp1,path_par,domain,grid_name,yelmo1%grd%nx,yelmo1%grd%ny)
    
    ! Initialize surface mass balance model (bnd%smb, bnd%T_srf)
    call smbpal_init(smbpal1,path_par,x=yelmo1%grd%xc,y=yelmo1%grd%yc,lats=yelmo1%grd%lat)
    
    ! Initialize marine melt model (bnd%bmb_shlf)
    call marshelf_init(mshlf1,path_par,"marine_shelf",yelmo1%grd%nx,yelmo1%grd%ny,domain,grid_name,yelmo1%bnd%regions,yelmo1%bnd%basins)
    
    ! Load other constant boundary variables (bnd%H_sed, bnd%Q_geo)
    call sediments_init(sed1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,grid_name)
    call geothermal_init(gthrm1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,grid_name)

    ! Initialize isostasy using present-day topography 
    ! values to calibrate the reference rebound
    call isos_init_state(isos1,z_bed=yelmo1%bnd%z_bed,z_bed_ref=yelmo1%bnd%z_bed_ref, &                !mmr
                           H_ice_ref=yelmo1%bnd%H_ice_ref,z_sl=yelmo1%bnd%z_sl*0.0,time=time)    !mmr


    ! === Update external modules and pass variables to yelmo boundaries =======

    call sealevel_update(sealev,year_bp=time_bp)
    yelmo1%bnd%z_sl  = sealev%z_sl 
    yelmo1%bnd%H_sed = sed1%now%H 

    ! Update snapclim
    call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=time_bp,domain=domain)

    ! Equilibrate snowpack for itm
    if (trim(smbpal1%par%abl_method) .eq. "itm") then 
        call smbpal_update_monthly_equil(smbpal1,snp1%now%tas,snp1%now%pr, &
                               yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time_bp,time_equil=100.0)
    end if 
     
    call smbpal_update_monthly(smbpal1,snp1%now%tas,snp1%now%pr, &
                               yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time_bp) 
    yelmo1%bnd%smb   = smbpal1%ann%smb*conv_we_ie*1e-3    ! [mm we/a] => [m ie/a]
    yelmo1%bnd%T_srf = smbpal1%ann%tsrf 

    call marshelf_update_shelf(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                        yelmo1%bnd%basins,yelmo1%bnd%z_sl,yelmo1%grd%dx,snp1%now%depth, &
                        snp1%now%to_ann,snp1%now%so_ann,dto_ann=snp1%now%to_ann-snp1%clim0%to_ann)

    call marshelf_update(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                         yelmo1%bnd%basins,yelmo1%bnd%z_sl,dx=yelmo1%grd%dx)

    yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
    yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  

    yelmo1%bnd%Q_geo    = gthrm1%now%ghf 
    
    call yelmo_print_bound(yelmo1%bnd)

    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base)
    call yelmo_init_state(yelmo1,time=time,thrm_method="robin-cold")

    if (yelmo1%par%use_restart) then 
        ! If using restart file, set boundary module variables 
        ! equal to restarted value as needed 
         
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

        ! Start timing 
        call yelmo_cpu_time(cpu_start_time)
    
        ! Run yelmo alone for several years with constant boundary conditions and topo
        ! to equilibrate thermodynamics and dynamics
        call yelmo_update_equil(yelmo1,time,time_tot=10.0_wp,     dt=1.0_wp,topo_fixed=.FALSE.)
        call yelmo_update_equil(yelmo1,time,time_tot=ctl%time_equil,dt=1.0_wp,topo_fixed=.TRUE.)

        write(*,*) "Initial equilibration complete."

        ! Initialize output files for checking progress 

        ! 2D file 
        call yelmo_write_init(yelmo1,file2D,time_init=time,units="years")  
        call write_step_2D_combined(yelmo1,isos1,snp1,mshlf1,smbpal1,file2D,time=time)

        ! 1D file 
        call yelmo_write_reg_init(yelmo1,file1D,time_init=time,units="years",mask=yelmo1%bnd%ice_allowed)
        call yelmo_write_reg_step(yelmo1,file1D,time=time) 
        
        ! Next perform 'coupled' model simulations for desired time
        do n = 1, ceiling((ctl%time_end-ctl%time_init)/ctl%dtt)

            ! Get current time 
            time    = ctl%time_init + n*ctl%dtt
            time_bp = ctl%time_const - 1950.0_wp 

            ! == SEA LEVEL ==========================================================
            call sealevel_update(sealev,year_bp=time_bp)
            yelmo1%bnd%z_sl  = sealev%z_sl 

            ! == Yelmo ice sheet ===================================================
            if (ctl%with_ice_sheet) call yelmo_update(yelmo1,time)

            ! == ISOSTASY ==========================================================
            call isos_update(isos1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_sl,time) 
            yelmo1%bnd%z_bed = isos1%now%z_bed


            ! == CLIMATE (ATMOSPHERE AND OCEAN) ====================================
            if (mod(time,ctl%dtt)==0) then
                ! Update snapclim (for elevation changes, but keep time=time_init)
                call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=time_bp,domain=domain)
            end if 

            ! == SURFACE MASS BALANCE ==============================================

            call smbpal_update_monthly(smbpal1,snp1%now%tas,snp1%now%pr, &
                                       yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time_bp) 
            yelmo1%bnd%smb   = smbpal1%ann%smb*conv_we_ie*1e-3       ! [mm we/a] => [m ie/a]
            yelmo1%bnd%T_srf = smbpal1%ann%tsrf 

            ! == MARINE AND TOTAL BASAL MASS BALANCE ===============================
            call marshelf_update_shelf(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                            yelmo1%bnd%basins,yelmo1%bnd%z_sl,yelmo1%grd%dx,snp1%now%depth, &
                            snp1%now%to_ann,snp1%now%so_ann,dto_ann=snp1%now%to_ann-snp1%clim0%to_ann)

            call marshelf_update(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                                 yelmo1%bnd%basins,yelmo1%bnd%z_sl,dx=yelmo1%grd%dx)

            yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
            yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  

            ! == MODEL OUTPUT ===================================

            if (mod(nint(time*100),nint(ctl%dt2D_out*100))==0) then
                call write_step_2D_combined(yelmo1,isos1,snp1,mshlf1,smbpal1,file2D,time)
            end if

            if (mod(nint(time*100),nint(ctl%dt1D_out*100))==0) then
                call yelmo_write_reg_step(yelmo1,file1D,time=time)
                 
            end if 

            if (mod(time,10.0)==0 .and. (.not. yelmo_log)) then
                write(*,"(a,f14.4)") "yelmo:: time = ", time
            end if 
            
        end do 

        write(*,*)
        write(*,*) "Spinup complete."
        write(*,*)

        ! Write the restart file for the end of the simulation
        call yelmo_restart_write(yelmo1,file_restart,time=time_bp) 

        ! Stop timing 
        call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)

        write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
        write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(ctl%time_end-ctl%time_init))/(cpu_dtime/3600.0), " kiloyears / hr"
    
    case("spinup_ismip6")
        ! Model can start from no spinup or equilibration (using restart file), 
        ! here it is run under constant boundary conditions to spinup 
        ! desired state. 

        write(*,*)
        write(*,*) "Performing spinup_ismip6."
        write(*,*) 

        ! Initialize variables inside of ismip6 object 
        call ismip6_forcing_init(ismp1,trim(outfldr)//"/ismip6.nml",gcm="noresm",scen=trim(ctl%scenario), &
                                                domain="Antarctica",grid_name="ANT-32KM")

        ! Update ismip6 forcing to current time
        call ismip6_forcing_update(ismp1,time)

        ! Initialize duplicate climate/smb/mshlf objects for use with ismip data
        
        snp2    = snp1 
        smbpal2 = smbpal1
        mshlf2  = mshlf1
        
        ! Make sure that tf is prescribed externally
        mshlf2%par%tf_method = 0 
        

        ! ===== basal friction optimization ======
        if (ctl%optimize) then 
            
            ! Ensure that cf_ref will be optimized (cb_method == set externally) 
            yelmo1%dyn%par%cb_method = -1  

            ! If not using restart, prescribe cf_ref to initial guess 
            if (.not. yelmo1%par%use_restart) then
                yelmo1%dyn%now%cf_ref = opt%cf_init 
            end if 

            ! Set tf_corr to zero initially 
            mshlf1%now%tf_corr = 0.0_wp 
            mshlf2%now%tf_corr = 0.0_wp 
            
        end if 
        ! ========================================

        ! Start timing 
        call yelmo_cpu_time(cpu_start_time)

if (ctl%with_ice_sheet) then 
        ! Run yelmo alone for several years with constant boundary conditions and topo
        ! to equilibrate thermodynamics and dynamics
        call yelmo_update_equil(yelmo1,time,time_tot=10.0_wp,       dt=1.0_wp,topo_fixed=.FALSE.)
        call yelmo_update_equil(yelmo1,time,time_tot=ctl%time_equil,dt=1.0_wp,topo_fixed=.TRUE.)
end if 

        write(*,*) "Initial equilibration complete."

        ! Initialize output files for checking progress 

        ! 2D file 
        call yelmo_write_init(yelmo1,file2D,time_init=time,units="years")  
        call write_step_2D_combined(yelmo1,isos1,snp1,mshlf1,smbpal1,file2D,time=time)

        ! 1D file 
        call yelmo_write_reg_init(yelmo1,file1D,time_init=time,units="years",mask=yelmo1%bnd%ice_allowed)
        call yelmo_write_reg_step(yelmo1,file1D,time=time) 
        
        ! Next perform 'coupled' model simulations for desired time
        do n = 1, ceiling((ctl%time_end-ctl%time_init)/ctl%dtt)

            ! Get current time 
            time    = ctl%time_init + n*ctl%dtt
            time_bp = ctl%time_const - 1950.0_wp 


            ! ===== basal friction optimization ==================
            if (ctl%optimize) then 

                ! === Optimization parameters =========
                
                ! Update model relaxation time scale and error scaling (in [m])
                opt%rel_tau = get_opt_param(time,time1=opt%rel_time1,time2=opt%rel_time2, &
                                                p1=opt%rel_tau1,p2=opt%rel_tau2,m=opt%rel_m)
                
                ! Set model tau, and set yelmo relaxation switch (2: gl-line and shelves relaxing; 0: no relaxation)
                yelmo1%tpo%par%topo_rel_tau = opt%rel_tau 
                yelmo1%tpo%par%topo_rel     = 2
                if (time .gt. opt%rel_time2) yelmo1%tpo%par%topo_rel = 0 

                ! === Optimization update step =========

                ! Update cf_ref based on error metric(s) 
                call update_cf_ref_errscaling_l21(yelmo1%dyn%now%cf_ref,yelmo1%tpo%now%H_ice, &
                                    yelmo1%tpo%now%dHicedt,yelmo1%bnd%z_bed,yelmo1%bnd%z_sl,yelmo1%dyn%now%ux_s,yelmo1%dyn%now%uy_s, &
                                    yelmo1%dta%pd%H_ice,yelmo1%dta%pd%uxy_s,yelmo1%dta%pd%H_grnd.le.0.0_prec, &
                                    yelmo1%tpo%par%dx,opt%cf_min,opt%cf_max,opt%sigma_err,opt%sigma_vel,opt%tau_c,opt%H0, &
                                    fill_dist=80.0_prec,dt=ctl%dtt)

                ! Update tf_corr based on error metric(s) 
                call update_tf_corr_l21(mshlf2%now%tf_corr,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%H_grnd,yelmo1%tpo%now%dHicedt, &
                                        yelmo1%dta%pd%H_ice,yelmo1%bnd%basins,opt%H_grnd_lim, &
                                        opt%tau_m,opt%m_temp,opt%tf_min,opt%tf_max,dt=ctl%dtt)

            end if 
            ! ====================================================


            ! == SEA LEVEL ==========================================================
            call sealevel_update(sealev,year_bp=time_bp)
            yelmo1%bnd%z_sl  = sealev%z_sl 

            ! == Yelmo ice sheet ===================================================
            if (ctl%with_ice_sheet) call yelmo_update(yelmo1,time)

            ! == ISOSTASY ==========================================================
            call isos_update(isos1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_sl,time) 
            yelmo1%bnd%z_bed = isos1%now%z_bed


            ! Update ismip6 forcing to current time
            call ismip6_forcing_update(ismp1,ctl%time_const, &
                                    use_ref_atm=.TRUE.,use_ref_ocn=.TRUE.)

            ! Set climate to present day 
            snp2%now = snp2%clim0

            ! == SURFACE MASS BALANCE ==============================================

            ! Calculate smb for present day 
            call smbpal_update_monthly(smbpal2,snp2%now%tas,snp2%now%pr, &
                                       yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time) 

            ! Apply ISMIP6 anomalies
            ! (apply to climate too, just for consistency)

            smbpal2%ann%smb  = smbpal2%ann%smb  + ismp1%smb%var(:,:,1,1)*1.0/(conv_we_ie*1e-3) ! [m ie/yr] => [mm we/a]
            smbpal2%ann%tsrf = smbpal2%ann%tsrf + ismp1%ts%var(:,:,1,1)

            do m = 1,12
                snp2%now%tas(:,:,m) = snp2%now%tas(:,:,m) + ismp1%ts%var(:,:,1,1)
                snp2%now%pr(:,:,m)  = snp2%now%pr(:,:,m)  + ismp1%pr%var(:,:,1,1)/365.0 ! [mm/yr] => [mm/d]
            end do 

            snp2%now%ta_ann = sum(snp2%now%tas,dim=3) / 12.0_wp 
            if (trim(domain) .eq. "Antarctica") then 
                snp2%now%ta_sum = sum(snp2%now%tas(:,:,[12,1,2]),dim=3)/3.0     ! Antarctica summer
            else 
                snp2%now%ta_sum = sum(snp2%now%tas(:,:,[6,7,8]),dim=3)/3.0      ! NH summer 
            end if 
            snp2%now%pr_ann = sum(snp2%now%pr,dim=3)  / 12.0 * 365.0     ! [mm/d] => [mm/a]
            
            ! == MARINE AND TOTAL BASAL MASS BALANCE ===============================
            call marshelf_update_shelf(mshlf2,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                            yelmo1%bnd%basins,yelmo1%bnd%z_sl,yelmo1%grd%dx,-ismp1%to%lev, &
                            ismp1%to%var(:,:,:,1),ismp1%so%var(:,:,:,1), &
                            dto_ann=ismp1%to%var(:,:,:,1)-ismp1%to_ref%var(:,:,:,1), &
                            tf_ann=ismp1%tf%var(:,:,:,1))

            ! Update temperature forcing field with tf_corr 
            mshlf2%now%tf_shlf = mshlf2%now%tf_shlf + mshlf2%now%tf_corr
            
            call marshelf_update(mshlf2,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                                 yelmo1%bnd%basins,yelmo1%bnd%z_sl,dx=yelmo1%grd%dx)

            ! Overwrite original mshlf and snp with ismip6 derived ones 
            snp1    = snp2
            smbpal1 = smbpal2  
            mshlf1  = mshlf2 

            yelmo1%bnd%smb      = smbpal1%ann%smb*conv_we_ie*1e-3   ! [mm we/a] => [m ie/a]
            yelmo1%bnd%T_srf    = smbpal1%ann%tsrf 

            yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
            yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  

            ! == MODEL OUTPUT ===================================

            if (mod(nint(time*100),nint(ctl%dt2D_out*100))==0) then
                call write_step_2D_combined(yelmo1,isos1,snp1,mshlf1,smbpal1,file2D,time)
            end if

            if (mod(nint(time*100),nint(ctl%dt1D_out*100))==0) then
                call yelmo_write_reg_step(yelmo1,file1D,time=time)
                 
            end if 

            if (mod(time,10.0)==0 .and. (.not. yelmo_log)) then
                write(*,"(a,f14.4)") "yelmo:: time = ", time
            end if 
            
        end do 

        write(*,*)
        write(*,*) "spinup_ismip6 complete."
        write(*,*)

        ! Write the restart file for the end of the simulation
        call yelmo_restart_write(yelmo1,file_restart,time=time_bp) 

        ! Stop timing 
        call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)

        write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
        write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(ctl%time_end-ctl%time_init))/(cpu_dtime/3600.0), " kiloyears / hr"
        
    case("transient_lgm_to_proj","transient_proj")
        ! Here it is assumed that the model has gone through spinup 
        ! and is ready for transient simulations 

        write(*,*)
        write(*,*) "Performing transient."
        write(*,*) 
 
        ! Initialize variables inside of ismip6 object 
        call ismip6_forcing_init(ismp1,trim(outfldr)//"/ismip6.nml",gcm="noresm",scen=trim(ctl%scenario), &
                                                domain="Antarctica",grid_name="ANT-32KM")

        ! Initialize duplicate climate/smb/mshlf objects for use with ismip data
        
        snp2    = snp1 
        smbpal2 = smbpal1
        mshlf2  = mshlf1
        
        ! Make sure that tf is prescribed externally
        mshlf2%par%tf_method = 0 

        ! Additionally make sure isostasy is update every timestep 
        isos1%par%dt = 1.0_wp 

        ! Start timing 
        call yelmo_cpu_time(cpu_start_time)
        
        ! Get current time 
        time    = ctl%time_init
        time_bp = time - 1950.0_wp 

        init_output_files = .TRUE. 

        ! Perform 'coupled' model simulations for desired time
        do n = 0, ceiling((ctl%time_end-ctl%time_init)/ctl%dtt)

            ! Get current time 
            time    = ctl%time_init + n*ctl%dtt
            time_bp = time - 1950.0_wp 
 
            ! == SEA LEVEL ==========================================================
            call sealevel_update(sealev,year_bp=0.0_wp)
            yelmo1%bnd%z_sl  = sealev%z_sl

            ! == Yelmo ice sheet ===================================================
            if (ctl%with_ice_sheet) call yelmo_update(yelmo1,time)
 
            ! == ISOSTASY ==========================================================
            call isos_update(isos1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_sl,time) 
            yelmo1%bnd%z_bed = isos1%now%z_bed

            if (time .le. ctl%time1) then 

                ! == CLIMATE (ATMOSPHERE AND OCEAN) ====================================

                if (mod(time,ctl%dtt)==0) then !mmr - this gives problems with restart when dtt is small if (mod(time,2.0)==0) then
                    ! Update snapclim
                    call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=time_bp,domain=domain)
                end if 

                ! == SURFACE MASS BALANCE ==============================================

                call smbpal_update_monthly(smbpal1,snp1%now%tas,snp1%now%pr, &
                                           yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time) 
                
                ! == MARINE AND TOTAL BASAL MASS BALANCE ===============================
                call marshelf_update_shelf(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                                yelmo1%bnd%basins,yelmo1%bnd%z_sl,yelmo1%grd%dx,snp1%now%depth, &
                                snp1%now%to_ann,snp1%now%so_ann,dto_ann=snp1%now%to_ann-snp1%clim0%to_ann)

                call marshelf_update(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                                     yelmo1%bnd%basins,yelmo1%bnd%z_sl,dx=yelmo1%grd%dx)

                
            end if 

            if (time .ge. ctl%time0) then 
                ! ISMIP6 forcing 

                ! Update ismip6 forcing to current time
                call ismip6_forcing_update(ismp1,time)

                ! Set climate to present day 
                snp2%now = snp2%clim0

                ! == SURFACE MASS BALANCE ==============================================

                ! Calculate smb for present day 
                call smbpal_update_monthly(smbpal2,snp2%now%tas,snp2%now%pr, &
                                           yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time) 
                
                ! Apply ISMIP6 anomalies
                ! (apply to climate just for consistency)

                smbpal2%ann%smb  = smbpal2%ann%smb  + ismp1%smb%var(:,:,1,1)*1.0/(conv_we_ie*1e-3) ! [m ie/yr] => [mm we/a]
                smbpal2%ann%tsrf = smbpal2%ann%tsrf + ismp1%ts%var(:,:,1,1)

                do m = 1,12
                    snp2%now%tas(:,:,m) = snp2%now%tas(:,:,m) + ismp1%ts%var(:,:,1,1)
                    snp2%now%pr(:,:,m)  = snp2%now%pr(:,:,m)  + ismp1%pr%var(:,:,1,1)/365.0 ! [mm/yr] => [mm/d]
                end do 

                snp2%now%ta_ann = sum(snp2%now%tas,dim=3) / 12.0_wp 
                if (trim(domain) .eq. "Antarctica") then 
                    snp2%now%ta_sum  = sum(snp2%now%tas(:,:,[12,1,2]),dim=3)/3.0  ! Antarctica summer
                else 
                    snp2%now%ta_sum  = sum(snp2%now%tas(:,:,[6,7,8]),dim=3)/3.0  ! NH summer 
                end if 
                snp2%now%pr_ann = sum(snp2%now%pr,dim=3)  / 12.0 * 365.0     ! [mm/d] => [mm/a]
                
                ! == MARINE AND TOTAL BASAL MASS BALANCE ===============================
                call marshelf_update_shelf(mshlf2,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                                yelmo1%bnd%basins,yelmo1%bnd%z_sl,yelmo1%grd%dx,-ismp1%to%lev, &
                                ismp1%to%var(:,:,:,1),ismp1%so%var(:,:,:,1), &
                                dto_ann=ismp1%to%var(:,:,:,1)-ismp1%to_ref%var(:,:,:,1), &
                                tf_ann=ismp1%tf%var(:,:,:,1))

                call marshelf_update(mshlf2,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                                     yelmo1%bnd%basins,yelmo1%bnd%z_sl,dx=yelmo1%grd%dx)

            end if 

            ! Determine which forcing to use based on time period 
            ! LGM to time0 CE      == snapclim 
            ! time0 CE to time1 CE == linear transition from snapclim to ismip6 
            ! time1 CE to future   == ismip6 

            if (time .le. ctl%time0) then 
                ! Only snapclim-based forcing 

                yelmo1%bnd%smb      = smbpal1%ann%smb*conv_we_ie*1e-3       ! [mm we/a] => [m ie/a]
                yelmo1%bnd%T_srf    = smbpal1%ann%tsrf 

                yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
                yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  

            else if (time .gt. ctl%time0 .and. time .le. ctl%time1) then 
                ! Linear-weighted average between snapclim and ismip6 forcing 

                time_wt = (time-ctl%time0) / (ctl%time1 - ctl%time0)

                yelmo1%bnd%smb      = (time_wt*smbpal2%ann%smb  + (1.0-time_wt)*smbpal1%ann%smb)  * conv_we_ie*1e-3
                yelmo1%bnd%T_srf    =  time_wt*smbpal2%ann%tsrf + (1.0-time_wt)*smbpal1%ann%tsrf

                yelmo1%bnd%bmb_shlf = time_wt*mshlf2%now%bmb_shlf + (1.0-time_wt)*mshlf1%now%bmb_shlf 
                yelmo1%bnd%T_shlf   = time_wt*mshlf2%now%T_shlf   + (1.0-time_wt)*mshlf1%now%T_shlf  

            else    ! time > ctl%time1
                ! Only ISMIP6 forcing 

                ! Overwrite original mshlf and snp with ismip6 derived ones 
                snp1    = snp2
                smbpal1 = smbpal2  
                mshlf1  = mshlf2 
                
                yelmo1%bnd%smb      = smbpal2%ann%smb*conv_we_ie*1e-3
                yelmo1%bnd%T_srf    = smbpal2%ann%tsrf 

                yelmo1%bnd%bmb_shlf = mshlf2%now%bmb_shlf  
                yelmo1%bnd%T_shlf   = mshlf2%now%T_shlf  

            
            end if
            
            ! == MODEL OUTPUT ===================================

            if (init_output_files) then 
                ! Initialize output files 
                call yelmo_write_init(yelmo1,file2D,time_init=time,units="years")
                call yelmo_write_reg_init(yelmo1,file1D,time_init=time,units="years",mask=yelmo1%bnd%ice_allowed) 
                init_output_files = .FALSE. 
            end if

            if (mod(nint(time*100),nint(ctl%dt2D_out*100))==0) then
                call write_step_2D_combined(yelmo1,isos1,snp1,mshlf1,smbpal1, &
                                                  file2D,time,snp2,mshlf2,smbpal2)
            end if

            if (mod(nint(time*100),nint(ctl%dt1D_out*100))==0) then
                call yelmo_write_reg_step(yelmo1,file1D,time=time)
                 
            end if 

            if (mod(time,10.0)==0 .and. (.not. yelmo_log)) then
                write(*,"(a,f14.4)") "yelmo:: time = ", time
            end if 
            
            if (time == ctl%time0) then 
                ! Write restart file at start of transition period
                call yelmo_restart_write(yelmo1,file_restart_trans,time=time) 
            end if 

            if (time == 1950.0_wp) then 
                ! Write restart file at start of hist period
                call yelmo_restart_write(yelmo1,file_restart_hist,time=time) 
            end if 

        end do 

        write(*,*)
        write(*,*) "Transient complete."
        write(*,*)

        ! Write the restart file for the end of the transient simulation
        call yelmo_restart_write(yelmo1,file_restart,time=time) 

        ! Stop timing 
        call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)

        write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
        write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(ctl%time_end-ctl%time_init))/(cpu_dtime/3600.0), " kiloyears / hr"
        
    end select

    ! Finalize program
    call yelmo_end(yelmo1,time=time)

contains
    
    subroutine write_step_2D_combined(ylmo,isos,snp,mshlf,srf,filename,time,snp2,mshlf2,srf2)

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

        type(snapclim_class),   intent(IN), optional :: snp2 
        type(marshelf_class),   intent(IN), optional :: mshlf2 
        type(smbpal_class),     intent(IN), optional :: srf2 
        
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
        
        call nc_write(filename,"dHicedt",ylmo%tpo%now%dHicedt,units="m/yr",long_name="Ice thickness rate of change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"mb_applied",ylmo%tpo%now%mb_applied,units="m",long_name="Applied net mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"H_grnd",ylmo%tpo%now%H_grnd,units="m",long_name="Ice thickness overburden", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="bar",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"calv",ylmo%tpo%now%calv,units="m/a ice equiv.",long_name="Calving rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"f_grnd",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",ylmo%tpo%now%f_ice,units="1",long_name="Ice fraction in grid cell", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"dist_grline",ylmo%tpo%now%dist_grline,units="km", &
!                       long_name="Distance to nearest grounding-line point", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"cf_ref",ylmo%dyn%now%cf_ref,units="--",long_name="Bed friction scalar", &
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

!         call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically-averaged velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        !call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
        !              dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        !call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
        !              dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

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

        call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW/m^2",long_name="Geothermal heat flux", &
                      dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
        
        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a ice equiv.",long_name="Basal mass balance", &
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
        call nc_write(filename,"S_shlf",mshlf%now%S_shlf,units="PSU",long_name="Shelf salinity", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dT_shlf",mshlf%now%dT_shlf,units="K",long_name="Shelf temperature anomaly", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_fp_shlf",mshlf%now%T_fp_shlf,units="K",long_name="Shelf freezing temperature", &
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

if (.FALSE. .and. present(snp2) .and. present(mshlf2)) then 
    ! Output values from second set of climate objects 
        call nc_write(filename,"Ta_ann_2",snp2%now%ta_ann,units="K",long_name="Near-surface air temperature (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Ta_sum_2",snp2%now%ta_sum,units="K",long_name="Near-surface air temperature (sum)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Pr_ann_2",snp2%now%pr_ann*1e-3,units="m/a water equiv.",long_name="Precipitation (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"T_shlf_2",mshlf2%now%T_shlf,units="K",long_name="Shelf temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"S_shlf_2",mshlf2%now%S_shlf,units="PSU",long_name="Shelf salinity", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dT_shlf_2",mshlf2%now%dT_shlf,units="K",long_name="Shelf temperature anomaly", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        if (trim(mshlf2%par%bmb_method) .eq. "pico") then 
            call nc_write(filename,"d_shlf_2",mshlf2%pico%now%d_shlf,units="km",long_name="Shelf distance to grounding line", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"d_if_2",mshlf2%pico%now%d_if,units="km",long_name="Shelf distance to ice front", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"boxes_2",mshlf2%pico%now%boxes,units="",long_name="Shelf boxes", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)      
            call nc_write(filename,"r_shlf_2",mshlf2%pico%now%r_shlf,units="",long_name="Ratio of ice shelf", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"T_box_2",mshlf2%pico%now%T_box,units="K?",long_name="Temperature of boxes", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"S_box_2",mshlf2%pico%now%S_box,units="PSU",long_name="Salinity of boxes", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"A_box_2",mshlf2%pico%now%A_box*1e-6,units="km2",long_name="Box area of ice shelf", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        else 
            call nc_write(filename,"tf_basin_2",mshlf2%now%tf_basin,units="K",long_name="Mean basin thermal forcing", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"tf_shlf_2",mshlf2%now%tf_shlf,units="K",long_name="Shelf thermal forcing", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"tf_corr_2",mshlf2%now%tf_corr,units="K",long_name="Shelf thermal forcing applied correction factor", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"slope_base_2",mshlf2%now%slope_base,units="",long_name="Shelf-base slope", &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        end if 

        call nc_write(filename,"dTa_ann_2",snp2%now%ta_ann-snp%clim0%ta_ann,units="K",long_name="Near-surface air temperature anomaly (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dTa_sum_2",snp2%now%ta_sum-snp%clim0%ta_sum,units="K",long_name="Near-surface air temperature anomaly (sum)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dPr_ann_2",(snp2%now%pr_ann-snp%clim0%pr_ann)*1e-3,units="m/a water equiv.",long_name="Precipitation anomaly (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
end if 

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
       
        call nc_write(filename,"ux_s",ylmo%dyn%now%ux_s,units="m/a",long_name="Surface velocity (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_s",ylmo%dyn%now%uy_s,units="m/a",long_name="Surface velocity (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

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
        
                ! Strain-rate and stress tensors 
        if (.TRUE.) then

            ! call nc_write(filename,"de",ylmo%mat%now%strn%de,units="a^-1",long_name="Effective strain rate", &
            !           dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
            ! call nc_write(filename,"te",ylmo%mat%now%strs%te,units="Pa",long_name="Effective stress", &
            !           dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
            ! call nc_write(filename,"visc_int",ylmo%mat%now%visc_int,units="Pa a m",long_name="Depth-integrated effective viscosity (SSA)", &
            !           dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
         
            ! call nc_write(filename,"de2D",ylmo%mat%now%strn2D%de,units="yr^-1",long_name="Effective strain rate", &
            !           dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"div2D",ylmo%mat%now%strn2D%div,units="yr^-1",long_name="Divergence strain rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            ! call nc_write(filename,"te2D",ylmo%mat%now%strs2D%te,units="Pa",long_name="Effective stress", &
            !           dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            
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

end program yelmox_ismip6



