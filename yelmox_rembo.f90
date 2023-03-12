

program yelmox

    use ncio 
    use yelmo 
    use ice_optimization 

    ! External libraries
    use sealevel 
    use isostasy  
    
    use rembo_sclimate 
    use snapclim
    use marine_shelf 
    use sediments 
    use geothermal
    
    use hyster 

    implicit none 

    type(yelmo_class)      :: yelmo1 
    type(snapclim_class)   :: snp1 
    type(sealevel_class)   :: sealev 
    type(marshelf_class)   :: mshlf1 
    type(sediments_class)  :: sed1 
    type(geothermal_class) :: gthrm1
    type(isos_class)       :: isos1
    
    type(hyster_class)     :: hyst1 

    character(len=256) :: outfldr, file1D, file2D, file1D_hyst, file_restart, domain 
    character(len=256) :: file_rembo
    character(len=512) :: path_par, path_const  
    real(wp) :: time_init, time_end, time_equil, time, dtt, dt1D_out, dt2D_out, dt_restart   
    real(wp) :: time_elapsed
    integer  :: n
    logical  :: calc_transient_climate
    
    logical :: use_hyster
    logical :: write_restart 
    real(4) :: convert_km3_Gt, var 
    real(4) :: dTa, dT_summer

    ! No-ice mask (to impose additional melting)
    logical, allocatable :: mask_noice(:,:)  
    logical :: lim_pd_ice 
    logical :: with_ice_sheet 
    logical :: optimize 

    type(ice_opt_params)  :: opt 

    real(8) :: cpu_start_time, cpu_end_time, cpu_dtime  
    
    ! Start timing 
    call yelmo_cpu_time(cpu_start_time)
    
    ! Assume program is running from the output folder
    outfldr = "./"
    
    ! Determine the parameter file from the command line 
    !call yelmo_load_command_line_args(path_par)
    path_par = trim(outfldr)//"yelmo_Greenland_rembo.nml"
    
    ! Timing and other parameters 
    call nml_read(path_par,"ctrl","time_init",      time_init)              ! [yr] Starting time
    call nml_read(path_par,"ctrl","time_end",       time_end)               ! [yr] Ending time
    call nml_read(path_par,"ctrl","time_equil",     time_equil)             ! [yr] Years to equilibrate first
    call nml_read(path_par,"ctrl","dtt",            dtt)                    ! [yr] Main loop time step 
    call nml_read(path_par,"ctrl","dt1D_out",       dt1D_out)               ! [yr] Frequency of 1D output 
    call nml_read(path_par,"ctrl","dt2D_out",       dt2D_out)               ! [yr] Frequency of 2D output 
    call nml_read(path_par,"ctrl","write_restart",  write_restart)
    call nml_read(path_par,"ctrl","transient",      calc_transient_climate) ! Calculate transient climate? 
    call nml_read(path_par,"ctrl","use_hyster",     use_hyster)             ! Use hyster?
    call nml_read(path_par,"ctrl","dT",             dT_summer)              ! Initial summer temperature anomaly
    call nml_read(path_par,"ctrl","lim_pd_ice",     lim_pd_ice)             ! Limit to pd ice extent (apply extra melting outside mask)
    call nml_read(path_par,"ctrl","with_ice_sheet", with_ice_sheet)         ! Active ice sheet? 
    call nml_read(path_par,"ctrl","optimize",       optimize)               ! Optimize basal friction?
    
    if (optimize) then 
        ! Load optimization parameters 

        call optimize_par_load(opt,path_par,"opt")

    end if 
    
    ! Define input and output locations 
    path_const   = trim(outfldr)//"yelmo_const_Earth.nml"
    file1D       = trim(outfldr)//"yelmo1D.nc"
    file2D       = trim(outfldr)//"yelmo2D.nc"
    file_restart = trim(outfldr)//"yelmo_restart.nc"          
    file1D_hyst  = trim(outfldr)//"yelmo1D_hyst.nc" 

    file_rembo   = trim(outfldr)//"yelmo-rembo.nc"

    ! How often to write a restart file 
    dt_restart   = 20e3                 ! [yr] 


    time = time_init 

    ! === Initialize ice sheet model =====

    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)

    ! Initialize data objects and load initial topography
    call yelmo_init(yelmo1,filename=path_par,grid_def="file",time=time_init)

    ! === Initialize external models (forcing for ice sheet) ======

    ! Store domain name as a shortcut 
    domain = yelmo1%par%domain 



    ! Ensure optimization fields are allocated and preassigned
    allocate(opt%cf_min(yelmo1%grd%nx,yelmo1%grd%ny))
    allocate(opt%cf_max(yelmo1%grd%nx,yelmo1%grd%ny))
    
    opt%cf_min = opt%cf_min_par 
    opt%cf_max = yelmo1%dyn%par%till_cf_ref

    ! Define no-ice mask from present-day data
    allocate(mask_noice(yelmo1%grd%nx,yelmo1%grd%ny))
    mask_noice = .FALSE. 
    where(yelmo1%dta%pd%H_ice .le. 0.0) mask_noice = .TRUE. 



    
    ! Initialize global sea level model (bnd%z_sl)
    call sealevel_init(sealev,path_par)

    ! Initialize bedrock model 
    call isos_init(isos1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,yelmo1%grd%dx)

    ! Initialize the climate model REMBO, including loading parameters from options_rembo 
    call rembo_init(real(time_init,8))

    ! Initialize hysteresis module for transient forcing experiments 
    call hyster_init(hyst1,path_par,time_init) 
    convert_km3_Gt = yelmo1%bnd%c%rho_ice *1e-3

    ! Initialize "climate" model (here for ocean forcing)
    call snapclim_init(snp1,path_par,domain,yelmo1%par%grid_name,yelmo1%grd%nx,yelmo1%grd%ny,yelmo1%bnd%basins)
    
    ! Initialize marine melt model (bnd%bmb_shlf)
    call marshelf_init(mshlf1,path_par,"marine_shelf",yelmo1%grd%nx,yelmo1%grd%ny, &
                        domain,yelmo1%par%grid_name,yelmo1%bnd%regions,yelmo1%bnd%basins)
    
    ! Load other constant boundary variables (bnd%H_sed, bnd%Q_geo)
    call sediments_init(sed1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,yelmo1%par%grid_name)
    call geothermal_init(gthrm1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,yelmo1%par%grid_name)
    ! === Update initial boundary conditions for current time and yelmo state =====
    ! ybound: z_bed, z_sl, H_sed, H_w, smb, T_srf, bmb_shlf , Q_geo

    ! Initialize isostasy using present-day topography 
    ! values to calibrate the reference rebound
    call isos_init_state(isos1,z_bed=yelmo1%bnd%z_bed,H_ice=yelmo1%tpo%now%H_ice, &
                                    z_sl=yelmo1%bnd%z_sl,z_bed_ref=yelmo1%bnd%z_bed_ref, &
                                    H_ice_ref=yelmo1%bnd%H_ice_ref, &
                                    z_sl_ref=yelmo1%bnd%z_sl*0.0,time=time)
                                    
    call sealevel_update(sealev,year_bp=time_init)
    yelmo1%bnd%z_sl  = sealev%z_sl 
    yelmo1%bnd%H_sed = sed1%now%H 
    
    if (use_hyster) then
        ! Update hysteresis variable 
        call hyster_calc_forcing(hyst1,time=time,var=yelmo1%reg%V_ice*convert_km3_Gt)
        dT_summer = hyst1%f_now 
    end if 

    ! Update REMBO, with correct topography, let it equilibrate for several years 
    ! do n = 1, 100
    !     time = time_init + real(n,8)    
    !     call rembo_update(real(time,8),real(dT_summer,8),real(yelmo1%tpo%now%z_srf,8),real(yelmo1%tpo%now%H_ice,8))
    ! end do 
    ! rembo_ann%time_emb = time_init 
    ! rembo_ann%time_smb = time_init  
    call rembo_update(real(time_init,8),real(dT_summer,8),real(yelmo1%tpo%now%z_srf,8),real(yelmo1%tpo%now%H_ice,8))
    
    ! Update surface mass balance and surface temperature from REMBO
    yelmo1%bnd%smb   = rembo_ann%smb    *yelmo1%bnd%c%conv_we_ie*1e-3       ! [mm we/a] => [m ie/a]
    yelmo1%bnd%T_srf = rembo_ann%T_srf
    
    ! Special treatment for Greenland
    if (lim_pd_ice) then 
        
        ! Impose additional negative mass balance to no ice points of 2 [m.i.e./a] melting
        where(mask_noice) yelmo1%bnd%smb = yelmo1%bnd%smb - 4.0 

    end if 
    
    ! Update snapclim
    call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=time_init,domain=domain,dx=yelmo1%grd%dx,basins=yelmo1%bnd%basins)

    call marshelf_update_shelf(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                        yelmo1%bnd%basins,yelmo1%bnd%z_sl,yelmo1%grd%dx,snp1%now%depth, &
                        snp1%now%to_ann,snp1%now%so_ann,dto_ann=snp1%now%to_ann-snp1%clim0%to_ann)

    call marshelf_update(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                         yelmo1%bnd%regions,yelmo1%bnd%basins,yelmo1%bnd%z_sl,dx=yelmo1%grd%dx)

    yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
    yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  

    yelmo1%bnd%Q_geo    = gthrm1%now%ghf 
    
    call yelmo_print_bound(yelmo1%bnd)

    time = time_init 

    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base)
    call yelmo_init_state(yelmo1,time=time_init,thrm_method="robin-cold")

    ! ===== basal friction optimization ======
    if (optimize) then 
        
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
    
    if (yelmo1%par%use_restart) then 
        ! If using restart file, set boundary module variables equal to restarted value 

        ! Set boundary module variables equal to restarted value         
        isos1%now%z_bed  = yelmo1%bnd%z_bed

    else 
        ! Run yelmo for several years with constant boundary conditions and topo
        ! to equilibrate thermodynamics and dynamics

        if (with_ice_sheet) then 
            call yelmo_update_equil(yelmo1,time,time_tot=10.0_prec, dt=1.0_prec,topo_fixed=.FALSE.)
            call yelmo_update_equil(yelmo1,time,time_tot=time_equil,dt=dtt,topo_fixed=.TRUE.)
        end if 

    end if 

    ! 2D file 
    if (dt2D_out .ne. 0.0) then 
        call yelmo_write_init(yelmo1,file2D,time_init=time,units="years")
        call write_step_2D_combined(yelmo1,rembo_ann,isos1,mshlf1,file2D,time=time)
    end if 

    ! 2D small file 
    ! call yelmo_write_init(yelmo1,file2D_small,time_init=time,units="years")
    ! call write_step_2D_combined_small(yelmo1,isos1,snp1,mshlf1,smbpal1,file2D_small,time=time)
    
    ! 1D file 
    ! call yelmo_write_reg_init(yelmo1,file1D,time_init=time,units="years",mask=yelmo1%bnd%ice_allowed)
    ! call yelmo_write_reg_step(yelmo1,file1D,time=time)

    call write_yelmo_init_combined(yelmo1,file_rembo,time_init=time,units="years", &
                    mask=yelmo1%bnd%ice_allowed,dT_min=hyst1%par%f_min,dT_max=hyst1%par%f_max)
    call write_step_2D_combined_small(yelmo1,hyst1,rembo_ann,isos1,mshlf1,file_rembo,time)

    ! Advance timesteps
    do n = 1, ceiling((time_end-time_init)/dtt)

        ! Get current time 
        time         = time_init + n*dtt
        time_elapsed = time - time_init

        ! == SEA LEVEL ==========================================================
        call sealevel_update(sealev,year_bp=time)
        yelmo1%bnd%z_sl  = sealev%z_sl 

        ! == Yelmo ice sheet ===================================================
        if (with_ice_sheet) then

            if (optimize) then 

                ! === Optimization update step =========

                if (opt%opt_cf .and. &
                        (time_elapsed .ge. opt%cf_time_init .and. time_elapsed .le. opt%cf_time_end) ) then  
                    ! Perform cf_ref optimization
                
                    ! Update cb_ref based on error metric(s) 
                    call optimize_cb_ref(yelmo1%dyn%now%cb_ref,yelmo1%tpo%now%H_ice, &
                                        yelmo1%tpo%now%dHidt,yelmo1%bnd%z_bed,yelmo1%bnd%z_sl,yelmo1%dyn%now%ux_s,yelmo1%dyn%now%uy_s, &
                                        yelmo1%dta%pd%H_ice,yelmo1%dta%pd%uxy_s,yelmo1%dta%pd%H_grnd, &
                                        opt%cf_min,opt%cf_max,yelmo1%tpo%par%dx,opt%sigma_err,opt%sigma_vel,opt%tau_c,opt%H0, &
                                        dt=dtt,fill_method=opt%fill_method,fill_dist=opt%sigma_err, &
                                        cb_tgt=yelmo1%dyn%now%cb_tgt)

                end if

                if (opt%opt_tf .and. &
                        (time_elapsed .ge. opt%tf_time_init .and. time_elapsed .le. opt%tf_time_end) ) then
                    ! Perform tf_corr optimization

                    call optimize_tf_corr(mshlf1%now%tf_corr,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%H_grnd,yelmo1%tpo%now%dHidt, &
                                            yelmo1%dta%pd%H_ice,yelmo1%dta%pd%H_grnd,opt%H_grnd_lim,opt%tau_m,opt%m_temp, &
                                            opt%tf_min,opt%tf_max,yelmo1%tpo%par%dx,sigma=opt%tf_sigma,dt=dtt)
                    ! call optimize_tf_corr_basin(mshlf1%now%tf_corr,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%H_grnd,yelmo1%tpo%now%dHidt, &
                    !                         yelmo1%dta%pd%H_ice,yelmo1%bnd%basins,opt%H_grnd_lim, &
                    !                         opt%tau_m,opt%m_temp,opt%tf_min,opt%tf_max,opt%tf_basins,dt=dtt)
                
                end if 

            
            end if 

            ! Update ice sheet to current time 
            call yelmo_update(yelmo1,time)

        end if 

        ! == ISOSTASY ==========================================================
        call isos_update(isos1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_sl,time,yelmo1%bnd%dzbdt_corr) 
        yelmo1%bnd%z_bed = isos1%now%z_bed


if (calc_transient_climate) then 
        ! == CLIMATE (ATMOSPHERE AND OCEAN) ====================================
        
        if (use_hyster) then
            ! Update forcing based on hysteresis module
            call hyster_calc_forcing(hyst1,time=time,var=yelmo1%reg%V_ice*convert_km3_Gt)
            write(*,*) "hyst: ", time, hyst1%dt, hyst1%dv_dt, hyst1%df_dt*1e6, hyst1%f_now 
            
            dT_summer = hyst1%f_now 
        end if 

        ! call REMBO1     
        call rembo_update(real(time,8),real(dT_summer,8),real(yelmo1%tpo%now%z_srf,8),real(yelmo1%tpo%now%H_ice,8))
        
        ! Update surface mass balance and surface temperature from REMBO
        yelmo1%bnd%smb   = rembo_ann%smb    *yelmo1%bnd%c%conv_we_ie*1e-3       ! [mm we/a] => [m ie/a]
        yelmo1%bnd%T_srf = rembo_ann%T_srf
         
        ! Special treatment for Greenland
        if (lim_pd_ice) then 
        
            ! Impose additional negative mass balance to no ice points of 4 [m.i.e./a] melting
            where(mask_noice) yelmo1%bnd%smb = yelmo1%bnd%smb - 4.0 

        end if 

        ! == MARINE AND TOTAL BASAL MASS BALANCE ===============================
        
        ! Update snapclim
        call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=time,domain=domain,dx=yelmo1%grd%dx,basins=yelmo1%bnd%basins) 

        call marshelf_update_shelf(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                        yelmo1%bnd%basins,yelmo1%bnd%z_sl,yelmo1%grd%dx,snp1%now%depth, &
                        snp1%now%to_ann,snp1%now%so_ann,dto_ann=snp1%now%to_ann-snp1%clim0%to_ann)

        call marshelf_update(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                             yelmo1%bnd%regions,yelmo1%bnd%basins,yelmo1%bnd%z_sl,dx=yelmo1%grd%dx)

end if 

        yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
        yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  

        ! == MODEL OUTPUT =======================================================

        if (dt2D_out .ne. 0.0 .and. mod(time,dt2D_out)==0) then 
            call write_step_2D_combined(yelmo1,rembo_ann,isos1,mshlf1,file2D,time=time)
        end if 

        if (mod(time,dt1D_out)==0) then 
            ! call yelmo_write_reg_step(yelmo1,file1D,time=time) 
            !call write_step_1D_combined(yelmo1,hyst1,file1D_hyst,time=time)
            call write_step_2D_combined_small(yelmo1,hyst1,rembo_ann,isos1,mshlf1,file_rembo,time)
        end if 

        if (write_restart .and. mod(time,dt_restart)==0) then 
            call yelmo_restart_write(yelmo1,file_restart,time=time) 
        end if 

        if (mod(time,10.0)==0) then
            write(*,"(a,f14.4)") "yelmo::       time = ", time
        end if 

        if (use_hyster .and. hyst1%kill) then 
            write(*,"(a,f12.3,a,f12.3)") "hyster:: kill switch activated. [time, f_now] = ", &
                                                        time, ", ", hyst1%f_now 
            exit  
        end if 

    end do 

    ! Write the restart file for the end of the simulation
    if (write_restart) then 
        call yelmo_restart_write(yelmo1,file_restart,time=time) 
    end if

!     ! Let's see if we can read a restart file 
!     call yelmo_restart_read(yelmo1,file_restart,time=time)

    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)

    write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(time-time_init))/(cpu_dtime/3600.0), " kiloyears / hr"

contains
    
    subroutine write_step_2D_combined_small(ylmo,hyst,rembo,isos,mshlf,filename,time)

        implicit none 
        
        type(yelmo_class),    intent(IN) :: ylmo
        type(hyster_class),   intent(IN) :: hyst
        type(rembo_class),    intent(IN) :: rembo
        type(isos_class),     intent(IN) :: isos 
        type(marshelf_class), intent(IN) :: mshlf 
        character(len=*),     intent(IN) :: filename
        real(wp), intent(IN) :: time

        ! Local variables
        integer  :: ncid, n, k
        real(wp) :: time_prev 
        real(wp) :: npmb, ntot, aar, smb_tot  
        real(wp) :: dT_axis(1000)
        type(yregions_class) :: reg

        ! Assume region to write is the global region of yelmo 
        reg = ylmo%reg 

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
        call nc_write(filename,"hyst_df_dt",hyst%df_dt*1e6,units="K/(1e6 yr)",long_name="hyst: forcing rate of change", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"hyst_dv_dt",hyst%dv_dt,units="Gt/yr",long_name="hyst: volume rate of change", &
                      dim1="time",start=[n],ncid=ncid)

        ! Write volume in volume-dT phase space
        call nc_read(filename,"dT_axis",dT_axis) 
        k = minloc(abs(dT_axis-hyst%f_now),dim=1)
        call nc_write(filename,"V_dT",reg%V_ice*1e-6,units="1e6 km^3",long_name="Ice volume", &
                      dim1="time",start=[k],ncid=ncid)
        
        ! == yelmo metrics ==

        ! Write model metrics (model speed, dt, eta)
        call yelmo_write_step_model_metrics(filename,ylmo,n,ncid)

        ! Write present-day data metrics (rmse[H],etc)
        ! call yelmo_write_step_pd_metrics(filename,ylmo,n,ncid)
        
        ! == 1D Variables ==
        
        call nc_write(filename,"V_ice",reg%V_ice*1e-6,units="1e6 km^3",long_name="Ice volume", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"A_ice",reg%A_ice*1e-6,units="1e6 km^2",long_name="Ice area", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"V_sle",reg%V_sle,units="m sle",long_name="Sea-level equivalent volume", &
                      dim1="time",start=[n],ncid=ncid)
        
        ! == yelmo_topography ==

        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/yr",long_name="Surface velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! == rembo climate == 
        call nc_write(filename,"ta_ann",rembo%T_ann,units="K",long_name="REMBO Near-surface air temperature (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ta_sum",rembo%T_jja,units="K",long_name="REMBO Near-surface air temperature (sum)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"pr_ann",rembo%pr*1e-3,units="m/a water equiv.",long_name="REMBO Precipitation (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"smb",rembo%smb*1e-3,units="m/yr water equiv.",long_name="REMBO Surface mass balance (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        

        ! == ice-sheet wide metrics == 

        ! Get integrated metrics (smb_tot [Gt/yr] and aar [unitless])
        ntot = count(ylmo%tpo%now%H_ice .gt. 0.0)

        if (ntot .gt. 0.0) then 
            npmb = count(ylmo%tpo%now%H_ice .gt. 0.0 .and. rembo%smb .gt. 0.0)
            aar  = real(npmb,prec) / real(ntot,prec) 

            smb_tot = (ylmo%tpo%par%dx**2)*sum(rembo%smb*1e-3,mask=ylmo%tpo%now%H_ice .gt. 0.0)

            ! Convert from m^3/yr => Gt/yr
            ! [m^3/yr] * [1000 kg/m^3] * [1e-12 Gt/kg] == [Gt/yr]
            smb_tot = smb_tot * (1000) *1e-12 
        else
            aar = 0.0
        end if  

        call nc_write(filename,"dT_jja",hyst%f_now,units="K",long_name="hyst: forcing value", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"smb_mean",smb_tot,units="Gt/yr",long_name="Mean smb over the ice sheet", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"aar",aar,units="1",long_name="Accumulation area ratio", &
                      dim1="time",start=[n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D_combined_small

    subroutine write_step_2D_combined(ylmo,rembo,isos,mshlf,filename,time)

        implicit none 
        
        type(yelmo_class),    intent(IN) :: ylmo
        type(rembo_class),    intent(IN) :: rembo
        type(isos_class),     intent(IN) :: isos 
        type(marshelf_class), intent(IN) :: mshlf 
        character(len=*),     intent(IN) :: filename
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


        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
        !               dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! call nc_write(filename,"Q_strn",ylmo%thrm%now%Q_strn/(ylmo%bnd%c%rho_ice*ylmo%thrm%now%cp),units="K a-1",long_name="Strain heating", &
        !               dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="J a-1 m-2",long_name="Basal frictional heating", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"bmb_grnd",ylmo%thrm%now%bmb_grnd,units="m/a ice equiv.",long_name="Basal mass balance (grounded)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_w",ylmo%thrm%now%H_w,units="m",long_name="Basal water layer thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ATT",ylmo%mat%now%ATT,units="a^-1 Pa^-3",long_name="Rate factor", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

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
        ! call nc_write(filename,"smb_errpd",ylmo%bnd%smb-ylmo%dta%pd%smb,units="m/a ice equiv.",long_name="Surface mass balance error wrt present day", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
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
 
        call nc_write(filename,"dT_shlf",mshlf%now%dT_shlf,units="K",long_name="Shelf temperature anomaly", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Comparison with present-day 
        call nc_write(filename,"H_ice_pd_err",ylmo%dta%pd%err_H_ice,units="m",long_name="Ice thickness error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf_pd_err",ylmo%dta%pd%err_z_srf,units="m",long_name="Surface elevation error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s_pd_err",ylmo%dta%pd%err_uxy_s,units="m/a",long_name="Surface velocity error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        

        ! === REMBO ANNUAL FIELDS ===

        call nc_write(filename,"Ta_ann",rembo%T_ann,units="K",long_name="REMBO Near-surface air temperature (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Ta_sum",rembo%T_jja,units="K",long_name="REMBO Near-surface air temperature (sum)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"pr_ann",rembo%pr*1e-3,units="m/a water equiv.",long_name="REMBO Precipitation (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"smb_ann",rembo%smb*1e-3,units="m/a water equiv.",long_name="REMBO Surface mass balance (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! call nc_write(filename,"dTa_ann",snp%now%ta_ann-snp%clim0%ta_ann,units="K",long_name="Near-surface air temperature anomaly (ann)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"dTa_sum",snp%now%ta_sum-snp%clim0%ta_sum,units="K",long_name="Near-surface air temperature anomaly (sum)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"dPr_ann",(snp%now%pr_ann-snp%clim0%pr_ann)*1e-3,units="m/a water equiv.",long_name="Precipitation anomaly (ann)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D_combined

    subroutine write_step_1D_combined(ylm,hyst,filename,time)

        implicit none 
        
        type(yelmo_class),  intent(IN) :: ylm
        type(hyster_class), intent(IN) :: hyst 
        character(len=*),   intent(IN) :: filename
        real(wp),           intent(IN) :: time

        ! Local variables
        integer  :: ncid, n, k 
        real(wp) :: time_prev 
        real(wp) :: dT_axis(1000) 
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
        call nc_write(filename,"hyst_df_dt",hyst%df_dt*1e6,units="K/(1e6 yr)",long_name="hyst: forcing rate of change", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"hyst_dv_dt",hyst%dv_dt,units="Gt/yr",long_name="hyst: volume rate of change", &
                      dim1="time",start=[n],ncid=ncid)

        ! Write volume in volume-dT phase space
        call nc_read(filename,"dT_axis",dT_axis) 
        k = minloc(abs(dT_axis-hyst%f_now),dim=1)
        call nc_write(filename,"V_dT",reg%V_ice*1e-6,units="1e6 km^3",long_name="Ice volume", &
                      dim1="time",start=[k],ncid=ncid)
        
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

    end subroutine write_step_1D_combined

    subroutine write_yelmo_init_combined(dom,filename,time_init,units,mask,dT_min,dT_max)

        implicit none 

        type(yelmo_class), intent(IN) :: dom 
        character(len=*),  intent(IN) :: filename, units 
        real(wp),          intent(IN) :: time_init
        logical,           intent(IN) :: mask(:,:) 
        real(wp),          intent(IN) :: dT_min 
        real(wp),          intent(IN) :: dT_max

        ! Local variables
        real(wp), allocatable :: dT_axis(:) 

        ! Initialize netcdf file and dimensions
        call nc_create(filename)
        call nc_write_dim(filename,"xc",        x=dom%grd%xc*1e-3,      units="kilometers")
        call nc_write_dim(filename,"yc",        x=dom%grd%yc*1e-3,      units="kilometers")
        call nc_write_dim(filename,"zeta",      x=dom%par%zeta_aa,      units="1")
        call nc_write_dim(filename,"time",      x=time_init,dx=1.0_wp,nx=1,units=trim(units),unlimited=.TRUE.)
        
        !============================================================
        ! Add temperature axis to 1D hysteresis file 
        allocate(dT_axis(1000))
        do n = 1, 1000 
            dT_axis(n) = dT_min + (dT_max-dT_min)*(n-1)/real(1000-1,prec)
        end do 
        call nc_write_dim(filename,"dT_axis",x=dT_axis,units="degC")
        
        ! Populate variable with missing values for now to initialize it
        dT_axis = missing_value 
        call nc_write(filename,"V_dT",dT_axis,dim1="dT_axis")
        !============================================================
        
        ! Static information
        call nc_write(filename,"mask", mask,  units="1",long_name="Region mask",dim1="xc",dim2="yc")
        
        return

    end subroutine write_yelmo_init_combined

end program yelmox



