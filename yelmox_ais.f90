

program yelmox

    use ncio 
    use yelmo 
    
    ! External libraries
    use sealevel 
    use isostasy  
    
    use snapclim
    use marine_shelf 
    use smbpal   
    use basal_hydrology 
    use sediments 
    use geothermal
    
    use hyster 

    implicit none 

    type(yelmo_class)      :: yelmo1 
    
    type(sealevel_class)   :: sealev 
    type(snapclim_class)   :: snp1 
    type(marshelf_class)   :: mshlf1 
    type(smbpal_class)     :: smbpal1 
    type(hydro_class)      :: hyd1  
    type(sediments_class)  :: sed1 
    type(geothermal_class) :: gthrm1
    type(isos_class)       :: isos1
    
    type(hyster_class)     :: hyst1 

    character(len=256) :: outfldr, file1D, file2D, file_restart, domain
    character(len=256) :: fAPIS1D, fWAIS1D, fEAIS1D
    character(len=512) :: path_par  
    real(prec) :: time_init, time_end, time_equil, time, dtt, dt1D_out, dt2D_out   
    integer    :: n
    logical    :: calc_transient_climate
    real(4) :: cpu_start_time, cpu_end_time 

    logical, parameter :: use_hyster = .FALSE. 
    real(4) :: conv_km3_Gt, var 
    real(4) :: dTa 

    ! Start timing 
    call cpu_time(cpu_start_time)

    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)

    ! Timing and other parameters 
    call nml_read(path_par,"control","time_init",    time_init)                 ! [yr] Starting time
    call nml_read(path_par,"control","time_end",     time_end)                  ! [yr] Ending time
    call nml_read(path_par,"control","time_equil",   time_equil)                ! [yr] Years to equilibrate first
    call nml_read(path_par,"control","dtt",          dtt)                       ! [yr] Main loop time step 
    call nml_read(path_par,"control","dt1D_out",     dt1D_out)                  ! [yr] Frequency of 1D output 
    call nml_read(path_par,"control","dt2D_out",     dt2D_out)                  ! [yr] Frequency of 2D output 
    call nml_read(path_par,"control","transient",    calc_transient_climate)    ! Calculate transient climate? 

    ! Assume program is running from the output folder
    outfldr = "./"

    ! Define input and output locations 
    file1D     = trim(outfldr)//"yelmo1D.nc"
    file2D     = trim(outfldr)//"yelmo2D.nc"
    ! jablasco
    fWAIS1D = trim(outfldr)//"WAIS1D.nc"
    fEAIS1D = trim(outfldr)//"EAIS1D.nc"
    fAPIS1D = trim(outfldr)//"APIS1D.nc"

    ! === Initialize ice sheet model =====

    ! Initialize data objects and load initial topography
    call yelmo_init(yelmo1,filename=path_par,grid_def="file",time=time_init)

    ! === Initialize external models (forcing for ice sheet) ======

    ! Store domain name as a shortcut 
    domain = yelmo1%par%domain 

    ! Initialize global sea level model (bnd%z_sl)
    call sealevel_init(sealev,path_par)

    ! Initialize bedrock model 
    call isos_init(isos1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,yelmo1%grd%dx)

    ! Initialize "climate" model (climate and ocean forcing)
    call snapclim_init(snp1,path_par,domain,yelmo1%par%grid_name,yelmo1%grd%nx,yelmo1%grd%ny)

    ! Initialize hysteresis module for transient forcing experiments 
    call hyster_init(hyst1,path_par,time_init) 
    conv_km3_Gt = rho_ice *1e-3

    ! Initialize surface mass balance model (bnd%smb, bnd%T_srf)
    call smbpal_init(smbpal1,path_par,x=yelmo1%grd%xc,y=yelmo1%grd%yc,lats=yelmo1%grd%lat)
    
    ! Initialize marine melt model (bnd%bmb_shlf)
    call marshelf_init(mshlf1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,yelmo1%par%grid_name,yelmo1%bnd%basins)
    
    ! Initialize basal hydrology (bnd%H_w)
    call hydro_init(hyd1,path_par,yelmo1%grd%nx,yelmo1%grd%ny)

    ! Load other constant boundary variables (bnd%H_sed, bnd%Q_geo)
    call sediments_init(sed1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,yelmo1%par%grid_name)
    call geothermal_init(gthrm1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,yelmo1%par%grid_name)
    ! === Update initial boundary conditions for current time and yelmo state =====
    ! ybound: z_bed, z_sl, H_sed, H_w, smb, T_srf, bmb_shlf , Q_geo



!mmr    ! This was using as reference H_ref and z_bed the initial H_ref and z_bed, which 
!mmr    ! is only correct when initializing from PD values; for a restart these were not
!mmr    ! set, z_bed_ref was zero and the bedrock bounced instantaneously

!mmr    ! Initialize isostasy - note this should use present-day topography values to 
!mmr    ! calibrate the reference rebound   
!mmr    call isos_init_state(isos1,z_bed=yelmo1%bnd%z_bed,z_bed_ref=yelmo1%bnd%z_bed, &
!mmr          H_ice_ref=yelmo1%tpo%now%H_ice,z_sl=yelmo1%bnd%z_sl*0.0,time=time_init)

!mmr   ! Initialize isostasy using present-day topography values to 
!mmr   ! calibrate the reference rebound
       call isos_init_state(isos1,z_bed=yelmo1%bnd%z_bed,z_bed_ref=yelmo1%bnd%z_bed_ref, &                !mmr
                               H_ice_ref=yelmo1%bnd%H_ice_ref,z_sl=yelmo1%bnd%z_sl*0.0,time=time_init)    !mmr


    call sealevel_update(sealev,year_bp=time_init)
    yelmo1%bnd%z_sl  = sealev%z_sl 
    yelmo1%bnd%H_sed = sed1%now%H 
    
    call hydro_init_state(hyd1,yelmo1%tpo%now%H_ice,time=time_init)

    if (use_hyster) then
        ! snapclim call using anomaly from the hyster package 
        call hyster_calc_forcing(hyst1,time=time,var=yelmo1%reg%V_ice*conv_km3_Gt)
        call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=time_init,domain=domain,dTa=hyst1%f_now)
    else
        ! Normal snapclim call 
        call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=time_init,domain=domain)

    end if 

    ! Equilibrate snowpack for itm
    if (trim(smbpal1%par%abl_method) .eq. "itm") then 
        call smbpal_update_monthly_equil(smbpal1,snp1%now%tas,snp1%now%pr, &
                               yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time_init,time_equil=100.0)
    end if 

! Testing related to present-day surface mass balance
!     snp1%now%tas = snp1%now%tas + 0.0 
!     snp1%now%pr  = snp1%now%pr*exp(0.05*(0.0))

!     call smbpal_update_monthly(smbpal1,snp1%now%tas,snp1%now%pr, &
!                                yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time_init, &
!                                file_out=trim(outfldr)//"smbpal.nc",write_now=.TRUE.,write_init=.TRUE.) 

!     stop 
    call smbpal_update_monthly(smbpal1,snp1%now%tas,snp1%now%pr, &
                               yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time_init) 
    yelmo1%bnd%smb   = smbpal1%ann%smb*conv_we_ie*1e-3    ! [mm we/a] => [m ie/a]
    yelmo1%bnd%T_srf = smbpal1%ann%tsrf 

!     yelmo1%bnd%smb   = yelmo1%dta%pd%smb
!     yelmo1%bnd%T_srf = yelmo1%dta%pd%t2m
    
    call marshelf_calc_Tshlf(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                         yelmo1%bnd%z_sl,depth=snp1%now%depth,to_ann=snp1%now%to_ann, &
                         dto_ann=snp1%now%to_ann - snp1%clim0%to_ann)

    call marshelf_update(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                         yelmo1%bnd%z_sl,regions=yelmo1%bnd%regions,dx=yelmo1%grd%dx)

    yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
    yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  

    yelmo1%bnd%Q_geo    = gthrm1%now%ghf 
    
    call yelmo_print_bound(yelmo1%bnd)

    time = time_init 

    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base)
    call yelmo_init_state(yelmo1,path_par,time=time_init,thrm_method="robin-cold")

    if (yelmo1%par%use_restart) then 
        ! If using restart file, set boundary module variables equal to restarted value 

        hyd1%now%H_water = yelmo1%bnd%H_w 
        isos1%now%z_bed  = yelmo1%bnd%z_bed

    end if 

    ! Run yelmo for several years with constant boundary conditions and topo
    ! to equilibrate thermodynamics and dynamics
    call yelmo_update_equil(yelmo1,time,time_tot=time_equil,topo_fixed=.FALSE.,dt=1.0,ssa_vel_max=500.0)
    
    ! 2D file 
    call yelmo_write_init(yelmo1,file2D,time_init=time,units="years")
!     call write_step_2D(yelmo1,file2D,time=time)  
!     call write_step_2D_small(yelmo1,file2D,time=time)  
    call write_step_2D_combined(yelmo1,isos1,snp1,mshlf1,smbpal1,hyd1,file2D,time=time)
    
    ! 1D file 
    call write_yreg_init(yelmo1,file1D,time_init=time_init,units="years",mask=yelmo1%bnd%ice_allowed)
    call write_yreg_step(yelmo1%reg,file1D,time=time)
    call write_yreg_init(yelmo1,fAPIS1D,time_init=time_init,units="years",mask=yelmo1%bnd%ice_apis)
    call calc_write_yreg_step(yelmo1,fAPIS1D,time=time,mask=yelmo1%bnd%ice_apis)
    call write_yreg_init(yelmo1,fWAIS1D,time_init=time_init,units="years",mask=yelmo1%bnd%ice_wais)
    call calc_write_yreg_step(yelmo1,fWAIS1D,time=time,mask=yelmo1%bnd%ice_wais)
    call write_yreg_init(yelmo1,fEAIS1D,time_init=time_init,units="years",mask=yelmo1%bnd%ice_eais)
    call calc_write_yreg_step(yelmo1,fEAIS1D,time=time,mask=yelmo1%bnd%ice_eais)
 

!     ! Testing smb ======
!     stop "Done."

!     ! Testing hyster ======
!     var = 0.0 
!     do n = 1, 200
!         time = time_init + n*dtt 
!         var  = var + max(0.0,real(100-n,prec)*dtt)
!         if (use_hyster) call hyster_calc_forcing(hyst1,time=time,var=var)
!         write(*,*) "hyst: ", time, hyst1%dv_dt, hyst1%df_dt, hyst1%f_now 
!     end do 

!     stop 

    ! Advance timesteps
    do n = 1, ceiling((time_end-time_init)/dtt)

        ! Get current time 
        time = time_init + n*dtt

        ! == SEA LEVEL ==========================================================
        call sealevel_update(sealev,year_bp=time)
        yelmo1%bnd%z_sl  = sealev%z_sl 

        ! == Yelmo ice sheet ===================================================
        call yelmo_update(yelmo1,time)

        ! == ISOSTASY ==========================================================
        call isos_update(isos1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_sl,time) 
        yelmo1%bnd%z_bed = isos1%now%z_bed


if (calc_transient_climate) then 
        ! == CLIMATE (ATMOSPHERE AND OCEAN) ====================================
        if (mod(time,dtt)==0) then !mmr - this gives problems with restart when dtt is small if (mod(time,2.0)==0) then
            if (use_hyster) then
                ! snapclim call using anomaly from the hyster package 
                call hyster_calc_forcing(hyst1,time=time,var=yelmo1%reg%V_ice*conv_km3_Gt)
                write(*,*) "hyst: ", time, hyst1%time(hyst1%n)-hyst1%time(hyst1%n-1), &
                                                        hyst1%dv_dt, hyst1%df_dt*1e6, hyst1%f_now 
                call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=time_init,domain=domain,dTa=hyst1%f_now)
            else
                ! Normal snapclim call 
                call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=time,domain=domain)
            end if 
        end if 

        ! == SURFACE MASS BALANCE ==============================================

        call smbpal_update_monthly(smbpal1,snp1%now%tas,snp1%now%pr, &
                                   yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time) 
        yelmo1%bnd%smb   = smbpal1%ann%smb*conv_we_ie*1e-3       ! [mm we/a] => [m ie/a]
        yelmo1%bnd%T_srf = smbpal1%ann%tsrf 

!         yelmo1%bnd%smb   = yelmo1%dta%pd%smb
!         yelmo1%bnd%T_srf = yelmo1%dta%pd%t2m
    
        ! == MARINE AND TOTAL BASAL MASS BALANCE ===============================
        call marshelf_calc_Tshlf(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                         yelmo1%bnd%z_sl,depth=snp1%now%depth,to_ann=snp1%now%to_ann, &
                         dto_ann=snp1%now%to_ann - snp1%clim0%to_ann)

        call marshelf_update(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                         yelmo1%bnd%z_sl,regions=yelmo1%bnd%regions,dx=yelmo1%grd%dx*1e-3)

end if 

        yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
        yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  

        ! == BASAL HYDROLOGY ====================================================
        call hydro_update(hyd1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%z_srf,yelmo1%bnd%z_sl, &
                          yelmo1%tpo%now%bmb,yelmo1%tpo%now%f_grnd,yelmo1%tpo%par%dx,yelmo1%tpo%par%dx,time)
        yelmo1%bnd%H_w   = hyd1%now%H_water

        ! == MODEL OUTPUT =======================================================

        if (mod(time,dt2D_out)==0) then 
!             call write_step_2D(yelmo1,file2D,time=time)
!             call write_step_2D_small(yelmo1,file2D,time=time)
            call write_step_2D_combined(yelmo1,isos1,snp1,mshlf1,smbpal1,hyd1,file2D,time=time)
        end if 

        if (mod(time,dt1D_out)==0) then 
            call write_yreg_step(yelmo1%reg,file1D,time=time)
            call calc_write_yreg_step(yelmo1,fAPIS1D,time=time,mask=yelmo1%bnd%ice_apis)
            call calc_write_yreg_step(yelmo1,fWAIS1D,time=time,mask=yelmo1%bnd%ice_wais)
            call calc_write_yreg_step(yelmo1,fEAIS1D,time=time,mask=yelmo1%bnd%ice_eais)
        end if 

        if (mod(time,10.0)==0) then
            write(*,"(a,f14.4)") "yelmo::       time = ", time
        end if 

        if (use_hyster .and. hyst1%kill) then 
            write(*,*) "hyster kill switch activated, f_now = ", hyst1%f_now 
            stop 
        end if 

    end do 

    ! Write a restart file for the end of the simulation
    file_restart = trim(outfldr)//"yelmo_restart.nc"          
    call yelmo_restart_write(yelmo1,file_restart,time=time) 

!     ! Let's see if we can read a restart file 
!     call yelmo_restart_read(yelmo1,file_restart,time=time)

    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call cpu_time(cpu_end_time)

    write(*,"(a,f12.3,a)") "Time  = ",(cpu_end_time-cpu_start_time)/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(time_end-time_init))/((cpu_end_time-cpu_start_time)/3600.0), "kiloyears / hr"

contains

    subroutine write_step_2D_small(ylmo,filename,time)

        implicit none 
        
        type(yelmo_class),      intent(IN) :: ylmo
!         type(snapclim_class),   intent(IN) :: snp 
!         type(marshelf_class),   intent(IN) :: mshlf 
!         type(smbpal_class),     intent(IN) :: srf 
!         type(hydro_class),      intent(IN) :: hyd  
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

    subroutine write_step_2D_combined(ylmo,isos,snp,mshlf,srf,hyd,filename,time)

        implicit none 
        
        type(yelmo_class),      intent(IN) :: ylmo
        type(isos_class),       intent(IN) :: isos 
        type(snapclim_class),   intent(IN) :: snp 
        type(marshelf_class),   intent(IN) :: mshlf 
        type(smbpal_class),     intent(IN) :: srf 
        type(hydro_class),      intent(IN) :: hyd  
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

        ! Write model speed 
        call nc_write(filename,"speed",ylmo%par%model_speed,units="kyr/hr",long_name="Model speed (Yelmo only)", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        
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

        call nc_write(filename,"dist_grline",ylmo%tpo%now%dist_grline,units="km", &
                      long_name="Distance to nearest grounding-line point", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"C_bed",ylmo%dyn%now%C_bed,units="Pa",long_name="Bed friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m-1",long_name="Basal friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff",ylmo%dyn%now%visc_eff,units="Pa a m",long_name="Effective viscosity (SSA)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically integrated velocity magnitude", &
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
        call nc_write(filename,"H_w",ylmo%bnd%H_w,units="m",long_name="Basal water layer thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb",ylmo%bnd%smb,units="m/a ice equiv.",long_name="Surface mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb_errpd",ylmo%bnd%smb-ylmo%dta%pd%smb,units="m/a ice equiv.",long_name="Surface mass balance error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_srf",ylmo%bnd%T_srf,units="K",long_name="Surface temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/a ice equiv.",long_name="Basal mass balance (shelf)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW/m^2",long_name="Geothermal heat flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a ice equiv.",long_name="Basal mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! External data
        call nc_write(filename,"dzbdt",isos%now%dzbdt,units="m/a",long_name="Bedrock uplift rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
 
        call nc_write(filename,"dT_shlf",mshlf%now%dT_shlf,units="K",long_name="Shelf temperature anomaly", &
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

        ! Ice thickness comparison with present-day 
        call nc_write(filename,"H_ice_errpd",ylmo%tpo%now%H_ice-ylmo%dta%pd%H_ice,units="m",long_name="Ice thickness error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s_errpd",ylmo%dyn%now%uxy_s-ylmo%dta%pd%uxy_s,units="m",long_name="Ice thickness error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ssa_mask_acx",ylmo%dyn%now%ssa_mask_acx,units="1",long_name="SSA mask (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_mask_acy",ylmo%dyn%now%ssa_mask_acy,units="1",long_name="SSA mask (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dzsdx",ylmo%tpo%now%dzsdx,units="m/m",long_name="Surface slope", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dzsdy",ylmo%tpo%now%dzsdy,units="m/m",long_name="Surface slope", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"f_grnd_acx",ylmo%tpo%now%f_grnd_acx,units="1",long_name="Grounded fraction (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_grnd_acy",ylmo%tpo%now%f_grnd_acy,units="1",long_name="Grounded fraction (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"taub_acx",ylmo%dyn%now%taub_acx,units="Pa",long_name="Basal stress (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taub_acy",ylmo%dyn%now%taub_acy,units="Pa",long_name="Basal stress (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"taud_acx",ylmo%dyn%now%taud_acx,units="Pa",long_name="Driving stress (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud_acy",ylmo%dyn%now%taud_acy,units="Pa",long_name="Driving stress (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux_b",ylmo%dyn%now%ux_b,units="m/a",long_name="Basal sliding velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_b",ylmo%dyn%now%uy_b,units="m/a",long_name="Basal sliding velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ux_i_bar",ylmo%dyn%now%ux_i_bar,units="m/a",long_name="Internal shear velocity (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_i_bar",ylmo%dyn%now%uy_i_bar,units="m/a",long_name="Internal shear velocity (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ux_bar",ylmo%dyn%now%ux_bar,units="m/a",long_name="Vertically averaged velocity (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_bar",ylmo%dyn%now%uy_bar,units="m/a",long_name="Vertically averaged velocity (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"beta_acx",ylmo%dyn%now%beta_acx,units="Pa a m-1",long_name="Dragging coefficient (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta_acy",ylmo%dyn%now%beta_acy,units="Pa a m-1",long_name="Dragging coefficient (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D_combined

    subroutine write_step_2D(ylmo,filename,time)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: time

        ! Local variables
        integer    :: ncid, n, i, j, nx, ny  
        real(prec) :: time_prev 

        nx = ylmo%tpo%par%nx 
        ny = ylmo%tpo%par%ny 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"dzsdx",ylmo%tpo%now%dzsdx,units="m/m",long_name="Surface slope", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dzsdy",ylmo%tpo%now%dzsdy,units="m/m",long_name="Surface slope", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"dzsrfdt",ylmo%tpo%now%dzsrfdt,units="m/a",long_name="Surface elevation change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dHicedt",ylmo%tpo%now%dHicedt,units="m/a",long_name="Ice thickness change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"f_grnd",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",ylmo%tpo%now%f_ice,units="1",long_name="Ice-covered fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="bar",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"calv",ylmo%tpo%now%calv,units="m/a",long_name="Calving rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_thermodynamics ==
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
        
        ! == yelmo_material ==
        call nc_write(filename,"visc",ylmo%mat%now%visc,units="Pa a",long_name="Viscosity", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"visc_int",ylmo%mat%now%visc_int,units="Pa a m",long_name="Vertically integerated viscosity", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        !call nc_write(filename,"ATT",ylmo%mat%now%ATT,units="a^-1 Pa^-3",long_name="Rate factor", &
        !              dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"enh",ylmo%mat%now%enh,units="1",long_name="Enhancement factor", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"enh_bar",ylmo%mat%now%enh_bar,units="1",long_name="Vertically averaged enhancement factor", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_dynamics ==

        call nc_write(filename,"ssa_mask_acx",ylmo%dyn%now%ssa_mask_acx,units="1",long_name="SSA mask (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_mask_acy",ylmo%dyn%now%ssa_mask_acy,units="1",long_name="SSA mask (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"C_bed",ylmo%dyn%now%C_bed,units="Pa",long_name="Bed friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m^-1",long_name="Basal friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff",ylmo%dyn%now%visc_eff,units="Pa a m",long_name="Effective viscosity (SSA)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"sigma_horiz_sq",ylmo%dyn%now%sigma_horiz_sq,units="1",long_name="Horizontal stress components squared", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"lhs_x",ylmo%dyn%now%lhs_x,units="Pa",long_name="Shear reduction (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"lhs_y",ylmo%dyn%now%lhs_y,units="Pa",long_name="Shear reduction (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"lhs_xy",ylmo%dyn%now%lhs_xy,units="Pa",long_name="Shear reduction magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"duxdz",ylmo%dyn%now%duxdz,units="1/a",long_name="Vertical shear (x)", &
                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"duydz",ylmo%dyn%now%duydz,units="1/a",long_name="Vertical shear (y)", &
                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"taud_acx",ylmo%dyn%now%taud_acx,units="Pa",long_name="Driving stress (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud_acy",ylmo%dyn%now%taud_acy,units="Pa",long_name="Driving stress (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"ux_i_bar",ylmo%dyn%now%ux_i_bar,units="m/a",long_name="Internal shear velocity (x)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_i_bar",ylmo%dyn%now%uy_i_bar,units="m/a",long_name="Internal shear velocity (y)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dd_ab",ylmo%dyn%now%dd_ab,units="m2 a-1",long_name="Diffusivity", &
                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"dd_ab_bar",ylmo%dyn%now%dd_ab_bar,units="m2 a-1",long_name="Diffusivity", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"ux_b",ylmo%dyn%now%ux_b,units="m/a",long_name="Basal sliding velocity (x)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_b",ylmo%dyn%now%uy_b,units="m/a",long_name="Basal sliding velocity (y)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"ux_bar",ylmo%dyn%now%ux_bar,units="m/a",long_name="Vertically integrated velocity (x)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_bar",ylmo%dyn%now%uy_bar,units="m/a",long_name="Vertically integrated velocity (y)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically integrated velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux",ylmo%dyn%now%ux,units="m/a",long_name="Horizontal velocity (x)", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uy",ylmo%dyn%now%uy,units="m/a",long_name="Horizontal velocity (y)", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uxy",ylmo%dyn%now%uxy,units="m/a",long_name="Horizontal velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uz",ylmo%dyn%now%uz,units="m/a",long_name="Vertical velocity", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"f_vbvs",ylmo%dyn%now%f_vbvs,units="1",long_name="Basal to surface velocity fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"f_shear_bar",ylmo%mat%now%f_shear_bar,units="1",long_name="Vertically averaged shearing fraction", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"de",ylmo%mat%now%strn%de,units="",long_name="Vertically averaged strain rate", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

!         call nc_write(filename,"de_bar",ylmo%mat%now%strn2D%de,units="",long_name="Vertically averaged strain rate", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"tau_mx",ylmo%dyn%now%tau_mx,units="Pa",long_name="Basal friction (x)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"tau_my",ylmo%dyn%now%tau_my,units="Pa",long_name="Basal friction (y)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! == yelmo_bound ==
        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"H_sed",ylmo%bnd%H_sed,units="m",long_name="Sediment thickness", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"H_w",ylmo%bnd%H_w,units="m",long_name="Basal water pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"smb",ylmo%bnd%smb,units="m/a",long_name="Surface mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_srf",ylmo%bnd%T_srf,units="K",long_name="Surface temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/a",long_name="Basal mass balance (shelf)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW/m^2",long_name="Geothermal heat flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        

        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a",long_name="Basal mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D

end program yelmox 



