

program yelmox

    use ncio 
    use yelmo 
    use yelmo_boundaries
    use basal_dragging 

    ! External libraries
    use sealevel 
    use isostasy  
    
    use snapclim
    use marine_shelf 
    use smbpal 
    use sediments 
    use geothermal
    
    use ice_optimization 
    use ice_enhancement 

    implicit none 

    type(yelmo_class)      :: yelmo1 
    
    type(sealevel_class)   :: sealev 
    type(ice_enh_class)    :: ice_enh  
    type(snapclim_class)   :: snp1 
    type(marshelf_class)   :: mshlf1 
    type(smbpal_class)     :: smbpal1  
    type(sediments_class)  :: sed1 
    type(geothermal_class) :: gthrm1
    type(isos_class)       :: isos1
    
    character(len=256) :: outfldr, file1D, file2D, file_restart_init, file_restart, domain 
    character(len=256) :: file_cf_ref 
    character(len=512) :: path_par, path_const  
    real(prec) :: time_init, time_end, time_equil, time, dtt, dt1D_out, dt2D_out, dt_restart 
    real(prec) :: time_init_final   
    integer    :: k, n
    logical    :: calc_transient_climate, load_cf_ref 
    real(prec) :: f_cf, f_cf_lim, f_hol, dtas_hol, dpr_hol, dsmb_negis   
    logical    :: write2D_now 

    real(prec), allocatable :: dsmb_now(:,:) 
    real(prec), allocatable :: dtas_now(:,:) 
    real(prec), allocatable :: dpr_now(:,:) 

    logical, parameter :: check_init     = .FALSE. 
    logical, parameter :: calc_ice_sheet = .TRUE. 

! === optimization ===
    
    logical                 :: optimize_cf
    logical                 :: overwrite_files
    real(prec)              :: err_scale 
    real(prec)              :: fill_dist 
    real(prec), allocatable :: cf_ref_dot(:,:) 
    integer                 :: n_iter 
    character(len=12)       :: iter_str
    integer                 :: iter 

    type(yelmo_class)       :: yelmo0
    type(smbpal_class)      :: smbpal0  
    type(isos_class)        :: isos0
    
! ====================
    
    real(8) :: cpu_start_time, cpu_end_time, cpu_dtime  
    
    ! Start timing 
    call yelmo_cpu_time(cpu_start_time)
    
    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)

    ! Timing and other parameters 
    call nml_read(path_par,"ctrl","time_init",    time_init)                 ! [yr] Starting time
    call nml_read(path_par,"ctrl","time_init_final",time_init_final)         ! [yr] Starting time for final iteration after optimization
    call nml_read(path_par,"ctrl","time_end",     time_end)                  ! [yr] Ending time
    call nml_read(path_par,"ctrl","time_equil",   time_equil)                ! [yr] Years to equilibrate first
    call nml_read(path_par,"ctrl","dtt",          dtt)                       ! [yr] Main loop time step 
    call nml_read(path_par,"ctrl","dt1D_out",     dt1D_out)                  ! [yr] Frequency of 1D output 
    call nml_read(path_par,"ctrl","dt2D_out",     dt2D_out)                  ! [yr] Frequency of 2D output 
    call nml_read(path_par,"ctrl","transient",    calc_transient_climate)    ! Calculate transient climate? 

    call nml_read(path_par,"ctrl","load_cf_ref",  load_cf_ref)               ! Load cf_ref from file? Otherwise define from cf_stream + inline tuning
    call nml_read(path_par,"ctrl","file_cf_ref",  file_cf_ref)               ! Filename holding cf_ref to load 
    call nml_read(path_par,"ctrl","optimize_cf",  optimize_cf)               ! Run optimization iterations on cf_ref?
    call nml_read(path_par,"ctrl","overwrite_files",overwrite_files)         ! Overwrite files from optimization iterations?
    call nml_read(path_par,"ctrl","n_iter",       n_iter)                    ! Number of optimization iterations
    
    call nml_read(path_par,"ctrl","f_cf",         f_cf)                      ! Friction scaling parameter 
    call nml_read(path_par,"ctrl","f_cf_lim",     f_cf_lim)                  ! Friction scaling parameter limit
    call nml_read(path_par,"ctrl","f_hol",        f_hol)                     ! Holocene index scaling parameter 
    call nml_read(path_par,"ctrl","dtas_hol",     dtas_hol)                  ! Anomaly to apply to default climate during the Holocene
    call nml_read(path_par,"ctrl","dpr_hol",      dpr_hol)                   ! Anomaly to apply to default climate during the Holocene
    call nml_read(path_par,"ctrl","dsmb_negis",   dsmb_negis)                ! Anomaly to apply to default climate during the Holocene
    call nml_read(path_par,"ctrl","fill_dist",    fill_dist)                 ! [km] Distance to fill cf_ref with minimum value instead of neighborhood mean

    ! Assume program is running from the output folder
    outfldr = "./"

    ! Define input and output locations 
    path_const   = trim(outfldr)//"yelmo_const_Earth.nml"
    file1D       = trim(outfldr)//"yelmo1D.nc"
    file2D       = trim(outfldr)//"yelmo2D.nc"
    file_restart = trim(outfldr)//"yelmo_restart.nc"          
    file_restart_init = trim(outfldr)//"yelmo_restart_init.nc"          
    
    ! How often to write a restart file 
    dt_restart   = 20e3                 ! [yr] 


    if (check_init) time_init = -7e3

    ! === opt ======

    ! Optimization ice-thickness scaling parameter (higher value=>slower but more stable optimization)
    err_scale = 1000.0 

    ! If not optimizing set iterations to 1, and eliminate optimization related choices
    if (.not. optimize_cf) then 
        n_iter = 1 
        time_init_final = time_init 
        overwrite_files = .TRUE. 
    end if 

    ! === opt ======

    ! === Initialize ice sheet model =====

    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)

    ! Initialize data objects and load initial topography
    call yelmo_init(yelmo1,filename=path_par,grid_def="file",time=time_init)

    ! Allocate smb, temp and precip anomaly field for writing 
    allocate(dsmb_now(yelmo1%grd%nx,yelmo1%grd%ny))
    dsmb_now = 0.0_prec 

    allocate(dtas_now(yelmo1%grd%nx,yelmo1%grd%ny))
    dtas_now = 0.0_prec 

    allocate(dpr_now(yelmo1%grd%nx,yelmo1%grd%ny))
    dpr_now = 0.0_prec 

    ! Define cf_ref_dot for later use 
    allocate(cf_ref_dot(yelmo1%grd%nx,yelmo1%grd%ny))
    cf_ref_dot = 0.0 
    
    ! === Initialize external models (forcing for ice sheet) ======

    ! Store domain name as a shortcut 
    domain = yelmo1%par%domain 

    ! Initialize global sea level model (bnd%z_sl)
    call sealevel_init(sealev,path_par)

    ! Initialize paleo enhancement factor time series (input to bnd%enh_srf)
    call ice_enh_init(ice_enh,path_par)

    ! Initialize bedrock model 
    call isos_init(isos1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,yelmo1%grd%dx)

    ! Initialize "climate" model (climate and ocean forcing)
    call snapclim_init(snp1,path_par,domain,yelmo1%par%grid_name,yelmo1%grd%nx,yelmo1%grd%ny)

    ! Ensure that global f_hol is used if running with 1ind_new method 
    if (trim(snp1%par%atm_type) .ne. "hybrid") snp1%hybrid%f_hol = f_hol 

    ! Initialize surface mass balance model (bnd%smb, bnd%T_srf)
    call smbpal_init(smbpal1,path_par,x=yelmo1%grd%xc,y=yelmo1%grd%yc,lats=yelmo1%grd%lat)
    
    ! Initialize marine melt model (bnd%bmb_shlf)
    call marshelf_init(mshlf1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,yelmo1%par%grid_name,yelmo1%bnd%regions,yelmo1%bnd%basins)
    
    ! Load other constant boundary variables (bnd%H_sed, bnd%Q_geo)
    call sediments_init(sed1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,yelmo1%par%grid_name)
    call geothermal_init(gthrm1,path_par,yelmo1%grd%nx,yelmo1%grd%ny,domain,yelmo1%par%grid_name)
    ! === Update initial boundary conditions for current time and yelmo state =====
    ! ybound: z_bed, z_sl, H_sed, H_w, smb, T_srf, bmb_shlf , Q_geo

    ! Initialize isostasy - note this should use present-day topography values to 
    ! calibrate the reference rebound!
    call isos_init_state(isos1,z_bed=yelmo1%bnd%z_bed,z_bed_ref=yelmo1%bnd%z_bed_ref, &
                               H_ice_ref=yelmo1%bnd%H_ice_ref,z_sl=yelmo1%bnd%z_sl*0.0,time=time_init)

    call sealevel_update(sealev,year_bp=time_init)
    yelmo1%bnd%z_sl  = sealev%z_sl 
    yelmo1%bnd%H_sed = sed1%now%H 
    
    ! Get initial boundary enhancement factor value
    call ice_enh_update(ice_enh,time=time_init)
    yelmo1%bnd%enh_srf = ice_enh%enh 

    ! Normal snapclim call 
    call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=time_init,domain=domain)
    
    ! Modify tas and precip
    call modify_tas(snp1%now%tas,snp1%now%ta_ann,snp1%now%ta_sum,dtas_now,dtas_hol,yelmo1%grd,time_init)
    call modify_pr(snp1%now%pr,snp1%now%pr_ann,dpr_now,dpr_hol,yelmo1%grd,time_init)
    
    ! Equilibrate snowpack for itm
    if (trim(smbpal1%par%abl_method) .eq. "itm") then 
        call smbpal_update_monthly_equil(smbpal1,snp1%now%tas,snp1%now%pr, &
                               yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time_init,time_equil=100.0)
    end if 

    ! Update snowpack again 
    call smbpal_update_monthly(smbpal1,snp1%now%tas,snp1%now%pr, &
                               yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time_init) 
    yelmo1%bnd%smb   = smbpal1%ann%smb*conv_we_ie*1e-3    ! [mm we/a] => [m ie/a]
    yelmo1%bnd%T_srf = smbpal1%ann%tsrf 

    ! Impose flux correction to smb 
    call modify_smb(yelmo1%bnd%smb,dsmb_now,dsmb_negis,yelmo1%bnd,yelmo1%grd,time_init)

!     yelmo1%bnd%smb   = yelmo1%dta%pd%smb
!     yelmo1%bnd%T_srf = yelmo1%dta%pd%t2m
    
    call marshelf_calc_Tshlf(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                         yelmo1%bnd%z_sl,depth=snp1%now%depth,to_ann=snp1%now%to_ann, &
                         dto_ann=snp1%now%to_ann - snp1%clim0%to_ann)

    call marshelf_update(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                         yelmo1%bnd%z_sl,dx=yelmo1%grd%dx)

    yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
    yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  
    yelmo1%bnd%Q_geo    = gthrm1%now%ghf 
    
    call yelmo_print_bound(yelmo1%bnd)

    time = time_init 
    
    ! ============================================================
    ! Load or define cf_ref (this will be overwritten if loading from restart file)

    if (yelmo1%dyn%par%cb_method .eq. -1) then 
    
        if (load_cf_ref) then 

            ! Parse filename with grid information
            call yelmo_parse_path(file_cf_ref,yelmo1%par%domain,yelmo1%par%grid_name)

            ! Load cf_ref from specified file 
            call nc_read(file_cf_ref,"cf_ref",yelmo1%dyn%now%cf_ref)

            ! Additionally modify cf_ref 
            !call modify_cf_ref(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,yelmo1%grd,domain,f_cf,f_cf_lim)
            
            ! Make sure minimum is consistent with parameter settings
            where (yelmo1%dyn%now%cf_ref .lt. yelmo1%dyn%par%cb_min) yelmo1%dyn%now%cf_ref = yelmo1%dyn%par%cb_min

        else
            ! Define cf_ref inline 

            !call set_cf_ref(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,yelmo1%grd,domain)

            call set_cf_ref_new(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,yelmo1%grd,domain,yelmo1%par%grid_name,f_cf)

        end if 

    end if 

    ! ============================================================

    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base)
    call yelmo_init_state(yelmo1,time=time,thrm_method="robin-cold")

    ! Update boundary conditions again as necessary, in case of change after generating restart file 
    
    ! Update the ice_allowed mask based on domain definition 
    call ybound_define_ice_allowed(yelmo1%bnd,yelmo1%par%domain)
        
    ! Get initial boundary enhancement factor value
    call ice_enh_update(ice_enh,time=time)
    yelmo1%bnd%enh_srf = ice_enh%enh 
    
    ! Prescribe initial enhancement factor everywhere too 
    do k = 1, yelmo1%mat%par%nz_aa 
        yelmo1%mat%now%enh(:,:,k) = yelmo1%bnd%enh_srf
    end do 

    ! GHF 
    yelmo1%bnd%Q_geo = gthrm1%now%ghf 
    
    ! SMB and T_srf
    yelmo1%bnd%smb   = smbpal1%ann%smb*conv_we_ie*1e-3    ! [mm we/a] => [m ie/a]
    yelmo1%bnd%T_srf = smbpal1%ann%tsrf 

    ! Impose flux correction to smb 
    call modify_smb(yelmo1%bnd%smb,dsmb_now,dsmb_negis,yelmo1%bnd,yelmo1%grd,time)
    
    ! ============================================================
    ! Load or define cf_ref again, in case of restart file

    if (yelmo1%dyn%par%cb_method .eq. -1) then 
    
        if (load_cf_ref) then 

            ! Parse filename with grid information
            call yelmo_parse_path(file_cf_ref,yelmo1%par%domain,yelmo1%par%grid_name)

            ! Load cf_ref from specified file 
            call nc_read(file_cf_ref,"cf_ref",yelmo1%dyn%now%cf_ref)

            ! Additionally modify cf_ref 
            !call modify_cf_ref(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,yelmo1%grd,domain,f_cf,f_cf_lim)
            
            ! Make sure minimum is consistent with parameter settings
            where (yelmo1%dyn%now%cf_ref .lt. yelmo1%dyn%par%cb_min) yelmo1%dyn%now%cf_ref = yelmo1%dyn%par%cb_min
            
        else
            ! Define cf_ref inline 

            !call set_cf_ref(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,yelmo1%grd,domain)

            call set_cf_ref_new(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,yelmo1%grd,domain,yelmo1%par%grid_name,f_cf)
            
        end if 

    end if 

    ! ============================================================

    if (.not. yelmo1%par%use_restart) then 
        ! Run initialization steps 

        ! Run yelmo for several years with constant boundary conditions and topo
        ! to equilibrate thermodynamics and dynamics
        call yelmo_update_equil(yelmo1,time,time_tot=20e3,dt=dtt,topo_fixed=.FALSE.,dyn_solver="sia")
        call yelmo_update_equil(yelmo1,time,time_tot=10e3,dt=dtt,topo_fixed=.FALSE.)

        ! Write a restart file 
        call yelmo_restart_write(yelmo1,file_restart_init,time)
!         stop "**** Done ****"

    end if 

if (.not. check_init) then 
    call yelmo_update_equil(yelmo1,time,time_tot=1e3,dt=dtt,topo_fixed=.FALSE.)
end if 

    ! Reset dep_time to time_init everywhere to avoid complications related to equilibration 
    yelmo1%mat%now%dep_time = time_init 

!=== opt ===

    yelmo0  = yelmo1 
    smbpal0 = smbpal1 
    isos0   = isos1 

do iter = 1, n_iter

    ! Set time_init to production run value for last iteration
    if (iter .eq. n_iter) time_init = time_init_final  

    smbpal1 = smbpal0 
    isos1   = isos0 
    isos1%par%time = time_init 

    yelmo0%dyn%now%cf_ref = yelmo1%dyn%now%cf_ref
    yelmo1  = yelmo0 
    call yelmo_set_time(yelmo1,time_init)   ! For safety

    ! Reset dep_time to time_init everywhere to avoid complications related to equilibration 
    yelmo1%mat%now%dep_time = time_init 

    ! Now set current time to initial time 
    time    = time_init 

    write(iter_str,*) iter 
    iter_str = adjustl(iter_str)

    if ( iter .lt. n_iter .and. (.not. overwrite_files) ) then
        ! Optimization filenames 
        file1D       = trim(outfldr)//"yelmo1D-"//trim(iter_str)//".nc"
        file2D       = trim(outfldr)//"yelmo2D-"//trim(iter_str)//".nc"
        file_restart = trim(outfldr)//"yelmo_restart-"//trim(iter_str)//".nc"  
    else
        ! Production run filenames
        file1D       = trim(outfldr)//"yelmo1D.nc"
        file2D       = trim(outfldr)//"yelmo2D.nc"
        file_restart = trim(outfldr)//"yelmo_restart.nc" 
    end if 

!=== opt ===

    ! 2D file 
    call yelmo_write_init(yelmo1,file2D,time_init=time,units="years") 
    call write_step_2D_combined(yelmo1,isos1,snp1,mshlf1,smbpal1,dsmb_now,dtas_now,dpr_now,file2D,time=time)
    
    ! 1D file 
    call yelmo_write_reg_init(yelmo1,file1D,time_init=time,units="years",mask=yelmo1%bnd%ice_allowed)
    call yelmo_write_reg_step(yelmo1,file1D,time=time)  

if (check_init) stop 

    ! Advance timesteps
    do n = 1, ceiling((time_end-time_init)/dtt)

        ! Get current time 
        time = time_init + n*dtt

        ! == SEA LEVEL ==========================================================
        call sealevel_update(sealev,year_bp=time)
        yelmo1%bnd%z_sl  = sealev%z_sl 

        ! Update boundary enhancement factor 
        call ice_enh_update(ice_enh,time=time)
        yelmo1%bnd%enh_srf = ice_enh%enh 
        
if (calc_ice_sheet) then 

        ! Update cf_ref if desired
        if (yelmo1%dyn%par%cb_method .eq. -1 .and. (.not. load_cf_ref) ) then
!             call set_cf_ref(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,yelmo1%grd,domain)
            call set_cf_ref_new(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,yelmo1%grd,domain,yelmo1%par%grid_name,f_cf)
        end if 
        
        if (yelmo1%dyn%par%cb_method .eq. -1 .and. load_cf_ref .and. iter .lt. n_iter) then
            ! If using tuned cf, perform additional optimization here 

            !if (time .ge. -2e3 .and. time .lt. 0.0 .and. mod(time,500.0)==0 ) then
                ! Update cf_ref every 500yr for the last 2000 yrs, except for year 0.

                !err_scale = get_opt_param(time,time1=-2e3,time2=0.0,p1=2000.0,p2=500.0,m=1.0)
            
            if (time .eq. 0.0) then 
                ! Update cf_ref at present day 

                call update_cf_ref_errscaling(yelmo1%dyn%now%cf_ref,cf_ref_dot,yelmo1%tpo%now%H_ice, &
                                        yelmo1%bnd%z_bed,yelmo1%dyn%now%ux_s,yelmo1%dyn%now%uy_s,yelmo1%dta%pd%H_ice, &
                                            yelmo1%dta%pd%uxy_s,yelmo1%dta%pd%H_grnd.le.0.0_prec,yelmo1%grd%dx, &
                                            cf_min=yelmo1%dyn%par%cb_min,cf_max=1.0,sigma_err=1.0,sigma_vel=200.0, &
                                            err_scale=err_scale,fill_dist=fill_dist,optvar="ice")

            end if 


        end if 


        ! == Yelmo ice sheet ===================================================
        call yelmo_update(yelmo1,time)

        ! == ISOSTASY ==========================================================
        call isos_update(isos1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_sl,time)
        yelmo1%bnd%z_bed = isos1%now%z_bed

end if 

if (calc_transient_climate) then 
        ! == CLIMATE (ATMOSPHERE AND OCEAN) ====================================
        if (mod(time,5.0)==0) then
            ! Normal snapclim call 
            call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,time=time,domain=domain)

            ! Modify tas and precip
            call modify_tas(snp1%now%tas,snp1%now%ta_ann,snp1%now%ta_sum,dtas_now,dtas_hol,yelmo1%grd,time)
            call modify_pr(snp1%now%pr,snp1%now%pr_ann,dpr_now,dpr_hol,yelmo1%grd,time)
    
        end if 

        ! == SURFACE MASS BALANCE ==============================================

        call smbpal_update_monthly(smbpal1,snp1%now%tas,snp1%now%pr, &
                                   yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time) 
        yelmo1%bnd%smb   = smbpal1%ann%smb*conv_we_ie*1e-3       ! [mm we/a] => [m ie/a]
        yelmo1%bnd%T_srf = smbpal1%ann%tsrf 

        ! Impose flux correction to smb 
        call modify_smb(yelmo1%bnd%smb,dsmb_now,dsmb_negis,yelmo1%bnd,yelmo1%grd,time)

!         yelmo1%bnd%smb   = yelmo1%dta%pd%smb
!         yelmo1%bnd%T_srf = yelmo1%dta%pd%t2m
    
        ! == MARINE AND TOTAL BASAL MASS BALANCE ===============================
        call marshelf_calc_Tshlf(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                         yelmo1%bnd%z_sl,depth=snp1%now%depth,to_ann=snp1%now%to_ann, &
                         dto_ann=snp1%now%to_ann - snp1%clim0%to_ann)

        call marshelf_update(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                         yelmo1%bnd%z_sl,dx=yelmo1%grd%dx*1e-3)

end if 

        yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
        yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  

        ! == MODEL OUTPUT =======================================================

        write2D_now = .FALSE. 
        if (iter .lt. n_iter) then
            ! If inside of optimization iterations, then write output with constant frequency as
            ! defined by dt2D_out  
            if (mod(time,dt2D_out)==0) write2D_now = .TRUE. 
        else 
            ! For final transient simulation use variable 2D output writing 
            ! via function check_out2D_time (to reduce output when not needed)
            write2D_now = check_out2D_time(time)
        end if 

        if (write2D_now) then  
            call write_step_2D_combined(yelmo1,isos1,snp1,mshlf1,smbpal1,dsmb_now,dtas_now,dpr_now,file2D,time=time)
        end if 

        if (mod(time,dt1D_out)==0) then 
            call yelmo_write_reg_step(yelmo1,file1D,time=time) 
        end if 

        if (mod(time,dt_restart)==0) then 
            call yelmo_restart_write(yelmo1,file_restart,time=time) 
        end if 

        if (mod(time,100.0)==0) then
            write(*,"(a,i5,f14.4)") "yelmo:: iter, time = ", iter, time
        end if  

        if (abs(time).lt.1e-5) then 
            ! At year==0, write some statistics 
            write(*,*) "PDstats:", iter, yelmo1%dta%pd%rmse_H, yelmo1%dta%pd%rmse_uxy, yelmo1%dta%pd%rmse_loguxy 
        end if 
        
    end do 

    ! Write the restart file for the end of the simulation
    call yelmo_restart_write(yelmo1,file_restart,time=time) 

end do ! end iter 


    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)

    write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(time_end-time_init))/(cpu_dtime/3600.0), " kiloyears / hr"

contains

    subroutine write_step_2D_combined(ylmo,isos,snp,mshlf,srf,dsmb_now,dtas_now,dpr_now,filename,time)

        implicit none 
        
        type(yelmo_class),      intent(IN) :: ylmo
        type(isos_class),       intent(IN) :: isos 
        type(snapclim_class),   intent(IN) :: snp 
        type(marshelf_class),   intent(IN) :: mshlf 
        type(smbpal_class),     intent(IN) :: srf 
        real(prec),             intent(IN) :: dsmb_now(:,:) 
        real(prec),             intent(IN) :: dtas_now(:,:)  
        real(prec),             intent(IN) :: dpr_now(:,:)  
        !type(sediments_class),  intent(IN) :: sed 
        !type(geothermal_class), intent(IN) :: gthrm
        !type(isos_class),       intent(IN) :: isos
        
        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: time

        ! Local variables
        integer    :: ncid, n
        real(prec) :: time_prev 

        real(prec) :: uxy_rmse, H_rmse, zsrf_rmse, loguxy_rmse 
        real(prec), allocatable :: tmp(:,:) 
        real(prec), allocatable :: tmp1(:,:) 
        
        allocate(tmp(ylmo%grd%nx,ylmo%grd%ny))
        allocate(tmp1(ylmo%grd%nx,ylmo%grd%ny))

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
        
        ! Write temp and precip anomaly field
        call nc_write(filename,"dsmb_now",dsmb_now,units="m/a",long_name="smb scaling field", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid) 
        call nc_write(filename,"dtas_now",dtas_now,units="K",long_name="Temp. scaling field", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid) 
        call nc_write(filename,"dpr_now",dpr_now,units="m/a",long_name="Precip scaling field", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"mb_applied",ylmo%tpo%now%mb_applied,units="m",long_name="Applied net mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"H_grnd",ylmo%tpo%now%H_grnd,units="m",long_name="Ice thickness overburden", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="bar",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"calv",ylmo%tpo%now%calv,units="m/a ice equiv.",long_name="Calving rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"f_grnd",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"f_ice",ylmo%tpo%now%f_ice,units="1",long_name="Ice fraction in grid cell", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"dist_grline",ylmo%tpo%now%dist_grline,units="km", &
!                       long_name="Distance to nearest grounding-line point", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"cf_ref",ylmo%dyn%now%cf_ref,units="--",long_name="Bed friction scalar", &
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
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically integrated velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"Q_strn",ylmo%thrm%now%Q_strn/(rho_ice*ylmo%thrm%now%cp),units="K a-1",long_name="Strain heating", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="mW m-2",long_name="Basal frictional heating", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"bmb_grnd",ylmo%thrm%now%bmb_grnd,units="m/a ice equiv.",long_name="Basal mass balance (grounded)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_w",ylmo%thrm%now%H_w,units="m",long_name="Basal water layer thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ATT",ylmo%mat%now%ATT,units="a^-1 Pa^-3",long_name="Rate factor", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"dep_time",ylmo%mat%now%dep_time*1e-3,units="kyr",long_name="Deposition time", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"depth_iso",ylmo%mat%now%depth_iso,units="m",long_name="Isochronal layer depth", &
                      dim1="xc",dim2="yc",dim3="age_iso",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"enh_bar",ylmo%mat%now%enh_bar,units="1",long_name="Vertically averaged enhancement factor", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"visc_int",ylmo%mat%now%visc_int,units="Pa a m",long_name="Vertically integrated viscosity", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"enh",ylmo%mat%now%enh,units="",long_name="Enhancement factor", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"enh_srf",ylmo%bnd%enh_srf,units="",long_name="Enhancement factor boundary condition", &
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
        call nc_write(filename,"T_srf",ylmo%bnd%T_srf,units="K",long_name="Surface temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/a ice equiv.",long_name="Basal mass balance (shelf)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW m-2",long_name="Geothermal heat flux", &
                      dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
        
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
       
!         call nc_write(filename,"pr",snp%now%pr*1e-3*365.0_prec,units="m/a water equiv.",long_name="Precipitation (monthly mean)", &
!                       dim1="xc",dim2="yc",dim3="month",dim4="time",start=[1,1,1,n],ncid=ncid)
              
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
!         call nc_write(filename,"z_srf_pd_err",ylmo%dta%pd%err_z_srf,units="m",long_name="Surface elevation error wrt present day", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s_pd_err",ylmo%dta%pd%err_uxy_s,units="m/a",long_name="Surface velocity error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"depth_iso_pd_err",ylmo%dta%pd%err_depth_iso,units="m", &
                      long_name="Isochronal layer depth error wrt present day", &
                      dim1="xc",dim2="yc",dim3="age_iso",dim4="time",start=[1,1,1,n],ncid=ncid)

!         call nc_write(filename,"smb_errpd",ylmo%bnd%smb-ylmo%dta%pd%smb,units="m/a ice equiv.",long_name="Surface mass balance error wrt present day", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

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

        call nc_write(filename,"ux_s",ylmo%dyn%now%ux_s,units="m/a",long_name="Surface velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_s",ylmo%dyn%now%uy_s,units="m/a",long_name="Surface velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"uz",ylmo%dyn%now%uz,units="m/a",long_name="Vertical velocity", &
                      dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",start=[1,1,1,n],ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D_combined

    subroutine modify_tas(tas,tas_ann,tas_sum,dtas_now,dtas_hol,grd,time)

        implicit none 

        real(prec),         intent(INOUT) :: tas(:,:,:)     ! [K]
        real(prec),         intent(INOUT) :: tas_ann(:,:)   ! [K]
        real(prec),         intent(INOUT) :: tas_sum(:,:)   ! [K]
        real(prec),         intent(OUT)   :: dtas_now(:,:)  ! [K]
        real(prec),         intent(IN)    :: dtas_hol       ! [K] Holocene temp anomaly
        type(ygrid_class),  intent(IN)    :: grd
        real(prec),         intent(IN)    :: time 

        ! Local variables
        real(prec), allocatable :: dtas_0kyr(:,:)
        real(prec), allocatable :: dtas_12kyr(:,:) 
        
        integer    :: k, nx, ny
        real(prec) :: t0, t1, t2, t3, wt 

        allocate(dtas_0kyr(grd%nx,grd%ny))
        allocate(dtas_12kyr(grd%nx,grd%ny))
        
        dtas_12kyr = 0.0_prec
        dtas_0kyr  = 0.0_prec 
        
        t0 = -12e3
        t1 =  -8e3 
        t2 =  -6e3
        t3 =  -2e3

        dtas_now = 0.0_prec 

        if (time .lt. t0) then 

            dtas_now = 0.0_prec 

        else if (time .ge. t0 .and. time .lt. t1) then 
            
            wt      = (time - t0) / (t1-t0)
            dtas_now = (1.0-wt)*dtas_12kyr + wt*dtas_hol 
            
        else if (time .ge. t1 .and. time .lt. t2) then 

            dtas_now = dtas_hol
            
        else if (time .ge. t2 .and. time .lt. t3) then 

            wt      = (time - t2) / (t3-t2)
            dtas_now = (1.0-wt)*dtas_hol + wt*dtas_0kyr 
        
        else ! time==t3 

            dtas_now = dtas_0kyr 

        end if 

        do k = 1, size(tas,3)
            tas(:,:,k) = tas(:,:,k) + dtas_now
        end do 

        ! Update annual and summer mean too for consistency
        tas_ann = tas_ann + dtas_now 
        tas_sum = tas_sum + dtas_now 

        return 

    end subroutine modify_tas 

    subroutine modify_pr(pr,pr_ann,dpr_now,dpr_hol,grd,time)

        implicit none 

        real(prec),         intent(INOUT) :: pr(:,:,:)      ! [mm/d]
        real(prec),         intent(INOUT) :: pr_ann(:,:)    ! [mm/a]
        real(prec),         intent(OUT)   :: dpr_now(:,:)   ! [m/a]
        real(prec),         intent(IN)    :: dpr_hol        ! [m/a] Holocene precip anomaly
        type(ygrid_class),  intent(IN)    :: grd
        real(prec),         intent(IN)    :: time 

        ! Local variables
        real(prec), allocatable :: dpr_0kyr(:,:) 
        real(prec), allocatable :: dpr_12kyr(:,:) 
        
        integer    :: k, nx, ny, nslice
        real(prec) :: t0, t1, t2, t3, wt 

        allocate(dpr_0kyr(grd%nx,grd%ny))
        allocate(dpr_12kyr(grd%nx,grd%ny))
        
        dpr_12kyr = 0.0_prec
        dpr_0kyr  = 0.0_prec 
        
        t0 = -12e3
        t1 =  -8e3 
        t2 =  -6e3
        t3 =   0e3

        dpr_now = 0.0_prec 

        if (time .lt. t0) then 

            dpr_now = 0.0_prec 

        else if (time .ge. t0 .and. time .lt. t1) then 
            
            wt      = (time - t0) / (t1-t0)
            dpr_now = (1.0-wt)*dpr_12kyr + wt*dpr_hol 
            
        else if (time .ge. t1 .and. time .lt. t2) then 

            dpr_now = dpr_hol
            
        else if (time .ge. t2 .and. time .lt. t3) then 

            wt      = (time - t2) / (t3-t2)
            dpr_now = (1.0-wt)*dpr_hol + wt*dpr_0kyr 
        
        else ! time==t3 

            dpr_now = dpr_0kyr 

        end if 

        do k = 1, size(pr,3)
            pr(:,:,k) = pr(:,:,k) + dpr_now * (1e3/365.0_prec)      ! dpr [m/a] => [mm/d] 
        end do 

        ! Update annual mean too for consistency
        pr_ann = sum(pr,dim=3)/12.0_prec * 365.0_prec               ! [mm/d] => [mm/a]
        
        return 

    end subroutine modify_pr 

        subroutine modify_smb(smb,dsmb_now,dsmb_negis,bnd,grd,time)

        implicit none 

        real(prec),         intent(INOUT) :: smb(:,:) 
        real(prec),         intent(OUT)   :: dsmb_now(:,:) 
        real(prec),         intent(IN)    :: dsmb_negis
        type(ygrid_class),  intent(IN)    :: grd 
        type(ybound_class), intent(IN)    :: bnd 
        real(prec),         intent(IN)    :: time 

        ! Local variables
        real(prec), allocatable :: dsmb_0kyr(:,:) 
        real(prec), allocatable :: dsmb_hol(:,:) 
        real(prec), allocatable :: dsmb_12kyr(:,:) 
        
        real(prec) :: t0, t1, t2, t3, t4, t5, t6 
        real(prec) :: wt 

        allocate(dsmb_12kyr(grd%nx,grd%ny))
        allocate(dsmb_hol(grd%nx,grd%ny))
        allocate(dsmb_0kyr(grd%nx,grd%ny))
        
        dsmb_12kyr = 0.0_prec

        dsmb_hol = 0.0_prec 
        call scale_cf_gaussian(dsmb_hol,-1.0,x0= 500.0, y0=-1300.0,sigma=80.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dsmb_hol,-1.0,x0= 600.0, y0=-1500.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dsmb_hol,-1.0,x0= 600.0, y0=-1800.0,sigma=80.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dsmb_hol,-2.0,x0= 600.0, y0=-1900.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dsmb_hol,-2.0,x0= 600.0, y0=-2000.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dsmb_hol,-2.0,x0= 600.0, y0=-2100.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        
        call scale_cf_gaussian(dsmb_hol,-5.0,x0= 500.0, y0=-2300.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dsmb_hol,-5.0,x0= 330.0, y0=-2600.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dsmb_hol,-5.0,x0= 240.0, y0=-2700.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        
        ! NEGIS
        call scale_cf_gaussian(dsmb_hol,dsmb_negis,x0= 420.0, y0=-1150.0,sigma=150.0,xx=grd%x*1e-3,yy=grd%y*1e-3)


        dsmb_0kyr = 0.0_prec 
        call scale_cf_gaussian(dsmb_0kyr,-1.0,x0= 500.0, y0=-1300.0,sigma=80.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dsmb_0kyr,-1.0,x0= 600.0, y0=-1500.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dsmb_0kyr,-0.5,x0= 500.0, y0=-1700.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dsmb_0kyr,-1.0,x0= 600.0, y0=-1800.0,sigma=80.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dsmb_0kyr,-1.0,x0= 550.0, y0=-1900.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dsmb_0kyr,-1.0,x0= 600.0, y0=-2000.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dsmb_0kyr,-1.0,x0= 600.0, y0=-2100.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        
        call scale_cf_gaussian(dsmb_0kyr,-2.0,x0= 500.0, y0=-2300.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dsmb_0kyr,-2.0,x0= 330.0, y0=-2600.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dsmb_0kyr,-2.0,x0= 240.0, y0=-2700.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        
        call scale_cf_gaussian(dsmb_0kyr,-0.5,x0=-300.0, y0=-2650.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        
        ! NEGIS
        call scale_cf_gaussian(dsmb_0kyr,dsmb_negis,x0= 420.0, y0=-1150.0,sigma=150.0,xx=grd%x*1e-3,yy=grd%y*1e-3)

        ! Ensure more negative mass balance for ice-free points during Holocene
        ! (helps avoid too much growth out of NEGIS)
!         where(bnd%H_ice_ref .eq. 0.0_prec .and. bnd%z_bed_ref .lt. 0.0_prec) 
!             dsmb_0kyr = -2.0_prec
!             dsmb_hol  = -2.0_prec 
!         end where 
        
        t0 =  0e3 
        t1 =  0e3 
        t2 = -6e3 
        t3 = -8e3 
        t4 = -12e3 
        t5 = -18e3 
        t6 = -100e3 

        if (time .lt. t6) then 

            dsmb_now = 0.0_prec 

        else if (time .ge. t6 .and. time .lt. t5) then 

            ! For 'glacial times', reduce negative smb by eg 80%
            dsmb_now = 0.0_prec 
            where (smb .lt. 0.0_prec) dsmb_now = -smb !*0.80_prec 

        else if (time .ge. t5 .and. time .lt. t4) then

            dsmb_now = 0.0_prec 

            ! For 'deglacial times', reduce negative smb by eg 50%
            where (smb .lt. 0.0_prec) dsmb_now = -smb*0.5_prec 

        else if (time .ge. t4 .and. time .lt. t3) then

            wt       = (time - t4) / (t3-t4)
            dsmb_now = (1.0-wt)*dsmb_12kyr + wt*dsmb_hol 

        else if (time .ge. t3 .and. time .lt. t2) then 

            dsmb_now = dsmb_hol

        else if (time .ge. t2 .and. time .lt. t1) then 
   
            wt = (time - t2) / (t1-t2)
            dsmb_now = (1.0-wt)*dsmb_hol + wt*dsmb_0kyr 
        
        else 

            dsmb_now = dsmb_0kyr

        end if 


        ! Apply to smb field
        smb = smb + dsmb_now 

        return 

    end subroutine modify_smb 

        subroutine set_cf_ref_new(dyn,tpo,thrm,bnd,grd,domain,grid_name,f_cf)
        ! Set cf_ref [unitless] with location specific tuning 

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd  
        type(ygrid_class),  intent(IN)    :: grd
        character(len=*),   intent(IN)    :: domain 
        character(len=*),   intent(IN)    :: grid_name 
        real(prec),         intent(IN)    :: f_cf 

        ! Local variables 
        character(len=512) :: file_vel 
        real(prec), allocatable :: uxy_srf(:,:) 

        allocate(uxy_srf(grd%nx,grd%ny))

        call nml_read(path_par,"ctrl","file_vel",  file_vel)               ! Filename holding PD uxy_srf field
        call yelmo_parse_path(file_vel,domain,grid_name)

        call nc_read(file_vel,"uxy_srf",uxy_srf,missing_value=MV)

        ! Initial value everywhere 
        dyn%now%cf_ref = 0.4 
        where (uxy_srf .gt.  10.0) dyn%now%cf_ref = 0.4
        where (uxy_srf .gt.  20.0) dyn%now%cf_ref = 0.2
        where (uxy_srf .gt.  50.0) dyn%now%cf_ref = 0.08
        where (uxy_srf .gt. 100.0) dyn%now%cf_ref = 0.002
        where (uxy_srf .gt. 200.0) dyn%now%cf_ref = 0.001
        where (uxy_srf .eq. MV) dyn%now%cf_ref = dyn%par%cb_min 


        ! Additional tuning 
        call scale_cf_gaussian(dyn%now%cf_ref,0.2  ,x0=-350.0, y0=-1450.0,sigma= 50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dyn%now%cf_ref,0.2  ,x0=-250.0, y0=-1600.0,sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dyn%now%cf_ref,0.2  ,x0= -50.0, y0=-1900.0,sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dyn%now%cf_ref,0.1  ,x0=-250.0, y0=-2000.0,sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        
        call scale_cf_gaussian(dyn%now%cf_ref,0.2,  x0=-150.0, y0=-2850.0,sigma= 50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        
        call scale_cf_gaussian(dyn%now%cf_ref,0.1,  x0= -50.0, y0=-1000.0,sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dyn%now%cf_ref,0.05, x0=-200.0, y0=-1200.0,sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        
        call scale_cf_gaussian(dyn%now%cf_ref,0.005,x0= 450.0, y0=-1150.0,sigma= 50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dyn%now%cf_ref,0.005,x0= 330.0, y0=-1250.0,sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dyn%now%cf_ref,0.005,x0= 300.0, y0=-1400.0,sigma= 50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                     

        ! Finally multiply the whole thing by f_cf to scale field up or down 
        dyn%now%cf_ref = f_cf * dyn%now%cf_ref 

        ! Ensure minimum value for PD ice-free points 
        where (bnd%H_ice_ref .eq. 0.0 .and. bnd%z_bed_ref .lt. 0.0) dyn%now%cf_ref = dyn%par%cb_min
        
        ! Eliminate extreme values 
        where (dyn%now%cf_ref .lt. dyn%par%cb_min) dyn%now%cf_ref = dyn%par%cb_min

        return 

    end subroutine set_cf_ref_new

    subroutine modify_cf_ref(dyn,tpo,thrm,bnd,grd,domain,f_cf,f_cf_lim)
        ! Set cf_ref [unitless] with location specific tuning 

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd  
        type(ygrid_class),  intent(IN)    :: grd
        character(len=*),   intent(IN)    :: domain 
        real(prec),         intent(IN)    :: f_cf 
        real(prec),         intent(IN)    :: f_cf_lim 

        ! First reduce maximum value everywhere 
        where (dyn%now%cf_ref .gt. 0.3) dyn%now%cf_ref = 0.3 

        call scale_cf_gaussian(dyn%now%cf_ref,0.05, x0=-150.0, y0=-1700.0,sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dyn%now%cf_ref,0.05, x0= -50.0, y0=-1800.0,sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)

        call scale_cf_gaussian(dyn%now%cf_ref,0.1,  x0=   0.0, y0=-1200.0,sigma=200.0,xx=grd%x*1e-3,yy=grd%y*1e-3)

        call scale_cf_gaussian(dyn%now%cf_ref,0.01, x0=-150.0, y0= -950.0,sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dyn%now%cf_ref,0.01, x0= 300.0, y0=-1000.0,sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
    
        call scale_cf_gaussian(dyn%now%cf_ref,0.01, x0=-100.0, y0=-2200.0,sigma= 80.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
        call scale_cf_gaussian(dyn%now%cf_ref,0.05, x0=-240.0, y0=-2400.0,sigma= 50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
    

        ! Modify cf_ref
        if (trim(domain) .eq. "Greenland") then

!             call scale_cf_gaussian(dyn%now%cf_ref,0.30, x0=-100.0, y0=-1400.0,sigma= 50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.30, x0=   0.0, y0=-1500.0,sigma= 50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.30, x0=   0.0, y0=-1700.0,sigma= 70.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.30, x0= 400.0, y0=-1800.0,sigma= 80.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.30, x0= 100.0, y0=-1900.0,sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.30, x0= 400.0, y0=-2000.0,sigma= 80.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.30, x0= 150.0, y0=-2000.0,sigma= 80.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.30, x0= 100.0, y0=-2200.0,sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
            
!             call scale_cf_gaussian(dyn%now%cf_ref,0.30, x0=  80.0, y0=-2400.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.30, x0=  50.0, y0=-2550.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.30, x0=   0.0, y0=-2700.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.30, x0= -50.0, y0=-2800.0,sigma=50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
            
!             call scale_cf_gaussian(dyn%now%cf_ref,0.10, x0= 300.0, y0=-1000.0,sigma= 50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.10, x0= 300.0, y0=-2200.0,sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
            
!             call scale_cf_gaussian(dyn%now%cf_ref,0.05, x0=   0.0, y0=-1300.0,sigma=300.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.05, x0= 400.0, y0=-2200.0,sigma= 80.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.03, x0=-300.0, y0=-1650.0,sigma= 50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.03, x0= 400.0, y0=-1650.0,sigma= 50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.03, x0= 450.0, y0=-1650.0,sigma= 50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.03, x0= 450.0, y0=-1950.0,sigma= 80.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.02, x0= -80.0, y0=-1200.0,sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
            
!             call scale_cf_gaussian(dyn%now%cf_ref,0.01, x0= -80.0, y0=-1000.0,sigma= 80.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.01, x0=-250.0, y0=-1050.0,sigma= 80.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
            
!             call scale_cf_gaussian(dyn%now%cf_ref,0.005,x0= 450.0, y0=-1150.0,sigma= 50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.005,x0= 400.0, y0=-1250.0,sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!             call scale_cf_gaussian(dyn%now%cf_ref,0.005,x0= 450.0, y0=-2250.0,sigma= 50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
            
!             call scale_cf_gaussian(dyn%now%cf_ref,0.002,x0= 500.0, y0=-2300.0,sigma= 50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
            
        end if 

        ! Finally multiply the whole thing by f_cf to scale field up or down, for 
        ! cf_ref values above a threshold 
        where (dyn%now%cf_ref .gt. f_cf_lim) dyn%now%cf_ref = f_cf * dyn%now%cf_ref 

        ! Ensure minimum value for PD ice-free points 
        !where (bnd%H_ice_ref .eq. 0.0 .and. bnd%z_bed_ref .lt. 0.0) dyn%now%cf_ref = dyn%par%cb_min
        
        ! Eliminate extreme values 
        where (dyn%now%cf_ref .lt. dyn%par%cb_min) dyn%now%cf_ref = dyn%par%cb_min
        
        return 

    end subroutine modify_cf_ref

    subroutine set_cf_ref(dyn,tpo,thrm,bnd,grd,domain)
        ! Set cf_ref [unitless] with location specific tuning 

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd  
        type(ygrid_class),  intent(IN)    :: grd
        character(len=*),   intent(IN)    :: domain 
        
        ! Local variables
        integer :: i, j, nx, ny 
        integer :: i1, i2, j1, j2 
        real(prec) :: f_scale 
        real(prec), allocatable :: lambda1(:,:) 
        real(prec), allocatable :: lambda_bed(:,:)  

        nx = grd%nx 
        ny = grd%ny 

        allocate(lambda1(nx,ny))
        allocate(lambda_bed(nx,ny))
        
        lambda1 = 1.0_prec 

            ! Additionally modify cf_ref
            if (trim(domain) .eq. "Greenland") then

                ! Reduction
                call scale_cf_gaussian(lambda1,0.1,x0= 600.0, y0=-1800.0,sigma=200.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cf_gaussian(lambda1,0.1,x0= 600.0, y0=-1900.0,sigma=200.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cf_gaussian(lambda1,0.1,x0= 600.0, y0=-2000.0,sigma=200.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cf_gaussian(lambda1,0.1,x0= 600.0, y0=-2100.0,sigma=200.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                
            end if 

            ! =============================================================================
            ! Step 2: calculate lambda functions to scale cf_ref from default value 
            
            !------------------------------------------------------------------------------
            ! lambda_bed: scaling as a function of bedrock elevation

            select case(trim(dyn%par%cb_scale))

                case("lin_zb")
                    ! Linear scaling function with bedrock elevation
                    
                    lambda_bed = calc_lambda_bed_lin(bnd%z_bed,dyn%par%cb_z0,dyn%par%cb_z1)

                case("exp_zb")
                    ! Exponential scaling function with bedrock elevation
                    
                    lambda_bed = calc_lambda_bed_exp(bnd%z_bed,dyn%par%cb_z0,dyn%par%cb_z1)

                case("till_const")
                    ! Constant till friction angle

                    lambda_bed = calc_lambda_till_const(dyn%par%till_phi_const)

                case("till_zb")
                    ! Linear till friction angle versus elevation

                    lambda_bed = calc_lambda_till_linear(bnd%z_bed,bnd%z_sl,dyn%par%till_phi_min,dyn%par%till_phi_max, &
                                                            dyn%par%till_phi_zmin,dyn%par%till_phi_zmax)

                case DEFAULT
                    ! No scaling

                    lambda_bed = 1.0_prec

            end select 
            
            ! Ensure lambda_bed is not below lower limit [default range 0:1] 
            where (lambda_bed .lt. dyn%par%cb_min) lambda_bed = dyn%par%cb_min


            ! =============================================================================
            ! Step 3: calculate cf_ref [non-dimensional]
            
            dyn%now%cf_ref = dyn%par%cf_stream*lambda1*lambda_bed
            
        return 

    end subroutine set_cf_ref

    subroutine scale_cf_gaussian(cf_ref,cf_new,x0,y0,sigma,xx,yy)

        implicit none 

        real(prec), intent(INOUT) :: cf_ref(:,:)
        real(prec), intent(IN) :: cf_new
        real(prec), intent(IN) :: x0
        real(prec), intent(IN) :: y0
        real(prec), intent(IN) :: sigma
        real(prec), intent(IN) :: xx(:,:)
        real(prec), intent(IN) :: yy(:,:)

        ! Local variables 
        integer :: nx, ny 
        real(prec), allocatable :: wts(:,:)
        
        nx = size(cf_ref,1)
        ny = size(cf_ref,2)

        allocate(wts(nx,ny))

        ! Get Gaussian weights 
        wts = 1.0/(2.0*pi*sigma**2)*exp(-((xx-x0)**2+(yy-y0)**2)/(2.0*sigma**2))
        wts = wts / maxval(wts)

        ! Scale cf_ref
        cf_ref = cf_ref*(1.0-wts) + cf_new*wts

        return 

    end subroutine scale_cf_gaussian

    function check_out2D_time(time) result(write_now)

        implicit none 

        real(prec), intent(IN) :: time 
        logical :: write_now 

        if (time .le. -130e3) then 
            dt2D_out = 5e3 
        else if (time .gt. -130e3 .and. time .le. -118e3) then 
            dt2D_out = 500.0 
        else if (time .gt. -118e3 .and. time .le. -22e3) then
            dt2D_out = 5e3 
        else if (time .gt.  -22e3 .and. time .le. -12e3) then
            dt2D_out = 1e3 
        else ! time .gt. -12e3 
            dt2D_out = 500.0 
        end if 
        
        if (mod(time,dt2D_out)==0) then
            write_now = .TRUE. 
        else
            write_now = .FALSE. 
        end if 

        return 

    end function check_out2D_time


end program yelmox 



