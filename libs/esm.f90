module esm
    ! This module contains routines that help with performing the esm suite
    ! of experiments. 
    
    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    use nml  
    use ncio 
    use varslice

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the working precision of the library (sp,dp)
    integer,  parameter :: wp = sp 

    ! Define default missing value 
    real(wp), parameter :: mv = -9999.0_wp 

    ! Class for holding ice-forcing data from esm archives
    type esm_forcing_class
        
        ! Experiment information
        character(len=256)     :: gcm 
        character(len=256)     :: scenario
        character(len=256)     :: experiment 
        character(len=256)     :: domain 
        character(len=256)     :: grid_name 
        character(len=256)     :: ctrl_run_type

        ! Climatologies
        ! Atmosphere
        type(varslice_class)   :: ts_ref
        type(varslice_class)   :: pr_ref
        ! Ocean
        type(varslice_class)   :: to_ref
        type(varslice_class)   :: so_ref

        ! Atmospheric fields
        type(varslice_class)   :: ts
        type(varslice_class)   :: pr
        !type(varslice_class)   :: dts
        !type(varslice_class)   :: dpr

        ! Oceanic fields 
        type(varslice_class)   :: to
        type(varslice_class)   :: so
        !type(varslice_class)   :: dto
        !type(varslice_class)   :: dso
        type(varslice_class)   :: tf_cor

        ! General fields 
        type(varslice_class)   :: basins
        type(varslice_class)   :: zs_ref
        type(varslice_class)   :: zs_hist
        type(varslice_class)   :: zs_proj

        ! ESM fields
        ! Atmospheric fields
        type(varslice_class)   :: ts_ref 
        type(varslice_class)   :: pr_ref 

        type(varslice_class)   :: ts_hist 
        type(varslice_class)   :: pr_hist 

        type(varslice_class)   :: ts_proj
        type(varslice_class)   :: pr_proj

        ! Oceanic fields 
        type(varslice_class)   :: to_ref
        type(varslice_class)   :: so_ref

        type(varslice_class)   :: to_hist
        type(varslice_class)   :: so_hist

        type(varslice_class)   :: to_proj
        type(varslice_class)   :: so_proj
        
    end type

    type esm_ice_var_class
        character(len=56)  :: name 
        character(len=128) :: long_name
        character(len=12)  :: var_type
        character(len=128) :: standard_name 
        character(len=128) :: units_in
        character(len=128) :: units_out
        real(wp) :: unit_scale 
        real(wp) :: unit_offset
    end type

    type esm_experiment_class
        character(len=56)   :: expname
        character(len=56)   :: group
        character(len=56)   :: model
        character(len=256)  :: experiment
        character(len=256)  :: file_suffix
    end type
        
    ! Class for holding ice output for writing to standard formats...
    type esm_ice_class
        type(esm_ice_var_class), allocatable :: vars(:)
    end type 


    private
    public :: esm_forcing_class
    public :: esm_experiment_class
    public :: esm_ice_class

    ! General routines
    public :: esm_experiment_def
    public :: esm_forcing_init
    public :: esm_forcing_update
    
    public :: esm_write_init
    public :: esm_write_step
    
contains
    
    subroutine esm_experiment_def(ie,expname,filename,group,model,domain)

        implicit none

        type(esm_experiment_class), intent(OUT) :: ie
        character(len=*), intent(IN) :: expname
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: group 
        character(len=*), intent(IN) :: model
        character(len=*), intent(IN) :: domain 

        ! Save the experiment name (ctrlAE, expAE01, etc)
        ie%expname = trim(expname)
        ie%group   = trim(group)
        ie%model   = trim(model)
        ie%model   = trim(domain)

        ! Load parameters associated with this experiment
        call nml_read(filename,ie%expname,"experiment",   ie%experiment)
        
        ! Define the ouput filename according to protocol
        ie%file_suffix = trim(domain)//"_"//trim(group)//"_"//trim(model)//"_"//trim(ie%expname)//".nc"

        return
        
    end subroutine esm_experiment_def

    subroutine esm_forcing_init(esm,filename,domain,grid_name,gcm,scenario,experiment)

        implicit none 
    
        type(esm_forcing_class), intent(INOUT) :: esm
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: domain 
        character(len=*), intent(IN) :: grid_name 
        character(len=*), intent(IN), optional :: gcm
        character(len=*), intent(IN), optional :: scenario
        character(len=*), intent(IN), optional :: experiment
    
        ! Local variables 
        character(len=256) :: gcm_now
        character(len=256) :: scenario_now
            
        character(len=256) :: group_prefix 
    
        ! Surface elevation
        character(len=256) :: grp_zs_esm

        ! Atmosphere
        character(len=256) :: grp_ts_clm  
        character(len=256) :: grp_pr_clm
        character(len=256) :: grp_ts_ref  
        character(len=256) :: grp_pr_ref 
        character(len=256) :: grp_ts_hist 
        character(len=256) :: grp_pr_hist 
        character(len=256) :: grp_ts_proj 
        character(len=256) :: grp_pr_proj 
            
        ! Ocean
        character(len=256) :: grp_to_clm 
        character(len=256) :: grp_so_clm
        character(len=256) :: grp_to_ref 
        character(len=256) :: grp_so_ref 
        character(len=256) :: grp_to_hist 
        character(len=256) :: grp_so_hist 
        character(len=256) :: grp_to_proj 
        character(len=256) :: grp_so_proj 
        character(len=256) :: grp_tf_cor
            
        integer  :: iloc, k 
        real(wp) :: tmp
        real(wp) :: time_par_clm(3),time_par_hist(3),time_par_proj(3) 
    
        ! First determine whether gcm+scenario provided or experiment
        ! obtain valid values for gcm and scenario to start.
        if (present(experiment)) then 
            ! Experiment provided, obtain gcm+scenario from it 
    
            if (trim(experiment) .eq. "ctrl" .or. trim(experiment) .eq. "ctrl0") then 
                gcm_now      = trim(experiment)
                scenario_now = trim(experiment)
            else 
                iloc = index(experiment,"_")
                if (iloc == 0) then
                    write(error_unit,*) "esm_ant_forcing_init:: Error: experiment &
                    &argument must be defined as gcm_scenario."
                    write(error_unit,*) "experiment = ", trim(experiment)
                    stop
                end if
    
                gcm_now = experiment(1:iloc-1)
                scenario_now = experiment(iloc+1:len_trim(experiment))
    
            end if 
    
        else if (present(gcm) .and. present(scenario)) then 
            ! gcm and scenario provided go forward as usual 
    
            gcm_now      = trim(gcm)
            scenario_now = trim(scenario)
    
        else
            ! Arguments not provided 
            write(error_unit,*) "esm_forcing_init:: Error: gcm+scenario or experiment must be provided &
            &as arguments."
            stop
        end if
    
        ! Define the current experiment characteristics
        esm%gcm        = trim(gcm_now)
        esm%scenario   = trim(scenario_now) 
        esm%experiment = trim(esm%gcm)//"_"//trim(esm%scenario) 
    
        write(*,*)
        write(*,*) "esm_forcing_init:: summary"
        write(*,*) "ctrl_run_type: ", trim(esm%ctrl_run_type)
        write(*,*) "gcm:           ", trim(esm%gcm)
        write(*,*) "scenario:      ", trim(esm%scenario)
        write(*,*) "experiment:    ", trim(esm%experiment)
        write(*,*) 
    
        group_prefix = "gcm_"
    
        grp_ts_ref   = trim(group_prefix)//"ts_ref"
        grp_pr_ref   = trim(group_prefix)//"pr_ref"
        grp_ts_hist  = trim(group_prefix)//"ts_hist"
        grp_pr_hist  = trim(group_prefix)//"pr_hist"
        grp_ts_proj  = trim(group_prefix)//"ts_proj"
        grp_pr_proj  = trim(group_prefix)//"pr_proj"
                    
        grp_to_ref   = trim(group_prefix)//"to_ref"
        grp_so_ref   = trim(group_prefix)//"so_ref"
        grp_to_hist  = trim(group_prefix)//"to_hist"
        grp_so_hist  = trim(group_prefix)//"so_hist"
        grp_to_proj  = trim(group_prefix)//"to_proj"
        grp_so_proj  = trim(group_prefix)//"so_proj"
        grp_tf_cor   = "tf_cor"
 
        ! Set negative values to time_par so that values are used directly from the file
        ! [time_init,time_end,timestep]
        time_par_ref  = [-1.0_wp,-1.0_wp,-1.0_wp]
        time_par_hist = [-1.0_wp,-1.0_wp,-1.0_wp]
        time_par_proj = [-1.0_wp,-1.0_wp,-1.0_wp]
                    
        ! Initialize all variables from namelist entries 
    
        ! General fields 
        call varslice_init_nml_esm(esm%basins,  filename,"imbie_basins",domain,grid_name,esm%gcm,esm%scenario)
            
        ! Surface elevation (check if different fields)
        call varslice_init_nml_esm(esm%z_srf,   filename,trim(grp_zs_esm), domain,grid_name,esm%gcm,esm%scenario)

        ! Amospheric fields
        call varslice_init_nml_esm(esm%ts_ref,  filename,trim(grp_ts_ref), domain,grid_name,esm%gcm,esm%scenario,time_par_ref)
        call varslice_init_nml_esm(esm%pr_ref,  filename,trim(grp_pr_ref), domain,grid_name,esm%gcm,esm%scenario,time_par_ref)
            
        call varslice_init_nml_esm(esm%ts_hist, filename,trim(grp_ts_hist), domain,grid_name,esm%gcm,esm%scenario,time_par_hist)
        call varslice_init_nml_esm(esm%pr_hist, filename,trim(grp_pr_hist), domain,grid_name,esm%gcm,esm%scenario,time_par_hist)
    
        call varslice_init_nml_esm(esm%ts_proj, filename,trim(grp_ts_proj), domain,grid_name,esm%gcm,esm%scenario,time_par_proj)
        call varslice_init_nml_esm(esm%pr_proj, filename,trim(grp_pr_proj), domain,grid_name,esm%gcm,esm%scenario,time_par_proj)
    
        ! Oceanic fields
        call varslice_init_nml_esm(esm%to_ref,  filename,trim(grp_to_ref),domain,grid_name,esm%gcm,esm%scenario)
        call varslice_init_nml_esm(esm%so_ref,  filename,trim(grp_so_ref),domain,grid_name,esm%gcm,esm%scenario)
    
        call varslice_init_nml_esm(esm%to_hist, filename,trim(grp_to_hist),domain,grid_name,esm%gcm,esm%scenario,time_par_hist)
        call varslice_init_nml_esm(esm%so_hist, filename,trim(grp_so_hist),domain,grid_name,esm%gcm,esm%scenario,time_par_hist)
    
        call varslice_init_nml_esm(esm%to_proj, filename,trim(grp_to_proj),domain,grid_name,esm%gcm,esm%scenario,time_par_proj)
        call varslice_init_nml_esm(esm%so_proj, filename,trim(grp_so_proj),domain,grid_name,esm%gcm,esm%scenario,time_par_proj)
     
        call varslice_init_nml_esm(esm%tf_cor,  filename,trim(grp_tf_cor),domain,grid_name,esm%gcm,esm%scenario)
    
        ! Load time-independent fields 
    
        ! Surface elevation
        call varslice_update(esm%zs_esm)

        ! Amospheric fields 
    
        ! Oceanic fields
        call varslice_update(esm%to_ref)
        call varslice_update(esm%so_ref)
        call varslice_update(esm%tf_cor)
    
        ! Remove missing values from the ocean, if possible
        do k = 1, size(esm%to_ref%var,3)
            if (count(esm%to_ref%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = minval(esm%to_ref%var(:,:,k,1),mask=esm%to_ref%var(:,:,k,1) .ne. mv)
                where(esm%to_ref%var(:,:,k,1) .eq. mv) 
                    esm%to_ref%var(:,:,k,1) = tmp
                end where 
            end if
            if (count(esm%so_ref%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = maxval(esm%so_ref%var(:,:,k,1),mask=esm%so_ref%var(:,:,k,1) .ne. mv)
                where(esm%so_ref%var(:,:,k,1) .eq. mv) 
                    esm%so_ref%var(:,:,k,1) = tmp
                end where 
            end if
        end do
    
        return 
    
    end subroutine esm_forcing_init
    
    subroutine esm_forcing_update(esm,time,use_ref_atm,use_ref_ocn)
        ! Update climatic fields. These will be used as bnd conditions for Yelmo.
        ! Output are anomaly fields with respect to a reference field

        implicit none 

        type(esm_forcing_class), intent(INOUT) :: esm
        real(wp), intent(IN) :: time
        logical,  intent(IN), optional :: use_ref_atm 
        logical,  intent(IN), optional :: use_ref_ocn 

        ! Local variables 
        integer  :: k 
        real(wp) :: tmp 
        character(len=56) :: slice_method 

        ! Get slices for current time
        slice_method = "extrap" 
        
        select case(trim(esm%ctrl_run_type))
        
            case("ctrl","opt")
                ! If ctrl or opt, run only climatology?
                ! check
                esm%dts = 0.0_wp
                esm%dpr = 1.0_wp
                esm%dto = 0.0_wp
                esm%dso = 0.0_wp 
        
            case("transient")
                ! To do: what reference do i use for the historical/projection period?
                ! Historical period
                if (time .lt. time_hist_end) then
                    ! === Atmospheric fields === 
                    call varslice_update(esm%ts_hist,[time],method=slice_method)
                    call varslice_update(esm%pr_hist,[time],method=slice_method)
                    esm%ts = esm%ts_hist
                    esm%pr = esm%pr_hist
                    ! ===   Oceanic fields   ===
                    call varslice_update(esm%to_hist,[time],method=slice_method)
                    call varslice_update(esm%so_hist,[time],method=slice_method)
                    esm%to = esm%to_hist
                    esm%so = esm%so_hist
                else if (time .gt. time_hist_end .and. time .le. time_proj_end) then
                    ! === Atmospheric fields ===
                    call varslice_update(esm%ts_proj, [time],method=slice_method)
                    call varslice_update(esm%pr_proj, [time],method=slice_method)
                    esm%ts = esm%ts_proj
                    esm%pr = esm%pr_proj
                    ! ===   Oceanic fields   ===
                    call varslice_update(esm%to_proj,[time],method=slice_method)
                    call varslice_update(esm%so_proj,[time],method=slice_method)
                    esm%to = esm%to_proj
                    esm%so = esm%so_proj
                ! If longer than projection period, take the mean of the last 30 years
                else
                    ! === Atmospheric fields ===
                    call varslice_update(esm%ts_proj, [time_proj_end-30.0,time_proj_end],method="range_mean")
                    call varslice_update(esm%pr_proj, [time_proj_end-30.0,time_proj_end],method="range_mean")
                    esm%ts  = esm%ts_proj
                    esm%pr  = esm%pr_proj
                    ! ===   Oceanic fields   ===
                    call varslice_update(esm%to_proj, [time_proj_end-30.0,time_proj_end],method="range_mean")
                    call varslice_update(esm%so_proj, [time_proj_end-30.0,time_proj_end],method="range_mean")
                    esm%to  = esm%to_proj
                    esm%so  = esm%so_proj        
                end if

            case DEFAULT
                write(*,*) "esm_forcing_update:: Error: ctrl_run_type not recognized: "//trim(esm%ctrl_run_type)
                stop 
 
        end select
       
        ! Remove missing values from the ocean, if possible
        do k = 1, size(esm%to%var,3)
            if (count(esm%to%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = minval(esm%to%var(:,:,k,1),mask=esm%to%var(:,:,k,1) .ne. mv)
                where(esm%to%var(:,:,k,1) .eq. mv) 
                    esm%to%var(:,:,k,1) = tmp
                end where 
            end if
            if (count(esm%so%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = maxval(esm%so%var(:,:,k,1),mask=esm%so%var(:,:,k,1) .ne. mv)
                where(esm%so%var(:,:,k,1) .eq. mv) 
                    esm%so%var(:,:,k,1) = tmp
                end where 
            end if
        end do
        
        return 

    end subroutine esm_forcing_update

    subroutine esm_clim_update(esm,time,domain)
        ! Routine to compute the esm reference climatology
        ! This is a monthly file

        implicit none

        type(esm_forcing_class), intent(INOUT) :: esm
        real(wp),                intent(IN)    :: time
        character(len=*),        intent(IN)    :: domain

        ! Local variables 
        integer  :: year_rand, month, year 
        real(wp) :: rand, tmp 
        character(len=56) :: slice_method 

        ! Get slices for current time
        slice_method = "extrap" 

        ! Obtain reference climatologies
        if (esm%clim_var) then
            ! if climate variability is true, we select a random year from reference period
            call random_number(rand)
            !year_rand = NINT(year_clim_init + (year_clim_end-year_clim_init)*rand)
            year_rand = NINT((esm%time_eq(2)-esm%time_eq(1))*rand)
            ! === Atmospheric fields ===
            ! jablasco: this has to be monthly data
            do month = 1, 12 
                call varslice_update(esm%ts_ref[:,:,month],[year_rand+month],method=slice_method)
                call varslice_update(esm%pr_ref[:,:,month],[year_rand+month],method=slice_method)
            end do
            ! ===   Oceanic fields   ===
            ! jablasco: to do, ocean with variability
            !call varslice_update(esm%to_ref,[year_rand],method=slice_method)
            !call varslice_update(esm%so_ref,[year_rand],method=slice_method)
        else
            ! TO DO, mean of every month
            ! === Atmospheric fields ===
            call varslice_update(esm%ts_ref, [esm%time_eq(1),esm%time_eq(2)],method="range_mean")
            call varslice_update(esm%pr_ref, [esm%time_eq(1),esm%time_eq(2)],method="range_mean")
            !esm%ts_clim = esm%ts_ref
            !esm%pr_clim = esm%pr_ref
            ! ===   Oceanic fields   ===
            ! jablasco: first constant ocean
            !call varslice_update(esm%to_ref,[time_ref_1,time_ref_2],method="range_mean")
            !call varslice_update(esm%so_ref,[time_ref_1,time_ref_2],method="range_mean")
            !esm%to_ref
            !esm%so_ref        
        end if

        return

    end subroutine esm_clim_update

    ! === varslice wrapper routines with esm specific options ===================

    subroutine varslice_init_nml_esm(vs,filename,group,domain,grid_name,gcm,scenario,time_par)
        ! Routine to load information related to a given 
        ! transient variable, so that it can be processed properly.

        implicit none 

        type(varslice_class),   intent(INOUT) :: vs
        character(len=*),       intent(IN)    :: filename
        character(len=*),       intent(IN)    :: group
        character(len=*),       intent(IN)    :: domain
        character(len=*),       intent(IN)    :: grid_name
        character(len=*),       intent(IN)    :: gcm
        character(len=*),       intent(IN)    :: scenario
        real(wp), optional,     intent(IN)    :: time_par(3)

        ! First load parameters from nml file 
        call varslice_par_load_esm(vs%par,filename,group,domain,grid_name,gcm,scenario,verbose=.TRUE.)

        if (present(time_par)) then 
            if (minval(time_par) .ge. 0.0) then
                ! Use time_par option provided as an argument
                vs%par%time_par = time_par 
            end if
        end if

        ! Perform remaining init operations 
        call varslice_init_data(vs) 

        return 

    end subroutine varslice_init_nml_esm

    subroutine varslice_par_load_esm(par,filename,group,domain,grid_name,gcm,scenario,verbose)
        ! Wrapper to routine varslice::varslice_par_load() that includes
        ! additional parsing arguments of esm gcm and scenario. 

        type(varslice_param_class), intent(OUT) :: par 
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: group
        character(len=*), intent(IN) :: domain
        character(len=*), intent(IN) :: grid_name  
        character(len=*), intent(IN) :: gcm
        character(len=*), intent(IN) :: scenario
        logical :: verbose 

        ! Local variables
        logical  :: init_pars 
        real(wp) :: time_par(3) 
        logical  :: print_summary 

        init_pars     = .FALSE.
        print_summary = verbose 

        call nml_read(filename,group,"filename",       par%filename,     init=init_pars)
        call nml_read(filename,group,"name",           par%name,         init=init_pars)
        call nml_read(filename,group,"units_in",       par%units_in,     init=init_pars)
        call nml_read(filename,group,"units_out",      par%units_out,    init=init_pars)
        call nml_read(filename,group,"unit_scale",     par%unit_scale,   init=init_pars)   
        call nml_read(filename,group,"unit_offset",    par%unit_offset,  init=init_pars)   
        call nml_read(filename,group,"with_time",      par%with_time,    init=init_pars)
        call nml_read(filename,group,"monthly",        par%use_monthly,  init=init_pars)   
        call nml_read(filename,group,"time_par",       par%time_par,     init=init_pars)   
        
        ! Parse filename as needed
        call parse_path(par%filename,domain,grid_name)
        call parse_path_esm(par%filename,gcm,scenario)

        ! Make sure time parameters are consistent time_par=[x0,x1,dx]
        if (par%time_par(3) .eq. 0) par%time_par(2) = par%time_par(1) 

        ! Summary 
        if (print_summary) then  
            write(*,*) "Loading: ", trim(filename), ":: ", trim(group)
            write(*,*) "filename    = ", trim(par%filename)
            write(*,*) "name        = ", trim(par%name)
            write(*,*) "units_in    = ", trim(par%units_in)
            write(*,*) "units_out   = ", trim(par%units_out)
            write(*,*) "unit_scale  = ", par%unit_scale
            write(*,*) "unit_offset = ", par%unit_offset
            write(*,*) "with_time   = ", par%with_time
            if (par%with_time) then
                write(*,*) "monthly     = ", par%use_monthly
                write(*,*) "time_par    = ", par%time_par
            end if
        end if 

        return

    end subroutine varslice_par_load_esm
    
    subroutine parse_path_esm(path,gcm,scenario)

        implicit none

        character(len=*), intent(INOUT) :: path
        character(len=*), intent(IN)    :: gcm
        character(len=*), intent(IN)    :: scenario

        call nml_replace(path,"{gcm}",        trim(gcm))
        call nml_replace(path,"{scenario}",   trim(scenario))

        return

    end subroutine parse_path_esm
    
    subroutine parse_path(path,domain,grid_name)

        implicit none 

        character(len=*), intent(INOUT) :: path 
        character(len=*), intent(IN)    :: domain, grid_name 

        call nml_replace(path,"{domain}",   trim(domain))
        call nml_replace(path,"{grid_name}",trim(grid_name))
        
        return 

    end subroutine parse_path
    
    ! === ESM OUTPUT ROUTINES ==========

    subroutine esm_write_init(filename,xc,yc,time,lon,lat,area,map_name,lambda,phi)

        implicit none 

        character(len=*),   intent(IN) :: filename
        real(wp),           intent(IN) :: xc(:)
        real(wp),           intent(IN) :: yc(:)
        real(wp),           intent(IN) :: time
        real(wp),           intent(IN) :: lon(:,:)
        real(wp),           intent(IN) :: lat(:,:)
        real(wp),           intent(IN) :: area(:,:)
        character(len=*),   intent(IN) :: map_name
        real(wp),           intent(IN) :: lambda
        real(wp),           intent(IN) :: phi 

        ! Local variables 
        character(len=12) :: xnm 
        character(len=12) :: ynm 
        
        xnm = "xc"
        ynm = "yc" 

        ! === Initialize netcdf file and dimensions =========

        ! Create the netcdf file 
        call nc_create(filename)

        ! Add grid axis variables to netcdf file
        call nc_write_dim(filename,xnm,x=xc*1e-3,units="kilometers")
        call nc_write_attr(filename,xnm,"_CoordinateAxisType","GeoX")

        call nc_write_dim(filename,ynm,x=yc*1e-3,units="kilometers")
        call nc_write_attr(filename,ynm,"_CoordinateAxisType","GeoY")
        
        ! Add time axis with current value 
        call nc_write_dim(filename,"time", x=time,dx=1.0_wp,nx=1,units="years",unlimited=.TRUE.)
        
        ! Projection information 
        call nc_write_map(filename,map_name,dble(lambda),phi=dble(phi))

        ! Lat-lon information
        call nc_write(filename,"lon2D",lon,dim1=xnm,dim2=ynm,grid_mapping=map_name)
        call nc_write_attr(filename,"lon2D","_CoordinateAxisType","Lon")
        call nc_write(filename,"lat2D",lat,dim1=xnm,dim2=ynm,grid_mapping=map_name)
        call nc_write_attr(filename,"lat2D","_CoordinateAxisType","Lat")

        call nc_write(filename,"area",  area*1e-6,  dim1=xnm,dim2=ynm,grid_mapping=map_name,units="km^2")
        call nc_write_attr(filename,"area","coordinates","lat2D lon2D")
        

        return

    end subroutine esm_write_init

    subroutine esm_write_step(filename,file_nml,time)

        implicit none 

        character(len=*),   intent(IN) :: filename
        character(len=*),   intent(IN) :: file_nml 
        real(wp),           intent(IN) :: time 

        ! Local variables 
        integer    :: ncid, n
        real(wp) :: time_prev 
        type(esm_ice_class) :: esm 

        if (.FALSE.) then    
            ! Open the file for writing
            call nc_open(filename,ncid,writable=.TRUE.)

            ! Determine current writing time step 
            n = nc_size(filename,"time",ncid)
            call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
            if (abs(time-time_prev).gt.1e-5) n = n+1 

        end if 
        ! Load up output variable meta information 
        call esm_load_ice_var_info(esm,file_nml,verbose=.TRUE.)
 
        return

    end subroutine esm_write_step

    subroutine esm_load_ice_var_info(esmi,filename,verbose)

        implicit none 

        type(esm_ice_class), intent(OUT) :: esmi 
        character(len=*),    intent(IN)  :: filename 
        logical,             intent(IN)  :: verbose 

        ! Local variables
        integer :: n  
        integer, parameter :: n_variables = 38

        type(esm_ice_var_class) :: v 

        ! First initialize esm object to hold variable meta data 
        if (allocated(esmi%vars)) deallocate(esmi%vars)
        allocate(esmi%vars(n_variables))

        ! Load individual variables by namelist group 
        call ice_var_par_load(esmi%vars(1), filename,var_name="lithk")
        call ice_var_par_load(esmi%vars(2), filename,var_name="orog")
        call ice_var_par_load(esmi%vars(3), filename,var_name="base")
        call ice_var_par_load(esmi%vars(4), filename,var_name="topg")
        call ice_var_par_load(esmi%vars(5), filename,var_name="hfgeoubed")
        call ice_var_par_load(esmi%vars(6), filename,var_name="acabf")
        call ice_var_par_load(esmi%vars(7), filename,var_name="libmassbfgr")
        call ice_var_par_load(esmi%vars(8), filename,var_name="libmassbffl")
        call ice_var_par_load(esmi%vars(9), filename,var_name="dlithkdt")
        call ice_var_par_load(esmi%vars(10),filename,var_name="xvelsurf")
        call ice_var_par_load(esmi%vars(11),filename,var_name="yvelsurf")
        call ice_var_par_load(esmi%vars(12),filename,var_name="zvelsurf")
        call ice_var_par_load(esmi%vars(13),filename,var_name="xvelbase")
        call ice_var_par_load(esmi%vars(14),filename,var_name="yvelbase")
        call ice_var_par_load(esmi%vars(15),filename,var_name="zvelbase")
        call ice_var_par_load(esmi%vars(16),filename,var_name="xvelmean")
        call ice_var_par_load(esmi%vars(17),filename,var_name="yvelmean")
        call ice_var_par_load(esmi%vars(18),filename,var_name="litemptop")
        call ice_var_par_load(esmi%vars(19),filename,var_name="litempbotgr")
        call ice_var_par_load(esmi%vars(20),filename,var_name="litempbotfl")
        call ice_var_par_load(esmi%vars(21),filename,var_name="strbasemag")
        call ice_var_par_load(esmi%vars(22),filename,var_name="licalvf")
        call ice_var_par_load(esmi%vars(23),filename,var_name="lifmassbf")
        call ice_var_par_load(esmi%vars(24),filename,var_name="lifmassbf")
        call ice_var_par_load(esmi%vars(25),filename,var_name="ligroundf")
        call ice_var_par_load(esmi%vars(26),filename,var_name="sftgif")
        call ice_var_par_load(esmi%vars(27),filename,var_name="sftgrf")
        call ice_var_par_load(esmi%vars(28),filename,var_name="sftflf")

        call ice_var_par_load(esmi%vars(29),filename,var_name="lim")
        call ice_var_par_load(esmi%vars(30),filename,var_name="limnsw")
        call ice_var_par_load(esmi%vars(31),filename,var_name="iareagr")
        call ice_var_par_load(esmi%vars(32),filename,var_name="iareafl")
        call ice_var_par_load(esmi%vars(33),filename,var_name="tendacabf")
        call ice_var_par_load(esmi%vars(34),filename,var_name="tendlibmassbf")
        call ice_var_par_load(esmi%vars(35),filename,var_name="tendlibmassbffl")
        call ice_var_par_load(esmi%vars(36),filename,var_name="tendlicalvf")
        call ice_var_par_load(esmi%vars(37),filename,var_name="tendlifmassbf")
        call ice_var_par_load(esmi%vars(38),filename,var_name="tendligroundf")

        if (verbose) then 

            ! === Print summary =========

            write(*,"(a40,a8,a65,a15)") &
                                    "Variable name",    &
                                    "Type",             &
                                    "Standard name",    &
                                    "Unit"
            
            do n = 1, n_variables 
                v = esm%vars(n)
                write(*,"(a40,a8,a65,a15)") &
                                    trim(v%long_name),      &
                                    trim(v%var_type),       &
                                    trim(v%standard_name),  &
                                    trim(v%units_out) 
            end do 


        end if 

        return 

    end subroutine esm_load_ice_var_info

    subroutine ice_var_par_load(esmv,filename,var_name)
        ! Load parmaeters associated with a given ice variable

        implicit none 

        type(esm_ice_var_class),    intent(OUT) :: esmv 
        character(len=*),           intent(IN)  :: filename 
        character(len=*),           intent(IN)  :: var_name

        ! Local variables 
        character(len=56) :: group 

        group = "esm_out_"//trim(var_name)

        call nml_read(filename,group,"name",            esmv%name)
        call nml_read(filename,group,"long_name",       esmv%long_name)
        call nml_read(filename,group,"var_type",        esmv%var_type)
        call nml_read(filename,group,"standard_name",   esmv%standard_name)
        call nml_read(filename,group,"units_in",        esmv%units_in)
        call nml_read(filename,group,"units_out",       esmv%units_out)
        call nml_read(filename,group,"unit_scale",      esmv%unit_scale)
        call nml_read(filename,group,"unit_offset",     esmv%unit_offset)

        return 

    end subroutine ice_var_par_load

end module esm
