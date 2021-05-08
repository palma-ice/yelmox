module transclim

    use ncio 
    use nml 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the working precision of the library (sp,dp)
    integer,  parameter :: wp = sp 

    ! Define default missing value 
    real(wp), parameter :: mv = -9999.0_wp 

    type transclim_param_class

        character(len=1024) :: filename
        character(len=56)   :: name 
        character(len=56)   :: units_in
        character(len=56)   :: units_out
        real(wp) :: unit_scale 
        real(wp) :: unit_offset
        logical  :: with_time 

        ! Internal parameters
        integer, allocatable :: dim(:)  
        integer  :: ndim 
        
    end type

    type transclim_class 

        type(transclim_param_class) :: par 

        real(wp) :: time_now
        
        real(wp), allocatable :: time(:)  
        real(wp), allocatable :: var(:,:,:) 
    end type 

    private 
    public :: transclim_class
    public :: transclim_update
    public :: transclim_init_nml 
    public :: transclim_init_arg

contains

    
    subroutine transclim_update(tclim,time)
        ! Routine to update transient climate forcing to match 
        ! current `time`. 

        implicit none 

        type(transclim_class),  intent(INOUT) :: tclim
        real(wp), optional,     intent(IN)    :: time       ! [yr] Current time 

        ! Local variables 
        integer :: k_now, k0, k1, nt  
        type(transclim_param_class) :: par 
        logical :: with_time 
        real(wp) :: time_now 

        ! Define shortcuts
        par = tclim%par 
        with_time = par%with_time 

        if (with_time) then 

            if (.not. present(time)) then 
                write(*,*) "transclim_update:: Error: current time must be given as an argument."
                stop 
            end if 

        end if 

        if (present(time)) then 
            time_now = time 
        else 
            time_now = tclim%time_now 
        end if 

        if (with_time .and. tclim%time_now .eq. time_now) then 

            ! Do nothing, the transclim object is already up to date 
            ! fo the current time. 

        else 

            ! 1. Determine the index of the current time, if needed
            if (with_time) then 
                nt = size(tclim%time)
                do k_now = 1, nt 
                    if (tclim%time(k_now) .eq. time_now) exit 
                end do
            end if 

            ! 2. Read variable and convert units as needed
            select case(par%ndim)

                case(1)

                    if (with_time) then 
                        ! 0D (point) variable plus time dimension 
                        call nc_read(par%filename,par%name,tclim%var(1,1,1),missing_value=mv, &
                                start=[k_now],count=[1])
                    else 
                        ! 1D variable
                        call nc_read(par%filename,par%name,tclim%var(:,1,1),missing_value=mv)
                    end if 

                case(2)

                    if (with_time) then 
                        ! 1D variable plus time dimension 
                        call nc_read(par%filename,par%name,tclim%var(:,1,1),missing_value=mv, &
                                start=[1,k_now],count=[par%dim(1),1])
                    else 
                        ! 2D variable 
                        call nc_read(par%filename,par%name,tclim%var(:,:,1),missing_value=mv)
                    end if 

                case(3)

                    if (with_time) then 
                        ! 2D variable plus time dimension 
                        call nc_read(par%filename,par%name,tclim%var(:,:,1),missing_value=mv, &
                                start=[1,1,k_now],count=[par%dim(1),par%dim(2),1])
                    else 
                        ! 3D variable 
                        call nc_read(par%filename,par%name,tclim%var,missing_value=mv)
                    end if 

                case(4)

                    if (with_time) then 
                        ! 3D variable plus time dimension 
                        call nc_read(par%filename,par%name,tclim%var,missing_value=mv, &
                                start=[1,1,1,k_now],count=[par%dim(1),par%dim(2),par%dim(3),1])
                    else 
                        ! 4D variable 
                        write(*,*) "transclim_update:: Error: 4D variable without time dimension &
                        &is not yet supported."
                        stop 
                    end if 

                case DEFAULT 

                    write(*,*) "transclim_update:: ndim not allowed."
                    write(*,*) "ndim = ", par%ndim 
                    stop 

            end select

            ! Apply scaling 
            where (tclim%var .ne. mv) 
                tclim%var = tclim%var*par%unit_scale + par%unit_offset
            end where 
            
        end if 

        return 

    end subroutine transclim_update

    subroutine transclim_init_nml(tclim,filename,group,domain,grid_name)
        ! Routine to load information related to a given 
        ! transient variable, so that it can be processed properly.

        implicit none 

        type(transclim_class), intent(INOUT) :: tclim
        character(len=*),      intent(IN)    :: filename
        character(len=*),      intent(IN)    :: group
        character(len=*),      intent(IN)    :: domain
        character(len=*),      intent(IN)    :: grid_name  
        
        ! Local variables 
        
        ! First load parameters from nml file 
        call transclim_par_load(tclim%par,filename,group,domain,grid_name)

        ! Perform remaining init operations 
        call transclim_init_data(tclim) 

        return 

    end subroutine transclim_init_nml

    subroutine transclim_init_arg(tclim,filename,name,units_in,units_out,scale,offset,with_time)
        ! Routine to load information related to a given 
        ! transient variable, so that it can be processed properly.

        implicit none 

        type(transclim_class), intent(INOUT) :: tclim
        character(len=*),      intent(IN)    :: filename
        character(len=*),      intent(IN)    :: name
        character(len=*),      intent(IN)    :: units_in 
        character(len=*),      intent(IN)    :: units_out 
        real(wp), optional,    intent(IN)    :: scale 
        real(wp), optional,    intent(IN)    :: offset 
        logical,  optional,    intent(IN)    :: with_time 

        ! Local variables 
        
        ! Define parameters from subroutine arguments 
        tclim%par%filename  = trim(filename) 
        tclim%par%name      = trim(name) 
        tclim%par%units_in  = trim(units_in) 
        tclim%par%units_out = trim(units_out) 
        
        tclim%par%unit_scale = 1.0_wp 
        if (present(scale)) tclim%par%unit_scale = scale 

        tclim%par%unit_offset = 0.0_wp 
        if (present(offset)) tclim%par%unit_offset = scale 
        
        tclim%par%with_time = .TRUE. 
        if (present(with_time)) tclim%par%with_time = with_time 

        ! Perform remaining init operations 
        call transclim_init_data(tclim) 
        
        return 

    end subroutine transclim_init_arg

    subroutine transclim_init_data(tclim)

        implicit none 

        type(transclim_class), intent(INOUT) :: tclim 

        ! Local variables  
        character(len=12), allocatable :: dim_names(:) 
        logical :: with_time 

        ! Local shortcut
        with_time = tclim%par%with_time 

        ! First make sure all data objects are deallocated 
        if (allocated(tclim%time)) deallocate(tclim%time)
        if (allocated(tclim%var))  deallocate(tclim%var)

        ! Get information from netcdf file 
        call nc_dims(tclim%par%filename,tclim%par%name,dim_names,tclim%par%dim)
        tclim%par%ndim = size(tclim%par%dim,1)

        ! Allocate data variables 
        select case(tclim%par%ndim)

            case(1)
                if (with_time) then 
                    allocate(tclim%time(tclim%par%dim(1)))
                    allocate(tclim%var(1,1,1))
                else 
                    allocate(tclim%var(tclim%par%dim(1),1,1))
                end if 
            case(2)
                if (with_time) then 
                    allocate(tclim%time(tclim%par%dim(2)))
                    allocate(tclim%var(tclim%par%dim(1),1,1))
                else 
                    allocate(tclim%var(tclim%par%dim(1),tclim%par%dim(2),1))
                end if 

            case(3)
                if (with_time) then
                    allocate(tclim%time(tclim%par%dim(3)))
                    allocate(tclim%var(tclim%par%dim(1),tclim%par%dim(2),1))
                else 
                    allocate(tclim%var(tclim%par%dim(1),tclim%par%dim(2),tclim%par%dim(3)))
                end if 
            case(4)
                if (with_time) then 
                    allocate(tclim%time(tclim%par%dim(4)))
                    allocate(tclim%var(tclim%par%dim(1),tclim%par%dim(2),tclim%par%dim(3)))
                else 
                    write(*,*) "transclim_init_data:: 4D array without time dimension is not yet supported."
                    stop 
                end if
            case DEFAULT 
                write(*,*) "transclim_init_data:: ndim not allowed."
                write(*,*) "ndim = ", tclim%par%ndim 
                stop 

        end select

        return 

    end subroutine transclim_init_data

    subroutine transclim_par_load(par,filename,group,domain,grid_name,init)

        type(transclim_param_class), intent(OUT) :: par 
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: group
        character(len=*), intent(IN) :: domain
        character(len=*), intent(IN) :: grid_name  
        logical, optional :: init 

        ! Local variables
        logical :: init_pars 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 

        call nml_read(filename,group,"filename",       par%filename,     init=init_pars)
        call nml_read(filename,group,"name",           par%name,         init=init_pars)
        call nml_read(filename,group,"units_in",       par%units_in,     init=init_pars)
        call nml_read(filename,group,"units_out",      par%units_out,    init=init_pars)
        call nml_read(filename,group,"unit_scale",     par%unit_scale,   init=init_pars)   
        call nml_read(filename,group,"unit_offset",    par%unit_offset,  init=init_pars)   
        call nml_read(filename,group,"with_time",      par%with_time,    init=init_pars)   
        
        ! Parse filename as needed
        call parse_path(par%filename,domain,grid_name)
        
        return

    end subroutine transclim_par_load

    subroutine parse_path(path,domain,grid_name)

        implicit none 

        character(len=*), intent(INOUT) :: path 
        character(len=*), intent(IN)    :: domain, grid_name 

        call nml_replace(path,"{domain}",   trim(domain))
        call nml_replace(path,"{grid_name}",trim(grid_name))
        
        return 

    end subroutine parse_path
    
    subroutine axis_init(x,x0,dx)

        implicit none 

        real(wp) :: x(:)
        real(wp), optional :: x0, dx
        real(wp) :: dx_tmp 
        integer :: i, nx  

        nx = size(x) 

        do i = 1, nx 
            x(i) = real(i-1,wp)
        end do 

        dx_tmp = 1.d0 
        if (present(dx)) dx_tmp = dx 
        
        x = x*dx_tmp  

        if (present(x0)) then 
            x = x + x0 
        else
            x = x + (-(nx-1.0)/2.0*dx_tmp)
        end if 

        return 
    end subroutine axis_init

end module transclim
