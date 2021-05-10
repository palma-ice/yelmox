module varslice

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

    type varslice_param_class

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
        real(wp) :: time_par(3) 

    end type

    type varslice_class 

        type(varslice_param_class) :: par 

        real(wp) :: time_now
        
        real(wp), allocatable :: x(:) 
        real(wp), allocatable :: y(:)
        real(wp), allocatable :: lev(:)  

        real(wp), allocatable :: time(:)  
        real(wp), allocatable :: var(:,:,:) 
    end type 

    private 
    public :: varslice_class
    public :: varslice_update
    public :: varslice_init_nml 
    public :: varslice_init_arg
    public :: varslice_end 

contains

    
    subroutine varslice_update(vs,time)
        ! Routine to update transient climate forcing to match 
        ! current `time`. 

        implicit none 

        type(varslice_class),  intent(INOUT) :: vs
        real(wp), optional,     intent(IN)   :: time       ! [yr] Current time 

        ! Local variables 
        integer :: k_now, k0, k1, nt  
        type(varslice_param_class) :: par 
        logical :: with_time 
        real(wp) :: time_now 

        ! Define shortcuts
        par = vs%par 
        with_time = par%with_time 

        if (with_time) then 

            if (.not. present(time)) then 
                write(*,*) "varslice_update:: Error: current time must be given as an argument."
                stop 
            end if 

        end if 

        if (present(time)) then 
            time_now = time 
        else 
            time_now = vs%time_now 
        end if 

        if (with_time .and. vs%time_now .eq. time_now) then 

            ! Do nothing, the varslice object is already up to date 
            ! fo the current time. 

        else 

            ! 1. Determine the index of the current time, if needed
            if (with_time) then 
                nt = size(vs%time)
                do k_now = 1, nt 
                    if (vs%time(k_now) .eq. time_now) exit 
                end do
            end if 

            ! 2. Read variable and convert units as needed
            select case(par%ndim)

                case(1)

                    if (with_time) then 
                        ! 0D (point) variable plus time dimension 
                        call nc_read(par%filename,par%name,vs%var(1,1,1),missing_value=mv, &
                                start=[k_now],count=[1])
                    else 
                        ! 1D variable
                        call nc_read(par%filename,par%name,vs%var(:,1,1),missing_value=mv)
                    end if 

                case(2)

                    if (with_time) then 
                        ! 1D variable plus time dimension 
                        call nc_read(par%filename,par%name,vs%var(:,1,1),missing_value=mv, &
                                start=[1,k_now],count=[par%dim(1),1])
                    else 
                        ! 2D variable 
                        call nc_read(par%filename,par%name,vs%var(:,:,1),missing_value=mv)
                    end if 

                case(3)

                    if (with_time) then 
                        ! 2D variable plus time dimension 
                        call nc_read(par%filename,par%name,vs%var(:,:,1),missing_value=mv, &
                                start=[1,1,k_now],count=[par%dim(1),par%dim(2),1])
                    else 
                        ! 3D variable 
                        call nc_read(par%filename,par%name,vs%var,missing_value=mv)
                    end if 

                case(4)

                    if (with_time) then 
                        ! 3D variable plus time dimension 
                        call nc_read(par%filename,par%name,vs%var,missing_value=mv, &
                                start=[1,1,1,k_now],count=[par%dim(1),par%dim(2),par%dim(3),1])
                    else 
                        ! 4D variable 
                        write(*,*) "varslice_update:: Error: 4D variable without time dimension &
                        &is not yet supported."
                        stop 
                    end if 

                case DEFAULT 

                    write(*,*) "varslice_update:: ndim not allowed."
                    write(*,*) "ndim = ", par%ndim 
                    stop 

            end select

            ! Make sure crazy values have been set to missing (for safety)
            where (abs(vs%var) .ge. 1e10) vs%var = mv 

            ! Apply scaling 
            where (vs%var .ne. mv) 
                vs%var = vs%var*par%unit_scale + par%unit_offset
            end where 
            
        end if 

        return 

    end subroutine varslice_update

    subroutine varslice_init_nml(vs,filename,group,domain,grid_name,verbose)
        ! Routine to load information related to a given 
        ! transient variable, so that it can be processed properly.

        implicit none 

        type(varslice_class),   intent(INOUT) :: vs
        character(len=*),       intent(IN)    :: filename
        character(len=*),       intent(IN)    :: group
        character(len=*),       intent(IN), optional :: domain
        character(len=*),       intent(IN), optional :: grid_name  
        logical,                intent(IN), optional :: verbose 
        ! Local variables 
        
        ! First load parameters from nml file 
        call varslice_par_load(vs%par,filename,group,domain,grid_name,verbose)

        ! Perform remaining init operations 
        call varslice_init_data(vs) 

        return 

    end subroutine varslice_init_nml

    subroutine varslice_init_arg(vs,filename,name,units_in,units_out,scale,offset, &
                                                            with_time,time_par)
        ! Routine to load information related to a given 
        ! transient variable, so that it can be processed properly.

        implicit none 

        type(varslice_class), intent(INOUT) :: vs
        character(len=*),      intent(IN)   :: filename
        character(len=*),      intent(IN)   :: name
        character(len=*),      intent(IN)   :: units_in 
        character(len=*),      intent(IN)   :: units_out 
        real(wp), optional,    intent(IN)   :: scale 
        real(wp), optional,    intent(IN)   :: offset 
        logical,  optional,    intent(IN)   :: with_time 
        real(wp), optional,    intent(IN)   :: time_par(:) 

        ! Local variables 
        
        ! Define parameters from subroutine arguments 
        vs%par%filename  = trim(filename) 
        vs%par%name      = trim(name) 
        vs%par%units_in  = trim(units_in) 
        vs%par%units_out = trim(units_out) 
        
        vs%par%unit_scale = 1.0_wp 
        if (present(scale)) vs%par%unit_scale = scale 

        vs%par%unit_offset = 0.0_wp 
        if (present(offset)) vs%par%unit_offset = offset 
        
        vs%par%with_time = .TRUE. 
        if (present(with_time)) vs%par%with_time = with_time 

        vs%par%time_par = [0.0,0.0,0.0]
        if (present(time_par)) vs%par%time_par = time_par 

        ! Perform remaining init operations 
        call varslice_init_data(vs) 
        
        return 

    end subroutine varslice_init_arg

    subroutine varslice_init_data(vs)

        implicit none 

        type(varslice_class), intent(INOUT) :: vs 

        ! Local variables  
        character(len=12), allocatable :: dim_names(:) 
        logical :: with_time 

        ! Local shortcut
        with_time = vs%par%with_time 

        ! First make sure all data objects are deallocated 
        call varslice_end(vs)

        ! Get information from netcdf file 
        call nc_dims(vs%par%filename,vs%par%name,dim_names,vs%par%dim)
        vs%par%ndim = size(vs%par%dim,1)

        if (with_time) then

            ! Initialize time vector from user parameters 
            call axis_init(vs%time,x0=vs%par%time_par(1), &
                                   x1=vs%par%time_par(2), &
                                   dx=vs%par%time_par(3))

            ! Check to make sure time vector matches netcdf file length 
            if (size(vs%time,1) .ne. vs%par%dim(vs%par%ndim)) then 
                write(*,*) "varslice_init_data:: Error: generated time coordinate &
                &does not match the length of the time dimension in the netcdf file."
                write(*,*) "time_par:    ", vs%par%time_par 
                write(*,*) "size(time):  ", size(vs%time,1)
                write(*,*) "nt (netcdf): ", vs%par%dim(vs%par%ndim)
                write(*,*) "filename:    ", trim(vs%par%filename)
                stop 
            end if 

        end if 

        ! Allocate coordinate and data variables,
        ! load coordinates too 
        select case(vs%par%ndim)

            case(1)
                if (with_time) then 
                    allocate(vs%var(1,1,1))
                else 
                    allocate(vs%var(vs%par%dim(1),1,1))
                end if 
            case(2)
                if (with_time) then 
                    allocate(vs%x(vs%par%dim(1)))
                    allocate(vs%var(vs%par%dim(1),1,1))

                    if (nc_exists_var(vs%par%filename,dim_names(1))) then 
                        call nc_read(vs%par%filename,dim_names(1),vs%x)
                    else
                        call axis_init(vs%x,nx=vs%par%dim(1))
                    end if
                else 
                    allocate(vs%x(vs%par%dim(1)))
                    allocate(vs%y(vs%par%dim(2)))
                    allocate(vs%var(vs%par%dim(1),vs%par%dim(2),1))

                    if (nc_exists_var(vs%par%filename,dim_names(1))) then 
                        call nc_read(vs%par%filename,dim_names(1),vs%x)
                    else
                        call axis_init(vs%x,nx=vs%par%dim(1))
                    end if
                    if (nc_exists_var(vs%par%filename,dim_names(2))) then 
                        call nc_read(vs%par%filename,dim_names(2),vs%y)
                    else
                        call axis_init(vs%y,nx=vs%par%dim(2))
                    end if
                    
                end if 

            case(3)
                if (with_time) then
                    allocate(vs%x(vs%par%dim(1)))
                    allocate(vs%y(vs%par%dim(2)))
                    allocate(vs%var(vs%par%dim(1),vs%par%dim(2),1))

                    if (nc_exists_var(vs%par%filename,dim_names(1))) then 
                        call nc_read(vs%par%filename,dim_names(1),vs%x)
                    else
                        call axis_init(vs%x,nx=vs%par%dim(1))
                    end if
                    if (nc_exists_var(vs%par%filename,dim_names(2))) then 
                        call nc_read(vs%par%filename,dim_names(2),vs%y)
                    else
                        call axis_init(vs%y,nx=vs%par%dim(2))
                    end if
                    
                else 
                    allocate(vs%x(vs%par%dim(1)))
                    allocate(vs%y(vs%par%dim(2)))
                    allocate(vs%lev(vs%par%dim(3)))
                    allocate(vs%var(vs%par%dim(1),vs%par%dim(2),vs%par%dim(3)))

                    if (nc_exists_var(vs%par%filename,dim_names(1))) then 
                        call nc_read(vs%par%filename,dim_names(1),vs%x)
                    else
                        call axis_init(vs%x,nx=vs%par%dim(1))
                    end if
                    if (nc_exists_var(vs%par%filename,dim_names(2))) then 
                        call nc_read(vs%par%filename,dim_names(2),vs%y)
                    else
                        call axis_init(vs%y,nx=vs%par%dim(2))
                    end if
                    if (nc_exists_var(vs%par%filename,dim_names(3))) then 
                        call nc_read(vs%par%filename,dim_names(3),vs%lev)
                    else
                        call axis_init(vs%lev,nx=vs%par%dim(3))
                    end if
                    
                end if

                
                 
            case(4)
                if (with_time) then 
                    allocate(vs%x(vs%par%dim(1)))
                    allocate(vs%y(vs%par%dim(2)))
                    allocate(vs%lev(vs%par%dim(3)))
                    allocate(vs%var(vs%par%dim(1),vs%par%dim(2),vs%par%dim(3)))

                    if (nc_exists_var(vs%par%filename,dim_names(1))) then 
                        call nc_read(vs%par%filename,dim_names(1),vs%x)
                    else
                        call axis_init(vs%x,nx=vs%par%dim(1))
                    end if
                    if (nc_exists_var(vs%par%filename,dim_names(2))) then 
                        call nc_read(vs%par%filename,dim_names(2),vs%y)
                    else
                        call axis_init(vs%y,nx=vs%par%dim(2))
                    end if
                    if (nc_exists_var(vs%par%filename,dim_names(3))) then 
                        call nc_read(vs%par%filename,dim_names(3),vs%lev)
                    else
                        call axis_init(vs%lev,nx=vs%par%dim(3))
                    end if
                    
                else 
                    write(*,*) "varslice_init_data:: 4D array without time dimension is not yet supported."
                    stop 
                end if

                    
            case DEFAULT 
                write(*,*) "varslice_init_data:: ndim not allowed."
                write(*,*) "ndim = ", vs%par%ndim 
                stop 

        end select

        return 

    end subroutine varslice_init_data

    subroutine varslice_end(vs)
        ! Deallocate all variables

        implicit none 

        type(varslice_class), intent(INOUT) :: vs 

        if (allocated(vs%par%dim))  deallocate(vs%par%dim)
        if (allocated(vs%x))        deallocate(vs%x)
        if (allocated(vs%y))        deallocate(vs%y)
        if (allocated(vs%lev))      deallocate(vs%lev)
        if (allocated(vs%time))     deallocate(vs%time)
        if (allocated(vs%var))      deallocate(vs%var)
        
        return 

    end subroutine varslice_end

    subroutine varslice_par_load(par,filename,group,domain,grid_name,init,verbose)

        type(varslice_param_class), intent(OUT) :: par 
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: group
        character(len=*), intent(IN), optional :: domain
        character(len=*), intent(IN), optional :: grid_name  
        logical, optional :: init 
        logical, optional :: verbose 

        ! Local variables
        logical  :: init_pars 
        real(wp) :: time_par(3) 
        logical  :: print_summary 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 

        print_summary = .TRUE. 
        if (present(verbose)) print_summary = verbose 

        call nml_read(filename,group,"filename",       par%filename,     init=init_pars)
        call nml_read(filename,group,"name",           par%name,         init=init_pars)
        call nml_read(filename,group,"units_in",       par%units_in,     init=init_pars)
        call nml_read(filename,group,"units_out",      par%units_out,    init=init_pars)
        call nml_read(filename,group,"unit_scale",     par%unit_scale,   init=init_pars)   
        call nml_read(filename,group,"unit_offset",    par%unit_offset,  init=init_pars)   
        call nml_read(filename,group,"with_time",      par%with_time,    init=init_pars)   
        call nml_read(filename,group,"time_par",       par%time_par,     init=init_pars)   
        
        ! Parse filename as needed
        if (present(domain) .and. present(grid_name)) then
            call parse_path(par%filename,domain,grid_name)
        end if 

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
                write(*,*) "time_par    = ", par%time_par
            end if
        end if 

        return

    end subroutine varslice_par_load

    subroutine parse_path(path,domain,grid_name)

        implicit none 

        character(len=*), intent(INOUT) :: path 
        character(len=*), intent(IN)    :: domain, grid_name 

        call nml_replace(path,"{domain}",   trim(domain))
        call nml_replace(path,"{grid_name}",trim(grid_name))
        
        return 

    end subroutine parse_path
    
    subroutine axis_init(x,x0,x1,dx,nx)

        implicit none 

        real(wp), allocatable, intent(OUT) :: x(:)
        real(wp), optional :: x0
        real(wp), optional :: x1
        real(wp), optional :: dx
        integer,  optional :: nx 

        ! Local variables 
        integer :: i  
        real(wp) :: x0_now
        real(wp) :: x1_now
        real(wp) :: dx_now
        integer  :: nx_now 

        dx_now = 1.0_wp 
        if (present(dx)) dx_now = dx 

        x0_now = 0.0_wp 
        if (present(x0)) x0_now = x0 
        
        ! Note: if x1 is present, nx is ignored 
        if (present(x1)) then 
            x1_now = x1 
        else if (present(nx)) then 
            x1_now = (nx-1)*dx_now 
        else 
            write(*,*) "axis_init:: Error: either x1 or nx must be present."
            stop 
        end if 

        if (allocated(x)) deallocate(x)

        nx_now = (x1_now-x0_now)/dx_now + 1
        allocate(x(nx_now))

        do i = 1, nx_now 
            x(i) = x0_now + dx_now*real(i-1,wp)
        end do 

        return

    end subroutine axis_init

end module varslice
