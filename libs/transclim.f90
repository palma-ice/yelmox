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

        character(len=56) :: name 
        character(len=56) :: filename
        character(len=56) :: units_in
        character(len=56) :: units_out
        real(wp) :: unit_scale 
        real(wp) :: unit_offset
        
        ! Internal parameters
        integer, allocatable :: dim(:)  
        integer  :: ndim 
        
    end type

    type transclim_class 

        type(transclim_param_class) :: par 

        real(wp), allocatable :: time_now
        
        real(wp), allocatable :: time(:) 
        
        real(wp) :: var0D  
        real(wp), allocatable :: var1D(:) 
        real(wp), allocatable :: var2D(:,:) 
        real(wp), allocatable :: var3D(:,:,:) 
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
        real(wp),               intent(IN)    :: time       ! [yr] Current time 

        ! Local variables 
        integer :: k_now, k0, k1, nt  
        type(transclim_param_class) :: par 

        ! Define shortcut to par 
        par = tclim%par 

        if (tclim%time_now .eq. time) then 

            ! Do nothing, the transclim object is already up to date 
            ! fo the current time. 

        else 

            ! 1. Determine the index of the current time 
            nt = size(tclim%time)
            do k_now = 1, nt 
                if (tclim%time(k_now) .eq. time) exit 
            end do

            ! 2. Read variable and convert units as needed
            select case(par%ndim)

                case(1)
                    ! 0D (point) variable plus time dimension 
                    call nc_read(par%filename,par%name,tclim%var0D,missing_value=mv, &
                            start=[k_now],count=[1])

                    if (tclim%var0D .ne. mv) then 
                        tclim%var0D = tclim%var0D*par%unit_scale + par%unit_offset
                    end if 

                case(2)
                    ! 1D variable plus time dimension 
                    call nc_read(par%filename,par%name,tclim%var1D,missing_value=mv, &
                            start=[1,k_now],count=[par%dim(1),1])

                    where (tclim%var1D .ne. mv) 
                        tclim%var1D = tclim%var1D*par%unit_scale + par%unit_offset
                    end where 

                case(3)
                    ! 2D variable plus time dimension 
                    call nc_read(par%filename,par%name,tclim%var2D,missing_value=mv, &
                            start=[1,1,k_now],count=[par%dim(1),par%dim(2),1])

                    where (tclim%var2D .ne. mv) 
                        tclim%var2D = tclim%var2D*par%unit_scale + par%unit_offset
                    end where 
                    
                case(4)
                    ! 3D variable plus time dimension 
                    call nc_read(par%filename,par%name,tclim%var3D,missing_value=mv, &
                            start=[1,1,1,k_now],count=[par%dim(1),par%dim(2),par%dim(3),1])
                    
                    where (tclim%var3D .ne. mv) 
                        tclim%var3D = tclim%var3D*par%unit_scale + par%unit_offset
                    end where 
                    
                case DEFAULT 

                    write(*,*) "transclim_update:: ndim not allowed."
                    write(*,*) "ndim = ", par%ndim 
                    stop 

            end select

            
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

    subroutine transclim_init_arg(tclim,filename,name,units_in,units_out,scale,offset)
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
        
        ! Local variables 
        
        ! Define parameters from subroutine arguments 
        tclim%par%filename  = trim(filename) 
        tclim%par%name      = trim(name) 
        tclim%par%units_in  = trim(units_in) 
        tclim%par%units_out = trim(units_out) 
        
        tclim%par%unit_scale = 1.0_wp 
        if (present(scale)) tclim%par%unit_scale = scale 

        tclim%par%unit_offset = 1.0_wp 
        if (present(offset)) tclim%par%unit_offset = scale 
        

        ! Perform remaining init operations 
        call transclim_init_data(tclim) 
        
        return 

    end subroutine transclim_init_arg

    subroutine transclim_init_data(tclim)

        implicit none 

        type(transclim_class), intent(INOUT) :: tclim 

        ! Local variables  
        character(len=12), allocatable :: dim_names(:) 

        ! First make sure all data objects are deallocated 
        if (allocated(tclim%time)) deallocate(tclim%time)
        if (allocated(tclim%var1D)) deallocate(tclim%var1D)
        if (allocated(tclim%var2D)) deallocate(tclim%var2D)
        if (allocated(tclim%var3D)) deallocate(tclim%var3D)

        ! Get information from netcdf file 
        call nc_dims(tclim%par%filename,tclim%par%name,dim_names,tclim%par%dim)
        tclim%par%ndim = size(tclim%par%dim,1)

        ! Allocate data variables 
        select case(tclim%par%ndim)

            case(1)
                allocate(tclim%time(tclim%par%dim(1)))
            case(2)
                allocate(tclim%time(tclim%par%dim(2)))
                allocate(tclim%var1D(tclim%par%dim(1)))
            case(3)
                allocate(tclim%time(tclim%par%dim(3)))
                allocate(tclim%var2D(tclim%par%dim(1),tclim%par%dim(2)))
            case(4)
                allocate(tclim%time(tclim%par%dim(4)))
                allocate(tclim%var3D(tclim%par%dim(1),tclim%par%dim(2),tclim%par%dim(3)))
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
    
end module transclim
