

module sediments

    use nml 
    use ncio 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: prec = sp 

    type sediments_param_class
        logical            :: use_obs
        character(len=512) :: obs_path
        character(len=56)  :: obs_name 
        real(prec)         :: sed_limit
    end type 

    type sediments_state_class 
        real(prec), allocatable :: H(:,:)
    end type 

    type sediments_class
        type(sediments_param_class) :: par 
        type(sediments_state_class) :: now 
    end type

    private
    public :: sediments_class
    public :: sediments_init 
    public :: sediments_update
    public :: sediments_end 

contains 

    subroutine sediments_init(sed,filename,nx,ny,domain,grid_name)
        ! Eventually interface should use 'sed' object directly
        ! instead of input/output arrays 

        implicit none 

        type(sediments_class), intent(INOUT) :: sed  
        character(len=*),      intent(IN)    :: filename 
        integer,               intent(IN)    :: nx, ny 
        character(len=*),      intent(IN)    :: domain, grid_name


        ! Load sediments parameters
        call sediments_par_load(sed%par,filename,domain,grid_name,init=.TRUE.)

        ! Allocate the sediments object 
        call sediments_allocate(sed%now,nx,ny)

        ! ====================================
        !
        ! Read in the sediment thickness
        !
        ! ====================================
        
        ! Read in array 
        if (sed%par%use_obs) then 
            ! Load sediment data from file 
            call nc_read(sed%par%obs_path,sed%par%obs_name,sed%now%H)
            write(*,*) "sediments_init:: Sediments loaded from: "
            write(*,*) trim(sed%par%obs_path)//" : "//trim(sed%par%obs_name)
        
        else 
            ! Set sediment thickness to zero by default 
            sed%now%H = 0.0 

        end if 

        ! ====================================
        !
        ! Additional calculations
        !
        ! ====================================

!         ! Limit the sediments to values above 0.1 m (for stability?)       
!         sed%now%H = max(sed%now%H,0.1)        

        ! ====================================
        !
        ! Summary
        !
        ! ====================================

        write(*,*) "sediments_init:: range H: ", minval(sed%now%H), maxval(sed%now%H)

        return 

    end subroutine sediments_init


    subroutine sediments_update(sed,year_bp)

        implicit none 

        type(sediments_class) :: sed 

        real(prec) :: year_bp 

        ! === To do ====


        return 

    end subroutine sediments_update


    subroutine sediments_end(sed)

        implicit none 

        type(sediments_class) :: sed

        ! Deallocate sediments state object
        call sediments_deallocate(sed%now)

        return 

    end subroutine sediments_end


    subroutine sediments_par_load(par,filename,domain,grid_name,init)

        type(sediments_param_class), intent(OUT) :: par
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: domain, grid_name 
        logical, optional :: init 
        logical :: init_pars 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 

        call nml_read(filename,"sediments","use_obs",    par%use_obs,     init=init_pars)
        call nml_read(filename,"sediments","obs_path",   par%obs_path,    init=init_pars)
        call nml_read(filename,"sediments","obs_name",   par%obs_name,    init=init_pars)

        ! Replace gridding template values from path
        call parse_path(par%obs_path,domain,grid_name)
        
        return

    end subroutine sediments_par_load


    ! =======================================================
    !
    ! sediments memory management
    !
    ! =======================================================

    subroutine sediments_allocate(now,nx,ny)

        implicit none 

        type(sediments_state_class) :: now 
        integer :: nx, ny 

        ! Make object is deallocated
        call sediments_deallocate(now)

        ! Allocate sediments 
        allocate(now%H(nx,ny))
        
        ! Initialize to zero
        now%H = 0.0 
        
        return

    end subroutine sediments_allocate

    subroutine sediments_deallocate(now)

        implicit none 

        type(sediments_state_class) :: now 

        ! Allocate state objects
        if (allocated(now%H))    deallocate(now%H)

        return

    end subroutine sediments_deallocate

    subroutine parse_path(path,domain,grid_name)

        implicit none 

        character(len=*), intent(INOUT) :: path 
        character(len=*), intent(IN)    :: domain, grid_name 

        call nml_replace(path,"{domain}",   trim(domain))
        call nml_replace(path,"{grid_name}",trim(grid_name))
        
        return 

    end subroutine parse_path
    

end module sediments 


