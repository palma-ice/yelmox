module dimarray

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

    type dimarray_param_class

        character(len=56) :: varnm 
        character(len=56) :: filename
        character(len=56) :: units
        integer, allocatable :: dim(:)  
        integer  :: ndim
        integer  :: dim_time 
        real(wp) :: unit_conv

    end type

    type dimarray_class 

        type(dimarray_param_class) :: par 

        real(wp), allocatable :: time(:) 
        real(wp), allocatable :: layer(:) 
        
        real(wp), allocatable :: time_now
        real(wp) :: var0D  
        real(wp), allocatable :: var1D(:) 
        real(wp), allocatable :: var2D(:,:) 
        real(wp), allocatable :: var3D(:,:,:) 
    end type 

    private 
    public :: dimarray_class
    public :: dimarray_update
    public :: dimarray_init

contains

    
    subroutine dimarray_update(da,time)
        ! Routine to update transient climate forcing to match 
        ! current `time`. 

        implicit none 

        type(dimarray_class),  intent(INOUT) :: da
        real(wp),               intent(IN)    :: time       ! [yr] Current time 

        ! Local variables 
        integer :: k_now, k0, k1, nt  
        type(dimarray_param_class) :: par 

        ! Define shortcut to par 
        par = da%par 

        if (da%time_now .eq. time) then 

            ! Do nothing, the dimarray object is already up to date 
            ! fo the current time. 

        else 

            ! 1. Determine the index of the current time 
            nt = size(da%time)
            do k_now = 1, nt 
                if (da%time(k_now) .eq. time) exit 
            end do

            ! 2. Read variable and convert units as needed
            select case(par%ndim)

                case(1)
                    ! 0D (point) variable plus time dimension 
                    call nc_read(par%filename,par%varnm,da%var0D,missing_value=mv, &
                            start=[k_now],count=[1])

                    !da%var0D
                case(2)
                    ! 1D variable plus time dimension 
                    call nc_read(par%filename,par%varnm,da%var1D,missing_value=mv, &
                            start=[1,k_now],count=[par%dim(1),1])

                case(3)
                    ! 2D variable plus time dimension 
                    call nc_read(par%filename,par%varnm,da%var2D,missing_value=mv, &
                            start=[1,1,k_now],count=[par%dim(1),par%dim(2),1])

                case(4)
                    ! 3D variable plus time dimension 
                    call nc_read(par%filename,par%varnm,da%var3D,missing_value=mv, &
                            start=[1,1,1,k_now],count=[par%dim(1),par%dim(2),par%dim(3),1])
                
                case DEFAULT 

                    write(*,*) "dimarray_update:: ndim not allowed."
                    write(*,*) "ndim = ", par%ndim 
                    stop 

            end select

            
        end if 

        return 

    end subroutine dimarray_update

    subroutine dimarray_init(da)
        ! Routine to load information related to a given 
        ! transient variable, so that it can be processed properly.

        implicit none 

        type(dimarray_class),  intent(INOUT) :: da
        
        ! Local variables 
        !


        
        return 

    end subroutine dimarray_init



end module dimarray
