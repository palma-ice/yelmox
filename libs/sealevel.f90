module sealevel

    use nml 
    use ncio 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: wp = sp 

    real(wp), parameter :: sec_year  = 365.0*24.0*60.0*60.0   ! [s/a]
    real(wp), parameter :: pi        = 3.14159265359

    type series_type
        character(len=512) :: filename 
        real(wp), allocatable :: time(:), var(:), sigma(:)
    end type 

    type sealevel_class 
        type(series_type) :: series 
        integer :: method 
        real(wp) :: z_sl_const 
        real(wp) :: time, z_sl, sigma  
    end type 

    private
    public :: sealevel_class 
    public :: sealevel_init 
    public :: sealevel_update 

contains 

    subroutine sealevel_init(sl,filename)

        implicit none 

        type(sealevel_class), intent(OUT) :: sl 
        character(len=*),     intent(IN)  :: filename 

        ! Local variables 
        character(len=56) :: varname 
        logical           :: use_nc 
        integer           :: n 

        call nml_read(filename,"sealevel","method",sl%method,init=.TRUE.)
        
        select case(sl%method)

            case(0)
                ! Load constant sea level value from parameter file 

                call nml_read(filename,"sealevel","z_sl_const",sl%z_sl_const,init=.TRUE.)
                

            case(1) 
                ! Load transient sea-level time series 

                ! Determine filename from which to load sea level time series 
                call nml_read(filename,"sealevel","sl_path",sl%series%filename,init=.TRUE.)
        
                use_nc = .FALSE. 
                n = len_trim(sl%series%filename)
                if (sl%series%filename(n-1:n) .eq. "nc") use_nc = .TRUE. 

                if (use_nc) then 

                    ! Get the variable name of interest
                    call nml_read(filename,"sealevel","sl_name",varname,init=.TRUE.)
                    
                    ! Read the time series from netcdf file 
                    call read_series_nc(sl%series,sl%series%filename,varname)

                else

                    ! Read the time series from ascii file
                    call read_series(sl%series,sl%series%filename)

                end if 
        
                ! Also set z_sl_const to zero for safety 
                sl%z_sl_const = 0.0_wp 

            case DEFAULT 

                write(*,*) ""
                write(*,*) "sealevel_init:: Error: method not recognized."
                write(*,*) "method = ", sl%method
                stop 

        end select

        return 

    end subroutine sealevel_init

    subroutine sealevel_update(sl,year_bp,z_sl_now)

        implicit none 

        type(sealevel_class), intent(INOUT) :: sl 
        real(wp),             intent(IN)    :: year_bp 
        real(wp), optional,   intent(IN)    :: z_sl_now 

        select case(sl%method)

            case(0)
                ! Assign sea-level constant 
                sl%time  = year_bp 
                sl%z_sl  = sl%z_sl_const
                sl%sigma = 0.0 

            case(1) 
                sl%time  = year_bp 
                sl%z_sl  = series_interp(sl%series,year_bp)
                sl%sigma = 0.0 
            
            case DEFAULT 

                write(*,*) ""
                write(*,*) "sealevel_update:: Error: method not recognized."
                write(*,*) "method = ", sl%method
                stop 

        end select

        if (present(z_sl_now)) then 
            ! Overwrite any options and use argument value 

            sl%z_sl = z_sl_now 

        end if 

        return 

    end subroutine sealevel_update

    subroutine read_series(series,filename)
        ! This subroutine will read a time series of
        ! two columns [time,var] from an ascii file.
        ! Header should be commented by "#" or "!"
        implicit none 

        type(series_type) :: series 
        character(len=*)  :: filename 

        integer, parameter :: f = 190
        integer, parameter :: nmax = 10000

        integer :: i, stat, n 
        character(len=256) :: str, str1 
        real(wp) :: x(nmax), y(nmax) 

        ! Open file for reading 
        open(f,file=filename,status="old")

        ! Read the header in the first line: 
        read(f,*,IOSTAT=stat) str

        do i = 1, nmax 
            read(f,'(a100)',IOSTAT=stat) str 

            ! Exit loop if the end-of-file is reached 
            if(IS_IOSTAT_END(stat)) exit 

            str1 = adjustl(trim(str))
!            str1=str
            if ( len(trim(str1)) .gt. 0 ) then 
                if ( .not. (str1(1:1) == "!" .or. &
                            str1(1:1) == "#") ) then 
                    read(str1,*) x(i), y(i) 
                end if
            end if  
        end do 


        ! Close the file
        close(f) 

        if (i .eq. nmax) then 
            write(*,*) "read_series:: warning: "// &
                       "Maximum length of time series reached, ", nmax
            write(*,*) "Time series in the file may be longer: ", trim(filename)
        end if 

        ! Allocate the time series object and store output data
        n = i-1 
        call series_allocate(series,n)

        series%time = x(1:n) 
        series%var  = y(1:n) 

        write(*,*) "read_series:: Time series read from file: "//trim(filename)
        write(*,*) "    range time: ",minval(series%time), maxval(series%time)
        write(*,*) "    range var : ",minval(series%var),  maxval(series%var)
        
        return 

    end subroutine read_series

    subroutine read_series_nc(series,filename,varname)
        ! This subroutine will read a time series of
        ! sea level from a netcdf file.

        implicit none 

        type(series_type) :: series 
        character(len=*)  :: filename 
        character(len=*)  :: varname 

        integer :: n 

        ! Allocate the time series object and store output data
        n = nc_size(filename,"time")

        call series_allocate(series,n)

        call nc_read(filename,"time",series%time)
        call nc_read(filename,varname,series%var)

        ! Convert from kiloyears to years
        ! ajr: this should not be hard coded!! 
        series%time = series%time*1e3 

        write(*,*) "read_series:: Time series read from file: "//trim(filename)
        write(*,*) "    range time: ",minval(series%time), maxval(series%time)
        write(*,*) "    range var : ",minval(series%var),  maxval(series%var)
        
        return 

    end subroutine read_series_nc

    function series_interp(series,year_bp) result(var)
        ! Wrapper for simple `interp_linear` function
        ! for series_types. 
        implicit none 

        type(series_type) :: series 
        real(wp) :: year_bp 
        real(wp) :: var 
        integer  :: nt, i 

        ! Interpolate series object
        var = interp_linear(series%time,series%var,xout=year_bp)

        return 

    end function series_interp

    function interp_linear(x,y,xout) result(yout)
        ! Simple linear interpolation of a point

        implicit none 
 
        real(wp), dimension(:), intent(IN) :: x, y
        real(wp), intent(IN) :: xout
        real(wp) :: yout 
        integer  :: i, j, n, nout 
        real(wp) :: alph

        n    = size(x) 

        if (xout .lt. x(1)) then
            yout = y(1)
        else if (xout .gt. x(n)) then
            yout = y(n)
        else
            do j = 1, n 
                if (x(j) .ge. xout) exit 
            end do

            if (j .eq. 1) then 
                yout = y(1) 
            else if (j .eq. n+1) then 
                yout = y(n)
            else 
                alph = (xout - x(j-1)) / (x(j) - x(j-1))
                yout = y(j-1) + alph*(y(j) - y(j-1))
            end if 
        end if 

        return 

    end function interp_linear
    
    subroutine series_allocate(series,nt)

        implicit none 

        type(series_type) :: series 
        integer :: nt 

        if (allocated(series%time))  deallocate(series%time)
        if (allocated(series%var))   deallocate(series%var)
        if (allocated(series%sigma)) deallocate(series%sigma)

        allocate(series%time(nt))
        allocate(series%var(nt))
        allocate(series%sigma(nt))
        
        ! Initialize variables to zero
        series%time  = 0.0
        series%var   = 0.0
        series%sigma = 0.0
        
        return 

    end subroutine series_allocate

end module sealevel
