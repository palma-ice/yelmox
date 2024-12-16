module timestepping

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    use nml
    use ncio 

    implicit none

    integer,  parameter :: sp  = kind(1.0)
    integer,  parameter :: dp  = kind(1.0d0)
    integer,  parameter :: wp  = sp 

    real(wp), parameter :: MV  = -9999.0_wp 

    real(wp), parameter :: time_tol = 1e-6 
    logical,  parameter :: verbose  = .TRUE. 

    type tstep_class
        character(len=56) :: units      ! External time units
        character(len=56) :: method     ! method of timekeeping
        logical  :: use_const_cal
        logical  :: use_const_bp
        real(wp) :: time
        real(wp) :: time_elapsed
        real(wp) :: time_cal 
        real(wp) :: time_bp
        real(wp) :: time_init
        real(wp) :: time_const_bp
        real(wp) :: time_const_cal
        integer  :: n

        ! Internal values
        real(wp) :: time_pd
        real(dp) :: comp_elapsed
        real(dp) :: comp_cal
        real(dp) :: comp_bp
        
    end type

    private
    public :: tstep_class
    public :: tstep_init
    public :: tstep_update
    public :: tstep_print

contains

    subroutine tstep_init(ts,time_init,method,units,const_cal,const_bp)
        ! method = "const": time_elapsed evolves, fixed time_cal and time_bp
        ! method = "cal"  : time_cal evolves, time_bp is set relative to it
        ! method = "bp"   : time_bp evolves, time_cal is set relative to it

        implicit none

        type(tstep_class), intent(INOUT) :: ts
        real(wp),          intent(IN)    :: time_init
        character(len=*),  intent(IN)    :: method
        character(len=*),  intent(IN)    :: units
        real(wp),          intent(IN), optional :: const_cal
        real(wp),          intent(IN), optional :: const_bp

        ! Local variables
        real(wp) :: time_years

        ! Load parameters
        ts%method = trim(method)
        ts%units  = trim(units)

        if (present(const_cal)) then
            ts%time_const_cal = convert_time_from_units(const_cal,ts%units)
            ts%use_const_cal  = .TRUE.
        else
            ts%time_const_cal = -1e8
            ts%use_const_cal  = .FALSE.
        end if
        
        if (present(const_bp)) then
            ts%time_const_bp = convert_time_from_units(const_bp,ts%units)
            ts%use_const_bp  = .TRUE.
        else
            ts%time_const_bp = -1e8
            ts%use_const_bp  = .FALSE.
        end if
        
        if (trim(ts%method) .eq. "const") then
            ! In this case, set both constant values to true independent of arguments provided
            ts%use_const_cal = .TRUE.
            ts%use_const_bp  = .TRUE.
        end if

        ! Set initial time
        ts%time_init    = convert_time_from_units(time_init,ts%units)

        ! Set time PD (calendar year) - used for converting between cal and bp timescales
        ts%time_pd = 1950.0 

        ! Set interations to zero
        ts%n = 0

        ! Set time keepers based on method
        select case(trim(ts%method))
            case("cal")
                
                ts%time_elapsed = 0.0
                ts%time_cal = ts%time_init
                ts%time_bp  = ts%time_cal - ts%time_pd
                ts%time     = convert_time_to_units(ts%time_cal,ts%units)

                if (ts%use_const_bp) ts%time_bp = ts%time_const_bp

                ! Consistency check
                if (ts%use_const_cal) then
                    write(error_unit,*) "tstep_init:: Error: a constant calendar time (const_cal) &
                    &has been specified with method='cal'. These cannot be used together."
                    write(error_unit,*) "const_cal = ", const_cal
                    stop
                end if

            case("bp")

                ts%time_elapsed = 0.0
                ts%time_bp  = ts%time_init
                ts%time_cal = ts%time_bp + ts%time_pd
                ts%time     = convert_time_to_units(ts%time_bp,ts%units)

                if (ts%use_const_cal) ts%time_cal = ts%time_const_cal
                
                ! Consistency check
                if (ts%use_const_bp) then
                    write(error_unit,*) "tstep_init:: Error: a constant before present time (const_bp) &
                    &has been specified with method='bp'. These cannot be used together."
                    write(error_unit,*) "const_bp = ", const_bp
                    stop
                end if
            
            case("const")

                ts%time_elapsed = ts%time_init      ! only case where time_elapsed can start at non-zero
                ts%time_cal = ts%time_const_cal
                ts%time_bp  = ts%time_const_bp
                ts%time     = convert_time_to_units(ts%time_elapsed,ts%units)

            case DEFAULT
                write(error_unit,*) "tstep_init:: Error: method not recognized."
                write(error_unit,*) "ts%method = ", trim(ts%method)
                stop
        end select

        ! Set summation compensation values to zero for each time keeper
        ts%comp_elapsed = 0.0
        ts%comp_cal     = 0.0
        ts%comp_bp      = 0.0

        return

    end subroutine tstep_init

    subroutine tstep_update(ts,dt)

        implicit none

        type(tstep_class), intent(INOUT) :: ts
        real(wp),          intent(IN)    :: dt 

        ! Local variables
        real(wp) :: dt_year

        ! Get timestep in years
        dt_year = convert_time_from_units(dt,ts%units)

        ! Update each time keeper and round for errors

        ts%n = ts%n + 1 

        call kahan_sum(ts%time_elapsed, ts%comp_elapsed, dt_year)

        if (ts%use_const_cal) then
            ts%time_cal = ts%time_const_cal
        else
            call kahan_sum(ts%time_cal, ts%comp_cal, dt_year)
        end if

        if (ts%use_const_bp) then
            ts%time_bp = ts%time_const_bp
        else
            call kahan_sum(ts%time_bp, ts%comp_bp, dt_year)
        end if
        
        ! Set output time based on method
        select case(trim(ts%method))
            case("cal")
                ts%time = convert_time_to_units(ts%time_cal,ts%units)
            case("bp")
                ts%time = convert_time_to_units(ts%time_bp,ts%units)
            case("const")
                ts%time = convert_time_to_units(ts%time_elapsed,ts%units)
        end select

        return

    end subroutine tstep_update

    subroutine tstep_print(ts)

        implicit none

        type(tstep_class), intent(IN) :: ts

        write(*,*) "ts: ", ts%n, ts%time, ts%time_cal, ts%time_bp, ts%time_elapsed

        return

    end subroutine tstep_print

    function round_time(time) result (rtime)

        implicit none

        real(wp), intent(IN) :: time
        real(wp) :: rtime

        rtime = nint(time*1e2)*1e-2_wp

        return

    end function round_time

    subroutine kahan_sum(time, comp, dt)
        ! Use Kahan Summation to avoid accumulation of summing
        ! errors over time.
        ! https://en.wikipedia.org/wiki/Kahan_summation_algorithm

        implicit none

        real(wp), intent(inout) :: time  ! Accumulated total
        real(dp), intent(inout) :: comp  ! Compensation for rounding errors
        real(wp), intent(in)    :: dt    ! Current increment
        
        ! Local variables
        real(8) :: y, t

        y = dble(dt) - comp    ! Compensate for lost low-order bits
        t = time + y           ! Temporarily add the compensated value
        comp = (t - time) - y  ! Compute new compensation
        time = t               ! Update the accumulated total

        return

    end subroutine kahan_sum

    function convert_time_to_units(time,units) result (time_units)

        implicit none

        real(wp),         intent(IN) :: time
        character(len=*), intent(IN) :: units
        real(wp) :: time_units

        select case(trim(units))
            case("years","yrs","year","yr")
                ! Internal units match external units
                time_units = time
            case("kiloyears","kyr","ka")
                time_units = time*1e-3
            case DEFAULT
                write(error_unit,*) "convert_time_to_units:: Error: units not recognized."
                write(error_unit,*) "units = ", trim(units)
        end select

        return

    end function convert_time_to_units
    
    function convert_time_from_units(time,units) result (time_yr)

        implicit none

        real(wp),         intent(IN) :: time
        character(len=*), intent(IN) :: units
        real(wp) :: time_yr

        select case(trim(units))
            case("years","yrs","year","yr")
                ! Internal units match external units
                time_yr = time
            case("kiloyears","kyr","ka")
                time_yr = time*1e3
            case DEFAULT
                write(error_unit,*) "convert_time_from_units:: Error: units not recognized."
                write(error_unit,*) "units = ", trim(units)
        end select

        return

    end function convert_time_from_units
    
end module timestepping