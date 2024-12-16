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
        logical  :: trans_cal
        logical  :: trans_bp
        real(wp) :: time
        real(wp) :: time_elapsed
        real(wp) :: time_cal 
        real(wp) :: time_bp
        real(wp) :: time_pd
        real(wp) :: time_init
        integer  :: n
    end type

    private
    public :: tstep_class
    public :: tstep_init
    public :: tstep_update

contains

    subroutine tstep_init(ts,time_init,time_pd,method,units,trans_cal,trans_bp)

        implicit none

        type(tstep_class), intent(INOUT) :: ts
        real(wp),          intent(IN)    :: time_init
        real(wp),          intent(IN)    :: time_pd 
        character(len=*),  intent(IN)    :: method
        character(len=*),  intent(IN)    :: units
        logical,           intent(IN), optional :: trans_cal
        logical,           intent(IN), optional :: trans_bp

        ! Local variables
        real(wp) :: time_years

        ! Load parameters
        ts%method = trim(method)
        ts%units  = trim(units)

        ts%trans_cal = .TRUE.
        ts%trans_bp  = .TRUE.
        if (present(trans_cal)) ts%trans_cal = trans_cal
        if (present(trans_bp))  ts%trans_bp  = trans_bp 

        ! Set initial time
        ts%time_init    = convert_time_from_units(time_init,ts%units)
        ts%time_pd      = convert_time_from_units(time_pd,ts%units)
        
        ! Set output time based on method
        select case(trim(ts%method))
            case("cal")
                ts%time_cal = ts%time_init
                ts%time_bp  = ts%time_cal - ts%time_pd
                ts%time     = convert_time_to_units(ts%time_cal,ts%units)
            case("bp")
                ts%time_bp  = ts%time_init
                ts%time_cal = ts%time_bp + ts%time_pd
                ts%time     = convert_time_to_units(ts%time_bp,ts%units)
            case DEFAULT
                write(error_unit,*) "tstep_update:: Error: method not recognized."
                write(error_unit,*) "ts%method = ", trim(ts%method)
                stop
        end select

        ! Set elapsed time and interations too
        ts%time_elapsed = 0.0
        ts%n = 0

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
        ts%time_elapsed = round_time(ts%time_elapsed + dt)
        
        if (ts%trans_cal) ts%time_cal = round_time(ts%time_cal + dt)
        if (ts%trans_bp)  ts%time_bp  = round_time(ts%time_bp  + dt)

        ts%n = ts%n + 1 

        ! Set output time based on method
        select case(trim(ts%method))
            case("cal")
                ts%time = convert_time_to_units(ts%time_cal,ts%units)
            case("bp")
                ts%time = convert_time_to_units(ts%time_bp,ts%units)
        end select

        return

    end subroutine tstep_update

    function round_time(time) result (rtime)

        implicit none

        real(wp), intent(IN) :: time
        real(wp) :: rtime

        rtime = nint(time*1e3)*1e-3_wp

        return

    end function round_time

    function convert_time_to_units(time,units) result (time_units)

        implicit none

        real(wp),         intent(IN) :: time
        character(len=*), intent(IN) :: units
        real(wp) :: time_units

        select case(trim(units))
            case("years","yrs","year","yr")
                ! Pass, internal units match external units
            case("kiloyears","kyr","ka")
                time_units = time*1e-3
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
                ! Pass, internal units match external units
            case("kiloyears","kyr","ka")
                time_yr = time*1e3
        end select

        return

    end function convert_time_from_units
    
end module timestepping