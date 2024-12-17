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
        character(len=56) :: units      ! Time units
        character(len=56) :: method     ! method of timekeeping
        real(wp) :: time_init
        real(wp) :: time_end
        logical  :: use_const_cal
        logical  :: use_const_rel
        real(wp) :: time
        real(wp) :: time_elapsed
        real(wp) :: time_cal 
        real(wp) :: time_rel
        real(wp) :: time_const_rel
        real(wp) :: time_const_cal
        integer  :: n

        ! Internal values
        real(wp) :: time_ref
        real(dp) :: comp_elapsed
        real(dp) :: comp_cal
        real(dp) :: comp_rel
        logical  :: is_transient
        logical  :: is_finished
    end type

    private
    public :: tstep_class
    public :: tstep_init
    public :: tstep_update
    public :: tstep_print
    public :: tstep_print_header

contains

    subroutine tstep_init(ts,time_init,time_end,method,units,time_ref,const_cal,const_rel)
        ! method = "const": time_elapsed evolves, fixed time_cal and time_rel
        ! method = "cal"  : time_cal evolves, time_rel is set relative to it
        ! method = "rel","sp","bp"  : time_rel evolves, time_cal is set relative to it

        implicit none

        type(tstep_class), intent(INOUT) :: ts
        real(wp),          intent(IN)    :: time_init
        real(wp),          intent(IN)    :: time_end
        character(len=*),  intent(IN)    :: method
        character(len=*),  intent(IN)    :: units
        real(wp),          intent(IN), optional :: time_ref
        real(wp),          intent(IN), optional :: const_cal
        real(wp),          intent(IN), optional :: const_rel

        ! Local variables
        real(wp) :: time_years

        ! Load parameters
        ts%method = trim(method)
        ts%units  = trim(units)
        
        if (present(time_ref)) then
            ts%time_ref = time_ref
        else
            ts%time_ref = 0.0
        end if

        if (present(const_cal)) then
            ts%time_const_cal = const_cal
            ts%use_const_cal  = .TRUE.
        else
            ts%time_const_cal = -1e8
            ts%use_const_cal  = .FALSE.
        end if
        
        if (present(const_rel)) then
            ts%time_const_rel = const_rel
            ts%use_const_rel  = .TRUE.
        else
            ts%time_const_rel = -1e8
            ts%use_const_rel  = .FALSE.
        end if
        
        if (trim(ts%method) .eq. "const") then
            ! In this case, set both constant values to true independent of arguments provided
            ts%use_const_cal  = .TRUE.
            ts%use_const_rel  = .TRUE.
        end if

        ! Set initial and end time
        ts%time_init = time_init
        ts%time_end  = time_end

        ! Initialize is_finished to .false.
        if (ts%time_end .gt. ts%time_init) then
            ts%is_finished = .FALSE.
        else
            write(error_unit,*) "tstep_init:: Error: time_init must be earlier than time_end."
            write(error_unit,*) "time_init = ", ts%time_init
            write(error_unit,*) "time_end  = ", ts%time_end
            stop
        end if

        ! Set interations to zero
        ts%n = 0

        ! Set time keepers based on method
        select case(trim(ts%method))
            case("cal")
                
                ts%time_elapsed = 0.0
                ts%time_cal     = ts%time_init
                ts%time_rel     = ts%time_cal - ts%time_ref
                ts%time         = ts%time_cal

                if (ts%use_const_rel) ts%time_rel = ts%time_const_rel

                ! Consistency check
                if (ts%use_const_cal) then
                    write(error_unit,*) "tstep_init:: warning: a constant calendar time (const_cal) &
                    &has been specified with method='cal'. The constant value will be ignored."
                    write(error_unit,*) "const_cal = ", const_cal
                    ts%use_const_cal = .FALSE.
                end if

            case("rel","sp","bp") ! Synonyms

                ts%time_elapsed = 0.0
                ts%time_rel     = ts%time_init
                ts%time_cal     = ts%time_rel + ts%time_ref
                ts%time         = ts%time_rel

                if (ts%use_const_cal) ts%time_cal = ts%time_const_cal
                
                ! Consistency check
                if (ts%use_const_rel) then
                    write(*,*) "tstep_init:: warning: a constant before present time (const_rel) &
                    &has been specified with one of method=['rel','sp','bp']. The constant value will be ignored."
                    write(*,*) "const_rel = ", const_rel
                    ts%use_const_rel = .FALSE.
                end if
            
            case("const")

                ts%time_elapsed = ts%time_init      ! only case where time_elapsed can start at non-zero
                ts%time_cal     = ts%time_const_cal
                ts%time_rel     = ts%time_const_rel
                ts%time         = ts%time_elapsed

            case DEFAULT
                write(error_unit,*) "tstep_init:: Error: method not recognized."
                write(error_unit,*) "ts%method = ", trim(ts%method)
                stop
        end select

        ! Set summation compensation values to zero for each time keeper
        ts%comp_elapsed = 0.0
        ts%comp_cal     = 0.0
        ts%comp_rel     = 0.0

        return
    
    end subroutine tstep_init

    subroutine tstep_update(ts,dt)

        implicit none

        type(tstep_class), intent(INOUT) :: ts
        real(wp),          intent(IN)    :: dt 

        if (ts%n .gt. 0) then

            ! Update each time keeper and round for errors

            call kahan_sum(ts%time_elapsed, ts%comp_elapsed, dt)

            if (ts%use_const_cal) then
                ts%time_cal = ts%time_const_cal
            else
                call kahan_sum(ts%time_cal, ts%comp_cal, dt)
            end if

            if (ts%use_const_rel) then
                ts%time_rel = ts%time_const_rel
            else
                call kahan_sum(ts%time_rel, ts%comp_rel, dt)
            end if
            
            ! Set output time based on method
            select case(trim(ts%method))
                case("cal")
                    ts%time = ts%time_cal
                case("rel","sp","bp")
                    ts%time = ts%time_rel
                case("const")
                    ts%time = ts%time_elapsed
            end select

        end if 

        ! Advance number of iterations
        ts%n = ts%n + 1 

        ! Finally check if time stepping is finished
        if (ts%time .ge. ts%time_end) then
            ts%is_finished = .TRUE.
        end if

        return

    end subroutine tstep_update

    subroutine tstep_print(ts)

        implicit none

        type(tstep_class), intent(IN) :: ts

        write(*,"(a5,i6,4f17.5)") "ts: ", ts%n, ts%time, ts%time_cal, ts%time_rel, ts%time_elapsed

        return

    end subroutine tstep_print

    subroutine tstep_print_header(ts)

        implicit none

        type(tstep_class), intent(IN) :: ts

        write(*,"(a5,a6,4a17)") "ts: ", "iter", "time", "time_cal", "time_rel", "time_elapsed"

        return

    end subroutine tstep_print_header

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
    
end module timestepping