module timer

    implicit none 

    type timer_class
        real(8) :: dtime_cpu
        real(8) :: dtime_mod
        real(8) :: time_cpu1, time_cpu2
        real(8) :: time_mod1, time_mod2

        character(len=128) :: str_dtime 
        character(len=128) :: str_rate  
    end type 

    interface timer_step
        module procedure    timer_step_none
        module procedure    timer_step_flt
        module procedure    timer_step_dble
    end interface

    private
    public :: timer_class 
    public :: timer_step 
    public :: timer_str
    public :: timer_str_comprate

contains

    subroutine timer_step_none(timer,step)

        implicit none 

        type(timer_class), intent(INOUT) :: timer 
        integer, intent(IN) :: step

        call timer_step_dble(timer,step,real(0.0,8))

        return

    end subroutine timer_step_none

    subroutine timer_step_flt(timer,step,time_mod)

        implicit none 

        type(timer_class), intent(INOUT) :: timer 
        integer, intent(IN) :: step 
        real(4), intent(IN) :: time_mod

        call timer_step_dble(timer,step,real(time_mod,8))

        return

    end subroutine timer_step_flt

    subroutine timer_step_dble(timer,step,time_mod)

        implicit none 

        type(timer_class), intent(INOUT) :: timer 
        integer, intent(IN) :: step 
        real(8), intent(IN) :: time_mod

        if (step .eq. 1) then 

            ! Get current cpu [s] and model [units_mod] times
            call timer_cpu_time(timer%time_cpu1)
            timer%time_mod1 = time_mod
        
        else if (step .eq. 2) then 
            
            ! Get current cpu [s] and model [units_mod] times
            call timer_cpu_time(timer%time_cpu2) 
            timer%time_mod2 = time_mod

            ! Get time differences [s], [units_mod]
            timer%dtime_cpu = timer%time_cpu2 - timer%time_cpu1 
            timer%dtime_mod = timer%time_mod2 - timer%time_mod1 

        else
            write(*,*) "timer_step:: unknown step: ", step 
            stop 
        end if 

        return 

    end subroutine timer_step_dble

    subroutine timer_str(timer,units,label)

        implicit none 

        type(timer_class), intent(INOUT) :: timer
        character(len=1),  intent(IN)    :: units 
        character(len=*),  intent(IN), optional :: label 

        ! Local variables
        character(len=512) :: time_label
        real(8) :: dtime 

        ! Get dtime in desired units
        dtime = convert_dtime(timer%dtime_cpu,units)

        time_label = "calc_time"
        if (present(label)) time_label = trim(label)
        
        write(timer%str_dtime,"(a,1x,f12.3,1x,a)") trim(time_label), dtime, "["//trim(units)//"]"

        return 

    end subroutine timer_str

    subroutine timer_str_comprate(timer,units,units_mod,label)

        implicit none

        type(timer_class), intent(INOUT) :: timer
        character(len=*),  intent(IN)    :: units 
        character(len=*),  intent(IN)    :: units_mod
        character(len=*),  intent(IN), optional :: label 

        ! Local variables
        character(len=512) :: time_label 
        character(len=512) :: units_rate
        real(8) :: dtime, mtime 
        real(8) :: rate 

        ! Get dtime in desired units
        dtime = convert_dtime(timer%dtime_cpu,units)

        ! Get model time in model units
        mtime = timer%dtime_mod 

        if (dtime .ne. 0.0) then
            
            rate = mtime / dtime

            units_rate = trim(units_mod)//"/"//trim(units)
            
            time_label = "calc_rate"
            if (present(label)) time_label = trim(label)

            write(timer%str_rate,"(a,1x,f12.3,1x,a)") trim(time_label), rate, "["//trim(units_rate)//"]"

        end if

        return

    end subroutine timer_str_comprate

    ! === INTERNAL FUNCTIONS ===

    subroutine timer_cpu_time(time)
        ! Calculate time intervals using system_clock.

        ! Note: for mulithreading, cpu_time() won't work properly.
        ! Instead, system_clock() should be used as it is here, 
        ! unless use_cpu_time=.TRUE. 

        !$ use omp_lib

        implicit none 

        real(8), intent(OUT) :: time 

        ! Local variables
        logical    :: using_omp

        ! Check openmp status - do not use global switch, since it may not have been initialized yet
        using_omp = .FALSE. 
        !$ using_omp = .TRUE.

        if (using_omp) then 
            ! --------------------------------------
            ! omp_get_wtime must be used for multithread openmp execution to get timing on master thread 
            ! The following lines will overwrite time with the result from omp_get_wtime on the master thread 

            !$ time = omp_get_wtime()

            ! --------------------------------------
            
        else 

            ! cpu_time can be used for serial execution to get timing on 1 processor
            call cpu_time(time)

        end if 

        return 

    end subroutine timer_cpu_time


    function convert_dtime(dtime,units) result(dtime_conv)
        ! Convert from [s] to [units]

        implicit none 
        
        real(8),          intent(IN) :: dtime 
        character(len=*), intent(IN) :: units 
        real(8) :: dtime_conv 

        select case(units)
            case("s")
                dtime_conv = dtime
            case("m")
                dtime_conv = dtime / 60.d0 
            case("h") 
                dtime_conv = dtime / (60.d0*60.d0)
            case DEFAULT
                write(*,*) "timer_print:: error: units not recognized: "//units
        end select

        return
        
    end function convert_dtime 

end module timer 
