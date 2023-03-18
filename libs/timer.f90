module timer

    implicit none 

    type timer_class
        
        real(8) :: time_cpu1
        real(8) :: time_cpu2
        real(8) :: time_mod1
        real(8) :: time_mod2
        
        character(len=56) :: label(100)
        real(8) :: dtime_cpu(100)
        real(8) :: dtime_mod(100)
        real(8) :: dtime_cpu_tot(100)
        real(8) :: dtime_mod_tot(100)
        integer :: nstep_tot(100)
        
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

    subroutine timer_step_none(timer,comp,label,reset)

        implicit none 

        type(timer_class), intent(INOUT) :: timer 
        integer,           intent(IN)    :: comp
        character(len=*),  intent(IN), optional :: label
        logical,           intent(IN), optional :: reset
        
        call timer_step_dble(timer,comp,real(0.0,8),label,reset)

        return

    end subroutine timer_step_none

    subroutine timer_step_flt(timer,comp,time_mod,label,reset)

        implicit none 

        type(timer_class), intent(INOUT) :: timer 
        integer,           intent(IN)    :: comp 
        real(4),           intent(IN)    :: time_mod
        character(len=*),  intent(IN), optional :: label
        logical,           intent(IN), optional :: reset
        
        call timer_step_dble(timer,comp,real(time_mod,8),label,reset)

        return

    end subroutine timer_step_flt

    subroutine timer_step_dble(timer,comp,time_mod,label,reset)

        implicit none 

        type(timer_class), intent(INOUT) :: timer 
        integer,           intent(IN)    :: comp 
        real(8),           intent(IN)    :: time_mod
        character(len=*),  intent(IN), optional :: label
        logical,           intent(IN), optional :: reset
        
        ! Local variables
        integer :: n

        if (comp .eq. 0) then 
            ! Initialize time1

            ! Get current cpu [s] and model [units_mod] times
            call timer_cpu_time(timer%time_cpu1)
            timer%time_mod1 = time_mod

            if (present(reset)) then 
                if (reset) then 
                    timer%dtime_cpu_tot = 0.0
                    timer%dtime_mod_tot = 0.0
                    timer%nstep_tot     = 0
                end if 
            end if
        else
            ! Update timer for component number `comp`. 

            ! Get current cpu [s] and model [units_mod] times
            call timer_cpu_time(timer%time_cpu2) 
            timer%time_mod2 = time_mod

            ! Get time differences [s], [units_mod]
            timer%dtime_cpu(comp) = timer%time_cpu2 - timer%time_cpu1 
            timer%dtime_mod(comp) = timer%time_mod2 - timer%time_mod1 

            ! Add it to the total time 
            timer%dtime_cpu_tot(comp) = timer%dtime_cpu_tot(comp) + timer%dtime_cpu(comp)
            timer%dtime_mod_tot(comp) = timer%dtime_mod_tot(comp) + timer%dtime_mod(comp)

            ! Increment the total steps too 
            timer%nstep_tot(comp) = timer%nstep_tot(comp) + 1 

            ! Assign label to this component if available
            if (present(label)) timer%label(comp) = trim(label) 
            
            ! Finally, store time2 in time1 to be able to compute the next time interval as needed
            ! (alternatively, timer_step(...,comp=0) can be called to reset the timer)
            timer%time_cpu1 = timer%time_cpu2 
            timer%time_mod1 = timer%time_mod2 

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
        real(8) :: dtime(100)

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
        real(8) :: dtime(100), mtime(100)
        real(8) :: rate(100)

        ! Get dtime in desired units
        dtime = convert_dtime(timer%dtime_cpu,units)

        ! Get model time in model units
        mtime = timer%dtime_mod 

        where (dtime .ne. 0.0)
            rate = mtime / dtime
        elsewhere
            rate = 0.0
        end where

        units_rate = trim(units_mod)//"/"//trim(units)
        
        time_label = "calc_rate"
        if (present(label)) time_label = trim(label)

        write(timer%str_rate,"(a,1x,f12.3,1x,a)") trim(time_label), rate, "["//trim(units_rate)//"]"

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


    elemental function convert_dtime(dtime,units) result(dtime_conv)
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
                dtime_conv = -1.0d0
        end select

        return
        
    end function convert_dtime 

end module timer 
