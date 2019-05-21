module timing

    implicit none 

    type timer_class
        real (8) :: dtime1, dtime2 
    end type 

contains

    subroutine timer_step(timer,step)

        implicit none 

        type(timer_class) :: timer 
        integer :: step 

        if (step .eq. 1) then 
            call cpu_time(timer%dtime1)           ! get current time in seconds   
        else if (step .eq. 2) then 
            call cpu_time(timer%dtime2)           ! get current time in seconds   
        else
            write(*,*) "timer_step:: unknown step: ", step 
            stop 
        end if 

        return 

    end subroutine timer_step 

    subroutine timer_print(timer,units,label)

        implicit none 

        type(timer_class) :: timer
        character(len=1)  :: units 
        character(len=*), optional :: label 
        character(len=512) :: time_label 

        double precision :: dtime 

        dtime = timer%dtime2 - timer%dtime1 

        select case(units)

            case("s")
                dtime = dtime 

            case("m")
                dtime = dtime / 60.d0 

            case("h") 
                dtime = dtime / (60.d0*60.d0)

            case DEFAULT
                write(*,*) "timer_print:: error: units not recognized: "//units

        end select

        time_label = "Calculation time: "
        if (present(label)) time_label = trim(label)//": "

        write(*,"(a,f18.3,1x,a1)") trim(time_label), dtime, units

        return 

    end subroutine timer_print




end module timing 
