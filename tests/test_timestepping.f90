program test_timestepping

    use nml
    use ncio
    use timestepping
    use yelmo_defs, only : wp
    

    type(tstep_class) :: ts
    integer  :: n, ntot
    real(wp) :: time_init, time_end
    real(wp) :: dtt
    logical  :: timesteps_complete

    time_init = -100.0
    time_end  =  1.0
    dtt       =  1.0

    call tstep_init(ts,time_init,method="bp",units="kyr")
    call tstep_print(ts)

    ! Advance timesteps
    timesteps_complete = .FALSE.

    do while (.not. timesteps_complete)

        ! Update timestep
        call tstep_update(ts,dtt)
        call tstep_print(ts)

        ! Check timesteps
        if (ts%time .ge. time_end) timesteps_complete = .TRUE. 

    end do

end program
