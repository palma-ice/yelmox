program test_timestepping

    use nml
    use ncio
    use timestepping
    use yelmo_defs, only : wp
    

    type(tstep_class) :: ts
    integer :: n, ntot
    real(wp) :: time_init, time_end
    real(wp) :: dtt
    logical  :: timesteps_complete

    time_init = 1800.0
    time_end  = 3000.0
    dtt       = 100.0

    call tstep_init(ts,time_init,time_pd=1950.0_wp,method="cal",units="years",trans_cal=.TRUE.,trans_bp=.FALSE.)
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
