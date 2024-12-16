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

    call tstep_init(ts,time_init,time_end,method="bp",units="kyr")
    call tstep_print(ts)

    ! Advance timesteps
    do while (.not. ts%is_finished)

        ! Update timestep
        call tstep_update(ts,dtt)
        call tstep_print(ts)

    end do

end program
