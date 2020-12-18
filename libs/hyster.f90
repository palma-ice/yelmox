
module hyster 

    use nml 

    implicit none 

    integer,  parameter :: sp  = kind(1.0)
    integer,  parameter :: dp  = kind(1.0d0)
    integer,  parameter :: wp  = sp 

    real(wp), parameter :: MV  = -9999.0_wp 

    type hyster_par_class  
        character(len=56) :: label 
        real(wp) :: dt_ave 
        real(wp) :: df_sign 
        real(wp) :: dv_dt_scale
        real(wp) :: df_dt_min
        real(wp) :: df_dt_max
        real(wp) :: f_min
        real(wp) :: f_max   
    end type 

    type hyster_class

        type(hyster_par_class) :: par 

        ! variables 
        real(wp), allocatable :: time(:)
        real(wp), allocatable :: var(:)
        real(wp) :: dv_dt
        real(wp) :: df_dt
        
        real(wp) :: f_now
        logical  :: kill
    end type 

    private
    public :: wp
    public :: hyster_class
    public :: hyster_init 
    public :: hyster_calc_forcing

contains

    subroutine hyster_init(hyst,filename,time,label)

        type(hyster_class), intent(INOUT) :: hyst 
        character(len=*),   intent(IN)    :: filename 
        real(wp),           intent(IN)    :: time 
        character(len=*),   intent(IN), optional :: label 
        
        integer :: ntot 
        character(len=56) :: par_label 

        par_label = "hyster"
        if (present(label)) par_label = trim(par_label)//"_"//trim(label)

        ! Load parameters 
        call nml_read(filename,trim(par_label),"dt_ave",        hyst%par%dt_ave)
        call nml_read(filename,trim(par_label),"df_sign",       hyst%par%df_sign)
        call nml_read(filename,trim(par_label),"dv_dt_scale",   hyst%par%dv_dt_scale)
        call nml_read(filename,trim(par_label),"df_dt_min",     hyst%par%df_dt_min)
        call nml_read(filename,trim(par_label),"df_dt_max",     hyst%par%df_dt_max)
        
        call nml_read(filename,trim(par_label),"f_min",         hyst%par%f_min)
        call nml_read(filename,trim(par_label),"f_max",         hyst%par%f_max)
        
        ! Make sure sign is only +1/-1 
        hyst%par%df_sign = sign(1.0_wp,hyst%par%df_sign)

        ! Define label for this hyster object 
        hyst%par%label = "hyster" 
        if (present(label)) hyst%par%label = trim(hyst%par%label)//"_"//trim(label)

        ! (Re)initialize hyster vectors to a large value 
        ! to store many timesteps.
        ntot = 2000
        if (allocated(hyst%time)) deallocate(hyst%time)
        if (allocated(hyst%var))  deallocate(hyst%var)
        allocate(hyst%time(ntot))
        allocate(hyst%var(ntot))

        ! Initialize variable values
        hyst%time  = MV
        hyst%var   = MV 
        hyst%dv_dt = 0.0_wp 
        hyst%df_dt = 0.0_wp

        ! Initialize values of forcing 
        if (hyst%par%df_sign .gt. 0) then 
            hyst%f_now = hyst%par%f_min 
        else 
            hyst%f_now = hyst%par%f_max 
        end if 

        ! Set kill switch to false to start 
        hyst%kill = .FALSE. 

        return 

    end subroutine hyster_init

  
    subroutine hyster_calc_forcing (hyst,time,var)
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Subroutine :  d t T r a n s 1
        ! Author     :  Alex Robinson
        ! Purpose    :  Generate correct T_warming for gradual changes over
        !               time (continuous stability diagram!)
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  

        type(hyster_class), intent(INOUT) :: hyst 
        real(wp),           intent(IN)    :: time
        real(wp),           intent(IN)    :: var 

        ! Local variables 
        real(wp), allocatable :: dv_dt(:)
        real(wp) :: dv_dt_now 
        real(wp) :: f_scale 
        integer  :: ntot, kmin, kmax, nk, k 
        real(wp) :: dt_tot, dt_now 

        ! Get size of hyst vectors 
        ntot = size(hyst%time,1) 

        ! Get current timestep 
        dt_now = time - hyst%time(ntot) 

        ! Remove oldest point from beginning and add current one to the end
        hyst%time = eoshift(hyst%time,1,boundary=time)
        hyst%var  = eoshift(hyst%var, 1,boundary=var)

        ! Determine range of indices of times within our 
        ! time-averaging window. 
        kmin = findloc(hyst%time .ge. time - hyst%par%dt_ave,value=.TRUE., &
                                                     dim=1,mask=hyst%time.ne.MV)
        kmax = ntot 

        ! Determine currently available time window
        ! Note: do not use kmin here, in case time step does not match dt_ave,
        ! rather, use all available times in the vector to see if enough 
        ! time has passed. 
        dt_tot = hyst%time(kmax) - minval(hyst%time,mask=hyst%time.ne.MV)

        if (dt_tot .lt. hyst%par%dt_ave) then 
            ! Not enough time has passed, maintain constant forcing 
            ! (to avoid affects of noisy derivatives)

            hyst%df_dt = 0.0_wp 

        else 
            ! Calculate forcing 

            ! Get current number of averaging points 
            nk = kmax - kmin + 1 
            
            allocate(dv_dt(nk-1)) 

            ! Calculate mean rate of change over time steps of interest
            dv_dt = (hyst%var(kmin+1:kmax)-hyst%var(kmin:kmax-1)) / &
                      (hyst%time(kmin+1:kmax)-hyst%time(kmin:kmax-1))
            hyst%dv_dt = sum(dv_dt) / real(nk,wp)

            ! Calculate the current df_dt
            ! BASED ON EXPONENTIAL (sharp transition, tuneable)
            ! Returns scalar in range [0-1], 0.6 at dv_dt==dv_dt_scale
            f_scale = exp(-abs(hyst%dv_dt)/hyst%par%dv_dt_scale)

            ! Get forcing rate of change in [f/1e6 a]
            hyst%df_dt = hyst%par%df_sign * ( hyst%par%df_dt_min + f_scale*(hyst%par%df_dt_max-hyst%par%df_dt_min) )

            ! Convert [f/1e6 a] => [f/a]
            hyst%df_dt = hyst%df_dt *1e-6 
            
            ! Avoid underflow errors 
            if (abs(hyst%df_dt) .lt. 1e-10) hyst%df_dt = 0.0 

            ! Once the rate is available, update the current forcing value 
            hyst%f_now = hyst%f_now + hyst%df_dt * dt_now 
            
        end if 

        ! Check if kill should be activated 
        if (hyst%f_now .lt. hyst%par%f_min .or. hyst%f_now .gt. hyst%par%f_max) then 
            hyst%kill = .TRUE. 
        end if 

        return

    end subroutine hyster_calc_forcing

  
end module hyster
