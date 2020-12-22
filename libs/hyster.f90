
module hyster 

    use nml 

    implicit none 

    integer,  parameter :: sp  = kind(1.0)
    integer,  parameter :: dp  = kind(1.0d0)
    integer,  parameter :: wp  = sp 

    real(wp), parameter :: MV  = -9999.0_wp 

    type hyster_par_class  
        character(len=56) :: label 
        character(len=56) :: method 
        real(wp) :: dt_ave 
        real(wp) :: df_sign 
        real(wp) :: dv_dt_scale
        real(wp) :: df_dt_min
        real(wp) :: df_dt_max
        real(wp) :: df_dt_const 
        real(wp) :: f_range(2)  
        real(wp) :: pc_eps 

        ! Internal parameters 
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
        
        real(wp) :: pc_df(3)
        real(wp) :: pc_eta(3)

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
        call nml_read(filename,trim(par_label),"method",        hyst%par%method)
         
        call nml_read(filename,trim(par_label),"dt_ave",        hyst%par%dt_ave)
        call nml_read(filename,trim(par_label),"df_sign",       hyst%par%df_sign)
        call nml_read(filename,trim(par_label),"dv_dt_scale",   hyst%par%dv_dt_scale)
        call nml_read(filename,trim(par_label),"df_dt_min",     hyst%par%df_dt_min)
        call nml_read(filename,trim(par_label),"df_dt_max",     hyst%par%df_dt_max)
        call nml_read(filename,trim(par_label),"df_dt_const",   hyst%par%df_dt_const)
        call nml_read(filename,trim(par_label),"f_range",       hyst%par%f_range)
        call nml_read(filename,trim(par_label),"pc_eps",        hyst%par%pc_eps)
        
        ! Make sure sign is only +1/-1 
        hyst%par%df_sign = sign(1.0_wp,hyst%par%df_sign)

        ! Get f_min and f_max values 
        hyst%par%f_min = hyst%par%f_range(1) 
        hyst%par%f_max = hyst%par%f_range(2) 
        
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

        hyst%pc_df  = hyst%par%df_dt_min 
        hyst%pc_eta = hyst%par%pc_eps

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
        real(wp) :: dvdt_fac 
        real(wp) :: pc_df_now 

        ! Since dv_dt is typically calculated over an averaging period,
        ! assume second-order PI controller parameters are needed. 
        integer, parameter :: pc_order = 2

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
            hyst%dv_dt = sum(dv_dt,mask= (hyst%time(kmin+1:kmax) .ne. MV) .and. &
                                         (hyst%time(kmin:kmax-1) .ne. MV)) / real(nk,wp)


            select case(trim(hyst%par%method))

                case("const") 
                    ! Apply a constant rate of change, independent of dv_dt
                    ! Set magnitude of df_dt right now in [f/(1e6 yr)] 
                    hyst%df_dt = hyst%par%df_dt_const

                case("exp")

                    ! Calculate the current df_dt
                    ! BASED ON EXPONENTIAL (sharp transition, tuneable)
                    ! Returns scalar in range [0-1], 0.6 at dv_dt==dv_dt_scale
                    ! Note: apply limit to dvdt_fac of a maximum value of 10, so 
                    ! that exp function doesn't explode (exp(-10)=>0.0)
                    dvdt_fac = min(abs(hyst%dv_dt)/hyst%par%dv_dt_scale,10.0_wp)
                    f_scale  = exp(-dvdt_fac)

                    ! Get forcing rate of change magnitude in [f/1e6 a]
                    hyst%df_dt = ( hyst%par%df_dt_min + f_scale*(hyst%par%df_dt_max-hyst%par%df_dt_min) )

                case("PI42","H312b","H312PID","H321PID","PID1")

                    ! Calculate adaptive dfdt value using proportional-integral (PI) methods
                    call set_adaptive_timestep_pc(pc_df_now,hyst%pc_df,hyst%pc_eta,hyst%par%pc_eps, &
                                        hyst%par%df_dt_min,hyst%par%df_dt_max,pc_order,hyst%par%method)

                    ! Remove oldest point from the end and insert latest point in beginning
                    hyst%pc_eta = eoshift(hyst%pc_eta,-1,boundary=abs(hyst%dv_dt))
                    hyst%pc_df  = eoshift(hyst%pc_df, -1,boundary=abs(pc_df_now))

                    ! Apply limits to eta so that algorithm works well. 
                    ! pc_eta should be greater than zero
                    hyst%pc_eta(1) = max(hyst%pc_eta(1),1e-3)

                    ! Get forcing rate of change magnitude in [f/1e6 a]
                    hyst%df_dt = hyst%pc_df(1) 

            end select 

            ! Convert [f/(1e6 yr)] => [f/yr] and apply sign of change
            hyst%df_dt = hyst%par%df_sign*hyst%df_dt *1e-6 

            ! Avoid underflow errors 
            if (abs(hyst%df_dt) .lt. 1e-8) hyst%df_dt = 0.0 

            ! Once the rate is available, update the current forcing value 
            hyst%f_now = hyst%f_now + (hyst%df_dt*dt_now) 
            
        end if 

        ! Check if kill should be activated 
        if (hyst%f_now .lt. hyst%par%f_min .or. hyst%f_now .gt. hyst%par%f_max) then 
            hyst%kill = .TRUE. 
        end if 

        return

    end subroutine hyster_calc_forcing
            
    subroutine set_adaptive_timestep_pc(dt_new,dt,eta,eps,dtmin,dtmax,pc_k,controller)
        ! Calculate the timestep following algorithm for 
        ! a general predictor-corrector (pc) method.
        ! Implemented followig Cheng et al (2017, GMD)

        implicit none 

        real(wp), intent(OUT) :: dt_new               ! [yr]   Timestep (n+1)
        real(wp), intent(IN)  :: dt(:)                ! [yr]   Timesteps (n:n-2)
        real(wp), intent(IN)  :: eta(:)               ! [X/yr] Maximum truncation error (n:n-2)
        real(wp), intent(IN)  :: eps                  ! [--]   Tolerance value (eg, eps=1e-4)
        real(wp), intent(IN)  :: dtmin                ! [yr]   Minimum allowed timestep, must be > 0
        real(wp), intent(IN)  :: dtmax                ! [yr]   Maximum allowed timestep
        integer,    intent(IN)  :: pc_k                 ! pc_k gives the order of the timestepping scheme (pc_k=2 for FE-SBE, pc_k=3 for AB-SAM)
        character(len=*), intent(IN) :: controller      ! Adaptive controller to use [PI42, H312b, H312PID]

        ! Local variables
        real(wp) :: dt_n, dt_nm1, dt_nm2          ! [yr]   Timesteps (n:n-2)
        real(wp) :: eta_n, eta_nm1, eta_nm2       ! [X/yr] Maximum truncation error (n:n-2)
        real(wp) :: rho_n, rho_nm1, rho_nm2
        real(wp) :: rhohat_n 
        real(wp) :: dt_adv 
        real(wp) :: dtmax_now
        real(wp) :: k_i 
        real(wp) :: k_p, k_d 

        ! Smoothing parameter; Söderlind and Wang (2006) method, Eq. 10
        ! Values on the order of [0.7,2.0] are reasonable. Higher kappa slows variation in dt
        real(wp), parameter :: kappa = 2.0_wp 
        
        ! Step 1: Save information needed for adapative controller algorithms 

        ! Save dt from several timesteps (potentially more available)
        dt_n    = max(dt(1),dtmin) 
        dt_nm1  = max(dt(2),dtmin) 
        dt_nm2  = max(dt(3),dtmin)

        ! Save eta from several timesteps (potentially more available)
        eta_n   = eta(1)
        eta_nm1 = eta(2)
        eta_nm2 = eta(3)

        ! Calculate rho from several timesteps 
        rho_nm1 = (dt_n   / dt_nm1) 
        rho_nm2 = (dt_nm1 / dt_nm2) 

        ! Step 2: calculate scaling for the next timestep (dt,n+1)
        select case(trim(controller))

            case("PI42")
                ! Söderlind and Wang, 2006; Cheng et al., 2017
                ! Deeper discussion in Söderlind, 2002. Note for example, 
                ! that Söderlind (2002) recommends:
                ! k*k_i >= 0.3 
                ! k*k_p >= 0.2 
                ! k*k_i + k*k_p <= 0.7 
                ! However, the default values of Söderlind and Wang (2006) and Cheng et al (2017)
                ! are outside of these bounds. 

                ! Default parameter values 
                k_i = 2.0_wp / (pc_k*5.0_wp)
                k_p = 1.0_wp / (pc_k*5.0_wp)
                
                ! Improved parameter values (reduced oscillations)
!                 k_i = 4.0_wp / (pc_k*10.0_wp)
!                 k_p = 3.0_wp / (pc_k*10.0_wp)

                ! Experimental parameter values (minimal oscillations, does not access largest timesteps)
!                 k_i = 0.5_wp / (pc_k*10.0_wp)
!                 k_p = 6.5_wp / (pc_k*10.0_wp)

                ! Default parameter values
                rho_n = calc_pi_rho_pi42(eta_n,eta_nm1,rho_nm1,eps,k_i,k_p,alpha_2=0.0_wp)

            case("H312b") 
                ! Söderlind (2003) H312b, Eq. 31+ (unlabeled) 
                
                rho_n = calc_pi_rho_H312b(eta_n,eta_nm1,eta_nm2,rho_nm1,rho_nm2,eps,k=real(pc_k,wp),b=8.0_wp)

            case("H312PID") 
                ! Söderlind (2003) H312PD, Eq. 38
                ! Note: Suggested k_i =(2/9)*1/pc_k, but lower value gives more stable solution

                !k_i = (2.0_wp/9.0_wp)*1.0_wp/real(pc_k,wp)
                k_i = 0.08_wp/real(pc_k,wp)

                rho_n = calc_pi_rho_H312PID(eta_n,eta_nm1,eta_nm2,eps,k_i)

            case("H321PID")

                k_i = 0.1  / real(pc_k,wp)
                k_p = 0.45 / real(pc_k,wp) 

                rho_n = calc_pi_rho_H321PID(eta_n,eta_nm1,eta_nm2,dt_n,dt_nm1,eps,k_i,k_p)
                
            case("PID1")

                k_i = 0.175 
                k_p = 0.075
                k_d = 0.01 

                rho_n = calc_pi_rho_PID1(eta_n,eta_nm1,eta_nm2,eps,k_i,k_p,k_d)
                
            case DEFAULT 

                write(*,*) "set_adaptive_timestep_pc:: Error: controller not recognized."
                write(*,*) "controller = ", trim(controller) 
                stop 

        end select 

        ! Scale rho_n for smoothness 
        rhohat_n = rho_n
        !rhohat_n = min(rho_n,1.1)
        !rhohat_n = 1.0_wp + kappa * atan((rho_n-1.0_wp)/kappa) ! Söderlind and Wang, 2006, Eq. 10
        
        ! Step 3: calculate the next time timestep (dt,n+1)
        dt_new = rhohat_n * dt_n

        ! Ensure timestep is also within parameter limits 
        dt_new = max(dtmin,dt_new)  ! dt >= dtmin
        dt_new = min(dtmax,dt_new)  ! dt <= dtmax

        return 

    end subroutine set_adaptive_timestep_pc

    function calc_pi_rho_pi42(eta_n,eta_nm1,rho_nm1,eps,k_i,k_p,alpha_2) result(rho_n)

        implicit none 

        real(wp), intent(IN) :: eta_n 
        real(wp), intent(IN) :: eta_nm1 
        real(wp), intent(IN) :: rho_nm1 
        real(wp), intent(IN) :: eps 
        real(wp), intent(IN) :: k_i 
        real(wp), intent(IN) :: k_p
        real(wp), intent(IN) :: alpha_2 
        real(wp) :: rho_n 

        ! Söderlind and Wang, 2006; Cheng et al., 2017
        ! Original formulation: Söderlind, 2002, Eq. 3.12:
        rho_n   = (eps/eta_n)**(k_i+k_p) * (eps/eta_nm1)**(-k_p) * rho_nm1**(-alpha_2)

        return 

    end function calc_pi_rho_pi42

    function calc_pi_rho_H312b(eta_n,eta_nm1,eta_nm2,rho_nm1,rho_nm2,eps,k,b) result(rho_n)

        implicit none 

        real(wp), intent(IN) :: eta_n 
        real(wp), intent(IN) :: eta_nm1 
        real(wp), intent(IN) :: eta_nm2 
        real(wp), intent(IN) :: rho_nm1
        real(wp), intent(IN) :: rho_nm2
        real(wp), intent(IN) :: eps 
        real(wp), intent(IN) :: k 
        real(wp), intent(IN) :: b 
        real(wp) :: rho_n 

        ! Local variables 
        real(wp) :: beta_1, beta_2, beta_3 
        real(wp) :: alpha_2, alpha_3 

        beta_1  =  1.0_wp / (k*b)
        beta_2  =  2.0_wp / (k*b)
        beta_3  =  1.0_wp / (k*b)
        alpha_2 = -3.0_wp / b 
        alpha_3 = -1.0_wp / b 

        ! Söderlind (2003) H312b, Eq. 31+ (unlabeled) 
        rho_n   = (eps/eta_n)**beta_1 * (eps/eta_nm1)**beta_2 * (eps/eta_nm2)**beta_3 &
                            * rho_nm1**alpha_2 * rho_nm2**alpha_3 

        return 

    end function calc_pi_rho_H312b

    function calc_pi_rho_H312PID(eta_n,eta_nm1,eta_nm2,eps,k_i) result(rho_n)

        implicit none 

        real(wp), intent(IN) :: eta_n 
        real(wp), intent(IN) :: eta_nm1 
        real(wp), intent(IN) :: eta_nm2 
        real(wp), intent(IN) :: eps 
        real(wp), intent(IN) :: k_i  
        real(wp) :: rho_n 

        ! Local variables 
        real(wp) :: k_i_1, k_i_2, k_i_3

        k_i_1   = k_i / 4.0_wp 
        k_i_2   = k_i / 2.0_wp 
        k_i_3   = k_i / 4.0_wp 

        ! Söderlind (2003) H312PID, Eq. 38
        rho_n   = (eps/eta_n)**k_i_1 * (eps/eta_nm1)**k_i_2 * (eps/eta_nm2)**k_i_3

        return 

    end function calc_pi_rho_H312PID

    function calc_pi_rho_H321PID(eta_n,eta_nm1,eta_nm2,dt_n,dt_nm1,eps,k_i,k_p) result(rho_n)

        implicit none 

        real(wp), intent(IN) :: eta_n 
        real(wp), intent(IN) :: eta_nm1 
        real(wp), intent(IN) :: eta_nm2 
        real(wp), intent(IN) :: dt_n 
        real(wp), intent(IN) :: dt_nm1 
        real(wp), intent(IN) :: eps 
        real(wp), intent(IN) :: k_i  
        real(wp), intent(IN) :: k_p
        real(wp) :: rho_n 

        ! Local variables 
        real(wp) :: k_i_1, k_i_2, k_i_3

        k_i_1   =   0.75_wp*k_i + 0.50_wp*k_p 
        k_i_2   =   0.50_wp*k_i 
        k_i_3   = -(0.25_wp*k_i + 0.50_wp*k_p)

        ! Söderlind (2003) H321PID, Eq. 42
        rho_n   = (eps/eta_n)**k_i_1 * (eps/eta_nm1)**k_i_2 * (eps/eta_nm2)**k_i_3 * (dt_n / dt_nm1)

        return 

    end function calc_pi_rho_H321PID

    function calc_pi_rho_PID1(eta_n,eta_nm1,eta_nm2,eps,k_i,k_p,k_d) result(rho_n)

        implicit none 

        real(wp), intent(IN) :: eta_n 
        real(wp), intent(IN) :: eta_nm1 
        real(wp), intent(IN) :: eta_nm2 
        real(wp), intent(IN) :: eps 
        real(wp), intent(IN) :: k_i  
        real(wp), intent(IN) :: k_p
        real(wp), intent(IN) :: k_d
        real(wp) :: rho_n 

        ! https://www.mathematik.uni-dortmund.de/~kuzmin/cfdintro/lecture8.pdf
        ! Page 20 (theoretical basis unclear/unknown)
        rho_n   = (eps/eta_n)**k_i * (eta_nm1/eta_n)**k_p * (eta_nm1**2/(eta_n*eta_nm2))**k_d

        return 

    end function calc_pi_rho_PID1

end module hyster
