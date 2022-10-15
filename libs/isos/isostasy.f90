module isostasy 
    ! Isostasy (Greek ísos "equal", stásis "standstill") is the state of 
    ! gravitational equilibrium between Earth's crust and mantle such that 
    ! the crust "floats" at an elevation that depends on its thickness and density.
    ! -- https://en.wikipedia.org/wiki/Isostasy

    use nml 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: wp = sp 

    real(wp), parameter :: g  = 9.81                ! [m/s^2]
    real(wp), parameter :: pi = 3.14159265359

    real(wp), parameter :: rho_ice  = 910.0         ! [kg/m^3]
    real(wp), parameter :: rho_sw   = 1028.0        ! [kg/m^3] 
    real(wp), parameter :: rho_w    = 1000.0        ! [kg/m^3] 
    real(wp), parameter :: rho_a    = 3300.0        ! [kg/m^3] 3370 used by Coulon et al (2021)

    real(wp), parameter :: nu       = 0.25          ! [-]   Poisson's ratio, Coulon et al (2021) 
    real(wp), parameter :: E        = 100.0         ! [GPa] Young's modulus, Coulon et al (2021) 
    
    real(wp), parameter :: r_earth  = 6.378e6       ! [m]  Earth's radius, Coulon et al (2021) 
    real(wp), parameter :: m_earth  = 5.972e24      ! [kg] Earth's mass,   Coulon et al (2021) 
    
    type isos_param_class 
        integer            :: method            ! Type of isostasy to use
        character(len=512) :: fname_kelvin      ! File containing precalculated zero-order Kelvin function values
        real(wp)           :: dt                ! [yr] Timestep to recalculate bedrock uplift rate
        real(wp)           :: tau               ! [yr] Asthenospheric relaxation constant
        real(wp)           :: He_lith           ! [km] Effective elastic thickness of lithosphere
        real(wp)           :: D_lith            ! [N-m] Lithosphere flexural rigidity
        
        ! Internal parameters 
        real(wp) :: L_w         ! [m] flexural length scale
        integer  :: nrad        ! [-] Radius of neighborhood for convolution, in number of grid points        
        real(wp) :: time 

    end type 

    type isos_state_class 
        
        real(wp), allocatable :: z_bed(:,:)         ! Bedrock elevation         [m]
        real(wp), allocatable :: dzbdt(:,:)         ! Rate of bedrock uplift    [m/a]
        real(wp), allocatable :: z_bed_ref(:,:)     ! Reference (unweighted) bedrock 

        real(wp), allocatable :: kei(:,:)           ! Kelvin function filter values 
        real(wp), allocatable :: G0(:,:)            ! Green's function values

        real(wp), allocatable :: tau(:,:)           ! [yr] Asthenospheric relaxation timescale field
        real(wp), allocatable :: He_lith(:,:)       ! [m]  Effective elastic thickness of the lithosphere
        real(wp), allocatable :: D_lith(:,:)        ! [m]  [N-m] Lithosphere flexural rigidity
        
        real(wp), allocatable :: we(:,:)            ! enfoncement du socle autour
                                                    ! d'une load unitaire, sera
                                                    ! dimensionne dans le module 
                                                    ! isostasie_mod, routine init_iso a 
                                                    ! we(nrad:LBOC,-nrad:nrad)

        real(wp), allocatable :: w0(:,:)          ! Reference displacement
        real(wp), allocatable :: w1(:,:)          ! New displacement          

        
    end type 

    type isos_class
        type(isos_param_class) :: par
        type(isos_state_class) :: now 

    end type

    private
    public :: isos_class 
    public :: isos_init, isos_init_state 
    public :: isos_update 
    public :: isos_end  

    public :: isos_set_field

contains 

    subroutine isos_init(isos,filename,nx,ny,dx)

        implicit none 

        type(isos_class), intent(OUT) :: isos 
        character(len=*), intent(IN)  :: filename 
        integer,  intent(IN) :: nx, ny 
        real(wp), intent(IN) :: dx 

        ! Local variables
        integer :: n 

        ! Load parameters
        call isos_par_load(isos%par,filename,init=.TRUE.)
        
        ! Calculate the flexural length scale
        ! (Coulon et al, 2021, Eq. in text after Eq. 3)
        ! Note: should be on the order of 100km
        isos%par%L_w = (isos%par%D_lith / (rho_a*g))**0.25 

        ! nrad, 400 km on each side
        ! ajr: 400km set internally, since the radius should be smaller.
        ! This needs thorough revision and code refactoring. 
        ! See Greve and Blatter (2009), Chpt 8, page 192 for methodology 
        ! and Le Muer and Huybrechts (1996). It seems that this value
        ! should be larger to capture the forebuldge at 5-6x radius of relative stiffness
        !isos%par%nrad = int((400000.0-0.1)/dx)+1
        isos%par%nrad = int(4.0*isos%par%L_w/dx)+1
        
        ! Initialize isos variables 
        call isos_allocate(isos%now,nx,ny,nrad=isos%par%nrad)
        
        ! Intially ensure all variables are zero 
        isos%now%z_bed_ref  = 0.0
        isos%now%z_bed      = 0.0 
        isos%now%dzbdt      = 0.0 

        isos%now%we         = 0.0

        isos%now%w0         = 0.0   
        isos%now%w1         = 0.0   
        
        ! Set time to very large value in the future 
        isos%par%time       = 1e10 

        ! Calculate the Kelvin function filter 
        call calc_kei_filter_2D(isos%now%kei,dx=dx,dy=dx, &
                        L_w=isos%par%L_w,filename=isos%par%fname_kelvin)

        ! Calculate the Green's function values
        call calc_greens_function(isos%now%G0,isos%now%kei,isos%par%L_w,isos%par%D_lith)

        ! Store initial values of parameters as constant fields
        isos%now%He_lith    = isos%par%He_lith
        isos%now%D_lith     = isos%par%D_lith
        isos%now%tau        = isos%par%tau 
        

        ! Load lithospheric table of Kelvin function values 
        call read_tab_litho(isos%now%we,filename=trim(isos%par%fname_kelvin), &
                            RL=isos%par%L_w,DL=isos%par%D_lith,dx=dx,dy=dx)

        write(*,*) "isos_init:: range(we):  ", minval(isos%now%we),     maxval(isos%now%we)
        !write(*,*) "isos_init:: range(G0):  ", minval(isos%now%G0),     maxval(isos%now%G0)
        !stop 

        return 

    end subroutine isos_init

    subroutine isos_init_state(isos,z_bed,H_ice,z_sl,z_bed_ref,H_ice_ref,z_sl_ref,time)

        implicit none 

        type(isos_class), intent(INOUT) :: isos 
        real(wp), intent(IN) :: z_bed(:,:)            ! [m] Current bedrock elevation 
        real(wp), intent(IN) :: H_ice(:,:)            ! [m] Current ice thickness  
        real(wp), intent(IN) :: z_sl(:,:)             ! [m] Current sea level 
        real(wp), intent(IN) :: z_bed_ref(:,:)        ! [m] Reference bedrock elevation (with known load)
        real(wp), intent(IN) :: H_ice_ref(:,:)        ! [m] Reference ice thickness (associated with reference z_bed)
        real(wp), intent(IN) :: z_sl_ref(:,:)         ! [m] Reference sea level (associated with reference z_bed)
        real(wp), intent(IN) :: time                  ! [a] Initial time 
        
        ! Store reference bedrock field
        isos%now%z_bed_ref = z_bed_ref 
        isos%now%dzbdt    = 0.0 

        select case(isos%par%method)

            case(0,1) 
                ! 0: Steady-state lithospheric depression 
                ! 1: LLRA - Local lithosphere, relaxing Asthenosphere

                call calc_litho_local(isos%now%w0,z_bed_ref,H_ice_ref,z_sl_ref)
                
            case(2)
                ! 2: ELRA - Elastic lithosphere, relaxing Asthenosphere

                ! Local lithosphere (LL)
                call calc_litho_regional(isos%now%w0,z_bed_ref,H_ice_ref,z_sl_ref,isos%now%we)

        end select 

        ! Define initial time of isostasy model 
        isos%par%time = time 

        ! Store initial bedrock field 
        isos%now%z_bed = z_bed 

        ! Call isos_update to diagnose rate of change
        ! (no change to z_bed will be applied since isos%par%time==time)
        call isos_update(isos,H_ice,z_sl,time)

        write(*,*) "isos_init_state:: "
        write(*,*) "  Initial time:  ", isos%par%time 
        write(*,*) "  range(w0):     ", minval(isos%now%w0), maxval(isos%now%w0)
        write(*,*) "  range(w1):     ", minval(isos%now%w1), maxval(isos%now%w1)
        write(*,*) "  range(z_bed):  ", minval(isos%now%z_bed), maxval(isos%now%z_bed)
        
        return 

    end subroutine isos_init_state

    subroutine isos_update(isos,H_ice,z_sl,time)

        implicit none 

        type(isos_class), intent(INOUT) :: isos 
        real(wp), intent(IN) :: H_ice(:,:)        ! [m] Current ice thickness 
        real(wp), intent(IN) :: z_sl(:,:)         ! [m] Current sea level 
        real(wp), intent(IN) :: time              ! [a] Current time 

        ! Local variables 
        real(wp) :: dt 

        ! Step 0: determine current timestep 
        dt = time - isos%par%time 

        ! Step 1: diagnose rate of bedrock uplift
        
        select case(isos%par%method)

            case(0)
                ! Steady-state lithosphere

                isos%now%w1    = isos%now%w0
                isos%now%dzbdt = 0.0 

            case(1)
                ! Local lithosphere, relaxing asthenosphere (LLRA)

                ! Local lithosphere (LL)
                call calc_litho_local(isos%now%w1,isos%now%z_bed,H_ice,z_sl)

                ! Relaxing asthenosphere (RA)
                isos%now%dzbdt = ((isos%now%z_bed_ref-isos%now%z_bed) - (isos%now%w1-isos%now%w0))/isos%now%tau

            case(2)
                ! Elastic lithosphere, relaxing asthenosphere (ELRA)
                
                ! Regional elastic lithosphere (EL)
                call calc_litho_regional(isos%now%w1,isos%now%z_bed,H_ice,z_sl,isos%now%we)

                ! Relaxing asthenosphere (RA)
                isos%now%dzbdt = ((isos%now%z_bed_ref-isos%now%z_bed) - (isos%now%w1-isos%now%w0))/isos%now%tau
            
            case(3) 
                ! Elementary GIA model (spatially varying ELRA with geoid - to do!)

                ! Regional elastic lithosphere (EL)
                call calc_litho_regional(isos%now%w1,isos%now%z_bed,H_ice,z_sl,isos%now%we)

                ! Asthenosphere timescale field 
                isos%now%tau = isos%par%tau 

                ! Relaxing asthenosphere (RA)
                call calc_uplift_relax(isos%now%dzbdt,isos%now%z_bed,isos%now%z_bed_ref, &
                                                w_b=isos%now%w1-isos%now%w0,tau=isos%now%tau)

        end select 

        ! Step 2: update bedrock elevation (every timestep > 0)
        if (dt .ge. isos%par%dt) then 

            isos%now%z_bed = isos%now%z_bed + isos%now%dzbdt*dt 

            isos%par%time  = time 

        end if 
            
        return 

    end subroutine isos_update

    subroutine isos_end(isos)

        implicit none 

        type(isos_class), intent(INOUT) :: isos 

        ! Deallocate isos object 
        call isos_deallocate(isos%now)

        return 

    end subroutine isos_end

    subroutine isos_par_load(par,filename,init)

        implicit none

        type(isos_param_class), intent(OUT) :: par
        character(len=*),       intent(IN)  :: filename 
        logical, optional        :: init

        ! Local variables 
        logical :: init_pars 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 
        
        call nml_read(filename,"isostasy","method",         par%method,         init=init_pars)
        call nml_read(filename,"isostasy","fname_kelvin",   par%fname_kelvin,   init=init_pars)
        call nml_read(filename,"isostasy","dt",             par%dt,             init=init_pars)
        call nml_read(filename,"isostasy","He_lith",        par%He_lith,        init=init_pars)
        call nml_read(filename,"isostasy","D_lith",         par%D_lith,         init=init_pars)
        call nml_read(filename,"isostasy","tau",            par%tau,            init=init_pars)
        
        return

    end subroutine isos_par_load

    subroutine isos_allocate(now,nx,ny,nrad)

        implicit none 

        type(isos_state_class), intent(INOUT) :: now 
        integer, intent(IN) :: nx, ny, nrad  

        ! Local variables
        integer :: nfilt 

        nfilt = 2*nrad+1 

        ! First ensure arrays are not allocated
        call isos_deallocate(now)

        ! Allocate arrays

        allocate(now%kei(nfilt,nfilt))
        allocate(now%G0(nfilt,nfilt))
        
        allocate(now%z_bed(nx,ny))
        allocate(now%dzbdt(nx,ny))
        allocate(now%z_bed_ref(nx,ny))
        
        
        allocate(now%tau(nx,ny))
        allocate(now%He_lith(nx,ny))
        allocate(now%D_lith(nx,ny))

        allocate(now%we(nfilt,nfilt))
        
        allocate(now%w0(nx,ny))
        allocate(now%w1(nx,ny))
        
        return 

    end subroutine isos_allocate

    subroutine isos_deallocate(now)

        implicit none 

        type(isos_state_class), intent(INOUT) :: now 

        if (allocated(now%kei))         deallocate(now%kei)
        if (allocated(now%G0))          deallocate(now%G0)
        
        if (allocated(now%z_bed))       deallocate(now%z_bed)
        if (allocated(now%dzbdt))       deallocate(now%dzbdt)
        if (allocated(now%z_bed_ref))   deallocate(now%z_bed_ref)
        
        if (allocated(now%tau))         deallocate(now%tau)
        if (allocated(now%He_lith))     deallocate(now%He_lith)
        if (allocated(now%D_lith))      deallocate(now%D_lith)
        
        if (allocated(now%we))          deallocate(now%we)
        
        if (allocated(now%w0))          deallocate(now%w0)
        if (allocated(now%w1))          deallocate(now%w1)
        
        
        return 

    end subroutine isos_deallocate


    subroutine isos_set_field(var,var_values,mask_values,mask,dx,sigma)
        ! Impose multiple var values according to the corresponding
        ! locations provided in a mask. 
        ! Additionally impose Gaussian smoothing via the
        ! smoothing radius sigma.

        implicit none 

        real(wp), intent(OUT) :: var(:,:) 
        real(wp), intent(IN)  :: var_values(:)
        real(wp), intent(IN)  :: mask_values(:)
        real(wp), intent(IN)  :: mask(:,:) 
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: sigma 

        ! Local variables
        integer :: j, n 

        ! Safety check 
        if (sigma .le. dx) then 
            write(*,*) "isos_set_field:: Error: sigma must be larger than dx."
            write(*,*) "dx    = ", dx 
            write(*,*) "sigma = ", sigma 
            stop 
        end if 

        ! Determine how many values should be assigned
        n = size(var_values,1)

        ! Initially set var=0 everywhere
        var = 0.0_wp 

        ! Loop over unique var values and assign them
        ! to correct regions of domain.
        do j = 1, n 

            where(mask .eq. mask_values(j)) var = var_values(j)

        end do

        ! Apply Gaussian smoothing as desired
        call smooth_gauss_2D(var,dx=dx,sigma=sigma)
        
        return

    end subroutine isos_set_field

    subroutine smooth_gauss_2D(var,dx,sigma,mask_apply,mask_use)
        ! Smooth out a field to avoid noise 
        ! mask_apply designates where smoothing should be applied 
        ! mask_use   designates which points can be considered in the smoothing filter 

        implicit none

        real(wp),   intent(INOUT) :: var(:,:)      ! [nx,ny] 2D variable
        real(wp),   intent(IN)    :: dx 
        real(wp),   intent(IN)    :: sigma  
        logical,    intent(IN), optional :: mask_apply(:,:) 
        logical,    intent(IN), optional :: mask_use(:,:) 

        ! Local variables
        integer  :: i, j, nx, ny, n, n2, k 
        integer  :: imx, ipx, jmx, jpx 
        real(wp), allocatable :: filter0(:,:), filter(:,:) 
        real(wp), allocatable :: var_old(:,:) 
        logical,  allocatable :: mask_apply_local(:,:) 
        logical,  allocatable :: mask_use_local(:,:) 

        real(wp), allocatable :: var_ext(:,:), var_ref_ext(:,:) 

        nx    = size(var,1)
        ny    = size(var,2)

        ! Determine half-width of filter as 2-sigma
        n2 = ceiling(2.0_wp*sigma / dx)

        ! Get total number of points for filter window in each direction
        n = 2*n2+1
        
        allocate(var_old(nx,ny))
        allocate(mask_apply_local(nx,ny))
        allocate(mask_use_local(nx,ny))
        allocate(filter0(n,n))
        allocate(filter(n,n))

        allocate(var_ext(-n2:nx+n2,-n2:ny+n2))
        allocate(var_ref_ext(-n2:nx+n2,-n2:ny+n2))
        
        ! Check whether mask_apply is available 
        if (present(mask_apply)) then 
            ! use mask_use to define neighborhood points
            
            mask_apply_local = mask_apply 

        else
            ! Assume that everywhere should be smoothed

            mask_apply_local = .TRUE.
        
        end if

        ! Check whether mask_use is available 
        if (present(mask_use)) then 
            ! use mask_use to define neighborhood points
            
            mask_use_local = mask_use 

        else
            ! Assume that mask_apply also gives the points to use for smoothing 

            mask_use_local = mask_apply_local
        
        end if

        ! Calculate default 2D Gaussian smoothing kernel
        filter0 = gauss_values(dx,dx,sigma=sigma,n=n)

        var_old = var 

        var_ref_ext = -9999.0 

        var_ref_ext(1:nx,1:ny) = var 
        do i = 0, -n2, -1
            k = -i+1
            var_ref_ext(i,:) = var_ref_ext(k,:) 
        end do 
        do i = nx+1, nx+n2 
            k = nx + ((nx+1)-i)
            var_ref_ext(i,:) = var_ref_ext(k,:) 
        end do 

        do j = 0, -n2, -1
            k = -j+1
            var_ref_ext(:,j) = var_ref_ext(:,k) 
        end do 
        do j = ny+1, ny+n2 
            k = ny + ((ny+1)-j)
            var_ref_ext(:,j) = var_ref_ext(:,k) 
        end do 

        if (count(var_ref_ext .eq. -9999.0) .gt. 0) then 
            write(*,*) "Missing points!"
            stop 
        end if 

        do j = 1, ny
        do i = 1, nx


            !if (mask_apply_local(i,j)) then 
                ! Apply smoothing to this point 

                ! Limit filter input to neighbors of interest
                filter = filter0 
                !where(.not. mask_use_local(i-n2:i+n2,j-n2:j+n2) ) filter = 0.0

                ! If neighbors are available, normalize and perform smoothing  
                if (sum(filter) .gt. 0.0) then 
                    filter = filter/sum(filter)
                    var_ext(i,j) = sum(var_ref_ext(i-n2:i+n2,j-n2:j+n2)*filter) 
                end if  

            !end if 

        end do 
        end do 

        ! Get variable on normal grid 
        var = var_ext(1:nx,1:ny)

        return 

    end subroutine smooth_gauss_2D

    function gauss_values(dx,dy,sigma,n) result(filt)
        ! Calculate 2D Gaussian smoothing kernel
        ! https://en.wikipedia.org/wiki/Gaussian_blur

        implicit none 

        real(wp), intent(IN) :: dx 
        real(wp), intent(IN) :: dy 
        real(wp), intent(IN) :: sigma 
        integer,  intent(IN) :: n 
        real(wp) :: filt(n,n) 

        ! Local variables 
        real(wp) :: x, y  
        integer  :: n2, i, j, i1, j1  

        if (mod(n,2) .ne. 1) then 
            write(*,*) "gauss_values:: error: n can only be odd."
            write(*,*) "n = ", n 
        end if 

        n2 = (n-1)/2 

        do j = -n2, n2 
        do i = -n2, n2 
            x = i*dx 
            y = j*dy 

            i1 = i+1+n2 
            j1 = j+1+n2 
            filt(i1,j1) = 1.0/(2.0*pi*sigma**2)*exp(-(x**2+y**2)/(2*sigma**2))

        end do 
        end do 
        
        ! Normalize to ensure sum to 1
        filt = filt / sum(filt)

        return 

    end function gauss_values

    subroutine calc_greens_function(G0,kei2D,L_w,D_lith)

        implicit none

        real(wp), intent(OUT) :: G0(:,:) 
        real(wp), intent(IN)  :: kei2D(:,:) 
        real(wp), intent(IN)  :: L_w 
        real(wp), intent(IN)  :: D_lith 

        G0 = -L_w**2 / (2.0*pi*D_lith) * kei2D

        return

    end subroutine calc_greens_function


    subroutine calc_kei_filter_2D(filt,dx,dy,L_w,filename)
        ! Calculate 2D Kelvin function (kei) smoothing kernel

        implicit none 

        real(wp), intent(OUT) :: filt(:,:) 
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dy 
        real(wp), intent(IN)  :: L_w  
        character(len=*), intent(IN) :: filename 

        ! Local variables 
        integer  :: i, j, i1, j1, n, n2
        real(wp) :: x, y, r

        real(wp), allocatable :: rn_vals(:) 
        real(wp), allocatable :: kei_vals(:) 
        
        ! Get size of filter array and half-width
        n  = size(filt,1) 
        n2 = (n-1)/2 

        ! Safety check
        if (size(filt,1) .ne. size(filt,2)) then 
            write(*,*) "calc_kei_filt:: error: array 'filt' must be square [n,n]."
            write(*,*) "size(filt): ", size(filt,1), size(filt,2)
            stop
        end if 

        ! Safety check
        if (mod(n,2) .ne. 1) then 
            write(*,*) "calc_kei_filt:: error: n can only be odd."
            write(*,*) "n = ", n
            stop  
        end if 

        ! Load tabulated kei values from file
        call load_kei_values(rn_vals,kei_vals,filename)

        ! Loop over filter array in two dimensions,
        ! calculate the distance from the center, normalized by L_w
        ! and impose correct Kelvin function value. 

        do j = -n2, n2 
        do i = -n2, n2

            x  = i*dx 
            y  = j*dy 
            r  = sqrt(x**2+y**2)

            ! Get actual index of array
            i1 = i+1+n2 
            j1 = j+1+n2 

            ! Get correct kei value for this point
            filt(i1,j1) = get_kei_value(r,L_w,rn_vals,kei_vals)

        end do 
        end do 
        
        return 

    end subroutine calc_kei_filter_2D

    function get_kei_value(r,L_w,rn_vals,kei_vals) result(kei)

        implicit none

        real(wp), intent(IN) :: r           ! [m] Radius from point load 
        real(wp), intent(IN) :: L_w         ! [m] Flexural length scale
        real(wp), intent(IN) :: rn_vals(:)  ! [-] Tabulated normalised radius values (r/L_w)
        real(wp), intent(IN) :: kei_vals(:) ! [-] Tabulated normalised Kelvin function values
        real(wp) :: kei

        ! Local variables 
        integer :: k, n 
        real(wp) :: rn_now 
        real(wp) :: wt 

        n = size(rn_vals,1) 

        ! Get current normalized radius from point load
        rn_now = r / L_w 

        if (rn_now .gt. rn_vals(n)) then 

            kei = kei_vals(n)

        else 

            do k = 1, n-1
                if (rn_now .ge. rn_vals(k) .and. rn_now .lt. rn_vals(k+1)) exit
            end do 

            ! Linear interpolation to get current kei value
            kei = kei_vals(k) &
                + (rn_now-rn_vals(k))/(rn_vals(k+1)-rn_vals(k))*(kei_vals(k+1)-kei_vals(k))   

        end if 

        ! Diagnostics (only works if function is changed to a subroutine!)
        !write(*,*) "get_kei_value: ", L_w, r, rn_now, k, rn_vals(k), kei_vals(k) 

        return

    end function get_kei_value

    subroutine load_kei_values(rn,kei,filename)

        implicit none

        real(wp), allocatable,  intent(OUT):: rn(:)
        real(wp), allocatable,  intent(OUT):: kei(:)
        character(len=*),       intent(IN) :: filename

        ! Local variables 
        integer :: k 
        integer, parameter :: ntot    = 1001
        integer, parameter :: filenum = 177

        if (allocated(rn))  deallocate(rn)
        if (allocated(kei)) deallocate(kei) 

        allocate(rn(ntot))
        allocate(kei(ntot)) 

        ! fonction de kelvin
        ! lecture de la table kei qui est tous les 0.01 entre 0 et 10
        ! STEPK=100=1/ecart 
        !cdc modification du chemin maintenant fonction de dir_inp
        ! trim(dir_inp)//'kelvin.res'
        open(filenum,file=trim(filename))
        read(filenum,*)  ! Skip first line
        do k = 1, ntot
            read(filenum,*) rn(k),kei(k)
        end do
        close(filenum)

        return

    end subroutine load_kei_values

    function calc_kei_value(r,L_w) result(kei)
        ! This function is based on Greve and Blatter (2009), Eq. 8.34.
        ! Combined with 8.11, it should be possible to obtain
        ! an analytical expression for the Kelvin function (kei). 
        ! This could eventually be used as a replacement for loading 
        ! the tabulated values from a file. 

        ! So far, this approach is not giving the right answer. 
        
        implicit none 

        real(wp), intent(IN) :: r       ! [m] Radius from point load 
        real(wp), intent(IN) :: L_w     ! [m] Flexural length scale
        real(wp) :: kei 

        ! Local variables
        real(wp) :: alpha 
        real(wp) :: f_now 
        real(wp) :: fac 

        alpha = r / (sqrt(2.0)*L_w)
        f_now = exp(-alpha)*(cos(alpha)+sin(alpha))

        fac = sqrt(2.0)**3 * L_w * (pi/4.0)
        
        kei = fac * f_now 

        ! Note: this doesn't give the right values of kei!!!

        return

    end function calc_kei_value

    subroutine read_tab_litho(we,filename,rl,dl,dx,dy) 
        !        subroutine qui donne la repartition des enfoncements
        !        en fonction de la rigidite de la lithosphere.
        !        definition du tableau
        !
        !   variables en entree-----------
        !      ROM masse volumique du manteau
        !      rl  rayon de rigidite relative
        !      dl  rigidite flexural
        !  
        !   variables en sortie------------
        !      we  deflection due a une load unitaire      
        !
        !
        !

        implicit none 

        real(wp), intent(INOUT) :: we(:,:) 
        character(len=*), intent(IN) :: filename 
        real(wp), intent(IN) :: rl, dl, dx, dy 

        ! Local variables
        integer, parameter :: nk = 1000
        integer :: i, j, k, nr 
        real(wp) :: kei(0:nk)
        real(wp) :: stepk, al, xl, dist
        integer :: num_kelvin = 177 

        real(wp), allocatable :: we00(:,:)

        ! Size of neighborhood 
        nr = (size(we,1)-1)/2 

        ! Allocate local WE object with proper indexing
        allocate(we00(-nr:nr,-nr:nr))

        ! pour la lithosphere
        stepk = 100.0
        al    = -rl*rl/(2.0*pi*dl)*dx*dy      

        ! fonction de kelvin
        ! lecture de la table kei qui est tous les 0.01 entre 0 et 10
        ! STEPK=100=1/ecart 
        !cdc modification du chemin maintenant fonction de dir_inp
        ! trim(dir_inp)//'kelvin.res'
        open(num_kelvin,file=trim(filename))
        read(num_kelvin,*)  ! Skip first line
        do k=1,nk
            read(num_kelvin,*) xl, kei(k)
        end do
        close(num_kelvin)

        do i = -nr, nr    
        do j = -nr, nr
            dist = dx*sqrt(1.0*(i*i+j*j))                                
            xl   = dist/rl*stepk                                          
            k    = int(xl)+1
            if ((k.gt.834).or.(dist.gt.dx*nr)) then                 
                we00(i,j) = 0.0                                              
            else 
                we00(i,j)=kei(k)+(kei(k+1)-kei(k))*(xl-k)
                if (K.eq.834) we00(i,j) = max(we00(i,j),0.0)                   
            endif
            we00(i,j) = we00(i,j)*al 

            ! write(*,*) i, j, dist, xl, k, we00(i,j) 
            ! write(*,*) kei(k), kei(k+1), xl-k
            ! if (i .eq. 0 .and. j .eq. 0) stop 

        end do
        end do

        ! normalisation and scaling
        we00  = we00/sum(we00) / (rho_a*g)

        ! Return solution to external object
        we = we00(-nr:nr,-nr:nr)

        ! Make sure too small values are eliminated to avoid underflow errors
        where(abs(we) .lt. 1e-15) we = 0.0_wp 

        return
    
    end subroutine read_tab_litho

    ! === isos physics routines ======================================

    elemental subroutine calc_litho_local(w1,z_bed,H_ice,z_sl)
        ! Calculate the local lithospheric loading from ice or ocean weight 
        ! in units of [N]

        implicit none 

        real(wp), intent(INOUT) :: w1
        real(wp), intent(IN)    :: z_bed, H_ice, z_sl 

        if (rho_ice*H_ice.ge.rho_sw*(z_sl-z_bed)) then
            ! Ice or land

            w1 = rho_ice*g*H_ice

        else
            ! Ocean

            w1 = rho_sw*g*(z_sl-z_bed)

        end if

        return 

    end subroutine calc_litho_local

    subroutine calc_litho_regional(w1,z_bed,H_ice,z_sl,we)
        ! Calculate the load on the lithosphere as
        ! distributed on an elastic plate. 

        implicit none

        real(wp), intent(INOUT) :: w1(:,:)      ! [N] Lithospheric load
        real(wp), intent(IN)    :: z_bed(:,:)   ! [m] Bed elevation
        real(wp), intent(IN)    :: H_ice(:,:)   ! [m] Ice thickness 
        real(wp), intent(IN)    :: z_sl(:,:)    ! [m] Sea level 
        real(wp), intent(IN)    :: we(:,:)      ! Regional filter function
        
        ! Calculate local lithospheric load first
        call calc_litho_local(w1,z_bed,H_ice,z_sl)

        ! Convolve the estimated point load with the regional
        ! filter to obtain the distributed load w1. 
        call convolve_load_elastic_plate(w1,we)

        return

    end subroutine calc_litho_regional

    subroutine convolve_load_elastic_plate(w1,we)
        ! litho-0.3.f            10 Novembre 1999             *     
        !
        ! Petit routine qui donne la repartition des enfoncements
        ! en fonction de la rigidite de la lithosphere.
        !
        ! En entree 
        !      ------------
        !     we(-nrad:nrad,-nrad:nrad)  : deflection due a une load unitaire 
        !                          defini dans tab-litho

        implicit none

        real(wp), intent(INOUT) :: w1(:,:)      ! [N] Lithospheric loading
        real(wp), intent(IN)    :: we(:,:)      ! Regional scaling filter

        ! Local variables
        integer :: IP,JP,LPX,LPY,II,SOM1,SOM2
        real(wp), allocatable :: w1_ext(:,:)
        real(wp), allocatable :: w_reg(:,:)
        real(wp), allocatable :: croix(:)

        integer :: i, j, nx ,ny, nr

        nx = size(w1,1)
        ny = size(w1,2)

        ! Size of regional neighborhood 
        nr = (size(we,1)-1)/2 

        ! Populate load on extended grid
        allocate(w1_ext(1-nr:nx+nr,1-nr:ny+nr))

        ! First fill in main grid points with current point load
        w1_ext(1:nx,1:ny) = w1 

        ! Populate the extended grid points
        do i = 1, nx
            w1_ext(i,1-nr:0)=w1_ext(i,1)
            w1_ext(i,ny+1:ny+nr)=w1_ext(i,ny)
        end do
        do j = 1, ny
            w1_ext(1-nr:0,j)=w1_ext(1,j)
            w1_ext(NX+1:NX+nr,j)=w1_ext(nx,j)
        end do
        
        ! Populate the extended grid corner points     
        w1_ext(1-nr:0,1-nr:0)         = w1_ext(1,1)
        w1_ext(1-nr:0,ny+1:ny+nr)     = w1_ext(1,ny)
        w1_ext(nx+1:nx+nr,1-nr:0)     = w1_ext(nx,1)
        w1_ext(nx+1:nx+nr,ny+1:ny+nr) = w1_ext(nx,ny)

        ! ----- allocation de w_reg  et de croix -----------

        allocate(w_reg(-nr:nr,-nr:nr))
        allocate(croix(0:nr))

        ! calcul de la deflexion par sommation des voisins appartenant
        ! au bloc de taille 2*nr
        som1 = 0.0
        som2 = 0.0

        ! On somme aussi les contributions des points exterieurs au domaine
        ! lorsque la charge est due a l'ocean. On suppose alors  que
        ! ces points ont la meme charge que les limites

        do j = 1, ny
        do i = 1, nx

            w1(i,j) = 0.0 
            ii      = 0

            ! Apply the neighborhood weighting to the load
            w_reg = we * w1_ext(i-nr:i+nr,j-nr:j+nr)

            ! Sum of all effects (tentative to avoid rounding errors)
            w1(i,j) = w_reg(0,0)

            ! dans croix on calcul la somme des effets WLOC puis on met le resultat dans W1   
            do lpx = 1, nr

                lpy = 0
                croix(lpy) = (w_reg(lpx,lpy)+w_reg(-lpx,lpy))

                croix(lpy) = (w_reg(lpy,lpx)+w_reg(lpy,-lpx)) + croix(lpy) 
                lpy = lpx
                croix(lpy) = ((w_reg(lpx,lpy)+w_reg(lpx,-lpy)) &
                              +(w_reg(-lpx,lpy)+w_reg(-lpx,-lpy))) 
                do lpy=1,lpx-1
                    croix(lpy) = (((w_reg(lpx,lpy)+w_reg(lpx,-lpy)) &
                                  +(w_reg(-lpx,lpy)+w_reg(-lpx,-lpy))) &
                                  +((w_reg(lpy,lpx)+w_reg(lpy,-lpx)) &
                                  +(w_reg(-lpy,lpx)+w_reg(-lpy,-lpx))))
                end do

                do lpy = 0, lpx 
                    w1(i,j) = w1(i,j) + croix(lpy) ! sommation dans W1
                end do
     
            end do

            ! --- FIN DE L'INTEGRATION SUR LE PAVE nr
            som1 = som1 + w1(i,j)
            som2 = som2 - w1_ext(i,j)/(rho_a*g)

        end do
        end do
   
        return

    end subroutine convolve_load_elastic_plate
    
    elemental subroutine calc_uplift_relax(dzbdt,z_bed,z_bed_ref,w_b,tau)
        ! Calculate rate of change of vertical bedrock height
        ! from a relaxing asthenosphere.

        implicit none

        real(wp), intent(OUT) :: dzbdt 
        real(wp), intent(IN)  :: z_bed 
        real(wp), intent(IN)  :: z_bed_ref
        real(wp), intent(IN)  :: w_b        ! w_b = w1-w0
        real(wp), intent(IN)  :: tau

        dzbdt = -((z_bed-z_bed_ref) - w_b) / tau
        
        return

    end subroutine calc_uplift_relax


end module isostasy
