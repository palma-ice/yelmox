module pico
    ! PICO module

    use nml 
    use ncio 
    use pico_geometry
    use pico_physics

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: prec = sp 

    ! Physical constants 
    real(prec), parameter :: rho_ice =  917.d0       ! Density ice           [kg/m^3] 
    real(prec), parameter :: rho_w   = 1000.d0       ! Density water         [kg/m^3] 
    real(prec), parameter :: rho_sw  = 1028.d0       ! Density seawater      [kg/m^3] 
    real(prec), parameter :: g       = 9.81d0        ! Gravitational accel.  [m/s^2]
    real(prec), parameter :: cp_o    = 3974.d0       ! Specific heat capacity of ocean mixed layer [J/kg*ºC]
    real(prec), parameter :: L_ice   = 3.34e5        ! Latent heat of fusion of ice [J/kg] 
 
    real(prec), parameter :: year_to_sec = 365.0*24.0*60.0*60.0
    real(prec), parameter :: lambda      = L_ice/cp_o ! K
    real(prec), parameter :: rho_ice_sw  = rho_ice / rho_sw 

    type pico_param_class

        logical    :: use_pico
        real(prec) :: n_box
        real(prec) :: a_pico, b_pico, c_pico
        real(prec) :: alpha_pico, beta_pico
        real(prec) :: rho_star
        real(prec) :: gamma_tstar
        real(prec) :: C_over
        real(prec) :: c_deep, depth_deep
        logical    :: find_ocean

        character(len=512) :: domain

    end type

    type pico_state_class 
        
        real(prec), allocatable   :: d_shlf(:,:)            ! Shelf distance to grl
        real(prec), allocatable   :: d_if(:,:)              ! Shelf distance to ice front
        real(prec), allocatable   :: r_shlf(:,:)            ! Relative ice shelf position
        real(prec), allocatable   :: n_shlf(:,:)            ! Max number of boxes per basin
        real(prec), allocatable   :: A_box(:,:)             ! Area of box
        real(prec), allocatable   :: boxes(:,:)             ! Box number of ice shelf
        real(prec), allocatable   :: bmb_shlf(:,:)          ! Shelf basal mass balance [m/a]
        real(prec), allocatable   :: T_box(:,:)             ! Temperature of box [ºC]
        real(prec), allocatable   :: S_box(:,:)             ! Salinity of box [PSU]
        real(prec), allocatable   :: CC(:,:)                ! Overtuning strength

        integer,    allocatable   :: mask_ocn_ref(:,:) 
        integer,    allocatable   :: mask_ocn(:,:) 
       
    end type 

    type pico_class
        type(pico_param_class) :: par 
        type(pico_state_class) :: now 
    end type

    private
    public :: pico_class
    public :: pico_init
    public :: pico_calc_geometry
    public :: pico_update_physics
    public :: pico_end

contains
    ! =================
    ! 
    ! Parameters load
    ! 
    ! =================

    subroutine pico_init(pico,filename,nx,ny,domain,grid_name,regions)

        implicit none

        type(pico_class), intent(OUT)     :: pico
        character(len=*), intent(IN)      :: filename
        integer, intent(IN)               :: nx, ny 
        character(len=*), intent(IN)      :: domain, grid_name
        real(prec), intent(IN)            :: regions(:,:)
        ! Load parameters
        call pico_par_load(pico%par,domain,filename)
        ! Allocate the object 
        call pico_allocate(pico%now,nx,ny)

        ! Define mask_ocn_ref based on regions mask 
        ! (these definitions should work for all North and Antarctica domains)
        pico%now%mask_ocn_ref = 0 
        where (regions .eq. 1.0_prec) pico%now%mask_ocn_ref = 1   
        where (regions .eq. 2.0_prec) pico%now%mask_ocn_ref = 1 


        select case(trim(pico%par%domain))

            case("Greenland") 
                ! Greenland specific ocean kill regions

                where (regions .ne. 1.3) pico%now%mask_ocn_ref = 2 

                ! ajr: not used anymore now with mask_ocn_ref formulation,
                ! keeping code here just to see if Greenland domain still calculated well.
                ! ! Kill regions that should not be calculated (for now)
                ! ! North America and Ellesmere Island (1.1,1.11)
                ! ! Svalbard (1.2,1.23)
                ! ! Iceland (1.31)
                ! ! open sea (1.0)
                ! where (regions .eq. 1.1 .or. regions .eq. 1.11) is_c_deep = .TRUE.
                ! where (regions .eq. 1.2 .or. regions .eq. 1.23) is_c_deep = .TRUE.
                ! where (regions .eq. 1.31) is_c_deep = .TRUE.
                ! where (regions .eq. 1.0)  is_c_deep = .TRUE.

            case("Antarctica") 
                ! Antarctica specific ocean kill regions

                ! Omit regions==2.11 which means c_deep is not applied to deep points within continental shelf
                where (regions .ne. 2.11) pico%now%mask_ocn_ref = 2 

            ! case("North") 
            !     ! North specific ocean kill regions

            !     ! Apply only in purely open-ocean regions  
            !     where (regions .eq. 1.0) pico%now%mask_ocn_ref = 2 

            case DEFAULT 
                ! Other domains: c_deep potentially applied everywhere 
                ! with deep ocean points 

                pico%now%mask_ocn_ref = 2 

        end select 

        ! ====================================
        ! 
        ! Summary
        ! 
        ! ====================================
        
        write(*,*) "range mask_ocn_ref: ", minval(pico%now%mask_ocn_ref), maxval(pico%now%mask_ocn_ref)
        
        ! Initialize variables 
        pico%now%bmb_shlf      = 0.0  
        pico%now%T_box         = 0.0
        pico%now%S_box         = 0.0

        return

    end subroutine pico_init

    ! Load PICO variables
    subroutine pico_par_load(par,domain,filename,init)

        type(pico_param_class), intent(OUT) :: par
        character(len=*), intent(IN) :: domain,filename
        logical, optional :: init
        logical :: init_pars

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE.

        call nml_read(filename,"pico","use_pico",    par%use_pico,    init=init_pars)
        call nml_read(filename,"pico","n_box",       par%n_box,       init=init_pars)
        call nml_read(filename,"pico","a_pico",      par%a_pico,      init=init_pars)
        call nml_read(filename,"pico","b_pico",      par%b_pico,      init=init_pars)
        call nml_read(filename,"pico","c_pico",      par%c_pico,      init=init_pars)
        call nml_read(filename,"pico","alpha_pico",  par%alpha_pico,  init=init_pars)
        call nml_read(filename,"pico","beta_pico",   par%beta_pico,   init=init_pars)
        call nml_read(filename,"pico","rho_star",    par%rho_star,    init=init_pars)
        call nml_read(filename,"pico","gamma_tstar", par%gamma_tstar, init=init_pars)
        call nml_read(filename,"pico","C_over",      par%C_over,      init=init_pars)
        call nml_read(filename,"pico","c_deep",      par%c_deep,      init=init_pars)
        call nml_read(filename,"pico","depth_deep",  par%depth_deep,  init=init_pars)
        call nml_read(filename,"pico","find_ocean",  par%find_ocean,  init=init_pars)

        ! Determine derived parameters
        par%domain = trim(domain)

        ! convert T_fp from C to K
        par%b_pico     = par%b_pico + 273.15
        ! 1Sv = 10^6 m3/s. Then 1/s to 1/yr
        par%C_over     = par%C_over*year_to_sec*10**6

        return

    end subroutine pico_par_load

    subroutine pico_end(pico)

        implicit none 

        type(pico_class) :: pico 

        ! Deallocate pico state object
        call pico_deallocate(pico%now)

        return 

    end subroutine pico_end

    subroutine pico_allocate(now,nx,ny)

        implicit none 

        type(pico_state_class) :: now 
        integer :: nx, ny 

        ! Make object is deallocated
        call pico_deallocate(now)

        ! Allocate pico 
        allocate(now%T_box(nx,ny))
        allocate(now%S_box(nx,ny))
        allocate(now%d_shlf(nx,ny))
        allocate(now%n_shlf(nx,ny))
        allocate(now%d_if(nx,ny))
        allocate(now%r_shlf(nx,ny))
        allocate(now%A_box(nx,ny))
        allocate(now%bmb_shlf(nx,ny))
        allocate(now%CC(nx,ny))
        allocate(now%boxes(nx,ny))

        allocate(now%mask_ocn_ref(nx,ny))
        allocate(now%mask_ocn(nx,ny))

        ! Initialize variables 
        now%T_box    = 0.0
        now%S_box    = 0.0
        now%d_shlf   = 0.0
        now%d_if     = 0.0
        now%r_shlf   = 0.0
        now%n_shlf   = 0.0
        now%A_box    = 0.0
        now%bmb_shlf = 0.0
        now%CC       = 0.0
        now%boxes    = 0.0

        ! By default set ocean points everywhere
        now%mask_ocn_ref = 1
        now%mask_ocn     = 1

        return

    end subroutine pico_allocate

    subroutine pico_deallocate(now)

        implicit none 

        type(pico_state_class) :: now 

        ! Allocate state objects
        if (allocated(now%T_box))    deallocate(now%T_box)
        if (allocated(now%S_box))    deallocate(now%S_box)
        if (allocated(now%d_shlf))   deallocate(now%d_shlf)
        if (allocated(now%d_if))     deallocate(now%d_if)
        if (allocated(now%r_shlf))   deallocate(now%r_shlf)
        if (allocated(now%n_shlf))   deallocate(now%n_shlf)
        if (allocated(now%A_box))    deallocate(now%A_box)
        if (allocated(now%bmb_shlf))    deallocate(now%bmb_shlf)
        if (allocated(now%CC))       deallocate(now%CC)
        if (allocated(now%boxes))    deallocate(now%boxes)        

        return

    end subroutine pico_deallocate

    subroutine parse_path(path,domain,grid_name)

        implicit none 

        character(len=*), intent(INOUT) :: path 
        character(len=*), intent(IN)    :: domain, grid_name 

        call nml_replace(path,"{domain}",   trim(domain))
        call nml_replace(path,"{grid_name}",trim(grid_name))
        
        return 

    end subroutine parse_path

    subroutine pico_calc_geometry(pico,H_ice,f_grnd,basins,dx)

        implicit none

        type(pico_class), intent(INOUT) :: pico
        real(prec), intent(IN) :: H_ice(:,:), f_grnd(:,:), basins(:,:)
        real(prec), intent(IN) :: dx

        ! Local variables
        logical, allocatable :: is_grline(:,:), is_margin(:,:)
        real(prec), allocatable :: boxes_basin(:,:)
        real(prec) :: d_max, d_max_basin
        integer :: m

        allocate(is_grline(size(f_grnd,1),size(f_grnd,2)))
        allocate(is_margin(size(f_grnd,1),size(f_grnd,2)))
        allocate(boxes_basin(size(f_grnd,1),size(f_grnd,2)))
        is_grline = calc_grline(f_grnd,H_ice)
        is_margin = calc_margin(H_ice)
        d_max = 0.0
        d_max_basin = 0.0

        ! compute ice shelf distance to margin and grounding-line
        call pico_calc_shelf_extent(pico%now%d_shlf,pico%now%d_if,pico%now%r_shlf,H_ice,f_grnd,is_grline,is_margin,dx)
        d_max = maxval(pico%now%d_shlf)

        do m = 1 , maxval(basins)
            d_max_basin = maxval(pico%now%d_shlf,basins .eq. m)
            boxes_basin = 0.0
            where(basins .eq. m)
                ! Compute boxes and area
                ! jablasco: create first indep. box basin mask seems to work
                boxes_basin = pico_calc_shelf_boxes(H_ice,f_grnd,pico%now%d_shlf,pico%now%r_shlf,pico%par%n_box,d_max,d_max_basin)
                pico%now%boxes = boxes_basin
                pico%now%A_box = pico_calc_area_boxes(boxes_basin,f_grnd,pico%par%n_box,dx)
            end where
            d_max_basin = 0.0
        end do        

    end subroutine pico_calc_geometry

    subroutine pico_update_physics(pico,to,so,H_ice,z_bed,f_grnd,z_sl)
        
        implicit none
        
        type(pico_class), intent(INOUT) :: pico
        real(prec), intent(IN) :: to(:,:),so(:,:),H_ice(:,:) 
        real(prec), intent(IN) :: z_bed(:,:) 
        real(prec), intent(IN) :: f_grnd(:,:)
        real(prec), intent(IN) :: z_sl(:,:) 

        ! Local variables
        integer :: m, i, j, k, l, nx, ny, ngr

        real(prec) :: bmb_floating, pm_point, pm_point_box0
        real(prec), allocatable :: H_ocn(:,:)
        real(prec), allocatable :: T_star(:,:)
        logical,    allocatable :: is_grline(:,:)

        nx = size(f_grnd,1)
        ny = size(f_grnd,2) 

        allocate(H_ocn(nx,ny)) 
        allocate(is_grline(nx,ny))
        allocate(T_star(nx,ny))

        ! Determine location of grounding line 
        is_grline = calc_grline(f_grnd,H_ice)
 
        ! Initialize Tstar, Tbox, Sbox
        T_star = 0.0
        pico%now%T_box  = to
        pico%now%S_box  = so

        ! PISM condition in PICO. Box 0 temperature cannot be below the pressure-melting point.
        do l = 1, ny
        do k = 1, nx
            pm_point_box0 = calc_theta_pm(pico%now%S_box(k,l),pico%par%a_pico,pico%par%b_pico,pico%par%c_pico,H_ice(k,l))
            if (pico%now%T_box(k,l) .lt. pm_point_box0) pico%now%T_box(k,l) = pm_point_box0 + 0.001
        end do
        end do

        where(f_grnd .eq. 1.0)
            pico%now%T_box = 0.0
            pico%now%S_box = 0.0
        end where

        if (pico%par%find_ocean) then 
            ! Determine which floating or partially floating points
            ! are connected to the open-ocean

            call find_open_ocean(pico%now%mask_ocn,f_grnd,pico%now%mask_ocn_ref)

        else 
            ! Set any floating or partially floating points to ocean points

            pico%now%mask_ocn = 0 
            where (f_grnd .lt. 1.0) pico%now%mask_ocn = 1
            where (f_grnd .lt. 1.0 .and. pico%now%mask_ocn_ref .eq. 2) pico%now%mask_ocn = 2

        end if

        ! 1. Calculate current ice shelf bmb field (grounded-ice bmb is
        ! calculated in ice-sheet model separately) ========
        
        do j = 1, ny
        do i = 1, nx

            ! Calculate ocean depth
            H_ocn(i,j) = max( (z_sl(i,j)-z_bed(i,j))-(rho_ice_sw*H_ice(i,j)),0.0 )
             
            if (pico%now%mask_ocn(i,j) .gt. 0) then
                do m = 1, pico%par%n_box
                    ! jablasco: test box 1 -> .eq. -> .ge. la m
                    if (m .eq. 1.0 .and. pico%now%boxes(i,j) .ge. 1.0) then
                        ! 1. Compute box 1
                        T_star(i,j) = calc_Tstar(to(i,j),so(i,j),H_ice(i,j),pico%par%a_pico,pico%par%b_pico,pico%par%c_pico)
                        pico%now%T_box(i,j) = calc_shelf_Tbox_1(to(i,j),so(i,j),pico%now%A_box(i,j),T_star(i,j), &
                                                                pico%par%C_over,pico%par%rho_star,pico%par%gamma_tstar, &
                                                                pico%par%alpha_pico,pico%par%beta_pico)
                        pico%now%S_box(i,j) = calc_shelf_Sbox_1(to(i,j),so(i,j),pico%now%T_box(i,j))
                        ! Compute overtuning
                        pico%now%CC(i,j) = calc_overtuning(to(i,j),so(i,j),pico%now%T_box(i,j),pico%now%S_box(i,j), &
                                                           pico%par%C_over,pico%par%rho_star,pico%par%alpha_pico,pico%par%beta_pico)             

                    ! jablasco:test box 1. 1.0 -> n_box
                    else if (m .gt. 1.0 .and. pico%now%boxes(i,j) .ge. m) then
                        ! Compute rest boxes
                        T_star(i,j) = calc_Tstar(pico%now%T_box(i,j),pico%now%S_box(i,j),H_ice(i,j), &
                                                 pico%par%a_pico,pico%par%b_pico,pico%par%c_pico)
                        call calc_shelf_TS_box_n(pico%now%T_box(i,j), pico%now%S_box(i,j), pico%now%A_box(i,j), &
                                                 T_star(i,j), pico%now%CC(i,j), pico%par%a_pico, pico%par%gamma_tstar)

                    end if

                end do

                ! Compute melting
                pm_point = calc_theta_pm(pico%now%S_box(i,j),pico%par%a_pico,pico%par%b_pico,pico%par%c_pico,H_ice(i,j))
                bmb_floating = calc_melt_rate_pico(pico%now%T_box(i,j),pm_point,pico%par%gamma_tstar)
               
                ! Apply melting to purely floating points or ocean
                if(f_grnd(i,j) .eq. 0.0) pico%now%bmb_shlf(i,j) = bmb_floating

                ! Ensure that ice accretion only occurs where ice exists 
                if (pico%now%bmb_shlf(i,j) .gt. 0.0 .and. H_ice(i,j) .eq. 0.0) pico%now%bmb_shlf(i,j) = 0.0 

                ! No accretion allowed in  Box 1
                if (pico%now%bmb_shlf(i,j) .gt. 0.0 .and. pico%now%boxes(i,j) .eq. 1.0) pico%now%bmb_shlf(i,j) = 0.0

                ! Note: the following condition is a good idea, but in grisli-ucm the grounding line
                ! is defined as the last *grounded* point, so this limit sets bmb at grounding line to zero
                ! ! Ensure that refreezing rate is less than the available water depth within one time step 
                ! pico%now%bmb_shlf = max(pico%now%bmb_shlf,-H_ocn/dt)
                
            else 
                ! Grounded point, or floating point not connected to the ocean 
                ! Set bmb_shlf to zero 

                pico%now%bmb_shlf(i,j) = 0.0
                
            end if 
 
        end do
        end do  

        call apply_c_deep(pico%now%bmb_shlf,pico%par%c_deep,pico%par%depth_deep,pico%now%mask_ocn,z_bed,z_sl,n_smth=0) 

        return
        
    end subroutine pico_update_physics

    subroutine find_open_ocean(mask,f_grnd,mask_ref)
        ! Brute-force routine to find all ocean points 
        ! (ie, when f_grnd < 1) connected to
        ! the open ocean as defined in mask_ref. 

        implicit none 

        integer,    intent(OUT) :: mask(:,:) 
        real(prec), intent(IN)  :: f_grnd(:,:)
        integer,    intent(IN)  :: mask_ref(:,:) 

        ! Local variables 
        integer :: i, j, q, nx, ny
        integer :: im1, ip1, jm1, jp1 
        integer :: n_unfilled 

        integer, parameter :: qmax = 1000 

        nx = size(mask,1)
        ny = size(mask,2) 

        ! First specify land points (mask==0)
        ! and assume all floating points are 'closed ocean' points (mask==-1)
        mask = 0 
        where (f_grnd .lt. 1.0) mask = -1 

        ! Now populate our ocean mask to be consistent 
        ! with known open ocean points that can be 
        ! normal open ocean (1) or deep ocean (2) 
        ! For now set all points to 1 for easier looping
        where (mask .eq. -1 .and. mask_ref .gt. 0) mask = 1 
         
        ! Iteratively fill in open-ocean points that are found
        ! next to known open-ocean points
        do q = 1, qmax 

            n_unfilled = 0 

            do j = 1, ny 
            do i = 1, nx 

                ! Get neighbor indices 
                im1 = max(i-1,1)
                ip1 = min(i+1,nx)
                jm1 = max(j-1,1)
                jp1 = min(j+1,ny)
                
                if (mask(i,j) .eq. 1) then 
                    ! This is an open-ocean point 
                    ! Define any neighbor ocean points as open-ocean points
                    
                    if (mask(im1,j) .eq. -1) then 
                        mask(im1,j) = 1
                        n_unfilled = n_unfilled + 1
                    end if 

                    if (mask(ip1,j) .eq. -1) then 
                        mask(ip1,j) = 1 
                        n_unfilled = n_unfilled + 1
                    end if

                    if (mask(i,jm1) .eq. -1) then 
                        mask(i,jm1) = 1 
                        n_unfilled = n_unfilled + 1
                    end if 

                    if (mask(i,jp1) .eq. -1) then 
                        mask(i,jp1) = 1 
                        n_unfilled = n_unfilled + 1
                    end if 

                end if
                    
            end do 
            end do  

            !write(*,*) q, n_unfilled, count(mask .eq. -1) 

            ! Exit loop if no more open-ocean points are found 
            if (n_unfilled .eq. 0) exit 

        end do 

        ! Finally populate deep ocean points 
        where (mask .gt. 0 .and. mask_ref .eq. 2) mask = 2 
        
        return 

    end subroutine find_open_ocean

end module pico
