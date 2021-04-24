

module marine_shelf
    ! Module to simulate the marine-shelf interface:
    ! Calculates the basal mass balance of an ice shelf (bmb_shlf)

    use nml 
    use ncio 

!     use yelmo_defs, only : sp, dp, wp, rho_ice, rho_w, rho_sw, g, parse_path 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the working precision of the library (sp,dp)
    integer,  parameter :: wp = sp 

    ! Physical constants 
    real(wp), parameter :: rho_ice =  917.d0     ! Density ice           [kg/m^3] 
    real(wp), parameter :: rho_w   = 1000.d0     ! Density water         [kg/m^3] 
    real(wp), parameter :: rho_sw  = 1028.d0     ! Density seawater      [kg/m^3] 
    real(wp), parameter :: g       = 9.81d0      ! Gravitational accel.  [m/s^2]
    real(wp), parameter :: cp_o    = 3974.d0     ! Specific heat capacity of ocean mixed layer [J/kg*ºC]
    real(wp), parameter :: L_ice   = 3.34e5      ! Latent heat of fusion of ice [J/kg] 
 
    real(wp), parameter :: lambda     = L_ice/cp_o
    real(wp), parameter :: rho_ice_sw = rho_ice / rho_sw 


    type marshelf_param_class

        character(len=56)  :: bmb_shlf_method, T_shlf_method
        logical            :: use_obs
        character(len=512) :: obs_path
        character(len=56)  :: obs_name 
        real(wp)         :: obs_f, obs_lim 
        character(len=56)  :: basin_name(50)
        real(wp)         :: basin_bmb(50)
        real(wp) :: c_shlf, kappa_shlf, f_grz_shlf
        real(wp) :: c_grz, kappa_grz, grz_length
        real(wp) :: gamma_lin, gamma_quad, gamma_quad_nl 

        real(wp) :: c_deep, depth_deep
        real(wp) :: T_fp 
        real(wp) :: depth_min, depth_max

        logical :: find_ocean 

        character(len=512) :: domain   
        real(wp) :: rho_ice, rho_sw

    end type 

    type marshelf_state_class 
        real(wp), allocatable :: bmb_shlf(:,:)          ! Shelf basal mass balance [m/a]
        real(wp), allocatable :: bmb_obs(:,:)           ! observed ice shelf melting [m/a]
        real(wp), allocatable :: bmb_ref(:,:)           ! Basal mass balance reference field
        real(wp), allocatable :: T_shlf(:,:)            ! Boundary ocean temps. for forcing the ice shelves
        real(wp), allocatable :: T_basin(:,:)           ! Basin average boundary ocean temps. for forcing the ice shelves
        real(wp), allocatable :: dT_shlf(:,:)           ! Boundary ocean temp anomalies for forcing the ice shelves
        real(wp), allocatable :: dT_basin(:,:)          ! Basin average boundary ocean temp anomalies for forcing the ice shelves
        real(wp), allocatable :: S_shlf(:,:)            ! Boundary salinity for forcing the ice shelves
        real(wp), allocatable :: kappa(:,:)             ! Shelf-melt coefficient [m/a/K]

        real(wp), allocatable :: z_base(:,:)            ! Ice-shelf base elevation (relative to sea level)
        real(wp), allocatable :: slope_base(:,:)        ! Ice-shelf base slope (slope=sin(theta)=length/hypotenuse)

        integer,  allocatable :: mask_ocn_ref(:,:) 
        integer,  allocatable :: mask_ocn(:,:) 
        
    end type 

    type marshelf_class
        type(marshelf_param_class) :: par 
        type(marshelf_state_class) :: now 
    end type

    private
    public :: marshelf_class
    public :: marshelf_init
    public :: marshelf_calc_Tshlf
    public :: marshelf_calc_Sshlf
    public :: marshelf_set_Tshlf 
    public :: marshelf_update
    public :: marshelf_end 

contains 

    subroutine marshelf_init(mshlf,filename,nx,ny,domain,grid_name,regions,basins)

        implicit none 

        type(marshelf_class), intent(OUT) :: mshlf
        character(len=*), intent(IN)      :: filename
        integer, intent(IN)               :: nx, ny 
        character(len=*), intent(IN)      :: domain, grid_name
        real(wp), intent(IN)            :: regions(:,:)
        real(wp), intent(IN)            :: basins(:,:)
        
        ! Local variables
        integer :: j 
        integer :: num

        ! Load parameters
        call marshelf_par_load(mshlf%par,filename,domain,grid_name)

        ! Allocate the object 
        call marshelf_allocate(mshlf%now,nx,ny)
        
        ! Decide whether to load observations or use psuedo-observations
        if (mshlf%par%use_obs) then

            ! Load observed field from file
            call nc_read(mshlf%par%obs_path,mshlf%par%obs_name,mshlf%now%bmb_obs) 

            ! Make negative since bmb is defined with melt negative, but observed field has melt positive
            mshlf%now%bmb_obs = -mshlf%now%bmb_obs

            ! Scale the obs as needed 
            mshlf%now%bmb_obs = mshlf%now%bmb_obs*mshlf%par%obs_f 

            ! Limit the max obs values as desired 
            where (mshlf%now%bmb_obs >  mshlf%par%obs_lim) mshlf%now%bmb_obs =  mshlf%par%obs_lim 
            where (mshlf%now%bmb_obs < -mshlf%par%obs_lim) mshlf%now%bmb_obs = -mshlf%par%obs_lim 
        
        else 

            ! Set default psuedo-observation value since none are available and make them negative (as above) 
            mshlf%now%bmb_obs = -1.0*mshlf%par%obs_f 

        end if 

        ! Make domain specific initializations
        select case(trim(mshlf%par%domain))

            case("Antarctica")

                ! Modify bmb_obs in specific basins according to parameter values 
                do j = 1, size(mshlf%par%basin_name)
                    
                    select case(trim(mshlf%par%basin_name(j)))

                        case("ronne")
                            where(basins .ge.  1.0 .and. basins .le. 2.0) &
                                mshlf%now%bmb_obs = mshlf%par%basin_bmb(j)
                        case("queen")
                            where(basins .ge.  3.0 .and. basins .le.  5.0) & 
                                mshlf%now%bmb_obs = mshlf%par%basin_bmb(j)
                        case("amery")
                            where(basins .ge.  6.0 .and. basins .le. 7.0) & 
                                mshlf%now%bmb_obs = mshlf%par%basin_bmb(j)
                        case("wilkes")
                            where(basins .ge. 8.0 .and. basins .le. 10.0) & 
                                mshlf%now%bmb_obs = mshlf%par%basin_bmb(j)
                        case("ross")
                            where(basins .ge. 11.0 .and. basins .le. 12.0) & 
                                mshlf%now%bmb_obs = mshlf%par%basin_bmb(j)
                        case("pine")
                            where(basins .ge. 13.0 .and. basins .le. 15.0) & 
                                mshlf%now%bmb_obs = mshlf%par%basin_bmb(j)
                        case("peninsula")
                            where(basins .ge. 16.0 .and. basins .le. 19.0) & 
                                mshlf%now%bmb_obs = mshlf%par%basin_bmb(j)
                        case("west")
                            where((basins .ge. 12.0 .and. basins .le. 19.0) .or. (basins .eq. 1.0)) & 
                                mshlf%now%bmb_obs = mshlf%par%basin_bmb(j)
                        case("east")
                            where(basins .ge.  2.0 .and. basins .le. 11.0) & 
                                mshlf%now%bmb_obs = mshlf%par%basin_bmb(j)
                        
                        case DEFAULT 

                            ! Basin name not recognized, do nothing 

                    end select 

                end do 
                
            case("Greenland") 

                ! Modify specific basins according to parameter values 
                do j = 1, size(mshlf%par%basin_name)
                    
                    select case(trim(mshlf%par%basin_name(j)))

                        case("east")
                            where(basins .eq. 3.0) &
                                mshlf%now%bmb_obs = mshlf%par%basin_bmb(j)

                        case("northeast")
                            where(basins .eq. 2.0) &
                                mshlf%now%bmb_obs = mshlf%par%basin_bmb(j)

                        case("northwest")
                            where(basins .eq. 1.0 .or. basins .eq. 8.0) &
                                mshlf%now%bmb_obs = mshlf%par%basin_bmb(j)

                        case DEFAULT 

                            ! Basin name not recognized, do nothing 

                    end select 

                end do 
                
            case DEFAULT 

                ! Pass - no basins defined for other domains yet 

        end select 

        ! ==============================================
        ! Generate reference ocean mask 
        ! (0: land, 1: open ocean, 2: deep ocean) 

        ! Define mask_ocn_ref based on regions mask 
        ! (these definitions should work for all North and Antarctica domains)
        mshlf%now%mask_ocn_ref = 0 
        where (regions .eq. 1.0_wp) mshlf%now%mask_ocn_ref = 1   
        where (regions .eq. 2.0_wp) mshlf%now%mask_ocn_ref = 1 


        select case(trim(mshlf%par%domain))

            case("Greenland") 
                ! Greenland specific ocean kill regions

                where (regions .ne. 1.3) mshlf%now%mask_ocn_ref = 2 

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
                where (regions .ne. 2.11) mshlf%now%mask_ocn_ref = 2 

            ! case("North") 
            !     ! North specific ocean kill regions

            !     ! Apply only in purely open-ocean regions  
            !     where (regions .eq. 1.0) mshlf%now%mask_ocn_ref = 2 

            case DEFAULT 
                ! Other domains: c_deep potentially applied everywhere 
                ! with deep ocean points 

                mshlf%now%mask_ocn_ref = 2 

        end select 

        ! ==============================================


        ! ====================================
        !
        ! Summary
        !
        ! ====================================
        
        write(*,*) "range bmb_obs:      ", minval(mshlf%now%bmb_obs), maxval(mshlf%now%bmb_obs)
        write(*,*) "range mask_ocn_ref: ", minval(mshlf%now%mask_ocn_ref), maxval(mshlf%now%mask_ocn_ref)
        
        ! Initialize variables 
        mshlf%now%bmb_shlf      = 0.0  
        mshlf%now%bmb_ref       = 0.0 
        mshlf%now%kappa         = 0.0 
        mshlf%now%T_shlf        = 0.0
        mshlf%now%T_basin       = 0.0 
        mshlf%now%dT_shlf       = 0.0
        mshlf%now%dT_basin      = 0.0 

        return 

    end subroutine marshelf_init
    
    subroutine marshelf_calc_Tshlf(mshlf,H_ice,z_bed,f_grnd,basins,z_sl,dx,depth,to_ann,dto_ann)
        ! Calculate the 2D field of T_shlf (or dT_shlf) from 3D ocean temperature fields

        implicit none 

        type(marshelf_class), intent(INOUT) :: mshlf
        real(wp), intent(IN) :: H_ice(:,:) 
        real(wp), intent(IN) :: z_bed(:,:) 
        real(wp), intent(IN) :: f_grnd(:,:), basins(:,:) 
        real(wp), intent(IN) :: z_sl(:,:)
        real(wp), intent(IN) :: dx
        real(wp), intent(IN) :: depth(:),to_ann(:,:,:),dto_ann(:,:,:)
       
        ! 0. Calculate water temps at depths of interest ============================

        select case(trim(mshlf%par%T_shlf_method))

            case("mean") 
                ! Computes the temperature for a mean value
                call calc_shelf_temperature_mean(mshlf%now%T_shlf,mshlf%now%dT_shlf,depth,to_ann,dto_ann, &
                                                 depth_range=[mshlf%par%depth_min,mshlf%par%depth_max])
            
            case("layer")
                ! Takes the nearest layer for the temperature
                call calc_shelf_temperature_layer(mshlf%now%T_shlf,mshlf%now%dT_shlf,depth,to_ann,dto_ann,H_ice)

           case("interp")
               ! Interpolation from the two nearest layers 
               call  calc_shelf_temperature_depth(mshlf%now%T_shlf,mshlf%now%dT_shlf,depth,to_ann,dto_ann,H_ice)
 
            case DEFAULT
                write(*,*) "marshelf_update:: error: T_shlf_method not recognized: ", mshlf%par%T_shlf_method
                write(*,*) "Must be one of [mean,layer,interp]"
                stop
        
        end select
       
        ! jablasco: if method basin compute temperature by basin 
        select case(trim(mshlf%par%bmb_shlf_method))
            case("quad-nl")

                mshlf%now%T_basin = calc_shelf_basin(f_grnd,basins,H_ice,mshlf%now%T_shlf)
 
        end select

        return 

    end subroutine marshelf_calc_Tshlf

    subroutine marshelf_calc_Sshlf(mshlf,H_ice,depth,so)
        ! Calculate the 2D field of S_shlf from 3D ocean salinity
        ! fields
        ! TO DO: update S_shlf method as T_shlf
        ! Right now: layer for PICO -> simplest case

        implicit none

        type(marshelf_class), intent(INOUT) :: mshlf
        real(wp), intent(IN) :: H_ice(:,:)
        real(wp), intent(IN) :: depth(:),so(:,:,:)

        ! 0. Calculate water salinity at depths of interest
        ! ============================

        call calc_shelf_salinity_layer(mshlf%now%S_shlf,depth,so,H_ice)

        select case(trim(mshlf%par%T_shlf_method))

            case("mean")
                ! Computes the salinity for a mean value
                call calc_shelf_salinity_mean(mshlf%now%S_shlf,depth,so, &
                                                 depth_range=[mshlf%par%depth_min,mshlf%par%depth_max])

            case("layer")
                ! Takes the nearest layer for the salinity
                call calc_shelf_salinity_layer(mshlf%now%S_shlf,depth,so,H_ice)

           case("interp")
               ! Interpolation from the two nearest layers
               call  calc_shelf_salinity_depth(mshlf%now%S_shlf,depth,so,H_ice)

        end select

        return

    end subroutine marshelf_calc_Sshlf

    subroutine marshelf_set_Tshlf(mshlf,to_ann,dto_ann)
        ! Specify the 2D T_shlf/dT_shlf field according to a spatially constant value 

        implicit none 

        type(marshelf_class), intent(INOUT) :: mshlf
        real(wp),           intent(IN)    :: to_ann 
        real(wp),           intent(IN)    :: dto_ann 
        
        mshlf%now%T_shlf  = to_ann 
        mshlf%now%dT_shlf = dto_ann 

        return 

    end subroutine marshelf_set_Tshlf

    subroutine marshelf_update(mshlf,H_ice,z_bed,f_grnd,basins,z_sl,dx)
        
        implicit none
        
        type(marshelf_class), intent(INOUT) :: mshlf
        real(wp), intent(IN) :: H_ice(:,:) 
        real(wp), intent(IN) :: z_bed(:,:) 
        real(wp), intent(IN) :: f_grnd(:,:)
        real(wp), intent(IN) :: basins(:,:)
        real(wp), intent(IN) :: z_sl(:,:) 
        !real(wp), intent(IN) :: depth(:),to_ann(:,:,:),dto_ann(:,:,:)
        real(wp), intent(IN) :: dx   ! grid resolution [m]

        ! Local variables
        integer :: i, j, nx, ny, ngr 

        real(wp) :: bmb_floating, bmb_grline
        real(wp), allocatable :: H_ocn(:,:)
        logical,    allocatable :: is_grline(:,:)   
        real(wp) :: grz_wt

        nx = size(f_grnd,1)
        ny = size(f_grnd,2) 

        allocate(H_ocn(nx,ny)) 
        allocate(is_grline(nx,ny))

        ! Step 1: calculate geometry ===========================

        do j = 1, ny
        do i = 1, nx

            if (f_grnd(i,j) .eq. 0.0) then 
                ! Floating ice shelves
                
                ! Calculate height of ice-shelf base relative to sea level 
                mshlf%now%z_base(i,j) = z_sl(i,j) - (H_ice(i,j)*rho_ice_sw)
                
                ! Calculate ocean depth
                H_ocn(i,j) = max(mshlf%now%z_base(i,j) - z_bed(i,j),0.0) 

            else 
                ! Grounded ice, define for completeness
            
                mshlf%now%z_base(i,j)       = z_bed(i,j) 
                mshlf%now%slope_base(i,j)   = 0.0 
                H_ocn(i,j)                  = 0.0 

            end if 
            
        end do 
        end do 

        ! Calculate slope of ice-shelf base 
        call calc_shelf_slope(mshlf%now%slope_base,mshlf%now%z_base,f_grnd)
        
        ! Determine location of grounding line 
        is_grline = calc_grline(f_grnd)

        ! Determine resolution scaling at the grounding line 
        grz_wt = min( (mshlf%par%grz_length*1e3) / dx, 1.0) 
        

        if (mshlf%par%find_ocean) then 
            ! Determine which floating or partially floating points
            ! are connected to the open-ocean

            call find_open_ocean(mshlf%now%mask_ocn,f_grnd,mshlf%now%mask_ocn_ref)

        else 
            ! Set any floating or partially floating points to ocean points

            mshlf%now%mask_ocn = 0 
            where (f_grnd .lt. 1.0) mshlf%now%mask_ocn = 1
            where (f_grnd .lt. 1.0 .and. mshlf%now%mask_ocn_ref .eq. 2) mshlf%now%mask_ocn = 2

        end if

        ! 2. Define the reference bmb field for floating ice =================

        if (trim(mshlf%par%bmb_shlf_method) .eq. "anom") then
            ! Modify the basic bmb_obs fields by scalars

            mshlf%now%bmb_ref  = mshlf%now%bmb_obs

        else
            ! method == "abs", absolute method does not use reference melt term

            mshlf%now%bmb_ref = 0.0

            ! Redefine dT_shlf as the temperature anomaly wrt the freezing
            ! temperature
            mshlf%now%dT_shlf  = mshlf%now%T_shlf  - mshlf%par%T_fp
            mshlf%now%dT_basin = mshlf%now%T_basin - mshlf%par%T_fp

            ! Following Lipscomb et al (2021, Eq. 6 text), ensure dT_shlf >= 0 
            where (mshlf%now%dT_shlf  .lt. 0.0) mshlf%now%dT_shlf  = 0.0 
            where (mshlf%now%dT_basin .lt. 0.0) mshlf%now%dT_basin = 0.0 
            
        end if

        ! 3. Calculate current ice shelf bmb field (grounded-ice bmb is
        ! calculated in ice-sheet model separately) ========

        do j = 1, ny
        do i = 1, nx

            ! Calculate ocean depth
            H_ocn(i,j) = max( (z_sl(i,j)-z_bed(i,j))-(rho_ice_sw*H_ice(i,j)),0.0 )

            if (mshlf%now%mask_ocn(i,j) .gt. 0) then
                ! Floating ice points, including deep ocean (mask_ocn==2)

                ! Calculate the default floating bmb value
                if (trim(mshlf%par%bmb_shlf_method) .eq. "anom") then
                        bmb_floating = calc_bmb_shelf_anom(mshlf%now%bmb_ref(i,j),mshlf%now%dT_shlf(i,j), &
                                                        mshlf%par%c_shlf,mshlf%par%kappa_shlf)
                else if (trim(mshlf%par%bmb_shlf_method) .eq. "lin") then
                        bmb_floating = calc_bmb_shelf_linear(mshlf%now%dT_shlf(i,j),mshlf%par%gamma_lin)
                else if (trim(mshlf%par%bmb_shlf_method) .eq. "quad") then
                        bmb_floating = calc_bmb_shelf_quad(mshlf%now%dT_shlf(i,j),mshlf%par%gamma_quad)
                else if (trim(mshlf%par%bmb_shlf_method) .eq. "quad-nl") then
                        bmb_floating = calc_bmb_shelf_basin(mshlf%now%dT_shlf(i,j),mshlf%now%dT_basin(i,j),mshlf%par%gamma_quad_nl)
                end if

                if (is_grline(i,j)) then
                    ! Scale grounding line points by grounding line melt factor
                    ! accounting for the resolution dependence
                    if (trim(mshlf%par%bmb_shlf_method) .eq. "anom") then
                        bmb_grline = calc_bmb_shelf_anom(mshlf%now%bmb_ref(i,j),mshlf%now%dT_shlf(i,j), &
                                                         mshlf%par%c_grz,mshlf%par%kappa_grz)
                    else if (trim(mshlf%par%bmb_shlf_method) .eq. "lin") then
                        bmb_grline = calc_bmb_shelf_linear(mshlf%now%dT_shlf(i,j),mshlf%par%gamma_lin)
                        bmb_grline = mshlf%par%f_grz_shlf*bmb_grline
                    else if (trim(mshlf%par%bmb_shlf_method) .eq. "quad") then
                        bmb_grline = calc_bmb_shelf_quad(mshlf%now%dT_shlf(i,j),mshlf%par%gamma_quad)
                        bmb_grline = mshlf%par%f_grz_shlf*bmb_grline
                    else if (trim(mshlf%par%bmb_shlf_method) .eq. "quad-nl") then
                        bmb_grline = calc_bmb_shelf_basin(mshlf%now%dT_shlf(i,j),mshlf%now%dT_basin(i,j),mshlf%par%gamma_quad_nl)
                        bmb_grline = mshlf%par%f_grz_shlf*bmb_grline
                    end if

                    mshlf%now%bmb_shlf(i,j) = (1.0-grz_wt)*bmb_floating + grz_wt*bmb_grline

                else 
                    ! Only apply floating bmb value 

                    mshlf%now%bmb_shlf(i,j) = bmb_floating

                end if 
        

                ! Ensure that ice accretion is only 0% of melting
                ! Note: ajr: this should be a parameter!!
                ! jablasco: allow accretion 
                !if (mshlf%now%bmb_shlf(i,j) .gt. 0.0) mshlf%now%bmb_shlf(i,j) = mshlf%now%bmb_shlf(i,j)*0.0 

                ! Ensure that ice accretion only occurs where ice exists 
                if (mshlf%now%bmb_shlf(i,j) .gt. 0.0 .and. H_ice(i,j) .eq. 0.0) mshlf%now%bmb_shlf(i,j) = 0.0 

                ! Note: the following condition is a good idea, but in grisli-ucm the grounding line
                ! is defined as the last *grounded* point, so this limit sets bmb at grounding line to zero
                ! ! Ensure that refreezing rate is less than the available water depth within one time step 
                ! mshlf%now%bmb_shlf = max(mshlf%now%bmb_shlf,-H_ocn/dt)
                
            else 
                ! Grounded point, or floating point not connected to the ocean 
                ! Set bmb_shlf to zero 

                mshlf%now%bmb_shlf(i,j) = 0.0
                
            end if 

        end do 
        end do  

        ! Apply c_deep melt value to deep ocean points
        ! n_smth can be zero when c_deep is a relatively small value (like -1 m/a). For
        ! c_deep = -50 m/a, n_smth=2 is more appropriate to ensure a smooth transition to 
        ! high melt rates. 
        call apply_c_deep(mshlf%par,mshlf%now%bmb_shlf,mshlf%now%mask_ocn,z_bed,z_sl,n_smth=0) 

        return
        
    end subroutine marshelf_update


    subroutine marshelf_end(mshlf)

        implicit none 

        type(marshelf_class) :: mshlf 

        ! Deallocate mshlf state object
        call marshelf_deallocate(mshlf%now)

        return 

    end subroutine marshelf_end


    subroutine marshelf_par_load(par,filename,domain,grid_name,init)

        type(marshelf_param_class), intent(OUT) :: par 
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: domain, grid_name  
        logical, optional :: init 
        logical :: init_pars 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 

        call nml_read(filename,"marine_shelf","bmb_shlf_method",par%bmb_shlf_method,init=init_pars)
        call nml_read(filename,"marine_shelf","T_shlf_method",  par%T_shlf_method,  init=init_pars)
        call nml_read(filename,"marine_shelf","use_obs",        par%use_obs,        init=init_pars)
        call nml_read(filename,"marine_shelf","obs_path",       par%obs_path,       init=init_pars)
        call nml_read(filename,"marine_shelf","obs_name",       par%obs_name,       init=init_pars)
        call nml_read(filename,"marine_shelf","obs_f",          par%obs_f,          init=init_pars)
        call nml_read(filename,"marine_shelf","obs_lim",        par%obs_lim,        init=init_pars)
        call nml_read(filename,"marine_shelf","basin_name",     par%basin_name,     init=init_pars)
        call nml_read(filename,"marine_shelf","basin_bmb",      par%basin_bmb,      init=init_pars)
        call nml_read(filename,"marine_shelf","c_shlf",         par%c_shlf,         init=init_pars)
        call nml_read(filename,"marine_shelf","kappa_shlf",     par%kappa_shlf,     init=init_pars)
        call nml_read(filename,"marine_shelf","f_grz_shlf",     par%f_grz_shlf,     init=init_pars)
        call nml_read(filename,"marine_shelf","c_grz",          par%c_grz,          init=init_pars)
        call nml_read(filename,"marine_shelf","kappa_grz",      par%kappa_grz,      init=init_pars)
        call nml_read(filename,"marine_shelf","gamma_lin",      par%gamma_lin,      init=init_pars)
        call nml_read(filename,"marine_shelf","gamma_quad",     par%gamma_quad,     init=init_pars)
        call nml_read(filename,"marine_shelf","gamma_quad_nl",  par%gamma_quad_nl,  init=init_pars)
        call nml_read(filename,"marine_shelf","grz_length",     par%grz_length,     init=init_pars)
        call nml_read(filename,"marine_shelf","c_deep",         par%c_deep,         init=init_pars)
        call nml_read(filename,"marine_shelf","depth_deep",     par%depth_deep,     init=init_pars)
        call nml_read(filename,"marine_shelf","T_fp",           par%T_fp,           init=init_pars)
        call nml_read(filename,"marine_shelf","depth_min",      par%depth_min,      init=init_pars)
        call nml_read(filename,"marine_shelf","depth_max",      par%depth_max,      init=init_pars)    
        call nml_read(filename,"marine_shelf","find_ocean",     par%find_ocean,     init=init_pars)   
        
        ! Determine derived parameters
        call parse_path(par%obs_path,domain,grid_name)
        par%domain = trim(domain)

        if (par%f_grz_shlf .eq. 0.0) then 
            write(*,*) "marshelf_par_load:: Error: f_grz_shlf cannot be zero."
            stop 
        end if 

        ! Define some parameters internally based on choices 
        par%kappa_shlf = par%kappa_grz * (1.0/par%f_grz_shlf)
        par%c_shlf     = par%c_grz     * (1.0/par%f_grz_shlf)
        
        ! jablasco: convert T_fp from C to K
        par%T_fp       = par%T_fp + 273.15

        return

    end subroutine marshelf_par_load

    ! =======================================================
    !
    ! marine_shelf interface physics
    !
    ! =======================================================
    
    elemental function calc_bmb_shelf_anom(bmb_ref,dT,c,kappa) result(bmb)
        ! Calculate basal mass balance of shelf (floating) ice [m/a] 
        ! as a function of a reference state bmb_ref and 
        ! a heat flux. 
        ! The variable dT can either represent the shelf temperature 
        ! relative to the water freezing point, or the temporal anomaly
        ! relative to a reference state

        implicit none

        real(wp), intent(IN)    :: bmb_ref, dT, c, kappa
        real(wp) :: bmb
        
        ! Use temperature anomaly relative to a reference state
        bmb = c*bmb_ref - kappa*dT
        
        return

    end function calc_bmb_shelf_anom

    elemental function calc_bmb_shelf_linear(dT_shlf,gamma_lin) result(bmb)
        ! Calculate basal mass balance of shelf (floating) ice [m/a] 
        ! as a linear law following Beckmann and Goosse (2003).
       ! The linear, local dependencyon thermal forcing assumes a balance 
        ! between vertical diffusive heat flux across the oceancavity top 
        ! boundary layer and latent heat due to melting and freezing. 

        implicit none

        real(wp), intent(IN)    :: dT_shlf, gamma_lin
        real(wp) :: bmb

        ! jablasco: careful!!
        ! -1.0 because bmb is negative
        bmb = -1.0*gamma_lin*(dT_shlf*(rho_sw*cp_o)/(rho_ice*L_ice))

        return

    end function calc_bmb_shelf_linear

    elemental function calc_bmb_shelf_quad(dT_shlf,gamma_quad) result(bmb)
        ! Calculate basal mass balance of shelf (floating) ice [m/a] 
        ! as a quadratic law following (Holland et al., 2008) with local ocean T.
        ! This accounts for positive feedback between the sub-shelf melting and the circulation in the cavity.

        implicit none

        real(wp), intent(IN)    :: dT_shlf, gamma_quad
        real(wp) :: bmb

        ! jablasco: careful!!
        ! -1.0 because bmb is negative
        bmb = -1.0*gamma_quad*(dT_shlf*(rho_sw*cp_o)/(rho_ice*L_ice))**2

        return

    end function calc_bmb_shelf_quad

        ! jablasco: quadratic function with mean shelf temperature
    elemental function calc_bmb_shelf_basin(dT_shlf,dT_basin,gamma_quad_nl) result(bmb)
        ! Calculate basal mass balance of shelf (floating) ice [m/a] 
        ! as a quadratic law following (Holland et al., 2008) but with local ocean T and mean shelf T.
        ! This accounts for positive feedback between the sub-shelf melting and the circulation in the cavity.

        implicit none

        real(wp), intent(IN)    :: dT_shlf, dT_basin, gamma_quad_nl
        real(wp) :: bmb

        ! jablasco: careful!!
        ! -1.0 because bmb is negative
        bmb = -1.0*gamma_quad_nl*(dT_shlf*dT_basin)*((rho_sw*cp_o)/(rho_ice*L_ice))**2

        return

    end function calc_bmb_shelf_basin

    subroutine apply_c_deep(par,bmb,mask_ocn,z_bed,z_sl,n_smth)
        ! Apply c_deep for killing ice in the deep ocean 
        
        implicit none 

        type(marshelf_param_class) :: par 
        real(wp), intent(INOUT)  :: bmb(:,:)  
        integer,    intent(IN)     :: mask_ocn(:,:) 
        real(wp), intent(IN)     :: z_bed(:,:) 
        real(wp), intent(IN)     :: z_sl(:,:)
        integer,    intent(IN)     :: n_smth       ! Smoothing neighborhood radius in grid points

        ! Local variables 
        integer :: i, j, nx, ny 
        real(wp) :: n_pts 
        logical,    allocatable :: is_c_deep(:,:) 
        real(wp), allocatable :: bmb_tmp(:,:) 

        nx = size(bmb,1)
        ny = size(bmb,2) 

        allocate(is_c_deep(nx,ny))
        allocate(bmb_tmp(nx,ny))

        if (n_smth .lt. 0 .or. n_smth .gt. 2) then 
            write(*,*) "apply_c_deep:: Error: n_smth should be 0, 1, or 2."
            write(*,*) "n_smth = ", n_smth 
            stop 
        end if 

        ! Apply c_deep to appropriate regions or keep bmb, whichever is more negative

        where(mask_ocn .eq. 2 .and. z_sl-z_bed .ge. par%depth_deep) bmb = par%c_deep 
        ! jablasco: parche
            !bmb = min(bmb,par%c_deep)

            !is_c_deep = .TRUE. 

        !elsewhere 

            !is_c_deep = .FALSE. 

        !end where 


        ! Make sure bmb transitions smoothly to c_deep value from neighbors 
        
        if (n_smth .gt. 0) then 
            
            ! Get size of neighborhood
            n_pts = real( (2*n_smth+1)**2, wp) 

            bmb_tmp = bmb 
            do i = 1+n_smth, nx-n_smth
            do j = 1+n_smth, ny-n_smth 

                ! Apply a neighborhood average value, if deep ocean is in neighborhood of shallow ocean
                if (is_c_deep(i,j) .and. &
                    count(.not. is_c_deep(i-n_smth:i+n_smth,j-n_smth:j+n_smth)) .gt. 0) then

                    bmb(i,j) = sum(bmb_tmp(i-n_smth:i+n_smth,j-n_smth:j+n_smth))/n_pts

                end if 

            end do 
            end do 

        end if 

        return 

    end subroutine apply_c_deep

    subroutine calc_shelf_slope(slope,z_base,f_grnd)
        ! Calculate the slope of the ice shelf base centered on aa-node, 
        ! where slope = sin(theta) = height/hypotenuse

        implicit none 

        real(wp), intent(OUT) :: slope(:,:) 
        real(wp), intent(IN)  :: z_base(:,:) 
        real(wp), intent(IN)  :: f_grnd(:,:) 

        ! Local variables 
        integer :: i, j, ip1, im1, jp1, jm1 

        
        return 

    end subroutine calc_shelf_slope

    function calc_grline(f_grnd) result(is_grline)
        ! Determine the grounding line given the grounded fraction f_grnd
        ! ie, is_grline is true for a floating point or partially floating  
        ! point with grounded neighbors

        implicit none 

        real(wp), intent(IN) :: f_grnd(:,:)
        logical :: is_grline(size(f_grnd,1),size(f_grnd,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(is_grline,1)
        ny = size(is_grline,2)

        is_grline = .FALSE. 
        do i = 2, nx-1 
        do j = 2, ny-1 
            
            ! Floating point or partially floating point with grounded neighbors
            if (f_grnd(i,j) .lt. 1.0 .and. &
                (f_grnd(i-1,j) .eq. 1.0 .or. f_grnd(i+1,j) .eq. 1.0 .or. &
                 f_grnd(i,j-1) .eq. 1.0 .or. f_grnd(i,j+1) .eq. 1.0) ) then 
                is_grline(i,j) = .TRUE. 

            end if 
            
        end do 
        end do 

        return 

    end function calc_grline

    function calc_shelf_basin(f_grnd,basins,H_ice,T_shlf) result(T_basin)
        ! jablasco: function to compute mean ocean temperature of ice-shelf basin
        implicit none

        real(wp), intent(IN)  :: f_grnd(:,:), basins(:,:), H_ice(:,:), T_shlf(:,:)
        real(wp) :: T_basin(size(f_grnd,1),size(f_grnd,2))

        ! Local variables
        integer :: i, j, m
        real(wp) :: n_pts, t_mean

        ! Assign to shelf values real ocean values
        T_basin(:,:) = T_shlf(:,:) 

        do m=1, int(maxval(basins))
            n_pts  = 0.0
            t_mean = 0.0
            do i = 1, size(T_shlf,1)
                do j = 1, size(T_shlf,2)
                    ! Floating ice point
                    if (f_grnd(i,j) .lt. 1.0 .and. H_ice(i,j) .gt. 0.0 .and. basins(i,j) .eq. m) then
                        !t_mean = t_mean+(1.0-f_grnd(i,j))*T_shlf(i,j)
                        t_mean = t_mean + T_shlf(i,j)
                        !n_pts = n_pts + (1.0-f_grnd(i,j))
                        n_pts = n_pts + 1.0
                    end if
                end do
            end do

            ! Assign shelf value
            do i = 1, size(T_shlf,1)
                do j = 1, size(T_shlf,2)
                    ! Basin floating ice point
                    if (f_grnd(i,j) .lt. 1.0 .and. H_ice(i,j) .gt. 0.0 .and. basins(i,j) .eq. m) then
                        T_basin(i,j) = t_mean / n_pts
                    end if
                end do
            end do

        end do

        return

    end function calc_shelf_basin

    subroutine calc_shelf_temperature_mean(Tshlf,dTshlf,depth,T_ocn,dT_ocn,depth_range)
        ! Calculate water temp for a given range of depths
        
        implicit none 
        
        real(wp), intent(OUT)    :: Tshlf(:,:), dTshlf(:,:)
        real(wp), intent(IN)     :: depth(:), T_ocn(:,:,:), dT_ocn(:,:,:)
        real(wp), intent(IN)     :: depth_range(:)

        ! Local variables 
        integer :: k0, k1 

        ! Note: this requires that k1 > k0, and it should be
        ! weighted by the thickness of the layers (to do!)
        ! Note: depth is z-coordinate, ie positive below sea level  
        k0 = minloc(abs(depth-depth_range(1)),dim=1)
        k1 = minloc(abs(depth-depth_range(2)),dim=1)

        if (k1 < k0) then 
            write(*,*) "calc_shelf_temperature_mean:: error in depth_range calculation."
            write(*,*) "depth_min, depth_max: ", depth_range 
            write(*,*) "indices(k0,k1): ", k0, k1
            write(*,*) "depths(k0,k1):  ", depth(k0), depth(k1) 
            stop 
        end if 

        ! Get the mean water temperature for these depths
        Tshlf  = sum(T_ocn(:,:,k0:k1),dim=3)  / (k1-k0+1)
        dTshlf = sum(dT_ocn(:,:,k0:k1),dim=3) / (k1-k0+1)

        return 

    end subroutine calc_shelf_temperature_mean

    subroutine calc_shelf_temperature_layer(Tshlf,dTshlf,depth,T_ocn,dT_ocn,H)
        ! Calculates the water temperature at the depth of the ice shelf
        ! It assigns the temperature of the nearest layer
        
        implicit none

        real(wp), intent(OUT)   :: Tshlf(:,:), dTshlf(:,:)
        real(wp), intent(IN)    :: depth(:), T_ocn(:,:,:), dT_ocn(:,:,:), H(:,:)

        ! Local variables
        integer :: k0, i, j

        do i = 1, size(Tshlf,1)
            do j = 1, size(Tshlf,2)
                ! The depth of the ice shelve is a 90% of the ice thickness
                k0 = minloc(abs(depth-rho_ice_sw*H(i,j)),dim=1)
                Tshlf(i,j)  = T_ocn(i,j,k0)
                dTshlf(i,j) = dT_ocn(i,j,k0)
            end do
        end do 

        return

    end subroutine calc_shelf_temperature_layer

    subroutine calc_shelf_temperature_depth(Tshlf,dTshlf,depth,T_ocn,dT_ocn,H)
        ! Calculates the water temperature from linear interpolation of vertical profile
        
        implicit none

        real(wp), intent(OUT)   :: Tshlf(:,:), dTshlf(:,:)
        real(wp), intent(IN)    :: depth(:), T_ocn(:,:,:), dT_ocn(:,:,:), H(:,:)

        ! Local variables
        integer :: i, j
        real(wp) :: depth_H 

        do i = 1, size(Tshlf,1)
            do j = 1, size(Tshlf,2)

                ! Depth of ice shelf
                depth_H = rho_ice_sw*H(i,j) 
                    
                ! Linearly interpolate to depth of ice shelf 
                Tshlf(i,j)  = interp_linear(depth,T_ocn(i,j,:), xout=depth_H)
                dTshlf(i,j) = interp_linear(depth,dT_ocn(i,j,:),xout=depth_H)

            end do
        end do

        return

    end subroutine calc_shelf_temperature_depth

    subroutine calc_shelf_salinity_layer(Sshlf,depth,S_ocn,H)
        ! Calculates the water salinity at the depth of the ice shelf
        ! It assigns the salinity of the nearest layer

        implicit none

        real(wp), intent(OUT)   :: Sshlf(:,:)
        real(wp), intent(IN)    :: depth(:), S_ocn(:,:,:), H(:,:)

        ! Local variables
        integer :: k0, i, j

        do i = 1, size(Sshlf,1)
            do j = 1, size(Sshlf,2)
                ! The depth of the ice shelve is a 90% of the ice thickness
                k0 = minloc(abs(depth-rho_ice_sw*H(i,j)),dim=1)
                Sshlf(i,j)  = S_ocn(i,j,k0)
            end do
        end do

        return

    end subroutine calc_shelf_salinity_layer

    subroutine calc_shelf_salinity_depth(Sshlf,depth,S_ocn,H)

        implicit none

        real(wp), intent(OUT)   :: Sshlf(:,:)
        real(wp), intent(IN)    :: depth(:), S_ocn(:,:,:), H(:,:)

        ! Local variables
        integer :: i, j
        real(wp) :: depth_H

        do i = 1, size(Sshlf,1)
            do j = 1, size(Sshlf,2)

                ! Depth of ice shelf
                depth_H = rho_ice_sw*H(i,j)

                ! Linearly interpolate to depth of ice shelf
                Sshlf(i,j)  = interp_linear(depth,S_ocn(i,j,:), xout=depth_H)

            end do
        end do

        return

    end subroutine calc_shelf_salinity_depth

    subroutine calc_shelf_salinity_mean(Sshlf,depth,S_ocn,depth_range)
        ! Calculate water temp for a given range of depths

        implicit none

        real(wp), intent(OUT)    :: Sshlf(:,:)
        real(wp), intent(IN)     :: depth(:), S_ocn(:,:,:)
        real(wp), intent(IN)     :: depth_range(:)

        ! Local variables
        integer :: k0, k1

        ! Note: this requires that k1 > k0, and it should be
        ! weighted by the thickness of the layers (to do!)
        ! Note: depth is z-coordinate, ie positive below sea level
        k0 = minloc(abs(depth-depth_range(1)),dim=1)
        k1 = minloc(abs(depth-depth_range(2)),dim=1)

        if (k1 < k0) then
            write(*,*) "calc_shelf_salinity_mean:: error in depth_range calculation."
            write(*,*) "depth_min, depth_max: ", depth_range
            write(*,*) "indices(k0,k1): ", k0, k1
            write(*,*) "depths(k0,k1):  ", depth(k0), depth(k1)
            stop
        end if

        ! Get the mean water temperature for these depths
        Sshlf  = sum(S_ocn(:,:,k0:k1),dim=3)  / (k1-k0+1)

        return

    end subroutine calc_shelf_salinity_mean

    ! =======================================================
    !
    ! marshelf memory management
    !
    ! =======================================================

    subroutine marshelf_allocate(now,nx,ny)

        implicit none 

        type(marshelf_state_class) :: now 
        integer :: nx, ny 

        ! Make object is deallocated
        call marshelf_deallocate(now)

        ! Allocate marshelf 
        allocate(now%bmb_shlf(nx,ny))
        allocate(now%bmb_obs(nx,ny))
        allocate(now%bmb_ref(nx,ny))
        allocate(now%T_shlf(nx,ny))
        allocate(now%T_basin(nx,ny))
        allocate(now%dT_shlf(nx,ny))
        allocate(now%dT_basin(nx,ny))
        allocate(now%S_shlf(nx,ny))
        allocate(now%kappa(nx,ny))
        
        allocate(now%z_base(nx,ny))
        allocate(now%slope_base(nx,ny))

        allocate(now%mask_ocn_ref(nx,ny))
        allocate(now%mask_ocn(nx,ny))

        ! Initialize variables 
        now%bmb_shlf        = 0.0  
        now%bmb_obs         = 0.0
        now%bmb_ref         = 0.0  
        now%T_shlf          = 0.0
        now%T_basin         = 0.0
        now%dT_shlf         = 0.0
        now%dT_basin        = 0.0
        now%S_shlf          = 0.0
        now%kappa           = 0.0
        
        now%z_base          = 0.0 
        now%slope_base      = 0.0 
        
        ! By default set ocean points everywhere
        now%mask_ocn_ref    = 1
        now%mask_ocn        = 1

        return

    end subroutine marshelf_allocate

    subroutine marshelf_deallocate(now)

        implicit none 

        type(marshelf_state_class) :: now 

        ! Allocate state objects
        if (allocated(now%bmb_shlf))    deallocate(now%bmb_shlf)
        if (allocated(now%bmb_obs))     deallocate(now%bmb_obs)
        if (allocated(now%bmb_ref))     deallocate(now%bmb_ref)
        if (allocated(now%T_shlf))      deallocate(now%T_shlf)
        if (allocated(now%T_basin))     deallocate(now%T_basin)
        if (allocated(now%dT_shlf))     deallocate(now%dT_shlf)
        if (allocated(now%dT_basin))    deallocate(now%dT_basin)
        if (allocated(now%kappa))       deallocate(now%kappa)
        
        if (allocated(now%z_base))      deallocate(now%z_base)
        if (allocated(now%slope_base))  deallocate(now%slope_base)
        
        if (allocated(now%mask_ocn_ref))    deallocate(now%mask_ocn_ref)
        if (allocated(now%mask_ocn))        deallocate(now%mask_ocn)
        
        return

    end subroutine marshelf_deallocate

    subroutine parse_path(path,domain,grid_name)

        implicit none 

        character(len=*), intent(INOUT) :: path 
        character(len=*), intent(IN)    :: domain, grid_name 

        call nml_replace(path,"{domain}",   trim(domain))
        call nml_replace(path,"{grid_name}",trim(grid_name))
        
        return 

    end subroutine parse_path
    
    function interp_linear(x,y,xout) result(yout)
        ! Simple linear interpolation of a point

        implicit none 
 
        real(wp), dimension(:), intent(IN) :: x, y
        real(wp), intent(IN) :: xout
        real(wp) :: yout 
        integer :: i, j, n, nout 
        real(wp) :: alph

        n    = size(x) 

        if (xout .lt. x(1)) then
            yout = y(1)
        else if (xout .gt. x(n)) then
            yout = y(n)
        else
            do j = 1, n 
                if (x(j) .ge. xout) exit 
            end do

            if (j .eq. 1) then 
                yout = y(1) 
            else if (j .eq. n+1) then 
                yout = y(n)
            else 
                alph = (xout - x(j-1)) / (x(j) - x(j-1))
                yout = y(j-1) + alph*(y(j) - y(j-1))
            end if 
        end if 

        return 

    end function interp_linear
    
    subroutine find_open_ocean(mask,f_grnd,mask_ref)
        ! Brute-force routine to find all ocean points 
        ! (ie, when f_grnd < 1) connected to
        ! the open ocean as defined in mask_ref. 

        implicit none 

        integer,    intent(OUT) :: mask(:,:) 
        real(wp), intent(IN)  :: f_grnd(:,:)
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

end module marine_shelf


