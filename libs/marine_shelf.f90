

module marine_shelf
    ! Module to simulate the marine-shelf interface:
    ! Calculates the basal mass balance of an ice shelf (bmb_shlf)

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    use nml 
    use ncio 

    use pico 

!     use yelmo_defs, only : sp, dp, wp, rho_ice, rho_w, rho_sw, g, parse_path 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the working precision of the library (sp,dp)
    integer,  parameter :: wp = sp 

    ! Physical constants 
    real(wp), parameter :: rho_ice =  918.d0    ! Density ice           [kg/m^3] 
    real(wp), parameter :: rho_w   = 1000.d0    ! Density water         [kg/m^3] 
    real(wp), parameter :: rho_sw  = 1028.d0    ! Density seawater      [kg/m^3] 
    real(wp), parameter :: g       = 9.81d0     ! Gravitational accel.  [m/s^2]
    real(wp), parameter :: cp_o    = 3974.d0    ! [J/(kg K)] Specific heat capacity of ocean mixed layer 
    real(wp), parameter :: L_ice   = 3.34e5     ! [J/kg] Latent heat of fusion of ice 
    
    real(wp), parameter :: T0      = 273.15     ! [K] Reference freezing temp 

    real(wp), parameter :: lambda     = L_ice/cp_o
    real(wp), parameter :: rho_ice_sw = rho_ice / rho_sw 
    real(wp), parameter :: omega      = (rho_sw*cp_o) / (rho_ice*L_ice)     ! [1/K]

    ! Global mask value definitions
    integer, parameter :: mask_val_land           = 0 
    integer, parameter :: mask_val_grounding_line = 1
    integer, parameter :: mask_val_floating_line  = 2
    integer, parameter :: mask_val_ocean          = 3
    integer, parameter :: mask_val_deep_ocean     = 4
    integer, parameter :: mask_val_lake           = 5
        
    type marshelf_param_class

        character(len=56)   :: bmb_method
        integer             :: tf_method 
        character(len=56)   :: interp_method 
        character(len=56)   :: interp_depth 
        logical             :: find_ocean 
        logical             :: use_obs
        character(len=512)  :: obs_path
        character(len=56)   :: obs_name 
        real(wp)            :: obs_scale, obs_lim
        character(len=56)   :: corr_method 
        real(wp)            :: basin_number(50)
        real(wp)            :: basin_bmb_corr(50)
        real(wp)            :: basin_tf_corr(50)
        logical             :: tf_correction
        character(len=512)  :: tf_path
        character(len=56)   :: tf_name

        real(wp)            :: bmb_max  
        real(wp)            :: c_deep
        real(wp)            :: depth_deep
        real(wp)            :: depth_const
        real(wp)            :: depth_min
        real(wp)            :: depth_max

        real(wp)            :: lambda1, lambda2, lambda3 
        real(wp)            :: gamma_lin, gamma_quad, gamma_quad_nl 
        real(wp)            :: gamma_prime 

        real(wp)            :: c_grz, kappa_grz, f_grz_shlf, grz_length
        
        character(len=512) :: domain   
        real(wp) :: rho_ice, rho_sw

    end type 

    type marshelf_state_class 
        real(wp), allocatable :: bmb_shlf(:,:)          ! Shelf basal mass balance [m/a]
        real(wp), allocatable :: bmb_ref(:,:)           ! Basal mass balance reference field
        real(wp), allocatable :: bmb_corr(:,:)          ! Basal mass balance correction field
        
        real(wp), allocatable :: T_shlf(:,:)            ! [K] Shelf temperature
        real(wp), allocatable :: dT_shlf(:,:)           ! [K] Shelf temperature anomaly relative to ref. state
        real(wp), allocatable :: S_shlf(:,:)            ! [K] Shelf temperature
        real(wp), allocatable :: T_fp_shlf(:,:)         ! [K] Shelf freezing-point temperature
        
        real(wp), allocatable :: tf_shlf(:,:)           ! Thermal forcing at the ice-shelf interface
        real(wp), allocatable :: tf_basin(:,:)          ! Basin-average thermal forcing at the ice-shelf interface
        real(wp), allocatable :: tf_corr(:,:)           ! Thermal correction at the ice-shelf interface by data (ismip6)
        real(wp), allocatable :: tf_corr_basin(:,:)     ! Thermal correction at the ice-shelf interface by basin
        
        real(wp), allocatable :: z_base(:,:)            ! Ice-shelf base elevation (relative to sea level)
        real(wp), allocatable :: slope_base(:,:)        ! Ice-shelf base slope (slope=sin(theta)=length/hypotenuse)

        integer,  allocatable :: mask_ocn_ref(:,:) 
        integer,  allocatable :: mask_ocn(:,:) 
        
    end type 

    type marshelf_class
        type(marshelf_param_class) :: par 
        type(marshelf_state_class) :: now
        type(pico_class)           :: pico  
    end type

    private
    public :: marshelf_class
    public :: marshelf_update_shelf
    public :: marshelf_update
    public :: marshelf_init
    public :: marshelf_end 

contains 
    
    subroutine marshelf_update_shelf(mshlf,H_ice,z_bed,f_grnd,basins,z_sl,dx, &
                                    depth,to_ann,so_ann,dto_ann,tf_ann)
        ! Calculate various 2D fields from 3D ocean fields representative 
        ! for the ice-shelf interface: T_shlf, dT_shlf, S_shlf 

        implicit none 

        type(marshelf_class), intent(INOUT) :: mshlf
        real(wp), intent(IN) :: H_ice(:,:) 
        real(wp), intent(IN) :: z_bed(:,:) 
        real(wp), intent(IN) :: f_grnd(:,:)
        real(wp), intent(IN) :: basins(:,:) 
        real(wp), intent(IN) :: z_sl(:,:)
        real(wp), intent(IN) :: dx
        real(wp), intent(IN) :: depth(:)
        real(wp), intent(IN) :: to_ann(:,:,:)
        real(wp), intent(IN) :: so_ann(:,:,:)
        real(wp), intent(IN) :: dto_ann(:,:,:)
        real(wp), intent(IN), optional :: tf_ann(:,:,:)

        ! Local variables
        integer :: i, j, nx, ny, nz 
        real(wp), allocatable :: depth_shlf 
        real(wp), allocatable :: wt_shlf(:) 

        nx = size(H_ice,1) 
        ny = size(H_ice,2) 
        nz = size(depth,1) 

        allocate(wt_shlf(nz)) 
        wt_shlf = 0.0 

        ! Loop over domain and update variables at each point (vertical interpolation)
        do j = 1, ny 
        do i = 1, nx 

            ! 1. Calculate the depth of the current shelf base 

            select case(trim(mshlf%par%interp_depth))

                case("shlf")
                    ! Assign the depth to the shelf depth 

                    if(H_ice(i,j) .gt. 0.0 .and. f_grnd(i,j) .lt. 1.0) then
                        ! Floating ice, depth == z_ice_base
                        depth_shlf = H_ice(i,j)*rho_ice_sw
                    else if(H_ice(i,j) .gt. 0.0 .and. f_grnd(i,j) .eq. 1.0) then
                        ! Grounded ice, depth == H_ocn = z_sl-z_bed 
                        depth_shlf = z_sl(i,j) - z_bed(i,j)
                    else 
                        ! Open ocean, depth == constant value, eg 2000 m. 
                        depth_shlf = mshlf%par%depth_const
                    end if

                case("bed")
                    ! Assign the depth corresponding to the bedrock 

                    depth_shlf = z_sl(i,j) - z_bed(i,j)

                case("const")

                    depth_shlf = mshlf%par%depth_const 

                case DEFAULT

                    write(*,*) "marshelf_update_shelf:: Error: interp_depth method not recognized."
                    write(*,*) "interp_depth = ", trim(mshlf%par%interp_depth)

            end select
            
            ! 2. Calculate weighting function for vertical depths ===========================

            select case(trim(mshlf%par%interp_method))

                case("mean") 
                    ! Equal weighting of layers within a specified depth range

                    call calc_shelf_variable_mean(wt_shlf,depth, &
                                depth_range=[mshlf%par%depth_min,mshlf%par%depth_max])

                case("layer")
                    ! All weight given to the nearest layer to depth of shelf

                    call calc_shelf_variable_layer(wt_shlf,depth,depth_shlf)
                    
                case("interp")
                    ! Interpolation weights from the two nearest layers to depth of shelf

                    call  calc_shelf_variable_depth(wt_shlf,depth,depth_shlf)

                case DEFAULT
                    write(*,*) "marshelf_update:: error: interp_method not recognized: ", mshlf%par%interp_method
                    write(*,*) "Must be one of [mean,layer,interp]"
                    stop
                
            end select
            
            ! Normalize weighting function 
            if (sum(wt_shlf) .gt. 0.0_wp) then 
                wt_shlf = wt_shlf / sum(wt_shlf) 
            else 
                write(*,*) "marshelf_update_shelf:: Error: weighting should be > 0."
                stop 
            end if 

            ! 3. Calculate water properties at depths of interest ============================

            mshlf%now%T_shlf(i,j)  = sum(to_ann(i,j,:) *wt_shlf)
            mshlf%now%S_shlf(i,j)  = sum(so_ann(i,j,:) *wt_shlf)
            mshlf%now%dT_shlf(i,j) = sum(dto_ann(i,j,:)*wt_shlf)

            if (present(tf_ann)) then 
                mshlf%now%tf_shlf(i,j) = sum(tf_ann(i,j,:)*wt_shlf)
            end if 

        end do 
        end do  

        return 

    end subroutine marshelf_update_shelf

    subroutine marshelf_update(mshlf,H_ice,z_bed,f_grnd,regions,basins,z_sl,dx)
        
        implicit none
       
        ! jablasco: added regions 
        type(marshelf_class), intent(INOUT) :: mshlf
        real(wp), intent(IN) :: H_ice(:,:) 
        real(wp), intent(IN) :: z_bed(:,:) 
        real(wp), intent(IN) :: f_grnd(:,:)
        real(wp), intent(IN) :: regions(:,:)
        real(wp), intent(IN) :: basins(:,:)
        real(wp), intent(IN) :: z_sl(:,:) 
        !real(wp), intent(IN) :: depth(:),to_ann(:,:,:),dto_ann(:,:,:)
        real(wp), intent(IN) :: dx   ! grid resolution [m]

        ! Local variables
        integer :: i, j, nx, ny, ngr 

        logical,  allocatable :: is_grline(:,:)   
        real(wp), allocatable :: gamma2D(:,:)
        logical :: with_slope 

        nx = size(f_grnd,1)
        ny = size(f_grnd,2) 

        allocate(is_grline(nx,ny))
        allocate(gamma2D(nx,ny))

        ! Step 1: calculate geometry ===========================

        do j = 1, ny
        do i = 1, nx

            if (f_grnd(i,j) .lt. 1.0) then 
                ! Floating ice shelves
                ! (even assume floating at grounding line to allow
                ! reasonable calculations for non-binary f_grnd)
                
                ! Calculate height of ice-shelf base relative to sea level 
                mshlf%now%z_base(i,j) = z_sl(i,j) - (H_ice(i,j)*rho_ice_sw)
                
            else 
                ! Grounded ice, define for completeness

                mshlf%now%z_base(i,j) = z_bed(i,j)

            end if 
            
        end do 
        end do 

        ! Calculate slope of ice-shelf base 
        call calc_shelf_slope(mshlf%now%slope_base,mshlf%now%z_base,dx)
        
        ! Define the ocean mask (land, grounding line, floating line, ocean, deep ocean, lakes)
        call set_ocean_mask(mshlf%now%mask_ocn,f_grnd,mshlf%now%mask_ocn_ref,mshlf%par%find_ocean)


        ! 2. Calculate various fields of interest ============================

        ! Calculate ocean water freezing point [K]
        call calc_freezing_point(mshlf%now%T_fp_shlf,mshlf%now%S_shlf,mshlf%now%z_base, &
                                    mshlf%par%lambda1,mshlf%par%lambda2,mshlf%par%lambda3, &
                                    T_ref=T0)

        ! Calculate the thermal forcing, if desired 
        ! (not used for pico, but good to diagnose anyway and needed for other methods)
        
        select case(mshlf%par%tf_method)

            case(0)
                ! tf_shlf defined externally, do nothing 

            case(1) 
                ! Determine tf_shlf via T_shlf - T_fp

                ! Consistency check 
                if (trim(mshlf%par%bmb_method) .eq. "anom") then 
                    write(error_unit,*) "marshelf_update:: Error: for bmb_method='anom', &
                    &tf_method=1 cannot be used."
                    write(error_unit,*) "bmb_method = ", trim(mshlf%par%bmb_method)
                    write(error_unit,*) "tf_method  = ", mshlf%par%tf_method
                    stop 
                end if 

                mshlf%now%tf_shlf = (mshlf%now%T_shlf - mshlf%now%T_fp_shlf) &
                            + mshlf%now%tf_corr + mshlf%now%tf_corr_basin

            case(2)
                ! Calculate the thermal forcing specific to the anom method 
                ! (assume that tf_shlf == dT_shlf == temp anomaly relative to a reference state)
                ! Note: for bmb_method='anom', tf_method automatically set to 2. 

                mshlf%now%tf_shlf = mshlf%now%dT_shlf + mshlf%now%tf_corr + mshlf%now%tf_corr_basin

            case DEFAULT 

                write(*,*) "marshelf_update:: Error: tf_method not recognized."
                write(*,*) "tf_method = ", mshlf%par%tf_method 
                stop 

        end select 
        

        ! Check consistency of tf_shlf field here. Make sure that 
        ! warmer water inland cannot be accessed below grounded ice. 

        ! To do: maybe this is not necessary?? ajr, 2022-02-10



        ! 3. Calculate current ice shelf bmb field (grounded-ice bmb is
        ! calculated in ice-sheet model separately) ========

        select case(trim(mshlf%par%bmb_method))

            case("pico") 
                ! Calculate bmb_shlf using the PICO box model 

                call pico_update(mshlf%pico,mshlf%now%T_shlf,mshlf%now%S_shlf, &
                                    H_ice,z_bed,f_grnd,z_sl,basins,mshlf%now%mask_ocn,dx)

                ! jablasco: to avoid ice shelves growing at the margin lets impose an averaged melt in region 2.1
                select case(trim(mshlf%par%domain))
                    case("Antarctica")
                        where (regions .eq. 2.1) mshlf%now%bmb_shlf = 0.5*mshlf%now%bmb_shlf+0.5*mshlf%par%c_deep
                end select 

            case("lin","quad","quad-nl","lin-slope", &
                    "quad-slope","quad-nl-slope","anom") 
                ! Calculate bmb_shlf using other available parameterizations 

                ! Check whether slope is being applied here 
                if (index(trim(mshlf%par%bmb_method),"slope") .gt. 0) then 
                    with_slope = .TRUE. 
                else 
                    with_slope = .FALSE. 
                end if 

                ! For simplicity, first calculate everywhere (ocn and grounded points)
                select case(trim(mshlf%par%bmb_method))

                    case("lin","lin-slope")

                        if (with_slope) then 
                            ! Scale gamma by gamma_prime*sin(theta)
                            gamma2D = mshlf%par%gamma_lin * &
                                        (mshlf%par%gamma_prime*mshlf%now%slope_base)
                        else 
                            gamma2D = mshlf%par%gamma_lin
                        end if 

                        call calc_bmb_linear(mshlf%now%bmb_shlf,mshlf%now%tf_shlf, &
                                                                        gamma2D,omega)
                            
                    case("quad","quad-slope")

                        if (with_slope) then 
                            ! Scale gamma by gamma_prime*sin(theta)
                            gamma2D = mshlf%par%gamma_quad * &
                                        (mshlf%par%gamma_prime*mshlf%now%slope_base)
                        else 
                            gamma2D = mshlf%par%gamma_quad
                        end if 
                        
                        call calc_bmb_quad(mshlf%now%bmb_shlf,mshlf%now%tf_shlf, &
                                                                        gamma2D,omega)
                        
                    case("quad-nl","quad-nl-slope")

                        if (with_slope) then 
                            ! Scale gamma by gamma_prime*sin(theta)
                            gamma2D = mshlf%par%gamma_quad_nl * &
                                        (mshlf%par%gamma_prime*mshlf%now%slope_base)
                        else 
                            gamma2D = mshlf%par%gamma_quad_nl
                        end if 
                        
                        ! Calculate basin-average thermal forcing 
                        call calc_variable_basin(mshlf%now%tf_basin,mshlf%now%tf_shlf, &
                                                        f_grnd,basins,H_ice,mshlf%now%mask_ocn)

                        ! Ensure tf_basin is non-negative following Lipscomb et al (2021)
                        where(mshlf%now%tf_basin .lt. 0.0_wp) mshlf%now%tf_basin = 0.0_wp 

                        call calc_bmb_quad_nl(mshlf%now%bmb_shlf,mshlf%now%tf_shlf,mshlf%now%tf_basin, &
                            mshlf%par%gamma_quad_nl,omega)
                    
                    case("anom")

                        call calc_bmb_anom(mshlf%now%bmb_shlf,mshlf%now%tf_shlf,mshlf%now%bmb_ref, &
                                    mshlf%par%kappa_grz,mshlf%par%c_grz,mshlf%par%f_grz_shlf, &
                                    mshlf%now%mask_ocn,mshlf%par%grz_length,dx)

                end select 

                ! Apply basin correction
                mshlf%now%bmb_shlf = mshlf%now%bmb_shlf + mshlf%now%bmb_corr
                
            case DEFAULT 

                write(*,*) "marshelf_update:: Error: bmb_method not recognized."
                write(*,*) "bmb_method = ", trim(mshlf%par%bmb_method)
                stop 

        end select 

        ! === Apply limitations ====

        ! The above routines calculate bmb_shlf everywhere with no limitations
        ! (ie, also for grounded ice and lakes if tf_shlf is defined). 
        ! Below we should apply specific limitations. 

        ! Apply maximum refreezing rate 
        where (mshlf%now%bmb_shlf .gt. mshlf%par%bmb_max)   &
                                mshlf%now%bmb_shlf = mshlf%par%bmb_max

        ! ajr: disable for testing!!
        ! ! Set bmb to zero for grounded and lake points 
        ! where (mshlf%now%mask_ocn .eq. mask_val_land .or. &
        !        mshlf%now%mask_ocn .eq. mask_val_lake)     &
        !                                     mshlf%now%bmb_shlf = 0.0 

        ! Ensure that ice accretion only occurs where ice exists 
        where (mshlf%now%bmb_shlf .gt. 0.0 .and. H_ice .eq. 0.0) &
                                            mshlf%now%bmb_shlf = 0.0 
        
        ! Apply c_deep melt value to deep ocean points
        ! n_smth can be zero when c_deep is a relatively small value (like -1 m/a). For
        ! c_deep = -50 m/a, n_smth=2 is more appropriate to ensure a smooth transition to 
        ! high melt rates. 
        call apply_c_deep(mshlf%par,mshlf%now%bmb_shlf,mshlf%now%mask_ocn,z_bed,z_sl,n_smth=0) 

        return
        
    end subroutine marshelf_update

    subroutine marshelf_init(mshlf,filename,group,nx,ny,domain,grid_name,regions,basins)

        implicit none 

        type(marshelf_class), intent(OUT) :: mshlf
        character(len=*), intent(IN)      :: filename
        character(len=*), intent(IN)      :: group
        integer, intent(IN)               :: nx, ny 
        character(len=*), intent(IN)      :: domain, grid_name
        real(wp), intent(IN)            :: regions(:,:)
        real(wp), intent(IN)            :: basins(:,:)
        
        ! Local variables
        integer  :: j 
        integer  :: num
        character(len=56) :: group_now
        real(wp) :: basin_number_now
        real(wp) :: tf_corr_now 

        ! Load parameters
        call marshelf_par_load(mshlf%par,filename,group,domain,grid_name)

        ! Allocate the object 
        call marshelf_allocate(mshlf%now,nx,ny)
        
        ! Decide whether to load observations or use psuedo-observations
        if (mshlf%par%use_obs) then

            ! Load observed field from file
            call nc_read(mshlf%par%obs_path,mshlf%par%obs_name,mshlf%now%bmb_ref) 

            ! Make negative since bmb is defined with melt negative, but observed field has melt positive
            mshlf%now%bmb_ref = -mshlf%now%bmb_ref

            ! Scale the obs as needed 
            mshlf%now%bmb_ref = mshlf%now%bmb_ref*mshlf%par%obs_scale 

            ! Limit the max obs values as desired 
            where (mshlf%now%bmb_ref >  mshlf%par%obs_lim) mshlf%now%bmb_ref =  mshlf%par%obs_lim 
            where (mshlf%now%bmb_ref < -mshlf%par%obs_lim) mshlf%now%bmb_ref = -mshlf%par%obs_lim 
        
        else 

            ! Set default psuedo-observation value since none are available and make them negative (as above) 
            mshlf%now%bmb_ref = -1.0*mshlf%par%obs_scale 

        end if 

        ! Initialize basin-wide correction fields
        mshlf%now%bmb_corr      = 0.0
        mshlf%now%tf_corr_basin = 0.0

        ! Define basin-wide corrections as needed     
        select case(trim(mshlf%par%corr_method))
            
            case("bmb")
                ! Modify specific basins according to parameter values 

                do j = 1, size(mshlf%par%basin_number)
                    where(basins .eq. mshlf%par%basin_number(j)) &
                        mshlf%now%bmb_corr = mshlf%par%basin_bmb_corr(j)
                end do 

            case("tf")
                ! Modify specific basins according to parameter values 

                do j = 1, size(mshlf%par%basin_number)

                    where(basins .eq. mshlf%par%basin_number(j)) &
                        mshlf%now%tf_corr_basin = mshlf%par%basin_tf_corr(j)

                end do

            case("tf-grl") 
                ! Modify specific basins according to parameter values 
                ! as defined in the parameter section 'tf_corr_grl'

                group_now = "tf_corr_grl" 

                ! ne = northeast
                call nml_read(filename,group_now,"ne",tf_corr_now)
                call apply_tf_corr_by_basin(mshlf%now%tf_corr_basin,basins,tf_corr_now, &
                                                basin_numbers=[2.0_wp])
                
                ! e = east
                call nml_read(filename,group_now,"e",tf_corr_now)
                call apply_tf_corr_by_basin(mshlf%now%tf_corr_basin,basins,tf_corr_now, &
                                                basin_numbers=[3.0_wp])
                
                
                ! se = southeast 
                call nml_read(filename,group_now,"se",tf_corr_now)
                call apply_tf_corr_by_basin(mshlf%now%tf_corr_basin,basins,tf_corr_now, &
                                                basin_numbers=[4.0_wp])

                
                ! w = west
                call nml_read(filename,group_now,"w",tf_corr_now)
                call apply_tf_corr_by_basin(mshlf%now%tf_corr_basin,basins,tf_corr_now, &
                                                basin_numbers=[6.0_wp,7.0_wp,8.0_wp])
                

            case("tf-ant") 
                ! Modify specific basins according to parameter values 
                ! as defined in the parameter section 'tf_corr_ant'

                group_now = "tf_corr_ant" 

                ! Ronne 
                call nml_read(filename,group_now,"ronne",tf_corr_now)
                call apply_tf_corr_by_basin(mshlf%now%tf_corr_basin,basins,tf_corr_now, &
                                                basin_numbers=[1.0_wp])
                
                ! Ross
                call nml_read(filename,group_now,"ross",tf_corr_now)
                call apply_tf_corr_by_basin(mshlf%now%tf_corr_basin,basins,tf_corr_now, &
                                                basin_numbers=[12.0_wp])
                
                ! Pine Island 
                call nml_read(filename,group_now,"pine",tf_corr_now)
                call apply_tf_corr_by_basin(mshlf%now%tf_corr_basin,basins,tf_corr_now, &
                                                basin_numbers=[14.0_wp])
                
                ! Abbott 
                call nml_read(filename,group_now,"abbott",tf_corr_now)
                call apply_tf_corr_by_basin(mshlf%now%tf_corr_basin,basins,tf_corr_now, &
                                                basin_numbers=[15.0_wp])
                
            case DEFAULT ! eg, "none", "None", "zero"

                ! DO NOTHING

        end select

        ! Load tf_corr field from file if desired 
        if (mshlf%par%tf_correction) then 

            call nc_read(mshlf%par%tf_path,mshlf%par%tf_name,mshlf%now%tf_corr)

        end if  

        ! ==============================================
        ! Generate reference ocean mask 
        ! (0: land, 1: open ocean, 2: deep ocean) 

        ! Define mask_ocn_ref based on regions mask 
        ! (these definitions should work for all North and Antarctica domains)
        mshlf%now%mask_ocn_ref = mask_val_land 
        where (regions .eq. 1.0_wp) mshlf%now%mask_ocn_ref = mask_val_ocean
        where (regions .eq. 2.0_wp) mshlf%now%mask_ocn_ref = mask_val_ocean

        select case(trim(mshlf%par%domain))

            case("Greenland") 
                ! Greenland specific ocean kill regions

                where (regions .ne. 1.3) mshlf%now%mask_ocn_ref = mask_val_deep_ocean

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
                where (regions .ne. 2.11) mshlf%now%mask_ocn_ref = mask_val_deep_ocean

            ! case("North") 
            !     ! North specific ocean kill regions

            !     ! Apply only in purely open-ocean regions  
            !     where (regions .eq. 1.0) mshlf%now%mask_ocn_ref = mask_val_deep_ocean

            case DEFAULT 
                ! Other domains: c_deep potentially applied everywhere 
                ! with deep ocean points 

                mshlf%now%mask_ocn_ref = mask_val_deep_ocean

        end select 

        ! ==============================================
        ! PICO 

        if (trim(mshlf%par%bmb_method) .eq. "pico") then 

            ! Initialize pico too 

            call pico_init(mshlf%pico,filename,nx,ny,domain)

        end if

        ! ====================================
        !
        ! Summary
        !
        ! ====================================
        
        write(*,*) "range bmb_ref:      ", minval(mshlf%now%bmb_ref), maxval(mshlf%now%bmb_ref)
        write(*,*) "range mask_ocn_ref: ", minval(mshlf%now%mask_ocn_ref), maxval(mshlf%now%mask_ocn_ref)
        
        ! Initialize variables 
        mshlf%now%bmb_shlf      = 0.0
        mshlf%now%tf_shlf       = 0.0
        mshlf%now%tf_basin      = 0.0 

        return 

    end subroutine marshelf_init

    subroutine apply_tf_corr_by_basin(tf_corr,basins,tf_corr_now,basin_numbers)
        ! Apply the value of tf_corr_now in the basins that correspond
        ! to the given basin_numbers of interest. 

        implicit none

        real(wp), intent(INOUT) :: tf_corr(:,:) 
        real(wp), intent(IN)    :: basins(:,:) 
        real(wp), intent(IN)    :: tf_corr_now
        real(wp), intent(IN)    :: basin_numbers(:) 

        ! Local variables 
        integer :: b, nb 
        real(wp) :: basin_number_now

        nb = size(basin_numbers)

        do b = 1, nb 
            basin_number_now = basin_numbers(b) 
            where(basins .eq. basin_number_now) tf_corr = tf_corr_now
        end do 

        return

    end subroutine apply_tf_corr_by_basin

    subroutine marshelf_end(mshlf)

        implicit none 

        type(marshelf_class) :: mshlf 

        ! Deallocate mshlf state object
        call marshelf_deallocate(mshlf%now)

        return 

    end subroutine marshelf_end


    subroutine marshelf_par_load(par,filename,group,domain,grid_name,init)

        type(marshelf_param_class), intent(OUT) :: par 
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: group
        character(len=*), intent(IN) :: domain, grid_name  
        logical, optional :: init 
        logical :: init_pars 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 

        call nml_read(filename,group,"bmb_method",     par%bmb_method,     init=init_pars)
        call nml_read(filename,group,"tf_method",      par%tf_method,      init=init_pars)
        call nml_read(filename,group,"interp_depth",   par%interp_depth,   init=init_pars)
        call nml_read(filename,group,"interp_method",  par%interp_method,  init=init_pars)
        call nml_read(filename,group,"find_ocean",     par%find_ocean,     init=init_pars)
        call nml_read(filename,group,"corr_method",    par%corr_method,    init=init_pars)   
        call nml_read(filename,group,"basin_number",   par%basin_number,   init=init_pars)
        call nml_read(filename,group,"basin_bmb_corr", par%basin_bmb_corr, init=init_pars)
        call nml_read(filename,group,"basin_tf_corr",  par%basin_tf_corr,  init=init_pars)       
        call nml_read(filename,group,"tf_correction",  par%tf_correction,  init=init_pars)
        call nml_read(filename,group,"tf_path",        par%tf_path,        init=init_pars)
        call nml_read(filename,group,"tf_name",        par%tf_name,        init=init_pars)
 
        call nml_read(filename,group,"bmb_max",        par%bmb_max,        init=init_pars)
        call nml_read(filename,group,"c_deep",         par%c_deep,         init=init_pars)
        call nml_read(filename,group,"depth_deep",     par%depth_deep,     init=init_pars)
        call nml_read(filename,group,"depth_const",    par%depth_const,    init=init_pars)
        call nml_read(filename,group,"depth_min",      par%depth_min,      init=init_pars)
        call nml_read(filename,group,"depth_max",      par%depth_max,      init=init_pars)    
        
        ! Freezing point 
        call nml_read(filename,group,"lambda1",        par%lambda1,        init=init_pars)
        call nml_read(filename,group,"lambda2",        par%lambda2,        init=init_pars)
        call nml_read(filename,group,"lambda3",        par%lambda3,        init=init_pars)
        
        ! bmb_method == [lin,quad,quad-nl]
        call nml_read(filename,group,"gamma_lin",      par%gamma_lin,      init=init_pars)
        call nml_read(filename,group,"gamma_quad",     par%gamma_quad,     init=init_pars)
        call nml_read(filename,group,"gamma_quad_nl",  par%gamma_quad_nl,  init=init_pars)
        call nml_read(filename,group,"gamma_prime",    par%gamma_prime,     init=init_pars)
        
        ! bmb_method == anom
        call nml_read(filename,group,"kappa_grz",      par%kappa_grz,      init=init_pars)
        call nml_read(filename,group,"c_grz",          par%c_grz,          init=init_pars)
        call nml_read(filename,group,"f_grz_shlf",     par%f_grz_shlf,     init=init_pars)
        call nml_read(filename,group,"grz_length",     par%grz_length,     init=init_pars)
        
        call nml_read(filename,group,"use_obs",        par%use_obs,        init=init_pars)
        call nml_read(filename,group,"obs_path",       par%obs_path,       init=init_pars)
        call nml_read(filename,group,"obs_name",       par%obs_name,       init=init_pars)
        call nml_read(filename,group,"obs_scale",      par%obs_scale,      init=init_pars)
        call nml_read(filename,group,"obs_lim",        par%obs_lim,        init=init_pars)
       
        ! Determine derived parameters
        call parse_path(par%obs_path,domain,grid_name)
        call parse_path(par%tf_path,domain,grid_name)
        par%domain = trim(domain)

        ! Consistency checks 
        if (par%f_grz_shlf .eq. 0.0) then 
            write(*,*) "marshelf_par_load:: Error: f_grz_shlf cannot be zero."
            stop 
        end if 

        if (trim(par%bmb_method) .eq. "pico") then 
            ! pico far-field ocean fields should be interpolated at the ocean bed (z_bed)
            par%interp_depth  = "bed"
            par%interp_method = "interp"
        end if 

        return

    end subroutine marshelf_par_load

    subroutine apply_c_deep(par,bmb,mask_ocn,z_bed,z_sl,n_smth)
        ! Apply c_deep for killing ice in the deep ocean 
        
        implicit none 

        type(marshelf_param_class) :: par 
        real(wp), intent(INOUT)    :: bmb(:,:)  
        integer,  intent(IN)       :: mask_ocn(:,:) 
        real(wp), intent(IN)       :: z_bed(:,:) 
        real(wp), intent(IN)       :: z_sl(:,:)
        integer,  intent(IN)       :: n_smth       ! Smoothing neighborhood radius in grid points

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

        where(mask_ocn .eq. mask_val_deep_ocean .and. z_sl-z_bed .ge. par%depth_deep) bmb = par%c_deep 

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

    subroutine calc_shelf_slope(slope,z_base,dx)
        ! Calculate the slope of the ice shelf base centered on aa-node, 
        ! where slope = sin(theta) = height/hypotenuse

        implicit none 

        real(wp), intent(OUT) :: slope(:,:) 
        real(wp), intent(IN)  :: z_base(:,:) 
        real(wp), intent(IN)  :: dx 

        ! Local variables 
        integer :: i, j, nx, ny, ip1, im1, jp1, jm1 
        real(wp) :: dz, dx_2, hypot

        nx = size(slope,1) 
        ny = size(slope,2) 

        dx_2 = dx*2.0_wp 

        do j = 1, ny 
        do i = 1, nx 
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1)
            jp1 = min(j+1,ny)

            ! Get mean vertical change in distance between x and y directions
            dz =  0.5_wp*(z_base(ip1,j)-z_base(im1,j)) &
                + 0.5_wp*(z_base(i,jp1)-z_base(i,jm1))

            ! Get hypotenuse 
            hypot = sqrt(dz**2 + dx_2**2) 

            ! Get slope == sin(theta) magnitude, independent of direction
            slope(i,j) = abs(dz) / hypot 

        end do 
        end do
        
        return 

    end subroutine calc_shelf_slope

    subroutine set_ocean_mask(mask_ocn,f_grnd,mask_ocn_ref,find_ocean)
        ! Define a mask of topographic properties
        ! 0: Land
        ! 1: Grounding line (grounded side, including partially floating points)
        ! 2: Floating line  (first fully floating points adjacent to grounding line)
        ! 3: Open ocean 
        ! 4: Lakes 

        implicit none

        integer,  intent(OUT) :: mask_ocn(:,:) 
        real(wp), intent(IN)  :: f_grnd(:,:)
        integer,  intent(IN)  :: mask_ocn_ref(:,:)
        logical,  intent(IN)  :: find_ocean         ! Should open-ocean and lake points be determined?

        ! Local variables 
        integer :: i, j, nx, ny  
        integer :: im1, ip1, jm1, jp1 

        integer, allocatable :: mask_open_ocn(:,:) 

        nx = size(mask_ocn,1)
        ny = size(mask_ocn,2) 

        ! First, by default set everything to open ocean
        mask_ocn = mask_val_ocean 

        do j = 1, ny 
        do i = 1, nx 
        
            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            ! Determine cases (land, grounding line, floating line)

            if (f_grnd(i,j) .gt. 0.0 .and. &
                (f_grnd(im1,j) .eq. 0.0 .or. f_grnd(ip1,j) .eq. 0.0 .or. &
                 f_grnd(i,jm1) .eq. 0.0 .or. f_grnd(i,jp1) .eq. 0.0) ) then 
                ! Grounded point or partially floating point with fully floating neighbors
                
                mask_ocn(i,j) = mask_val_grounding_line

            else if (f_grnd(i,j) .eq. 0.0 .and. &
                (f_grnd(im1,j) .eq. 1.0 .or. f_grnd(ip1,j) .eq. 1.0 .or. &
                 f_grnd(i,jm1) .eq. 1.0 .or. f_grnd(i,jp1) .eq. 1.0) ) then 
                ! Fully floating point with partially floating or grounded neighbors
                
                mask_ocn(i,j) = mask_val_floating_line

            else if (f_grnd(i,j) .eq. 1.0) then 
                ! Land point 

                mask_ocn(i,j) = mask_val_land

            else if (f_grnd(i,j) .eq. 0.0 .and. mask_ocn_ref(i,j) .eq. mask_val_deep_ocean) then

                mask_ocn(i,j) = mask_val_deep_ocean

            else

                mask_ocn(i,j) = mask_val_ocean 

            end if 
            
        end do 
        end do 

        ! Diagnose further cases

            
        if (find_ocean) then 
            ! If desired, determine which floating or partially 
            ! floating points are connected to the open-ocean
            ! (otherwise all ocean points are assumed to be open-ocean points)

            call find_open_ocean_and_lakes(mask_ocn,mask_ocn_ref)

        end if 

        return

    end subroutine set_ocean_mask


    subroutine find_open_ocean_and_lakes(mask_ocn,mask_ocn_ref)
        ! Brute-force routine to find all ocean points 
        ! (ie, when f_grnd < 1) connected to
        ! the open ocean as defined in mask_ocn_ref. 

        ! Mask returned should be the same as original mask_ocn,
        ! except with floating points and ground-line points that are not
        ! connected to the open ocean set to lakes. 

        implicit none 

        integer, intent(INOUT) :: mask_ocn(:,:)
        integer, intent(IN)    :: mask_ocn_ref(:,:) 

        ! Local variables 
        integer :: i, j, q, nx, ny
        integer :: im1, ip1, jm1, jp1 
        integer :: n_unfilled 

        integer, allocatable :: mask(:,:) 
        integer, allocatable :: mask0(:,:) 

        integer, parameter :: qmax = 1000 

        nx = size(mask_ocn,1)
        ny = size(mask_ocn,2) 

        ! Allocate local mask object for diagnosing points of interest
        allocate(mask(nx,ny))
        allocate(mask0(nx,ny))

        ! Initially assume all points are 'closed ocean' points (mask==-1)
        mask = -1 

        ! Set purely land points to zero 
        where (mask_ocn .eq. mask_val_land) mask = 0 

        ! Set grounding-line points to specific value too (mask==-2)
        where (mask_ocn .eq. mask_val_grounding_line)  mask = -2

        ! Now populate our ocean mask to be consistent 
        ! with known open ocean points that can be 
        ! normal open ocean (1) or deep ocean (2) 
        where (mask .eq. -1 .and. &
            (mask_ocn_ref .eq. mask_val_ocean .or. &
             mask_ocn_ref .eq. mask_val_deep_ocean) )  mask = 1 
         
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

        ! Update the ocean mask to reflect open ocean and lakes
        ! (deep ocean remains deep ocean always)
        where(mask .eq. -1.0 .and. mask_ocn .eq. mask_val_ocean) &
                                            mask_ocn = mask_val_lake

        where(mask .eq. -1.0 .and. mask_ocn .eq. mask_val_floating_line) &
                                            mask_ocn = mask_val_lake

        ! Finally, check which grounding-line points should also
        ! be connected to the open ocean.

        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices 
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1)
            jp1 = min(j+1,ny)
            
            if (mask(i,j) .eq. -2) then 
                ! This is a grounding-line point.
                ! If no neighbors are open ocean, set it to a lake point

                if ( mask(im1,j) .ne. 1 .and. &
                     mask(ip1,j) .ne. 1 .and. &
                     mask(i,jm1) .ne. 1 .and. &
                     mask(i,jp1) .ne. 1 ) then 
                    ! No neighbors are open ocean, set
                    ! this point to 'closed ocean' (lake)

                    mask_ocn(i,j) = mask_val_lake

                end if 
                
            end if 

        end do
        end do

        return 

    end subroutine find_open_ocean_and_lakes

    subroutine calc_variable_basin(var_basin,var2D,f_grnd,basins,H_ice,mask_ocn)
        ! Compute mean ocean temperature of ice-shelf basin

        implicit none

        real(wp), intent(OUT) :: var_basin(:,:) 
        real(wp), intent(IN)  :: var2D(:,:) 
        real(wp), intent(IN)  :: f_grnd(:,:)
        real(wp), intent(IN)  :: basins(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        integer,  intent(IN)  :: mask_ocn(:,:)

        ! Local variables
        integer :: i, j, m, nx, ny, npts 
        real(wp) :: var_mean

        nx = size(var_basin,1)
        ny = size(var_basin,2) 

        ! Initially assign to shelf values real ocean values
        var_basin = var2D 

        ! Loop over each basin
        do m=1, int(maxval(basins))
            
            ! First calculate the basin-wide average variable value,
            ! limited to floating ice shelf points

            var_mean = 0.0
            npts     = 0
            
            do j = 1, ny
            do i = 1, nx
                
                ! Grounding-line or floating ice point, or if none available,
                ! then any floating line points available in the basin.
                if (basins(i,j) .eq. m .and. f_grnd(i,j) .lt. 1.0 .and.  &
                     (H_ice(i,j) .gt. 0.0 .or. &
                      mask_ocn(i,j) .eq. mask_val_floating_line .or. &
                      mask_ocn(i,j) .eq. mask_val_grounding_line) ) then
                    
                    var_mean = var_mean + var2D(i,j)
                    npts     = npts + 1

                end if
            
            end do
            end do

             
            if (npts .gt. 0) then
                ! Get mean basin value
                
                var_mean = var_mean / real(npts,wp)

            else
                ! Set to zero for now and print a warning for safety.
                ! This case is unlikely to happen.

                var_mean = 0.0_wp 
                
                write(*,*) "Warning:: calc_variable_basin:: no floating ice points &
                &available in this basin for calculating the basin-wide mean."
                write(*,*) "basin = ", m 
                write(*,*) "n_flt = ", count(basins .eq. m .and. &
                            (f_grnd .lt. 1.0 .and. H_ice .gt. 0.0) )
                write(*,*) "n_ice = ", count(basins .eq. m .and. &
                            (H_ice .gt. 0.0) )
            end if 

            ! Assign shelf value to all points in the basin
            do j = 1, ny
            do i = 1, nx
                
                ! Basin point
                if (basins(i,j) .eq. m) then
                    var_basin(i,j) = var_mean
                end if
            
            end do
            end do

        end do

        return

    end subroutine calc_variable_basin

    subroutine calc_shelf_variable_mean(wt_shlf,depth,depth_range)
        ! Calculate average variable value for a given range of depths
        ! at a specific point (x,y). 

        implicit none 
        
        real(wp), intent(OUT)   :: wt_shlf(:)
        real(wp), intent(IN)    :: depth(:)
        real(wp), intent(IN)    :: depth_range(:)

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

        ! Get index weights to produce mean water temperature for these depths
        wt_shlf        = 0.0 
        wt_shlf(k0:k1) = 1.0 

        return 

    end subroutine calc_shelf_variable_mean

    subroutine calc_shelf_variable_layer(wt_shlf,depth,depth_shlf)
        ! Calculates the water temperature at the depth of the ice shelf
        ! It assigns the temperature of the nearest layer
        
        implicit none

        real(wp), intent(OUT)   :: wt_shlf(:)
        real(wp), intent(IN)    :: depth(:)
        real(wp), intent(IN)    :: depth_shlf

        ! Local variables
        integer :: k0, i, j

        ! Determine layer closest to target depth 
        k0 = minloc(abs(depth-depth_shlf),dim=1)

        ! Assign weighting function 
        wt_shlf     = 0.0 
        wt_shlf(k0) = 1.0 

        return

    end subroutine calc_shelf_variable_layer

    subroutine calc_shelf_variable_depth(wt_shlf,depth,depth_shlf)
        ! Calculates the water temperature from linear interpolation of vertical profile
        
        implicit none

        real(wp), intent(OUT)   :: wt_shlf(:)
        real(wp), intent(IN)    :: depth(:)
        real(wp), intent(IN)    :: depth_shlf

        ! Local variables 
        integer  :: j, n 
        real(wp) :: alpha 

        n    = size(depth) 

        do j = 1, n 
            if (depth(j) .ge. depth_shlf) exit 
        end do

        if (j .eq. 1) then 
            wt_shlf    = 0.0 
            wt_shlf(1) = 1.0 
        else if (j .eq. n+1) then 
            wt_shlf    = 0.0 
            wt_shlf(n) = 1.0 
        else 
            alpha = (depth_shlf - depth(j-1)) / (depth(j) - depth(j-1))
            wt_shlf      = 0.0 
            wt_shlf(j-1) = 1.0-alpha 
            wt_shlf(j)   = alpha
        end if  

        return

    end subroutine calc_shelf_variable_depth

    elemental subroutine calc_freezing_point(to_fp,so,z_base,lambda1,lambda2,lambda3,T_ref)
        ! Calculate the water freezing point following 
        ! Favier et al (2019), Eq. 3
        ! z_base is ice-base elevation, which is negative below sea level

        implicit none 

        real(wp), intent(OUT) :: to_fp
        real(wp), intent(IN)  :: so 
        real(wp), intent(IN)  :: z_base 
        real(wp), intent(IN)  :: lambda1 
        real(wp), intent(IN)  :: lambda2 
        real(wp), intent(IN)  :: lambda3 
        real(wp), intent(IN)  :: T_ref 

        to_fp = lambda1*so + lambda2 + lambda3*z_base + T_ref 

        return 

    end subroutine calc_freezing_point

    elemental subroutine calc_bmb_linear(bmb,tf,gamma,omega)
        ! Calculate the basal mass balance, linear method
        ! Favier et al (2019), Eq. 2

        ! Note: omega = (rho_sw*cp_o) / (rho_ice*L_ice)

        implicit none 

        real(wp), intent(OUT) :: bmb            ! [m i.e./yr] Basal mass balance 
        real(wp), intent(IN)  :: tf             ! [K] Thermal forcing 
        real(wp), intent(IN)  :: gamma          ! [m yr-1] Heat exchange velocity 
        real(wp), intent(IN)  :: omega          ! [1/K] Physical scaling constant 

        bmb = -gamma * omega * tf 

        return 

    end subroutine calc_bmb_linear

    elemental subroutine calc_bmb_quad(bmb,tf,gamma,omega)
        ! Calculate the basal mass balance, (local) quadratic method
        ! Favier et al (2019), Eq. 4

        implicit none 

        real(wp), intent(OUT) :: bmb            ! [m i.e./yr] Basal mass balance 
        real(wp), intent(IN)  :: tf             ! [K] Thermal forcing 
        real(wp), intent(IN)  :: gamma          ! [m yr-1] Heat exchange velocity 
        real(wp), intent(IN)  :: omega          ! [1/K] Physical scaling constant 

        bmb = -gamma * (omega*omega) * (tf*tf)

        return 

    end subroutine calc_bmb_quad
    
    elemental subroutine calc_bmb_quad_nl(bmb,tf,tf_basin,gamma,omega)
        ! Calculate the basal mass balance, quadratic-nonlocal method
        ! Favier et al (2019), Eq. 5

        implicit none 

        real(wp), intent(OUT) :: bmb            ! [m i.e./yr] Basal mass balance 
        real(wp), intent(IN)  :: tf             ! [K] Thermal forcing 
        real(wp), intent(IN)  :: tf_basin       ! [K] Basin thermal forcing 
        real(wp), intent(IN)  :: gamma          ! [m yr-1] Heat exchange velocity 
        real(wp), intent(IN)  :: omega          ! [1/K] Physical scaling constant 

        bmb = -gamma * (omega*omega) * (tf*tf_basin)

        return 

    end subroutine calc_bmb_quad_nl
    
    elemental subroutine calc_bmb_anom(bmb,tf,bmb_ref,kappa_grz,c_grz,f_grz_shlf,mask_ocn,grz_length,dx)
        ! Calculate the ice-shelf basal mass balance following an anomaly approach:
        !
        ! bmb = bmb_ref + kappa*dT 
        !
        ! First calculate bmb assuming grounding line parameters everywhere, then
        ! for the open shelf away from the grounding line, scale bmb according to 
        ! the prescribed ratio between grounding line melt and shelf melt (f_grz_shlf) 

        implicit none 

        real(wp), intent(OUT) :: bmb            ! [m i.e./yr] Basal mass balance 
        real(wp), intent(IN)  :: tf             ! [K] Thermal forcing 
        real(wp), intent(IN)  :: bmb_ref        ! [m i.e./yr] Reference basal mass balance 
        real(wp), intent(IN)  :: kappa_grz      ! [m yr-1 / K] Heat flux coeff. grounding zone 
        real(wp), intent(IN)  :: c_grz          ! [-] Scalar coefficient, grounding zone  
        real(wp), intent(IN)  :: f_grz_shlf     ! [-] Ratio of gl to shelf melt rate
        integer,  intent(IN)  :: mask_ocn       ! [-] Ocean mask
        real(wp), intent(IN)  :: grz_length     ! [km] Grounding-zone length scale
        real(wp), intent(IN)  :: dx             ! [m] Grid resolution 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(wp) :: bmb_grline 
        real(wp) :: bmb_floating 
        real(wp) :: f_shlf_grz 
        real(wp) :: grz_wt

        ! Get ratio of shelf melt to gr-line melt (inverse of parameter):
        f_shlf_grz = 1.0_wp / f_grz_shlf 

        ! Determine resolution scaling at the grounding line 
        grz_wt = min( (grz_length*1e3) / dx, 1.0) 
        
        ! Calculate grounding line bmb:
        bmb_grline = c_grz*bmb_ref - kappa_grz*tf

        ! Calculate floating shelf bmb: 
        bmb_floating = bmb_grline*f_shlf_grz

        if (mask_ocn .eq. mask_val_grounding_line .or. mask_ocn .eq. mask_val_floating_line) then 
            ! Calculate gl-bmb as weighted average between gl-line value and shelf value 
            ! to account for resolution dependence 
            
            bmb = (1.0-grz_wt)*bmb_floating + grz_wt*bmb_grline
        
        else
            ! Set bmb to floating value 

            bmb = bmb_floating 

        end if 

        return 

    end subroutine calc_bmb_anom
    
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
        allocate(now%bmb_ref(nx,ny))
        allocate(now%bmb_corr(nx,ny))
        
        allocate(now%T_shlf(nx,ny))
        allocate(now%dT_shlf(nx,ny))
        allocate(now%S_shlf(nx,ny))
        allocate(now%T_fp_shlf(nx,ny))
        
        allocate(now%tf_shlf(nx,ny))
        allocate(now%tf_basin(nx,ny))
        allocate(now%tf_corr(nx,ny))
        allocate(now%tf_corr_basin(nx,ny))

        allocate(now%z_base(nx,ny))
        allocate(now%slope_base(nx,ny))

        allocate(now%mask_ocn_ref(nx,ny))
        allocate(now%mask_ocn(nx,ny))

        ! Initialize variables 
        now%bmb_shlf        = 0.0  
        now%bmb_ref         = 0.0 
        now%bmb_corr        = 0.0       
 
        now%T_shlf          = 0.0
        now%dT_shlf         = 0.0
        now%S_shlf          = 0.0
        now%T_fp_shlf       = 0.0

        now%tf_shlf         = 0.0
        now%tf_basin        = 0.0
        now%tf_corr         = 0.0
        now%tf_corr_basin   = 0.0     

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
        if (allocated(now%bmb_shlf))        deallocate(now%bmb_shlf)
        if (allocated(now%bmb_ref))         deallocate(now%bmb_ref)
        if (allocated(now%bmb_corr))        deallocate(now%bmb_corr)       
 
        if (allocated(now%T_shlf))          deallocate(now%T_shlf)
        if (allocated(now%dT_shlf))         deallocate(now%dT_shlf)
        if (allocated(now%S_shlf))          deallocate(now%S_shlf)
        if (allocated(now%T_fp_shlf))       deallocate(now%T_fp_shlf)
        
        if (allocated(now%tf_shlf))         deallocate(now%tf_shlf)
        if (allocated(now%tf_basin))        deallocate(now%tf_basin)
        if (allocated(now%tf_corr))         deallocate(now%tf_corr)
        if (allocated(now%tf_corr_basin))   deallocate(now%tf_corr_basin)       
 
        if (allocated(now%z_base))          deallocate(now%z_base)
        if (allocated(now%slope_base))      deallocate(now%slope_base)
        
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
    
end module marine_shelf


