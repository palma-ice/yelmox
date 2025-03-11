module ice2ocean

    use ncio 

    implicit none
    
    private
    public :: calc_fwf

    integer,  parameter :: wp = kind(1.0)


contains

function calc_fwf(rho_water,rho_ice,sec_year,H_ice,dHidt,f_grnd,dx,dy,mask,hemisphere) result(fwf)
    ! This subroutine is based on yelmo_regions.f90 calc_yregions(...)

    implicit none 

    real(wp) :: fwf

    ! Local variables
    integer  :: nx, ny
    real(wp) :: npts_flt

    real(wp) :: m3_km3 = 1e-9 
    real(wp) :: conv_km3a_Sv

    real(wp)   :: rho_water, rho_ice ! bnd%c%rho_w, bnd%c%rho_ice
    real(wp)   :: sec_year ! bnd%c%sec_year
    real(wp), allocatable   :: H_ice(:,:), dHidt(:,:) ! tpo%now%H_ice, tpo%now%dHidt
    real(wp), allocatable   :: f_grnd(:,:)
    real(wp)           :: dx, dy ! tpo%par%dx, tpo%par%dy
    real(wp), allocatable :: mask(:,:) 
    character(len=*) :: hemisphere

    logical, allocatable :: mask_flt(:,:), basin_mask(:,:), extended_north_basin_mask(:,:)
    real(wp), allocatable :: H_af(:,:) 
    logical, parameter :: use_extended_north=.TRUE.

    character(len=*), parameter :: fwf_def= "mask_and_floating_ice" !"mask_and_floating_ice" 

    ! Conversion parameter 
    conv_km3a_Sv = 1e-6*(1e9*rho_water/rho_ice)/sec_year   

    ! Grid size 
    nx = size(H_ice,1)
    ny = size(H_ice,2)

    allocate(mask_flt(nx,ny))
    allocate(H_af(nx,ny)) 
    allocate(extended_north_basin_mask(nx,ny))

    ! Define masks
    mask_flt  = (H_ice .gt. 0.0 .and. f_grnd .eq. 0.0)
    npts_flt  = real(count(mask_flt),wp)

    ! ===== Compute fwf =====
    select case(fwf_def)
        case("floating_ice") 
            fwf = -sum(dHidt)*dx*dy*m3_km3*conv_km3a_Sv ! mirar a calcular la derivada de V
        case("mask_and_floating_ice")
            ! Mask values:
            !   1 = atlantic ocean
            !   2 = arctic ocean
            !   3 = subarctic ocean
            !   4 = indian ocean
            !   5 = pacific ocean
            !   6 = southern ocean

            select case(hemisphere)
                case("north")
                    if (use_extended_north) then
                        extended_north_basin_mask = (mask .eq. 1 .or. mask .eq. 3)
                        basin_mask = (mask_flt .and. extended_north_basin_mask)   ! floating points in the atlantic ocean area
                    else
                        basin_mask = (mask_flt .and. mask .eq. 1)   ! floating points in the atlantic ocean area
                    end if

                    if (real(count(basin_mask), wp) .gt. 0) then
                        fwf = max(-sum(dHidt, mask=basin_mask)*dx*dy*m3_km3*conv_km3a_Sv, 0.0)
                    else
                        fwf = 0.0
                    end if
                case("south")
                    basin_mask = (mask_flt .and. mask .eq. 6)   ! floating points in the southern ocean area
                    
                    if (real(count(basin_mask), wp) .gt. 0) then
                        fwf = max(-sum(dHidt, mask=basin_mask)*dx*dy*m3_km3*conv_km3a_Sv, 0.0) * 1/3    ! fwf is multiplied for the proportion of atlantic ocean in the domain
                    else
                        fwf = 0.0
                    end if
                case DEFAULT
                    basin_mask = mask_flt
                    if (real(count(basin_mask), wp) .gt. 0) then
                        fwf = max(-sum(dHidt, mask=basin_mask)*dx*dy*m3_km3*conv_km3a_Sv, 0.0)
                    else
                        fwf = 0.0
                    end if
            end select
           
        case DEFAULT
            fwf = max(-sum(dHidt)*dx*dy*m3_km3*conv_km3a_Sv, 0.0)
        end select
    return

end function calc_fwf

end module ice2ocean
