

module lithosphere 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: prec = sp 

    ! Densities
    real(prec), parameter :: rho_a   = 3300.d0   ! Density asthenosphere [kg/m^3]
    real(prec), parameter :: rho_ice =  910.d0   ! Density ice           [kg/m^3] 
    real(prec), parameter :: rho_w   = 1000.d0   ! Density water         [kg/m^3] 
    real(prec), parameter :: rho_sw  = 1028.d0   ! Density seawater      [kg/m^3] 
    
    real(prec), parameter :: sec_year = 365.d0*24.d0*60.d0*60.d0

    private 
    public :: calc_bedrock_llra_relaxed, calc_bedrock_llra

contains 

    subroutine calc_bedrock_llra_relaxed(zb0,zb,H_w,H_sw,H_ice)

        implicit none 

        real(prec), dimension(:,:), intent(out)   :: zb0 
        real(prec), dimension(:,:), intent(in)    :: zb, H_w, H_sw, H_ice
        
        ! Local variables
        real(prec), dimension(:,:), allocatable :: H_load 

        ! Allocate local array
        allocate(H_load(size(zb,1),size(zb,2)))

        ! Calculate the total load in thickness on top of aesthenosphere [m a.e.] 
        ! (ie, meters aesthenosphere equivalent)
        H_load = (rho_w/rho_a)*H_w + (rho_sw/rho_a)*H_sw + (rho_ice/rho_a)*H_ice 

        ! Calculate the relaxed lithosphere if the loading were removed 
        zb0 = zb + H_load 

        return 

    end subroutine calc_bedrock_llra_relaxed

    subroutine calc_bedrock_llra(dzbdt,zb0,zb,H_w,H_sw,H_ice,dt,tau,fr)

        ! Isostatic bedrock adjustment with local
        ! lithosphere and relaxing asthenosphere (LLRA)
        ! Calculates the rate of uplift of current bedrock elevation [m/a]
        ! given ice (H), basal water (H_w) and oceanic (H_sw) loading 

        implicit none 

        real(prec), dimension(:,:), intent(inout) :: dzbdt
        real(prec), dimension(:,:), intent(in)    :: zb0, zb, H_w, H_sw, H_ice
        real(prec), intent(in)                    :: dt, tau, fr 
        
        ! Local variables
        real(prec), dimension(:,:), allocatable :: H_load, zb_new 

        ! Check parameter ranges 
        if (fr .lt. 0.d0 .or. fr .gt. 1.d0) then 
            write(*,*) "lithosphere:: error: llra_f must be between 0 and 1."
            stop 
        end if 

        ! Allocate local arrays
        allocate(H_load(size(zb,1),size(zb,2)))
        allocate(zb_new(size(zb,1),size(zb,2)))

        ! Calculate the total load in thickness on top of aesthenosphere [m a.e.] 
        ! (ie, meters aesthenosphere equivalent)
        H_load = (rho_w/rho_a)*H_w + (rho_sw/rho_a)*H_sw + (rho_ice/rho_a)*H_ice 

        ! Calculate new bedrock elevation wrt loading (ice+ocean) [m]
        zb_new = ( zb*tau + dt*(zb0 - fr*H_load) ) / (dt+tau)

        ! Calculate bedrock uplift rate [m/a] 
        dzbdt = (zb_new - zb) / dt

        return 

    end subroutine calc_bedrock_llra

    
end module lithosphere

