module nautilus
    ! ============================================================================================= !
    ! Nautilus Ocean Box Model
    ! ---------------------------------------------------------------------------------------------
    ! Description:
    !   * Fortran 90 code based on Zickfield and Rahmstorf (2003, Ocean Dynamics) model 
    !   * The original model was adapted for paleoclimate studies by Jorge Alvarez-Solas (2019)
    !   * Now I (Sergio Pérez-Montero) am generalizing the code for the coupling with Yelmo (2024)
    !
    ! ====================== !

    use ncio
    use nml
    use obm_defs

    implicit none

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: prec = sp
    integer,  parameter :: wp   = sp 

    real(wp), parameter :: sec_year  = 365.0*24.0*60.0*60.0   ! [s/a]
    real(wp), parameter :: pi        = 3.14159265359

    private
    public :: nautilus_init 
    public :: nautilus_update 
    
    contains

    subroutine nautilus_init(ntls, filename, group)

        implicit none 

        type(obm_class), intent(INOUT) :: ntls
        character(len=*),  intent(IN) :: filename, group 

        ! Read parameters from namelist file
        call nml_read(filename, group, "specific_heat_capacity", ntls%par%specific_heat_capacity)
        call nml_read(filename, group, "density_of_seawater", ntls%par%density_of_seawater)
        call nml_read(filename, group, "therm_exp_coeff", ntls%par%therm_exp_coeff)
        call nml_read(filename, group, "haline_exp_coeff", ntls%par%haline_exp_coeff)
        call nml_read(filename, group, "reference_salinity", ntls%par%reference_salinity)
        call nml_read(filename, group, "thermal_coupling_constant", ntls%par%thermal_coupling_constant)
        call nml_read(filename, group, "emp_flow_k", ntls%par%emp_flow_k)
        call nml_read(filename, group, "vol_n", ntls%par%vol_n)
        call nml_read(filename, group, "vol_t", ntls%par%vol_t)
        call nml_read(filename, group, "vol_td", ntls%par%vol_td)
        call nml_read(filename, group, "vol_s", ntls%par%vol_s)
        call nml_read(filename, group, "depth_n", ntls%par%depth_n)
        call nml_read(filename, group, "depth_t", ntls%par%depth_t)
        call nml_read(filename, group, "depth_td", ntls%par%depth_td)
        call nml_read(filename, group, "depth_s", ntls%par%depth_s)

        ! Initialize prognostic states of the model
        call nml_read(filename, group, "tn_init", ntls%tn)
        call nml_read(filename, group, "tt_init", ntls%tt)
        call nml_read(filename, group, "ttd_init", ntls%ttd)
        call nml_read(filename, group, "ts_init", ntls%ts)
        call nml_read(filename, group, "sn_init", ntls%sn)
        call nml_read(filename, group, "st_init", ntls%st)
        call nml_read(filename, group, "std_init", ntls%std)
        call nml_read(filename, group, "ss_init", ntls%ss)
        call nml_read(filename, group, "fn_init", ntls%fn)
        call nml_read(filename, group, "ft_init", ntls%ft) 
        call nml_read(filename, group, "fs_init", ntls%fs) 
        call nml_read(filename, group, "phin_init", ntls%phin)
        call nml_read(filename, group, "phit_init", ntls%phit) 
        call nml_read(filename, group, "thetan_init", ntls%thetan)
        call nml_read(filename, group, "thetat_init", ntls%thetat) 
        call nml_read(filename, group, "thetas_init", ntls%thetas)
        call nml_read(filename, group, "lambdan_init", ntls%lambdan)
        call nml_read(filename, group, "lambdat_init", ntls%lambdat)
        call nml_read(filename, group, "lambdatd_init", ntls%lambdatd)  
        call nml_read(filename, group, "lambdas_init", ntls%lambdas)

        return

    end subroutine nautilus_init

    subroutine nautilus_update(ntls,dt)

        implicit none 

        type(obm_class), intent(INOUT) :: ntls  
        real(prec), intent(IN)    :: dt 
        real(wp), parameter :: factor = 365*24*60*60 * 1e6

        ! Local variables
        real(prec) :: dTsdt, dTndt, dTtdt, dTtddt 
        real(prec) :: dSsdt, dSndt, dStdt, dStddt 
        real(prec) :: fn, ft, fs

        ! Transform external fluxes (Sv) to m3a-1 (1 Sv = 1e6 m3s-1 = 1e6 * 365*24*60*60 m3yr-1)
        fn = ntls%fn * factor
        ft = ntls%ft * factor
        fs = ntls%fs * factor

        ! Update diagnostic variables
        ntls%m = calc_overturning(ntls%par%emp_flow_k, ntls%par%haline_exp_coeff, ntls%par%therm_exp_coeff, ntls%ts, ntls%ss, ntls%tn, ntls%sn)
        call calc_all_lambda(ntls)

        ! Compute time derivatives, box loop: [...-> South --> Tropics --> North --> Tropics Deep -...]
        dTsdt = calc_temperature_derivative(ntls%par%vol_s, ntls%lambdas, ntls%m, -ntls%ts, ntls%ttd, ntls%thetas)               
        dTtdt = calc_temperature_derivative(ntls%par%vol_t, ntls%lambdat, ntls%m, -ntls%tt, ntls%ts, ntls%thetat)
        dTndt = calc_temperature_derivative(ntls%par%vol_n, ntls%lambdan, ntls%m, -ntls%tn, ntls%tt, ntls%thetan)                                          
        dTtddt = calc_temperature_derivative(ntls%par%vol_td, 0.0, ntls%m, -ntls%ttd, ntls%tn, 0.0)   
        
        dSsdt = calc_salinity_derivative(ntls%par%vol_s, ntls%par%reference_salinity, ntls%m, -ntls%ss, ntls%std, ntls%phit, -fs)
        dStdt = calc_salinity_derivative(ntls%par%vol_t, ntls%par%reference_salinity, ntls%m, -ntls%st, ntls%ss, ntls%phin, -ntls%phit)
        dSndt = calc_salinity_derivative(ntls%par%vol_n, ntls%par%reference_salinity, ntls%m, -ntls%sn, ntls%st, -fn, -ntls%phin)
        dStddt = calc_salinity_derivative(ntls%par%vol_td, 0.0, ntls%m, -ntls%std, ntls%sn, 0.0, 0.0)

        ! Update prognostic variables
        ntls%ts = update_state(ntls%ts, dTsdt, dt)
        ntls%tt = update_state(ntls%tt, dTtdt, dt)
        ntls%tn = update_state(ntls%tn, dTndt, dt)
        ntls%ttd = update_state(ntls%ttd, dTtddt, dt)

        ntls%ss = update_state(ntls%ss, dSsdt, dt)
        ntls%st = update_state(ntls%st, dStdt, dt)
        ntls%sn = update_state(ntls%sn, dSndt, dt)
        ntls%std = update_state(ntls%std, dStddt, dt)

        return

    end subroutine nautilus_update

    subroutine calc_all_lambda(ntls)
        implicit none 

        type(obm_class), intent(INOUT) :: ntls 

        ntls%lambdas = calc_individual_thermal_coupling_constant(ntls%par%thermal_coupling_constant, &
                                                                 ntls%par%specific_heat_capacity, ntls%par%density_of_seawater, &
                                                                 ntls%par%depth_s)
        ntls%lambdat = calc_individual_thermal_coupling_constant(ntls%par%thermal_coupling_constant, &
                                                                 ntls%par%specific_heat_capacity, ntls%par%density_of_seawater, &
                                                                 ntls%par%depth_t)
        ntls%lambdan = calc_individual_thermal_coupling_constant(ntls%par%thermal_coupling_constant, &
                                                                 ntls%par%specific_heat_capacity, ntls%par%density_of_seawater, &
                                                                 ntls%par%depth_n)                                  

        return

    end subroutine calc_all_lambda

    function update_state(u, dudt, deltat) result(newu)
        real(prec) :: u, dudt, deltat, newu
        newu = u + dudt * deltat
        return 
    end function update_state

    function calc_overturning(k, beta, alpha, t1, s1, t2, s2) result(m)
        real(prec) :: k, beta, alpha, t1, s1, t2, s2, m
        m = k * (beta * (s2 - s1) - alpha * (t2 - t1))
        return
    end function calc_overturning

    function calc_temperature_derivative(vol1, lambda1, m, t1, t2, theta1) result(dTdt)
        ! WARNING: make sure you provide the correct sign to the flux
        real(prec) :: vol1, lambda1, m, t1, t2, theta1, dTdt
        dTdt = m / vol1 * (t2 + t1) + lambda1 * (theta1 + t1)
        return
    end function calc_temperature_derivative

    function calc_salinity_derivative(vol1, s0, m, s1, s2, f1, f2) result(dSdt)
        ! WARNING: make sure you provide the correct sign to the flux
        real(prec) :: vol1, s0, f1, f2, m, s1, s2, dSdt
        dSdt = m / vol1 * (s2 + s1) + s0 / vol1 * (f1 + f2)
        return
    end function calc_salinity_derivative

    function calc_individual_thermal_coupling_constant(l, c, rho, z) result(lambda)
        real(preci) :: l, c, rho, z, lambda
        lambda = l / (c * rho * z)
        return
    end function calc_individual_thermal_coupling_constant

end module nautilus