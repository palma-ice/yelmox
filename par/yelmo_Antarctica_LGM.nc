&ctrl
    time_init       = 0.0               ! [yr] Starting time 
    time_end        = 100e3              ! [yr] Ending time
    time_equil      = 500.0              ! [yr] Equilibration time 
    time_const      = 1950.0            ! [yr CE] Assumed constant time for non-transient climate (-19050 == 21 kyr ago)
    time_lgm_step   = 15e3              ! [yr] Time to switch time_const => -19050 (21 kyr ago)
    time_pd_step    = 30e3              ! [yr] Time to switch time_const => 1950.0 (PD, 0 kyr ago)
    dtt             = 5.0               ! [yr] Main loop timestep 
    dt1D_out        = 100.0             ! [yr] Frequency of writing 1D output files 
    dt2D_out        = 20e3              ! [yr] Frequency of writing 2D output files
    dt2D_small_out  = 0.0               ! [yr] Frequency of writing 2D output files
    dt_restart      = 20e3              ! [yr] Frequency of writing restart files
    transient_clim  = True              ! Calculate transient climate? 
    use_lgm_step    = False             ! Activate LGM step in forcing at time=time_lgm_step?
    use_pd_step     = False             ! Activate PD step-return in forcing at time=time_pd_step?
    with_ice_sheet  = True              ! Is the ice sheet active?
    equil_method    = "none"             ! "none", "relax", "opt"
/

&opt_L21
    cf_init         = 0.2               ! Initial value of cf_ref everywhere
    cf_min          = 1e-4              ! [--] Minimum allowed cf value
    cf_max          = 1.0               ! [--] Maximum allowed cf value
    tau_c           = 500.0             ! [yr] L21: Optimization relaxation timescale 
    H0              = 100.0             ! [m]  L21: Optimization ice-thickness error scaling 
    
    sigma_err       = 1.5               ! [--] Smoothing radius for error to calculate correction in cf_ref (in multiples of dx)
    sigma_vel       = 200.0             ! [m/a] Speed at which smoothing diminishes to zero
    fill_method     = "cf_min"          ! nearest, analog, cf_min
    
    rel_tau1        = 10.0              ! [yr] Relaxation timescale in time period 1
    rel_tau2        = 1000.0            ! [yr] Relaxation timescale in time period 2
    rel_time1       = 5e3               ! [yr] Time limit for time period 1
    rel_time2       = 10e3              ! [yr] Time limit for time period 2
    rel_m           = 2.0               ! [--] Non-linear exponent to scale interpolation of rel_tau between time1 and time2 

    opt_tf          = False             ! Should thermal forcing be optimized too?
    H_grnd_lim      = 500.0             ! [m] Distance from grounding line to consider in terms of thickness above flotation
    tau_m           = 100.0             ! [yr] Optimization relxation timescale 
    m_temp          = 10.0              ! [m yr^−1 degC^−1] Optimization scaling value
    tf_min          = -2.0              ! [degC] Minimum allowed tf_corr value 
    tf_max          =  5.0              ! [degC] Maximum allowed tf_corr value 
    
/

&yelmo
    domain          = "Antarctica"
    grid_name       = "ANT-16KM"
    grid_path       = "ice_data/{domain}/{grid_name}/{grid_name}_REGIONS.nc"
    experiment      = "None"            ! "None", "EISMINT1", "MISMIP3D"  to apply special boundary conditions
    restart         = "/home/javier/yelmo-ucm/yelmox/restart_16km_pd/yelmo_restart.nc"
    restart_z_bed   = True              ! Take z_bed (and z_bed_sd) from restart file
    restart_H_ice   = True              ! Take H_ice from restart file
    log_timestep    = False 
    disable_kill    = False             ! Disable automatic kill if unstable
    zeta_scale      = "exp"             ! "linear", "exp", "tanh"
    zeta_exp        = 2.0  
    nz_aa           = 10                ! Vertical resolution in ice
    dt_method       = 2                 ! 0: no internal timestep, 1: adaptive, cfl, 2: adaptive, pc
    dt_min          = 0.01              ! [a] Minimum timestep 
    cfl_max         = 0.5               ! Maximum value is 1.0, lower will be more stable
    cfl_diff_max    = 0.12              ! Bueler et al (2007), Eq. 25  
    pc_method       = "AB-SAM"          ! "FE-SBE", "AB-SAM", "HEUN"
    pc_controller   = "PI42"            ! PI42, H312b, H312PID, H321PID, PID1
    pc_use_H_pred   = True              ! Use predicted H_ice instead of corrected H_ice 
    pc_filter_vel   = True              ! Use mean vel. solution of current and previous timestep
    pc_corr_vel     = False             ! Calculate dynamics again for corrected H_ice
    pc_n_redo       = 5                 ! How many times can the same iteration be repeated (when high error exists)
    pc_tol          = 10.0              ! [m/a] Tolerance threshold to redo timestep
    pc_eps          = 3.0               ! Predictor-corrector tolerance 
    
/

&ytopo
    solver              = "impl-upwind"     ! "expl", "impl-upwind"
    calv_flt_method     = "vm-l19"          ! "zero", "simple", "flux", "vm-l19", "kill"
    calv_grnd_method    = "stress-b12"      ! "zero", "stress-b12"
    bmb_gl_method       = "pmp"             ! "fcmp": flotation; "fmp": full melt; "pmp": partial melt; "nmp": no melt
    fmb_method          = 1                 ! 0: prescribe boundary field fmb_shlf; 1: calculate proportional fmb~bmb_shlf.
    surf_gl_method      = 0                 ! 0: binary (max grnd/flt elevation), 1: subgrid average elevation
    margin_flt_subgrid  = True              ! Allow fractional ice area for floating margins
    margin2nd           = False             ! Apply second-order upwind approximation to gradients at the margin
    use_bmb             = True              ! Use basal mass balance in mass conservation equation
    topo_fixed          = False             ! Keep ice thickness fixed, perform other ytopo calculations
    topo_rel            = 0                 ! 0: No relaxation; 1: relax shelf; 2: relax shelf + gl; 3: all points 
    topo_rel_tau        = 10.0              ! [a] Time scale for relaxation 
    topo_rel_field      = "H_ref"           ! "H_ref" or "H_ice_n"
    calv_H_lim          = 200.0             ! [m] Calving limit in ice thickness (thinner ice calves)
    calv_tau            = 1.0               ! [a] Characteristic calving time
    calv_thin           = 30.0              ! [m/yr] Calving rate to apply to very thin ice
    H_min_grnd          = 5.0               ! [m] Minimum ice thickness at grounded margin (thinner ice is ablated) - helps with stability
    H_min_flt           = 0.0               ! [m] Minimum ice thickness at floating margin (thinner ice is ablated) - helps with stability
    sd_min              = 100.0             ! [m] calv_grnd(z_bed_sd <  = sd_min)     = 0.0 
    sd_max              = 500.0             ! [m] calv_grnd(z_bed_sd >  = sd_max)     = calv_max  
    calv_max            = 10.0              ! [m/a] Maximum grounded calving rate from high stdev(z_bed)
    grad_lim            = 0.1               ! [m/m] Maximum allowed sloped in gradient calculations (dz/dx,dH/dx)
    dist_grz            = 200.0             ! [km] Radius to consider part of "grounding-line zone" (grz)
    gl_sep              = 1                 ! 1: Linear f_grnd_acx/acy and binary f_grnd, 2: area f_grnd, average to acx/acy
    gl_sep_nx           = 15                ! [-] Number of interpolation points (nx*nx) to calculate grounded area at grounding line
    diffuse_bmb_shlf    = False             ! Allow bmb_shlf to permeate inland at the grounding line 
    fmb_scale           = 1.0               ! Scaling of fmb ~ scale*bmb, scale=10 suggested by Pollard and DeConto (2016)
    kt                  = 0.0025            ! [m yr-1 Pa-1] vm-l19 calving scaling parameter (L21 use 0.0025)
    w2                  = 25                ! [-] vm-l19 calving eigenvalue weighting coefficient
    k2                  = 3e9               ! [m yr] eigen calving scaling factor (Albrecht et al, 2021 recommend 1e17 m s == 3.2e9 m yr)
/

&ydyn
    solver              = "diva"          ! "fixed", "sia", "ssa", "hybrid", "diva", "diva-noslip", "l1l2", "l1l2-noslip"
    visc_method         = 1               ! 0: constant visc=visc_const, 1: dynamic viscosity
    visc_const          = 1e7             ! [Pa a] Constant value for viscosity (if visc_method=0)
    beta_method         = 3               ! 0: constant beta; 1: linear (beta=cb/u0); 2: psuedo-plastic-power; 3: Regularized Coulomb
    beta_const          = 1e3             ! [Pa a m−1] Constant value of basal friction coefficient to be used
    beta_q              = 0.2             ! Dragging law exponent 
    beta_u0             = 100.0           ! [m/a] Regularization term for regularized Coulomb law (beta_method=3)
    beta_gl_scale       = 2               !  0: beta*beta_gl_f, 1: H_grnd linear scaling, 2: Zstar scaling 
    beta_gl_stag        = 3               !  0: simple staggering, 1: Upstream beta at gl, 2: downstream, 3: f_grnd_ac scaling, 4: vel weighting 
    beta_gl_f           = 1.0             ! [-] Scaling of beta at the grounding line (for beta_gl_scale=0)
    taud_gl_method      = 0               !  0: binary, no subgrid, 1: Two-sided gradient
    H_grnd_lim          = 500.0           ! [m] For beta_gl_scale=1, reduce beta linearly between H_grnd=0 and H_grnd_lim 
    n_sm_beta           = 0               ! [-] Standard deviation in gridpoints for Gaussian smoothing of beta (0==no smoothing)
    beta_min            = 100.0           ! [Pa a m-1] Minimum value of beta allowed for grounded ice (for stability)
    eps_0               = 1e-6            ! [1/a] Regularization term for effective viscosity - minimum strain rate
    ssa_lis_opt         = "-i minres -p jacobi -maxiter 100 -tol 1.0e-2 -initx_zeros false"  ! See Lis library !-omp_num_threads 2
    ssa_lat_bc          = "floating"      ! "all", "marine", "floating", "none", "slab"
    ssa_beta_max        = 1e20            ! [Pa a m-1] Maximum value of beta for which ssa should be calculated 
    ssa_vel_max         = 5000.0          ! [m a-1] SSA velocity limit to avoid spurious results 
    ssa_iter_max        = 5               ! Number of maximum allowed iterations over ssa to converge on vel. solution
    ssa_iter_rel        = 0.7             ! [--] Relaxation fraction [0:1] to stabilize ssa iterations
    ssa_iter_conv       = 1e-2            ! [--] L2 relative error convergence limit to exit ssa iterations
    taud_lim            = 2e5             ! [Pa] Maximum allowed driving stress 
    cb_sia              = 0.0             ! [m a-1 Pa-1] Bed roughness coefficient for SIA sliding
    
/

&ytill
    method          = 1                 ! -1: set externally; 1: calculate cb_ref online  
    scale           = "exp"             ! "none", "lin", or "exp" : scaling with elevation 
    is_angle        = False             ! cb_ref is a till strength angle?
    z0              = -400.0            ! [m] Bedrock rel. to sea level, lower limit
    z1              =    0.0            ! [m] Bedrock rel. to sea level, upper limit
    cf_min          =  0.001            ! [-- or deg] Minimum value of cf
    cf_ref          =  0.25             ! [-- or deg] Reference/const/max value of cf
/

&yneff 
    method          = 2                 ! -1: external N_eff, 0: neff_const, 1: overburden pressure, 2: Leguy param., 3: Till pressure
    const           = 1e7               ! == rho_ice*g*(eg 1000 m ice thickness)
    p               = 1.0               ! *neff_method=2* marine connectivity exponent (0: none, 1: full)
    set_water       = False             ! *neff_method=3* Prescribe H_w = H_w_max for temperate ice instead of using H_w field?
    N0              = 1000.0            ! [Pa] *neff_method=3* Reference effective pressure 
    delta           = 0.04              ! [--] *neff_method=3* Fraction of overburden pressure for saturated till
    e0              = 0.69              ! [--] *neff_method=3* Reference void ratio at N0 
    Cc              = 0.12              ! [--] *neff_method=3* Till compressibility    
/

&ymat
    flow_law                = "glen"        ! Only "glen" is possible right now
    rf_method               = 1             ! -1: set externally; 0: rf_const everywhere; 1: standard function 
    rf_const                = 1e-18         ! [Pa^-3 a^-1]
    rf_use_eismint2         = False         ! Only applied for rf_method=1
    rf_with_water           = False         ! Only applied for rf_method=1, scale rf by water content?
    n_glen                  = 3.0           ! Glen flow law exponent
    visc_min                = 1e3           ! [Pa a] Minimum allowed viscosity 
    de_max                  = 0.5           ! [a-1]  Maximum allowed effective strain rate
    enh_method              = "shear3D"     ! "simple","shear2D", "shear3D", "paleo-shear" 
    enh_shear               = 4.0
    enh_stream              = 1.0
    enh_shlf                = 0.7
    enh_umin                = 50.0          ! [m/yr] Minimum transition velocity to enh_stream (enh_method="paleo-shear")
    enh_umax                = 500.0         ! [m/yr] Maximum transition velocity to enh_stream (enh_method="paleo-shear")
    calc_age                = False         ! Calculate age tracer field?
    age_iso                 = 0.0
    tracer_method           = "expl"        ! "expl", "impl": used for age and 'paleo-shear' enh fields
    tracer_impl_kappa       = 1.5           ! [m2 a-1] Artificial diffusion term for implicit tracer solving 
    
/

&ytherm
    method          = "temp"            ! "fixed", "robin", "temp", "enth"
    dt_method       = "FE"              ! "FE", "AB", "SAM"
    solver_advec    = "impl-upwind"     ! "expl", "impl-upwind"
    gamma           = 1.0               ! [K] Scalar for the pressure melting point decay function 
    use_strain_sia  = False             ! True: calculate strain heating from SIA approx.
    n_sm_qstrn      = 0                 ! [-] Standard deviation in gridpoints for Gaussian smoothing of strain heating (0==no smoothing)
    n_sm_qb         = 0                 ! [-] Standard deviation in gridpoints for Gaussian smoothing of basal heating (0==no smoothing)
    use_const_cp    = False             ! Use specified constant value of heat capacity?
    const_cp        = 2009.0            ! [J kg-1 K-1] Specific heat capacity 
    use_const_kt    = False             ! Use specified constant value of heat conductivity?
    const_kt        = 6.62e7            ! [J a-1 m-1 K-1] Thermal conductivity [W m-1 K-1 * sec_year] => [J a-1 m-1 K-1]
    enth_cr         = 1e-3              ! [--] Conductivity ratio for temperate ice (kappa_temp     = enth_cr*kappa_cold)
    omega_max       = 0.01              ! [--] Maximum allowed water content fraction 
    till_rate       = 1e-3              ! [m/a] Basal water till drainage rate (water equiv.)
    H_w_max         = 2.0               ! [m] Maximum limit to basal water layer thickness (water equiv.)

    rock_method     = "active"          ! "equil" (not active bedrock) or "active"
    nzr_aa          = 5                 ! Number of vertical points in bedrock 
    zeta_scale_rock = "exp-inv"         ! "linear", "exp-inv"
    zeta_exp_rock   = 2.0  
    H_rock          = 2000.0            ! [m] Lithosphere thickness 
    cp_rock         = 1000.0            ! [J kg-1 K-1] Specific heat capacity of bedrock
    kt_rock         = 6.3e7             ! [J a-1 m-1 K-1] Thermal conductivity of bedrock [W m-1 K-1 * sec_year] => [J a-1 m-1 K-1]
/

&yelmo_masks
    basins_load     = True 
    basins_path     = "ice_data/{domain}/{grid_name}/{grid_name}_BASINS-nasa.nc" 
    basins_nms      = "basin_reese" "basin_mask"
    regions_load    = True 
    regions_path    = "ice_data/{domain}/{grid_name}/{grid_name}_BASINS-nasa.nc"
    regions_nms     = "mask_regions" "None"
/

&yelmo_init_topo
    init_topo_load    = True 
    init_topo_path    = "ice_data/{domain}/{grid_name}/{grid_name}_TOPO-BedMachine.nc"
    init_topo_names   = "H_ice" "z_bed" "z_bed_sd" "z_srf"      ! Ice thickness, bedrock elevation, bedrock noise, surface elevation
    init_topo_state   = 0                                       ! 0: from file, 1: ice-free, 2: ice-free, rebounded 
    z_bed_f_sd        = 0.0                                     ! Scaling fraction to modify z_bed = z_bed + f_sd*z_bed_sd
    smooth_H_ice      = 0.0                                     ! Smooth ice thickness field at loading time, with sigma=N*dx
    smooth_z_bed      = 0.0                                     ! Smooth bedrock field at loading time, with sigma=N*dx
    
/

&yelmo_data 
    pd_topo_load      = True 
    pd_topo_path      = "ice_data/{domain}/{grid_name}/{grid_name}_TOPO-BedMachine.nc"
    pd_topo_names     = "H_ice" "z_bed" "z_bed_sd" "z_srf"  ! Ice thickness, bedrock elevation, bedrock noise, surface elevation
    pd_tsrf_load      = True 
    pd_tsrf_path      = "ice_data/{domain}/{grid_name}/{grid_name}_RACMO23-ERAINT-HYBRID_1981-2010.nc"
    pd_tsrf_name      = "T_srf"                      ! Surface temperature (or near-surface temperature)
    pd_tsrf_monthly   = True
    pd_smb_load       = True 
    pd_smb_path       = "ice_data/{domain}/{grid_name}/{grid_name}_RACMO23-ERA-INTERIM_monthly_1981-2010.nc"
    pd_smb_name       = "smb"                        ! Surface mass balance 
    pd_smb_monthly    = True 
    pd_vel_load       = True 
    pd_vel_path       = "ice_data/{domain}/{grid_name}/{grid_name}_VEL-R11-2.nc"
    pd_vel_names      = "ux_srf" "uy_srf"            ! Surface velocity 
    pd_age_load       = False 
    pd_age_path       = "input/{domain}/{grid_name}/{grid_name}_STRAT-M15.nc"
    pd_age_names      = "age_iso" "depth_iso"        ! ages of isochrones; depth of isochrones 

/

! ===== EXTERNAL LIBRARIES =====

&sealevel
    method              = 0                             ! 0: use z_sl_const, 1: read time series from sl_path
    z_sl_const          = -130.0                          ! [m] Sea-level relative to present day
    sl_path             = "input/sl_deglaciation.dat"   ! Input time series filename
    sl_name             = "none"                        ! Variable name in netcdf file, if relevant
/

&snap
    atm_type           = "anom"
    ocn_type           = "anom"
    fname_at           = "input/ones.dat"
    fname_ao           = "input/ones.dat"
    fname_ap           = "input/ones.dat"
    fname_as           = "input/ones.dat"
    fname_bt           = "input/ones.dat"
    fname_bo           = "input/ones.dat"
    fname_bp           = "input/ones.dat"
    fname_bs           = "input/ones.dat"
    lapse              = 0.0080 0.0065    ! lapse_ann, lapse_sum 
    dTa_const          = -10.0  
    dTo_const          = -2.5
    dSo_const          = 3.0 
    f_to               = 0.25 
    f_p                = 0.05
    f_stdev            = 0.0  
/

&snap_hybrid
    hybrid_path        = "ice_data/Greenland/paleo300ka_hybrid_aicc2012.nc"
    f_eem              = 1.0 
    f_glac             = 1.0 
    f_hol              = 1.0 
    f_seas             = 1.0 
    f_to               = 0.2
/

&snap_clim0
    clim_path       = "ice_data/{domain}/{grid_name}/{grid_name}_RACMO23-ERAINT-HYBRID_1981-2010.nc"
    clim_names      = "zs" "T_srf" "pr" ""                                                   ! Default (if monthly = True) : "zs" "T_srf" "pr" "" (if False) "z_srf" "t2m_ann" "t2m_djf" "pr_ann"
    clim_monthly    = True                                                                   ! If True and no snowfall (sf) then the third variable is precip and the fourth one irrelevant
    clim_time       = 0
    clim_stdev_path = "None"
    clim_stdev_name = "pr" 
    ocn_path        = "ice_data/{domain}/{grid_name}/{grid_name}_OCEAN_ISMIP6_J20.nc"
    ocn_names       = "z" "mask_ocn" "to" "so"
    ocn_monthly     = False
    ocn_time        = 0
/

&snap_clim1
    clim_path       = "ice_data/PlioMIP2/{domain}/{grid_name}/Atmosphere/{grid_name}_IPSLCM5A2_E280.nc"
    clim_names      = "z_srf" "t2m" "pr" ""
    clim_monthly    = True
    clim_time       = 0
    clim_stdev_path = "None"
    clim_stdev_name = "pr" 
    ocn_path        = "none"
    ocn_names       = "depth" "mask_ocn" "to" "so"
    ocn_monthly     = False
    ocn_time        = 0
/

&snap_clim2
    clim_path       = "ice_data/PlioMIP2/{domain}/{grid_name}/Atmosphere/{grid_name}_IPSLCM5A2_Eoi400.nc"
    clim_names      = "z_srf" "t2m" "pr" ""
    clim_monthly    =  True
    clim_time       = -21000
    clim_stdev_path = "None"
    clim_stdev_name = "pr"
    ocn_path        = "none"
    ocn_names       = "depth" "mask_ocn" "to" "so"
    ocn_monthly     = False
    ocn_time        = -21000
/

&snap_clim3
    clim_path       = "ice_data/{domain}/{grid_name}/Montoya2008_sig250km/{grid_name}_lgm_1p7strong.nc"
    clim_names      = "zs" "t2m_ann" "t2m_sum" "pr_ann"
    clim_monthly    = False
    clim_time       = -21000
    clim_stdev_path = "None"
    clim_stdev_name = "pr" 
    ocn_path        = "ice_data/{domain}/{grid_name}/Montoya2008_sig100km/{grid_name}_lgm_1p7strong_ocean.nc"
    ocn_names       = "depth" "mask_ocn" "to"
    ocn_monthly     = False
    ocn_time        = -21000
/

&smbpal
    const_insol = True
    const_kabp  = 0.0
    insol_fldr  = "input"
    abl_method  = "itm"   ! "itm" or "pdd"
    sigma_snow  = 5.0     ! Standard deviation of temperature [K] over ice and snow
    sigma_melt  = 5.0     ! Standard deviation of temperature [K] in ablation zone
    sigma_land  = 5.0     ! Standard deviation of temperature [K] over land 
    
    sf_a        = 0.273   ! Snow fraction function (parameter a)
    sf_b        = 273.6   ! Snow fraction function (parameter b)
    firn_fac    = 0.0266  ! Reduction in surface temp versus melting [K/mm w.e.]

    mm_snow = 3.0         ! PDD snow melt factor [mm w.e./K-d]
    mm_ice  = 8.0         ! PDD ice melt factor  [mm w.e./K-d]
    
/

&itm
        
    H_snow_max      =  5.0e3         ! Maximum allowed snowpack thickness [mm w.e.]
    Pmaxfrac        =  0.6           ! Refreezing factor 

    ! ## ITM 
    trans_a         =  0.46 
    trans_b         =  6e-5
    trans_c         =  0.01
    itm_c           = -55.0
    itm_t           =  10.0
    itm_b           =  3.0 
    itm_lat0        =  -60.0 

    ! ## Surface albedo
    alb_snow_dry       =  0.8
    alb_snow_wet       =  0.65
    alb_ice            =  0.4
    alb_land           =  0.2
    alb_forest         =  0.1
    alb_ocean          =  0.1
    H_snow_crit_desert =  10.0
    H_snow_crit_forest = 100.0
    melt_crit          =   0.5      ! Critical melt rate to reduce snow albedo [mm/day]

/

&marine_shelf
    bmb_method      = "quad-nl"     ! "pico", "anom", "lin", "quad", "quad-nl", "*-slope"
    tf_method       = 1             ! 0: tf defined externally, 1: calculate as T_shlf-t_fp_shlf, 2: use dT_shlf
    interp_depth    = "shlf"        ! "shlf", "bed", "const"
    interp_method   = "interp"      ! "mean", "layer", "interp" ! Same method is used for salinity
    find_ocean      = True          ! Find which points are connected to open ocean (ie, ignore subglacial lakes)
    
    corr_method     = "tf-ant"      ! "tf": corrected tf basin, "bmb": corrected bmb basin, "tf-ant": use values in tf_corr_ant parameter section
    basin_number    = 1             ! 1 2 3
    basin_bmb_corr  = 0.0           ! -1.0  -1.0 -1.0 [corrected bmb m/yr; -: melting, +: accretion]
    basin_tf_corr   = 0.0           ! 1.0   2.0  3.0  
    tf_correction   = True
    tf_path         = "ice_data/{domain}/{grid_name}/{grid_name}_OCEAN_ISMIP6_J20.nc"
    tf_name         = "dT_nl"
 
    bmb_max         = 0.0           ! [m/a] Maximum allowed bmb value (eg, for no refreezing, set bmb_max=0.0) 
    
    ! Deep ocean bmb 
    c_deep          = -1.0          ! [m/a] Melt rate in deep ocean, should be negative but not huge
    depth_deep      = 2000.0        ! [m] Depth at which deep ocean is assumed
    depth_const     = 2000.0        ! [m] Constant shelf depth for interp_depth=const
    depth_min       = 1100.0        ! [m] interp_method=mean, min depth to include
    depth_max       = 1500.0        ! [m] interp_method=mean, max depth to include
    
    ! Freezing point coefficients
    lambda1         = -0.0575       ! [degC PSU−1] Liquidus slope
    lambda2         = 0.0832        ! [degC] Liquidus intercept 
    lambda3         = 7.59e-4       ! [degC m−1] Liquidus pressure coefficient

    ! bmb_method == [lin,quad,quad-nl] 
    gamma_lin       = 630.0         ! [m/yr] Linear heat excange velocity. Favier et al., (2019) suggest 630.
    gamma_quad      = 11100.0       ! [m/yr] Quadratic heat exchange velocity. Jourdain et al. (2020) suggest 11100.
    gamma_quad_nl   = 14500.0       ! [m/yr] Quadratic non-local heat exchange velocity. Jourdain et al. (2020) suggest 14500.
    gamma_prime     = 100.0         ! [-] Scalar when '*-slope' method is used

    ! bmb_method == anom
    kappa_grz       = 10.0          ! [m/a K-1] Heat flux coefficient [kappa*dT]
    c_grz           = 1.0           ! [-] bmb_ref scalar [c_grz*bmb_ref]
    f_grz_shlf      = 10.0          ! [-] Ratio of grz melt to shelf melt
    grz_length      = 16.0          ! [km] Length scale of assumed grounding zone

    ! bmb_method == anom, bmb_ref parameters 
    use_obs         = False          ! Load obs for bmb_ref
    obs_path        = "ice_data/{domain}/{grid_name}/{grid_name}_BMELT-R13.nc"
    obs_name        = "bm_eq_reese"
    obs_scale       = 1.0           ! [-] bmb_ref = obs_scale*bmb_obs
    obs_lim         = 100.0         ! [m/yr] Maximum magnitude of bmb_ref
    basin_name      = "none"        !"northwest" "northeast" "east"
    basin_bmb       = 0.0           !        -1.0       -1.0   -1.0
    
/

&tf_corr_ant
    ronne   = 0.2                   ! [K] tf_corr for ronne basin (basin=1)
    ross    = 0.3                   ! [K] tf_corr for ronne ross  (basin=12)
    pine    = -1.6                   ! [K] tf_corr for ronne pine  (basin=14)
    abbott  = 0.0                   ! [K] tf_corr for Abbott glacier (basin=15)
/

&pico
    n_box          = 5.0        ! Maximum number of boxes in PICO.  
    a_pico         = -0.0572    ! Salinity coefficient of freezing equation [ºC/PSU]. Reese et al., (2018) suggest −0.0572.
    b_pico         = 0.0788     ! Constant coefficient of freezing equation [ºC]. Reese et al., (2018) suggest 0.0788.
    c_pico         = 7.77e-8    ! Pressure coefficient of freezing equation [ºC/Pa]. Reese et al., (2018) suggest 7.77e-8.
    alpha_pico     = 7.5e-5     ! Thermal expansion coefficien [1/ºC]. Reese et al., (2018) suggest 7.5e-5.
    beta_pico      = 7.7e-4     ! Salt contraction coefficient [1/PSU]. Reese et al., (2018) suggest 7.7e-4.
    rho_star       = 1033.0     ! Reference density [kg/m3]. Reese et al., (2018) suggest 1033.0.
    gamma_tstar    = 630.0      ! Effective turbulent temperature exchange velocity [m/yr]. Reese et al., (2018) uses 630.0 and suggest 600-1000.
    C_over         = 1.0        ! Overturning strength [Sv m3/kg]. Reese et al., (2018) suggest 1.0.

/

&sediments
    use_obs   = False
    obs_path  = "ice_data/{domain}/{grid_name}/{grid_name}_SED-L97.nc"
    obs_name  = "z_sed"

/

&geothermal
    use_obs   = True 
    obs_path  = "ice_data/{domain}/{grid_name}/{grid_name}_GHF-D13.nc"
    obs_name  = "ghf_mean"
    ghf_const = 50.0           ! [mW/m^2]

/

&isostasy 
    method       = 2                    ! 0: constant lithosphere depression; 1: LLRA (not tested); 2: ELRA
    fname_kelvin = "input/kelvin.res"   ! File containing precalculated zero-order Kelvin function values 
    dt           = 50.0                 ! [yr] Time step to recalculate bedrock uplift rate 
    tau          = 3000.0               ! [yr] Relaxation time 
    D_lith       = 9.87e24              ! [N-m] Flexural rigidity 
    R_lith       = 131910.0             ! [m] Radius of relative stiffness
    
/
