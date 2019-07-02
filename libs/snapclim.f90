

module snapclim 

    !use yelmo_defs !,only :: sp, dp, prec, pi, yelmo_parse_path
    use ncio 
    use nml 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: prec = sp 

    real(prec), parameter :: sec_year  = 365.0*24.0*60.0*60.0   ! [s/a]
    real(prec), parameter :: pi        = 3.14159265359

    type series_type
        character(len=512) :: filename 
        real(prec), allocatable :: time(:), var(:), sigma(:)
    end type 

    type series_2D_type
        character(len=512) :: filename 
        real(prec), allocatable :: time(:), var(:,:), sigma(:,:)
    end type     

    type snapshot_param_class 
        character(len=512) :: clim_path 
        character(len=56)  :: clim_names(4)
        logical            :: clim_monthly 
        real(prec)         :: clim_year_bp 

        character(len=512) :: ocn_path 
        character(len=56)  :: ocn_names(3)
        logical            :: ocn_monthly 
        real(prec)         :: ocn_year_bp 
        
    end type 

    type hybrid_class 
        character(len=512)   :: hybrid_path 
        real(prec)           :: f_eem, f_glac, f_hol, f_seas, f_to
        type(series_2D_type) :: dTmon
    
    end type 

    type snapclim_param_class
        integer    :: nx, ny
        character(len=56)  :: clim_type
        character(len=512) :: fname_at, fname_ap, fname_ao 
        character(len=512) :: fname_bt, fname_bp, fname_bo   
        logical    :: spatial_const_clim
        real(prec) :: ta_clim1, ta_clim2, ta_clim3
        logical    :: spatial_const_ocn 
        real(prec) :: to_clim1, to_clim2, to_clim3
        logical    :: artificial_acc
        real(prec) :: acc_ross, acc_ronne
        real(prec) :: lapse(2)
        real(prec) :: dTa_const 
        real(prec) :: dTo_const 
        real(prec) :: f_p

    end type

    type snapclim_state_class 

        type(snapshot_param_class) :: par 

        ! Climate variables 
        real(prec), allocatable :: mask(:,:)     
        real(prec), allocatable :: z_srf(:,:)

        real(prec), allocatable :: tas(:,:,:)
        real(prec), allocatable :: pr(:,:,:)
        real(prec), allocatable :: sf(:,:,:)

        real(prec), allocatable :: ta_ann(:,:)
        real(prec), allocatable :: ta_sum(:,:)
        real(prec), allocatable :: pr_ann(:,:)
        
        real(prec), allocatable :: tsl(:,:,:)
        real(prec), allocatable :: prcor(:,:,:)
        real(prec), allocatable :: tsl_ann(:,:)
        real(prec), allocatable :: tsl_sum(:,:)
        real(prec), allocatable :: prcor_ann(:,:)
        
        ! Oceanic variables
        integer :: nzo
        real(prec), allocatable :: depth(:) 
        real(prec), allocatable :: mask_ocn(:,:,:) 
        real(prec), allocatable :: to_ann(:,:,:) 
        
        ! Indices 
        real(prec) :: at, ap, ao, bt, bp, bo

    end type 

    type snapclim_class

        type(snapclim_param_class) :: par 

        type(snapclim_state_class) :: now  
        type(snapclim_state_class) :: clim0, clim1, clim2, clim3  
        type(series_type)          :: at, ap, ao 
        type(series_type)          :: bt, bp, bo  
        
        type(hybrid_class)         :: hybrid 
    end type 

    private
    public :: snapclim_class 
    public :: snapclim_init 
    public :: snapclim_update 
!     public :: snapclim_end 

    ! === NOTES ===

    ! clim_type indicates the type of climate forcing aplied
    ! =1 | No climate_forcing per se ==> only climatologic fields (or any other reference climate)
    ! =2 | Anomalies methode by default (the anomaly from two climate fields interpolated by 1 index+ climatologies)
    ! =3 | Anomalies methode two index (two anomalies from three climate fields interpolated by 2 index+ climatologies)
    ! =4 | Absolute method 1 index (two absolute climates weighted in time by 1 index | no climatologies)
    ! =5 | Absolute method 2 index (three absolute climates weighted in time by 2 index | no climatologies) 
    !    | More cases TO BE CONSIDERED (theoretical domains...)

    ! climX_monthly
    ! =False | The climate_forcing files ONLY provide Tann and Tjja and annual Precip
    ! =True  | The climate_forcing files provide monthly temperatures and precipitations

contains 

    subroutine snapclim_init(snp,filename,domain,grid_name,nx,ny)
        ! This subroutine will initialize four climate snapshots
        ! (clim0,clim1,clim2,clim3) which will be used for temporal
        ! interpolation to determine the current climate forcing. 
        ! Note: snapclim_update should be called after initialiation

        implicit none 

        type(snapclim_class), intent(INOUT) :: snp 
        character(len=*),     intent(IN)    :: filename 
        character(len=*),     intent(IN)    :: domain, grid_name
        integer,    intent(IN) :: nx, ny  

        ! Local variables 
        logical :: load_clim1, load_clim2, load_clim3 
        integer :: k, nzo  
        real(prec), allocatable :: depth(:) 

        ! Load parameters 
        call snapclim_par_load(snp%par,snp%hybrid,filename)
        snp%par%nx = nx 
        snp%par%ny = ny 

        ! Determine which snapshots should be loaded
        select case(trim(snp%par%clim_type))

            case("const","hybrid","anom")

                load_clim1 = .FALSE.
                load_clim2 = .FALSE. 
                load_clim3 = .FALSE. 

            case("anom_1ind","abs_1nd")
                load_clim1 = .TRUE.
                load_clim2 = .TRUE. 
                load_clim3 = .FALSE. 

            case DEFAULT 
                load_clim1 = .TRUE.
                load_clim2 = .TRUE. 
                load_clim3 = .TRUE. 
 
        end select 

        ! Initialize output ocean elevation (ie, -depth) and ocean temp arrays
        nzo = 23 
        allocate(depth(nzo))

        ! Intialize ocean depths between 0 and 3000 m 
        do k = 1, nzo-2 
            depth(k) = (k-1)*2000.0/(nzo-3)
        end do 
        depth(nzo-1:nzo) = [2500.0,3000.0]

        write(*,"(a)") "snapclim_init:: ocean depths initialized."
        write(*,"(a,i5)")      "nzo =   ", nzo 
        write(*,"(a,50f10.1)") "depth = ", depth 

        ! ================================================================
        !
        ! Step 2: Read in all necessary input climate data from
        !         snapshots and time series files 
        !
        ! ================================================================

        ! == clim0: reference climate (eg, present day) ==

        call snapshot_par_load(snp%clim0%par,filename,"snapclim_clim0",domain,grid_name,init=.TRUE.)
        call read_climate_snapshot(snp%clim0,nx,ny,snp%par%lapse,snp%par%f_p,domain)
        call read_ocean_snapshot(snp%clim0,nx,ny,depth=depth)
            
        if (load_clim1) then
            ! == clim1: snapshot 1 (eg, present day from model) == 

            call snapshot_par_load(snp%clim1%par,filename,"snapclim_clim1",domain,grid_name,init=.TRUE.)
                
            if (snp%par%spatial_const_clim) then
                call set_climate_snapshot(snp%clim1,snp%clim0,nx,ny,year_bp=0.0,tas=snp%par%ta_clim1,f_p=snp%par%f_p)
            else 
                call read_climate_snapshot(snp%clim1,nx,ny,snp%par%lapse,snp%par%f_p,domain)
            end if 

            if (snp%par%spatial_const_ocn) then
                call set_ocn_snapshot(snp%clim1,nx,ny,depth=depth,year_bp=0.0,t_ocn=snp%par%to_clim1)
            else 
                call read_ocean_snapshot(snp%clim1,nx,ny,depth=depth)
            end if 
        end if 

        if (load_clim2) then
            ! == clim2: snapshot 2 (eg, LGM with strong AMOC) == 

            call snapshot_par_load(snp%clim2%par,filename,"snapclim_clim2",domain,grid_name,init=.TRUE.)
                
            if (snp%par%spatial_const_clim) then
                call set_climate_snapshot(snp%clim2,snp%clim0,nx,ny,year_bp=0.0,tas=snp%par%ta_clim2,f_p=snp%par%f_p)
            else 
                call read_climate_snapshot(snp%clim2,nx,ny,snp%par%lapse,snp%par%f_p,domain)
            end if 

            if (snp%par%spatial_const_ocn) then
                call set_ocn_snapshot(snp%clim2,nx,ny,depth=depth,year_bp=0.0,t_ocn=snp%par%to_clim2)
            else 
                call read_ocean_snapshot(snp%clim2,nx,ny,depth=depth)
            end if 

        end if 

        if (load_clim3) then
            ! == clim3: snapshot 3 (eg, LGM with weak AMOC) == 

            call snapshot_par_load(snp%clim3%par,filename,"snapclim_clim3",domain,grid_name,init=.TRUE.)
                
            if (snp%par%spatial_const_clim) then
                call set_climate_snapshot(snp%clim3,snp%clim0,nx,ny,year_bp=0.0,tas=snp%par%ta_clim3,f_p=snp%par%f_p)
            else 
                call read_climate_snapshot(snp%clim3,nx,ny,snp%par%lapse,snp%par%f_p,domain)
            end if 

            if (snp%par%spatial_const_ocn) then
                call set_ocn_snapshot(snp%clim3,nx,ny,depth=depth,year_bp=0.0,t_ocn=snp%par%to_clim3)
            else 
                call read_ocean_snapshot(snp%clim3,nx,ny,depth=depth)
            end if 

        end if 


        ! In case we are using the hybrid monthly anomaly time series to generate climate forcing
        if (trim(snp%par%clim_type) .eq. "hybrid") then  
            call read_series_hybrid(snp%hybrid%dTmon,snp%hybrid%hybrid_path, &
                        snp%hybrid%f_eem,snp%hybrid%f_glac,snp%hybrid%f_hol,snp%hybrid%f_seas)
        end if 

        ! Read in forcing time series
        ! Make sure the time series files do not conatin empty lines at the end
        call read_series(snp%at, snp%par%fname_at)
        call read_series(snp%ap, snp%par%fname_ap)
        call read_series(snp%ao, snp%par%fname_ao)
        call read_series(snp%bt, snp%par%fname_bt)
        call read_series(snp%bp, snp%par%fname_bp)
        call read_series(snp%bo, snp%par%fname_bo)
   
        ! Also check consistency that in the case of absolute forcing clim0 = clim1,
        ! which is the reference climate 
        if (trim(snp%par%clim_type) .eq. "abs_1ind" .or. trim(snp%par%clim_type) .eq. "abs_2ind") then 
            write(*,*) "snapclim_init:: overriding user choice for clim0, now clim0=clim1 for consistency."
            snp%clim0 = snp%clim1 
        end if 

        ! Initialize the current climate state to the reference climate
        ! (this also allocates the now variables).
        snp%now = snp%clim0 

        return 

    end subroutine snapclim_init

    subroutine snapclim_update(snp,z_srf,year_bp,domain,dTa,dTo)

        implicit none 

        type(snapclim_class), intent(INOUT) :: snp
        real(prec), intent(IN)    :: z_srf(:,:) 
        real(prec), intent(IN)    :: year_bp    ! Current simulation year
        character(len=*), intent(IN) :: domain 
        real(prec), intent(IN), optional :: dTa   ! For clim_type='anom'
        real(prec), intent(IN), optional :: dTo   ! For clim_type='anom'
        
        ! Local variables
        real(prec) :: at, ap, ao, bt, bp, bo
        real(prec) :: dTa_now, dTo_now  
        real(prec) :: dT(12) 
        logical :: south 
        integer :: m 
        integer :: nx, ny, i, j

        ! Determine the current values of various indices 
        at     = series_interp(snp%at, year_bp)
        ap     = series_interp(snp%ap, year_bp)
        ao     = series_interp(snp%ao, year_bp)
        bt     = series_interp(snp%bt, year_bp)
        bp     = series_interp(snp%bp, year_bp)
        bo     = series_interp(snp%bo, year_bp)

        !write(*,"(6f12.2)") year_bp, at, ap, ao, bt, bp, bo

        ! Determine whether the domain is in the south or not, by checking
        ! for the substrings ANT/ant in the domain name 
        if (trim(domain) .eq. "Antarctica") then 
            south = .TRUE. 
        else 
            south = .FALSE.
        end if 

        ! Call the appropriate forcing subroutine to obtain anomalies 
        select case(trim(snp%par%clim_type))
            case("const") 
                ! Constant steady-state climate
                call forclim_noforcing(snp%now,snp%clim0)

            case("hybrid") 
                ! Use hybrid input forcing curve (predetermined - no snapshots)
                dT = series_2D_interp(snp%hybrid%dTmon, year_bp)
                call forclim_hybrid(snp%now,snp%clim0,dT,to_fac=snp%hybrid%f_to,f_p=snp%par%f_p,is_south=south)

                ! Update external indices 
                at = sum(dT)/real(size(dT,1))   ! Annual mean 
                ao = at*snp%hybrid%f_to 

            case("anom")
                ! Apply a simple anomaly to clim0 

                ! First calculate the current anomaly 
                dTa_now = at*snp%par%dTa_const
                dTo_now = ao*snp%par%dTo_const

                ! Overwrite with input arguments if available
                if (present(dTa)) dTa_now = dTa 
                if (present(dTo)) dTo_now = dTo  

                call forclim_anom(snp%now,snp%clim0,dTa_now,dTo_now,snp%par%f_p)

            case("anom_1ind")
                call forclim_anom_1ind(snp%now,snp%clim0,snp%clim1,snp%clim2,at,ap,ao)

            case("anom_2ind")
                call forclim_anom_2ind(snp%now,snp%clim0,snp%clim1,snp%clim2,snp%clim3,at,ap,ao,bt,bp,bo)

            case("abs_1ind")
                call forclim_abs_1ind(snp%now,snp%clim1,snp%clim2,at,ap,ao)

            case("abs_2ind")
                call forclim_abs_2ind(snp%now,snp%clim1,snp%clim2,snp%clim3,at,ap,ao,bt,bp,bo) 
            
            case DEFAULT 
                write(*,*) "snapclim_update:: error: "// &
                           "clim_type not recognized: "// trim(snp%par%clim_type) 
                stop 

        end select 

        ! Correct sea-level temps and precip to current surface elevation
        snp%now%ta_ann = snp%now%tsl_ann - snp%par%lapse(1)*z_srf
        snp%now%ta_sum = snp%now%tsl_sum - snp%par%lapse(2)*z_srf
        snp%now%pr_ann = snp%now%prcor_ann*exp(snp%par%f_p*(snp%now%ta_ann-snp%now%tsl_ann))
        snp%now%to_ann = snp%now%to_ann ! no correction to sea level needed in the ocean

        ! Monthly values 
        do m = 1, 12 
            if (south) then 
                snp%now%tas(:,:,m) = snp%now%tsl(:,:,m) - &
                    z_srf*(snp%par%lapse(1)+(snp%par%lapse(2)-snp%par%lapse(1))*cos(2*pi*(m*30.0-30.0)/360.0))
            else 
                snp%now%tas(:,:,m) = snp%now%tsl(:,:,m) - &
                    z_srf*(snp%par%lapse(1)+(snp%par%lapse(1)-snp%par%lapse(2))*cos(2*pi*(m*30.0-30.0)/360.0))
            end if 
        end do 

        snp%now%pr = snp%now%prcor*exp(snp%par%f_p*(snp%now%tas-snp%now%tsl))

        ! jablasco: artificial accumulation
        nx = size(snp%now%pr_ann,1)
        ny = size(snp%now%pr_ann,2)
        if (snp%par%artificial_acc) then
           do i=2,nx-1
              do j=2,ny-1
                   if(i.eq.69 .and. j.eq.47) then !Ross=[68,46]
                      snp%now%pr_ann(i,j) = snp%par%acc_ross*(12000.0/365.0)/32.0
                   else if(i.eq.54 .and. j.eq.88) then !Ronne=[53,87]
                      snp%now%pr_ann(i,j) = snp%par%acc_ronne*(12000.0/365.0)/32.0
                   end if
              end do
           end do
        end if

        ! Store current index values in object for output 
        snp%now%at = at 
        snp%now%ap = ap 
        snp%now%ao = ao 
        snp%now%bt = bt 
        snp%now%bp = bp
        snp%now%bo = bo 
        
        return 

    end subroutine snapclim_update


    ! ======================================================================
    !
    ! The routines below return the current 'sea-level' annual mean temp,
    ! summer temp and annual precip. based on temporally changing weighting
    ! coefficients:
    ! 1ind = 1 index
    ! 2ind = 2 indices 
    ! anom = anomaly based forcing
    ! abs  = absolute based forcing 
    !
    ! ======================================================================

    subroutine forclim_noforcing(clim_now,clim0)
        ! Current climate equals clim0 

        implicit none

        type(snapclim_state_class), intent(OUT) :: clim_now 
        type(snapclim_state_class), intent(IN)  :: clim0

        clim_now = clim0 

        return

    end  subroutine forclim_noforcing

    subroutine forclim_anom(clim_now,clim0,dTa,dTo,f_p)
        ! Calculate current climate from clim0 plus anomaly 

        implicit none 

        type(snapclim_state_class), intent(OUT) :: clim_now
        type(snapclim_state_class), intent(IN)  :: clim0
        real(prec), intent(IN) :: dTa, dTo    ! Current anomaly value for atm and ocn
        real(prec), intent(IN) :: f_p         ! [frac/K] Precip scalar (eg f_p=5)

        clim_now = clim0 
        clim_now%tsl_ann   = clim0%tsl_ann + dTa 
        clim_now%tsl_sum   = clim0%tsl_sum + dTa
        clim_now%to_ann    = clim0%to_ann  + dTo
        clim_now%prcor_ann = clim0%prcor_ann*exp(f_p*dTa) 

        ! Monthly values 
        clim_now%tsl       = clim0%tsl + dTa 
        clim_now%prcor     = clim0%prcor*exp(f_p*dTa)
        
!         clim_now%prcor_ann = clim0%prcor_ann*(1+ap*0.01)

        return 
    end subroutine forclim_anom

    subroutine forclim_anom_1ind(clim_now,clim0,clim1,clim2,at,ap,ao)

        implicit none

        type(snapclim_state_class), intent(INOUT) :: clim_now
        type(snapclim_state_class), intent(IN)    :: clim0, clim1, clim2 
        real(prec) :: at, ap, ao   

        clim_now%tsl_ann = clim0%tsl_ann + (1.0-at)*(clim2%tsl_ann-clim1%tsl_ann)
        clim_now%tsl_sum = clim0%tsl_sum + (1.0-at)*(clim2%tsl_sum-clim1%tsl_sum)
        clim_now%to_ann  = clim0%to_ann  + (1.0-ao)*(clim2%to_ann-clim1%to_ann)

        where (clim1%prcor_ann.ne.0.0)
            clim_now%prcor_ann = clim0%prcor_ann*((1.0-ap)*(clim2%prcor_ann/clim1%prcor_ann) + ap)
        elsewhere 
            clim_now%prcor_ann = clim0%prcor_ann
        endwhere

        ! Monthly values 
        clim_now%tsl = clim0%tsl + (1.0-at)*(clim2%tsl-clim1%tsl)
        clim_now%prcor = clim0%prcor 
        where (clim1%prcor.ne.0.0) &
            clim_now%prcor = clim0%prcor*((1.0-ap)*(clim2%prcor/clim1%prcor) + ap)

        return

    end subroutine forclim_anom_1ind

    subroutine forclim_anom_2ind(clim_now,clim0,clim1,clim2,clim3,at,ap,ao,bt,bp,bo)

        implicit none

        type(snapclim_state_class), intent(INOUT) :: clim_now
        type(snapclim_state_class), intent(IN)    :: clim0, clim1, clim2, clim3 
        real(prec) :: at, ap, ao, bt, bp, bo

        ! option 1 : strictly correct

!       clim_now%tsl_ann = clim0%tsl_ann + (1.0-at)*(((1.0-bt)*(clim2%tsl_ann)+(bt*clim3%tsl_ann))-clim1%tsl_ann)
!       clim_now%tsl_sum = clim0%tsl_sum + (1.0-at)*(((1.0-bt)*(clim2%tsl_sum)+(bt*clim3%tsl_sum))-clim1%tsl_sum)
!       clim_now%to_ann  = clim0%to_ann  + (1.0-ao)*(((1.0-bo)*(clim2%to_ann) +(bo*clim3%to_ann)) -clim1%to_ann)

        ! option 2 : deconvoluted
                    
         clim_now%tsl_ann = clim0%tsl_ann + &
                     ((1.0-at)*(clim2%tsl_ann-clim1%tsl_ann))+(bt*(clim3%tsl_ann-clim2%tsl_ann))
         clim_now%tsl_sum = clim0%tsl_sum + &
                     ((1.0-at)*(clim2%tsl_sum-clim1%tsl_sum))+(bt*(clim3%tsl_sum-clim2%tsl_sum))

         clim_now%to_ann  = clim0%to_ann  + &
                     ((1.0-ao)*(clim2%to_ann-clim1%to_ann)) +(bo*(clim3%to_ann-clim2%to_ann))

        where (clim1%prcor_ann(:,:).ne.0)
            clim_now%prcor_ann = clim0%prcor_ann* &
                        ((1.0-ap)*(((1.0-bp)*(clim2%prcor_ann)+(bp*clim3%prcor_ann))/clim1%prcor_ann)+ap)
        elsewhere 
            clim_now%prcor_ann = clim0%prcor_ann
        endwhere

        ! Monthly values 
        clim_now%tsl = clim0%tsl + &
            ((1.0-at)*(clim2%tsl-clim1%tsl))+(bt*(clim3%tsl-clim2%tsl))
        clim_now%prcor = clim0%prcor 
        where (clim1%prcor.ne.0.0) &
            clim_now%prcor = clim0%prcor* &
                ((1.0-ap)*(((1.0-bp)*(clim2%prcor)+(bp*clim3%prcor))/clim1%prcor)+ap)

        return

    end subroutine forclim_anom_2ind

    subroutine forclim_abs_1ind(clim_now,clim1,clim2,at,ap,ao)

        implicit none

        type(snapclim_state_class), intent(INOUT) :: clim_now
        type(snapclim_state_class), intent(IN)    :: clim1, clim2
        real(prec) :: at, ap, ao

        clim_now%tsl_ann   = (1.0-at)*clim2%tsl_ann   + at*clim1%tsl_ann    
        clim_now%tsl_sum   = (1.0-at)*clim2%tsl_sum   + at*clim1%tsl_sum  
        clim_now%to_ann    = (1.0-ao)*clim2%to_ann    + ao*clim1%to_ann  
        clim_now%prcor_ann = (1.0-ap)*clim2%prcor_ann + ap*clim1%prcor_ann

        ! Monthly values 
        clim_now%tsl   = (1.0-at)*clim2%tsl   + at*clim1%tsl 
        clim_now%prcor = (1.0-ap)*clim2%prcor + ap*clim1%prcor 

        return

    end subroutine forclim_abs_1ind

    subroutine forclim_abs_2ind(clim_now,clim1,clim2,clim3,at,ap,ao,bt,bp,bo)

        implicit none

        type(snapclim_state_class), intent(INOUT) :: clim_now
        type(snapclim_state_class), intent(IN)    :: clim1, clim2, clim3
        real(prec) :: at, ap, ao, bt, bp, bo 

        clim_now%tsl_ann    = (1.0-at)*((1.0-bt)*(clim2%tsl_ann)  + bt*clim3%tsl_ann)   + at*clim1%tsl_ann
        clim_now%tsl_sum    = (1.0-at)*((1.0-bt)*(clim2%tsl_sum)  + bt*clim3%tsl_sum)   + at*clim1%tsl_sum
        clim_now%to_ann     = (1.0-ao)*((1.0-bo)*(clim2%to_ann)   + bo*clim3%to_ann)    + ao*clim1%to_ann
        clim_now%prcor_ann  = (1.0-ap)*((1.0-bp)*(clim2%prcor_ann)+ bp*clim3%prcor_ann) + ap*clim1%prcor_ann

        ! Monthly values 
        clim_now%tsl    = (1.0-at)*((1.0-bt)*(clim2%tsl)  + bt*clim3%tsl)   + at*clim1%tsl
        clim_now%prcor  = (1.0-ap)*((1.0-bp)*(clim2%prcor)+ bp*clim3%prcor) + ap*clim1%prcor

        return

    end subroutine forclim_abs_2ind

    subroutine forclim_hybrid(clim_now,clim0,dT,to_fac,f_p,is_south)

        implicit none 

        type(snapclim_state_class), intent(INOUT) :: clim_now
        type(snapclim_state_class), intent(IN)    :: clim0
        real(prec),                 intent(IN)    :: dT(:)
        real(prec),                 intent(IN)    :: to_fac 
        real(prec),                 intent(IN)    :: f_p    
        logical,                    intent(IN)    :: is_south 

        ! Local variables 
        real(prec) :: dT_ann, dT_sum 
        integer :: m 

        dT_ann = sum(dT) / real(size(dT,1))

        if (is_south) then 
            dT_sum = sum(dT([1,2,12]))/3.0
        else 
            dT_sum = sum(dT([6,7,8]))/3.0 
        end if 
        
        !write(*,*) "forclim_hybrid:: ", dT 

        ! Get annual and summer anomalies
        clim_now%tsl_ann = clim0%tsl_ann + dT_ann 
        clim_now%tsl_sum = clim0%tsl_sum + dT_sum 

        ! No change to oceanic temps and precip scaling 
        clim_now%to_ann    = clim0%to_ann + dT_ann*to_fac 
        
        ! Clausius-Clapeyron scaling of precipitation based on temp change
        clim_now%prcor_ann = clim0%prcor_ann*exp(f_p*dT_ann)

        ! Monthly values 
        do m = 1, 12 
            clim_now%tsl(:,:,m)   = clim0%tsl(:,:,m) + dT(m) 
            clim_now%prcor(:,:,m) = clim0%prcor(:,:,m)*exp(f_p*dT(m)) 
        end do 

        return 

    end subroutine forclim_hybrid

    subroutine snapclim_allocate(clim,nx,ny)

        implicit none 

        type(snapclim_state_class) :: clim 
        integer :: nx, ny

        ! Climatic variables
        if (allocated(clim%tas))       deallocate(clim%tas)
        if (allocated(clim%tsl))       deallocate(clim%tsl)
        if (allocated(clim%pr))        deallocate(clim%pr)
        if (allocated(clim%prcor))     deallocate(clim%prcor)
        if (allocated(clim%sf))        deallocate(clim%sf)


        if (allocated(clim%ta_ann))    deallocate(clim%ta_ann)
        if (allocated(clim%tsl_ann))   deallocate(clim%tsl_ann)
        if (allocated(clim%ta_sum))    deallocate(clim%ta_sum)
        if (allocated(clim%tsl_sum))   deallocate(clim%tsl_sum)
        if (allocated(clim%pr_ann))    deallocate(clim%pr_ann)
        if (allocated(clim%prcor_ann)) deallocate(clim%prcor_ann)
        if (allocated(clim%z_srf))     deallocate(clim%z_srf)
        if (allocated(clim%mask))      deallocate(clim%mask)
        
        allocate(clim%tas(nx,ny,12))
        allocate(clim%tsl(nx,ny,12))
        allocate(clim%pr(nx,ny,12))
        allocate(clim%prcor(nx,ny,12))
        allocate(clim%sf(nx,ny,12))
        
        allocate(clim%ta_ann(nx,ny))
        allocate(clim%tsl_ann(nx,ny))
        allocate(clim%ta_sum(nx,ny))
        allocate(clim%tsl_sum(nx,ny))
        allocate(clim%pr_ann(nx,ny))
        allocate(clim%prcor_ann(nx,ny))
        allocate(clim%z_srf(nx,ny))
        allocate(clim%mask(nx,ny))
        
        ! Initialize to zero
        clim%tas   = 0.0
        clim%tsl   = 0.0
        clim%pr    = 0.0
        clim%prcor = 0.0
        clim%sf    = 0.0
        
        clim%ta_ann    = 0.0
        clim%tsl_ann   = 0.0
        clim%ta_sum    = 0.0
        clim%tsl_sum   = 0.0
        clim%pr_ann    = 0.0
        clim%prcor_ann = 0.0
        clim%z_srf     = 0.0
        clim%mask      = 0.0
        
        
        return 

    end subroutine snapclim_allocate 
    
    subroutine snapclim_par_load(par,hpar,filename,init)
        ! == Parameters from namelist file ==

        implicit none 

        type(snapclim_param_class), intent(OUT) :: par
        type(hybrid_class),         intent(OUT) :: hpar 
        character(len=*), intent(IN) :: filename 
        logical, optional :: init 
        logical :: init_pars 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 

        call nml_read(filename,"snapclim","clim_type",          par%clim_type,      init=init_pars)
        call nml_read(filename,"snapclim","fname_at",           par%fname_at,      init=init_pars)
        call nml_read(filename,"snapclim","fname_ap",           par%fname_ap,      init=init_pars)
        call nml_read(filename,"snapclim","fname_ao",           par%fname_ao,      init=init_pars)
        call nml_read(filename,"snapclim","fname_bt",           par%fname_bt,      init=init_pars)
        call nml_read(filename,"snapclim","fname_bp",           par%fname_bp,      init=init_pars)
        call nml_read(filename,"snapclim","fname_bo",           par%fname_bo,      init=init_pars)
        call nml_read(filename,"snapclim","spatial_const_clim", par%spatial_const_clim,init=init_pars)
        call nml_read(filename,"snapclim","ta_clim1",           par%ta_clim1,       init=init_pars)
        call nml_read(filename,"snapclim","ta_clim2",           par%ta_clim2,       init=init_pars)
        call nml_read(filename,"snapclim","ta_clim3",           par%ta_clim3,       init=init_pars)
        call nml_read(filename,"snapclim","spatial_const_ocn",  par%spatial_const_ocn, init=init_pars)
        call nml_read(filename,"snapclim","to_clim1",           par%to_clim1,       init=init_pars)
        call nml_read(filename,"snapclim","to_clim2",           par%to_clim2,       init=init_pars)
        call nml_read(filename,"snapclim","to_clim3",           par%to_clim3,       init=init_pars)
        call nml_read(filename,"snapclim","artificial_acc",     par%artificial_acc, init=init_pars)
        call nml_read(filename,"snapclim","acc_ross",           par%acc_ross,       init=init_pars)
        call nml_read(filename,"snapclim","acc_ronne",          par%acc_ronne,      init=init_pars)
        call nml_read(filename,"snapclim","lapse",              par%lapse,          init=init_pars)
        call nml_read(filename,"snapclim","dTa_const",          par%dTa_const,      init=init_pars)
        call nml_read(filename,"snapclim","dTo_const",          par%dTo_const,      init=init_pars)
        call nml_read(filename,"snapclim","f_p",                par%f_p,            init=init_pars)
        
        call nml_read(filename,"snapclim_hybrid","hybrid_path", hpar%hybrid_path,  init=init_pars)
        call nml_read(filename,"snapclim_hybrid","f_eem",       hpar%f_eem,        init=init_pars)
        call nml_read(filename,"snapclim_hybrid","f_glac",      hpar%f_glac,       init=init_pars)
        call nml_read(filename,"snapclim_hybrid","f_hol",       hpar%f_hol,        init=init_pars)
        call nml_read(filename,"snapclim_hybrid","f_seas",      hpar%f_seas,       init=init_pars)
        call nml_read(filename,"snapclim_hybrid","f_to",        hpar%f_to,         init=init_pars)
        
        return 

    end subroutine snapclim_par_load

    subroutine snapshot_par_load(par,filename,name,domain,grid_name,init)

        implicit none 

        type(snapshot_param_class), intent(OUT) :: par
        character(len=*), intent(IN) :: filename, name
        character(len=*), intent(IN) :: domain, grid_name 
        logical, optional :: init 
        logical :: init_pars 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 

        call nml_read(filename,name,"clim_path",    par%clim_path,     init=init_pars)
        call nml_read(filename,name,"clim_names",   par%clim_names,    init=init_pars)
        call nml_read(filename,name,"clim_monthly", par%clim_monthly,  init=init_pars)
        call nml_read(filename,name,"clim_year_bp", par%clim_year_bp,  init=init_pars)
        call nml_read(filename,name,"ocn_path",     par%ocn_path,      init=init_pars)
        call nml_read(filename,name,"ocn_names",    par%ocn_names,     init=init_pars)
        call nml_read(filename,name,"ocn_monthly",  par%ocn_monthly,   init=init_pars)
        call nml_read(filename,name,"ocn_year_bp",  par%ocn_year_bp,   init=init_pars)
        
        ! Subsitute domain/grid_name (equivalent to yelmo_parse_path)
        call parse_path(par%clim_path,domain,grid_name)
        call parse_path(par%ocn_path,domain,grid_name)
        
        return 

    end subroutine snapshot_par_load

    subroutine read_climate_snapshot(clim,nx,ny,lapse,f_p,domain)
        ! `names` is a vector of names in the netcdf file that 
        ! correspond to the fields to be read in:
        ! (1) 2D elevation field
        ! (2) 2D [nx,ny] or 3D [nx,ny,nmonth] near-surface temperature field 
        ! (3) 2D [nx,ny] or 3D [nx,ny,nmonth] precipitation, or snowfall field 
        ! (4) 2D [nx,ny] or 3D [nx,ny,nmonth] rainfall field, in case snowfall
        !     was provided in field 3.
        implicit none 

        type(snapclim_state_class), intent(INOUT) :: clim
        integer,          intent(IN) :: nx, ny  
        real(prec),       intent(IN) :: lapse(2)
        real(prec),       intent(IN) :: f_p 
        character(len=*), intent(IN) :: domain   
        
        ! Local variables
        real(prec) :: lapse_mon(12)   
        real(prec), allocatable :: rf(:,:,:)
        character(len=56)  :: tmp_str 
        integer :: m 
        logical :: south 

        south = .FALSE. 
        if (trim(domain).eq."Antarctica") south = .TRUE. 

        ! (Re)allocate the clim object
        call clim_allocate(clim,nx,ny)

        ! Allocate helper array
        allocate(rf(nx,ny,12))

        if (trim(clim%par%clim_path) .eq. "None") then 
            ! Simply set fields to zero 

            clim%par%clim_year_bp = 0.0 

            clim%z_srf     = 0.0 
            clim%mask      = 0.0 
            clim%tas       = 0.0 
            clim%pr        = 0.0 
            clim%sf        = 0.0 
            clim%ta_ann    = 0.0 
            clim%ta_sum    = 0.0 
            
            clim%tsl       = 0.0 
            clim%tsl_ann   = 0.0 
            clim%tsl_sum   = 0.0 
            clim%prcor     = 0.0 
            clim%prcor_ann = 0.0 
            
        else 
            ! Read data using parameters

            ! Read in the elevation
            call nc_read(clim%par%clim_path,clim%par%clim_names(1),clim%z_srf)

            ! Define mask based on elevations 
            clim%mask = 0.0 
            where(clim%z_srf .gt. 0.0) clim%mask = 1.0 

            if (clim%par%clim_monthly) then 
                ! Read in monthly climate fields, then get the averages

                call nc_read(clim%par%clim_path,clim%par%clim_names(2),clim%tas)
                call nc_read_attr(clim%par%clim_path,clim%par%clim_names(2),"units",tmp_str)

                if (minval(clim%tas) .lt. 100.0) then
                    clim%tas=clim%tas+273.15
                end if

                call nc_read(clim%par%clim_path,clim%par%clim_names(3),clim%sf)
                ! Check the rainfall field also needs to be loaded, or if 
                ! precip was provided as the second name 
                if (trim(clim%par%clim_names(3)) .eq. "sf") then 
                    call nc_read(clim%par%clim_path,clim%par%clim_names(4),rf)
                    clim%pr = clim%sf+rf        ! [mm/d]
                else
                    clim%pr = clim%sf           ! [mm/d]
                end if 

                ! Calculate derived fields from monthly fields...

                ! Annual mean temperature [K]
                clim%ta_ann = sum(clim%tas,dim=3)
                clim%ta_ann = clim%ta_ann/12.0 

                ! Summer (July or Jan) temperature [K]
                if (south) then 
                    clim%ta_sum = sum(clim%tas(:,:,[12,1,2]),dim=3)
                else 
                    clim%ta_sum = sum(clim%tas(:,:,[6,7,8]),dim=3)
                end if 
                clim%ta_sum = clim%ta_sum/3.0 

                ! Annual mean precipitation [mm/a]
                clim%pr_ann = sum(clim%pr,dim=3)/12.0 *365.0  ! [mm/a]

            else
                ! Read in specific derived climate fields 

                call nc_read(clim%par%clim_path,clim%par%clim_names(2),clim%ta_ann)
                call nc_read(clim%par%clim_path,clim%par%clim_names(3),clim%ta_sum)
                call nc_read_attr(clim%par%clim_path,clim%par%clim_names(2),"units",tmp_str)
                
                if (minval(clim%ta_ann) .lt. 100.0) then
                   clim%ta_ann = clim%ta_ann + 273.15
                   clim%ta_sum = clim%ta_sum + 273.15
                end if

                call nc_read(clim%par%clim_path,clim%par%clim_names(4),clim%pr_ann)
       
                clim%pr_ann = clim%pr_ann *365.0  ! [mm/d] => [mm/a]

                ! Translate to monthly fields to have them available 
                do m = 1, 12

                    ! Summer (July or Jan) temperature [degrees C]
                    if (south) then 
                        clim%tas(:,:,m) = clim%ta_ann + &
                            (clim%ta_sum-clim%ta_ann)*cos(2*pi*(m*30.0-30.0)/360.0)
                    else 
                        clim%tas(:,:,m) = clim%ta_ann - &
                            (clim%ta_sum-clim%ta_ann)*cos(2*pi*(m*30.0-30.0)/360.0)
                    
                    end if 

                    ! Precip distributed evenly throughout the year 
                    clim%pr(:,:,m) = clim%pr_ann/365.0             ! [mm/d]
                    clim%sf(:,:,m) = clim%pr(:,:,m)/365.0          ! [mm/d]
                end do 
        
            end if 

            ! ====================================
            !
            ! Additional calculations
            !
            ! ====================================

            ! Calculate sea-level temperatures 
            clim%tsl_ann = clim%ta_ann + lapse(1)*clim%z_srf 
            clim%tsl_sum = clim%ta_sum + lapse(2)*clim%z_srf 

            ! Calculate sea-level corrected precip 
            clim%prcor_ann = clim%pr_ann/exp(f_p*(clim%ta_ann-clim%tsl_ann))

            do m = 1, 12

                ! Monthly lapse rates from ann and sum (jan or jul)
                if (south) then 
                    clim%tsl(:,:,m) = clim%tas(:,:,m) + &
                        clim%z_srf*(lapse(1)+(lapse(2)-lapse(1))*cos(2*pi*(m*30.0-30.0)/360.0))
                else 
                    clim%tsl(:,:,m) = clim%tas(:,:,m) + &
                        clim%z_srf*(lapse(1)+(lapse(1)-lapse(2))*cos(2*pi*(m*30.0-30.0)/360.0))
                end if 

                ! Precip
                clim%prcor(:,:,m) = clim%pr(:,:,m) / exp(f_p*(clim%tas(:,:,m)-clim%tsl(:,:,m)))   ! [m/a]

            end do 
        
        end if 

        ! ====================================
        !
        ! Summary
        !
        ! ====================================

        write(*,"(a)") "read_climate_snapshot:: Climate snapshot loaded: "//trim(clim%par%clim_path)
        write(*,"(4x,a,2f12.2)") "range: ta_ann     [K]    = ", minval(clim%ta_ann),    maxval(clim%ta_ann)
        write(*,"(4x,a,2f12.2)") "range: tsl_ann    [K]    = ", minval(clim%tsl_ann),   maxval(clim%tsl_ann)
        write(*,"(4x,a,2f12.2)") "range: ta_sum     [K]    = ", minval(clim%ta_sum),    maxval(clim%ta_sum)
        write(*,"(4x,a,2f12.2)") "range: tsl_sum    [K]    = ", minval(clim%tsl_sum),   maxval(clim%tsl_sum)
        write(*,"(4x,a,3f12.2)") "range: pr_ann     [mm/a] = ", minval(clim%pr_ann),    maxval(clim%pr_ann), sum(clim%pr_ann)/count(clim%pr_ann.gt.-999.9)
        write(*,"(4x,a,3f12.2)") "range: prcor_ann  [mm/a] = ", minval(clim%prcor_ann), maxval(clim%prcor_ann), sum(clim%prcor_ann)/count(clim%prcor_ann.gt.-999.9)
        write(*,"(4x,a,2f12.2)") "range: z_srf      [m]    = ", minval(clim%z_srf),     maxval(clim%z_srf)
        write(*,"(a,f12.2)") "snapshot year_bp =", clim%par%clim_year_bp       
 
        return 

    end subroutine read_climate_snapshot

    subroutine read_ocean_snapshot(ocn,nx,ny,depth)
        ! Read oceanic snapshot from netcdf file
        ! Parameters loaded earlier via snapshot_par_load()
        ! The oceanic snapshot needs the following information:
        ! (1) depth dimension of 3D/4D field (ie, nz of [nx,ny,nz])
        ! (2) 3D land-ocean mask specifying original boundaries of the simulation
        ! (3) 3D [nx,ny,nz0] or 4D [nx,ny,nz0,nmonth] ocean temperature field with no missing values
        !     (ie, interpolated even where there was land to be safe) 
        !
        ! Subroutine will output the ocn (climate) object with the input field [nx,ny,nz0]
        ! interpolated vertically to the desired depth levels [nx,ny,nz]

        implicit none 

        type(snapclim_state_class), intent(INOUT) :: ocn 
        integer,    intent(IN) :: nx, ny
        real(prec), intent(IN) :: depth(:)

        ! Local helping variables 
        real(prec), allocatable :: depth0(:) 
        integer,    allocatable :: mask3D(:,:,:)
        real(prec), allocatable :: tocn3D(:,:,:)
        real(prec), allocatable :: tocn4D(:,:,:,:)
        integer :: nz, nz0, m 
        integer :: i, j, k 

        if (trim(ocn%par%ocn_path) .eq. "None") then 
            ! Set all fields to zero 

            ocn%nzo      = size(depth)
            call ocn_allocate(ocn,nx,ny,ocn%nzo)

            ocn%depth    = 0.0
            ocn%to_ann   = 0.0 
            ocn%mask_ocn = 0.0 

            ocn%par%ocn_year_bp = 0.0 

        else 
            ! Read data from parameters 

            ! Determine the length of the input depth dimension 
            nz0 = nc_size(ocn%par%ocn_path,ocn%par%ocn_names(1))
            allocate(depth0(nz0))
            allocate(mask3D(nx,ny,nz0))
            allocate(tocn3D(nx,ny,nz0))

            ! Read in the input ocean mask and depth
            call nc_read(ocn%par%ocn_path,ocn%par%ocn_names(1),depth0)
            call nc_read(ocn%par%ocn_path,ocn%par%ocn_names(2),mask3D)

            if (ocn%par%ocn_monthly) then 
                ! Read in monthly climate fields, then get the averages

                allocate(tocn4D(nx,ny,nz0,12))
                call nc_read(ocn%par%ocn_path,ocn%par%ocn_names(3),tocn4D)

                ! Calculate annual mean temperature [degrees C]
                tocn3D = sum(tocn4D,dim=4)
                tocn3D = tocn3D/12.0 

                ! Delete working array to be safe
                deallocate(tocn4D) 

            else
                ! Read in specific climate fields 

                call nc_read(ocn%par%ocn_path,ocn%par%ocn_names(3),tocn3D)
        
            end if 

            ! Ensure units are in Kelvin
            if (minval(tocn3D) .lt. 100.0) then
                tocn3D = tocn3D + 273.15 
            end if 

            ! Determine the length of output depth dimension
            ! and allocate the ocn object
            ocn%nzo   = size(depth)
            call ocn_allocate(ocn,nx,ny,ocn%nzo)
            ocn%depth = depth 

            ! Interpolate vertical levels to those needed by the model
            do k = 1, size(ocn%depth)
                do j = 1, ny 
                    do i = 1, nx 

                        ! Interpolate ocean temp to current depth level (negative to convert from depth to height)
                        ocn%to_ann(i,j,k) = interp_linear(depth0,tocn3D(i,j,:),xout=ocn%depth(k))

                        ! Get nearest mask value 
                        if (ocn%depth(k) .gt. maxval(depth0)) then
                            ! If the model ocean depth does not exist, set mask to land
                            ocn%mask_ocn(i,j,k)=0

                        else 
                            ! Interpolate ocean mask to current depth 
                            ocn%mask_ocn(i,j,k) = &
                                nint(interp_linear(depth0,real(mask3D(i,j,:)),xout=ocn%depth(k)))
                        
                        end if 
                        
                    end do 
                end do
            end do

        end if


        ! ====================================
        !
        ! Summary
        !
        ! ====================================

        write(*,"(a)") "read_ocean_snapshot:: Ocean snapshot loaded: "//trim(ocn%par%ocn_path)
        write(*,"(4x,a,2f12.2)") "depth  range = ", minval(ocn%depth),  maxval(ocn%depth)
        write(*,"(4x,a,2f12.2)") "to_ann range = ", minval(ocn%to_ann), maxval(ocn%to_ann)
        write(*,"(a,f12.2)") "snapshot year_bp =", ocn%par%ocn_year_bp 

        return 

    end subroutine read_ocean_snapshot

    subroutine set_climate_snapshot(clim,clim0,nx,ny,year_bp,tas,f_p)

        implicit none 

        type(snapclim_state_class), intent(INOUT) :: clim
        type(snapclim_state_class), intent(IN)    :: clim0 

        integer, intent(IN) :: nx, ny 
        real(prec), intent(IN) :: year_bp 
        real(prec), intent(IN) :: tas 
        real(prec), intent(IN) :: f_p 
        integer :: m 

        ! (Re)allocate the clim object
        call clim_allocate(clim,nx,ny)
        
        clim%mask      = 0.0 
        clim%z_srf     = 0.0

        clim%tas       = tas 
        clim%ta_ann    = tas 
        clim%ta_sum    = tas 
        clim%tsl_ann   = tas 
        clim%tsl_sum   = tas 

        clim%prcor_ann = clim0%prcor_ann*exp(f_p*(clim%ta_ann-clim0%tsl_ann))
        clim%pr_ann    = clim%prcor_ann 
        
        do m = 1, 12 
            clim%pr(:,:,m) = clim%prcor_ann
            clim%sf(:,:,m) = clim%prcor_ann
        end do 

        clim%par%clim_year_bp = year_bp 

        return 

    end subroutine set_climate_snapshot

    subroutine set_ocn_snapshot(ocn,nx,ny,depth,year_bp,t_ocn)

        implicit none 

        type(snapclim_state_class), intent(INOUT) :: ocn 
        integer, intent(IN) :: nx, ny
        real(prec), intent(IN) :: depth(:)
        real(prec), intent(IN) :: year_bp 
        real(prec), intent(IN) :: t_ocn  

        ! Determine the length of output depth dimension
        ! and allocate the ocn object
        ocn%nzo      = size(depth)
        call ocn_allocate(ocn,nx,ny,ocn%nzo)
        ocn%depth    = depth 

        ocn%to_ann   = t_ocn 
        ocn%mask_ocn = 0 

        ocn%par%ocn_year_bp = year_bp 

        return 

    end subroutine set_ocn_snapshot

    subroutine read_series(series,filename)
        ! This subroutine will read a time series of
        ! two columns [time,var] from an ascii file.
        ! Header should be commented by "#" or "!"
        implicit none 

        type(series_type) :: series 
        character(len=*)  :: filename 

        integer, parameter :: f = 190
        integer, parameter :: nmax = 10000

        integer :: i, stat, n 
        character(len=256) :: str, str1 
        real(prec) :: x(nmax), y(nmax) 

        ! Open file for reading 
        open(f,file=filename,status="old")

        ! Read the header in the first line: 
        read(f,*,IOSTAT=stat) str

        do i = 1, nmax 
            read(f,'(a100)',IOSTAT=stat) str 

            ! Exit loop if the end-of-file is reached 
            if(IS_IOSTAT_END(stat)) exit 

            str1 = adjustl(trim(str))
!            str1=str
            if ( len(trim(str1)) .gt. 0 ) then 
                if ( .not. (str1(1:1) == "!" .or. &
                            str1(1:1) == "#") ) then 
                    read(str1,*) x(i), y(i) 
                end if
            end if  
        end do 


        ! Close the file
        close(f) 

        if (i .eq. nmax) then 
            write(*,*) "read_series:: warning: "// &
                       "Maximum length of time series reached, ", nmax
            write(*,*) "Time series in the file may be longer: ", trim(filename)
        end if 

        ! Allocate the time series object and store output data
        n = i-1 
        call series_allocate(series,n)

        series%time = x(1:n) 
        series%var  = y(1:n) 

        write(*,*) "read_series:: Time series read from file: "//trim(filename)
        write(*,*) "    range time: ",minval(series%time), maxval(series%time)
        write(*,*) "    range var : ",minval(series%var), maxval(series%var)
        
        return 

    end subroutine read_series 

    subroutine read_series_hybrid(series,filename,f_eem,f_glac,f_hol,f_seas)

        implicit none 

        type(series_2D_type), intent(OUT) :: series 
        character(len=*),     intent(IN)  :: filename 
        real(prec),              intent(IN)  :: f_eem,f_glac,f_hol,f_seas

        ! Local variables
        integer :: nm, nt, i  
        real(prec), allocatable :: dTann(:), dataset(:), dTseas(:,:)

        ! Initialize dimensions and data object 
        nm = 12 
        nt = nc_size(filename,"time")
        call series_2D_allocate(series,nt,nm)
        allocate(dTann(nt),dataset(nt))
        allocate(dTseas(nm,nt))

        call nc_read(filename,"time",series%time)
        series%time = series%time*1d3               ! ka => a 

        call nc_read(filename,"dTann",dTann)
        call nc_read(filename,"dTseas",dTseas)
        call nc_read(filename,"dataset",dataset)
        
        ! Scale temperature anomalies based on parameter scaling factors
        where(dataset .eq. 1.0 .and. dTann .gt. 0.0) dTann = dTann * f_hol 
        where(dataset .eq. 2.0 .and. dTann .lt. 0.0) dTann = dTann * f_glac

        ! Make sure the Eemian factor is applied to all positive temperatures during the Eemian only!
        where(series%time .ge. -130e3 .and. series%time .le. -115e3 .and. dTann .gt. 0.0) 
            dTann = dTann * f_eem 
        end where

        do i = 1,nm
            series%var(i,:) = dTann + (dTseas(i,:) * f_seas) 
        end do

        write(*,"(a)") "read_series_hybrid:: read anomalies from "//trim(filename)
        write(*,"(a15,2f12.2)") "range time: ", minval(series%time), maxval(series%time)
        write(*,"(a15,2f12.2)") "range dT: ", minval(series%var), maxval(series%var)

        return 

    end subroutine read_series_hybrid

    function series_2D_interp(series,year_bp) result(var)
        ! Wrapper for simple `interp_linear` function
        ! for series_2D_types. 
        implicit none 

        type(series_2D_type) :: series 
        real(prec) :: year_bp 
        real(prec) :: var(size(series%var,1))
        integer :: nt, nm, i 

        ! Interpolate series object
        nm = size(series%var,1)
        do i = 1, nm 
            var(i) = interp_linear(series%time,series%var(i,:),xout=year_bp)
        end do 

        return 

    end function series_2D_interp 

    function series_interp(series,year_bp) result(var)
        ! Wrapper for simple `interp_linear` function
        ! for series_types. 
        implicit none 

        type(series_type) :: series 
        real(prec) :: year_bp 
        real(prec) :: var 
        integer :: nt, i 

        ! Interpolate series object
        var = interp_linear(series%time,series%var,xout=year_bp)

        return 

    end function series_interp 

    function interp_linear(x,y,xout) result(yout)
        ! Simple linear interpolation of a point

        implicit none 
 
        real(prec), dimension(:), intent(IN) :: x, y
        real(prec), intent(IN) :: xout
        real(prec) :: yout 
        integer :: i, j, n, nout 
        real(prec) :: alph

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
        
    ! ==================================================
    !
    ! Basic object allocations, etc. below...
    !
    ! ==================================================

    subroutine clim_allocate(clim,nx,ny)

        implicit none 

        type(snapclim_state_class) :: clim 
        integer :: nx, ny

        ! Climatic variables
        if (allocated(clim%tas))       deallocate(clim%tas)
        if (allocated(clim%tsl))       deallocate(clim%tsl)
        if (allocated(clim%pr))        deallocate(clim%pr)
        if (allocated(clim%prcor))     deallocate(clim%prcor)
        if (allocated(clim%sf))        deallocate(clim%sf)


        if (allocated(clim%ta_ann))    deallocate(clim%ta_ann)
        if (allocated(clim%tsl_ann))   deallocate(clim%tsl_ann)
        if (allocated(clim%ta_sum))    deallocate(clim%ta_sum)
        if (allocated(clim%tsl_sum))   deallocate(clim%tsl_sum)
        if (allocated(clim%pr_ann))    deallocate(clim%pr_ann)
        if (allocated(clim%prcor_ann)) deallocate(clim%prcor_ann)
        if (allocated(clim%z_srf))     deallocate(clim%z_srf)
        if (allocated(clim%mask))      deallocate(clim%mask)
        
        allocate(clim%tas(nx,ny,12))
        allocate(clim%tsl(nx,ny,12))
        allocate(clim%pr(nx,ny,12))
        allocate(clim%prcor(nx,ny,12))
        allocate(clim%sf(nx,ny,12))
        
        allocate(clim%ta_ann(nx,ny))
        allocate(clim%tsl_ann(nx,ny))
        allocate(clim%ta_sum(nx,ny))
        allocate(clim%tsl_sum(nx,ny))
        allocate(clim%pr_ann(nx,ny))
        allocate(clim%prcor_ann(nx,ny))
        allocate(clim%z_srf(nx,ny))
        allocate(clim%mask(nx,ny))

        return 

    end subroutine clim_allocate 

    subroutine ocn_allocate(ocn,nx,ny,nz)

        implicit none 

        type(snapclim_state_class) :: ocn 
        integer :: nx, ny, nz 

        ! Oceanic variables
        if (allocated(ocn%depth))  deallocate(ocn%depth)
        if (allocated(ocn%to_ann)) deallocate(ocn%to_ann)
        if (allocated(ocn%mask_ocn))   deallocate(ocn%mask_ocn)
        
        allocate(ocn%depth(nz))
        allocate(ocn%to_ann(nx,ny,nz))
        allocate(ocn%mask_ocn(nx,ny,nz))
        
        if (.not. (nx .eq. size(ocn%ta_ann,1) .and. &
                   ny .eq. size(ocn%ta_ann,2)) ) then 

            write(*,*) "ocn_allocate:: error: ocn [nx,ny] dimensions do not &
                       &match clim [nx,ny] dimensions."
            write(*,*) "clim nx,ny: ", size(ocn%ta_ann,1), size(ocn%ta_ann,2)
            write(*,*) "ocn  nx,ny: ", nx, ny 
            stop 
        end if 

        return 

    end subroutine ocn_allocate 

    subroutine series_allocate(series,nt)

        implicit none 

        type(series_type) :: series 
        integer :: nt 

        if (allocated(series%time))  deallocate(series%time)
        if (allocated(series%var))   deallocate(series%var)
        if (allocated(series%sigma)) deallocate(series%sigma)

        allocate(series%time(nt))
        allocate(series%var(nt))
        allocate(series%sigma(nt))

        return 

    end subroutine series_allocate

    subroutine series_2D_allocate(series,nt,nm)

        implicit none 

        type(series_2D_type) :: series 
        integer :: nt, nm 

        if (allocated(series%time))  deallocate(series%time)
        if (allocated(series%var))   deallocate(series%var)
        if (allocated(series%sigma)) deallocate(series%sigma)

        allocate(series%time(nt))
        allocate(series%var(nm,nt))
        allocate(series%sigma(nm,nt))

        return 

    end subroutine series_2D_allocate

    subroutine parse_path(path,domain,grid_name)

        implicit none 

        character(len=*), intent(INOUT) :: path 
        character(len=*), intent(IN)    :: domain, grid_name 

        call nml_replace(path,"{domain}",   trim(domain))
        call nml_replace(path,"{grid_name}",trim(grid_name))
        
        return 

    end subroutine parse_path
    


! ajr: to move into its own forcing module
!             case("hyst")
!                 ! Hysteresis 
                
!                 if (hyst_atm%par%kill) then 
!                     write(*,*) "hyst_atm kill condition reached - stopping simulation!"
!                     stop 
!                 end if 

!                 ! Get constant climate from climatology 
!                 call forclim_noforcing(clim_now,clim0)

!                 ! Calculate current ice volume [m^3 ice] => [km^3] == [Gt]
!                 Vtot     = sum(H_ice,mask=H_ice.gt.0.0)*dx*dy *(1e3/910.0) *1e-9
!                 Vtot_flt = sum(H_ice,mask=H_ice.gt.0.0 .and. f_grnd.eq.0.0)*dx*dy *(1e3/910.0) *1e-9

!                 ! Determine rate of forcing change and new forcing value
!                 if (snp%hyst_atm%par%use_hyster) then 
!                     call hyster_calc_rate(snp%hyst_atm,time=dble(year_bp),var=dble(Vtot),verbose=.TRUE.)
!                     clim_now%tsl_ann = clim_now%tsl_ann + snp%hyst_atm%f_now
!                     clim_now%tsl_sum = clim_now%tsl_sum + snp%hyst_atm%f_now
!                     at = snp%hyst_atm%f_now 
!                 end if

!                 if (snp%hyst_ocn%par%use_hyster) then 
!                     call hyster_calc_rate(snp%hyst_ocn,time=dble(year_bp),var=dble(Vtot),verbose=.TRUE.)
!                     clim_now%to_ann  = clim_now%to_ann  + snp%hyst_ocn%f_now  
!                     ao = snp%hyst_ocn%f_now 
!                 end if 


end module snapclim 
