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

    real(wp), parameter :: G  = 9.81                  ! [m/s^2]
    real(wp), parameter :: pi = 3.14159265359

    real(wp), parameter :: rho_ice    = 910.0                 ! [kg/m^3]
    real(wp), parameter :: rho_ocean  = 1028.0                ! [kg/m^3] 
    real(wp), parameter :: rho_water  = 1000.0                ! [kg/m^3] 
    real(wp), parameter :: rho_mantle = 3300.0                ! [kg/m^3] 

    type isos_param_class 
        integer             :: method           ! Type of isostasy to use
        character(len=512)  :: fname_kelvin     ! File containing precalculated zero-order Kelvin function values
        real(wp)          :: dt               ! [yr] Timestep to recalculate bedrock uplift rate
        real(wp)          :: tau              ! [yr] Asthenospheric relaxation constant
        real(wp)          :: DL               ! [N-m] lithosphere flexural rigidity (D_lith)
        real(wp)          :: RL               ! [m] radius of relative stiffness (R_lith)
        
        ! Internal parameters 

        integer    :: lbloc     ! distance en noeuds sur laquelle on calcule
                                ! l'influence d'une charge sur la deflexion de 
                                ! la lithosphere.  

        
        real(wp) :: time 

    end type 

    type isos_state_class 
        
        real(wp), allocatable :: z_bed(:,:)       ! Bedrock elevation         [m]
        real(wp), allocatable :: dzbdt(:,:)       ! Rate of bedrock uplift    [m/a]
        real(wp), allocatable :: z_bed_ref(:,:)   ! Reference (unweighted) bedrock 

        real(wp), allocatable :: tau(:,:)           ! [yr] Asthenospheric relaxation timescale field
        real(wp), allocatable :: He_lith(:,:)       ! [m]  Effective elastic thickness of the lithosphere
        real(wp), allocatable :: D_lith(:,:)        ! [m]  [N-m] Lithosphere flexural rigidity
        
        real(wp), allocatable :: we(:,:)          ! enfoncement du socle autour
                                                    ! d'une charge unitaire, sera
                                                    ! dimensionne dans le module 
                                                    ! isostasie_mod, routine init_iso a 
                                                    ! WE(LBLOC:LBOC,-LBLOC:LBLOC)

        ! Relaxation
        real(wp), allocatable :: charge(:,:)      ! charge sur une maille = RO G H
                                                    ! unite : 
                                                    ! sera dimensionne dans
                                                    ! isostasie_mod 
                                                    ! CHARGE(1-LBLOC:NX+LBLOC,1-LBLOC,NY+LBLOC)

        real(wp), allocatable :: w0(:,:)          ! Current weighting
        real(wp), allocatable :: w1(:,:)          ! New weighting          

        
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

contains 

    subroutine isos_init(isos,filename,nx,ny,dx)

        implicit none 

        type(isos_class), intent(OUT) :: isos 
        character(len=*), intent(IN)  :: filename 
        integer, intent(IN) :: nx, ny 
        real(wp), intent(IN) :: dx

        integer :: LBLOC 

        ! Load parameters
        call isos_par_load(isos%par,filename,init=.TRUE.)
        
        ! LBLOC, 400 km on each side
        ! ajr: 400km set internally, since the radius should be smaller.
        ! This needs thorough revision and code refactoring. 
        ! See Greve and Blatter (2009), Chpt 8, page 192 for methodology 
        ! and Le Muer and Huybrechts (1996). It seems that this value
        ! should be larger to capture the forebuldge at 5-6x radius of relative stiffness
        isos%par%lbloc = int((400000.0-0.1)/dx)+1
        lbloc = isos%par%lbloc
        
        ! Initialize isos variables 
        call isos_allocate(isos%now,nx,ny,nrad=lbloc)
    
        ! Intially ensure all variables are zero 
        isos%now%we         = 0.0 
        isos%now%charge     = 0.0  
        isos%now%z_bed_ref  = 0.0
        isos%now%z_bed      = 0.0 
        isos%now%dzbdt      = 0.0 

        ! Set time to very large value in the future 
        isos%par%time       = 1e10 

        ! Load lithospheric table of Kelvin function values 
        call read_tab_litho(isos%now%we,filename=trim(isos%par%fname_kelvin), &
                            RL=isos%par%RL,DL=isos%par%DL,dx=dx,dy=dx)

        write(*,*) "isos_init:: range(WE):  ", minval(isos%now%we),     maxval(isos%now%we)

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

        ! Store initial values of parameters
        isos%now%tau = isos%par%tau 
        isos%now%He_lith = 1.0_wp   ! ajr: to do!
        isos%now%D_lith  = 1.0_wp   ! ajr: to do!


        ! Initialize the charge field to match reference topography
        call init_charge(isos%now%charge,H_ice_ref,z_bed_ref,z_sl_ref,isos%par%lbloc)

        select case(isos%par%method)

            case(0)
                ! Steady-state lithospheric depression 

                isos%now%w0 = isos%now%charge/(rho_mantle*G)
                isos%now%w1 = isos%now%w0

            case(1,2) 
                ! 1: LLRA - Local lithosphere, relaxing Aesthenosphere
                ! 2: ELRA - Elastic lithosphere, relaxing Aesthenosphere

                call litho(isos%now%w1,isos%now%charge,isos%now%we,isos%par%LBLOC)
                isos%now%w0 = isos%now%w1 

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
        write(*,*) "  range(charge): ", minval(isos%now%charge), maxval(isos%now%charge)
        write(*,*) "  range(w0):     ", minval(isos%now%w0), maxval(isos%now%w0)
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

                isos%now%w1     = isos%now%w0
                isos%now%dzbdt = 0.0 

            case(1)
                ! Local lithosphere, relaxing aesthenosphere (LLRA)

                ! Local lithosphere (LL)
                call calc_litho_local(isos%now%w1,isos%now%z_bed,H_ice,z_sl)

                ! Relaxing aesthenosphere (RA)
                isos%now%dzbdt = ((isos%now%z_bed_ref-isos%now%z_bed) - (isos%now%w1-isos%now%w0))/isos%par%tau

            case(2)
                ! Elastic lithosphere, relaxing aesthenosphere (ELRA)
                
                ! Regional elastic lithosphere (EL)
                call calc_litho_regional(isos%now%w1,isos%now%charge,isos%now%we,isos%now%z_bed,H_ice,z_sl)

                ! Relaxing aesthenosphere (RA)
                isos%now%dzbdt = ((isos%now%z_bed_ref-isos%now%z_bed) - (isos%now%w1-isos%now%w0))/isos%par%tau
            
            case(3) 
                ! Elementary GIA model (spatially varying ELRA with geoid - to do!)

                ! Regional elastic lithosphere (EL)
                call calc_litho_regional(isos%now%w1,isos%now%charge,isos%now%we,isos%now%z_bed,H_ice,z_sl)

                ! Aesthenosphere timescale field 
                isos%now%tau = isos%par%tau 

                ! Relaxing aesthenosphere (RA)
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
        call nml_read(filename,"isostasy","tau",            par%tau,            init=init_pars)
        call nml_read(filename,"isostasy","D_lith",         par%DL,             init=init_pars)
        call nml_read(filename,"isostasy","R_lith",         par%RL,             init=init_pars)
        
        return

    end subroutine isos_par_load

    subroutine isos_allocate(now,nx,ny,nrad)

        implicit none 

        type(isos_state_class), intent(INOUT) :: now 
        integer, intent(IN) :: nx, ny, nrad  

        ! First ensure arrays are not allocated
        call isos_deallocate(now)

        ! Allocate arrays

        allocate(now%z_bed(nx,ny))
        allocate(now%dzbdt(nx,ny))
        allocate(now%z_bed_ref(nx,ny))
        
        allocate(now%tau(nx,ny))
        allocate(now%He_lith(nx,ny))
        allocate(now%D_lith(nx,ny))

!         allocate(now%we(nx,ny))
!         allocate(now%charge(nx,ny))
!         allocate(now%w0(nx,ny))
!         allocate(now%w1(nx,ny))
        
!         allocate(now%we(-nrad:nrad,-nrad:nrad))
!         allocate(now%charge(1-nrad:nx+nrad,1-nrad:ny+nrad))
        allocate(now%we(nrad*2+1,nrad*2+1))
        allocate(now%charge(1:nx+2*nrad,1:ny+2*nrad))

        allocate(now%w0(nx,ny))
        allocate(now%w1(nx,ny))
        
        return 

    end subroutine isos_allocate

    subroutine isos_deallocate(now)

        implicit none 

        type(isos_state_class), intent(INOUT) :: now 

        if (allocated(now%z_bed))       deallocate(now%z_bed)
        if (allocated(now%dzbdt))       deallocate(now%dzbdt)
        if (allocated(now%z_bed_ref))   deallocate(now%z_bed_ref)
        
        if (allocated(now%tau))         deallocate(now%tau)
        if (allocated(now%He_lith))     deallocate(now%He_lith)
        if (allocated(now%D_lith))      deallocate(now%D_lith)
        
        if (allocated(now%we))          deallocate(now%we)
        if (allocated(now%charge))      deallocate(now%charge)
        if (allocated(now%w0))          deallocate(now%w0)
        if (allocated(now%w1))          deallocate(now%w1)
        
        
        return 

    end subroutine isos_deallocate


    ! === ISOS physics routines ======================================

    subroutine init_charge(charge,H0,BSOC0,sealevel,lbloc)
        !******** initialisation de CHARGE ***********
        ! pour calcul du socle initial on calcule la charge
        ! avec l'etat initial suppose en equilibre S0, H0, Bsoc0

        implicit none 

        real(wp), intent(INOUT) :: charge(:,:)
        real(wp), intent(IN) :: H0(:,:)
        real(wp), intent(IN) :: BSOC0(:,:)
        real(wp), intent(IN) :: sealevel(:,:)
        integer, intent(IN) :: lbloc 

        integer :: i, j, nx, ny, nrad 
        real(wp), allocatable :: charge_local(:,:)

        ! Size of neighborhood 
        nrad = lbloc

        nx = size(BSOC0,1)
        ny = size(BSOC0,2)

        
        allocate(charge_local(1-nrad:nx+nrad,1-nrad:ny+nrad))
        charge_local(1-nrad:nx+nrad,1-nrad:ny+nrad) = charge 

        
        do i = 1, nx
        do j = 1, ny
            if (rho_ice*H0(i,j).ge.rho_ocean*(sealevel(i,j)-BSOC0(i,j))) then
                ! glace ou terre
                
                charge_local(i,j) = (rho_ice*G)*H0(i,j)

            else
                ! ocean
                
                charge_local(i,j) = (rho_ocean*G)*(sealevel(i,j)-BSOC0(i,j))

            endif

        end do
        end do

        do j = 1, ny  ! parties de charge_local a l'exterieure de la grille
            charge_local(1-lbloc:0,j)     = charge_local(1,j)
            charge_local(nx+1:nx+lbloc,j) = charge_local(nx,j)
        end do
        do i = 1, nx
            charge_local(i,1-lbloc:0)     = charge_local(i,1)
            charge_local(i,ny+1:ny+lbloc) = charge_local(i,ny)
        end do

        do j = 1, ny  ! parties de charge_local a l'exterieure de la grille
            charge_local(1-lbloc:0,j)     = charge_local(1,j)
            charge_local(nx+1:nx+lbloc,j) = charge_local(nx,j)
        end do
        do i = 1, nx
            charge_local(i,1-lbloc:0)     = charge_local(i,1)
            charge_local(i,ny+1:ny+lbloc) = charge_local(i,ny)
        end do

        charge_local(1-lbloc:0,1-lbloc:0)         = charge_local(1,1) !valeurs aux quatres coins
        charge_local(1-lbloc:0,ny+1:ny+lbloc)     = charge_local(1,ny) ! exterieurs au domaine
        charge_local(nx+1:nx+lbloc,1-lbloc:0)     = charge_local(nx,1)
        charge_local(nx+1:nx+lbloc,ny+1:ny+lbloc) = charge_local(nx,ny)
        
        ! Return charge_local to the external variable (to match indices from 1:nx+2*nrad)
        charge = charge_local(1-nrad:nx+nrad,1-nrad:ny+nrad)

        return 

    end subroutine init_charge

    subroutine read_tab_litho(WE,filename,RL,DL,DX,DY) 
        !        subroutine qui donne la repartition des enfoncements
        !        en fonction de la rigidite de la lithosphere.
        !        definition du tableau
        !
        !   variables en entree-----------
        !      ROM masse volumique du manteau
        !      RL  rayon de rigidite relative
        !      DL  rigidite flexural
        !  
        !   variables en sortie------------
        !      WE  deflection due a une charge unitaire      
        !
        !
        !

        implicit none 

        real(wp), intent(INOUT) :: WE(:,:) 
        character(len=*), intent(IN) :: filename 
        real(wp), intent(IN) :: RL, DL, DX, DY 

        ! Local variables
        integer, parameter :: nk = 1000
        integer :: i, j, k 
        integer :: LBLOC 
        real(wp) :: kei(0:nk)
        real(wp) :: stepk, AL, XL, DIST
        real(wp) :: som
        integer :: num_kelvin = 177 

        real(wp), allocatable :: WE00(:,:)

        ! Size of neighborhood 
        LBLOC = (size(WE,1)-1)/2 

        ! Allocate local WE object with proper indexing
        allocate(WE00(-LBLOC:LBLOC,-LBLOC:LBLOC))

        ! pour la lithosphere
        STEPK = 100.0
        AL    = -RL*RL/(2.0*pi*DL)*DX*DY      

        ! fonction de kelvin
        ! lecture de la table kei qui est tous les 0.01 entre 0 et 10
        ! STEPK=100=1/ecart 
        !cdc modification du chemin maintenant fonction de dir_inp
        ! trim(dir_inp)//'kelvin.res'
        open(num_kelvin,file=trim(filename))
        read(num_kelvin,*)  ! Skip first line
        do K=1,NK
            read(num_kelvin,*) XL,kei(K)
        end do
        close(num_kelvin)

        do i = -LBLOC, LBLOC    
        do j = -LBLOC, LBLOC
            DIST = dx*sqrt(1.0*(i*i+j*j))                                
            XL   = DIST/RL*STEPK                                          
            K    = int(XL)
            if ((K.gt.834).or.(DIST.gt.dx*LBLOC)) then                 
                WE00(i,j) = 0.0                                              
            else 
                WE00(i,j)=kei(k)+(kei(k+1)-kei(k))*(XL-K)
                if (K.eq.834) WE00(i,j) = max(WE00(i,j),0.0)                   
            endif
            WE00(i,j) = WE00(i,j)*AL 
        end do
        end do

        
        ! normalisation
        som = SUM(WE00)
        WE00  = WE00/(som*rho_mantle*G)

        ! Return solution to external object
        WE = WE00(-LBLOC:LBLOC,-LBLOC:LBLOC)

        ! Make sure too small values are eliminated 
        where(abs(WE) .lt. 1e-12) WE = 0.0 

        return
    
    end subroutine read_tab_litho

    subroutine litho(W1,CHARGE,WE,LBLOC)
        ! litho-0.3.f            10 Novembre 1999             *     
        !
        ! Petit routine qui donne la repartition des enfoncements
        ! en fonction de la rigidite de la lithosphere.
        !
        ! En entree 
        !      ------------
        !     WE(-LBLOC:LBLOC,-LBLOC:LBLOC)  : deflection due a une charge unitaire 
        !                          defini dans tab-litho
        !     LBLOC : relie a la distance : distance en noeud autour de laquelle 
        !     la flexure de la lithosphere est calculee
        !
        !     CHARGE(1-LBLOC:NX+LBLOC,1-LBLOC:NY+LBLOC) : poids par unite de surface
        !             (unite ?)        au temps time, calcule avant  'appel a litho 
        !                     dans taubed ou initial2 
         
        implicit none

        real(wp), intent(INOUT) :: W1(:,:)       ! enfoncement courant
        real(wp), intent(IN)    :: CHARGE(:,:)
        real(wp), intent(IN)    :: WE(:,:)
        integer, intent(IN)    :: LBLOC 

        ! Local variables
        integer :: IP,JP,LPX,LPY,II,SOM1,SOM2
        real(wp), allocatable :: WLOC(:,:)
        real(wp), allocatable :: croix(:)

        integer :: i, j, nx ,ny
        integer :: nrad  
        real(wp), allocatable :: charge_local(:,:)

        nx = size(W1,1)
        ny = size(W1,2)

        nrad = LBLOC 
        allocate(charge_local(1-nrad:nx+nrad,1-nrad:ny+nrad))
        charge_local(1-nrad:nx+nrad,1-nrad:ny+nrad) = charge 

        ! ----- allocation de WLOC  et de croix -----------

        if (allocated(WLOC)) deallocate(WLOC)
        allocate(WLOC(-LBLOC:LBLOC,-LBLOC:LBLOC))

        if (allocated(croix)) deallocate(croix)
        allocate(croix(0:LBLOC))

        ! calcul de la deflexion par sommation des voisins appartenant
        ! au bloc de taille 2*LBLOC
        som1 = 0.0
        som2 = 0.0

        ! On somme aussi les contributions des points exterieurs au domaine
        ! lorsque la charge est due a l'ocean. On suppose alors  que
        ! ces points ont la meme charge que les limites

        do j = 1, ny
        do i = 1, nx

            W1(i,j) = 0.0 
            ii      = 0

            ! Apply the neighborhood weighting to the charge
            WLOC = WE * charge_local(I-LBLOC:I+LBLOC,J-LBLOC:J+LBLOC)

            ! sommation de tous les effets (tentative pour
            ! eviter les erreurs d'arrondi)
            W1(i,j)=WLOC(0,0)

            ! dans croix on calcul la somme des effets WLOC puis on met le resultat dans W1   
            do LPX=1,LBLOC
                LPY=0
                CROIX(LPY)=(WLOC(LPX,LPY)+WLOC(-LPX,LPY))

                CROIX(LPY)=(WLOC(LPY,LPX)+WLOC(LPY,-LPX)) + CROIX(LPY) 
                LPY=LPX
                CROIX(LPY)= ((WLOC(LPX,LPY)+WLOC(LPX,-LPY)) &
                             +(WLOC(-LPX,LPY)+WLOC(-LPX,-LPY))) 
                do LPY=1,LPX-1
                    CROIX(LPY)= (((WLOC(LPX,LPY)+WLOC(LPX,-LPY)) &
                                +(WLOC(-LPX,LPY)+WLOC(-LPX,-LPY))) &
                                +((WLOC(LPY,LPX)+WLOC(LPY,-LPX)) &
                                +(WLOC(-LPY,LPX)+WLOC(-LPY,-LPX))))
                end do

                do LPY=0,LPX 
                    W1(i,j)=W1(i,j)+CROIX(LPY) ! sommation dans W1
                end do
     
            end do

            ! --- FIN DE L'INTEGRATION SUR LE PAVE LBLOC
            som1 = som1 + W1(i,j)
            som2 = som2 - charge_local(i,j)/(rho_mantle*G)

        end do
        end do
   
        return

    end subroutine litho
        
    subroutine calc_litho_regional(w1,charge,we,z_bed,H_ice,z_sl)
        ! Previously known as `taubed`
        ! Routine qui calcul la charge en chaque point de la grille
        ! puis appel la routine litho pour calculer la contribution 
        ! de chaque point a la deflexion de la lithosphere
        !    En sortie
        !   ------------
        !       CHARGE(1-LBLOC:NX+LBLOC,1-LBLOC:NY+LBLOC) : poids par unite de surface
        !               (unite ?)   Elle est calculee initialement dans initial2
        !               Poids de la colonne d'eau ou de la colonne de glace.
        !               a l'exterieur du domaine : 1-LBLOC:1 et NX+1:NX+LBLOC
        !               on donne les valeurs des bords de la grille
        !               CHARGE est utilise par litho uniquement
        !
        !       W1(NX,NY) est l'enfoncement courant, c'est le resultat 
        !               de la routine litho
        !               W1 peut etre calcule de plusieurs facons selon le modele 
        !               d'isostasie utilise
        !

        implicit none

        real(wp), intent(INOUT) :: w1(:,:) 
        real(wp), intent(INOUT) :: charge(:,:) 
        real(wp), intent(IN) :: we(:,:) 
        real(wp), intent(IN) :: z_bed(:,:) 
        real(wp), intent(IN) :: H_ice(:,:) 
        real(wp), intent(IN) :: z_sl(:,:) 
        
        ! Local variables 
        integer :: nrad 
        integer :: i, j, nx, ny 
        real(wp), allocatable :: charge_local(:,:)

        nx = size(W1,1)
        ny = size(W1,2)

        ! Size of neighborhood 
        nrad = (size(we,1)-1)/2 

        allocate(charge_local(1-nrad:nx+nrad,1-nrad:ny+nrad))
        charge_local(1-nrad:nx+nrad,1-nrad:ny+nrad) = charge 

        ! ********* calcul de W1 l'enfoncement d'equilibre au temps t
        ! NLITH est defini dans isostasie et permet le choix du modele d'isostasie

        ! avec rigidite de la lithosphere

        do j = 1, ny 
        do i = 1, nx
            if (rho_ice*H_ice(i,j).ge.rho_ocean*(z_sl(i,j)-z_bed(i,j))) then
                ! glace ou terre 
                charge_local(i,j)=(rho_ice*G)*H_ice(i,j)
            else
                ! ocean
                charge_local(i,j)=(rho_ocean*G)*(z_sl(i,j)-z_bed(i,j))
            endif
        end do
        end do


        ! il faut remplir charge_local dans les parties a l'exterieur de la grille :
        ! a l'exterieur de la grille charge_local est egale a la valeur sur le bord

        do j = 1, ny
            charge_local(1-nrad:0,j)=charge_local(1,j)      ! bord W
            charge_local(NX+1:NX+nrad,j)=charge_local(nx,j) ! bord E
        end do
        do i = 1, nx
            charge_local(i,1-nrad:0)=charge_local(i,1)      ! bord S
            charge_local(i,ny+1:ny+nrad)=charge_local(i,ny) ! bord N
        end do

        ! valeur dans les quatres coins exterieurs au domaine       
        charge_local(1-nrad:0,1-nrad:0)         = charge_local(1,1)   ! coin SW
        charge_local(1-nrad:0,ny+1:ny+nrad)     = charge_local(1,ny)  ! coin NW
        charge_local(nx+1:nx+nrad,1-nrad:0)     = charge_local(nx,1)  ! coin SE
        charge_local(nx+1:nx+nrad,ny+1:ny+nrad) = charge_local(nx,ny) ! coin NE

        call litho(w1,charge_local,we,nrad)

        ! Return charge to the external variable (to match indices from 1:nx+2*nrad)
        charge = charge_local(1-nrad:nx+nrad,1-nrad:ny+nrad)

        return

    end subroutine calc_litho_regional

    elemental subroutine calc_litho_local(w1,z_bed,H_ice,z_sl)
        ! Previously the else-statement in the routine `taubed`

        implicit none 

        real(wp), intent(INOUT) :: w1
        real(wp), intent(IN)    :: z_bed, H_ice, z_sl 

        if (rho_ice*H_ice.ge.rho_ocean*(z_sl-z_bed)) then
            ! Ice or land 
            w1 = rho_ice/rho_mantle*H_ice

        else
            ! Ocean
            w1 = rho_ocean/rho_mantle*(z_sl-z_bed)

        end if

        return 

    end subroutine calc_litho_local

    elemental subroutine calc_uplift_relax(dzbdt,z_bed,z_bed_ref,w_b,tau)

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
