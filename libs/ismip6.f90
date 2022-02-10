module ismip6
    ! This module contains routines that help with performing the ISMIP6 suite
    ! of experiments. 
    
    use nml  
    use ncio 
    use varslice

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the working precision of the library (sp,dp)
    integer,  parameter :: wp = sp 

    ! Define default missing value 
    real(wp), parameter :: mv = -9999.0_wp 

    ! Class for holding ice-forcing data from ISMIP6 archives
    type ismip6_forcing_class
        
        ! Experiment information 
        character(len=256)     :: gcm 
        character(len=256)     :: scen 
        character(len=256)     :: experiment 
        
        ! Current state:

        ! Atmospheric fields
        type(varslice_class)   :: ts
        type(varslice_class)   :: pr
        type(varslice_class)   :: smb

        ! Oceanic fields 
        type(varslice_class)   :: to
        type(varslice_class)   :: so
        type(varslice_class)   :: tf

        ! Resources: 

        ! General fields 
        type(varslice_class)   :: basins

        ! Atmospheric fields
        type(varslice_class)   :: ts_ref 
        type(varslice_class)   :: pr_ref 
        type(varslice_class)   :: smb_ref

        type(varslice_class)   :: ts_hist 
        type(varslice_class)   :: pr_hist 
        type(varslice_class)   :: smb_hist

        type(varslice_class)   :: ts_proj
        type(varslice_class)   :: pr_proj
        type(varslice_class)   :: smb_proj


        ! Oceanic fields 
        type(varslice_class)   :: to_ref
        type(varslice_class)   :: so_ref
        type(varslice_class)   :: tf_ref
        type(varslice_class)   :: tf_cor

        type(varslice_class)   :: to_hist
        type(varslice_class)   :: so_hist
        type(varslice_class)   :: tf_hist

        type(varslice_class)   :: to_proj
        type(varslice_class)   :: so_proj
        type(varslice_class)   :: tf_proj
        
    end type

    type ismip6_ice_var_class
        character(len=56)  :: name 
        character(len=128) :: long_name
        character(len=12)  :: var_type
        character(len=128) :: standard_name 
        character(len=128) :: units_in
        character(len=128) :: units_out
        real(wp) :: unit_scale 
        real(wp) :: unit_offset
    end type

        
    ! Class for holding ice output for writing to standard formats...
    type ismip6_ice_class
        type(ismip6_ice_var_class), allocatable :: vars(:)
    end type 


    private
    public :: ismip6_forcing_class
    public :: ismip6_ice_class
    public :: ismip6_forcing_init
    public :: ismip6_forcing_update

    public :: ismip6_write_init
    public :: ismip6_write_step

contains
    

    subroutine ismip6_forcing_init(ism,filename,gcm,scen,domain,grid_name)

        implicit none 

        type(ismip6_forcing_class), intent(INOUT) :: ism
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: gcm
        character(len=*), intent(IN) :: scen
        character(len=*), intent(IN), optional :: domain 
        character(len=*), intent(IN), optional :: grid_name 

        ! Local variables 
        character(len=256) :: group_prefix 
        character(len=256) :: grp_ts_ref 
        character(len=256) :: grp_pr_ref 
        character(len=256) :: grp_smb_ref 
        character(len=256) :: grp_ts_hist 
        character(len=256) :: grp_pr_hist 
        character(len=256) :: grp_smb_hist 
        character(len=256) :: grp_ts_proj 
        character(len=256) :: grp_pr_proj 
        character(len=256) :: grp_smb_proj 
        
        character(len=256) :: grp_to_ref 
        character(len=256) :: grp_so_ref 
        character(len=256) :: grp_tf_ref 
        character(len=256) :: grp_tf_cor
        character(len=256) :: grp_to_hist 
        character(len=256) :: grp_so_hist 
        character(len=256) :: grp_tf_hist 
        character(len=256) :: grp_to_proj 
        character(len=256) :: grp_so_proj 
        character(len=256) :: grp_tf_proj 
        
        integer  :: k 
        real(wp) :: tmp 
        
        ! Define the current experiment characteristics
        ism%gcm        = trim(gcm)
        ism%scen       = trim(scen) 
        ism%experiment = trim(ism%gcm)//"_"//trim(ism%scen) 

        select case(trim(ism%experiment))

            case("noresm_rcp85","noresm_ctrl")
                ! Control and RCP85 scenarios use the same files 
                ! since ctrl specific forcing adapted in update step 

                group_prefix = "noresm_rcp85_"

                grp_ts_ref   = trim(group_prefix)//"ts_ref"
                grp_pr_ref   = trim(group_prefix)//"pr_ref"
                grp_smb_ref  = trim(group_prefix)//"smb_ref"
                grp_ts_hist  = trim(group_prefix)//"ts_hist"
                grp_pr_hist  = trim(group_prefix)//"pr_hist"
                grp_smb_hist = trim(group_prefix)//"smb_hist"
                grp_ts_proj  = trim(group_prefix)//"ts_proj"
                grp_pr_proj  = trim(group_prefix)//"pr_proj"
                grp_smb_proj = trim(group_prefix)//"smb_proj"
                
                grp_to_ref   = "to_ref"
                grp_so_ref   = "so_ref"
                grp_tf_ref   = "tf_ref"
                grp_tf_cor   = "tf_cor"
                grp_to_hist  = trim(group_prefix)//"to_hist"
                grp_so_hist  = trim(group_prefix)//"so_hist"
                grp_tf_hist  = trim(group_prefix)//"tf_hist"
                grp_to_proj  = trim(group_prefix)//"to_proj"
                grp_so_proj  = trim(group_prefix)//"so_proj"
                grp_tf_proj  = trim(group_prefix)//"tf_proj"
                
            case DEFAULT 

                write(*,*) "ismip6_forcing_init:: Error: experiment not recognized."
                write(*,*) "experiment = ", trim(ism%experiment) 
                stop 

        end select

        ! Initialize all variables from namelist entries 

        ! General fields 
        call varslice_init_nml(ism%basins,   filename,group="imbie_basins",domain=domain,grid_name=grid_name)
        
        ! Amospheric fields
        call varslice_init_nml(ism%ts_ref,   filename,group=trim(grp_ts_ref), domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%pr_ref,   filename,group=trim(grp_pr_ref), domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%smb_ref,  filename,group=trim(grp_smb_ref),domain=domain,grid_name=grid_name)
        
        call varslice_init_nml(ism%ts_hist,  filename,group=trim(grp_ts_hist), domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%pr_hist,  filename,group=trim(grp_pr_hist), domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%smb_hist, filename,group=trim(grp_smb_hist),domain=domain,grid_name=grid_name)

        call varslice_init_nml(ism%ts_proj,  filename,group=trim(grp_ts_proj), domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%pr_proj,  filename,group=trim(grp_pr_proj), domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%smb_proj, filename,group=trim(grp_smb_proj),domain=domain,grid_name=grid_name)

        ! Oceanic fields
        call varslice_init_nml(ism%to_ref,   filename,group=trim(grp_to_ref),domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%so_ref,   filename,group=trim(grp_so_ref),domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%tf_ref,   filename,group=trim(grp_tf_ref),domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%tf_cor,   filename,group=trim(grp_tf_cor),domain=domain,grid_name=grid_name)

        call varslice_init_nml(ism%to_hist,  filename,group=trim(grp_to_hist),domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%so_hist,  filename,group=trim(grp_so_hist),domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%tf_hist,  filename,group=trim(grp_tf_hist),domain=domain,grid_name=grid_name)

        call varslice_init_nml(ism%to_proj,  filename,group=trim(grp_to_proj),domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%so_proj,  filename,group=trim(grp_so_proj),domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%tf_proj,  filename,group=trim(grp_tf_proj),domain=domain,grid_name=grid_name)

        ! Load time-independent fields

        ! Amospheric fields 
        call varslice_update(ism%ts_ref)
        call varslice_update(ism%pr_ref)
        call varslice_update(ism%smb_ref)

        ! Oceanic fields
        call varslice_update(ism%to_ref)
        call varslice_update(ism%so_ref)
        call varslice_update(ism%tf_ref)
        call varslice_update(ism%tf_cor)

        ! Remove missing values from the ocean, if possible
        do k = 1, size(ism%to_ref%var,3)
            if (count(ism%to_ref%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = minval(ism%to_ref%var(:,:,k,1),mask=ism%to_ref%var(:,:,k,1) .ne. mv)
                where(ism%to_ref%var(:,:,k,1) .eq. mv) 
                    ism%to_ref%var(:,:,k,1) = tmp
                end where 
            end if
            if (count(ism%so_ref%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = maxval(ism%so_ref%var(:,:,k,1),mask=ism%so_ref%var(:,:,k,1) .ne. mv)
                where(ism%so_ref%var(:,:,k,1) .eq. mv) 
                    ism%so_ref%var(:,:,k,1) = tmp
                end where 
            end if
            if (count(ism%tf_ref%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = maxval(ism%tf_ref%var(:,:,k,1),mask=ism%tf_ref%var(:,:,k,1) .ne. mv)
                where(ism%tf_ref%var(:,:,k,1) .eq. mv) 
                    ism%tf_ref%var(:,:,k,1) = tmp
                end where 
            end if
        end do

        return 

    end subroutine ismip6_forcing_init


    subroutine ismip6_forcing_update(ism,time,use_ref_atm,use_ref_ocn)

        implicit none 

        type(ismip6_forcing_class), intent(INOUT) :: ism
        real(wp), intent(IN) :: time
        logical,  intent(IN), optional :: use_ref_atm 
        logical,  intent(IN), optional :: use_ref_ocn 

        ! Local variables 
        integer  :: k 
        real(wp) :: tmp 
        character(len=56) :: slice_method 

        ! Get slices for current time

        slice_method = "interp" 

        ! === Atmospheric fields ==================================
        
        if (time .lt. 1995) then 

            ! call varslice_update(ism%ts_hist, [1950.0_wp,1980.0_wp],method="range_mean")
            ! call varslice_update(ism%pr_hist, [1950.0_wp,1980.0_wp],method="range_mean")
            ! call varslice_update(ism%smb_hist,[1950.0_wp,1980.0_wp],method="range_mean")

            ! ism%ts  = ism%ts_hist 
            ! ism%pr  = ism%pr_hist 
            ! ism%smb = ism%smb_hist
            
            ! Prior to 1995, no anomalies 
            ! Set atm fields to reference values (ie, zero anomaly) 
            ism%ts  = ism%ts_ref 
            ism%pr  = ism%pr_ref 
            ism%smb = ism%smb_ref 

            ! Since atm fields are anomalies, also set actual variable to zero 
            ism%ts%var  = 0.0_wp 
            ism%pr%var  = 0.0_wp 
            ism%smb%var = 0.0_wp 
            
        else if (time .ge. 1995 .and. time .le. 2014) then 
            
            if ( (trim(ism%gcm) .eq. "noresm" .and. time .ge. 1995) &
                    .or. time .ge. 2005 ) then 
                ! noresm hist file only goes until 1994, other gcms hist file goes to 2004
                call varslice_update(ism%ts_proj, [time],method=slice_method)
                call varslice_update(ism%pr_proj, [time],method=slice_method)
                call varslice_update(ism%smb_proj,[time],method=slice_method)
                
                ism%ts  = ism%ts_proj 
                ism%pr  = ism%pr_proj 
                ism%smb = ism%smb_proj 
                
            else 
                ! Load hist variable as normal 
                call varslice_update(ism%ts_hist, [time],method=slice_method)
                call varslice_update(ism%pr_hist, [time],method=slice_method)
                call varslice_update(ism%smb_hist,[time],method=slice_method)
                
                ism%ts  = ism%ts_hist 
                ism%pr  = ism%pr_hist 
                ism%smb = ism%smb_hist 
                
            end if 

                
        else if (time .ge. 2015 .and. time .le. 2100) then 

            call varslice_update(ism%ts_proj, [time],method=slice_method)
            call varslice_update(ism%pr_proj, [time],method=slice_method)
            call varslice_update(ism%smb_proj,[time],method=slice_method) 

            ism%ts  = ism%ts_proj
            ism%pr  = ism%pr_proj
            ism%smb = ism%smb_proj
            
        else ! time .gt. 2100

            call varslice_update(ism%ts_proj, [2090.0_wp,2100.0_wp],method="range_mean")
            call varslice_update(ism%pr_proj, [2090.0_wp,2100.0_wp],method="range_mean")
            call varslice_update(ism%smb_proj,[2090.0_wp,2100.0_wp],method="range_mean")

            ism%ts  = ism%ts_proj
            ism%pr  = ism%pr_proj
            ism%smb = ism%smb_proj
            
        end if

        ! === Oceanic fields ==================================

        if (time .lt. 1995) then 
            ! Prehistoric 

            ! Oceanic fields 
            ! call varslice_update(ism%to_hist,[1950.0_wp,1980.0_wp],method="range_mean")
            ! call varslice_update(ism%so_hist,[1950.0_wp,1980.0_wp],method="range_mean")
            ! call varslice_update(ism%tf_hist,[1950.0_wp,1980.0_wp],method="range_mean")

            ! ism%to = ism%to_hist
            ! ism%so = ism%so_hist
            ! ism%tf = ism%tf_hist

            ! Set reference oceanic fields
            ism%to = ism%to_ref
            ism%so = ism%so_ref
            ism%tf = ism%tf_ref

        else if (time .ge. 1995 .and. time .le. 2014) then 
            ! Historical period 

            if ( (trim(ism%gcm) .eq. "noresm" .and. time .ge. 1995) &
                    .or. time .ge. 2005) then 
                ! noresm hist file only goes until 1994, other hist files go until 2004
                call varslice_update(ism%to_proj, [time],method=slice_method)
                call varslice_update(ism%so_proj, [time],method=slice_method)
                call varslice_update(ism%tf_proj, [time],method=slice_method)
                
                ism%to = ism%to_proj
                ism%so = ism%so_proj
                ism%tf = ism%tf_proj
                
            else 
                ! Load hist variable as normal 
                call varslice_update(ism%to_hist, [time],method=slice_method)
                call varslice_update(ism%so_hist, [time],method=slice_method)
                call varslice_update(ism%tf_hist, [time],method=slice_method)
                
                ism%to = ism%to_hist
                ism%so = ism%so_hist
                ism%tf = ism%tf_hist
                
            end if 
            
                
        else if (time .ge. 2015 .and. time .le. 2100) then 
            ! Projection period 1 

            call varslice_update(ism%to_proj,[time],method=slice_method)
            call varslice_update(ism%so_proj,[time],method=slice_method)
            call varslice_update(ism%tf_proj,[time],method=slice_method)

            ism%to = ism%to_proj
            ism%so = ism%so_proj
            ism%tf = ism%tf_proj
               
        else ! time .gt. 2100
            ! Projection period 2 

            call varslice_update(ism%to_proj,[2090.0_wp,2100.0_wp],method="range_mean")
            call varslice_update(ism%so_proj,[2090.0_wp,2100.0_wp],method="range_mean")
            call varslice_update(ism%tf_proj,[2090.0_wp,2100.0_wp],method="range_mean")

            ism%to = ism%to_proj
            ism%so = ism%so_proj
            ism%tf = ism%tf_proj 
            
        end if

        ! === Additional calculations ======================

        if (trim(ism%scen) .eq. "ctrl") then 
            ! For control scenario, override above choices and 
            ! set control atm and ocn 

            ! Set atm fields to reference values (ie, zero anomaly) 
            ism%ts  = ism%ts_ref 
            ism%pr  = ism%pr_ref 
            ism%smb = ism%smb_ref 

            ! Since atm fields are anomalies, also set actual variable to zero 
            ism%ts%var  = 0.0_wp 
            ism%pr%var  = 0.0_wp 
            ism%smb%var = 0.0_wp 
            
            ! Set ocn fields to reference values 
            ism%to = ism%to_ref 
            ism%so = ism%so_ref 
            ism%tf = ism%tf_ref 

        end if 
        
        ! If desired, override other choices and use reference atm fields
        if (present(use_ref_atm)) then 
        if (use_ref_atm) then  

            ism%ts  = ism%ts_ref 
            ism%pr  = ism%pr_ref 
            ism%smb = ism%smb_ref 

            ! Since atm fields are anomalies, also set actual variable to zero 
            ism%ts%var  = 0.0_wp 
            ism%pr%var  = 0.0_wp 
            ism%smb%var = 0.0_wp 
            
        end if 
        end if

        ! If desired, override other choices and use reference ocean fields
        if (present(use_ref_ocn)) then 
        if (use_ref_ocn) then  

            ism%to = ism%to_ref 
            ism%so = ism%so_ref 
            ism%tf = ism%tf_ref 

        end if 
        end if

        ! Remove missing values from the ocean, if possible
        do k = 1, size(ism%to%var,3)
            if (count(ism%to%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = minval(ism%to%var(:,:,k,1),mask=ism%to%var(:,:,k,1) .ne. mv)
                where(ism%to%var(:,:,k,1) .eq. mv) 
                    ism%to%var(:,:,k,1) = tmp
                end where 
            end if
            if (count(ism%so%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = maxval(ism%so%var(:,:,k,1),mask=ism%so%var(:,:,k,1) .ne. mv)
                where(ism%so%var(:,:,k,1) .eq. mv) 
                    ism%so%var(:,:,k,1) = tmp
                end where 
            end if
            if (count(ism%tf%var(:,:,k,1) .ne. mv) .gt. 0) then
                tmp = maxval(ism%tf%var(:,:,k,1),mask=ism%tf%var(:,:,k,1) .ne. mv)
                where(ism%tf%var(:,:,k,1) .eq. mv) 
                    ism%tf%var(:,:,k,1) = tmp
                end where 
            end if
        end do

        ! Apply oceanic correction factor to each depth level
        do k = 1, size(ism%to%var,3)
            where(ism%to%var(:,:,k,1) .ne. mv .and. ism%tf_cor%var(:,:,1,1) .ne. mv)   
                ism%to%var(:,:,k,1) = ism%to%var(:,:,k,1) + ism%tf_cor%var(:,:,1,1)
            end where
            where(ism%tf%var(:,:,k,1) .ne. mv .and. ism%tf_cor%var(:,:,1,1) .ne. mv)
                ism%tf%var(:,:,k,1) = ism%tf%var(:,:,k,1) + ism%tf_cor%var(:,:,1,1) 
            end where
        end do 
        
        return 

    end subroutine ismip6_forcing_update

    ! === ISMIP6 OUTPUT ROUTINES ==========

    subroutine ismip6_write_init(filename,xc,yc,time,lon,lat,area,map_name,lambda,phi)

        implicit none 

        character(len=*),   intent(IN) :: filename
        real(wp),           intent(IN) :: xc(:)
        real(wp),           intent(IN) :: yc(:)
        real(wp),           intent(IN) :: time
        real(wp),           intent(IN) :: lon(:,:)
        real(wp),           intent(IN) :: lat(:,:)
        real(wp),           intent(IN) :: area(:,:)
        character(len=*),   intent(IN) :: map_name
        real(wp),           intent(IN) :: lambda
        real(wp),           intent(IN) :: phi 

        ! Local variables 
        character(len=12) :: xnm 
        character(len=12) :: ynm 
        
        xnm = "xc"
        ynm = "yc" 

        ! === Initialize netcdf file and dimensions =========

        ! Create the netcdf file 
        call nc_create(filename)

        ! Add grid axis variables to netcdf file
        call nc_write_dim(filename,xnm,x=xc*1e-3,units="kilometers")
        call nc_write_attr(filename,xnm,"_CoordinateAxisType","GeoX")

        call nc_write_dim(filename,ynm,x=yc*1e-3,units="kilometers")
        call nc_write_attr(filename,ynm,"_CoordinateAxisType","GeoY")
        
        ! Add time axis with current value 
        call nc_write_dim(filename,"time", x=time,dx=1.0_wp,nx=1,units="years",unlimited=.TRUE.)
        
        ! Projection information 
        call nc_write_map(filename,map_name,dble(lambda),phi=dble(phi))

        ! Lat-lon information
        call nc_write(filename,"lon2D",lon,dim1=xnm,dim2=ynm,grid_mapping=map_name)
        call nc_write_attr(filename,"lon2D","_CoordinateAxisType","Lon")
        call nc_write(filename,"lat2D",lat,dim1=xnm,dim2=ynm,grid_mapping=map_name)
        call nc_write_attr(filename,"lat2D","_CoordinateAxisType","Lat")

        call nc_write(filename,"area",  area*1e-6,  dim1=xnm,dim2=ynm,grid_mapping=map_name,units="km^2")
        call nc_write_attr(filename,"area","coordinates","lat2D lon2D")
        

        return

    end subroutine ismip6_write_init


    subroutine ismip6_write_step(filename,file_nml,time)

        implicit none 

        character(len=*),   intent(IN) :: filename
        character(len=*),   intent(IN) :: file_nml 
        real(wp),           intent(IN) :: time 

        ! Local variables 
        integer    :: ncid, n
        real(wp) :: time_prev 
        type(ismip6_ice_class) :: ismp 

if (.FALSE.) then    
        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

end if 
        ! Load up output variable meta information 
        call ismip6_load_ice_var_info(ismp,file_nml,verbose=.TRUE.)

        ! Proceed to calculating and writing each variable 

        ! Variables to output every five years 
        if (time - time_prev .ge. 5.0_wp) then 
            ! Five years have passed, proceed 




        end if 


        ! Next, variables that should be output every year 
        return

    end subroutine ismip6_write_step

    subroutine ismip6_load_ice_var_info(ismp,filename,verbose)

        implicit none 

        type(ismip6_ice_class), intent(OUT) :: ismp 
        character(len=*),       intent(IN)  :: filename 
        logical,                intent(IN)  :: verbose 

        ! Local variables
        integer :: n  
        integer, parameter :: n_variables = 38

        type(ismip6_ice_var_class) :: v 

        ! First initialize ismp object to hold variable meta data 
        if (allocated(ismp%vars)) deallocate(ismp%vars)
        allocate(ismp%vars(n_variables))

        ! Load individual variables by namelist group 
        call ice_var_par_load(ismp%vars(1), filename,var_name="lithk")
        call ice_var_par_load(ismp%vars(2), filename,var_name="orog")
        call ice_var_par_load(ismp%vars(3), filename,var_name="base")
        call ice_var_par_load(ismp%vars(4), filename,var_name="topg")
        call ice_var_par_load(ismp%vars(5), filename,var_name="hfgeoubed")
        call ice_var_par_load(ismp%vars(6), filename,var_name="acabf")
        call ice_var_par_load(ismp%vars(7), filename,var_name="libmassbfgr")
        call ice_var_par_load(ismp%vars(8), filename,var_name="libmassbffl")
        call ice_var_par_load(ismp%vars(9), filename,var_name="dlithkdt")
        call ice_var_par_load(ismp%vars(10),filename,var_name="xvelsurf")
        call ice_var_par_load(ismp%vars(11),filename,var_name="yvelsurf")
        call ice_var_par_load(ismp%vars(12),filename,var_name="zvelsurf")
        call ice_var_par_load(ismp%vars(13),filename,var_name="xvelbase")
        call ice_var_par_load(ismp%vars(14),filename,var_name="yvelbase")
        call ice_var_par_load(ismp%vars(15),filename,var_name="zvelbase")
        call ice_var_par_load(ismp%vars(16),filename,var_name="xvelmean")
        call ice_var_par_load(ismp%vars(17),filename,var_name="yvelmean")
        call ice_var_par_load(ismp%vars(18),filename,var_name="litemptop")
        call ice_var_par_load(ismp%vars(19),filename,var_name="litempbotgr")
        call ice_var_par_load(ismp%vars(20),filename,var_name="litempbotfl")
        call ice_var_par_load(ismp%vars(21),filename,var_name="strbasemag")
        call ice_var_par_load(ismp%vars(22),filename,var_name="licalvf")
        call ice_var_par_load(ismp%vars(23),filename,var_name="lifmassbf")
        call ice_var_par_load(ismp%vars(24),filename,var_name="lifmassbf")
        call ice_var_par_load(ismp%vars(25),filename,var_name="ligroundf")
        call ice_var_par_load(ismp%vars(26),filename,var_name="sftgif")
        call ice_var_par_load(ismp%vars(27),filename,var_name="sftgrf")
        call ice_var_par_load(ismp%vars(28),filename,var_name="sftflf")

        call ice_var_par_load(ismp%vars(29),filename,var_name="lim")
        call ice_var_par_load(ismp%vars(30),filename,var_name="limnsw")
        call ice_var_par_load(ismp%vars(31),filename,var_name="iareagr")
        call ice_var_par_load(ismp%vars(32),filename,var_name="iareafl")
        call ice_var_par_load(ismp%vars(33),filename,var_name="tendacabf")
        call ice_var_par_load(ismp%vars(34),filename,var_name="tendlibmassbf")
        call ice_var_par_load(ismp%vars(35),filename,var_name="tendlibmassbffl")
        call ice_var_par_load(ismp%vars(36),filename,var_name="tendlicalvf")
        call ice_var_par_load(ismp%vars(37),filename,var_name="tendlifmassbf")
        call ice_var_par_load(ismp%vars(38),filename,var_name="tendligroundf")

        if (verbose) then 

            ! === Print summary =========

            write(*,"(a40,a8,a65,a15)") &
                                    "Variable name",    &
                                    "Type",             &
                                    "Standard name",    &
                                    "Unit"
            
            do n = 1, n_variables 
                v = ismp%vars(n)
                write(*,"(a40,a8,a65,a15)") &
                                    trim(v%long_name),      &
                                    trim(v%var_type),       &
                                    trim(v%standard_name),  &
                                    trim(v%units_out) 
            end do 


        end if 

        return 

    end subroutine ismip6_load_ice_var_info

    subroutine ice_var_par_load(ismpv,filename,var_name)
        ! Load parmaeters associated with a given ice variable

        implicit none 

        type(ismip6_ice_var_class), intent(OUT) :: ismpv 
        character(len=*),           intent(IN)  :: filename 
        character(len=*),           intent(IN)  :: var_name

        ! Local variables 
        character(len=56) :: group 

        group = "ismip6_out_"//trim(var_name)

        call nml_read(filename,group,"name",            ismpv%name)
        call nml_read(filename,group,"long_name",       ismpv%long_name)
        call nml_read(filename,group,"var_type",        ismpv%var_type)
        call nml_read(filename,group,"standard_name",   ismpv%standard_name)
        call nml_read(filename,group,"units_in",        ismpv%units_in)
        call nml_read(filename,group,"units_out",       ismpv%units_out)
        call nml_read(filename,group,"unit_scale",      ismpv%unit_scale)
        call nml_read(filename,group,"unit_offset",     ismpv%unit_offset)

        return 

    end subroutine ice_var_par_load

end module ismip6