module ismip6
    ! This module contains routines that help with performing the ISMIP6 suite
    ! of experiments. 
    
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
        
        ! Current state:

        ! Atmospheric fields
        type(varslice_class)   :: ts
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
        type(varslice_class)   :: smb_ref

        type(varslice_class)   :: ts_hist 
        type(varslice_class)   :: smb_hist

        type(varslice_class)   :: ts_proj
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

    ! Class for holding ice output for writing to standard formats...
    type ismip6_ice_class

    end type 


    private
    public :: ismip6_forcing_class
    public :: ismip6_ice_class
    public :: ismip6_forcing_init
    public :: ismip6_forcing_update

contains
    

    subroutine ismip6_forcing_init(ism,filename,experiment,domain,grid_name)

        implicit none 

        type(ismip6_forcing_class), intent(INOUT) :: ism
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: experiment 
        character(len=*), intent(IN), optional :: domain 
        character(len=*), intent(IN), optional :: grid_name 

        ! Local variables 
        character(len=256) :: group_prefix 


        select case(trim(experiment))

            case("noresm_rcp85")

                group_prefix = "noresm_rcp85_"

            case DEFAULT 

                write(*,*) "ismip6_forcing_init:: Error: experiment not recognized."
                write(*,*) "experiment = ", trim(experiment) 
                stop 

        end select

        ! Initialize all variables from namelist entries 

        ! General fields 
        call varslice_init_nml(ism%basins,   filename,group="imbie_basins",domain=domain,grid_name=grid_name)
        
        ! Amospheric fields
        call varslice_init_nml(ism%ts_ref,   filename,group=trim(group_prefix)//"ts_ref",domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%smb_ref,  filename,group=trim(group_prefix)//"smb_ref",domain=domain,grid_name=grid_name)
        
        call varslice_init_nml(ism%ts_hist,  filename,group=trim(group_prefix)//"ts_hist",domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%smb_hist, filename,group=trim(group_prefix)//"smb_hist",domain=domain,grid_name=grid_name)

        call varslice_init_nml(ism%ts_proj,  filename,group=trim(group_prefix)//"ts_proj",domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%smb_proj, filename,group=trim(group_prefix)//"smb_proj",domain=domain,grid_name=grid_name)

        ! Oceanic fields
        call varslice_init_nml(ism%to_ref,   filename,group="to_ref",domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%so_ref,   filename,group="so_ref",domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%tf_ref,   filename,group="tf_ref",domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%tf_cor,   filename,group="tf_cor",domain=domain,grid_name=grid_name)

        call varslice_init_nml(ism%to_hist,  filename,group=trim(group_prefix)//"to_hist",domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%so_hist,  filename,group=trim(group_prefix)//"so_hist",domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%tf_hist,  filename,group=trim(group_prefix)//"tf_hist",domain=domain,grid_name=grid_name)

        call varslice_init_nml(ism%to_proj,  filename,group=trim(group_prefix)//"to_proj",domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%so_proj,  filename,group=trim(group_prefix)//"so_proj",domain=domain,grid_name=grid_name)
        call varslice_init_nml(ism%tf_proj,  filename,group=trim(group_prefix)//"tf_proj",domain=domain,grid_name=grid_name)

        ! Load time-independent fields

        ! Amospheric fields 
        call varslice_update(ism%ts_ref)
        call varslice_update(ism%smb_ref)

        ! Oceanic fields
        call varslice_update(ism%to_ref)
        call varslice_update(ism%so_ref)
        call varslice_update(ism%tf_ref)
        call varslice_update(ism%tf_cor)

        return 

    end subroutine ismip6_forcing_init


    subroutine ismip6_forcing_update(ism,time)

        implicit none 

        type(ismip6_forcing_class), intent(INOUT) :: ism
        real(wp), intent(IN) :: time

        ! Local variables 
        integer :: k 

        ! Get slices for current time

        ! === Atmospheric fields ==================================
        
        if (time .lt. 1950) then 

            call varslice_update(ism%ts_hist,1950.0_wp)
            call varslice_update(ism%smb_hist,1950.0_wp)

            ism%ts  = ism%ts_hist 
            ism%smb = ism%smb_hist 
            
        else if (time .ge. 1950 .and. time .le. 1994) then 

            call varslice_update(ism%ts_hist,time)
            call varslice_update(ism%smb_hist,time)

            ism%ts  = ism%ts_hist 
            ism%smb = ism%smb_hist 
            
        else if (time .ge. 1995 .and. time .le. 2100) then 

            call varslice_update(ism%ts_proj,time)
            call varslice_update(ism%smb_proj,time)

            ism%ts  = ism%ts_proj
            ism%smb = ism%smb_proj
            
        else ! time .gt. 2100

            call varslice_update(ism%ts_proj,2100.0_wp)
            call varslice_update(ism%smb_proj,2100.0_wp)

            ism%ts  = ism%ts_proj
            ism%smb = ism%smb_proj
            
        end if

        ! === Oceanic fields ==================================

        if (time .lt. 1850) then 

            ! Oceanic fields 
            call varslice_update(ism%to_hist,1850.0_wp)
            call varslice_update(ism%so_hist,1850.0_wp)
            call varslice_update(ism%tf_hist,1850.0_wp)

            ism%to = ism%to_hist
            ism%so = ism%so_hist
            ism%tf = ism%tf_hist
            
        else if (time .ge. 1850 .and. time .le. 1994) then 

            ! Oceanic fields 
            call varslice_update(ism%to_hist,time)
            call varslice_update(ism%so_hist,time)
            call varslice_update(ism%tf_hist,time)

            ism%to = ism%to_hist
            ism%so = ism%so_hist
            ism%tf = ism%tf_hist
            
        else if (time .ge. 1995 .and. time .le. 2100) then 

            ! Oceanic fields 
            call varslice_update(ism%to_proj,time)
            call varslice_update(ism%so_proj,time)
            call varslice_update(ism%tf_proj,time)

            ism%to = ism%to_proj
            ism%so = ism%so_proj
            ism%tf = ism%tf_proj
            
        else ! time .gt. 2100

            ! Oceanic fields 
            call varslice_update(ism%to_proj,2100.0_wp)
            call varslice_update(ism%so_proj,2100.0_wp)
            call varslice_update(ism%tf_proj,2100.0_wp)

            ism%to = ism%to_proj
            ism%so = ism%so_proj
            ism%tf = ism%tf_proj
            
        end if

            

        ! === Additional calculations ======================

        ! Apply oceanic correction factor to each depth level
        do k = 1, size(ism%to%var,3)   
            ism%to%var(:,:,k) = ism%to%var(:,:,k) + ism%tf_cor%var(:,:,1) 
            ism%tf%var(:,:,k) = ism%tf%var(:,:,k) + ism%tf_cor%var(:,:,1) 
        end do 

        return 

    end subroutine ismip6_forcing_update

end module ismip6