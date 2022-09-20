module varslice

    use ncio 
    use nml 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the working precision of the library (sp,dp)
    integer,  parameter :: wp = sp 

    ! Define default missing value 
    real(wp), parameter :: mv = -9999.0_wp 

    type varslice_param_class

        character(len=1024) :: filename
        character(len=56)   :: name 
        character(len=56)   :: units_in
        character(len=56)   :: units_out
        real(wp) :: unit_scale 
        real(wp) :: unit_offset
        logical  :: with_time 

        ! Internal parameters
        integer, allocatable :: dim(:)  
        integer  :: ndim 
        real(wp) :: time_par(3) 

    end type

    type varslice_class 

        type(varslice_param_class) :: par 

        ! Parameters defined during update call
        real(wp)          :: time_range(2)
        character(len=56) :: slice_method 
        integer           :: range_rep 
        
        ! Variable information
        real(wp), allocatable :: x(:) 
        real(wp), allocatable :: y(:)
        real(wp), allocatable :: lev(:)  

        real(wp), allocatable :: time(:)  
        real(wp), allocatable :: var(:,:,:,:) 
    end type 

    private 
    public :: varslice_class
    public :: varslice_update
    public :: varslice_init_nml 
    public :: varslice_init_arg
    public :: varslice_end 

contains

    
    subroutine varslice_update(vs,time,method,rep)
        ! Routine to update transient climate forcing to match 
        ! current `time`. 

        implicit none 

        type(varslice_class),       intent(INOUT) :: vs
        real(wp),         optional, intent(IN)    :: time(:)    ! [yr] Current time, or time range 
        character(len=*), optional, intent(IN)    :: method     ! slice_method (only if with_time==True)
        integer,          optional, intent(IN)    :: rep        ! Only if with_time==True, and slice_method==range_*
        
        ! Local variables 
        integer :: k, k0, k1, nt, nt_now, nt_out
        integer :: n1, n2, i, j, l    
        type(varslice_param_class) :: par 
        logical  :: with_time 
        real(wp) :: time_range(2) 
        real(wp) :: time_wt(2)
        character(len=56) :: slice_method
        character(len=56) :: vec_method 
        integer  :: range_rep 
        integer,  allocatable :: kk(:) 
        real(wp), allocatable :: var(:,:,:,:) 

        ! Define shortcuts
        par = vs%par 
        with_time = par%with_time 

        if (with_time) then 

            if (.not. present(time)) then 
                write(*,*) "varslice_update:: Error: current time or time range &
                            &must be given as an argument (1D array)."
                stop 
            end if 

            ! Consistency check 
            if (size(time,1) .eq. 2) then 
                if (time(2) .lt. time(1)) then 
                    write(*,*) "varslice_update:: Error: time(2) should be >= time(1)."
                    write(*,*) "time = ", time
                    stop 
                end if
            end if 

        end if 

        slice_method = "exact"
        if (present(method)) slice_method = trim(method)

        range_rep = 1
        if (present(rep)) range_rep = rep 


        if (present(time)) then

            if (size(time) .eq. 2) then 
                time_range = time
            else 
                time_range(1:2) = time(1) 
            end if 

        else 

            time_range = vs%time_range 

        end if 

        
        if ( with_time .and. trim(slice_method) .eq. trim(vs%slice_method) &
                .and. range_rep .eq. vs%range_rep &
                .and. vs%time_range(1) .eq. time_range(1) &
                .and. vs%time_range(2) .eq. time_range(2) ) then 

            ! Do nothing, the varslice object is already up to date 
            ! fo the current time and method. 

        else 

            ! Set parameters in varslice object 
            vs%slice_method = trim(slice_method)
            vs%range_rep    = range_rep 
            vs%time_range   = time_range 

            ! 2. Read variable and convert units as needed

            if (.not. with_time) then 
                ! Handle cases that do not have time dimension (simpler)

                select case(par%ndim)

                    case(1)

                        ! 1D variable
                        call nc_read(par%filename,par%name,vs%var(:,1,1,1),missing_value=mv)

                    case(2)

                        ! 2D variable 
                        call nc_read(par%filename,par%name,vs%var(:,:,1,1),missing_value=mv)

                    case(3)

                        ! 3D variable 
                        call nc_read(par%filename,par%name,vs%var(:,:,:,1),missing_value=mv)

                    case DEFAULT 

                        write(*,*) "varslice_update:: ndim >= 4 with no time dimension is not allowed."
                        write(*,*) "ndim = ", par%ndim 
                        stop 

                end select

            else 
                ! Cases with a time dimension (more complicated)

                if (trim(slice_method) .eq. "interp" .or. &
                    trim(slice_method) .eq. "extrap") then 
                    ! Update time range for interp/extrap methods 

                    ! Time range should cover all available data since 
                    ! bracketing indices to the desired time will be found
                    vs%time_range(1) = minval(vs%time) 
                    vs%time_range(2) = maxval(vs%time) 

                    ! Additional consistency check 
                    if (size(time,1) .ne. 1) then 
                        write(*,*) "varslice_update:: Error: to use slice_method=['interp','extrap'], &
                        &only one time should be provided as an argument."
                        write(*,*) "time = ", time 
                        stop 
                    end if 

                end if 

                ! Determine indices of data to load 

                call get_indices(k0,k1,vs%time,vs%time_range,trim(slice_method),time(1))



                if (k0 .gt. 0 .and. k1 .gt. 0) then 
                    ! Dimension range is available for loading, proceed 

                    ! Get size of time dimension needed for loading
                    nt_now = k1-k0+1

                    ! write(*,*) "time: ", vs%time_range, k0, k1 
                    ! write(*,*) "      ", vs%time(k0), vs%time(k1)

                    if (allocated(var)) deallocate(var) 

                    select case(par%ndim)

                        case(1)

                            ! Allocate local var to the right size 
                            allocate(var(nt_now,1,1,1))

                            ! 0D (point) variable plus time dimension 
                            call nc_read(par%filename,par%name,var,missing_value=mv, &
                                    start=[k0],count=[nt_now])


                        case(2)

                            ! Allocate local var to the right size 
                            allocate(var(par%dim(1),nt_now,1,1))

                            ! 1D variable plus time dimension 
                            call nc_read(par%filename,par%name,var,missing_value=mv, &
                                    start=[1,k0],count=[par%dim(1),nt_now])

                        case(3)
        
                            ! Allocate local var to the right size 
                            allocate(var(par%dim(1),par%dim(2),nt_now,1))

                            ! 2D variable plus time dimension 
                            call nc_read(par%filename,par%name,var,missing_value=mv, &
                                    start=[1,1,k0],count=[par%dim(1),par%dim(2),nt_now])

                        case(4)

                            ! Allocate local var to the right size 
                            allocate(var(par%dim(1),par%dim(2),par%dim(3),nt_now))

                            ! 3D variable plus time dimension 
                            call nc_read(par%filename,par%name,var,missing_value=mv, &
                                    start=[1,1,1,k0],count=[par%dim(1),par%dim(2),par%dim(3),nt_now]) 

                        case DEFAULT 

                            write(*,*) "varslice_update:: ndim > 4 with time dimension not allowed."
                            write(*,*) "ndim = ", par%ndim 
                            stop 

                    end select


                    ! At this point, the local var variable has been defined 
                    ! with data from the file for the appropriate time indices 

                    ! Next, we need to allocate the vs%var variable to the 
                    ! appropriate size and perform any calculations on the time 
                    ! indices of the local var variable as needed. 

                    ! Handle special case: if only one time is available 
                    ! for interp/extrap methods, then change method 
                    ! to exact 
                    if ( (trim(vs%slice_method) .eq. "interp" .or. & 
                          trim(vs%slice_method) .eq. "extrap") .and. &
                          vs%time(k0) .eq. vs%time(k1)) then 
                        ! Same time is given for upper and lower bound

                        slice_method = "exact"     
                    end if 


                    ! Now, allocate the vs%var variable to the right size 

                    select case(trim(slice_method)) 

                        case("exact","range")
                            ! Allocate vs%var to the same size as var 
                            ! and store all values 

                            if (size(vs%var,1) .eq. size(var,1) .and. &
                                size(vs%var,2) .eq. size(var,2) .and. &
                                size(vs%var,3) .eq. size(var,3) .and. &
                                size(vs%var,4) .eq. size(var,4) ) then 

                                ! No allocation needed size is the same 

                            else 

                                deallocate(vs%var)
                                allocate(vs%var(size(vs%var,1),size(vs%var,2),size(vs%var,3),size(vs%var,4)))

                            end if 

                            ! Store data in vs%var 
                            vs%var = var 

                        case("interp","extrap")
                            ! var should have two time dimensions to interpolate
                            ! between. Allocate vs%var to appropriate size and 
                            ! perform interpolation 

                            if (nt_now .ne. 2) then 
                                write(*,*) "varslice_update:: Error: something went wrong during &
                                &interpolation. More than 2 time slices are available to interpolate &
                                &between. Check!"
                                write(*,*) "nt_now = ", nt_now 
                                write(*,*) "time   = ", time 
                                write(*,*) "indices k0, k1: ", k0, k1 
                                write(*,*) "times : ", vs%time(k0:k1) 
                                stop 
                            end if 

                            ! Calculate time weighting between two extremes
                            time_wt(2) = (time(1)-vs%time(k0)) / (vs%time(k1) - vs%time(k0))
                            time_wt(1) = 1.0_wp - time_wt(2) 

                            if (minval(time_wt) .lt. 0.0_wp .or. maxval(time_wt) .gt. 1.0_wp) then 
                                write(*,*) "varslice_update:: Error: interpolation weights are incorrect."
                                write(*,*) "time_wt  = ", time_wt 
                                write(*,*) "time     = ", time 
                                write(*,*) "time(k0) = ", vs%time(k0)
                                write(*,*) "time(k1) = ", vs%time(k1) 
                                stop
                            end if 
                            
                            ! Note: slice_method='interp' and 'extrap' use the same method, since 
                            ! the indices determine interpolation weights (ie, 
                            ! for slice_method='interp', if the time of interest lies
                            ! outside of the bounds of the data, then k0=k1=-1 and 
                            ! output data are set to missing values) 

                            vec_method = "mean"     ! interp methods use the (weighted) mean
                            nt_out     = 1          ! We expect one time out 

                            deallocate(vs%var)

                            select case(par%ndim)

                                case(1)
                                    allocate(vs%var(nt_out,1,1,1))

                                    ! Calculate the vector value desired (wtd mean)
                                    call calc_vec_value(vs%var(1,1,1,1),var([1,2],1,1,1),vec_method,mv,wt=time_wt)

                                case(2)
                                    allocate(vs%var(size(var,1),nt_out,1,1))

                                    do i = 1, size(vs%var,1)
                                        ! Calculate the vector value desired (mean,sd,...)
                                        call calc_vec_value(vs%var(i,1,1,1),var(i,[1,2],1,1),vec_method,mv,wt=time_wt)
                                    end do 
                                    
                                case(3)
                                    allocate(vs%var(size(var,1),size(var,2),nt_out,1))
                                    
                                    do j = 1, size(vs%var,2)
                                    do i = 1, size(vs%var,1)
                                        ! Calculate the vector value desired (mean,sd,...)
                                        call calc_vec_value(vs%var(i,j,1,1),var(i,j,[1,2],1),vec_method,mv,wt=time_wt)
                                    end do
                                    end do 
                                    
                                case(4)
                                    allocate(vs%var(size(var,1),size(var,2),size(var,3),nt_out))

                                    do l = 1, size(vs%var,3)
                                    do j = 1, size(vs%var,2)
                                    do i = 1, size(vs%var,1)
                                        ! Calculate the vector value desired (mean,sd,...)
                                        call calc_vec_value(vs%var(i,j,l,1),var(i,j,l,[1,2]),vec_method,mv,wt=time_wt)
                                    end do
                                    end do
                                    end do 
                                    
                            end select

                        case("range_mean","range_sd","range_min","range_max","range_sum")
                            ! Allocate vs%var to match desired output size, 
                            ! and calculate output values 

                            ! Define 'vec_method'
                            n1 = index(slice_method,"_")
                            n2 = len_trim(slice_method)
                            vec_method = slice_method(n1+1:n2)

                            ! Size of dimension out is the size of the 
                            ! repitition desired. Ie, range_rep = 1 means 
                            ! to apply mean/sd/etc to all values along dimension
                            ! with the result of calculating 1 value. 
                            ! range_rep = 12 means apply calculation to every 12th 
                            ! value, resulting in 12 values along dimension.

                            nt_out = vs%range_rep 
                            
                            ! Make sure that var has at least as many values as we expect 
                            if (nt_out .gt. nt_now) then 
                                write(*,*) "varslice_update:: Error: the specified time range &
                                    & does not provide enough data points to be consistent with &
                                    & the specified value of range_rep."
                                write(*,*) "time_range      = ", vs%time_range 
                                write(*,*) "nt (time_range) = ", nt_now 
                                write(*,*) "range_rep       = ", vs%range_rep 
                                write(*,*) "range_rep must be <= nt."
                                stop 
                            end if 

                            deallocate(vs%var)

                            select case(par%ndim)

                                case(1)
                                    allocate(vs%var(nt_out,1,1,1))

                                    ! Calculate each slice 
                                    do k = 1, nt_out 

                                        ! Get indices for current repitition
                                        call get_rep_indices(kk,i0=k,i1=nt_now,nrep=vs%range_rep)

                                        ! Calculate the vector value desired (mean,sd,...)
                                        call calc_vec_value(vs%var(k,1,1,1),var(kk,1,1,1),vec_method,mv)

                                    end do 

                                case(2)
                                    allocate(vs%var(size(var,1),nt_out,1,1))

                                    ! Calculate each slice 
                                    do k = 1, nt_out 

                                        ! Get indices for current repitition
                                        call get_rep_indices(kk,i0=k,i1=nt_now,nrep=vs%range_rep)

                                        do i = 1, size(vs%var,1)
                                            ! Calculate the vector value desired (mean,sd,...)
                                            call calc_vec_value(vs%var(i,k,1,1),var(i,kk,1,1),vec_method,mv)
                                        end do 

                                    end do 
                                    
                                case(3)
                                    allocate(vs%var(size(var,1),size(var,2),nt_out,1))
                                
                                    ! Calculate each slice 
                                    do k = 1, nt_out 

                                        ! Get indices for current repitition
                                        call get_rep_indices(kk,i0=k,i1=nt_now,nrep=vs%range_rep)

                                        do j = 1, size(vs%var,2)
                                        do i = 1, size(vs%var,1)
                                            ! Calculate the vector value desired (mean,sd,...)
                                            call calc_vec_value(vs%var(i,j,k,1),var(i,j,kk,1),vec_method,mv)
                                        end do
                                        end do 
                                        
                                    end do 
                                    
                                case(4)
                                    allocate(vs%var(size(var,1),size(var,2),size(var,3),nt_out))

                                    ! Calculate each slice 
                                    do k = 1, nt_out 

                                        ! Get indices for current repitition
                                        call get_rep_indices(kk,i0=k,i1=nt_now,nrep=vs%range_rep)

                                        do l = 1, size(vs%var,3)
                                        do j = 1, size(vs%var,2)
                                        do i = 1, size(vs%var,1)
                                            ! Calculate the vector value desired (mean,sd,...)
                                            call calc_vec_value(vs%var(i,j,l,k),var(i,j,l,kk),vec_method,mv)
                                        end do
                                        end do
                                        end do 
                                        
                                    end do 
                                    
                            end select


                    end select

                else 
                    ! Dimension range was not found, set variable to missing values 

                    vs%var = mv 

                end if

            end if 

            ! Make sure crazy values have been set to missing (for safety)
            where (abs(vs%var) .ge. 1e10) vs%var = mv 

            ! Apply scaling 
            where (vs%var .ne. mv) 
                vs%var = vs%var*par%unit_scale + par%unit_offset
            end where 
            
        end if 

        return 

    end subroutine varslice_update


    subroutine get_indices(k0,k1,x,xrange,slice_method,x_interp)
        ! Get the indices k0 and k1 that 
        ! correspond to the lower and upper bound 
        ! range of xmin <= x <= xmax. 

        ! Resulting indices should either match 
        ! the range exactly, or bracket the range of interest 

        ! Note: routine assumes xmin <= xmax! 

        implicit none 

        integer,  intent(OUT) :: k0 
        integer,  intent(OUT) :: k1 
        real(wp), intent(IN)  :: x(:) 
        real(wp), intent(IN)  :: xrange(2)
        character(len=*), intent(IN) :: slice_method 
        real(wp), intent(IN)  :: x_interp               ! only used for interp methods

        ! Local variables 
        integer  :: k, nk 
        real(wp) :: xmin, xmax 

        xmin = xrange(1)
        xmax = xrange(2) 

        nk = size(x,1) 

        ! Get lower bound 
        k0 = 1 
        do k = 1, nk 
            if (x(k) .gt. xmin) exit 
            k0 = k 
        end do 

        if (xmax .eq. xmin) then 

            k1 = k0 

        else 

            ! Get upper bound 
            k1 = nk 
            do k = nk, k0, -1 
                if (x(k) .lt. xmax) exit 
                k1 = k 
            end do 

        end if 

        ! Make sure indices work for slice_method of choice
        select case(trim(slice_method))

            case("exact") 

                ! If index wasn't found, set idices to -1
                if (x(k0) .ne. xmin) then 
                    k0 = -1 
                    k1 = -1 
                end if 
     
            case("interp") 

                ! If xmin/xmax are not found in range, set indices to -1
                if (x_interp .gt. xmax .or. x_interp .lt. xmin) then 
                    k0 = -1
                    k1 = -1 

                else 

                    ! Redo indices to get nearest bracketing points in time_range

                    ! Get lower bound 
                    k0 = 1 
                    do k = 1, nk 
                        if (x(k) .gt. x_interp) exit 
                        k0 = k 
                    end do 

                    ! Get upper bound 
                    k1 = nk 
                    do k = nk, k0, -1 
                        if (x(k) .lt. x_interp) exit 
                        k1 = k 
                    end do 

                end if 

            case("extrap") 

                ! If xmin/xmax are not found in range, 
                ! set indices to extreme bound
                if (minval(x) .gt. x_interp) then 

                    k0 = 1
                    k1 = 1 
                
                else if (maxval(x) .lt. x_interp) then 

                    k0 = nk
                    k1 = nk

                else 

                    ! Redo indices to get nearest bracketing points in time_range

                    ! Get lower bound 
                    k0 = 1 
                    do k = 1, nk 
                        if (x(k) .gt. x_interp) exit 
                        k0 = k 
                    end do 

                    ! Get upper bound 
                    k1 = nk 
                    do k = nk, k0, -1 
                        if (x(k) .lt. x_interp) exit 
                        k1 = k 
                    end do 

                end if 

            case("range","range_mean","range_sd","range_min","range_max") 
                 
                ! If xmin/xmax are not found in range, set indices to -1
                if (minval(x) .gt. xmax .or. maxval(x) .lt. xmin) then 
                    k0 = -1
                    k1 = -1 
                end if 

            ! No default case
        end select 

        return 

    end subroutine get_indices

    subroutine get_rep_indices(ii,i0,i1,nrep)
        ! Given starting value i0 and final value i1, 
        ! and number of values to skip nrep, generate 
        ! an vector of indices ii. 
        ! eg, i0 = 1, i1 = 36, nrep = 12
        ! => ii = [1,13,25]
        ! eg, i0 = 2, i1 = 36, nrep = 12
        ! => ii = [2,14,26]
        
        implicit none 

        integer, allocatable, intent(OUT) :: ii(:) 
        integer, intent(IN) :: i0 
        integer, intent(IN) :: i1 
        integer, intent(IN) :: nrep 

        ! Local variables   
        integer :: i, ni, ntot 
        integer :: jj(10000)

        ni = i1-i0+1 

        jj = 0  

        do i = 1, ni
            jj(i) = i0 + nrep*(i-1) 
            if (jj(i) .gt. i1) then 
                jj(i) = 0
                exit 
            end if 
        end do 

        ntot = count(jj .gt. 0)
        if (allocated(ii)) deallocate(ii)
        allocate(ii(ntot)) 

        ii = jj(1:ntot) 

        return 

    end subroutine get_rep_indices

    subroutine calc_vec_value(val,var,method,mv,wt)

        implicit none 

        real(wp),         intent(OUT) :: val 
        real(wp),         intent(IN)  :: var(:) 
        character(len=*), intent(IN)  :: method 
        real(wp),         intent(IN)  :: mv 
        real(wp),         intent(IN), optional :: wt(:) 

        ! Local variables 
        integer  :: ntot 
        real(wp) :: mean, variance 
        real(wp), allocatable :: wt_now(:) 

        ntot = count(var .ne. mv)
        
        if (ntot .gt. 0) then 
            ! Values available for calculations 

            allocate(wt_now(size(var)))
            if (present(wt)) then 
                wt_now = wt 
            else 
                wt_now = 1.0_wp
            end if 

            where(var .eq. mv) wt_now = 0.0_wp 

            select case(trim(method))

                case("sum","mean")
                    ! Calculate a weighted sum/mean with the right weights for each case 

                    ! Normalize weights to sum to one
                    if (trim(method) .eq. "mean") wt_now = wt_now / sum(wt_now)
                    
                    val = sum(wt_now*var)

                case("sd")
                    ! Calculate the weighted mean and then weighted standard deviation 

                    ! Normalize weights to sum to one
                    wt_now = wt_now / sum(wt_now)

                    if (ntot .ge. 2) then

                        mean     = sum(wt_now*var)
                        variance = real(ntot/(ntot-1),wp) * sum( wt_now*(var-mean)**2 )
                        val      = sqrt(variance)

                    else 

                        val = mv 

                    end if 

                case("min") 

                    val = minval(var,mask=var.ne.mv)

                case("max")

                    val = maxval(var,mask=var.ne.mv)

                case DEFAULT 

                    write(*,*) "calc_vec_value:: Error: method not recognized."
                    write(*,*) "method = ", trim(method) 
                    stop 

            end select

        else 
            ! No values available in vector, set to missing value 

            val = mv 

        end if 

        return 

    end subroutine calc_vec_value

    ! jablasco: add gcm and scenario
    subroutine varslice_init_nml(vs,filename,group,domain,grid_name,gcm,verbose)
        ! Routine to load information related to a given 
        ! transient variable, so that it can be processed properly.

        implicit none 

        type(varslice_class),   intent(INOUT) :: vs
        character(len=*),       intent(IN)    :: filename
        character(len=*),       intent(IN)    :: group
        character(len=*),       intent(IN), optional :: domain
        character(len=*),       intent(IN), optional :: grid_name
        character(len=*),       intent(IN), optional :: gcm
        logical,                intent(IN), optional :: verbose 
        ! Local variables 
        
        ! First load parameters from nml file 
        call varslice_par_load(vs%par,filename,group,domain,grid_name,gcm,verbose)

        ! Perform remaining init operations 
        call varslice_init_data(vs) 

        return 

    end subroutine varslice_init_nml

    subroutine varslice_init_arg(vs,filename,name,units_in,units_out,scale,offset, &
                                                            with_time,time_par)
        ! Routine to load information related to a given 
        ! transient variable, so that it can be processed properly.

        implicit none 

        type(varslice_class), intent(INOUT) :: vs
        character(len=*),      intent(IN)   :: filename
        character(len=*),      intent(IN)   :: name
        character(len=*),      intent(IN)   :: units_in 
        character(len=*),      intent(IN)   :: units_out 
        real(wp), optional,    intent(IN)   :: scale 
        real(wp), optional,    intent(IN)   :: offset 
        logical,  optional,    intent(IN)   :: with_time 
        real(wp), optional,    intent(IN)   :: time_par(:) 

        ! Local variables 
        
        ! Define parameters from subroutine arguments 
        vs%par%filename  = trim(filename) 
        vs%par%name      = trim(name) 
        vs%par%units_in  = trim(units_in) 
        vs%par%units_out = trim(units_out) 
        
        vs%par%unit_scale = 1.0_wp 
        if (present(scale)) vs%par%unit_scale = scale 

        vs%par%unit_offset = 0.0_wp 
        if (present(offset)) vs%par%unit_offset = offset 
        
        vs%par%with_time = .TRUE. 
        if (present(with_time)) vs%par%with_time = with_time 

        vs%par%time_par = [0.0,0.0,0.0]
        if (present(time_par)) vs%par%time_par = time_par 

        ! Perform remaining init operations 
        call varslice_init_data(vs) 
        
        return 

    end subroutine varslice_init_arg

    subroutine varslice_init_data(vs)

        implicit none 

        type(varslice_class), intent(INOUT) :: vs 

        ! Local variables  
        character(len=12), allocatable :: dim_names(:) 
        logical :: with_time 

        ! Local shortcut
        with_time = vs%par%with_time 

        ! First make sure all data objects are deallocated 
        call varslice_end(vs)

        ! Get information from netcdf file 
        call nc_dims(vs%par%filename,vs%par%name,dim_names,vs%par%dim)
        vs%par%ndim = size(vs%par%dim,1)

        if (with_time) then

            ! Initialize time vector from user parameters 
            call axis_init(vs%time,x0=vs%par%time_par(1), &
                                   x1=vs%par%time_par(2), &
                                   dx=vs%par%time_par(3))

            ! Check to make sure time vector matches netcdf file length 
            if (size(vs%time,1) .ne. vs%par%dim(vs%par%ndim)) then 
                write(*,*) "varslice_init_data:: Error: generated time coordinate &
                &does not match the length of the time dimension in the netcdf file."
                write(*,*) "time_par:    ", vs%par%time_par 
                write(*,*) "size(time):  ", size(vs%time,1)
                write(*,*) "nt (netcdf): ", vs%par%dim(vs%par%ndim)
                write(*,*) "filename:    ", trim(vs%par%filename)
                stop 
            end if 

        end if 

        ! Allocate coordinate and data variables,
        ! load coordinates too 
        select case(vs%par%ndim)

            case(1)
                if (with_time) then 
                    allocate(vs%var(1,1,1,1))
                else 
                    allocate(vs%var(vs%par%dim(1),1,1,1))
                end if 
            case(2)
                if (with_time) then 
                    allocate(vs%x(vs%par%dim(1)))
                    allocate(vs%var(vs%par%dim(1),1,1,1))

                    if (nc_exists_var(vs%par%filename,dim_names(1))) then 
                        call nc_read(vs%par%filename,dim_names(1),vs%x)
                    else
                        call axis_init(vs%x,nx=vs%par%dim(1))
                    end if
                else 
                    allocate(vs%x(vs%par%dim(1)))
                    allocate(vs%y(vs%par%dim(2)))
                    allocate(vs%var(vs%par%dim(1),vs%par%dim(2),1,1))

                    if (nc_exists_var(vs%par%filename,dim_names(1))) then 
                        call nc_read(vs%par%filename,dim_names(1),vs%x)
                    else
                        call axis_init(vs%x,nx=vs%par%dim(1))
                    end if
                    if (nc_exists_var(vs%par%filename,dim_names(2))) then 
                        call nc_read(vs%par%filename,dim_names(2),vs%y)
                    else
                        call axis_init(vs%y,nx=vs%par%dim(2))
                    end if
                    
                end if 

            case(3)
                if (with_time) then
                    allocate(vs%x(vs%par%dim(1)))
                    allocate(vs%y(vs%par%dim(2)))
                    allocate(vs%var(vs%par%dim(1),vs%par%dim(2),1,1))

                    if (nc_exists_var(vs%par%filename,dim_names(1))) then 
                        call nc_read(vs%par%filename,dim_names(1),vs%x)
                    else
                        call axis_init(vs%x,nx=vs%par%dim(1))
                    end if
                    if (nc_exists_var(vs%par%filename,dim_names(2))) then 
                        call nc_read(vs%par%filename,dim_names(2),vs%y)
                    else
                        call axis_init(vs%y,nx=vs%par%dim(2))
                    end if
                    
                else 
                    allocate(vs%x(vs%par%dim(1)))
                    allocate(vs%y(vs%par%dim(2)))
                    allocate(vs%lev(vs%par%dim(3)))
                    allocate(vs%var(vs%par%dim(1),vs%par%dim(2),vs%par%dim(3),1))

                    if (nc_exists_var(vs%par%filename,dim_names(1))) then 
                        call nc_read(vs%par%filename,dim_names(1),vs%x)
                    else
                        call axis_init(vs%x,nx=vs%par%dim(1))
                    end if
                    if (nc_exists_var(vs%par%filename,dim_names(2))) then 
                        call nc_read(vs%par%filename,dim_names(2),vs%y)
                    else
                        call axis_init(vs%y,nx=vs%par%dim(2))
                    end if
                    if (nc_exists_var(vs%par%filename,dim_names(3))) then 
                        call nc_read(vs%par%filename,dim_names(3),vs%lev)
                    else
                        call axis_init(vs%lev,nx=vs%par%dim(3))
                    end if
                    
                end if
                  
            case(4)
                if (with_time) then 
                    allocate(vs%x(vs%par%dim(1)))
                    allocate(vs%y(vs%par%dim(2)))
                    allocate(vs%lev(vs%par%dim(3)))
                    allocate(vs%var(vs%par%dim(1),vs%par%dim(2),vs%par%dim(3),1))

                    if (nc_exists_var(vs%par%filename,dim_names(1))) then 
                        call nc_read(vs%par%filename,dim_names(1),vs%x)
                    else
                        call axis_init(vs%x,nx=vs%par%dim(1))
                    end if
                    if (nc_exists_var(vs%par%filename,dim_names(2))) then 
                        call nc_read(vs%par%filename,dim_names(2),vs%y)
                    else
                        call axis_init(vs%y,nx=vs%par%dim(2))
                    end if
                    if (nc_exists_var(vs%par%filename,dim_names(3))) then 
                        call nc_read(vs%par%filename,dim_names(3),vs%lev)
                    else
                        call axis_init(vs%lev,nx=vs%par%dim(3))
                    end if
                    
                else 
                    write(*,*) "varslice_init_data:: 4D array without time dimension is not yet supported."
                    stop 
                end if

                    
            case DEFAULT 
                write(*,*) "varslice_init_data:: ndim > 4 not allowed."
                write(*,*) "ndim = ", vs%par%ndim 
                stop 

        end select

        return 

    end subroutine varslice_init_data

    subroutine varslice_end(vs)
        ! Deallocate all variables

        implicit none 

        type(varslice_class), intent(INOUT) :: vs 

        if (allocated(vs%par%dim))  deallocate(vs%par%dim)
        if (allocated(vs%x))        deallocate(vs%x)
        if (allocated(vs%y))        deallocate(vs%y)
        if (allocated(vs%lev))      deallocate(vs%lev)
        if (allocated(vs%time))     deallocate(vs%time)
        if (allocated(vs%var))      deallocate(vs%var)
        
        return 

    end subroutine varslice_end

    ! jablasco added gcm scen
    subroutine varslice_par_load(par,filename,group,domain,grid_name,gcm,init,verbose)

        type(varslice_param_class), intent(OUT) :: par 
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: group
        character(len=*), intent(IN), optional :: domain
        character(len=*), intent(IN), optional :: grid_name  
        character(len=*), intent(IN), optional :: gcm
        logical, optional :: init 
        logical, optional :: verbose 

        ! Local variables
        logical  :: init_pars 
        real(wp) :: time_par(3) 
        logical  :: print_summary 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 

        print_summary = .TRUE. 
        if (present(verbose)) print_summary = verbose 

        call nml_read(filename,group,"filename",       par%filename,     init=init_pars)
        call nml_read(filename,group,"name",           par%name,         init=init_pars)
        call nml_read(filename,group,"units_in",       par%units_in,     init=init_pars)
        call nml_read(filename,group,"units_out",      par%units_out,    init=init_pars)
        call nml_read(filename,group,"unit_scale",     par%unit_scale,   init=init_pars)   
        call nml_read(filename,group,"unit_offset",    par%unit_offset,  init=init_pars)   
        call nml_read(filename,group,"with_time",      par%with_time,    init=init_pars)   
        call nml_read(filename,group,"time_par",       par%time_par,     init=init_pars)   
        
        ! Parse filename as needed
        if (present(domain) .and. present(grid_name)) then
            call parse_path(par%filename,domain,grid_name)
        end if 

        ! Parse filename as needed
        if (present(gcm)) then
            call parse_path_ismip6(par%filename,gcm)
        end if

        ! Make sure time parameters are consistent time_par=[x0,x1,dx]
        if (par%time_par(3) .eq. 0) par%time_par(2) = par%time_par(1) 

        ! Summary 
        if (print_summary) then  
            write(*,*) "Loading: ", trim(filename), ":: ", trim(group)
            write(*,*) "filename    = ", trim(par%filename)
            write(*,*) "name        = ", trim(par%name)
            write(*,*) "units_in    = ", trim(par%units_in)
            write(*,*) "units_out   = ", trim(par%units_out)
            write(*,*) "unit_scale  = ", par%unit_scale
            write(*,*) "unit_offset = ", par%unit_offset
            write(*,*) "with_time   = ", par%with_time
            if (par%with_time) then
                write(*,*) "time_par    = ", par%time_par
            end if
        end if 

        return

    end subroutine varslice_par_load

    subroutine parse_path(path,domain,grid_name)

        implicit none 

        character(len=*), intent(INOUT) :: path 
        character(len=*), intent(IN)    :: domain, grid_name 

        call nml_replace(path,"{domain}",   trim(domain))
        call nml_replace(path,"{grid_name}",trim(grid_name))
        
        return 

    end subroutine parse_path
   
    ! jablasco
    subroutine parse_path_ismip6(path,gcm)

        implicit none

        character(len=*), intent(INOUT) :: path
        character(len=*), intent(IN)    :: gcm

        call nml_replace(path,"{gcm}",   trim(gcm))

        return

    end subroutine parse_path_ismip6
 
    subroutine axis_init(x,x0,x1,dx,nx)

        implicit none 

        real(wp), allocatable, intent(OUT) :: x(:)
        real(wp), optional :: x0
        real(wp), optional :: x1
        real(wp), optional :: dx
        integer,  optional :: nx 

        ! Local variables 
        integer :: i  
        real(wp) :: x0_now
        real(wp) :: x1_now
        real(wp) :: dx_now
        integer  :: nx_now 

        dx_now = 1.0_wp 
        if (present(dx)) dx_now = dx 

        x0_now = 0.0_wp 
        if (present(x0)) x0_now = x0 
        
        ! Note: if x1 is present, nx is ignored 
        if (present(x1)) then 
            x1_now = x1 
        else if (present(nx)) then 
            x1_now = (nx-1)*dx_now 
        else 
            write(*,*) "axis_init:: Error: either x1 or nx must be present."
            stop 
        end if 

        if (allocated(x)) deallocate(x)

        nx_now = (x1_now-x0_now)/dx_now + 1
        allocate(x(nx_now))

        do i = 1, nx_now 
            x(i) = x0_now + dx_now*real(i-1,wp)
        end do 

        return

    end subroutine axis_init

end module varslice
