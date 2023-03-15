module timeout

    use nml 
    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    implicit none

    integer,  parameter :: sp  = kind(1.0)
    integer,  parameter :: dp  = kind(1.0d0)
    integer,  parameter :: wp  = sp 

    real(wp), parameter :: MV  = -9999.0_wp 

    logical,  parameter :: verbose = .TRUE. 

    type timeout_class
        character(len=56)     :: method
        real(wp), allocatable :: times(:)
    end type

    private
    public :: timeout_class
    public :: timeout_init 

contains

    subroutine timeout_init(tm,filename,group,time_init,time_end)

        implicit none 

        type(timeout_class), intent(INOUT) :: tm
        character(len=*),    intent(IN)    :: filename 
        character(len=*),    intent(IN)    :: group 
        real(wp),            intent(IN)    :: time_init 
        real(wp),            intent(IN)    :: time_end
        
        ! Local variables
        integer  :: k, n, k0, k1 
        real(wp) :: dt_const
        character(len=512) :: timeout_file 

        integer, parameter :: nmax = 100000
        real(wp) :: times(nmax)

        ! Load parameters
        call nml_read(filename,group,"method",tm%method)
        
        times = MV 

        select case(trim(tm%method))

            case("const")
                
                call nml_read(filename,group,"dt",dt_const)

                if (dt_const .gt. 0.0) then
                    
                    n = ceiling( (time_end-time_init)/dt_const ) + 1
                    
                    if (n .gt. nmax) then 
                        write(error_unit,*) "timeout_init:: Error: too many output timesteps desired. &
                        &Maximum value limited to nmax = ", nmax 
                        write(error_unit,*) "time_init = ", time_init 
                        write(error_unit,*) "time_end  = ", time_end 
                        write(error_unit,*) "dt        = ", dt_const 
                        write(error_unit,*) "n         = ", n 
                        stop
                    end if

                    do k = 1, n 
                        times(k) = time_init + (k-1)*dt_const 
                    end do

                    k0 = 1 
                    k1 = k0+n-1

                end if

            case("file","times")

                if (trim(tm%method) .eq. "file") then 
                    ! Load time information from an input file 

                    call nml_read(filename,group,"file",timeout_file)

                    ! Get a vector of times 
                    call timeout_parse_file(times,timeout_file)

                else
                    ! Load time information from a parameter
                
                    call nml_read(filename,group,"times",times,init=.FALSE.)

                end if 

                k0 = -1 
                k1 = -1 

                do k = 1, nmax 
                    if (k0 .lt. 0 .and. times(k) .ge. time_init .and. times(k) .ne. MV) then 
                        k0 = k
                    end if 
                    if (k0 .gt. 0 .and. times(k) .le. time_end  .and. times(k) .ne. MV) then
                        k1 = k 
                    end if
                end do

                if (k0 .lt. 0 .or. k1 .lt. 0) then 
                    ! No proper range of indices was found, 
                    ! prescribe first and list times as output times.
                    k0 = 1
                    k1 = 2
                    times = MV 
                    times(1) = time_init
                    times(2) = time_end 
                end if

                n = k1-k0+1

            case DEFAULT

                write(error_unit,*) "timeout_init:: Error: timeout method not recognized."
                write(error_unit,*) "timeout.method = ", trim(tm%method)
                stop 
            
        end select 

        if (allocated(tm%times)) deallocate(tm%times)
        allocate(tm%times(n))
        tm%times(1:n) = times(k0:k1)


        if (verbose) then 
            ! Print summary
            write(*,*) "timeout: ", trim(tm%method)
            write(*,*) "  time_init = ", time_init 
            write(*,*) "  time_end  = ", time_end 
            write(*,*) "  n         = ", n 

            do k = 1, n, 5
                k0 = k 
                k1 = min(k+5-1,n)
                write(*,"(5g12.3)") tm%times(k0:k1)
            end do 

        end if 

        return

    end subroutine timeout_init

    subroutine timeout_parse_file(times,filename)

        implicit none 

        real(wp),         intent(INOUT) :: times(:) 
        character(len=*), intent(IN)    :: filename 

        ! Local variables 
        integer :: q, k, nmax, io, stat  
        character(len=256) :: line_str 

        integer, parameter :: nmax_file = 10000

        open(newunit=io,file=filename,status="old",action="read")

        k = 0 

        do q = 1, nmax_file 

            ! Read line as a string first 
            read(io,"(a10000)",iostat=stat) line_str 
            line_str = trim(adjustl(line_str))

            if (stat .ne. 0) exit
            if (trim(line_str) .eq. "" .or. line_str(1:1) .eq. "#") then 
                ! Skip this line, do nothing, it is either a comment or an empty line
            else
                ! This line should contain information we want

                if (index(line_str,":") .gt. 0) then 
                    ! Parse this line as a set of times time1:dtime:time2


                else
                    ! Assume only one value is available

                    k = k+1
                    times(k) = string_to_wp(line_str)
                    !write(*,*) k, " ", times(k)

                end if

            end if

        end do 


        close(io) 

        return 

    end subroutine timeout_parse_file

    subroutine parse_time_vector(times,timestr)

        implicit none

        real(wp),         intent(INOUT) :: times(:)
        character(len=*), intent(IN)    :: timestr 

        return

    end subroutine parse_time_vector

    function timeout_check(tm,time) result(out_now)

        implicit none

        type(timeout_class), intent(IN) :: tm
        real(wp), intent(IN) :: time 
        real(wp) :: out_now 

        return

    end function timeout_check


    ! === Helper functions (borrowed from nml.f90) ===
    
    function string_to_wp(string) result(value)

        implicit none 

        character(len=*), intent(IN) :: string 
        real(wp) :: value 

        character(len=256) :: tmpstr 
        integer  :: stat, n
        real(wp) :: x 

        tmpstr = trim(adjustl(string))
        n      = len_trim(tmpstr)

        read(tmpstr(1:n),*,IOSTAT=stat) x

        value = 0
        if (stat .eq. 0) then 
            value = x 
        else
            n = len_trim(tmpstr)-1
            READ(tmpstr(1:n),*,IOSTAT=stat) x
            if (stat .ne. 0) then 
                write(*,*) "nml:: ","Error converting string to number!"
                write(*,*) "|",trim(tmpstr),"|",n,stat,x
            else
                value = x 
            end if 
        end if 

        return 

    end function string_to_wp

    subroutine string_to_vector(string,value)

        implicit none 

        character(len=*), intent(IN) :: string 
        character(len=*) :: value(:)
        character(len=256) :: tmpvec(size(value))
        character(len=256) :: tmpstr, fmt 
        integer :: stat, n, q, q1, q2, j 

        tmpstr = trim(adjustl(string))
        n      = len_trim(tmpstr)+2

        tmpvec(:) = "" 

        q1 = 1 
        do q = 1, size(tmpvec)
            q2 = index(tmpstr(q1:n)," ") + q1
            if (q2 .gt. q1 .and. q2 .le. n) then 
                tmpvec(q) = tmpstr(q1:q2-1)
                q1 = q2

                ! Make sure gaps of more than one space are properly handled
                do j = 1, 1000
                    if (tmpstr(q1:q1) == " ") q1 = q1+1
                    if (q1 .ge. n) exit 
                end do 

!                 ! Eliminate quotes
!                 q2 = len_trim(tmpvec(q))
!                 if (tmpvec(q)(1:1) == '"') tmpvec(q) = trim(adjustl(tmpvec(q)(2:q2)))
!                 q2 = len_trim(tmpvec(q))
!                 if (tmpvec(q)(q2:q2) == '"') tmpvec(q) = trim(tmpvec(q)(1:q2-1))
                ! Remove quotes around string if they exist 
                call remove_quotes_comma(tmpvec(q))
            
            end if 
        end do 
        
        value = tmpvec 

        return 

    end subroutine string_to_vector

    subroutine remove_quotes_comma(string)

        implicit none 
        character(len=*), intent(INOUT) :: string 
        integer :: i, n 

!         ! Eliminate quotes
!         n = len_trim(string)
!         if (n == 1 .and. trim(string) == '"') then 
!             string = ""
!         else if (n > 0) then 
!             if (string(1:1) == '"') string = trim(adjustl(string(2:n)))
!             n = len_trim(string)
!             if (n > 1  .and. string(n:n) == '"') string = trim(string(1:n-1))
!             if (n == 1 .and. string(n:n) == '"') string = ""
            
!         end if 

        ! Eliminate quotes
        n = len_trim(string)
        do i = 1,n 
            if (string(i:i) == '"' .or. string(i:i) == "'") string(i:i) = " "
        end do 
        string = trim(adjustl(string))

        ! Remove final comma too
        n = len_trim(string)
        if (n > 0) then 
            if (string(n:n) == ",") string(n:n) = " "
            string = trim(adjustl(string))
        end if 
        
        return 

    end subroutine remove_quotes_comma

end module timeout