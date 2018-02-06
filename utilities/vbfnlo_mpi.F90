! this  module provides helper variables/functions for mpi

! The overall MPI concept is
! * there is a master process (vbfnlo_mpi_masterid, assumed 0), which
!   * is responsible for input/output, to stdout and files
!   * always has the combined final data
! * there are worker processes (1:n)
!   * have output disabled
!   * run the full initialization, i.e. parse input, set couplings, ...
!   * run a fraction of the phasespace points
! * in general ALLREDUCE is used, so that also the workers have the combined information
!   * this allows them to update their grids themselves
!   * more code can be reused that assumes there is only 1 process running

#ifdef WITH_MPI
module vbfnlo_mpi
    use mpi
    use monaco, only: monaco_end_of_iteration, jump_next_numbers
    integer, public :: vbfnlo_mpi_myid
    integer, public :: vbfnlo_mpi_nprocs
    integer, private :: ierr
    integer, private :: comm
    integer, public, parameter :: vbfnlo_mpi_masterid = 0
    integer*8, public :: ncall_thisproc
    integer*8, private :: ncall_start

    logical, parameter, private :: ldebug = .false.

contains
    subroutine InitMPI()
        use globalvars, only: lglobalprint
        
        call MPI_INIT(ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, vbfnlo_mpi_myid, ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, vbfnlo_mpi_nprocs, ierr)
        print*, 'Started MPI job ', vbfnlo_mpi_myid+1, '/', vbfnlo_mpi_nprocs

        if (vbfnlo_mpi_myid /= vbfnlo_mpi_masterid) then
            lglobalprint = .false.
        endif
    end subroutine

    subroutine FinalMPI()
        call MPI_FINALIZE(ierr)
    end subroutine

    subroutine combine_grids()
        use monaco, only: grid, g, ngrid, ndimen

        type(grid) :: summedgrid
        allocate(summedgrid%d(ngrid, ndimen))
        allocate(summedgrid%di(ngrid, ndimen))

        call MPI_ALLREDUCE(g%fb, summedgrid%fb, 1, MPI_DOUBLE, MPI_SUM, comm, ierr)
        call MPI_ALLREDUCE(g%f2b, summedgrid%f2b, 1, MPI_DOUBLE, MPI_SUM, comm, ierr)
        call MPI_ALLREDUCE(g%d, summedgrid%d, ngrid*ndimen, MPI_DOUBLE, MPI_SUM, comm, ierr)
        call MPI_ALLREDUCE(g%di, summedgrid%di, ngrid*ndimen, MPI_DOUBLE, MPI_SUM, comm, ierr)
        call MPI_ALLREDUCE(g%ngtz, summedgrid%ngtz, 1, MPI_DOUBLE, MPI_SUM, comm, ierr)
        call MPI_ALLREDUCE(g%nltz, summedgrid%nltz, 1, MPI_DOUBLE, MPI_SUM, comm, ierr)
        call MPI_ALLREDUCE(g%nez, summedgrid%nez, 1, MPI_DOUBLE, MPI_SUM, comm, ierr)
        call MPI_ALLREDUCE(g%ncalls, summedgrid%ncalls, 1, MPI_DOUBLE, MPI_SUM, comm, ierr)
        call MPI_ALLREDUCE(g%famax, summedgrid%famax, 1, MPI_DOUBLE, MPI_MAX, comm, ierr)
        call MPI_ALLREDUCE(g%frmax, summedgrid%frmax, 1, MPI_DOUBLE, MPI_MAX, comm, ierr)
        ! print*, 'MPI_REDUCE', vbfnlo_mpi_myid, ' of ', vbfnlo_mpi_nprocs
        g = summedgrid
    end subroutine

    subroutine combine_histograms()
        use hist_stor

        include 'hist.inc'

        integer hist_in
        integer hlen
        integer hlenx, hleny
        type(hist_store), dimension(histmax) :: total_hist
        type(hist_store2d), dimension(histmax) :: total_hist2d

        do hist_in = 1, 2*NLOoffset
            if (hstored(hist_in)) then
                hlen = size(hist(hist_in)%hsum)
                allocate(total_hist(hist_in)%hsum(hlen))
                allocate(total_hist(hist_in)%hsquare(hlen))
                allocate(total_hist(hist_in)%hcount(hlen))
                call MPI_ALLREDUCE(hist(hist_in)%hsum, total_hist(hist_in)%hsum, hlen, MPI_DOUBLE, MPI_SUM, comm, ierr)
                call MPI_ALLREDUCE(hist(hist_in)%hsquare, total_hist(hist_in)%hsquare, hlen, MPI_DOUBLE, MPI_SUM, comm, ierr)
                call MPI_ALLREDUCE(hist(hist_in)%hcount, total_hist(hist_in)%hcount, hlen, MPI_DOUBLE, MPI_SUM, comm, ierr)
                hist(hist_in) = total_hist(hist_in)
            endif
        enddo

        do hist_in = 1, 2*NLOoffset2d
            if (h2dstored(hist_in)) then
                hlenx = size(hist2d(hist_in)%hsum, 1)
                hleny = size(hist2d(hist_in)%hsum, 2)
                allocate(total_hist2d(hist_in)%hsum(hlenx, hleny))
                allocate(total_hist2d(hist_in)%hsquare(hlenx, hleny))
                allocate(total_hist2d(hist_in)%hcount(hlenx, hleny))
                call MPI_ALLREDUCE(hist2d(hist_in)%hsum, total_hist2d(hist_in)%hsum, hlenx*hleny, MPI_DOUBLE, MPI_SUM, comm, ierr)
                call MPI_ALLREDUCE(hist2d(hist_in)%hsquare, total_hist2d(hist_in)%hsquare, hlenx*hleny, MPI_DOUBLE, MPI_SUM, comm, ierr)
                call MPI_ALLREDUCE(hist2d(hist_in)%hcount, total_hist2d(hist_in)%hcount, hlenx*hleny, MPI_DOUBLE, MPI_SUM, comm, ierr)
                hist2d(hist_in) = total_hist2d(hist_in)
            endif
        enddo

    end subroutine

    subroutine mpi_prepare_loop(iteration, ncall)
        use monaco, only: rtype, RTYPE_XORSHIFT
        integer*8, intent(in) :: ncall
        integer, intent(in) :: iteration
        integer*8 ncall_worker
        integer i

        ncall_worker = ncall/vbfnlo_mpi_nprocs
        if (vbfnlo_mpi_myid == vbfnlo_mpi_masterid) then
            ! make sure the total number of calls adds up
            ncall_thisproc = ncall - (vbfnlo_mpi_nprocs-1) * ncall_worker
        else
            ncall_thisproc = ncall_worker
        endif

        if (rtype /= RTYPE_XORSHIFT) then
            ncall_start = (vbfnlo_mpi_myid - 1) * ncall_worker
            if (vbfnlo_mpi_myid == vbfnlo_mpi_masterid) then
                ! make sure the total number of calls adds up
                ncall_start =  ncall - ncall_thisproc
            endif
            call jump_next_numbers(ncall_start)
        endif

        if (rtype == RTYPE_XORSHIFT .and. iteration == 1) then ! only xorshift has a quick jump method
            do i=1,vbfnlo_mpi_myid ! maximum here: 2**64
                call xorshift_jump()
            enddo
        endif
        if (ldebug) then
            print*, 'MPI ID ', vbfnlo_mpi_myid, 'nstart', ncall_start, '+', ncall_thisproc, '/', ncall
        endif
        !                  print*, 'callstart', vbfnlo_mpi_myid, ncall_start, ncall_start+ncall_thisproc
    endsubroutine

    subroutine mpi_after_loop(ncall)
        use monaco_rng_mz, only: monran_print_state
        use monaco, only: rtype, RTYPE_XORSHIFT
        integer*8, intent(in) :: ncall
        call combine_grids()
        call monaco_end_of_iteration()
        ! print*, 'before jump'
        ! TODO: this jump is only to ensure we use exactly the same random numbers as for a single run; 
        ! it can probably be removed by using different seeds
        if (rtype /= RTYPE_XORSHIFT) then
            call jump_next_numbers(ncall - ncall_thisproc - ncall_start)
        endif
        if (ldebug) call monran_print_state()
        ! print*, 'after jump'
        ! call monran_print_state()
    endsubroutine

end module
#endif
