!   An extremely modified version of VEGAS - fundamental
!   structural differences, along with selection of
!   pseudo-random or quasi-random vector generation,
!   checking for special cases, and modified output
!   statements, among other things, distinguish the two.

!   Adam Duff, <duff@phenom.physics.wisc.edu>
!   Initial version:  1991 June 25
!   Last modified:  2000 Dec 5 by D.Zeppenfeld

!       new version has all real*8 defined as double precision
!       and allows to identify names for input and output grids as
!       variables in call of monaco_read and monaco_write

!       all integer (except for ndim in monaco_init... calls) declared as
!       integer*8 to allow more than 2^30 events. Only change in monaco use
!       is declararion of npt (=ncall) as integer*8

!   Subroutine performs n-dimensional monte carlo integ"n
!      - by g.p. lepage    sept 1976/(rev)aug 1979
!      - algorithm described in j comp phys 27,192(1978)

module monaco
    use globalvars, only: lglobalprint
   implicit none

! declare global parameters
   integer, parameter, public :: ngrid=48, ndi=24, nxi=ngrid*ndi

! declare external common variables
   integer :: nprn = 0            !default print cumulative information
   integer :: ndev = 6            !default output device
   double precision, dimension(ndi) :: xl = 0.0d0        !default integration lower bounds
   double precision, dimension(ndi) :: xu = 1.0d0        !default integration upper bounds
   double precision :: acc = -1.0d0            !default termination accuracy
   integer :: it = 0            !interations completed
   integer :: ndo = 1            !number of subdivisions per axis
   double precision :: si = 0.0d0            !sum S / sigma^2
   double precision :: swgt = 0.0d0            !sum 1 / sigma^2
   double precision :: schi = 0.0d0            !sum S^2 / sigma^2
   double precision, allocatable, dimension(:, :) :: xi  ! ngrid, ndimen       !location of i-th division on j-th axis
   double precision :: alph = 0.875d0         !default convergence parameter
   integer :: ndmx = 48            !default number of grid divisions ! TODO combine with ngrid
   integer :: mds = 0            !use stratified and/or importance sampling
   
   integer, public :: rtype  ! seed RTYPE_ defines
   integer, public, parameter :: RTYPE_MZ=0
   integer, public, parameter :: RTYPE_SOBOL=1
   integer, public, parameter :: RTYPE_BUILTIN=2
   integer, public, parameter :: RTYPE_XORSHIFT=3

   double precision calls
   double precision ti ! integral for one iteration
   double precision tsi ! sigma for one iteration
   double precision reffic
   double precision aeffic
   double precision fltz
   double precision fez
   double precision fgtz

! declare internal common variables
   double precision xjac
   integer nd ! =ndmx=48, number of grid cuts
   integer*8 ng ! =1
   integer*8 npg ! =ncall
   integer*8 ndm ! =nd-1
   double precision dv2g ! 1/(npg-1) for ng=1

   type :: grid
      double precision, allocatable, dimension(:, :) :: d !(ngrid,ndi) ! sum of f**2
      double precision, allocatable, dimension(:, :) :: di !(ngrid,ndi) ! sum of f
      double precision fb ! sum f
      double precision f2b ! sum f**2
      integer*8 nltz, nez, ngtz ! number of points less, equal and greater than 0
      integer*8 ncalls
      double precision frmax
      double precision famax
   end type
   type(grid) :: g ! all values that get filled by monaco_put

   integer*8, allocatable, dimension(:) :: kg !(ndi)
   integer, public :: ndimen = 0
   double precision, allocatable, dimension(:) :: dx !(ndi)
   double precision :: dxg, xnd
   double precision avgi !total integral over all iterations
   double precision sd ! total sigma over all iterations
   double precision chi2a ! chi^2/iteration

   integer, allocatable, dimension(:) :: ia ! ndimen

   save
   private

   public :: reweight_grid, monaco_init, monaco_init2, monaco_init3
   public :: monaco_get, monaco_put, monaco_read, monaco_write
   public :: monaco_result, monaco_end_of_iteration
   public :: grid, g
   public :: jump_next_numbers
   
contains

   subroutine monaco_init1(ndim, npt)
   integer*4, intent(in) :: ndim !in:  number of dimensions in hypercube 
   integer*8, intent(in) :: npt !in:  number of evaluations per iteration 
      call monaco_init(ndim, npt, oldinitid=1)
   end subroutine

   subroutine monaco_init2(ndim, npt)
   integer*4, intent(in) :: ndim !in:  number of dimensions in hypercube 
   integer*8, intent(in) :: npt !in:  number of evaluations per iteration 
      call monaco_init(ndim, npt, oldinitid=2)
   end subroutine
   
   subroutine monaco_init3(ndim, npt)
   integer*4, intent(in) :: ndim !in:  number of dimensions in hypercube 
   integer*8, intent(in) :: npt !in:  number of evaluations per iteration 
      call monaco_init(ndim, npt, oldinitid=3)
   end subroutine

   subroutine monaco_init(ndim, npt, oldinitid)
      use monaco_rng_mz, only: imonrn
      use monaco_rng_sob, only: imonso
      use globalvars, only: seed
      
   implicit none

   integer*4, intent(in) :: ndim !in:  number of dimensions in hypercube 
   integer*8, intent(in) :: npt !in:  number of evaluations per iteration 

   integer, intent(in), optional :: oldinitid

   integer initid
   integer k
   integer*4 seed1, seed2

! declare local variables
   integer*8 i, j
   double precision xin(ngrid), rc, xn, dr, xo

   if (present(oldinitid)) then
      initid = oldinitid
   else
      initid = 0
   endif

   if (ndim /= ndimen ) then
      ndimen = ndim
      if (allocated(xi)) deallocate(xi)
      allocate(xi(ngrid, ndimen))
      xi = 1.0

      if (allocated(dx)) deallocate(dx)
      allocate(dx(ndimen))

      if (allocated(ia)) deallocate(ia)
      allocate(ia(ndimen))

   endif

   if (.not. lglobalprint) nprn = -1

   if (initid < 1) then
   ! set the two seeds
      if (mod(seed,2).eq.0) then
         seed1 = seed / 2
         seed2 = seed / 2
      else
         seed1 = (seed + 1) / 2
         seed2 = (seed - 1) / 2
      endif

   ! test to see that ndimen is less than available from either the
   ! dimensionality of the grid, or the Sobol generator (40).

      if ( ndimen .gt. min( ndi, 40 ) ) then
         write(ndev,*)
         write(ndev,*) "MONACO called with ndim > ndi"
         write(ndev,*) "ndim =", ndimen, ",  ndi =", ndi
         stop
      end if


      if (rtype == RTYPE_MZ) then
         call imonrn( ndimen, seed1, seed2 )
      elseif (rtype == RTYPE_SOBOL) then
         call imonso( ndimen )
      elseif (rtype == RTYPE_BUILTIN) then
          ! TODO: put a seed here; handling depends on gfortran version, so dangerous
          call random_seed()
      elseif (rtype == RTYPE_XORSHIFT) then
         call xorshift_seed(ndimen, seed)
      else
          write(ndev,*)
          write(ndev,*) "MONACO:  invalid random sequence generator choice"
          write(ndev,*) "rtype =", rtype
          stop
      endif
      ndo = 1
      xi(1,1:ndimen) = 1.0d0
      avgi = 0.0d0
      sd = 0.0d0
      chi2a = 0.0d0
   endif

   if (initid < 2) then
      ! initializes cumulative variables, but not grid
      it = 0
      si = 0.0d0
      swgt = 0d0
      schi = 0d0
   endif


   if (initid < 3) then
   ! no initialisation

      nd = ndmx
      ng = 1
      ! stratified sampling:
      if ( mds .ne. 0 ) then
         ng = int( (dble( npt ) / 2.0d0)**(1.0d0 / dble( ndimen )) )
         mds = 1
         if ( (2 * ng - ndmx) .ge. 0 ) then
            mds = -1
            npg = ng / ndmx + 1
            nd = ng / npg
            ng = npg * nd
         end if
      end if
      k = ng**ndimen
      npg = npt / k
      if ( npg .lt. 2 ) then
         npg = 2
      end if
      calls = dble( npg * k )
      dxg = 1.0d0 / dble( ng )
      dv2g = (calls * dxg**ndimen)**2 / (dble( npg )**2 * dble( npg - 1 ))
      xnd = dble( nd )
      ndm = nd - 1
      dxg = dxg * xnd
      xjac = 1.0d0 / calls
      do j = 1, ndimen
         dx(j) = xu(j) - xl(j)
         xjac = xjac * dx(j)
      end do

   ! rebin, preserving bin density

      if ( nd .ne. ndo ) then
         rc = dble( ndo ) / xnd
         do j = 1, ndimen
            k = 0
            xn = 0.0d0
            dr = xn
            do i = k+1,ndm
               do while (rc > dr)
                  k = k + 1
                  dr = dr + 1.0d0
                  xo = xn
                  xn = xi(k,j)
               enddo
               dr = dr - rc
               xin(i) = xn - (xn - xo) * dr
            enddo
            do i = 1, ndm
               xi(i,j) = xin(i)
            end do
            xi(nd,j) = 1.0d0
         end do
         ndo = nd
      end if

   200  format( / " MONACO input parameters:" / &
              " ndim =", i3, ",  ncall =", i8,a, "   rtype =", i2 / &
         " nprn =", i3, ",  acc = ", 1pe7.1e1, ",  alph = ", 1pe8.2e1, &
              ",  nd =", i3, ",  mds =", i3   )   !    /
   !     &        " integration bounds: (lower, upper)" /
   !     &        ( " dimension", i3, ":   ( ", 1pe13.6, ", ", 1pe13.6,
   !     &        " )" ) )
      if ( nprn .ge. 0 ) then
         if (calls.gt.1d6) then
            write(ndev,200) ndimen, int( calls/1024**2 ),"M,", rtype, &
                         nprn, acc, alph, nd, mds    !,
   !     &                   (j, xl(j), xu(j), j=1, ndim)
         elseif (calls.gt.1d4) then
            write(ndev,200) ndimen, int( calls )/1024,"k,", rtype, &
                         nprn, acc, alph, nd, mds    !,
   !     &                   (j, xl(j), xu(j), j=1, ndim)
         else
            write(ndev,200) ndimen, int( calls ),",", rtype, &
                         nprn, acc, alph, nd, mds    !,
   !     &                   (j, xl(j), xu(j), j=1, ndim)
         endif
      end if
   endif
! setup main integration loop
! note that number of points per iteration remains as defined above

   call reset_iteration_variables(dimchange=.true.)
   end subroutine

!-------------------------------------------------------------------------------

   subroutine monaco_get(r, wgt)
      use monaco_rng_mz, only: monran
      use monaco_rng_sob, only: monsob
   implicit none

! declare output variables
   double precision, intent(out) :: r(ndimen) ! generated point in hypercube
   double precision, intent(out) :: wgt ! generated weight of point
   
! declare local variables
   integer*8 j
   double precision rand(ndimen), xn, xo, rc

! generate point, with accompanying weight

   if ( rtype .eq. RTYPE_MZ ) then
      call monran( rand )
   else if ( rtype .eq. RTYPE_SOBOL ) then
      call monsob( rand )
   else if ( rtype .eq. RTYPE_BUILTIN ) then
      call random_number(rand)
   else if ( rtype .eq. RTYPE_XORSHIFT ) then
      call xorshift_getrand(rand)
   end if
   wgt = xjac
   do j = 1, ndimen
      if ( rand(j) .eq. 0.0d0 ) then
         rand(j) = 1.0d-15
      end if
      xn = (dble( kg(j) ) - rand(j)) * dxg + 1.0d0
      ia(j) = min( int( xn ), ngrid )
      if ( ia(j) .le. 1 ) then
         xo = xi(ia(j),j)
         rc = (xn - dble( ia(j) )) * xo
      else
         xo = xi(ia(j),j) - xi(ia(j)-1,j)
         rc = xi(ia(j)-1,j) + (xn - dble( ia(j) )) * xo
      end if
      r(j) = xl(j) + rc * dx(j)
      wgt = wgt * xo * xnd
   end do

   end subroutine

   subroutine jump_next_numbers(jumpsize)
      use monaco_rng_mz, only: get_next_rn
      use monaco_rng_sob, only: monsob
      integer*8, intent(in) :: jumpsize
      integer*8 :: i
      double precision :: throw_away
      double precision :: ran_dump(ndimen)
      if (rtype == 0) then
         do i = 1, jumpsize*ndimen
            throw_away = get_next_rn()
         enddo
      elseif (rtype == 1) then
         do i = 1, jumpsize
            call monsob(ran_dump)
         enddo
      endif
   end subroutine 



!-------------------------------------------------------------------------------

   subroutine monaco_put(r, wgt, value, rew_grid)
   implicit none

! declare input variables
   double precision, intent(in), dimension(ndimen) :: r      !in:  generated point in hypercube 
   double precision, intent(in) :: wgt   !in:  weight as generated by "monaco_get" 
   double precision, intent(in) :: value   !in:  value of integrand 
   double precision, intent(in) :: rew_grid  !in:  reweighting factor for grid generation in order to improve statistics in certain PS regions 
   ! logical, optional, intent(in) :: combine
   ! logical :: combine_grid

! declare local variables
   integer*8 j
   double precision fret, f, f2, rc

   ! if (present(combine)) then
   !    combine_grid = combine
   ! else
   !    combine_grid = .true.
   ! endif

! add point
   g%ncalls = g%ncalls + 1
   fret = value
   g%frmax = max( g%frmax, abs( fret ))
   if ( fret .ne. 0.0d0 ) then
      if ( fret .gt. 0.0d0 ) then
         g%ngtz = g%ngtz + 1
      else
         g%nltz = g%nltz + 1
      end if
   else
      g%nez = g%nez + 1
   end if

   f = wgt * fret
   g%famax = max( g%famax, abs( f ))
   f2 = f**2
   g%fb = g%fb + f
   g%f2b = g%f2b + f2

   ! modify f and f2 by rew_grid in order to improve statistics in certain PS regions
   f = f * rew_grid
   f2 = f**2


   do j = 1,ndimen
       g%di(ia(j),j) = g%di(ia(j),j) + f
       if ( mds .ge. 0 ) then 
           g%d(ia(j),j) = g%d(ia(j),j) + f2
       end if
   enddo

   end subroutine

   subroutine monaco_end_of_iteration()
   implicit none
   integer*8 i, j, k
   character(len=500) fmt_iteration, fmt_griddata

   double precision wgt, ti2

   g%f2b = sqrt( g%f2b * dble( npg ))
   g%f2b = (g%f2b - g%fb) * (g%f2b + g%fb)
   ti = ti + g%fb
   tsi = tsi + g%f2b

   ! stratified sampling
   if ( mds .lt. 0 ) then
      do j = 1, ndimen
         g%d(ia(j),j) = g%d(ia(j),j) + g%f2b
      end do
   end if
   k = ndimen
   do while (k > 0)
      kg(k) = mod( kg(k), ng ) + 1
      ! if ( kg(k) .ne. 1 ) then
      !    g%fb = 0.0d0
      !    g%f2b = 0.0d0
      !    k = 0
      !    return
      ! end if
      k = k - 1
   enddo

   ! TODO move last point of iteration to separate function
! compute final results for this iteration

   if ( g%frmax .ne. 0.0d0 ) then
      reffic = abs( ti ) / g%frmax
   else
      reffic = 0.0d0
   end if
   if ( g%famax .ne. 0.0d0 ) then
      aeffic = abs( ti ) / (g%famax * calls)
   else
      aeffic = 0.0d0
   end if
   fltz = dble( g%nltz ) / calls
   fez = dble( g%nez ) / calls
   fgtz = dble( g%ngtz ) / calls

   tsi = tsi * dv2g
   ti2 = ti**2
   if ( tsi .ne. 0.0d0 ) then
      wgt = 1.0d0 / tsi
   else
      wgt = 0.0d0
   end if
   si = si + ti * wgt
   swgt = swgt + wgt
   schi = schi + ti2 * wgt
   if ( swgt .ne. 0.0d0 ) then
      avgi = si / swgt
   else
      avgi = 0.0d0
   end if
   chi2a = (schi - si * avgi) / (dble( it ) - 0.999999d0)
   if ( swgt .gt. 0.0d0 ) then
      sd = sqrt( 1.0d0 / swgt )
   else
      sd = 0.0d0
   end if

   it = it + 1
   fmt_iteration = '( / " iteration", i3, ":" / &
           " integral = ", 1pe13.6, ",  sigma = ", 1pe9.3 / &
           " efficacy = ", 1pe9.3, ",  raw efficacy = ", 1pe9.3 / &
           " f_positive = ", 1pe8.3e1, ",  f_zero = ", 1pe8.3e1, &
           ",  f_negative = ", 1pe8.3e1 / &
           " accumulated statistics:" / &
           " integral = ", 1pe13.6, ",  sigma = ", 1pe9.3, &
           ",  chi^2/iteration = ", 1pe9.3 )'
   fmt_griddata = '( / " grid data for axis", i3, ":" / &
         & 5x, "x", 6x, "delta_i", 8x, "x", 6x, "delta_i", &
         & 8x, "x", 6x, "delta_i", 8x, "x", 6x, "delta_i", &
         & 8x, "x", 6x, "delta_i", 8x, "x", 6x, "delta_i" / &
         & (1x, 1pe9.3,1x, 1pe9.3, 3x, 1pe9.3,1x, 1pe9.3, 3x, &
         & 1pe9.3,1x, 1pe9.3, 3x, 1pe9.3,1x, 1pe9.3, 3x, &
         & 1pe9.3,1x, 1pe9.3, 3x, 1pe9.3,1x, 1pe9.3 ) )'
   if ( nprn .ge. 0 ) then
      if ( tsi .gt. 0.0d0 ) then
         tsi = sqrt( tsi )
      else
         tsi = 0.0d0
      end if
      write(ndev,fmt_iteration) it, &
         ti, tsi, &
         aeffic, reffic, &
         fgtz, fez, fltz, &
         avgi, sd, abs( chi2a )
      if ( nprn .gt. 0 ) then
         do j = 1, ndimen
            write(ndev,fmt_griddata) j, (xi(i,j), g%di(i,j), i = 1+nprn/2, nd)
         end do
      end if
   end if

   call refine_grid
   end subroutine

   subroutine refine_grid()
! declare local variables
   implicit none
   integer*8 i, j
   double precision dt(ndimen), xin(ngrid), ri(ngrid)
   double precision xo, xn, rc, dr
   integer k

   do j = 1, ndimen
      xo = g%d(1,j)
      xn = g%d(2,j)
      g%d(1,j) = (xo + xn) / 2.0d0
      dt(j) = g%d(1,j)
      do i = 2, ndm
         g%d(i,j) = xo + xn
         xo = xn
         xn = g%d(i+1,j)
         g%d(i,j) = (g%d(i,j) + xn) / 3.0d0
         dt(j) = dt(j) + g%d(i,j)
      end do
      g%d(nd,j) = (xn + xo) / 2.0d0
      dt(j) = dt(j) + g%d(nd,j)
   end do

   do j = 1, ndimen
      rc = 0.0d0
      do i = 1, nd
         ri(i) = 0.0d0
         if ((g%d(i,j).ge.0d0).and.(dt(j).ne.0d0)) then
            g%d(i,j) = max( 1.0d-30, g%d(i,j) )
            xo = dt(j) / g%d(i,j)
            ri(i) = ( (xo - 1.0d0) / (xo * log( xo )) )**alph
         else
            ri(i) = (  1d0 / log( 1d30 )  )**alph
         endif
         rc = rc + ri(i)
      end do
      rc = rc / xnd
      k = 0
      xn = 0.0d0
      dr = xn
      i = k
      do i=k+1,ndm
         do while (rc > dr)
            k = k + 1
            dr = dr + ri(k)
            xo = xn
            xn = xi(k,j)
         enddo
         dr = dr - rc
         xin(i) = xn - (xn - xo) * dr / ri(k)
      enddo
      do i = 1, ndm
         xi(i,j) = xin(i)
      end do
      xi(nd,j) = 1.0d0
   end do

   call reset_iteration_variables(dimchange=.false.)
   end subroutine

   subroutine reset_iteration_variables(dimchange)
      implicit none
      logical, intent(in) :: dimchange
      ti = 0.0d0
      tsi = 0.0d0
      g%nltz = 0
      g%nez = 0
      g%ngtz = 0
      g%frmax = 0.0d0
      g%famax = 0.0d0

      if (dimchange) then
         ! ensure the arrays have correct size, especially when switching to real emission
         if (allocated(kg) .and. size(kg) /= ndimen) deallocate(kg)
         if (allocated(g%d) .and. size(g%d, 2) /= ndimen) deallocate(g%d)
         if (allocated(g%di) .and. size(g%di, 2) /= ndimen) deallocate(g%di)
      endif

      if (.not.allocated(kg))   allocate(kg(ndimen))
      if (.not.allocated(g%d))   allocate(g%d(ngrid, ndimen))
      if (.not.allocated(g%di))   allocate(g%di(ngrid, ndimen))
      kg(1:ndimen) = 1
      g%d(1:nd,1:ndimen) = 0d0
      g%di(1:nd,1:ndimen) = 0d0

      g%fb = 0.0d0
      g%f2b = 0.0d0
      g%ncalls = 0
   end subroutine


   subroutine monaco_result(mean, sdev, chi2)
      implicit none

      ! declare output variables
      double precision, intent(out) :: mean   ! function integral 
      double precision, intent(out) :: sdev   ! standard deviation of integral 
      double precision, intent(out) :: chi2   ! chi^2 per degree of freedom 

      ! transfer values
      mean = avgi
      sdev = sd
      chi2 = abs( chi2a )
   end subroutine

!-------------------------------------------------------------------------------

   subroutine monaco_read(file_name)
   implicit none

   character*250, intent(in) :: file_name
   character(len=100) :: readformat

! declare local variables
   integer i, j, k
   integer ioerrread, ioerropen

! open file and read grid

   open(unit=15, file=file_name, iostat=ioerropen, status="old")
   if (ioerropen /= 0) then
      if (lglobalprint) write(6,'(/A)') " MONACO:  continuing with uniform grid"
      close( unit=15 )
      return
   endif

   readformat = '( 3( 1x, 1pd23.16 ) )'
   do j = 1, ndimen
      do i = 0, (ngrid-1)/3
         read(unit=15,fmt=readformat,iostat=ioerrread) ( xi(3*i+k,j), k=1, 3 )
      end do
   end do
   close( unit=15 )

   if (ioerrread /= 0 ) then
      if (lglobalprint) write(6,*) " MONACO:  read error on file unit 15"
      stop
   endif

   if (lglobalprint) write(6,'(/A, A)') " MONACO:  grid read from file ", file_name
   return

   end subroutine

!-------------------------------------------------------------------------------

   subroutine monaco_write(file_name)
   implicit none

   character*250, intent(in) :: file_name
   integer ioerroropen, ioerrorwrite

   character(len=100) :: writeformat

! declare local variables
   integer*8 i, j, k

! open file and write grid
   open( unit=16, file=file_name, iostat=ioerroropen, status="unknown" )
   writeformat = '( 3( 1x, 1pd23.16 ) )'
   if (ioerroropen /= 0) then
      write(6,'(/A)') " MONACO:  open error on file unit " // trim(file_name)
      return
   endif
   do j = 1, ndimen
      do i = 0, (ngrid-1)/3
         write(unit=16, fmt=writeformat, iostat=ioerrorwrite) ( xi(3*i+k,j), k=1, 3 )
      end do
   end do
   do j = ndimen+1, ndi
      do i = 0, (ngrid-1)/3
         write(unit=16, fmt=writeformat, iostat=ioerrorwrite) (/ 1d0, 1d0, 1d0 /)
      end do
   end do
   close( unit=16 )
   if(ioerrorwrite /= 0) then
      write(6,'(/A)') " MONACO:  write error on file unit " // trim(file_name)
      return
   endif

   if (lglobalprint) write(6,'(/A)') " MONACO:  grid written to file " // trim(file_name)
   end subroutine


! **************************************************************************
   double precision function reweight_grid (p, v)
! **************************************************************************

! Function which allows for a kinematics-dependent reweighting within the generation
! of grids for VEGAS, which can improve statistics in subdominant parts of distributions.
! Since VEGAS optimises the integration with respect to minimizing the error on the
! total cross section such regions are per default sampled poorly.
! In order to improve the error in this part of distributions it might be useful
! to define a function here which is large in this part of the PS and small in the
! other regions. One example for enhancing the statistics in high-mass tails
! of multi-boson production is provided below.
! While improving distributions, using a non-constant reweight_grid usually
! leads to larger errors on total cross sections. This feature should therefore
! be used only for generating distributions.

   implicit none

#include "VBFNLO/utilities/global.inc"

   real*8 p(0:3,max_p)
   real*8 v(0:3,max_v)
   real*8 s, s0, qvvv(0:3), mass2
   external mass2

   ! integer mu, i, expo

! -----------------------------------------

   ! default
   reweight_grid = 1d0

! -----------------------------------------

!      ! enhance high-mass tails of m_VV(V) for VV(V) and VBF V(V)jj production
!      ! useful e.g. for anomalous couplings analyses

!      s0 = 300d0**2     ! energy at which the reweighting sets in
!      expo = 4          ! exponent of the enhancement (s/s0) .
!                        ! for anomalous couplings runs a smaller exponent
!                        ! than for pure standard model runs might be useful

!      do mu = 0,3
!         qvvv(mu) = 0d0
!         do i = 1, n_v
!            qvvv(mu) = qvvv(mu) + v(mu,i)
!         enddo
!      enddo
!      s = mass2(qvvv)
!      if (s.gt.s0) then
!         reweight_grid = (s/s0)**expo
!      else
!         reweight_grid = 1d0
!      endif

! -----------------------------------------

   return
   end function

end module

