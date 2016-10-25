c
c	MONACO.F
c	An extremely modified version of VEGAS - fundamental
c	structural differences, along with selection of
c	pseudo-random or quasi-random vector generation,
c	checking for special cases, and modified output
c	statements, among other things, distinguish the two.
c
c	Adam Duff, <duff@phenom.physics.wisc.edu>
c	Initial version:  1991 June 25
c	Last modified:  2000 Dec 5 by D.Zeppenfeld
c
c       new version has all real*8 defined as double precision
c       and allows to identify names for input and output grids as
c       variables in call of monaco_read and monaco_write
c
c       all integer (except for ndim in monaco_init... calls) declared as 
c       integer*8 to allow more than 2^30 events. Only change in monaco use
c       is declararion of npt (=ncall) as integer*8
c
c	Subroutine performs n-dimensional monte carlo integ"n
c      - by g.p. lepage    sept 1976/(rev)aug 1979
c      - algorithm described in j comp phys 27,192(1978)
c
      subroutine monaco_init(
     &              ndim,	!in:  number of dimensions in hypercube
     &              npt,	!in:  number of evaluations per iteration
     &              seed        !in:  number for different seeds 
     &                      )
      implicit none
c
      integer*4 ndim
      integer*8 npt
      integer*4 seed1, seed2
      integer seed
c
c declare global parameters
c
      integer*4 ngrid, ndi, nxi
      parameter( ngrid=48, ndi=24, nxi=ngrid*ndi )
c
c declare external common variables
c
      integer*8 ncall, itmx, nprn, ndev
      double precision xl(ndi), xu(ndi), acc
      common /bveg1/ xl, xu, acc, ncall, itmx, nprn, ndev
c
      integer*8 it, ndo
      double precision si, swgt, schi, xi(ngrid,ndi)
      common /bveg2/ si, swgt, schi, xi, it, ndo
c
      integer*8 ndmx, mds
      double precision alph
      common /bveg3/ alph, ndmx, mds
c
      double precision calls, ti, tsi
      common /bveg4/ calls, ti, tsi
c
      integer*8 rtype
      common /bveg5/ rtype
c
      double precision reffic, aeffic, frmax, famax, fltz, fez, fgtz
      common /bveg6/ reffic, aeffic, frmax, famax,
     &               fgtz, fez, fltz
c
c declare internal common variables
c
      double precision xjac
      common /monaig/ xjac
c
      integer*8 nd, ng, npg, ndm, nltz, nez, ngtz
      double precision d(ngrid,ndi), di(ngrid,ndi)
      double precision dv2g, fb, f2b
      common /monaip/ d, di,
     &                dv2g, fb, f2b,
     &                nd, ng, npg, ndm, nltz, nez, ngtz
c
      integer*8 kg(ndi), k, ndimen
      double precision dx(ndi), dxg, xnd
      common /monall/ dx, dxg, xnd,
     &                kg, k, ndimen
c
      double precision avgi, sd, chi2a
      common /monsta/ avgi, sd, chi2a
c
c declare local variables
c
      integer*8 istart, i, j
      double precision xin(ngrid), rc, xn, dr, xo
c
c set default values
c
      data nprn /0/		      !default print cumulative information
      data ndev /6/		      !default output device
      data xl /ndi*0.0d0/	      !default integration lower bounds
      data xu /ndi*1.0d0/	      !default integration upper bounds
      data acc /-1.0d0/		      !default termination accuracy
      data it /0/		      !interations completed
      data ndo /1/		      !number of subdivisions per axis
      data si /0.0d0/		      !sum S / sigma^2
      data swgt /0.0d0/		      !sum 1 / sigma^2
      data schi /0.0d0/		      !sum S^2 / sigma^2
      data xi /nxi*1.0d0/	      !location of i-th division on j-th axis
      data alph /0.875d0/	      !default convergence parameter
      data ndmx /48/		      !default number of grid divisions
      data mds /0/		      !use stratified and/or importance sampling
      data rtype /0/		      !use Sobol quasi-random sequences
c
      data istart /0/
      save

c
c set the two seeds
      if (mod(seed,2).eq.0) then
         seed1 = seed / 2
         seed2 = seed / 2
      else
         seed1 = (seed + 1) / 2
         seed2 = (seed - 1) / 2
      endif

c
c test to see that ndim is less than available from either the
c dimensionality of the grid, or the Sobol generator (40).
c
      if ( ndim .gt. min( ndi, 40 ) ) then
         write(ndev,*)
         write(ndev,*) "MONACO called with ndim > ndi"
         write(ndev,*) "ndim =", ndim, ",  ndi =", ndi
         stop
      end if
c
      if ( istart .ne. 1 ) then
         call imonrn( ndi, seed1, seed2 )
         call imonso( ndi )
         istart = 1
      end if
      ndo = 1
      do j = 1, ndim
 1       xi(1,j) = 1.0d0
      end do
      avgi = 0.0d0
      sd = 0.0d0
      chi2a = 0.0d0
c
      entry monaco_init1(
     &         ndim,		!in:  number of dimensions in hypercube
     &         npt		!in:  number of evaluations per iteration
     &                  )
c
c initializes cumulative variables, but not grid
c
      it = 0
      si = 0.0d0
      swgt = si
      schi = si
c
      entry monaco_init2(
     &         ndim,		!in:  number of dimensions in hypercube
     &         npt		!in:  number of evaluations per iteration
     &                  )
c
c no initialisation
c
      ncall = npt
      nd = ndmx
      ng = 1
      if ( mds .ne. 0 ) then
         ng = int( (dble( ncall ) / 2.0d0)**(1.0d0 / dble( ndim )) )
         mds = 1
         if ( (2 * ng - ndmx) .ge. 0 ) then
            mds = -1
            npg = ng / ndmx + 1
            nd = ng / npg
            ng = npg * nd
         end if
      end if
      k = ng**ndim
      npg = ncall / k
      if ( npg .lt. 2 ) then
         npg = 2
      end if
      calls = dble( npg * k )
      dxg = 1.0d0 / dble( ng )
      dv2g = (calls * dxg**ndim)**2 / (dble( npg )**2 * dble( npg - 1 ))
      xnd = dble( nd )
      ndm = nd - 1
      dxg = dxg * xnd
      xjac = 1.0d0 / calls
      do j = 1, ndim
         dx(j) = xu(j) - xl(j)
         xjac = xjac * dx(j)
      end do
c
c rebin, preserving bin density
c
      if ( nd .ne. ndo ) then
         rc = dble( ndo ) / xnd
         do j = 1, ndim
            k = 0
            xn = 0.0d0
            dr = xn
            i = k
 10         k = k + 1
            dr = dr + 1.0d0
            xo = xn
            xn = xi(k,j)
 20         if ( rc .gt. dr ) goto 10
            i = i + 1
            dr = dr - rc
            xin(i) = xn - (xn - xo) * dr
            if ( i .lt. ndm ) goto 20
            do i = 1, ndm
               xi(i,j) = xin(i)
            end do
            xi(nd,j) = 1.0d0
         end do
         ndo = nd
      end if
c
 200  format( / " MONACO input parameters:" /
     &        " ndim =", i3, ",  ncall =", i8,a, "   rtype =", i2 /
     &   " nprn =", i3, ",  acc = ", 1pe7.1e1, ",  alph = ", 1pe7.2e1,
     &        ",  nd =", i3, ",  mds =", i3   )   !    /
c     &        " integration bounds: (lower, upper)" /
c     &        ( " dimension", i3, ":   ( ", 1pe13.6, ", ", 1pe13.6,
c     &        " )" ) )
      if ( nprn .ge. 0 ) then
         if (calls.gt.1d6) then
            write(ndev,200) ndim, int( calls/1024**2 ),"M,", rtype,
     &                   nprn, acc, alph, nd, mds    !,
c     &                   (j, xl(j), xu(j), j=1, ndim)
         elseif (calls.gt.1d4) then
            write(ndev,200) ndim, int( calls )/1024,"k,", rtype,
     &                   nprn, acc, alph, nd, mds    !,
c     &                   (j, xl(j), xu(j), j=1, ndim)
         else
            write(ndev,200) ndim, int( calls ),",", rtype,
     &                   nprn, acc, alph, nd, mds    !,
c     &                   (j, xl(j), xu(j), j=1, ndim)
         endif
      end if
c
      entry monaco_init3(
     &         ndim,		!in:  number of dimensions in hypercube
     &         npt		!in:  number of evaluations per iteration
     &                  )
c
c setup main integration loop
c note that number of points per iteration remains as defined above
c
      ti = 0.0d0
      tsi = 0.0d0
      nltz = 0
      nez = 0
      ngtz = 0
      frmax = 0.0d0
      famax = 0.0d0
c
      do j = 1, ndim
         kg(j) = 1
         do i = 1, nd
            d(i,j) = ti
            di(i,j) = ti
         end do
      end do
      ndimen = ndim
c
      fb = 0.0d0
      f2b = 0.0d0
      k = 0
c
c done
c
      return
      end
c
c-------------------------------------------------------------------------------
c
      subroutine monaco_get(
     &              r,		!out:  generated point in hypercube
     &              wgt      	!out:  generated weight of point
cerw     &              itn		!out:  iteration number
     &                     )
      implicit none
c
c declare output variables
c
cerw      integer*8 itn
      double precision r(*), wgt
c
c declare global parameters
c
      integer*8 ngrid, ndi, nxi
      parameter( ngrid=48, ndi=24, nxi=ngrid*ndi )
c
c declare external common variables
c
      integer*8 ncall, itmx, nprn, ndev
      double precision xl(ndi), xu(ndi), acc
      common /bveg1/ xl, xu, acc, ncall, itmx, nprn, ndev
c
      integer*8 it, ndo
      double precision si, swgt, schi, xi(ngrid,ndi)
      common /bveg2/ si, swgt, schi, xi, it, ndo
c
      integer*8 rtype
      common /bveg5/ rtype
c
c declare internal common variables
c
      double precision xjac
      common /monaig/ xjac
c
      integer*8 kg(ndi), k, ndimen
      double precision dx(ndi), dxg, xnd
      common /monall/ dx, dxg, xnd,
     &                kg, k, ndimen
c
c declare local variables
c
      integer*8 ia(ndi), j
      double precision rand(ndi), xn, xo, rc
      double precision rand_dum(40)
c
c generate point, with accompanying weight
c
      if ( rtype .eq. 0 ) then
         call monran( rand )
      else if ( rtype .eq. 1 ) then
         call monsob( rand_dum )
         do j = 1, ndi
            rand(j) = rand_dum(j)
         end do
      else
         write(ndev,*)
      write(ndev,*) "MONACO:  invalid random sequence generator choice"
         write(ndev,*) "rtype =", rtype
         stop
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
c
c done
c
      return
      end
c
c-------------------------------------------------------------------------------
c
      subroutine monaco_put(
     &              r,		!in:  generated point in hypercube
     &              wgt,	!in:  weight as generated by "monaco_get"
     &              value,	!in:  value of integrand
     &              rew_grid  !in:  reweighting factor for grid generation in order to improve statistics in certain PS regions
     &                     )
      implicit none
c
c declare input variables
c
      double precision r(*), wgt, value, rew_grid
c
c declare global parameters
c
      integer*8 ngrid, ndi, nxi
      parameter( ngrid=48, ndi=24, nxi=ngrid*ndi )
c
c declare external common variables
c
      integer*8 ncall, itmx, nprn, ndev
      double precision xl(ndi), xu(ndi), acc
      common /bveg1/ xl, xu, acc, ncall, itmx, nprn, ndev
c
      integer*8 it, ndo
      double precision si, swgt, schi, xi(ngrid,ndi)
      common /bveg2/ si, swgt, schi, xi, it, ndo
c
      integer*8 ndmx, mds
      double precision alph
      common /bveg3/ alph, ndmx, mds
c
      double precision calls, ti, tsi
      common /bveg4/ calls, ti, tsi
c
      double precision reffic, aeffic, frmax, famax, fltz, fez, fgtz
      common /bveg6/ reffic, aeffic, frmax, famax,
     &               fgtz, fez, fltz
c
c declare internal common variables
c
      integer*8 nd, ng, npg, ndm, nltz, nez, ngtz
      double precision d(ngrid,ndi), di(ngrid,ndi)
      double precision dv2g, fb, f2b
      common /monaip/ d, di,
     &                dv2g, fb, f2b,
     &                nd, ng, npg, ndm, nltz, nez, ngtz
c
      integer*8 kg(ndi), k, ndimen
      double precision dx(ndi), dxg, xnd
      common /monall/ dx, dxg, xnd,
     &                kg, k, ndimen
c
      double precision avgi, sd, chi2a
      common /monsta/ avgi, sd, chi2a
c
c declare local variables
c
      integer*8 ia(ndi), i, j
      double precision dt(ndi), xin(ngrid), ri(ngrid)
      double precision fret, f, f2, ti2, xo, xn, rc, dr
c
c add point
c
      k = k + 1
      fret = value
      frmax = max( frmax, abs( fret ))
      if ( fret .ne. 0.0d0 ) then
         if ( fret .gt. 0.0d0 ) then
            ngtz = ngtz + 1
         else
            nltz = nltz + 1
         end if
      else
         nez = nez + 1
      end if
c
      f = wgt * fret
      famax = max( famax, abs( f ))
      f2 = f * f
      fb = fb + f
      f2b = f2b + f2

      ! modify f and f2 by rew_grid in order to improve statistics in certain PS regions
      f = f * rew_grid
      f2 = f * f

c
      do j = 1, ndimen
         rc = ( r(j) - xl(j) ) / dx(j)
         ia(j) = 1
 10      if ( (rc .ge. xi(ia(j),j)) .and. (ia(j) .lt. ngrid) ) then
            ia(j) = ia(j) + 1
            goto 10
         end if
         di(ia(j),j) = di(ia(j),j) + f
         if ( mds .ge. 0 ) then
            d(ia(j),j) = d(ia(j),j) + f2
         end if
      end do
      if ( k .lt. npg ) then
         return
      end if
c
      f2b = sqrt( f2b * dble( npg ))
      f2b = (f2b - fb) * (f2b + fb)
      ti = ti + fb
      tsi = tsi + f2b
      if ( mds .lt. 0 ) then
         do j = 1, ndimen
            d(ia(j),j) = d(ia(j),j) + f2b
         end do
      end if
      k = ndimen
 20   kg(k) = mod( kg(k), ng ) + 1
      if ( kg(k) .ne. 1 ) then
         fb = 0.0d0
         f2b = 0.0d0
         k = 0
         return
      end if
      k = k - 1
      if ( k .gt. 0 ) goto 20
c
c compute final results for this iteration
c
      if ( frmax .ne. 0.0d0 ) then
         reffic = abs( ti ) / frmax
      else
         reffic = 0.0d0
      end if
      if ( famax .ne. 0.0d0 ) then
         aeffic = abs( ti ) / (famax * calls)
      else
         aeffic = 0.0d0
      end if
      fltz = dble( nltz ) / calls
      fez = dble( nez ) / calls
      fgtz = dble( ngtz ) / calls
c
      tsi = tsi * dv2g
      ti2 = ti * ti
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
c
      it = it + 1
 201  format( / " iteration", i3, ":" /
     &        " integral = ", 1pe13.6, ",  sigma = ", 1pe9.3 /
     &        " efficacy = ", 1pe9.3, ",  raw efficacy = ", 1pe9.3 /
     &        " f_positive = ", 1pe8.3e1, ",  f_zero = ", 1pe8.3e1,
     &        ",  f_negative = ", 1pe8.3e1 /
     &        " accumulated statistics:" /
     &        " integral = ", 1pe13.6, ",  sigma = ", 1pe9.3,
     &        ",  chi^2/iteration = ", 1pe9.3 )
 202  format( / " grid data for axis", i3, ":" /
     &        5x, "x", 6x, "delta_i", 8x, "x", 6x, "delta_i",
     &        8x, "x", 6x, "delta_i", 8x, "x", 6x, "delta_i",
     &        8x, "x", 6x, "delta_i", 8x, "x", 6x, "delta_i" /
     &        (1x, 1pe9.3,1x, 1pe9.3, 3x, 1pe9.3,1x, 1pe9.3, 3x,
     &        1pe9.3,1x, 1pe9.3, 3x, 1pe9.3,1x, 1pe9.3, 3x,
     &        1pe9.3,1x, 1pe9.3, 3x, 1pe9.3,1x, 1pe9.3 ) )
      if ( nprn .ge. 0 ) then
         if ( tsi .gt. 0.0d0 ) then
            tsi = sqrt( tsi )
         else
            tsi = 0.0d0
         end if
         write(ndev,201) it,
     &      ti, tsi,
     &      aeffic, reffic,
     &      fgtz, fez, fltz,
     &      avgi, sd, abs( chi2a )
         if ( nprn .gt. 0 ) then
            do j = 1, ndimen
               write(ndev,202) j, (xi(i,j), di(i,j), i = 1+nprn/2, nd)
            end do
         end if
      end if
c
c refine grid
c
      do j = 1, ndimen
         xo = d(1,j)
         xn = d(2,j)
         d(1,j) = (xo + xn) / 2.0d0
         dt(j) = d(1,j)
         do i = 2, ndm
            d(i,j) = xo + xn
            xo = xn
            xn = d(i+1,j)
            d(i,j) = (d(i,j) + xn) / 3.0d0
            dt(j) = dt(j) + d(i,j)
         end do
         d(nd,j) = (xn + xo) / 2.0d0
         dt(j) = dt(j) + d(nd,j)
      end do
c
      do j = 1, ndimen
         rc = 0.0d0
         do i = 1, nd
            ri(i) = 0.0d0
            if ((d(i,j).ge.0d0).and.(dt(j).ne.0d0)) then
               d(i,j) = max( 1.0d-30, d(i,j) )
               xo = dt(j) / d(i,j)
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
 30      k = k + 1
         dr = dr + ri(k)
         xo = xn
         xn = xi(k,j)
 40      if ( rc .gt. dr ) goto 30
         i = i + 1
         dr = dr - rc
         xin(i) = xn - (xn - xo) * dr / ri(k)
         if ( i .lt. ndm ) goto 40
         do i = 1, ndm
            xi(i,j) = xin(i)
         end do
         xi(nd,j) = 1.0d0
      end do
c
      ti = 0.0d0
      tsi = 0.0d0
      nltz = 0
      nez = 0
      ngtz = 0
      frmax = 0.0d0
      famax = 0.0d0
c
      do j = 1, ndimen
         kg(j) = 1
         do i = 1, nd
            d(i,j) = ti
            di(i,j) = ti
         end do
      end do
c
      fb = 0.0d0
      f2b = 0.0d0
      k = 0
c
c done
c
      return
      end
c
c-------------------------------------------------------------------------------
c
      subroutine monaco_result(
     &              mean,	!out:  function integral
     &              sdev,	!out:  standard deviation of integral
     &              chi2	!out:  chi^2 per degree of freedom
     &                        )
      implicit none
c
c declare output variables
c
      double precision mean, sdev, chi2
c
c declare external global parameters
c
      integer*8 ngrid, ndi, nxi
      parameter( ngrid=48, ndi=24, nxi=ngrid*ndi )
c
c declare internal global parameters
c
      double precision avgi, sd, chi2a
      common /monsta/ avgi, sd, chi2a
c
c declare common variables
c
      integer*8 ncall, itmx, nprn, ndev
      double precision xl(ndi), xu(ndi), acc
      common /bveg1/ xl, xu, acc, ncall, itmx, nprn, ndev
c
      integer*8 it, ndo
      double precision si, swgt, schi, xi(ngrid,ndi)
      common /bveg2/ si, swgt, schi, xi, it, ndo
c
      integer*8 ndmx, mds
      double precision alph
      common /bveg3/ alph, ndmx, mds
c
      double precision calls, ti, tsi
      common /bveg4/ calls, ti, tsi
c
      integer*8 rtype
      common /bveg5/ rtype
c
      double precision reffic, aeffic, frmax, famax, fltz, fez, fgtz
      common /bveg6/ reffic, aeffic, frmax, famax,
     &               fgtz, fez, fltz
c
c transfer values
c
      mean = avgi
      sdev = sd
      chi2 = abs( chi2a )
c
c done
c
      return
      end
c
c-------------------------------------------------------------------------------
c
      subroutine monaco_read(file_name)
      implicit none
c
c declare global parameters
c
      character*250 file_name
      integer*8 ngrid, ndi
      parameter( ngrid=48, ndi=24 )
c
c declare local variables
c
      integer*8 i, j, k
c
c declare common variables
c
      integer*8 it, ndo
      double precision si, swgt, schi, xi(ngrid,ndi)
      common /bveg2/ si, swgt, schi, xi, it, ndo
c
c open file and read grid
c
      open( unit=15, file=file_name, access="sequential", 
     *      err=20, status="old" )
 10   format( 3( 1x, 1pd23.16 ) )
      do j = 1, ndi
         do i = 0, (ngrid-1)/3
            read(unit=15,fmt=10,err=30) ( xi(3*i+k,j), k=1, 3 )
         end do
      end do
      close( unit=15 )
c
      write(6,*)
      write(6,*) "MONACO:  grid read from file ",file_name
      return
c
c otherwise, signal open error (non-fatal)
c
 20   write(6,*)
c      write(6,*) "MONACO:  open error on file unit 15" //
c     &           " - continuing with uniform grid"
      write(6,*) "MONACO:  continuing with uniform grid"
      return
c
c otherwise, signal read error (fatal)
c
 30   write(6,*)
      write(6,*) "MONACO:  read error on file unit 15"
      stop
c
c done
c
      end
c
c-------------------------------------------------------------------------------
c
      subroutine monaco_write(file_name)
      implicit none
c
c declare global parameters
c
      character*250 file_name
      integer*8 ngrid, ndi
      parameter( ngrid=48, ndi=24 )
c
c declare local variables
c
      integer*8 i, j, k
c
c declare common variables
c
      integer*8 it, ndo
      double precision si, swgt, schi, xi(ngrid,ndi)
      common /bveg2/ si, swgt, schi, xi, it, ndo
c
c open file and write grid
c
      open( unit=16, file=file_name, access="sequential", 
     *      err=20, status="unknown" )
 10   format( 3( 1x, 1pd23.16 ) )
      do j = 1, ndi
         do i = 0, (ngrid-1)/3
            write(unit=16,fmt=10,err=30) ( xi(3*i+k,j), k=1, 3 )
         end do
      end do
      close( unit=16 )
c
      write(6,*)
      write(6,*) "MONACO:  grid written to file ",
     &            file_name(1:index(file_name,"  ")-1)
      return
c
c otherwise, signal open error (non-fatal)
c
 20   write(6,*)
      write(6,*) "MONACO:  open error on file unit ",
     &            file_name(1:index(file_name,"  ")-1)
      return
c
c otherwise, signal write error (non-fatal)
c
 30   write(6,*)
      write(6,*) "MONACO:  write error on file unit ",
     &            file_name(1:index(file_name,"  ")-1)
      close( unit=16 )
      return
c
c done
c
      end
c
c-------------------------------------------------------------------------------
c
c	This subroutine is a "universal" random number generator
c	as proposed by Marsaglia and Zaman in report FSU-SCRI-87-50
c	Slightly modified by F. James, 1988, to generate a vector
c	of pseudorandom numbers "rvec" of length "len", and also
c	by A. Duff, 1991.
c
c	To get the values in the Marsaglia-Zaman paper, put
c	ij = 1802, kl = 9373
c
      subroutine imonrn( ndim, seed1, seed2)
c
      implicit none
      integer*4 seed1, seed2
      integer*4 ndim, nd
      integer*8 i, j, k, l, m, ij, kl, ii, jj, i1, i2
      double precision s, t, u(97), c, cd, cm
c
      common /monrnc/ u, c, cd, cm, i1,i2, nd
      save
c
c      data ij /1802/, kl /9373/
      ij = 1802 + seed1
      kl = 9373 + seed2

c
      nd = ndim
      i = mod( ij / 177, 177 ) + 2
      j = mod( ij, 177 ) + 2
      k = mod( kl / 169, 178 ) + 1
      l = mod( kl, 169 )
c
      do ii = 1, 97
         s = 0.0d0
         t = 0.5d0
         do jj = 1, 24
            m = mod( mod( i * j, 179 ) * k, 179 )
            i = j
            j = k
            k = m
            l = mod( 53 * l + 1, 169 )
            if ( mod( l * m, 64 ) .ge. 32 ) then
               s = s + t
            end if
            t = 0.5d0 * t
         end do
         u(ii) = s
      end do
c
      c = 362436.0d0 / 16777216.0d0
      cd = 7654321.0d0 / 16777216.0d0
      cm = 16777213.0d0 / 16777216.0d0
c
      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine monran( rvec )
c
      implicit none
      integer*8 i, j, ivec
      integer*4 ndim
      double precision rvec(*), uni
      double precision u(97), c, cd, cm
c
      common /monrnc/ u, c, cd, cm, i, j, ndim
      save
c
      data i /97/, j /33/
c
      do ivec = 1, ndim
         uni = u(i) - u(j)
         if ( uni .lt. 0.0d0 ) then
            uni = uni + 1.0d0
         end if
         u(i) = uni
         i = i - 1
         if ( i .le. 0 ) then
            i = 97
         end if
         j = j - 1
         if ( j .le. 0 ) then
            j = 97
         end if
         c = c - cd
         if ( c .lt. 0.0d0 ) then
            c = c + cm
         end if
         uni = uni - c
         if ( uni .lt. 0.0d0 ) then
            uni = uni + 1.0d0
         end if
         rvec(ivec) = uni
      end do
c
      return
      end
c
c--------- allow for restart of random numbers at predetermined point
c
      subroutine monran_set(id)
      implicit none
      integer*4 ndim, id, ii
      integer*8 i, j
      double precision u(97), c, cd, cm
      common /monrnc/ u, c, cd, cm, i, j, ndim
c
      integer*4 ndims
      integer*8 is, js
      double precision us(97), cs, cds, cms
c
      common /monrncsave/ us, cs, cds, cms, is, js, ndims
      if (id.eq.1) then   ! save to monrncsave
         do ii = 1,97 
            us(ii) = u(ii)
         enddo
         cs  = c
         cds = cd
         cms = cm
         ndims = ndim
         is = i
         js = j
         print*," random number setting in monaco saved "
      elseif(id.eq.2) then   ! restore monrnc from  monrncsave
         do ii = 1,97 
            u(ii) = us(ii)
         enddo
         c  = cs
         cd = cds
         cm = cms
         ndim = ndims
         i = is
         j = js
         print*," random number setting in monaco restored to saved "
      endif
      end
c
c---------------------------------------------------------------------------
c
c	This is a stripped down version of Algorithm 659
c	as appeared in ACM Trans. Math. Software, vol.14, no.1
c	March 1988, pgs. 88-100
c
      subroutine imonso( ndim )
c
      implicit none
      integer*8 v(40,30), s, maxcol, count, poly(40), lastq(40)
c_unused      integer*8 atmost, i, j, j2, k, l, m, newv, exor
      integer*8 atmost, i, j, j2, k, l, m, newv
      integer in
c_unused      integer*8 tau(13)
      integer*8 temp
      integer*4 ndim
c_unused      logical flag(2), includ(8)
      logical  includ(8)
      double precision recipd
      equivalence( poly, lastq )
c
      common /monsoc/ recipd, v, s, maxcol, count, poly
      save
c
      data atmost /1073741823/
      data poly /1, 3, 7, 11, 13, 19, 25, 37, 59, 47, 61, 55, 41,
     &           67, 97, 91, 109, 103, 115, 131, 193, 137, 145, 143,
     &           241, 157, 185, 167, 229, 171, 213, 191, 253, 203, 211,
     &           239, 247, 285, 369, 299/
      data (v(in,1), in=1, 40) /40*1/
      data (v(in,2), in=3, 40) /1, 3, 1, 3, 1, 3, 3, 1, 3, 1, 3, 1, 3,
     &                        1, 1, 3, 1, 3, 1, 3, 1, 3, 3, 1, 3, 1,
     &                        3, 1, 3, 1, 1, 3, 1, 3, 1, 3, 1, 3/
      data (v(in,3), in=4, 40) /7, 5, 1, 3, 3, 7, 5, 5, 7, 7, 1, 3, 3,
     &                        7, 5, 1, 1, 5, 3, 3, 1, 7, 5, 1, 3, 3,
     &                        7, 5, 1, 1, 5, 7, 7, 5, 1, 3, 3/
      data (v(in,4), in=6, 40) /1, 7, 9, 13, 11, 1, 3, 7, 9, 5, 13, 13,
     &                        11, 3, 15, 5, 3, 15, 7, 9, 13, 9, 1, 11,
     &                        7, 5, 15, 1, 15, 11, 5, 3, 1, 7, 9/
      data (v(in,5), in=8, 40) /9, 3, 27, 15, 29, 21, 23, 19, 11, 25, 7,
     &                        13, 17, 1, 25, 29, 3, 31, 11, 5, 23, 27,
     &                        19, 21, 5, 1, 17, 13, 7, 15, 9, 31, 9/
      data (v(in,6), in=14, 40) /37, 33, 7, 5, 11, 39, 63, 27, 17, 15,
     &                         23, 29, 3, 21, 13, 31, 25, 9, 49, 33,
     &                         19, 29, 11, 19, 27, 15, 25/
      data (v(in,7), in=20, 40) /13, 33, 115, 41, 79, 17, 29, 119, 75,
     &                         73, 105, 7, 59, 65, 21, 3, 113, 61, 89,
     &                         45, 107/
      data (v(in,8), in=38, 40) /7, 23, 39/
c      data tau /0, 0, 1, 3, 5, 8, 11, 15, 19, 23, 27, 31, 35/
c
      s = ndim
      i = atmost
      maxcol = 0
 10   maxcol = maxcol + 1
         i = i / 2
         if ( i .gt. 0 ) goto 10
c
      do i = 1, maxcol
         v(1,i) = 1
      end do
c
      do i = 2, s
         j = poly(i)
         m = 0
 20      j = j / 2
         if ( j .gt. 0 ) then
            m = m + 1
            goto 20
         end if
c
         j = poly(i)
         do k = m, 1, -1
            j2 = j / 2
            includ(k) = (j .ne. (2 * j2))
            j = j2
         end do
c
         do j = m + 1, maxcol
            newv = v(i,j-m)
            l = 1
            do k = 1, m
               l = 2 * l
               if ( includ(k) ) then
                  temp = l * v(i,j-k)
                  newv = ieor( newv, temp )	  !assume VAX xor function
c                 newv = exor( newv, temp )	  !only if necessary
               end if
            end do
            v(i,j) = newv
         end do
      end do
c
      l = 1
      do j = maxcol - 1, 1, -1
         l = 2 * l
         do i = 1, s
            v(i,j) = v(i,j) * l
         end do
      end do
c
      recipd = 1.0d0 / (2 * l)
c
      count = 0
      do i = 1, s
         lastq(i) = 0
      end do
c
      return
      end
c
c-------------------------------------------------------------------------------
c
      subroutine monsob( rvec )
c
      implicit none
      integer*8 v(40,30), s, maxcol, count, lastq(40)
c_unused      integer*8 i, i2, l, exor
      integer*8 i, i2, l
      double precision recipd, rvec(40)
c
      common /monsoc/ recipd, v, s, maxcol, count, lastq
      save
c
      l = 0
      i = count
 10   l = l + 1
      i2 = i / 2
      if ( i .ne. (2 * i2) ) then
         i = i2
         goto 10
      end if
c
      if ( l .gt. maxcol ) then
         stop "MONACO:  Sobol generator - too many calls."  !not too likely!
      end if
c
      do i = 1, s
         lastq(i) = ieor( lastq(i), v(i,l) )	  !again - VAX Fortran specific
c        lastq(i) = exor( lastq(i), v(i,l) )	  !use only if necessary
         rvec(i) = lastq(i) * recipd
      end do
c
      count = count + 1
c
      return
      end
c
c-------------------------------------------------------------------------------
c
      integer*8 function exor( iin,
     &                         jin )
c
      implicit none
      integer*8 i, j, k, l, iin, jin, i2, j2
c
      i = iin
      j = jin
      k = 0
      l = 1
c
 10   if ( (i .eq. 0) .and. (j .eq. 0) ) then
         exor = k
         return
      end if
c
      i2 = i / 2
      j2 = j / 2
      if ( (i .eq. (2 * i2)) .neqv. (j .eq. (2 * j2)) ) then
         k = k + l
      end if
      i = i2
      j = j2
      l = 2 * l
      goto 10
c
      end

c-------------------------------------------------------------
c
c	This subroutine is a "universal" random number generator
c	as proposed by Marsaglia and Zaman in report FSU-SCRI-87-50
c	Slightly modified by F. James, 1988, to generate a vector
c	of pseudorandom numbers "rvec" of length "len", and also
c	by A. Duff, 1991.
c
c	To get the values in the Marsaglia-Zaman paper, put
c	ij = 1802, kl = 9373
c
      subroutine iranmr( ndim, seed1, seed2 )
c
      implicit none
      integer*4 ndim, nd, seed1, seed2
      integer*8 i, j, k, l, m, ij, kl, ii, jj
      double precision s, t, u(97), c, cd, cm
c
      common /comrmr/ u, c, cd, cm, nd
      save /comrmr/
c
c     data ij /1802/, kl /9373/
      ij = 1802 + seed2
      kl = 9373 + seed1

      nd = ndim
      i = mod( ij / 177, 177 ) + 2
      j = mod( ij, 177 ) + 2
      k = mod( kl / 169, 178 ) + 1
      l = mod( kl, 169 )
c
      do ii = 1, 97
         s = 0.0d0
         t = 0.5d0
         do jj = 1, 24
            m = mod( mod( i * j, 179 ) * k, 179 )
            i = j
            j = k
            k = m
            l = mod( 53 * l + 1, 169 )
            if ( mod( l * m, 64 ) .ge. 32 ) then
               s = s + t
            end if
            t = 0.5d0 * t
         end do
         u(ii) = s
      end do
c
      c = 362436.0d0 / 16777216.0d0
      cd = 7654321.0d0 / 16777216.0d0
      cm = 16777213.0d0 / 16777216.0d0
c
      return
      end
c
c-------------------------------------------------------------------------------
c
      subroutine ranmar( rvec )
c
      implicit none
      integer*8 i, j, ivec
      integer*4 ndim
      double precision rvec(*), uni
      double precision u(97), c, cd, cm
c
      common /comrmr/ u, c, cd, cm, ndim
      save /comrmr/
c
      data i /97/, j /33/
c
      do ivec = 1, ndim
         uni = u(i) - u(j)
         if ( uni .lt. 0.0d0 ) then
            uni = uni + 1.0d0
         end if
         u(i) = uni
         i = i - 1
         if ( i .le. 0 ) then
            i = 97
         end if
         j = j - 1
         if ( j .le. 0 ) then
            j = 97
         end if
         c = c - cd
         if ( c .lt. 0.0d0 ) then
            c = c + cm
         end if
         uni = uni - c
         if ( uni .lt. 0.0d0 ) then
            uni = uni + 1.0d0
         end if
         rvec(ivec) = uni
      end do
c
      return
      end



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

      integer mu, i, expo

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
      end



