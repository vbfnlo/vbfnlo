

      subroutine vcartx(p,vmass,vwidth,ncart,nsv , vc)
c
c This subroutine computes an effective VECTOR wavefunction for an internal
c vector boson line. The propagator is inlcuded in the Feynman gauge.
c
c input:
c       real    p(0:3)         : four-momentum of vector boson
c       real    vmass          : mass          of vector boson
c       integer ncart = 0,1,2,3: cartesian polarization direction
c                                of vector boson
c       integer nsv  = -1 or 1 : +1 for final, -1 for initial
c
c output:
c       complex vc(6)          : vector wavefunction       epsilon^mu(v)
c     
      implicit none
      double complex vc(6), d
      double precision p(0:3),vmass,vwidth,q2
      integer ncart, mu, nsv

      q2 = p(0)**2-p(1)**2-p(2)**2-p(3)**2
      if (vmass.eq.0d0) then
         d = 1d0/q2
      else
         d = 1d0/dcmplx( q2-vmass**2, vmass*vwidth )
      endif
      do mu = 0,3
         vc(mu+1) = dcmplx(0d0,0d0)
      enddo
      if (ncart.eq.0) then
         vc(1) = d
      else
         vc(ncart+1) = -d
      endif
c
      vc(5) = dcmplx(p(0),p(3))*nsv
      vc(6) = dcmplx(p(1),p(2))*nsv
c
      return
      end 

********************************************************************************

      subroutine vcartxdp(p,ncart,nsv , vc)
c     
      implicit none
      double complex vc(6), d
      double precision p(0:3)
      integer ncart, mu, nsv

      d = 1d0
      do mu = 0,3
         vc(mu+1) = dcmplx(0d0,0d0)
      enddo
      if (ncart.eq.0) then
         vc(1) = d
      else
         vc(ncart+1) = -d
      endif
c
      vc(5) = dcmplx(p(0),p(3))*nsv
      vc(6) = dcmplx(p(1),p(2))*nsv
c
      return
      end
