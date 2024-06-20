      subroutine vvhaxx(v1,v2,sc,bos1,bos2,vertex)
c
c This subroutine computes an amplitude of the vector-vector-scalar
c coupling using anomalous Higgs couplings
c
c Adapted from vvsxxx by Michael Rauch, July 2011
c
c input:
c       complex v1(6)          : first  vector                        v1
c       complex v2(6)          : second vector                        v2
c       complex sc(3)          : input  scalar                        s
c       integer bos1, bos2     : type of bosons attached
c
c output:
c       complex vertex         : amplitude                gamma(v1,v2,s)
c     
      implicit none
      double complex v1(6),v2(6),sc(3),vertex
      integer bos1, bos2
c
      double precision p1(0:4), p2(0:4), ps(0:4), dotrr
      double complex dotrc, dotcc, epsccrr, I3point
      external dotrr, dotrc, dotcc, epsccrr, I3point
c     koppln common blocks
      double complex ahvv(3,4,4), ahvvL(3,4,4)
      common/tensorhvv/ ahvv, ahvvL
c     formfactors
      integer ffac
      double precision m2ff, mff 
      double complex formfactor
      logical lff
      common/ formfacmass/ m2ff,mff,ffac,lff


      p1(0) = dreal(v1(5))
      p1(1) = dreal(v1(6))
      p1(2) = dimag(v1(6))
      p1(3) = dimag(v1(5))
      p1(4) = dotrr(p1,p1)
      p2(0) = dreal(v2(5))
      p2(1) = dreal(v2(6))
      p2(2) = dimag(v2(6))
      p2(3) = dimag(v2(5))
      p2(4) = dotrr(p2,p2)
      ps(0) = dreal(sc(2))
      ps(1) = dreal(sc(3))
      ps(2) = dimag(sc(3))
      ps(3) = dimag(sc(2))
      ps(4) = dotrr(ps,ps)

      if(lff) then
         if(ffac.eq.1) then
            formfactor = M2ff/(-p1(4)+M2ff)*M2ff/(-p2(4)+M2ff)
         else if(ffac.eq.2) then
            formfactor = -2*m2ff*I3point(mff,p1(4),p2(4),ps(4))
         else
            write(*,*) "Please select ffac=1 or ffac=2 in anom_HVV.dat"
            stop
         endif
      else
        formfactor = 1d0
      endif

      vertex = sc(1)*formfactor*(
     &           ahvv(1,bos1,bos2)*dotcc(v1,v2)
     &         + ahvv(2,bos1,bos2)
     &           *( dotrr(p1,p2)*dotcc(v1,v2)
     &            - dotrc(p2,v1)*dotrc(p1,v2) )
     &         + ahvv(3,bos1,bos2)*epsccrr(v1,v2,p1,p2)
     &         )
c
      return
      end
