c #define QUAD_ENHANCED
c***************************************************************************
      real*8 function mjj2Q(p1,p2)
c***************************************************************************
c
c       This function calculates the invariant mass
c       square of the two four-momenta as given.
c       ie.,  mjj2 = (p1 + p2)^2
c
c       input: real*8 p1(0:3), p2(0:3) four momenta
c
c***************************************************************************
      implicit none
c declare input/output variables
      real*8 p1(0:3), p2(0:3)
c calculate
#ifdef QUAD_ENHANCED
      mjj2Q = ( real( p1(0) + p2(0),16 )**2
     &          -  real( p1(1) + p2(1),16 )**2
     &          -  real( p1(2) + p2(2),16 )**2
     &          -  real( p1(3) + p2(3),16 )**2 )
#else
      mjj2Q = ( ( p1(0) + p2(0) )**2
     &       - ( p1(1) + p2(1) )**2
     &       - ( p1(2) + p2(2) )**2
     &       - ( p1(3) + p2(3) )**2 )
      return
#endif
      end


c***************************************************************************
      real*8 function qjj2Q(p1,p2)
c***************************************************************************
c
c       This function calculates the invariant momentum
c       transfer square of the two four-momenta as given.
c       ie.,  qjj2 = (p1 - p2)^2
c
c       input: real*8 p1(0:3), p2(0:3) four momenta
c
c***************************************************************************
      implicit none

c declare input/output variables

      real*8 p1(0:3), p2(0:3)

c calculate
#ifdef QUAD_ENHANCED
      qjj2Q = ( real( p1(0) - p2(0),16 )**2
     &       - real( p1(1) - p2(1),16 )**2
     &       - real( p1(2) - p2(2),16 )**2
     &       - real( p1(3) - p2(3),16 )**2 )
#else
      qjj2Q = ( ( p1(0) - p2(0) )**2
     &       - ( p1(1) - p2(1) )**2
     &       - ( p1(2) - p2(2) )**2
     &       - ( p1(3) - p2(3) )**2 )
#endif
      return
      end


c***************************************************************************
      function dotrrQ(p1,p2)
c***************************************************************************
c
c     dotrr(p1,p2) = p1.p2
c
c***************************************************************************
      implicit none

      real*8 dotrrQ,p1(0:3),p2(0:3)
#ifdef QUAD_ENHANCED
      dotrrQ = real(p1(0),16)*real(p2(0),16) - real(p1(1),16)*real(p2(1),16) -
     $         real(p1(2),16)*real(p2(2),16) - real(p1(3),16)*real(p2(3),16)
#else
      dotrrQ = p1(0)*p2(0) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)
#endif
      end

c***************************************************************************
      real*8 FUNCTION mass2q(p)
c***************************************************************************
c     calculates p^2
c***************************************************************************
      real*8 p(0:3)
#ifdef QUAD_ENHANCED
      mass2q= reaL(p(0),16)**2-real(p(1),16)**2-real(p(2),16)**2-real(p(3),16)**2
#else
      mass2q= p(0)*p(0) - p(1)*p(1) - p(2)*p(2) - p(3)*p(3)
#endif
      return
      end

#ifdef WITH_QUAD
c***************************************************************************
      real*16 function qjj2Quad(p1,p2)
c***************************************************************************
c
c       This function calculates the invariant momentum
c       transfer square of the two four-momenta as given.
c       ie.,  qjj2 = (p1 - p2)^2
c
c       input: real*8 p1(0:3), p2(0:3) four momenta
c
c***************************************************************************
      implicit none

c declare input/output variables

      real*16 p1(0:3), p2(0:3)

c calculate

      qjj2Quad = ( real( p1(0) - p2(0),16 )**2
     &       - real( p1(1) - p2(1),16 )**2
     &       - real( p1(2) - p2(2),16 )**2
     &       - real( p1(3) - p2(3),16 )**2 )

      return
      end


c***************************************************************************
      function dotrrQuad(p1,p2)
c***************************************************************************
c
c     dotrr(p1,p2) = p1.p2
c
c***************************************************************************
      implicit none

      real*16 dotrrQuad,p1(0:3),p2(0:3)

      dotrrQuad = p1(0)*p2(0) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)

      end

c***************************************************************************
      complex*32 function dotccQuad(v1,v2)
c***************************************************************************
c
c     dotcc(v1,v2) = v1.v2
c
c***************************************************************************

      implicit none
      complex*32 v1(0:3), v2(0:3)
      dotccQuad = v1(0)*v2(0)-v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)
      end

c***************************************************************************
      complex*32 function dotcrQuad(v1,v2)
c***************************************************************************
c
c     dotrc(v1,v2) = v1.v2
c
c***************************************************************************

      implicit none
      real*16 v2(0:3)
      complex*32  v1(0:3)
      dotcrQuad = v1(0)*v2(0)-v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)
      end

c***************************************************************************
      complex*32 function dotrcquad(v1,v2)
c***************************************************************************
c
c     dotrc(v1,v2) = v1.v2
c
c***************************************************************************

      implicit none
      real*16 v1(0:3)
      complex*32  v2(0:3)
      dotrcquad = v1(0)*v2(0)-v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)
      end




c***************************************************************************
      real*16 FUNCTION mass2quad(p)
c***************************************************************************
c     calculates p^2
c***************************************************************************
      real*16 p(0:3)

      mass2quad= p(0)*p(0)-p(1)*p(1)-p(2)*p(2)-p(3)*p(3)

      return
      end



c***************************************************************************
      real*16 function mjj2Quad(p1,p2)
c***************************************************************************
c
c       This function calculates the invariant mass
c       square of the two four-momenta as given.
c       ie.,  mjj2 = (p1 + p2)^2
c
c       input: real*8 p1(0:3), p2(0:3) four momenta
c
c***************************************************************************
      implicit none
c declare input/output variables
      real*16 p1(0:3), p2(0:3)
c calculate
      mjj2QUAD = ( real( p1(0) + p2(0),16 )**2
     &          -  real( p1(1) + p2(1),16 )**2
     &          -  real( p1(2) + p2(2),16 )**2
     &          -  real( p1(3) + p2(3),16 )**2 )
      return
      end

#endif