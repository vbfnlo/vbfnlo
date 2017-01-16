c Author: Francisco Campanario
c Date: 21.03.2007
c Description: Input of divergent part of b,C,d,E functions

c #define UVPOLE_ONLY

      function B0finDiv(p1sq,musq,i)
      implicit none
      complex * 16 B0finDiv
      real * 8 p1sq,musq
      real * 8 tiny
      parameter (tiny=1d-6)
      Integer i
      if (musq.lt.0) then
         write(*,*) 'ERROR in B0fin: mu^2 MUST be a positive number'
         stop
      endif

#ifdef UVPOLE_ONLY
      if(i.eq.1) then
	B0finDiv = (1d0,0d0)
      else
	B0finDiv = 0d0
      endif
      return
#endif

      if (abs(p1sq).lt.tiny) then
      B0finDiv = (0d0,0d0)
      return
      ELSE
      If(i.eq.1)  then 

      B0finDiv = (1d0,0d0)
      RETURN
      else

      B0finDiv=(0d0,0d0)
      RETURN
      endif
      
      endif

      end


      function C0finDiv(p1sq,p2sq,s12,musq,i)
      implicit none
      complex*16 C0finDiv,C0finDiv_1m,C0finDiv_2m
      external C0finDiv_1m,C0finDiv_2m
      real*8 p1sq,p2sq,s12,musq
      integer i    
      real*8 tiny
      parameter (tiny=1d-6)

#ifdef UVPOLE_ONLY
      if(i.eq.1) then
	C0finDiv = (1d0,0d0)
      else
	C0finDiv = 0d0
      endif
      return
#endif

      if(abs(p1sq).lt.tiny) then
      if(abs(p2sq).lt.tiny) then
      if(abs(s12).lt.tiny) then
      C0finDiv=(0d0,0d0)
      return
      else
      C0finDiv=C0finDiv_1m(s12,musq,i) 
      return
      endif !s12
      else !p2sq
      if(abs(s12).lt.tiny) then
      C0finDiv=C0finDiv_1m(p2sq,musq,i)
      return 
      else
      C0finDiv=C0finDiv_2m(p2sq,s12,musq,i) 
      return
      endif !s12
      endif  !p2sq
      else
      if(abs(p2sq).lt.tiny) then
      if(abs(s12).lt.tiny) then
      C0finDiv=C0finDiv_1m(p1sq,musq,i) 
      return
      else
      C0finDiv=C0finDiv_2m(p1sq,s12,musq,i)
      return 
      endif !s12
      else 
      if(abs(s12).lt.tiny) then
      C0finDiv=C0finDiv_2m(p1sq,p2sq,musq,i) 
      return
      else
      C0finDiv=(0d0,0d0)
      return
      endif !s12
      endif  !p2sq
      endif

      end


      function C0finDiv_1m(p1sq,musq,i)
      implicit none
      complex*16 C0finDiv_1m
      real*8 p1sq,musq
      complex * 16 ipi,l1
      parameter (ipi=(0,3.141592653589793238462643383279502884197D0))     
      integer i

       if (i.eq.1) then

       l1=log(abs(p1sq)/musq)
      
       If(p1sq.gt.0.d0) l1=l1-ipi
      
       C0finDiv_1m=(-l1)/p1sq
       return
       else
      
       C0finDiv_1m=1/p1sq
       return
      
       endif

       End


      function C0finDiv_2m(p1sq,p2sq,musq,i)
      implicit none
      complex*16 C0finDiv_2m
      real*8 p1sq,p2sq,musq
      complex * 16 ipi,l1
      parameter (ipi=(0,3.141592653589793238462643383279502884197D0))     
      integer i


       if(i.eq.1) then

       l1=log(abs(p1sq)/abs(p2sq))

       if (p1sq.gt.0) then
       if (p2sq.lt.0) l1 = l1 - ipi
       else
       if (p2sq.gt.0) l1 = l1 + ipi
       endif

       C0finDiv_2m=l1/(p2sq-p1sq)
       return
       else
       C0finDiv_2m=(0d0,0d0)
       return
       endif
       End


c    ***********   D00m_fin(s,t,musq)   *************
c
c Scalar box with MASSLESS PROPAGATORS and with external kinematics: 
c
c                     k
c     p1 ->-----------<-------------<-- p4
c             |                 |
c             |                 |
c             |                 |  
c     p2 ->-------------------------<-- p3
c
c
c    s = (p1+p2)^2
c    t = (p2+p3)^2
c    p1^2 = 0, p2^2 = 0, p3^2 = 0, p4^2 = 0
c    musq = mu^2 = reference dimensional scale 
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+p3)^2 = 
c          N_ep * D00m_fin(s,t,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

      function D00m_finDiv(s,t,musq,i)
      implicit none
      complex * 16 D00m_finDiv
      real * 8 s,t,musq
      real * 8 ms,mt
      complex * 16 lnms,lnmt
      complex * 16 ipi
      parameter (ipi=(0d0,3.14159265358979323846264338328d0))
      real * 8 prefactor, theta
      external theta
      integer i

      
      if (musq.lt.0d0) then
         write(*,*) 
     #        'POSSIBLE ERROR IN D00m_fin: SCALE MUSQ LESS THAN ZERO!!'
      endif

      prefactor = 1d0/(s*t)

      if(i.eq.1) then
      ms = -s/musq
      mt = -t/musq
            
      lnms = log(abs(ms)) - ipi*theta(-ms)
      lnmt = log(abs(mt)) - ipi*theta(-mt)

      D00m_finDiv =-2d0*prefactor * (lnms+lnmt)
      return
      else
      D00m_finDiv = 4d0*prefactor 
      return
      endif

      end


      function D01m_finDiv(s,t,m4sq,musq,i)
      implicit none
      complex * 16 D01m_finDiv
      real * 8 s,t,m4sq,musq
      real * 8 ms,mt,mm4sq
      complex * 16 lnms,lnmt,lnmm4sq
      complex * 16 ipi
      parameter (ipi=(0d0,3.14159265358979323846264338328d0))
      real * 8 prefactor, theta
      integer i
      
      if (musq.lt.0d0) then
         write(*,*) 
     #        'POSSIBLE ERROR IN D01m_fin: SCALE MUSQ LESS THAN ZERO!!'
      endif

      prefactor = 1d0/(s*t)

      if(i.eq.1) then
      ms = -s/musq
      mt = -t/musq
      mm4sq = - m4sq/musq
      
      lnms = log(abs(ms)) - ipi*theta(-ms)
      lnmt = log(abs(mt)) - ipi*theta(-mt)
      lnmm4sq = log(abs(mm4sq)) - ipi*theta(-mm4sq)

      D01m_finDiv = -2d0*prefactor *(lnms+lnmt-lnmm4sq)
      return
      else

      D01m_finDiv = 2d0*prefactor 
      return
      endif


      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCC                SCALAR INTEGRALS                CCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c    ***********   D02m_fin(s,t,m3sq,m4sq,musq)   *************
c
c Scalar box with MASSLESS PROPAGATORS and with external kinematics: 
c
c                     k
c     p1 ->-----------<-------------<-- p4
c             |                 |
c             |                 |
c             |                 |  
c     p2 ->-------------------------<-- p3
c
c
c    s = (p1+p2)^2
c    t = (p2+p3)^2
c    p1^2 = 0, p2^2 = 0, p3^2 = m3sq <>0, p4^2 = m4sq <>0
c    musq = mu^2 = reference dimensional scale 
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+p3)^2 = 
c          N_ep * D02m_fin(s,t,m3sq,m4sq,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

      function D02m_finDiv(s,t,m3sq,m4sq,musq,i)
      implicit none
      complex * 16 D02m_finDiv
      real * 8 s,t,m3sq,m4sq,musq
      real * 8 ms,mt,mm3sq,mm4sq
      complex * 16 lnms,lnmt,lnmm3sq,lnmm4sq
       complex * 16 ipi
      parameter (ipi=(0d0,3.14159265358979323846264338328d0))
      real * 8 prefactor, theta
      integer i      

      if (musq.lt.0d0) then
         write(*,*) 
     #        'POSSIBLE ERROR IN D02m_fin: SCALE MUSQ LESS THAN ZERO!!'
      endif

      prefactor = 1d0/(s*t)

      if(i.eq.1) then
      ms = -s/musq
      mt = -t/musq
      mm3sq = - m3sq/musq
      mm4sq = - m4sq/musq
      
      lnms = log(abs(ms)) - ipi*theta(-ms)
      lnmt = log(abs(mt)) - ipi*theta(-mt)
      lnmm3sq = log(abs(mm3sq)) - ipi*theta(-mm3sq)
      lnmm4sq = log(abs(mm4sq)) - ipi*theta(-mm4sq)
 
        
      D02m_finDiv = prefactor *(lnmm3sq+lnmm4sq-lnms-2*lnmt)
      return
      else

       D02m_finDiv = prefactor
      return
      endif

      end


c    ***********   D02m_fine(s,t,m1sq,m3sq,musq)   *************
c
c Scalar box with MASSLESS PROPAGATORS and with external kinematics: 
c
c                     k
c     p1 ->-----------<-------------<-- p4
c             |                 |
c             |                 |
c             |                 |  
c     p2 ->-------------------------<-- p3
c
c
c    s = (p1+p2)^2
c    t = (p2+p3)^2
c    p1^2 = m1sq, p2^2 = 0, p3^2 = m3sq , p4^2 =0
c    musq = mu^2 = reference dimensional scale 
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+p3)^2 = 
c          N_ep * D02m_fin3(s,t,m1sq,m3sq,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

      function D02m_fineDiv(s,t,m1sq,m3sq,musq,i)
      implicit none
      complex * 16 D02m_fineDiv
      real * 8 s,t,m1sq,m3sq,musq
      real * 8 ms,mt,mm1sq,mm3sq
      complex * 16 lnms,lnmt,lnmm1sq,lnmm3sq
      complex * 16 ipi
      parameter (ipi=(0d0,3.14159265358979323846264338328d0))
      real * 8 prefactor, theta
      integer i       


      if (musq.lt.0d0) then
         write(*,*) 
     #        'POSSIBLE ERROR IN D02m_fine: SCALE MUSQ LESS THAN ZERO!!'
      endif

      if (i.eq.1) then
      prefactor = 1d0/(m1sq*m3sq- s*t)

      ms = -s/musq
      mt = -t/musq
      mm1sq = - m1sq/musq
      mm3sq = - m3sq/musq
      
      lnms = log(abs(ms)) - ipi*theta(-ms)
      lnmt = log(abs(mt)) - ipi*theta(-mt)
      lnmm1sq = log(abs(mm1sq)) - ipi*theta(-mm1sq)
      lnmm3sq = log(abs(mm3sq)) - ipi*theta(-mm3sq)


      D02m_fineDiv = 2.d0*prefactor*(lnms+lnmt-lnmm1sq-lnmm3sq)
      Return

      else

      D02m_fineDiv =(0d0,0d0)
      endif

      end


      function D03m_finDiv(s,t,m2sq,m3sq,m4sq,musq,i)
      implicit none
      complex * 16 D03m_finDiv
      real * 8 s,t,m2sq,m3sq,m4sq,musq
      real * 8 ms,mt,mm2sq,mm4sq
      complex * 16 lnms,lnmt,lnmm2sq,lnmm4sq
      complex * 16 ipi
      parameter (ipi=(0d0,3.14159265358979323846264338328d0))
      real * 8 prefactor, theta
      integer i

      if (musq.lt.0d0) then
         write(*,*) 
     #        'POSSIBLE ERROR IN D03m_fin: SCALE MUSQ LESS THAN ZERO!!'
      endif


      if(i.eq.1) then

      prefactor = 1d0/(s*t-m2sq*m4sq)
      ms = -s/musq
      mt = -t/musq
      mm2sq = - m2sq/musq

      mm4sq = - m4sq/musq

      lnms = log(abs(ms)) - ipi*theta(-ms)
      lnmt = log(abs(mt)) - ipi*theta(-mt)
      lnmm2sq = log(abs(mm2sq)) - ipi*theta(-mm2sq)
      lnmm4sq = log(abs(mm4sq)) - ipi*theta(-mm4sq)
  
      D03m_finDiv = prefactor *(lnmm2sq+lnmm4sq-lnms-lnmt)
      return
      else

      D03m_finDiv =(0d0,0d0)
      return
      endif
      end


      function D0finDiv(s,t,p1sq,p2sq,p3sq,p4sq,musq,i)
      implicit none
      complex * 16 D0finDiv
      real * 8 s,t,musq,p1sq,p2sq,p3sq,p4sq
      Complex*16 D01m_finDiv,D02m_finDiv,D02m_fineDiv,D03m_finDiv,D00m_finDiv
      External D01m_finDiv,D02m_finDiv,D02m_fineDiv,D03m_finDiv,D00m_finDiv
      real * 8 tiny
      parameter (tiny=1d-6)
      integer i

#ifdef UVPOLE_ONLY
      if(i.eq.1) then
	D0finDiv = (1d0,0d0)
      else
	D0finDiv = 0d0
      endif
      return
#endif

c D0finDiv(s,t,p1sq,p2sq,p3sq,p4sq,musq,i)
c Author: Francisco Campanario
c Date:2008-04-16
c Modified: 2009-01-14(Including s=0 or t=0)
c I have used Carlo Oleari D02m,D03m,D04m
c Bug solve in analytical continuation of D03m
c Need  my  D00m,D02me,D01m functions
c
c                     k
c     p1 ->-----------<-------------<-- p4
c             |                 |
c             |                 |
c             |                 |  
c     p2 ->-------------------------<-- p3
c
c
c    s=(p1+p2)^2
c    t=(p2+p3)^2
c    musq = mu^2 reference dimensional scale
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+q1)^2= 
c          N_ep * D0finDiv(s,t,q1^2,q2^2)
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)
c   the external momenta could be whaterver value      


      if(abs(p1sq).lt.tiny) then
          if(abs(p2sq).lt.tiny)then
               if(abs(p3sq).lt.tiny)then
                   if(abs(p4sq).lt.tiny) then
                      if(abs(s).lt.tiny) then
                         if(abs(t).lt.tiny) then   
c  D(0,0,0,0,0,0) 
                             D0finDiv=0d0
                               RETURN
                         else    
c  D(0,t,0,0,0,0)                         
                            D0finDiv=0d0
                               RETURN
                         endif  !t if
                      else
                         if(abs(t).lt.tiny) then
c  D(s,0,0,0,0,0)
                            D0finDiv=0d0
                               RETURN
                         else      
c  D(s,t,0,0,0,0)                                                      
                             D0finDiv=D00m_finDiv(s,t,musq,i)
                             RETURN
                         endif  !t if
                       endif    !s if   
                   else   ! p4sq
                      if(abs(s).lt.tiny) then
                         if(abs(t).lt.tiny) then   
c  D(0,0,0,0,0,p4sq) 
                             D0finDiv=0d0
                               RETURN
                         else    
c  D(0,t,0,0,0,p4sq)                         
                            D0finDiv=0d0
                               RETURN
                         endif  !t if
                      else
                         if(abs(t).lt.tiny) then
c  D(s,0,0,0,0,p4sq)
                            D0finDiv=0d0
                               RETURN
                         else      
c D(s,t,0,0,0,p4sq) 
                            D0finDiv=D01m_finDiv(s,t,p4sq,musq,i)
                              RETURN

                         endif  !t if
                       endif    !s if   
                   endif        !m4sq if
               else             !m3sq if
                   if(abs(p4sq).lt.tiny) then
                      if(abs(s).lt.tiny) then
                         if(abs(t).lt.tiny) then   
c  D(0,0,0,0,p3sq,0) 
                             D0finDiv=0d0
                               RETURN
                         else    
c  D(0,t,0,0,p3sq,0)                         
                            D0finDiv=0d0
                               RETURN
                         endif  !t if
                      else
                         if(abs(t).lt.tiny) then
c  D(s,0,0,0,p3sq,0)
                            D0finDiv=0d0
                               RETURN
                         else      
c  D(s,t,0,0,p3sq,0)                                                      
                             D0finDiv=D01m_finDiv(s,t,p3sq,musq,i) 
                             RETURN
                         endif  !t if
                       endif     
                   else   ! m4sq
                      if(abs(s).lt.tiny) then
                         if(abs(t).lt.tiny) then   
c  D(0,0,0,0,p3sq,p4sq) 
                             D0finDiv=0d0
                               RETURN
                         else    
c  D(0,t,0,0,p3sq,p4sq)                         
                            D0finDiv=0d0
                               RETURN
                         endif  !t if
                      else
                         if(abs(t).lt.tiny) then
c  D(s,0,0,0,p3sq,p4sq)
                            D0finDiv=0d0
                               RETURN
                         else      
c  D(s,t,0,0,p3sq,p4sq)                                                      
                         D0finDiv=D02m_finDiv(s,t,p3sq,p4sq,musq,i)     
                             RETURN
                         endif  !t if
                       endif    !s if   
                   endif  ! p4sq 
               endif      ! p3sq
          else          ! p2sq
                     if(abs(p3sq).lt.tiny)then
                           if(abs(p4sq).lt.tiny) then
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,0,p2sq,0,0) 
                                   D0finDiv=0d0
                                   RETURN
                                 else    
c  D(0,t,0,p2sq,0,0)                         
                                   D0finDiv=0d0
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,0,p2sq,0,0)
                                   D0finDiv=0d0
                                   RETURN
                                 else      
c  D(s,t,0,p2sq,0,0)                                                      
                                    D0finDiv=D01m_finDiv(s,t,p2sq,musq,i)
                                    RETURN
                                 endif  !t if
                              endif    !s if   
                           else   !p4sq
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,0,p2sq,0,p4sq) 
                                   D0finDiv=D00m_finDiv(p2sq,p4sq,musq,i)
                                   RETURN
                                 else    
c  D(0,t,0,p2sq,0,p4sq)                         
                                   D0finDiv=D01m_finDiv(p2sq,p4sq,t,musq,i)
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,0,p2sq,0,p4sq)
                                   D0finDiv=D01m_finDiv(p2sq,p4sq,s,musq,i)
                                   RETURN
                                 else      
c  D(s,t,0,p2sq,0,p4sq)                                                      
                                     D0finDiv=D02m_fineDiv(s,t,p2sq,p4sq,musq,i) 
                                    RETURN
                                 endif  !t if
                              endif    !s if   
                           endif  ! p4sq
                     else    !p3sq
                            if(abs(p4sq).lt.tiny) then
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,0,p2sq,p3sq,0) 
                                   D0finDiv=0d0
                                   RETURN
                                 else    
c  D(0,t,0,p2sq,p3sq,0)                         
                                   D0finDiv=0d0
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,0,p2sq,p3sq,0)
                                   D0finDiv=0d0
                                   RETURN
                                 else      
c  D(s,t,0,p2sq,p3sq,0)                                                      
                                    D0finDiv=D02m_finDiv(t,s,p2sq,p3sq,musq,i) 
                                    RETURN
                                 endif  !t if
                              endif    !s if   
                             else  ! p4sq
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,0,p2sq,p3sq,p4sq) 
                                   D0finDiv=D01m_finDiv(p2sq,p4sq,p3sq,musq,i)
                                   RETURN
                                 else    
c  D(0,t,0,p2sq,p3sq,p4sq)                         
                                   D0finDiv=D02m_finDiv(p2sq,p4sq,p3sq,t,musq,i)
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,0,p2sq,p3sq,p4sq)
                                   D0finDiv=D02m_finDiv(p4sq,p2sq,s,p3sq,musq,i)
                                   RETURN
                                 else      
c  D(s,t,0,p2sq,p3sq,p4sq)                                                      
                                   D0finDiv=D03m_finDiv(s,t,p2sq,p3sq,p4sq,musq,i) 
                                   RETURN
                                 endif  !t if
                              endif    !s if   
                             endif  ! p4sq
                      endif !p3sq
           endif   !p2sq
      else    ! p1sq
           if(abs(p2sq).lt.tiny)then
               if(abs(p3sq).lt.tiny)then
                   if(abs(p4sq).lt.tiny) then


                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,p1sq,0,0,0) 
                                   D0finDiv=0d0
                                   RETURN
                                 else    
c  D(0,t,0,p1sq,0,0)                         
                                   D0finDiv=0d0
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,0,p1sq,0,0)
                                   D0finDiv=0d0
                                   RETURN
                                 else      
c  D(s,t,0,p1sq,0,0)                                                      

                             D0finDiv=D01m_finDiv(s,t,p1sq,musq,i)
                             RETURN
                                 endif  !t if
                              endif    !s if   
                   else  ! p4sq
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,p1sq,0,0,p4sq) 
                                   D0finDiv=0d0
                                   RETURN
                                 else    
c  D(0,t,p1sq,0,0,p4sq)                         
                                   D0finDiv=0d0
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,p1sq,0,0,p4sq)
                                   D0finDiv=0d0
                                   RETURN
                                 else      
c  D(s,t,p1sq,0,0,p4sq)                                                      
                             D0finDiv=D02m_finDiv(t,s,p4sq,p1sq,musq,i) 
                             RETURN
                                 endif  !t if
                              endif    !s if   
                   endif !p4sq
               else  !p3sq
                   if(abs(p4sq).lt.tiny) then
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,p1sq,0,p3sq,0) 
                                   D0finDiv=D00m_finDiv(p1sq,p3sq,musq,i)
                                   RETURN
                                 else    
c  D(0,t,p1sq,0,p3sq,0)                         
                                   D0finDiv=D01m_finDiv(p1sq,p3sq,t,musq,i)
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,p1sq,0,p3sq,0)
                                   D0finDiv=D01m_finDiv(p1sq,p3sq,s,musq,i)
                                   RETURN
                                 else      
c  D(s,t,p1sq,0,p3sq,0)                                                      
                             D0finDiv=D02m_fineDiv(s,t,p1sq,p3sq,musq,i) 
                             RETURN
                                 endif  !t if
                              endif    !s if   
                   else   !p4sq
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,p1sq,0,p3sq,p4sq) 
                                   D0finDiv=D01m_finDiv(p1sq,p3sq,p4sq,musq,i)
                                   RETURN
                                 else    
c  D(0,t,p1sq,0,p3sq,p4sq)                         
                                   D0finDiv=D02m_finDiv(p1sq,p3sq,t,p4sq,musq,i)
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,p1sq,0,p3sq,p4sq)
                                   D0finDiv=D02m_finDiv(p3sq,p1sq,s,p4sq,musq,i)
                                   RETURN
                                 else      
c  D(s,t,p1sq,0,p3sq,p4sq)                                                      
                             D0finDiv=D03m_finDiv(t,s,p3sq,p4sq,p1sq,musq,i)
                            RETURN
                                 endif  !t if
                              endif    !s if   
                   endif  !p4sq
               endif   !p3sq
          else   !p2sq
               if(abs(p3sq).lt.tiny)then
                   if(abs(p4sq).lt.tiny) then
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,p1sq,p2sq,0,0) 
                                   D0finDiv=0d0
                                   RETURN
                                 else    
c  D(0,t,p1sq,p2sq,0,0)                         
                                   D0finDiv=0d0
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,p1sq,p2sq,0,0)
                                   D0finDiv=0d0
                                   RETURN
                                 else      
c  D(s,t,p1sq,p2sq,0,0) 
                             D0finDiv=D02m_finDiv(s,t,p1sq,p2sq,musq,i)
                             RETURN
                                 endif  !t if
                              endif    !s if   
                   else !p4sq
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,p1sq,p2sq,0,p4sq) 
                                   D0finDiv=D01m_finDiv(p2sq,p4sq,p1sq,musq,i)
                                   RETURN
                                 else    
c  D(0,t,p1sq,p2sq,0,p4sq)                         
                                   D0finDiv=D02m_finDiv(p4sq,p2sq,t,p1sq,musq,i)
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,p1sq,p2sq,0,p4sq)
                                   D0finDiv=D02m_finDiv(p2sq,p4sq,s,p1sq,musq,i)
                                   RETURN
                                 else      
c  D(s,t,p1sq,p2sq,0,p4sq)                                                      
                             D0finDiv=D03m_finDiv(s,t,p4sq,p1sq,p2sq,musq,i) 
                                   RETURN
                                 endif  !t if
                              endif    !s if   
                   endif !p4sq
               else   !p3sq
                   if(abs(p4sq).lt.tiny) then

                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,p1sq,p2sq,p3sq,0) 
                                   D0finDiv=D01m_finDiv(p1sq,p3sq,p2sq,musq,i)
                                   RETURN
                                 else    
c  D(0,t,p1sq,p2sq,p3sq,0)                         
                                   D0finDiv=D02m_finDiv(p3sq,p1sq,p2sq,t,musq,i)
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,p1sq,p2sq,p3sq,0)
                                   D0finDiv=D02m_finDiv(p1sq,p3sq,p2sq,s,musq,i)
                                   RETURN
                                 else      
c  D(s,t,p1sq,p2sq,p3sq,0) 
                             D0finDiv=D03m_finDiv(t,s,p1sq,p2sq,p3sq,musq,i) 
                             RETURN
                                 endif  !t if
                              endif    !s if   
                   else  !p4sq
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,p1sq,p2sq,p3sq,p4sq) 
                                   D0finDiv=D02m_fineDiv(p2sq,p4sq,p1sq,p3sq,musq,i)
                                   RETURN
                                 else    
c  D(0,t,p1sq,p2sq,p3sq,p4sq)                         
                                   D0finDiv=D03m_finDiv(p4sq,p2sq,p3sq,t,p1sq,musq,i)
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,p1sq,p2sq,p3sq,p4sq)
                                   D0finDiv=D03m_finDiv(p4sq,p2sq,p1sq,s,p3sq,musq,i)
                                   RETURN
                                 else      
c  D(s,t,p1sq,p2sq,p3sq,p4sq)                                                      
                             D0finDiv=(0d0,0d0)
                             RETURN
                                 endif  !t if
                              endif    !s if   
                   endif !p4sq
               endif  !p3sq
           endif  ! p2sq
         endif  !p1sq
         end


