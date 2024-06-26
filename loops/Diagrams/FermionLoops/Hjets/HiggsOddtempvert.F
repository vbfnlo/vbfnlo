       subroutine gluonverHiggsOdd(m,p1,p2,p3,mup1,mup2,comp,result)
c ************************************************************************************
c determine the pentagon of
c p1(gluon)    ---->------->--------->--------->---   p4(B) p4^2=M_h
c                $                   $              
c                $                   $              
c                $                   $              
c                V                   V              
c                $                   $              
c                $                   $              
c               p2                   p3             
c              mu_p2                 mu_p3           
c              p2^2=0                p3^2=0           
c Note: No couplings are assigned
c Note: To make it shorter in the promgram: mu_p2,...->mup2,... 
c Notation of External momenta: p1+p2+p3+p4=0 
c mu_p1,mu_p2,mu_p3,mu_p4 should be think as external current 
c musq is the renormalization scale energy  
c************************************************************************************
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       IMPLICIT NONE
      Real*8   p1(0:3),p2(0:3),p3(0:3)
       Real*8   p1sq, p1p2, p1p3, p2sq, p2p3, p3sq
       Real*8   s12, s23
       Complex*16   p1mup1, p1mup2, p2mup1, p2mup2,p3mup1, p3mup2
       Complex*16   mup1mup2
       Real*8 dotrr
       Complex*16 C0fin1M
       EXTERNAL dotrr, C0fin1M
       Complex*16 C0123
       Real*8 C0123R
       Real*8 C0123I
      Complex*16 F(1)
      Real*8 FI(1),FR(1)
      Complex*16 mup1(0:3),mup2(0:3)
      Complex*16 dotrc,dotcc,EPSCCRR,result
       Real*8 p5t,m,msq
      EXTERNAL   dotrc,dotcc,EPSCCRR
      Integer alpha
      Save F
      integer comp
c*******************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the scalar products. Not inlcueded the contraction of the
c       moments with the external currents  
c************************************************************************************
c************************************************************************************
c************************************************************************************

       msq=m*m
   

       p1sq = dotrr(p1,p1)
       p1p2 = dotrr(p1,p2)
       p2sq = dotrr(p2,p2)
       s12 = (p1sq +p2sq+ 2*p1p2) 

c************************************************************************************
c************************************************************************************
c************************************************************************************
c    Calling C_ij,D_ij Functions    
c************************************************************************************
c************************************************************************************
       if (comp.gt.0) then 

      C0123=C0fin1M(m,p1sq,p2sq,s12) 
      C0123R=Dble(C0123)
      C0123I=DIMAG(c0123)

c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************


       FR(1) =C0123R
       FI(1) =C0123I
       F(1)=DCMPLX(FR(1),FI(1))
       endif
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c               PART THAT DEPENDS ON THE EXTERNAL CURRENT
c************************************************************************************

c************************************************************************************
c       Amplitude                         
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************


       result =-4*msq*F(1)*EPSCCRR(mup1,mup2,p1,p2)
       result= result/m
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
