      subroutine GFHiggs_B0FinG_QUAD(m,p1sq,musq,B0)
      real*16 m,p1sq,musq
      complex*32 B0
      complex*32 B0tMQUAD
      external B0tMQUAD
      real*16 eps
      parameter(eps=1d-7)
      If(abs(m).gt.eps) then
         B0= B0tMQUAD(m,p1sq)
      else
         print*,"QUAD:B0 No program yet"
         stop
      endif
      end subroutine
      

c
c--------------  B0tM(m,qsq): regularized 2-point function --------------
c
      function B0tMQUAD(m,qsq) 
      use avh_olo_qp
      implicit none
      real*16 m, qsq, qsqn
      complex*32 B0tMQUAD
c evaluate scalar 2-point function for equal masses m on propagators 
c  
c    B0 = Int d^4k [k^2-m^2]^-1 [(k-q)^2-m^2]^-1 
c
c Subtracting the divergent piece, 1/eps - gamma + log(4pi mu^2/m^2),
c one obtains the modified scalar 2-point function B_0~ which is evaluated 
c here
c
c   B0tM(m,q^2) = - int_0^1 dx log[ 1 - q^2/(m^2-i eps) x(1-x) ]
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2000 April 7
c	Last modified:    2000 November 12
c  
c     For OneLOop 
      complex*32 reslt(0:2)
      real*16 mu2,m12,m22    
      real*16 phi, beta, srt, lnfac, re, im
      real*16 eps, pi
      parameter (pi=3.141 592 653 589 793 238 462 643 383 279 502 88q0)
      parameter (eps=5q-8)  ! limit of q^2/m^2 << 1 approximation 
      
      logical UseOneLOop
      common /OneLOopLogical/ UseOneLOop
      
      
    
      if(UseOneLOop) then
            m12=m*m
            m22=m*m
            mu2=1.0q0    
        call olo_scale(dble(sqrt(mu2)))
        call olo_onshell(dble(eps))
        call olo_b0(reslt,qsq,m12,m22)
        
        B0tMQUAD = reslt(0)- log(mu2/m12)*QCMPLX(1.0q0,0.0q0)
        
c        print*,'------'
c        print*, 'B02(oneloop)=',B0tMQUAD
        return
      endif


      qsqn = qsq/(m*m)
      if ( qsqn.lt.-eps ) then
         srt = sqrt(1q0-4q0/qsqn)
         lnfac = log( (srt-1q0)/(srt+1q0) )
         B0tMQUAD = 2q0 + srt*lnfac
      elseif ( abs(qsqn).le.eps )then
         B0tMQUAD = qsqn/6q0* ( 1q0+qsqn*0.1q0*(1q0+qsqn/7q0 *
     &        ( 1q0+qsqn/7q0*( 1q0+2q0/11q0*qsqn ) ) ) )
      elseif (qsqn.lt.4q0) then
         srt = sqrt(4q0/qsqn-1q0)
         phi = atan(1q0/srt)
         B0tMQUAD = 2q0 - 2q0*srt*phi
      elseif (qsqn.eq.4q0) then
         B0tMQUAD = 2q0
      else
         beta = sqrt(1q0-4q0/qsqn)
         lnfac = log( (1q0-beta)/(1q0+beta) )
         re = 2q0 + beta*lnfac
         im = pi*beta
         B0tMQUAD = QCMPLX( re, im )
      endif
c      print*, 'B02(internal)=',B0tMQUAD
      return
      end
