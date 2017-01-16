c.....Bubble correction: finite part.           
c.....musq = mu^2 is an external scale
c 
c.....int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2 = 
c          N_ep * B0fin(p1sq,musq)
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)
c
c.....ris(1) = finite part
c.....ris(2) = coeff of 1/ep
c.....ris(3) = coeff of 1/ep^2      

      function B0fin(p1sq,musq)
      implicit none
      complex * 16 B0fin
      real * 8 p1sq,musq
      complex * 16 ris(3)
      complex * 16 l1
      complex * 16 I
      parameter (I=(0,1))
      double precision pi
      parameter (pi=3.141592653589793238462643383279502884197D0)
      complex * 16 ipi
      parameter (ipi=(0,3.141592653589793238462643383279502884197D0))


      real * 8 tiny
      parameter (tiny=1d-7)

      if (musq.lt.0) then
         write(*,*) "ERROR in B0fin: mu^2 MUST be a positive number"
         stop
      endif

      if (abs(p1sq).lt.tiny) then
c     BUBBLE WITH ZERO EXTERNAL INVARIANTS
c      write(*,*) "Warning: B0fin called with external"//
c     #       " invariant equal to zero"
         ris(3) = 0.d0
         ris(2) = 0.d0
         ris(1) = 0.d0
         B0fin = ris(1)
         return
      endif
      l1 =  log(abs(p1sq/musq))
      if (p1sq.gt.0) then
         l1 = l1 -ipi
      endif
      ris(3) = 0d0
      ris(2) = 1d0
      ris(1) = 2-l1
      B0fin = ris(1)
      end


c --------------------------------------
c   B0finG3=B0finG(M1,M1,s,musq)                                 
c --------------------------------------
      complex*16 function B0t1(M1,s,musq)
      implicit none
      real*8 pi
      complex*16 Ipi, Ieps 
      parameter(pi=3.14159265358979324d0,
     & Ipi=(0d0,3.14159265358979324d0),Ieps=(0d0,1d-38))
      double precision M1,s,musq,msq
      complex*16 b,g1,g2

      msq = M1*M1
      if(abs(s).gt.1d-6) then
      b = Sqrt(s*s-4.d0*s*(msq-Ieps))
      g1 = .5d0*(s+b)/s
      g2 = .5d0*(s-b)/s

      B0t1 = -Log((s-Ieps)/musq) 
     &        + (g1*Log((g1-1.d0)/g1)-Log(g1-1))
     &        + (g2*Log((g2-1.d0)/g2)-Log(g2-1)) + 2.d0 


      else
      B0t1  = -dLog(msq/musq)
      
      endif


      end


c --------------------------------------
c   B0finG3=B0finG(M1,M1,s,musq)                                 
c --------------------------------------
      complex*16 function B0tMDiv(M1,s,i)
       double precision M1,s

      integer i
      if(i.eq.1) then
      B0tMDiv = 1.d0
      else
      B0tMDiv = 0d0
      endif
      end


c--------------  B0tM(m,qsq): regularized 2-point function --------------
c
      complex*16 function B0tM1(m,qsq) 
      implicit none
      double precision m, qsq, qsqn

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

      double precision phi, beta, srt, lnfac, re, im
      double precision eps, pi
      parameter (pi=3.14159 26535 89793d0)
      parameter (eps=5d-4)  ! limit of q^2/m^2 << 1 approximation 

      qsqn = qsq/m**2
      if ( qsqn.lt.-eps ) then
         srt = sqrt(1d0-4d0/qsqn)
         lnfac = log( (srt-1d0)/(srt+1d0) )
         B0tM1 = 2d0 + srt*lnfac
      elseif ( abs(qsqn).le.eps )then
         B0tM1 = qsqn/6d0* ( 1d0+qsqn*0.1d0*(1d0+qsqn/7d0 *
     &        ( 1d0+qsqn/7d0*( 1d0+2d0/11d0*qsqn ) ) ) )
      elseif (qsqn.lt.4d0) then
         srt = sqrt(4d0/qsqn-1d0)
         phi = atan(1d0/srt)
         B0tM1 = 2d0 - 2d0*srt*phi
      elseif (qsqn.eq.4d0) then
         B0tM1 = 2d0
      else
         beta = sqrt(1d0-4d0/qsqn)
         lnfac = log( (1d0-beta)/(1d0+beta) )
         re = 2d0 + beta*lnfac
         im = pi*beta
         B0tM1 = cmplx( re, im )
      endif
      return
      end


c --------------------------------------
c   B0finG1=B0finG(0,0,s,musq)                                 
c --------------------------------------
      complex*16 function B0finG1(s,musq)
      implicit none
      real*8 pi
      complex*16 Ipi, Ieps 
      parameter(pi=3.14159265358979324d0,
     & Ipi=(0d0,3.14159265358979324d0),Ieps=(0d0,1d-38))
      double precision s,musq

      if (s.gt.0.d0) then 
         B0finG1 = -(dLog(s/musq)-Ipi)+2.d0
      else
         B0finG1 = -dLog(-s/musq)+2.d0
      endif
      end

c --------------------------------------
c   B0finG2=B0finG(M1,0,s,musq)                                 
c --------------------------------------
      complex*16 function B0finG2(M1,s,musq)
      implicit none
      real*8 pi
      complex*16 Ipi, Ieps 
      parameter(pi=3.14159265358979324d0,
     & Ipi=(0d0,3.14159265358979324d0),Ieps=(0d0,1d-38))
      double precision M1,s,musq,msq

      msq = M1*M1
      if (abs(s).lt.1d-7) then
        B0finG2 = -dLog(msq/musq)+ 1.d0
      elseif (s.gt.msq) then 
         B0finG2 = -dLog(msq/musq)
     &          +(msq-s)/s*(dLog((s-msq)/msq)-Ipi)+2.d0
      else
         B0finG2 = -dLog(msq/musq)
     &          +(msq-s)/s*dLog(-(s-msq)/msq)+2.d0
      endif 
      end
c
