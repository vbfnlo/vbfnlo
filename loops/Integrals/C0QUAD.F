      function C0finG_QUAD(m0,m1,m2,s1,s2,s3,musq)
      IMPLICIT NONE
      complex*32 C0finG_QUAD
      real*16 m0,m1,m2,s1,s2,s3,musq
      complex*32 I3pointQUAD1!,C0finG1M,C0finG1MPrime,C0finG21,C0finG2M
      External I3pointQUAD1! ,C0finG1M,C0finG1MPrime,C0finG21,C0finG2M,I3point
      real*16 eps
c      complex*16 rslt(0:2)
      parameter (eps=1d-7)
      logical UseAVH, UseFF
      real*16 xpi(6)
      integer ier
      parameter (UseAVH=.false.)
      parameter (UseFF=.true.)
cccccccccccccccccccccccccccccccccccccccccccccccc
c Input
c --------------------------------------
c   C0finG1M=C0finG(0d0,0d0,m2,s1,s2,s3,musq)                             
c --------------------------------------
c --------------------------------------
c   C0finG1MPrime=C0finG(0d0,M1,0d0,s1,s2,s3,musq)                        
c --------------------------------------
c
c  C0finG2(M1,M2,M3,s1,s2,s3,musq) only s3=0
c-------------------------------------------
c
c      C0fin(0,0,0,s1,s2,s3)
c----------------------------------------------
c      if(UseAVH) then
cFC#ifdef WITH_ONELOOP
cFC         call avh_olo_mu_set(sqrt(musq))
cFC         call avh_olo_onshell(eps)
cFC         call avh_olo_c0m(rslt,s1,s2,s3,m0**2,m1**2,m2**2)
cFC         C0finG = rslt(0) 
cFC#else
cFC         call avh_error()
cFC#endif
c      else 
        if (UseFF) then
        if ( (abs(m0).lt.eps) .and.
     &       (abs(m1).lt.eps) .and.
     &       (abs(m2).lt.eps) ) then
           call C0fin_QUAD(s1,s2,s3,musq,C0finG_QUAD)
        else if((m0.eq.m1).and.(m0.eq.m2)) then 
            C0finG_QUAD=I3pointQUAD1(m0,s1,s2,s3)
            return
        else
          xpi(1) = m0**2
          xpi(2) = m1**2
          xpi(3) = m2**2
          xpi(4) = s1
          xpi(5) = s2
          xpi(6) = s3
          ier=0
          call vbfquadffxc0(C0finG_QUAD,xpi,ier)
          return
        endif
      else
c                   M0,M1,M3
        if ( (abs(m0).lt.eps) .and.
     &       (abs(m1).lt.eps) .and.
     &       (abs(m2).lt.eps) ) then
           call C0fin_QUAD(s1,s2,s3,musq,C0finG_QUAD)
                return
          endif
          if((m0.eq.m1).and.(m0.eq.m2)) then 
                    C0finG_QUAD=I3pointQUAD1(m0,s1,s2,s3)
                    return
               endif

          print*, "C: Not debug"
          stop
        
      endif

      end

c     C0fin_QUAD(p1sq,p2sq,p3sq,musq)
c
c             
c              /|===== p3^2
c             / | 
c            /  |
c           /   |
c p1^2 =====    | 
c           \   |
c            \  |
c             \ |
c              \|===== p2^2
c             
c
c    musq = mu^2 e' una scala dimensionale esterna!!
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2 = 
c          N_ep * C0fin_QUAD(p1sq,p2sq,p3sq,musq)
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)
c      complex*32 function C0fin_QUAD(p1sq,p2sq,p3sq,musq)
      subroutine C0fin_QUAD(p1sq,p2sq,p3sq,musq,C0_fin_QUAD)
      implicit none
      complex*32 C0_fin_QUAD
      real*16  p1sq,p2sq,p3sq,musq
c     ris(1) = finite part, ris(2) = coeff of 1/ep, 
c     ris(3) = coeff of 1/ep^2      
      complex*32 ris(3)
      real*16  qsq(3),tmp(3)
      real*16  dilog_QUAD
      complex*32 lr,l1,l2,lr2,lr3
      complex*32 I
      parameter (I=(0q0,1q0))
      real*16 pi,pi2
      parameter (pi=3.14159265358979323846264338327950288q0,
     &     pi2 =9.86960440108935861883449099987615114q0)
      complex*32 ipi
      parameter (ipi=(0,3.14159265358979323846264338327950288q0))
      integer j,offshellleg,imax,ii
c      complex*32 C03
      real*16  r3,r2,detsq,det,x,y,max,lomx,lomy
c      complex*16 C0t

      real*16  tiny
      parameter (tiny=1q-15)
      External dilog_QUAD

c      complex*32 C0

c      common/debug_B0C0D0/debugB0C0D0

c      if (debugB0C0D0) then
c         C0fin_QUAD = C0(p1sq,p2sq,p3sq,musq)
c         return
c      endif


      if (musq.lt.0q0) then
      write(*,*) 'ERROR in C0fin_QUAD: mu^2 MUST be a positive number'
      stop
      endif
      
c      c0fin_QUAD = C0t(p1sq,p2sq,p3sq,musq)
c      return
      
c      write(*,*) 'C0 with args',p1sq,p2sq,p3sq,musq

      offshellleg = 0q0
      do j=1,3 
         qsq(j) = 0.q0
      enddo
c      if (p1sq.ne.0) then
      if (abs(p1sq).gt.tiny) then
         offshellleg = offshellleg +1
         qsq(offshellleg) = p1sq
      endif
c      if (p2sq.ne.0) then
      if (abs(p2sq).gt.tiny) then
         offshellleg = offshellleg +1
         qsq(offshellleg) = p2sq
      endif
c      if (p3sq.ne.0) then
      if (abs(p3sq).gt.tiny) then
         offshellleg = offshellleg +1
         qsq(offshellleg) = p3sq
      endif
c      Print*, 'offshellleg',offshellleg
      if (offshellleg.eq.1) then
c     TRIANGLE WITH ONLY ONE EXTERNAL INVARIANT      
         l1 =  Log(abs(qsq(1)/musq))
         if (qsq(1).lt.0q0) then
c     do nothing
         else
            l1 = l1 -ipi
         endif
         ris(3) = 1.q0
         ris(2) = -l1
         ris(1) = 1.q0/2q0*l1*l1 - pi2/6q0         
         do j=1,3
            ris(j) = ris(j)/qsq(1)
         enddo

      elseif (offshellleg.eq.2) then
c TRIANGLE WITH TWO EXTERNAL INVARIANTS  
         lr = Log(abs(qsq(2)/qsq(1)))
         l1 = Log(abs(qsq(1)/musq))
         l2 = Log(abs(qsq(2)/musq))
         if ((qsq(1).lt.0q0).and.(qsq(2).lt.0q0)) then
c     do nothing            
         elseif ((qsq(1).gt.0q0).and.(qsq(2).lt.0q0)) then
            lr = lr + ipi
            l1 = l1 - ipi
         elseif ((qsq(1).lt.0q0).and.(qsq(2).gt.0q0)) then
            lr = lr - ipi
            l2 = l2 - ipi
         elseif ((qsq(1).gt.0q0).and.(qsq(2).gt.0q0)) then
            l1 = l1 - ipi
            l2 = l2 - ipi
         endif
         ris(3) = 0.q0
         ris(2) = lr
         ris(1) = 1.q0/2q0*(l1*l1 - l2*l2)
         do j=1,3
            ris(j) = ris(j)/(qsq(1)-qsq(2))
         enddo
      elseif (offshellleg.eq.3) then         
c TRIANGLE WITH THREE EXTERNAL INVARIANTS           
         ris(3) = 0.q0
         ris(2) = 0.q0

c     ris(1) = C03(qsq(1),qsq(2),qsq(3),0d0,0d0,0d0)
         
         
c order the qsq(i) with absolute max value in the first position
         max = 0q0
         do ii=1,3
            if (abs(qsq(ii)).gt.max) then
               max = abs(qsq(ii))
               imax = ii
            endif
         enddo
         tmp(1) = qsq(imax)

         j = 2
         do ii=1,3
            if (ii.ne.imax) then               
               tmp(j) = qsq(ii)
               j = j+1
            endif
         enddo
         do ii=1,3
            qsq(ii) = tmp(ii)
         enddo

c         write(*,*) 'qsq(i)  ============> ', qsq

         r3 = qsq(3)/qsq(1)
         r2 = qsq(2)/qsq(1)

c         if (sqrt(r2)+sqrt(r3).gt.1d0) then 
c            write(*,*) 'NOT YET IMPLEMENTED'
c            stop
c         endif
         
         detsq = (1q0-r2-r3)*(1q0-r2-r3) - 4q0*r2*r3

c         print*,'detsq',detsq

         if (detsq.lt.0q0) then
c            write(*,*) qsq,detsq
c            write(*,*) r2,r3
c            write(*,*) 'WARNING: this case has NOT yet been implemented',detsq
c            write(*,*) 'RETURN 0 from C0fin_QUAD function'
            C0_fin_QUAD = 0q0
            detsq = 0q0
            RETURN
         endif
         
         det = sqrt(detsq)
         x = 1.q0/2q0/r2*(r3+r2-1q0+det)
         y = 1.q0/2q0/r3*(r3+r2-1q0+det)

c         write(*,*) 'x, y ==> ',x,y

c         write(*,*) 'det ==> ',det

         if ((x.gt.1q0).or.(y.gt.1q0)) then             
            write(*,*) 'ERROR in C0fin_QUAD: x and/or y have values 
     #           bigger than one', x, y
            stop
         endif
        
         lomx = Log(1q0-x)
         lomy = Log(1q0-y)
         lr2 = Log(abs(r2))
         lr3 = Log(abs(r3))
         if (r2.lt.0q0) then
            if (qsq(1).lt.0q0) then
               lr2 = lr2 - ipi
            else
               lr2 = lr2 + ipi
            endif
         endif
         if (r3.lt.0q0) then
            if (qsq(1).lt.0q0) then
               lr3 = lr3 - ipi
            else
               lr3 = lr3 + ipi
            endif
         endif

c       write(*,*) 'x,dilog ==>', x,dilog_QUAD(x)     
         
         ris(1) = 1q0/qsq(1)*(1q0-x)*(1q0-y)/(1q0-x*y)*(2q0*dilog_QUAD(x)
     &   +2q0*dilog_QUAD(y)+
     &        (lomx+lomy)**2+2q0*lr2*lomy+2q0*lr3*lomx+lr3*lr2+pi2/3q0)
 
     
      else
c TRIANGLE WITH ZERO EXTERNAL INVARIANTS
         write(*,*) 'Warning: C0fin_QUAD called with external '//
     #        'invariants equal to zero'
         write(*,*) 'dot prods',p1sq,p2sq,p3sq,musq
         ris(3) = 0.q0
         ris(2) = 0.q0
         ris(1) = 0.q0
      endif

      C0_fin_QUAD = ris(1)
      return
      end


c------------- C(m,q1^2,q2^2,P^2)  3-point fuction ---------------------------
c
c      complex*16 function I3point(m,q1sq,q2sq,Psq)
      complex*32 function I3pointQUAD1(mt,q1sqt,q2sqt,Psqt)
      implicit none
      real*16 m, q1sq, q2sq, psq 
      real*16 mt,q1sqt,q2sqt,psqt,eps1

c evaluate scalar 3-point function for equal masses m on propagators 
c  
c  I3 = 1/(i*pi^2) * Int d^4k [k^2-m^2]^-1 [(k+q1)^2-m^2]^-1 [(k-q2)^2-m^2]^-1
c
c     = i/pi^2 C("t Hooft,Veltman) = -C0(Passarino,Veltman)
c
c in terms of 12 dilogarithms. P = q1+q2 is assumed to be space-like. q1 
c and q2 may be space-, time- or lightlike
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2000 November 12
c	Last modified:    
c  

      real*16 lambda, rtlam, q1q2, q1p, q2p, fmsq
      real*16 rp(3), rm(3), eps, precise, pi

      complex*32 sp(3), sm(3), fmt2, ieps, fac, z, di, quadcdilog
c      external li2
      parameter (eps=1q-36, ieps = eps*(0q0,1q0), precise=1q-9 )
      parameter (pi=3.141 592 653 589 793 238 462 643 383 279 502 88q0, fac=-1q0/2q0)


      eps1=1d-10
      if(abs(mt).le.eps1) then
         m=0q0
      else
         m=mt
      endif
      if(abs(q1sqt).le.eps1) then
         q1sq=0q0
      else
         q1sq=q1sqt
      endif
      if(abs(q2sqt).le.eps1) then
         q2sq=0q0
      else
         q2sq=q2sqt
      endif
      if(abs(Psqt).le.eps1) then
         Psq=0q0
      else
         Psq=Psqt
      endif

c determine dot products of input momenta
      q1q2 = 0.5q0*(psq - q1sq - q2sq)
      q1p  = 0.5q0*(psq + q1sq - q2sq)
      q2p  = 0.5q0*(psq - q1sq + q2sq)
      lambda = q1q2**2-q1sq*q2sq
      
c$$$         print*,'m',m
c$$$         print*,'q1sq',q1sq
c$$$         print*,'q2sq',q2sq
c$$$         print*,'Psq',Psq
c$$$         print*,''
c$$$         print*,'q1q2',q1q2
c$$$         print*,'q1p',q1p
c$$$         print*,'q2p',q2p
c$$$         print*,'lambda ',lambda 
c$$$         print*,'**********************'
c$$$         print*,'lambda/lambda+',lambda/(q1q2**2+q1sq*q2sq)
c$$$c         STOP

      if (lambda.le.0q0) then
         print*,'----------- in I3point ----------'
         print*,'m',m
         print*,'q1sq',q1sq
         print*,'q2sq',q2sq
         print*,'Psq',Psq
         print*,''
         print*,'q1q2',q1q2
c         print*,'q1p',q1p
c         print*,'q2p',q2p
         print*,'lambda ',lambda 
         print*,'lambda/lambda+',lambda/(q1q2**2+q1sq*q2sq)
         print*,'------------------'
         write(*,*) " singular lambda in 3-point function. Reset to 0 "
         I3pointQUAD1 = 0q0
c         STOP
         return
      else
         rtlam = sqrt(lambda)
      endif


c determine factors for call of spence functions
      fmsq = 4q0*m**2
      fmt2 = fmsq - ieps
      di = q1sq*q2sq*psq + fmt2*lambda 
      if ( abs(REAL(di)).gt.eps )then
         di = 1q0/di
      else
         print*," singular point called in I3point "
         stop
      endif

      rp(1) = q2p + rtlam
      rm(1) = q2p - rtlam
      if (abs(rm(1)).lt.abs(rp(1))*precise ) then
         rm(1) = q2sq*psq/rp(1)
      elseif (abs(rp(1)).lt.abs(rm(1))*precise ) then
         rp(1) = q2sq*psq/rm(1)
      endif

      rp(2) = q1q2 + rtlam
      rm(2) = q1q2 - rtlam
      if (abs(rm(2)).lt.abs(rp(2))*precise ) then
         rm(2) = q1sq*q2sq/rp(2)
      elseif (abs(rp(2)).lt.abs(rm(2))*precise ) then
         rp(2) = q1sq*q2sq/rm(2)
      endif

      rp(3) = q1p + rtlam
      rm(3) = q1p - rtlam
      if (abs(rm(3)).lt.abs(rp(3))*precise ) then
         rm(3) = q1sq*psq/rp(3)
      elseif (abs(rp(3)).lt.abs(rm(3))*precise ) then
         rp(3) = q1sq*psq/rm(3)
      endif

      if (q1sq.ne.0q0) then
         z = rtlam*sqrt((q1sq-fmt2)*q1sq) 
         sp(1) = di*( q2p*q1sq + z )
         sm(1) = di*( q2p*q1sq - z )
      else
         sp(1) = 0q0
         sm(1) = 0q0
      endif

      if (psq.ne.0q0) then
         z = rtlam*sqrt((psq-fmt2)*psq) 
         sp(2) = di*( q1q2*psq + z )
         sm(2) = di*( q1q2*psq - z )
      else
         sp(2) = 0q0
         sm(2) = 0q0
      endif

      if (q2sq.ne.0q0) then
         z = rtlam*sqrt((q2sq-fmt2)*q2sq) 
         sp(3) = di*( q1p*q2sq + z )
         sm(3) = di*( q1p*q2sq - z )
      else
         sp(3) = 0q0
         sm(3) = 0q0
      endif

      z = 0
      if (q1sq.lt.0q0 .or. q1sq.gt.fmsq) then
         z = z + quadcdilog(rm(1)*sp(1)) + quadcdilog(rm(1)*sm(1)) - 
     &           quadcdilog(rp(1)*sm(1)) - quadcdilog(rp(1)*sp(1))
      elseif (q1sq.gt.0q0) then
         z = z + 2q0*REAL( quadcdilog(rm(1)*sp(1)) ) - 
     &           2q0*REAL( quadcdilog(rp(1)*sp(1)) )
      endif
      if (psq.lt.0q0 .or. psq.gt.fmsq) then
         z = z - quadcdilog(rm(2)*sp(2)) - quadcdilog(rm(2)*sm(2)) + 
     &           quadcdilog(rp(2)*sm(2)) + quadcdilog(rp(2)*sp(2))
      elseif (psq.gt.0q0) then
         z = z - 2q0*REAL( quadcdilog(rm(2)*sp(2)) ) + 
     &           2q0*REAL( quadcdilog(rp(2)*sp(2)) )
      endif
      if (q2sq.lt.0q0 .or. q2sq.gt.fmsq) then
         z = z + quadcdilog(rm(3)*sp(3)) + quadcdilog(rm(3)*sm(3)) - 
     &           quadcdilog(rp(3)*sm(3)) - quadcdilog(rp(3)*sp(3))
      elseif (q2sq.gt.0q0) then
         z = z + 2q0*REAL( quadcdilog(rm(3)*sp(3)) ) - 
     &           2q0*REAL( quadcdilog(rp(3)*sp(3)) )
      endif

      I3pointQUAD1 = z*fac/rtlam

!      Print*, 'I3point',I3point

      return
       end


c
c-------------  spence function  li2(z) ------------------------
c
      complex*32 function QUADcdilog(zin)
      implicit none
      complex*32 zin, z, ans, zext,test
      real*16 r2, fac
c
c determine the value of the dilogarithm 
c
c    li2(z) = - int_0^1  log(1-zt)/t dt  with cut along the positive 
c                                        real axis, z>1
c
c       Francisco Campanario
c  

      real*16 c0,c1,c2,c4,c6,c8,c10,c12,c14,c16,c18,c20,c22
      real*16  c24,c26,c28,c30,c32,c34,c36,c38
      real*16 b0,b1,b2,b4,b6,b8,b10,b12,b14,b16,b18,b20,b22
      real*16  b24,b26,b28,b30,b32,b34,b36,b38
      real*16 d0,d1,d2,d4,d6,d8,d10,d12,d14,d16,d18,d20,d22
      real*16  d24,d26,d28,d30,d32,d34,d36,d38
     
      real*16  I2,I4,I9,I16,I25,I36,I49,I64,I81,I100,I121,I144
      real*16 I169
      complex*32  u,u2,u4,u8,u16,u32
     

      parameter (b0=1q0,            d0 =1q0,      c0= b0/d0,I2=1q0/2q0, I4=1/4q0 )
      parameter (b1=-1q0/2q0,       d1 =2q0,      c1= b1/d1,I9=1q0/9q0)
      parameter (b2= 1q0/6q0,       d2 =6q0,      c2= b2/d2,I16=1q0/16q0)
      parameter (b4=-1q0/30q0,      d4 =12q1,     c4= b4/d4,I25=1q0/25q0)
      parameter (b6=1q0/42q0,       d6 =504q1,    c6= b6/d6,I36=1q0/36q0)
      parameter (b8=-1q0/30q0,      d8 =36288q1,  c8= b8/d8,I49=1q0/49q0)
      parameter (b10=5q0/66q0,      d10=399168q2,  c10=b10/d10,I64=1q0/64q0)
      parameter (b12=-691q0/2730q0, d12=62270208q2, c12=b12/d12,I81=1q0/81q0)
      parameter (b14=7q0/6q0,       d14=1307674368q3, c14=b14/d14,I100=1q0/100q0)
      parameter (b16=-3617q0/510q0, d16=355687428096q3, c16=b16/d16,I121=1q0/121q0)
      parameter (b18=43867q0/798q0, d18=121645100408832q3, c18=b18/d18,I144=1q0/144q0)
      parameter (b20=-174611q0/330q0,d20=5109094217170944q4, c20=b20/d20)
      parameter (b22=854513q0/138q0,d22=25 852 016 738 884 976 64q4,  c22=b22/d22)

      parameter (b24=-(236364091q0/2730q0),d24=15 511 210 043 330 985 984q6,  c24=b24/d24)
      parameter (b26=8553103q0/6q0,d26=10 888 869 450 418 352 160 768q6, c26=b26/d26)
      parameter (b28=-23749461029q0/870q0
     - ,d28=8 841 761 993 739 701 954 543 616q6, c28=b28/d28)
       parameter (b30=8615841276005q0/14322q0,
     - d30=8 222 838 654 177 922 817 725 562 88q7,  c30=b30/d30)
      parameter (b32=-7709321041217q0/510q0,
     -   d32=8 683 317 618 811 886 495 518 194 401 28q7,  c32=b32/d32)
      parameter (b34=2577687858367q0/6q0,
     -  d34=10 333 147 966 386 144 929 666 651 337 523 2q8, c34=b34/d34)
      parameter (b36=-26315271553053477373q0/1919190q0
     - ,d36=13 763 753 091 226 345 046 315 979 581 580 902 4q8,  
     -  c36=b36/d36)
      parameter (b38=2929993913841559q0/6q0,
     - d38=20 397 882 081 197 443 358 640 281 739 902 897 356 8q8,  
     - c38=b38/d38)
      parameter (I169=1q0/169q0)

      real*16 eps,eps1, epst, pi, pi2o6,pi2o3
      parameter (eps=1q-8, epst=1q-6,eps1=1q-70)
      parameter (pi=3.141 592 653 589 793 238 462 643 383 279 502 88q0
     &  , pi2o6=1.644 934 066 848 226 436 472 415 166 646 025 19q0, 
     &    pi2o3=3.289 868 133 696 452 872 944 830 333 292 050 38q0)

      real*16 zreal,zimag
      complex*32 zlog,zlog1
c
c debug information
      logical ldebug
      parameter (ldebug=.false.)

      z = zin
c      print*,' li2 call with z = ',z
      u = z*z
      zreal=REAL(z)
      zimag=IMAG(z)

      r2 = zreal*zreal + zimag*zimag

      if (r2.lT.eps) then
         QUADcdilog = ((((I25*z+I16)*z+I9)*z+I4)*z+1q0)*z 
         return
      elseif (r2.lE.epst) then

         ans=((((((((((((I169*z+I144)*z+I121)*z+I100)*z+I81)*z+I64)*z+I49)*z+I36)
     c     *z+I25)*z+I16)*z+I9)*z+I4)*z+1q0)*z 
         QUADcdilog=ans
         return
      endif
      if (zreal.ge.1q0 .and. zimag.eq.0q0 ) then
         z = z + (0q0,1q0)*eps1
      endif
c
c use z-->1/z and z--> 1-z mappings of the spence function to restrict 
c agument to unit circle in the complex plane with Re(z) <= 0.5
c
      zext = (0q0,0q0)
      fac = 1q0

      if(zREAL.lt.0q0) then  
       ZLOG=Log(-z)
       zlog1=Log(1Q0-z)
       zext = zlog1*(zlog1*I2-zlog) -pi2o6
       z = 1q0/(1q0-z)
       zreal=REAL(z)
       zimag=IMAG(z)
      r2 = zreal*zreal + zimag*zimag
      endif


        if (sqrt(r2).gt.1q0) then     ! map z ---> 1/z
         fac = -fac
         ZLOG=Log(-z)
         zext = zext  -(pi2o6 + I2*zlog*zlog)
         z = 1q0/z
        zreal=REAL(z)
        zimag=IMAG(z)
      endif


      if (zreal.gt.I2) then     ! map new z ---> 1-z
         zlog=log(z)
         zlog1=log(1q0-z)
         zext = zext + fac*(pi2o6-zlog*zlog1)
         fac = -fac
         z = 1q0-z
      endif
c
c now use t = 1 - exp(-u) mapping to write Li(z) in terms of Bernoulli 
c numbers
c
      u = - log(1q0-z)
      u2 = u*u
      u4=u2*u2
      u8=u4*u4
      u16=u8*u8
      u32=u16*u16

      ans = ((c2*u+c1)*u+c0)*u
      ans = (((((c14*u2+c12)*u2+c10)*u2+c8)*u2+c6)*u2+c4)*u4*u+ans
      test=c16*u16*u
      ans=((((((c30*u2+c28)*u2+c26)*u2+c24)*u2+c22)*u2+c20)*u2+c18)*u2*u*u16+
     c test+ans
      test=c32*u32*u
      ans=((c38*u2+c36)*u2+c34)*u2*u32*u+test+ans
      QUADcdilog = fac*ans + zext
      Return
      end



