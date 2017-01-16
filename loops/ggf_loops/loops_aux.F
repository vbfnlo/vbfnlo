c
c---------------  tens3 = 3point tensors ------------------------------
c
      subroutine tens3( mt, p1sq, p2sq, psq, B0, C0, 
     &                  Cij)
      implicit none
      double precision mt, p1sq, p2sq, psq
      complex*16 B0(3), C0, Cij(2,4)
c
c  determine the Passarino-Veltman tensor decomposition for the 
c  three-point tensor integrals
c
c                   1                 1;  k_mu;   k_mu k_nu
c C0;C_mu;C_munu =------Int d^4k --------------------------------------
c                 i pi^2        [-k^2+mt^2][-(k+p1)^2+mt^2][-(k+p1+p2)^2+mt^2] 
c
c with
c
c   C_mu = p1_mu C11  +  p_2_mu C12
c
c  C_munu = p1_mu p1_nu C21 + p2_mu p2_nu C22 + 
c           (p1_mu p2_nu + p1_nu p2_mu) C23  -  g_munu C24
c
c  for notation see Passarino&Veltman, NP B160 (1979) 151 and my notes
c  on gg-->Hgg, pages 4.1--4.8
c
C INPUT:  mt, p1sq, p2sq, psq=(p1+p2)^2     mass and 4-momenta**2
C         B0, C0                   4 scalar integrals; the 3 B0 are, 
c                                  in PV notation:
c         B0(1) = B0(1,2)          B_0 function with subtraction of 
c         B0(2) = B0(2,3)          divergent term
c         B0(3) = B0(1,3)
c
c OUTPUT: Cij(n,m) = C_nm          form factors in the tensor integrals
c         n=1,2; m=1,2,3,4         a la PV
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2000 November
c	Last modified:    

      real*8 f1, f2, Xi(2,2), p1p2, deti
      complex*16 b13, b12, b23, r1,r2,r3,r4,r5,r6
      complex*16 C11, C12, C21, C22, C23, C24, C23p
      integer i,j, iscount
      logical ldebug
      parameter (ldebug=.false.)
      data iscount /0/
      save iscount

      real*8 srt
      external srt

      p1p2 = (psq-p1sq-p2sq)/2d0
      f1 = p1sq
      f2 = psq - p1sq
      if (psq.eq.0d0) then
         deti = 0.25d0*(p1sq-p2sq)**2
      elseif (p1sq.eq.0d0) then
         deti = 0.25d0*(psq-p2sq)**2
      elseif (p2sq.eq.0d0) then
         deti = 0.25d0*(psq-p1sq)**2
      else
         deti = p1p2**2-p1sq*p2sq
         if (abs(deti).lt.1d-30) then
            do i=1,2
               do j = 1,4
                  Cij(i,j) = 0
               enddo
            enddo
            iscount = iscount+1
            if (iscount.le.1000) then
               print*," Warning: singular point in TENS3 "
               return
            else
               print*," stop: too many singular points in TENS3 "
               stop
            endif
         endif
c         if (abs(deti).lt.1d-6*max(abs(p1p2**2),abs(p1sq*p2sq))) then
c            print*," singular Cij. deti suppressed by ",
c     1           abs(deti)/max(abs(p1p2**2),abs(p1sq*p2sq))
c            print*," p1 =  ",srt(p1sq)
c            print*," p2 =  ",srt(p2sq)
c            print*," p12 = ",srt(psq)
c         endif         
      endif
      deti = 1d0/deti
      Xi(1,1) = p2sq*deti
      Xi(1,2) = -p1p2*deti
      Xi(2,1) = -p1p2*deti
      Xi(2,2) = p1sq*deti
      b12 = B0(1)
      b23 = B0(2)
      b13 = B0(3)
c
      R1 = (f1*C0+b13-b23)/2d0
      R2 = (f2*C0+b12-b13)/2d0
      C11 = Xi(1,1)*R1 + Xi(1,2)*R2
      C12 = Xi(2,1)*R1 + Xi(2,2)*R2
      Cij(1,1) = C11
      Cij(1,2) = C12
      Cij(1,3) = 0
      Cij(1,4) = 0
c
      C24 = -0.5d0*mt**2*C0 + 0.25d0*(1d0+b23-f1*C11-f2*C12)
      R3 = -b13/4d0 + (b23+f1*C11)/2d0 - C24
      R4 = -(b13-b23)/4d0 + f1*C12/2d0
      R5 = -(b12-b13)/4d0 + f2*C11/2d0
      R6 = b13/4d0 + f2*C12/2d0 - C24

      Cij(2,1) = Xi(1,1)*R3 + Xi(1,2)*R5
      Cij(2,3) = Xi(2,1)*R3 + Xi(2,2)*R5
      Cij(2,2) = Xi(2,1)*R4 + Xi(2,2)*R6
      Cij(2,4) = C24

      if (ldebug) then
         c21 = cij(2,1)
         c22 = cij(2,2)
         c23 = cij(2,3)
         C23p= Xi(1,1)*R4 + Xi(1,2)*R6
         print*," C23 comparison ",C23, C23p,C23/(C23p+1d-30)
         print*,"r1:",p1sq*C11+p1p2*c12+0.5d0*(p1sq*c0+b13-b23)
         print*,"r2:",p1p2*C11+p2sq*C12+0.5d0*((p2sq+2*p1p2)*c0+
     &                b12-b13)
c         print*," c24 = ",c24
c         r1 = p1sq*C21+p2sq*C22+2*p1p2*C23-4*c24+b23-mt**2*C0
c         print*,r1,r1/(c24+1d-60)
         print*,"r4f",p1p2*c22+p1sq*c23+(b13-b23)/4d0+p1sq*c12/2d0
         print*,"r4:",p1p2*c22+p1sq*c23-(b13-b23)/4d0+p1sq*c12/2d0
         print*,"r5f",p1p2*c21+p2sq*c23+(b12-b13)/4d0+(2*p1p2+p2sq)*c11/2d0
         print*,"r5:",p1p2*c21+p2sq*c23-(b12-b13)/4d0+(2*p1p2+p2sq)*c11/2d0
      endif

      return
      end

c
c---------------  tens4 = 4point tensors -----------------------------
c
      subroutine tens4( mt, p1sq, p2sq, p3sq, p1p2, p1p3, p2p3, 
     &                  C0123,  C0124,  C0134,  C0234,
     &                  Cij123, Cij124, Cij134, Cij234, 
     &                  D0, Dij)
      implicit none
      double precision mt, msq, p1sq, p2sq, p3sq, p1p2, p1p3, p2p3
      complex*16 C0123,  C0124,  C0134,  C0234
      complex*16 Cij123(2,4), Cij124(2,4), Cij134(2,4), Cij234(2,4)
      complex*16 D0, Dij(3,13)
c
c  determine the Passarino-Veltman tensor decomposition for the 
c  four-point tensor integrals
c
c                                       1
c   D0; D_mu; D_mu,nu; D_mu,nu,rho = ------- Int d^4k
c                                     i pi^2 
c
c           1;  k_mu;   k_mu k_nu; k_mu k_nu k_rho
c   -------------------------------------------------------------------
c   [-k^2+mt^2][-(k+p1)^2+mt^2][-(k+p1+p2)^2+mt^2][-(k+p1+p2+p3)^2+m^2]
c
c  with
c
c   D_mu = p1_mu D11  +  p2_mu D12  +  p3_mu D13
c
c   D_munu = p1_mu p1_nu D21 + p2_mu p2_nu D22 + ...
c
c  for notation see Passarino&Veltman, NP B160 (1979) 151 and my notes
c  on gg-->Hgg, pages 4.9--4.19. Note: since I use Bjorken Drell g_munu
c  one needs to replace pi.pj --> -pi.pj, delta_munu = -g_munu in PV
c  to get the results below
c
C INPUT:  mt, p1sq, p1p2,...                mass and 4-momenta**2
C         C0123 = C0(1,2,3) = C0(p1,p2)     scalar three point 
C         C0124 = C0(1,2,4) = C0(p1,p2+p3)  functions in PV notation
C         C0134 = C0(1,3,4) = C0(p1+p2,p3)
C         C0234 = C0(2,3,4) = C0(p2,p3)
C         Cij123(n,m) = C_nm(1,2,3) ....    higher C_nm form factors
C                                           as in tens3
c         D0 = D0(p1,p2,p3)                 scalar four point function
c
c OUTPUT: Dij(n,m) = D_nm                   form factors in the tensor 
c                                           integrals a la PV
c         nm = 11, 12, 13                   ff"s for D_mu
c         nm = 21, 22, 23, 24, 25, 26, 27   ff"s for D_munu
c         nm = 31, 32, 33, ..., 39, 310, 311, 312  ff"s for D_mu,nu,rho
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2000 December 
c	Last modified:    2002 October 29
c	Last modified: 19.09.07, Michael Kubocz   
c       Additional D00ij and D0000 coefficients for Denner-Dittmaier-tensor 
c       decomposition (used in SUBROUTINE DDtens51m in ggf_amp_aux.F) 

      real*8 f1, f2, f3, Xi(3,3), X(3,3), deti
      complex*16 r(20:58), Dijp(3,13), t(9), d310pp
      integer i,j, iscount
      logical ldebug, ldb0
      parameter (ldebug=.false.)
      common /bdebug/ ldb0
      real*8 srt !, dtexm
      external srt
      data iscount /0/
      save iscount

      f1 = p1sq
      f2 = p2sq+2*p1p2
      f3 = p3sq+2*(p1p3+p2p3)
      deti = -p1sq*p2sq*p3sq-2*p1p2*p2p3*p1p3
     &       +p1sq*p2p3**2+p2sq*p1p3**2+p3sq*p1p2**2
      do i=1,3
         do j = 1,13
            Dij(i,j) = 0
         enddo
      enddo
      if (abs(deti).lt.1d-30) then
         iscount = iscount+1
         if (iscount.le.1000) then
            print*," Warning: singular point in TENS4 "
            return
         else
            print*," stop: too many singular points in TENS4 "
            stop
         endif
      endif
      deti = 1/deti
c
      Xi(1,1) = (p2sq*p3sq-p2p3**2)*deti
      Xi(1,2) = (p1p3*p2p3-p3sq*p1p2)*deti
      Xi(1,3) = (p1p2*p2p3-p2sq*p1p3)*deti

      Xi(2,1) = Xi(1,2)
      Xi(2,2) = (p1sq*p3sq-p1p3**2)*deti
      Xi(2,3) = (p1p2*p1p3-p1sq*p2p3)*deti

      Xi(3,1) = Xi(1,3)
      Xi(3,2) = Xi(2,3)
      Xi(3,3) = (p1sq*p2sq-p1p2**2)*deti

      if (ldb0) then
         print("(g15.4)")," det(-X) = ",1/deti
         print*," 1,2,3 ",srt(p1sq),srt(p2sq),srt(p3sq)
c         if (p1sq.ne.0) then
c            dtexm = (p1p2+p2p3)**2*(p1p3+p2p3)**2/
c     1              (p1sq+p2sq+p3sq+2*(p1p2+p1p3+p2p3))
c         elseif (p2sq.ne.0) then
c            dtexm = (p1p2+p1p3)**2*(p1p3+p2p3)**2/
c     1              (p1sq+p2sq+p3sq+2*(p1p2+p1p3+p2p3))
c         else
c            dtexm = (p1p2+p1p3)**2*(p1p2+p2p3)**2/
c     1              (p1sq+p2sq+p3sq+2*(p1p2+p1p3+p2p3))
c         endif
c         print*," detx ratio ",1/deti/dtexm
         print*," 12 ",srt(p1sq+2*p1p2+p2sq)," 13 ",
     1          srt(p1sq+2*p1p3+p3sq)," 23 ",srt(p2sq+2*p2p3+p3sq)
         print*," check inversion of X matrix for n = 3 "
         X(1,1) = -p1sq
         X(1,2) = -p1p2
         X(1,3) = -p1p3
         X(2,1) = -p1p2
         X(2,2) = -p2sq
         X(2,3) = -p2p3
         X(3,1) = -p1p3
         X(3,2) = -p2p3
         X(3,3) = -p3sq
         do i = 1,3
            do j = 1,3
               print*,"X*X^-1(",i,j,") = ", 
     &                X(i,1)*Xi(1,j)+X(i,2)*Xi(2,j)+X(i,3)*Xi(3,j)
            enddo
         enddo
         print*," Xi = "
         do i = 1,3
            print 13,(Xi(i,j),j=1,3)
         enddo
         print*," X = "
         do i = 1,3
            print 13,(X(i,j),j=1,3)
         enddo
 13      format (3g14.4)
      endif
c
      if (ldb0) then
         if (.false.) then
            C0123 = 0
            c0124 = 0
            c0234 = 0
            c0134 = 0
            do i=1,2
               do j = 1,4
                  cij123(i,j) = 0
                  cij124(i,j) = 0
                  cij234(i,j) = 0
                  cij134(i,j) = 0
               enddo
            enddo
         endif
         print*," input for D1j "
         print*," D0 = ",d0
         print*," C0123 = ",C0123
         print*," C0124 = ",C0124
         print*," C0234 = ",C0234
         print*," C0134 = ",C0134
      endif
      R(20) = (f1*D0+C0134-C0234)/2d0
      R(21) = (f2*D0+C0124-C0134)/2d0
      R(22) = (f3*D0+C0123-C0124)/2d0
      Dij(1,1) = Xi(1,1)*R(20) + Xi(1,2)*R(21) + Xi(1,3)*R(22)
      Dij(1,2) = Xi(2,1)*R(20) + Xi(2,2)*R(21) + Xi(2,3)*R(22)
      Dij(1,3) = Xi(3,1)*R(20) + Xi(3,2)*R(21) + Xi(3,3)*R(22)
c      do j = 4,13
c         Dij(1,j) = 0
c      enddo
ccccccc
      if (ldebug) then
         print*," r20 = ",R(20)
         print*," r21 = ",R(21)
         print*," r22 = ",R(22)
         print*,"r20c ",
     &           -(p1sq*dij(1,1)+p1p2*dij(1,2)+p1p3*dij(1,3))/R(20)
         print*,"r21c ",
     &           -(p1p2*dij(1,1)+p2sq*dij(1,2)+p2p3*dij(1,3))/R(21)
         print*,"r22c ",
     &           -(p1p3*dij(1,1)+p2p3*dij(1,2)+p3sq*dij(1,3))/R(22)
      endif
ccccccc
      Dij(2,7) = -mt**2*D0 
     &           - (f1*dij(1,1)+f2*dij(1,2)+f3*dij(1,3)-C0234)/2d0
c
      R(30) = (f1*Dij(1,1)+Cij134(1,1)+C0234)/2d0 - Dij(2,7)
      R(33) = (f1*Dij(1,2)+Cij134(1,1)-Cij234(1,1))/2d0
      R(36) = (f1*Dij(1,3)+Cij134(1,2)-Cij234(1,2))/2d0

      R(31) = (f2*Dij(1,1)+Cij124(1,1)-Cij134(1,1))/2d0
      R(34) = (f2*Dij(1,2)+Cij124(1,2)-Cij134(1,1))/2d0 - Dij(2,7)
      R(37) = (f2*Dij(1,3)+Cij124(1,2)-Cij134(1,2))/2d0

      R(32) = (f3*Dij(1,1)+Cij123(1,1)-Cij124(1,1))/2
      R(35) = (f3*Dij(1,2)+Cij123(1,2)-Cij124(1,2))/2
      R(38) = (f3*Dij(1,3)-Cij124(1,2))/2d0 - Dij(2,7)

      Dij(2,1) = Xi(1,1)*R(30) + Xi(1,2)*R(31) + Xi(1,3)*R(32)
      Dij(2,4) = Xi(2,1)*R(30) + Xi(2,2)*R(31) + Xi(2,3)*R(32)
      Dij(2,5) = Xi(3,1)*R(30) + Xi(3,2)*R(31) + Xi(3,3)*R(32)

      Dij(2,2) = Xi(2,1)*R(33) + Xi(2,2)*R(34) + Xi(2,3)*R(35)
      Dij(2,6) = Xi(3,1)*R(33) + Xi(3,2)*R(34) + Xi(3,3)*R(35)

      Dij(2,3) = Xi(3,1)*R(36) + Xi(3,2)*R(37) + Xi(3,3)*R(38)
cccccccccc
      if (ldebug) then
         Dijp(2,4) =Xi(1,1)*R(33) + Xi(1,2)*R(34) + Xi(1,3)*R(35)
         Dijp(2,5) = Xi(1,1)*R(36) + Xi(1,2)*R(37) + Xi(1,3)*R(38)
         Dijp(2,6) = Xi(2,1)*R(36) + Xi(2,2)*R(37) + Xi(2,3)*R(38)
         print*," ratios of 2 versions for D24, D25, D26 "
         print*," D24 ",dij(2,4),dij(2,4)/dijp(2,4)
         print*," D25 ",dij(2,5),dij(2,5)/dijp(2,5)
         print*," D26 ",dij(2,6),dij(2,6)/dijp(2,6)
      endif
cccccccccc
c      do j = 8,13
c         Dij(2,j) = 0
c      enddo
c
      Dij(3,11) = -mt**2/2d0*dij(1,1) 
     &            - (f1*Dij(2,1)+f2*dij(2,4)+f3*dij(2,5)+C0234)/4d0
      Dij(3,12) = -mt**2/2d0*dij(1,2) 
     &            - (f1*Dij(2,4)+f2*dij(2,2)+f3*dij(2,6)-Cij234(1,1))/4d0
      Dij(3,13) = -mt**2/2d0*dij(1,3) 
     &            - (f1*Dij(2,5)+f2*dij(2,6)+f3*dij(2,3)-Cij234(1,2))/4d0
c
      R(41) = (f1*Dij(2,1)+Cij134(2,1)-C0234)/2d0 - 2*Dij(3,11)
      R(42) = (f2*Dij(2,1)+Cij124(2,1)-Cij134(2,1))/2d0
      R(43) = (f3*Dij(2,1)+Cij123(2,1)-Cij124(2,1))/2d0

      R(44) = (f1*Dij(2,4)+Cij134(2,1)+Cij234(1,1))/2d0 - Dij(3,12)
      R(50) = (f1*Dij(2,2)+Cij134(2,1)-Cij234(2,1))/2d0 
      R(56) = (f1*Dij(2,3)+Cij134(2,2)-Cij234(2,2))/2d0

      R(45) = (f2*Dij(2,4)+Cij124(2,3)-Cij134(2,1))/2d0 - Dij(3,11)
      R(51) = (f2*Dij(2,2)+Cij124(2,2)-Cij134(2,1))/2d0 - 2*Dij(3,12)
      R(57) = (f2*Dij(2,3)+Cij124(2,2)-Cij134(2,2))/2d0 

      R(46) = (f3*Dij(2,4)+Cij123(2,3)-Cij124(2,3))/2d0
      R(52) = (f3*Dij(2,2)+Cij123(2,2)-Cij124(2,2))/2d0
      R(58) = (f3*Dij(2,3)            -Cij124(2,2))/2d0 - 2*Dij(3,13) 

      Dij(3,1) = Xi(1,1)*R(41) + Xi(1,2)*R(42) + Xi(1,3)*R(43)
      Dij(3,4) = Xi(2,1)*R(41) + Xi(2,2)*R(42) + Xi(2,3)*R(43)
      Dij(3,5) = Xi(3,1)*R(41) + Xi(3,2)*R(42) + Xi(3,3)*R(43)

      Dij(3,6) = Xi(1,1)*R(50) + Xi(1,2)*R(51) + Xi(1,3)*R(52)
      Dij(3,2) = Xi(2,1)*R(50) + Xi(2,2)*R(51) + Xi(2,3)*R(52)
      Dij(3,8) = Xi(3,1)*R(50) + Xi(3,2)*R(51) + Xi(3,3)*R(52)

      Dij(3,7) = Xi(1,1)*R(56) + Xi(1,2)*R(57) + Xi(1,3)*R(58)
      Dij(3,9) = Xi(2,1)*R(56) + Xi(2,2)*R(57) + Xi(2,3)*R(58)
      Dij(3,3) = Xi(3,1)*R(56) + Xi(3,2)*R(57) + Xi(3,3)*R(58)

      Dij(3,10) = Xi(3,1)*R(44) + Xi(3,2)*R(45) + Xi(3,3)*R(46)
*----------------------------------------------------------------------
c WARNING: Previously unused Dij(1,7-13) variables are used
c to define D00ij and D0000 functions for DD-E-coefficients in 
c amplitudes/ggf/ggf_amp_aux.F --> subroutine DDtens51m
c In PV notation we have:
c Dij(1,7)=D416, Dij(1,8)=D417, Dij(1,9)=D418
c Dij(1,10)=D419, Dij(1,11)=D420, Dij(1,12)=D421, D(1,13)=D422
      msq=mt*mt
      Dij(1,7)=(-C0234+f1*Dij(3,1)+f2*Dij(3,4)+f3*Dij(3,5)+2d0
     &        *msq*Dij(2,1))/6d0
      Dij(1,8)=(-Cij234(2,1)+f1*Dij(3,6)+f2*Dij(3,2)+f3*Dij(3,8)+2d0
     &        *msq*Dij(2,2))/6d0
      Dij(1,9)=(-Cij234(2,2)+f1*Dij(3,7)+f2*Dij(3,9)+f3*Dij(3,3)+2d0
     &        *msq*Dij(2,3))/6d0
      Dij(1,10)=(Cij234(1,1)+f1*Dij(3,4)+f2*Dij(3,6)+f3*Dij(3,10)+2d0
     &         *msq*Dij(2,4))/6d0
      Dij(1,11)=(Cij234(1,2)+f1*Dij(3,5)+f2*Dij(3,10)+f3*Dij(3,7)+2d0
     &         *msq*Dij(2,5))/6d0
      Dij(1,12)=(-Cij234(2,3)+f1*Dij(3,10)+f2*Dij(3,8)+f3*Dij(3,9)+2d0
     &         *msq*Dij(2,6))/6d0
      Dij(1,13)=(Cij234(2,4)+1d0/12d0-p1sq*Dij(1,7)-p2sq*Dij(1,8)-p3sq
     &         *Dij(1,9)-2d0*p1p2*Dij(1,10)-2d0*p1p3*Dij(1,11)-2d0*p2p3
     &         *Dij(1,12)-msq*Dij(2,7))/6d0 
*----------------------------------------------------------------------
      if (ldb0) then
         print*," R_20/R_22 = ",r(20)/r(22)
         print*," R_30/R_32 = ",r(30)/r(32)
         print*," R_33/R_35 = ",r(33)/r(35)
         print*," R_36/R_38 = ",r(36)/r(38)

         print*," R_41/R_43 = ",r(41)/r(43)
         print*," R_50/R_52 = ",r(50)/r(52)
         print*," R_56/R_58 = ",r(56)/r(58)
         print*," R_44/R_46 = ",r(44)/r(46)
      endif
      if (ldebug) then
         Dijp(3,4) = Xi(1,1)*R(44) + Xi(1,2)*R(45) + Xi(1,3)*R(46)
         Dijp(3,6) = Xi(2,1)*R(44) + Xi(2,2)*R(45) + Xi(2,3)*R(46)
         print*," ratios of 2 versions for D3j "
         print*," D34 ",dij(3,4),dij(3,4)/dijp(3,4)
         print*," D36 ",dij(3,6),dij(3,6)/dijp(3,6)
         t(1) = (cij134(2,3)+f1*Dij(2,5)+cij234(1,2))/2-dij(3,13)
         t(2) = (cij134(2,3)+f1*Dij(2,6)-cij234(2,3))/2
         t(3) = (cij134(2,4)+f1*Dij(2,7)-cij234(2,4))/2

         t(4) = (cij124(2,3)-cij134(2,3)+f2*Dij(2,5))/2
         t(5) = (cij124(2,2)-cij134(2,3)+f2*Dij(2,6))/2-dij(3,13)
         t(6) = (cij124(2,4)-cij134(2,4)+f2*Dij(2,7))/2

         t(7) = (           -cij124(2,3)+f3*Dij(2,5))/2-dij(3,11)
         t(8) = (           -cij124(2,2)+f3*Dij(2,6))/2-dij(3,12)
         t(9) = (cij123(2,4)-cij124(2,4)+f3*Dij(2,7))/2

         Dijp(3,5) = Xi(1,1)*t(1) + Xi(1,2)*t(4) + Xi(1,3)*t(7)
         Dijp(3,10)= Xi(2,1)*t(1) + Xi(2,2)*t(4) + Xi(2,3)*t(7)
         Dijp(3,7) = Xi(3,1)*t(1) + Xi(3,2)*t(4) + Xi(3,3)*t(7)

         D310pp    = Xi(1,1)*t(2) + Xi(1,2)*t(5) + Xi(1,3)*t(8)
         Dijp(3,8) = Xi(2,1)*t(2) + Xi(2,2)*t(5) + Xi(2,3)*t(8)
         Dijp(3,9) = Xi(3,1)*t(2) + Xi(3,2)*t(5) + Xi(3,3)*t(8)
         
         Dijp(3,11)= Xi(1,1)*t(3) + Xi(1,2)*t(6) + Xi(1,3)*t(9)
         Dijp(3,12)= Xi(2,1)*t(3) + Xi(2,2)*t(6) + Xi(2,3)*t(9)
         Dijp(3,13)= Xi(3,1)*t(3) + Xi(3,2)*t(6) + Xi(3,3)*t(9)

         print*," D35 ",dij(3,5),dij(3,5)/dijp(3,5)
         print*," D37 ",dij(3,7),dij(3,7)/dijp(3,7)
         print*," D38 ",dij(3,8),dij(3,8)/dijp(3,8)
         print*," D39 ",dij(3,9),dij(3,9)/dijp(3,9)
         print*," D310 ",dij(3,10),dij(3,10)/dijp(3,10)
         print*," D310pp ",d310pp,d310pp/dij(3,10)
         print*," D311 ",dij(3,11),dij(3,11)/dijp(3,11)
         print*," D312 ",dij(3,12),dij(3,12)/dijp(3,12)
         print*," D313 ",dij(3,13),dij(3,13)/dijp(3,13)
      endif
cccccccccccc
      return
      end


c----  B0t(qsq): regularized 2-point function for massless loop  -----
c
      complex*16 function B0t(qsq,musq)
      implicit none
      double precision qsq, musq

c evaluate the finite piece of the massless B0 function
c  
c             d^d k 
c    B0 = Int------- [-k^2+idelta]^-1 [-(k-q)^2+idelta]^-1 
c             i pi^2
c
c       = pi^-eps Gamma(1+eps) musq^-eps [1/eps + B0t(qsq,musq)]
c
c i.e. B0t(qsq,musq) = 2 - log( (-qsq-idelta)/musq )
c
c subroutine requires positive musq
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2002 September 13
c	Last modified:    2008 by Guiseppe Bozzi
c  
      complex*16 ipi
      parameter (ipi=(0d0,3.14159 26535 89793d0))
      logical ldebug
      parameter (ldebug=.true.)
      real *8 tiny
      parameter (tiny=1d-7)
c
      if (ldebug) then
         if (musq.le.0d0) then
            print*," Unacceptable mu^2 in B0t(q^2,mu^2): mu^2 = ",
     1             musq," ==> set B0t = 0 "
            b0t=0
            return
         endif
      endif
c
      if (abs(qsq).lt.tiny) then
         b0t = 0d0
      elseif (qsq.lt.0d0) then
         b0t = 2d0 - log(-qsq/musq)
      else
         b0t = 2d0 - log(qsq/musq) + ipi
      endif
      end


c
c------------- C0(q1^2,q2^2,P^2)  massless 3-point fuction ------------
c
      complex*16 function C0t(q1sq,q2sq,Psq,musq)
      implicit none
      double precision q1sq, q2sq, psq, musq

c evaluate the finite part of the scalar 3-point function for zero mass
c propagators 
c  
c  C0 = 1/(i*pi^2) * Int d^dk [-k^2-i0]^-1 [-(k+q1)^2-i0]^-1 [-(k-q2)^2-i0]^-1
c
c       = pi^-eps Gamma(1+eps) musq^-eps [IR + C0t(q1sq,q2sq,psq)]
c
c where IR represents the potentially divergent terms when one or more
c of the arguments vanish. P = q1+q2. The divergent terms are
c 
c (1) IR = 0              all 3 of q1sq,q2sq,psq nonzero
c 
c (2) IR = 1/(t-qsq)*log(t/qsq)*1/eps    for 1 vanishing argument, the
c                                        other two being t and qsq
c
c (3) IR = -1/qsq*(1/eps**2 - 1/eps log(qsq/musq))  when exactly one 
c                                                   argument, qsq,
c                                                   is nonzero
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2002 September 13
c	Last modified:    2008 by Guiseppe Bozzi
c  
      complex*16 ipi, z, res, vli2
      parameter (ipi=(0d0,3.14159 26535 89793d0))
      external vli2
      double precision pi,pi2o3,pi2o6,pi2o2,tiny
      parameter (pi=3.14159 26535 89793d0,pi2o3=pi**2/3d0)
      parameter (pi2o6=pi**2/6d0,pi2o2=pi**2/2d0,tiny=1d-3)
      logical ldebug, ltaylor, first_call
      parameter (ldebug=.false.)
      double precision x, flambda, tsqlami, qsq(3), 
     1                 sum, rz, logz, ln1, ln2,
     2                 musq_t, asum, fac, resr, ser, fb(8), xn
      integer izero, inz, i, j
      data first_call /.true./
      save first_call, fb
c
      if (ldebug) then
         if (musq.le.0d0) then
            print*," Unacceptable mu^2 in C0t(q1sq,q2sq,psq,mu^2): ",
     1         "mu^2 = ",musq," ==> set C0t = 0 "
            C0t=0
            return
         endif
      endif
c
      inz = 0
      izero = 0
      do i=1,3
         qsq(i)=0
      enddo
c      if (q1sq.eq.0d0) then
c         izero = izero+1
c      else
c         inz = inz+1
c         qsq(inz)=q1sq
c      endif
c      if (q2sq.eq.0d0) then
c         izero = izero+1
c      else
c         inz = inz+1
c         qsq(inz)=q2sq
c      endif
c      if (psq.eq.0d0) then
c         izero = izero+1
c      else
c         inz = inz+1
c         qsq(inz)=psq
c      endif
      if (abs(q1sq).lt.tiny) then
         izero = izero+1
      else
         inz = inz+1
         qsq(inz)=q1sq
      endif
      if (abs(q2sq).lt.tiny) then
         izero = izero+1
      else
         inz = inz+1
         qsq(inz)=q2sq
      endif
      if (abs(psq).lt.tiny) then
         izero = izero+1
      else
         inz = inz+1
         qsq(inz)=psq
      endif
      if (ldebug) then
         print("(a,3g12.4,a,g12.4)"), "qsq_i = ",q1sq,q2sq,psq,
     1                                " mu^2 = ",musq
         print*," number of q_i^2 = 0  ",izero
         print*," number of q_i^2 ne 0 ",inz," Sum = ",inz+izero
         if ((inz+izero).ne.3) print*," WARNING: inz+izero .ne. 3 "
      endif
c
      if (izero.eq.0) then   ! all qsq(i) nonzero, this is full C0
         flambda = (qsq(1)**2+qsq(2)**2+qsq(3)**2) - 
     1           2*(qsq(1)*qsq(2)+qsq(1)*qsq(3)+qsq(2)*qsq(3))
         sum = qsq(1)+qsq(2)+qsq(3)
         asum = abs(qsq(1))+abs(qsq(2))+abs(qsq(3))
         ltaylor = abs(abs(sum)/asum-1).lt.1d-13      ! same sign for all qsq
         asum=min(abs(sum-2*qsq(1)),abs(sum-2*qsq(2)),
     1            abs(sum-2*qsq(3)))
         ltaylor = ltaylor .and. abs(flambda).lt.(0.01*asum**2)
c
         res = (0d0,0d0)
         if (ldebug) then
            print*," lambda = ",flambda/4d0," x_taylor = ",
     1                          flambda/asum**2
            ltaylor=.false.
         endif
         if (flambda.gt.0d0 .and. .not.ltaylor) then
            tsqlami = 1d0/sqrt(flambda)
            do i=1,3
               x = (sum-2*qsq(i))*tsqlami
               if (x.gt.1d0) then
                  rz = (x-1)/(x+1)
                  z = rz
                  logz = log(rz)
                  res = res + 2*Vli2(z)+0.5d0*logz**2-pi2o3 
                  if (qsq(i).gt.0d0) then
                     res = res + ipi*logz
                  else
                     res = res - ipi*logz
                  endif
               elseif (x.gt.0d0) then
                  rz = (x-1)/(x+1)
                  z = rz
                  logz = log(-rz)
                  res = res + 2*Vli2(z)+0.5d0*logz**2+pi2o6
               elseif (x.gt.-1d0) then
                  rz = (x+1)/(x-1)
                  z = rz
                  logz = -log(-rz)
                  res = res - 2*Vli2(z)-0.5d0*logz**2-pi2o6
               else 
                  rz = (x+1)/(x-1)
                  z = rz
                  logz = -log(rz)
                  res = res - 2*Vli2(z)-0.5d0*logz**2+pi2o3 
                  if (qsq(i).gt.0d0) then
                     res = res + ipi*logz
                  else
                     res = res - ipi*logz
                  endif
               endif
               if (ldebug) print*," i,x, z = ",i,x,rz
            enddo
            C0t = res*tsqlami
         elseif (flambda.lt.0d0 .and. .not.ltaylor) then  ! the lambda < 0 case
            tsqlami = 1d0/sqrt(-flambda)         ! + imag part of 2 sqrt(lambda)
            do i=1,3
               x = (sum-2*qsq(i))*tsqlami
               z = dcmplx(x,-1d0)/dcmplx(x,1d0)
               if (ldebug) then  
                  print*," i,x, z = ",i,x,z
               endif
               res = res + dimag(Vli2(z))
            enddo
            C0t = res*2*tsqlami
         endif
         if (ldebug) then
            print*," C0t via Vli2 calls: ",C0t
            ltaylor = abs(flambda).lt.(0.5*asum**2)
         endif
         if (ltaylor) then
            if (first_call) then
               fb(1) = 2d0/9d0
               fb(2) = 16d0/75d0
               fb(3) = 142d0/735d0
               fb(4) = 496d0/2835d0
               fb(5) = 6086d0/38115d0
               fb(6) = 86048d0/585585d0
               fb(7) = 92054d0/675675d0
               fb(8) = 1655008d0/13018005d0 
c               fb(9) = 32976682d0/276441165d0
               first_call = .false.
            endif
            musq_t = abs(sum/3d0)
            resr = 0
            do i = 1,3
               fac = 1/(sum-2*qsq(i))
               x = flambda*fac**2
c               if (ldebug) print*," i = ",i," x = ",x
               ln1 = log(abs(qsq(i))/musq_t)
               ser = ln1
               j = 1
               xn = x
c               if (ldebug) print*," log = ",ln1
               do while (abs(xn).gt.1d-15 .and. j.lt.9)
                  ser = ser + xn*(ln1/(2*j+1)-fb(j))
c                  if (ldebug) then
c                     print*," n = ",j," a_n = ",ln1/(2*j+1)-fb(j)
c                  endif
                  xn = xn*x
                  j = j+1
               enddo
               if (ldebug) print*," C_i = ",2*ser*fac
               resr = resr + 2*ser*fac
            enddo
            if (ldebug) then
               print*," C0t via Li2 calls: ",C0t
               print*," C0t via Taylor exp.",dcmplx(resr,0d0)
            endif
            C0t = dcmplx(resr,0d0)
         endif
      elseif (izero.eq.1) then
         ln1 = log(abs(qsq(1))/musq)
         ln2 = log(abs(qsq(2))/musq)
         res = 0.5d0*(ln1**2-ln2**2)
         if (qsq(1).gt.0d0) then
            res = res - ipi*ln1 - pi2o2
         endif
         if (qsq(2).gt.0d0) then
            res = res + ipi*ln2 + pi2o2
         endif
         C0t = res/(qsq(2)-qsq(1))
      elseif (izero.eq.2) then
         ln1 = log(abs(qsq(1))/musq)
         res = 0.5d0*ln1**2-pi2o6
         if (qsq(1).gt.0d0) then
            res = res - ipi*ln1 - pi2o2
         endif
         C0t = -res/qsq(1)
      else
         print*," WARNING: C0t called with 3 zero arguments "
         print*," C0t reset to 0 "
         C0t = 0
      endif
      if (ldebug) then
         print*," output: C0t = ",C0t
      endif
      return 
      end
c


********************************************************************************
c--------------  B0tM(m,qsq): regularized 2-point function --------------

      complex*16 function B0tM(m,qsq) 
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
         srt = dsqrt(1d0-4d0/qsqn)
         lnfac = dlog( (srt-1d0)/(srt+1d0) )
         B0tM = 2d0 + srt*lnfac
      elseif ( abs(qsqn).le.eps )then
         B0tM = qsqn/6d0* ( 1d0+qsqn*0.1d0*(1d0+qsqn/7d0 *
     &        ( 1d0+qsqn/7d0*( 1d0+2d0/11d0*qsqn ) ) ) )
      elseif (qsqn.lt.4d0) then
         srt = dsqrt(4d0/qsqn-1d0)
         phi = atan(1d0/srt)
         B0tM = 2d0 - 2d0*srt*phi
      elseif (qsqn.eq.4d0) then
         B0tM = 2d0
      else
         beta = dsqrt(1d0-4d0/qsqn)
         lnfac = dlog( (1d0-beta)/(1d0+beta) )
         re = 2d0 + beta*lnfac
         im = pi*beta
         B0tM = dcmplx( re, im )
      endif
      return
      end


********************************************************************************
********************************************************************************


