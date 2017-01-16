****************************************************************************
*  file for dot-products and tensor contractions
****************************************************************************
*   This file contains functions for dot-products of real and complex 
*   4-vectors and routines for contractions with the Levi-Civita 
*   epsilon-tensor and arbitrary tensors.
c***************************************************************************
c   LIST OF ALL FUNCTIONS AND SUBROUTINES IN THIS FILE:
c
c   double precision FUNCTION dotrr(p1,p2)
c   double complex   FUNCTION dotcc(v1,v2)
c   double complex   FUNCTION dotrc(v1,v2)
c   double precision FUNCTION dot0p(p1,p2)
c   double complex   FUNCTION dotqj(j1)
c   real*8           FUNCTION mjj2(p1,p2)
c   real*8           FUNCTION qjj2(p1,p2)
c   real*8           FUNCTION mass2(p)
c   double precision FUNCTION psumsq(p1,p2)
c   double complex   FUNCTION contract_Tjj(T,j1,j2)
c   complex*16       FUNCTION contract_TJJr(T,j1,j2)
c   real*8           FUNCTION contract_Trjjr(T,j1,j2)
c   SUBROUTINE EPSRRR(RES,A,B,C)
c   SUBROUTINE EPSCRR(RES,A,B,C)
c   SUBROUTINE EPSCCR(RES,A,B,C)
c   SUBROUTINE EPSCCC(RES,A,B,C)
c   SUBROUTINE EPSRR (T,A,B)
c   SUBROUTINE contract_CPO(J1cJ2,J1,J2,k1,k2)
c   SUBROUTINE contract_CPE(J1cJ2,J1,J2,k1,k2)
c   SUBROUTINE contract_T1j(T,jc,jout)
c   SUBROUTINE contract_T2j(T,jc,jout)
c   SUBROUTINE contract_Tr1j(T,jc,jout)
c   SUBROUTINE contract_T1jr(T,jr,jout)
c   SUBROUTINE contract_Tr2j(T,jc,jout)
c   SUBROUTINE contract_T2jr(T,jr,jout)
c   SUBROUTINE contract_Tr3j(T,jc,i,Tout)
c   SUBROUTINE contract_Tr3jr(T,jr,i,Tout)
c   SUBROUTINE contract_T3j(T,jc,i,Tout)
c   SUBROUTINE contract_T3jr(T,jr,i,Tout)
c   SUBROUTINE contract_Tr4j(T,jc,i,Tout)
c   SUBROUTINE contract_T4j(T,jc,i,Tout)
c   SUBROUTINE contract_Tr4jr(T,jr,i,Tout)
cMK sign convention of eps
cMK The convention is eps_abcd v^a v^b v^c v^d with eps_0123 = +1
cMK this differs compared to the conventions in FeynCalc and FormCalc, where eps^0123 = +1
c   double complex   FUNCTION EPSCRRR(a,b,c,d)
c   double complex   FUNCTION EPSCCRR(a,b,c,d)
c   double complex   FUNCTION EPSCCCR(a,b,c,d)
c   double complex   FUNCTION EPSCCCC(a,b,c,d)
c   double precision FUNCTION EPSRRRR(a,b,c,d)
c   double complex   FUNCTION EPSCCRMU(a,b,c,mu)
c   double complex   FUNCTION EPSCRRMU(a,b,c,mu)
c   double precision FUNCTION EPSRRRMU(a,b,c,mu)
c   double complex   FUNCTION EPSCRMUNU(a,b,mu,nu)
c   double precision FUNCTION EPSRRMUNU(a,b,mu,nu)
c   double precision FUNCTION EPSRMUNURHO(a,mu,nu,rho)
c   double real*8    FUNCTION EPSMUNURHOSIGMA(mu,nu,rho,sigma)
c
c***************************************************************************
*   Last modified: 06.07.2006
****************************************************************************


c***************************************************************************
      function dotrr(p1,p2)
c***************************************************************************
c
c     dotrr(p1,p2) = p1.p2
c
c***************************************************************************
      implicit none

      double precision dotrr,p1(0:3),p2(0:3)

      dotrr = p1(0)*p2(0) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)

      end


c***************************************************************************
      double complex function dotcc(v1,v2)
c***************************************************************************
c
c     dotcc(v1,v2) = v1.v2
c
c***************************************************************************
      implicit none
      double complex v1(0:3), v2(0:3)
      dotcc = v1(0)*v2(0)-v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)
      end


c***************************************************************************
      double complex function dotcr(v1,v2)
c***************************************************************************
c
c     dotrc(v1,v2) = v1.v2
c
c***************************************************************************
      implicit none
      double precision v2(0:3)
      double complex  v1(0:3)
      dotcr = v1(0)*v2(0)-v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)
      end


c***************************************************************************
      double complex function dotrc(v1,v2)
c***************************************************************************
c
c     dotrc(v1,v2) = v1.v2
c
c***************************************************************************
      implicit none
      double precision v1(0:3)
      double complex  v2(0:3)
      dotrc = v1(0)*v2(0)-v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)
      end


c***************************************************************************
      function dot0p(p1,p2)
c***************************************************************************
c
c     dot0p(p1,p2) = p1.p2
c
c***************************************************************************
      implicit none
      double precision dot0p,p1(0:3),p2(0:3)
      dot0p = p1(0)*p2(0) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)
c$$$      if (abs(dot0p).lt.0.01*abs(p1(0)*p2(0))) then !section for small p1.p2
c$$$*         print*," old dotp = ",dot0p," ratio = ",dot0p/p1(0)/p2(0)
c$$$         fac = p2(0)/p1(0)
c$$$         pt(1)=p2(1)-fac*p1(1)
c$$$         pt(2)=p2(2)-fac*p1(2)
c$$$         pt(3)=p2(3)-fac*p1(3)
c$$$         dot0p = (pt(1)**2+pt(2)**2+pt(3)**2)/(2*fac)
c$$$*         print*," new dotp = ",dot0p
c$$$      end if

      end


c***************************************************************************
      double complex function dotqj(v2)
c***************************************************************************
c
c     dotqj(j1) = j1.p1, where p1 is four-momentum associated with current j1
c
c***************************************************************************
      implicit none
      double precision v1(0:3)
      double complex  v2(0:5)
      v1(0) = dreal(v2(4))
      v1(1) = dreal(v2(5))
      v1(2) = dimag(v2(5))
      v1(3) = dimag(v2(4))
      dotqj = v1(0)*v2(0)-v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)
      end


c***************************************************************************
            real*8 function mjj2(p1,p2)
c***************************************************************************
c
c	This function calculates the invariant mass
c	square of the two four-momenta as given.
c	ie.,  mjj2 = (p1 + p2)^2
c
c       input: real*8 p1(0:3), p2(0:3) four momenta
c
c***************************************************************************
      implicit none
c declare input/output variables
      real*8 p1(0:3), p2(0:3)
c calculate
      mjj2 = ( ( p1(0) + p2(0) )**2
     &       - ( p1(1) + p2(1) )**2
     &       - ( p1(2) + p2(2) )**2
     &       - ( p1(3) + p2(3) )**2 )
      return
      end


c***************************************************************************
            real*8 function qjj2(p1,p2)
c***************************************************************************
c
c	This function calculates the invariant momentum
c	transfer square of the two four-momenta as given.
c	ie.,  qjj2 = (p1 - p2)^2
c
c       input: real*8 p1(0:3), p2(0:3) four momenta
c
c***************************************************************************
      implicit none
c declare input/output variables
      real*8 p1(0:3), p2(0:3)
c calculate
      qjj2 = ( ( p1(0) - p2(0) )**2
     &       - ( p1(1) - p2(1) )**2
     &       - ( p1(2) - p2(2) )**2
     &       - ( p1(3) - p2(3) )**2 )
      return
      end


c***************************************************************************
      real*8 FUNCTION mass2(p)
c***************************************************************************
c     calculates p^2
c***************************************************************************
      implicit none
      real*8 p(0:3)
      mass2= p(0)*p(0)-p(1)*p(1)-p(2)*p(2)-p(3)*p(3)
      return
      end


c***************************************************************************
      double precision FUNCTION psumsq(p1,p2)
c***************************************************************************
c     calculates (p1 + p2)^2
c***************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3)
      psumsq = (p1(0)+p2(0))**2 - (p1(1)+p2(1))**2 -
     1         (p1(2)+p2(2))**2 - (p1(3)+p2(3))**2
      return
      end


c***************************************************************************
      SUBROUTINE EPSRRR(RES,A,B,C)
c***************************************************************************
c
c     This subroutine calculates the four vector
c     RES^mu = EPS^mu_alpha_beta_rho * A^alpha * B^beta * C^rho
c     where EPS is the totally antisymmetric tensor in Bjorken-Drell
c     convention: EPS_0123=+1
c
c     input:  double precision A(0:3), B(0:3), C(0:3)
c     output: double precision res(0:3) (= EPS[mu,A,B,C])
c
c***************************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION  RES(0:3), A(0:3)
      DOUBLE PRECISION  B(0:3),C(0:3),BC12,BC23,BC31,BC01,BC02,BC03

      BC12 = B(1)*C(2) - B(2)*C(1)
      BC23 = B(2)*C(3) - B(3)*C(2)
      BC31 = B(3)*C(1) - B(1)*C(3)
      BC01 = B(0)*C(1) - B(1)*C(0)
      BC02 = B(0)*C(2) - B(2)*C(0)
      BC03 = B(0)*C(3) - B(3)*C(0)
      RES(0) = A(1)*BC23 + A(2)*BC31 + A(3)*BC12
      RES(1) = -A(2)*BC03 + A(3)*BC02 + A(0)*BC23
      RES(2) = A(1)*BC03 + A(0)*BC31 - A(3)*BC01
      RES(3) = -A(1)*BC02 + A(2)*BC01 + A(0)*BC12

      RETURN
      END


c***************************************************************************
      SUBROUTINE EPSCRR(RES,A,B,C)
c***************************************************************************
c
c     This subroutine calculates the four vector
c     RES^mu = EPS^mu_alpha_beta_rho * A^alpha * B^beta * C^rho
c     where EPS is the totally antisymmetric tensor in Bjorken-Drell
c     convention: EPS_0123=+1
c
c     input:  double complex A(0:3)
c             double precision B(0:3), C(0:3)
c     output: double complex res(0:3) (= EPS[mu,A,B,C])
c
c***************************************************************************
      IMPLICIT NONE

      DOUBLE COMPLEX  RES(0:3), A(0:3)
      DOUBLE PRECISION  B(0:3),C(0:3),BC12,BC23,BC31,BC01,BC02,BC03

      BC12 = B(1)*C(2) - B(2)*C(1)
      BC23 = B(2)*C(3) - B(3)*C(2)
      BC31 = B(3)*C(1) - B(1)*C(3)
      BC01 = B(0)*C(1) - B(1)*C(0)
      BC02 = B(0)*C(2) - B(2)*C(0)
      BC03 = B(0)*C(3) - B(3)*C(0)
      RES(0) = A(1)*BC23 + A(2)*BC31 + A(3)*BC12
      RES(1) = -A(2)*BC03 + A(3)*BC02 + A(0)*BC23
      RES(2) = A(1)*BC03 + A(0)*BC31 - A(3)*BC01
      RES(3) = -A(1)*BC02 + A(2)*BC01 + A(0)*BC12

      RETURN
      END


c***************************************************************************
      SUBROUTINE EPSCCR(RES,A,B,C)
c***************************************************************************
c
c     This subroutine calculates the four vector
c     RES^mu = EPS^mu_alpha_beta_rho * A^alpha * B^beta * C^rho
c     where EPS is the totally antisymmetric tensor in Bjorken-Drell
c     convention: EPS_0123=+1
c
c     input:  double complex A(0:3), B(0:3)
c             double precision C(0:3)
c     output: double complex res(0:3) (= EPS[mu,A,B,C])
c
c***************************************************************************
      IMPLICIT NONE

      DOUBLE COMPLEX  RES(0:3), A(0:3), B(0:3),BC12,BC23,BC31,BC01,BC02,BC03
      DOUBLE PRECISION  C(0:3)

      BC12 = B(1)*C(2) - B(2)*C(1)
      BC23 = B(2)*C(3) - B(3)*C(2)
      BC31 = B(3)*C(1) - B(1)*C(3)
      BC01 = B(0)*C(1) - B(1)*C(0)
      BC02 = B(0)*C(2) - B(2)*C(0)
      BC03 = B(0)*C(3) - B(3)*C(0)
      RES(0) = A(1)*BC23 + A(2)*BC31 + A(3)*BC12
      RES(1) = -A(2)*BC03 + A(3)*BC02 + A(0)*BC23
      RES(2) = A(1)*BC03 + A(0)*BC31 - A(3)*BC01
      RES(3) = -A(1)*BC02 + A(2)*BC01 + A(0)*BC12

      RETURN
      END


c***************************************************************************
      SUBROUTINE EPSCCC(RES,A,B,C)
c***************************************************************************
c
c     This subroutine calculates the four vector
c     RES^mu = EPS^mu_alpha_beta_rho * A^alpha * B^beta * C^rho
c     where EPS is the totally antisymmetric tensor in Bjorken-Drell
c     convention: EPS_0123=+1
c
c     input:  double complex A(0:3), B(0:3), C(0:3)
c     output: double complex res(0:3) (= EPS[mu,A,B,C])
c
c***************************************************************************
      IMPLICIT NONE

      DOUBLE COMPLEX  RES(0:3), A(0:3), B(0:3), C(0:3)
      DOUBLE COMPLEX  BC12, BC23, BC31, BC01, BC02, BC03

      BC12 = B(1)*C(2) - B(2)*C(1)
      BC23 = B(2)*C(3) - B(3)*C(2)
      BC31 = B(3)*C(1) - B(1)*C(3)
      BC01 = B(0)*C(1) - B(1)*C(0)
      BC02 = B(0)*C(2) - B(2)*C(0)
      BC03 = B(0)*C(3) - B(3)*C(0)
      RES(0) = A(1)*BC23 + A(2)*BC31 + A(3)*BC12
      RES(1) = -A(2)*BC03 + A(3)*BC02 + A(0)*BC23
      RES(2) = A(1)*BC03 + A(0)*BC31 - A(3)*BC01
      RES(3) = -A(1)*BC02 + A(2)*BC01 + A(0)*BC12

      RETURN
      END


c***************************************************************************
c      double precision function EPSRRRR(d,a,b,c)
c***************************************************************************
c     This subroutine calculates the complete contraction
c     EPS_mu_alpha_beta_rho * D^mu * A^alpha * B^beta * C^rho
c     where EPS is the totally antisymmetric tensor in Bjorken-Drell
c     convention: EPS_0123=+1
c
c***************************************************************************
c      implicit none
c      double precision a(0:3), b(0:3), c(0:3), d(0:3), res(0:3)
c      DOUBLE PRECISION BC12,BC23,BC31,BC01,BC02,BC03
cc
c      BC12 = B(1)*C(2) - B(2)*C(1)
c      BC23 = B(2)*C(3) - B(3)*C(2)
c      BC31 = B(3)*C(1) - B(1)*C(3)
c      BC01 = B(0)*C(1) - B(1)*C(0)
c      BC02 = B(0)*C(2) - B(2)*C(0)
c      BC03 = B(0)*C(3) - B(3)*C(0)
c      RES(0) = A(1)*BC23 + A(2)*BC31 + A(3)*BC12
c      RES(1) = -A(2)*BC03 + A(3)*BC02 + A(0)*BC23
c      RES(2) = A(1)*BC03 + A(0)*BC31 - A(3)*BC01
c      RES(3) = -A(1)*BC02 + A(2)*BC01 + A(0)*BC12
ccc
c      epsrrrr = d(0)*res(0)-d(1)*res(1)-d(2)*res(2)-d(3)*res(3)

c      END


c***************************************************************************
      SUBROUTINE  EPSRR (T,A,B)
c***************************************************************************
c     This subroutine calculates the second rank tensor      
c     T^mu^nu = EPS^mu^nu_alpha_beta * A^alpha * B^beta *
c     where EPS is the totally antisymmetric tensor in Bjorken-Drell
c     convention: EPS_0123=+1
c     input:  real*8 A(0:3), B(0:3)
c     output: real*8 T(0:3,0:3)
*
c***************************************************************************

      implicit none
      
      real*8 T(0:3,0:3), a(0:3), b(0:3)
      integer i
      
      do i = 0,3
         t(i,i) = 0
      enddo
      
      t(0,1) = a(2) * b(3) - a(3) * b(2)
      t(0,2) = a(3) * b(1) - a(1) * b(3)
      t(0,3) = a(1) * b(2) - a(2) * b(1)
      t(1,2) = - (a(0) * b(3) - a(3) * b(0))
      t(1,3) = - (a(2) * b(0) - a(0) * b(2))
      t(2,3) = - (a(0) * b(1) - a(1) * b(0))
      t(1,0) = a(3) * b(2) -a(2) * b(3)
      t(2,0) = a(1) * b(3) -a(3) * b(1)
      t(3,0) = a(2) * b(1) -a(1) * b(2)
      t(2,1) = - (a(3) * b(0) -a(0) * b(3))
      t(3,1) = - (a(0) * b(2) -a(2) * b(0))
      t(3,2) = - (a(1) * b(0) -a(0) * b(1))
      
      return
      end 


c***************************************************************************
      subroutine contract_CPO(J1cJ2,J1,J2,k1,k2)
c***************************************************************************
c
c     This subroutine computes 
c     J1cJ2 = Eps[mu,nu,rho,sigma] k1_rho k2_sigma J1_mu J2_nu
c
c     input:  double complex J1(0:3),J2(0:3)  4-current
c             real*8 k1(0:3),k2(0:3) 4-momenta
c     output: double complex J1cJ2 (= Eps[J1,J2,k1,k2])
c     
c***************************************************************************
      implicit none

      DOUBLE COMPLEX  J1(0:3),J2(0:3),J2p(0:3),J2pJ1
      double complex J1cJ2
      real*8 k1(0:3),k2(0:3)

c local variables


      call epscrr(J2p,J2,k1,k2)
      J2pJ1 = J2p(0)*J1(0)-J2p(1)*J1(1)-J2p(2)*J1(2)-J2p(3)*J1(3)
c
c     compute Eps[mu,nu.k1,k2] J1_mu J2_nu 
      J1cJ2 =  J2pJ1
      
      return
      end


c***************************************************************************
      subroutine contract_CPE(J1cJ2,J1,J2,k1,k2)
c***************************************************************************
c
c     This subroutine computes 
c     J1cJ2 = (g^(mu,nu) k1.k2 - k1^nu k2^mu) J1_mu J2_nu
c
c     input:  double complex J1(0:3),J2(0:3) 4-current
c             real*8 k1(0:3),k2(0:3) 4-momenta
c     output: double complex J1cJ2 (= (k1.k2) (J1.J2) - (k1.J2) (k2.J1))
c     
c***************************************************************************
      implicit none

      DOUBLE COMPLEX  J1(0:3),J2(0:3),J12,J2k1,J1k2
      double complex J1cJ2
      real*8 k1(0:3),k2(0:3),k12
 
c local variables


c     dot products
      k12 = k1(0)*k2(0)-k1(1)*k2(1)-k1(2)*k2(2)-k1(3)*k2(3)

      J12 = J1(0)*J2(0)-J1(1)*J2(1)-J1(2)*J2(2)-J1(3)*J2(3)

      J2k1 = J2(0)*k1(0)-J2(1)*k1(1)-J2(2)*k1(2)-J2(3)*k1(3)

      J1k2 = J1(0)*k2(0)-J1(1)*k2(1)-J1(2)*k2(2)-J1(3)*k2(3)

      J1cJ2 = (J12 * k12 - J2k1 * J1k2)
      
      return
      end


c***************************************************************************
      double complex function contract_Tjj(T,j1,j2)
c***************************************************************************
c     
c     contract complex rank 2 tensor T^{mu,nu} with 
c     two complex vectors j1 and j2
c
c     input: double complex rank 2 tensor T(0:3,0:3)
c            double complex four vectors j1(0:3),j2(0:3)
c
c***************************************************************************
      implicit none
      complex*16 T(0:3,0:3), j1(0:3), j2(0:3), resv(0:3)
      integer mu
      do mu = 0,3
         resv(mu) = T(mu,0)*j2(0) - T(mu,1)*j2(1) 
     &            - T(mu,2)*j2(2) - T(mu,3)*j2(3)
      enddo
      contract_Tjj = resv(0)*j1(0) - resv(1)*j1(1) 
     &             - resv(2)*j1(2) - resv(3)*j1(3)
      return
      end


c *************************************************************************     
      function contract_Tjjr(T,j1,j2)
c *************************************************************************     
      implicit none
c contract complex rank 2 tensor T^{mu,nu} with two real vectors j1 and j2
      complex*16 T(0:3,0:3),  resv(0:3),contract_Tjjr
      real*8 j1(0:3), j2(0:3)
      integer mu
      do mu = 0,3
         resv(mu) = T(mu,0)*j2(0) - T(mu,1)*j2(1) 
     &            - T(mu,2)*j2(2) - T(mu,3)*j2(3)
      enddo
      contract_Tjjr = resv(0)*j1(0) - resv(1)*j1(1) 
     &             - resv(2)*j1(2) - resv(3)*j1(3)
c      return
      end      
c ************************************************************************ 


c ************************************************************************
      function contract_Trjjr(T,j1,j2)
c ************************************************************************      
      implicit none
c contract real rank 2 tensor T^{mu,nu} with two real vectors j1 and j2
      real*8 T(0:3,0:3),  resv(0:3),contract_Trjjr
      real*8 j1(0:3), j2(0:3)
      integer mu
      do mu = 0,3
         resv(mu) = T(mu,0)*j2(0) - T(mu,1)*j2(1) 
     &            - T(mu,2)*j2(2) - T(mu,3)*j2(3)
      enddo
      contract_Trjjr = resv(0)*j1(0) - resv(1)*j1(1) 
     &             - resv(2)*j1(2) - resv(3)*j1(3)
c      return
      end     
c ************************************************************************      

      
c***************************************************************************
      subroutine contract_T1j(T,jc,jout)
c***************************************************************************
c
c     contract first index of complex rank 2 tensor T^{mu,nu} with 
c     complex vector jc; jout(mu) = T^{nu,mu} jc(nu) 
c
c     input:  double complex rank 2 tensor T(0:3,0:3)
c             double complex four vectors jc(0:3)
c     output: double complex four vectors jout(0:3)
c
c***************************************************************************
      implicit none
      complex*16 T(0:3,0:3), jc(0:3), jout(0:3)
      integer mu
      do mu = 0,3
         jout(mu) = T(0,mu)*jc(0) - T(1,mu)*jc(1) 
     &            - T(2,mu)*jc(2) - T(3,mu)*jc(3)
      enddo
      end


c***************************************************************************
      subroutine contract_T2j(T,jc,jout)
c***************************************************************************
c
c     contract second index of complex rank 2 tensor T^{mu,nu} with 
c     complex vector jc; jout(mu) = T^{mu,nu} jc(nu) 
c
c     input:  double complex rank 2 tensor T(0:3,0:3)
c             double complex four vectors jc(0:3)
c     output: double complex four vectors jout(0:3)
c
c***************************************************************************
      implicit none
      complex*16 T(0:3,0:3), jc(0:3), jout(0:3)
      integer mu
      do mu = 0,3
         jout(mu) = T(mu,0)*jc(0) - T(mu,1)*jc(1) 
     &            - T(mu,2)*jc(2) - T(mu,3)*jc(3)
      enddo
      end


c***************************************************************************
c **************************************************************************    
      subroutine contract_Tr1j(T,jc,jout)
c *************************************************************************     
      implicit none
c contract first index of real rank 2 tensor T^{mu,nu} with complex vector jc
c jout(mu) = T^{nu,mu} jc(nu) 
      complex*16  jc(0:3), jout(0:3)
      double precision T(0:3,0:3)
      integer mu
      do mu = 0,3
         jout(mu) = T(0,mu)*jc(0) - T(1,mu)*jc(1) 
     &            - T(2,mu)*jc(2) - T(3,mu)*jc(3)
      enddo
      end
c **************************************************************************    
      

c **************************************************************************
      subroutine contract_T1jr(T,jr,jout)
c **************************************************************************    
      implicit none
c contract first index of complex rank 2 tensor T^{mu,nu} with real vector jr
c jout(mu) = T^{nu,mu} jc(nu) 
      complex*16 T(0:3,0:3), jout(0:3)
      real*8 jr(0:3)
      integer mu
      do mu = 0,3
         jout(mu) = T(0,mu)*jr(0) - T(1,mu)*jr(1) 
     &            - T(2,mu)*jr(2) - T(3,mu)*jr(3)
      enddo
      end
c **************************************************************************


c **************************************************************************
      subroutine contract_Tr2j(T,jc,jout)
c **************************************************************************    
      implicit none
c contract second index of real rank 2 tensor T^{mu,nu} with complex vector jc
c jout(mu) = T^{mu,nu} jc(nu) 
      complex*16  jc(0:3), jout(0:3)
      double precision T(0:3,0:3)
      integer mu
      do mu = 0,3
         jout(mu) = T(mu,0)*jc(0) - T(mu,1)*jc(1) 
     &            - T(mu,2)*jc(2) - T(mu,3)*jc(3)
      enddo
      end
c ************************************************************************** 
   
      
c **************************************************************************
      subroutine contract_T2jr(T,jr,jout)
c **************************************************************************    
      implicit none
c contract second index of complex rank 2 tensor T^{mu,nu} with real vector jr
c jout(mu) = T^{mu,nu} jr(nu) 
      complex*16 T(0:3,0:3), jout(0:3)
      real*8 jr(0:3)
      integer mu
      do mu = 0,3
         jout(mu) = T(mu,0)*jr(0) - T(mu,1)*jr(1) 
     &            - T(mu,2)*jr(2) - T(mu,3)*jr(3)
      enddo
      end
c **************************************************************************


c **************************************************************************    
      subroutine contract_Tr3j(T,jc,i,Tout)
c **************************************************************************
c contract i-th index of real rank 3 tensor T^(mu,nu,rho) with complex vector jc
      implicit none
      double precision T(0:3,0:3,0:3)
      double complex jc(0:3),Tout(0:3,0:3)
      integer i,mu,nu
      if (i.eq.1) then
         do mu=0,3
            do nu=0,3
               Tout(mu,nu)= T(0,mu,nu)*jc(0)-T(1,mu,nu)*jc(1)-
     -              T(2,mu,nu)*jc(2)-T(3,mu,nu)*jc(3)
            enddo
         enddo
      endif
      if(i.eq.2) then
         do mu=0,3
            do nu=0,3
               Tout(mu,nu)=T(mu,0,nu)*jc(0)-T(mu,1,nu)*jc(1)-
     -              T(mu,2,nu)*jc(2)-T(mu,3,nu)*jc(3)
            enddo
         enddo
      endif
      if(i.eq.3) then
         do mu=0,3
            do nu=0,3
               Tout(mu,nu)=T(mu,nu,0)*jc(0)-T(mu,nu,1)*jc(1)-
     -              T(mu,nu,2)*jc(2)-T(mu,nu,3)*jc(3)
            enddo
         enddo
      endif
      end
c **************************************************************************
      

c **************************************************************************    
      subroutine contract_Tr3jr(T,jr,i,Tout)
c **************************************************************************    
c contract i-th index of real rank 3 tensor T^(mu,nu,rho) with real vector jr
      implicit none
      double precision T(0:3,0:3,0:3),jr(0:3),Tout(0:3,0:3)
      integer i,mu,nu
      if (i.eq.1) then
         do mu=0,3
            do nu=0,3
               Tout(mu,nu)= T(0,mu,nu)*jr(0)-T(1,mu,nu)*jr(1)-
     -              T(2,mu,nu)*jr(2)-T(3,mu,nu)*jr(3)
            enddo
         enddo
      endif
      if(i.eq.2) then
         do mu=0,3
            do nu=0,3
               Tout(mu,nu)=T(mu,0,nu)*jr(0)-T(mu,1,nu)*jr(1)-
     -              T(mu,2,nu)*jr(2)-T(mu,3,nu)*jr(3)
            enddo
         enddo
      endif
      if(i.eq.3) then
         do mu=0,3
            do nu=0,3
               Tout(mu,nu)=T(mu,nu,0)*jr(0)-T(mu,nu,1)*jr(1)-
     -              T(mu,nu,2)*jr(2)-T(mu,nu,3)*jr(3)
            enddo
         enddo
      endif
      end
c **************************************************************************

      
c **************************************************************************    
      subroutine contract_T3j(T,jc,i,Tout)
c **************************************************************************    
c contract i-th index of complex rank 3 tensor T^(mu,nu,rho) with complex vector jc      
      implicit none
      double complex T(0:3,0:3,0:3),jc(0:3),Tout(0:3,0:3)
      integer i,mu,nu
      if (i.eq.1) then
         do mu=0,3
            do nu=0,3
               Tout(mu,nu)= T(0,mu,nu)*jc(0)-T(1,mu,nu)*jc(1)-
     -              T(2,mu,nu)*jc(2)-T(3,mu,nu)*jc(3)
            enddo
         enddo
      endif
      if(i.eq.2) then
         do mu=0,3
            do nu=0,3
               Tout(mu,nu)=T(mu,0,nu)*jc(0)-T(mu,1,nu)*jc(1)-
     -              T(mu,2,nu)*jc(2)-T(mu,3,nu)*jc(3)
            enddo
         enddo
      endif
      if(i.eq.3) then
         do mu=0,3
            do nu=0,3
               Tout(mu,nu)=T(mu,nu,0)*jc(0)-T(mu,nu,1)*jc(1)-
     -              T(mu,nu,2)*jc(2)-T(mu,nu,3)*jc(3)
            enddo
         enddo
      endif
      end
c **************************************************************************


c **************************************************************************    
      subroutine contract_T3jr(T,jr,i,Tout)
c **************************************************************************    
c contract i-th index of complex rank 3 tensor T^(mu,nu,rho) with complex vector jc      
      implicit none
      double complex T(0:3,0:3,0:3),Tout(0:3,0:3)
      double precision jr(0:3)
      integer i,mu,nu
      if (i.eq.1) then
         do mu=0,3
            do nu=0,3
               Tout(mu,nu)= T(0,mu,nu)*jr(0)-T(1,mu,nu)*jr(1)-
     -              T(2,mu,nu)*jr(2)-T(3,mu,nu)*jr(3)
            enddo
         enddo
      endif
      if(i.eq.2) then
         do mu=0,3
            do nu=0,3
               Tout(mu,nu)=T(mu,0,nu)*jr(0)-T(mu,1,nu)*jr(1)-
     -              T(mu,2,nu)*jr(2)-T(mu,3,nu)*jr(3)
            enddo
         enddo
      endif
      if(i.eq.3) then
         do mu=0,3
            do nu=0,3
               Tout(mu,nu)=T(mu,nu,0)*jr(0)-T(mu,nu,1)*jr(1)-
     -              T(mu,nu,2)*jr(2)-T(mu,nu,3)*jr(3)
            enddo
         enddo
      endif
      end
c **************************************************************************
      
      
c **************************************************************************    
      subroutine contract_Tr4j(T,jc,i,Tout)
c **************************************************************************    
c contracts i-th index of real rank 4 tensor T^(mu,nu,rho,sigma) with complex vector jc
      implicit none
      double precision T(0:3,0:3,0:3,0:3)
      double complex jc(0:3),Tout(0:3,0:3,0:3)
      integer i,mu,nu,rho
      if ( i.eq. 1) then
         do mu=0,3
            do nu=0,3
               do rho=0,3
                  Tout(mu,nu,rho)=T(0,mu,nu,rho)*jc(0)-T(1,mu,nu,rho)*jc(1)
     &                 -T(2,mu,nu,rho)*jc(2)-T(3,mu,nu,rho)*jc(3)
               enddo
            enddo
         enddo
      endif
      if ( i.eq.2) then
         do mu=0,3
            do nu=0,3
               do rho=0,3
                  Tout(mu,nu,rho)=T(mu,0,nu,rho)*jc(0)-T(mu,1,nu,rho)*jc(1)
     &                 -T(mu,2,nu,rho)*jc(2)-T(mu,3,nu,rho)*jc(3)
               enddo
            enddo
         enddo
      endif
      if (i.eq.3) then
         do mu=0,3
            do nu=0,3
               do rho=0,3
                  Tout(mu,nu,rho)=T(mu,nu,0,rho)*jc(0)-T(mu,nu,1,rho)*jc(1)
     &                 -T(mu,nu,2,rho)*jc(2)-T(mu,nu,3,rho)*jc(3)
               enddo
            enddo
         enddo
      endif
      if (i.eq. 4) then
         do mu=0,3
            do nu=0,3
               do rho=0,3
                  Tout(mu,nu,rho)=T(mu,nu,rho,0)*jc(0)-T(mu,nu,rho,1)*jc(1)
     &                 -T(mu,nu,rho,2)*jc(2)-T(mu,nu,rho,3)*jc(3)
               enddo
            enddo
         enddo
      endif
      end
c **************************************************************************    
      

c **************************************************************************
      subroutine contract_T4j(T,jc,i,Tout)
c **************************************************************************   
c contracts i-th index of complex rank 4 tensor T^(mu,nu,rho,sigma) with complex vector jc
      implicit none
      double complex T(0:3,0:3,0:3,0:3)
      double complex jc(0:3),Tout(0:3,0:3,0:3)
      integer i,mu,nu,rho
      if ( i.eq. 1) then
         do mu=0,3
            do nu=0,3
               do rho=0,3
                  Tout(mu,nu,rho)=T(0,mu,nu,rho)*jc(0)-T(1,mu,nu,rho)*jc(1)
     &                 -T(2,mu,nu,rho)*jc(2)-T(3,mu,nu,rho)*jc(3)
               enddo
            enddo
         enddo
      endif
      if ( i.eq.2) then
         do mu=0,3
            do nu=0,3
               do rho=0,3
                  Tout(mu,nu,rho)=T(mu,0,nu,rho)*jc(0)-T(mu,1,nu,rho)*jc(1)
     &                 -T(mu,2,nu,rho)*jc(2)-T(mu,3,nu,rho)*jc(3)
               enddo
            enddo
         enddo
      endif
      if (i.eq.3) then
         do mu=0,3
            do nu=0,3
               do rho=0,3
                  Tout(mu,nu,rho)=T(mu,nu,0,rho)*jc(0)-T(mu,nu,1,rho)*jc(1)
     &                 -T(mu,nu,2,rho)*jc(2)-T(mu,nu,3,rho)*jc(3)
               enddo
            enddo
         enddo
      endif
      if (i.eq. 4) then
         do mu=0,3
            do nu=0,3
               do rho=0,3
                  Tout(mu,nu,rho)=T(mu,nu,rho,0)*jc(0)-T(mu,nu,rho,1)*jc(1)
     &                 -T(mu,nu,rho,2)*jc(2)-T(mu,nu,rho,3)*jc(3)
               enddo
            enddo
         enddo
      endif
      end
c **************************************************************************

      
c **************************************************************************    
      subroutine contract_Tr4jr(T,jr,i,Tout)
c **************************************************************************
c contracts i-th index of real rank 4 tensor T^(mu,nu,rho,sigma) with real vector jr
      implicit none
      double precision T(0:3,0:3,0:3,0:3),jr(0:3),Tout(0:3,0:3,0:3)      
      integer i,mu,nu,rho
      if ( i.eq. 1) then
         do mu=0,3
            do nu=0,3
               do rho=0,3
                  Tout(mu,nu,rho)=T(0,mu,nu,rho)*jr(0)-T(1,mu,nu,rho)*jr(1)
     &                 -T(2,mu,nu,rho)*jr(2)-T(3,mu,nu,rho)*jr(3)
               enddo
            enddo
         enddo
      endif
      if ( i.eq.2) then
         do mu=0,3
            do nu=0,3
               do rho=0,3
                  Tout(mu,nu,rho)=T(mu,0,nu,rho)*jr(0)-T(mu,1,nu,rho)*jr(1)
     &                 -T(mu,2,nu,rho)*jr(2)-T(mu,3,nu,rho)*jr(3)
               enddo
            enddo
         enddo
      endif
      if (i.eq.3) then
         do mu=0,3
            do nu=0,3
               do rho=0,3
                  Tout(mu,nu,rho)=T(mu,nu,0,rho)*jr(0)-T(mu,nu,1,rho)*jr(1)
     &                 -T(mu,nu,2,rho)*jr(2)-T(mu,nu,3,rho)*jr(3)
               enddo
            enddo
         enddo
      endif
      if (i.eq. 4) then
         do mu=0,3
            do nu=0,3
               do rho=0,3
                  Tout(mu,nu,rho)=T(mu,nu,rho,0)*jr(0)-T(mu,nu,rho,1)*jr(1)
     &                 -T(mu,nu,rho,2)*jr(2)-T(mu,nu,rho,3)*jr(3)
               enddo
            enddo
         enddo
      endif
      end
c **************************************************************************

c***************************************************************************
      function EPSCRRR(a,b,c,d)
c****************************************************************************   
c EPSCRRR calculates the the value of a 4 dimensional totally antisymmetric tensor
c EPS-mu-nu-rho-sigma which is contracted with 1 complex current A(mu) and 3 real momenta
c B(nu),C(rho),D(sigma)      
cMK The convention is eps_abcd v^a v^b v^c v^d with eps_0123 = +1
cMK this differs compared to the conventions in FeynCalc and FormCalc, where eps^0123 = +1
      IMPLICIT NONE
      double complex a(0:3),epscrrr
      double precision b(0:3),c(0:3),d(0:3)
      epscrrr= a(0)*b(1)*c(2)*d(3)-a(0)*b(1)*c(3)*d(2)-a(0)*b(2)*c(1)*d(3)
     &        +a(0)*b(2)*c(3)*d(1)+a(0)*b(3)*c(1)*d(2)-a(0)*b(3)*c(2)*d(1)
     &        -a(1)*b(0)*c(2)*d(3)+a(1)*b(0)*c(3)*d(2)+a(1)*b(2)*c(0)*d(3)
     &        -a(1)*b(2)*c(3)*d(0)-a(1)*b(3)*c(0)*d(2)+a(1)*b(3)*c(2)*d(0)
     &        +a(2)*b(0)*c(1)*d(3)-a(2)*b(0)*c(3)*d(1)-a(2)*b(1)*c(0)*d(3)
     &        +a(2)*b(1)*c(3)*d(0)+a(2)*b(3)*c(0)*d(1)-a(2)*b(3)*c(1)*d(0)
     &        -a(3)*b(0)*c(1)*d(2)+a(3)*b(0)*c(2)*d(1)+a(3)*b(1)*c(0)*d(2)
     &        -a(3)*b(1)*c(2)*d(0)-a(3)*b(2)*c(0)*d(1)+a(3)*b(2)*c(1)*d(0)
      END
c***************************************************************************

      
c***************************************************************************
      function EPSCCRR(a,b,c,d)
c***************************************************************************    
c EPSCCRR calculates the the value of a 4 dimensional totally antisymmetric tensor
c EPS-mu-nu-rho-sigma which is contracted with 2 complex currents A(mu) and B(nu) and 2 real 
c momenta C(rho) and D(sigma)      
cMK The convention is eps_abcd v^a v^b v^c v^d with eps_0123 = +1
cMK this differs compared to the conventions in FeynCalc and FormCalc, where eps^0123 = +1
      IMPLICIT NONE
      double complex a(0:3),b(0:3),epsccrr
      double precision c(0:3),d(0:3)
      epsccrr= a(0)*b(1)*c(2)*d(3)-a(0)*b(1)*c(3)*d(2)-a(0)*b(2)*c(1)*d(3)
     &        +a(0)*b(2)*c(3)*d(1)+a(0)*b(3)*c(1)*d(2)-a(0)*b(3)*c(2)*d(1)
     &        -a(1)*b(0)*c(2)*d(3)+a(1)*b(0)*c(3)*d(2)+a(1)*b(2)*c(0)*d(3)
     &        -a(1)*b(2)*c(3)*d(0)-a(1)*b(3)*c(0)*d(2)+a(1)*b(3)*c(2)*d(0)
     &        +a(2)*b(0)*c(1)*d(3)-a(2)*b(0)*c(3)*d(1)-a(2)*b(1)*c(0)*d(3)
     &        +a(2)*b(1)*c(3)*d(0)+a(2)*b(3)*c(0)*d(1)-a(2)*b(3)*c(1)*d(0)
     &        -a(3)*b(0)*c(1)*d(2)+a(3)*b(0)*c(2)*d(1)+a(3)*b(1)*c(0)*d(2)
     &        -a(3)*b(1)*c(2)*d(0)-a(3)*b(2)*c(0)*d(1)+a(3)*b(2)*c(1)*d(0)
      END      
c***************************************************************************    


c***************************************************************************
      function EPSRRCC(a,b,c,d)
c***************************************************************************      
c EPSCCRR calculates the the value of a 4 dimensional totally antisymmetric tensor
c EPS-mu-nu-rho-sigma which is contracted with 2 complex currents A(mu) and B(nu) and 2 real c momenta C(rho) and D(sigma)      
cMK The convention is eps_abcd v^a v^b v^c v^d with eps_0123 = +1
cMK this differs compared to the conventions in FeynCalc and FormCalc, where eps^0123 = +1
      IMPLICIT NONE
      double complex c(0:3),d(0:3),epsrrcc
      double precision a(0:3),b(0:3)
      epsrrcc= a(0)*b(1)*c(2)*d(3)-a(0)*b(1)*c(3)*d(2)-a(0)*b(2)*c(1)*d(3)
     &        +a(0)*b(2)*c(3)*d(1)+a(0)*b(3)*c(1)*d(2)-a(0)*b(3)*c(2)*d(1)
     &        -a(1)*b(0)*c(2)*d(3)+a(1)*b(0)*c(3)*d(2)+a(1)*b(2)*c(0)*d(3)
     &        -a(1)*b(2)*c(3)*d(0)-a(1)*b(3)*c(0)*d(2)+a(1)*b(3)*c(2)*d(0)
     &        +a(2)*b(0)*c(1)*d(3)-a(2)*b(0)*c(3)*d(1)-a(2)*b(1)*c(0)*d(3)
     &        +a(2)*b(1)*c(3)*d(0)+a(2)*b(3)*c(0)*d(1)-a(2)*b(3)*c(1)*d(0)
     &        -a(3)*b(0)*c(1)*d(2)+a(3)*b(0)*c(2)*d(1)+a(3)*b(1)*c(0)*d(2)
     &        -a(3)*b(1)*c(2)*d(0)-a(3)*b(2)*c(0)*d(1)+a(3)*b(2)*c(1)*d(0)
      END      
c***************************************************************************      

c***************************************************************************
      function EPSCCCC(a,b,c,d)
c***************************************************************************      
c EPSCCRR calculates the the value of a 4 dimensional totally antisymmetric tensor
c EPS-mu-nu-rho-sigma which is contracted with 2 complex currents A(mu) and B(nu) and 2 real c momenta C(rho) and D(sigma)      
cMK The convention is eps_abcd v^a v^b v^c v^d with eps_0123 = +1
cMK this differs compared to the conventions in FeynCalc and FormCalc, where eps^0123 = +1
      IMPLICIT NONE
      double complex a(0:3),b(0:3),c(0:3),d(0:3),epscccc
      epscccc= a(0)*b(1)*c(2)*d(3)-a(0)*b(1)*c(3)*d(2)-a(0)*b(2)*c(1)*d(3)
     &        +a(0)*b(2)*c(3)*d(1)+a(0)*b(3)*c(1)*d(2)-a(0)*b(3)*c(2)*d(1)
     &        -a(1)*b(0)*c(2)*d(3)+a(1)*b(0)*c(3)*d(2)+a(1)*b(2)*c(0)*d(3)
     &        -a(1)*b(2)*c(3)*d(0)-a(1)*b(3)*c(0)*d(2)+a(1)*b(3)*c(2)*d(0)
     &        +a(2)*b(0)*c(1)*d(3)-a(2)*b(0)*c(3)*d(1)-a(2)*b(1)*c(0)*d(3)
     &        +a(2)*b(1)*c(3)*d(0)+a(2)*b(3)*c(0)*d(1)-a(2)*b(3)*c(1)*d(0)
     &        -a(3)*b(0)*c(1)*d(2)+a(3)*b(0)*c(2)*d(1)+a(3)*b(1)*c(0)*d(2)
     &        -a(3)*b(1)*c(2)*d(0)-a(3)*b(2)*c(0)*d(1)+a(3)*b(2)*c(1)*d(0)
      END      
c***************************************************************************

c***************************************************************************
      function EPSCCCR(a,b,c,d)
c***************************************************************************      
c EPSCCRR calculates the the value of a 4 dimensional totally antisymmetric tensor
c EPS-mu-nu-rho-sigma which is contracted with 2 complex currents A(mu) and B(nu) and 2 real c momenta C(rho) and D(sigma)      
cMK The convention is eps_abcd v^a v^b v^c v^d with eps_0123 = +1
cMK this differs compared to the conventions in FeynCalc and FormCalc, where eps^0123 = +1
      IMPLICIT NONE
      double complex a(0:3),c(0:3),b(0:3),epscccr
      double precision d(0:3)
      epscccr= a(0)*b(1)*c(2)*d(3)-a(0)*b(1)*c(3)*d(2)-a(0)*b(2)*c(1)*d(3)
     &        +a(0)*b(2)*c(3)*d(1)+a(0)*b(3)*c(1)*d(2)-a(0)*b(3)*c(2)*d(1)
     &        -a(1)*b(0)*c(2)*d(3)+a(1)*b(0)*c(3)*d(2)+a(1)*b(2)*c(0)*d(3)
     &        -a(1)*b(2)*c(3)*d(0)-a(1)*b(3)*c(0)*d(2)+a(1)*b(3)*c(2)*d(0)
     &        +a(2)*b(0)*c(1)*d(3)-a(2)*b(0)*c(3)*d(1)-a(2)*b(1)*c(0)*d(3)
     &        +a(2)*b(1)*c(3)*d(0)+a(2)*b(3)*c(0)*d(1)-a(2)*b(3)*c(1)*d(0)
     &        -a(3)*b(0)*c(1)*d(2)+a(3)*b(0)*c(2)*d(1)+a(3)*b(1)*c(0)*d(2)
     &        -a(3)*b(1)*c(2)*d(0)-a(3)*b(2)*c(0)*d(1)+a(3)*b(2)*c(1)*d(0)
      END      
c***************************************************************************      

c***************************************************************************
      function EPSRCCC(a,b,c,d)
c***************************************************************************      
c EPSCCRR calculates the the value of a 4 dimensional totally antisymmetric tensor
c EPS-mu-nu-rho-sigma which is contracted with 2 complex currents A(mu) and B(nu) and 2 real c momenta C(rho) and D(sigma)      
cMK The convention is eps_abcd v^a v^b v^c v^d with eps_0123 = +1
cMK this differs compared to the conventions in FeynCalc and FormCalc, where eps^0123 = +1
      IMPLICIT NONE
      double complex d(0:3),c(0:3),b(0:3),epsrccc
      double precision a(0:3)
      epsrccc= a(0)*b(1)*c(2)*d(3)-a(0)*b(1)*c(3)*d(2)-a(0)*b(2)*c(1)*d(3)
     &        +a(0)*b(2)*c(3)*d(1)+a(0)*b(3)*c(1)*d(2)-a(0)*b(3)*c(2)*d(1)
     &        -a(1)*b(0)*c(2)*d(3)+a(1)*b(0)*c(3)*d(2)+a(1)*b(2)*c(0)*d(3)
     &        -a(1)*b(2)*c(3)*d(0)-a(1)*b(3)*c(0)*d(2)+a(1)*b(3)*c(2)*d(0)
     &        +a(2)*b(0)*c(1)*d(3)-a(2)*b(0)*c(3)*d(1)-a(2)*b(1)*c(0)*d(3)
     &        +a(2)*b(1)*c(3)*d(0)+a(2)*b(3)*c(0)*d(1)-a(2)*b(3)*c(1)*d(0)
     &        -a(3)*b(0)*c(1)*d(2)+a(3)*b(0)*c(2)*d(1)+a(3)*b(1)*c(0)*d(2)
     &        -a(3)*b(1)*c(2)*d(0)-a(3)*b(2)*c(0)*d(1)+a(3)*b(2)*c(1)*d(0)
      END      
c***************************************************************************      

c***************************************************************************
      function EPSRRRR(a,b,c,d)
c***************************************************************************    
c EPSRRRR calculates the the value of a 4 dimensional totally antisymmetric tensor
c EPS-mu-nu-rho-sigma which is contracted with 4 real momenta
c A(mu),B(nu),C(rho),D(sigma)      
cMK The convention is eps_abcd v^a v^b v^c v^d with eps_0123 = +1
cMK this differs compared to the conventions in FeynCalc and FormCalc, where eps^0123 = +1
      IMPLICIT NONE      
      double precision a(0:3),b(0:3),c(0:3),d(0:3),epsrrrr
      epsrrrr= a(0)*b(1)*c(2)*d(3)-a(0)*b(1)*c(3)*d(2)-a(0)*b(2)*c(1)*d(3)
     &        +a(0)*b(2)*c(3)*d(1)+a(0)*b(3)*c(1)*d(2)-a(0)*b(3)*c(2)*d(1)
     &        -a(1)*b(0)*c(2)*d(3)+a(1)*b(0)*c(3)*d(2)+a(1)*b(2)*c(0)*d(3)
     &        -a(1)*b(2)*c(3)*d(0)-a(1)*b(3)*c(0)*d(2)+a(1)*b(3)*c(2)*d(0)
     &        +a(2)*b(0)*c(1)*d(3)-a(2)*b(0)*c(3)*d(1)-a(2)*b(1)*c(0)*d(3)
     &        +a(2)*b(1)*c(3)*d(0)+a(2)*b(3)*c(0)*d(1)-a(2)*b(3)*c(1)*d(0)
     &        -a(3)*b(0)*c(1)*d(2)+a(3)*b(0)*c(2)*d(1)+a(3)*b(1)*c(0)*d(2)
     &        -a(3)*b(1)*c(2)*d(0)-a(3)*b(2)*c(0)*d(1)+a(3)*b(2)*c(1)*d(0)
      END
c****************************************************************************   
      

c****************************************************************************
      function EPSCCRMU(a,b,c,mu)
c****************************************************************************      
c EPSCCRMU calculates the remainung 4 Vector when the 4 dimensional totally antisymmetric
c tensor EPS-mu-nu-rho-sigma is contracted with 3 vectors A(mu), B(nu), C(rho), where
c A and B are assumed to be complex to include the current and C is assumed to be real to
c include the momenta      
      IMPLICIT NONE
      integer mu
      double complex epsccrmu,a(0:3),b(0:3)
      double complex ab01,ab02,ab03,ab12,ab13,ab23
      double precision c(0:3)
      ab01= a(0)*b(1)-a(1)*b(0)
      ab02= a(0)*b(2)-a(2)*b(0)
      ab03= a(0)*b(3)-a(3)*b(0)
      ab12= a(1)*b(2)-a(2)*b(1)
      ab13= a(1)*b(3)-a(3)*b(1)
      ab23= a(2)*b(3)-a(3)*b(2)
      if (mu.eq.0) then
       epsccrmu= -c(1)*ab23 +c(2)*ab13-c(3)*ab12
      endif
      if (mu.eq.1) then
       epsccrmu= -c(0)*ab23 +c(2)*ab03 -c(3)*ab02
      endif
      if (mu.eq.2) then
       epsccrmu= c(0)*ab13 -c(1)*ab03 +c(3)*ab01
      endif
      if (mu.eq.3) then
       epsccrmu= -c(0)*ab12 +c(1)*ab02 -c(2)*ab01
      endif
      END
c***************************************************************************


c****************************************************************************
      function EPSCRRMU(a,b,c,mu)
c****************************************************************************   
c EPSCRRMU calculates the remainung 4 Vector when the 4 dimensional totally antisymmetric
c tensor EPS-mu-nu-rho-sigma is contracted with 3 vectors A(mu), B(nu), C(rho), where
c A is assumed to be complex to include the current and B and C are assumed to be real to
c include the momenta      
cMK The convention is eps_abcd v^a v^b v^c v^d with eps_0123 = +1
cMK this differs compared to the conventions in FeynCalc and FormCalc, where eps^0123 = +1
      IMPLICIT NONE
      integer mu
      double complex epscrrmu,a(0:3)
      double complex ab01,ab02,ab03,ab12,ab13,ab23
      double precision b(0:3),c(0:3)
      ab01= a(0)*b(1)-a(1)*b(0)
      ab02= a(0)*b(2)-a(2)*b(0)
      ab03= a(0)*b(3)-a(3)*b(0)
      ab12= a(1)*b(2)-a(2)*b(1)
      ab13= a(1)*b(3)-a(3)*b(1)
      ab23= a(2)*b(3)-a(3)*b(2)
      if (mu.eq.0) then
       epscrrmu= -c(1)*ab23 +c(2)*ab13-c(3)*ab12
      endif
      if (mu.eq.1) then
       epscrrmu= -c(0)*ab23 +c(2)*ab03 -c(3)*ab02
      endif
      if (mu.eq.2) then
       epscrrmu= c(0)*ab13 -c(1)*ab03 +c(3)*ab01
      endif
      if (mu.eq.3) then
       epscrrmu= -c(0)*ab12 +c(1)*ab02 -c(2)*ab01
      endif
      END
c***************************************************************************
 

c***************************************************************************
      function EPSRRRMU(a,b,c,mu)
c***************************************************************************
c EPSRRRMU calculates the same as epscrrmu in the case you have no complex current but 3
c real momenta to contract      
cMK The convention is eps_abcd v^a v^b v^c v^d with eps_0123 = +1
cMK this differs compared to the conventions in FeynCalc and FormCalc, where eps^0123 = +1
      IMPLICIT NONE
      integer mu
      double precision epsrrrmu
      double precision ab01,ab02,ab03,ab12,ab13,ab23
      double precision a(0:3),b(0:3),c(0:3)
      ab01= a(0)*b(1)-a(1)*b(0)
      ab02= a(0)*b(2)-a(2)*b(0)
      ab03= a(0)*b(3)-a(3)*b(0)
      ab12= a(1)*b(2)-a(2)*b(1)
      ab13= a(1)*b(3)-a(3)*b(1)
      ab23= a(2)*b(3)-a(3)*b(2)
      if (mu.eq.0) then
       epsrrrmu= -c(1)*ab23 +c(2)*ab13-c(3)*ab12
      endif
      if (mu.eq.1) then
       epsrrrmu= -c(0)*ab23 +c(2)*ab03 -c(3)*ab02
      endif
      if (mu.eq.2) then
       epsrrrmu= c(0)*ab13 -c(1)*ab03 +c(3)*ab01
      endif
      if (mu.eq.3) then
       epsrrrmu= -c(0)*ab12 +c(1)*ab02 -c(2)*ab01
      endif
      END
c***************************************************************************
      
      
c***************************************************************************
      function EPSCRMUNU(a,b,mu,nu)
c***************************************************************************
c EPSCRMUNU calculates the remaining 4x4 Matrix when the 4 dimensional totally antisymmetric
c tensor EPS-mu-nu-rho-sigma is contracted with 2 vectors A(mu), B(nu), where A is
c assumed to be complex to include the current and B is assumed to be real to include 
c the momenta.
c COMMENT: note that when you have a product of time- and spatial components, e.g. ab01
c you have to change the sign because the spatial component has a minus      
      IMPLICIT NONE
      integer mu,nu
      Double complex epscrmunu,a(0:3)
      double complex ab01,ab02,ab03,ab12,ab13,ab23
      double precision  b(0:3)      
      ab01= a(0)*b(1)-a(1)*b(0)
      ab02= a(0)*b(2)-a(2)*b(0)
      ab03= a(0)*b(3)-a(3)*b(0)
      ab12= a(1)*b(2)-a(2)*b(1)
      ab13= a(1)*b(3)-a(3)*b(1)
      ab23= a(2)*b(3)-a(3)*b(2)
      if (mu .eq. nu) then
       epscrmunu =0
      endif
      if (mu .eq.0 .and. nu.eq.1) then
       epscrmunu= -ab23
      endif
      if (mu .eq.0 .and. nu.eq.2) then
       epscrmunu= +ab13
      endif
      if (mu .eq.0 .and. nu.eq.3) then
       epscrmunu= -ab12
      endif
      if (mu .eq.1 .and. nu.eq.0) then
       epscrmunu= +ab23
      endif
      if (mu .eq.1 .and. nu.eq.2) then
       epscrmunu= +ab03
      endif
      if (mu .eq.1 .and. nu.eq.3) then
       epscrmunu= -ab02
      endif
      if (mu .eq.2 .and. nu.eq.0) then
       epscrmunu= -ab13
      endif
      if (mu .eq.2 .and. nu.eq.1) then
       epscrmunu= -ab03
      endif
      if (mu .eq.2 .and. nu.eq.3) then
       epscrmunu= +ab01
      endif
      if (mu .eq.3 .and. nu.eq.0) then
       epscrmunu= +ab12
      endif
      if (mu .eq.3 .and. nu.eq.1) then
       epscrmunu= +ab02
      endif
      if (mu .eq.3 .and. nu.eq.2) then
       epscrmunu= -ab01
      endif      
      END
c***********************************************************************
      

c***********************************************************************
      function EPSRRMUNU(a,b,mu,nu)
c***********************************************************************
c EPSRRMUNU calculates the remaining 4x4 Matrix when the 4 dimensional totally antisymmetric
c tensor EPS-mu-nu-rho-sigma is contracted with 2 vectors A(mu), B(nu), where A and
c B are assumed to be real to include the momenta.
c COMMENT: note that when you have a product of time- and spatial components, e.g. ab01
c you have to change the sign because the spatial component has a minus      
      IMPLICIT NONE
      integer mu,nu
      Double precision epsrrmunu,a(0:3),b(0:3)
      double precision ab01,ab02,ab03,ab12,ab13,ab23            
      ab01= a(0)*b(1)-a(1)*b(0)
      ab02= a(0)*b(2)-a(2)*b(0)
      ab03= a(0)*b(3)-a(3)*b(0)
      ab12= a(1)*b(2)-a(2)*b(1)
      ab13= a(1)*b(3)-a(3)*b(1)
      ab23= a(2)*b(3)-a(3)*b(2)
      if (mu .eq. nu) then
       epsrrmunu =0
      endif
      if (mu .eq.0 .and. nu.eq.1) then
       epsrrmunu= -ab23
      endif
      if (mu .eq.0 .and. nu.eq.2) then
       epsrrmunu= +ab13
      endif
      if (mu .eq.0 .and. nu.eq.3) then
        epsrrmunu= -ab12
      endif
      if (mu .eq.1 .and. nu.eq.0) then
       epsrrmunu= +ab23
      endif
      if (mu .eq.1 .and. nu.eq.2) then
       epsrrmunu= +ab03
      endif
      if (mu .eq.1 .and. nu.eq.3) then
       epsrrmunu= -ab02
      endif
      if (mu .eq.2 .and. nu.eq.0) then
       epsrrmunu= -ab13
      endif
      if (mu .eq.2 .and. nu.eq.1) then
       epsrrmunu= -ab03
      endif
      if (mu .eq.2 .and. nu.eq.3) then
       epsrrmunu= +ab01
      endif
      if (mu .eq.3 .and. nu.eq.0) then
       epsrrmunu= +ab12
      endif
      if (mu .eq.3 .and. nu.eq.1) then
       epsrrmunu= +ab02
      endif
      if (mu .eq.3 .and. nu.eq.2) then
       epsrrmunu= -ab01
      endif      
      END
c******************************************************************************

c******************************************************************************
      function epsrmunurho(a,mu,nu,rho)
c******************************************************************************
c epsrmunurho calculates remaining rank 3 Tensor, when a totally antisymmetric
c epsilon tensor eps-Mu-Nu-Rho-Sigma is contractes with a real vector a(      
      implicit none
      double precision a(0:3),epsrmunurho
      integer mu,nu,rho 
      epsrmunurho=0.d0 
      if (mu .eq. 0 .and. nu .eq. 1 .and. rho .eq. 2) then
       epsrmunurho= -a(3)
      endif 
      if (mu .eq. 0 .and. nu .eq. 1 .and. rho .eq. 3) then
       epsrmunurho =+a(2)
      endif
      if (mu .eq. 0 .and. nu .eq. 2 .and. rho .eq. 1) then
       epsrmunurho = +a(3)
      endif
      if (mu .eq. 0 .and. nu .eq. 2 .and. rho .eq. 3) then
       epsrmunurho = -a(1)
      endif
      if (mu .eq. 0 .and. nu .eq. 3 .and. rho .eq. 1) then
       epsrmunurho = -a(2)
      endif
      if (mu .eq. 0 .and. nu .eq. 3 .and. rho .eq. 2) then
       epsrmunurho = +a(1)
      endif
      if (mu .eq. 1 .and. nu .eq. 0 .and. rho .eq. 2) then
       epsrmunurho = +a(3)
      endif
      if (mu .eq. 1 .and. nu .eq. 0 .and. rho .eq. 3) then
       epsrmunurho = -a(2)
      endif
      if (mu .eq. 1 .and. nu .eq. 2 .and. rho .eq. 0) then
       epsrmunurho = -a(3)
      endif
      if (mu .eq. 1 .and. nu .eq. 2 .and. rho .eq. 3) then
       epsrmunurho = -a(0)
      endif
      if (mu .eq. 1 .and. nu .eq. 3 .and. rho .eq. 0) then
       epsrmunurho = +a(2)
      endif
      if (mu .eq. 1 .and. nu .eq. 3 .and. rho .eq. 2) then
       epsrmunurho = +a(0)
      endif
      if (mu .eq. 2 .and. nu .eq. 0 .and. rho .eq. 1) then
       epsrmunurho = -a(3)
      endif
      if (mu .eq. 2 .and. nu .eq. 0 .and. rho .eq. 3) then
       epsrmunurho = +a(1)
      endif
      if (mu .eq. 2 .and. nu .eq. 1 .and. rho .eq. 0) then
       epsrmunurho = +a(3)
      endif
      if (mu .eq. 2 .and. nu .eq. 1 .and. rho .eq. 3) then
       epsrmunurho = +a(0)
      endif
      if (mu .eq. 2 .and. nu .eq. 3 .and. rho .eq. 0) then
       epsrmunurho = -a(1)
      endif
      if (mu .eq. 2 .and. nu .eq. 3 .and. rho .eq. 1) then
       epsrmunurho = -a(0)
      endif
      if (mu .eq. 3 .and. nu .eq. 0 .and. rho .eq. 1) then
       epsrmunurho = a(2)
      endif
      if (mu .eq. 3 .and. nu .eq. 0 .and. rho .eq. 2) then
       epsrmunurho = -a(1)
      endif
      if (mu .eq. 3 .and. nu .eq. 1 .and. rho .eq. 0) then
       epsrmunurho = -a(2)
      endif
      if (mu .eq. 3 .and. nu .eq. 1 .and. rho .eq. 2) then
       epsrmunurho = -a(0)
      endif
      if (mu .eq. 3 .and. nu .eq. 2 .and. rho .eq. 0) then
       epsrmunurho = +a(1)
      endif
      if (mu .eq. 3 .and. nu .eq. 2 .and. rho .eq. 1) then
       epsrmunurho = +a(0)
      endif
      end
c******************************************************************************
      
c******************************************************************************
      function epsmunurhosigma(mu,nu,rho,sigma)
c******************************************************************************
c epsmunurhosigma calculates the totally antisymmetric epsilon tensor EPS-mu-nu-rho-sigma
c when no index is contracted      
      IMPLICIT NONE
      integer mu,nu,rho,sigma
      real*8 epsmunurhosigma
      epsmunurhosigma=0.d0
      if (mu .eq. 0 .and. nu .eq. 1 .and. rho .eq. 2 .and. sigma .eq. 3) then
      epsmunurhosigma= -1.
      endif
      if (mu .eq. 0 .and. nu .eq. 1 .and. rho .eq. 3 .and. sigma .eq. 2) then
      epsmunurhosigma= 1.
      endif
      if (mu .eq. 0 .and. nu .eq. 2 .and. rho .eq. 1 .and. sigma .eq. 3) then
      epsmunurhosigma= 1.
      endif
      if (mu .eq. 0 .and. nu .eq. 2 .and. rho .eq. 3 .and. sigma .eq. 1) then
      epsmunurhosigma= -1.
      endif
      if (mu .eq. 0 .and. nu .eq. 3 .and. rho .eq. 1 .and. sigma .eq. 2) then
      epsmunurhosigma= -1.
      endif
      if (mu .eq. 0 .and. nu .eq. 3 .and. rho .eq. 2 .and. sigma .eq. 1) then
      epsmunurhosigma= 1.
      endif
      if (mu .eq. 1 .and. nu .eq. 0 .and. rho .eq. 2 .and. sigma .eq. 3) then
      epsmunurhosigma= 1.
      endif
      if (mu .eq. 1 .and. nu .eq. 0 .and. rho .eq. 3 .and. sigma .eq. 2) then
      epsmunurhosigma= -1.
      endif
      if (mu .eq. 1 .and. nu .eq. 2 .and. rho .eq. 0 .and. sigma .eq. 3) then
      epsmunurhosigma= -1.
      endif
      if (mu .eq. 1 .and. nu .eq. 2 .and. rho .eq. 3 .and. sigma .eq. 0) then
      epsmunurhosigma= 1.
      endif
      if (mu .eq. 1 .and. nu .eq. 3 .and. rho .eq. 0 .and. sigma .eq. 2) then
      epsmunurhosigma= 1.
      endif
      if (mu .eq. 1 .and. nu .eq. 3 .and. rho .eq. 2 .and. sigma .eq. 0) then
      epsmunurhosigma= -1.
      endif
      if (mu .eq. 2 .and. nu .eq. 0 .and. rho .eq. 1 .and. sigma .eq. 3) then
      epsmunurhosigma= -1.
      endif
      if (mu .eq. 2 .and. nu .eq. 0 .and. rho .eq. 3 .and. sigma .eq. 1) then
      epsmunurhosigma= 1.
      endif
      if (mu .eq. 2 .and. nu .eq. 1 .and. rho .eq. 0 .and. sigma .eq. 3) then
      epsmunurhosigma= 1.
      endif
      if (mu .eq. 2 .and. nu .eq. 1 .and. rho .eq. 3 .and. sigma .eq. 0) then
      epsmunurhosigma= -1.
      endif
      if (mu .eq. 2 .and. nu .eq. 3 .and. rho .eq. 0 .and. sigma .eq. 1) then
      epsmunurhosigma= -1.
      endif
      if (mu .eq. 2 .and. nu .eq. 3 .and. rho .eq. 1 .and. sigma .eq. 0) then
      epsmunurhosigma= 1.
      endif
      if (mu .eq. 3 .and. nu .eq. 0 .and. rho .eq. 1 .and. sigma .eq. 2) then
      epsmunurhosigma= 1.
      endif
      if (mu .eq. 3 .and. nu .eq. 0 .and. rho .eq. 2 .and. sigma .eq. 1) then
      epsmunurhosigma= -1.
      endif
      if (mu .eq. 3 .and. nu .eq. 1 .and. rho .eq. 0 .and. sigma .eq. 2) then
      epsmunurhosigma= -1.
      endif
      if (mu .eq. 3 .and. nu .eq. 1 .and. rho .eq. 2 .and. sigma .eq. 0) then
      epsmunurhosigma= 1.
      endif
      if (mu .eq. 3 .and. nu .eq. 2 .and. rho .eq. 0 .and. sigma .eq. 1) then
      epsmunurhosigma= 1.
      endif
      if (mu .eq. 3 .and. nu .eq. 2 .and. rho .eq. 1 .and. sigma .eq. 0) then
      epsmunurhosigma= -1.
      endif                  
      end      


c***************************************************************************
      double complex function contract_Trjj(T,j1,j2)
c***************************************************************************
c     
c     contract complex rank 2 tensor T^{mu,nu} with 
c     two complex vectors j1 and j2
c
c     input: double complex rank 2 tensor T(0:3,0:3)
c            double complex four vectors j1(0:3),j2(0:3)
c
c***************************************************************************
      implicit none
      double complex j1(0:3), j2(0:3), resv(0:3)
      double precision T(0:3,0:3)
      integer mu
      do mu = 0,3
         resv(mu) = T(mu,0)*j2(0) - T(mu,1)*j2(1) 
     &            - T(mu,2)*j2(2) - T(mu,3)*j2(3)
      enddo
      contract_Trjj = resv(0)*j1(0) - resv(1)*j1(1) 
     &             - resv(2)*j1(2) - resv(3)*j1(3)
      return
      end

c***************************************************************************
      double precision function contract_Trjcj(T,j1)
c***************************************************************************
c     
c     contract double precision rank 2 tensor T^{mu,nu} with 
c     two complex vectors j1 and j1*
c
c     input: double pecision rank 2 tensor T(0:3,0:3)
c            double complex four vectors j1(0:3)
c
c***************************************************************************

      implicit none
      double complex j1(0:3), j2(0:3), resv(0:3)
      double precision T(0:3,0:3)
      integer mu
      do mu=0,3
      j2(mu)=dconjg(j1(mu))
      enddo
      do mu = 0,3
         resv(mu) = T(mu,0)*j2(0) - T(mu,1)*j2(1) 
     &            - T(mu,2)*j2(2) - T(mu,3)*j2(3)
      enddo
      contract_Trjcj = dreal(resv(0)*j1(0) - resv(1)*j1(1) 
     &             - resv(2)*j1(2) - resv(3)*j1(3))
      return
      end

c***************************************************************************
      double precision function contract_grjcj(j1)
c***************************************************************************
c     
c     contract metric g^{mu,nu} with 
c     two complex vectors j1 and j1*
c
c     input: double complex four vectors j1(0:3)
c
c***************************************************************************

      implicit none
      double complex j1(0:3)
      integer mu

      contract_grjcj = j1(0)*dconjg(j1(0))
      do mu = 1,3
         contract_grjcj = contract_grjcj - j1(mu)*dconjg(j1(mu))
      enddo
      return
      end

