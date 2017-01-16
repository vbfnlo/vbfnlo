       subroutine tens_red3_new_Re_Com_R4_G(m0,m1,m2,p1sq,p2sq,s12,B0R_23,B0R_13,B0R_12, 
     & BijR_23,BijR_13,BijR_12,B0I_23,B0I_13,B0I_12,BijI_23,BijI_13,BijI_12,
     & C0,C0r,C0I,Cijr,CijI)
C                tens_red3 = 3-point tensors
C     Francisco Campanario
C    Email: Francam@particle.uni-karlsruhe.de
C Date: 25/02/2010
c
c  determine the Passarino-Veltman tensor decomposition for the three-point
c  tensor integrals
c
c                       d^4k k_mu1 k_mu2 k_mu3,...,k_mu1 k_mu2 k_mu3 k_mu4 k_mu5
c C0;C_mu;C_munu =Int -------  -----------------------------------------
c                     (2pi)^4    [k^2-m^2][(k+p1)^2-m^2][(k+p1+p2)^2-m^2] 
c
c with
c
c    C_mu = p1_mu C11  +  p_2_mu C12
c
c  C_munu = p1_mu p1_nu C21 + p2_mu p2_nu C22 + 
c           (p1_mu p2_nu + p1_nu p2_mu) C23  -  g_munu C24
c
c  for notation see Passarino&Veltman, NP B160 (1979) 151 and my notes
c
C INPUT:  p1sq, p2sq, s12          external invariants: p1^2, p2^2, s12=(p1+p2)^2
C         B0, C0                   4 scalar integrals; the 3 B0 are, 
c                                  in PV notation:
c         B0(1) = B0(1,2) = B0(p1)  B_0 function with subtraction of 
c         B0(2) = B0(2,3) = B0(p2)  divergent term
c         B0(3) = B0(1,3) = B0(s12)
c
c OUTPUT: Cij(n,m) = C_nm          form factors in the tensor integrals
c          n=1,2,3,4; n=1,2        a la PV
c

      implicit none
      real * 8  p1sq, p2sq, s12
      complex*16 C0,ct
      real*8 r1, r2r1, det,p1p2,r10,r21

      real*8 BijR_12(2,3),BijR_13(2,3),BijR_23(2,3)
      real*8 BijI_12(2,3),BijI_13(2,3),BijI_23(2,3)

      real*8 B1r_12,B1r_13,B1r_23,Cijr(9,4)
      real*8 B1I_12,B1I_13,B1I_23,CijI(9,4)

      real*8 B0r_23, B0r_13, B0r_12, C0r
      real*8 B0I_23, B0I_13, B0I_12, C0I 

      real*8 z11,z12,z21,z22,iz11,iz22

      real*8 Rr(2),RI(2),PRr(2),PRI(2)

      real*8 B21r_12,B21r_13,B21r_23
      real*8 B21I_12,B21I_13,B21I_23
      real*8 B22r_12,B22r_13,B22r_23
      real*8 B22I_12,B22I_13,B22I_23

      real*8 m0,m0sq,m1,m2,m1sq,m2sq
      real*8 deter,detAbs
      integer i,j
      Logical Singular,singular1
      common/DetCountML/Singular
c New for Gram
      common/GramSwitch/GramCS,GramDS
      logical GramCS,GRAMDS
      Integer jjinit

      real*8 tempjj,tempkl,IX,Zmax
      Common/Decide/tempjj,tempkl,IX,ZMax,jjinit
      Save/Decide/
c New variables to compare        
      real*8 musq
      Common/musqInv/musq

      complex*16 A0finG
      external A0finG
      real*8 accuracyC(0:4,5),AccuracyD(0:5,4)

      Common/Accuracy/AccuracyC,AccuracyD
      real*8 ratio
       logical printmy
       common/mprint/printmy
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ratio=1d-6 ! use at the end

      singular1=.false.
      p1p2 = (s12 - p1sq - p2sq)*0.5d0

      m0sq=m0*m0
      m1sq=m1*m1
      m2sq=m2*m2

       r10 = p1sq-m1sq+m0sq
       r21 = s12-m2sq -(p1sq-m1sq)
       r1 = r10
       r2r1 = r21


       deter = abs(2.d0*(p1sq*p2sq - p1p2*p1p2))
       detAbs=  abs(2.d0*(abs(p1sq*p2sq)+abs(p1p2*p1p2)))
      
c       If( (deter/detAbs).le.1d-5) Singular=.true.

       If( abs(deter)/abs(detAbs).lt.4.6d-3) then
         Singular=.true.
         singular1=.true.
c         print*, "Here1",abs(deter)/abs(detAbs),Singular, GramCS
         do i=1,4 
            do j=1,9
               Cijr(j,i)=0d0
               CijI(j,i)=0d0
               C0r=0d0
               C0I=0d0
            enddo
         enddo
c         stop
c         return
         endif

 
       C0r=Dble(C0)
       C0I=DImag(C0)


       B1r_12 = BijR_12(1,1)
       B1r_13 = BijR_13(1,1)  
       B1r_23 = BijR_23(1,1) 
       
       B1I_12 = BijI_12(1,1) 
       B1I_13 = BijI_13(1,1) 
       B1I_23 = BijI_23(1,1) 
   
       B21r_12= BijR_12(1,2)
       B22r_12= BijR_12(2,2)
       B21I_12= BijI_12(1,2)
       B22I_12= BijI_12(2,2)

       B21r_13= BijR_13(1,2)
       B22r_13= BijR_13(2,2)
       B21I_13= BijI_13(1,2)
       B22I_13= BijI_13(2,2)


       B21r_23= BijR_23(1,2)
       B22r_23= BijR_23(2,2)
       B21I_23= BijI_23(1,2)
       B22I_23= BijI_23(2,2)


      If(abs(p1sq).gt.abs(p1p2)) then
          z11=2d0*p1sq
          iz11=1d0/z11
          z12=2d0*p1p2 
          z21=z12*iz11
          z22=2d0*p2sq-z12*z21
          iz22=1d0/z22
c          iorder(1)=1
c          iorder(2)=2
          det=z11*z22 

c 1-2
       PRr(1) = (B0r_13 - B0r_23 - C0r*r1)
       PRr(2) = (B0r_12 - B0r_13 - C0r*r2r1)
   
       Rr(2)=(PRr(2)-z21*PRr(1))*iz22 
       Rr(1)=(PRr(1)-z12*Rr(2))*iz11 

       Cijr(1,1) = Rr(1)
       Cijr(2,1) = Rr(2)

       PRI(1) = (B0I_13 - B0I_23 - C0I*r1)
       PRI(2) = (B0I_12 - B0I_13 - C0I*r2r1)
      
       RI(2)=(PRI(2)-z21*PRI(1))*iz22 
       RI(1)=(PRI(1)-z12*RI(2))*iz11 

       CijI(1,1) = RI(1)
       CijI(2,1) = RI(2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c C00
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Cijr(4,2) = ( B0r_23 +2d0*m0sq*C0r + Cijr(1,1)*r1 +Cijr(2,1)*r2r1 +1.d0)*0.25d0
      CijI(4,2) = ( B0I_23 +2d0*m0sq*C0I + CijI(1,1)*r1 +CijI(2,1)*r2r1)*0.25d0    
c 3-5
       PRr(1) = (B1r_13 + B0r_23 - Cijr(1,1)*r1 - Cijr(4,2)*2.d0)
       PRr(2) = (B1r_12 - B1r_13 - Cijr(1,1)*r2r1)

       Rr(2)=(PRr(2)-z21*PRr(1))*iz22 
       Rr(1)=(PRr(1)-z12*Rr(2))*iz11 

       Cijr(1,2) =Rr(1)
       Cijr(3,2) =Rr(2)

       PRI(1) = (B1I_13 + B0I_23 - CijI(1,1)*r1 - CijI(4,2)*2.d0)
       PRI(2) = (B1I_12 - B1I_13 - CijI(1,1)*r2r1)

       RI(2)=(PRI(2)-z21*PRI(1))*iz22 
       RI(1)=(PRI(1)-z12*RI(2))*iz11 

       CijI(1,2) =RI(1)
       CijI(3,2) =RI(2)

c 4-6
       PRr(1) = (  B1r_13 - B1r_23 - Cijr(2,1)*r1)
       PRr(2) = (- B1r_13          -Cijr(2,1)*r2r1 -Cijr(4,2)*2.d0)
      
       Rr(2)=(PRr(2)-z21*PRr(1))*iz22 
c      Rr(1)=(PRr(1)-z12*Rr(2))*iz11 

       Cijr(2,2) = Rr(2)
c      Cijr(3,2) = Rr(1)


       PRI(1) = (  B1I_13 - B1I_23 - CijI(2,1)*r1)
       PRI(2) = (- B1I_13          - CijI(2,1)*r2r1 -CijI(4,2)*2.d0)

       RI(2)=(PRI(2)-z21*PRI(1))*iz22 
c      RI(1)=(PRI(1)-z12*RI(2))*iz11 

       CijI(2,2) = RI(2)
c      CijI(3,2) = RI(1)


c  001,002
       PRr(1) = B22R_13 - B22r_23 - r1 *Cijr(4, 2)
       PRr(2) = B22r_12 - B22r_13 - r2r1 *Cijr(4, 2)

       Rr(2)=(PRr(2)-z21*PRr(1))*iz22 
       Rr(1)=(PRr(1)-z12*Rr(2))*iz11 

       Cijr(5,3)=Rr(1)
       Cijr(6,3)=Rr(2)


       PrI(1) = B22I_13 - B22I_23 - r1 *CijI(4, 2)
       PrI(2) = B22I_12 - B22I_13 - r2r1 *CijI(4, 2)

       RI(2)=(PrI(2)-z21*PrI(1))*iz22 
       RI(1)=(PrI(1)-z12*RI(2))*iz11 

       CijI(5,3)=RI(1)
       CijI(6,3)=RI(2)


c 111,211

       PRr(1) = -B0r_23 + B21r_13 - r1* Cijr(1,2) - 4d0* Cijr(5,3) 
       PRr(2) = B21r_12 - B21r_13 - r2r1* Cijr(1,2)
   
       Rr(2)=(PRr(2)-z21*PRr(1))*iz22 
       Rr(1)=(PRr(1)-z12*Rr(2))*iz11 

       Cijr(1,3) = Rr(1)
       Cijr(3,3) = Rr(2)

       PrI(1) = -B0I_23 + B21I_13 - r1* CijI(1,2) - 4d0* CijI(5,3) 
       PrI(2) = B21I_12 - B21I_13 - r2r1* CijI(1,2)
   
       RI(2)=(PrI(2)-z21*PrI(1))*iz22 
       RI(1)=(PrI(1)-z12*RI(2))*iz11 

       CijI(1,3) = RI(1)
       CijI(3,3) = RI(2)

c  121,221

       PRr(1)=B21r_13 + B1r_23 - r1* Cijr(3, 2) - 2d0* Cijr(6, 3)
       PRr(2)=-B21r_13 - r2r1* Cijr(3, 2) - 2d0* Cijr(5, 3)
       
       Rr(2)=(PRr(2)-z21*PRr(1))*iz22 
c       Rr(1)=(PRr(1)-z12*Rr(2))*iz11 

c       Cijr(1,1) = Rr(1)
       Cijr(4,3) = Rr(2)
******
       PRI(1)=B21I_13 + B1I_23 - r1* CijI(3, 2) - 2d0* CijI(6, 3)
       PRI(2)=-B21I_13 - r2r1* CijI(3, 2) - 2d0* CijI(5, 3)
       
       RI(2)=(PRI(2)-z21*PRI(1))*iz22 
c       RI(1)=(PRI(1)-z12*RI(2))*iz11 

c       CijI(1,1) = RI(1)
       CijI(4,3) = RI(2)

c  122,222
       PRr(1)=B21r_13 - B21r_23 - r1*Cijr(2, 2)
       PRr(2)=-B21r_13 - r2r1* Cijr(2, 2) - 4d0*Cijr(6, 3)

       Rr(2)=(PRr(2)-z21*PRr(1))*iz22 
c       Rr(1)=(PRr(1)-z12*Rr(2))*iz11 

c       Cijr(1,1) = Rr(1)
       Cijr(2,3) = Rr(2)


       PrI(1)=B21I_13 - B21I_23  - r1*CijI(2, 2)
       PrI(2)=-B21I_13 - r2r1* CijI(2, 2) - 4d0*CijI(6, 3)

       RI(2)=(PRI(2)-z21*PRI(1))*iz22 
c       RI(1)=(PRI(1)-z12*RI(2))*iz11 

c       CijI(1,1) = RI(1)
       CijI(2,3) = RI(2)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 0000

       Ct=(4*m0sq + 4*m1sq + 4*m2sq - p1sq - 
     &     p2sq - s12)/192.d0 + 
     &  (DCMPLX(B22R_23,B22I_23) 
     & +2*m0sq*DCMPLX(CijR(4,2),CijI(4,2)) + 
     &     r10*DCMPLX(CijR(5,3),CijI(5,3)) 
     &   + r21*DCMPLX(CijR(6,3),CijI(6,3))
     &     )/8.d0
     
       CijR(9,4)=Dble(ct)
       CijI(9,4)=DIMAG(ct)


c 0011,0021
       PrR(1)=BijR_13(2,3) + BijR_23(2,2) - r10*CijR(5,3) - 2*CijR(9,4)
       PrR(2)=BijR_12(2,3) - BijR_13(2,3) - r21*CijR(5,3) 

       Rr(2)=(PrR(2)-z21*PrR(1))*iz22 
       Rr(1)=(PrR(1)-z12*Rr(2))*iz11 
       CijR(6,4)=Rr(1)   
       CijR(8,4)=Rr(2)

       PrI(1)=BijI_13(2,3) + BijI_23(2,2) - r10*CijI(5,3) - 2*CijI(9,4)
       PrI(2)=BijI_12(2,3) - BijI_13(2,3) - r21*CijI(5,3) 
       RI(2)=(PrI(2)-z21*PrI(1))*iz22 
       RI(1)=(PrI(1)-z12*RI(2))*iz11 

       CijI(6,4)=RI(1)   
       CijI(8,4)=RI(2)


c0012,0022
       PrR(1)=BijR_13(2,3) - BijR_23(2,3) - r10*CijR(6,3)
       PrR(2)=-BijR_13(2,3) - r21*CijR(6,3) - 2*CijR(9,4)  
       Rr(2)=(PrR(2)-z21*PrR(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c      Cij(8,4)=Rr(1)
       CijR(7,4)=Rr(2)


       PrI(1)=BijI_13(2,3) - BijI_23(2,3) - r10*CijI(6,3)
       PrI(2)=-BijI_13(2,3) - r21*CijI(6,3) - 2*CijI(9,4)  
       RI(2)=(PrI(2)-z21*PrI(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c      Cij(8,4)=Rr(1)
       CijI(7,4)=RI(2)


       return
ccccccccccccccccc
cccccccccccccccc
ccccccccccccccccc
       else
		  z11=2d0*p1p2
		  iz11=1d0/z11
		  z21=2d0*p1sq*iz11
		  z12=2d0*p2sq
		  z22=z11-z12*z21
		  iz22=1d0/z22
c          iorder(1)=2
c          iorder(2)=1
          det=-z11*z22

c 1-2
       PRr(1) = (B0r_13 - B0r_23 - C0r*r1)
       PRr(2) = (B0r_12 - B0r_13 - C0r*r2r1)
   
       Rr(2)=(PRr(1)-z21*PRr(2))*iz22 
       Rr(1)=(PRr(2)-z12*Rr(2))*iz11 
 
       Cijr(1,1) = Rr(1)
       Cijr(2,1) = Rr(2)
 
       PRI(1) = (B0I_13 - B0I_23 - C0I*r1)
       PRI(2) = (B0I_12 - B0I_13 - C0I*r2r1)
      
       RI(2)=(PRI(1)-z21*PRI(2))*iz22 
       RI(1)=(PRI(2)-z12*RI(2))*iz11 
 
       CijI(1,1) = RI(1)
       CijI(2,1) = RI(2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c C00
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Cijr(4,2) = ( B0r_23 +2d0*m0sq*C0r + Cijr(1,1)*r1 +Cijr(2,1)*r2r1 +1.d0)*0.25d0
      CijI(4,2) = ( B0I_23 +2d0*m0sq*C0I + CijI(1,1)*r1 +CijI(2,1)*r2r1)*0.25d0
c 3-5
       PRr(1) = (B1r_13 + B0r_23 - Cijr(1,1)*r1 - Cijr(4,2)*2.d0)
       PRr(2) = (B1r_12 - B1r_13 - Cijr(1,1)*r2r1)
       
       Rr(2)=(PRr(1)-z21*PRr(2))*iz22 
       Rr(1)=(PRr(2)-z12*Rr(2))*iz11 
       
       Cijr(1,2) =Rr(1)
       Cijr(3,2) =Rr(2)
       
       PRI(1) = (B1I_13 + B0I_23 - CijI(1,1)*r1 - CijI(4,2)*2.d0)
       PRI(2) = (B1I_12 - B1I_13 - CijI(1,1)*r2r1)
       
       RI(2)=(PRI(1)-z21*PRI(2))*iz22 
       RI(1)=(PRI(2)-z12*RI(2))*iz11 
       
       CijI(1,2) =RI(1)
       CijI(3,2) =RI(2)

c 4-6
       PRr(1) = (  B1r_13 - B1r_23 - Cijr(2,1)*r1)
       PRr(2) = (- B1r_13          -Cijr(2,1)*r2r1 -Cijr(4,2)*2.d0)
       
       Rr(2)=(PRr(1)-z21*PRr(2))*iz22 
c      Rr(1)=(PRr(2)-z12*Rr(2))*iz11 

       Cijr(2,2) = Rr(2)
c      Cijr(3,2) = Rr(1)


       PRI(1) = (  B1I_13 - B1I_23 - CijI(2,1)*r1)
       PRI(2) = (- B1I_13          - CijI(2,1)*r2r1 -CijI(4,2)*2.d0)
       
       RI(2)=(PRI(1)-z21*PRI(2))*iz22 
c      RI(1)=(PRI(2)-z12*RI(2))*iz11 


       CijI(2,2) = RI(2)
c      CijI(3,2) = RI(1)


c  001,002
       PRr(1) = B22R_13 - B22r_23 - r1 *Cijr(4, 2)
       PRr(2) = B22r_12 - B22r_13 - r2r1 *Cijr(4, 2)

       Rr(2)=(PRr(1)-z21*PRr(2))*iz22 
       Rr(1)=(PRr(2)-z12*Rr(2))*iz11 

       Cijr(5,3)=Rr(1)
       Cijr(6,3)=Rr(2)


       PrI(1) = B22I_13 - B22I_23 - r1 *CijI(4, 2)
       PrI(2) = B22I_12 - B22I_13 - r2r1 *CijI(4, 2)

       RI(2)=(PrI(1)-z21*PrI(2))*iz22 
       RI(1)=(PrI(2)-z12*RI(2))*iz11 

       CijI(5,3)=RI(1)
       CijI(6,3)=RI(2)
c 111,211
       PRr(1) = -B0r_23 + B21r_13 - r1* Cijr(1,2) - 4d0* Cijr(5,3) 
       PRr(2) = B21r_12 - B21r_13 - r2r1* Cijr(1,2)
   
       Rr(2)=(PRr(1)-z21*PRr(2))*iz22 
       Rr(1)=(PRr(2)-z12*Rr(2))*iz11 

       Cijr(1,3) = Rr(1)
       Cijr(3,3) = Rr(2)

       PrI(1) = -B0I_23 + B21I_13 - r1* CijI(1,2) - 4d0* CijI(5,3) 
       PrI(2) = B21I_12 - B21I_13 - r2r1* CijI(1,2)
   
       RI(2)=(PrI(1)-z21*PrI(2))*iz22 
       RI(1)=(PrI(2)-z12*RI(2))*iz11 

       CijI(1,3) = RI(1)
       CijI(3,3) = RI(2)
c  121,221
       PRr(1)=B21r_13 + B1r_23 - r1* Cijr(3, 2) - 2d0* Cijr(6, 3)
       PRr(2)=-B21r_13 - r2r1* Cijr(3, 2) - 2d0* Cijr(5, 3)
       
       Rr(2)=(PRr(1)-z21*PRr(2))*iz22 
c       Rr(1)=(PRr(2)-z12*Rr(2))*iz11 

c       Cijr(1,1) = Rr(1)
       Cijr(4,3) = Rr(2)
******
       PRI(1)=B21I_13 + B1I_23 - r1* CijI(3, 2) - 2d0* CijI(6, 3)
       PRI(2)=-B21I_13 - r2r1* CijI(3, 2) - 2d0* CijI(5, 3)
       
       RI(2)=(PRI(1)-z21*PRI(2))*iz22 
c       RI(1)=(PRI(2)-z12*RI(2))*iz11 

c       CijI(1,1) = RI(1)
       CijI(4,3) = RI(2)

c  122,222
       PRr(1)=B21r_13 - B21r_23 - r1*Cijr(2, 2)
       PRr(2)=-B21r_13 - r2r1* Cijr(2, 2) - 4d0*Cijr(6, 3)

       Rr(2)=(PRr(1)-z21*PRr(2))*iz22 
c       Rr(1)=(PRr(2)-z12*Rr(2))*iz11 

c       Cijr(1,1) = Rr(1)
       Cijr(2,3) = Rr(2)


       PrI(1)=B21I_13 - B21I_23  - r1*CijI(2, 2)
       PrI(2)=-B21I_13 - r2r1* CijI(2, 2) - 4d0*CijI(6, 3)

       RI(2)=(PRI(1)-z21*PRI(2))*iz22 
c       RI(1)=(PRI(2)-z12*RI(2))*iz11 

c       CijI(1,1) = RI(1)
       CijI(2,3) = RI(2)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 0000
c       Cij(9,4)=(4*m0sq + 4*m1sq + 4*m2sq - p1sq - 
c     &     p2sq - s12)/192.d0 + 
c     &  (B23(2,2) + 2*m0sq*Cij(4,2) + 
c     &     r10*Cij(5,3) + r21*Cij(6,3)
c     &     )/8.d0

       Ct=(4*m0sq + 4*m1sq + 4*m2sq - p1sq - 
     &     p2sq - s12)/192.d0 + 
     &  (DCMPLX(B22R_23,B22I_23) 
     & +2*m0sq*DCMPLX(CijR(4,2),CijI(4,2)) + 
     &     r10*DCMPLX(CijR(5,3),CijI(5,3)) 
     &   + r21*DCMPLX(CijR(6,3),CijI(6,3))
     &     )/8.d0
     
       CijR(9,4)=Dble(ct)
       CijI(9,4)=DIMAG(ct)  

c 0011,0021
c       Pr(1)=B13(2,3) + B23(2,2) - r10*Cij(5,3) - 2*Cij(9,4)
c       Pr(2)=B12(2,3) - B13(2,3) - r21*Cij(5,3) 
c       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c       Cij(6,4)=Rr(1)   
c       Cij(8,4)=Rr(2)

       PrR(1)=BijR_13(2,3) + BijR_23(2,2) - r10*CijR(5,3) - 2*CijR(9,4)
       PrR(2)=BijR_12(2,3) - BijR_13(2,3) - r21*CijR(5,3) 

       Rr(2)=(PrR(1)-z21*PrR(2))*iz22 
       Rr(1)=(PrR(2)-z12*Rr(2))*iz11 
       CijR(6,4)=Rr(1)   
       CijR(8,4)=Rr(2)

       PrI(1)=BijI_13(2,3) + BijI_23(2,2) - r10*CijI(5,3) - 2*CijI(9,4)
       PrI(2)=BijI_12(2,3) - BijI_13(2,3) - r21*CijI(5,3) 
       RI(2)=(PrI(1)-z21*PrI(2))*iz22 
       RI(1)=(PrI(2)-z12*RI(2))*iz11 

       CijI(6,4)=RI(1)   
       CijI(8,4)=RI(2)
c0012,0022
c       Pr(1)=B13(2,3) - B23(2,3) - r10*Cij(6,3)
c       Pr(2)=-B13(2,3) - r21*Cij(6,3) - 2*Cij(9,4)  
c       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c      Cij(8,4)=Rr(1)
c       Cij(7,4)=Rr(2)

       PrR(1)=BijR_13(2,3) - BijR_23(2,3) - r10*CijR(6,3)
       PrR(2)=-BijR_13(2,3) - r21*CijR(6,3) - 2*CijR(9,4)  
       Rr(2)=(PrR(1)-z21*PrR(2))*iz22 
c       Rr(1)=(Pr(2)-z12*Rr(2))*iz11 
c      Cij(8,4)=Rr(1)
       CijR(7,4)=Rr(2)


       PrI(1)=BijI_13(2,3) - BijI_23(2,3) - r10*CijI(6,3)
       PrI(2)=-BijI_13(2,3) - r21*CijI(6,3) - 2*CijI(9,4)  
       RI(2)=(PrI(1)-z21*PrI(2))*iz22 
c       Rr(1)=(Pr(2)-z12*Rr(2))*iz11 
c      Cij(8,4)=Rr(1)
       CijI(7,4)=RI(2)


       endif


         singular=.false.
         If( abs(deter)/abs(detAbs).lt.ratio) then
         Singular=.true.
c         singular1=.true.
c         print*, "Here2HereCCRank4",abs(deter)/abs(detAbs),Singular, GramCS
         do i=1,4 
            do j=1,9
               Cijr(j,i)=0d0
               CijI(j,i)=0d0
               C0r=0d0
               C0I=0d0
            enddo
         enddo
c         stop
         return
         endif

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
