       subroutine tens_red3_new_Re_Com_R4_G_QUAD(m0,m1,m2,p1sq,p2sq,s12,B0R_23,B0R_13,B0R_12, 
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
      real*16  p1sq, p2sq, s12
      complex*32 C0,ct
      real*16 r1, r2r1, det,p1p2,r10,r21

      real*16 BijR_12(2,3),BijR_13(2,3),BijR_23(2,3)
      real*16 BijI_12(2,3),BijI_13(2,3),BijI_23(2,3)

      real*16 B1r_12,B1r_13,B1r_23,Cijr(9,4)
      real*16 B1I_12,B1I_13,B1I_23,CijI(9,4)

      real*16 B0r_23, B0r_13, B0r_12, C0r
      real*16 B0I_23, B0I_13, B0I_12, C0I 

      real*16 z11,z12,z21,z22,iz11,iz22

      real*16 Rr(2),RI(2),PRr(2),PRI(2)

      real*16 B21r_12,B21r_13,B21r_23
      real*16 B21I_12,B21I_13,B21I_23
      real*16 B22r_12,B22r_13,B22r_23
      real*16 B22I_12,B22I_13,B22I_23

      real*16 m0,m0sq,m1,m2,m1sq,m2sq
      real*16 deter,detAbs
      integer i,j
      Logical Singular,singular1
      common/DetCountML/Singular
c New for Gram
      common/GramSwitch/GramCS,GramDS
      logical GramCS,GRAMDS


cFC      real*16 tempjj,tempkl,IX,Zmax
c      Common/Decide/tempjj,tempkl,IX,ZMax,jjinit
c      Save/Decide/
c New variables to compare        
c      real*16 musq
c      Common/musqInv/musq
cFC      complex*32 Cij(30,9),C0123
cFC      real*16 zz11,zz12,zz22
cFC      real*16 cijrtemp(9,4),cijitemp(9,4),C0ITEMP,C0RTEMP
cFC      complex*32 test2g,test2p1p2,ten2g, ten2g_Dble,ten2p1p2
cFC      complex*32 ten2p1p2_Dble,test2p1p2_dble
cFC      complex*32 test3gp2,ten3gp2,ten3gp2_dble
cFC      complex*32 test2g_dble,test3gp2_dble

cFC      complex*32 A0finG,A02,A03
cFC      external A0finG
cFC      real*16 accuracyC(0:4,5),AccuracyD(0:5,4)
cFC      real*16 accuracyCR(9,0:4,5)
cFC      Common/Accuracy/AccuracyC,AccuracyD
      real*16 ratio
       logical printmy
       common/mprint/printmy
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ratio=1q-6 ! use at the end

      singular1=.false.
      p1p2 = (s12 - p1sq - p2sq)*0.5q0

      m0sq=m0*m0
      m1sq=m1*m1
      m2sq=m2*m2

       r10 = p1sq-m1sq+m0sq
       r21 = s12-m2sq -(p1sq-m1sq)
       r1 = r10
       r2r1 = r21


       deter = abs(2.q0*(p1sq*p2sq - p1p2*p1p2))
       detAbs=  abs(2.q0*(abs(p1sq*p2sq)+abs(p1p2*p1p2)))
      
c       If( (deter/detAbs).le.1q-5) Singular=.true.

       If( abs(deter)/abs(detAbs).lt.4.6q-9) then
         Singular=.true.
         singular1=.true.
c         print*, "Here1",abs(deter)/abs(detAbs),Singular, GramCS
         do i=1,4 
            do j=1,9
               Cijr(j,i)=0q0
               CijI(j,i)=0q0
               C0r=0q0
               C0I=0q0
            enddo
         enddo
c         stop
c         return
         endif

 
       C0r=REAL(C0)
       C0I=IMAG(C0)


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
          z11=2q0*p1sq
          iz11=1q0/z11
          z12=2q0*p1p2 
          z21=z12*iz11
          z22=2q0*p2sq-z12*z21
          iz22=1q0/z22
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

      Cijr(4,2) = ( B0r_23 +2q0*m0sq*C0r + Cijr(1,1)*r1 +Cijr(2,1)*r2r1 +1.q0)*0.25q0
      CijI(4,2) = ( B0I_23 +2q0*m0sq*C0I + CijI(1,1)*r1 +CijI(2,1)*r2r1)*0.25q0    
c 3-5
       PRr(1) = (B1r_13 + B0r_23 - Cijr(1,1)*r1 - Cijr(4,2)*2.q0)
       PRr(2) = (B1r_12 - B1r_13 - Cijr(1,1)*r2r1)

       Rr(2)=(PRr(2)-z21*PRr(1))*iz22 
       Rr(1)=(PRr(1)-z12*Rr(2))*iz11 

       Cijr(1,2) =Rr(1)
       Cijr(3,2) =Rr(2)

       PRI(1) = (B1I_13 + B0I_23 - CijI(1,1)*r1 - CijI(4,2)*2.q0)
       PRI(2) = (B1I_12 - B1I_13 - CijI(1,1)*r2r1)

       RI(2)=(PRI(2)-z21*PRI(1))*iz22 
       RI(1)=(PRI(1)-z12*RI(2))*iz11 

       CijI(1,2) =RI(1)
       CijI(3,2) =RI(2)

c 4-6
       PRr(1) = (  B1r_13 - B1r_23 - Cijr(2,1)*r1)
       PRr(2) = (- B1r_13          -Cijr(2,1)*r2r1 -Cijr(4,2)*2.q0)
      
       Rr(2)=(PRr(2)-z21*PRr(1))*iz22 
c      Rr(1)=(PRr(1)-z12*Rr(2))*iz11 

       Cijr(2,2) = Rr(2)
c      Cijr(3,2) = Rr(1)


       PRI(1) = (  B1I_13 - B1I_23 - CijI(2,1)*r1)
       PRI(2) = (- B1I_13          - CijI(2,1)*r2r1 -CijI(4,2)*2.q0)

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

       PRr(1) = -B0r_23 + B21r_13 - r1* Cijr(1,2) - 4q0* Cijr(5,3) 
       PRr(2) = B21r_12 - B21r_13 - r2r1* Cijr(1,2)
   
       Rr(2)=(PRr(2)-z21*PRr(1))*iz22 
       Rr(1)=(PRr(1)-z12*Rr(2))*iz11 

       Cijr(1,3) = Rr(1)
       Cijr(3,3) = Rr(2)

       PrI(1) = -B0I_23 + B21I_13 - r1* CijI(1,2) - 4q0* CijI(5,3) 
       PrI(2) = B21I_12 - B21I_13 - r2r1* CijI(1,2)
   
       RI(2)=(PrI(2)-z21*PrI(1))*iz22 
       RI(1)=(PrI(1)-z12*RI(2))*iz11 

       CijI(1,3) = RI(1)
       CijI(3,3) = RI(2)

c  121,221

       PRr(1)=B21r_13 + B1r_23 - r1* Cijr(3, 2) - 2q0* Cijr(6, 3)
       PRr(2)=-B21r_13 - r2r1* Cijr(3, 2) - 2q0* Cijr(5, 3)
       
       Rr(2)=(PRr(2)-z21*PRr(1))*iz22 
c       Rr(1)=(PRr(1)-z12*Rr(2))*iz11 

c       Cijr(1,1) = Rr(1)
       Cijr(4,3) = Rr(2)
******
       PRI(1)=B21I_13 + B1I_23 - r1* CijI(3, 2) - 2q0* CijI(6, 3)
       PRI(2)=-B21I_13 - r2r1* CijI(3, 2) - 2q0* CijI(5, 3)
       
       RI(2)=(PRI(2)-z21*PRI(1))*iz22 
c       RI(1)=(PRI(1)-z12*RI(2))*iz11 

c       CijI(1,1) = RI(1)
       CijI(4,3) = RI(2)

c  122,222
       PRr(1)=B21r_13 - B21r_23 - r1*Cijr(2, 2)
       PRr(2)=-B21r_13 - r2r1* Cijr(2, 2) - 4q0*Cijr(6, 3)

       Rr(2)=(PRr(2)-z21*PRr(1))*iz22 
c       Rr(1)=(PRr(1)-z12*Rr(2))*iz11 

c       Cijr(1,1) = Rr(1)
       Cijr(2,3) = Rr(2)


       PrI(1)=B21I_13 - B21I_23  - r1*CijI(2, 2)
       PrI(2)=-B21I_13 - r2r1* CijI(2, 2) - 4q0*CijI(6, 3)

       RI(2)=(PRI(2)-z21*PRI(1))*iz22 
c       RI(1)=(PRI(1)-z12*RI(2))*iz11 

c       CijI(1,1) = RI(1)
       CijI(2,3) = RI(2)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 0000

       Ct=(4*m0sq + 4*m1sq + 4*m2sq - p1sq - 
     &     p2sq - s12)/192.q0 + 
     &  (CMPLX(B22R_23,B22I_23,16) 
     & +2*m0sq*CMPLX(CijR(4,2),CijI(4,2),16) + 
     &     r10*CMPLX(CijR(5,3),CijI(5,3),16) 
     &   + r21*CMPLX(CijR(6,3),CijI(6,3),16)
     &     )/8.q0
     
       CijR(9,4)=REAL(ct)
       CijI(9,4)=IMAG(ct)


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


c       if(Singular1)goto 300


       return
ccccccccccccccccc
cccccccccccccccc
ccccccccccccccccc
       else
		  z11=2q0*p1p2
		  iz11=1q0/z11
		  z21=2q0*p1sq*iz11
		  z12=2q0*p2sq
		  z22=z11-z12*z21
		  iz22=1q0/z22
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

      Cijr(4,2) = ( B0r_23 +2q0*m0sq*C0r + Cijr(1,1)*r1 +Cijr(2,1)*r2r1 +1.q0)*0.25q0
      CijI(4,2) = ( B0I_23 +2q0*m0sq*C0I + CijI(1,1)*r1 +CijI(2,1)*r2r1)*0.25q0
c 3-5
       PRr(1) = (B1r_13 + B0r_23 - Cijr(1,1)*r1 - Cijr(4,2)*2.q0)
       PRr(2) = (B1r_12 - B1r_13 - Cijr(1,1)*r2r1)
       
       Rr(2)=(PRr(1)-z21*PRr(2))*iz22 
       Rr(1)=(PRr(2)-z12*Rr(2))*iz11 
       
       Cijr(1,2) =Rr(1)
       Cijr(3,2) =Rr(2)
       
       PRI(1) = (B1I_13 + B0I_23 - CijI(1,1)*r1 - CijI(4,2)*2.q0)
       PRI(2) = (B1I_12 - B1I_13 - CijI(1,1)*r2r1)
       
       RI(2)=(PRI(1)-z21*PRI(2))*iz22 
       RI(1)=(PRI(2)-z12*RI(2))*iz11 
       
       CijI(1,2) =RI(1)
       CijI(3,2) =RI(2)

c 4-6
       PRr(1) = (  B1r_13 - B1r_23 - Cijr(2,1)*r1)
       PRr(2) = (- B1r_13          -Cijr(2,1)*r2r1 -Cijr(4,2)*2.q0)
       
       Rr(2)=(PRr(1)-z21*PRr(2))*iz22 
c      Rr(1)=(PRr(2)-z12*Rr(2))*iz11 

       Cijr(2,2) = Rr(2)
c      Cijr(3,2) = Rr(1)


       PRI(1) = (  B1I_13 - B1I_23 - CijI(2,1)*r1)
       PRI(2) = (- B1I_13          - CijI(2,1)*r2r1 -CijI(4,2)*2.q0)
       
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
       PRr(1) = -B0r_23 + B21r_13 - r1* Cijr(1,2) - 4q0* Cijr(5,3) 
       PRr(2) = B21r_12 - B21r_13 - r2r1* Cijr(1,2)
   
       Rr(2)=(PRr(1)-z21*PRr(2))*iz22 
       Rr(1)=(PRr(2)-z12*Rr(2))*iz11 

       Cijr(1,3) = Rr(1)
       Cijr(3,3) = Rr(2)

       PrI(1) = -B0I_23 + B21I_13 - r1* CijI(1,2) - 4q0* CijI(5,3) 
       PrI(2) = B21I_12 - B21I_13 - r2r1* CijI(1,2)
   
       RI(2)=(PrI(1)-z21*PrI(2))*iz22 
       RI(1)=(PrI(2)-z12*RI(2))*iz11 

       CijI(1,3) = RI(1)
       CijI(3,3) = RI(2)
c  121,221
       PRr(1)=B21r_13 + B1r_23 - r1* Cijr(3, 2) - 2q0* Cijr(6, 3)
       PRr(2)=-B21r_13 - r2r1* Cijr(3, 2) - 2q0* Cijr(5, 3)
       
       Rr(2)=(PRr(1)-z21*PRr(2))*iz22 
c       Rr(1)=(PRr(2)-z12*Rr(2))*iz11 

c       Cijr(1,1) = Rr(1)
       Cijr(4,3) = Rr(2)
******
       PRI(1)=B21I_13 + B1I_23 - r1* CijI(3, 2) - 2q0* CijI(6, 3)
       PRI(2)=-B21I_13 - r2r1* CijI(3, 2) - 2q0* CijI(5, 3)
       
       RI(2)=(PRI(1)-z21*PRI(2))*iz22 
c       RI(1)=(PRI(2)-z12*RI(2))*iz11 

c       CijI(1,1) = RI(1)
       CijI(4,3) = RI(2)

c  122,222
       PRr(1)=B21r_13 - B21r_23 - r1*Cijr(2, 2)
       PRr(2)=-B21r_13 - r2r1* Cijr(2, 2) - 4q0*Cijr(6, 3)

       Rr(2)=(PRr(1)-z21*PRr(2))*iz22 
c       Rr(1)=(PRr(2)-z12*Rr(2))*iz11 

c       Cijr(1,1) = Rr(1)
       Cijr(2,3) = Rr(2)


       PrI(1)=B21I_13 - B21I_23  - r1*CijI(2, 2)
       PrI(2)=-B21I_13 - r2r1* CijI(2, 2) - 4q0*CijI(6, 3)

       RI(2)=(PRI(1)-z21*PRI(2))*iz22 
c       RI(1)=(PRI(2)-z12*RI(2))*iz11 

c       CijI(1,1) = RI(1)
       CijI(2,3) = RI(2)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 0000
c       Cij(9,4)=(4*m0sq + 4*m1sq + 4*m2sq - p1sq - 
c     &     p2sq - s12)/192.q0 + 
c     &  (B23(2,2) + 2*m0sq*Cij(4,2) + 
c     &     r10*Cij(5,3) + r21*Cij(6,3)
c     &     )/8.q0

       Ct=(4*m0sq + 4*m1sq + 4*m2sq - p1sq - 
     &     p2sq - s12)/192.q0 + 
     &  (CMPLX(B22R_23,B22I_23,16) 
     & +2*m0sq*CMPLX(CijR(4,2),CijI(4,2),16) + 
     &     r10*CMPLX(CijR(5,3),CijI(5,3),16) 
     &   + r21*CMPLX(CijR(6,3),CijI(6,3),16)
     &     )/8.q0
     
       CijR(9,4)=REAL(ct)
       CijI(9,4)=IMAG(ct)  

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


c       if(Singular1) goto 300  
       endif


cFC 300   if(GramCS.and.singular1) then
cFCc          print*,"ActiveR4"
cFC            if(m0*m1*m2.eq.0q0) then
cFC               
cFC               A02=0q0
cFC               A03=0q0
cFC               C0123=C0
cFC               B0_23=DCMPLX(B0R_23,B0I_23)
cFC               B0_12=DCMPLX(B0R_12,B0I_12)
cFC               B0_13=DCMPLX(B0R_13,B0I_13)
cFCc               print*, p1sq,p2sq,s12,musq
cFCc               print*, "masses",m0,m1,m2
cFCc               print*, "B0_23",B0_23
cFCc               print*, "B0_13",B0_13
cFCc               print*, "B0_12",B0_12
cFCc               print*,  "C123",C0123
cFC           call ten_red3_Gram1(p1sq,p2sq,s12,musq,B0_23,B0_13,B0_12,C0123,Cij)
cFC           elseif(abs(m0**2/(m1*m2)-1q0).lt.1d-3) then
cFC               C0123=C0
cFC               B0_23=DCMPLX(B0R_23,B0I_23)
cFC               B0_12=DCMPLX(B0R_12,B0I_12)
cFC               B0_13=DCMPLX(B0R_13,B0I_13)
cFC
cFC               A02=A0finG(M1,musq)  
cFC               A03=A0finG(M2,musq)  
cFCc                 print*, "masses",m0,m1,m2          
cFC                call ten_red3_Gram2_M(m0,p1sq,p2sq,s12,musq,B0_23,B0_13,B0_12,C0123,Cij)
cFC          
cFC          endif
cFC
cFC
cFC         if(printmy) then
cFC           print*, "c0",c0
cFC           PRINT*, "c0123gram",c0123
cFC           PRINT*
cFC           PRINT*,  "ix",IX,TEMPJJ,TEMPKL
cFC           do i=1,5
cFC         print*, "AccuracyC", i,AccuracyC(0,i),AccuracyC(1,i),AccuracyC(2,i),AccuracyC(3,i),AccuracyC(4,i)
cFC         enddo
cFC         endif
cFC
cFC
cFCc              print*,  "C123",C0123
cFC
cFC       index(1)=2
cFC       index(2)=4
cFC       index(3)=6
cFC       index(4)=9
cFC
cFC       do j1=1,4
cFC          do j2= 1,index(j1)
cFC             CijRtemp(j2,j1)=Dble(Cij(j2,j1))
cFC             CijItemp(j2,j1)=imag(Cij(j2,j1))     
cFCc             print*, j2,j1, CijRtemp(j2,j1),CijItemp(j2,j1)
cFCc             print*, "Complex",Cij(j2,j1)
cFCc             print*, "Dble", DCMPLX(CijR(j2,j1),CijI(j2,j1))
cFCc             print*, j2,j1, DCMPLX(CijRtemp(j2,j1),CijItemp(j2,j1))/(DCMPLX(CijR(j2,j1),CijI(j2,j1)))
cFC          enddo
cFC       enddo
cFC        C0rTEMP=Dble(C0123)
cFC        C0ITEMP=DImag(C0123)
cFC     
cFC	
cFC	zz11=2q0*p1sq
cFC	zz12=2q0*p1p2
cFC	zz22=2q0*p2sq
cFC 
cFCc Checking accuracy 
cFC       test2g=2q0*(B0_23+ C0123*m0sq)
cFC       test2g_Dble=2q0*(B0_23+ C0*m0sq)
cFCc       print*, "C0123",C0123,p1sq,p2sq,s12,m0sq,m0sq,m0sq
cFCc       print*, "C0",C0,p1sq,p2sq,s12,m0sq,m0sq,m0sq
cFC       
cFCc Only valid for equal masses of massless case:
cFC
cFC       test2p1p2=-B0_13*(ZZ12+ZZ22)/2q0 +(r2r1)*B0_23
cFC     1   -r1*(B0_12-B0_13-(r2r1)*C0123)
cFC
cFC       test2p1p2_dble=-B0_13*(ZZ12+ZZ22)/2q0 +(r2r1)*B0_23
cFC     1   -r1*(B0_12-B0_13-(r2r1)*C0)
cFC
cFC
cFCc General Case:
cFCc       test3gp2 =-(ZZ12+Z22/2q0)*B0_23
cFC
cFC       test3gp2= A02 + (-2*A03 + 2*B0_12*m0sq - 2*B0_13*m0sq - 
cFC     -     (B0_23 + C0123*m0sq)*
cFC     -      (2*(m1sq - m2sq + ZZ12) + ZZ22))/2.q0
cFC
cFC       test3gp2_dble= A02 + (-2*A03 + 2*B0_12*m0sq - 2*B0_13*m0sq - 
cFC     -     (B0_23 + C0*m0sq)*
cFC     -      (2*(m1sq - m2sq + ZZ12) + ZZ22))/2.q0
cFC
cFC
cFC
cFCcFC        ZZ12/3. + ZZ22/6. + 
cFCcFC     -  (ZZ11*ZZ12*Cij123(1,3))/2. + 
cFCcFC     -  (ZZ22**2*Cij123(2,3))/2. + 
cFCcFC     -  ZZ12**2*Cij123(3,3) + 
cFCcFC     -  (ZZ11*ZZ22*Cij123(3,3))/2. + 
cFCcFC     -  (3*ZZ12*ZZ22*Cij123(4,3))/2. + 
cFCcFC     -  6*ZZ12*Cij123(5,3) + 
cFCcFC     -  6*ZZ22*Cij123(6,3)
cFC
cFC       ten3gp2 = (zz11*zz12*DCMPLX(CijRtemp(1,3),CijItemp(1,3)))/2.q0 +
cFC     -  (ZZ22**2*DCMPLX(CijRtemp(2,3),CijItemp(2,3)))/2.q0 + 
cFC     -  zz12**2*DCMPLX(CijRtemp(3,3),CijItemp(3,3)) + 
cFC     -  (zz11*ZZ22*DCMPLX(CijRtemp(3,3),CijItemp(3,3)))/2.q0 +
cFC     -   (3*zz12*ZZ22*DCMPLX(CijRtemp(4,3),CijItemp(4,3)))/2.q0 
cFC     -  + 6*zz12*DCMPLX(CijRtemp(5,3),CijItemp(5,3)) + 
cFC     -  6*ZZ22*DCMPLX(CijRtemp(6,3),CijItemp(6,3))+
cFC     -   ZZ12/3q0  + ZZ22/6q0
cFC
cFC       ten3gp2_Dble = (zz11*zz12*DCMPLX(CijR(1,3),CijI(1,3)))/2.q0 + 
cFC     -  (ZZ22**2*DCMPLX(CijR(2,3),CijI(2,3)))/2.q0 + 
cFC     -  zz12**2*DCMPLX(CijR(3,3),CijI(3,3)) + 
cFC     -  (zz11*ZZ22*DCMPLX(CijR(3,3),CijI(3,3)))/2.q0 +
cFC     -   (3*zz12*ZZ22*DCMPLX(CijR(4,3),CijI(4,3)))/2.q0 
cFC     - + 6*zz12*DCMPLX(CijR(5,3),CijI(5,3)) + 
cFC     -  6*ZZ22*DCMPLX(CijR(6,3),CijI(6,3))+
cFC     -   ZZ12/3q0  + ZZ22/6q0
cFC
cFCc       print*, "test3gp2",test3gp2,test3gp2_dble
cFCc       print*, "ten3gp2",ten3gp2,ten3gp2/test3gp2-1q0
cFCc       print*, "ten3gp2_dle",ten3gp2_dble,ten3gp2_dble/test3gp2_dble-1q0
cFCc       print*, "ratio Dble_GRAM", ten3gp2_dble/ten3gp2
cFCc
cFC           ten2g=ZZ11*DCMPLX(CijRtemp(1,2),CijItemp(1,2))+ZZ22*DCMPLX(CijRtemp(2,2),CijItemp(2,2))+ 
cFC     1  2q0*(ZZ12*DCMPLX(CijRtemp(3,2),CijItemp(3,2)) + 4q0*DCMPLX(CijRtemp(4,2),CijItemp(4,2)))
cFC     2 -1q0
cFC   
cFC      ten2g_Dble=ZZ11*DCMPLX(CijR(1,2),CijI(1,2))+ZZ22*DCMPLX(CijR(2,2),CijI(2,2))+ 
cFC     1  2q0*(ZZ12*DCMPLX(CijR(3,2),CijI(3,2)) + 4q0*DCMPLX(CijR(4,2),CijI(4,2)))
cFC     2 -1q0
cFC   
cFC   
cFC      
cFC           ten2p1p2=ZZ11*(ZZ12*DCMPLX(CijRtemp(1,2),CijItemp(1,2))+ZZ22*DCMPLX(CijRtemp(3,2),CijItemp(3,2))) + 
cFC     1  ZZ12*(ZZ22*DCMPLX(CijRtemp(2,2),CijItemp(2,2)) + ZZ12*DCMPLX(CijRtemp(3,2),CijItemp(3,2))+2*DCMPLX(CijRtemp(4,2),CijItemp(4,2)))
cFCc          
cFC 
cFC      ten2p1p2_Dble=ZZ11*(ZZ12*DCMPLX(CijR(1,2),CijI(1,2))+ZZ22*DCMPLX(CijR(3,2),CijI(3,2))) + 
cFC     1  ZZ12*(ZZ22*DCMPLX(CijR(2,2),CijI(2,2)) + ZZ12*DCMPLX(CijR(3,2),CijI(3,2))+2*DCMPLX(CijR(4,2),CijI(4,2)))
cFC
cFC
cFCc        print*, "test2g", test2g, test2g_dble
cFCc        print*, "ten2p1p2",test2p1p2,test2p1p2_dble
cFCc        print*, 'ten2g',ten2g
cFCc        print*, 'ten2p1p2',ten2p1p2
cFCc        print*, 'ten2g_dble',ten2g_Dble
cFCc        print*, 'ten2p1p2_Dble',ten2p1p2_Dble	  
cFC
cFC
cFC	 If(abs(test2g).gt.1d-6) then
cFC         ten2g     =abs(ten2g/test2g-1q0)
cFC	 ten2g_Dble=abs(ten2g_Dble/test2g_dble-1q0)
cFC	  else
cFC	 ten2g     =abs(ten2g-test2g)
cFC	 ten2g_Dble=abs(ten2g_Dble-test2g_dble)
cFC	 endif
cFC 	  
cFC	 If(abs(test3gp2).gt.1d-6) then
cFC         ten3gp2     =abs(ten3gp2/test3gp2-1q0)
cFC	 ten3gp2_Dble=abs(ten3gp2_Dble/test3gp2_dble-1q0)
cFC	  else
cFC	 ten3gp2     =abs(ten3gp2-test3gp2)
cFC	 ten3gp2_Dble=abs(ten3gp2_Dble-test3gp2_dble)
cFC	 endif
cFC
cFC
cFC
cFC
cFC	  If(abs(test2p1p2_dble).gt.1d-6) then
cFC         ten2p1p2     =abs(ten2p1p2/test2p1p2-1q0)
cFC	 ten2p1p2_Dble=abs(ten2p1p2_Dble/test2p1p2_dble-1q0)
cFC	  else
cFC	 ten2p1p2     =abs(ten2p1p2-test2p1p2)
cFC	 ten2p1p2_Dble=abs(ten2p1p2_Dble-test2p1p2_dble)
cFC	 endif
cFC
cFCc        print*, 'ten2g',ten2g
cFCc        print*, 'ten2p1p2',ten2p1p2
cFC
cFCc        print*, 'ten2g_dble',ten2g_Dble
cFCc        print*, 'ten2p1p2_Dble',ten2p1p2_Dble	  
cFC
cFC
cFC	ten2g_dble=Max(abs(ten2g_dble),abs(ten2p1p2_Dble),abs(ten3gp2_dble)) 
cFC	ten2g=Max(abs(ten2g),abs(ten2p1p2),abs(ten3gp2)) 
cFC
cFC 
cFC        
cFCc        print*, 'ten2g',ten2g
cFCc        print*, 'ten2g_dble',ten2g_Dble
cFC
cFC         if(abs(ten2g).lt.abs(ten2g_dble)) then
cFCc         write(*,'(A12,9E12.3)'),"GOOD RESCUED",abs(ten2g),abs(ten2g_Dble),abs(deter)/abs(detAbs),abs(IX*deter),accuracyC(4,5),accuracyC(3,5),accuracyC(2,5),accuracyC(1,5),accuracyC(0,5)
cFC         else
cFCc         write(*,'(A12,9E12.3)'),"BAD RESCUED",abs(ten2g),abs(ten2g_Dble),abs(deter)/abs(detAbs),abs(IX*deter),accuracyC(4,5),accuracyC(3,5),accuracyC(2,5),accuracyC(1,5),accuracyC(0,5)
cFC         endif
cFC
cFC         if(abs(ten2g).lt.(abs(ten2g_dble))) then
cFC
cFC         do j1=1,4
cFC          do j2= 1,index(j1)
cFC             CijR(j2,j1)=CijRtemp(j2,j1)
cFC             CijI(j2,j1)=CijItemp(j2,j1)
cFC          enddo
cFC        enddo
cFC        C0r=C0rTEMP
cFC        C0I=C0ITEMP
cFC        C0=C0123
cFC
cFC       test2g_Dble=2q0*(B0_23+ C0*m0sq)
cFC       test2p1p2_dble=-B0_13*(ZZ12+ZZ22)/2q0 +(r2r1)*B0_23
cFC     1   -r1*(B0_12-B0_13-(r2r1)*C0)
cFCc General Case:
cFC       test3gp2_dble= A02 + (-2*A03 + 2*B0_12*m0sq - 2*B0_13*m0sq - 
cFC     -     (B0_23 + C0*m0sq)*
cFC     -      (2*(m1sq - m2sq + ZZ12) + ZZ22))/2.q0
cFC       ten3gp2_Dble = (zz11*zz12*DCMPLX(CijR(1,3),CijI(1,3)))/2.q0 + 
cFC     -  (ZZ22**2*DCMPLX(CijR(2,3),CijI(2,3)))/2.q0 + 
cFC     -  zz12**2*DCMPLX(CijR(3,3),CijI(3,3)) + 
cFC     -  (zz11*ZZ22*DCMPLX(CijR(3,3),CijI(3,3)))/2.q0 +
cFC     -   (3*zz12*ZZ22*DCMPLX(CijR(4,3),CijI(4,3)))/2.q0 
cFC     - + 6*zz12*DCMPLX(CijR(5,3),CijI(5,3)) + 
cFC     -  6*ZZ22*DCMPLX(CijR(6,3),CijI(6,3))+
cFC     -   ZZ12/3q0  + ZZ22/6q0
cFC
cFC      ten2g_Dble=ZZ11*DCMPLX(CijR(1,2),CijI(1,2))+ZZ22*DCMPLX(CijR(2,2),CijI(2,2))+ 
cFC     1  2q0*(ZZ12*DCMPLX(CijR(3,2),CijI(3,2)) + 4q0*DCMPLX(CijR(4,2),CijI(4,2)))
cFC     2 -1q0
cFC      ten2p1p2_Dble=ZZ11*(ZZ12*DCMPLX(CijR(1,2),CijI(1,2))+ZZ22*DCMPLX(CijR(3,2),CijI(3,2))) + 
cFC     1  ZZ12*(ZZ22*DCMPLX(CijR(2,2),CijI(2,2)) + ZZ12*DCMPLX(CijR(3,2),CijI(3,2))+2*DCMPLX(CijR(4,2),CijI(4,2)))
cFC
cFC	 If(abs(test2g_dble).gt.1d-6) then
cFC	 ten2g_Dble=abs(ten2g_Dble/test2g_dble-1q0)
cFC	  else
cFC	 ten2g_Dble=abs(ten2g_Dble-test2g_dble)
cFC	 endif
cFC 	  
cFC	 If(abs(test3gp2_dble).gt.1d-6) then
cFC	 ten3gp2_Dble=abs(ten3gp2_Dble/test3gp2_dble-1q0)
cFC	  else
cFC	 ten3gp2_Dble=abs(ten3gp2_Dble-test3gp2_dble)
cFC	 endif
cFC
cFC	  If(abs(test2p1p2_dble).gt.1d-6) then
cFC	 ten2p1p2_Dble=abs(ten2p1p2_Dble/test2p1p2_dble-1q0)
cFC	  else
cFC	 ten2p1p2_Dble=abs(ten2p1p2_Dble-test2p1p2_dble)
cFC	 endif
cFC
cFC	ten2g_dble=Max(abs(ten2g_dble),abs(ten2p1p2_Dble),abs(ten3gp2_dble)) 
cFCc        print*, "here",ten2g_dble/ten2g
cFCc        pause
cFC
cFC        Singular=.false.
cFC
cFC         return
cFC         else
cFC         Singular=.false.
cFC         If( abs(deter)/abs(detAbs).lt.ratio) then
cFC         Singular=.true.
cFCc         singular1=.true.
cFCc         print*, "Here1HereCCRank4",abs(deter)/abs(detAbs),Singular, GramCS
cFC         do i=1,4 
cFC            do j=1,9
cFC               Cijr(j,i)=0q0
cFC               CijI(j,i)=0q0
cFC               C0r=0q0
cFC               C0I=0q0
cFC            enddo
cFC         enddo
cFCc         stop
cFC         return
cFC         endif
cFC         endif
cFC
cFC
cFC         else
cFC         singular=.false.
cFC         If( abs(deter)/abs(detAbs).lt.ratio) then
cFC         Singular=.true.
cFCc         singular1=.true.
cFCc         print*, "Here2HereCCRank4",abs(deter)/abs(detAbs),Singular, GramCS
cFC         do i=1,4 
cFC            do j=1,9
cFC               Cijr(j,i)=0q0
cFC               CijI(j,i)=0q0
cFC               C0r=0q0
cFC               C0I=0q0
cFC            enddo
cFC         enddo
cFCc         stop
cFC         return
cFC         endif
cFC         endif

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
