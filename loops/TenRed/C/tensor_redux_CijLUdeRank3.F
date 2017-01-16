      subroutine tens_red3_new_Re_ComRank3(p1sq,p2sq,s12,B0R_23,B0R_13,B0R_12, 
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
c                       d^4k           1;  k_mu;   k_mu k_nu
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
      complex*16 B0_23, B0_13, B0_12, C0
      real*8 r1, r2r1, det,p1p2
      real*8 B1r_12,B1r_13,B1r_23,Cijr(6,3)
      real*8 B1I_12,B1I_13,B1I_23,CijI(6,3)
      real*8 B0r_23, B0r_13, B0r_12, C0r
      real*8 B0I_23, B0I_13, B0I_12, C0I 
      real*8 z11,z12,z21,z22,iz11,iz22
      real*8 Rr(2),RI(2),PRr(2),PRI(2)
      REAL*8 detabs,deter
      Integer i,j

      real*8 BijR_12(2,2),BijR_13(2,2),BijR_23(2,2)
      real*8 BijI_12(2,2),BijI_13(2,2),BijI_23(2,2)

      real*8 B21r_12,B21r_13,B21r_23
      real*8 B21I_12,B21I_13,B21I_23
      real*8 B22r_12,B22r_13,B22r_23
      real*8 B22I_12,B22I_13,B22I_23
      logical Singular,singular1
       common/DetCountML/Singular
c New for Gram
      common/GramSwitch/GramCS,GramDS
      logical GramCS,GramDS
      Integer jjinit
      real*8 para
      real*8 tempjj,tempkl,IX,Zmax
      Common/Decide/tempjj,tempkl,IX,ZMax,jjinit
      Save/Decide/
c New variables to compare        
      real*8 musq
      Common/musqInv/musq
      Complex*16 Cij(30,9),C0123
      real*8 zz11,zz12,zz22
      real*8 cijrtemp(9,4),cijitemp(9,4),C0ITEMP,C0RTEMP
      complex*16 test2g,test2p1p2,ten2g, ten2g_Dble,ten2p1p2
      complex*16 ten2p1p2_Dble,test2p1p2_dble
      complex*16 test3gp2,ten3gp2,ten3gp2_dble
      complex*16 test2g_dble,test3gp2_dble
      integer j1,j2,index(4)
      complex*16 A0finG,A02,A03
      external A0finG
      real*8 accuracyC(0:4,5),AccuracyD(0:5,4)

      Common/Accuracy/AccuracyC,AccuracyD
      real*8 ratio
      real*8 m0sq,m1sq,m2sq
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ratio=1d-6 ! use at the end

      singular1=.false.
CCCCCCCCC
       para=1.3d-3
CCCCCCC

      p1p2 = (s12 - p1sq - p2sq)*0.5d0

      r1 = p1sq
      r2r1 = s12 - r1
      
       deter = abs(2.d0*(p1sq*p2sq - p1p2*p1p2))
       detAbs=  abs(2.d0*(abs(p1sq*p2sq)+abs(p1p2*p1p2)))


       If( abs(deter)/abs(detAbs).lt.4.6d-3) then
         Singular=.true.
         singular1=.true.
         do i=1,3 
            do j=1,6
               Cijr(j,i)=0d0
               CijI(j,i)=0d0
               C0r=0d0
               C0I=0d0
            enddo
         enddo
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
*********************************************


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
       Cijr(4,2) = ( B0r_23 + Cijr(1,1)*r1 +Cijr(2,1)*r2r1 +1.d0)*0.25d0
       CijI(4,2) = ( B0I_23 + CijI(1,1)*r1 +CijI(2,1)*r2r1)*0.25d0
      
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

       if(Singular1)goto 300
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
       Cijr(4,2) = ( B0r_23 + Cijr(1,1)*r1 +Cijr(2,1)*r2r1 +1.d0)*0.25d0
       CijI(4,2) = ( B0I_23 + CijI(1,1)*r1 +CijI(2,1)*r2r1)*0.25d0
      
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

       if(Singular1) goto 300  


       endif


 300   if(GramCS.and.singular1) then
c          print*,"ActiveR30"
        m0sq=0d0
        m1sq=0d0
        m2sq=0d0
c            if(m*m*m.eq.0d0) then
               
               A02=0d0
               A03=0d0
               C0123=C0
               B0_23=DCMPLX(B0R_23,B0I_23)
               B0_12=DCMPLX(B0R_12,B0I_12)
               B0_13=DCMPLX(B0R_13,B0I_13)
c               print*, p1sq,p2sq,s12,musq
c               print*, "masses",m0,m1,m2
c               print*, "B0_23",B0_23
c               print*, "B0_13",B0_13
c               print*, "B0_12",B0_12
c               print*,  "C123",C0123
           call ten_red3_Gram1(p1sq,p2sq,s12,musq,B0_23,B0_13,B0_12,C0123,Cij)
c           elseif(abs(m**2/(m*m)-1d0).lt.1d-3) then
c               C0123=C0
c               B0_23=DCMPLX(B0R_23,B0I_23)
c               B0_12=DCMPLX(B0R_12,B0I_12)
c               B0_13=DCMPLX(B0R_13,B0I_13)
c
c               A02=A0finG(M,musq)  
c               A03=A0finG(M,musq)  
c                 print*, "masses",m0,m1,m2          
c                call ten_red3_Gram2_M(m,p1sq,p2sq,s12,musq,B0_23,B0_13,B0_12,C0123,Cij)
c          
c          endif


c              print*,  "C123",C0123

       index(1)=2
       index(2)=4
       index(3)=6
       index(4)=9

       do j1=1,3
          do j2= 1,index(j1)
             CijRtemp(j2,j1)=Dble(Cij(j2,j1))
             CijItemp(j2,j1)=DImag(Cij(j2,j1))     
c             print*, j2,j1, CijRtemp(j2,j1),CijItemp(j2,j1)
c             print*, "Complex",Cij(j2,j1)
c             print*, "Dble", DCMPLX(CijR(j2,j1),CijI(j2,j1))
c             print*, j2,j1, DCMPLX(CijRtemp(j2,j1),CijItemp(j2,j1))/(DCMPLX(CijR(j2,j1),CijI(j2,j1)))
          enddo
       enddo
        C0rTEMP=Dble(C0123)
        C0ITEMP=DImag(C0123)
     
	
	zz11=2d0*p1sq
	zz12=2d0*p1p2
	zz22=2d0*p2sq
 
c Checking accuracy 
       test2g=2d0*(B0_23+ C0123*m0sq)
       test2g_Dble=2d0*(B0_23+ C0*m0sq)
c       print*, "C0123",C0123,p1sq,p2sq,s12,m0sq,m0sq,m0sq
c       print*, "C0",C0,p1sq,p2sq,s12,m0sq,m0sq,m0sq
       
c Only valid for equal masses of massless case:

       test2p1p2=-B0_13*(ZZ12+ZZ22)/2d0 +(r2r1)*B0_23
     1   -r1*(B0_12-B0_13-(r2r1)*C0123)

       test2p1p2_dble=-B0_13*(ZZ12+ZZ22)/2d0 +(r2r1)*B0_23
     1   -r1*(B0_12-B0_13-(r2r1)*C0)


c General Case:
c       test3gp2 =-(ZZ12+Z22/2d0)*B0_23

       test3gp2= A02 + (-2*A03 + 2*B0_12*m0sq - 2*B0_13*m0sq - 
     -     (B0_23 + C0123*m0sq)*
     -      (2*(m1sq - m2sq + ZZ12) + ZZ22))/2.d0

       test3gp2_dble= A02 + (-2*A03 + 2*B0_12*m0sq - 2*B0_13*m0sq - 
     -     (B0_23 + C0*m0sq)*
     -      (2*(m1sq - m2sq + ZZ12) + ZZ22))/2.d0


cFC        ZZ12/3. + ZZ22/6. + 
cFC     -  (ZZ11*ZZ12*Cij123(1,3))/2. + 
cFC     -  (ZZ22**2*Cij123(2,3))/2. + 
cFC     -  ZZ12**2*Cij123(3,3) + 
cFC     -  (ZZ11*ZZ22*Cij123(3,3))/2. + 
cFC     -  (3*ZZ12*ZZ22*Cij123(4,3))/2. + 
cFC     -  6*ZZ12*Cij123(5,3) + 
cFC     -  6*ZZ22*Cij123(6,3)

       ten3gp2 = (zz11*zz12*DCMPLX(CijRtemp(1,3),CijItemp(1,3)))/2.d0 +
     -  (ZZ22**2*DCMPLX(CijRtemp(2,3),CijItemp(2,3)))/2.d0 + 
     -  zz12**2*DCMPLX(CijRtemp(3,3),CijItemp(3,3)) + 
     -  (zz11*ZZ22*DCMPLX(CijRtemp(3,3),CijItemp(3,3)))/2.d0 +
     -   (3*zz12*ZZ22*DCMPLX(CijRtemp(4,3),CijItemp(4,3)))/2.d0 
     -  + 6*zz12*DCMPLX(CijRtemp(5,3),CijItemp(5,3)) + 
     -  6*ZZ22*DCMPLX(CijRtemp(6,3),CijItemp(6,3))+
     -   ZZ12/3d0  + ZZ22/6d0

       ten3gp2_Dble = (zz11*zz12*DCMPLX(CijR(1,3),CijI(1,3)))/2.d0 + 
     -  (ZZ22**2*DCMPLX(CijR(2,3),CijI(2,3)))/2.d0 + 
     -  zz12**2*DCMPLX(CijR(3,3),CijI(3,3)) + 
     -  (zz11*ZZ22*DCMPLX(CijR(3,3),CijI(3,3)))/2.d0 +
     -   (3*zz12*ZZ22*DCMPLX(CijR(4,3),CijI(4,3)))/2.d0 
     - + 6*zz12*DCMPLX(CijR(5,3),CijI(5,3)) + 
     -  6*ZZ22*DCMPLX(CijR(6,3),CijI(6,3))+
     -   ZZ12/3d0  + ZZ22/6d0

c       print*, "test3gp2",test3gp2,test3gp2_dble
c       print*, "ten3gp2",ten3gp2,ten3gp2/test3gp2-1d0
c       print*, "ten3gp2_dle",ten3gp2_dble,ten3gp2_dble/test3gp2_dble-1d0
c       print*, "ratio Dble_GRAM", ten3gp2_dble/ten3gp2
c
           ten2g=ZZ11*DCMPLX(CijRtemp(1,2),CijItemp(1,2))+ZZ22*DCMPLX(CijRtemp(2,2),CijItemp(2,2))+ 
     1  2d0*(ZZ12*DCMPLX(CijRtemp(3,2),CijItemp(3,2)) + 4d0*DCMPLX(CijRtemp(4,2),CijItemp(4,2)))
     2 -1d0
   
      ten2g_Dble=ZZ11*DCMPLX(CijR(1,2),CijI(1,2))+ZZ22*DCMPLX(CijR(2,2),CijI(2,2))+ 
     1  2d0*(ZZ12*DCMPLX(CijR(3,2),CijI(3,2)) + 4d0*DCMPLX(CijR(4,2),CijI(4,2)))
     2 -1d0


           ten2p1p2=ZZ11*(ZZ12*DCMPLX(CijRtemp(1,2),CijItemp(1,2))+ZZ22*DCMPLX(CijRtemp(3,2),CijItemp(3,2))) + 
     1  ZZ12*(ZZ22*DCMPLX(CijRtemp(2,2),CijItemp(2,2)) 
     1  + ZZ12*DCMPLX(CijRtemp(3,2),CijItemp(3,2))+2*DCMPLX(CijRtemp(4,2),CijItemp(4,2)))
c          
 
      ten2p1p2_Dble=ZZ11*(ZZ12*DCMPLX(CijR(1,2),CijI(1,2))+ZZ22*DCMPLX(CijR(3,2),CijI(3,2))) + 
     1  ZZ12*(ZZ22*DCMPLX(CijR(2,2),CijI(2,2)) + ZZ12*DCMPLX(CijR(3,2),CijI(3,2))+2*DCMPLX(CijR(4,2),CijI(4,2)))


c        print*, "test2g", test2g, test2g_dble
c        print*, "ten2p1p2",test2p1p2,test2p1p2_dble
c        print*, 'ten2g',ten2g
c        print*, 'ten2p1p2',ten2p1p2
c        print*, 'ten2g_dble',ten2g_Dble
c        print*, 'ten2p1p2_Dble',ten2p1p2_Dble	  


	 If(abs(test2g).gt.1d-6) then
         ten2g     =abs(ten2g/test2g-1d0)
	 ten2g_Dble=abs(ten2g_Dble/test2g_dble-1d0)
	  else
	 ten2g     =abs(ten2g-test2g)
	 ten2g_Dble=abs(ten2g_Dble-test2g_dble)
	 endif
 	  
	 If(abs(test3gp2).gt.1d-6) then
         ten3gp2     =abs(ten3gp2/test3gp2-1d0)
	 ten3gp2_Dble=abs(ten3gp2_Dble/test3gp2_dble-1d0)
	  else
	 ten3gp2     =abs(ten3gp2-test3gp2)
	 ten3gp2_Dble=abs(ten3gp2_Dble-test3gp2_dble)
	 endif


	  If(abs(test2p1p2_dble).gt.1d-6) then
         ten2p1p2     =abs(ten2p1p2/test2p1p2-1d0)
	 ten2p1p2_Dble=abs(ten2p1p2_Dble/test2p1p2_dble-1d0)
	  else
	 ten2p1p2     =abs(ten2p1p2-test2p1p2)
	 ten2p1p2_Dble=abs(ten2p1p2_Dble-test2p1p2_dble)
	 endif

c        print*, 'ten2g',ten2g
c        print*, 'ten2p1p2',ten2p1p2

c        print*, 'ten2g_dble',ten2g_Dble
c        print*, 'ten2p1p2_Dble',ten2p1p2_Dble


c Check NAN
        if(ten2g+1d0.eq.ten2g)ten2g= 1d99
        if(ten2p1p2+1d0.eq.ten2p1p2)ten2p1p2= 1d99
        if(ten3gp2+1d0.eq.ten3gp2)ten3gp2 = 1d99
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	ten2g_dble=Max(abs(ten2g_dble),abs(ten2p1p2_Dble),abs(ten3gp2_dble))
         ten2g=Max(abs(ten2g),abs(ten2p1p2),abs(ten3gp2))


c        print*, 'ten2g',ten2g
c        print*, 'ten2g_dble',ten2g_Dble

         if(abs(ten2g).lt.abs(ten2g_dble)) then
c         write(*,'(A12,9E12.3)'),"GOOD RESCUED",abs(ten2g),abs(ten2g_Dble),abs(deter)/abs(detAbs),abs(IX*deter),accuracyC(4,5),accuracyC(3,5),accuracyC(2,5),accuracyC(1,5),accuracyC(0,5)
         else
c         write(*,'(A12,9E12.3)'),"BAD RESCUED",abs(ten2g),abs(ten2g_Dble),abs(deter)/abs(detAbs),abs(IX*deter),accuracyC(4,5),accuracyC(3,5),accuracyC(2,5),accuracyC(1,5),accuracyC(0,5)
         endif

         if(abs(ten2g).lt.(abs(ten2g_dble))) then

         do j1=1,3
          do j2= 1,index(j1)
             CijR(j2,j1)=CijRtemp(j2,j1)
             CijI(j2,j1)=CijItemp(j2,j1)
          enddo
        enddo
        C0r=C0rTEMP
        C0I=C0ITEMP
        C0=C0123

       test2g_Dble=2d0*(B0_23+ C0*m0sq)
       test2p1p2_dble=-B0_13*(ZZ12+ZZ22)/2d0 +(r2r1)*B0_23
     1   -r1*(B0_12-B0_13-(r2r1)*C0)
c General Case:
       test3gp2_dble= A02 + (-2*A03 + 2*B0_12*m0sq - 2*B0_13*m0sq - 
     -     (B0_23 + C0*m0sq)*
     -      (2*(m1sq - m2sq + ZZ12) + ZZ22))/2.d0
       ten3gp2_Dble = (zz11*zz12*DCMPLX(CijR(1,3),CijI(1,3)))/2.d0 + 
     -  (ZZ22**2*DCMPLX(CijR(2,3),CijI(2,3)))/2.d0 + 
     -  zz12**2*DCMPLX(CijR(3,3),CijI(3,3)) + 
     -  (zz11*ZZ22*DCMPLX(CijR(3,3),CijI(3,3)))/2.d0 +
     -   (3*zz12*ZZ22*DCMPLX(CijR(4,3),CijI(4,3)))/2.d0 
     - + 6*zz12*DCMPLX(CijR(5,3),CijI(5,3)) + 
     -  6*ZZ22*DCMPLX(CijR(6,3),CijI(6,3))+
     -   ZZ12/3d0  + ZZ22/6d0

      ten2g_Dble=ZZ11*DCMPLX(CijR(1,2),CijI(1,2))+ZZ22*DCMPLX(CijR(2,2),CijI(2,2))+ 
     1  2d0*(ZZ12*DCMPLX(CijR(3,2),CijI(3,2)) + 4d0*DCMPLX(CijR(4,2),CijI(4,2)))
     2 -1d0
      ten2p1p2_Dble=ZZ11*(ZZ12*DCMPLX(CijR(1,2),CijI(1,2))+ZZ22*DCMPLX(CijR(3,2),CijI(3,2))) + 
     1  ZZ12*(ZZ22*DCMPLX(CijR(2,2),CijI(2,2)) + ZZ12*DCMPLX(CijR(3,2),CijI(3,2))+2*DCMPLX(CijR(4,2),CijI(4,2)))

	 If(abs(test2g_dble).gt.1d-6) then
	 ten2g_Dble=abs(ten2g_Dble/test2g_dble-1d0)
	  else
	 ten2g_Dble=abs(ten2g_Dble-test2g_dble)
	 endif
 	  
	 If(abs(test3gp2_dble).gt.1d-6) then
	 ten3gp2_Dble=abs(ten3gp2_Dble/test3gp2_dble-1d0)
	  else
	 ten3gp2_Dble=abs(ten3gp2_Dble-test3gp2_dble)
	 endif

	  If(abs(test2p1p2_dble).gt.1d-6) then
	 ten2p1p2_Dble=abs(ten2p1p2_Dble/test2p1p2_dble-1d0)
	  else
	 ten2p1p2_Dble=abs(ten2p1p2_Dble-test2p1p2_dble)
	 endif

	ten2g_dble=Max(abs(ten2g_dble),abs(ten2p1p2_Dble),abs(ten3gp2_dble)) 
c        print*, "here",ten2g_dble/ten2g
c        pause

        Singular=.false.

         return
         else
         Singular=.false.
         If( abs(deter)/abs(detAbs).lt.ratio) then
         Singular=.true.
c         singular1=.true.
c         print*, "HereCCRank31",abs(deter)/abs(detAbs),Singular, GramCS
         do i=1,3 
            do j=1,6
               Cijr(j,i)=0d0
               CijI(j,i)=0d0
               C0r=0d0
               C0I=0d0
            enddo
         enddo
c         stop
         return
         endif
         endif


         else
        Singular=.false.
         If( abs(deter)/abs(detAbs).lt.ratio) then
         Singular=.true.
c         singular1=.true.
c         print*, "Here2HereCCRank3",abs(deter)/abs(detAbs),Singular, GramCS
         do i=1,3 
            do j=1,6
               Cijr(j,i)=0d0
               CijI(j,i)=0d0
               C0r=0d0
               C0I=0d0
            enddo
         enddo
c         stop
         return
         endif
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
