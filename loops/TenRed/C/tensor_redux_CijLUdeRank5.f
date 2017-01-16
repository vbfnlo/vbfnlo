       subroutine ten_red_LUdecomR5(p1sq,p2sq,s12,B023,B013,B012, 
     &                     Cijr,CijI,Cij)
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
      real*8 r10, r21, det,p1p2
      real*8 Cijr(4,2)
      real*8 CijI(4,2)
       real*8 z11,z12,z21,z22,iz11,iz22
ccccccccccccccccc
      complex*16 Rr(2),PR(2)
      complex*16 B023, B013, B012
      complex*16 B12(6,11),B13(6,11),B23(6,11)
      complex*16 Cij(30,9)

      p1p2 = (s12 - p1sq - p2sq)*0.5d0

      r10 = p1sq
      r21 = s12 - r10

      call ten_red2_forGram(p1sq,B012,B12)
      call ten_red2_forGram(s12,B013,B13)
      call ten_red2_forGram(p2sq,B023,B23)
  
      Cij(1,1)=DCMPLX(CijR(1,1),CijI(1,1))
      Cij(2,1)=DCMPLX(CijR(2,1),CijI(2,1))

      Cij(1,2)=DCMPLX(CijR(1,2),CijI(1,2))
      Cij(2,2)=DCMPLX(CijR(2,2),CijI(2,2))
      Cij(3,2)=DCMPLX(CijR(3,2),CijI(3,2))
      Cij(4,2)=DCMPLX(CijR(4,2),CijI(4,2))


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
c  001,002
       Pr(1) = B13(2, 2) - B23(2, 2) - r10 *Cij(4, 2)
       Pr(2) = B12(2, 2) - B13(2, 2) - r21 *Cij(4, 2)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 

       Cij(5,3)=Rr(1)
       Cij(6,3)=Rr(2)

c 111,211
       Pr(1) = -B023 + B13(1,2) - r10* Cij(1,2) - 4* Cij(5,3) 
       Pr(2) = B12(1,2) - B13(1,2) - r21* Cij(1,2)
   
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 

       Cij(1,3) = Rr(1)
       Cij(3,3) = Rr(2)

c  121,221

       Pr(1)=B13(1, 2) + B23(1, 1) - r10* Cij(3, 2) - 2* Cij(6, 3)
       Pr(2)=-B13(1, 2) - r21* Cij(3, 2) - 2* Cij(5, 3)
       
       Rr(2)=(PR(2)-z21*Pr(1))*iz22 
c       Rr(1)=(PR(1)-z12*Rr(2))*iz11 

c       Cij(3,3) = Rr(1)
       Cij(4,3) = Rr(2)

c  122,222
       Pr(1)=B13(1, 2) - B23(1, 2) - r10*Cij(2, 2)
       Pr(2)=-B13(1, 2) - r21* Cij(2, 2) - 4*Cij(6, 3)

       Rr(2)=(PR(2)-z21*Pr(1))*iz22 
c       Rr(1)=(PR(1)-z12*Rr(2))*iz11 

c       Cij(4,1) = Rr(1)
       Cij(2,3) = Rr(2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 0000

       Cij(9,4)=(-p1sq - p2sq - s12)/192d0 + 
     &  (B23(2,2) + r10*Cij(5,3) +r21*Cij(6,3))/8d0
       
c 0011,0021
       Pr(1)=B13(2,3) + B23(2,2) - r10*Cij(5,3) - 2*Cij(9,4)
       Pr(2)=B12(2,3) - B13(2,3) - r21*Cij(5,3) 
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
       Cij(6,4)=Rr(1)   
       Cij(8,4)=Rr(2)
c0012,0022
       Pr(1)=B13(2,3) - B23(2,3) - r10*Cij(6,3)
       Pr(2)=-B13(2,3) - r21*Cij(6,3) - 2*Cij(9,4)  
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c      Cij(8,4)=Rr(1)
       Cij(7,4)=Rr(2)

c 1111,2111
       Pr(1)=B023 + B13(1,3) - r10*Cij(1,3) - 6*Cij(6,4)
       Pr(2)=B12(1,3) - B13(1,3) - r21*Cij(1,3) 
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
       Cij(1,4)=Rr(1)
       Cij(3,4)=Rr(2)

c 1211 2211
       Pr(1)=B13(1,3) - B23(1,1) - r10*Cij(3,3) - 4*Cij(8,4)
       Pr(2)=-B13(1,3) - r21*Cij(3,3) - 2*Cij(6,4) 
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(3,4)=Rr(1)
       Cij(4,4)=Rr(2)

c 1221 2221
       Pr(1)=B13(1,3) + B23(1,2) - r10*Cij(4,3) - 2*Cij(7,4)
       Pr(2)=-B13(1,3) - r21*Cij(4,3) - 4*Cij(8,4) 
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(4,4)=Rr(1)
       Cij(5,4)=Rr(2)

c 1222,2222
       Pr(1)=B13(1,3) - B23(1,3) - r10*Cij(2,3)
       Pr(2)=-B13(1,3) - r21*Cij(2,3) - 6*Cij(7,4)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(5,4)=Rr(1)
       Cij(2,4)=Rr(2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 00001,00002
       Pr(1)=B13(3,4) - B23(3,4) - r10*Cij(9,4)
       Pr(2)=B12(3,4) - B13(3,4) - r21*Cij(9,4)   
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
       Cij(11,5)=Rr(1)
       Cij(12,5)=Rr(2)

CCCCCC

c 00111,00211
       Pr(1)=B13(2,4) - B23(2,2) - r10*Cij(6,4) - 4*Cij(11,5)
       Pr(2)=B12(2,4) - B13(2,4) - r21*Cij(6,4)  
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
       Cij(7,5)=Rr(1)
       Cij(9,5)=Rr(2)
c 00121,00221
       Pr(1)=B13(2,4) + B23(2,3) - r10*Cij(8,4) - 2*Cij(12,5)
       Pr(2)=-B13(2,4) - r21*Cij(8,4) - 2*Cij(11,5)  
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(9,5)=Rr(1)
       Cij(10,5)=Rr(2)
c 00122.00222
       Pr(1)=B13(2,4) - B23(2,4) - r10*Cij(7,4)
       Pr(2)=-B13(2,4) - r21*Cij(7,4) - 4*Cij(12,5)   
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(10,5)=Rr(1)
       Cij(8,5)=Rr(2)
CCCCC

c 11111,21111
       Pr(1)=-B023 + B13(1,4) - r10*Cij(1,4) - 8*Cij(7,5)
       Pr(2)= B12(1,4) - B13(1,4) - r21*Cij(1,4)  
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
       Cij(1,5)=Rr(1)
       Cij(3,5)=Rr(2)
c 12111,22111
       Pr(1)=B13(1,4) + B23(1,1) - r10*Cij(3,4) - 6*Cij(9,5)
       Pr(2)=-B13(1,4) - r21*Cij(3,4) - 2*Cij(7,5) 
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(3,5)=Rr(1)
       Cij(4,5)=Rr(2)
c 12211,22211
       Pr(1)=B13(1,4) - B23(1,2) - r10*Cij(4,4) - 4*Cij(10,5)
       Pr(2)=-B13(1,4) - r21*Cij(4,4) - 4*Cij(9,5)  
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(4,5)=Rr(1)
       Cij(5,5)=Rr(2)
c 12221,22221
       Pr(1)=B13(1,4) + B23(1,3) - r10*Cij(5,4) - 2*Cij(8,5)
       Pr(2)=-B13(1,4) - r21*Cij(5,4) - 6*Cij(10,5)   
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(5,5)=Rr(1)
       Cij(6,5)=Rr(2)
c 12222,22222
       Pr(1)=B13(1,4) - B23(1,4) - r10*Cij(2,4)
       Pr(2)=-B13(1,4) - r21*Cij(2,4) - 8*Cij(8,5)  
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(6,5)=Rr(1)
       Cij(2,5)=Rr(2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 000000
      
       Rr(1)=(B23(3,4)+r10*Cij(11,5)+r21*Cij(12,5))/12d0
     &  +(p1sq*p1sq+(p2sq+s12)*p1sq+p2sq*p2sq+s12*s12+p2sq*s12)/8640d0

       Cij(16,6)=Rr(1)
       

c 000011 000021
       Pr(1)=B13(3,5) + B23(3,4) - r10*Cij(11,5) - 2*Cij(16,6)
       Pr(2)=B12(3,5) - B13(3,5) - r21*Cij(11,5)


       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 

       Cij(13,6)=Rr(1)
       Cij(15,6)=Rr(2)

c 000012 000022
       Pr(1)=B13(3,5) - B23(3,5) - r10*Cij(12,5)
       Pr(2)=-B13(3,5) - r21*Cij(12,5) - 2*Cij(16,6)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 

c      Cij(15,6)=Rr(1)
        Cij(14,6)=Rr(2)

c 001111 002111
        
       Pr(1)=B13(2,5) + B23(2,2) - r10*Cij(7,5) - 6*Cij(13,6)
       Pr(2)=B12(2,5) - B13(2,5) - r21*Cij(7,5)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
        Cij(8,6)=Rr(1)
        Cij(10,6)=Rr(2)


c 001211 002211

       Pr(1)=B13(2,5) - B23(2,3) - r10*Cij(9,5) - 4*Cij(15,6)
       Pr(2)=-B13(2,5) - r21*Cij(9,5) - 2*Cij(13,6)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(10,6)=Rr(1)
        Cij(11,6)=Rr(2)

c 001221 002221

       Pr(1)=B13(2,5) + B23(2,4) - r10*Cij(10,5) - 2*Cij(14,6)
       Pr(2)=-B13(2,5) - r21*Cij(10,5) - 4*Cij(15,6)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(11,6)=Rr(1)
        Cij(12,6)=Rr(2)

c 001222 002222
 
       Pr(1)=B13(2,5) - B23(2,5) - r10*Cij(8,5)
       Pr(2)=-B13(2,5) - r21*Cij(8,5) - 6*Cij(14,6)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(12,6)=Rr(1)
        Cij(9,6)=Rr(2)


c 111111 211111

       Pr(1)=B023 + B13(1,5) - r10*Cij(1,5) - 10*Cij(8,6)
       Pr(2)=B12(1,5) - B13(1,5) - r21*Cij(1,5)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
        Cij(1,6)=Rr(1)
        Cij(3,6)=Rr(2)

c 121111 221111

       Pr(1)=B13(1,5) - B23(1,1) - r10*Cij(3,5) - 8*Cij(10,6)
       Pr(2)=-B13(1,5) - r21*Cij(3,5) - 2*Cij(8,6)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(3,6)=Rr(1)
        Cij(4,6)=Rr(2)


c 122111 222111

       Pr(1)=B13(1,5) + B23(1,2) - r10*Cij(4,5) - 6*Cij(11,6)
       Pr(2)=-B13(1,5) - r21*Cij(4,5) - 4*Cij(10,6)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(4,6)=Rr(1)
        Cij(5,6)=Rr(2)


c 122211  222211

       Pr(1)=B13(1,5) - B23(1,3) - r10*Cij(5,5) - 4*Cij(12,6)
       Pr(2)=-B13(1,5) - r21*Cij(5,5) - 6*Cij(11,6)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(5,6)=Rr(1)
        Cij(6,6)=Rr(2)

c 122221  222221

       Pr(1)=B13(1,5) + B23(1,4) - r10*Cij(6,5) - 2*Cij(9,6)
       Pr(2)=-B13(1,5) - r21*Cij(6,5) - 8*Cij(12,6)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(6,6)=Rr(1)
        Cij(7,6)=Rr(2)

c 122222  222222

       Pr(1)=B13(1,5) - B23(1,5) - r10*Cij(2,5)
       Pr(2)=-B13(1,5) - r21*Cij(2,5) - 10*Cij(9,6)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(7,6)=Rr(1)
        Cij(2,6)=Rr(2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 0000001 0000002

       Pr(1)=B13(4,6) - B23(4,6) - r10*Cij(16,6)
       Pr(2)=B12(4,6) - B13(4,6) - r21*Cij(16,6)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
        Cij(19,7)=Rr(1)
        Cij(20,7)=Rr(2)
c 0000111 0000211
       Pr(1)=B13(3,6) - B23(3,4) - r10*Cij(13,6) - 4*Cij(19,7)
       Pr(2)=B12(3,6) - B13(3,6) - r21*Cij(13,6)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
        Cij(15,7)=Rr(1)
        Cij(17,7)=Rr(2)
c 0000121 0000221
       Pr(1)=B13(3,6) + B23(3,5) - r10*Cij(15,6) - 2*Cij(20,7)
       Pr(2)=-B13(3,6) - r21*Cij(15,6) - 2*Cij(19,7)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(17,7)=Rr(1)
        Cij(18,7)=Rr(2)
c 0000122 0000222
       Pr(1)=B13(3,6) - B23(3,6) - r10*Cij(14,6)
       Pr(2)=-B13(3,6) - r21*Cij(14,6) - 4*Cij(20,7)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(18,7)=Rr(1)
        Cij(16,7)=Rr(2)
c 0011111 0021111
       Pr(1)=B13(2,6) - B23(2,2) - r10*Cij(8,6) - 8*Cij(15,7)
       Pr(2)=B12(2,6) - B13(2,6) - r21*Cij(8,6)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 

        Cij(9,7)=Rr(1)
        Cij(11,7)=Rr(2)

c 0012111 0022111
       Pr(1)=B13(2,6) + B23(2,3) - r10*Cij(10,6) - 6*Cij(17,7)
       Pr(2)=-B13(2,6) - r21*Cij(10,6) - 2*Cij(15,7)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(11,7)=Rr(1)
        Cij(12,7)=Rr(2)

c 0012211 0022211
       Pr(1)=B13(2,6) - B23(2,4) - r10*Cij(11,6) - 4*Cij(18,7)
       Pr(2)=-B13(2,6) - r21*Cij(11,6) - 4*Cij(17,7)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(12,7)=Rr(1)
        Cij(13,7)=Rr(2)
c 0012221 0022221

       Pr(1)=B13(2,6) + B23(2,5) - r10*Cij(12,6) - 2*Cij(16,7)
       Pr(2)=-B13(2,6) - r21*Cij(12,6) - 6*Cij(18,7)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(13,7)=Rr(1)
        Cij(14,7)=Rr(2)
c 0012222 0022222
       Pr(1)=B13(2,6) - B23(2,6) - r10*Cij(9,6)
       Pr(2)=-B13(2,6) - r21*Cij(9,6) - 8*Cij(16,7)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(14,7)=Rr(1)
        Cij(10,7)=Rr(2)

c 1111111 2111111
       Pr(1)=-B023 + B13(1,6) - r10*Cij(1,6) - 12*Cij(9,7)
       Pr(2)=B12(1,6) - B13(1,6) - r21*Cij(1,6)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
        Cij(1,7)=Rr(1)
        Cij(3,7)=Rr(2)

c 1211111 2211111
       Pr(1)=B13(1,6) + B23(1,1) - r10*Cij(3,6) - 10*Cij(11,7)
       Pr(2)=-B13(1,6) - r21*Cij(3,6) - 2*Cij(9,7)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(3,7)=Rr(1)
        Cij(4,7)=Rr(2)

c 1221111 2221111
       Pr(1)=B13(1,6) - B23(1,2) - r10*Cij(4,6) - 8*Cij(12,7)
       Pr(2)=-B13(1,6) - r21*Cij(4,6) - 4*Cij(11,7)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(4,7)=Rr(1)
        Cij(5,7)=Rr(2)


c 1222111 2222111
       Pr(1)=B13(1,6) + B23(1,3) - r10*Cij(5,6) - 6*Cij(13,7)
       Pr(2)=-B13(1,6) - r21*Cij(5,6) - 6*Cij(12,7)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(5,7)=Rr(1)
        Cij(6,7)=Rr(2)


c 1222211 2222211
       Pr(1)=B13(1,6) - B23(1,4) - r10*Cij(6,6) - 4*Cij(14,7)
       Pr(2)=-B13(1,6) - r21*Cij(6,6) - 8*Cij(13,7)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(6,7)=Rr(1)
        Cij(7,7)=Rr(2)
c 1222221 2222221
       Pr(1)=B13(1,6) + B23(1,5) - r10*Cij(7,6) - 2*Cij(10,7)
       Pr(2)=-B13(1,6) - r21*Cij(7,6) - 10*Cij(14,7)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(7,7)=Rr(1)
        Cij(8,7)=Rr(2)
c 1222222 2222222 
       Pr(1)=B13(1,6) - B23(1,6) - r10*Cij(2,6)
       Pr(2)=-B13(1,6) - r21*Cij(2,6) - 12*Cij(10,7)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(8,7)=Rr(1)
        Cij(2,7)=Rr(2)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 00000000

        Rr(1)=(-3*p1sq*p1sq*p1sq - 3*p1sq*p1sq*(p2sq + s12) - 
     &     3*(p2sq + s12)*(p2sq*p2sq + s12*s12) - 
     &     p1sq*(3*p2sq*p2sq + 4*p2sq*s12 + 3*s12*s12))/1.29024d6 + 
     &  (B23(4,6) + r10*Cij(19,7) + r21*Cij(20,7))/16d0

        Cij(25,8)=Rr(1)

c 00000011 00000021

       Pr(1)=B13(4,7) + B23(4,6) - r10*Cij(19,7) - 2*Cij(25,8)
       Pr(2)=B12(4,7) - B13(4,7) - r21*Cij(19,7)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 

        Cij(22,8)=Rr(1)
        Cij(24,8)=Rr(2)

c 00000012 00000022
       Pr(1)=B13(4,7) - B23(4,7) - r10*Cij(20,7)
       Pr(2)=-B13(4,7) - r21*Cij(20,7) - 2*Cij(25,8)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(24,8)=Rr(1)
        Cij(23,8)=Rr(2)


c 00001111 00002111

       Pr(1)=B13(3,7) + B23(3,4) - r10*Cij(15,7) - 6*Cij(22,8)
       Pr(2)=B12(3,7) - B13(3,7) - r21*Cij(15,7)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
        Cij(17,8)=Rr(1)
        Cij(19,8)=Rr(2)

c 00001211 00002211

       Pr(1)=B13(3,7) - B23(3,5) - r10*Cij(17,7) - 4*Cij(24,8)
       Pr(2)=-B13(3,7) - r21*Cij(17,7) - 2*Cij(22,8)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(19,8)=Rr(1)
        Cij(20,8)=Rr(2)

c 00001221 00002221
       Pr(1)=B13(3,7) + B23(3,6) - r10*Cij(18,7) - 2*Cij(23,8)
       Pr(2)=-B13(3,7) - r21*Cij(18,7) - 4*Cij(24,8)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(20,8)=Rr(1)
        Cij(21,8)=Rr(2)

c 00001222 00002222
       Pr(1)=B13(3,7) - B23(3,7) - r10*Cij(16,7)
       Pr(2)=-B13(3,7) - r21*Cij(16,7) - 6*Cij(23,8)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(21,8)=Rr(1)
        Cij(18,8)=Rr(2)


c 00111111 00211111

       Pr(1)=B13(2,7) + B23(2,2) - r10*Cij(9,7) - 10*Cij(17,8)
       Pr(2)=B12(2,7) - B13(2,7) - r21*Cij(9,7)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
        Cij(10,8)=Rr(1)
        Cij(12,8)=Rr(2)

c 00121111 00221111

       Pr(1)=B13(2,7) - B23(2,3) - r10*Cij(11,7) - 8*Cij(19,8)
       Pr(2)=-B13(2,7) - r21*Cij(11,7) - 2*Cij(17,8)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(12,8)=Rr(1)
        Cij(13,8)=Rr(2)

c 00122111 00222111
       Pr(1)=B13(2,7) + B23(2,4) - r10*Cij(12,7) - 6*Cij(20,8)
       Pr(2)=-B13(2,7) - r21*Cij(12,7) - 4*Cij(19,8)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(13,8)=Rr(1)
        Cij(14,8)=Rr(2)

c 00122211 00222211
       Pr(1)=B13(2,7) - B23(2,5) - r10*Cij(13,7) - 4*Cij(21,8)
       Pr(2)=-B13(2,7) - r21*Cij(13,7) - 6*Cij(20,8)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(14,8)=Rr(1)
        Cij(15,8)=Rr(2)

c 00122221  00222221
       Pr(1)=B13(2,7) + B23(2,6) - r10*Cij(14,7) - 2*Cij(18,8)
       Pr(2)=-B13(2,7) - r21*Cij(14,7) - 8*Cij(21,8)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(15,8)=Rr(1)
        Cij(16,8)=Rr(2)
c 00122222  00222222
       Pr(1)=B13(2,7) - B23(2,7) - r10*Cij(10,7)
       Pr(2)=-B13(2,7) - r21*Cij(10,7) - 10*Cij(18,8)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(16,8)=Rr(1)
        Cij(11,8)=Rr(2)

c 11111111 21111111

       Pr(1)=B023 + B13(1,7) - r10*Cij(1,7) - 14*Cij(10,8)
       Pr(2)=B12(1,7) - B13(1,7) - r21*Cij(1,7)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
        Cij(1,8)=Rr(1)
        Cij(3,8)=Rr(2)

c 12111111 22111111

       Pr(1)=B13(1,7) - B23(1,1) - r10*Cij(3,7) - 12*Cij(12,8)
       Pr(2)=-B13(1,7) - r21*Cij(3,7) - 2*Cij(10,8)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(3,8)=Rr(1)
        Cij(4,8)=Rr(2)


c 12211111 22211111
       Pr(1)=B13(1,7) + B23(1,2) - r10*Cij(4,7) - 10*Cij(13,8)
       Pr(2)=-B13(1,7) - r21*Cij(4,7) - 4*Cij(12,8)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(4,8)=Rr(1)
        Cij(5,8)=Rr(2)


c 12221111 22221111
       Pr(1)=B13(1,7) - B23(1,3) - r10*Cij(5,7) - 8*Cij(14,8)
       Pr(2)=-B13(1,7) - r21*Cij(5,7) - 6*Cij(13,8)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(5,8)=Rr(1)
        Cij(6,8)=Rr(2)


c 12222111  22222111

       Pr(1)=B13(1,7) + B23(1,4) - r10*Cij(6,7) - 6*Cij(15,8)
       Pr(2)=-B13(1,7) - r21*Cij(6,7) - 8*Cij(14,8)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(6,8)=Rr(1)
        Cij(7,8)=Rr(2)

c 12222211  22222211
       Pr(1)=B13(1,7) - B23(1,5) - r10*Cij(7,7) - 4*Cij(16,8)
       Pr(2)=-B13(1,7) - r21*Cij(7,7) - 10*Cij(15,8)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(7,8)=Rr(1)
        Cij(8,8)=Rr(2)


c 12222221  22222221

       Pr(1)=B13(1,7) + B23(1,6) - r10*Cij(8,7) - 2*Cij(11,8)
       Pr(2)=-B13(1,7) - r21*Cij(8,7) - 12*Cij(16,8)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(8,8)=Rr(1)
        Cij(9,8)=Rr(2)


c 12222222  22222222

       Pr(1)=B13(1,7) - B23(1,7) - r10*Cij(2,7)
       Pr(2)=-B13(1,7) - r21*Cij(2,7) - 14*Cij(11,8)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c        Cij(9,8)=Rr(1)
        Cij(2,8)=Rr(2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c 000000001 000000002
       Pr(1)=B13(5,8) - B23(5,8) - r10*Cij(25,8)
       Pr(2)=B12(5,8) - B13(5,8) - r21*Cij(25,8)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
       Cij(29,9)=Rr(1)
       Cij(30,9)=Rr(2)


c 000000111 000000211
       Pr(1)=B13(4,8) - B23(4,6) - r10*Cij(22,8) - 4*Cij(29,9)
       Pr(2)=B12(4,8) - B13(4,8) - r21*Cij(22,8)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
       Cij(25,9)=Rr(1)
       Cij(27,9)=Rr(2)

c 000000121 000000221
       Pr(1)=B13(4,8) + B23(4,7) - r10*Cij(24,8) - 2*Cij(30,9)
       Pr(2)=-B13(4,8) - r21*Cij(24,8) - 2*Cij(29,9)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(27,9)=Rr(1)
       Cij(28,9)=Rr(2)


c 000000122 000000222
       Pr(1)=B13(4,8) - B23(4,8) - r10*Cij(23,8)
       Pr(2)=-B13(4,8) - r21*Cij(23,8) - 4*Cij(30,9)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(28,9)=Rr(1)
       Cij(26,9)=Rr(2)


c 000011111 000021111
       Pr(1)=B13(3,8) - B23(3,4) - r10*Cij(17,8) - 8*Cij(25,9)
       Pr(2)=B12(3,8) - B13(3,8) - r21*Cij(17,8)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
       Cij(19,9)=Rr(1)
       Cij(21,9)=Rr(2)


c 000012111 000022111
       Pr(1)=B13(3,8) + B23(3,5) - r10*Cij(19,8) - 6*Cij(27,9)
       Pr(2)=-B13(3,8) - r21*Cij(19,8) - 2*Cij(25,9)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(21,9)=Rr(1)
       Cij(22,9)=Rr(2)


c 000012211 000022211
       Pr(1)=B13(3,8) - B23(3,6) - r10*Cij(20,8) - 4*Cij(28,9)
       Pr(2)=-B13(3,8) - r21*Cij(20,8) - 4*Cij(27,9)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(22,9)=Rr(1)
       Cij(23,9)=Rr(2)

c 000012221 000022221
       Pr(1)=B13(3,8) + B23(3,7) - r10*Cij(21,8) - 2*Cij(26,9)
       Pr(2)=-B13(3,8) - r21*Cij(21,8) - 6*Cij(28,9)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(23,9)=Rr(1)
       Cij(24,9)=Rr(2)

c 000012222 000022222
       Pr(1)=B13(3,8) - B23(3,8) - r10*Cij(18,8)
       Pr(2)=-B13(3,8) - r21*Cij(18,8) - 8*Cij(26,9)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(24,9)=Rr(1)
       Cij(20,9)=Rr(2)


c 001111111 002111111
       Pr(1)=B13(2,8) - B23(2,2) - r10*Cij(10,8) - 12*Cij(19,9)
       Pr(2)=B12(2,8) - B13(2,8) - r21*Cij(10,8)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
       Cij(11,9)=Rr(1)
       Cij(13,9)=Rr(2)
       
         
c 001211111 002211111
       Pr(1)=B13(2,8) + B23(2,3) - r10*Cij(12,8) - 10*Cij(21,9)
       Pr(2)=-B13(2,8) - r21*Cij(12,8) - 2*Cij(19,9)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(13,9)=Rr(1)
       Cij(14,9)=Rr(2)
         
         
c 001221111 002221111
       Pr(1)=B13(2,8) - B23(2,4) - r10*Cij(13,8) - 8*Cij(22,9)
       Pr(2)=-B13(2,8) - r21*Cij(13,8) - 4*Cij(21,9)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(14,9)=Rr(1)
       Cij(15,9)=Rr(2)
         
c 001222111 002222111
       Pr(1)=B13(2,8) + B23(2,5) - r10*Cij(14,8) - 6*Cij(23,9)
       Pr(2)=-B13(2,8) - r21*Cij(14,8) - 6*Cij(22,9)

       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(15,9)=Rr(1)
       Cij(16,9)=Rr(2)
         
c 001222211 002222211
       Pr(1)=B13(2,8) - B23(2,6) - r10*Cij(15,8) - 4*Cij(24,9)
       Pr(2)=-B13(2,8) - r21*Cij(15,8) - 8*Cij(23,9)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(16,9)=Rr(1)
       Cij(17,9)=Rr(2)

c 001222221 002222221
       Pr(1)=B13(2,8) + B23(2,7) - r10*Cij(16,8) - 2*Cij(20,9)
       Pr(2)=-B13(2,8) - r21*Cij(16,8) - 10*Cij(24,9)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(17,9)=Rr(1)
       Cij(18,9)=Rr(2)

c 001222222 002222222
       Pr(1)=B13(2,8) - B23(2,8) - r10*Cij(11,8)
       Pr(2)=-B13(2,8) - r21*Cij(11,8) - 12*Cij(20,9)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(18,9)=Rr(1)
       Cij(12,9)=Rr(2)

c 111111111 211111111
       Pr(1)=-B023 + B13(1,8) - r10*Cij(1,8) - 16*Cij(11,9)
       Pr(2)=B12(1,8) - B13(1,8) - r21*Cij(1,8)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
       Cij(1,9)=Rr(1)
       Cij(3,9)=Rr(2)
c 121111111 221111111
       Pr(1)=B13(1,8) + B23(1,1) - r10*Cij(3,8) - 14*Cij(13,9)
       Pr(2)=-B13(1,8) - r21*Cij(3,8) - 2*Cij(11,9)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(3,9)=Rr(1)
       Cij(4,9)=Rr(2)       
       
c 122111111 222111111
       Pr(1)=B13(1,8) - B23(1,2) - r10*Cij(4,8) - 12*Cij(14,9)
       Pr(2)=-B13(1,8) - r21*Cij(4,8) - 4*Cij(13,9)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(4,9)=Rr(1)
       Cij(5,9)=Rr(2)       
c 122211111 222211111
       Pr(1)=B13(1,8) + B23(1,3) - r10*Cij(5,8) - 10*Cij(15,9)
       Pr(2)=-B13(1,8) - r21*Cij(5,8) - 6*Cij(14,9)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(5,9)=Rr(1)
       Cij(6,9)=Rr(2)       
c 122221111 222221111
       Pr(1)=B13(1,8) - B23(1,4) - r10*Cij(6,8) - 8*Cij(16,9)
       Pr(2)=-B13(1,8) - r21*Cij(6,8) - 8*Cij(15,9)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(6,9)=Rr(1)
       Cij(7,9)=Rr(2)
c 122222111 222222111
       Pr(1)=B13(1,8) + B23(1,5) - r10*Cij(7,8) - 6*Cij(17,9)
       Pr(2)=-B13(1,8) - r21*Cij(7,8) - 10*Cij(16,9)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(7,9)=Rr(1)
       Cij(8,9)=Rr(2)
c 122222211 222222211
       Pr(1)=B13(1,8) - B23(1,6) - r10*Cij(8,8) - 4*Cij(18,9)
       Pr(2)=-B13(1,8) - r21*Cij(8,8) - 12*Cij(17,9)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(8,9)=Rr(1)
       Cij(9,9)=Rr(2)
c 122222221 222222221
       Pr(1)=B13(1,8) + B23(1,7) - r10*Cij(9,8) - 2*Cij(12,9)
       Pr(2)=-B13(1,8) - r21*Cij(9,8) - 14*Cij(18,9)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(9,9)=Rr(1)
       Cij(10,9)=Rr(2)


c 122222222 222222222
       Pr(1)=B13(1,8) - B23(1,8) - r10*Cij(2,8)
       Pr(2)=-B13(1,8) - r21*Cij(2,8) - 16*Cij(12,9)
       Rr(2)=(Pr(2)-z21*Pr(1))*iz22 
c       Rr(1)=(Pr(1)-z12*Rr(2))*iz11 
c       Cij(10,9)=Rr(1)
       Cij(2,9)=Rr(2)


cFCc                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
cFC       Cij(36,10)=(2*p1sq**4+2*p1sq**3*(p2sq+s12)+p1sq**2*(2*p2sq**2+3*p
cFC     &  2sq*s12+2*s12**2)+p1sq*(2*p2sq**3+3*p2sq**2*s12+3*p2sq*s12**2+2*
cFC     &  s12**3)+2*(p2sq**4+p2sq**3*s12+p2sq**2*s12**2+p2sq*s12**3+s12**4
cFC     &  ))/4.8384d7+(B23(5,8)+r10*Cij(29,9)+(r21)*Cij(30,9))/20.d0
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, ii1, 1}
cFC       Pr(1)=B13(5,9)+B23(5,8)-r10*Cij(29,9)-2*Cij(36,10)
cFC       Pr(2)=B12(5,9)-B13(5,9)-r21*Cij(29,9)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFC       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFC       Cij(33,10)=Rr(1)
cFC       Cij(35,10)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, ii1, 2}
cFC       Pr(1)=B13(5,9)-B23(5,9)-r10*Cij(30,9)
cFC       Pr(2)=-B13(5,9)-r21*Cij(30,9)-2*Cij(36,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(35,10)=Rr(1)
cFC       Cij(34,10)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 1, 1, 1}
cFC       Pr(1)=B13(4,9)+B23(4,6)-r10*Cij(25,9)-6*Cij(33,10)
cFC       Pr(2)=B12(4,9)-B13(4,9)-r21*Cij(25,9)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFC       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFC       Cij(28,10)=Rr(1)
cFC       Cij(30,10)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 1, 1}
cFC       Pr(1)=B13(4,9)-B23(4,7)-r10*Cij(27,9)-4*Cij(35,10)
cFC       Pr(2)=-B13(4,9)-r21*Cij(27,9)-2*Cij(33,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(30,10)=Rr(1)
cFC       Cij(31,10)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 2, 1}
cFC       Pr(1)=B13(4,9)+B23(4,8)-r10*Cij(28,9)-2*Cij(34,10)
cFC       Pr(2)=-B13(4,9)-r21*Cij(28,9)-4*Cij(35,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(31,10)=Rr(1)
cFC       Cij(32,10)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 2, 2}
cFC       Pr(1)=B13(4,9)-B23(4,9)-r10*Cij(26,9)
cFC       Pr(2)=-B13(4,9)-r21*Cij(26,9)-6*Cij(34,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(32,10)=Rr(1)
cFC       Cij(29,10)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(3,9)+B23(3,4)-r10*Cij(19,9)-10*Cij(28,10)
cFC       Pr(2)=B12(3,9)-B13(3,9)-r21*Cij(19,9)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFC       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFC       Cij(21,10)=Rr(1)
cFC       Cij(23,10)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(3,9)-B23(3,5)-r10*Cij(21,9)-8*Cij(30,10)
cFC       Pr(2)=-B13(3,9)-r21*Cij(21,9)-2*Cij(28,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(23,10)=Rr(1)
cFC       Cij(24,10)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(3,9)+B23(3,6)-r10*Cij(22,9)-6*Cij(31,10)
cFC       Pr(2)=-B13(3,9)-r21*Cij(22,9)-4*Cij(30,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(24,10)=Rr(1)
cFC       Cij(25,10)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(3,9)-B23(3,7)-r10*Cij(23,9)-4*Cij(32,10)
cFC       Pr(2)=-B13(3,9)-r21*Cij(23,9)-6*Cij(31,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(25,10)=Rr(1)
cFC       Cij(26,10)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(3,9)+B23(3,8)-r10*Cij(24,9)-2*Cij(29,10)
cFC       Pr(2)=-B13(3,9)-r21*Cij(24,9)-8*Cij(32,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(26,10)=Rr(1)
cFC       Cij(27,10)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(3,9)-B23(3,9)-r10*Cij(20,9)
cFC       Pr(2)=-B13(3,9)-r21*Cij(20,9)-10*Cij(29,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(27,10)=Rr(1)
cFC       Cij(22,10)=Rr(2)
cFCc                {0, 0, ii1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,9)+B23(2,2)-r10*Cij(11,9)-14*Cij(21,10)
cFC       Pr(2)=B12(2,9)-B13(2,9)-r21*Cij(11,9)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFC       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFC       Cij(12,10)=Rr(1)
cFC       Cij(14,10)=Rr(2)
cFCc                {0, 0, ii1, 2, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,9)-B23(2,3)-r10*Cij(13,9)-12*Cij(23,10)
cFC       Pr(2)=-B13(2,9)-r21*Cij(13,9)-2*Cij(21,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(14,10)=Rr(1)
cFC       Cij(15,10)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,9)+B23(2,4)-r10*Cij(14,9)-10*Cij(24,10)
cFC       Pr(2)=-B13(2,9)-r21*Cij(14,9)-4*Cij(23,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(15,10)=Rr(1)
cFC       Cij(16,10)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,9)-B23(2,5)-r10*Cij(15,9)-8*Cij(25,10)
cFC       Pr(2)=-B13(2,9)-r21*Cij(15,9)-6*Cij(24,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(16,10)=Rr(1)
cFC       Cij(17,10)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(2,9)+B23(2,6)-r10*Cij(16,9)-6*Cij(26,10)
cFC       Pr(2)=-B13(2,9)-r21*Cij(16,9)-8*Cij(25,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(17,10)=Rr(1)
cFC       Cij(18,10)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(2,9)-B23(2,7)-r10*Cij(17,9)-4*Cij(27,10)
cFC       Pr(2)=-B13(2,9)-r21*Cij(17,9)-10*Cij(26,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(18,10)=Rr(1)
cFC       Cij(19,10)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(2,9)+B23(2,8)-r10*Cij(18,9)-2*Cij(22,10)
cFC       Pr(2)=-B13(2,9)-r21*Cij(18,9)-12*Cij(27,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(19,10)=Rr(1)
cFC       Cij(20,10)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(2,9)-B23(2,9)-r10*Cij(12,9)
cFC       Pr(2)=-B13(2,9)-r21*Cij(12,9)-14*Cij(22,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(20,10)=Rr(1)
cFC       Cij(13,10)=Rr(2)
cFCc                {ii1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B023+B13(1,9)-r10*Cij(1,9)-18*Cij(12,10)
cFC       Pr(2)=B12(1,9)-B13(1,9)-r21*Cij(1,9)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFC       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFC       Cij(1,10)=Rr(1)
cFC       Cij(3,10)=Rr(2)
cFCc                {ii1, 2, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,9)-B23(1,1)-r10*Cij(3,9)-16*Cij(14,10)
cFC       Pr(2)=-B13(1,9)-r21*Cij(3,9)-2*Cij(12,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(3,10)=Rr(1)
cFC       Cij(4,10)=Rr(2)
cFCc                {ii1, 2, 2, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,9)+B23(1,2)-r10*Cij(4,9)-14*Cij(15,10)
cFC       Pr(2)=-B13(1,9)-r21*Cij(4,9)-4*Cij(14,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(4,10)=Rr(1)
cFC       Cij(5,10)=Rr(2)
cFCc                {ii1, 2, 2, 2, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,9)-B23(1,3)-r10*Cij(5,9)-12*Cij(16,10)
cFC       Pr(2)=-B13(1,9)-r21*Cij(5,9)-6*Cij(15,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(5,10)=Rr(1)
cFC       Cij(6,10)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,9)+B23(1,4)-r10*Cij(6,9)-10*Cij(17,10)
cFC       Pr(2)=-B13(1,9)-r21*Cij(6,9)-8*Cij(16,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(6,10)=Rr(1)
cFC       Cij(7,10)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,9)-B23(1,5)-r10*Cij(7,9)-8*Cij(18,10)
cFC       Pr(2)=-B13(1,9)-r21*Cij(7,9)-10*Cij(17,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(7,10)=Rr(1)
cFC       Cij(8,10)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(1,9)+B23(1,6)-r10*Cij(8,9)-6*Cij(19,10)
cFC       Pr(2)=-B13(1,9)-r21*Cij(8,9)-12*Cij(18,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(8,10)=Rr(1)
cFC       Cij(9,10)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(1,9)-B23(1,7)-r10*Cij(9,9)-4*Cij(20,10)
cFC       Pr(2)=-B13(1,9)-r21*Cij(9,9)-14*Cij(19,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(9,10)=Rr(1)
cFC       Cij(10,10)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(1,9)+B23(1,8)-r10*Cij(10,9)-2*Cij(13,10)
cFC       Pr(2)=-B13(1,9)-r21*Cij(10,9)-16*Cij(20,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(10,10)=Rr(1)
cFC       Cij(11,10)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(1,9)-B23(1,9)-r10*Cij(2,9)
cFC       Pr(2)=-B13(1,9)-r21*Cij(2,9)-18*Cij(13,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(11,10)=Rr(1)
cFC       Cij(2,10)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ii1}
cFC       Pr(1)=B13(6,10)-B23(6,10)-r10*Cij(36,10)
cFC       Pr(2)=B12(6,10)-B13(6,10)-r21*Cij(36,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFC       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFC       Cij(41,11)=Rr(1)
cFC       Cij(42,11)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, ii1, 1, 1}
cFC       Pr(1)=B13(5,10)-B23(5,8)-r10*Cij(33,10)-4*Cij(41,11)
cFC       Pr(2)=B12(5,10)-B13(5,10)-r21*Cij(33,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFC       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFC       Cij(37,11)=Rr(1)
cFC       Cij(39,11)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, ii1, 2, 1}
cFC       Pr(1)=B13(5,10)+B23(5,9)-r10*Cij(35,10)-2*Cij(42,11)
cFC       Pr(2)=-B13(5,10)-r21*Cij(35,10)-2*Cij(41,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(39,11)=Rr(1)
cFC       Cij(40,11)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, ii1, 2, 2}
cFC       Pr(1)=B13(5,10)-B23(5,10)-r10*Cij(34,10)
cFC       Pr(2)=-B13(5,10)-r21*Cij(34,10)-4*Cij(42,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(40,11)=Rr(1)
cFC       Cij(38,11)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 1, 1, 1, 1}
cFC       Pr(1)=B13(4,10)-B23(4,6)-r10*Cij(28,10)-8*Cij(37,11)
cFC       Pr(2)=B12(4,10)-B13(4,10)-r21*Cij(28,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFC       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFC       Cij(31,11)=Rr(1)
cFC       Cij(33,11)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 1, 1, 1}
cFC       Pr(1)=B13(4,10)+B23(4,7)-r10*Cij(30,10)-6*Cij(39,11)
cFC       Pr(2)=-B13(4,10)-r21*Cij(30,10)-2*Cij(37,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(33,11)=Rr(1)
cFC       Cij(34,11)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 2, 1, 1}
cFC       Pr(1)=B13(4,10)-B23(4,8)-r10*Cij(31,10)-4*Cij(40,11)
cFC       Pr(2)=-B13(4,10)-r21*Cij(31,10)-4*Cij(39,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(34,11)=Rr(1)
cFC       Cij(35,11)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 2, 2, 1}
cFC       Pr(1)=B13(4,10)+B23(4,9)-r10*Cij(32,10)-2*Cij(38,11)
cFC       Pr(2)=-B13(4,10)-r21*Cij(32,10)-6*Cij(40,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(35,11)=Rr(1)
cFC       Cij(36,11)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 2, 2, 2}
cFC       Pr(1)=B13(4,10)-B23(4,10)-r10*Cij(29,10)
cFC       Pr(2)=-B13(4,10)-r21*Cij(29,10)-8*Cij(38,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(36,11)=Rr(1)
cFC       Cij(32,11)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(3,10)-B23(3,4)-r10*Cij(21,10)-12*Cij(31,11)
cFC       Pr(2)=B12(3,10)-B13(3,10)-r21*Cij(21,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFC       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFC       Cij(23,11)=Rr(1)
cFC       Cij(25,11)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(3,10)+B23(3,5)-r10*Cij(23,10)-10*Cij(33,11)
cFC       Pr(2)=-B13(3,10)-r21*Cij(23,10)-2*Cij(31,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(25,11)=Rr(1)
cFC       Cij(26,11)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(3,10)-B23(3,6)-r10*Cij(24,10)-8*Cij(34,11)
cFC       Pr(2)=-B13(3,10)-r21*Cij(24,10)-4*Cij(33,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(26,11)=Rr(1)
cFC       Cij(27,11)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(3,10)+B23(3,7)-r10*Cij(25,10)-6*Cij(35,11)
cFC       Pr(2)=-B13(3,10)-r21*Cij(25,10)-6*Cij(34,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(27,11)=Rr(1)
cFC       Cij(28,11)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(3,10)-B23(3,8)-r10*Cij(26,10)-4*Cij(36,11)
cFC       Pr(2)=-B13(3,10)-r21*Cij(26,10)-8*Cij(35,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(28,11)=Rr(1)
cFC       Cij(29,11)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(3,10)+B23(3,9)-r10*Cij(27,10)-2*Cij(32,11)
cFC       Pr(2)=-B13(3,10)-r21*Cij(27,10)-10*Cij(36,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(29,11)=Rr(1)
cFC       Cij(30,11)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(3,10)-B23(3,10)-r10*Cij(22,10)
cFC       Pr(2)=-B13(3,10)-r21*Cij(22,10)-12*Cij(32,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(30,11)=Rr(1)
cFC       Cij(24,11)=Rr(2)
cFCc                {0, 0, ii1, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,10)-B23(2,2)-r10*Cij(12,10)-16*Cij(23,11)
cFC       Pr(2)=B12(2,10)-B13(2,10)-r21*Cij(12,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFC       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFC       Cij(13,11)=Rr(1)
cFC       Cij(15,11)=Rr(2)
cFCc                {0, 0, ii1, 2, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,10)+B23(2,3)-r10*Cij(14,10)-14*Cij(25,11)
cFC       Pr(2)=-B13(2,10)-r21*Cij(14,10)-2*Cij(23,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(15,11)=Rr(1)
cFC       Cij(16,11)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,10)-B23(2,4)-r10*Cij(15,10)-12*Cij(26,11)
cFC       Pr(2)=-B13(2,10)-r21*Cij(15,10)-4*Cij(25,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(16,11)=Rr(1)
cFC       Cij(17,11)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,10)+B23(2,5)-r10*Cij(16,10)-10*Cij(27,11)
cFC       Pr(2)=-B13(2,10)-r21*Cij(16,10)-6*Cij(26,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(17,11)=Rr(1)
cFC       Cij(18,11)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,10)-B23(2,6)-r10*Cij(17,10)-8*Cij(28,11)
cFC       Pr(2)=-B13(2,10)-r21*Cij(17,10)-8*Cij(27,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(18,11)=Rr(1)
cFC       Cij(19,11)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(2,10)+B23(2,7)-r10*Cij(18,10)-6*Cij(29,11)
cFC       Pr(2)=-B13(2,10)-r21*Cij(18,10)-10*Cij(28,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(19,11)=Rr(1)
cFC       Cij(20,11)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(2,10)-B23(2,8)-r10*Cij(19,10)-4*Cij(30,11)
cFC       Pr(2)=-B13(2,10)-r21*Cij(19,10)-12*Cij(29,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(20,11)=Rr(1)
cFC       Cij(21,11)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(2,10)+B23(2,9)-r10*Cij(20,10)-2*Cij(24,11)
cFC       Pr(2)=-B13(2,10)-r21*Cij(20,10)-14*Cij(30,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(21,11)=Rr(1)
cFC       Cij(22,11)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(2,10)-B23(2,10)-r10*Cij(13,10)
cFC       Pr(2)=-B13(2,10)-r21*Cij(13,10)-16*Cij(24,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(22,11)=Rr(1)
cFC       Cij(14,11)=Rr(2)
cFCc                {ii1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=-B023+B13(1,10)-r10*Cij(1,10)-20*Cij(13,11)
cFC       Pr(2)=B12(1,10)-B13(1,10)-r21*Cij(1,10)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFC       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFC       Cij(1,11)=Rr(1)
cFC       Cij(3,11)=Rr(2)
cFCc                {ii1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,10)+B23(1,1)-r10*Cij(3,10)-18*Cij(15,11)
cFC       Pr(2)=-B13(1,10)-r21*Cij(3,10)-2*Cij(13,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(3,11)=Rr(1)
cFC       Cij(4,11)=Rr(2)
cFCc                {ii1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,10)-B23(1,2)-r10*Cij(4,10)-16*Cij(16,11)
cFC       Pr(2)=-B13(1,10)-r21*Cij(4,10)-4*Cij(15,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(4,11)=Rr(1)
cFC       Cij(5,11)=Rr(2)
cFCc                {ii1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,10)+B23(1,3)-r10*Cij(5,10)-14*Cij(17,11)
cFC       Pr(2)=-B13(1,10)-r21*Cij(5,10)-6*Cij(16,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(5,11)=Rr(1)
cFC       Cij(6,11)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,10)-B23(1,4)-r10*Cij(6,10)-12*Cij(18,11)
cFC       Pr(2)=-B13(1,10)-r21*Cij(6,10)-8*Cij(17,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(6,11)=Rr(1)
cFC       Cij(7,11)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,10)+B23(1,5)-r10*Cij(7,10)-10*Cij(19,11)
cFC       Pr(2)=-B13(1,10)-r21*Cij(7,10)-10*Cij(18,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(7,11)=Rr(1)
cFC       Cij(8,11)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,10)-B23(1,6)-r10*Cij(8,10)-8*Cij(20,11)
cFC       Pr(2)=-B13(1,10)-r21*Cij(8,10)-12*Cij(19,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(8,11)=Rr(1)
cFC       Cij(9,11)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(1,10)+B23(1,7)-r10*Cij(9,10)-6*Cij(21,11)
cFC       Pr(2)=-B13(1,10)-r21*Cij(9,10)-14*Cij(20,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(9,11)=Rr(1)
cFC       Cij(10,11)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(1,10)-B23(1,8)-r10*Cij(10,10)-4*Cij(22,11)
cFC       Pr(2)=-B13(1,10)-r21*Cij(10,10)-16*Cij(21,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(10,11)=Rr(1)
cFC       Cij(11,11)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(1,10)+B23(1,9)-r10*Cij(11,10)-2*Cij(14,11)
cFC       Pr(2)=-B13(1,10)-r21*Cij(11,10)-18*Cij(22,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(11,11)=Rr(1)
cFC       Cij(12,11)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(1,10)-B23(1,10)-r10*Cij(2,10)
cFC       Pr(2)=-B13(1,10)-r21*Cij(2,10)-20*Cij(14,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(12,11)=Rr(1)
cFC       Cij(2,11)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
cFC       Cij(49,12)=(-5*p1sq**5-5*p1sq**4*(p2sq+s12)-p1sq**3*(5*p2sq**2+8*
cFC     &  p2sq*s12+5*s12**2)-p1sq**2*(5*p2sq**3+9*p2sq**2*s12+9*p2sq*s12**
cFC     &  2+5*s12**3)-p1sq*(5*p2sq**4+8*p2sq**3*s12+9*p2sq**2*s12**2+8*p2s
cFC     &  q*s12**3+5*s12**4)-5*(p2sq**5+p2sq**4*s12+p2sq**3*s12**2+p2sq**2
cFC     &  *s12**3+p2sq*s12**4+s12**5))/7.6640256d9+(B23(6,10)+r10*Cij(41
cFC     &  ,11)+(r21)*Cij(42,11))/24.d0
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ii1, 1}
cFC       Pr(1)=B13(6,11)+B23(6,10)-r10*Cij(41,11)-2*Cij(49,12)
cFC       Pr(2)=B12(6,11)-B13(6,11)-r21*Cij(41,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFC       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFC       Cij(46,12)=Rr(1)
cFC       Cij(48,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ii1, 2}
cFC       Pr(1)=B13(6,11)-B23(6,11)-r10*Cij(42,11)
cFC       Pr(2)=-B13(6,11)-r21*Cij(42,11)-2*Cij(49,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(48,12)=Rr(1)
cFC       Cij(47,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, ii1, 1, 1, 1}
cFC       Pr(1)=B13(5,11)+B23(5,8)-r10*Cij(37,11)-6*Cij(46,12)
cFC       Pr(2)=B12(5,11)-B13(5,11)-r21*Cij(37,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFC       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFC       Cij(41,12)=Rr(1)
cFC       Cij(43,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, ii1, 2, 1, 1}
cFC       Pr(1)=B13(5,11)-B23(5,9)-r10*Cij(39,11)-4*Cij(48,12)
cFC       Pr(2)=-B13(5,11)-r21*Cij(39,11)-2*Cij(46,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(43,12)=Rr(1)
cFC       Cij(44,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, ii1, 2, 2, 1}
cFC       Pr(1)=B13(5,11)+B23(5,10)-r10*Cij(40,11)-2*Cij(47,12)
cFC       Pr(2)=-B13(5,11)-r21*Cij(40,11)-4*Cij(48,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(44,12)=Rr(1)
cFC       Cij(45,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, ii1, 2, 2, 2}
cFC       Pr(1)=B13(5,11)-B23(5,11)-r10*Cij(38,11)
cFC       Pr(2)=-B13(5,11)-r21*Cij(38,11)-6*Cij(47,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(45,12)=Rr(1)
cFC       Cij(42,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(4,11)+B23(4,6)-r10*Cij(31,11)-10*Cij(41,12)
cFC       Pr(2)=B12(4,11)-B13(4,11)-r21*Cij(31,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFC       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFC       Cij(34,12)=Rr(1)
cFC       Cij(36,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(4,11)-B23(4,7)-r10*Cij(33,11)-8*Cij(43,12)
cFC       Pr(2)=-B13(4,11)-r21*Cij(33,11)-2*Cij(41,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(36,12)=Rr(1)
cFC       Cij(37,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(4,11)+B23(4,8)-r10*Cij(34,11)-6*Cij(44,12)
cFC       Pr(2)=-B13(4,11)-r21*Cij(34,11)-4*Cij(43,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(37,12)=Rr(1)
cFC       Cij(38,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(4,11)-B23(4,9)-r10*Cij(35,11)-4*Cij(45,12)
cFC       Pr(2)=-B13(4,11)-r21*Cij(35,11)-6*Cij(44,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(38,12)=Rr(1)
cFC       Cij(39,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(4,11)+B23(4,10)-r10*Cij(36,11)-2*Cij(42,12)
cFC       Pr(2)=-B13(4,11)-r21*Cij(36,11)-8*Cij(45,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(39,12)=Rr(1)
cFC       Cij(40,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(4,11)-B23(4,11)-r10*Cij(32,11)
cFC       Pr(2)=-B13(4,11)-r21*Cij(32,11)-10*Cij(42,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(40,12)=Rr(1)
cFC       Cij(35,12)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(3,11)+B23(3,4)-r10*Cij(23,11)-14*Cij(34,12)
cFC       Pr(2)=B12(3,11)-B13(3,11)-r21*Cij(23,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFC       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFC       Cij(25,12)=Rr(1)
cFC       Cij(27,12)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(3,11)-B23(3,5)-r10*Cij(25,11)-12*Cij(36,12)
cFC       Pr(2)=-B13(3,11)-r21*Cij(25,11)-2*Cij(34,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(27,12)=Rr(1)
cFC       Cij(28,12)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(3,11)+B23(3,6)-r10*Cij(26,11)-10*Cij(37,12)
cFC       Pr(2)=-B13(3,11)-r21*Cij(26,11)-4*Cij(36,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(28,12)=Rr(1)
cFC       Cij(29,12)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(3,11)-B23(3,7)-r10*Cij(27,11)-8*Cij(38,12)
cFC       Pr(2)=-B13(3,11)-r21*Cij(27,11)-6*Cij(37,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(29,12)=Rr(1)
cFC       Cij(30,12)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(3,11)+B23(3,8)-r10*Cij(28,11)-6*Cij(39,12)
cFC       Pr(2)=-B13(3,11)-r21*Cij(28,11)-8*Cij(38,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(30,12)=Rr(1)
cFC       Cij(31,12)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(3,11)-B23(3,9)-r10*Cij(29,11)-4*Cij(40,12)
cFC       Pr(2)=-B13(3,11)-r21*Cij(29,11)-10*Cij(39,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(31,12)=Rr(1)
cFC       Cij(32,12)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(3,11)+B23(3,10)-r10*Cij(30,11)-2*Cij(35,12)
cFC       Pr(2)=-B13(3,11)-r21*Cij(30,11)-12*Cij(40,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(32,12)=Rr(1)
cFC       Cij(33,12)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(3,11)-B23(3,11)-r10*Cij(24,11)
cFC       Pr(2)=-B13(3,11)-r21*Cij(24,11)-14*Cij(35,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(33,12)=Rr(1)
cFC       Cij(26,12)=Rr(2)
cFCc                {0, 0, ii1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,11)+B23(2,2)-r10*Cij(13,11)-18*Cij(25,12)
cFC       Pr(2)=B12(2,11)-B13(2,11)-r21*Cij(13,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFC       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFC       Cij(14,12)=Rr(1)
cFC       Cij(16,12)=Rr(2)
cFCc                {0, 0, ii1, 2, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,11)-B23(2,3)-r10*Cij(15,11)-16*Cij(27,12)
cFC       Pr(2)=-B13(2,11)-r21*Cij(15,11)-2*Cij(25,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(16,12)=Rr(1)
cFC       Cij(17,12)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,11)+B23(2,4)-r10*Cij(16,11)-14*Cij(28,12)
cFC       Pr(2)=-B13(2,11)-r21*Cij(16,11)-4*Cij(27,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(17,12)=Rr(1)
cFC       Cij(18,12)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,11)-B23(2,5)-r10*Cij(17,11)-12*Cij(29,12)
cFC       Pr(2)=-B13(2,11)-r21*Cij(17,11)-6*Cij(28,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(18,12)=Rr(1)
cFC       Cij(19,12)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,11)+B23(2,6)-r10*Cij(18,11)-10*Cij(30,12)
cFC       Pr(2)=-B13(2,11)-r21*Cij(18,11)-8*Cij(29,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(19,12)=Rr(1)
cFC       Cij(20,12)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,11)-B23(2,7)-r10*Cij(19,11)-8*Cij(31,12)
cFC       Pr(2)=-B13(2,11)-r21*Cij(19,11)-10*Cij(30,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(20,12)=Rr(1)
cFC       Cij(21,12)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(2,11)+B23(2,8)-r10*Cij(20,11)-6*Cij(32,12)
cFC       Pr(2)=-B13(2,11)-r21*Cij(20,11)-12*Cij(31,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(21,12)=Rr(1)
cFC       Cij(22,12)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(2,11)-B23(2,9)-r10*Cij(21,11)-4*Cij(33,12)
cFC       Pr(2)=-B13(2,11)-r21*Cij(21,11)-14*Cij(32,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(22,12)=Rr(1)
cFC       Cij(23,12)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(2,11)+B23(2,10)-r10*Cij(22,11)-2*Cij(26,12)
cFC       Pr(2)=-B13(2,11)-r21*Cij(22,11)-16*Cij(33,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(23,12)=Rr(1)
cFC       Cij(24,12)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(2,11)-B23(2,11)-r10*Cij(14,11)
cFC       Pr(2)=-B13(2,11)-r21*Cij(14,11)-18*Cij(26,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(24,12)=Rr(1)
cFC       Cij(15,12)=Rr(2)
cFCc                {ii1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B023+B13(1,11)-r10*Cij(1,11)-22*Cij(14,12)
cFC       Pr(2)=B12(1,11)-B13(1,11)-r21*Cij(1,11)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFC       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFC       Cij(1,12)=Rr(1)
cFC       Cij(3,12)=Rr(2)
cFCc                {ii1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,11)-B23(1,1)-r10*Cij(3,11)-20*Cij(16,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(3,11)-2*Cij(14,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(3,12)=Rr(1)
cFC       Cij(4,12)=Rr(2)
cFCc                {ii1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,11)+B23(1,2)-r10*Cij(4,11)-18*Cij(17,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(4,11)-4*Cij(16,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(4,12)=Rr(1)
cFC       Cij(5,12)=Rr(2)
cFCc                {ii1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,11)-B23(1,3)-r10*Cij(5,11)-16*Cij(18,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(5,11)-6*Cij(17,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(5,12)=Rr(1)
cFC       Cij(6,12)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,11)+B23(1,4)-r10*Cij(6,11)-14*Cij(19,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(6,11)-8*Cij(18,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(6,12)=Rr(1)
cFC       Cij(7,12)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,11)-B23(1,5)-r10*Cij(7,11)-12*Cij(20,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(7,11)-10*Cij(19,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(7,12)=Rr(1)
cFC       Cij(8,12)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,11)+B23(1,6)-r10*Cij(8,11)-10*Cij(21,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(8,11)-12*Cij(20,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(8,12)=Rr(1)
cFC       Cij(9,12)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,11)-B23(1,7)-r10*Cij(9,11)-8*Cij(22,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(9,11)-14*Cij(21,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(9,12)=Rr(1)
cFC       Cij(10,12)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(1,11)+B23(1,8)-r10*Cij(10,11)-6*Cij(23,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(10,11)-16*Cij(22,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(10,12)=Rr(1)
cFC       Cij(11,12)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(1,11)-B23(1,9)-r10*Cij(11,11)-4*Cij(24,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(11,11)-18*Cij(23,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(11,12)=Rr(1)
cFC       Cij(12,12)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(1,11)+B23(1,10)-r10*Cij(12,11)-2*Cij(15,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(12,11)-20*Cij(24,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(12,12)=Rr(1)
cFC       Cij(13,12)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(1,11)-B23(1,11)-r10*Cij(2,11)
cFC       Pr(2)=-B13(1,11)-r21*Cij(2,11)-22*Cij(15,12)
cFC       Rr(2)=(Pr(2)-z21*Pr(1))*iz22
cFCc       Rr(1)=(Pr(1)-  z12*Rr(2))*iz11
cFCc       Cij(13,12)=Rr(1)
cFC       Cij(2,12)=Rr(2)


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

c  001,002
       Pr(1) = B13(2, 2) - B23(2, 2) - r10*Cij(4, 2)
       Pr(2) = B12(2, 2) - B13(2, 2) - r21*Cij(4, 2)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11 

       Cij(5,3)=Rr(1)
       Cij(6,3)=Rr(2)

c 111,211
       Pr(1) = -B023 + B13(1,2) - r10* Cij(1,2) - 4* Cij(5,3) 
       Pr(2) = B12(1,2) - B13(1,2) - r21* Cij(1,2)
   
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11 

       Cij(1,3) = Rr(1)
       Cij(3,3) = Rr(2)

c  121,221

       Pr(1)=B13(1, 2) + B23(1, 1) - r10* Cij(3, 2) - 2* Cij(6, 3)
       Pr(2)=-B13(1, 2) - r21* Cij(3, 2) - 2* Cij(5, 3)
       
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 

c       Cij(3,3) = Rr(1)
       Cij(4,3) = Rr(2)

c  122,222
       Pr(1)=B13(1, 2) - B23(1, 2) - r10* Cij(2, 2)
       Pr(2)=-B13(1, 2) - r21* Cij(2, 2) - 4* Cij(6, 3)
   
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 

c       Cij(1,1) = Rr(1)
       Cij(2,3) = Rr(2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 0000
       Cij(9,4)=(-p1sq - p2sq - s12)/192d0 + 
     &  (B23(2,2) + r10*Cij(5,3) +r21*Cij(6,3))/8d0
  

c 0011,0021
       Pr(1)=B13(2,3) + B23(2,2) - r10*Cij(5,3) - 2*Cij(9,4)
       Pr(2)=B12(2,3) - B13(2,3) - r21*Cij(5,3) 
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11 

       Cij(6,4)=Rr(1)   
       Cij(8,4)=Rr(2)
c0012,0022
       Pr(1)=B13(2,3) - B23(2,3) - r10*Cij(6,3)
       Pr(2)=-B13(2,3) - r21*Cij(6,3) - 2*Cij(9,4)  
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 

c      Cij(8,4)=Rr(1)
       Cij(7,4)=Rr(2)

c 1111,2111
       Pr(1)=B023 + B13(1,3) - r10*Cij(1,3) - 6*Cij(6,4)
       Pr(2)=B12(1,3) - B13(1,3) - r21*Cij(1,3) 
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11 

       Cij(1,4)=Rr(1)
       Cij(3,4)=Rr(2)

c 1211 2211
       Pr(1)=B13(1,3) - B23(1,1) - r10*Cij(3,3) - 4*Cij(8,4)
       Pr(2)=-B13(1,3) - r21*Cij(3,3) - 2*Cij(6,4) 
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 

c       Cij(3,4)=Rr(1)
       Cij(4,4)=Rr(2)

c 1221 2221
       Pr(1)=B13(1,3) + B23(1,2) - r10*Cij(4,3) - 2*Cij(7,4)
       Pr(2)=-B13(1,3) - r21*Cij(4,3) - 4*Cij(8,4) 
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c      Rr(1)=(PR(2)-z12*Rr(2))*iz11 

c       Cij(4,4)=Rr(1)
       Cij(5,4)=Rr(2)

c 1222,2222
       Pr(1)=B13(1,3) - B23(1,3) - r10*Cij(2,3)
       Pr(2)=-B13(1,3) - r21*Cij(2,3) - 6*Cij(7,4) 
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 

c       Cij(5,4)=Rr(1)
       Cij(2,4)=Rr(2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 00001,00002
       Pr(1)=B13(3,4) - B23(3,4) - r10*Cij(9,4)
       Pr(2)=B12(3,4) - B13(3,4) - r21*Cij(9,4)   
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11 

       Cij(11,5)=Rr(1)
       Cij(12,5)=Rr(2)

CCCCCC
c 00111,00211
       Pr(1)=B13(2,4) - B23(2,2) - r10*Cij(6,4) - 4*Cij(11,5)
       Pr(2)=B12(2,4) - B13(2,4) - r21*Cij(6,4)  
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11 

       Cij(7,5)=Rr(1)
       Cij(9,5)=Rr(2)
c 00121,00221
       Pr(1)=B13(2,4) + B23(2,3) - r10*Cij(8,4) - 2*Cij(12,5)
       Pr(2)=-B13(2,4) - r21*Cij(8,4) - 2*Cij(11,5)  
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 

c       Cij(9,5)=Rr(1)
       Cij(10,5)=Rr(2)
c 00122.00222
       Pr(1)=B13(2,4) - B23(2,4) - r10*Cij(7,4)
       Pr(2)=-B13(2,4) - r21*Cij(7,4) - 4*Cij(12,5)   
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 

c       Cij(10,5)=Rr(1)
       Cij(8,5)=Rr(2)
CCCCC

c 11111,21111
       Pr(1)=-B023 + B13(1,4) - r10*Cij(1,4) - 8*Cij(7,5)
       Pr(2)= B12(1,4) - B13(1,4) - r21*Cij(1,4)  
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11 

       Cij(1,5)=Rr(1)
       Cij(3,5)=Rr(2)
c 12111,22111
       Pr(1)=B13(1,4) + B23(1,1) - r10*Cij(3,4) - 6*Cij(9,5)
       Pr(2)=-B13(1,4) - r21*Cij(3,4) - 2*Cij(7,5) 
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 

c       Cij(3,5)=Rr(1)
       Cij(4,5)=Rr(2)
c 12211,22211
       Pr(1)=B13(1,4) - B23(1,2) - r10*Cij(4,4) - 4*Cij(10,5)
       Pr(2)=-B13(1,4) - r21*Cij(4,4) - 4*Cij(9,5)  
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 

c       Cij(4,5)=Rr(1)
       Cij(5,5)=Rr(2)
c 12221,22221
       Pr(1)=B13(1,4) + B23(1,3) - r10*Cij(5,4) - 2*Cij(8,5)
       Pr(2)=-B13(1,4) - r21*Cij(5,4) - 6*Cij(10,5)   
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 

c       Cij(5,5)=Rr(1)
       Cij(6,5)=Rr(2)
c 12222,22222
       Pr(1)=B13(1,4) - B23(1,4) - r10*Cij(2,4)
       Pr(2)=-B13(1,4) - r21*Cij(2,4) - 8*Cij(8,5)  
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 

c       Cij(6,5)=Rr(1)
       Cij(2,5)=Rr(2)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 000000
      
       Rr(1)=(B23(3,4)+r10*Cij(11,5)+r21*Cij(12,5))/12d0
     &  +(p1sq*p1sq+(p2sq+s12)*p1sq+p2sq*p2sq+s12*s12+p2sq*s12)/8640d0

       Cij(16,6)=Rr(1)
       

c 000011 000021
       Pr(1)=B13(3,5) + B23(3,4) - r10*Cij(11,5) - 2*Cij(16,6)
       Pr(2)=B12(3,5) - B13(3,5) - r21*Cij(11,5)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11 

       Cij(13,6)=Rr(1)
       Cij(15,6)=Rr(2)

c 000012 000022

       Pr(1)=B13(3,5) - B23(3,5) - r10*Cij(12,5)
       Pr(2)=-B13(3,5) - r21*Cij(12,5) - 2*Cij(16,6)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c       Cij(15,6)=Rr(1)
        Cij(14,6)=Rr(2)

c 001111 002111
        
       Pr(1)=B13(2,5) + B23(2,2) - r10*Cij(7,5) - 6*Cij(13,6)
       Pr(2)=B12(2,5) - B13(2,5) - r21*Cij(7,5)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
        Cij(8,6)=Rr(1)
        Cij(10,6)=Rr(2)


c 001211 002211

       Pr(1)=B13(2,5) - B23(2,3) - r10*Cij(9,5) - 4*Cij(15,6)
       Pr(2)=-B13(2,5) - r21*Cij(9,5) - 2*Cij(13,6)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(10,6)=Rr(1)
        Cij(11,6)=Rr(2)

c 001221 002221

       Pr(1)=B13(2,5) + B23(2,4) - r10*Cij(10,5) - 2*Cij(14,6)
       Pr(2)=-B13(2,5) - r21*Cij(10,5) - 4*Cij(15,6)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(11,6)=Rr(1)
        Cij(12,6)=Rr(2)

c 001222 002222
 
       Pr(1)=B13(2,5) - B23(2,5) - r10*Cij(8,5)
       Pr(2)=-B13(2,5) - r21*Cij(8,5) - 6*Cij(14,6)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c       Cij(12,6)=Rr(1)
        Cij(9,6)=Rr(2)


c 111111 211111

       Pr(1)=B023 + B13(1,5) - r10*Cij(1,5) - 10*Cij(8,6)
       Pr(2)=B12(1,5) - B13(1,5) - r21*Cij(1,5)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
        Cij(1,6)=Rr(1)
        Cij(3,6)=Rr(2)

c 121111 221111

       Pr(1)=B13(1,5) - B23(1,1) - r10*Cij(3,5) - 8*Cij(10,6)
       Pr(2)=-B13(1,5) - r21*Cij(3,5) - 2*Cij(8,6)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(3,6)=Rr(1)
        Cij(4,6)=Rr(2)


c 122111 222111

       Pr(1)=B13(1,5) + B23(1,2) - r10*Cij(4,5) - 6*Cij(11,6)
       Pr(2)=-B13(1,5) - r21*Cij(4,5) - 4*Cij(10,6)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(4,6)=Rr(1)
        Cij(5,6)=Rr(2)


c 122211  222211

       Pr(1)=B13(1,5) - B23(1,3) - r10*Cij(5,5) - 4*Cij(12,6)
       Pr(2)=-B13(1,5) - r21*Cij(5,5) - 6*Cij(11,6)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(5,6)=Rr(1)
        Cij(6,6)=Rr(2)

c 122221  222221

       Pr(1)=B13(1,5) + B23(1,4) - r10*Cij(6,5) - 2*Cij(9,6)
       Pr(2)=-B13(1,5) - r21*Cij(6,5) - 8*Cij(12,6)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(6,6)=Rr(1)
        Cij(7,6)=Rr(2)

c 122222  222222

       Pr(1)=B13(1,5) - B23(1,5) - r10*Cij(2,5)
       Pr(2)=-B13(1,5) - r21*Cij(2,5) - 10*Cij(9,6)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(7,6)=Rr(1)
        Cij(2,6)=Rr(2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 0000001 0000002

       Pr(1)=B13(4,6) - B23(4,6) - r10*Cij(16,6)
       Pr(2)=B12(4,6) - B13(4,6) - r21*Cij(16,6)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
        Cij(19,7)=Rr(1)
        Cij(20,7)=Rr(2)
c 0000111 0000211
       Pr(1)=B13(3,6) - B23(3,4) - r10*Cij(13,6) - 4*Cij(19,7)
       Pr(2)=B12(3,6) - B13(3,6) - r21*Cij(13,6)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
        Cij(15,7)=Rr(1)
        Cij(17,7)=Rr(2)
c 0000121 0000221
       Pr(1)=B13(3,6) + B23(3,5) - r10*Cij(15,6) - 2*Cij(20,7)
       Pr(2)=-B13(3,6) - r21*Cij(15,6) - 2*Cij(19,7)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(17,7)=Rr(1)
        Cij(18,7)=Rr(2)
c 0000122 0000222
       Pr(1)=B13(3,6) - B23(3,6) - r10*Cij(14,6)
       Pr(2)=-B13(3,6) - r21*Cij(14,6) - 4*Cij(20,7)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(18,7)=Rr(1)
        Cij(16,7)=Rr(2)
c 0011111 0021111
       Pr(1)=B13(2,6) - B23(2,2) - r10*Cij(8,6) - 8*Cij(15,7)
       Pr(2)=B12(2,6) - B13(2,6) - r21*Cij(8,6)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
        Cij(9,7)=Rr(1)
        Cij(11,7)=Rr(2)

c 0012111 0022111
       Pr(1)=B13(2,6) + B23(2,3) - r10*Cij(10,6) - 6*Cij(17,7)
       Pr(2)=-B13(2,6) - r21*Cij(10,6) - 2*Cij(15,7)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(11,7)=Rr(1)
        Cij(12,7)=Rr(2)

c 0012211 0022211
       Pr(1)=B13(2,6) - B23(2,4) - r10*Cij(11,6) - 4*Cij(18,7)
       Pr(2)=-B13(2,6) - r21*Cij(11,6) - 4*Cij(17,7)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(12,7)=Rr(1)
        Cij(13,7)=Rr(2)
c 0012221 0022221

       Pr(1)=B13(2,6) + B23(2,5) - r10*Cij(12,6) - 2*Cij(16,7)
       Pr(2)=-B13(2,6) - r21*Cij(12,6) - 6*Cij(18,7)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(13,7)=Rr(1)
        Cij(14,7)=Rr(2)
c 0012222 0022222
       Pr(1)=B13(2,6) - B23(2,6) - r10*Cij(9,6)
       Pr(2)=-B13(2,6) - r21*Cij(9,6) - 8*Cij(16,7)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(14,7)=Rr(1)
        Cij(10,7)=Rr(2)

c 1111111 2111111
       Pr(1)=-B023 + B13(1,6) - r10*Cij(1,6) - 12*Cij(9,7)
       Pr(2)=B12(1,6) - B13(1,6) - r21*Cij(1,6)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
        Cij(1,7)=Rr(1)
        Cij(3,7)=Rr(2)

c 1211111 2211111
       Pr(1)=B13(1,6) + B23(1,1) - r10*Cij(3,6) - 10*Cij(11,7)
       Pr(2)=-B13(1,6) - r21*Cij(3,6) - 2*Cij(9,7)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(3,7)=Rr(1)
        Cij(4,7)=Rr(2)

c 1221111 2221111
       Pr(1)=B13(1,6) - B23(1,2) - r10*Cij(4,6) - 8*Cij(12,7)
       Pr(2)=-B13(1,6) - r21*Cij(4,6) - 4*Cij(11,7)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(4,7)=Rr(1)
        Cij(5,7)=Rr(2)


c 1222111 2222111
       Pr(1)=B13(1,6) + B23(1,3) - r10*Cij(5,6) - 6*Cij(13,7)
       Pr(2)=-B13(1,6) - r21*Cij(5,6) - 6*Cij(12,7)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(5,7)=Rr(1)
        Cij(6,7)=Rr(2)


c 1222211 2222211
       Pr(1)=B13(1,6) - B23(1,4) - r10*Cij(6,6) - 4*Cij(14,7)
       Pr(2)=-B13(1,6) - r21*Cij(6,6) - 8*Cij(13,7)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(6,7)=Rr(1)
        Cij(7,7)=Rr(2)
c 1222221 2222221
       Pr(1)=B13(1,6) + B23(1,5) - r10*Cij(7,6) - 2*Cij(10,7)
       Pr(2)=-B13(1,6) - r21*Cij(7,6) - 10*Cij(14,7)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(7,7)=Rr(1)
        Cij(8,7)=Rr(2)
c 1222222 2222222 
       Pr(1)=B13(1,6) - B23(1,6) - r10*Cij(2,6)
       Pr(2)=-B13(1,6) - r21*Cij(2,6) - 12*Cij(10,7)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(8,7)=Rr(1)
        Cij(2,7)=Rr(2)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 00000000

        Rr(1)=(-3*p1sq*p1sq*p1sq - 3*p1sq*p1sq*(p2sq + s12) - 
     &     3*(p2sq + s12)*(p2sq*p2sq + s12*s12) - 
     &     p1sq*(3*p2sq*p2sq + 4*p2sq*s12 + 3*s12*s12))/1.29024d6 + 
     &  (B23(4,6) + r10*Cij(19,7) + r21*Cij(20,7))/16d0

        Cij(25,8)=Rr(1)

c 00000011 00000021

       Pr(1)=B13(4,7) + B23(4,6) - r10*Cij(19,7) - 2*Cij(25,8)
       Pr(2)=B12(4,7) - B13(4,7) - r21*Cij(19,7)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
        Cij(22,8)=Rr(1)
        Cij(24,8)=Rr(2)

c 00000012 00000022
       Pr(1)=B13(4,7) - B23(4,7) - r10*Cij(20,7)
       Pr(2)=-B13(4,7) - r21*Cij(20,7) - 2*Cij(25,8)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(24,8)=Rr(1)
        Cij(23,8)=Rr(2)


c 00001111 00002111

       Pr(1)=B13(3,7) + B23(3,4) - r10*Cij(15,7) - 6*Cij(22,8)
       Pr(2)=B12(3,7) - B13(3,7) - r21*Cij(15,7)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
        Cij(17,8)=Rr(1)
        Cij(19,8)=Rr(2)

c 00001211 00002211

       Pr(1)=B13(3,7) - B23(3,5) - r10*Cij(17,7) - 4*Cij(24,8)
       Pr(2)=-B13(3,7) - r21*Cij(17,7) - 2*Cij(22,8)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(19,8)=Rr(1)
        Cij(20,8)=Rr(2)

c 00001221 00002221
       Pr(1)=B13(3,7) + B23(3,6) - r10*Cij(18,7) - 2*Cij(23,8)
       Pr(2)=-B13(3,7) - r21*Cij(18,7) - 4*Cij(24,8)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(20,8)=Rr(1)
        Cij(21,8)=Rr(2)

c 00001222 00002222
       Pr(1)=B13(3,7) - B23(3,7) - r10*Cij(16,7)
       Pr(2)=-B13(3,7) - r21*Cij(16,7) - 6*Cij(23,8)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(21,8)=Rr(1)
        Cij(18,8)=Rr(2)


c 00111111 00211111

       Pr(1)=B13(2,7) + B23(2,2) - r10*Cij(9,7) - 10*Cij(17,8)
       Pr(2)=B12(2,7) - B13(2,7) - r21*Cij(9,7)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
        Cij(10,8)=Rr(1)
        Cij(12,8)=Rr(2)

c 00121111 00221111

       Pr(1)=B13(2,7) - B23(2,3) - r10*Cij(11,7) - 8*Cij(19,8)
       Pr(2)=-B13(2,7) - r21*Cij(11,7) - 2*Cij(17,8)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(12,8)=Rr(1)
        Cij(13,8)=Rr(2)

c 00122111 00222111
       Pr(1)=B13(2,7) + B23(2,4) - r10*Cij(12,7) - 6*Cij(20,8)
       Pr(2)=-B13(2,7) - r21*Cij(12,7) - 4*Cij(19,8)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(13,8)=Rr(1)
        Cij(14,8)=Rr(2)

c 00122211 00222211
       Pr(1)=B13(2,7) - B23(2,5) - r10*Cij(13,7) - 4*Cij(21,8)
       Pr(2)=-B13(2,7) - r21*Cij(13,7) - 6*Cij(20,8)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(14,8)=Rr(1)
        Cij(15,8)=Rr(2)

c 00122221  00222221
       Pr(1)=B13(2,7) + B23(2,6) - r10*Cij(14,7) - 2*Cij(18,8)
       Pr(2)=-B13(2,7) - r21*Cij(14,7) - 8*Cij(21,8)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(15,8)=Rr(1)
        Cij(16,8)=Rr(2)
c 00122222  00222222
       Pr(1)=B13(2,7) - B23(2,7) - r10*Cij(10,7)
       Pr(2)=-B13(2,7) - r21*Cij(10,7) - 10*Cij(18,8)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(16,8)=Rr(1)
        Cij(11,8)=Rr(2)

c 11111111 21111111

       Pr(1)=B023 + B13(1,7) - r10*Cij(1,7) - 14*Cij(10,8)
       Pr(2)=B12(1,7) - B13(1,7) - r21*Cij(1,7)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
        Cij(1,8)=Rr(1)
        Cij(3,8)=Rr(2)

c 12111111 22111111

       Pr(1)=B13(1,7) - B23(1,1) - r10*Cij(3,7) - 12*Cij(12,8)
       Pr(2)=-B13(1,7) - r21*Cij(3,7) - 2*Cij(10,8)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(3,8)=Rr(1)
        Cij(4,8)=Rr(2)


c 12211111 22211111
       Pr(1)=B13(1,7) + B23(1,2) - r10*Cij(4,7) - 10*Cij(13,8)
       Pr(2)=-B13(1,7) - r21*Cij(4,7) - 4*Cij(12,8)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(4,8)=Rr(1)
        Cij(5,8)=Rr(2)


c 12221111 22221111
       Pr(1)=B13(1,7) - B23(1,3) - r10*Cij(5,7) - 8*Cij(14,8)
       Pr(2)=-B13(1,7) - r21*Cij(5,7) - 6*Cij(13,8)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(5,8)=Rr(1)
        Cij(6,8)=Rr(2)


c 12222111  22222111

       Pr(1)=B13(1,7) + B23(1,4) - r10*Cij(6,7) - 6*Cij(15,8)
       Pr(2)=-B13(1,7) - r21*Cij(6,7) - 8*Cij(14,8)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(6,8)=Rr(1)
        Cij(7,8)=Rr(2)

c 12222211  22222211
       Pr(1)=B13(1,7) - B23(1,5) - r10*Cij(7,7) - 4*Cij(16,8)
       Pr(2)=-B13(1,7) - r21*Cij(7,7) - 10*Cij(15,8)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(7,8)=Rr(1)
        Cij(8,8)=Rr(2)


c 12222221  22222221

       Pr(1)=B13(1,7) + B23(1,6) - r10*Cij(8,7) - 2*Cij(11,8)
       Pr(2)=-B13(1,7) - r21*Cij(8,7) - 12*Cij(16,8)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(8,8)=Rr(1)
        Cij(9,8)=Rr(2)


c 12222222  22222222

       Pr(1)=B13(1,7) - B23(1,7) - r10*Cij(2,7)
       Pr(2)=-B13(1,7) - r21*Cij(2,7) - 14*Cij(11,8)

       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11 
c        Cij(9,8)=Rr(1)
        Cij(2,8)=Rr(2)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c 000000001 000000002
       Pr(1)=B13(5,8) - B23(5,8) - r10*Cij(25,8)
       Pr(2)=B12(5,8) - B13(5,8) - r21*Cij(25,8)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11
       Cij(29,9)=Rr(1)
       Cij(30,9)=Rr(2)


c 000000111 000000211
       Pr(1)=B13(4,8) - B23(4,6) - r10*Cij(22,8) - 4*Cij(29,9)
       Pr(2)=B12(4,8) - B13(4,8) - r21*Cij(22,8)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11
       Cij(25,9)=Rr(1)
       Cij(27,9)=Rr(2)

c 000000121 000000221
       Pr(1)=B13(4,8) + B23(4,7) - r10*Cij(24,8) - 2*Cij(30,9)
       Pr(2)=-B13(4,8) - r21*Cij(24,8) - 2*Cij(29,9)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11
c       Cij(27,9)=Rr(1)
       Cij(28,9)=Rr(2)


c 000000122 000000222
       Pr(1)=B13(4,8) - B23(4,8) - r10*Cij(23,8)
       Pr(2)=-B13(4,8) - r21*Cij(23,8) - 4*Cij(30,9)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11
c       Cij(28,9)=Rr(1)
       Cij(26,9)=Rr(2)

c 000011111 000021111
       Pr(1)=B13(3,8) - B23(3,4) - r10*Cij(17,8) - 8*Cij(25,9)
       Pr(2)=B12(3,8) - B13(3,8) - r21*Cij(17,8)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11
       Cij(19,9)=Rr(1)
       Cij(21,9)=Rr(2)

c 000012111 000022111
       Pr(1)=B13(3,8) + B23(3,5) - r10*Cij(19,8) - 6*Cij(27,9)
       Pr(2)=-B13(3,8) - r21*Cij(19,8) - 2*Cij(25,9)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11
c       Cij(21,9)=Rr(1)
       Cij(22,9)=Rr(2)

c 000012211 000022211
       Pr(1)=B13(3,8) - B23(3,6) - r10*Cij(20,8) - 4*Cij(28,9)
       Pr(2)=-B13(3,8) - r21*Cij(20,8) - 4*Cij(27,9)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11
c       Cij(22,9)=Rr(1)
       Cij(23,9)=Rr(2)

c 000012221 000022221
       Pr(1)=B13(3,8) + B23(3,7) - r10*Cij(21,8) - 2*Cij(26,9)
       Pr(2)=-B13(3,8) - r21*Cij(21,8) - 6*Cij(28,9)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11
c       Cij(23,9)=Rr(1)
       Cij(24,9)=Rr(2)

c 000012222 000022222
       Pr(1)=B13(3,8) - B23(3,8) - r10*Cij(18,8)
       Pr(2)=-B13(3,8) - r21*Cij(18,8) - 8*Cij(26,9)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11
c       Cij(24,9)=Rr(1)
       Cij(20,9)=Rr(2)


c 001111111 002111111
       Pr(1)=B13(2,8) - B23(2,2) - r10*Cij(10,8) - 12*Cij(19,9)
       Pr(2)=B12(2,8) - B13(2,8) - r21*Cij(10,8)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11
       Cij(11,9)=Rr(1)
       Cij(13,9)=Rr(2)


c 001211111 002211111
       Pr(1)=B13(2,8) + B23(2,3) - r10*Cij(12,8) - 10*Cij(21,9)
       Pr(2)=-B13(2,8) - r21*Cij(12,8) - 2*Cij(19,9)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11
c       Cij(13,9)=Rr(1)
       Cij(14,9)=Rr(2)
         
         
c 001221111 002221111
       Pr(1)=B13(2,8) - B23(2,4) - r10*Cij(13,8) - 8*Cij(22,9)
       Pr(2)=-B13(2,8) - r21*Cij(13,8) - 4*Cij(21,9)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11
c       Cij(14,9)=Rr(1)
       Cij(15,9)=Rr(2)
         
c 001222111 002222111
       Pr(1)=B13(2,8) + B23(2,5) - r10*Cij(14,8) - 6*Cij(23,9)
       Pr(2)=-B13(2,8) - r21*Cij(14,8) - 6*Cij(22,9)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11
c       Cij(15,9)=Rr(1)
       Cij(16,9)=Rr(2)
         
c 001222211 002222211
       Pr(1)=B13(2,8) - B23(2,6) - r10*Cij(15,8) - 4*Cij(24,9)
       Pr(2)=-B13(2,8) - r21*Cij(15,8) - 8*Cij(23,9)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11
c       Cij(16,9)=Rr(1)
       Cij(17,9)=Rr(2)

c 001222221 002222221
       Pr(1)=B13(2,8) + B23(2,7) - r10*Cij(16,8) - 2*Cij(20,9)
       Pr(2)=-B13(2,8) - r21*Cij(16,8) - 10*Cij(24,9)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11
c       Cij(17,9)=Rr(1)
       Cij(18,9)=Rr(2)

c 001222222 002222222
       Pr(1)=B13(2,8) - B23(2,8) - r10*Cij(11,8)
       Pr(2)=-B13(2,8) - r21*Cij(11,8) - 12*Cij(20,9)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11
c       Cij(18,9)=Rr(1)
       Cij(12,9)=Rr(2)

c 111111111 211111111
       Pr(1)=-B023 + B13(1,8) - r10*Cij(1,8) - 16*Cij(11,9)
       Pr(2)=B12(1,8) - B13(1,8) - r21*Cij(1,8)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
       Rr(1)=(PR(2)-z12*Rr(2))*iz11
       Cij(1,9)=Rr(1)
       Cij(3,9)=Rr(2)
c 121111111 221111111
       Pr(1)=B13(1,8) + B23(1,1) - r10*Cij(3,8) - 14*Cij(13,9)
       Pr(2)=-B13(1,8) - r21*Cij(3,8) - 2*Cij(11,9)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11
c       Cij(3,9)=Rr(1)
       Cij(4,9)=Rr(2)       
       
c 122111111 222111111
       Pr(1)=B13(1,8) - B23(1,2) - r10*Cij(4,8) - 12*Cij(14,9)
       Pr(2)=-B13(1,8) - r21*Cij(4,8) - 4*Cij(13,9)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11
c       Cij(4,9)=Rr(1)
       Cij(5,9)=Rr(2)       
c 122211111 222211111
       Pr(1)=B13(1,8) + B23(1,3) - r10*Cij(5,8) - 10*Cij(15,9)
       Pr(2)=-B13(1,8) - r21*Cij(5,8) - 6*Cij(14,9)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11
c       Cij(5,9)=Rr(1)
       Cij(6,9)=Rr(2)       
c 122221111 222221111
       Pr(1)=B13(1,8) - B23(1,4) - r10*Cij(6,8) - 8*Cij(16,9)
       Pr(2)=-B13(1,8) - r21*Cij(6,8) - 8*Cij(15,9)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11
c       Cij(6,9)=Rr(1)
       Cij(7,9)=Rr(2)
c 122222111 222222111
       Pr(1)=B13(1,8) + B23(1,5) - r10*Cij(7,8) - 6*Cij(17,9)
       Pr(2)=-B13(1,8) - r21*Cij(7,8) - 10*Cij(16,9)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11
c       Cij(7,9)=Rr(1)
       Cij(8,9)=Rr(2)
c 122222211 222222211
       Pr(1)=B13(1,8) - B23(1,6) - r10*Cij(8,8) - 4*Cij(18,9)
       Pr(2)=-B13(1,8) - r21*Cij(8,8) - 12*Cij(17,9)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11
c       Cij(8,9)=Rr(1)
       Cij(9,9)=Rr(2)
c 122222221 222222221
       Pr(1)=B13(1,8) + B23(1,7) - r10*Cij(9,8) - 2*Cij(12,9)
       Pr(2)=-B13(1,8) - r21*Cij(9,8) - 14*Cij(18,9)
       Rr(2)=(PR(1)-z21*PR(2))*iz22 
c       Rr(1)=(PR(2)-z12*Rr(2))*iz11
c       Cij(9,9)=Rr(1)
       Cij(10,9)=Rr(2)
cFCc 122222222 222222222
cFC       Pr(1)=B13(1,8) - B23(1,8) - r10*Cij(2,8)
cFC       Pr(2)=-B13(1,8) - r21*Cij(2,8) - 16*Cij(12,9)
cFC       Rr(2)=(PR(1)-z21*PR(2))*iz22 
cFCc       Rr(1)=(PR(2)-z12*Rr(2))*iz11
cFCc       Cij(10,9)=Rr(1)
cFC       Cij(2,9)=Rr(2)
cFC
cFC
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
cFC       Cij(36,10)=(2*p1sq**4+2*p1sq**3*(p2sq+s12)+p1sq**2*(2*p2sq**2+3*p
cFC     &  2sq*s12+2*s12**2)+p1sq*(2*p2sq**3+3*p2sq**2*s12+3*p2sq*s12**2+2*
cFC     &  s12**3)+2*(p2sq**4+p2sq**3*s12+p2sq**2*s12**2+p2sq*s12**3+s12**4
cFC     &  ))/4.8384d7+(B23(5,8)+r10*Cij(29,9)+(r21)*Cij(30,9))/20.d0
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, ii1, 1}
cFC       Pr(1)=B13(5,9)+B23(5,8)-r10*Cij(29,9)-2*Cij(36,10)
cFC       Pr(2)=B12(5,9)-B13(5,9)-r21*Cij(29,9)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFC       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFC       Cij(33,10)=Rr(1)
cFC       Cij(35,10)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, ii1, 2}
cFC       Pr(1)=B13(5,9)-B23(5,9)-r10*Cij(30,9)
cFC       Pr(2)=-B13(5,9)-r21*Cij(30,9)-2*Cij(36,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(35,10)=Rr(1)
cFC       Cij(34,10)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 1, 1, 1}
cFC       Pr(1)=B13(4,9)+B23(4,6)-r10*Cij(25,9)-6*Cij(33,10)
cFC       Pr(2)=B12(4,9)-B13(4,9)-r21*Cij(25,9)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFC       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFC       Cij(28,10)=Rr(1)
cFC       Cij(30,10)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 1, 1}
cFC       Pr(1)=B13(4,9)-B23(4,7)-r10*Cij(27,9)-4*Cij(35,10)
cFC       Pr(2)=-B13(4,9)-r21*Cij(27,9)-2*Cij(33,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(30,10)=Rr(1)
cFC       Cij(31,10)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 2, 1}
cFC       Pr(1)=B13(4,9)+B23(4,8)-r10*Cij(28,9)-2*Cij(34,10)
cFC       Pr(2)=-B13(4,9)-r21*Cij(28,9)-4*Cij(35,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(31,10)=Rr(1)
cFC       Cij(32,10)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 2, 2}
cFC       Pr(1)=B13(4,9)-B23(4,9)-r10*Cij(26,9)
cFC       Pr(2)=-B13(4,9)-r21*Cij(26,9)-6*Cij(34,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(32,10)=Rr(1)
cFC       Cij(29,10)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(3,9)+B23(3,4)-r10*Cij(19,9)-10*Cij(28,10)
cFC       Pr(2)=B12(3,9)-B13(3,9)-r21*Cij(19,9)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFC       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFC       Cij(21,10)=Rr(1)
cFC       Cij(23,10)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(3,9)-B23(3,5)-r10*Cij(21,9)-8*Cij(30,10)
cFC       Pr(2)=-B13(3,9)-r21*Cij(21,9)-2*Cij(28,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(23,10)=Rr(1)
cFC       Cij(24,10)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(3,9)+B23(3,6)-r10*Cij(22,9)-6*Cij(31,10)
cFC       Pr(2)=-B13(3,9)-r21*Cij(22,9)-4*Cij(30,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(24,10)=Rr(1)
cFC       Cij(25,10)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(3,9)-B23(3,7)-r10*Cij(23,9)-4*Cij(32,10)
cFC       Pr(2)=-B13(3,9)-r21*Cij(23,9)-6*Cij(31,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(25,10)=Rr(1)
cFC       Cij(26,10)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(3,9)+B23(3,8)-r10*Cij(24,9)-2*Cij(29,10)
cFC       Pr(2)=-B13(3,9)-r21*Cij(24,9)-8*Cij(32,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(26,10)=Rr(1)
cFC       Cij(27,10)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(3,9)-B23(3,9)-r10*Cij(20,9)
cFC       Pr(2)=-B13(3,9)-r21*Cij(20,9)-10*Cij(29,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(27,10)=Rr(1)
cFC       Cij(22,10)=Rr(2)
cFCc                {0, 0, ii1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,9)+B23(2,2)-r10*Cij(11,9)-14*Cij(21,10)
cFC       Pr(2)=B12(2,9)-B13(2,9)-r21*Cij(11,9)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFC       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFC       Cij(12,10)=Rr(1)
cFC       Cij(14,10)=Rr(2)
cFCc                {0, 0, ii1, 2, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,9)-B23(2,3)-r10*Cij(13,9)-12*Cij(23,10)
cFC       Pr(2)=-B13(2,9)-r21*Cij(13,9)-2*Cij(21,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(14,10)=Rr(1)
cFC       Cij(15,10)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,9)+B23(2,4)-r10*Cij(14,9)-10*Cij(24,10)
cFC       Pr(2)=-B13(2,9)-r21*Cij(14,9)-4*Cij(23,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(15,10)=Rr(1)
cFC       Cij(16,10)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,9)-B23(2,5)-r10*Cij(15,9)-8*Cij(25,10)
cFC       Pr(2)=-B13(2,9)-r21*Cij(15,9)-6*Cij(24,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(16,10)=Rr(1)
cFC       Cij(17,10)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(2,9)+B23(2,6)-r10*Cij(16,9)-6*Cij(26,10)
cFC       Pr(2)=-B13(2,9)-r21*Cij(16,9)-8*Cij(25,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(17,10)=Rr(1)
cFC       Cij(18,10)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(2,9)-B23(2,7)-r10*Cij(17,9)-4*Cij(27,10)
cFC       Pr(2)=-B13(2,9)-r21*Cij(17,9)-10*Cij(26,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(18,10)=Rr(1)
cFC       Cij(19,10)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(2,9)+B23(2,8)-r10*Cij(18,9)-2*Cij(22,10)
cFC       Pr(2)=-B13(2,9)-r21*Cij(18,9)-12*Cij(27,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(19,10)=Rr(1)
cFC       Cij(20,10)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(2,9)-B23(2,9)-r10*Cij(12,9)
cFC       Pr(2)=-B13(2,9)-r21*Cij(12,9)-14*Cij(22,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(20,10)=Rr(1)
cFC       Cij(13,10)=Rr(2)
cFCc                {ii1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B023+B13(1,9)-r10*Cij(1,9)-18*Cij(12,10)
cFC       Pr(2)=B12(1,9)-B13(1,9)-r21*Cij(1,9)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFC       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFC       Cij(1,10)=Rr(1)
cFC       Cij(3,10)=Rr(2)
cFCc                {ii1, 2, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,9)-B23(1,1)-r10*Cij(3,9)-16*Cij(14,10)
cFC       Pr(2)=-B13(1,9)-r21*Cij(3,9)-2*Cij(12,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(3,10)=Rr(1)
cFC       Cij(4,10)=Rr(2)
cFCc                {ii1, 2, 2, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,9)+B23(1,2)-r10*Cij(4,9)-14*Cij(15,10)
cFC       Pr(2)=-B13(1,9)-r21*Cij(4,9)-4*Cij(14,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(4,10)=Rr(1)
cFC       Cij(5,10)=Rr(2)
cFCc                {ii1, 2, 2, 2, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,9)-B23(1,3)-r10*Cij(5,9)-12*Cij(16,10)
cFC       Pr(2)=-B13(1,9)-r21*Cij(5,9)-6*Cij(15,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(5,10)=Rr(1)
cFC       Cij(6,10)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,9)+B23(1,4)-r10*Cij(6,9)-10*Cij(17,10)
cFC       Pr(2)=-B13(1,9)-r21*Cij(6,9)-8*Cij(16,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(6,10)=Rr(1)
cFC       Cij(7,10)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,9)-B23(1,5)-r10*Cij(7,9)-8*Cij(18,10)
cFC       Pr(2)=-B13(1,9)-r21*Cij(7,9)-10*Cij(17,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(7,10)=Rr(1)
cFC       Cij(8,10)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(1,9)+B23(1,6)-r10*Cij(8,9)-6*Cij(19,10)
cFC       Pr(2)=-B13(1,9)-r21*Cij(8,9)-12*Cij(18,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(8,10)=Rr(1)
cFC       Cij(9,10)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(1,9)-B23(1,7)-r10*Cij(9,9)-4*Cij(20,10)
cFC       Pr(2)=-B13(1,9)-r21*Cij(9,9)-14*Cij(19,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(9,10)=Rr(1)
cFC       Cij(10,10)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(1,9)+B23(1,8)-r10*Cij(10,9)-2*Cij(13,10)
cFC       Pr(2)=-B13(1,9)-r21*Cij(10,9)-16*Cij(20,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(10,10)=Rr(1)
cFC       Cij(11,10)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(1,9)-B23(1,9)-r10*Cij(2,9)
cFC       Pr(2)=-B13(1,9)-r21*Cij(2,9)-18*Cij(13,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(11,10)=Rr(1)
cFC       Cij(2,10)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ii1}
cFC       Pr(1)=B13(6,10)-B23(6,10)-r10*Cij(36,10)
cFC       Pr(2)=B12(6,10)-B13(6,10)-r21*Cij(36,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFC       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFC       Cij(41,11)=Rr(1)
cFC       Cij(42,11)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, ii1, 1, 1}
cFC       Pr(1)=B13(5,10)-B23(5,8)-r10*Cij(33,10)-4*Cij(41,11)
cFC       Pr(2)=B12(5,10)-B13(5,10)-r21*Cij(33,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFC       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFC       Cij(37,11)=Rr(1)
cFC       Cij(39,11)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, ii1, 2, 1}
cFC       Pr(1)=B13(5,10)+B23(5,9)-r10*Cij(35,10)-2*Cij(42,11)
cFC       Pr(2)=-B13(5,10)-r21*Cij(35,10)-2*Cij(41,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(39,11)=Rr(1)
cFC       Cij(40,11)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, ii1, 2, 2}
cFC       Pr(1)=B13(5,10)-B23(5,10)-r10*Cij(34,10)
cFC       Pr(2)=-B13(5,10)-r21*Cij(34,10)-4*Cij(42,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(40,11)=Rr(1)
cFC       Cij(38,11)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 1, 1, 1, 1}
cFC       Pr(1)=B13(4,10)-B23(4,6)-r10*Cij(28,10)-8*Cij(37,11)
cFC       Pr(2)=B12(4,10)-B13(4,10)-r21*Cij(28,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFC       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFC       Cij(31,11)=Rr(1)
cFC       Cij(33,11)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 1, 1, 1}
cFC       Pr(1)=B13(4,10)+B23(4,7)-r10*Cij(30,10)-6*Cij(39,11)
cFC       Pr(2)=-B13(4,10)-r21*Cij(30,10)-2*Cij(37,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(33,11)=Rr(1)
cFC       Cij(34,11)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 2, 1, 1}
cFC       Pr(1)=B13(4,10)-B23(4,8)-r10*Cij(31,10)-4*Cij(40,11)
cFC       Pr(2)=-B13(4,10)-r21*Cij(31,10)-4*Cij(39,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(34,11)=Rr(1)
cFC       Cij(35,11)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 2, 2, 1}
cFC       Pr(1)=B13(4,10)+B23(4,9)-r10*Cij(32,10)-2*Cij(38,11)
cFC       Pr(2)=-B13(4,10)-r21*Cij(32,10)-6*Cij(40,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(35,11)=Rr(1)
cFC       Cij(36,11)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 2, 2, 2}
cFC       Pr(1)=B13(4,10)-B23(4,10)-r10*Cij(29,10)
cFC       Pr(2)=-B13(4,10)-r21*Cij(29,10)-8*Cij(38,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(36,11)=Rr(1)
cFC       Cij(32,11)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(3,10)-B23(3,4)-r10*Cij(21,10)-12*Cij(31,11)
cFC       Pr(2)=B12(3,10)-B13(3,10)-r21*Cij(21,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFC       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFC       Cij(23,11)=Rr(1)
cFC       Cij(25,11)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(3,10)+B23(3,5)-r10*Cij(23,10)-10*Cij(33,11)
cFC       Pr(2)=-B13(3,10)-r21*Cij(23,10)-2*Cij(31,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(25,11)=Rr(1)
cFC       Cij(26,11)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(3,10)-B23(3,6)-r10*Cij(24,10)-8*Cij(34,11)
cFC       Pr(2)=-B13(3,10)-r21*Cij(24,10)-4*Cij(33,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(26,11)=Rr(1)
cFC       Cij(27,11)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(3,10)+B23(3,7)-r10*Cij(25,10)-6*Cij(35,11)
cFC       Pr(2)=-B13(3,10)-r21*Cij(25,10)-6*Cij(34,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(27,11)=Rr(1)
cFC       Cij(28,11)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(3,10)-B23(3,8)-r10*Cij(26,10)-4*Cij(36,11)
cFC       Pr(2)=-B13(3,10)-r21*Cij(26,10)-8*Cij(35,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(28,11)=Rr(1)
cFC       Cij(29,11)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(3,10)+B23(3,9)-r10*Cij(27,10)-2*Cij(32,11)
cFC       Pr(2)=-B13(3,10)-r21*Cij(27,10)-10*Cij(36,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(29,11)=Rr(1)
cFC       Cij(30,11)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(3,10)-B23(3,10)-r10*Cij(22,10)
cFC       Pr(2)=-B13(3,10)-r21*Cij(22,10)-12*Cij(32,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(30,11)=Rr(1)
cFC       Cij(24,11)=Rr(2)
cFCc                {0, 0, ii1, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,10)-B23(2,2)-r10*Cij(12,10)-16*Cij(23,11)
cFC       Pr(2)=B12(2,10)-B13(2,10)-r21*Cij(12,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFC       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFC       Cij(13,11)=Rr(1)
cFC       Cij(15,11)=Rr(2)
cFCc                {0, 0, ii1, 2, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,10)+B23(2,3)-r10*Cij(14,10)-14*Cij(25,11)
cFC       Pr(2)=-B13(2,10)-r21*Cij(14,10)-2*Cij(23,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(15,11)=Rr(1)
cFC       Cij(16,11)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,10)-B23(2,4)-r10*Cij(15,10)-12*Cij(26,11)
cFC       Pr(2)=-B13(2,10)-r21*Cij(15,10)-4*Cij(25,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(16,11)=Rr(1)
cFC       Cij(17,11)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,10)+B23(2,5)-r10*Cij(16,10)-10*Cij(27,11)
cFC       Pr(2)=-B13(2,10)-r21*Cij(16,10)-6*Cij(26,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(17,11)=Rr(1)
cFC       Cij(18,11)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,10)-B23(2,6)-r10*Cij(17,10)-8*Cij(28,11)
cFC       Pr(2)=-B13(2,10)-r21*Cij(17,10)-8*Cij(27,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(18,11)=Rr(1)
cFC       Cij(19,11)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(2,10)+B23(2,7)-r10*Cij(18,10)-6*Cij(29,11)
cFC       Pr(2)=-B13(2,10)-r21*Cij(18,10)-10*Cij(28,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(19,11)=Rr(1)
cFC       Cij(20,11)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(2,10)-B23(2,8)-r10*Cij(19,10)-4*Cij(30,11)
cFC       Pr(2)=-B13(2,10)-r21*Cij(19,10)-12*Cij(29,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(20,11)=Rr(1)
cFC       Cij(21,11)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(2,10)+B23(2,9)-r10*Cij(20,10)-2*Cij(24,11)
cFC       Pr(2)=-B13(2,10)-r21*Cij(20,10)-14*Cij(30,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(21,11)=Rr(1)
cFC       Cij(22,11)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(2,10)-B23(2,10)-r10*Cij(13,10)
cFC       Pr(2)=-B13(2,10)-r21*Cij(13,10)-16*Cij(24,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(22,11)=Rr(1)
cFC       Cij(14,11)=Rr(2)
cFCc                {ii1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=-B023+B13(1,10)-r10*Cij(1,10)-20*Cij(13,11)
cFC       Pr(2)=B12(1,10)-B13(1,10)-r21*Cij(1,10)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFC       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFC       Cij(1,11)=Rr(1)
cFC       Cij(3,11)=Rr(2)
cFCc                {ii1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,10)+B23(1,1)-r10*Cij(3,10)-18*Cij(15,11)
cFC       Pr(2)=-B13(1,10)-r21*Cij(3,10)-2*Cij(13,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(3,11)=Rr(1)
cFC       Cij(4,11)=Rr(2)
cFCc                {ii1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,10)-B23(1,2)-r10*Cij(4,10)-16*Cij(16,11)
cFC       Pr(2)=-B13(1,10)-r21*Cij(4,10)-4*Cij(15,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(4,11)=Rr(1)
cFC       Cij(5,11)=Rr(2)
cFCc                {ii1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,10)+B23(1,3)-r10*Cij(5,10)-14*Cij(17,11)
cFC       Pr(2)=-B13(1,10)-r21*Cij(5,10)-6*Cij(16,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(5,11)=Rr(1)
cFC       Cij(6,11)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,10)-B23(1,4)-r10*Cij(6,10)-12*Cij(18,11)
cFC       Pr(2)=-B13(1,10)-r21*Cij(6,10)-8*Cij(17,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(6,11)=Rr(1)
cFC       Cij(7,11)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,10)+B23(1,5)-r10*Cij(7,10)-10*Cij(19,11)
cFC       Pr(2)=-B13(1,10)-r21*Cij(7,10)-10*Cij(18,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(7,11)=Rr(1)
cFC       Cij(8,11)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,10)-B23(1,6)-r10*Cij(8,10)-8*Cij(20,11)
cFC       Pr(2)=-B13(1,10)-r21*Cij(8,10)-12*Cij(19,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(8,11)=Rr(1)
cFC       Cij(9,11)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(1,10)+B23(1,7)-r10*Cij(9,10)-6*Cij(21,11)
cFC       Pr(2)=-B13(1,10)-r21*Cij(9,10)-14*Cij(20,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(9,11)=Rr(1)
cFC       Cij(10,11)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(1,10)-B23(1,8)-r10*Cij(10,10)-4*Cij(22,11)
cFC       Pr(2)=-B13(1,10)-r21*Cij(10,10)-16*Cij(21,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(10,11)=Rr(1)
cFC       Cij(11,11)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(1,10)+B23(1,9)-r10*Cij(11,10)-2*Cij(14,11)
cFC       Pr(2)=-B13(1,10)-r21*Cij(11,10)-18*Cij(22,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(11,11)=Rr(1)
cFC       Cij(12,11)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(1,10)-B23(1,10)-r10*Cij(2,10)
cFC       Pr(2)=-B13(1,10)-r21*Cij(2,10)-20*Cij(14,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(12,11)=Rr(1)
cFC       Cij(2,11)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
cFC       Cij(49,12)=(-5*p1sq**5-5*p1sq**4*(p2sq+s12)-p1sq**3*(5*p2sq**2+8*
cFC     &  p2sq*s12+5*s12**2)-p1sq**2*(5*p2sq**3+9*p2sq**2*s12+9*p2sq*s12**
cFC     &  2+5*s12**3)-p1sq*(5*p2sq**4+8*p2sq**3*s12+9*p2sq**2*s12**2+8*p2s
cFC     &  q*s12**3+5*s12**4)-5*(p2sq**5+p2sq**4*s12+p2sq**3*s12**2+p2sq**2
cFC     &  *s12**3+p2sq*s12**4+s12**5))/7.6640256d9+(B23(6,10)+r10*Cij(41
cFC     &  ,11)+(r21)*Cij(42,11))/24.d0
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ii1, 1}
cFC       Pr(1)=B13(6,11)+B23(6,10)-r10*Cij(41,11)-2*Cij(49,12)
cFC       Pr(2)=B12(6,11)-B13(6,11)-r21*Cij(41,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFC       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFC       Cij(46,12)=Rr(1)
cFC       Cij(48,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ii1, 2}
cFC       Pr(1)=B13(6,11)-B23(6,11)-r10*Cij(42,11)
cFC       Pr(2)=-B13(6,11)-r21*Cij(42,11)-2*Cij(49,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(48,12)=Rr(1)
cFC       Cij(47,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, ii1, 1, 1, 1}
cFC       Pr(1)=B13(5,11)+B23(5,8)-r10*Cij(37,11)-6*Cij(46,12)
cFC       Pr(2)=B12(5,11)-B13(5,11)-r21*Cij(37,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFC       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFC       Cij(41,12)=Rr(1)
cFC       Cij(43,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, ii1, 2, 1, 1}
cFC       Pr(1)=B13(5,11)-B23(5,9)-r10*Cij(39,11)-4*Cij(48,12)
cFC       Pr(2)=-B13(5,11)-r21*Cij(39,11)-2*Cij(46,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(43,12)=Rr(1)
cFC       Cij(44,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, ii1, 2, 2, 1}
cFC       Pr(1)=B13(5,11)+B23(5,10)-r10*Cij(40,11)-2*Cij(47,12)
cFC       Pr(2)=-B13(5,11)-r21*Cij(40,11)-4*Cij(48,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(44,12)=Rr(1)
cFC       Cij(45,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, 0, 0, ii1, 2, 2, 2}
cFC       Pr(1)=B13(5,11)-B23(5,11)-r10*Cij(38,11)
cFC       Pr(2)=-B13(5,11)-r21*Cij(38,11)-6*Cij(47,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(45,12)=Rr(1)
cFC       Cij(42,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(4,11)+B23(4,6)-r10*Cij(31,11)-10*Cij(41,12)
cFC       Pr(2)=B12(4,11)-B13(4,11)-r21*Cij(31,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFC       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFC       Cij(34,12)=Rr(1)
cFC       Cij(36,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(4,11)-B23(4,7)-r10*Cij(33,11)-8*Cij(43,12)
cFC       Pr(2)=-B13(4,11)-r21*Cij(33,11)-2*Cij(41,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(36,12)=Rr(1)
cFC       Cij(37,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(4,11)+B23(4,8)-r10*Cij(34,11)-6*Cij(44,12)
cFC       Pr(2)=-B13(4,11)-r21*Cij(34,11)-4*Cij(43,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(37,12)=Rr(1)
cFC       Cij(38,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(4,11)-B23(4,9)-r10*Cij(35,11)-4*Cij(45,12)
cFC       Pr(2)=-B13(4,11)-r21*Cij(35,11)-6*Cij(44,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(38,12)=Rr(1)
cFC       Cij(39,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(4,11)+B23(4,10)-r10*Cij(36,11)-2*Cij(42,12)
cFC       Pr(2)=-B13(4,11)-r21*Cij(36,11)-8*Cij(45,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(39,12)=Rr(1)
cFC       Cij(40,12)=Rr(2)
cFCc                {0, 0, 0, 0, 0, 0, ii1, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(4,11)-B23(4,11)-r10*Cij(32,11)
cFC       Pr(2)=-B13(4,11)-r21*Cij(32,11)-10*Cij(42,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(40,12)=Rr(1)
cFC       Cij(35,12)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(3,11)+B23(3,4)-r10*Cij(23,11)-14*Cij(34,12)
cFC       Pr(2)=B12(3,11)-B13(3,11)-r21*Cij(23,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFC       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFC       Cij(25,12)=Rr(1)
cFC       Cij(27,12)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(3,11)-B23(3,5)-r10*Cij(25,11)-12*Cij(36,12)
cFC       Pr(2)=-B13(3,11)-r21*Cij(25,11)-2*Cij(34,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(27,12)=Rr(1)
cFC       Cij(28,12)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(3,11)+B23(3,6)-r10*Cij(26,11)-10*Cij(37,12)
cFC       Pr(2)=-B13(3,11)-r21*Cij(26,11)-4*Cij(36,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(28,12)=Rr(1)
cFC       Cij(29,12)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(3,11)-B23(3,7)-r10*Cij(27,11)-8*Cij(38,12)
cFC       Pr(2)=-B13(3,11)-r21*Cij(27,11)-6*Cij(37,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(29,12)=Rr(1)
cFC       Cij(30,12)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(3,11)+B23(3,8)-r10*Cij(28,11)-6*Cij(39,12)
cFC       Pr(2)=-B13(3,11)-r21*Cij(28,11)-8*Cij(38,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(30,12)=Rr(1)
cFC       Cij(31,12)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(3,11)-B23(3,9)-r10*Cij(29,11)-4*Cij(40,12)
cFC       Pr(2)=-B13(3,11)-r21*Cij(29,11)-10*Cij(39,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(31,12)=Rr(1)
cFC       Cij(32,12)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(3,11)+B23(3,10)-r10*Cij(30,11)-2*Cij(35,12)
cFC       Pr(2)=-B13(3,11)-r21*Cij(30,11)-12*Cij(40,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(32,12)=Rr(1)
cFC       Cij(33,12)=Rr(2)
cFCc                {0, 0, 0, 0, ii1, 2, 2, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(3,11)-B23(3,11)-r10*Cij(24,11)
cFC       Pr(2)=-B13(3,11)-r21*Cij(24,11)-14*Cij(35,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(33,12)=Rr(1)
cFC       Cij(26,12)=Rr(2)
cFCc                {0, 0, ii1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,11)+B23(2,2)-r10*Cij(13,11)-18*Cij(25,12)
cFC       Pr(2)=B12(2,11)-B13(2,11)-r21*Cij(13,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFC       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFC       Cij(14,12)=Rr(1)
cFC       Cij(16,12)=Rr(2)
cFCc                {0, 0, ii1, 2, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,11)-B23(2,3)-r10*Cij(15,11)-16*Cij(27,12)
cFC       Pr(2)=-B13(2,11)-r21*Cij(15,11)-2*Cij(25,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(16,12)=Rr(1)
cFC       Cij(17,12)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,11)+B23(2,4)-r10*Cij(16,11)-14*Cij(28,12)
cFC       Pr(2)=-B13(2,11)-r21*Cij(16,11)-4*Cij(27,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(17,12)=Rr(1)
cFC       Cij(18,12)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,11)-B23(2,5)-r10*Cij(17,11)-12*Cij(29,12)
cFC       Pr(2)=-B13(2,11)-r21*Cij(17,11)-6*Cij(28,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(18,12)=Rr(1)
cFC       Cij(19,12)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,11)+B23(2,6)-r10*Cij(18,11)-10*Cij(30,12)
cFC       Pr(2)=-B13(2,11)-r21*Cij(18,11)-8*Cij(29,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(19,12)=Rr(1)
cFC       Cij(20,12)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(2,11)-B23(2,7)-r10*Cij(19,11)-8*Cij(31,12)
cFC       Pr(2)=-B13(2,11)-r21*Cij(19,11)-10*Cij(30,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(20,12)=Rr(1)
cFC       Cij(21,12)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(2,11)+B23(2,8)-r10*Cij(20,11)-6*Cij(32,12)
cFC       Pr(2)=-B13(2,11)-r21*Cij(20,11)-12*Cij(31,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(21,12)=Rr(1)
cFC       Cij(22,12)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(2,11)-B23(2,9)-r10*Cij(21,11)-4*Cij(33,12)
cFC       Pr(2)=-B13(2,11)-r21*Cij(21,11)-14*Cij(32,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(22,12)=Rr(1)
cFC       Cij(23,12)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(2,11)+B23(2,10)-r10*Cij(22,11)-2*Cij(26,12)
cFC       Pr(2)=-B13(2,11)-r21*Cij(22,11)-16*Cij(33,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(23,12)=Rr(1)
cFC       Cij(24,12)=Rr(2)
cFCc                {0, 0, ii1, 2, 2, 2, 2, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(2,11)-B23(2,11)-r10*Cij(14,11)
cFC       Pr(2)=-B13(2,11)-r21*Cij(14,11)-18*Cij(26,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(24,12)=Rr(1)
cFC       Cij(15,12)=Rr(2)
cFCc                {ii1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B023+B13(1,11)-r10*Cij(1,11)-22*Cij(14,12)
cFC       Pr(2)=B12(1,11)-B13(1,11)-r21*Cij(1,11)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFC       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFC       Cij(1,12)=Rr(1)
cFC       Cij(3,12)=Rr(2)
cFCc                {ii1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,11)-B23(1,1)-r10*Cij(3,11)-20*Cij(16,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(3,11)-2*Cij(14,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(3,12)=Rr(1)
cFC       Cij(4,12)=Rr(2)
cFCc                {ii1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,11)+B23(1,2)-r10*Cij(4,11)-18*Cij(17,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(4,11)-4*Cij(16,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(4,12)=Rr(1)
cFC       Cij(5,12)=Rr(2)
cFCc                {ii1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,11)-B23(1,3)-r10*Cij(5,11)-16*Cij(18,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(5,11)-6*Cij(17,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(5,12)=Rr(1)
cFC       Cij(6,12)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,11)+B23(1,4)-r10*Cij(6,11)-14*Cij(19,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(6,11)-8*Cij(18,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(6,12)=Rr(1)
cFC       Cij(7,12)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,11)-B23(1,5)-r10*Cij(7,11)-12*Cij(20,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(7,11)-10*Cij(19,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(7,12)=Rr(1)
cFC       Cij(8,12)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,11)+B23(1,6)-r10*Cij(8,11)-10*Cij(21,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(8,11)-12*Cij(20,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(8,12)=Rr(1)
cFC       Cij(9,12)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1}
cFC       Pr(1)=B13(1,11)-B23(1,7)-r10*Cij(9,11)-8*Cij(22,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(9,11)-14*Cij(21,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(9,12)=Rr(1)
cFC       Cij(10,12)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1}
cFC       Pr(1)=B13(1,11)+B23(1,8)-r10*Cij(10,11)-6*Cij(23,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(10,11)-16*Cij(22,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(10,12)=Rr(1)
cFC       Cij(11,12)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1}
cFC       Pr(1)=B13(1,11)-B23(1,9)-r10*Cij(11,11)-4*Cij(24,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(11,11)-18*Cij(23,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(11,12)=Rr(1)
cFC       Cij(12,12)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1}
cFC       Pr(1)=B13(1,11)+B23(1,10)-r10*Cij(12,11)-2*Cij(15,12)
cFC       Pr(2)=-B13(1,11)-r21*Cij(12,11)-20*Cij(24,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(12,12)=Rr(1)
cFC       Cij(13,12)=Rr(2)
cFCc                {ii1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}
cFC       Pr(1)=B13(1,11)-B23(1,11)-r10*Cij(2,11)
cFC       Pr(2)=-B13(1,11)-r21*Cij(2,11)-22*Cij(15,12)
cFC       Rr(2)=(Pr(1)-z21*Pr(2))*iz22
cFCc       Rr(1)=(Pr(2)-  z12*Rr(2))*iz11
cFCc       Cij(13,12)=Rr(1)
cFC       Cij(2,12)=Rr(2)


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
