CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                tens_red3 = 3-point tensors
C     Francisco Campanario
C    Email: Francam@particle.uni-karlsruhe.de
C Date: 4/02/2008
C  Tensor reduction a la Passarino-Vetlman for the tensor
C In this file, we should put as an input the 3 scalar invariants, 
C the 3 B0, and C0.
C Notation: Example CR. the "R" means the real part of C0
C                   CI, the "I" means the imaginary part of C0
C It speeds up the code. 
C OUTPUT:  is given in C0r,C0i,Cijr,CijI. "r":  real part, "I" imaginary part
C NOTE!!!! The code seems to be 3 times faster than the old one..
C The indices are pertumted in comparison with the old subroutines.
C the rank of the tensor integral is store in the second argument.
C C_old(1,2)-> C(2,1)(rank 1 coefficient 2) it is faster in this way. 
C I was not able to find this recursion relations in the literature.
C I did it myself, so, one has to look into my notes if want to generalize them.
C I have as well the case of a massive loop with one mass with no restrictions in the invariants.
C I have more version which only needs the input in complex form
C and the output is given as a complex array. The subroutine still makes the separation
C between imaginary part and real part, but it is done inside of it. Here I am saving doing that since 
C I need and use this CR and CI in my penline and boxline subroutines.
C If you need other case ask me!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine tens_red3_new_Re_Com_Div(p1sq,p2sq,s12,B0_23,B0_13,B0_12, 
     #                     C0,!Cij,
     # C0r,C0I,Cijr,CijI)
      implicit none
      real * 8  p1sq, p2sq, s12
      complex*16 B0_23, B0_13, B0_12, C0!, Cij(4,2)
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
      real*8 r1, r2r1, deter,Ideter,p1p2,IZ(2,2)
      real*8 B1r_12,B1r_13,B1r_23,Rr(2),Cijr(4,2)
      real*8 B1I_12,B1I_13,B1I_23,RI(2),CijI(4,2)
      real*8 B0r_23, B0r_13, B0r_12, C0r
      real*8 B0I_23, B0I_13, B0I_12, C0I 

      p1p2 = (s12 - p1sq - p2sq)*0.5d0

      r1 = p1sq
      r2r1 = s12 - r1
      

      deter = 2.d0*(p1sq*p2sq - p1p2*p1p2)
      Ideter=1.d0/deter
     
      IZ(1,1) =  p2sq*Ideter
      IZ(2,1) = -p1p2*Ideter
      IZ(1,2) =  IZ(2,1)
      IZ(2,2) =  p1sq*Ideter

     
      B0r_12=Dble(B0_12)
      B0r_13=Dble(B0_13)
      B0r_23=Dble(B0_23)
      C0r=Dble(C0)

      B0I_12=DImag(B0_12)
      B0I_13=DImag(B0_13)
      B0I_23=DImag(B0_23)
      C0I=DImag(C0)

      B1r_12 = -B0r_12*0.5d0
      B1r_13 = -B0r_13*0.5d0
      B1r_23 = -B0r_23*0.5d0
      
      B1I_12 = -B0I_12*0.5d0
      B1I_13 = -B0I_13*0.5d0
      B1I_23 = -B0I_23*0.5d0


c 1-2
      Rr(1) = (B0r_13 - B0r_23 - C0r*r1)
      Rr(2) = (B0r_12 - B0r_13 - C0r*r2r1)
      
      RI(1) = (B0I_13 - B0I_23 - C0I*r1)
      RI(2) = (B0I_12 - B0I_13 - C0I*r2r1)
      

      Cijr(1,1) = (Rr(1)*IZ(1,1) + Rr(2)*IZ(2,1))
      Cijr(2,1) = (Rr(1)*IZ(1,2) + Rr(2)*IZ(2,2))
      
      CijI(1,1) = (RI(1)*IZ(1,1) + RI(2)*IZ(2,1))
      CijI(2,1) = (RI(1)*IZ(1,2) + RI(2)*IZ(2,2))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c C00
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Cijr(4,2) = ( B0r_23 + Cijr(1,1)*r1 +Cijr(2,1)*r2r1 
cccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc
c Commented out for divergent part
c      +1.d0
ccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccc
c
     -  )*0.25d0

      CijI(4,2) = ( B0I_23 + CijI(1,1)*r1 +CijI(2,1)*r2r1)*0.25d0
      
c 3-5
      Rr(1) = (B1r_13 + B0r_23 - Cijr(1,1)*r1 - Cijr(4,2)*2.d0)
      Rr(2) = (B1r_12 - B1r_13 - Cijr(1,1)*r2r1)
     
      RI(1) = (B1I_13 + B0I_23 - CijI(1,1)*r1 - CijI(4,2)*2.d0)
      RI(2) = (B1I_12 - B1I_13 - CijI(1,1)*r2r1)
      
      Cijr(1,2) =(Rr(1)*IZ(1,1)+Rr(2)*IZ(2,1))
      CijI(1,2) =(RI(1)*IZ(1,1)+RI(2)*IZ(2,1))


c 4-6
      Rr(1) = (  B1r_13 - B1r_23 - Cijr(2,1)*r1)
      Rr(2) = (- B1r_13          -Cijr(2,1)*r2r1 -Cijr(4,2)*2.d0)
      
      RI(1) = (  B1I_13 - B1I_23 - CijI(2,1)*r1)
      RI(2) = (- B1I_13          - CijI(2,1)*r2r1 -CijI(4,2)*2.d0)

      Cijr(2,2) = (Rr(1)*IZ(1,2)+Rr(2)*IZ(2,2))
      Cijr(3,2) = (Rr(1)*IZ(1,1)+Rr(2)*IZ(2,1))

      CijI(2,2) = (RI(1)*IZ(1,2)+RI(2)*IZ(2,2))
      CijI(3,2) = (RI(1)*IZ(1,1)+RI(2)*IZ(2,1))
     
c$$$      Cij(1,1)=DCMPLX(Cijr(1,1),CijI(1,1))
c$$$      Cij(2,1)=DCMPLX(Cijr(2,1),CijI(2,1))
c$$$
c$$$      Cij(1,2)=DCMPLX(Cijr(1,2),CijI(1,2))
c$$$      Cij(2,2)=DCMPLX(Cijr(2,2),CijI(2,2))
c$$$      Cij(3,2)=DCMPLX(Cijr(3,2),CijI(3,2))
c$$$      Cij(4,2)=DCMPLX(Cijr(4,2),CijI(4,2))


      end



