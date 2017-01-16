CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                tens_red4 = 4-point tensors
C     Francisco Campanario
C    Email: Francam@particle.uni-karlsruhe.de
C Date: 4/02/2008
C  Tensor reduction a la Passarino-Vetlman for the tensor
C In this file, we should put as an input the 6 scalar invariants,
C the 4 C0, and the Cij got it from tens_red3_new_Re_Com and D0.
C Notation: Example CR234. the "R" means the real part of C0234
C                   CI234, the "I" means the imaginary part of C0234
C This separation is done already in tens_red3_....
C It speeds up the code.
C OUTPUT:  is given in D0r,D0i,Dijr,DijI. "r":  real part, "I" imaginary part
C NOTE!!!! The code seems to be 3 to 4 times faster than the old one..
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
      subroutine tens_red4_new_Re_Com_Div(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,
     &                     C0r_234, C0r_134, C0r_124, C0r_123,
     &                     Cijr_234, Cijr_134, Cijr_124, Cijr_123,
     &                     C0I_234, C0I_134, C0I_124, C0I_123,
     &                     CijI_234, CijI_134, CijI_124, CijI_123,
     &                     D0,!Dij,
     &   D0r, D0I,Dijr,DijI)
      implicit none
      real * 8 p1sq, p2sq, p3sq, p1p2, p1p3,p2p3
c$$$      complex*16 C0_234, C0_134, C0_124, C0_123
c$$$      complex*16 Cij_234(4,2), Cij_134(4,2), Cij_124(4,2), Cij_123(4,2)
      complex*16 D0!, Dij(13,3)
c, Dijp(13,3), Dijpp(13,3)
c      complex*16 t(9)
c
c  determine the Passarino-Veltman tensor decomposition for the four-point
c  tensor integrals
c
c                                          d^4k
c   D0; D_mu; D_mu,nu; D_mu,nu,rho =  Int ------
c                                         (2pi)^4
c
c              1;  k_mu;   k_mu k_nu; k_mu k_nu k_rho
c      -------------------------------------------------------------------
c         [k^2-m^2][(k+p1)^2-m^2][(k+p1+p2)^2-m^2][(k+p1+p2+p3)^2-m^2]
c with
c
c   D_mu = p1_mu D11  +  p2_mu D12  +  p3_mu D13
c
c   D_munu = p1_mu p2_nu D21 + p2_mu p2_nu D22 + ...
c
c  for notation see Passarino&Veltman, NP B160 (1979) 151
c
C INPUT:  psq, pq,...                        kinematics invariants
C         C0_123 = C0(1,2,3) = C0(p1,p2)     scalar three point
C         C0_124 = C0(1,2,4) = C0(p1,p2+p3)  functions in PV notation
C         C0_134 = C0(1,3,4) = C0(p1+p2,p3)
C         C0_234 = C0(2,3,4) = C0(p2,p3)
C         Cij_123(n,m) = C_nm(1,2,3) ....    higher C_nm form factors
C                                            as in tens_red3
c         D0 = D0(p,q,l)                  scalar four point function
c
c OUTPUT: Dij(n,m) = D_nm                    form factors in the tensor
c                                            integrals a la PV
c         nm = 11, 21, 31                    ff's for D_mu
c         nm = 12, 22, 32, 42, 52, 62, 72    ff's for D_munu
c         nm = 13, 23, 33, ..., 93, 103, 113, 123  ff's for D_mu,nu,rho
c
      real*8 r1, r2r1, r3r2, IZ(3,3), deter,Ideter,Iv
      real*8 Rr(1:3), RI(1:3)
      real*8 C0r_234, C0r_134, C0r_124, C0r_123
      real*8 Cijr_234(4,2), Cijr_134(4,2), Cijr_124(4,2), Cijr_123(4,2)
      real*8 C0I_234, C0I_134, C0I_124, C0I_123
      real*8 CijI_234(4,2), CijI_134(4,2), CijI_124(4,2), CijI_123(4,2)
      real*8 D0r, Dijr(13,3)
      real*8 D0I, DijI(13,3)


      r1   = p1sq
      r2r1 = p2sq+2.d0*p1p2
      r3r2 = 2.d0*(p1p3+p2p3)+p3sq

      deter =-2.d0*(-2.d0*p1p2*p1p3*p2p3 + p1p3*p1p3*p2sq + p1p2*p1p2*p3sq
     -  + p1sq*(p2p3*p2p3 - p2sq*p3sq))

      Ideter=1.d0/deter

      IZ(1,1) = (p2sq*p3sq - p2p3*p2p3)*Ideter
      IZ(2,1) = (p1p3*p2p3 - p1p2*p3sq)*Ideter
      IZ(3,1) = (p1p2*p2p3 - p1p3*p2sq)*Ideter
      IZ(1,2) = IZ(2,1)
      IZ(2,2) = (p1sq*p3sq-p1p3*p1p3)*Ideter
      IZ(3,2) = (p1p2*p1p3 - p1sq*p2p3)*Ideter
      IZ(1,3) = IZ(3,1)
      IZ(2,3) = IZ(3,2)
      IZ(3,3) = (p1sq*p2sq - p1p2*p1p2)*Ideter

      D0r=Dble(D0)
      D0I=DImag(D0)
c
      Rr(1) = (C0r_134 - C0r_234 - D0r*r1   )
      Rr(2) = (C0r_124 - C0r_134 - D0r*r2r1 )
      Rr(3) = (C0r_123 - C0r_124 - D0r*r3r2 )
c
      Dijr(1,1) =( Rr(1)*IZ(1,1) + Rr(2)*IZ(2,1) + Rr(3)*IZ(3,1)  )
      Dijr(2,1) =( Rr(1)*IZ(1,2) + Rr(2)*IZ(2,2) + Rr(3)*IZ(3,2)  )
      Dijr(3,1) =( Rr(1)*IZ(1,3) + Rr(2)*IZ(2,3) + Rr(3)*IZ(3,3)  )

c
      RI(1) = (C0I_134 - C0I_234 - D0I*r1   )
      RI(2) = (C0I_124 - C0I_134 - D0I*r2r1 )
      RI(3) = (C0I_123 - C0I_124 - D0I*r3r2 )
c
      DijI(1,1) =( RI(1)*IZ(1,1) + RI(2)*IZ(2,1) + RI(3)*IZ(3,1)  )
      DijI(2,1) =( RI(1)*IZ(1,2) + RI(2)*IZ(2,2) + RI(3)*IZ(3,2)  )
      DijI(3,1) =( RI(1)*IZ(1,3) + RI(2)*IZ(2,3) + RI(3)*IZ(3,3)  )


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c D00
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Dijr(7,2) = (C0r_234 + Dijr(1,1)*r1 + Dijr(2,1)*r2r1
     &     +Dijr(3,1)*r3r2)*0.5d0

      DijI(7,2) = (C0I_234 + DijI(1,1)*r1 + DijI(2,1)*r2r1
     &     +DijI(3,1)*r3r2)*0.5d0


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c 30-31-32

c
      Rr(1) = (Cijr_134(1,1) + C0r_234 - Dijr(1,1)*r1 - Dijr(7,2)*2.d0)
      Rr(2) = (Cijr_124(1,1) - Cijr_134(1,1) - Dijr(1,1)*r2r1        )
      Rr(3) = (Cijr_123(1,1) - Cijr_124(1,1) - Dijr(1,1)*r3r2        )
c
      Dijr(1,2) = ( Rr(1)*IZ(1,1) + Rr(2)*IZ(2,1) + Rr(3)*IZ(3,1)  )
      Dijr(4,2) = ( Rr(1)*IZ(1,2) + Rr(2)*IZ(2,2) + Rr(3)*IZ(3,2)  )
      Dijr(5,2) = ( Rr(1)*IZ(1,3) + Rr(2)*IZ(2,3) + Rr(3)*IZ(3,3)  )
c
      RI(1) = (CijI_134(1,1) + C0I_234  - DijI(1,1)*r1 - DijI(7,2)*2.d0)
      RI(2) = (CijI_124(1,1) - CijI_134(1,1) - DijI(1,1)*r2r1        )
      RI(3) = (CijI_123(1,1) - CijI_124(1,1) - DijI(1,1)*r3r2        )
c
      DijI(1,2) = ( RI(1)*IZ(1,1) + RI(2)*IZ(2,1) + RI(3)*IZ(3,1)  )
      DijI(4,2) = ( RI(1)*IZ(1,2) + RI(2)*IZ(2,2) + RI(3)*IZ(3,2)  )
      DijI(5,2) = ( RI(1)*IZ(1,3) + RI(2)*IZ(2,3) + RI(3)*IZ(3,3)  )

c 33-34-35

      Rr(1) = (Cijr_134(1,1) -Cijr_234(1,1) - Dijr(2,1)*r1          )
      Rr(2)=(Cijr_124(2,1)-Cijr_134(1,1)-Dijr(2,1)*r2r1-Dijr(7,2)*2.d0)
      Rr(3) = (Cijr_123(2,1) -Cijr_124(2,1) - Dijr(2,1)*r3r2        )

c      Dijr(4,2) = ( Rr(1)*IZ(1,1) + Rr(2)*IZ(2,1) + Rr(3)*IZ(3,1)  )
      Dijr(2,2) = ( Rr(1)*IZ(1,2) + Rr(2)*IZ(2,2) + Rr(3)*IZ(3,2)  )
      Dijr(6,2) = ( Rr(1)*IZ(1,3) + Rr(2)*IZ(2,3) + Rr(3)*IZ(3,3)  )

      RI(1) = (CijI_134(1,1) -CijI_234(1,1) - DijI(2,1)*r1          )
      RI(2)=(CijI_124(2,1)-CijI_134(1,1)-DijI(2,1)*r2r1-DijI(7,2)*2.d0)
      RI(3) = (CijI_123(2,1) -CijI_124(2,1) - DijI(2,1)*r3r2        )

c      DijI(4,2) = ( RI(1)*IZ(1,1) + RI(2)*IZ(2,1) + RI(3)*IZ(3,1)  )
      DijI(2,2) = ( RI(1)*IZ(1,2) + RI(2)*IZ(2,2) + RI(3)*IZ(3,2)  )
      DijI(6,2) = ( RI(1)*IZ(1,3) + RI(2)*IZ(2,3) + RI(3)*IZ(3,3)  )

c 36-37-38

      Rr(1) = (Cijr_134(2,1) - Cijr_234(2,1) - Dijr(3,1)*r1      )
      Rr(2) = (Cijr_124(2,1) - Cijr_134(2,1) - Dijr(3,1)*r2r1    )
      Rr(3) = (       - Cijr_124(2,1) - Dijr(3,1)*r3r2 -Dijr(7,2)*2.d0)

c      Dijr(5,2) = ( Rr(1)*IZ(1,1) + Rr(2)*IZ(2,1) + Rr(3)*IZ(3,1)  )
c      Dijr(6,2) = ( Rr(1)*IZ(1,2) + Rr(2)*IZ(2,2) + Rr(3)*IZ(3,2)  )
      Dijr(3,2) = ( Rr(1)*IZ(1,3) + Rr(2)*IZ(2,3) + Rr(3)*IZ(3,3)  )

      RI(1) = (CijI_134(2,1) - CijI_234(2,1) - DijI(3,1)*r1      )
      RI(2) = (CijI_124(2,1) - CijI_134(2,1) - DijI(3,1)*r2r1    )
      RI(3) = (       - CijI_124(2,1) - DijI(3,1)*r3r2 -DijI(7,2)*2.d0)

c      DijI(5,2) = ( RI(1)*IZ(1,1) + RI(2)*IZ(2,1) + RI(3)*IZ(3,1)  )
c      DijI(6,2) = ( RI(1)*IZ(1,2) + RI(2)*IZ(2,2) + RI(3)*IZ(3,2)  )
      DijI(3,2) = ( RI(1)*IZ(1,3) + RI(2)*IZ(2,3) + RI(3)*IZ(3,3)  )


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      Dijr(11,3)=0.25d0*(-C0r_234+Dijr(1,2)*r1+Dijr(4,2)*r2r1
     -      +Dijr(5,2)*r3r2)

      Dijr(12,3)=0.25d0*(Cijr_234(1,1)+Dijr(4,2)*r1+Dijr(2,2)*r2r1
     -                  +Dijr(6,2)*r3r2)
      Dijr(13,3)=0.25d0*(Cijr_234(2,1)+Dijr(5,2)*r1 + Dijr(6,2)*r2r1
     -                   + Dijr(3,2)*r3r2)


      DijI(11,3)=0.25d0*(-C0I_234+DijI(1,2)*r1+DijI(4,2)*r2r1
     -  +DijI(5,2)*r3r2)

      DijI(12,3)=0.25d0*(CijI_234(1,1)+DijI(4,2)*r1+DijI(2,2)*r2r1
     -                  +DijI(6,2)*r3r2)
      DijI(13,3)=0.25d0*(CijI_234(2,1)+DijI(5,2)*r1 + DijI(6,2)*r2r1
     -                   + DijI(3,2)*r3r2)


c 41-42-43


      Rr(1) =(Cijr_134(1,2) - C0r_234 - Dijr(1,2)*r1  - Dijr(11,3)*4.d0)
      Rr(2) =(Cijr_124(1,2) - Cijr_134(1,2) - Dijr(1,2)*r2r1       )
      Rr(3) =(Cijr_123(1,2) - Cijr_124(1,2) - Dijr(1,2)*r3r2       )

      Dijr(1,3) = ( Rr(1)*IZ(1,1) + Rr(2)*IZ(2,1) + Rr(3)*IZ(3,1)  )
      Dijr(4,3) = ( Rr(1)*IZ(1,2) + Rr(2)*IZ(2,2) + Rr(3)*IZ(3,2)  )
      Dijr(5,3) = ( Rr(1)*IZ(1,3) + Rr(2)*IZ(2,3) + Rr(3)*IZ(3,3)  )

      RI(1) =(CijI_134(1,2) - C0I_234 - DijI(1,2)*r1  - DijI(11,3)*4.d0)
      RI(2) =(CijI_124(1,2) - CijI_134(1,2) - DijI(1,2)*r2r1       )
      RI(3) =(CijI_123(1,2) - CijI_124(1,2) - DijI(1,2)*r3r2       )

      DijI(1,3) = ( RI(1)*IZ(1,1) + RI(2)*IZ(2,1) + RI(3)*IZ(3,1)  )
      DijI(4,3) = ( RI(1)*IZ(1,2) + RI(2)*IZ(2,2) + RI(3)*IZ(3,2)  )
      DijI(5,3) = ( RI(1)*IZ(1,3) + RI(2)*IZ(2,3) + RI(3)*IZ(3,3)  )

c 50-51-52

      Rr(1) =  (Cijr_134(1,2) - Cijr_234(1,2) - Dijr(2,2)*r1)
      Rr(2) =  (Cijr_124(2,2) - Cijr_134(1,2) - Dijr(2,2)*r2r1
     &     -Dijr(12,3)*4.d0)
      Rr(3) =  (Cijr_123(2,2) - Cijr_124(2,2) - Dijr(2,2)*r3r2)


      Dijr(6,3) = ( Rr(1)*IZ(1,1) + Rr(2)*IZ(2,1) + Rr(3)*IZ(3,1)  )
      Dijr(2,3) = ( Rr(1)*IZ(1,2) + Rr(2)*IZ(2,2) + Rr(3)*IZ(3,2)  )
      Dijr(8,3) = ( Rr(1)*IZ(1,3) + Rr(2)*IZ(2,3) + Rr(3)*IZ(3,3)  )

      RI(1) =  (CijI_134(1,2) - CijI_234(1,2) - DijI(2,2)*r1)
      RI(2) =  (CijI_124(2,2) - CijI_134(1,2) - DijI(2,2)*r2r1
     &     -DijI(12,3)*4.d0)
      RI(3) =  (CijI_123(2,2) - CijI_124(2,2) - DijI(2,2)*r3r2)

      DijI(6,3) = ( RI(1)*IZ(1,1) + RI(2)*IZ(2,1) + RI(3)*IZ(3,1)  )
      DijI(2,3) = ( RI(1)*IZ(1,2) + RI(2)*IZ(2,2) + RI(3)*IZ(3,2)  )
      DijI(8,3) = ( RI(1)*IZ(1,3) + RI(2)*IZ(2,3) + RI(3)*IZ(3,3)  )

c
c 56-57-58
c
c
      Rr(1) = (Cijr_134(2,2) - Cijr_234(2,2) - Dijr(3,2)*r1)
      Rr(2) = (Cijr_124(2,2) - Cijr_134(2,2) - Dijr(3,2)*r2r1)
      Rr(3) = (     -Cijr_124(2,2) - Dijr(3,2)*r3r2 - Dijr(13,3)*4.d0)
c
      Dijr(7,3) = ( Rr(1)*IZ(1,1) + Rr(2)*IZ(2,1) + Rr(3)*IZ(3,1)  )
      Dijr(9,3) = ( Rr(1)*IZ(1,2) + Rr(2)*IZ(2,2) + Rr(3)*IZ(3,2)  )
      Dijr(3,3) = ( Rr(1)*IZ(1,3) + Rr(2)*IZ(2,3) + Rr(3)*IZ(3,3)  )
c
      RI(1) = (CijI_134(2,2) - CijI_234(2,2) - DijI(3,2)*r1)
      RI(2) = (CijI_124(2,2) - CijI_134(2,2) - DijI(3,2)*r2r1)
      RI(3) = (     -CijI_124(2,2) - DijI(3,2)*r3r2 - DijI(13,3)*4.d0)
c
      DijI(7,3) = ( RI(1)*IZ(1,1) + RI(2)*IZ(2,1) + RI(3)*IZ(3,1)  )
      DijI(9,3) = ( RI(1)*IZ(1,2) + RI(2)*IZ(2,2) + RI(3)*IZ(3,2)  )
      DijI(3,3) = ( RI(1)*IZ(1,3) + RI(2)*IZ(2,3) + RI(3)*IZ(3,3)  )
c 44-45-46
      Rr(1) = (Cijr_134(1,2) + Cijr_234(1,1) - Dijr(4,2)*r1
     &     -Dijr(12,3)*2.d0)
      Rr(2) = (Cijr_124(3,2) - Cijr_134(1,2) - Dijr(4,2)*r2r1
     &     - Dijr(11,3)*2.d0)
      Rr(3) = (Cijr_123(3,2) - Cijr_124(3,2)  -Dijr(4,2)*r3r2)
c
C      Dijr(4,3) = ( Rr(1)*IZ(1,1) + Rr(2)*IZ(2,1) + Rr(3)*IZ(3,1)  )
C      Dijr(6,3) = ( Rr(1)*IZ(1,2) + Rr(2)*IZ(2,2) + Rr(3)*IZ(3,2)  )
      Dijr(10,3) = ( Rr(1)*IZ(1,3) + Rr(2)*IZ(2,3) + Rr(3)*IZ(3,3)  )
c
c
      RI(1) = (CijI_134(1,2) + CijI_234(1,1) - DijI(4,2)*r1
     &     -DijI(12,3)*2.d0)
      RI(2) = (CijI_124(3,2) - CijI_134(1,2) - DijI(4,2)*r2r1
     &     - DijI(11,3)*2.d0)
      RI(3) = (CijI_123(3,2) - CijI_124(3,2)  -DijI(4,2)*r3r2)

C      DijI(4,3) = ( RI(1)*IZ(1,1) + RI(2)*IZ(2,1) + RI(3)*IZ(3,1)  )
C      DijI(6,3) = ( RI(1)*IZ(1,2) + RI(2)*IZ(2,2) + RI(3)*IZ(3,2)  )
      DijI(10,3) = ( RI(1)*IZ(1,3) + RI(2)*IZ(2,3) + RI(3)*IZ(3,3)  )
c
c to define D00ij and D0000 functions!
c In PV notation we have:
c Dij(1,7)=D416, Dij(1,8)=D417, Dij(1,9)=D418
c Dij(1,10)=D419, Dij(1,11)=D420, Dij(1,12)=D421, D(1,13)=D422
c
      Iv=1.d0/6.d0
c
c
      Dijr(7,1)=(C0r_234+Dijr(1,3)*r1+Dijr(4,3)*r2r1+Dijr(5,3)*r3r2 )*Iv
c
      Dijr(8,1)=(Cijr_234(1,2)+Dijr(6,3)*r1+Dijr(2,3)*r2r1
     - +Dijr(8,3)*r3r2)*Iv
c
      Dijr(9,1)=(Cijr_234(2,2)+Dijr(7,3)*r1+Dijr(9,3)*r2r1
     - +Dijr(3,3)*r3r2)*Iv
c
      Dijr(10,1)=(-Cijr_234(1,1)+Dijr(4,3)*r1+Dijr(6,3)*r2r1
     - +Dijr(10,3)*r3r2)*Iv
c
      Dijr(11,1)=(-Cijr_234(2,1)+Dijr(5,3)*r1+Dijr(10,3)*r2r1
     - +Dijr(7,3)*r3r2)*Iv
c
      Dijr(12,1)=(Cijr_234(3,2)+Dijr(10,3)*r1+Dijr(8,3)*r2r1
     - +Dijr(9,3)*r3r2)*Iv
c
      Dijr(13,1)=(Cijr_234(4,2)+Dijr(11,3)*r1+Dijr(12,3)*r2r1
     -      + Dijr(13,3)*r3r2
ccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccc
c    Comented out for Divergent part
c      +Iv
ccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc
     -   )*Iv
cc
cc
c
      DijI(7,1)=(C0I_234+DijI(1,3)*r1+DijI(4,3)*r2r1+DijI(5,3)*r3r2 )*Iv
c
      DijI(8,1)=(CijI_234(1,2)+DijI(6,3)*r1+DijI(2,3)*r2r1
     - +DijI(8,3)*r3r2)*Iv
c
      DijI(9,1)=(CijI_234(2,2)+DijI(7,3)*r1+DijI(9,3)*r2r1
     -+DijI(3,3)*r3r2)*Iv
c
      DijI(10,1)=(-CijI_234(1,1)+DijI(4,3)*r1+DijI(6,3)*r2r1
     - +DijI(10,3)*r3r2)*Iv
c
      DijI(11,1)=(-CijI_234(2,1)+DijI(5,3)*r1+DijI(10,3)*r2r1
     - +DijI(7,3)*r3r2)*Iv
c
      DijI(12,1)=(CijI_234(3,2)+DijI(10,3)*r1+DijI(8,3)*r2r1
     - +DijI(9,3)*r3r2)*Iv
c
      DijI(13,1)=(CijI_234(4,2)+DijI(11,3)*r1+DijI(12,3)*r2r1
     -      + DijI(13,3)*r3r2)*Iv

c$$$      Dij(1,1)=DCMPLX(Dijr(1,1),DijI(1,1))
c$$$      Dij(2,1)=DCMPLX(Dijr(2,1),DijI(2,1))
c$$$      Dij(3,1)=DCMPLX(Dijr(3,1),DijI(3,1))
c$$$c
c$$$      Dij(7,1)=DCMPLX(Dijr(7,1),DijI(7,1))
c$$$      Dij(8,1)=DCMPLX(Dijr(8,1),DijI(8,1))
c$$$      Dij(9,1)=DCMPLX(Dijr(9,1),DijI(9,1))
c$$$      Dij(10,1)=DCMPLX(Dijr(10,1),DijI(10,1))
c$$$      Dij(11,1)=DCMPLX(Dijr(11,1),DijI(11,1))
c$$$      Dij(12,1)=DCMPLX(Dijr(12,1),DijI(12,1))
c$$$      Dij(13,1)=DCMPLX(Dijr(13,1),DijI(13,1))
c$$$c
c$$$      Dij(1,2)=DCMPLX(Dijr(1,2),DijI(1,2))
c$$$      Dij(2,2)=DCMPLX(Dijr(2,2),DijI(2,2))
c$$$      Dij(3,2)=DCMPLX(Dijr(3,2),DijI(3,2))
c$$$      Dij(4,2)=DCMPLX(Dijr(4,2),DijI(4,2))
c$$$      Dij(5,2)=DCMPLX(Dijr(5,2),DijI(5,2))
c$$$      Dij(6,2)=DCMPLX(Dijr(6,2),DijI(6,2))
c$$$      Dij(7,2)=DCMPLX(Dijr(7,2),DijI(7,2))
c$$$c
c$$$      Dij(1,3)=DCMPLX(Dijr(1,3),DijI(1,3))
c$$$      Dij(2,3)=DCMPLX(Dijr(2,3),DijI(2,3))
c$$$      Dij(3,3)=DCMPLX(Dijr(3,3),DijI(3,3))
c$$$      Dij(4,3)=DCMPLX(Dijr(4,3),DijI(4,3))
c$$$      Dij(5,3)=DCMPLX(Dijr(5,3),DijI(5,3))
c$$$      Dij(6,3)=DCMPLX(Dijr(6,3),DijI(6,3))
c$$$      Dij(7,3)=DCMPLX(Dijr(7,3),DijI(7,3))
c$$$      Dij(8,3)=DCMPLX(Dijr(8,3),DijI(8,3))
c$$$      Dij(9,3)=DCMPLX(Dijr(9,3),DijI(9,3))
c$$$      Dij(10,3)=DCMPLX(Dijr(10,3),DijI(10,3))
c$$$      Dij(11,3)=DCMPLX(Dijr(11,3),DijI(11,3))
c$$$      Dij(12,3)=DCMPLX(Dijr(12,3),DijI(12,3))
c$$$      Dij(13,3)=DCMPLX(Dijr(13,3),DijI(13,3))
c
      end

