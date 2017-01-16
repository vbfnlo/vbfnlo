CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                tens_red4 = 4-point tensors
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine tens_red4(psq, qsq, lsq, pq, pl, ql,
     #                     C0_234, C0_134, C0_124, C0_123,
     #                     Cij_234, Cij_134, Cij_124, Cij_123,
     #                     D0,
     #                     Dij)
      implicit none
      real * 8 psq, qsq, lsq, pq, pl, ql
      complex*16 C0_234, C0_134, C0_124, C0_123
      complex*16 Cij_234(2,4), Cij_134(2,4), Cij_124(2,4), Cij_123(2,4)
      complex*16 D0, Dij(3,13), Dijp(3,13), Dijpp(3,13)
      complex*16 t(9)
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
c         [k^2-m^2][(k+p)^2-m^2][(k+p+q)^2-m^2][(k+p+q+l)^2-m^2]
c with
c
c   D_mu = p_mu D11  +  q_mu D12  +  l_mu D13
c
c   D_munu = p_mu p_nu D21 + q_mu q_nu D22 + ...
c
c  for notation see Passarino&Veltman, NP B160 (1979) 151
c
C INPUT:  psq, pq,...                        kinematics invariants
C         C0_123 = C0(1,2,3) = C0(p,q)     scalar three point
C         C0_124 = C0(1,2,4) = C0(p,q+l)  functions in PV notation
C         C0_134 = C0(1,3,4) = C0(p+q,l)
C         C0_234 = C0(2,3,4) = C0(q,l)
C         Cij_123(n,m) = C_nm(1,2,3) ....    higher C_nm form factors
C                                            as in tens_red3
c         D0 = D0(p,q,l)                  scalar four point function
c
c OUTPUT: Dij(n,m) = D_nm                    form factors in the tensor
c                                            integrals a la PV
c         nm = 11, 12, 13                    ff"s for D_mu
c         nm = 21, 22, 23, 24, 25, 26, 27    ff"s for D_munu
c         nm = 31, 32, 33, ..., 39, 310, 311, 312  ff"s for D_mu,nu,rho
c
      real*8 f1, f2, f3, Xm1(3,3), deter, m
      complex*16 R(1:3),res(1:3)

      m = 0d0
      f1 = -psq
      f2 = -qsq-2*pq
      f3 = -lsq-2*ql-2*pl

      deter = psq*qsq*lsq-psq*ql**2-pq**2*lsq+2*pq*pl*ql-pl**2*qsq
      Xm1(1,1) = (qsq*lsq-ql**2)/deter
      Xm1(1,2) = -(pq*lsq-pl*ql)/deter
      Xm1(1,3) = (pq*ql-pl*qsq)/deter
      Xm1(2,1) = Xm1(1,2)
      Xm1(2,2) = (psq*lsq-pl**2)/deter
      Xm1(2,3) = (-psq*ql+pq*pl)/deter
      Xm1(3,1) = Xm1(1,3)
      Xm1(3,2) = Xm1(2,3)
      Xm1(3,3) = -(-psq*qsq+pq**2)/deter


      R(1) = 0.5d0*(C0_134-C0_234+f1*D0)
      R(2) = 0.5d0*(C0_124-C0_134+f2*D0)
      R(3) = 0.5d0*(C0_123-C0_124+f3*D0)
      call prod_mat_col3(Xm1,R(1),res)
      Dij(1,1) = res(1)
      Dij(1,2) = res(2)
      Dij(1,3) = res(3)


      Dij(2,7) = (C0_234+2*m**2*D0-f1*Dij(1,1)-f2*Dij(1,2)
     #     -f3*Dij(1,3))/2

c 30-31-32
      R(1) = 0.5d0*(Cij_134(1,1)+C0_234+f1*Dij(1,1))-Dij(2,7)
      R(2) = 0.5d0*(Cij_124(1,1)-Cij_134(1,1)+f2*Dij(1,1))
      R(3) = 0.5d0*(Cij_123(1,1)-Cij_124(1,1)+f3*Dij(1,1))
      call prod_mat_col3(Xm1,R(1),res)
      Dij(2,1) = res(1)
      Dij(2,4) = res(2)
      Dij(2,5) = res(3)

c 33-34-35
      R(1) = 0.5d0*(Cij_134(1,1)-Cij_234(1,1)+f1*Dij(1,2))
      R(2) = 0.5d0*(Cij_124(1,2)-Cij_134(1,1)+f2*Dij(1,2))-Dij(2,7)
      R(3) = 0.5d0*(Cij_123(1,2)-Cij_124(1,2)+f3*Dij(1,2))
      call prod_mat_col3(Xm1,R(1),res)
      Dijp(2,4) = res(1)
      Dij(2,2) = res(2)
      Dij(2,6) = res(3)


c 36-37-38
      R(1) = 0.5d0*(Cij_134(1,2)-Cij_234(1,2)+f1*Dij(1,3))
      R(2) = 0.5d0*(Cij_124(1,2)-Cij_134(1,2)+f2*Dij(1,3))
      R(3) = 0.5d0*(-Cij_124(1,2)+f3*Dij(1,3))-Dij(2,7)
      call prod_mat_col3(Xm1,R(1),res)
      Dijp(2,5) = res(1)
      Dijp(2,6) = res(2)
      Dij(2,3) = res(3)


      Dij(3,11) = 0.5d0*(-0.5d0*C0_234 + m**2*Dij(1,1)
     #     - 0.5d0*(f1*Dij(2,1)+f2*Dij(2,4)+f3*Dij(2,5)))
      Dij(3,12) = 0.5d0*(0.5d0*Cij_234(1,1) + m**2*Dij(1,2)
     #     - 0.5d0*(f2*Dij(2,2)+f1*Dij(2,4)+f3*Dij(2,6)))
      Dij(3,13) = 0.5d0*(0.5d0*Cij_234(1,2) + m**2*Dij(1,3)
     #     - 0.5d0*(f1*Dij(2,5)+f2*Dij(2,6)+f3*Dij(2,3)))


c 41-42-43
      R(1) = 0.5d0*(Cij_134(2,1) - C0_234 + f1*Dij(2,1))-2*Dij(3,11)
      R(2) = 0.5d0*(Cij_124(2,1) - Cij_134(2,1) + f2*Dij(2,1))
      R(3) = 0.5d0*(Cij_123(2,1) - Cij_124(2,1) + f3*Dij(2,1))
      call prod_mat_col3(Xm1,R(1),res)
      Dij(3,1) = res(1)
      Dij(3,4) = res(2)
      Dij(3,5) = res(3)

c 50-51-52
      R(1) = 0.5d0*(Cij_134(2,1) - Cij_234(2,1) + f1*Dij(2,2))
      R(2) = 0.5d0*(Cij_124(2,2) - Cij_134(2,1) + f2*Dij(2,2))
     #     -2*Dij(3,12)
      R(3) = 0.5d0*(Cij_123(2,2) - Cij_124(2,2) + f3*Dij(2,2))
      call prod_mat_col3(Xm1,R(1),res)
      Dij(3,6) = res(1)
      Dij(3,2) = res(2)
      Dij(3,8) = res(3)

c 56-57-58
      R(1) = 0.5d0*(Cij_134(2,2) - Cij_234(2,2) + f1*Dij(2,3))
      R(2) = 0.5d0*(Cij_124(2,2) - Cij_134(2,2) + f2*Dij(2,3))
      R(3) = 0.5d0*(-Cij_124(2,2) + f3*Dij(2,3))-2*Dij(3,13)
      call prod_mat_col3(Xm1,R(1),res)
      Dij(3,7) = res(1)
      Dij(3,9) = res(2)
      Dij(3,3) = res(3)

c 44-45-46
      R(1) = 0.5d0*(Cij_134(2,1) + Cij_234(1,1) + f1*Dij(2,4))
     #     - Dij(3,12)
      R(2) = 0.5d0*(Cij_124(2,3) - Cij_134(2,1) + f2*Dij(2,4))
     #     - Dij(3,11)
      R(3) = 0.5d0*(Cij_123(2,3) - Cij_124(2,3) + f3*Dij(2,4))
      call prod_mat_col3(Xm1,R(1),res)
      Dijp(3,4) = res(1)
      Dijp(3,6) = res(2)
      Dij(3,10) = res(3)


      t(1) = (Cij_134(2,3)+f1*Dij(2,5)+Cij_234(1,2))/2-Dij(3,13)
      t(2) = (Cij_134(2,3)+f1*Dij(2,6)-Cij_234(2,3))/2
      t(3) = (Cij_134(2,4)+f1*Dij(2,7)-Cij_234(2,4))/2

      t(4) = (Cij_124(2,3)-Cij_134(2,3)+f2*Dij(2,5))/2
      t(5) = (Cij_124(2,2)-Cij_134(2,3)+f2*Dij(2,6))/2-Dij(3,13)
      t(6) = (Cij_124(2,4)-Cij_134(2,4)+f2*Dij(2,7))/2

      t(7) = (           -Cij_124(2,3)+f3*Dij(2,5))/2-Dij(3,11)
      t(8) = (           -Cij_124(2,2)+f3*Dij(2,6))/2-Dij(3,12)
      t(9) = (Cij_123(2,4)-Cij_124(2,4)+f3*Dij(2,7))/2

      Dijp(3,5) = Xm1(1,1)*t(1) + Xm1(1,2)*t(4) + Xm1(1,3)*t(7)
      Dijp(3,10)= Xm1(2,1)*t(1) + Xm1(2,2)*t(4) + Xm1(2,3)*t(7)
      Dijp(3,7) = Xm1(3,1)*t(1) + Xm1(3,2)*t(4) + Xm1(3,3)*t(7)

      Dijpp(3,10) = Xm1(1,1)*t(2) + Xm1(1,2)*t(5) + Xm1(1,3)*t(8)
      Dijp(3,8) = Xm1(2,1)*t(2) + Xm1(2,2)*t(5) + Xm1(2,3)*t(8)
      Dijp(3,9) = Xm1(3,1)*t(2) + Xm1(3,2)*t(5) + Xm1(3,3)*t(8)

      Dijp(3,11)= Xm1(1,1)*t(3) + Xm1(1,2)*t(6) + Xm1(1,3)*t(9)
      Dijp(3,12)= Xm1(2,1)*t(3) + Xm1(2,2)*t(6) + Xm1(2,3)*t(9)
      Dijp(3,13)= Xm1(3,1)*t(3) + Xm1(3,2)*t(6) + Xm1(3,3)*t(9)

c (giuseppe bozzi)
c WARNING: Previously unused Dij(1,7-13) variables are used
c to define D00ij and D0000 functions!
c In PV notation we have:
c Dij(1,7)=D416, Dij(1,8)=D417, Dij(1,9)=D418
c Dij(1,10)=D419, Dij(1,11)=D420, Dij(1,12)=D421, D(1,13)=D422
      Dij(1,7)=(C0_234-f1*Dij(3,1)-f2*Dij(3,4)-f3*Dij(3,5))/6d0
      Dij(1,8)=(Cij_234(2,1)-f1*Dij(3,6)-f2*Dij(3,2)-f3*Dij(3,8))/6d0
      Dij(1,9)=(Cij_234(2,2)-f1*Dij(3,7)-f2*Dij(3,9)-f3*Dij(3,3))/6d0
      Dij(1,10)=-(Cij_234(1,1)+f1*Dij(3,4)+f2*Dij(3,6)+f3*Dij(3,10))/6d0
      Dij(1,11)=-(Cij_234(1,2)+f1*Dij(3,5)+f2*Dij(3,10)+f3*Dij(3,7))/6d0
      Dij(1,12)=(Cij_234(2,3)-f1*Dij(3,10)-f2*Dij(3,8)-f3*Dij(3,9))/6d0
      Dij(1,13)=(Cij_234(2,4)+1/12d0-psq*Dij(1,7)-qsq*Dij(1,8)-
     \     lsq*Dij(1,9)-2d0*pq*Dij(1,10)-2d0*pl*Dij(1,11)-
     \     2d0*ql*Dij(1,12))/6d0

      if (.false.) then
         print*," D34 ",Dij(3,4),Dijp(3,4),Dij(3,4)/Dijp(3,4)
         print*," D36 ",Dij(3,6),Dijp(3,6),Dij(3,6)/Dijp(3,6)
         print*," D35 ",Dij(3,5),Dijp(3,5),Dij(3,5)/Dijp(3,5)
         print*," D37 ",Dij(3,7),Dijp(3,7),Dij(3,7)/Dijp(3,7)
         print*," D38 ",Dij(3,8),Dijp(3,8),Dij(3,8)/Dijp(3,8)
         print*," D39 ",Dij(3,9),Dijp(3,9),Dij(3,9)/Dijp(3,9)
         print*," D310 ",Dij(3,10),Dijp(3,10),Dij(3,10)/Dijp(3,10)
         print*," D310pp ",Dij(3,10),Dijpp(3,10),Dijpp(3,10)/Dij(3,10)
         print*," D311 ",Dij(3,11),Dijp(3,11),Dij(3,11)/Dijp(3,11)
         print*," D312 ",Dij(3,12),Dijp(3,12),Dij(3,12)/Dijp(3,12)
         print*," D313 ",Dij(3,13),Dijp(3,13),Dij(3,13)/Dijp(3,13)
      endif

      end


      subroutine prod_mat_col3(M,col,res)
      implicit none
      real * 8 M(1:3,1:3)
      complex * 16 col(1:3), res(1:3)
      integer i,j
      do i=1,3
         res(i) = 0d0
      enddo
      do i=1,3
         do j=1,3
            res(i) = res(i) + M(i,j)*col(j)
         enddo
      enddo
      end
