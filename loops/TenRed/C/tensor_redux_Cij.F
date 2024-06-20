CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                tens_red3 = 3-point tensors
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine tens_red3(psq, qsq, pq_sq, B0_23, B0_13, B0_12, 
     #                     C0, Cij)
      implicit none
      real * 8  psq, qsq, pq_sq
      complex*16 B0_23, B0_13, B0_12, C0, Cij(2,4), Cijp(2,4)
c
c  determine the Passarino-Veltman tensor decomposition for the three-point
c  tensor integrals
c
c                       d^4k           1;  k_mu;   k_mu k_nu
c C0;C_mu;C_munu =Int -------  -----------------------------------------
c                     (2pi)^4    [k^2-m^2][(k+p)^2-m^2][(k+p+q)^2-m^2] 
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
C INPUT:  psq, qsq, pq_sq          external invariants: p^2, q^2, (p+q)^2
C         B0, C0                   4 scalar integrals; the 3 B0 are, 
c                                  in PV notation:
c         B0(1) = B0(1,2) = B0(p)  B_0 function with subtraction of 
c         B0(2) = B0(2,3) = B0(q)  divergent term
c         B0(3) = B0(1,3) = B0(p+q)
c
c OUTPUT: Cij(n,m) = C_nm          form factors in the tensor integrals
c         n=1,2; m=1,2,3,4         a la PV
c
      real*8 f1, f2, Xm1(2,2), deter, pq
      complex*16 B1_12,B1_13,B1_23,R(2), res(2)
      real * 8 m

      pq = (pq_sq - psq - qsq)/2
      m = 0d0
      f1 = -psq
      f2 = -qsq-2*pq
      

      deter = psq*qsq-pq**2
      Xm1(1,1) = qsq/deter
      Xm1(1,2) = -pq/deter
      Xm1(2,1) = Xm1(1,2)
      Xm1(2,2) = psq/deter

c      B0_12 = B0(1)
c      B0_23 = B0(2)
c      B0_13 = B0(3)
      
      B1_12 = -0.5d0*B0_12
      B1_13 = -0.5d0*B0_13
      B1_23 = -0.5d0*B0_23
      
c 1-2
      R(1) = 0.5d0*(B0_13-B0_23+f1*C0)
      R(2) = 0.5d0*(B0_12-B0_13+f2*C0)
      call prod_mat_col2(Xm1,R(1),res)
      Cij(1,1) = res(1)
      Cij(1,2) = res(2)
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c c"e" da sistemare per poli!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Cij(2,4) = 0.5d0*( + 0.5d0*B0_23 + m**2*C0
     #     - 0.5d0*f1*Cij(1,1) - 0.5d0*f2*Cij(1,2)) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 1/4 comes from keeping only the FINITE parts
     #     +1d0/4
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c 3-5
      R(1) = 0.5d0*(B1_13 + B0_23 + f1*Cij(1,1))-Cij(2,4)
      R(2) = 0.5d0*(B1_12 - B1_13 + f2*Cij(1,1))
      call prod_mat_col2(Xm1,R(1),res)
      Cij(2,1) = res(1)
      Cij(2,3) = res(2)

c 4-6
      R(1) = 0.5d0*(B1_13 - B1_23 + f1*Cij(1,2))
      R(2) = 0.5d0*( - B1_13 + f2*Cij(1,2)) - Cij(2,4)
      call prod_mat_col2(Xm1,R(1),res)
      Cijp(2,3) = res(1)
      Cij(2,2) = res(2)

      if (.false.) then
         print*," C23 ",Cij(2,3),Cijp(2,3),Cij(2,3)/Cijp(2,3)
      endif

      end


      subroutine prod_mat_col2(M,col,res)
      implicit none
      real * 8 M(1:2,1:2)
      complex * 16 col(1:2), res(1:2)
      integer i,j
      do i=1,2 
         res(i) = 0d0
      enddo
      do i=1,2
         do j=1,2
            res(i) = res(i) + M(i,j)*col(j)
         enddo
      enddo
      end


