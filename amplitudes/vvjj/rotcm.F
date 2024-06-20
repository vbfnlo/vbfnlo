cccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine BOOSTCM(P1,P2,P3,P4, T, Tinv)
c    Boost from the lab frame to the CM frame and a rotation to the z axis.
c    input:
c    Pi(0:3)      : vector-bosons momenta in the lab frame
c         i = 1,2 : particles before scattering
c         i = 3,4 : particles after scattering
c
c   output:
c     T(0:3,0:3)    : Transformation matrix (Boost and rotation)
c     Tinv(0:3,0:3) : Inverse matrix of T 
c

       implicit none
       DOUBLE PRECISION rZero, rOne
       PARAMETER( rZero=0.0d0, rOne=1.0d0)       
       
       REAL*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3)
       REAL*8 T(0:3,0:3), Tinv(0:3,0:3)
       
       REAL*8 Q1(0:3),Q2(0:3),Q3(0:3),Q4(0:3)
       REAL*8 R(0:3), T2(0:3,0:3)     
       REAL*8 X, Y1, Y2, R_SQ
       REAL*8 BETA(3), N1(3), N3(3), O1(3)       
       REAL*8 GAMMA                  
       REAL*8 L(0:3,0:3), RA(0:3,0:3), Linv(0:3,0:3), RAinv(0:3,0:3)      
       REAL*8 Q1_SQ, Q3_SQ, Q1_INV, Q3_INV, W1, W2, O_INV
       
       INTEGER I, J, mu ,K

       logical dswitch, gswitch
       parameter(dswitch=.false. ,gswitch=.false.) ! or .true., it will help you to test things :)

       REAL*8 F1(0:3), F2(0:3), F3(0:3), F4(0:3)
       REAL*8 UNI(0:3,0:3), UNI2(0:3,0:3), uni3(0:3,0:3)       
       REAL*8 QC1(0:3),QC2(0:3),QC3(0:3),QC4(0:3)
       REAL*8 Z1, Z2, Z3,Z4
       REAL*8 ZB1, ZB2,ZB3,ZB4
       REAL*8 ZBB1,ZBB2,ZBB3,ZBB4
       
#include "gmunu.inc"

       X = 0D0
       L(0:3,0:3) = 0.0d0
       RA(0:3,0:3) = 0.0d0
       T(0:3,0:3) = 0.0d0
       Tinv(0:3,0:3) = 0.0d0
       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCC     BOOST    CCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     
       R(0) = p1(0)+p2(0)
       DO K = 1,3
          R(k) = - (p1(k)+p2(k))          
       ENDDO
       
       DO I = 1,3
          BETA(I) = R(I)/R(0)
          X = X + BETA(I)**2
       ENDDO
    
       if (X.GT.1D-16.AND.X.LT.(1D0-1D-12)) then
          GAMMA = 1D0/DSQRT(1D0-X)
          Y1 = 1D0/X       
          Y2 = GAMMA-1D0

          L = RESHAPE( (/gamma , -gamma*beta(1) , -gamma*beta(2) , -gamma*beta(3),
     &       -gamma*beta(1) , 1d0 + Y1*Y2*beta(1)**2, Y1*Y2*beta(1)*beta(2), Y1*Y2*beta(1)*beta(3),
     &       -gamma*beta(2) , Y1*Y2*beta(1)*beta(2), 1d0 + Y1*Y2*beta(2)**2, Y1*Y2*beta(2)*beta(3),
     &       -gamma*beta(3) , Y1*Y2*beta(1)*beta(3), Y1*Y2*beta(2)*beta(3), 1d0 + Y1*Y2*beta(3)**2/),(/4,4/) )     
                
          do i = 0,3
           q1(i) = L(i,0)*p1(0)-L(i,1)*p1(1)-L(i,2)*p1(2)-L(i,3)*p1(3)
           q2(i) = L(i,0)*p2(0)-L(i,1)*p2(1)-L(i,2)*p2(2)-L(i,3)*p2(3)
           q3(i) = L(i,0)*p3(0)-L(i,1)*p3(1)-L(i,2)*p3(2)-L(i,3)*p3(3)
           q4(i) = L(i,0)*p4(0)-L(i,1)*p4(1)-L(i,2)*p4(2)-L(i,3)*p4(3)
          enddo

       else
          DO i = 0,3
             q1(i) = p1(i)
             q2(i) = p2(i)
             q3(i) = p3(i)
             q4(i) = p4(i)
          ENDDO
       R_sq = R(0)**2-R(1)**2-R(2)**2-R(3)**2
       If(R_sq.le.0d0)
     *    write(6,1000) R, R_sq
          print*, "W1 =", P1(0:3)
          print*, "W2 = ", P2(0:3)
       endif
 1000  FORMAT (" The reference vector ",4G12.3," is not timelike. "/
     1        " R1**2 = ",G12.3)

       if (dswitch) then
C dot product calculation (invariance)
          z1 = p1(0)**2 - (p1(1)**2+p1(2)**2+p1(3)**2)
          z2 = p2(0)**2 - (p2(1)**2+p2(2)**2+p2(3)**2)
          z3 = p3(0)**2 - (p3(1)**2+p3(2)**2+p3(3)**2)
          z4 = p4(0)**2 - (p4(1)**2+p4(2)**2+p4(3)**2)
          
          zb1 = q1(0)**2 - (q1(1)**2+q1(2)**2+q1(3)**2)
          zb2 = q2(0)**2 - (q2(1)**2+q2(2)**2+q2(3)**2)
          zb3 = q3(0)**2 - (q3(1)**2+q3(2)**2+q3(3)**2)
          zb4 = q4(0)**2 - (q4(1)**2+q4(2)**2+q4(3)**2) 
          
          zbb1 = (p1(0)+p2(0))**2 - ((p1(1)+p2(1))**2+(p1(2)+p2(2))**2+(p1(3)+p2(3))**2)
          zbb2 = (p3(0)+p4(0))**2 - ((p3(1)+p4(1))**2+(p3(2)+p4(2))**2+(p3(3)+p4(3))**2)
          zbb3 = (q1(0)+q2(0))**2 - ((q1(1)+q2(1))**2+(q1(2)+q2(2))**2+(q1(3)+q2(3))**2)
          zbb4 = (q3(0)+q4(0))**2 - ((q3(1)+q4(1))**2+(q3(2)+q4(2))**2+(q3(3)+q4(3))**2)

          print*, "before boosted momenta, invariant sq initial = ", ZBB1
          print*, "before boosted momenta, invariant sq final = ", ZBB2           
          print*, "boosted momenta, invariant sq initial = ", ZBB3
          print*, "boosted momenta, invariant sq final = ", ZBB4           

          print*, "boosted momenta, invariant sq initial for 1 = ", z1
          print*, "boosted momenta, invariant sq final = ", zb1 
          print*, "boosted momenta, invariant sq initial for 2 = ", z2
          print*, "boosted momenta, invariant sq final = ", zb2 
          print*, "boosted momenta, invariant sq initial for 3 = ", z3
          print*, "boosted momenta, invariant sq final = ", zb3 
          print*, "boosted momenta, invariant sq initial for 4 = ", z4
          print*, "boosted momenta, invariant sq final = ", zb4           
          print*, "  ========        ========       ======" 
          print*, "p1 before boost = ", p1(0:3)
          print*, "p2 before boost = ", p2(0:3)
          print*, "p3 before boost = ", p3(0:3)
          print*, "p4 before boost = ", p4(0:3)        
          print*, "p1 after boost  = ", q1(0:3)
          print*, "p2 after boost  = ", q2(0:3)
          print*, "p3 after boost  = ", q3(0:3)
          print*, "p4 after boost  = ", q4(0:3)          
       endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCC   ROTATION   CCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      q1_sq = 0d0
      q3_sq = 0d0   
      
      DO i = 1,3
        q1_sq = q1_sq + q1(i)**2
        q3_sq = q3_sq + q3(i)**2      
      ENDDO

      q1_inv = 1d0/DSQRT(q1_sq)
      q3_inv = 1d0/DSQRT(q3_sq)
      
      DO i = 1,3
        n1(i) = q1(i)*q1_inv
        n3(i) = q3(i)*q3_inv
      ENDDO
      
      W1 = 0d0
      W2 = 0d0
      DO i = 1,3
         W1 = W1 + n1(i)*n3(i)
      ENDDO

      DO i = 1,3
         o1(i) = n3(i) - n1(i)*W1
         W2 = W2 + o1(i)**2
      ENDDO

      o_inv = 1d0/DSQRT(W2)
     
      Ra = RESHAPE( (/rOne ,rZero,rZero,rZero,
     &           rZero, o1(1)*o_inv, o1(3)*o_inv*n1(2)-o1(2)*o_inv*n1(3),n1(1),
     &           rZero, o1(2)*o_inv, o1(1)*o_inv*n1(3)-o1(3)*o_inv*n1(1),n1(2),
     &           rZero, o1(3)*o_inv, o1(2)*o_inv*n1(1)-o1(1)*o_inv*n1(2),n1(3)/), (/4,4/))
   
           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCC TRANSFORMATION MATRIX   CCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
c
C The total transformation matrix is given by
c   T = Ra(rot matrix)*LT(Lorentz matrix)
C those were defined before
c T^(mu nu)

      DO i = 0,3
         DO k = 0,3
            T(i,k) = 0
            DO j = 0,3
              T(i,k) = T(i,k) + Ra(i,j)*L(j,k)
            ENDDO
         ENDDO
      ENDDO

          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC TRANSFORMATION MATRIX INVERSE   CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      

      RAinv = Transpose(Ra)
      
      T2 = Transpose(L)
      
      Linv(0,1) = - T2(0,1)
      Linv(1,0) = - T2(1,0)
      Linv(0,2) = - T2(0,2)
      Linv(2,0) = - T2(2,0)
      Linv(0,3) = - T2(0,3)
      Linv(3,0) = - T2(3,0)
      Linv(0,0) = T2(0,0)
      
      do i = 1,3     
       do k = 1,3
        Linv(i,k) = T2(i,k)
       enddo
      enddo

      DO i = 0,3
         DO k = 0,3
            Tinv(i,k) = 0            
            DO j = 0,3
               Tinv(i,k) = Tinv(i,k) + Linv(i,j)*Rainv(j,k)               
            ENDDO
         ENDDO
      ENDDO
                         
      if(gswitch) then 
      
       do i = 0,3
         do j = 0,3
           uni(i,j) = 0
           do k = 0,3
              uni(i,j) = ra(i,0)*RAinv(0,j)+ra(i,1)*RAinv(1,j)+ra(i,2)*rainv(2,j)+ra(i,3)*rainv(3,j)
              uni2(i,j) = l(i,0)*linv(0,j)+l(i,1)*linv(1,j)+l(i,2)*linv(2,j)+l(i,3)*linv(3,j)  
              uni3(i,j) = tinv(i,0)*t(0,j)+tinv(i,1)*t(1,j)+tinv(i,2)*t(2,j)+tinv(i,3)*t(3,j)  
           enddo
         enddo
       enddo 
       
       do i = 0,3
           f1(i) = T(i,0)*p1(0)-T(i,1)*p1(1)-T(i,2)*p1(2)-T(i,3)*p1(3)
           f2(i) = T(i,0)*p2(0)-T(i,1)*p2(1)-T(i,2)*p2(2)-T(i,3)*p2(3)
           f3(i) = T(i,0)*p3(0)-T(i,1)*p3(1)-T(i,2)*p3(2)-T(i,3)*p3(3)
           f4(i) = T(i,0)*p4(0)-T(i,1)*p4(1)-T(i,2)*p4(2)-T(i,3)*p4(3)
       enddo   



       qc1(0) = +Tinv(0,0)*f1(0)+Tinv(0,1)*f1(1)+Tinv(0,2)*f1(2)+Tinv(0,3)*f1(3)
       qc2(0) = +Tinv(0,0)*f2(0)+Tinv(0,1)*f2(1)+Tinv(0,2)*f2(2)+Tinv(0,3)*f2(3)
       qc3(0) = +Tinv(0,0)*f3(0)+Tinv(0,1)*f3(1)+Tinv(0,2)*f3(2)+Tinv(0,3)*f3(3)
       qc4(0) = +Tinv(0,0)*f4(0)+Tinv(0,1)*f4(1)+Tinv(0,2)*f4(2)+Tinv(0,3)*f4(3)          
      
       do i = 1,3
         qc1(i) = -Tinv(i,0)*f1(0)-Tinv(i,1)*f1(1)-Tinv(i,2)*f1(2)-Tinv(i,3)*f1(3)
         qc2(i) = -Tinv(i,0)*f2(0)-Tinv(i,1)*f2(1)-Tinv(i,2)*f2(2)-Tinv(i,3)*f2(3)
         qc3(i) = -Tinv(i,0)*f3(0)-Tinv(i,1)*f3(1)-Tinv(i,2)*f3(2)-Tinv(i,3)*f3(3)
         qc4(i) = -Tinv(i,0)*f4(0)-Tinv(i,1)*f4(1)-Tinv(i,2)*f4(2)-Tinv(i,3)*f4(3)          
       enddo       

       print*, "p1 after boost and rotation = ", f1(0:3)
       print*, "p2 after boost and rotation = ", f2(0:3)
       print*, "p3 after boost and rotation = ", f3(0:3)
       print*, "p4 after boost and rotation = ", f4(0:3)
       print*, "      "
       print*, "p1 before boost = ", p1(0:3)
       print*, "p2 before boost = ", p2(0:3)
       print*, "p3 before boost = ", p3(0:3)
       print*, "p4 before boost = ", p4(0:3)        
       print*, "      "       
       print*, "p1 after boosted to cm and boosted back", qc1(0:3)
       print*, "p2 after boosted to cm and boosted back", qc2(0:3)
       print*, "p3 after boosted to cm and boosted back", qc3(0:3)
       print*, "p4 after boosted to cm and boosted back", qc4(0:3)           
          
      print*, "rotation unit matrix"
      print*, uni(0,:)
      print*, uni(1,:)
      print*, uni(2,:)
      print*, uni(3,:)
      
      print*, "lorentz unit matrix"
      print*, uni2(0,:)
      print*, uni2(1,:)
      print*, uni2(2,:)
      print*, uni2(3,:) 
      
      print*, "unit matrix"
      print*, uni3(0,:)
      print*, uni3(1,:)
      print*, uni3(2,:)
      print*, uni3(3,:)              
       
      endif 
     
      END
       
