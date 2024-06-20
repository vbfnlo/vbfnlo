************************************************************************
* Operations with matrices:
* The LU Decomposition, Determinant and inverse of a matrix using LU decomposition are M. Rauch's subroutines.
* The QR col pivot decomposition, Householder vector and all other subroutines required, including the explicit inverse 2x2
* and 3x3 matrices are G. Perez' subroutines.
* P.S.: The QR col pivot decomposition works, but I couldn't manage to return back the pivot. So, I took it out.
****************************************************************         
        SUBROUTINE MAXEigenvalue(A,m,maxEigenv) !9-->m
        
        IMPLICIT NONE
        
        INTEGER I,J,K,L
        INTEGER N, M, P
        PARAMETER(N=4)
        DOUBLE COMPLEX A(M,M)
        DOUBLE COMPLEX maxEigenv
        
        DOUBLE COMPLEX Q(M,M), R(M,M)
        DOUBLE COMPLEX Aeigen(M,M,0:N)
        DOUBLE COMPLEX A1(M,M)
        DOUBLE PRECISION AArray(M)
        DOUBLE PRECISION REST
        
        INTEGER ier
        DOUBLE PRECISION ACC
        
        Aeigen(1:M,1:M,0) = A(1:M,1:M) 
        A1(1:M,1:M) = A(1:M,1:M)
        Aeigen(1:M,1:M,1:N) = (0.0D0,0.0D0)
        maxEigenv = (0.0D0,0.0D0)
        
        k = 0        
        do while(k.lt.N)  !optimal 10         
         Q(1:M,1:M) = (0.0D0,0.0D0)
         R(1:M,1:M) = (0.0D0,0.0D0)
        
         CALL CQRpivot(Aeigen(1,1,k),M,M,Q,R)
         k = k+1
         do i =1,M
          do j = 1,M
           do l = 1,M
             Aeigen(i,j,k) = Aeigen(i,j,k) + R(i,l)*Q(l,j)
           enddo
          enddo 
         enddo     
****correct this please
!          rest = 0.0d0                  
!          do i = 1,9
!           do j = 1,9
!            if (real(Aeigen(i,j,k-1)).lt.1.0d-6) then
!             rest = rest + abs(real(Aeigen(i,j,k)) - real(Aeigen(i,j,k-1)))
!            else                      
!             rest = rest + (1.0d0 - (real(Aeigen(i,j,k))/ real(Aeigen(i,j,k-1))))         
!            endif
!           enddo
!          enddo 
         
!          if(abs(rest/81.0d0).le.1.0d-2.or.k.eq.N) then
!            maxEigenv = max(Real(Aeigen(1,1,k)),Real(Aeigen(2,2,k)),Real(Aeigen(3,3,k)),Real(Aeigen(4,4,k)),Real(Aeigen(5,5,k)),Real(Aeigen(6,6,k)),Real(Aeigen(7,7,k)),Real(Aeigen(8,8,k)),Real(Aeigen(9,9,k)))         
            do P=1,M
             AArray(P) = Real(Aeigen(P,P,k))
            enddo
            maxEigenv = maxval(AArray)
!             exit
!          endif         
        enddo         

        END
****************************************************************         
        SUBROUTINE settoCU(n,m,A)
        IMPLICIT NONE
        
        INTEGER N,M
        DOUBLE COMPLEX A(N,M)
        INTEGER I,J
        DOUBLE COMPLEX cZERO, cONE
        PARAMETER(cZERO = (0.0D0,0.0D0), cONE=(1.0D0,0.0D0))
        
        DO I = 1,N
         DO J = 1,M
          A(I,J) = cZERO
         ENDDO
         A(I,I) = cONE
        ENDDO
        
        END
        
****************************************************************         
        SUBROUTINE settoRU(n,m,A)
        IMPLICIT NONE
        
        INTEGER N,M
        DOUBLE PRECISION A(N,M)
        INTEGER I,J
        DOUBLE PRECISION rZERO, rONE
        PARAMETER(rZERO = 0.0D0, rONE=1.0D0)
        
        DO I = 1,N
         DO J = 1,M
          A(I,J) = rZERO
         ENDDO
         A(I,I) = rONE
        ENDDO
        
        END        
        
****************************************************************         
        SUBROUTINE CQREq(A,B,n,m,x)
* It calculates x for x.A = B, for complex values.
* A(n,m), B(1,m), x(1,n)
* and A is decomposed using QR decomposition with Householder vector.        
        IMPLICIT NONE
        
        INTEGER N,M
        DOUBLE COMPLEX A(N,M)
        DOUBLE COMPLEX B(M), X(N), Y(N)
        DOUBLE COMPLEX Q(N,N), R(N,M),Qdag(N,N)
        
        DOUBLE COMPLEX TEMP 
        
        DOUBLE COMPLEX cZERO,cONE, cIMAG
        PARAMETER(cZERO = (0.0D0,0.0D0), cONE=(1.0D0,0.0D0), cIMAG=(0.0D0,1.0D0))
        DOUBLE PRECISION rZERO, rONE
        PARAMETER(rZERO=0.0D0, rONE=1.0D0)
        
        INTEGER I,J,K,L
        INTEGER NN,MM
        
******************        
        y(1:N) = cZero
        x(1:N) = cZero
        temp = cZero
        
******************                
* Does a QR decomposition: A = QR
* where Q is a unitary matrix and R and upper triangular matrix.
        CALL CQRpivot(A,N,M,Q,R)
* 
* xA = B  is calculated as xQR=B, then xQ=Y and YR=B        
*
* Y calculation, using YR=B
        y(1) = b(1)/R(1,1)        
        do k = 2, n
         do i = 1, k-1
          y(k) = y(k)+ y(i)*R(i,k)
         enddo
         y(k) = (b(k)-y(k))/R(k,k)
        enddo        
*
* X calculation, using xQ=Y, x=YQdag
        Qdag = transpose(conjg(Q))
        
        do i = 1, n
         do j = 1,n
          x(j) = x(j)+y(i)*Qdag(i,j)
         enddo
        enddo
        
        END
        
************************************************************************        
        SUBROUTINE CQRpivot(A,N,M,Q,R)!,Piv)  ! complex QR calculation
*   A(N,M) matrix to decompose
*   Q(N,N) unitary matrix, Qdagger*Q = 1
*   R(N,M) upper triangular matrix 
*   A(n,m) P(m,m) =  Q(n,n,) R(n,m)

        IMPLICIT NONE
    
        INTEGER N,M   ! N=rows, M=Columns
        DOUBLE COMPLEX A(N,M)
        DOUBLE COMPLEX Q(N,N), R(N,M)
        DOUBLE COMPLEX X(N,M), RR(N,M),QQ(N,N)
        DOUBLE COMPLEX H(N,N), Hi(N,N)
        DOUBLE COMPLEX V(N), W, CUnit(N,N)

        DOUBLE COMPLEX TEMP
    
        INTEGER I,J,K,L
        INTEGER NN,MM
        
        Logical tswitch
        PARAMETER(tswitch = .false.)
        DOUBLE COMPLEX T(N,M), Htest(N,N)
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!         INTEGER PERM, PIV(N)
!         DOUBLE COMPLEX  P(M,M), PTRANS(M,M)        
!         DOUBLE PRECISION RNORM(M),MAXNORM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

*************        
        call settoCU(n,n,Q)
        call settoCU(n,n,CUnit)
        do i =1,n
         do j =1,m
         R(i,j) = A(i,j)
         enddo
        enddo               
*************         
        
* NB: The Householder transformation will set to zero one column at time
* therefore, R = Hm*Hm-1*....*H1*A and Q = H1*H2*...Hm
        j = 1  
        do while (j.lt.m)
         mm = m-j+1
         nn = n-j+1      

* Get Householder transformation matrix: H,
* and apply the transformation to R from the left and Q from the right.  
         X(1:N,1:M) = (0.0d0,0.0d0)         
         Hi(1:n,1:n) = (0.0d0,0.0d0)

         X(1:nn,j) = R(j:n,j)            
         call CHouse(X(1,j),nn,Hi(1:nn,1:nn))
         
         H(1:n,1:n) = (0.0d0,0.0d0)         
         if(j.eq.1) then
          H(1:n,1:n) = Hi(1:n,1:n)
         else
          do k = 1, j-1
          H(k,1:n) = (0.0d0,0.0d0)
          H(j-1:n,k) = (0.0d0,0.0d0)
          H(k,k) = (1.0d0,0.0d0)
          enddo                    
          H(j:n,j:n) = Hi(1:nn,1:nn)
         endif        

         RR(1:N,1:M) = (0.0d0,0.0d0)
         QQ(1:N,1:N) = (0.0d0,0.0d0)
         call Prod_Matrix(H,R,n,n,m,RR)
         call Prod_Matrix(Q,H,n,n,n,QQ)               
         
         R(1:N,1:M) = RR(1:N,1:M)
         Q(1:N,1:N) = QQ(1:N,1:N)      
         
         j = j+1
       enddo     

       if(tswitch) then
       print*, "          "
       print*, "R real = "
       do i = 1,N
        print*, real(R(i,:))
       enddo
       print*, "R imaginary= "
       do i = 1,N
        print*, dimag(R(i,:))
       enddo
       print*, "          "
       print*, "Q real= "
       do i = 1,N
        print*, real(Q(i,:))
       enddo
       print*, "Q imaginary= "
       do i = 1,N
        print*, dimag(Q(i,:))
       enddo
       
       call Prod_Matrix(Q,transpose(conjg(Q)),n,n,n,Htest)
       do i = 1,n
        write(*,*) "QQdag real ", real(Htest(i,:))
        write(*,*) "QQdag complex", dimag(Htest(i,:))
       enddo

       T(1:N,1:M) = (0.0d0,0.0d0)
       call Prod_Matrix(Q,R,n,n,m,T)                   
       do i = 1,n
        do j = 1,m
          write(*,*) "original real", real(A(i,j))," | ", "decomposition real", real(T(i,j))
          write(*,*) "original complex", dimag(A(i,j))," | ", "decomposition complex", dimag(T(i,j))          
        enddo
       enddo
       print*, " ================== "         
       endif

       END 
        
****************************************************************         
        SUBROUTINE CHouse(x,n, H)
** Calculates the Householder vector for a complex vector
** and its matrix Hv(n,n)        
        IMPLICIT NONE

        INTEGER N        
        DOUBLE COMPLEX X(N)
        DOUBLE PRECISION XNORM,VNorm, S, Xscale
        DOUBLE COMPLEX H(N,N), UNIT(N,N)
        DOUBLE COMPLEX V(N),U(N)
        DOUBLE COMPLEX W,Wup,Wdown,B
        
        DOUBLE COMPLEX cZERO
        PARAMETER(cZERO=(0.0D0,0.0D0))
        INTEGER I,J,K
        
        LOGICAL nzero
        nzero=.false.

        call settoCU(n,n,unit)
        call settoCU(n,n,H)        
        Xscale = 0.0d0
        B = cZero
        Wup = cZero
        Wdown = cZero
        
        do i = 2,N  !It sets nzero to true, if there's at least one component different from zero
         xscale = max(xscale,abs(x(i)))
        enddo 
        if(xscale.ne.0) then
          nzero=.true.
        endif
              
        if(.not.nzero) then
         v(1:N) = cZero
         H(1:N,1:N) = Unit(1:N,1:N)
        else   
         call NORM(n,x,xnorm)  
         
         if( x(1).ne.cZero) then
          if (Real(x(1)).ne.0.0d0) then
           s = - sign(1.0d0,Real(X(1)))
          else if (Real(x(1)).eq.0.0d0) then
           s = - sign(1.0d0,Dimag(X(1)))
          endif  
          v(1) = X(1) - S*XNORM*(X(1)/abs(X(1)))
          v(2:N) = X(2:N)
         else
          v(1) = XNORM
          v(2:N) = X(2:N)
         endif    
         
         do i =1,n
          B = B + conjg(V(i))*V(i)
         enddo
         do i = 1,n
          Wup = Wup + conjg(x(i))*v(i)
         enddo
         do i = 1,n
          Wdown = Wdown + conjg(v(i))*x(i)
         enddo
         W = (1.0d0 + (Wup/Wdown))/B         
         
         do i = 1,n
          do j = 1,n
           H(i,j) = unit(i,j) -W*v(i)*conjg(v(j))
          enddo
         enddo
        endif
        
        END 
        

****************************************************************         
        SUBROUTINE NORM(N,V,NORMV)
        IMPLICIT NONE
        
        INTEGER N
        DOUBLE COMPLEX V(N),SUMA
        DOUBLE PRECISION NORMV
        INTEGER I,J
        
        suma = (0.0d0,0.0d0)
        
        do i = 1,n
          suma = suma + V(i)*CONJG(V(i))
        enddo
        NORMV = SQRT((suma))
        
        END
        
****************************************************************         
        SUBROUTINE Prod_Matrix(A,B,N,M,L,PRO)
        IMPLICIT NONE
        
        INTEGER N,M,L        
        DOUBLE COMPLEX A(N,M), B(M,L), PRO(N,L)
        DOUBLE COMPLEX SUMA
        DOUBLE PRECISION Preal, Pim
        
        INTEGER I,J,K

        Pro(1:N,1:L) = (0.0d0,0.0d0)        
        
        do i = 1,N
         do j = 1,L
          Preal = 0.0d0
          Pim = 0.0d0
           do k = 1,M
           Pro(i,j) = Pro(i,j) + a(i,k)*b(k,j)          
           enddo
           if(abs(real(Pro(i,j))).lt.1.0d-12) Pro(i,j) = dcmplx(0.0d0,dimag(Pro(i,j)))
           if(abs(dimag(Pro(i,j))).lt.1.0d-12) Pro(i,j) = dcmplx(real(Pro(i,j)),0.0d0)
         enddo
        enddo              

        END
        
************************************************************************
        subroutine LUInverse(A,Inv) !why not n?
        
        IMPLICIT NONE
        
        DOUBLE COMPLEX cIMAG, cZERO, cONE
        PARAMETER( cIMAG=(0.0D0,1.0D0), cZERO=(0.0D0, 0.0D0), cONE=(1.0D0 ,0.0D0))        
        
        integer n, perm(9)
        parameter(n = 9 )
        double complex A(9,9), Inv(9,9)
        double complex LU(9,9)
        double complex K(9,9)
        double complex deter
        double complex cdet
        external cdet
        
        integer i,j
        
        LU(1:9,1:9) = cZero
        Inv(1:9,1:9) = cZero
       
        LU=A
        call Inverse(LU,Inv,n,perm)     
        end

************************************************************************
************************************************************************
* Inverse computes the inverse of a matrix.
* Input:
*   A: n-by-n matrix A
*   n: dimension of A
* Output:
*   A: mangled LU decomposition of A
*   Inv: inverse of A
*   perm: permutation vector

        subroutine Inverse(A, Inv, n, perm)
        implicit none
        integer n, perm(*)
        double complex A(n,*), Inv(n,*)

        integer i,j

        call CLUDecomp(A, n, perm)

        do i=1,n
          do j=1,n
            Inv(j,i) = 0
          enddo
          Inv(i,i) = 1.d0
          call CLUBackSubst(A, n, perm, Inv(1,i))
        enddo
        
        end

#define MAXDIM 9

************************************************************************
* LUDecomp computes the LU decomposition of the n-by-n matrix A
* by Gaussian Elimination with partial pivoting;
* compact (in situ) storage scheme
* Input:
*   A: n-by-n matrix to LU-decompose
*   n: dimension of A
* Output:
*   A: mangled LU decomposition of A in the form
*     ( y11 y12 ... y1m )
*     ( x21 y22 ... y2m )
*     ( x31 x32 ... y3m )
*     ( ............... )
*     ( xm1 xm2 ... ymm )
*   where 
*     (   1   0 ...   0 )  ( y11 y12 ... y1m )
*     ( x21   1 ...   0 )  (   0 y22 ... y2m )
*     ( x31 x32 ...   0 )  (   0   0 ... y3m )  =  Permutation(A)
*     ( ............... )  ( ............... )
*     ( xm1 xm2 ...   1 )  (   0   0 ... ymm ) 
*   perm: permutation vector

        subroutine LUDecomp(A, n, perm)
        implicit none
        integer n, perm(*)
        double precision A(n,*)

        integer i, j, k, imax
        double precision tmp
        double precision Amax

        do j = 1, n
* do U part (minus diagonal one)
          do i = 1, j - 1
            do k = 1, i - 1
              A(i,j) = A(i,j) - A(i,k)*A(k,j)
            enddo
          enddo

* do L part (plus diagonal from U case)
          Amax = 0
          do i = j, n
            do k = 1, j - 1
              A(i,j) = A(i,j) - A(i,k)*A(k,j)
            enddo

* do partial pivoting ...
* find the pivot
            if( abs(A(i,j)) .gt. Amax ) then
              Amax = abs(A(i,j))
              imax = i
            endif
          enddo

* exchange rows
          perm(j) = imax
          do k = 1, n
            tmp = A(j,k)
            A(j,k) = A(imax,k)
            A(imax,k) = tmp
          enddo

* division by the pivot element
          if( A(j,j) .eq. 0 ) then
            tmp = 1D123
          else
            tmp = 1/A(j,j)
          endif
          do i = j + 1, n
            A(i,j) = A(i,j)*tmp
          enddo
        enddo
        end

        subroutine CLUDecomp(A, n, perm)
        implicit none
        integer n, perm(*)
        double complex A(n,*)

        integer i, j, k, imax
        double complex tmp
        double precision Amax

        do j = 1, n
* do U part (minus diagonal one)
          do i = 1, j - 1
            do k = 1, i - 1
              A(i,j) = A(i,j) - A(i,k)*A(k,j)
            enddo
          enddo

* do L part (plus diagonal from U case)
          Amax = 0
          do i = j, n
            do k = 1, j - 1
              A(i,j) = A(i,j) - A(i,k)*A(k,j)
            enddo

* do partial pivoting ...
* find the pivot
            if( abs(A(i,j)) .gt. Amax ) then
              Amax = abs(A(i,j))
              imax = i
            endif
          enddo

* exchange rows
          perm(j) = imax
          do k = 1, n
            tmp = A(j,k)
            A(j,k) = A(imax,k)
            A(imax,k) = tmp
          enddo

* division by the pivot element
          if( A(j,j) .eq. 0 ) then
            tmp = 1D123
          else
            tmp = 1/A(j,j)
          endif
          do i = j + 1, n
            A(i,j) = A(i,j)*tmp
          enddo
        enddo
        end

************************************************************************
* LUBacksubst computes the x in A.x = b from the LU-decomposed A.
* Input:
*   A: LU-decomposed n-by-n matrix A
*   b: input vector b in A.x = b
*   n: dimension of A
*   p: Permutation vector from LU decomposition
* Output:
*   b: solution vector x in A.x = b

        subroutine LUBackSubst(A, n, p, b)
        implicit none

        integer n, p(*)
        double precision A(n,*)
        double complex b(*)

        integer i, j
        double precision tmp

* permute b 
        do i = 1, n
          tmp = b(i)
          b(i) = b(p(i))
          b(p(i)) = tmp
        enddo

* forward substitution L.Y = B
        do i = 1, n
          do j = 1, i - 1
            b(i) = b(i) - A(i,j)*b(j)
          enddo
        enddo

* backward substitution U.X = Y
        do i = n, 1, -1
          do j = i + 1, n
            b(i) = b(i) - A(i,j)*b(j)
          enddo
          b(i) = b(i)/A(i,i)
        enddo
        end


        subroutine CLUBackSubst(A, n, p, b)
        implicit none

        integer n, p(*)
        double complex A(n,*)
        double complex b(*)

        integer i, j
        double complex tmp

* permute b 
        do i = 1, n
          tmp = b(i)
          b(i) = b(p(i))
          b(p(i)) = tmp
        enddo

* forward substitution L.Y = B
        do i = 1, n
          do j = 1, i - 1
            b(i) = b(i) - A(i,j)*b(j)
          enddo
        enddo

* backward substitution U.X = Y
        do i = n, 1, -1
          do j = i + 1, n
            b(i) = b(i) - A(i,j)*b(j)
          enddo
          b(i) = b(i)/A(i,i)
        enddo
        end

************************************************************************
* Det computes the determinant of a matrix.
* Input:
*   A: n-by-n matrix A
*   n: dimension of A
* Output:
*   determinant of A
* Warning: A is overwritten

        double precision function Det(A, n) 

        implicit none
        integer n
        double precision A(n,*)

        integer i, perm(n)

        Det = 1
        call LUDecomp(A, n, perm)
        do i = 1, n
          Det = Det*A(i,i)
          if( perm(i) .ne. i ) Det = -Det
        enddo
        end

        double complex function CDet(A, n)
        implicit none
        integer n
        double complex A(n,*)

        integer i, perm(n)

        CDet = 1
        call CLUDecomp(A, n, perm)
        do i = 1, n
          CDet = CDet*A(i,i)
          if( perm(i) .ne. i ) CDet = -CDet
        enddo
        end

************************************************************************
** Empirical inverse matrix calculation for 2x2 and 3x3        

* Inverse matrix 2x2
** This subroutine calculates the inverse of a 2x2 matrix
      SUBROUTINE matinv2(A,B)
      
      IMPLICIT NONE
      
      DOUBLE COMPLEX A(2,2) !original matrix
      DOUBLE COMPLEX B(2,2) !inverse matrix
      DOUBLE COMPLEX detinv !determinant 
      
! Calculate the inverse determinant of the matrix
      detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))
! Calculate the inverse of the matrix
      B(1,1) = +detinv * A(2,2)
      B(2,1) = -detinv * A(2,1)
      B(1,2) = -detinv * A(1,2)
      B(2,2) = +detinv * A(1,1)      
      
      END
      
* Inverse matrix 3x3
** This subroutine calculates the inverse of a 3x3 matrix
      SUBROUTINE matinv3(A,B)
      
      IMPLICIT NONE
      
      DOUBLE COMPLEX A(3,3) !original matrix
      DOUBLE COMPLEX B(3,3) !inverse matrix
      DOUBLE COMPLEX detinv !determinant 
      
! Calculate the inverse determinant of the matrix
      detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)
     &          - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)
     &         + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

! Calculate the inverse of the matrix
      B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
      B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
      B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
      B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
      B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
      B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
      B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
      B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
      B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))      
      
      END 
