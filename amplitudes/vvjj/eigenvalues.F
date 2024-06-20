ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Calculate eigenvalues of a complex n x n matrix
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Reference: 
c Folkmar Bornemann, course notes to Numerische Mathematik 1&2, 
c Technische Universitaet Muenchen, Germany
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Author: Michael Rauch, last changed: Apr 24, 2013
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine hessenberg(A, m, ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c turns A into Hessenberg form via Householder reflections
c Input:  A(1..m,1..m): complex m x m matrix 
c         m: dimension of A
c Output: A(1..m,1..m): transformed A
c         ier: error flag
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       implicit none

       integer m, ier
       double complex A(m,m)

       integer i,j,k, ip
       double complex rho, vecv(m), vTv, vecvTA(m)
       double precision absa

       ier = 0

       do i=1,m-1
         ip=i+1
c get Householder vector         
c -- leave first line untouched
         absa=0d0
         do j=ip,m
           absa = absa + abs(A(j,i))**2
         enddo
         if (A(ip,i) .ne. 0d0) then
           rho = - A(ip,i)/abs(A(ip,i)) * sqrt(absa)
         else
           rho = - sqrt(absa)
         endif
         if (rho.ne.0d0) then
c otherwise all entries are zero 
c -> has already required form
c -> next index
           do j=ip,m
             vecv(j) = -A(j,i)/rho
           enddo
           vecv(ip) = vecv(ip) + 1d0
           if (dreal(vecv(ip)).eq.0d0) then
c shouldn't happen (?)
             ier = 1
             return
           endif
c substitute it in
           do k=i,m
             vecvTA(k) = 0d0
             do j=ip,m
               vecvTA(k) = vecvTA(k) + dconjg(vecv(j))*A(j,k)
             enddo
             do j=ip,m
               A(j,k) = A(j,k) - vecv(j)/dreal(vecv(ip))*vecvTA(k)
             enddo
           enddo
c and now H^dagger from right to preserve unitarity
           do j=1,m
             vecvTA(j) = 0d0
             do k=ip,m
               vecvTA(j) = vecvTA(j) + A(j,k)*vecv(k)
             enddo
             do k=ip,m
               A(j,k) = A(j,k) 
     &                  - vecvTA(j)/dreal(vecv(ip))*dconjg(vecv(k))

             enddo
           enddo
         endif
       enddo

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine eigenvaluestep(A, m, mdo)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computes a single step via QR decomposition using Givens rotation
c Input:  A(1..m,1..m): complex m x m matrix 
c         m: dimension of A
c         mdo: A(1..mdo,1..mdo) submatrix to perform step on
c Output: A(1..m,1..m): transformed A
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       

       implicit none

       integer m, mdo
       double complex A(m,m)

       integer i,j,k, ip
       double complex cg(mdo), sg(mdo), tg, tmpAij
       double complex wshift, ddiff, sdiff, wnorm
c Wilkinson shift
c from eigenvalues of lower-left 2x2 matrix
       ddiff = (A(mdo-1,mdo-1)-A(mdo,mdo))/2d0
       if (ddiff.eq.0d0) then
         sdiff = 1d0
       else
         sdiff = ddiff/abs(ddiff)
       endif
       wshift = A(mdo,mdo) 
       wnorm = (abs(ddiff)+sqrt(ddiff**2+A(mdo-1,mdo)*A(mdo,mdo-1)))
       if (wnorm .ne. 0d0) then
         wshift = wshift - sdiff*A(mdo-1,mdo)*A(mdo,mdo-1)/wnorm
       endif
       do i=1,mdo
         A(i,i) = A(i,i) - wshift
       enddo

c make tri-diagonal by Givens rotation
       do i=1,mdo-1
         if (A(i+1,i) .eq. 0d0) then
c already zero -- nothing to do
           cg(i) = 1d0
           sg(i) = 0d0
         else
           if ( abs(A(i+1,i)) .ge. abs(A(i,i)) ) then
             tg = - A(i,i)/A(i+1,i)
             sg(i) = 1d0/sqrt(1d0+tg**2)
             cg(i) = sg(i)*tg
           else
             tg = - A(i+1,i)/A(i,i)
             cg(i) = 1d0/sqrt(1d0+tg**2)
             sg(i) = cg(i)*tg
           endif
c rotate
           do j=i,mdo
             tmpAij   = cg(i)*A(i,j)   - sg(i)*A(i+1,j) 
             A(i+1,j) = cg(i)*A(i+1,j) + sg(i)*A(i,j) 
             A(i  ,j) = tmpAij
           enddo
         endif
       enddo
c now A has become R and Q is stored in cg,sg
c Calculate R.Q
       do i=1,mdo-1
         do j=1,i+1
           tmpAij   = cg(i)*A(j,i)   - sg(i)*A(j,i+1) 
           A(j,i+1) = cg(i)*A(j,i+1) + sg(i)*A(j,i) 
           A(j,i  ) = tmpAij
         enddo
       enddo

c undo Wilkinson shift
       do i=1,mdo
         A(i,i) = A(i,i) + wshift
       enddo

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine eigenvalues(A, m, acc, maxiter, ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computes eigenvalues via QR decomposition
c Input:  A(1..m,1..m): complex m x m matrix 
c         m: dimension of A
c         acc: requested accuracy 
c              defined as sum|(i+1,i)-off-diag|/sum|diag|
c Output: A(1..m,1..m): transformed A with eigenvalues on diagonal
c         acc: achieved accuracy
c         ier: error flag
c              ier>0: internal errors 
c                     (4: maximal number of iterations reached)
c              ier<0: - number of iterations used
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       

       implicit none

       integer m, maxiter, ier
       double complex A(m,m)
       double precision acc

       integer i, iter, mloc
       double precision myacc, sumoff, sumdiag
       double complex acc_nom, acc_denom
       logical errorinunitarization
c transform A into Hessenberg form by successive Householder reflections
       errorinunitarization =.false.
       sumdiag = 0d0
       do i=1,m
         sumdiag = sumdiag + abs(A(i,i))
       enddo
!       if(isnan(dble(sumdiag)))then
!       print*, 'sumdiag is nan'
!       endif
       if(sumdiag.lt.acc/m/m)then
        return
       endif

       call hessenberg(A,m,ier)
c check for errors
       if (ier.ne.0) then
         return
       endif

       mloc = m

       do iter=1,maxiter
         ier = -iter
c call a step
         call eigenvaluestep(A, m, mloc)


c check whether accuracy is ok
c compare last off-diagonal element to last two diagonal ones
         myacc = 0d0
         do while ( (mloc.gt.1) .and. (myacc.lt.acc) )
          if(abs(A(mloc,mloc-1)).lt.sumdiag*acc/m/m) then
            A(mloc,mloc-1) =0d0
            A(mloc-1,mloc) =0d0
          endif
           acc_nom   = 2*A(mloc,mloc-1)
           acc_denom = (abs(A(mloc-1,mloc-1))+abs(A(mloc,mloc)))
           if (abs(acc_nom).eq.0d0 .and. abs(acc_denom).eq.0d0) then
             myacc = 0d0
           elseif (abs(acc_nom).ne.0d0 .and. abs(acc_denom).eq.0d0) then
             myacc = acc * 1000d0
           else
             myacc = abs(acc_nom/acc_denom)
!             if(isnan(dble(myacc)))then
!              print*, 'myacc is nan'
!             endif
           endif
           if (myacc.lt.acc) mloc = mloc-1
         enddo


c got it
         if (mloc.le.1) then
          sumoff = 0d0
          do i=1,m-1
            sumoff = sumoff + abs(A(i+1,i))
          enddo
           acc = sumoff/sumdiag

           return 
         endif

       enddo

c could not achieve accuracy
       sumoff = 0d0
       do i=1,m-1
         sumoff = sumoff + abs(A(i+1,i))
       enddo
       acc = sumoff/sumdiag
       ier = 4
       if(ier.eq.4) then
        print*, 'could not achieve accuracy'
        print*, acc
        print*, 'after', iter ,'iterations'
        print*, 'with size', mloc
        print*, sumoff
        print*, sumdiag
        errorinunitarization =.true.
       endif

       return
       end


