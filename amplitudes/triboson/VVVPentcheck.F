
      subroutine VVVPentcheck(psi, p, z1, qz1, z2, qz2, z3, qz3, 
     1                                      j5z1z2)

c****************************************************************
c
c this subroutine computes the born and virtual matrix elements with 
c Z3, A3 replaced by q_3 and compares it with the analytic result
c which consists of a sum of diagrams with 2 vector bosons
c in the final state.
c
c***************************************************************

      implicit none

c  input / output

      double complex psi(2,-1:1,2)
      double precision p(0:3,8)
      double complex j5z1z2(0:3,3)
      double precision qz1(0:4), qz2(0:4), qz3(0:4)
      double complex z1(6), z2(6), z3(6)
      integer number

c  local variables

      integer i, mu, k, isig
      double complex m1(3), m2(3), m3(3)
      double complex mv1(3), mv2(3), mv3(3)
      double complex mvirt(3)
      double complex mvcheck(3)
      
      double complex psiz1(2,-1:1,2), psiz2(2,-1:1,2)

      double precision fqz1(0:4,2), fqz2(0:4,2)

      double complex jz1(0:5,-1:1,2), jz2(0:5,-1:1,2)
      double complex jvz1(0:5,-1:1,2), jvz2(0:5,-1:1,2)

      integer icount, bcount, lbox
      data bcount /0/
      data icount /0/
      data lbox /0/
      save icount, bcount, lbox

      double complex zero
      parameter (zero=(0d0,0d0))

      logical ldebug
      parameter (ldebug = .false.)

c  external functions

      double complex s1c, dotcc, dotrc
      external s1c, dotcc, dotrc

c--------------------------- Begin Code ----------------------------

c calculate bras and kets
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qz1(0),z1(1),
     &           psiz1(1,-1,1),fqz1(0,1)) 
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qz2(0),z2(1),
     &           psiz2(1,-1,2),fqz2(0,2))

c different currents
      call curr6(-1,psiz2(1,-1,2),fqz2(0,2),psi(1,-1,1),p(0,1),jz2(0,-1,2))
      call curr6(-1,psi(1,-1,2),p(0,2),psiz1(1,-1,1),fqz1(0,1),jz1(0,-1,1))

      call TRIboxlinec(-1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
     1     .true., z1(1), 2, jz1(0,-1,1), jvz1(0,-1,1) )
      call TRIboxlinec(-1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
     1     .true., z2(1), 1, jz2(0,-1,2), jvz2(0,-1,2) )

c q qbar -> VVV, z3 -> q(Z_3)
      m1(3) = -dotcc(z2(1),jz1(0,-1,1))    !Z1 Z2 Z3
      m2(3) = -dotcc(z1(1),jz2(0,-1,2))    !Z1 Z3 Z2
     1                +dotcc(z2(1),jz1(0,-1,1))
      m3(3) = +dotcc(z1(1),jz2(0,-1,2))    !Z3 Z1 Z2

c Virtual diagrams only consist of boxes * cvirt Mborn
      mv1(3) = -dotcc(z2(1),jvz1(0,-1,1))
      mv2(3) = +dotcc(z2(1),jvz1(0,-1,1))
     1                 -dotcc(z1(1),jvz2(0,-1,2)) 
      mv3(3) = +dotcc(z1(1),jvz2(0,-1,2)) 

      mvcheck(1) = + mv3(3)
      mvcheck(2) = + mv2(3)
      mvcheck(3) = + mv1(3)

      do i = 1,3
         mvirt(i) = -dotrc(qz3(0),j5z1z2(0,i))
      enddo

c      if (abs(mvcheck(1)/mvirt(1)-1d0).gt.1d-4) then
c      write(*,*)'1', abs(mvcheck(1)/mvirt(1)-1d0)
c      elseif (abs(mvcheck(2)/mvirt(2)-1d0).gt.1d-4) then
c      write(*,*)'2', abs(mvcheck(2)/mvirt(2)-1d0)
c      elseif (abs(mvcheck(3)/mvirt(3)-1d0).gt.1d-4) then
c      write(*,*)'3', abs(mvcheck(3)/mvirt(3)-1d0)
c      endif

      if( abs(mvcheck(1)/mvirt(1)-1d0).gt.1d-1 ) then
         if (mvcheck(1).eq.0) lbox = lbox+1
         icount = icount+1
         do mu = 0,3
            j5z1z2(mu,1) = (0.0d0,0.0d0)
            j5z1z2(mu,2) = (0.0d0,0.0d0)
            j5z1z2(mu,3) = (0.0d0,0.0d0)
         enddo
      else if( abs(mvcheck(2)/mvirt(2)-1d0).gt.1d-1 ) then
         if (mvcheck(2).eq.0) lbox = lbox+1
         icount = icount+1
         do mu = 0,3
            j5z1z2(mu,1) = (0.0d0,0.0d0)
            j5z1z2(mu,2) = (0.0d0,0.0d0)
            j5z1z2(mu,3) = (0.0d0,0.0d0)
         enddo
      else if( abs(mvcheck(3)/mvirt(3)-1d0).gt.1d-1 ) then
         if (mvcheck(3).eq.0) lbox = lbox+1
         icount = icount+1
         do mu = 0,3
            j5z1z2(mu,1) = (0.0d0,0.0d0)
            j5z1z2(mu,2) = (0.0d0,0.0d0)
            j5z1z2(mu,3) = (0.0d0,0.0d0)
         enddo
      else
         bcount = bcount+1
      endif

      do i = 10,10000000,1000
      if (bcount.eq.1d2*i) print*,'bcount,icount,lbox = ',bcount,icount,lbox
      enddo

      return
      end
