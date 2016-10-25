
      subroutine ZZZPentcheck(psi, p, z1, qz1, z2, qz2, z3, qz3, 
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
      double complex j5z1z2(0:3,-1:1,3)
      double precision qz1(0:4), qz2(0:4), qz3(0:4)
      double complex z1(6), z2(6), z3(6)

c  local variables
      integer i, mu, k, isig
      double complex m1(3,-1:1), m2(3,-1:1), m3(3,-1:1)
      double complex mv1(3,-1:1), mv2(3,-1:1), mv3(3,-1:1)
      double complex mvirt(3,-1:1)
      double complex mvcheck(3,-1:1)
      
      double complex psiz1(2,-1:1,2), psiz2(2,-1:1,2)
      double precision fqz1(0:4,2), fqz2(0:4,2)
      double complex jz1(0:5,-1:1,2), jz2(0:5,-1:1,2)
      double complex jvz1(0:5,-1:1,2), jvz2(0:5,-1:1,2)

      integer icount, bcount
      data bcount /0/
      data icount /0/
      save icount, bcount

      double complex zero
      parameter (zero=(0d0,0d0))

      logical ldebug
      parameter (ldebug = .false.)

c  external functions
      double complex s1c, dotcc, dotrc
      external s1c, dotcc, dotrc

c--------------------------- Begin Code ----------------------------

c calculate bras and kets
      do isig= -1,1,2
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qz1(0),z1(1),
     &              psiz1(1,isig,1),fqz1(0,1)) 
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qz2(0),z2(1),
     &              psiz2(1,isig,2),fqz2(0,2))
      enddo

c different currents
      call curr6(1,psiz2(1,-1,2),fqz2(0,2),psi(1,-1,1),p(0,1),jz2(0,-1,2))
      call curr6(1,psi(1,-1,2),p(0,2),psiz1(1,-1,1),fqz1(0,1),jz1(0,-1,1))

      call TRIboxlinec(1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
     1     .true., z1(1), 2, jz1(0,-1,1), jvz1(0,-1,1) )
      call TRIboxlinec(1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
     1     .true., z2(1), 1, jz2(0,-1,2), jvz2(0,-1,2) )

c q qbar -> VVV, z3 -> q(Z_3)
      do isig = -1,1,2
         m1(3,isig) = -dotcc(z2(1),jz1(0,isig,1))    !Z1 Z2 Z3
         m2(3,isig) = -dotcc(z1(1),jz2(0,isig,2))    !Z1 Z3 Z2
     1                +dotcc(z2(1),jz1(0,isig,1))
         m3(3,isig) = +dotcc(z1(1),jz2(0,isig,2))    !Z3 Z1 Z2

c Virtual diagrams only consist of boxes * cvirt Mborn
         mv1(3,isig) = -dotcc(z2(1),jvz1(0,isig,1))
         mv2(3,isig) = +dotcc(z2(1),jvz1(0,isig,1))
     1                 -dotcc(z1(1),jvz2(0,isig,2)) 
         mv3(3,isig) = +dotcc(z1(1),jvz2(0,isig,2)) 

         mvcheck(1,isig) = + mv3(3,isig)
         mvcheck(2,isig) = + mv2(3,isig)
         mvcheck(3,isig) = + mv1(3,isig)

         do i = 1,3
            mvirt(i,isig) = -dotrc(qz3(0),j5z1z2(0,isig,i))
         enddo
      
         if( abs(mvirt(1,isig)/mvcheck(1,isig)-1d0).gt.1d-1 .or.
     1       abs(mvcheck(1,isig)/mvirt(1,isig)-1d0).gt.1d-1 .or.
     2       abs(mvirt(2,isig)/mvcheck(2,isig)-1d0).gt.1d-1 .or.
     3       abs(mvcheck(2,isig)/mvirt(2,isig)-1d0).gt.1d-1 .or.
     4       abs(mvirt(3,isig)/mvcheck(3,isig)-1d0).gt.1d-1 .or.
     5       abs(mvcheck(3,isig)/mvirt(3,isig)-1d0).gt.1d-1 ) then
            icount = icount+1
            do mu = 0,3
               j5z1z2(mu,isig,1) = (0.0d0,0.0d0)
               j5z1z2(mu,isig,2) = (0.0d0,0.0d0)
               j5z1z2(mu,isig,3) = (0.0d0,0.0d0)
            enddo
         else if ( mvcheck(1,isig).eq.0 .or.
     1             mvcheck(2,isig).eq.0 .or.
     2             mvcheck(3,isig).eq.0 ) then
            icount = icount+1
            do mu = 0,3
               j5z1z2(mu,isig,1) = (0.0d0,0.0d0)
               j5z1z2(mu,isig,2) = (0.0d0,0.0d0)
               j5z1z2(mu,isig,3) = (0.0d0,0.0d0)
            enddo
         else
            bcount = bcount+1
         endif      

      enddo !isig

      do i = 10,10000000,1000
      if (bcount.eq.1d2*i) print*,"ratio= ",icount*1.0d0/(bcount+icount)
      enddo

      return
      end
