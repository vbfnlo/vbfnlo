      subroutine qqVVgg(p,v,selc,m2s)

      implicit none

      include 'VBFNLO/utilities/global.inc'
      include 'VBFNLO/utilities/lha.inc'
      include 'VBFNLO/utilities/scales.inc'
      include 'VBFNLO/utilities/coupl.inc'
      include 'tensor.inc'

      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)

      integer j,i,k,l,mu,isig,l1,l2
      integer sign(4),polmap(2)

      double precision p(0:3,max_p), v(0:3,max_v), m2s(2)
      double precision pbar(0:3,4), qbar(0:3,4),
     &       gluvecb(0:4,2), gluvec(0:4,2), glumom(0:4)
      double precision pmomv(0:4,2,5),pmomvv(0:4,2,5,5),
     &       pmomvvv(0:4,2,5,5,5)

      double complex mat(6,-1:1,3), qqgg(2)
      double complex psi(2,-1:1,4),
     &       psiv(2,-1:1,2,5),psivv(2,-1:1,2,5,5),
     &       psivvv(2,-1:1,2,5,5,5)
      double complex effgluhel(0:3)
      double complex jv(1:6,-1:1,5:5,2), jvv(1:6,-1:1,5,5,3), 
     &       jvvv(1:6,-1:1,5,5,5,4)

      logical ldebug
      parameter (ldebug=.false.)
    
      double complex dotcc, dotrc
      double precision dotrr
      external dotrr, dotrc, dotcc

      DOUBLE PRECISION fpials(2:3), fpi
      PARAMETER (fpi=4d0*pi)

      double complex temp,temp1
      integer selc,cs

      double precision pall(0:3,8), mg2uures, mg2ddres
      integer fsign(8)

      do i=1,2
        m2s(i)=0d0
      enddo

c fix strong coupling gs**2 for the two quarks:

      fpials(2) = fpi*als(1,1)
      fpials(3) = fpi*als(2,1)
c      fpials(2) = fpi*alfas

!
! Get information about different subprocesses and calculate external quark spinors

      Select Case (selc)

!-----------------------------------------------------
      Case (1) 
!  q  q~ -> e+ ve mu- vm~ g g  

      sign(1)=1
      sign(2)=-1
      polmap(1)=1
      polmap(2)=1

      do mu=0,3
       pbar(mu,1)=p(mu,1)
       pbar(mu,2)=p(mu,2)
       gluvecb(mu,1)=p(mu,3)
       gluvecb(mu,2)=p(mu,4)
       gluvec(mu,1)=polmap(1)*gluvecb(mu,1)
       gluvec(mu,2)=polmap(2)*gluvecb(mu,2)
       glumom(mu)=gluvec(mu,1)+gluvec(mu,2)
       do i=1,2
         qbar(mu,i)=sign(i)*pbar(mu,i)
       enddo
      enddo

!-----------------------------------------------------

      Case (2) 
!  q  g -> e+ ve mu- vm~ q  g  

      sign(1)=1
      sign(2)=1
      polmap(1)=-1
      polmap(2)=1

      do mu=0,3
       pbar(mu,1)=p(mu,1)
       pbar(mu,2)=p(mu,3)
       gluvecb(mu,1)=p(mu,2)
       gluvecb(mu,2)=p(mu,4)
       gluvec(mu,1)=polmap(1)*gluvecb(mu,1)
       gluvec(mu,2)=polmap(2)*gluvecb(mu,2)
       glumom(mu)=gluvec(mu,1)+gluvec(mu,2)
       do i=1,2
         qbar(mu,i)=sign(i)*pbar(mu,i)
       enddo
      enddo

!-----------------------------------------------------

      Case (3) 
!  q~ g -> e+ ve mu- vm~ q~ g

      sign(1)=-1
      sign(2)=-1
      polmap(1)=-1
      polmap(2)=1

      do mu=0,3
       pbar(mu,1)=p(mu,3)
       pbar(mu,2)=p(mu,1)
       gluvecb(mu,1)=p(mu,2)
       gluvecb(mu,2)=p(mu,4)
       gluvec(mu,1)=polmap(1)*gluvecb(mu,1)
       gluvec(mu,2)=polmap(2)*gluvecb(mu,2)
       glumom(mu)=gluvec(mu,1)+gluvec(mu,2)
       do i=1,2
         qbar(mu,i)=sign(i)*pbar(mu,i)
       enddo
      enddo

!-----------------------------------------------------

      Case (4) 
!  g g -> e+ ve mu- vm~ q q~  

      sign(1)=-1
      sign(2)=1
      polmap(1)=-1
      polmap(2)=-1

      do mu=0,3
       pbar(mu,1)=p(mu,4)
       pbar(mu,2)=p(mu,3)
       gluvecb(mu,1)=p(mu,1)
       gluvecb(mu,2)=p(mu,2)
       gluvec(mu,1)=polmap(1)*gluvecb(mu,1)
       gluvec(mu,2)=polmap(2)*gluvecb(mu,2)
       glumom(mu)=gluvec(mu,1)+gluvec(mu,2)
       do i=1,2
         qbar(mu,i)=sign(i)*pbar(mu,i)
       enddo
      enddo

      end select

      gluvec(4,1)=0d0
      gluvec(4,2)=0d0
      glumom(4) = glumom(0)**2-glumom(1)**2-glumom(2)**2-glumom(3)**2

#ifdef WITH_MADGRAPH
!
! Write settings in variables for MadGraph
!
      if (ldebug) then
        fsign(1)=sign(1)
        fsign(2)=sign(2)
        fsign(7)=polmap(1)
        fsign(8)=polmap(2)
        fsign(3)=-1
        fsign(4)=+1
        fsign(5)=+1
        fsign(6)=-1
        do mu=0,3
          pall(mu,1) = pbar(mu,1)
          pall(mu,2) = pbar(mu,2)
          pall(mu,7) = gluvecb(mu,1)
          pall(mu,8) = gluvecb(mu,2)
          pall(mu,3) = v(mu,2)
          pall(mu,4) = v(mu,1)
          pall(mu,5) = v(mu,3)
          pall(mu,6) = v(mu,4)
        enddo
        call SUU_WWGG(pall,fsign,mg2uures)
        call SDD_WWGG(pall,fsign,mg2ddres)
      endif
#endif

!-----------------------------------------------------
!-----------------------------------------------------

      call psi0m(2,pbar(0,1),sign(1),psi)

      do l1=-1,1,2       ! Loop over gluon (1) helicity
       do l2 = -1,1,2    ! Loop over gluon (2) helicity

        do j=1,3
        do i=1,6
          do isig=-1,1,2
            mat(i,isig,j)=0d0
          enddo
        enddo
        enddo

        call helvec(gluvecb(0,1),polmap(1),l1,glupolc(0,1,l1))
        call helvec(gluvecb(0,2),polmap(2),l2,glupolc(0,2,l2))

!
!  Add from one to three gauge bosons to the bra
!
!  in psig() etc: g1=1, g2=2, w+=3, w-=4
!  order in arguments is as in math. expression <bra|V...   ...V|ket>
!

       do isig=-1,1,2
! first V to <bra|
        call bra2c(psi(1,isig,2),.true.,qbar(0,2),isig,
     &    gluvec(0,1),glupolc(0,1,l1),psiv(1,isig,2,1),pmomv(0,2,1))
        call bra2c(psi(1,isig,2),.true.,qbar(0,2),isig,
     &    gluvec(0,2),glupolc(0,2,l2),psiv(1,isig,2,2),pmomv(0,2,2))
        call bra2c(psi(1,isig,2),.true.,qbar(0,2),isig,
     &    qp(0,1),wp(1,1),psiv(1,isig,2,3),pmomv(0,2,3))
        call bra2c(psi(1,isig,2),.true.,qbar(0,2),isig,
     &    qm(0,1),wm(1,1),psiv(1,isig,2,4),pmomv(0,2,4))

! second V to <bra|V
        do i=1,4
          if (i.ne.1)
     &      call bra2c(psiv(1,isig,2,i),.false.,pmomv(0,2,i),isig,
     &    gluvec(0,1),glupolc(0,1,l1),psivv(1,isig,2,i,1),pmomvv(0,2,i,1))
          if (i.ne.2)
     &      call bra2c(psiv(1,isig,2,i),.false.,pmomv(0,2,i),isig,
     &    gluvec(0,2),glupolc(0,2,l2),psivv(1,isig,2,i,2),pmomvv(0,2,i,2))
          if (i.ne.3)
     &      call bra2c(psiv(1,isig,2,i),.false.,pmomv(0,2,i),isig,
     &    qp(0,1),wp(1,1),psivv(1,isig,2,i,3),pmomvv(0,2,i,3))
          if (i.ne.4)
     &      call bra2c(psiv(1,isig,2,i),.false.,pmomv(0,2,i),isig,
     &    qm(0,1),wm(1,1),psivv(1,isig,2,i,4),pmomvv(0,2,i,4))
        enddo

! third V to <bra|VV
        do i=1,4
        do j=1,4
          if ((i.ne.1).and.(j.ne.1))
     &      call bra2c(psivv(1,isig,2,i,j),.false.,pmomvv(0,2,i,j),isig,
     &    gluvec(0,1),glupolc(0,1,l1),psivvv(1,isig,2,i,j,1),pmomvvv(0,2,i,j,1))
          if ((i.ne.2).and.(j.ne.2))
     &      call bra2c(psivv(1,isig,2,i,j),.false.,pmomvv(0,2,i,j),isig,
     &    gluvec(0,2),glupolc(0,2,l2),psivvv(1,isig,2,i,j,2),pmomvvv(0,2,i,j,2))
          if ((i.ne.3).and.(j.ne.3))
     &      call bra2c(psivv(1,isig,2,i,j),.false.,pmomvv(0,2,i,j),isig,
     &    qp(0,1),wp(1,1),psivvv(1,isig,2,i,j,3),pmomvvv(0,2,i,j,3))
          if ((i.ne.4).and.(j.ne.4))
     &      call bra2c(psivv(1,isig,2,i,j),.false.,pmomvv(0,2,i,j),isig,
     &    qm(0,1),wm(1,1),psivvv(1,isig,2,i,j,4),pmomvvv(0,2,i,j,4))
        enddo
        enddo

!
!  Add from one to four gauge bosons to the ket
!

! first V to |ket>
        call ket2c(psi(1,isig,1),.true.,qbar(0,1),isig,
     &    gluvec(0,1),glupolc(0,1,l1),psiv(1,isig,1,1),pmomv(0,1,1))
        call ket2c(psi(1,isig,1),.true.,qbar(0,1),isig,
     &    gluvec(0,2),glupolc(0,2,l2),psiv(1,isig,1,2),pmomv(0,1,2))
        call ket2c(psi(1,isig,1),.true.,qbar(0,1),isig,
     &    qp(0,1),wp(1,1),psiv(1,isig,1,3),pmomv(0,1,3))
        call ket2c(psi(1,isig,1),.true.,qbar(0,1),isig,
     &    qm(0,1),wm(1,1),psiv(1,isig,1,4),pmomv(0,1,4))

! second V to V|ket>
        do i=1,4
          if (i.ne.1)
     &      call ket2c(psiv(1,isig,1,i),.false.,pmomv(0,1,i),isig,
     &    gluvec(0,1),glupolc(0,1,l1),psivv(1,isig,1,1,i),pmomvv(0,1,1,i))
          if (i.ne.2)
     &      call ket2c(psiv(1,isig,1,i),.false.,pmomv(0,1,i),isig,
     &    gluvec(0,2),glupolc(0,2,l2),psivv(1,isig,1,2,i),pmomvv(0,1,2,i))
          if (i.ne.3)
     &      call ket2c(psiv(1,isig,1,i),.false.,pmomv(0,1,i),isig,
     &    qp(0,1),wp(1,1),psivv(1,isig,1,3,i),pmomvv(0,1,3,i))
          if (i.ne.4)
     &      call ket2c(psiv(1,isig,1,i),.false.,pmomv(0,1,i),isig,
     &    qm(0,1),wm(1,1),psivv(1,isig,1,4,i),pmomvv(0,1,4,i))
        enddo

! third V to VV|ket>
        do i=1,4
        do j=1,4
          if ((i.ne.1).and.(j.ne.1))
     &      call ket2c(psivv(1,isig,1,i,j),.false.,pmomvv(0,1,i,j),isig,
     &    gluvec(0,1),glupolc(0,1,l1),psivvv(1,isig,1,1,i,j),pmomvvv(0,1,1,i,j))
          if ((i.ne.2).and.(j.ne.2))
     &      call ket2c(psivv(1,isig,1,i,j),.false.,pmomvv(0,1,i,j),isig,
     &    gluvec(0,2),glupolc(0,2,l2),psivvv(1,isig,1,2,i,j),pmomvvv(0,1,2,i,j))
          if ((i.ne.3).and.(j.ne.3))
     &      call ket2c(psivv(1,isig,1,i,j),.false.,pmomvv(0,1,i,j),isig,
     &    qp(0,1),wp(1,1),psivvv(1,isig,1,3,i,j),pmomvvv(0,1,3,i,j))
          if ((i.ne.4).and.(j.ne.4))
     &      call ket2c(psivv(1,isig,1,i,j),.false.,pmomvv(0,1,i,j),isig,
     &    qm(0,1),wm(1,1),psivvv(1,isig,1,4,i,j),pmomvvv(0,1,4,i,j))
        enddo
        enddo

      enddo ! isig

!
! Now cycle through position of g1 and build currents
! always need both gluons attached
! 
      do i=2,4
      do j=2,4
        if (i.eq.j) cycle
      do k=2,4
        if ((i.eq.k).or.(j.eq.k)) cycle
! all 4 bosons must be attached 
        call curr6(-1,psivvv(1,-1,2,i,j,k),pmomvvv(0,2,i,j,k),psi(1,-1,1),qbar(0,1),jvvv(1,-1,i,j,k,4))
        call curr6(-1,psivv(1,-1,2,i,j),pmomvv(0,2,i,j),psiv(1,-1,1,k),pmomv(0,1,k),jvvv(1,-1,i,j,k,3))
        call curr6(-1,psiv(1,-1,2,i),pmomv(0,2,i),psivv(1,-1,1,j,k),pmomvv(0,1,j,k),jvvv(1,-1,i,j,k,2))
        call curr6(-1,psi(1,-1,2),qbar(0,2),psivvv(1,-1,1,i,j,k),pmomvvv(0,1,i,j,k),jvvv(1,-1,i,j,k,1))
      enddo
      enddo
      enddo

!
! Now neutral currents attached to fermion line
!
      do i=1,2
        j=3-i
        call curr6(+1,psivv(1,-1,2,i,j),pmomvv(0,2,i,j),psi(1,-1,1),qbar(0,1),jvv(1,-1,i,j,3))
        call curr6(+1,psiv(1,-1,2,i),pmomv(0,2,i),psiv(1,-1,1,j),pmomv(0,1,j),jvv(1,-1,i,j,2))
        call curr6(+1,psi(1,-1,2),qbar(0,2),psivv(1,-1,1,i,j),pmomvv(0,1,i,j),jvv(1,-1,i,j,1))
      enddo

!
!  Contract currents with the leptonic tensors
!

      do i=2,4
      do j=2,4
        if (i.eq.j) cycle
      do k=2,4
        if ((i.eq.k).or.(j.eq.k)) cycle
      do l=1,4
! check for order of g's for colour structure
        if (((i.eq.2).and.(l.le.1)).or. 
     &      ((j.eq.2).and.(l.le.2)).or. 
     &      ((k.eq.2).and.(l.le.3))) then
          cs = 1  ! 1-2
        else
          cs = 2  ! 2-1
        endif
! check for order of W's
! for u-quarks only 4-3 is allowed, for d only 3-4
        if (((k.eq.3).and.((j.eq.4).or.(i.eq.4)))
     &      .or.((j.eq.3).and.(i.eq.4)))
     &    mat(cs,-1,3) = mat(cs,-1,3) 
     &             + dotcc(jvvv(1,-1,i,j,k,l),glupolc(0,1,l1))
        if (((k.eq.4).and.((j.eq.3).or.(i.eq.3)))
     &      .or.((j.eq.4).and.(i.eq.3)))
     &    mat(3+cs,-1,3) = mat(3+cs,-1,3) 
     &             + dotcc(jvvv(1,-1,i,j,k,l),glupolc(0,1,l1))
      enddo
      enddo
      enddo
      enddo

      do isig=-1,1,2
      do i=1,2
         j=3-i
      do k=1,3
! check for order of g's for colour structure
        if ((i.eq.1).and.(j.eq.2)) then
          cs = 1  ! 1-2
        else
          cs = 2  ! 2-1
        endif
        mat(cs,isig,2) = mat(cs,isig,2) 
     &             - dotcc(jvv(1,isig,i,j,k),zww(0,1))
        mat(cs,isig,1) = mat(cs,isig,1) 
     &             - dotcc(jvv(1,isig,i,j,k),aww(0,1))
      enddo
      enddo
      enddo
      

! 
!   Get effective g(1)g(2) -> g polarization
!
        do mu=0,3
        effgluhel(mu) = (2d0*dotrc(gluvec(0,1),glupolc(0,2,l2))*glupolc(mu,1,l1)
     &       -2d0*dotrc(gluvec(0,2),glupolc(0,1,l1))*glupolc(mu,2,l2)
     &       +dotcc(glupolc(0,1,l1),glupolc(0,2,l2))*(gluvec(mu,2)-gluvec(mu,1)))
     &       /(glumom(4))
        enddo

!
! Now add again the vector bosons to the <bra| and |ket>
!
      do isig=-1,1,2
! first V to <bra|
        call bra2c(psi(1,isig,2),.true.,qbar(0,2),isig,
     &    glumom(0),effgluhel(0),psiv(1,isig,2,5),pmomv(0,2,5)) !5

! first V to |ket>
        call ket2c(psi(1,isig,1),.true.,qbar(0,1),isig,
     &    glumom(0),effgluhel(0),psiv(1,isig,1,5),pmomv(0,1,5)) !5
       enddo ! isig


!
! Now cycle through position of g(5) and build currents
! 
      do i=3,4
        j=7-i
        call curr6(-1,psivv(1,-1,2,i,j),pmomvv(0,2,i,j),psi(1,-1,1),qbar(0,1),jvv(1,-1,i,j,3))
        call curr6(-1,psiv(1,-1,2,i),pmomv(0,2,i),psiv(1,-1,1,j),pmomv(0,1,j),jvv(1,-1,i,j,2))
        call curr6(-1,psi(1,-1,2),qbar(0,2),psivv(1,-1,1,i,j),pmomvv(0,1,i,j),jvv(1,-1,i,j,1))
      enddo

      do i=3,4
        j=7-i
! check for order of W's
! for u-quarks only 4-3 is allowed, for d only 3-4
        if ((j.eq.3).and.(i.eq.4)) then
          do k=1,3
            mat(3,-1,3) = mat(3,-1,3) 
     &                + dotcc(jvv(1,-1,i,j,k),effgluhel(0))
          enddo
        elseif ((j.eq.4).and.(i.eq.3)) then
          do k=1,3
            mat(6,-1,3) = mat(6,-1,3) 
     &                + dotcc(jvv(1,-1,i,j,k),effgluhel(0))
          enddo
        endif
      enddo

      call curr6(+1,psiv(1,-1,2,5),pmomv(0,2,5),psi(1,-1,1),qbar(0,1),jv(1,-1,5,2))
      call curr6(+1,psi(1,-1,2),qbar(0,2),psiv(1,-1,1,5),pmomv(0,1,5),jv(1,-1,5,1))

      do isig=-1,1,2
      do i=1,2
        mat(3,isig,2) = mat(3,isig,2) 
     &             - dotcc(jv(1,isig,5,i),zww(0,1))
        mat(3,isig,1) = mat(3,isig,1) 
     &             - dotcc(jv(1,isig,5,i),aww(0,1))
      enddo
      enddo

!
! Multiply with the remaining e-w couplings
!

      do isig=-1,1,2
        do i=1,6
          mat(i,isig,3) = mat(i,isig,3) * clr(3,3,isig)**2
        enddo
        do i=1,3
! order, as we modify in-place
          mat(i+3,isig,2) = mat(i,isig,2) * clr(4,2,isig)
          mat(i+3,isig,1) = mat(i,isig,1) * clr(4,1,isig)
          mat(i  ,isig,2) = mat(i,isig,2) * clr(3,2,isig)
          mat(i  ,isig,1) = mat(i,isig,1) * clr(3,1,isig)
        enddo
      enddo

#ifdef WITH_MADGRAPH
!
! Compare with MadGraph
!
      if (ldebug) then
         print*,"MadGraph check qqWWgg"
         print*,"Permutation (selc) & gluon helicities: ", selc, l1, l2
         do j=1,3
         do i=1,6
         do isig=-1,1,2
           print*,"mat(",j,",",i,",",isig,",",l1,",",l2,") = ",mat(i,isig,j)
           temp = dreal(mat(i,isig,j))**2+dimag(mat(i,isig,j))**2
           print*,"mg2(",j,",",i,",",isig,",",l1,",",l2,") = ",
     &       -tamprr(i+6*(j-1),isig,l1*polmap(1),l2*polmap(2),1)
           temp1 = 
     &        dreal(tamprr(i+6*(j-1),isig,l1*polmap(1),l2*polmap(2),1))**2
     &       +dimag(tamprr(i+6*(j-1),isig,l1*polmap(1),l2*polmap(2),1))**2
           print*,"temp/temp1 = ",temp/temp1
         enddo
         enddo
         enddo
         read(*,*)
      endif
#endif

!
! Sum everything up
!


      do isig=-1,1,2

        do i=1,2
          qqgg(i) = 0d0
        enddo
        do i=1,3
          qqgg(1) = qqgg(1) + mat(1,isig,i) + mat(3,isig,i) 
          qqgg(2) = qqgg(2) + mat(2,isig,i) - mat(3,isig,i) 
        enddo
        m2s(1) = m2s(1) 
     &           + 16d0/3d0*(abs(qqgg(1))**2+abs(qqgg(2))**2)
     &                - 4d0/3d0*dreal(qqgg(1)*dconjg(qqgg(2)))

        do i=1,2
          qqgg(i) = 0d0
        enddo
        do i=1,3
          qqgg(1) = qqgg(1) + mat(4,isig,i) + mat(6,isig,i) 
          qqgg(2) = qqgg(2) + mat(5,isig,i) - mat(6,isig,i) 
        enddo
        m2s(2) = m2s(2) 
     &           + 16d0/3d0*(abs(qqgg(1))**2+abs(qqgg(2))**2)
     &                - 4d0/3d0*dreal(qqgg(1)*dconjg(qqgg(2)))

      enddo

      enddo     !Close Gluon 2 helicity loop HELICITY
      enddo     !Close Gluon 1 helicity loop HELICITY

      do i=1,2
        m2s(i) = m2s(i) * fpials(2)**2
      enddo

#ifdef WITH_MADGRAPH
      if (ldebug) then
         print*,"Total partonic amplitude (selc=",selc,")"
         print*,"  uubar:"
         print*,"    VBFNLO:  ", m2s(1)
         print*,"    MadGraph:", mg2uures/g**4*fpials(2)**2
         print*,"    Ratio:   ", m2s(1)/mg2uures*g**4/fpials(2)**2
         print*,"  ddbar:"
         print*,"    VBFNLO:  ", m2s(2)
         print*,"    MadGraph:", mg2ddres/g**4*fpials(2)**2
         print*,"    Ratio:   ", m2s(2)/mg2ddres*g**4/fpials(2)**2
         read(*,*)
      endif
#endif

      end


