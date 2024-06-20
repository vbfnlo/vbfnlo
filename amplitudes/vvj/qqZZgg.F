      subroutine qqZAZAgg(p,v,selc,bos,m2s)

      implicit none

      include 'VBFNLO/utilities/global.inc'
      include 'VBFNLO/utilities/lha.inc'
      include 'VBFNLO/utilities/scales.inc'
      include 'VBFNLO/utilities/coupl.inc'
      include 'tensor.inc'

      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)

      integer bos,j,i,k,l,mu,isig,l1,l2
      integer sign(4),polmap(2)

      double precision p(0:3,max_p), v(0:3,max_v), m2s(2)
      double precision pbar(0:3,4), qbar(0:3,4),
     &       gluvecb(0:4,2), gluvec(0:4,2), glumom(0:4)
      double precision pmomv(0:4,2,5),pmomvv(0:4,2,5,5),
     &       pmomvvv(0:4,2,5,5,5)

      double complex mat(6,-1:1,6), qqgg(2)
      double complex psi(2,-1:1,2),
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

c  helicity selection
      INTEGER h
      COMMON /hcount / h

c variables for ZZ,ZA,AA 
      integer typ, typl
      integer cpp(4), cpm(4)
      double complex cepsp(6,4), cepsm(6,4)
      logical lhaszvv, linit
      data linit /.true./
      save linit

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
        call SUU_ZZGG(pall,fsign,mg2uures)
        call SDD_ZZGG(pall,fsign,mg2ddres)
      endif
#endif

!-----------------------------------------------------

c map polarisation vectors and indices onto the arrays, so we can use do-loops
c makes debugging way easier
c typ: 1:AA, 2: ZA, 3: AZ, 4: ZZ (so can just stop at right index for AA and ZA)
      if (linit) then
        cpp(1) = 1
        cpp(2) = 2
        cpp(3) = 1
        cpp(4) = 2
        cpm(1) = 1
        cpm(2) = 1
        cpm(3) = 2
        cpm(4) = 2
        if (bos.eq.22) then
          typ = 4
          lhaszvv = .true.
        else if (bos.eq.21) then
          typ = 2
          lhaszvv = .true.
        else if (bos.eq.11) then
          typ = 1
          lhaszvv = .false.
        else
          write(*,*) "Error in ggZZg: wrong bos ",bos
          stop
        endif
c        if (narrowwidth.eq.0) then
          typl = 1
c        else
c          typl = typ
c        endif
        linit = .false.
      endif
      do mu=1,6
        cepsp(mu,1) = atau1(mu,1)
        cepsp(mu,2) = wp(mu,1)
        cepsp(mu,3) = atau1(mu,1)
        cepsp(mu,4) = wp(mu,1)
        cepsm(mu,1) = atau2(mu,1)
        cepsm(mu,2) = atau2(mu,1)
        cepsm(mu,3) = wm(mu,1)
        cepsm(mu,4) = wm(mu,1)
      enddo
      if (bos.eq.11) then
        do mu=0,3
          qp(mu,1) = qa1(mu,1)
          qm(mu,1) = qa2(mu,1)
        enddo
      endif

!-----------------------------------------------------

      call psi0m(2,pbar(0,1),sign(1),psi)

      do l1=-1,1,2       ! Loop over gluon (1) helicity
       do l2 = -1,1,2    ! Loop over gluon (2) helicity

        do j=1,6
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
!  in psig() etc: g1=1, g2=2, A_mu=3, Z_mu=4
!  order in arguments is as in math. expression <bra|V...   ...V|ket>
!

       do isig=-1,1,2
! first V to <bra|
        call bra2c(psi(1,isig,2),.true.,qbar(0,2),isig,
     &    gluvec(0,1),glupolc(0,1,l1),psiv(1,isig,2,1),pmomv(0,2,1))
        call bra2c(psi(1,isig,2),.true.,qbar(0,2),isig,
     &    gluvec(0,2),glupolc(0,2,l2),psiv(1,isig,2,2),pmomv(0,2,2))
        call bra2c(psi(1,isig,2),.true.,qbar(0,2),isig,
     &    qm(0,1),atau2(1,1),psiv(1,isig,2,3),pmomv(0,2,3))
        if (typ.ge.3) 
     &    call bra2c(psi(1,isig,2),.true.,qbar(0,2),isig,
     &      qm(0,1),wm(1,1),psiv(1,isig,2,4),pmomv(0,2,4))

! second V to <bra|V
        do i=1,4
          if (i.ne.1)
     &      call bra2c(psiv(1,isig,2,i),.false.,pmomv(0,2,i),isig,
     &        gluvec(0,1),glupolc(0,1,l1),psivv(1,isig,2,i,1),pmomvv(0,2,i,1))
          if (i.ne.2)
     &      call bra2c(psiv(1,isig,2,i),.false.,pmomv(0,2,i),isig,
     &        gluvec(0,2),glupolc(0,2,l2),psivv(1,isig,2,i,2),pmomvv(0,2,i,2))
          if ((i.ne.3).and.(i.ne.4))
     &      call bra2c(psiv(1,isig,2,i),.false.,pmomv(0,2,i),isig,
     &        qm(0,1),atau2(1,1),psivv(1,isig,2,i,3),pmomvv(0,2,i,3))
          if ((i.ne.3).and.(i.ne.4).and.(typ.ge.3))
     &      call bra2c(psiv(1,isig,2,i),.false.,pmomv(0,2,i),isig,
     &        qm(0,1),wm(1,1),psivv(1,isig,2,i,4),pmomvv(0,2,i,4))
        enddo

! third V to <bra|VV
        do i=1,4
        do j=1,4
          if ((i.ne.1).and.(j.ne.1))
     &      call bra2c(psivv(1,isig,2,i,j),.false.,pmomvv(0,2,i,j),isig,
     &        gluvec(0,1),glupolc(0,1,l1),psivvv(1,isig,2,i,j,1),pmomvvv(0,2,i,j,1))
          if ((i.ne.2).and.(j.ne.2))
     &      call bra2c(psivv(1,isig,2,i,j),.false.,pmomvv(0,2,i,j),isig,
     &        gluvec(0,2),glupolc(0,2,l2),psivvv(1,isig,2,i,j,2),pmomvvv(0,2,i,j,2))
          if ((i.ne.3).and.(j.ne.3).and.(i.ne.4).and.(j.ne.4))
     &      call bra2c(psivv(1,isig,2,i,j),.false.,pmomvv(0,2,i,j),isig,
     &        qm(0,1),atau2(1,1),psivvv(1,isig,2,i,j,3),pmomvvv(0,2,i,j,3))
          if ((i.ne.3).and.(j.ne.3).and.(i.ne.4).and.(j.ne.4).and.(typ.ge.3))
     &      call bra2c(psivv(1,isig,2,i,j),.false.,pmomvv(0,2,i,j),isig,
     &        qm(0,1),wm(1,1),psivvv(1,isig,2,i,j,4),pmomvvv(0,2,i,j,4))
        enddo
        enddo

!
!  Add from one to three gauge bosons to the ket
!

! first V to |ket>
        call ket2c(psi(1,isig,1),.true.,qbar(0,1),isig,
     &    gluvec(0,1),glupolc(0,1,l1),psiv(1,isig,1,1),pmomv(0,1,1))
        call ket2c(psi(1,isig,1),.true.,qbar(0,1),isig,
     &    gluvec(0,2),glupolc(0,2,l2),psiv(1,isig,1,2),pmomv(0,1,2))
        call ket2c(psi(1,isig,1),.true.,qbar(0,1),isig,
     &    qm(0,1),atau2(1,1),psiv(1,isig,1,3),pmomv(0,1,3))
        if (typ.ge.3) 
     &    call ket2c(psi(1,isig,1),.true.,qbar(0,1),isig,
     &      qm(0,1),wm(1,1),psiv(1,isig,1,4),pmomv(0,1,4))

! second V to V|ket>
        do i=1,4
          if (i.ne.1)
     &      call ket2c(psiv(1,isig,1,i),.false.,pmomv(0,1,i),isig,
     &        gluvec(0,1),glupolc(0,1,l1),psivv(1,isig,1,1,i),pmomvv(0,1,1,i))
          if (i.ne.2)
     &      call ket2c(psiv(1,isig,1,i),.false.,pmomv(0,1,i),isig,
     &        gluvec(0,2),glupolc(0,2,l2),psivv(1,isig,1,2,i),pmomvv(0,1,2,i))
          if ((i.ne.3).and.(i.ne.4))
     &      call ket2c(psiv(1,isig,1,i),.false.,pmomv(0,1,i),isig,
     &        qm(0,1),atau2(1,1),psivv(1,isig,1,3,i),pmomvv(0,1,3,i))
          if ((i.ne.3).and.(i.ne.4).and.(typ.ge.3))
     &      call ket2c(psiv(1,isig,1,i),.false.,pmomv(0,1,i),isig,
     &        qm(0,1),wm(1,1),psivv(1,isig,1,4,i),pmomvv(0,1,4,i))
        enddo

! third V to VV|ket>
        do i=1,4
        do j=1,4
          if ((i.ne.1).and.(j.ne.1))
     &      call ket2c(psivv(1,isig,1,i,j),.false.,pmomvv(0,1,i,j),isig,
     &        gluvec(0,1),glupolc(0,1,l1),psivvv(1,isig,1,1,i,j),pmomvvv(0,1,1,i,j))
          if ((i.ne.2).and.(j.ne.2))
     &      call ket2c(psivv(1,isig,1,i,j),.false.,pmomvv(0,1,i,j),isig,
     &        gluvec(0,2),glupolc(0,2,l2),psivvv(1,isig,1,2,i,j),pmomvvv(0,1,2,i,j))
          if ((i.ne.3).and.(j.ne.3).and.(i.ne.4).and.(j.ne.4))
     &      call ket2c(psivv(1,isig,1,i,j),.false.,pmomvv(0,1,i,j),isig,
     &        qm(0,1),atau2(1,1),psivvv(1,isig,1,3,i,j),pmomvvv(0,1,3,i,j))
          if ((i.ne.3).and.(j.ne.3).and.(i.ne.4).and.(j.ne.4).and.(typ.ge.3))
     &      call ket2c(psivv(1,isig,1,i,j),.false.,pmomvv(0,1,i,j),isig,
     &        qm(0,1),wm(1,1),psivvv(1,isig,1,4,i,j),pmomvvv(0,1,4,i,j))
        enddo
        enddo

      enddo ! isig

!
! Now cycle through position of V_mu and build currents
! always need both gluons attached
! 
      do i=1,4
        if ((i.eq.4).and.(typ.lt.3)) cycle
      do j=1,4
        if ((j.eq.4).and.(typ.lt.3)) cycle
        if (i.eq.j) cycle
        if ((i.eq.3).and.(j.eq.4)) cycle
        if ((i.eq.4).and.(j.eq.3)) cycle
      do k=1,4
        if ((k.eq.4).and.(typ.lt.3)) cycle
        if ((i.eq.k).or.(j.eq.k)) cycle
        if ((i.eq.3).and.(k.eq.4)) cycle
        if ((i.eq.4).and.(k.eq.3)) cycle
        if ((j.eq.3).and.(k.eq.4)) cycle
        if ((j.eq.4).and.(k.eq.3)) cycle
! all 4 bosons must be attached 
        call curr6(+1,psivvv(1,-1,2,i,j,k),pmomvvv(0,2,i,j,k),psi(1,-1,1),qbar(0,1),jvvv(1,-1,i,j,k,4))
        call curr6(+1,psivv(1,-1,2,i,j),pmomvv(0,2,i,j),psiv(1,-1,1,k),pmomv(0,1,k),jvvv(1,-1,i,j,k,3))
        call curr6(+1,psiv(1,-1,2,i),pmomv(0,2,i),psivv(1,-1,1,j,k),pmomvv(0,1,j,k),jvvv(1,-1,i,j,k,2))
        call curr6(+1,psi(1,-1,2),qbar(0,2),psivvv(1,-1,1,i,j,k),pmomvvv(0,1,i,j,k),jvvv(1,-1,i,j,k,1))
      enddo
      enddo
      enddo

!
! Now (VV) attached to fermion line
!
      if (lhaszvv) then
        do i=1,2
          j=3-i
          call curr6(+1,psivv(1,-1,2,i,j),pmomvv(0,2,i,j),psi(1,-1,1),qbar(0,1),jvv(1,-1,i,j,3))
          call curr6(+1,psiv(1,-1,2,i),pmomv(0,2,i),psiv(1,-1,1,j),pmomv(0,1,j),jvv(1,-1,i,j,2))
          call curr6(+1,psi(1,-1,2),qbar(0,2),psivv(1,-1,1,i,j),pmomvv(0,1,i,j),jvv(1,-1,i,j,1))
        enddo
      endif

!
!  Contract currents with the leptonic tensors
!

      do i=1,4
        if ((i.eq.4).and.(typ.lt.3)) cycle
      do j=1,4
        if ((j.eq.4).and.(typ.lt.3)) cycle
        if (i.eq.j) cycle
        if ((i.eq.3).and.(j.eq.4)) cycle
        if ((i.eq.4).and.(j.eq.3)) cycle
      do k=1,4
        if ((k.eq.4).and.(typ.lt.3)) cycle
        if ((i.eq.k).or.(j.eq.k)) cycle
        if ((i.eq.3).and.(k.eq.4)) cycle
        if ((i.eq.4).and.(k.eq.3)) cycle
        if ((j.eq.3).and.(k.eq.4)) cycle
        if ((j.eq.4).and.(k.eq.3)) cycle
      do l=1,4 ! position for inserting V_e
! check for order of g's for colour structure
        if ((i.eq.1).or. 
     &      ((i.eq.1).and.(j.eq.2)).or. 
     &      (k.eq.2)) then
          cs = 1  ! 1-2
        else
          cs = 2  ! 2-1
        endif

! build matrix elements
        do isig=-1,1,2
          if ((i.eq.3).or.(j.eq.3).or.(k.eq.3)) then 
c A_mu
cc + A_e
            mat(cs,isig,1) = mat(cs,isig,1) 
     &               + dotcc(jvvv(1,isig,i,j,k,l),atau1(1,1))
            if (typ.ge.2) 
cc + Z_e
     &        mat(cs,isig,2) = mat(cs,isig,2) 
     &               + dotcc(jvvv(1,isig,i,j,k,l),wp(1,1))
          elseif (typ.ge.3) then
c Z_mu
cc + A_e
            mat(cs,isig,3) = mat(cs,isig,3) 
     &               + dotcc(jvvv(1,isig,i,j,k,l),atau1(1,1))
cc + Z_e
            mat(cs,isig,4) = mat(cs,isig,4) 
     &               + dotcc(jvvv(1,isig,i,j,k,l),wp(1,1))
          endif
        enddo

      enddo
      enddo
      enddo
      enddo

      if (lhaszvv) then
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
          mat(cs,isig,5) = mat(cs,isig,5) 
     &               - dotcc(jvv(1,isig,i,j,k),wz(1,1))
          mat(cs,isig,6) = mat(cs,isig,6) 
     &               - dotcc(jvv(1,isig,i,j,k),wa(1,1))
        enddo
        enddo
        enddo
      endif
      

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
! eff gg to <bra|
        call bra2c(psi(1,isig,2),.true.,qbar(0,2),isig,
     &    glumom(0),effgluhel(0),psiv(1,isig,2,5),pmomv(0,2,5)) !5
! V_mu to <bra|(gg)
        call bra2c(psiv(1,isig,2,5),.false.,pmomv(0,2,5),isig,
     &    qm(0,1),atau2(1,1),psivv(1,isig,2,5,3),pmomvv(0,2,5,3))
        if (typ.ge.3)
     &    call bra2c(psiv(1,isig,2,5),.false.,pmomv(0,2,5),isig,
     &      qm(0,1),wm(1,1),psivv(1,isig,2,5,4),pmomvv(0,2,5,4))
! eff gg to <bra|V_mu
        call bra2c(psiv(1,isig,2,3),.false.,pmomv(0,2,3),isig,
     &    glumom(0),effgluhel(0),psivv(1,isig,2,3,5),pmomvv(0,2,3,5))
        if (typ.ge.3)
     &    call bra2c(psiv(1,isig,2,4),.false.,pmomv(0,2,4),isig,
     &      glumom(0),effgluhel(0),psivv(1,isig,2,4,5),pmomvv(0,2,4,5))

! eff gg to |ket>
        call ket2c(psi(1,isig,1),.true.,qbar(0,1),isig,
     &    glumom(0),effgluhel(0),psiv(1,isig,1,5),pmomv(0,1,5)) !5
! V_mu to (gg)|ket>
        call ket2c(psiv(1,isig,1,5),.false.,pmomv(0,1,5),isig,
     &    qm(0,1),atau2(1,1),psivv(1,isig,1,3,5),pmomvv(0,1,3,5))
        if (typ.ge.3)
     &    call ket2c(psiv(1,isig,1,5),.false.,pmomv(0,1,5),isig,
     &      qm(0,1),wm(1,1),psivv(1,isig,1,4,5),pmomvv(0,1,4,5))
! eff gg to V_mu|ket>
        call ket2c(psiv(1,isig,1,3),.false.,pmomv(0,1,3),isig,
     &    glumom(0),effgluhel(0),psivv(1,isig,1,5,3),pmomvv(0,1,5,3))
        if (typ.ge.3)
     &    call ket2c(psiv(1,isig,1,4),.false.,pmomv(0,1,4),isig,
     &      glumom(0),effgluhel(0),psivv(1,isig,1,5,4),pmomvv(0,1,5,4))

       enddo ! isig


!
! Now cycle through position of V_e and build currents
! 
      do i=3,5,2
        j=8-i
        call curr6(+1,psivv(1,-1,2,i,j),pmomvv(0,2,i,j),psi(1,-1,1),qbar(0,1),jvv(1,-1,i,j,3))
        call curr6(+1,psiv(1,-1,2,i),pmomv(0,2,i),psiv(1,-1,1,j),pmomv(0,1,j),jvv(1,-1,i,j,2))
        call curr6(+1,psi(1,-1,2),qbar(0,2),psivv(1,-1,1,i,j),pmomvv(0,1,i,j),jvv(1,-1,i,j,1))
      enddo
      if (typ.ge.3) then
        do i=4,5,1
          j=9-i
          call curr6(+1,psivv(1,-1,2,i,j),pmomvv(0,2,i,j),psi(1,-1,1),qbar(0,1),jvv(1,-1,i,j,3))
          call curr6(+1,psiv(1,-1,2,i),pmomv(0,2,i),psiv(1,-1,1,j),pmomv(0,1,j),jvv(1,-1,i,j,2))
          call curr6(+1,psi(1,-1,2),qbar(0,2),psivv(1,-1,1,i,j),pmomvv(0,1,i,j),jvv(1,-1,i,j,1))
        enddo
      endif

      do isig=-1,1,2
        do k=1,3
c A_mu
          do i=3,5,2
            j=8-i
cc + A_e
            mat(3,isig,1) = mat(3,isig,1) 
     &               + dotcc(jvv(1,isig,i,j,k),atau1(1,1))
            if (typ.ge.2) 
cc + Z_e
     &        mat(3,isig,2) = mat(3,isig,2) 
     &               + dotcc(jvv(1,isig,i,j,k),wp(1,1))
          enddo
          if (typ.ge.3) then
c Z_mu
            do i=4,5,1
              j=9-i
cc + A_e
              mat(3,isig,3) = mat(3,isig,3) 
     &                 + dotcc(jvv(1,isig,i,j,k),atau1(1,1))
cc + Z_e
              mat(3,isig,4) = mat(3,isig,4) 
     &               + dotcc(jvv(1,isig,i,j,k),wp(1,1))
            enddo
          endif

        enddo
      enddo

      if (lhaszvv) then
        call curr6(+1,psiv(1,-1,2,5),pmomv(0,2,5),psi(1,-1,1),qbar(0,1),jv(1,-1,5,2))
        call curr6(+1,psi(1,-1,2),qbar(0,2),psiv(1,-1,1,5),pmomv(0,1,5),jv(1,-1,5,1))

        do isig=-1,1,2
        do i=1,2
          mat(3,isig,5) = mat(3,isig,5) 
     &               - dotcc(jv(1,isig,5,i),wz(1,1))
          mat(3,isig,6) = mat(3,isig,6) 
     &               - dotcc(jv(1,isig,5,i),wa(1,1))
        enddo
        enddo
      endif

!
! Multiply with the remaining e-w couplings
!

      do isig=-1,1,2
        do i=1,3
! order, as we modify in-place
          do j=1,4
            mat(i+3,isig,j) = mat(i,isig,j) * clr(4,cpp(j),isig) * clr(4,cpm(j),isig) 
            mat(i  ,isig,j) = mat(i,isig,j) * clr(3,cpp(j),isig) * clr(3,cpm(j),isig) 
          enddo
          mat(i+3,isig,5) = mat(i,isig,5) * clr(4,2,isig)
          mat(i+3,isig,6) = mat(i,isig,6) * clr(4,1,isig)
          mat(i  ,isig,5) = mat(i,isig,5) * clr(3,2,isig)
          mat(i  ,isig,6) = mat(i,isig,6) * clr(3,1,isig)
        enddo
      enddo

#ifdef WITH_MADGRAPH
!
! Compare with MadGraph
!
      if (ldebug) then
         print*,"MadGraph check qqZZgg"
         print*,"Permutation (selc) & gluon helicities: ", selc, l1, l2
         do j=1,6
         do i=1,6
         do isig=-1,1,2
           print*,"mat(",j,",",i,",",isig,",",l1,",",l2,") = ",mat(i,isig,j)
           temp = dreal(mat(i,isig,j))**2+dimag(mat(i,isig,j))**2
           print*,"mg2(",j,",",i,",",isig,",",l1,",",l2,") = ",
     &       -tamprr(i+6*(j-1),isig,l1*polmap(1),l2*polmap(2),h)
           temp1 = 
     &        dreal(tamprr(i+6*(j-1),isig,l1*polmap(1),l2*polmap(2),h))**2
     &       +dimag(tamprr(i+6*(j-1),isig,l1*polmap(1),l2*polmap(2),h))**2
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
        do i=1,6
          qqgg(1) = qqgg(1) + mat(1,isig,i) + mat(3,isig,i) 
          qqgg(2) = qqgg(2) + mat(2,isig,i) - mat(3,isig,i) 
        enddo
        m2s(1) = m2s(1) 
     &           + 16d0/3d0*(abs(qqgg(1))**2+abs(qqgg(2))**2)
     &                - 4d0/3d0*dreal(qqgg(1)*dconjg(qqgg(2)))

        do i=1,2
          qqgg(i) = 0d0
        enddo
        do i=1,6
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


