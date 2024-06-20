      subroutine qqWHgg(p, v, selc, m2s)

      implicit none

      include 'VBFNLO/utilities/global.inc'
      include 'VBFNLO/utilities/lha.inc'
      include 'VBFNLO/utilities/scales.inc'
      include 'VBFNLO/utilities/coupl.inc'
      include 'tensor.inc'

      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1), XM2(6), XMG(6), B(6,6,6)

      integer j, i, k, l, mu, ep, isig1, isig2
      integer sign(4), polmap(2)

      double precision p(0:3, max_p), v(0:3, max_v), m2s, cpl(5)
      double precision pbar(0:3,4),
     &       qbar(0:3,4), gluvec(0:4,2), gluvecb(0:4,2), glumom(0:4)

      double complex psi(2,-1:1,4),
     &       ubdgg(3), psiv(2,-1:1,2,5), psivv(2,-1:1,2,5,5),
     &       jv(1:6,-1:1,5:5,2), jvv(1:6,-1:1,5,5,3), effgluhel(0:3) 
      double precision pmomv(0:4,2,5), pmomvv(0:4,2,5,5) 
      logical ldebug, ldebugm
      parameter (ldebug=.false.)
      parameter (ldebugm=.false.)
    
      double complex dotcc, temp, dotrc, temp1

      DOUBLE PRECISION fpials(2:3), fpi
      PARAMETER (fpi=4d0*pi)

      integer selc
      double precision pall(0:3,8), mg2res
      integer fsign(8)

      m2s=0d0

c fix strong coupling gs**2 for the two quarks:

      fpials(2) = fpi*als(1,1)
      fpials(3) = fpi*als(2,1)
c      fpials(2) = fpi*alfas

!
! Get information about different subprocesses and calculate external quark spinors

      Select Case (selc)

!-----------------------------------------------------
      Case (1) 
!  u  d~ -> e+ ve  H g g  
!( u~ d  -> e- ve~ H g g )

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
         qbar(mu, i)=sign(i)*pbar(mu, i)
       enddo
      enddo

!-----------------------------------------------------

      Case (2) 
!  u  g -> e+ ve  H d  g  
!( u~ g -> e- ve~ H d~ g )

      sign(1)=1
      sign(2)=1
      polmap(1)=-1
      polmap(2)=1

      do mu=0,3
       pbar(mu,1)=p(mu,1)
       pbar(mu,2)=p(mu,4)
       gluvecb(mu,1)=p(mu,2)
       gluvecb(mu,2)=p(mu,3)
       gluvec(mu,1)=polmap(1)*gluvecb(mu,1)
       gluvec(mu,2)=polmap(2)*gluvecb(mu,2)
       glumom(mu)=gluvec(mu,1)+gluvec(mu,2)
       do i=1,2
         qbar(mu, i)=sign(i)*pbar(mu, i)
       enddo
      enddo

!-----------------------------------------------------

      Case (3) 
!  d~ g -> e+ ve  H u~ g  
!( d  g -> e- ve~ H u  g )

      sign(1)=-1
      sign(2)=-1
      polmap(1)=-1
      polmap(2)=1

      do mu=0,3
       pbar(mu,1)=p(mu,4)
       pbar(mu,2)=p(mu,1)
       gluvecb(mu,1)=p(mu,2)
       gluvecb(mu,2)=p(mu,3)
       gluvec(mu,1)=polmap(1)*gluvecb(mu,1)
       gluvec(mu,2)=polmap(2)*gluvecb(mu,2)
       glumom(mu)=gluvec(mu,1)+gluvec(mu,2)
       do i=1,2
         qbar(mu, i)=sign(i)*pbar(mu, i)
       enddo
      enddo

!-----------------------------------------------------

      Case (4) 
!  g g -> e+ ve  H u~ d   
!( g g -> e- ve~ H u  d~ )

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
         qbar(mu, i)=sign(i)*pbar(mu, i)
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
        fsign(3)=+1
        fsign(4)=-1
        fsign(5)=+1
        do mu=0,3
          pall(mu,1) = pbar(mu,1)
          pall(mu,2) = pbar(mu,2)
          pall(mu,7) = gluvecb(mu,1)
          pall(mu,8) = gluvecb(mu,2)
          pall(mu,3) = wnumom(mu)
          pall(mu,4) = wemmom(mu)
          pall(mu,5) = higmom(mu) 
        enddo
        
        call SUD_WHGG(pall, fsign, mg2res)
      endif
#endif

!-----------------------------------------------------
!-----------------------------------------------------

      call psi0m(2, pbar(0,1), sign(1), psi)

      do isig1=-1,1,2       ! Loop over gluon (1) helicity
       do isig2 = -1,1,2    ! Loop over gluon (2) helicity

        ubdgg(1)=0d0
        ubdgg(2)=0d0
        ubdgg(3)=0d0

        call helvec(gluvecb(0,1), polmap(1), isig1, glupolc(0,1, isig1))
        call helvec(gluvecb(0,2), polmap(2), isig2, glupolc(0,2, isig2))

!
!  Add from one to four gauge bosons to the bra
!
!  in psig() etc: g1=1, g2=2, a1=3, a2=4
!  order in arguments is as in math. expression <bra|V...   ...V|ket>


! first V to <bra|
        call bra2c(psi(1,-1,2),.true., qbar(0,2),-1,
     &    gluvec(0,1), glupolc(0,1, isig1), psiv(1,-1,2,1), pmomv(0,2,1))
        call bra2c(psi(1,-1,2),.true., qbar(0,2),-1,
     &    gluvec(0,2), glupolc(0,2, isig2), psiv(1,-1,2,2), pmomv(0,2,2))

! second V to <bra|V
        do i=1,2 ! instead 4
          if (i.ne.1)
     &      call bra2c(psiv(1,-1,2, i),.false., pmomv(0,2, i),-1,
     &    gluvec(0,1), glupolc(0,1, isig1), psivv(1,-1,2, i,1), pmomvv(0,2, i,1))
          if (i.ne.2)
     &      call bra2c(psiv(1,-1,2, i),.false., pmomv(0,2, i),-1,
     &    gluvec(0,2), glupolc(0,2, isig2), psivv(1,-1,2, i,2), pmomvv(0,2, i,2))
        enddo


!
!  Add from one to four gauge bosons to the ket
!

! first V to |ket>
        call ket2c(psi(1,-1,1),.true., qbar(0,1),-1,
     &    gluvec(0,1), glupolc(0,1, isig1), psiv(1,-1,1,1), pmomv(0,1,1))
        call ket2c(psi(1,-1,1),.true., qbar(0,1),-1,
     &    gluvec(0,2), glupolc(0,2, isig2), psiv(1,-1,1,2), pmomv(0,1,2))

! second V to V|ket>
        do i=1,2 ! instead 4
          if (i.ne.1)
     &      call ket2c(psiv(1,-1,1, i),.false., pmomv(0,1, i),-1,
     &    gluvec(0,1), glupolc(0,1, isig1), psivv(1,-1,1,1, i), pmomvv(0,1,1, i))
          if (i.ne.2)
     &      call ket2c(psiv(1,-1,1, i),.false., pmomv(0,1, i),-1,
     &    gluvec(0,2), glupolc(0,2, isig2), psivv(1,-1,1,2, i), pmomvv(0,1,2, i))
        enddo


!
! Now cycle through position of W and build currents
! always need both gluons attached
! 
      do i=1,4
         do j=1,4
            if (i.eq.j) cycle
            if ((i.ne.1).and.(j.ne.1)) cycle
            if ((i.ne.2).and.(j.ne.2)) cycle
               ! i, j = 1,2 or 2,1
            call curr6(-1, psivv(1,-1,2, i, j), pmomvv(0,2, i, j), psi(1,-1,1), qbar(0,1), jvv(1,-1, i, j,3))
            call curr6(-1, psiv(1,-1,2, i), pmomv(0,2, i), psiv(1,-1,1, j), pmomv(0,1, j), jvv(1,-1, i, j,2))
            call curr6(-1, psi(1,-1,2), qbar(0,2), psivv(1,-1,1, i, j), pmomvv(0,1, i, j), jvv(1,-1, i, j,1))
         enddo
      enddo

!
!  Contract currents with the leptonic tensors
!  ubdgg: 1 = order g(1) g(2), 2 = vice versa
!

      do i=1,4
         do j=1,4
            if (i.eq.j) cycle
            if ((i.ne.1).and.(j.ne.1)) cycle
            if ((i.ne.2).and.(j.ne.2)) cycle
            k=0
            l=0
            call checkorder(i, j, k, l, ep, cpl)
            do mu=1,3
               temp = dotcc(jvv(1,-1, i, j, mu), wtowh(0,1))
               ubdgg(ep) = ubdgg(ep) + cpl(mu)*temp
            enddo
         enddo
      enddo
      

! 
!   Get effective g(1)g(2) -> g polarization
!
        do mu=0,3
        effgluhel(mu) = (2d0*dotrc(gluvec(0,1), glupolc(0,2, isig2))*glupolc(mu,1, isig1)
     &       -2d0*dotrc(gluvec(0,2), glupolc(0,1, isig1))*glupolc(mu,2, isig2)
     &       +dotcc(glupolc(0,1, isig1), glupolc(0,2, isig2))*(gluvec(mu,2)-gluvec(mu,1)))
     &       /(glumom(4))
        enddo

!
! Now add again the vector bosons to the <bra| and |ket>
!
! first V to <bra|
        call bra2c(psi(1,-1,2),.true., qbar(0,2),-1,
     &    glumom(0), effgluhel(0), psiv(1,-1,2,5), pmomv(0,2,5)) !5

        do i=3,4
          call bra2c(psiv(1,-1,2, i),.false., pmomv(0,2, i),-1,
     &      glumom(0), effgluhel(0), psivv(1,-1,2, i,5), pmomvv(0,2, i,5)) !i5
        enddo
!

! first V to |ket>
        call ket2c(psi(1,-1,1),.true., qbar(0,1),-1,
     &    glumom(0), effgluhel(0), psiv(1,-1,1,5), pmomv(0,1,5)) !5

!! second V to V|ket>
        do i=3,4
          call ket2c(psiv(1,-1,1, i),.false., pmomv(0,1, i),-1,
     &      glumom(0), effgluhel(0), psivv(1,-1,1,5, i), pmomvv(0,1,5, i)) !5i
        enddo
!

!
! Now cycle through position of W and build currents
! always need three-gluon-vertex (5) attached
! 
      do i=3,5
        if (i.ne.5) cycle
        call curr6(-1, psiv(1,-1,2, i), pmomv(0,2, i), psi(1,-1,1), qbar(0,1), jv(1,-1, i,2))
        call curr6(-1, psi(1,-1,2), qbar(0,2), psiv(1,-1,1, i), pmomv(0,1, i), jv(1,-1, i,1))
      enddo

      do i=3,5
         if (i.ne.5) cycle
         do mu=1,2
            cpl(mu) = clr(1,3,-1)
            temp = dotcc(jv(1,-1, i, mu), wtowh(0,1))
            ubdgg(3) = ubdgg(3) + cpl(mu)*temp
         enddo
      enddo

#ifdef WITH_MADGRAPH
!
! Compare with MadGraph
!
      if (ldebugm) then
         do i=1,3
           temp = dreal(ubdgg(i))**2+dimag(ubdgg(i))**2
           temp1 = 
     &        dreal(tamprr(i,-1, isig1*polmap(1), isig2*polmap(2),1))**2
     &       +dimag(tamprr(i,-1, isig1*polmap(1), isig2*polmap(2),1))**2
           if(abs(temp/temp1 -1).ge.1.0d-6) then
              print*,"MadGraph check qqWHgg"
              print*,"Type & gluon helicities: ", selc, isig1, isig2
              print*,"mat(", i,",", isig1,",", isig2,") = ", ubdgg(i)
              print*,"mg2(", i,",", isig1,",", isig2,") = ",
     &       -tamprr(i,-1, isig1*polmap(1), isig2*polmap(2),1)
              print*,"temp/temp1 = ", temp/temp1
              print*, "matrix elements qqWHgg.F are only within", abs(temp/temp1-1) 
           endif
         enddo
!         read(*,*)
      endif
#endif

!
! Sum everything up
!

      ubdgg(1) = ubdgg(1) + ubdgg(3)
      ubdgg(2) = ubdgg(2) - ubdgg(3)

      m2s=m2s+16d0/3d0*(abs(ubdgg(1))**2+abs(ubdgg(2))**2)
     &              - 4d0/3d0*dreal(ubdgg(1)*dconjg(ubdgg(2)))

      enddo     !Close Gluon 2 helicity loop HELICITY
      enddo     !Close Gluon 1 helicity loop HELICITY

      m2s = m2s * fpials(2)**2

      end


