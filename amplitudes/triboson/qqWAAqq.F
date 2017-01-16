      subroutine qqWAAqq(p,v,selc,m2su,m2sd)

      implicit none

      include 'VBFNLO/utilities/global.inc'
      include 'VBFNLO/utilities/lha.inc'
      include 'VBFNLO/utilities/scales.inc'
      include 'VBFNLO/utilities/coupl.inc'
      include 'tensor.inc'

      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)

      integer j,i,k,mu,chi

      integer sign(4),iflavu,iflavd


      double precision p(0:3,max_p), v(0:3,max_v),m2su(3), m2sd(3)
      double precision pbar(0:3,4),qbar(0:3,4)


      double complex psi(2,-1:1,4),
     &       uuud(2,-1:1,-1:1), uddd(2,-1:1,-1:1), 
     &       psiv(2,-1:1,4,5),psivv(2,-1:1,4,5,5),
     &       psivvw(2,-1:1,4), psivw(2,-1:1,4,3:4), psinovw(2,-1:1,4),
     &       jnovw(1:6,-1:1,3), jvw(1:6,-1:1,3,3:4), jvvw(1:6,-1:1,3),
     &       gvv14(1:6,-1:1), gvv32(1:6,-1:1), gvv34(1:6,-1:1),
     &       gv14(1:6,-1:1,3:4), gv32(1:6,-1:1,3:4), gv34(1:6,-1:1,3:4),
     &       gnov14(1:6,-1:1), gnov32(1:6,-1:1), gnov34(1:6,-1:1)

      double precision pmomv(0:4,4,5),pmomvv(0:4,4,5,5),
     &       gmomnov14(0:4,-1:1),gmomnov32(0:4,-1:1),gmomnov34(0:4,-1:1),
     &       gmomv14(0:4,-1:1,3:4),gmomv32(0:4,-1:1,3:4),gmomv34(0:4,-1:1,3:4),
     &       gmomvv14(0:4,-1:1),gmomvv32(0:4,-1:1),gmomvv34(0:4,-1:1),
     &       pmomvvw(0:4,-1:1,4),pmomvw(0:4,-1:1,4,3:4),pmomnovw(0:4,-1:1,4)

      double complex im
      logical ldebug
      parameter (ldebug=.false.)
      parameter (im = (0d0,1d0))
    
      double complex dotcc

      double complex temp, temp1


      DOUBLE PRECISION fpials(2:3), fpi
      PARAMETER (fpi=4d0*pi)

      integer selc

      double precision pall(0:3,8), mg2res
      integer fsign(8)

      do i=1,3
        m2su(i)=0d0
        m2sd(i)=0d0
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
!  u  u  > ve e+  a a d  u
!  u  d  > ve e+  a a d  d
!( u~ u~ > e- ve~ a a d~ u~ )
!( u~ d~ > e- ve~ a a d~ d~ )
!
      iflavu=3
      iflavd=4
      sign(1)=+1
      sign(2)=+1
      sign(3)=+1
      sign(4)=+1
      do mu=0,3
      pbar(mu,1)=p(mu,1)   ! 1-> ket
      pbar(mu,2)=p(mu,3)   ! 2-> bra
      pbar(mu,3)=p(mu,2)   ! 3-> ket
      pbar(mu,4)=p(mu,4)   ! 4-> bra
      enddo
!
!-----------------------------------------------------
      Case (2)  
!  u  u~ > ve e+  a a d  u~
!  u  d~ > ve e+  a a d  d~
!( u~ u  > e- ve~ a a d~ u  )
!( u~ d  > e- ve~ a a d~ d  )
!
      iflavu=3
      iflavd=4
      sign(1)=1
      sign(2)=1  
      sign(3)=-1
      sign(4)=-1
      do mu=0,3
      pbar(mu,1)=p(mu,1)   ! 1-> ket
      pbar(mu,2)=p(mu,3)   ! 2-> bra
      pbar(mu,3)=p(mu,4)   ! 3-> ket
      pbar(mu,4)=p(mu,2)   ! 4-> bra
      enddo

!
!-----------------------------------------------------
      Case (3)  
!  u  d~ > ve e+  a a u~ u  
!  d  d~ > ve e+  a a u~ d  
!( u~ d  > e- ve~ a a u  u~ )
!( d~ d  > e- ve~ a a u  d~ )
!
      iflavu=3
      iflavd=4
      sign(1)=-1
      sign(2)=-1
      sign(3)=1
      sign(4)=1
      do mu=0,3
      pbar(mu,1)=p(mu,3)   ! 1-> ket  
      pbar(mu,2)=p(mu,2)   ! 2-> bra
      pbar(mu,3)=p(mu,1)   ! 3-> ket
      pbar(mu,4)=p(mu,4)   ! 4-> bra
      enddo
!
!-----------------------------------------------------
      Case (4)  
!  u~ d~ > ve e+  a a u~ u~
!  d~ d~ > ve e+  a a u~ d~
!( u  d  > e- ve~ a a u  u  )
!( d  d  > e- ve~ a a u  d  )
!
      iflavu=3
      iflavd=4
      sign(1)=-1
      sign(2)=-1
      sign(3)=-1
      sign(4)=-1
      do mu=0,3
      pbar(mu,1)=p(mu,3)   ! 1-> ket
      pbar(mu,2)=p(mu,2)   ! 2-> bra
      pbar(mu,3)=p(mu,4)   ! 3-> ket
      pbar(mu,4)=p(mu,1)   ! 4-> bra
      enddo

      end select

      do mu=0,3
       do i=1,4
         qbar(mu,i)=sign(i)*pbar(mu,i)
       enddo
      enddo

#ifdef WITH_MADGRAPH
!
! Write settings in variables for MadGraph
!
      if (ldebug) then
        fsign(1)= sign(1)
        fsign(2)= sign(2)
        fsign(7)= sign(4)
        fsign(8)= sign(3)
        fsign(3)=+1
        fsign(4)=-1
        fsign(5)=+1
        fsign(6)=+1
        do mu=0,3
          pall(mu,1) = pbar(mu,1)
          pall(mu,2) = pbar(mu,2)
          pall(mu,7) = pbar(mu,4)
          pall(mu,8) = pbar(mu,3)
          pall(mu,3) = wnumom(mu)
          pall(mu,4) = wemmom(mu)
          pall(mu,5) = phot1fourvec(mu)
          pall(mu,6) = phot2fourvec(mu)
        enddo
        call SUD_WAAUU(pall,fsign,mg2res)
        call SUD_WAADD(pall,fsign,mg2res)
      endif
#endif

!-----------------------------------------------------
!-----------------------------------------------------

      call psi0m(4,pbar(0,1),sign(1),psi)
      
      do i=1,2
      do j=-1,1,2
      do k=-1,1,2
        uuud(i,j,k)=0d0
        uddd(i,j,k)=0d0
      enddo
      enddo
      enddo

!
!  First construct the fermion line which has no W attached
!
!  Add up to two photons to psi
!
!  in psiv() etc: a1=3, a2=4
!  order in arguments is as in math. expression <bra|V...   ...V|ket>
!

        do chi=-1,1,2
          do i=1,3,2
            j=i+1
! first V to <bra| or |ket>
            call bra2c(psi(1,chi,j),.true.,qbar(0,j),chi,
     &        phot1fourvec(0),atau1(1,1),psiv(1,chi,j,3),pmomv(0,j,3))
            call ket2c(psi(1,chi,i),.true.,qbar(0,i),chi,
     &        phot1fourvec(0),atau1(1,1),psiv(1,chi,i,3),pmomv(0,i,3))
            call bra2c(psi(1,chi,j),.true.,qbar(0,j),chi,
     &        phot2fourvec(0),atau2(1,1),psiv(1,chi,j,4),pmomv(0,j,4))
            call ket2c(psi(1,chi,i),.true.,qbar(0,i),chi,
     &        phot2fourvec(0),atau2(1,1),psiv(1,chi,i,4),pmomv(0,i,4))

! second V to <bra|V or V|ket>
            call bra2c(psiv(1,chi,j,4),.false.,pmomv(0,j,4),chi,
     &        phot1fourvec(0),atau1(1,1),psivv(1,chi,j,4,3),pmomvv(0,j,4,3))
            call bra2c(psiv(1,chi,j,3),.false.,pmomv(0,j,3),chi,
     &        phot2fourvec(0),atau2(1,1),psivv(1,chi,j,3,4),pmomvv(0,j,3,4))
            call ket2c(psiv(1,chi,i,4),.false.,pmomv(0,i,4),chi,
     &        phot1fourvec(0),atau1(1,1),psivv(1,chi,i,3,4),pmomvv(0,i,3,4))
            call ket2c(psiv(1,chi,i,3),.false.,pmomv(0,i,3),chi,
     &        phot2fourvec(0),atau2(1,1),psivv(1,chi,i,4,3),pmomvv(0,i,4,3))

          enddo
        enddo

!
! Now cycle through position of gluon and build currents
! <2|1> is the <d|u> quark line in the common combination
! at the end add gluon propagator
! 
! initialize gvv14, gvv32 and gvv34 first
        do i=-1,1,2
        do mu=1,4
          gvv14(mu,i) = 0d0
          gvv32(mu,i) = 0d0
          gvv34(mu,i) = 0d0
        enddo
          gvv14(5,i) = dcmplx(
     &      qbar(0,1)-qbar(0,4)-phot1fourvec(0)-phot2fourvec(0),
     &      qbar(3,1)-qbar(3,4)-phot1fourvec(3)-phot2fourvec(3))
          gvv14(6,i) = dcmplx(
     &      qbar(1,1)-qbar(1,4)-phot1fourvec(1)-phot2fourvec(1),
     &      qbar(2,1)-qbar(2,4)-phot1fourvec(2)-phot2fourvec(2))
          gvv32(5,i) = dcmplx(
     &      qbar(0,3)-qbar(0,2)-phot1fourvec(0)-phot2fourvec(0),
     &      qbar(3,3)-qbar(3,2)-phot1fourvec(3)-phot2fourvec(3))
          gvv32(6,i) = dcmplx(
     &      qbar(1,3)-qbar(1,2)-phot1fourvec(1)-phot2fourvec(1),
     &      qbar(2,3)-qbar(2,2)-phot1fourvec(2)-phot2fourvec(2))
          gvv34(5,i) = dcmplx(
     &      qbar(0,3)-qbar(0,4)-phot1fourvec(0)-phot2fourvec(0),
     &      qbar(3,3)-qbar(3,4)-phot1fourvec(3)-phot2fourvec(3))
          gvv34(6,i) = dcmplx(
     &      qbar(1,3)-qbar(1,4)-phot1fourvec(1)-phot2fourvec(1),
     &      qbar(2,3)-qbar(2,4)-phot1fourvec(2)-phot2fourvec(2))
        enddo
        do i=3,4
          j = 7-i
! <4|1>
          call curr6add(+1,psivv(1,-1,4,i,j),pmomvv(0,4,i,j),psi(1,-1,1),qbar(0,1),gvv14(1,-1))
          call curr6add(+1,psiv(1,-1,4,i),pmomv(0,4,i),psiv(1,-1,1,j),pmomv(0,1,j),gvv14(1,-1))
          call curr6add(+1,psi(1,-1,4),qbar(0,4),psivv(1,-1,1,i,j),pmomvv(0,1,i,j),gvv14(1,-1))
! <2|3>
          call curr6add(+1,psivv(1,-1,2,i,j),pmomvv(0,2,i,j),psi(1,-1,3),qbar(0,3),gvv32(1,-1))
          call curr6add(+1,psiv(1,-1,2,i),pmomv(0,2,i),psiv(1,-1,3,j),pmomv(0,3,j),gvv32(1,-1))
          call curr6add(+1,psi(1,-1,2),qbar(0,2),psivv(1,-1,3,i,j),pmomvv(0,3,i,j),gvv32(1,-1))
! <4|3>
          call curr6add(+1,psivv(1,-1,4,i,j),pmomvv(0,4,i,j),psi(1,-1,3),qbar(0,3),gvv34(1,-1))
          call curr6add(+1,psiv(1,-1,4,i),pmomv(0,4,i),psiv(1,-1,3,j),pmomv(0,3,j),gvv34(1,-1))
          call curr6add(+1,psi(1,-1,4),qbar(0,4),psivv(1,-1,3,i,j),pmomvv(0,3,i,j),gvv34(1,-1))
! <4|1>
          call curr6   (+1,psiv(1,-1,4,i),pmomv(0,4,i),psi(1,-1,1),qbar(0,1),gv14(1,-1,i))
          call curr6add(+1,psi(1,-1,4),qbar(0,4),psiv(1,-1,1,i),pmomv(0,1,i),gv14(1,-1,i))
! <2|3>
          call curr6   (+1,psiv(1,-1,2,i),pmomv(0,2,i),psi(1,-1,3),qbar(0,3),gv32(1,-1,i))
          call curr6add(+1,psi(1,-1,2),qbar(0,2),psiv(1,-1,3,i),pmomv(0,3,i),gv32(1,-1,i))
! <4|3>
          call curr6   (+1,psiv(1,-1,4,i),pmomv(0,4,i),psi(1,-1,3),qbar(0,3),gv34(1,-1,i))
          call curr6add(+1,psi(1,-1,4),qbar(0,4),psiv(1,-1,3,i),pmomv(0,3,i),gv34(1,-1,i))
        enddo
! <4|1>
        call curr6(+1,psi(1,-1,4),qbar(0,4),psi(1,-1,1),qbar(0,1),gnov14(1,-1))
! <2|3>
        call curr6(+1,psi(1,-1,2),qbar(0,2),psi(1,-1,3),qbar(0,3),gnov32(1,-1))
! <4|3>
        call curr6(+1,psi(1,-1,4),qbar(0,4),psi(1,-1,3),qbar(0,3),gnov34(1,-1))
! gluon propagator
        call propagate(+1,gvv14(1,-1),gmomvv14(0,-1))
        call propagate(+1,gvv32(1,-1),gmomvv32(0,-1))
        call propagate(+1,gvv34(1,-1),gmomvv34(0,-1))
        do i=3,4
          call propagate(+1,gv14(1,-1,i),gmomv14(0,-1,i))
          call propagate(+1,gv32(1,-1,i),gmomv32(0,-1,i))
          call propagate(+1,gv34(1,-1,i),gmomv34(0,-1,i))
        enddo
        call propagate(+1,gnov14(1,-1),gmomnov14(0,-1))
        call propagate(+1,gnov32(1,-1),gmomnov32(0,-1))
        call propagate(+1,gnov34(1,-1),gmomnov34(0,-1))

!
! Now build the other fermion line, attaching the photons and the W-current as necessary
!
! First add again the vector bosons to the <bra| and |ket>
! Order is important to avoid double counting !!!
!
! first all <bra| and |ket> with 0 A
!! 0 A via W current
        call bra2c(psi(1,-1,2),.true.,qbar(0,2),-1,
     &    wemnumom(0),w0(1,1),
     &    psinovw(1,-1,2),pmomnovw(0,-1,2)) !w0 to <2|
        call bra2c(psi(1,-1,4),.true.,qbar(0,4),-1,
     &    wemnumom(0),w0(1,1),
     &    psinovw(1,-1,4),pmomnovw(0,-1,4)) !w0 to <4|
        call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,
     &    wemnumom(0),w0(1,1),
     &    psinovw(1,-1,1),pmomnovw(0,-1,1)) !w0 to |1>
        call ket2c(psi(1,-1,3),.true.,qbar(0,3),-1,
     &    wemnumom(0),w0(1,1),
     &    psinovw(1,-1,3),pmomnovw(0,-1,3)) !w0 to |3>

! second all <bra| and |ket> with 1 A
!! 1 A via W current
        call bra2c(psi(1,-1,2),.true.,qbar(0,2),-1,
     &    wemnuga1mom(0),wpwa1(0,1),
     &    psivw(1,-1,2,3),pmomvw(0,-1,2,3)) !w3 to <2|
        call bra2c(psi(1,-1,2),.true.,qbar(0,2),-1,
     &    wemnuga2mom(0),wpwa2(0,1),
     &    psivw(1,-1,2,4),pmomvw(0,-1,2,4)) !w4 to <2|
        call bra2c(psi(1,-1,4),.true.,qbar(0,4),-1,
     &    wemnuga1mom(0),wpwa1(0,1),
     &    psivw(1,-1,4,3),pmomvw(0,-1,4,3)) !w3 to <4|
        call bra2c(psi(1,-1,4),.true.,qbar(0,4),-1,
     &    wemnuga2mom(0),wpwa2(0,1),
     &    psivw(1,-1,4,4),pmomvw(0,-1,4,4)) !w4 to <4|
        call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,
     &    wemnuga1mom(0),wpwa1(0,1),
     &    psivw(1,-1,1,3),pmomvw(0,-1,1,3)) !w3 to |1>
        call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,
     &    wemnuga2mom(0),wpwa2(0,1),
     &    psivw(1,-1,1,4),pmomvw(0,-1,1,4)) !w4 to |1>
        call ket2c(psi(1,-1,3),.true.,qbar(0,3),-1,
     &    wemnuga1mom(0),wpwa1(0,1),
     &    psivw(1,-1,3,3),pmomvw(0,-1,3,3)) !w3 to |3>
        call ket2c(psi(1,-1,3),.true.,qbar(0,3),-1,
     &    wemnuga2mom(0),wpwa2(0,1),
     &    psivw(1,-1,3,4),pmomvw(0,-1,3,4)) !w4 to |3>
!! 1 A attached before (0A)-W
        do i=3,4
          call bra2cadd(psiv(1,-1,2,i),.false.,pmomv(0,2,i),-1,
     &      wemnumom(0),w0(1,1),
     &      psivw(1,-1,2,i),pmomvw(0,-1,2,i),clr(iflavd,1,-1)) !w0 to <2|i
          call bra2cadd(psiv(1,-1,4,i),.false.,pmomv(0,4,i),-1,
     &      wemnumom(0),w0(1,1),
     &      psivw(1,-1,4,i),pmomvw(0,-1,4,i),clr(iflavd,1,-1)) !w0 to <4|i
          call ket2cadd(psiv(1,-1,1,i),.false.,pmomv(0,1,i),-1,
     &      wemnumom(0),w0(1,1),
     &      psivw(1,-1,1,i),pmomvw(0,-1,1,i),clr(iflavu,1,-1)) !w0 to i|1>
          call ket2cadd(psiv(1,-1,3,i),.false.,pmomv(0,3,i),-1,
     &      wemnumom(0),w0(1,1),
     &      psivw(1,-1,3,i),pmomvw(0,-1,3,i),clr(iflavu,1,-1)) !w0 to i|3>
        enddo !i
!! 1 A attached after (0A)-W
        call bra2cadd(psinovw(1,-1,2),.false.,pmomnovw(0,-1,2),-1,
     &    phot1fourvec(0),atau1(1,1),
     &    psivw(1,-1,2,3),pmomvw(0,-1,2,3),clr(iflavu,1,-1)) !3 to <2|w0
        call bra2cadd(psinovw(1,-1,2),.false.,pmomnovw(0,-1,2),-1,
     &    phot2fourvec(0),atau2(1,1),
     &    psivw(1,-1,2,4),pmomvw(0,-1,2,4),clr(iflavu,1,-1)) !4 to <2|w0
        call bra2cadd(psinovw(1,-1,4),.false.,pmomnovw(0,-1,4),-1,
     &    phot1fourvec(0),atau1(1,1),
     &    psivw(1,-1,4,3),pmomvw(0,-1,4,3),clr(iflavu,1,-1)) !3 to <4|w0
        call bra2cadd(psinovw(1,-1,4),.false.,pmomnovw(0,-1,4),-1,
     &    phot2fourvec(0),atau2(1,1),
     &    psivw(1,-1,4,4),pmomvw(0,-1,4,4),clr(iflavu,1,-1)) !4 to <4|w0
        call ket2cadd(psinovw(1,-1,1),.false.,pmomnovw(0,-1,1),-1,
     &    phot1fourvec(0),atau1(1,1),
     &    psivw(1,-1,1,3),pmomvw(0,-1,1,3),clr(iflavd,1,-1)) !3 to w0|1>
        call ket2cadd(psinovw(1,-1,1),.false.,pmomnovw(0,-1,1),-1,
     &    phot2fourvec(0),atau2(1,1),
     &    psivw(1,-1,1,4),pmomvw(0,-1,1,4),clr(iflavd,1,-1)) !4 to w0|1>
        call ket2cadd(psinovw(1,-1,3),.false.,pmomnovw(0,-1,3),-1,
     &    phot1fourvec(0),atau1(1,1),
     &    psivw(1,-1,3,3),pmomvw(0,-1,3,3),clr(iflavd,1,-1)) !3 to w0|3>
        call ket2cadd(psinovw(1,-1,3),.false.,pmomnovw(0,-1,3),-1,
     &    phot2fourvec(0),atau2(1,1),
     &    psivw(1,-1,3,4),pmomvw(0,-1,3,4),clr(iflavd,1,-1)) !4 to w0|3>

! third all <bra| and |ket> with 2 A
!! 2 A via W current
        call bra2c(psi(1,-1,2),.true.,qbar(0,2),-1,
     &    wemnugagamom(0),wpwaa(0,1),
     &    psivvw(1,-1,2),pmomvvw(0,-1,2)) !w34 to <2|
        call bra2c(psi(1,-1,4),.true.,qbar(0,4),-1,
     &    wemnugagamom(0),wpwaa(0,1),
     &    psivvw(1,-1,4),pmomvvw(0,-1,4)) !w34 to <4|
        call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,
     &    wemnugagamom(0),wpwaa(0,1),
     &    psivvw(1,-1,1),pmomvvw(0,-1,1)) !w34 to |1>
        call ket2c(psi(1,-1,3),.true.,qbar(0,3),-1,
     &    wemnugagamom(0),wpwaa(0,1),
     &    psivvw(1,-1,3),pmomvvw(0,-1,3)) !w34 to |3>
!! 1 A attached before (1A)-W
        call bra2cadd(psiv(1,-1,2,4),.false.,pmomv(0,2,4),-1,
     &    wemnuga1mom(0),wpwa1(0,1),
     &    psivvw(1,-1,2),pmomvvw(0,-1,2),clr(iflavd,1,-1)) !w3 to <2|4
        call bra2cadd(psiv(1,-1,2,3),.false.,pmomv(0,2,3),-1,
     &    wemnuga2mom(0),wpwa2(0,1),
     &    psivvw(1,-1,2),pmomvvw(0,-1,2),clr(iflavd,1,-1)) !w4 to <2|3
        call bra2cadd(psiv(1,-1,4,4),.false.,pmomv(0,4,4),-1,
     &    wemnuga1mom(0),wpwa1(0,1),
     &    psivvw(1,-1,4),pmomvvw(0,-1,4),clr(iflavd,1,-1)) !w3 to <4|4
        call bra2cadd(psiv(1,-1,4,3),.false.,pmomv(0,4,3),-1,
     &    wemnuga2mom(0),wpwa2(0,1),
     &    psivvw(1,-1,4),pmomvvw(0,-1,4),clr(iflavd,1,-1)) !w4 to <4|3
        call ket2cadd(psiv(1,-1,1,4),.false.,pmomv(0,1,4),-1,
     &    wemnuga1mom(0),wpwa1(0,1),
     &    psivvw(1,-1,1),pmomvvw(0,-1,1),clr(iflavu,1,-1)) !w3 to 4|1>
        call ket2cadd(psiv(1,-1,1,3),.false.,pmomv(0,1,3),-1,
     &    wemnuga2mom(0),wpwa2(0,1),
     &    psivvw(1,-1,1),pmomvvw(0,-1,1),clr(iflavu,1,-1)) !w4 to 3|1>
        call ket2cadd(psiv(1,-1,3,4),.false.,pmomv(0,3,4),-1,
     &    wemnuga1mom(0),wpwa1(0,1),
     &    psivvw(1,-1,3),pmomvvw(0,-1,3),clr(iflavu,1,-1)) !w3 to 4|3>
        call ket2cadd(psiv(1,-1,3,3),.false.,pmomv(0,3,3),-1,
     &    wemnuga2mom(0),wpwa2(0,1),
     &    psivvw(1,-1,3),pmomvvw(0,-1,3),clr(iflavu,1,-1)) !w4 to 3|3>
!! 1 A attached after (1A)-W
        call bra2cadd(psivw(1,-1,2,3),.false.,pmomvw(0,-1,2,3),-1,
     &    phot2fourvec(0),atau2(1,1),
     &    psivvw(1,-1,2),pmomvvw(0,-1,2),clr(iflavu,1,-1)) !4 to <2|w3
        call bra2cadd(psivw(1,-1,2,4),.false.,pmomvw(0,-1,2,4),-1,
     &    phot1fourvec(0),atau1(1,1),
     &    psivvw(1,-1,2),pmomvvw(0,-1,2),clr(iflavu,1,-1)) !3 to <2|w4
        call bra2cadd(psivw(1,-1,4,3),.false.,pmomvw(0,-1,4,3),-1,
     &    phot2fourvec(0),atau2(1,1),
     &    psivvw(1,-1,4),pmomvvw(0,-1,4),clr(iflavu,1,-1)) !4 to <4|w3
        call bra2cadd(psivw(1,-1,4,4),.false.,pmomvw(0,-1,4,4),-1,
     &    phot1fourvec(0),atau1(1,1),
     &    psivvw(1,-1,4),pmomvvw(0,-1,4),clr(iflavu,1,-1)) !3 to <4|w4
        call ket2cadd(psivw(1,-1,1,3),.false.,pmomvw(0,-1,1,3),-1,
     &    phot2fourvec(0),atau2(1,1),
     &    psivvw(1,-1,1),pmomvvw(0,-1,1),clr(iflavd,1,-1)) !4 to w3|1>
        call ket2cadd(psivw(1,-1,1,4),.false.,pmomvw(0,-1,1,4),-1,
     &    phot1fourvec(0),atau1(1,1),
     &    psivvw(1,-1,1),pmomvvw(0,-1,1),clr(iflavd,1,-1)) !3 to w4|1>
        call ket2cadd(psivw(1,-1,3,3),.false.,pmomvw(0,-1,3,3),-1,
     &    phot2fourvec(0),atau2(1,1),
     &    psivvw(1,-1,3),pmomvvw(0,-1,3),clr(iflavd,1,-1)) !4 to w3|3>
        call ket2cadd(psivw(1,-1,3,4),.false.,pmomvw(0,-1,3,4),-1,
     &    phot1fourvec(0),atau1(1,1),
     &    psivvw(1,-1,3),pmomvvw(0,-1,3),clr(iflavd,1,-1)) !3 to w4|3>
!! 2 A attached before (0A)-W
        do i=3,4
          j = 7-i
          call bra2cadd(psivv(1,-1,2,i,j),.false.,pmomvv(0,2,i,j),-1,
     &      wemnumom(0),w0(1,1),
     &      psivvw(1,-1,2),pmomvvw(0,-1,2),clr(iflavd,1,-1)**2) !w0 to <2|ij
          call bra2cadd(psivv(1,-1,4,i,j),.false.,pmomvv(0,4,i,j),-1,
     &      wemnumom(0),w0(1,1),
     &      psivvw(1,-1,4),pmomvvw(0,-1,4),clr(iflavd,1,-1)**2) !w0 to <4|ij
          call ket2cadd(psivv(1,-1,1,i,j),.false.,pmomvv(0,1,i,j),-1,
     &      wemnumom(0),w0(1,1),
     &      psivvw(1,-1,1),pmomvvw(0,-1,1),clr(iflavu,1,-1)**2) !w0 to ij|1>
          call ket2cadd(psivv(1,-1,3,i,j),.false.,pmomvv(0,3,i,j),-1,
     &      wemnumom(0),w0(1,1),
     &      psivvw(1,-1,3),pmomvvw(0,-1,3),clr(iflavu,1,-1)**2) !w0 to ij|3>
        enddo !i

!
! Now cycle through position of g and build currents
! 

! 2 A in current
        call curr6   (-1,psivvw(1,-1,2),pmomvvw(0,-1,2),
     &         psi(1,-1,1),qbar(0,1),jvvw(1,-1,1)) !<2|WAAg|1>
        call curr6   (-1,psivvw(1,-1,2),pmomvvw(0,-1,2),
     &         psi(1,-1,3),qbar(0,3),jvvw(1,-1,2)) !<2|WAAg|3>
        call curr6   (-1,psivvw(1,-1,4),pmomvvw(0,-1,4),
     &         psi(1,-1,1),qbar(0,1),jvvw(1,-1,3)) !<4|WAAg|1>
        call curr6add(-1,psi(1,-1,2),qbar(0,2),
     &         psivvw(1,-1,1),pmomvvw(0,-1,1),jvvw(1,-1,1)) !<2|gWAA|1>
        call curr6add(-1,psi(1,-1,2),qbar(0,2),
     &         psivvw(1,-1,3),pmomvvw(0,-1,3),jvvw(1,-1,2)) !<2|gWAA|3>
        call curr6add(-1,psi(1,-1,4),qbar(0,4),
     &         psivvw(1,-1,1),pmomvvw(0,-1,1),jvvw(1,-1,3)) !<4|gWAA|1>
        do i=3,4
          j = 7-i
          call curr6addm(-1,psivw(1,-1,2,i),pmomvw(0,-1,2,i),
     &           psiv(1,-1,1,j),pmomv(0,1,j),jvvw(1,-1,1),
     &           clr(iflavu,1,-1),clr(iflavu,1,+1)) !<2|WAgA|1>
          call curr6addm(-1,psivw(1,-1,2,i),pmomvw(0,-1,2,i),
     &           psiv(1,-1,3,j),pmomv(0,3,j),jvvw(1,-1,2),
     &           clr(iflavu,1,-1),clr(iflavu,1,+1)) !<2|WAgA|3>
          call curr6addm(-1,psivw(1,-1,4,i),pmomvw(0,-1,4,i),
     &           psiv(1,-1,1,j),pmomv(0,1,j),jvvw(1,-1,3),
     &           clr(iflavu,1,-1),clr(iflavu,1,+1)) !<4|WAgA|1>
          call curr6addm(-1,psiv(1,-1,2,i),pmomv(0,2,i),
     &           psivw(1,-1,1,j),pmomvw(0,-1,1,j),jvvw(1,-1,1),
     &           clr(iflavd,1,-1),clr(iflavd,1,+1)) !<2|AgWA|1>
          call curr6addm(-1,psiv(1,-1,2,i),pmomv(0,2,i),
     &           psivw(1,-1,3,j),pmomvw(0,-1,3,j),jvvw(1,-1,2),
     &           clr(iflavd,1,-1),clr(iflavd,1,+1)) !<2|AgWA|3>
          call curr6addm(-1,psiv(1,-1,4,i),pmomv(0,4,i),
     &           psivw(1,-1,1,j),pmomvw(0,-1,1,j),jvvw(1,-1,3),
     &           clr(iflavd,1,-1),clr(iflavd,1,+1)) !<4|AgWA|1>
          call curr6addm(-1,psivv(1,-1,2,i,j),pmomvv(0,2,i,j),
     &           psinovw(1,-1,1),pmomnovw(0,-1,1),jvvw(1,-1,1),
     &           clr(iflavd,1,-1)**2,clr(iflavd,1,+1)**2) !<2|AAgW|1>
          call curr6addm(-1,psivv(1,-1,2,i,j),pmomvv(0,2,i,j),
     &           psinovw(1,-1,3),pmomnovw(0,-1,3),jvvw(1,-1,2),
     &           clr(iflavd,1,-1)**2,clr(iflavd,1,+1)**2) !<2|AAgW|3>
          call curr6addm(-1,psivv(1,-1,4,i,j),pmomvv(0,4,i,j),
     &           psinovw(1,-1,1),pmomnovw(0,-1,1),jvvw(1,-1,3),
     &           clr(iflavd,1,-1)**2,clr(iflavd,1,+1)**2) !<4|AAgW|1>
          call curr6addm(-1,psinovw(1,-1,2),pmomnovw(0,-1,2),
     &           psivv(1,-1,1,i,j),pmomvv(0,1,i,j),jvvw(1,-1,1),
     &           clr(iflavu,1,-1)**2,clr(iflavu,1,+1)**2) !<2|WgAA|1>
          call curr6addm(-1,psinovw(1,-1,2),pmomnovw(0,-1,2),
     &           psivv(1,-1,3,i,j),pmomvv(0,3,i,j),jvvw(1,-1,2),
     &           clr(iflavu,1,-1)**2,clr(iflavu,1,+1)**2) !<2|WgAA|3>
          call curr6addm(-1,psinovw(1,-1,4),pmomnovw(0,-1,4),
     &           psivv(1,-1,1,i,j),pmomvv(0,1,i,j),jvvw(1,-1,3),
     &           clr(iflavu,1,-1)**2,clr(iflavu,1,+1)**2) !<4|WgAA|1>
        enddo
        do chi=-1,1,2
          temp = dotcc(jvvw(1,-1,1),gnov34(1,chi))
          uuud(1,-1,chi) = uuud(1,-1,chi) + temp
          uddd(1,-1,chi) = uddd(1,-1,chi) + temp
          temp = dotcc(jvvw(1,-1,2),gnov14(1,chi))
          uuud(2,-1,chi) = uuud(2,-1,chi) + temp
          temp = dotcc(jvvw(1,-1,3),gnov32(1,chi))
          uddd(2,-1,chi) = uddd(2,-1,chi) + temp
        enddo !chi
! 1 A in current 
        do i=3,4
          j = 7-i
          call curr6   (-1,psivw(1,-1,2,i),pmomvw(0,-1,2,i),
     &           psi(1,-1,1),qbar(0,1),jvw(1,-1,1,i)) !<2|WAg|1>
          call curr6   (-1,psivw(1,-1,2,i),pmomvw(0,-1,2,i),
     &           psi(1,-1,3),qbar(0,3),jvw(1,-1,2,i)) !<2|WAg|3>
          call curr6   (-1,psivw(1,-1,4,i),pmomvw(0,-1,4,i),
     &           psi(1,-1,1),qbar(0,1),jvw(1,-1,3,i)) !<4|WAg|1>
          call curr6add(-1,psi(1,-1,2),qbar(0,2),
     &           psivw(1,-1,1,i),pmomvw(0,-1,1,i),jvw(1,-1,1,i)) !<2|gWA|1>
          call curr6add(-1,psi(1,-1,2),qbar(0,2),
     &           psivw(1,-1,3,i),pmomvw(0,-1,3,i),jvw(1,-1,2,i)) !<2|gWA|3>
          call curr6add(-1,psi(1,-1,4),qbar(0,4),
     &           psivw(1,-1,1,i),pmomvw(0,-1,1,i),jvw(1,-1,3,i)) !<4|gWA|1>
          call curr6addm(-1,psiv(1,-1,2,i),pmomv(0,2,i),
     &           psinovw(1,-1,1),pmomnovw(0,-1,1),jvw(1,-1,1,i),
     &           clr(iflavd,1,-1),clr(iflavd,1,+1)) !<2|AgW|1>
          call curr6addm(-1,psiv(1,-1,2,i),pmomv(0,2,i),
     &           psinovw(1,-1,3),pmomnovw(0,-1,3),jvw(1,-1,2,i),
     &           clr(iflavd,1,-1),clr(iflavd,1,+1)) !<2|AgW|3>
          call curr6addm(-1,psiv(1,-1,4,i),pmomv(0,4,i),
     &           psinovw(1,-1,1),pmomnovw(0,-1,1),jvw(1,-1,3,i),
     &           clr(iflavd,1,-1),clr(iflavd,1,+1)) !<4|AgW|1>
          call curr6addm(-1,psinovw(1,-1,2),pmomnovw(0,-1,2),
     &           psiv(1,-1,1,i),pmomv(0,1,i),jvw(1,-1,1,i),
     &           clr(iflavu,1,-1),clr(iflavu,1,+1)) !<2|WgA|1>
          call curr6addm(-1,psinovw(1,-1,2),pmomnovw(0,-1,2),
     &           psiv(1,-1,3,i),pmomv(0,3,i),jvw(1,-1,2,i),
     &           clr(iflavu,1,-1),clr(iflavu,1,+1)) !<2|WgA|3>
          call curr6addm(-1,psinovw(1,-1,4),pmomnovw(0,-1,4),
     &           psiv(1,-1,1,i),pmomv(0,1,i),jvw(1,-1,3,i),
     &           clr(iflavu,1,-1),clr(iflavu,1,+1)) !<4|WgA|1>
          do chi=-1,1,2
            temp = dotcc(jvw(1,-1,1,i),gv34(1,chi,j))
            uuud(1,-1,chi) = uuud(1,-1,chi) + temp*clr(iflavu,1,-1)
            uddd(1,-1,chi) = uddd(1,-1,chi) + temp*clr(iflavd,1,-1)
            temp = dotcc(jvw(1,-1,2,i),gv14(1,chi,j))
            uuud(2,-1,chi) = uuud(2,-1,chi) + temp*clr(iflavu,1,-1)
            temp = dotcc(jvw(1,-1,3,i),gv32(1,chi,j))
            uddd(2,-1,chi) = uddd(2,-1,chi) + temp*clr(iflavd,1,-1)
          enddo !chi
        enddo
! 0 A in current
        call curr6   (-1,psinovw(1,-1,2),pmomnovw(0,-1,2),
     &         psi(1,-1,1),qbar(0,1),jnovw(1,-1,1))
        call curr6   (-1,psinovw(1,-1,2),pmomnovw(0,-1,2),
     &         psi(1,-1,3),qbar(0,3),jnovw(1,-1,2))
        call curr6   (-1,psinovw(1,-1,4),pmomnovw(0,-1,4),
     &         psi(1,-1,1),qbar(0,1),jnovw(1,-1,3))
        call curr6add(-1,psi(1,-1,2),qbar(0,2),
     &         psinovw(1,-1,1),pmomnovw(0,-1,1),jnovw(1,-1,1))
        call curr6add(-1,psi(1,-1,2),qbar(0,2),
     &         psinovw(1,-1,3),pmomnovw(0,-1,3),jnovw(1,-1,2))
        call curr6add(-1,psi(1,-1,4),qbar(0,4),
     &         psinovw(1,-1,1),pmomnovw(0,-1,1),jnovw(1,-1,3))
        do chi=-1,1,2
          temp = dotcc(jnovw(1,-1,1),gvv34(1,chi))
          uuud(1,-1,chi) = uuud(1,-1,chi) + temp*clr(iflavu,1,-1)**2
          uddd(1,-1,chi) = uddd(1,-1,chi) + temp*clr(iflavd,1,-1)**2
          temp = dotcc(jnovw(1,-1,2),gvv14(1,chi))
          uuud(2,-1,chi) = uuud(2,-1,chi) + temp*clr(iflavu,1,-1)**2
          temp = dotcc(jnovw(1,-1,3),gvv32(1,chi))
          uddd(2,-1,chi) = uddd(2,-1,chi) + temp*clr(iflavd,1,-1)**2
        enddo !chi

#ifdef WITH_MADGRAPH
!
! Compare with MadGraph
!
      if (ldebug) then
         print*,"MadGraph check qqWAAqq"
         print*,"Type: ", selc
         do chi=-1,1,2
         do i=1,2
c          uduu
           print*,"uu -> du WAA and crossed"
           print*,"mat(",i,",-1,",chi,") = ",uuud(i,-1,chi)
           temp = dreal(uuud(i,-1,chi))**2+dimag(uuud(i,-1,chi))**2
           print*,"mg2(",i,",-1,",chi,") = ",
     &       -tampr(i,-1,chi,1)
           temp1 = 
     &        dreal(tampr(i,-1,chi,1))**2
     &       +dimag(tampr(i,-1,chi,1))**2
           print*,"temp/temp1 = ",temp/temp1
c          uddd
           print*,"ud -> dd WAA and crossed"
           print*,"mat(",i,",-1,",chi,") = ",uddd(i,-1,chi)
           temp = dreal(uddd(i,-1,chi))**2+dimag(uddd(i,-1,chi))**2
           print*,"mg2(",i,",-1,",chi,") = ",
     &       -tampr(i,-1,chi,2)
           temp1 = 
     &        dreal(tampr(i,-1,chi,2))**2
     &       +dimag(tampr(i,-1,chi,2))**2
           print*,"temp/temp1 = ",temp/temp1
         enddo
         enddo
         read(*,*)
      endif
#endif

!
! Sum everything up
!

      do chi=-1,1,2
        m2su(1)=m2su(1)
     &           +2d0*abs(uuud(1,-1,chi))**2
     &           +2d0*abs(uuud(2,-1,chi))**2
     &           +4d0/3d0*dreal(uuud(1,-1,chi)*dconjg(uuud(2,-1,chi)))
        m2su(2)=m2su(2)
     &           +2d0*abs(uuud(2,-1,chi))**2
        m2su(3)=m2su(3)
     &           +2d0*abs(uuud(1,-1,chi))**2
        m2sd(1)=m2sd(1)
     &           +2d0*abs(uddd(1,-1,chi))**2
     &           +2d0*abs(uddd(2,-1,chi))**2
     &           +4d0/3d0*dreal(uddd(1,-1,chi)*dconjg(uddd(2,-1,chi)))
        m2sd(2)=m2sd(2)
     &           +2d0*abs(uddd(2,-1,chi))**2
        m2sd(3)=m2sd(3)
     &           +2d0*abs(uddd(1,-1,chi))**2
      enddo

!
! Multiply with remaining coupling constants
!
      do i=1,3
        m2su(i) = m2su(i) * fpials(2)**2 * clr(1,3,-1)**2
        m2sd(i) = m2sd(i) * fpials(2)**2 * clr(1,3,-1)**2
      enddo

      end


