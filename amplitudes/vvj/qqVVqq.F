      subroutine qqVVqq(p,v,selc,m2s)

      implicit none

      include 'VBFNLO/utilities/global.inc'
      include 'VBFNLO/utilities/lha.inc'
      include 'VBFNLO/utilities/scales.inc'
      include 'VBFNLO/utilities/coupl.inc'
      include 'tensor.inc'

      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)

      integer j,i,k,mu,isig,isig2
      integer sign(4)

      double precision p(0:3,max_p), v(0:3,max_v), m2s(7,2,2)
      double precision pbar(0:3,4), qbar(0:3,4)
     
      double complex psi(2,-1:1,4),psiv(2,-1:1,4,4),psivv(2,-1:1,4,4,4)
      double complex gvv12(1:6,-1:1,3:4), gvv32(1:6,-1:1,3:4),
     &               gv12(1:6,-1:1,3:4), gv32(1:6,-1:1,3:4),
     &               gx12(1:6,-1:1), gx32(1:6,-1:1),
     &               jx14(1:6,-1:1), jx34(1:6,-1:1), 
     &               jv14(1:6,-1:1,3:4), jv34(1:6,-1:1,3:4),
     &               jvv14(1:6,-1:1,3:4), jvv34(1:6,-1:1,3:4)
      double precision pmomv(0:4,4,4),pmomvv(0:4,4,4,4)
      double precision gmomvv12(0:4,-1:1,3:4), gmomvv32(0:4,-1:1,3:4),
     &               gmomv12(0:4,-1:1,3:4), gmomv32(0:4,-1:1,3:4),
     &               gmomx12(0:4,-1:1), gmomx32(0:4,-1:1)

      double complex mat(4,3:4,3:4,-1:1,-1:1)

      double complex im
      logical ldebug
      parameter (ldebug=.false.)
      parameter (im = (0d0,1d0))

      double complex dotcc
      external dotcc
    
      DOUBLE PRECISION fpials(2:3), fpi
      PARAMETER (fpi=4d0*pi)

      double complex temp,temp1
      integer selc

      double precision pall(0:3,8)
      double precision mg2uures(3), mg2ddres(3)
      double precision mg2udres(3), mg2dures(3)
      integer fsign(8)

      do i=1,7
      do j=1,2
      do k=1,2
        m2s(i,j,k)=0d0
      enddo
      enddo
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
!  q  q  > e+ ve mu- vm~ q  q
!
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
!  q  q~ > e+ ve mu- vm~ q  q~
!
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
!  q~ q  > e+ ve mu- vm~ q  q~ 
!
      sign(1)=1
      sign(2)=-1
      sign(3)=-1
      sign(4)=1
      do mu=0,3
      pbar(mu,1)=p(mu,2)   ! 1-> ket  
      pbar(mu,2)=p(mu,1)   ! 2-> bra
      pbar(mu,3)=p(mu,4)   ! 3-> ket
      pbar(mu,4)=p(mu,3)   ! 4-> bra
      enddo
!
!-----------------------------------------------------
      Case (4)  
!  q~ q~ > e+ ve mu- vm~ q~ q~
!
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
        fsign(3)=-1
        fsign(4)=+1
        fsign(5)=+1
        fsign(6)=-1
        do mu=0,3
          pall(mu,1) = pbar(mu,1)
          pall(mu,2) = pbar(mu,2)
          pall(mu,7) = pbar(mu,4)
          pall(mu,8) = pbar(mu,3)
          pall(mu,3) = v(mu,2)
          pall(mu,4) = v(mu,1)
          pall(mu,5) = v(mu,3)
          pall(mu,6) = v(mu,4)
        enddo
        call SUU_WWUU(pall,fsign,mg2uures(1),1)
        call SUU_WWUU(pall,fsign,mg2uures(2),2)
        call SUU_WWUU(pall,fsign,mg2uures(3),3)
        call SDD_WWDD(pall,fsign,mg2ddres(1),1)
        call SDD_WWDD(pall,fsign,mg2ddres(2),2)
        call SDD_WWDD(pall,fsign,mg2ddres(3),3)
        call SUD_WWUD(pall,fsign,mg2udres(1),1)
        call SUD_WWUD(pall,fsign,mg2udres(2),2)
        call SUD_WWUD(pall,fsign,mg2udres(3),3)
        call SDU_WWDU(pall,fsign,mg2dures(1),1)
        call SDU_WWDU(pall,fsign,mg2dures(2),2)
        call SDU_WWDU(pall,fsign,mg2dures(3),3)
      endif
#endif

!-----------------------------------------------------
!-----------------------------------------------------

      call psi0m(4,pbar(0,1),sign(1),psi)
      
      do i=1,4
      do j=3,4
      do k=3,4
      do isig=-1,1,2
      do isig2=-1,1,2
        mat(i,j,k,isig,isig2) = 0d0
      enddo
      enddo
      enddo
      enddo
      enddo

!
!  First construct the fermion line which has everything attached
!  (1-2 or 3-2)
!
!  Add particles to psi
!
!  in psiv() etc: a=1, z=2, w+=3, w-=4
!  order in arguments is as in math. expression <bra|V...   ...V|ket>
!

        do isig=-1,1,2
          do i=1,3,2
            j=i+1
! first V to <bra| or |ket>
            call bra2c(psi(1,isig,j),.true.,qbar(0,j),isig,
     &        qww(0,1),aww(0,1),psiv(1,isig,j,1),pmomv(0,j,1))
            call bra2c(psi(1,isig,j),.true.,qbar(0,j),isig,
     &        qww(0,1),zww(0,1),psiv(1,isig,j,2),pmomv(0,j,2))
            call ket2c(psi(1,isig,i),.true.,qbar(0,i),isig,
     &        qww(0,1),aww(0,1),psiv(1,isig,i,1),pmomv(0,i,1))
            call ket2c(psi(1,isig,i),.true.,qbar(0,i),isig,
     &        qww(0,1),zww(0,1),psiv(1,isig,i,2),pmomv(0,i,2))
            if (isig.eq.-1) then ! W only to left
              call bra2c(psi(1,isig,j),.true.,qbar(0,j),isig,
     &          qp(0,1),wp(1,1),psiv(1,isig,j,3),pmomv(0,j,3))
              call bra2c(psi(1,isig,j),.true.,qbar(0,j),isig,
     &          qm(0,1),wm(1,1),psiv(1,isig,j,4),pmomv(0,j,4))
              call ket2c(psi(1,isig,i),.true.,qbar(0,i),isig,
     &          qp(0,1),wp(1,1),psiv(1,isig,i,3),pmomv(0,i,3))
              call ket2c(psi(1,isig,i),.true.,qbar(0,i),isig,
     &          qm(0,1),wm(1,1),psiv(1,isig,i,4),pmomv(0,i,4))

! second V to <bra|V or V|ket> ! only W's
              call bra2c(psiv(1,isig,j,4),.false.,pmomv(0,j,4),isig,
     &          qp(0,1),wp(1,1),psivv(1,isig,j,4,3),pmomvv(0,j,4,3))
              call bra2c(psiv(1,isig,j,3),.false.,pmomv(0,j,3),isig,
     &          qm(0,1),wm(1,1),psivv(1,isig,j,3,4),pmomvv(0,j,3,4))
              call ket2c(psiv(1,isig,i,4),.false.,pmomv(0,i,4),isig,
     &          qp(0,1),wp(1,1),psivv(1,isig,i,3,4),pmomvv(0,i,3,4))
              call ket2c(psiv(1,isig,i,3),.false.,pmomv(0,i,3),isig,
     &          qm(0,1),wm(1,1),psivv(1,isig,i,4,3),pmomvv(0,i,4,3))
            endif

          enddo
        enddo

!
! Now cycle through position of gluon and build currents
! 
! initialize g*12 and g*32 
        do isig=-1,1,2
        do i=3,4
        do mu=1,4
          gvv12(mu,isig,i) = 0d0
          gvv32(mu,isig,i) = 0d0
          gv12(mu,isig,i) = 0d0
          gv32(mu,isig,i) = 0d0
        enddo
          gvv12(5,isig,i) = dcmplx(
     &      qbar(0,1)-qbar(0,2)-qp(0,1)-qm(0,1),
     &      qbar(3,1)-qbar(3,2)-qp(3,1)-qm(3,1))
          gvv12(6,isig,i) = dcmplx(
     &      qbar(1,1)-qbar(1,2)-qp(1,1)-qm(1,1),
     &      qbar(2,1)-qbar(2,2)-qp(2,1)-qm(2,1))
          gvv32(5,isig,i) = dcmplx(
     &      qbar(0,3)-qbar(0,2)-qp(0,1)-qm(0,1),
     &      qbar(3,3)-qbar(3,2)-qp(3,1)-qm(3,1))
          gvv32(6,isig,i) = dcmplx(
     &      qbar(1,3)-qbar(1,2)-qp(1,1)-qm(1,1),
     &      qbar(2,3)-qbar(2,2)-qp(2,1)-qm(2,1))
        enddo
          gv12(5,isig,3) = dcmplx(
     &      qbar(0,1)-qbar(0,2)-qm(0,1),
     &      qbar(3,1)-qbar(3,2)-qm(3,1))
          gv12(6,isig,3) = dcmplx(
     &      qbar(1,1)-qbar(1,2)-qm(1,1),
     &      qbar(2,1)-qbar(2,2)-qm(2,1))
          gv12(5,isig,4) = dcmplx(
     &      qbar(0,1)-qbar(0,2)-qp(0,1),
     &      qbar(3,1)-qbar(3,2)-qp(3,1))
          gv12(6,isig,4) = dcmplx(
     &      qbar(1,1)-qbar(1,2)-qp(1,1),
     &      qbar(2,1)-qbar(2,2)-qp(2,1))
          gv32(5,isig,3) = dcmplx(
     &      qbar(0,3)-qbar(0,2)-qm(0,1),
     &      qbar(3,3)-qbar(3,2)-qm(3,1))
          gv32(6,isig,3) = dcmplx(
     &      qbar(1,3)-qbar(1,2)-qm(1,1),
     &      qbar(2,3)-qbar(2,2)-qm(2,1))
          gv32(5,isig,4) = dcmplx(
     &      qbar(0,3)-qbar(0,2)-qp(0,1),
     &      qbar(3,3)-qbar(3,2)-qp(3,1))
          gv32(6,isig,4) = dcmplx(
     &      qbar(1,3)-qbar(1,2)-qp(1,1),
     &      qbar(2,3)-qbar(2,2)-qp(2,1))
        enddo

! add currents 
c both W
cc u quark
ccc <2|1>
cccc W
        call curr6addm(-1,psivv(1,-1,2,4,3),pmomvv(0,2,4,3),psi(1,-1,1),qbar(0,1),
     &    gvv12(1,-1,3),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psiv(1,-1,2,4),pmomv(0,2,4),psiv(1,-1,1,3),pmomv(0,1,3),
     &    gvv12(1,-1,3),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psi(1,-1,2),qbar(0,2),psivv(1,-1,1,4,3),pmomvv(0,1,4,3),
     &    gvv12(1,-1,3),clr(3,3,-1)**2,clr(3,3,+1)**2)
cccc A
        call curr6addm(+1,psiv(1,-1,2,1),pmomv(0,2,1),psi(1,-1,1),qbar(0,1),
     &    gvv12(1,-1,3),-clr(3,1,-1),-clr(3,1,+1))
        call curr6addm(+1,psi(1,-1,2),qbar(0,2),psiv(1,-1,1,1),pmomv(0,1,1),
     &    gvv12(1,-1,3),-clr(3,1,-1),-clr(3,1,+1))
cccc Z
        call curr6addm(+1,psiv(1,-1,2,2),pmomv(0,2,2),psi(1,-1,1),qbar(0,1),
     &    gvv12(1,-1,3),-clr(3,2,-1),-clr(3,2,+1))
        call curr6addm(+1,psi(1,-1,2),qbar(0,2),psiv(1,-1,1,2),pmomv(0,1,2),
     &    gvv12(1,-1,3),-clr(3,2,-1),-clr(3,2,+1))
ccc <2|3>
cccc W
        call curr6addm(-1,psivv(1,-1,2,4,3),pmomvv(0,2,4,3),psi(1,-1,3),qbar(0,3),
     &    gvv32(1,-1,3),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psiv(1,-1,2,4),pmomv(0,2,4),psiv(1,-1,3,3),pmomv(0,3,3),
     &    gvv32(1,-1,3),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psi(1,-1,2),qbar(0,2),psivv(1,-1,3,4,3),pmomvv(0,3,4,3),
     &    gvv32(1,-1,3),clr(3,3,-1)**2,clr(3,3,+1)**2)
cccc A
        call curr6addm(+1,psiv(1,-1,2,1),pmomv(0,2,1),psi(1,-1,3),qbar(0,3),
     &    gvv32(1,-1,3),-clr(3,1,-1),-clr(3,1,+1))
        call curr6addm(+1,psi(1,-1,2),qbar(0,2),psiv(1,-1,3,1),pmomv(0,3,1),
     &    gvv32(1,-1,3),-clr(3,1,-1),-clr(3,1,+1))
cccc Z
        call curr6addm(+1,psiv(1,-1,2,2),pmomv(0,2,2),psi(1,-1,3),qbar(0,3),
     &    gvv32(1,-1,3),-clr(3,2,-1),-clr(3,2,+1))
        call curr6addm(+1,psi(1,-1,2),qbar(0,2),psiv(1,-1,3,2),pmomv(0,3,2),
     &    gvv32(1,-1,3),-clr(3,2,-1),-clr(3,2,+1))
cc d quark
ccc <2|1>
cccc W
        call curr6addm(-1,psivv(1,-1,2,3,4),pmomvv(0,2,3,4),psi(1,-1,1),qbar(0,1),
     &    gvv12(1,-1,4),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psiv(1,-1,2,3),pmomv(0,2,3),psiv(1,-1,1,4),pmomv(0,1,4),
     &    gvv12(1,-1,4),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psi(1,-1,2),qbar(0,2),psivv(1,-1,1,3,4),pmomvv(0,1,3,4),
     &    gvv12(1,-1,4),clr(3,3,-1)**2,clr(3,3,+1)**2)
cccc A
        call curr6addm(+1,psiv(1,-1,2,1),pmomv(0,2,1),psi(1,-1,1),qbar(0,1),
     &    gvv12(1,-1,4),-clr(4,1,-1),-clr(4,1,+1))
        call curr6addm(+1,psi(1,-1,2),qbar(0,2),psiv(1,-1,1,1),pmomv(0,1,1),
     &    gvv12(1,-1,4),-clr(4,1,-1),-clr(4,1,+1))
cccc Z
        call curr6addm(+1,psiv(1,-1,2,2),pmomv(0,2,2),psi(1,-1,1),qbar(0,1),
     &    gvv12(1,-1,4),-clr(4,2,-1),-clr(4,2,+1))
        call curr6addm(+1,psi(1,-1,2),qbar(0,2),psiv(1,-1,1,2),pmomv(0,1,2),
     &    gvv12(1,-1,4),-clr(4,2,-1),-clr(4,2,+1))
ccc <2|3>
cccc W
        call curr6addm(-1,psivv(1,-1,2,3,4),pmomvv(0,2,3,4),psi(1,-1,3),qbar(0,3),
     &    gvv32(1,-1,4),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psiv(1,-1,2,3),pmomv(0,2,3),psiv(1,-1,3,4),pmomv(0,3,4),
     &    gvv32(1,-1,4),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psi(1,-1,2),qbar(0,2),psivv(1,-1,3,3,4),pmomvv(0,3,3,4),
     &    gvv32(1,-1,4),clr(3,3,-1)**2,clr(3,3,+1)**2)
cccc A
        call curr6addm(+1,psiv(1,-1,2,1),pmomv(0,2,1),psi(1,-1,3),qbar(0,3),
     &    gvv32(1,-1,4),-clr(4,1,-1),-clr(4,1,+1))
        call curr6addm(+1,psi(1,-1,2),qbar(0,2),psiv(1,-1,3,1),pmomv(0,3,1),
     &    gvv32(1,-1,4),-clr(4,1,-1),-clr(4,1,+1))
cccc Z
        call curr6addm(+1,psiv(1,-1,2,2),pmomv(0,2,2),psi(1,-1,3),qbar(0,3),
     &    gvv32(1,-1,4),-clr(4,2,-1),-clr(4,2,+1))
        call curr6addm(+1,psi(1,-1,2),qbar(0,2),psiv(1,-1,3,2),pmomv(0,3,2),
     &    gvv32(1,-1,4),-clr(4,2,-1),-clr(4,2,+1))

c single W attached (fermion 2 designates u/d)
cc <2|1>
ccc u quark
        call curr6addm(-1,psiv(1,-1,2,4),pmomv(0,2,4),psi(1,-1,1),qbar(0,1),
     &    gv12(1,-1,3),clr(3,3,-1),clr(3,3,+1))
        call curr6addm(-1,psi(1,-1,2),qbar(0,2),psiv(1,-1,1,4),pmomv(0,1,4),
     &    gv12(1,-1,3),clr(3,3,-1),clr(3,3,+1))
ccc d quark
        call curr6addm(-1,psiv(1,-1,2,3),pmomv(0,2,3),psi(1,-1,1),qbar(0,1),
     &    gv12(1,-1,4),clr(3,3,-1),clr(3,3,+1))
        call curr6addm(-1,psi(1,-1,2),qbar(0,2),psiv(1,-1,1,3),pmomv(0,1,3),
     &    gv12(1,-1,4),clr(3,3,-1),clr(3,3,+1))
cc <2|3>
ccc u quark
        call curr6addm(-1,psiv(1,-1,2,4),pmomv(0,2,4),psi(1,-1,3),qbar(0,3),
     &    gv32(1,-1,3),clr(3,3,-1),clr(3,3,+1))
        call curr6addm(-1,psi(1,-1,2),qbar(0,2),psiv(1,-1,3,4),pmomv(0,3,4),
     &    gv32(1,-1,3),clr(3,3,-1),clr(3,3,+1))
ccc d quark
        call curr6addm(-1,psiv(1,-1,2,3),pmomv(0,2,3),psi(1,-1,3),qbar(0,3),
     &    gv32(1,-1,4),clr(3,3,-1),clr(3,3,+1))
        call curr6addm(-1,psi(1,-1,2),qbar(0,2),psiv(1,-1,3,3),pmomv(0,3,3),
     &    gv32(1,-1,4),clr(3,3,-1),clr(3,3,+1))

c zero W
        call curr6(+1,psi(1,-1,2),qbar(0,2),psi(1,-1,1),qbar(0,1),gx12(1,-1))
        call curr6(+1,psi(1,-1,2),qbar(0,2),psi(1,-1,3),qbar(0,3),gx32(1,-1))

! gluon propagator
        do i=3,4
          call propagate(+1,gvv12(1,-1,i),gmomvv12(0,-1,i))
          call propagate(+1,gvv32(1,-1,i),gmomvv32(0,-1,i))
          call propagate(+1,gv12(1,-1,i),gmomv12(0,-1,i))
          call propagate(+1,gv32(1,-1,i),gmomv32(0,-1,i))
        enddo
        call propagate(+1,gx12(1,-1),gmomx12(0,-1))
        call propagate(+1,gx32(1,-1),gmomx32(0,-1))

!
! Now build the other fermion line, attaching the remaining particles as necessary
!
! initialize j*14 and j*34 
        do isig=-1,1,2
        do i=3,4
        do mu=1,4
          jvv14(mu,isig,i) = 0d0
          jvv34(mu,isig,i) = 0d0
          jv14(mu,isig,i) = 0d0
          jv34(mu,isig,i) = 0d0
        enddo
          jvv14(5,isig,i) = dcmplx(
     &      qbar(0,1)-qbar(0,4)-qp(0,1)-qm(0,1),
     &      qbar(3,1)-qbar(3,4)-qp(3,1)-qm(3,1))
          jvv14(6,isig,i) = dcmplx(
     &      qbar(1,1)-qbar(1,4)-qp(1,1)-qm(1,1),
     &      qbar(2,1)-qbar(2,4)-qp(2,1)-qm(2,1))
          jvv34(5,isig,i) = dcmplx(
     &      qbar(0,3)-qbar(0,4)-qp(0,1)-qm(0,1),
     &      qbar(3,3)-qbar(3,4)-qp(3,1)-qm(3,1))
          jvv34(6,isig,i) = dcmplx(
     &      qbar(1,3)-qbar(1,4)-qp(1,1)-qm(1,1),
     &      qbar(2,3)-qbar(2,4)-qp(2,1)-qm(2,1))
        enddo
          jv14(5,isig,3) = dcmplx(
     &      qbar(0,1)-qbar(0,4)-qm(0,1),
     &      qbar(3,1)-qbar(3,4)-qm(3,1))
          jv14(6,isig,3) = dcmplx(
     &      qbar(1,1)-qbar(1,4)-qm(1,1),
     &      qbar(2,1)-qbar(2,4)-qm(2,1))
          jv14(5,isig,4) = dcmplx(
     &      qbar(0,1)-qbar(0,4)-qp(0,1),
     &      qbar(3,1)-qbar(3,4)-qp(3,1))
          jv14(6,isig,4) = dcmplx(
     &      qbar(1,1)-qbar(1,4)-qp(1,1),
     &      qbar(2,1)-qbar(2,4)-qp(2,1))
          jv34(5,isig,3) = dcmplx(
     &      qbar(0,3)-qbar(0,4)-qm(0,1),
     &      qbar(3,3)-qbar(3,4)-qm(3,1))
          jv34(6,isig,3) = dcmplx(
     &      qbar(1,3)-qbar(1,4)-qm(1,1),
     &      qbar(2,3)-qbar(2,4)-qm(2,1))
          jv34(5,isig,4) = dcmplx(
     &      qbar(0,3)-qbar(0,4)-qp(0,1),
     &      qbar(3,3)-qbar(3,4)-qp(3,1))
          jv34(6,isig,4) = dcmplx(
     &      qbar(1,3)-qbar(1,4)-qp(1,1),
     &      qbar(2,3)-qbar(2,4)-qp(2,1))
        enddo

! add currents 
c both W
cc u quark
ccc <4|1>
cccc W
        call curr6addm(-1,psivv(1,-1,4,4,3),pmomvv(0,4,4,3),psi(1,-1,1),qbar(0,1),
     &    jvv14(1,-1,3),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psiv(1,-1,4,4),pmomv(0,4,4),psiv(1,-1,1,3),pmomv(0,1,3),
     &    jvv14(1,-1,3),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psi(1,-1,4),qbar(0,4),psivv(1,-1,1,4,3),pmomvv(0,1,4,3),
     &    jvv14(1,-1,3),clr(3,3,-1)**2,clr(3,3,+1)**2)
cccc A
        call curr6addm(+1,psiv(1,-1,4,1),pmomv(0,4,1),psi(1,-1,1),qbar(0,1),
     &    jvv14(1,-1,3),-clr(3,1,-1),-clr(3,1,+1))
        call curr6addm(+1,psi(1,-1,4),qbar(0,4),psiv(1,-1,1,1),pmomv(0,1,1),
     &    jvv14(1,-1,3),-clr(3,1,-1),-clr(3,1,+1))
cccc Z
        call curr6addm(+1,psiv(1,-1,4,2),pmomv(0,4,2),psi(1,-1,1),qbar(0,1),
     &    jvv14(1,-1,3),-clr(3,2,-1),-clr(3,2,+1))
        call curr6addm(+1,psi(1,-1,4),qbar(0,4),psiv(1,-1,1,2),pmomv(0,1,2),
     &    jvv14(1,-1,3),-clr(3,2,-1),-clr(3,2,+1))
ccc <4|3>
cccc W
        call curr6addm(-1,psivv(1,-1,4,4,3),pmomvv(0,4,4,3),psi(1,-1,3),qbar(0,3),
     &    jvv34(1,-1,3),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psiv(1,-1,4,4),pmomv(0,4,4),psiv(1,-1,3,3),pmomv(0,3,3),
     &    jvv34(1,-1,3),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psi(1,-1,4),qbar(0,4),psivv(1,-1,3,4,3),pmomvv(0,3,4,3),
     &    jvv34(1,-1,3),clr(3,3,-1)**2,clr(3,3,+1)**2)
cccc A
        call curr6addm(+1,psiv(1,-1,4,1),pmomv(0,4,1),psi(1,-1,3),qbar(0,3),
     &    jvv34(1,-1,3),-clr(3,1,-1),-clr(3,1,+1))
        call curr6addm(+1,psi(1,-1,4),qbar(0,4),psiv(1,-1,3,1),pmomv(0,3,1),
     &    jvv34(1,-1,3),-clr(3,1,-1),-clr(3,1,+1))
cccc Z
        call curr6addm(+1,psiv(1,-1,4,2),pmomv(0,4,2),psi(1,-1,3),qbar(0,3),
     &    jvv34(1,-1,3),-clr(3,2,-1),-clr(3,2,+1))
        call curr6addm(+1,psi(1,-1,4),qbar(0,4),psiv(1,-1,3,2),pmomv(0,3,2),
     &    jvv34(1,-1,3),-clr(3,2,-1),-clr(3,2,+1))
cc d quark
ccc <4|1>
cccc W
        call curr6addm(-1,psivv(1,-1,4,3,4),pmomvv(0,4,3,4),psi(1,-1,1),qbar(0,1),
     &    jvv14(1,-1,4),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psiv(1,-1,4,3),pmomv(0,4,3),psiv(1,-1,1,4),pmomv(0,1,4),
     &    jvv14(1,-1,4),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psi(1,-1,4),qbar(0,4),psivv(1,-1,1,3,4),pmomvv(0,1,3,4),
     &    jvv14(1,-1,4),clr(3,3,-1)**2,clr(3,3,+1)**2)
cccc A
        call curr6addm(+1,psiv(1,-1,4,1),pmomv(0,4,1),psi(1,-1,1),qbar(0,1),
     &    jvv14(1,-1,4),-clr(4,1,-1),-clr(4,1,+1))
        call curr6addm(+1,psi(1,-1,4),qbar(0,4),psiv(1,-1,1,1),pmomv(0,1,1),
     &    jvv14(1,-1,4),-clr(4,1,-1),-clr(4,1,+1))
cccc Z
        call curr6addm(+1,psiv(1,-1,4,2),pmomv(0,4,2),psi(1,-1,1),qbar(0,1),
     &    jvv14(1,-1,4),-clr(4,2,-1),-clr(4,2,+1))
        call curr6addm(+1,psi(1,-1,4),qbar(0,4),psiv(1,-1,1,2),pmomv(0,1,2),
     &    jvv14(1,-1,4),-clr(4,2,-1),-clr(4,2,+1))
ccc <4|3>
cccc W
        call curr6addm(-1,psivv(1,-1,4,3,4),pmomvv(0,4,3,4),psi(1,-1,3),qbar(0,3),
     &    jvv34(1,-1,4),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psiv(1,-1,4,3),pmomv(0,4,3),psiv(1,-1,3,4),pmomv(0,3,4),
     &    jvv34(1,-1,4),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psi(1,-1,4),qbar(0,4),psivv(1,-1,3,3,4),pmomvv(0,3,3,4),
     &    jvv34(1,-1,4),clr(3,3,-1)**2,clr(3,3,+1)**2)
cccc A
        call curr6addm(+1,psiv(1,-1,4,1),pmomv(0,4,1),psi(1,-1,3),qbar(0,3),
     &    jvv34(1,-1,4),-clr(4,1,-1),-clr(4,1,+1))
        call curr6addm(+1,psi(1,-1,4),qbar(0,4),psiv(1,-1,3,1),pmomv(0,3,1),
     &    jvv34(1,-1,4),-clr(4,1,-1),-clr(4,1,+1))
cccc Z
        call curr6addm(+1,psiv(1,-1,4,2),pmomv(0,4,2),psi(1,-1,3),qbar(0,3),
     &    jvv34(1,-1,4),-clr(4,2,-1),-clr(4,2,+1))
        call curr6addm(+1,psi(1,-1,4),qbar(0,4),psiv(1,-1,3,2),pmomv(0,3,2),
     &    jvv34(1,-1,4),-clr(4,2,-1),-clr(4,2,+1))

c single W attached (fermion 4 designates u/d)
cc <4|1>
ccc u quark
        call curr6addm(-1,psiv(1,-1,4,4),pmomv(0,4,4),psi(1,-1,1),qbar(0,1),
     &    jv14(1,-1,3),clr(3,3,-1),clr(3,3,+1))
        call curr6addm(-1,psi(1,-1,4),qbar(0,4),psiv(1,-1,1,4),pmomv(0,1,4),
     &    jv14(1,-1,3),clr(3,3,-1),clr(3,3,+1))
ccc d quark
        call curr6addm(-1,psiv(1,-1,4,3),pmomv(0,4,3),psi(1,-1,1),qbar(0,1),
     &    jv14(1,-1,4),clr(3,3,-1),clr(3,3,+1))
        call curr6addm(-1,psi(1,-1,4),qbar(0,4),psiv(1,-1,1,3),pmomv(0,1,3),
     &    jv14(1,-1,4),clr(3,3,-1),clr(3,3,+1))
cc <4|3>
ccc u quark
        call curr6addm(-1,psiv(1,-1,4,4),pmomv(0,4,4),psi(1,-1,3),qbar(0,3),
     &    jv34(1,-1,3),clr(3,3,-1),clr(3,3,+1))
        call curr6addm(-1,psi(1,-1,4),qbar(0,4),psiv(1,-1,3,4),pmomv(0,3,4),
     &    jv34(1,-1,3),clr(3,3,-1),clr(3,3,+1))
ccc d quark
        call curr6addm(-1,psiv(1,-1,4,3),pmomv(0,4,3),psi(1,-1,3),qbar(0,3),
     &    jv34(1,-1,4),clr(3,3,-1),clr(3,3,+1))
        call curr6addm(-1,psi(1,-1,4),qbar(0,4),psiv(1,-1,3,3),pmomv(0,3,3),
     &    jv34(1,-1,4),clr(3,3,-1),clr(3,3,+1))

c zero W
        call curr6(+1,psi(1,-1,4),qbar(0,4),psi(1,-1,1),qbar(0,1),jx14(1,-1))
        call curr6(+1,psi(1,-1,4),qbar(0,4),psi(1,-1,3),qbar(0,3),jx34(1,-1))

!
! Now cycle through position of g and build matrix elements
! 

        do isig=-1,1,2
        do isig2=-1,1,2
        do i=3,4
        do j=3,4
          mat(1,i,j,isig,isig2) = -dotcc(gvv12(1,isig,i),jx34(1,isig2))
     &                            -dotcc(gx12(1,isig),jvv34(1,isig2,j))
          mat(2,i,j,isig,isig2) = dotcc(gvv32(1,isig,i),jx14(1,isig2))
     &                           +dotcc(gx32(1,isig),jvv14(1,isig2,j))
! Only fill the allowed combinations
          if (i.ne.j) then
            mat(3,i,j,isig,isig2) = -dotcc(gv12(1,isig,i),jv34(1,isig2,j))
            mat(4,i,j,isig,isig2) = +dotcc(gv32(1,isig,i),jv14(1,isig2,j))
          endif
        enddo
        enddo
        enddo 
        enddo 

#ifdef WITH_MADGRAPH
!
! Compare with MadGraph
!
      if (ldebug) then
         print*,"MadGraph check qqWWqq"
         print*,"Type: ", selc
         do i=3,4
         do j=3,4
         do isig=-1,1,2
         do isig2=-1,1,2
         do k=1,4
           print*,"mat(",k,",",i,",",j,",",isig,",",isig2,") = ",mat(k,i,j,isig,isig2)
           temp = dreal(mat(k,i,j,isig,isig2))**2+dimag(mat(k,i,j,isig,isig2))**2
           print*,"mg2(",k,",",i,",",j,",",isig,",",isig2,") = ",
     &       -tampr(k+4*(i-3)+8*(j-3),isig,isig2,1)
           temp1 = 
     &        dreal(tampr(k+4*(i-3)+8*(j-3),isig,isig2,1))**2
     &       +dimag(tampr(k+4*(i-3)+8*(j-3),isig,isig2,1))**2
           print*,"temp/temp1 = ",temp/temp1
         enddo
         enddo
         enddo
         read(*,*)
         enddo
         enddo
      endif
#endif

!
! Sum everything up
!

c index of m2s should be understood (almost) bitwise:
c 1: <2|1><4|3> fermion lines (mat(1,...))
c 2: <2|3><4|1> fermion lines (mat(2,...))
c 3: 1+2
c 4: <2|W|1><4|W|3> fermion line (mat(3,...))
c 5: <2|W|3><4|W|1> fermion line (mat(4,...))
c 6: 1+5
c 7: 2+4

      do i=3,4
      do j=3,4
      do isig=-1,1,2
      do isig2=-1,1,2
        m2s(1,i-2,j-2)=m2s(1,i-2,j-2)
     &           +2d0*abs(mat(1,i,j,isig,isig2))**2
        m2s(2,i-2,j-2)=m2s(2,i-2,j-2)
     &           +2d0*abs(mat(2,i,j,isig,isig2))**2
        m2s(4,i-2,j-2)=m2s(4,i-2,j-2)
     &           +2d0*abs(mat(3,i,j,isig,isig2))**2
        m2s(5,i-2,j-2)=m2s(5,i-2,j-2)
     &           +2d0*abs(mat(4,i,j,isig,isig2))**2
        if (isig.eq.isig2) then ! interference only possible if same helicity
          m2s(3,i-2,j-2)=m2s(3,i-2,j-2)
     &           -4d0/3d0*dreal(mat(1,i,j,isig,isig2)*dconjg(mat(2,i,j,isig,isig2)))
          m2s(6,i-2,j-2)=m2s(6,i-2,j-2)
     &           -4d0/3d0*dreal(mat(1,i,j,isig,isig2)*dconjg(mat(4,i,j,isig,isig2)))
          m2s(7,i-2,j-2)=m2s(7,i-2,j-2)
     &           -4d0/3d0*dreal(mat(2,i,j,isig,isig2)*dconjg(mat(3,i,j,isig,isig2)))
        endif
      enddo
      enddo
        m2s(3,i-2,j-2)=m2s(3,i-2,j-2)+m2s(1,i-2,j-2)+m2s(2,i-2,j-2)
        m2s(6,i-2,j-2)=m2s(6,i-2,j-2)+m2s(1,i-2,j-2)+m2s(5,i-2,j-2)
        m2s(7,i-2,j-2)=m2s(7,i-2,j-2)+m2s(2,i-2,j-2)+m2s(4,i-2,j-2)
      enddo
      enddo

!
! Multiply with remaining alpha_s coupling constants
!
      do i=1,7
      do j=1,2
      do k=1,2
        m2s(i,j,k) = m2s(i,j,k) * fpials(2)**2 
      enddo
      enddo
      enddo

c Compare with MadGraph
#ifdef WITH_MADGRAPH
      if (ldebug) then
         print*,"Total partonic amplitude (selc=",selc,")"
         print*,"  uu 1-2:"
         print*,"    VBFNLO:  ", m2s(1,1,1)
         print*,"    MadGraph:", mg2uures(1)/g**4*fpials(2)**2
         print*,"    Ratio:   ", m2s(1,1,1)/mg2uures(1)*g**4/fpials(2)**2
         print*,"  uu 3-2:"
         print*,"    VBFNLO:  ", m2s(2,1,1)
         print*,"    MadGraph:", mg2uures(2)/g**4*fpials(2)**2
         print*,"    Ratio:   ", m2s(2,1,1)/mg2uures(2)*g**4/fpials(2)**2
         print*,"  uu 1-2+3-2:"
         print*,"    VBFNLO:  ", m2s(3,1,1)
         print*,"    MadGraph:", mg2uures(3)/g**4*fpials(2)**2
         print*,"    Ratio:   ", m2s(3,1,1)/mg2uures(3)*g**4/fpials(2)**2
         print*,"  dd 1-2:"
         print*,"    VBFNLO:  ", m2s(1,2,2)
         print*,"    MadGraph:", mg2ddres(1)/g**4*fpials(2)**2
         print*,"    Ratio:   ", m2s(1,2,2)/mg2ddres(1)*g**4/fpials(2)**2
         print*,"  dd 3-2:"
         print*,"    VBFNLO:  ", m2s(2,2,2)
         print*,"    MadGraph:", mg2ddres(2)/g**4*fpials(2)**2
         print*,"    Ratio:   ", m2s(2,2,2)/mg2ddres(2)*g**4/fpials(2)**2
         print*,"  dd 1-2+3-2:"
         print*,"    VBFNLO:  ", m2s(3,2,2)
         print*,"    MadGraph:", mg2ddres(3)/g**4*fpials(2)**2
         print*,"    Ratio:   ", m2s(3,2,2)/mg2ddres(3)*g**4/fpials(2)**2
         print*,"  ud 1-2:"
         print*,"    VBFNLO:  ", m2s(4,2,1)
         print*,"    MadGraph:", mg2udres(1)/g**4*fpials(2)**2
         print*,"    Ratio:   ", m2s(4,2,1)/mg2udres(1)*g**4/fpials(2)**2
         print*,"  ud 3-2:"
         print*,"    VBFNLO:  ", m2s(2,2,1)
         print*,"    MadGraph:", mg2udres(2)/g**4*fpials(2)**2
         print*,"    Ratio:   ", m2s(2,2,1)/mg2udres(2)*g**4/fpials(2)**2
         print*,"  ud 1-2+3-2:"
         print*,"    VBFNLO:  ", m2s(7,2,1)
         print*,"    MadGraph:", mg2udres(3)/g**4*fpials(2)**2
         print*,"    Ratio:   ", m2s(7,2,1)/mg2udres(3)*g**4/fpials(2)**2
         print*,"  du 1-2:"
         print*,"    VBFNLO:  ", m2s(4,1,2)
         print*,"    MadGraph:", mg2dures(1)/g**4*fpials(2)**2
         print*,"    Ratio:   ", m2s(4,1,2)/mg2dures(1)*g**4/fpials(2)**2
         print*,"  du 3-2:"
         print*,"    VBFNLO:  ", m2s(2,1,2)
         print*,"    MadGraph:", mg2dures(2)/g**4*fpials(2)**2
         print*,"    Ratio:   ", m2s(2,1,2)/mg2dures(2)*g**4/fpials(2)**2
         print*,"  du 1-2+3-2:"
         print*,"    VBFNLO:  ", m2s(7,1,2)
         print*,"    MadGraph:", mg2dures(3)/g**4*fpials(2)**2
         print*,"    Ratio:   ", m2s(7,1,2)/mg2dures(3)*g**4/fpials(2)**2
         read(*,*)
      endif
#endif

      end


