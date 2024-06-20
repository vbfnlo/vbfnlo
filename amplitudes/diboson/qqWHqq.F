      subroutine qqWHqq(p, v, selc, m2su, m2sd)

      implicit none

      include 'VBFNLO/utilities/global.inc'
      include 'VBFNLO/utilities/lha.inc'
      include 'VBFNLO/utilities/scales.inc'
      include 'VBFNLO/utilities/coupl.inc'
      include 'tensor.inc'

      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1), XM2(6), XMG(6), B(6,6,6)

      integer j, i, k, mu, chi
      integer sign(4), iflavu, iflavd

      double precision p(0:3, max_p), v(0:3, max_v), m2su(3), m2sd(3)
      double precision pbar(0:3,4), qbar(0:3,4)

      double complex psi(2,-1:1,4),
     &       uuud(2,-1:1,-1:1), uddd(2,-1:1,-1:1), 
     &       psiwh(2,-1:1,4), jwh(1:6,-1:1,3),
     &       gnov14(1:6,-1:1), gnov32(1:6,-1:1), gnov34(1:6,-1:1)

      double precision gmomnov14(0:4,-1:1), gmomnov32(0:4,-1:1), 
     &       gmomnov34(0:4,-1:1),
     &       pmomwh(0:4,-1:1,4) 

      logical ldebug, ldebugm
      ! enable general debugging
      parameter (ldebug=.false.) 
      ! enable matrixelement checking with madgraph
      parameter (ldebugm=.false.) 
    
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
!  u  u  > ve e+  H d  u
!  u  d  > ve e+  H d  d
!( u~ u~ > e- ve~ H d~ u~ )
!( u~ d~ > e- ve~ H d~ d~ )

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

!-----------------------------------------------------
      Case (2)  
!  u  u~ > ve e+  H d  u~
!  u  d~ > ve e+  H d  d~
!( u~ u  > e- ve~ H d~ u  )
!( u~ d  > e- ve~ H d~ d  )
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
!  u  d~ > ve e+  H u~ u  
!  d  d~ > ve e+  H u~ d  
!( u~ d  > e- ve~ H u  u~ )
!( d~ d  > e- ve~ H u  d~ )
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
!  u~ d~ > ve e+  H u~ u~
!  d~ d~ > ve e+  H u~ d~
!( u  d  > e- ve~ H u  u  )
!( d  d  > e- ve~ H u  d  )
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
            qbar(mu, i)=sign(i)*pbar(mu, i)
         enddo
      enddo

#ifdef WITH_MADGRAPH
!
! Write settings in variables for MadGraph
!
      if (ldebugm) then
         fsign(1)= sign(1)
         fsign(2)= sign(2)
         fsign(7)= sign(4)
         fsign(8)= sign(3)
         fsign(3)=+1
         fsign(4)=-1
         fsign(5)=+1
         do mu=0,3
            pall(mu,1) = pbar(mu,1)
            pall(mu,2) = pbar(mu,2)
            pall(mu,7) = pbar(mu,4)
            pall(mu,8) = pbar(mu,3)
            pall(mu,3) = wnumom(mu)
            pall(mu,4) = wemmom(mu)
            pall(mu,5) = higmom(mu) 
         enddo
         call SUD_WHUU(pall, fsign, mg2res)
         call SUD_WHDD(pall, fsign, mg2res)
      endif
#endif

!-----------------------------------------------------
!-----------------------------------------------------

      call psi0m(4, pbar(0,1), sign(1), psi)

      do i=1,2
         do j=-1,1,2
            do k=-1,1,2
               uuud(i, j, k)=0d0
               uddd(i, j, k)=0d0
            enddo
         enddo
      enddo

!
! Now cycle through position of gluon and build currents
! <2|1> is the <d|u> quark line in the common combination
! at the end add gluon propagator
! 


! <4|1>
        call curr6(+1, psi(1,-1,4), qbar(0,4), psi(1,-1,1), qbar(0,1), gnov14(1,-1))
! <2|3>
        call curr6(+1, psi(1,-1,2), qbar(0,2), psi(1,-1,3), qbar(0,3), gnov32(1,-1))
! <4|3>
        call curr6(+1, psi(1,-1,4), qbar(0,4), psi(1,-1,3), qbar(0,3), gnov34(1,-1))
! gluon propagator
!        enddo
        call propagate(+1, gnov14(1,-1), gmomnov14(0,-1))
        call propagate(+1, gnov32(1,-1), gmomnov32(0,-1))
        call propagate(+1, gnov34(1,-1), gmomnov34(0,-1))

!
! Now build the other fermion line, attaching the WH-current as necessary
!
! First add again the vector bosons to the <bra| and |ket>
! Order is important to avoid double counting !!!


! all <bra| and |ket> needed
        call bra2c(psi(1,-1,2),.true., qbar(0,2),-1,
     &    wemnuhmom(0), wtowh(0,1),
     &    psiwh(1,-1,2), pmomwh(0,-1,2)) !w34 to <2|
        call bra2c(psi(1,-1,4),.true., qbar(0,4),-1,
     &    wemnuhmom(0), wtowh(0,1),
     &    psiwh(1,-1,4), pmomwh(0,-1,4)) !w34 to <4|
        call ket2c(psi(1,-1,1),.true., qbar(0,1),-1,
     &    wemnuhmom(0), wtowh(0,1),
     &    psiwh(1,-1,1), pmomwh(0,-1,1)) !w34 to |1>
        call ket2c(psi(1,-1,3),.true., qbar(0,3),-1,
     &    wemnuhmom(0), wtowh(0,1),
     &    psiwh(1,-1,3), pmomwh(0,-1,3)) !w34 to |3>


!
! Now cycle through position of g and build currents
! 

        call curr6   (-1, psiwh(1,-1,2), pmomwh(0,-1,2),
     &         psi(1,-1,1), qbar(0,1), jwh(1,-1,1)) !<2|WH g|1>
        call curr6   (-1, psiwh(1,-1,2), pmomwh(0,-1,2),
     &         psi(1,-1,3), qbar(0,3), jwh(1,-1,2)) !<2|WH g|3>
        call curr6   (-1, psiwh(1,-1,4), pmomwh(0,-1,4),
     &         psi(1,-1,1), qbar(0,1), jwh(1,-1,3)) !<4|WH g|1>
        call curr6add(-1, psi(1,-1,2), qbar(0,2),
     &         psiwh(1,-1,1), pmomwh(0,-1,1), jwh(1,-1,1)) !<2|g WH|1>
        call curr6add(-1, psi(1,-1,2), qbar(0,2),
     &         psiwh(1,-1,3), pmomwh(0,-1,3), jwh(1,-1,2)) !<2|g WH|3>
        call curr6add(-1, psi(1,-1,4), qbar(0,4),
     &         psiwh(1,-1,1), pmomwh(0,-1,1), jwh(1,-1,3)) !<4|g WH|1>

        do chi=-1,1,2
          temp = dotcc(jwh(1,-1,1), gnov34(1, chi))
          uuud(1,-1, chi) = uuud(1,-1, chi) + temp
          uddd(1,-1, chi) = uddd(1,-1, chi) + temp

          temp = dotcc(jwh(1,-1,2), gnov14(1, chi))
          uuud(2,-1, chi) = uuud(2,-1, chi) + temp
          temp = dotcc(jwh(1,-1,3), gnov32(1, chi))
          uddd(2,-1, chi) = uddd(2,-1, chi) + temp
        enddo !chi


#ifdef WITH_MADGRAPH
!
! Compare with MadGraph
!
        if (ldebugm) then
           do chi=-1,1,2
              do i=1,2
c          uduu
                 temp = dreal(uuud(i,-1, chi))**2+dimag(uuud(i,-1, chi))**2
                 temp1 = dreal(tampr(i,-1, chi,1))**2 
     &                  +dimag(tampr(i,-1, chi,1))**2
                 if(abs(temp/temp1 -1).ge.1.0d-5) then
                    print*,"MadGraph check qqWHqq"
                    print*,"Type: ", selc
                    print*,"uu -> du WH and crossed"
                    print*,"mat(", i,",-1,", chi,") = ", uuud(i,-1, chi)
                    print*,"mg2(", i,",-1,", chi,") = ", tampr(i,-1, chi,1)
                    print*,"temp/temp1 = ", temp/temp1
                    print*, "matrix elements qqWHqq.F are only within", abs(temp/temp1-1) 
                 endif      
c          uddd
                 temp = dreal(uddd(i,-1, chi))**2+dimag(uddd(i,-1, chi))**2
                 temp1 = dreal(tampr(i,-1, chi,2))**2
     &                  +dimag(tampr(i,-1, chi,2))**2
                 if(abs(temp/temp1 -1).ge.1.0d-5) then
                    print*,"MadGraph check qqWHqq"
                    print*,"Type: ", selc
                    print*,"ud -> dd WH and crossed"
                    print*,"mat(", i,",-1,", chi,") = ", uddd(i,-1, chi)
                    print*,"mg2(", i,",-1,", chi,") = ", -tampr(i,-1, chi,2)
                    ! uddd has a sign flipped compared to madgraph
                    print*,"temp/temp1 = ", temp/temp1
                    print*, "matrix elements qqWHqq.F are only within", abs(temp/temp1-1) 
                 endif      
              enddo
           enddo
!           pause to check me by hand           
!           read(*,*)
        endif
#endif

!
! Sum everything up
!

      do chi=-1,1,2
        m2su(1)=m2su(1)
     &           +2d0*abs(uuud(1,-1, chi))**2
     &           +2d0*abs(uuud(2,-1, chi))**2
     &           +4d0/3d0*dreal(uuud(1,-1, chi)*dconjg(uuud(2,-1, chi)))
        m2su(2)=m2su(2)
     &           +2d0*abs(uuud(2,-1, chi))**2
        m2su(3)=m2su(3)
     &           +2d0*abs(uuud(1,-1, chi))**2
        m2sd(1)=m2sd(1)
     &           +2d0*abs(uddd(1,-1, chi))**2
     &           +2d0*abs(uddd(2,-1, chi))**2
     &           +4d0/3d0*dreal(uddd(1,-1, chi)*dconjg(uddd(2,-1, chi)))
        m2sd(2)=m2sd(2)
     &           +2d0*abs(uddd(2,-1, chi))**2
        m2sd(3)=m2sd(3)
     &           +2d0*abs(uddd(1,-1, chi))**2
      enddo

!
! Multiply with remaining coupling constants
!
      do i=1,3
        m2su(i) = m2su(i) * fpials(2)**2 * clr(1,3,-1)**2
        m2sd(i) = m2sd(i) * fpials(2)**2 * clr(1,3,-1)**2
      enddo

      end


