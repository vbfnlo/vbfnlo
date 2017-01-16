! This subroutine returns the real emission contribution to the order alpha_s 
! hadronic matrix element squared for the process p p --> e- ve gamma gamma j j
! jet. It is also used as part of the NLO calculation.

      SUBROUTINE m2s_qqWAAjj(bos,nlo,lokt,xi,p,v,rn,xuz,ps_number,NW,m2s)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "tensor.inc"

      integer ps_number,NW
 
      integer i,j,k,mu,
     &  nlo,signpdf1,signpdf2,
     &  bos
      logical lokt(11),fincolcalc,pdfchange
      logical ldebug
      parameter (ldebug=.false.)

      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)

      double precision xi(nx), p(0:3,max_p,max_kin), 
     &   v(0:3,max_v,max_kin), m2s_qqwgaglujfey,
     &   m2s_qqwgaglujf(4,3),pch(0:3,max_p,max_kin),
     &   
     &   
     &   muu(3),mdd(3),muuch(3),mddch(3),m2s_qqwgaglujch(4,3),
     &   q_sf,pdf(-6:6,2),amp,m2s_qqwgaglujdd(4,3),dotrr,
     &   m2s_qqwgaglujddch(4,3),
     &   helastestmom(0:3,7),
     &   mgg(4),mggch(4),m2s(0:11),
     &   
     &   fincol(12),fincolampsq(12),coll3,coll4,
     &   fincolconv,coll7,coll8,pdftmp(-6:6,2),
     &   rn(1)

      double precision dipole(10),
     &   xuz(2,2:11),kinem(11),dipolev(6,6),dipolea(3),dipoleb(3),borncache147(2,2),
     &   borncache(9,4),borncache145(3,2)
      double precision kininv(4,4), ratio(5)
      logical kintrue

      double precision zero,inv2,four,two
      parameter (zero=0d0,inv2=1d0/2d0,four=4d0,two=2d0)

      double complex zm2i(2:3), dotrc
      save zm2i

      common/BORNCACHE/borncache145,borncache147

      double complex help1c,help2c,help3c,
     &   Contract_Trjj
      external dotrr,Contract_Trjj

      real*8 partons(0:7,max_jets,max_kin)
      integer npartons(max_kin)
      common /partons/partons,npartons

      double precision maxeventweight
      common/maxevents/maxeventweight

      parameter(fincolcalc=.true.)

      real*8 rsep
      external rsep

      INTEGER init/0/
      SAVE init

c  helicity selection
      INTEGER h
      COMMON /hcount / h

      if ( init .eq. 0 ) then
        maxeventweight=0d0

        zm2i(2) = 1d0/dcmplx(xm2(2),-xmg(2))
        zm2i(3) = 1d0/dcmplx(xm2(3),-xmg(3))

        if (bos .eq. 411) then
           pdfchange = .false.
           write(6,*) " "
           write(6,*) "W-AAjj amplitude square information:"
           write(6,*) "---------------------------------------------"
           write(6,*) " "
        elseif (bos .eq. 311) then
           pdfchange = .true.
           write(6,*) " "
           write(6,*) "W+AAjj amplitude square information:"
           write(6,*) "---------------------------------------------"
           write(6,*) " "
        else 
           write(6,*) " "
           write(6,*) "m2s_qqWAAjj.F called with invalid process."
           write(6,*) 'bos=',bos
           write(6,*) " "
        endif

        if (ldebug .and. fincolcalc) then
           write(6,*) "Calculation of finite collinear terms included"
        endif

      endif

      init=1

! initialize
      do i=1,11
      kinem(i)=zero
      enddo
      do i=1,6
      do j=1,6
      dipolev(i,j)=zero
      enddo
      enddo
      do j=1,3
      dipolea(j)=zero
      dipoleb(j)=zero
      enddo
      do i=1,10
      dipole(i)=zero
      enddo
      do i=1,3
      muu(i)=zero
      mdd(i)=zero
      muuch(i)=zero
      mddch(i)=zero
      enddo
      do i=1,4
      mgg(i)=zero
      mggch(i)=zero
      enddo
      do i=0,11
      m2s(i)=zero
      enddo
      m2s_qqwgaglujfey=zero
      amp=zero
      coll3=zero
      coll4=zero
      coll7=zero
      coll8=zero

! clear cache
!
      do mu=1,9
       do i=1,4
         borncache(mu,i)=zero
       enddo
      enddo
      do mu=1,3
       do i=1,2
         borncache145(mu,i)=zero
       enddo
      enddo
      do mu=1,2
       do i=1,2
         borncache147(mu,i)=zero
       enddo
      enddo

! exchange incoming particles

      do mu=0,3
      pch(mu,1,1)=p(mu,2,1)
      pch(mu,2,1)=p(mu,1,1)
      pch(mu,3,1)=p(mu,3,1)
      pch(mu,4,1)=p(mu,4,1)
      enddo

!  call topologies for crossed kinematics

      call getyourscalesready_trib(xi,p(0,1,1),v(0,1,1),pdfchange,pdf,nlo,1)

c for WAA precalculate leptonic tensors
c lepton spinors and polarisation vectors
c select helicity: h = random number for lepton & photon helicity (h=1:4)
      h = mod(h,4) + 1
      ie = (-1)**(h+1)
      iu = (-1)**((h-1)/2)

      do i=11,1,-1
        if (lokt(i)) then
! effective W currents
          do mu=0,3
            phot1fourvec(mu)=v(mu,3,i)
            phot2fourvec(mu)=v(mu,4,i)
            wnumom(mu)=v(mu,1,i)
            wemmom(mu)=v(mu,2,i)
            wemnumom(mu)=v(mu,1,i)+v(mu,2,i)
            wemnuga1mom(mu)=v(mu,1,i)+v(mu,2,i)+v(mu,3,i)
            wemnuga2mom(mu)=v(mu,1,i)+v(mu,2,i)+v(mu,4,i)
            wemnugagamom(mu)=v(mu,1,i)+v(mu,2,i)+v(mu,3,i)+v(mu,4,i)
            qw0(mu,i)=wemnumom(mu)
            qa1(mu,i)=phot1fourvec(mu)
            qa2(mu,i)=phot2fourvec(mu)
          enddo
          phot1fourvec(4)=zero
          phot2fourvec(4)=zero
          wemnumom(4)=dotrr(wemnumom,wemnumom)
          wemnuga1mom(4)=dotrr(wemnuga1mom,wemnuga1mom)
          wemnuga2mom(4)=dotrr(wemnuga2mom,wemnuga2mom)
          wemnugagamom(4)=dotrr(wemnugagamom,wemnugagamom)

c use symmetry and always reduce to W+AA case
          if (bos.eq.311) then
            CALL OXXXXX(v(0,1,i),ZERO ,-1,1,wve(1,i))           !ve 
            CALL IXXXXX(v(0,2,i),ZERO ,1,-1,wep(1,i))           !e+ 
          elseif (bos.eq.411) then
            CALL IXXXXX(v(0,1,i),ZERO ,1,-1,wep(1,i))           !e-  -> e+
            CALL OXXXXX(v(0,2,i),ZERO ,-1,1,wve(1,i))           !ve~ -> ve
          endif
          CALL VXXXXX(phot1fourvec(0),ZERO ,ie,1,atau1(1,i))          !A1
          CALL VXXXXX(phot2fourvec(0),ZERO ,iu,1,atau2(1,i))          !A2
          CALL JIOXXX(wep(1,i),wve(1,i),GWF,WMASS,WWIDTH, w0(1,i))!W+- 

          if (with_anom) then ! anomalous gauge boson couplings
c            using global form factor for all tensors of one phase space point
c            this ensures proper cancellations for anomalous contributions
c            energy scale is invariant WAA mass
             call anomal_formfactor(wemnugagamom(0),wemnumom(0),phot1fourvec(0),phot2fourvec(0))

             call wptowanew_anomal(v(0,1,i),i,h,wpwa1(0,i),1,NW)!W+ -> W+ A -> ve e+ a
             call wptowanew_anomal(v(0,1,i),i,h,wpwa2(0,i),2,NW)!W+ -> W+ A -> ve e+ a
             call wptowaa_anomal(v(0,1,i),i,h,wpwaa(0,i),0,NW)  !W+ -> ve e+ a a

          else !SM

             call wptowanew(v(0,1,i),i,h,wpwa1(0,i),1,NW)!W+ -> W+ A -> ve e+ a
             call wptowanew(v(0,1,i),i,h,wpwa2(0,i),2,NW)!W+ -> W+ A -> ve e+ a
             call wptowaa(v(0,1,i),i,h,wpwaa(0,i),0,NW)  !W+ -> ve e+ a a

          endif ! anom or SM

c define epsilon~ = espilon - k * k.epsilon/M**2 for wpwa1, wpwa2 and wpwaa
          help1c = -zm2i(3) * dotrc(wemnuga1mom(0),wpwa1(0,i))
          help2c = -zm2i(3) * dotrc(wemnuga2mom(0),wpwa2(0,i))
          help3c = -zm2i(3) * dotrc(wemnugagamom(0),wpwaa(0,i))
          do mu=0,3
            wpwa1(mu,i) = wpwa1(mu,i) + wemnuga1mom(mu) * help1c
            wpwa2(mu,i) = wpwa2(mu,i) + wemnuga2mom(mu) * help2c
            wpwaa(mu,i) = wpwaa(mu,i) + wemnugagamom(mu) * help3c
          enddo
        endif
      enddo

!
!  Born - Kinematics
!
      if (lokt(1)) then

      do i=1,4
        if((Loops_sub_NLO.eq.1).or.(sub_number.eq.1)) then
          call qqWAAgg(p,v,i,mgg(i))  ! Processes 1-4
          if (i.ne.4) then
            call qqWAAgg(pch,v,i,mggch(i))  !Processes 5-7
          endif
        endif

        if((Loops_sub_NLO.eq.1).or.(sub_number.eq.2)) then
          call qqWAAqq(p,v,i,muu,mdd)   ! Processes 8-11
          call qqWAAqq(pch,v,i,muuch,mddch)  ! Processes 12-15
        endif
        do k=1,3
          m2s_qqwgaglujdd(i,k)=mdd(k)
          m2s_qqwgaglujddch(i,k)=mddch(k)
          m2s_qqwgaglujf(i,k)=muu(k)
          m2s_qqwgaglujch(i,k)=muuch(k)
        enddo
      enddo

!                                                                                     PDF CHECK AGAINST DIPOLES
      amp= zero
     &   +m2s_qqwgaglujf(1,1)   *(pdf(-2,1)*pdf(-2,2)+pdf(-4,1)*pdf(-4,2))   !u~u~ > u~d~   !pdf check   ! 1
     &   +m2s_qqwgaglujf(2,1)   *(pdf(-2,1)*pdf( 2,2)+pdf(-4,1)*pdf( 4,2))   !u~u  > u d~   !pdf check   ! 2-3
     &   +m2s_qqwgaglujch(2,1)  *(pdf( 2,1)*pdf(-2,2)+pdf( 4,1)*pdf(-4,2))                  !pdf check 
     &   +m2s_qqwgaglujf(3,1)   *(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2)) !u~d  > u u~   !pdf check	  ! 4-5
     &   +m2s_qqwgaglujch(3,1)  *(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                  !pdf check
     &   +m2s_qqwgaglujf(4,1)*inv2   *(pdf( 2,1)*pdf( 1,2)+pdf( 4,1)*pdf( 3,2)) !u d  > u u  !pdf check   ! 6-7
     &   +m2s_qqwgaglujch(4,1)*inv2  *(pdf( 1,1)*pdf( 2,2)+pdf( 3,1)*pdf( 4,2))              !pdf check
     &   +m2s_qqwgaglujddch(4,1)*(pdf( 1,1)*pdf( 1,2)+pdf( 3,1)*pdf( 3,2))   !d d  > d u    !pdf check   ! 8
     &   +m2s_qqwgaglujddch(1,1)*inv2  *(pdf(-1,1)*pdf(-2,2)+pdf(-3,1)*pdf(-4,2)) !d~u~ > d~d~ !pdf check   ! 9-10
     &   +m2s_qqwgaglujdd(1,1)*inv2*(pdf(-2,1)*pdf(-1,2)+pdf(-4,1)*pdf(-3,2))              !pdf check 
     &   +m2s_qqwgaglujdd(2,1)  *(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))   !u~d  > d d~   !pdf check   ! 11-12
     &   +m2s_qqwgaglujddch(2,1)*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                  !pdf check
     &   +m2s_qqwgaglujdd(3,1)  *(pdf(-1,1)*pdf( 1,2)+pdf(-3,1)*pdf( 3,2))   !d~d > u d~    !pdf check   ! 13-14
     &   +m2s_qqwgaglujddch(3,1)*(pdf( 1,1)*pdf(-1,2)+pdf( 3,1)*pdf(-3,2))                  !pdf check 
!!!!!!!!!!!!
! ! ! u~u > s~c, c~c > d~u  ! no dipole!!!!!!!!!!!!
     &   +m2s_qqwgaglujf(2,2)   *(pdf(-2,1)*pdf( 2,2)+pdf(-4,1)*pdf( 4,2))                  !pdf check   ! 15-16
     &   +m2s_qqwgaglujch(2,2)  *(pdf( 2,1)*pdf(-2,2)+pdf( 4,1)*pdf(-4,2))                  !pdf check
! ! ! u~d > ss~, c~s > dd~ 
     &   +m2s_qqwgaglujdd(2,2)  *(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))                  !pdf check   ! 17-18
     &   +m2s_qqwgaglujddch(2,2)*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                  !pdf check
! ! ! u~d > cc~, c~s > u u~ 
     &   +m2s_qqwgaglujf(3,2)  *(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))                   !pdf check   ! 19-20
     &   +m2s_qqwgaglujch(3,2) *(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                   !pdf check
! ! ! u~s > u~c, c~d > c~u
     &   +m2s_qqwgaglujf(3,3) *(pdf(-2,1)*pdf( 3,2)+pdf(-4,1)*pdf( 1,2))                    !pdf check   ! 21-22
     &   +m2s_qqwgaglujch(3,3)*(pdf( 3,1)*pdf(-2,2)+pdf( 1,1)*pdf(-4,2))                    !pdf check 
! ! ! u~s~> d~s~, c~d~ > s~d~
     &   +m2s_qqwgaglujddch(1,3)  *(pdf(-2,2)*pdf(-3,1)+pdf(-4,2)*pdf(-1,1))                  !pdf check   ! 23-24
     &   +m2s_qqwgaglujdd(1,3)*(pdf(-3,2)*pdf(-2,1)+pdf(-1,2)*pdf(-4,1))                  !pdf check
! ! ! u~s > d~s, c~d > s~d 
     &   +m2s_qqwgaglujdd(2,3)  *(pdf(-2,1)*pdf( 3,2)+pdf(-4,1)*pdf( 1,2))                  !pdf check   ! 25-26
     &   +m2s_qqwgaglujddch(2,3)*(pdf( 3,1)*pdf(-2,2)+pdf( 1,1)*pdf(-4,2))                  !pdf check
! ! ! u~c~> u~s~, c~u~>c~d~
     &   +m2s_qqwgaglujf(1,2)  *(pdf(-2,1)*pdf(-4,2)+pdf(-4,1)*pdf(-2,2))                   !pdf check   ! 27-28
     &   +m2s_qqwgaglujch(1,2)*(pdf(-4,1)*pdf(-2,2)+pdf(-2,1)*pdf(-4,2))                    !pdf check
! ! ! ! u~c> d~c, c~u > s~u
     &   +m2s_qqwgaglujf(2,3)  *(pdf(-2,1)*pdf( 4,2)+pdf(-4,1)*pdf( 2,2))                   !pdf check   ! 29-30 x
     &   +m2s_qqwgaglujch(2,3)*(pdf( 4,1)*pdf(-2,2)+pdf( 2,1)*pdf(-4,2))                    !pdf check
! ! ! u s> uc, cd > c u
     &   +m2s_qqwgaglujf(4,3)  *(pdf( 2,1)*pdf( 3,2)+pdf( 4,1)*pdf( 1,2))                   !pdf check   ! 31-32 
     &   +m2s_qqwgaglujch(4,3)*(pdf( 3,1)*pdf( 2,2)+pdf( 1,1)*pdf( 4,2))                    !pdf check
! ! ! d~d> s~c, s~s > d~u ! no dipole
     &   +m2s_qqwgaglujdd(3,2)  *(pdf(-1,1)*pdf( 1,2)+pdf(-3,1)*pdf( 3,2))                  !pdf check   ! 33-34
     &   +m2s_qqwgaglujddch(3,2)*(pdf( 1,1)*pdf(-1,2)+pdf( 3,1)*pdf(-3,2))                  !pdf check 
! ! ! d~s> d~c, s~d > s~u
     &   +m2s_qqwgaglujdd(3,3)  *(pdf(-1,1)*pdf( 3,2)+pdf(-3,1)*pdf( 1,2))                  !pdf check   ! 35-36 x
     &   +m2s_qqwgaglujddch(3,3)*(pdf( 3,1)*pdf(-1,2)+pdf( 1,1)*pdf(-3,2))                  !pdf check
! ! ! d s> u s, s d > c d   
     &   +m2s_qqwgaglujddch(4,3)*(pdf( 1,1)*pdf( 3,2)+pdf( 3,1)*pdf( 1,2))                  !pdf check   ! 37-38
     &   +m2s_qqwgaglujdd(4,3)  *(pdf( 3,1)*pdf( 1,2)+pdf( 1,1)*pdf( 3,2))                  !pdf check
! 
      m2s_qqwgaglujfey=amp/36d0
! 
!       
      amp=zero
! ! ! ! u~d> g g, c~s > g g
     &   +(mgg(1)  *(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))                              !pdf check    !39-40
     &   + mggch(1)*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                              !pdf check
     &    )/72d0
! ! u~g> d~g, c~g > s~g
     &   +(mgg(2)  *(pdf(-2,1)*pdf( 0,2)+pdf(-4,1)*pdf( 0,2))                              !pdf check    !41-42
     &   + mggch(2)*(pdf( 0,1)*pdf(-2,2)+pdf( 0,1)*pdf(-4,2))                              !pdf check
     &    )/96d0                        
! ! ! ! d g> u g, s g > c g
     &   +(mgg(3) *(pdf( 1,1)*pdf( 0,2)+pdf( 3,1)*pdf( 0,2))                                !pdf check    !43-44 x
     &   +mggch(3)*(pdf( 0,1)*pdf( 1,2)+pdf( 0,1)*pdf( 3,2))                               !pdf check
     &   )/96d0
! ! ! g g> u d~, g g > c s~
     &   + 2d0*mgg(4)/256d0  *(pdf(0,1)*pdf( 0,2))                                               !pdf check    !45 x

      m2s_qqwgaglujfey=m2s_qqwgaglujfey+amp

!       do mu=0,3
!       helastestmom(mu,1)=p(mu,1,1)
!       helastestmom(mu,2)=p(mu,2,1)
!       helastestmom(mu,3)=v(mu,1,1)
!       helastestmom(mu,4)=v(mu,2,1)
!       helastestmom(mu,5)=v(mu,3,1)
!       helastestmom(mu,6)=p(mu,3,1)
!       helastestmom(mu,7)=p(mu,4,1)
!       enddo
! 
!       CALL SCBDSBD(helastestmom,help1)
! 
!       print*, help1,m2s_qqwgaglujdd(3,3)/36d0/help1
!       STOP
!       endif

      endif

! PDF(5, 4, 3, 2, 1, 0, -1, ......, -5)
! for(b, c, s, u, d, g, d_bar, ..., b_bar)
! (switched u <-> d comp. to CTEQ)


!       do i=2,11
!       do mu=0,3
!       print*, v(mu,3,1)-v(mu,3,i),i
!       enddo
!       print*,"---------------------------_"
!       enddo
!       STOP

! Kinematics 1

      if (lokt(2)) then !.and.dotrr(p(0,1,1),p(0,4,1)).le.1d0) then   ! Dipole 1,43

      call getyourscalesready_trib(xi,p(0,1,2),v(0,1,2),pdfchange,pdftmp,nlo,2)

      Call daisjdipole67_trib(2,1,1,xuz(1,2),p(0,4,1),p(0,3,1),p(0,1,1),p(0,1,2),v(0,1,2),dipolea(1),ps_number)
      Call daisjdipole67_trib(2,2,1,xuz(1,2),p(0,4,1),p(0,3,1),p(0,1,1),p(0,1,2),v(0,1,2),dipoleb(1),ps_number)
      Call daisjdipole65_trib(2,1,1,xuz(1,2),p(0,4,1),p(0,3,1),p(0,1,1),p(0,1,2),v(0,1,2),dipolev(1,1),ps_number)
      Call daisjdipole65_trib(2,1,2,xuz(1,2),p(0,4,1),p(0,3,1),p(0,1,1),p(0,1,2),v(0,1,2),dipolev(1,2),ps_number)

      Call daisjdipole65_trib(2,2,1,xuz(1,2),p(0,4,1),p(0,3,1),p(0,1,1),p(0,1,2),v(0,1,2),dipolev(1,3),ps_number)

      if((Loops_sub_NLO.ne.1).and.(sub_number.eq.1)) then
      dipole(1)=zero
      dipole(2)=zero
      else
      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)
      endif

      kinem(1)= zero 
     &        + dipole(1)*(pdf(-2,1)*pdf(-2,2)+pdf(-4,1)*pdf(-4,2))    !u~u~ > u~d~  !pdf check   !C11
     &        + dipole(1)*(pdf( 2,1)*pdf(-2,2)+pdf( 4,1)*pdf(-4,2))    !u~u  > u d~  !pdf check   !C21
     &        + dipole(2)*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))    !u~d  > u u~  !pdf check   !C31 
     &        + dipole(2)*inv2*(pdf( 2,1)*pdf( 1,2)+pdf( 4,1)*pdf( 3,2))!u d  > u u   !pdf check  ! 6-7 !C41
     &        + dipole(2)*(pdf( 1,1)*pdf( 1,2)+pdf( 3,1)*pdf( 3,2))    !d d  > d u   !pdf check  ! 8 !C42
     &        + dipole(1)*inv2*(pdf(-1,1)*pdf(-2,2)+pdf(-3,1)*pdf(-4,2))!d~u~ > d~d~  !pdf check   !C12
     &        + dipole(1)*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))    !u~d  > d d~  !pdf check    !C22
     &        + dipole(2)*(pdf(-1,1)*pdf( 1,2)+pdf(-3,1)*pdf( 3,2))    !d~d  > u d~  !pdf check    !C32
!!!!!!-
     &        + dipole(2)*(pdf(-2,1)*pdf( 3,2)+pdf(-4,1)*pdf( 1,2))    !u~s  > u~c , c~d  > c~u   !pdf check
     &        + dipole(1)*(pdf(-2,2)*pdf(-3,1)+pdf(-4,2)*pdf(-1,1))    !u~s~ > d~s~, c~d~ > s~d~  !pdf check
     &        + dipole(1)*(pdf( 3,1)*pdf(-2,2)+pdf( 1,1)*pdf(-4,2))    !u~s  > d~s , c~d  > s~d   !pdf check    !25-26
     &        + dipole(1)*(pdf(-2,1)*pdf(-4,2)+pdf(-4,1)*pdf(-2,2))    !u~c~ > u~s~, c~u~ > c~d~  !pdf check
     &        + dipole(1)*(pdf( 4,1)*pdf(-2,2)+pdf( 2,1)*pdf(-4,2))    !u~c  > d~c , c~u  > s~u   !pdf check   !29-30
     &        + dipole(2)*(pdf( 2,1)*pdf( 3,2)+pdf( 4,1)*pdf( 1,2))    !u s  > u c , c d  > c u   !pdf check
     &        + dipole(2)*(pdf(-1,1)*pdf( 3,2)+pdf(-3,1)*pdf( 1,2))    !d~s  > d~c , s~d  > s~u   !pdf check   !35-36
     &        + dipole(2)*(pdf( 3,1)*pdf( 1,2)+pdf( 1,1)*pdf( 3,2))    !d s  > u s , s d  > c d   !pdf check    !37-38
!! !-------------------------------------------------------------------------------------------------------------
     &        + dipolev(1,1)*inv2*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))   !u~d  > g g , c~s > g g  !pdf check
     &        + dipolev(1,2)*inv2*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2)) !pdf check
!! !-------------------------------------------------------------------------------------------------------------
     &        + dipolev(4,2)*(pdf( 0,1)*pdf(-2,2)+pdf( 0,1)*pdf(-4,2))   !u~g> d g, c~g > s g     !pdf check  41-42
!! !-------------------------------------------------------------------------------------------------------------
     &        + dipolev(4,1)*(pdf( 0,1)*pdf( 1,2)+pdf( 0,1)*pdf( 3,2))   !d g  > u g              !pdf check   43-44
!! !-------------------------------------------------------------------------------------------------------------
     &        + 2d0*dipolev(5,3)*(pdf(0,1)*pdf( 0,2))                    !g g > u d~              !pdf check 

      endif

! Kinematics 2

      if (lokt(3)) then ! .and.dotrr(p(0,2,1),p(0,4,1)).le.1d0) then   ! Dipole 2,43

      call getyourscalesready_trib(xi,p(0,1,3),v(0,1,3),pdfchange,pdftmp,nlo,3)

      Call daisjdipole67_trib(3,1,1,xuz(1,3),p(0,4,1),p(0,3,1),p(0,2,1),p(0,1,3),v(0,1,3),dipolea(1),ps_number)
      Call daisjdipole67_trib(3,2,1,xuz(1,3),p(0,4,1),p(0,3,1),p(0,2,1),p(0,1,3),v(0,1,3),dipoleb(1),ps_number)
      Call daisjdipole65_trib(3,1,2,xuz(1,3),p(0,4,1),p(0,3,1),p(0,2,1),p(0,1,3),v(0,1,3),dipolev(1,1),ps_number)
      Call daisjdipole65_trib(3,1,1,xuz(1,3),p(0,4,1),p(0,3,1),p(0,2,1),p(0,1,3),v(0,1,3),dipolev(1,2),ps_number)
      Call daisjdipole65_trib(3,2,1,xuz(1,3),p(0,4,1),p(0,3,1),p(0,2,1),p(0,1,3),v(0,1,3),dipolev(1,3),ps_number)

      if((Loops_sub_NLO.ne.1).and.(sub_number.eq.1)) then
      dipole(1)=zero
      dipole(2)=zero
      else
      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)
      endif

      kinem(2)= zero
     &        + dipole(1)*(pdf(-2,1)*pdf(-2,2)+pdf(-4,1)*pdf(-4,2))   !u~u~ > u~d~  !pdf check  !C11
     &        + dipole(1)*(pdf(-2,1)*pdf( 2,2)+pdf(-4,1)*pdf( 4,2))   !u~u  > u d~  !pdf check  !C21 
     &        + dipole(2)*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))   !u~d  > u u~  !pdf check  !C31
     &        + dipole(2)*inv2*(pdf( 1,1)*pdf( 2,2)+pdf( 3,1)*pdf( 4,2)) !u d   > u u  !pdf check !6-7 !C41
     &        + dipole(2)*(pdf( 1,1)*pdf( 1,2)+pdf( 3,1)*pdf( 3,2))   !d d  > d u   !pdf check   ! 8  !C42
     &        + dipole(1)*inv2*(pdf(-2,1)*pdf(-1,2)+pdf(-4,1)*pdf(-3,2))   !d~u~ > d~d~  !pdf check !C12
     &        + dipole(1)*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))   !u~d  > d d~  !pdf check !C22
     &        + dipole(2)*(pdf( 1,1)*pdf(-1,2)+pdf( 3,1)*pdf(-3,2))   !d~d  > u d~  !pdf check !C32
!!!!!!!
     &        + dipole(2)*(pdf( 3,1)*pdf(-2,2)+pdf( 1,1)*pdf(-4,2))   !u~s  > u~c , c~d  > c~u   !pdf check !
     &        + dipole(1)*(pdf(-3,2)*pdf(-2,1)+pdf(-1,2)*pdf(-4,1))   !u~s~ > d~s~, c~d~ > s~d~  !pdf check
     &        + dipole(1)*(pdf(-2,1)*pdf( 3,2)+pdf(-4,1)*pdf( 1,2))   !u~s  > d~s , c~d  > s~d   !pdf check !25-26
     &        + dipole(1)*(pdf(-4,1)*pdf(-2,2)+pdf(-2,1)*pdf(-4,2))   !u~c~ > u~s~, c~u~ > c~d~  !pdf check
     &        + dipole(1)*(pdf(-2,1)*pdf( 4,2)+pdf(-4,1)*pdf( 2,2))   !u~c  > d~c , c~u  > s~u   !pdf check !29-30
     &        + dipole(2)*(pdf( 3,1)*pdf( 2,2)+pdf( 1,1)*pdf( 4,2))   !u s  > u c , c d  > c u   !pdf check
     &        + dipole(2)*(pdf( 3,1)*pdf(-1,2)+pdf( 1,1)*pdf(-3,2))   !d~s  > d~c , s~d  > s~u   !pdf check !35-36
     &        + dipole(2)*(pdf( 1,1)*pdf( 3,2)+pdf( 3,1)*pdf( 1,2))   !d s  > u s , s d  > c d   !pdf check  !37-38
!! !-------------------------------------------------------------------------------------------------------------
     &        + dipolev(1,1)*inv2*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))  !u~d  > g g , c~s > g g    !pdf check
     &        + dipolev(1,2)*inv2*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                             !pdf check
!! !-------------------------------------------------------------------------------------------------------------
     &        + dipolev(4,1)*(pdf(-2,1)*pdf( 0,2)+pdf(-4,1)*pdf( 0,2))  !u~g  > d g , c~g > s g    !pdf check  41-42
!! !-------------------------------------------------------------------------------------------------------------
     &        + dipolev(4,2)*(pdf( 1,1)*pdf( 0,2)+pdf( 3,1)*pdf( 0,2))  !d g  > u g               !pdf check 43-44
!! !-------------------------------------------------------------------------------------------------------------
     &        +  2d0*dipolev(5,3)*(pdf(0,1)*pdf( 0,2))                    !g g > u d~              !pdf check 

      endif

! Kinematics 3

      if (lokt(4)) then !.and.dotrr(p(0,1,1),p(0,4,1)).le.1d0) then   ! Dipole 12,4

      call getyourscalesready_trib(xi,p(0,1,4),v(0,1,4),pdfchange,pdftmp,nlo,4)
 
      Call daibdipole147_trib(4,1,1,xuz(1,4),p(0,4,1),p(0,1,1),p(0,2,1),p(0,1,4),v(0,1,4),dipolea(1),ps_number)
      Call daibdipole147_trib(4,2,1,xuz(1,4),p(0,4,1),p(0,1,1),p(0,2,1),p(0,1,4),v(0,1,4),dipoleb(1),ps_number)
      Call daibdipole145_trib(4,1,1,xuz(1,4),p(0,4,1),p(0,1,1),p(0,2,1),p(0,1,4),v(0,1,4),dipolev(1,1),ps_number)
      Call daibdipole145_trib(4,1,2,xuz(1,4),p(0,4,1),p(0,1,1),p(0,2,1),p(0,1,4),v(0,1,4),dipolev(1,2),ps_number)

      if((Loops_sub_NLO.ne.1).and.(sub_number.eq.1)) then
      dipole(1)=zero
      dipole(2)=zero
      else
      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)
      endif

      Call daibdipole145_trib(4,2,1,xuz(1,4),p(0,4,1),p(0,1,1),p(0,2,1),p(0,1,4),v(0,1,4),dipolev(1,3),ps_number)

      kinem(3)= zero 
     &        + dipole(1)*(pdf(-2,1)*pdf(-2,2)+pdf(-4,1)*pdf(-4,2))    !u~u~ > u~d~  !pdf check
     &        + dipole(1)*(pdf( 2,1)*pdf(-2,2)+pdf( 4,1)*pdf(-4,2))    !u~u  > u d~  !pdf check
     &        + dipole(2)*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))    !u~d  > u u~  !pdf check
     &        + dipole(2)*inv2*(pdf( 2,1)*pdf( 1,2)+pdf( 4,1)*pdf( 3,2))    !u d  > u u   !pdf check !6-7
     &        + dipole(2)*(pdf( 1,1)*pdf( 1,2)+pdf( 3,1)*pdf( 3,2))    !d d  > d u   !pdf check    ! 8
     &        + dipole(1)*inv2*(pdf(-1,1)*pdf(-2,2)+pdf(-3,1)*pdf(-4,2)) !d~u~ > d~d~  !pdf check
     &        + dipole(1)*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))    !u~d  > d d~  !pdf check
     &        + dipole(2)*(pdf(-1,1)*pdf( 1,2)+pdf(-3,1)*pdf( 3,2))    !d~d  > u d~  !pdf check
!!!!!!!!!
     &        + dipole(2)*(pdf(-2,1)*pdf( 3,2)+pdf(-4,1)*pdf( 1,2))    !u~s  > u~c , c~d  > c~u   !pdf check
     &        + dipole(1)*(pdf(-2,2)*pdf(-3,1)+pdf(-4,2)*pdf(-1,1))    !u~s~ > d~s~, c~d~ > s~d~  !pdf check
     &        + dipole(1)*(pdf( 3,1)*pdf(-2,2)+pdf( 1,1)*pdf(-4,2))    !u~s  > d~s , c~d  > s~d   !pdf check  !25-26
     &        + dipole(1)*(pdf(-2,1)*pdf(-4,2)+pdf(-4,1)*pdf(-2,2))    !u~c~ > u~s~, c~u~ > c~d~  !pdf check
     &        + dipole(1)*(pdf( 4,1)*pdf(-2,2)+pdf( 2,1)*pdf(-4,2))    !u~c  > d~c , c~u  > s~u   !pdf check   !29-30
     &        + dipole(2)*(pdf( 2,1)*pdf( 3,2)+pdf( 4,1)*pdf( 1,2))    !u s  > u c , c d  > c u   !pdf check
     &        + dipole(2)*(pdf(-1,1)*pdf( 3,2)+pdf(-3,1)*pdf( 1,2))    !d~s  > d~c , s~d  > s~u   !pdf check  !35-36
     &        + dipole(2)*(pdf( 3,1)*pdf( 1,2)+pdf( 1,1)*pdf( 3,2))    !d s  > u s , s d  > c d   !pdf check   !37-38
!!! !-------------------------------------------------------------------------------------------------------------
     &        + dipolev(1,1)*inv2*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))   !u~d  > g g , c~s > g g    !pdf check
     &        + dipolev(1,2)*inv2*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                              !pdf check
!! !-------------------------------------------------------------------------------------------------------------
     &        + dipolev(4,2)*(pdf( 0,1)*pdf(-2,2)+pdf( 0,1)*pdf(-4,2))      !u~g  > d g , c~g > s g    !pdf check 41-42
!!! !-------------------------------------------------------------------------------------------------------------
     &        + dipolev(4,1)*(pdf( 0,1)*pdf( 1,2)+pdf( 0,1)*pdf( 3,2))      !d g  > u g              !pdf check 43-44
!! !-------------------------------------------------------------------------------------------------------------
     &        + 2d0*dipolev(5,3) *(pdf(0,1)*pdf( 0,2))                        !g g > u d~              !pdf check  45

!
!     Fill Cache
!
      borncache(1,1)=borncache147(1,1)
      borncache(2,1)=borncache147(1,2)
      borncache(3,1)=borncache147(2,1)
      borncache(4,1)=borncache147(2,2)
      borncache(5,1)=borncache145(1,1)
      borncache(6,1)=borncache145(1,2)
      borncache(7,1)=borncache145(2,1)
      borncache(8,1)=borncache145(2,2)
      borncache(9,1)=borncache145(3,1)


      endif

! Kinematics 4

      if (lokt(5)) then !.and.dotrr(p(0,2,1),p(0,4,1)).le.1d0) then   ! Dipole 21,4

      call getyourscalesready_trib(xi,p(0,1,5),v(0,1,5),pdfchange,pdftmp,nlo,5)

      Call daibdipole147_trib(5,1,1,xuz(1,5),p(0,4,1),p(0,2,1),p(0,1,1),p(0,1,5),v(0,1,5),dipolea(1),ps_number)
      Call daibdipole147_trib(5,2,1,xuz(1,5),p(0,4,1),p(0,2,1),p(0,1,1),p(0,1,5),v(0,1,5),dipoleb(1),ps_number)
      Call daibdipole145_trib(5,1,2,xuz(1,5),p(0,4,1),p(0,2,1),p(0,1,1),p(0,1,5),v(0,1,5),dipolev(1,1),ps_number)
      Call daibdipole145_trib(5,1,1,xuz(1,5),p(0,4,1),p(0,2,1),p(0,1,1),p(0,1,5),v(0,1,5),dipolev(1,2),ps_number)

      Call daibdipole145_trib(5,2,1,xuz(1,5),p(0,4,1),p(0,2,1),p(0,1,1),p(0,1,5),v(0,1,5),dipolev(1,3),ps_number)

      if((Loops_sub_NLO.ne.1).and.(sub_number.eq.1)) then
      dipole(1)=zero
      dipole(2)=zero
      else
      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)
      endif

      kinem(4)= zero 
     &        + dipole(1)*(pdf(-2,1)*pdf(-2,2)+pdf(-4,1)*pdf(-4,2))     !u~u~  > u~d~  !pdf check
     &        + dipole(1)*(pdf(-2,1)*pdf( 2,2)+pdf(-4,1)*pdf( 4,2))     !u~u   > u d~  !pdf check
     &        + dipole(2)*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))     !u~d   > u u~  !pdf check
     &        + dipole(2)*inv2*(pdf( 1,1)*pdf( 2,2)+pdf( 3,1)*pdf( 4,2))     !u d   > u u   !pdf check !6-7
     &        + dipole(2)*(pdf( 1,1)*pdf( 1,2)+pdf( 3,1)*pdf( 3,2))     !d d   > d u   !pdf check  ! 8
     &        + dipole(1)*inv2*(pdf(-2,1)*pdf(-1,2)+pdf(-4,1)*pdf(-3,2))     !d~u~  > d~d~  !pdf check
     &        + dipole(1)*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))     !u~d   > d d~  !pdf check
     &        + dipole(2)*(pdf( 1,1)*pdf(-1,2)+pdf( 3,1)*pdf(-3,2))     !d~d   > u d~  !pdf check
!!!!!
     &        + dipole(2)*(pdf( 3,1)*pdf(-2,2)+pdf( 1,1)*pdf(-4,2))     !u~s   > u~c, c~d  > c~u   !pdf check 
     &        + dipole(1)*(pdf(-3,2)*pdf(-2,1)+pdf(-1,2)*pdf(-4,1))     !u~s~ > d~s~, c~d~ > s~d~  !pdf check
     &        + dipole(1)*(pdf(-2,1)*pdf( 3,2)+pdf(-4,1)*pdf( 1,2))     !u~s  > d~s , c~d  > s~d   !pdf check !25-26
     &        + dipole(1)*(pdf(-4,1)*pdf(-2,2)+pdf(-2,1)*pDf(-4,2))     !u~c~ > u~s~, c~u~ > c~d~  !pdf check
     &        + dipole(1)*(pdf(-2,1)*pdf( 4,2)+pdf(-4,1)*pdf( 2,2))     !u~c  > d~c , c~u  > s~u   !pdf check  !29-30
     &        + dipole(2)*(pdf( 3,1)*pdf( 2,2)+pdf( 1,1)*pdf( 4,2))     !u s  > u c , c d  > c u   !pdf check
     &        + dipole(2)*(pdf( 3,1)*pdf(-1,2)+pdf( 1,1)*pdf(-3,2))     !d~s  > d~c , s~d  > s~u   !pdf check  !35-36
     &        + dipole(2)*(pdf( 1,1)*pdf( 3,2)+pdf( 3,1)*pdf( 1,2))     !d s  > u s , s d  > c d   !pdf check   !37 -38
! !-------------------------------------------------------------------------------------------------------------
     &     + dipolev(1,1)*inv2*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))     !u~d  > g g , c~s  > g g    !pdf check
     &     + dipolev(1,2)*inv2*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                                 !pdf check  
!! !!-------------------------------------------------------------------------------------------------------------
     &     + dipolev(4,1)*(pdf(-2,1)*pdf( 0,2)+pdf(-4,1)*pdf( 0,2))     !u~g  > d g , c~g  > s g    !pdf check   41-42
!! !!-------------------------------------------------------------------------------------------------------------
     &     + dipolev(4,2)*(pdf( 1,1)*pdf( 0,2)+pdf( 3,1)*pdf( 0,2))    !d g  > u g                  !pdf check 43-44
!! !!-------------------------------------------------------------------------------------------------------------
     &     + 2d0*dipolev(5,3) *(pdf(0,1)*pdf( 0,2))                     !g g > u d~              !pdf check 

!
!     Fill Cache
!
      borncache(1,3)=borncache147(1,1)
      borncache(2,3)=borncache147(1,2)
      borncache(3,3)=borncache147(2,1)
      borncache(4,3)=borncache147(2,2)
      borncache(5,3)=borncache145(1,1)
      borncache(6,3)=borncache145(1,2)
      borncache(7,3)=borncache145(2,1)
      borncache(8,3)=borncache145(2,2)
      borncache(9,3)=borncache145(3,1)

      endif

! Kinematics 5

      if (lokt(6)) then !.and.dotrr(p(0,1,1),p(0,3,1)).le.1d0) then   ! Dipole 1,34

      call getyourscalesready_trib(xi,p(0,1,6),v(0,1,6),pdfchange,pdftmp,nlo,6)

      Call daisjdipole67_trib(6,1,1,xuz(1,6),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipolea(1),ps_number)
      Call daisjdipole67_trib(6,2,1,xuz(1,6),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipoleb(1),ps_number)
      Call daisjdipole65_trib(6,1,1,xuz(1,6),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipolev(1,1),ps_number)
      Call daisjdipole65_trib(6,1,2,xuz(1,6),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipolev(1,2),ps_number)

      if((Loops_sub_NLO.ne.1).and.(sub_number.eq.1)) then
      dipole(1)=zero
      dipole(2)=zero
      else
      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)
      endif
      if((Loops_sub_NLO.ne.1).and.(sub_number.eq.2)) then
      dipole(4)=zero
      dipole(6)=zero
      else
      dipole(4)=dipolea(2)
      dipole(6)=dipoleb(2)
      endif

      Call daisjdipole67_trib(6,1,2,xuz(1,6),p(0,3,1), p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipolea(1),ps_number)
      Call daisjdipole67_trib(6,2,2,xuz(1,6),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipoleb(1),ps_number)
      Call daisjdipole65_trib(6,3,1,xuz(1,6),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipolev(1,3),ps_number)

      dipole(3)=dipolea(3)  
      dipole(5)=dipoleb(3)
 
      kinem(5)= zero 
     &        + dipole(2)*inv2*(pdf( 2,1)*pdf( 1,2)+pdf( 4,1)*pdf( 3,2))    !u d  > u u   !pdf check !6-7
     &        + dipole(1)*inv2*(pdf(-1,1)*pdf(-2,2)+pdf(-3,1)*pdf(-4,2))    !d~u~ > d~d~  !pdf check??
! !-------------------------------------------------------------------------------------------------------------
     &        + dipolev(1,1)*inv2*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))   !u~d  > g g , c~s > g g    !pdf check
     &        + dipolev(1,2)*inv2*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                              !pdf check
! !-------------------------------------------------------------------------------------------------------------
     &        + dipole(3)*(pdf(-2,1)*pdf( 0,2)+pdf(-4,1)*pdf( 0,2))      !u~g  > d g , c~g > s g    !pdf check   41-42
     &        + dipole(4)*(pdf( 0,1)*pdf(-2,2)+pdf( 0,1)*pdf(-4,2))                                 !pdf check
! !-------------------------------------------------------------------------------------------------------------
     &        + dipole(5)*(pdf( 1,1)*pdf( 0,2)+pdf( 3,1)*pdf( 0,2))      !d g  > u g              !pdf check 43-44
     &        + dipole(6)*(pdf( 0,1)*pdf( 1,2)+pdf( 0,1)*pdf( 3,2))                         !pdf check 
! !-------------------------------------------------------------------------------------------------------------
     &        + 2d0*dipolev(5,3) *(pdf(0,1)*pdf( 0,2))                       !g g > u d~              !pdf check 

      endif

! Kinematics 6

      if (lokt(7)) then !.and.dotrr(p(0,2,1),p(0,3,1)).le.1d0) then   ! Dipole 2,34

      call getyourscalesready_trib(xi,p(0,1,7),v(0,1,7),pdfchange,pdftmp,nlo,7)

      Call daisjdipole67_trib(7,1,1,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipolea(1),ps_number)
      Call daisjdipole67_trib(7,2,1,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipoleb(1),ps_number)
      Call daisjdipole65_trib(7,1,2,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipolev(1,1),ps_number)
      Call daisjdipole65_trib(7,1,1,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipolev(1,2),ps_number)

      if((Loops_sub_NLO.ne.1).and.(sub_number.eq.1)) then
      dipole(1)=zero
      dipole(2)=zero
      else
      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)
      endif
      if((Loops_sub_NLO.ne.1).and.(sub_number.eq.2)) then
      dipole(3)=zero
      dipole(5)=zero
      else
      dipole(3)=dipolea(2)
      dipole(5)=dipoleb(2)
      endif

      Call daisjdipole67_trib(7,1,2,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipolea(1),ps_number)
      Call daisjdipole67_trib(7,2,2,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipoleb(1),ps_number)

      Call daisjdipole65_trib(7,3,1,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipolev(1,3),ps_number)


      dipole(4)=dipolea(3)
      dipole(6)=dipoleb(3)

      kinem(6)= zero 
     &        + dipole(2)*inv2*(pdf( 1,1)*pdf( 2,2)+pdf( 3,1)*pdf( 4,2))   !u d  > u u    !pdf check !6-7
     &        + dipole(1)*inv2*(pdf(-2,1)*pdf(-1,2)+pdf(-4,1)*pdf(-3,2))   !d~u~ > d~d~   !pdf check
! !-------------------------------------------------------------------------------------------------------------
     &        + dipolev(1,1)*inv2*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))  !u~d  > g g , c~s > g g    !pdf check
     &        + dipolev(1,2)*inv2*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                             !pdf check
! !-------------------------------------------------------------------------------------------------------------
     &        + dipole(3)*(pdf(-2,1)*pdf( 0,2)+pdf(-4,1)*pdf( 0,2)) !u~g  > d~g , c~g > s g    !pdf check   41-42
     &        + dipole(4)*(pdf( 0,1)*pdf(-2,2)+pdf( 0,1)*pdf(-4,2)) !pdf check
! !-------------------------------------------------------------------------------------------------------------
     &        + dipole(5)*(pdf( 1,1)*pdf( 0,2)+pdf( 3,1)*pdf( 0,2))   !d g  > u g              !pdf check 43-44
     &        + dipole(6)*(pdf( 0,1)*pdf( 1,2)+pdf( 0,1)*pdf( 3,2))                            !pdf check
! !-------------------------------------------------------------------------------------------------------------
     &        + 2d0*dipolev(5,3) *(pdf(0,1)*pdf( 0,2))                 !g g > u d~              !pdf check 


      endif

! Kinematics 7

      if (lokt(8)) then !.and.dotrr(p(0,1,1),p(0,3,1)).le.1d0) then   ! Dipole 12,3

      call getyourscalesready_trib(xi,p(0,1,8),v(0,1,8),pdfchange,pdftmp,nlo,8)

      Call daibdipole147_trib(8,1,1,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipolea(1),ps_number)
      Call daibdipole147_trib(8,2,1,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipoleb(1),ps_number)
      Call daibdipole145_trib(8,1,1,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipolev(1,1),ps_number)
      Call daibdipole145_trib(8,1,2,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipolev(1,2),ps_number)

      if((Loops_sub_NLO.ne.1).and.(sub_number.eq.1)) then
      dipole(1)=zero
      dipole(2)=zero
      else
      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)
      endif
      if((Loops_sub_NLO.ne.1).and.(sub_number.eq.2)) then
      dipole(4)=zero
      dipole(6)=zero
      else
      dipole(4)=dipolea(2)
      dipole(6)=dipoleb(2)
      endif

      Call daibdipole147_trib(8,1,2,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipolea(1),ps_number)
      Call daibdipole147_trib(8,2,2,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipoleb(1),ps_number)
      Call daibdipole145_trib(8,3,1,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipolev(1,3),ps_number)

      dipole(3)=dipolea(3)
      dipole(5)=dipoleb(3)

      kinem(7)= zero 
     &        + dipole(2)*inv2*(pdf( 2,1)*pdf( 1,2)+pdf( 4,1)*pdf( 3,2))    !u d  > u u  !pdf check !6-7
     &        + dipole(1)*inv2*(pdf(-1,1)*pdf(-2,2)+pdf(-3,1)*pdf(-4,2))    !d~u~ > d~d~ !pdf check??
! !-------------------------------------------------------------------------------------------------------------
     &        + dipolev(1,1)*inv2*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))   !u~d  > g g  !pdf check
     &        + dipolev(1,2)*inv2*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                !pdf check
! !-------------------------------------------------------------------------------------------------------------
     &        + dipole(3)*(pdf(-2,1)*pdf( 0,2)+pdf(-4,1)*pdf( 0,2))     !u~g  > d~g , c~g > s g   !pdf check   41-42
     &        + dipole(4)*(pdf( 0,1)*pdf(-2,2)+pdf( 0,1)*pdf(-4,2))                               !pdf check
! !-------------------------------------------------------------------------------------------------------------
     &        + dipole(5)*(pdf( 1,1)*pdf( 0,2)+pdf( 3,1)*pdf( 0,2))     !d g  > u g              !pdf check 43-44
     &        + dipole(6)*(pdf( 0,1)*pdf( 1,2)+pdf( 0,1)*pdf( 3,2))                        !pdf check XX
! !-------------------------------------------------------------------------------------------------------------
     &        + 2d0*dipolev(5,3) *(pdf(0,1)*pdf( 0,2))


!
!     Fill Cache
!
      borncache(1,2)=borncache147(1,1)
      borncache(2,2)=borncache147(1,2)
      borncache(3,2)=borncache147(2,1)
      borncache(4,2)=borncache147(2,2)
      borncache(5,2)=borncache145(1,1)
      borncache(6,2)=borncache145(1,2)
      borncache(7,2)=borncache145(2,1)
      borncache(8,2)=borncache145(2,2)
      borncache(9,2)=borncache145(3,2)

      endif

! Kinematics 8

      if (lokt(9)) then !.and.dotrr(p(0,2,1),p(0,3,1)).le.1d0) then   ! Dipole 21,3

      call getyourscalesready_trib(xi,p(0,1,9),v(0,1,9),pdfchange,pdftmp,nlo,9)

      Call daibdipole147_trib(9,1,1,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipolea(1),ps_number)
      Call daibdipole147_trib(9,2,1,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipoleb(1),ps_number)
      Call daibdipole145_trib(9,1,2,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipolev(1,1),ps_number)
      Call daibdipole145_trib(9,1,1,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipolev(1,2),ps_number)

      if((Loops_sub_NLO.ne.1).and.(sub_number.eq.1)) then
      dipole(1)=zero
      dipole(2)=zero
      else
      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)
      endif
      if((Loops_sub_NLO.ne.1).and.(sub_number.eq.2)) then
      dipole(3)=zero
      dipole(5)=zero
      else
      dipole(3)=dipolea(2)
      dipole(5)=dipoleb(2)
      endif

      Call daibdipole147_trib(9,1,2,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipolea(1),ps_number)
      Call daibdipole147_trib(9,2,2,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipoleb(1),ps_number)
      Call daibdipole145_trib(9,3,1,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipolev(1,3),ps_number)

      dipole(4)=dipolea(3)
      dipole(6)=dipoleb(3)

      kinem(8)= zero 
     &        + dipole(2)*inv2*(pdf( 1,1)*pdf( 2,2)+pdf( 3,1)*pdf( 4,2))     !u d  > u u  !pdf check !6-7
     &        + dipole(1)*inv2*(pdf(-2,1)*pdf(-1,2)+pdf(-4,1)*pdf(-3,2))     !d~u~ > d~d~ !pdf check
! !-------------------------------------------------------------------------------------------------------------
     &        + dipolev(1,1)*inv2*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))  !u~d  > g g  !pdf check
     &        + dipolev(1,2)*inv2*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))               !pdf check
! !-------------------------------------------------------------------------------------------------------------
     &        + dipole(3)*(pdf(-2,1)*pdf( 0,2)+pdf(-4,1)*pdf( 0,2)) !u~g  > d~g , c~g > s g   !pdf check  41-42
     &        + dipole(4)*(pdf( 0,1)*pdf(-2,2)+pdf( 0,1)*pdf(-4,2))                           !pdf check
! !-------------------------------------------------------------------------------------------------------------
     &        + dipole(5)*(pdf( 1,1)*pdf( 0,2)+pdf( 3,1)*pdf( 0,2))   !d g  > u g               !pdf check 43-44
     &        + dipole(6)*(pdf( 0,1)*pdf( 1,2)+pdf( 0,1)*pdf( 3,2))                             !pdf check
! !-------------------------------------------------------------------------------------------------------------
     &        + 2d0*dipolev(5,3) *(pdf(0,1)*pdf( 0,2))                 !g g > u d~              !pdf check 


!
!     Fill Cache
!
      borncache(1,4)=borncache147(1,1)
      borncache(2,4)=borncache147(1,2)
      borncache(3,4)=borncache147(2,1)
      borncache(4,4)=borncache147(2,2)
      borncache(5,4)=borncache145(1,1)
      borncache(6,4)=borncache145(1,2)
      borncache(7,4)=borncache145(2,1)
      borncache(8,4)=borncache145(2,2)
      borncache(9,4)=borncache145(3,2)

      endif


! Kinematics 9

      if (lokt(10)) then !.and.dotrr(p(0,4,1),p(0,3,1)).le.1d0) then   ! Dipole 1,34

      call getyourscalesready_trib(xi,p(0,1,10),v(0,1,10),pdfchange,pdftmp,nlo,10)

      CALL dasijdipole40_trib(10,1,1,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1),ps_number)
      CALL dasijdipole40_trib(10,2,1,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1),ps_number)

      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)

!       CALL dasijdipole40_trib(1,2,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1),ps_number)
!       CALL dasijdipole40_trib(2,2,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1),ps_number)
! 
!       dipole(3)=dipolea(2)
!       dipole(4)=dipoleb(1)

      CALL test40_trib(10,1,2,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1),ps_number)
      CALL test40_trib(10,2,2,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1),ps_number)

      dipole(3)=dipolea(2)
      dipole(4)=dipoleb(1)

!       CALL dasijdipole40_trib(1,3,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1),ps_number)
!       CALL dasijdipole40_trib(2,3,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1),ps_number)

      CALL test40_trib(10,1,3,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1),ps_number)
      CALL test40_trib(10,2,3,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1),ps_number)

      dipole(5)=dipolea(2)
      dipole(6)=dipoleb(1)

      CALL dasijdipole40_trib(10,1,4,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1),ps_number)
      CALL dasijdipole40_trib(10,2,4,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1),ps_number)

      dipole(7)=dipolea(1)
      dipole(8)=dipoleb(1)


      kinem(9)= zero
     &        + dipole(7)*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))   !u~d  > u u~    !pdf check
     &        + dipole(8)*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                   !pdf check
     &        + dipole(7)*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))   !u~d  > d d~    !pdf check
     &        + dipole(8)*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                   !pdf check
!!!!!!
     &        + dipole(7)*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))   !u~d  > s s~    !pdf check
     &        + dipole(8)*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                   !pdf check
     &        + dipole(7)*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))  !u~d  > c c~    !pdf check
     &        + dipole(8)*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                  !pdf check
! !-------------------------------------------------------------------------------------------------------------
     &        + dipole(1)*inv2*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))  !u~d  > g g     !pdf check
     &        + dipole(2)*inv2*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                  !pdf check
! !-------------------------------------------------------------------------------------------------------------
     &        + dipole(3)*(pdf(-2,1)*pdf( 0,2)+pdf(-4,1)*pdf( 0,2))  !u~g  > d~g , c~g > s g   !pdf check  41-42
     &        + dipole(4)*(pdf( 0,1)*pdf(-2,2)+pdf( 0,1)*pdf(-4,2))                            !pdf check  
! !-------------------------------------------------------------------------------------------------------------
     &        + dipole(5)*(pdf( 1,1)*pdf( 0,2)+pdf( 3,1)*pdf( 0,2))  !d g  > u g              !pdf check 43-44
     &        + dipole(6)*(pdf( 0,1)*pdf( 1,2)+pdf( 0,1)*pdf( 3,2))                           !pdf check


      endif


! Kinematics 10

      if (lokt(11)) then !.and.dotrr(p(0,4,1),p(0,3,1)).le.1d0) then   ! Dipole 2,34

      call getyourscalesready_trib(xi,p(0,1,11),v(0,1,11),pdfchange,pdftmp,nlo,11)

      CALL dasijdipole40_trib(11,1,1,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipolea(1),ps_number)
      CALL dasijdipole40_trib(11,2,1,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipoleb(1),ps_number)

      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)

!       CALL dasijdipole40_trib(1,2,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipolea(1),ps_number)
!       CALL dasijdipole40_trib(2,2,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipoleb(1),ps_number)

      CALL test40_trib(11,1,2,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipolea(1),ps_number)
      CALL test40_trib(11,2,2,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipoleb(1),ps_number)
      dipole(3)=dipolea(1)
      dipole(4)=dipoleb(2)

!       CALL dasijdipole40_trib(1,3,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipolea(1),ps_number)
!       CALL dasijdipole40_trib(2,3,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipoleb(1),ps_number)

      CALL test40_trib(11,1,3,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipolea(1),ps_number)
      CALL test40_trib(11,2,3,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipoleb(1),ps_number)
! 
      dipole(5)=dipolea(1)
      dipole(6)=dipoleb(2)

      CALL dasijdipole40_trib(11,1,4,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipolea(1),ps_number)
      CALL dasijdipole40_trib(11,2,4,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipoleb(1),ps_number)

      dipole(7)=dipolea(1)
      dipole(8)=dipoleb(1)

      kinem(10)= zero 
     &        + dipole(7)*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))   !u~d  > u u~  !pdf check
     &        + dipole(8)*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                 !pdf check
     &        + dipole(7)*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))   !u~d  > d d~  !pdf check
     &        + dipole(8)*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                 !pdf check
!!!!!
     &        + dipole(7)*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))   !u~d  > s s~  !pdf check
     &        + dipole(8)*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                 !pdf check
     &        + dipole(7)*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))   !u~d  > c c~  !pdf check
     &        + dipole(8)*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                 !pdf check 
! !-------------------------------------------------------------------------------------------------------------
     &        + dipole(1)*inv2*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))  !u~d  > g g   !pdf check
     &        + dipole(2)*inv2*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                !pdf check
! !-------------------------------------------------------------------------------------------------------------
     &         + dipole(3)*(pdf(-2,1)*pdf( 0,2)+pdf(-4,1)*pdf( 0,2))    !u~g  > d g, c~g  > s g  !pdf check   41-42
     &         + dipole(4)*(pdf( 0,1)*pdf(-2,2)+pdf( 0,1)*pdf(-4,2))                             !pdf check 
! !-------------------------------------------------------------------------------------------------------------
     &         + dipole(5)*(pdf( 1,1)*pdf( 0,2)+pdf( 3,1)*pdf( 0,2))    !d g  > u g              !pdf check 43-44
     &         + dipole(6)*(pdf( 0,1)*pdf( 1,2)+pdf( 0,1)*pdf( 3,2))                             !pdf check

      endif

!       if (m2s(0).ne.0d0) then
!       write(66,*), dotrr(p(0,1,1),p(0,4,1)), m2s(1),m2s(0)-m2s(1)
!       endif


      fincolconv=two*pi/dotrr(p(0,1,1),p(0,2,1))

      if (pdfchange) then
      signpdf1=-1
      signpdf2=-1
      else
      signpdf1=1
      signpdf2=1
      endif
      signpdf1=signpdf1*sign1
      signpdf2=signpdf2*sign2

!
!  Finally calculate finite collinear terms
!
      if((Loops_sub_NLO.eq.1).or.(sub_number.eq.1)) then
      if (lokt(4).and.fincolcalc) then !this is one of the initial-initial kinematics, e.g. eq 5.147 in Catani-Seymour
!                                                                                       first splits
      call getyourscalesready_trib(xi,p(0,1,4),v(0,1,4),pdfchange,pdf,nlo,4)
      q_sf=dsqrt(mufsq(1,4))
 
      CALL finitecollqqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,-2,fincol(1),signpdf1)  ! ub d
!                                                  =p(0,2,1)
      CALL finitecollqqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,1,fincol(2),signpdf1)   ! d  ub
!                                                  =p(0,2,1)
      CALL finitecollqgREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,-2,fincol(3),signpdf1)  ! ub g
!                                                  =p(0,2,1)
      CALL finitecollgqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,fincol(4),signpdf1)     ! g  ub 
!                                                  =p(0,2,1)
!       CALL finitecollqgREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,1,fincol(5)) ! d  g
! !                                                  =p(0,2,1)
!       CALL finitecollgqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,fincol(6))   ! g  d 
! !                                                  =p(0,2,1)

      CALL finitecollqqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,-4,fincol(7),signpdf1)  ! cb s
!                                                  =p(0,2,1)
      CALL finitecollqqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,3,fincol(8),signpdf1)   ! s  cb
!                                                  =p(0,2,1)
      CALL finitecollqgREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,-4,fincol(9),signpdf1)  ! cb g
!                                                  =p(0,2,1)
!      CALL finitecollgqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,fincol(10))   ! g  cb 
!                                                  =p(0,2,1)
      fincol(10)=fincol(4)
!       CALL finitecollqgREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,3,fincol(11))  ! s  g
! !                                                  =p(0,2,1)
! !      CALL finitecollgqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,fincol(12))   ! g  s 
! !                                                  =p(0,2,1)
!       fincol(12)=fincol(6)


!        do k=0,3
!           helastestmom(k,1)=p(k,1,4)
!           helastestmom(k,2)=p(k,2,4)
!           helastestmom(k,3)=p(k,3,4)
!           helastestmom(k,4)=v(k,3,4)
!           helastestmom(k,5)=v(k,1,4)
!           helastestmom(k,6)=v(k,2,4)
!        enddo
! 
!       CAll SUDGENUAA(helastestmom,fincolampsq(1))
      fincolampsq(1)=borncache(5,1) 
!       Call SUGDENUAA(helastestmom,fincolampsq(3))
      fincolampsq(3)=borncache(7,1) 

!       do k=0,3
!          helastestmom(k,1)=p(k,2,4)
!          helastestmom(k,2)=p(k,1,4)
!          helastestmom(k,3)=p(k,3,4)
!          helastestmom(k,4)=v(k,3,4)
!          helastestmom(k,5)=v(k,1,4)
!          helastestmom(k,6)=v(k,2,4)
!       enddo

!       CAll SUDGENUAA(helastestmom,fincolampsq(2))
      fincolampsq(2)=borncache(6,1)
!       Call SUGDENUAA(helastestmom,fincolampsq(4))
      fincolampsq(4)=borncache(1,1)

      amp=(fincol(1)*pdf(1,2)+fincol(7)*pdf(3,2))   *fincolampsq(1)
     &   +(fincol(2)*pdf(-2,2)+fincol(8)*pdf(-4,2)) *fincolampsq(2)
     &   +(fincol(3)*pdf(0,2)+fincol(9)*pdf(0,2))   *fincolampsq(3)
     &   +(fincol(4)*pdf(-2,2)+fincol(10)*pdf(-4,2))*fincolampsq(4)


!       amp=2d0/(1d0-xuz(1,4))*(xuz(1,4)/(1d0-xi(1)*xuz(1,4)))*fincolampsq(1)!*pdf(-2,1)*pdf(1,2)
!      &    *8.*pi**2/2./dotrr(p(0,1,4),p(0,2,4))

!       help1=xuz(1,4)*xi(1)
!       call pdfproton(help1,q_sf,pdftmp(-6,1))
!       do mu=-6,6
!       pdf(mu,1)=pdftmp(mu,1)/help1
!       enddo
!       amp=1d0/(1d0-xuz(1,4))*(xuz(1,4)/(1d0-xi(1)*xuz(1,4)))*fincolampsq(1)
!       amp=1d0/(1d0-xuz(1,4))*(xuz(1,4)
!      &    /(1d0-xi(1)*xuz(1,4)))*pdf(-2,1)*pdf(1,2)*fincolampsq(1)
!       coll3=amp*8d0*pi**2/2d0/dotrr(p(0,1,4),p(0,2,4))*2d0
!       print*, alfas, gg(1), q_sf

      coll3=2d0*amp*fincolconv*alfas/xuz(1,4)

!       do mu=1,4
!       print*,fincolampsq(mu) , mu
!       enddo
!       print*,"-------------"
!       do mu=1,4
!       print*,fincol(mu) , mu
!       enddo
!       print*,"-------------"
!       STOP
!       do mu=1,9
!       print*,borncache(mu,1) , mu
!       enddo
!       print*,"-------------"
! 
!       STOP


      endif


      if (lokt(8).and.fincolcalc) then

      call getyourscalesready_trib(xi,p(0,1,8),v(0,1,8),pdfchange,pdf,nlo,8)
      q_sf=dsqrt(mufsq(1,8))

      CALL finitecollqgREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,1,fincol(5),signpdf1)   ! d  g
!                                                  =p(0,2,1)
      CALL finitecollgqREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,fincol(6),signpdf1)     ! g  d 
!                                                  =p(0,2,1)
      CALL finitecollqgREM(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,3,fincol(11),signpdf1)  ! s  g
!                                                  =p(0,2,1)
!      CALL finitecollgqREM(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,fincol(12))   ! g  s 
!                                                  =p(0,2,1)
      fincol(12)=fincol(6)

!       do k=0,3
! !          helastestmom(k,1)=p(k,1,8)
!          helastestmom(k,2)=p(k,2,8)
!          helastestmom(k,3)=p(k,3,8)
!          helastestmom(k,4)=v(k,3,8)
!          helastestmom(k,5)=v(k,1,8)
!          helastestmom(k,6)=v(k,2,8)
!       enddo
! 
!       Call SDGUENUAA(helastestmom,fincolampsq(5))
      fincolampsq(5)=borncache(4,2)

!       do k=0,3
!          helastestmom(k,1)=p(k,2,8)
!          helastestmom(k,2)=p(k,1,8)
!          helastestmom(k,3)=p(k,3,8)
!          helastestmom(k,4)=v(k,3,8)
!          helastestmom(k,5)=v(k,1,8)
!          helastestmom(k,6)=v(k,2,8)
!       enddo
! 
!       Call SDGUENUAA(helastestmom,fincolampsq(6))
       fincolampsq(6)=borncache(3,2)

!       do mu=5,6
!       print*,fincolampsq(mu) , mu
!       enddo
!       print*,"-------------"
!       do mu=1,9
!       print*,borncache(mu,2) , mu
!       enddo
!       print*,"-------------"
!       print*, fincolampsq(5)/borncache(4,2)
!       print*, fincolampsq(6)/borncache(3,2)
!       STOP

      amp=(fincol(5)*pdf(0,2)+fincol(11)*pdf(0,2))  *fincolampsq(5)!
     &   +(fincol(6)*pdf(1,2)+fincol(12)*pdf(3,2))  *fincolampsq(6)

      coll7=amp*2d0*fincolconv*alfas/xuz(1,8)

!       help1=xuz(1,8)*xi(1)
!       call pdfproton(help1,q_sf,pdftmp(-6,1))
!       do mu=-6,6
!       pdf(mu,1)=pdftmp(mu,1)/help1
!       enddo
! 
! 
!       amp=1d0/(1d0-xuz(1,8))*(xuz(1,8)
!      &    /(1d0-xi(1)*xuz(1,8)))*fincolampsq(5)*pdf(1,1)*pdf(0,2)

!       amp=(fincol(5)*pdf(0,2))*fincolampsq(5)*alfas/2./pi
!       coll7=amp*8d0*pi**2/2d0/dotrr(p(0,1,8),p(0,2,8))*2d0

      endif


!      
      if (lokt(5).and.fincolcalc) then  ! this is one of the initial-initial kinematics, e.g. eq 5.147 in Catani-Seymour
!
      call getyourscalesready_trib(xi,p(0,1,5),v(0,1,5),pdfchange,pdf,nlo,5)
      q_sf=dsqrt(mufsq(1,5))
!                                                                                  first splits
      CALL finitecollqqREM(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,5),q_sf,-2,fincol(1),signpdf2)  ! ub d
!                                                  =p(0,2,1)
      CALL finitecollqqREM(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,5),q_sf,1,fincol(2),signpdf2)   ! d  ub
!                                                  =p(0,2,1)
      CALL finitecollqgREM(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,5),q_sf,-2,fincol(3),signpdf2)  ! ub g
!                                                  =p(0,2,1)
      CALL finitecollgqREM(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,5),q_sf,fincol(4),signpdf2)     ! g  ub 
!                                                  =p(0,2,1)
!       CALL finitecollqgREM(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,5),q_sf,1,fincol(5))   ! d  g
! !                                                  =p(0,2,1)
!       CALL finitecollgqREM(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,5),q_sf,fincol(6))     ! g  d 
! !                                                  =p(0,2,1)
      CALL finitecollqqREM(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,5),q_sf,-4,fincol(7),signpdf2)  ! cb s
!                                                  =p(0,2,1)
      CALL finitecollqqREM(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,5),q_sf,3,fincol(8),signpdf2)   ! s  cb
!                                                  =p(0,2,1)
      CALL finitecollqgREM(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,5),q_sf,-4,fincol(9),signpdf2)  ! cb g
!                                                  =p(0,2,1)
!      CALL finitecollgqREM(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,5),q_sf,fincol(10))     ! g  cb 
!                                                  =p(0,2,1)
      fincol(10)=fincol(4)
!       CALL finitecollqgREM(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,5),q_sf,3,fincol(11))   ! s  g
! !                                                  =p(0,2,1)
! !      CALL finitecollgqREM(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,4),q_sf,fincol(12))     ! g  s 
! !                                                  =p(0,2,1)
!       fincol(12)=fincol(6)


      do k=0,3
         helastestmom(k,1)=p(k,1,5)
         helastestmom(k,2)=p(k,2,5)
         helastestmom(k,3)=p(k,3,5)
         helastestmom(k,4)=v(k,3,5)
         helastestmom(k,5)=v(k,1,5)
         helastestmom(k,6)=v(k,2,5)
      enddo

!       CAll SUDGENUAA(helastestmom,fincolampsq(1))
      fincolampsq(1)=borncache(5,3)
!       Call SUGDENUAA(helastestmom,fincolampsq(3))
      fincolampsq(3)=borncache(7,3)


!       do k=0,3
!          helastestmom(k,1)=p(k,2,5)
!          helastestmom(k,2)=p(k,1,5)
!          helastestmom(k,3)=p(k,3,5)
!          helastestmom(k,4)=v(k,3,5)
!          helastestmom(k,5)=v(k,1,5)
!          helastestmom(k,6)=v(k,2,5)
!       enddo
! 
!       CAll SUDGENUAA(helastestmom,fincolampsq(2))
      fincolampsq(2)=borncache(6,3)
!       Call SUGDENUAA(helastestmom,fincolampsq(4))
      fincolampsq(4)=borncache(1,3)
!       Call SDGUENUAA(helastestmom,fincolampsq(6))
      

!       do mu=1,4
!       print*,fincolampsq(mu) , mu
!       enddo
!       print*,"-------------"
!       do mu=1,9
!       print*,borncache(mu,3) , mu
!       enddo
!       print*,"-------------"
!       print*,fincolampsq(1)/borncache(5,3)
!       print*,fincolampsq(2)/borncache(6,3)
!       print*,fincolampsq(3)/borncache(7,3)
!       print*,fincolampsq(4)/borncache(1,3)
!       STOP

      amp=(fincol(1)*pdf( 1,1)+fincol(7)*pdf( 3,1)) *fincolampsq(1)
     &   +(fincol(2)*pdf(-2,1)+fincol(8)*pdf(-4,1)) *fincolampsq(2)
     &   +(fincol(3)*pdf(0,1)+fincol(9)*pdf(0,1))   *fincolampsq(3)
     &   +(fincol(4)*pdf(-2,1)+fincol(10)*pdf(-4,1))*fincolampsq(4)

!       help1=xuz(1,5)*xi(2)
!       call pdfproton(help1,q_sf,pdftmp(-6,2))
!       do mu=-6,6
!       pdf(mu,2)=pdftmp(mu,2)/help1
!       enddo
!       amp=1d0/(1d0-xuz(1,5))*(xuz(1,5)/(1d0-xi(2)*xuz(1,5)))*fincolampsq(1)*pdf(-2,2)*pdf(1,1) !(fincol(1)*pdf(1,1))   

      coll4=amp*2d0*fincolconv*alfas/xuz(1,5)

      endif


      if (lokt(9).and.fincolcalc) then

      call getyourscalesready_trib(xi,p(0,1,9),v(0,1,9),pdfchange,pdf,nlo,9)
      q_sf=dsqrt(mufsq(1,9))

      CALL finitecollqgREM(xi(2),xuz(1,9),p(0,1,9),p(0,2,9),p(0,3,9),q_sf,1,fincol(5),signpdf2)   ! d  g
!                                                  =p(0,2,1)
      CALL finitecollgqREM(xi(2),xuz(1,9),p(0,1,9),p(0,2,9),p(0,3,9),q_sf,fincol(6),signpdf2)     ! g  d 
!                                                  =p(0,2,1)
      CALL finitecollqgREM(xi(2),xuz(1,9),p(0,1,9),p(0,2,9),p(0,3,9),q_sf,3,fincol(11),signpdf2)  ! s  g
!                                                  =p(0,2,1)
!      CALL finitecollgqREM(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,4),q_sf,fincol(12))   ! g  s 
!                                                  =p(0,2,1)
      fincol(12)=fincol(6)


!       do k=0,3
!          helastestmom(k,1)=p(k,1,9)
!          helastestmom(k,2)=p(k,2,9)
!          helastestmom(k,3)=p(k,3,9)
!          helastestmom(k,4)=v(k,3,9)
!          helastestmom(k,5)=v(k,1,9)
!          helastestmom(k,6)=v(k,2,9)
!       enddo
! 
!       Call SDGUENUAA(helastestmom,fincolampsq(5))
      fincolampsq(5)=borncache(4,4)


!       do k=0,3
!          helastestmom(k,1)=p(k,2,9)
!          helastestmom(k,2)=p(k,1,9)
!          helastestmom(k,3)=p(k,3,9)
!          helastestmom(k,4)=v(k,3,9)
!          helastestmom(k,5)=v(k,1,9)
!          helastestmom(k,6)=v(k,2,9)
!       enddo
! 
!       Call SDGUENUAA(helastestmom,fincolampsq(6))
      fincolampsq(6)=borncache(3,4)

!       do mu=5,6
!       print*,fincolampsq(mu) , mu
!       enddo
!       print*,"-------------"
!       do mu=1,9
!       print*,borncache(mu,4) , mu
!       enddo
!       print*,"-------------"
!       print*,fincolampsq(5)/borncache(4,4)
!       print*,fincolampsq(6)/borncache(3,4)
!       STOP

      amp=(fincol(5)*pdf(0,1)+fincol(11)*pdf(0,1))  *fincolampsq(5)
     &   +(fincol(6)*pdf(1,1)+fincol(12)*pdf(3,1))  *fincolampsq(6)

      coll8=amp*2d0*fincolconv*alfas/xuz(1,9)

      endif
      endif

      m2s(1) = m2s_qqwgaglujfey
      m2s(2) =-kinem(1)
      m2s(3) =-kinem(2)
      m2s(4) =+coll3-kinem(3)
      m2s(5) =+coll4-kinem(4)
      m2s(6) =-kinem(5)
      m2s(7) =-kinem(6)
      m2s(8) =+coll7-kinem(7)
      m2s(9) =+coll8-kinem(8)
      m2s(10)=-kinem(9)
      m2s(11)=-kinem(10)

!  Sum the amplitude
!
      m2s(0)=m2s(1)+m2s(2)+m2s(3)+m2s(4)+m2s(5)+m2s(6)+m2s(7)+m2s(8)+m2s(9)+m2s(10)+m2s(11)
      do i=0,11
      m2s(i)=m2s(i)*4d0 ! random photon helicity summation
      enddo

!      print*, m2s
!      print*, "------"
!      print*, gg,gw,gwwa
!      STOP

!       if (.false.) then
! 
!       do mu=0,3
!         do i=1,11
!          ptelwtmp(mu,i)=v(mu,1,i)+v(mu,2,i)+v(mu,3,i)
!         enddo
!       enddo
!       do i=1,11
!       ptelw(i)=dsqrt(ptelwtmp(1,i)**2+ptelwtmp(2,i)**2)
!       enddo
! 
!       if(abs(maxeventweight).le.abs(m2s(0))) then
!       maxeventweight=m2s(0)
!       print*, "new maxweight=", m2s(0), " ptelw=", ptelw(1), " lokt=", lokt(1)
!       print*, "p2p4=",dotrr(p(0,2,1),p(0,4,1))
!       print*, "p1p4=",dotrr(p(0,1,1),p(0,4,1))
!       print*, "p2p3=",dotrr(p(0,2,1),p(0,3,1))
!       print*, "p1p3=",dotrr(p(0,1,1),p(0,3,1))
!       print*, "p1p2=",dotrr(p(0,1,1),p(0,2,1))
!       print*, "p3p4=",dotrr(p(0,3,1),p(0,4,1))
!       print*, "p1v3=",dotrr(p(0,1,1),v(0,3,1))
!       print*, "p2v3=",dotrr(p(0,2,1),v(0,3,1))
!       print*, "p3v3=",dotrr(p(0,3,1),v(0,3,1))
!       print*, "p4v3=",dotrr(p(0,4,1),v(0,3,1))
!       print*, "pt photon =", dsqrt(v(1,3,1)**2+v(2,3,1)**2)
!       print*, "pt lepton =", dsqrt(v(1,1,1)**2+v(2,1,1)**2)
!       print*, "--------"
!       print*, "lokt 2", lokt(2), " lokt 4", lokt(4), " Dipole p1p4 ", m2s(1)/(m2s(2)+m2s(4)), partons(5,1,2), partons(5,1,4), ptelw(2), ptelw(4),  partons(6,1,2), partons(6,1,4)
!       print*, "lokt 3", lokt(3), " lokt 5", lokt(5), " Dipole p2p4 ", m2s(1)/(m2s(3)+m2s(5)), partons(5,1,3), partons(5,1,5), ptelw(3), ptelw(5),  partons(6,1,3), partons(6,1,5)
!       print*, "lokt 6", lokt(6), " lokt 8", lokt(8), " Dipole p1p3 ", m2s(1)/(m2s(6)+m2s(8)), partons(5,1,6), partons(5,1,8), ptelw(6), ptelw(8),  partons(6,1,6), partons(6,1,8)
!       print*, "lokt 7", lokt(7), " lokt 9", lokt(9), " Dipole p2p3 ", m2s(1)/(m2s(7)+m2s(9)), partons(5,1,7), partons(5,1,9), ptelw(7), ptelw(9),  partons(6,1,7), partons(6,1,9)
!       print*, "lokt10", lokt(10), " lokt11", lokt(11), " Dipole p3p4", m2s(1)/(m2s(10)+m2s(11)), partons(5,1,10), partons(5,1,11), ptelw(10), ptelw(11),  partons(6,1,10), partons(6,1,11)
!       print*,  "--------"
!       print*,  "--------"
!       print*,  "--------"
!       if (ptelw(1).le.40.) pause
!       endif
!       endif

!       if (m2s(0).ne.0d0) then
!       if (dotrr(p(0,1,1),p(0,4,1)).le.2d0.or.dotrr(p(0,2,1),p(0,4,1)).le.2d0)
!      & write(66,*), dotrr(p(0,1,1),p(0,4,1)),dotrr(p(0,2,1),p(0,4,1)), m2s(1),m2s(0)-m2s(1)
!       if (dotrr(p(0,1,1),p(0,3,1)).le.2d0.or.dotrr(p(0,2,1),p(0,3,1)).le.2d0)
!      & write(67,*), dotrr(p(0,1,1),p(0,3,1)),dotrr(p(0,2,1),p(0,3,1)), m2s(1),m2s(0)-m2s(1)
!       if (dotrr(p(0,3,1),p(0,4,1)).le.2d0) then
!       write(68,*), dotrr(p(0,3,1),p(0,4,1)), m2s(1),m2s(0)-m2s(1)
!       print*, dotrr(p(0,3,1),p(0,4,1))
!       endif
!       endif

      if (.false.) then
      if (abs(dotrr(p(0,2,1),p(0,4,1))).le.1e-1.or.abs(dotrr(p(0,1,1),p(0,4,1))).le.1e-1.or.
     &    abs(dotrr(p(0,2,1),p(0,3,1))).le.1e-1.or.abs(dotrr(p(0,1,1),p(0,3,1))).le.1e-1.or.
     &    abs(dotrr(p(0,3,1),p(0,4,1))).le.1e-1) then
!      if (abs(m2s(1)/(m2s(10)+m2s(11))).ge.0.9.and.abs(m2s(1)/(m2s(3)+m2s(5))).le.1.1) then
      print*, "p2p4=",dotrr(p(0,2,1),p(0,4,1))
      print*, "p1p4=",dotrr(p(0,1,1),p(0,4,1))
      print*, "p2p3=",dotrr(p(0,2,1),p(0,3,1))
      print*, "p1p3=",dotrr(p(0,1,1),p(0,3,1))
      print*, "p1p2=",dotrr(p(0,1,1),p(0,2,1))
      print*, "p3p4=",dotrr(p(0,3,1),p(0,4,1))
      print*, "p1v3=",dotrr(p(0,1,1),v(0,3,1))
      print*, "p2v3=",dotrr(p(0,2,1),v(0,3,1))
      print*, "p3v3=",dotrr(p(0,3,1),v(0,3,1))
      print*, "p4v3=",dotrr(p(0,4,1),v(0,3,1))
      print*, "drive integration via ", m2s(0), " Fract ", m2s(0)/m2s(1)
      print*, "real emission ME      ", m2s(1) , lokt(1)
      print*, "Dipole 1              ", 4d0*kinem(1) , lokt(2)
      print*, "Dipole 2              ", 4d0*kinem(2) , lokt(3)
      print*, "Dipole 3              ", 4d0*kinem(3) , lokt(4)
      print*, "Dipole 4              ", 4d0*kinem(4) , lokt(5)
      print*, "Dipole 5              ", 4d0*kinem(5) , lokt(6)
      print*, "Dipole 6              ", 4d0*kinem(6) , lokt(7)
      print*, "Dipole 7              ", 4d0*kinem(7) , lokt(8)
      print*, "Dipole 8              ", 4d0*kinem(8) , lokt(9)
      print*, "Dipole 9              ", 4d0*kinem( 9) , lokt(10)
      print*, "Dipole 10             ", 4d0*kinem(10) , lokt(11)
      print*, "REM/Dipole 1+3        ", m2s(1)/(4d0*kinem(1)+4d0*kinem(3))
      print*, "REM/Dipole 2+4        ", m2s(1)/(4d0*kinem(2)+4d0*kinem(4))
      print*, "REM/Dipole 5+7        ", m2s(1)/(4d0*kinem(5)+4d0*kinem(7))
      print*, "REM/Dipole 6+8        ", m2s(1)/(4d0*kinem(6)+4d0*kinem(8))
      print*, "REM/Dipole 9+10       ", m2s(1)/(4d0*kinem(9)+4d0*kinem(10))
      print*, "REM/Dipole tot        ", m2s(1)/(4d0*(kinem(1)+kinem(2)+
     &   kinem(3)+kinem(4)+kinem(5)+kinem(6)+kinem(7)+kinem(8)+kinem(9)+kinem(10)))
      print*, "---------------------------------"
c      read(*,*)
      endif
      endif

      if (.false.) then
        open(unit=31,file="irfinite.chk",ACCESS='APPEND')
        kintrue=.false.
        do i=1,4
          do j=i,4
            if (i.eq.j) then ! soft
              kininv(i,i) = p(0,i,1)
            else !collinear
              kininv(i,j) = sqrt(abs(dotrr(p(0,i,1),p(0,j,1))))
            endif
            if (kininv(i,j) .lt. 1d0) kintrue=.true.
          enddo
        enddo
        ratio(1) = m2s_qqwgaglujfey/(kinem( 1)+kinem( 3))
        ratio(2) = m2s_qqwgaglujfey/(kinem( 2)+kinem( 4))
        ratio(3) = m2s_qqwgaglujfey/(kinem( 5)+kinem( 7))
        ratio(4) = m2s_qqwgaglujfey/(kinem( 6)+kinem( 8))
        ratio(5) = m2s_qqwgaglujfey/(kinem( 9)+kinem(10))
c        do i=1,5
c          if ( (ratio(i).gt.0) .and. (ratio(i).lt.5) ) kintrue = .true.
c        enddo
        if (kintrue) then ! write it out
          do i=1,4
            do j=i,4
              write(31,'(g16.8)') kininv(i,j)
            enddo
          enddo
          do i=1,5
            write(31,'(g16.8)') ratio(i)
          enddo
          write(31,*) 
        endif
        close(31)
      endif

      end


      subroutine getyourscalesready_trib(xi,p,v,pdfchange,pdf,nlo,l)

      implicit none

      include 'VBFNLO/utilities/global.inc'
      include 'VBFNLO/utilities/scales.inc'
      include 'VBFNLO/utilities/coupl.inc'

      double precision q_sf,pdf(-6:6,2),pdftmp(-6:6,2),xi(nx),four,
     &  p(0:3,max_p),v(0:3,max_v)
      integer k,nlo,pdfsign
      logical pdfchange
      integer l

      parameter(four=4d0)

      q_sf = dsqrt(mufsq(1,l))

      call pdfproton( xi(1), q_sf, pdftmp(-6,1) )        ! f_a(x1)

      q_sf = dsqrt(mufsq(2,l))
      call pdfproton( xi(2), q_sf, pdftmp(-6,2) )        ! f_b(x2)

! change sign for W+ case
      if (pdfchange) then
        pdfsign = -1
      else
        pdfsign = +1
      endif

      do k=-6,6
        pdf(k,1)=pdftmp(pdfsign*sign1*k,1)/xi(1)
        pdf(k,2)=pdftmp(pdfsign*sign2*k,2)/xi(2)
      enddo


! running coupling
      alfas  =  als(1,l) 
      G      =  dsqrt(four*pi*alfas)
      GG(1)  =  -G
      GG(2)  =  -G

!      print*, "scales routine ",xi, xi2

      end


! THIS SUBROUTINE CALCULATES THE DIPOLES ACCORDING TO CATANI-SEYMOUR

      subroutine daisjdipole67_trib(n,sel,kin,xuz,pmi,pk,pa,p,v,dipoleres, 
     &     ps_number)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"

      integer ps_number

      integer mu,nu,sel,kin,n
c INPUT

c
      double precision pmi(0:3),pk(0:3),pa(0:3)
      double precision pis,CF,CA
      parameter (pis=3.141592653589793238d0,CF=4d0/3d0,CA=3d0)
      logical ldebug
      parameter (ldebug=.false.)

c OUTPUT
c
      double precision dipole(0:3,0:3),p(0:3,3),v(0:3,4)
     &    ,pipk,pkpa,pipa,xuz(2),dipoleres(3),bornmom(0:3,7),color,symmetry
c
      double complex borncurr(0:3)
      double precision Contract_Trjcj,borntest1,borntest2,metric(0:3,0:3)

      double precision xika,ui
      double precision dotrr
      external dotrr,Contract_Trjcj

      do mu=1,3
      dipoleres(mu) = 0
      enddo
      if(Loops_sub_NLO.ne.1) then
c     kin1 -> both, kin2 -> g
      if((sub_number.eq.2).and.(kin.eq.2)) return
      endif
      do mu=0,3
      do nu=0,mu-1
      metric(mu,nu)=0d0
      metric(nu,mu)=0d0
      enddo
      enddo
      metric(0,0)=1d0
      metric(1,1)=-1d0
      metric(2,2)=-1d0
      metric(3,3)=-1d0

      pipk=dotrr(pmi,pk)
      pkpa=dotrr(pk,pa)
      pipa=dotrr(pmi,pa)

      xika= xuz(1)
      ui=xuz(2)

      if (ldebug) then
      if (abs(xika-(pkpa+pipa-pipk)/(pkpa+pipa)).ge.1e-10) then
      PRINT*, "error related to xika in dipole 67"
      STOP
      endif
      if (abs(ui-pipa/(pipa+pkpa)).ge.1e-10) then
      PRINT*, "error related to ui in dipole 67"
      STOP
      endif
      do mu=0,3
      if (abs(p(mu,1)-xika*pa(mu)).ge.1e-10) then
      print*, "wrong ptilde_a(",mu," ) in 67"
      STOP
      endif
      if (abs(p(mu,3)-(pk(mu)+pmi(mu)-(1d0-xika)*pa(mu))).ge.1e-10) then
      print*, "wrong ptilde_k(",mu," ) in 67"
      STOP
      endif
      enddo
      endif

      do mu=0,3
      if (kin.eq.1) then
      bornmom(mu,1)=p(mu,2)!!!
      bornmom(mu,2)=p(mu,1)
      elseif (kin.eq.2) then
      bornmom(mu,1)=p(mu,1)!!!
      bornmom(mu,2)=p(mu,2)
      endif
      bornmom(mu,3)=p(mu,3)
      bornmom(mu,4)=v(mu,1)
      bornmom(mu,5)=v(mu,2)
      bornmom(mu,6)=v(mu,3)
      bornmom(mu,7)=v(mu,4)
      enddo

      if (sel.eq.1) then
      call ugdAAWj_curr(bornmom,borncurr,n,311)
      color=4d0      !born matrix element colorfactor =CF*Trace(1_3)
     &     *(-1d0)/2d0 !Dipole colorfactor = -1/2 CA/CA 
      symmetry=1d0/96d0  ! 1/3*1/8 (inital state colors) * 1/2 *1/2 (initial state helicities)

      if (ldebug) then
      CAll UGDAAWj_msq(bornmom,borntest1,n,311,ps_number) !u~ g > d a ve
      borntest2=-((contract_Trjcj(metric,borncurr)
     &           )*4./96.)
      if (abs(borntest1-borntest2).ge.1e-10) then
      Print*, "error in matrixelement calculation 67:", (BORNTEST1-BORNTEST2)/BORNTEST1*100d0, "%"
      STOP
      endif
      endif

      elseif (sel.eq.2) then
      Call dguAAWj_curr(bornmom,borncurr,n,311)
      color=4d0      !born matrix element colorfactor =CF*Trace(1_3)
     &     *(-1d0)/2d0 !Dipole colorfactor = -1/2 CA/CA 
      symmetry=1d0/96d0  ! 1/3*1/8 (inital state colors) * 1/2 *1/2 (initial state helicities)

      if (ldebug) then      
      CAll DGUAAWj_msq(bornmom,borntest1,n,311,ps_number)
      borntest2=-((contract_Trjcj(metric,borncurr)
     &           )*4./96.)
      if (abs(borntest1-borntest2).ge.1e-10) then
      Print*, "error in matrixelement calculation 67:", (BORNTEST1-BORNTEST2)/BORNTEST1*100d0, "%"
      STOP
      endif
      endif

      endif

      do mu=0,3
      do nu=0,mu
      dipole(mu,nu)=-1d0/2d0/pipa/xika*8d0*pis*alfas*CF*
     &  ( (1d0-xika)/xika*2d0*ui*(1d0-ui)/pipk*(pmi(mu)/ui-pk(mu)/(1d0-ui))*(pmi(nu)/ui-pk(nu)/(1d0-ui)))
      if (mu.eq.nu.and.mu.eq.0) then
      dipole(mu,nu)=dipole(mu,nu)+1d0/2d0/pipa/xika*8d0*pis*alfas*CF*xika
      elseif (mu.eq.nu.and.mu.ne.0) then
      dipole(mu,nu)=dipole(mu,nu)-1d0/2d0/pipa/xika*8d0*pis*alfas*CF*xika
      endif
      if (mu.ne.nu) then
      dipole(nu,mu)=dipole(mu,nu)
      endif
      enddo
      enddo

      dipoleres(1)=(Contract_Trjcj(dipole(0,0),borncurr)
     &             )*color*symmetry

      do mu=0,3
      do nu=0,mu
      dipole(mu,nu)=-1d0/2d0/pipa/xika*16d0*pis*alfas*CA*
     &  ( (1d0-xika)/xika*ui*(1d0-ui)/pipk*(pmi(mu)/ui-pk(mu)/(1d0-ui))*(pmi(nu)/ui-pk(nu)/(1d0-ui)))
      if (mu.eq.nu.and.mu.eq.0) then
      dipole(mu,nu)=dipole(mu,nu)+1d0/2d0/pipa/xika*16d0*pis*alfas*CA*(1./(1.-xika+ui)-1.+xika*(1.-xika))
      elseif (mu.eq.nu.and.mu.ne.0) then
      dipole(mu,nu)=dipole(mu,nu)-1d0/2d0/pipa/xika*16d0*pis*alfas*CA*(1./(1.-xika+ui)-1.+xika*(1.-xika))
      endif
      if (mu.ne.nu) then
      dipole(nu,mu)=dipole(mu,nu)
      endif
      enddo
      enddo

      dipoleres(2)=(Contract_Trjcj(dipole(0,0),borncurr)
     &             )*color*symmetry


      dipole(0,0)=-1d0/2d0/pipa/xika*8.*pis*alfas*CF*(2./(1.-xika+ui)-(1.+xika))

      color=(CA-2.*CF)/2./CF*4.
      dipoleres(3)=-(Contract_Trjcj(metric(0,0),borncurr)
     &              )*color*symmetry*dipole(0,0)


      end


      subroutine daibdipole147_trib(n,sel,kin,xuz,pmi,pa,pb,p,v,dipoleres,
     &     ps_number)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"

      integer ps_number

      integer mu,nu,sel,kin,n
c INPUT
c
      double precision pmi(0:3),pb(0:3),pa(0:3),p(0:3,3),v(0:3,4)
      double precision pis,CF,CA
      parameter (pis=3.141592653589793238d0,CF=4d0/3d0,CA=3d0)
      logical ldebug
      parameter (ldebug=.false.)
c OUTPUT
c
      double precision dipole(0:3,0:3),xuz(2),pges(0:3),
     & pipb,papb,pipa,bornmom(0:3,7),dipoleres(3),color,symmetry,
     & borncache147(2,2),borncache145(3,2)

      common/BORNCACHE/borncache145,borncache147 

      double complex borncurr(0:3)
      double precision Contract_Trjcj,borntest1,borntest2,metric(0:3,0:3)

c
      double precision xiab,ui
      double precision dotrr
      external dotrr,Contract_Trjcj

      do mu=1,3
      dipoleres(mu) = 0
      enddo
      if(Loops_sub_NLO.ne.1) then
c     kin1 -> both, kin2 -> g
      if((sub_number.eq.2).and.(kin.eq.2)) return
      endif
      do mu=0,3
      do nu=0,mu-1
      metric(mu,nu)=0d0
      metric(nu,mu)=0d0
      enddo
      enddo
      metric(0,0)=1d0
      metric(1,1)=-1d0
      metric(2,2)=-1d0
      metric(3,3)=-1d0

      pipb=dotrr(pmi,pb)
      papb=dotrr(pa,pb)
      pipa=dotrr(pmi,pa)

      xiab= xuz(1)
      ui=xuz(2)

      if (ldebug) then
      if (abs(xiab-(papb-pipa-pipb)/papb).ge.1e-10) then
      PRINT*, "error related to xika in dipole 147"
      STOP
      endif
      do mu=0,3
      if (abs(p(mu,1)-xiab*pa(mu)).ge.1e-10) then
      print*, "wrong ptilde_a(",mu," ) in 147"
      STOP
      endif
      if (abs(p(mu,2)-pb(mu)).ge.1e-10) then
      print*, "wrong pb(",mu," ) in 147"
      STOP
      endif
      enddo
      endif

      do mu=0,3
      if (kin.eq.1) then
      bornmom(mu,1)=p(mu,2)!!!
      bornmom(mu,2)=p(mu,1)
      elseif (kin.eq.2) then
      bornmom(mu,1)=p(mu,1)!!!
      bornmom(mu,2)=p(mu,2)
      endif
      bornmom(mu,3)=p(mu,3)   !d~
      bornmom(mu,4)=v(mu,1)   !e-
      bornmom(mu,5)=v(mu,2)   !ve~
      bornmom(mu,6)=v(mu,3)   !a1
      bornmom(mu,7)=v(mu,4)   !a2
      enddo

      if (sel.eq.1) then
      Call ugdAAWj_curr(bornmom,borncurr,n,311)
      color=4d0      !born matrix element colorfactor =CF*Trace(1_3)
     &     *(-1d0)/2d0 !Dipole colorfactor = -1/2 CA/CA 
      symmetry=1d0/96d0  ! 1/2*1/8 (inital state colors) * 1/2 *1/2 (initial state helicities)
     
      if (ldebug) then
      CAll UGDAAWj_msq(bornmom,borntest1,n,311,ps_number)
      borntest2=-((contract_Trjcj(metric,borncurr)
     &           )*4./96.)
      if (abs(borntest1-borntest2).ge.1e-10) then
      Print*, "error in matrixelement calculation 147:", BORNTEST1, BORNTEST2
      STOP
      endif
      do mu=0,3
      pges(mu)=bornmom(mu,1)+bornmom(mu,2)
      enddo
!       Write(67,*) dsqrt(dotrr(pges,pges)),abs(borntest1-borntest2)/borntest1, borntest1, borntest2
      endif

      elseif (sel.eq.2) then
      Call dguAAWj_curr(bornmom,borncurr,n,311)
      color=4d0      !born matrix element colorfactor =CF*Trace(1_3)
     &     *(-1d0)/2d0 !Dipole colorfactor = -1/2 CA/CA 
      symmetry=1d0/96d0  ! 1/2*1/8 (inital state colors) * 1/2 *1/2 (initial state helicities)

    
      if (ldebug) then
      CAll DGUAAWj_msq(bornmom,borntest1,n,311,ps_number)
      borntest2=-((contract_Trjcj(metric,borncurr)
     &            )*4./96.)
      if (abs(borntest1-borntest2).ge.1e-10) then
      Print*, "error in matrixelement calculation 147:", BORNTEST1, BORNTEST2
      STOP
      endif
      endif

      endif

      do mu=0,3
      do nu=0,mu
      dipole(mu,nu)=-1d0/2d0/pipa/xiab*8d0*pis*alfas*CF*
     &  ( (1d0-xiab)/xiab*2d0*papb/(pipa*pipb)*(pmi(mu)-pipa/papb*pb(mu))*(pmi(nu)-pipa/papb*pb(nu)))
      if (mu.eq.nu.and.mu.eq.0) then
      dipole(mu,nu)=dipole(mu,nu)+1d0/2d0/pipa/xiab*8d0*pis*alfas*CF*xiab
      elseif (mu.eq.nu.and.mu.ne.0) then
      dipole(mu,nu)=dipole(mu,nu)-1d0/2d0/pipa/xiab*8d0*pis*alfas*CF*xiab
      endif
      if (mu.ne.nu) then
      dipole(nu,mu)=dipole(mu,nu)
      endif
      enddo
      enddo

      dipoleres(1)=(Contract_Trjcj(dipole(0,0),borncurr)
     &             )*color*symmetry


      do mu=0,3
      do nu=0,mu
      dipole(mu,nu)=-1d0/2d0/pipa/xiab*16d0*pis*alfas*CA*
     &  ( (1d0-xiab)/xiab*papb/pipa/pipb*(pmi(mu)-pipa/papb*pb(mu))*(pmi(nu)-pipa/papb*pb(nu)))
      if (mu.eq.nu.and.mu.eq.0) then
      dipole(mu,nu)=dipole(mu,nu)+1d0/2d0/pipa/xiab*16d0*pis*alfas*CA*(xiab/(1.-xiab)+xiab*(1.-xiab))
      elseif (mu.eq.nu.and.mu.ne.0) then
      dipole(mu,nu)=dipole(mu,nu)-1d0/2d0/pipa/xiab*16d0*pis*alfas*CA*(xiab/(1.-xiab)+xiab*(1.-xiab))
      endif
      if (mu.ne.nu) then
      dipole(nu,mu)=dipole(mu,nu)
      endif
      enddo
      enddo

      dipoleres(2)=(Contract_Trjcj(dipole(0,0),borncurr)
     &             )*color*symmetry


      dipole(0,0)=-1d0/2d0/pipa/xiab*8.*pis*alfas*CF*(2./(1.-xiab)-(1.+xiab))  
      color=-CA/2./CF*4.
      dipoleres(3)=-(Contract_Trjcj(metric(0,0),borncurr)
     &              )*color*symmetry*dipole(0,0)

      borncache147(sel,kin)=-((contract_Trjcj(metric,borncurr)
     &                       )*4./96.)

      end


      subroutine daisjdipole65_trib(n,sel,kin,xuz,pmi,pk,pa,p,v,dipoleres,
     &     ps_number)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"

      integer ps_number

      integer mu,sel,kin,n
c INPUT

      double precision pmi(0:3),pk(0:3),pa(0:3)
      double precision pis,CF,CA,Tr
      parameter (pis=3.141592653589793238d0,CF=4d0/3d0,CA=3d0,Tr=1d0/2d0)
      logical ldebug
      parameter (ldebug=.false.)

c OUTPUT
c
      double precision dipole,p(0:3,3),v(0:3,4)
     &    ,pipk,pkpa,pipa,xuz(2),dipoleres(6),bornmom(0:3,7),color(3)
c
      double precision born
      double precision Contract_Trjcj

      double precision xika,ui
      double precision dotrr
      external dotrr,Contract_Trjcj

      do mu=1,6
      dipoleres(mu) = 0
      enddo
      if((Loops_sub_NLO.ne.1).and.(sub_number.eq.2)) return
      pipk=dotrr(pmi,pk)
      pkpa=dotrr(pk,pa)
      pipa=dotrr(pmi,pa)

      xika= xuz(1)
      ui=xuz(2)

      if (ldebug) then
      if (abs(xika-(pkpa+pipa-pipk)/(pkpa+pipa)).ge.1e-10) then
      PRINT*, "error related to xika in dipole 65",xika,(pkpa+pipa-pipk)/(pkpa+pipa)
      STOP
      endif
      if (abs(ui-pipa/(pipa+pkpa)).ge.1e-10) then
      PRINT*, "error related to ui in dipole 65"
      STOP
      endif
      do mu=0,3
      if (abs(p(mu,1)-xika*pa(mu)).ge.1e-10) then
      print*, "wrong ptilde_a(",mu," ) in 65"
      STOP
      endif
      if (abs(p(mu,3)-(pk(mu)+pmi(mu)-(1d0-xika)*pa(mu))).ge.1e-10) then
      print*, "wrong ptilde_k(",mu," ) in 65"
      STOP
      endif
      enddo
      endif

      do mu=0,3
      if (kin.eq.1) then
      bornmom(mu,1)=p(mu,1)  ! u~
      bornmom(mu,2)=p(mu,2)  ! d
      elseif (kin.eq.2) then
      bornmom(mu,1)=p(mu,2)  ! u~
      bornmom(mu,2)=p(mu,1)  ! d
      endif
      bornmom(mu,3)=p(mu,3)  ! g
      bornmom(mu,4)=v(mu,1)  ! e-
      bornmom(mu,5)=v(mu,2)  ! ve~
      bornmom(mu,6)=v(mu,3)  ! a1
      bornmom(mu,7)=v(mu,4)  ! a2
      enddo

      color(1)=-CA/2./CF
      color(2)=(CA-2.*CF)/2./CF
      color(3)=0.

      if (sel.eq.1) then
      Call UDGAAWj_msq(bornmom,born,n,311, ps_number)
      elseif (sel.eq.2) then
      call UGDAAWj_msq(bornmom,born,n,311,ps_number)
      elseif (sel.eq.3) then
      Call DGUAAWj_msq(bornmom,born,n,311,ps_number)
      endif


      dipole=-1./2./pipa/xika*8.*pis*alfas*CF*(2./(1.-xika+ui)-(1.+xika))

      do mu=1,3
      dipoleres(mu)=born*dipole*color(mu)
      enddo

      dipole=-1./2./pipa/xika*8.*pis*alfas*TR*(1.-2.*xika*(1.-xika))

      do mu=4,6
      dipoleres(mu)=born*dipole*color(mu-3)
      enddo

      end


      subroutine daibdipole145_trib(n,sel,kin,xuz,pmi,pa,pb,p,v,dipoleres,
     &     ps_number)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"

      integer ps_number

      integer mu,sel,kin,n
c INPUT
c
      double precision pmi(0:3),pb(0:3),pa(0:3),p(0:3,3),v(0:3,4)
      double precision pis,CF,CA,TR
      parameter (pis=3.141592653589793238d0,CF=4d0/3d0,CA=3d0,TR=1d0/2d0)!
      logical ldebug
      parameter (ldebug=.false.)
c OUTPUT
c
      double precision dipole,xuz(2),
     &    pipb,papb,pipa,bornmom(0:3,7),dipoleres(6),color(3),
     &    borncache145(3,2),borncache147(2,2)

      common/BORNCACHE/borncache145,borncache147

      double precision born
      double precision Contract_Trjcj

c
      double precision xiab,ui
      double precision dotrr
      external dotrr,Contract_Trjcj

      do mu=1,6
      dipoleres(mu) = 0
      enddo
      if((Loops_sub_NLO.ne.1).and.(sub_number.eq.2)) return
      pipb=dotrr(pmi,pb)
      papb=dotrr(pa,pb)
      pipa=dotrr(pmi,pa)

      xiab= xuz(1)
      ui=xuz(2)

      if (ldebug) then
      if (abs(xiab-(papb-pipa-pipb)/papb).ge.1e-10) then
      PRINT*, "error related to xika in dipole 145"
      STOP
      endif
      do mu=0,3
      if (abs(p(mu,1)-xiab*pa(mu)).ge.1e-10) then
      print*, "wrong ptilde_a(",mu," ) in 145"
      STOP
      endif
      if (abs(p(mu,2)-pb(mu)).ge.1e-10) then
      print*, "wrong pb(",mu," ) in 145"
      STOP
      endif
      enddo
      endif

      do mu=0,3
      if (kin.eq.1) then
      bornmom(mu,1)=p(mu,1)  ! u~
      bornmom(mu,2)=p(mu,2)  ! d
      elseif (kin.eq.2) then
      bornmom(mu,1)=p(mu,2)  ! u~
      bornmom(mu,2)=p(mu,1)  ! d
      endif
      bornmom(mu,3)=p(mu,3)   !g
      bornmom(mu,4)=v(mu,1)   !e-
      bornmom(mu,5)=v(mu,2)   !ve~
      bornmom(mu,6)=v(mu,3)   !a1
      bornmom(mu,7)=v(mu,4)   !a2
      enddo

      color(1)=(CA-2d0*CF)/2d0/CF
      color(2)=-CA/2d0/CF
      color(3)=0d0

      if (sel.eq.1) then
      Call UDGAAWj_msq(bornmom,born,n,311, ps_number)
      elseif (sel.eq.2) then
      Call UGDAAWj_msq(bornmom,born,n,311,ps_number)
      elseif (sel.eq.3) then
      Call DGUAAWj_msq(bornmom,born,n,311,ps_number)
      endif

      dipole=-1d0/2d0/pipa/xiab*8d0*pis*alfas*CF*(2d0/(1d0-xiab)-(1d0+xiab))

      do mu=1,3
      dipoleres(mu)=born*dipole*color(mu)
      enddo

      dipole=-1d0/2d0/pipa/xiab*8d0*pis*alfas*TR*(1d0-2d0*xiab*(1d0-xiab))

      do mu=4,6
      dipoleres(mu)=born*dipole*color(mu-3)
      enddo

      borncache145(sel,kin)=born

      end


      subroutine dasijdipole40_trib(n,sel,kin,xuz,pmi,pj,pa,p,v,dipoleres,
     &     ps_number)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"

      integer ps_number

      integer mu,nu,sel,kin,n
c INPUT

c
      double precision pmi(0:3),pj(0:3),pa(0:3)
      double precision pis,CF,CA,TR
      parameter (pis=3.141592653589793238d0,CF=4d0/3d0,CA=3d0,TR=1d0/2d0)
      logical ldebug
      parameter (ldebug=.false.)

c OUTPUT
c
      double precision dipole(0:3,0:3),p(0:3,3),v(0:3,4)
     &    ,pipj,pjpa,pipa,xuz(2),dipoleres(3),bornmom(0:3,7),color,symmetry,borntmp,
     &    dipolesc
c
      double complex borncurr(0:3)
      double precision Contract_Trjcj

      double precision xija,zi,zj
      double precision dotrr
      external dotrr,Contract_Trjcj

      do mu=1,3
      dipoleres(mu) = 0
      enddo
      if(Loops_sub_NLO.ne.1) then
      if((sub_number.eq.1).and.(kin.ne.1)) return
      if((sub_number.eq.2).and.(kin.eq.1)) return
      endif
      borntmp=0d0
      dipolesc=0d0

      pipj=dotrr(pmi,pj)
      pjpa=dotrr(pj,pa)
      pipa=dotrr(pmi,pa)

      xija=xuz(1)
      zi=xuz(2)

      zj=1d0-zi

      if (ldebug) then
      if (abs(xija-(pipa+pjpa-pipj)/(pipa+pjpa)).ge.1e-10) then
      PRINT*, "error related to xija in dipole 40", xija, (pipa+pjpa-pipj)/(pipa+pjpa)
      STOP
      endif
      if (abs(zi-pipa/(pipa+pjpa)).ge.1e-10) then
      PRINT*, "error related to zi in dipole 40"
      STOP
      endif
      do mu=0,3
!       if (sel.eq.2) then
!       if (abs(p(mu,2)-xija*pa(mu)).ge.1e-10) then
!       print*, "kin 1 wrong ptilde_a(",mu," ) in 40", p(mu,2),xija*pa(mu)
!       STOP
!       endif
!       elseif (sel.eq.1) then
!       if (abs(p(mu,1)-xija*pa(mu)).ge.1e-10) then
!       print*, "kin 1 wrong ptilde_a(",mu," ) in 40", p(mu,2),xija*pa(mu)
!       STOP
!       endif
!       endif
      if (abs(p(mu,3)-(pj(mu)+pmi(mu)-(1d0-xija)*pa(mu))).ge.1e-10) then
      print*, "wrong ptilde_ij(",mu," ) in 40"
      STOP
      endif
      enddo

      endif

      do mu=0,3
      if (sel.eq.1) then
      bornmom(mu,1)=p(mu,1)
      bornmom(mu,2)=p(mu,2)
      elseif (sel.eq.2) then
      bornmom(mu,1)=p(mu,2)
      bornmom(mu,2)=p(mu,1)
      endif
      bornmom(mu,3)=p(mu,3)
      bornmom(mu,4)=v(mu,1)
      bornmom(mu,5)=v(mu,2)
      bornmom(mu,6)=v(mu,3)
      bornmom(mu,7)=v(mu,4)
      enddo

      if (kin.eq.1.or.kin.eq.4) then
      call udgAAWj_curr(bornmom,borncurr,n,311)
      color=4d0      !born matrix element colorfactor =CF*Trace(1_3)
     &     *(-1d0)/2d0 !Dipole colorfactor = -1/2 CA/CA 
      symmetry=1d0/36d0  ! 1/3*1/3 (inital state colors) * 1/2 *1/2 (initial state helicities)

      elseif (kin.eq.2) then
      CALL UGDAAWj_msq(bornmom,borntmp,n,311,ps_number)

      elseif (kin.eq.3) then
      CALL DGUAAWj_msq(bornmom,borntmp,n,311,ps_number)
      endif

!       if (ldebug) then
!       do mu=0,3
!       do nu=0,mu-1
!       metric(mu,nu)=0d0
!       metric(nu,mu)=0d0
!       enddo
!       enddo
!       metric(0,0)=1d0
!       metric(1,1)=-1d0
!       metric(2,2)=-1d0
!       metric(3,3)=-1d0
!       CAll SUDGENUAA(bornmom,borntest1)
!       borntest2=-((contract_Trjcj(metric,borncurr(0,-1))+contract_Trjcj(metric,borncurr(0,1)))*4./36.)
!       if (abs(borntest1-borntest2).ge.1e-7) then
!       Print*, "error in matrixelement calculation 40:", (BORNTEST1-BORNTEST2)/BORNTEST1*100d0, "%"
!       Print*, borntest1/borntest2
!       STOP
!       endif
!       endif


      if (kin.eq.1) then

      do mu=0,3
      do nu=0,mu
      dipole(mu,nu)=-1d0/2d0/pipj/xija*16d0*pis*alfas*CA*
     &  1d0/pipj*(zi*pmi(mu)-zj*pj(mu))*(zi*pmi(nu)-zj*pj(nu))
      if (mu.eq.nu.and.mu.eq.0) then
      dipole(mu,nu)=dipole(mu,nu)+1d0/2d0/pipj/xija*16d0*pis*alfas*CA*
     &  (1d0/(1d0-zi+(1d0-xija))+1d0/(1d0-zj+(1d0-xija))-2d0)
      elseif (mu.eq.nu.and.mu.ne.0) then
      dipole(mu,nu)=dipole(mu,nu)-1d0/2d0/pipj/xija*16d0*pis*alfas*CA*
     &  (1d0/(1d0-zi+(1d0-xija))+1d0/(1d0-zj+(1d0-xija))-2d0)
      endif
      if (mu.ne.nu) then
      dipole(nu,mu)=dipole(mu,nu)
      endif
      enddo
      enddo


      dipoleres(1)=(Contract_Trjcj(dipole(0,0),borncurr))
     &           *color*symmetry

      dipoleres(2)=0d0
      dipoleres(3)=0d0

      elseif (kin.eq.2.or.kin.eq.3) then

      dipolesc=-1d0/2d0/pipj/xija*8d0*pis*alfas*CF*(2d0/(1d0-zi+(1d0-xija))-(1d0+zi))

      color=-CA/2d0/CF
      dipoleres(1)=borntmp*color*dipolesc

      color=(CA-2d0*CF)/2d0/CF
      dipoleres(2)=borntmp*color*dipolesc
      color=-1d0/2d0
      dipoleres(3)=borntmp*color*dipolesc

      elseif (kin.eq.4) then
      color=-1./2.*4.

      do mu=0,3
      do nu=0,mu
      dipole(mu,nu)=-1d0/2d0/pipj/xija*8d0*pis*alfas*TR*
     &  (-2d0)/pipj*(zi*pmi(mu)-zj*pj(mu))*(zi*pmi(nu)-zj*pj(nu))
      if (mu.eq.nu.and.mu.eq.0) then
      dipole(mu,nu)=dipole(mu,nu)+1d0/2d0/pipj/xija*8d0*pis*alfas*TR
      elseif (mu.eq.nu.and.mu.ne.0) then
      dipole(mu,nu)=dipole(mu,nu)-1d0/2d0/pipj/xija*8d0*pis*alfas*TR
      endif
      if (mu.ne.nu) then
      dipole(nu,mu)=dipole(mu,nu)
      endif
      enddo
      enddo

      dipoleres(1)=(Contract_Trjcj(dipole(0,0),borncurr))
     &           *color*symmetry

      dipoleres(2)=0d0
      dipoleres(3)=0d0

      endif

      end


      subroutine test40_trib(n,sel,kin,xuz,pj,pmi,pa,p,v,dipoleres,ps_number)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"

      integer ps_number

      integer mu,nu,sel,kin,n
c INPUT

c
      double precision pmi(0:3),pj(0:3),pa(0:3)
      double precision pis,CF,CA,TR
      parameter (pis=3.141592653589793238d0,CF=4d0/3d0,CA=3d0,TR=1d0/2d0)
      logical ldebug
      parameter (ldebug=.false.)

c OUTPUT
c
      double precision p(0:3,3),v(0:3,4)
     &    ,pipj,pjpa,pipa,xuz(2),dipoleres(3),bornmom(0:3,7),color,borntmp,
     &    dipolesc
c
      double complex borncurr(0:3)
      double precision Contract_Trjcj,borntest2,metric(0:3,0:3),dptmp

      double precision xija,zi,ui,zj
      double precision dotrr
      external dotrr,Contract_Trjcj

      do mu=1,3
      dipoleres(mu) = 0
      enddo
      if((Loops_sub_NLO.ne.1).and.(sub_number.eq.2)) return
      borntmp=0d0
      dipolesc=0d0

      pipj=dotrr(pmi,pj)
      pjpa=dotrr(pj,pa)
      pipa=dotrr(pmi,pa)

      xija=xuz(1)
      ui=xuz(2)
      zi=pipa/(pipa+pjpa)
      zj=1d0-zi

      do mu=0,3
      if (sel.eq.1) then
      bornmom(mu,1)=p(mu,1)
      bornmom(mu,2)=p(mu,2)
      elseif (sel.eq.2) then
      bornmom(mu,1)=p(mu,2)
      bornmom(mu,2)=p(mu,1)
      endif
      bornmom(mu,3)=p(mu,3)
      bornmom(mu,4)=v(mu,1)
      bornmom(mu,5)=v(mu,2)
      bornmom(mu,6)=v(mu,3)
      bornmom(mu,7)=v(mu,4)
      enddo

      dptmp=-1./2./xija*8.*pis*alfas*4./3.*(2./(1.-zi+1.-xija)-1.-zi)!/pipj


      if (kin.eq.2) then
      CALL UGDAAWj_msq(bornmom,borntmp,n,311,ps_number)
      elseif (kin.eq.3) then
      CALL DGUAAWj_msq(bornmom,borntmp,n,311,ps_number)
      endif


      do mu=0,3
      do nu=0,mu-1
      metric(mu,nu)=0d0
      metric(nu,mu)=0d0
      enddo
      enddo
      metric(0,0)=1d0
      metric(1,1)=-1d0
      metric(2,2)=-1d0
      metric(3,3)=-1d0
      CAll ugdAAWj_curr(bornmom,borncurr,n,311)
      borntest2=-((contract_Trjcj(metric,borncurr)
     &           )*4./96.)


      dipolesc=-1./2./pipj/xija*8.*pis*alfas*CF*(2d0/(1d0-zi+(1d0-xija))-(1d0+zi))

      color=-CA/2d0/CF
      dipoleres(1)=borntmp*color*dipolesc

      color=(CA-2d0*CF)/2d0/CF
      dipoleres(2)=borntmp*color*dipolesc
      color=-1d0/2d0
      dipoleres(3)=borntmp*color*dipolesc

!       one=1.
!       two=2.
!       omx=one-xija
!       z=zi
!       omz=one-z
!       x=xija
!       print*, "pipa=",pipa," pjpa=",pjpa, " pipj=", pipj
! 
!       test1=4.*pis*alfas/pipj/x*(two/(omz+omx)-one-z)
!       test2=1./2./pipj/xija*8.*pis*alfas*(2d0/(1d0-zi+(1d0-xija))-(1d0+zi))
! 
!       print*, "My Dipole cons", test1, test2
!       print*, "Born of My Dipole", borntmp
!       print*, "My Dipole", borntmp*dipolesc

      end


!       subroutine dijkdipole(id,pmi,pj,pk,ptildeij,ptildek,dipole)
! 
!       implicit none
! 
!       include 'VBFNLO/utilities/scales.inc'
!       include 'VBFNLO/utilities/coupl.inc'
! 
!       integer mu,id
! c INPUT
! c
!       double precision pmi(0:3),pj(0:3),pk(0:3)
!       double precision pis,CF,CA,TR
!       parameter (pis=3.141592653589793d0,CF=4d0/3d0,CA=3d0,TR=1d0/2d0)
! c OUTPUT
! c
!       double precision ptildeij(0:3),ptildek(0:3),dipole
! c
!       double precision yijk,ztildei,ztildej
!       double precision dotrr
!       external dotrr
! 
!       yijk= dotrr(pmi,pj)/(dotrr(pmi,pj)+dotrr(pj,pk)+dotrr(pk,pmi))
! 
!       do mu=0,3
!       ptildek(mu) =1d0/(1d0-yijk)*pk(mu)
!       ptildeij(mu)=pmi(mu)+pj(mu)+pk(mu)-ptildek(mu)
!       enddo
!       ztildei=dotrr(pmi,ptildek)/dotrr(ptildeij,ptildek)
!       ztildej=1d0-ztildei
! 
!       if(id.eq.1) then
!       dipole=-1d0/dotrr(pmi,pj)*8d0*pis*alfas*CF*(2d0/(1d0-ztildei*(1d0-yijk))-(1d0+ztildei))  ! Eq 5.7 (qi,gj),k
!       elseif(id.eq.2) then
!       dipole=8d0*pis*alfas*TR*(2d0/(1d0-ztildei*(1d0-yijk))-(1d0+ztildei))  ! Eq 5.7 (qi,gj),k
!       endif
! 
! 
!       end
! 
! 
!       subroutine dasijdipole(id,pmi,pj,pa,ptildeij,ptildea,dipole)
! 
!       implicit none
! 
!       include 'VBFNLO/utilities/scales.inc'
!       include 'VBFNLO/utilities/coupl.inc'
! 
!       integer mu,id
! c INPUT
! c
!       double precision pmi(0:3),pj(0:3),pa(0:3)
!       double precision pis,CF,CA
!       parameter (pis=3.141592653589793d0,CF=4d0/3d0,CA=3d0)
! c OUTPUT
! c
!       double precision ptildeij(0:3),ptildea(0:3),dipole
! c
!       double precision xija,ztildei,ztildej
!       double precision dotrr
!       external dotrr
! 
!       xija= (dotrr(pmi,pa)+dotrr(pj,pa)-dotrr(pmi,pj))/(dotrr(pmi,pa)+dotrr(pj,pa))
! 
!       do mu=0,3
!       ptildea(mu) =xija*pa(mu)
!       ptildeij(mu)=pmi(mu)+pj(mu)-(1d0-xija)*ptildea(mu)
!       enddo
!       ztildei=dotrr(pmi,ptildea)/dotrr(ptildeij,ptildea)
!       ztildej=1d0-ztildei
! 
!       dipole=-1d0/2d0/dotrr(pmi,pj)/xija*8d0*pis*alfas*CF*(2d0/(1d0-ztildei*(1d0-xija))-(1d0+ztildei))
! 
! 
!       end
!
! 
!       subroutine daisjdipole(id,pmi,pk,pa,ptildeai,ptildek,dipole)
! 
!       implicit none
! 
!       include 'VBFNLO/utilities/scales.inc'
!       include 'VBFNLO/utilities/coupl.inc'
! 
!       integer mu,id
! c INPUT
! c
!       double precision pmi(0:3),pk(0:3),pa(0:3)
!       double precision pis,CF,CA
!       parameter (pis=3.141592653589793d0,CF=4d0/3d0,CA=3d0)
! c OUTPUT
! c
!       double precision ptildeai(0:3),ptildek(0:3),dipole
! c
!       double precision xika,ui
!       double precision dotrr
!       external dotrr
! 
!       xika= (dotrr(pk,pa)+dotrr(pmi,pa)-dotrr(pmi,pk))/(dotrr(pk,pa)+dotrr(pmi,pa))
! 
!       do mu=0,3
!       ptildeai(mu)=xika*pa(mu)
!       ptildek(mu) =pk(mu)+pmi(mu)-(1d0-xika)*pa(mu)
!       enddo
!       ui=dotrr(pmi,pa)/(dotrr(pmi,pa)+dotrr(pk,pa))
! 
!       dipole=-1d0/2d0/dotrr(pa,pmi)/xika*8d0*pis*alfas*CF*(2d0/(1d0-xika+ui)-(1d0+xika))
! 
! 
!       end


!       if (m2s(0).ne.0d0) then
!       write(66,*), dotrr(p(0,1,1),p(0,4,1)), m2s(1),m2s(0)-m2s(1)
!       endif

!       if (abs(dotrr(p(0,2,1),p(0,4,1))).le.1e-2.or.abs(dotrr(p(0,1,1),p(0,4,1))).le.1e-2.or.
!     &    abs(dotrr(p(0,2,1),p(0,3,1))).le.1e-2.or.abs(dotrr(p(0,1,1),p(0,3,1))).le.1e-2.or.
!     &    abs(dotrr(p(0,3,1),p(0,4,1))).le.1e-2) then
! !      if (abs(m2s(1)/(m2s(10)+m2s(11))).ge.0.9.and.abs(m2s(1)/(m2s(3)+m2s(5))).le.1.1) then
!       print*, "p2p4=",dotrr(p(0,2,1),p(0,4,1))
!       print*, "p1p4=",dotrr(p(0,1,1),p(0,4,1))
!       print*, "p2p3=",dotrr(p(0,2,1),p(0,3,1))
!       print*, "p1p3=",dotrr(p(0,1,1),p(0,3,1))
!       print*, "p1p2=",dotrr(p(0,1,1),p(0,2,1))
!       print*, "p3p4=",dotrr(p(0,3,1),p(0,4,1))
!       print*, "p1v3=",dotrr(p(0,1,1),v(0,3,1))
!       print*, "p2v3=",dotrr(p(0,2,1),v(0,3,1))
!       print*, "p3v3=",dotrr(p(0,3,1),v(0,3,1))
!       print*, "p4v3=",dotrr(p(0,4,1),v(0,3,1))
!       print*, "drive integration via ", m2s(0), " Fract ", m2s(0)/m2s(1)
!       print*, "real emission ME      ", m2s(1) , lokt(1)
!       print*, "Dipole 1              ", m2s(2) , lokt(2)
!       print*, "Dipole 2              ", m2s(3) , lokt(3)
!       print*, "Dipole 3              ", m2s(4) , lokt(4)
!       print*, "Dipole 4              ", m2s(5) , lokt(5)
!       print*, "Dipole 5              ", m2s(6) , lokt(6)
!       print*, "Dipole 6              ", m2s(7) , lokt(7)
!       print*, "Dipole 7              ", m2s(8) , lokt(8)
!       print*, "Dipole 8              ", m2s(9) , lokt(9)
!       print*, "Dipole 9              ", m2s(10) , lokt(10)
!       print*, "Dipole 10             ", m2s(11) , lokt(11)
!       print*, "REM/Dipole 1+3        ", m2s(1)/(m2s(2)+m2s(4))
!       print*, "REM/Dipole 2+4        ", m2s(1)/(m2s(3)+m2s(5))
!       print*, "REM/Dipole 5+7        ", m2s(1)/(m2s(6)+m2s(8))
!       print*, "REM/Dipole 6+8        ", m2s(1)/(m2s(7)+m2s(9))
!       print*, "REM/Dipole 9+10       ", m2s(1)/(m2s(10)+m2s(11))
!       print*, "---------------------------------"
! !!     pause
!       endif
! ! 
! !       print*, "ending scales", mufsq(1,1),mufsq(2,1),als(1,1)#
