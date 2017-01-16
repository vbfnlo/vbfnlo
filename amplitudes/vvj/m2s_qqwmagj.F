! This subroutine returns the real emission contribution to the order alpha_s 
! hadronic matrix element squared for the process p p --> e- ve gamma j
! jet. It is also used as part of the NLO calculation.

      subroutine m2s_qqwgagluj(xuz,nlo,pdfchange,lokt,xi,p,v,m2s)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "tensorlep.inc"
 
      integer i,j,k,mu,nlo,signpdf1,signpdf2

      logical lokt(11),fincolcalc,pdfchange

      integer lokstear(11)
      common/whywhat/lokstear

      double precision xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin),m2s_qqwgaglujfey,
     &   m2s_qqwgaglujf(4,3),pch(0:3,max_p,max_kin),
     &   muu(3),mdd(3),muuch(3),mddch(3),m2s_qqwgaglujch(4,3),
     &   q_sf,pdf(-6:6,2),amp,m2s_qqwgaglujdd(4,3),dotrr,m2s_qqwgaglujddch(4,3),
     &   helastestmom(0:3,7),mgg(4),mggch(4),m2s(0:11),
     &   fincol(12),fincolampsq(12),coll3,coll4,
     &   fincolconv,coll7,coll8,pdftmp(-6:6,2)

      double precision dipole(10),
     &   xuz(2,2:max_kin),kinem(11),dipolev(6,6),dipolea(3),dipoleb(3),borncache147(2,2),
     &   borncache(9,4),borncache145(3,2)

      double precision zero,inv2,four,two

      parameter (zero=0d0,inv2=1d0/2d0,four=4d0,two=2d0)

      common/BORNCACHE/borncache145,borncache147

      double complex Contract_Trjj
      external dotrr,Contract_Trjj

      parameter(fincolcalc=.true.)

      real*8 rsep
      external rsep

      INTEGER init/0/
      SAVE init

      if (.not.pdfchange) then

      if ( init .eq. 0 ) then
         write(6,*) ' '
         write(6,*) 'W-Ajj amplitude square information:'
         write(6,*) '------------------------------------------------'
         write(6,*) ' '
         write(6,*) ' '
         write(6,*) '------------------------------------------------'

         call printnfl(.true.)
         lbquark = 0
      endif

      elseif (pdfchange) then

      if ( init .eq. 0 ) then
         write(6,*) ' '
         write(6,*) 'W+Ajj amplitude square information:'
         write(6,*) '-----------------------------------------------'
         write(6,*) ' '
         write(6,*) ' '
         write(6,*) '-----------------------------------------------'

         call printnfl(.true.)
         lbquark = 0
      endif

      endif

      init=1

! initialize
      do i=1,11
      kinem(i)=zero
      enddo
      do i=1,6
      do j=1,3
      dipolev(i,j)=zero
      dipolea(j)=zero
      dipoleb(j)=zero
      enddo
      enddo
      do i=1,10
      dipole(i)=zero
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
! nlo!!!!!!
      call getyourscalesready(xi,p(0,1,1),v(0,1,1),pdfchange,pdf,nlo,1)

      if (lokt(1)) then
!       print*, "here 1"


! effective W currents
      do mu=0,3
      photfourvec(mu)=v(mu,3,1)
      wemnumom(mu)=v(mu,1,1)+v(mu,2,1)
      wemnugamom(mu)=v(mu,1,1)+v(mu,2,1)+v(mu,3,1)
      enddo
      photfourvec(4)=zero
      wemnugamom(4)=dotrr(wemnugamom,wemnugamom)
      CALL OXXXXX(v(0,1,1   ),zero ,-1, 1, wem ) !  e-
      CALL IXXXXX(v(0,2,1   ),zero , 1,-1, wve ) !  ve~
      Call JIOXXX(wve,wem,gwf,wmass,wwidth,wmcurr)
      do mu=0,3
      CALL VCARTX(wemnumom,WMASS , wwidth ,mu,-1, wavew)  
      CALL IOVXXX(wve,wem,wavew,GWF ,currlepw(mu))
      enddo
      do mu=-1,1,2
      Call HELVEC(photfourvec(0),1,mu,wafpolc(0,mu))
      CALL VXXXXX(v(0,3,1   ),zero ,mu,+1, waf(1,mu)) 
      Call SWTOEVEA(mu,currlepwga(0,mu))
      enddo

      do i=1,4
      call calcgg(p,v,i,mgg(i))
      if (i.ne.4) then
      call calcgg(pch,v,i,mggch(i))
      endif
      call calcuu(p,v,i,muu)
      call calcdd(p,v,i,mdd)
      call calcuu(pch,v,i,muuch)
      call calcdd(pch,v,i,mddch)
        do k=1,3
        m2s_qqwgaglujdd(i,k)=mdd(k)
        m2s_qqwgaglujddch(i,k)=mddch(k)
        m2s_qqwgaglujf(i,k)=muu(k)
        m2s_qqwgaglujch(i,k)=muuch(k)
        enddo
      enddo
!
!  Born - Kinematics
!
!     if (lokt(1)) then
!                                              PDF CHECK AGAINST DIPOLES
      amp= zero
     &   +m2s_qqwgaglujf(1,1)   *(pdf(-2,1)*pdf(-2,2)+pdf(-4,1)*pdf(-4,2))   !u~u~ > u~d~   !pdf check   ! 1
     &   +m2s_qqwgaglujf(2,1)   *(pdf(-2,1)*pdf( 2,2)+pdf(-4,1)*pdf( 4,2))   !u~u  > u d~   !pdf check   ! 2-3
     &   +m2s_qqwgaglujch(2,1)  *(pdf( 2,1)*pdf(-2,2)+pdf( 4,1)*pdf(-4,2))                  !pdf check 
     &   +m2s_qqwgaglujf(3,1)   *(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2)) !u~d  > u u~   !pdf check	  ! 4-5
     &   +m2s_qqwgaglujch(3,1)  *(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                  !pdf check
     &   +m2s_qqwgaglujf(4,1)*inv2   *(pdf( 2,1)*pdf( 1,2)+pdf( 4,1)*pdf( 3,2)) !u d  > u u  !pdf check   ! 6-7
     &   +m2s_qqwgaglujch(4,1)*inv2  *(pdf( 1,1)*pdf( 2,2)+pdf( 3,1)*pdf( 4,2))              !pdf check
     &   +m2s_qqwgaglujdd(1,1)  *(pdf( 1,1)*pdf( 1,2)+pdf( 3,1)*pdf( 3,2))   !d d  > d u    !pdf check   ! 8
     &   +m2s_qqwgaglujdd(2,1)*inv2  *(pdf(-1,1)*pdf(-2,2)+pdf(-3,1)*pdf(-4,2)) !d~u~ > d~d~ !pdf check   ! 9-10
     &   +m2s_qqwgaglujddch(2,1)*inv2*(pdf(-2,1)*pdf(-1,2)+pdf(-4,1)*pdf(-3,2))              !pdf check 
     &   +m2s_qqwgaglujdd(3,1)  *(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))   !u~d  > d d~   !pdf check   ! 11-12
     &   +m2s_qqwgaglujddch(3,1)*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                  !pdf check
     &   +m2s_qqwgaglujdd(4,1)  *(pdf(-1,1)*pdf( 1,2)+pdf(-3,1)*pdf( 3,2))   !d~d > u d~    !pdf check   ! 13-14
     &   +m2s_qqwgaglujddch(4,1)*(pdf( 1,1)*pdf(-1,2)+pdf( 3,1)*pdf(-3,2))                  !pdf check 
!!!!!!!!!!!!
! ! ! u~u > s~c, c~c > d~u  ! no dipole!!!!!!!!!!!!
     &   +m2s_qqwgaglujf(2,2)   *(pdf(-2,1)*pdf( 2,2)+pdf(-4,1)*pdf( 4,2))                  !pdf check   ! 15-16
     &   +m2s_qqwgaglujch(2,2)  *(pdf( 2,1)*pdf(-2,2)+pdf( 4,1)*pdf(-4,2))                  !pdf check
! ! ! u~d > ss~, c~s > dd~ 
     &   +m2s_qqwgaglujdd(3,2)  *(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))                  !pdf check   ! 17-18
     &   +m2s_qqwgaglujddch(3,2)*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                  !pdf check
! ! ! u~d > cc~, c~s > u u~ 
     &   +m2s_qqwgaglujf(3,2)  *(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))                   !pdf check   ! 19-20
     &   +m2s_qqwgaglujch(3,2) *(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))                   !pdf check
! ! ! u~s > u~c, c~d > c~u
     &   +m2s_qqwgaglujf(3,3) *(pdf(-2,1)*pdf( 3,2)+pdf(-4,1)*pdf( 1,2))                    !pdf check   ! 21-22
     &   +m2s_qqwgaglujch(3,3)*(pdf( 3,1)*pdf(-2,2)+pdf( 1,1)*pdf(-4,2))                    !pdf check 
! ! ! u~s~> d~s~, c~d~ > s~d~
     &   +m2s_qqwgaglujdd(2,3)  *(pdf(-2,2)*pdf(-3,1)+pdf(-4,2)*pdf(-1,1))                  !pdf check   ! 23-24
     &   +m2s_qqwgaglujddch(2,3)*(pdf(-3,2)*pdf(-2,1)+pdf(-1,2)*pdf(-4,1))                  !pdf check
! ! ! u~s > d~s, c~d > s~d 
     &   +m2s_qqwgaglujdd(3,3)  *(pdf(-2,1)*pdf( 3,2)+pdf(-4,1)*pdf( 1,2))                  !pdf check   ! 25-26
     &   +m2s_qqwgaglujddch(3,3)*(pdf( 3,1)*pdf(-2,2)+pdf( 1,1)*pdf(-4,2))                  !pdf check
! ! ! u~c~> u~s~, c~u~>c~d~
     &   +m2s_qqwgaglujf(1,2)  *(pdf(-2,1)*pdf(-4,2)+pdf(-4,1)*pdf(-2,2))                   !pdf check   ! 27-28
     &   +m2s_qqwgaglujch(1,2)*(pdf(-4,1)*pdf(-2,2)+pdf(-2,1)*pdf(-4,2))                    !pdf check
! ! ! ! u~c> d~c, c~u > s~u
     &   +m2s_qqwgaglujf(2,3)  *(pdf(-2,1)*pdf( 4,2)+pdf(-4,1)*pdf( 2,2))                   !pdf check   ! 29-30 x
     &   +m2s_qqwgaglujch(2,3)*(pdf( 4,1)*pdf(-2,2)+pdf( 2,1)*pdf(-4,2))                    !pdf check
! ! ! u s> uc, cd > c u
     &   +m2s_qqwgaglujf(4,2)  *(pdf( 2,1)*pdf( 3,2)+pdf( 4,1)*pdf( 1,2))                   !pdf check   ! 31-32 
     &   +m2s_qqwgaglujch(4,2)*(pdf( 3,1)*pdf( 2,2)+pdf( 1,1)*pdf( 4,2))                    !pdf check
! ! ! d~d> s~c, s~s > d~u ! no dipole
     &   +m2s_qqwgaglujdd(4,2)  *(pdf(-1,1)*pdf( 1,2)+pdf(-3,1)*pdf( 3,2))                  !pdf check   ! 33-34
     &   +m2s_qqwgaglujddch(4,2)*(pdf( 1,1)*pdf(-1,2)+pdf( 3,1)*pdf(-3,2))                  !pdf check 
! ! ! d~s> d~c, s~d > s~u
     &   +m2s_qqwgaglujdd(4,3)  *(pdf(-1,1)*pdf( 3,2)+pdf(-3,1)*pdf( 1,2))                  !pdf check   ! 35-36 x
     &   +m2s_qqwgaglujddch(4,3)*(pdf( 3,1)*pdf(-1,2)+pdf( 1,1)*pdf(-3,2))                  !pdf check
! ! ! d s> u s, s d > c d   
     &   +m2s_qqwgaglujdd(1,3)  *(pdf( 1,1)*pdf( 3,2)+pdf( 3,1)*pdf( 1,2))                  !pdf check   ! 37-38
     &   +m2s_qqwgaglujddch(1,3)*(pdf( 3,1)*pdf( 1,2)+pdf( 1,1)*pdf( 3,2))                  !pdf check
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

      call getyourscalesready(xi,p(0,1,2),v(0,1,2),pdfchange,pdftmp,nlo,2)

      Call daisjdipole67(1,1,xuz(1,2),p(0,4,1),p(0,3,1),p(0,1,1),p(0,1,2),v(0,1,2),dipolea(1))
      Call daisjdipole67(2,1,xuz(1,2),p(0,4,1),p(0,3,1),p(0,1,1),p(0,1,2),v(0,1,2),dipoleb(1))
      Call daisjdipole65(1,1,xuz(1,2),p(0,4,1),p(0,3,1),p(0,1,1),p(0,1,2),v(0,1,2),dipolev(1,1))
      Call daisjdipole65(1,2,xuz(1,2),p(0,4,1),p(0,3,1),p(0,1,1),p(0,1,2),v(0,1,2),dipolev(1,2))

      Call daisjdipole65(2,1,xuz(1,2),p(0,4,1),p(0,3,1),p(0,1,1),p(0,1,2),v(0,1,2),dipolev(1,3))

      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)

      kinem(1)= zero 
     &        + dipole(1)*(pdf(-2,1)*pdf(-2,2)+pdf(-4,1)*pdf(-4,2))    !u~u~ > u~d~  !pdf check
     &        + dipole(1)*(pdf( 2,1)*pdf(-2,2)+pdf( 4,1)*pdf(-4,2))    !u~u  > u d~  !pdf check
     &        + dipole(2)*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))    !u~d  > u u~  !pdf check
     &        + dipole(2)*inv2*(pdf( 2,1)*pdf( 1,2)+pdf( 4,1)*pdf( 3,2))!u d  > u u   !pdf check  ! 6-7
     &        + dipole(2)*(pdf( 1,1)*pdf( 1,2)+pdf( 3,1)*pdf( 3,2))    !d d  > d u   !pdf check  ! 8
     &        + dipole(1)*inv2*(pdf(-1,1)*pdf(-2,2)+pdf(-3,1)*pdf(-4,2))!d~u~ > d~d~  !pdf check
     &        + dipole(1)*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))    !u~d  > d d~  !pdf check
     &        + dipole(2)*(pdf(-1,1)*pdf( 1,2)+pdf(-3,1)*pdf( 3,2))    !d~d  > u d~  !pdf check
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

      call getyourscalesready(xi,p(0,1,3),v(0,1,3),pdfchange,pdftmp,nlo,3)

      Call daisjdipole67(1,1,xuz(1,3),p(0,4,1),p(0,3,1),p(0,2,1),p(0,1,3),v(0,1,3),dipolea(1))
      Call daisjdipole67(2,1,xuz(1,3),p(0,4,1),p(0,3,1),p(0,2,1),p(0,1,3),v(0,1,3),dipoleb(1))
      Call daisjdipole65(1,2,xuz(1,3),p(0,4,1),p(0,3,1),p(0,2,1),p(0,1,3),v(0,1,3),dipolev(1,1))
      Call daisjdipole65(1,1,xuz(1,3),p(0,4,1),p(0,3,1),p(0,2,1),p(0,1,3),v(0,1,3),dipolev(1,2))
      Call daisjdipole65(2,1,xuz(1,3),p(0,4,1),p(0,3,1),p(0,2,1),p(0,1,3),v(0,1,3),dipolev(1,3))


      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)

      kinem(2)= zero
     &        + dipole(1)*(pdf(-2,1)*pdf(-2,2)+pdf(-4,1)*pdf(-4,2))   !u~u~ > u~d~  !pdf check
     &        + dipole(1)*(pdf(-2,1)*pdf( 2,2)+pdf(-4,1)*pdf( 4,2))   !u~u  > u d~  !pdf check
     &        + dipole(2)*(pdf( 1,1)*pdf(-2,2)+pdf( 3,1)*pdf(-4,2))   !u~d  > u u~  !pdf check
     &        + dipole(2)*inv2*(pdf( 1,1)*pdf( 2,2)+pdf( 3,1)*pdf( 4,2)) !u d   > u u  !pdf check !6-7
     &        + dipole(2)*(pdf( 1,1)*pdf( 1,2)+pdf( 3,1)*pdf( 3,2))   !d d  > d u   !pdf check   ! 8
     &        + dipole(1)*inv2*(pdf(-2,1)*pdf(-1,2)+pdf(-4,1)*pdf(-3,2))   !d~u~ > d~d~  !pdf check
     &        + dipole(1)*(pdf(-2,1)*pdf( 1,2)+pdf(-4,1)*pdf( 3,2))   !u~d  > d d~  !pdf check
     &        + dipole(2)*(pdf( 1,1)*pdf(-1,2)+pdf( 3,1)*pdf(-3,2))   !d~d  > u d~  !pdf check
!!!!!!!
     &        + dipole(2)*(pdf( 3,1)*pdf(-2,2)+pdf( 1,1)*pdf(-4,2))   !u~s  > u~c , c~d  > c~u   !pdf check
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

      call getyourscalesready(xi,p(0,1,4),v(0,1,4),pdfchange,pdftmp,nlo,4)
 
      Call daibdipole147(1,1,xuz(1,4),p(0,4,1),p(0,1,1),p(0,2,1),p(0,1,4),v(0,1,4),dipolea(1))
      Call daibdipole147(2,1,xuz(1,4),p(0,4,1),p(0,1,1),p(0,2,1),p(0,1,4),v(0,1,4),dipoleb(1))
      Call daibdipole145(1,1,xuz(1,4),p(0,4,1),p(0,1,1),p(0,2,1),p(0,1,4),v(0,1,4),dipolev(1,1))
      Call daibdipole145(1,2,xuz(1,4),p(0,4,1),p(0,1,1),p(0,2,1),p(0,1,4),v(0,1,4),dipolev(1,2))

      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)

      Call daibdipole145(2,1,xuz(1,4),p(0,4,1),p(0,1,1),p(0,2,1),p(0,1,4),v(0,1,4),dipolev(1,3))

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

      call getyourscalesready(xi,p(0,1,5),v(0,1,5),pdfchange,pdftmp,nlo,5)

      Call daibdipole147(1,1,xuz(1,5),p(0,4,1),p(0,2,1),p(0,1,1),p(0,1,5),v(0,1,5),dipolea(1))
      Call daibdipole147(2,1,xuz(1,5),p(0,4,1),p(0,2,1),p(0,1,1),p(0,1,5),v(0,1,5),dipoleb(1))
      Call daibdipole145(1,2,xuz(1,5),p(0,4,1),p(0,2,1),p(0,1,1),p(0,1,5),v(0,1,5),dipolev(1,1))
      Call daibdipole145(1,1,xuz(1,5),p(0,4,1),p(0,2,1),p(0,1,1),p(0,1,5),v(0,1,5),dipolev(1,2))

      Call daibdipole145(2,1,xuz(1,5),p(0,4,1),p(0,2,1),p(0,1,1),p(0,1,5),v(0,1,5),dipolev(1,3))


      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)
      dipole(3)=dipolea(2)

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

      call getyourscalesready(xi,p(0,1,6),v(0,1,6),pdfchange,pdftmp,nlo,6)

      Call daisjdipole67(1,1,xuz(1,6),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipolea(1))
      Call daisjdipole67(2,1,xuz(1,6),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipoleb(1))
      Call daisjdipole65(1,1,xuz(1,6),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipolev(1,1))
      Call daisjdipole65(1,2,xuz(1,6),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipolev(1,2))

      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)
      dipole(4)=dipolea(2)
      dipole(6)=dipoleb(2)

      Call daisjdipole67(1,2,xuz(1,6),p(0,3,1), p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipolea(1))
      Call daisjdipole67(2,2,xuz(1,6),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipoleb(1))
      Call daisjdipole65(3,1,xuz(1,6),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipolev(1,3))

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

      call getyourscalesready(xi,p(0,1,7),v(0,1,7),pdfchange,pdftmp,nlo,7)

      Call daisjdipole67(1,1,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipolea(1))
      Call daisjdipole67(2,1,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipoleb(1))
      Call daisjdipole65(1,2,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipolev(1,1))
      Call daisjdipole65(1,1,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipolev(1,2))

      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)
      dipole(3)=dipolea(2)
      dipole(5)=dipoleb(2)

      Call daisjdipole67(1,2,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipolea(1))
      Call daisjdipole67(2,2,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipoleb(1))

      Call daisjdipole65(3,1,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipolev(1,3))


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

      call getyourscalesready(xi,p(0,1,8),v(0,1,8),pdfchange,pdftmp,nlo,8)

      Call daibdipole147(1,1,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipolea(1))
      Call daibdipole147(2,1,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipoleb(1))
      Call daibdipole145(1,1,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipolev(1,1))
      Call daibdipole145(1,2,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipolev(1,2))

      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)
      dipole(4)=dipolea(2)
      dipole(6)=dipoleb(2)

      Call daibdipole147(1,2,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipolea(1))
      Call daibdipole147(2,2,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipoleb(1))
      Call daibdipole145(3,1,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipolev(1,3))

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

      call getyourscalesready(xi,p(0,1,9),v(0,1,9),pdfchange,pdftmp,nlo,9)

      Call daibdipole147(1,1,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipolea(1))
      Call daibdipole147(2,1,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipoleb(1))
      Call daibdipole145(1,2,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipolev(1,1))
      Call daibdipole145(1,1,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipolev(1,2))

      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)
      dipole(3)=dipolea(2)
      dipole(5)=dipoleb(2)

      Call daibdipole147(1,2,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipolea(1))
      Call daibdipole147(2,2,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipoleb(1))
      Call daibdipole145(3,1,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipolev(1,3))

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

      call getyourscalesready(xi,p(0,1,10),v(0,1,10),pdfchange,pdftmp,nlo,10)

      CALL dasijdipole40(1,1,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1))
      CALL dasijdipole40(2,1,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1))

      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)

!       CALL dasijdipole40(1,2,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1))
!       CALL dasijdipole40(2,2,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1))
! 
!       dipole(3)=dipolea(2)
!       dipole(4)=dipoleb(1)

      CALL test40(1,2,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1))
      CALL test40(2,2,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1))

      dipole(3)=dipolea(2)
      dipole(4)=dipoleb(1)

!       CALL dasijdipole40(1,3,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1))
!       CALL dasijdipole40(2,3,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1))

      CALL test40(1,3,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1))
      CALL test40(2,3,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1))

      dipole(5)=dipolea(2)
      dipole(6)=dipoleb(1)

      CALL dasijdipole40(1,4,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1))
      CALL dasijdipole40(2,4,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1))

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

      call getyourscalesready(xi,p(0,1,11),v(0,1,11),pdfchange,pdftmp,nlo,11)

      CALL dasijdipole40(1,1,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipolea(1))
      CALL dasijdipole40(2,1,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipoleb(1))

      dipole(1)=dipolea(1)
      dipole(2)=dipoleb(1)

!       CALL dasijdipole40(1,2,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipolea(1))
!       CALL dasijdipole40(2,2,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipoleb(1))

      CALL test40(1,2,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipolea(1))
      CALL test40(2,2,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipoleb(1))
      dipole(3)=dipolea(1)
      dipole(4)=dipoleb(2)

!       CALL dasijdipole40(1,3,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipolea(1))
!       CALL dasijdipole40(2,3,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipoleb(1))

      CALL test40(1,3,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipolea(1))
      CALL test40(2,3,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipoleb(1))
! 
      dipole(5)=dipolea(1)
      dipole(6)=dipoleb(2)

      CALL dasijdipole40(1,4,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipolea(1))
      CALL dasijdipole40(2,4,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipoleb(1))

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
      if (lokt(4).and.fincolcalc) then !this is one of the initial-initial kinematics, e.g. eq 5.147 in Catani-Seymour
!                                                                                       first splits
      call getyourscalesready(xi,p(0,1,4),v(0,1,4),pdfchange,pdf,nlo,4)
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
!       CAll SUDGAENU(helastestmom,fincolampsq(1))
      fincolampsq(1)=borncache(5,1) 
!       Call SUGDAEVE(helastestmom,fincolampsq(3))
      fincolampsq(3)=borncache(7,1) 

!       do k=0,3
!          helastestmom(k,1)=p(k,2,4)
!          helastestmom(k,2)=p(k,1,4)
!          helastestmom(k,3)=p(k,3,4)
!          helastestmom(k,4)=v(k,3,4)
!          helastestmom(k,5)=v(k,1,4)
!          helastestmom(k,6)=v(k,2,4)
!       enddo

!      CAll SUDGAENU(helastestmom,fincolampsq(2))
      fincolampsq(2)=borncache(6,1)
!      Call SUGDAEVE(helastestmom,fincolampsq(4))
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

      call getyourscalesready(xi,p(0,1,8),v(0,1,8),pdfchange,pdf,nlo,8)
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
!       Call SDGUAEVE(helastestmom,fincolampsq(5))
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
!       Call SDGUAEVE(helastestmom,fincolampsq(6))
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
      call getyourscalesready(xi,p(0,1,5),v(0,1,5),pdfchange,pdf,nlo,5)
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

!       CAll SUDGAENU(helastestmom,fincolampsq(1))
      fincolampsq(1)=borncache(5,3)
!       Call SUGDAEVE(helastestmom,fincolampsq(3))
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
!       CAll SUDGAENU(helastestmom,fincolampsq(2))
      fincolampsq(2)=borncache(6,3)
!       Call SUGDAEVE(helastestmom,fincolampsq(4))
      fincolampsq(4)=borncache(1,3)
!       Call SDGUAEVE(helastestmom,fincolampsq(6))
      

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

      call getyourscalesready(xi,p(0,1,9),v(0,1,9),pdfchange,pdf,nlo,9)
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
!       Call SDGUAEVE(helastestmom,fincolampsq(5))
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
!       Call SDGUAEVE(helastestmom,fincolampsq(6))
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
      if (abs(dotrr(p(0,2,1),p(0,4,1))).le.1e-2.or.abs(dotrr(p(0,1,1),p(0,4,1))).le.1e-2.or.
     &    abs(dotrr(p(0,2,1),p(0,3,1))).le.1e-2.or.abs(dotrr(p(0,1,1),p(0,3,1))).le.1e-2.or.
     &    abs(dotrr(p(0,3,1),p(0,4,1))).le.1e-2) then
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
      print*, "Dipole 1              ", m2s(2) , lokt(2)
      print*, "Dipole 2              ", m2s(3) , lokt(3)
      print*, "Dipole 3              ", m2s(4) , lokt(4)
      print*, "Dipole 4              ", m2s(5) , lokt(5)
      print*, "Dipole 5              ", m2s(6) , lokt(6)
      print*, "Dipole 6              ", m2s(7) , lokt(7)
      print*, "Dipole 7              ", m2s(8) , lokt(8)
      print*, "Dipole 8              ", m2s(9) , lokt(9)
      print*, "Dipole 9              ", m2s(10) , lokt(10)
      print*, "Dipole 10             ", m2s(11) , lokt(11)
      print*, "REM/Dipole 1+3        ", m2s(1)/(m2s(2)+m2s(4))
      print*, "REM/Dipole 2+4        ", m2s(1)/(m2s(3)+m2s(5))
      print*, "REM/Dipole 5+7        ", m2s(1)/(m2s(6)+m2s(8))
      print*, "REM/Dipole 6+8        ", m2s(1)/(m2s(7)+m2s(9))
      print*, "REM/Dipole 9+10       ", m2s(1)/(m2s(10)+m2s(11))
      print*, "---------------------------------"
!      pause
      endif
      endif

      end


      subroutine getyourscalesready(xi,p,v,pdfchange,pdf,nlo,l)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"

      double precision q_sf,pdf(-6:6,2),pdftmp(-6:6,2),xi(nx),four,
     &  p(0:3,max_p),v(0:3,max_v)
      integer k,nlo
      logical pdfchange
      integer l

      parameter(four=4d0)


      q_sf = dsqrt(mufsq(1,l))

      call pdfproton( xi(1), q_sf, pdf(-6,1) )        ! f_a(x1)

      q_sf = dsqrt(mufsq(2,l))
      call pdfproton( xi(2), q_sf, pdf(-6,2) )        ! f_b(x2)

      do k=-6,6,1
        pdftmp(k,1)=pdf(k,1)/xi(1)
        pdftmp(k,2)=pdf(k,2)/xi(2)
      enddo
      do k=-6,6,1
        pdf(k,1)=pdftmp(sign1*k,1)
        pdf(k,2)=pdftmp(sign2*k,2)
      enddo

! change pdfs for W+ case
      if (pdfchange) then
      do k=-6,6,1
        pdftmp(k,1)=pdf(k,1)
        pdftmp(k,2)=pdf(k,2)
      enddo
      do k=-6,6,1
        pdf(k,1)=pdftmp(-k,1)
        pdf(k,2)=pdftmp(-k,2)
      enddo
      endif


! running coupling
      alfas  =  als(1,l) 
      G      =  dsqrt(four*pi*alfas)
      GG(1)  =  -G
      GG(2)  =  -G

!      print*, "scales routine ",xi, xi2

      end


      subroutine finitecollqqREM(xa,x,patilde,pb,pin,mufac,quark,result,pdfchange)

      implicit none

#include "VBFNLO/utilities/global.inc"

      integer quark
      integer pdfchange

      double precision xa,x,patilde(0:3),pb(0:3),pin(0:3),result,dotrr,
     &   xxa,pdfqxa,pdfgxa,pdfqxxa,pdfgxxa,logomx, logomxxa,logx,omx,omxxa,
     &   Logmupapb2,logmupapi2,NF,mufacsq,mufac,pdftmp(-6:6),one,inv2

      parameter (nf=5d0,one=1d0,inv2=1d0/2d0)

      external dotrr

      mufacsq=mufac**2
      xxa=x*xa
      omxxa=one-xxa
      omx=one-x
!       print*, xa,xxa, " in rems"

      logx=Log(x)
!       if (x.ge.0.999995) logx=0d0

      logomx=Log(omx)
!       if (x.le.0.000005) logomx=0d0

      logomxxa=Log(omxxa)
!       if (xxa.le.0.000005) logomxxa=0d0

      
      logmupapb2=Log(mufacsq*inv2/dotrr(patilde,pb))
      logmupapi2=Log(mufacsq*inv2/dotrr(patilde,pin))

      call pdfproton( xa, mufac, pdftmp(-6) )
      pdfqxa=pdftmp(pdfchange*quark)/xa
      pdfgxa=pdftmp(0)/xa
      call pdfproton( xxa, mufac, pdftmp(-6) )   
      pdfqxxa=pdftmp(pdfchange*quark)/xxa
      pdfgxxa=pdftmp(0)/xxa
  
      result=(logmupapb2*pdfgxa)/(16.*omx) - (9.*logmupapi2*pdfgxa)/(16.*omx) + (7.*Logomx*pdfgxa)/(16.*omx) - 
     -  (Logx*pdfgxa)/(2.*omx) - (11.*pdfqxa)/(4.*omx**2) + (logmupapb2*pdfqxa)/(3.*omx**2) - 
     -  (3.*logmupapi2*pdfqxa)/omx**2 + (7.*Logomx*pdfqxa)/(3.*omx**2) - (4.*Logx*pdfqxa)/omx**2 + 
     -  (Nf*pdfqxa)/(6.*omx**2) + (4.*pdfqxa)/(3.*omx) - (logmupapb2*pdfqxa)/(6.*omx) + 
     -  (3.*logmupapi2*pdfqxa)/(2.*omx) - (7.*Logomx*pdfqxa)/(6.*omx) + (8.*Logx*pdfqxa)/(3.*omx) + (pdfgxa*x)/omx - 
     -  (logmupapb2*pdfgxa*x)/(8.*omx) + (9.*logmupapi2*pdfgxa*x)/(8.*omx) - (7.*Logomx*pdfgxa*x)/(8.*omx) + 
     -  (Logx*pdfgxa*x)/omx - (4.*pdfqxa*x)/(3.*omx) - (logmupapb2*pdfqxa*x)/(6.*omx) + 
     -  (3.*logmupapi2*pdfqxa*x)/(2.*omx) - (7.*Logomx*pdfqxa*x)/(6.*omx) + (8.*Logx*pdfqxa*x)/(3.*omx) + 
     -  (11.*pdfqxxa*x)/(4.*omx**2) - (logmupapb2*pdfqxxa*x)/(3.*omx**2) + (3.*logmupapi2*pdfqxxa*x)/omx**2 - 
     -  (7.*Logomx*pdfqxxa*x)/(3.*omx**2) - (Nf*pdfqxxa*x)/(6.*omx**2) - (113.*pdfqxxa*x)/(12.*omx*omxxa) + 
     -  (logmupapb2*pdfqxxa*x)/(4.*omx*omxxa) - (9.*logmupapi2*pdfqxxa*x)/(4.*omx*omxxa) - 
     -  (11.*Logomxxa*pdfqxxa*x)/(4.*omx*omxxa) + (logmupapb2*Logomxxa*pdfqxxa*x)/(3.*omx*omxxa) - 
     -  (3.*logmupapi2*Logomxxa*pdfqxxa*x)/(omx*omxxa) + (7.*Logomxxa**2*pdfqxxa*x)/(6.*omx*omxxa) + 
     -  (Nf*pdfqxxa*x)/(6.*omx*omxxa) + (Logomxxa*Nf*pdfqxxa*x)/(6.*omx*omxxa) + 
     -  (17.*pdfqxxa*Pi**2*x)/(18.*omx*omxxa) - (pdfgxa*x**2)/omx + (logmupapb2*pdfgxa*x**2)/(8.*omx) - 
     -  (9.*logmupapi2*pdfgxa*x**2)/(8.*omx) + (7.*Logomx*pdfgxa*x**2)/(8.*omx) - (Logx*pdfgxa*x**2)/omx + 
     -  (4.*Logx*pdfqxa*x**2)/(3.*omx**2)


      end


      subroutine finitecollgqREM(xa,x,patilde,pb,pin,mufac,result,pdfchange)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"

      integer quark
      integer pdfchange
      integer maxflav

      double precision xa,x,patilde(0:3),pb(0:3),pin(0:3),result,dotrr,
     &   xxa,pdfqxa,pdfgxa,pdfqxxa,pdfgxxa,logomx, logomxxa,logx,omx,omxxa,
     &   Logmupapb2,logmupapi2,NF,mufacsq,mufac,pdftmp(-6:6),one,inv2,zero
      parameter (nf=5d0,one=1d0,inv2=1d0/2d0,zero=0d0)

      external dotrr

      if (lbquark .eq. 0) then
        maxflav = 4
      else
        maxflav = 5
      endif

      mufacsq=mufac**2
      xxa=x*xa
      omxxa=one-xxa
      omx=one-x

      logx=Log(x)
!       if (x.ge.0.999995) logx=0d0

      logomx=Log(omx)
!       if (x.le.0.000005) logomx=0d0

      logomxxa=Log(omxxa)
!       if (xxa.le.0.000005) logomxxa=0d0

      logmupapb2=Log(mufacsq*inv2/dotrr(patilde,pb))
      logmupapi2=Log(mufacsq*inv2/dotrr(patilde,pin))

      pdfqxa=zero
      pdfqxxa=zero
      
      call pdfproton( xa, mufac, pdftmp(-6) )
      pdfgxa=pdftmp(0)/xa

      do quark=1,maxflav  ! sum over quark gluon splittings
      pdfqxa=pdfqxa+pdftmp(pdfchange*quark)/xa
      enddo
      do quark=-maxflav,-1
      pdfqxa=pdfqxa+pdftmp(pdfchange*quark)/xa
      enddo

      call pdfproton( xxa, mufac, pdftmp(-6) )    
      pdfgxxa=pdftmp(0)/xxa
      do quark=1,maxflav  ! sum over quark gluon splittings
      pdfqxxa=pdfqxxa+pdftmp(pdfchange*quark)/xxa
      enddo
      do quark=-maxflav,-1
      pdfqxxa=pdfqxxa+pdftmp(pdfchange*quark)/xxa
      enddo

      result=(-9.*pdfgxa)/(4.*omx**2) - (3.*logmupapb2*pdfgxa)/omx**2 - (3.*logmupapi2*pdfgxa)/omx**2 + 
     -  (9.*Logomx*pdfgxa)/omx**2 - (6.*Logx*pdfgxa)/omx**2 + (6.*logmupapb2*pdfgxa)/omx + (6.*logmupapi2*pdfgxa)/omx - 
     -  (18.*Logomx*pdfgxa)/omx + (12.*Logx*pdfgxa)/omx + (4.*logmupapb2*pdfqxa)/(3.*omx) + 
     -  (4.*logmupapi2*pdfqxa)/(3.*omx) - (4.*Logomx*pdfqxa)/omx + (8.*Logx*pdfqxa)/(3.*omx) - 
     -  (3.*logmupapb2*pdfgxa)/(omx*x) - (3.*logmupapi2*pdfgxa)/(omx*x) + (9.*Logomx*pdfgxa)/(omx*x) - 
     -  (6.*Logx*pdfgxa)/(omx*x) - (4.*logmupapb2*pdfqxa)/(3.*omx*x) - (4.*logmupapi2*pdfqxa)/(3.*omx*x) + 
     -  (4.*Logomx*pdfqxa)/(omx*x) - (8.*Logx*pdfqxa)/(3.*omx*x) - (3.*logmupapb2*pdfgxa*x)/omx - 
     -  (3.*logmupapi2*pdfgxa*x)/omx + (9.*Logomx*pdfgxa*x)/omx - (6.*Logx*pdfgxa*x)/omx + (9.*pdfgxxa*x)/(4.*omx**2) + 
     -  (3.*logmupapb2*pdfgxxa*x)/omx**2 + (3.*logmupapi2*pdfgxxa*x)/omx**2 - (9.*Logomx*pdfgxxa*x)/omx**2 - 
     -  (227.*pdfgxxa*x)/(12.*omx*omxxa) - (11.*logmupapb2*pdfgxxa*x)/(4.*omx*omxxa) - 
     -  (11.*logmupapi2*pdfgxxa*x)/(4.*omx*omxxa) - (9.*Logomxxa*pdfgxxa*x)/(4.*omx*omxxa) - 
     -  (3.*logmupapb2*Logomxxa*pdfgxxa*x)/(omx*omxxa) - (3.*logmupapi2*Logomxxa*pdfgxxa*x)/(omx*omxxa) + 
     -  (9.*Logomxxa**2*pdfgxxa*x)/(2.*omx*omxxa) + (8.*NF*pdfgxxa*x)/(9.*omx*omxxa) + 
     -  (logmupapb2*NF*pdfgxxa*x)/(6.*omx*omxxa) + (logmupapi2*NF*pdfgxxa*x)/(6.*omx*omxxa) + (4.*pdfqxa*x)/(3.*omx) - 
     -  (2.*logmupapb2*pdfqxa*x)/(3.*omx) - (2.*logmupapi2*pdfqxa*x)/(3.*omx) + (2.*Logomx*pdfqxa*x)/omx - 
     -  (4.*Logx*pdfqxa*x)/(3.*omx) + (3.*pdfgxxa*Pi**2*x)/(2.*omx*omxxa) + (3.*logmupapb2*pdfgxa*x**2)/omx + 
     -  (3.*logmupapi2*pdfgxa*x**2)/omx - (9.*Logomx*pdfgxa*x**2)/omx +
     -  (6.*Logx*pdfgxa*x**2)/omx

      end


      subroutine finitecollqgREM(xa,x,patilde,pb,pin,mufac,quark,result,pdfchange)

      implicit none

#include "VBFNLO/utilities/global.inc"


      integer quark
      integer pdfchange

      double precision xa,x,patilde(0:3),pb(0:3),pin(0:3),result,dotrr,
     &   xxa,pdfqxa,pdfgxa,pdfqxxa,pdfgxxa,logomx, logomxxa,logx,omx,omxxa,
     &   Logmupapb2,logmupapi2,NF,mufacsq,mufac,pdftmp(-6:6),one,inv2,zero
      parameter (nf=5d0,one=1d0,inv2=1d0/2d0,zero=0d0)

      external dotrr

      mufacsq=mufac**2
      xxa=x*xa

      omxxa=one-xxa
      omx=one-x

      logx=Log(x)
!       if (x.ge.0.999995) logx=0d0

      logomx=Log(omx)
!       if (x.le.0.000005) logomx=0d0

      logomxxa=Log(omxxa)
!       if (xxa.le.0.000005) logomxxa=0d0

      logmupapb2=Log(mufacsq*inv2/dotrr(patilde,pb))
      logmupapi2=Log(mufacsq*inv2/dotrr(patilde,pin))

      pdfqxa=zero
      pdfqxxa=zero

      call pdfproton( xa, mufac, pdftmp(-6) )
      pdfgxa=pdftmp(0)/xa
      pdfqxa=pdftmp(pdfchange*quark)/xa
         
      call pdfproton( xxa, mufac, pdftmp(-6) )    
      pdfgxxa=pdftmp(0)/xxa
      pdfqxxa=pdftmp(pdfchange*quark)/xxa


      result=(-9.*logmupapb2*pdfgxa)/(16.*omx) + (logmupapi2*pdfgxa)/(16.*omx) + (17.*Logomx*pdfgxa)/(16.*omx) - 
     -  (Logx*pdfgxa)/(2.*omx) + (4.*pdfqxa)/3. + pdfqxa/(4.*omx**2) - (3.*logmupapb2*pdfqxa)/(2.*omx**2) + 
     -  (logmupapi2*pdfqxa)/(6.*omx**2) + (17.*Logomx*pdfqxa)/(6.*omx**2) - (4.*Logx*pdfqxa)/(3.*omx**2) + pdfgxa*x + 
     -  (9.*logmupapb2*pdfgxa*x)/8. - (logmupapi2*pdfgxa*x)/8. - (17.*Logomx*pdfgxa*x)/8. + Logx*pdfgxa*x - 
     -  (pdfqxxa*x)/(4.*omx**2) + (3.*logmupapb2*pdfqxxa*x)/omx**2 - (logmupapi2*pdfqxxa*x)/(3.*omx**2) - 
     -  (17.*Logomx*pdfqxxa*x)/(3.*omx**2) - (77.*pdfqxxa*x)/(12.*omx*omxxa) - (9.*logmupapb2*pdfqxxa*x)/(4.*omx*omxxa) + 
     -  (logmupapi2*pdfqxxa*x)/(4.*omx*omxxa) + (Logomxxa*pdfqxxa*x)/(4.*omx*omxxa) - 
     -  (3.*logmupapb2*Logomxxa*pdfqxxa*x)/(omx*omxxa) + (logmupapi2*Logomxxa*pdfqxxa*x)/(3.*omx*omxxa) + 
     -  (17.*Logomxxa**2*pdfqxxa*x)/(6.*omx*omxxa) + (7.*pdfqxxa*Pi**2*x)/(18.*omx*omxxa) - 
     -  (3.*logmupapb2*pdfqxa*x**2)/(2.*omx**2) + (logmupapi2*pdfqxa*x**2)/(6.*omx**2) + 
     -  (17.*Logomx*pdfqxa*x**2)/(6.*omx**2) - (4.*Logx*pdfqxa*x**2)/(3.*omx**2)


      end


      subroutine calcgg (p,v,selc,m2s_qqwgagluj)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "tensorlep.inc"

      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)

      integer i,mu,isig1,isig2,isig3

      integer sign(4),polmap(2)

      double precision p(0:3,max_p), v(0:3,max_v),m2s_qqwgagluj
      double precision pbar(0:3,4),
     &       fermpphot(0:4,4),
     &       qbar(0:3,4),gluvec(0:4,2),gluvecb(0:4,2),glumom(0:4)

      double precision fermpg(0:4,4),fermpg17(0:4),fermpg16(0:4),fermpg15(0:4),fermpg27(0:4),
     &       fermpg25(0:4),fermpg26(0:4)

      double complex psi(2,-1:1,4),
     &       jqqggg(0:5,-1:1,8),
     &       fourgbcurr(0:3,2),ubdgg(2),jqqgggtmp(0:5,-1:1,11),jqqgggoff(0:5,-1:1,8),
     &       psig17(2,-1:1),psig176(2,-1:1),psig1765(2,-1:1),psig175(2,-1:1),psig1756(2,-1:1),
     &       psig16(2,-1:1),psig167(2,-1:1),psig1675(2,-1:1),psig165(2,-1:1),psig1657(2,-1:1),
     &       psig15(2,-1:1),psig157(2,-1:1),psig1576(2,-1:1),psig156(2,-1:1),psig1567(2,-1:1),
     &       psig25(2,-1:1),psig26(2,-1:1),psig27(2,-1:1),psig167g(2,-1:1),psig167g5(2,-1:1),
     &       psig1567g(2,-1:1),
     &       psig276(2,-1:1),psig2765(2,-1:1),psig275(2,-1:1),psig2756(2,-1:1),
     &       psig267(2,-1:1),psig2675(2,-1:1),psig265(2,-1:1),psig2657(2,-1:1),
     &       psig257(2,-1:1),psig2576(2,-1:1),psig256(2,-1:1),psig2567(2,-1:1),
     &       psig267g(2,-1:1),psig267g5(2,-1:1),
     &       psig2567g(2,-1:1),
     &       effgluhel(0:3)

      logical ldebug
      parameter (ldebug=.false.)
      double complex dotcc

      double complex temp,dotrc


      integer selc

      m2s_qqwgagluj=0d0
!
! Get information about different subprocesses and calculate external quark spinors

      Select Case (selc)

!-----------------------------------------------------
      Case (1) !   u~ d -> e- ve~ a g g

      sign(1)=1
      sign(2)=-1
      polmap(1)=1
      polmap(2)=1

      do mu=0,3
       pbar(mu,2)=p(mu,1)
       pbar(mu,1)=p(mu,2)
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

      Case (2) !   u~ g -> e- ve~ a d~ g

      sign(1)=-1
      sign(2)=-1
      polmap(1)=-1
      polmap(2)=1

      do mu=0,3
       pbar(mu,2)=p(mu,1)
       pbar(mu,1)=p(mu,4)
       gluvecb(mu,1)=p(mu,2)
       gluvecb(mu,2)=p(mu,3)
       gluvec(mu,1)=polmap(1)*gluvecb(mu,1)
       gluvec(mu,2)=polmap(2)*gluvecb(mu,2)
       glumom(mu)=gluvec(mu,1)+gluvec(mu,2)
       do i=1,2
         qbar(mu,i)=sign(i)*pbar(mu,i)
       enddo
      enddo

!-----------------------------------------------------

      Case (3) !   d g -> e- ve~ a u g

      sign(1)=1
      sign(2)=1
      polmap(1)=-1
      polmap(2)=1

      do mu=0,3
       pbar(mu,2)=p(mu,4)
       pbar(mu,1)=p(mu,1)
       gluvecb(mu,1)=p(mu,2)
       gluvecb(mu,2)=p(mu,3)
       gluvec(mu,1)=polmap(1)*gluvecb(mu,1)
       gluvec(mu,2)=polmap(2)*gluvecb(mu,2)
       glumom(mu)=gluvec(mu,1)+gluvec(mu,2)
       do i=1,2
         qbar(mu,i)=sign(i)*pbar(mu,i)
       enddo
      enddo

!-----------------------------------------------------

      Case (4) !   g g -> e- ve~ a u d~

      sign(1)=-1
      sign(2)=1
      polmap(1)=-1
      polmap(2)=-1

      do mu=0,3
       pbar(mu,2)=p(mu,4)
       pbar(mu,1)=p(mu,3)
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


!-----------------------------------------------------
!-----------------------------------------------------

      do isig2=-1,1,2       ! Loop over gluon (1) helicity
       do isig3 = -1,1,2    ! Loop over gluon (2) helicity
         do isig1= -1,1,2   ! Loop over photon helicity

      ubdgg(1)=0d0
      ubdgg(2)=0d0

       
      call psi0m(2,pbar(0,1),sign(1),psi)

      call helvec(gluvecb(0,1),polmap(1),isig2,glupolc(0,1,isig2))
      call helvec(gluvecb(0,2),polmap(2),isig3,glupolc(0,2,isig3))
!
!  Emitt three gauge bosons from the quark line before W
!  convention see my notes (MadGraph Code u~ d > e- ve~ a g(7=1) g(6=2)
!

!   UPPER CASE

!  g(1) g(2) a
        call bra2c(psi(1,-1,2),.true.,qbar(0,2),-1,gluvec(0,1),glupolc(0,1,isig2),psig17(1,-1),fermpg17(0))
        call bra2c(psig17(1,-1),.false.,fermpg17(0),-1,gluvec(0,2) ,glupolc(0,2,isig3),psig176(1,-1),fermpg(0,2))
        call bra2c(psig176(1,-1),.false.,fermpg(0,2),-1,photfourvec(0),wafpolc(0,isig1),psig1765(1,-1),fermpg(0,3))
        call curr6(-1,psig1765(1,-1),fermpg(0,3),psi(1,-1,1),qbar(0,1),jqqggg(0,-1,1))
!  g(2) g(1) a
        call bra2c(psi(1,-1,2),.true.,qbar(0,2),-1,gluvec(0,2),glupolc(0,2,isig3),psig16(1,-1),fermpg16(0))
        call bra2c(psig16(1,-1),.false.,fermpg16(0),-1,gluvec(0,1),glupolc(0,1,isig2),psig167(1,-1),fermpg(0,2))
        call bra2c(psig167(1,-1),.false.,fermpg(0,2),-1,photfourvec(0),wafpolc(0,isig1),psig1675(1,-1),fermpg(0,3))
        call curr6(-1,psig1675(1,-1),fermpg(0,3),psi(1,-1,1),qbar(0,1),jqqggg(0,-1,2))
!  g(1) a g(2)
        call bra2c(psig17(1,-1),.false.,fermpg17(0),-1,photfourvec(0),wafpolc(0,isig1),psig175(1,-1),fermpg(0,2))
        call bra2c(psig175(1,-1),.false.,fermpg(0,2),-1,gluvec(0,2),glupolc(0,2,isig3),psig1756(1,-1),fermpg(0,3))
        call curr6(-1,psig1756(1,-1),fermpg(0,3),psi(1,-1,1),qbar(0,1),jqqggg(0,-1,3))
!  g(2) a g(1)
        call bra2c(psig16(1,-1),.false.,fermpg16(0),-1,photfourvec(0),wafpolc(0,isig1),psig165(1,-1),fermpg(0,2))
        call bra2c(psig165(1,-1),.false.,fermpg(0,2),-1,gluvec(0,1),glupolc(0,1,isig2),psig1657(1,-1),fermpg(0,3))
        call curr6(-1,psig1657(1,-1),fermpg(0,3),psi(1,-1,1),qbar(0,1),jqqggg(0,-1,4))
! a g(1) g(2)
        call bra2c(psi(1,-1,2),.true.,qbar(0,2),-1,photfourvec(0),wafpolc(0,isig1),psig15(1,-1),fermpg15(0))
        call bra2c(psig15(1,-1),.false.,fermpg15(0),-1,gluvec(0,1),glupolc(0,1,isig2),psig157(1,-1),fermpg(0,2))
        call bra2c(psig157(1,-1),.false.,fermpg(0,2),-1,gluvec(0,2),glupolc(0,2,isig3),psig1576(1,-1),fermpg(0,3))
        call curr6(-1,psig1576(1,-1),fermpg(0,3),psi(1,-1,1),qbar(0,1),jqqggg(0,-1,5))
! a g(2) g(1)
        call bra2c(psig15(1,-1),.false.,fermpg15(0),-1,gluvec(0,2),glupolc(0,2,isig3),psig156(1,-1),fermpg(0,2))
        call bra2c(psig156(1,-1),.false.,fermpg(0,2),-1,gluvec(0,1),glupolc(0,1,isig2),psig1567(1,-1),fermpg(0,3))
        call curr6(-1,psig1567(1,-1),fermpg(0,3),psi(1,-1,1),qbar(0,1),jqqggg(0,-1,6))

      do mu=0,3
        fourgbcurr(mu,1)=jqqggg(mu,-1,1)+jqqggg(mu,-1,3)+jqqggg(mu,-1,5)
        fourgbcurr(mu,2)=jqqggg(mu,-1,2)+jqqggg(mu,-1,4)+jqqggg(mu,-1,6)
      enddo

!
!  Contract to get amplitude
!
        temp=dotcc(fourgbcurr(0,1),currlepw(0))
        ubdgg(1)=ubdgg(1)+temp*clr(1,3,-1)*clr(3,1,-1)*gg(1)**2

        temp=dotcc(fourgbcurr(0,2),currlepw(0))
        ubdgg(2)=ubdgg(2)+temp*clr(1,3,-1)*clr(3,1,-1)*gg(1)**2

!
!  Emitt three gauge bosons from the quark line before W
!
!       LOWER CASE

!  g(1) g(2) a
        call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,gluvec(0,1),glupolc(0,1,isig2),psig27(1,-1),fermpg27(0))
        call ket2c(psig27(1,-1),.false.,fermpg27(0),-1,gluvec(0,2) ,glupolc(0,2,isig3),psig276(1,-1),fermpg(0,2))
        call ket2c(psig276(1,-1),.false.,fermpg(0,2),-1,photfourvec(0),wafpolc(0,isig1),psig2765(1,-1),fermpg(0,3))
        call curr6(-1,psi(1,-1,2),qbar(0,2),psig2765(1,-1),fermpg(0,3),jqqggg(0,-1,1))
!  g(2) g(1) a
        call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,gluvec(0,2),glupolc(0,2,isig3),psig26(1,-1),fermpg26(0))
        call ket2c(psig26(1,-1),.false.,fermpg26(0),-1,gluvec(0,1),glupolc(0,1,isig2),psig267(1,-1),fermpg(0,2))
        call ket2c(psig267(1,-1),.false.,fermpg(0,2),-1,photfourvec(0),wafpolc(0,isig1),psig2675(1,-1),fermpg(0,3))
        call curr6(-1,psi(1,-1,2),qbar(0,2),psig2675(1,-1),fermpg(0,3),jqqggg(0,-1,2))
!  g(1) a g(2)
        call ket2c(psig27(1,-1),.false.,fermpg27(0),-1,photfourvec(0),wafpolc(0,isig1),psig275(1,-1),fermpg(0,2))
        call ket2c(psig275(1,-1),.false.,fermpg(0,2),-1,gluvec(0,2),glupolc(0,2,isig3),psig2756(1,-1),fermpg(0,3))
        call curr6(-1,psi(1,-1,2),qbar(0,2),psig2756(1,-1),fermpg(0,3),jqqggg(0,-1,3))
!  g(2) a g(1)
        call ket2c(psig26(1,-1),.false.,fermpg26(0),-1,photfourvec(0),wafpolc(0,isig1),psig265(1,-1),fermpg(0,2))
        call ket2c(psig265(1,-1),.false.,fermpg(0,2),-1,gluvec(0,1),glupolc(0,1,isig2),psig2657(1,-1),fermpg(0,3))
        call curr6(-1,psi(1,-1,2),qbar(0,2),psig2657(1,-1),fermpg(0,3),jqqggg(0,-1,4))
! a g(1) g(2)
        call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,photfourvec(0),wafpolc(0,isig1),psig25(1,-1),fermpg25(0))
        call ket2c(psig25(1,-1),.false.,fermpg25(0),-1,gluvec(0,1),glupolc(0,1,isig2),psig257(1,-1),fermpg(0,2))
        call ket2c(psig257(1,-1),.false.,fermpg(0,2),-1,gluvec(0,2),glupolc(0,2,isig3),psig2576(1,-1),fermpg(0,3))
        call curr6(-1,psi(1,-1,2),qbar(0,2),psig2576(1,-1),fermpg(0,3),jqqggg(0,-1,5))
! a g(2) g(1)
        call ket2c(psig25(1,-1),.false.,fermpg25(0),-1,gluvec(0,2),glupolc(0,2,isig3),psig256(1,-1),fermpg(0,2))
        call ket2c(psig256(1,-1),.false.,fermpg(0,2),-1,gluvec(0,1),glupolc(0,1,isig2),psig2567(1,-1),fermpg(0,3))
        call curr6(-1,psi(1,-1,2),qbar(0,2),psig2567(1,-1),fermpg(0,3),jqqggg(0,-1,6))

      do mu=0,3
        fourgbcurr(mu,2)=jqqggg(mu,-1,1)+jqqggg(mu,-1,3)+jqqggg(mu,-1,5)
        fourgbcurr(mu,1)=jqqggg(mu,-1,2)+jqqggg(mu,-1,4)+jqqggg(mu,-1,6)
      enddo

!
!  Contract to get amplitude
!
        temp=dotcc(fourgbcurr(0,1),currlepw(0))
        ubdgg(1)=ubdgg(1)+temp*clr(1,3,-1)*clr(4,1,-1)*gg(1)**2
        temp=dotcc(fourgbcurr(0,2),currlepw(0))
        ubdgg(2)=ubdgg(2)+temp*clr(1,3,-1)*clr(4,1,-1)*gg(1)**2

!
!  Emitt two gauge bosons from the quark line before W
!
!       UPPER CASE
!
!         call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,gluvec(0,1),glupolc(0,1,isig2),psig27(1,-1),fermpg(0,1))
!         call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,gluvec(0,2),glupolc(0,2,isig3),psig26(1,-1),fermpg(0,2))
!         call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,photfourvec(0),wafpolc(0,isig1),psig25(1,-1),fermpg(0,3))

        call curr6(-1,psig165(1,-1),fermpphot(0,1),psig27(1,-1),fermpg27(0),jqqggg(0,-1,1))
        call curr6(-1,psig175(1,-1),fermpphot(0,1),psig26(1,-1),fermpg26(0),jqqggg(0,-1,2))
        call curr6(-1,psig157(1,-1),fermpphot(0,1),psig26(1,-1),fermpg26(0),jqqggg(0,-1,3))
        call curr6(-1,psig156(1,-1),fermpphot(0,1),psig27(1,-1),fermpg27(0),jqqggg(0,-1,4))
        call curr6(-1,psig176(1,-1),fermpphot(0,1),psig25(1,-1),fermpg25(0),jqqggg(0,-1,5))
        call curr6(-1,psig167(1,-1),fermpphot(0,1),psig25(1,-1),fermpg25(0),jqqggg(0,-1,6))

      do mu=0,3
        fourgbcurr(mu,2)=(jqqggg(mu,-1,1)+jqqggg(mu,-1,4))*clr(3,1,-1)+jqqggg(mu,-1,6)*clr(4,1,-1)
        fourgbcurr(mu,1)=(jqqggg(mu,-1,2)+jqqggg(mu,-1,3))*clr(3,1,-1)+jqqggg(mu,-1,5)*clr(4,1,-1)
      enddo

!
!  Contract to get amplitude
!
        temp=dotcc(fourgbcurr(0,1),currlepw(0))
        ubdgg(1)=ubdgg(1)+temp*clr(1,3,-1)*gg(1)**2
        temp=dotcc(fourgbcurr(0,2),currlepw(0))
        ubdgg(2)=ubdgg(2)+temp*clr(1,3,-1)*gg(1)**2

!
!  Emitt two gauge bosons from the quark line before W
!
!       LOWER CASE
!
!         call bra2r(psi(1,-1,2),.true.,qbar(0,2),-1,gluvec(0,1),glupol(0,1,isig2),psig17(1,-1),fermpg(0,1))
!         call bra2r(psi(1,-1,2),.true.,qbar(0,2),-1,gluvec(0,2),glupol(0,2,isig3),psig16(1,-1),fermpg(0,2))
!         call bra2r(psi(1,-1,2),.true.,qbar(0,2),-1,photfourvec(0),wafpol(0,isig1),psig15(1,-1),fermpg(0,3))

        call curr6(-1,psig17(1,-1),fermpg17(0),psig265(1,-1),fermpphot(0,1),jqqggg(0,-1,1))
        call curr6(-1,psig16(1,-1),fermpg16(0),psig275(1,-1),fermpphot(0,1),jqqggg(0,-1,2))
        call curr6(-1,psig16(1,-1),fermpg16(0),psig257(1,-1),fermpphot(0,1),jqqggg(0,-1,3))
        call curr6(-1,psig17(1,-1),fermpg17(0),psig256(1,-1),fermpphot(0,1),jqqggg(0,-1,4))
        call curr6(-1,psig15(1,-1),fermpg15(0),psig276(1,-1),fermpphot(0,1),jqqggg(0,-1,5))
        call curr6(-1,psig15(1,-1),fermpg15(0),psig267(1,-1),fermpphot(0,1),jqqggg(0,-1,6))

      do mu=0,3
        fourgbcurr(mu,1)=(jqqggg(mu,-1,1)+jqqggg(mu,-1,4))*clr(4,1,-1)+jqqggg(mu,-1,6)*clr(3,1,-1)
        fourgbcurr(mu,2)=(jqqggg(mu,-1,2)+jqqggg(mu,-1,3))*clr(4,1,-1)+jqqggg(mu,-1,5)*clr(3,1,-1)
      enddo

!
!  Contract to get amplitude
!
        temp=dotcc(fourgbcurr(0,1),currlepw(0))
        ubdgg(1)=ubdgg(1)+temp*clr(1,3,-1)*gg(1)**2
        temp=dotcc(fourgbcurr(0,2),currlepw(0))
        ubdgg(2)=ubdgg(2)+temp*clr(1,3,-1)*gg(1)**2

! 
!   Get effective g(1)g(2) -> g polarization
!

        do mu=0,3
        effgluhel(mu) = (2d0*dotrc(gluvec(0,1),glupolc(0,2,isig3))*glupolc(mu,1,isig2)
     &       -2d0*dotrc(gluvec(0,2),glupolc(0,1,isig2))*glupolc(mu,2,isig3)
     &       +dotcc(glupolc(0,1,isig2),glupolc(0,2,isig3))*(gluvec(mu,2)-gluvec(mu,1)))
     &       /(glumom(4))
        enddo
!
!   three gauge bosons emitted from quark line, one is W -> e- ve~, other is g -> g g
!
!
!       UPPER CASE
!
! effglu a
        call bra2c(psi(1,-1,2),.true.,qbar(0,2),-1,glumom(0),effgluhel(0),psig167g(1,-1),fermpg(0,1))
        call bra2c(psig167g(1,-1),.false.,fermpg(0,1),-1,photfourvec(0),wafpolc(0,isig1),psig167g5(1,-1),fermpg(0,3))
        call curr6(-1,psig167g5(1,-1),fermpg(0,3),psi(1,-1,1),qbar(0,1),jqqggg(0,-1,1))
! a effglu
        call bra2c(psi(1,-1,2),.true.,qbar(0,2),-1,photfourvec(0),wafpolc(0,isig1),psig15(1,-1),fermpg(0,1))
        call bra2c(psig15(1,-1),.false.,fermpg(0,1),-1,glumom(0),effgluhel(0),psig1567g(1,-1),fermpg(0,2))
        call curr6(-1,psig1567g(1,-1),fermpg(0,2),psi(1,-1,1),qbar(0,1),jqqggg(0,-1,2))

      do mu=0,3
        fourgbcurr(mu,1)=jqqggg(mu,-1,1)+jqqggg(mu,-1,2)
        fourgbcurr(mu,2)=-fourgbcurr(mu,1)
      enddo

!
!  Contract to get amplitude
!
        temp=dotcc(fourgbcurr(0,1),currlepw(0))
        ubdgg(1)=ubdgg(1)+temp*clr(1,3,-1)*clr(3,1,-1)*gg(1)**2
        temp=dotcc(fourgbcurr(0,2),currlepw(0))
        ubdgg(2)=ubdgg(2)+temp*clr(1,3,-1)*clr(3,1,-1)*gg(1)**2

!
!       LOWER CASE
!
! effglu a
        call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,glumom(0),effgluhel(0),psig267g(1,-1),fermpg(0,1))
        call ket2c(psig267g(1,-1),.false.,fermpg(0,1),-1,photfourvec(0),wafpolc(0,isig1),psig267g5(1,-1),fermpg(0,3))
        call curr6(-1,psi(1,-1,2),qbar(0,2),psig267g5(1,-1),fermpg(0,3),jqqggg(0,-1,1))
! a effglu
        call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,photfourvec(0),wafpolc(0,isig1),psig25(1,-1),fermpg(0,1))
        call ket2c(psig25(1,-1),.false.,fermpg(0,1),-1,glumom(0),effgluhel(0),psig2567g(1,-1),fermpg(0,2))
        call curr6(-1,psi(1,-1,2),qbar(0,2),psig2567g(1,-1),fermpg(0,2),jqqggg(0,-1,2))

      do mu=0,3
        fourgbcurr(mu,1)=jqqggg(mu,-1,1)+jqqggg(mu,-1,2)
        fourgbcurr(mu,2)=-fourgbcurr(mu,1)
      enddo

!
!  Contract to get amplitude
!
        temp=dotcc(fourgbcurr(0,1),currlepw(0))
        ubdgg(1)=ubdgg(1)+temp*clr(1,3,-1)*clr(4,1,-1)*gg(1)**2
        temp=dotcc(fourgbcurr(0,2),currlepw(0))
        ubdgg(2)=ubdgg(2)+temp*clr(1,3,-1)*clr(4,1,-1)*gg(1)**2

!
!  All remaining W -> e- ve~,  A W (gg)_eff
!
!
        call bra2c(psi(1,-1,2),.true.,qbar(0,2),-1,photfourvec(0),wafpolc(0,isig1),psig15(1,-1),fermpg(0,1))
        call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,glumom(0),effgluhel(0),psig267g(1,-1),fermpg(0,2))
        call curr6(-1,psig15(1,-1),fermpg(0,1),psig267g(1,-1),fermpg(0,2),jqqggg(0,-1,1))

        call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,photfourvec(0),wafpolc(0,isig1),psig25(1,-1),fermpg(0,2))
        call bra2c(psi(1,-1,2),.true.,qbar(0,2),-1,glumom(0),effgluhel(0),psig167g(1,-1),fermpg(0,1))
        call curr6(-1,psig167g(1,-1),fermpg(0,2),psig25(1,-1),fermpg(0,1),jqqggg(0,-1,2))

      do mu=0,3
        fourgbcurr(mu,1)=jqqggg(mu,-1,1)*clr(3,1,-1)+jqqggg(mu,-1,2)*clr(4,1,-1)
        fourgbcurr(mu,2)=-fourgbcurr(mu,1)
      enddo

!
!  Contract to get amplitude
!
        temp=dotcc(fourgbcurr(0,1),currlepw(0))
        ubdgg(1)=ubdgg(1)+temp*clr(1,3,-1)*gg(1)**2
        temp=dotcc(fourgbcurr(0,2),currlepw(0))
        ubdgg(2)=ubdgg(2)+temp*clr(1,3,-1)*gg(1)**2


!     UPPER CASE

      Call curr6(-1,psig176(1,-1),fermpphot(0,1),psi(1,-1,1),qbar(0,1),jqqggg(0,-1,1))
      Call curr6(-1,psig17(1,-1),fermpphot(0,1),psi(1,-1,1),qbar(0,1),jqqgggtmp(0,-1,1))
      Call curr6(-1,psig167(1,-1),fermpphot(0,1),psi(1,-1,1),qbar(0,1),jqqggg(0,-1,2))
      Call curr6(-1,psig16(1,-1),fermpphot(0,1),psi(1,-1,1),qbar(0,1),jqqgggtmp(0,-1,2))
      Call curr6(-1,psig167g(1,-1),fermpphot(0,1),psi(1,-1,1),qbar(0,1),jqqggg(0,-1,3))
      Call curr6(-1,psi(1,-1,2),qbar(0,2),psi(1,-1,1),qbar(0,1),jqqgggtmp(0,-1,3))


!     LOWER CASE

      Call curr6(-1,psi(1,-1,2),qbar(0,2),psig267(1,-1),fermpphot(0,1),jqqggg(0,-1,4))
      Call curr6(-1,psi(1,-1,2),qbar(0,2),psig26(1,-1),fermpphot(0,1),jqqgggtmp(0,-1,4))
      Call curr6(-1,psi(1,-1,2),qbar(0,2),psig276(1,-1),fermpphot(0,1),jqqggg(0,-1,5))
      Call curr6(-1,psi(1,-1,2),qbar(0,2),psig27(1,-1),fermpphot(0,1),jqqgggtmp(0,-1,5))
      Call curr6(-1,psi(1,-1,2),qbar(0,2),psig267g(1,-1),fermpphot(0,1),jqqggg(0,-1,6))


!     MIDDLE CASE

      Call curr6(-1,psig16(1,-1),fermpphot(0,1),psig27(1,-1),fermpphot(0,1),jqqggg(0,-1,7))
      Call curr6(-1,psig16(1,-1),fermpphot(0,1),psi(1,-1,1),qbar(0,1),jqqgggtmp(0,-1,8))
      Call curr6(-1,psi(1,-1,2),qbar(0,2),psig27(1,-1),fermpphot(0,1),jqqgggtmp(0,-1,9))

      Call curr6(-1,psig17(1,-1),fermpphot(0,1),psig26(1,-1),fermpphot(0,1),jqqggg(0,-1,8))
      Call curr6(-1,psig17(1,-1),fermpphot(0,1),psi(1,-1,1),qbar(0,1),jqqgggtmp(0,-1,10))
      Call curr6(-1,psi(1,-1,2),qbar(0,2),psig26(1,-1),fermpphot(0,1),jqqgggtmp(0,-1,11))

      do mu=0,3
      jqqgggoff(mu,-1,1) = jqqggg(mu,-1,1)
     &  +dotcc(jqqgggtmp(0,-1,1),glupolc(0,2,isig3))*wemnugamom(mu)/dcmplx(xm2(3),-xmg(3))
      jqqgggoff(mu,-1,2) = jqqggg(mu,-1,2)
     &  +dotcc(jqqgggtmp(0,-1,2),glupolc(0,1,isig2))*wemnugamom(mu)/dcmplx(xm2(3),-xmg(3))
      jqqgggoff(mu,-1,3) = jqqggg(mu,-1,3)
     &  +dotcc(jqqgggtmp(0,-1,3),effgluhel(0))*wemnugamom(mu)/dcmplx(xm2(3),-xmg(3))
      jqqgggoff(mu,-1,4) = jqqggg(mu,-1,4)
     &  -dotcc(jqqgggtmp(0,-1,4),glupolc(0,1,isig2))*wemnugamom(mu)/dcmplx(xm2(3),-xmg(3))
      jqqgggoff(mu,-1,5) = jqqggg(mu,-1,5)
     &  -dotcc(jqqgggtmp(0,-1,5),glupolc(0,2,isig3))*wemnugamom(mu)/dcmplx(xm2(3),-xmg(3))
      jqqgggoff(mu,-1,6) = jqqggg(mu,-1,6)
     &  -dotcc(jqqgggtmp(0,-1,3),effgluhel(0))*wemnugamom(mu)/dcmplx(xm2(3),-xmg(3))
      jqqgggoff(mu,-1,7)=jqqggg(mu,-1,7)
     &      -dotcc(jqqgggtmp(0,-1,8),glupolc(0,1,isig2))*wemnugamom(mu)/dcmplx(xm2(3),-xmg(3))
     &      +dotcc(jqqgggtmp(0,-1,9),glupolc(0,2,isig3))*wemnugamom(mu)/dcmplx(xm2(3),-xmg(3))
      jqqgggoff(mu,-1,8)=jqqggg(mu,-1,8)
     &      -dotcc(jqqgggtmp(0,-1,10),glupolc(0,2,isig3))*wemnugamom(mu)/dcmplx(xm2(3),-xmg(3))
     &      +dotcc(jqqgggtmp(0,-1,11),glupolc(0,1,isig2))*wemnugamom(mu)/dcmplx(xm2(3),-xmg(3))
      enddo

!
!  Contract to get amplitude
!
      do mu=0,3
        fourgbcurr(mu,1)=-jqqgggoff(mu,-1,8)-jqqgggoff(mu,-1,1)-jqqgggoff(mu,-1,3)
     &            -jqqgggoff(mu,-1,4)-jqqgggoff(mu,-1,6)
        fourgbcurr(mu,2)=-jqqgggoff(mu,-1,7)-jqqgggoff(mu,-1,2)+jqqgggoff(mu,-1,3)
     &            -jqqgggoff(mu,-1,5)+jqqgggoff(mu,-1,6)
      enddo

        temp=dotcc(fourgbcurr(0,1),currlepwga(0,isig1))
        ubdgg(1)=ubdgg(1)+temp*clr(1,3,-1)*gg(1)**2
        temp=dotcc(fourgbcurr(0,2),currlepwga(0,isig1))
        ubdgg(2)=ubdgg(2)+temp*clr(1,3,-1)*gg(1)**2


      m2s_qqwgagluj=m2s_qqwgagluj+16d0/3d0*(abs(ubdgg(1))**2+
     &       abs(ubdgg(2))**2) - 4d0/3d0*DBLE(ubdgg(1)*
     &       dconjg(ubdgg(2)))

! !----------------------------------------------------------------------------------
! Check
!       do mu=0,3
! !     ------------------------------------
!       helastestmom(mu,1)=p(mu,1)
! !     ------------------------------------
!       helastestmom(mu,2)=p(mu,2)
! !     ------------------------------------
!       helastestmom(mu,3)=v(mu,1)
! !     ------------------------------------
!       helastestmom(mu,4)=v(mu,2)
! !     ------------------------------------
!       helastestmom(mu,5)=v(mu,3)
! !     ------------------------------------
!       helastestmom(mu,6)=p(mu,4)
! !     ------------------------------------
!       helastestmom(mu,7)=p(mu,3)
! !     ------------------------------------
!       enddo
! 
!       do i=1,7
!       ipart(i)=1
!       enddo
! 
!       ihel(1)=1
!       ihel(2)=-1
!       ihel(3)=-1
!       ihel(4)=+1
!       ihel(5)=isig1
!       ihel(6)=isig3
!       ihel(7)=isig2
! 
!       call Subdgg(helastestmom,ihel,ipart,help1)
!       print*, "HELS   ", ihel(1),ihel(2),ihel(6),ihel(7)
!       print*, "HELAS  ", help1
!       print*, "My Code", m2s_qqwgagluj
!       print*, "Ratio", help1/ m2s_qqwgagluj
!       print*, "-------------------"
! 
!       STOP
! 
!-------------------------------------------------------------------

      enddo     !Close Photon helicity loop HELICITY
      enddo     !Close Gluon 2 helicity loop HELICITY
      enddo     !Close Gluon 1 helicity loop HELICITY


      end


      subroutine calcuu (p,v,selc,m2s_qqwgagluj)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "tensorlep.inc"

      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)

      integer j,i,k,l,mu,isig1,isig2,ihel(7),ipart(7),chi

      integer sign(4),iflavu,iflavd

      double precision p(0:3,max_p), v(0:3,max_v),m2s_qqwgagluj(3)
      double precision Q2(0:3),pbar(0:3,4),
     &       fermpphot(0:4,4),fermpphotw(0:4,8),
     &       Q2a(0:3),
     &       Q2uza(0:3),Q2dza(0:3),qbar(0:3,4),
     &       fermpphotww(0:4,4),fermpphotwa(0:4,6)
      double complex psi(2,-1:1,4), jqqu(0:5,-1:1,2), jqqd(0:5,-1:1,2),
     &       psiwmaz(2,-1:1,4),
     &       ubububdb(2,-1:1,-1:1,-1:1,-1:1),
     &       jqqza(0:5,-1:1,2),jqqaz(0:5,-1:1,2),
     &       psiww(2,-1:1,4),jqqvw(0:5,-1:1,4),propwa(2),
     &       psiwwa(2,-1:1,8),jqqwwa(0:5,-1:1,7),
     &       jqqwaa(0:5,-1:1,20),psiwaa(2,-1:1,6),propa(2),jqqwaaoff(0:5,-1:1,4)

      logical ldebug
      parameter (ldebug=.true.)
    
      double precision helastestmom(0:3,7)
      double precision m2s_qqwgaglu,m2s_wgaglu
      double complex contract_Tjj,dotcc
      double precision dotrr
      double complex temp,dotrc
      double complex temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8

      integer selc
     
      external m2s_qqwgaglu,m2s_wgaglu,contract_Tjj,dotcc,dotrc,dotrr


      m2s_qqwgagluj(1)=0d0
      m2s_qqwgagluj(2)=0d0
      m2s_qqwgagluj(3)=0d0
!
      
! We lay out a methode similiar to the "vbfnlo-approach", i.e. we calculate subtopologies using HELAS.
! Afterwards they are dotted with currents calclated the Zeppenfeld-Hagiwara-way

!  the final state down-quark has momentum p(.,3)
!  the final state up-quark has momentum p(.,4)

!  SUBROCESSES OF TYPE
!
! u~ u~ -> e- ve~ a u~ d~

      Select Case (selc)
!     to identify subprocess from crossing
      Case (1)  ! u~ u~ > e- ve~ a u~ d~
!
      iflavu=3
      iflavd=4
      sign(1)=-1
      sign(2)=-1  
      sign(3)=-1
      sign(4)=-1
      do mu=0,3
      pbar(mu,1)=p(mu,3)   ! 1-> ket
      pbar(mu,2)=p(mu,1)   ! 2-> bra
      pbar(mu,3)=p(mu,4)   ! 3-> ket
      pbar(mu,4)=p(mu,2)   ! 4-> bra
      enddo
!
      Case (2)  ! u~ u > e- ve~ a u d~
!
      iflavu=3
      iflavd=4
      sign(1)=-1
      sign(2)=-1  
      sign(3)=1
      sign(4)=1
      do mu=0,3
      pbar(mu,1)=p(mu,3)   ! 1-> ket
      pbar(mu,2)=p(mu,1)   ! 2-> bra
      pbar(mu,3)=p(mu,2)   ! 3-> ket
      pbar(mu,4)=p(mu,4)   ! 4-> bra
      enddo

!
      Case (3)  ! u~ d > e- ve~ a u u~
!
      iflavu=3
      iflavd=4
      sign(1)=1
      sign(2)=1
      sign(3)=-1
      sign(4)=-1
      do mu=0,3
      pbar(mu,1)=p(mu,2)   ! 1-> ket  !changed p(.,3) <-> p(.,4) for dipoles!
      pbar(mu,2)=p(mu,3)   ! 2-> bra
      pbar(mu,3)=p(mu,4)   ! 3-> ket
      pbar(mu,4)=p(mu,1)   ! 4-> bra
      enddo
!
      Case (4)  ! u d > e- ve~ a u u
!
      iflavu=3
      iflavd=4
      sign(1)=1
      sign(2)=1
      sign(3)=1
      sign(4)=1
      do mu=0,3
      pbar(mu,1)=p(mu,2)   ! 1-> ket
      pbar(mu,2)=p(mu,4)   ! 2-> bra
      pbar(mu,3)=p(mu,1)   ! 3-> ket
      pbar(mu,4)=p(mu,3)   ! 4-> bra
      enddo

      END SELECT

      call psi0m(4,pbar(0,1),sign(1),psi)

      do mu=0,3
       do i=1,4
         qbar(mu,i)=sign(i)*pbar(mu,i)
       enddo
      enddo

      do mu=0,3
      Q2(mu)= -qbar(mu,4)+qbar(mu,3)!p(mu,2)-p(mu,4)
      Q2a(mu)=-qbar(mu,2)+qbar(mu,3)!p(mu,1)-p(mu,4)
      Q2uza(mu)=-qbar(mu,2)+qbar(mu,3)-photfourvec(mu)!p(mu,1)-p(mu,4)-v(mu,3)
      Q2dza(mu)=-qbar(mu,4)+qbar(mu,3)-photfourvec(mu)!p(mu,2)-p(mu,4)-v(mu,3)
      enddo
      propwa(1)=1d0/(dotrr(q2dza,q2dza))
      propwa(2)=1d0/(dotrr(q2uza,q2uza))
      propa(1)=1d0/dotrr(q2,q2)
      propa(2)=1d0/dotrr(q2a,q2a)

!
! get the f-fbar currents
!
      call curr6(1,psi(1,-1,2),qbar(0,2),psi(1,-1,1),qbar(0,1),jqqu(0,-1,1))
      call curr6(1,psi(1,-1,4),qbar(0,4),psi(1,-1,3),qbar(0,3),jqqu(0,-1,2))

      call curr6(1,psi(1,-1,4),qbar(0,4),psi(1,-1,1),qbar(0,1),jqqd(0,-1,2))
      call curr6(1,psi(1,-1,2),qbar(0,2),psi(1,-1,3),qbar(0,3),jqqd(0,-1,1))   

!------------------------------------------------------------
!------------------------------------------------------------
      do isig1=-1,1,2 ! sum over external photon helicity   +
!------------------------------------------------------------
!------------------------------------------------------------
      do i=-1,1,2
       do j=-1,1,2
        do k=-1,1,2
         do l=-1,1,2
!                   u u d u
         ubububdb(1,i,j,k,l)=dcmplx(0d0,0d0)
         ubububdb(2,i,j,k,l)=dcmplx(0d0,0d0)
         enddo
        enddo
       enddo
      enddo


      if(ldebug) then

      do mu=0,3
      helastestmom(mu,1)=p(mu,1)
      helastestmom(mu,2)=p(mu,2)
      helastestmom(mu,3)=v(mu,1)
      helastestmom(mu,4)=v(mu,2)
      helastestmom(mu,5)=v(mu,3)
      helastestmom(mu,6)=p(mu,4)
      helastestmom(mu,7)=p(mu,3)
      enddo

      do k=1,7
      ipart(k)=1
      enddo

      ihel(1)=1
      ihel(2)=1
      ihel(6)=1     !QUARKS
      ihel(7)=1
!------------------------
      ihel(3)=-1
      ihel(4)=1
      ihel(5)=isig1
!------------------------
      endif


       do chi=-1,1,2

!  photon radiation off the quark legs

        call bra2c(psi(1,chi,2),.true.,qbar(0,2),chi,photfourvec(0),wafpolc(0,isig1),
     &              psiwmaz(1,chi,2),fermpphot(0,2))
        call ket2c(psi(1,chi,1),.true.,qbar(0,1),chi,photfourvec(0),wafpolc(0,isig1),
     &              psiwmaz(1,chi,1),fermpphot(0,1))
        call bra2c(psi(1,chi,4),.true.,qbar(0,4),chi,photfourvec(0),wafpolc(0,isig1),
     &              psiwmaz(1,chi,4),fermpphot(0,4))
        call ket2c(psi(1,chi,3),.false.,qbar(0,3),chi,photfourvec(0),wafpolc(0,isig1),
     &              psiwmaz(1,chi,3),fermpphot(0,3))
        enddo

! form currents
!
        call curr6(1,psi(1,-1,2),qbar(0,2),psiwmaz(1,-1,3),fermpphot(0,3),jqqaz(0,-1,1))
        call curr6(1,psiwmaz(1,-1,2),fermpphot(0,2),psi(1,-1,3),qbar(0,3),jqqza(0,-1,1))

        call curr6(1,psi(1,-1,4),qbar(0,4),psiwmaz(1,-1,3),fermpphot(0,3),jqqaz(0,-1,2))
        call curr6(1,psiwmaz(1,-1,4),fermpphot(0,4),psi(1,-1,3),qbar(0,3),jqqza(0,-1,2))

! W radiation off the quark legs, subsequent leptonic decay

      call bra2c(psi(1,-1,4),.true.,qbar(0,4),-1,wemnumom(0),currlepw(0),
     &              psiww(1,-1,4),fermpphotww(0,4))       
      call bra2c(psi(1,-1,2),.true.,qbar(0,2),-1,wemnumom(0),currlepw(0),
     &              psiww(1,-1,2),fermpphotww(0,2))
      call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,wemnumom(0),currlepw(0),
     &              psiww(1,-1,1),fermpphotww(0,1))

      Call curr6(-1,psi(1,-1,2),qbar(0,2),psiww(1,-1,1),fermpphotww(0,1),jqqvw(0,-1,1))
      Call curr6(-1,psiww(1,-1,2),fermpphotww(0,2),psi(1,-1,1),qbar(0,1),jqqvw(0,-1,2))

      Call curr6(-1,psi(1,-1,4),qbar(0,4),psiww(1,-1,1),fermpphotww(0,1),jqqvw(0,-1,3))
      Call curr6(-1,psiww(1,-1,4),fermpphotww(0,4),psi(1,-1,1),qbar(0,1),jqqvw(0,-1,4))

      do isig2=-1,1,2
!
!  W A / W G  case
!
!  lower line GA
!
      temp1=dotcc(jqqvw(0,-1,1),jqqaz(0,isig2,2))
      temp2=dotcc(jqqvw(0,-1,1),jqqza(0,isig2,2))
      temp3=dotcc(jqqvw(0,-1,2),jqqaz(0,isig2,2))
      temp4=dotcc(jqqvw(0,-1,2),jqqza(0,isig2,2))

      temp=(temp1+temp2+temp3+temp4)*gg(1)**2*clr(1,3,-1)*clr(3,1,isig2)*propwa(1)

      ubububdb(1,-1,isig2,-1,isig2)=ubububdb(1,-1,isig2,-1,isig2)+temp


!  upper line GA

      temp1=dotcc(jqqvw(0,-1,3),jqqaz(0,isig2,1))
      temp2=dotcc(jqqvw(0,-1,3),jqqza(0,isig2,1))
      temp3=dotcc(jqqvw(0,-1,4),jqqaz(0,isig2,1))
      temp4=dotcc(jqqvw(0,-1,4),jqqza(0,isig2,1))

      temp=(temp1+temp2+temp3+temp4)*gg(1)**2*clr(1,3,-1)*clr(iflavu,1,isig2)*propwa(2)

      ubububdb(2,isig2,-1,-1,isig2)=ubububdb(2,-1,isig2,-1,isig2)+temp


      enddo


!  finally emit up to three gauge bosons from a single quark line

! upper case

      Call bra2c(psiww(1,-1,4),.false.,fermpphotww(0,4),-1,photfourvec(0),wafpolc(0,isig1),
     &           psiwwa(1,-1,1),fermpphotw(0,1))

      Call bra2c(psiww(1,-1,2),.false.,fermpphotww(0,2),-1,photfourvec(0),wafpolc(0,isig1),
     &           psiwwa(1,-1,3),fermpphotw(0,3))

      Call bra2c(psiwmaz(1,-1,2),.false.,fermpphot(0,2),-1,wemnumom(0),currlepw(0),
     &           psiwwa(1,-1,4),fermpphotw(0,4))

      call bra2c(psi(1,-1,2),.false.,qbar(0,2),-1,wemnugamom(0),currlepwga(0,isig1),
     &           psiwwa(1,-1,5),fermpphotw(0,5))

      call curr6(-1,psiwwa(1,-1,5),fermpphotw(0,5),psi(1,-1,3),qbar(0,3),jqqwwa(0,-1,7))

      Call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,wemnumom(0),currlepw(0),
     &           psiwaa(1,-1,1),fermpphotwa(0,1))

      Call curr6(-1,psiwmaz(1,-1,2),fermpphot(0,2),psiwaa(1,-1,1),fermpphotwa(0,1),jqqwaa(0,-1,1))
      Call curr6(-1,psiwwa(1,-1,4),fermpphotw(0,4),psi(1,-1,1),qbar(0,1),jqqwaa(0,-1,2))
      Call curr6(-1,psiww(1,-1,2),fermpphotww(0,2),psiwmaz(1,-1,1),fermpphot(0,1),jqqwaa(0,-1,3))
      Call curr6(-1,psiwwa(1,-1,3),fermpphotw(0,3),psi(1,-1,1),qbar(0,1),jqqwaa(0,-1,4))

      Call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,wemnugamom(0),currlepwga(0,isig1),
     &           psiwaa(1,-1,2),fermpphotwa(0,2))

      Call curr6(-1,psi(1,-1,2),qbar(0,2),psiwaa(1,-1,2),fermpphotwa(0,2),jqqwaa(0,-1,5))
      Call curr6(-1,psiwwa(1,-1,5),fermpphotw(0,5),psi(1,-1,1),qbar(0,1),jqqwaa(0,-1,6))

      do mu=0,3
      jqqwaaoff(mu,-1,1)=jqqwaa(mu,-1,5)+jqqu(mu,-1,1)*
     &              dotrc(wemnugamom,currlepwga(0,isig1))/dcmplx(xm2(3),-xmg(3))
      jqqwaaoff(mu,-1,2)=jqqwaa(mu,-1,6)-jqqu(mu,-1,1)*
     &              dotrc(wemnugamom,currlepwga(0,isig1))/dcmplx(xm2(3),-xmg(3))
      enddo

      Call ket2c(psiwmaz(1,-1,1),.false.,fermpphot(0,1),-1,wemnumom(0),currlepw(0),
     &           psiwaa(1,-1,3),fermpphotwa(0,3))
      Call curr6(-1,psi(1,-1,2),qbar(0,2),psiwaa(1,-1,3),fermpphotwa(0,3),jqqwaa(0,-1,7))

      Call ket2c(psiwaa(1,-1,1),.false.,fermpphotwa(0,1),-1,photfourvec(0),wafpolc(0,isig1),
     &           psiwaa(1,-1,4),fermpphotwa(0,4))
      Call curr6(-1,psi(1,-1,2),qbar(0,2),psiwaa(1,-1,4),fermpphotwa(0,4),jqqwaa(0,-1,8))

      do isig2=-1,1,2
      temp1= dotcc(jqqwaa(0,-1,1),jqqu(0,isig2,2))
      temp2= dotcc(jqqwaa(0,-1,2),jqqu(0,isig2,2))
      temp3= dotcc(jqqwaa(0,-1,3),jqqu(0,isig2,2))
      temp4= dotcc(jqqwaa(0,-1,4),jqqu(0,isig2,2))
      temp5= -dotcc(jqqwaaoff(0,-1,1),jqqu(0,isig2,2)) 
      temp6= -dotcc(jqqwaaoff(0,-1,2),jqqu(0,isig2,2))
      temp7= dotcc(jqqwaa(0,-1,7),jqqu(0,isig2,2))
      temp8= dotcc(jqqwaa(0,-1,8),jqqu(0,isig2,2))

      ubububdb(1,-1,isig2,-1,isig2)= ubububdb(1,-1,isig2,-1,isig2)+ 
     &  ( temp1*(clr(3,1,-1))
     &   +temp2*(clr(3,1,-1))
     &   +temp3*(clr(4,1,-1))
     &   +temp4*(clr(4,1,-1))
     &   +temp5
     &   +temp6
     &   +temp7*(clr(4,1,-1))
     &   +temp8*(clr(3,1,-1)))*propa(1)*clr(1,3,-1)*gg(1)**2

!        test2= (
!      &   temp1*(clr(3,1,-1))
!      &   +temp2*(clr(3,1,-1))
!      &   +temp3*(clr(4,1,-1))
!      &   +temp4*(clr(4,1,-1))
!      &   +temp5
!      &   +temp6
!      &   +temp7*(clr(4,1,-1))
!      &   +temp8*(clr(3,1,-1))
!      &   )*propa(1)*clr(1,3,-1)*gg(1)**2
! 
       enddo

!       ihel(1)=1
!       ihel(2)=-1
!       ihel(6)=-1     !QUARKS
!       ihel(7)=1
!       print*,ubububdb(1,-ihel(1),-ihel(2),-ihel(7),-ihel(6))
!       CALL Sreal(helastestmom,ihel,ipart,help1c)
!       print*, help1c
!       print*, ubububdb(1,-ihel(1),-ihel(2),-ihel(7),-ihel(6))/help1c
! 
!       STOP


!
!  lower case
!
      Call curr6(-1,psiwmaz(1,-1,4),fermpphot(0,4),psiwaa(1,-1,1),fermpphotwa(0,1),jqqwaa(0,-1,9))
      Call curr6(-1,psiwwa(1,-1,1),fermpphotw(0,1),psi(1,-1,1),qbar(0,1),jqqwaa(0,-1,10))
      Call curr6(-1,psiww(1,-1,4),fermpphotww(0,4),psiwmaz(1,-1,1),fermpphot(0,1),jqqwaa(0,-1,11))
      Call curr6(-1,psi(1,-1,4),qbar(0,4),psiwaa(1,-1,2),fermpphotwa(0,2),jqqwaa(0,-1,12))

      call bra2c(psi(1,-1,4),.false.,qbar(0,4),-1,wemnugamom(0),currlepwga(0,isig1),
     &           psiwwa(1,-1,8),fermpphotw(0,8))

      Call curr6(-1,psiwwa(1,-1,8),fermpphotw(0,8),psi(1,-1,1),qbar(0,1),jqqwaa(0,-1,13))

      do mu=0,3
      jqqwaaoff(mu,-1,3)=jqqwaa(mu,-1,12)+jqqd(mu,-1,2)*
     &              dotrc(wemnugamom,currlepwga(0,isig1))/dcmplx(xm2(3),-xmg(3))
      jqqwaaoff(mu,-1,4)=jqqwaa(mu,-1,13)-jqqd(mu,-1,2)*
     &              dotrc(wemnugamom,currlepwga(0,isig1))/dcmplx(xm2(3),-xmg(3))
      enddo

      Call bra2c(psiwmaz(1,-1,4),.false.,fermpphot(0,4),-1,wemnumom(0),currlepw(0),
     &           psiwaa(1,-1,5),fermpphotwa(0,5))

      Call curr6(-1,psiwaa(1,-1,5),fermpphotwa(0,5),psi(1,-1,1),qbar(0,1),jqqwaa(0,-1,14))
      Call curr6(-1,psi(1,-1,4),qbar(0,4),psiwaa(1,-1,3),fermpphotwa(0,3),jqqwaa(0,-1,15))
      Call curr6(-1,psi(1,-1,4),qbar(0,4),psiwaa(1,-1,4),fermpphotwa(0,4),jqqwaa(0,-1,16))

      do isig2=-1,1,2
      temp1= dotcc(jqqwaa(0,-1,9),jqqd(0,isig2,1))
      temp2= dotcc(jqqwaa(0,-1,10),jqqd(0,isig2,1))
      temp3= dotcc(jqqwaa(0,-1,11),jqqd(0,isig2,1))
      temp4= dotcc(jqqwaa(0,-1,16),jqqd(0,isig2,1))
      temp5= -dotcc(jqqwaaoff(0,-1,3),jqqd(0,isig2,1))
      temp6= -dotcc(jqqwaaoff(0,-1,4),jqqd(0,isig2,1))
      temp7= +dotcc(jqqwaa(0,-1,14),jqqd(0,isig2,1))
      temp8= +dotcc(jqqwaa(0,-1,15),jqqd(0,isig2,1))

      ubububdb(2,isig2,-1,-1,isig2)= ubububdb(2,isig2,-1,-1,isig2)+
     & ( temp1*(clr(iflavu,1,-1))
     &  +temp2*(clr(iflavd,1,-1))
     &  +temp3*(clr(iflavd,1,-1))
     &  +temp4*(clr(iflavu,1,-1))
     &  +temp5
     &  +temp6
     &  +temp7*(clr(iflavu,1,-1))
     &  +temp8*(clr(iflavd,1,-1)))*clr(1,3,-1)*propa(2)*gg(1)**2

      enddo

! MATRIXELEMENT

      do i=-1,1,2
       do j=-1,1,2
        do k=-1,1,2
         do l=-1,1,2
           m2s_qqwgagluj(1)=m2s_qqwgagluj(1)
     &              +2d0*abs(ubububdb(1,i,j,k,l))**2
     &              +2d0*abs(ubububdb(2,i,j,k,l))**2
     &              +4d0/3d0*DBLE(ubububdb(1,i,j,k,l)*dconjg(ubububdb(2,i,j,k,l)))
           m2s_qqwgagluj(2)=m2s_qqwgagluj(2)
     &              +2d0*abs(ubububdb(2,i,j,k,l))**2
           m2s_qqwgagluj(3)=m2s_qqwgagluj(3)
     &              +2d0*abs(ubububdb(1,i,j,k,l))**2
         enddo
        enddo
       enddo
      enddo

! ----------------------------------------------------------------------#+
! CHECK
! ! ----------------------------------------------------------------------#+
!       test1=0d0
!       help2=0d0
! 
!       do i=-1,1,2
!        do j=-1,1,2
!         do k=-1,1,2
!          do l=-1,1,2
! 
!       ihel(1)=i
!       ihel(2)=j
!       ihel(7)=k
!       ihel(6)=l
!       CALL Suduu(helastestmom,ihel,ipart,help1)!C2
!       help2=+abs(help1c)+help2	
!       test1=0d0! test1
!      &              +2d0*abs(ubububdb(1,i,j,k,l))**2
!      &              +2d0*abs(ubububdb(2,i,j,k,l))**2
!      &              +2d0/3d0*2d0*DBLE(ubububdb(1,i,j,k,l)*dconjg(ubububdb(2,i,j,k,l)))
! 
! 
!       PRINT*, "My Code ", test1
!       Print*, i,j,k,l
!       PRINT*, "HELAS   ", help1
!       print*, ihel(1),ihel(2),ihel(3),ihel(4)
!       PRINT*, "Fraction", test1/help1
!       Print*, "-------------------------"
!          enddo
!         enddo
!        enddo
!       enddo
! 
!       STOP
! 
!       print*, test1/36d0
!       CAll STESTSQUBUBUBDB(helastestmom,help1)
!       print*, test1/36d0,help1,help1/test1*36d0
! ! 
! ! 
!       STOP
!       elseif (selid.eq.2) then
! 
!       do i=-1,1,2
!        do j=-1,1,2
!         do k=-1,1,2
!          do l=-1,1,2
!            help1=help1
!      &              +2d0*abs(ubububdb(2,i,j,k,l))**2
! 
!          enddo
!         enddo
!        enddo
!       enddo

!       do i=-1,1,2
!        do j=-1,1,2
!         do k=-1,1,2
!          do l=-1,1,2
! 
!       ihel(1)=i
!       ihel(2)=j
!       ihel(7)=k
!       ihel(6)=l
!       CALL SUBSUBC(helastestmom,ihel,ipart,help1)!C2
!       test1=
! !      &              +2d0*abs(ubububdb(1,i,j,k,l))**2
!      &              +2d0*abs(ubububdb(1,i,j,k,l))**2
! !      &              +2d0/3d0*2d0*DBLE(ubububdb(1,i,j,k,l)*dconjg(ubububdb(2,i,j,k,l)))
! 
! 
!       PRINT*, "My Code ", test1
!       Print*, i,j,k,l
!       PRINT*, "HELAS   ", help1
!       print*, ihel(1),ihel(2),ihel(3),ihel(4)
!       PRINT*, "Fraction", abs(help1c)/test1
!       Print*, "-------------------------"
!          enddo
!         enddo
!        enddo
!       enddo
! 
!       STOP
! 
!       else
! 
!       do i=-1,1,2
!        do j=-1,1,2
!         do k=-1,1,2
!          do l=-1,1,2
!            m2s_qqwgagluj=m2s_qqwgagluj
!      &              +2d0*abs(ubububdb(1,i,j,k,l))**2
! 
!          enddo
!         enddo
!        enddo
!       enddo
! 
!       do i=-1,1,2
!        do j=-1,1,2
!         do k=-1,1,2
!          do l=-1,1,2
! 
!       ihel(1)=-i
!       ihel(2)=-j
!       ihel(7)=-k
!       ihel(6)=-l
!       CALL  SUSUC(helastestmom,ihel,ipart,help1c)!C2
!       test1=
! !      &              +2d0*abs(ubububdb(1,i,j,k,l))**2
!      &              +2d0*abs(ubububdb(1,i,j,k,l))**2
! !      &              +2d0/3d0*2d0*DBLE(ubububdb(1,i,j,k,l)*dconjg(ubububdb(2,i,j,k,l)))
! 
! 
!       PRINT*, "My Code ", test1
!       Print*, i,j,k,l
!       PRINT*, "HELAS   ", abs(help1c)
!       print*, ihel(1),ihel(2),ihel(3),ihel(4)
!       PRINT*, "Fraction", abs(help1c)/test1
!       Print*, "-------------------------"
!          enddo
!         enddo
!        enddo
!       enddo
! 
! 
!       endif
! 
! 
!       STOP
! 
      enddo ! close loop over external photon helicity

      end


      subroutine calcdd(p,v,selc,m2s_qqwgaglujdd)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "tensorlep.inc"

      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)

      integer j,i,k,l,mu,isig1,isig2,ihel(7),ipart(7),chi

      integer sign(4)

      double precision p(0:3,max_p), v(0:3,max_v),m2s_qqwgaglujdd(3)
      double precision Q2(0:3),pbar(0:3,4),
     &       fermpphot(0:4,4),fermpphotw(0:4,8),
     &       Q2a(0:3),
     &       Q2uza(0:3),Q2dza(0:3),qbar(0:3,4),
     &       fermpphotww(0:4,4),fermpphotwa(0:4,6)
      double complex psi(2,-1:1,4), jqqu(0:5,-1:1,2), jqqd(0:5,-1:1,2),
     &       psiwmaz(2,-1:1,4),
     &       ubububdb(2,-1:1,-1:1,-1:1,-1:1),
     &       jqqza(0:5,-1:1,2),jqqaz(0:5,-1:1,2),
     &       psiww(2,-1:1,4),jqqvw(0:5,-1:1,4),propwa(2),
     &       psiwwa(2,-1:1,8),
     &       jqqwaa(0:5,-1:1,20),psiwaa(2,-1:1,6),propa(2),jqqwaaoff(0:5,-1:1,4)

      logical ldebug
      parameter (ldebug=.true.)
    
      double precision helastestmom(0:3,7)
      double precision m2s_qqwgaglu,m2s_wgaglu
      double complex contract_Tjj,dotcc
      double precision dotrr
      double complex temp,dotrc
      double complex temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8

      integer selc
     
      external m2s_qqwgaglu,m2s_wgaglu,contract_Tjj,dotcc,dotrc,dotrr


      m2s_qqwgaglujdd(1)=0d0
      m2s_qqwgaglujdd(2)=0d0
      m2s_qqwgaglujdd(3)=0d0
!

      
! We lay out a methode similiar to the "vbfnlo-approach", i.e. we calculate subtopologies using HELAS.
! Afterwards they are dotted with currents calclated the Zeppenfeld-Hagiwara-way

!  the final state down-quark has momentum p(.,3)
!  the final state up-quark has momentum p(.,4)

!  SUBROCESSES OF TYPE
!
! d d -> e- ve~ a d u

!

      SELECT CASE (selc)

      CASE (1)  ! d d -> e- ve~ a d u

      sign(1)=1
      sign(2)=1  
      sign(3)=1
      sign(4)=1
      do mu=0,3
      pbar(mu,1)=p(mu,1)   ! 1-> ket
      pbar(mu,2)=p(mu,3)   ! 2-> bra
      pbar(mu,3)=p(mu,2)   ! 3-> ket
      pbar(mu,4)=p(mu,4)   ! 4-> bra
      enddo

      CASE (2)  ! d~ u~ -> e- ve~ a d~ d~

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

      CASE (3)  ! u~ d -> e- ve~ a d d~

      sign(1)=-1
      sign(2)=-1
      sign(3)=1
      sign(4)=1
      do mu=0,3
      pbar(mu,1)=p(mu,3)   ! 1-> ket
      pbar(mu,2)=p(mu,1)   ! 2-> bra
      pbar(mu,3)=p(mu,2)   ! 3-> ket
      pbar(mu,4)=p(mu,4)   ! 4-> bra
      enddo

      CASE (4)  ! d~ d -> e- ve~ a u d~

      sign(1)=1
      sign(2)=1
      sign(3)=-1
      sign(4)=-1
      do mu=0,3
      pbar(mu,1)=p(mu,2)   ! 1-> ket 
      pbar(mu,2)=p(mu,3)   ! 2-> bra
      pbar(mu,3)=p(mu,4)   ! 3-> ket
      pbar(mu,4)=p(mu,1)   ! 4-> bra
      enddo

      END SELECT

      call psi0m(4,pbar(0,1),sign(1),psi)

      do mu=0,3
       do i=1,4
         qbar(mu,i)=sign(i)*pbar(mu,i)
       enddo
      enddo

      do mu=0,3
      Q2(mu)= -qbar(mu,4)+qbar(mu,3)!p(mu,2)-p(mu,4)
      Q2a(mu)=-qbar(mu,4)+qbar(mu,1)!p(mu,1)-p(mu,4)
      Q2uza(mu)=-qbar(mu,4)+qbar(mu,1)-photfourvec(mu)!p(mu,1)-p(mu,4)-v(mu,3)
      Q2dza(mu)=-qbar(mu,4)+qbar(mu,3)-photfourvec(mu)!p(mu,2)-p(mu,4)-v(mu,3)
      enddo
      propwa(1)=1d0/(dotrr(q2dza,q2dza))
      propwa(2)=1d0/(dotrr(q2uza,q2uza))
      propa(1)=1d0/dotrr(q2,q2)
      propa(2)=1d0/dotrr(q2a,q2a)

!
! get the f-fbar currents
!
      call curr6(1,psi(1,-1,2),qbar(0,2),psi(1,-1,1),qbar(0,1),jqqu(0,-1,1))
      call curr6(1,psi(1,-1,4),qbar(0,4),psi(1,-1,3),qbar(0,3),jqqu(0,-1,2))

      call curr6(1,psi(1,-1,2),qbar(0,4),psi(1,-1,3),qbar(0,1),jqqd(0,-1,2))
      call curr6(1,psi(1,-1,4),qbar(0,2),psi(1,-1,1),qbar(0,3),jqqd(0,-1,1))   

!------------------------------------------------------------
!------------------------------------------------------------
      do isig1=-1,1,2 ! sum over external photon helicity   +
!------------------------------------------------------------
!------------------------------------------------------------

      do i=-1,1,2
       do j=-1,1,2
        do k=-1,1,2
         do l=-1,1,2
!                   u u d u
         ubububdb(1,i,j,k,l)=dcmplx(0d0,0d0)
         ubububdb(2,i,j,k,l)=dcmplx(0d0,0d0)
         enddo
        enddo
       enddo
      enddo


      if(ldebug) then

      do mu=0,3
      helastestmom(mu,1)=p(mu,2)
      helastestmom(mu,2)=p(mu,1)
      helastestmom(mu,3)=v(mu,1)
      helastestmom(mu,4)=v(mu,2)
      helastestmom(mu,5)=v(mu,3)
      helastestmom(mu,6)=p(mu,3)
      helastestmom(mu,7)=p(mu,4)
      enddo

      do k=1,7
      ipart(k)=1
      enddo

      ihel(1)=-1
      ihel(2)=-1
      ihel(6)=-1     !QUARKS
      ihel(7)=-1
!------------------------
      ihel(3)=-1
      ihel(4)=1
      ihel(5)=isig1
!------------------------
      endif


       do chi=-1,1,2

!  photon radiation off the quark legs

        call bra2c(psi(1,chi,2),.true.,qbar(0,2),chi,photfourvec(0),wafpolc(0,isig1),
     &              psiwmaz(1,chi,2),fermpphot(0,2))
        call ket2c(psi(1,chi,1),.true.,qbar(0,1),chi,photfourvec(0),wafpolc(0,isig1),
     &              psiwmaz(1,chi,1),fermpphot(0,1))
        call bra2c(psi(1,chi,4),.true.,qbar(0,4),chi,photfourvec(0),wafpolc(0,isig1),
     &              psiwmaz(1,chi,4),fermpphot(0,4))
        call ket2c(psi(1,chi,3),.true.,qbar(0,3),chi,photfourvec(0),wafpolc(0,isig1),
     &              psiwmaz(1,chi,3),fermpphot(0,3))
        enddo

! form currents
!
        call curr6(1,psi(1,-1,4),qbar(0,4),psiwmaz(1,-1,1),fermpphot(0,1),jqqaz(0,-1,1))
        call curr6(1,psiwmaz(1,-1,4),fermpphot(0,2),psi(1,-1,1),qbar(0,3),jqqza(0,-1,1))

        call curr6(1,psi(1,-1,4),qbar(0,4),psiwmaz(1,-1,3),fermpphot(0,3),jqqaz(0,-1,2))
        call curr6(1,psiwmaz(1,-1,4),fermpphot(0,4),psi(1,-1,3),qbar(0,3),jqqza(0,-1,2))

! W radiation off the quark legs, subsequent leptonic decay

      call bra2c(psi(1,-1,4),.true.,qbar(0,4),-1,wemnumom(0),currlepw(0),
     &              psiww(1,-1,4),fermpphotww(0,4))       
      call bra2c(psi(1,-1,2),.true.,qbar(0,2),-1,wemnumom(0),currlepw(0),
     &              psiww(1,-1,2),fermpphotww(0,2))
      call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,wemnumom(0),currlepw(0),
     &              psiww(1,-1,1),fermpphotww(0,1))
      call ket2c(psi(1,-1,3),.true.,qbar(0,3),-1,wemnumom(0),currlepw(0),
     &              psiww(1,-1,3),fermpphotww(0,3))


      Call curr6(-1,psi(1,-1,2),qbar(0,2),psiww(1,-1,1),fermpphotww(0,1),jqqvw(0,-1,1))
      Call curr6(-1,psiww(1,-1,2),fermpphotww(0,2),psi(1,-1,1),qbar(0,1),jqqvw(0,-1,2))

      Call curr6(-1,psi(1,-1,2),qbar(0,2),psiww(1,-1,3),fermpphotww(0,3),jqqvw(0,-1,3))
      Call curr6(-1,psiww(1,-1,2),fermpphotww(0,2),psi(1,-1,3),qbar(0,3),jqqvw(0,-1,4))

      do isig2=-1,1,2
!
!  W A / W G  case
!
!  lower line GA
!
      temp1=dotcc(jqqvw(0,-1,1),jqqaz(0,isig2,2))
      temp2=dotcc(jqqvw(0,-1,1),jqqza(0,isig2,2))
      temp3=dotcc(jqqvw(0,-1,2),jqqaz(0,isig2,2))
      temp4=dotcc(jqqvw(0,-1,2),jqqza(0,isig2,2))

      temp=(temp1+temp2+temp3+temp4)*gg(1)**2*clr(1,3,-1)*clr(4,1,isig2)*propwa(1)

      ubububdb(1,-1,isig2,-1,isig2)=ubububdb(1,-1,isig2,-1,isig2)+temp

!       test2=(temp1+temp2+temp3+temp4)*gg(1)**2*clr(1,3,-1)*clr(4,1,isig2)*propwa(1)
!       print*, "MY Code", test2
!       CALL Stestdd(helastestmom,ihel,ipart,help1c)
!       print*, "HELAS  ", help1c
!       print*, "Frac   ", help1c/test2
! 
!       STOP


!  upper line GA

      temp1=dotcc(jqqvw(0,-1,3),jqqaz(0,isig2,1))
      temp2=dotcc(jqqvw(0,-1,3),jqqza(0,isig2,1))
      temp3=dotcc(jqqvw(0,-1,4),jqqaz(0,isig2,1))
      temp4=dotcc(jqqvw(0,-1,4),jqqza(0,isig2,1))

      temp=(temp1+temp2+temp3+temp4)*gg(1)**2*clr(1,3,-1)*clr(4,1,isig2)*propwa(2)

      ubububdb(2,isig2,-1,-1,isig2)=ubububdb(2,-1,isig2,-1,isig2)+temp

!       test2=ubububdb(2,-1,-1,-1,-1)
!       print*, "MY Code", test2
!       CALL Stestdd(helastestmom,ihel,ipart,help1c)
!       print*, "HELAS  ", help1c
!       print*, "Frac   ", help1c/test2
! 
!       STOP

      enddo

!  finally emit up to three gauge bosons from a single quark line

! upper case

      Call bra2c(psiww(1,-1,4),.false.,fermpphotww(0,4),-1,photfourvec(0),wafpolc(0,isig1),
     &           psiwwa(1,-1,1),fermpphotw(0,1))

      Call bra2c(psiww(1,-1,2),.false.,fermpphotww(0,2),-1,photfourvec(0),wafpolc(0,isig1),
     &           psiwwa(1,-1,3),fermpphotw(0,3))

      Call bra2c(psiwmaz(1,-1,2),.false.,fermpphot(0,2),-1,wemnumom(0),currlepw(0),
     &           psiwwa(1,-1,4),fermpphotw(0,4))

      call bra2c(psi(1,-1,2),.true.,qbar(0,2),-1,wemnugamom(0),currlepwga(0,isig1),
     &           psiwwa(1,-1,5),fermpphotw(0,5))

!      call curr6(-1,psiwwa(1,-1,5),fermpphotw(0,5),psi(1,-1,3),qbar(0,3),jqqwwa(0,-1,7))

      Call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,wemnumom(0),currlepw(0),
     &           psiwaa(1,-1,1),fermpphotwa(0,1))

      Call curr6(-1,psiwmaz(1,-1,2),fermpphot(0,2),psiwaa(1,-1,1),fermpphotwa(0,1),jqqwaa(0,-1,1))
      Call curr6(-1,psiwwa(1,-1,4),fermpphotw(0,4),psi(1,-1,1),qbar(0,1),jqqwaa(0,-1,2))
      Call curr6(-1,psiww(1,-1,2),fermpphotww(0,2),psiwmaz(1,-1,1),fermpphot(0,1),jqqwaa(0,-1,3))
      Call curr6(-1,psiwwa(1,-1,3),fermpphotw(0,3),psi(1,-1,1),qbar(0,1),jqqwaa(0,-1,4))

      Call ket2c(psi(1,-1,1),.true.,qbar(0,1),-1,wemnugamom(0),currlepwga(0,isig1),
     &           psiwaa(1,-1,2),fermpphotwa(0,2))

      Call curr6(-1,psi(1,-1,2),qbar(0,2),psiwaa(1,-1,2),fermpphotwa(0,2),jqqwaa(0,-1,5))
      Call curr6(-1,psiwwa(1,-1,5),fermpphotw(0,5),psi(1,-1,1),qbar(0,1),jqqwaa(0,-1,6))

      do mu=0,3
      jqqwaaoff(mu,-1,1)=jqqwaa(mu,-1,5)+jqqu(mu,-1,1)*
     &              dotrc(wemnugamom,currlepwga(0,isig1))/dcmplx(xm2(3),-xmg(3))
      jqqwaaoff(mu,-1,2)=jqqwaa(mu,-1,6)-jqqu(mu,-1,1)*
     &              dotrc(wemnugamom,currlepwga(0,isig1))/dcmplx(xm2(3),-xmg(3))
      enddo

      Call ket2c(psiwmaz(1,-1,1),.false.,fermpphot(0,1),-1,wemnumom(0),currlepw(0),
     &           psiwaa(1,-1,3),fermpphotwa(0,3))
      Call curr6(-1,psi(1,-1,2),qbar(0,2),psiwaa(1,-1,3),fermpphotwa(0,3),jqqwaa(0,-1,7))

      Call ket2c(psiwaa(1,-1,1),.false.,fermpphotwa(0,1),-1,photfourvec(0),wafpolc(0,isig1),
     &           psiwaa(1,-1,4),fermpphotwa(0,4))
      Call curr6(-1,psi(1,-1,2),qbar(0,2),psiwaa(1,-1,4),fermpphotwa(0,4),jqqwaa(0,-1,8))

      do isig2=-1,1,2
      temp1= dotcc(jqqwaa(0,-1,1),jqqu(0,isig2,2))
      temp2= dotcc(jqqwaa(0,-1,2),jqqu(0,isig2,2))
      temp3= dotcc(jqqwaa(0,-1,3),jqqu(0,isig2,2))
      temp4= dotcc(jqqwaa(0,-1,4),jqqu(0,isig2,2))
      temp5= -dotcc(jqqwaaoff(0,-1,1),jqqu(0,isig2,2)) 
      temp6= -dotcc(jqqwaaoff(0,-1,2),jqqu(0,isig2,2))
      temp7= dotcc(jqqwaa(0,-1,7),jqqu(0,isig2,2))
      temp8= dotcc(jqqwaa(0,-1,8),jqqu(0,isig2,2))

      ubububdb(1,-1,isig2,-1,isig2)= ubububdb(1,-1,isig2,-1,isig2)+ 
     &  ( temp1*(clr(3,1,-1))
     &   +temp2*(clr(3,1,-1))
     &   +temp3*(clr(4,1,-1))
     &   +temp4*(clr(4,1,-1))
     &   +temp5
     &   +temp6
     &   +temp7*(clr(4,1,-1))
     &   +temp8*(clr(3,1,-1))
     &     )*propa(1)*clr(1,3,-1)*gg(1)**2

!       print*, "MY Code",ubububdb(1,-1,-1,-1,-1)
!       CALL Stestdd(helastestmom,ihel,ipart,help1c)
!       print*, "HELAS  ", help1c
!       print*, "Frac   ", help1c/ubububdb(1,-1,-1,-1,-1)
! 
!       STOP
       enddo

!  lower case
!
      Call curr6(-1,psiwmaz(1,-1,2),fermpphot(0,2),psiww(1,-1,3),fermpphotww(0,3),jqqwaa(0,-1,9))
      Call curr6(-1,psiwwa(1,-1,3),fermpphotw(0,3),psi(1,-1,3),qbar(0,3),jqqwaa(0,-1,10))
      Call curr6(-1,psiww(1,-1,2),fermpphotww(0,2),psiwmaz(1,-1,3),fermpphot(0,3),jqqwaa(0,-1,11))

      call ket2c(psi(1,-1,3),.true.,qbar(0,3),-1,wemnugamom(0),currlepwga(0,isig1),
     &           psiwaa(1,-1,4),fermpphotw(0,4))
      Call curr6(-1,psi(1,-1,2),qbar(0,3),psiwaa(1,-1,4),fermpphotwa(0,4),jqqwaa(0,-1,12))


      call bra2c(psi(1,-1,2),.true.,qbar(0,2),-1,wemnugamom(0),currlepwga(0,isig1),
     &           psiwwa(1,-1,8),fermpphotw(0,8))

      Call curr6(-1,psiwwa(1,-1,8),fermpphotw(0,8),psi(1,-1,3),qbar(0,3),jqqwaa(0,-1,13))

      do mu=0,3
      jqqwaaoff(mu,-1,3)=jqqwaa(mu,-1,12)+jqqd(mu,-1,2)*
     &              dotrc(wemnugamom,currlepwga(0,isig1))/dcmplx(xm2(3),-xmg(3))
      jqqwaaoff(mu,-1,4)=jqqwaa(mu,-1,13)-jqqd(mu,-1,2)*
     &              dotrc(wemnugamom,currlepwga(0,isig1))/dcmplx(xm2(3),-xmg(3))
      enddo

      Call bra2c(psiwmaz(1,-1,2),.false.,fermpphot(0,2),-1,wemnumom(0),currlepw(0),
     &           psiwaa(1,-1,5),fermpphotwa(0,5))

      Call curr6(-1,psiwaa(1,-1,5),fermpphotwa(0,5),psi(1,-1,3),qbar(0,3),jqqwaa(0,-1,14))

      Call ket2c(psiwmaz(1,-1,3),.false.,fermpphot(0,3),-1,wemnumom(0),currlepw(0),
     &           psiwaa(1,-1,3),fermpphotwa(0,3))
      Call ket2c(psiww(1,-1,3),.false.,fermpphotww(0,3),-1,photfourvec,wafpolc(0,isig1),
     &           psiwaa(1,-1,4),fermpphotwa(0,4))

      Call curr6(-1,psi(1,-1,2),qbar(0,2),psiwaa(1,-1,3),fermpphotwa(0,3),jqqwaa(0,-1,15))
      Call curr6(-1,psi(1,-1,2),qbar(0,4),psiwaa(1,-1,4),fermpphotwa(0,4),jqqwaa(0,-1,16))

      do isig2=-1,1,2
      temp1= dotcc(jqqwaa(0,-1,9),jqqd(0,isig2,1))
      temp2= dotcc(jqqwaa(0,-1,10),jqqd(0,isig2,1))
      temp3= dotcc(jqqwaa(0,-1,11),jqqd(0,isig2,1))
      temp4= dotcc(jqqwaa(0,-1,16),jqqd(0,isig2,1))
      temp5= -dotcc(jqqwaaoff(0,-1,3),jqqd(0,isig2,1))
      temp6= -dotcc(jqqwaaoff(0,-1,4),jqqd(0,isig2,1))
      temp7= +dotcc(jqqwaa(0,-1,14),jqqd(0,isig2,1))
      temp8= +dotcc(jqqwaa(0,-1,15),jqqd(0,isig2,1))

      ubububdb(2,isig2,-1,-1,isig2)= ubububdb(2,isig2,-1,-1,isig2)+
     & ( temp1*(clr(3,1,-1))
     &  +temp2*(clr(4,1,-1))
     &  +temp3*(clr(4,1,-1))
     &  +temp4*(clr(3,1,-1))
     &  +temp5
     &  +temp6
     &  +temp7*(clr(3,1,-1))
     &  +temp8*(clr(4,1,-1))
     &   )*clr(1,3,-1)*propa(2)*gg(1)**2

!       print*, "MY Code", ubububdb(2,-1,-1,-1,-1)
!       CALL Stestdd(helastestmom,ihel,ipart,help1c)
!       print*, "HELAS  ", help1c
!       print*, "Frac   ", help1c/ubububdb(2,-1,-1,-1,-1)
! 
!       STOP

      enddo

! MATRIXELEMENT

      do i=-1,1,2
       do j=-1,1,2
        do k=-1,1,2
         do l=-1,1,2
           m2s_qqwgaglujdd(1)=m2s_qqwgaglujdd(1)
     &              +2d0*abs(ubububdb(1,i,j,k,l))**2
     &              +2d0*abs(ubububdb(2,i,j,k,l))**2
     &              +4d0/3d0*DBLE(ubububdb(1,i,j,k,l)*dconjg(ubububdb(2,i,j,k,l)))
           m2s_qqwgaglujdd(2)=m2s_qqwgaglujdd(2)
     &              +2d0*abs(ubububdb(2,i,j,k,l))**2
           m2s_qqwgaglujdd(3)=m2s_qqwgaglujdd(3)
     &              +2d0*abs(ubububdb(1,i,j,k,l))**2
         enddo
        enddo
       enddo
      enddo

!----------------------------------------------------------------------#+
!CHECK
!----------------------------------------------------------------------#+
! 
!       if (selid.eq.1) then
! 
!       do i=-1,1,2
!        do j=-1,1,2
!         do k=-1,1,2
!          do l=-1,1,2
! 
!       ihel(1)=-i
!       ihel(2)=-j
!       ihel(7)=-k
!       ihel(6)=-l
!       CALL STESTREALUBD(helastestmom,ihel,ipart,help1c)!C2
!       test1=
!      &              +2d0*abs(ubububdb(1,i,j,k,l))**2
!      &              +2d0*abs(ubububdb(2,i,j,k,l))**2
!      &              +2d0/3d0*2d0*DBLE(ubububdb(1,i,j,k,l)*dconjg(ubububdb(2,i,j,k,l)))
! 
! 
!       PRINT*, "My Code ", test1
!       Print*, i,j,k,l
!       PRINT*, "HELAS   ", abs(help1c)
!       print*, ihel(1),ihel(2),ihel(3),ihel(4)
!       PRINT*, "Fraction", abs(help1c)/test1
!       Print*, "-------------------------"
!          enddo
!         enddo
!        enddo
!       enddo
! 
!       elseif (selid.eq.2) then
! 
!       do i=-1,1,2
!        do j=-1,1,2
!         do k=-1,1,2
!          do l=-1,1,2
!            m2s_qqwgaglujdd=m2s_qqwgaglujdd
!      &              +2d0*abs(ubububdb(2,i,j,k,l))**2
! 
!          enddo
!         enddo
!        enddo
!       enddo
! 
!       do i=-1,1,2
!        do j=-1,1,2
!         do k=-1,1,2
!          do l=-1,1,2
! 
!       ihel(1)=-i
!       ihel(2)=-j
!       ihel(7)=-k
!       ihel(6)=-l
!       CALL Subsbdbsb(helastestmom,ihel,ipart,help1)!C2
!       test1=   +2d0*abs(ubububdb(1,i,j,k,l))**2
! 
! 
!       PRINT*, "My Code ", test1
!       Print*, i,j,k,l
!       PRINT*, "HELAS   ", help1
!       print*, ihel(1),ihel(2),ihel(3),ihel(4)
!       PRINT*, "Fraction", help1/test1
!       Print*, "-------------------------"
!          enddo
!         enddo
!        enddo
!       enddo
!       STOP
! 
!       else
!  
!       do i=-1,1,2
!        do j=-1,1,2
!         do k=-1,1,2
!          do l=-1,1,2
!            m2s_qqwgaglujdd=m2s_qqwgaglujdd
!      &              +2d0*abs(ubububdb(1,i,j,k,l))**2
! 
!          enddo
!         enddo
!        enddo
!       enddo
! 
!       do i=-1,1,2
!        do j=-1,1,2
!         do k=-1,1,2
!          do l=-1,1,2
! 
!       ihel(1)=-i
!       ihel(2)=-j
!       ihel(7)=-k
!       ihel(6)=-l
!       CALL SDSUS(helastestmom,ihel,ipart,help1c)!C2
!       test1=  +2d0*abs(ubububdb(1,i,j,k,l))**2
! 
!       PRINT*, "My Code ", test1
!       Print*, i,j,k,l
!       PRINT*, "HELAS   ", abs(help1c)
!       print*, ihel(1),ihel(2),ihel(3),ihel(4)
!       PRINT*, "Fraction", abs(help1c)/test1
!       Print*, "-------------------------"
!          enddo
!         enddo
!        enddo
!       enddo
! 
! 
!       endif
! 
! 
!       STOP

      enddo ! close loop over external photon helicity

      end

! THIS SUBROUTINE CALCULATES THE DIPOLES ACCORDING TO CATANI-SEYMOUR

      subroutine daisjdipole67(sel,kin,xuz,pi,pk,pa,p,v,dipoleres)

      implicit none

#include "VBFNLO/utilities/coupl.inc"

      integer mu,nu,sel,kin
c INPUT

c
      double precision pi(0:3),pk(0:3),pa(0:3)
      double precision pis,CF,CA
      parameter (pis=3.141592653589793238d0,CF=4d0/3d0,CA=3d0)
      logical ldebug
      parameter (ldebug=.false.)

c OUTPUT
c
      double precision dipole(0:3,0:3),p(0:3,3),v(0:3,3)
     &    ,pipk,pkpa,pipa,xuz(2),dipoleres(3),bornmom(0:3,6),color,symmetry
c
      double complex borncurr(0:3,-1:1)
      double precision Contract_Trjcj,borntest1,borntest2,metric(0:3,0:3)

      double precision xika,ui
      double precision dotrr
      external dotrr,Contract_Trjcj

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

      pipk=dotrr(pi,pk)
      pkpa=dotrr(pk,pa)
      pipa=dotrr(pi,pa)

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
      if (abs(p(mu,3)-(pk(mu)+pi(mu)-(1d0-xika)*pa(mu))).ge.1e-10) then
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
      bornmom(mu,4)=v(mu,3)
      bornmom(mu,5)=v(mu,1)
      bornmom(mu,6)=v(mu,2)
      enddo

      if (sel.eq.1) then
      Call UGDAEVEM(bornmom,borncurr)
      color=4d0      !born matrix element colorfactor =CF*Trace(1_3)
     &     *(-1d0)/2d0 !Dipole colorfactor = -1/2 CA/CA 
      symmetry=1d0/96d0  ! 1/3*1/8 (inital state colors) * 1/2 *1/2 (initial state helicities)

      if (ldebug) then
      CAll SUGDAEVE(bornmom,borntest1) !u~ g > d a ve
      borntest2=-((contract_Trjcj(metric,borncurr(0,-1))+contract_Trjcj(metric,borncurr(0,1)))*4/96)
      if (abs(borntest1-borntest2).ge.1e-10) then
      Print*, "error in matrixelement calculation 67:", (BORNTEST1-BORNTEST2)/BORNTEST1*100d0, "%"
      STOP
      endif
      endif

      elseif (sel.eq.2) then
      Call DGUAEVEM(bornmom,borncurr)
      color=4d0      !born matrix element colorfactor =CF*Trace(1_3)
     &     *(-1d0)/2d0 !Dipole colorfactor = -1/2 CA/CA 
      symmetry=1d0/96d0  ! 1/3*1/8 (inital state colors) * 1/2 *1/2 (initial state helicities)

      if (ldebug) then      
      CAll SDGUAEVE(bornmom,borntest1)
      borntest2=-1d0*((contract_Trjcj(metric,borncurr(0,-1))+contract_Trjcj(metric,borncurr(0,1)))*4./96.)
      if (abs(borntest1-borntest2).ge.1e-10) then
      Print*, "error in matrixelement calculation 67:", (BORNTEST1-BORNTEST2)/BORNTEST1*100d0, "%"
      STOP
      endif
      endif

      endif

      do mu=0,3
      do nu=0,mu
      dipole(mu,nu)=-1d0/2d0/pipa/xika*8d0*pis*alfas*CF*
     &  ( (1d0-xika)/xika*2d0*ui*(1d0-ui)/pipk*(pi(mu)/ui-pk(mu)/(1d0-ui))*(pi(nu)/ui-pk(nu)/(1d0-ui)))
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

      dipoleres(1)=(Contract_Trjcj(dipole(0,0),borncurr(0,-1))+Contract_Trjcj(dipole(0,0),borncurr(0,1)))
     &           *color*symmetry

      do mu=0,3
      do nu=0,mu
      dipole(mu,nu)=-1d0/2d0/pipa/xika*16d0*pis*alfas*CA*
     &  ( (1d0-xika)/xika*ui*(1d0-ui)/pipk*(pi(mu)/ui-pk(mu)/(1d0-ui))*(pi(nu)/ui-pk(nu)/(1d0-ui)))
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

      dipoleres(2)=(Contract_Trjcj(dipole(0,0),borncurr(0,-1))+Contract_Trjcj(dipole(0,0),borncurr(0,1)))
     &           *color*symmetry


      dipole(0,0)=-1d0/2d0/pipa/xika*8.*pis*alfas*CF*(2./(1.-xika+ui)-(1.+xika))

      color=(CA-2.*CF)/2./CF*4.
      dipoleres(3)=-(Contract_Trjcj(metric(0,0),borncurr(0,-1))+Contract_Trjcj(metric(0,0),borncurr(0,1)))
     &           *color*symmetry*dipole(0,0)


      end


      subroutine daibdipole147(sel,kin,xuz,pi,pa,pb,p,v,dipoleres)

      implicit none

#include "VBFNLO/utilities/coupl.inc"

      integer mu,nu,sel,kin
c INPUT
c
      double precision pi(0:3),pb(0:3),pa(0:3),p(0:3,3),v(0:3,3)
      double precision pis,CF,CA
      parameter (pis=3.141592653589793238d0,CF=4d0/3d0,CA=3d0)
      logical ldebug
      parameter (ldebug=.false.)
c OUTPUT
c
      double precision dipole(0:3,0:3),xuz(2),pges(0:3),
     & pipb,papb,pipa,bornmom(0:3,6),dipoleres(3),color,symmetry,
     & borncache147(2,2),borncache145(3,2)

      common/BORNCACHE/borncache145,borncache147 

      double complex borncurr(0:3,-1:1)
      double precision Contract_Trjcj,borntest1,borntest2,metric(0:3,0:3)

c
      double precision xiab,ui
      double precision dotrr
      external dotrr,Contract_Trjcj

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

      pipb=dotrr(pi,pb)
      papb=dotrr(pa,pb)
      pipa=dotrr(pi,pa)

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
      bornmom(mu,4)=v(mu,3)   !a
      bornmom(mu,5)=v(mu,1)   !e-
      bornmom(mu,6)=v(mu,2)   !ve~
      enddo

      if (sel.eq.1) then
      Call UGDAEVEM(bornmom,borncurr)
      color=4d0      !born matrix element colorfactor =CF*Trace(1_3)
     &     *(-1d0)/2d0 !Dipole colorfactor = -1/2 CA/CA 
      symmetry=1d0/96d0  ! 1/2*1/8 (inital state colors) * 1/2 *1/2 (initial state helicities)
     
      if (ldebug) then
      CAll SUGDAEVE(bornmom,borntest1)
      borntest2=-((contract_Trjcj(metric,borncurr(0,-1))+contract_Trjcj(metric,borncurr(0,1)))*4./96.)
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
      Call DGUAEVEM(bornmom,borncurr)
      color=4d0      !born matrix element colorfactor =CF*Trace(1_3)
     &     *(-1d0)/2d0 !Dipole colorfactor = -1/2 CA/CA 
      symmetry=1d0/96d0  ! 1/2*1/8 (inital state colors) * 1/2 *1/2 (initial state helicities)

    
      if (ldebug) then
      CAll SDGUAEVE(bornmom,borntest1)
      borntest2=-((contract_Trjcj(metric,borncurr(0,-1))+contract_Trjcj(metric,borncurr(0,1)))*4./96.)
      if (abs(borntest1-borntest2).ge.1e-10) then
      Print*, "error in matrixelement calculation 147:", BORNTEST1, BORNTEST2
      STOP
      endif
      endif

      endif

      do mu=0,3
      do nu=0,mu
      dipole(mu,nu)=-1d0/2d0/pipa/xiab*8d0*pis*alfas*CF*
     &  ( (1d0-xiab)/xiab*2d0*papb/(pipa*pipb)*(pi(mu)-pipa/papb*pb(mu))*(pi(nu)-pipa/papb*pb(nu)))
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

      dipoleres(1)=(Contract_Trjcj(dipole(0,0),borncurr(0,-1))+Contract_Trjcj(dipole(0,0),borncurr(0,1)))
     &           *color*symmetry


      do mu=0,3
      do nu=0,mu
      dipole(mu,nu)=-1d0/2d0/pipa/xiab*16d0*pis*alfas*CA*
     &  ( (1d0-xiab)/xiab*papb/pipa/pipb*(pi(mu)-pipa/papb*pb(mu))*(pi(nu)-pipa/papb*pb(nu)))
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

      dipoleres(2)=(Contract_Trjcj(dipole(0,0),borncurr(0,-1))+Contract_Trjcj(dipole(0,0),borncurr(0,1)))
     &           *color*symmetry


      dipole(0,0)=-1d0/2d0/pipa/xiab*8.*pis*alfas*CF*(2./(1.-xiab)-(1.+xiab))  
      color=-CA/2./CF*4.
      dipoleres(3)=-(Contract_Trjcj(metric(0,0),borncurr(0,-1))+Contract_Trjcj(metric(0,0),borncurr(0,1)))
     &           *color*symmetry*dipole(0,0)

      borncache147(sel,kin)=-((contract_Trjcj(metric,borncurr(0,-1))+contract_Trjcj(metric,borncurr(0,1)))*4./96.)

!       print*, "147 ",sel,kin,-(Contract_Trjcj(metric(0,0),borncurr(0,-1))+Contract_Trjcj(metric(0,0),borncurr(0,1)))

      end


      subroutine daisjdipole65(sel,kin,xuz,pi,pk,pa,p,v,dipoleres)

      implicit none

#include "VBFNLO/utilities/coupl.inc"

      integer mu,sel,kin
c INPUT

      double precision pi(0:3),pk(0:3),pa(0:3)
      double precision pis,CF,CA,Tr
      parameter (pis=3.141592653589793238d0,CF=4d0/3d0,CA=3d0,Tr=1d0/2d0)
      logical ldebug
      parameter (ldebug=.false.)

c OUTPUT
c
      double precision dipole,p(0:3,3),v(0:3,3)
     &    ,pipk,pkpa,pipa,xuz(2),dipoleres(6),bornmom(0:3,6),color(3)
c
      double precision born
      double precision Contract_Trjcj

      double precision xika,ui
      double precision dotrr
      external dotrr,Contract_Trjcj

      pipk=dotrr(pi,pk)
      pkpa=dotrr(pk,pa)
      pipa=dotrr(pi,pa)

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
      if (abs(p(mu,3)-(pk(mu)+pi(mu)-(1d0-xika)*pa(mu))).ge.1e-10) then
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
      bornmom(mu,4)=v(mu,3)  ! a
      bornmom(mu,5)=v(mu,1)  ! e-
      bornmom(mu,6)=v(mu,2)  ! ve~
      enddo

      color(1)=-CA/2./CF
      color(2)=(CA-2.*CF)/2./CF
      color(3)=0.

      if (sel.eq.1) then
      Call SUDGAENU(bornmom,born)
      elseif (sel.eq.2) then
      Call SUGDAEVE(bornmom,born)
      elseif (sel.eq.3) then
      Call SDGUAEVE(bornmom,born)
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


      subroutine daibdipole145(sel,kin,xuz,pi,pa,pb,p,v,dipoleres)

      implicit none

#include "VBFNLO/utilities/coupl.inc"

      integer mu,sel,kin
c INPUT
c
      double precision pi(0:3),pb(0:3),pa(0:3),p(0:3,3),v(0:3,3)
      double precision pis,CF,CA,TR
      parameter (pis=3.141592653589793238d0,CF=4d0/3d0,CA=3d0,TR=1d0/2d0)!
      logical ldebug
      parameter (ldebug=.false.)
c OUTPUT
c
      double precision dipole,xuz(2),
     &    pipb,papb,pipa,bornmom(0:3,6),dipoleres(6),color(3),
     &    borncache145(3,2),borncache147(2,2)

      common/BORNCACHE/borncache145,borncache147

      double precision born
      double precision Contract_Trjcj

c
      double precision xiab,ui
      double precision dotrr
      external dotrr,Contract_Trjcj

      pipb=dotrr(pi,pb)
      papb=dotrr(pa,pb)
      pipa=dotrr(pi,pa)

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
      bornmom(mu,4)=v(mu,3)   !a
      bornmom(mu,5)=v(mu,1)   !e-
      bornmom(mu,6)=v(mu,2)   !ve~
      enddo

      color(1)=(CA-2d0*CF)/2d0/CF
      color(2)=-CA/2d0/CF
      color(3)=0d0

      if (sel.eq.1) then
      Call SUDGAENU(bornmom,born)
      elseif (sel.eq.2) then
      Call SUGDAEVE(bornmom,born)
      elseif (sel.eq.3) then
      Call SDGUAEVE(bornmom,born)
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


      subroutine dasijdipole40(sel,kin,xuz,pi,pj,pa,p,v,dipoleres)

      implicit none

#include "VBFNLO/utilities/coupl.inc"

      integer mu,nu,sel,kin
c INPUT

c
      double precision pi(0:3),pj(0:3),pa(0:3)
      double precision pis,CF,CA,TR
      parameter (pis=3.141592653589793238d0,CF=4d0/3d0,CA=3d0,TR=1d0/2d0)
      logical ldebug
      parameter (ldebug=.false.)

c OUTPUT
c
      double precision dipole(0:3,0:3),p(0:3,3),v(0:3,3)
     &    ,pipj,pjpa,pipa,xuz(2),dipoleres(3),bornmom(0:3,6),color,symmetry,borntmp,
     &    dipolesc
c
      double complex borncurr(0:3,-1:1)
      double precision Contract_Trjcj

      double precision xija,zi,zj
      double precision dotrr
      external dotrr,Contract_Trjcj

      borntmp=0d0
      dipolesc=0d0

      pipj=dotrr(pi,pj)
      pjpa=dotrr(pj,pa)
      pipa=dotrr(pi,pa)

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
      if (abs(p(mu,3)-(pj(mu)+pi(mu)-(1d0-xija)*pa(mu))).ge.1e-10) then
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
      bornmom(mu,4)=v(mu,3)
      bornmom(mu,5)=v(mu,1)
      bornmom(mu,6)=v(mu,2)
      enddo

      if (kin.eq.1.or.kin.eq.4) then
      Call UDGAENUM(bornmom,borncurr)
      color=4d0      !born matrix element colorfactor =CF*Trace(1_3)
     &     *(-1d0)/2d0 !Dipole colorfactor = -1/2 CA/CA 
      symmetry=1d0/36d0  ! 1/3*1/3 (inital state colors) * 1/2 *1/2 (initial state helicities)

      elseif (kin.eq.2) then
      CALL SUGDAEVE(bornmom,borntmp)

      elseif (kin.eq.3) then
      CALL SDGUAEVE(bornmom,borntmp)
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
!       CAll SUDGAENU(bornmom,borntest1)
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
     &  1d0/pipj*(zi*pi(mu)-zj*pj(mu))*(zi*pi(nu)-zj*pj(nu))
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


      dipoleres(1)=(Contract_Trjcj(dipole(0,0),borncurr(0,-1))+Contract_Trjcj(dipole(0,0),borncurr(0,1)))
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
     &  (-2d0)/pipj*(zi*pi(mu)-zj*pj(mu))*(zi*pi(nu)-zj*pj(nu))
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

      dipoleres(1)=(Contract_Trjcj(dipole(0,0),borncurr(0,-1))+Contract_Trjcj(dipole(0,0),borncurr(0,1)))
     &           *color*symmetry

      dipoleres(2)=0d0
      dipoleres(3)=0d0

      endif

      end


      subroutine test40(sel,kin,xuz,pj,pi,pa,p,v,dipoleres)

      implicit none

#include "VBFNLO/utilities/coupl.inc"

      integer mu,nu,sel,kin
c INPUT

c
      double precision pi(0:3),pj(0:3),pa(0:3)
      double precision pis,CF,CA,TR
      parameter (pis=3.141592653589793238d0,CF=4d0/3d0,CA=3d0,TR=1d0/2d0)
      logical ldebug
      parameter (ldebug=.false.)

c OUTPUT
c
      double precision p(0:3,3),v(0:3,3)
     &    ,pipj,pjpa,pipa,xuz(2),dipoleres(3),bornmom(0:3,6),color,borntmp,
     &    dipolesc
c
      double complex borncurr(0:3,-1:1)
      double precision Contract_Trjcj,borntest2,metric(0:3,0:3),dptmp

      double precision xija,zi,ui,zj
      double precision dotrr
      external dotrr,Contract_Trjcj

      borntmp=0d0
      dipolesc=0d0

      pipj=dotrr(pi,pj)
      pjpa=dotrr(pj,pa)
      pipa=dotrr(pi,pa)

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
      bornmom(mu,4)=v(mu,3)
      bornmom(mu,5)=v(mu,1)
      bornmom(mu,6)=v(mu,2)
      enddo

      dptmp=-1./2./xija*8.*pis*alfas*4./3.*(2./(1.-zi+1.-xija)-1.-zi)!/pipj


      if (kin.eq.2) then
      CALL SUGDAEVE(bornmom,borntmp)
      elseif (kin.eq.3) then
      CALL SDGUAEVE(bornmom,borntmp)
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
      CAll UGDAEVEM(bornmom,borncurr)
      borntest2=-((contract_Trjcj(metric,borncurr(0,-1))+contract_Trjcj(metric,borncurr(0,1)))*4./96.)


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



