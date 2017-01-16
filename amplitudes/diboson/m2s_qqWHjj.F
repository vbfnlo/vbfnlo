! This subroutine returns the real emission contribution to the order alpha_s 
! hadronic matrix element squared for the process p p --> e- ve H j j
! jet. It is also used as part of the NLO calculation.

      SUBROUTINE m2s_qqWHjj(bos,nlo,lokt,xi,p,v,rn,xuz,NW,NWH,m2s)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/lha.inc"
#include "tensor.inc"


      ! Narrowwidth flags: 0 off (full width), 1 on (force onshell)
      ! NW: switches final W boson
      ! NWH: switches Higgs boson, default in amplitudes.F is 0
      integer NW,NWH
 
      integer i,j,k,mu,
     &  nlo,signpdf1,signpdf2,
     &   bos
      logical lokt(11),fincolcalc,pdfchange

      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)

      double precision xi(nx), p(0:3,max_p,max_kin), 
     &   v(0:3,max_v,max_kin), m2s_fey,
     &   m2s_f(4,3),pch(0:3,max_p,max_kin),
     &   muu(3),mdd(3),muuch(3),mddch(3),m2s_ch(4,3),
     &   q_sf,pdf(-6:6,2),pdfckm(-6:6,2),
     &   amp,m2s_dd(4,3),dotrr,
     &   m2s_ddch(4,3),
     &   mgg(4),mggch(4),m2s(0:11),
     &   fincol(12),fincolampsq(12),coll3,coll4,
     &   fincolconv,coll7,coll8,pdftmp(-6:6,2),
     &   rn(1)


      double precision dipole(10),
     &   xuz(2,2:11),kinem(11),dipolev(6,6),dipolea(3),dipoleb(3),borncache147(2,2),
     &   borncache(9,4),borncache145(3,2)
   
      ! stores individual matrix elements
      ! size should be determined by running and outputting maximum k
      double precision res(114)

      ! variables for family/flavor loop
      integer fam1, fl2, fam2, u, d, q, fam2end

      double precision zero,inv2,four,two,zero66(6,6)
      parameter (zero=0d0,inv2=1d0/2d0,four=4d0,two=2d0)

      common/BORNCACHE/borncache145,borncache147

      double complex Contract_Trjj
      external dotrr,Contract_Trjj

      real*8 partons(0:7,max_jets,max_kin)
      integer npartons(max_kin)
      common /partons/partons,npartons

      double precision maxeventweight
      common/maxevents/maxeventweight

      parameter(fincolcalc=.true.)

      integer FL_VVgg
      external FL_VVgg

c variables for leshouches interface
      DOUBLE PRECISION weight,rnumb,RandomNumber

      ! for inclusion of ckm effects
      double precision getpdf, getpdfckm, getpdfckm1, getpdfckm2

      ! switch between WHjj and Wjj
      logical hashiggs

      ! higgs propagator factor
      double precision fac
      ! constants required for higgs propagator factor are already included at BKOPOU

      real*8 rsep
      external rsep

      INTEGER init/0/
      SAVE init

      logical doSubprocesses, doNLOCalc  
      ! true if sub_number is used to do different subprocesses
      
      logical ldebug
      parameter (ldebug=.false.)
      logical lprintmomenta
      parameter (lprintmomenta=.false.)

c  helicity selection
c     INTEGER h
c     COMMON /hcount / h

      doNLOCalc = nlo.ne.0

      if ( init .eq. 0 ) then
         maxeventweight=0d0

         write(6,*) " "
         if (bos.eq.46) then
            pdfchange = .false.
            write(6,*) "W-Hjj amplitude square information:"
            hashiggs = .true.
         elseif (bos.eq.36) then
            pdfchange = .true.
            write(6,*) "W+Hjj amplitude square information:"
            hashiggs = .true.
         elseif (bos.eq.4) then
            pdfchange = .false.
            write(6,*) "W-jj amplitude square information:"
            hashiggs = .false.
         elseif (bos.eq.3) then
            pdfchange = .true.
            write(6,*) "W+jj amplitude square information:"
            hashiggs = .false.
         else 
            write(6,*) "m2s_qqWHjj.F called with invalid process."
            write(6,*) 'bos=',bos
         endif
         write(6,*) "---------------------------------------------"
         write(6,*) " "

         if (ldebug) then
            if (fincolcalc) then
               write(6,*) "Calculation of finite collinear terms included"
            else
!              write(6,*) "Calculation of finite collinear terms excluded"
               print*, "WARNING: finite collinear terms excluded!"
            endif
         endif

         call printnfl(.false.) 
         init = 1

      endif


! initialize
      amp=zero
      do i=1,6
         do j=1,6
            zero66(i,j)=zero
         enddo
      enddo
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
      do i=1,12
         fincol(i)=zero
      enddo
      m2s_fey=zero
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

! exchange incoming particles to calculate permutations simultaneously
      do mu=0,3
         pch(mu,1,1)=p(mu,2,1)
         pch(mu,2,1)=p(mu,1,1)
         pch(mu,3,1)=p(mu,3,1)
         pch(mu,4,1)=p(mu,4,1)
      enddo

      if (ldebug .and. .false.) then
         do i=1,1
            print*, "mom i"
            print*, "n_v", n_v
            do k=1,n_v
               print*, v(0:3,k,i)
            enddo
            print*, "n_p", n_p
            do k=1,n_p
               print*, p(0:3,k,i)
            enddo
         enddo
      endif

      doSubprocesses =  doNLOCalc .and. (Loops_sub_NLO.gt.1) 
     &       .or. (.not. doNLOCalc .and. (Loops_sub_LO.gt.1))
      ! sub_number splits different subprocesses
      ! sub_number 1: qqqq 
      ! sub_number 2: qqgg
      ! for details see implementation below 

!  call topologies for crossed kinematics

      call getyourscalesready_dibos(xi,pdfchange,pdf,pdfckm,1)

c for WH precalculate leptonic tensors
c lepton spinors and polarisation vectors

c select helicity: h = random number for lepton & photon helicity (h=1:4)
c     h = mod(h,4) + 1
c     ie = (-1)**(h+1)
c     iu = (-1)**((h-1)/2)

      do i=11,1,-1
         if (lokt(i)) then
! effective W currents
            do mu=0,3
               !TODO: are there common electroweak parts in all kinematics?
               wnumom(mu)=v(mu,1,i) !TODO: remove unused variables here, NLO case
               wemmom(mu)=v(mu,2,i)
               wemnumom(mu)=v(mu,1,i)+v(mu,2,i)
               qw0(mu,i)=wemnumom(mu)
               if (hashiggs) then
                  higmom(mu) = 0d0
                  do j=3,n_v ! sum all higgs decay products
                     higmom(mu) = higmom(mu) + v(mu,j,i)
                  enddo
                  wemnuhmom(mu)=v(mu,1,i)+v(mu,2,i)+higmom(mu)
                  qh(mu,i)=higmom(mu)
               else
                  qh(mu,i)=0d0 ! needed for dipoles
                  wemnuhmom(mu)=v(mu,1,i)+v(mu,2,i)
               endif
            enddo
            if (hashiggs) then
               higmom(4)=dotrr(higmom,higmom)
               wemnumom(4)=dotrr(wemnumom,wemnumom)
            endif
            wemnuhmom(4)=dotrr(wemnuhmom,wemnuhmom)


            if(ldebug .and. i.eq.1 .and. .false.) then
               print*, 'momnames'
               print*, 'higmom', higmom
               print*, 'wemnuhmom', wemnuhmom
               print*, 'qw0', qw0(:,i)
               print*, 'qh', qh(:,i)
            endif


            if (bos.eq.3 .or. bos.eq.36) then
               CALL OXXXXX(v(0,1,i),ZERO ,-1,1,wve(1,i))           !ve 
               CALL IXXXXX(v(0,2,i),ZERO ,1,-1,wep(1,i))           !e+ 
            elseif (bos.eq.4 .or. bos.eq.46) then
               CALL IXXXXX(v(0,1,i),ZERO ,1,-1,wep(1,i))           !e-  -> e+
               CALL OXXXXX(v(0,2,i),ZERO ,-1,1,wve(1,i))           !ve~ -> ve
            endif

c          CALL IXXXXX(v(0,1,i),ZERO ,1,-1,wep(1,i))           !e-  -> e+
c          CALL OXXXXX(v(0,2,i),ZERO ,-1,1,wve(1,i))           !ve~ -> ve

            if (hashiggs) then
               CALL SXXXXX(qh(0,i),1,hig(1,i))      
            endif


            if (hashiggs) then
               if (with_anom) then ! anomalous gauge boson couplings
                  if (bos.eq.36) then
                     call jwpmtowh_an(v(0,1,i),i,wtowh(0,i),3)  !W+ -> ve e+ h
                  elseif (bos.eq.46) then
                     call jwpmtowh_an(v(0,1,i),i,wtowh(0,i),4)  !W- -> ve e- h
                  endif
               else !SM
                  if (bos.eq.36) then
                     call jwpmtowh(v(0,1,i),i,wtowh(0,i),3)
                  else if (bos.eq.46) then
                     call jwpmtowh(v(0,1,i),i,wtowh(0,i),4)
                  endif
               endif !with_anom
            else ! no hashiggs
               CALL JIOXXX(wep(1,i),wve(1,i),GWF,WMASS,WWIDTH,wtowh(0,i))!W+- 
            endif

         endif !lokt
      enddo ! loop over kinematics


      ! todo: cache this somewhere (qh(4,i))
      ! masses independent of kinemativ
      do mu=0,3
         higmom(mu) = 0d0
         do j=3,n_v ! sum all higgs decay products
            higmom(mu) = higmom(mu) + v(mu,j,1)
         enddo
      enddo
      higmom(4)=dotrr(higmom,higmom)
! from amplitudes/vvjj/qqhqq.F      
c
!      if(ldebug) print*, 'hmass: ', higmom(4)

      if(hashiggs) then
         fac = 1.d0 / ( (higmom(4)-xm2(6))**2 + xmg(6)**2 ) !propagator^2 only
         fac = fac *16*pi*xmg(6)
      else
         fac = 1d0
      endif
c
c multiply by 16 pi*M Gamma to compensate for the phase space integral in
c LIPSN 
! (here LISVHj)
c
c   Int dq^2/2pi * dOmega/(32 pi^2) * 1/(q^2-M*2)^2+(MGamma)^2) 
c
c   = 1/(16 pi MGamma)
!   
! relative error due to finete integral-limits and approximation M**2 >> MGamma: 1e-5
c
c (where the propagator factor is fac above)


!
!  Born - Kinematics
!
      if (lokt(1)) then


      do i=1,4
         if(.not. doSubprocesses .or. (sub_number.eq.1) )  then
            call qqWHqq(p,v,i,muu,mdd)   ! Processes 8-11
            call qqWHqq(pch,v,i,muuch,mddch)  ! Processes 12-15
         endif
         do k=1,3
            m2s_dd(i,k)=mdd(k)
            m2s_ddch(i,k)=mddch(k)
            m2s_f(i,k)=muu(k)
            m2s_ch(i,k)=muuch(k)
         enddo

         if(.not. doSubprocesses .or. (sub_number.eq.2) )  then
            ! in qg->qg: p3 is gluon, p4 is quark
            call qqWHgg(p,v,i,mgg(i))  ! Processes 1-4
            if (i.ne.4) then
               call qqWHgg(pch,v,i,mggch(i))  !Processes 5-7
            endif
         endif

      enddo

       ! reset result counter
       k=FL_VVgg((/0,0,0,0/),-1)
       if(.not. doSubprocesses .or. (sub_number.eq.1) )  then
         !diagrams qqqq

         do fam1=1,2  !coupling to w
            do fl2=1,2 !observer
               if (fl2.eq.1 .and. (lincludebb .or. lincludebbfinaldo)) then
                  fam2end=2
               else
                  fam2end=1
               endif
               do fam2=0,fam2end !here add b quark
                  ! the two quarks on coupling to the w
                  u = 2*fam1           !u,c quarks
                  d = u - 1            !d,s quarks
                  ! the observer quark
                  q = fl2 + 2*fam2     ! u,d,s,c quarks

!                 print*, 'fam'
!                 print*, fam1,fl2,fam2
!                 print*, u,d,q,k
                  if (q.eq.u) then
!                     print*, q, u

                     ! 1/36d0 : (N_c**2) * (2**2) (average over quark color and spin)
                     ! interference terms uuud
                     ! TODO: Colorstructure
                     ! flavor output does not include ckm effects here
                     !   ubar ubar -> ubar dbar could also be ubar sbar
                     !   only final state d quarks effected (in uuud)
                     !   (or u in uddd)
                     ! in this section with the pdfs the label q is not
                     ! always the observer but symbolises a quark that
                     ! should not get ckm changes
                     k=FL_VVgg((/-u,-u,-d,-u/),0); res(k)=m2s_f(1,1) *getpdf(pdf,-q,-q)/36d0     
                     k=FL_VVgg((/-u, u,-d, u/),0); res(k)=m2s_f(2,1) *getpdf(pdf,-q, q)/36d0
                     k=FL_VVgg((/ u,-u,-d, u/),0); res(k)=m2s_ch(2,1)*getpdf(pdf, q,-q)/36d0
                     k=FL_VVgg((/-u, d, u,-u/),0); res(k)=m2s_f(3,1) *getpdfckm2(pdf,pdfckm,-q, d)/36d0
                     k=FL_VVgg((/ d,-u, u,-u/),0); res(k)=m2s_ch(3,1)*getpdfckm1(pdf,pdfckm, d,-q)/36d0
                     k=FL_VVgg((/ u, d, u, u/),0); res(k)=m2s_f(4,1) *getpdfckm2(pdf,pdfckm, q, d)/36d0*0.5d0
                     k=FL_VVgg((/ d, u, u, u/),0); res(k)=m2s_ch(4,1)*getpdfckm1(pdf,pdfckm, d, q)/36d0*0.5d0
                     ! additional 0.5 for same-quark final state (no interchange possible)
                  elseif (q.eq.d) then                             
                     ! interference terms uddd
                     k=FL_VVgg((/-u,-d,-d,-d/),0); res(k)=m2s_dd(1,1)  *getpdfckm1(pdf,pdfckm,-u,-q)/36d0*0.5d0
                     k=FL_VVgg((/-d,-u,-d,-d/),0); res(k)=m2s_ddch(1,1)*getpdfckm2(pdf,pdfckm,-q,-u)/36d0*0.5d0
                     k=FL_VVgg((/-u, d,-d, d/),0); res(k)=m2s_dd(2,1)  *getpdfckm1(pdf,pdfckm,-u, q)/36d0
                     k=FL_VVgg((/ d,-u,-d, d/),0); res(k)=m2s_ddch(2,1)*getpdfckm2(pdf,pdfckm, q,-u)/36d0
                     k=FL_VVgg((/-d, d,-u, d/),0); res(k)=m2s_dd(3,1)  *getpdf(pdf,-q, q)/36d0
                     k=FL_VVgg((/ d,-d,-u, d/),0); res(k)=m2s_ddch(3,1)*getpdf(pdf, q,-q)/36d0
                     k=FL_VVgg((/ d, d, u, d/),0); res(k)=m2s_dd(4,1)  *getpdf(pdf, q, q)/36d0
                  else
                     ! non-interference terms
                     ! here: only ckm for u,d not q
                     if (fl2.eq.1) then !down-type observer

                        if (fam2.ne.2 .or. lincludebbfinaldo .or. lincludebb) then
                        ! final state observer
                        ! take bb ohne if includebb or includebbfinal set
                      k=FL_VVgg((/-u, d,-q, q/),0); res(k)=  m2s_dd(2,2)  *getpdfckm1(pdf,pdfckm,-u, d)/36d0
                      k=FL_VVgg((/ d,-u,-q, q/),0); res(k)=  m2s_ddch(2,2)*getpdfckm2(pdf,pdfckm, d,-u)/36d0
                        endif

                        if ((fam2.eq.2).and. .not. (lincludebbfinaldo .or. lincludebb)) then
                           cycle !don't do the rest of the contributions (which include initial b quarks)
                        else
                           if (fam2.eq.2) then
                              ! b bar -> u d is double counted, halfed here by hand
                              ! remember to increase size of res if anything gets changed here
                      k=fl_vvgg((/-q, q, u,-d/),0); res(k)=  m2s_dd(3,2)  *getpdf(pdf,-q, q)/36d0*0.5d0
                      k=fl_vvgg((/ q,-q, u,-d/),0); res(k)=  m2s_ddch(3,2)*getpdf(pdf, q,-q)/36d0*0.5d0
                           else
                      k=fl_vvgg((/-q, q, u,-d/),0); res(k)=  m2s_dd(3,2)  *getpdf(pdf,-q, q)/36d0  
                      k=fl_vvgg((/ q,-q, u,-d/),0); res(k)=  m2s_ddch(3,2)*getpdf(pdf, q,-q)/36d0
                           endif
                      k=FL_VVgg((/-u,-q,-d,-q/),0); res(k)=  m2s_dd(1,3)  *getpdfckm1(pdf,pdfckm,-u,-q)/36d0
                      k=FL_VVgg((/-q,-u,-d,-q/),0); res(k)=  m2s_ddch(1,3)*getpdfckm2(pdf,pdfckm,-q,-u)/36d0
                      k=FL_VVgg((/-u, q,-d, q/),0); res(k)=  m2s_dd(2,3)  *getpdfckm1(pdf,pdfckm,-u, q)/36d0
                      k=FL_VVgg((/ q,-u,-d, q/),0); res(k)=  m2s_ddch(2,3)*getpdfckm2(pdf,pdfckm, q,-u)/36d0
                      k=FL_VVgg((/-q, d, u,-q/),0); res(k)=  m2s_dd(3,3)  *getpdf(pdf,-q, d)/36d0
                      k=FL_VVgg((/ d,-q, u,-q/),0); res(k)=  m2s_ddch(3,3)*getpdf(pdf, d,-q)/36d0
                      k=FL_VVgg((/ q, d, u, d/),0); res(k)=  m2s_dd(4,3)  *getpdf(pdf, q, d)/36d0
                      k=FL_VVgg((/ d, q, u, d/),0); res(k)=  m2s_ddch(4,3)*getpdf(pdf, d, q)/36d0
                        endif
                     else if (fl2.eq.2) then !up-type observer             
                      k=FL_VVgg((/-q, q,-d, u/),0); res(k)=  m2s_f(2,2) *getpdf(pdf,-q, q)/36d0
                      k=FL_VVgg((/ q,-q,-d, u/),0); res(k)=  m2s_ch(2,2)*getpdf(pdf, q,-q)/36d0
                      k=FL_VVgg((/-q,-u,-d,-u/),0); res(k)=  m2s_f(1,2) *getpdf(pdf,-q,-u)/36d0
                      k=FL_VVgg((/-u,-q,-d,-u/),0); res(k)=  m2s_ch(1,2)*getpdf(pdf,-u,-q)/36d0
                      k=FL_VVgg((/-u, q,-d, q/),0); res(k)=  m2s_f(2,3) *getpdf(pdf,-u, q)/36d0
                      k=FL_VVgg((/ q,-u,-d, q/),0); res(k)=  m2s_ch(2,3)*getpdf(pdf, q,-u)/36d0
                      k=FL_VVgg((/-u, d, q,-q/),0); res(k)=  m2s_f(3,2) *getpdfckm2(pdf,pdfckm,-u, d)/36d0
                      k=FL_VVgg((/ d,-u, q,-q/),0); res(k)=  m2s_ch(3,2)*getpdfckm1(pdf,pdfckm, d,-u)/36d0
                      k=FL_VVgg((/-q, d, u,-q/),0); res(k)=  m2s_f(3,3) *getpdfckm2(pdf,pdfckm,-q, d)/36d0
                      k=FL_VVgg((/ d,-q, u,-q/),0); res(k)=  m2s_ch(3,3)*getpdfckm1(pdf,pdfckm, d,-q)/36d0
                      k=FL_VVgg((/ q, d, q, u/),0); res(k)=  m2s_f(4,3) *getpdfckm2(pdf,pdfckm, q, d)/36d0
                      k=FL_VVgg((/ d, q, q, u/),0); res(k)=  m2s_ch(4,3)*getpdfckm1(pdf,pdfckm, d, q)/36d0
                     endif
                  endif
               enddo
            enddo
         enddo

      endif

      if(.not. doSubprocesses .or. (sub_number.eq.2) )  then
         ! diagrams qqgg

         do fam1=1,2  !coupling to w
            ! the observer is a gluon
            ! the two quarks on coupling to the w
            u = 2*fam1           !u,c quarks
            d = u - 1            !d,s quarks

            k=FL_VVgg((/-u, d, 0, 0/),0); res(k)=  mgg(1)  *getpdfckm1(pdf,pdfckm,-u, d)/72d0 
            k=FL_VVgg((/ d,-u, 0, 0/),0); res(k)=  mggch(1)*getpdfckm2(pdf,pdfckm, d,-u)/72d0
            !(N_c**2) * (2**2) *2 (average over quark color and spin) and
            !interchange gluons
            k=FL_VVgg((/-u, 0,-d, 0/),0); res(k)=  mgg(2)  *getpdf(pdf,-u, 0)/96d0
            k=FL_VVgg((/ 0,-u,-d, 0/),0); res(k)=  mggch(2)*getpdf(pdf, 0,-u)/96d0
            !(N_c) * 2)  (average over quark color and spin) 
            !*(N_c**2 - 1) * 2 (average over gluon and spin) 
            k=FL_VVgg((/ d, 0, u, 0/),0); res(k)=  mgg(3)  *getpdf(pdf, d, 0)/96d0
            k=FL_VVgg((/ 0, d, u, 0/),0); res(k)=  mggch(3)*getpdf(pdf, 0, d)/96d0
            !(N_c) * 2)  (average over quark color and spin) 
            !*(N_c**2 - 1) * 2 (average over gluon and spin) 
            k=FL_VVgg((/ 0, 0, u,-d/),0); res(k)=  mgg(4)  *getpdf(pdf, 0, 0)/256d0
            !(N_c**2 - 1)**2 * 2**2 (average over gluon and spin) 

            ! remember to increase size of res if anything gets changed here
         enddo


      endif ! sub_number 2

      !print*, 'max k, size of res?',k
      amp = zero
      do k=1,FL_VVgg((/0,0,0,0/),-2) !gets number of subprocesses
         amp = amp+res(k)
      enddo
      m2s_fey=amp

      endif !lokt(1)

c...Les Houches interface - the most propable subprocess   
      if ((lha.or.hepmc).and..not.doNLO) then
         i=0
         weight=0.d0
         rnumb=RandomNumber()
         ! k now is number of different processes
         do while((i.le.k).and.(weight.le.rnumb*m2s_fey))
            i=i+1
            weight=weight+res(i)
            iprocess=i
         enddo
      endif

! PDF(5, 4, 3, 2, 1, 0, -1, ......, -5)
! for(b, c, s, u, d, g, d_bar, ..., b_bar)
! (switched u <-> d comp. to CTEQ)


! Kinematics 1

      if (lokt(2)) then !.and.dotrr(p(0,1,1),p(0,4,1)).le.1d0) then   ! Dipole 1,43

         call getyourscalesready_dibos(xi,pdfchange,pdftmp,pdftmp,2)

         if(.not. doSubprocesses .or. (sub_number.eq.1) )  then
         Call daisjdipole67_dib(2,1,1,xuz(1,2),p(0,4,1),p(0,3,1),p(0,1,1),p(0,1,2),v(0,1,2),dipolea(1))
         Call daisjdipole67_dib(2,2,1,xuz(1,2),p(0,4,1),p(0,3,1),p(0,1,1),p(0,1,2),v(0,1,2),dipoleb(1))

         dipole(1)=dipolea(1)
         dipole(2)=dipoleb(1)

         do fam1=1,2 
            do fl2=1,2
               if (fl2.eq.1 .and. lincludebb) then
                  fam2end=2
               else
                  fam2end=1
               endif
               do fam2=0,fam2end 
                  u = 2*fam1           !u,c quarks
                  d = u - 1            !d,s quarks
                  q = fl2 + 2*fam2     !d,u,s,c,b quarks
                  if (q.eq.u) then
                     kinem(1) = kinem(1)
     &               + dipole(1)*(getpdf(pdf,-q,-q)+getpdf(pdf, q,-q))
     &               + dipole(2)*(getpdfckm2(pdf,pdfckm, q, d)*0.5d0+getpdf(pdf,-q, d))
                  else if (q.eq.d) then
                     kinem(1) = kinem(1)
     &               + dipole(1)*(getpdfckm2(pdf,pdfckm,-q,-u)*0.5d0+getpdf(pdf, q,-u))
     &               + dipole(2)*(getpdf(pdf, q, q)+getpdf(pdf,-q, q))
                  else 
                     kinem(1) = kinem(1)
     &               + dipole(1)*(getpdf(pdf,-q,-u)+getpdf(pdf, q,-u))
     &               + dipole(2)*(getpdf(pdf, q, d)+getpdf(pdf,-q, d))
                  endif
               enddo
            enddo
         enddo

         endif

         if(.not. doSubprocesses .or. (sub_number.eq.2) )  then

         Call daisjdipole65_dib(2,1,1,xuz(1,2),p(0,4,1),p(0,3,1),p(0,1,1),p(0,1,2),v(0,1,2),dipolev(1,1))
         Call daisjdipole65_dib(2,1,2,xuz(1,2),p(0,4,1),p(0,3,1),p(0,1,1),p(0,1,2),v(0,1,2),dipolev(1,2))
         Call daisjdipole65_dib(2,2,1,xuz(1,2),p(0,4,1),p(0,3,1),p(0,1,1),p(0,1,2),v(0,1,2),dipolev(1,3))

         do fam1=1,2 
            u = 2*fam1           !u,c quarks
            d = u - 1            !d,s quarks
            kinem(1) = kinem(1)
     &       + dipolev(1,1)*0.5d0*getpdfckm1(pdf,pdfckm,-u, d)
     &       + dipolev(1,2)*0.5d0*getpdfckm2(pdf,pdfckm, d,-u)
     &       + dipolev(4,2)*getpdf(pdf, 0,-u)
     &       + dipolev(4,1)*getpdf(pdf, 0, d)
     &       + dipolev(5,3)*getpdf(pdf, 0, 0)
         enddo

         endif
      endif

! Kinematics 2

      if (lokt(3)) then ! .and.dotrr(p(0,2,1),p(0,4,1)).le.1d0) then   ! Dipole 2,43

         call getyourscalesready_dibos(xi,pdfchange,pdftmp,pdftmp,3)

         if(.not. doSubprocesses .or. (sub_number.eq.1) )  then
         Call daisjdipole67_dib(3,1,1,xuz(1,3),p(0,4,1),p(0,3,1),p(0,2,1),p(0,1,3),v(0,1,3),dipolea(1))
         Call daisjdipole67_dib(3,2,1,xuz(1,3),p(0,4,1),p(0,3,1),p(0,2,1),p(0,1,3),v(0,1,3),dipoleb(1))

         dipole(1)=dipolea(1)
         dipole(2)=dipoleb(1)

         do fam1=1,2 
            do fl2=1,2
               if (fl2.eq.1 .and. lincludebb) then
                  fam2end=2
               else
                  fam2end=1
               endif
               do fam2=0,fam2end 
                  u = 2*fam1           !u,c quarks
                  d = u - 1            !d,s quarks
                  q = fl2 + 2*fam2     !d,u,s,c,b quarks
                  if (q.eq.u) then
                     kinem(2) = kinem(2)
     &               + dipole(1)*(getpdf(pdf,-q,-q)+getpdf(pdf,-q, q))
     &               + dipole(2)*(getpdfckm1(pdf,pdfckm, d, q)*0.5d0+getpdf(pdf, d,-q))
                  else if (q.eq.d) then
                     kinem(2) = kinem(2)
     &               + dipole(1)*(getpdfckm1(pdf,pdfckm,-u,-q)*0.5d0+getpdf(pdf,-u, q))
     &               + dipole(2)*(getpdf(pdf, q, q)+getpdf(pdf, q,-q))
                  else 
                     kinem(2) = kinem(2)
     &               + dipole(1)*(getpdf(pdf,-u,-q)+getpdf(pdf,-u, q))
     &               + dipole(2)*(getpdf(pdf, d, q)+getpdf(pdf, d,-q))
                  endif
               enddo
            enddo
         enddo
         endif

         if(.not. doSubprocesses .or. (sub_number.eq.2) )  then
         Call daisjdipole65_dib(3,1,2,xuz(1,3),p(0,4,1),p(0,3,1),p(0,2,1),p(0,1,3),v(0,1,3),dipolev(1,1))
         Call daisjdipole65_dib(3,1,1,xuz(1,3),p(0,4,1),p(0,3,1),p(0,2,1),p(0,1,3),v(0,1,3),dipolev(1,2))
         Call daisjdipole65_dib(3,2,1,xuz(1,3),p(0,4,1),p(0,3,1),p(0,2,1),p(0,1,3),v(0,1,3),dipolev(1,3))

         do fam1=1,2 
            u = 2*fam1           !u,c quarks
            d = u - 1            !d,s quarks
            kinem(2) = kinem(2)
     &       + dipolev(1,1)*0.5d0*getpdfckm1(pdf,pdfckm,-u, d)
     &       + dipolev(1,2)*0.5d0*getpdfckm2(pdf,pdfckm, d,-u)
     &       + dipolev(4,1)*getpdf(pdf,-u, 0)
     &       + dipolev(4,2)*getpdf(pdf, d, 0)
     &       + dipolev(5,3)*getpdf(pdf, 0, 0)
         enddo
         endif
      endif

! Kinematics 3

      if (lokt(4)) then !.and.dotrr(p(0,1,1),p(0,4,1)).le.1d0) then   ! Dipole 12,4

         call getyourscalesready_dibos(xi,pdfchange,pdftmp,pdftmp,4)

         Call daibdipole147_dib(4,1,1,xuz(1,4),p(0,4,1),p(0,1,1),p(0,2,1),p(0,1,4),v(0,1,4),dipolea(1))
         Call daibdipole147_dib(4,2,1,xuz(1,4),p(0,4,1),p(0,1,1),p(0,2,1),p(0,1,4),v(0,1,4),dipoleb(1))
         if(.not. doSubprocesses .or. (sub_number.eq.1) )  then

         dipole(1)=dipolea(1)
         dipole(2)=dipoleb(1)

         do fam1=1,2 
            do fl2=1,2
               if (fl2.eq.1 .and. lincludebb) then
                  fam2end=2
               else
                  fam2end=1
               endif
               do fam2=0,fam2end 
                  u = 2*fam1           !u,c quarks
                  d = u - 1            !d,s quarks
                  q = fl2 + 2*fam2     !d,u,s,c,b quarks
                  if (q.eq.u) then
                     kinem(3) = kinem(3)
     &               + dipole(1)*(getpdf(pdf,-q,-q)+getpdf(pdf, q,-q))
     &               + dipole(2)*(getpdfckm2(pdf,pdfckm, q, d)*0.5d0+getpdf(pdf,-q, d))
                  else if (q.eq.d) then
                     kinem(3) = kinem(3)
     &               + dipole(1)*(getpdfckm2(pdf,pdfckm,-q,-u)*0.5d0+getpdf(pdf, q,-u))
     &               + dipole(2)*(getpdf(pdf, q, q)+getpdf(pdf,-q, q))
                  else 
                     kinem(3) = kinem(3)
     &               + dipole(1)*(getpdf(pdf,-q,-u)+getpdf(pdf, q,-u))
     &               + dipole(2)*(getpdf(pdf, q, d)+getpdf(pdf,-q, d))
                  endif
               enddo
            enddo
         enddo

         endif

         if(.not. doSubprocesses .or. (sub_number.eq.2) )  then
         Call daibdipole145_dib(4,1,1,xuz(1,4),p(0,4,1),p(0,1,1),p(0,2,1),p(0,1,4),v(0,1,4),dipolev(1,1))
         Call daibdipole145_dib(4,1,2,xuz(1,4),p(0,4,1),p(0,1,1),p(0,2,1),p(0,1,4),v(0,1,4),dipolev(1,2))
         Call daibdipole145_dib(4,2,1,xuz(1,4),p(0,4,1),p(0,1,1),p(0,2,1),p(0,1,4),v(0,1,4),dipolev(1,3))

         do fam1=1,2 
            u = 2*fam1           !u,c quarks
            d = u - 1            !d,s quarks
            kinem(3) = kinem(3)
     &       + dipolev(1,2)*0.5d0*getpdfckm2(pdf,pdfckm, d,-u)
     &       + dipolev(1,1)*0.5d0*getpdfckm1(pdf,pdfckm,-u, d)
     &       + dipolev(4,2)*getpdf(pdf, 0,-u)
     &       + dipolev(4,1)*getpdf(pdf, 0, d)
     &       + dipolev(5,3)*getpdf(pdf, 0, 0)
         enddo
         endif

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

         call getyourscalesready_dibos(xi,pdfchange,pdftmp,pdftmp,5)

         Call daibdipole147_dib(5,1,1,xuz(1,5),p(0,4,1),p(0,2,1),p(0,1,1),p(0,1,5),v(0,1,5),dipolea(1))
         Call daibdipole147_dib(5,2,1,xuz(1,5),p(0,4,1),p(0,2,1),p(0,1,1),p(0,1,5),v(0,1,5),dipoleb(1))
         if(.not. doSubprocesses .or. (sub_number.eq.1) )  then
            dipole(1)=dipolea(1)
            dipole(2)=dipoleb(1)

         do fam1=1,2 
            do fl2=1,2
               if (fl2.eq.1 .and. lincludebb) then
                  fam2end=2
               else
                  fam2end=1
               endif
               do fam2=0,fam2end 
                  u = 2*fam1           !u,c quarks
                  d = u - 1            !d,s quarks
                  q = fl2 + 2*fam2     !d,u,s,c,b quarks
                  if (q.eq.u) then
                     kinem(4) = kinem(4)
     &               + dipole(1)*(getpdf(pdf,-q,-q)+getpdf(pdf,-q, q))
     &               + dipole(2)*(getpdfckm1(pdf,pdfckm, d, q)*0.5d0+getpdf(pdf, d,-q))
                  else if (q.eq.d) then
                     kinem(4) = kinem(4)
     &               + dipole(1)*(getpdfckm1(pdf,pdfckm,-u,-q)*0.5d0+getpdf(pdf,-u, q))
     &               + dipole(2)*(getpdf(pdf, q, q)+getpdf(pdf, q,-q))
                  else 
                     kinem(4) = kinem(4)
     &               + dipole(1)*(getpdf(pdf,-u,-q)+getpdf(pdf,-u, q))
     &               + dipole(2)*(getpdf(pdf, d, q)+getpdf(pdf, d,-q))
                  endif
               enddo
            enddo
         enddo
         endif
         if(.not. doSubprocesses .or. (sub_number.eq.2) )  then
         Call daibdipole145_dib(5,1,2,xuz(1,5),p(0,4,1),p(0,2,1),p(0,1,1),p(0,1,5),v(0,1,5),dipolev(1,1))
         Call daibdipole145_dib(5,1,1,xuz(1,5),p(0,4,1),p(0,2,1),p(0,1,1),p(0,1,5),v(0,1,5),dipolev(1,2))
         Call daibdipole145_dib(5,2,1,xuz(1,5),p(0,4,1),p(0,2,1),p(0,1,1),p(0,1,5),v(0,1,5),dipolev(1,3))

         do fam1=1,2 
            u = 2*fam1           !u,c quarks
            d = u - 1            !d,s quarks
            kinem(4) = kinem(4)
     &       + dipolev(1,1)*0.5d0*getpdfckm1(pdf,pdfckm,-u, d)
     &       + dipolev(1,2)*0.5d0*getpdfckm2(pdf,pdfckm, d,-u)
     &       + dipolev(4,1)*getpdf(pdf,-u, 0)
     &       + dipolev(4,2)*getpdf(pdf, d, 0)
     &       + dipolev(5,3)*getpdf(pdf, 0, 0)
         enddo
         endif

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

         call getyourscalesready_dibos(xi,pdfchange,pdftmp,pdftmp,6)

         Call daisjdipole67_dib(6,1,1,xuz(1,6),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipolea(1))
         Call daisjdipole67_dib(6,2,1,xuz(1,6),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipoleb(1))

         dipole(1)=dipolea(1)
         dipole(2)=dipoleb(1)
         dipole(4)=dipolea(2)
         dipole(6)=dipoleb(2)

         if(.not. doSubprocesses .or. (sub_number.eq.1) )  then

         do fam1=1,2 
            u = 2*fam1           !u,c quarks
            d = u - 1            !d,s quarks
            kinem(5) = kinem(5)
     &               + dipole(2)*(getpdfckm2(pdf,pdfckm, u, d)*0.5d0) !u d  -> u u 
     &               + dipole(1)*(getpdfckm2(pdf,pdfckm,-d,-u)*0.5d0) !d~u~ -> d~d~
         enddo
! !--------------------------------------------------------------------------------------------------------------
         endif
         if(.not. doSubprocesses .or. (sub_number.eq.2) )  then
         Call daisjdipole65_dib(6,1,1,xuz(1,6),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipolev(1,1))
         Call daisjdipole65_dib(6,1,2,xuz(1,6),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipolev(1,2))
         Call daisjdipole65_dib(6,3,1,xuz(1,6),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipolev(1,3))
         
         Call daisjdipole67_dib(6,1,2,xuz(1,6),p(0,3,1), p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipolea(1))
         Call daisjdipole67_dib(6,2,2,xuz(1,6),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,6),v(0,1,6),dipoleb(1))

            dipole(3)=dipolea(3)  
            dipole(5)=dipoleb(3)

         do fam1=1,2 
            u = 2*fam1           !u,c quarks
            d = u - 1            !d,s quarks
            kinem(5) = kinem(5)
     &       + dipolev(1,1)*0.5d0*getpdfckm1(pdf,pdfckm,-u, d)
     &       + dipolev(1,2)*0.5d0*getpdfckm2(pdf,pdfckm, d,-u)
     &       + dipole(3)*getpdf(pdf,-u, 0)
     &       + dipole(4)*getpdf(pdf, 0,-u)
     &       + dipole(5)*getpdf(pdf, d, 0)
     &       + dipole(6)*getpdf(pdf, 0, d)
     &       + dipolev(5,3)*getpdf(pdf, 0, 0)
         enddo
         endif

      endif

! Kinematics 6

      if (lokt(7)) then !.and.dotrr(p(0,2,1),p(0,3,1)).le.1d0) then   ! Dipole 2,34

         call getyourscalesready_dibos(xi,pdfchange,pdftmp,pdftmp,7)

         Call daisjdipole67_dib(7,1,1,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipolea(1))
         Call daisjdipole67_dib(7,2,1,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipoleb(1))
         dipole(1)=dipolea(1)
         dipole(2)=dipoleb(1)
         dipole(3)=dipolea(2)
         dipole(5)=dipoleb(2)
         
         if(.not. doSubprocesses .or. (sub_number.eq.1)) then

         do fam1=1,2 
            u = 2*fam1           !u,c quarks
            d = u - 1            !d,s quarks
            kinem(6) = kinem(6)
     &               + dipole(2)*(getpdfckm1(pdf,pdfckm, d, u)*0.5d0) ! u d  -> u u 
     &               + dipole(1)*(getpdfckm1(pdf,pdfckm,-u,-d)*0.5d0) ! d~u~ -> d~d~
         enddo
! !-------------------------------------------------------------------------------------------------------------
         endif
         if(.not. doSubprocesses .or. (sub_number .eq. 2)) then
         Call daisjdipole65_dib(7,1,2,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipolev(1,1))
         Call daisjdipole65_dib(7,1,1,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipolev(1,2))
         Call daisjdipole65_dib(7,3,1,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipolev(1,3))
         
         Call daisjdipole67_dib(7,1,2,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipolea(1))
         Call daisjdipole67_dib(7,2,2,xuz(1,7),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,7),v(0,1,7),dipoleb(1))

         dipole(4)=dipolea(3)
         dipole(6)=dipoleb(3)

         do fam1=1,2 
            u = 2*fam1           !u,c quarks
            d = u - 1            !d,s quarks
            kinem(6) = kinem(6)
     &       + dipolev(1,1)*0.5d0*getpdfckm1(pdf,pdfckm,-u, d)
     &       + dipolev(1,2)*0.5d0*getpdfckm2(pdf,pdfckm, d,-u)
     &       + dipole(3)*getpdf(pdf,-u, 0)
     &       + dipole(4)*getpdf(pdf, 0,-u)
     &       + dipole(5)*getpdf(pdf, d, 0)
     &       + dipole(6)*getpdf(pdf, 0, d)
     &       + dipolev(5,3)*getpdf(pdf, 0, 0)
         enddo
         endif
      endif

! Kinematics 7

      if (lokt(8)) then !.and.dotrr(p(0,1,1),p(0,3,1)).le.1d0) then   ! Dipole 12,3

         call getyourscalesready_dibos(xi,pdfchange,pdftmp,pdftmp,8)

         Call daibdipole147_dib(8,1,1,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipolea(1))
         Call daibdipole147_dib(8,2,1,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipoleb(1))

         dipole(1)=dipolea(1)
         dipole(2)=dipoleb(1)
         dipole(4)=dipolea(2)
         dipole(6)=dipoleb(2)

         if(.not. doSubprocesses .or. (sub_number .eq. 1)) then
         do fam1=1,2 
            u = 2*fam1           !u,c quarks
            d = u - 1            !d,s quarks
            kinem(7) = kinem(7)
     &               + dipole(2)*(getpdfckm2(pdf,pdfckm, u, d)*0.5d0) ! u d  -> u u 
     &               + dipole(1)*(getpdfckm2(pdf,pdfckm,-d,-u)*0.5d0) ! d~u~ -> d~d~
         enddo
! !-------------------------------------------------------------------------------------------------------------
         endif
         if(.not. doSubprocesses .or. (sub_number .eq. 2)) then
         Call daibdipole147_dib(8,1,2,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipolea(1))
         Call daibdipole147_dib(8,2,2,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipoleb(1))

         dipole(3)=dipolea(3)
         dipole(5)=dipoleb(3)
         
         Call daibdipole145_dib(8,1,1,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipolev(1,1))
         Call daibdipole145_dib(8,1,2,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipolev(1,2))
         Call daibdipole145_dib(8,3,1,xuz(1,8),p(0,3,1),p(0,1,1),p(0,2,1),p(0,1,8),v(0,1,8),dipolev(1,3))
         do fam1=1,2 
            u = 2*fam1           !u,c quarks
            d = u - 1            !d,s quarks
            kinem(7) = kinem(7)
     &       + dipolev(1,1)*0.5d0*getpdfckm1(pdf,pdfckm,-u, d)
     &       + dipolev(1,2)*0.5d0*getpdfckm2(pdf,pdfckm, d,-u)
     &       + dipole(3)*getpdf(pdf,-u, 0)
     &       + dipole(4)*getpdf(pdf, 0,-u)
     &       + dipole(5)*getpdf(pdf, d, 0)
     &       + dipole(6)*getpdf(pdf, 0, d)
     &       + dipolev(5,3)*getpdf(pdf, 0, 0)
         enddo

         endif
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

         call getyourscalesready_dibos(xi,pdfchange,pdftmp,pdftmp,9)

         Call daibdipole147_dib(9,1,1,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipolea(1))
         Call daibdipole147_dib(9,2,1,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipoleb(1))

         dipole(1)=dipolea(1)
         dipole(2)=dipoleb(1)
         dipole(3)=dipolea(2)
         dipole(5)=dipoleb(2)


         if(.not. doSubprocesses .or. (sub_number .eq. 1)) then
         do fam1=1,2 
            u = 2*fam1           !u,c quarks
            d = u - 1            !d,s quarks
            kinem(8) = kinem(8)
     &               + dipole(2)*(getpdfckm1(pdf,pdfckm, d, u)*0.5d0) ! u d  -> u u 
     &               + dipole(1)*(getpdfckm1(pdf,pdfckm,-u,-d)*0.5d0) ! d~u~ -> d~d~
         enddo
! !-------------------------------------------------------------------------------------------------------------
         endif
         if(.not. doSubprocesses .or. (sub_number .eq. 2)) then
         Call daibdipole145_dib(9,1,2,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipolev(1,1))
         Call daibdipole145_dib(9,1,1,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipolev(1,2))

         Call daibdipole147_dib(9,1,2,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipolea(1))
         Call daibdipole147_dib(9,2,2,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipoleb(1))
         Call daibdipole145_dib(9,3,1,xuz(1,9),p(0,3,1),p(0,2,1),p(0,1,1),p(0,1,9),v(0,1,9),dipolev(1,3))

         dipole(4)=dipolea(3)
         dipole(6)=dipoleb(3)

         do fam1=1,2 
            u = 2*fam1           !u,c quarks
            d = u - 1            !d,s quarks
            kinem(8) = kinem(8)
     &       + dipolev(1,1)*0.5d0*getpdfckm(pdf,pdfckm,-u, d)
     &       + dipolev(1,2)*0.5d0*getpdfckm(pdf,pdfckm, d,-u)
     &       + dipole(3)*getpdf(pdf,-u, 0)
     &       + dipole(4)*getpdf(pdf, 0,-u)
     &       + dipole(5)*getpdf(pdf, d, 0)
     &       + dipole(6)*getpdf(pdf, 0, d)
     &       + dipolev(5,3)*getpdf(pdf, 0, 0)
         enddo
         endif

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

         call getyourscalesready_dibos(xi,pdfchange,pdftmp,pdftmp,10)

         if(.not. doSubprocesses .or. (sub_number .eq. 1)) then
         CALL dasijdipole40_dib(10,1,4,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1))
         CALL dasijdipole40_dib(10,2,4,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1))

         dipole(1)=dipolea(1)
         dipole(2)=dipoleb(1)
         
         kinem(9)=zero
         do fam1=1,2 
            u = 2*fam1           !u,c quarks
            d = u - 1            !d,s quarks
            kinem(9) = kinem(9) 
     &        + dipole(1)*getpdfckm1(pdf,pdfckm,-u, d)
     &        + dipole(2)*getpdfckm2(pdf,pdfckm, d,-u)
         enddo
         if(lincludebb .or. lincludebbfinaldo) then
            kinem(9) = kinem(9) * 5d0 ! number of q qbar final states (with b)
         else
            kinem(9) = kinem(9) * 4d0 !(without b)
         endif
         endif
         if(.not. doSubprocesses .or. (sub_number .eq. 2)) then
         CALL dasijdipole40_dib(10,1,1,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1))
         CALL dasijdipole40_dib(10,2,1,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1))

         dipole(7)=dipolea(1)
         dipole(8)=dipoleb(1)

! original:         
                 CALL test40_dib(10,1,2,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1))
                 CALL test40_dib(10,2,2,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1))

!         CALL dasijdipole40_dib(10,1,2,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1))
!         CALL dasijdipole40_dib(10,2,2,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1))

         dipole(3)=dipolea(2)
         dipole(4)=dipoleb(1)

! those two calls are the same (momenta 3,4 must be interchanged!)        
!        CALL dasijdipole40_dib(10,1,3,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1))
!               CALL test40_dib(10,1,3,xuz(1,10),p(0,4,1),p(0,3,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1))

!        CALL dasijdipole40_dib(10,2,3,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1))
!               CALL test40_dib(10,2,3,xuz(1,10),p(0,4,1),p(0,3,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1))

!original:         
         CALL test40_dib(10,1,3,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1))
         CALL test40_dib(10,2,3,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1))

!        CALL dasijdipole40_dib(10,1,3,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipolea(1))
!        CALL dasijdipole40_dib(10,2,3,xuz(1,10),p(0,3,1),p(0,4,1),p(0,1,1),p(0,1,10),v(0,1,10),dipoleb(1))

         dipole(5)=dipolea(2)
         dipole(6)=dipoleb(1)

         do fam1=1,2 
            u = 2*fam1           !u,c quarks
            d = u - 1            !d,s quarks
            kinem(9) = kinem(9)
     &       + dipole(3)*getpdf(pdf,-u, 0)
     &       + dipole(4)*getpdf(pdf, 0,-u)
     &       + dipole(5)*getpdf(pdf, d, 0)
     &       + dipole(6)*getpdf(pdf, 0, d)
     &       + dipole(7)*getpdfckm(pdf,pdfckm,-u, d)*0.5d0
     &       + dipole(8)*getpdfckm(pdf,pdfckm, d,-u)*0.5d0
         enddo
         endif

      endif


! Kinematics 10

      if (lokt(11)) then !.and.dotrr(p(0,4,1),p(0,3,1)).le.1d0) then   ! Dipole 2,34

         call getyourscalesready_dibos(xi,pdfchange,pdftmp,pdftmp,11)

         if(.not. doSubprocesses .or. (sub_number .eq. 1)) then
         CALL dasijdipole40_dib(11,1,4,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipolea(1))
         CALL dasijdipole40_dib(11,2,4,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipoleb(1))

         dipole(1)=dipolea(1)
         dipole(2)=dipoleb(1)

         kinem(10)=zero
         do fam1=1,2 
            u = 2*fam1           !u,c quarks
            d = u - 1            !d,s quarks
            kinem(10) = kinem(10) 
     &        + dipole(1)*getpdfckm1(pdf,pdfckm,-u, d)
     &        + dipole(2)*getpdfckm2(pdf,pdfckm, d,-u)
         enddo
         if(lincludebb .or. lincludebbfinaldo) then
            kinem(10) = kinem(10) * 5d0 ! number of q qbar final states (with b)
         else
            kinem(10) = kinem(10) * 4d0 !(without b)
         endif
         endif
         if(.not. doSubprocesses .or. (sub_number .eq. 2)) then
         CALL dasijdipole40_dib(11,1,1,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipolea(1))
         CALL dasijdipole40_dib(11,2,1,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipoleb(1))

         dipole(7)=dipolea(1)
         dipole(8)=dipoleb(1)

         ! see kinem(9) for alternative calls
!       CALL dasijdipole40_dib(11,1,2,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipolea(1))
!       CALL dasijdipole40_dib(11,2,2,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipoleb(1))

               CALL test40_dib(11,1,2,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipolea(1))
               CALL test40_dib(11,2,2,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipoleb(1))
         dipole(3)=dipolea(1)
         dipole(4)=dipoleb(2)

!       CALL dasijdipole40_dib(11,1,3,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipolea(1))
!       CALL dasijdipole40_dib(11,2,3,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipoleb(1))

         CALL test40_dib(11,1,3,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipolea(1))
         CALL test40_dib(11,2,3,xuz(1,11),p(0,3,1),p(0,4,1),p(0,2,1),p(0,1,11),v(0,1,11),dipoleb(1))
! 
         dipole(5)=dipolea(1)
         dipole(6)=dipoleb(2)

         do fam1=1,2 
            u = 2*fam1           !u,c quarks
            d = u - 1            !d,s quarks
            kinem(10) = kinem(10)
     &       + dipole(3)*getpdf(pdf,-u, 0)
     &       + dipole(4)*getpdf(pdf, 0,-u)
     &       + dipole(5)*getpdf(pdf, d, 0)
     &       + dipole(6)*getpdf(pdf, 0, d)
     &       + dipole(7)*getpdfckm(pdf,pdfckm,-u, d)*0.5d0
     &       + dipole(8)*getpdfckm(pdf,pdfckm, d,-u)*0.5d0
         enddo
         endif

      endif

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
      if((.not. doSubprocesses) .or. (sub_number .eq. 2)) then

!         print*, 'borncache'
!         print*, lokt(4), borncache(5,1), borncache(7,1), borncache(6,1), borncache(1,1)
!         print*, lokt(8), borncache(4,2), borncache(3,2)
!         print*, lokt(5), borncache(5,3), borncache(7,3), borncache(6,3), borncache(1,3) 
!         print*, lokt(9), borncache(4,4), borncache(3,4)
      if (lokt(4).and.fincolcalc) then !this is one of the initial-initial kinematics, e.g. eq 5.147 in Catani-Seymour
!                                                                                       first splits
      call getyourscalesready_dibos(xi,pdfchange,pdf,pdfckm,4)
      q_sf=dsqrt(mufsq(1,4))

         CALL finitecollqqREMWH(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,-2,fincol(1),signpdf1)  ! ub d
!                                                  =p(0,2,1)
         CALL finitecollqqREMWH(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,1,fincol(2),signpdf1)   ! d  ub
!                                                  =p(0,2,1)
         CALL finitecollqqREMWH(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,-4,fincol(7),signpdf1)  ! cb s
!                                                  =p(0,2,1)
         CALL finitecollqqREMWH(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,3,fincol(8),signpdf1)   ! s  cb
!                                                  =p(0,2,1)
         CALL finitecollqgREMWH(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,-2,fincol(3),signpdf1)  ! ub g
!                                                  =p(0,2,1)
         CALL finitecollgqREMWH(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,fincol(4),signpdf1)     ! g  ub 
!                                                  =p(0,2,1)
         CALL finitecollqgREMWH(xi(1),xuz(1,4),p(0,1,4),p(0,2,4),p(0,3,4),q_sf,-4,fincol(9),signpdf1)  ! cb g
!                                                  =p(0,2,1)
      ! print*, 'xi callpart', xi(1), xuz(1,4)
      fincol(10)=fincol(4) 
      ! included bbar, same in lokt(5)
      ! let quarks count till 5 in finitecollgqREM, q in fincolgq is observer quark
      fincolampsq(1)=borncache(5,1) 
      fincolampsq(3)=borncache(7,1) 
      fincolampsq(2)=borncache(6,1)
      fincolampsq(4)=borncache(1,1)

      ! in finite collinears
      ! ckm effects on all quark pdfs in m2s_qqWHjj
      ! observer quarks only in finitecollgq (there use nfl to do 4flav/5flav)
      ! every other quark is coupled to w
      amp=(fincol(1)*pdf( 1,2)+fincol( 7)*pdf( 3,2))   *fincolampsq(1)
     &   +(fincol(2)*pdf(-2,2)+fincol( 8)*pdf(-4,2))   *fincolampsq(2)
     &   +(fincol(3)*pdf( 0,2)+fincol( 9)*pdf( 0,2))   *fincolampsq(3)
     &   +(fincol(4)*pdf(-2,2)+fincol(10)*pdf(-4,2))   *fincolampsq(4)

      !print*, fincolampsq(1:4)
      !print*, fincol(1:10)

      coll3=2d0*amp*fincolconv*alfas/xuz(1,4)
      endif


!      
      if (lokt(5).and.fincolcalc) then  ! this is one of the initial-initial kinematics, e.g. eq 5.147 in Catani-Seymour
!
      call getyourscalesready_dibos(xi,pdfchange,pdf,pdfckm,5)
      q_sf=dsqrt(mufsq(1,5))
!                                                                                  first splits
      CALL finitecollqqREMWH(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,5),q_sf,-2,fincol(1),signpdf2)  ! ub d
!                                                  =p(0,2,1)
      CALL finitecollqqREMWH(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,5),q_sf,1,fincol(2),signpdf2)   ! d  ub
!                                                  =p(0,2,1)
      CALL finitecollqqREMWH(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,5),q_sf,-4,fincol(7),signpdf2)  ! cb s
!                                                  =p(0,2,1)
      CALL finitecollqqREMWH(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,5),q_sf,3,fincol(8),signpdf2)   ! s  cb
!                                                  =p(0,2,1)
      CALL finitecollqgREMWH(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,5),q_sf,-2,fincol(3),signpdf2)  ! ub g
!                                                  =p(0,2,1)
      CALL finitecollgqREMWH(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,5),q_sf,fincol(4),signpdf2)     ! g  ub 
!                                                  =p(0,2,1)
      CALL finitecollqgREMWH(xi(2),xuz(1,5),p(0,1,5),p(0,2,5),p(0,3,5),q_sf,-4,fincol(9),signpdf2)  ! cb g
!                                                  =p(0,2,1)
      fincol(10)=fincol(4)
      fincolampsq(1)=borncache(5,3)
      fincolampsq(3)=borncache(7,3)
      fincolampsq(2)=borncache(6,3)
      fincolampsq(4)=borncache(1,3)

      amp=(fincol(1)*pdf( 1,1)+fincol(7)*pdf( 3,1)) *fincolampsq(1)
     &   +(fincol(2)*pdf(-2,1)+fincol(8)*pdf(-4,1)) *fincolampsq(2)
     &   +(fincol(3)*pdf(0,1)+fincol(9)*pdf(0,1))   *fincolampsq(3)
     &   +(fincol(4)*pdf(-2,1)+fincol(10)*pdf(-4,1))*fincolampsq(4)

      coll4=amp*2d0*fincolconv*alfas/xuz(1,5)

      endif

      if (lokt(8).and.fincolcalc) then
      call getyourscalesready_dibos(xi,pdfchange,pdf,pdfckm,8)
      q_sf=dsqrt(mufsq(1,8))

      CALL finitecollqgREMWH(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,1,fincol(5),signpdf1)   ! d  g
!                                                  =p(0,2,1)
      CALL finitecollgqREMWH(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,fincol(6),signpdf1)     ! g  d 
!                                                  =p(0,2,1)
      CALL finitecollqgREMWH(xi(1),xuz(1,8),p(0,1,8),p(0,2,8),p(0,3,8),q_sf,3,fincol(11),signpdf1)  ! s  g
!                                                  =p(0,2,1)
      fincol(12)=fincol(6)
      fincolampsq(5)=borncache(4,2)
      fincolampsq(6)=borncache(3,2)

      amp=(fincol(5)*pdf(0,2)+fincol(11)*pdf(0,2))  *fincolampsq(5)!
     &   +(fincol(6)*pdf(1,2)+fincol(12)*pdf(3,2))  *fincolampsq(6)

      coll7=amp*2d0*fincolconv*alfas/xuz(1,8)
      endif


      if (lokt(9).and.fincolcalc) then
      call getyourscalesready_dibos(xi,pdfchange,pdf,pdfckm,9)
      q_sf=dsqrt(mufsq(1,9))

      CALL finitecollqgREMWH(xi(2),xuz(1,9),p(0,1,9),p(0,2,9),p(0,3,9),q_sf,1,fincol(5),signpdf2)   ! d  g
!                                                  =p(0,2,1)
      CALL finitecollgqREMWH(xi(2),xuz(1,9),p(0,1,9),p(0,2,9),p(0,3,9),q_sf,fincol(6),signpdf2)     ! g  d 
!                                                  =p(0,2,1)
      CALL finitecollqgREMWH(xi(2),xuz(1,9),p(0,1,9),p(0,2,9),p(0,3,9),q_sf,3,fincol(11),signpdf2)  ! s  g
!                                                  =p(0,2,1)
      fincol(12)=fincol(6)
      fincolampsq(5)=borncache(4,4)
      fincolampsq(6)=borncache(3,4)

      amp=(fincol(5)*pdf(0,1)+fincol(11)*pdf(0,1))  *fincolampsq(5)
     &   +(fincol(6)*pdf(1,1)+fincol(12)*pdf(3,1))  *fincolampsq(6)

      coll8=amp*2d0*fincolconv*alfas/xuz(1,9)

      endif
      endif ! sub_number 

      m2s(1) = m2s_fey
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
      m2s(0)=m2s(1)+m2s(2)+m2s(3)+m2s(4)+m2s(5)+m2s(6)+m2s(7)+m2s(8)+m2s(9)+m2s(10)+m2s(11)


!     do i=0,11
!        m2s(i)=m2s(i) !*4d0 ! random photon helicity summation
!                      ! so there is no helicity summation in the code
!                      ! before this line
!     enddo

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


      if(lprintmomenta) then
            print*, "momenta"
            print*, "n_v", n_v
            do k=1,n_v
               print*, v(0:3,k,i)
            enddo
            print*, "n_p", n_p
            do k=1,n_p
               print*, p(0:3,k,i)
            enddo
      endif

      if (ldebug) then
         if ( abs(dotrr(p(0,2,1),p(0,4,1))).le.1e-1 
     &    .or.abs(dotrr(p(0,1,1),p(0,4,1))).le.1e-1
     &    .or.abs(dotrr(p(0,2,1),p(0,3,1))).le.1e-1
     &    .or.abs(dotrr(p(0,1,1),p(0,3,1))).le.1e-1
     &    .or.abs(dotrr(p(0,3,1),p(0,4,1))).le.1e-1) then
!       if (abs(m2s(1)/(m2s(10)+m2s(11))).ge.0.9.and.abs(m2s(1)/(m2s(3)+m2s(5))).le.1.1) then
!      print*, "=============================="
!      print*, "p1p2=",dotrr(p(0,1,1),p(0,2,1))
!      print*, "=============================="
!      print*, "p1v3=",dotrr(p(0,1,1),v(0,3,1))
!      print*, "p2v3=",dotrr(p(0,2,1),v(0,3,1))
!      print*, "p3v3=",dotrr(p(0,3,1),v(0,3,1))
!      print*, "p4v3=",dotrr(p(0,4,1),v(0,3,1))
         print*, "=============================="

         print*, "p2p4=",dotrr(p(0,2,1),p(0,4,1))
         print*, "p1p4=",dotrr(p(0,1,1),p(0,4,1))
         print*, "p2p3=",dotrr(p(0,2,1),p(0,3,1))
         print*, "p1p3=",dotrr(p(0,1,1),p(0,3,1))
         print*, "p3p4=",dotrr(p(0,3,1),p(0,4,1))
         print*, "=============================="

         print*, "sub_number            ", sub_number
         print*, "lincludebbfinaldo     ", lincludebbfinaldo
         print*, "lincludebb            ", lincludebb
!         print*, "bbfinaldo  1,10,11    ", lincludebbfinaldo(1), lincludebbfinaldo(10), lincludebbfinaldo(11)
         print*, "coll3                 ", coll3
         print*, "coll4                 ", coll4
         print*, "coll7                 ", coll7
         print*, "coll8                 ", coll8
         print*, "higgs propagator fac  ", fac , "not included"
         print*, "drive integration via ", m2s(0), " Fract ", m2s(0)/m2s(1)
         print*, "real emission ME      ", m2s(1) , lokt(1)
         print*, "Dipole 1              ", kinem(1) , lokt(2)
         print*, "Dipole 2              ", kinem(2) , lokt(3)
         print*, "Dipole 3              ", kinem(3) , lokt(4)
         print*, "Dipole 4              ", kinem(4) , lokt(5)
         print*, "Dipole 5              ", kinem(5) , lokt(6)
         print*, "Dipole 6              ", kinem(6) , lokt(7)
         print*, "Dipole 7              ", kinem(7) , lokt(8)
         print*, "Dipole 8              ", kinem(8) , lokt(9)
         print*, "Dipole 9              ", kinem( 9) , lokt(10)
         print*, "Dipole 10             ", kinem(10) , lokt(11)
         print*, "REM/Dipole 1+3        ", m2s(1)/(kinem(1)+kinem(3))
         print*, "REM/Dipole 2+4        ", m2s(1)/(kinem(2)+kinem(4))
         print*, "REM/Dipole 5+7        ", m2s(1)/(kinem(5)+kinem(7))
         print*, "REM/Dipole 6+8        ", m2s(1)/(kinem(6)+kinem(8))
         print*, "REM/Dipole 9+10       ", m2s(1)/(kinem(9)+kinem(10))
         print*, "REM/Dipole tot        ", m2s(1)/((kinem(1)+kinem(2)+
     &   kinem(3)+kinem(4)+kinem(5)+kinem(6)+kinem(7)+kinem(8)+kinem(9)+kinem(10)))
         print*, "---------------------------------"
         read(*,*)
         endif
      endif !ldebug

! include higgs propagator factor
      ! factor independent of kinematic
      m2s=m2s*fac

      end


      subroutine getyourscalesready_dibos(xi,pdfchange,pdf,pdfckm,l)

      implicit none

      include 'VBFNLO/utilities/global.inc'
      include 'VBFNLO/utilities/scales.inc'
      include 'VBFNLO/utilities/coupl.inc'

      double precision q_sf
      real*8, intent(out) :: pdf(-6:6,2),pdfckm(-6:6,2)
      real*8 pdftmp(-6:6,2)
      real*8, intent(in) :: xi(nx)
      integer k,pdfsign
      logical pdfchange
      integer, intent(in) :: l

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
      call pdfckmtransform(pdf(-6,1),pdfckm(-6,1))


! running coupling
      alfas  =  als(1,l) 
      G      =  dsqrt(4d0*pi*alfas)
      GG(1)  =  -G
      GG(2)  =  -G

!      print*, "scales routine ",xi, xi2

      end


      subroutine finitecollqqREMWH(xa,x,patilde,pb,pin,mufac,quark,result,pdfchange)

      implicit none

#include "VBFNLO/utilities/global.inc"

      integer quark
      integer pdfchange

      double precision xa,x,patilde(0:3),pb(0:3),pin(0:3),result,dotrr,
     &   xxa,pdfqxa,pdfgxa,pdfqxxa,logomx, logomxxa,logx,omx,omxxa,
     &   Logmupapb2,logmupapi2,NF,mufacsq,mufac,pdftmp(-6:6),one,inv2
      double precision pdfckm(-6:6)

      parameter (one=1d0,inv2=1d0/2d0)
      external dotrr

      nf=dble(nfl)

      mufacsq=mufac**2
      xxa=x*xa
      omxxa=one-xxa
      omx=one-x
!      print*, xa,xxa, " in rems"

      logx=Log(x)
!       if (x.ge.0.999995) logx=0d0

      logomx=Log(omx)
!       if (x.le.0.000005) logomx=0d0

      logomxxa=Log(omxxa)
!       if (xxa.le.0.000005) logomxxa=0d0

      
      logmupapb2=Log(mufacsq*inv2/dotrr(patilde,pb))
      logmupapi2=Log(mufacsq*inv2/dotrr(patilde,pin))

      ! no change here to include bbar
      ! change here to include ckm
      call pdfproton( xa, mufac, pdftmp(-6) )
      call pdfckmtransformsingle(pdftmp,pdfckm)
      pdfqxa=pdfckm(pdfchange*quark)/xa !pdf before splitting !ckm!
      pdfgxa=pdftmp(0)/xa
      call pdfproton( xxa, mufac, pdftmp(-6) )   
      call pdfckmtransformsingle(pdftmp,pdfckm)
      pdfqxxa=pdfckm(pdfchange*quark)/xxa  !pdf after splitting !ckm!
!      print*, pdfqxxa, xxa

      ! first pdf is before splitting and only shows up in q -> q g
      ! second pdf is after splitting and calculation shows it also only
      ! shows up in q -> q g and not in g -> q q
      ! therefore ckm is included in both pdfs

!      pdfgxxa=pdftmp(0)/xxa !not used here
  
!      print*, 'in fincolqq'
!      print*, x, xa, xxa, 'pdfs', pdfqxa, pdfgxa, pdfqxxa, 'logs',  logmupapb2, logmupapi2, logomx, logomxxa

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


      subroutine finitecollgqREMWH(xa,x,patilde,pb,pin,mufac,result,pdfchange)

      implicit none

#include "VBFNLO/utilities/global.inc"

      integer quark
      integer pdfchange

      double precision xa,x,patilde(0:3),pb(0:3),pin(0:3),result,dotrr,
     &   xxa,pdfqxa,pdfgxa,pdfqxxa,pdfgxxa,logomx, logomxxa,logx,omx,omxxa,
     &   Logmupapb2,logmupapi2,NF,mufacsq,mufac,pdftmp(-6:6),one,inv2,zero
      parameter (one=1d0,inv2=1d0/2d0,zero=0d0)
      integer nfext
      external dotrr

      !nfl: global variable for number of quark flavors to include
      ! nf: double version of nfl which shows up in p,k terms of catani/seymour
      !nfext: number of external flavors used, only 5 if nfl=5 and flag set to include external bbbar
      nf  = dble(nfl)
      if (nfl.eq.5 .and. lincludebb) then
         nfext = 5
      else if(nfl.lt.5) then
         if (.not.lincludebb) then
            nfext = 4
         else 
            print*, 'input file inconsistent: including b bbar but not using 5flavor scheme'
            stop
         endif
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

      ! include bbar by nf 4->5
      ! no change here to include ckm
      do quark=1,nfext  ! sum over quark gluon splittings
      pdfqxa=pdfqxa+pdftmp(pdfchange*quark)/xa
      enddo
      do quark=-nfext,-1
      pdfqxa=pdfqxa+pdftmp(pdfchange*quark)/xa
      enddo

      call pdfproton( xxa, mufac, pdftmp(-6) )    
      pdfgxxa=pdftmp(0)/xxa

      !not used here:
!     do quark=1,nfext  ! sum over quark gluon splittings
!     pdfqxxa=pdfqxxa+pdftmp(pdfchange*quark)/xxa
!     enddo
!     do quark=-nfext,-1
!     pdfqxxa=pdfqxxa+pdftmp(pdfchange*quark)/xxa
!     enddo

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


      subroutine finitecollqgREMWH(xa,x,patilde,pb,pin,mufac,quark,result,pdfchange)

      implicit none

#include "VBFNLO/utilities/global.inc"


      integer quark
      integer pdfchange

      double precision xa,x,patilde(0:3),pb(0:3),pin(0:3),result,dotrr,
     &   xxa,pdfqxa,pdfgxa,pdfqxxa,logomx, logomxxa,logx,omx,omxxa,
     &   Logmupapb2,logmupapi2,NF,mufacsq,mufac,pdftmp(-6:6),one,inv2,zero
      parameter (one=1d0,inv2=1d0/2d0,zero=0d0)
      external dotrr

      nf=dble(nfl)

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

      ! no change here to include ckm
      ! no change here to include bbar
      call pdfproton( xa, mufac, pdftmp(-6) )
      pdfgxa=pdftmp(0)/xa
      pdfqxa=pdftmp(pdfchange*quark)/xa
         
      call pdfproton( xxa, mufac, pdftmp(-6) )    
!      pdfgxxa=pdftmp(0)/xxa !not used here
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


! THIS SUBROUTINE CALCULATES THE DIPOLES ACCORDING TO CATANI-SEYMOUR

      subroutine daisjdipole67_dib(n,sel,kin,xuz,pmi,pk,pa,p,v,dipoleres)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"


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
      double precision dipole(0:3,0:3),p(0:3,3),
     &                 v(0:3,max_v),pipk,pkpa,pipa,xuz(2),dipoleres(3),
     &                 bornmom(0:3,6),color,symmetry
c
      double complex borncurr(0:3)
      double precision Contract_Trjcj,borntest1,borntest2,metric(0:3,0:3)

      double precision xika,ui
      double precision dotrr
      external dotrr,Contract_Trjcj

      do mu=1,3
         dipoleres(mu) = 0
      enddo
! Loops_sub_NLO handeled in calling part      
!     if(Loops_sub_NLO.ne.1) then
c     kin1 -> both, kin2 -> g
!         if((sub_number.eq.2).and.(kin.eq.2)) return
!      endif
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
!         bornmom(mu,6)=v(mu,3)
      enddo

      if (sel.eq.1) then
         call ugdHWj_curr(bornmom,borncurr,n,3)
         color=4d0      !born matrix element colorfactor =CF*Trace(1_3)
     &     *(-1d0)/2d0 !Dipole colorfactor = -1/2 CA/CA 
         symmetry=1d0/96d0  ! 1/3*1/8 (inital state colors) * 1/2 *1/2 (initial state helicities)

         if (ldebug) then
            CAll UGDHWj_msq(bornmom,borntest1,n,3) !u~ g > d a ve
            borntest2=-((contract_Trjcj(metric,borncurr)
     &           )*4./96.)
            if (abs(borntest1-borntest2).ge.1e-10) then
               Print*, "error in matrixelement calculation 67:", (BORNTEST1-BORNTEST2)/BORNTEST1*100d0, "%"
               STOP
            endif
         endif

      elseif (sel.eq.2) then
         Call dguHWj_curr(bornmom,borncurr,n,3)
         color=4d0      !born matrix element colorfactor =CF*Trace(1_3)
     &     *(-1d0)/2d0 !Dipole colorfactor = -1/2 CA/CA 
         symmetry=1d0/96d0  ! 1/3*1/8 (inital state colors) * 1/2 *1/2 (initial state helicities)

         if (ldebug) then      
            CAll DGUHWj_msq(bornmom,borntest1,n,3)
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


      subroutine daibdipole147_dib(n,sel,kin,xuz,pmi,pa,pb,p,v,dipoleres)

         implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"

         integer mu,nu,sel,kin,n
c INPUT
c
         double precision pmi(0:3),pb(0:3),pa(0:3),p(0:3,3),v(0:3,max_v)
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

         double complex borncurr(0:3)
         double precision Contract_Trjcj,borntest1,borntest2,metric(0:3,0:3)

c
         double precision xiab,ui
         double precision dotrr
         external dotrr,Contract_Trjcj

         do mu=1,3
            dipoleres(mu) = 0
         enddo
! Loops_sub_NLO handeled in calling part         
!         if(Loops_sub_NLO.ne.1) then
!c     kin1 -> both, kin2 -> g
!            if((sub_number.eq.2).and.(kin.eq.2)) return
!         endif
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
!            bornmom(mu,6)=v(mu,3)   !h
         enddo

         if (sel.eq.1) then
            Call ugdHWj_curr(bornmom,borncurr,n,3)
            color=4d0      !born matrix element colorfactor =CF*Trace(1_3)
     &     *(-1d0)/2d0 !Dipole colorfactor = -1/2 CA/CA 
            symmetry=1d0/96d0  ! 1/2*1/8 (inital state colors) * 1/2 *1/2 (initial state helicities)

            if (ldebug) then
               CAll UGDHWj_msq(bornmom,borntest1,n,3)
               borntest2=-((contract_Trjcj(metric,borncurr)
     &           )*4./96.)
!               print*, "m2s_qqWHjj borntest ", borntest1, borntest2
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
            Call dguHWj_curr(bornmom,borncurr,n,3)
            color=4d0      !born matrix element colorfactor =CF*Trace(1_3)
     &     *(-1d0)/2d0 !Dipole colorfactor = -1/2 CA/CA 
            symmetry=1d0/96d0  ! 1/2*1/8 (inital state colors) * 1/2 *1/2 (initial state helicities)


            if (ldebug) then
               CAll DGUHWj_msq(bornmom,borntest1,n,3)
               borntest2=-((contract_Trjcj(metric,borncurr)
     &            )*4./96.)
!               print*, "borntest ", borntest1, borntest2
               if (abs(borntest1-borntest2).ge.1e-10) then
                  Print*, "error in matrixelement calculation 147:", BORNTEST1, BORNTEST2
                  STOP
               endif
            endif

         endif

!            if (kin.eq.1) then
!               bornmom(mu,1)=p(mu,2)!!!
!               bornmom(mu,2)=p(mu,1)
!            elseif (kin.eq.2) then
!               bornmom(mu,1)=p(mu,1)!!!
!              bornmom(mu,2)=p(mu,2)
!            endif
!            bornmom(mu,3)=p(mu,3)   !d~
!            bornmom(mu,4)=v(mu,1)   !e-
!            bornmom(mu,5)=v(mu,2)   !ve~
!            bornmom(mu,6)=v(mu,3)   !h
!            if (ldebug) then
!               CAll SUD_wmhgL(bornmom,borntest1)
!               borntest2=-((contract_Trjcj(metric,borncurr(0,-1))+contract_Trjcj(metric,borncurr(0,1)))*4./36.)
!               if (abs(borntest1-borntest2).ge.1e-7) then
!                  Print*, "zz error in matrixelement calculation 147:", (BORNTEST1-BORNTEST2)/BORNTEST1*100d0, "%"
!                  Print*, borntest1,borntest2
!!                  STOP
!               endif
!            endif
!            if (ldebug) then
!               CAll SUD_wphgL(bornmom,borntest1)
!               borntest2=-((contract_Trjcj(metric,borncurr(0,-1))+contract_Trjcj(metric,borncurr(0,1)))*4./36.)
!               if (abs(borntest1-borntest2).ge.1e-7) then
!                  Print*, "xx error in matrixelement calculation 147:", (BORNTEST1-BORNTEST2)/BORNTEST1*100d0, "%"
!                  Print*, borntest1,borntest2
!!                  STOP
!               endif
!            endif

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


      subroutine daisjdipole65_dib(n,sel,kin,xuz,pmi,pk,pa,p,v,dipoleres)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"

      integer mu,sel,kin,n
c INPUT

      double precision pmi(0:3),pk(0:3),pa(0:3)
      double precision pis,CF,CA,Tr
      parameter (pis=3.141592653589793238d0,CF=4d0/3d0,CA=3d0,Tr=1d0/2d0)
      logical ldebug
      parameter (ldebug=.false.)

c OUTPUT
c
      double precision dipole,p(0:3,3),v(0:3,max_v),
     &    pipk,pkpa,pipa,xuz(2),dipoleres(6),bornmom(0:3,6),color(3)

      double precision born
      double precision Contract_Trjcj

      double precision xika,ui
      double precision dotrr
      external dotrr,Contract_Trjcj

      do mu=1,6
         dipoleres(mu) = 0
      enddo
! Loops_sub_NLO handeled in calling part      
 !     if((Loops_sub_NLO.ne.1).and.(sub_number.eq.2)) return
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
!         bornmom(mu,6)=v(mu,3)  ! h
      enddo

      color(1)=-CA/2./CF
      color(2)=(CA-2.*CF)/2./CF
      color(3)=0.

      if (sel.eq.1) then
         Call UDGHWj_msq(bornmom,born,n,3)
      elseif (sel.eq.2) then
         call UGDHWj_msq(bornmom,born,n,3)
      elseif (sel.eq.3) then
         Call DGUHWj_msq(bornmom,born,n,3)
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


      subroutine daibdipole145_dib(n,sel,kin,xuz,pmi,pa,pb,p,v,dipoleres)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"

      integer mu,sel,kin,n
c INPUT
c
      double precision pmi(0:3),pb(0:3),pa(0:3),p(0:3,3),v(0:3,max_v)
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

      do mu=1,6
         dipoleres(mu) = 0
      enddo
! Loops_sub_NLO handeled in calling part      
!      if((Loops_sub_NLO.ne.1).and.(sub_number.eq.2)) return
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
!         bornmom(mu,6)=v(mu,3)   !h
      enddo

      color(1)=(CA-2d0*CF)/2d0/CF
      color(2)=-CA/2d0/CF
      color(3)=0d0

      if (sel.eq.1) then
         Call UDGHWj_msq(bornmom,born,n,3)
      elseif (sel.eq.2) then
         Call UGDHWj_msq(bornmom,born,n,3)
      elseif (sel.eq.3) then
         Call DGUHWj_msq(bornmom,born,n,3)
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


      subroutine dasijdipole40_dib(n,sel,kin,xuz,pmi,pj,pa,p,v,dipoleres)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"

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
      double precision dipole(0:3,0:3),p(0:3,3),v(0:3,max_v)
     &    ,pipj,pjpa,pipa,xuz(2),dipoleres(3),bornmom(0:3,6),color,symmetry,borntmp,
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
! Loops_sub_NLO handeled in calling part      
!      if(Loops_sub_NLO.ne.1) then
!         if((sub_number.eq.1).and.(kin.ne.1)) return
!         if((sub_number.eq.2).and.(kin.eq.1)) return
!      endif
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
            read*
         endif
         if (abs(zi-pipa/(pipa+pjpa)).ge.1e-10) then
            PRINT*, "error related to zi in dipole 40"
            read*
         endif
!         print*, "kin 1 ptilde_a(",mu," ) check in 40"
!         read*,
         do mu=0,3
!           if (sel.eq.2) then
!              if (abs(p(mu,2)-xija*pa(mu)).ge.1e-10) then
!                 print*, "kin 1 wrong ptilde_a(",mu," ) in 40", p(mu,2),xija*pa(mu)
!                 read*
!              endif
!           elseif (sel.eq.1) then
!              if (abs(p(mu,1)-xija*pa(mu)).ge.1e-10) then
!                 print*, "kin 1 wrong ptilde_a(",mu," ) in 40", p(mu,2),xija*pa(mu)
!                 read*
!              endif
!           endif
            if (abs(p(mu,3)-(pj(mu)+pmi(mu)-(1d0-xija)*pa(mu))).ge.1e-10) then
               print*, "wrong ptilde_ij(",mu," ) in 40"
               read*
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
!         bornmom(mu,6)=v(mu,3)
      enddo

      if (kin.eq.1.or.kin.eq.4) then
         call udgHWj_curr(bornmom,borncurr,n,3)
         color=4d0      !born matrix element colorfactor =CF*Trace(1_3)
     &     *(-1d0)/2d0 !Dipole colorfactor = -1/2 CA/CA 
         symmetry=1d0/36d0  ! 1/3*1/3 (inital state colors) * 1/2 *1/2 (initial state helicities)

      elseif (kin.eq.2) then
         CALL UGDHWj_msq(bornmom,borntmp,n,3)

      elseif (kin.eq.3) then
         CALL DGUHWj_msq(bornmom,borntmp,n,3)
      endif

!       if (kin.eq.1 .or. kin.eq.4 .and. ldebug) then
!          do mu=0,3
!             do nu=0,mu-1
!                metric(mu,nu)=0d0
!                metric(nu,mu)=0d0
!             enddo
!          enddo
!          metric(0,0)=1d0
!          metric(1,1)=-1d0
!          metric(2,2)=-1d0
!          metric(3,3)=-1d0
!          !CAll SUD_wphgL(bornmom,testsign,borntest1)
!          CAll UDGHWj_msq(bornmom,borntest1,n,36)
!          borntest2=-(contract_Trjcj(metric,borncurr))*4./36.
!          print*, "trying born test", borntest1, borntest2
!          if (abs(borntest1-borntest2).ge.1e-7) then
!             Print*, "error in matrixelement calculation 40:", (BORNTEST1-BORNTEST2)/BORNTEST1*100d0, "%"
!             Print*, borntest1/borntest2
!             read*
!          endif
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


      subroutine test40_dib(n,sel,kin,xuz,pj,pmi,pa,p,v,dipoleres)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"

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
      double precision p(0:3,3),v(0:3,max_v),
     &    pipj,pjpa,pipa,xuz(2),dipoleres(3),bornmom(0:3,7),color,borntmp,
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
! Loops_sub_NLO handeled in calling part      
!      if((Loops_sub_NLO.ne.1).and.(sub_number.eq.2)) return
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
!         bornmom(mu,6)=v(mu,3)
      enddo

      dptmp=-1./2./xija*8.*pis*alfas*4./3.*(2./(1.-zi+1.-xija)-1.-zi)!/pipj


      if (kin.eq.2) then
         CALL UGDHWj_msq(bornmom,borntmp,n,3)
      elseif (kin.eq.3) then
         CALL DGUHWj_msq(bornmom,borntmp,n,3)
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
      CAll ugdHWj_curr(bornmom,borncurr,n,3)
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



