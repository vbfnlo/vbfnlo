      SUBROUTINE mgCheck()

      implicit none
#include "ggfloMG.inc"
#include "ggf_process.inc"

      integer i,j,MGfsign(6)
      real*8 MGrh!(0:3)
      integer isig,ipol,NHEL(6),isig1,isig3,i1pol,i2pol,i3pol,isig2,isig4
      complex*16 m1kb(-1:1,12),m2kb(-1:1,12),m1kb2(12),m2kb2(12)
      real*8 ampHel1,ampHel2,ampHel3,ANS
      integer aMG,bMG,cMG,dMG,eMG,fMG
      REAL*8 ALFAS,XMT,ALFA,XMZ,XMW,SIN2W,XMH
      COMMON /BKOPIN/ ALFAS,XMT,ALFA,XMZ,XMW,SIN2W,XMH
c... MadGraph variables
      double precision G,tau,series_t,series_p,scalarf,axialf,PI,V,WM

      REAL*8 MATRIXucucgH2,MATRIXuuuugH2,MATRIXAUD,MATRIXuguggH2,
     &     MATRIXHgg,MATRIXgggggH2,MATRIXgggggA2,MATRIXuuxccxgH2,
     &     MATRIXuuxuuxgH2,MATRIXguuuuxH2,MATRIXguuccxH2
      EXTERNAL MATRIXucucgH2,MATRIXuuuugH2,MATRIXAUD,MATRIXuguggH2,
     &     MATRIXuuxccxgH2,MATRIXuuxuuxgH2,MATRIXguuuuxH2,
     &     MATRIXguuccxH2
cfc New MG5 subroutine
      Real*8 MATRIXuuuugH3,MATRIXucucgH3,MATRIXubuHubug3
      External MATRIXuuuugH3,MATRIXucucgH3,MATRIXubuHubug3
      Real*8 MATRIXubuHdbdg3,MATRIXguuuuH3,MATRIXguuddH3
      External MATRIXubuHdbdg3,MATRIXguuuuH3,MATRIXguuddH3
cfc
      common/colorna/colorN
      complex*16 colorN(24),colorM(24)

      ampHel1=0.0d0
      ampHel2=0.0d0
      ampHel3=0.0d0

c... couplings from koppeln.F (old VBFNLO) needed for coupl.inc
c... needed in MadGraph subroutines
c$$$      gh(1)=(5.08415494103491686d-005,0.0d0)
c$$$      gh(2)=(0.0d0,0.0d0)
c$$$      ga(1)=(0.0d0,0.0d0)
c$$$      ga(2)=(7.62623241155237496d-005,0.0d0)

c... from MG Source/MODELS/couplings.f
      PI=3.14159265358979323846d0
      WM=sqrt(XMZ**2/2d0+
     &     sqrt(XMZ**4/4d0-Pi/(1.414213562d0)*alfa/1.16639d-5*XMZ**2))
      V= 2d0*WM*DSQRT(SIN2W)/DSQRT(ALFA*4d0*PI) 
      G = DSQRT(4d0*PI*ALFAS)
      scalarf = 1d0
      axialf  = 0d0
      tau = XMH**2/(4d0*m_top**2)
cfc I have comment in the values of series_t and series_p
cfc      series_t = 1d0
cfc      series_p = 1d0
      series_t = 1d0 + tau*7d0/30d0 + tau**2*2d0/21d0 + tau**3*26d0/525d0
      series_p = 1d0 + tau/3d0      + tau**2*8d0/45d0 + tau**3*4d0/35d0
      gh(1) = dcmplx( scalarf*g**2/4d0/PI/(3d0*PI*V)*series_t, 0d0)
      gh(2) = dcmplx( axialf *g**2/4d0/PI/(2d0*PI*V)*series_p, 0d0)
c Pseudo-scalar Higgs coupling:
      ga(1) = dcmplx( 0d0                  , 0d0)
      ga(2) = dcmplx( g**2/4d0/PI/(2d0*PI*V)*series_p, 0d0)
c   Couplings for the 4-vertices VVVS
      gh4(1) = G*gh(1)
      gh4(2) = G*gh(2)
      ga4(1) = G*ga(1)
      ga4(2) = G*ga(2)
c      print*,'gh4(1) ',gh4(1) 
c      print*,'gh4(2) ',gh4(2) 
c      print*,'ga4(1) ',ga4(1) 
c      print*,'ga4(2) ',ga4(2) 


c      do i=1,6
c         print*,'MGpbar',i,'=',(MGpbar(j,i),j=0,3)
c      enddo

*---------------------------------------------------------------
      if (subprQQ) then
         MGfsign(1)=1
         MGfsign(2)=1 
         MGfsign(3)=1
         MGfsign(4)=1
         MGfsign(5)=1
         MGfsign(6)=1
c         NHEL(1)=-1
c         NHEL(2)=-1
c         NHEL(3)=-1
c         NHEL(4)=-1
         NHEL(5)=-1
         ipol =1
         isig1=1
         isig3=1
         isig2=1
         isig4=1

cmk      do ipol =1,2
cmk      do isig1=-1,1,2
cmk      do isig3=-1,1,2
cmk      do isig2=-1,1,2
cmk      do isig4=-1,1,2
c      do i=1,12     !for u d 
*--------- q q --> q q g H ---------
         NHEL(1)=isig1
         NHEL(2)=isig2 !isig3
         NHEL(3)=isig3 !isig1
         NHEL(4)=isig4 !isig3
         if (ipol.eq.1) then
            NHEL(5)=-1
         else
            NHEL(5)= 1
         endif
         
c... for check comment out the 2nd,3rd,4th Higgs3jQQ subroutine in m2s_ggf3jh.F       
c... check: uu --> uugH & uc --> ucgH
         if(mguu) then
            if(equal_f) then
               If(MGdebug1) then
               print*,'qq --> qqgH'
               endif
               ampHel1=ampHel1+MATRIXuuuugH2(MGpbar,NHEL,MGfsign) 
cfc New MG5 subroutine. Comment line or not
cmk               ampHel1=ampHel1+MATRIXuuuugH3(MGpbar,NHEL,MGfsign) 
cmk               print*,'ampHel1',ampHel1
            else
               If(MGdebug1) then
               print*,'qQ --> qQgH'
               endif
               ampHel1=ampHel1+MATRIXucucgH2(MGpbar(0,1),NHEL,MGfsign) 
cfc New MG5 subroutine.  Comment line or not
cmk               ampHel1=ampHel1+MATRIXucucgH3(MGpbar(0,1),NHEL,MGfsign) 
cmk               print*,'ampHel1',ampHel1
            endif
         endif
cfc check: uub :> u ub g H & u db > u db g H
         if(mguub) then
            MGfsign(1)=1 
            MGfsign(2)=1  
            MGfsign(3)=-1 
            MGfsign(4)=-1 
            MGfsign(5)=1
            if(equal_f) then
               If(MGdebug1) then
               print*,'qq~ --> qq~gH'
               endif
               ampHel1=ampHel1+MATRIXuuxuuxgH2(MGpbar,NHEL,MGfsign) 
cfc New MG5 subroutine. Comment line or not
cmk            ampHel1=ampHel1+MATRIXubuHubug3(MGpbar,NHEL,MGfsign) 
cfc               print*,'ampHel1',ampHel1
            else
               If(MGdebug1) then
               print*,'qq~ --> QQ~gH'
               endif
               ampHel1=ampHel1+MATRIXuuxccxgH2(MGpbar,NHEL,MGfsign) 
cfc New MG5 subroutine.  Comment line or not
cmk            ampHel1=ampHel1+MATRIXubuHdbdg3(MGpbar,NHEL,MGfsign) 
cfc               print*,'ampHel1',ampHel1
            endif
         endif
cfc check: gu > u uub H && g u > u ddb H
         if(mggu) then

            if (ipol.eq.1) then
               NHEL(1)=-1
            else
               NHEL(1)= 1
            endif
            NHEL(2)=isig1
            NHEL(3)=isig2       !isig3
            NHEL(4)=isig3
            NHEL(5)=isig4       !isig3
 
            MGfsign(1)=-1 
            MGfsign(2)=1
            MGfsign(3)=1 !1
            MGfsign(4)=1 !1
            MGfsign(5)=-1

            if(equal_f) then
               If(MGdebug1) then
               print*,'gq --> qqq~H'
               endif
               ampHel1=ampHel1+ MATRIXguuuuxH2(MGpbar,NHEL,MGfsign)
cfc New MG5 subroutine. Comment line or not
cmk               ampHel1=ampHel1+ MATRIXguuuuH3(MGpbar,NHEL,MGfsign)
cmk            do i=1,24        !for u u 
cmk               print*,'GGFLOamp1',GGFLOamp1(isig1,isig3,ipol,i)
cmk               print*,'MGamp1   ',MGamp1(NHEL(1),NHEL(2),NHEL(5),i)
cmk               print*,'ratio',i,'',abs(abs(GGFLOamp1(isig1,isig3,ipol,i))
cmk  &                 /abs(MGamp1(NHEL(1),NHEL(2),NHEL(5),i))-1d0)
cmk            enddo
cmk            print*,''
cmk            print*,'ampHel1',ampHel1
            else
               If(MGdebug1) then
               print*,'gq --> qQQ~H'
               endif
               ampHel1=ampHel1+MATRIXguuccxH2(MGpbar,NHEL,MGfsign) 
cfc New MG5 subroutine.  Comment line or not
cmk               ampHel1=ampHel1+MATRIXguuddH3(MGpbar,NHEL,MGfsign) 
cmk             do i=12,24       !for gu
cmk                print*,'GGFLOamp1',GGFLOamp1(isig1,isig3,ipol,i)
cmk                print*,'MGamp1   ',MGamp1(NHEL(1),NHEL(2),NHEL(5),i)
cmk                print*,'ratio',i,'',abs(abs(GGFLOamp1(isig1,isig3,ipol,i))
cmk   &                 /abs(MGamp1(NHEL(1),NHEL(2),NHEL(5),i))-1d0)
cmk             enddo
cmk             print*,'ampHel1',ampHel1
            endif
         endif

c... check: uc --> ucgA & uu --> uugA
c         ampHel1=MATRIXAUD(MGpbar(0,1),NHEL,MGfsign) 
cmk      do i=1,24              !for u u 
cmk         print*,'GGFLOamp1',GGFLOamp1(isig1,isig3,ipol,i)
cmk         print*,'MGamp1   ',MGamp1(NHEL(1),NHEL(2),NHEL(5),i)
cmk         print*,'ratio',i,'',abs(abs(GGFLOamp1(isig1,isig3,ipol,i))
cmk  &        /abs(MGamp1(NHEL(1),NHEL(2),NHEL(5),i))-1d0)
cmk      enddo
         If(MGdebug1) then
         print*,''
         do i=1,8
            print*,'mgJAMP    (',i,')',mgJAMP(i)
            print*,'GGFLOcolor(',i,')',GGFLOcolor(i)
            print*,'ratio',abs(abs(mgJAMP(i))/abs(GGFLOcolor(i))-1d0)
         enddo
         print*,''
         print*,'GGFLOampSq1',GGFLOampSq1
         print*,'MGampsq1   ',MGampsq1
         print*,'ratio',abs(abs(MGampsq1/GGFLOampSq1)-1d0)
         if (abs(abs(MGampsq1/GGFLOampSq1)-1d0).lt.1d-1) then
            Print*, ipol, isig1,isig3,isig2,isig4
cmk            pause
         endif
         print*,isig1,isig2,isig3,isig4,ipol
         do i=1,24 
            print*,'GGFLOamp1',GGFLOamp1(isig1,isig3,ipol,i)
            print*,'MGamp1   ',MGamp1(NHEL(1),NHEL(2),NHEL(5),i)
            print*,'ratio',i,'',abs(abs(GGFLOamp1(isig1,isig3,ipol,i))
     &           /abs(MGamp1(NHEL(1),NHEL(2),NHEL(5),i))-1d0)
         enddo
         endif
cmk      enddo
cmk      enddo
cmk      enddo
cmk      enddo
cmk      enddo
cmk         if(abs(abs(ampHel1/GGFLOampSq1)-1d0).gt.1d-9) then
cmk            print*,'GGFLOampSq1',GGFLOampSq1
cmk            print*,'ampHel1   ',ampHel1
cmk            print*,'ratio',abs(abs(ampHel1/GGFLOampSq1)-1d0)
cmk         endif
      endif
*--------- ug ---------
      if (subprQG) then

c$$$      CALL SMATRIXuguggh(MGpbar,ANS)
c$$$      print*,'ANS',ANS
c$$$      STOP
c... CHECK possible with res(8) in m2s_ggf3jh.F !!!!!!!!
      i1pol = 1
      i2pol = 1
      i3pol = 1

      MGfsign(1)=1
      MGfsign(2)=-1 
      MGfsign(3)=1
      MGfsign(4)=1
      MGfsign(5)=1
      MGfsign(6)=1
 
      NHEL(1)=-1
c      NHEL(2)=-1
      NHEL(3)=-1
c      NHEL(4)=-1
c      NHEL(5)=-1

      if (i1pol.eq.1) then
         NHEL(2)=-1
      else
         NHEL(2)= 1
      endif
      if (i2pol.eq.1) then
         NHEL(4)=-1
      else
         NHEL(4)= 1
      endif
      if (i3pol.eq.1) then
         NHEL(5)=-1
      else
         NHEL(5)= 1
      endif

c... for check comment only the 2nd Higgs3jQG subroutine in m2s_ggf3jh.F     
      ampHel2=MATRIXuguggH2(MGpbar(0,1),NHEL,MGfsign)
      If(MGdebug1) then
      do i=1,84
c         if (abs(abs(GGFLOamp2(i1pol,i2pol,i3pol,i))/
c     &        abs(MGamp2(NHEL(2),NHEL(4),NHEL(5),i))-1d0).gt.1d-5) then
         print*,'GGFLOamp2',i,'=',GGFLOamp2(i1pol,i2pol,i3pol,i)
         print*,'MGamp2   ',i,'=',MGamp2(NHEL(2),NHEL(4),NHEL(5),i)
         print*,'ratio ',i,'',abs(GGFLOamp2(i1pol,i2pol,i3pol,i))/
     &             abs(MGamp2(NHEL(2),NHEL(4),NHEL(5),i))-1d0
c         print*,''
c         endif
      enddo
      do i=1,6
         print*,'mgJAMP    (',i,')',mgJAMP(i)
         print*,'GGFLOcolor(',i,')',GGFLOcolor(i)
         print*,'ratio',abs(abs(mgJAMP(i))/abs(GGFLOcolor(i))-1d0)
      enddo
      print*,''
      print*,'GGFLOampSq2',GGFLOampSq2
      print*,'MGampSq2',MGampSq2
      print*,'qg ratio2',GGFLOampSq2/MGampSq2-1d0
      STOP
      endif
      endif

*--------- gg ---------
c... for gg
c$$$         gpol(1) = 1
c$$$         gpol(2) = 1
c$$$         gpol(3) = 1
c$$$         gpol(4) = 1
c$$$         gpol(5) = 1
      if (subprGG) then

      MGfsign(1)=-1
      MGfsign(2)=-1 
      MGfsign(3)=1
      MGfsign(4)=1
      MGfsign(5)=1
      MGfsign(6)=1

      NHEL(1)=-1
      NHEL(2)=-1
      NHEL(3)=-1
      NHEL(4)=-1
      NHEL(5)=-1

      if(Heff.or.HloopTop.or.HloopBottom) then
c          ampHel3=MATRIXHgg(MGpbar(0,1),NHEL,MGfsign)
         ampHel3=MATRIXgggggH2(MGpbar(0,1),NHEL,MGfsign)
      endif
      if(Aeff.or.AloopTop.or.AloopBottom) then
         ampHel3=MATRIXgggggA2(MGpbar(0,1),NHEL,MGfsign)
      endif

      If(MGdebug1) then
      do i=1,390
c         if(abs(abs(GGFLOamp3(i))/abs(MGamp3(i))-1d0).gt.1d-3) then
            print*,'GGFLOamp3(',i,')',GGFLOamp3(i)
            print*,'MGamp3   (',i,')',MGamp3(i)
            print*,'ratio gg: GGFLO/MG',i,'',abs(abs(GGFLOamp3(i))
     &           /abs(MGamp3(i))-1d0)
c         endif
         enddo
      endif

      
      print*,''
      do i=1,24
         print*,'mgJAMP    (',i,')',mgJAMP(i)
         print*,'GGFLOcolor(',i,')',GGFLOcolor(i)
         print*,'ratio',abs(abs(mgJAMP(i))/abs(GGFLOcolor(i))-1d0)
      enddo

c      print*,''
      print*,'GGFLOampSq3',GGFLOampSq3
      print*,'MGampSq3   ',MGampSq3
      print*,'gg ratio3     ',abs(GGFLOampSq3/MGampSq3-1d0)
      print*,''
      
      colorM(1)=((MGamp3(101)+MGamp3(361))
     &     -(MGamp3(363)+MGamp3(389))-(MGamp3(23)+MGamp3(283))
     &     +(MGamp3(25)-MGamp3(103))+(MGamp3(285)+MGamp3(387)))*2d0
      colorM(2)= ((MGamp3(23)-MGamp3(153))
     &     +(-MGamp3(25)+MGamp3(103))+(-MGamp3(101)+MGamp3(231))
     &     +(-MGamp3(155)+MGamp3(257))-(MGamp3(233)+MGamp3(259)))*2d0
      colorM(3)=((MGamp3(23)-MGamp3(153))
     &     -(MGamp3(25)+MGamp3(129))+(-MGamp3(75)+MGamp3(205))
     &     +(MGamp3(77)+MGamp3(127))-(MGamp3(155)+MGamp3(207)))*2d0
      colorM(4)=(-(MGamp3(23)+MGamp3(283))
     &     +(MGamp3(25)+MGamp3(129))+(MGamp3(49)+MGamp3(309))
     &     +(MGamp3(51)-MGamp3(127))+(MGamp3(285)-MGamp3(311)))*2d0
      colorM(5)=((MGamp3(23)+MGamp3(283))
     &     -(MGamp3(25)+MGamp3(129))-(MGamp3(49)+MGamp3(309))
     &     +(-MGamp3(51)+MGamp3(127))+(-MGamp3(285)+MGamp3(311)))*2d0
      colorM(6)=((-MGamp3(23)+MGamp3(153))
     &     +(MGamp3(25)+MGamp3(129))+(MGamp3(75)-MGamp3(205))
     &     -(MGamp3(77)+MGamp3(127))+(MGamp3(155)+MGamp3(207)))*2d0
      colorM(7)=((-MGamp3(23)+MGamp3(153))
     &     +(MGamp3(25)-MGamp3(103))+(MGamp3(101)-MGamp3(231))
     &     +(MGamp3(155)-MGamp3(257))+(MGamp3(233)+MGamp3(259)))*2d0
      colorM(8)=((MGamp3(23)+MGamp3(283))
     &     +(-MGamp3(25)+MGamp3(103))-(MGamp3(101)+MGamp3(361))
     &     -(MGamp3(285)+MGamp3(387))+(MGamp3(363)+MGamp3(389)))*2d0
      colorM(9)=((MGamp3(49)+MGamp3(309))
     &     +(MGamp3(51)+MGamp3(77))-(MGamp3(75)+MGamp3(335))
     &     -(MGamp3(311)+MGamp3(387))+(MGamp3(337)+MGamp3(389)))*2d0
      colorM(10)=((-MGamp3(49)+MGamp3(179))
     &     +(-MGamp3(51)+MGamp3(127))+(MGamp3(101)-MGamp3(231))
     &     +(-MGamp3(103)-MGamp3(129))+(-MGamp3(181)+MGamp3(233)))*2d0
      colorM(11)=((MGamp3(49)-MGamp3(179))
     &     +(MGamp3(51)+MGamp3(77))+(-MGamp3(75)+MGamp3(205))
     &     +(MGamp3(181)+MGamp3(259))+(-MGamp3(207)-MGamp3(257)))*2d0
      colorM(12)=((-MGamp3(49)+MGamp3(179))
     &     +(-MGamp3(51)-MGamp3(77))+(MGamp3(75)-MGamp3(205))
     &     +(-MGamp3(181)-MGamp3(259))+(MGamp3(207)+MGamp3(257)))*2d0
      colorM(13)=((MGamp3(49)-MGamp3(179))
     &     +(MGamp3(51)-MGamp3(127))+(-MGamp3(101)+MGamp3(231))
     &     +(MGamp3(103)+MGamp3(129))+(MGamp3(181)-MGamp3(233)))*2d0
      colorM(14)=((-MGamp3(49)-MGamp3(309))
     &     +(-MGamp3(51)-MGamp3(77))+(MGamp3(75)+MGamp3(335))
     &     +(MGamp3(311)+MGamp3(387))+(-MGamp3(337)-MGamp3(389)))*2d0
       colorM(15)=((MGamp3(75)+MGamp3(335))
     &     +(-MGamp3(77)-MGamp3(127))+(-MGamp3(101)-MGamp3(361))
     &     +(MGamp3(103)+MGamp3(129))+(-MGamp3(337)+MGamp3(363)))*2d0
      colorM(16)=((MGamp3(101)+MGamp3(361))
     &     +(MGamp3(337)-MGamp3(363))+(MGamp3(77)+MGamp3(127))
     &     +(-MGamp3(103)-MGamp3(129))+(-MGamp3(75)-MGamp3(335)))*2d0
      colorM(17)=((MGamp3(153)+MGamp3(283))
     &     +(MGamp3(155)+MGamp3(207))+(-MGamp3(205)-MGamp3(335))
     &     +(-MGamp3(285)-MGamp3(387))+(MGamp3(337)+MGamp3(389)))*2d0
      colorM(18)=((MGamp3(153)+MGamp3(283))
     &     +(MGamp3(155)-MGamp3(257))+(-MGamp3(179)-MGamp3(309))
     &     +(MGamp3(181)+MGamp3(259))+(-MGamp3(285)+MGamp3(311)))*2d0
      colorM(19)=((-MGamp3(181)-MGamp3(259))
     &     +(-MGamp3(153)-MGamp3(283))+(-MGamp3(155)+MGamp3(257))
     &     +(MGamp3(179)+MGamp3(309))+(MGamp3(285)-MGamp3(311)))*2d0
      colorM(20)=((-MGamp3(153)-MGamp3(283))
     &     +(-MGamp3(155)-MGamp3(207))+(MGamp3(205)+MGamp3(335))
     &     +(MGamp3(285)+MGamp3(387))+(-MGamp3(337)-MGamp3(389)))*2d0
      colorM(21)=((-MGamp3(179)-MGamp3(309))
     &     +(MGamp3(181)-MGamp3(233))+(MGamp3(231)+MGamp3(361))
     &     +(MGamp3(311)+MGamp3(387))+(-MGamp3(363)-MGamp3(389)))*2d0
      colorM(22)=((MGamp3(179)+MGamp3(309))
     &     +(-MGamp3(181)+MGamp3(233))+(-MGamp3(231)-MGamp3(361))
     &     +(-MGamp3(311)-MGamp3(387))+(MGamp3(363)+MGamp3(389)))*2d0
      colorM(23)=((MGamp3(205)+MGamp3(335))
     &     +(-MGamp3(207)-MGamp3(257))+(-MGamp3(231)-MGamp3(361))
     &     +(MGamp3(233)+MGamp3(259))+(-MGamp3(337)+MGamp3(363)))*2d0
      colorM(24)=((-MGamp3(205)-MGamp3(335))
     &     +(MGamp3(207)+MGamp3(257))+(MGamp3(231)+MGamp3(361))
     &     +(-MGamp3(233)-MGamp3(259))+(MGamp3(337)-MGamp3(363)))*2d0

c... color(1)
         colorM(1)=((MGamp3(101)+MGamp3(361))
     &     -(MGamp3(363)+MGamp3(389))-(MGamp3(23)+MGamp3(283))
     &     +(MGamp3(25)-MGamp3(103))+(MGamp3(285)+MGamp3(387)))*2d0
        colorM(2)= ((MGamp3(23)-MGamp3(153))
     &     +(-MGamp3(25)+MGamp3(103))+(-MGamp3(101)+MGamp3(231))
     &     +(-MGamp3(155)+MGamp3(257))-(MGamp3(233)+MGamp3(259)))*2d0
        colorM(3)=((MGamp3(23)-MGamp3(153))
     &     -(MGamp3(25)+MGamp3(129))+(-MGamp3(75)+MGamp3(205))
     &     +(MGamp3(77)+MGamp3(127))-(MGamp3(155)+MGamp3(207)))*2d0
         colorM(4)=(-(MGamp3(23)+MGamp3(283))
     &     +(MGamp3(25)+MGamp3(129))+(MGamp3(49)+MGamp3(309))
     &     +(MGamp3(51)-MGamp3(127))+(MGamp3(285)-MGamp3(311)))*2d0
         colorM(5)=((MGamp3(23)+MGamp3(283))
     &     -(MGamp3(25)+MGamp3(129))-(MGamp3(49)+MGamp3(309))
     &     +(-MGamp3(51)+MGamp3(127))+(-MGamp3(285)+MGamp3(311)))*2d0
          colorM(6)=((-MGamp3(23)+MGamp3(153))
     &     +(MGamp3(25)+MGamp3(129))+(MGamp3(75)-MGamp3(205))
     &     -(MGamp3(77)+MGamp3(127))+(MGamp3(155)+MGamp3(207)))*2d0
         colorM(7)=((-MGamp3(23)+MGamp3(153))
     &     +(MGamp3(25)-MGamp3(103))+(MGamp3(101)-MGamp3(231))
     &     +(MGamp3(155)-MGamp3(257))+(MGamp3(233)+MGamp3(259)))*2d0
          colorM(8)=((MGamp3(23)+MGamp3(283))
     &     +(-MGamp3(25)+MGamp3(103))-(MGamp3(101)+MGamp3(361))
     &     -(MGamp3(285)+MGamp3(387))+(MGamp3(363)+MGamp3(389)))*2d0
          colorM(9)=((MGamp3(49)+MGamp3(309))
     &     +(MGamp3(51)+MGamp3(77))-(MGamp3(75)+MGamp3(335))
     &     -(MGamp3(311)+MGamp3(387))+(MGamp3(337)+MGamp3(389)))*2d0
          colorM(10)=((-MGamp3(49)+MGamp3(179))
     &     +(-MGamp3(51)+MGamp3(127))+(MGamp3(101)-MGamp3(231))
     &     +(-MGamp3(103)-MGamp3(129))+(-MGamp3(181)+MGamp3(233)))*2d0

      colorM(11)=((MGamp3(49)-MGamp3(179))
     &     +(MGamp3(51)+MGamp3(77))+(-MGamp3(75)+MGamp3(205))
     &     +(MGamp3(181)+MGamp3(259))+(-MGamp3(207)-MGamp3(257)))*2d0
      colorM(12)=((-MGamp3(49)+MGamp3(179))
     &     +(-MGamp3(51)-MGamp3(77))+(MGamp3(75)-MGamp3(205))
     &     +(-MGamp3(181)-MGamp3(259))+(MGamp3(207)+MGamp3(257)))*2d0
      colorM(13)=((MGamp3(49)-MGamp3(179))
     &     +(MGamp3(51)-MGamp3(127))+(-MGamp3(101)+MGamp3(231))
     &     +(MGamp3(103)+MGamp3(129))+(MGamp3(181)-MGamp3(233)))*2d0
      colorM(14)=((-MGamp3(49)-MGamp3(309))
     &     +(-MGamp3(51)-MGamp3(77))+(MGamp3(75)+MGamp3(335))
     &     +(MGamp3(311)+MGamp3(387))+(-MGamp3(337)-MGamp3(389)))*2d0
       colorM(15)=((MGamp3(75)+MGamp3(335))
     &     +(-MGamp3(77)-MGamp3(127))+(-MGamp3(101)-MGamp3(361))
     &     +(MGamp3(103)+MGamp3(129))+(-MGamp3(337)+MGamp3(363)))*2d0
      colorM(16)=((MGamp3(101)+MGamp3(361))
     &     +(MGamp3(337)-MGamp3(363))+(MGamp3(77)+MGamp3(127))
     &     +(-MGamp3(103)-MGamp3(129))+(-MGamp3(75)-MGamp3(335)))*2d0
      colorM(17)=((MGamp3(153)+MGamp3(283))
     &     +(MGamp3(155)+MGamp3(207))+(-MGamp3(205)-MGamp3(335))
     &     +(-MGamp3(285)-MGamp3(387))+(MGamp3(337)+MGamp3(389)))*2d0
      colorM(18)=((MGamp3(153)+MGamp3(283))
     &     +(MGamp3(155)-MGamp3(257))+(-MGamp3(179)-MGamp3(309))
     &     +(MGamp3(181)+MGamp3(259))+(-MGamp3(285)+MGamp3(311)))*2d0
      colorM(19)=((-MGamp3(181)-MGamp3(259))
     &     +(-MGamp3(153)-MGamp3(283))+(-MGamp3(155)+MGamp3(257))
     &     +(MGamp3(179)+MGamp3(309))+(MGamp3(285)-MGamp3(311)))*2d0
      colorM(20)=((-MGamp3(153)-MGamp3(283))
     &     +(-MGamp3(155)-MGamp3(207))+(MGamp3(205)+MGamp3(335))
     &     +(MGamp3(285)+MGamp3(387))+(-MGamp3(337)-MGamp3(389)))*2d0
      colorM(21)=((-MGamp3(179)-MGamp3(309))
     &     +(MGamp3(181)-MGamp3(233))+(MGamp3(231)+MGamp3(361))
     &     +(MGamp3(311)+MGamp3(387))+(-MGamp3(363)-MGamp3(389)))*2d0
      colorM(22)=((MGamp3(179)+MGamp3(309))
     &     +(-MGamp3(181)+MGamp3(233))+(-MGamp3(231)-MGamp3(361))
     &     +(-MGamp3(311)-MGamp3(387))+(MGamp3(363)+MGamp3(389)))*2d0
      colorM(23)=((MGamp3(205)+MGamp3(335))
     &     +(-MGamp3(207)-MGamp3(257))+(-MGamp3(231)-MGamp3(361))
     &     +(MGamp3(233)+MGamp3(259))+(-MGamp3(337)+MGamp3(363)))*2d0

          colorM(24)=((-MGamp3(205)-MGamp3(335))
     &     +(MGamp3(207)+MGamp3(257))+(MGamp3(231)+MGamp3(361))
     &     +(-MGamp3(233)-MGamp3(259))+(MGamp3(337)-MGamp3(363)))*2d0

c$$$      print*,'(MGamp3(101)+MGamp3(361))',(MGamp3(101)+MGamp3(361))*2d0
c$$$      print*,'-(MGamp3(363)+MGamp3(389))',-(MGamp3(363)+MGamp3(389))*2d0
c$$$      print*,'-(MGamp3(23)+MGamp3(283))',-(MGamp3(23)+MGamp3(283))*2d0
c$$$      print*,'(MGamp3(25)-MGamp3(103))',(MGamp3(25)-MGamp3(103))*2d0
c$$$      print*,'(MGamp3(285)+MGamp3(387))',(MGamp3(285)+MGamp3(387))*2d0
c$$$      print*,'color(1) pent sum',((MGamp3(101)+MGamp3(361))
c$$$     &     -(MGamp3(363)+MGamp3(389))-(MGamp3(23)+MGamp3(283))
c$$$     &     +(MGamp3(25)-MGamp3(103))+(MGamp3(285)+MGamp3(387)))*2d0
c$$$c... color(2)
c$$$      print*,'(-MGamp3(101)+MGamp3(231))',(-MGamp3(101)+MGamp3(231))*2d0
c$$$      print*,'-(MGamp3(233)+MGamp3(259))',-(MGamp3(233)+MGamp3(259))*2d0
c$$$      print*,'(MGamp3(23)-MGamp3(153))',(MGamp3(23)-MGamp3(153))*2d0
c$$$      print*,'(-MGamp3(25)+MGamp3(103))',(-MGamp3(25)+MGamp3(103))*2d0
c$$$      print*,'(-MGamp3(155)+MGamp3(257))',(-MGamp3(155)+MGamp3(257))*2d0
c$$$      print*,'color(2) pent sum',((MGamp3(23)-MGamp3(153))
c$$$     &     +(-MGamp3(25)+MGamp3(103))+(-MGamp3(101)+MGamp3(231))
c$$$     &     +(-MGamp3(155)+MGamp3(257))-(MGamp3(233)+MGamp3(259)))*2d0
c$$$c... color(3)
c$$$      print*,'(MGamp3(77)+MGamp3(127))',(MGamp3(77)+MGamp3(127))*2d0
c$$$      print*,'(MGamp3(23)-MGamp3(153))',(MGamp3(23)-MGamp3(153))*2d0
c$$$      print*,'-(MGamp3(25)+MGamp3(129))',-(MGamp3(25)+MGamp3(129))*2d0
c$$$      print*,'-(MGamp3(155)+MGamp3(207))',-(MGamp3(155)+MGamp3(207))*2d0
c$$$      print*,'(-MGamp3(75)+MGamp3(205))',(-MGamp3(75)+MGamp3(205))*2d0
c$$$      print*,'color(3) pent sum',((MGamp3(23)-MGamp3(153))
c$$$     &     -(MGamp3(25)+MGamp3(129))+(-MGamp3(75)+MGamp3(205))
c$$$     &     +(MGamp3(77)+MGamp3(127))-(MGamp3(155)+MGamp3(207)))*2d0
c$$$c... color(4)
c$$$      print*,'(MGamp3(51)-MGamp3(127))',(MGamp3(51)-MGamp3(127))*2d0
c$$$      print*,'-(MGamp3(23)+MGamp3(283))',-(MGamp3(23)+MGamp3(283))*2d0
c$$$      print*,'(MGamp3(25)+MGamp3(129))',(MGamp3(25)+MGamp3(129))*2d0
c$$$      print*,'(MGamp3(49)+MGamp3(309))',(MGamp3(49)+MGamp3(309))*2d0
c$$$      print*,'(MGamp3(285)-MGamp3(311))',(MGamp3(285)-MGamp3(311))*2d0
c$$$      print*,'color(4) pent sum',(-(MGamp3(23)+MGamp3(283))
c$$$     &     +(MGamp3(25)+MGamp3(129))+(MGamp3(49)+MGamp3(309))
c$$$     &     +(MGamp3(51)-MGamp3(127))+(MGamp3(285)-MGamp3(311)))*2d0
c$$$c... color(5)
c$$$      print*,'(-MGamp3(51)+MGamp3(127))',(-MGamp3(51)+MGamp3(127))*2d0
c$$$      print*,'(MGamp3(23)+MGamp3(283))',(MGamp3(23)+MGamp3(283))*2d0
c$$$      print*,'-(MGamp3(25)+MGamp3(129))',-(MGamp3(25)+MGamp3(129))*2d0
c$$$      print*,'-(MGamp3(49)+MGamp3(309))',-(MGamp3(49)+MGamp3(309))*2d0
c$$$      print*,'(-MGamp3(285)+MGamp3(311))',(-MGamp3(285)+MGamp3(311))*2d0
c$$$      print*,'color(5) pent sum',((MGamp3(23)+MGamp3(283))
c$$$     &     -(MGamp3(25)+MGamp3(129))-(MGamp3(49)+MGamp3(309))
c$$$     &     +(-MGamp3(51)+MGamp3(127))+(-MGamp3(285)+MGamp3(311)))*2d0
c$$$c... color(6)
c$$$      print*,'-(MGamp3(77)+MGamp3(127))',-(MGamp3(77)+MGamp3(127))*2d0
c$$$      print*,'(-MGamp3(23)+MGamp3(153))',(-MGamp3(23)+MGamp3(153))*2d0
c$$$      print*,'(MGamp3(25)+MGamp3(129))',(MGamp3(25)+MGamp3(129))*2d0
c$$$      print*,'(MGamp3(155)+MGamp3(207))',(MGamp3(155)+MGamp3(207))*2d0
c$$$      print*,'(MGamp3(75)-MGamp3(205))',(MGamp3(75)-MGamp3(205))*2d0
c$$$      print*,'color(6) pent sum',((-MGamp3(23)+MGamp3(153))
c$$$     &     +(MGamp3(25)+MGamp3(129))+(MGamp3(75)-MGamp3(205))
c$$$     &     -(MGamp3(77)+MGamp3(127))+(MGamp3(155)+MGamp3(207)))*2d0
c$$$c... color(7)
c$$$      print*,'(MGamp3(101)-MGamp3(231))',(MGamp3(101)-MGamp3(231))*2d0
c$$$      print*,'(MGamp3(233)+MGamp3(259))',(MGamp3(233)+MGamp3(259))*2d0
c$$$      print*,'(-MGamp3(23)+MGamp3(153))',(-MGamp3(23)+MGamp3(153))*2d0
c$$$      print*,'(MGamp3(25)-MGamp3(103))',(MGamp3(25)-MGamp3(103))*2d0
c$$$      print*,'(MGamp3(155)-MGamp3(257))',(MGamp3(155)-MGamp3(257))*2d0
c$$$      print*,'color(7) pent sum',((-MGamp3(23)+MGamp3(153))
c$$$     &     +(MGamp3(25)-MGamp3(103))+(MGamp3(101)-MGamp3(231))
c$$$     &     +(MGamp3(155)-MGamp3(257))+(MGamp3(233)+MGamp3(259)))*2d0
c$$$c... color(8)
c$$$      print*,'-(MGamp3(101)+MGamp3(361))',-(MGamp3(101)+MGamp3(361))*2d0
c$$$      print*,'(MGamp3(363)+MGamp3(389))',(MGamp3(363)+MGamp3(389))*2d0
c$$$      print*,'(MGamp3(23)+MGamp3(283))',(MGamp3(23)+MGamp3(283))*2d0
c$$$      print*,'(-MGamp3(25)+MGamp3(103))',(-MGamp3(25)+MGamp3(103))*2d0
c$$$      print*,'-(MGamp3(285)+MGamp3(387))',-(MGamp3(285)+MGamp3(387))*2d0
c$$$      print*,'color(8) pent sum',((MGamp3(23)+MGamp3(283))
c$$$     &     +(-MGamp3(25)+MGamp3(103))-(MGamp3(101)+MGamp3(361))
c$$$     &     -(MGamp3(285)+MGamp3(387))+(MGamp3(363)+MGamp3(389)))*2d0
c$$$c... color(9)
c$$$      print*,'(MGamp3(337)+MGamp3(389))',(MGamp3(337)+MGamp3(389))*2d0
c$$$      print*,'(MGamp3(51)+MGamp3(77))',(MGamp3(51)+MGamp3(77))*2d0
c$$$      print*,'(MGamp3(49)+MGamp3(309))',(MGamp3(49)+MGamp3(309))*2d0
c$$$      print*,'-(MGamp3(311)+MGamp3(387))',-(MGamp3(311)+MGamp3(387))*2d0
c$$$      print*,'-(MGamp3(75)+MGamp3(335))',-(MGamp3(75)+MGamp3(335))*2d0
c$$$      print*,'color(9) pent sum',((MGamp3(49)+MGamp3(309))
c$$$     &     +(MGamp3(51)+MGamp3(77))-(MGamp3(75)+MGamp3(335))
c$$$     &     -(MGamp3(311)+MGamp3(387))+(MGamp3(337)+MGamp3(389)))*2d0
c$$$c... color(10)
cFC      print*,'(MGamp3(101)+MGamp3(-231))',(MGamp3(101)-MGamp3(231))*2d0
cFC      print*,'(-MGamp3(181)+MGamp3(233))',(-MGamp3(181)+MGamp3(233))*2d0
cFC      print*,'(-MGamp3(51)+MGamp3(127))',(-MGamp3(51)+MGamp3(127))*2d0
cFC      print*,'(-MGamp3(103)-MGamp3(129))',(-MGamp3(103)-MGamp3(129))*2d0
cFC      print*,'(-MGamp3(49)+MGamp3(179))',(-MGamp3(49)+MGamp3(179))*2d0
cFC      print*,'color(10) pent sum',((-MGamp3(49)+MGamp3(179))
cFC     &     +(-MGamp3(51)+MGamp3(127))+(MGamp3(101)-MGamp3(231))
cFC     &     +(-MGamp3(103)-MGamp3(129))+(-MGamp3(181)+MGamp3(233)))*2d0
c$$$c... color(11)
c$$$      print*,'(MGamp3(181)+MGamp3(259))',(MGamp3(181)+MGamp3(259))*2d0
c$$$      print*,'(MGamp3(51)+MGamp3(77))',(MGamp3(51)+MGamp3(77))*2d0
c$$$      print*,'(-MGamp3(207)-MGamp3(257))',(-MGamp3(207)-MGamp3(257))*2d0
c$$$      print*,'(MGamp3(49)-MGamp3(179))',(MGamp3(49)-MGamp3(179))*2d0
c$$$      print*,'(-MGamp3(75)+MGamp3(205))',(-MGamp3(75)+MGamp3(205))*2d0
c$$$      print*,'color(11) pent sum',((MGamp3(49)-MGamp3(179))
c$$$     &     +(MGamp3(51)+MGamp3(77))+(-MGamp3(75)+MGamp3(205))
c$$$     &     +(MGamp3(181)+MGamp3(259))+(-MGamp3(207)-MGamp3(257)))*2d0
c$$$c... color(12)
c$$$      print*,'(-MGamp3(181)-MGamp3(259))',(-MGamp3(181)-MGamp3(259))*2d0
c$$$      print*,'(-MGamp3(51)-MGamp3(77))',(-MGamp3(51)-MGamp3(77))*2d0
c$$$      print*,'(MGamp3(207)+MGamp3(257))',(MGamp3(207)+MGamp3(257))*2d0
c$$$      print*,'(-MGamp3(49)+MGamp3(179))',(-MGamp3(49)+MGamp3(179))*2d0
c$$$      print*,'(MGamp3(75)-MGamp3(205))',(MGamp3(75)-MGamp3(205))*2d0
c$$$      print*,'color(12) pent sum',((-MGamp3(49)+MGamp3(179))
c$$$     &     +(-MGamp3(51)-MGamp3(77))+(MGamp3(75)-MGamp3(205))
c$$$     &     +(-MGamp3(181)-MGamp3(259))+(MGamp3(207)+MGamp3(257)))*2d0
c$$$c... color(13)
c$$$      print*,'(-MGamp3(101)+MGamp3(231))',(-MGamp3(101)+MGamp3(231))*2d0
c$$$      print*,'(MGamp3(181)-MGamp3(233))',(MGamp3(181)-MGamp3(233))*2d0
c$$$      print*,'(MGamp3(51)-MGamp3(127))',(MGamp3(51)-MGamp3(127))*2d0
c$$$      print*,'(MGamp3(103)+MGamp3(129))',(MGamp3(103)+MGamp3(129))*2d0
c$$$      print*,'(MGamp3(49)-MGamp3(179))',(MGamp3(49)-MGamp3(179))*2d0
c$$$      print*,'color(13) pent sum',((MGamp3(49)-MGamp3(179))
c$$$     &     +(MGamp3(51)-MGamp3(127))+(-MGamp3(101)+MGamp3(231))
c$$$     &     +(MGamp3(103)+MGamp3(129))+(MGamp3(181)-MGamp3(233)))*2d0
c$$$c... color(14)
c$$$      print*,'(-MGamp3(337)-MGamp3(389))',(-MGamp3(337)-MGamp3(389))*2d0
c$$$      print*,'(-MGamp3(51)-MGamp3(77))',(-MGamp3(51)-MGamp3(77))*2d0
c$$$      print*,'(-MGamp3(49)-MGamp3(309))',(-MGamp3(49)-MGamp3(309))*2d0
c$$$      print*,'(MGamp3(311)+MGamp3(387)',(MGamp3(311)+MGamp3(387))*2d0
c$$$      print*,'(MGamp3(75)+MGamp3(335))',(MGamp3(75)+MGamp3(335))*2d0
c$$$      print*,'color(14) pent sum',((-MGamp3(49)-MGamp3(309))
c$$$     &     +(-MGamp3(51)-MGamp3(77))+(MGamp3(75)+MGamp3(335))
c$$$     &     +(MGamp3(311)+MGamp3(387))+(-MGamp3(337)-MGamp3(389)))*2d0
c$$$c... color(15)
c$$$      print*,'(-MGamp3(101)-MGamp3(361))',(-MGamp3(101)-MGamp3(361))*2d0
c$$$      print*,'(-MGamp3(337)+MGamp3(363))',(-MGamp3(337)+MGamp3(363))*2d0
c$$$      print*,'(-MGamp3(77)-MGamp3(127))',(-MGamp3(77)-MGamp3(127))*2d0
c$$$      print*,'(MGamp3(103)+MGamp3(129))',(MGamp3(103)+MGamp3(129))*2d0
c$$$      print*,'(MGamp3(75)+MGamp3(335))',(MGamp3(75)+MGamp3(335))*2d0
c$$$      print*,'color(15) pent sum',((MGamp3(75)+MGamp3(335))
c$$$     &     +(-MGamp3(77)-MGamp3(127))+(-MGamp3(101)-MGamp3(361))
c$$$     &     +(MGamp3(103)+MGamp3(129))+(-MGamp3(337)+MGamp3(363)))*2d0
c$$$c... color(16)
c$$$      print*,'(MGamp3(101)+MGamp3(361))',(MGamp3(101)+MGamp3(361))*2d0
c$$$      print*,'(MGamp3(337)-MGamp3(363))',(MGamp3(337)-MGamp3(363))*2d0
c$$$      print*,'(MGamp3(77)+MGamp3(127))',(MGamp3(77)+MGamp3(127))*2d0
c$$$      print*,'(-MGamp3(103)-MGamp3(129))',(-MGamp3(103)-MGamp3(129))*2d0
c$$$      print*,'(-MGamp3(75)-MGamp3(335))',(-MGamp3(75)-MGamp3(335))*2d0
c$$$      print*,'color(16) pent sum',((MGamp3(101)+MGamp3(361))
c$$$     &     +(MGamp3(337)-MGamp3(363))+(MGamp3(77)+MGamp3(127))
c$$$     &     +(-MGamp3(103)-MGamp3(129))+(-MGamp3(75)-MGamp3(335)))*2d0
c$$$c... color(17)
c$$$      print*,'(MGamp3(337)+MGamp3(389))',(MGamp3(337)+MGamp3(389))*2d0
c$$$      print*,'(MGamp3(153)+MGamp3(283))',(MGamp3(153)+MGamp3(283))*2d0
c$$$      print*,'(MGamp3(155)+MGamp3(207))',(MGamp3(155)+MGamp3(207))*2d0
c$$$      print*,'(-MGamp3(285)-MGamp3(387))',(-MGamp3(285)-MGamp3(387))*2d0
c$$$      print*,'(-MGamp3(205)-MGamp3(335))',(-MGamp3(205)-MGamp3(335))*2d0
c$$$      print*,'color(17) pent sum',((MGamp3(153)+MGamp3(283))
c$$$     &     +(MGamp3(155)+MGamp3(207))+(-MGamp3(205)-MGamp3(335))
c$$$     &     +(-MGamp3(285)-MGamp3(387))+(MGamp3(337)+MGamp3(389)))*2d0
c$$$c... color(18)
c$$$      print*,'(MGamp3(181)+MGamp3(259))',(MGamp3(181)+MGamp3(259))*2d0
c$$$      print*,'(MGamp3(153)+MGamp3(283))',(MGamp3(153)+MGamp3(283))*2d0
c$$$      print*,'(MGamp3(155)-MGamp3(257))',(MGamp3(155)-MGamp3(257))*2d0
c$$$      print*,'(-MGamp3(179)-MGamp3(309))',(-MGamp3(179)-MGamp3(309))*2d0
c$$$      print*,'(-MGamp3(285)+MGamp3(311))',(-MGamp3(285)+MGamp3(311))*2d0
c$$$      print*,'color(18) pent sum',((MGamp3(153)+MGamp3(283))
c$$$     &     +(MGamp3(155)-MGamp3(257))+(-MGamp3(179)-MGamp3(309))
c$$$     &     +(MGamp3(181)+MGamp3(259))+(-MGamp3(285)+MGamp3(311)))*2d0
c$$$c... color(19)
c$$$      print*,'(-MGamp3(181)-MGamp3(259))',(-MGamp3(181)-MGamp3(259))*2d0
c$$$      print*,'(-MGamp3(153)-MGamp3(283))',(-MGamp3(153)-MGamp3(283))*2d0
c$$$      print*,'(-MGamp3(155)+MGamp3(257))',(-MGamp3(155)+MGamp3(257))*2d0
c$$$      print*,'(MGamp3(179)+MGamp3(309)',(MGamp3(179)+MGamp3(309))*2d0
c$$$      print*,'(MGamp3(285)-MGamp3(311))',(MGamp3(285)-MGamp3(311))*2d0
c$$$      print*,'color(19) pent sum',((-MGamp3(181)-MGamp3(259))
c$$$     &     +(-MGamp3(153)-MGamp3(283))+(-MGamp3(155)+MGamp3(257))
c$$$     &     +(MGamp3(179)+MGamp3(309))+(MGamp3(285)-MGamp3(311)))*2d0
c$$$c... color(20)
c$$$      print*,'(-MGamp3(337)-MGamp3(389))',(-MGamp3(337)-MGamp3(389))*2d0
c$$$      print*,'(-MGamp3(153)-MGamp3(283))',(-MGamp3(153)-MGamp3(283))*2d0
c$$$      print*,'(-MGamp3(155)-MGamp3(207))',(-MGamp3(155)-MGamp3(207))*2d0
c$$$      print*,'(MGamp3(285)+MGamp3(387))',(MGamp3(285)+MGamp3(387))*2d0
c$$$      print*,'(MGamp3(205)+MGamp3(335))',(MGamp3(205)+MGamp3(335))*2d0
c$$$      print*,'color(20) pent sum',((-MGamp3(153)-MGamp3(283))
c$$$     &     +(-MGamp3(155)-MGamp3(207))+(MGamp3(205)+MGamp3(335))
c$$$     &     +(MGamp3(285)+MGamp3(387))+(-MGamp3(337)-MGamp3(389)))*2d0
c$$$c... color(21)
c$$$      print*,'(MGamp3(231)+MGamp3(361))',(MGamp3(231)+MGamp3(361))*2d0
c$$$      print*,'(-MGamp3(363)-MGamp3(389))',(-MGamp3(363)-MGamp3(389))*2d0
c$$$      print*,'(MGamp3(181)-MGamp3(233))',(MGamp3(181)-MGamp3(233))*2d0
c$$$      print*,'(-MGamp3(179)-MGamp3(309))',(-MGamp3(179)-MGamp3(309))*2d0
c$$$      print*,'(MGamp3(311)+MGamp3(387))',(MGamp3(311)+MGamp3(387))*2d0
c$$$      print*,'color(21) pent sum',((-MGamp3(179)-MGamp3(309))
c$$$     &     +(MGamp3(181)-MGamp3(233))+(MGamp3(231)+MGamp3(361))
c$$$     &     +(MGamp3(311)+MGamp3(387))+(-MGamp3(363)-MGamp3(389)))*2d0
c$$$c... color(22)
c$$$      print*,'(-MGamp3(231)-MGamp3(361))',(-MGamp3(231)-MGamp3(361))*2d0
c$$$      print*,'(MGamp3(363)+MGamp3(389))',(MGamp3(363)+MGamp3(389))*2d0
c$$$      print*,'(-MGamp3(181)+MGamp3(233))',(-MGamp3(181)+MGamp3(233))*2d0
c$$$      print*,'(MGamp3(179)+MGamp3(309))',(MGamp3(179)+MGamp3(309))*2d0
c$$$      print*,'(-MGamp3(311)-MGamp3(387))',(-MGamp3(311)-MGamp3(387))*2d0
c$$$      print*,'color(22) pent sum',((MGamp3(179)+MGamp3(309))
c$$$     &     +(-MGamp3(181)+MGamp3(233))+(-MGamp3(231)-MGamp3(361))
c$$$     &     +(-MGamp3(311)-MGamp3(387))+(MGamp3(363)+MGamp3(389)))*2d0
c$$$c... color(23)
c$$$      print*,'(-MGamp3(231)-MGamp3(361))',(-MGamp3(231)-MGamp3(361))*2d0
c$$$      print*,'(-MGamp3(337)+MGamp3(363))',(-MGamp3(337)+MGamp3(363))*2d0
c$$$      print*,'(MGamp3(233)+MGamp3(259))',(MGamp3(233)+MGamp3(259))*2d0
c$$$      print*,'(-MGamp3(207)-MGamp3(257))',(-MGamp3(207)-MGamp3(257))*2d0
c$$$      print*,'(MGamp3(205)-MGamp3(335))',(MGamp3(205)+MGamp3(335))*2d0
c$$$      print*,'color(23) pent sum',((MGamp3(205)+MGamp3(335))
c$$$     &     +(-MGamp3(207)-MGamp3(257))+(-MGamp3(231)-MGamp3(361))
c$$$     &     +(MGamp3(233)+MGamp3(259))+(-MGamp3(337)+MGamp3(363)))*2d0
c$$$c... color(24)
cFC      print*,'(MGamp3(231)+MGamp3(361))',(MGamp3(231)+MGamp3(361))*2d0
cFC      print*,'(MGamp3(337)-MGamp3(363))',(MGamp3(337)-MGamp3(363))*2d0
cFC      print*,'(-MGamp3(233)-MGamp3(259))',(-MGamp3(233)-MGamp3(259))*2d0
cFC      print*,'(MGamp3(207)+MGamp3(257))',(MGamp3(207)+MGamp3(257))*2d0
cFC      print*,'(-MGamp3(205)-MGamp3(335))',(-MGamp3(205)-MGamp3(335))*2d0
cFC      print*,'color(24) pent sum',((-MGamp3(205)-MGamp3(335))
cFC     &     +(MGamp3(207)+MGamp3(257))+(MGamp3(231)+MGamp3(361))
cFC     &     +(-MGamp3(233)-MGamp3(259))+(MGamp3(337)-MGamp3(363)))*2d0
cFC
         iF(HloopTop) THEN
      print*
      do i=1,24
         if(abs(colorN(i)/colorM(i) -1d0).gt.2d-6) then
            print*,"colorM(",i,")=", colorM(i)
            print*,"colorN(",i,")=", colorN(i)
            print*,"ratio colorM/colorN -1d0", abs(colorN(i)/colorM(i)
     &           -1d0), colorN(i)/colorM(i) -1d0
            print*
         endif
      enddo
      ENDIF
      PAUSE
      
      endif
*---------------------------------------------------------------

      END
