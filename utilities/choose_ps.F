
****************************************************************************
*   This file contains a function to decide which phase space should be
*   taken in case of multi-channels.
c***************************************************************************
c   LIST OF ALL FUNCTIONS AND SUBROUTINES IN THIS FILE:
c
c     LOGICAL FUNCTION Choose_PS(PS_number, v, L)
c     LOGICAL FUNCTION VA_Choose(PS_number, mom, L, bos1)
c     LOGICAL FUNCTION VVA_Choose(PS_number, mom, L, bos1, bos2)
c     LOGICAL FUNCTION WWA_Choose_anom(PS_number, mom, L, bos1, bos2)
c     LOGICAL FUNCTION WZA_Choose_anom(PS_number, mom, L, bos1, bos2)
c     LOGICAL FUNCTION ZZA_Choose_anom(PS_number, mom, L, bos1, bos2)
c     LOGICAL FUNCTION VAA_Choose(PS_number, mom, L, bos1)
c     LOGICAL FUNCTION WAA_Choose_anom(PS_number, mom, L, bos1)
c     LOGICAL FUNCTION ZAA_Choose_anom(PS_number, mom, L, bos1)
c     LOGICAL FUNCTION ZnAA_Choose_anom(PS_number, mom, L, bos1)
c     LOGICAL FUNCTION AAA_Choose_anom(PS_number, mom, L, bos1)
c
c***************************************************************************
*   Vera Hankele, <vera@particle.uni-karlsruhe.de>
*   Initial version:  09 Februar 2009
*   Last modified: 18.07.2013 by Nicolas Kaiser
****************************************************************************

c*************************************************************************
      LOGICAL FUNCTION Choose_PS(PS_number, v, L)
c*************************************************************************
c     This is the general routine to choose a certain phase space generator.
c*************************************************************************
c     INPUT 
c     PS_number:  Number of phase space generator used to generate this point
c     v    :      4-momenta of the particles, which are no jets. 
c     L    :      id for momentum set
c
c     OUTPUT
c     Cuts   : .true. or .false. depending on whether this phase space
c              generator is selected or not.
c*************************************************************************

      IMPLICIT NONE

#include "global.inc"
#include "process.inc"

      integer ps_number

c     input variables
      double precision v(0:3,max_v,max_kin) 
      integer L

c     external functions
      logical VA_Choose, VVA_Choose, VAA_Choose, WWA_Choose_anom
      logical WZA_Choose_anom, ZZA_Choose_anom, WAA_Choose_anom
      logical ZAA_Choose_anom, ZnAA_Choose_anom, AAA_Choose_anom
      external VA_Choose, VVA_Choose, VAA_Choose, WWA_Choose_anom
      external WZA_Choose_anom, ZZA_Choose_anom, WAA_Choose_anom 
      external ZAA_Choose_anom, ZnAA_Choose_anom, AAA_Choose_anom

c For most processes no multi-channel approach is used -> .TRUE.
      Choose_PS = .true.

      SELECT CASE(procID)
      CASE(WPA,WMA,WPAJ,WMAJ,WPAJJ,WMAJJ,QCDWPAJJ,QCDWMAJJ,QCDZAJJ_l)! charge of W does not matter
        if(PS_loops.eq.1) then
          Choose_PS = .true.
        else
          Choose_PS = VA_Choose(PS_number, v(0,1,L), L , 3)
        endif
      CASE(QCDZAJJ_n)! charge of W does not matter
          Choose_PS = .true.
      CASE(ZA,ZAjj)
         Choose_PS = VA_Choose(PS_number, v(0,1,L), L , 2)
      CASE(WWA, WPhadWMA, WPWMhadA)
         if(with_anom) then
            Choose_PS = WWA_Choose_anom(PS_number, v(0,1,L), L , 3, 3)
         else
            Choose_PS = VVA_Choose(PS_number, v(0,1,L), L , 3, 3)
         endif
      CASE(ZZA, ZZhadA, ZZnA)
         if(with_anom) then
            Choose_PS = ZZA_Choose_anom(PS_number, v(0,1,L), L , 2, 2)
         else
            Choose_PS = VVA_Choose(PS_number, v(0,1,L), L , 2, 2)
         endif
      CASE(WPZA,WMZA, WPhadZA, WPZhadA, WMhadZA, WMZhadA)  ! charge of W does not matter
         if(with_anom) then
            Choose_PS = WZA_Choose_anom(PS_number, v(0,1,L), L , 3, 2)
         else
            Choose_PS = VVA_Choose(PS_number, v(0,1,L), L , 3, 2)
         endif
      CASE(WPAA,WMAA,WPAAj,WMAAj)  ! charge of W does not matter
         if(with_anom) then
            Choose_PS = WAA_Choose_anom(PS_number, v(0,1,L), L , 3)
         else
            Choose_PS = VAA_Choose(PS_number, v(0,1,L), L , 3)
         endif
      CASE(ZAA)  
         if(with_anom) then
            Choose_PS = ZAA_Choose_anom(PS_number, v(0,1,L), L , 2)
         else
            Choose_PS = VAA_Choose(PS_number, v(0,1,L), L , 2)
         endif
      CASE(ZnAA)  
         if(with_anom) then
            Choose_PS = ZnAA_Choose_anom(PS_number, v(0,1,L), L , 2)
         else
            Choose_PS = .true.
         endif
      CASE(AAA)  
         if(with_anom) then
            Choose_PS = AAA_Choose_anom(PS_number, v(0,1,L), L , 1)
         else
            Choose_PS = .true.
         endif
      END SELECT

      RETURN
      END


c*************************************************************************
      LOGICAL FUNCTION VA_Choose(PS_number, mom, L, bos1)
c*************************************************************************
c     This routine selects the right phase space generator for VA
c     processes
c     Note: PS1 is (must be!) the three-body decay here!
c           GF-ZA only has one phase-space, which is that one.
c*************************************************************************
      IMPLICIT NONE

c     input variables
      integer PS_number, L
      double precision mom(0:3,3)
      integer bos1

c local variables
      double precision tempvector(0:4)
      double precision tempva(0:4)
      double precision choosev1,chooseva
      integer mu, myps

c     external functions
      double precision dotrr
      external dotrr


      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)


      do mu=0,3
         tempvector(mu)=mom(mu,1)+mom(mu,2)
         tempva(mu)=tempvector(mu)+mom(mu,3)
      enddo
      tempvector(4) = dotrr(tempvector,tempvector)
      tempva(4) = dotrr(tempva,tempva)

      choosev1=abs(sqrt(tempvector(4))-sqrt(xm2(bos1)))
      chooseva=abs(sqrt(tempva(4))-sqrt(xm2(bos1)))

      if(chooseva.le.choosev1) then
        myps = 1
      else 
        myps = 2
      endif

      VA_Choose = myps.eq.ps_number
      
      RETURN
      END

c*************************************************************************
      LOGICAL FUNCTION VVA_Choose(PS_number, mom, L, bos1, bos2)
c*************************************************************************
c     This routine selects the right phase space generator for VVA
c     processes
c     phase spaces: 1: photon from continuum
c                   2: photon from bos1
c                   3: photon from bos2
c*************************************************************************
      IMPLICIT NONE

#include "global.inc"

      integer ps_number

c     input variables
      integer L
      double precision mom(0:3,5)
      integer bos1, bos2

c local variables
      double precision tempp(0:3),tempm(0:3)
      double precision temppa(0:3),tempma(0:3)
      double precision choosep,choosem,choosepa,choosema
      integer mu, myps
      logical VVA_ps

c     external functions
      double precision dotrr
      external dotrr

      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

      VVA_ps = .true.

      do mu=0,3
         tempp(mu)=mom(mu,1)+mom(mu,2)
         tempm(mu)=mom(mu,3)+mom(mu,4)
         temppa(mu)=tempp(mu)+mom(mu,5)
         tempma(mu)=tempm(mu)+mom(mu,5)
      enddo

      choosep=abs(sqrt(dotrr(tempp,tempp))-sqrt(xm2(bos1)))
      choosem=abs(sqrt(dotrr(tempm,tempm))-sqrt(xm2(bos2)))
      choosepa=abs(sqrt(dotrr(temppa,temppa))-sqrt(xm2(bos1)))
      choosema=abs(sqrt(dotrr(tempma,tempma))-sqrt(xm2(bos2)))

      
      if ( (choosepa.lt.choosep).and.(choosepa.lt.choosema) ) then
        myps = 2
      elseif ( (choosema.lt.choosem).and.(choosema.lt.choosepa) ) then
        myps = 3
      else
        myps = 1
      endif
      if (ps_number.ne.myps) VVA_ps=.false.

      VVA_Choose = VVA_ps
      
      RETURN
      END

c*************************************************************************
      LOGICAL FUNCTION WWA_Choose_anom(PS_number, mom, L, bos1, bos2)
c*************************************************************************
c     This routine selects the right phase space generator for WWA
c     production with anomalous couplings
c     phase spaces: 1: photon from continuum, Higgs in WW
c                   2: photon from bos1
c                   3: photon from bos2
c*************************************************************************
      IMPLICIT NONE

#include "global.inc"

      integer ps_number

c     input variables
      integer L
      double precision mom(0:3,5)
      integer bos1, bos2

c local variables
      double precision tempp(0:3),tempm(0:3)
      double precision temppa(0:3),tempma(0:3),temppm(0:3)
      double precision choosep,choosem,choosepa,choosema,chooseh
      double precision gammap,gammam,gammah
      integer mu, myps
      logical VVA_ps

c     external functions
      double precision dotrr
      external dotrr

      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

      VVA_ps = .true.

      do mu=0,3
         tempp(mu)=mom(mu,1)+mom(mu,2)
         tempm(mu)=mom(mu,3)+mom(mu,4)
         temppm(mu)=mom(mu,1)+mom(mu,2)+mom(mu,3)+mom(mu,4)
         temppa(mu)=tempp(mu)+mom(mu,5)
         tempma(mu)=tempm(mu)+mom(mu,5)
      enddo

      gammap = xmg(bos1)/sqrt(xm2(bos1))
      gammam = xmg(bos2)/sqrt(xm2(bos2))
      gammah = min(xmg(6)/sqrt(xm2(6)),10.) ! avoid that this matches too often due to large Higgs width

      choosep=abs(sqrt(dotrr(tempp,tempp))-sqrt(xm2(bos1)))
      choosem=abs(sqrt(dotrr(tempm,tempm))-sqrt(xm2(bos2)))
      choosepa=abs(sqrt(dotrr(temppa,temppa))-sqrt(xm2(bos1)))
      choosema=abs(sqrt(dotrr(tempma,tempma))-sqrt(xm2(bos2)))

      chooseh=abs(sqrt(dotrr(temppm,temppm))-sqrt(xm2(6)))/gammah

      
      if ( (choosepa.lt.choosep).and.(choosepa.lt.choosema) .and. chooseh.gt.min(choosepa/gammap,choosema/gammam) ) then
        myps = 2
      elseif ( (choosema.lt.choosem).and.(choosema.lt.choosepa)  .and. chooseh.gt.min(choosepa/gammap,choosema/gammam) ) then
        myps = 3
      else
        myps = 1
      endif
      if (ps_number.ne.myps) VVA_ps=.false.

      WWA_Choose_anom = VVA_ps
      
      RETURN
      END

c*************************************************************************
      LOGICAL FUNCTION WZA_Choose_anom(PS_number, mom, L, bos1, bos2)
c*************************************************************************
c     This routine selects the right phase space generator for the WZA
c     process with anomalous couplings 
c*************************************************************************
      IMPLICIT NONE

#include "global.inc"
#include "coupl.inc"

      integer ps_number

c     input variables
      integer L
      double precision mom(0:3,5)
      integer bos1, bos2

c local variables
      double precision tempp(0:3),tempm(0:3)
      double precision temppa(0:3),tempma(0:3)
      double precision choosep,choosem,choosepa,choosema,chooseh
      double precision gammap,gammam,gammah
      integer mu, myps
      logical WZA_ps

c     external functions
      double precision dotrr
      external dotrr

      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

      WZA_ps = .true.

c widths of particles
      gammap = xmg(bos1)/sqrt(xm2(bos1))
      gammam = xmg(bos2)/sqrt(xm2(bos2))
      gammah = min(xmg(6)/sqrt(xm2(6)),10.) ! avoid that this matches too often due to large Higgs width

      do mu=0,3
         tempp(mu)=mom(mu,1)+mom(mu,2)
         tempm(mu)=mom(mu,3)+mom(mu,4)
         temppa(mu)=tempp(mu)+mom(mu,5)
         tempma(mu)=tempm(mu)+mom(mu,5)
      enddo

      choosep=abs(sqrt(dotrr(tempp,tempp))-sqrt(xm2(bos1)))/gammap
      choosem=abs(sqrt(dotrr(tempm,tempm))-sqrt(xm2(bos2)))/gammam
      choosepa=abs(sqrt(dotrr(temppa,temppa))-sqrt(xm2(bos1)))/gammap
      choosema=abs(sqrt(dotrr(tempma,tempma))-sqrt(xm2(bos2)))/gammam
      chooseh=abs(sqrt(dotrr(tempma,tempma))-sqrt(xm2(6)))/gammah

      if ( (choosepa.lt.choosep).and.(choosepa.lt.choosema).and.(choosepa.lt.chooseh) ) then
        myps = 2
      elseif ( ((choosema.lt.choosem).and.(choosema.lt.choosepa))
     &    .or. ((chooseh.lt.choosema).and.(chooseh.lt.choosepa)) ) then
        myps = 3
      else
        myps = 1
      endif
      if (ps_number.ne.myps) WZA_ps=.false.

      WZA_Choose_anom = WZA_ps
      
      RETURN
      END


c*************************************************************************
      LOGICAL FUNCTION ZZA_Choose_anom(PS_number, mom, L, bos1, bos2)
c*************************************************************************
c     This routine selects the right phase space generator for ZZA
c     with anomalous couplings
c*************************************************************************
      IMPLICIT NONE

#include "global.inc"
#include "coupl.inc"

      integer ps_number

c     input variables
      integer L
      double precision mom(0:3,5)
      integer bos1, bos2

c local variables
      double precision tempp(0:3),tempm(0:3),temppm(0:3)
      double precision temppa(0:3),tempma(0:3)
      double precision choosep,choosem,choosepa,choosema,chooseh1,chooseh2,chooseh3
      double precision gammap,gammam,gammah
      integer mu, myps
      logical ZZA_ps

c     external functions
      double precision dotrr
      external dotrr

      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

      ZZA_ps = .true.

c widths of particles
      gammap = xmg(bos1)/sqrt(xm2(bos1))
      gammam = xmg(bos2)/sqrt(xm2(bos2))
      gammah = min(xmg(6)/sqrt(xm2(6)),10.) ! avoid that this matches too often due to large Higgs width

      do mu=0,3
         tempp(mu)=mom(mu,1)+mom(mu,2)
         tempm(mu)=mom(mu,3)+mom(mu,4)
         temppm(mu)=tempp(mu)+tempm(mu)
         temppa(mu)=tempp(mu)+mom(mu,5)
         tempma(mu)=tempm(mu)+mom(mu,5)       
      enddo

      choosep=abs(sqrt(dotrr(tempp,tempp))-sqrt(xm2(bos1)))/gammap
      choosem=abs(sqrt(dotrr(tempm,tempm))-sqrt(xm2(bos2)))/gammam
      choosepa=abs(sqrt(dotrr(temppa,temppa))-sqrt(xm2(bos1)))/gammap
      choosema=abs(sqrt(dotrr(tempma,tempma))-sqrt(xm2(bos2)))/gammam
      chooseh1=abs(sqrt(dotrr(temppm,temppm))-sqrt(xm2(6)))/gammah
      chooseh2=abs(sqrt(dotrr(tempp,tempp))-sqrt(xm2(6)))/gammah
      chooseh3=abs(sqrt(dotrr(tempm,tempm))-sqrt(xm2(6)))/gammah

      if ( ((choosepa.lt.choosep).and.(choosepa.lt.choosema)) 
     &  .or. (chooseh2.lt.min(choosep,chooseh1,chooseh3)) ) then
         myps=2
      elseif ( ((choosema.lt.choosem).and.(choosema.lt.choosepa)) 
     &  .or. (chooseh3.lt.min(choosem,chooseh1,chooseh2)) ) then
         myps=3
      else
         myps=1 
      endif 

      if (ps_number.ne.myps) ZZA_ps=.false.

      ZZA_Choose_anom = ZZA_ps
      
      RETURN
      END

c*************************************************************************
      LOGICAL FUNCTION VAA_Choose(PS_number, mom, L, bos1)
c*************************************************************************
c     This routine selects the right phase space generator for VAA
c     processes
c*************************************************************************
      IMPLICIT NONE

c     input variables
      integer PS_number, L
      double precision mom(0:3,4)
      integer bos1

c local variables
      double precision tempvector(0:4)
      double precision tempv2(0:4), tempv3(0:4), tempva(0:4)
      double precision choosev1,choosev2,choosev3,chooseva
      integer mu, myps
      logical VAA_ps

c     external functions
      double precision dotrr
      external dotrr


      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)


      VAA_ps = .true.

      do mu=0,3
         tempvector(mu)=mom(mu,1)+mom(mu,2)
         tempv2(mu)=tempvector(mu)+mom(mu,3)
         tempv3(mu)=tempvector(mu)+mom(mu,4)
         tempva(mu)=tempv2(mu)+mom(mu,4)
      enddo
      tempvector(4) = dotrr(tempvector,tempvector)
      tempv2(4) = dotrr(tempv2,tempv2)
      tempv3(4) = dotrr(tempv3,tempv3)
      tempva(4) = dotrr(tempva,tempva)

      choosev1=abs(sqrt(tempvector(4))-sqrt(xm2(bos1)))
      choosev2=abs(sqrt(tempv2(4))-sqrt(xm2(bos1)))
      choosev3=abs(sqrt(tempv3(4))-sqrt(xm2(bos1)))
      chooseva=abs(sqrt(tempva(4))-sqrt(xm2(bos1)))

      if(choosev2.le.min(choosev1,choosev3,chooseva)) then
        myps = 2
      else if(choosev3.le.min(choosev1,choosev2,chooseva)) then
        myps = 3
      else if(chooseva.le.min(choosev1,choosev2,choosev3)) then
        myps = 4
      else 
        myps = 1
      endif
      if (ps_number.ne.myps) VAA_ps=.false.

      VAA_Choose = VAA_ps
      
      RETURN
      END

c*************************************************************************
      LOGICAL FUNCTION WAA_Choose_anom(PS_number, mom, L, bos1)
c*************************************************************************
c     This routine selects the right phase space generator for WAA
c     with anomalous couplings
c*************************************************************************
      IMPLICIT NONE

c     input variables
      integer PS_number, L
      double precision mom(0:3,4)
      integer bos1

c local variables
      double precision tempvector(0:4)
      double precision tempv2(0:4), tempv3(0:4), tempva(0:4), tempaa(0:4)
      double precision choosev1,choosev2,choosev3,chooseva,chooseaa
      integer mu, myps
      logical WAA_ps

c     external functions
      double precision dotrr
      external dotrr


      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)


      WAA_ps = .true.

      do mu=0,3
         tempvector(mu)=mom(mu,1)+mom(mu,2)
         tempv2(mu)=tempvector(mu)+mom(mu,3)
         tempv3(mu)=tempvector(mu)+mom(mu,4)
         tempva(mu)=tempv2(mu)+mom(mu,4)
         tempaa(mu)=mom(mu,3)+mom(mu,4)
      enddo
      tempvector(4) = dotrr(tempvector,tempvector)
      tempv2(4) = dotrr(tempv2,tempv2)
      tempv3(4) = dotrr(tempv3,tempv3)
      tempva(4) = dotrr(tempva,tempva)
      tempaa(4) = dotrr(tempaa,tempaa)

      choosev1=abs(sqrt(tempvector(4))-sqrt(xm2(bos1)))
      choosev2=abs(sqrt(tempv2(4))-sqrt(xm2(bos1)))
      choosev3=abs(sqrt(tempv3(4))-sqrt(xm2(bos1)))
      chooseva=abs(sqrt(tempva(4))-sqrt(xm2(bos1)))
      chooseaa=abs(sqrt(tempaa(4))-sqrt(xm2(6))) !for H -> A A

      if((choosev2.le.min(choosev1,choosev3,chooseva)).and.(chooseaa.gt.7.0D-2)) then
        myps = 2
      else if((choosev3.le.min(choosev1,choosev2,chooseva)).and.(chooseaa.gt.7.0D-2)) then
        myps = 3
      else if((chooseva.le.min(choosev1,choosev2,choosev3)).and.(chooseaa.gt.7.0D-2)) then
        myps = 4
      else if(chooseaa.le.7.0D-2) then
        myps = 5
      else 
        myps = 1
      endif
      if (ps_number.ne.myps) WAA_ps=.false.

      WAA_Choose_anom = WAA_ps
      
      RETURN
      END

c*************************************************************************
      LOGICAL FUNCTION ZAA_Choose_anom(PS_number, mom, L, bos1)
c*************************************************************************
c     This routine selects the right phase space generator for ZAA
c     with anomalous couplings
c*************************************************************************
      IMPLICIT NONE

c     input variables
      integer PS_number, L
      double precision mom(0:3,4)
      integer bos1

c local variables
      double precision tempvector(0:4)
      double precision tempv2(0:4), tempv3(0:4), tempva(0:4), tempaa(0:4)
      double precision choosev1,choosev2,choosev3,chooseva,choosehaa,choosehv2,choosehv3
      integer mu, myps
      logical ZAA_ps

c     external functions
      double precision dotrr
      external dotrr


      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)


      ZAA_ps = .true.

      do mu=0,3
         tempvector(mu)=mom(mu,1)+mom(mu,2)
         tempv2(mu)=tempvector(mu)+mom(mu,3)
         tempv3(mu)=tempvector(mu)+mom(mu,4)
         tempva(mu)=tempv2(mu)+mom(mu,4)
         tempaa(mu)=mom(mu,3)+mom(mu,4)
      enddo
      tempvector(4) = dotrr(tempvector,tempvector)
      tempv2(4) = dotrr(tempv2,tempv2)
      tempv3(4) = dotrr(tempv3,tempv3)
      tempva(4) = dotrr(tempva,tempva)
      tempaa(4) = dotrr(tempaa,tempaa)

      choosev1=abs(sqrt(tempvector(4))-sqrt(xm2(bos1)))
      choosev2=abs(sqrt(tempv2(4))-sqrt(xm2(bos1)))
      choosev3=abs(sqrt(tempv3(4))-sqrt(xm2(bos1)))
      chooseva=abs(sqrt(tempva(4))-sqrt(xm2(bos1)))
      choosehaa=abs(sqrt(tempaa(4))-sqrt(xm2(6))) ! for H -> AA
      choosehv2=abs(sqrt(tempv2(4))-sqrt(xm2(6))) ! for H -> ZA1
      choosehv3=abs(sqrt(tempv3(4))-sqrt(xm2(6))) ! for H -> ZA2

      if(((choosev2.le.min(choosev1,choosev3,chooseva)).and.(min(choosehaa,choosehv2,choosehv3).gt.7.0D-2))
     1   .or.((choosehv2.le.7.0D-2).and.(choosehv2.lt.min(choosehaa,choosehv3)))) then
        myps = 2
      else if(((choosev3.le.min(choosev1,choosev2,chooseva)).and.(min(choosehaa,choosehv2,choosehv3).gt.7.0D-2))
     1       .or.((choosehv3.le.7.0D-2).and.(choosehv3.lt.min(choosehaa,choosehv2)))) then
        myps = 3
      else if((chooseva.le.min(choosev1,choosev2,choosev3)).and.(min(choosehaa,choosehv2,choosehv3).gt.7.0D-2)) then
        myps = 4
      else if((choosehaa.le.7.0D-2).and.(choosehaa.le.min(choosehv2,choosehv3))) then
        myps = 5
      else 
        myps = 1
      endif
      if (ps_number.ne.myps) ZAA_ps=.false.

      ZAA_Choose_anom = ZAA_ps
      
      RETURN
      END


c*************************************************************************
      LOGICAL FUNCTION ZnAA_Choose_anom(PS_number, mom, L, bos1)
c*************************************************************************
c     This routine selects the right phase space generator for ZAA (Z->nu nu~)
c     with anomalous couplings
c*************************************************************************
      IMPLICIT NONE

c     input variables
      integer PS_number, L
      double precision mom(0:3,4)
      integer bos1

c local variables
      double precision tempvector(0:4)
      double precision tempv2(0:4), tempv3(0:4), tempva(0:4), tempaa(0:4)
      double precision choosev1,choosehaa,choosehv2,choosehv3
      integer mu, myps
      logical ZnAA_ps

c     external functions
      double precision dotrr
      external dotrr


      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)


      ZnAA_ps = .true.

      do mu=0,3
         tempvector(mu)=mom(mu,1)+mom(mu,2)
         tempv2(mu)=tempvector(mu)+mom(mu,3)
         tempv3(mu)=tempvector(mu)+mom(mu,4)
         tempva(mu)=tempv2(mu)+mom(mu,4)
         tempaa(mu)=mom(mu,3)+mom(mu,4)
      enddo
      tempvector(4) = dotrr(tempvector,tempvector)
      tempv2(4) = dotrr(tempv2,tempv2)
      tempv3(4) = dotrr(tempv3,tempv3)
      tempva(4) = dotrr(tempva,tempva)
      tempaa(4) = dotrr(tempaa,tempaa)

      choosev1=abs(sqrt(tempvector(4))-sqrt(xm2(bos1)))
      choosehaa=abs(sqrt(tempaa(4))-sqrt(xm2(6))) ! for H -> AA
      choosehv2=abs(sqrt(tempv2(4))-sqrt(xm2(6))) ! for H -> ZA1
      choosehv3=abs(sqrt(tempv3(4))-sqrt(xm2(6))) ! for H -> ZA2

      if(choosehv2.le.min(choosehaa,choosehv3).and.choosehv2.le.7.0D-2) then
        myps = 2
      else if(choosehv3.le.min(choosehaa,choosehv2).and.choosehv3.le.7.0D-2) then
        myps = 3
      else if(choosehaa.le.min(choosehv2,choosehv3).and.choosehaa.le.7.0D-2) then
        myps = 4
      else 
        myps = 1
      endif
      if (ps_number.ne.myps) ZnAA_ps=.false.

      ZnAA_Choose_anom = ZnAA_ps
      
      RETURN
      END

c*************************************************************************
      LOGICAL FUNCTION AAA_Choose_anom(PS_number, mom, L, bos1)
c*************************************************************************
c     This routine selects the right phase space generator for AAA 
c     with anomalous couplings
c*************************************************************************
      IMPLICIT NONE

c     input variables
      integer PS_number, L
      double precision mom(0:3,3)
      integer bos1

c local variables
      double precision tempv12(0:4), tempv13(0:4), tempv23(0:4)
      double precision chooseh12,chooseh13,chooseh23
      integer mu, myps
      logical AAA_ps

c     external functions
      double precision dotrr
      external dotrr


      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)


      AAA_ps = .true.

      do mu=0,3
         tempv12(mu)=mom(mu,1)+mom(mu,2)
         tempv13(mu)=mom(mu,1)+mom(mu,3)
         tempv23(mu)=mom(mu,2)+mom(mu,3)
      enddo
      tempv12(4) = dotrr(tempv12,tempv12)
      tempv13(4) = dotrr(tempv13,tempv13)
      tempv23(4) = dotrr(tempv23,tempv23)


      chooseh12=abs(sqrt(tempv12(4))-sqrt(xm2(6)))! for H -> A1 A2
      chooseh13=abs(sqrt(tempv13(4))-sqrt(xm2(6))) ! for H -> A1 A3
      chooseh23=abs(sqrt(tempv23(4))-sqrt(xm2(6))) ! for H -> A2 A3


      if(chooseh12.le.min(chooseh13,chooseh23).and.chooseh12.le.7.0D-2) then
        myps = 2
      else if(chooseh13.le.min(chooseh12,chooseh23).and.chooseh13.le.7.0D-2) then
        myps = 3
      else if(chooseh23.le.min(chooseh12,chooseh13).and.chooseh23.le.7.0D-2) then
        myps = 4
      else 
        myps = 1
      endif
      if (ps_number.ne.myps) AAA_ps=.false.

      AAA_Choose_anom = AAA_ps
      
      RETURN
      END
