c
c phase space generator for pp --> V1 + V2 + (np-nx) jets, 
c                                  V1 --> l1 l2, V2 --> l3 l4
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  1995 March 25
c	Last modified:  2005 Dec.        
c
      subroutine def_ps_vvnj(bos0,bos1,bos2, iwidth0, iwidth1, iwidth2)
          use globalvars, only: ldoblha
      implicit none 
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/cuts.inc"
#include "VBFNLO/utilities/spin2coupl.inc"
#include "VBFNLO/utilities/process.inc"
      integer N1, bos0, bos1, bos2, bos(0:2), i
      integer iwidth0, iwidth1, iwidth2 
      parameter ( N1=max_p-nx )
c
c set up variables in common block BLIPSIvvnj for call of LIPSNVV
c 
      DOUBLE PRECISION RM2,RMG,RM2MIN,RM2MAX, S,M2MIN, YJMIN,YJMAX,PTJMIN,EJMIN
      INTEGER INFOJ, iwidth
      COMMON /BLIPSIvvnj/ RM2(0:2),RMG(0:2),RM2MIN(0:2),RM2MAX(0:2),
     1                S,M2MIN,YJMIN(N1),YJMAX(N1),PTJMIN(N1),EJMIN(N1),
     2                INFOJ(N1), iwidth(0:2)
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
c definitions for second Higgs
      double precision BH2WW,BH2ZZ,BH2GG,BH2TT,BH2BB,BH2CC,BH2TAU,BH2MU,
     &                 BH2GAM, BH2GAMZ, XMH2, XGH2, 
     &                 sinba, cosba, ch2ww, ch2zz
      COMMON /BRANCH2/ BH2WW,BH2ZZ,BH2GG,BH2TT,BH2BB,BH2CC,BH2TAU,BH2MU,
     &                 BH2GAM, BH2GAMZ, XMH2, XGH2, 
     &                 sinba, cosba, ch2ww, ch2zz
      
      DOUBLE PRECISION xg ! parameter for simulation of VV distribution
      integer bos_var
      
      common /bosin/ xg,bos_var,bos

      double precision M2_i(max_res),MG_i(max_res),F_i(max_res)         ! for resonances of intermediate particles
      double precision base_M2,base_MG
      integer N_res
      
c only relevant for spin-2 scenario
      double precision xmassx, xwidthx      !mass and width of spin-2 resonance
c
      s = ecm**2
C Choose cutoff on diparton invariant mass and Q^2 for initial-final 
C collinear cutoff: for numerical purposes only
      m2min = 0.1d0**2

      bos(0) = bos0
      bos(1) = bos1
      bos(2) = bos2
      
      bos_var = bos0      
      xg = 0d0	! for initialization only (overruled for bos0=-6)
      
c set iwidth(i) = 0 to generate distributions in zero width approximation.
c the amplitude squared should have the normal Breit-Wigner factor included, 
c which then gets evaluted at q^2 = m^2, i.e. the program expects a factor 
c 1/(mGamma)^2 * |M_decay|^2 from the matrix element squared
      iwidth(0) = iwidth0
      iwidth(1) = iwidth1
      iwidth(2) = iwidth2
C set resonances: 
c       bos0 = 0    ==>  broad V1 + V2 mass distribution
c       bos0 = 6    ==>  sharp Higgs resonance for H-->V1 V2 decay
c       bos0 = -6   ==>  sharp Higgs resonance for H-->V1 V2 decay plus broad distribution
c       bos0 = 39   ==>  sharp Spin-2 resonance for Spin-2-->V1 V2 decay


! settings for intermediate (higgs / spin2 / W / Z / continuum) resonances

      rm2min(0) = 1d-3              ! min/max values may be changed in each event
      rm2max(0) = s/2d0
      rm2(0) = (sqrt(xm2(bos(1)))+sqrt(xm2(bos(2))))**2     ! intermediate particle, continuum
      rmg(0) = sqrt(rm2(0))*500d0
      N_res=0
      if (bos0.eq.6) then
         rm2(0) = xm2(6)            ! no continuum, Higgs resonance is base resonance
         rmg(0) = xmg(6)
         rm2min(0) = sqrt(rm2(0))*(1d0-150*rmg(0)/rm2(0))
         rm2min(0) = max( rm2min(0), 1d-3 )**2
         rm2max(0) = min(rm2(0)*(1d0+150*rmg(0)/rm2(0))**2,s/2d0)
      elseif (bos0.eq.39) then
         if ((abs(f1)+abs(f2)+abs(f3)+abs(f4)+abs(f5)).GE.
     &       (abs(f6)+abs(f7))) then
           xmassx=sp2mass
           xwidthx=sp2width
         else
           xmassx=msp2tripn
           xwidthx=sp2tripnwidth
         endif
         rm2(0) = xmassx**2         ! no continuum, Spin2 resonance is base resonance
         rmg(0) = xmassx*xwidthx
         rm2min(0) = sqrt(rm2(0))*(1d0-150*rmg(0)/rm2(0))
         rm2min(0) = max( rm2min(0), 1d-3 )**2
         rm2max(0) = min(rm2(0)*(1d0+150*rmg(0)/rm2(0))**2,s/2d0)
      elseif (bos0.ne.0 .and. bos0.ne.39 .and. abs(bos0).ne.6) then
         print*," bos0 = ",bos0," not implemented in ps_vvnj.f "
         stop
      endif
      base_M2 = rm2(0)
      base_MG = rmg(0)
      if (bos0.eq.-6) then
        N_res = N_res + 1     ! add higgs resonance to continuum
        M2_i(N_res) = xm2(6)
        MG_i(N_res) = xmg(6)
        F_i(N_res)  = 0.5d0
      endif
      if ( (bos0.eq.-6 .or. bos0.eq.0) .and. (bos1.eq.2 .and. bos2.eq.2) ) then   ! ZZ
        N_res = N_res + 1     ! add Z->4l resonance to continuum
        M2_i(N_res) = xm2(2)
        MG_i(N_res) = xmg(2)
        F_i(N_res)  = 0.25d0
      endif
      if ( (bos0.eq.-6 .or. bos0.eq.0) .and. ( ((bos1.eq.4 .or. bos1.eq.3) .and. bos2.eq.2) 
     &            .or. ((bos2.eq.4 .or. bos2.eq.3) .and. bos1.eq.2) ) ) then     ! WZ
        N_res = N_res + 1     ! add W->3l resonance to continuum
        M2_i(N_res) = xm2(3)
        MG_i(N_res) = xmg(3)
        F_i(N_res)  = 0.25d0
      endif
      if (model.eq.3) then
        N_res = N_res + 1     ! add W->3l resonance to continuum
        M2_i(N_res) = xmh2**2
        MG_i(N_res) = xmh2*xgh2
        F_i(N_res)  = 0.5d0
      endif

      if (with_spin2) then
         Select Case (procID)
         Case(ZZjj_lnu, ZZjj_ll, WPWMjj, WPZjj, WMZjj)
           if ((abs(f1)+abs(f2)+abs(f3)+abs(f4)+abs(f5)).GE.
     &       (abs(f6)+abs(f7))) then
             xmassx=sp2mass
             xwidthx=sp2width
           else
             Select Case (procID)
             Case(ZZjj_lnu, ZZjj_ll, WPWMjj)
               xmassx=msp2tripn
               xwidthx=sp2tripnwidth
             Case(WPZjj, WMZjj)
               xmassx=msp2trippm
               xwidthx=sp2trippmwidth
             END SELECT
           endif
           N_res = N_res + 1     ! add spin-2 resonance to continuum
           M2_i(N_res) = xmassx**2
           MG_i(N_res) = xmassx*xwidthx
           F_i(N_res)  = 0.5d0
         end select
      endif ! with_spin2

      call MultiResonancesInit(1, base_m2, base_mg, 1d0, M2_i, MG_i, F_i, N_res)


! settings for vector bosons

      do i = 1,2
         rm2(i) = xm2(abs(bos(i)))
         rmg(i) = xmg(abs(bos(i)))
         rm2min(i) = 1d-3
         rm2max(i) = s/2d0
      enddo


      do i = 1,N1
         yjmin(i) = 0d0
         yjmax(i) = ymax_def
         ptjmin(i) = max(ptmin_def(i),10d0)
         ejmin(i) = 0d0
         infoj(i) = -1
      enddo
 60   format ( " machine energy = ", f10.1, "  GeV" )
 70   format ( " integration range for boson ",i2," = ",f10.3,
     1         "  to ",f12.3," GeV" )
 72   format ( " boson ",i2," has mass ",f8.2," GeV and width ",
     1         f10.4," GeV" )
 80   format ( " Minimal mom.transfer or mass for parton pair: ",f8.2,
     1         " MeV" )
      if (.not.ldoblha) then
        write(6,*)
        write(6,*) "       phase space generator information  "
        write(6,*) "-----------------------------------------------"
        write(6,*)
        do i = 0,2
           write(6,72) bos(i),sqrt(rm2(i)),rmg(i)/sqrt(rm2(i))
           write(6,70) bos(i),sqrt(rm2min(i)),sqrt(rm2max(i))
        enddo
        write(6,*)
c        write(6,10) ptjmin
c        write(6,20) ejmin
c        write(6,30) yjmin
c        write(6,40) yjmax
c        write(6,50) infoj
        write(6,60) sqrt(s)
        write(6,80) sqrt(m2min)*1000.d0
        write(6,*)
        write(6,*) "-----------------------------------------------"
        write(6,*)
      endif
      
      end
C
C**********************  LIPSNVV   ****************************************
C
      SUBROUTINE LIPSNVV(N,RD,RN,
     1                   K1,K2,D,P,
     2                   X1,X2,W)
      IMPLICIT NONE
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/spin2coupl.inc"
#include "VBFNLO/utilities/process.inc"
c
      INTEGER N,NRANDLG,NRANDL,N1
      DOUBLE PRECISION TPI
      PARAMETER (N1=max_p-nx, NRANDLG=N1*3, nrandl=(N1+3)*3, TPI=2D0*PI )
      DOUBLE PRECISION R(NRANDL),RD(*),RN,
     1       K1(0:3),K2(0:3),Q(0:4,2),D(0:3,max_v),P(0:3,N1),
     2       X1,X2,W
C
C Phase space generator for production of boson + N jets at hadron colliders
C
C  INPUT:   RD        random number array (between 0 and 1)
C           RN        additional random number which will not be driven by
C                     VEGAS
C
C  OUTPUT:  K1,K2     4-vectors for 2 incoming partons
C           X1,X2     corresponding values of Feynman x
C           Q(0:4,2)  4-vectors for the vector boson resonances V1 and V2,
C                     the 4th component contains the mass**2
C           D(0:3,4)  4-vectors for the 4 decay leptons, assumed massless
C                     V1(Q(*,1)) --> D(*,1)+D(*,2)
C                     V2(Q(*,2)) --> D(*,3)+D(*,4)
C           P(*,N)    4-vectors for the N jets, assumed massless
C           W         the weight of the phase space point
C
C  For a final state photon the ENTRY LIPSV has to be used, which generates a
C  massless four-vector Q
C
C  INPUT INFORMATION on the phase space region is transferred via the 
C    common block /BLIPSIvvnj/
C

      DOUBLE PRECISION RM2,RMG,RM2MIN,RM2MAX, S,M2MIN, YJMIN,YJMAX,PTJMIN,EJMIN
      INTEGER INFOJ, iwidth
C           S         center of mass energy**2 of the machine
      COMMON /BLIPSIvvnj/ RM2(0:2),RMG(0:2),RM2MIN(0:2),RM2MAX(0:2),
     1                S,M2MIN,YJMIN(N1),YJMAX(N1),PTJMIN(N1),EJMIN(N1),
     2                INFOJ(N1), iwidth(0:2)
C
C   RM2            resonance mass squared for resonance
C   RMG            resonance total width * mass for resonance 
C   RM2_MIN;_MAX   min & max values of the inv. mass**2 of resonance
C   YJMIN;YJMAX    rapidity range of the individual jets:
C                  YJMIN(I) < ABS(Y(I)) < YJMAX(I)
C   INFOJ(I)=1     PT of jet I is integration variable
C           =2     energy of jet I is integration variable
C   PTJMIN         minimal PT values
C   EJMIN          minimal energy values of individual jets
C
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      
      DOUBLE PRECISION xg
      integer bos0, bos(0:2) 		
      common /bosin/ xg,bos0,bos	! info on VV decay 

C  OUTPUT INFORMATION is transferred via the common block /BLIPSO/
C
      DOUBLE PRECISION Y,PT,PHI,PNJ
      COMMON /BLIPSO/ Y(N1), PT(N1), PHI(N1), PNJ(0:4)
C
C  local variables
      DOUBLE PRECISION RS,FAC,DELY,XMIN,XMAX,X,DELX,COSHY,SINHY,M2,dotrr,
     1       YCM,YMIN,YMAX, SUMPST,Q0ST,RSHAT,PCM(0:3),QL(0:4,0:2),
     2       STH,CTH,PHID, XMJ,QV,YBAR,PJT2,SHA,ACM, 
     3       xmass, xmvv, ev1, pv1
     
      DOUBLE PRECISION rx

      logical resonance, nonresonance, resonance2, multiresonance

      external dotrr,xmass, resonance, nonresonance, resonance2, multiresonance
      INTEGER I,J,II,JJ,MU,ic
      data ic /1/
      save ic
      logical ldebug
      parameter (ldebug=.false.)
      

c only relevant for spin-2 scenario
      double precision xwidthx   !mass and width of the spin-2 resonance
      double precision hw             !mass and width of the higgs resonance
      double precision ovp2         !OverLapPoint of spin-2 and higgs resonance for two cases

    
      if (ldebug) print*," enter LIPSNVV"
      DO I = 1,(N+3)*3
         R(I) = RD(I)
      ENDDO

      FAC = 1/(16*PI**3)
      RS = SQRT(S)
      W = TPI
C  Reset local variables and incoming momenta
      DO 10 I = 0,3
         PCM(I) = 0
         PNJ(I) = 0
         K1(I) = 0
         K2(I) = 0
10    CONTINUE
      PNJ(4) = 0

C
C  Get 4-momenta of the N massless jets
      DO 20 I = 1,N
         W = W*FAC
         J = 3*(I+2)
C
C  generate p_T first for INFOJ=1. Use p_T * coshy = E to limit 
C  the rapidity range in the second step
         IF(INFOJ(I).EQ.1) THEN
            XMIN = 2D0/RS
            XMAX = 1/PTJMIN(I)
            DELX = XMAX-XMIN
            X = XMIN + R(J+2)*DELX
            PT(I) = 1/X
            YMAX = RS/(2*PT(I))
            YMAX = LOG(YMAX+SQRT(YMAX**2-1))
            W = W*DELX*PT(I)**3
         ELSEIF(INFOJ(I).EQ.2) THEN
            XMIN = 2.001D0/RS
            XMAX = 1/EJMIN(I)
            DELX = XMAX-XMIN
            X = XMIN + R(J+2)*DELX
            P(0,I) = 1/X
            IF(PTJMIN(I).GT.0) THEN
              YMAX = P(0,I)/PTJMIN(I)
              YMAX = LOG(YMAX+SQRT(YMAX**2-1))
            ELSE 
              YMAX = 1000
            ENDIF
         elseif(infoj(i).eq.-1) then
            if (r(j+2).gt.0.2) then
               XMIN = 2D0/RS
               XMAX = 1/PTJMIN(I)
               DELX = XMAX-XMIN
               X = XMIN + (R(J+2)-0.2)/0.8d0*DELX
               PT(I) = 1/X
               W = W/0.8d0*DELX*PT(I)**3
            else
               xmin = 0d0
               xmax = ptjmin(i)
               DELX = XMAX-XMIN
               PT(i) = XMIN + R(J+2)/0.2d0*DELX
               W = W/0.2d0*DELX*PT(I)
            endif
            YMAX = RS/(2*PT(I))
            YMAX = LOG(YMAX+SQRT(YMAX**2-1))
         ELSE
            W = 0
            RETURN
         ENDIF

c         YMAX = MIN(YMAX,YJMAX(I))
         IF (YMAX.LE.YJMIN(I)) THEN
            W = 0
            RETURN
         ENDIF
         DELY = 2*(YMAX-YJMIN(I))
C         IF (I.EQ.1) THEN
         IF (I.LE.3) THEN
            Y(I) = -YMAX+R(J+1)*DELY
            IF (Y(I).GT.-YJMIN(I)) Y(I) = 2*YJMIN(I) + Y(I)
         ELSE
C  fix random number to 0 or 1 for delta-y = 0 with previous jet: R=0 or 1 
C  corresponds to lowest dijet invariant mass
            IF (Y(I-1).LT.0D0) THEN
               IF (Y(I-1).LT.-YMAX) THEN
                  Y(I) = -YMAX + R(J+1)*DELY
                  IF (Y(I).GT.-YJMIN(I)) Y(I) = 2*YJMIN(I) + Y(I)
               ELSEIF (Y(I-1).LT.-YJMIN(I)) THEN
                  Y(I) = Y(I-1) + R(J+1)*DELY
                  IF (Y(I).GT.-YJMIN(I)) Y(I) = 2*YJMIN(I) + Y(I)
                  IF (Y(I).GT.YMAX) Y(I) = -2*YMAX + Y(I)
               ELSE
                  Y(I) = -YJMIN(I) - R(J+1)*DELY
                  IF (Y(I).LT.-YMAX) Y(I) = 2*YMAX + Y(I)
               ENDIF
            ELSE
               IF (Y(I-1).GT.YMAX) THEN
                  Y(I) = YMAX - R(J+1)*DELY
                  IF (Y(I).LT.YJMIN(I)) Y(I) = -2*YJMIN(I) + Y(I)
               ELSEIF (Y(I-1).GT.YJMIN(I)) THEN
                  Y(I) = Y(I-1) - R(J+1)*DELY
                  IF (Y(I).LT.YJMIN(I)) Y(I) = -2*YJMIN(I) + Y(I)
                  IF (Y(I).LT.-YMAX) Y(I) = 2*YMAX + Y(I)
               ELSE
                  Y(I) = YJMIN(I) + R(J+1)*DELY
                  IF (Y(I).GT.YMAX) Y(I) = -2*YMAX + Y(I)
               ENDIF
            ENDIF
         ENDIF
         SINHY = SINH(Y(I))
         COSHY = SQRT(1+SINHY**2)
         W = W*DELY

         IF(abs(INFOJ(I)).EQ.1) THEN
            P(0,I) = PT(I)*COSHY
         ELSEIF(INFOJ(I).EQ.2) THEN
            PT(I) = P(0,I)/COSHY
            W = W*DELX*P(0,I)**3/COSHY**2
         ENDIF

         IF (I.EQ.1) THEN
            PHI(I) = PI*(2*RN-1)
         ELSE
            PHI(I) = PHI(I-1) + TPI*R(J+3)
            if (phi(i).gt.pi) phi(i) = phi(i)-tpi
         ENDIF

         W = W*TPI

         P(1,I) = PT(I)*COS(PHI(I))
         P(2,I) = PT(I)*SIN(PHI(I))
         P(3,I) = PT(I)*SINHY
         if (ldebug) then
            print*," P_",i," = ",(p(mu,i),mu=0,3)
         endif
         DO 30 J = 0,3
            PNJ(J) = PNJ(J)+P(J,I)
30       CONTINUE
20    CONTINUE
C  determine N-jet invariant mass
      IF (N.GT.1) THEN
         PNJ(4) = PNJ(0)**2-PNJ(1)**2-PNJ(2)**2-PNJ(3)**2
      ELSE
         PNJ(4) = 1D-18
      ENDIF
      IF (PNJ(4).GE.S) THEN
         W=0
         RETURN
      ENDIF

c  check invariant masses of all parton pairs: eliminate events which 
c  may cause numerical problems 
      do ii = 1,N-1
         do jj = ii+1,N
            m2 = 2*dotrr(p(0,ii),p(0,jj))
            if (m2.lt.m2min) then
               W = 0
               return
            endif
         enddo
      enddo

c-------------------------------------------------------------------
c
C  Generate Q**2 = Q(4) according to Breit Wigner distribution
      ic = 3-ic      ! oscillate between ic = 1 and ic = 2
      if (ldebug) print*, " ic = ",ic
      do ii = 0,2
         i = ii*(2-ii)*ic + (ii/2)*(3-ic) ! = 0 for ii=0, ic for ii=1, 3-ic for ii=2
         if (ldebug) print*," ic, ii,i = ",ic,ii,i 
         if (iwidth(i).eq.0) then
            ql(4,i) = rm2(i)
            w = w*rmg(i)/2
            goto 40
         endif
         if (ii.eq.1) then
            rm2max(i) = ql(4,0)
         elseif (ii.eq.2) then
            rm2max(i) = (sqrt(ql(4,0))-sqrt(ql(4,ic)))**2
            if ( bos0.ne.0 .or. (bos(1).eq.bos(2)) ) then      ! W+W-, W+W+, W-W-, ZZ
              if (rm2min(1).eq.rm2min(2)) then
                 rm2max(i) = min( rm2max(i), ql(4,ic) )
                 w = w*2d0
              endif
            endif
         endif

         if ( bos0.ne.0 .or. (bos(1).eq.bos(2)) ) then      ! W+W-, W+W+, W-W-, ZZ production
           rx = R(1+3*ii)
         else                                               ! WZ production
           rx = R(1+3*i)
         endif
         if ( i.gt.0 ) then
         ! possible different choices here:
         ! * increase nonresonance-part of lighter V to 0.5:
         !   works better for ZZ -> lnu, but less efficient for other processes
         ! * check for "RM2(i)-8d0*RMG(i)", else only nonresonant part:
         !   improves total number of events, but increases runtime even more, due to more events passing PS-generator,
         !   but not mapped optimal on random number
            if (                                                                    ! split rnd for Q2 of bosons in the following cases:
     &            bos0.eq.-6 .or.                                                   ! WW / ZZ production (higgs (+Sp2) + continuum)
     &            ( (bos0.eq.6 .or. bos0.eq.39) .and. bos(1).eq.2 .and. ii.eq.1 )   ! (light) Higgs / Spin2 -> ZZ, heavier Z
     &         ) then
               if (ii.eq.1) then                ! heavier V
                  if (rx .lt. 0.3d0) then
                     W = W / 0.3d0
                     if (.not.NonResonance(RM2MIN(i), min(RM2(i)-8d0*RMG(i),RM2MAX(i)), rx/0.3d0, W, Ql(4,i))) RETURN 
                  else
                     W = W / 0.7d0
                     if (.not.Resonance2(RM2(i), RMG(i), RM2(i)-8d0*RMG(i), RM2MAX(i), 8d0, (rx-0.3d0)/0.7d0, W, Ql(4,i))) RETURN
                  endif
               else                             ! lighter V
                  if (rx .lt. 0.3d0) then
                     W = W / 0.3d0
                     if (.not.NonResonance(RM2MIN(i), min(RM2(i)-8d0*RMG(i),RM2MAX(i)), rx/0.3d0, W, Ql(4,i))) RETURN 
                  else
                     W = W / 0.7d0
                     if (.not.Resonance2(RM2(i), RMG(i), RM2(i)-8d0*RMG(i), RM2MAX(i), 8d0, (rx-0.3d0)/0.7d0, W, Ql(4,i))) RETURN
                  endif
               endif
            elseif (bos0.eq.0 .and. .not.(bos(1).eq.bos(2))) then                   ! WZ production
               if (bos(i).eq.2) then
                  if (rx .lt. 0.3d0) then
                     W = W / 0.3d0
                     if (.not.NonResonance(RM2MIN(i), min(RM2(i)-8d0*RMG(i),RM2MAX(i)), rx/0.3d0, W, Ql(4,i))) RETURN 
                  else
                     W = W / 0.7d0
                     if (.not.Resonance2(RM2(i), RMG(i), RM2(i)-8d0*RMG(i), RM2MAX(i), 8d0, (rx-0.3d0)/0.7d0, W, Ql(4,i))) RETURN
                  endif
               else
                 if (.not.Resonance2(RM2(i), RMG(i), RM2MIN(i), RM2MAX(i), 8d0, rx, W, Ql(4,i))) RETURN
               endif
            else                                                                    ! W+W+, W-W-, lighter Z from Higgs/Spin2
               if (.not.Resonance(RM2(i), RMG(i), RM2MIN(i), RM2MAX(i), rx, W, Ql(4,i))) RETURN
            endif
         else
            if (.not. MultiResonance(1, rx, rm2min(0), rm2max(0), W, Ql(4,0))) return
         endif

 40      continue
      enddo
c
c-------------------------------------------------------------------

      if (ldebug) then
         print*," mass of boson 0 ",sqrt(ql(4,0))
         print*," mass of boson 1 ",sqrt(ql(4,1))
         print*," mass of boson 2 ",sqrt(ql(4,2))
      endif
      x = (sqrt(ql(4,1))+sqrt(ql(4,2)))**2
      if (x.ge.ql(4,0)) then
         w = 0
         if (ldebug) print*," mass of boson 0 too small "
         return
      endif
      
      
C  Next get the center of mass rapidity YCM
      QV = SQRT(Ql(4,0))
      XMJ = SQRT(PNJ(4))
      YBAR = 0.5D0*LOG( (PNJ(0)+PNJ(3))/(PNJ(0)-PNJ(3)) )
      PJT2 = PNJ(1)**2+PNJ(2)**2
      SHA = (S-(QV+XMJ)**2-4*PJT2) / (4*(PNJ(4)+PJT2)) 
      IF (SHA.LE.0D0) THEN
        W=0
         RETURN
      ENDIF
      SHA = SQRT(SHA)
      ACM = LOG( SHA+SQRT(SHA**2+1) ) 
      YMIN = YBAR - ACM
      YMAX = YBAR + ACM

      DELY = YMAX-YMIN
      YCM = YMIN + R(12)*DELY
      SINHY = SINH(YCM)
      COSHY = SQRT(1+SINHY**2)
      W = W*DELY

C
C  now get Q and SHAT etc.
      SUMPST = PNJ(1)**2+PNJ(2)**2+(PNJ(3)*COSHY-PNJ(0)*SINHY)**2
      Q0ST = SQRT(Ql(4,0)+SUMPST)
      RSHAT = Q0ST + SQRT(PNJ(4)+SUMPST)
      PCM(0) = RSHAT*COSHY
      PCM(3) = RSHAT*SINHY

      DO mu = 0,3
        Ql(mu,0) = PCM(mu)-PNJ(mu)
      ENDDO

      X1 = (PCM(0)+PCM(3))/RS
      X2 = (PCM(0)-PCM(3))/RS
      IF (X1.GE.1.0D0 .OR. X2.GE.1.0D0) THEN
         W=0
         RETURN
      ENDIF
      K1(0) =  X1*RS/2
      K1(3) =  K1(0)
      K2(0) =  X2*RS/2
      K2(3) = -K2(0)
c  check initial-final momentum transfer for all partons: eliminate events 
c  which may cause numerical problems 
      do ii = 1,N         
         m2 = 2*min( p(0,ii)*k1(0)-p(3,ii)*k1(3),
     1               p(0,ii)*k2(0)-p(3,ii)*k2(3) )
         if (m2.lt.m2min) then
            W = 0
            return
         endif
      enddo

      QV = Ql(0,0)**2-Ql(1,0)**2-Ql(2,0)**2-Ql(3,0)**2
      IF(ldebug .and. ABS(QV/Ql(4,0)-1D0).GT.1D-6) THEN
         PRINT*," "
         PRINT*," MASS OF V WRONG: ",QV, qv/ql(4,0)
         PRINT*," Q = ",(Ql(i,0),i=0,4)
         PRINT*," "
      ENDIF

C  insert Jacobian for dtau * delta(PCM(0)-Ql(0,0)-PNJ(0))/(2 Ql(0,0))
      W = W * RSHAT/(S*Q0ST)

C
C now calculate the decay of the VV cluster, given by Ql(*,0), into 4 leptons
C first the decay cluster --> VV
      CTH = -1+2*R(2)
      STH = SQRT(1-CTH**2)
      PHID = PHI(1) + TPI*R(3)
      xmvv = sqrt(ql(4,0))
      ev1 = (ql(4,0)-ql(4,2)+ql(4,1))/(2*xmvv)
      if ((ev1**2-ql(4,1)) .le. 1D-9) then
         w = 0D0
         return
      end if
      pv1 = sqrt(ev1**2-ql(4,1))
      QL(0,1) = ev1
      QL(1,1) = pv1*STH*COS(PHID)
      QL(2,1) = pv1*STH*SIN(PHID)
      QL(3,1) = pv1*CTH
c boost q1 from the q0 rest frame to the lab frame and get q2 from mom.cons.
      call boostn(ql(0,1),ql(0,0),q(0,1))
      do mu = 0,3
         q(mu,2) = ql(mu,0) - q(mu,1)
      enddo
      q(4,1) = ql(4,1)
      q(4,2) = ql(4,2)
c and multiply weight factor by beta/8 pi
      W = W * 2*pv1/(xmvv*8*PI) 
C  finally let V1,V2 decay into 2 massless leptons each
C  generate D1,D3 in V rest frame
      do i = 1,2
         j = i*3
         ii = 2*i-1
         CTH = -1+2*R(j+2)
         if ((1D0-CTH**2) .le. 1D-8) then
            W = 0D0
            return
         end if
         STH = SQRT(1-CTH**2)
         PHID = PHI(1) + TPI*R(j+3)
         D(0,ii) = SQRT(Q(4,i))/2
         D(1,ii) = D(0,ii)*STH*COS(PHID)
         D(2,ii) = D(0,ii)*STH*SIN(PHID)
         D(3,ii) = D(0,ii)*CTH
C
C  boost D1 to its lab frame value and get D2 from momentum conservation
         CALL BOOSTN(D(0,ii),Q(0,i),D(0,ii))  
         do mu = 0,3
            D(mu,ii+1) = Q(mu,i)-D(mu,ii)
         enddo
C  phase space factor for the decay
         W = W /(8*PI)                        ! *2*TPI /(32*PI**2)

         if (ldebug) then
c  check mass of Vi and 2 decay bosons
            print 301," mass of boson ",i,sqrt(q(4,i)),
     1             (q(0,i)**2-q(1,i)**2-q(2,i)**2-q(3,i)**2)/q(4,i)
            print 301," mass of lepton ",ii,xmass(d(0,ii)),
     1             xmass(d(0,ii))/d(0,ii)
            print 301," mass of lepton ",ii+1,xmass(d(0,ii+1)),
     1             xmass(d(0,ii+1))/d(0,ii+1)
         endif
 301     format(a,i5,2f15.8)

      enddo
C
C  finally put in flux factor for the parton collisions and conversion to fb
      W = W / (2*RSHAT**2) * 3.89379304D+11

c check momentum conservation
      if (ldebug .or. lwarn) then
         do mu=0,3
            qv = k1(mu)+k2(mu)
            do i = 1,N
               qv = qv - p(mu,i)
            enddo
            do i = 1,4
               qv = qv - d(mu,i)
            enddo
            if (ldebug) then
               print 301," mom.cons. for mu = ",mu,qv,qv/(k1(0)+k2(0))
            endif
            if (lwarn .and. abs(qv).gt.1d-10*(k1(0)+k2(0)) ) then
               print*," WARNING from LIPSNVV "
               print 301," no mom.cons. for mu = ",
     1                    mu,qv,qv/(k1(0)+k2(0))
            endif
         enddo
      endif

      END

c-----------------------------------------------------------

      DOUBLE PRECISION FUNCTION xmass(p)

      IMPLICIT NONE

      DOUBLE PRECISION p(0:3), xm


      xm = p(0)**2-p(1)**2-p(2)**2-p(3)**2
      if (xm.gt.0) then
         xmass = sqrt(xm)
      elseif (xm.lt.0) then
         xmass = -sqrt(-xm)
      else
         xmass = 0
      endif
      end
