c
c phase space generator for pp --> V + (np-nx) jets, V --> l1 l2
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  1995 March 25
c	Last modified:  2003 March 20
c
      subroutine def_ps_vnj(bos)
       use globalvars, only: lglobalprint, ldoblha
      implicit none 
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/cuts.inc"
#include "VBFNLO/utilities/spin2coupl.inc"
#include "VBFNLO/utilities/process.inc"
      integer*4 N1, bos, i 
      parameter ( N1=max_p-nx )
c
c set up variables in common block BLIPSIvnj for call of LIPSN
c 
      DOUBLE PRECISION RM2,RMG,RM2MIN,RM2MAX, S,M2MIN, YJMIN,YJMAX,PTJMIN,EJMIN
      INTEGER INFOJ
      COMMON /BLIPSIvnj/ RM2,RMG,RM2MIN,RM2MAX,S,M2MIN,
     1    YJMIN(N1),YJMAX(N1),PTJMIN(N1),EJMIN(N1),INFOJ(N1)

      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

c only relevant for spin-2 scenario
      double precision xmassx, xwidthx      !mass and width of spin-2 resonance
c
      s = ecm**2
C Choose cutoff on diparton invariant mass and Q^2 for initial-final 
C collinear cutoff: for numerical purposes only
      m2min = 0.1d0**2
      if (bos.eq.39) then   !only relevant for spin-2 scenario
        if ((abs(f1)+abs(f2)+abs(f3)+abs(f4)+abs(f5)).GE.
     &       (abs(f6)+abs(f7))) then
           xmassx=sp2mass
           xwidthx=sp2width
        else
           xmassx=msp2tripn
           xwidthx=sp2tripnwidth
        endif
        rm2 = xmassx**2
        rmg = xmassx*xwidthx
      else
         rm2 = xm2(bos)
         rmg = xmg(bos)
      endif
      if (bos.eq.1) then 
         rm2min = 1d-3
         rm2max = s/2d0
      elseif (bos.eq.3.or.bos.eq.4)then
             rm2min = 1d-3     
             rm2max = s/2d0 
      elseif(bos.eq.6) then
        SELECT CASE(procID)
          CASE(AAjj)     !pp-> higgs instead of spin-2 jj -> aa jj
           rm2min = 1d-3
           rm2max = s/2d0
          CASE(Hjj,Hjjj)
            if (higgsscheme .eq. 0) then
              rm2min = sqrt(rm2)*( 1.d0 - 15.d1*rmg/rm2 )
              if ( rm2min.lt.0.d0 ) then
                 rm2min = 1d-3
              else
                 rm2min = rm2min**2
              end if
              rm2max = rm2*( 1.d0 + 15.d1*rmg/rm2 )**2
              rm2max = min(rm2max,(sqrt(rm2) + 200)**2)
            else
              rm2min = 1d-3
              rm2max = s/2d0
            endif
          CASE DEFAULT
           rm2min = sqrt(rm2)*( 1.d0 - 15.d1*rmg/rm2 )
           if ( rm2min.lt.0.d0 ) then
              rm2min = 1d-3
           else
              rm2min = rm2min**2
           end if
           rm2max = rm2*( 1.d0 + 15.d1*rmg/rm2 )**2
           rm2max = min(rm2max,(sqrt(rm2) + 200)**2)
        END SELECT
      else
           rm2min = 1d-3
           rm2max = s/2d0
      endif
      do i = 1,N1
         yjmin(i) = 0d0
         yjmax(i) = ymax_def
         ptjmin(i) = max(ptmin_def(i),10d0)
         ejmin(i) = 00.0
         infoj(i) = -1
      enddo

 10   format ( " pTjmin = ", 5(f6.1,2x) )
 20   format ( " Ejmin  = ", 5(f6.1,2x) )
 30   format ( " yjmin  = ", 5(f6.1,2x) )
 40   format ( " yjmax  = ", 5(f6.1,2x) )
 50   format ( " infoj  = ", 5(i6,2x) )
 60   format ( " machine energy = ", f10.1, "  GeV" )
 70   format ( " integration range for boson ",i2," = ",f10.3,
     1         "  to ",f12.3," GeV" )
 72   format ( " boson ",i2," has mass ",f6.2," GeV and width ",
     1         f7.4," GeV" )
 73   format ( " boson ",i2," has mass ",f6.2," GeV" )
 80   format ( " Minimal mom.transfer or mass for parton pair: ",f8.2,
     1         " MeV" )
      if (.not.ldoblha .and. lglobalprint) then
        write(6,*)
        write(6,*) "       phase space generator information  "
        write(6,*) "-----------------------------------------------"
        write(6,*)
        if (sqrt(rm2).ne.0d0) then 
           write(6,72) bos,sqrt(rm2),rmg/sqrt(rm2)
        else
           write(6,73) bos,sqrt(rm2)
        endif
        write(6,70) bos,sqrt(rm2min),sqrt(rm2max)
        write(6,*)
        write(6,10) ptjmin
        write(6,20) ejmin
        write(6,30) yjmin
        write(6,40) yjmax
        write(6,50) infoj
        write(6,60) sqrt(s)
        write(6,80) sqrt(m2min)*1000.d0
        write(6,*)
        write(6,*) "-----------------------------------------------"
        write(6,*)
      endif
      
      end

C
C**********************  LIPSN0   ****************************************
C
      SUBROUTINE LIPSN0(N,RD,RN,
     1                 K1,K2,Q,P,
     2                 X1,X2,W)
      IMPLICIT NONE
#include "VBFNLO/utilities//global.inc"
c
      INTEGER N,NRANDLG,NRANDL,N1
      DOUBLE PRECISION TPI
      PARAMETER (N1=max_p-nx, NRANDLG=N1*3, nrandl=(N1+1)*3, TPI=2D0*PI )
      DOUBLE PRECISION R(NRANDL),RD(*),RN,
     1       K1(0:3),K2(0:3),Q(0:4),P(0:3,N1),
     2       X1,X2,W
C
C Phase space generator for production of boson + N jets at hadron colliders
C
C  INPUT:   RD,RG     random number array (between 0 and 1)
C           RN        additional random number which will not be driven by
C                     VEGAS
C
C  OUTPUT:  K1,K2     4-vectors for 2 incoming partons
C           X1,X2     corresponding values of Feynman x
C           Q(0:4)    4-vector for the vector boson resonance V, the 4th
C                     component contains the mass**2
C           P(*,N)    4-vectors for the N jets, assumed massless
C           W         the weight of the phase space point
C
C  INPUT INFORMATION on the phase space region is transferred via the 
C    common block /BLIPSIvnj/
C

      DOUBLE PRECISION RM2,RMG,RM2MIN,RM2MAX, S,M2MIN, YJMIN,YJMAX,PTJMIN,EJMIN
      INTEGER INFOJ
C           S         center of mass energy**2 of the machine
      COMMON /BLIPSIvnj/ RM2,RMG,RM2MIN,RM2MAX,S,M2MIN,
     1    YJMIN(N1),YJMAX(N1),PTJMIN(N1),EJMIN(N1),INFOJ(N1)
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
C  OUTPUT INFORMATION is transferred via the common block /BLIPSO/
C
      DOUBLE PRECISION Y,PT,PHI,PNJ
      COMMON /BLIPSO/ Y(N1), PT(N1), PHI(N1), PNJ(0:4)
C
C  local variables
      DOUBLE PRECISION RS,FAC,DELY,XMIN,XMAX,X,DELX,COSHY,SINHY,M2,dotrr,
     1       YCM,YMIN,YMAX, SUMPST,Q0ST,RSHAT,PCM(0:3),
     2       XMJ,QV,YBAR,PJT2,SHA,ACM
      external dotrr
      INTEGER I,J,II,JJ
      DO 1 I = 1,N*3
         R(I) = RD(I)
 1    CONTINUE
      
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
         J = 3*(I-1)
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
               XMAX = 1/max(PTJMIN(I),10d0)
               DELX = XMAX-XMIN
               X = XMIN + (R(J+2)-0.2)/0.8d0*DELX
               PT(I) = 1/X
               W = W/0.8d0*DELX*PT(I)**3
            else
               xmin = 0
               xmax = max(ptjmin(i),10d0)
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
            PHI(I) = PHI(I-1) + TPI*(R(J+3)+ 0.25 )
            if (phi(i).gt.pi) phi(i) = phi(i)-tpi
            if (phi(i).gt.pi) phi(i) = phi(i)-tpi
         ENDIF
         
         W = W*TPI

         P(1,I) = PT(I)*COS(PHI(I))
         P(2,I) = PT(I)*SIN(PHI(I))
         P(3,I) = PT(I)*SINHY
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
    
      Q(4) = RM2
      QV = SQRT(Q(4))
    
C  Next get the center of mass rapidity YCM
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
      YCM = YMIN + R(3)*DELY !!!!!!!
      SINHY = SINH(YCM)
      COSHY = SQRT(1+SINHY**2)
      W = W*DELY

C
C  now get Q and SHAT etc.
      SUMPST = PNJ(1)**2+PNJ(2)**2+(PNJ(3)*COSHY-PNJ(0)*SINHY)**2
      Q0ST = SQRT(Q(4)+SUMPST)
      RSHAT = Q0ST + SQRT(PNJ(4)+SUMPST)
      PCM(0) = RSHAT*COSHY
      PCM(3) = RSHAT*SINHY

      DO 40 I = 0,3
        Q(I) = PCM(I)-PNJ(I)
40    CONTINUE

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

C
C  insert Jacobian for dtau * delta(PCM(0)-Q(0)-PNJ(0))/(2 Q(0))
      W = W * RSHAT/(S*Q0ST)

C
C  finally put in flux factor for the parton collisions and conversion to fb
      W = W / (2*RSHAT**2) * 3.89379304D+11

      END


C
C**********************  LIPSN   ****************************************
C
      SUBROUTINE LIPSN(N,RD,RN,
     1                 K1,K2,Q,D1,D2,P,
     2                 X1,X2,W, NW)
      use globalvars, only: ldoblha
      IMPLICIT NONE
#include "VBFNLO/utilities/global.inc"
c
      INTEGER N,NRANDLG,NRANDL,N1, NW
      DOUBLE PRECISION TPI
      PARAMETER (N1=max_p-nx, NRANDLG=N1*3, nrandl=(N1+1)*3, TPI=2D0*PI )
      DOUBLE PRECISION R(NRANDL),RD(*),RN,
     1       K1(0:3),K2(0:3),Q(0:4),D1(0:3),D2(0:3),P(0:3,N1),
     2       X1,X2,W
C
C Phase space generator for production of boson + N jets at hadron colliders
C
C  INPUT:   RD,RG     random number array (between 0 and 1)
C           RN        additional random number which will not be driven by
C                     VEGAS
C
C  OUTPUT:  K1,K2     4-vectors for 2 incoming partons
C           X1,X2     corresponding values of Feynman x
C           Q(0:4)    4-vector for the vector boson resonance V, the 4th
C                     component contains the mass**2
C           D1,D2     4-vectors for the 2 decay leptons of V, assumed massless
C           P(*,N)    4-vectors for the N jets, assumed massless
C           W         the weight of the phase space point
C           NW        Q^2 (Higgs) is generated in narrow width approx.
C                       = 0 = OS
C                       = 1 = BW
C
C  For a final state photon the ENTRY LIPSV has to be used, which generates a
C  massless four-vector Q
C
C  INPUT INFORMATION on the phase space region is transferred via the 
C    common block /BLIPSIvnj/
C

      DOUBLE PRECISION RM2,RMG,RM2MIN,RM2MAX, S,M2MIN, YJMIN,YJMAX,PTJMIN,EJMIN
      INTEGER INFOJ
C           S         center of mass energy**2 of the machine
      COMMON /BLIPSIvnj/ RM2,RMG,RM2MIN,RM2MAX,S,M2MIN,
     1    YJMIN(N1),YJMAX(N1),PTJMIN(N1),EJMIN(N1),INFOJ(N1)
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
C  OUTPUT INFORMATION is transferred via the common block /BLIPSO/
C
      DOUBLE PRECISION Y,PT,PHI,PNJ
      COMMON /BLIPSO/ Y(N1), PT(N1), PHI(N1), PNJ(0:4)
C
C  local variables
      DOUBLE PRECISION RS,FAC,DELY,XMIN,XMAX,X,DELX,COSHY,SINHY,M2,dotrr,
     1       YCM,YMIN,YMAX, SUMPST,Q0ST,RSHAT,PCM(0:3),
     2       STH,CTH,PHID, XMJ,QV,YBAR,PJT2,SHA,ACM
      external dotrr
      INTEGER I,J,II,JJ,BOS

** rescaling momentum for BW distribution when calculating EW cor'ns
      double precision sh, sinc, shin, disc, scafac
      integer mu


      DO 1 I = 1,(N+1)*3
         R(I) = RD(I)
 1    CONTINUE
      BOS = 2
c      GOTO 5
c
c      ENTRY LIPSV(RG,RN,
c     1            K1,K2,Q,P,
c     2            X1,X2,W)
c      DO 2 I = 1,N*3
c         R(I) = RG(I)
c 2    CONTINUE
c      BOS = 1
c 5    CONTINUE
      
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
         J = 3*I
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
               XMAX = 1/max(PTJMIN(I),10d0)
               DELX = XMAX-XMIN
               X = XMIN + (R(J+2)-0.2)/0.8d0*DELX
               PT(I) = 1/X
               W = W/0.8d0*DELX*PT(I)**3
            else
               xmin = 0
               xmax = max(ptjmin(i),10d0)
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
           if (ldoblha .and. NW.eq.1 .and. I.eq.N) then
             PHI(I) = PHI(I-1) + TPI*(R(1)+0.25d0)
           else
             PHI(I) = PHI(I-1) + TPI*(R(J+3)+0.25d0)
           endif
            if (phi(i).gt.pi) phi(i) = phi(i)-tpi
            if (phi(i).gt.pi) phi(i) = phi(i)-tpi
         ENDIF
         
         W = W*TPI

         P(1,I) = PT(I)*COS(PHI(I))
         P(2,I) = PT(I)*SIN(PHI(I))
         P(3,I) = PT(I)*SINHY
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

      IF (BOS.EQ.1) THEN
         Q(4) = 0
         QV = 0
      ELSE
      if (NW.eq.1) then ! narrow width approximation
         Q(4) = rm2
         QV = SQRT(Q(4))
         W = W * rmg/2.0d0
      else if (NW.eq.0) then
C  Generate Q**2 = Q(4) according to Breit Wigner distribution
         XMAX = (RM2MAX-RM2)/RMG
         XMAX = ATAN(XMAX) 
         XMIN = (RM2MIN-RM2)/RMG
         XMIN = ATAN(XMIN)
         DELX = XMAX-XMIN
         X   = XMIN + R(1)*DELX
         Q(4) = RM2+RMG*TAN(X)

         IF (Q(4).LE.0) THEN
            W = 0
            RETURN
         ENDIF
         QV = SQRT(Q(4))
         W = W*( (Q(4)-RM2)**2+RMG**2 )*DELX/(TPI*RMG)
      endif !NW
       
C  ... this factor to the weight W corresponds to dq**2/2pi
      ENDIF

C  Next get the center of mass rapidity YCM
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
      YCM = YMIN + R(6)*DELY
      SINHY = SINH(YCM)
      COSHY = SQRT(1+SINHY**2)
      W = W*DELY

C
C  now get Q and SHAT etc.
      SUMPST = PNJ(1)**2+PNJ(2)**2+(PNJ(3)*COSHY-PNJ(0)*SINHY)**2
      Q0ST = SQRT(Q(4)+SUMPST)
      RSHAT = Q0ST + SQRT(PNJ(4)+SUMPST)
      PCM(0) = RSHAT*COSHY
      PCM(3) = RSHAT*SINHY

      DO 40 I = 0,3
        Q(I) = PCM(I)-PNJ(I)
40    CONTINUE

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

      if (bos.eq.2) then
        QV = Q(0)**2-Q(1)**2-Q(2)**2-Q(3)**2
        IF(ABS(QV/Q(4)-1D0).GT.1D-6) THEN
          PRINT*," "
          PRINT*," MASS OF V WRONG: ",QV, qv/q(4)
          PRINT*," Q = ",Q
          PRINT*," "
        ENDIF
      endif
C
C  insert Jacobian for dtau * delta(PCM(0)-Q(0)-PNJ(0))/(2 Q(0))
      W = W * RSHAT/(S*Q0ST)

      IF (BOS.EQ.2) THEN
C  let V decay into 2 massless leptons
C  generate D1 in V rest frame
         CTH = -1+2*R(2)
         STH = SQRT(1-CTH**2)
         PHID = PHI(1) + TPI*R(3)
         D1(0) = SQRT(Q(4))/2
         D1(1) = D1(0)*STH*COS(PHID)
         D1(2) = D1(0)*STH*SIN(PHID)
         D1(3) = D1(0)*CTH
C
C  boost D1 to its lab frame value and get D2 from momentum conservation
         CALL BOOSTN(D1,Q,D1)      
         DO 50 I = 0,3
            D2(I) = Q(I)-D1(I)
50       CONTINUE
C  phase space factor for the decay
         W = W /(8*PI)                        ! *2*TPI /(32*PI**2)
      ENDIF


c$$$      write(*,*)'INITIAL!'
c$$$      do mu = 0, 3
c$$$         write(*,*)'mu, Higgs =', mu, q(mu)
c$$$         write(*,*)'mu, inc 1 =', mu, k1(mu)
c$$$         write(*,*)'mu, inc 2 =', mu, k2(mu)
c$$$         write(*,*)'mu, out 1 =', mu, p(mu,1)
c$$$         write(*,*)'mu, out 2 =', mu, p(mu,2)
c$$$         write(*,*)'mu, out 3 =', mu, p(mu,3)
c$$$         write(*,*)'mu, sum =', mu, k1(mu) + k2(mu) - q(mu) - p(mu,1) -
c$$$     -        p(mu,2) - p(mu,3)
c$$$      end do         
c$$$      write(*,*)'ini sh =', q(0)*q(0) - q(1)*q(1) - q(2)*q(2) - q(3)*q(3)   

** For ew corns, need to rescale Higgs momentum OS if Higgs is BW distributed
      if (nw .eq. 7) then
         sh = q(0)*q(0) 
         sinc = (k1(0)+k2(0))*(k1(0)+k2(0))
         shin = (k1(0) + k2(0))*q(0)
         do mu = 1, 3
            sh = sh - q(mu)*q(mu)
            sinc = sinc - (k1(mu)+k2(mu))*(k1(mu)+k2(mu))
            shin = shin - (k1(mu) + k2(mu))*q(mu)
         end do
         disc = shin**2 - sinc*sh + rm2*(sinc - 2d0*shin + sh)
         scafac = (sinc-rm2)/(sinc-shin+sqrt(disc))
         do mu = 0, 3
            p(mu,1) = p(mu,1)*scafac
            p(mu,2) = p(mu,2)*scafac
            p(mu,3) = p(mu,3)*scafac
            q(mu) = k1(mu) + k2(mu) - p(mu,1) - p(mu,2) - p(mu,3) 
         end do
         if (scafac .lt. 0d0) then
            w = 0d0
            return
         end if
      end if
** end of rescaling

c$$$      if (scafac .lt. 0d0) then
c$$$         write(*,*)'FINAL - OS!'
c$$$         do mu = 0, 3
c$$$            write(*,*)'mu, Higgs =', mu, q(mu)
c$$$            write(*,*)'mu, inc 1 =', mu, k1(mu)
c$$$            write(*,*)'mu, inc 2 =', mu, k2(mu)
c$$$            write(*,*)'mu, out 1 =', mu, p(mu,1)
c$$$            write(*,*)'mu, out 2 =', mu, p(mu,2)
c$$$            write(*,*)'mu, out 3 =', mu, p(mu,3)
c$$$            write(*,*)'mu, sum =', mu, k1(mu) + k2(mu) - q(mu) - p(mu,1) -
c$$$     -           p(mu,2) - p(mu,3)
c$$$         end do        
c$$$         
c$$$         write(*,*)'ofshell sh =', sh
c$$$         write(*,*)'os sh =', q(0)*q(0) - q(1)*q(1) - q(2)*q(2) - q(3)*q(3)
c$$$      end if

C
C  finally put in flux factor for the parton collisions and conversion to fb
      W = W / (2*RSHAT**2) * 3.89379304D+11

      END


C
C**********************  LIPSN   ****************************************
C
      SUBROUTINE LIPSNthreebod(N,R,RN,
     1                 K1,K2,Q,D1,D2,D3,P,
     2                 X1,X2,W)
      IMPLICIT NONE
#include "VBFNLO/utilities/global.inc"
c
      INTEGER N,NRANDLG,NRANDL,N1,mu
      REAL*8 TPI
      PARAMETER (N1=max_p-nx, NRANDLG=N1*3, nrandl=(N1+1)*3, TPI=2D0*PI )
      REAL*8 R(NRANDL),RG,RN,
     1       K1(0:3),K2(0:3),Q(0:4),D1(0:3),D2(0:3),D3(0:3),P(0:3,N1),
     2       X1,X2,W
C
C Phase space generator for production of boson + N jets at hadron colliders
C
C  INPUT:   RD,RG     random number array (between 0 and 1)
C           RN        additional random number which will not be driven by
C                     VEGAS
C
C  OUTPUT:  K1,K2     4-vectors for 2 incoming partons
C           X1,X2     corresponding values of Feynman x
C           Q(0:4)    4-vector for the vector boson resonance V, the 4th
C                     component contains the mass**2
C           D1,D2,D3  4-vectors for the 3 decay leptons of V, assumed massless
C           P(*,N)    4-vectors for the N jets, assumed massless
C           W         the weight of the phase space point
C
C  For a final state photon the ENTRY LIPSV has to be used, which generates a
C  massless four-vector Q
C
C  INPUT INFORMATION on the phase space region is transferred via the 
C    common block /BLIPSIvnj/
C

      REAL*8 RM2,RMG,RM2MIN,RM2MAX, S,M2MIN, YJMIN,YJMAX,PTJMIN,EJMIN
      INTEGER INFOJ
C           S         center of mass energy**2 of the machine
      COMMON /BLIPSIvnj/ RM2,RMG,RM2MIN,RM2MAX,S,M2MIN,
     1    YJMIN(N1),YJMAX(N1),PTJMIN(N1),EJMIN(N1),INFOJ(N1)
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
C  OUTPUT INFORMATION is transferred via the common block /BLIPSO/
C
      REAL*8 Y,PT,PHI,PNJ
      COMMON /BLIPSO/ Y(N1), PT(N1), PHI(N1), PNJ(0:4)
C
C  local variables
      REAL*8 RS,FAC,DELY,XMIN,XMAX,X,DELX,COSHY,SINHY,M2,dotrr,
     1       YCM,YMIN,YMAX, SUMPST,Q0ST,RSHAT,PCM(0:3),
     2       XMJ,QV,YBAR,PJT2,SHA,ACM,WTEMP

      INTEGER I,J,II,JJ,BOS


      external dotrr

      BOS = 2
      GOTO 5

!       ENTRY LIPSV(RG,RN,
!      1            K1,K2,Q,P,
!      2            X1,X2,W)
      DO 2 I = 1,N*3
         R(I) = RG(I)
 2           CONTINUE
      BOS = 1
 5        CONTINUE
      
      FAC = 1/(16*PI**3)
      RS = SQRT(S)
      W = TPI
C  Reset local variables and incoming momenta
      DO 10 I = 0,3
         PCM(I) = 0
         PNJ(I) = 0
         K1(I) = 0
         K2(I) = 0
 10          CONTINUE
      PNJ(4) = 0

C
C  Get 4-momenta of the N massless jets
      DO 20 I = 1,N
         W = W*FAC
         J = 3*(I-1)
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
               XMAX = 1/max(PTJMIN(I),10d0)
               DELX = XMAX-XMIN
               X = XMIN + (R(J+2)-0.2)/0.8d0*DELX
               PT(I) = 1/X
               W = W/0.8d0*DELX*PT(I)**3
            else
               xmin = 0
               xmax = max(ptjmin(i),10d0)
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
            PHI(I) = PHI(I-1) + TPI*(R(J+3)+0.25d0)
            if (phi(i).gt.pi) phi(i) = phi(i)-tpi
            if (phi(i).gt.pi) phi(i) = phi(i)-tpi
         ENDIF
         
         W = W*TPI

         P(1,I) = PT(I)*COS(PHI(I))
         P(2,I) = PT(I)*SIN(PHI(I))
         P(3,I) = PT(I)*SINHY
         DO 30 J = 0,3
            PNJ(J) = PNJ(J)+P(J,I)
 30                CONTINUE
 20                    CONTINUE
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

      IF (BOS.EQ.1) THEN
         Q(4) = 0
         QV = 0
      ELSE
C  Generate Q**2 = Q(4) according to Breit Wigner distribution
         XMAX = (RM2MAX-RM2)/RMG
         XMAX = ATAN(XMAX) 
         XMIN = (RM2MIN-RM2)/RMG
         XMIN = ATAN(XMIN)
         DELX = XMAX-XMIN
         X   = XMIN + R(3*N+1)*DELX
         Q(4) = RM2+RMG*TAN(X)

         IF (Q(4).LE.0) THEN
            W = 0
            RETURN
         ENDIF
         QV = SQRT(Q(4))
         W = W*( (Q(4)-RM2)**2+RMG**2 )*DELX/(TPI*RMG)
       
C  ... this factor to the weight W corresponds to dq**2/2pi
      ENDIF


C  Next get the center of mass rapidity YCM
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

      YCM = YMIN + R(3)*DELY
      SINHY = SINH(YCM)
      COSHY = SQRT(1+SINHY**2)
      W = W*DELY

C
C  now get Q and SHAT etc.
      SUMPST = PNJ(1)**2+PNJ(2)**2+(PNJ(3)*COSHY-PNJ(0)*SINHY)**2
      Q0ST = SQRT(Q(4)+SUMPST)
      RSHAT = Q0ST + SQRT(PNJ(4)+SUMPST)
      PCM(0) = RSHAT*COSHY
      PCM(3) = RSHAT*SINHY

      DO 40 I = 0,3
        Q(I) = PCM(I)-PNJ(I)
 40         CONTINUE

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

      if (bos.eq.2) then
        QV = Q(0)**2-Q(1)**2-Q(2)**2-Q(3)**2
        IF(ABS(QV/Q(4)-1D0).GT.1D-6) THEN
          PRINT*,' '
          PRINT*,' MASS OF V WRONG: ',QV, qv/q(4)
          PRINT*,' Q = ',Q
          PRINT*,' '
        ENDIF
      endif
C
C  insert Jacobian for dtau * delta(PCM(0)-Q(0)-PNJ(0))/(2 Q(0))
      W = W * RSHAT/(S*Q0ST)
      WTEMP=W

c
c  threebody decay
c
      do mu=0,3
      d1(mu)=0d0
      d2(mu)=0d0
      d3(mu)=0d0
      enddo

      Call ThreeBodyDecay0(R(3*N+2), Q(0),Q(4), d1, d2, d3, wtemp)

      w=wtemp/(2*RSHAT**2)

      w=w*3.89379304D+11

      END
