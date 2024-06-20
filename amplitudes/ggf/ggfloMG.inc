      real*8 MGpbar(0:3,6)
      complex*16 GGFLOamp1(-1:1,-1:1,2,24),MGamp1(-1:1,-1:1,-1:1,24),
     &           GGFLOamp2(2,2,2,90),MGamp2(-1:1,-1:1,-1:1,90),
     &           GGFLOamp3(400),MGamp3(400),
     &		 mgJAMP(24),GGFLOcolor(24)
      real*8 MGampsq1,GGFLOampSq1,MGampsq2,GGFLOampSq2,
     &       MGampsq3,GGFLOampSq3,GGFLOampSq1sumhel

      COMMON /MadGcheck/ GGFLOamp1,MGamp1,GGFLOampSq1,MGampSq1,
     &                 GGFLOamp2,MGamp2,GGFLOampSq2,MGampSq2,
     &                 GGFLOamp3,MGamp3,GGFLOampSq3,MGampSq3,
     &		       MGpbar,mgJAMP,GGFLOcolor

      COMMON /MadGcheck1/GGFLOampSq1sumhel

       logical MGdebug,lgauge2
        common/MGdebugC/MGdebug,lgauge2
       logical MGdebug1,MGPen,MgHex,PRINT_AMP
       common/MGdebugC1/MGdebug1,MGPen,MGHex,PRINT_AMP
c      PARAMETER (MGdebug = .true.) 
c  Replace polarization -> momenta
c       PARAMETER (lgauge2 = .false.)
c higgs effective theory couplings

      double complex      gh(2) , ga(2) , gphph(2), gh4(2),ga4(2)
      double precision    ghh(2), gha(2), gaa(2), dumm
      common /COUPL_HEFT/ gh    , ga    , gphph ,
     &                    ghh   , gha   , gaa   , dumm,
     &			  gh4   , ga4
 

