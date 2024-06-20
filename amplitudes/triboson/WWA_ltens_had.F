c ------------------------------------------------------------

      subroutine AZ_TOWWA_had(P,N,h,a_wwa,z_wwa, NW)

      IMPLICIT NONE

      DOUBLE PRECISION P(0:3,2:7)
      COMPLEX*16 a_wwa(0:5), z_wwa(0:5)
      INTEGER lsign(2:7,2),ll,h,N
      INTEGER NW ! = 0: all diagrams, = 1: no non resonant diagrams

      ! possible lepton helicities:
      DATA (lsign(ll,1),ll=2,7) / -1, -1, -1, -1, -1, 99/  
      DATA (lsign(ll,2),ll=2,7) / -1, -1, -1, +1, +1, 99/  

      call AZ_TOWWA1_had(P,N,lsign(2,h),a_wwa(0),z_wwa(0),NW)

      ! helicity combi nr.h is selected

      end
      
c ------------------------------------------------------------

      SUBROUTINE AZ_TOWWA1_had(P,N,lsign,a_wwa,z_wwa,NW)
C  
C  adapted from ZTOWWA1 and ATOWWA1
C 
C 
C  
C FOR PROCESS : z -> ve e+ vm~ mu- a  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=  26)
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  43)
      DOUBLE PRECISION     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      DOUBLE PRECISION P(0:3,2:7), Q(0:3)
      INTEGER N
      COMPLEX*16 a_wwa(0:5), z_wwa(0:5)
      INTEGER lsign(2:7)
      integer NW ! = 0: all diagrams, = 1: no non resonant diagrams
C  
C LOCAL VARIABLES 
C  
      INTEGER I, mu
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)

c  lepton spinors and W+-, Z/A polarization vectors
#include "tensor.inc"

      logical lswitch
      parameter (lswitch=.true.)

C
C GLOBAL VARIABLES
C
#include "VBFNLO/utilities/coupl.inc"      
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"

C                                                                          
C ----------
C BEGIN CODE
C ----------
      if (lswitch) then      ! read in lepton an W states from common block
         do i=1,6
            w(i,2)  = wve(i,N)
            w(i,3)  = wep(i,N)
            w(i,4)  = wvm(i,N)
            w(i,5)  = wmu(i,N)
            w(i,6)  = atau(i,N)
            w(i,10) = wp(i,N)
            w(i,14) = wm(i,N)
         enddo
      endif

      do mu = 0,3
         q(mu) = p(mu,2)+p(mu,3)+p(mu,4)+p(mu,5)+p(mu,6)+p(mu,7)
      enddo


c calculate spinors needed for ATOWWA and ZTOWWA

      CALL FVIXXX(W(1,3   ),W(1,6   ),GWP_ADL ,ZERO     ,ZERO    ,W(1,7   ))
      CALL JIOXXX(W(1,7   ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,8   ))
      CALL FVIXXX(W(1,4   ),W(1,8  ),GWF ,ZERO     ,ZERO    ,W(1,9  ))
      CALL FVIXXX(W(1,4   ),W(1,10  ),GWF ,ZERO     ,ZERO    ,W(1,11  ))
      CALL FVOXXX(W(1,5  ),W(1,6   ),GWM_ADL ,ZERO     ,ZERO    ,W(1,12  ))
      CALL jwmwpa_anomal3(W(1,10  ),W(1,6   ),W(1,13  ))
      CALL FVIXXX(W(1,4   ),W(1,13  ),GWF ,ZERO     ,ZERO    ,W(1,15  ))
      CALL FVIXXX(W(1,11   ), W(1,6   ), GWM_ADL, ZERO, ZERO, W(1,16  ))
      CALL JIOXXX(W(1,4   ),W(1,12   ),GWF ,WMASS   ,WWIDTH  ,W(1,17   ))
      CALL jwpawm_anomal3(W(1,6  ),W(1,14   ),W(1,18  ))
      CALL FVOXXX(W(1,2   ),W(1,14  ),GWF ,ZERO     ,ZERO    ,W(1,19  ))
      CALL FVOXXX(W(1,2   ),W(1,17  ),GWF ,ZERO     ,ZERO    ,W(1,20  ))
      CALL FVOXXX(W(1,19   ), W(1,6   ), GWP_ADL, ZERO, ZERO, W(1,21  ))
      CALL FVOXXX(W(1,2   ),W(1,18  ),GWF ,ZERO     ,ZERO    ,W(1,22  ))

      CALL FVOXXX(W(1,5   ),W(1,8  ),GWF ,ZERO     ,ZERO     ,W(1,23  ))
      CALL FVOXXX(W(1,12   ),W(1,10  ),GWF ,ZERO     ,ZERO   ,W(1,24  ))
      CALL FVOXXX(W(1,5   ),W(1,13  ),GWF ,ZERO     ,ZERO    ,W(1,25  ))
      CALL FVIXXX(W(1,7   ),W(1,14  ),GWF ,ZERO     ,ZERO    ,W(1,26  ))
      CALL FVIXXX(W(1,3   ),W(1,17  ),GWF ,ZERO     ,ZERO    ,W(1,27  ))
      CALL FVIXXX(W(1,3   ),W(1,18  ),GWF ,ZERO     ,ZERO    ,W(1,28  ))

      CALL jhww_anomal3(W(1,14  ),W(1,10  ),W(1,29  ))  !for new higgs graphs   

      CALL FVOXXX(W(1,2  ),W(1,6   ),GWP_AUN ,ZERO     ,ZERO    ,W(1,30  ))
      CALL JIOXXX(W(1,3  ),W(1,30  ),GWF ,WMASS    ,WWIDTH  ,W(1,31  ))
      CALL FVIXXX(W(1,4  ),W(1,31  ),GWF ,ZERO     ,ZERO    ,W(1,32  ))
      CALL FVOXXX(W(1,30 ),W(1,14  ),GWF ,ZERO     ,ZERO    ,W(1,33  ))
      CALL FVIXXX(W(1,3  ),W(1,14  ),GWF ,ZERO     ,ZERO    ,W(1,34  ))
      CALL FVIXXX(W(1,34 ),W(1,6   ),GWP_AUN ,ZERO     ,ZERO    ,W(1,35  ))

      CALL FVIXXX(W(1,4  ),W(1,6   ),GWM_AUN ,ZERO     ,ZERO    ,W(1,36  ))
      CALL FVIXXX(W(1,36 ),W(1,10  ),GWF ,ZERO     ,ZERO    ,W(1,37  ))
      CALL FVOXXX(W(1,5  ),W(1,10  ),GWF ,ZERO     ,ZERO    ,W(1,38  ))
      CALL FVOXXX(W(1,38 ),W(1,6   ),GWM_AUN ,ZERO     ,ZERO    ,W(1,39  ))
      CALL JIOXXX(W(1,36 ),W(1,5   ),GWF ,WMASS    ,WWIDTH  ,W(1,40  ))
      CALL FVOXXX(W(1,2  ),W(1,40  ),GWF ,ZERO     ,ZERO    ,W(1,41  ))


c calculate ampiltudes for ATOWWA
      do mu = 0,3
          CALL TRIBOSVCARTX(Q(0),ZERO ,ZERO, MU,-1,W(1,1   ))
          CALL IOVXXX(W(1,9   ),W(1,5   ),W(1,1   ),GWM_ADL ,AMP(1   ))
          CALL IOVXXX(W(1,11   ),W(1,12   ),W(1,1   ),GWM_ADL ,AMP(2   ))
          CALL IOVXXX(W(1,15   ),W(1,5   ),W(1,1   ),GWM_ADL ,AMP(3   ))
          CALL IOVXXX(W(1,16   ),W(1,5   ),W(1,1   ),GWM_ADL ,AMP(4   ))
          CALL wwa_anomal3(W(1,14  ),W(1,8   ),W(1,1   ),AMP(5   ))
          CALL wwa_anomal3(W(1,17  ),W(1,10   ),W(1,1   ),AMP(6   ))
          CALL wwaa_anomal4(W(1,14  ),W(1,6  ),W(1,10   ),W(1,1  ), AMP(7+2   ))
          CALL wwa_anomal3(W(1,14  ),W(1,13   ),W(1,1   ),AMP(8-1   ))
          CALL wwa_anomal3(W(1,18  ),W(1,10   ),W(1,1   ),AMP(9-1   ))
          CALL IOVXXX(W(1,7   ),W(1,19   ),W(1,1   ),GWP_ADL ,AMP(10   ))
          CALL IOVXXX(W(1,3   ),W(1,20   ),W(1,1   ),GWP_ADL ,AMP(11   ))
          CALL IOVXXX(W(1,3   ),W(1,21   ),W(1,1   ),GWP_ADL ,AMP(12   ))
          CALL IOVXXX(W(1,3   ),W(1,22   ),W(1,1   ),GWP_ADL ,AMP(13   ))

c     new higgs-coupling: aah
          CALL aah_anomal3(W(1,1   ),W(1,6   ),W(1,29),AMP(14   ))    

          a_wwa(mu) = AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)
     &             +AMP(  14)  ! Higgs contribution

          if (procid.eq.WPhadWMA) then
          ! additional graphs for W+ decay
            CALL IOVXXX(W(1,32  ),W(1,5    ),W(1,1   ),GAL ,AMP(15   ))
            CALL wwa_anomal3(W(1,14  ),W(1,31  ),W(1,1   ),AMP(16  ))
            CALL IOVXXX(W(1,3   ),W(1,33   ),W(1,1   ),GWP_ADL ,AMP(17   ))
            CALL IOVXXX(W(1,34  ),W(1,30   ),W(1,1   ),GWP_AUN ,AMP(18   ))
            CALL IOVXXX(W(1,26  ),W(1,2    ),W(1,1   ),GWP_AUN ,AMP(19   ))
            CALL IOVXXX(W(1,27  ),W(1,2    ),W(1,1   ),GWP_AUN ,AMP(20   ))
            CALL IOVXXX(W(1,35  ),W(1,2    ),W(1,1   ),GWP_AUN ,AMP(21   ))
            CALL IOVXXX(W(1,28  ),W(1,2    ),W(1,1   ),GWP_AUN ,AMP(22   ))
            a_wwa(mu) = a_wwa(mu) + AMP(  15)+AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)
     &                            + AMP(  20)+AMP(  21)+AMP(  22)
          endif

          if (procid.eq.WPWMhadA) then
          ! additional graphs for W- decay
            CALL IOVXXX(W(1,37  ),W(1,5    ),W(1,1   ),GWM_ADL ,AMP(15   ))
            CALL wwa_anomal3(W(1,40  ),W(1,10  ),W(1,1   ),AMP(16  ))
            CALL IOVXXX(W(1,4   ),W(1,23   ),W(1,1   ),GWM_AUN ,AMP(17   ))
            CALL IOVXXX(W(1,36  ),W(1,38   ),W(1,1   ),GWM_AUN ,AMP(18   ))
            CALL IOVXXX(W(1,4   ),W(1,24   ),W(1,1   ),GWM_AUN ,AMP(19   ))
            CALL IOVXXX(W(1,4   ),W(1,39   ),W(1,1   ),GWM_AUN ,AMP(20   ))
            CALL IOVXXX(W(1,4   ),W(1,25   ),W(1,1   ),GWM_AUN ,AMP(21   ))
            CALL IOVXXX(W(1,3   ),W(1,41   ),W(1,1   ),GAL ,AMP(22   ))
            a_wwa(mu) = a_wwa(mu) + AMP(  15)+AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)
     &                            + AMP(  20)+AMP(  21)+AMP(  22)
          endif

      enddo


c calculate extra spinors for ZTOWWA

      CALL FVOXXX(W(1,5  ),W(1,31  ),GWF ,ZERO     ,ZERO    ,W(1,42  ))
      CALL FVIXXX(W(1,3  ),W(1,40  ),GWF ,ZERO     ,ZERO    ,W(1,43  ))

        
c calculate amplitudes for ZTOWWA
      do mu = 0,3
          CALL TRIBOSVCARTX(Q(0),ZMASS ,ZWIDTH, MU,-1,W(1,1   ))
          CALL IOVXXX(W(1,9   ),W(1,5   ),W(1,1   ),GWM_ZDL ,AMP(1   ))
          CALL IOVXXX(W(1,11   ),W(1,12   ),W(1,1   ),GWM_ZDL ,AMP(2   ))
          CALL IOVXXX(W(1,15   ),W(1,5   ),W(1,1   ),GWM_ZDL ,AMP(3   ))
          CALL IOVXXX(W(1,16   ),W(1,5   ),W(1,1   ),GWM_ZDL ,AMP(4   ))
          CALL wwz_anomal3(W(1,14  ),W(1,8   ),W(1,1   ),AMP(5   ))
          CALL wwz_anomal3(W(1,17  ),W(1,10   ),W(1,1   ),AMP(6   ))
          CALL wwza_anomal4(W(1,14  ),W(1,1  ),W(1,10   ),W(1,6  ), AMP(7+2   ))
          CALL wwz_anomal3(W(1,14  ),W(1,13   ),W(1,1   ),AMP(8-1   ))
          CALL wwz_anomal3(W(1,18  ),W(1,10   ),W(1,1   ),AMP(9-1   ))
          CALL IOVXXX(W(1,4   ),W(1,23   ),W(1,1   ),GWM_ZUN ,AMP(10   ))
          CALL IOVXXX(W(1,4   ),W(1,24   ),W(1,1   ),GWM_ZUN ,AMP(11   ))
          CALL IOVXXX(W(1,4   ),W(1,25   ),W(1,1   ),GWM_ZUN ,AMP(12   ))
          CALL IOVXXX(W(1,7   ),W(1,19   ),W(1,1   ),GWP_ZDL ,AMP(13   )) ! 10 -> 13
          CALL IOVXXX(W(1,3   ),W(1,20   ),W(1,1   ),GWP_ZDL ,AMP(14   )) ! 11 -> 14
          CALL IOVXXX(W(1,3   ),W(1,21   ),W(1,1   ),GWP_ZDL ,AMP(15   )) ! 12 -> 15
          CALL IOVXXX(W(1,3   ),W(1,22   ),W(1,1   ),GWP_ZDL ,AMP(16   )) ! 13 -> 16
          CALL IOVXXX(W(1,26   ),W(1,2   ),W(1,1   ),GWP_ZUN ,AMP(17   ))
          CALL IOVXXX(W(1,27   ),W(1,2   ),W(1,1   ),GWP_ZUN ,AMP(18   ))
          CALL IOVXXX(W(1,28   ),W(1,2   ),W(1,1   ),GWP_ZUN ,AMP(19   ))

c     new higgs-coupling: azh
          CALL azh_anomal3(W(1,6   ),W(1,1   ),W(1,29  ),AMP(20   ))    

          z_wwa(mu) = AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &           +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &           +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &           +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)
     &           +AMP(  20)  ! Higgs contribution

          if (procid.eq.WPhadWMA) then
          ! additional graphs for W+ decay
            CALL IOVXXX(W(1,32  ),W(1,5    ),W(1,1   ),GZL ,AMP(21   ))
            CALL wwz_anomal3(W(1,14  ),W(1,31  ),W(1,1   ), AMP(22  ))
            CALL IOVXXX(W(1,4   ),W(1,42   ),W(1,1   ),GZN ,AMP(23   ))
            CALL IOVXXX(W(1,3   ),W(1,33   ),W(1,1   ),GWP_ZDL ,AMP(24   ))
            CALL IOVXXX(W(1,34  ),W(1,30   ),W(1,1   ),GWP_ZUN ,AMP(25   ))
            CALL IOVXXX(W(1,35  ),W(1,2    ),W(1,1   ),GWP_ZUN ,AMP(26   ))
            z_wwa(mu) = z_wwa(mu) + AMP(  21)+AMP(  22)+AMP(  23)+AMP(  24)+AMP(  25)+AMP(  26)
          endif

          if (procid.eq.WPWMhadA) then
          ! additional graphs for W- decay
            CALL IOVXXX(W(1,37  ),W(1,5    ),W(1,1   ),GWM_ZDL ,AMP(21   ))
            CALL wwz_anomal3(W(1,40  ),W(1,10  ),W(1,1   ),AMP(22  ))
            CALL IOVXXX(W(1,36  ),W(1,38   ),W(1,1   ),GWM_ZUN ,AMP(23   ))
            CALL IOVXXX(W(1,4   ),W(1,39   ),W(1,1   ),GWM_ZUN ,AMP(24   ))
            CALL IOVXXX(W(1,3   ),W(1,41   ),W(1,1   ),GZL ,AMP(25   ))
            CALL IOVXXX(W(1,43  ),W(1,2    ),W(1,1   ),GZN ,AMP(26   ))
            z_wwa(mu) = z_wwa(mu) + AMP(  21)+AMP(  22)+AMP(  23)+AMP(  24)+AMP(  25)+AMP(  26)
          endif

      enddo
        
      a_wwa(4) = dcmplx(q(0),q(3))
      a_wwa(5) = dcmplx(q(1),q(2))
      z_wwa(4) = dcmplx(q(0),q(3))
      z_wwa(5) = dcmplx(q(1),q(2))
      
      END

c ------------------------------------------------------------


c ------------------------------------------------------------

      SUBROUTINE WMTOWA_had(v,N,h,wmatens,NW)

      IMPLICIT NONE

      DOUBLE PRECISION P(0:3,2:5), v(0:3,6)
      COMPLEX*16 wmatens(0:5)
      INTEGER lsign(2:5,2),ll,mu,h,N
      INTEGER NW ! = 0: all diagrams, = 1: no non resonant diagrams

      !possible lepton helicities: 
      DATA (lsign(ll,1),ll=2,5) / -1, -1, -1, -1/ !h=1
      DATA (lsign(ll,2),ll=2,5) / -1, -1, +1, +1/ !h=2

      do mu = 0,3
         P(mu,2) = v(mu,3)
         P(mu,3) = v(mu,4)
         P(mu,4) = v(mu,5)
         P(mu,5) = v(mu,6)
      enddo

      call WMTOWA1_had(P,N,lsign(2,h),wmatens,NW)

      ! helicity combi nr.h is selected

      end

c ------------------------------------------------------------

      SUBROUTINE WMTOWA1_had(P,N,lsign,wmatens,NW)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : w- -> vm~ mu- a  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=   3) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=   7) 
      DOUBLE PRECISION     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      DOUBLE PRECISION P(0:3,2:5), Q(0:3)
      INTEGER N
      COMPLEX*16 wmatens(0:5)
      INTEGER lsign(2:5)
      INTEGER NW ! = 0: all diagrams, = 1: no non resonant diagrams

C  
C LOCAL VARIABLES 
C  
      INTEGER I, mu
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)

c  lepton spinors and W+-, Z/A polarization vectors
#include "tensor.inc"

      LOGICAL lswitch
      PARAMETER (lswitch=.true.)

C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"                       
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
C                                                                          
C ----------
C BEGIN CODE
C ----------
      if (lswitch) then      ! read in lepton an W states from common block
         do i=1,6
            w(i,2)  = wvm(i,N)
            w(i,3)  = wmu(i,N)
            w(i,4)  = atau(i,N)
         enddo
      endif

      do mu = 0,3
         q(mu) = p(mu,2)+p(mu,3)+p(mu,4)+p(mu,5)
      enddo

      if (.not.lswitch) then     
         CALL OXXXXX(P(0,2   ),ZERO ,  lsign(2),+1, W(1,2   )) ! mu-       
         CALL IXXXXX(P(0,3   ),ZERO , -lsign(3),-1, W(1,3   )) ! vm~       
         CALL VXXXXX(P(0,4)+P(0,5),ZERO ,  +1,  +1, W(1,4   )) ! A
         print*," mu  ",(wmu(i,N)/w(i,2),i=1,6)
         print*," vm  ",(wvm(i,N)/w(i,3),i=1,6)
         print*," A   ",(atau(i,N)/w(i,4),i=1,6)
      endif
      do mu = 0,3
      CALL TRIBOSVCARTX(Q(0),WMASS ,WWIDTH, MU,-1,W(1,1   ))
      CALL JIOXXX(W(1,2   ),W(1,3   ),GWF ,WMASS   ,WWIDTH  ,W(1,5   ))  
      CALL wwa_anomal3(W(1,5   ),W(1,1   ),W(1,4   ),AMP(1   )) 
      CALL FVIXXX(W(1,2   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,4   ),GWM_ADL ,AMP(2   ))            

      wmatens(mu) = +AMP(   1)+AMP(   2)

          if (procid.eq.WPWMhadA) then
          ! additional graphs for W- decay
            CALL FVIXXX(W(1,2   ),W(1,4   ),GWM_AUN ,ZERO    ,ZERO    ,W(1,7   ))    
            CALL IOVXXX(W(1,7   ),W(1,3   ),W(1,1   ),GWF ,AMP(3   ))            
            wmatens(mu) = wmatens(mu) + AMP(  3)
          endif

      enddo
      wmatens(4) = dcmplx(q(0),q(3))
      wmatens(5) = dcmplx(q(1),q(2))
c
      return
      END


c----------------------------------------------------------------------

      SUBROUTINE WPTOWA_had(v,N,h,wpatens,NW)

      IMPLICIT NONE

      DOUBLE PRECISION v(0:3,6)
      DOUBLE PRECISION P(0:3,2:5)
      
      COMPLEX*16 wpatens(0:5)
      INTEGER lsign(2:5,2),ll,mu,h,N
      INTEGER NW ! = 0: all diagrams, = 1: no non resonant diagrams

      ! possible lepton helicities:
      DATA (lsign(ll,1),ll=2,5) / -1, -1, -1, -1/  !h=1
      DATA (lsign(ll,2),ll=2,5) / -1, -1,  1,  1/  !h=2

      do mu = 0,3
         P(mu,2) = v(mu,1)
         P(mu,3) = v(mu,2)
         P(mu,4) = v(mu,5)
         P(mu,5) = v(mu,6)
      enddo

      call WPTOWA1_had(P,N,lsign(2,h),wpatens,NW)

      ! helicity combi nr.h is selected

      end

c ------------------------------------------------------------

      SUBROUTINE WPTOWA1_had(P,N,lsign,wpatens,NW)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : w+ -> ve e+ a  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=   3)
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=   7)
      DOUBLE PRECISION     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      DOUBLE PRECISION P(0:3,2:5), Q(0:3)
      INTEGER N
      COMPLEX*16 wpatens(0:5)
      INTEGER lsign(2:5)
      integer NW ! = 0: all diagrams, = 1: no non resonant diagrams
C  
C LOCAL VARIABLES 
C  
      INTEGER I, mu
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)

c  lepton spinors and W+-, Z/A polarization vectors
#include "tensor.inc"

  
C GLOBAL VARIABLES
  
#include "VBFNLO/utilities/coupl.inc"       
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"

      logical lswitch
      parameter (lswitch=.true.)
                       
C                                                                          
C ----------
C BEGIN CODE
C ----------
      if (lswitch) then      ! read in lepton an W states from common block
         do i=1,6
            w(i,2)  = wve(i,N)
            w(i,3)  = wep(i,N)
            w(i,4)  = atau(i,N)
         enddo
      endif

      do mu = 0,3
         q(mu) = p(mu,2)+p(mu,3)+p(mu,4)+p(mu,5)
      enddo

      if (.not.lswitch) then     
         CALL OXXXXX(P(0,2   ),ZERO ,  lsign(2),+1,W(1,2   )) ! ve       
         CALL IXXXXX(P(0,3   ),ZERO , -lsign(3),-1,W(1,3   )) ! e+       
         CALL VXXXXX(P(0,4)+P(0,5),ZERO ,  +1,  +1, W(1,4   )) ! A
         print*,' ve  ',(wve(i,N)/w(i,2),i=1,6)
         print*,' e+  ',(wep(i,N)/w(i,3),i=1,6)
         print*,' A  ',(atau(i,N)/w(i,4),i=1,6)
      endif
      do mu = 0,3
      CALL TRIBOSVCARTX(Q(0),WMASS ,WWIDTH, MU,-1,W(1,1   ))
      CALL JIOXXX(W(1,3   ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,5   )) 
      CALL wwa_anomal3(W(1,1   ),W(1,5   ),W(1,4   ),AMP(1   ))  
      CALL FVOXXX(W(1,2   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,3   ),W(1,6   ),W(1,4   ),GWP_ADL ,AMP(2   ))            

      wpatens(mu) = AMP(   1)+AMP(   2)

          if (procid.eq.WPhadWMA) then
          ! additional graphs for W+ decay
            CALL FVOXXX(W(1,2   ),W(1,4   ),GWP_AUN ,ZERO    ,ZERO    ,W(1,7   ))    
            CALL IOVXXX(W(1,3   ),W(1,7   ),W(1,1   ),GWF ,AMP(3   ))            
            wpatens(mu) = wpatens(mu) + AMP(  3)
          endif

      enddo
      wpatens(4) = dcmplx(q(0),q(3))
      wpatens(5) = dcmplx(q(1),q(2))
c
      return
      END


      subroutine WWA_ZtoWW_had(P,N,z,NW)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : z -> ve e+ mu- vm~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=   5) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  13) 
      DOUBLE PRECISION     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      INTEGER N
      DOUBLE PRECISION Q(0:3), P(0:3,2:5)
      complex*16 z(0:5)
      integer NW ! = 0: all diagrams, = 1: no non resonant diagrams
C  
C LOCAL VARIABLES 
C  
      INTEGER I, mu
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
      INTEGER NHEL(2:5), ihel
      DATA (NHEL(IHEL),IHEL=2,5) / -1,+1,-1,+1/
      Save NHEL

c  lepton spinors and W+- polarization vectors

#include "tensor.inc"

      logical lswitch
      parameter (lswitch=.true.)

C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
C                                                                          
C ----------
C BEGIN CODE
C ----------    
  
      if (lswitch) then      ! read in lepton an W states from common block
         do i =1,6
            w(i,2) = wve(i,N)
            w(i,3) = wep(i,N)
            w(i,4) = wmu(i,N)
            w(i,5) = wvm(i,N)
            w(i,6) = wp(i,N)
         enddo
      endif

      do mu = 0,3
         q(mu) = p(mu,2)+p(mu,3)+p(mu,4)+p(mu,5)
      enddo

!      if (.not.lswitch) then
!         CALL OXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1,W(1,2   ))        
!         CALL IXXXXX(P(0,3   ),ZERO ,NHEL(3   ),-1,W(1,3   ))        
!         CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1,W(1,4   ))        
!         CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1,W(1,5   ))        
!         CALL JIOXXX(W(1,3   ),W(1,2),GWF ,WMASS,WWIDTH,W(1,6))    
!      endif
      do mu = 0,3
      CALL TRIBOSVCARTX(Q(0     ),ZMASS ,Zwidth,MU,-1,W(1,1   )) 
      CALL FVOXXX(W(1,4   ),W(1,6   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,1   ),GWM_ZUN ,AMP(1   ))            ! W-
      CALL jwmwpz_anomal3(W(1,6   ),W(1,1   ),W(1,8   ))
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,8   ),GWF ,AMP(2   ))            
      CALL FVOXXX(W(1,4   ),W(1,1   ),GWM_ZDL ,ZERO    ,ZERO    ,W(1,9   ))    ! W-
      CALL IOVXXX(W(1,5   ),W(1,9   ),W(1,6   ),GWF ,AMP(3   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GWP_ZDL ,ZERO    ,ZERO    ,W(1,10  ))    ! W+
      CALL JIOXXX(W(1,10  ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,11  ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,11  ),GWF ,AMP(4   ))            
      CALL FVOXXX(W(1,2   ),W(1,1   ),GWP_ZUN ,ZERO    ,ZERO    ,W(1,12  ))    ! W+
      CALL JIOXXX(W(1,3   ),W(1,12  ),GWF ,WMASS   ,WWIDTH  ,W(1,13  ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,13  ),GWF ,AMP(5   ))    

         z(mu) = +AMP(   2)
     &           +AMP(   1)+AMP(   3)+AMP(   4)+AMP(   5) !non resonant

      enddo
      z(4) = dcmplx(q(0),q(3))
      z(5) = dcmplx(q(1),q(2))
c
      return
      END

c---------------------------------------------------------------------

      subroutine WWA_AtoWW_had(P,N,a, NW)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : a -> ve e+ mu- vm~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=   4) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  12) 
      DOUBLE PRECISION     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      INTEGER N
      DOUBLE PRECISION Q(0:3), P(0:3,2:5)
      complex*16 a(0:5)
      INTEGER NW ! = 0: all diagrams, = 1: no non resonant diagrams
C  
C LOCAL VARIABLES 
C
      INTEGER I, mu
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
      INTEGER NHEL(2:5), ihel
      DATA (NHEL(IHEL),IHEL=2,5) /-1,+1,-1,+1/
      Save NHEL
c  lepton spinors and W+- polarization vectors
#include "tensor.inc"

      logical lswitch
      parameter (lswitch=.true.)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
C                                                                          
C ----------
C BEGIN CODE
C ----------
      if (lswitch) then      ! read in lepton an W states from common block
         do i=1,6
            w(i,2) = wve(i,N)
            w(i,3) = wep(i,N)
            w(i,4) = wmu(i,N)
            w(i,5) = wvm(i,N)
            w(i,6) = wp(i,N)
         enddo
      endif

      do mu = 0,3
         q(mu) = p(mu,2)+p(mu,3)+p(mu,4)+p(mu,5)
      enddo

!      if (.not.lswitch) then
!         CALL OXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1,W(1,2   ))        
!         CALL IXXXXX(P(0,3   ),ZERO ,NHEL(3   ),-1,W(1,3   ))        
!         CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1,W(1,4   ))        
!         CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1,W(1,5   ))        
!         CALL JIOXXX(W(1,3   ),W(1,2),GWF ,WMASS,WWIDTH,W(1,6))    
!      endif
      do mu = 0,3
      CALL TRIBOSVCARTX(Q(0),ZERO ,ZERO, MU,-1,W(1,1   ))
      CALL jwmwpa_anomal3(W(1,6   ),W(1,1   ),W(1,7   ))
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,7   ),GWF ,AMP(1   ))            
      CALL FVOXXX(W(1,4   ),W(1,1   ),GWM_ADL ,ZERO    ,ZERO    ,W(1,8   ))    ! W-
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,6   ),GWF ,AMP(2   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GWP_ADL ,ZERO    ,ZERO    ,W(1,9   ))    ! W+
      CALL JIOXXX(W(1,9   ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,10  ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,10  ),GWF ,AMP(3   ))            

         a(mu) = +AMP(   1)
     &        +AMP(   2)+AMP(   3) !non resonant

      if (procid .eq. WPhadWMA) then
        ! new graphs for W+ decay
        CALL FVOXXX(W(1,2   ),W(1,1   ),GAU ,ZERO    ,ZERO    ,W(1,11  ))    
        CALL JIOXXX(W(1,3   ),W(1,11  ),GWF ,WMASS   ,WWIDTH  ,W(1,12  ))    
        CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,12  ),GWF ,AMP(4   ))            
        a(mu) = a(mu) + AMP(   4)
      endif

      if (procid .eq. WPWMhadA) then
        ! new graphs for W- decay
        CALL JIOXXX(W(1,3   ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,11  ))    
        CALL FVOXXX(W(1,4   ),W(1,11  ),GWF ,ZERO    ,ZERO    ,W(1,12  ))    
        CALL IOVXXX(W(1,5   ),W(1,12  ),W(1,1   ),GAU ,AMP(4   ))            
        a(mu) = a(mu) + AMP(   4)
      endif

      enddo
      a(4) = dcmplx(q(0),q(3))
      a(5) = dcmplx(q(1),q(2))
c
      END

