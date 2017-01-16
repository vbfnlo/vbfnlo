c ------------------------------------------------------------

      subroutine AZ_TOWWA_anomal(P,N,h,a_wwa,z_wwa, NW)

      IMPLICIT NONE

      DOUBLE PRECISION P(0:3,2:7)
      COMPLEX*16 a_wwa(0:5), z_wwa(0:5)
      INTEGER lsign(2:7,2),ll,h,N
      INTEGER NW ! = 0: all diagrams, = 1: no non resonant diagrams

      ! possible lepton helicities:
      DATA (lsign(ll,1),ll=2,7) / -1, -1, -1, -1, -1, 99/  
      DATA (lsign(ll,2),ll=2,7) / -1, -1, -1, +1, +1, 99/  

      call AZ_TOWWA1_anomal(P,N,lsign(2,h),a_wwa(0),z_wwa(0),NW)

      ! helicity combi nr.h is selected

      end
      
c ------------------------------------------------------------

      SUBROUTINE AZ_TOWWA1_anomal(P,N,lsign,a_wwa,z_wwa,NW)
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
      PARAMETER (NGRAPHS=  20)
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  29)
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
      
c
      if (.not.lswitch) then
         CALL OXXXXX(P(0,2   ),ZERO , lsign(2),+1,W(1,2   )) !ve
         CALL IXXXXX(P(0,3   ),ZERO ,-lsign(3),-1,W(1,3   )) !e+
         CALL OXXXXX(P(0,4   ),ZERO , lsign(4),+1,W(1,4   )) !mu-
         CALL IXXXXX(P(0,5   ),ZERO ,-lsign(5),-1,W(1,5   )) !vm~
         CALL VXXXXX(P(0,6)+P(0,7),ZERO , +1  ,+1,W(1,6   )) !A
         CALL JIOXXX(W(1,3   ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,10  ))    !W+
         CALL JIOXXX(W(1,4   ),W(1,5   ),GWF ,WMASS   ,WWIDTH  ,W(1,14  ))    !W-
c         print*,' ve ',(wve(i,N)/w(i,2),i=1,6)
c         print*,' ep ',(wep(i,N)/w(i,3),i=1,6)
c         print*,' mu ',(wmu(i,N)/w(i,4),i=1,6)
c         print*,' vm ',(wvm(i,N)/w(i,5),i=1,6)
c         print*,' A ',( atau(i,N)/w(i,6),i=1,6)
c         print*,' wp ',( wp(i,N)/w(i,10),i=1,6)
c         print*,' wm ',( wm(i,N)/w(i,14),i=1,6)
      endif
      
      
c calculate spinors needed for ATOWWA and ZTOWWA

      CALL FVIXXX(W(1,3   ),W(1,6   ),GAL ,ZERO     ,ZERO    ,W(1,7   ))
      CALL JIOXXX(W(1,7   ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,8   ))
      CALL FVIXXX(W(1,4   ),W(1,8  ),GWF ,ZERO     ,ZERO    ,W(1,9  ))
      CALL FVIXXX(W(1,4   ),W(1,10  ),GWF ,ZERO     ,ZERO    ,W(1,11  ))
      CALL FVOXXX(W(1,5  ),W(1,6   ),GAL ,ZERO     ,ZERO    ,W(1,12  ))

c      CALL JVVXXX(W(1,10  ),W(1,6   ),GWWA ,WMASS   ,WWIDTH  ,W(1,13  ))
      CALL jwmwpa_anomal3(W(1,10  ),W(1,6   ),W(1,13  ))

      CALL FVIXXX(W(1,4   ),W(1,13  ),GWF ,ZERO     ,ZERO    ,W(1,15  ))
      CALL FVIXXX(W(1,11   ), W(1,6   ), GAL, ZERO, ZERO, W(1,16  ))
      CALL JIOXXX(W(1,4   ),W(1,12   ),GWF ,WMASS   ,WWIDTH  ,W(1,17   ))

c      CALL JVVXXX(W(1,6  ),W(1,14   ),GWWA ,WMASS   ,WWIDTH  ,W(1,18  ))
      CALL jwpawm_anomal3(W(1,6  ),W(1,14   ),W(1,18  ))

      CALL FVOXXX(W(1,2   ),W(1,14  ),GWF ,ZERO     ,ZERO    ,W(1,19  ))
      CALL FVOXXX(W(1,2   ),W(1,17  ),GWF ,ZERO     ,ZERO    ,W(1,20  ))
      CALL FVOXXX(W(1,19   ), W(1,6   ), GAL, ZERO, ZERO, W(1,21  ))
      CALL FVOXXX(W(1,2   ),W(1,18  ),GWF ,ZERO     ,ZERO    ,W(1,22  ))

      CALL jhww_anomal3(W(1,14  ),W(1,10  ),W(1,29  ))  !for new higgs graphs   
        
c calculate ampiltudes for ATOWWA
      do mu = 0,3
          CALL TRIBOSVCARTX(Q(0),ZERO ,ZERO, MU,-1,W(1,1   ))
          CALL IOVXXX(W(1,9   ),W(1,5   ),W(1,1   ),GAL ,AMP(1   ))
          CALL IOVXXX(W(1,11   ),W(1,12   ),W(1,1   ),GAL ,AMP(2   ))
          CALL IOVXXX(W(1,15   ),W(1,5   ),W(1,1   ),GAL ,AMP(3   ))
          CALL IOVXXX(W(1,16   ),W(1,5   ),W(1,1   ),GAL ,AMP(4   ))

c          CALL VVVXXX(W(1,14  ),W(1,8   ),W(1,1   ),GWWA ,AMP(5   ))
          CALL wwa_anomal3(W(1,14  ),W(1,8   ),W(1,1   ),AMP(5   ))

c          CALL VVVXXX(W(1,17  ),W(1,10   ),W(1,1   ),GWWA ,AMP(6   ))
          CALL wwa_anomal3(W(1,17  ),W(1,10   ),W(1,1   ),AMP(6   ))

c          CALL W3W3XX(W(1,14  ),W(1,6  ),W(1,10   ),W(1,1  ),GWWA, GWWA, AMP(7+2   ))
          CALL wwaa_anomal4(W(1,14  ),W(1,6  ),W(1,10   ),W(1,1  ), AMP(7+2   ))

c          CALL VVVXXX(W(1,14  ),W(1,13   ),W(1,1   ),GWWA ,AMP(8-1   ))
          CALL wwa_anomal3(W(1,14  ),W(1,13   ),W(1,1   ),AMP(8-1   ))

c          CALL VVVXXX(W(1,18  ),W(1,10   ),W(1,1   ),GWWA ,AMP(9-1   ))
          CALL wwa_anomal3(W(1,18  ),W(1,10   ),W(1,1   ),AMP(9-1   ))

          CALL IOVXXX(W(1,7   ),W(1,19   ),W(1,1   ),GAL ,AMP(10   ))
          CALL IOVXXX(W(1,3   ),W(1,20   ),W(1,1   ),GAL ,AMP(11   ))
          CALL IOVXXX(W(1,3   ),W(1,21   ),W(1,1   ),GAL ,AMP(12   ))
          CALL IOVXXX(W(1,3   ),W(1,22   ),W(1,1   ),GAL ,AMP(13   ))

c     new higgs-coupling: aah
          CALL aah_anomal3(W(1,1   ),W(1,6   ),W(1,29),AMP(20   ))    


          a_wwa(mu) = AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)
     &             +AMP(  20)  ! Higgs contribution

      enddo
        
c calculate additional spinors for ZTOWWA
        
      CALL FVOXXX(W(1,5   ),W(1,8  ),GWF ,ZERO     ,ZERO    ,W(1,23  ))
      CALL FVOXXX(W(1,12   ),W(1,10  ),GWF ,ZERO     ,ZERO    ,W(1,24  ))
      CALL FVOXXX(W(1,5   ),W(1,13  ),GWF ,ZERO     ,ZERO    ,W(1,25  ))
      CALL FVIXXX(W(1,7   ),W(1,14  ),GWF ,ZERO     ,ZERO    ,W(1,26  ))
      CALL FVIXXX(W(1,3   ),W(1,17  ),GWF ,ZERO     ,ZERO    ,W(1,27  ))
      CALL FVIXXX(W(1,3   ),W(1,18  ),GWF ,ZERO     ,ZERO    ,W(1,28  ))
        
c calculate amplitudes for ZTOWWA
      do mu = 0,3
          CALL TRIBOSVCARTX(Q(0),ZMASS ,ZWIDTH, MU,-1,W(1,1   ))
          CALL IOVXXX(W(1,9   ),W(1,5   ),W(1,1   ),GZL ,AMP(1   ))
          CALL IOVXXX(W(1,11   ),W(1,12   ),W(1,1   ),GZL ,AMP(2   ))
          CALL IOVXXX(W(1,15   ),W(1,5   ),W(1,1   ),GZL ,AMP(3   ))
          CALL IOVXXX(W(1,16   ),W(1,5   ),W(1,1   ),GZL ,AMP(4   ))

c          CALL VVVXXX(W(1,14  ),W(1,8   ),W(1,1   ),GWWZ ,AMP(5   ))
          CALL wwz_anomal3(W(1,14  ),W(1,8   ),W(1,1   ),AMP(5   ))

c          CALL VVVXXX(W(1,17  ),W(1,10   ),W(1,1   ),GWWZ ,AMP(6   ))
          CALL wwz_anomal3(W(1,17  ),W(1,10   ),W(1,1   ),AMP(6   ))

c          CALL W3W3XX(W(1,14  ),W(1,1  ),W(1,10   ),W(1,6  ),GWWZ, GWWA, AMP(7+2   ))
          CALL wwza_anomal4(W(1,14  ),W(1,1  ),W(1,10   ),W(1,6  ), AMP(7+2   ))

c          CALL VVVXXX(W(1,14  ),W(1,13   ),W(1,1   ),GWWZ ,AMP(8-1   ))
          CALL wwz_anomal3(W(1,14  ),W(1,13   ),W(1,1   ),AMP(8-1   ))

c          CALL VVVXXX(W(1,18  ),W(1,10   ),W(1,1   ),GWWZ ,AMP(9-1   ))
          CALL wwz_anomal3(W(1,18  ),W(1,10   ),W(1,1   ),AMP(9-1   ))

          CALL IOVXXX(W(1,4   ),W(1,23   ),W(1,1   ),GZN ,AMP(10   ))
          CALL IOVXXX(W(1,4   ),W(1,24   ),W(1,1   ),GZN ,AMP(11   ))
          CALL IOVXXX(W(1,4   ),W(1,25   ),W(1,1   ),GZN ,AMP(12   ))
          CALL IOVXXX(W(1,7   ),W(1,19   ),W(1,1   ),GZL ,AMP(13   )) ! 10 -> 13
          CALL IOVXXX(W(1,3   ),W(1,20   ),W(1,1   ),GZL ,AMP(14   )) ! 11 -> 14
          CALL IOVXXX(W(1,3   ),W(1,21   ),W(1,1   ),GZL ,AMP(15   )) ! 12 -> 15
          CALL IOVXXX(W(1,3   ),W(1,22   ),W(1,1   ),GZL ,AMP(16   )) ! 13 -> 16
          CALL IOVXXX(W(1,26   ),W(1,2   ),W(1,1   ),GZN ,AMP(17   ))
          CALL IOVXXX(W(1,27   ),W(1,2   ),W(1,1   ),GZN ,AMP(18   ))
          CALL IOVXXX(W(1,28   ),W(1,2   ),W(1,1   ),GZN ,AMP(19   ))

c     new higgs-coupling: azh
          CALL azh_anomal3(W(1,6   ),W(1,1   ),W(1,29  ),AMP(20   ))    

          z_wwa(mu) = AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &           +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &           +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &           +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)
     &           +AMP(  20)  ! Higgs contribution

      enddo
        
      a_wwa(4) = dcmplx(q(0),q(3))
      a_wwa(5) = dcmplx(q(1),q(2))
      z_wwa(4) = dcmplx(q(0),q(3))
      z_wwa(5) = dcmplx(q(1),q(2))
      
      END

c ------------------------------------------------------------
