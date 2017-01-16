c ------------------------------------------------------------

      subroutine AZ_TOWWA(P,N,h,a_wwa,z_wwa, NW)

      IMPLICIT NONE

      DOUBLE PRECISION P(0:3,2:7)
      COMPLEX*16 a_wwa(0:5), z_wwa(0:5)
      INTEGER lsign(2:7,2),ll,h,N
      INTEGER NW ! = 0: all diagrams, = 1: no non resonant diagrams

      ! possible lepton helicities:
      DATA (lsign(ll,1),ll=2,7) / -1, -1, -1, -1, -1, 99/  
      DATA (lsign(ll,2),ll=2,7) / -1, -1, -1, +1, +1, 99/  

      call AZ_TOWWA1(P,N,lsign(2,h),a_wwa(0),z_wwa(0),NW)

      ! helicity combi nr.h is selected

      end
      
c ------------------------------------------------------------

      SUBROUTINE AZ_TOWWA1(P,N,lsign,a_wwa,z_wwa,NW)
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
      PARAMETER (NGRAPHS=  19)
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
         CALL JIOXXX(W(1,3   ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,10  ))    
         CALL JIOXXX(W(1,4   ),W(1,5   ),GWF ,WMASS   ,WWIDTH  ,W(1,14  ))    
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
      CALL JVVXXX(W(1,10  ),W(1,6   ),GWWA ,WMASS   ,WWIDTH  ,W(1,13  ))
      CALL FVIXXX(W(1,4   ),W(1,13  ),GWF ,ZERO     ,ZERO    ,W(1,15  ))
      CALL FVIXXX(W(1,11   ), W(1,6   ), GAL, ZERO, ZERO, W(1,16  ))
      CALL JIOXXX(W(1,4   ),W(1,12   ),GWF ,WMASS   ,WWIDTH  ,W(1,17   ))
      CALL JVVXXX(W(1,6  ),W(1,14   ),GWWA ,WMASS   ,WWIDTH  ,W(1,18  ))
      CALL FVOXXX(W(1,2   ),W(1,14  ),GWF ,ZERO     ,ZERO    ,W(1,19  ))
      CALL FVOXXX(W(1,2   ),W(1,17  ),GWF ,ZERO     ,ZERO    ,W(1,20  ))
      CALL FVOXXX(W(1,19   ), W(1,6   ), GAL, ZERO, ZERO, W(1,21  ))
      CALL FVOXXX(W(1,2   ),W(1,18  ),GWF ,ZERO     ,ZERO    ,W(1,22  ))
        
c calculate ampiltudes for ATOWWA
      do mu = 0,3
          CALL TRIBOSVCARTX(Q(0),ZERO ,ZERO, MU,-1,W(1,1   ))
          CALL IOVXXX(W(1,9   ),W(1,5   ),W(1,1   ),GAL ,AMP(1   ))
          CALL IOVXXX(W(1,11   ),W(1,12   ),W(1,1   ),GAL ,AMP(2   ))
          CALL IOVXXX(W(1,15   ),W(1,5   ),W(1,1   ),GAL ,AMP(3   ))
          CALL IOVXXX(W(1,16   ),W(1,5   ),W(1,1   ),GAL ,AMP(4   ))
          CALL VVVXXX(W(1,14  ),W(1,8   ),W(1,1   ),GWWA ,AMP(5   ))
          CALL VVVXXX(W(1,17  ),W(1,10   ),W(1,1   ),GWWA ,AMP(6   ))
          CALL W3W3XX(W(1,14  ),W(1,6  ),W(1,10   ),W(1,1  ),GWWA, GWWA, AMP(7+2   ))
          CALL VVVXXX(W(1,14  ),W(1,13   ),W(1,1   ),GWWA ,AMP(8-1   ))
          CALL VVVXXX(W(1,18  ),W(1,10   ),W(1,1   ),GWWA ,AMP(9-1   ))
          CALL IOVXXX(W(1,7   ),W(1,19   ),W(1,1   ),GAL ,AMP(10   ))
          CALL IOVXXX(W(1,3   ),W(1,20   ),W(1,1   ),GAL ,AMP(11   ))
          CALL IOVXXX(W(1,3   ),W(1,21   ),W(1,1   ),GAL ,AMP(12   ))
          CALL IOVXXX(W(1,3   ),W(1,22   ),W(1,1   ),GAL ,AMP(13   ))

          a_wwa(mu) = AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)

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
          CALL VVVXXX(W(1,14  ),W(1,8   ),W(1,1   ),GWWZ ,AMP(5   ))
          CALL VVVXXX(W(1,17  ),W(1,10   ),W(1,1   ),GWWZ ,AMP(6   ))
          CALL W3W3XX(W(1,14  ),W(1,1  ),W(1,10   ),W(1,6  ),GWWZ, GWWA, AMP(7+2   ))
          CALL VVVXXX(W(1,14  ),W(1,13   ),W(1,1   ),GWWZ ,AMP(8-1   ))
          CALL VVVXXX(W(1,18  ),W(1,10   ),W(1,1   ),GWWZ ,AMP(9-1   ))
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

          z_wwa(mu) = AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &           +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &           +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &           +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)

      enddo
        
      a_wwa(4) = dcmplx(q(0),q(3))
      a_wwa(5) = dcmplx(q(1),q(2))
      z_wwa(4) = dcmplx(q(0),q(3))
      z_wwa(5) = dcmplx(q(1),q(2))
      
      END

c ------------------------------------------------------------

      subroutine ATOWWA(P,N,h,wwa,NW)

      IMPLICIT NONE

      DOUBLE PRECISION P(0:3,2:7)
      COMPLEX*16 wwa(0:5)
      INTEGER lsign(2:7,2),ll,h,N
      INTEGER NW ! = 0: all diagrams, = 1: no non resonant diagrams

      ! possible lepton helicities:
      DATA (lsign(ll,1),ll=2,7) / -1, -1, -1, -1, -1, -1/  !h=1
      DATA (lsign(ll,2),ll=2,7) / -1, -1, -1, -1, +1, +1/  !h=2

      call ATOWWA1(P,N,lsign(2,h),wwa(0),NW)

      ! helicity combi nr.h is selected

      end

c ------------------------------------------------------------

      SUBROUTINE ATOWWA1(P,N,lsign,wwa,NW)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : a -> ve e+ vm~ mu- a  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=  13)
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  22)
      DOUBLE PRECISION     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      DOUBLE PRECISION P(0:3,2:7), Q(0:3)
      INTEGER N
      COMPLEX*16 wwa(0:5)
      INTEGER lsign(2:7)
      INTEGER NW ! = 0: all diagrams, = 1: no non resonant diagrams
C  
C LOCAL VARIABLES 
C  
      INTEGER I, mu
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
c  lepton spinors and W+-, Z/A polarization vector
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
         CALL JIOXXX(W(1,3   ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,10  ))    
         CALL JIOXXX(W(1,4   ),W(1,5   ),GWF ,WMASS   ,WWIDTH  ,W(1,14  ))    
c         print*,' ve ',(wve(i,N)/w(i,2),i=1,6)
c         print*,' ep ',(wep(i,N)/w(i,3),i=1,6)
c         print*,' mu ',(wmu(i,N)/w(i,4),i=1,6)
c         print*,' vm ',(wvm(i,N)/w(i,5),i=1,6)
c         print*,' A ',( atau(i,N)/w(i,6),i=1,6)
c         print*,' wp ',( wp(i,N)/w(i,10),i=1,6)
c         print*,' wm ',( wm(i,N)/w(i,14),i=1,6)
      endif
      do mu = 0,3
      CALL TRIBOSVCARTX(Q(0),ZERO ,ZERO, MU,-1,W(1,1   )) 
      CALL FVOXXX(W(1,5   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL FVIXXX(W(1,3   ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,8   ))    
      CALL JIOXXX(W(1,8   ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,9   ))    
      CALL IOVXXX(W(1,4   ),W(1,7   ),W(1,9   ),GWF ,AMP(1   ))            
      CALL FVIXXX(W(1,4   ),W(1,10  ),GWF ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL FVIXXX(W(1,11  ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL IOVXXX(W(1,12  ),W(1,5   ),W(1,6   ),GAL ,AMP(2   ))            
      CALL JVVXXX(W(1,10  ),W(1,6   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     13  ))                                                          
      CALL IOVXXX(W(1,4   ),W(1,7   ),W(1,13  ),GWF ,AMP(3   ))            
      CALL IOVXXX(W(1,11  ),W(1,7   ),W(1,6   ),GAL ,AMP(4   ))            
      CALL VVVXXX(W(1,14  ),W(1,9   ),W(1,1   ),GWWA ,AMP(5   ))           
      CALL JVVXXX(W(1,10  ),W(1,1   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     15  ))                                                          
      CALL FVIXXX(W(1,4   ),W(1,15  ),GWF ,ZERO    ,ZERO    ,W(1,16  ))    
      CALL IOVXXX(W(1,16  ),W(1,5   ),W(1,6   ),GAL ,AMP(6   ))            
      CALL VVVXXX(W(1,14  ),W(1,13  ),W(1,1   ),GWWA ,AMP(7   ))           
      CALL VVVXXX(W(1,14  ),W(1,15  ),W(1,6   ),GWWA ,AMP(8   ))           
      CALL JW3WXX(W(1,6   ),W(1,10  ),W(1,1   ),GWWA ,GWWA ,WMASS   ,      
     &     WWIDTH  ,W(1,17  ))                                             
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,17  ),GWF ,AMP(9   ))            
      CALL FVIXXX(W(1,8   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL IOVXXX(W(1,18  ),W(1,2   ),W(1,14  ),GWF ,AMP(10  ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GAL ,ZERO    ,ZERO    ,W(1,19  ))    
      CALL JIOXXX(W(1,19  ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,20  ))    
      CALL FVIXXX(W(1,4   ),W(1,20  ),GWF ,ZERO    ,ZERO    ,W(1,21  ))    
      CALL IOVXXX(W(1,21  ),W(1,5   ),W(1,6   ),GAL ,AMP(11  ))            
      CALL FVIXXX(W(1,19  ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,22  ))    
      CALL IOVXXX(W(1,22  ),W(1,2   ),W(1,14  ),GWF ,AMP(12  ))            
      CALL VVVXXX(W(1,14  ),W(1,20  ),W(1,6   ),GWWA ,AMP(13  ))           
c      wwa(mu) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
c     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
c     &             -AMP(  11)-AMP(  12)-AMP(  13)
      wwa(mu) = AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)
      enddo
      wwa(4) = dcmplx(q(0),q(3))
      wwa(5) = dcmplx(q(1),q(2))
c
      END

c ------------------------------------------------------------

      subroutine ZTOWWA(P,N,h,wwa,NW)

      IMPLICIT NONE

      DOUBLE PRECISION P(0:3,2:7)
      COMPLEX*16 wwa(0:5)
      INTEGER lsign(2:7,2),ll,h,N
      INTEGER NW ! = 0: all diagrams, = 1: no non resonant diagrams

      ! possible lepton helicities:
      DATA (lsign(ll,1),ll=2,7) / -1, -1, -1, -1, -1, -1/  !h=1
      DATA (lsign(ll,2),ll=2,7) / -1, -1, -1, -1, +1, +1/  !h=2

      call ZTOWWA1(P,N,lsign(2,h),wwa(0),NW)

      ! helicity combi nr.h is selected

      end

c ------------------------------------------------------------


      SUBROUTINE ZTOWWA1(P,N,lsign,wwa,NW)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : z -> ve e+ vm~ mu- a  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=  19)
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  28)
      DOUBLE PRECISION     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      DOUBLE PRECISION P(0:3,2:7), Q(0:3)
      INTEGER N
      COMPLEX*16 wwa(0:5)
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
         CALL JIOXXX(W(1,3   ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,10  ))    
         CALL JIOXXX(W(1,4   ),W(1,5   ),GWF ,WMASS   ,WWIDTH  ,W(1,14  ))    
c         print*,' ve ',(wve(i,N)/w(i,2),i=1,6)
c         print*,' ep ',(wep(i,N)/w(i,3),i=1,6)
c         print*,' mu ',(wmu(i,N)/w(i,4),i=1,6)
c         print*,' vm ',(wvm(i,N)/w(i,5),i=1,6)
c         print*,' A ',( atau(i,N)/w(i,6),i=1,6)
c         print*,' wp ',( wp(i,N)/w(i,10),i=1,6)
c         print*,' wm ',( wm(i,N)/w(i,14),i=1,6)
      endif
      do mu = 0,3
      CALL TRIBOSVCARTX(Q(0),ZMASS ,ZWIDTH, MU,-1,W(1,1   )) 
      CALL FVOXXX(W(1,5   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL FVIXXX(W(1,3   ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,8   ))    
      CALL JIOXXX(W(1,8   ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,9   ))    
      CALL IOVXXX(W(1,4   ),W(1,7   ),W(1,9   ),GWF ,AMP(1   ))            
      CALL FVIXXX(W(1,4   ),W(1,10  ),GWF ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL FVIXXX(W(1,11  ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL IOVXXX(W(1,12  ),W(1,5   ),W(1,6   ),GAL ,AMP(2   ))            
      CALL JVVXXX(W(1,10  ),W(1,6   ),GWWA ,WMASS   ,WWIDTH  ,W(1,         
     &     13  ))                                                          
      CALL IOVXXX(W(1,4   ),W(1,7   ),W(1,13  ),GWF ,AMP(3   ))            
      CALL IOVXXX(W(1,11  ),W(1,7   ),W(1,6   ),GAL ,AMP(4   ))            
      CALL VVVXXX(W(1,14  ),W(1,9   ),W(1,1   ),GWWZ ,AMP(5   ))           
      CALL JVVXXX(W(1,10  ),W(1,1   ),GWWZ ,WMASS   ,WWIDTH  ,W(1,         
     &     15  ))                                                          
      CALL FVIXXX(W(1,4   ),W(1,15  ),GWF ,ZERO    ,ZERO    ,W(1,16  ))    
      CALL IOVXXX(W(1,16  ),W(1,5   ),W(1,6   ),GAL ,AMP(6   ))            
      CALL VVVXXX(W(1,14  ),W(1,13  ),W(1,1   ),GWWZ ,AMP(7   ))           
      CALL VVVXXX(W(1,14  ),W(1,15  ),W(1,6   ),GWWA ,AMP(8   ))           
      CALL JW3WXX(W(1,1   ),W(1,10  ),W(1,6   ),GWWZ ,GWWA ,WMASS   ,      
     &     WWIDTH  ,W(1,17  ))                                             
      CALL IOVXXX(W(1,4   ),W(1,5   ),W(1,17  ),GWF ,AMP(9   ))            
      CALL FVIXXX(W(1,4   ),W(1,1   ),GZN ,ZERO    ,ZERO    ,W(1,18  ))    
      CALL JIOXXX(W(1,18  ),W(1,5   ),GWF ,WMASS   ,WWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,8   ),W(1,2   ),W(1,19  ),GWF ,AMP(10  ))            
      CALL FVIXXX(W(1,18  ),W(1,10  ),GWF ,ZERO    ,ZERO    ,W(1,20  ))    
      CALL IOVXXX(W(1,20  ),W(1,5   ),W(1,6   ),GAL ,AMP(11  ))            
      CALL VVVXXX(W(1,19  ),W(1,10  ),W(1,6   ),GWWA ,AMP(12  ))           
      CALL FVIXXX(W(1,8   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,21  ))    
      CALL IOVXXX(W(1,21  ),W(1,2   ),W(1,14  ),GWF ,AMP(13  ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GZL ,ZERO    ,ZERO    ,W(1,22  ))    
      CALL JIOXXX(W(1,22  ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,23  ))    
      CALL FVIXXX(W(1,4   ),W(1,23  ),GWF ,ZERO    ,ZERO    ,W(1,24  ))    
      CALL IOVXXX(W(1,24  ),W(1,5   ),W(1,6   ),GAL ,AMP(14  ))            
      CALL FVIXXX(W(1,22  ),W(1,6   ),GAL ,ZERO    ,ZERO    ,W(1,25  ))    
      CALL IOVXXX(W(1,25  ),W(1,2   ),W(1,14  ),GWF ,AMP(15  ))            
      CALL VVVXXX(W(1,14  ),W(1,23  ),W(1,6   ),GWWA ,AMP(16  ))           
      CALL FVOXXX(W(1,2   ),W(1,1   ),GZN ,ZERO    ,ZERO    ,W(1,26  ))    
      CALL IOVXXX(W(1,8   ),W(1,26  ),W(1,14  ),GWF ,AMP(17  ))            
      CALL JIOXXX(W(1,3   ),W(1,26  ),GWF ,WMASS   ,WWIDTH  ,W(1,27  ))    
      CALL FVIXXX(W(1,4   ),W(1,27  ),GWF ,ZERO    ,ZERO    ,W(1,28  ))    
      CALL IOVXXX(W(1,28  ),W(1,5   ),W(1,6   ),GAL ,AMP(18  ))            
      CALL VVVXXX(W(1,14  ),W(1,27  ),W(1,6   ),GWWA ,AMP(19  ))           
c      wwa(mu) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
c     &             -AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
c     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)-AMP(  15)
c     &             -AMP(  16)-AMP(  17)-AMP(  18)-AMP(  19)
      wwa(mu) = AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)+AMP(  15)
     &             +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)
      enddo
      wwa(4) = dcmplx(q(0),q(3))
      wwa(5) = dcmplx(q(1),q(2))
c
      END
