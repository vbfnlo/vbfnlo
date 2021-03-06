c here are 2 subroutines DIBOSwptowpzora DIBOSwmtowmzora

      subroutine DIBOSwptowpz_had(P,N,wplus)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : w+ -> ve e+ mu- mu+   
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=   7)
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  15) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      INTEGER N
      REAL*8 Q(0:3), P(0:3,2:5)
      complex*16 wplus(0:5)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
C  
C LOCAL VARIABLES 
C  
      INTEGER I, mu
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
      INTEGER NHEL(2:5), ihel
      DATA (NHEL(IHEL),IHEL=2,5) / +1,-1,-1,+1/
      Save NHEL
c  lepton spinors and W+- polarization vectors
#include "tensor.inc"

      logical lswitch
      parameter (lswitch=.true.)
C                                                                          
C ----------
C BEGIN CODE
C ----------
      if (lswitch) then      ! read in lepton an W states from common block
         do i=1,6
            w(i,2) = wve(i,N)
            w(i,3) = wep(i,N)
            w(i,4) = wmu(i,N)
            w(i,5) = wmubar(i,N)
         enddo
      endif

      do mu = 0,3
         q(mu) = p(mu,2)+p(mu,3)+p(mu,4)+p(mu,5)
      enddo

      call anomal_formfactor(q,q,q,q) ! only first entry relevant

      if (.not.lswitch) then
         CALL IXXXXX(P(0,3   ),ZERO ,NHEL(2   ),-1,W(1,2   ))
         CALL OXXXXX(P(0,2   ),ZERO ,NHEL(3   ),+1,W(1,3   ))        
         CALL OXXXXX(P(0,5   ),ZERO ,NHEL(4   ),+1,W(1,4   ))        
         CALL IXXXXX(P(0,4   ),ZERO ,NHEL(5   ),-1,W(1,5   ))        
      endif  

      do mu = 0,3
      CALL DIBOSVCARTX(Q(0     ),WMASS ,Wwidth,MU,-1,W(1,1   )) 
      CALL JIOXXX(W(1,3   ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,6   ))    
      CALL FVOXXX(W(1,4   ),W(1,6   ),GWF ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,1   ),GWF ,AMP(1   ))            
      CALL jaww_anomal3(W(1,1   ),W(1,6   ),W(1,8   ))                        
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,8   ),GZ_AF ,AMP(2   ))            ! Z decay
      CALL jzww_anomal3(W(1,1   ),W(1,6   ), W(1, 9   )) 
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,9   ),GZ_ZF ,AMP(3   ))            ! Z decay
      CALL FVIXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL JIOXXX(W(1,10  ),W(1,2   ),GWP_ZUN ,ZMASS   ,ZWIDTH  ,W(1,11  ))  ! W decay
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,11  ),GZ_ZF ,AMP(4   ))            ! Z decay
      CALL FVOXXX(W(1,2   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL JIOXXX(W(1,3   ),W(1,12  ),GWP_ADL ,ZERO    ,AWIDTH  ,W(1,13  ))  ! W decay
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,13  ),GZ_AF ,AMP(5   ))            ! Z decay
      CALL JIOXXX(W(1,3   ),W(1,12  ),GWP_ZDL ,ZMASS   ,ZWIDTH  ,W(1,14  ))  ! W decay
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,14  ),GZ_ZF ,AMP(6   ))            ! Z decay

      wplus(mu) = -AMP( 1)-AMP( 2)-AMP( 3)-AMP( 4)-AMP( 5)-AMP(6)

      if (procid .eq. WPhadZ) then
        ! new graph for W decay
        CALL JIOXXX(W(1,10  ),W(1,2   ),GWP_AUN ,ZERO   ,ZERO  ,W(1,15  ))    ! W decay
        CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,15  ),GAL ,AMP(7   ))            
        wplus(mu) = wplus(mu) -AMP( 7)
      endif

      if (procid.eq.WPZhad .and. Ztouu) then
        ! replace graph 1
      CALL FVIXXX(W(1,5   ),W(1,6   ),GWF ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,4   ),W(1,1   ),GWF ,AMP(7   ))            
        wplus(mu) = wplus(mu) - AMP( 7) + AMP( 1)
      endif

      enddo
      wplus(4) = dcmplx(q(0),q(3))
      wplus(5) = dcmplx(q(1),q(2))
c
      return
      END


      subroutine DIBOSwmtowmz_had(P,N,wminus)
C  
C Generated by MadGraph II Version 3.0. Updated 02/19/04                
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : w- -> e- ve~ mu- mu+   
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=   7) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  15) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      INTEGER N
      REAL*8 Q(0:3), P(0:3,2:5)
      complex*16 wminus(0:5)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"
C  
C LOCAL VARIABLES 
C  
      INTEGER I, mu
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
      INTEGER NHEL(2:5), ihel
      DATA (NHEL(IHEL),IHEL=2,5) / +1,-1,-1,+1/
      Save NHEL
c  lepton spinors and W+- polarization vectors
#include "tensor.inc"

      logical lswitch
      parameter (lswitch=.true.)
C                                                                          
C ----------
C BEGIN CODE
C ----------
      if (lswitch) then      ! read in lepton an W states from common block
         do i=1,6
            w(i,2) = wve(i,N)
            w(i,3) = wep(i,N)
            w(i,4) = wmu(i,N)
            w(i,5) = wmubar(i,N)
         enddo

      endif

      do mu = 0,3
         q(mu) = p(mu,2)+p(mu,3)+p(mu,4)+p(mu,5)
      enddo

      call anomal_formfactor(q,q,q,q) ! only first entry relevant

      if (.not.lswitch) then
         CALL IXXXXX(P(0,3   ),ZERO ,NHEL(2   ),-1,W(1,2   ))
         CALL OXXXXX(P(0,2   ),ZERO ,NHEL(3   ),+1,W(1,3   ))        
         CALL OXXXXX(P(0,5   ),ZERO ,NHEL(4   ),+1,W(1,4   ))        
         CALL IXXXXX(P(0,4   ),ZERO ,NHEL(5   ),-1,W(1,5   ))        
      endif  

      do mu = 0,3
      CALL DIBOSVCARTX(Q(0     ),WMASS ,Wwidth,MU,-1,W(1,1   )) 
      CALL JIOXXX(W(1,3   ),W(1,2   ),GWF ,WMASS   ,WWIDTH  ,W(1,6   ))
      CALL jaww_anomal3(W(1,6   ),W(1,1   ),W(1,7   ))    
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,7   ),GZ_AF ,AMP(1   ))            ! Z decay
      CALL jzww_anomal3(W(1,6   ),W(1,1   ), W(1, 8   )) 
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,8   ),GZ_ZF ,AMP(2   ))            ! Z decay
      CALL FVOXXX(W(1,4   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,5   ),W(1,9   ),W(1,6   ),GWF ,AMP(3   ))            
      CALL FVIXXX(W(1,3   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL JIOXXX(W(1,10  ),W(1,2   ),GWM_ADL ,ZERO    ,AWIDTH  ,W(1,11  ))  ! W decay
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,11  ),GZ_AF ,AMP(4   ))            ! Z decay
      CALL JIOXXX(W(1,10  ),W(1,2   ),GWM_ZDL ,ZMASS   ,ZWIDTH  ,W(1,12  ))  ! W decay
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,12  ),GZ_ZF ,AMP(5   ))            ! Z decay
      CALL FVOXXX(W(1,2   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL JIOXXX(W(1,3   ),W(1,13  ),GWM_ZUN ,ZMASS   ,ZWIDTH  ,W(1,14  ))  ! W decay
      CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,14  ),GZ_ZF ,AMP(6   ))            ! Z decay

      wminus(mu) = -AMP( 1)-AMP( 2)-AMP( 3)-AMP( 4)-AMP( 5)-AMP(6)

      if (procid .eq. WMhadZ) then
        ! new graph for W decay
        CALL JIOXXX(W(1,3   ),W(1,13  ),GWM_AUN ,ZERO   ,ZERO  ,W(1,15  ))    ! W decay
        CALL IOVXXX(W(1,5   ),W(1,4   ),W(1,15  ),GAL ,AMP(7   ))            
        wminus(mu) = wminus(mu) -AMP( 7)
      endif

      if (procid.eq.WMZhad .and. Ztouu) then
        ! replace graph 3
      CALL FVIXXX(W(1,5   ),W(1,1   ),GWF ,ZERO    ,ZERO    ,W(1,15  ))    
      CALL IOVXXX(W(1,15  ),W(1,4   ),W(1,6   ),GWF ,AMP(7   ))            
        wminus(mu) = wminus(mu) - AMP( 7) + AMP( 3)
      endif

      enddo


      wminus(4) = dcmplx(q(0),q(3))
      wminus(5) = dcmplx(q(1),q(2))
c
      return
      END

             
