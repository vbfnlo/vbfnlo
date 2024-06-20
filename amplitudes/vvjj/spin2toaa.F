c  for V1V2 --> spin-2 particle --> gamma gamma
c (process 240: AAjj)

c by Jessica Frank
       
      subroutine wwsp2toaa(q1,q2,P,lsign,wwsp2aa)

C FOR PROCESS : w+ w- (-> spin-2 particle) -> a a
C  
      IMPLICIT NONE
 
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=  4) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  8) 

      REAL*8  Q1(0:3), Q2(0:3), P(0:3,3:4)
      complex*16 wwsp2aa(0:3,0:3)

      INTEGER  mu,nu,lsign(3:4),i
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(18,NWAVEFUNCS)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)

      double precision absmomsp2sq, absmomwpsq, absmomwmsq, formfac

#include "tensorsp2aa.inc"
#include "VBFNLO/utilities/spin2coupl.inc"

C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
! C                                                                          
C ----------
C BEGIN CODE
C ----------

      if (lswitch) then
      	do i = 1,6

         w(i,3) = sp2au(i)
         w(i,4) = sp2ae(i)

      	enddo
      else 
c final photons
      CALL VXXXXX(P(0,3),ZERO ,lsign(3),1,W(1,3))          !A1
      CALL VXXXXX(P(0,4),ZERO ,lsign(4),1,W(1,4))          !A2
      endif
      
      DO MU = 0,3
         DO NU = 0,3

c initial vector bosons      
      CALL VCARTX(Q1(0),WMASS ,wwidth, MU,-1,W(1,1   ))        ! W+
      CALL VCARTX(Q2(0),WMASS ,wwidth, NU,-1,W(1,2   ))        ! W- 

c singlet diagram

      if ((f1.EQ.0.0d0) .and. (f2.EQ.0.0d0) .and. (f3.EQ.0.0d0) .and.
     &    (f4.EQ.0.0d0) .and. (f5.EQ.0.0d0)) then

      AMP(1  )=0.0d0

      else

      CALL UWWXXX_SING(W(1,1  ),W(1,2 ) , W(1,5  ))
      CALL AATXXX_SING(W(1,3  ),W(1,4  ),W(1,5  ) , AMP(1   ))

      endif

c uncharged triplet diagram

      if ((f6.EQ.0.0d0) .and. (f7.EQ.0.0d0)) then

      AMP(2  )=0.0d0
      AMP(3  )=0.0d0
      AMP(4  )=0.0d0

      else

      CALL UWWXXX_TRIPN(W(1,1  ),W(1,2 ) , W(1,6  ))
      CALL AATXXX_TRIPN(W(1,3  ),W(1,4  ),W(1,6  ) , AMP(2   ))

c charged triplet diagrams

!       CALL UWAXXX_TRIPPM(W(1,1  ),W(1,3 ) , W(1,7  ))
!       CALL WATXXX_TRIPPM(W(1,2  ),W(1,4  ),W(1,7  ) , AMP(3   ))
!       CALL UWAXXX_TRIPPM(W(1,1  ),W(1,4 ) , W(1,8  ))
!       CALL WATXXX_TRIPPM(W(1,2  ),W(1,3  ),W(1,8  ) , AMP(4   ))

      endif

c formfactor
      absmomwpsq=abs((dreal(W(5,1)))**2-(dreal(W(6,1)))**2
     &               -(dimag(W(6,1)))**2-(dimag(W(5,1)))**2)

      absmomwmsq=abs((dreal(W(5,2)))**2-(dreal(W(6,2)))**2
     &               -(dimag(W(6,2)))**2-(dimag(W(5,2)))**2)

      if ((f1.EQ.0.0d0) .and. (f2.EQ.0.0d0) .and. (f3.EQ.0.0d0) .and.
     &    (f4.EQ.0.0d0) .and. (f5.EQ.0.0d0)) then

      absmomsp2sq=abs((dreal(W(17,6)))**2-(dreal(W(18,6)))**2
     &                -(dimag(W(18,6)))**2-(dimag(W(17,6)))**2)

      else

      absmomsp2sq=abs((dreal(W(17,5)))**2-(dreal(W(18,5)))**2
     &                -(dimag(W(18,5)))**2-(dimag(W(17,5)))**2)

      endif

      formfac=(lambdaff**2/(absmomwpsq+lambdaff**2)*
     &         lambdaff**2/(absmomwmsq+lambdaff**2)*
     &         lambdaff**2/(absmomsp2sq+lambdaff**2))**nff


!            wwsp2aa(mu,nu) = (AMP(1   )+AMP(2   )+AMP(3   )+AMP(4   ))
!      &                      *formfac

           wwsp2aa(mu,nu) = (AMP(1   )+AMP(2   ))
     &                      *formfac

	 enddo
        enddo

       end
	

      subroutine zzsp2toaa(q1,q2,P,lsign,zzsp2aa)

C FOR PROCESS : z z (-> spin-2 particle) -> a a  
C  
      IMPLICIT NONE
 
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=  6) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  10) 

      REAL*8  Q1(0:3), Q2(0:3), P(0:3,3:4)
      complex*16 zzsp2aa(0:3,0:3)

      INTEGER  mu,nu,lsign(3:4),i
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(18,NWAVEFUNCS)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)

      double precision absmomsp2sq, absmomz1sq, absmomz2sq, formfac

#include "tensorsp2aa.inc"
#include "VBFNLO/utilities/spin2coupl.inc"

C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
! C                                                                          
C ----------
C BEGIN CODE
C ----------
     
      if (lswitch) then
      	do i = 1,6

         w(i,3) = sp2au(i)
         w(i,4) = sp2ae(i)

      	enddo
      else 
c final photons
      CALL VXXXXX(P(0,3),ZERO ,lsign(3),1,W(1,3))          !A1
      CALL VXXXXX(P(0,4),ZERO ,lsign(4),1,W(1,4))          !A2
      endif
       
      DO MU = 0,3
         DO NU = 0,3

c initial vector bosons      
      CALL VCARTX(Q1(0),ZMASS ,zwidth, MU,-1,W(1,1   ))    
      CALL VCARTX(Q2(0),ZMASS ,zwidth, NU,-1,W(1,2   ))    

c singlet diagrams 

      if ((f1.EQ.0.0d0) .and. (f2.EQ.0.0d0) .and. (f3.EQ.0.0d0) .and.
     &    (f4.EQ.0.0d0) .and. (f5.EQ.0.0d0)) then

      AMP(1  )=0.0d0
      AMP(2  )=0.0d0
      AMP(3  )=0.0d0

      else

      CALL UZZXXX_SING(W(1,1  ),W(1,2 ) , W(1,5 ))
      CALL AATXXX_SING(W(1,3  ),W(1,4  ),W(1,5 ) , AMP(1   ))
!       CALL UAZXXX_SING(W(1,3  ),W(1,1  ) , W(1,6 ))
!       CALL AZTXXX_SING(W(1,4  ),W(1,2  ),W(1,6 ) , AMP(2   ))
!       CALL UAZXXX_SING(W(1,4  ),W(1,1  ) , W(1,7 ))
!       CALL AZTXXX_SING(W(1,3  ),W(1,2  ),W(1,7 ) , AMP(3   ))

      endif

c triplet diagrams (only uncharged possible)

      if ((f6.EQ.0.0d0) .and. (f7.EQ.0.0d0)) then

      AMP(4  )=0.0d0
      AMP(5  )=0.0d0
      AMP(6  )=0.0d0

      else

      CALL UZZXXX_TRIPN(W(1,1  ),W(1,2 ) , W(1,8 ))
      CALL AATXXX_TRIPN(W(1,3  ),W(1,4  ),W(1,8 ) , AMP(4   ))
!       CALL UAZXXX_TRIPN(W(1,3  ),W(1,1  ) , W(1,9 ))
!       CALL AZTXXX_TRIPN(W(1,4  ),W(1,2  ),W(1,9 ) , AMP(5   ))
!       CALL UAZXXX_TRIPN(W(1,4  ),W(1,1  ) , W(1,10 ))
!       CALL AZTXXX_TRIPN(W(1,3  ),W(1,2  ),W(1,10 ) , AMP(6  ))

      endif

c formfactor
      absmomz1sq=abs((dreal(W(5,1)))**2-(dreal(W(6,1)))**2
     &               -(dimag(W(6,1)))**2-(dimag(W(5,1)))**2)

      absmomz2sq=abs((dreal(W(5,2)))**2-(dreal(W(6,2)))**2
     &               -(dimag(W(6,2)))**2-(dimag(W(5,2)))**2)

      if ((f1.EQ.0.0d0) .and. (f2.EQ.0.0d0) .and. (f3.EQ.0.0d0) .and.
     &    (f4.EQ.0.0d0) .and. (f5.EQ.0.0d0)) then

      absmomsp2sq=abs((dreal(W(17,8)))**2-(dreal(W(18,8)))**2
     &                -(dimag(W(18,8)))**2-(dimag(W(17,8)))**2)

      else

      absmomsp2sq=abs((dreal(W(17,5)))**2-(dreal(W(18,5)))**2
     &                -(dimag(W(18,5)))**2-(dimag(W(17,5)))**2)

      endif

      formfac=(lambdaff**2/(absmomz1sq+lambdaff**2)*
     &         lambdaff**2/(absmomz2sq+lambdaff**2)*
     &         lambdaff**2/(absmomsp2sq+lambdaff**2))**nff


!            zzsp2aa(mu,nu) = (AMP(1   )+AMP(2   )+AMP(3   )
!      &                      +AMP(4   )+AMP(5   )+AMP(6   ))
!      &                      *formfac


           zzsp2aa(mu,nu) = (AMP(1   )+AMP(4   ))
     &                      *formfac


	 enddo
        enddo

       end


      subroutine azsp2toaa(q1,q2,P,lsign,azsp2aa)

C FOR PROCESS : a z (-> spin-2 particle) -> a a  
C  
      IMPLICIT NONE
 
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=  6) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  10) 

      REAL*8  Q1(0:3), Q2(0:3), P(0:3,3:4)
      complex*16 azsp2aa(0:3,0:3)

      INTEGER  mu,nu,lsign(3:4),i
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(18,NWAVEFUNCS)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)

      double precision absmomsp2sq, absmomasq, absmomzsq, formfac

#include "tensorsp2aa.inc"
#include "VBFNLO/utilities/spin2coupl.inc"
c
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
! C                                                                          
C ----------
C BEGIN CODE
C ----------
     
      if (lswitch) then
      	do i = 1,6

         w(i,3) = sp2au(i)
         w(i,4) = sp2ae(i)

      	enddo
      else 
c final photons
      CALL VXXXXX(P(0,3),ZERO ,lsign(3),1,W(1,3))          !A1
      CALL VXXXXX(P(0,4),ZERO ,lsign(4),1,W(1,4))          !A2
      endif
       
      DO MU = 0,3
         DO NU = 0,3

c initial vector bosons      
      CALL VCARTX(Q1(0),zero ,zero, MU,-1,W(1,1   ))    
      CALL VCARTX(Q2(0),ZMASS ,zwidth, NU,-1,W(1,2   ))   

c singlet diagrams  

      if ((f1.EQ.0.0d0) .and. (f2.EQ.0.0d0) .and. (f3.EQ.0.0d0) .and.
     &    (f4.EQ.0.0d0) .and. (f5.EQ.0.0d0)) then

      AMP(1  )=0.0d0
      AMP(2  )=0.0d0
      AMP(3  )=0.0d0

      else

c save time if couplings are equal and therefore T-gamma-Z-vertex is zero
        if ((f1.EQ.f2) .and. (f3.EQ.f4)) then   
         AMP(1   )=0.0d0
         W(17,5 )= W(5,1  )+W(5,2  )
         W(18,5 )= W(6,1  )+W(6,2  )
        else

        CALL UAZXXX_SING(W(1,1  ),W(1,2 ) , W(1,5 ))
        CALL AATXXX_SING(W(1,3  ),W(1,4  ),W(1,5 ) , AMP(1   ))

        endif

!       CALL UAAXXX_SING(W(1,3  ),W(1,1  ) , W(1,6 ))
!       CALL AZTXXX_SING(W(1,4  ),W(1,2  ),W(1,6 ) , AMP(2   ))
!       CALL UAAXXX_SING(W(1,4  ),W(1,1  ) , W(1,7 ))
!       CALL AZTXXX_SING(W(1,3  ),W(1,2  ),W(1,7 ) , AMP(3   ))

      endif

c triplet diagrams (only uncharged possible)

      if ((f6.EQ.0.0d0) .and. (f7.EQ.0.0d0)) then

      AMP(4  )=0.0d0
      AMP(5  )=0.0d0
      AMP(6  )=0.0d0

      else

      CALL UAZXXX_TRIPN(W(1,1  ),W(1,2 ) , W(1,8 ))
      CALL AATXXX_TRIPN(W(1,3  ),W(1,4  ),W(1,8 ) , AMP(4   ))
!       CALL UAAXXX_TRIPN(W(1,3  ),W(1,1  ) , W(1,9 ))
!       CALL AZTXXX_TRIPN(W(1,4  ),W(1,2  ),W(1,9 ) , AMP(5   ))
!       CALL UAAXXX_TRIPN(W(1,4  ),W(1,1  ) , W(1,10 ))
!       CALL AZTXXX_TRIPN(W(1,3  ),W(1,2  ),W(1,10 ) , AMP(6  ))

      endif

c formfactor
      absmomasq=abs((dreal(W(5,1)))**2-(dreal(W(6,1)))**2
     &               -(dimag(W(6,1)))**2-(dimag(W(5,1)))**2)

      absmomzsq=abs((dreal(W(5,2)))**2-(dreal(W(6,2)))**2
     &               -(dimag(W(6,2)))**2-(dimag(W(5,2)))**2)

      if ((f1.EQ.0.0d0) .and. (f2.EQ.0.0d0) .and. (f3.EQ.0.0d0) .and.
     &    (f4.EQ.0.0d0) .and. (f5.EQ.0.0d0)) then

      absmomsp2sq=abs((dreal(W(17,8)))**2-(dreal(W(18,8)))**2
     &                -(dimag(W(18,8)))**2-(dimag(W(17,8)))**2)

      else

      absmomsp2sq=abs((dreal(W(17,5)))**2-(dreal(W(18,5)))**2
     &                -(dimag(W(18,5)))**2-(dimag(W(17,5)))**2)

      endif

      formfac=(lambdaff**2/(absmomasq+lambdaff**2)*
     &         lambdaff**2/(absmomzsq+lambdaff**2)*
     &         lambdaff**2/(absmomsp2sq+lambdaff**2))**nff


!             azsp2aa(mu,nu) = (AMP(1   )+AMP(2   )+AMP(3   )
!      &                       +AMP(4   )+AMP(5   )+AMP(6   ))
!      &                       *formfac

            azsp2aa(mu,nu) = (AMP(1   )+AMP(4   ))
     &                       *formfac


	 enddo
        enddo

       end


      subroutine aasp2toaa(q1,q2,P,lsign,aasp2aa)

C FOR PROCESS : a a (-> spin-2 particle) -> a a  
C  
      IMPLICIT NONE
 
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=  6) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  10) 

      REAL*8  Q1(0:3), Q2(0:3), P(0:3,3:4)
      complex*16 aasp2aa(0:3,0:3)

      INTEGER  mu,nu,lsign(3:4),i
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(18,NWAVEFUNCS)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)

      double precision absmomsp2sq, absmoma1sq, absmoma2sq, formfac

#include "tensorsp2aa.inc"

C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/spin2coupl.inc"
   
      logical lswitch
      parameter (lswitch=.true.)
! C                                                                          
C ----------
C BEGIN CODE
C ----------
     
      if (lswitch) then
      	do i = 1,6

         w(i,3) = sp2au(i)
         w(i,4) = sp2ae(i)
         
      	enddo
      else 
c final photons
      CALL VXXXXX(P(0,3),ZERO ,lsign(3),1,W(1,3))          !A1
      CALL VXXXXX(P(0,4),ZERO ,lsign(4),1,W(1,4))          !A2
      endif
       
      DO MU = 0,3
         DO NU = 0,3

c initial vector bosons      
      CALL VCARTX(Q1(0),zero ,zero, MU,-1,W(1,1   ))    
      CALL VCARTX(Q2(0),zero ,zero, NU,-1,W(1,2   )) 

c singlet diagrams 

      if ((f1.EQ.0.0d0) .and. (f2.EQ.0.0d0) .and. (f3.EQ.0.0d0) .and.
     &    (f4.EQ.0.0d0) .and. (f5.EQ.0.0d0)) then

      AMP(1  )=0.0d0
      AMP(2  )=0.0d0
      AMP(3  )=0.0d0

      else

      CALL UAAXXX_SING(W(1,1  ),W(1,2 ) , W(1,5 ))
      CALL AATXXX_SING(W(1,3  ),W(1,4  ),W(1,5 ) , AMP(1   ))
!       CALL UAAXXX_SING(W(1,3  ),W(1,1  ) , W(1,6 ))
!       CALL AATXXX_SING(W(1,4  ),W(1,2  ),W(1,6 ) , AMP(2   ))
!       CALL UAAXXX_SING(W(1,4  ),W(1,1  ) , W(1,7 ))
!       CALL AATXXX_SING(W(1,3  ),W(1,2  ),W(1,7 ) , AMP(3   ))

      endif

c triplet diagrams (only uncharged possible)

      if ((f6.EQ.0.0d0) .and. (f7.EQ.0.0d0)) then

      AMP(4  )=0.0d0
      AMP(5  )=0.0d0
      AMP(6  )=0.0d0

      else

      CALL UAAXXX_TRIPN(W(1,1  ),W(1,2 ) , W(1,8 ))
      CALL AATXXX_TRIPN(W(1,3  ),W(1,4  ),W(1,8 ) , AMP(4   ))
!       CALL UAAXXX_TRIPN(W(1,3  ),W(1,1  ) , W(1,9 ))
!       CALL AATXXX_TRIPN(W(1,4  ),W(1,2  ),W(1,9 ) , AMP(5   ))
!       CALL UAAXXX_TRIPN(W(1,4  ),W(1,1  ) , W(1,10 ))
!       CALL AATXXX_TRIPN(W(1,3  ),W(1,2  ),W(1,10 ) , AMP(6  ))

      endif

c formfactor
      absmoma1sq=abs((dreal(W(5,1)))**2-(dreal(W(6,1)))**2
     &               -(dimag(W(6,1)))**2-(dimag(W(5,1)))**2)

      absmoma2sq=abs((dreal(W(5,2)))**2-(dreal(W(6,2)))**2
     &               -(dimag(W(6,2)))**2-(dimag(W(5,2)))**2)

      if ((f1.EQ.0.0d0) .and. (f2.EQ.0.0d0) .and. (f3.EQ.0.0d0) .and.
     &    (f4.EQ.0.0d0) .and. (f5.EQ.0.0d0)) then

      absmomsp2sq=abs((dreal(W(17,8)))**2-(dreal(W(18,8)))**2
     &                -(dimag(W(18,8)))**2-(dimag(W(17,8)))**2)

      else

      absmomsp2sq=abs((dreal(W(17,5)))**2-(dreal(W(18,5)))**2
     &                -(dimag(W(18,5)))**2-(dimag(W(17,5)))**2)

      endif

      formfac=(lambdaff**2/(absmoma1sq+lambdaff**2)*
     &         lambdaff**2/(absmoma2sq+lambdaff**2)*
     &         lambdaff**2/(absmomsp2sq+lambdaff**2))**nff


!       aasp2aa(mu,nu) = (AMP(1   )+AMP(2   )+AMP(3   )
!      &                 +AMP(4   )+AMP(5   )+AMP(6   ))
!      &                 *formfac

      aasp2aa(mu,nu) = (AMP(1   )+AMP(4   ))
     &                 *formfac


	 enddo
        enddo

       RETURN
	end
	
				
! c   for Higgs instead of spin-2 particle (result should be nearly the same for process 240 with this
! c   (and in phasespace.F in call def_ps_vnj( ) 6 instead of 39) and for process 101) 
! c   (difference: effective coupling constant)
!        
!       subroutine wwsp2toaa(q1,q2,P,lsign,wwsp2aa)
! 
! C FOR PROCESS : w+ w- (-> higgs) -> a a  
! C  
!       IMPLICIT NONE
!  
!       INTEGER    NGRAPHS
!       PARAMETER (NGRAPHS=  1) 
!       INTEGER    NWAVEFUNCS
!       PARAMETER (NWAVEFUNCS=  5) 
! 
!       REAL*8  Q1(0:3), Q2(0:3), P(0:3,3:4)
!       complex*16 wwsp2aa(0:3,0:3)
! 
!       INTEGER  mu,nu,lsign(3:4),i
!       COMPLEX*16 AMP(NGRAPHS)
!       COMPLEX*16 W(6,NWAVEFUNCS)
!       REAL*8     ZERO
!       PARAMETER (ZERO=0D0)
! 
!       double precision v1, paraw, parat
! 
!       double complex cone, gc(2)
!       parameter( cone = ( 1.0d0, 0.0d0 ) )
! 
! 
! #include "tensorsp2aa.inc"
! 
! C GLOBAL VARIABLES
! C  
! #include "VBFNLO/utilities/coupl.inc"
!    
!       logical lswitch
!       parameter (lswitch=.true.)
! ! C                                                                          
! C ----------
! C BEGIN CODE
! C ----------
! 
! c calculate effective h gamma gamma coupling
!       v1 = 2.0d0*wmass/gw
!       paraw=(hmass/2.0d0/wmass)
!       parat=(hmass/2.0d0/tmass)
!       gc(1) = alpha/(3.1415926535897931d0*v1)
!      & *47.0d0/18.0d0*(1.0d0+(paraw)**2*66.0d0/235.0d0
!      & +(paraw)**4*228.0d0/1645.0d0
!      & +(paraw)**6*696.0d0/8225.0d0+
!      & (paraw)**8*5248.0d0/90475.0d0+(paraw)**10*1280.0d0/29939.0d0
!      & +(paraw)**12*54528.0d0/1646645.0d0
!      & -(parat)**2*56.0d0/705.0d0
!      & -(parat)**4*32.0d0/987.0d0)
!       
!       gc(2) = (0.0d0, 0.0d0)
! 
!       if (lswitch) then
!       	do i = 1,6
! 
!          w(i,3) = sp2au(i)
!          w(i,4) = sp2ae(i)
! 
!       	enddo
!       else 
! c final photons
!       CALL VXXXXX(P(0,3),ZERO ,lsign(3),1,W(1,3))          !A1
!       CALL VXXXXX(P(0,4),ZERO ,lsign(4),1,W(1,4))          !A2
!       endif
!        
!       DO MU = 0,3
!          DO NU = 0,3
! 
! c initial vector bosons      
!       CALL VCARTX(Q1(0),WMASS ,wwidth, MU,-1,W(1,1   ))        ! W+
!       CALL VCARTX(Q2(0),WMASS ,wwidth, NU,-1,W(1,2   ))        ! W-
! 
!       CALL HVVXXX(W(1,1  ),W(1,2   ),GWWH ,HMASS   ,HWIDTH  ,W(1,         
!      &     5  ))                                                          
!       CALL VVSHXX(W(1,3  ),W(1,4   ),W(1,5  ),gc ,AMP(1   ))        
!       
!       wwsp2aa(mu,nu) = AMP(1   )
! 
! 
! 	 enddo
!         enddo
! 
!        end
! 	
! 
!       subroutine zzsp2toaa(q1,q2,P,lsign,zzsp2aa)
! 
! C FOR PROCESS : z z (-> higgs) -> a a  
! C  
!       IMPLICIT NONE
!  
!       INTEGER    NGRAPHS
!       PARAMETER (NGRAPHS=  3) 
!       INTEGER    NWAVEFUNCS
!       PARAMETER (NWAVEFUNCS=  7) 
! 
!       REAL*8  Q1(0:3), Q2(0:3), P(0:3,3:4)
!       complex*16 zzsp2aa(0:3,0:3)
! 
!       INTEGER  mu,nu,lsign(3:4),i
!       COMPLEX*16 AMP(NGRAPHS)
!       COMPLEX*16 W(6,NWAVEFUNCS)
!       REAL*8     ZERO
!       PARAMETER (ZERO=0D0)
!       double precision v1, paraw, parat
! 
!       double complex cone, gc(2)
!       parameter( cone = ( 1.0d0, 0.0d0 ) )
! 
! 
! #include "tensorsp2aa.inc"
! 
! C GLOBAL VARIABLES
! C  
! #include "VBFNLO/utilities/coupl.inc"
!    
!       logical lswitch
!       parameter (lswitch=.true.)
! ! C                                                                          
! C ----------
! C BEGIN CODE
! C ----------
! c calculate effective h gamma gamma coupling
!       v1 = 2.0d0*wmass/gw
!       paraw=(hmass/2.0d0/wmass)
!       parat=(hmass/2.0d0/tmass)
!       gc(1) = alpha/(3.1415926535897931d0*v1)
!      & *47.0d0/18.0d0*(1.0d0+(paraw)**2*66.0d0/235.0d0
!      & +(paraw)**4*228.0d0/1645.0d0
!      & +(paraw)**6*696.0d0/8225.0d0+
!      & (paraw)**8*5248.0d0/90475.0d0+(paraw)**10*1280.0d0/29939.0d0
!      & +(paraw)**12*54528.0d0/1646645.0d0
!      & -(parat)**2*56.0d0/705.0d0
!      & -(parat)**4*32.0d0/987.0d0) 
! 
!       gc(2) = (0.0d0, 0.0d0)
! 
!      
!       if (lswitch) then
!       	do i = 1,6
! 
!          w(i,3) = sp2au(i)
!          w(i,4) = sp2ae(i)
! 
!       	enddo
!       else 
! c final photons
!       CALL VXXXXX(P(0,3),ZERO ,lsign(3),1,W(1,3))          !A1
!       CALL VXXXXX(P(0,4),ZERO ,lsign(4),1,W(1,4))          !A2
!       endif
!        
!       DO MU = 0,3
!          DO NU = 0,3
! 
! c initial vector bosons      
!       CALL VCARTX(Q1(0),ZMASS ,zwidth, MU,-1,W(1,1   ))    
!       CALL VCARTX(Q2(0),ZMASS ,zwidth, NU,-1,W(1,2   ))     
! 
!       CALL HVVXXX(W(1,1  ),W(1,2   ),GZZH ,HMASS   ,HWIDTH  ,W(1,         
!      &     5  ))                                                          
!       CALL VVSHXX(W(1,3  ),W(1,4   ),W(1,5  ),gc ,AMP(1   ))    
! c       CALL HVVXXX(W(1,1  ),W(1,3   ),cone ,HMASS   ,HWIDTH  ,W(1,         
! c      &     6  ))                                                          
! c       CALL VVSXXX(W(1,2  ),W(1,4   ),W(1,6  ),cone ,AMP(2   ))   
! c       CALL HVVXXX(W(1,1  ),W(1,4   ),cone ,HMASS   ,HWIDTH  ,W(1,         
! c      &     7  ))                                                          
! c       CALL VVSXXX(W(1,2  ),W(1,3   ),W(1,7  ),cone ,AMP(3   ))    
! 
!       
!       zzsp2aa(mu,nu) = AMP(1   )
! c     &                 +AMP(2   )+AMP(3   )
! 
! 
! 	 enddo
!            enddo
!  
!        end
! 
! 
!       subroutine azsp2toaa(q1,q2,P,lsign,azsp2aa)
! 
! C FOR PROCESS : a z (-> higgs) -> a a  
! C  
!       IMPLICIT NONE
!  
!       INTEGER    NGRAPHS
!       PARAMETER (NGRAPHS=  3) 
!       INTEGER    NWAVEFUNCS
!       PARAMETER (NWAVEFUNCS=  7) 
! 
!       REAL*8  Q1(0:3), Q2(0:3), P(0:3,3:4)
!       complex*16 azsp2aa(0:3,0:3)
! 
!       INTEGER  mu,nu,lsign(3:4),i
!       COMPLEX*16 AMP(NGRAPHS)
!       COMPLEX*16 W(6,NWAVEFUNCS)
!       REAL*8     ZERO
!       PARAMETER (ZERO=0D0)
!       double complex cone
!       parameter( cone = ( 1.0d0, 0.0d0 ) )
! 
! #include "tensorsp2aa.inc"
! c
! C GLOBAL VARIABLES
! C  
! #include "VBFNLO/utilities/coupl.inc"
!    
!       logical lswitch
!       parameter (lswitch=.true.)
! ! C                                                                          
! C ----------
! C BEGIN CODE
! C ----------
!      
!       if (lswitch) then
!       	do i = 1,6
! 
!          w(i,3) = sp2au(i)
!          w(i,4) = sp2ae(i)
! 
!       	enddo
!       else 
! c final photons
!       CALL VXXXXX(P(0,3),ZERO ,lsign(3),1,W(1,3))          !A1
!       CALL VXXXXX(P(0,4),ZERO ,lsign(4),1,W(1,4))          !A2
!       endif
!        
!       DO MU = 0,3
!          DO NU = 0,3
! 
! c initial vector bosons      
!       CALL VCARTX(Q1(0),zero ,zero, MU,-1,W(1,1   ))    
!       CALL VCARTX(Q2(0),ZMASS ,zwidth, NU,-1,W(1,2   ))  
! 
! c       CALL HVVXXX(W(1,1  ),W(1,2   ),cone ,HMASS   ,HWIDTH  ,W(1,         
! c      &     5  ))                                                          
! c       CALL VVSXXX(W(1,3  ),W(1,4   ),W(1,5  ),cone ,AMP(1   ))    
! c       CALL HVVXXX(W(1,1  ),W(1,3   ),cone ,HMASS   ,HWIDTH  ,W(1,         
! c      &     6  ))                                                          
! c       CALL VVSXXX(W(1,2  ),W(1,4   ),W(1,6  ),cone ,AMP(2   ))   
! c       CALL HVVXXX(W(1,1  ),W(1,4   ),cone ,HMASS   ,HWIDTH  ,W(1,         
! c      &     7  ))                                                          
! c       CALL VVSXXX(W(1,2  ),W(1,3   ),W(1,7  ),cone ,AMP(3   ))   
! 
!       
! c      azsp2aa(mu,nu) = AMP(1   )+AMP(2   )+AMP(3   )
!       azsp2aa(mu,nu) = ZERO
! 
! 	 enddo
!         enddo
! 
!        end
! 
! 
!       subroutine aasp2toaa(q1,q2,P,lsign,aasp2aa)
! 
! C FOR PROCESS : a a (-> higgs) -> a a  
! C  
!       IMPLICIT NONE
!  
!       INTEGER    NGRAPHS
!       PARAMETER (NGRAPHS=  3) 
!       INTEGER    NWAVEFUNCS
!       PARAMETER (NWAVEFUNCS=  7) 
! 
!       REAL*8  Q1(0:3), Q2(0:3), P(0:3,3:4)
!       complex*16 aasp2aa(0:3,0:3)
! 
!       INTEGER  mu,nu,lsign(3:4),i
!       COMPLEX*16 AMP(NGRAPHS)
!       COMPLEX*16 W(6,NWAVEFUNCS)
!       REAL*8     ZERO
!       PARAMETER (ZERO=0D0)
!       double complex cone
!       parameter( cone = ( 1.0d0, 0.0d0 ) )
! 
! #include "tensorsp2aa.inc"
! 
! C GLOBAL VARIABLES
! C  
! #include "VBFNLO/utilities/coupl.inc"
!    
!       logical lswitch
!       parameter (lswitch=.true.)
! ! C                                                                          
! C ----------
! C BEGIN CODE
! C ----------
!      
!       if (lswitch) then
!       	do i = 1,6
! 
!          w(i,3) = sp2au(i)
!          w(i,4) = sp2ae(i)
! 
!       	enddo
!       else 
! c final photons
!       CALL VXXXXX(P(0,3),ZERO ,lsign(3),1,W(1,3))          !A1
!       CALL VXXXXX(P(0,4),ZERO ,lsign(4),1,W(1,4))          !A2
!       endif
!        
!       DO MU = 0,3
!          DO NU = 0,3
! 
! c initial vector bosons      
!       CALL VCARTX(Q1(0),zero ,zero, MU,-1,W(1,1   ))    
!       CALL VCARTX(Q2(0),zero ,zero, NU,-1,W(1,2   ))     
! 
! 
! c       CALL HVVXXX(W(1,1  ),W(1,2   ),cone ,HMASS   ,HWIDTH  ,W(1,         
! c      &     5  ))                                                          
! c       CALL VVSXXX(W(1,3  ),W(1,4   ),W(1,5  ),cone ,AMP(1   ))    
! c       CALL HVVXXX(W(1,1  ),W(1,3   ),cone ,HMASS   ,HWIDTH  ,W(1,         
! c      &     6  ))                                                          
! c       CALL VVSXXX(W(1,2  ),W(1,4   ),W(1,6  ),cone ,AMP(2   ))   
! c       CALL HVVXXX(W(1,1  ),W(1,4   ),cone ,HMASS   ,HWIDTH  ,W(1,         
! c      &     7  ))                                                          
! c       CALL VVSXXX(W(1,2  ),W(1,3   ),W(1,7  ),cone ,AMP(3   ))   
! 
!       
! c      aasp2aa(mu,nu) = AMP(1   )+AMP(2   )+AMP(3   )
!        aasp2aa(mu,nu) = ZERO 
! 
! 	 enddo
!         enddo
! 
!        RETURN
! 	end



