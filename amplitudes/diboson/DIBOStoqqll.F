C*******************************************************
C
C 	subroutine for processes:
C 	Z -> u~ + d + e+ + v_e
C   &   Z -> u + d~ + e- + v_e~
C
C	returns summed, squared amplitudes 
C
C*******************************************************
C
      subroutine DIBOSZtoQQll(P, N, z)
C input P: momenta of outgoing fermions
C input N: NLO
C input z: return value including the amplitude
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=   5) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS= 13)  
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      INTEGER N
      REAL*8 Q(0:3), P(0:3,2:5)
      complex*16 z(0:5)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J, mu
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS) !wavefunctions 
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

C reading states of q1, q2, l & v_l from common block

      if (lswitch) then      
         do i=1,6
C           w(i,1)            !  Z
            w(i,2) = wep(i,N) !  e+
            w(i,3) = wve(i,N) !  ve
            w(i,4) = wmu(i,N) !  d
            w(i,5) = wvm(i,N) !  u~
            w(i,6) = wp(i,N)  !  W+
            w(i,7) = wm(i,N)  !  W-
         enddo
      endif

      do mu = 0,3
         q(mu) = p(mu,2)+p(mu,3)+p(mu,4)+p(mu,5)
      enddo

c      do mu = 0,3
c      CALL DIBOSVCARTX(q(0), ZMASS, Zwidth, mu, -1, w(1,1)) 

C       Z ---> W+ W-;     (W+ ---> e+ + ve);   (W- ---> d + u~)
c       CALL VVVXXX(w(1,7), w(1,6), w(1,1), gwwz, amp(1))
      CALL JVVXXX(w(1,7), w(1,6), gwwz, Zmass, Zwidth, w(1,9))

C       Z ---> u~ + u1;   u1 ---> d + W+;   (W+ ---> e + ve)
      CALL FVOXXX(w(1,4), w(1,6), gwf, zero, zero, w(1,8))    ! w(1,8) = u1
c      CALL IOVXXX(w(1,5), w(1,8), w(1,1), GWP_ZUN, amp(2))
      CALL JIOXXX(w(1,5), w(1,8), GWP_ZUN, Zmass, Zwidth, w(1,10))

C       Z ---> d + d1~;   d1~ ---> u~ + W+;   (W+ ---> e + ve)
      CALL FVIXXX(w(1,5), w(1,6), gwf, zero, zero, w(1,8))    ! w(1,8) = d1~
c      CALL IOVXXX(w(1,8), w(1,4), w(1,1), GWP_ZDL, amp(3))
      CALL JIOXXX(w(1,8), w(1,4), GWP_ZDL, Zmass, Zwidth, w(1,11))

C       Z ---> e + e1;    e1 ---> ve + W-;  (W- ---> d + u~)
      CALL FVOXXX(w(1,3), w(1,7), gwf, zero, zero, w(1,8))    ! w(1,8)= e1
c      CALL IOVXXX(w(1,2), w(1,8), w(1,1), GWM_ZDL, amp(4))
      CALL JIOXXX(w(1,2), w(1,8), GWM_ZDL, Zmass, Zwidth, w(1,12))

C       Z ---> ve + ve1;  ve1 ---> e + W-;    (W- ---> d + u~)
      CALL FVIXXX(w(1,2), w(1,7), gwf, zero, zero, w(1,8))       ! w(1,8)= ve1
c      CALL IOVXXX(w(1,8), w(1,3), w(1,1), GWM_ZUN, amp(5))
      CALL JIOXXX(w(1,8), w(1,3), GWM_ZUN, Zmass, Zwidth, w(1,13))


    
c      z(mu) = -amp(1)
c      if (narrowwidth.eq.0) then
c        z(mu) = z(mu) - amp(2) - amp(3) - amp(4) - amp(5)
c      endif

c      enddo

      z(:) = -w(:,9)-w(:,10)-w(:,11)-w(:,12)-w(:,13)
      

      z(4) = dcmplx(q(0),q(3))
      z(5) = dcmplx(q(1),q(2))
c
      return
      END



C*******************************************************
C
C 	subroutine for process:
C       A -> u~ + d + e+ + v_e
C   &   A -> u + d~ + e- + v_e~
C
C	returns summed, squared amplitudes 
C   	
C*******************************************************

      subroutine DIBOSAtoQQll(P,N,a)
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=   4) 
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=  13) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      INTEGER N
      REAL*8 Q(0:3), P(0:3,2:5)
      complex*16 a(0:5)
C  
C GLOBAL VARIABLES
C  
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"          
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
C                                                                          
C  
C LOCAL VARIABLES 
C
      INTEGER I,J, mu
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(6,NWAVEFUNCS)
      INTEGER NHEL(2:5), ihel
      DATA (NHEL(IHEL),IHEL=2,5) /+1,-1,-1,+1/
      Save NHEL
c  lepton spinors and W+- polarization vectors
#include "tensor.inc"

      logical lswitch
      parameter (lswitch=.true.)
C ----------
C BEGIN CODE
C ----------

C reading states of q1, q2, l & v_l from common block
      if (lswitch) then      
         do i=1,6
C           w(i,1)            !  A
            w(i,2) = wep(i,N) !  e+
            w(i,3) = wve(i,N) !  ve
            w(i,4) = wmu(i,N) !  d
            w(i,5) = wvm(i,N) !  u~
            w(i,6) = wp(i,N)  !  W+
            w(i,7) = wm(i,N)  !  W-
         enddo
      endif

      do mu = 0,3
         q(mu) = p(mu,2)+p(mu,3)+p(mu,4)+p(mu,5)
      enddo


c      do mu = 0,3
c      CALL DIBOSVCARTX(Q(0), ZERO, ZERO, mu, -1, W(1,1   ))

C       A ---> W+ W-;     (W+ ---> e+ + ve);   (W- ---> d + u~)
c       CALL VVVXXX(w(1,7), w(1,6), w(1,1), gwwa, amp(1))
      CALL JVVXXX(w(1,7), w(1,6), gwwa, zero, zero, w(1,9))

C       A ---> u~ + u1;   u1 ---> d + W+;   (W+ ---> e+ + ve) GWP_AUN
      CALL FVOXXX(w(1,4), w(1,6), gwf, zero, zero, w(1,8))    ! w(1,8) = u1
c      CALL IOVXXX(w(1,5), w(1,8), w(1,1), GWP_AUN, amp(2))
      CALL JIOXXX(w(1,5), w(1,8), GWP_AUN, zero, zero, w(1,10))

C       A ---> d + d1~;   d1~ ---> u~ + W+;   (W+ ---> e+ + ve) GWP_ADL
      CALL FVIXXX(w(1,5), w(1,6), gwf, zero, zero, w(1,8))    ! w(1,8) = d1~
c      CALL IOVXXX(w(1,8), w(1,4), w(1,1), GWP_ADL, amp(3))
      CALL JIOXXX(w(1,8), w(1,4), GWP_ADL, zero, zero, w(1,11))

C       A ---> e+ + e1;    e1 ---> ve + W-;  (W- ---> d + u~) GWM_ADL
      CALL FVOXXX(w(1,3), w(1,7), gwf, zero, zero, w(1,8))    ! w(1,8)= e1
c      CALL IOVXXX(w(1,2), w(1,8), w(1,1), GWM_ADL, amp(4))
      CALL JIOXXX(w(1,2), w(1,8), GWM_ADL, zero, zero, w(1,12))    

C      !!!! Forbidden process: GWM_AUN is zero here.                !!!!
c      !!!! In charged conjugated final states, ve is given as u    !!!!
C       A ---> ve + ve1; ve1 ---> e+ + W-;   (W- ---> d + u~) GWM_ADL
      CALL FVIXXX(w(1,2), w(1,7), gwf, zero, zero, w(1,8))       ! w(1,8)= ve1
c      CALL IOVXXX(w(1,8), w(1,3), w(1,1), GWM_AUN, amp(5))
      CALL JIOXXX(w(1,8), w(1,3), GWM_AUN, zero, zero, w(1,13))

c      a(mu) = -amp(1)
c      if (narrowwidth.eq.0) then
c        a(mu) = a(mu) - amp(2) - amp(3) - amp(4) - amp(5)
c      endif
c      enddo

      a(:) = -w(:,9)-w(:,10)-w(:,11)-w(:,12)-w(:,13)

      a(4) = dcmplx(q(0),q(3))
      a(5) = dcmplx(q(1),q(2))
c
      END

