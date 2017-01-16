*******************************************************************
*  file for the cached D-functions 
*******************************************************************
*  Content: the cached D0- & Dij-Function and some service routines
*
*     SUBROUTINE D_1M_INIT_CACHE()
*     SUBROUTINE D_1M_CLEAR_CACHE()
*     FUNCTION D_1M_GET(mass, q1sq, q2sq, q3sq, q1q2, q1q3, q2q3)
*     FUNCTION D_1M_VAL(coefficient, key)
*     FUNCTION D0_1M(mass, q1sq, q2sq, q3sq, q1q2, q1q3, q2q3)
*     FUNCTION Dij_1M(coefficient, mass, q1sq, q2sq, q3sq, q1q2, q1q3, q2q3)
*     FUNCTION D_1M_ALL(mass, q1sq, q2sq, q3sq, q1q2, q1q3, q2q3, Dij)
*
*******************************************************************
*  Main dependencies: 
*     amp_aux.f      |  D04, tens4 
*     c_functions.f  |  C_1M_GET, C_1M_VAL 
*******************************************************************
*  Author: Gunnar Klaemke
*  Date: 16.10.07
*  Modified: Michel Kubocz
*  Date: 19.09.07
*  Additional D00ij and D0000 coefficients for Denner-Dittmaier-tensor decomp.,
*  new SUBROUTINE D_1M_ALL2 (used in SUBROUTINE DDtens51m in ggf_amp_aux.F) 
*  with different ordering of tensor coefficients in comparison to D_1M_ALL 
*******************************************************************


**************************************************************************
      SUBROUTINE D_1M_INIT_CACHE()
**************************************************************************
*     Initializes the cache for the D_1M functions.
**************************************************************************
*     Author: Gunnar Klaemke
*     Date: 16.10.07
**************************************************************************

      IMPLICIT NONE
      
#include "d_functions.inc"

      INTEGER i 

      Cache_index=1
      Cache_count=0
      CallCounter=0
      
      do i=1,LT_Size
         LookupTable(i) = 0
      enddo
      
c      print *,""
c      print *," * Initialized cache for D_1M function"
c      print *,""
      RETURN
      END


**************************************************************************
      SUBROUTINE D_1M_CLEAR_CACHE()
**************************************************************************
*     Clears the cache for the D_1M functions.
**************************************************************************
*     Author: Gunnar Klaemke
*     Date: 16.10.7
**************************************************************************

      IMPLICIT NONE
      
#include "d_functions.inc"

      Cache_index=1
      Cache_count=0
      CallCounter=0

      RETURN
      END


**************************************************************************
      FUNCTION D_1M_GET(mass, q1sq, q2sq, q3sq, q1q2, q1q3, q2q3)
**************************************************************************
*     Calculates the Passarino-Veltman D0-function and the Dij tensor 
*     coefficients for 1 mass in the loop. The results are stored in a 
*     cache, unless the results have already been calculated. Also 
*     calculates the results for the permutation q1<->q3 and stores them 
*     in the cache. The function returns the Index at which the arguments
*     and results are stored as a key.
*     Input: 
*     mass            (REAL*8) the mass of the particle in the loop
*     q1sq,q2sq,q3sq  (REAL*8) q1^2, q2^2 and q3^2 of external momenta
*     q1q2,q1q3,q2q3  (REAL*8) dot-products of external momenta
*     Output:
*     D_1m_get       (INTEGER) a key for D_1M_VAL 
**************************************************************************
*     Author: Gunnar Klaemke
*     Date: 16.10.07
*     Modified: Michael Kubocz
*     Date: 30.01.2011      
*     Additional D00ij and D0000 coefficients for Denner-Dittmaier-tensor decomp.
**************************************************************************

      IMPLICIT NONE
      
      INTEGER D_1M_GET

      REAL*8 mass, q1sq, q2sq, q3sq, q1q2, q1q3, q2q3
      REAL*8 f1,f2,f3,msq
      
      REAL*8 arg(7),qsq
      CHARACTER*(8*7) arg_data !org
      INTEGER*4 arg_bits(2*7)
*     use strings for faster comparison of keys with arguments
      EQUIVALENCE (arg_data,arg)  
      EQUIVALENCE (arg_bits,arg)  
*     for numerical computation of the d0-function
      COMPLEX*16 D04, C_1M_VAL
      INTEGER C_1M_GET
      EXTERNAL D04, C_1M_VAL, C_1M_GET
      
#include "d_functions.inc"
#include "c_coefficients.inc"
#include "d_coefficients.inc"

      INTEGER keyq1q2, keyq1q23, keyq12q3, keyq2q3, keyq1q2DD
      REAL*8 q12sq, q23sq, q123sq
      COMPLEX*16 C0q1q2,  C0q1q23,  C0q12q3,  C0q2q3
      COMPLEX*16 Cijq1q2(2,4), Cijq1q23(2,4), Cijq12q3(2,4), Cijq2q3(2,4),
     &           Cijq1q2DD(2,4) 
      real*8 C0q1q2r,C0q1q23r,C0q12q3r,C0q2q3r, 
     &       Cijq1q2r(4,2),Cijq1q23r(4,2),Cijq12q3r(4,2),Cijq2q3r(4,2)
      real*8 C0q1q2i,C0q1q23i,C0q12q3i,C0q2q3i, 
     &       Cijq1q2i(4,2),Cijq1q23i(4,2),Cijq12q3i(4,2),Cijq2q3i(4,2)

      COMPLEX*16 D0, Dij(13,3)
      real*8 D0r,D0i,Dij_Re(13,3),Dij_Im(13,3)

      INTEGER l_ch, ll
*     this is the position in the cache, where the last key was found
      INTEGER last_lookup
      SAVE last_lookup
*     the number of ignored bits defines the number of binary digits that 
*     shall be ignored in comparisons
      INTEGER*4 IGNORED_BITS
      PARAMETER (IGNORED_BITS=2**10-1)   ! 10 bits = 3 digits  !org

*     fill input arguments into an array
      arg(1)=mass
      arg(2)=q1sq
      arg(3)=q2sq
      arg(4)=q3sq
      arg(5)=q1q2
      arg(6)=q1q3
      arg(7)=q2q3
*     overwrite the least significant bits of each argument with 0
*     this allows to disregard the last digits in comparisons,
*     because 2 arguments might differ only in their last digits (for numerical reasons).
*     A naive comparison would consider those as not equal, though they actually are quite equal.
      arg_bits( 1)=AND(arg_bits( 1),NOT(IGNORED_BITS))
      arg_bits( 3)=AND(arg_bits( 3),NOT(IGNORED_BITS))
      arg_bits( 5)=AND(arg_bits( 5),NOT(IGNORED_BITS))
      arg_bits( 7)=AND(arg_bits( 7),NOT(IGNORED_BITS))
      arg_bits( 9)=AND(arg_bits( 9),NOT(IGNORED_BITS))
      arg_bits(11)=AND(arg_bits(11),NOT(IGNORED_BITS))
      arg_bits(13)=AND(arg_bits(13),NOT(IGNORED_BITS))

      CallCounter = CallCounter+1 ! count all calls to this subroutine
      if (Cache_count.ne.0) then   ! only if cache is not empty

*     first, check if the position stored in the lookup table is right
         if (CallCounter.lt.LT_Size) then 
            ll = LookupTable(CallCounter)
         else
            ll = 0
         endif
         
         if ((ll.ne.0).and.(ll.lt.Cache_count)) then ! the values might be at this position
            if (Cache_keys_data(ll).eq.arg_data) then ! arguments are already stored in cache
               D_1M_GET = ll    ! return index
               RETURN
            endif  
         endif
         
*     the lookup table doesn't give the right answer -> the cache has to be searched:

*     begin search at the position where the last key was found ('last_lookup')

         DO l_ch=1,Cache_count       ! loop over all entries
            if (Cache_keys_data(last_lookup).eq.arg_data) then ! arguments are already stored in cache
               LookupTable(CallCounter) = last_lookup ! save position in look-up table
               D_1M_GET = last_lookup ! return index
               RETURN
            endif  
            if (last_lookup.eq.Cache_count) last_lookup=0   ! jump to start, if end is reached
            last_lookup=last_lookup+1
         ENDDO

      endif    

*     the arguments were not found in the cache.
*     => evaluate dij coefficients:

*     calculate further momentum^2
      q12sq =q1sq +q2sq+2*q1q2
      q23sq =q2sq +q3sq+2*q2q3
      q123sq=q12sq+q3sq+2*q2q3+2*q1q3
*     get keys for the C0s and Cijs
      keyq1q2 =C_1M_GET(mass, q1sq, q2sq, q12sq)
      keyq1q23=C_1M_GET(mass, q1sq,q23sq,q123sq)
      keyq12q3=C_1M_GET(mass,q12sq, q3sq,q123sq)
      keyq2q3 =C_1M_GET(mass, q2sq, q3sq, q23sq)
*     get the C0s
      C0q1q2 =C_1M_VAL(cc0, keyq1q2)
      C0q1q23=C_1M_VAL(cc0, keyq1q23)
      C0q12q3=C_1M_VAL(cc0, keyq12q3)
      C0q2q3 =C_1M_VAL(cc0, keyq2q3)

      C0q1q2r =-Dble(C0q1q2)
      C0q1q23r=-Dble(C0q1q23)
      C0q12q3r=-Dble(C0q12q3)
      C0q2q3r =-Dble(C0q2q3)

      C0q1q2i =-DImag(C0q1q2)
      C0q1q23i=-DImag(C0q1q23)
      C0q12q3i=-DImag(C0q12q3)
      C0q2q3i =-DImag(C0q2q3)
*     get the Cijs
*     keyq1q2:
      Cijq1q2(1,1)=C_1M_VAL(cc11, keyq1q2)
      Cijq1q2(1,2)=C_1M_VAL(cc12, keyq1q2)
      Cijq1q2(2,1)=C_1M_VAL(cc21, keyq1q2)
      Cijq1q2(2,2)=C_1M_VAL(cc22, keyq1q2)
      Cijq1q2(2,3)=C_1M_VAL(cc23, keyq1q2)
      Cijq1q2(2,4)=C_1M_VAL(cc24, keyq1q2)

      Cijq1q2r(1,1)=-Dble(Cijq1q2(1,1))
      Cijq1q2r(2,1)=-Dble(Cijq1q2(1,2))
      Cijq1q2r(1,2)=-Dble(Cijq1q2(2,1))
      Cijq1q2r(2,2)=-Dble(Cijq1q2(2,2))
      Cijq1q2r(3,2)=-Dble(Cijq1q2(2,3))
      Cijq1q2r(4,2)= Dble(Cijq1q2(2,4))

      Cijq1q2i(1,1)=-DImag(Cijq1q2(1,1))
      Cijq1q2i(2,1)=-DImag(Cijq1q2(1,2))
      Cijq1q2i(1,2)=-DImag(Cijq1q2(2,1))
      Cijq1q2i(2,2)=-DImag(Cijq1q2(2,2))
      Cijq1q2i(3,2)=-DImag(Cijq1q2(2,3))
      Cijq1q2i(4,2)= DImag(Cijq1q2(2,4))
*     keyq1q23:
      Cijq1q23(1,1)=C_1M_VAL(cc11, keyq1q23)
      Cijq1q23(1,2)=C_1M_VAL(cc12, keyq1q23)
      Cijq1q23(2,1)=C_1M_VAL(cc21, keyq1q23)
      Cijq1q23(2,2)=C_1M_VAL(cc22, keyq1q23)
      Cijq1q23(2,3)=C_1M_VAL(cc23, keyq1q23)
      Cijq1q23(2,4)=C_1M_VAL(cc24, keyq1q23)

      Cijq1q23r(1,1)=-Dble(Cijq1q23(1,1))
      Cijq1q23r(2,1)=-Dble(Cijq1q23(1,2))
      Cijq1q23r(1,2)=-Dble(Cijq1q23(2,1))
      Cijq1q23r(2,2)=-Dble(Cijq1q23(2,2))
      Cijq1q23r(3,2)=-Dble(Cijq1q23(2,3))
      Cijq1q23r(4,2)= Dble(Cijq1q23(2,4))

      Cijq1q23i(1,1)=-DImag(Cijq1q23(1,1))
      Cijq1q23i(2,1)=-DImag(Cijq1q23(1,2))
      Cijq1q23i(1,2)=-DImag(Cijq1q23(2,1))
      Cijq1q23i(2,2)=-DImag(Cijq1q23(2,2))
      Cijq1q23i(3,2)=-DImag(Cijq1q23(2,3))
      Cijq1q23i(4,2)= DImag(Cijq1q23(2,4))
*     keyq12q3:
      Cijq12q3(1,1)=C_1M_VAL(cc11, keyq12q3)
      Cijq12q3(1,2)=C_1M_VAL(cc12, keyq12q3)
      Cijq12q3(2,1)=C_1M_VAL(cc21, keyq12q3)
      Cijq12q3(2,2)=C_1M_VAL(cc22, keyq12q3)
      Cijq12q3(2,3)=C_1M_VAL(cc23, keyq12q3)
      Cijq12q3(2,4)=C_1M_VAL(cc24, keyq12q3)

      Cijq12q3r(1,1)=-Dble(Cijq12q3(1,1))
      Cijq12q3r(2,1)=-Dble(Cijq12q3(1,2))
      Cijq12q3r(1,2)=-Dble(Cijq12q3(2,1))
      Cijq12q3r(2,2)=-Dble(Cijq12q3(2,2))
      Cijq12q3r(3,2)=-Dble(Cijq12q3(2,3))
      Cijq12q3r(4,2)= Dble(Cijq12q3(2,4))

      Cijq12q3i(1,1)=-DImag(Cijq12q3(1,1))
      Cijq12q3i(2,1)=-DImag(Cijq12q3(1,2))
      Cijq12q3i(1,2)=-DImag(Cijq12q3(2,1))
      Cijq12q3i(2,2)=-DImag(Cijq12q3(2,2))
      Cijq12q3i(3,2)=-DImag(Cijq12q3(2,3))
      Cijq12q3i(4,2)= DImag(Cijq12q3(2,4))
*     keyq2q3:
      Cijq2q3(1,1)=C_1M_VAL(cc11, keyq2q3)
      Cijq2q3(1,2)=C_1M_VAL(cc12, keyq2q3)
      Cijq2q3(2,1)=C_1M_VAL(cc21, keyq2q3)
      Cijq2q3(2,2)=C_1M_VAL(cc22, keyq2q3)
      Cijq2q3(2,3)=C_1M_VAL(cc23, keyq2q3)
      Cijq2q3(2,4)=C_1M_VAL(cc24, keyq2q3)

      Cijq2q3r(1,1)=-Dble(Cijq2q3(1,1))
      Cijq2q3r(2,1)=-Dble(Cijq2q3(1,2))
      Cijq2q3r(1,2)=-Dble(Cijq2q3(2,1))
      Cijq2q3r(2,2)=-Dble(Cijq2q3(2,2))
      Cijq2q3r(3,2)=-Dble(Cijq2q3(2,3))
      Cijq2q3r(4,2)= Dble(Cijq2q3(2,4))

      Cijq2q3i(1,1)=-DImag(Cijq2q3(1,1))
      Cijq2q3i(2,1)=-DImag(Cijq2q3(1,2))
      Cijq2q3i(1,2)=-DImag(Cijq2q3(2,1))
      Cijq2q3i(2,2)=-DImag(Cijq2q3(2,2))
      Cijq2q3i(3,2)=-DImag(Cijq2q3(2,3))
      Cijq2q3i(4,2)= DImag(Cijq2q3(2,4))

*     calculate D0
      D0=d04(q1sq,q2sq,q3sq,q123sq,q12sq,q23sq,
     &             mass,mass,mass,mass)
      D0r=Dble(D0)
      D0i=DImag(D0)

*     calculate Dijs
      CALL tens_red4_new_Re_Com_1M(mass,q1sq,q2sq,q3sq,q1q2,q1q3,q2q3,
     &         C0q2q3r,C0q12q3r,C0q1q23r,C0q1q2r,
     &         Cijq2q3r,Cijq12q3r,Cijq1q23r,Cijq1q2r,
     &         C0q2q3i,C0q12q3i,C0q1q23i,C0q1q2i,
     &         Cijq2q3i,Cijq12q3i,Cijq1q23i,Cijq1q2i,
     &         D0,D0r,D0i,Dij_Re,Dij_Im)

c... combine real and imaginery parts
      Dij(1,1)=DCMPLX(Dij_Re(1,1),Dij_Im(1,1))
      Dij(2,1)=DCMPLX(Dij_Re(2,1),Dij_Im(2,1))
      Dij(3,1)=DCMPLX(Dij_Re(3,1),Dij_Im(3,1))
 
      Dij(7,1)=DCMPLX(Dij_Re(7,1),Dij_Im(7,1))
      Dij(8,1)=DCMPLX(Dij_Re(8,1),Dij_Im(8,1))
      Dij(9,1)=DCMPLX(Dij_Re(9,1),Dij_Im(9,1))
      Dij(10,1)=DCMPLX(Dij_Re(10,1),Dij_Im(10,1))
      Dij(11,1)=DCMPLX(Dij_Re(11,1),Dij_Im(11,1))
      Dij(12,1)=DCMPLX(Dij_Re(12,1),Dij_Im(12,1))
      Dij(13,1)=DCMPLX(Dij_Re(13,1),Dij_Im(13,1))
      
      Dij(1,2)=DCMPLX(Dij_Re(1,2),Dij_Im(1,2))
      Dij(2,2)=DCMPLX(Dij_Re(2,2),Dij_Im(2,2))
      Dij(3,2)=DCMPLX(Dij_Re(3,2),Dij_Im(3,2))
      Dij(4,2)=DCMPLX(Dij_Re(4,2),Dij_Im(4,2))
      Dij(5,2)=DCMPLX(Dij_Re(5,2),Dij_Im(5,2))
      Dij(6,2)=DCMPLX(Dij_Re(6,2),Dij_Im(6,2))
      Dij(7,2)=DCMPLX(Dij_Re(7,2),Dij_Im(7,2))
      
      Dij(1,3)=DCMPLX(Dij_Re(1,3),Dij_Im(1,3))
      Dij(2,3)=DCMPLX(Dij_Re(2,3),Dij_Im(2,3))
      Dij(3,3)=DCMPLX(Dij_Re(3,3),Dij_Im(3,3))
      Dij(4,3)=DCMPLX(Dij_Re(4,3),Dij_Im(4,3))
      Dij(5,3)=DCMPLX(Dij_Re(5,3),Dij_Im(5,3))
      Dij(6,3)=DCMPLX(Dij_Re(6,3),Dij_Im(6,3))
      Dij(7,3)=DCMPLX(Dij_Re(7,3),Dij_Im(7,3))
      Dij(8,3)=DCMPLX(Dij_Re(8,3),Dij_Im(8,3))
      Dij(9,3)=DCMPLX(Dij_Re(9,3),Dij_Im(9,3))
      Dij(10,3)=DCMPLX(Dij_Re(10,3),Dij_Im(10,3))
      Dij(11,3)=DCMPLX(Dij_Re(11,3),Dij_Im(11,3))
      Dij(12,3)=DCMPLX(Dij_Re(12,3),Dij_Im(12,3))
      Dij(13,3)=DCMPLX(Dij_Re(13,3),Dij_Im(13,3))

*     fill key & values in the cache
      Cache_keys_data(Cache_index)=arg_data
      Cache_values(dd0  ,Cache_index)=D0

      Cache_values(dd11 ,Cache_index)=Dij(1,1)
      Cache_values(dd12 ,Cache_index)=Dij(2,1)
      Cache_values(dd13 ,Cache_index)=Dij(3,1)

      Cache_values(dd21 ,Cache_index)=Dij(1,2)
      Cache_values(dd22 ,Cache_index)=Dij(2,2)
      Cache_values(dd23 ,Cache_index)=Dij(3,2)
      Cache_values(dd24 ,Cache_index)=Dij(4,2)
      Cache_values(dd25 ,Cache_index)=Dij(5,2)
      Cache_values(dd26 ,Cache_index)=Dij(6,2)
      Cache_values(dd27 ,Cache_index)=-Dij(7,2)  

      Cache_values(dd31 ,Cache_index)=Dij(1,3)
      Cache_values(dd32 ,Cache_index)=Dij(2,3)
      Cache_values(dd33 ,Cache_index)=Dij(3,3)
      Cache_values(dd34 ,Cache_index)=Dij(4,3)
      Cache_values(dd35 ,Cache_index)=Dij(5,3)
      Cache_values(dd36 ,Cache_index)=Dij(6,3)
      Cache_values(dd37 ,Cache_index)=Dij(7,3)
      Cache_values(dd38 ,Cache_index)=Dij(8,3)
      Cache_values(dd39 ,Cache_index)=Dij(9,3)
      Cache_values(dd310,Cache_index)=Dij(10,3)
      Cache_values(dd311,Cache_index)=-Dij(11,3) 
      Cache_values(dd312,Cache_index)=-Dij(12,3) 
      Cache_values(dd313,Cache_index)=-Dij(13,3) 

c Denner-Dittmaier coefficients stored in unused Dij(7,1)-Dij(13,1)
      Cache_values(dd416,Cache_index)=Dij(7,1)  
      Cache_values(dd417,Cache_index)=Dij(8,1)
      Cache_values(dd418,Cache_index)=Dij(9,1)
      Cache_values(dd419,Cache_index)=Dij(10,1)
      Cache_values(dd420,Cache_index)=Dij(11,1)
      Cache_values(dd421,Cache_index)=Dij(12,1)
      Cache_values(dd422,Cache_index)=Dij(13,1)

*     return this index as the functions result
      D_1M_GET = Cache_index
*     reset last look_up
      last_lookup = Cache_index  
*     fill index into lookup table
      LookupTable(CallCounter) = Cache_index

*     increment number of entries in cache
      if (Cache_count.lt.Cache_Size) Cache_count=Cache_count+1      
*     increment current cache position, jump to start if end is reached
      if (Cache_index.eq.Cache_Size) then 
         Cache_index=0
         print *,"  WARNING: Cache overflow in D-functions! "
      endif
      Cache_index=Cache_index+1

c$$$      RETURN  check!!!!

*     repeat this for the permutation q1<->q3
      qsq=arg(4)                ! permute q1sq<->q3sq
      arg(4)=arg(2)
      arg(2)=qsq
      qsq=arg(7)                ! permute q1q2<->q2q3
      arg(7)=arg(5)
      arg(5)=qsq

*     fill key & values in the cache, use symmetry relations
      Cache_keys_data(Cache_index)=arg_data
      Cache_values(dd0  ,Cache_index)= D0
      Cache_values(dd11 ,Cache_index)=-Dij(3,1)-D0
      Cache_values(dd12 ,Cache_index)=-Dij(2,1)-D0
      Cache_values(dd13 ,Cache_index)=-Dij(1,1)-D0
      Cache_values(dd21 ,Cache_index)= Dij(3,2)+2*Dij(3,1)+D0
      Cache_values(dd22 ,Cache_index)= Dij(2,2)+2*Dij(2,1)+D0
      Cache_values(dd23 ,Cache_index)= Dij(1,2)+2*Dij(1,1)+D0
      Cache_values(dd24 ,Cache_index)= Dij(6,2)+Dij(3,1)+Dij(2,1)+D0
      Cache_values(dd25 ,Cache_index)= Dij(5,2)+Dij(3,1)+Dij(1,1)+D0
      Cache_values(dd26 ,Cache_index)= Dij(4,2)+Dij(2,1)+Dij(1,1)+D0
      Cache_values(dd27 ,Cache_index)= -Dij(7,2)
      Cache_values(dd31 ,Cache_index)=-3*Dij(3,1)-3*Dij(3,2)-Dij(3,3)-D0
      Cache_values(dd32 ,Cache_index)=-3*Dij(2,1)-3*Dij(2,2)-Dij(2,3)-D0
      Cache_values(dd33 ,Cache_index)=-3*Dij(1,1)-3*Dij(1,2)-Dij(1,3)-D0
      Cache_values(dd34 ,Cache_index)=-2*Dij(3,1)-Dij(2,1)-2*Dij(6,2)-Dij(9,3)-Dij(3,2)-D0
      Cache_values(dd35 ,Cache_index)=-2*Dij(3,1)-Dij(1,1)-2*Dij(5,2)-Dij(7,3)-Dij(3,2)-D0
      Cache_values(dd36 ,Cache_index)=-Dij(3,1)-2*Dij(2,1)-2*Dij(6,2)-Dij(8,3)-Dij(2,2)-D0
      Cache_values(dd37 ,Cache_index)=-Dij(3,1)-2*Dij(1,1)-2*Dij(5,2)-Dij(5,3)-Dij(1,2)-D0
      Cache_values(dd38 ,Cache_index)=-2*Dij(2,1)-Dij(1,1)-2*Dij(4,2)-Dij(6,3)-Dij(2,2)-D0
      Cache_values(dd39 ,Cache_index)=-Dij(2,1)-2*Dij(1,1)-2*Dij(4,2)-Dij(4,3)-Dij(1,2)-D0
      Cache_values(dd310,Cache_index)=-Dij(3,1)-Dij(2,1)-Dij(1,1)-Dij(6,2)
     &                                -Dij(10,3)-Dij(5,2)-Dij(4,2)-D0
      Cache_values(dd311,Cache_index)=-(-Dij(7,2)-Dij(13,3))
      Cache_values(dd312,Cache_index)=-(-Dij(7,2)-Dij(12,3))
      Cache_values(dd313,Cache_index)=-(-Dij(7,2)-Dij(11,3))

c D00ij and D0000 functions for Denner-Dittmaier E-coefficients !! 
      msq=mass*mass
      f1=arg(2)
      f2=arg(3)+2d0*arg(5)
      f3=arg(4)+2d0*(arg(6)+arg(7))

      keyq1q2DD=C_1M_GET(mass,q2sq,q1sq,q12sq)
      Cijq1q2DD(1,1)=C_1M_VAL(cc11,keyq1q2DD)
      Cijq1q2DD(1,2)=C_1M_VAL(cc12,keyq1q2DD)
      Cijq1q2DD(2,1)=C_1M_VAL(cc21,keyq1q2DD)
      Cijq1q2DD(2,2)=C_1M_VAL(cc22,keyq1q2DD)
      Cijq1q2DD(2,3)=C_1M_VAL(cc23,keyq1q2DD)
      Cijq1q2DD(2,4)=C_1M_VAL(cc24,keyq1q2DD)

      Cache_values(dd416,Cache_index)=
     & (-C0q1q2+f1*Cache_values(dd31,Cache_index)+f2
     & *Cache_values(dd34 ,Cache_index)+f3*Cache_values(dd35,Cache_index)
     & +2d0*msq*Cache_values(dd21 ,Cache_index))/6d0

      Cache_values(dd417,Cache_index)=
     & (-Cijq1q2DD(2,1)+f1*Cache_values(dd36,Cache_index)+f2
     & *Cache_values(dd32,Cache_index)+f3*Cache_values(dd38,Cache_index)
     & +2d0*msq*Cache_values(dd22,Cache_index))/6d0

      Cache_values(dd418,Cache_index)=
     & (-Cijq1q2DD(2,2)+f1*Cache_values(dd37,Cache_index)+f2
     & *Cache_values(dd39,Cache_index)+f3*Cache_values(dd33 ,Cache_index)
     & +2d0*msq*Cache_values(dd23,Cache_index))/6d0

      Cache_values(dd419,Cache_index)=
     & (Cijq1q2DD(1,1)+f1*Cache_values(dd34,Cache_index)+f2
     & *Cache_values(dd36,Cache_index)+f3*Cache_values(dd310,Cache_index)
     & +2d0*msq*Cache_values(dd24,Cache_index))/6d0

      Cache_values(dd420,Cache_index)=
     & (Cijq1q2DD(1,2)+f1*Cache_values(dd35,Cache_index)+f2
     & *Cache_values(dd310,Cache_index)+f3*Cache_values(dd37,Cache_index)
     & +2d0*msq*Cache_values(dd25,Cache_index))/6d0

      Cache_values(dd421,Cache_index)=
     & (-Cijq1q2DD(2,3)+f1*Cache_values(dd310,Cache_index)+f2
     & *Cache_values(dd38,Cache_index)+f3*Cache_values(dd39,Cache_index)
     & +2d0*msq*Cache_values(dd26,Cache_index))/6d0

      Cache_values(dd422,Cache_index)=
     & (Cijq1q2DD(2,4)+1d0/12d0-arg(2)*Cache_values(dd416,Cache_index)
     & -arg(3)*Cache_values(dd417,Cache_index)-arg(4)
     & *Cache_values(dd418,Cache_index)-2d0*arg(5)
     & *Cache_values(dd419,Cache_index)-2d0*arg(6)
     & *Cache_values(dd420,Cache_index)-2d0*arg(7)
     & *Cache_values(dd421,Cache_index)-msq
     & *Cache_values(dd27,Cache_index))/6d0

*     increment number of entries in cache
      if (Cache_count.lt.Cache_Size) Cache_count=Cache_count+1      
*     increment current cache position, jump to start if end is reached
      if (Cache_index.eq.Cache_Size) then 
         Cache_index=0
         print *,"  WARNING: Cache overflow in D-functions! "
      endif
      Cache_index=Cache_index+1
      
*     good bye
      RETURN
      END


**************************************************************************
      FUNCTION D_1M_VAL(coefficient, key)
**************************************************************************
*     Returns the Passarino-Veltman 4-point coefficients for a given key.
*     Use 'D_1M_GET' to obtain a key for some set of momenta.
*     Use the following constants for the coefficient:
*     dd0 : D0 scalar Integral
*     dd11, dd12, dd13,
*     dd21, dd22, dd23, dd24, dd25, dd26, dd27, 
*     dd31, dd32, dd33, dd34, dd35, dd36, dd37, dd38, dd39,  
*     dd310, dd311, dd312, dd313  : Dij tensor coefficients 
*
*     Use this function in the following way:
*     ...
*     mykey =D_1M_GET(my_mass,my_q1sq,my_q2sq,my_q3sq,...)
*     my_D0 =D_1M_VAL(dd0 ,key)
*     my_D12=D_1M_VAL(dd12,key)
*     my_D37=D_1M_VAL(dd37,key)
*     ...
*     Important: Use different keys for different sets of momenta:
*     ...
*     key_A=D_1M_GET(mass_A,q1sq_A,q2sq_A,q3sq_A,...)
*     key_B=D_1M_GET(mass_B,q1sq_B,q2sq_B,q3sq_B,...)
*     ...
**************************************************************************
*     Author: Gunnar Klaemke
*     Date: 22.05.06
**************************************************************************
      
      IMPLICIT NONE
      COMPLEX*16 D_1M_VAL
      INTEGER coefficient, key

#include "d_functions.inc"
#include "d_coefficients.inc"
      
      D_1M_VAL=Cache_values(coefficient,key)

      RETURN
      END


**************************************************************************
      FUNCTION D0_1M(mass, q1sq, q2sq, q3sq, q1q2, q1q3, q2q3)
**************************************************************************
*     Returns the Passarino-Veltman D0 function for one mass in the loop.
**************************************************************************
*     Author: Gunnar Klaemke
*     Date: 22.05.06
**************************************************************************

      IMPLICIT NONE
      COMPLEX*16 D0_1M
      REAL*8 mass, q1sq, q2sq, q3sq, q1q2, q1q3, q2q3
      INTEGER D_1M_GET
      COMPLEX*16 D_1M_VAL
      INTEGER key

#include "d_coefficients.inc"

      key  =D_1M_GET(mass, q1sq, q2sq, q3sq, q1q2, q1q3, q2q3)
      D0_1M=D_1M_VAL(dd0, key)
      
      RETURN
      END


**************************************************************************
      FUNCTION Dij_1M(coefficient, mass, q1sq, q2sq, q3sq, q1q2, q1q3, q2q3)
**************************************************************************
*     Returns the Passarino-Veltman Dij 4-point tensor coefficients for 
*     one mass in the loop.
*     Possible coefficients: dd11, dd12, dd13, dd21, dd22, dd23, dd24, dd25, 
*                            dd26, dd27, dd31, dd32, dd33, dd34, dd35, dd36, 
*                            dd37, dd38, dd39, dd310,dd311,dd312,dd313
**************************************************************************
*     Author: Gunnar Klaemke
*     Date: 22.05.06
**************************************************************************

      IMPLICIT NONE
      COMPLEX*16 Dij_1M
      REAL*8 mass, q1sq, q2sq, q3sq, q1q2, q1q3, q2q3
      INTEGER coefficient
      INTEGER D_1M_GET
      COMPLEX*16 D_1M_VAL
      INTEGER key

#include "d_coefficients.inc"

      key   =D_1M_GET(mass, q1sq, q2sq, q3sq, q1q2, q1q3, q2q3)
      Dij_1M=D_1M_VAL(coefficient, key)
      
      RETURN
      END


**************************************************************************
      SUBROUTINE D_1M_ALL(mass, q1sq, q2sq, q3sq, q1q2, q1q3, q2q3, DDij)
**************************************************************************
*     Returns all the Passarino-Veltman Dij 4-point tensor coefficients for 
*     one mass in the loop.
*
*     Input: 
*     mass            (REAL*8) the mass of the particle in the loop
*     q1sq,q2sq,q3sq  (REAL*8) q1^2, q2^2 and q3^2 of external momenta
*     q1q2,q1q3,q2q3  (REAL*8) dot-products of external momenta
*     Output:
*     DDij(1:24)  (COMPLEX*15) D0=DDij(dd0) and all 23 Dij coefficients. 
**************************************************************************
*     Author: Gunnar Klaemke
*     Date: 07.06.06
**************************************************************************

      IMPLICIT NONE
      REAL*8 mass, q1sq, q2sq, q3sq, q1q2, q1q3, q2q3
      COMPLEX*16 DDij(1:24)
      INTEGER D_1M_GET
      INTEGER key,i

#include "d_functions.inc"
      
      key = D_1M_GET(mass, q1sq, q2sq, q3sq, q1q2, q1q3, q2q3)
      
      do i=1,24
         DDij(i)=Cache_values(i,key)
      enddo
      
      RETURN
      END

**************************************************************************
      SUBROUTINE D_1M_ALL2(mass,q1sq,q2sq,q3sq,q1q2,q1q3,q2q3,Dij)
**************************************************************************
*     Returns the Passarino-Veltman Dij 4-point tensor coefficients for 
*     one mass in old convetion --> DDij(3,13) without D0
*
*     Input: 
*     mass            (REAL*8) the mass of the particle in the loop
*     q1sq,q2sq,q3sq  (REAL*8) q1^2, q2^2 and q3^2 of external momenta
*     q1q2,q1q3,q2q3  (REAL*8) dot-products of external momenta
*     Output:
*     DDij(3,13)  (COMPLEX*16) all 23 Dij coefficients + 7 coeff for DD 
**************************************************************************
*     Author: Michael Kubocz
*     Date:   20.09.07
**************************************************************************

      IMPLICIT NONE

#include "d_functions.inc"

      real*8 mass,q1sq,q2sq,q3sq,q1q2,q1q3,q2q3
      complex*16 Dij(3,13)

      integer D_1M_GET
      integer key,i,j

      key=D_1M_GET(mass,q1sq,q2sq,q3sq,q1q2,q1q3,q2q3)
      
      do i=1,3
         do j = 1,13
            Dij(i,j)=DCMPLX(0d0,0d0)
c            print*,'Dij(i,j)',Dij(i,j)
         enddo
c         print*,''
      enddo

      Dij(1,1) =Cache_values(2 ,key)
      Dij(1,2) =Cache_values(3 ,key)
      Dij(1,3) =Cache_values(4 ,key)
      Dij(2,1) =Cache_values(5 ,key)
      Dij(2,2) =Cache_values(6 ,key)
      Dij(2,3) =Cache_values(7 ,key)
      Dij(2,4) =Cache_values(8 ,key)
      Dij(2,5) =Cache_values(9 ,key)
      Dij(2,6) =Cache_values(10,key)
      Dij(2,7) =-Cache_values(11,key)
      Dij(3,1) =Cache_values(12,key)
      Dij(3,2) =Cache_values(13,key)
      Dij(3,3) =Cache_values(14,key)
      Dij(3,4) =Cache_values(15,key)
      Dij(3,5) =Cache_values(16,key)
      Dij(3,6) =Cache_values(17,key)
      Dij(3,7) =Cache_values(18,key)
      Dij(3,8) =Cache_values(19,key)
      Dij(3,9) =Cache_values(20,key)
      Dij(3,10)=Cache_values(21,key)
      Dij(3,11)=-Cache_values(22,key)
      Dij(3,12)=-Cache_values(23,key)
      Dij(3,13)=-Cache_values(24,key)
c  Denner-Dittmaier coefficients
      Dij(1,7) =Cache_values(25,key) 
      Dij(1,8) =Cache_values(26,key)
      Dij(1,9) =Cache_values(27,key)
      Dij(1,10)=Cache_values(28,key)
      Dij(1,11)=Cache_values(29,key)
      Dij(1,12)=Cache_values(30,key)
      Dij(1,13)=Cache_values(31,key)

      RETURN
      END

**************************************************************************
      SUBROUTINE D_1M_ALL3(mass,q1sq,q2sq,q3sq,q1q2,q1q3,q2q3,Dij)
**************************************************************************
*     Returns the Passarino-Veltman Dij 4-point tensor coefficients for 
*     one mass in old convetion, but with mirrored arrays --> DDij(13,3) 
*     without D0
*
*     Input: 
*     mass            (REAL*8) the mass of the particle in the loop
*     q1sq,q2sq,q3sq  (REAL*8) q1^2, q2^2 and q3^2 of external momenta
*     q1q2,q1q3,q2q3  (REAL*8) dot-products of external momenta
*     Output:
*     DDij(13,3)  (COMPLEX*16) all 23 Dij coefficients + 7 coeff for DD 
**************************************************************************
*     Author: Michael Kubocz
*     Date:   20.09.07
**************************************************************************

      IMPLICIT NONE

#include "d_functions.inc"

      real*8 mass,q1sq,q2sq,q3sq,q1q2,q1q3,q2q3
      complex*16 Dij(13,3)

      integer D_1M_GET
      integer key,i,j

      key=D_1M_GET(mass,q1sq,q2sq,q3sq,q1q2,q1q3,q2q3)

      do i=1,13
         do j = 1,3
            Dij(i,j)=DCMPLX(0d0,0d0)
         enddo
      enddo

      Dij(1,1) =Cache_values(2 ,key)
      Dij(2,1) =Cache_values(3 ,key)
      Dij(3,1) =Cache_values(4 ,key)
      Dij(1,2) =Cache_values(5 ,key)
      Dij(2,2) =Cache_values(6 ,key)
      Dij(3,2) =Cache_values(7 ,key)
      Dij(4,2) =Cache_values(8 ,key)
      Dij(5,2) =Cache_values(9 ,key)
      Dij(6,2) =Cache_values(10,key)
      Dij(7,2) =-Cache_values(11,key)
      Dij(1,3) =Cache_values(12,key)
      Dij(2,3) =Cache_values(13,key)
      Dij(3,3) =Cache_values(14,key)
      Dij(4,3) =Cache_values(15,key)
      Dij(5,3) =Cache_values(16,key)
      Dij(6,3) =Cache_values(17,key)
      Dij(7,3) =Cache_values(18,key)
      Dij(8,3) =Cache_values(19,key)
      Dij(9,3) =Cache_values(20,key)
      Dij(10,3)=Cache_values(21,key)
      Dij(11,3)=-Cache_values(22,key)
      Dij(12,3)=-Cache_values(23,key)
      Dij(13,3)=-Cache_values(24,key)
c  Denner-Dittmaier coefficients
      Dij(7,1) =Cache_values(25,key) 
      Dij(8,1) =Cache_values(26,key)
      Dij(9,1) =Cache_values(27,key)
      Dij(10,1)=Cache_values(28,key)
      Dij(11,1)=Cache_values(29,key)
      Dij(12,1)=Cache_values(30,key)
      Dij(13,1)=Cache_values(31,key)

      END

