*******************************************************************
*  file for the cached C-functions (improved version)
*******************************************************************
*  Content: the cached C0- & Cij-Function and some service routines
*
*     SUBROUTINE C_1M_INIT_CACHE()
*     SUBROUTINE C_1M_CLEAR_CACHE()
*     C_1M_GET(mass,q1sq,q2sq,q12sq)
*     C_1M_VAL(coefficient, key)
*     C0_1M(mass,q1sq,q2sq,q12sq)
*     Cij_1M(coefficient,mass,q1sq,q2sq,q12sq)
*
*******************************************************************
*  Main dependencies: 
*     amp_aux.f     | I3point, tens3, B0tM 
*******************************************************************
*  Author: Gunnar Klaemke
*  Date: 08.10.07
*******************************************************************

**************************************************************************
      SUBROUTINE C_1M_INIT_CACHE()
**************************************************************************
*     Initializes the cache for the C_1M functions.
**************************************************************************
*     Author: Gunnar Klaemke
*     Date: 08.10.07
**************************************************************************

      IMPLICIT NONE
      
#include "c_functions.inc"

      INTEGER i 

      Cache_index=1
      Cache_count=0
      CallCounter=0
      
      do i=1,LT_Size
         LookupTable(i) = 0
      enddo

c      print *,""
c      print *," * Initialized cache for C_1M function"
c      print *,""
      RETURN
      END


**************************************************************************
      SUBROUTINE C_1M_CLEAR_CACHE()
**************************************************************************
*     Resets the cache for the C_1M functions.
**************************************************************************
*     Author: Gunnar Klaemke
*     Date: 08.10.07
**************************************************************************

      IMPLICIT NONE
      
#include "c_functions.inc"

      Cache_index=1
      Cache_count=0
      CallCounter=0

      RETURN
      END


**************************************************************************
      FUNCTION C_1M_GET(mass,q1sq,q2sq,q12sq)
**************************************************************************
*     Calculates the Passarino-Veltman C0-function and the Cij tensor 
*     coefficients for 1 mass in the loop. The results are stored in a 
*     cache, unless the results have already been calculated. Also 
*     calculates the results for the permutation q1<->q2 and stores them 
*     in the cache. The function returns the Index at which the arguments
*     and results are stored as a key.
*     Uses 'I3point', 'B0tM' and 'tens3' for numerical evaluations.
*     Input: 
*     mass        (REAL*8) the mass of the particle in the loop
*     q1sq,q2sq   (REAL*8) q1^2 and q2^2 of external momenta
*     q12sq       (REAL*8) q12^2=(q1+q2)^2 of external momenta
*     Output:
*     c_1m_get:  (INTEGER) a key for C_1M_VAL 
**************************************************************************
*     Author: Gunnar Klaemke
*     Date: 08.10.07
**************************************************************************

      IMPLICIT NONE
      
      INTEGER C_1M_GET
      REAL*8 mass, q1sq, q2sq, q12sq 

      REAL*8 arg(4),qsq
      CHARACTER*(8*4) arg_data
      INTEGER*4 arg_bits(2*4)
*     use strings for faster comparison of keys with arguments
      EQUIVALENCE (arg_data,arg)  
      EQUIVALENCE (arg_bits,arg)  
*     for numerical computation of the c0-function
      COMPLEX*16 I3point, B0tM
      EXTERNAL I3point, B0tM
*     this is the position in the cache, where the last key was found in a search loop
      INTEGER last_lookup
      SAVE last_lookup

#include "c_functions.inc"
#include "c_coefficients.inc"

      complex*16 B0(3), C0, Cij(4,2)
      real*8 C0r,C0I,Cij_Re(4,2),Cij_Im(4,2)
      INTEGER l_ch, ll
*     the number of ignored bits defines the number of binary digits that 
*     shall be ignored in comparisons
      INTEGER*4 IGNORED_BITS
      PARAMETER (IGNORED_BITS=2**10-1)   ! 10 bits = 3 digits


*  - START -

*     fill input arguments into an array
      arg(1)=mass
      arg(2)=q1sq
      arg(3)=q2sq
      arg(4)=q12sq
*     overwrite the least significant bits of each argument with 0
*     this allows to disregard the last digits in comparisons,
*     because 2 arguments might differ only in their last digits (for numerical reasons).
*     A naive comparison would consider those as not equal, though they actually are quite equal.
      arg_bits(1)=AND(arg_bits(1),NOT(IGNORED_BITS))
      arg_bits(3)=AND(arg_bits(3),NOT(IGNORED_BITS))
      arg_bits(5)=AND(arg_bits(5),NOT(IGNORED_BITS))
      arg_bits(7)=AND(arg_bits(7),NOT(IGNORED_BITS))

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
               C_1M_GET = ll    ! return index
               RETURN
            endif  
         endif
         
*     the lookup table doesn't give the right answer -> the cache has to be searched:
        
*     begin search at the position where the last key was found ('last_lookup')

         DO l_ch=1,Cache_count       ! loop over all entries
            if (Cache_keys_data(last_lookup).eq.arg_data) then ! arguments are already stored in cache
               LookupTable(CallCounter) = last_lookup ! save position in look-up table
               C_1M_GET = last_lookup ! return index
               RETURN
            endif  
            if (last_lookup.eq.Cache_count) last_lookup=0   ! jump to start, if end is reached
            last_lookup=last_lookup+1
         ENDDO
        
      endif    

*     the arguments were not found in the cache:
*     evaluate cij coefficients
      
      B0(1)=B0tM(mass,q1sq)
      B0(2)=B0tM(mass,q2sq)
      B0(3)=B0tM(mass,q12sq)

      C0 = I3point(mass,q1sq,q2sq,q12sq) 
      CALL tens_red3_new_Re_Com_1M(mass,q1sq,q2sq,q12sq,B0(2),B0(3),B0(1),
     &                               C0,C0r,C0I,Cij_Re,Cij_Im)

c... combine real and imaginary parts of Cij
      Cij(1,1)=DCMPLX(Cij_Re(1,1),Cij_Im(1,1))
      Cij(2,1)=DCMPLX(Cij_Re(2,1),Cij_Im(2,1))

      Cij(1,2)=DCMPLX(Cij_Re(1,2),Cij_Im(1,2))
      Cij(2,2)=DCMPLX(Cij_Re(2,2),Cij_Im(2,2))
      Cij(3,2)=DCMPLX(Cij_Re(3,2),Cij_Im(3,2))
      Cij(4,2)=DCMPLX(Cij_Re(4,2),Cij_Im(4,2))

*     fill key & values in the cache
      Cache_keys_data(Cache_index)=arg_data
      Cache_values(cc0 ,Cache_index)=-C0
      Cache_values(cc11,Cache_index)=-Cij(1,1)
      Cache_values(cc12,Cache_index)=-Cij(2,1)

      Cache_values(cc21,Cache_index)=-Cij(1,2)
      Cache_values(cc22,Cache_index)=-Cij(2,2)
      Cache_values(cc23,Cache_index)=-Cij(3,2)
      Cache_values(cc24,Cache_index)=Cij(4,2)

*     return this index as the functions result
      C_1M_GET = Cache_index
*     reset last look_up
      last_lookup = Cache_index
*     fill index into lookup table
      LookupTable(CallCounter) = Cache_index
      
*     increment number of entries in cache
      if (Cache_count.lt.Cache_Size) Cache_count=Cache_count+1      
*     increment current cache position, jump to start if end is reached
      if (Cache_index.eq.Cache_Size) then 
         Cache_index=0
         print *," Warning: Cache overflow in C-functions!!"
      endif
      Cache_index=Cache_index+1

      RETURN  ! ?????????????????????
*     repeat this for the permutation q1<->q2
      qsq=arg(3) ! permute
      arg(3)=arg(2)
      arg(2)=qsq

*     fill key & values in the cache, use symmetry relations
      Cache_keys_data(Cache_index)=arg_data
      Cache_values(cc0 ,Cache_index)=C0
      Cache_values(cc11,Cache_index)=-Cij(2,1)-C0
      Cache_values(cc12,Cache_index)=-Cij(1,1)-C0
      Cache_values(cc21,Cache_index)=Cij(2,2) + 2*Cij(2,1) + C0
      Cache_values(cc22,Cache_index)=Cij(1,2) + 2*Cij(1,1) + C0
      Cache_values(cc23,Cache_index)=Cij(3,2) + Cij(1,2) + Cij(1,1) + C0
      Cache_values(cc24,Cache_index)=Cij(4,2)

*     increment number of entries in cache
      if (Cache_count.lt.Cache_Size) Cache_count=Cache_count+1      
*     increment current cache position, jump to start if end is reached
      if (Cache_index.eq.Cache_Size) then 
         Cache_index=0
         print *," Warning: Cache overflow in C-functions!!"
      endif
      Cache_index=Cache_index+1
      
*     finished
      RETURN
      END


**************************************************************************
      FUNCTION C_1M_VAL(coefficient, key)
**************************************************************************
*     Returns the Passarino-Veltman 3-point coefficients for a given key.
*     Use 'C_1M_GET' to obtain a key for some set of momenta.
*     Use the following constants for the coefficient:
*     cc0 : C0 scalar Integral
*     cc11, cc12, cc21, cc22, cc23, cc24 : Cij tensor coefficients 
*
*     Use the function in the following way:
*     ...
*     mykey =C_1M_GET(my_mass,my_q1sq,my_q2sq,my_q12sq)
*     my_C0 =C_1M_VAL(cc0 ,key)
*     my_C12=C_1M_VAL(cc12,key)
*     my_C24=C_1M_VAL(cc24,key)
*     ...
*     Important: Use different keys for different sets of momenta:
*     ...
*     key_A=C_1M_GET(mass_A,q1sq_A,q2sq_A,q12sq_A)
*     key_B=C_1M_GET(mass_B,q1sq_B,q2sq_B,q12sq_B)
*     ...
*
*     The key contains a validation tag, which is checked first. If a key
*     has become invalid (because of overwritten cache-entries), the cache
*     is searched again. (the first 16 bits of the key are the index, and
*     the next 16 bits are for the validation tag)
*
**************************************************************************
*     Author: Gunnar Klaemke
*     Date: 08.10.07
**************************************************************************
      
      IMPLICIT NONE
      COMPLEX*16 C_1M_VAL
      INTEGER coefficient, key

#include "c_functions.inc"
#include "c_coefficients.inc"

      C_1M_VAL=Cache_values(coefficient,key)

      RETURN
      END


**************************************************************************
      FUNCTION C0_1M(mass,q1sq,q2sq,q12sq)
**************************************************************************
*     Returns the Passarino-Veltman C0 function for one mass in the loop.
**************************************************************************
*     Author: Gunnar Klaemke
*     Date: 08.10.07
**************************************************************************

      IMPLICIT NONE
      COMPLEX*16 C0_1M
      REAL*8 mass, q1sq, q2sq, q12sq 
      INTEGER C_1M_GET
      COMPLEX*16 C_1M_VAL
      INTEGER key

#include "c_coefficients.inc"

      key=C_1M_GET(mass,q1sq,q2sq,q12sq)
      C0_1M=C_1M_VAL(cc0, key)
      
      RETURN
      END


**************************************************************************
      FUNCTION Cij_1M(coefficient,mass,q1sq,q2sq,q12sq)
**************************************************************************
*     Returns the Passarino-Veltman Cij 3-point tensor coefficients for 
*     one mass in the loop.
*     Possible coefficients:  cc11, cc12, cc21, cc22, cc23, cc24
**************************************************************************
*     Author: Gunnar Klaemke
*     Date: 08.10.07
**************************************************************************

      IMPLICIT NONE
      COMPLEX*16 Cij_1M
      REAL*8 mass, q1sq, q2sq, q12sq 
      INTEGER coefficient
      INTEGER C_1M_GET
      COMPLEX*16 C_1M_VAL
      INTEGER key

#include "c_coefficients.inc"

      key=C_1M_GET(mass,q1sq,q2sq,q12sq)
      Cij_1M=C_1M_VAL(coefficient, key)
      
      RETURN
      END
