********************************************************************************
********************************************************************************
*** FC_tools.F                                                               ***
*** 27 August 2012                                                           ***
*** sophy@particle.uni-karlsruhe.de                                          ***
***                                                                          ***
*** This set of subroutine and functions is stolen and adapted from FormCalc ***
*** and includes the kinetic functions necessary to use FC matrix elements   ***
*** in VBFNLO.                                                               ***
***   (FormCalc is by Thomas Hahn, www.formcalc.de)                          ***
***                                                                          ***
********************************************************************************
********************************************************************************

*** This subroutine fulfils the same role as VecSet does in FormCalc - i.e.
*** it calculates the values of the necessary momentum and polarization vectors

      subroutine FC_kin(i, mass, mom, momx, momy, momz)

      implicit none


*** ARGUMENTS:
***    i = particle number (FC labelling!)
***    m = mass of particle i
***    p = absolute value of 3-momentum of particle i
***    ex, ey, ez = unit 3-vector of momentum of particle i
      integer i
      double precision mass, mom, momx, momy, momz


      double precision sqrt2
      parameter (sqrt2 = 1.4142135623730950488016887242096981D0)


** Vector K, from VecSet.F and associated dummy variables
      double complex KK(5,2,2)
      double precision q0, deltaq, sinth, sumq
      double precision oneQez, oneMez
      double complex expIphi

** Spinor and vector values and labels.  Final components of spi/vec are spinor 
** IDs
      double complex spi(2,98)
      double complex vec(2,2,41)
      integer x, y

      common /FCvectors/ spi, vec

** Vector label
      integer k
      external k


** Do loop variable, mu = spacetime


** Debugging parameter
      integer debug


** Setting KK(i,x,y) and vec(x,y,k(i))  - see VecSet.F  
      if( mass**2 .lt. 1D-14 ) then
         q0 = mom
         deltaq = 0
      else
         q0 = sqrt(mom**2 + mass**2)
         deltaq = (mass**2)/(q0+mom)
      end if
      
      sinth = momx**2 + momy**2
      oneQez = 1 + momz
      if(oneQez .lt. 0.5D0) then
         oneQez = sinth/(1 - momz)
      endif
      oneMez = 1 - momz
      if(oneMez .lt. 0.5D0) then 
         oneMez = sinth/(1 + momz)      
      endif
      

      KK(i,1,1) = q0*oneQez - deltaq*momz
      KK(i,2,2) = q0*oneMez + deltaq*momz
      KK(i,2,1) = mom*DCMPLX(momx, momy)
      KK(i,1,2) = DCONJG(KK(i,2,1))
      
      do x = 1, 2
         do y = 1, 2
            vec(x,y,k(i)) = KK(i,x,y)
         end do
      end do
         

** Now the spinors and their dummy variables.  This is utterly stolen from 
** VecSet.F
      if( sinth .lt. 1D-14 ) then
         expIphi = 1
      else
         sinth = sqrt(sinth)
         expIphi = DCMPLX(momx, momy)/sinth
      endif

      expIphi = 1/sqrt2*sqrt(oneQez)*expIphi
      sinth = 1/sqrt2*sqrt(oneMez)
      sumq = sqrt(q0 + mom)
      deltaq = sqrt(deltaq)


** The spinors themselves.  For labelling conventions, see notebook 6, page 110
*        if (i .eq. 1 .or. i .eq. 2) then

      spi(1,9+(16*i)) = deltaq*sinth
      spi(2,9+(16*i)) = -deltaq*expIphi
      
      spi(1,10+(16*i)) = sumq*sinth
      spi(2,10+(16*i)) = -sumq*expIphi
      
      spi(1,13+(16*i)) = DCONJG(sumq*expIphi)
      spi(2,13+(16*i)) = sumq*sinth
      
      spi(1,14+(16*i)) = DCONJG(deltaq*expIphi)
      spi(2,14+(16*i)) = deltaq*sinth

*        end if

** Dotted spinors 
*        if (i .eq. 3 .or. i .eq. 5) then

      spi(1,11+(16*i)) = deltaq*sinth
      spi(2,11+(16*i)) = DCONJG(-deltaq*expIphi)
      
      spi(1,12+(16*i)) = sumq*sinth
      spi(2,12+(16*i)) = DCONJG(-sumq*expIphi)
      
      spi(1,15+(16*i)) = sumq*expIphi
      spi(2,15+(16*i)) = sumq*sinth
      
      spi(1,16+(16*i)) = deltaq*expIphi
      spi(2,16+(16*i)) = deltaq*sinth
      
*        end if


** Debugging messages
      if (debug .eq. 2) then
         
         write(*,*)'     '
         write(*,*)'i = ', i

         write(*,*)'mom =',  mom
         write(*,*)'mom e =', momx, momy, momz

         write(*,*)'KK(i,1,1) =',  KK(i,1,1)
         write(*,*)'KK(i,1,2) =',  KK(i,1,2)
         write(*,*)'KK(i,2,1) =',  KK(i,2,1)
         write(*,*)'KK(i,2,2) =',  KK(i,2,2)
      endif

      
      end


*******************************************************************************
*******************************************************************************

      subroutine FCinvariants(SInv, TInv)

      implicit none

#include "global.inc"

** Spinor and vector values and labels
      double complex spi(2,98)
      double complex vec(2,2,41)

      common /FCvectors/ spi, vec   


** do loops - run over particle numbers
      integer i, j

** The invariants SInv and TInv, and their associated dummy variables
      double precision S11, S22, T11, T22
      double complex S12, T12
      double precision SInv(n_p+n_v,n_p+n_v), TInv(n_p+n_v,n_p+n_v)

      integer k
      external k


** The invariants
      do i = 1, (n_p+n_v)
         do j = 1, (n_p)+(n_v)

** Now the SInvariant - see SInvariant.F (S_ij = (k_i + k_j)^2
            S11 = DBLE(vec(1,1,k(i))) + DBLE(vec(1,1,k(j)))
            S22 = DBLE(vec(2,2,k(i))) + DBLE(vec(2,2,k(j)))
            S12 = vec(1,2,k(i)) + vec(1,2,k(j))

            SInv(i,j) = S11*S22 - DBLE(S12)**2 - DIMAG(S12)**2

** And the TInvariant - see TInvariant.F (T_ij = (k_i - k_j)^2)
            T11 = DBLE(vec(1,1,k(i))) - DBLE(vec(1,1,k(j)))
            T22 = DBLE(vec(2,2,k(i))) - DBLE(vec(2,2,k(j)))
            T12 = vec(1,2,k(i)) - vec(1,2,k(j))
      
            TInv(i,j) = T11*T22 - DBLE(T12)**2 - DIMAG(T12)**2  

         end do
      end do


** Debugging messages
c$$$      do i = 1, 5
c$$$         write(*,*)'     '
c$$$         write(*,*)'i, k(i) = ', i, k(i)
c$$$
c$$$         write(*,*)'vec(1,1,k(i)) =', vec(1,1,k(i))
c$$$         write(*,*)'vec(1,2,k(i)) =', vec(1,2,k(i))
c$$$         write(*,*)'vec(2,2,k(i)) =', vec(2,2,k(i))
c$$$
c$$$         do j = 1, 5
c$$$            write(*,*) 'j =', j
c$$$
c$$$            write(*,*)'vec(1,1,k(i)) =', vec(1,1,k(i))
c$$$            write(*,*)'vec(1,2,k(i)) =', vec(1,2,k(i))
c$$$            write(*,*)'vec(2,2,k(i)) =', vec(2,2,k(i))
c$$$         end do
c$$$      end do
c$$$
c$$$      do i = 1, 5
c$$$         write(*,*)'     '
c$$$         write(*,*)'i = ', i
c$$$         
c$$$         do j = 1, 5
c$$$            write(*,*) 'j =', j
c$$$            write(*,*)'SInv(i,j) =', SInv(i,j)
c$$$            write(*,*)'TInv(i,j) =', TInv(i,j)
c$$$         end do
c$$$      end do


      end


*******************************************************************************
*******************************************************************************

*** This function returns the integer id of a particular vector.  See util.h of
*** FormCalc.

      integer function k(i)
      implicit none

      integer i

      k = 8*i + 1

      end


*******************************************************************************
*******************************************************************************

*** This function returns the integer id of a (polarization) vector.  See 
*** util.h of FormCalc.

      integer function e(i)
      implicit none

      integer Hel(5)
      common /FChelicities/ Hel

      integer i

      e = 8*i + 3 + Hel(i)

      end


*******************************************************************************
*******************************************************************************

*** This function returns the integer id of a conujugate (polarization) vector. 
*** See util.h of FormCalc.

      integer function ec(i)
      implicit none

      integer Hel(5)
      common /FChelicities/ Hel

      integer i

      ec = 8*i + 3 - Hel(i)

      end


******************************************************************************
******************************************************************************

*** This function calculates the scalar product of two 4-vectors.  

      double complex function vPair(a,b)

      implicit none

** The arguments a, b, of the function specify the vectors' position in the 
** vector array 'vec'
      integer a, b

** Spinor and vector values and labels.  Final components of spi/vec are spinor 
** IDs
      double complex spi(2,98)
      double complex vec(2,2,41)

      common /FCvectors/ spi, vec

** Dummy variables
      integer x, y
      double complex AA(2,2), BB(2,2)


      do x = 1, 2
         do y = 1, 2
            AA(x,y) = vec(x,y, a)
            BB(x,y) = vec(x,y, b)
         end do
      end do

      vPair = 0.5D0*(AA(1,1)*BB(2,2) + AA(2,2)*BB(1,1) -
     &     AA(1,2)*BB(2,1) - AA(2,1)*BB(1,2) )
      

      end


*******************************************************************************
*******************************************************************************

*** This function returns the integer id of a particular dotted spinor.  See 
*** util.h of FormCalc.

      integer function DottedSpinor(i, s, om)
      implicit none

      integer i, s, om

      integer Hel(5)
      common /FChelicities/ Hel

      DottedSpinor = 2*s*Hel(i) + (16*i) + om + 7

      end


*******************************************************************************
*******************************************************************************

*** This function returns the integer id of a particular spinor.  See util.h of
*** FormCalc.

      integer function Spinor(i, s, om)
      implicit none

      integer i, s, om

      integer Hel(5)
      common /FChelicities/ Hel

      Spinor = 2*s*Hel(i) + (16*i) + om + 5

      end


*******************************************************************************
*******************************************************************************

*** This function computes the value of Spinor.eps.Spinor.  Stolen from 
*** util/kin/SeS.F of FormCalc.

      double complex function vSeS(a, b)
      implicit none

** Labels - a and b are the integer id's of the two spinors.
      integer a, b, x

** Spinor and vector values and labels.  Final components of spi/vec are spinor 
** IDs
      double complex spi(2,98)
      double complex vec(2,2,41)

      common /FCvectors/ spi, vec

** Dummy variables
      double complex AA(2), BB(2)


      do x = 1,2
         AA(x) = spi(x, a)
         BB(x) = spi(x, b)
      end do
      
      vSeS = AA(1)*BB(2) - AA(2)*BB(1)

      
      end


*******************************************************************************
*******************************************************************************

*** This function computes the value of Spinor.Spinor.  Stolen from 
*** util/kin/SxS.F of FormCalc.

      double complex function vSxS(a, b)
      implicit none

** Labels - a and b are the integer id's of the two spinors.
      integer a, b, x

** Spinor and vector values and labels.  Final components of spi/vec are spinor 
** IDs
      double complex spi(2,98)
      double complex vec(2,2,41)

      common /FCvectors/ spi, vec

** Dummy variables
      double complex AA(2), BB(2)


      do x = 1, 2
         AA(x) = spi(x, a)
         BB(x) = spi(x, b)
      end do

      vSxS = AA(1)*BB(1) + AA(2)*BB(2)
      
      
      end


*******************************************************************************
*******************************************************************************

*** This function computes the value of Vector.eps.Spinor by defining a new 
*** spinor.  Stolen from util/kin/VeS.F of FormCalc.

      integer function vVeS(a, b)
      implicit none

** Labels - a and b are the integer id's of the two spinors.
      integer a, b, x, y

** Spinor and vector values and labels.  Final components of spi/vec are spinor 
** IDs
      double complex spi(2,98)
      double complex vec(2,2,41)

      common /FCvectors/ spi, vec

** Dummy variables
      double complex AA(2,2), BB(2)


      do x = 1, 2
         do y = 1, 2
            AA(x,y) = vec(x,y, a)
         end do
         BB(x) = spi(x, b)
      end do
      
      spi(1, 1) = AA(1,1)*BB(2) - AA(1,2)*BB(1)
      spi(2, 1) = AA(2,1)*BB(2) - AA(2,2)*BB(1)
      vVeS = 1
      
      
      end


*******************************************************************************
*******************************************************************************

*** This function computes the value of Vector.Spinor by defining a new spinor.
*** Stolen from util/kin/VxS.F of FormCalc.

      integer function vVxS(a, b)
      implicit none

** Labels - a and b are the integer id's of the two spinors.
      integer a, b, x, y

** Spinor and vector values and labels.  Final components of spi/vec are spinor 
** IDs
      double complex spi(2,98)
      double complex vec(2,2,41)

      common /FCvectors/ spi, vec

** Dummy variables
      double complex AA(2,2), BB(2)


      do x = 1, 2
         do y = 1, 2
            AA(x,y) = vec(x,y, a)
         end do
         BB(x) = spi(x, b)
      end do

      spi(1, 1) = AA(1,1)*BB(1) + AA(1,2)*BB(2)
      spi(2, 1) = AA(2,1)*BB(1) + AA(2,2)*BB(2)
      vVxS = 1
      

      end


*******************************************************************************
*******************************************************************************

*** This function computes the value of Bar[Vector].eps.Spinor by defining a 
*** new spinor.  Stolen from util/kin/BeS.F of FormCalc. 

      integer function vBeS(a, b)
      implicit none

** Labels - a and b are the integer id's of the two spinors.
      integer a, b, x, y

** Spinor and vector values and labels.  Final components of spi/vec are spinor 
** IDs
      double complex spi(2,98)
      double complex vec(2,2,41)

      common /FCvectors/ spi, vec

** Dummy variables
      double complex AA(2,2), BB(2)


      do x = 1, 2
         do y = 1, 2
            AA(x,y) = vec(x,y, a)
         end do
         BB(x) = spi(x, b)
      end do

      spi(1, 2) = AA(1,2)*BB(1) + AA(2,2)*BB(2)
      spi(2, 2) = -AA(1,1)*BB(1) - AA(2,1)*BB(2)
      vBeS = 2
      

      end


*******************************************************************************
*******************************************************************************

*** This function computes the value of Bar[Vector].Spinor by defining a new 
*** spinor.  Stolen from util/kin/BxS.F of FormCalc.

      integer function vBxS(a, b)
      implicit none

** Labels - a and b are the integer id's of the two spinors.
      integer a, b, x, y

** Spinor and vector values and labels.  Final components of spi/vec are spinor 
** IDs
      double complex spi(2,98)
      double complex vec(2,2,41)

      common /FCvectors/ spi, vec

** Dummy variables
      double complex AA(2,2), BB(2)

        
      do x = 1, 2
         do y = 1, 2
            AA(x,y) =  vec(x,y, a)
         end do
         BB(x) = spi(x, b)
      end do

      spi(1, 2) = AA(2,2)*BB(1) - AA(1,2)*BB(2)
      spi(2, 2) = AA(1,1)*BB(2) - AA(2,1)*BB(1)
      vBxS = 2


      end


*******************************************************************************
*******************************************************************************

*** Function stolen from SquaredME.F of FormCalc to calculate tree*loop.

      double precision function sumup(MatSUN, CCloop, CCtree)
      implicit none

      double complex CCloop(1), CCtree(1), m
      
      double precision MatSUN
      
      sumup = 0
      
      m = 0
      m = m + CCloop(1)*MatSUN
      sumup = sumup + DBLE(DCONJG(CCtree(1))*m)
      
      end


*******************************************************************************
*******************************************************************************

*** Function stolen from SquaredME.F of FormCalc to calculate tree*loop.

      double precision function sumup4(MatSUN, CCloop, CCtree)
      implicit none
      
      integer jSUNloop, jSUNtree
      double complex CCloop(2), CCtree(2)
      double complex m
      double complex MatSUN(2,2)

      sumup4 = 0d0

      do jSUNtree = 1, 2
         m = (0d0,0d0)
         do jSUNloop = 1, 2
            m = m + CCloop(jSUNloop)*MatSUN(jSUNloop,jSUNtree)
         enddo
         sumup4 = sumup4 + DBLE(DCONJG(CCtree(jSUNtree))*m)
      enddo

      end


*******************************************************************************
