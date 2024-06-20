c***********************************************************************************
!****************subroutine contract_TZZ****************************************
!   Tinghua Chen, <txchen2@shockers.wichita.edu>
!   Last modified: 2023, February
!
!  These subroutines calcluate the tensor structure of the HZZ.
!  The input should be the current J1,J2, momentum k1,k2,
!  output is the contraction of J1TJ2 with anomalous couplings a2,a3 
!  and the corresponding propagator
!  reference: Phys.Rev.D 74 (2006) 095001

        subroutine contract_TmnZZ(J1,J2,k1,k2,J1TJ2Z)
        implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/koppln.inc"   
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

! declaration of dummy arguments
        double complex J1(0:3), J2(0:3)  !input 
        double precision k1(0:4), k2(0:4) !output
        double complex J1TJ2Z
        
! local         
        double complex propzz, J1CJ2_E, J1CJ2_O 
! anomalous couplings, use for a2 ,a3
        double complex ahvv(3,4,4), ahvvL(3,4,4)
        common/tensorhvv/ ahvv, ahvvL
! external
        double precision  clr, xm2, xmg, b, v, a
        COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     &                  V(4,5),A(4,5)

    ! use the cpe and cpo subroutine 
        call contract_CPE(J1CJ2_E,J1,J2,k1,k2) 
        call contract_CPO(J1CJ2_O,J1,J2,k1,k2)
        propzz = (1.0d0
     &       /dcmplx(k1(4)-xm2(2),xmg(2))
     &       /dcmplx(k2(4)-xm2(2),xmg(2)))

        J1TJ2Z = (ahvv(2,2,2) * J1CJ2_E + ahvv(3,2,2) * J1CJ2_O) * propzz
        return
        end

c***********************************************************************************
!****************subroutine contract_TWW****************************************
!   
!  These subroutines calcluate the tensor structure of the HWW.
!  The input should be the current J1,J2, momentum k1,k2,
!  output is the contraction of J1TJ2 with anomalous couplings a2,a3 
!  and the corresponding propagator
!  reference: Phys.Rev.D 74 (2006) 095001

        subroutine contract_TmnWW(J1,J2,k1,k2,J1TJ2W)
        implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/koppln.inc"   
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

! local
        double complex J1(0:3), J2(0:3), J1CJ2_E,J1CJ2_O !input
        double complex J1TJ2W
        double precision k1(0:4), k2(0:4) 
        double complex propww
! anomalous couplings, use for a2 ,a3
        double complex ahvv(3,4,4), ahvvL(3,4,4)
        common/tensorhvv/ ahvv, ahvvL
! external
        double precision  clr, xm2, xmg, b, v, a
        COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     &                  V(4,5),A(4,5)

! use the cpe and cpo subroutine 
        call contract_CPE(J1CJ2_E,J1,J2,k1,k2) 
        call contract_CPO(J1CJ2_O,J1,J2,k1,k2)
         propww = (1.0d0
     &      /dcmplx(k1(4)-xm2(3),xmg(3))
     &      /dcmplx(k2(4)-xm2(3),xmg(3)))

        J1TJ2W = (ahvv(2,3,4) * J1CJ2_E + ahvv(3,3,4) * J1CJ2_O) *propww
        return
        end
c****************************hjjja_virt*****************************************
c        This subroutine calcualte the a2 term for the anomalous hvv 
c        J1T = (g^(mu,nu))k1.k2- k1^nu k2^mu)J1^mu
c        input: doule complex J1 (0:3) 4-current, 
c        double precision k1(0:3), k2(0:3) 4 momenta   
c        output: double complex J1T = (k1.k2)J1-k1.(k2.J1)
         
        subroutine contract_CPEJ1(J1T,J1,k1,k2)
                implicit none
     
        double complex J1(0:3), J1k2, J1T(0:3)
        double precision k1(0:3), k2(0:3),k12
c local variables
        integer mu,nu         
c dot products
        k12 = k1(0)*k2(0)-k1(1)*k2(1)-k1(2)*k2(2)-k1(3)*k2(3)
        J1k2 = J1(0)*k2(0)-J1(1)*k2(1)-J1(2)*k2(2)-J1(3)*k2(3)
     
        J1T(0) = k12*J1(0) - k1(0)*J1k2
        J1T(1) = k12*J1(1) - k1(1)*J1k2
        J1T(2) = k12*J1(2) - k1(2)*J1k2
        J1T(3) = k12*J1(3) - k1(3)*J1k2
!        J1T(:) = J1T(:)-k1*J1k2
        return
        end


c****************************hjjja_virt*****************************************
c        This subroutine calcualte the a2 term for the anomalous hvv 
c        J2T = (g^(mu,nu))k1.k2- k1^nu k2^mu)J1^mu
c        input: doule complex J1 (0:3) 4-current, 
c        double precision k1(0:3), k2(0:3) 4 momenta   
c        output: double complex J1T = (k1.k2)J1-k1.(k2.J1)
         
        subroutine contract_CPEJ2(J2T,J2,k1,k2)
        implicit none

        double complex J2(0:3), J2k1, J2T(0:3)
        double precision k1(0:3), k2(0:3),k12
c local variables
        integer mu,nu         
c dot products
        k12 = k1(0)*k2(0)-k1(1)*k2(1)-k1(2)*k2(2)-k1(3)*k2(3)
        J2k1 = J2(0)*k1(0)-J2(1)*k1(1)-J2(2)*k1(2)-J2(3)*k1(3)

        J2T(0) = k12*J2(0) - k2(0)*J2k1
        J2T(1) = k12*J2(1) - k2(1)*J2k1
        J2T(2) = k12*J2(2) - k2(2)*J2k1
        J2T(3) = k12*J2(3) - k2(3)*J2k1
!        J1T(:) = J1T(:)-k1*J1k2
        return
        end
