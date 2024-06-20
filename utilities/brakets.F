****************************************************************************
*  file for spinors, currents and polarization vectors
****************************************************************************
*   This file contains functions for calculations with spinors (bra
*   and ket), current 4- and 6-vectors and polarization vectors. 
c***************************************************************************
c   LIST OF ALL FUNCTIONS AND SUBROUTINES IN THIS FILE:
c
c     SUBROUTINE CHIOM(NF,PBAR,SIGMA,CHI,OMFAC)
c     SUBROUTINE PSI0M( NF,PBAR,SIGN,PSI )
c     subroutine curr( sigmax,psibar,psi,j )
c     subroutine curr6( sigmax,psibar,pb,psi,p,j )
c     SUBROUTINE BRA2R(CHI,CHREAL,P,SIGMAP,K,REPS,RESULT,PPLUSK)
c     SUBROUTINE BRA2C( CHI,CHREAL,P,SIGMAP,K,CEPS,RESULT,PPLUSK )    (as ENTRY)
c     SUBROUTINE KET2R( CHI,CHREAL,P,SIGMAP,K,REPS,RESULT,PMINK )
c     SUBROUTINE KET2C( CHI,CHREAL,P,SIGMAP,K,CEPS,RESULT,PMINK )     (as ENTRY)
c     FUNCTION   S1R( BRA,RVEC,TIMEEX,SIGMA,KET )
c     FUNTCION   S1C( BRA,CVEC,TIMEEX,SIGMA,KET )     (as ENTRY)
c     SUBROUTINE POLVEC(KBAR,LAMBDA,EPSCAR)
c     SUBROUTINE HELVEC(KBAR,SIGMAK,LAMBDA,EPSHEL)    (as ENTRY)
c     FUNCTION   SC3(CHII,A1,A2,A3,CHIF,ALPHA)
c
c***************************************************************************
*   Last modified: 12.07.2006
****************************************************************************


C-------------------------------------------------------------------
C
c
      SUBROUTINE CHIOM(NF,PBAR,SIGMA,CHI,OMFAC)

      implicit none
C
C  Calculates the bra and ket vectors for external fermions as given
C  in Hagiwara, Zeppenfeld  MAD/PH/402, Appendix A, Eq.A.9
C         checked on May 24., 88
C
C  INPUT:
C  ------
C
C  NF            Number of fermions to calculate bras or kets for.
C                Fermions with odd index generate kets, even indexed
C                ones generate bras.
C
C  PBAR(0:3,NF)  Physical momenta of the NF fermions. All fermions
C                are assumed to be massless.
C
C  SIGMA(NF)     unbared helicity labels for the NF fermions
C                  SIGMA =  helicity*2 (+1 or -1) for fermions
C                  SIGMA = -helicity*2 (+1 or -1) for antifermions
C
C
C  OUTPUT:
C  -------
C
C  CHI(2,NF)     2 dimensional complex vector containing the ket |i>
C                for odd and the bra <i| for even numbered fermions
C
C  OMFAC         product of sqrt(2*pbar(0)) of all fermions entering
C                the normalization factor for amplitudes with massless
C                fermions
C
C
      INTEGER IF, NF, SIGMA(NF)
      real*8  PBAR(0:3,NF), OMFAC
      real*8  PABS, PBARX, PBARY, PBARZ, PAPZ, NORMAL, EPS
      complex*16    CHI(2,NF)
      PARAMETER (EPS=1.0d-30)
C
      OMFAC = 1.0d0
      DO 10 IF = 1,NF
        OMFAC = OMFAC* 2D0 *PBAR(0,IF)
        PABS  = PBAR(0,IF)
        PBARX = PBAR(1,IF)
        IF (MOD(IF,2).EQ.1) THEN
          PBARY = PBAR(2,IF)
        ELSE
          PBARY = -PBAR(2,IF)
        ENDIF
        PBARZ = PBAR(3,IF)
        IF (PBARZ.GT. 0.0d0 ) THEN
           PAPZ = PABS + PBARZ
        ELSE
           PAPZ = (PBARX**2+PBARY**2)/(PABS - PBARZ)
        ENDIF
C
C  Treat the pbar along negative z-axis case first
        IF(PAPZ .LE. EPS*PABS) THEN
          IF (SIGMA(IF).EQ.-1) THEN
            CHI(1,IF) = -1
            CHI(2,IF) = 0
          ELSEIF (SIGMA(IF).EQ.1) THEN
            CHI(1,IF) = 0
            CHI(2,IF) = 1
          ELSE
            WRITE(*,*) " Illegal fermion helicity SIGMA = ",SIGMA(IF)
          ENDIF
        ELSE
C
C And now the general case
          NORMAL = 1.0d0 / SQRT(2D0 *PABS*PAPZ)
          IF (SIGMA(IF).EQ.-1) THEN
            CHI(1,IF) = dcmplx(-NORMAL*PBARX,NORMAL*PBARY)
            CHI(2,IF) = NORMAL * PAPZ
          ELSEIF (SIGMA(IF).EQ.1) THEN
            CHI(1,IF) = NORMAL * PAPZ
            CHI(2,IF) = dcmplx(NORMAL*PBARX,NORMAL*PBARY)
          ELSE
            WRITE(*,*) " Illegal fermion helicity SIGMA = ",SIGMA(IF)
          ENDIF
        ENDIF
10    CONTINUE
      OMFAC = SQRT(OMFAC)
      RETURN
      END


C********************************************************************
C
C Calculates the bra and ket vectors for external massless fermions as
C given in Hagiwara, Zeppenfeld  MAD/PH/402, Appendix A, Eq.A.9
C checked on May 24., 88
C
C Modification on Aug. 14, 1992: 
C
C     The omega factors, i.e. sqrt(2E) are multiplied in already at 
C     the Chi level and returned are both helicity Weyl spinors
C     S*sqrt(2E)*chi_{+} and S*sqrt(2E)*chi_{-} 
C
C     In addition the subroutine was transformed to real*8
C
C  INPUT:
C  ------
C
C  NF            Number of fermions to calculate bras or kets for.
C                Fermions with odd index generate kets, even indexed
C                ones generate bras.
C
C  PBAR(0:3,NF)  Physical momenta of the NF fermions. All fermions
C                are assumed to be massless.
C
C  SIGN(NF)      The sign factors distinguishing fermion (S=+1) from
C                antifermion (S=-1) for all NF fermions
C
C
C  OUTPUT:
C  -------
C
C  PSI(2,-1:1,NF) 2 dimensional complex vector containing the ket |i>
C                 for odd and the bra <i| for even numbered fermions
C                 time S*sqrt(2E)
C
C********************************************************************
C
      SUBROUTINE PSI0M( NF,PBAR,SIGN,PSI )
      IMPLICIT NONE
      INTEGER  IF, NF, SIGN(NF)
      DOUBLE COMPLEX  PSI(2,-1:1,NF)
      DOUBLE PRECISION  PBAR(0:3,NF), PABS, PBARX, PBARY, PBARZ, PAPZ
      DOUBLE PRECISION  NORMAL, NX, NY, NZ, EPS
      PARAMETER ( EPS=1.d-30 )
CC
      DO IF = 1,NF
        PABS  = PBAR(0,IF)
        PBARX = PBAR(1,IF)
        IF ( MOD(IF,2).EQ.1 ) THEN
          PBARY =  PBAR(2,IF)
        ELSE
          PBARY = -PBAR(2,IF)
        ENDIF
        PBARZ = PBAR(3,IF)
        IF ( PBARZ.GT.0.d0 ) THEN
           PAPZ = PABS + PBARZ
        ELSE
           PAPZ = (PBARX**2 + PBARY**2)/(PABS - PBARZ)
        ENDIF
c
c treat the pbar along negative z-axis case first
c
        IF ( PAPZ.LE.EPS*PABS ) THEN
           NORMAL = SIGN(IF)*SQRT(2.d0*PBAR(0,IF))
           PSI(1,-1,IF) = -NORMAL 
           PSI(2,-1,IF) = 0.d0
           PSI(1,1,IF)  = 0.d0
           PSI(2,1,IF)  = NORMAL
        ELSE
c
c and now the general case
c
           NORMAL = SIGN(IF)/SQRT(PAPZ)
           NX = NORMAL*PBARX
           NY = NORMAL*PBARY
           NZ = NORMAL*PAPZ
           PSI(1,-1,IF) = DCMPLX(-NX,NY)
           PSI(2,-1,IF) = NZ
           PSI(1,1,IF)  = NZ
           PSI(2,1,IF)  = DCMPLX(NX,NY)
        ENDIF
      ENDDO
ccc
      END
C
C************************************  CURR  *************************
C
C  CURR calculates the current 
C
C        PSIBAR gamm^mu PSI
C 
C  of Hagiwara, Zeppenfeld, Nucl.Phys.B313 (1989) 560, eq. 2.20 for 
C  the two possible polarizations of the Weyl spinors PSIBAR and PSI
C
C  INPUT:
C  ------
C
C    PSIBAR(2,-1:1)  two Weyl spinors for 2 helicity states each
C    PSI(2,-1:1)
C
C  OUTPUT:
C  -------
C
C  J(0:3,-1:1)       the current as defined in eq.2.27 for two possible
C                    helicity combinations for massless fermions
C
      subroutine curr( sigmax,psibar,psi,j )

      implicit none

      double complex  j(0:3,-1:1), psibar(2,-1:1), psi(2,-1:1), zj2
      double complex  z1, z2, z3, z4
      integer  sigmax,sig
cc
      do sig = -1,sigmax,2
         z1 = psibar(1,sig) * psi(1,sig)
         z2 = psibar(2,sig) * psi(2,sig)
         z3 = psibar(1,sig) * psi(2,sig)
         z4 = psibar(2,sig) * psi(1,sig)
         j(0,sig) = z1 + z2
         if (sig.eq.-1) then
            j(1,sig) = -(z3+z4)
            zj2 = z3-z4
            j(2,sig) = dcmplx(-dimag(zj2),dreal(zj2))
            j(3,sig) = z2-z1
         else
            j(1,sig) = z3+z4
            zj2 = z4-z3
            j(2,sig) = dcmplx(-dimag(zj2),dreal(zj2))
            j(3,sig) = z1-z2
         endif
      enddo
ccc
      end
C
C************************************  CURR6  *************************
C
C  CURR6 calculates the current 
C
C        PSIBAR gamm^mu PSI
C 
C  of Hagiwara, Zeppenfeld, Nucl.Phys.B313 (1989) 560, eq. 2.20 for 
C  the two possible polarizations of the Weyl spinors PSIBAR and PSI
c  and stores the current momentum in the last 2 complex entries of the 
c  current, as in HELAS
C
C  INPUT:
C  ------
C
C    PSIBAR(2,-1:1)  two Weyl spinors for 2 helicity states each
C    PSI(2,-1:1)
c    pb(0:3)         the outflowing momentum of the bra spinor psibar
c    p(0:3)          the inflowing momentum of the ket spinor psi
C
C  OUTPUT:
C  -------
C
C  J(0:3,-1:1)       the current as defined in eq.2.27 for two possible
C                    helicity combinations for massless fermions
C  J(4:5,-1:1)       the momentum q = p-pb of the current, flowing away from 
c                    the fermion line stored as
c                    J(4) = dcmplx(q0,qz), J(5) = dcmplx(qx,qy)
C
      subroutine curr6( sigmax,psibar,pb,psi,p,j )

      implicit none

      double complex  j(0:5,-1:1), psibar(2,-1:1), psi(2,-1:1), zj2
      double precision pb(0:3), p(0:3)
      double complex  z1, z2, z3, z4, q4, q5
      integer  sigmax,sig
cc
      q4 = dcmplx(p(0)-pb(0),p(3)-pb(3))
      q5 = dcmplx(p(1)-pb(1),p(2)-pb(2))
      do sig = -1,sigmax,2
         z1 = psibar(1,sig) * psi(1,sig)
         z2 = psibar(2,sig) * psi(2,sig)
         z3 = psibar(1,sig) * psi(2,sig)
         z4 = psibar(2,sig) * psi(1,sig)
         j(0,sig) = z1 + z2
         if (sig.eq.-1) then
            j(1,sig) = -(z3+z4)
            zj2 = z3-z4
            j(2,sig) = dcmplx(-dimag(zj2),dreal(zj2))
            j(3,sig) = z2-z1
         else
            j(1,sig) = z3+z4
            zj2 = z4-z3
            j(2,sig) = dcmplx(-dimag(zj2),dreal(zj2))
            j(3,sig) = z1-z2
         endif
         j(4,sig) = q4
         j(5,sig) = q5
      enddo
ccc
      end
C
C
C******************************** bra2r, bra2c ************************
C
C   Calculates the bra vector <i,k| as given in
C      MAD/PH/402, Appendix A, Eqn A.10, top eqn.
C      checked on May 26, 88
C
C   Modified to double precision and nonzero time component of 
C   polarization vector Aug. 92
C
C   INPUT:
C            chi      double complex array(2)
C                        bra <i| as given by eq. A.9 if chreal = .true.
C                        any bra <...| if chreal = .false.
C
C            chreal  logical
C                        .true. if one component of chi is real
C                        .false. otherwise
C
C            p      double precision  array(0:3)
C                        standard four momentum for fermion
C                        p(0:3) = [E, px, py, pz]
C
C            sigmap      integer
C                        chirality/helicity factor for fermion
C                        allowed values: +1,-1
C
C            k      double precision  array(0:4)
C                        standard four momentum for boson and mass**2
C                        k(0:4) = [E, kx, ky, kz, mass**2]
C
C            reps     double precision  array(0:3)
C            ceps     double complex array(0:3)
C                     real or complex polarisation vector of the boson.
C
C   OUTPUT:
C            result      double complex array(1:2)
C                        two component row vector <i,k| on Eqn. A.10
C
C            pplusk      double precision  array(0:4)
C                        four momentum p+k. The fourth component
C                        contains the square (p+k)**2
C
C
      SUBROUTINE BRA2R(CHI,CHREAL,P,SIGMAP,K,REPS,RESULT,PPLUSK)

      implicit none
c
c arguments
c
      DOUBLE COMPLEX  CHI(2), CEPS(0:3), RESULT(2)
      DOUBLE PRECISION  P(0:3), K(0:4), REPS(0:3), PPLUSK(0:4)
      INTEGER  SIGMAP
      LOGICAL  CHREAL
c
c local variables
c
      DOUBLE PRECISION  A0,A1,A2,A3,B0,B1,B2,B3,C,CR,CI,
     &                  D0,D1,D2,D3,CCR,CCI
      DOUBLE PRECISION  T1R,T1I,T2R,T2I
      DOUBLE PRECISION  PROP
      INTEGER  INDEX1, INDEX2
      LOGICAL  DEXIST
c
c Compute product of the three matrices in A.10 explicitly in terms 
C of components.
c If one component of Chi is real it is denoted by "c", the other is 
C complex and separated as cr + i*ci. If both are complex replace 
C c --> ccr + i*cci. The 2 possible helicity indices are combined by 
C using some Pauli matrix algebra.
c
c Computations are reduced by computing (chi*eps)*(p+k).
c The (chi*eps) temporary is stored in the row vector:
c   ( t1r + I*t1i, t2r + I*t2i )
c
c Since the 4-vector eps can be real or complex, the chi*eps part is 
C done with two different entries into the subroutine: bra2r and bra2c
c
      A0 = REPS(0)
      IF ( SIGMAP.EQ.1 ) THEN
         A1 = REPS(1)
         INDEX1 = 1
         INDEX2 = 2
      ELSE
         A1 = -REPS(1)
         INDEX1 = 2
         INDEX2 = 1
      ENDIF
      A2 = REPS(2)
      A3 = REPS(3)
      DEXIST = .FALSE.
      GOTO 1
c
c entry for complex polarization vector
c
      ENTRY BRA2C( CHI,CHREAL,P,SIGMAP,K,CEPS,RESULT,PPLUSK )
      A0 = DREAL(CEPS(0))
      D0 = DIMAG(CEPS(0))
      IF ( SIGMAP.EQ.1 ) THEN
         A1 = DREAL(CEPS(1))
         D1 = DIMAG(CEPS(1))
         INDEX1 = 1
         INDEX2 = 2
      ELSEIF ( SIGMAP.EQ.-1 ) THEN
         A1 = -DREAL(CEPS(1))
         D1 = -DIMAG(CEPS(1))
         INDEX1 = 2
         INDEX2 = 1
      ELSE
c
c unrecognised value of simgap
c
        WRITE (*,*) "Invalid Sigmap in BRA2 : Sigmap = ",SIGMAP
        RESULT(1) = 0.d0
        RESULT(2) = 0.d0
        RETURN
      ENDIF
      A2 = DREAL(CEPS(2))
      D2 = DIMAG(CEPS(2))
      A3 = DREAL(CEPS(3))
      D3 = DIMAG(CEPS(3))
      DEXIST = .TRUE.

  1   CONTINUE
c
c compute the sum (p+k) and store in pplusk:
c   (p+k)(0:3) = pplusk(0:3), (p+k)**2 = pplusk(4)
c
      PPLUSK(0) = P(0) + K(0)
      PPLUSK(1) = P(1) + K(1)
      PPLUSK(2) = P(2) + K(2)
      PPLUSK(3) = P(3) + K(3)
      PPLUSK(4) = PPLUSK(0)**2 - PPLUSK(1)**2 - 
     &            PPLUSK(2)**2 - PPLUSK(3)**2
      PROP = 1.d0/PPLUSK(4)

      B0 = PPLUSK(0)
      IF ( SIGMAP.EQ.1 ) THEN
            B1 = PPLUSK(1)
      ELSE
            B1 = - PPLUSK(1)
      ENDIF
      B2 = PPLUSK(2)
      B3 = PPLUSK(3)
c
c now calculate Chi*eps*prop
c
      IF ( CHREAL ) THEN
         C  = DREAL(CHI(INDEX1))*PROP
         CR = DREAL(CHI(INDEX2))*PROP
         CI = DIMAG(CHI(INDEX2))*PROP

         IF ( DEXIST ) THEN
            T1R = -C*(A3-A0) + CR*(D2-A1) + CI*(A2+D1)
            T1I = -C*(D3-D0) + (D2-A1)*CI - (A2+D1)*CR
            T2R =  CR*(A3+A0) - C*(A1+D2) - CI*(D3+D0)
            T2I =  C*(A2-D1) + CI*(A3+A0) + CR*(D3+D0)
         ELSE
            T1R = -C*(A3-A0) - CR*A1 + CI*A2
            T1I = -A1*CI - A2*CR
            T2R = CR*(A3+A0) - C*A1
            T2I = C*A2 + CI*(A3+A0)
         ENDIF
      ELSE
         CCR = DREAL(CHI(INDEX1))*PROP
         CCI = DIMAG(CHI(INDEX1))*PROP
         CR  = DREAL(CHI(INDEX2))*PROP
         CI  = DIMAG(CHI(INDEX2))*PROP

         IF ( DEXIST ) THEN
            T1R = -CCR*(A3-A0) + CCI*(D3-D0) + CR*(D2-A1) + CI*(A2+D1)
            T1I = -CCR*(D3-D0) - CCI*(A3-A0) + (D2-A1)*CI - (A2+D1)*CR
            T2R =  CR*(A3+A0) - CCR*(A1+D2) - CCI*(A2-D1) - CI*(D3+D0)
            T2I = CCR*(A2-D1) - CCI*(A1+D2) + CI*(A3+A0) + CR*(D3+D0)
         ELSE
            T1R = -CCR*(A3-A0) - CR*A1 + CI*A2
            T1I = -CCI*(A3-A0) - A1*CI - A2*CR
            T2R = CR*(A3+A0) - CCR*A1 - CCI*A2
            T2I = CCR*A2 - CCI*A1 + CI*(A3+A0)
         ENDIF
      ENDIF

      RESULT(INDEX1) = DCMPLX( T1R*(B0+B3) + T2R*B1 - T2I*B2,
     &                         T1I*(B0+B3) + T2I*B1 + T2R*B2 )
      RESULT(INDEX2) = DCMPLX( T1R*B1 + T1I*B2 + T2R*(B0-B3),
     &                         T1I*B1 - T1R*B2 + T2I*(B0-B3) )
ccc
      END
C
C********************* ket2r, ket2c **********************************
C
C   Calculates the ket vector |k,i> as given in
C      MAD/PH/402, Appendix A, Eqn A.10, bottom eqn.
C      checked on May 26, 88
C
C   Modified to double precision and nonzero time component of 
C   polarization vector Aug. 92
C
C   INPUT:
C            chi      double complex array(2)
C                        ket |i> as given by eq. A.9 if chreal = .true.
C                        any ket |...> if chreal = .false.
C
C            chreal  logical
C                        .true. if one component of chi is real
C                        .false. otherwise
C
C            p      double precision  array(0:3)
C                        standard four momentum for fermion
C                        p(0:3) = [E, px, py, pz]
C
C            sigmap      integer
C                        chirality/helicity factor for fermion
C                        allowed values: +1,-1
C
C            k      double precision  array(0:4)
C                        standard four momentum for boson and mass**2
C                        k(0:4) = [E, kx, ky, kz, mass**2]
C
C            reps   double precision  array(0:3)
C            ceps   double complex array(0:3)
C                   real or complex polarisation vector of the boson. 
C
C   OUTPUT:
C            result      double complex array(1:2)
C                        two component column vector |k,i> on Eqn. A.10
C
C            pmink      double precision  array(0:4)
C                        four momentum p-k. The fourth component
C                        contains the square (p-k)**2
C
      SUBROUTINE KET2R( CHI,CHREAL,P,SIGMAP,K,REPS,RESULT,PMINK )

      implicit none
c
c arguments
c
      DOUBLE COMPLEX  CHI(2), CEPS(0:3), RESULT(2)
      DOUBLE PRECISION  P(0:3), K(0:4), REPS(0:3), PMINK(0:4)
      INTEGER  SIGMAP
      LOGICAL  CHREAL
c
c local variables
c
      DOUBLE PRECISION   A0,A1,A2,A3,B0,B1,B2,B3,C,CR,CI,
     &                   D0,D1,D2,D3,CCI,CCR
      DOUBLE PRECISION   T1R,T1I,T2R,T2I
      DOUBLE PRECISION   PROP
      INTEGER  INDEX1, INDEX2
      LOGICAL  DEXIST
c
c Compute product of the three matrices in A.10 explicitly in terms 
C of components.
c If one component of Chi is real it is denoted by "c", the other is 
C complex and separated as cr + i*ci. If both are complex replace 
C c --> ccr + i*cci. The 2 possible helicity indices are combined by 
C using some Pauli matrix algebra.
c
c Computations are reduced by computing (p-k)*(eps*chi).
c The (eps*chi) temporary is stored in the vector:
c     ( t1r + I*t1i, t2r + I*t2i )
c
c Since the 4-vector eps can be real or complex, the chi*eps part is 
C done with two different entries into the subroutine: bra2r and bra2c
c
      A0 = REPS(0)
      IF ( SIGMAP.EQ.1 ) THEN
         A1 = REPS(1)
         INDEX1 = 1
         INDEX2 = 2
      ELSE
         A1 = -REPS(1)
         INDEX1 = 2
         INDEX2 = 1
      ENDIF
      A2 = REPS(2)
      A3 = REPS(3)
      DEXIST = .FALSE.
      GOTO 1
c
c entry for complex polarization vector
c
      ENTRY KET2C( CHI,CHREAL,P,SIGMAP,K,CEPS,RESULT,PMINK )
      A0 = DREAL(CEPS(0))
      D0 = DIMAG(CEPS(0))
      IF ( SIGMAP.EQ.1 ) THEN
         A1 = DREAL(CEPS(1))
         D1 = DIMAG(CEPS(1))
         INDEX1 = 1
         INDEX2 = 2
      ELSEIF ( SIGMAP.EQ.-1 ) THEN
         A1 = -DREAL(CEPS(1))
         D1 = -DIMAG(CEPS(1))
         INDEX1 = 2
         INDEX2 = 1
      ELSE
c
c unrecognised value of simgap
c
         WRITE(6,*) "Invalid Sigmap in BRA2 : Sigmap = ",SIGMAP
         RESULT(1) = 0.d0
         RESULT(2) = 0.d0
         RETURN
      ENDIF
      A2 = DREAL(CEPS(2))
      D2 = DIMAG(CEPS(2))
      A3 = DREAL(CEPS(3))
      D3 = DIMAG(CEPS(3))
      DEXIST = .TRUE.

  1   CONTINUE
c
c compute the sum (p-k) and store in pmink:
c   (p-k)(0:3) = pmink(0:3), (p-k)**2 = pmink(4)
c
      PMINK(0) = P(0) - K(0)
      PMINK(1) = P(1) - K(1)
      PMINK(2) = P(2) - K(2)
      PMINK(3) = P(3) - K(3)
      
      PMINK(4) = PMINK(0)**2 - PMINK(1)**2 - PMINK(2)**2 - PMINK(3)**2
      PROP = 1.d0/PMINK(4)

      B0 = PMINK(0)
      IF ( SIGMAP.EQ.1 ) THEN
            B1 = PMINK(1)
      ELSE
            B1 = -PMINK(1)
      ENDIF
      B2 = PMINK(2)
      B3 = PMINK(3)
c
c now calculate Chi*eps*prop
c
        IF ( CHREAL ) THEN
           C  = DREAL(CHI(INDEX1))*PROP
           CR = DREAL(CHI(INDEX2))*PROP
           CI = DIMAG(CHI(INDEX2))*PROP

           IF ( DEXIST ) THEN
              T1R = -C*(A3-A0) - CR*(D2+A1) - CI*(A2-D1)
              T1I = -C*(D3-D0) - CI*(D2+A1) + CR*(A2-D1)
              T2R =  C*(D2-A1) + CR*(A3+A0) - CI*(D3+D0)
              T2I = -C*(A2+D1) + CI*(A3+A0) + CR*(D3+D0)
           ELSE
              T1R = -C*(A3-A0) - CR*A1 - CI*A2
              T1I = -CI*A1 + CR*A2
              T2R = -C*A1 + CR*(A3+A0)
              T2I = -C*A2 + CI*(A3+A0)
           ENDIF
        ELSE
           CCR = DREAL(CHI(INDEX1))*PROP
           CCI = DIMAG(CHI(INDEX1))*PROP
           CR  = DREAL(CHI(INDEX2))*PROP
           CI  = DIMAG(CHI(INDEX2))*PROP

           IF ( DEXIST ) THEN
              T1R = -CCR*(A3-A0)+ CCI*(D3-D0) - CR*(D2+A1) - CI*(A2-D1)
              T1I = -CCR*(D3-D0)- CCI*(A3-A0) - CI*(D2+A1) + CR*(A2-D1)
              T2R =  CCR*(D2-A1)+ CCI*(A2+D1) + CR*(A3+A0) - CI*(D3+D0)
              T2I =  CCI*(D2-A1)- CCR*(A2+D1) + CI*(A3+A0) + CR*(D3+D0)
           ELSE
              T1R = -CCR*(A3-A0) - CR*A1 - CI*A2
              T1I = -CCI*(A3-A0) - CI*A1 + CR*A2
              T2R = -CCR*A1 + CCI*A2 + CR*(A3+A0)
              T2I = -CCR*A2 - CCI*A1 + CI*(A3+A0)
           ENDIF
        ENDIF

        RESULT(INDEX1) = DCMPLX( T1R*(B0+B3) + T2R*B1 + T2I*B2,
     &                           T1I*(B0+B3) + T2I*B1 - T2R*B2 )
        RESULT(INDEX2) = DCMPLX( T1R*B1 - T1I*B2 + T2R*(B0-B3),
     &                           T1I*B1 + T1R*B2 + T2I*(B0-B3) )
ccc
      END
C
C**********************  S1R,S1C  *************************************
C  checked on 31.5.88
C
C  S1R and S1C calculate the complex number 
C             <bra|(rvec-slash)_sigma|ket>
C  in the notation of HZ
C
C  INPUT:
C  ------
C
C  BRA(2)           double complex bra vector
C
C  RVEC(0:3)        double precision 4-vector
C  CVEC(0:3)        double complex 4-vector
C
C  TIMEEX           logical variable:
C                      .true. if time component of RVEC/CVEC is nonzero
C                      .false. if time component of RVEC/CVEC vanishes
C
C  SIGMA            integer index of (rvec-slash)_sigma.
C                      Possible values: +1,-1
C
C  KET(2)           complex ket vector
C
C  OUTPUT:
C  -------
C
C  S1R,S1C          complex function value
C
      FUNCTION S1R( BRA,RVEC,TIMEEX,SIGMA,KET )

      implicit none

      DOUBLE COMPLEX  BRA(2), CVEC(0:3), KET(2), S1R, S1C
      DOUBLE PRECISION  RVEC(0:3)
      INTEGER  SIGMA
      LOGICAL  TIMEEX
c
c local variables
c
      DOUBLE PRECISION  A0,A1,A2,A3, B0,B1,B2,B3
      DOUBLE PRECISION  SL11R,SL11I,SL12R,SL12I,SL21R,SL21I, 
     &                  SL22R,SL22I
      DOUBLE PRECISION  K1R,K1I,K2R,K2I, KA1R,KA1I,KA2R,KA2I
cc
      K1R = DREAL(KET(1))
      K1I = DIMAG(KET(1))
      K2R = DREAL(KET(2))
      K2I = DIMAG(KET(2))

      A1 = RVEC(1)
      A2 = RVEC(2)
      A3 = RVEC(3)
      IF ( TIMEEX ) THEN
        A0 = RVEC(0)
        IF ( SIGMA.EQ.1 ) THEN
          SL11R =  A3 - A0
          SL22R = -A3 - A0
        ELSE
          SL11R =  A3 + A0
          SL22R = -A3 + A0
        ENDIF
      ELSE
        SL11R =  A3
        SL22R = -A3
      ENDIF

      KA1R = SL11R*K1R + A1*K2R + A2*K2I
      KA1I = SL11R*K1I - A2*K2R + A1*K2I
      KA2R = SL22R*K2R + A1*K1R - A2*K1I
      KA2I = SL22R*K2I + A1*K1I + A2*K1R

      K1R = DREAL(BRA(1))
      K1I = DIMAG(BRA(1))
      K2R = DREAL(BRA(2))
      K2I = DIMAG(BRA(2))

      IF (SIGMA.EQ.-1) THEN
        S1R = DCMPLX( K1R*KA1R - K1I*KA1I + K2R*KA2R - K2I*KA2I,
     &                K1R*KA1I + K1I*KA1R + K2R*KA2I + K2I*KA2R )
      ELSE
        S1R = DCMPLX( -K1R*KA1R + K1I*KA1I - K2R*KA2R + K2I*KA2I,
     &                -K1R*KA1I - K1I*KA1R - K2R*KA2I - K2I*KA2R )
      ENDIF

      RETURN
c
c entry for complex 4-vector
c
      ENTRY S1C( BRA,CVEC,TIMEEX,SIGMA,KET )

      K1R = DREAL(KET(1))
      K1I = DIMAG(KET(1))
      K2R = DREAL(KET(2))
      K2I = DIMAG(KET(2))

      A1 = DREAL(CVEC(1))
      A2 = DREAL(CVEC(2))
      A3 = DREAL(CVEC(3))

      B1 = DIMAG(CVEC(1))
      B2 = DIMAG(CVEC(2))
      B3 = DIMAG(CVEC(3))

      IF ( TIMEEX ) THEN
        A0 = DREAL(CVEC(0))
        B0 = DIMAG(CVEC(0))
        IF ( SIGMA.EQ.1 ) THEN
          SL11R =  A3 - A0
          SL22R = -A3 - A0
          SL11I =  B3 - B0
          SL22I = -B3 - B0
        ELSE
          SL11R =  A3 + A0
          SL22R = -A3 + A0
          SL11I =  B3 + B0
          SL22I = -B3 + B0
        ENDIF
      ELSE
        SL11R =  A3
        SL22R = -A3
        SL11I =  B3
        SL22I = -B3
      ENDIF

      SL12R = A1 + B2
      SL12I = B1 - A2
      SL21R = A1 - B2
      SL21I = A2 + B1

      KA1R = SL11R*K1R - SL11I*K1I + SL12R*K2R - SL12I*K2I
      KA1I = SL11R*K1I + SL11I*K1R + SL12R*K2I + SL12I*K2R
      KA2R = SL21R*K1R - SL21I*K1I + SL22R*K2R - SL22I*K2I
      KA2I = SL21R*K1I + SL21I*K1R + SL22R*K2I + SL22I*K2R

      K1R = DREAL(BRA(1))
      K1I = DIMAG(BRA(1))
      K2R = DREAL(BRA(2))
      K2I = DIMAG(BRA(2))

      IF ( SIGMA.EQ.-1 ) THEN
        S1R = DCMPLX( K1R*KA1R - K1I*KA1I + K2R*KA2R - K2I*KA2I,
     &                K1R*KA1I + K1I*KA1R + K2R*KA2I + K2I*KA2R )
      ELSE
        S1R = DCMPLX( -K1R*KA1R + K1I*KA1I - K2R*KA2R + K2I*KA2I,
     &                -K1R*KA1I - K1I*KA1R - K2R*KA2I - K2I*KA2R )
      ENDIF
ccc
      RETURN
      END


C
C*********** polarisation vector (POLVEC) ************
C*********** entry HELVEC ****************************
C
C   Calculates the polarisation vector eps_mu(k,lambda) as given in
C   MAD/PH/402, Appendix A, Eqn A.12 and A.13
C   checked on May 24, 88
C   bug in HELVEC for KT=0 corrected Nov.4, 1993
C
C   INPUT:
C
C            kbar      real*8  array(0:4)
C                        physical four momentum for boson and mass
C                        kbar(0:4) = [Ebar, kbarx, kbary, kbarz, mass**2]
C
C            sigmak      integer
C                        sign factor for vector boson
C                        = +1 outgoing
C                        = -1 incoming
C
C            lambda      integer
C                        rectangular polarization index as used
C                        in Eqn A.12. Allowed values= 1,2,3.
C                        or vector boson helicity (allowed values
C                        -1,0,1) when entry HELVEC is used
C
C
C   OUTPUT:
C            epscar  real*8  array(0:3)
C                        polarisation vector epsilon_mu(lambda)
C                               in eq. A.12 (cartesian basis)
C
C            epshel  complex*16  array(0:3)
C                        polarisation vector epsilon_mu(lambda)
C                               in eq. A.12/13 (helicity basis)
C
C                  The polarization vector for helicity -1,1 or
C                  polarization index 1,2 has vanishing time component.
C                  For longitudinal polarisation, current conservation
C                  is used to eliminate the time component. This is
C                  only possible for MASSLESS FERMIONS. For massive
C                  fermions the polarization vectors defined in eq.A.12
C                  must be used!
C
      SUBROUTINE POLVEC(KBAR,LAMBDA,EPSCAR)

      implicit none
C
C   arguments:
C
      double precision  KBAR(0:4)
      double precision  EPSCAR(0:3)
      double complex EPSHEL(0:3)
      INTEGER SIGMAK, LAMBDA
C
C   local variables:
C
      double precision  A1,A2,A3,B1,B2
      double precision  KT,KMAGNI,NORMAL
C
C   code:
C
C
C            Compute polarization vector:
C            Eqn A.12
C
      EPSCAR(0) = 0.0d0
      IF (LAMBDA .EQ. 1) THEN
         KT = SQRT(KBAR(1)**2+KBAR(2)**2)
         IF (KT .GT. 0) THEN
            IF (KBAR(4) .EQ. 0 ) THEN
               KMAGNI = KBAR(0)
            ELSE
               KMAGNI = SQRT(KBAR(1)**2+KBAR(2)**2+KBAR(3)**2)
            ENDIF
            NORMAL = 1.0d0 / (KMAGNI*KT)
            EPSCAR(1) = KBAR(1)*KBAR(3)*NORMAL
            EPSCAR(2) = KBAR(2)*KBAR(3)*NORMAL
            EPSCAR(3) =          -KT**2*NORMAL
         ELSE
            EPSCAR(1) = dsign(1.0d0, KBAR(3))
            EPSCAR(2) = 0.0d0
            EPSCAR(3) = 0.0d0
         ENDIF
      ELSE IF (LAMBDA .EQ. 2) THEN
         KT = SQRT(KBAR(1)**2+KBAR(2)**2)
         IF (KT .GT. 0) THEN
            EPSCAR(1) = -KBAR(2)/KT
            EPSCAR(2) =  KBAR(1)/KT
            EPSCAR(3) = 0.0d0
         ELSE
            EPSCAR(1) = 0.0d0
            EPSCAR(2) = 1.0d0
            EPSCAR(3) = 0.0d0
         ENDIF
      ELSE IF (LAMBDA .EQ. 3) THEN
         IF (KBAR(4) .EQ. 0) THEN
C
C                  Mass = 0 not allowed for longitudinal polarization
C
            WRITE(*,*) "Mass = 0 for Lambda = 3 in POLVEC"
            EPSCAR(1) = 0.0d0
            EPSCAR(2) = 0.0d0
            EPSCAR(3) = 0.0d0
            RETURN
         END IF
            print*,'HELVEC(longitudinal polarisation) is buggy ?!'
!             stop
         KMAGNI = SQRT(KBAR(1)**2+KBAR(2)**2+KBAR(3)**2)
         NORMAL = KBAR(0)/(KBAR(4)*KMAGNI)
         EPSCAR(0) = KMAGNI*KMAGNI/KBAR(0)*NORMAL
         EPSCAR(1) = KBAR(1)*NORMAL
         EPSCAR(2) = KBAR(2)*NORMAL
         EPSCAR(3) = KBAR(3)*NORMAL
      ELSE
C
C                  Unrecognised value of Lambda
C
         WRITE (*,*) "Invalid Lambda in POLVEC: Lambda = ",LAMBDA
         EPSCAR(1) = 0.0d0
         EPSCAR(2) = 0.0d0
         EPSCAR(3) = 0.0d0
      END IF
      RETURN
C
      ENTRY HELVEC(KBAR,SIGMAK,LAMBDA,EPSHEL)
C
C   code:
C
C
C            Compute polarization vector:
C            Use Eqn A.12 and combine according to A.13
C            For outgoing bosons (sigmak=+1) the complex
C            conjugate polarisation vector is returned.
C
      EPSHEL(0) = 0.0d0
      IF (ABS(LAMBDA) .EQ. 1) THEN
          KT = SQRT(2D0 *(KBAR(1)**2+KBAR(2)**2))
            IF (KT .GT. 0) THEN
            IF (KBAR(4) .EQ. 0 ) THEN
                  KMAGNI = KBAR(0)
            ELSE
                  KMAGNI = SQRT(KBAR(1)**2+KBAR(2)**2+KBAR(3)**2)
            ENDIF
            NORMAL = -LAMBDA/(KMAGNI*KT)
            A1 = KBAR(1)*KBAR(3)*NORMAL
            A2 = KBAR(2)*KBAR(3)*NORMAL
            A3 =    -KT**2 * 0.5d0 *NORMAL
            B1 = -KBAR(2)/KT
            B2 =  KBAR(1)/KT
            EPSHEL(1) = dcmplx(A1, SIGMAK * B1)
            EPSHEL(2) = dcmplx(A2, SIGMAK * B2)
            EPSHEL(3) = A3
          ELSE
            NORMAL = -LAMBDA / sqrt( 2D0 )
            EPSHEL(1) = dsign(1.0d0, KBAR(3)) * NORMAL
            EPSHEL(2) = dcmplx(0.0d0, -LAMBDA * SIGMAK * NORMAL)
            EPSHEL(3) = 0.0d0
          ENDIF
      ELSE IF (LAMBDA .EQ. 0) THEN
            IF (KBAR(4) .EQ. 0) THEN
C
C                  Mass = 0 not allowed for longitudinal polarization
C
                  WRITE(*,*) "Mass = 0 for Lambda = 0 in HELVEC"
                  EPSHEL(1) = 0.0d0
                  EPSHEL(2) = 0.0d0
                  EPSHEL(3) = 0.0d0
                  RETURN
                  END IF
            print*,'HELVEC(longitudinal polarisation) is buggy ?!'
!             stop
            KMAGNI = SQRT(KBAR(1)**2+KBAR(2)**2+KBAR(3)**2)
            NORMAL = KBAR(0)/(sqrt(abs(KBAR(4)))*KMAGNI)
            EPSHEL(0) = KMAGNI*KMAGNI/KBAR(0)*NORMAL
            EPSHEL(1) = KBAR(1)*NORMAL
            EPSHEL(2) = KBAR(2)*NORMAL
            EPSHEL(3) = KBAR(3)*NORMAL
            if(KBAR(4).lt.0d0) EPSHEL(0) = EPSHEL(0)*(0d0,1d0)
      ELSE
C
C                  Unrecognised value of Lambda
C
            WRITE (*,*) "Invalid Lambda in HELVEC: Lambda = ",LAMBDA
                  EPSHEL(1) = 0.0d0
                  EPSHEL(2) = 0.0d0
                  EPSHEL(3) = 0.0d0
            END IF
      RETURN
      END

      FUNCTION SC3(CHII,A1,A2,A3,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC3
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A3(0:3), AUX(0:3,3)
      REAL*8  A2(0:3)
C
      N = 3
      DO I = 0,3
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC3  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END

c----------------------------------------------------------------------

      subroutine propagate(chi,curr,mom)

      implicit none

      integer chi ! maximum chirality
      double complex curr(1:6,-1:1) ! current
      double complex currout(1:6,-1:1) ! current
      double precision mom(0:4,-1:1) ! momentum

      integer sig, i
      double precision dotrr
      external dotrr

      do sig = -1,chi,2

        mom(0,sig) = dreal(curr(5,sig))
        mom(1,sig) = dreal(curr(6,sig))
        mom(2,sig) = dimag(curr(6,sig))
        mom(3,sig) = dimag(curr(5,sig))
        mom(4,sig) = dotrr(mom(0,sig),mom(0,sig))
      
        do i=1,4
          curr(i,sig) = -curr(i,sig)/mom(4,sig)   ! -g^\mu_\nu/q^2 (Feynman gauge)
        enddo

      enddo

      end

c----------------------------------------------------------------------

      subroutine BRA2CADD( CHI,CHREAL,P,SIGMAP,K,CEPS,RESULT,PPLUSK,cpl)
      implicit none
      DOUBLE COMPLEX  CHI(2), CEPS(0:3), RESULT(2)
      DOUBLE PRECISION  P(0:3), K(0:4), REPS(0:3), PPLUSK(0:4), cpl
      INTEGER  SIGMAP
      LOGICAL  CHREAL

      double complex tmpres(2)
      double precision tmpp(0:4)
      integer i
      logical lassert
      parameter (lassert=.false.)

      call BRA2C( CHI,CHREAL,P,SIGMAP,K,CEPS,tmpres,tmpp )
      goto 1

      entry BRA2RADD( CHI,CHREAL,P,SIGMAP,K,REPS,RESULT,PPLUSK,cpl)
      call BRA2R( CHI,CHREAL,P,SIGMAP,K,REPS,tmpres,tmpp )

1     continue
      if (lassert) then
        do i=0,4
          if ( dabs(pplusk(i)-tmpp(i)) .gt. 
     &        (1d-8*max(dabs(pplusk(i)),1d0)) ) then
            print*,"Momentum mismatch in BRA2CADD, entry", i
            print*,"expected     value:", pplusk(i)
            print*,"recalculated value:", tmpp(i)
            stop
          endif
        enddo
      endif
      do i=1,2
        result(i) = result(i) + tmpres(i)*cpl
      enddo

      end

c----------------------------------------------------------------------

      subroutine KET2CADD( CHI,CHREAL,P,SIGMAP,K,CEPS,RESULT,PPLUSK,cpl)
      implicit none
      DOUBLE COMPLEX  CHI(2), CEPS(0:3), RESULT(2)
      DOUBLE PRECISION  P(0:3), K(0:4), REPS(0:3), PPLUSK(0:4), cpl
      INTEGER  SIGMAP
      LOGICAL  CHREAL

      double complex tmpres(2)
      double precision tmpp(0:4)
      integer i
      logical lassert
      parameter (lassert=.false.)

      call KET2C( CHI,CHREAL,P,SIGMAP,K,CEPS,tmpres,tmpp )
      goto 1

      entry KET2RADD( CHI,CHREAL,P,SIGMAP,K,REPS,RESULT,PPLUSK,cpl)
      call KET2R( CHI,CHREAL,P,SIGMAP,K,REPS,tmpres,tmpp )

1     continue
      if (lassert) then
        do i=0,4
          if ( dabs(pplusk(i)-tmpp(i)) .gt. 
     &        (1d-8*max(dabs(pplusk(i)),1d0)) ) then
            print*,"Momentum mismatch in KET2CADD, entry", i
            print*,"expected     value:", pplusk(i)
            print*,"recalculated value:", tmpp(i)
            stop
          endif
        enddo
      endif
      do i=1,2
        result(i) = result(i) + tmpres(i)*cpl
      enddo

      end

c----------------------------------------------------------------------

      subroutine curr6addm( sigmax,psibar,pb,psi,p,j,cplm,cplp)
      double complex  j(0:5,-1:1), psibar(2,-1:1), psi(2,-1:1)
      double precision pb(0:3), p(0:3), cplm, cplp

      integer  sigmax

      double complex tmpj(0:5,-1:1)
      double precision cpl(-1:1)
      integer sig, i
      logical lassert
      parameter (lassert=.false.)

      cpl(-1) = cplm
      cpl(+1) = cplp
      goto 1

      entry curr6add( sigmax,psibar,pb,psi,p,j)
      cpl(-1) = 1d0
      cpl(+1) = 1d0

1     continue
      call curr6( sigmax,psibar,pb,psi,p,tmpj )
      if (lassert) then
        do i=4,5
          if ( abs(j(i,-1)-tmpj(i,-1)) .gt. 
     &        (1d-8*max(abs(j(i,-1)),1d0)) ) then
            print*,"Momentum mismatch in curr6add, entry", i
            print*,"expected     value:", j(i,-1)
            print*,"recalculated value:", tmpj(i,-1)
            stop
          endif
        enddo
      endif
      do sig=-1,sigmax,2
        do i=0,3
          j(i,sig) = j(i,sig) + tmpj(i,sig)*cpl(sig)
        enddo
      enddo

      end

c----------------------------------------------------------------------

      subroutine checkorder(i,j,k,l,pos,cpl)

      implicit none

      include 'VBFNLO/utilities/coupl.inc'

      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)

      integer c,i,j,k,l,pos
      integer pos3
      integer pos4
      double precision cpl(5)

! get order of 1 and 2
      if ((i.eq.2).and.((j.eq.1).or.(k.eq.1).or.(l.eq.1)) .or.
     &    (j.eq.2).and.((k.eq.1).or.(l.eq.1)) .or. 
     &    (k.eq.2).and.((l.eq.1)) ) then 
        pos = 2
      else
        pos = 1
      endif
      pos3 = 0
      pos4 = 0
      if (i.eq.3) pos3 = 1
      if (j.eq.3) pos3 = 2
      if (k.eq.3) pos3 = 3
      if (l.eq.3) pos3 = 4
      if (i.eq.4) pos4 = 1
      if (j.eq.4) pos4 = 2
      if (k.eq.4) pos4 = 3
      if (l.eq.4) pos4 = 4

      do c=1,5
        cpl(c) = clr(1,3,-1)
        if ((pos3.gt.0) .and. (pos3.lt.c)) cpl(c) = cpl(c)*clr(4,1,-1)
        if (pos3.ge.c) cpl(c) = cpl(c)*clr(3,1,-1)
        if ((pos4.gt.0) .and. (pos4.lt.c)) cpl(c) = cpl(c)*clr(4,1,-1)
        if (pos4.ge.c) cpl(c) = cpl(c)*clr(3,1,-1)
      enddo
      if (l.eq.0) then
        if (pos3.eq.0) l = 3
        if (pos4.eq.0) l = 4
      endif
      if (k.eq.0) then
        if (pos4.eq.0) k = 4
        if (pos3.eq.0) k = 3
      endif

      end
