!======================================================
!
!  Define common block containing KOPPLN inputs.
!
!  Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
!  Dave Rainwater,    <rain@pheno.physics.wisc.edu>
!  Initial version:  1997 Mar 3
!
!======================================================
!

      double precision  alfas, xmt, alfa, xmz, xmw, sin2w, xmh, gf, alfas_lo
      common /bkopin/   alfas, xmt, alfa, xmz, xmw, sin2w, xmh, gf, alfas_lo

!  xmb = mb pole mass,   BMASS_H = mb(m_H) etc.
      double precision xmb, BMASS_H, xmc, CMASS_H, xmtau 
      common /quarkmasses/ xmb, BMASS_H, xmc, CMASS_H, xmtau 
      
! name of the SLHA file which is read-in
      character*250 slhaFILE
      common /slha/ slhaFILE
      
