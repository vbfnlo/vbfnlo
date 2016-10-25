c======================================================
c
c  Define common block containing KOPPLN inputs.
c
c  Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c  Dave Rainwater,    <rain@pheno.physics.wisc.edu>
c  Initial version:  1997 Mar 3
c
c======================================================
c

      double precision  alfas, xmt, alfa, xmz, xmw, sin2w, xmh, gf, alfas_lo
      common /bkopin/   alfas, xmt, alfa, xmz, xmw, sin2w, xmh, gf, alfas_lo

c  xmb = mb pole mass,   BMASS_H = mb(m_H) etc.
      double precision xmb, BMASS_H, xmc, CMASS_H, xmtau 
      common /quarkmasses/ xmb, BMASS_H, xmc, CMASS_H, xmtau 
      
c name of the SLHA file which is read-in
      character*250 slhaFILE
      common /slha/ slhaFILE
      
