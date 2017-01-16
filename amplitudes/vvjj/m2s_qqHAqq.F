**sophy: 5th January 2010: Altered slightly to account for the PDF arrays
** new file, 27-10-2009
** Calculation of real photon corrections.  By Terrance Figy

c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2003 February 12
c	Last modified:    2003 February 27
c
      subroutine m2s_wbfHA(
     &                   bos,   !in:  Boson identifier, 2 = Z, 6=H
     &                   nlo,   !in:  NLO=1: create subtraction term; LO = 0
     &                   lok,   !in:  need to calculate m2s(1:3)? T or F 
     &                   xi,	!in:  Feynman x parameters of incoming fermions
     &                   p,	!in:  momenta associated with external fermions
     &                   v,	!in:  momenta associated with WW decay fermions
     &                   xuz,   !in:  x, u, z for subtraction terms
     &                   m2s    !out: |M|^2*pdf1*pdf2 
     &                        )
      implicit none
c
c declare input/output variables
c
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
c      include 'hepup.inc'

      integer bos
      real*8 xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), 
     1       rn(6), xuz(2,2:max_kin), m2s(0:max_kin)
      logical lok(max_kin)
c
c declare external functions
c
      integer FL_ZHg
      integer DecayMode
      external FL_ZHg
      external DecayMode


      real*8 RandomNumber
      external RandomNumber


c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"


c
c declare local variables
c
      real*8 q_sf

      integer init/0/, I, J
      save init
c
c declare variables for summation process
c
      INTEGER ires

      integer FSIGN(6), gsign, physToDiag(5), nlo, nmin, nmax

      save nmin, nmax
c
c store contributions from subprocess ip in res(ip,ID) where
c  ID = 1  : the real emission |M|^2 * pdf
c     = 2,3: sutraction terms for emision off 12 or 34 line
c     = 0  : subtracted result which drives integration, i.e
c res(*,0) = res(*,1)-res(*,2)-res(*,3)
      real*8 res(maxnumsubproc,0:max_kin)
c
c declare parton distribution variables
c
      real*8 pdf(-6:6,2)
c
c define program switches
c
      logical ldebug
      parameter (ldebug=.false.)
      logical lInPhoton,lQuarkIni
      parameter (lInPhoton=.false.,lQuarkIni=.true.)

c      data ldebug /.false./		!output debug information
c
c if first iteration, output acceptance cuts and fix constant input 
c parameters for call of qqZWqqji.f
c

** sophy: include file for photon PDF, and new PDF array which includes photons
      double precision Qpdf(-7:7,2)
#include "VBFNLO/utilities/Apdf.inc"


      do i = 1, 6
         rn(i) = RandomNumber()
      end do


      if ( init .eq. 1 ) then
         do i = 1,nmax
            do j=0,max_kin
               res(i,j) = 0
            enddo
         enddo
      else
         if (bos.eq.6) then
            write(6,*) ' '
            write(6,*) 'Hjj + photon amplitude square information:'
            write(6,*) '-----------------------------------------------'
            write(6,*) ' '
            if(lInPhoton) then
               write(6,*) 'initial photons are on'
            else
               write(6,*) 'initial photons are off '
            endif
            if(nlo.eq.1) write(6,*)' this is the subtracted real MEsq'
         endif
         call printnfl(.true.)
         print*,' '
        
         init = 1
         do i = 1,maxnumsubproc
            do j = 0,n_qed
               res(i,j) = 0
            enddo
         enddo

      endif

c define scales
c      call calc_qsq(p,v,q12,q34,3) ! modify calc_qsq 
          
c
c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c RECALL Qpdf(\pm 7) is the photon
c
      q_sf = sqrt(mufsq(1,1))

c      print*,'q_sf=',q_sf
c      read(*,*)
      call pdfproton(xi(1), q_sf, pdf(-6,1) ) 
      Qpdf(-7,1) = GAMpdf
      Qpdf(7,1) = GAMpdf
      call pdfproton(xi(2), q_sf, pdf(-6,2) )  
      Qpdf(-7,2) = GAMpdf
      Qpdf(7,2) = GAMpdf

* sophy: setting up special PDF array to include photons
      do i = 1, 2
         do j = -6, 6
            Qpdf(j,i) = pdf(j,i)
         end do
      end do
         
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Now sum over the subprocesses contributing to Z production

      nmax = FL_ZHg(fsign,-1)            !reset counter for subprocesses to 0
      if(lQuarkIni) then
C*******************  q1 q3 ---> q2 q4 g Z   **********************
         
c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon
C NOTE: for call of wbf_ha2j it is important that p(*,1,*) and p(*,3,*)
c correspond to 1-2 fermion line ALWAYS, i.e physToDiag(1/2)={1,3} and 
c similarly physToDiag(3/4)={2,4} for the 3-4 fermion line
      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1
      gsign    = 1

*      write(*,*)'first, p =', p

*      write(*,*)'pre first wbf_ha2j'

      call wbf_ha2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,
     1     xuz,Qpdf,res,nmin,nmax)

*      write(*,*)'post first wbf_ha2j'


      if (ldebug) call resprint(nmin,nmax,res)

C*******************  q1 qb4 ---> q2 qb3 g Z   **********************

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon

      fsign(3) = -1
      fsign(4) = -1

      call wbf_ha2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,
     1     xuz,Qpdf,res,nmin,nmax)
     
      if (ldebug) call resprint(nmin,nmax,res)

C*******************  qbar2 q3 ---> qbar1 q4 g Z   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = 1
      
      call wbf_ha2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,
     1     xuz,Qpdf,res,nmin,nmax)

      if (ldebug) call resprint(nmin,nmax,res)

C*******************  qbar2 qb4 ---> qbar1 qb3 g Z   ******************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1 
      fsign(4) = -1
         
      call wbf_ha2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,
     1     xuz,Qpdf,res,nmin,nmax)

      if (ldebug) call resprint(nmin,nmax,res)
      endif
      if(lInPhoton) then
C*******************  g q3 ---> qb1 q2 q4 Z   **********************

      physToDiag(1)=5
      physToDiag(2)=3
      physToDiag(3)=1
      physToDiag(4)=4
      physToDiag(5)=2

      fsign(1) = -1
      fsign(2) =  1
      fsign(3) =  1
      fsign(4) =  1
      gsign    = -1

      call wbf_ha2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,
     1     xuz,Qpdf,res,nmin,nmax)

      if (ldebug) call resprint(nmin,nmax,res)

C*******************  q1 g ---> q2 qb3 q4 Z   **********************

      physToDiag(1)=1
      physToDiag(2)=5
      physToDiag(3)=2
      physToDiag(4)=3
      physToDiag(5)=4

      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1

      call wbf_ha2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,
     $     xuz,Qpdf,res,nmin,nmax)

      if (ldebug) call resprint(nmin,nmax,res)

C*******************  g qbar4 ---> qbar3 qb1 q2 Z   **********************

      physToDiag(1)=5
      physToDiag(2)=4
      physToDiag(3)=1
      physToDiag(4)=3
      physToDiag(5)=2

      fsign(1) = -1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) = -1

      call wbf_ha2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,
     1     xuz,Qpdf,res,nmin,nmax)

      if (ldebug) call resprint(nmin,nmax,res)

C*******************  qbar2 g ---> qbar1 qb3 q4 Z   **********************

      physToDiag(1)=2
      physToDiag(2)=5
      physToDiag(3)=1
      physToDiag(4)=3
      physToDiag(5)=4
      
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1

      call wbf_ha2j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,
     $     xuz,Qpdf,res,nmin,nmax)

      if (ldebug) call resprint(nmin,nmax,res)
      endif ! end initial phtn section 
c*****************  end of process evaluation  *********************

      do j=0,n_qed
         m2s(j) = 0
         DO IRES = 1,NMAX
            m2s(j) = m2s(j) + RES(IRES,j)
         ENDDO         
      enddo

c$$$      print*,'m2s(0)=',m2s(0)
c$$$      print*,'lok=',lok(1)
      RETURN
      END


c******************************************************************
c
c   begin subroutine wbf_ha2j
c
c*****************************************************************
      subroutine wbf_ha2j(xi,p,v,physToDiag,fsign,gsign,bos, 
     1                    nlo,lok,xuz,
     2                    Qpdf,res,nmin,nmax)
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/mssm.inc"
#include "VBFNLO/utilities/lha.inc"
      real*8 p(0:3,max_p,max_kin), v(0:3,max_v,max_kin)
      real*8 xi(nx), xuz(2,2:max_kin)
      real*8 Qpdf(-7:7,2), res(maxnumsubproc,0:max_kin)
      integer physToDiag(5), fsign(6), gsign, bos, nlo, nmin, nmax
      logical lok(max_kin)
      integer isoft             !0 for gluon and 7 for photon
      parameter(isoft=7)
c
c wbf_ha2j calls the amplitude square routines 
c             qqzqqj    for qq-->qqZ g   for bos=2      and
c             qqhqq     for qq-->qqH     for bos=6
c for the subtraction terms for the NLO cross section calculation
c
c  INPUT:  p(0:3,5,3)      external physical parton momenta
c          v(0:3,nv)       Z decay momenta, for Higgs production only
c                          the sum q(mu) = v(mu,1)+...+v(mu,nv) is needed
c          physToDiag(5)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign,gsign     sign factors for amplitude calls; see qqZqq(j)
c          nlo             nlo = 0: calculate LO only (i.e. no subtraction)
c                          nlo = 1: calculate full NLO subtraction
c          lok(3)          lok(ID)=T means momenta set ID passes acceptance
c                          cuts, i.e res(k,ID) needs to be calculated
c  OUTPUT:
c          uucc(ID)   output in format as in qqZqq(j), but corresponding
c          etc.            to res(*,ID) in m2s_qqZqq
c                          ID = 1  : the real emission |M|^2 * Qpdf
c                             = 2,3: sutraction terms for emision off 12 or 
c     or   res(k,ID)                 34 line
c                             = 0  : subtracted result which drives 
c                                    integration, i.e 
c                                    res(*,0) = res(*,1)+res(*,2)+res(*,3)
c
c  In and Output
c    nmin, nmax            range of process numbers for this call
c
      real*8 N                              ! color factors
      parameter(N=3d0)

c Note: Factor 9 from color sum included in qqbqq. I am using color summed
c amplitudes here, while CS use color averages. 
c
      real*8 uuccb(-1:1,-1:1,2:max_kin), ddccb(-1:1,-1:1,2:max_kin),
     1       uussb(-1:1,-1:1,2:max_kin), ddssb(-1:1,-1:1,2:max_kin),
     2       udscb(-1:1,-1:1,2:max_kin), ducsb(-1:1,-1:1,2:max_kin),
     3       uucc(0:max_kin), uuss(0:max_kin), ddcc(0:max_kin),
     4       ddss(0:max_kin), udsc(0:max_kin), ducs(0:max_kin)
      real*8 pbar(0:3,6), qbar(0:4),gsub(2:max_kin), 
     1       e_in(2), dotrr, polcolq, polcolg

      double precision NCborn(0:1,0:1,max_kin),tree(6)
      double precision CCborn(0:1,max_kin)
      double precision NCreal(0:1,0:1)
      double precision CCreal(0:1)
      double precision q2,q3,v2,v3
      integer iflav(5), diagToPhys(5), FL_ZHg
      external dotrr, FL_ZHg

      logical ChargedCurrent, sametype, oneAntiparticle
      logical loldp, ldebug,ldebug1,ldebug2
      parameter (ldebug=.false.)
      parameter (ldebug1=.false.,ldebug2=.false.)
      data e_in/0d0,0d0/
      save e_in,loldp, uuccb,uussb,ddccb,ddssb,udscb,ducsb,
     1     polcolq,polcolg
      integer i,if1,if2,icc1,j,k,mu, fc(6), fcb(6),nd
c
      double precision mompart(0:3,1:6),momphtn(0:3),qsq
c
      double precision QuarkECharge(1:5) !d,u,s,c,b relative electric chargese 
      data QuarkECharge /-1.0D0,2.0D0,-1.0D0,2.0D0,-1.0D0/
      double precision e2QiQj, sigmaf(1:4),sum_dipoles,e2,e3
c
c     sigma(1) = fsign(1) = +1 ( incoming quark)
c     sigma(2) = -fsign(2) = -1 (outgoing quark)
c
c     this is the mapping to the dipole subtraction terms
c     this only is relevant for the final state photon radiation
c
      integer ndip 
      integer idip(2:13),jdip(2:13)
      data idip /3,4,3,4,3,4,1,1,2,2,1,2/ ! emitter
      data jdip /4,3,1,1,2,2,3,4,3,4,2,1/ ! spectator


      do j = 0,max_kin
         uucc(j) = 0d0
         uuss(j) = 0d0
         ddcc(j) = 0d0
         ddss(j) = 0d0
         udsc(j) = 0d0
         ducs(j) = 0d0
      enddo
      nmin = nmax+1
      do i = 1,5
         diagToPhys(physToDiag(i)) = i
      enddo


      loldp = e_in(1).eq.p(0,3,1) .and. e_in(2).eq.p(0,4,1)
c  reset the LO amplitude to 0 to avoid wrong subtraction for initial gluon
      if (.not.loldp) then            ! this is a new phase space point
         e_in(1) = p(0,3,1)
         e_in(2) = p(0,4,1)
         do i = -1,1
            do j = -1,1
               do k = 2,max_kin
                  uuccb(i,j,k) = 0d0
                  uussb(i,j,k) = 0d0
                  ddccb(i,j,k) = 0d0
                  ddssb(i,j,k) = 0d0
                  udscb(i,j,k) = 0d0
                  ducsb(i,j,k) = 0d0
               enddo
            enddo
         enddo
c
c  determine factor for 
c  polarization and spin average
         
         polcolq = 1d0/(4d0*N**2*xi(1)*xi(2))
         if(isoft.eq.7) then
            polcolg = 1d0/(4d0*N*xi(1)*xi(2)) ! initial photon
         elseif(isoft.eq.0) then
            polcolg = 1d0/(4d0*N*(N**2-1)*xi(1)*xi(2)) ! initial gluon
         endif
      endif
      
c      print*,polcolg,polcolq
         
c get the real emission amplitude squared, store it in uucs(-1,1) etc.

*      do i = 1, 5
*         write(*,*)'i, physToDiag(i) =', i, physToDiag(i)
*      end do

      do mu = 0,3
*         write(*,*)'mu, p(mu,5,1) =', mu, p(mu,5,1)
         do i = 1,5
            pbar(mu,physToDiag(i))=p(mu,i,1)
         enddo
         qbar(mu) = pbar(mu,5)
*         write(*,*)'pbar(mu,5) =', pbar(mu,5)
      enddo
      qbar(4)=0

      if (bos.eq.6) then
         do mu = 0,3
            pbar(mu,5) = 0      ! dummy momentum
            pbar(mu,6) = 0
            do i = 1,n_v
               pbar(mu,6) = pbar(mu,6) + v(mu,i,1) ! Higgs momentum
            enddo
         enddo
         fsign(5) = 0
         fsign(6) = 1

*         write(*,*)'m2s_qqHAqq, qbar =', qbar
*         write(*,*)'m2s_qqHAqq, pbar =', pbar

         if(isoft.eq.7) then    ! for photons
            if (lok(1) .or. nlo.eq.0) then
               call qqHqqPHT(pbar,fsign,qbar,gsign,
     $           uucc(1),uuss(1),ddcc(1),ddss(1),udsc(1),ducs(1))
c$$$               print*,'-------'
c$$$               print*,'udsc=',udsc(1)
c$$$               print*,'uucc=',uucc(1)
c$$$               
c$$$               do i=1,6
c$$$                  tree(i)=0d0
c$$$               enddo
c$$$
c$$$               call qqHAqq_mg(pbar,fsign,qbar,gsign,
c$$$     $              tree(1),tree(2),tree(3),tree(4),tree(5),tree(6))
c$$$               print*,'tree(1)=',tree(1)
c$$$               print*,'ratio=',tree(1)/uucc(1)
c$$$               print*,'tree(5)=',tree(5)
c$$$               print*,'tree(6)=',tree(6)
c$$$               print*,'ratio=',udsc(1)/tree(5)
c$$$               print*,'fsign:',fsign
c$$$               print*,'gsign:',gsign
c$$$               print*,'---------'
            endif
         endif
c$$$         read(*,*)
      else
         print*,' Invalid entry BOS = ',bos
         print*,' in call of wbf_ha2j. Must be 2 or 6. '
         stop
      endif
c
c the real matrix elements squareds
c
      NCreal(0,0)=uucc(1)
      NCreal(0,1)=uuss(1)
      NCreal(1,0)=ddcc(1)
      NCreal(1,1)=ddss(1)
      CCreal(0)=udsc(1)
      CCreal(1)=ducs(1)

c      print*,'ncreal=',ncreal
c      print*,'ccreal=',ccreal

c for the NLO case get the subtraction terms; 
c first the case with a final state gluon
      if (gsign.eq.1) then
c$$$         if (ldebug) then
c$$$            print*,' final state photon section in wbf_ha2j '
c$$$            print 101,' fsign = ',fsign
c$$$ 101        format(a,6i5,a,i5)
c$$$         endif
         if (nlo.eq.1) then
            do j = 2,n_qed            ! loop of dipoles
               do mu = 0,3           ! 
                  do i = 1,4
                     pbar(mu,physToDiag(i))=p(mu,i,j)
                  enddo
               enddo
               
               do mu = 0,3
                  pbar(mu,5) = 0 ! dummy momentum
                  pbar(mu,6) = 0
                  do i = 1,max_v
                     pbar(mu,6) = pbar(mu,6) + v(mu,i,j) ! Higgs momentum
                  enddo
               enddo
               fsign(5) = 0
               fsign(6) = 1

               if (bos.eq.6 .and. lok(j)) then
                  call qqhqqi( 0 )
                  call qqHqq(pbar,fsign,0,
     1                       uucc(j),uuss(j),ddcc(j),
     2                       ddss(j),udsc(j),ducs(j),tree)
c     compute dipole factors
c     The photon is p5
c
                  do mu = 0,3
                     mompart(mu,1) = p(mu,1,1) ! quark
                     mompart(mu,2) = p(mu,2,1) ! quark
                     mompart(mu,3) = p(mu,3,1) ! quark
                     mompart(mu,4) = p(mu,4,1) ! quark
                     mompart(mu,5) = v(mu,1,1) ! lepton/b quark
                     mompart(mu,6) = v(mu,2,1) ! lepton/b quark
                     momphtn(mu)   = p(mu,5,1) ! photon
                  enddo
                  call get_qed_dipole(mompart,momphtn,idip(j),
     $                 jdip(j),xuz(1,j),gsub(j))  
               else
                  uucc(j) = 0d0
                  uuss(j) = 0d0
                  ddcc(j) = 0d0
                  ddss(j) = 0d0
                  udsc(j) = 0d0
                  ducs(j) = 0d0

                  gsub(j) = 0d0
c
               endif

               NCborn(0,0,j)=uucc(j)
               NCborn(0,1,j)=uuss(j)
               NCborn(1,0,j)=ddcc(j)
               NCborn(1,1,j)=ddss(j)
               CCborn(0,j)=udsc(j)
               CCborn(1,j)=ducs(j)
c
c We only need two dipole factors for the initial 
c save matrix elements for later use with initial photons
c
               uuccb(fsign(1),fsign(3),j) = uucc(j)
               uussb(fsign(1),fsign(3),j) = uuss(j)
               ddccb(fsign(1),fsign(3),j) = ddcc(j)
               ddssb(fsign(1),fsign(3),j) = ddss(j)
               udscb(fsign(1),fsign(3),j) = udsc(j)
               ducsb(fsign(1),fsign(3),j) = ducs(j)

c              
            enddo
         endif
c
         iflav(5) = 21           ! final state photon id
         do if1=1,nflVBF           !(nfl/2)*2
            iflav(1)=if1*fsign(physToDiag(1))
            iflav(3)=if1*fsign(physToDiag(3))
            sigmaf(1) = fsign(physToDiag(1)) ! charge flow sign
            sigmaf(3) = -fsign(physToDiag(3)) ! charge flow sign
            do if2=1,nflVBF       !(nfl/2)*2
               iflav(2)=if2*fsign(physToDiag(2))
               iflav(4)=if2*fsign(physToDiag(4))
               sigmaf(2) = fsign(physToDiag(2)) ! charge flow sign
               sigmaf(4) = -fsign(physToDiag(4)) ! charge flow sign
               k=FL_ZHg(iflav,2)
               if (lok(1)) then
                  res(k,1)=Qpdf(sign1*iflav(1),1)
     &                 *Qpdf(sign2*iflav(2),2)
     &                 *NCreal(mod(if1,2),mod(if2,2))*polcolq
               else
                  res(k,1) = 0
               endif

               if(nlo.eq.1) then
                  res(k,0) = res(k,1) ! initialize to real
                  do j=2,n_qed
                     if(lok(j)) then  ! check to see if event passes cuts
c                        print*,'---------------'
c get dipole factor (done above)                       
                        res(k,j) = Qpdf(sign1*iflav(1),1)
     &                 *Qpdf(sign2*iflav(2),2)*
     &                       NCborn(mod(if1,2),mod(if2,2),j)*polcolq

                        e2QiQj = QuarkECharge(abs(iflav(jdip(j))))*
     $                       QuarkECharge(abs(iflav(idip(j))))/9.0D0*
     $                       dble(sigmaf(idip(j))*sigmaf(jdip(j)))*EL**2
c
c                        print*,'e2QiQj=',e2QiQj
c
                        res(k,j) = res(k,j) * gsub(j) * e2QiQj
                        
c compute subtraction term                        
                     else
                        res(k,j) = 0
                     endif

                     
                     res(k,0) = res(k,0) + res(k,j)
c                     print*,'-------------------------'
                  enddo
               else
                  res(k,0) = res(k,1)
                  do nd = 2,n_qed
                     res(k,nd) = 0
                  enddo
               endif ! 
c check subtraction here   
c     
               if(ldebug) then
c               qsq = dotrr(p(0,3,1),p(0,5,1)) ! check initial state singularities
                  qsq = p(0,5,1) ! the soft singularity 
               if(qsq.lt.1.0d-2) then
c
                  print*,'------------'
                  print*,'NC processes'
                  print*,' final state photon section in wbf_ha2j '
                  print 101,' fsign = ',fsign
 101              format(a,6i5,a,i5)
c     
                  sum_dipoles = 0
                  do j=2,n_qed
                     sum_dipoles = sum_dipoles + res(k,j)
                     
                     if (res(k,1).gt.0.0d0) then
                        print*,'born/real (',idip(j),',',jdip(j),')=',
     $                       res(k,j)/res(k,1)                       
                     else
                        print*,'born(',idip(j),',',jdip(j),')=',res(k,j)
c                        print*,'real(',idip(j),',',jdip(j),')=',res(k,1) 
                     endif    
                  enddo
c
                  if (res(k,1).gt.0.0d0) then                     
                     print*,'sum_dipoles/real=',
     $                    sum_dipoles/res(k,1)
                  else
                     print*,'real=',res(k,1)
                     print*,'sum_dipoles=',sum_dipoles
                  endif
                  
               endif            !debug
               endif
                 
               
C Now check if there is a CC contribution for this choice of initial state
C flavors (i,j). First: check whether initial uu or dd, not ud
               sametype=(mod(if1,2)).eq.(mod(if2,2))
               oneAntiparticle=
     &            ((fsign(physToDiag(1))*fsign(physToDiag(2))).eq.-1) 
c true if only one particle is an antiparticle
               ChargedCurrent=(oneAntiparticle.and.sametype) .or.
     &                 ( (.not.oneAntiparticle).and.(.not.sametype) )
               if ( ChargedCurrent .and. if1.le.4 .and. if2.le.4 ) then  
! never allow external bottom quarks in charged currents
c change 1<-->2 and 3<-->4 in outgoing quark flavors
                  iflav(3)=(if1+2*mod(if1,2)-1)*fsign(physToDiag(3))
                  iflav(4)=(if2+2*mod(if2,2)-1)*fsign(physToDiag(4))
                  icc1 = abs(iflav(diagtophys(1)))

                  k=FL_ZHg(iflav,2)

                  if (lok(1)) then
                     res(k,1)=Qpdf(sign1*iflav(1),1)
     &                    *Qpdf(sign2*iflav(2),2)
     &                    *CCreal(mod(icc1,2))*polcolq
                  else
                     res(k,1) = 0
                  endif
                  
                  if(nlo.eq.1) then
                     res(k,0) = res(k,1) ! initialize to real
                     do j=2,n_qed
                        if(lok(j)) then ! check to see if event passes cuts 
                           res(k,j) = Qpdf(sign1*iflav(1),1)
     &                 *Qpdf(sign2*iflav(2),2)*
     &                       CCborn(mod(icc1,2),j)*polcolq

                        e2QiQj = QuarkECharge(abs(iflav(jdip(j))))*
     $                       QuarkECharge(abs(iflav(idip(j))))/9.0d0*
     $                       dble(sigmaf(idip(j))*sigmaf(jdip(j)))*EL**2  
c     compute subtraction term                        
                        res(k,j) = res(k,j) * gsub(j) * e2QiQj
c
                        else
                           res(k,j) = 0
                        endif
                        res(k,0) = res(k,0) + res(k,j)
                     enddo
                  else
                     res(k,0) = res(k,1)
                     do nd = 2,n_qed
                        res(k,nd) = 0
                     enddo
                  endif         ! 
c     you will put the subtraction stuff here  
c                  
                  if(ldebug) then
                     if(qsq.lt.1.0d-2) then
c     
                        print*,'------------'
                        print*,'CC processes'
c     
                        print*,'NC processes'
                        print*,' final state photon section in wbf_ha2j'
                        print 101,' fsign = ',fsign
                        sum_dipoles = 0
                        do j=2,n_qed
                           sum_dipoles = sum_dipoles + res(k,j)
                           
                           if (res(k,1).gt.0.0d0) then
                              print*,'born/real (',idip(j),',',jdip(j),
     &                             ')=',res(k,j)/res(k,1)                       
                           else
                              print*,'born(',idip(j),',',jdip(j),')=',res(k,j)
c     print*,'real(',idip(j),',',jdip(j),')=',res(k,1)                        
                           endif    
                        enddo
c     
                        if (res(k,1).gt.0.0d0) then                     
                           print*,'sum_dipoles/real=',
     $                          sum_dipoles/res(k,1)
                        else
                           print*,'real=',res(k,1)
                           print*,'sum_dipoles=',sum_dipoles
                        endif
                        
                     endif      !debug
                  endif
c                 
                  iflav(3)=if1*fsign(physToDiag(3))
               endif            ! endif(CC)
            enddo
         enddo 

      elseif (gsign.eq.-1) then           !initial photon section
         j = 0
         if (nlo.eq.0) then
            if (fsign(1).eq.-fsign(2)) then 
               j=2              ! j=2: emission off 1,2 line
            elseif(fsign(3).eq.-fsign(4)) then 
               j=3              ! j=3: emission off 3,4 line
            endif
         elseif(nlo.eq.1) then
            do i=1,6
               fc(i)=fsign(i)
               fcb(i)=fsign(i)
            enddo
        
            if (fsign(1).eq.-fsign(2)) then
               j=2              ! j=2: emission off 1,2 line
               ndip = 8          !dipole factor does this work?
               fc(2)=-fsign(2)
               fcb(1)=-fsign(1)

               q2 = dotrr(qbar,pbar(0,2))*xuz(1,ndip)
               q3 = dotrr(qbar,pbar(0,1))*xuz(1,ndip)
            elseif(fsign(3).eq.-fsign(4)) then
               j=3              ! j=3: emission off 3,4 line
               ndip = 11          !dipole factor
               fc(4)=-fsign(4)
               fcb(3)=-fsign(3)
               q2 = dotrr(qbar,pbar(0,4))*xuz(1,ndip)
               q3 = dotrr(qbar,pbar(0,3))*xuz(1,ndip)
            endif

            v2 = EL**2/9.0D0*( (1.0d0 - xuz(1,ndip))**2 + 
     $           xuz(1,ndip)**2 )
            v3 = v2

c
c     get born matrix elements    
c
c     photon + Q -> q qbar Q 
c
            NCborn(0,0,2)=uuccb(fc(1),fc(3),ndip)
            NCborn(0,1,2)=uussb(fc(1),fc(3),ndip)
            NCborn(1,0,2)=ddccb(fc(1),fc(3),ndip)
            NCborn(1,1,2)=ddssb(fc(1),fc(3),ndip)
            CCborn(0,2)=udscb(fc(1),fc(3),ndip)
            CCborn(1,2)=ducsb(fc(1),fc(3),ndip)

            NCborn(0,0,3)=uuccb(fcb(1),fcb(3),ndip)
            NCborn(0,1,3)=uussb(fcb(1),fcb(3),ndip)
            NCborn(1,0,3)=ddccb(fcb(1),fcb(3),ndip)
            NCborn(1,1,3)=ddssb(fcb(1),fcb(3),ndip)
            CCborn(0,3)=udscb(fcb(1),fcb(3),ndip)
            CCborn(1,3)=ducsb(fcb(1),fcb(3),ndip)

c            print*,'udscb',udscb(1,1,7),udscb(1,1,10)
c
c     dipole factor
c

         endif
c
c         print*,'j=',j

         do if1=1,nflVBF            !(nfl/2)*2
            do if2=1,nflVBF         !(nfl/2)*2
               iflav(j-1) = isoft !0
               iflav(4-j)=if1*fsign(physToDiag(4-j))
               iflav(6-j)=if1*fsign(physToDiag(6-j))
               iflav(j+1)=if2*fsign(physToDiag(j+1))
               iflav(5)  =if2*fsign(physToDiag(5))
               k=FL_ZHg(iflav,j)
               if (lok(1)) then
                  res(k,1)=Qpdf(sign1*iflav(1),1)    
     &                    *Qpdf(sign2*iflav(2),2)     
     &                 *NCreal(mod(if1,2),mod(if2,2))*polcolg
               else
                  res(k,1) = 0
               endif
c
c               print*,'ndip=',ndip
c
               if(nlo.eq.1) then
                  res(k,0) = res(k,1) ! initialize to real
c
                  if(j.eq.3) e2QiQj = QuarkECharge(if2)**2 
                  if(j.eq.2) e2QiQj = QuarkECharge(if1)**2 
c
                  do nd=2,n_qed
                     if(lok(nd) .and. nd.eq.ndip) then ! check to see if event passes cuts                   
c     compute subtraction term 
                        res(k,nd) = 
     $                       Qpdf(sign1*iflav(1),1)    
     &                       *Qpdf(sign2*iflav(2),2)     
     &                       *(NCborn(mod(if1,2),mod(if2,2),2)*v2/q2+
     $                       NCborn(mod(if1,2),mod(if2,2),3)*v3/q3)*
     $                       polcolg * e2QiQj !QuarkECharge(if2)**2 
c                       
c                        print*,'e2Q2=',QuarkECharge(if2)
c     
                     else
                        res(k,nd) = 0
                  endif
                  res(k,0) = res(k,0) - res(k,nd)
               enddo
            else
               res(k,0) = res(k,1)
               do nd = 2,n_qed
                  res(k,nd) = 0
               enddo
            endif
c     res(k,0) = res(k,1)
c     
            if(ldebug1.and.lok(1)) then
c               qsq = dotrr(qbar,pbar(0,4))
c               qsq = dotrr(qbar,pbar(0,4))
c               qsq = dotrr(qbar,pbar(0,1))
c               qsq = dotrr(qbar,pbar(0,3))
c               print*,'here'
c     qsq = p(0,5,1) ! the soft singularity 
               qsq = 1d-2
               if((q3.lt.qsq .or. q2.lt.qsq) .and. res(k,1).gt.1.0d-2) then
c     
                  print*,'------------'
                  print*,'NC processes'
                  print*,' final state photon section in wbf_ha2j '
                  print 101,' fsign = ',fsign
                  print 101,' iflav = ',iflav
                  print*,'charge=',e2QiQj
c     101              format(a,6i5,a,i5)
c     
                  sum_dipoles = 0
                  do nd=2,n_qed
                     sum_dipoles = sum_dipoles + res(k,nd)
                     
                     if (res(k,1).gt.0.0d0) then
                        print*,'born/real (',idip(nd),',',jdip(nd),')=',
     $                       res(k,nd)/res(k,1)                       
                     else
                        print*,'born(',idip(nd),',',jdip(nd),')=',res(k,nd)
c     print*,'real(',idip(j),',',jdip(j),')=',res(k,1)                        
                     endif    
                  enddo
c     
                  if (res(k,1).gt.0.0d0) then                     
                     print*,'sum_dipoles/real=',
     $                    sum_dipoles/res(k,1)
                     print*,'real=',res(k,1)
                  else
                     print*,'real=',res(k,1)
                     print*,'sum_dipoles=',sum_dipoles
                  endif
                  
               endif            !debug
            endif
 
c     perform subtraction here

         enddo
      enddo
C Next the CC contributions; no b quarks allowed here
         iflav(diagToPhys(5)) = isoft !0
         do if1 = 1,(nfl/2)*2    !4    sum over all flavors for quark 1 and 2
            iflav(diagToPhys(1)) = fsign(1)*if1
            iflav(diagToPhys(2)) = fsign(2)*(if1-1+2*mod(if1,2))
            do if2 = mod(if1,2)+1,(nfl/2)*2,2 ! flavor of q3 set by q2 (mod 2)
               iflav(diagToPhys(3)) = fsign(3)*if2
               iflav(diagToPhys(4)) = fsign(4)*(if2-1+2*mod(if2,2))
               k=FL_ZHg(iflav,j)
               if (lok(1)) then
                  res(k,1)=Qpdf(sign1*iflav(1),1)
     &                    *Qpdf(sign2*iflav(2),2)
     &                    *CCreal(mod(if1,2))*polcolg
               else
                  res(k,1) = 0
               endif

               if(nlo.eq.1) then
                  res(k,0) = res(k,1) ! initialize to real
                  
                  do nd=2,n_qed
                     if(lok(nd) .and. nd.eq.ndip) then ! check to see if event passes cuts                   

                        if(j.eq.3) then
                           e2 = QuarkECharge(
     $                          abs(iflav(diagToPhys(1))))**2
                           e3 = QuarkECharge(
     $                          abs(iflav(diagToPhys(2))))**2
                        elseif(j.eq.2) then
                           e2 = QuarkECharge(
     $                          abs(iflav(diagToPhys(3))))**2
                           e3 = QuarkECharge(
     $                          abs(iflav(diagToPhys(4))))**2
                        endif

c     compute subtraction term 
                        res(k,nd) = Qpdf(sign1*iflav(1),1)    
     &                       *Qpdf(sign2*iflav(2),2)     
     &                       *(CCborn(mod(if1,2),2)*v2*e2/q2+
     $                       CCborn(mod(if1,2),3)*v3*e3/q3)*polcolg
     $                       
                     else
                        res(k,nd) = 0
                     endif
                     res(k,0) = res(k,0) - res(k,nd)
                  enddo
               else
                  res(k,0) = res(k,1)
                  do nd = 2,n_qed
                     res(k,nd) = 0
                  enddo
               endif
               
               
c     perform subtraction here
c     res(k,0) = res(k,1)
               
               if(ldebug2.and.lok(1)) then
                  qsq = 0.1d0*xuz(1,ndip)
c                  if(xuz(1,ndip).lt.qsq) then
                  if((q3.lt.qsq .or. q2.lt.qsq) .and. res(k,1).gt.1.0d-2) then
c     
                     print*,'------------'
                     print*,'CC processes'
                     print*,' final state photon section in wbf_ha2j '
                     print 101,' fsign = ',fsign
                     print 101,' iflav = ',iflav
                     print*,'e2,e3:',e2,e3,j

c     
                     sum_dipoles = 0
                     do nd=2,n_qed
                        sum_dipoles = sum_dipoles + res(k,nd)
                        
                        if (res(k,1).gt.0.0d0) then
                           print*,'born/real (',idip(nd),',',jdip(nd),
     $                          ')=',res(k,nd)/res(k,1)                       
                        else
                           print*,'born(',idip(nd),',',jdip(nd),')=',res(k,nd)
c     print*,'real(',idip(j),',',jdip(j),')=',res(k,1)                        
                        endif    
                     enddo
c     
                     if (res(k,1).gt.0.0d0) then                     
                        print*,'sum_dipoles/real=',
     $                       sum_dipoles/res(k,1)
                        print*,'real=',res(k,1)
                     else
                        print*,'real=',res(k,1)
                        print*,'sum_dipoles=',sum_dipoles
                     endif
                     
                  endif         !debug
               endif            !debug
               
            enddo
         enddo 
      endif                     ! gsign 
      
      nmax = FL_ZHg(iflav,-2)

      end
c******************************************************************
c
c   end subroutine wbf_ha2j
c
c*****************************************************************
