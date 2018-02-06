********************************************************************************
********************************************************************************
***                                                                          ***
***  ew_qqV_eff.F                                                            ***
***  sophy@particle.uni-karlsruhe.de                                         ***
***                                                                          ***
***  This set of subroutines calculates an effective qqV coupling based on   ***
***  the loop and counterterm diagrams.                                      ***
***                                                                          ***
********************************************************************************
********************************************************************************

*** The 'initialisation' and 'finalisation' subroutine: setting the masses, 
*** charges, couplings, momenta in the correct order and then combining to
*** give us a coupling in VBFNLO format.

      subroutine qqV_eff(pbar,fsign,fcplCT,debug)

      implicit none

** Debugging flag
      integer debug


#include "mssm.inc"
#include "global.inc"


** The effective qqV coupling and its tree level component
** Labelling: wbf(cur, qk, gen)
** cur       =   1    u + c -> u + H + c
**               2    u + s -> u + H + s
**               3    d + c -> d + H + c
**               4    d + s -> d + H + s
**               5    u + s -> d + H + c
**               6    d + c -> u + H + s
** qk        =   1    q + q        ->   q + H + q
**               2    qbar + q     ->   qbar + H + q
**               3    q + qbar     ->   q + H + qbar
**               4    qbar + qbar  ->   qbar + H + qbar
** NOTE: Be careful of charged current cur-qk labels
** gen       =   1    upper line, first generation
**               2    lower line, second generation    
      double complex wbfL(6,4,2), wbfR(6,4,2)
      double complex wbfLtr(6,4,2), wbfRtr(6,4,2)


** Coupling for use in VBFNLO.  This equals:
**        upperCT*lowerTREE + upperTREE*lowerCT
** Labelling: (isig, cur)
** isig = 1     LH * LH
**        2     LH * RH
**        3     RH * LH
**        4     RH * RH
      double complex fcplCT(4,6)


** pbar and fsign are the inputs from the main vbfnlo program - the momentum
** and the particle/antiparticle label.  The final entry in each array is
** particle number i.  p is the physical momentum and mom2 is the square of the
** momenta.  For pbar, i does not necessarily correspond to the particle label 
** i - see notebook 7, page 12.
** Labelling of p: (spacetime index,particle number i)
** Labelling of mom2: (particle number i)
**                 i = 1         upper incoming quark
**                   = 2         upper outgoing quark
**                   = 3         lower incoming quark
**                   = 4         lower outgoing quark
**                   = 5         upper vector boson
**                   = 6         lower vector boson
      double precision pbar(0:3,4+max_v), p(0:3,6), mom2(6)
      integer fsign(4+max_v)

** Do loop parameters
      integer qk, mu, i, cur


** Clearing the values of fcplCT
      do cur = 1, 6
         do qk = 1, 4
            fcplCT(qk,cur) = (0D0,0D0)
         end do
      end do
      
      
** Inputs - finding the quark pattern from fsign.  See notebook 7, page 1.
      if (fsign(1) .gt. 0) then
         if (fsign(3) .gt. 0) then
            qk = 1
         else
            qk = 3
         end if 
      else
         if (fsign(3) .gt. 0) then
            qk = 2
         else
            qk = 4
         end if
      end if

      if (debug .ne. 0) then
         write(*,*)'  '
         write(*,*)'in qqV, qk =', qk
      end if


** Setting the physical momentum.  I realise the momentum notation is
** confusing - see Notebook 7, page 12.
      do i = 1, 4
         do mu = 0, 3
            p(mu,i) = pbar(mu,i)
         end do
      end do

      if (fsign(1) .eq. -1) then
         do mu = 0,3
            p(mu,1) = pbar(mu,2)
            p(mu,2) = pbar(mu,1)
         end do
      end if      
      if (fsign(3) .eq. -1) then
         do mu = 0,3
            p(mu,3) = pbar(mu,4)
            p(mu,4) = pbar(mu,3)
         end do
      end if

** vector bosons' momenta
      do mu = 0, 3
         p(mu,5) = p(mu,1) - p(mu,2)
         p(mu,6) = p(mu,3) - p(mu,4)
      end do

** momentum squared calculation
      do i = 1, 6
         mom2(i) = p(0,i)**2 - p(1,i)**2 - p(2,i)**2 - p(3,i)**2
*         if (i .eq. 1) then 
*            write(*,*)'  '
*            write(*,*)'qk =', qk
*         end if
*         write(*,*)'i, mom2(i) =', i, mom2(i)
*         do mu = 0, 3
*            write(*,*)'mu, p(mu,i) =', mu, p(mu,i)
*         end do
      end do


** Setting the mass and charge for each qqV current when all externals are 
** quarks (rather than antiquarks).  
** labelling: qqV current, particle number, generation
** qqV current   1   u -> u (NC),  u -> d (CC)
**               2   d -> d (NC),  d -> u (CC)
** particle number   1 -> 2
** generation  1 = upper line, 2 = lower line
      qqVmass(1,1,1) = MsU
      qqVmass(1,2,1) = MD
      qqVmass(2,1,1) = MD
      qqVmass(2,2,1) = MsU

      qqVmass(1,1,2) = MC
      qqVmass(1,2,2) = MS
      qqVmass(2,1,2) = MS
      qqVmass(2,2,2) = MC

** labelling: qqV current, particle number
      qqVcharge(1,1) = Qu
      qqVcharge(1,2) = Qd
      qqVcharge(2,1) = Qd
      qqVcharge(2,2) = Qu


      if (qk .eq. 1) then

         if (debug .ne. 0) then
            write(*,*)'  '
            write(*,*)'Calculating effective qqV couplings'
         end if

** Subroutine that calculates the effective coupling values
         call qqV(mom2,model,sector,wbfLtr,wbfRtr,wbfL,wbfR,debug)


** Writing out the couplings ... 
         if (debug .ne. 0) then
            do cur = 1, 6
               if (cur .eq. 1) then
                  write(*,*) '             '
                  write(*,*)'And the results ...'
               end if
               write(*,*)'CURRENT =', cur
** Here, i runs over qk pattern, NOT particle number
               do i = 1, 4
                  write(*,*)'qk =', i
                  write(*,*)'upper cplg (L/R) = ', wbfL(cur,i,1),
     -                 wbfR(cur,i,1)
                  write(*,*)'lower cplg (L/R) = ', wbfL(cur,i,2),
     -                 wbfR(cur,i,2)
                  write(*,*)'  '
               end do
               write(*,*)'  '
            end do            
         end if
      end if


** Calculating fcplCT.  This equals:
**        upperCT*lowerTREE + upperTREE*lowerCT
** Labelling: fcplCT(isig, wbf current)
** isig = 1     LH * LH
**        2     LH * RH
**        3     RH * LH
**        4     RH * RH
** Labelling: wbf(cur,qk,gen)
      do cur = 1, 6

** upper + lower
         fcplCT(1,cur) = wbfL(cur,qk,1)*wbfLtr(cur,qk,2) + 
     -        wbfLtr(cur,qk,1)*wbfL(cur,qk,2)
         if (cur .le. 4) then
            fcplCT(2,cur) = wbfL(cur,qk,1)*wbfRtr(cur,qk,2) + 
     -           wbfLtr(cur,qk,1)*wbfR(cur,qk,2)
            fcplCT(3,cur) = wbfR(cur,qk,1)*wbfLtr(cur,qk,2) + 
     -           wbfRtr(cur,qk,1)*wbfL(cur,qk,2)
            fcplCT(4,cur) = wbfR(cur,qk,1)*wbfRtr(cur,qk,2) + 
     -           wbfRtr(cur,qk,1)*wbfR(cur,qk,2)
         end if

*** For the (very important) purposes of debugging, can include loop 
*** corrections for only upper or lower vertices if desired
** upper
*         fcplCT(1,cur) = wbfL(cur,qk,1)*wbfLtr(cur,qk,2) 
*         if (cur .le. 4) then
*            fcplCT(2,cur) = wbfL(cur,qk,1)*wbfRtr(cur,qk,2) 
*            fcplCT(3,cur) = wbfR(cur,qk,1)*wbfLtr(cur,qk,2) 
*            fcplCT(4,cur) = wbfR(cur,qk,1)*wbfRtr(cur,qk,2) 
*         end if

** lower
*         fcplCT(1,cur) = wbfLtr(cur,qk,1)*wbfL(cur,qk,2)
*         if (cur .le. 4) then
*            fcplCT(2,cur) = wbfLtr(cur,qk,1)*wbfR(cur,qk,2)
*            fcplCT(3,cur) = wbfRtr(cur,qk,1)*wbfL(cur,qk,2)
*            fcplCT(4,cur) = wbfRtr(cur,qk,1)*wbfR(cur,qk,2)
*         end if


         if (qk .eq. 2 .or. qk .eq. 3) then
            do i = 1, 4
               fcplCT(i,cur) = -fcplCT(i,cur)
            end do
         end if

         if (debug .ne. 0) then
            write(*,*)'CURRENT =', cur
            write(*,*)'fcplCT(1,cur) =', fcplCT(1,cur)
            write(*,*)'fcplCT(2,cur) =', fcplCT(2,cur)
            write(*,*)'fcplCT(3,cur) =', fcplCT(3,cur)
            write(*,*)'fcplCT(4,cur) =', fcplCT(4,cur)
            write(*,*)'  '

            if (qk .eq. 4) then
               stop
            end if
         end if
      end do

*      do cur = 1, 6
*         write(*,*)'  '
*         do i = 1, 4
*            write(*,*)'QQV: k, isig =', cur, i
*            write(*,*)'QQV: fcplCT(isig, k) =', fcplCT(i,cur) 
*         end do
*      end do
*      stop

     
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates and combined the contributions to the qqV 
*** effective couplings in the SM or MSSM

      subroutine qqV(mom2,model,sector,wbfLtr,wbfRtr,wbfL,wbfR,debug)

      implicit none


** Debugging flag
      integer debug

** The effective qqV coupling and tree level contribution
** Labelling: (cur, qk, gen)
** cur       =   1    u + c -> u + H + c
**               2    u + s -> u + H + s
**               3    d + c -> d + H + c
**               4    d + s -> d + H + s
**               5    u + s -> d + H + c
**               6    d + c -> u + H + s
** qk        =   1    q + q        ->   q + H + q
**               2    qbar + q     ->   qbar + H + q
**               3    q + qbar     ->   q + H + qbar
**               4    qbar + qbar  ->   qbar + H + qbar
** NOTE: Be careful of charged current cur-qk labels
** gen       =   1    upper line, first generation
**               2    lower line, second generation        
      double complex wbfL(6,4,2), wbfR(6,4,2)
      double complex wbfLtr(6,4,2), wbfRtr(6,4,2)


** Parts of the effective qqV coupling: tree, loop and counterterm
** Labelling: (qqVcur,gen,ncC)
** qqVcur    =   1    incoming particle is up-type
**               2    incoming particle is down-type
** ncC       =   1    neutral current, quarks
**               2    neutral current, anti-quarks
**               3    charged current, quarks
**               4    charged current, anti-quarks
      double complex cplLtr(2,2,4), cplRtr(2,2,4)
      double complex cplLlp(2,2,4), cplRlp(2,2,4)
      double complex cplLct(2,2,4), cplRct(2,2,4)


** mom2(i) is the momenta^2 of all the quarks and vector bosons, labelling:
**    i     =     1    upper incoming quark
**                2    upper outgoing quark
**                3    lower incoming quark
**                4    lower outgoing quark
**                5    upper vector boson
**                6    lower vector boson
** p12, p22 and p32 are the momenta squared of the quarks and vector boson
** for the qqV current line we're considering (i.e. upper or lower)
      double precision mom2(6)
      double precision p12, p22, p32


** Do loop parameters for running over generation (upper/lower line) and 
** qqV current, quark number (i) , wbf current, quark pattern and neutral /
** charged quark /antiquark qqV pattern (ncC).  Also dummy parameters for
** qqVcur and ncC.  Finally, do loop parameters for neutralinos, charginos
** and sfermions.
      integer gen, qqVcur, i, cur, qk, ncC
      integer qqVcurD, ncCD
      integer neu1, cha1, sfe1, sfe2

** Charge, mass, fermion field renormalisation constants and couplings for 
** considered qqV current.  Labelling: (particle number i)
**    i      =   1   incoming quark
**               2   outgoing quark
      double precision Qq(2), Qmass(2), Qmass2(2)
      double complex dZfL(2), dZfR(2)
      double precision QgL(2), QgR(2), dQgL(2), dQgR(2)

** SUSY couplings and masses for the qqV current being considered
** Labelling: final entry = particle number (i)
**            first entry = neutralino/chargino index
      double complex SNeu1(4,2), SCha1(2,2), SCha2(2,2)
      double precision SNeu2(2), SNeu3(2)
** chargino-neutralino-W coupling 
** labelling: (cnw, chargino, neutralino)
      double complex cnwQ(4,2,4)

** sfermion masses and mixing matrices
** Labelling: sfe1, (sfe2), particle number i
      double precision massSF(2,2), mass2SF(2,2)
      double complex mUS(2,2,2), mUSC(2,2,2)

      integer sector, model


** Parameter and looptools declarations
#include "mssm.inc"


** Changing the sign convention of the SU(2) derivative when working in 
** the MSSM
      if (model .eq. 2) then
         SW = -SW
      end if


** Clearing the total effective coupling
      do cur = 1, 6
         do qk = 1, 4
            do gen = 1, 2
               wbfL(cur,qk,gen) = (0D0,0D0)
               wbfR(cur,qk,gen) = (0D0,0D0)
            end do
         end do
      end do

** Running over qqV current and generation (upper/lower line)
      do qqVcur = 1, 2
         do gen = 1, 2


** Clearing the pieces of the coupling values ...
            do ncC = 1, 4
               cplLtr(qqVcur,gen,ncC) = (0D0,0D0)
               cplRtr(qqVcur,gen,ncC) = (0D0,0D0)
               cplLct(qqVcur,gen,ncC) = (0D0,0D0)
               cplRct(qqVcur,gen,ncC) = (0D0,0D0)
               cplLlp(qqVcur,gen,ncC) = (0D0,0D0)
               cplRlp(qqVcur,gen,ncC) = (0D0,0D0)
            end do


** Setting the momenta squared
            if (gen .eq. 1) then
** upper quark line
               p12 = mom2(1)
               p22 = mom2(2)
               p32 = mom2(5)
            else if (gen .eq. 2) then
** lower quark line
               p12 = mom2(3)
               p22 = mom2(4)
               p32 = mom2(6)
            end if

            p12 = massQ**2
            p22 = massQ**2


            if (debug .eq. 1) then
               write(*,*)'p12, p22, p32 =', p12, p22, p32
            end if

** Setting the values of mass etc for each current: running over particle 
** number
            do i = 1, 2
               Qmass(i) = qqVmass(qqVcur,i,gen)
               Qmass2(i) = Qmass(i)**2
               Qq(i) = qqVcharge(qqVcur,i)

               if (Qq(i) .eq. Qu) then
                  QgL(i) = gL3
                  QgR(i) = gR3
                  dQgL(i) = dgL3
                  dQgR(i) = dgR3
                  if (gen .eq. 1) then
                     dZfL(i) = dZfLu
                     dZfR(i) = dZfRu
                  else if (gen .eq. 2) then
                     dZfL(i) = dZfLc
                     dZfR(i) = dZfRc
                  end if

** sfermion parameters
                  do sfe1 = 1, 2
                     massSF(sfe1,i) = MSf(sfe1,3,gen)
                     mass2SF(sfe1,i) = MSf2(sfe1,3,gen)
                     do sfe2 = 1, 2
                        mUS(sfe1,sfe2,i) = USf(sfe1,sfe2,3,gen)
                        mUSC(sfe1,sfe2,i) = USfC(sfe1,sfe2,3,gen)
                     end do
                  end do
                  
** neutralino couplings
                  do neu1 = 1, 4
                     SNeu1(neu1,i) = SNeu1u(neu1,gen)
                  end do
                  SNeu2(i) = SNeu2u
                  SNeu3(i) = SNeu3u
                  
** chargino couplings
                  do cha1 = 1, 2
                     SCha1(cha1,i) = SCha1u(cha1,gen)
                     SCha2(cha1,i) = SCha2u(cha1)
                  end do

               else if (Qq(i) .eq. Qd) then
                  QgL(i) = gL4
                  QgR(i) = gR4
                  dQgL(i) = dgL4
                  dQgR(i) = dgR4
                  if (gen .eq. 1) then
                     dZfL(i) = dZfLd
                     dZfR(i) = dZfRd
                  else if (gen .eq. 2) then
                     dZfL(i) = dZfLs
                     dZfR(i) = dZfRs
                  end if

** sfermion parameters
                  do sfe1 = 1, 2
                     massSF(sfe1,i) = MSf(sfe1,4,gen)
                     mass2SF(sfe1,i) = MSf2(sfe1,4,gen)
                     do sfe2 = 1, 2
                        mUS(sfe1,sfe2,i) = USf(sfe1,sfe2,4,gen)
                        mUSC(sfe1,sfe2,i) = USfC(sfe1,sfe2,4,gen)
                     end do
                  end do
                  
** neutralino couplings
                  do neu1 = 1, 4
                     SNeu1(neu1,i) = SNeu1d(neu1,gen)
                  end do
                  SNeu2(i) = SNeu2d
                  SNeu3(i) = SNeu3d

** chargino couplings
                  do cha1 = 1, 2
                     SCha1(cha1,i) = SCha1d(cha1,gen)
                     SCha2(cha1,i) = SCha2d(cha1)
                  end do

               end if
            end do

** chargino-neutralino-W couplings
            if (qqVcur .eq. 1) then
               do cha1 = 1, 2
                  do neu1 = 1, 4
                     cnwQ(1,cha1,neu1) = cnw(1,cha1,neu1)
                     cnwQ(2,cha1,neu1) = cnw(2,cha1,neu1)
                     cnwQ(3,cha1,neu1) = cnw(3,cha1,neu1)
                     cnwQ(4,cha1,neu1) = cnw(4,cha1,neu1)
                  end do
               end do
            else if (qqVcur .eq. 2) then
               do cha1 = 1, 2
                  do neu1 = 1, 4
                     cnwQ(1,cha1,neu1) = cnw(4,cha1,neu1)
                     cnwQ(2,cha1,neu1) = cnw(3,cha1,neu1)
                     cnwQ(3,cha1,neu1) = cnw(2,cha1,neu1)
                     cnwQ(4,cha1,neu1) = cnw(1,cha1,neu1)
                  end do
               end do
            end if

** Debugging output: labelling            
            if (debug .ne. 0) then
               if (gen .eq. 1) then
                  write(*,*)'  '
               end if
               write(*,*)'  '
               write(*,*)'qqV CURRENT = ', qqVcur
               write(*,*)'generation (upper/lower line) = ', gen   
c$$$               do i = 1, 2
c$$$                  if (i .eq. 1) then
c$$$                     write(*,*) 'upper line'
c$$$                  else
c$$$                     write(*,*) 'lower line'
c$$$                  end if                     
c$$$                  write(*,*)'Qmass, Qq =', Qmass(i), Qq(i)
c$$$                  write(*,*)'QgL, QgR =', QgL(i), QgR(i)
c$$$               end do
            end if
               

*** FIRSTLY: The tree level
            call tree(qqVcur,gen,cplLtr,cplRtr,QgL,QgR)

** Debugging output: tree level components
            if (debug .ne. 0) then
               write(*,*)'  '
               write(*,*)'Tree level couplings ..'
               write(*,*)'qqZ (L/R) = ', cplLtr(qqVcur,gen,1),
     -              cplRtr(qqVcur,gen,1)
               write(*,*)'qqZbar (L/R) = ', cplLtr(qqVcur,gen,2),
     -              cplRtr(qqVcur,gen,2)
               write(*,*)'qqW (L/R) = ', cplLtr(qqVcur,gen,3),
     -              cplRtr(qqVcur,gen,3)
               write(*,*)'qqWbar (L/R) = ', cplLtr(qqVcur,gen,4),
     -              cplRtr(qqVcur,gen,4)
            end if


*** NEXT : The counterterms
            call counter(qqVcur,gen,cplLct,cplRct,Qq,QgL,QgR,dQgL,dQgR,
     -           Qmass,Qmass2,dZfL,dZfR)

** Debugging output: counterterm components
            if (debug .ne. 0) then
               write(*,*)'  '
               write(*,*)'Counterterm couplings ..'
               write(*,*)'qqZ (L/R) = ', cplLct(qqVcur,gen,1),
     -              cplRct(qqVcur,gen,1)
               write(*,*)'qqZbar (L/R) = ', cplLct(qqVcur,gen,2),
     -              cplRct(qqVcur,gen,2)
               write(*,*)'qqW (L/R) = ', cplLct(qqVcur,gen,3),
     -              cplRct(qqVcur,gen,3)
               write(*,*)'qqWbar (L/R) = ', cplLct(qqVcur,gen,4),
     -              cplRct(qqVcur,gen,4)
            end if

      
*** FINALLY the lopp contributions, which runs iff all particles are included  
*** in the loops. (This subroutine should only run if non-(s)fermions are
*** considered, but you can't be too careful.
            if (sector .ge. 3) then

               call loop(qqVcur,gen,p12,p22,p32,cplLlp,cplRlp,QgL,QgR,
     -              Qq,Qmass,Qmass2,SNeu1,SNeu2,SNeu3,SCha1,SCha2,cnwQ,
     -              massSF,mass2SF,mUS,mUSC)
                     

** Debugging output: loop components
               if (debug .ne. 0) then
                  write(*,*)'  '
                  write(*,*)'Loop couplings ..'
                  write(*,*)'qqZ (L/R) = ', cplLlp(qqVcur,gen,1),
     -                 cplRlp(qqVcur,gen,1)
                  write(*,*)'qqZbar (L/R) = ', cplLlp(qqVcur,gen,2),
     -                 cplRlp(qqVcur,gen,2)
                  write(*,*)'qqW (L/R) = ', cplLlp(qqVcur,gen,3),
     -                 cplRlp(qqVcur,gen,3)
                  write(*,*)'qqWbar (L/R) = ', cplLlp(qqVcur,gen,4),
     -                 cplRlp(qqVcur,gen,4)
               end if


            end if

         end do
      end do

      
** Calculating the total effective coupling
** RECALL:  labelling of cplL is (cur, qk, gen)
**          labelling of cplLtr etc. is (qqVcur, gen, ncC)
      do gen = 1, 2
         do qk = 1, 4
            do cur = 1, 6
               if (cur .le. 4) then
                  if (gen .eq. 1) then
                     if (cur .le. 2) then
                        qqVcurD = 1
                     else
                        qqVcurD = 2
                     end if
                     if (qk .eq. 1 .or. qk .eq. 3) then
                        ncCD = 1
                     else
                        ncCD = 2
                     end if
                  else
                     if (cur .eq. 1 .or. cur .eq. 3) then
                        qqVcurD = 1
                     else
                        qqVcurD = 2
                     end if
                     if (qk .le. 2) then
                        ncCD = 1
                     else
                        ncCD = 2
                     end if
                  end if
               else
                  if (gen .eq. 1) then
                     if (qk .eq. 1 .or. qk .eq. 3) then
                        ncCD = 3
                     else
                        ncCD = 4
                     end if
                  else
                     if (qk .le. 2) then
                        ncCD = 3
                     else
                        ncCD = 4
                     end if
                  end if
                  qqVcurD = cur - 4
                  if (qk .eq. 2 .or. (qk .eq. 1 .and. gen .eq. 2) 
     -                 .or. (qk .eq. 4 .and. gen .eq. 1)) then
                     if (cur .eq. 5) then
                        qqVcurD = 2
                     else
                        qqVcurD = 1
                     end if
                  end if
               end if
               
               wbfLtr(cur,qk,gen) = cplLtr(qqVcurD,gen,ncCD)
               wbfRtr(cur,qk,gen) = cplRtr(qqVcurD,gen,ncCD)
               
               wbfL(cur,qk,gen) = cplLlp(qqVcurD,gen,ncCD) +
     -              cplLct(qqVcurD,gen,ncCD)
               wbfR(cur,qk,gen) = cplRlp(qqVcurD,gen,ncCD) +
     -              cplRct(qqVcurD,gen,ncCD)
               
            end do
         end do
      end do


** Correcting for the sign convention of the SU(2) derivative when working in 
** the MSSM
      if (model .eq. 2) then
         SW = -SW
      end if


      end


********************************************************************************
********************************************************************************

*** This subroutine calculates the tree level contribution to the qqV coupling

      subroutine tree(qqVcur,gen,cplLtr,cplRtr,QgL,QgR)

      implicit none


** These are the tree parts of the effective right / left couplings for uuZ,
** ddZ, udW, .. Labelling: (qqVcur,gen,ncC)
** qqVcur    =   1    incoming particle is up-type
**               2    incoming particle is down-type
** gen       =   1    upper line, first generation
**               2    lower line, second generation
** ncC       =   1    neutral current, quarks
**               2    neutral current, anti-quarks
**               3    charged current, quarks
**               4    charged current, anti-quarks
      double complex cplLtr(2,2,4), cplRtr(2,2,4)


** L/R coupling for the qqV current
      double precision QgL(2), QgR(2)


** Labels of the generation (line number) and qqV current
      integer gen, qqVcur

** Debugging flag
      integer debug

** Parameter declaration
#include "mssm.inc"


** qqZ
      cplLtr(qqVcur,gen,1) = EL*QgL(1)
      cplRtr(qqVcur,gen,1) = EL*QgR(1)
** qqZ bar
      cplLtr(qqVcur,gen,2) = -EL*QgL(1)
      cplRtr(qqVcur,gen,2) = -EL*QgR(1)

** qqW
      cplLtr(qqVcur,gen,3) = EL/(SW*Sqrt(2.D0))
      cplRtr(qqVcur,gen,3) = (0D0,0D0)
** qqW bar
      cplLtr(qqVcur,gen,4) = -EL/(SW*Sqrt(2.D0))
      cplRtr(qqVcur,gen,4) = (0D0,0D0)


      end


********************************************************************************
********************************************************************************

*** This subroutine calls and appropriately assigns the loop components for the
*** qqV current

      subroutine loop(qqVcur,gen,q12,q22,q32,cplLlp,cplRlp,QgL,QgR,Qq,
     -     Qmass,Qmass2,SNeu1,SNeu2,SNeu3,SCha1,SCha2,cnwQ,massSF,
     -     mass2SF,mUS,mUSC)

      implicit none


** These are the loop parts of the effective right / left couplings for uuZ,
** ddZ, udW, .. Labelling: (qqVcur, gen,ncC)
** qqVcur    =   1    incoming particle is up-type
**               2    incoming particle is down-type
** gen       =   1    upper line, first generation
**               2    lower line, second generation
** ncC       =   1    neutral current, quarks
**               2    neutral current, anti-quarks
**               3    charged current, quarks
**               4    charged current, anti-quarks
      double complex cplLlp(2,2,4), cplRlp(2,2,4)


** Effective couplings, labelling (L/R)
      double complex qqZlp(2), qqWlp(2)
      double complex qqZbarlp(2), qqWbarlp(2)


** Squared momenta of the quarks and vector boson
      double precision q12, q22, q32


** Values of coupling, charges, mass for a specific current.  Labelling:
** (particle number = incoming 1 / outgoing 2)
      double precision QgL(2), QgR(2), Qq(2), Qmass(2), Qmass2(2)

** SUSY couplings and masses for the qqV current being considered
** Labelling: final entry = particle number (i)
**            first entry = neutralino/chargino index
      double complex SNeu1(4,2), SCha1(2,2), SCha2(2,2)
      double precision SNeu2(2), SNeu3(2)
** chargino-neutralino-W coupling 
** labelling: (cnw, chargino, neutralino)
      double complex cnwQ(4,2,4)

** sfermion masses and mixing matrices
** Labelling: sfe1, (sfe2), particle number i
      double precision massSF(2,2), mass2SF(2,2)
      double complex mUS(2,2,2), mUSC(2,2,2)

** qqV current and generation (line) label
      integer qqVcur, gen

** Debugging flag
      integer debug

#include "global.inc"
#include "mssm.inc"


** Clearing the loop values
      qqZlp(1) = (0D0,0D0)
      qqZlp(2) = (0D0,0D0)
      qqWlp(1) = (0D0,0D0)
      qqWlp(2) = (0D0,0D0)


      call qqZ(qqVcur,qqZlp,q12,q22,q32,QgL,QgR,Qq,Qmass,Qmass2)
      call qqW(qqVcur,qqWlp,q12,q22,q32,QgL,QgR,Qq,Qmass,Qmass2)

** The SUSY diagrams
      if ((model .eq. 2) .and. (sector .ge. 5)) then
         SW = -SW
         call qqZsusy(qqVcur,qqZlp,q12,q22,q32,QgL,QgR,Qmass,Qmass2,
     -        SNeu1,SNeu2,SNeu3,SCha1,SCha2,massSF,mass2SF,mUS,mUSC)
         call qqWsusy(qqVcur,qqWlp,q12,q22,q32,QgL,QgR,Qmass,Qmass2,
     -        SNeu1,SNeu2,SNeu3,SCha1,SCha2,cnwQ,massSF,mass2SF,mUS,
     -        mUSC)
         SW = -SW
      end if

      
      cplLlp(qqVcur,gen,1) = qqZlp(1)
      cplRlp(qqVcur,gen,1) = qqZlp(2)

      cplLlp(qqVcur,gen,2) = -qqZlp(1)
      cplRlp(qqVcur,gen,2) = -qqZlp(2)

      cplLlp(qqVcur,gen,3) = qqWlp(1)
      cplRlp(qqVcur,gen,3) = qqWlp(2)

      cplLlp(qqVcur,gen,4) = -qqWlp(1)
      cplRlp(qqVcur,gen,4) = -qqWlp(2)


      end


********************************************************************************
********************************************************************************

*** This subroutine calculates the loop contribution to the effective L/R 
*** coupling for the qqZ vertex in the SM

      subroutine qqZ(qqVcur,qqZlp,q12,q22,q32,QgL,QgR,Qq,Qmass,Qmass2)

      implicit none


** Effective coupling from loops for qqZ coupling Labelling (L/R, gen)
      double complex qqZlp(2)


** Couplings, charges and masses for particle number (i) for the qqV current
** being considered
      double precision QgL(2), QgR(2), Qq(2), Qmass(2), Qmass2(2)


** Value of separate loop contributions
      double complex loopNL, loopNR, loopPL, loopPR

** Integral labels
      integer iint1, iint2, iint3, iint4, iint5, iint6, iint7

** EM coupling
      double precision Alfa, Alfa2


** Squared momenta of the quarks and vector boson
      double precision q12, q22, q32


** Label for the qqV current
      integer qqVcur

** Debugging flag
      integer debug


** Parameter and looptools declarations
#include "looptools.h"
#include "mssm.inc"
#include "global.inc"


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** Clearing loop values
      loopNL = (0D0,0D0)
      loopPL = (0D0,0D0)
      loopNR = (0D0,0D0)
      loopPR = (0D0,0D0)


** Integral values
      iint1 = Bget(q32,Qmass2(1),Qmass2(1))
      iint2 = Bget(q32,Qmass2(2),Qmass2(2))   
      iint3 = Cget(Qmass2(1),q32,Qmass2(1),0.D0,Qmass2(1),Qmass2(1))
      iint4 = Cget(Qmass2(1),q32,Qmass2(1),MW2,Qmass2(2),Qmass2(2))
      iint5 = Cget(Qmass2(1),q32,Qmass2(1),MZ2,Qmass2(1),Qmass2(1))
      iint6 = Bget(q32,MW2,MW2)
      iint7 = Cget(q32,Qmass2(1),Qmass2(1),MW2,MW2,Qmass2(2))


** First: left handed!
      loopNL = loopNL + (Alfa*CW)/(8.D0*Pi*SW2)*(EL)/SW 

      loopNL = loopNL - 
     -     (Alfa*CW)/(4.D0*Pi*SW2)*(EL*Bval(bb0,iint6))/SW 

      loopNL = loopNL - 
     -     (Alfa*CW)/(2.D0*Pi*SW2)*(EL*Cval(cc00,iint7))/SW 

      loopNL = loopNL + 
     -     1/(8.D0*Pi*SW2)*(Cval(cc2,iint7)*
     -     (5*Alfa*CW*(EL*Qmass2(1)) - 
     -     2*Alfa*CW*(EL*q32)))/SW 

      loopNL = loopNL + 
     -     1/(8.D0*Pi*SW2)*(Cval(cc0,iint7)*
     -     (5*Alfa*CW*(EL*Qmass2(1)) + 
     -     Alfa*CW*(EL*(-2*Qmass2(2) - 2*q32))))/SW
      
      loopPL = loopPL + 
     -     Cval(cc0,iint4)*(Alfa/(2.D0*Pi*SW2)*(EL*QgL(2)*Qmass2(1)) - 
     -     Alfa/(4.D0*Pi*SW2)*
     -     (EL*(QgR(2)*Qmass2(2) + QgL(2)*(q32 - MW2)))) + 
     -     Cval(cc0,iint5)*(Alfa/Pi*(EL*QgL(1)**3*Qmass2(1)) + 
     -     Alfa/(2.D0*Pi)*(EL*QgR(1)**3*Qmass2(1)) - 
     -     Alfa/(2.D0*Pi)*(EL*QgL(1)**2*
     -     (QgR(1)*Qmass2(1) + QgL(1)*(q32 - MZ2)))) + 
     -     Cval(cc0,iint3)*(Alfa/Pi*
     -     (EL*QgL(1)*Qmass2(1)*Qq(1)**2) + 
     -     Alfa/(2.D0*Pi)*(EL*QgR(1)*Qmass2(1)*Qq(1)**2) - 
     -     Alfa/(2.D0*Pi)*(EL*(QgR(1)*Qmass2(1) + QgL(1)*q32)*
     -     Qq(1)**2)) + 
     -     Alfa/(4.D0*Pi*SW2)*(EL*QgL(2)*Bval(bb0,iint2)) - 
     -     Alfa/(2.D0*Pi*SW2)*(EL*QgL(2)*Cval(cc00,iint4)) + 
     -     Alfa*(EL*Cval(cc00,iint5)*
     -     (-(1/Pi*(QgL(1)**3)))) + 
     -     1/(4.D0*Pi*SW2)*(Cval(cc1,iint4)*
     -     (3*Alfa*(EL*QgL(2)*Qmass2(1)) - 
     -     Alfa*(EL*QgL(2)*q32))) + 
     -     Alfa*(EL*Bval(bb0,iint1)*
     -     (QgL(1)*(1/(2.D0*Pi)*QgL(1)**2 + 1/(2.D0*Pi)*Qq(1)**2))) - 
     -     Alfa/(8.D0*Pi*SW2)*(EL*
     -     (QgL(2) + QgL(1)*(2*SW2*QgL(1)**2 + 2*SW2*Qq(1)**2))) + 
     -     Alfa*(EL*Cval(cc00,iint3)*
     -     (-(1/Pi*(QgL(1)*Qq(1)**2)))) + 
     -     Cval(cc1,iint3)*((3*Alfa)/(2.D0*Pi)*(EL*QgL(1)*
     -     Qmass2(1)*Qq(1)**2) + 
     -     Alfa/(2.D0*Pi)*(EL*QgR(1)*Qmass2(1)*Qq(1)**2) +
     -     Alfa*(EL*(-(1/(2.D0*Pi)*(QgL(1)*q32*Qq(1)**2)))))
      
      loopPL = loopPL + 
     -     1/(4.D0*Pi*SW2)*(Cval(cc2,iint4)*
     -     (3*Alfa*(EL*QgL(2)*Qmass2(1)) - 
     -     Alfa*(EL*QgL(2)*q32))) + 
     -     Cval(cc1,iint5)*((3*Alfa)/(2.D0*Pi)*(EL*
     -     QgL(1)**3*Qmass2(1)) + 
     -     Alfa/(2.D0*Pi)*(EL*QgR(1)**3*Qmass2(1)) + 
     -     Alfa*(EL*(-(1/(2.D0*Pi)*(QgL(1)**3*q32))))) + 
     -     Cval(cc2,iint5)*((3*Alfa)/(2.D0*Pi)*(EL*QgL(1)**3*Qmass2(1)) + 
     -     Alfa/(2.D0*Pi)*(EL*QgR(1)**3*Qmass2(1)) + 
     -     Alfa*(EL*(-(1/(2.D0*Pi)*(QgL(1)**3*q32))))) + 
     -     Cval(cc2,iint3)*((3*Alfa)/(2.D0*Pi)*(EL*QgL(1)*Qmass2(1)*Qq(1)**2) + 
     -     Alfa/(2.D0*Pi)*(EL*QgR(1)*Qmass2(1)*Qq(1)**2) +
     -     Alfa*(EL*(-(1/(2.D0*Pi)*(QgL(1)*q32*Qq(1)**2)))))
      
** Right handed
      loopNR = loopNR + 
     -     1/(8.D0*Pi*SW2)*(Cval(cc2,iint7)*
     -     (5*Alfa*CW*((3/5.D0)*EL*Qmass2(1))))/SW 

      loopNR = loopNR + 
     -     1/(8.D0*Pi*SW2)*(Cval(cc0,iint7)*
     -     (5*Alfa*CW*((3/5.D0)*EL*Qmass2(1))))/SW
      
      loopPR = loopPR + 
     -     Cval(cc0,iint4)*(Alfa/(2.D0*Pi*SW2)*((1/2.D0)*
     -     EL*QgL(2)*Qmass2(1))) + 
     -     Cval(cc0,iint5)*(Alfa/Pi*((1/2.D0)*EL*QgL(1)**3*Qmass2(1)) + 
     -     Alfa/(2.D0*Pi)*(2*EL*QgR(1)**3*Qmass2(1)) + 
     -     Alfa/(2.D0*Pi)*(EL*QgR(1)**2*
     -     (-(QgL(1)*Qmass2(1)) + QgR(1)*(-q32 + MZ2)))) + 
     -     Cval(cc0,iint3)*(Alfa/Pi*
     -     ((1/2.D0)*EL*QgL(1)*Qmass2(1)*Qq(1)**2) + 
     -     Alfa/(2.D0*Pi)*(2*EL*QgR(1)*Qmass2(1)*Qq(1)**2) - Alfa/(2.D0*Pi)*
     -     (EL*(QgL(1)*Qmass2(1) + QgR(1)*q32)*Qq(1)**2)) - 
     -     Alfa/(4.D0*Pi)*(EL*QgR(1)*(QgR(1)**2 + Qq(1)**2)) + 
     -     Alfa*(EL*Cval(cc00,iint5)*
     -     (-1/Pi*(QgR(1)**3))) + 
     -     1/(4.D0*Pi*SW2)*(Cval(cc1,iint4)*
     -     (3*Alfa*((1/3.D0)*EL*QgL(2)*Qmass2(1)))) + 
     -     Alfa*(EL*Bval(bb0,iint1)*
     -     (QgR(1)*(1/(2.D0*Pi)*QgR(1)**2 + 1/(2.D0*Pi)*Qq(1)**2))) + 
     -     Alfa*(EL*Cval(cc00,iint3)*
     -     (-1/Pi*(QgR(1)*Qq(1)**2))) + 
     -     Cval(cc1,iint3)*((3*Alfa)/(2.D0*Pi)*((1/3.D0)*EL*QgL(1)*
     -     Qmass2(1)*Qq(1)**2) + 
     -     Alfa/(2.D0*Pi)*(3*EL*QgR(1)*Qmass2(1)*Qq(1)**2) +
     -     Alfa*(EL*(-1/(2.D0*Pi)*(QgR(1)*q32*Qq(1)**2))))
      
      loopPR = loopPR + 
     -     1/(4.D0*Pi*SW2)*(Cval(cc2,iint4)*
     -     (3*Alfa*((1/3.D0)*EL*QgL(2)*Qmass2(1)))) + 
     -     Cval(cc1,iint5)*((3*Alfa)/(2.D0*Pi)*((1/3.D0)*EL*
     -     QgL(1)**3*Qmass2(1)) + 
     -     Alfa/(2.D0*Pi)*(3*EL*QgR(1)**3*Qmass2(1)) + 
     -     Alfa*(EL*(-1/(2.D0*Pi)*(QgR(1)**3*q32)))) + 
     -     Cval(cc2,iint5)*((3*Alfa)/(2.D0*Pi)*((1/3.D0)*
     -     EL*QgL(1)**3*Qmass2(1)) + 
     -     Alfa/(2.D0*Pi)*(3*EL*QgR(1)**3*Qmass2(1)) + 
     -     Alfa*(EL*(-1/(2.D0*Pi)*(QgR(1)**3*q32)))) + 
     -     Cval(cc2,iint3)*((3*Alfa)/(2.D0*Pi)*((1/3.D0)*
     -     EL*QgL(1)*Qmass2(1)*Qq(1)**2) + 
     -     Alfa/(2.D0*Pi)*(3*EL*QgR(1)*Qmass2(1)*Qq(1)**2) +
     -     Alfa*(EL*(-1/(2.D0*Pi)*(QgR(1)*q32*Qq(1)**2))))
      
      
      if (qqVcur .eq. 1) then
         loopNL = -loopNL
         loopNR = -loopNR
      end if


** Set individual parts to zero for debugging
*      loopNL = (0D0,0D0)
*      loopPL = (0D0,0D0)
*      loopNR = (0D0,0D0)
*      loopPR = (0D0,0D0)

         
      qqZlp(1) = loopNL + loopPL
      qqZlp(2) = loopNR + loopPR      


      end


********************************************************************************
********************************************************************************

*** This subroutine calculates the loop contribution to the effective L/R 
*** coupling for the qqW vertex in the SM

      subroutine qqW(qqVcur,qqWlp,q12,q22,q32,QgL,QgR,Qq,Qmass,Qmass2)

      implicit none

** Effective coupling from loops for qqW coupling Labelling (L/R)
      double complex qqWlp(2)


** Couplings, charges and masses for particle number (i) for the qqV current
** being considered
      double precision QgL(2), QgR(2), Qq(2), Qmass(2), Qmass2(2)


** Value of separate loop contributions
      double complex loopNL, loopNR, loopPL, loopPR

** Integral labels
      integer iint1, iint2, iint3, iint4, iint5
      integer iint6, iint7, iint8, iint9

** EM coupling
      double precision Alfa, Alfa2

** Squared momenta of the quarks and vector boson
      double precision q12, q22, q32

      double complex dup1, dup2


** Label for the qqV current
      integer qqVcur

** Debugging flag
      integer debug


** Looptools and parameters declarations
#include "looptools.h"
#include "mssm.inc"
#include "global.inc"


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** Clearing loop values
      loopNL = (0D0,0D0)
      loopPL = (0D0,0D0)
      loopNR = (0D0,0D0)
      loopPR = (0D0,0D0)


** Integer values
      iint1 = Bget(q32,Qmass2(1),Qmass2(2))
      iint2 = Cget(Qmass2(1),q32,Qmass2(2),0.D0,Qmass2(1),Qmass2(2))
      iint3 = Cget(Qmass2(1),q32,Qmass2(2),MZ2,Qmass2(1),Qmass2(2))
      iint4 = Bget(q32,0.D0,MW2)
      iint5 = Bget(q32,MW2,MZ2)
      iint6 = Cget(q32,Qmass2(1),Qmass2(2),0.D0,MW2,Qmass2(2))
      iint7 = Cget(q32,Qmass2(2),Qmass2(1),MW2,MZ2,Qmass2(2))
      iint8 = Cget(q32,Qmass2(1),Qmass2(2),MW2,MZ2,Qmass2(1))
      iint9 = Cget(q32,Qmass2(2),Qmass2(1),0.D0,MW2,Qmass2(1))


** First: left handed!     
      loopNL = loopNL + 
     -     Cval(cc1,iint7)*(-(Alfa*CW)/(4.D0*Pi*SW2)*
     -     (EL*QgL(2)*(-Qmass2(1) + Qmass2(2)))/Sqrt(2.D0)) + 
     -     Cval(cc1,iint6)*(Alfa/(4.D0*Pi)*(EL*(-Qmass2(1) + 
     -     Qmass2(2))*Qq(2))/
     -     (SW*Sqrt(2.D0))) + 
     -     Cval(cc0,iint6)*(-Alfa/(2.D0*Pi)*(EL*
     -     (-Qmass2(1) + Qmass2(2) + q32)*
     -     Qq(2))/(SW*Sqrt(2.D0))) - 
     -     Alfa/(2.D0*Pi)*(EL*Qq(2)*Bval(bb0,iint4))/(SW*Sqrt(2.D0)) - 
     -     (Alfa*CW)/(2.D0*Pi*SW2)*
     -     (EL*QgL(2)*Bval(bb0,iint5))/Sqrt(2.D0) - 
     -     Alfa/Pi*(EL*Qq(2)*Cval(cc00,iint6))/(SW*Sqrt(2.D0)) - 
     -     (Alfa*CW)/(Pi*SW2)*
     -     (EL*QgL(2)*Cval(cc00,iint7))/Sqrt(2.D0) + 
     -     Cval(cc2,iint6)*(-Alfa/(2.D0*Pi)*(EL*Qq(2)*
     -     (-Qmass2(1) + q32 - 3/2.D0*Qmass2(2)))/
     -     (SW*Sqrt(2.D0))) + 
     -     Cval(cc0,iint7)*(-(3*Alfa*CW)/(4.D0*Pi*SW2)*
     -     (EL*(QgR(2)*Qmass2(2) + 
     -     QgL(2)*(-Qmass2(1) + 2/3.D0*q32)))/Sqrt(2.D0)) + 
     -     Alfa/(4.D0*Pi*SW2)*(EL*(SW2*Qq(2) + CW*(QgL(2)*SW)))/
     -     (SW*Sqrt(2.D0)) + 
     -     Cval(cc2,iint7)*(-(Alfa*CW)/(2.D0*Pi*SW2)*
     -     (EL*QgL(2)*(-Qmass2(2) + q32 - 3/2.D0*Qmass2(1)))/
     -     Sqrt(2.D0))


      loopNL = loopNL + 
     -     Cval(cc1,iint8)*(-(Alfa*CW)/(4.D0*Pi*SW2)*
     -     (EL*QgL(1)*(-Qmass2(1) + Qmass2(2)))/Sqrt(2.D0)) + 
     -     Cval(cc1,iint9)*(Alfa/(4.D0*Pi)*(EL*(-Qmass2(1) + 
     -     Qmass2(2))*Qq(1))/
     -     (SW*Sqrt(2.D0))) + 
     -     Alfa/(2.D0*Pi)*(EL*Qq(1)*Bval(bb0,iint4))/(SW*Sqrt(2.D0)) + 
     -     (Alfa*CW)/(2.D0*Pi*SW2)*
     -     (EL*QgL(1)*Bval(bb0,iint5))/Sqrt(2.D0) + 
     -     (Alfa*CW)/(Pi*SW2)*
     -     (EL*QgL(1)*Cval(cc00,iint8))/Sqrt(2.D0) + 
     -     Alfa/Pi*(EL*Qq(1)*Cval(cc00,iint9))/(SW*Sqrt(2.D0)) + 
     -     Cval(cc0,iint9)*(Alfa/(2.D0*Pi)*(EL*(Qmass2(1) - Qmass2(2) + q32)*
     -     Qq(1))/(SW*Sqrt(2.D0))) + 
     -     Cval(cc2,iint8)*((Alfa*CW)/(2.D0*Pi*SW2)*
     -     (EL*QgL(1)*(-Qmass2(1) + q32 - 3/2.D0*Qmass2(2)))/
     -     Sqrt(2.D0)) + Cval(cc0,iint8)*
     -     ((3*Alfa*CW)/(4.D0*Pi*SW2)*
     -     (EL*(QgR(1)*Qmass2(1) + 
     -     QgL(1)*(-Qmass2(2) + 2/3.D0*q32)))/Sqrt(2.D0)) - 
     -     Alfa/(4.D0*Pi*SW2)*(EL*(SW2*Qq(1) + CW*(QgL(1)*SW)))/
     -     (SW*Sqrt(2.D0)) + 
     -     Cval(cc2,iint9)*(Alfa/(2.D0*Pi)*(EL*Qq(1)*
     -     (-Qmass2(2) + q32 - 3/2.D0*Qmass2(1)))/
     -     (SW*Sqrt(2.D0)))


      dup1 = -Qmass2(1) - 2*Qmass2(2) + q32
      dup2 = -2*Qmass2(1) - Qmass2(2) + q32

      loopPL = loopPL + 
     -     Cval(cc0,iint3)*(Alfa/(2.D0*Pi)*(EL*0*
     -     (-(QgL(1)*QgL(2)) + QgR(1)*QgR(2))*Qmass(1)*
     -     Qmass(2))/(SW*Sqrt(2.D0)) + 
     -     Alfa/(2.D0*Pi)*(EL*QgL(1)*QgL(2)*
     -     (Qmass2(1) + Qmass2(2) - q32 + MZ2))/
     -     (SW*Sqrt(2.D0))) - 
     -     Alfa/(4.D0*Pi)*(EL*(QgL(1)*QgL(2) + Qq(1)*Qq(2)))/(SW*Sqrt(2.D0)) + 
     -     Alfa/(2.D0*Pi)*(EL*(QgL(1)*QgL(2) + Qq(1)*Qq(2))*Bval(bb0,iint1))/
     -     (SW*Sqrt(2.D0)) - 
     -     Alfa/(2.D0*Pi)*(EL*(-Qmass2(1) - Qmass2(2) + q32)*Qq(1)*
     -     Qq(2)*Cval(cc0,iint2))/(SW*Sqrt(2.D0)) - 
     -     Alfa/Pi*(EL*Qq(1)*Qq(2)*Cval(cc00,iint2))/(SW*Sqrt(2.D0)) - 
     -     Alfa/Pi*(EL*QgL(1)*QgL(2)*Cval(cc00,iint3))/(SW*Sqrt(2.D0)) + 
     -     1/(2.D0*Pi)*(Cval(cc2,iint2)*
     -     (Alfa*(-(dup1*EL*Qq(1)*Qq(2)))
     -     ))/(SW*Sqrt(2.D0)) + 
     -     1/(2.D0*Pi)*(Cval(cc2,iint3)*
     -     (Alfa*(-(dup1*EL*QgL(1)*QgL(2)))))/(SW*Sqrt(2.D0)) + 
     -     1/(2.D0*Pi)*(Cval(cc1,iint3)*
     -     (Alfa*(-(dup2*EL*QgL(1)*QgL(2)))))/(SW*Sqrt(2.D0)) + 
     -     1/(2.D0*Pi)*(Cval(cc1,iint2)*
     -     (Alfa*(-(dup2*EL*Qq(1)*Qq(2)))))/(SW*Sqrt(2.D0))


** Right handed
      loopNR = loopNR + 
     -     Cval(cc2,iint6)*((3*Alfa)/(4.D0*Pi)*
     -     (EL*Qmass(1)*Qmass(2)*Qq(2))/(SW*Sqrt(2.D0))) + 
     -     Cval(cc0,iint7)*(-(3*Alfa*CW)/(4.D0*Pi*SW2)*
     -     (EL*(-QgL(2) + QgR(2))*Qmass(1)*Qmass(2))/Sqrt(2.D0)) + 
     -     Cval(cc2,iint7)*((3*Alfa*CW)/(4.D0*Pi*SW2)*
     -     (EL*QgL(2)*Qmass(1)*Qmass(2))/Sqrt(2.D0))

      loopNR = loopNR + 
     -     Cval(cc2,iint8)*(-((3*Alfa*CW)/(4.D0*Pi*SW2)*
     -     (EL*QgL(1)*Qmass(1)*Qmass(2))/Sqrt(2.D0))) + Cval(cc0,iint8)*
     -     ((3*Alfa*CW)/(4.D0*Pi*SW2)*
     -     (EL*(-QgL(1) + QgR(1))*Qmass(1)*Qmass(2))/Sqrt(2.D0)) - 
     -     Alfa/(4.D0*Pi*SW2)*(EL*0*(SW2*Qq(1) + CW*(QgL(1)*SW)))/
     -     (SW*Sqrt(2.D0)) + 
     -     Cval(cc2,iint9)*(-(3*Alfa)/(4.D0*Pi)*
     -     (EL*Qmass(1)*Qmass(2)*Qq(1))/(SW*Sqrt(2.D0)))
      
      
      loopPR = loopPR + 
     -     Cval(cc0,iint3)*(Alfa/(2.D0*Pi)*(EL*
     -     (-(QgL(1)*QgL(2)) + QgR(1)*QgR(2))*Qmass(1)*
     -     Qmass(2))/(SW*Sqrt(2.D0))) + 
     -     1/(2.D0*Pi)*(Cval(cc2,iint2)*
     -     (Alfa*(EL*Qmass(1)*
     -     Qmass(2)*Qq(1)*Qq(2))
     -     ))/(SW*Sqrt(2.D0)) + 
     -     1/(2.D0*Pi)*(Cval(cc2,iint3)*
     -     (Alfa*(EL*QgL(1)*QgL(2)*Qmass(1)*Qmass(2))))/(SW*Sqrt(2.D0)) + 
     -     1/(2.D0*Pi)*(Cval(cc1,iint3)*
     -     (Alfa*(EL*QgL(1)*QgL(2)*Qmass(1)*Qmass(2))))/(SW*Sqrt(2.D0)) + 
     -     1/(2.D0*Pi)*(Cval(cc1,iint2)*
     -     (Alfa*(EL*Qmass(1)*Qmass(2)*Qq(1)*Qq(2))))/(SW*Sqrt(2.D0))


** Set individual parts to zero for debugging
*      loopNL = (0D0,0D0)
*      loopPL = (0D0,0D0)
*      loopNR = (0D0,0D0)
*      loopPR = (0D0,0D0)
         
      if (qqVcur .eq. 2) then
         loopNL = -loopNL
         loopNR = -loopNR
      end if

         
      qqWlp(1) = loopNL + loopPL
      qqWlp(2) = loopNR + loopPR  


      end


********************************************************************************
********************************************************************************

*** This subroutine calculates counterterm contributions to the effective
*** coupling for all qqV currents

      subroutine counter(qqVcur,gen,cplLct,cplRct,Qq,QgL,QgR,dQgL,dQgR,
     -     Qmass,Qmass2,dZfL,dZfR)

      implicit none


** Effective coupling from counterterms. 
** Labelling: (qqVcur,gen,ncC)
** qqVcur    =   1    incoming particle is up-type
**               2    incoming particle is down-type
** gen       =   1    upper line, first generation
**               2    lower line, second generation  
** ncC       =   1    neutral current, quarks
**               2    neutral current, anti-quarks
**               3    charged current, quarks
**               4    charged current, anti-quarks
      double complex cplLct(2,2,4), cplRct(2,2,4)


** Labels for qqV current and generation (upper / lower line)
      integer qqVcur, gen


** Charge, mass, fermion field renormalisation constants and couplings for 
** considered qqV current.  Labelling: (particle number i)
**    i      =   1   incoming quark
**               2   outgoing quark
      double precision Qq(2), QgL(2), QgR(2), dQgL(2), dQgR(2)
      double complex dZfL(2), dZfR(2)
      double precision Qmass(2), Qmass2(2)

** Debugging flag
      integer debug


** Parameter declarations
#include "mssm.inc" 
#include "global.inc"


** qqZ
      cplLct(qqVcur,gen,1) = 1/2.D0*(EL*QgL(1)*DCONJG(dZfL(1))) +1/2.D0*
     -     (EL*((dZfL(1)*QgL(1) + dZZZ1*QgL(1) -dZAZ1*Qq(1)+2*dQgL(1))))

      cplRct(qqVcur,gen,1) = 1/2.D0*(EL*QgR(1)*DCONJG(dZfR(1))) +1/2.D0*
     -     (EL*((dZZZ1*QgR(1) - dZAZ1*Qq(1) +2*dQgR(1)+QgR(1)*dZfR(1))))

** qqZ bar
      cplLct(qqVcur,gen,2) = -cplLct(qqVcur,gen,1)
      cplRct(qqVcur,gen,2) = -cplRct(qqVcur,gen,1)

** qqW
      cplLct(qqVcur,gen,3) = (EL*(1/2.D0*DCONJG(dZfL(2))))/(SW*
     -     Sqrt(2.D0)) + 1/2.D0*(EL*dZfL(1))/(SW*Sqrt(2.D0)) + 1/(2.D0*
     -     SW2)*(EL*(-2*dSW1*SW+2*SW2*dZet1+SW2*dZW1))/(SW*Sqrt(2.D0))

      cplRct(qqVcur,gen,3) = (0D0,0D0)

** qqW bar 
      cplLct(qqVcur,gen,4) = -cplLct(qqVcur,gen,3)
      cplRct(qqVcur,gen,4) = (0D0,0D0)


      end


********************************************************************************
********************************************************************************

*** This subroutine calculates the SUSY loop contribution to the effective L/R 
*** coupling for the qqZ vertex.  NOTE that we can split the loop diagrams into
*** SUSY and SM parts because we're pretending that external particles have no
*** mass and thus do not couple to the higgs

      subroutine qqZsusy(qqVcur,qqZlp,q12,q22,q32,QgL,QgR,Qmass,Qmass2,
     -     SNeu1,SNeu2,SNeu3,SCha1,SCha2,massSF,mass2SF,mUS,mUSC)

      implicit none


** Effective coupling from loops for qqZ coupling Labelling (L/R, gen)
      double complex qqZlp(2)


** SM couplings and masses for particle number (i) for the qqV current
** being considered
      double precision QgL(2), QgR(2), Qmass(2), Qmass2(2)

** SUSY couplings and masses  for the qqV current being considered
** Labelling: final entry = particle number (i)
**            first entry = neutralino/chargino index
      double complex SNeu1(4,2), SCha1(2,2), SCha2(2,2)
      double precision SNeu2(2), SNeu3(2)

** sfermion masses and mixing matrices
** Labelling: sfe1, (sfe2), particle number i
      double precision massSF(2,2), mass2SF(2,2)
      double complex mUS(2,2,2), mUSC(2,2,2)

** do loop variables
      integer sfe1, sfe2, neu1, neu2, cha1, cha2


** Value of separate loop contributions
      double complex loopNL, loopNR, loopPL, loopPR

** Integral labels
      integer iint1(4,2,2), iint2(4,4,2), iint3(4,4), iint4(2,2,2)
      integer iint5(2,2,2), iint6(2,2)

** EM coupling
      double precision Alfa, Alfa2


** do loop variables
      integer Sfe4, Sfe5, Neu4, Neu5, Cha4, Cha5

** dummy parameters
      double complex dup1, dup2, dup3, dup4, dup5, dup6, dup7
      double complex dup8, dup9,  dup10, dup11, dup12, dup13, dup14
      double complex dup15, dup16, dup17, dup18, dup19, dup20, dup21
      double complex dup22, dup23, dup24, dup25, dup26, dup27, dup28
      double complex dup29, dup30
      double complex tmp1, tmp2, tmp3


** Squared momenta of the quarks and vector boson
      double precision q12, q22, q32


** Label for the qqV current
      integer qqVcur


      double precision IndexDelta(4,4)

** Debugging flag
      integer debug


** Parameter and looptools declarations
#include "looptools.h"
#include "mssm.inc"
#include "global.inc"


      Alfa = AlfaQED
      Alfa2 = AlfaQED2

      do neu1 = 1, 4
         do neu2 = 1, 4
            if (neu1 .eq. neu2) then
               IndexDelta(neu1,neu2) = (1D0,0D0)
            else
               IndexDelta(neu1,neu2) = (0D0,0D0)
            end if
         end do
      end do


** Clearing loop values
      loopNL = (0D0,0D0)
      loopPL = (0D0,0D0)
      loopNR = (0D0,0D0)
      loopPR = (0D0,0D0)


** Integral values
      do Cha5 = 1,2
         do Cha4 = 1,2
            iint6(Cha4,Cha5) = Bget(q32,MCha2(Cha4),MCha2(Cha5))
         enddo
      enddo
      do Neu5 = 1,4
         do Neu4 = 1,4
            iint3(Neu4,Neu5) = Bget(q32,MNeu2(Neu4),MNeu2(Neu5))
         enddo
      enddo
      do Sfe4 = 1,2
         do Cha5 = 1,2
            do Cha4 = 1,2
               iint5(Cha4,Cha5,Sfe4) = Cget(q32,Qmass2(1),Qmass2(1),
     -              MCha2(Cha4),MCha2(Cha5),mass2SF(Sfe4,2))
            enddo
         enddo
      enddo
      do Sfe5 = 1,2
         do Sfe4 = 1,2
            do Cha4 = 1,2
               iint4(Cha4,Sfe4,Sfe5) = Cget(Qmass2(1),q32,Qmass2(1),
     -              MCha2(Cha4),mass2SF(Sfe4,2),mass2SF(Sfe5,2))
            enddo
         enddo
      enddo
      do Sfe4 = 1,2
         do Neu5 = 1,4
            do Neu4 = 1,4
               iint2(Neu4,Neu5,Sfe4) = Cget(q32,Qmass2(1),Qmass2(1),
     -              MNeu2(Neu4),MNeu2(Neu5),mass2SF(Sfe4,1))
            enddo
         enddo
      enddo
      do Sfe5 = 1,2
         do Sfe4 = 1,2
            do Neu4 = 1,4
               iint1(Neu4,Sfe4,Sfe5) = Cget(Qmass2(1),q32,Qmass2(1),
     -              MNeu2(Neu4),mass2SF(Sfe4,1),mass2SF(Sfe5,1))
            enddo
         enddo
      enddo


** Calculating the L/R couplings
      do Sfe4 = 1,2
         do Cha5 = 1,2
            do Cha4 = 1,2

               tmp1 = (0D0,0D0)
               tmp2 = (0D0,0D0)
               tmp3 = (0D0,0D0)

               dup28 = -Qmass2(1) + mass2SF(Sfe4,2)
               dup29 = SW2*IndexDelta(Cha4,Cha5) - 1/2.D0*(UCha(Cha4,2)*
     -              UChaC(Cha5,2)) - UCha(Cha4,1)*UChaC(Cha5,1)
               dup30 = SW2*IndexDelta(Cha4,Cha5) - 1/2.D0*(VCha(Cha5,2)*
     -              VChaC(Cha4,2)) - VCha(Cha5,1)*VChaC(Cha4,1)

               if (qqVcur .eq. 2) then
                  dup30 = SW2*IndexDelta(Cha4,Cha5) - 1/2.D0*
     -                 (UCha(Cha5,2)*UChaC(Cha4,2)) - 
     -                 UCha(Cha5,1)*UChaC(Cha4,1)
                  dup29 = SW2*IndexDelta(Cha4,Cha5) - 
     -                 1/2.D0*(VCha(Cha4,2)*VChaC(Cha5,2)) - 
     -                 VCha(Cha4,1)*VChaC(Cha5,1)
               end if


               tmp1 =1/(8.D0*CW*MW*MW2*Pi*SW2)*(Alfa*(dup30*EL*Qmass(1)*
     -              (2*MW2*(DCONJG(SCha1(Cha5,1))*DCONJG(SCha2(Cha4,1))*
     -              MCha(Cha4)*mUS(Sfe4,1,2)*mUSC(Sfe4,1,2)) - 
     -              MW*(DCONJG(SCha1(Cha4,2))*DCONJG(SCha1(Cha5,1))*
     -              Sqrt(2.D0)*MCha(Cha4)*mUS(Sfe4,2,2)*
     -              mUSC(Sfe4,1,2)))) + Alfa*(dup28*(dup29*EL*(-(2*MW2*
     -              (DCONJG(SCha2(Cha4,1))*mUS(Sfe4,1,2)*
     -              mUSC(Sfe4,2,2)*SCha1(Cha5,2))) + MW*
     -              (DCONJG(SCha1(Cha4,2))*Sqrt(2.D0)*mUS(Sfe4,2,2)*
     -              mUSC(Sfe4,2,2)*SCha1(Cha5,2)))) + dup30*EL*(2*MW2*
     -              (DCONJG(SCha2(Cha4,1))*MCha(Cha4)*MCha(Cha5)*
     -              mUS(Sfe4,1,2)*mUSC(Sfe4,2,2)*SCha1(Cha5,2)) - 
     -              MW*(DCONJG(SCha1(Cha4,2))*Sqrt(2.D0)*MCha(Cha4)*
     -              MCha(Cha5)*mUS(Sfe4,2,2)*mUSC(Sfe4,2,2)*
     -              SCha1(Cha5,2)))))/(SW*Sqrt(2.D0))
               
               tmp1 = tmp1 + 1/(8.D0*CW*MW*MW2*Pi*SW2)*(Alfa*(dup28*
     -              dup29*EL*(2*MW*MW2*(DCONJG(SCha2(Cha4,1))*
     -              Sqrt(2.D0)*mUS(Sfe4,1,2)*mUSC(Sfe4,1,2)*
     -              SCha2(Cha5,1)) - 2*MW2*(DCONJG(SCha1(Cha4,2))*
     -              mUS(Sfe4,2,2)*mUSC(Sfe4,1,2)*SCha2(Cha5,1))) + 
     -              dup30*EL*(-(2*MW*MW2*(DCONJG(SCha2(Cha4,1))*
     -              Sqrt(2.D0)*MCha(Cha4)*MCha(Cha5)*mUS(Sfe4,1,2)*
     -              mUSC(Sfe4,1,2)*SCha2(Cha5,1))) + 2*MW2*
     -              (DCONJG(SCha1(Cha4,2))*MCha(Cha4)*MCha(Cha5)*
     -              mUS(Sfe4,2,2)*mUSC(Sfe4,1,2)*SCha2(Cha5,1)))) + 
     -              Alfa*(dup29*EL*Qmass(1)*(-(MW*(Sqrt(2.D0)*
     -              MCha(Cha4)*mUS(Sfe4,1,2)*mUSC(Sfe4,2,2)*
     -              SCha1(Cha4,1)*SCha1(Cha5,2))) + 2*MW2*
     -              (MCha(Cha4)*mUS(Sfe4,1,2)*mUSC(Sfe4,1,2)*
     -              SCha1(Cha4,1)*SCha2(Cha5,1)))))/(SW*Sqrt(2.D0))
               
               tmp2 = 1/(8.D0*CW*MW*MW2*Pi*SW2)*(Alfa*((dup30*EL*
     -              Qmass(1)*(2*MW2*(DCONJG(SCha1(Cha5,1))*
     -              DCONJG(SCha2(Cha4,1))*MCha(Cha4)*mUS(Sfe4,1,2)*
     -              mUSC(Sfe4,1,2)) + MW*(-(DCONJG(SCha1(Cha4,2))*
     -              DCONJG(SCha1(Cha5,1))*Sqrt(2.D0)*MCha(Cha4)*
     -              mUS(Sfe4,2,2)*mUSC(Sfe4,1,2))+Sqrt(2.D0)*MCha(Cha5)*
     -              mUS(Sfe4,1,2)*mUSC(Sfe4,2,2)*SCha1(Cha4,1)*
     -              SCha1(Cha5,2))) + dup29*EL*Qmass(1)*(MW*
     -              (DCONJG(SCha1(Cha4,2))*DCONJG(SCha1(Cha5,1))*
     -              Sqrt(2.D0)*MCha(Cha5)*mUS(Sfe4,2,2)*mUSC(Sfe4,1,2) - 
     -              Sqrt(2.D0)*MCha(Cha4)*mUS(Sfe4,1,2)*mUSC(Sfe4,2,2)*
     -              SCha1(Cha4,1)*SCha1(Cha5,2)) + MW2*(-2*
     -              DCONJG(SCha1(Cha5,1))*DCONJG(SCha2(Cha4,1))*
     -              MCha(Cha5)*mUS(Sfe4,1,2)*mUSC(Sfe4,1,2) + 2*
     -              MCha(Cha4)*mUS(Sfe4,1,2)*mUSC(Sfe4,1,2)*
     -              SCha1(Cha4,1)*SCha2(Cha5,1))))))/(SW*Sqrt(2.D0))
               
               tmp2 = tmp2- Alfa/(4.D0*CW*MW*Pi*SW2)*(dup30*EL*Qmass(1)*
     -              MCha(Cha5)*mUS(Sfe4,1,2)*mUSC(Sfe4,1,2)*
     -              SCha1(Cha4,1)*SCha2(Cha5,1))/(SW*Sqrt(2.D0))

               tmp3 = 1/(8.D0*CW*MW*MW2*Pi*SW2)*(Alfa*((MW*(dup30*EL*
     -              Qmass(1)**2*DCONJG(SCha1(Cha5,1))*Sqrt(2.D0)*
     -              mUS(Sfe4,1,2)*mUSC(Sfe4,1,2)*SCha1(Cha4,1)) + 
     -              dup29*EL*Qmass(1)**2*(2*MW2*(DCONJG(SCha2(Cha4,1))*
     -              mUS(Sfe4,1,2)*mUSC(Sfe4,2,2)*SCha1(Cha5,2)) - MW*
     -              (DCONJG(SCha1(Cha4,2))*Sqrt(2.D0)*mUS(Sfe4,2,2)*
     -              mUSC(Sfe4,2,2)*SCha1(Cha5,2))))) + Alfa*(dup29*EL*
     -              (-MW*(Qmass(1)*Sqrt(2.D0)*MCha(Cha4)*mUS(Sfe4,1,2)*
     -              mUSC(Sfe4,2,2)*SCha1(Cha4,1)*SCha1(Cha5,2)))+dup30*
     -              EL*(MW*(Qmass(1)*Sqrt(2.D0)*MCha(Cha5)*
     -              mUS(Sfe4,1,2)*mUSC(Sfe4,2,2)*SCha1(Cha4,1)*
     -              SCha1(Cha5,2)))))/(SW*Sqrt(2.D0))
               
               tmp3 = tmp3 + 1/(4.D0*CW*MW*Pi*SW2)*(Alfa*(dup29*EL*
     -              Qmass(1)**2*(-(MW*(DCONJG(SCha2(Cha4,1))*Sqrt(2.D0)*
     -              mUS(Sfe4,1,2)*mUSC(Sfe4,1,2)*SCha2(Cha5,1))) + 
     -              DCONJG(SCha1(Cha4,2))*mUS(Sfe4,2,2)*mUSC(Sfe4,1,2)*
     -              SCha2(Cha5,1))) +Alfa*(dup29*EL*Qmass(1)*MCha(Cha4)*
     -              mUS(Sfe4,1,2)*mUSC(Sfe4,1,2)*SCha1(Cha4,1)*
     -              SCha2(Cha5,1) - dup30*EL*Qmass(1)*MCha(Cha5)*
     -              mUS(Sfe4,1,2)*mUSC(Sfe4,1,2)*SCha1(Cha4,1)*
     -              SCha2(Cha5,1)))/(SW*Sqrt(2.D0))
               
               loopNL = loopNL + tmp1*Cval(cc0,iint5(Cha4,Cha5,Sfe4)) + 
     -              tmp2*Cval(cc1,iint5(Cha4,Cha5,Sfe4)) + tmp3*
     -              Cval(cc2,iint5(Cha4,Cha5,Sfe4)) + Alfa/(8.D0*CW*MW*
     -              MW2*Pi*SW2)*(Bval(bb0,iint6(Cha4,Cha5))*(dup29*EL*
     -              (MW*(2*MW2*(DCONJG(SCha2(Cha4,1))*Sqrt(2.D0)*
     -              mUS(Sfe4,1,2)*mUSC(Sfe4,1,2)*SCha2(Cha5,1)) + 
     -              DCONJG(SCha1(Cha4,2))*Sqrt(2.D0)*mUS(Sfe4,2,2)*
     -              mUSC(Sfe4,2,2)*SCha1(Cha5,2))+MW2*(-2*
     -              DCONJG(SCha2(Cha4,1))*mUS(Sfe4,1,2)*mUSC(Sfe4,2,2)*
     -              SCha1(Cha5,2)-2*DCONJG(SCha1(Cha4,2))*mUS(Sfe4,2,2)*
     -              mUSC(Sfe4,1,2)*SCha2(Cha5,1)))))/(SW*Sqrt(2.D0)) + 
     -              Alfa/(4.D0*CW*MW*MW2*Pi*SW2)*
     -              (Cval(cc00,iint5(Cha4,Cha5,Sfe4))*(dup29*EL*(MW*
     -              (-(2*MW2*(DCONJG(SCha2(Cha4,1))*Sqrt(2.D0)*
     -              mUS(Sfe4,1,2)*mUSC(Sfe4,1,2)*SCha2(Cha5,1))) - 
     -              DCONJG(SCha1(Cha4,2))*Sqrt(2.D0)*mUS(Sfe4,2,2)*
     -              mUSC(Sfe4,2,2)*SCha1(Cha5,2)) + MW2*(2*
     -              DCONJG(SCha2(Cha4,1))*mUS(Sfe4,1,2)*mUSC(Sfe4,2,2)*
     -              SCha1(Cha5,2)+2*DCONJG(SCha1(Cha4,2))*mUS(Sfe4,2,2)*
     -              mUSC(Sfe4,1,2)*SCha2(Cha5,1)))))/(SW*Sqrt(2.D0))

               tmp1 = (0D0,0D0)
               tmp2 = (0D0,0D0)
               tmp3 = (0D0,0D0)

               tmp1 = 1/(8.D0*CW*MW*MW2*Pi*SW2)*(Alfa*(dup30*EL*
     -              Qmass(1)*(2*MW2*(DCONJG(SCha1(Cha5,1))*
     -              DCONJG(SCha2(Cha4,1))*MCha(Cha4)*mUS(Sfe4,1,2)*
     -              mUSC(Sfe4,1,2)) - MW*(DCONJG(SCha1(Cha4,2))*
     -              DCONJG(SCha1(Cha5,1))*Sqrt(2.D0)*MCha(Cha4)*
     -              mUS(Sfe4,2,2)*mUSC(Sfe4,1,2)))) + Alfa*(-(MW*(dup29*
     -              EL*DCONJG(SCha1(Cha5,1))*Sqrt(2.D0)*MCha(Cha4)*
     -              MCha(Cha5)*mUS(Sfe4,1,2)*mUSC(Sfe4,1,2)*
     -              SCha1(Cha4,1))) + dup28*(MW*(dup30*EL*
     -              DCONJG(SCha1(Cha5,1))*Sqrt(2.D0)*mUS(Sfe4,1,2)*
     -              mUSC(Sfe4,1,2)*SCha1(Cha4,1)))))/(SW*Sqrt(2.D0))
               
               tmp1 = tmp1 + 1/(8.D0*CW*MW*MW2*Pi*SW2)*(Alfa*(dup29*EL*
     -              Qmass(1)*(-(MW*(Sqrt(2.D0)*MCha(Cha4)*mUS(Sfe4,1,2)*
     -              mUSC(Sfe4,2,2)*SCha1(Cha4,1)*SCha1(Cha5,2)))+ 2*MW2*
     -              (MCha(Cha4)*mUS(Sfe4,1,2)*mUSC(Sfe4,1,2)*
     -              SCha1(Cha4,1)*SCha2(Cha5,1)))))/(SW*Sqrt(2.D0))
               
               tmp2 = 1/(8.D0*CW*MW*MW2*Pi*SW2)*(Alfa*((dup30*EL*
     -              Qmass(1)*(2*MW2*(DCONJG(SCha1(Cha5,1))*
     -              DCONJG(SCha2(Cha4,1))*MCha(Cha4)*mUS(Sfe4,1,2)*
     -              mUSC(Sfe4,1,2)) + MW*(-(DCONJG(SCha1(Cha4,2))*
     -              DCONJG(SCha1(Cha5,1))*Sqrt(2.D0)*MCha(Cha4)*
     -              mUS(Sfe4,2,2)*mUSC(Sfe4,1,2))+Sqrt(2.D0)*MCha(Cha5)*
     -              mUS(Sfe4,1,2)*mUSC(Sfe4,2,2)*SCha1(Cha4,1)*
     -              SCha1(Cha5,2))) + dup29*EL*Qmass(1)*(MW*
     -              (DCONJG(SCha1(Cha4,2))*DCONJG(SCha1(Cha5,1))*
     -              Sqrt(2.D0)*MCha(Cha5)*mUS(Sfe4,2,2)*mUSC(Sfe4,1,2) - 
     -              Sqrt(2.D0)*MCha(Cha4)*mUS(Sfe4,1,2)*mUSC(Sfe4,2,2)*
     -              SCha1(Cha4,1)*SCha1(Cha5,2)) + MW2*(-2*
     -              DCONJG(SCha1(Cha5,1))*DCONJG(SCha2(Cha4,1))*
     -              MCha(Cha5)*mUS(Sfe4,1,2)*mUSC(Sfe4,1,2) + 2*
     -              MCha(Cha4)*mUS(Sfe4,1,2)*mUSC(Sfe4,1,2)*
     -              SCha1(Cha4,1)*SCha2(Cha5,1))))))/(SW*Sqrt(2.D0))
               
               tmp2 = tmp2- Alfa/(4.D0*CW*MW*Pi*SW2)*(dup30*EL*Qmass(1)*
     -              MCha(Cha5)*mUS(Sfe4,1,2)*mUSC(Sfe4,1,2)*
     -              SCha1(Cha4,1)*SCha2(Cha5,1))/(SW*Sqrt(2.D0))

               tmp3 = 1/(8.D0*CW*MW*MW2*Pi*SW2)*(Alfa*((-1)*(MW*(dup30*
     -              EL*Qmass2(1)*DCONJG(SCha1(Cha5,1))*Sqrt(2.D0)*
     -              mUS(Sfe4,1,2)*mUSC(Sfe4,1,2)*SCha1(Cha4,1)) + dup29*
     -              EL*Qmass2(1)*(2*MW2*(DCONJG(SCha2(Cha4,1))*
     -              mUS(Sfe4,1,2)*mUSC(Sfe4,2,2)*SCha1(Cha5,2)) - MW*
     -              (DCONJG(SCha1(Cha4,2))*Sqrt(2.D0)*mUS(Sfe4,2,2)*
     -              mUSC(Sfe4,2,2)*SCha1(Cha5,2))))) + Alfa*(dup29*EL*
     -              (Qmass(1)*(-(2*MW2*(DCONJG(SCha1(Cha5,1))*
     -              DCONJG(SCha2(Cha4,1))*MCha(Cha5)*mUS(Sfe4,1,2)*
     -              mUSC(Sfe4,1,2))) + MW*(DCONJG(SCha1(Cha4,2))*
     -              DCONJG(SCha1(Cha5,1))*Sqrt(2.D0)*MCha(Cha5)*
     -              mUS(Sfe4,2,2)*mUSC(Sfe4,1,2)))) +dup30*EL*(Qmass(1)*
     -              (2*MW2*(DCONJG(SCha1(Cha5,1))*DCONJG(SCha2(Cha4,1))*
     -              MCha(Cha4)*mUS(Sfe4,1,2)*mUSC(Sfe4,1,2)) - MW*
     -              (DCONJG(SCha1(Cha4,2))*DCONJG(SCha1(Cha5,1))*
     -              Sqrt(2.D0)*MCha(Cha4)*mUS(Sfe4,2,2)*
     -              mUSC(Sfe4,1,2))))))/(SW*Sqrt(2.D0))
               
               tmp3 = tmp3 + 1/(4.D0*CW*MW*Pi*SW2)*(Alfa*((-1)*dup29*EL*
     -              Qmass2(1)*(-(MW*(DCONJG(SCha2(Cha4,1))*Sqrt(2.D0)*
     -              mUS(Sfe4,1,2)*mUSC(Sfe4,1,2)*SCha2(Cha5,1))) + 
     -              DCONJG(SCha1(Cha4,2))*mUS(Sfe4,2,2)*mUSC(Sfe4,1,2)*
     -              SCha2(Cha5,1))))/(SW*Sqrt(2.D0))

               loopNR = loopNR + tmp1*Cval(cc0,iint5(Cha4,Cha5,Sfe4)) + 
     -              tmp2*Cval(cc1,iint5(Cha4,Cha5,Sfe4)) + tmp3*
     -              Cval(cc2,iint5(Cha4,Cha5,Sfe4)) + Alfa/(8.D0*CW*MW*
     -              MW2*Pi*SW2)*(Bval(bb0,iint6(Cha4,Cha5))*(MW*(dup30*
     -              EL*DCONJG(SCha1(Cha5,1))*Sqrt(2.D0)*mUS(Sfe4,1,2)*
     -              mUSC(Sfe4,1,2)*SCha1(Cha4,1))))/(SW*Sqrt(2.D0)) + 
     -              Alfa/(4.D0*CW*MW*MW2*Pi*SW2)*
     -              (Cval(cc00,iint5(Cha4,Cha5,Sfe4))*(-(MW*(dup30*EL*
     -              DCONJG(SCha1(Cha5,1))*Sqrt(2.D0)*mUS(Sfe4,1,2)*
     -              mUSC(Sfe4,1,2)*SCha1(Cha4,1)))))/(SW*Sqrt(2.D0))
               
*               call loop3_Cha4_Cha5_Sfe4
            enddo
         enddo
      enddo

      do Sfe5 = 1,2
         do Sfe4 = 1,2
            do Cha4 = 1,2
               dup27 = QgL(2)*mUS(Sfe5,1,2)*mUSC(Sfe4,1,2) + 
     -              QgR(2)*mUS(Sfe5,2,2)*mUSC(Sfe4,2,2)

               loopPL = loopPL + Alfa/(4.D0*MW*MW2*Pi*SW2)*(dup27*EL*
     -              Cval(cc00,iint4(Cha4,Sfe4,Sfe5))*((MW*(2*MW2*
     -              (DCONJG(SCha2(Cha4,1))*Sqrt(2.D0)*mUS(Sfe4,1,2)*
     -              mUSC(Sfe5,1,2)*SCha2(Cha4,1))+DCONJG(SCha1(Cha4,2))*
     -              Sqrt(2.D0)*mUS(Sfe4,2,2)*mUSC(Sfe5,2,2)*
     -              SCha1(Cha4,2))+MW2*(-2*DCONJG(SCha2(Cha4,1))*
     -              mUS(Sfe4,1,2)*mUSC(Sfe5,2,2)*SCha1(Cha4,2) - 
     -              2*DCONJG(SCha1(Cha4,2))*mUS(Sfe4,2,2)*
     -              mUSC(Sfe5,1,2)*SCha2(Cha4,1)))))/Sqrt(2.D0)


               loopPR = loopPR + Alfa/(4.D0*MW*MW2*Pi*SW2)*
     -              (dup27*EL*Cval(cc00,iint4(Cha4,Sfe4,Sfe5))*
     -              (MW*(DCONJG(SCha1(Cha4,1))*Sqrt(2.D0)*mUS(Sfe4,1,2)*
     -              mUSC(Sfe5,1,2)*SCha1(Cha4,1))))/Sqrt(2.D0)
*               call loop2_Cha4_Sfe4_Sfe5
            enddo
         enddo
      enddo

      do Sfe4 = 1,2
         do Neu5 = 1,4
            do Neu4 = 1,4
               dup12 = 3*CW*(mUS(Sfe4,1,1)*SNeu1(Neu4,1)) - 
     -              2*MW*(SW*mUS(Sfe4,2,1)*SNeu3(1)*ZNeu(Neu4,1))
               dup13 = -(3*CW*(mUS(Sfe4,1,1)*SNeu1(Neu4,1))) + 
     -              2*MW*(SW*mUS(Sfe4,2,1)*SNeu3(1)*ZNeu(Neu4,1))
               dup14 = 3*CW*(SNeu2(1)*ZNeu(Neu5,2)) + SW*ZNeu(Neu5,1)
               dup15 = 3*CW*(SNeu2(1)*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
               dup16 = ZNeu(Neu5,3)*ZNeuC(Neu4,3) - ZNeu(Neu5,4)*
     -              ZNeuC(Neu4,4)
               dup17 = -(ZNeu(Neu5,3)*ZNeuC(Neu4,3)) + 
     -              ZNeu(Neu5,4)*ZNeuC(Neu4,4)
               dup18 = ZNeu(Neu4,3)*ZNeuC(Neu5,3) - 
     -              ZNeu(Neu4,4)*ZNeuC(Neu5,4)
               dup19 = MW*(dup15*mUS(Sfe4,1,1)) + 
     -              3*CW*(DCONJG(SNeu1(Neu4,1))*mUS(Sfe4,2,1))
               dup20 = 3*CW*(DCONJG(SNeu1(Neu5,1))*mUSC(Sfe4,1,1)) - 
     -              2*MW*(SW*mUSC(Sfe4,2,1)*SNeu3(1)*ZNeuC(Neu5,1))
               dup21 = -(3*CW*(DCONJG(SNeu1(Neu5,1))*mUSC(Sfe4,1,1))) + 
     -              2*MW*(SW*mUSC(Sfe4,2,1)*SNeu3(1)*ZNeuC(Neu5,1))
               dup22 = dup17*dup19*dup21 + dup13*dup18*(MW*(dup14*
     -              mUSC(Sfe4,1,1))+3*CW*(mUSC(Sfe4,2,1)*SNeu1(Neu5,1)))
               dup24 = dup12*dup16*dup20 + dup18*dup19*(MW*(dup14*
     -              mUSC(Sfe4,1,1))+3*CW*(mUSC(Sfe4,2,1)*SNeu1(Neu5,1)))


               dup26 = dup18*dup19*(MW*(dup14*mUSC(Sfe4,1,1)) + 
     -              3*CW*(mUSC(Sfe4,2,1)*SNeu1(Neu5,1)))

               loopPL = loopPL + Alfa/(144.D0*CW*CW2*MW2*Pi*SW2)*(dup26*
     -              EL*Bval(bb0,iint3(Neu4,Neu5)))/SW - Alfa/(72.D0*CW*
     -              CW2*MW2*Pi*SW2)*(dup26*EL*
     -              Cval(cc00,iint2(Neu4,Neu5,Sfe4)))/SW + 
     -              Cval(cc1,iint2(Neu4,Neu5,Sfe4))*(Alfa/(144.D0*CW*
     -              CW2*MW2*Pi*SW2)*(dup22*EL*Qmass(1)*MNeu(Neu4))/SW - 
     -              Alfa/(144.D0*CW*CW2*MW2*Pi*SW2)*(EL*Qmass(1)*(dup18*
     -              dup19*dup21+dup13*dup17*(MW*(dup14*mUSC(Sfe4,1,1)) + 
     -              3*CW*(mUSC(Sfe4,2,1)*SNeu1(Neu5,1))))*MNeu(Neu5))/
     -              SW)+ Cval(cc2,iint2(Neu4,Neu5,Sfe4))*(-Alfa/(144.D0*
     -              CW*CW2*MW2*Pi*SW2)*(dup24*EL*Qmass(1)**2)/SW + Alfa/
     -              (144.D0*CW*CW2*MW2*Pi*SW2)*(EL*Qmass(1)*(dup13*
     -              dup18*(MW*(dup14*mUSC(Sfe4,1,1)) + 3*CW*
     -              (mUSC(Sfe4,2,1)*SNeu1(Neu5,1))))*MNeu(Neu4))/SW - 
     -              Alfa/(144.D0*CW*CW2*MW2*Pi*SW2)*(EL*Qmass(1)*(dup13*
     -              dup17*(MW*(dup14*mUSC(Sfe4,1,1)) + 3*CW*
     -              (mUSC(Sfe4,2,1)*SNeu1(Neu5,1))))*MNeu(Neu5))/SW)
               
               loopPL = loopPL + Cval(cc0,iint2(Neu4,Neu5,Sfe4))*(-Alfa/
     -              (144.D0*CW*CW2*MW2*Pi*SW2)*(dup26*EL*(Qmass(1)**2 - 
     -              mass2SF(Sfe4,1)))/SW + Alfa/(144.D0*CW*CW2*MW2*Pi*
     -              SW2)*(dup22*EL*Qmass(1)*MNeu(Neu4))/SW - Alfa/
     -              (144.D0*CW*CW2*MW2*Pi*SW2)*(EL*(dup17*dup19*(MW*
     -              (dup14*mUSC(Sfe4,1,1)) + 3*CW*(mUSC(Sfe4,2,1)*
     -              SNeu1(Neu5,1))))*MNeu(Neu4)*MNeu(Neu5))/SW)


               dup26 = dup13*dup17*dup21
               
               loopPR = loopPR + Alfa/(144.D0*CW*CW2*MW2*Pi*SW2)*(dup26*
     -              EL*Bval(bb0,iint3(Neu4,Neu5)))/SW - Alfa/(72.D0*CW* 
     -              CW2*MW2*Pi*SW2)*(dup26*EL*
     -              Cval(cc00,iint2(Neu4,Neu5,Sfe4)))/SW + 
     -              Cval(cc1,iint2(Neu4,Neu5,Sfe4))*(Alfa/(144.D0*CW*
     -              CW2*MW2*Pi*SW2)*(dup22*EL*Qmass(1)*MNeu(Neu4))/SW - 
     -              Alfa/(144.D0*CW*CW2*MW2*Pi*SW2)*(EL*Qmass(1)*(dup18*
     -              dup19*dup21+dup13*dup17*(MW*(dup14*mUSC(Sfe4,1,1)) + 
     -              3*CW*(mUSC(Sfe4,2,1)*SNeu1(Neu5,1))))*MNeu(Neu5))/
     -              SW) + Cval(cc2,iint2(Neu4,Neu5,Sfe4))*(-Alfa/
     -              (144.D0*CW*CW2*MW2*Pi*SW2)*((-1)*dup24*EL*
     -              Qmass2(1))/SW + Alfa/(144.D0*CW*CW2*MW2*Pi*SW2)*
     -              (EL*Qmass(1)*(dup17*dup19*dup21)*MNeu(Neu4))/SW - 
     -              Alfa/(144.D0*CW*CW2*MW2*Pi*SW2)*(EL*Qmass(1)*(dup18*
     -              dup19*dup21)*MNeu(Neu5))/SW)
               
               loopPR = loopPR + Cval(cc0,iint2(Neu4,Neu5,Sfe4))*(-Alfa/
     -              (144.D0*CW*CW2*MW2*Pi*SW2)*(dup26*EL*(Qmass2(1) - 
     -              mass2SF(Sfe4,1)))/SW + Alfa/(144.D0*CW*CW2*MW2*Pi*
     -              SW2)*(dup22*EL*Qmass(1)*MNeu(Neu4))/SW - Alfa/
     -              (144.D0*CW*CW2*MW2*Pi*SW2)*(EL*(dup12*dup18*dup20)*
     -              MNeu(Neu4)*MNeu(Neu5))/SW)
*               call loop1_Neu4_Neu5_Sfe4
            enddo
         enddo
      enddo

      do Sfe5 = 1,2
         do Sfe4 = 1,2
            do Neu4 = 1,4
               dup1 = QgL(1)*mUS(Sfe5,1,1)*mUSC(Sfe4,1,1) + 
     -              QgR(1)*mUS(Sfe5,2,1)*mUSC(Sfe4,2,1)
               dup2 = 3*CW*(mUS(Sfe4,1,1)*SNeu1(Neu4,1)) - 
     -              2*MW*(SW*mUS(Sfe4,2,1)*SNeu3(1)*ZNeu(Neu4,1))
               dup3 = -(3*CW*(mUS(Sfe4,1,1)*SNeu1(Neu4,1))) + 
     -              2*MW*(SW*mUS(Sfe4,2,1)*SNeu3(1)*ZNeu(Neu4,1))
               dup4 = 3*CW*(SNeu2(1)*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)
               dup5 = 3*CW*(SNeu2(1)*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
               dup6 = MW*(dup5*mUS(Sfe4,1,1)) + 
     -              3*CW*(DCONJG(SNeu1(Neu4,1))*mUS(Sfe4,2,1))
               dup7 = 3*CW*(DCONJG(SNeu1(Neu4,1))*mUSC(Sfe5,1,1)) - 
     -              2*MW*(SW*mUSC(Sfe5,2,1)*SNeu3(1)*ZNeuC(Neu4,1))
               dup8 = -(3*CW*(DCONJG(SNeu1(Neu4,1))*mUSC(Sfe5,1,1))) + 
     -              2*MW*(SW*mUSC(Sfe5,2,1)*SNeu3(1)*ZNeuC(Neu4,1))


               loopPL = loopPL + Alfa/(36.D0*CW2*MW2*Pi*SW2)*(dup1*EL*
     -              Cval(cc00,iint1(Neu4,Sfe4,Sfe5))*(dup6*(MW*(dup4*
     -              mUSC(Sfe5,1,1)) + 3*CW*(mUSC(Sfe5,2,1)*
     -              SNeu1(Neu4,1))))) 


               loopPR = loopPR + Alfa/(36.D0*CW2*MW2*Pi*SW2)*(dup1*EL*
     -              Cval(cc00,iint1(Neu4,Sfe4,Sfe5))*(dup2*dup7)) 
*               call loop1_Neu4_Sfe4_Sfe5
            enddo
         enddo
      enddo


** Set individual parts to zero for debugging
*      loopNL = (0D0,0D0)
*      loopPL = (0D0,0D0)
*      loopNR = (0D0,0D0)
*      loopPR = (0D0,0D0)


      if (qqVcur .eq. 2) then
         loopNL = -loopNL
         loopNR = -loopNR
      end if

         
      qqZlp(1) = qqZlp(1) + loopNL + loopPL
      qqZlp(2) = qqZlp(2) + loopNR + loopPR     


      end 


********************************************************************************
********************************************************************************

*** Subroutine calculating the SUSY loop contributions to the qqW vertex.  NOTE 
*** that we can split the loop diagrams into SUSY and SM parts because we're 
*** pretending that external particles have no mass and thus do not couple to 
*** the higgs

      subroutine qqWsusy(qqVcur,qqWlp,q12,q22,q32,QgL,QgR,Qmass,Qmass2,
     -     SNeu1,SNeu2,SNeu3,SCha1,SCha2,cnwQ,massSF,mass2SF,mUS,mUSC)

      implicit none


** Effective coupling from loops for qqZ coupling Labelling (L/R, gen)
      double complex qqWlp(2)


** SM couplings and masses for particle number (i) for the qqV current
** being considered
      double precision QgL(2), QgR(2), Qmass(2), Qmass2(2)

** SUSY couplings and masses for the qqV current being considered
** Labelling: final entry = particle number (i)
**            first entry = neutralino/chargino index
      double complex SNeu1(4,2), SCha1(2,2), SCha2(2,2)
      double precision SNeu2(2), SNeu3(2)
** chargino-neutralino-W coupling 
** labelling: (cnw, chargino, neutralino)
      double complex cnwQ(4,2,4)

** sfermion masses and mixing matrices
** Labelling: sfe1, (sfe2), particle number i
      double precision massSF(2,2), mass2SF(2,2)
      double complex mUS(2,2,2), mUSC(2,2,2)

** do loop variables
      integer sfe1, sfe2, neu1, neu2, cha1, cha2


** Value of separate loop contributions
      double complex loopNL, loopNR, loopPL, loopPR

** Integral labels
      integer iint1(2,4), iint2(2,4,2), iint3(2,4,2), iint4(4,2,2)

** EM coupling
      double precision Alfa, Alfa2


** do loop variables
      integer Sfe4, Sfe5, Neu4, Cha4

** dummy parameters
      double complex dup1, dup2, dup3, dup4, dup5, dup6, dup7 
      double complex dup8, dup9, dup10, dup11, dup12, dup13, dup14
      double complex dup15, dup16, dup17, dup18, dup19, dup20, dup21
      double complex tmp1, tmp2


** Squared momenta of the quarks and vector boson
      double precision q12, q22, q32


** Label for the qqV current
      integer qqVcur

** Debugging flag
      integer debug


** Parameter and looptools declarations
#include "looptools.h"
#include "mssm.inc"
#include "global.inc"


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** Clearing loop values
      loopNL = (0D0,0D0)
      loopPL = (0D0,0D0)
      loopNR = (0D0,0D0)
      loopPR = (0D0,0D0)


** Setting the integrals
      do Neu4 = 1,4
         do Cha4 = 1,2
            iint1(Cha4,Neu4) = Bget(q32,MCha2(Cha4),MNeu2(Neu4))
         enddo
      enddo
      do Sfe4 = 1,2
         do Neu4 = 1,4
            do Cha4 = 1,2
               iint2(Cha4,Neu4,Sfe4) = Cget(q32,Qmass2(1),Qmass2(2),
     -              MCha2(Cha4),MNeu2(Neu4),mass2SF(Sfe4,1))
               iint3(Cha4,Neu4,Sfe4) = Cget(q32,Qmass2(2),Qmass2(1),
     -              MCha2(Cha4),MNeu2(Neu4),mass2SF(Sfe4,2))
            enddo
         enddo
      enddo
      do Sfe5 = 1,2
         do Sfe4 = 1,2
            do Neu4 = 1,4
               iint4(Neu4,Sfe4,Sfe5) = Cget(Qmass2(1),q32,Qmass2(2),
     -              MNeu2(Neu4),mass2SF(Sfe4,1),mass2SF(Sfe5,2))
            enddo
         enddo
      enddo


** The L/R components
      do Sfe4 = 1,2
         do Neu4 = 1,4
            do Cha4 = 1,2

               dup1 = -(cnwQ(2,Cha4,Neu4)*MCha(Cha4)) + 
     -              cnwQ(1,Cha4,Neu4)*MNeu(Neu4)
               dup2 = cnwQ(1,Cha4,Neu4)*MCha(Cha4) - 
     -              cnwQ(2,Cha4,Neu4)*MNeu(Neu4)
               dup3 = -(3*CW*(mUS(Sfe4,1,1)*SNeu1(Neu4,1))) + 
     -              2*MW*(SW*mUS(Sfe4,2,1)*SNeu3(1)*ZNeu(Neu4,1))
               dup4 = 3*CW*(SNeu2(1)*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
               dup5 = -(Qmass2(2)*cnwQ(1,Cha4,Neu4)) + 
     -              cnwQ(1,Cha4,Neu4)*mass2SF(Sfe4,1) - 
     -              cnwQ(2,Cha4,Neu4)*MCha(Cha4)*MNeu(Neu4)
               dup6 = MW*(dup4*mUS(Sfe4,1,1)) + 
     -              3*CW*(DCONJG(SNeu1(Neu4,1))*mUS(Sfe4,2,1))
               dup7 = -(cnwQ(2,Cha4,Neu4)*MCha(Cha4)) + 
     -              cnwQ(1,Cha4,Neu4)*MNeu(Neu4)
               dup8 = cnwQ(1,Cha4,Neu4)*MCha(Cha4) - 
     -              cnwQ(2,Cha4,Neu4)*MNeu(Neu4)
               dup9 = 3*CW*(SNeu2(2)*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)
               dup10 = Qmass2(1)*cnwQ(2,Cha4,Neu4) - 
     -              cnwQ(2,Cha4,Neu4)*mass2SF(Sfe4,2) + 
     -              cnwQ(1,Cha4,Neu4)*MCha(Cha4)*MNeu(Neu4)
               dup11 = -(3*CW*(DCONJG(SNeu1(Neu4,2))*mUSC(Sfe4,1,2))) + 
     -              2*MW*(SW*mUSC(Sfe4,2,2)*SNeu3(2)*ZNeuC(Neu4,1))

               
               loopNL = loopNL + Alfa/(24.D0*CW*MW*MW2*Pi*SW2)*
     -              (Bval(bb0,iint1(Cha4,Neu4))*(dup6*EL*(MW*
     -              (cnwQ(1,Cha4,Neu4)*Sqrt(2.D0)*mUSC(Sfe4,2,1)*
     -              SCha1(Cha4,1)) - 2*MW2*(cnwQ(1,Cha4,Neu4)*
     -              mUSC(Sfe4,1,1)*SCha2(Cha4,2)))))/(SW*Sqrt(2.D0)) + 
     -              Alfa/(12.D0*CW*MW*MW2*Pi*SW2)*
     -              (Cval(cc00,iint2(Cha4,Neu4,Sfe4))*(dup6*EL*(-(MW*
     -              (cnwQ(1,Cha4,Neu4)*Sqrt(2.D0)*mUSC(Sfe4,2,1)*
     -              SCha1(Cha4,1))) + 2*MW2*(cnwQ(1,Cha4,Neu4)*
     -              mUSC(Sfe4,1,1)*SCha2(Cha4,2)))))/(SW*Sqrt(2.D0))
               
               loopNL = loopNL + Cval(cc1,iint2(Cha4,Neu4,Sfe4))*(-Alfa/
     -              (24.D0*CW*MW2*Pi*SW2)*(dup2*dup6*EL*Qmass(2)*
     -              DCONJG(SCha1(Cha4,2))*mUSC(Sfe4,1,1))/SW - Alfa/
     -              (24.D0*CW*MW2*Pi*SW2)*(dup1*dup3*EL*Qmass(1)*
     -              mUSC(Sfe4,2,1)*SCha1(Cha4,1))/SW - Alfa/(24.D0*CW*
     -              MW2*Pi*SW2)*(dup6*EL*(-Qmass2(1) + Qmass2(2))*
     -              cnwQ(1,Cha4,Neu4)*mUSC(Sfe4,2,1)*SCha1(Cha4,1))/SW + 
     -              Alfa/(12.D0*CW*MW*Pi*SW2)*(dup1*dup3*EL*Qmass(1)*
     -              mUSC(Sfe4,1,1)*SCha2(Cha4,2))/(SW*Sqrt(2.D0)) + 
     -              Alfa/(12.D0*CW*MW*Pi*SW2)*(dup6*EL*(-Qmass2(1) + 
     -              Qmass2(2))*cnwQ(1,Cha4,Neu4)*mUSC(Sfe4,1,1)*
     -              SCha2(Cha4,2))/(SW*Sqrt(2.D0)))
               
               loopNL = loopNL + 1/(24.D0*CW*MW*MW2*Pi*SW2)*
     -              (Cval(cc2,iint2(Cha4,Neu4,Sfe4))*(Alfa*(-(MW*(dup2*
     -              dup6*EL*Qmass(2)*DCONJG(SCha1(Cha4,2))*Sqrt(2.D0)*
     -              mUSC(Sfe4,1,1))) + dup3*EL*(-(MW*(Qmass(1)*Qmass(2)*
     -              cnwQ(2,Cha4,Neu4)*DCONJG(SCha1(Cha4,2))*Sqrt(2.D0)*
     -              mUSC(Sfe4,1,1)))) + dup6*EL*(Qmass2(2)*(-(MW*
     -              (cnwQ(1,Cha4,Neu4)*Sqrt(2.D0)*mUSC(Sfe4,2,1)*
     -              SCha1(Cha4,1))) + 2*MW2*(cnwQ(1,Cha4,Neu4)*
     -              mUSC(Sfe4,1,1)*SCha2(Cha4,2)))))))/(SW*Sqrt(2.D0))
               
               loopNL = loopNL + Cval(cc0,iint2(Cha4,Neu4,Sfe4))*(-Alfa/
     -              (24.D0*CW*MW2*Pi*SW2)*(dup6*EL*Qmass(2)*
     -              cnwQ(1,Cha4,Neu4)*DCONJG(SCha1(Cha4,2))*MCha(Cha4)*
     -              mUSC(Sfe4,1,1))/SW + Alfa/(24.D0*CW*MW2*Pi*SW2)*
     -              (dup5*dup6*EL*mUSC(Sfe4,2,1)*SCha1(Cha4,1))/SW + 
     -              Alfa/(24.D0*CW*MW2*Pi*SW2)*(dup3*EL*Qmass(1)*
     -              cnwQ(2,Cha4,Neu4)*MCha(Cha4)*mUSC(Sfe4,2,1)*
     -              SCha1(Cha4,1))/SW - Alfa/(12.D0*CW*MW*Pi*SW2)*(dup5*
     -              dup6*EL*mUSC(Sfe4,1,1)*SCha2(Cha4,2))/(SW*
     -              Sqrt(2.D0)) - Alfa/(12.D0*CW*MW*Pi*SW2)*(dup3*EL*
     -              Qmass(1)*cnwQ(2,Cha4,Neu4)*MCha(Cha4)*
     -              mUSC(Sfe4,1,1)*SCha2(Cha4,2))/(SW*Sqrt(2.D0)))

               
               loopNR = loopNR + Alfa/(24.D0*CW*MW*MW2*Pi*SW2)*
     -              (Bval(bb0,iint1(Cha4,Neu4))*(-(MW*(dup3*EL*
     -              cnwQ(2,Cha4,Neu4)*DCONJG(SCha1(Cha4,2))*Sqrt(2.D0)*
     -              mUSC(Sfe4,1,1)))))/(SW*Sqrt(2.D0)) + Alfa/(12.D0*CW*
     -              MW*MW2*Pi*SW2)*(Cval(cc00,iint2(Cha4,Neu4,Sfe4))*
     -              (MW*(dup3*EL*cnwQ(2,Cha4,Neu4)*
     -              DCONJG(SCha1(Cha4,2))*Sqrt(2.D0)*mUSC(Sfe4,1,1))))/
     -              (SW*Sqrt(2.D0))
               
               loopNR = loopNR + Cval(cc1,iint2(Cha4,Neu4,Sfe4))*(-Alfa/
     -              (24.D0*CW*MW2*Pi*SW2)*(dup2*dup6*EL*Qmass(1)*
     -              DCONJG(SCha1(Cha4,2))*mUSC(Sfe4,1,1))/SW + Alfa/
     -              (24.D0*CW*MW2*Pi*SW2)*(dup3*EL*(-Qmass2(1) + 
     -              Qmass2(2))*cnwQ(2,Cha4,Neu4)*DCONJG(SCha1(Cha4,2))*
     -              mUSC(Sfe4,1,1))/SW-Alfa/(24.D0*CW*MW2*Pi*SW2)*(dup1*
     -              dup3*EL*Qmass(2)*mUSC(Sfe4,2,1)*SCha1(Cha4,1))/
     -              SW + Alfa/(12.D0*CW*MW*Pi*SW2)*(dup1*dup3*EL*
     -              Qmass(2)*mUSC(Sfe4,1,1)*SCha2(Cha4,2))/(SW*
     -              Sqrt(2.D0)))
               
               loopNR = loopNR + 1/(24.D0*CW*MW*MW2*Pi*SW2)*
     -              (Cval(cc2,iint2(Cha4,Neu4,Sfe4))*(Alfa*(dup3*EL*(MW*
     -              (Qmass2(2)*cnwQ(2,Cha4,Neu4)*DCONJG(SCha1(Cha4,2))*
     -              Sqrt(2.D0)*mUSC(Sfe4,1,1))) + dup1*dup3*EL*Qmass(2)*
     -              (-(MW*(Sqrt(2.D0)*mUSC(Sfe4,2,1)*SCha1(Cha4,1))) + 
     -              2*MW2*(mUSC(Sfe4,1,1)*SCha2(Cha4,2))) + dup6*EL*
     -              (Qmass(1)*Qmass(2)*(MW*(cnwQ(1,Cha4,Neu4)*
     -              Sqrt(2.D0)*mUSC(Sfe4,2,1)*SCha1(Cha4,1)) - 2*MW2*
     -              (cnwQ(1,Cha4,Neu4)*mUSC(Sfe4,1,1)*
     -              SCha2(Cha4,2)))))))/(SW*Sqrt(2.D0))
               
               loopNR = loopNR + Cval(cc0,iint2(Cha4,Neu4,Sfe4))*(-Alfa/
     -              (24.D0*CW*MW2*Pi*SW2)*(dup6*EL*Qmass(1)*
     -              cnwQ(1,Cha4,Neu4)*DCONJG(SCha1(Cha4,2))*MCha(Cha4)*
     -              mUSC(Sfe4,1,1))/SW + Alfa/(24.D0*CW*MW2*Pi*SW2)*
     -              (dup3*EL*DCONJG(SCha1(Cha4,2))*(Qmass2(2)*
     -              cnwQ(2,Cha4,Neu4)-cnwQ(2,Cha4,Neu4)*mass2SF(Sfe4,1)+ 
     -              cnwQ(1,Cha4,Neu4)*MCha(Cha4)*MNeu(Neu4))*
     -              mUSC(Sfe4,1,1))/SW + Alfa/(24.D0*CW*MW2*Pi*SW2)*
     -              (dup3*EL*Qmass(2)*cnwQ(2,Cha4,Neu4)*MCha(Cha4)*
     -              mUSC(Sfe4,2,1)*SCha1(Cha4,1))/SW -Alfa/(12.D0*CW*MW*
     -              Pi*SW2)*(dup3*EL*Qmass(2)*cnwQ(2,Cha4,Neu4)*
     -              MCha(Cha4)*mUSC(Sfe4,1,1)*SCha2(Cha4,2))/(SW*
     -              Sqrt(2.D0)))
*               call loop1_Cha4_Neu4_Sfe4


               tmp1 = -Alfa/(12.D0*CW*MW*Pi*SW2)*(EL*Qmass2(1)*
     -              cnwQ(2,Cha4,Neu4)*DCONJG(SCha2(Cha4,1))*(MW*(dup9*
     -              mUSC(Sfe4,1,2)) + 3*CW*(mUSC(Sfe4,2,2)*
     -              SNeu1(Neu4,2)))*mUS(Sfe4,1,2))/(SW*Sqrt(2.D0)) + 
     -              Alfa/(24.D0*CW*MW2*Pi*SW2)*(EL*Qmass2(1)*
     -              cnwQ(2,Cha4,Neu4)*DCONJG(SCha1(Cha4,2))*(MW*(dup9*
     -              mUSC(Sfe4,1,2)) + 3*CW*(mUSC(Sfe4,2,2)*
     -              SNeu1(Neu4,2)))*mUS(Sfe4,2,2))/SW + Alfa/(24.D0*CW*
     -              MW2*Pi*SW2)*(dup11*EL*Qmass(1)*Qmass(2)*
     -              cnwQ(1,Cha4,Neu4)*mUS(Sfe4,1,2)*SCha1(Cha4,1))/SW - 
     -              Alfa/(24.D0*CW*MW2*Pi*SW2)*(dup7*EL*Qmass(1)*(MW*
     -              (dup9*mUSC(Sfe4,1,2)) + 3*CW*(mUSC(Sfe4,2,2)*
     -              SNeu1(Neu4,2)))*mUS(Sfe4,1,2)*SCha1(Cha4,1))/SW
               
               tmp2 = -(Alfa/(12.D0*CW*MW*Pi*SW2)*(dup10*EL*
     -              DCONJG(SCha2(Cha4,1))*(MW*(dup9*mUSC(Sfe4,1,2)) + 
     -              3*CW*(mUSC(Sfe4,2,2)*SNeu1(Neu4,2)))*mUS(Sfe4,1,2))/
     -              (SW*Sqrt(2.D0))) + Alfa/(12.D0*CW*MW*Pi*SW2)*
     -              (dup11*EL*Qmass(2)*cnwQ(1,Cha4,Neu4)*
     -              DCONJG(SCha2(Cha4,1))*MCha(Cha4)*mUS(Sfe4,1,2))/
     -              (SW*Sqrt(2.D0)) - Alfa/(24.D0*CW*MW2*Pi*SW2)*
     -              (dup11*EL*Qmass(2)*cnwQ(1,Cha4,Neu4)*
     -              DCONJG(SCha1(Cha4,2))*MCha(Cha4)*mUS(Sfe4,2,2))/SW 
               
               tmp2 = tmp2 + Alfa/(24.D0*CW*MW2*Pi*SW2)*(dup10*EL*
     -              DCONJG(SCha1(Cha4,2))*(MW*(dup9*mUSC(Sfe4,1,2)) + 
     -              3*CW*(mUSC(Sfe4,2,2)*SNeu1(Neu4,2)))*mUS(Sfe4,2,2))/
     -              SW + Alfa/(24.D0*CW*MW2*Pi*SW2)*(EL*Qmass(1)*
     -              cnwQ(2,Cha4,Neu4)*(MW*(dup9*mUSC(Sfe4,1,2)) + 3*CW*
     -              (mUSC(Sfe4,2,2)*SNeu1(Neu4,2)))*MCha(Cha4)*
     -              mUS(Sfe4,1,2)*SCha1(Cha4,1))/SW
               
               loopNL = loopNL + tmp2*Cval(cc0,iint3(Cha4,Neu4,Sfe4)) + 
     -              tmp1*Cval(cc2,iint3(Cha4,Neu4,Sfe4)) + 
     -              Bval(bb0,iint1(Cha4,Neu4))*(Alfa/(12.D0*CW*MW*Pi*
     -              SW2)*(EL*cnwQ(2,Cha4,Neu4)*DCONJG(SCha2(Cha4,1))*
     -              (MW*(dup9*mUSC(Sfe4,1,2)) + 3*CW*(mUSC(Sfe4,2,2)*
     -              SNeu1(Neu4,2)))*mUS(Sfe4,1,2))/(SW*Sqrt(2.D0)) - 
     -              Alfa/(24.D0*CW*MW2*Pi*SW2)*(EL*cnwQ(2,Cha4,Neu4)*
     -              DCONJG(SCha1(Cha4,2))*(MW*(dup9*mUSC(Sfe4,1,2)) + 
     -              3*CW*(mUSC(Sfe4,2,2)*SNeu1(Neu4,2)))*mUS(Sfe4,2,2))/
     -              SW)+ Cval(cc00,iint3(Cha4,Neu4,Sfe4))*(-(Alfa/(6.D0*
     -              CW*MW*Pi*SW2)*(EL*cnwQ(2,Cha4,Neu4)*
     -              DCONJG(SCha2(Cha4,1))*(MW*(dup9*mUSC(Sfe4,1,2)) + 3*
     -              CW*(mUSC(Sfe4,2,2)*SNeu1(Neu4,2)))*mUS(Sfe4,1,2))/
     -              (SW*Sqrt(2.D0))) + Alfa/(12.D0*CW*MW2*Pi*SW2)*(EL*
     -              cnwQ(2,Cha4,Neu4)*DCONJG(SCha1(Cha4,2))*(MW*(dup9*
     -              mUSC(Sfe4,1,2)) + 3*CW*(mUSC(Sfe4,2,2)*
     -              SNeu1(Neu4,2)))*mUS(Sfe4,2,2))/SW) 

               loopNL = loopNL + Cval(cc1,iint3(Cha4,Neu4,Sfe4))*
     -              (Alfa/(12.D0*CW*MW*Pi*SW2)*(dup11*dup8*EL*Qmass(2)*
     -              DCONJG(SCha2(Cha4,1))*mUS(Sfe4,1,2))/(SW*
     -              Sqrt(2.D0)) + Alfa/(12.D0*CW*MW*Pi*SW2)*(EL*
     -              (-Qmass2(1) + Qmass2(2))*cnwQ(2,Cha4,Neu4)*
     -              DCONJG(SCha2(Cha4,1))*(MW*(dup9*mUSC(Sfe4,1,2)) + 3*
     -              CW*(mUSC(Sfe4,2,2)*SNeu1(Neu4,2)))*mUS(Sfe4,1,2))/
     -              (SW*Sqrt(2.D0)) - Alfa/(24.D0*CW*MW2*Pi*SW2)*
     -              (dup11*dup8*EL*Qmass(2)*DCONJG(SCha1(Cha4,2))*
     -              mUS(Sfe4,2,2))/SW - Alfa/(24.D0*CW*MW2*Pi*SW2)*(EL*
     -              (-Qmass2(1) + Qmass2(2))*cnwQ(2,Cha4,Neu4)*
     -              DCONJG(SCha1(Cha4,2))*(MW*(dup9*mUSC(Sfe4,1,2)) + 3*
     -              CW*(mUSC(Sfe4,2,2)*SNeu1(Neu4,2)))*mUS(Sfe4,2,2))/
     -              SW - Alfa/(24.D0*CW*MW2*Pi*SW2)*(dup7*EL*Qmass(1)*
     -              (MW*(dup9*mUSC(Sfe4,1,2)) + 3*CW*(mUSC(Sfe4,2,2)*
     -              SNeu1(Neu4,2)))*mUS(Sfe4,1,2)*SCha1(Cha4,1))/SW)


               tmp1 = Alfa/(12.D0*CW*MW*Pi*SW2)*(dup11*dup8*EL*Qmass(1)*
     -              DCONJG(SCha2(Cha4,1))*mUS(Sfe4,1,2))/(SW*
     -              Sqrt(2.D0))+ Alfa/(12.D0*CW*MW*Pi*SW2)*(EL*Qmass(1)*
     -              Qmass(2)*cnwQ(2,Cha4,Neu4)*DCONJG(SCha2(Cha4,1))*
     -              (MW*(dup9*mUSC(Sfe4,1,2)) + 3*CW*(mUSC(Sfe4,2,2)*
     -              SNeu1(Neu4,2)))*mUS(Sfe4,1,2))/(SW*Sqrt(2.D0)) - 
     -              Alfa/(24.D0*CW*MW2*Pi*SW2)*(dup11*dup8*EL*Qmass(1)*
     -              DCONJG(SCha1(Cha4,2))*mUS(Sfe4,2,2))/SW - Alfa/
     -              (24.D0*CW*MW2*Pi*SW2)*(EL*Qmass(1)*Qmass(2)*
     -              cnwQ(2,Cha4,Neu4)*DCONJG(SCha1(Cha4,2))*(MW*(dup9*
     -              mUSC(Sfe4,1,2)) + 3*CW*(mUSC(Sfe4,2,2)*
     -              SNeu1(Neu4,2)))*mUS(Sfe4,2,2))/SW - Alfa/(24.D0*CW*
     -              MW2*Pi*SW2)*(dup11*EL*Qmass2(1)*cnwQ(1,Cha4,Neu4)*
     -              mUS(Sfe4,1,2)*SCha1(Cha4,1))/SW 
               
               tmp2 = Alfa/(12.D0*CW*MW*Pi*SW2)*(dup11*EL*Qmass(1)*
     -              cnwQ(1,Cha4,Neu4)*DCONJG(SCha2(Cha4,1))*MCha(Cha4)*
     -              mUS(Sfe4,1,2))/(SW*Sqrt(2.D0)) - Alfa/(24.D0*CW*MW2*
     -              Pi*SW2)*(dup11*EL*Qmass(1)*cnwQ(1,Cha4,Neu4)*
     -              DCONJG(SCha1(Cha4,2))*MCha(Cha4)*mUS(Sfe4,2,2))/SW + 
     -              Alfa/(24.D0*CW*MW2*Pi*SW2)*(dup11*EL*(-(Qmass2(1)*
     -              cnwQ(1,Cha4,Neu4)) + cnwQ(1,Cha4,Neu4)*
     -              mass2SF(Sfe4,2) - cnwQ(2,Cha4,Neu4)*MCha(Cha4)*
     -              MNeu(Neu4))*mUS(Sfe4,1,2)*SCha1(Cha4,1))/SW
               
               tmp2 = tmp2 + Alfa/(24.D0*CW*MW2*Pi*SW2)*(EL*Qmass(2)*
     -              cnwQ(2,Cha4,Neu4)*(MW*(dup9*mUSC(Sfe4,1,2)) + 3*CW*
     -              (mUSC(Sfe4,2,2)*SNeu1(Neu4,2)))*MCha(Cha4)*
     -              mUS(Sfe4,1,2)*SCha1(Cha4,1))/SW 
               
               loopNR = loopNR + tmp2*Cval(cc0,iint3(Cha4,Neu4,Sfe4)) + 
     -              tmp1*Cval(cc2,iint3(Cha4,Neu4,Sfe4)) + 
     -              Bval(bb0,iint1(Cha4,Neu4))*(Alfa/(24.D0*CW*MW2*Pi*
     -              SW2)*(dup11*EL*cnwQ(1,Cha4,Neu4)*mUS(Sfe4,1,2)*
     -              SCha1(Cha4,1))/SW)+Cval(cc00,iint3(Cha4,Neu4,Sfe4))*
     -              (-Alfa/(12.D0*CW*MW2*Pi*SW2)*(dup11*EL*
     -              cnwQ(1,Cha4,Neu4)*mUS(Sfe4,1,2)*SCha1(Cha4,1))/SW) 
               
               loopNR = loopNR + Cval(cc1,iint3(Cha4,Neu4,Sfe4))*(Alfa/
     -              (12.D0*CW*MW*Pi*SW2)*(dup11*dup8*EL*Qmass(1)*
     -              DCONJG(SCha2(Cha4,1))*mUS(Sfe4,1,2))/(SW*
     -              Sqrt(2.D0))- Alfa/(24.D0*CW*MW2*Pi*SW2)*(dup11*dup8*
     -              EL*Qmass(1)*DCONJG(SCha1(Cha4,2))*mUS(Sfe4,2,2))/
     -              SW + Alfa/(24.D0*CW*MW2*Pi*SW2)*(dup11*EL*
     -              (-Qmass2(1) + Qmass2(2))*cnwQ(1,Cha4,Neu4)*
     -              mUS(Sfe4,1,2)*SCha1(Cha4,1))/SW -Alfa/(24.D0*CW*MW2*
     -              Pi*SW2)*(dup7*EL*Qmass(2)*(MW*(dup9*mUSC(Sfe4,1,2))+ 
     -              3*CW*(mUSC(Sfe4,2,2)*SNeu1(Neu4,2)))*mUS(Sfe4,1,2)*
     -              SCha1(Cha4,1))/SW)
*               call loop2_Cha4_Neu4_Sfe4
            enddo
         enddo
      enddo

      do Sfe5 = 1,2
         do Sfe4 = 1,2
            do Neu4 = 1,4
               dup12 = 3*CW*(mUS(Sfe4,1,1)*SNeu1(Neu4,1)) - 
     -              2*MW*(SW*mUS(Sfe4,2,1)*SNeu3(1)*ZNeu(Neu4,1))
               dup13 = -(3*CW*(mUS(Sfe4,1,1)*SNeu1(Neu4,1))) + 
     -              2*MW*(SW*mUS(Sfe4,2,1)*SNeu3(1)*ZNeu(Neu4,1))
               dup14 = 3*CW*(SNeu2(2)*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)
               dup15 = 3*CW*(SNeu2(1)*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
               dup16 = MW*(dup15*mUS(Sfe4,1,1)) + 
     -              3*CW*(DCONJG(SNeu1(Neu4,1))*mUS(Sfe4,2,1))
               dup17 = 3*CW*(DCONJG(SNeu1(Neu4,2))*mUSC(Sfe5,1,2)) - 
     -              2*MW*(SW*mUSC(Sfe5,2,2)*SNeu3(2)*ZNeuC(Neu4,1))
               dup18 = -(3*CW*(DCONJG(SNeu1(Neu4,2))*mUSC(Sfe5,1,2))) + 
     -              2*MW*(SW*mUSC(Sfe5,2,2)*SNeu3(2)*ZNeuC(Neu4,1))


               loopPL = loopPL - Alfa/(36.D0*CW2*MW2*Pi*SW2)*(EL*
     -              Cval(cc00,iint4(Neu4,Sfe4,Sfe5))*(dup16*(MW*(dup14*
     -              mUSC(Sfe5,1,2)) + 3*CW*(mUSC(Sfe5,2,2)*
     -              SNeu1(Neu4,2))))*mUS(Sfe5,1,2)*mUSC(Sfe4,1,1))/(SW*
     -              Sqrt(2.D0)) 


               loopPR = loopPR - Alfa/(36.D0*CW2*MW2*Pi*SW2)*(EL*
     -              Cval(cc00,iint4(Neu4,Sfe4,Sfe5))*(dup12*dup17)*
     -              mUS(Sfe5,1,2)*mUSC(Sfe4,1,1))/(SW*Sqrt(2.D0)) 
*               call loop3_Neu4_Sfe4_Sfe5
            enddo
         enddo
      enddo


** Set individual parts to zero for debugging
*      loopNL = (0D0,0D0)
*      loopPL = (0D0,0D0)
*      loopNR = (0D0,0D0)
*      loopPR = (0D0,0D0)

         
      if (qqVcur .eq. 2) then
         loopNL = -loopNL
         loopNR = -loopNR
      end if

         
      qqWlp(1) = qqWlp(1) + loopNL + loopPL
      qqWlp(2) = qqWlp(2) + loopNR + loopPR  


      end


********************************************************************************
