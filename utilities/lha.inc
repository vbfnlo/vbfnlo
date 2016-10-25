c...output LHA event file name
      character*250 lhaFile, hepmcfile
      common/fileLHA/lhaFile,hepmcfile
c...unweigthing numbers
      integer nevall,nevacc,prenevunw,desiredEventCount
      common/unnumbers/ prenevunw,nevall,nevacc,desiredEventCount
c...unweigthing related
      real*8 maxweight_PS, sumwgt
      integer*8 nevraw
      integer nevovr, nevnewmax, prenevacc, estimatedEvents
      logical partialUnweight
      common/unwevents/ maxweight_PS, sumwgt, 
     &                  nevraw, 
     &                  estimatedEvents, nevovr, nevnewmax, prenevacc,
     &                  partialUnweight
c...multichannel-unweighting related
      integer maxContribs
      parameter(maxContribs=100) 
      real*8 contribXS(maxContribs), contribXSFraction(maxContribs), contribMaxWeight(maxContribs)
      integer contributionNumber, numContribs, contribEstimatedEvents(maxContribs)
      logical multiChannelUnweighting
      common/unwmultch/ contribXS, contribXSFraction, contribMaxWeight, 
     &                  contributionNumber, numContribs, contribEstimatedEvents,
     &                  multiChannelUnweighting
c...number of bosons
      integer n_b
      common/numbers/ n_b	   
c...input for Les Houches interface
      integer numParticles
      integer maxNumSubProc
      parameter(numParticles=8) 
      parameter(maxNumSubProc=3000) 
c...the most subprobable process for LHA interface 
      integer iprocess
      common /process/ iprocess	
c...helicities
      integer helicity
      common /helicities/helicity(6)
c...Les Houches switch
      logical lha,hepmc,doNLO,unweight,unweighting,taumass_switch
      common /LesHouches/ lha,hepmc,doNLO,unweight,unweighting,taumass_switch
c...cross section, error, etc. information
      double precision xsection, sdev, eventweight, maxweight,wtmax
      integer numcolor,numdecay	
      common /localHEPRUP/ xsection, sdev, eventweight, wtmax,
     &                     maxweight, numcolor, numdecay
c...local version of HEPEUP
      integer LMAXNUP
      parameter ( LMAXNUP = 20 )
      integer LNUP, LIDPRUP, LISTUP, LIDUP, LMOTHUP, LICOLUP, LSPINUP
      double precision LXWGTUP
      common /localHEPEUP/ LXWGTUP(maxNumSubProc),
     &               LNUP(maxNumSubProc), LIDPRUP(maxNumSubProc),
     &               LISTUP(LMAXNUP,maxNumSubProc),
     &               LIDUP(LMAXNUP,maxNumSubProc),
     &               LMOTHUP(2,LMAXNUP,maxNumSubProc),
     &               LSPINUP(LMAXNUP,maxNumSubProc),
     &               LICOLUP(2,LMAXNUP,maxNumSubProc)


C...User process event common block...      
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP 
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP),NUP,IDPRUP,IDUP(MAXNUP),ICOLUP(2,MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP)
      SAVE/HEPEUP/

C...User process initialization common block...
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/EBMUP(2),
     &XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &LPRUP(MAXPUP),IDBMUP(2),PDFGUP(2),PDFSUP(2),IDWTUP,NPRUP
      SAVE /HEPRUP/
