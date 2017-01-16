!...output LHA event file name
      character*250 lhaFile, hepmcfile
      common/fileLHA/lhaFile,hepmcfile
!...unweigthing numbers
      integer nevall,nevacc,prenevunw,desiredEventCount
      common/unnumbers/ prenevunw,nevall,nevacc,desiredEventCount
!...unweigthing related
      real*8 maxweight_PS, sumwgt
      integer*8 nevraw
      integer nevovr, nevnewmax, prenevacc, estimatedEvents
      logical partialUnweight

      common/unwevents/ maxweight_PS, sumwgt
      common/unwevents/ nevraw
      common/unwevents/ estimatedEvents, nevovr, nevnewmax, prenevacc
      common/unwevents/ partialUnweight
!...multichannel-unweighting related
      integer maxContribs
      parameter(maxContribs=100)
      real*8 contribXS(maxContribs), contribXSFraction(maxContribs), contribMaxWeight(maxContribs)
      integer contributionNumber, numContribs, contribEstimatedEvents(maxContribs)
      logical multiChannelUnweighting

      common/unwmultch/ contribXS, contribXSFraction, contribMaxWeight
      common/unwmultch/ contributionNumber, numContribs, contribEstimatedEvents
      common/unwmultch/ multiChannelUnweighting
!...number of bosons
      integer n_b
      common/numbers/ n_b
!...input for Les Houches interface
      integer numParticles
      integer maxNumSubProc
      parameter(numParticles=8)
      parameter(maxNumSubProc=4000)
!...the most subprobable process for LHA interface
      integer iprocess
      common /process/ iprocess
!...helicities
      integer helicity
      common /helicities/helicity(6)
!...Les Houches switch
      logical lha,hepmc,doNLO,unweight,unweighting,taumass_switch
      common /LesHouches/ lha,hepmc,doNLO,unweight,unweighting,taumass_switch
!...cross section, error, etc. information
      double precision xsection, sdev, eventweight, maxweight,wtmax
      integer numcolor,numdecay
      common /localHEPRUP/ xsection, sdev, eventweight, wtmax
      common /localHEPRUP/ maxweight, numcolor, numdecay
!...local version of HEPEUP
      integer LMAXNUP
      parameter ( LMAXNUP = 20 )
      integer LNUP, LIDPRUP, LISTUP, LIDUP, LMOTHUP, LICOLUP, LSPINUP
      double precision LXWGTUP
      common /localHEPEUP/ LXWGTUP(maxNumSubProc)
      common /localHEPEUP/ LNUP(maxNumSubProc), LIDPRUP(maxNumSubProc)
      common /localHEPEUP/ LISTUP(LMAXNUP,maxNumSubProc)
      common /localHEPEUP/ LIDUP(LMAXNUP,maxNumSubProc)
      common /localHEPEUP/ LMOTHUP(2,LMAXNUP,maxNumSubProc)
      common /localHEPEUP/ LSPINUP(LMAXNUP,maxNumSubProc)
      common /localHEPEUP/ LICOLUP(2,LMAXNUP,maxNumSubProc)


!...User process event common block...
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/ XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP(5,MAXNUP)
      COMMON/HEPEUP/ VTIMUP(MAXNUP),SPINUP(MAXNUP),NUP,IDPRUP,IDUP(MAXNUP),ICOLUP(2,MAXNUP)
      COMMON/HEPEUP/ ISTUP(MAXNUP),MOTHUP(2,MAXNUP)
      SAVE/HEPEUP/

!...User process initialization common block...
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/ EBMUP(2)
      COMMON/HEPRUP/ XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP)
      COMMON/HEPRUP/ LPRUP(MAXPUP),IDBMUP(2),PDFGUP(2),PDFSUP(2),IDWTUP,NPRUP
      SAVE /HEPRUP/
