! variables set in the input file
      LOGICAL writeROOT, writeTOP, writeGNU, writeDATA !, writePAW
      LOGICAL replaceFiles, smearbins, calcerror1d, calcerror2d, calcerrorGnuplot, calcerrorRoot
      character*50 rootfile, topfile, gnufile, datafile ! , pawfile
      double precision smearamount
      common/histuser/ writeROOT, writeTOP, writeGNU, writeDATA, replaceFiles
      common/histuser/ rootfile, topfile, gnufile, datafile !, writePAW, pawfile
      common/histsmear/ smearamount, smearbins, calcerror1d, calcerror2d, calcerrorGnuplot, calcerrorRoot

! properties of the histograms
      INTEGER HISTMAX, BINMAX, TOTAL

      PARAMETER (HISTMAX = 600, BINMAX = 128, TOTAL = HISTMAX * BINMAX)
      integer  nbin(HISTMAX)
      double precision hlo(HISTMAX), hup(HISTMAX), rwid(HISTMAX), dbin(HISTMAX)
      character*64 htitles(HISTMAX)

      integer NLOoffset
! list of histograms to be stored in file
      LOGICAL hstored(1:HISTMAX)
      common /histcb/ hlo, hup, rwid, dbin, nbin, NLOoffset, hstored, htitles

! properties of the two dimensional histograms
      INTEGER HIST2dMAX, BINMAXX, BINMAXY, NLOoffset2d

      PARAMETER (HIST2dMAX = 150, BINMAXX = 128, BINMAXY = 128)
      integer nbinx(HIST2dMAX), nbiny(HIST2dMAX)
      double precision hloX(HIST2dMAX), hupX(HIST2dMAX), hloY(HIST2dMAX), hupY(HIST2dMAX)
      double precision rwidX(HIST2dMAX), rwidY(HIST2dMAX), dbinX(HIST2dMAX), dbinY(HIST2dMAX)
      character*64 h2dtitles(HIST2dMAX)
      character*64 xlabels(HIST2dMAX), ylabels(HIST2dMAX)

      LOGICAL h2dstored(1:HIST2dMAX)
      common /hist2dcb/ hloX, hloY, hupX, hupY, rwidX, rwidY, dbinX, dbinY, nbinX, nbinY
      common /hist2dcb/ NLOoffset2d, h2dstored, h2dtitles, xlabels, ylabels

