****************************************************************************
* This file contains the subroutines for the histogram handling.
* There is an internal representation of the histograms, which is filled
* by calls of "FillHist". In the end of the program, this data is transferred
* to other histogram representations (e.g. hbook, top-drawer).
* Each histogram has an ID, which is used for the filling.
c***************************************************************************
c   LIST OF ALL FUNCTIONS AND SUBROUTINES IN THIS FILE:
c
c     SUBROUTINE InitHistograms
c     SUBROUTINE HistogramEvent(weight,dsig,dps, rand, p, x, v, jets, n_jets, 
c                                   NLO)
c     SUBROUTINE FillHist(ID, value, weight)
c     SUBROUTINE FillHist2d(ID, X, Y, weight, NLO)
c     SUBROUTINE CreateHist(ID, title, bins, min, max)
c     SUBROUTINE Create2dHist(ID, title, xlabel, ylabel, xbins, ybins, xmin,  
c                                  xmax, ymin, ymax)
c     SUBROUTINE WriteHistograms
c     SUBROUTINE WriteRootHists
c     SUBROUTINE WritePawHists
c     SUBROUTINE WriteTopHists
c     SUBROUTINE WriteGnuHists
c     SUBROUTINE WriteDataHists
c     SUBROUTINE SelectFile
c
c***************************************************************************
c
c      to add a histogram
c     --------------------
c     1: put a "CreateHist" in "InitHistograms". Use unique IDs for each 
c         histogram
c     2: put a proper "FillHist" in "HistogramEvents"
c
c***************************************************************************

c*************************************************************************
      SUBROUTINE InitHistograms
c*************************************************************************
c     Set up histograms
c*************************************************************************
          use globalvars, only: lglobalprint
          use readinput
          use hist_stor

      implicit none

#include "hist.inc"
#include "global.inc"
#include "lha.inc"

      integer id
      real*8 hwid
      character*5 ext

      ! No "reallocation" case needed, since InitHistograms is called only once.
      ! Furthermore, HISTMAX and HIST2dMAX are fixed parameters
      if (.not.allocated(hist))   allocate(hist(HISTMAX))
      if (.not.allocated(hist2d)) allocate(hist2d(HIST2dMAX))

      NLOoffset = HISTMAX/3
      NLOoffset2d = HIST2dMAX/3
c
c reset histogram contents to zero     
c 
      nbin(1:histmax) = 100
      hlo(1:histmax) = 0d0
      hup(1:histmax) = 1d0
      htitles(1:histmax) = ""
      hstored(1:histmax) = .false.

c
c do the same for the two-dimensional histograms
c
      nbinX(1:hist2dmax) = 100
      nbinY(1:hist2dmax) = 100
      hloX(1:hist2dmax) = 0d0
      hloY(1:hist2dmax) = 0d0
      hupX(1:hist2dmax) = 1d0
      hupY(1:hist2dmax) = 1d0
      h2dtitles(1:hist2dmax) = ""
      xlabels(1:hist2dmax) = ""
      ylabels(1:hist2dmax) = ""
      h2dstored(1:hist2dmax) = .false.

****************************
**    HISTOGRAM SETUP     **
****************************

c     loadfile has to be called before all CreateHist calls
      if (lglobalprint) then
         print *,""
         print *,"         information on histogram config"
         print *,"-----------------------------------------------"
         print *,""
      endif
      call loadfile("histograms.dat",.false.)
      call read_logical("SMEARING",smearbins,.true.)
      if (lglobalprint) write(*,*) "SMEARING           =     ", smearbins
      if (smearbins) then
        call read_Real("SMEAR_VALUE",smearamount,0.1d0)
        if (smearamount.gt.0.499 .or. smearamount.lt.1d-6) then
           print*, "SMEAR_VALUE is outside of allowed range"//
     &          " (10^-6 to 0.499)!"
           print*, "Default (0.1) is used."
           smearamount = 0.1
        endif
        write(*,*) "SMEAR_VALUE        =     ", smearamount
      endif
      call read_logical("CALC_ERROR_GNUPLOT",calcerrorGnuplot,.false.)
      if (lglobalprint) write(*,*) "CALC_ERROR_GNUPLOT =     ", calcerrorGnuplot
      call read_logical("CALC_ERROR_ROOT",calcerrorRoot,.false.)
      if (lglobalprint) write(*,*) "CALC_ERROR_ROOT    =     ", calcerrorRoot
      call read_logical("CALC_ERROR_1D",calcerror1d,.true.)
      if (lglobalprint) write(*,*) "CALC_ERROR_1D      =     ", calcerror1d
      call read_logical("CALC_ERROR_2D",calcerror2d,.false.)
      if (lglobalprint) write(*,*) "CALC_ERROR_2D      =     ", calcerror2d
      if (lglobalprint) print *,""

*     call CreateHist(ID, title, #bins, min, max)
      call CreateHist(1, "dS/dpT_j (fb/GeV)",100,0d0,250d0)    !pT of tagging jets
      call CreateHist(2, "dS/dpTmax_j (fb/GeV)",100,0d0,500d0) !pT of tagging jet with higher pT
      call CreateHist(3, "dS/dpTmin_j (fb/GeV)",100,0d0,200d0) !pT of tagging jet with lower pT
      call CreateHist(4, "dS/dy_j (fb)",100,-5d0,5d0)        !y of tagging jets
      call CreateHist(5, "dS/dy_j1 (fb)",100,-5d0,5d0)       !y of tagging jet with higher pT
      call CreateHist(6, "dS/dy_j2 (fb)",100,-5d0,5d0)       !y of tagging jet with lower pT
      call CreateHist(7, "dS/dpTmax_l (fb/GeV)",100,0d0,500d0)     !pT_max of all leptons
      call CreateHist(8, "dS/dpTmin_l (fb/GeV)",100,0d0,200d0)     !pT_min of all leptons
      call CreateHist(9, "dS/d|eta|max_l (fb)",100,0d0,5d0)  !|eta|_max of all leptons
      call CreateHist(10, "dS/d|eta|min_l (fb)",100,0d0,5d0)  !|eta|_min of all leptons
      call CreateHist(11, "dS/dPhi_jj (fb)", 100, -180d0, 180d0)  !azimuthal angle of tagging jets

c     example for 2d-hist
c     (ID, title, x-axis label, y-axis label, x #bins, y #bins, xmin, xmax, ymin, ymax)
      call Create2dHist(1, "d2S/dy_jj dm_jj (fb/GeV)", 
     -     "eta_jj", "m_jj", 25, 25, 0d0, 6d0, 0d0, 800d0)


* If you do NOT wish a histogram to be saved to file, switch it off here.
!      hstored( 1) = .false.  
!      hstored( 2) = .false.
!      hstored( 3) = .false.
!      hstored( 4) = .false.
!      hstored( 5) = .false.
!      hstored( 6) = .false.
!      hstored( 7) = .false.
!      hstored( 8) = .false.
!      hstored( 9) = .false.
!      hstored(10) = .false.
!      hstored(11) = .false.
!      h2dstored( 1) = .false.

c     closefile has to be called after all CreateHist calls
      call closefile

      if (lglobalprint) then
      print *,""
      print *,"-----------------------------------------------"
      print *,""
      print *,""
      print *,""
      endif

***********************************************************************
c     
c select also the NLO histograms to be stored in a file
c
      if ((hasNLO.and.doNLO) .or. ewcor_switch) then
         do id = 1, NLOoffset
            hstored(id+NLOoffset)=hstored(id)
            hstored(id+2*NLOoffset)=hstored(id)
         enddo

         do id = 1, NLOoffset2d
            h2dstored(id+NLOoffset2d)=h2dstored(id)
            h2dstored(id+2*NLOoffset2d)=h2dstored(id)
         enddo
       endif
c   
c calculate inverse bin widths for each histogram
c
      do id = 1, histmax
         hwid = hup(id) - hlo(id)
         
         if ( nbin(id).gt.0 .and. hwid.ne.0d0 ) then
            rwid(id) = real( nbin(id) ) / hwid
            dbin(id) = 1/rwid(id)
         else
            dbin(id) = 0d0
            rwid(id) = 0d0
         end if
      end do
c
c the same for the 2d histograms
c
      do id = 1, hist2dmax
         hwid = hupX(id) - hloX(id)
         if ( nbinX(id).gt.0 .and. hwid.ne.0d0 ) then
            rwidX(id) = real( nbinX(id) ) / hwid
            dbinX(id) = 1/rwidX(id)
         else
            dbinX(id) = 0d0
            rwidX(id) = 0d0
         endif
         hwid = hupY(id) - hloY(id)
         if ( nbinY(id).gt.0 .and. hwid.ne.0d0 ) then
            rwidY(id) = real( nbinY(id) ) / hwid
            dbinY(id) = 1/rwidY(id)
         else
            dbinY(id) = 0d0
            rwidY(id) = 0d0
         endif         
      end do
c
c read input file from "vbfnlo.dat"
c
      call loadfile("vbfnlo.dat",.false.)
      call read_string("XSECFILE",XSECfile,"xsection")
      call read_logical("ROOT",writeROOT,.false.)
*      call read_logical("HBOOK",writePAW,.false.)
      call read_logical("TOP",writeTOP,.false.)
      call read_logical("GNU",writeGNU,.true.)
      call read_logical("DATA",writeDATA,.true.)
      call read_logical("REPLACE",replaceFiles,.true.)
      call read_string("ROOTFILE",rootfile,"histograms")
*      call read_string("PAWFILE",PAWfile,"histograms")
      call read_string("TOPFILE",TOPfile,"histograms")
      call read_string("GNUFILE",GNUfile,"histograms")
      call read_string("DATAFILE",DATAfile,"histograms")
      call closefile
      if (.not.replacefiles) then
         ext = ".out"
         call selectfile(xsecfile,ext)
         ext = ".root"
         call selectfile(rootfile,ext)
*         ext = ".paw"
*         call selectfile(pawfile,ext)
         ext = ".top"
         call selectfile(topfile,ext)
         ext = ".gp"
         call selectfile(gnufile,ext)
         ext = ".dir"
         call selectfile(datafile,ext)
      else
        xsecfile= trim(xsecfile)//".out"
        rootfile= trim(rootfile)//".root"
*        pawfile= pawfile(1:index(pawfile," ")-1)//".paw"
        topfile= trim(topfile)//".top"        
        gnufile= trim(gnufile)//".gp"        
        datafile= trim(datafile)//".dir"        
      endif
#ifndef WITH_ROOT
      if (writeROOT) then
        write(*,*)'Sorry, ROOT output not included.'
        write(*,*)'Please reconfigure with --with-root or'//
     &       ' switch off ROOT histogram output!'
        writeROOT = .false.
        stop
      endif
#endif

c
c initialize variables for real emission error estimation
c
      ! this is not needed and messes up allocation as number of kinematics is not known
!      call InitRealHist(1)

#ifdef WITH_ROOT
      if (writeRoot) then 
        call openrootsession(rootfile)
        call InitRootHists()
      endif
#endif

      RETURN
      END


c*************************************************************************
      SUBROUTINE HistogramEvent(weight,dsig,dps, rand, p, x, v,
     &                           jets, n_jets, leptons, nl,
     &                           invisible,ninv, photons, nphotons, NLO)
c*************************************************************************
c     Fill histograms
c RECALL: jets/leptons/invisible labelling (j,i)                           
c             i = jet/lepton label number (1->max_jets/max_v)              
c             j = 0-3 = 4-momentum                                         
c                     = 4   = mass                                         
c                     = 5   = pT                                           
c                     = 6   = rapidity y                                   
c                     = 7   = azimuthal angle phi                          
c                     = 8   = (leptons only) particle ID                   
c*************************************************************************
         use hist_stor

      implicit none

#include "global.inc"
#include "hist.inc"
#include "cuts.inc"

      real*8 weight, dsig, dps
      integer NLO 
      real*8 p(0:3,max_p), v(0:3,max_v), x(nx)
      real*8 rand(max_PS_dim)
      real*8 jets(0:7,max_jets), leptons(0:8,max_v), photons(0:7,max_v)
      real*8 invisible(0:8,max_v)
      integer n_jets, nl, nphotons,i,ninv
      real*8 pTmaxlep, pTminlep, ymaxlep, yminlep

      real*8 dw
      real*8 deltaY, deltaPhi, mjj2, mjj

      double precision RSEPS

      external RSEPS
      external mjj2

* variables for real emission error estimation

#include "hist_real.inc"

      if (NLO.lt.-1) return

c bookkeeping for real emission histograms
      if (nlo .ge. 1) then
        do i = NLOoffset+1,2*NLOoffset
           hsave(i)%HIST_CALLS=0
        enddo
        do i = NLOoffset2d+1,2*NLOoffset2d
           hsave2d(i)%hist_calls=0
        enddo
      endif

      if (dsig.eq.0d0) return

c W-W-jj is obtained from W+W+jj:  Momenta should be flipped
      if (procID .eq. 260) call momswitch(p,v,jets,leptons,invisible)

c getting the IDs of the tagging jets
      if (vbfprocess .or. semilepAlwaysVBFCuts) call DefineTaggingJets(jets, n_jets)

      dw = dsig * weight


c fill histograms

**********************
* pT of tagging jets *
**********************
      if (tag1.gt.0 .and. tag2.gt.0) then
         call FillHist( 1, jets(5,tag1), dw*0.5d0, NLO)
         call FillHist( 1, jets(5,tag2), dw*0.5d0, NLO)
         call FillHist( 2, jets(5,tag1), dw, NLO)
         call FillHist( 3, jets(5,tag2), dw, NLO)
      else
         if (n_jets.ge.1) then
            call FillHist( 1, jets(5,1), dw*0.5d0, NLO)
            call FillHist( 2, jets(5,1), dw, NLO)
         endif
         if (n_jets.ge.2) then   
            call FillHist( 1, jets(5,2), dw*0.5d0, NLO)
            call FillHist( 3, jets(5,2), dw, NLO)
         end if
      endif

***********************
* eta of tagging jets *
***********************
      if (tag1.gt.0 .and. tag2.gt.0) then
         call FillHist( 4, jets(6,tag1), dw*0.5d0, NLO)
         call FillHist( 4, jets(6,tag2), dw*0.5d0, NLO)
         call FillHist( 5, jets(6,tag1), dw, NLO)
         call FillHist( 6, jets(6,tag2), dw, NLO)
      else
         if (n_jets.ge.1) then
            call FillHist( 4, jets(6,1), dw*0.5d0, NLO)
            call FillHist( 5, jets(6,1), dw, NLO)
         endif
         if (n_jets.ge.2) then   
            call FillHist( 4, jets(6,2), dw*0.5d0, NLO)
            call FillHist( 6, jets(6,2), dw, NLO)
         end if
      endif

**********************
* pT of leptons      *
**********************
      pTmaxlep = 0
      do i = 1, nl
         pTmaxlep = max(pTmaxlep,leptons(5,i))
      enddo
      pTminlep = 1.0d6
      do i = 1, nl
         pTminlep = min(pTminlep,leptons(5,i))
      enddo
      if (nl.gt.0) then  
         call FillHist( 7, pTmaxlep, dw, NLO)
         call FillHist( 8, pTminlep, dw, NLO)
      endif

*************************
* |eta| of leptons      *
*************************
      ymaxlep = 0
      do i = 1, nl
         ymaxlep = max(ymaxlep,abs(leptons(6,i)))
      enddo
      yminlep = 1.0d6
      do i = 1, nl
         yminlep = min(yminlep,abs(leptons(6,i)))
      enddo
      if (nl.gt.0) then  
         call FillHist( 9, ymaxlep, dw, NLO)
         call FillHist(10, yminlep, dw, NLO)
      endif

******************
*  Delta Phi_jj  *
******************
      if (n_jets.ge.2) then
         if (tag1.gt.0 .and. tag2.gt.0) then
            if (jets(6,tag1).gt.jets(6,tag2)) then
               deltaPhi = jets(7,tag1) - jets(7,tag2)
            else
               deltaPhi = jets(7,tag2) - jets(7,tag1)
            endif
         else
            if (jets(6,1).gt.jets(6,2)) then
               deltaPhi = jets(7,1) - jets(7,2)
            else
               deltaPhi = jets(7,2) - jets(7,1)
            endif
         endif
         if (deltaPhi.gt.Pi) then
            deltaPhi = deltaPhi - 2d0*pi
         else if (deltaPhi.lt.-Pi) then
            deltaPhi = deltaPhi + 2d0*pi
         endif
         deltaPhi = deltaPhi * 180d0 / Pi
         call FillHist(11, deltaPhi, dw, NLO)
      endif

****************************
* example for 2d-histogram *
****************************
      if (n_jets.ge.2) then
         deltay=abs(jets(6,1)-jets(6,2))
         mjj=sqrt(mjj2(jets(0,1),jets(0,2)))
         call FillHist2d( 1,deltay,mjj, dw, NLO)
      endif

#ifdef WITH_ROOT
      if (writeRoot) then 
        call FillRootHists(dw,jets,n_jets,leptons,nl,
     &         invisible,ninv,photons,nphotons,NLO)
      endif
#endif

      RETURN
      END


c*************************************************************************
      SUBROUTINE FillHist(ID, value, weight, NLO)
c*************************************************************************
c     Fill value in histogram. 
c     Smear over 2 bins if value is within fraction drb of bin boundary
c*************************************************************************
      use hist_stor

      implicit none

#include "hist.inc"
#include "global.inc"

      Integer ID, NLO
      real*8 value, weight
      Integer hist_in, BIN, binl, binh
      real*8 X,Y 
      real*8 xbin, yl, yh
      real*8 drb, tdrbi, rat

c variables for real emission error estimation

#include "hist_real.inc"

      integer i
      integer mem, mem_L, mem_H

      if (value .ne. value) return

      if (smearbins) then
        drb   = smearamount
        tdrbi = 0.5d0/drb
      else
        drb   = 0.0d0
        tdrbi = 1.0d0
      endif

      hist_in = ID
      if (NLO.ne.0) hist_in=hist_in+NLOoffset
      x = value
      y = weight*rwid(hist_in)

      if (.not.hstored(HIST_in)) return
      
      xbin = (X-HLO(hist_in)) * RWID(hist_in)
      if (xbin.gt.drb) then
         bin = min( INT(xbin), nbin(hist_in) ) 
         if ( (xbin-bin).lt.drb ) then ! x in low end of bin
            binl = bin+1
            binh = bin+2
            rat = (xbin-bin)*tdrbi + 0.5
            yh = y*rat
            yl = y-yh
            bin = 0
         elseif( (xbin-bin).lt.(1.0-drb) ) then ! x in main part of bin
            bin = bin+2
C            elseif( (xbin-bin).lt.1.0 ) then   ! x in high part of bin
         elseif( (xbin-bin).lt.1.0.and.xbin.lt.nbin(hist_in) ) then ! x in high part of bin
            binl = bin+2
            binh = bin+3
            rat = (xbin-bin-1.0)*tdrbi + 0.5
            yh = y*rat
            yl = y-yh
            bin = 0
         else                   ! x more than drp above upper bin
            bin = nbin(hist_in)+2
         endif
      elseif(xbin.lt.-drb) then
         bin = 1
      else                      ! xbin in [-drb,drb]
         binl = 1
         binh = 2
         rat = xbin * tdrbi + 0.5 ! in [0,1]
         yh = y*rat
         yl = y-yh
         bin = 0
      endif
      if (hist_in.eq.-1) then
         print*," "
         print*," input to hf1dp: x,y = ",value,weight
         print*," bin = ",bin
         if (bin.eq.0) then
            print*," bin_or. = ",min( INT(xbin), nbin(hist_in) ),
     1           " binl,h = ",binl,binh
            print*," xbin = ",xbin," rat = ",rat
            print*," yl = ",yl," yh = ",yh
         endif
      endif

      if (nlo .ge. 1 .and. (calcerror1d.or.calcerrorGnuplot.or.calcerrorRoot)) then
        hsave(hist_in)%HIST_CALLS = hsave(hist_in)%HIST_CALLS + 1
        hsave(hist_in)%HIST_CALLS_MAX = max(hsave(hist_in)%HIST_CALLS_MAX,hsave(hist_in)%HIST_CALLS)
        if (hsave(hist_in)%HIST_CALLS.gt.max_calls) then
            print*, "**************************************************************"
           print*, "The same histogram was called too often per phase space point!"
           print*, "Please increase max_calls in utilities/hist_real.inc,"
           print*, 'followed by a "make all install"!'
           print*, "**************************************************************"
           STOP
        endif
      endif

      if (bin.eq.0) then
        if (calcerror1d.or.calcerrorGnuplot.or.calcerrorRoot) then
         if (nlo .le. 0) then  ! LO + virt
          call inc(hist(hist_in)%hsum(BINL), YL)
          call inc(hist(hist_in)%hsum(BINH), YH)
          call inc(hist(hist_in)%hSQUARE(BINL), YL*YL)
          call inc(hist(hist_in)%hCOUNT(BINL), 1d0)
          call inc(hist(hist_in)%hSQUARE(BINH), YH*YH)
          call inc(hist(hist_in)%hCOUNT(BINH), 1d0)
         endif
         if (nlo .ge. 1) then ! reals
           mem_L=0
           mem_H=0
           do i=1,2*(nlo-1)
             if (hsave(HIST_in)%SAVE_BIN(i,hsave(HIST_in)%HIST_CALLS) .eq. BINL) mem_L = i
             if (hsave(HIST_in)%SAVE_BIN(i,hsave(HIST_in)%HIST_CALLS) .eq. BINH) mem_H = i
           enddo
           if (mem_L.eq.0) then
             mem_L = 2*nlo - 1
             hsave(HIST_in)%SAVE_BIN(mem_L,hsave(HIST_in)%HIST_CALLS) = BINL
           endif
           if (mem_H.eq.0) then
             mem_H = 2*nlo
             hsave(HIST_in)%SAVE_BIN(mem_H,hsave(HIST_in)%HIST_CALLS) = BINH
           endif
           call inc(hsave(HIST_in)%SAVE_Y(mem_L,hsave(HIST_in)%HIST_CALLS), YL)
           call inc(hsave(HIST_in)%SAVE_Y(mem_H,hsave(HIST_in)%HIST_CALLS), YH)
         endif
        else ! no error calculation needed
            call inc(hist(hist_in)%hsum(BINL), YL)
            call inc(hist(hist_in)%hsum(BINH), YH)
        endif
      else
        if (calcerror1d.or.calcerrorGnuplot.or.calcerrorRoot) then
         if (nlo .le. 0) then 
           call inc(hist(hist_in)%hSUM(BIN), Y)
           call inc(hist(hist_in)%hSQUARE(BIN), Y*Y)
           call inc(hist(hist_in)%hCOUNT(BIN), 1d0)
         endif
         if (nlo .ge. 1) then
           mem=0
           do i=1,2*(nlo-1)
             if (hsave(HIST_in)%SAVE_BIN(i,hsave(HIST_in)%HIST_CALLS) .eq. BIN) mem = i
           enddo
           if (mem.eq.0) then
             mem = 2*nlo
             hsave(HIST_in)%SAVE_BIN(mem,hsave(HIST_in)%HIST_CALLS) = BIN
           endif
           call inc(hsave(HIST_in)%SAVE_Y(mem,hsave(HIST_in)%HIST_CALLS), Y)
         endif
        else
           call inc(hist(hist_in)%hSUM(BIN), Y)
        endif
      endif
      
      RETURN
      END


c*************************************************************************
      SUBROUTINE FillHist2d(ID, X, Y, weight, NLO)
c*************************************************************************
c     Fill value in 2d-histogram. 
c     Smear over bins if value is within fraction drb of bin boundary
c*************************************************************************
         use hist_stor

      implicit none

#include "hist.inc"
#include "global.inc"

      Integer ID, NLO
      real*8 X, Y, weight
      Integer hist_in, BINX, binlx, binhx, BINY, binly, binhy
      real*8 xb, yb, wll, whl, wlh, whh, w, wh, wl, w_h, w_l
      real*8 drb, tdrbi, ratx, raty

c variables for real emission error estimation

#include "hist_real.inc"

      integer i,j
      integer mem_x, mem_y
      integer mem_x_L, mem_x_H, mem_y_L, mem_y_H

      if (weight .ne. weight) return

      if (smearbins) then
        drb   = smearamount
        tdrbi = 0.5d0/drb
      else
        drb   = 0.0d0
        tdrbi = 1.0d0
      endif

      hist_in = ID
      if (NLO.ne.0) hist_in=hist_in+NLOoffset2d
      w=weight*rwidX(hist_in)*rwidY(hist_in)
      xb = (X-HLOX(hist_in)) * RWIDX(hist_in)
      yb = (Y-HLOY(hist_in)) * RWIDY(hist_in)

      if (.not.h2dstored(HIST_in)) return

c     X-direction
      if (xb.gt.drb) then
         binx = min( INT(xb), nbinX(hist_in) ) 
         if ( (xb-binx).lt.drb ) then ! x in low end of bin    
            binlx = binx+1
            binhx = binx+2
            ratx = (xb-binx)*tdrbi + 0.5
            wh = w*ratx
            wl = w-wh
            binx = 0
         elseif( (xb-binx).lt.(1.0-drb) ) then ! x in main part of bin
            binx = binx+2
         elseif( (xb-binx).lt.1.0.and.xb.lt.nbinX(hist_in) ) then ! x in high part of bin
            binlx = binx+2
            binhx = binx+3
            ratx = (xb-binx-1.0)*tdrbi + 0.5
            wh = w*ratx
            wl = w-wh
            binx = 0
         else                   ! x more than drp above upper bin
            binx = nbinX(hist_in)+2
         endif
      elseif(xb.lt.-drb) then
         binx = 1
      else
         binlx = 1
         binhx = 2
         ratx = xb * tdrbi + 0.5 ! in [0,1]
         wh = w*ratx
         wl = w-wh
         binx = 0
      endif
c     Y-direction
      if (yb.gt.drb) then
         biny = min( INT(yb), nbinY(hist_in) ) 
         if ( (yb-biny).lt.drb ) then ! y in low end of bin    
            binly = biny+1
            binhy = biny+2
            raty = (yb-biny)*tdrbi + 0.5
            whh = wh*raty
            whl = wh-whh
            wlh = wl*raty
            wll = wl-wlh
            w_h = w*raty
            w_l = w-w_h
            biny = 0
         elseif( (yb-biny).lt.(1.0-drb) ) then ! y in main part of bin
            biny = biny+2
         elseif( (yb-biny).lt.1.0.and.yb.lt.nbinY(hist_in) ) then ! y in high part of bin
            binly = biny+2
            binhy = biny+3
            raty = (yb-biny-1.0)*tdrbi + 0.5
            whh = wh*raty
            whl = wh-whh
            wlh = wl*raty
            wll = wl-wlh
            w_h = w*raty
            w_l = w-w_h
            biny = 0
         else                   ! y more than drp above upper bin
            biny = nbinY(hist_in)+2
         endif
      elseif(yb.lt.-drb) then
         biny = 1
      else
         binly = 1
         binhy = 2
         raty = yb * tdrbi + 0.5 ! in [0,1]
         whh = wh*raty
         whl = wh-whh
         wlh = wl*raty
         wll = wl-wlh
         w_h = w*raty
         w_l = w-w_h
         biny = 0
      endif

      if (nlo .ge. 1 .and. (calcerror2d.or.calcerrorRoot)) then
        hsave2d(hist_in)%hist_calls = hsave2d(hist_in)%hist_calls + 1
        hsave2d(hist_in)%hist_calls_MAX = max(hsave2d(hist_in)%hist_calls_MAX,hsave2d(hist_in)%hist_calls)
        if (hsave2d(hist_in)%hist_calls.gt.max_calls) then
           print*, "**************************************************************"
           print*, "The same histogram was called too often per phase space point!"
           print*, "Please increase max_calls in utilities/hist_real.inc,"
           print*, 'followed by a "make all install"!'
           print*, "**************************************************************"
           STOP
        endif
      endif

      if (binx.eq.0) then
         if (biny.eq.0) then
           if (calcerror2d.or.calcerrorRoot) then
            if (nlo .le. 0) then 
             call inc(hist2d(HIST_IN)%hSUM(BINLX,BINLY), WLL)
             call inc(hist2d(HIST_in)%hSQUARE(BINLX,BINLY), WLL*WLL)
             call inc(hist2d(HIST_in)%hCOUNT(BINLX,BINLY), 1d0)

             call inc(hist2d(HIST_IN)%hSUM(BINHX,BINLY), WHL)
             call inc(hist2d(HIST_in)%hSQUARE(BINHX,BINLY), WHL*WHL)
             call inc(hist2d(HIST_in)%hCOUNT(BINHX,BINLY), 1d0)

             call inc(hist2d(HIST_IN)%hSUM(BINLX,BINHY), WLH)
             call inc(hist2d(HIST_in)%hSQUARE(BINLX,BINHY), WLH*WLH)
             call inc(hist2d(HIST_in)%hCOUNT(BINLX,BINHY), 1d0)

             call inc(hist2d(HIST_IN)%hSUM(BINHX,BINHY), WHH)
             call inc(hist2d(HIST_in)%hSQUARE(BINHX,BINHY), WHH*WHH)
             call inc(hist2d(HIST_in)%hCOUNT(BINHX,BINHY), 1d0)
            endif
            if (nlo .ge. 1) then
              mem_x_L = 0
              mem_x_H = 0
              mem_y_L = 0
              mem_y_H = 0
              do i=1,2*(nlo-1)
                if (hsave2d(HIST_in)%SAVE_BIN_X(i,hsave2d(HIST_in)%hist_calls) .eq. BINLX) mem_x_L = i
                if (hsave2d(HIST_in)%SAVE_BIN_X(i,hsave2d(HIST_in)%hist_calls) .eq. BINHX) mem_x_H = i
              enddo
              do j=1,2*(nlo-1)
                if (hsave2d(HIST_in)%SAVE_BIN_Y(j,hsave2d(HIST_in)%hist_calls) .eq. BINLY) mem_y_L = j
                if (hsave2d(HIST_in)%SAVE_BIN_Y(j,hsave2d(HIST_in)%hist_calls) .eq. BINHY) mem_y_H = j
              enddo
              if (mem_x_L.eq.0) then
                mem_x_L=2*nlo -1
                hsave2d(HIST_in)%SAVE_BIN_X(mem_x_L,hsave2d(HIST_in)%hist_calls) = BINLX
              endif
              if (mem_x_H.eq.0) then
                mem_x_H=2*nlo
                hsave2d(HIST_in)%SAVE_BIN_X(mem_x_H,hsave2d(HIST_in)%hist_calls) = BINHX
              endif
              if (mem_y_L.eq.0) then
                mem_y_L=2*nlo - 1
                hsave2d(HIST_in)%SAVE_BIN_Y(mem_y_L,hsave2d(HIST_in)%hist_calls) = BINLY
              endif
              if (mem_y_H.eq.0) then
                mem_y_H=2*nlo
                hsave2d(HIST_in)%SAVE_BIN_Y(mem_y_H,hsave2d(HIST_in)%hist_calls) = BINHY
              endif
              call inc(hsave2d(HIST_in)%save_y(mem_x_L,mem_y_L,hsave2d(HIST_in)%hist_calls), WLL) 
              call inc(hsave2d(HIST_in)%save_y(mem_x_L,mem_y_H,hsave2d(HIST_in)%hist_calls), WLH)
              call inc(hsave2d(HIST_in)%save_y(mem_x_H,mem_y_L,hsave2d(HIST_in)%hist_calls), WHL)
              call inc(hsave2d(HIST_in)%save_y(mem_x_H,mem_y_H,hsave2d(HIST_in)%hist_calls), WHH)
            endif
           else
            call inc(hist2d(HIST_IN)%hSUM(BINLX,BINLY), WLL)
            call inc(hist2d(HIST_IN)%hSUM(BINHX,BINLY), WHL)
            call inc(hist2d(HIST_IN)%hSUM(BINLX,BINHY), WLH)
            call inc(hist2d(HIST_IN)%hSUM(BINHX,BINHY), WHH)
           endif
         else
           if (calcerror2d.or.calcerrorRoot) then
            if (nlo .le. 0) then 
             call inc(hist2d(HIST_IN)%hSUM(BINLX,BINY), WL)
             call inc(hist2d(HIST_in)%hSQUARE(BINLX,BINY), WL*WL)
             call inc(hist2d(HIST_in)%hCOUNT(BINLX,BINY), 1d0)

             call inc(hist2d(HIST_IN)%hSUM(BINHX,BINY), WH)
             call inc(hist2d(HIST_in)%hSQUARE(BINHX,BINY), WH*WH)
             call inc(hist2d(HIST_in)%hCOUNT(BINHX,BINY), 1d0)
            endif
            if (nlo .ge. 1) then
              mem_x_L = 0
              mem_x_H = 0
              mem_y = 0
              do i=1,2*(nlo-1)
                if (hsave2d(HIST_in)%SAVE_BIN_X(i,hsave2d(HIST_in)%hist_calls) .eq. BINLX) mem_x_L = i
                if (hsave2d(HIST_in)%SAVE_BIN_X(i,hsave2d(HIST_in)%hist_calls) .eq. BINHX) mem_x_H = i
              enddo
              do j=1,2*(nlo-1)
                if (hsave2d(HIST_in)%SAVE_BIN_Y(j,hsave2d(HIST_in)%hist_calls) .eq. BINY) mem_y = j
              enddo
              if (mem_x_L.eq.0) then
                mem_x_L=2*nlo -1
                hsave2d(HIST_in)%SAVE_BIN_X(mem_x_L,hsave2d(HIST_in)%hist_calls) = BINLX
              endif
              if (mem_x_H.eq.0) then
                mem_x_H=2*nlo
                hsave2d(HIST_in)%SAVE_BIN_X(mem_x_H,hsave2d(HIST_in)%hist_calls) = BINHX
              endif
              if (mem_y.eq.0) then
                mem_y=2*nlo
                hsave2d(HIST_in)%SAVE_BIN_Y(mem_y,hsave2d(HIST_in)%hist_calls) = BINY
              endif
              call inc(hsave2d(HIST_in)%save_y(mem_x_L,mem_y,hsave2d(HIST_in)%hist_calls), WL)
              call inc(hsave2d(HIST_in)%save_y(mem_x_H,mem_y,hsave2d(HIST_in)%hist_calls), WH)
            endif
        else
            call inc(hist2d(HIST_IN)%hSUM(BINLX,BINY), WL)
            call inc(hist2d(HIST_IN)%hSUM(BINHX,BINY), WH)
           endif
         endif
      else
         if (biny.eq.0) then
           if (calcerror2d.or.calcerrorRoot) then
            if (nlo .le. 0) then 
             call inc(hist2d(HIST_IN)%hSUM(BINX,BINLY), W_L)
             call inc(hist2d(HIST_in)%hSQUARE(BINX,BINLY), W_L*W_L)
             call inc(hist2d(HIST_in)%hCOUNT(BINX,BINLY), 1d0)

             call inc(hist2d(HIST_IN)%hSUM(BINX,BINHY), W_H)
             call inc(hist2d(HIST_in)%hSQUARE(BINX,BINHY), W_H*W_H)
             call inc(hist2d(HIST_in)%hCOUNT(BINX,BINHY), 1d0)
            endif
            if (nlo .ge. 1) then
              mem_x = 0
              mem_y_L = 0
              mem_y_H = 0
              do i=1,2*(nlo-1)
                if (hsave2d(HIST_in)%SAVE_BIN_X(i,hsave2d(HIST_in)%hist_calls) .eq. BINX) mem_x = i
              enddo
              do j=1,2*(nlo-1)
                if (hsave2d(HIST_in)%SAVE_BIN_Y(j,hsave2d(HIST_in)%hist_calls) .eq. BINLY) mem_y_L = j
                if (hsave2d(HIST_in)%SAVE_BIN_Y(j,hsave2d(HIST_in)%hist_calls) .eq. BINHY) mem_y_H = j
              enddo
              if (mem_x.eq.0) then
                mem_x=2*nlo
                hsave2d(HIST_in)%SAVE_BIN_X(mem_x,hsave2d(HIST_in)%hist_calls) = BINX
              endif
              if (mem_y_L.eq.0) then
                mem_y_L=2*nlo - 1
                hsave2d(HIST_in)%SAVE_BIN_Y(mem_y_L,hsave2d(HIST_in)%hist_calls) = BINLY
              endif
              if (mem_y_H.eq.0) then
                mem_y_H=2*nlo
                hsave2d(HIST_in)%SAVE_BIN_Y(mem_y_H,hsave2d(HIST_in)%hist_calls) = BINHY
              endif
              call inc(hsave2d(HIST_in)%save_y(mem_x,mem_y_L,hsave2d(HIST_in)%hist_calls), W_L)
              call inc(hsave2d(HIST_in)%save_y(mem_x,mem_y_H,hsave2d(HIST_in)%hist_calls), W_H)
            endif
        else
            call inc(hist2d(HIST_IN)%hSUM(BINX,BINLY), W_L)
            call inc(hist2d(HIST_IN)%hSUM(BINX,BINHY), W_H)
           endif
         else
           if (calcerror2d.or.calcerrorRoot) then
            if (nlo .le. 0) then 
              call inc(hist2d(HIST_IN)%hSUM(BINX,BINY), W)
              call inc(hist2d(HIST_in)%hSQUARE(BINX,BINY), W*W)
              call inc(hist2d(HIST_in)%hCOUNT(BINX,BINY), 1d0)
            endif
            if (nlo .ge. 1) then
              mem_x = 0
              mem_y = 0
              do i=1,2*(nlo-1)
                if (hsave2d(HIST_in)%SAVE_BIN_X(i,hsave2d(HIST_in)%hist_calls) .eq. BINX) mem_x = i
              enddo
              do j=1,2*(nlo-1)
                if (hsave2d(HIST_in)%SAVE_BIN_Y(j,hsave2d(HIST_in)%hist_calls) .eq. BINY) mem_y = j
              enddo
              if (mem_x.eq.0) then
                mem_x=2*nlo
                hsave2d(HIST_in)%SAVE_BIN_X(mem_x,hsave2d(HIST_in)%hist_calls) = BINX
              endif
              if (mem_y.eq.0) then
                mem_y=2*nlo
                hsave2d(HIST_in)%SAVE_BIN_Y(mem_y,hsave2d(HIST_in)%hist_calls) = BINY
              endif
              call inc(hsave2d(HIST_in)%save_y(mem_x,mem_y,hsave2d(HIST_in)%hist_calls), W)
            endif
        else
              call inc(hist2d(HIST_IN)%hSUM(BINX,BINY), W)
           endif
         endif
      endif

      RETURN
      END


      SUBROUTINE InitRealHist(kinematics)
         use hist_stor

      implicit none

#include "hist.inc"
#include "global.inc"
#include "hist_real.inc"

      integer :: id ! id of histogram
      integer :: kinematics ! max number of kinematics
      integer :: i ! loops over lo,nlo,kfac

      real_kinematics = kinematics        ! real_kinematics = n_kin or n_qed, depending on the part currently calculated


      if (.not.allocated(hsave))   allocate(hsave(HISTMAX))
      if (.not.allocated(hsave2d)) allocate(hsave2d(HIST2dMAX))

      
      ! Since InitHistograms is called only once, no reallocation is needed: the only variable parameter in the allocation
      ! is hstored/h2dstored, which should be set only in InitHistograms


      do id = 1, NLOoffset
         if (hstored(id)) then
            do i = 0,1
               if (.not.allocated(hsave(id+i*NLOoffset)%save_y))   allocate(hsave(id+i*NLOoffset)%save_y(2*max_kin,max_calls))
               if (.not.allocated(hsave(id+i*NLOoffset)%save_bin)) allocate(hsave(id+i*NLOoffset)%save_bin(2*max_kin,max_calls))
               hsave(id+i*NLOoffset)%save_y = 0
               hsave(id+i*NLOoffset)%save_bin = 0
               hsave(id+i*NLOoffset)%hist_calls = 0
               hsave(id+i*NLOoffset)%hist_calls_max = 0
            enddo
         endif
      enddo


      do id = 1, NLOoffset2d
         if (h2dstored(id)) then
            do i = 0,1
               if (.not.allocated(hsave2d(id+i*NLOoffset2d)%save_y)) 
     &                     allocate(hsave2d(id+i*NLOoffset2d)%save_y(2*max_kin,2*max_kin,max_calls))
               if (.not.allocated(hsave2d(id+i*NLOoffset2d)%save_bin_x)) 
     &                     allocate(hsave2d(id+i*NLOoffset2d)%save_bin_x(2*max_kin,max_calls))
               if (.not.allocated(hsave2d(id+i*NLOoffset2d)%save_bin_y)) 
     &                     allocate(hsave2d(id+i*NLOoffset2d)%save_bin_y(2*max_kin,max_calls))
               hsave2d(id+i*NLOoffset2d)%save_y = 0
               hsave2d(id+i*NLOoffset2d)%save_bin_x = 0
               hsave2d(id+i*NLOoffset2d)%save_bin_y = 0
               hsave2d(id+i*NLOoffset2d)%hist_calls = 0
               hsave2d(id+i*NLOoffset2d)%hist_calls_max = 0
            enddo
         endif
      enddo


      RETURN
      END


      SUBROUTINE SaveRealHist()
         use hist_stor

      implicit none

#include "hist.inc"
#include "global.inc"
#include "hist_real.inc"

      integer calls,HID,i,j

      if (calcerror1d.or.calcerrorGnuplot.or.calcerrorRoot) then

      do HID = NLOoffset+1,2*NLOoffset
        if (hstored(HID)) then
          do calls=1,hsave(HID)%HIST_CALLS_MAX
            do i=1,2*real_kinematics
              if (hsave(HID)%SAVE_BIN(i,calls).ne.0 .and. hsave(HID)%SAVE_Y(i,calls).ne.0.d0) then
                call inc(hist(HID)%hSUM(hsave(HID)%SAVE_BIN(i,calls)), hsave(HID)%SAVE_Y(i,calls))
                call inc(hist(HID)%hSQUARE(hsave(HID)%SAVE_BIN(i,calls)), hsave(HID)%SAVE_Y(i,calls)**2)
                call inc(hist(HID)%hCOUNT(hsave(HID)%SAVE_BIN(i,calls)), 1d0)
              endif
              hsave(HID)%SAVE_BIN(i,calls)=0
              hsave(HID)%SAVE_Y(i,calls)=0.
            enddo
          enddo
        endif
        hsave(HID)%HIST_CALLS=0
        hsave(HID)%HIST_CALLS_MAX=0
      enddo

      endif

      if (calcerror2d.or.calcerrorRoot) then

      do HID = NLOoffset2d+1,2*NLOoffset2d
        if (h2dstored(HID)) then
          do calls=1,hsave2d(HID)%hist_calls_MAX
            do i=1,2*real_kinematics
             do j=1,2*real_kinematics
              if (hsave2d(HID)%SAVE_BIN_X(i,calls).ne.0 .and. 
     &               hsave2d(HID)%SAVE_BIN_Y(j,calls).ne.0 .and. 
     &               hsave2d(HID)%save_y(i,j,calls).ne.0.d0) then
                  call inc(hist2d(HID)%hSUM(hsave2d(HID)%SAVE_BIN_X(i,calls),hsave2d(HID)%SAVE_BIN_Y(j,calls)),
     &             hsave2d(HID)%save_y(i,j,calls))
                  call inc(hist2d(HID)%hSQUARE(hsave2d(HID)%SAVE_BIN_X(i,calls),hsave2d(HID)%SAVE_BIN_Y(j,calls)),
     &             hsave2d(HID)%save_y(i,j,calls)**2)
                  call inc(hist2d(HID)%hCOUNT(hsave2d(HID)%SAVE_BIN_X(i,calls),hsave2d(HID)%SAVE_BIN_Y(j,calls)), 1d0)
              endif
              hsave2d(HID)%save_y(i,j,calls)=0.
             enddo
             hsave2d(HID)%SAVE_BIN_X(i,calls)=0
            enddo
            do j=1,2*real_kinematics
             hsave2d(HID)%SAVE_BIN_Y(j,calls)=0
            enddo
          enddo
        endif
        hsave2d(HID)%hist_calls=0
        hsave2d(HID)%hist_calls_MAX=0
      enddo

      endif

      RETURN
      END


c*************************************************************************
      SUBROUTINE CreateHist(ID, title, bins, min, max)
c*************************************************************************
c     Set up the histogram for the internal representation.
c*************************************************************************
       use globalvars, only: lglobalprint
       use readinput, only: read_HistList
      use hist_stor

      implicit none

#include "hist.inc"

      integer ID
      character*(*) title
      integer bins
      real*8 min, max
      character*14 idstring
      character*7 idstringnr
      real*8 histconf(2), histconfdefault(2)
      character*7 myadjustl
      external myadjustl

      ! loop index
      integer i

      histconfdefault(1) = min
      histconfdefault(2) = max

      write( idstringnr, '(i3)' )  ID
      idstringnr = myadjustl(trim(idstringnr))
      idstring = "HIST_ID"//idstringnr

      
c check if internal arrays are large enough
      if (ID .gt. NLOoffset) then
        print*,"Error: Histogram id ",ID," larger than maximum:",HISTMAX,"/3"
        print*,"       Increase HISTMAX in utilities/hist.inc"
        stop
      endif
      if (bins .gt. BINMAX) then
        print*,"Error: Number of bins ",bins," for histogram ",ID,
     -        "larger than maximum:",BINMAX
        print*,"       Increase BINMAX in utilities/hist.inc"
        stop
      endif

** Reading in user-defined histogram limits
      call read_HistList(idstring,histconf,histconfdefault)
      if (lglobalprint) write(*,'(2A,2F10.2,2A)') idstring, "=", histconf, "  # ", title
      
      htitles(ID) = trim(title)
      write(htitles(ID+NLOoffset),*) title," (NLO)"
      write(htitles(ID+2*NLOoffset),*) title," (K)"


      do i=0,1
         ! check for reallocation: the user might have called CreateHist twice with the same histogram ID
         if (allocated(hist(id+i*NLOoffset)%hsum))    deallocate(hist(id+i*NLOoffset)%hsum)
         allocate(hist(id+i*NLOoffset)%hsum(bins+2))
         if (allocated(hist(id+i*NLOoffset)%hsquare)) deallocate(hist(id+i*NLOoffset)%hsquare)
         allocate(hist(id+i*NLOoffset)%hsquare(bins+2))
         if (allocated(hist(id+i*NLOoffset)%hcount))  deallocate(hist(id+i*NLOoffset)%hcount)
         allocate(hist(id+i*NLOoffset)%hcount(bins+2))
         ! todo check for enough memory
         ! stat=allocatestatus
   !      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
         hist(id+i*NLOoffset)%hsum = 0
         hist(id+i*NLOoffset)%hsquare = 0
         hist(id+i*NLOoffset)%hcount = 0
      enddo

      do i=0,2
         nbin(ID+i*NLOoffset) = bins
         hlo(ID+i*NLOoffset) = histconf(1)!min
         hup(ID+i*NLOoffset) = histconf(2)!max

      enddo

      hstored(ID) = .true.
     
      END


c*************************************************************************
      SUBROUTINE Create2dHist(ID, title, xlabel, ylabel, xbins, ybins, 
     -     xmin, xmax, ymin, ymax)
c*************************************************************************
c     Set up a 2d-histogram for the internal representation.
c*************************************************************************
       use globalvars, only: lglobalprint
       use readinput, only: read_Hist2dList
         use hist_stor

      implicit none

#include "hist.inc"

      integer ID
      character*(*) title, xlabel, ylabel
      integer xbins, ybins
      real*8 xmin, xmax, ymin, ymax
      
      character*14 idstring
      character*5 idstringnr
      real*8 histconf(4), histconfdefault(4)
      character*5 myadjustl
      external myadjustl
      
      ! loop index
      integer i
      
      histconfdefault(1) = xmin
      histconfdefault(2) = xmax
      histconfdefault(3) = ymin
      histconfdefault(4) = ymax
      
      write( idstringnr, '(i3)' )  ID
      idstringnr = myadjustl(trim(idstringnr))
      idstring = "2DHIST_ID"//idstringnr
      
c check if internal arrays are large enough
      if (ID .gt. NLOoffset2d) then
        print*,"Error: 2d-Histogram id ",ID," larger than maximum:",HIST2dMAX,"/3"
        print*,"       Increase HIST2dMAX in utilities/hist.inc"
        stop
      endif
      if (xbins .gt. BINMAXX) then
        print*,"Error: Number of x-bins ",xbins," for 2d-histogram ",ID,"larger than maximum:",BINMAXX
        print*,"       Increase BINMAXX in utilities/hist.inc"
        stop
      endif
      if (ybins .gt. BINMAXY) then
        print*,"Error: Number of y-bins ",ybins," for 2d-histogram ",ID,"larger than maximum:",BINMAXY
        print*,"       Increase BINMAXY in utilities/hist.inc"
        stop
      endif

** Reading in user-defined axis limits
      call read_Hist2dList(idstring,histconf,histconfdefault)
      if (lglobalprint) write(*,'(2A,4F10.2,2A)') idstring, "=", histconf, "  # ", title
      
      h2dtitles(ID) = title
      write(h2dtitles(ID+NLOoffset2d),*) title," (NLO)"
      write(h2dtitles(ID+2*NLOoffset2d),*) title," (K)"

      do i=0,2 ! loop over LO, NLO, kfac
         xlabels(ID+i*NLOoffset2d) = xlabel
         ylabels(ID+i*NLOoffset2d) = ylabel

         nbinX(ID+i*NLOoffset2d) = xbins
         hloX(ID+i*NLOoffset2d) = histconf(1) !xmin
         hupX(ID+i*NLOoffset2d) = histconf(2) !xmax
         nbinY(ID+i*NLOoffset2d) = ybins
         hloY(ID+i*NLOoffset2d) = histconf(3) !ymin
         hupY(ID+i*NLOoffset2d) = histconf(4) !ymax      
      enddo

      do i=0,1 ! loop over LO, NLO

         ! check for reallocation: the user might have called CreateHist twice with the same histogram ID
         if (allocated(hist2d(id+i*NLOoffset2d)%hsum))    deallocate(hist2d(id+i*NLOoffset2d)%hsum)
         allocate(hist2d(id+i*NLOoffset2d)%hsum(xbins+2,ybins+2))
         if (allocated(hist2d(id+i*NLOoffset2d)%hsquare)) deallocate(hist2d(id+i*NLOoffset2d)%hsquare)
         allocate(hist2d(id+i*NLOoffset2d)%hsquare(xbins+2,ybins+2))
         if (allocated(hist2d(id+i*NLOoffset2d)%hcount))  deallocate(hist2d(id+i*NLOoffset2d)%hcount)
         allocate(hist2d(id+i*NLOoffset2d)%hcount(xbins+2,ybins+2))
         ! todo check for enough memory
         ! stat=allocatestatus
   !      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
         hist2d(id+i*NLOoffset2d)%hsum = 0
         hist2d(id+i*NLOoffset2d)%hsquare = 0
         hist2d(id+i*NLOoffset2d)%hcount =0
      enddo

      h2dstored(ID) = .true.
     
      RETURN
      END


c*************************************************************************
      SUBROUTINE WriteHistograms
c*************************************************************************
c     write out histograms in various formats
c*************************************************************************

      implicit none

#include "hist.inc"    

      print*,' '

#ifdef WITH_ROOT
      if (writeRoot) then 
        call CloseRootHists
        call WriteRootHists
      endif
#endif

* #ifdef WITH_CERN
*      if (writePaw) call WritePawHists
* #endif

      if (writeTop) call WriteTopHists

      if (writegnu) call WriteGnuHists

      if (writedata) then
        call WriteDataHists
        call WriteGnuCommandFile
      endif

      RETURN
      END


#ifdef WITH_ROOT
c*************************************************************************
      SUBROUTINE WriteRootHists
c*************************************************************************
c     Transfers the histogram data to the Root format
c*************************************************************************
         use hist_stor

      implicit none

#include "hist.inc"

      integer i,ID,TID,histdp,bin,binx,biny
      double precision x,y,z, error
      character*4 htype
      
      DO TID = 1,NLOoffset
         IF (hstored(TID)) THEN
            do i=0,1
              ID = TID+i*NLOoffset
              if (i.eq.0) then
                htype = "LO"
              else if (i.eq.1) then
                htype = "NLO"
              endif
              call createroothisto(TID,ID,htitles(ID),
     &               nbin(ID),hlo(ID),hup(ID),htype)
              histdp = ID 
              do bin = 1,nbin(ID)+2
                 ! if (count(histdp,bin).gt.0d0) then
                 if (hist(histdp)%hsum(bin).ne.0d0) then
                    x = hlo(ID) + dbin(ID)*(bin-1.5d0)
                    y = hist(histdp)%hsum(bin)
                    if (calcerrorRoot .and. hist(histdp)%hcount(bin) .gt. 0d0) then
                      error = sqrt(hist(histdp)%hsquare(bin) - hist(histdp)%hsum(bin)**2/hist(histdp)%hcount(bin))
                    else
                      error = 0d0
                    endif
                    call fillroothisto(ID,x,y,error)
                 endif
              enddo
            enddo
         ENDIF
      ENDDO

      DO TID = 1,NLOoffset2d
         IF (h2dstored(TID)) THEN
            do i=0,1
              ID = TID+i*NLOoffset2d
              if (i.eq.0) then
                htype = "LO"
              else if (i.eq.1) then
                htype = "NLO"
              endif
              call createroothisto2d(TID,ID,htitles(ID),
     &               nbinx(ID),hlox(ID),hupx(ID),
     &               nbiny(ID),hloy(ID),hupy(ID),htype)
              histdp = ID 
              do binx = 1,nbinx(ID)+2
                 do biny = 1,nbiny(ID)+2
                    ! if (count2d(histdp,binx,biny).gt.1d-3) then
                    if (hist2d(histdp)%hsum(binx,biny).ne.0d0) then
                       x = hlox(ID) + dbinx(ID)*(binx-1.5d0)
                       y = hloy(ID) + dbiny(ID)*(biny-1.5d0)
                       z = hist2d(histdp)%hsum(binx,biny)
                       if (calcerrorRoot .and. hist2d(histdp)%hcount(binx,biny) .gt. 0d0) then
                         error = sqrt(hist2d(histdp)%hsquare(binx,biny) 
     &                         - hist2d(histdp)%hsum(binx,biny)**2/hist2d(histdp)%hcount(binx,biny))
                       else
                         error = 0d0
                       endif
                       call fillroothisto2d(ID,x,y,z,error)
                    endif
                 enddo
              enddo
            enddo
         ENDIF
      ENDDO

      call closerootsession
      print *," created root file for histograms : ", rootfile
      
      RETURN
      END
#endif


c*************************************************************************
      SUBROUTINE WriteTopHists
c*************************************************************************
c     Transfers the histogram data to the top-drawer format
c*************************************************************************
         use hist_stor

      implicit none

#include "hist.inc"

      integer ID, ibin, i1

      double precision delx,xres(0:BINMAX+2,2), xmin,xmax,ymin,ymax,ymaxK,yminK

      OPEN(UNIT=20,FILE=topfile,STATUS="REPLACE")

      write (20,"(A)") " set size 11 by 8.5 "

      DO ID = 1,HISTMAX
         IF (hstored(ID)) THEN
            write (20,*)
            write (20,"(a)") " new frame "

            write (20,"(A)") " set intensity 4 "
            write (20,"(A)") " set font duplex " 
            write (20,"(A)") " set tick size .05 "
            write (20,"(A)") " set title size 2.0 "
            
            delx = (hup(id)-hlo(id))/(real(nbin(id))+1e-10)
            XRES(0,1) = hlo(id) - DELX/2. ! stores underflow
            XRES(1,1) = hlo(id) + DELX/2.
            DO  IBIN = 2,nbin(id)+1
               XRES(IBIN,1) = XRES(IBIN-1,1) + DELX
            END DO
            xmin = hlo(id)
            xmax = hup(id)
            ymin = 1d33
            ymax = -1d33
            yminK = 0
            ymaxK = 2
            if (id .le. 2*NLOoffset) then
                do ibin = 0,nbin(id)+1  !ibin=nbin+1 for overflow
                   xres(ibin,2) = hist(id)%hsum(ibin+1)               
                enddo
                ymin = minval(xres(0:nbin(id)+1, 2))
                ymax = maxval(xres(0:nbin(id)+1, 2))
            else !kfac
                do ibin = 0,nbin(id)+1  !ibin=nbin+1 for overflow
                   xres(ibin,2) = 
     &               hist(id-NLOoffset)%hsum(ibin+1)/hist(id-2*NLOoffset)%hsum(ibin+1)
                enddo
                ymin = minval(xres(0:nbin(id)+1, 2))
                ymax = maxval(xres(0:nbin(id)+1, 2))
            endif
            do i1 = 0,nbin(id)+1
               if (xres(i1,2).eq.0.) xres(i1,2) = 1e-20
            end do
            write (20,"(A,i3,a)") " title 2.0 7.5 ""Fig.",
     1           id,": "//htitles(id)
            write (20,"(A)") " set window  x 1.2 to 5.2 y 3 to 7.0 "
            write (20,202) " set limits x ",xmin," to ",xmax,
     1                             " y ",ymin," to ",ymax
 202        format(2(a,f8.2),2(a,f12.5)) 
            write(20,*) " ( "
            DO  IBIN = 1,NBIN(ID)
               WRITE (20,"(2g20.6)")  XRES(ibin,1), xres(ibin,2)
            END DO
            write (20,"(A)") " join solid "
                       
         ENDIF
      ENDDO

      CLOSE(20)

      print *," created top-drawer file for histograms : ", topfile

      RETURN
      END


c*************************************************************************
      SUBROUTINE WriteGnuHists
c*************************************************************************
c     Transfers the histogram data to the gnuplot format
c*************************************************************************
         use hist_stor

      implicit none

#include "hist.inc"

      integer ID, ibin, jbin

      double precision xres(1:BINMAX+2,3)
      real*8, dimension(:,:) :: xres2d(1:BINMAXX+2,1:BINMAXY+2)

      OPEN(UNIT=20,FILE=gnufile,STATUS="REPLACE")

      write (20,"(A)") 'set output "'//
     &gnufile(1:len_trim(gnufile)-3)//'.ps"'
      write (20,"(A)") "set terminal postscript color"
      write (20,"(A)") "set style data histep"
      write (20,"(A)") "# set style data lines"
      write (20,"(A)") "set key off"

c     the 1d-histograms
      DO ID = 1,HISTMAX
         IF (hstored(ID)) THEN
            write (20,*)
            write (20,*) "#",ID,": ", trim(htitles(id))
            write (20,"(A)") 'set title "'// trim(htitles(id))//'"'
            if (calcerrorGnuplot) then
              write (20,"(A)") 'plot "-", "-" with yerrorbars lt 1 pt 0'
            else
              write (20,"(A)") 'plot "-" '
            endif
            
            DO  IBIN = 1,nbin(id)+2
               XRES(IBIN,1) = hlo(id) + dbin(id)*(IBIN-1.5d0)
            END DO

            if (id .le. 2*NLOoffset) then
c LO and NLO
              do ibin = 1,nbin(id)+2  !ibin=nbin+2 for overflow
                 xres(ibin,2) = hist(id)%hsum(ibin)               
                 if (calcerrorGnuplot .and. hist(id)%hcount(ibin) .gt. 0d0) then
                   xres(ibin,3) = hist(id)%hsquare(ibin) - 
     &                             hist(id)%hsum(ibin)**2/hist(id)%hcount(ibin)
                 else
                   xres(ibin,3) = 0d0
                 endif
              enddo
            else
c K-factor
              do ibin = 1,nbin(id)+2  !ibin=nbin+2 for overflow
                 if (hist(id-2*NLOoffset)%hsum(ibin).ne.0) then
                   xres(ibin,2) = 
     &               hist(id-NLOoffset)%hsum(ibin)/hist(id-2*NLOoffset)%hsum(ibin)
                   if (calcerrorGnuplot .and. hist(id-2*NLOoffset)%hsum(ibin).ne.0) then
                     xres(ibin,3) = 
     &                hist(id-2*NLOoffset)%hsquare(ibin)
     &                /hist(id-2*NLOoffset)%hsum(ibin)**2 
     &                - 1d0/hist(id-2*NLOoffset)%hcount(ibin)
                   else
                     xres(ibin,3) = 0
                   endif
                   if (calcerrorGnuplot .and. hist(id-NLOoffset)%hsum(ibin).ne.0) then
                     xres(ibin,3) = xres(ibin,3) + 
     &                  hist(id-NLOoffset)%hsquare(ibin)
     &                  /hist(id-NLOoffset)%hsum(ibin)**2 
     &                  - 1d0/hist(id-NLOoffset)%hcount(ibin)
                   endif
                   xres(ibin,3) = xres(ibin,3) * xres(ibin,2)**2 
                 else
                   xres(ibin,2) = -1
                   xres(ibin,3) = 0
                 endif
              enddo
            endif

            DO  IBIN = 2,NBIN(ID)+1    ! without under-/overflow bins
               WRITE (20,"(2g20.6)")  XRES(ibin,1), xres(ibin,2)
            END DO
            write (20,"(A)") "e"
            if (calcerrorGnuplot) then
              DO  IBIN = 2,NBIN(ID)+1
                 WRITE (20,"(3g20.6)")  XRES(ibin,1), xres(ibin,2), sqrt(xres(ibin,3))
              END DO
              write (20,"(A)") "e"
            endif
          ENDIF
      ENDDO

c     the 2d-histograms
      DO ID = 1, HIST2dMAX
         IF (h2dstored(ID)) THEN
            write (20,*)
            write (20,*) "#",ID,": ", trim(h2dtitles(id))
            write (20,"(A)") "set view map"
            write (20,"(A)") "set size 3.9/5.0,3.4/3.0"
            write (20,"(A)") "unset surface"
            write (20,"(A)") "set style data pm3d"
            write (20,"(A)") "set style function pm3d"
            write (20,"(A)") "set pm3d implicit at b"            
            write (20,"(A)") "set pm3d map corners2color c1"
            write (20,"(A)") 'set title "'// trim(h2dtitles(id))//'"'
            write (20,"(A)") 'set xlabel "'// trim(xlabels(id))//'"'
            write (20,"(A)") 'set ylabel "'// trim(ylabels(id))//'"'
c            write (20,"(A)") "set palette rgbformulae 7, 5, 15 negative"
            write (20,"(A)") "set palette rgbformulae 30,31,32 negative"

            write (20,"(A)") 'splot "-" '
            
            xres2d = 0      ! overflow bins should be 0 for proper colour scaling
            if (id .le. 2*NLOoffset2d) then
              DO IBIN = 1, NBINX(ID)+1
                 DO JBIN = 1, NBINY(ID)+1
                   xres2d(ibin,jbin) = hist2d(id)%hsum(ibin,jbin)
                 ENDDO
              ENDDO
            else
              DO IBIN = 1, NBINX(ID)+1
                 DO JBIN = 1, NBINY(ID)+1
                   if (hist2d(id-2*NLOoffset2d)%hsum(ibin,jbin).ne.0) then
                     xres2d(ibin,jbin) = 
     &           hist2d(id-NLOoffset2d)%hsum(ibin,jbin)/hist2d(id-2*NLOoffset2d)%hsum(ibin,jbin)
                   else
                     xres2d(ibin,jbin) = -1
                   endif
                 ENDDO
              ENDDO
            endif
                  
            DO IBIN = 2, NBINX(ID)+2       ! corners2color c1 is usually lower left corner -> write this out here
               DO JBIN = 2, NBINY(ID)+2    ! gnuplot needs overflow bin as well
                 WRITE (20,"(3g20.6)")  hlox(id)+dbinx(ID)*(IBIN-2d0), 
     &                  hloy(id)+dbiny(ID)*(JBIN-2d0), xres2d(ibin,jbin)
               ENDDO
               write (20,"(A)") " "
            ENDDO
            write (20,"(A)") "e"
         ENDIF
      ENDDO

      CLOSE(20)
      print *," created GNUplot file for histograms : ", gnufile
      
      RETURN
      END


c*************************************************************************
      SUBROUTINE WriteGnuCommandFile
c*************************************************************************
c     Write a Gnuplot command file into the data folder
c*************************************************************************

      implicit none

#include "hist.inc"

      integer ID
      character*255 tmpstr
      character*50 charid
      character*50 myadjustl
      external myadjustl

      write(tmpstr,*) 'mkdir -p '//trim(datafile)
      CALL system(tmpstr)

      write(tmpstr,'(A)') trim(datafile)//"/hist_noerror.gp"
      OPEN(UNIT=20,FILE=tmpstr,STATUS="REPLACE")

      write (20,"(A)") 'set output "hist_noerror.ps"'
      write (20,"(A)") "set terminal postscript color"
      write (20,"(A)") "set style data histep"
      write (20,"(A)") "# set style data lines"
      write (20,"(A)") "set key off"

      if (calcerror1d) then
         write(tmpstr,'(A)') trim(datafile)//"/hist_witherror_1d.gp"
         OPEN(UNIT=21,FILE=tmpstr,STATUS="REPLACE")
         write (21,"(A)") 'set output "hist_witherror_1d.ps"'
         write (21,"(A)") "set terminal postscript color"
         write (21,"(A)") "set style data histep"
         write (21,"(A)") "# set style data lines"
         write (21,"(A)") "set key off"
      endif

c     the 1d-histograms
      DO ID = 1,HISTMAX
         IF (hstored(ID)) THEN
            write (20,*)
            write (20,*) "#",ID,": ", trim(htitles(id))
            write (20,"(A)") 'set title "'//trim(htitles(id))//'"'
            write (20,*) "set xrange [",hlo(id),":",hup(id),"]"
            if (calcerror1d) then
               write (21,*)
               write (21,*) "#",ID,": ", trim(htitles(id))
               write (21,"(A)") 'set title "'//trim(htitles(id))//'"'
               write (21,*) "set xrange [",hlo(id),":",hup(id),"]"
            endif
            if (id .le. NLOoffset) then
               ! LO
               write(charid,*) ID
               charid = myadjustl(charid)
               write (20,"(A)") 'plot "LO/hist.'//trim(charid)//'.dat" using (($1+$2)/2):3'
               if (calcerror1d)
     &            write (21,"(A)") 'plot "LO/hist.'//trim(charid)//'.dat" using (($1+$2)/2):3, '//
     &                                  '"LO/hist.'//trim(charid)//'.dat" using (($1+$2)/2):3:4 with yerrorbars lt 1 pt 0'
            elseif (id .le. 2*NLOoffset) then
               ! NLO
               write(charid,*) ID-NLOoffset
               charid = myadjustl(charid)
               write (20,"(A)") 'plot "NLO/hist.'//trim(charid)//'.dat" using (($1+$2)/2):3'
               if (calcerror1d)
     &            write (21,"(A)") 'plot "NLO/hist.'//trim(charid)//'.dat" using (($1+$2)/2):3, '//
     &                                  '"NLO/hist.'//trim(charid)//'.dat" using (($1+$2)/2):3:4 with yerrorbars lt 1 pt 0'
            else
               ! K-factor
               write(charid,*) ID-2*NLOoffset
               charid = myadjustl(charid)
               write (20,"(A)") 'plot "Kfac/hist.'//trim(charid)//'.dat" using (($1+$2)/2):3'
               if (calcerror1d)
     &            write (21,"(A)") 'plot "Kfac/hist.'//trim(charid)//'.dat" using (($1+$2)/2):3, '//
     &                                  '"Kfac/hist.'//trim(charid)//'.dat" using (($1+$2)/2):3:4 with yerrorbars lt 1 pt 0'
            endif
          ENDIF
      ENDDO

c     the 2d-histograms
      DO ID = 1, HIST2dMAX
         IF (h2dstored(ID)) THEN
            write (20,*)
            write (20,*) "#",ID,": ", trim(h2dtitles(id))
            write (20,"(A)") "set view map"
            write (20,"(A)") "set size 3.9/5.0,3.4/3.0"
            write (20,"(A)") "unset surface"
            write (20,"(A)") "set style data pm3d"
            write (20,"(A)") "set style function pm3d"
            write (20,"(A)") "set pm3d implicit at b"            
            write (20,"(A)") "set pm3d map corners2color c1"
            write (20,"(A)") 'set title "'//trim(h2dtitles(id))//'"'
            write (20,"(A)") 'set xlabel "'//trim(xlabels(id))//'"'
            write (20,"(A)") 'set ylabel "'//trim(ylabels(id))//'"'
            write (20,"(A)") "set palette rgbformulae 30,31,32 negative"
            write (20,*) "set xrange [",hlox(id),":",hupx(id),"]"
            write (20,*) "set yrange [",hloy(id),":",hupy(id),"]"
            if (id .le. NLOoffset2d) then
               ! LO
               write(charid,*) ID
               charid = myadjustl(charid)
               write (20,"(A)") 'splot "LO/hist2.'//trim(charid)//'.dat" using 1:3:5'
            elseif (id .le. 2*NLOoffset2d) then
               ! NLO
               write(charid,*) ID-NLOoffset2d
               charid = myadjustl(charid)
               write (20,"(A)") 'splot "NLO/hist2.'//trim(charid)//'.dat" using 1:3:5'
            else
               ! K-factor
               write(charid,*) ID-2*NLOoffset2d
               charid = myadjustl(charid)
               write (20,"(A)") 'splot "Kfac/hist2.'//trim(charid)//'.dat" using 1:3:5'
            endif
         ENDIF
      ENDDO

      CLOSE(20)
      if (calcerror1d) CLOSE(21)
      print *," created GNUplot command file for in the histogram data folder"
      
      RETURN
      END


c*************************************************************************
      SUBROUTINE WriteDataHists
c*************************************************************************
c     Transfers the histogram data to a general data format
c*************************************************************************
         use hist_stor

      implicit none

#include "hist.inc"

      integer ID, ibin, jbin


      real*8, dimension(:,:) :: xres(0:BINMAX+2,3) 

      real*8, dimension(:,:,:) :: xres2d(0:BINMAXX+2,0:BINMAXY+2,3) 
      character*255 tmpstr
      character*50  charid
      character*50 myadjustl
      external myadjustl

      write(tmpstr,*) 'mkdir -p '//trim(datafile)//"/LO"
      CALL system(tmpstr)
      write(tmpstr,*) 'mkdir -p '//trim(datafile)//"/NLO"
      CALL system(tmpstr)
      write(tmpstr,*) 'mkdir -p '//trim(datafile)//"/Kfac"
      CALL system(tmpstr)

c     the 1d-histograms
      DO ID = 1,HISTMAX
         IF (hstored(ID)) THEN

            if (id .le. NLOoffset) then
               write(charid,*) ID
               charid = myadjustl(charid)
               write(tmpstr,'(A)') trim(datafile)//
     &                             "/LO/hist."//
     &                             trim(charid)//".dat"
            else if (id .le. 2*NLOoffset) then
               write(charid,*) ID-NLOoffset
               charid = myadjustl(charid)
               write(tmpstr,'(A)') trim(datafile)//
     &                             "/NLO/hist."//
     &                             trim(charid)//".dat"
            else
               write(charid,*) ID-2*NLOoffset
               charid = myadjustl(charid)
               write(tmpstr,'(A)') trim(datafile)//
     &                             "/Kfac/hist."//
     &                             trim(charid)//".dat"
            endif
            OPEN(UNIT=20,FILE=tmpstr,STATUS="REPLACE")
            write (20,*) "#",ID,": ", trim(htitles(id))

            DO  IBIN = 1,nbin(id)+2
               XRES(IBIN,1) = hlo(id) + dbin(id)*(IBIN-1.5d0)
            END DO

            if (id .le. 2*NLOoffset) then
c LO and NLO
              do ibin = 1,nbin(id)+2  !ibin=nbin+2 for overflow
                 xres(ibin,2) = hist(id)%hsum(ibin)               
                 if (calcerrorGnuplot .and. hist(id)%hcount(ibin) .gt. 0d0) then
                   xres(ibin,3) = hist(id)%hsquare(ibin) - 
     &                             hist(id)%hsum(ibin)**2/hist(id)%hcount(ibin)
                 else
                   xres(ibin,3) = 0d0
                 endif
              enddo
            else
c K-factor
              do ibin = 1,nbin(id)+2  !ibin=nbin+2 for overflow
                 if (hist(id-2*NLOoffset)%hsum(ibin).ne.0) then
                   xres(ibin,2) = 
     &               hist(id-NLOoffset)%hsum(ibin)/hist(id-2*NLOoffset)%hsum(ibin)
                   if (calcerrorGnuplot .and. hist(id-2*NLOoffset)%hsum(ibin).ne.0) then
                     xres(ibin,3) = 
     &                hist(id-2*NLOoffset)%hsquare(ibin)
     &                /hist(id-2*NLOoffset)%hsum(ibin)**2 
     &                - 1d0/hist(id-2*NLOoffset)%hcount(ibin)
                   else
                     xres(ibin,3) = 0
                   endif
                   if (calcerrorGnuplot .and. hist(id-NLOoffset)%hsum(ibin).ne.0) then
                     xres(ibin,3) = xres(ibin,3) + 
     &                  hist(id-NLOoffset)%hsquare(ibin)
     &                  /hist(id-NLOoffset)%hsum(ibin)**2 
     &                  - 1d0/hist(id-NLOoffset)%hcount(ibin)
                   endif
                   xres(ibin,3) = xres(ibin,3) * xres(ibin,2)**2 
                 else
                   xres(ibin,2) = -1
                   xres(ibin,3) = 0
                 endif
              enddo
            endif
            do ibin = 1,nbin(id)+2
               if (xres(ibin,3).lt.0d0) xres(ibin,3) = 0d0
            end do

            if (calcerror1d) then
              DO  IBIN = 1,NBIN(ID)+2
                 WRITE (20,"(4g20.6)")  hlo(id)+dbin(id)*(IBIN-2), 
     &                                  hlo(id)+dbin(id)*(IBIN-1), 
     &                                  xres(ibin,2), sqrt(xres(ibin,3))
              END DO
            else
              DO  IBIN = 1,NBIN(ID)+2
                 WRITE (20,"(3g20.6)")  hlo(id)+dbin(id)*(IBIN-2), 
     &                                  hlo(id)+dbin(id)*(IBIN-1), 
     &                                  xres(ibin,2)
              END DO
            endif

            close(20)
            
          ENDIF
      ENDDO

c     the 2d-histograms
      DO ID = 1, HIST2dMAX
         IF (h2dstored(ID)) THEN
            if (id .le. NLOoffset2d) then
               write(charid,*) ID
               charid = myadjustl(charid)
               write(tmpstr,'(A)') trim(datafile)// "/LO/hist2."// trim(charid)//".dat"
            else if (id .le. 2*NLOoffset2d) then
               write(charid,*) ID-NLOoffset2d
               charid = myadjustl(charid)
               write(tmpstr,'(A)') trim(datafile)// "/NLO/hist2."// trim(charid)//".dat"
            else
               write(charid,*) ID-2*NLOoffset2d
               charid = myadjustl(charid)
               write(tmpstr,'(A)') trim(datafile)// "/Kfac/hist2."//trim(charid)//".dat"
            endif
            OPEN(UNIT=20,FILE=tmpstr,STATUS="REPLACE")
            write (20,*) "#",ID,": ", trim(h2dtitles(id))
            
            if (id .le. 2*NLOoffset2d) then
              DO IBIN = 1, NBINX(ID)+2
                 DO JBIN = 1, NBINY(ID)+2
                   xres2d(ibin,jbin,2) = hist2d(id)%hsum(ibin,jbin)
                   if (calcerror2d .and. hist2d(id)%hcount(ibin,jbin) .gt. 0d0) then
                     xres2d(ibin,jbin,3) = 
     &                 hist2d(id)%hsquare(ibin,jbin) - 
     &                 hist2d(id)%hsum(ibin,jbin)**2/hist2d(id)%hcount(ibin,jbin)
                   else
                     xres2d(ibin,jbin,3) = 0d0
                   endif
                 ENDDO
              ENDDO
            else
              DO IBIN = 1, NBINX(ID)+2
                 DO JBIN = 1, NBINY(ID)+2
                  if (hist2d(id-2*NLOoffset2d)%hsum(ibin,jbin).ne.0) then
                    xres2d(ibin,jbin,2) = 
     &        hist2d(id-NLOoffset2d)%hsum(ibin,jbin)/hist2d(id-2*NLOoffset2d)%hsum(ibin,jbin)
                    if (calcerror2d .and. hist2d(id-2*NLOoffset2d)%hsum(ibin,jbin).ne.0) then
                      xres2d(ibin,jbin,3) = 
     &                  hist2d(id-2*NLOoffset2d)%hsquare(ibin,jbin)
     &                  /hist2d(id-2*NLOoffset2d)%hsum(ibin,jbin)**2
     &                  - 1d0/hist2d(id-2*NLOoffset2d)%hcount(ibin,jbin)
                    else
                      xres2d(ibin,jbin,3) = 0
                    endif
                    if (calcerror2d .and. hist2d(id-NLOoffset2d)%hsum(ibin,jbin).ne.0) then
                      xres2d(ibin,jbin,3) = xres2d(ibin,jbin,3) + 
     &                  hist2d(id-NLOoffset2d)%hsquare(ibin,jbin)
     &                  /hist2d(id-NLOoffset2d)%hsum(ibin,jbin)**2
     &                  - 1d0/hist2d(id-NLOoffset2d)%hcount(ibin,jbin)
                    endif
                    xres2d(ibin,jbin,3) = xres2d(ibin,jbin,3) * 
     &                                    xres2d(ibin,jbin,2)**2 
                  else
                    xres2d(ibin,jbin,2) = -1
                    xres2d(ibin,jbin,3) = 0
                  endif
                 ENDDO
              ENDDO
            endif
            do ibin = 1,nbinx(id)+2
              do jbin = 1,nbiny(id)+2
               if (xres2d(ibin,jbin,3).lt.0d0) xres2d(ibin,jbin,3) = 0d0
              end do
            end do
                  
            if (calcerror2d) then
              DO IBIN = 1, NBINX(ID)+2
                 DO JBIN = 1, NBINY(ID)+2
                   WRITE (20,"(6g20.6)")  
     &               hlox(id)+dbinx(ID)*(IBIN-2d0), 
     &               hlox(id)+dbinx(ID)*(IBIN-1d0), 
     &               hloy(id)+dbiny(ID)*(JBIN-2d0), 
     &               hloy(id)+dbiny(ID)*(JBIN-1d0), 
     &               xres2d(ibin,jbin,2),
     &               sqrt(xres2d(ibin,jbin,3))
                 ENDDO
                 write (20,"(A)") " "
              ENDDO
            else
              DO IBIN = 1, NBINX(ID)+2
                 DO JBIN = 1, NBINY(ID)+2
                   WRITE (20,"(5g20.6)")  
     &               hlox(id)+dbinx(ID)*(IBIN-2d0), 
     &               hlox(id)+dbinx(ID)*(IBIN-1d0), 
     &               hloy(id)+dbiny(ID)*(JBIN-2d0), 
     &               hloy(id)+dbiny(ID)*(JBIN-1d0), 
     &               xres2d(ibin,jbin,2)
                 ENDDO
                 write (20,"(A)") " "
              ENDDO
            endif
            CLOSE(20)
         ENDIF
      ENDDO

      print *," created data files for histograms : ", datafile
      
      RETURN
      END


ccccccccccccccccccccc begin subroutine SelectFile ccccccccccccccc
c
c  This subroutine selects a file name which does not exist yet
c  Terrance Figy <terrance@pheno.physics.wisc.edu
c  August 17, 2001  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine SelectFile(name,extens)

      implicit none

      character*50 name,lname,lname1
      character*1 digit(0:9)
      character*5 extens
      integer ii,i,j
      logical FileAlreadyExists
      data DIGIT /"0","1","2","3","4","5","6","7","8","9"/
 
      lname1 = name
      do 20 ii = 1,99
     
         if((ii.ge.1).and.(ii.le.9)) then
            
            lname = 
     &           trim(lname1)//"."//DIGIT(ii)//extens
            inquire(file=lname,exist=FileAlreadyExists)
            if(FileAlreadyExists) then
               goto 20
            else
               name = lname(1:index(lname," ")-1)
               return
            endif
         elseif(ii.gt.9) then
            i = (ii/10)
            j = mod(ii,10)
            lname = trim(lname1)//"."//
     &              DIGIT(i)//DIGIT(j)//extens
            inquire(file=lname,exist=FileAlreadyExists)
            if(FileAlreadyExists) then
               goto 20
            else
               name = trim(lname)
               return
            endif
         endif

         name = trim(name)
         
 20      continue
         write(*,*) "-------- fatal crash in SelectFile ----------"
         stop
         end


ccccccccccccccccccccc function myadjustl ccccccccccccccc
c
c implementation of FORTRAN90 intrinsic adjustl
c returns character string left-aligned
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function myadjustl(instr)
      implicit none
      character(*) instr, myadjustl

      myadjustl = instr
      do while (myadjustl(1:1) == ' ') 
        myadjustl = myadjustl(2:len_trim(myadjustl))
      enddo

      return
      end
