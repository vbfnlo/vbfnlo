********************************************************************************
********************************************************************************
*** mstw2008.F                                                               ***
*** 7th June 2011                                                            ***
*** sophy@particle.uni-karlsruhe.de                                          ***
***                                                                          ***
*** Interface for MSTW 2008 PDFs.  Almost entirely from the MSTW2008 package ***
*** from HEPFORGE:  http://projects.hepforge.org/mstwpdf/code/code.html      ***
***                                                                          ***
*** INFORMATION FROM MSTW:                                                   ***
***   reference: A. D. Martin, W. J. Stirling, R. S. Thorne and G. Watt,     ***
***               "Parton distributions for the LHC",                        ***
***               Eur. Phys. J. C63 (2009) 189-285                           ***
***               [arXiv:0901.0002 [hep-ph]].                                ***
***   Fortran interpolation code for MSTW PDFs, building on existing         ***
***   MRST Fortran code and Jeppe Andersen's C++ code.                       ***
***     user interface:                                                      ***
***        call GetAllPDFsAlt(prefix,ih,x,q,xpdf,xphoton)                    ***
***   Comments to Graeme Watt <Graeme.Watt(at)cern.ch>.                      ***
***   Updated 25/06/2010: Enlarge allowed range for m_c and m_b.             ***
***   Updated 25/01/2011: Fix "NaN" bug for q <= m_c when m_c^2 < 1.25 GeV^2 ***
***                                                                          ***
********************************************************************************
********************************************************************************

      subroutine init_mstw(nlo)

      implicit none

      integer nlo

      integer mstw_nlo
      common /mstw/ mstw_nlo


      mstw_nlo = nlo

** Initialisation statements
      if (nlo .eq. 0) then
         print*,"Initialize PDF set MSTW2008 LO"
      else
         print*,"Initialize PDF set MSTW2008 NLO"
      end if

c$$$** Initialising alphas
c$$$      if (nlo .eq. 0) then
c$$$         call mstwINITALPHAS(nlo,1D0,1D0,0.68183D0,1.40D0,4.75D0,1.D10)
c$$$      else
c$$$         call mstwINITALPHAS(nlo,1D0,1D0,0.49128D0,1.40D0,4.75D0,1.D10)
c$$$      end if

      
      end


********************************************************************************
********************************************************************************

      subroutine mstw2008(x, q, pdf)

      implicit none

      double precision x, q, pdf(-6:6)
      double precision xphoton

      integer mstw_nlo
      common /mstw/ mstw_nlo

* declarations of path variables
      character*250 GetPdfsetsPath
      character prefix*150
      character*250 path
      character*40 datafile
      external GetPdfsetsPath


** Setting prefix for grid files
      path=GetPdfsetsPath()
      if (mstw_nlo .eq. 0) then
         datafile = 'mstw2008lo'
         prefix = trim(path)//"/"//datafile
      else if (mstw_nlo .eq. 1) then
         datafile = 'mstw2008nlo'
         prefix = trim(path)//"/"//datafile
      end if


*   An LHAPDF-like interface where the PDFs
*   are returned in an array pdf(-6:6).  Indices given by PDG
*   convention (apart from the gluon is 0, not 21):
*                -6,  -5,  -4,  -3,  -2,  -1,0,1,2,3,4,5,6
*            = tbar,bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b,t.
      call mstwGetAllPDFsAlt(prefix,0,x,q,pdf,xphoton)


      end


********************************************************************************
********************************************************************************

C--   Alternative LHAPDF-like interface: return PDFs in an array.

      subroutine mstwGetAllPDFsAlt(prefix,ih,x,q,xpdf,xphoton)

      implicit none

      integer ih,f
      double precision x,q,xpdf(-6:6),xphoton,xvalence,mstwGetOnePDF
      character*(*) prefix

      do f = 1, 6
         xpdf(f) = mstwGetOnePDF(prefix,ih,x,q,f) ! quarks
         xvalence = mstwGetOnePDF(prefix,ih,x,q,f+6) ! valence quarks
         xpdf(-f) = xpdf(f) - xvalence ! antiquarks
      end do
      xpdf(0) = mstwGetOnePDF(prefix,ih,x,q,0) ! gluon
      xphoton = mstwGetOnePDF(prefix,ih,x,q,13) ! photon
      
      return

      end


********************************************************************************
********************************************************************************

C--   Get only one parton flavour 'f', using PDG notation
C--   (apart from gluon has f=0, not 21):
C--    f =   -6,  -5,  -4,  -3,  -2,  -1,0,1,2,3,4,5,6
C--      = tbar,bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b,t.
C--   Can also get valence quarks directly:
C--    f =  7, 8, 9,10,11,12.
C--      = dv,uv,sv,cv,bv,tv.
C--   Photon: f = 13.

      double precision function mstwGetOnePDF(prefix,ih,x,q,f)

      implicit none

      logical warn,fatal
      parameter(warn=.false.,fatal=.false.)
C--   Set warn=.true. to turn on warnings when extrapolating.
C--   Set fatal=.false. to return zero instead of terminating when
C--    invalid input values of x and q are used.

      integer ih,f,nhess,nx,nq,np,nqc0,nqb0,n,m,ip,io,
     &     alphaSorder,alphaSnfmax,nExtraFlavours,mstw_lentrim
      double precision x,q,xmin,xmax,qsqmin,qsqmax,mc2,mb2,eps,
     &     dummy,qsq,xlog,qsqlog,res,res1,anom,mstwExtrapolatePDF,
     &     mstwInterpolatePDF,distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ
      parameter(nx=64,nq=48,np=12)
      parameter(xmin=1d-6,xmax=1d0,qsqmin=1d0,qsqmax=1d9,eps=1d-6)
      parameter(nhess=2*20)
      character set*2,prefix*(*),filename*200,oldprefix(0:nhess)*200
      character dummyChar,dummyWord*50
      double precision ff(np,nx,nq)
      double precision qqorig(nq),qq(nq),xx(nx),cc(np,0:nhess,nx,nq,4,4)
      double precision xxl(nx),qql(nq)

C--   Store distance along each eigenvector, tolerance,
C--   heavy quark masses and alphaS parameters in COMMON block.
      common/mstwCommonS/distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ,alphaSorder,alphaSnfmax
      save

      data xx/1d-6,2d-6,4d-6,6d-6,8d-6,
     &     1d-5,2d-5,4d-5,6d-5,8d-5,
     &     1d-4,2d-4,4d-4,6d-4,8d-4,
     &     1d-3,2d-3,4d-3,6d-3,8d-3,
     &     1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     &     .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     &     .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     &     .5d0,.525d0,.55d0,.575d0,.6d0,.625d0,.65d0,.675d0,
     &     .7d0,.725d0,.75d0,.775d0,.8d0,.825d0,.85d0,.875d0,
     &     .9d0,.925d0,.95d0,.975d0,1d0/
      data qqorig/1.d0,
     &     1.25d0,1.5d0,0.d0,0.d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,
     &     1d1,1.2d1,0.d0,0.d0,2.6d1,4d1,6.4d1,1d2,
     &     1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     &     1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     &     1.8d6,3.2d6,5.6d6,1d7,1.8d7,3.2d7,5.6d7,1d8,
     &     1.8d8,3.2d8,5.6d8,1d9/

      if (f.lt.-6.or.f.gt.13) then
         print *,"Error: invalid parton flavour = ",f
         stop
      end if

      if (ih.lt.0.or.ih.gt.nhess) then
         print *,"Error: invalid eigenvector number = ",ih
         stop
      end if

C--   Check if the requested parton set is already in memory.
      if (oldprefix(ih).ne.prefix) then

C--   Start of initialisation for eigenvector set "i" ...
C--   Do this only the first time the set "i" is called,
C--   OR if the prefix has changed from the last time.

C--   Check that the character arrays "oldprefix" and "filename"
C--   are large enough.
         if (mstw_lentrim(prefix).gt.len(oldprefix(ih))) then
            print *,"Error in mstwGetOnePDF: increase size of oldprefix"
            stop
         else if (mstw_lentrim(prefix)+7.gt.len(filename)) then
            print *,"Error in mstwGetOnePDF: increase size of filename"
            stop
         end if

         write(set,'(I2.2)') ih  ! convert integer to string
C--   Remove trailing blanks from prefix before assigning filename.
         filename = prefix(1:mstw_lentrim(prefix))//'.'//set//'.dat'
C--   Line below can be commented out if you don't want this message.
         print *,"Reading PDF grid from ",filename(1:mstw_lentrim(filename))
         open(unit=33,file=filename,iostat=io,status='old')
         if (io.ne.0) then
            print *,"Error in mstwGetOnePDF: can't open ",
     &           filename(1:mstw_lentrim(filename))
            stop
         end if

C--   Read header containing heavy quark masses and alphaS values.
         read(33,*) 
         read(33,*)
         read(33,*) dummyChar,dummyWord,dummyWord,dummyChar,
     &        distance,tolerance
         read(33,*) dummyChar,dummyWord,dummyChar,mCharm
         read(33,*) dummyChar,dummyWord,dummyChar,mBottom
         read(33,*) dummyChar,dummyWord,dummyChar,alphaSQ0
         read(33,*) dummyChar,dummyWord,dummyChar,alphaSMZ
         read(33,*) dummyChar,dummyWord,dummyWord,dummyChar,
     &        alphaSorder,alphaSnfmax
         read(33,*) dummyChar,dummyWord,dummyChar,nExtraFlavours
         read(33,*)
         read(33,*)
         mc2=mCharm**2
         mb2=mBottom**2

C--   Check that the heavy quark masses are sensible.
C--   Redistribute grid points if not in usual range.
         do m=1,nq
            qq(m) = qqorig(m)
         end do
         if (mc2.le.qq(1).or.mc2+eps.ge.qq(8)) then
            print *,"Error in mstwGetOnePDF: invalid mCharm = ",mCharm
            stop
         else if (mc2.lt.qq(2)) then
            nqc0=2
            qq(4)=qq(2)
            qq(5)=qq(3)
         else if (mc2.lt.qq(3)) then
            nqc0=3
            qq(5)=qq(3)
         else if (mc2.lt.qq(6)) then
            nqc0=4
         else if (mc2.lt.qq(7)) then
            nqc0=5
            qq(4)=qq(6)
         else
            nqc0=6
            qq(4)=qq(6)
            qq(5)=qq(7)
         end if
         if (mb2.le.qq(12).or.mb2+eps.ge.qq(17)) then
            print *,"Error in mstwGetOnePDF: invalid mBottom = ",mBottom
            stop
         else if (mb2.lt.qq(13)) then
            nqb0=13
            qq(15)=qq(13)
         else if (mb2.lt.qq(16)) then
            nqb0=14
         else
            nqb0=15
            qq(14)=qq(16)
         end if
         qq(nqc0)=mc2
         qq(nqc0+1)=mc2+eps
         qq(nqb0)=mb2
         qq(nqb0+1)=mb2+eps

C--   The nExtraFlavours variable is provided to aid compatibility
C--   with future grids where, for example, a photon distribution
C--   might be provided (cf. the MRST2004QED PDFs).
         if (nExtraFlavours.lt.0.or.nExtraFlavours.gt.1) then
            print *,"Error in mstwGetOnePDF: invalid nExtraFlavours = ",
     &           nExtraFlavours
            stop
         end if

C--   Now read in the grids from the grid file.
         do n=1,nx-1
            do m=1,nq
               if (nExtraFlavours.gt.0) then
                  if (alphaSorder.eq.2) then ! NNLO
                     read(33,'(12(1pe12.4))',iostat=io)
     &                    (ff(ip,n,m),ip=1,12)
                  else          ! LO or NLO
                     ff(10,n,m) = 0.d0 ! = chm-cbar
                     ff(11,n,m) = 0.d0 ! = bot-bbar
                     read(33,'(10(1pe12.4))',iostat=io)
     &                    (ff(ip,n,m),ip=1,9),ff(12,n,m)
                  end if
               else             ! nExtraFlavours = 0
                  if (alphaSorder.eq.2) then ! NNLO
                     ff(12,n,m) = 0.d0 ! = photon
                     read(33,'(11(1pe12.4))',iostat=io)
     &                 (ff(ip,n,m),ip=1,11)
                  else          ! LO or NLO
                     ff(10,n,m) = 0.d0 ! = chm-cbar
                     ff(11,n,m) = 0.d0 ! = bot-bbar
                     ff(12,n,m) = 0.d0 ! = photon
                     read(33,'(9(1pe12.4))',iostat=io)
     &                    (ff(ip,n,m),ip=1,9)
                  end if
               end if
               if (io.ne.0) then
                  print *,"Error in mstwGetOnePDF reading ",filename
                  stop
               end if
            enddo
         enddo

C--   Check that ALL the file contents have been read in.
         read(33,*,iostat=io) dummy
         if (io.eq.0) then
            print *,"Error in mstwGetOnePDF: not at end of ",filename
            stop
         end if
         close(unit=33)

C--   PDFs are identically zero at x = 1.
         do m=1,nq
            do ip=1,np
               ff(ip,nx,m)=0d0
            enddo
         enddo

         do n=1,nx
            xxl(n)=log10(xx(n))
         enddo
         do m=1,nq
            qql(m)=log10(qq(m))
         enddo

C--   Initialise all parton flavours.
         do ip=1,np
            call mstwInitialisePDF(ip,np,ih,nhess,nx,nq,nqc0,nqb0,
     &           xxl,qql,ff,cc)
         enddo

         oldprefix(ih) = prefix

C--   ... End of initialisation for eigenvector set "ih".

      end if                    ! oldprefix(ih).ne.prefix

*--------------------------------------------------------------------

      qsq=q*q
C--   If mc2 < qsq < mc2+eps, then qsq = mc2+eps.
      if (qsq.gt.qq(nqc0).and.qsq.lt.qq(nqc0+1)) qsq = qq(nqc0+1)
C--   If mb2 < qsq < mb2+eps, then qsq = mb2+eps.
      if (qsq.gt.qq(nqb0).and.qsq.lt.qq(nqb0+1)) qsq = qq(nqb0+1)
      
      xlog=log10(x)
      qsqlog=log10(qsq)

      res = 0.d0

      if (f.eq.0) then          ! gluon
         ip = 1
      else if (f.ge.1.and.f.le.5) then ! quarks
         ip = f+1
      else if (f.le.-1.and.f.ge.-5) then ! antiquarks
         ip = -f+1
      else if (f.ge.7.and.f.le.11) then ! valence quarks
         ip = f
      else if (f.eq.13) then    ! photon
         ip = 12
      else if (abs(f).ne.6.and.f.ne.12) then
         if (warn.or.fatal) print *,"Error in mstwGetOnePDF: f = ",f
         if (fatal) stop
      end if
      
      if (x.le.0.d0.or.x.gt.xmax.or.q.le.0.d0) then

         if (warn.or.fatal) print *,"Error in mstwGetOnePDF: x,qsq = ",
     &        x,qsq
         if (fatal) stop

      else if (abs(f).eq.6.or.f.eq.12) then ! set top quarks to zero
         
         res = 0.d0

      else if (qsq.lt.qsqmin) then ! extrapolate to low Q^2

         if (warn) then
            print *, "Warning in mstwGetOnePDF, extrapolating: f = ",f,
     &           ", x = ",x,", q = ",q
         end if

         if (x.lt.xmin) then    ! extrapolate to low x

            res = mstwExtrapolatePDF(ip,np,ih,nhess,xlog,
     &           log10(qsqmin),nx,nq,xxl,qql,cc)
            res1 = mstwExtrapolatePDF(ip,np,ih,nhess,xlog,
     &           log10(1.01D0*qsqmin),nx,nq,xxl,qql,cc)
            if (f.le.-1.and.f.ge.-5) then ! antiquark = quark - valence
               res = res - mstwExtrapolatePDF(ip+5,np,ih,nhess,xlog,
     &              log10(qsqmin),nx,nq,xxl,qql,cc)
               res1 = res1 - mstwExtrapolatePDF(ip+5,np,ih,nhess,xlog,
     &              log10(1.01D0*qsqmin),nx,nq,xxl,qql,cc)
            end if
            
         else                   ! do usual interpolation
            
            res = mstwInterpolatePDF(ip,np,ih,nhess,xlog,
     &           log10(qsqmin),nx,nq,xxl,qql,cc)
            res1 = mstwInterpolatePDF(ip,np,ih,nhess,xlog,
     &           log10(1.01D0*qsqmin),nx,nq,xxl,qql,cc)
            if (f.le.-1.and.f.ge.-5) then ! antiquark = quark - valence
               res = res - mstwInterpolatePDF(ip+5,np,ih,nhess,xlog,
     &              log10(qsqmin),nx,nq,xxl,qql,cc)
               res1 = res1 - mstwInterpolatePDF(ip+5,np,ih,nhess,xlog,
     &              log10(1.01D0*qsqmin),nx,nq,xxl,qql,cc)
            end if
            
         end if

C--   Calculate the anomalous dimension, dlog(xf)/dlog(qsq),
C--   evaluated at qsqmin.  Then extrapolate the PDFs to low
C--   qsq < qsqmin by interpolating the anomalous dimenion between
C--   the value at qsqmin and a value of 1 for qsq << qsqmin.
C--   If value of PDF at qsqmin is very small, just set
C--   anomalous dimension to 1 to prevent rounding errors.
C--   Impose minimum anomalous dimension of -2.5.
         if (abs(res).ge.1.D-5) then
            anom = max( -2.5D0, (res1-res)/res/0.01D0 )
         else
            anom = 1.D0
         end if
         res = res*(qsq/qsqmin)**(anom*qsq/qsqmin+1.D0-qsq/qsqmin)

      else if (x.lt.xmin.or.qsq.gt.qsqmax) then ! extrapolate

         if (warn) then
            print *, "Warning in mstwGetOnePDF, extrapolating: f = ",f,
     &           ", x = ",x,", q = ",q
         end if

         res = mstwExtrapolatePDF(ip,np,ih,nhess,xlog,
     &        qsqlog,nx,nq,xxl,qql,cc)
         
         if (f.le.-1.and.f.ge.-5) then ! antiquark = quark - valence
            res = res - mstwExtrapolatePDF(ip+5,np,ih,nhess,xlog,
     &           qsqlog,nx,nq,xxl,qql,cc)
         end if

      else                      ! do usual interpolation
         
         res = mstwInterpolatePDF(ip,np,ih,nhess,xlog,
     &        qsqlog,nx,nq,xxl,qql,cc)

         if (f.le.-1.and.f.ge.-5) then ! antiquark = quark - valence
            res = res - mstwInterpolatePDF(ip+5,np,ih,nhess,xlog,
     &           qsqlog,nx,nq,xxl,qql,cc)
         end if
            
      end if
      
      mstwGetOnePDF = res

      return
      end


********************************************************************************
********************************************************************************

      subroutine mstwInitialisePDF(ip,np,ih,nhess,nx,my,myc0,myb0,
     &     xx,yy,ff,cc)

      implicit none

      integer nhess,ih,nx,my,myc0,myb0,j,k,l,m,n,ip,np
      double precision xx(nx),yy(my),ff(np,nx,my),
     &     ff1(nx,my),ff2(nx,my),ff12(nx,my),ff21(nx,my),
     &     yy0(4),yy1(4),yy2(4),yy12(4),z(16),
     &     cl(16),cc(np,0:nhess,nx,my,4,4),iwt(16,16),
     &     mstw_polderiv1,mstw_polderiv2,mstw_polderiv3,d1,d2,d1d2,xxd

      data iwt/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &     0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
     &     -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
     &     2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,
     &     0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
     &     0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
     &     0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
     &     0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,
     &     -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
     &     0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
     &     9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
     &     -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,
     &     2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
     &     0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
     &     -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
     &     4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1/

      do m=1,my
         ff1(1,m)=mstw_polderiv1(xx(1),xx(2),xx(3),
     &        ff(ip,1,m),ff(ip,2,m),ff(ip,3,m))
         ff1(nx,m)=mstw_polderiv3(xx(nx-2),xx(nx-1),xx(nx),
     &        ff(ip,nx-2,m),ff(ip,nx-1,m),ff(ip,nx,m))
         do n=2,nx-1
            ff1(n,m)=mstw_polderiv2(xx(n-1),xx(n),xx(n+1),
     &           ff(ip,n-1,m),ff(ip,n,m),ff(ip,n+1,m))
         enddo
      enddo

C--   Calculate the derivatives at qsq=mc2,mc2+eps,mb2,mb2+eps
C--   in a similar way as at the endpoints qsqmin and qsqmax.
      do n=1,nx
         do m=1,my
            if (myc0.eq.2.and.m.eq.1) then
               ff2(n,m)=(ff(ip,n,m+1)-ff(ip,n,m))/(yy(m+1)-yy(m))
            else if (myc0.eq.2.and.m.eq.2) then
               ff2(n,m)=(ff(ip,n,m)-ff(ip,n,m-1))/(yy(m)-yy(m-1))
            else if (m.eq.1.or.m.eq.myc0+1.or.m.eq.myb0+1) then
               ff2(n,m)=mstw_polderiv1(yy(m),yy(m+1),yy(m+2),
     &              ff(ip,n,m),ff(ip,n,m+1),ff(ip,n,m+2))
            else if (m.eq.my.or.m.eq.myc0.or.m.eq.myb0) then
               ff2(n,m)=mstw_polderiv3(yy(m-2),yy(m-1),yy(m),
     &              ff(ip,n,m-2),ff(ip,n,m-1),ff(ip,n,m))
            else
               ff2(n,m)=mstw_polderiv2(yy(m-1),yy(m),yy(m+1),
     &              ff(ip,n,m-1),ff(ip,n,m),ff(ip,n,m+1))
            end if
         end do
      end do

C--   Calculate the cross derivatives (d/dx)(d/dy).
      do m=1,my
         ff12(1,m)=mstw_polderiv1(xx(1),xx(2),xx(3),
     &        ff2(1,m),ff2(2,m),ff2(3,m))
         ff12(nx,m)=mstw_polderiv3(xx(nx-2),xx(nx-1),xx(nx),
     &        ff2(nx-2,m),ff2(nx-1,m),ff2(nx,m))
         do n=2,nx-1
            ff12(n,m)=mstw_polderiv2(xx(n-1),xx(n),xx(n+1),
     &           ff2(n-1,m),ff2(n,m),ff2(n+1,m))
         enddo
      enddo

C--   Calculate the cross derivatives (d/dy)(d/dx).
      do n=1,nx
         do m = 1, my
            if (myc0.eq.2.and.m.eq.1) then
               ff21(n,m)=(ff1(n,m+1)-ff1(n,m))/(yy(m+1)-yy(m))
            else if (myc0.eq.2.and.m.eq.2) then
               ff21(n,m)=(ff1(n,m)-ff1(n,m-1))/(yy(m)-yy(m-1))
            else if (m.eq.1.or.m.eq.myc0+1.or.m.eq.myb0+1) then
               ff21(n,m)=mstw_polderiv1(yy(m),yy(m+1),yy(m+2),
     &              ff1(n,m),ff1(n,m+1),ff1(n,m+2))
            else if (m.eq.my.or.m.eq.myc0.or.m.eq.myb0) then
               ff21(n,m)=mstw_polderiv3(yy(m-2),yy(m-1),yy(m),
     &              ff1(n,m-2),ff1(n,m-1),ff1(n,m))
            else
               ff21(n,m)=mstw_polderiv2(yy(m-1),yy(m),yy(m+1),
     &              ff1(n,m-1),ff1(n,m),ff1(n,m+1))
            end if
         end do
      end do

C--   Take the average of (d/dx)(d/dy) and (d/dy)(d/dx).
      do n=1,nx
         do m = 1, my
            ff12(n,m)=0.5*(ff12(n,m)+ff21(n,m))
         end do
      end do

      do n=1,nx-1
         do m=1,my-1
            d1=xx(n+1)-xx(n)
            d2=yy(m+1)-yy(m)
            d1d2=d1*d2
            
            yy0(1)=ff(ip,n,m)
            yy0(2)=ff(ip,n+1,m)
            yy0(3)=ff(ip,n+1,m+1)
            yy0(4)=ff(ip,n,m+1)
            
            yy1(1)=ff1(n,m)
            yy1(2)=ff1(n+1,m)
            yy1(3)=ff1(n+1,m+1)
            yy1(4)=ff1(n,m+1)
            
            yy2(1)=ff2(n,m)
            yy2(2)=ff2(n+1,m)
            yy2(3)=ff2(n+1,m+1)
            yy2(4)=ff2(n,m+1)
            
            yy12(1)=ff12(n,m)
            yy12(2)=ff12(n+1,m)
            yy12(3)=ff12(n+1,m+1)
            yy12(4)=ff12(n,m+1)
            
            do k=1,4
               z(k)=yy0(k)
               z(k+4)=yy1(k)*d1
               z(k+8)=yy2(k)*d2
               z(k+12)=yy12(k)*d1d2
            enddo
            
            do l=1,16
               xxd=0.d0
               do k=1,16
                  xxd=xxd+iwt(k,l)*z(k)
               enddo
               cl(l)=xxd
            enddo
            l=0
            do k=1,4
               do j=1,4
                  l=l+1
                  cc(ip,ih,n,m,k,j)=cl(l)
               enddo
            enddo
         enddo
      enddo
      return

      end


********************************************************************************
********************************************************************************

      double precision function mstwInterpolatePDF(ip,np,ih,nhess,x,y,
     &     nx,my,xx,yy,cc)

      implicit none

      integer ih,nx,my,nhess,mstw_locx8,l,m,n,ip,np
      double precision xx(nx),yy(my),cc(np,0:nhess,nx,my,4,4),
     &     x,y,z,t,u

      n=mstw_locx8(xx,nx,x)
      m=mstw_locx8(yy,my,y)
      
      t=(x-xx(n))/(xx(n+1)-xx(n))
      u=(y-yy(m))/(yy(m+1)-yy(m))
      
      z=0.d0
      do l=4,1,-1
         z=t*z+((cc(ip,ih,n,m,l,4)*u+cc(ip,ih,n,m,l,3))*u
     &        +cc(ip,ih,n,m,l,2))*u+cc(ip,ih,n,m,l,1)
      enddo

      mstwInterpolatePDF = z

      return

      end


********************************************************************************
********************************************************************************

      double precision function mstwExtrapolatePDF(ip,np,ih,nhess,x,y,
     &     nx,my,xx,yy,cc)

      implicit none

      integer ih,nx,my,nhess,mstw_locx8,n,m,ip,np
      double precision xx(nx),yy(my),cc(np,0:nhess,nx,my,4,4),
     &     x,y,z,f0,f1,z0,z1,mstwInterpolatePDF
      
      n=mstw_locx8(xx,nx,x)           ! 0: below xmin, nx: above xmax
      m=mstw_locx8(yy,my,y)           ! 0: below qsqmin, my: above qsqmax
      
C--   If extrapolation in small x only:
      if (n.eq.0.and.m.gt.0.and.m.lt.my) then
         f0 = mstwInterpolatePDF(ip,np,ih,nhess,xx(1),y,nx,my,xx,yy,cc)
         f1 = mstwInterpolatePDF(ip,np,ih,nhess,xx(2),y,nx,my,xx,yy,cc)
         if (f0.gt.1.d-3.and.f1.gt.1.d-3) then
            z = exp(log(f0)+(log(f1)-log(f0))/(xx(2)-xx(1))*(x-xx(1)))
         else
            z = f0+(f1-f0)/(xx(2)-xx(1))*(x-xx(1))
         end if
C--   If extrapolation into large q only:
      else if (n.gt.0.and.m.eq.my) then
         f0 = mstwInterpolatePDF(ip,np,ih,nhess,x,yy(my),nx,my,xx,yy,cc)
         f1 = mstwInterpolatePDF(ip,np,ih,nhess,x,yy(my-1),nx,my,xx,yy,cc)
         if (f0.gt.1.d-3.and.f1.gt.1.d-3) then
            z = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*
     &           (y-yy(my)))
         else
            z = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
         end if
C--   If extrapolation into large q AND small x:
      else if (n.eq.0.and.m.eq.my) then
         f0 = mstwInterpolatePDF(ip,np,ih,nhess,xx(1),yy(my),nx,my,xx,yy,cc)
         f1 = mstwInterpolatePDF(ip,np,ih,nhess,xx(1),yy(my-1),nx,my,xx,yy,
     &        cc)
         if (f0.gt.1.d-3.and.f1.gt.1.d-3) then
            z0 = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*
     &           (y-yy(my)))
         else
            z0 = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
         end if
         f0 = mstwInterpolatePDF(ip,np,ih,nhess,xx(2),yy(my),nx,my,xx,yy,cc)
         f1 = mstwInterpolatePDF(ip,np,ih,nhess,xx(2),yy(my-1),nx,my,xx,yy,
     &        cc)
         if (f0.gt.1.d-3.and.f1.gt.1.d-3) then
            z1 = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*
     &           (y-yy(my)))
         else
            z1 = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
         end if
         if (z0.gt.1.d-3.and.z1.gt.1.d-3) then
            z = exp(log(z0)+(log(z1)-log(z0))/(xx(2)-xx(1))*(x-xx(1)))
         else
            z = z0+(z1-z0)/(xx(2)-xx(1))*(x-xx(1))
         end if
      else
         print *,"Error in mstwExtrapolatePDF"
         stop
      end if

      mstwExtrapolatePDF = z      

      return

      end


********************************************************************************
********************************************************************************

      integer function mstw_locx8(xx,nx,x)

C--   returns an integer j such that x lies inbetween xx(j) and xx(j+1).
C--   nx is the length of the array with xx(nx) the highest element.

      implicit none

      integer nx,jl,ju,jm
      double precision x,xx(nx)

      if(x.eq.xx(1)) then
         mstw_locx8=1
         return
      endif
      if(x.eq.xx(nx)) then
         mstw_locx8=nx-1  
         return
      endif
      ju=nx+1
      jl=0
    1 if((ju-jl).le.1) go to 2
      jm=(ju+jl)/2
      if(x.ge.xx(jm)) then
         jl=jm
      else
         ju=jm
      endif
      go to 1
    2 mstw_locx8=jl
      return

      end


********************************************************************************
********************************************************************************

      double precision function mstw_polderiv1(x1,x2,x3,y1,y2,y3)

C--   returns the estimate of the derivative at x1 obtained by a
C--   polynomial interpolation using the three points (x_i,y_i).

      implicit none

      double precision x1,x2,x3,y1,y2,y3

      mstw_polderiv1=(x3*x3*(y1-y2)+2.d0*x1*(x3*(-y1+y2)+x2*(y1-y3))
     &     +x2*x2*(-y1+y3)+x1*x1*(-y2+y3))/((x1-x2)*(x1-x3)*(x2-x3))

      return

      end


********************************************************************************
********************************************************************************

      double precision function mstw_polderiv2(x1,x2,x3,y1,y2,y3)

C--   returns the estimate of the derivative at x2 obtained by a
C--   polynomial interpolation using the three points (x_i,y_i).

      implicit none

      double precision x1,x2,x3,y1,y2,y3

      mstw_polderiv2=(x3*x3*(y1-y2)-2.d0*x2*(x3*(y1-y2)+x1*(y2-y3))
     &     +x2*x2*(y1-y3)+x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3))

      return

      end


********************************************************************************
********************************************************************************

      double precision function mstw_polderiv3(x1,x2,x3,y1,y2,y3)

C--   returns the estimate of the derivative at x3 obtained by a
C--   polynomial interpolation using the three points (x_i,y_i).

      implicit none

      double precision x1,x2,x3,y1,y2,y3

      mstw_polderiv3=(x3*x3*(-y1+y2)+2.d0*x2*x3*(y1-y3)+x1*x1*(y2-y3)
     &     +x2*x2*(-y1+y3)+2.d0*x1*x3*(-y2+y3))/
     &     ((x1-x2)*(x1-x3)*(x2-x3))

      return

      end

********************************************************************************
********************************************************************************

C--   G.W. 05/07/2010 Avoid using Fortran 90 intrinsic function
C--   "len_trim" since not supported by all Fortran 77 compilers.

      integer function mstw_lentrim(s)

      implicit none
      character*(*) s

      do mstw_lentrim = len(s), 1, -1
         if (s(mstw_lentrim:mstw_lentrim) .ne. ' ') return
      end do
      return

      end


********************************************************************************
