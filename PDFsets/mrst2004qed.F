*sophy: 05-01-2011
** Altered such that the photon PDF is separated from the quark/gluon
*sophy: 14-10-2009
** New file by Terrance Figy: for mrst2004qed PDF (needed for photon induced 
**  processes.
** Altered open statements in subroutine mrst1 and mrst2
C--   Alternative LHAPDF-like interface: return PDFs in an array.
      subroutine mrst2004qed(mode,x,q,xpdf,GAMpdf)
      implicit none
      integer mode
      real*8 x, q, xpdf(-6:6), GAMpdf
      real*8 upv,dnv,usea,dsea,str,chm,bot,glu,phot

      call mrstqed(x,q,mode,
     &     upv,dnv,usea,dsea,str,chm,bot,glu,phot)

      xpdf(-6) = 0              ! tbar
      xpdf(-5) = bot            ! bbar
      xpdf(-4) = chm            ! cbar
      xpdf(-3) = str            ! sbar
      xpdf(-2) = usea           ! ubar
      xpdf(-1) = dsea           ! dbar
      xpdf(0) = glu             ! g
      xpdf(1) = dnv+dsea        ! d
      xpdf(2) = upv+usea        ! u
      xpdf(3) = str             ! s
      xpdf(4) = chm             ! c
      xpdf(5) = bot             ! b
      xpdf(6) = 0               ! t

      GAMpdf = phot   ! photon


c      print*,'photon=',phot
      end  

      subroutine mrstqed(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,
     x    glu,phot)
C***************************************************************C
C								C
C  This is a package for the new MRST 2004-QED parton           C
C  distributions.                                               C     
C  Reference: A.D. Martin, R.G. Roberts, W.J. Stirling and      C
C  R.S. Thorne, hep-ph/0411040                                  C
C                                                               C
C  There are 2 pdf sets:  proton distributions (mode = 1)       C
C                        neutron distributions (mode = 2)       C
C                                                               C
C  These sets read grids qed6-10gridp.dat and qed6-10gridn.dat  C
C  respectively.                                                C
C  Note that the neutron grid was replaced in January 2007 to   C
C  correct a bug in the original version that affected the      C
C  ubar and dbar distributions at high x.                       C
C							        C
C								C
C  Note the extra argument "phot" for the photon distribution   C
C  in the proton/neutron. As always, the quantity returned      C
C  is x (f,x,Q^2).                                              C
C								C
C   This subroutine uses an improved interpolation procedure    C 
C   for extracting values of the pdf's from the grid            C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C***************************************************************C
      implicit none
      double precision x, q, q2, qsqmax, qsqmin, xmin, xmax
      double precision upv, dnv, usea, dsea, str, chm, bot, glu, phot
      integer mode
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      q2=q*q
cMR: Remove excessive warnings
c      if(q2.lt.qsqmin.or.q2.gt.qsqmax) print 99,q2
c      if(x.lt.xmin.or.x.gt.xmax)       print 98,x
          if(mode.eq.1) then
        call mrst1(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu,phot) 
      elseif(mode.eq.2) then
        call mrst2(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu,phot) 
      endif 
      end

      subroutine mrst1(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
      implicit none
      double precision upv,dnv,usea,dsea,str,chm,bot,glu,phot
      double precision qsq, q2save, qsqlog, qsqmax, qsqmin
      double precision x, xlog, xmax, xmin, xsave
      double precision emb2, emc2
      integer nx, nq, nqb, nqc, init, m, n, nqb0, nqc0, np
      parameter(nx=49,nq=37,np=9,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),f9(nx,nq),fc(nx,nqc),
     .fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .cc9(nx,nq,4,4),ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)

*sophy: declarations of path variables
      character*250 GetPdfsetsPath
      character*250 path
      character datafile*40
      external GetPdfsetsPath

      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
*sophy: Initialize statement
        print*,"Initialize PDF set mrst2004qed"

*sophy: changed open statement
        path=GetPdfsetsPath()
        datafile = 'qed6-10gridp.dat'
        path=trim(path)//"/"//datafile
        open(unit=33,file=path,status='old')

        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m),f9(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea 9=photon
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
      f9(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppeV1(nx,nq,xxl,qql,f1,cc1)
      call jeppeV1(nx,nq,xxl,qql,f2,cc2)
      call jeppeV1(nx,nq,xxl,qql,f3,cc3)
      call jeppeV1(nx,nq,xxl,qql,f4,cc4)
      call jeppeV1(nx,nq,xxl,qql,f6,cc6)
      call jeppeV1(nx,nq,xxl,qql,f8,cc8)
      call jeppeV1(nx,nq,xxl,qql,f9,cc9)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppeV1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppeV1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppeV2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppeV2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppeV2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppeV2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppeV2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppeV2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)
      call jeppeV2(xlog,qsqlog,nx,nq,xxl,qql,cc9,phot)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppeV2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppeV2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(9f10.5)
      end

      subroutine mrst2(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
      implicit none
      double precision upv,dnv,usea,dsea,str,chm,bot,glu,phot, emb2,emc2
      double precision qsq, q2save, qsqlog, qsqmin, qsqmax
      double precision x, xlog, xmin, xmax, xsave
      integer nx, nq, nqb, nqc, init, m, n, nqb0, nqc0, np
*      Implicit Double Precision (A-H,O-Z)
      parameter(nx=49,nq=37,np=9,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),f9(nx,nq),fc(nx,nqc),
     .fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .cc9(nx,nq,4,4),ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)

*sophy: declarations of path variables
      character*250 GetPdfsetsPath
      character*250 path
      character datafile*40
      external GetPdfsetsPath

      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10

*sophy: Initialize statement
        print*,"Initialize PDF set mrst2004qed"

*sophy: changed open statement
        path=GetPdfsetsPath()
        datafile = 'qed6-10gridn.dat'
        path=trim(path)//"/"//datafile
        open(unit=34,file=path,status='old')

        do 20 n=1,nx-1
        do 20 m=1,nq
        read(34,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m),f9(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea 9=photon
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
      f9(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppeV1(nx,nq,xxl,qql,f1,cc1)
      call jeppeV1(nx,nq,xxl,qql,f2,cc2)
      call jeppeV1(nx,nq,xxl,qql,f3,cc3)
      call jeppeV1(nx,nq,xxl,qql,f4,cc4)
      call jeppeV1(nx,nq,xxl,qql,f6,cc6)
      call jeppeV1(nx,nq,xxl,qql,f8,cc8)
      call jeppeV1(nx,nq,xxl,qql,f9,cc9)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppeV1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppeV1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppeV2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppeV2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppeV2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppeV2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppeV2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppeV2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)
      call jeppeV2(xlog,qsqlog,nx,nq,xxl,qql,cc9,phot)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppeV2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppeV2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(9f10.5)
      end

      subroutine jeppeV1(nx,my,xx,yy,ff,cc)
      implicit none
      integer nx, my, iwt, j, k, l, m, n
      double precision xx, yy, cc, cl, d1, d2, d1d2, dx, dy, z
      double precision ff, ff1, ff2, ff12, xxd, yy0, yy1, yy2, yy12
      double precision polderiv
      dimension xx(nx),yy(my),ff(nx,my),ff1(nx,my),ff2(nx,my),
     -     ff12(nx,my),yy0(4),yy1(4),yy2(4),yy12(4),z(16),
     -     cl(16),cc(nx,my,4,4),iwt(16,16)

      data iwt/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
     x		  -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
     x		  2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,
     x		  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
     x		  0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
     x		  0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,
     x		  -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
     x		  9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
     x		  -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,
     x		  2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
     x		  -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
     x		  4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1/


      do 42 m=1,my
      dx=xx(2)-xx(1)
      ff1(1,m)=(ff(2,m)-ff(1,m))/dx
      dx=xx(nx)-xx(nx-1)
      ff1(nx,m)=(ff(nx,m)-ff(nx-1,m))/dx
      do 41 n=2,nx-1
      ff1(n,m)=polderiv(xx(n-1),xx(n),xx(n+1),ff(n-1,m),ff(n,m),
     xff(n+1,m))
   41 continue
   42 continue

      do 44 n=1,nx
      dy=yy(2)-yy(1)
      ff2(n,1)=(ff(n,2)-ff(n,1))/dy
      dy=yy(my)-yy(my-1)
      ff2(n,my)=(ff(n,my)-ff(n,my-1))/dy
      do 43 m=2,my-1
      ff2(n,m)=polderiv(yy(m-1),yy(m),yy(m+1),ff(n,m-1),ff(n,m),
     xff(n,m+1))
   43 continue
   44 continue

      do 46 m=1,my
      dx=xx(2)-xx(1)
      ff12(1,m)=(ff2(2,m)-ff2(1,m))/dx
      dx=xx(nx)-xx(nx-1)
      ff12(nx,m)=(ff2(nx,m)-ff2(nx-1,m))/dx
      do 45 n=2,nx-1
      ff12(n,m)=polderiv(xx(n-1),xx(n),xx(n+1),ff2(n-1,m),ff2(n,m),
     xff2(n+1,m))
   45 continue
   46 continue

      do 53 n=1,nx-1
      do 52 m=1,my-1
      d1=xx(n+1)-xx(n)
      d2=yy(m+1)-yy(m)
      d1d2=d1*d2

      yy0(1)=ff(n,m)
      yy0(2)=ff(n+1,m)
      yy0(3)=ff(n+1,m+1)
      yy0(4)=ff(n,m+1)

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

      do 47 k=1,4
      z(k)=yy0(k)
      z(k+4)=yy1(k)*d1
      z(k+8)=yy2(k)*d2
      z(k+12)=yy12(k)*d1d2
   47 continue

      do 49 l=1,16
      xxd=0.
      do 48 k=1,16
      xxd=xxd+iwt(k,l)*z(k)
   48 continue
      cl(l)=xxd
   49 continue
      l=0
      do 51 k=1,4
      do 50 j=1,4
      l=l+1
      cc(n,m,k,j)=cl(l)
   50 continue
   51 continue
   52 continue
   53 continue
      return
      end

      subroutine jeppeV2(x,y,nx,my,xx,yy,cc,z)
C--   G.W. 02/07/2007 Allow extrapolation to small x and large q.
      implicit none
      double precision x, y, xx, yy, cc, z, t, u, z0, z1
      integer locx
      double precision f0, f1
      integer nx, my, l, m, n
      dimension xx(nx),yy(my),cc(nx,my,4,4)      
      
      n=locx(xx,nx,x)
      m=locx(yy,my,y)
      
      if (n.gt.0.and.n.lt.nx.and.m.gt.0.and.m.lt.my) then
C--   Do usual interpolation.
         t=(x-xx(n))/(xx(n+1)-xx(n))
         u=(y-yy(m))/(yy(m+1)-yy(m))
         z=0.d0
         do l=4,1,-1
            z=t*z+((cc(n,m,l,4)*u+cc(n,m,l,3))*u
     &           +cc(n,m,l,2))*u+cc(n,m,l,1)
         enddo
         
      else if (n.eq.0.and.m.gt.0.and.m.lt.my) then
C--   Extrapolate to small x.
         call jeppeV3(xx(1),y,nx,my,xx,yy,cc,f0)
         call jeppeV3(xx(2),y,nx,my,xx,yy,cc,f1)
         if (f0.gt.0.d0.and.f1.gt.0.d0) then
            z = exp(log(f0)+(log(f1)-log(f0))/(xx(2)-xx(1))*(x-xx(1)))
         else
            z = f0+(f1-f0)/(xx(2)-xx(1))*(x-xx(1))
         end if
         
      else if (n.gt.0.and.m.eq.my) then
C--   Extrapolate to large q.
         call jeppeV3(x,yy(my),nx,my,xx,yy,cc,f0)
         call jeppeV3(x,yy(my-1),nx,my,xx,yy,cc,f1)
         if (f0.gt.0.d0.and.f1.gt.0.d0) then
            z = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*
     &           (y-yy(my)))
         else
            z = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
         end if
         
      else if (n.eq.0.and.m.eq.my) then
C--   Extrapolate to small x AND large q.
         call jeppeV3(xx(1),yy(my),nx,my,xx,yy,cc,f0)
         call jeppeV3(xx(1),yy(my-1),nx,my,xx,yy,cc,f1)
         if (f0.gt.0.d0.and.f1.gt.0.d0) then
            z0 = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*
     &           (y-yy(my)))
         else
            z0 = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
         end if
         call jeppeV3(xx(2),yy(my),nx,my,xx,yy,cc,f0)
         call jeppeV3(xx(2),yy(my-1),nx,my,xx,yy,cc,f1)
         if (f0.gt.0.d0.and.f1.gt.0.d0) then
            z1 = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*
     &           (y-yy(my)))
         else
            z1 = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
         end if
         if (z0.gt.0.d0.and.z1.gt.0.d0) then
            z = exp(log(z0)+(log(z1)-log(z0))/(xx(2)-xx(1))*(x-xx(1)))
         else
            z = z0+(z1-z0)/(xx(2)-xx(1))*(x-xx(1))
         end if

      else
C--   Set parton distribution to zero otherwise.
         z = 0.d0

      end if
      
      return
      end

C--   G.W. 02/07/2007 Copy of the original jeppeV2,
C--   only used for extrapolation.
      subroutine jeppeV3(x,y,nx,my,xx,yy,cc,z)
      implicit none
      double precision x, y, xx, yy, cc, z, t, u
      integer nx, my, l, m, n
      integer locx
      dimension xx(nx),yy(my),cc(nx,my,4,4)      
      n=locx(xx,nx,x)
      m=locx(yy,my,y)
      t=(x-xx(n))/(xx(n+1)-xx(n))
      u=(y-yy(m))/(yy(m+1)-yy(m))
      z=0.d0
      do l=4,1,-1
         z=t*z+((cc(n,m,l,4)*u+cc(n,m,l,3))*u
     &        +cc(n,m,l,2))*u+cc(n,m,l,1)
      enddo
      return
      end
      
      integer function locx(xx,nx,x)
      implicit none
      double precision xx, x
      integer nx, jl, jm, ju
      dimension xx(nx)
c$$$      if(x.le.xx(1)) then
      if(x.eq.xx(1)) then ! G.W. 02/07/2007
      locx=1
      return
      endif
c$$$      if(x.ge.xx(nx)) then 
      if(x.eq.xx(nx)) then ! G.W. 02/07/2007
      locx=nx-1  
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
    2 locx=jl
      return
      end


      double precision function polderiv(x1,x2,x3,y1,y2,y3)
      implicit none
      double precision x1, x2, x3, y1, y2, y3 
      polderiv=(x3*x3*(y1-y2)-2.0*x2*(x3*(y1-y2)+x1*
     .(y2-y3))+x2*x2*(y1-y3)+x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3))
      return
      end
      

