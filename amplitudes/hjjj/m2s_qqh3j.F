c
c	Terrance Figy, <terrance@physics.wisc.edu>
c	Initial version:  2005 July
c	Last modified:    2008 Oct 1
c   
c     |Born|^2 + 2 Re[Born*cg(virt_finite)],   nlo = 1 old comment line
c     |Born|^2 + 2 Re[Born*cg(virt)]+<I(eps)>, nlo = 1
c     |born|^2                                 nlo = 0
c
      subroutine m2s_wbfhj(
     &                   bos,   !in:  Boson identifier, 2 = Z, 6=H
     &                   nlo,   !in:  NLO=1: brn+ virt; NLO = 0:brn
     &                   xi,	!in:  Feynman x parameters of incoming fermions
     &                   p,	!in:  momenta associated with external fermions
     &                   v,	!in:  momenta associated with WW decay fermions
     &                   rn,     !in:  random number for z integration 
     &                   m2s  )  !out: |M|^2*pdf1*pdf2 
          use globalvars, only: ldoblha
          use VBFNLO_HIGGSTO, only: VBFNLO_HIGGSTO_gridHt
      implicit none
c
c declare input/output variables
c
#include "VBFNLO/utilities/global.inc"
      integer bos
      real*8 xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin),m2s,rn
c
c declare external functions
c
      integer FL_H3j
      external FL_H3j
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
#include "VBFNLO/utilities/process.inc"
c
c declare local variables
c
c
      real*8 q_sf,alsotpi,x1,x2,J1,J2
      real*8 Ax(2),Bx(3),Cx(3),Dxi(3)
c
      integer init/0/, I, J, mu
      save init
c
c declare variables for summation process
c
      INTEGER ires

      integer FSIGN(4+max_v), gsign, physToDiag(5), 
     $     nlo, nmin, nmax,nproc(8)
      save nmin, nmax,nproc
c
c store contributions from subprocess ip in res(ip) where
c  
      real*8 res(maxnumsubproc)
      real*8 sum,dot0p
      external dot0p
c
c declare parton distribution variables
c
      real*8 pdf(-6:6,2,5)          !pdf's
c
c define program switches
c
      logical ldebug,ldebug1,linfo
      parameter (ldebug=.false.)
      parameter(ldebug1=.true.,linfo=.false.)

      real*8 weight,rnumb,RandomNumber

c CPS scheme
      DOUBLE PRECISION CLR,XM2,XMG,Bcoupl,Vcoupl,Acoupl
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),Bcoupl(6,6,6),Vcoupl(4,5),Acoupl(4,5)
      double precision qinv, qvv(0:3)
      double precision qgammaq
      common /VBFNLO_HIGGS_CPS/ qgammaq

c just a little setup for BLHA mode
      if (ldoblha) then
        physToDiag(1)=blha_physdiagmap(1,blha_cursubproc,blha_curproc)
        physToDiag(2)=blha_physdiagmap(2,blha_cursubproc,blha_curproc)
        physToDiag(3)=blha_physdiagmap(3,blha_cursubproc,blha_curproc)
        physToDiag(4)=blha_physdiagmap(4,blha_cursubproc,blha_curproc)
        physToDiag(5)=blha_physdiagmap(5,blha_cursubproc,blha_curproc)
        fsign(1) = blha_fsign(1,blha_cursubproc,blha_curproc)
        fsign(2) = blha_fsign(2,blha_cursubproc,blha_curproc)
        fsign(3) = blha_fsign(3,blha_cursubproc,blha_curproc)
        fsign(4) = blha_fsign(4,blha_cursubproc,blha_curproc)
        gsign    = blha_fsign(5,blha_cursubproc,blha_curproc)

        call wbf_h3j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,
     1    pdf,res,nmin,nmax) 

        return
      endif

c else not BLHA mode:
c
c if first iteration,
c
      if ( init .ne. 0 ) then
         do i = 1,nmax
               res(i) = 0
         enddo
      else     
         if (bos.eq.6) then
            write(6,*) " "
            write(6,*) "H3j amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         endif
         call printnfl(.true.)
         print*," "
c
         init = 1
         do i = 1,maxnumsubproc
            res(i) = 0
         enddo
      endif
      
      if(higgsscheme.eq.2) then    ! higgsscheme=2=CPS
        SELECT CASE(procid)
        CASE(Hjjj)
          do mu=0,3
            qvv(mu) = v(mu,1,1)+v(mu,2,1)+v(mu,3,1)+v(mu,4,1)
          enddo
          qinv = sqrt(qvv(0)**2-qvv(1)**2-qvv(2)**2-qvv(3)**2)
          call VBFNLO_HIGGSTO_gridHt(qinv,qgammaq)
          qgammaq = qinv*qgammaq
        CASE DEFAULT
          qgammaq = xmg(6)
        END SELECT
      else
        qgammaq = xmg(6)
      endif

c
c define scales 
c    
c
c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c 
c
c     initialize
      do j=1,5
      do i = -6,6
         pdf(i,1,j) = 0.0d0
         pdf(i,2,j) = 0.0d0
      enddo
      enddo
      
      q_sf = sqrt(mufsq(1,1))
      call pdfproton( xi(1), q_sf, pdf(-6,1,1)) !f_a/p(xi(1))

      q_sf = sqrt(mufsq(2,1))
      call pdfproton( xi(2), q_sf, pdf(-6,2,1)) !f_b/p(xi(2))
  
#ifdef WITH_NLO       
      if (nlo.eq.1) then 

         x1 = (1.0d0 - xi(1))*RN**2 + xi(1)
         x1 = xi(1)/x1
         J1 = (1.0d0 - xi(1))*2.0D0*RN
         call pdfproton(x1, q_sf, pdf(-6,1,2)) ! f_a/p(xi(1)/z)
         do i=-6,6
            do j=2,4
                  pdf(i,1,j+1) = pdf(i,1,2)
            enddo
         enddo

         x2 = (1.0d0 - xi(2))*RN**2 + xi(2)
         x2 = xi(2)/x2
         J2 = (1.0d0 - xi(2))*2.0d0*RN
     
         call pdfproton(x2, q_sf, pdf(-6,2,2)) ! f_b/p(xi(2)/z)
         do i=-6,6
            do j=2,4
               pdf(i,2,j+1) = pdf(i,2,2)
            enddo
         enddo

         alsotpi = als(1,1)/(2d0*pi)
 
         call fincollinear(Ax,Bx,Cx,Dxi,xi,rn,p,3,1) !3,1 color correlation
   
c     sum over quark pdf's
         sum = 0.0d0


         do i=1,nflVBF   !(nfl/2)*2
            sum = sum + pdf(i,1,2)+pdf(-i,1,2)
         enddo

c     quarks
         do i=1,5
            pdf(i,1,2) = alsotpi*
     $           (J1*((pdf(i,1,2)-pdf(i,1,1))*
     $           Bx(1)+pdf(i,1,2)*Cx(1)+pdf(0,1,2)*Ax(2))+
     $           pdf(i,1,1)*Dxi(1))
            pdf(-i,1,2) = alsotpi*
     $           (J1*((pdf(-i,1,2)-pdf(-i,1,1))*
     $           Bx(1)+pdf(-i,1,2)*Cx(1)+pdf(0,1,2)*Ax(2))+
     $           pdf(-i,1,1)*Dxi(1))
         enddo

c gluon 
         pdf(0,1,2) = alsotpi*
     $        (J1*((pdf(0,1,2)-pdf(0,1,1))*
     $        Bx(3)+pdf(0,1,2)*Cx(3)+sum*Ax(1))+
     $        pdf(0,1,1)*Dxi(3))
     
         call fincollinear(Ax,Bx,Cx,Dxi,xi,rn,p,5,1) !5,1 color correlation

         do i=1,5
            pdf(i,1,3) = alsotpi*
     $           (J1*((pdf(i,1,3)-pdf(i,1,1))*
     $           Bx(2)+pdf(i,1,3)*Cx(2)+pdf(0,1,3)*Ax(2))+
     $           pdf(i,1,1)*Dxi(2))
            pdf(-i,1,3) = alsotpi*
     $           (J1*((pdf(-i,1,3)-pdf(-i,1,1))*
     $           Bx(2)+pdf(-i,1,3)*Cx(2)+pdf(0,1,3)*Ax(2))+
     $           pdf(-i,1,1)*Dxi(2))
         enddo

c     gluon
         pdf(0,1,3) = alsotpi*
     $        (J1*((pdf(0,1,3)-pdf(0,1,1))*
     $        Bx(3)+pdf(0,1,3)*Cx(3)+sum*Ax(1))+
     $        pdf(0,1,1)*Dxi(3))
     
         call fincollinear(Ax,Bx,Cx,Dxi,xi,rn,p,4,2) !4,2 color correlation
     
c     sum over quark pdfs
         sum = 0.0d0
         do i=1,nflVBF   !(nfl/2)*2
            sum = sum + pdf(i,2,4)+pdf(-i,2,4)
         enddo

c     quarks       
         do i=1,5
            pdf(i,2,4) = alsotpi*
     $           (J2*((pdf(i,2,4)-pdf(i,2,1))*
     $           Bx(1)+pdf(i,2,4)*Cx(1)+pdf(0,2,4)*Ax(2))
     $           +pdf(i,2,1)*Dxi(1))
            pdf(-i,2,4) = alsotpi*
     $           (J2*((pdf(-i,2,4)-pdf(-i,2,1))*
     $           Bx(1)+pdf(-i,2,4)*Cx(1)+pdf(0,2,4)*Ax(2))+
     $           pdf(-i,2,1)*Dxi(1))
         enddo

c     gluon
         pdf(0,2,4) = alsotpi*
     $        (J2*((pdf(0,2,4)-pdf(0,2,1))*
     $        Bx(3)+pdf(0,2,4)*Cx(3)+sum*Ax(1))+
     $        pdf(0,2,1)*Dxi(3))
     
         call fincollinear(Ax,Bx,Cx,Dxi,xi,rn,p,5,2) !5,2 color correlation

         do i=1,5
            pdf(i,2,5) = alsotpi*
     $           (J2*((pdf(i,2,5)-pdf(i,2,1))*
     $           Bx(2)+pdf(i,2,5)*Cx(2)+pdf(0,2,5)*Ax(2))+
     $           pdf(i,2,1)*Dxi(2))
            pdf(-i,2,5) = alsotpi*
     $           (J2*((pdf(-i,2,5)-pdf(-i,2,1))*
     $           Bx(2)+pdf(-i,2,5)*Cx(1)+pdf(0,2,5)*Ax(2))+
     $           pdf(-i,2,1)*Dxi(2))
         enddo

c     gluon
         pdf(0,2,5) = alsotpi*
     $        (J2*((pdf(0,2,5)-pdf(0,2,1))*
     $        Bx(3)+pdf(0,2,5)*Cx(3)+sum*Ax(1))+
     $        pdf(0,2,1)*Dxi(3))

      endif !nlo.eq.1
#endif /* WITH_NLO */
c     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Now sum over the subprocesses contributing to H production

      nmax = FL_H3j(fsign,-1)            !reset counter for subprocesses to 0

C*******************  q1 q3 ---> q2 q4 g H   *********************
c
c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon
C NOTE: for call of wbf_zh4j it is important that p(*,1,*) and p(*,3,*)
c correspond to 1-2 fermion line ALWAYS, i.e physToDiag(1/2)={1,3} and 
c similarly physToDiag(3/4)={2,4} for the 3-4 fermion line
      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1
      gsign = 1

      call wbf_h3j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,
     1     pdf,res,nmin,nmax) 
c
      if(linfo) print*,"nmin,nmax",nmin,nmax 
      if (init.eq.1) nproc(1) = nmax      

C*******************  q1 qb4 ---> q2 qb3 g H   **********************

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon

      fsign(3) = -1
      fsign(4) = -1

      call wbf_h3j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,
     1              pdf,res,nmin,nmax)
      if(linfo) print*,"nmin,nmax",nmin,nmax 
      if (init.eq.1) nproc(2) = nmax
     
C*******************  qbar2 q3 ---> qbar1 q4 g H   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = 1
      
      call wbf_h3j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,
     1              pdf,res,nmin,nmax)
      if(linfo) print*,"nmin,nmax",nmin,nmax 
      if (init.eq.1) nproc(3) = nmax

C*******************  qbar2 qb4 ---> qbar1 qb3 g H   ******************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1 
      fsign(4) = -1
         
      call wbf_h3j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,
     1              pdf,res,nmin,nmax)
      if(linfo) print*,"nmin,nmax",nmin,nmax 
      if (init.eq.1) nproc(4) = nmax

C*******************  g q3 ---> qb1 q2 q4 H   **********************

      physToDiag(1)=5
      physToDiag(2)=3
      physToDiag(3)=2
      physToDiag(4)=4
      physToDiag(5)=1

      fsign(1) = -1
      fsign(2) =  1
      fsign(3) =  1
      fsign(4) =  1
      gsign = -1

      call wbf_h3j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,
     1              pdf,res,nmin,nmax)
      if(linfo) print*,"nmin,nmax",nmin,nmax 
      if (init.eq.1) nproc(5) = nmax
 
C*******************  g qbar4 ---> qbar3 qb1 q2 H   **********************

      physToDiag(1)=5
      physToDiag(2)=4
      physToDiag(3)=2
      physToDiag(4)=3
      physToDiag(5)=1

      fsign(1) = -1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) = -1

      call wbf_h3j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,
     1              pdf,res,nmin,nmax)
      if(linfo) print*,"nmin,nmax",nmin,nmax        
      if (init.eq.1) nproc(6) = nmax
 
C*******************  q1 g ---> q2 qb3 q4 H   **********************

      physToDiag(1)=1
      physToDiag(2)=5
      physToDiag(3)=2
      physToDiag(4)=4
      physToDiag(5)=3

      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      gsign = -1

      call wbf_h3j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,
     1              pdf,res,nmin,nmax)
      if(linfo) print*,"nmin,nmax",nmin,nmax  
      if (init.eq.1) nproc(7) = nmax
      
C*******************  qb2 g ---> qb1 qb3 q4 H   **********************

      physToDiag(1)=2
      physToDiag(2)=5
      physToDiag(3)=1
      physToDiag(4)=4
      physToDiag(5)=3

      fsign(1) =  -1
      fsign(2) =  -1
      fsign(3) = -1
      fsign(4) =  1
      gsign = -1          

      call wbf_h3j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,
     1              pdf,res,nmin,nmax)
      if(linfo) print*,"nmin,nmax",nmin,nmax 
      if (init.eq.1) nproc(8) = nmax
c*****************  end of process evaluation  ********************
c
      if (init.eq.1) then
         init = init+1
         if (lwarn) print 199," proc #s for Hjjj are ",nproc
 199     format(a,8i5)
      endif

      m2s = 0d0
      DO IRES = 1,NMAX
         m2s = m2s + RES(IRES)
      ENDDO
c
c... Les Houches interface - the most probable subprocess 3jets at LO
      if ((lha.or.hepmc).and..not.doNLO) then
         i=0
         weight=0.d0
         rnumb=RandomNumber()
         do while((i.le.nmax).and.(weight.le.rnumb*m2s))
            i=i+1
            weight=weight+res(i)
            iprocess=i
         enddo
      endif  
c... Les Houches interface       
c 
      RETURN
      END


c******************************************************************
c
c   begin subroutine wbf_h3j
c
c*****************************************************************
      subroutine wbf_h3j(xi,p,v,physToDiag,fsign,gsign,bos,nlo,
     1     pdf,res,nmin,nmax)
          use globalvars, only: ldoblha
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
      real*8 p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), xi(nx)

      real*8 pdf(-6:6,2,5), res(maxnumsubproc)
      integer physToDiag(5), fsign(4+max_v), gsign, 
     $     bos, nlo, nmin, nmax
c
      real*8 uuccv(3),uussv(3),ddccv(3),ddssv(3),
     $     udscv(3),ducsv(3)          !finite virtual 
      real*8 uuccb(3),uussb(3),ddccb(3),ddssb(3),
     $     udscb(3),ducsb(3)          !born 
      real*8 uucc(3),uuss(3),ddcc(3),ddss(3),
     $     udsc(3),ducs(3)          
c wbf_h3j calls the amplitude square routines 
c             qqhqq    for qq-->qqgH     for bos=6
c for the subtraction terms for the NLO cross section calculation
c
c  INPUT:  p(0:3,np,mxd)      external physical parton momenta
c          v(0:3,nv,mxd)      for Higgs production only
c                          the sum q(mu) = v(mu,1)+...+v(mu,nv) is needed
c          physToDiag(5)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign,gsign     sign factors for amplitude calls; see qqZqq(j)
c          nlo             nlo = 0: calculate LO only (born)
c                          nlo = 1: calculate full NLO (born+virt) 
c  OUTPUT:
c  res(k)            
c                            
c  In and Output
c    nmin, nmax            range of process numbers for this call
c
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
#include "color.inc"
c
      real*8 pbar(0:3,4+max_v),qbar(0:4)
      real*8 polcolgq,polcolqq
      double precision NCmatrix(0:1,0:1,3)
      double precision CCmatrix(0:1,3)
      double precision NCbornmsq(0:1,0:1,3)
      double precision CCbornmsq(0:1,3)
      double precision colcorr(5),born(3)
c     variables for subtractions
      real*8 lns21,lnt21,lnu21,lns43,lnt43,lnu43
      real*8 s21,t21,u21,s43,t43,u43,asotpi
      real*8 musq,l21,l43,fac(3)
c     variables for convolution of pdf's
      integer iflav(5),diagToPhys(5),FL_H3j
      real*8 dot0p
      external FL_H3j
      external dot0p
      logical ChargedCurrent, sametype, oneAntiparticle
      logical ldebug
      parameter (ldebug=.false.)
      integer i,if1,if2,icc1,j,k,mu,jc   
      
! note about nflav: number of flavours in dipole should always be 5 as that's
! what we use for alfas running.  THe number of external PDFs isn't what we're
! concerned with here
      double precision nflav
      parameter (nflav = 5d0)  


      polcolqq = 1d0/(NC**2 * 4d0 * xi(1)*xi(2)) 
      polcolgq = 1d0/(NC * (NC**2-1d0) * 4d0 * xi(1)*xi(2))
c     ave over polarizations and colors
      if(ldebug) then
         print*,"momentum"
         do i = 1,n_p
            print*,"p(mu,i,1)=",(p(mu,i,1),mu=0,3)
         enddo
      endif

      do i=1,3
         uucc(i) = 0d0
         uuss(i) = 0d0
         ddcc(i) = 0d0
         ddss(i) = 0d0
         udsc(i) = 0d0
         ducs(i) = 0d0
      enddo
      
      nmin = nmax+1
      do i = 1,5
         diagToPhys(physToDiag(i)) = i
      enddo
         
c get the born + virt amplitude squared, 
      do mu = 0,3
         do i = 1,5
            pbar(mu,physToDiag(i))=p(mu,i,1)
         enddo
         qbar(mu) = pbar(mu,5)
      enddo

      do mu = 0,3
         pbar(mu,5) = 0         ! dummy momentum
         pbar(mu,6) = 0
         pbar(mu,7) = 0
         pbar(mu,8) = 0
         do i = 1,n_v
            pbar(mu,6) = pbar(mu,6) + v(mu,i,1) ! Higgs momentum
         enddo
      enddo
      fsign(5) = 0
      fsign(6) = 1
      fsign(7) = 0
      fsign(8) = 0      

c     get born + virt_finite  
      call qqHqqj_c_virt(nlo,pbar,fsign,qbar,gsign,0,
     1     uuccv,uussv,ddccv,ddssv,udscv,ducsv)  

      if (ldoblha) then
        return
      endif

c     set default values
      do i=1,3
         uucc(i) = uuccv(i) 
         uuss(i) = uussv(i) 
         ddcc(i) = ddccv(i) 
         ddss(i) = ddssv(i) 
         udsc(i) = udscv(i) 
         ducs(i) = ducsv(i)
      enddo

c     get contribution from 2Re[M_v(div) M_Born] + <I(eps)>
c     if nlo = 1,2,-2
c     For nlo = -2 do not include fin.coll and  
c     2Re[M_v(div) M_Born] + <I(eps)> for this case only 
c     2Re[M_v(boxline) M_Born] is to be included
c

#ifdef WITH_NLO
      if (nlo.eq.1) then

         call qqHqqj_c_virt(0,pbar,fsign,qbar,gsign,0,
     1        uuccb,uussb,ddccb,ddssb,udscb,ducsb) 

c     compute kinematic factors
c     21 quark current
c     43 quark current 
         musq = mursq(1,1)      !scales.inc
         asotpi = als(1,1)/(2d0*pi)
         l21 = als(1,1)/(2d0*pi)
         l43 = als(2,1)/(2d0*pi)

         if (gsign.eq.-1) then
            if (fsign(1).eq.-fsign(2) .and. 
     1          fsign(3).eq.fsign(4) ) then
               l43 = 0            !  initial gluon attached to 1-2 line only
            elseif (fsign(1).eq.fsign(2) .and. 
     1              fsign(3).eq.-fsign(4) ) then
               l21 = 0            !  initial gluon attached to 3-4 line only
            endif
         endif

c     mandlestam variables
         s21=2d0*dot0p(pbar(0,1),pbar(0,2)) !the sign doesn't matter here
         t21=2d0*dot0p(pbar(0,1),qbar(0))
         u21=2d0*dot0p(pbar(0,2),qbar(0))
      
         s43=2d0*dot0p(pbar(0,3),pbar(0,4)) !the sign doesn't matter here
         t43=2d0*dot0p(pbar(0,3),qbar(0))
         u43=2d0*dot0p(pbar(0,4),qbar(0))

c     logs
         lns21 = log(abs(s21/musq))
         lnt21 = log(abs(t21/musq))
         lnu21 = log(abs(u21/musq))

         lns43 = log(abs(s43/musq))
         lnt43 = log(abs(t43/musq))
         lnu43 = log(abs(u43/musq))
     
c     below are the factors that result from
c     2 Re[born*virt_div]+<I(eps)>
! note about nflav: number of flavours in dipole should always be 5 as that's
! what we use for alfas running.  THe number of external PDFs isn't what we're
! concerned with here
         fac(1) = -asotpi*2d0/9d0 *(-25.0d0*CAc-45.0d0*CF+
!     $        8.0d0*(nfl/2)*2*Tf+3.0d0*(CAc+2.0d0*CF)*pi**2) !terms that multiply born
     $        8.0d0*nflav*Tf+3.0d0*(CAc+2.0d0*CF)*pi**2) !terms that multiply born
     $     + asotpi*(-4.0d0*pi**2/3.0d0 + 10.0d0)*CF  ! 
         fac(2) = l21*( -1d0/2d0*(CAc-2.0d0*CF)*lns21**2+
     $        3d0/2d0*(CAc-2d0*CF)*lns21+CAc/2d0*(lnu21**2+lnt21**2)+
!     $        1d0/18d0*(6d0*NFL*Tf-30d0*CAc)*(lnu21+lnt21))!multiply born_21
     $        1d0/18d0*(6d0*nflav*Tf-30d0*CAc)*(lnu21+lnt21))!multiply born_21
         fac(3) = l43*(-1d0/2d0*(CAc-2.0d0*CF)*lns43**2+
     $        3d0/2d0*(CAc-2d0*CF)*lns43+CAc/2d0*(lnu43**2+lnt43**2)+
!     $        1d0/18d0*(6d0*NFL*Tf-30d0*CAc)*(lnu43+lnt43)) !multiply born_43
     $        1d0/18d0*(6d0*nflav*Tf-30d0*CAc)*(lnu43+lnt43)) !multiply born_43

         do i = 2,3 
            uucc(i) = uuccv(i) + (fac(i)+fac(1)) * uuccb(i)
            uuss(i) = uussv(i) + (fac(i)+fac(1)) * uussb(i)
            ddcc(i) = ddccv(i) + (fac(i)+fac(1)) * ddccb(i)
            ddss(i) = ddssv(i) + (fac(i)+fac(1)) * ddssb(i)
            udsc(i) = udscv(i) + (fac(i)+fac(1)) * udscb(i)
            ducs(i) = ducsv(i) + (fac(i)+fac(1)) * ducsb(i)
         enddo

      endif !nlo.eq.1
#endif /* WITH_NLO */

    
c     create NCmatrix and CCmatrix born + virt 
      do i=2,3
         NCmatrix(0,0,i)=uucc(i) 
         NCmatrix(0,1,i)=uuss(i)
         NCmatrix(1,0,i)=ddcc(i)
         NCmatrix(1,1,i)=ddss(i)
         CCmatrix(0,i)=udsc(i)
         CCmatrix(1,i)=ducs(i)
      enddo

c     create NCmatrix and CCmatrix born 
c     used for the finite collinear terms
c     j = 1 born, j = 2 born_21, j=3 born_43

#ifdef WITH_NLO
      if (nlo.eq.1) then

         do j = 2,3
            NCbornmsq(0,0,j)=uuccb(j)
            NCbornmsq(0,1,j)=uussb(j)
            NCbornmsq(1,0,j)=ddccb(j)
            NCbornmsq(1,1,j)=ddssb(j)
            CCbornmsq(0,j)=udscb(j)
            CCbornmsq(1,j)=ducsb(j)
         enddo

      endif !nlo.eq.1
#endif /* WITH_NLO */

c first the case with one  final state gluon
      if((gsign.eq.1)) then
         if (ldebug) then
            print*," final state gluon section in wbf_H3j "            
            print 101," fsign = ",fsign
 101        format(a,6i5,a,i5)
         endif
c          
         iflav(5) = 0           ! final state gluon id
         do if1=1,nflVBF      !(nfl/2)*2
            do if2=1,nflVBF   !(nfl/2)*2
               iflav(1)=if1*fsign(physToDiag(1))!
               iflav(3)=if1*fsign(physToDiag(3))!
               iflav(2)=if2*fsign(physToDiag(2))
               iflav(4)=if2*fsign(physToDiag(4))

               do jc = 2,3      !loop over upper and lower
                  k=FL_H3j(iflav,jc)
                  res(k)=pdf(sign1*iflav(1),1,1) !1
     &                 *pdf(sign2*iflav(2),2,1) !2 
     &                 *NCmatrix(mod(if1,2),mod(if2,2),jc)*polcolqq

#ifdef WITH_NLO
                  if (nlo.eq.1) then       

c     include finite collinear terms 
c     get born
                     do j=2,3
                        born(j)= NCbornmsq(mod(if1,2),mod(if2,2),j)
                     enddo
c     
c     compute color correlations: colcorr(5) 
c     1 is the msq, 2 is msq^{1a},3 is msq^{3a},
c     4 is msq^{2b}, 5 is msq^{3b}
                     if(jc.eq.2) then !upper line
                        colcorr(2) = (CF-CAc/2.0d0)/CF*born(2)
                        colcorr(3) = CAc/2.0d0/CF*born(2)
                        colcorr(4) = born(2)

                        res(k) = res(k) + pdf(sign1*iflav(1),1,2) !1(1a)
     &                       *pdf(sign2*iflav(2),2,1)* !2
     $                       colcorr(2)*polcolqq 
                        res(k)= res(k)+ pdf(sign1*iflav(1),1,3) !1(3a)
     &                       *pdf(sign2*iflav(2),2,1)* !2
     $                       colcorr(3)*polcolqq
                        res(k) = res(k) +  pdf(sign1*iflav(1),1,1) !1(2b)
     &                       *pdf(sign2*iflav(2),2,4)* !2
     $                       colcorr(4)*polcolqq    
                       
                        
                     elseif(jc.eq.3) then !lower line
                        colcorr(2) = born(3)
                        colcorr(4) = (CF-CAc/2.0d0)/CF*born(3)
                        colcorr(5) = CAc/2.0d0/CF*born(3)

                        res(k) = res(k) + pdf(sign1*iflav(1),1,2) !1(1a)
     &                       *pdf(sign2*iflav(2),2,1)* !2
     $                       colcorr(2)*polcolqq 
                        
                        res(k) = res(k) +  pdf(sign1*iflav(1),1,1) !1(2b)
     &                       *pdf(sign2*iflav(2),2,4)* !2
     $                       colcorr(4)*polcolqq    
                        res(k) = res(k) +  pdf(sign1*iflav(1),1,1) !1(3b)
     &                       *pdf(sign2*iflav(2),2,5)* !2
     $                       colcorr(5)*polcolqq                          
                     endif

                  endif !nlo.eq.1
#endif /* WITH_NLO */

               enddo            !loop over upper and lower line
     
C     Now check if there is a CC contribution for this choice of initial state
C     flavors (i,j). First: check whether initial uu or dd, not ud
               sametype=(mod(if1,2)).eq.(mod(if2,2))
               oneAntiparticle=
     &              ((fsign(physToDiag(1))*fsign(physToDiag(2))).eq.-1) 
c     true if only one particle is an antiparticle
               ChargedCurrent=(oneAntiparticle.and.sametype) .or.
     &              ( (.not.oneAntiparticle).and.(.not.sametype) )
               if( ChargedCurrent .and. if1.le.4 .and. if2.le.4 ) then   
c     change 1<-->2 and 3<-->4 in outgoing quark flavors
                  iflav(3)=(if1+2*mod(if1,2)-1)*fsign(physToDiag(3))
                  iflav(4)=(if2+2*mod(if2,2)-1)*fsign(physToDiag(4))

                  icc1 = abs(iflav(diagtophys(1)))
     
                  do jc=2,3  !loop over upper and lower 
                     k=FL_H3j(iflav,jc)
                     res(k)=pdf(sign1*iflav(1),1,1)
     &                    *pdf(sign2*iflav(2),2,1)
     &                    *CCmatrix(mod(icc1,2),jc)*polcolqq


#ifdef WITH_NLO
                     if (nlo.eq.1) then

c     include finite collinear terms 
c     get born
                        do j=2,3
                           born(j)= CCbornmsq(mod(icc1,2),j)
                        enddo

c     compute color correlations: colcorr(5) 
c     1 is the msq, 2 is msq^{1a},3 is msq^{3a},
c     4 is msq^{2b}, 5 is msq^{3b}
                        if(jc.eq.2) then !upper line
                           colcorr(2) = (CF-CAc/2.0d0)/CF*born(2)
                           colcorr(3) = CAc/2.0d0/CF*born(2)
                           colcorr(4) = born(2)
                              
                           res(k) = res(k) + pdf(sign1*iflav(1),1,2) !1(1a)
     &                             *pdf(sign2*iflav(2),2,1)* !2
     $                             colcorr(2)*polcolqq 
                           res(k)= res(k)+ pdf(sign1*iflav(1),1,3) !1(3a)
     &                             *pdf(sign2*iflav(2),2,1)* !2
     $                             colcorr(3)*polcolqq
                           res(k) = res(k) +  pdf(sign1*iflav(1),1,1) !1(2b)
     &                             *pdf(sign2*iflav(2),2,4)* !2
     $                             colcorr(4)*polcolqq    
                              
                              
                        elseif(jc.eq.3) then !lower line
                           colcorr(2) = born(3)
                           colcorr(4) = (CF-CAc/2.0d0)/CF*born(3)
                           colcorr(5) = CAc/2.0d0/CF*born(3)
                              
                           res(k) = res(k) + pdf(sign1*iflav(1),1,2) !1(1a)
     &                             *pdf(sign2*iflav(2),2,1)* !2
     $                             colcorr(2)*polcolqq 
                              
                           res(k) = res(k) +  pdf(sign1*iflav(1),1,1) !1(2b)
     &                             *pdf(sign2*iflav(2),2,4)* !2
     $                             colcorr(4)*polcolqq    
                           res(k) = res(k) +  pdf(sign1*iflav(1),1,1) !1(3b)
     &                             *pdf(sign2*iflav(2),2,5)* !2
     $                             colcorr(5)*polcolqq                    
                        endif
                           
                     endif !nlo.eq.1
#endif /* WITH_NLO */

                  enddo      !loop over upper and lower line
               endif         !endif(CC)     
            enddo            !if2
         enddo               !if1
     
      elseif((gsign.eq.-1)) then ! one final state gluon        

         l43 = 1.0d0
         l21 = 1.0d0
     
         jc = 0                 !color structure
         if (fsign(1).eq.-fsign(2) .and. 
     1        fsign(3).eq.fsign(4) ) then
            l43 = 0             !  initial gluon attached to 1-2 line only
            jc = 2
         elseif (fsign(1).eq.fsign(2) .and. 
     1           fsign(3).eq.-fsign(4) ) then
            l21 = 0             !  initial gluon attached to 3-4 line only
            jc = 3
         endif

c     neutral current section
         iflav(diagtophys(5)) = 0 !gluon
         do if1=1,nflVBF     !(nfl/2)*2           !21-line
            do if2=1,nflVBF  !(nfl/2)*2        !43-line
               iflav(diagToPhys(1))=if1*fsign(1)
               iflav(diagToPhys(2))=if1*fsign(2)
               iflav(diagToPhys(3))=if2*fsign(3)
               iflav(diagToPhys(4))=if2*fsign(4) 
               k=FL_H3j(iflav,jc)
c
               res(k) = pdf(sign1*iflav(1),1,1) !1 
     &              *pdf(sign2*iflav(2),2,1) !2  
     &              *NCmatrix(mod(if1,2),mod(if2,2),jc)*polcolgq          
c     
#ifdef WITH_NLO
               if (nlo.eq.1) then
                  do j=2,3
                     born(j)= NCbornmsq(mod(if1,2),mod(if2,2),j)
                  enddo
c     
                  if(jc.eq.2) then
                     
c     born(jc)= NCbornmsq(mod(if1,2),mod(if2,2),jc)
                     
                     res(k)=res(k)+l21*pdf(sign1*iflav(1),1,2) !1,(1a) 
     &                    *pdf(sign2*iflav(2),2,1)* !2
     $                    born(2)*polcolgq/2.0d0
                     res(k)=res(k)+l21*pdf(sign1*iflav(1),1,3) !1,(3a) 
     &                    *pdf(sign2*iflav(2),2,1)* !2
     $                    born(2)*polcolgq/2.0d0
                     res(k)=res(k)+l21*pdf(sign1*iflav(1),1,1) !1,(2b) 
     &                    *pdf(sign2*iflav(2),2,4)* !2
     $                    born(2)*polcolgq
c     
                  elseif(jc.eq.3) then
                     
c                     born(jc)= NCbornmsq(mod(if1,2),mod(if2,2),jc) 
                     res(k)=res(k)+l43*pdf(sign1*iflav(1),1,1) !1,(2b) 
     &                    *pdf(sign2*iflav(2),2,4)* !2
     $                    born(3)*polcolgq/2.0d0
                     res(k)=res(k)+l43*pdf(sign1*iflav(1),1,1) !1,(3b) 
     &                    *pdf(sign2*iflav(2),2,5)* !2
     $                    born(3)*polcolgq/2.0d0
                     res(k)=res(k)+l43*pdf(sign1*iflav(1),1,2) !1,(1a) 
     &                    *pdf(sign2*iflav(2),2,1)* !2
     $                    born(3)*polcolgq
                  endif         !jc

               endif !nlo.eq.1
#endif /* WITH_NLO */

            enddo               !if2
         enddo                  !if1
c     
C     Next the CC contributions; no b quarks allowed here
         iflav(diagToPhys(5)) = 0 !gluon
         do if1 = 1,(nflVBF/2)*2    !4           ! sum over all flavors for quark 1 and 2 CHARGED CURRRENT!
            iflav(diagToPhys(1)) = fsign(1)*if1
            iflav(diagToPhys(2)) = fsign(2)*(if1-1+2*mod(if1,2))
            do if2 = mod(if1,2)+1,(nflVBF/2)*2,2!4,2 ! flavor of q3 set by q2 (mod 2)
               iflav(diagToPhys(3)) = fsign(3)*if2
               iflav(diagToPhys(4)) = fsign(4)*(if2-1+2*mod(if2,2))
               k=FL_H3j(iflav,jc)
               res(k)=pdf(sign1*iflav(1),1,1)
     &                    *pdf(sign2*iflav(2),2,1)
     &                    *CCmatrix(mod(if1,2),jc)*polcolgq
#ifdef WITH_NLO
               if (nlo.eq.1) then

                  do j=2,3
                     born(j)= CCbornmsq(mod(if1,2),j)
                  enddo
c     
                  if(jc.eq.2) then
c     gQ3->qb1q2 Q4
c     born(jc)= CCbornmsq(mod(if1,2),jc)
                     res(k)=res(k)+l21*pdf(sign1*iflav(1),1,2) !1,(1a) 
     &                    *pdf(sign2*iflav(2),2,1)* !2
     $                    born(2)*polcolgq/2.0d0
                     res(k)=res(k)+l21*pdf(sign1*iflav(1),1,3) !1,(3a) 
     &                    *pdf(sign2*iflav(2),2,1)* !2
     $                    born(2)*polcolgq/2.0d0
                     res(k)=res(k)+l21*pdf(sign1*iflav(1),1,1) !1,(2b) 
     &                    *pdf(sign2*iflav(2),2,4)* !2
     $                    born(2)*polcolgq
                  elseif(jc.eq.3) then
c     q1g->Qb3Q4 q2
c     born(jc)= CCbornmsq(mod(if1,2),jc)
                     res(k)=res(k)+l43*pdf(sign1*iflav(1),1,1) !1,(2b) 
     &                    *pdf(sign2*iflav(2),2,4)* !2
     $                    born(3)*polcolgq/2.0d0
                     res(k)=res(k)+l43*pdf(sign1*iflav(1),1,1) !1,(3b) 
     &                    *pdf(sign2*iflav(2),2,5)* !2
     $                    born(3)*polcolgq/2.0d0
                     res(k)=res(k)+l43*pdf(sign1*iflav(1),1,2) !1,(1a) 
     &                    *pdf(sign2*iflav(2),2,1)* !2
     $                    born(3)*polcolgq
                  endif         !jc

               endif !nlo.eq.1
#endif /* WITH_NLO */

            enddo
         enddo
     
      else
         stop
      endif

      nmax = FL_H3j(iflav,-2)
      return
      end
      
c******************************************************************
c
c   end subroutine wbf_h3j
c
c*****************************************************************


#ifdef WITH_NLO
c     finite collinear coefficients
      subroutine fincollinear(A,B,C,D,xi,RN,p,ip,ap)
c     This subroutine computes the coefficients
c     A(z),B(z),C(z), and D(x)
c     see notes
c
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"
#include "color.inc"
C input:
      integer ndip              !dipole mom. conf. label
      integer ip,ap             !integer label of for pi.pa
c     ip and ap are index labels 
c     by using a map one can use the same notation at CS
      real*8 p(0:3,max_p,max_kin)      !4-momentum
      real*8 xi(2)              !Feynman-x
      real*8 RN                 !random number for z integration

C output:
      real*8 A(2),B(3),C(3),D(3)
c functions
      real*8 dot0p
      external dot0p
c     local variables
      real*8 qsq,lnQomu !pi.pa,ln(pi.pa/muf^2)
      real*8 z,zmin,x,ln1mx,lnrat,omz,lnz
      real*8 gammag,gammaq,omx
      integer i 
      parameter(zmin=0.999995d0) !cutoff on z for num.stability
      
! note about nflav: number of flavours in dipole should always be 5 as that's
! what we use for alfas running.  The number of external PDFs isn't what we're
! concerned with here
      double precision nflav
      parameter (nflav = 5d0)

c     splitting coefficients are labeled as:
C     For A(n)
c     n = 1 is qg
c     n = 2 is gq
c     For B,C,D(n)
c     n = 1 is qq with gamma_q,lp is a quark
c     n = 2 is qq with gamma_g,lp is a gluon
c     n = 3 is gg
c
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     gamma_g and gamma_q
! note about nflav: number of flavours in dipole should always be 5 as that's
! what we use for alfas running.  THe number of external PDFs isn't what we're
! concerned with here
!      gammag = 11d0/6d0*CAc-2d0/3d0*NFL*Tf
      gammag = 11d0/6d0*CAc-2d0/3d0*nflav*Tf
      gammaq = 1.5d0*CF
c      print*,'gammag=',gammag
c      print*,'gammaq=',gammaq
c      print*,'creal=',creal
c      
      x = xi(ap)
c
c      print*,'rn=',RN
      z = (1.0d0 - x)*RN**2 + x
c      z = RN
c
      omx = 1.0d0 - x
      omz = 1.0d0 - z
c
      qsq = 2.0d0*dot0p(p(0,ip,1),p(0,ap,1))
c      print*,'qsq=',qsq
      lnQomu = log(qsq/mufsq(1,1))
c      print*,'lnQomu=',lnQomu
c
      ln1mx = log(omx)
c0.999995d0
      if(z.lt.zmin) then
         lnrat = lnQomu + log((1.0d0 - z)/z)
c         print*,'lnrat=',lnrat
         lnz = log(z)
c
         A(1) = CF*((1.0d0+(1.0d0 - z)**2)/z*lnrat + z) !A_qg^{l,a}
         A(2) = Tf*((z**2+(1-z)**2)*lnrat + 2.0d0*z*(1.0d0-z)) !A_gq^{l,a}
c         A(1)=0.0d0
c         A(2)=0.0d0
c
         B(1) = CF*(2.0d0/omz * (lnrat + lnz)-gammaq/CF/omz) !B_qq^{l,a}
         B(2) = CF*(2.0d0/omz * (lnrat + lnz)-gammag/CAc/omz)
         B(3) = CAc*(2.0d0/omz * (lnrat + lnz)-gammaq/CF/omz) !B_gg^{l,a}
c
c         B(3)=0.0d0
         C(1) = CF*(-(1.0d0+z)*lnrat - 2.0d0*lnz/omz + omz) !C_qq 
         C(2) = C(1)            !C_qq
c factor of 2 missing 
c corrected on Feb 6, 2007 
c see notes
         C(3) = 2.0D0*CAc*(-lnz/omz + (omz/z - 1.0d0 + z*omz)*lnrat) !C_gg
c         C(3) = CAc*(-2.0d0*lnz/omz + (omz/z - 1.0d0 + z*omz)*lnrat) !C_gg
c         C(3)=0.0d0
c     
      else
         A(1) = 0.0d0
         A(2) = 0.0d0
         do i=1,3
            B(i)=0.0d0
            C(i)=0.0d0
         enddo
      endif
c
! note about nflav: number of flavours in dipole should always be 5 as that's
! what we use for alfas running.  THe number of external PDFs isn't what we're
! concerned with here
      D(1) = CF*(crealH3j + 1.5d0 -(gammaq/CF)-(gammaq/CF)*ln1mx
     $     + ln1mx**2 + 1.5d0*lnQomu + 2.0d0*ln1mx * lnQomu)!D_qq
      D(2) = CF*(crealH3j + 1.5d0 -(gammag/CAc)-(gammag/CAc)*ln1mx
     $     + ln1mx**2 + 1.5d0*lnQomu + 2.0d0*ln1mx * lnQomu)
      D(3) = CAc*(2.0d0*ln1mx * lnQomu + (gammag/CAc)*lnQomu +
!     $     crealH3j - 5d0/9d0 + ln1mx**2 + 16d0/9d0*NFL*Tf/CAc -
     $     crealH3j - 5d0/9d0 + ln1mx**2 + 16d0/9d0*nflav*Tf/CAc -
     $     3d0/2d0*ln1mx)       !D_gg
c      do i=1,3
c         D(i)=0.0d0
c      enddo
c     
      end
#endif
ccccccccccccccccccccccccc

c******************************************************************
c
c   begin function FL_H3j
c
c*****************************************************************
      integer function FL_H3j(iflav,colstruc)
      implicit none
      integer iflav(5),colstruc    ! input for color and flavor assignment
c  fill hepup color and flavor assignments for WBF processes with
c  one attached external gluon and count subprocesses. 
c  There are 2 color structures depending into which fermion line the 
c  gluon is inserted
c
c  colstruc = 2:upper
c             3:lower fermion line
c
c  cases can be distinguished according to which iflav(i)=0, i.e. corresponds
c  to the gluon

c  id1,2 are flavor identifiers for incoming quarks
c  id3,4 are flavor identifiers for outgoing quarks
c  gsign=+1: outgoing gluon, gsign=-1 for incoming gluon
c  colstruc =2,3 determines whether the gluon is coupled to the 
c  uppper or lower fermion line:

c  Note that colstruc is also a flag:
c  if colstruc=-1, then we are resetting.

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
      
      integer lkup
      common /localkup/ lkup(numParticles,7,maxNumSubProc)
 
      integer listposition
      integer numresets
      save listposition
      save numresets

      data numresets /-1/
      data listposition /0/

      if (colstruc.eq.-1) then  ! we are resetting
c       write(*,*) "we are resetting"
        listposition=0
        numresets=numresets+1
      else if (colstruc.eq.-2) then
c       don"t do anything, just return the number of subprocesses.
      else
c increment the counter regardless of whether or not it"s the 1st time through.
        listposition=listposition+1

c       we fill info for this subprocess,if required
        if(numresets.eq.0) then
          call fillColoredPartons_H3j(iflav,
     1                                listposition+numdecay,colstruc)
        endif       ! numresets.eq.0
      endif     
c      print*,'listposition=',listposition
      FL_H3j=listposition

      end
c******************************************************************************
c
c   end function FL_H3j
c
c******************************************************************************


c*****************************************************************************
c
c    begin  subroutine fillColoredPartons_H3j
c
c*****************************************************************************
      subroutine fillColoredPartons_H3j(iflav,listposition,colstruc)
c  assigns values to the variables in the common block localHEPUP
c  in particular, this subroutine assigns values to those variables that 
c  will be stored in the lookup tables generated by writeHEPUPtable.  
c  As the name suggests, this routine only stores the information for the 
c  colored partons.  Particles without color will be dealt with in the 
c  subroutine fillColorless.

      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/process.inc"

      integer iflav(5),id1,id2,id3,id4,id5,listposition,i,i1
      integer iflavour(5),colstruc

      logical ldebug
      parameter (ldebug=.false.)


      do i1=1,5
         iflavour(i1)=iflav(i1) 
         if(iflavour(i1).eq.0) iflavour(i1)=21   
      enddo

c...flavours of 4 quarks and gluon

      id1=iflavour(1)
      id2=iflavour(2)
      id3=iflavour(3)
      id4=iflavour(4)
      id5=iflavour(5)

      select case (process)
      case(Hjjj_WW, Hjjj_ZZ_ll, Hjjj_ZZ_lnu)
         lnup(listposition)=numParticles+3 
      case(Hjjj)
         lnup(listposition)=numParticles-2
      case(Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_bbar)
         lnup(listposition)=numParticles            
      end select
      
      listup(1,listposition)=-1 !incoming quarks
      listup(2,listposition)=-1
      listup(3,listposition)=1  !outgoing quarks
      listup(4,listposition)=1
      listup(5,listposition)=1
 
c...four quarks and gluon
      lidup(1,listposition)=id1    
      lidup(2,listposition)=id2
      lidup(3,listposition)=id3
      lidup(4,listposition)=id4
      lidup(5,listposition)=id5  
      
      do i=3,n_p
         lmothup(1,i,listposition)=1
         lmothup(2,i,listposition)=2
      enddo

c...final state gluon
      if(id5.eq.21)then
         if(colstruc.eq.2)then
            if(id1.gt.0) then
               licolup(1,1,listposition)=501
               licolup(1,3,listposition)=502
               licolup(1,5,listposition)=501
               licolup(2,1,listposition)=0
               licolup(2,3,listposition)=0
               licolup(2,5,listposition)=502
            else
               licolup(1,1,listposition)=0
               licolup(1,3,listposition)=0
               licolup(1,5,listposition)=502
               licolup(2,1,listposition)=501
               licolup(2,3,listposition)=502
               licolup(2,5,listposition)=501               
            endif            
            if(id2.gt.0) then
               licolup(1,2,listposition)=503
               licolup(1,4,listposition)=503
               licolup(2,2,listposition)=0
               licolup(2,4,listposition)=0
            else
               licolup(1,2,listposition)=0
               licolup(1,4,listposition)=0
               licolup(2,2,listposition)=503
               licolup(2,4,listposition)=503
            endif
         elseif(colstruc.eq.3)then
            if(id1.gt.0) then
               licolup(1,1,listposition)=501
               licolup(1,3,listposition)=501
               licolup(2,1,listposition)=0
               licolup(2,3,listposition)=0
            else
               licolup(1,1,listposition)=0
               licolup(1,3,listposition)=0
               licolup(2,1,listposition)=501
               licolup(2,3,listposition)=501
            endif         
            if(id2.gt.0) then
               licolup(1,2,listposition)=502
               licolup(1,4,listposition)=503
               licolup(1,5,listposition)=502
               licolup(2,2,listposition)=0
               licolup(2,4,listposition)=0
               licolup(2,5,listposition)=503
            else
               licolup(1,2,listposition)=0
               licolup(1,4,listposition)=0
               licolup(1,5,listposition)=503
               licolup(2,2,listposition)=502
               licolup(2,4,listposition)=503
               licolup(2,5,listposition)=502
            endif
         endif
      endif

c...initial particles gluon and quark/antiquark 
      if(id1.eq.21)then 
         if(id3.gt.0)then
            licolup(1,1,listposition)=502
            licolup(1,3,listposition)=502
            licolup(2,1,listposition)=501
            licolup(2,3,listposition)=0
         else
            licolup(1,1,listposition)=501
            licolup(1,3,listposition)=0
            licolup(2,1,listposition)=502
            licolup(2,3,listposition)=502
         endif 
         if(id5.gt.0)then
            licolup(1,5,listposition)=501
            licolup(2,5,listposition)=0
         else
            licolup(1,5,listposition)=0
            licolup(2,5,listposition)=501            
         endif
         if(id2.gt.0) then
            licolup(1,2,listposition)=503
            licolup(1,4,listposition)=503
            licolup(2,2,listposition)=0
            licolup(2,4,listposition)=0
         else
            licolup(1,2,listposition)=0
            licolup(1,4,listposition)=0
            licolup(2,2,listposition)=503
            licolup(2,4,listposition)=503
         endif         
      endif
      
c...initial particles quark/antiquark and gluon 
      if(id2.eq.21)then 
         if(id1.gt.0)then
            licolup(1,1,listposition)=501
            licolup(1,3,listposition)=501
            licolup(2,1,listposition)=0
            licolup(2,3,listposition)=0
         else
            licolup(1,1,listposition)=0
            licolup(1,3,listposition)=0
            licolup(2,1,listposition)=501
            licolup(2,3,listposition)=501
         endif 
         if(id5.gt.0)then
            licolup(1,5,listposition)=502
            licolup(2,5,listposition)=0
         else
            licolup(1,5,listposition)=0
            licolup(2,5,listposition)=502            
         endif
         if(id4.gt.0) then
            licolup(1,2,listposition)=503
            licolup(1,4,listposition)=503
            licolup(2,2,listposition)=502
            licolup(2,4,listposition)=0
         else
            licolup(1,2,listposition)=502
            licolup(1,4,listposition)=0
            licolup(2,2,listposition)=503
            licolup(2,4,listposition)=503
         endif
      endif

      if(ldebug)then
         print*,"====================================="
         print("(1P,5I6)"),id1,id2,id3,id4,id5
         print*,"-------------------------------------"
         print("(1P,5I6)"),(licolup(1,i1,listposition),i1=1,5)
         print("(1P,5I6)"),(licolup(2,i1,listposition),i1=1,5)
         print*,"====================================="
         print*,""
      endif
    
      end  ! fillColoredPartons_H3j
c*****************************************************************************
c
c    end subroutine fillColoredPartons_H3j
c
c*****************************************************************************
