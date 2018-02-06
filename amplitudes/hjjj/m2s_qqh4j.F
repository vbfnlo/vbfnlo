c
c	Terrance Figy, <terrance@physics.wisc.edu>
c	Initial version:  2005 July
c	Last modified:    2017 February, Michael Rauch
c   
c     
      subroutine m2s_wbfh4j(
     1                   bos,   !in:  Boson identifier, 2 = Z, 6=H
     3                   nlo,   !in:  NLO=1: create subtraction term; LO = 0
     4                   lok,   !in:  need to calculate m2s(1:max_kin) T or F 
     5                   xi,	!in:  Feynman x parameters of incoming fermions
     6                   p,	!in:  momenta associated with external fermions
     7                   v,	!in:  momenta associated with WW decay fermions
     8                   xuz,   !in:  x, u, z for subtraction terms
     9                   m2s  )  !out: |M|^2*pdf1*pdf2 
          use globalvars, only: ldoblha
          use VBFNLO_HIGGSTO, only: VBFNLO_HIGGSTO_gridHt
      implicit none
c
c declare input/output variables
c
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
      integer bos
      double precision xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), 
     1     xuz(2,2:max_kin), m2s(0:max_kin)
      logical lok(max_kin)

c
c declare external functions
c
      integer FL_H4j
      external FL_H4j

c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"

c declare local variables
      double precision q_sf
      integer init/0/, I, J,k, mu
      save init

c declare variables for summation process
      INTEGER ires

      integer FSIGN(4+max_v), gsign(2), physToDiag(6), 
     $     nlo, nd,fs(6+max_v)
      integer nmin,nmax,nproc(13)
      save nmin, nmax, nproc

c store contributions from subprocess ip in res(ip,ID) where
c  ID = 1  : the real emission |M|^2 * pdf
c     = 2,3: sutraction terms for emision off 12 or 34 line
c     = 0  : subtracted result which drives integration, i.e
c res(*,0) = res(*,1)+res(*,2)+res(*,3)+...+res(*,max_kin)
      double precision res(maxnumsubproc,0:max_kin)


      double precision sum,delta
      double precision dot0p
      external dot0p

c declare parton distribution variables
      double precision pdf(-6:6,2,5) 

c define program switches
      logical lskip,ldebug1,linfo
      logical lchsb,lcoll,ldebugc
      parameter (lchsb=.false.,lcoll=.false.)
      parameter (lskip=.false.)
      parameter(ldebug1=.false.,linfo=.false.)
      parameter(ldebugc=.false.)

      double precision weight,rnumb,RandomNumber

c CPS scheme
      DOUBLE PRECISION CLR,XM2,XMG,Bcoupl,Vcoupl,Acoupl
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),Bcoupl(6,6,6),Vcoupl(4,5),Acoupl(4,5)
      double precision qinv, qvv(0:3)
      double precision qgammaq
      common /VBFNLO_HIGGS_CPS/ qgammaq

c if first iteration,
      if ( init .eq. 0 .and. .not.ldoblha ) then
         if (bos.eq.6) then
            write(6,*) " "
            write(6,*) "H4j amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         endif
         call printnfl(.true.)
         print*," "

         init = 1
         do i = 1,maxnumsubproc
            do j = 0,max_kin
               res(i,j) = 0.0d0
            enddo
         enddo
      else if ( init .eq. 0 .and. ldoblha ) then
         init = 1
      else
         do i = 1,nmax
            do j=0,max_kin
               res(i,j) = 0
            enddo
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

       if (ldoblha) then

         physToDiag(1)=blha_physdiagmap(1,blha_cursubproc,blha_curproc)
         physToDiag(2)=blha_physdiagmap(2,blha_cursubproc,blha_curproc)
         physToDiag(3)=blha_physdiagmap(3,blha_cursubproc,blha_curproc)
         physToDiag(4)=blha_physdiagmap(4,blha_cursubproc,blha_curproc)
         physToDiag(5)=blha_physdiagmap(5,blha_cursubproc,blha_curproc)
         physToDiag(6)=blha_physdiagmap(6,blha_cursubproc,blha_curproc)
         fsign(1) = blha_fsign(1,blha_cursubproc,blha_curproc)
         fsign(2) = blha_fsign(2,blha_cursubproc,blha_curproc)
         fsign(3) = blha_fsign(3,blha_cursubproc,blha_curproc)
         fsign(4) = blha_fsign(4,blha_cursubproc,blha_curproc)
         fsign(5) = blha_fsign(5,blha_cursubproc,blha_curproc)
         fsign(6) = blha_fsign(6,blha_cursubproc,blha_curproc)
         gsign(1) = blha_fsign(5,blha_cursubproc,blha_curproc)
         gsign(2) = blha_fsign(6,blha_cursubproc,blha_curproc)

         if (mod(blha_idsubproc(blha_cursubproc,blha_curproc)/100000,10) .eq. 1) then
           ! 4q2g
           call wbf_h4jg(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     &       pdf,res,nmin,nmax)
         else if (mod(blha_idsubproc(blha_cursubproc,blha_curproc)/100000,10) .eq. 2) then
           ! 6q
           call wbf_h4jq(xi,p,v,physToDiag,fsign,bos,nlo,lok,xuz,
     &       pdf,res,nmin,nmax)
         endif
         return
       endif

c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
      pdf = 0

      call ctrans
      
      q_sf = sqrt(mufsq(1,1))
      call pdfproton( xi(1), q_sf, pdf(-6,1,1)) !f_a/p(x1/x_{41,a})

      q_sf = sqrt(mufsq(2,1))
      call pdfproton( xi(2), q_sf, pdf(-6,2,1)) !f_b/p(x2/x_{42,b})         
     
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Now sum over the subprocesses contributing to H production
      nmax = FL_H4j(fsign,fsign,-1)   !reset counter for subprocesses to 0

C*******************  q1 q3 ---> q2 q4 g g H   **********************
c         
c   physToDiag(ext.momentum label) = Feynman diagram label
c
      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon
      physToDiag(6)=6   ! gluon
C NOTE: for call of wbf_zh4j it is important that p(*,1,*) and p(*,3,*)
c correspond to 1-2 fermion line ALWAYS, i.e physToDiag(1/2)={1,3} and 
c similarly physToDiag(3/4)={2,4} for the 3-4 fermion line
      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1
      gsign(1) = 1
      gsign(2) = 1

      call wbf_h4jg(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1     pdf,res,nmin,nmax)   

      if(linfo) print*,"nmin,nmax1",nmin,nmax     
      if (init.eq.1) nproc(1) = nmax 

C*******************  q1 qb4 ---> q2 qb3 g g H   **********************

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon
      physToDiag(6)=6   ! gluon
      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = -1
      fsign(4) = -1
      gsign(1) = 1
      gsign(2) = 1

      call wbf_h4jg(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

      if(linfo) print*,"nmin,nmax",nmin,nmax
      if (init.eq.1) nproc(2) = nmax 

C*******************  qbar2 q3 ---> qbar1 q4 g g H   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon
      physToDiag(6)=6   ! gluon

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = 1
      gsign(1) = 1
      gsign(2) = 1
c      
      call wbf_h4jg(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

      if(linfo) print*,"nmin,nmax3",nmin,nmax
      if (init.eq.1) nproc(3) = nmax 

C*******************  qbar2 qb4 ---> qbar1 qb3 g g H   ******************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon
      physToDiag(6)=6   ! gluon

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1 
      fsign(4) = -1
      gsign(1) = 1
      gsign(2) = 1
         
      call wbf_h4jg(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

      if(linfo) print*,"nmin,nmax",nmin,nmax
      if (init.eq.1) nproc(4) = nmax 

C******************* g5 q3 ---> qb1 q2 q4 g6 H   **********************
c                    p1 p2      p5  p3 p4 p6           
c     
      physToDiag(1)=5
      physToDiag(2)=3
      physToDiag(3)=2
      physToDiag(4)=4
      physToDiag(5)=1
      physToDiag(6)=6

      fsign(1) = -1
      fsign(2) =  1
      fsign(3) =  1
      fsign(4) =  1
      gsign(1) = -1
      gsign(2) = 1

      call wbf_h4jg(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

      if(linfo) print*,"nmin,nmax5",nmin,nmax
      if (init.eq.1) nproc(5) = nmax 

C*******************  g qbar4 ---> qbar3 qb1 q2 g H   **********************

      physToDiag(1)=5
      physToDiag(2)=4
      physToDiag(3)=2
      physToDiag(4)=3
      physToDiag(5)=1
      physToDiag(6)=6

      fsign(1) = -1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) = -1
      gsign(1) = -1
      gsign(2) = 1

      call wbf_h4jg(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

      if(linfo) print*,"nmin,nmax",nmin,nmax
      if (init.eq.1) nproc(6) = nmax 

C*******************  q1 g ---> q2 qb3 q4 g H   **********************
c                     p1 p2     p3 p5  p4 p6
      physToDiag(1)=1
      physToDiag(2)=5
      physToDiag(3)=2
      physToDiag(4)=4
      physToDiag(5)=3
      physToDiag(6)=6

      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      gsign(1) = -1
      gsign(2) =  1

      call wbf_h4jg(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

      if(linfo) print*,"nmin,nmax7",nmin,nmax
      if (init.eq.1) nproc(7) = nmax 

C*******************  qb2 g ---> qb1 qb3 q4 g H   *************

      physToDiag(1)=2
      physToDiag(2)=5
      physToDiag(3)=1
      physToDiag(4)=4
      physToDiag(5)=3
      physToDiag(6)=6

      fsign(1) =  -1
      fsign(2) =  -1
      fsign(3) = -1
      fsign(4) =  1
      gsign(1) = -1
      gsign(2) =  1            

      call wbf_h4jg(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

      if(linfo) print*,"nmin,nmax",nmin,nmax
      if (init.eq.1) nproc(8) = nmax 

c******************g g -> qb1 q2 qb3 q4 
c                 p1 p2   p5  p3 p6  p4 
 
      physToDiag(1)=5
      physToDiag(2)=6
      physToDiag(3)=2!2
      physToDiag(4)=4!4
      physToDiag(5)=1!1
      physToDiag(6)=3!3

      fsign(1) = -1
      fsign(2) = 1
      fsign(3) = -1
      fsign(4) = 1
      gsign(1) = -1
      gsign(2) = -1            

      call wbf_h4jg(xi,p,v,physToDiag,fsign,gsign,bos,nlo,lok,xuz,
     1              pdf,res,nmin,nmax)

      if(linfo) print*,"nmin,nmax9",nmin,nmax
      if (init.eq.1) nproc(9) = nmax 
     
c***************** q1 q3 -> q2 q4 qb5 q6
c                  p1 p2    p3 p4  p5 p6
      phystodiag(1)=1           !1/2
      phystodiag(2)=3           !3/4
      phystodiag(3)=2           !1/2
      phystodiag(4)=4           !3/4
      phystodiag(5)=5           !anti-quark
      phystodiag(6)=6           !quark
      
      fs(1) = 1
      fs(2) = 1
      fs(3) = 1
      fs(4) = 1
      fs(5) = -1
      fs(6) = 1

      call wbf_h4jq(xi,p,v,physToDiag,fs,bos,nlo,lok,xuz,
     1     pdf(-6,1,1),res,nmin,nmax)

      if(linfo) print*,"nmin,nmax 10",nmin,nmax
      if (init.eq.1) nproc(10) = nmax 
  
c****************qb2 q3 -> qb1 q4 qb5 q6
c                p1 p2     p3  p4 p5  p6
      phystodiag(1)=2           !1/2
      phystodiag(2)=3           !3/4
      phystodiag(3)=1           !1/2
      phystodiag(4)=4           !3/4
      phystodiag(5)=5           !anti-quark
      phystodiag(6)=6           !quark
      
      fs(1) = -1
      fs(2) = -1
      fs(3) = 1
      fs(4) = 1
      fs(5) = -1
      fs(6) = 1
      call wbf_h4jq(xi,p,v,physToDiag,fs,bos,nlo,lok,xuz,
     1     pdf(-6,1,1),res,nmin,nmax)

      if(linfo) print*,"nmin,nmax 11",nmin,nmax
      if (init.eq.1) nproc(11) = nmax 
     
c****************q1 qb4 -> q2 qb3 qb5 q6
c                p1  p2    p3  p4 p5  p6
      phystodiag(1)=1           !1/2
      phystodiag(2)=4           !3/4
      phystodiag(3)=2           !1/2
      phystodiag(4)=3           !3/4
      phystodiag(5)=5           !anti-quark
      phystodiag(6)=6           !quark
      
      fs(1) = 1
      fs(2) = 1
      fs(3) = -1
      fs(4) = -1
      fs(5) = -1
      fs(6) = 1

      call wbf_h4jq(xi,p,v,physToDiag,fs,bos,nlo,lok,xuz,
     1     pdf(-6,1,1),res,nmin,nmax)
           if(linfo) print*,"nmin,nmax 12",nmin,nmax
      if (init.eq.1) nproc(12) = nmax 
      
c**********************qb2 qb4 -> qb1 qb3 qb5 q6
c                      p1  p2     p3  p4 p5  p6

      phystodiag(1)=2           !1/2
      phystodiag(2)=4           !3/4
      phystodiag(3)=1           !1/2
      phystodiag(4)=3           !3/4
      phystodiag(5)=5           !anti-quark
      phystodiag(6)=6           !quark
      
      fs(1) = -1
      fs(2) = -1
      fs(3) = -1
      fs(4) = -1
      fs(5) = -1
      fs(6) = 1
c
      call wbf_h4jq(xi,p,v,physToDiag,fs,bos,nlo,lok,xuz,
     1     pdf(-6,1,1),res,nmin,nmax)

      if(linfo) print*,"nmin,nmax13",nmin,nmax
      if (init.eq.1) nproc(13) = nmax 
     
c*****************  end of process evaluation  ********************

      do nd=0,n_kin
         m2s(nd) = 0.0d0
         DO IRES = 1,NMAX            
            m2s(nd) = m2s(nd) + RES(IRES,nd)
         ENDDO 
c
         if(nd.eq.0) then
c...  Les Houches interface - the most probable subprocess 3jets at LO
            if ((lha.or.hepmc).and..not.doNLO) then
               i=0
               weight=0.d0
               rnumb=RandomNumber()
               do while((i.le.nmax).and.(weight.le.rnumb*m2s(0)))
                  i=i+1
                  weight=weight+res(i,0)
                  iprocess=i
               enddo
            endif  
c...  Les Houches interface       
c 
         endif
      enddo 
c
#ifdef WITH_NLO
c     check subtraction here
      if(lchsb) then
c     take the limit of pi.pj -> 0               
         delta = dot0p(p(0,4,1),p(0,5,1))
         if(abs(delta).lt.5d0.and.lok(1)) then
            print*,"----------------"
            sum = 0.0d0
            do nd = 2,n_kin
               if(lok(nd)) then 
                  print*,"D(",nd,")=",m2s(nd),lok(nd)
               endif
               sum = sum + m2s(nd)
            enddo
            print*,"m2real=",m2s(1)
            print*,"dip_sum=",sum
            print*,"m2s_real=",m2s(1),lok(1)
            print*,"m2s_sub=",m2s(0)
c            print*,"All subprocesses"
c            print*,"alpha_s=",als(1,1)
            print*,"delta=",delta
            print*,"---------------------------"
            if(abs(sum).gt.0d0) then
               print*,"m2real/dip_sum",abs(m2s(1)/sum)
            endif
            print*,"---------------------------"
            print*,m2s(1)+sum,"=?",m2s(0)
         endif    
      endif                     !debug
#endif /* WITH_NLO */

      RETURN
      END

c******************************************************************
c
c   begin subroutine wbf_h4jg
c
c*****************************************************************
      subroutine wbf_h4jg(xi,p,v,physToDiag,fsign,gsign,bos, 
     1     nlo,lok,xuz,pdf,res,nmin,nmax)
      use globalvars, only: ldoblha
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
      double precision p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), xi(nx), 
     $     xuz(2,2:max_kin)
      double precision pdf(-6:6,2,5), res(maxnumsubproc,0:max_kin)
      integer physToDiag(6), fsign(4+max_v), gsign(2), 
     $     bos, nlo, nmin, nmax
      logical lok(max_kin)
c
      double precision uucc_c(3),uuss_c(3),ddcc_c(3),ddss_c(3),
     $     udsc_c(3),ducs_c(3)
      double precision uucc(0:6),uuss(0:6),ddcc(0:6),ddss(0:6),
     $     udsc(0:6),ducs(0:6)
      double precision udscbb1     !testing purposes
c wbf_h4jg calls the amplitude square routines 
c             qqhqqgg_mg    for qq-->qqggH     for bos=6
c for the subtraction terms for the NLO cross section calculation
c
c  INPUT:  p(0:3,np,max_kin)      external physical parton momenta
c          v(0:3,nv,max_kin)      for Higgs production only
c                          the sum q(mu) = v(mu,1)+...+v(mu,nv) is needed
c          physToDiag(6)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign,gsign     sign factors for amplitude calls; see qqZqq(j)
c          nlo             nlo = 0: calculate LO only (i.e. no subtraction)
c                          nlo = 1: calculate full NLO subtraction
c          lok(max_kin)          lok(ID)=T means momenta set ID passes acceptance
c                          cuts, i.e res(k,ID) needs to be calculated
c
c     need to put in nlo switch
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
      integer ip(2:max_kin),kp(2:max_kin),jp(2:max_kin),nd,Nmx
      double precision pbar(0:3,4+max_v),qbar(0:4,2) 
      double precision polcolgg,polcolgq,polcolqq,polcolab 
      double precision NCmatrix_r(0:1,0:1,0:6)
      double precision CCmatrix_r(0:1,0:6)
      double precision NCmatrix_b(0:1,0:1,2,2:max_kin)
      double precision CCmatrix_b(0:1,2,2:max_kin)
      double precision NCmatrix_v(0:1,0:1,2,2:max_kin)
      double precision CCmatrix_v(0:1,2,2:max_kin)
      double precision sub(4,2:max_kin),subv(2:max_kin),qvec(0:3,2:max_kin)
      double precision msq_b(0:2,2:max_kin),msq_v(0:2,2:max_kin),dip(2:max_kin)

      integer iflav(6), diagToPhys(6)
      integer ICOL
      double precision dot0p
      integer FL_H4j
      external FL_H4j
      external dot0p
      logical ChargedCurrent, sametype, oneAntiparticle
      logical ldebug,ldebug1
      logical ldebug2,ldebug1n
      parameter (ldebug=.false.) 
      parameter(ldebug1=.false.,ldebug2=.false.)
      parameter(ldebug1n=.false.)
      logical ldebugmg
      parameter (ldebugmg=.false.)
      logical lcoll
      parameter(lcoll=.false.)  !switch of collinear fin piece
      integer i,if1,if2,icc1,j,k,mu, fc(4+max_v), gc
      double precision delta,sum,newdelta
      integer deltapair
      integer numdip21(12),numdip43(12),jd        !map to dipoles
      integer numdip4321(4) !only first 4 used 
c      data numdip21 /2,4,5,7,8,10,15,18,23,24,9,11/
c      data numdip43 /5,8,12,14,15,17,18,19,20,21,22,25/
c      data numdip4321 /22,18,23,9/
c      save numdip21,numdip43,numdip4321
c     set up map
      numdip21(1) =  8
      numdip21(2) = 10
      numdip21(3) = 18
      numdip21(4) = 15
      numdip21(5) = 5
      numdip21(6) = 7
      numdip21(7) = 2 
      numdip21(8) = 4
c     
      numdip21(9) = 9           !qQ
      numdip21(10)= 11
      numdip21(11) = 23         !q~Q
      numdip21(12) = 24
c     
      numdip43(1) = 5 
      numdip43(2) = 8
      numdip43(3) = 12
      numdip43(4) = 14
      numdip43(5) = 15
      numdip43(6) = 17
      numdip43(7) = 18
      numdip43(8) = 20

      numdip43(9) = 19          !qQ
      numdip43(10) = 21
      numdip43(11) = 22         !qQ~
      numdip43(12) = 25

      numdip4321(1)=22 ! 3,gQ~ 
      numdip4321(2)=18 ! 3,gQ
      numdip4321(3)=23 ! 3,q~g
      numdip4321(4)=9  ! 3,qg

      polcolqq = 1d0/(NC**2 * 8d0 * xi(1)*xi(2)) 
      polcolgq = 1d0/(NC * (NC**2-1d0) * 4d0 * xi(1)*xi(2))
      polcolgg = 1d0/((NC**2-1d0)**2 * 4d0 * xi(1)*xi(2))
c     ave over polarizations and colors
c     symmetry factor for two gluons

      if(ldebug) then
         print*,"momentum"
         do i = 1,n_p
            print*,"p(mu,i,1)=",(p(mu,i,1),mu=0,3)
         enddo
      endif
      
      do i = 0,6                !loop over color flow
         uucc(i) = 0d0
         uuss(i) = 0d0
         ddcc(i) = 0d0
         ddss(i) = 0d0
         udsc(i) = 0d0
         ducs(i) = 0d0
      enddo
c      
      nmin = nmax+1
c
      do i = 1,6
         diagToPhys(physToDiag(i)) = i
      enddo
         
c get the real emission amplitude squared, 
      do mu = 0,3
         do i = 1,6
            pbar(mu,physToDiag(i))=p(mu,i,1)
         enddo
         qbar(mu,1) = pbar(mu,5)
         qbar(mu,2) = pbar(mu,6)
      enddo
c
c
c
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
c
c              
      if (lok(1)) then
         call qqh2q2g(pbar,fsign,qbar,gsign,NCmatrix_r,CCmatrix_r)
#if WITH_MADGRAPH
         if (ldebugmg) then
           call qqhqqgg_mg_c(pbar,fsign,qbar,gsign,
     $       uucc,uuss,ddcc,ddss,udsc,ducs)
           write(*,*) "Hjjj 2q2g Comparison:"
           do i=0,6                  !loop over color flow 
             write(*,*) "  colour structure: c = ", i
             write(*,*) "  NC(0,0,c)", NCmatrix_r(0,0,i), uucc(i), abs(1-uucc(i)/NCmatrix_r(0,0,i))
             write(*,*) "  NC(0,1,c)", NCmatrix_r(0,1,i), uuss(i), abs(1-uuss(i)/NCmatrix_r(0,1,i))
             write(*,*) "  NC(1,0,c)", NCmatrix_r(1,0,i), ddcc(i), abs(1-ddcc(i)/NCmatrix_r(1,0,i))
             write(*,*) "  NC(1,1,c)", NCmatrix_r(1,1,i), ddss(i), abs(1-ddss(i)/NCmatrix_r(1,1,i))
             write(*,*) "  CC(0,c)", CCmatrix_r(0,i), udsc(i), abs(1-udsc(i)/CCmatrix_r(0,i))
             write(*,*) "  CC(1,c)", CCmatrix_r(1,i), ducs(i), abs(1-udsc(i)/CCmatrix_r(0,i))
           enddo
         endif
#endif
      endif

      if (ldoblha) then
        return
      endif

c     
c
c
c for the NLO case get the subtraction terms; 
c first the case with two  final state gluon
      if((gsign(1).eq.1).and.(gsign(2).eq.1)) then
         if (ldebug) then
            print*," final state gluon section in wbf_H4j "            
            print 101," fsign = ",fsign
 101        format(a,6i5,a,i5)
         endif

#ifdef WITH_NLO
         if (nlo.eq.1) then

c     create dipole subtraction terms
            call create_ikj(ip,kp,jp)
            do nd = 2,n_kin
               call dips(nd,p,xuz,ip(nd),kp(nd),jp(nd),sub(1,nd),
     $              subv(nd),qvec(0,nd))
            enddo
          
c     two gluons are always in the final state !
            do nd = 2,n_kin
               do mu = 0,3          
                  do i = 1,4
                     pbar(mu,physToDiag(i))=p(mu,i,nd)
                  enddo
                  qbar(mu,1)=p(mu,5,nd) !gluon momenta
               enddo

               if (lok(nd)) then
     
                  call qqhqqj_c_virt(0,pbar(0,1),fsign(1),qbar(0,1),gsign(1),nd,
     $                 uucc_c,uuss_c,ddcc_c,
     $                 ddss_c,udsc_c,ducs_c)
   
c     j = 1 is gluon emission of 21 line of born
c     j = 2 is gluon emission of 43 line of born
     
                  do j = 2,3
                     NCmatrix_b(0,0,j-1,nd)=uucc_c(j) !|M_b|^2
                     NCmatrix_b(0,1,j-1,nd)=uuss_c(j)
                     NCmatrix_b(1,0,j-1,nd)=ddcc_c(j)
                     NCmatrix_b(1,1,j-1,nd)=ddss_c(j)
                     CCmatrix_b(0,j-1,nd)=udsc_c(j)
                     CCmatrix_b(1,j-1,nd)=ducs_c(j)
                  enddo
c     for spin correlations                  
                  call qqHqqj_spcor(qvec(0,nd),pbar(0,1),fsign(1),
     $               qbar(0,1),gsign(1),nd,uucc_c,uuss_c,
     $               ddcc_c,ddss_c,udsc_c,ducs_c)
                  do j = 2,3
                     NCmatrix_v(0,0,j-1,nd)=uucc_c(j) !q_mu q_nu M_b^mu M_b*^nu
                     NCmatrix_v(0,1,j-1,nd)=uuss_c(j)
                     NCmatrix_v(1,0,j-1,nd)=ddcc_c(j)
                     NCmatrix_v(1,1,j-1,nd)=ddss_c(j)
                     CCmatrix_v(0,j-1,nd)=udsc_c(j)
                     CCmatrix_v(1,j-1,nd)=ducs_c(j)
                  enddo

               endif !nd
            enddo !dipole loop (nd)         
         endif !nlo.eq.1
#endif /* WITH_NLO */

         iflav(6) = 0           ! final state gluon id
         iflav(5) = 0           ! final state gluon id
         do if1=1,nflVBF        !(nfl/2)*2
            do if2=1,nflVBF     !(nfl/2)*2
               iflav(1)=if1*fsign(physToDiag(1))
               iflav(3)=if1*fsign(physToDiag(3))
               iflav(2)=if2*fsign(physToDiag(2))
               iflav(4)=if2*fsign(physToDiag(4))
c
               if(nlo.eq.1) then
#ifdef WITH_NLO
                  k=FL_H4j(iflav,diagToPhys,2)
                  if (lok(1)) then
                     res(k,1)=pdf(sign1*iflav(1),1,1) !1
     &                    *pdf(sign2*iflav(2),2,1) !2 
     &                    *NCmatrix_r(mod(if1,2),mod(if2,2),0)*polcolqq
                  else
                     res(k,1) = 0
                  endif
                  res(k,0) = res(k,1)
#endif
               else             !loop over color flow for LO 
                  do ICOL = 1,6
                     k=FL_H4j(iflav,diagToPhys,ICOL)
                     if (lok(1)) then
                        res(k,1)=pdf(sign1*iflav(1),1,1) !1
     &                       *pdf(sign2*iflav(2),2,1) !2 
     &                       *NCmatrix_r(mod(if1,2),mod(if2,2),ICOL)*
     &                       polcolqq
                     else
                        res(k,1) = 0
                     endif
                     res(k,0) = res(k,1)
                  enddo
               endif !nlo

#ifdef WITH_NLO
               if (nlo.eq.1) then
c     get born 
                  do nd = 2,n_kin
                     msq_b(1,nd) = NCmatrix_b(mod(if1,2),
     $                    mod(if2,2),1,nd)
                     msq_b(2,nd) = NCmatrix_b(mod(if1,2),
     $                    mod(if2,2),2,nd)
                     msq_v(1,nd) = NCmatrix_v(mod(if1,2),
     $                    mod(if2,2),1,nd)
                     msq_v(2,nd) = NCmatrix_v(mod(if1,2),
     $                    mod(if2,2),2,nd)
                  enddo
                  call dipsub(gsign,fsign,msq_b,msq_v,sub,subv,dip) !put in
c  gsign to control which dipoles are computed
               
                  do nd = 2,n_kin
                     if (lok(nd)) then
                        res(k,nd) = -pdf(sign1*iflav(1),1,1)
     &                       *pdf(sign2*iflav(2),2,1)*dip(nd)*polcolqq
c     put in collinear piece"s
                        res(k,0) = res(k,0)+res(k,nd)             
                     else
                        res(k,nd) = 0
                     endif       
                  enddo            !dipoles
               
                  if(ldebug1) then
c     take the limit of pi.pj -> 0
                    delta = dot0p(p(0,5,1),p(0,6,1))
                    deltapair=56
                    newdelta = dot0p(p(0,3,1),p(0,5,1))
                    if (newdelta.lt.delta) then
                      delta=newdelta
                      deltapair=35
                    endif
                    newdelta = dot0p(p(0,3,1),p(0,6,1))
                    if (newdelta.lt.delta) then
                      delta=newdelta
                      deltapair=36
                    endif
                    newdelta = dot0p(p(0,1,1),p(0,6,1))
                    if (newdelta.lt.delta) then
                      delta=newdelta
                      deltapair=16
                    endif
                    newdelta = dot0p(p(0,1,1),p(0,5,1))
                    if (newdelta.lt.delta) then
                      delta=newdelta
                      deltapair=15
                    endif
                    newdelta = dot0p(p(0,4,1),p(0,6,1))
                    if (newdelta.lt.delta) then
                      delta=newdelta
                      deltapair=46
                    endif
                    newdelta = dot0p(p(0,4,1),p(0,5,1))
                    if (newdelta.lt.delta) then
                      delta=newdelta
                      deltapair=45
                    endif
                    newdelta = dot0p(p(0,2,1),p(0,6,1))
                    if (newdelta.lt.delta) then
                      delta=newdelta
                      deltapair=26
                    endif
                    newdelta = dot0p(p(0,2,1),p(0,5,1))
                    if (newdelta.lt.delta) then
                      delta=newdelta
                      deltapair=25
                    endif
                    if(delta.lt.1d0) then
                      if(res(k,1).ne.0d0) then
                        print*,"NC",iflav
                        sum = 0
                        do nd = 2,n_kin
                          sum = sum + res(k,nd)
                          if(res(k,1).ne.0d0) then
                             print*,"D(",ip(nd),kp(nd),jp(nd),
     -                            ")/|m2s|",nd,
     $                            res(k,nd)/res(k,1),lok(nd)
                          else
                             print*,"D(",ip(nd),kp(nd),jp(nd),")",nd
     $                            ,res(k,nd),lok(nd)
                          endif
                        enddo
                        print*,"m2real=",res(k,1),lok(1)
                        print*,"dip_sum=",sum
                        print*,"---------------------------------"
                        print*,"dip_sum/res(k,1)",abs(sum/res(k,1))
                        print*,"----------------------------------"
                        print*,"momentum"
                        do i = 1,n_p
                        print*,"p(mu,",i,")=",(p(mu,i,1),mu=0,3)
                        enddo
                        print*,"H(mu)=",(v(mu,2,1),mu=0,3)
                        print*,"SCATTERXXX NCfg ",deltapair,delta,abs(sum/res(k,1)+1)
                        else
!                          print*,"SCATTERXXX NCfg ",deltapair,delta,"NaN"
                      endif
                    endif
                  endif            !debug
               endif !nlo.eq.1
#endif /* WITH_NLO */

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

                  if(nlo.eq.1) then
#ifdef WITH_NLO
                     k=FL_H4j(iflav,diagToPhys,2)
                     if (lok(1)) then
                        res(k,1)=pdf(sign1*iflav(1),1,1)
     &                       *pdf(sign2*iflav(2),2,1)
     &                       *CCmatrix_r(mod(icc1,2),0)*polcolqq
                     else
                        res(k,1) = 0
                     endif
                     res(k,0) = res(k,1)
#endif /* WITH_NLO */
                  else          !loop over color flow
                     do ICOL = 1,6
                        k=FL_H4j(iflav,diagToPhys,ICOL)
                        if (lok(1)) then
                           res(k,1)=pdf(sign1*iflav(1),1,1)
     &                          *pdf(sign2*iflav(2),2,1)
     &                          *CCmatrix_r(mod(icc1,2),ICOL)*polcolqq
                        else
                           res(k,1) = 0
                        endif
                        res(k,0) = res(k,1)
                     enddo
                  endif

#ifdef WITH_NLO
                  if (nlo.eq.1) then
                     do nd = 2,n_kin
                        msq_b(1,nd) = CCmatrix_b(mod(icc1,2),1,nd) !m21_b
                        msq_b(2,nd) = CCmatrix_b(mod(icc1,2),2,nd) !m43_b
                        msq_v(1,nd) = CCmatrix_v(mod(icc1,2),1,nd)
                        msq_v(2,nd) = CCmatrix_v(mod(icc1,2),2,nd)
                     enddo

                     call dipsub(gsign,fsign,msq_b,msq_v,sub,subv,dip) !put in
c gsign to control which dipoles are computed
                  
                     do nd = 2,n_kin
                        if (lok(nd)) then

                           res(k,nd) = -pdf(sign1*iflav(1),1,1)
     $                          *pdf(sign2*iflav(2),2,1)*
     $                          dip(nd)*polcolqq
                  
                           res(k,0) = res(k,0)+res(k,nd)
                        else
                           res(k,nd) = 0
                        endif
                     enddo         !dipole

                     if(ldebug1) then
c     take the limit of pi.pj -> 0
                       delta = dot0p(p(0,5,1),p(0,6,1))
                       deltapair=56
                       newdelta = dot0p(p(0,3,1),p(0,5,1))
                       if (newdelta.lt.delta) then
                         delta=newdelta
                         deltapair=35
                       endif
                       newdelta = dot0p(p(0,3,1),p(0,6,1))
                       if (newdelta.lt.delta) then
                         delta=newdelta
                         deltapair=36
                       endif
                       newdelta = dot0p(p(0,1,1),p(0,6,1))
                       if (newdelta.lt.delta) then
                         delta=newdelta
                         deltapair=16
                       endif
                       newdelta = dot0p(p(0,1,1),p(0,5,1))
                       if (newdelta.lt.delta) then
                         delta=newdelta
                         deltapair=15
                       endif
                       newdelta = dot0p(p(0,4,1),p(0,6,1))
                       if (newdelta.lt.delta) then
                         delta=newdelta
                         deltapair=46
                       endif
                       newdelta = dot0p(p(0,4,1),p(0,5,1))
                       if (newdelta.lt.delta) then
                         delta=newdelta
                         deltapair=45
                       endif
                       newdelta = dot0p(p(0,2,1),p(0,6,1))
                       if (newdelta.lt.delta) then
                         delta=newdelta
                         deltapair=26
                       endif
                       newdelta = dot0p(p(0,2,1),p(0,5,1))
                       if (newdelta.lt.delta) then
                         delta=newdelta
                         deltapair=25
                       endif
                       if(delta.lt.1d0) then
                         if(res(k,1).ne.0d0) then
                           print*,"CC",iflav
                           sum = 0
                           do nd = 2,n_kin
                             sum = sum + res(k,nd)
                             if(res(k,1).ne.0d0) then
                                print*,"D(",ip(nd),kp(nd),jp(nd),
     -                               ")/|m2s|",nd,
     $                               res(k,nd)/res(k,1),lok(nd)
                             else
                                print*,"D(",ip(nd),kp(nd),jp(nd),")",nd
     $                               ,res(k,nd),lok(nd)
                             endif
                           enddo
                           print*,"m2real=",res(k,1),lok(1)
                           print*,"dip_sum=",sum
                           print*,"---------------------------------"
                           print*,"dip_sum/res(k,1)",abs(sum/res(k,1))
                           print*,"----------------------------------"
                           print*,"momentum"
                           do i = 1,n_p
                           print*,"p(mu,",i,")=",(p(mu,i,1),mu=0,3)
                           enddo
                           print*,"H(mu)=",(v(mu,2,1),mu=0,3)
                           print*,"SCATTERXXX CCfg ",deltapair,delta,abs(sum/res(k,1)+1)
                           else
!                             print*,"SCATTERXXX CCfg ",deltapair,delta,"NaN"
                         endif
                       endif
                     endif         !debugging   
                  endif !nlo.eq.1
#endif /* WITH_NLO */

               endif            !cc section
            enddo               !if2
         enddo                  !if1
     
      elseif((gsign(1).eq.-1)) then   
c     This should also work for two initial gluons.
c     create dipole subtraction terms
c
         if(gsign(2).eq.1) then !one i.s. gluon
            polcolab = polcolgq
            Nmx = 12
         elseif(gsign(2).eq.-1) then !two i.s. gluons
            polcolab = polcolgg
            Nmx = 4
         endif
c
         
#ifdef WITH_NLO
         if (nlo.eq.1) then
            do nd = 2,n_kin          !initialize to zero
               do j = 2,3
                  NCmatrix_b(0,0,j-1,nd)=0d0 !|M_b|^2
                  NCmatrix_b(0,1,j-1,nd)=0d0
                  NCmatrix_b(1,0,j-1,nd)=0d0
                  NCmatrix_b(1,1,j-1,nd)=0d0
                  CCmatrix_b(0,j-1,nd)=0d0
                  CCmatrix_b(1,j-1,nd)=0d0

                  NCmatrix_v(0,0,j-1,nd)=0d0 !|M_v|^2
                  NCmatrix_v(0,1,j-1,nd)=0d0
                  NCmatrix_v(1,0,j-1,nd)=0d0
                  NCmatrix_v(1,1,j-1,nd)=0d0
                  CCmatrix_v(0,j-1,nd)=0d0
                  CCmatrix_v(1,j-1,nd)=0d0
               enddo
            enddo                 
         
            call create_ikj(ip,kp,jp)
            do nd = 2,n_kin          !initialize sub,subv
               call dips(nd,p,xuz,ip(nd),kp(nd),jp(nd),sub(1,nd),
     $              subv(nd),qvec(0,nd))
            enddo
 
            do jd = 1,Nmx          !Nmx dipoles set up a map numdip(j)
               if((fsign(1).eq.-fsign(2)).and.(fsign(3).eq.fsign(4))) then
                  nd = numdip21(jd) !define jd gluon emission of 21 line
c     define sign factors for born matrix elements squared
c     in the process of editing
c     for this to work q2 is the quark and q1 is the antiquark
c     
c     decay momentum
                  do mu = 0,3
                     pbar(mu,5) = 0 ! dummy momentum
                     pbar(mu,6) = 0
                     pbar(mu,7) = 0
                     pbar(mu,8) = 0
                     do i = 1,n_v
                        pbar(mu,6) = pbar(mu,6) + v(mu,i,nd) ! Higgs moment.
                     enddo
                  enddo
                  fc(5) = 0
                  fc(6) = 1
                  fc(7) = 0
                  fc(8) = 0

                  if(jd.le.8) then !gQ-> (q2q1~)QH 
                     gc = gsign(1) !g1 is always the i.s. gluon
                     do i=1,4
                        fc(i) = fsign(i)
                     enddo

                     do mu =0,3
                        do i=1,5   !6 is the f.s. soft/collinear gluon
                           pbar(mu,phystodiag(i)) = p(mu,i,nd)
                        enddo
                        qbar(mu,1) = pbar(mu,5)
                     enddo
                            
                  elseif(jd.le.10) then ! qQ->qQgH (f1f(3/4) -> f2f(4/3))
                     gc=gsign(2) !g2 is always the f.s.gluon
                     do mu = 0,3
                        pbar(mu,phystodiag(2)) = p(mu,2,nd) !43 quark line
                        pbar(mu,phystodiag(4)) = p(mu,4,nd) !43 quark line
                        qbar(mu,1) = p(mu,5,nd) !hard gluon is the last entry
                        pbar(mu,1) = p(mu,1,nd) !q1 is i.s. quark
                        pbar(mu,phystodiag(3)) = p(mu,3,nd) !q2 is f.s. quark
                     enddo

                     do i=1,4
                        fc(i)=fsign(i)
                     enddo
                     fc(1)=-fsign(1)
                  
                  elseif(jd.le.12) then ! q~Q ->q~QgH (f2f(3/4)->f1f(4/3))
                     gc=gsign(2)   !g2 is always the f.s.gluon
                     do mu = 0,3
                        pbar(mu,phystodiag(2)) = p(mu,2,nd) !43 quark line
                        pbar(mu,phystodiag(4)) = p(mu,3,nd) !43 quark line
                        qbar(mu,1) = p(mu,5,nd) !hard gluon is the last entry
                        pbar(mu,2) = p(mu,1,nd) !q2 is i.s. quark
                        pbar(mu,phystodiag(5)) = p(mu,4,nd) !q1 is f.s. quark
                     enddo          
 
                     do i=1,4
                        fc(i)=fsign(i)
                     enddo
                     fc(2)=-fsign(2)
                  else
                     stop
                  endif
               
               
               elseif((fsign(3).eq.-fsign(4)).and.(fsign(1).eq.fsign(2))) then 
                  nd = numdip43(jd) !gluon emission of 43 line
c  decay momenta
                  do mu = 0,3
                     pbar(mu,5) = 0 ! dummy momentum
                     pbar(mu,6) = 0
                     pbar(mu,7) = 0
                     pbar(mu,8) = 0
                     do i = 1,n_v
                        pbar(mu,6) = pbar(mu,6) + v(mu,i,nd) ! Higgs moment.
                     enddo
                  enddo
                  fc(5) = 0
                  fc(6) = 1
                  fc(7) = 0
                  fc(8) = 0

                  if(jd.le.8) then !qg-> q Q3Q4~H 
                     gc = gsign(1) !g1 is always the i.s. gluon
                     do i=1,4
                        fc(i) = fsign(i)
                     enddo

                     do mu =0,3
                        do i=1,5   !6 is the f.s. soft/collinear gluon
                           pbar(mu,phystodiag(i)) = p(mu,i,nd)
                        enddo
                        qbar(mu,1) = pbar(mu,5)                   
                     enddo
                  
                  elseif(jd.le.10) then ! qQ->qQgH (f(1/2)f3 -> f(2/1)f4)
                     gc=gsign(2) !g2 is always the f.s.gluon
                     do mu = 0,3
                        pbar(mu,phystodiag(1)) = p(mu,1,nd) !21 quark line
                        pbar(mu,phystodiag(3)) = p(mu,3,nd) !21 quark line
                        qbar(mu,1) = p(mu,5,nd) !hard gluon is the last entry
                        pbar(mu,3) = p(mu,2,nd) !q3 is i.s. quark
                        pbar(mu,phystodiag(4)) = p(mu,4,nd) !q4 is f.s. quark
                     enddo

                     do i=1,4
                        fc(i)=fsign(i)
                     enddo
                     fc(3)=-fsign(3)
                  
                  elseif(jd.le.12) then ! qQ~ ->qQ~gH (f(2/1)f4->f(1/2)f3)
                     gc=gsign(2)   !g2 is always the f.s.gluon
                     do mu = 0,3
                        pbar(mu,phystodiag(1)) = p(mu,1,nd) !21 quark line
                        pbar(mu,phystodiag(3)) = p(mu,3,nd) !21 quark line
                        qbar(mu,1) = p(mu,5,nd) !f.s. gluon is the last entry
                        pbar(mu,phystodiag(5)) = p(mu,4,nd) !q3 is f.s. quark
                        pbar(mu,4) = p(mu,2,nd) !q4 is i.s. quark
                     enddo

                     do i=1,4
                        fc(i)=fsign(i)
                     enddo
                     fc(4)=-fsign(4)
                  else
                     stop
                  endif !jd

               elseif((fsign(1).eq.-fsign(2)).and.(fsign(3).eq.-fsign(4))) then
c     gluons attached to both 21 and 43 line
                  nd = numdip4321(jd)
c     real emmission graphs is
c     g(p1)+g(p2)->q(p3)+Q(p4)+q~(p5)+Q~(p4)+H(P)
c     decay 
                  do mu = 0,3
                     pbar(mu,5) = 0 ! dummy momentum
                     pbar(mu,6) = 0
                     pbar(mu,7) = 0
                     pbar(mu,8) = 0
                     do i = 1,n_v
                        pbar(mu,6) = pbar(mu,6) + v(mu,i,nd) ! Higgs mom.
                     enddo
                  enddo
                  fc(5) = 0
                  fc(6) = 1
                  fc(7) = 0
                  fc(8) = 0

                  if(nd.eq.9) then ! qg
                     gc = gsign(2)
                     do i=1,4
                        fc(i)=fsign(i)
                     enddo
                     fc(1)=-fsign(1) !q1 on i.s.
                     do mu=0,3
                        qbar(mu,1) = p(mu,2,nd)
                        pbar(mu,phystodiag(4))=p(mu,4,nd)
                        pbar(mu,phystodiag(3))=p(mu,3,nd)
                        pbar(mu,phystodiag(6))=p(mu,5,nd)
                        pbar(mu,1)=p(mu,1,nd)
                     enddo

                  elseif(nd.eq.18) then !gQ
                     gc = gsign(1)
                     do i=1,4
                        fc(i)=fsign(i)
                     enddo
                     fc(3)=-fsign(3) !q3 on i.s.
                     do mu=0,3
                        qbar(mu,1) = p(mu,1,nd)
                        pbar(mu,phystodiag(4))=p(mu,4,nd)
                        pbar(mu,phystodiag(3))=p(mu,3,nd)
                        pbar(mu,phystodiag(5))=p(mu,5,nd)
                        pbar(mu,3)=p(mu,2,nd)
                     enddo

                  elseif(nd.eq.22) then !gQ~ 
                     gc = gsign(1)
                     do i=1,4
                        fc(i)=fsign(i)
                     enddo
                     fc(4)=-fsign(4) !q4 on i.s.
                     do mu=0,3
                        qbar(mu,1) = p(mu,1,nd)
                        pbar(mu,phystodiag(3))=p(mu,3,nd)
                        pbar(mu,phystodiag(5))=p(mu,4,nd)
                        pbar(mu,phystodiag(6))=p(mu,5,nd)
                        pbar(mu,4)=p(mu,2,nd)
                     enddo
            
                  elseif(nd.eq.23) then !q~g
                     gc = gsign(2)
                     do i=1,4
                        fc(i)=fsign(i)
                     enddo
                     fc(2)=-fsign(2) !q2 on i.s.
                     do mu=0,3
                        qbar(mu,1) = p(mu,2,nd)
                        pbar(mu,phystodiag(5))=p(mu,4,nd)
                        pbar(mu,phystodiag(4))=p(mu,3,nd)
                        pbar(mu,phystodiag(6))=p(mu,5,nd)
                        pbar(mu,2)=p(mu,1,nd)
                     enddo
                  else
                     stop
                  endif !nd

               endif !fsign

               if (lok(nd)) then

                  call qqhqqj_c_virt(0,pbar(0,1),fc,qbar(0,1),gc,nd,
     $                 uucc_c,uuss_c,ddcc_c,
     $                 ddss_c,udsc_c,ducs_c)
                  if(ldebug) then
                     print*,"---------------"
                     print*,"uucc_c=",uucc_c
                     print*,"gsign=",gc
                     print*,"fsign=",fc
                     do i=1,6   
                        print*,"pbar(",i,")=",(pbar(mu,i),mu=0,3)
                     enddo
                     print*,"qbar=",(qbar(mu,1),mu=0,3)
                     print*,"jd=",jd
                  endif
            
c     j = 1 is gluon emission of 21 line of born
c     j = 2 is gluon emission of 43 line of born    
                  do j = 2,3
                     NCmatrix_b(0,0,j-1,nd)=uucc_c(j) !|M_b|^2
                     NCmatrix_b(0,1,j-1,nd)=uuss_c(j)
                     NCmatrix_b(1,0,j-1,nd)=ddcc_c(j)
                     NCmatrix_b(1,1,j-1,nd)=ddss_c(j)
                     CCmatrix_b(0,j-1,nd)=udsc_c(j)
                     CCmatrix_b(1,j-1,nd)=ducs_c(j)
                  enddo
c     for spin correlations                  
                  call qqHqqj_spcor(qvec(0,nd),pbar(0,1),fc,
     $                 qbar(0,1),gc,nd,uucc_c,uuss_c,
     $                 ddcc_c,ddss_c,udsc_c,ducs_c)
                  do j = 2,3
                     NCmatrix_v(0,0,j-1,nd)=uucc_c(j) !q_mu q_nu M_b^mu M_b*^nu
                     NCmatrix_v(0,1,j-1,nd)=uuss_c(j)
                     NCmatrix_v(1,0,j-1,nd)=ddcc_c(j)
                     NCmatrix_v(1,1,j-1,nd)=ddss_c(j)
                     CCmatrix_v(0,j-1,nd)=udsc_c(j)
                     CCmatrix_v(1,j-1,nd)=ducs_c(j)
                  enddo
               endif !lok(nd)               
            enddo !dipoles
         endif !nlo
#endif /* WITH_NLO */

c     neutral current section
      iflav(diagtophys(5)) = 0  !gluon
      iflav(diagtophys(6)) = 0  !gluon
      do if1=1,nflVBF      !(nfl/2)*2              !21-line
         iflav(diagToPhys(1))=if1*fsign(1)
         iflav(diagToPhys(2))=if1*fsign(2)
         do if2=1,nflVBF   !(nfl/2)*2           !43-line
            iflav(diagToPhys(3))=if2*fsign(3)
            iflav(diagToPhys(4))=if2*fsign(4) 

            if(nlo.eq.1) then
#ifdef WITH_NLO
               k=FL_H4j(iflav,diagToPhys,2)
               if (lok(1)) then
                  res(k,1)=pdf(sign1*iflav(1),1,1) !1 
     &              *pdf(sign2*iflav(2),2,1) !2  
     &                 *NCmatrix_r(mod(if1,2),mod(if2,2),0)*polcolab
               else
                  res(k,1) = 0
               endif   
               res(k,0) = res(k,1)
#endif
            else
               do ICOL=1,6      !loop over color flows
                  k=FL_H4j(iflav,diagToPhys,ICOL)
                  if (lok(1)) then
                     res(k,1)=pdf(sign1*iflav(1),1,1) !1 
     &                    *pdf(sign2*iflav(2),2,1) !2  
     &                    *NCmatrix_r(mod(if1,2),mod(if2,2),ICOL)*
     $                    polcolab
                  else
                     res(k,1) = 0
                  endif   
                  res(k,0) = res(k,1)
               enddo
            endif !nlo

c     dipoles
c     get born and dipoles

#ifdef WITH_NLO
            if (nlo.eq.1) then
               do nd=2,n_kin
                  msq_b(1,nd) = NCmatrix_b(mod(if1,2),
     $                 mod(if2,2),1,nd)
                  msq_b(2,nd) = NCmatrix_b(mod(if1,2),
     $                 mod(if2,2),2,nd)
                  msq_v(1,nd) = NCmatrix_v(mod(if1,2),
     $                 mod(if2,2),1,nd)
                  msq_v(2,nd) = NCmatrix_v(mod(if1,2),
     $                 mod(if2,2),2,nd)  
               enddo
               call dipsub(gsign,fsign,msq_b,msq_v,sub,subv,dip)       

               do nd = 2,n_kin
            
                  if(lok(nd)) then
  
                     if(ldebug) then
                        print*,"NC"
                        print*, "msq_b1(",nd,")=",msq_b(1,nd) 
                        print*, "msq_b2(",nd,")=",msq_b(2,nd)
                        print*,  "msq_v1(",nd,")=",msq_v(1,nd) 
                        print*,  "msq_v2(",nd,")=",msq_v(2,nd)
                        print*,"dip(",nd,")=",dip(nd)
                        print*,"gsign=",gsign
                        print*,"fsign=",fsign
                        print*,"-------------------"
                     endif

                     res(k,nd) = -pdf(sign1*iflav(1),1,1)
     &                    *pdf(sign2*iflav(2),2,1)*dip(nd)/(xi(1)*xi(2))
     
                     res(k,0) = res(k,0)+res(k,nd)             
                  else     
                     res(k,nd) = 0
                  endif !lok(nd)
               enddo !dipoles (nd)

c     debugging 
c     check subtraction
c     c     take the limit of pi.pj -> 0  
               if(ldebug2) then
c     take the limit of pi.pj -> 0
                 delta = 1d99
                 newdelta = abs(dot0p(p(0,2,1),p(0,6,1)))
                 if (newdelta.lt.delta) then
                   delta=newdelta
                   deltapair=26
                 endif
                 if((fsign(1).eq.-fsign(2)).and.(fsign(3).eq.-fsign(4))) then
                   ! gg ->
                   newdelta = abs(dot0p(p(0,1,1),p(0,5,1)))
                   if (newdelta.lt.delta) then
                     delta=newdelta
                     deltapair=15
                   endif
                   newdelta = abs(dot0p(p(0,1,1),p(0,3,1)))
                   if (newdelta.lt.delta) then
                     delta=newdelta
                     deltapair=13
                   endif
                   newdelta = abs(dot0p(p(0,2,1),p(0,4,1)))
                   if (newdelta.lt.delta) then
                     delta=newdelta
                     deltapair=24
                   endif
                 else
                   ! gq/qg ->
                   newdelta = abs(dot0p(p(0,1,1),p(0,6,1)))
                   if (newdelta.lt.delta) then
                     delta=newdelta
                     deltapair=16
                   endif
                   newdelta = abs(dot0p(p(0,3,1),p(0,6,1)))
                   if (newdelta.lt.delta) then
                     delta=newdelta
                     deltapair=36
                   endif
                   newdelta = abs(dot0p(p(0,4,1),p(0,6,1)))
                   if (newdelta.lt.delta) then
                     delta=newdelta
                     deltapair=46
                   endif
                   newdelta = abs(dot0p(p(0,5,1),p(0,6,1)))
                   if (newdelta.lt.delta) then
                     delta=newdelta
                     deltapair=56
                   endif
                   if((fsign(1).eq.-fsign(2)).and.(fsign(3).eq.fsign(4))) then
                     ! gq ->
                     newdelta = abs(dot0p(p(0,1,1),p(0,3,1)))
                     if (newdelta.lt.delta) then
                       delta=newdelta
                       deltapair=13
                     endif
                     newdelta = abs(dot0p(p(0,1,1),p(0,5,1)))
                     if (newdelta.lt.delta) then
                       delta=newdelta
                       deltapair=15
                     endif
                   else
                     ! qg ->
                     newdelta = abs(dot0p(p(0,2,1),p(0,4,1)))
                     if (newdelta.lt.delta) then
                       delta=newdelta
                       deltapair=24
                     endif
                     newdelta = abs(dot0p(p(0,2,1),p(0,5,1)))
                     if (newdelta.lt.delta) then
                       delta=newdelta
                       deltapair=25
                     endif
                   endif
                 endif
                 if(delta.lt.1d0) then
                   if(res(k,1).ne.0d0) then
                     print*,"NC",iflav
                     sum = 0
                     do nd = 2,n_kin
                       sum = sum + res(k,nd)
                       if(res(k,1).ne.0d0) then
                          print*,"D(",ip(nd),kp(nd),jp(nd),
     -                         ")/|m2s|",nd,
     $                         res(k,nd)/res(k,1),lok(nd)
                       else
                          print*,"D(",ip(nd),kp(nd),jp(nd),")",nd
     $                         ,res(k,nd),lok(nd)
                       endif
                     enddo
                     print*,"m2real=",res(k,1),lok(1)
                     print*,"dip_sum=",sum
                     print*,"---------------------------------"
                     print*,"dip_sum/res(k,1)",abs(sum/res(k,1))
                     print*,"----------------------------------"
                     print*,"momentum"
                     do i = 1,n_p
                     print*,"p(mu,",i,")=",(p(mu,i,1),mu=0,3)
                     enddo
                     print*,"H(mu)=",(v(mu,2,1),mu=0,3)
                     print*,"SCATTERXXX NCig ",deltapair,delta,abs(sum/res(k,1)+1)
                    else
!                      print*,"SCATTERXXX NCig ",deltapair,delta,"NaN"
                   endif
                 endif
               endif            !debugging

            endif !nlo.eq.1
#endif /* WITH_NLO */

            enddo               !if2
         enddo                  !if1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
C     Next the CC contributions; no b quarks allowed here
         iflav(diagToPhys(5)) = 0 !gluon
         iflav(diagToPhys(6)) = 0 !gluon
         do if1 = 1,(nflVBF/2)*2   !4           ! sum over all flavors for quark 1 and 2
            iflav(diagToPhys(1)) = fsign(1)*if1
            iflav(diagToPhys(2)) = fsign(2)*(if1-1+2*mod(if1,2))
            do if2 = mod(if1,2)+1,(nflVBF/2)*2,2!4,2 ! flavor of q3 set by q2 (mod 2)
               iflav(diagToPhys(3)) = fsign(3)*if2
               iflav(diagToPhys(4)) = fsign(4)*(if2-1+2*mod(if2,2))

               if(nlo.eq.1) then
#ifdef WITH_NLO
                  k=FL_H4j(iflav,diagToPhys,2)
                  if (lok(1)) then
                     res(k,1)=pdf(sign1*iflav(1),1,1)
     &                   *pdf(sign2*iflav(2),2,1)
     &                    *CCmatrix_r(mod(if1,2),0)*polcolab
                  else
                     res(k,1) = 0.0d0
                  endif
                  res(k,0)=res(k,1)
#endif
               else
                  do ICOL=1,6   !loop over color flows
                     k=FL_H4j(iflav,diagToPhys,ICOL)
                     if (lok(1)) then
                        res(k,1)=pdf(sign1*iflav(1),1,1)
     &                       *pdf(sign2*iflav(2),2,1)
     &                 *CCmatrix_r(mod(if1,2),ICOL)*polcolab  
                     else
                        res(k,1) = 0.0d0
                     endif
                     res(k,0)=res(k,1)
                  enddo
               endif !nlo
                  
c     get born and dipoles

#ifdef WITH_NLO
               if (nlo.eq.1) then 
                  do nd=2,n_kin
                     msq_b(1,nd) = CCmatrix_b(mod(if1,2),1,nd) !m21_b
                     msq_b(2,nd) = CCmatrix_b(mod(if1,2),2,nd) !m43_b
                     msq_v(1,nd) = CCmatrix_v(mod(if1,2),1,nd)
                     msq_v(2,nd) = CCmatrix_v(mod(if1,2),2,nd)
                  enddo
                  call dipsub(gsign,fsign,msq_b,msq_v,sub,subv,dip) 

                  do nd = 2,n_kin
                  
                     if (lok(nd)) then
                     
                        if(ldebug) then
                           print*,"CC"
                           print*, "msq_b1(",nd,")=",msq_b(1,nd) 
                           print*, "msq_b2(",nd,")=",msq_b(2,nd)
                           print*,  "msq_v1(",nd,")=",msq_v(1,nd) 
                           print*,  "msq_v2(",nd,")=",msq_v(2,nd)
                           print*,"dip(",nd,")=",dip(nd)
                           print*,"-------------------"
                        endif
     
                        res(k,nd) = -pdf(sign1*iflav(1),1,1)
     &                       *pdf(sign2*iflav(2),2,1)*dip(nd)/
     &                       (xi(1)*xi(2))     
        
                        res(k,0) = res(k,0)+res(k,nd)
                     else
                        res(k,nd) = 0
                     endif
                  enddo            !dipoles

c     debugging
                  if(ldebug2) then
c     take the limit of pi.pj -> 0
                    delta = 1d99
                    newdelta = abs(dot0p(p(0,2,1),p(0,6,1)))
                    if (newdelta.lt.delta) then
                      delta=newdelta
                      deltapair=26
                    endif
                    if((fsign(1).eq.-fsign(2)).and.(fsign(3).eq.-fsign(4))) then
                      ! gg ->
                      newdelta = abs(dot0p(p(0,1,1),p(0,5,1)))
                      if (newdelta.lt.delta) then
                        delta=newdelta
                        deltapair=15
                      endif
                      newdelta = abs(dot0p(p(0,1,1),p(0,3,1)))
                      if (newdelta.lt.delta) then
                        delta=newdelta
                        deltapair=13
                      endif
                      newdelta = abs(dot0p(p(0,2,1),p(0,4,1)))
                      if (newdelta.lt.delta) then
                        delta=newdelta
                        deltapair=24
                      endif
                    else
                      ! gq/qg ->
                      newdelta = abs(dot0p(p(0,1,1),p(0,6,1)))
                      if (newdelta.lt.delta) then
                        delta=newdelta
                        deltapair=16
                      endif
                      newdelta = abs(dot0p(p(0,3,1),p(0,6,1)))
                      if (newdelta.lt.delta) then
                        delta=newdelta
                        deltapair=36
                      endif
                      newdelta = abs(dot0p(p(0,4,1),p(0,6,1)))
                      if (newdelta.lt.delta) then
                        delta=newdelta
                        deltapair=46
                      endif
                      newdelta = abs(dot0p(p(0,5,1),p(0,6,1)))
                      if (newdelta.lt.delta) then
                        delta=newdelta
                        deltapair=56
                      endif
                      if((fsign(1).eq.-fsign(2)).and.(fsign(3).eq.fsign(4))) then
                        ! gq ->
                        newdelta = abs(dot0p(p(0,1,1),p(0,3,1)))
                        if (newdelta.lt.delta) then
                          delta=newdelta
                          deltapair=13
                        endif
                        newdelta = abs(dot0p(p(0,1,1),p(0,5,1)))
                        if (newdelta.lt.delta) then
                          delta=newdelta
                          deltapair=15
                        endif
                      else
                        ! qg ->
                        newdelta = abs(dot0p(p(0,2,1),p(0,4,1)))
                        if (newdelta.lt.delta) then
                          delta=newdelta
                          deltapair=24
                        endif
                        newdelta = abs(dot0p(p(0,2,1),p(0,5,1)))
                        if (newdelta.lt.delta) then
                          delta=newdelta
                          deltapair=25
                        endif
                      endif
                    endif
                    if(delta.lt.1d0) then
                      if(res(k,1).ne.0d0) then
                        print*,"CC",iflav
                        sum = 0
                        do nd = 2,n_kin
                          sum = sum + res(k,nd)
                          if(res(k,1).ne.0d0) then
                             print*,"D(",ip(nd),kp(nd),jp(nd),
     -                            ")/|m2s|",nd,
     $                            res(k,nd)/res(k,1),lok(nd)
                          else
                             print*,"D(",ip(nd),kp(nd),jp(nd),")",nd
     $                            ,res(k,nd),lok(nd)
                          endif
                        enddo
                        print*,"m2real=",res(k,1),lok(1)
                        print*,"dip_sum=",sum
                        print*,"---------------------------------"
                        print*,"dip_sum/res(k,1)",abs(sum/res(k,1))
                        print*,"----------------------------------"
                        print*,"momentum"
                        do i = 1,n_p
                        print*,"p(mu,",i,")=",(p(mu,i,1),mu=0,3)
                        enddo
                        print*,"H(mu)=",(v(mu,2,1),mu=0,3)
                        print*,"SCATTERXXX CCig ",deltapair,delta,abs(sum/res(k,1)+1)
                      else
!                        print*,"SCATTERXXX CCig ",deltapair,delta,"NaN"
                      endif
                    endif
                  endif            !debugging

               endif !nlo.eq.1
#endif /* WITH_NLO */

            enddo               !if2
         enddo                  !if1   
         
      else
         stop
      endif

      nmax = FL_H4j(iflav,diagTophys,-2)
      return
      end

c******************************************************************
c
c   end subroutine wbf_h4jg
c
c*****************************************************************
c
c******************************************************************
c
c   begin function FL_H4j 
c
c*****************************************************************
      integer function FL_H4j(iflav,diagToPhys,colstruc)

      implicit none

      integer iflav(6),colstruc    ! input for color and flavor assignment
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
 
      integer listposition,diagToPhys(6)
      integer numresets
      save listposition
      save numresets

      data numresets /-1/
      data listposition /0/
      
!       print*, iflav,diagToPhys,colstruc

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
          call fillColoredPartons_H4j(iflav,diagToPhys,
     1                                listposition+numdecay,colstruc)
        endif       ! numresets.eq.0
      endif     
c      print*,"listposition=",listposition
      FL_H4j=listposition
      end
c******************************************************************************
c
c   end function FL_ZHg
c
c******************************************************************************
c*****************************************************************************
c
c    begin  subroutine fillColoredPartons_H4j
c
c*****************************************************************************
      subroutine fillColoredPartons_H4j(iflav,
     $                                  diagToPhys,
     $                                  listposition,
     $                                  colstruc)
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

      integer iflav(6),listposition,i,i1
      integer iflavour(6),colstruc

      integer cid(2,6),diagToPhys(6),Ngluon
      logical ldebug,IsFinalGluon,IsQuark
      parameter (ldebug=.false.)
c...  supply diagToPhys
c...

      Ngluon=0
      do i1=1,6
         iflavour(i1)=iflav(i1) 
         if(iflavour(i1).eq.0) then 
            iflavour(i1)=21
            Ngluon = Ngluon + 1
         endif
      enddo

      select case (process)
      case(HjjjjLO_WW, HjjjjLO_ZZ_ll, HjjjjLO_ZZ_lnu)
         lnup(listposition)=numParticles+4 
      case(HjjjjLO)
         lnup(listposition)=numParticles-1
      case(HjjjjLO_AA, HjjjjLO_mu, HjjjjLO_tau, HjjjjLO_bbar)
         lnup(listposition)=numParticles+1                   
      end select
      
      listup(1,listposition)=-1 !incoming quarks
      listup(2,listposition)=-1
      listup(3,listposition)=1  !outgoing quarks
      listup(4,listposition)=1
      listup(5,listposition)=1
      listup(6,listposition)=1
 
c...four quarks and two gluons
      do i =1,6
         lidup(i,listposition)=iflavour(i)
      enddo
c...
      do i=3,n_p
         lmothup(1,i,listposition)=1
         lmothup(2,i,listposition)=2
      enddo
c... Color flows : (1) T[2,1,5,6]T[4,3]  (2) T[2,1,5] T[4,3,6]
c...             : (3) T[2,1,6,5] T[4,3] (4) T[2,1,6] T[4,3,5]
c...             : (5) T[2,1] T[4,3,5,6] (6) T[2,1] T[4,3,6,5]
c...two final/initial state gluons
      if(Ngluon.eq.2)then

         do i=1,4               !loop over all particles (diagramatic labels)
            IsQuark = (iflavour(diagToPhys(i)).gt.0)
            
            if(IsQuark) then    !quark or anti-quark/ intial gluon 
               cid(1,i)=1       !color
               cid(2,i)=2       !anti-color
            else
               cid(2,i)=1
               cid(1,i)=2
            endif
         enddo
         
         do i=5,6               !loop over all particles (diagramatic labels)
            IsFinalGluon = (diagToPhys(i).gt.2)
         
         if(IsFinalGluon) then !quark or anti-quark/ intial gluon 
            cid(1,i)=1          !color
            cid(2,i)=2          !anti-color
         else
            cid(2,i)=1
            cid(1,i)=2
         endif
      enddo


         if(colstruc.eq.1) then
               licolup(cid(1,1),diagToPhys(1),listposition)=501
               licolup(cid(1,2),diagToPhys(2),listposition)=503
               licolup(cid(1,5),diagToPhys(5),listposition)=502
               licolup(cid(1,6),diagToPhys(6),listposition)=501
               licolup(cid(2,1),diagToPhys(1),listposition)=0
               licolup(cid(2,2),diagToPhys(2),listposition)=0
               licolup(cid(2,5),diagToPhys(5),listposition)=503
               licolup(cid(2,6),diagToPhys(6),listposition)=502
                        
               licolup(cid(1,3),diagToPhys(3),listposition)=504
               licolup(cid(1,4),diagToPhys(4),listposition)=504
               licolup(cid(2,3),diagToPhys(3),listposition)=0
               licolup(cid(2,4),diagToPhys(4),listposition)=0
        
         elseif(colstruc.eq.2) then
        
               licolup(cid(1,1),diagToPhys(1),listposition)=501
               licolup(cid(1,2),diagToPhys(2),listposition)=502
               licolup(cid(1,5),diagToPhys(5),listposition)=501
               licolup(cid(2,1),diagToPhys(1),listposition)=0
               licolup(cid(2,2),diagToPhys(2),listposition)=0
               licolup(cid(2,5),diagToPhys(5),listposition)=502
              
               licolup(cid(1,3),diagToPhys(3),listposition)=503
               licolup(cid(1,4),diagToPhys(4),listposition)=504
               licolup(cid(1,6),diagToPhys(6),listposition)=503
               licolup(cid(2,3),diagToPhys(3),listposition)=0
               licolup(cid(2,4),diagToPhys(4),listposition)=0
               licolup(cid(2,6),diagToPhys(6),listposition)=504

         elseif(colstruc.eq.3) then !
           
            licolup(cid(1,1),diagToPhys(1),listposition)=501
            licolup(cid(1,2),diagToPhys(2),listposition)=503
            licolup(cid(1,6),diagToPhys(6),listposition)=502
            licolup(cid(1,5),diagToPhys(5),listposition)=501
            licolup(cid(2,1),diagToPhys(1),listposition)=0
            licolup(cid(2,2),diagToPhys(2),listposition)=0
            licolup(cid(2,6),diagToPhys(6),listposition)=503
            licolup(cid(2,5),diagToPhys(5),listposition)=502
                     
            licolup(cid(1,3),diagToPhys(3),listposition)=504
            licolup(cid(1,4),diagToPhys(4),listposition)=504
            licolup(cid(2,3),diagToPhys(3),listposition)=0
            licolup(cid(2,4),diagToPhys(4),listposition)=0
         
         elseif(colstruc.eq.4) then !
         
            licolup(cid(1,1),diagToPhys(1),listposition)=501
            licolup(cid(1,2),diagToPhys(2),listposition)=502
            licolup(cid(1,6),diagToPhys(6),listposition)=501
            licolup(cid(2,1),diagToPhys(1),listposition)=0
            licolup(cid(2,2),diagToPhys(2),listposition)=0
            licolup(cid(2,6),diagToPhys(6),listposition)=502
            
            licolup(cid(1,3),diagToPhys(3),listposition)=503
            licolup(cid(1,4),diagToPhys(4),listposition)=504
            licolup(cid(1,5),diagToPhys(5),listposition)=503
            licolup(cid(2,3),diagToPhys(3),listposition)=0
            licolup(cid(2,4),diagToPhys(4),listposition)=0
            licolup(cid(2,5),diagToPhys(5),listposition)=504
            
         elseif(colstruc.eq.5) then !
           
            licolup(cid(1,3),diagToPhys(3),listposition)=501
            licolup(cid(1,4),diagToPhys(4),listposition)=503
            licolup(cid(1,5),diagToPhys(5),listposition)=502
            licolup(cid(1,6),diagToPhys(6),listposition)=501
            licolup(cid(2,3),diagToPhys(3),listposition)=0
            licolup(cid(2,4),diagToPhys(4),listposition)=0
            licolup(cid(2,5),diagToPhys(5),listposition)=503
            licolup(cid(2,6),diagToPhys(6),listposition)=502
           
            licolup(cid(1,1),diagToPhys(1),listposition)=504
            licolup(cid(1,2),diagToPhys(2),listposition)=504
            licolup(cid(2,1),diagToPhys(1),listposition)=0
            licolup(cid(2,2),diagToPhys(2),listposition)=0
            
         elseif(colstruc.eq.6) then
            licolup(cid(1,3),diagToPhys(3),listposition)=501
            licolup(cid(1,4),diagToPhys(4),listposition)=503
            licolup(cid(1,6),diagToPhys(6),listposition)=502
            licolup(cid(1,5),diagToPhys(5),listposition)=501
            licolup(cid(2,3),diagToPhys(3),listposition)=0
            licolup(cid(2,4),diagToPhys(4),listposition)=0
            licolup(cid(2,6),diagToPhys(6),listposition)=503
            licolup(cid(2,5),diagToPhys(5),listposition)=502
            
            licolup(cid(1,1),diagToPhys(1),listposition)=504
            licolup(cid(1,2),diagToPhys(2),listposition)=504
            licolup(cid(2,1),diagToPhys(1),listposition)=0
            licolup(cid(2,2),diagToPhys(2),listposition)=0
           
         endif

      elseif(Ngluon.eq.0) then  !6 quark 
         do i=1,6               !loop over all particles (diagramatic labels)
            IsQuark = (iflavour(diagToPhys(i)).gt.0)
         
            if(IsQuark) then    !quark or anti-quark/ intial gluon 
               cid(1,i)=1       !color
               cid(2,i)=2       !anti-color
            else
               cid(2,i)=1
               cid(1,i)=2
            endif
         enddo
         
         if(colstruc.eq.1) then ! 
            licolup(cid(1,1),diagToPhys(1),listposition)=501
            licolup(cid(1,2),diagToPhys(2),listposition)=502
            licolup(cid(1,6),diagToPhys(6),listposition)=501
            licolup(cid(1,5),diagToPhys(5),listposition)=502!0
c
            licolup(cid(2,1),diagToPhys(1),listposition)=0
            licolup(cid(2,2),diagToPhys(2),listposition)=0
            licolup(cid(2,6),diagToPhys(6),listposition)=0
            licolup(cid(2,5),diagToPhys(5),listposition)=0!502
            
            licolup(cid(1,3),diagToPhys(3),listposition)=503
            licolup(cid(1,4),diagToPhys(4),listposition)=503
            licolup(cid(2,3),diagToPhys(3),listposition)=0
            licolup(cid(2,4),diagToPhys(4),listposition)=0

         elseif(colstruc.eq.2) then

            licolup(cid(1,3),diagToPhys(3),listposition)=501
            licolup(cid(1,4),diagToPhys(4),listposition)=502
            licolup(cid(1,6),diagToPhys(6),listposition)=501
            licolup(cid(1,5),diagToPhys(5),listposition)=502!0
c
            licolup(cid(2,3),diagToPhys(3),listposition)=0
            licolup(cid(2,4),diagToPhys(4),listposition)=0
            licolup(cid(2,6),diagToPhys(6),listposition)=0
            licolup(cid(2,5),diagToPhys(5),listposition)=0!502
            
            licolup(cid(1,1),diagToPhys(1),listposition)=503
            licolup(cid(1,2),diagToPhys(2),listposition)=503
            licolup(cid(2,1),diagToPhys(1),listposition)=0
            licolup(cid(2,2),diagToPhys(2),listposition)=0

         endif

      endif


      if(ldebug)then
         print*,"colstruc = ", colstruc
         print*,"====================================="
         print("(1P,6I6)"),(iflavour(i1),i1=1,6)
         print*,"-------------------------------------"
         print("(1P,6I6)"),(licolup(1,i1,listposition),i1=1,6)
         print("(1P,6I6)"),(licolup(2,i1,listposition),i1=1,6)
         print*,"====================================="
         print*,""
      endif
    
      end  ! fillColoredPartons_H4j
c*****************************************************************************
c
c    end subroutine fillColoredPartons_H4j
c
c*****************************************************************************



