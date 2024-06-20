c     This subroutine computes ALL the virtual corrections along the 
c     quark line depicted below and returns it as a Lorentz vector. Only 
c     TWO contractions with external currents (among jmu1, jmu2, jmu3) are done
c     
c
c      q,jmu1            t,jmu3
c        |                 |
c       \|/               \|/
c        |                 |
c -->-------------------------->---
c   p            |             p"
c               /|\
c                |
c              l,jmu2
c
c     p^2 = p"^2 = 0
c     
c     The notation for the external momenta is
c     
c     p + q + l + t = p" = pp
c     musq = mu^2 = common external scale (renormalization scale, eventually)
c     
c  
c     Only TWO external currents are passed to the subroutine
c     
c     Carlo Oleari    13/6/2005
c
      subroutine J_virtual_tri_box_pent(
c     bra(p"), ket(p) and helicity (of the ket fermion)
     &     brapp,ketp,hel,
c     momenta and squared renormalization scale
     &     p,q,l,t,pp,musq,
c     input currents (any couple among [J1,J2,J3]). ORDER IS IMPORTANT!!
c     They can be  jmu1,jmu2  or  jmu1,jmu3   or   jmu2,jmu3
     &     jal,jbe,
c     uncontracted external current. It can be 1, 2 or 3
     &     num_ret_current,
c     returned current
     &     Jret)

c
c     modified by Guiseppe Bozzi and Vera Hankele to include WWA code
c

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"

      double complex brapp(2,-1:1),ketp(2,-1:1)
      integer hel
      double precision p(0:3),pp(0:3),q(0:3),l(0:3),t(0:3)
      double precision musq
      double complex jal(0:3),jbe(0:3)
      integer num_ret_current
      double complex jret(0:3)

      double complex jmu1(0:3),jmu2(0:3),jmu3(0:3)
      double precision psq, qsq, lsq, tsq,
     &     pq, pl, pt, ql, qt, lt, ppp, ppq, ppl
      double precision pq_sq,pql_sq,pqlt_sq,ql_sq,qlt_sq,lt_sq
      double complex C0fin,B0fin
      double complex E02
      double complex B0_12,B0_13,B0_14,B0_15,B0_23,B0_24,B0_25,B0_34,
     &     B0_35,B0_45
      double complex C0_345,C0_245,C0_235,C0_234,C0_145,C0_135,C0_134,
     &     C0_125,C0_124,C0_123
      double complex Cij_345(2,4),Cij_245(2,4),Cij_235(2,4),Cij_234(2,4),
     &     Cij_145(2,4),Cij_135(2,4),Cij_134(2,4),Cij_125(2,4),
     &     Cij_124(2,4),Cij_123(2,4)
      double complex D0_2345,D0_1345,D0_1245,D0_1235,D0_1234
      double complex Dij_2345(3,13),Dij_1345(3,13),Dij_1245(3,13),
     &     Dij_1235(3,13),Dij_1234(3,13),D0fin
      double complex E0_12345, Eij_12345(4,46)
      common/E_functions/E0_12345, Eij_12345
      common/D_functions/D0_2345,D0_1345,D0_1245,D0_1235,D0_1234,
     &     Dij_2345,Dij_1345,Dij_1245,Dij_1235,Dij_1234
      common/C_functions/
     &     C0_345,C0_245,C0_235,C0_234,C0_145,
     &     C0_135,C0_134,C0_125,C0_124,C0_123,
     &     Cij_345,Cij_245,Cij_235,Cij_234,Cij_145,
     &     Cij_135,Cij_134,Cij_125,Cij_124,Cij_123
      common/B_functions/B0_13,B0_14,B0_24,B0_34,B0_23,B0_45,B0_35
      integer compon
      double complex tri_box,virt_tri_box,pent,virt_pent


      logical debug
c      parameter (debug=.true.)
      parameter (debug=.false.)

      double complex jref(0:3,0:3)
      integer mu
      double precision dotrr
      external D0fin

      psq = 0d0                 ! dotrr(p,p)
      qsq = dotrr(q,q)
      lsq = dotrr(l,l)
      tsq = dotrr(t,t)

      pq = dotrr(p,q)
      pl = dotrr(p,l)
      ql = dotrr(q,l)
      qt = dotrr(q,t)
      lt = dotrr(l,t)

      pt = dotrr(p,t)

      ppp = dotrr(p,pp)
      ppq = dotrr(pp,q)
      ppl = dotrr(pp,l)

c (p+q)^2
      pq_sq = qsq+2*pq+psq
c (p+q+l)^2
      pql_sq = lsq+2*ql+2*pl+qsq+2*pq+psq
c (p+q+l+t)^2
      pqlt_sq =  0d0 ! tsq+2*lt+2*qt+2*pt+lsq+2*ql+2*pl+qsq+2*pq+psq
c (q+l)^2
      ql_sq = lsq+2*ql+qsq  
c (q+l+t)^2
      qlt_sq = tsq+2*lt+2*qt+lsq+2*ql+qsq
c (l+t)^2
      lt_sq = tsq+2*lt+lsq

      if (debug) then
         write(*,*) "more prod scal"
         write(*,*) psq, qsq, lsq, tsq, pq, pl, ql, qt, lt, pt, ppp, ppq, ppl
         write(*,*) pq_sq,pql_sq,pqlt_sq,ql_sq,qlt_sq,lt_sq
      endif

c B0(p)
      B0_12 = 0d0 ! B0fin(psq,musq)
c B0(p+q)
      B0_13 = B0fin(pq_sq,musq)
c B0(p+q+l)
      B0_14 = B0fin(pql_sq,musq)
c B0(p+q+l+t)
      B0_15 = B0fin(pqlt_sq,musq)
c B0(q)
      B0_23 = B0fin(qsq,musq)
c B0(q+l)
      B0_24 = B0fin(ql_sq,musq)
c B0(q+l+t)
      B0_25 = B0fin(qlt_sq,musq)
c B0(l)
      B0_34 = B0fin(lsq,musq)
c B0(l+t)
      B0_35 = B0fin(lt_sq,musq)
c B0(t)
      B0_45 = B0fin(tsq,musq)

c C0(l,t)
      C0_345 = C0fin(lsq,tsq,lt_sq,musq)
c C0(q+l,t)
      C0_245 = C0fin(ql_sq,tsq,qlt_sq,musq)
c C0(q,l+t)
      C0_235 = C0fin(qsq,lt_sq,qlt_sq,musq)
c C0(q,l)
      C0_234 = C0fin(qsq,lsq,ql_sq,musq)
c C0(p+q+l,t)
      C0_145 = C0fin(pql_sq,tsq,pqlt_sq,musq)
c C0(p+q,l+t)
      C0_135 = C0fin(pq_sq,lt_sq,pqlt_sq,musq)
c C0(p+q,l)
      C0_134 = C0fin(pq_sq,lsq,pql_sq,musq)
c C0(p,q+l+t)
      C0_125 = C0fin(psq,qlt_sq,pqlt_sq,musq)
c C0(p,q+l)
      C0_124 = C0fin(psq,ql_sq,pql_sq,musq)
c C0(p,q)
      C0_123 = C0fin(psq,qsq,pq_sq,musq)

c Cij_345
      call tens_red3(lsq,tsq,lt_sq, B0_45, B0_35, B0_34, 
     &     C0_345, Cij_345)
c Cij_245
      call tens_red3(ql_sq,tsq,qlt_sq, B0_45, B0_25, B0_24, 
     &     C0_245, Cij_245)
c Cij_235
      call tens_red3(qsq,lt_sq,qlt_sq, B0_35, B0_25, B0_23, 
     &     C0_235, Cij_235)
c Cij_234
      call tens_red3(qsq,lsq,ql_sq, B0_34, B0_24, B0_23, 
     &     C0_234, Cij_234)
c Cij_145
      call tens_red3(pql_sq,tsq,pqlt_sq, B0_45, B0_15, B0_14, 
     &     C0_145, Cij_145)
c Cij_135
      call tens_red3(pq_sq,lt_sq,pqlt_sq, B0_35, B0_15, B0_13, 
     &     C0_135, Cij_135)
c Cij_134
      call tens_red3(pq_sq,lsq,pql_sq, B0_34, B0_14, B0_13, 
     &     C0_134, Cij_134)
c Cij_125
      call tens_red3(psq,qlt_sq,pqlt_sq, B0_25, B0_15, B0_12, 
     &     C0_125, Cij_125)
c Cij_124
      call tens_red3(psq,ql_sq,pql_sq, B0_24, B0_14, B0_12, 
     &     C0_124, Cij_124)
c Cij_123
      call tens_red3(psq,qsq,pq_sq, B0_23, B0_13, B0_12, 
     &     C0_123, Cij_123)

c      call print_Cij(C0_123,Cij_123)
c      call print_Cij(C0_245,Cij_245)

c      SELECT CASE(procID)
c NEW NEW NEW NEW - has to be checked for every sigle process
c      CASE(WWA)
      D0_2345 = D0fin(ql_sq,lt_sq,qsq,lsq,tsq,qlt_sq,musq)
      D0_1345 = D0fin(pql_sq,lt_sq,pq_sq,lsq,tsq,pqlt_sq,musq)
      D0_1245 = D0fin(pql_sq,qlt_sq,psq,ql_sq,tsq,pqlt_sq,musq)
      D0_1235 = D0fin(pq_sq,qlt_sq,psq,qsq,lt_sq,pqlt_sq,musq)
      D0_1234 = D0fin(pq_sq,ql_sq,psq,qsq,lsq,pql_sq,musq)

c      CASE DEFAULT
c OLD OLD OLD - works for everything except WWA code
c      D0_2345 = D04m_fin(lt_sq,ql_sq,qlt_sq,qsq,lsq,tsq,musq)
c      D0_1345 = D03m_fin(lt_sq,pql_sq,pq_sq,lsq,tsq,musq)
c      D0_1245 = D02m_fin(qlt_sq,pql_sq,ql_sq,tsq,musq)
c      D0_1235 = D02m_fin(qlt_sq,pq_sq,qsq,lt_sq,musq)
c      D0_1234 = D03m_fin(ql_sq,pq_sq,pql_sq,lsq,qsq,musq)

c      print*,'D0_2345 = ',D0_2345

c      test1 = D04m_fin(lt_sq,ql_sq,qlt_sq,qsq,lsq,tsq,musq)
c      test2 = D03m_fin(lt_sq,pql_sq,pq_sq,lsq,tsq,musq)
c      test3 = D02m_fin(qlt_sq,pql_sq,ql_sq,tsq,musq)
c      test4 = D02m_fin(qlt_sq,pq_sq,qsq,lt_sq,musq)
c      test5 = D03m_fin(ql_sq,pq_sq,pql_sq,lsq,qsq,musq)

c      END SELECT

c      temp = dreal(test1)-dreal(D0_2345)
c      test = dimag(test1)-dimag(D0_2345)
c      if (temp.gt.1.0d-15 .or. test.gt.1.0d-13) then
c         print*,'D0_2345 (old)= ',test1
c         print*,'D0_2345 (new)= ',D0_2345
c         pause
c      endif
c      temp = dreal(test2)-dreal(D0_1345)
c      test = dimag(test2)-dimag(D0_1345)
c      if (temp.gt.1.0d-15 .or. test.gt.1.0d-13) then
c         print*,'D0_1345 (old)= ',test2
c         print*,'D0_1345 (new)= ',D0_1345
c         pause
c      endif
c      temp = dreal(test3)-dreal(D0_1245)
c      test = dimag(test3)-dimag(D0_1245)
c      if (temp.gt.1.0d-15 .or. test.gt.1.0d-13) then
c         print*,'D0_1245 (old)= ',test3
c         print*,'D0_1245 (new)= ',D0_1245
c         pause
c      endif
c      temp = dreal(test4)-dreal(D0_1235)
c      test = dimag(test4)-dimag(D0_1235)
c      if (temp.gt.1.0d-15 .or. test.gt.1.0d-13) then
c         print*,'D0_1235 (old)= ',test4
c         print*,'D0_1235 (new)= ',D0_1235
c         pause
c      temp = dreal(test5)-dreal(D0_1234)
c      test = dimag(test5)-dimag(D0_1234)
c      if (temp.gt.1.0d-15 .or. test.gt.1.0d-13) then
c      endif
c         print*,'D0_1234 (old)= ',test5
c         print*,'D0_1234 (new)= ',D0_1234
c         pause
c      endif

c$$$c D0_2345 = D00(q,l,t)
c$$$c      D0_2345 = D0(qsq,lsq,tsq,ql,qt,lt)
c$$$c the following two expressions are equivalent
c$$$      D0_2345 = D04m_fin(lt_sq,ql_sq,qlt_sq,qsq,lsq,tsq,musq)
c$$$c      write(*,*) "D0_2345",D0_2345
c$$$c      D0_2345 = D04m_fin(ql_sq,lt_sq,qsq,lsq,tsq,qlt_sq,musq)
c$$$c      write(*,*) "D0_2345",D0_2345
c$$$      
c$$$c D0_1345 = D00(p+q,l,t)         
c$$$c      D0_1345 = D0(pq_sq,lsq,tsq,ql+pl,qt+pt,lt)
c$$$      D0_1345 = D03m_fin(lt_sq,pql_sq,pq_sq,lsq,tsq,musq)
c$$$c      write(*,*) "D0_1345",D0_1345
c$$$
c$$$c D0_1245 = D00(p,q+l,t) 
c$$$c      D0_1245 = D0(psq,ql_sq,tsq,pl+pq,pt,lt+qt)
c$$$      D0_1245 = D02m_fin(qlt_sq,pql_sq,ql_sq,tsq,musq)
c$$$c      write(*,*) "D0_1245",D0_1245
c$$$
c$$$c D0_1235 = D00(p,q,l+t)  
c$$$c      D0_1235 = D0(psq,qsq,lt_sq,pq,pt+pl,qt+ql)
c$$$      D0_1235 = D02m_fin(qlt_sq,pq_sq,qsq,lt_sq,musq)
c$$$c      write(*,*) "D0_1235",D0_1235
c$$$
c$$$c D0_1234 = D00(p,q,l)
c$$$c      D0_1234 = D0(psq,qsq,lsq,pq,pl,ql)
c$$$      D0_1234 = D03m_fin(ql_sq,pq_sq,pql_sq,lsq,qsq,musq)
c$$$c      write(*,*) "D0_1234",D0_1234

c Dij_2345
      call tens_red4(qsq,lsq,tsq,ql,qt,lt,
     &     C0_345, C0_245, C0_235, C0_234,   
     &     Cij_345, Cij_245, Cij_235, Cij_234,   
     &     D0_2345, 
     &     Dij_2345)

c Dij_1345
      call tens_red4(pq_sq,lsq,tsq,ql+pl,qt+pt,lt,
     &     C0_345, C0_145, C0_135, C0_134,   
     &     Cij_345, Cij_145, Cij_135, Cij_134,   
     &     D0_1345, 
     &     Dij_1345)

c Dij_1245
      call tens_red4(psq,ql_sq,tsq,pl+pq,pt,lt+qt,
     &     C0_245, C0_145, C0_125, C0_124,   
     &     Cij_245, Cij_145, Cij_125, Cij_124,   
     &     D0_1245, 
     &     Dij_1245)

c Dij_1235
      call tens_red4(psq,qsq,lt_sq,pq,pt+pl,qt+ql,
     &     C0_235, C0_135, C0_125, C0_123,   
     &     Cij_235, Cij_135, Cij_125, Cij_123,   
     &     D0_1235, 
     &     Dij_1235)

c Dij_1234
      call tens_red4(psq,qsq,lsq,pq,pl,ql,
     &     C0_234, C0_134, C0_124, C0_123,   
     &     Cij_234, Cij_134, Cij_124, Cij_123,   
     &     D0_1234, 
     &     Dij_1234)

c      call print_Dij(D0_1234,Dij_1234)
c      call print_Dij(D0_2345,Dij_2345)      
c      call print_Dij(D0_1345,Dij_1345)      
c      call print_Dij(D0_1245,Dij_1245)      
c      call print_Dij(D0_1235,Dij_1235)      

c E0_12345 = E02(p,q,l,t)
      E0_12345 = E02(psq, qsq, lsq, tsq,
     &     pq, pl, pt, ql, qt, lt, 
     &     D0_2345,D0_1345,D0_1245,D0_1235,D0_1234)
c      write(*,*) "E02",E0_12345

c Eij_12345
      call tens_red5(psq, qsq, lsq, tsq,
     &     pq, pl, pt, ql, qt, lt, 
     &     D0_2345,D0_1345,D0_1245,D0_1235,D0_1234,
     &     Dij_2345,Dij_1345,Dij_1245,Dij_1235,Dij_1234,
     &     E0_12345, 
     &     Eij_12345)

      if (debug) then
         write(*,*) "E0_12345" 
         write(*,*) E0_12345
         write(*,*) "D0_2345,D0_1345,D0_1245,D0_1235,D0_1234"
         write(*,*) D0_2345,D0_1345,D0_1245,D0_1235,D0_1234
         write(*,*) "C0_345,C0_245,C0_235,C0_234,C0_145,"//
     &        "C0_135,C0_134,C0_125,C0_124,C0_123"
         write(*,*) C0_345,C0_245,C0_235,C0_234,C0_145,
     &        C0_135,C0_134,C0_125,C0_124,C0_123
         write(*,*) "B0_12,B0_13,B0_14,B0_15,B0_23,B0_24,B0_25,B0_34,"//
     &        "B0_35,B0_45"
         write(*,*) B0_12,B0_13,B0_14,B0_15,B0_23,B0_24,B0_25,B0_34,
     &        B0_35,B0_45
      endif
      
c      call print_Cij(C0_123,Cij_123)
c      call print_Dij(D0_2345,Dij_2345)
c      call print_Eij(E0_12345,Eij_12345)            
c      write(*,*) qsq,lsq,pq,pl,ppp,ppq,ppl,ql
      
      jref(0,0)=1d0 
      jref(0,1)=0d0
      jref(0,2)=0d0
      jref(0,3)=0d0
      
      jref(1,0)=0d0
      jref(1,1)=-1d0
      jref(1,2)=0d0
      jref(1,3)=0d0
      
      jref(2,0)=0d0
      jref(2,1)=0d0
      jref(2,2)=-1d0
      jref(2,3)=0d0

      jref(3,0)=0d0
      jref(3,1)=0d0
      jref(3,2)=0d0
      jref(3,3)=-1d0
      
      if  (num_ret_current.eq.1) then
         do mu=0,3
            jmu2(mu)=jal(mu)
            jmu3(mu)=jbe(mu)
         enddo
         do compon=0,3
            do mu=0,3
               jmu1(mu) = jref(compon,mu)
            enddo
            call compute_Jprod_gammaprod(
     &           brapp,ketp,hel,
     &           p,q,l,t,pp,
     &           jmu1,jmu2,jmu3)
            tri_box = virt_tri_box(qsq,lsq,pq,pl,ppp,ppq,ppl,ql)
            pent    =    virt_pent(qsq,lsq,pq,pl,ppp,ppq,ppl,ql)
            Jret(compon) = tri_box + pent
         enddo
      elseif  (num_ret_current.eq.2) then
         do mu=0,3
            jmu1(mu)=jal(mu)
            jmu3(mu)=jbe(mu)
         enddo
         do compon=0,3
            do mu=0,3
               jmu2(mu) = jref(compon,mu)
            enddo
            call compute_Jprod_gammaprod(
     &           brapp,ketp,hel,
     &           p,q,l,t,pp,
     &           jmu1,jmu2,jmu3)
            tri_box = virt_tri_box(qsq,lsq,pq,pl,ppp,ppq,ppl,ql)
            pent    =    virt_pent(qsq,lsq,pq,pl,ppp,ppq,ppl,ql)
            Jret(compon) = tri_box + pent
         enddo
      elseif  (num_ret_current.eq.3) then
         do mu=0,3
            jmu1(mu)=jal(mu)
            jmu2(mu)=jbe(mu)
         enddo
         do compon=0,3
            do mu=0,3
               jmu3(mu) = jref(compon,mu)
            enddo
            call compute_Jprod_gammaprod(
     &           brapp,ketp,hel,
     &           p,q,l,t,pp,
     &           jmu1,jmu2,jmu3)
            tri_box = virt_tri_box(qsq,lsq,pq,pl,ppp,ppq,ppl,ql)
            pent    =    virt_pent(qsq,lsq,pq,pl,ppp,ppq,ppl,ql)
            Jret(compon) = tri_box + pent
         enddo
      else
         write(*,*) "WRONG returned current",num_ret_current
         stop
      endif
      end


c     Compute J.p_i, where J"s are the external currents, and all the 
c     <bra|...|ket> gamma prods
c     The returned values are in the common blocks
c
c  Carlo Oleari    13/6/2005
c
      subroutine compute_Jprod_gammaprod(
c     bra(p"), ket(p) and helicity 
     &     bra,ket,hel,
c     momenta
     &     p,q,l,t,pp,
c     input currents
     &     jmu1,jmu2,jmu3)
      implicit none
      double complex bra(2,-1:1),ket(2,-1:1)
      integer hel
      double precision p(0:3),pp(0:3),q(0:3),l(0:3),t(0:3)
      double complex jmu1(0:3),jmu2(0:3),jmu3(0:3)

      double complex cq(0:3),cl(0:3)
      double complex dpmu1,dpmu2,dpmu3,dqmu1,dqmu2,dqmu3
      double complex dlmu1,dlmu2,dlmu3,dppmu1,dppmu2,dppmu3
      double complex dmu1mu2,dmu1mu3,dmu2mu3
      common/dotp_pol_vec/dpmu1,dpmu2,dpmu3,dqmu1,dqmu2,dqmu3,
     &     dlmu1,dlmu2,dlmu3,dppmu1,dppmu2,dppmu3,
     &     dmu1mu2,dmu1mu3,dmu2mu3
      double complex gp_mu3_l_mu2_q_mu1,gp_l,gp_mu2,gp_mu3_l_q
      double complex gp_mu2_q_mu1,gp_mu3_mu2_q,gp_q,gp_mu3_l_mu2
      double complex gp_l_mu2_q,gp_mu3_mu2_mu1,gp_mu3_l_mu1,gp_mu3
      double complex gp_mu1,gp_l_q_mu1,gp_mu3_q_mu1,gp_l_mu2_mu1
      common/gammaprods/gp_mu3_l_mu2_q_mu1,gp_l,gp_mu2,gp_mu3_l_q,
     &     gp_mu2_q_mu1,gp_mu3_mu2_q,gp_q,gp_mu3_l_mu2,
     &     gp_l_mu2_q,gp_mu3_mu2_mu1,gp_mu3_l_mu1,gp_mu3,
     &     gp_mu1,gp_l_q_mu1,gp_mu3_q_mu1,gp_l_mu2_mu1
c      
      double complex local_mom(0:3,5), dotrc, dotcc
      external dotcc, dotrc
      integer al
      
      logical debug
c      parameter (debug=.true.)
      parameter (debug=.false.)

      dpmu1 = dotrc(p,jmu1)
      dpmu2 = dotrc(p,jmu2)
      dpmu3 = dotrc(p,jmu3)

      dqmu1 = dotrc(q,jmu1)
      dqmu2 = dotrc(q,jmu2)
      dqmu3 = dotrc(q,jmu3)

      dlmu1 = dotrc(l,jmu1)
      dlmu2 = dotrc(l,jmu2)
      dlmu3 = dotrc(l,jmu3)

      dppmu1 = dotrc(pp,jmu1)
      dppmu2 = dotrc(pp,jmu2)
      dppmu3 = dotrc(pp,jmu3)

      dmu1mu2 = dotcc(jmu1,jmu2)
      dmu1mu3 = dotcc(jmu1,jmu3)
      dmu2mu3 = dotcc(jmu2,jmu3)
      
      do al=0,3
         cl(al) = l(al)
         cq(al) = q(al)
      enddo

      do al=0,3
         local_mom(al,1) = jmu3(al)
         local_mom(al,2) = cl(al)
         local_mom(al,3) = jmu2(al)
         local_mom(al,4) = cq(al)
         local_mom(al,5) = jmu1(al)
      enddo
      call bra_slash_ket(bra,ket,hel,local_mom,5,gp_mu3_l_mu2_q_mu1)

      do al=0,3
         local_mom(al,1) = jmu3(al)
         local_mom(al,2) = cl(al)
         local_mom(al,3) = cq(al)
      enddo
      call bra_slash_ket(bra,ket,hel,local_mom,3,gp_mu3_l_q)

      do al=0,3
         local_mom(al,1) = jmu2(al)
         local_mom(al,2) = cq(al)
         local_mom(al,3) = jmu1(al)
      enddo
      call bra_slash_ket(bra,ket,hel,local_mom,3,gp_mu2_q_mu1)

      do al=0,3
         local_mom(al,1) = jmu3(al)
         local_mom(al,2) = jmu2(al)
         local_mom(al,3) = cq(al)
      enddo
      call bra_slash_ket(bra,ket,hel,local_mom,3,gp_mu3_mu2_q)

      do al=0,3
         local_mom(al,1) = jmu3(al)
         local_mom(al,2) = cl(al)
         local_mom(al,3) = jmu2(al)
      enddo
      call bra_slash_ket(bra,ket,hel,local_mom,3,gp_mu3_l_mu2)

      do al=0,3
         local_mom(al,1) = cl(al)
         local_mom(al,2) = jmu2(al)
         local_mom(al,3) = cq(al)
      enddo
      call bra_slash_ket(bra,ket,hel,local_mom,3,gp_l_mu2_q)

      do al=0,3
         local_mom(al,1) = jmu3(al)
         local_mom(al,2) = jmu2(al)
         local_mom(al,3) = jmu1(al)
      enddo
      call bra_slash_ket(bra,ket,hel,local_mom,3,gp_mu3_mu2_mu1)

      do al=0,3
         local_mom(al,1) = jmu3(al)
         local_mom(al,2) = cl(al)
         local_mom(al,3) = jmu1(al)
      enddo
      call bra_slash_ket(bra,ket,hel,local_mom,3,gp_mu3_l_mu1)

      do al=0,3
         local_mom(al,1) = cl(al)
         local_mom(al,2) = cq(al)
         local_mom(al,3) = jmu1(al)
      enddo
      call bra_slash_ket(bra,ket,hel,local_mom,3,gp_l_q_mu1)

      do al=0,3
         local_mom(al,1) = jmu3(al)
         local_mom(al,2) = cq(al)
         local_mom(al,3) = jmu1(al)
      enddo
      call bra_slash_ket(bra,ket,hel,local_mom,3,gp_mu3_q_mu1)

      do al=0,3
         local_mom(al,1) = cl(al)
         local_mom(al,2) = jmu2(al)
         local_mom(al,3) = jmu1(al)
      enddo
      call bra_slash_ket(bra,ket,hel,local_mom,3,gp_l_mu2_mu1)

      do al=0,3
         local_mom(al,1) = cl(al)
      enddo
      call bra_slash_ket(bra,ket,hel,local_mom,1,gp_l)

      do al=0,3
         local_mom(al,1) = jmu2(al)
      enddo
      call bra_slash_ket(bra,ket,hel,local_mom,1,gp_mu2)

      do al=0,3
         local_mom(al,1) = cq(al)
      enddo
      call bra_slash_ket(bra,ket,hel,local_mom,1,gp_q)

      do al=0,3
         local_mom(al,1) = jmu3(al)
      enddo
      call bra_slash_ket(bra,ket,hel,local_mom,1,gp_mu3)

      do al=0,3
         local_mom(al,1) = jmu1(al)
      enddo
      call bra_slash_ket(bra,ket,hel,local_mom,1,gp_mu1)

      if (debug) then
         write(*,*) "gamma prods"
         write(*,*) gp_l,gp_q,gp_mu2,gp_mu3,gp_mu3_l_mu2,gp_mu3_l_q,
     &        gp_mu3_mu2_q,gp_mu3_l_mu2,gp_l_mu2_q
         
         write(*,*) 
         write(*,*) gp_mu3_l_mu2_q_mu1,gp_mu2_q_mu1,gp_mu3_mu2_mu1,
     &        gp_mu3_l_mu1,gp_l_q_mu1,gp_l_mu2_mu1,gp_mu3_q_mu1,gp_mu1
      endif

      if (debug) then
         write(*,*) gp_mu3_l_mu2_q_mu1,gp_l,gp_mu2,gp_mu3_l_q,
     &        gp_mu2_q_mu1,gp_mu3_mu2_q,gp_q,gp_mu3_l_mu2,
     &        gp_l_mu2_q,gp_mu3_mu2_mu1,gp_mu3_l_mu1,gp_mu3,
     &        gp_mu1,gp_l_q_mu1,gp_mu3_q_mu1,gp_l_mu2_mu1
      endif
      end



