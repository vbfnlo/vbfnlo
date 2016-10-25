c
c---------- BCD_fill(pb1,pb2,pj1,pj2,pl1,pl2) ---------------------
c
c      subroutine BCD_fill(pb1,pb2,pj1,pj2,pl1,pl2)
c      implicit none
c      double precision pb1(0:3),pb2(0:3),pj1(0:3),pj2(0:3),
c     1                 pl1(0:3),pl2(0:3)
c
c BCD_fill determines the finite parts of all B0, C_ij, and D_ij 
c functions for the weak boson fusion process qQ --> qQV, V-->l1+l2 
c from the input momenta
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2002 October 14
c	Last modified:    2002 October 31

c  INPUT:  pb1,pb2           beam momenta
c          pj1,pj2           (anti)quark-jet momenta
c          pl1,pl2           decay momenta of the vector boson V
c
c     pb1  -------------------------- pj1
c                 $         ~ q2                    q2 = pl1+pl2 always
c                 $           ~------ pl1
c                 $ q1         \
c                 $             ----- pl2
c                 $
c     pb2  -------------------------- pj2
c
c The principal output are the D_ij(k2,q2,q1) functions for the various
c crossed processes and associated form factor functions Teps1(q1^2,t),
c Teps2(q2^2,t), Tborn(q1^2,q^2,t) as defined on p.4.8 and 4.9 of my 
c notes. The Tborn defined here combines several terms which appear 
c as coefficients of the Born amplitude M_6:
c
c  Tborn(q1^2,q2^2,t) = T_b(q1^2,t)+T_b(q2^2,t)+B0t(t)-5+pi^2/3
c
c The D_ij etc. are calculated for ID=1,...,8 cases
c
c   ID = 1    q-line, t-channel, k1 = pb1, k2 = pj1, q1 = pj2-pb2
c                     t = t12 = (q2+pj1)^2
c      = 2    q-line, u-channel, k1 = pb1, k2 = pj1, q1 = pj2-pb2
c                     t = u12 = (q1+pj1)^2    i.e. q1 <--> q2
c      = 3    qbar-line, t-channel, k1 = -pj1, k2 = -pb1, q1 = pj2-pb2
c                     t = u12 = (q1+pj1)^2
c      = 4    qbar-line, u-channel, k1 = -pj1, k2 = -pb1, q1 = pj2-pb2
c                     t = u12 = (q2+pj1)^2    i.e. q1 <--> q2
c
c   ID = 5    Q-line, t-channel, k1 = pb2, k2 = pj2, q1 = pj1-pb1
c                     t = t12 = (q2+pj2)^2
c      = 6    Q-line, u-channel, k1 = pb2, k2 = pj2, q1 = pj1-pb1
c                     t = u12 = (q1+pj2)^2    i.e. q1 <--> q2
c      = 7    Qbar-line, t-channel, k1 = -pj2, k2 = -pb2, q1 = pj1-pb1
c                     t = u12 = (q1+pj2)^2
c      = 8    Qbar-line, u-channel, k1 = -pj2, k2 = -pb2, q1 = pj1-pb1
c                     t = u12 = (q2+pj2)^2    i.e. q1 <--> q2
c
c They are stored in D0v(ID) and Dijv(i,j,ID) etc. Similarly the 
c required Teps1 are stored in Teps1(ID) etc.
c
Common block for the output
c      complex*16 D0v(8), Dijv(3,13,8), Teps1(8), Teps2(8), Tborn(8),
c     1           Tg1(8), Tg2(8)
c      double precision s(8), t(8), q1sq(8), q2sq(8)
c      double precision q1(0:3,8), q2(0:3,8), k1(0:3,8), k2(0:3,8)
c      common /bcd_qqv/ s,t,q1sq,q2sq, D0v,Dijv, Teps1,Teps2,Tborn,
c     1                                          Tg1,Tg2
c      common /bcd_qqv_aux/ q1,q2,k1,k2
c
c local variables, external functions etc.
c      double precision p1p2, p1p3, p2p3
c      complex*16 C0123(8),  C0124(8),  C0134(8),  C0234(8)
c      complex*16 Cij123(2,4,8), Cij124(2,4,8), 
c     1           Cij134(2,4,8), Cij234(2,4,8)
c      complex*16 B0a(3), B0aq1(8), B0aq2(8), B0at(8)
c      complex*16 B0t, C0t, D0tmusq
c      external B0t, C0t, D0tmusq
c      double precision s12, t12, u12, s34, t34, u34
c      double precision x1,x2,x3,x4, tmq1i, tmq2i, q1v(0:3,2), q2v(0:3)
c      double precision dotrr, psumsq
c      external dotrr, psumsq
c      double precision pi,pi2o3
c      parameter (pi=3.14159 26535 89793d0,pi2o3=pi**2/3d0)
c      integer id, mu, id0
c      data id0/0/
c      logical ldebug, ldb0
c      parameter (ldebug=.false.)
c      common /bdebug/ ldb0
common block for debugging via bcd_fill_c
c      complex*16 D0c, Dijc(3,13), Teps1c, Teps2c, Tbornc,Tg1c,Tg2c
c      double precision sc, tc, q1sqc, q2sqc
c      common /bcd_qqv_c/ sc,tc,q1sqc,q2sqc, D0c,Dijc, 
c     1                   Teps1c,Teps2c,Tbornc,Tg1c,Tg2c
c
c set the kinematic variables for the 8 cases
c      do mu = 0,3
c         q1v(mu,1) = pj2(mu)-pb2(mu)             ! the actual q1 and q2
c         q1v(mu,2) = pj1(mu)-pb1(mu)
c         q2v(mu) = pl1(mu)+pl2(mu)
c and the effective q1 and q2 for individual subsets of diagrams
c         do id = 1,8
c            id0 = 1 + id/5
c            if (mod(id,2).eq.1) then             ! id=1,3,5,7
c               q1(mu,id) = q1v(mu,id0)
c               q2(mu,id) = q2v(mu)
c            else                                 ! id=2,4,6,8
c               q2(mu,id) = q1v(mu,id0)
c               q1(mu,id) = q2v(mu)
c            endif
c         enddo
c         do id = 1,2
c            k1(mu,id) = pb1(mu)
c            k2(mu,id) = pj1(mu)
c            k1(mu,id+2) = -pj1(mu)
c            k2(mu,id+2) = -pb1(mu)
c            k1(mu,id+4) = pb2(mu)
c            k2(mu,id+4) = pj2(mu)
c            k1(mu,id+6) = -pj2(mu)
c            k2(mu,id+6) = -pb2(mu)
c         enddo
c      enddo
c      s12 = -2*dotrr(pb1,pj1)
c      t12 = psumsq(pj1,q2(0,1))
c      u12 = psumsq(pj1,q1(0,1))
c      s34 = -2*dotrr(pb2,pj2)
c      t34 = psumsq(pj2,q2(0,5))
c      u34 = psumsq(pj2,q1(0,5))
c      do id = 1,4
c         s(id) = s12
c         s(id+4) = s34
c      enddo
c      t(1) = t12
c      t(2) = u12
c      t(5) = t34
c      t(6) = u34
c      q1sq(1) = dotrr(q1(0,1),q1(0,1))
c      q2sq(1) = dotrr(q2(0,1),q2(0,1))
c      q1sq(2) = q2sq(1)
c      q2sq(2) = q1sq(1)
c      q1sq(5) = dotrr(q1(0,5),q1(0,5))
c      q2sq(5) = q2sq(1)
c      q1sq(6) = q2sq(5)
c      q2sq(6) = q1sq(5)
c      do id = 3,4
c         t(id) = t(5-id)
c         q1sq(id) = q1sq(id-2)
c         q2sq(id) = q2sq(id-2)
c         t(id+4) = t(9-id)
c         q1sq(id+4) = q1sq(id+2)
c         q2sq(id+4) = q2sq(id+2)
c      enddo
c      if (ldebug) then
c         print*," table of s,t,q1sq,q2sq "
c         print*," s12, t12, u12 = ",s12,t12,u12
c         print*," s34, t34, u34 = ",s34,t34,u34
c         do id = 1,8
c            print("(i4,4f18.4)"),id,s(id),t(id),q1sq(id),q2sq(id)
c         enddo
c         print*," q1sq,q2sq check "
c         do id = 1,8
c            print*,id,q1sq(id)/dotrr(q1(0,id),q1(0,id)),
c     1                q2sq(id)/dotrr(q2(0,id),q2(0,id))
c         enddo
c      endif
c      if (ldebug .and. pb1(2).ne.0 .and. pb2(2).ne.0) then
c      if (ldebug) then
c         do id = 1,8
c            ldb0 = id.eq.5
c            call bcd_fill_c(k1(0,id),k2(0,id),q1(0,id),q2(0,id))
c            print*," D0: ",abs(d0v(id)/d0c-1)
c            print("(a,3g10.2)")," D1j ",abs(dijv(1,1,id)/dijc(1,1)-1),
c     1                               abs(dijv(1,2,id)/dijc(1,2)-1),
c     2                               abs(dijv(1,3,id)/dijc(1,3)-1)
c            print("(a,7g10.2)")," D2j ",abs(dijv(2,1,id)/dijc(2,1)-1),
c     1                               abs(dijv(2,2,id)/dijc(2,2)-1),
c     1                               abs(dijv(2,3,id)/dijc(2,3)-1),
c     1                               abs(dijv(2,4,id)/dijc(2,4)-1),
c     1                               abs(dijv(2,5,id)/dijc(2,5)-1),
c     1                               abs(dijv(2,6,id)/dijc(2,6)-1),
c     1                               abs(dijv(2,7,id)/dijc(2,7)-1)
c            print("(a,7g10.2)")," D3j ",abs(dijv(3,1,id)/dijc(3,1)-1),
c     1                               abs(dijv(3,2,id)/dijc(3,2)-1),
c     1                               abs(dijv(3,3,id)/dijc(3,3)-1),
c     1                               abs(dijv(3,4,id)/dijc(3,4)-1),
c     1                               abs(dijv(3,5,id)/dijc(3,5)-1),
c     1                               abs(dijv(3,6,id)/dijc(3,6)-1),
c     1                               abs(dijv(3,7,id)/dijc(3,7)-1)
c            print("(a,7g10.2)"),"     ",abs(dijv(3,8,id)/dijc(3,8)-1),
c     1                       abs(dijv(3,9,id)/dijc(3,9)-1),
c     1                       abs(dijv(3,10,id)/dijc(3,10)-1),
c     1                       abs(dijv(3,11,id)/dijc(3,11)-1),
c     1                       abs(dijv(3,12,id)/dijc(3,12)-1),
c     1                       abs(dijv(3,13,id)/dijc(3,13)-1)
c         enddo
c      endif
c
c now determine the required B0, C0 and D0 functions
c D0 = D0t(s,t,q1sq,q2sq) = D0t(s,t,q2sq,q1sq) with musq = -s(id)
c      do id = 1,7,2
c         D0v(id) = D0tmusq(s(id),t(id),q1sq(id),q2sq(id),-s(id))
c      enddo
c      do id = 2,6,4
c         D0v(id) = D0v(id+1)
c         D0v(id+2) = D0v(id-1)
c      enddo
c
c C0123 = C0(0,q2sq,t), C0124 = C0(0,s,0), C0134 = C0(t,q1sq,0),
c C0234 = C0(q2sq,q1sq,s)
c      do id = 1,5,4
c         C0124(id) = C0t(0d0,s(id),0d0,-s(id))
c         C0234(id) = C0t(q2sq(id),q1sq(id),s(id),-s(id))
c      enddo
c      do id = 2,4
c         C0124(id) = C0124(1)
c         C0234(id) = C0234(1)
c         C0124(id+4) = C0124(5)
c         C0234(id+4) = C0234(5)
c      enddo
c      do id = 1,8,2
c         C0123(id) = C0t(0d0,q2sq(id),t(id),-s(id))
c         C0134(id) = C0t(t(id),q1sq(id),0d0,-s(id))
c      enddo
c      do id = 2,6,4
c         C0123(id) = C0134(id+1)
c         C0134(id) = C0123(id+1)
c         C0123(id+2) = C0134(id-1)
c         C0134(id+2) = C0123(id-1)
c      enddo
c
c      if (ldebug.and..false.) then
c         print*," check of C0 functions for ID = 1...8 "
c         print*," ID,   C0123,    C0124,    C0134,     C0234/correct-1"
c         do id = 1,8
c            x1 = abs(C0123(id)/C0t(q2sq(id),t(id),0d0,-s(id))-1)
c            x2 = abs(C0124(id)/C0t(0d0,s(id),0d0,-s(id))-1)
c            x3 = abs(C0134(id)/C0t(q1sq(id),0d0,t(id),-s(id))-1)
c            x4 = abs(C0234(id)/C0t(q2sq(id),q1sq(id),s(id),-s(id))-1)
c            print("(i4,4g12.4)"),id,x1,x2,x3,x4
c         enddo
c         do id = 1,8
c            print*," id = ",id
c            print*," s,t,u=",s(id),t(id),q1sq(id)+q2sq(id)-s(id)-t(id)
c            print*," q1sq,q2sq = ",q1sq(id),q2sq(id)
c            print*," 2k2.q1 = ",q2sq(id)-s(id)-t(id)
c            print*," 2k2.q2 = ",t(id)-q2sq(id)
c         enddo
c      endif
c      if (ldebug) pause
c
c and now the B0 functions: B0t(q1sq), B0t(q2sq), B0t(t)
c Note: B0t(s,musq=s) = 2
c      do id = 1,5,4
c         B0aq1(id) = B0t(q1sq(id),-s(id))
c         B0aq2(id) = B0t(q2sq(id),-s(id))
c      enddo
c      do id = 2,4
c         B0aq1(id) = B0aq2(id-1)
c         B0aq2(id) = B0aq1(id-1)
c      enddo
c      do id = 6,8
c         B0aq1(id) = B0aq2(id-1)
c         B0aq2(id) = B0aq1(id-1)
c      enddo
c      do id = 1,2
c         B0at(id) = B0t(t(id),-s(id))
c         B0at(id+4) = B0t(t(id+4),-s(id+4))
c      enddo
c      do id = 3,7,4
c         B0at(id) = B0at(id-1)
c         B0at(id+1) = B0at(id-2)
c      enddo
c
c      if (ldebug) then
c         print*," check of B0 functions for ID = 1...8 "
c         print*," ID,   B0(q1sq),  B0(q2sq),  B0(t):|stored/correct-1|"
c         do id = 1,8
c            x1 = abs(B0aq1(id)/B0t(q1sq(id),-s(id))-1)
c            x2 = abs(B0aq2(id)/B0t(q2sq(id),-s(id))-1)
c            x3 = abs(B0at(id)/B0t(t(id),-s(id))-1)
c            x4 = abs(2/B0t(s(id),-s(id))-1)
c            print("(i4,4g12.4)"),id,x1,x2,x3,x4
c         enddo
c      endif
c
c Now everything is ready to calculate the C_ij and D_ij functions
c      do id = 1,8
c         B0a(1) = 0
c         B0a(2) = B0aq2(id)
c         B0a(3) = B0at(id)
c         call tens3(0d0,0d0,q2sq(id),t(id),B0a,C0123(id),
c     1              Cij123(1,1,id))
c         B0a(1) = 0d0
c         B0a(2) = 2
c         B0a(3) = 0d0
c         call tens3(0d0,0d0,s(id),0d0,B0a,C0124(id),
c     1              Cij124(1,1,id))
c         B0a(1) = B0at(id)
c         B0a(2) = B0aq1(id)
c         B0a(3) = 0d0
c         call tens3(0d0,t(id),q1sq(id),0d0,B0a,C0134(id),
c     1              Cij134(1,1,id))
c         B0a(1) = B0aq2(id)
c         B0a(2) = B0aq1(id)
c         B0a(3) = 2
c         call tens3(0d0,q2sq(id),q1sq(id),s(id),B0a,C0234(id),
c     1              Cij234(1,1,id))
c         p1p2 = (t(id)-q2sq(id))/2
c         p1p3 = (q2sq(id)-s(id)-t(id))/2
c         p2p3 = (s(id)-q1sq(id)-q2sq(id))/2
c         ldb0 = ldebug .and. id.eq.5
c         call tens4( 0d0, 0d0, q2sq(id), q1sq(id), p1p2, p1p3, p2p3, 
c     &               C0123(id),  C0124(id),  C0134(id),  C0234(id),
c     &               Cij123(1,1,id), Cij124(1,1,id), 
c     &               Cij134(1,1,id), Cij234(1,1,id), 
c     &               D0v(id), Dijv(1,1,id) )
c         if(ldebug) pause
c and the coefficients of the tensors in M_6:
c         tmq1i = 1d0/(t(id)-q1sq(id))
c         Teps1(id) = ( (B0at(id)-B0aq1(id))*(2*t(id)+3*q1sq(id))*tmq1i
c     &                 +2*B0aq1(id)+1-2*q1sq(id)*C0134(id) ) * tmq1i
c         Tg1(id) = (B0aq1(id)-B0at(id))*tmq1i
c
c         tmq2i = 1d0/(t(id)-q2sq(id))
c         Teps2(id) = ( (B0at(id)-B0aq2(id))*(2*t(id)+3*q2sq(id))*tmq2i
c     &                 +2*B0aq2(id)+1-2*q2sq(id)*C0123(id) ) * tmq2i
c         Tg2(id) = (B0aq2(id)-B0at(id))*tmq2i
c         Tborn(id) = ( 2*q1sq(id)*(B0at(id)-B0aq1(id)) + t(id)*B0at(id)
c     &               -q1sq(id)*B0aq1(id) )*tmq1i - 2*q1sq(id)*C0134(id)
c     &             + ( 2*q2sq(id)*(B0at(id)-B0aq2(id)) + t(id)*B0at(id)
c     &               -q2sq(id)*B0aq2(id) )*tmq2i - 2*q2sq(id)*C0123(id)
c     &             + B0at(id) - 5d0 + pi2o3
C NOTE: the finite term (-5+pi^2/3) in Tborn is fixed by EW gauge invariance!
C Since I have factored out M_born_total*(-2/eps^2 + ... + pi^2/3 - 7),
C and since essentially unrelated couplings will multiply the two sets 
C of boxline graphs, each set individually must vanish when replacing
C eps_i ---> q_i in subroutine boxline
c      enddo
c      if (ldebug) then
c         id0 = mod(id0,8)+1
c         write(11,*) " "
c         write(11,*) " id = ",id0
c         write(11,"(a,4g22.14)") " k1 = ",(k1(mu,id0),mu=0,3)
c         write(11,"(a,4g22.14)") " k2 = ",(k2(mu,id0),mu=0,3)
c         write(11,"(a,4g22.14)") " q1 = ",(q1(mu,id0),mu=0,3)
c         write(11,"(a,4g22.14)") " q2 = ",(q2(mu,id0),mu=0,3)
c         write(11,*) " s, t       = ",s(id0),t(id0)
c         write(11,*) " q1^2, q2^2 = ",q1sq(id0),q2sq(id0)
c         write(11,*) " D0(k2,q2,q1)  = ",d0v(id0)
c         write(11,*) " C0(k2,q2)     = ",C0123(id0)
c         write(11,*) " C0(k2,q1+q2)  = ",C0124(id0)
c         write(11,*) " C0(q2,q1)     = ",C0234(id0)
c         write(11,*) " C0(k2+q2,q1)  = ",C0134(id0)
c         do id = 1,8
c            call bcd_fill_c(k1(0,id),k2(0,id),q1(0,id),q2(0,id))
c            print*," D0: ",abs(d0v(id)/d0c-1)
c            print("(a,3g10.2)")," D1j ",abs(dijv(1,1,id)/dijc(1,1)-1),
c     1                               abs(dijv(1,2,id)/dijc(1,2)-1),
c     2                               abs(dijv(1,3,id)/dijc(1,3)-1)
c            print("(a,7g10.2)")," D2j ",abs(dijv(2,1,id)/dijc(2,1)-1),
c     1                               abs(dijv(2,2,id)/dijc(2,2)-1),
c     1                               abs(dijv(2,3,id)/dijc(2,3)-1),
c     1                               abs(dijv(2,4,id)/dijc(2,4)-1),
c     1                               abs(dijv(2,5,id)/dijc(2,5)-1),
c     1                               abs(dijv(2,6,id)/dijc(2,6)-1),
c     1                               abs(dijv(2,7,id)/dijc(2,7)-1)
c            print("(a,7g10.2)")," D3j ",abs(dijv(3,1,id)/dijc(3,1)-1),
c     1                               abs(dijv(3,2,id)/dijc(3,2)-1),
c     1                               abs(dijv(3,3,id)/dijc(3,3)-1),
c     1                               abs(dijv(3,4,id)/dijc(3,4)-1),
c     1                               abs(dijv(3,5,id)/dijc(3,5)-1),
c     1                               abs(dijv(3,6,id)/dijc(3,6)-1),
c     1                               abs(dijv(3,7,id)/dijc(3,7)-1)
c            print("(a,7g10.2)"),"     ",abs(dijv(3,8,id)/dijc(3,8)-1),
c     1                       abs(dijv(3,9,id)/dijc(3,9)-1),
c     1                       abs(dijv(3,10,id)/dijc(3,10)-1),
c     1                       abs(dijv(3,11,id)/dijc(3,11)-1),
c     1                       abs(dijv(3,12,id)/dijc(3,12)-1),
c     1                       abs(dijv(3,13,id)/dijc(3,13)-1)
c         enddo
c      endif
c      return
c      end

c
c---------------- BCD_fill_diboson(k1,k2,q1,q2) -------------------------------
c
      subroutine BCD_fill_diboson(k1,k2,q1,q2)
      implicit none
      double precision k1(0:3),k2(0:3),q1(0:3),q2(0:3)
c
c This subroutine does the same as BCD_fill for a single momentum
c assignment as in the diagram
c
c     k1  -->---------->------------- k2
c                 $         $
c                 $         $
c                 $ q1      $  q2
c                 $         $ 
c
c results are stored in Dijc(*,*) etc. i.e.
c
Common block for the output
      complex*16 D0c, Dijc(3,13), Teps1c, Teps2c, Tbornc,Tg1c,Tg2c
      double precision sc, tc, q1sqc, q2sqc
      common /bcd_qqv_c/ D0c,Dijc,Teps1c,Teps2c,Tbornc,Tg1c,Tg2c, 
     1                   sc,tc,q1sqc,q2sqc
c
c local variables, external functions etc.
#include "VBFNLO/utilities/scales.inc"

      double precision p1p2, p1p3, p2p3
      complex*16 C0123,  C0124,  C0134,  C0234
      complex*16 Cij123(2,4), Cij124(2,4), 
     1           Cij134(2,4), Cij234(2,4)
      complex*16 B0a(3), B0aq1, B0aq2, B0at, B0as
      complex*16 B0t, C0t, D0tmusq
      external B0t, C0t, D0tmusq
      double precision tmq1i, tmq2i
      double precision dotrr, psumsq
      external dotrr, psumsq
      double precision pi,pi2o3
      parameter (pi=3.14159 26535 89793d0,pi2o3=pi**2/3d0)
c
      sc = -2*dotrr(k1,k2)
      tc = psumsq(k2,q2)
      q1sqc = dotrr(q1,q1)
      q2sqc = dotrr(q2,q2)
c
c now determine the required B0, C0 and D0 functions
c D0 = D0t(s,t,q1sq,q2sq) = D0t(s,t,q2sq,q1sq) with musq = -sc
      D0c = D0tmusq(sc,tc,q2sqc,q1sqc,sc)
c
c C0123 = C0(0,q2sq,t), C0124 = C0(0,s,0), C0134 = C0(t,q1sq,0),
c C0234 = C0(q2sq,q1sq,s)
      C0124 = C0t(0d0,sc,0d0,sc)
      C0234 = C0t(q2sqc,q1sqc,sc,sc)
      C0123 = C0t(0d0,q2sqc,tc,sc)
      C0134 = C0t(tc,q1sqc,0d0,sc)
c
c and now the B0 functions: B0t(q1sq), B0t(q2sq), B0t(t)
c Note: B0t(-s,musq=s) = 2
      B0aq1 = B0t(q1sqc,sc)
      B0aq2 = B0t(q2sqc,sc)
      B0at  = B0t(tc,sc)
      B0as  = B0t(sc,sc)
c
c Now everything is ready to calculate the C_ij and D_ij functions
      B0a(1) = 0
      B0a(2) = B0aq2
      B0a(3) = B0at
      call tens3(0d0,0d0,q2sqc,tc,B0a,C0123,Cij123(1,1))
      B0a(1) = 0d0
      B0a(2) = B0as
      B0a(3) = 0d0
      call tens3(0d0,0d0,sc,0d0,B0a,C0124,Cij124(1,1))
      B0a(1) = B0at
      B0a(2) = B0aq1
      B0a(3) = 0d0
      call tens3(0d0,tc,q1sqc,0d0,B0a,C0134,Cij134(1,1))
      B0a(1) = b0aq2
      B0a(2) = B0aq1
      B0a(3) = B0as
      call tens3(0d0,q2sqc,q1sqc,sc,B0a,C0234,Cij234(1,1))
      p1p2 = dotrr(k2,q2)
      p1p3 = dotrr(k2,q1)
      p2p3 = dotrr(q2,q1)
      call tens4( 0d0, 0d0, q2sqc, q1sqc, p1p2, p1p3, p2p3, 
     &            C0123,  C0124,  C0134,  C0234,
     &            Cij123(1,1), Cij124(1,1), 
     &            Cij134(1,1), Cij234(1,1), 
     &            D0c, Dijc(1,1) )
c and the coefficients of the tensors in M_6:
      tmq1i = 1d0/(tc-q1sqc)
      Teps1c = ( (B0at-B0aq1)*(2*tc+3*q1sqc)*tmq1i
     &     +2*B0aq1+1-2*q1sqc*C0134 ) * tmq1i
      Tg1c = (B0aq1-B0at)*tmq1i
c
      tmq2i = 1d0/(tc-q2sqc)
      Teps2c = ( (B0at-b0aq2)*(2*tc+3*q2sqc)*tmq2i
     &     +2*b0aq2+1-2*q2sqc*C0123 ) * tmq2i
      Tg2c = (B0aq2-B0at)*tmq2i
      Tbornc = ( 2*q1sqc*(B0at-B0aq1) + tc*B0at
     &     -q1sqc*B0aq1 )*tmq1i - 2*q1sqc*C0134
     &     + ( 2*q2sqc*(B0at-b0aq2) + tc*B0at
     &     -q2sqc*b0aq2 )*tmq2i - 2*q2sqc*C0123
     &     + B0at - 5d0 + pi2o3
      return
      end

c------------- D0(s,t,q1sq,q2sq)  massless 4-point function -----------
c
      complex*16 function D0tmusq(s,t,q1sq,q2sq,musq)
      implicit none
      double precision s, t, q1sq, q2sq, musq
      double precision sl, tl, q1sql, q2sql
      double precision p1sq, p2sq, p3sq, p1p2, p1p3, p2p3

c evaluate the finite part of the scalar 4-point function for zero mass
c propagators 
c  
c  D0 = 1/(i*pi^2) * Int d^dk [-k^2-i0]^-1 [-(k+p1)^2-i0]^-1 
c                      [-(k+p1+p2)^2-i0]^-1 [-(k+p1+p2+p3)^2-i0]^-1
c
c       = pi^-eps Gamma(1+eps) musq^-eps [IR + D0t(s,t,q1sq,q2sq,musq)]
c
c where p1^2 = 0 = (p1+p2+p3)^2 is assumed. The arguments of D0t are
c 
c   s = (p2+p3)^2, t = (p1+p2)^2, q1sq = p3^3 and q2sq = p2^2
c
c IR represents the divergent terms and is given by
c 
c  IR = 1/(s*t)*[1/eps**2 + 1/eps (log(q1sq/t)+log(q2sq/t)+log(musq/-s)) ]
c
c Alternatively, arguments as in the call of tens4 can be used, i.e.
c
c  D0t(s,t,q1sq,q2sq) = D0t4(p1sq, p2sq, p3sq, p1p2, p1p3, p2p3) 
c
c This code agrees with the one written by Carlo. Checked 10/15/2002
c       Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2002 October 5
c	Last modified:    2002 October 14
c  
      complex*16 ipi, ieps, z1, z2, zlog, vli2, D0t4musq
      parameter (ipi=(0d0,3.14159 26535 89793d0),ieps=(0d0,1d-16))
      external vli2
      double precision pi,pi2o2,pi2o6
      parameter (pi=3.14159 26535 89793d0,pi2o2=pi**2/2d0)
      parameter (pi2o6=pi**2/6d0)
      logical ldebug
      parameter (ldebug=.false.)
      double precision rlog, rz1, rz2, rz3
c
      sl = s
      tl = t
      q1sql = q1sq
      q2sql = q2sq
      goto 10
c
      entry D0t4musq(p1sq, p2sq, p3sq, p1p2, p1p3, p2p3, musq)
      if (p1sq.ne.0d0 .or. 
     1    abs(2*(p1p2+p2p3+p1p3)/(p2sq+p3sq)+1).gt.1d-12) then
         print*," Dont use D0t4 for for less than 2 massless momenta "
         print*," (p1p2+p2p3+p1p3)/(p2sq+p3sq) = ",
     1           (p1p2+p2p3+p1p3)/(p2sq+p3sq)
         stop
      endif
      sl = p2sq+2*p2p3+p3sq
      tl = p2sq+2*p1p2
      q1sql = p3sq
      q2sql = p2sq
 10   continue
      if (q1sql.eq.0d0 .or. q2sql.eq.0d0) then
         print*," WARNING: zero Q^2 in D0t: ",q1sql,q2sql
         print*," Set D0t = 0 "
         D0tmusq = 0
         return 
      endif
      if (musq.le.0d0) then
         print*," WARNING: non-positive musq in D0t: ",musq
         print*," Set D0t = 0 "
         D0tmusq = 0
         return 
      endif
         
c
      rz1 = tl/q1sql
      rz2 = tl/q2sql
      rz3 = -sl/musq
      rlog = -log( abs(rz1*rz2*rz3) )
      zlog = rlog
      if (rz1.gt.0) then
         z1 = 1-rz1
      elseif (tl.lt.0d0) then        !   t<0, q1sq>0
         z1 = 1-rz1-ieps
         zlog = zlog - ipi
      else                           !   t>0, q1sq<0
         z1 = 1-rz1+ieps
         zlog = zlog + ipi
      endif
      if (rz2.gt.0) then
         z2 = 1-rz2
      elseif (tl.lt.0d0) then        !   t<0, q2sq>0
         z2 = 1-rz2-ieps
         zlog = zlog - ipi
      else                           !   t>0, q2sq<0
         z2 = 1-rz2+ieps
         zlog = zlog + ipi
      endif
      if (rz3.lt.0d0) then           !   s>0
         zlog = zlog + ipi
      endif
      D0tmusq = (0.5d0*zlog**2 + 2*(vli2(z1)+vli2(z2)) - pi2o6)/(sl*tl)
c
      if (ldebug) then
         print("(a,2g14.4,a,g14.4)"), 
     1         " t/qsq ratios: ",rz1, rz2," t = ",tl
         print("(a,f8.4,a,f8.4)"),
     1         " zlog = ",dreal(zlog),"+ i pi * ",dimag(zlog)/pi
         z1 = 0.5d0*zlog**2
         rz1 = (dreal(z1)-0.5d0*rlog**2)/pi**2
         if (rlog.ne.0d0) then
            rz2 = dimag(z1)/pi/rlog
         else
            rz2 = 0d0
         endif
         print("(2a,f8.4,a,f8.4,a)")," extra terms from log: ",
     1   " pi**2*(",rz1,")+ i pi rlog *(",rz2,")" 
         print*," D0t = ",d0tmusq
      endif
      return 
      end


