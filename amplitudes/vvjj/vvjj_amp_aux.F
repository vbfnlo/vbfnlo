c
c---------- BCD_fill(pb1,pb2,pj1,pj2,pl1,pl2) ---------------------
c
      subroutine BCD_fill(pb1,pb2,pj1,pj2,pl1,pl2)
      implicit none
      double precision pb1(0:3),pb2(0:3),pj1(0:3),pj2(0:3),
     1                 pl1(0:3),pl2(0:3)
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
      complex*16 D0v(8), Dijv(3,13,8), Teps1(8), Teps2(8), Tborn(8),
     1           Tg1(8), Tg2(8)
      double precision s(8), t(8), q1sq(8), q2sq(8)
      double precision q1(0:3,8), q2(0:3,8), k1(0:3,8), k2(0:3,8)
*      common /bcd_qqv/ s,t,q1sq,q2sq, D0v,Dijv, Teps1,Teps2,Tborn,
*     1                                          Tg1,Tg2
      common /bcd_qqv/ D0v,Dijv, Teps1,Teps2,Tborn,Tg1,Tg2,
     1                 s,t,q1sq,q2sq
      common /bcd_qqv_aux/ q1,q2,k1,k2
c
c local variables, external functions etc.
      double precision p1p2, p1p3, p2p3
      complex*16 C0123(8),  C0124(8),  C0134(8),  C0234(8)
      complex*16 Cij123(2,4,8), Cij124(2,4,8), 
     1           Cij134(2,4,8), Cij234(2,4,8)
      complex*16 B0a(3), B0aq1(8), B0aq2(8), B0at(8)
      complex*16 B0t, C0t, D0t
      external B0t, C0t, D0t
      double precision s12, t12, u12, s34, t34, u34
      double precision x1,x2,x3,x4, tmq1i, tmq2i, q1v(0:3,2), q2v(0:3)
      double precision dotrr, psumsq
      external dotrr, psumsq
      double precision pi,pi2o3
      parameter (pi=3.14159 26535 89793d0,pi2o3=pi**2/3d0)
      integer id, mu, id0
      data id0/0/
      logical ldebug, ldb0
      parameter (ldebug=.false.)
      common /bdebug/ ldb0
common block for debugging via bcd_fill_c
      complex*16 D0c, Dijc(3,13), Teps1c, Teps2c, Tbornc,Tg1c,Tg2c
      double precision sc, tc, q1sqc, q2sqc
      common /bcd_qqv_c/ D0c,Dijc,Teps1c,Teps2c,Tbornc,Tg1c,Tg2c, 
     1                   sc,tc,q1sqc,q2sqc
c
c set the kinematic variables for the 8 cases
      do mu = 0,3
         q1v(mu,1) = pj2(mu)-pb2(mu)             ! the actual q1 and q2
         q1v(mu,2) = pj1(mu)-pb1(mu)
         q2v(mu) = pl1(mu)+pl2(mu)
c and the effective q1 and q2 for individual subsets of diagrams
         do id = 1,8
            id0 = 1 + id/5
            if (mod(id,2).eq.1) then             ! id=1,3,5,7
               q1(mu,id) = q1v(mu,id0)
               q2(mu,id) = q2v(mu)
            else                                 ! id=2,4,6,8
               q2(mu,id) = q1v(mu,id0)
               q1(mu,id) = q2v(mu)
            endif
         enddo
         do id = 1,2
            k1(mu,id) = pb1(mu)
            k2(mu,id) = pj1(mu)
            k1(mu,id+2) = -pj1(mu)
            k2(mu,id+2) = -pb1(mu)
            k1(mu,id+4) = pb2(mu)
            k2(mu,id+4) = pj2(mu)
            k1(mu,id+6) = -pj2(mu)
            k2(mu,id+6) = -pb2(mu)
         enddo
      enddo
      s12 = -2*dotrr(pb1,pj1)
      t12 = psumsq(pj1,q2(0,1))
      u12 = psumsq(pj1,q1(0,1))
      s34 = -2*dotrr(pb2,pj2)
      t34 = psumsq(pj2,q2(0,5))
      u34 = psumsq(pj2,q1(0,5))
      do id = 1,4
         s(id) = s12
         s(id+4) = s34
      enddo
      t(1) = t12
      t(2) = u12
      t(5) = t34
      t(6) = u34
      q1sq(1) = dotrr(q1(0,1),q1(0,1))
      q2sq(1) = dotrr(q2(0,1),q2(0,1))
      q1sq(2) = q2sq(1)
      q2sq(2) = q1sq(1)
      q1sq(5) = dotrr(q1(0,5),q1(0,5))
      q2sq(5) = q2sq(1)
      q1sq(6) = q2sq(5)
      q2sq(6) = q1sq(5)
      do id = 3,4
         t(id) = t(5-id)
         q1sq(id) = q1sq(id-2)
         q2sq(id) = q2sq(id-2)
         t(id+4) = t(9-id)
         q1sq(id+4) = q1sq(id+2)
         q2sq(id+4) = q2sq(id+2)
      enddo
      if (ldebug) then
         print*," table of s,t,q1sq,q2sq "
         print*," s12, t12, u12 = ",s12,t12,u12
         print*," s34, t34, u34 = ",s34,t34,u34
         do id = 1,8
            print("(i4,4f18.4)"),id,s(id),t(id),q1sq(id),q2sq(id)
         enddo
         print*," q1sq,q2sq check "
         do id = 1,8
            print*,id,q1sq(id)/dotrr(q1(0,id),q1(0,id)),
     1                q2sq(id)/dotrr(q2(0,id),q2(0,id))
         enddo
      endif
c      if (ldebug .and. pb1(2).ne.0 .and. pb2(2).ne.0) then
      if (ldebug) then
         do id = 1,8
            ldb0 = id.eq.5
            call bcd_fill_c(k1(0,id),k2(0,id),q1(0,id),q2(0,id))
            print*," D0: ",abs(d0v(id)/d0c-1)
            print("(a,3g10.2)")," D1j ",abs(dijv(1,1,id)/dijc(1,1)-1),
     1                               abs(dijv(1,2,id)/dijc(1,2)-1),
     2                               abs(dijv(1,3,id)/dijc(1,3)-1)
            print("(a,7g10.2)")," D2j ",abs(dijv(2,1,id)/dijc(2,1)-1),
     1                               abs(dijv(2,2,id)/dijc(2,2)-1),
     1                               abs(dijv(2,3,id)/dijc(2,3)-1),
     1                               abs(dijv(2,4,id)/dijc(2,4)-1),
     1                               abs(dijv(2,5,id)/dijc(2,5)-1),
     1                               abs(dijv(2,6,id)/dijc(2,6)-1),
     1                               abs(dijv(2,7,id)/dijc(2,7)-1)
            print("(a,7g10.2)")," D3j ",abs(dijv(3,1,id)/dijc(3,1)-1),
     1                               abs(dijv(3,2,id)/dijc(3,2)-1),
     1                               abs(dijv(3,3,id)/dijc(3,3)-1),
     1                               abs(dijv(3,4,id)/dijc(3,4)-1),
     1                               abs(dijv(3,5,id)/dijc(3,5)-1),
     1                               abs(dijv(3,6,id)/dijc(3,6)-1),
     1                               abs(dijv(3,7,id)/dijc(3,7)-1)
            print("(a,7g10.2)"),"     ",abs(dijv(3,8,id)/dijc(3,8)-1),
     1                       abs(dijv(3,9,id)/dijc(3,9)-1),
     1                       abs(dijv(3,10,id)/dijc(3,10)-1),
     1                       abs(dijv(3,11,id)/dijc(3,11)-1),
     1                       abs(dijv(3,12,id)/dijc(3,12)-1),
     1                       abs(dijv(3,13,id)/dijc(3,13)-1)
         enddo
      endif
c
c now determine the required B0, C0 and D0 functions
c D0 = D0t(s,t,q1sq,q2sq) = D0t(s,t,q2sq,q1sq)
      do id = 1,7,2
         D0v(id) = D0t(s(id),t(id),q1sq(id),q2sq(id))
      enddo
      do id = 2,6,4
         D0v(id) = D0v(id+1)
         D0v(id+2) = D0v(id-1)
      enddo
c
c C0123 = C0(0,q2sq,t), C0124 = C0(0,s,0), C0134 = C0(t,q1sq,0),
c C0234 = C0(q2sq,q1sq,s)
      do id = 1,5,4
         C0124(id) = C0t(0d0,s(id),0d0,-s(id))
         C0234(id) = C0t(q2sq(id),q1sq(id),s(id),-s(id))
      enddo
      do id = 2,4
         C0124(id) = C0124(1)
         C0234(id) = C0234(1)
         C0124(id+4) = C0124(5)
         C0234(id+4) = C0234(5)
      enddo
      do id = 1,8,2
         C0123(id) = C0t(0d0,q2sq(id),t(id),-s(id))
         C0134(id) = C0t(t(id),q1sq(id),0d0,-s(id))
      enddo
      do id = 2,6,4
         C0123(id) = C0134(id+1)
         C0134(id) = C0123(id+1)
         C0123(id+2) = C0134(id-1)
         C0134(id+2) = C0123(id-1)
      enddo
c
      if (ldebug.and..false.) then
         print*," check of C0 functions for ID = 1...8 "
         print*," ID,   C0123,    C0124,    C0134,     C0234/correct-1"
         do id = 1,8
            x1 = abs(C0123(id)/C0t(q2sq(id),t(id),0d0,-s(id))-1)
            x2 = abs(C0124(id)/C0t(0d0,s(id),0d0,-s(id))-1)
            x3 = abs(C0134(id)/C0t(q1sq(id),0d0,t(id),-s(id))-1)
            x4 = abs(C0234(id)/C0t(q2sq(id),q1sq(id),s(id),-s(id))-1)
            print("(i4,4g12.4)"),id,x1,x2,x3,x4
         enddo
         do id = 1,8
            print*," id = ",id
            print*," s,t,u=",s(id),t(id),q1sq(id)+q2sq(id)-s(id)-t(id)
            print*," q1sq,q2sq = ",q1sq(id),q2sq(id)
            print*," 2k2.q1 = ",q2sq(id)-s(id)-t(id)
            print*," 2k2.q2 = ",t(id)-q2sq(id)
         enddo
      endif
      if (ldebug) read(*,*)
c
c and now the B0 functions: B0t(q1sq), B0t(q2sq), B0t(t)
c Note: B0t(s,musq=s) = 2
      do id = 1,5,4
         B0aq1(id) = B0t(q1sq(id),-s(id))
         B0aq2(id) = B0t(q2sq(id),-s(id))
      enddo
      do id = 2,4
         B0aq1(id) = B0aq2(id-1)
         B0aq2(id) = B0aq1(id-1)
      enddo
      do id = 6,8
         B0aq1(id) = B0aq2(id-1)
         B0aq2(id) = B0aq1(id-1)
      enddo
      do id = 1,2
         B0at(id) = B0t(t(id),-s(id))
         B0at(id+4) = B0t(t(id+4),-s(id+4))
      enddo
      do id = 3,7,4
         B0at(id) = B0at(id-1)
         B0at(id+1) = B0at(id-2)
      enddo
c
      if (ldebug) then
         print*," check of B0 functions for ID = 1...8 "
         print*," ID,   B0(q1sq),  B0(q2sq),  B0(t):|stored/correct-1|"
         do id = 1,8
            x1 = abs(B0aq1(id)/B0t(q1sq(id),-s(id))-1)
            x2 = abs(B0aq2(id)/B0t(q2sq(id),-s(id))-1)
            x3 = abs(B0at(id)/B0t(t(id),-s(id))-1)
            x4 = abs(2/B0t(s(id),-s(id))-1)
            print("(i4,4g12.4)"),id,x1,x2,x3,x4
         enddo
      endif
c
c Now everything is ready to calculate the C_ij and D_ij functions
      do id = 1,8
         B0a(1) = 0
         B0a(2) = B0aq2(id)
         B0a(3) = B0at(id)
         call tens3(0d0,0d0,q2sq(id),t(id),B0a,C0123(id),
     1              Cij123(1,1,id))
         B0a(1) = 0d0
         B0a(2) = 2
         B0a(3) = 0d0
         call tens3(0d0,0d0,s(id),0d0,B0a,C0124(id),
     1              Cij124(1,1,id))
         B0a(1) = B0at(id)
         B0a(2) = B0aq1(id)
         B0a(3) = 0d0
         call tens3(0d0,t(id),q1sq(id),0d0,B0a,C0134(id),
     1              Cij134(1,1,id))
         B0a(1) = B0aq2(id)
         B0a(2) = B0aq1(id)
         B0a(3) = 2
         call tens3(0d0,q2sq(id),q1sq(id),s(id),B0a,C0234(id),
     1              Cij234(1,1,id))
         p1p2 = (t(id)-q2sq(id))/2
         p1p3 = (q2sq(id)-s(id)-t(id))/2
         p2p3 = (s(id)-q1sq(id)-q2sq(id))/2
         ldb0 = ldebug .and. id.eq.5
         call tens4( 0d0, 0d0, q2sq(id), q1sq(id), p1p2, p1p3, p2p3, 
     &               C0123(id),  C0124(id),  C0134(id),  C0234(id),
     &               Cij123(1,1,id), Cij124(1,1,id), 
     &               Cij134(1,1,id), Cij234(1,1,id), 
     &               D0v(id), Dijv(1,1,id) )
         if(ldebug) read(*,*)
c and the coefficients of the tensors in M_6:
         tmq1i = 1d0/(t(id)-q1sq(id))
         Teps1(id) = ( (B0at(id)-B0aq1(id))*(2*t(id)+3*q1sq(id))*tmq1i
     &                 +2*B0aq1(id)+1-2*q1sq(id)*C0134(id) ) * tmq1i
         Tg1(id) = (B0aq1(id)-B0at(id))*tmq1i
c
         tmq2i = 1d0/(t(id)-q2sq(id))
         Teps2(id) = ( (B0at(id)-B0aq2(id))*(2*t(id)+3*q2sq(id))*tmq2i
     &                 +2*B0aq2(id)+1-2*q2sq(id)*C0123(id) ) * tmq2i
         Tg2(id) = (B0aq2(id)-B0at(id))*tmq2i
         Tborn(id) = ( 2*q1sq(id)*(B0at(id)-B0aq1(id)) + t(id)*B0at(id)
     &               -q1sq(id)*B0aq1(id) )*tmq1i - 2*q1sq(id)*C0134(id)
     &             + ( 2*q2sq(id)*(B0at(id)-B0aq2(id)) + t(id)*B0at(id)
     &               -q2sq(id)*B0aq2(id) )*tmq2i - 2*q2sq(id)*C0123(id)
     &             + B0at(id) - 5d0 + pi2o3
C NOTE: the finite term (-5+pi^2/3) in Tborn is fixed by EW gauge invariance!
C Since I have factored out M_born_total*(-2/eps^2 + ... + pi^2/3 - 7),
C and since essentially unrelated couplings will multiply the two sets 
C of boxline graphs, each set individually must vanish when replacing
C eps_i ---> q_i in subroutine boxline
      enddo
      if (ldebug) then
         id0 = mod(id0,8)+1
         write(11,*) " "
         write(11,*) " id = ",id0
         write(11,"(a,4g22.14)") " k1 = ",(k1(mu,id0),mu=0,3)
         write(11,"(a,4g22.14)") " k2 = ",(k2(mu,id0),mu=0,3)
         write(11,"(a,4g22.14)") " q1 = ",(q1(mu,id0),mu=0,3)
         write(11,"(a,4g22.14)") " q2 = ",(q2(mu,id0),mu=0,3)
         write(11,*) " s, t       = ",s(id0),t(id0)
         write(11,*) " q1^2, q2^2 = ",q1sq(id0),q2sq(id0)
         write(11,*) " D0(k2,q2,q1)  = ",d0v(id0)
         write(11,*) " C0(k2,q2)     = ",C0123(id0)
         write(11,*) " C0(k2,q1+q2)  = ",C0124(id0)
         write(11,*) " C0(q2,q1)     = ",C0234(id0)
         write(11,*) " C0(k2+q2,q1)  = ",C0134(id0)
         do id = 1,8
            call bcd_fill_c(k1(0,id),k2(0,id),q1(0,id),q2(0,id))
            print*," D0: ",abs(d0v(id)/d0c-1)
            print("(a,3g10.2)")," D1j ",abs(dijv(1,1,id)/dijc(1,1)-1),
     1                               abs(dijv(1,2,id)/dijc(1,2)-1),
     2                               abs(dijv(1,3,id)/dijc(1,3)-1)
            print("(a,7g10.2)")," D2j ",abs(dijv(2,1,id)/dijc(2,1)-1),
     1                               abs(dijv(2,2,id)/dijc(2,2)-1),
     1                               abs(dijv(2,3,id)/dijc(2,3)-1),
     1                               abs(dijv(2,4,id)/dijc(2,4)-1),
     1                               abs(dijv(2,5,id)/dijc(2,5)-1),
     1                               abs(dijv(2,6,id)/dijc(2,6)-1),
     1                               abs(dijv(2,7,id)/dijc(2,7)-1)
            print("(a,7g10.2)")," D3j ",abs(dijv(3,1,id)/dijc(3,1)-1),
     1                               abs(dijv(3,2,id)/dijc(3,2)-1),
     1                               abs(dijv(3,3,id)/dijc(3,3)-1),
     1                               abs(dijv(3,4,id)/dijc(3,4)-1),
     1                               abs(dijv(3,5,id)/dijc(3,5)-1),
     1                               abs(dijv(3,6,id)/dijc(3,6)-1),
     1                               abs(dijv(3,7,id)/dijc(3,7)-1)
            print("(a,7g10.2)"),"     ",abs(dijv(3,8,id)/dijc(3,8)-1),
     1                       abs(dijv(3,9,id)/dijc(3,9)-1),
     1                       abs(dijv(3,10,id)/dijc(3,10)-1),
     1                       abs(dijv(3,11,id)/dijc(3,11)-1),
     1                       abs(dijv(3,12,id)/dijc(3,12)-1),
     1                       abs(dijv(3,13,id)/dijc(3,13)-1)
         enddo
      endif
      return
      end

c
c---------------- BCD_fill_c(k1,k2,q1,q2) -------------------------------
c
      subroutine BCD_fill_c(k1,k2,q1,q2)
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
      double precision p1p2, p1p3, p2p3
      complex*16 C0123,  C0124,  C0134,  C0234
      complex*16 Cij123(2,4), Cij124(2,4), 
     1           Cij134(2,4), Cij234(2,4)
      complex*16 B0a(3), B0aq1, B0aq2, B0at, B0as
      complex*16 B0t, C0t, D0t
      external B0t, C0t, D0t
      double precision tmq1i, tmq2i
      double precision dotrr, psumsq
      external dotrr, psumsq
      double precision pi,pi2o3
      parameter (pi=3.14159 26535 89793d0,pi2o3=pi**2/3d0)

      double complex D0fin,C0fin,B0fin
      external D0fin,C0fin,B0fin
c
      sc = -2*dotrr(k1,k2)
      tc = psumsq(k2,q2)
      q1sqc = dotrr(q1,q1)
      q2sqc = dotrr(q2,q2)

c
c now determine the required B0, C0 and D0 functions
c D0 = D0t(s,t,q1sq,q2sq) = D0t(s,t,q2sq,q1sq)
c      D0c = D0t(sc,tc,q2sqc,q1sqc)
c      Print*, "D0c", D0c
      D0c = D0fin(tc,sc,0d0,q2sqc,q1sqc,0d0,-sc)
c     print*, "D0_general", D0fin(tc,sc,0d0,q2sqc,q1sqc,0d0,-sc)/D0c
c     
c C0123 = C0(0,q2sq,t), C0124 = C0(0,s,0), C0134 = C0(t,q1sq,0),
c C0234 = C0(q2sq,q1sq,s)
c      C0124 =    C0t(0d0,sc,0d0,-sc)
      C0124 = -C0fin(0d0,sc,0d0,-sc)
c      print*, "C0124", C0fin(0d0,sc,0d0,-sc)/C0124
c      C0234 =    C0t(q2sqc,q1sqc,sc,-sc)
      C0234=  -C0fin(q2sqc,q1sqc,sc,-sc)
c      print*, "C0234", C0fin(q2sqc,q1sqc,sc,-sc)/C0234
c      C0123 =   C0t(0d0,q2sqc,tc,-sc)
      c0123 =-C0fin(0d0,q2sqc,tc,-sc)
c      print*, "C0123", C0fin(0d0,q2sqc,tc,-sc)/C0123
c      C0134 =    C0t(tc,q1sqc,0d0,-sc)
      C0134 = -C0fin(tc,q1sqc,0d0,-sc)
c      print*, "C0134", C0fin(tc,q1sqc,0d0,-sc)/C0134
c
c and now the B0 functions: B0t(q1sq), B0t(q2sq), B0t(t)
c Note: B0t(s,musq=s) = 2
c      B0aq1 =   B0t(q1sqc,-sc)
      B0aq1 = B0fin(q1sqc,-sc)
c      print*, "B0aq1 ",B0fin(q1sqc,-sc)/B0aq1
c      B0aq2 =   B0t(q2sqc,-sc)
      B0aq2 = B0fin(q2sqc,-sc)
c      print*, "B0aq2 ",B0fin(q2sqc,-sc)/B0aq2
c      B0at  =   B0t(tc,-sc)
      B0at  = B0fin(tc,-sc)
c      print*, "B0at ",B0fin(tc,-sc)/B0at
c      B0as  =   B0t(sc,-sc)
      B0as  = B0fin(sc,-sc)
c      print*, "B0at ",B0fin(tc,-sc)/B0at
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
c----  B0t(qsq): regularized 2-point function for massless loop  -----
c
c      complex*16 function B0t(qsq,musq)
c      implicit none
c      double precision qsq, musq

c evaluate the finite piece of the massless B0 function
c  
c             d^d k 
c    B0 = Int------- [-k^2+idelta]^-1 [-(k-q)^2+idelta]^-1 
c             i pi^2
c
c       = pi^-eps Gamma(1+eps) musq^-eps [1/eps + B0t(qsq,musq)]
c
c i.e. B0t(qsq,musq) = 2 - log( (-qsq-idelta)/musq )
c
c subroutine requires positive musq
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2002 September 13
c	Last modified:    2002 September 13
c  
c      complex*16 ipi
c      parameter (ipi=(0d0,3.14159 26535 89793d0))
c      logical ldebug
c      parameter (ldebug=.true.)
c
c      if (ldebug) then
c         if (musq.le.0d0) then
c            print*," Unacceptable mu^2 in B0t(q^2,mu^2): mu^2 = ",
c     1             musq," ==> set B0t = 0 "
c            b0t=0
c            return
c         endif
c      endif
c
c      if (qsq.lt.0d0) then
c         b0t = 2 - log(-qsq/musq)
c      else
c         b0t = 2 - log(qsq/musq) + ipi
c      endif
c      end
c
c------------- C0(q1^2,q2^2,P^2)  massless 3-point fuction ------------
c
c      complex*16 function C0t(q1sq,q2sq,Psq,musq)
c      implicit none
c      double precision q1sq, q2sq, psq, musq

c evaluate the finite part of the scalar 3-point function for zero mass
c propagators 
c  
c  C0 = 1/(i*pi^2) * Int d^dk [-k^2-i0]^-1 [-(k+q1)^2-i0]^-1 [-(k-q2)^2-i0]^-1
c
c       = pi^-eps Gamma(1+eps) musq^-eps [IR + C0t(q1sq,q2sq,psq)]
c
c where IR represents the potentially divergent terms when one or more
c of the arguments vanish. P = q1+q2. The divergent terms are
c 
c (1) IR = 0              all 3 of q1sq,q2sq,psq nonzero
c 
c (2) IR = 1/(t-qsq)*log(t/qsq)*1/eps    for 1 vanishing argument, the
c                                        other two being t and qsq
c
c (3) IR = -1/qsq*(1/eps**2 - 1/eps log(qsq/musq))  when exactly one 
c                                                   argument, qsq,
c                                                   is nonzero
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2002 September 13
c	Last modified:    2002 October 25
c  
c      complex*16 ipi, z, res, vli2
c      parameter (ipi=(0d0,3.14159 26535 89793d0))
c      external vli2
c      double precision pi,pi2o3,pi2o6,pi2o2
c      parameter (pi=3.14159 26535 89793d0,pi2o3=pi**2/3d0)
c      parameter (pi2o6=pi**2/6d0,pi2o2=pi**2/2d0)
c      logical ldebug, ltaylor, first_call
c      parameter (ldebug=.false.)
c      double precision x, flambda, tsqlami, qsq(3), 
c$$$     1                 sum, rz, logz, ln1, ln2,
c$$$     2                 musq_t, asum, fac, resr, ser, fb(8), xn
c$$$      integer izero, inz, i, j
c$$$      data first_call /.true./
c$$$      save first_call, fb
c$$$c
c$$$      if (ldebug) then
c$$$         if (musq.le.0d0) then
c$$$            print*," Unacceptable mu^2 in C0t(q1sq,q2sq,psq,mu^2): ",
c$$$     1         "mu^2 = ",musq," ==> set C0t = 0 "
c$$$            C0t=0
c$$$            return
c$$$         endif
c$$$      endif
c$$$c
c$$$      inz = 0
c$$$      izero = 0
c$$$      do i=1,3
c$$$         qsq(i)=0
c$$$      enddo
c$$$      if (q1sq.eq.0d0) then
c$$$         izero = izero+1
c$$$      else
c$$$         inz = inz+1
c$$$         qsq(inz)=q1sq
c$$$      endif
c$$$      if (q2sq.eq.0d0) then
c$$$         izero = izero+1
c$$$      else
c$$$         inz = inz+1
c$$$         qsq(inz)=q2sq
c$$$      endif
c$$$      if (psq.eq.0d0) then
c$$$         izero = izero+1
c$$$      else
c$$$         inz = inz+1
c$$$         qsq(inz)=psq
c$$$      endif
c$$$      if (ldebug) then
c$$$         print("(a,3g12.4,a,g12.4)"), "qsq_i = ",q1sq,q2sq,psq,
c$$$     1                                " mu^2 = ",musq
c$$$         print*," number of q_i^2 = 0  ",izero
c$$$         print*," number of q_i^2 ne 0 ",inz," Sum = ",inz+izero
c$$$         if ((inz+izero).ne.3) print*," WARNING: inz+izero .ne. 3 "
c$$$      endif
c$$$c
c$$$      if (izero.eq.0) then   ! all qsq(i) nonzero, this is full C0
c$$$         flambda = (qsq(1)**2+qsq(2)**2+qsq(3)**2) - 
c$$$     1           2*(qsq(1)*qsq(2)+qsq(1)*qsq(3)+qsq(2)*qsq(3))
c$$$         sum = qsq(1)+qsq(2)+qsq(3)
c$$$         asum = abs(qsq(1))+abs(qsq(2))+abs(qsq(3))
c$$$         ltaylor = abs(abs(sum)/asum-1).lt.1d-13      ! same sign for all qsq
c$$$         asum=min(abs(sum-2*qsq(1)),abs(sum-2*qsq(2)),
c$$$     1            abs(sum-2*qsq(3)))
c$$$         ltaylor = ltaylor .and. abs(flambda).lt.(0.01*asum**2)
c$$$c
c$$$         res = (0d0,0d0)
c$$$         if (ldebug) then
c$$$            print*," lambda = ",flambda/4," x_taylor = ",
c$$$     1                          flambda/asum**2
c$$$            ltaylor=.false.
c$$$         endif
c$$$         if (flambda.gt.0d0 .and. .not.ltaylor) then
c$$$            tsqlami = 1/sqrt(flambda)
c$$$            do i=1,3
c$$$               x = (sum-2*qsq(i))*tsqlami
c$$$               if (x.gt.1d0) then
c$$$                  rz = (x-1)/(x+1)
c$$$                  z = rz
c$$$                  logz = log(rz)
c$$$                  res = res + 2*Vli2(z)+0.5d0*logz**2-pi2o3 
c$$$                  if (qsq(i).gt.0d0) then
c$$$                     res = res + ipi*logz
c$$$                  else
c$$$                     res = res - ipi*logz
c$$$                  endif
c$$$               elseif (x.gt.0d0) then
c$$$                  rz = (x-1)/(x+1)
c$$$                  z = rz
c$$$                  logz = log(-rz)
c$$$                  res = res + 2*Vli2(z)+0.5d0*logz**2+pi2o6
c$$$               elseif (x.gt.-1d0) then
c$$$                  rz = (x+1)/(x-1)
c$$$                  z = rz
c$$$                  logz = -log(-rz)
c$$$                  res = res - 2*Vli2(z)-0.5d0*logz**2-pi2o6
c$$$               else 
c$$$                  rz = (x+1)/(x-1)
c$$$                  z = rz
c$$$                  logz = -log(rz)
c$$$                  res = res - 2*Vli2(z)-0.5d0*logz**2+pi2o3 
c$$$                  if (qsq(i).gt.0d0) then
c$$$                     res = res + ipi*logz
c$$$                  else
c$$$                     res = res - ipi*logz
c$$$                  endif
c$$$               endif
c$$$               if (ldebug) print*," i,x, z = ",i,x,rz
c$$$            enddo
c$$$            C0t = res*tsqlami
c$$$         elseif (flambda.lt.0d0 .and. .not.ltaylor) then  ! the lambda < 0 case
c$$$            tsqlami = 1/sqrt(-flambda)         ! + imag part of 2 sqrt(lambda)
c$$$            do i=1,3
c$$$               x = (sum-2*qsq(i))*tsqlami
c$$$               z = dcmplx(x,-1d0)/dcmplx(x,1d0)
c$$$               if (ldebug) then  
c$$$                  print*," i,x, z = ",i,x,z
c$$$               endif
c$$$               res = res + dimag(Vli2(z))
c$$$            enddo
c$$$            C0t = res*2*tsqlami
c$$$         endif
c$$$         if (ldebug) then
c$$$            print*," C0t via Vli2 calls: ",C0t
c$$$            ltaylor = abs(flambda).lt.(0.5*asum**2)
c$$$         endif
c$$$         if (ltaylor) then
c$$$            if (first_call) then
c$$$               fb(1) = 2d0/9d0
c$$$               fb(2) = 16d0/75d0
c$$$               fb(3) = 142d0/735d0
c$$$               fb(4) = 496d0/2835d0
c$$$               fb(5) = 6086d0/38115d0
c$$$               fb(6) = 86048d0/585585d0
c$$$               fb(7) = 92054d0/675675d0
c$$$               fb(8) = 1655008d0/13018005d0 
c$$$c               fb(9) = 32976682d0/276441165d0
c$$$               first_call = .false.
c$$$            endif
c$$$            musq_t = abs(sum/3)
c$$$            resr = 0
c$$$            do i = 1,3
c$$$               fac = 1/(sum-2*qsq(i))
c$$$               x = flambda*fac**2
c$$$c               if (ldebug) print*," i = ",i," x = ",x
c$$$               ln1 = log(abs(qsq(i))/musq_t)
c$$$               ser = ln1
c$$$               j = 1
c$$$               xn = x
c$$$c               if (ldebug) print*," log = ",ln1
c$$$               do while (abs(xn).gt.1d-15 .and. j.lt.9)
c$$$                  ser = ser + xn*(ln1/(2*j+1)-fb(j))
c$$$c                  if (ldebug) then
c$$$c                     print*," n = ",j," a_n = ",ln1/(2*j+1)-fb(j)
c$$$c                  endif
c$$$                  xn = xn*x
c$$$                  j = j+1
c$$$               enddo
c$$$               if (ldebug) print*," C_i = ",2*ser*fac
c$$$               resr = resr + 2*ser*fac
c$$$            enddo
c$$$            if (ldebug) then
c$$$               print*," C0t via Vli2 calls: ",C0t
c$$$               print*," C0t via Taylor exp.",dcmplx(resr,0d0)
c$$$            endif
c$$$            C0t = dcmplx(resr,0d0)
c$$$         endif
c$$$      elseif (izero.eq.1) then
c$$$         ln1 = log(abs(qsq(1))/musq)
c$$$         ln2 = log(abs(qsq(2))/musq)
c$$$         res = 0.5d0*(ln1**2-ln2**2)
c$$$         if (qsq(1).gt.0d0) then
c$$$            res = res - ipi*ln1 - pi2o2
c$$$         endif
c$$$         if (qsq(2).gt.0d0) then
c$$$            res = res + ipi*ln2 + pi2o2
c$$$         endif
c$$$         C0t = res/(qsq(2)-qsq(1))
c$$$      elseif (izero.eq.2) then
c$$$         ln1 = log(abs(qsq(1))/musq)
c$$$         res = 0.5d0*ln1**2-pi2o6
c$$$         if (qsq(1).gt.0d0) then
c$$$            res = res - ipi*ln1 - pi2o2
c$$$         endif
c$$$         C0t = -res/qsq(1)
c$$$      else
c$$$         print*," WARNING: C0t called with 3 zero arguments "
c$$$         print*," C0t reset to 0 "
c$$$         C0t = 0
c$$$      endif
c$$$      if (ldebug) then
c$$$         print*," output: C0t = ",C0t
c$$$      endif
c$$$      return 
c$$$      end
c
c------------- D0(s,t,q1sq,q2sq)  massless 4-point function -----------
c
      complex*16 function D0t(s,t,q1sq,q2sq)
      implicit none
      double precision s, t, q1sq, q2sq
      double precision sl, tl, q1sql, q2sql
      double precision p1sq, p2sq, p3sq, p1p2, p1p3, p2p3

c evaluate the finite part of the scalar 4-point function for zero mass
c propagators 
c  
c  D0 = 1/(i*pi^2) * Int d^dk [-k^2-i0]^-1 [-(k+p1)^2-i0]^-1 
c                      [-(k+p1+p2)^2-i0]^-1 [-(k+p1+p2+p3)^2-i0]^-1
c
c       = pi^-eps Gamma(1+eps) (-s)^-eps [IR + D0t(s,t,q1sq,q2sq)]
c
c where p1^2 = 0 = (p1+p2+p3)^2 is assumed. The arguments of D0t are
c 
c   s = (p2+p3)^2, t = (p1+p2)^2, q1sq = p3^3 and q2sq = p2^2
c
c IR represents the divergent terms and is given by
c 
c  IR = 1/(s*t)*[1/eps**2 + 1/eps (log(q1sq/t)+log(q2sq/t)) ]
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
      complex*16 ipi, ieps, z1, z2, zlog, vli2, D0t4
      parameter (ipi=(0d0,3.14159 26535 89793d0),ieps=(0d0,1d-16))
      external vli2
      double precision pi,pi2o2,pi2o6
      parameter (pi=3.14159 26535 89793d0,pi2o2=pi**2/2d0)
      parameter (pi2o6=pi**2/6d0)
      logical ldebug
      parameter (ldebug=.false.)
      double precision rlog, rz1, rz2 
c
      sl = s
      tl = t
      q1sql = q1sq
      q2sql = q2sq
      goto 10
c
      entry D0t4(p1sq, p2sq, p3sq, p1p2, p1p3, p2p3)
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
         D0t = 0
         return 
      endif
c
      rz1 = tl/q1sql
      rz2 = tl/q2sql
      rlog = -log( abs(rz1*rz2) )
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
      D0t = (0.5d0*zlog**2 + 2*(vli2(z1)+vli2(z2)) - pi2o6)/(sl*tl)
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
         print*," D0t = ",d0t
      endif
      return 
      end
