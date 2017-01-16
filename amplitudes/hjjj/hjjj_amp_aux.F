c
c---------------- BCD_fill_v(k1,k2,q1,q2) -------------------------------
c
      subroutine BCD_fill_v(k1,k2,q1,q2)
      implicit none
      double precision k1(0:3),k2(0:3),q1(0:4),q2(0:4)
c
c     q1s = 0 for  the massless gluon
c
c This subroutine does the same as BCD_fill for a single momentum
c assignment as in the diagram
c
c     ID = 1  D0(k2,q2,q1),Dij(k2,q2,q1)
c     k1  -->---------->------------- k2
c                 $         $
c                 $         $
c                 $ q1      $  q2
c                 $         $ 
c               gluon     vector       
c 
c     ID = 2  D0(k2,q1,q2),Dij(k2,q1,q2)
c     k1  -->---------->------------- k2
c                 $         $
c                 $         $
c                 $ q2      $  q1
c                 $         $ 
c               vector     gluon      
c
c     ID = 3 D0(q1,k2,q2),Dij(q1,k2,q2)
c
c
c results are stored in Dijv(*,*,ID) etc. i.e.
c
c     renormalization scale
c
#include "VBFNLO/utilities/scales.inc"
c
Common block for the output
      complex*16 D0v(3), Dijv(3,13,3), Teps1(2), Teps2(2), Tborn(2,2),
     $     Tg1(2),Tg2(2),Ceps(2)
      common /bcd_qqvg/ D0v,Dijv,Teps1,Teps2,Tborn,Tg1,Tg2,Ceps
      double precision s, t, q1sq, q2sq,u
c
c local variables, external functions etc.
      double precision p1p2, p1p3, p2p3
      double precision musq     !renormalization scale
      complex*16 C0123(3),  C0124(3),  C0134(3),  C0234(3)
      complex*16 Cij123(2,4), Cij124(2,4), 
     1           Cij134(2,4), Cij234(2,4)
      complex*16 B0a(3), B0aq1, B0aq2, B0at, B0as,B0au
      complex*16 B0t, C0t, D0t1m
      external B0t, C0t, D0t1m
      double precision tmq2i,umq2i
      double precision dot0p, psumsq
      external dot0p, psumsq
      double precision pi,pi2o3
      parameter (pi=3.14159 26535 89793d0,pi2o3=pi**2/3d0)

c
      s = -2*dot0p(k1,k2)
      t = psumsq(k2,q2)
      u = psumsq(k2,q1)
      q1sq = dot0p(q1,q1)        !q1s = 0 
c      print*,'q1sq=',q1sq
      q2sq = dot0p(q2,q2)
c      print*,'q2sq=',q2sq
      musq = mursq(1,1)        !set in scales.inc
c
c now determine the required B0, C0 and D0 functions
c D0v(1) = D0t(s,t,0,q2sq,musq) = D0t(s,t,q2sq,0,musq)
c D0v(2) = D0t(s,u,0,q2sq,musq) = D0t(s,u,q2sq,0,musq)
      D0v(1) = D0t1m(s,t,q2sq,musq)
      D0v(2) = D0t1m(s,u,q2sq,musq)
      D0v(3) = D0t1m(t,u,q2sq,musq)
c ID=1
c C0123 = C0(0,q2sq,t), C0124 = C0(0,s,0), C0134 = C0(t,0,0),
c C0234 = C0(0,q2sq,s)
c ID=2
c C0123 = C0(0,0,u), C0124 =C0(0,s,0) , C0134 = C0(u,q2s,0),
c C0234 = C0(0,q2s,s)
c ID=3 
c  C0123 = C0(0,0,u), C0124 =C0(0,t,0) , C0134 = C0(u,q2s,0),
c C0234 = C0(0,q2s,t)
      
      C0124(1) = C0t(0d0,s,0d0,musq)
      C0234(1) = C0t(q2sq,0d0,s,musq)
      C0123(1) = C0t(0d0,q2sq,t,musq)
      C0134(1) = C0t(t,0d0,0d0,musq) !
c
      C0124(2) = C0t(0d0,s,0d0,musq) !=C0124(1)
      C0234(2) = C0t(0d0,q2sq,s,musq) !=C0234(1)
      C0123(2) = C0t(0d0,0d0,u,musq)
      C0134(2) = C0t(u,q2sq,0d0,musq)
c     
      C0124(3) = C0t(0d0,t,0d0,musq) !=C0134(1) 
      C0234(3) = C0t(0d0,q2sq,t,musq) !=C0123(1) 
      C0123(3) = C0t(0d0,0d0,u,musq) !=C0123(2)
      C0134(3) = C0t(u,q2sq,0d0,musq) !=C0134(2)    
c
c     and now the B0 functions: B0t(q2sq), B0t(t)
c     Note: B0t(s,musq=s) = 2
      B0aq1 = 0.0d0
      B0aq2 = B0t(q2sq,musq)
      B0at  = B0t(t,musq)
      B0as  = B0t(s,musq)
      B0au  = B0t(u,musq)
c
c Now everything is ready to calculate the C_ij and D_ij functions
c     ID = 1
c     p1 = k2,p2 = q2,p3 = q1
      B0a(1) = 0
      B0a(2) = B0aq2
      B0a(3) = B0at
      call tens3(0d0,0d0,q2sq,t,B0a,C0123(1),Cij123(1,1))
      B0a(1) = 0d0
      B0a(2) = B0as
      B0a(3) = 0d0
      call tens3(0d0,0d0,s,0d0,B0a,C0124(1),Cij124(1,1))
      B0a(1) = B0at
      B0a(2) = 0d0
      B0a(3) = 0d0
      call tens3(0d0,t,0d0,0d0,B0a,C0134(1),Cij134(1,1))
      B0a(1) = b0aq2
      B0a(2) = 0d0
      B0a(3) = B0as
      call tens3(0d0,q2sq,0d0,s,B0a,C0234(1),Cij234(1,1))
      p1p2 = dot0p(k2,q2)
      p1p3 = dot0p(k2,q1)
      p2p3 = dot0p(q2,q1)
      call tens4( 0d0, 0d0, q2sq, 0d0, p1p2, p1p3, p2p3, 
     &            C0123(1),  C0124(1),  C0134(1),  C0234(1),
     &            Cij123(1,1), Cij124(1,1), 
     &            Cij134(1,1), Cij234(1,1), 
     &            D0v(1), Dijv(1,1,1) )
c
c     ID = 2
c     p1 = k2,p2=q1,p3=q2 
c
      B0a(1) = 0
      B0a(2) = 0
      B0a(3) = B0au
      call tens3(0d0,0d0,0d0,u,B0a,C0123(2),Cij123(1,1))
      B0a(1) = 0d0
      B0a(2) = B0as
      B0a(3) = 0d0
      call tens3(0d0,0d0,s,0d0,B0a,C0124(2),Cij124(1,1))
      B0a(1) = B0au
      B0a(2) = B0aq2
      B0a(3) = 0d0
      call tens3(0d0,u,q2sq,0d0,B0a,C0134(2),Cij134(1,1))
      B0a(1) = 0d0
      B0a(2) = B0aq2
      B0a(3) = B0as
      call tens3(0d0,0d0,q2sq,s,B0a,C0234(2),Cij234(1,1))
      p1p2 = dot0p(k2,q1)
      p1p3 = dot0p(k2,q2)
      p2p3 = dot0p(q2,q1)
      call tens4( 0d0, 0d0, 0d0, q2sq, p1p2, p1p3, p2p3, 
     &            C0123(2),  C0124(2),  C0134(2),  C0234(2),
     &            Cij123(1,1), Cij124(1,1), 
     &            Cij134(1,1), Cij234(1,1), 
     &            D0v(2), Dijv(1,1,2) )
c
c      ID = 3
c     p1 = q1,p2=k2,p3=q2
c
      B0a(1) = 0
      B0a(2) = 0
      B0a(3) = B0au
      call tens3(0d0,0d0,0d0,u,B0a,C0123(3),Cij123(1,1))
      B0a(1) = 0d0
      B0a(2) = B0at
      B0a(3) = 0d0
      call tens3(0d0,0d0,t,0d0,B0a,C0124(3),Cij124(1,1))
      B0a(1) = B0au
      B0a(2) = B0aq2
      B0a(3) = 0d0
      call tens3(0d0,u,q2sq,0d0,B0a,C0134(3),Cij134(1,1))
      B0a(1) = 0d0
      B0a(2) = B0aq2
      B0a(3) = B0at
      call tens3(0d0,0d0,q2sq,t,B0a,C0234(3),Cij234(1,1))
      p1p2 = dot0p(k2,q1)
      p1p3 = dot0p(q2,q1)
      p2p3 = dot0p(k2,q2)
      call tens4( 0d0, 0d0, 0d0, q2sq, p1p2, p1p3, p2p3, 
     &            C0123(3),  C0124(3),  C0134(3),  C0234(3),
     &            Cij123(1,1), Cij124(1,1), 
     &            Cij134(1,1), Cij234(1,1), 
     &            D0v(3), Dijv(1,1,3) )
c     
c     
c and the coefficients of the tensors in M_6:
c     
c     check these formulas
      tmq2i = 1d0/(t-q2sq)
      umq2i = 1d0/(u-q2sq)
      
      Teps1(1) = (2*B0at+1)/t
      Teps2(1) = ( (B0at-b0aq2)*(2*t+3*q2sq)*tmq2i
     &     +2*b0aq2+1-2*q2sq*C0123(1) ) * tmq2i
      Tg1(1) = B0at/t
      Tg2(1) = (B0at-B0aq2)*tmq2i

      Ceps(1) = 2*(t*C0134(1)+1d0)/t
      Ceps(2) = 2*(u*C0123(2)+1d0)/u
c
      Tborn(1,1) =  ( 2*q2sq*(B0at-b0aq2) + t*B0at
     &     -q2sq*b0aq2 )*tmq2i - 2*q2sq*C0123(1)
     &     +2*B0at -5d0 +pi2o3 !-5d0 + pi2o3  
c
c      print*,'C0123-c0123=',C0123(1)-C0123(1)
      Teps1(2) = (2*B0au+1)/u
      Teps2(2) = ( (B0au-b0aq2)*(2*u+3*q2sq)*umq2i
     &     +2*b0aq2+1-2*q2sq*C0134(2) ) * umq2i
      Tg1(2) = B0au/u
      Tg2(2) = (B0au-B0aq2)*umq2i
      Tborn(2,1) =  ( 2*q2sq*(B0au-b0aq2) + u*B0au
     &     -q2sq*b0aq2 )*umq2i - 2*q2sq*C0134(2)
     &     +2*B0au -5d0 + pi2o3 !-5d0 + pi2o3 
c 
c     5-pi2o3 + cvirt = -3 see expression for div piece  
c     
      Tborn(1,2) = -2*(t*C0134(1)+1) +( 2*q2sq*(B0at-b0aq2) + t*B0at
     &     -q2sq*b0aq2 )*tmq2i - 2*q2sq*C0123(1)
     &     + B0at -5d0 + pi2o3 !-5d0 + pi2o3 
c
      Tborn(2,2) = -2*(u*C0123(2)+1) +( 2*q2sq*(B0au-b0aq2) + u*B0au
     &     -q2sq*b0aq2 )*umq2i - 2*q2sq*C0134(2)
     &     + B0au -5d0 + pi2o3 !-5d0 + pi2o3
c      print*,'bcd_fill'
c      print*,'Teps1=',Teps1
c      print*,'Teps2=',Teps2
c      print*,'Tg1=',Tg1
c      print*,'Tg2=',Tg2
c      print*,'Tborn(1,1)=',Tborn(1,1)
c      print*,'Tborn(2,1)=',Tborn(2,1)
c      print*,'Tborn(1,2)=',Tborn(1,2)
c      print*,'Tborn(2,2)=',Tborn(2,2)
      return
      end
c$$$c
c$$$c----  B0t(qsq): regularized 2-point function for massless loop  -----
c$$$c
c$$$      complex*16 function B0t(qsq,musq)
c$$$      implicit none
c$$$      double precision qsq, musq
c$$$
c$$$c evaluate the finite piece of the massless B0 function
c$$$c  
c$$$c             d^d k 
c$$$c    B0 = Int------- [-k^2+idelta]^-1 [-(k-q)^2+idelta]^-1 
c$$$c             i pi^2
c$$$c
c$$$c       = pi^-eps Gamma(1+eps) musq^-eps [1/eps + B0t(qsq,musq)]
c$$$c
c$$$c i.e. B0t(qsq,musq) = 2 - log( (-qsq-idelta)/musq )
c$$$c
c$$$c subroutine requires positive musq
c$$$c
c$$$c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c$$$c	Initial version:  2002 September 13
c$$$c	Last modified:    2002 September 13
c$$$c  
c$$$      complex*16 ipi
c$$$      parameter (ipi=(0d0,3.14159 26535 89793d0))
c$$$      logical ldebug
c$$$      parameter (ldebug=.true.)
c$$$c
c$$$      if (ldebug) then
c$$$         if (musq.le.0d0) then
c$$$            print*,' Unacceptable mu^2 in B0t(q^2,mu^2): mu^2 = ',
c$$$     1             musq,' ==> set B0t = 0 '
c$$$            b0t=0
c$$$            return
c$$$         endif
c$$$      endif
c$$$c
c$$$      if (qsq.lt.0d0) then
c$$$         b0t = 2 - log(-qsq/musq)
c$$$      else
c$$$         b0t = 2 - log(qsq/musq) + ipi
c$$$      endif
c$$$      end
c$$$c
c$$$c------------- C0(q1^2,q2^2,P^2)  massless 3-point fuction ------------
c$$$c
c$$$      complex*16 function C0t(q1sq,q2sq,Psq,musq)
c$$$      implicit none
c$$$      double precision q1sq, q2sq, psq, musq
c$$$
c$$$c evaluate the finite part of the scalar 3-point function for zero mass
c$$$c propagators 
c$$$c  
c$$$c  C0 = 1/(i*pi^2) * Int d^dk [-k^2-i0]^-1 [-(k+q1)^2-i0]^-1 [-(k-q2)^2-i0]^-1
c$$$c
c$$$c       = pi^-eps Gamma(1+eps) musq^-eps [IR + C0t(q1sq,q2sq,psq)]
c$$$c
c$$$c where IR represents the potentially divergent terms when one or more
c$$$c of the arguments vanish. P = q1+q2. The divergent terms are
c$$$c 
c$$$c (1) IR = 0              all 3 of q1sq,q2sq,psq nonzero
c$$$c 
c$$$c (2) IR = 1/(t-qsq)*log(t/qsq)*1/eps    for 1 vanishing argument, the
c$$$c                                        other two being t and qsq
c$$$c
c$$$c (3) IR = -1/qsq*(1/eps**2 - 1/eps log(qsq/musq))  when exactly one 
c$$$c                                                   argument, qsq,
c$$$c                                                   is nonzero
c$$$c
c$$$c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c$$$c	Initial version:  2002 September 13
c$$$c	Last modified:    2002 October 25
c$$$c  
c$$$      complex*16 ipi, z, res, vli2
c$$$      parameter (ipi=(0d0,3.14159 26535 89793d0))
c$$$      external vli2
c$$$      double precision pi,pi2o3,pi2o6,pi2o2
c$$$      parameter (pi=3.14159 26535 89793d0,pi2o3=pi**2/3d0)
c$$$      parameter (pi2o6=pi**2/6d0,pi2o2=pi**2/2d0)
c$$$      logical ldebug, ltaylor, first_call
c$$$      parameter (ldebug=.false.)
c$$$      double precision x, flambda, tsqlami, qsq(3), 
c$$$     1                 sum, rz, logz, ln1, ln2,
c$$$     2                 musq_t, asum, fac, resr, ser, fb(8), xn
c$$$      integer izero, inz, i, j
c$$$      data first_call /.true./
c$$$      save first_call, fb
c$$$c
c$$$      if (ldebug) then
c$$$         if (musq.le.0d0) then
c$$$            print*,' Unacceptable mu^2 in C0t(q1sq,q2sq,psq,mu^2): ',
c$$$     1         'mu^2 = ',musq,' ==> set C0t = 0 '
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
c$$$         print('(a,3g12.4,a,g12.4)'), 'qsq_i = ',q1sq,q2sq,psq,
c$$$     1                                ' mu^2 = ',musq
c$$$         print*,' number of q_i^2 = 0  ',izero
c$$$         print*,' number of q_i^2 ne 0 ',inz,' Sum = ',inz+izero
c$$$         if ((inz+izero).ne.3) print*,' WARNING: inz+izero .ne. 3 '
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
c$$$            print*,' lambda = ',flambda/4,' x_taylor = ',
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
c$$$               if (ldebug) print*,' i,x, z = ',i,x,rz
c$$$            enddo
c$$$            C0t = res*tsqlami
c$$$         elseif (flambda.lt.0d0 .and. .not.ltaylor) then  ! the lambda < 0 case
c$$$            tsqlami = 1/sqrt(-flambda)         ! + imag part of 2 sqrt(lambda)
c$$$            do i=1,3
c$$$               x = (sum-2*qsq(i))*tsqlami
c$$$               z = dcmplx(x,-1d0)/dcmplx(x,1d0)
c$$$               if (ldebug) then  
c$$$                  print*,' i,x, z = ',i,x,z
c$$$               endif
c$$$               res = res + dimag(Vli2(z))
c$$$            enddo
c$$$            C0t = res*2*tsqlami
c$$$         endif
c$$$         if (ldebug) then
c$$$            print*,' C0t via Vli2 calls: ',C0t
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
c$$$c               if (ldebug) print*,' i = ',i,' x = ',x
c$$$               ln1 = log(abs(qsq(i))/musq_t)
c$$$               ser = ln1
c$$$               j = 1
c$$$               xn = x
c$$$c               if (ldebug) print*,' log = ',ln1
c$$$               do while (abs(xn).gt.1d-15 .and. j.lt.9)
c$$$                  ser = ser + xn*(ln1/(2*j+1)-fb(j))
c$$$c                  if (ldebug) then
c$$$c                     print*,' n = ',j,' a_n = ',ln1/(2*j+1)-fb(j)
c$$$c                  endif
c$$$                  xn = xn*x
c$$$                  j = j+1
c$$$               enddo
c$$$               if (ldebug) print*,' C_i = ',2*ser*fac
c$$$               resr = resr + 2*ser*fac
c$$$            enddo
c$$$            if (ldebug) then
c$$$               print*,' C0t via Vli2 calls: ',C0t
c$$$               print*,' C0t via Taylor exp.',dcmplx(resr,0d0)
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
c$$$         print*,' WARNING: C0t called with 3 zero arguments '
c$$$         print*,' C0t reset to 0 '
c$$$         C0t = 0
c$$$      endif
c$$$      if (ldebug) then
c$$$         print*,' output: C0t = ',C0t
c$$$      endif
c$$$      return 
c$$$      end
c
c------------- D0(s,t,q1sq,q2sq)  massless 4-point function -----------
c
      complex*16 function D0t2m(s,t,q1sq,q2sq)
      implicit none
      double precision s, t, q1sq, q2sq,musq
      double precision sl, tl, q1sql, q2sql,qsql


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
C
c For the 1 mass box:
c D0 = 1/(i*pi^2) * Int d^dk [-k^2-i0]^-1 [-(k+p1)^2-i0]^-1 
c                      [-(k+p1+p2)^2-i0]^-1 [-(k+p1+p2+p3)^2-i0]^-1
c
c       = pi^-eps Gamma(1+eps) (musq)^-eps [IR + D0t(s,t,q1sq,musq)]
c divergent terms:
c
c  IR = 2/(s*t)*[1/eps**2 + 1/eps (log(qsq/t)+log(-musq/s)]
c
c Alternatively, arguments as in the call of tens4 can be used, i.e.
c
c  D0t(s,t,q1sq,q2sq) = D0t4(p1sq, p2sq, p3sq, p1p2, p1p3, p2p3) 
c
c This code agrees with the one written by Carlo. Checked 10/15/2002
c       Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2002 October 5
c	Last modified:    2005 June 27 by Terrance Figy
c
c       needs to checked numerically
c  
      complex*16 ipi, ieps, z1, z2, zlog, vli2,D0t1m
      complex*16 zlog1
      complex*16 zslog,ztlog,zqsqlog
      parameter (ipi=(0d0,3.14159 26535 89793d0),ieps=(0d0,1d-16))
      external vli2
      double precision pi,pi2o2,pi2o6,tpi2o3
      parameter (pi=3.14159 26535 89793d0,pi2o2=pi**2/2d0)
      parameter (pi2o6=pi**2/6d0,tpi2o3=2d0*pi**2/3d0)
      logical ldebug
      parameter (ldebug=.false.)
      double precision rlog, rz1, rz2,rz,rlog1
      double precision rslog,rtlog,rqsqlog,rsomu,rtomu,rqsqomu
c     new variables
      integer ioption
c
      
      sl = s
      tl = t
      q1sql = q1sq
      q2sql = q2sq
      ioption = 2
      goto 10
      entry D0t1m(s,t,q1sq,musq)
      sl=s
      tl=t
      qsql = q1sq
      ioption=1
      goto 10
c
c$$$*      entry D0t4(p1sq, p2sq, p3sq, p1p2, p1p3, p2p3)
c$$$      if (p1sq.ne.0d0 .or. 
c$$$     1    abs(2*(p1p2+p2p3+p1p3)/(p2sq+p3sq)+1).gt.1d-12) then
c$$$         print*,' Dont use D0t4 for for less than 2 massless momenta '
c$$$         print*,' (p1p2+p2p3+p1p3)/(p2sq+p3sq) = ',
c$$$     1           (p1p2+p2p3+p1p3)/(p2sq+p3sq)
c$$$         stop
c$$$      endif
c$$$      sl = p2sq+2*p2p3+p3sq
c$$$      tl = p2sq+2*p1p2
c$$$      q1sql = p3sq
c$$$      q2sql = p2sq
 10   continue

      if(ioption.eq.2) then ! 2 mass opposite box
         rz1 = tl/q1sql
         rz2 = tl/q2sql
         rlog = -log( abs(rz1*rz2) )
         zlog = rlog
         if (rz1.gt.0) then
            z1 = 1-rz1
         elseif (tl.lt.0d0) then !   t<0, q1sq>0
            z1 = 1-rz1-ieps
            zlog = zlog - ipi
         else                   !   t>0, q1sq<0
            z1 = 1-rz1+ieps
            zlog = zlog + ipi
         endif
         if (rz2.gt.0) then
            z2 = 1-rz2
         elseif (tl.lt.0d0) then !   t<0, q2sq>0
            z2 = 1-rz2-ieps
            zlog = zlog - ipi
         else                   !   t>0, q2sq<0
            z2 = 1-rz2+ieps
            zlog = zlog + ipi
         endif
         D0t2m = (0.5d0*zlog**2 + 2*(vli2(z1)+vli2(z2)) - pi2o6)/(sl*tl)
      elseif(ioption.eq.1) then ! 1-mass box
         if(musq.lt.0d0) then
            print*,'musq is less than 0'
            print*,'D0 is set to 1'
            D0t1m = 0d0
            return
         endif
         rz = tl/sl
         rz1 = qsql/sl
         rz2 = qsql/tl
         rsomu = -sl/musq
         rtomu = -tl/musq
         rqsqomu = -qsql/musq
c
         rlog = log(abs(rz))
         rlog1 = log(abs(rz1))
         zlog = rlog
         zlog1 = rlog1

         rslog = log(abs(rsomu))
         rtlog = log(abs(rtomu))
         rqsqlog = log(abs(rqsqomu))

         zslog =  rslog
         ztlog =  rtlog
         zqsqlog =  rqsqlog

         if(rsomu.lt.0d0) then
            zslog = zslog -ipi
         endif
         if(rtomu.lt.0d0) then
            ztlog = ztlog -ipi
         endif
         if(rqsqomu.lt.0d0) then
            zqsqlog = zqsqlog -ipi
         endif
         
         if(rz2.gt.0) then      
            z2 = 1-rz2
         elseif(tl.lt.0d0) then !q2sq >0,t<0
            z2=1-rz2 +ieps
         else                   !q2sq<0,t>0
            z2=1-rz2 -ieps
         endif
         if(rz1.gt.0) then      
            z1 = 1-rz1            
         elseif(sl.lt.0d0) then !q2sq >0,s<0
            zlog1 = zlog1 - ipi
            z1=1-rz1 +ieps
         else                   !q2sq<0,s>0
            zlog1 = zlog1 + ipi
            z1=1-rz1 -ieps
         endif
         if(rz.lt.0d0) then
            if(tl.gt.0d0) then       !s >0,t<0
               zlog = zlog - ipi
            else                !s<0,t>0
               zlog = zlog + ipi 
            endif
         endif
         if(ldebug) then
            D0t1m = (zlog**2 + zlog - zlog1**2 -2*vli2(z1)-
     1           2*vli2(z2)-tpi2o3)/(sl*tl)
            print*,'-----------------------'
            print*,'musq=',musq
            print*,'D0t1m with musq = -s fixed'
            print*,D0t1m 
         endif
         D0t1m = (zlog + zslog**2 + ztlog**2 - zqsqlog**2 
     1        - 2*Vli2(z1)-2*Vli2(z2)-tpi2o3)/(sl*tl)
         if(ldebug) then
            print*,'D0t1m with musq arbitrary'
            print*,D0t1m 
            print*,'-----------------------'
         endif
c
c     check this formula: it is vital
      else
         print*,'invalid selection in D0'
         stop
      endif
c
c      if (ldebug) then
c         print('(a,2g14.4,a,g14.4)'), 
c     1         ' t/qsq ratios: ',rz1, rz2,' t = ',tl
c         print('(a,f8.4,a,f8.4)'),
c     1         ' zlog = ',dreal(zlog),'+ i pi * ',dimag(zlog)/pi
c         z1 = 0.5d0*zlog**2
c         rz1 = (dreal(z1)-0.5d0*rlog**2)/pi**2
c         if (rlog.ne.0d0) then
c            rz2 = dimag(z1)/pi/rlog
c         else
c            rz2 = 0d0
c         endif
c         print('(2a,f8.4,a,f8.4,a)'),' extra terms from log: ',
c     1   ' pi**2*(',rz1,')+ i pi rlog *(',rz2,')' 
c         print*,' D0t = ',d0t
c      endif
      return 
      end
c$$$c
c$$$c-------------  spence function  vli2(z) ------------------------
c$$$c
c$$$      complex*16 function vli2(zin)
c$$$      implicit none
c$$$      complex*16 zin, z, u, u2, unpo, ans, zext
c$$$      double precision r, r2, r2n, fac
c$$$c
c$$$c determine the value of the dilogarithm 
c$$$c
c$$$c    vli2(z) = - int_0^1  log(1-zt)/t dt  with cut along the positive 
c$$$c                                        real axis, z>1
c$$$c
c$$$c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c$$$c	Initial version:  2000 November 6
c$$$c	Last modified:    2000 November 12
c$$$c  
c$$$      integer i
c$$$      double precision c0,c1,c2,c4,c6,c8,c10,c12,c14,c16,c18,c20,c22
c$$$      double precision b0,b1,b2,b4,b6,b8,b10,b12,b14,b16,b18,b20,b22
c$$$      double precision d0,d1,d2,d4,d6,d8,d10,d12,d14,d16,d18,d20,d22
c$$$      parameter (b0=1d0,            d0 =1d0,      c0= b0/d0)
c$$$      parameter (b1=-1d0/2d0,       d1 =d0*2d0,   c1= b1/d1)
c$$$      parameter (b2= 1d0/6d0,       d2 =d1*3d0,   c2= b2/d2)
c$$$      parameter (b4=-1d0/30d0,      d4 =d2*20d0,  c4= b4/d4)
c$$$      parameter (b6=1d0/42d0,       d6 =d4*42d0,  c6= b6/d6)
c$$$      parameter (b8=-1d0/30d0,      d8 =d6*72d0,  c8= b8/d8)
c$$$      parameter (b10=5d0/66d0,      d10=d8*110d0, c10=b10/d10)
c$$$      parameter (b12=-691d0/2730d0, d12=d10*156d0,c12=b12/d12)
c$$$      parameter (b14=7d0/6d0,       d14=d12*210d0,c14=b14/d14)
c$$$      parameter (b16=-3617d0/510d0, d16=d14*272d0,c16=b16/d16)
c$$$      parameter (b18=43867d0/798d0, d18=d16*342d0,c18=b18/d18)
c$$$      parameter (b20=-174611d0/330d0,d20=d18*420d0,c20=b20/d20)
c$$$      parameter (b22=854513d0/138d0,d22=d20*506d0,c22=b22/d22)
c$$$      double precision eps, epst, pi, pi2o6
c$$$      parameter (eps=1d-16, epst=1d-3)
c$$$      parameter (pi=3.14159 26535 89793d0, pi2o6=pi**2/6)
c$$$c
c$$$c debug information
c$$$      logical ldebug
c$$$      parameter (ldebug=.false.)
c$$$
c$$$      z = zin
c$$$c      print*,' vli2 call with z = ',z
c$$$      u = z**2
c$$$      r2 = dreal(z)**2+dimag(z)**2 
c$$$      if (r2.lt.eps) then
c$$$         vli2 = z + u/4d0
c$$$         return
c$$$      elseif (r2.lt.epst) then
c$$$         ans = z + u/4
c$$$         do i = 3,11
c$$$            u = u*z
c$$$            ans = ans + u/i**2
c$$$         enddo
c$$$         vli2 = ans
c$$$         return
c$$$      endif
c$$$      if (dreal(z).ge.1d0 .and. dimag(z).eq.0 ) then
c$$$         z = z + (0d0,1d0)*eps
c$$$      endif
c$$$c
c$$$c use z-->1/z and z--> 1-z mappings of the spence function to restrict 
c$$$c agument to unit circle in the complex plane with Re(z) <= 0.5
c$$$c
c$$$      zext = (0d0,0d0)
c$$$      fac = 1
c$$$      if (r2.gt.1d0) then     ! map z ---> 1/z
c$$$         fac = -fac
c$$$         zext = -pi2o6 - 0.5d0*(log(-z))**2
c$$$         z = 1/z
c$$$      endif
c$$$      if (dreal(z).gt.0.5d0) then     ! map new z ---> 1-z
c$$$         zext = zext + fac*(pi2o6-log(z)*log(1-z))
c$$$         fac = -fac
c$$$         z = 1-z
c$$$      endif
c$$$c
c$$$c now use t = 1 - exp(-u) mapping to write Li(z) in terms of Bernoulli 
c$$$c numbers
c$$$c
c$$$      u = - log(1-z)
c$$$      r2 = abs(u)**2
c$$$      u2 = u*u
c$$$      ans = u*(c0 + u*(c1+c2*u))
c$$$      r2n = r2*r2       !r^4
c$$$
c$$$      unpo = u2*u2*u
c$$$      ans = ans + c4*unpo
c$$$
c$$$      unpo = unpo*u2
c$$$      ans = ans + c6*unpo
c$$$
c$$$      r = r2n*r2n       !r^8
c$$$      unpo = unpo*u2
c$$$      ans = ans + c8*unpo
c$$$
c$$$      r2n = r*r2        !r^10 
c$$$      if ((r2n*c10).gt.eps) then
c$$$         unpo = unpo*u2
c$$$         ans = ans + c10*unpo
c$$$      else
c$$$         vli2 = fac * ans + zext
c$$$         if (ldebug) print*,' exit vli2s at n=8 '
c$$$         return
c$$$      endif
c$$$
c$$$      unpo = unpo*u2
c$$$      ans = ans + c12*unpo
c$$$
c$$$      unpo = unpo*u2
c$$$      ans = ans + c14*unpo
c$$$
c$$$      unpo = unpo*u2
c$$$      ans = ans + c16*unpo
c$$$
c$$$      r2n = r2n*r
c$$$      if ((r2n*c18).gt.eps) then
c$$$         unpo = unpo*u2
c$$$         ans = ans + c18*unpo
c$$$      else
c$$$         vli2 = fac * ans + zext
c$$$         if (ldebug) print*,' exit vli2s at n=16 '
c$$$         return
c$$$      endif
c$$$
c$$$      unpo = unpo*u2
c$$$      ans = ans + c20*unpo
c$$$
c$$$      unpo = unpo*u2
c$$$      ans = ans + c22*unpo
c$$$
c$$$      vli2 = fac * ans + zext
c$$$      end
c$$$
c$$$c
c$$$c----- psumsq(p1,p2) = (p1+p2).(p1+p2) -----------------
c$$$c
c$$$c      function psumsq(p1,p2)
c$$$c      implicit none
c$$$c      double precision psumsq,p1(0:3),p2(0:3)
c$$$c      psumsq = (p1(0)+p2(0))**2 - (p1(1)+p2(1))**2 -
c$$$c     1         (p1(2)+p2(2))**2 - (p1(3)+p2(3))**2
c$$$c      end
c$$$
c$$$c      real*8 function srt(x2)
c$$$c      real*8 x2
c$$$c      if (x2.gt.0d0) then
c$$$c         srt = sqrt(x2)
c$$$c      elseif (x2.lt.0d0) then
c$$$c         srt = -sqrt(-x2)
c$$$c      else
c$$$c         srt = 0d0
c$$$c      endif
c$$$c      end
c$$$
