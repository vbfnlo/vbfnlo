
c D0fin(s,t,p1sq,p2sq,p3sq,p4sq,musq)
c Author: Francisco Campanario
c Date:2008-04-16
c I have used Carlo Oleari D02m,D03m,D04m
c Bug solve in analytical continuation of D03m
c Need  my  D00m,D02me,D01m functions
c
c                     k
c     p1 ->-----------<-------------<-- p4
c             |                 |
c             |                 |
c             |                 |  
c     p2 ->-------------------------<-- p3
c
c
c    s=(p1+p2)^2
c    t=(p2+p3)^2
c    musq = mu^2 reference dimensional scale
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+q1)^2= 
c          N_ep * D0fin(s,t,q1^2,q2^2)
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)
c   the external momenta could be whaterver value


      function D0fin(s,t,p1sq,p2sq,p3sq,p4sq,musq)
      implicit none
      double complex D0fin
      double precision s,t,musq,p1sq,p2sq,p3sq,p4sq
      double complex D01m_fin,D02m_fin,D02m_fine,
     -     D03m_fin,D04m_fin,D00m_fin
      External D01m_fin,D02m_fin,D02m_fine,
     -     D03m_fin,D04m_fin,D00m_fin
      real * 8 tiny
      parameter (tiny=1d-6)
      

c      print*,'s = ',s
c      print*,'t = ',t
c      print*,'p1sq = ',p1sq
c      print*,'p2sq = ',p2sq
c      print*,'p3sq = ',p3sq
c      print*,'p4sq = ',p4sq
c      print*,'musq = ',musq

c      s =  -64972.6841d0
c      t =  -7028.19238d0
c      p1sq =   16848.6339d0
c      p2sq =  -6486.13781d0
c      p3sq =   6513.00115d0
c      p4sq =   0.0d0
c      musq =   3699.34781d0


      if(abs(p1sq).lt.tiny) then
         if(abs(p2sq).lt.tiny)then
            if(abs(p3sq).lt.tiny)then
               if(abs(p4sq).lt.tiny) then
c     D(s,t,0,0,0,0)                                                      
                  D0fin=D00m_fin(s,t,musq)
                  RETURN
               else
c     D(s,t,0,0,0,p4sq) 
                  D0fin=D01m_fin(s,t,p4sq,musq)
                  RETURN
               endif 
            else 
               if(abs(p4sq).lt.tiny) then
c     D(s,t,0,0,p3sq,0)
                  D0fin=D01m_fin(s,t,p3sq,musq)
                  RETURN 
               else
c     D(s,t,0,0,p3sq,p4sq)
                  D0fin=D02m_fin(s,t,p3sq,p4sq,musq)
                  RETURN
               endif
            endif
         else
            if(abs(p3sq).lt.tiny)then
               if(abs(p4sq).lt.tiny) then
c     D(s,t,0,p2sq,0,0)
                  D0fin=D01m_fin(s,t,p2sq,musq)
                  RETURN
               else
c     D(s,t,0,p2sq,0,p4sq)
                  D0fin=D02m_fine(s,t,p2sq,p4sq,musq) 
                  RETURN
               endif 
            else 
               if(abs(p4sq).lt.tiny) then
c     D(s,t,0,p2sq,p3sq,0)
                  D0fin=D02m_fin(t,s,p2sq,p3sq,musq) 
                  RETURN
               else 
c     D(s,t,0,p2sq,p3sq,p4sq)
                  D0fin=D03m_fin(s,t,p2sq,p3sq,p4sq,musq) 
               endif
            endif
         endif 
      else
         if(abs(p2sq).lt.tiny)then
            if(abs(p3sq).lt.tiny)then
               if(abs(p4sq).lt.tiny) then
c     D(s,t,p1sq,0,0,0) 
                  D0fin=D01m_fin(s,t,p1sq,musq)
                  RETURN
               else
c     D(s,t,p1sq,0,0,p4sq) 
                  D0fin=D02m_fin(t,s,p4sq,p1sq,musq) 
                  RETURN
               endif 
            else 
               if(abs(p4sq).lt.tiny) then
c     D(s,t,p1sq,0,p3sq,0) 
                  D0fin=D02m_fine(s,t,p1sq,p3sq,musq) 
                  RETURN
               else
c     D(s,t,p1sq,0,p3sq,p4sq) 
                  D0fin=D03m_fin(t,s,p3sq,p4sq,p1sq,musq)
                  RETURN
               endif
            endif
         else
            if(abs(p3sq).lt.tiny)then
               if(abs(p4sq).lt.tiny) then
c     D(s,t,p1sq,p2sq,0,0) 
                  D0fin=D02m_fin(s,t,p1sq,p2sq,musq)
                  RETURN
               else
c     D(s,t,p1sq,p2sq,0,p4sq) 
                  D0fin=D03m_fin(s,t,p4sq,p1sq,p2sq,musq) 
                  RETURN
               endif 
            else 
               if(abs(p4sq).lt.tiny) then
c     D(s,t,p1sq,p2sq,p3sq,0) 
                  D0fin=D03m_fin(t,s,p1sq,p2sq,p3sq,musq) 
                  RETURN
               else
c     D(s,t,p1sq,p2sq,p3sq,p4sq) 
                  D0fin=D04m_fin(s,t,p1sq,p2sq,p3sq,p4sq,musq)
                  RETURN
               endif
            endif
         endif
      endif
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCC                SCALAR INTEGRALS                CCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c    ***********   D00m_fin(s,t,musq)   *************
c
c Scalar box with MASSLESS PROPAGATORS and with external kinematics: 
c
c                     k
c     p1 ->-----------<-------------<-- p4
c             |                 |
c             |                 |
c             |                 |  
c     p2 ->-------------------------<-- p3
c
c
c    s = (p1+p2)^2
c    t = (p2+p3)^2
c    p1^2 = 0, p2^2 = 0, p3^2 = 0, p4^2 = 0
c    musq = mu^2 = reference dimensional scale 
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+p3)^2 = 
c          N_ep * D00m_fin(s,t,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

      function D00m_fin(s,t,musq)
      implicit none
      double complex D00m_fin
      double precision s,t,musq
      double precision ms,mt
      double complex lnms,lnmt,lnsot


      double complex ipi
      parameter (ipi=(0d0,3.14159265358979323846264338328d0))
      double precision pi,pi2, pi2o3t5
      parameter (pi=3.14159265358979323846264338328d0,
     -     pi2 = 9.86960440108935861883449099988d0,
     -   pi2o3t5=16.4493406684822643647241516664602519d0 )

      double complex ris
      double precision prefactor, theta, dilog, im_part
      external theta,im_part,dilog


      if (musq.lt.0d0) then
         write(*,*) 
     -        'POSSIBLE ERROR IN D00m_fin: SCALE MUSQ LESS THAN ZERO!!'
      endif

      prefactor = 1d0/(s*t)
      ms = -s/musq
      mt = -t/musq
      
      
      lnms = log(abs(ms)) - ipi*theta(-ms)
      lnmt = log(abs(mt)) - ipi*theta(-mt)
      lnsot = log(abs(s/t)) + ipi * im_part(s,t)

      ris = lnms**2+lnmt**2-lnsot**2-pi2o3t5

      D00m_fin = prefactor * ris
      end


c    ***********   D01m_fin(s,t,m4sq,musq)   *************
c
c Scalar box with MASSLESS PROPAGATORS and with external kinematics: 
c
c                     k
c     p1 ->-----------<-------------<-- p4
c             |                 |
c             |                 |
c             |                 |  
c     p2 ->-------------------------<-- p3
c
c
c    s = (p1+p2)^2
c    t = (p2+p3)^2
c    p1^2 = 0, p2^2 = 0, p3^2 = 0, p4^2 = m4sq <>0
c    musq = mu^2 = reference dimensional scale 
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+p3)^2 = 
c          N_ep * D01m_fin(s,t,m3sq,m4sq,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

      function D01m_fin(s,t,m4sq,musq)
      implicit none
      double complex D01m_fin
      double precision s,t,m4sq,musq
      double precision ms,mt,mm4sq
      double complex lnms,lnmt,lnmm4sq,lnsot
      double complex li2arg1,li2arg2
      double precision arg1,arg2
      double complex ipi
      parameter (ipi=(0d0,3.14159265358979323846264338328d0))
      double precision pi,pi2, pi2o3t2
      parameter (pi=3.14159265358979323846264338328d0,
     -     pi2 = 9.86960440108935861883449099988d0,
     -   pi2o3t2= 6.57973626739290574588966066658410076d0 )

      double complex ris
      double precision prefactor, theta, dilog, im_part
      external theta,im_part,dilog


      if (musq.lt.0d0) then
         write(*,*) 
     -        'POSSIBLE ERROR IN D01m_fin: SCALE MUSQ LESS THAN ZERO!!'
      endif

      prefactor = 1d0/(s*t)
      ms = -s/musq
      mt = -t/musq
      mm4sq = - m4sq/musq
      
      lnms = log(abs(ms)) - ipi*theta(-ms)
      lnmt = log(abs(mt)) - ipi*theta(-mt)
      lnmm4sq = log(abs(mm4sq)) - ipi*theta(-mm4sq)
      lnsot = log(abs(s/t)) + ipi * im_part(s,t)

      arg1 = 1-m4sq/s
      arg2 = 1-m4sq/t

c     (m4sq/s.lt.0d0)
      if (arg1.gt.1d0) then
         li2arg1 = -dilog(1d0/arg1) - log(arg1)**2/2+pi2/3 
     -        - ipi*log(arg1)*im_part(m4sq,s)
      else
         li2arg1 = dilog(arg1)
      endif

c     (m4sq/t.lt.0d0)
      if (arg2.gt.1d0)  then
         li2arg2 = -dilog(1d0/arg2) - log(arg2)**2/2+pi2/3 
     -        - ipi*log(arg2)*im_part(m4sq,t)
      else
         li2arg2 = dilog(arg2)
      endif

      ris = -lnmm4sq**2+lnms**2+lnmt**2-lnsot**2-2d0*(li2arg1+li2arg2)
     -   -pi2o3t2

      D01m_fin = prefactor * ris
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCC                SCALAR INTEGRALS                CCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c    ***********   D02m_fin(s,t,m3sq,m4sq,musq)   *************
c
c Scalar box with MASSLESS PROPAGATORS and with external kinematics: 
c
c                     k
c     p1 ->-----------<-------------<-- p4
c             |                 |
c             |                 |
c             |                 |  
c     p2 ->-------------------------<-- p3
c
c
c    s = (p1+p2)^2
c    t = (p2+p3)^2
c    p1^2 = 0, p2^2 = 0, p3^2 = m3sq <>0, p4^2 = m4sq <>0
c    musq = mu^2 = reference dimensional scale 
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+p3)^2 = 
c          N_ep * D02m_fin(s,t,m3sq,m4sq,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

      function D02m_fin(s,t,m3sq,m4sq,musq)
      implicit none
      double complex D02m_fin
      double precision s,t,m3sq,m4sq,musq
      double precision ms,mt,mm3sq,mm4sq
      double complex lnms,lnmt,lnmm3sq,lnmm4sq,lnsot
      double complex li2arg1,li2arg2
      double precision arg1,arg2
      double complex ipi
      parameter (ipi=(0d0,3.14159265358979323846264338328d0))
      double precision pi,pi2
      parameter (pi=3.14159265358979323846264338328d0,
     -     pi2 = 9.86960440108935861883449099988d0)

      double complex ris
      double precision prefactor, theta, dilog, im_part
      external theta,im_part,dilog
      if (musq.lt.0d0) then
         write(*,*) 
     -        'POSSIBLE ERROR IN D02m_fin: SCALE MUSQ LESS THAN ZERO!!'
      endif

      prefactor = 1d0/(s*t)
      ms = -s/musq
      mt = -t/musq
      mm3sq = - m3sq/musq
      mm4sq = - m4sq/musq
      
      lnms = log(abs(ms)) - ipi*theta(-ms)
      lnmt = log(abs(mt)) - ipi*theta(-mt)
      lnmm3sq = log(abs(mm3sq)) - ipi*theta(-mm3sq)
      lnmm4sq = log(abs(mm4sq)) - ipi*theta(-mm4sq)
      lnsot = log(abs(s/t)) + ipi * im_part(s,t)

      arg1 = 1-m3sq/t
      arg2 = 1-m4sq/t

c     (m3sq/t.lt.0d0)
      if (arg1.gt.1d0) then
         li2arg1 = -dilog(1d0/arg1) - log(arg1)**2/2+pi2/3 
     -        - ipi*log(arg1)*im_part(m3sq,t)
      else
         li2arg1 = dilog(arg1)
      endif

c     (m4sq/t.lt.0d0)
      if (arg2.gt.1d0)  then
         li2arg2 = -dilog(1d0/arg2) - log(arg2)**2/2+pi2/3 
     -        - ipi*log(arg2)*im_part(m4sq,t)
      else
         li2arg2 = dilog(arg2)
      endif

      ris = lnms**2/2+lnmt**2-lnmm3sq**2/2-lnmm4sq**2/2
     -     +lnmm3sq*lnmm4sq-(lnmm4sq+lnmm3sq-lnms)*lnms
     -     -2*li2arg1-2*li2arg2-lnsot**2 - pi2/6
        
      D02m_fin = prefactor * ris
      end


c    ***********   D02m_fine(s,t,m1sq,m3sq,musq)   *************
c
c Scalar box with MASSLESS PROPAGATORS and with external kinematics: 
c
c                     k
c     p1 ->-----------<-------------<-- p4
c             |                 |
c             |                 |
c             |                 |  
c     p2 ->-------------------------<-- p3
c
c
c    s = (p1+p2)^2
c    t = (p2+p3)^2
c    p1^2 = m1sq, p2^2 = 0, p3^2 = m3sq , p4^2 =0
c    musq = mu^2 = reference dimensional scale 
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+p3)^2 = 
c          N_ep * D02m_fine(s,t,m1sq,m3sq,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

      function D02m_fine(s,t,m1sq,m3sq,musq)
      implicit none
      double complex D02m_fine
      double precision s,t,m1sq,m3sq,musq
      double precision ms,mt,mm1sq,mm3sq
      double complex lnms,lnmt,lnmm1sq,lnmm3sq,lnsot,lnarg5
      double complex li2arg1,li2arg2,li2arg3,li2arg4,li2arg5
      double precision arg1,arg2,arg3,arg4,arg5
      double complex ipi
      parameter (ipi=(0d0,3.14159265358979323846264338328d0))
      double precision pi,pi2
      parameter (pi=3.14159265358979323846264338328d0,
     -     pi2 = 9.86960440108935861883449099988d0)

      double complex ris
      double precision prefactor, theta, dilog, im_part, sgn
      double precision im_a, im_b,im_ab, eta3

      external theta, dilog, im_part, sgn,eta3


      if (musq.lt.0d0) then
         write(*,*) 
     -        'POSSIBLE ERROR IN D02m_fine: SCALE MUSQ LESS THAN ZERO!!'
      endif

      prefactor = 1d0/(m1sq*m3sq- s*t)
      ms = -s/musq
      mt = -t/musq
      mm1sq = - m1sq/musq
      mm3sq = - m3sq/musq
      
      lnms = log(abs(ms)) - ipi*theta(-ms)
      lnmt = log(abs(mt)) - ipi*theta(-mt)
      lnmm1sq = log(abs(mm1sq)) - ipi*theta(-mm1sq)
      lnmm3sq = log(abs(mm3sq)) - ipi*theta(-mm3sq)
      lnsot = log(abs(s/t)) + ipi*im_part(s,t)
c      write(*,*) "lnsot", lnsot
 
      arg1 = 1-m1sq/t
      arg2 = 1-m3sq/t
      arg3 = 1-m1sq/s
      arg4 = 1-m3sq/s
      arg5= 1-(m1sq*m3sq)/(s*t)

c$$$      write(*,*) "arg1", arg1
c$$$      write(*,*) "arg2", arg2
c$$$      write(*,*) "arg3", arg3
c$$$      write(*,*) "arg4", arg4
c$$$      write(*,*) "arg5", arg5

c     (m1sq/t.lt.0d0)
      if (arg1.gt.1d0) then
         li2arg1 = -dilog(1d0/arg1) - log(arg1)**2/2+pi2/3 
     -        - ipi*log(arg1)*im_part(m1sq,t)
      else
         li2arg1 = dilog(arg1)
      endif
     

c     (m3sq/t.lt.0d0)
      if (arg2.gt.1d0)  then
         li2arg2 = -dilog(1d0/arg2) - log(arg2)**2/2+pi2/3 
     -        - ipi*log(arg2)*im_part(m3sq,t)
      else
         li2arg2 = dilog(arg2)
      endif


c     (m1sq/s.lt.0d0)
      if (arg3.gt.1d0) then
         li2arg3 = -dilog(1d0/arg3) - log(arg3)**2/2+pi2/3 
     -        - ipi*log(arg3)*im_part(m1sq,s)
      else
         li2arg3 = dilog(arg3)
      endif

c     (m3sq/s.lt.0d0)
      if (arg4.gt.1d0)  then
         li2arg4 = -dilog(1d0/arg4) - log(arg4)**2/2+pi2/3 
     -        - ipi*log(arg4)*im_part(m3sq,s)
      else
         li2arg4 = dilog(arg4)
      endif

      im_a = m1sq/s*(1d0/m1sq-1d0/s)
      im_b = m3sq/t*(1d0/m3sq-1d0/t)
      im_ab =m1sq*m3sq/(s*t)*(1d0/m1sq-1d0/s+1d0/m3sq-1d0/t)


c     m1sq*m3sq/(s*t)  < 0
      if (arg5.gt.1d0) then
         li2arg5 = -dilog(1d0/arg5) - log(arg5)**2/2+pi2/3 
     -        - ipi*log(arg5)*sgn(im_ab)
      else
         li2arg5 = dilog(arg5)
      endif

      if (arg5.lt.0d0) then
         lnarg5 = log(abs(arg5)) - ipi*sgn(im_ab)
      else
         lnarg5 = log(arg5)
      endif
      
      li2arg5 = li2arg5  + 2*ipi*eta3(im_a,im_b,im_ab)*lnarg5


c$$$      write(*,*) "li2arg1", li2arg1
c$$$      write(*,*) "li2arg2", li2arg2
c$$$      write(*,*) "li2arg3", li2arg3
c$$$      write(*,*) "li2arg4", li2arg4
c$$$      write(*,*) "li2arg5", li2arg5
c$$$      write(*,*) "********************"
c$$$      write(*,*) "eta3(im_a,im_b,im_ab)",  eta3(im_a,im_b,im_ab)
c$$$      write(*,*) "********************"
c$$$      write(*,*) "lnms", lnms
c$$$      write(*,*) "lnmt", lnmt
c$$$      write(*,*) "lnmms1sq", lnmm1sq
c$$$      write(*,*) "lnmms3sq", lnmm3sq
c$$$      write(*,*) "lnsot", lnsot
 
      ris = lnmm1sq**2+lnmm3sq**2-lnms**2-lnmt**2+lnsot**2+
     -  2.d0*(li2arg1+li2arg2+li2arg3+li2arg4-li2arg5)   

c      write(*,*) "ris", ris

      D02m_fine = prefactor*ris
      end


c    ***********   D03m_fin(s,t,m2sq,m3sq,m4sq,musq)   *************
c
c Scalar box with MASSLESS PROPAGATORS and with external kinematics: 
c
c
c                     k
c     p1 ->-----------<-------------<-- p4
c             |                 |
c             |                 |
c             |                 |  
c     p2 ->-------------------------<-- p3
c
c
c    s = (p1+p2)^2
c    t = (p2+p3)^2
c    p1^2 = 0, p2^2 = m2sq <>0, p3^2 = m3sq <>0, p4^2 = m4sq <>0
c    musq = mu^2 = reference dimensional scale 
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+p3)^2 =  
c          N_ep * D03m_fin(s,t,m2sq,m3sq,m4sq,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

      function D03m_fin(s,t,m2sq,m3sq,m4sq,musq)
      implicit none
      double complex D03m_fin
      double precision s,t,m2sq,m3sq,m4sq,musq
      double complex D03m_fin_BDK

      D03m_fin = D03m_fin_BDK(s,t,m2sq,m3sq,m4sq,musq)
c      D03m_finpaco = D03m_fin_DN(s,t,m2sq,m3sq,m4sq,musq)
      
      end


      function D03m_fin_BDK(s,t,m2sq,m3sq,m4sq,musq)
      implicit none
      double complex D03m_fin_BDK
      double precision s,t,m2sq,m3sq,m4sq,musq
      double precision ms,mt,mm2sq,mm3sq,mm4sq,sot
      double complex lnms,lnmt,lnmm2sq,lnmm3sq,lnmm4sq,lnsot,lnarg3
      double complex li2arg1,li2arg2,li2arg3
      double precision arg1,arg2,arg3
      double complex ipi
      parameter (ipi=(0d0,3.14159265358979323846264338328d0))
      double precision pi,pi2
      parameter (pi=3.14159265358979323846264338328d0,
     -     pi2 = 9.86960440108935861883449099988d0)

      double complex ris
      double precision prefactor, theta, dilog, sgn, im_part
      double precision im_a, im_b,im_ab, eta3

      external theta, dilog, sgn, im_part,eta3


      integer icount, bcount
      data icount /0/
      data bcount /0/
      save icount, bcount

      if (musq.lt.0d0) then
         write(*,*) 
     -        'POSSIBLE ERROR IN D03m_fin: SCALE MUSQ LESS THAN ZERO!!'
      endif


c  1/(ms*mt-mm2s*mm4s)*
c  (1/2*ln(ms)^2+1/2*ln(mt)^2-1/2*ln(mm2s)^2-1/2*ln(mm4s)^2-ln(ms/mt)^2
c  +ln(mm2s)*ln(mm3s)-(ln(mm3s)+ln(mm2s)-ln(mt))*ln(mt)
c  +ln(mm3s)*ln(mm4s)-(ln(mm4s)+ln(mm3s)-ln(ms))*ln(ms)
c  -2*li2(1-mm2s/ms)
c  -2*li2(1-mm4s/mt)
c  +2*li2(1-mm2s*mm4s/ms/mt))

      prefactor = 1d0/(s*t-m2sq*m4sq)
      ms = -s/musq
      mt = -t/musq
      mm2sq = - m2sq/musq
      mm3sq = - m3sq/musq
      mm4sq = - m4sq/musq
      sot = s/t
      
      lnms = log(abs(ms)) - ipi*theta(-ms)
      lnmt = log(abs(mt)) - ipi*theta(-mt)
      lnmm2sq = log(abs(mm2sq)) - ipi*theta(-mm2sq)
      lnmm3sq = log(abs(mm3sq)) - ipi*theta(-mm3sq)
      lnmm4sq = log(abs(mm4sq)) - ipi*theta(-mm4sq)
      lnsot = log(abs(sot)) + ipi * im_part(s,t)

     
      arg1 = 1-m2sq/s
      arg2 = 1-m4sq/t
      arg3 = 1-m2sq*m4sq/s/t
      

c     m2sq/s < 0 
      if (arg1.gt.1d0) then
         li2arg1 = -dilog(1d0/arg1) - log(arg1)**2/2+pi2/3 
     -        - ipi*log(arg1)*im_part(m2sq,s)
      else
         li2arg1 = dilog(arg1)
      endif

c     m4sq/t < 0
      if (arg2.gt.1d0) then
         li2arg2 = -dilog(1d0/arg2) - log(arg2)**2/2+pi2/3 
     -        - ipi*log(arg2)*im_part(m4sq,t)
      else
         li2arg2 = dilog(arg2)
      endif


      im_a = m2sq/s*(1d0/m2sq-1d0/s)
      im_b = m4sq/t*(1d0/m4sq-1d0/t)
      im_ab = m2sq*m4sq/s/t*(1d0/m2sq-1d0/s+1d0/m4sq-1d0/t)


c     m2sq*m4sq/(s*t)  < 0
c      print*,'arg3= ',arg3
      if (arg3.gt.1d0) then
         li2arg3 = -dilog(1d0/arg3) - log(arg3)**2/2+pi2/3 
     -        - ipi*log(arg3)*sgn(im_ab)
      else
         li2arg3 = dilog(arg3)
      endif

      if (arg3.lt.0d0) then
         lnarg3 = log(abs(arg3)) - ipi*sgn(im_ab)
      else
         lnarg3 = log(arg3)
      endif

c      print*,'li2arg3 = ',li2arg3      
      li2arg3 = li2arg3  + 2*ipi*eta3(im_a,im_b,im_ab)*lnarg3
c      print*,'li2arg3 = ',li2arg3  
c       bcount = bcount + 1
c      test = dimag(2*ipi*eta3(im_a,im_b,im_ab)*lnarg3)
c      if (test.ne.0.0d0) then
c         icount = icount + 1
c         print*,'icount/bcount = ',1.0d0*icount/bcount, icount, bcount
c      print*,'test = ',test
c      endif

      ris = lnms**2/2+lnmt**2/2-lnmm2sq**2/2-lnmm4sq**2/2
     -     +lnmm2sq*lnmm3sq
     -     -(lnmm3sq+lnmm2sq-lnmt)*lnmt+lnmm3sq*lnmm4sq
     -     -(lnmm4sq+lnmm3sq-lnms)*lnms      
      
      ris = ris - 2*li2arg1 - 2*li2arg2 + 2*li2arg3 
     -     - lnsot**2
        
      D03m_fin_BDK = prefactor * ris
      end


c    ***********   D04m_fin(s,t,m1sq,m2sq,m3sq,m4sq,musq)   *************
c
c Scalar box with MASSLESS PROPAGATORS and with external kinematics: 
c
c
c                     k
c     p1 ->-----------<-------------<-- p4
c             |                 |
c             |                 |
c             |                 |  
c     p2 ->-------------------------<-- p3
c
c
c    s = (p1+p2)^2
c    t = (p2+p3)^2
c    p1^2 = m1sq <>0, p2^2 = m2sq <>0, p3^2 = m3sq <>0, p4^2 = m4sq <>0
c    musq = mu^2 = reference dimensional scale  (NOT used here, since 
c                 this box is finite 
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+p3)^2 =  
c          N_ep * D04m_fin(s,t,m1sq,m2sq,m3sq,m4sq,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)
c

c	m2 ... arbitrary, but positive -> set m2 = 1d0


c
c
c     D0 with 4 nonzero external masses
c     and zero internal propagators      
c     use Eq.(41) from Denner, Nierste, Scharf		     
c

      function D04m_fin(s,t,m1sq,m2sq,m3sq,m4sq,musq)
      implicit none
      double complex D04m_fin
      double precision s,t,m1sq,m2sq,m3sq,m4sq,musq

      double complex ipi
      parameter (ipi=(0d0,3.14159265358979323846264338328d0))
      double precision pi,pi2
      parameter (pi=3.14159265358979323846264338328d0,
     -     pi2 = 9.86960440108935861883449099988d0)

      integer k

      double precision dilog, eta3, theta, sgn      

      double complex ln_mx(2),lnk12,lnk13,lnk14,lnk23,total

      double precision m2,k12,k13,k14,k23,k24,k34
      double precision a,b,c,d,sqrtdelta,x(2),qtemp
      double precision k34ok13,im_k34ok13,k24ok12,im_k24ok12
      double precision im_x(2),im_k34ok13_x(2),im_k24ok12_x(2),arg,logarg
      double complex eta_k34ok13(2),eta_k24ok12(2)
      double complex lik34ok13(2), logk34ok13(2)
      double complex lik24ok12(2), logk24ok12(2)


      external dilog, eta3, theta, sgn      

      m2 = 1d0                  ! > 0 ALWAYS!
      
      k12 = -m1sq/m2
      k13 = -s/m2
      k14 = -m4sq/m2
      k23 = -m2sq/m2
      k24 = -t/m2
      k34 = -m3sq/m2
      

      a = k24*k34
      b = k13*k24+k12*k34-k14*k23
      c = k12*k13
      d = k23
      

      if ((b**2-4*a*c).lt.0d0) then
         write(*,*) 'ERROR: delta less than zero!!'
         write(*,*) 'RETURN ZERO'
c use BDK version of D04m, since in this case it's REAL!!!
c         D04m_fin = D04m_fin_real(s,t,m1sq,m2sq,m3sq,m4sq,musq)
         D04m_fin = 0d0
c         write(*,*) 's,t,m1sq,m2sq,m3sq,m4sq'
c         write(*,*) s,t,m1sq,m2sq,m3sq,m4sq
         return         
c         stop
      else
         sqrtdelta = sqrt(b**2-4*a*c)
      endif


      ! the following schoolbook formula is numerically instable
      ! x(1) = (-b + sqrtdelta)/(2d0*a)
      ! x(2) = (-b - sqrtdelta)/(2d0*a)

      ! using formula from Numerical recipes, chapter 5.6
      qtemp = -0.5d0*(b + sign(sqrtdelta,b))
      if (b .lt. 0.d0) then
        x(1) = qtemp / a
        x(2) = c / qtemp
      else
        x(2) = qtemp / a
        x(1) = c / qtemp
      endif
            
      ! write(*,*) 'x(1) / x(1)_old, x(2) / x(2)_old',x(1) / ((-b + sqrtdelta)/(2d0*a)), x(2)/((-b - sqrtdelta)/(2d0*a))


      k34ok13 = k34/k13
      im_k34ok13 = (1d0-m3sq/s)/s
      
      k24ok12 = k24/k12
      im_k24ok12 = (1d0-t/m1sq)/m1sq

      im_x(1) = m2sq/sqrtdelta
      im_x(2) = -im_x(1)


      do k=1,2
         im_k34ok13_x(k) = im_k34ok13*x(k) + k34ok13*im_x(k)
         im_k24ok12_x(k) = im_k24ok12*x(k) + k24ok12*im_x(k)
         eta_k34ok13(k) = 
     -        2*ipi*eta3(-im_x(k),im_k34ok13,-im_k34ok13_x(k))
         eta_k24ok12(k) = 2*ipi*
     -        eta3(-im_x(k),im_k24ok12,-im_k24ok12_x(k))
         ln_mx(k) = log(abs(x(k))) - theta(x(k))*sgn(im_x(k))*ipi
      enddo

      lnk12 = log(abs(k12))-ipi*theta(-k12)
      lnk13 = log(abs(k13))-ipi*theta(-k13)
      lnk14 = log(abs(k14))-ipi*theta(-k14)
      lnk23 = log(abs(k23))-ipi*theta(-k23)

c      write(*,*) lnk12,lnk13,lnk14,lnk23

      do k=1,2
         arg = 1d0 + k34ok13*x(k)
         if (arg.gt.1d0) then
c     complex dilog
            logarg = log(arg)
            lik34ok13(k) = -dilog(1/arg) - 0.5d0 * logarg**2 + pi2/3
     -           + ipi * sgn(im_k34ok13_x(k)) * logarg
         else
            lik34ok13(k) = dilog(arg)
         endif
         logk34ok13(k) = log(abs(arg)) 
     -        + ipi*theta(-arg) * sgn(im_k34ok13_x(k))
c         write(*,*) lik34ok13(k),logk34ok13(k)
      enddo


      do k=1,2
         arg = 1d0 + k24ok12*x(k)
         if (arg.gt.1d0) then
c     complex dilog
            logarg = log(arg)
            lik24ok12(k) = -dilog(1/arg) - 0.5d0 * logarg**2 + pi2/3
     -           + ipi * sgn(im_k24ok12_x(k)) * logarg
         else
            lik24ok12(k) = dilog(arg)
         endif
         logk24ok12(k) = log(abs(arg)) 
     -        + ipi*theta(-arg) * sgn(im_k24ok12_x(k))
      enddo


c      write(*,*) '=========>', 
c     #     (lnk12+lnk13-lnk14-lnk23)*(-ln_mx(1)+ln_mx(2))

      total = 0d0
      do k=1,2
         total = total + (-1)**k*( - 0.5d0 *ln_mx(k)**2 
     -        - lik34ok13(k) - eta_k34ok13(k) * logk34ok13(k)
     -        - lik24ok12(k) - eta_k24ok12(k) * logk24ok12(k))
      enddo
      
      total = total + (lnk12+lnk13-lnk14-lnk23)*(-ln_mx(1)+ln_mx(2))

      D04m_fin = total/a/m2**2/(x(1)-x(2))
      
      end


********************************************************************************
********************************************************************************

c        FUNCTION D04(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4)
        FUNCTION D04(P1t,P2t,P3t,P4t,P12t,P23t,M1t,M2t,M3t,M4t)
************************************************************************
*  SCALAR 4-POINT FUNCTION WITH AT LEAST ONE MASS ZERO                 *
*  P1,P2,P3,P4 = SQUARED EXTERNAL MOMENTA			       *
*  P12 = (p1+p2)**2,  P23 = (p2+p3)**2				       *
*----------------------------------------------------------------------*
*  2.1.92  SD	         					       *
*  rearrangement to exploit massless external momenta   14.3.01  DZ    *
*  Modified: Michael Kubocz                                            *
*  Interception of NANs e.g. caused by log(0) etc. (see below)         *
*  Check, whether arguments of D04 are numerically "0" with eps1=1d-7  *
************************************************************************
        implicit none
	REAL*8 M(4),P(4,4),K(4,4)
        real*8 pi,eps,eps1
        real*8 im1,im2
        real*8 m1,m2,m3,m4
        real*8 m1t,m2t,m3t,m4t
        real*8 m02,m12,m22,m32,m42
        real*8 mm0,mm1,mm2,mm3,mm4 
        real*8 p1,p2,p3,p4,p12,p23
        real*8 p1t,p2t,p3t,p4t,p12t,p23t
        real*8 q0,q1,q2,q3,q4,q00,q12,q23
	COMPLEX*16 A1,A2,A3,A4,SWAP
	COMPLEX*16 SS(4), XX(2), X(2,4),RS(4,4)
	COMPLEX*16 S0(4),XX0(2),X0(2,4), R(4,4),G(2)
        COMPLEX*16 D04,VLI2,V_ETA,SQE,ETAS
	COMPLEX*16 AA,BB,CC,DD,IEPS,H,HH,L1,L2,L3,L4
	INTEGER I,J, i1,i2,i3,i4
        complex*16 D1,D2,D3,D4

        EXTERNAL VLI2, sqe, etas,V_ETA


        D1=DCMPLX(0d0,0d0)
        D2=DCMPLX(0d0,0d0)
        D3=DCMPLX(0d0,0d0)
        D4=DCMPLX(0d0,0d0)
        eps1=1d-7

        if(abs(P1t).le.eps1) then
           P1=0d0
        else
           P1=P1t
        endif
        if(abs(P2t).le.eps1) then
           P2=0d0
        else
           P2=P2t
        endif
        if(abs(P3t).le.eps1) then
           P3=0d0
        else
           P3=P3t
        endif
        if(abs(P4t).le.eps1) then
           P4=0d0
        else
           P4=P4t
        endif
        if(abs(P12t).le.eps1) then
           P12=0d0
        else
           P12=P12t
        endif
        if(abs(P23t).le.eps1) then
           P23=0d0
        else
           P23=P23t
        endif
        if(abs(M1t).le.eps1) then
           M1=0d0
        else
           M1=M1t
        endif
        if(abs(M2t).le.eps1) then
           M2=0d0
        else
           M2=M2t
        endif
        if(abs(M3t).le.eps1) then
           M3=0d0
        else
           M3=M3t
        endif
        if(abs(M4t).le.eps1) then
           M4=0d0
        else
           M4=M4t
        endif

        MM1=M1
        MM2=M2
        MM3=M3
        MM4=M4
        M12=M1*M1
        M22=M2*M2
        M32=M3*M3
        M42=M4*M4
        Q1=P1
        Q2=P2
        Q3=P3
	Q4=P4
        Q12=P12
        Q23=P23

C	IS AT LEAST ONE MASS ZERO ???
	IF (MM1*MM2*MM3*MM4.NE.0D0) GOTO 130

C	PERMUTATE UNTIL MM3=0D0
	GOTO 20
10	CONTINUE
	MM0=MM1
	MM1=MM2
	MM2=MM3
	MM3=MM4
	MM4=MM0
	M02=M12
	M12=M22
	M22=M32
	M32=M42
	M42=M02
	Q00=Q12
	Q12=Q23
	Q23=Q00
	Q0=Q1
	Q1=Q2
	Q2=Q3
	Q3=Q4
	Q4=Q0
20	IF (MM3.NE.0D0) GOTO 10
C	ONLY MM3 IS ZERO
	IF (MM1*MM2*MM4.NE.0D0) GOTO 30
C	ONLY MM3 AND MM4 ARE ZERO ==> 3->2, 4->3...
	IF ((MM1*MM2.NE.0D0).AND.(MM4.EQ.0D0)) GOTO 10
C	ONLY MM2 AND MM3 ARE ZERO
	IF ((MM1*MM4.NE.0D0).AND.(MM2.EQ.0D0)) GOTO 40
	WRITE(*,*)'CASE OF THIS SPECIAL D0-FUNCTION NOT IMPLEMENTED!'
	STOP

C	****** NO MASS EQUAL TO ZERO ******
130	CONTINUE
	EPS=1D-18
	IEPS=DCMPLX(0D0,EPS)
cMK check for configuration that caused problems
      if(p1.gt.0d0.and.p2.gt.0d0.and.p3.gt.0d0.and.p4.gt.0d0) then
        print*,"D0 function with all momenta light-like and non-zero massees not working"
        stop
      endif
c check for massless external momentum: excellent candidate for p13,
c leading to r13 >=1 and real.
        if (abs(q1*q2*q3*q4).le.eps1) then
           if (abs(q2).le.eps1) then
              I1 = 2
              I2 = 3
              I3 = 1
              I4 = 4
           elseif (abs(q1).le.eps1) then
              I1 = 1
              I2 = 3
              I3 = 2
              I4 = 4
           elseif (abs(q3).le.eps1) then
              I1 = 2
              I2 = 4
              I3 = 1
              I4 = 3
           else
              I1 = 1
              I2 = 4
              I3 = 2
              I4 = 3
           endif
           M(i1)=MM1
           M(i2)=MM2
           M(i3)=MM3
           M(i4)=MM4
           P(i1,i2)=Q1
           P(i3,i2)=Q2
           P(i3,i4)=Q3
           P(i1,i4)=Q4
           P(i1,i3)=Q12
           P(i3,i1)=Q12
           P(i2,i4)=Q23
           P(i4,i2)=Q23
	ELSEIF( ABS((MM1**2+MM3**2-Q12)/MM1/MM3).LT.2D0 ) THEN
C	R13 WOULD BE NOT REAL. -> PERMUTATION! -> R(2,4) IS NOT REAL.
	   M(1)=MM2
	   M(2)=MM3
	   M(3)=MM4
	   M(4)=MM1
	   P(1,2)=Q2
	   P(1,3)=Q23
	   P(1,4)=Q1
	   P(2,3)=Q3
	   P(2,4)=Q12
	   P(3,4)=Q4
	ELSE
C	R(1,3) IS REAL.
	   M(1)=MM1
	   M(2)=MM2
	   M(3)=MM3
	   M(4)=MM4
	   P(1,2)=Q1
	   P(1,3)=Q12
	   P(1,4)=Q4
	   P(2,3)=Q2
	   P(2,4)=Q23
	   P(3,4)=Q3
	ENDIF

	DO 11 J=2,4
	DO 11 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
        IF( DBLE(K(I,J)).LT.-2D0 ) THEN
c        IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
11	CONTINUE

	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	S0(1)=R(1,2)
	S0(2)=R(2,3)
	S0(3)=R(3,4)
	S0(4)=R(1,4)
	AA=K(3,4)/R(2,4)+R(1,3)*K(1,2)-K(1,4)*R(1,3)/R(2,4)-K(2,3)
	BB=(R(2,4)-1D0/R(2,4))*(R(1,3)-1D0/R(1,3))
     *		+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)/R(1,3)+R(2,4)*K(3,4)-K(1,4)*R(2,4)/R(1,3)-K(2,3)
	DD=K(2,3)-R(1,3)*K(1,2)-R(2,4)*K(3,4)+R(1,3)*R(2,4)*K(1,4)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	XX0(1)=SQE(AA,BB,CC)
	XX0(2)=CC/AA/XX0(1)
c	IF (ABS(DREAL(XX0(1)-XX(2))).LT.ABS(DREAL(XX0(1)-XX(1)))) THEN
	IF (ABS(XX0(1)-XX(2)).LT.ABS(XX0(1)-XX(1))) THEN
	  SWAP  =XX0(1)
	  XX0(1)=XX0(2)
	  XX0(2)=SWAP
	ENDIF

	DO 12 I=1,2
	G(I)  =SIGN( 1D0,DREAL(AA*(XX(I)-XX(3-I))) )
	 X(I,1)= XX(I)/R(2,4)
	X0(I,1)=XX0(I)/R(2,4)
	 X(I,2)= XX(I)/R(2,4)*R(1,3)
	X0(I,2)=XX0(I)/R(2,4)*R(1,3)
	 X(I,3)= XX(I)*R(1,3)
	X0(I,3)=XX0(I)*R(1,3)
	 X(I,4)= XX(I)
	X0(I,4)=XX0(I)
12	CONTINUE

	D04 = DCMPLX(0D0,0D0)
	DO 13 I=1,2
	DO 13 J=1,4
	A1 = 1D0+X0(I,J)*S0(J) + ABS(1D0+X0(I,J)*S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)*SS(J)))
	A2 = 1D0+X0(I,J)/S0(J) + ABS(1D0+X0(I,J)/S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)/SS(J)))
c  org. code:
c$$$           D04 = D04 + (-1D0)**(I+J)*(
c$$$     &          VLI2(A1)+V_ETA(-X(I,J),SS(J))*LOG(A1)
c$$$     &          +VLI2(A2)+V_ETA(-X(I,J),1D0/SS(J))*LOG(A2))

        if(abs(V_ETA(-X(I,J),SS(J))).ne.0d0) then
           D1=V_ETA(-X(I,J),SS(J))*LOG(A1)
        else
           D1=DCMPLX(0d0,0d0)
        endif
        if(abs(V_ETA(-X(I,J),1D0/SS(J))).ne.0d0) then
           D2=V_ETA(-X(I,J),1D0/SS(J))*LOG(A2)
        else
           D2=DCMPLX(0d0,0d0)
        endif
        D04=D04+(-1D0)**(I+J)*(VLI2(A1)+VLI2(A2)+D1+D2)
c   The enquiry avoids occurrence of NANs causing by LOG(A1) for A1=0 
c   and LOG(A2) for A2=0. At that points also V_ETA(-X(I,J),1D0/SS(J)) 
c   or V_ETA(-X(I,J),SS(J) are 0. (Michael Kubocz)

13	CONTINUE

c	IF( DIMAG(R(1,3)).EQ.0D0 ) THEN !org (makes troubles in squark pentagons)
	IF( abs(DIMAG(R(1,3))).le.eps1 ) THEN
	DO 14 I=1,2
	   A1 = (K(1,3)-2D0*R(1,3))/XX0(I)
     *		      -R(1,3)*K(1,4)+K(3,4)
     	   A2 = ((K(2,4)-2D0*R(2,4))*R(1,3)*XX0(I)
     *		      -R(2,4)*K(3,4)+K(2,3))/DD
	   A3 = (K(1,3)-2D0*R(1,3))*R(2,4)/XX0(I)
     *		      -R(1,3)*K(1,2)+K(2,3)
	   A4 = ((K(2,4)-2D0*R(2,4))*XX0(I)
     *		      -R(2,4)*K(1,4)+K(1,2))/DD
           if (A1.ne.0) then
	     L1 = LOG( A1-ABS(A1)*IEPS )
           else
             L1 = -1d99
           endif
           if (A2.ne.0) then
             L2 = LOG( A2+ABS(A2)*IEPS*G(I)*SIGN(1D0,DREAL(R(1,3))
     *                                            *DIMAG(RS(2,4))) )
           else
             L2 = -1d99
           endif
           if (A3.ne.0) then
             L3 = LOG( A3-ABS(A3)*IEPS )
           else
             L3 = -1d99
           endif
           if (A4.ne.0) then
             L4 = LOG( A4+ABS(A4)*IEPS*G(I)*SIGN(1D0,DIMAG(RS(2,4))) )
           else
             L4 = -1d99
           endif

c org. code:
c$$$	   D04 = D04 
c$$$     &         + (3D0-2D0*I)*(
c$$$     *		       ETAS( -XX(I),R(1,3),RS(1,3) )
c$$$     *		          *( LOG(R(1,3)*XX(I)) + L1 + L2 )
c$$$     *		     + ETAS( -XX(I),1D0/R(2,4),1D0/RS(2,4) )
c$$$     *		          *( LOG(XX(I)/R(2,4)) + L3 + L4 )
c$$$     *		     - ( ETAS( -XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4) )
c$$$     *		       + V_ETA( RS(1,3),1D0/RS(2,4) )                  )
c$$$     *		        *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 + L2 )
c$$$     *	  	     + V_ETA( RS(1,3),1D0/RS(2,4) )
c$$$     *		       *ETAS(-XX(I),-R(1,3)/R(2,4),-RS(1,3)/RS(2,4))   )

            if(abs(ETAS(-XX(I),R(1,3),RS(1,3))).ne.0d0) then
               D1=ETAS(-XX(I),R(1,3),RS(1,3))*(LOG(R(1,3)*XX(I))+L1+L2)
            else
               D1=DCMPLX(0d0,0d0)
            endif
            if(abs(ETAS(-XX(I),1D0/R(2,4),1D0/RS(2,4))).ne.0d0) then
               D2=ETAS(-XX(I),1D0/R(2,4),1D0/RS(2,4))*(LOG(XX(I)/R(2,4))
     &           +L3+L4)
            else
               D2=DCMPLX(0d0,0d0)
            endif
            if((abs(ETAS(-XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4))).ne.0d0)
     &          .or.(abs(V_ETA(RS(1,3),1D0/RS(2,4))).ne.0d0)) then 
               D3=-(ETAS(-XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4))
     &            +V_ETA(RS(1,3),1D0/RS(2,4)))*(LOG(XX(I)*R(1,3)/R(2,4))+L3+L2)
            else
               D3=DCMPLX(0d0,0d0)
            endif
            D4=V_ETA(RS(1,3),1D0/RS(2,4))*ETAS(-XX(I),-R(1,3)/R(2,4),-RS(1,3)/RS(2,4))
            D04=D04+(3D0-2D0*I)*(D1+D2+D3+D4)
c   The enquiry avoids occurrence of NANs causing by LOG(0). At that points 
c   also V_ETA(...) and ETAS(...) are 0. (Michael Kubocz)

14	CONTINUE
	ELSE
	DO 15 I=1,2
	   L1 = LOG( R(2,4)/XX0(I)+XX0(I)/R(2,4)+K(1,2)
     *		     -XX0(I)/R(2,4)*EPS*BB*G(I) )
	   L2 = LOG( R(1,3)*XX0(I)+1D0/XX0(I)/R(1,3)+K(3,4)
     *		     -XX0(I)*R(1,3)*EPS*BB*G(I) )
	   L3 = LOG( R(1,3)/R(2,4)*XX0(I)+R(2,4)/XX0(I)/R(1,3)+K(2,3)
     *		     -XX0(I)*R(1,3)/R(2,4)*EPS*BB*G(I) )

	   D04 = D04 
     &          + (3D0-2D0*I)*(
     *		     V_ETA(-XX(I),1D0/R(2,4))
     *		      *( LOG(XX(I)/R(2,4)) + L1 )
     *		    +V_ETA(-XX(I),R(1,3))
     *		      *( LOG(R(1,3)*XX(I)) + L2 )
     *		    -( V_ETA(-XX(I),R(1,3)/R(2,4))
     *		      +V_ETA(R(1,3),1D0/R(2,4)) )
     *		      *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 )
     *	  	    +V_ETA(R(1,3),1D0/R(2,4))
     *		      *V_ETA(-XX(I),-R(1,3)/R(2,4))
     *		      *(1D0-G(I)*SIGN(1D0,DREAL(BB)))	    )
15	CONTINUE
	ENDIF

	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN


C--->	***************** SPEZIELL ( --> T.SACK-PROMOTION )
C	D1=Q12-M12
C	D2=Q2 -M22
C	D3=Q3 -M42
C	IF ((D1*D2.LE.0D0).OR.(D2*D3.LE.0D0)) THEN
C	   WRITE(*,*) 'THE CASE OF DIFFERENT SIGNS OF THE D1,D2,D3'
C	   WRITE(*,*) 'IN D04(...) IS NOT IMPLEMENTED !!!'
C	   STOP
C	ENDIF
C	NM1=ABS(MM1/D1)
C	NM2=ABS(MM2/D2)
C	NM3=ABS(MM4/D3)
C	NP1=Q2/D2**2+Q12/D1**2+(Q1-Q2-Q12)/D1/D2
C	NP2=Q2/D2**2+ Q3/D3**2+(Q23-Q2-Q3)/D2/D3
C	NP3=Q3/D3**2+Q12/D1**2+(Q4-Q3-Q12)/D1/D3
C	D04=C04(NP1,NP2,NP3,NM1,NM2,NM3)/D1/D2/D3

C	*************** ALLGEMEIN


C	****** ONLY MM3 IS ZERO ******
30	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)
	M(1)=MM1
	M(2)=MM2
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 1 J=2,4
	DO 1 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
1	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(3,4)/R(2,4)-K(2,3)
	BB=K(1,3)*(1D0/R(2,4)-R(2,4))+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(1,3)*K(1,4)*R(2,4)+R(2,4)*K(3,4)-K(2,3)
	DD=K(2,3)-R(2,4)*K(3,4)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 2 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
2	CONTINUE
	D04 = DCMPLX(0D0,0D0)
	DO 3 I=1,2
	D04 = D04 + (2D0*I-3D0)*(
     *		VLI2(1D0+SS(4)*X(I,4))
     *	       -VLI2(1D0+SS(1)*X(I,1))
     *	       +VLI2(1D0+X(I,4)/SS(4))
     *	       -VLI2(1D0+X(I,1)/SS(1))
     *	       +V_ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       -V_ETA(-X(I,1),SS(1))*LOG(1D0+SS(1)*X(I,1))
     *	       +V_ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -V_ETA(-X(I,1),1D0/SS(1))*LOG(1D0+X(I,1)/SS(1))
     *	       -VLI2(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +VLI2(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	       -V_ETA(-X(I,4),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +V_ETA(-X(I,1),(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))   )
	IF (DIMAG(R(2,4)).NE.0D0) THEN
	   H=V_ETA(-1D0/XX(I),R(2,4))
	ELSE
	   H=DCMPLX(0D0,0D0)
	   IF (DREAL(R(2,4)).LT.0D0) THEN
	      HH=-1D0/XX(I)
	      IM1=DIMAG(HH)
	      IM2=DIMAG(RS(2,4))
              pi = 4.D0*datan(1.D0)
	      IF ((IM1.GT.0D0).AND.(IM2.GT.0D0)) THEN
	         H=-DCMPLX(0D0,2D0*PI)
	      ENDIF
	      IF ((IM1.LT.0D0).AND.(IM2.LT.0D0)) THEN
	         H=+DCMPLX(0D0,2D0*PI)
	      ENDIF
	   ENDIF
	ENDIF
	D04 = D04 + (2D0*I-3D0)*
     *	          H*( LOG( (K(1,2)-R(2,4)*K(1,4)
     *			  +XX(I)*(1D0/R(2,4)-R(2,4)))/DD )
     *		     +LOG(K(1,3)-IEPS) )
3	CONTINUE
	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN

C	****** ONLY MM2 AND MM3 ARE ZERO ******
40	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)

	M(1)=MM1
	M(2)=10D0
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 4 J=2,4
	DO 4 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.2) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.2) K(I,J)=K(I,J)-M(J)/M(I)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
4	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(2,4)*K(3,4)-K(2,3)
	BB=K(1,3)*K(2,4)+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(2,3)
	DD=K(2,3)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 5 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
5	CONTINUE
	D04 = DCMPLX(0D0,0D0)
	DO 6 I=1,2
	D04 = D04 + (2D0*I-3D0)*(
     *		VLI2(1D0+SS(4)*X(I,4))
     *	       +VLI2(1D0+X(I,4)/SS(4))
     *	       +V_ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       +V_ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -VLI2(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -VLI2(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       -V_ETA(-XX(I),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -V_ETA(-XX(I),(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       +LOG(-XX(I))*( LOG(K(1,2)-IEPS)
     *			     +LOG(K(1,3)-IEPS)-LOG(K(2,3)-IEPS) ) )
6	CONTINUE
	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))

	RETURN

	END

********************************************************************************


************************************************************************
