      function D0finG_QUAD(m0,m1,m2,m3,s12,s23,p1,p2,p3,p4,musq)
      IMPLICIT NONE
      complex*32 D0finG_QUAD
      real*16 m0,m1,m2,m3,p1,p2,p3,p4,s12,s23,musq

      complex*32 D04QUAD1 !D0fin!,D0t1,D0t2,D0t3,D04
      External D04QUAD1   !D0fin!,D0t1,D0t2,D0t3,D04
c      complex*32 rslt(0:2)
      real*16 eps
      parameter (eps=1d-6)
      logical UseFF
      real*16 xpi(13)
      integer ier
c      parameter (UseAVH=.false.)
      parameter (UseFF=.true.)
      
      real*8 m0_D,m1_D,m2_d,m3_D,s12_D,s23_D,p1_D,p2_D,p3_D,p4_D,musq_D

      external D04
      complex*16 D04, D0finG_dble

c MISSING THIS ONE!
c      complex*16 function D0t1(s1,s2,s3,s4,s,t,mi1s,
c     delta is the mass regulator 
c     1.) mi4s <>0, s3 <=> 0, s4 <=> 0, s <>0, t<>0
c     mi4 must be nonzero 
c-------------------------
c      complex*16 function D0t2(s1,s2,s3,s4,s,t,mi1s,
c     $     mi2s,mi3s,mi4s,mus,id)
c     delta is the mass regulator 
c     2.) mi4s <>0, s1 = 0, s2 <> 0 s3 <=> 0, s4 <> 0, s <>0, t<>0
c     mi4 must be nonzero mi1=mi2=mi3=0
c     mi3 = mi1 = mi2 = 0 
C--------------------------------
c      complex*16 function D0t3(s1,s2,s3,s4,s,t,mi1s,
c     $     mi2s,mi3s,mi4s,mus,id)
c     3.) mi4s <>0, mi3s <> 0, s1 = 0 s2 <> 0 s3 <> 0, s4 <> 0, s <>0, t<>0
c     mi4 must be nonzero mi1=mi2=0
c     mi3 must be nonzero

cfc      if(UseAVH) then
cfc#ifdef WITH_ONELOOP
cfc         m0s=m0*m0
cfc         m1s=m1*m1
cfc         m2s=m2*m2
cfc         m3s=m3*m3
cfc         call avh_olo_mu_set(sqrt(musq))
cfc         call avh_olo_onshell(eps)
cfc         call avh_olo_d0m(rslt,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s)
cfc         D0finG = rslt(0) 
cfc#else
cfc         call avh_error()
cfc#endif
cfc      else 
        if (UseFF) then
c take out mass combination which FF cannot yet handle
        if ( (abs(m0).lt.eps) .and.
     &       (abs(m1).lt.eps) .and.
     &       (abs(m2).lt.eps) .and.
     &       (abs(m3).lt.eps) ) then
c           print*, "D0_masses0"
           call D0fin_QUAD(s12,s23,p1,p2,p3,p4,musq,D0finG_QUAD)
         else if ( (m0.eq.m1) .and.
     &             (m0.eq.m2) .and.
     &             (m0.eq.m3) ) then
           D0finG_QUAD=D04QUAD1(p1,p2,p3,p4,s12,s23,m0,m1,m2,m3)
c           print*, "D02_equalmasses",D0finG_QUAD

           m0_D=   m0
           m1_D=   m1 
           m2_d=   m2 
           m3_D=   m3 
           s12_D=  s12
           s23_D=  s23
           p1_D=   p1 
           p2_D=   p2 
           p3_D=   p3 
           p4_D=   p4 
           musq_D= musq
         D0finG_dble=D04(p1_d,p2_d,p3_d,p4_d,s12_d,s23_d
     &             ,m0_d,m1_d,m2_d,m3_d)
         if(abs(D0finG_dble/D0fing_QUAD-1d0).gt.1d-7) then
            print*,"Attention in D0"
            print*,"dble",D0finG_dble
            print*,"quad",D0finG_quad
            print*, "ratio",abs(D0finG_dble/D0fing_QUAD-1d0)
         endif
c           print*, "D02_equalmasses",D0finG_dble
         else
c           print*, "D0_using ff"
           xpi( 1) = m0**2
           xpi( 2) = m1**2
           xpi( 3) = m2**2
           xpi( 4) = m3**2
           xpi( 5) = p1
           xpi( 6) = p2
           xpi( 7) = p3
           xpi( 8) = p4
           xpi( 9) = s12
           xpi(10) = s23
           ier=0
           call vbfquadffxd0(D0finG_QUAD,xpi,ier)
c           print*, "D0_using ff",D0finG_QUAD
cfc           xpi_d( 1) = m0_d**2
cfc           xpi_d( 2) = m1_d**2
cfc           xpi_d( 3) = m2_d**2
cfc           xpi_d( 4) = m3_d**2
cfc           xpi_d( 5) = p1_d
cfc           xpi_d( 6) = p2_d
cfc           xpi_d( 7) = p3_d
cfc           xpi_d( 8) = p4_d
cfc           xpi_d( 9) = s12_d
cfc           xpi_d(10) = s23_d
cfc           ier=0
c           call vbfffxd0(D0finG_dble,xpi_d,ier)
c           print*, "D0_using ff",D0finG_dble
c           D0finG_QUAD=D0finG_dble
c           pause

         endif
        else

        if ( (abs(m0).lt.eps) .and.
     &       (abs(m1).lt.eps) .and.
     &       (abs(m2).lt.eps) .and.
     &       (abs(m3).lt.eps) ) then
           call D0fin_QUAD(s12,s23,p1,p2,p3,p4,musq,D0finG_QUAD)
        else if ( (m0.eq.m1) .and.
     &             (m0.eq.m2) .and.
     &             (m0.eq.m3) ) then
           D0finG_QUAD=D04quad1(p1,p2,p3,p4,s12,s23,m0,m1,m2,m3)
         else
            print*,"D:Wrong"
            stop
        endif    
        return
        endif
        return
        end function


c      function D0fin_QUAD(s,t,p1sq,p2sq,p3sq,p4sq,musq)
       subroutine D0fin_QUAD(s,t,p1sq,p2sq,p3sq,p4sq,musq,D0_fin_QUAD)
       implicit none
       complex*32 D0_fin_QUAD
       real*16 s,t,musq,p1sq,p2sq,p3sq,p4sq
c      Complex*32 D01m_fin_QUAD,D02m_fin_QUAD,D02m_fin_QUADe,D03m_fin_QUAD
c     -   ,D04m_fin_QUAD,D00m_fin_QUAD
c      External D01m_fin_QUAD,D02m_fin_QUAD,D02m_fin_QUADe,D03m_fin_QUAD
c     -   ,D04m_fin_QUAD,D00m_fin_QUAD
       real*16 tiny
       parameter (tiny=1q-13)
  
c$$$      print*, s,t,p1sq,p2sq,p3sq,p4sq,musq

      if(abs(p1sq).lt.tiny) then
          if(abs(p2sq).lt.tiny)then
               if(abs(p3sq).lt.tiny)then
                   if(abs(p4sq).lt.tiny) then
c  D(s,t,0,0,0,0)                                                      
                             call D00mfin_QUAD(s,t,musq,D0_fin_QUAD)
                             RETURN
                   else
c D(s,t,0,0,0,p4sq) 
                            call D01mfin_QUAD(s,t,p4sq,musq,D0_fin_QUAD)
                              RETURN
                   endif 
               else 
                   if(abs(p4sq).lt.tiny) then
c D(s,t,0,0,p3sq,0)
                            call D01mfin_QUAD(s,t,p3sq,musq,D0_fin_QUAD)
                             RETURN 
                   else
c D(s,t,0,0,p3sq,p4sq)
                       call D02mfin_QUAD(s,t,p3sq,p4sq,musq,D0_fin_QUAD)
                             RETURN
                   endif
               endif
          else
                     if(abs(p3sq).lt.tiny)then
                           if(abs(p4sq).lt.tiny) then
c D(s,t,0,p2sq,0,0)
                            call D01mfin_QUAD(s,t,p2sq,musq,D0_fin_QUAD)
                             RETURN
                           else
c D(s,t,0,p2sq,0,p4sq)
                     call D02mfin_QUADe(s,t,p2sq,p4sq,musq,D0_fin_QUAD) 
                             RETURN
                           endif 
                     else 
                            if(abs(p4sq).lt.tiny) then
c D(s,t,0,p2sq,p3sq,0)
                       call D02mfin_QUAD(t,s,p2sq,p3sq,musq,D0_fin_QUAD) 
                             RETURN
                             else 
c D(s,t,0,p2sq,p3sq,p4sq)
                  call D03mfin_QUAD_BDK(s,t,p2sq,p3sq,p4sq,musq,D0_fin_QUAD)
c                  call D03mfin_QUAD(s,t,p2sq,p3sq,p4sq,musq,D0_fin_QUAD) 
c                              write(*,*)" D0_fin_QUAD", D0_fin_QUAD
                             endif
                      endif
           endif 
      else
           if(abs(p2sq).lt.tiny)then
               if(abs(p3sq).lt.tiny)then
                   if(abs(p4sq).lt.tiny) then
c  D(s,t,p1sq,0,0,0) 
                            call D01mfin_QUAD(s,t,p1sq,musq,D0_fin_QUAD)
                             RETURN
                   else
c  D(s,t,p1sq,0,0,p4sq)
                      call D02mfin_QUAD(t,s,p4sq,p1sq,musq,D0_fin_QUAD) 
                             RETURN
                   endif 
               else 
                   if(abs(p4sq).lt.tiny) then
c  D(s,t,p1sq,0,p3sq,0) 
                    call D02mfin_QUADe(s,t,p1sq,p3sq,musq,D0_fin_QUAD) 
                             RETURN
                   else
c  D(s,t,p1sq,0,p3sq,p4sq) 
              call D03mfin_QUAD_BDK(t,s,p3sq,p4sq,p1sq,musq,D0_fin_QUAD)
c                   call D03mfin_QUAD(t,s,p3sq,p4sq,p1sq,musq,D0_fin_QUAD)
                            RETURN
                   endif
               endif
          else
               if(abs(p3sq).lt.tiny)then
                   if(abs(p4sq).lt.tiny) then
c  D(s,t,p1sq,p2sq,0,0) 
                      call D02mfin_QUAD(s,t,p1sq,p2sq,musq,D0_fin_QUAD)
                             RETURN
                   else
c  D(s,t,p1sq,p2sq,0,p4sq) 
              call D03mfin_QUAD_BDK(s,t,p4sq,p1sq,p2sq,musq,D0_fin_QUAD)
c                  call D03mfin_QUAD(s,t,p4sq,p1sq,p2sq,musq,D0_fin_QUAD) 
                             RETURN
                   endif 
               else 
                   if(abs(p4sq).lt.tiny) then
c  D(s,t,p1sq,p2sq,p3sq,0) 
              call D03mfin_QUAD_BDK(t,s,p1sq,p2sq,p3sq,musq,D0_fin_QUAD)
c              call D03mfin_QUAD(t,s,p1sq,p2sq,p3sq,musq,D0_fin_QUAD) 
                             RETURN
                   else
c  D(s,t,p1sq,p2sq,p3sq,p4sq) 
             call D04mfin_QUAD(s,t,p1sq,p2sq,p3sq,p4sq,musq,D0_fin_QUAD)
                             RETURN
                   endif
               endif
           endif
         endif
         end subroutine


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCC                SCALAR INTEGRALS                CCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c    ***********   D00m_fin_QUAD(s,t,musq)   *************
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
c          N_ep * D00m_fin_QUAD(s,t,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

c      complex*32 function D00m_fin_QUAD(s,t,musq)
      subroutine  D00mfin_QUAD(s,t,musq,D00m_fin_QUAD)
      implicit none
      complex*32 D00m_fin_QUAD
      real*16  s,t,musq
      real*16  ms,mt
      complex*32 lnms,lnmt,lnsot


      complex*32 ipi
      parameter (ipi=(0q0,3.14159265358979323846264338327950288q0))
      real*16  pi,pi2, pi2o3t5
      parameter (pi=3.141 592 653 589 793 238 462 643 383 279 502 88q0,
     #     pi2 =9.86960440108935861883449099987615114q0,
     #  pi2o3t5=16.4493406684822643647241516664602519q0 )

      complex*32 ris
      real*16  prefactor, theta_QUAD, dilog_QUAD, im_part_QUAD
      External dilog_QUAD   
      if (musq.lt.0q0) then
         write(*,*) 
     #'POSSIBLE ERROR IN D00m_fin_QUAD: SCALE MUSQ LESS THAN ZERO!!'
      endif

      prefactor = 1q0/(s*t)
      ms = -s/musq
      mt = -t/musq
      
      lnms = Log(abs(ms)) - ipi*theta_QUAD(-ms)
      lnmt = Log(abs(mt)) - ipi*theta_QUAD(-mt)
      lnsot = Log(abs(s/t)) + ipi * im_part_QUAD(s,t)

c      ris = lnms**2+lnmt**2-lnsot**2-pi2o3t5
      ris = lnms*lnms+lnmt*lnmt-lnsot*lnsot-pi2o3t5

      D00m_fin_QUAD = prefactor * ris
      return
      end


c    ***********   D01m_fin_QUAD(s,t,m4sq,musq)   *************
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
c          N_ep * D01m_fin_QUAD(s,t,m3sq,m4sq,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

c      complex*32 function D01m_fin_QUAD(s,t,m4sq,musq)
      subroutine D01mfin_QUAD(s,t,m4sq,musq,D01m_fin_QUAD)
      implicit none
      complex*32 D01m_fin_QUAD
      real*16  s,t,m4sq,musq
      real*16  ms,mt,mm4sq
      complex*32 lnms,lnmt,lnmm4sq,lnsot
      complex*32 li2arg1,li2arg2
      real*16  arg1,arg2
      complex*32 ipi
      parameter (ipi=(0q0,3.14159265358979323846264338327950288q0))
      real*16  pi,pi2, pi2o3t2,pi2o3
      parameter (pi=3.14159265358979323846264338327950288q0,
     #     pi2 =9.86960440108935861883449099987615114q0,
     #  pi2o3  =3.28986813369645287294483033329205038q0 ,
     #  pi2o3t2=6.57973626739290574588966066658410076q0 )

      complex*32 ris
      real*16  prefactor, theta_QUAD, dilog_QUAD, im_part_QUAD
       real*16 logarg1,logarg2
      External dilog_QUAD
      if (musq.lt.0q0) then
         write(*,*) 
     # 'POSSIBLE ERROR IN D01m_fin_QUAD: SCALE MUSQ LESS THAN ZERO!!'
      endif

      prefactor = 1q0/(s*t)
      ms = -s/musq
      mt = -t/musq
      mm4sq = - m4sq/musq
      
      lnms = Log(abs(ms)) - ipi*theta_QUAD(-ms)
      lnmt = Log(abs(mt)) - ipi*theta_QUAD(-mt)
      lnmm4sq = Log(abs(mm4sq)) - ipi*theta_QUAD(-mm4sq)
      lnsot = Log(abs(s/t)) + ipi * im_part_QUAD(s,t)

      arg1 = 1q0-m4sq/s
      arg2 = 1q0-m4sq/t

c     (m4sq/s.lt.0d0)
      if (arg1.gt.1q0) then
         logarg1= Log(arg1)
         li2arg1 = -dilog_QUAD(1q0/arg1) - logarg1*logarg1/2q0+pi2o3 
     #        - ipi*logarg1*im_part_QUAD(m4sq,s)
      else
         li2arg1 = dilog_QUAD(arg1)
      endif

c     (m4sq/t.lt.0d0)
      if (arg2.gt.1q0)  then
         logarg2=Log(arg2)
         li2arg2 = -dilog_QUAD(1q0/arg2) - logarg2*logarg2/2q0+pi2o3 
     #        - ipi*logarg2*im_part_QUAD(m4sq,t)
      else
         li2arg2 = dilog_QUAD(arg2)
      endif

c      ris = -lnmm4sq**2+lnms**2+lnmt**2-lnsot**2-2d0*(li2arg1+li2arg2)
c     #   -pi2o3t2

      ris = -lnmm4sq*lnmm4sq+lnms*lnms+lnmt*lnmt-lnsot*lnsot
     # -2q0*(li2arg1+li2arg2)-pi2o3t2

      D01m_fin_QUAD = prefactor * ris
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCC                SCALAR INTEGRALS                CCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c    ***********   D02m_fin_QUAD(s,t,m3sq,m4sq,musq)   *************
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
c          N_ep * D02m_fin_QUAD(s,t,m3sq,m4sq,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

c      complex*32 function D02m_fin_QUAD(s,t,m3sq,m4sq,musq)
       subroutine D02mfin_QUAD(s,t,m3sq,m4sq,musq,D02m_fin_QUAD)
      implicit none
      complex*32 D02m_fin_QUAD
      real*16  s,t,m3sq,m4sq,musq
      real*16  ms,mt,mm3sq,mm4sq
      complex*32 lnms,lnmt,lnmm3sq,lnmm4sq,lnsot
      complex*32 li2arg1,li2arg2
      real*16  arg1,arg2
      complex*32 ipi
      parameter (ipi=(0q0,3.14159265358979323846264338327950288q0))
      real*16  pi,pi2,pi2o3,pi2o6
      parameter (pi=3.14159265358979323846264338327950288q0,
     #  pi2o3  =3.28986813369645287294483033329205038q0 ,
     #  pi2o6  =1.64493406684822643647241516664602519q0,
     #     pi2 =9.86960440108935861883449099987615114q0)

      complex*32 ris
      real*16  prefactor, theta_QUAD, dilog_QUAD, im_part_QUAD
      real*16 logarg1,logarg2    
      External dilog_QUAD
      if (musq.lt.0q0) then
         write(*,*) 
     #'POSSIBLE ERROR IN D02m_fin_QUAD: SCALE MUSQ LESS THAN ZERO!!'
      endif

      prefactor = 1q0/(s*t)
      ms = -s/musq
      mt = -t/musq
      mm3sq = - m3sq/musq
      mm4sq = - m4sq/musq
      
      lnms = Log(abs(ms)) - ipi*theta_QUAD(-ms)
      lnmt = Log(abs(mt)) - ipi*theta_QUAD(-mt)
      lnmm3sq = Log(abs(mm3sq)) - ipi*theta_QUAD(-mm3sq)
      lnmm4sq = Log(abs(mm4sq)) - ipi*theta_QUAD(-mm4sq)
      lnsot = Log(abs(s/t)) + ipi * im_part_QUAD(s,t)

      arg1 = 1q0-m3sq/t
      arg2 = 1q0-m4sq/t

c     (m3sq/t.lt.0d0)
      if (arg1.gt.1q0) then
         logarg1=Log(arg1)
         li2arg1 = -dilog_QUAD(1q0/arg1) - logarg1*logarg1/2q0+pi2o3 
     #        - ipi*logarg1*im_part_QUAD(m3sq,t)
      else
         li2arg1 = dilog_QUAD(arg1)
      endif

c     (m4sq/t.lt.0d0)
      if (arg2.gt.1q0)  then
         logarg2=Log(arg2)
         li2arg2 = -dilog_QUAD(1q0/arg2) - logarg2*logarg2/2q0+pi2o3 
     #        - ipi*logarg2*im_part_QUAD(m4sq,t)
      else
         li2arg2 = dilog_QUAD(arg2)
      endif

c      ris = lnms**2/2+lnmt**2-lnmm3sq**2/2-lnmm4sq**2/2
c     #     +lnmm3sq*lnmm4sq-(lnmm4sq+lnmm3sq-lnms)*lnms
c     #     -2*li2arg1-2*li2arg2-lnsot**2 - pi2/6
    
      ris = (lnms*lnms-lnmm3sq*lnmm3sq-lnmm4sq*lnmm4sq)/2q0+lnmt*lnmt
     #     +lnmm3sq*lnmm4sq-(lnmm4sq+lnmm3sq-lnms)*lnms
     #     -2*(li2arg1+li2arg2)-lnsot*lnsot - pi2o6

    
      D02m_fin_QUAD = prefactor * ris
      return
      end


c    ***********   D02m_fin_QUADe(s,t,m1sq,m3sq,musq)   *************
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
c          N_ep * D02m_fin_QUAD3(s,t,m1sq,m3sq,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

c      complex*32 function D02m_fin_QUADe(s,t,m1sq,m3sq,musq)
      subroutine D02mfin_QUADe(s,t,m1sq,m3sq,musq,D02m_fin_QUADe)
      implicit none
      complex*32 D02m_fin_QUADe
      real*16  s,t,m1sq,m3sq,musq
      real*16  ms,mt,mm1sq,mm3sq
      complex*32 lnms,lnmt,lnmm1sq,lnmm3sq,lnsot,lnarg5
      complex*32 li2arg1,li2arg2,li2arg3,li2arg4,li2arg5
      real*16  arg1,arg2,arg3,arg4,arg5
      complex*32 ipi
      parameter (ipi=(0q0,3.14159265358979323846264338327950288q0))
      real*16  pi,pi2,pi2o3
      parameter (pi=3.14159265358979323846264338327950288q0,
     #  pi2o3  =3.28986813369645287294483033329205038q0 ,
     #     pi2 =9.86960440108935861883449099987615114q0 )

      complex*32 ris
      real*16  prefactor, theta_QUAD, dilog_QUAD, im_part_QUAD, sgn_QUAD
      real*16  im_a, im_b,im_ab,  eta3_QUAD
      real*16 logarg1,logarg2,logarg3,logarg4,logarg5             
      External dilog_QUAD
      if (musq.lt.0q0) then
         write(*,*) 
     #'POSSIBLE ERROR IN D02m_fin_QUADe: SCALE MUSQ LESS THAN ZERO!!'
      endif

      prefactor = 1q0/(m1sq*m3sq- s*t)
      ms = -s/musq
      mt = -t/musq
      mm1sq = - m1sq/musq
      mm3sq = - m3sq/musq
      
      lnms = Log(abs(ms)) - ipi*theta_QUAD(-ms)
      lnmt = Log(abs(mt)) - ipi*theta_QUAD(-mt)
      lnmm1sq = Log(abs(mm1sq)) - ipi*theta_QUAD(-mm1sq)
      lnmm3sq = Log(abs(mm3sq)) - ipi*theta_QUAD(-mm3sq)
      lnsot = Log(abs(s/t)) + ipi*im_part_QUAD(s,t)
c      write(*,*) "lnsot", lnsot
 
      arg1 = 1q0-m1sq/t
      arg2 = 1q0-m3sq/t
      arg3 = 1q0-m1sq/s
      arg4 = 1q0-m3sq/s
      arg5=  1q0-(m1sq*m3sq)/(s*t)

c$$$      write(*,*) "arg1", arg1
c$$$      write(*,*) "arg2", arg2
c$$$      write(*,*) "arg3", arg3
c$$$      write(*,*) "arg4", arg4
c$$$      write(*,*) "arg5", arg5

c     (m1sq/t.lt.0d0)
      if (arg1.gt.1q0) then
         logarg1=Log(arg1)
         li2arg1 = -dilog_QUAD(1q0/arg1) - logarg1*logarg1/2q0+pi2o3 
     #        - ipi*logarg1*im_part_QUAD(m1sq,t)
      else
         li2arg1 = dilog_QUAD(arg1)
      endif


c     (m3sq/t.lt.0d0)
      if (arg2.gt.1q0)  then
         logarg2=Log(arg2)
         li2arg2 = -dilog_QUAD(1q0/arg2) - logarg2*logarg2/2q0+pi2o3 
     #        - ipi*logarg2*im_part_QUAD(m3sq,t)
      else
         li2arg2 = dilog_QUAD(arg2)
      endif

c     (m1sq/s.lt.0d0)
      if (arg3.gt.1q0)  then
         logarg3=Log(arg3)
         li2arg3 = -dilog_QUAD(1q0/arg3) - logarg3*logarg3/2q0+pi2o3 
     #        - ipi*logarg3*im_part_QUAD(m1sq,s)
      else
         li2arg3 = dilog_QUAD(arg3)
      endif
c     (m3sq/s.lt.0d0)
      if (arg4.gt.1q0)  then
         logarg4=Log(arg4)
         li2arg4 = -dilog_QUAD(1q0/arg4) - logarg4*logarg4/2q0+pi2o3 
     #        - ipi*logarg4*im_part_QUAD(m3sq,s)
      else
         li2arg4 = dilog_QUAD(arg4)
      endif

      im_a = m1sq/s*(1q0/m1sq-1q0/s)
      im_b = m3sq/t*(1q0/m3sq-1q0/t)
      im_ab =m1sq*m3sq/(s*t)*(1q0/m1sq-1q0/s+1q0/m3sq-1q0/t)
       

c     m1sq*m3sq/(s*t)  < 0
      if (arg5.gt.1q0) then
         logarg5=Log(arg5)
         li2arg5 = -dilog_QUAD(1q0/arg5) - logarg5* logarg5/2q0+pi2o3 
     #        - ipi*logarg5*sgn_QUAD(im_ab)
      else
         li2arg5 = dilog_QUAD(arg5)
      endif

      if (arg5.lt.0d0) then
         lnarg5 = log(abs(arg5)) - ipi*sgn_quad(im_ab)
      else
         lnarg5 = log(arg5)
      endif
      li2arg5 = li2arg5  + 2*ipi*eta3_quad(im_a,im_b,im_ab)*lnarg5
c      ris = lnmm1sq**2+lnmm3sq**2-lnms**2-lnmt**2+lnsot**2+
c     #  2.d0*(li2arg1+li2arg2+li2arg3+li2arg4-li2arg5)   

      ris=lnmm1sq*lnmm1sq+lnmm3sq*lnmm3sq-lnms*lnms-lnmt*lnmt
     #  +lnsot*lnsot+2.q0*(li2arg1+li2arg2+li2arg3+li2arg4-li2arg5)   

c      write(*,*) "ris", ris

      D02m_fin_QUADe = prefactor*ris
      return
      end


c    ***********   D03m_fin_QUAD(s,t,m2sq,m3sq,m4sq,musq)   *************
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
c          N_ep * D03m_fin_QUAD(s,t,m2sq,m3sq,m4sq,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

c      function D03m_fin_QUAD(s,t,m2sq,m3sq,m4sq,musq)
      subroutine D03mfin_QUAD(s,t,m2sq,m3sq,m4sq,musq,D03m_fin_QUAD)
      implicit none
      complex*32 D03m_fin_QUAD
      real*16  s,t,m2sq,m3sq,m4sq,musq
c      complex*32 D03m_fin_QUAD_BDK,D03m_fin_QUAD_DN
       call D03mfin_QUAD_BDK(s,t,m2sq,m3sq,m4sq,musq,D03m_fin_QUAD)
c      D03m_fin_QUAD = D03m_fin_QUAD_BDK(s,t,m2sq,m3sq,m4sq,musq)
c      D03m_fin_QUAD = D03m_fin_QUAD_DN(s,t,m2sq,m3sq,m4sq,musq)
c      write(*,*)" here 2 D03m_fin_QUAD",D03m_fin_QUAD
      return
      end


c      function D03m_fin_QUAD_BDK(s,t,m2sq,m3sq,m4sq,musq)
      subroutine D03mfin_QUAD_BDK(s,t,m2sq,m3sq,m4sq,musq,D03m_fin_QUAD_BDK)
      implicit none
      complex*32 D03m_fin_QUAD_BDK
c      complex*16 D03m_fin_QUAD_BDK_R
      real*16  s,t,m2sq,m3sq,m4sq,musq
      real*16  ms,mt,mm2sq,mm3sq,mm4sq,sot
      complex*32 lnms,lnmt,lnmm2sq,lnmm3sq,lnmm4sq,lnsot,lnarg3
      complex*32 li2arg1,li2arg2,li2arg3
      real*16  arg1,arg2,arg3
      complex*32 ipi
      parameter (ipi=(0q0,3.14159265358979323846264338327950288q0))
      real*16  pi,pi2,pi2o3
      parameter (pi=3.14159265358979323846264338327950288q0,
     #  pi2o3  = 3.28986813369645287294483033329205038q0 ,
     #     pi2 = 9.86960440108935861883449099987615114q0)
      complex*32 ris
      real*16  theta_QUAD, dilog_QUAD, sgn_QUAD, im_part_QUAD
      real*16  im_a, im_b,im_ab,  eta3_QUAD, prefactor
      real*16 logarg1,logarg2,logarg3  
      External dilog_QUAD
      if (musq.lt.0q0) then
         write(*,*) 
     #'POSSIBLE ERROR IN D03m_fin_QUAD: SCALE MUSQ LESS THAN ZERO!!'
      endif

      prefactor = 1.Q0/(s*t-m2sq*m4sq)
      ms = -s/musq
      mt = -t/musq
      mm2sq = - m2sq/musq
      mm3sq = - m3sq/musq
      mm4sq = - m4sq/musq
      sot = s/t
 
      lnms = Log(abs(ms)) - ipi*theta_QUAD(-ms)
      lnmt = Log(abs(mt)) - ipi*theta_QUAD(-mt)
      lnmm2sq = Log(abs(mm2sq)) - ipi*theta_QUAD(-mm2sq)
      lnmm3sq = Log(abs(mm3sq)) - ipi*theta_QUAD(-mm3sq)
      lnmm4sq = Log(abs(mm4sq)) - ipi*theta_QUAD(-mm4sq)
      lnsot = Log(abs(sot)) + ipi * im_part_QUAD(s,t)

     
      arg1 = 1.Q0-m2sq/s
      arg2 = 1.Q0-m4sq/t
      arg3 = 1.Q0-m2sq*m4sq/s/t
      

c     m2sq/s < 0 
      if (arg1.gt.1q0) then
         logarg1=Log(arg1)
         li2arg1 = -dilog_QUAD(1q0/arg1) - logarg1*logarg1/2q0+pi2o3 
     #        - ipi*logarg1*im_part_QUAD(m2sq,s)
      else
         li2arg1 = dilog_QUAD(arg1)
      endif

c     m4sq/t < 0
      if (arg2.gt.1q0)  then
         logarg2=Log(arg2)
         li2arg2 = -dilog_QUAD(1q0/arg2) - logarg2*logarg2/2q0+pi2o3 
     #        - ipi*logarg2*im_part_QUAD(m4sq,t)
      else
         li2arg2 = dilog_QUAD(arg2)
      endif

      im_a = m2sq/s*(1q0/m2sq-1q0/s)
      im_b = m4sq/t*(1q0/m4sq-1q0/t)
      im_ab = m2sq*m4sq/s/t*(1q0/m2sq-1q0/s+1q0/m4sq-1q0/t)

c     m2sq*m4sq/(s*t)  < 0
      if (arg3.gt.1q0) then
         logarg3=Log(arg3)
         li2arg3 = -dilog_QUAD(1q0/arg3) - logarg3*logarg3/2q0+pi2o3 
     #        - ipi*logarg3*sgn_QUAD(im_ab)
      else
         li2arg3 = dilog_QUAD(arg3)
      endif

      if (arg3.lt.0d0) then
         lnarg3 = log(abs(arg3)) - ipi*sgn_QUAD(im_ab)
      else
         lnarg3 = log(arg3)
      endif
c Not checked, this line, but without it, the results does not agree with
C QCD looptools      

      li2arg3 = li2arg3  + 2*ipi*eta3_QUAD(im_a,im_b,im_ab)*lnarg3


      ris = (lnms*lnms+lnmt*lnmt-lnmm2sq*lnmm2sq-lnmm4sq*lnmm4sq)/2q0
     #     +lnmm2sq*lnmm3sq
     #     -(lnmm3sq+lnmm2sq-lnmt)*lnmt+lnmm3sq*lnmm4sq
     #     -(lnmm4sq+lnmm3sq-lnms)*lnms

      ris = (ris + 2q0*(-li2arg1 - li2arg2 + li2arg3) 
     #     - lnsot*lnsot)*prefactor
      
    
      D03m_fin_QUAD_BDK = ris
c      D03m_fin_QUAD_BDK =QCMPLX(real(D03m_fin_QUAD_BDK),imag(D03m_fin_QUAD_BDK))


      return
    

      end


c    ***********   D04m_fin_QUAD(s,t,m1sq,m2sq,m3sq,m4sq,musq)   *************
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
c                 this box is fin_QUADite 
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+p3)^2 =  
c          N_ep * D04m_fin_QUAD(s,t,m1sq,m2sq,m3sq,m4sq,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)
c

c	m2 ... arbitrary, but positive -> set m2 = 1d0


c
c
c     D0 with 4 nonzero external masses
c     and zero internal propagators      
c     use Eq.(41) from Denner, Nierste, Scharf		     
c

c      complex*32 function D04m_fin_QUAD(s,t,m1sq,m2sq,m3sq,m4sq,musq)
      subroutine D04mfin_QUAD(s,t,m1sq,m2sq,m3sq,m4sq,musq,D04m_fin_QUAD)
      implicit none
      complex*32 D04m_fin_QUAD
      real*16  s,t,m1sq,m2sq,m3sq,m4sq,musq

      complex*32 ipi
      parameter (ipi=(0q0,3.14159265358979323846264338327950288q0))
      real*16  pi,pi2,pi2o3
      parameter (pi=3.14159265358979323846264338327950288q0,
     #  pi2o3  = 3.28986813369645287294483033329205038q0 ,
     #     pi2 =9.86960440108935861883449099987615114q0 )

      integer k

      real*16   dilog_QUAD,  eta3_QUAD, theta_QUAD, sgn_QUAD      
      complex*32 ln_mx(2),lnk12,lnk13,lnk14,lnk23,total

      real*16  m2,k12,k13,k14,k23,k24,k34
      real*16  a,b,c,d,sqrtdelta,x(2)
      real*16  k34ok13,im_k34ok13,k24ok12,im_k24ok12
      real*16  im_x(2),im_k34ok13_x(2),im_k24ok12_x(2),arg,logarg
      complex*32 eta_k34ok13(2),eta_k24ok12(2)
      complex*32 lik34ok13(2), logk34ok13(2)
      complex*32 lik24ok12(2), logk24ok12(2)

      real*16 delta,Ia2
      External dilog_QUAD
      m2 = 1q0                  ! > 0 ALWAYS!
      
      k12 = -m1sq/m2
      k13 = -s/m2
      k14 = -m4sq/m2
      k23 = -m2sq/m2
      k24 = -t/m2
      k34 = -m3sq/m2
      

      a = k24*k34
      Ia2=1q0/(2q0*a)
      b = k13*k24+k12*k34-k14*k23
      c = k12*k13
      d = k23
      delta=(b*b-4*a*c)


      if (delta.lt.0q0) then
         write(*,*) 'ERROR: delta less than zero!!'
         write(*,*) 'RETURN ZERO'
c use BDK version of D04m, since in this case it's REAL!!!
c         D04m_fin_QUAD = D04m_fin_QUAD_real(s,t,m1sq,m2sq,m3sq,m4sq,musq)
         D04m_fin_QUAD = 0q0
c         write(*,*) 's,t,m1sq,m2sq,m3sq,m4sq'
c         write(*,*) s,t,m1sq,m2sq,m3sq,m4sq
         return         
c         stop
      else
         sqrtdelta = sqrt(delta)
      endif


      x(1) = (-b + sqrtdelta)*Ia2
      x(2) = (-b - sqrtdelta)*Ia2

c      write(*,*) 'x(1) x(2)',x(1),x(2)

      k34ok13 = k34/k13
      im_k34ok13 = (1q0-m3sq/s)/s
      
      k24ok12 = k24/k12
      im_k24ok12 = (1q0-t/m1sq)/m1sq

      im_x(1) = m2sq/sqrtdelta
      im_x(2) = -im_x(1)


      do k=1,2
         im_k34ok13_x(k) = im_k34ok13*x(k) + k34ok13*im_x(k)
         im_k24ok12_x(k) = im_k24ok12*x(k) + k24ok12*im_x(k)
         eta_k34ok13(k) = 
     #        2*ipi* eta3_QUAD(-im_x(k),im_k34ok13,-im_k34ok13_x(k))
         eta_k24ok12(k) = 2q0*ipi*
     #         eta3_QUAD(-im_x(k),im_k24ok12,-im_k24ok12_x(k))
      ln_mx(k) = Log(abs(x(k))) - theta_QUAD(x(k))*sgn_QUAD(im_x(k))*ipi
      enddo

      lnk12 = Log(abs(k12))-ipi*theta_QUAD(-k12)
      lnk13 = Log(abs(k13))-ipi*theta_QUAD(-k13)
      lnk14 = Log(abs(k14))-ipi*theta_QUAD(-k14)
      lnk23 = Log(abs(k23))-ipi*theta_QUAD(-k23)

c      write(*,*) lnk12,lnk13,lnk14,lnk23

      do k=1,2
         arg = 1q0 + k34ok13*x(k)
         if (arg.gt.1q0) then
c     complex dilog_QUAD
            logarg = Log(arg)
       lik34ok13(k) = -dilog_QUAD(1q0/arg) - 0.5q0*logarg*logarg + pi2o3
     #           + ipi * sgn_QUAD(im_k34ok13_x(k))*logarg
         else
            lik34ok13(k) = dilog_QUAD(arg)
         endif
         logk34ok13(k) = Log(abs(arg)) 
     #        + ipi*theta_QUAD(-arg) * sgn_QUAD(im_k34ok13_x(k))
c         write(*,*) lik34ok13(k),logk34ok13(k)
      enddo


      do k=1,2
         arg = 1q0 + k24ok12*x(k)
         if (arg.gt.1q0) then
c     complex dilog_QUAD
            logarg = Log(arg)
       lik24ok12(k) = -dilog_QUAD(1q0/arg) - 0.5q0*logarg*logarg + pi2o3
     #           + ipi * sgn_QUAD(im_k24ok12_x(k)) * logarg
         else
            lik24ok12(k) = dilog_QUAD(arg)
         endif
         logk24ok12(k) = Log(abs(arg)) 
     #        + ipi*theta_QUAD(-arg) * sgn_QUAD(im_k24ok12_x(k))
      enddo


c      write(*,*) '=========>', 
c     #     (lnk12+lnk13-lnk14-lnk23)*(-ln_mx(1)+ln_mx(2))

      total = 0q0
      do k=1,2
         total = total + (-1)**k*( - 0.5q0 *ln_mx(k)*ln_mx(k) 
     #        - lik34ok13(k) - eta_k34ok13(k) * logk34ok13(k)
     #        - lik24ok12(k) - eta_k24ok12(k) * logk24ok12(k))
      enddo
      
      total = total + (lnk12+lnk13-lnk14-lnk23)*(-ln_mx(1)+ln_mx(2))

      D04m_fin_QUAD = total/a/m2**2/(x(1)-x(2))
      return
      end


cfc************************************************************************
cfcc        FUNCTION D04(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4)
cfc        FUNCTION D04QUAD1(P1t,P2t,P3t,P4t,P12t,P23t,M1t,M2t,M3t,M4t)
cfc************************************************************************
cfc*  SCALAR 4-POINT FUNCTION WITH AT LEAST ONE MASS ZERO                 *
cfc*  P1,P2,P3,P4 = SQUARED EXTERNAL MOMENTA			       *
cfc*  P12 = (p1+p2)**2,  P23 = (p2+p3)**2				       *
cfc*----------------------------------------------------------------------*
cfc*  2.1.92  SD	         					       *
cfc*  rearrangement to exploit massless external momenta   14.3.01  DZ    *
cfc*  Modified: Michael Kubocz                                            *
cfc*  Interception of NANs e.g. caused by Log(0) etc. (see below)         *
cfc************************************************************************
cfcc        IMPLICIT REAL*16 (A-Z)
cfc        implicit none
cfc        real*16 p1t,p2t,p3t,p4t,p12t,p23t
cfc        real*16 m1t,m2t,m3t,m4t
cfc        complex*32 D04QUAD1
cfc
cfc	REAL*16 M(4),P(4,4),K(4,4)
cfc        real*16 pi,eps,eps1
cfc        real*16 im1,im2
cfc        real*16 m1,m2,m3,m4
cfc    
cfc
cfc        real*16 m02,m12,m22,m32,m42
cfc        real*16 mm0,mm1,mm2,mm3,mm4 
cfc        real*16 p1,p2,p3,p4,p12,p23
cfc   
cfc
cfc        real*16 q0,q1,q2,q3,q4,q00,q12,q23
cfc	COMPLEX*32 A1,A2,A3,A4,SWAP
cfc	COMPLEX*32 SS(4), XX(2), X(2,4),RS(4,4)
cfc	COMPLEX*32 S0(4),XX0(2),X0(2,4), R(4,4),G(2)
cfc        
cfc	COMPLEX*32 AA,BB,CC,DD,IEPS,H,HH,L1,L2,L3,L4
cfc	INTEGER I,J, i1,i2,i3,i4
cfc        complex*32 D1,D2,D3,D4
cfc        
cfc        COMPLEX*32 QUADCDILOG,ETA_QUAD,SQE_QUAD,ETAS_QUAD
cfc        EXTERNAL QUADCDILOG,ETA_QUAD,SQE_QUAD,ETAS_QUAD
cfc
cfc        D1=QCMPLX(0q0,0q0)
cfc        D2=QCMPLX(0q0,0q0)
cfc        D3=QCMPLX(0q0,0q0)
cfc        D4=QCMPLX(0q0,0q0)
cfc        eps1=1d-7
cfc
cfc        if(abs(P1t).le.eps1) then
cfc           P1=0q0
cfc        else
cfc           P1=P1t
cfc        endif
cfc        if(abs(P2t).le.eps1) then
cfc           P2=0q0
cfc        else
cfc           P2=P2t
cfc        endif
cfc        if(abs(P3t).le.eps1) then
cfc           P3=0q0
cfc        else
cfc           P3=P3t
cfc        endif
cfc        if(abs(P4t).le.eps1) then
cfc           P4=0q0
cfc        else
cfc           P4=P4t
cfc        endif
cfc        if(abs(P12t).le.eps1) then
cfc           P12=0q0
cfc        else
cfc           P12=P12t
cfc        endif
cfc        if(abs(P23t).le.eps1) then
cfc           P23=0q0
cfc        else
cfc           P23=P23t
cfc        endif
cfc        if(abs(M1t).le.eps1) then
cfc           M1=0q0
cfc        else
cfc           M1=M1t
cfc        endif
cfc        if(abs(M2t).le.eps1) then
cfc           M2=0q0
cfc        else
cfc           M2=M2t
cfc        endif
cfc        if(abs(M3t).le.eps1) then
cfc           M3=0q0
cfc        else
cfc           M3=M3t
cfc        endif
cfc        if(abs(M4t).le.eps1) then
cfc           M4=0q0
cfc        else
cfc           M4=M4t
cfc        endif
cfc
cfc        MM1=M1
cfc        MM2=M2
cfc        MM3=M3
cfc        MM4=M4
cfc        M12=M1*M1
cfc        M22=M2*M2
cfc        M32=M3*M3
cfc        M42=M4*M4
cfc        Q1=P1
cfc        Q2=P2
cfc        Q3=P3
cfc	Q4=P4
cfc        Q12=P12
cfc        Q23=P23
cfc
cfcC	IS AT LEAST ONE MASS ZERO ???
cfc	IF (MM1*MM2*MM3*MM4.NE.0q0) GOTO 130
cfc
cfc
cfcC	****** NO MASS EQUAL TO ZERO ******
cfc130	CONTINUE
cfc	EPS=1q-18
cfc	IEPS=QCMPLX(0q0,EPS)
cfcc check for massless external momentum: excellent candidate for p13,
cfcc leading to r13 >=1 and real.
cfcc$$$        if (q1*q2*q3*q4.eq.0d0) then !org
cfc        if (abs(q1*q2*q3*q4).le.eps) then
cfcc$$$           if (q2.eq.0d0) then
cfc           if (abs(q2).le.eps) then
cfc              I1 = 2
cfc              I2 = 3
cfc              I3 = 1
cfc              I4 = 4
cfcc$$$           elseif (q1.eq.0d0) then
cfc           elseif (abs(q1).le.eps) then
cfc              I1 = 1
cfc              I2 = 3
cfc              I3 = 2
cfc              I4 = 4
cfcc$$$           elseif (q3.eq.0d0) then
cfc           elseif (abs(q3).le.eps) then
cfc              I1 = 2
cfc              I2 = 4
cfc              I3 = 1
cfc              I4 = 3
cfc           else
cfc              I1 = 1
cfc              I2 = 4
cfc              I3 = 2
cfc              I4 = 3
cfc           endif
cfc           M(i1)=MM1
cfc           M(i2)=MM2
cfc           M(i3)=MM3
cfc           M(i4)=MM4
cfc           P(i1,i2)=Q1
cfc           P(i3,i2)=Q2
cfc           P(i3,i4)=Q3
cfc           P(i1,i4)=Q4
cfc           P(i1,i3)=Q12
cfc           P(i3,i1)=Q12
cfc           P(i2,i4)=Q23
cfc           P(i4,i2)=Q23
cfc	ELSEIF( ABS((MM1**2+MM3**2-Q12)/MM1/MM3).LT.2q0 ) THEN
cfcC	R13 WOULD BE NOT REAL. -> PERMUTATION! -> R(2,4) IS NOT REAL.
cfc	   M(1)=MM2
cfc	   M(2)=MM3
cfc	   M(3)=MM4
cfc	   M(4)=MM1
cfc	   P(1,2)=Q2
cfc	   P(1,3)=Q23
cfc	   P(1,4)=Q1
cfc	   P(2,3)=Q3
cfc	   P(2,4)=Q12
cfc	   P(3,4)=Q4
cfc	ELSE
cfcC	R(1,3) IS REAL.
cfc	   M(1)=MM1
cfc	   M(2)=MM2
cfc	   M(3)=MM3
cfc	   M(4)=MM4
cfc	   P(1,2)=Q1
cfc	   P(1,3)=Q12
cfc	   P(1,4)=Q4
cfc	   P(2,3)=Q2
cfc	   P(2,4)=Q23
cfc	   P(3,4)=Q3
cfc	ENDIF
cfc
cfc	DO 11 J=2,4
cfc	DO 11 I=1,J-1
cfc	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
cfc	R(I,J) =SQE_QUAD(QCMPLX(1q0,0q0),QCMPLX(-K(I,J),0q0),
cfc     *	            QCMPLX(1q0,0q0))
cfc        IF( K(I,J).LT.-2q0 ) THEN
cfcc        IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
cfc	   RS(I,J)=SQE_QUAD(QCMPLX(1q0,0q0),QCMPLX(-K(I,J),EPS),
cfc     *	               QCMPLX(1q0,0q0))
cfc	ELSE
cfc	   RS(I,J)=R(I,J)
cfc	ENDIF
cfc11	CONTINUE
cfc 
cfc        
cfc
cfc	SS(1)=RS(1,2)
cfc	SS(2)=RS(2,3)
cfc	SS(3)=RS(3,4)
cfc	SS(4)=RS(1,4)
cfc	S0(1)=R(1,2)
cfc	S0(2)=R(2,3)
cfc	S0(3)=R(3,4)
cfc	S0(4)=R(1,4)
cfc	AA=K(3,4)/R(2,4)+R(1,3)*K(1,2)-K(1,4)*R(1,3)/R(2,4)-K(2,3)
cfc	BB=(R(2,4)-1q0/R(2,4))*(R(1,3)-1q0/R(1,3))
cfc     *		+K(1,2)*K(3,4)-K(1,4)*K(2,3)
cfc	CC=K(1,2)/R(1,3)+R(2,4)*K(3,4)-K(1,4)*R(2,4)/R(1,3)-K(2,3)
cfc	DD=K(2,3)-R(1,3)*K(1,2)-R(2,4)*K(3,4)+R(1,3)*R(2,4)*K(1,4)
cfc	XX(1)=SQE_QUAD(AA,BB,CC+IEPS*DD)
cfc	XX(2)=(CC+IEPS*DD)/AA/XX(1)
cfc	XX0(1)=SQE_QUAD(AA,BB,CC)
cfc	XX0(2)=CC/AA/XX0(1)
cfcc	IF (ABS(DREAL(XX0(1)-XX(2))).LT.ABS(DREAL(XX0(1)-XX(1)))) THEN
cfc	IF (ABS(XX0(1)-XX(2)).LT.ABS(XX0(1)-XX(1))) THEN
cfc	  SWAP  =XX0(1)
cfc	  XX0(1)=XX0(2)
cfc	  XX0(2)=SWAP
cfc	ENDIF
cfc
cfc	DO 12 I=1,2
cfc	G(I)  =SIGN( 1Q0,real(AA*(XX(I)-XX(3-I))) )
cfc	 X(I,1)= XX(I)/R(2,4)
cfc	X0(I,1)=XX0(I)/R(2,4)
cfc	 X(I,2)= XX(I)/R(2,4)*R(1,3)
cfc	X0(I,2)=XX0(I)/R(2,4)*R(1,3)
cfc	 X(I,3)= XX(I)*R(1,3)
cfc	X0(I,3)=XX0(I)*R(1,3)
cfc	 X(I,4)= XX(I)
cfc	X0(I,4)=XX0(I)
cfc12	CONTINUE
cfc
cfc	D04QUAD1 = QCMPLX(0q0,0q0)
cfc	DO 13 I=1,2
cfc	DO 13 J=1,4
cfc	A1 = 1q0+X0(I,J)*S0(J) + ABS(1q0+X0(I,J)*S0(J))*IEPS*
cfc     *				  SIGN(1q0,imag(X(I,J)*SS(J)))
cfc	A2 = 1q0+X0(I,J)/S0(J) + ABS(1q0+X0(I,J)/S0(J))*IEPS*
cfc     *				  SIGN(1q0,imag(X(I,J)/SS(J)))
cfc
cfcc  org. code:
cfcc$$$           D04 = D04 + (-1D0)**(I+J)*(
cfcc$$$     &          LI2(A1)+ETA(-X(I,J),SS(J))*LOG(A1)
cfcc$$$     &          +LI2(A2)+ETA(-X(I,J),1D0/SS(J))*LOG(A2))
cfc
cfc        if(abs(ETA_QUAD(-X(I,J),SS(J))).ne.0q0) then
cfc           D1=ETA_QUAD(-X(I,J),SS(J))*LOG(A1)
cfc        else
cfc           D1=QCMPLX(0q0,0q0)
cfc        endif
cfc        if(abs(ETA_QUAD(-X(I,J),1q0/SS(J))).ne.0q0) then
cfc           D2=ETA_QUAD(-X(I,J),1q0/SS(J))*LOG(A2)
cfc        else
cfc           D2=QCMPLX(0q0,0q0)
cfc        endif
cfc        D04QUAD1=D04QUAD1+(-1q0)**(I+J)*(QUADCDILOG(A1)+QUADCDILOG(A2)+D1+D2)
cfcc   The enquiry avoids occurrence of NANs causing by LOG(A1) for A1=0 
cfcc   and LOG(A2) for A2=0. At that points also ETA_QUAD(-X(I,J),1D0/SS(J)) 
cfcc   or ETA_QUAD(-X(I,J),SS(J) are 0. (Michael Kubocz)
cfc
cfc13	CONTINUE
cfc        
cfc        
cfc        
cfc
cfc
cfcc        print*,'DIMAG(R(1,3))',DIMAG(R(1,3))
cfcc	IF( DIMAG(R(1,3)).EQ.0D0 ) THEN !org (makes troubles in squark pentagons)
cfc	IF( abs(imag(R(1,3))).le.eps ) THEN
cfc	DO 14 I=1,2
cfc	   A1 = (K(1,3)-2q0*R(1,3))/XX0(I)
cfc     *		      -R(1,3)*K(1,4)+K(3,4)
cfc     	   A2 = ((K(2,4)-2q0*R(2,4))*R(1,3)*XX0(I)
cfc     *		      -R(2,4)*K(3,4)+K(2,3))/DD
cfc	   A3 = (K(1,3)-2q0*R(1,3))*R(2,4)/XX0(I)
cfc     *		      -R(1,3)*K(1,2)+K(2,3)
cfc	   A4 = ((K(2,4)-2q0*R(2,4))*XX0(I)
cfc     *		      -R(2,4)*K(1,4)+K(1,2))/DD
cfc	   L1 = LOG( A1-ABS(A1)*IEPS )
cfc     	   L2 = LOG( A2+ABS(A2)*IEPS*G(I)*SIGN(1q0,real(R(1,3))
cfc     *				        	  *imag(RS(2,4))) )
cfc	   L3 = LOG( A3-ABS(A3)*IEPS )
cfc	   L4 = LOG( A4+ABS(A4)*IEPS*G(I)*SIGN(1q0,imag(RS(2,4))) )
cfc
cfcc org. code:
cfcc$$$	   D04 = D04 
cfcc$$$     &         + (3D0-2D0*I)*(
cfcc$$$     *		       ETAS_QUAD( -XX(I),R(1,3),RS(1,3) )
cfcc$$$     *		          *( LOG(R(1,3)*XX(I)) + L1 + L2 )
cfcc$$$     *		     + ETAS_QUAD( -XX(I),1D0/R(2,4),1D0/RS(2,4) )
cfcc$$$     *		          *( LOG(XX(I)/R(2,4)) + L3 + L4 )
cfcc$$$     *		     - ( ETAS_QUAD( -XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4) )
cfcc$$$     *		       + ETA( RS(1,3),1D0/RS(2,4) )                  )
cfcc$$$     *		        *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 + L2 )
cfcc$$$     *	  	     + ETA( RS(1,3),1D0/RS(2,4) )
cfcc$$$     *		       *ETAS_QUAD(-XX(I),-R(1,3)/R(2,4),-RS(1,3)/RS(2,4))   )
cfc
cfc            if(abs(ETAS_QUAD(-XX(I),R(1,3),RS(1,3))).ne.0q0) then
cfc               D1=ETAS_QUAD(-XX(I),R(1,3),RS(1,3))*(LOG(R(1,3)*XX(I))+L1+L2)
cfc            else
cfc               D1=QCMPLX(0q0,0q0)
cfc            endif
cfc            if(abs(ETAS_QUAD(-XX(I),1q0/R(2,4),1q0/RS(2,4))).ne.0q0) then
cfc               D2=ETAS_QUAD(-XX(I),1q0/R(2,4),1q0/RS(2,4))*(LOG(XX(I)/R(2,4))
cfc     &           +L3+L4)
cfc            else
cfc               D2=QCMPLX(0q0,0q0)
cfc            endif
cfc            if((abs(ETAS_QUAD(-XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4))).ne.0q0)
cfc     &          .or.(abs(ETA_QUAD(RS(1,3),1q0/RS(2,4))).ne.0q0)) then 
cfc               D3=-(ETAS_QUAD(-XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4))
cfc     &            +ETA_QUAD(RS(1,3),1q0/RS(2,4)))*(LOG(XX(I)*R(1,3)/R(2,4))+L3+L2)
cfc            else
cfc               D3=QCMPLX(0q0,0q0)
cfc            endif
cfc            D4=ETA_QUAD(RS(1,3),1q0/RS(2,4))*ETAS_QUAD(-XX(I),-R(1,3)/R(2,4),-RS(1,3)/RS(2,4))
cfc            D04QUAD1=D04QUAD1+(3q0-2q0*I)*(D1+D2+D3+D4)
cfcc   The enquiry avoids occurrence of NANs causing by LOG(0). At that points 
cfcc   also ETA(...) and ETAS_QUAD(...) are 0. (Michael Kubocz)
cfc
cfc
cfc         
cfc
cfc         
cfc
cfc
cfc14	CONTINUE
cfc	ELSE
cfc	DO 15 I=1,2
cfc	   L1 = LOG( R(2,4)/XX0(I)+XX0(I)/R(2,4)+K(1,2)
cfc     *		     -XX0(I)/R(2,4)*EPS*BB*G(I) )
cfc	   L2 = LOG( R(1,3)*XX0(I)+1Q0/XX0(I)/R(1,3)+K(3,4)
cfc     *		     -XX0(I)*R(1,3)*EPS*BB*G(I) )
cfc	   L3 = LOG( R(1,3)/R(2,4)*XX0(I)+R(2,4)/XX0(I)/R(1,3)+K(2,3)
cfc     *		     -XX0(I)*R(1,3)/R(2,4)*EPS*BB*G(I) )
cfc
cfc	   D04QUAD1 = D04QUAD1 
cfc     &          + (3q0-2q0*I)*(
cfc     *		     ETA_QUAD(-XX(I),1q0/R(2,4))
cfc     *		      *( LOG(XX(I)/R(2,4)) + L1 )
cfc     *		    +ETA_QUAD(-XX(I),R(1,3))
cfc     *		      *( LOG(R(1,3)*XX(I)) + L2 )
cfc     *		    -( ETA_QUAD(-XX(I),R(1,3)/R(2,4))
cfc     *		      +ETA_QUAD(R(1,3),1Q0/R(2,4)) )
cfc     *		      *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 )
cfc     *	  	    +ETA_QUAD(R(1,3),1Q0/R(2,4))
cfc     *		      *ETA_QUAD(-XX(I),-R(1,3)/R(2,4))
cfc     *		      *(1q0-G(I)*SIGN(1q0,real(BB)))	    )
cfc
cfc                  
cfc15	CONTINUE
cfc	ENDIF
cfc
cfc	D04QUAD1 = D04QUAD1/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
cfc	RETURN
cfc 
cfcC--->	***************** SPEZIELL ( --> T.SACK-PROMOTION )
cfcC	D1=Q12-M12
cfcC	D2=Q2 -M22
cfcC	D3=Q3 -M42
cfcC	IF ((D1*D2.LE.0D0).OR.(D2*D3.LE.0D0)) THEN
cfcC	   WRITE(*,*) 'THE CASE OF DIFFERENT SIGNS OF THE D1,D2,D3'
cfcC	   WRITE(*,*) 'IN D04(...) IS NOT IMPLEMENTED !!!'
cfcC	   STOP
cfcC	ENDIF
cfcC	NM1=ABS(MM1/D1)
cfcC	NM2=ABS(MM2/D2)
cfcC	NM3=ABS(MM4/D3)
cfcC	NP1=Q2/D2**2+Q12/D1**2+(Q1-Q2-Q12)/D1/D2
cfcC	NP2=Q2/D2**2+ Q3/D3**2+(Q23-Q2-Q3)/D2/D3
cfcC	NP3=Q3/D3**2+Q12/D1**2+(Q4-Q3-Q12)/D1/D3
cfcC	D04=C04(NP1,NP2,NP3,NM1,NM2,NM3)/D1/D2/D3
cfc
cfcC	*************** ALLGEMEIN
cfc
cfc
cfcC	****** ONLY MM3 IS ZERO ******
cfc30	CONTINUE
cfc	EPS=1q-47
cfc	IEPS=QCMPLX(0Q0,EPS)
cfc	M(1)=MM1
cfc	M(2)=MM2
cfc	M(3)=10Q0
cfc	M(4)=MM4
cfc	P(1,2)=Q1
cfc	P(1,3)=Q12
cfc	P(1,4)=Q4
cfc	P(2,3)=Q2
cfc	P(2,4)=Q23
cfc	P(3,4)=Q3
cfc	DO 1 J=2,4
cfc	DO 1 I=1,J-1
cfc	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
cfc	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
cfc	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
cfc	R(I,J) =SQE_QUAD(QCMPLX(1Q0,0Q0),QCMPLX(-K(I,J),0Q0),
cfc     *	            QCMPLX(1Q0,0Q0))
cfc	IF( imag(R(I,J)).EQ.0Q0 ) THEN
cfc	   RS(I,J)=SQE_QUAD(QCMPLX(1Q0,0Q0),QCMPLX(-K(I,J),EPS),
cfc     *	               QCMPLX(1Q0,0Q0))
cfc	ELSE
cfc	   RS(I,J)=R(I,J)
cfc	ENDIF
cfc1	CONTINUE
cfc	SS(1)=RS(1,2)
cfc	SS(2)=RS(2,3)
cfc	SS(3)=RS(3,4)
cfc	SS(4)=RS(1,4)
cfc	AA=K(3,4)/R(2,4)-K(2,3)
cfc	BB=K(1,3)*(1Q0/R(2,4)-R(2,4))+K(1,2)*K(3,4)-K(1,4)*K(2,3)
cfc	CC=K(1,2)*K(1,3)-K(1,3)*K(1,4)*R(2,4)+R(2,4)*K(3,4)-K(2,3)
cfc	DD=K(2,3)-R(2,4)*K(3,4)
cfc	XX(1)=SQE_QUAD(AA,BB,CC+IEPS*DD)
cfc	XX(2)=(CC+IEPS*DD)/AA/XX(1)
cfc	DO 2 I=1,2
cfc	X(I,1)=XX(I)/R(2,4)
cfc	X(I,2)=XX(I)/R(2,4)*R(1,3)
cfc	X(I,3)=XX(I)*R(1,3)
cfc	X(I,4)=XX(I)
cfc2	CONTINUE
cfc	D04QUAD1 = QCMPLX(0Q0,0Q0)
cfc	DO 3 I=1,2
cfc	D04QUAD1 = D04QUAD1 + (2Q0*I-3Q0)*(
cfc     *		quadcdilog(1Q0+SS(4)*X(I,4))
cfc     *	       -quadcdilog(1Q0+SS(1)*X(I,1))
cfc     *	       +quadcdilog(1Q0+X(I,4)/SS(4))
cfc     *	       -quadcdilog(1Q0+X(I,1)/SS(1))
cfc     *	       +ETA_QUAD(-X(I,4),SS(4))*LOG(1Q0+SS(4)*X(I,4))
cfc     *	       -ETA_QUAD(-X(I,1),SS(1))*LOG(1Q0+SS(1)*X(I,1))
cfc     *	       +ETA_QUAD(-X(I,4),1Q0/SS(4))*LOG(1Q0+X(I,4)/SS(4))
cfc     *	       -ETA_QUAD(-X(I,1),1Q0/SS(1))*LOG(1Q0+X(I,1)/SS(1))
cfc     *	       -quadcdilog(1Q0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
cfc     *	       +quadcdilog(1Q0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))
cfc     *	       -ETA_QUAD(-X(I,4),(K(3,4)-IEPS)/(K(1,3)-IEPS))
cfc     *	           *LOG(1Q0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
cfc     *	       +ETA_QUAD(-X(I,1),(K(2,3)-IEPS)/(K(1,3)-IEPS))
cfc     *	           *LOG(1Q0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))   )
cfc	IF (imag(R(2,4)).NE.0D0) THEN
cfc	   H=ETA_QUAD(-1Q0/XX(I),R(2,4))
cfc	ELSE
cfc	   H=QCMPLX(0Q0,0Q0)
cfc	   IF (real(R(2,4)).LT.0D0) THEN
cfc	      HH=-1Q0/XX(I)
cfc	      IM1=imag(HH)
cfc	      IM2=imag(RS(2,4))
cfc              pi = 4.Q0*atan(1.Q0)
cfc	      IF ((IM1.GT.0Q0).AND.(IM2.GT.0Q0)) THEN
cfc	         H=-QCMPLX(0Q0,2Q0*PI)
cfc	      ENDIF
cfc	      IF ((IM1.LT.0Q0).AND.(IM2.LT.0Q0)) THEN
cfc	         H=+QCMPLX(0Q0,2Q0*PI)
cfc	      ENDIF
cfc	   ENDIF
cfc	ENDIF
cfc	D04QUAD1 = D04QUAD1 + (2Q0*I-3Q0)*
cfc     *	          H*( LOG( (K(1,2)-R(2,4)*K(1,4)
cfc     *			  +XX(I)*(1Q0/R(2,4)-R(2,4)))/DD )
cfc     *		     +LOG(K(1,3)-IEPS) )
cfc3	CONTINUE
cfc	D04QUAD1 = D04QUAD1/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
cfc	RETURN
cfc
cfcC	****** ONLY MM2 AND MM3 ARE ZERO ******
cfc40	CONTINUE
cfc	EPS=1D-17
cfc	IEPS=QCMPLX(0Q0,EPS)
cfc
cfc	M(1)=MM1
cfc	M(2)=10Q0
cfc	M(3)=10Q0
cfc	M(4)=MM4
cfc	P(1,2)=Q1
cfc	P(1,3)=Q12
cfc	P(1,4)=Q4
cfc	P(2,3)=Q2
cfc	P(2,4)=Q23
cfc	P(3,4)=Q3
cfc	DO 4 J=2,4
cfc	DO 4 I=1,J-1
cfc	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
cfc	IF (I.EQ.2) K(I,J)=K(I,J)-M(I)/M(J)
cfc	IF (J.EQ.2) K(I,J)=K(I,J)-M(J)/M(I)
cfc	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
cfc	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
cfc	R(I,J) =SQE_QUAD(QCMPLX(1Q0,0Q0),QCMPLX(-K(I,J),0Q0),
cfc     *	            QCMPLX(1Q0,0Q0))
cfc	IF( imag(R(I,J)).EQ.0Q0 ) THEN
cfc	   RS(I,J)=SQE_QUAD(QCMPLX(1Q0,0Q0),QCMPLX(-K(I,J),EPS),
cfc     *	               QCMPLX(1Q0,0Q0))
cfc	ELSE
cfc	   RS(I,J)=R(I,J)
cfc	ENDIF
cfc4	CONTINUE
cfc	SS(1)=RS(1,2)
cfc	SS(2)=RS(2,3)
cfc	SS(3)=RS(3,4)
cfc	SS(4)=RS(1,4)
cfc	AA=K(2,4)*K(3,4)-K(2,3)
cfc	BB=K(1,3)*K(2,4)+K(1,2)*K(3,4)-K(1,4)*K(2,3)
cfc	CC=K(1,2)*K(1,3)-K(2,3)
cfc	DD=K(2,3)
cfc	XX(1)=SQE_QUAD(AA,BB,CC+IEPS*DD)
cfc	XX(2)=(CC+IEPS*DD)/AA/XX(1)
cfc	DO 5 I=1,2
cfc	X(I,1)=XX(I)/R(2,4)
cfc	X(I,2)=XX(I)/R(2,4)*R(1,3)
cfc	X(I,3)=XX(I)*R(1,3)
cfc	X(I,4)=XX(I)
cfc5	CONTINUE
cfc	D04QUAD1 = QCMPLX(0Q0,0Q0)
cfc	DO 6 I=1,2
cfc	D04QUAD1 = D04QUAD1 + (2Q0*I-3Q0)*(
cfc     *		quadcdilog(1Q0+SS(4)*X(I,4))
cfc     *	       +quadcdilog(1Q0+X(I,4)/SS(4))
cfc     *	       +ETA_QUAD(-X(I,4),SS(4))*LOG(1Q0+SS(4)*X(I,4))
cfc     *	       +ETA_QUAD(-X(I,4),1Q0/SS(4))*LOG(1Q0+X(I,4)/SS(4))
cfc     *	       -quadcdilog(1Q0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
cfc     *	       -quadcdilog(1Q0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
cfc     *	       -ETA_QUAD(-XX(I),(K(3,4)-IEPS)/(K(1,3)-IEPS))
cfc     *	           *LOG(1Q0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
cfc     *	       -ETA_QUAD(-XX(I),(K(2,4)-IEPS)/(K(1,2)-IEPS))
cfc     *	           *LOG(1Q0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
cfc     *	       +LOG(-XX(I))*( LOG(K(1,2)-IEPS)
cfc     *			     +LOG(K(1,3)-IEPS)-LOG(K(2,3)-IEPS) ) )
cfc6	CONTINUE
cfc	D04QUAD1 = D04QUAD1/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
cfc
cfc	RETURN
cfc
cfc
cfc
cfc	END


***********************************************************************
        FUNCTION ETA_QUAD(C1,C2)                                            
***********************************************************************
*       COMPLEX ETA-FUNKTION                                           
*---------------------------------------------------------------------*
*       8.06.90    ANSGAR DENNER                                       
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*32 ETA_QUAD,C1,C2                                           
        REAL*16     PI,IM1,IM2,IM12                                     
                                                                       
        PI     = 4q0*ATAN(1q0)                                        
        IM1    = IMAG(C1)                                             
        IM2    = IMAG(C2)                                             
        IM12   = IMAG(C1*C2)                                          
                                                                       
        IF(IM1.LT.0q0.AND.IM2.LT.0q0.AND.IM12.GT.0q0) THEN             
            ETA_QUAD = CMPLX(0q0,2q0*PI,16)                                   
        ELSE IF (IM1.GT.0q0.AND.IM2.GT.0q0.AND.IM12.LT.0q0) THEN       
            ETA_QUAD = CMPLX(0q0,-2q0*PI,16)                                  
        ELSE                                                           
            ETA_QUAD = CMPLX(0q0,0q0,16)                                          
        END IF                                                         
        END                                                            

***********************************************************************
        FUNCTION ETAS_QUAD(Y,R,RS)                                            
***********************************************************************
*       MODIFIED ETA-FUNKTION                                           
*---------------------------------------------------------------------*
*       18.1.94   SD                                       
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*32 ETA_QUAD,ETAS_QUAD,Y,R,RS
        REAL*16     PI,IMY,IMRS
                                                                       
        PI     = 4q0*ATAN(1q0)                                        

	IF( IMAG(R).NE.0q0 ) THEN
	    ETAS_QUAD = ETA_QUAD(Y,R)
	ELSE	    
	    IF( REAL(R).GT.0q0 ) THEN
		ETAS_QUAD = CMPLX(0q0,0q0,16)
	    ELSE
	 	IMY  = IMAG(Y)
		IMRS = IMAG(RS)
		ETAS_QUAD = 2q0*CMPLX(0q0,PI,16)*(
     *			(1q0+SIGN(1q0,-IMY))*(1q0+SIGN(1q0,-IMRS))-
     *			(1q0+SIGN(1q0, IMY))*(1q0+SIGN(1q0, IMRS))
     *					  )/4q0
	    ENDIF
	ENDIF
        END                                                            

***********************************************************************
        FUNCTION SQE_QUAD(A,B,C)                                            
***********************************************************************
*       SOLUTION OF QUADRATIC EQUATION				      *
*---------------------------------------------------------------------*
*       13.1.92  SD						      *
***********************************************************************
        IMPLICIT REAL*16 (A-Z)                                        
        COMPLEX*32 A,B,C,SQE_QUAD,X1,X2

	X1=(-B+SQRT(B**2-4q0*A*C))/2q0/A
	X2=(-B-SQRT(B**2-4q0*A*C))/2q0/A

	IF (ABS(X1).GT.ABS(X2)) THEN
	   SQE_QUAD=X1
	ELSE
	   SQE_QUAD=X2
	ENDIF

        END


      real*16 function eta3_QUAD(im_a,im_b,im_ab)
      implicit none
c      real*16  eta3_QUAD
      real*16 im_a,im_b,im_ab
      real*16  theta_QUAD
       eta3_QUAD = theta_QUAD(-im_a)*theta_QUAD(-im_b)*theta_QUAD(im_ab)
     #     - theta_QUAD(im_a)*theta_QUAD(im_b)*theta_QUAD(-im_ab)
      end

      real*16 function theta_QUAD(x)
      implicit none
c      real*16  theta_QUAD
      real*16 x
      if (x.gt.0q0) then
         theta_QUAD = 1q0
      else
         theta_QUAD = 0q0
      endif
      end

      real*16 function sgn_QUAD(x)
      implicit none
c      real*16  sgn_QUAD
      real*16 x
      if (x.gt.0q0) then
         sgn_QUAD = +1q0
      else
         sgn_QUAD = -1q0
      endif
      end


      real*16 function im_part_QUAD(a,b)
      implicit none
      real*16  a,b!,im_part_QUAD
      if ((a/b).gt.0q0) then
         im_part_QUAD = 0q0
      else
         if (a.lt.0q0) then
            im_part_QUAD = +1q0
         else
            im_part_QUAD = -1q0
         endif
      endif
      
      end



