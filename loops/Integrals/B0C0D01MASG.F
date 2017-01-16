
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
      function A0finG(m,musq)
      IMPLICIT NONE
      complex*16 A0finG
      complex*16 A0finG1,A0finG2
      External A0finG1,A0finG2
      real*8 eps
      parameter (eps=1d-7)
      real*8 m,musq

      logical UseAVH, UseFF
      integer ier
      parameter (UseAVH=.false.)
      parameter (UseFF=.true.)
      
      if(UseAVH) then
#ifdef WITH_ONELOOP
         call avh_olo_mu_set(sqrt(musq))
         call avh_olo_onshell(eps)
         call avh_olo_a0m(rslt,m**2)
        
         A0finG = rslt(0)
#else
         call avh_error()
#endif
      else if (UseFF) then
         ier=0
         call vbfffxa0(A0finG,0d0,musq,m**2,ier)
      else
         if(m.ge.1d-7) then
            A0finG=A0finG2(m,musq)
         else
            A0finG=A0finG1(musq)
         endif
      endif

      return
      end


      function A0finGDiv(m,musq,id)
      IMPLICIT NONE
      complex*16 A0finGDiv
      complex*16 A0finG1Div,A0finG2Div
      External A0finG1Div,A0finG2Div
      real*8 eps
      parameter (eps=1d-7)
      real*8 m,musq

      integer id
      logical UseAVH, UseFF
      integer ier
      parameter (UseAVH=.false.)
      parameter (UseFF=.true.)

      if(UseAVH) then
#ifdef WITH_ONELOOP
         call avh_olo_mu_set(sqrt(musq))
         call avh_olo_onshell(eps)
         call avh_olo_a0m(rslt,m**2)
         
         A0finGDiv = rslt(id)
#else
         call avh_error()
#endif
      else if (UseFF) then
         ier=0
         if (id.eq.0) then
            call vbfffxa0(A0finGDiv,0d0,musq,m**2,ier)
         else
            call vbfff_diverror()
         endif
      else
         if(m.ge.1d-7) then
            A0finGDiv=A0finG2Div(m,musq,id)
         else
            A0finGDiv=A0finG1Div(musq,id)
         endif
      endif
         return
      end


      function B0finG(m0,m1,p1sq,musq)
      implicit none
      complex*16 B0finG,B0fin
      complex*16 B0finG1,B0finG2,B0finG3
      complex*16 B0finG4,B0finG5,B0finG6,B2i1e,B1i1e,B1i0e,B0t1
      External B0finG1,B0finG2,B0finG3,B0fin
      External B0finG4,B0finG5,B0finG6,B2i1e,B1i1e,B1i0e,B0t1
      real*8 m0,m1,p1sq,musq
      real*8 eps

      parameter (eps=1d-7)
      logical NotAVH
      parameter (NotAVH=.true.)
c      print*,"NotAVH",NotAVH
c --------------------------------------
c   B0finG2=B0finG(M1,0,s,musq)  s<=>0                               
c --------------------------------------
c   M1,0,s, 
c   M1,0,0
c   0,M,s                ?     = ?B1i1e(s,mis,mus,id)
c   0,M,0
c   B0finG2(m1,0d0,musq)     =   B1i0e(m02,musq,0)

c --------------------------------------
c   B0finG5=B0finG(M1,M2,s,musq)                                 
c --------------------------------------
c  M1,M2,s
c  M,M,s, (but there is another one for this case)
c
c   B0finG5(m0,m1,p1sq,musq) =   B2i1e(p1sq,m02,m12,musq,0)

c --------------------------------------
c   B0finG6=B0finG(M1,M1,0,musq)          
c --------------------------------------
c M,M, 0

c   B0finG6(m1,musq)         =   B0t1(m1,0d0,musq)     
   
c   M1,M2,S    v
c   M1,M2,0    not yet
c   M,M,S      v
c   M,M,0      v
c   M,0,S      v
c   M,0,0      v
c   0,M,S      v
c   0,M,0      v
C   0,0,0      v
c   0,0,s      v
  
      if(NotAVH) then
      if (m0.ge.eps) then
         if (m1.ge.eps) then
            if (abs(p1sq).ge.eps) then
               if(abs((m0-m1)).gt.eps) then
c B0finG(m0,m1,p1sq,musq)
               B0finG= B0finG5(m0,m1,p1sq,musq)
               return
               else
c B0finG(m0,m0,p1sq,musq)                  
               B0finG= B0finG4(m0,p1sq,musq)
               return
               endif
            else   ! p1s =0
c B0finG(m0,m1,0,musq)'
               if(abs((m0-m1)).gt.eps) then
c B0finG(m0,m1,p1sq,musq)
                Print*, 'No program yet'
                stop
               return
               else
c B0finG(m0,m0,0,musq)                  
               B0finG= B0finG6(m0,musq)
               return
               endif
               return
            endif
        else   ! m1s=0
c            if (abs(p1sq).ge.eps) then
c B0finG(m0,0,p1sq,musq)
c              B0finG=B0finG2(m0,p1sq,musq)
c               return
c             else
c B0finG(m0,0,0,musq)
             B0finG= B0finG2(m0,p1sq,musq)! B0finG2(m1,0d0,musq)= B1i0e(m02,musq,0)
               return
c            endif
        endif  ! m1s=0
c m0=0       
      else  
        if (m1.ge.eps) then
c            if (abs(p1sq).ge.eps) then
c B0finG(0,m1,p1sq,musq)                  
c               B0finG= B0finG2(m1,p1sq,musq)
c               return
c               else   ! p1s =0
c B0finG(0,m1,0,musq)'
c           print*, "HERE, HERE"
               B0finG= B0finG2(m1,p1sq,musq) ! B0finG2(m1,0d0,musq)= B1i0e(m02,musq,0)
               return
c            endif
        else   ! m1s=0
c            if (abs(p1sq).ge.eps) then
c B0finG(0,0,p1sq,musq)
              B0finG=B0fin(p1sq,musq)
             return
c             else
c B0finG(0,0,0,musq)
c             B0finG= 0d0
c             return
c            endif
        endif  ! m1s=0
       endif   ! m0s =0 
      return
      else
#ifdef WITH_ONELOOP
         call avh_olo_mu_set(sqrt(musq))
         call avh_olo_onshell(eps)
         call avh_olo_b0m(rslt,p1sq,m0**2,m1**2)
         B0finG = rslt(0)
         return
#else
         call avh_error()
#endif
      endif
    
      end


      function B0finGDiv(m0,m1,p1sq,musq,id)
      implicit none
      complex*16 B0finGDiv,B0finDiv
      complex*16 B0finG1Div,B0finG2Div,B0finG3Div
      complex*16 B0finG4Div,B0finG5Div,B0finG6Div,B2i1e,B1i1e,B1i0e,B0t1Div
      External B0finG1Div,B0finG2Div,B0finG3Div,B0finDiv
      External B0finG4Div,B0finG5Div,B0finG6Div,B2i1e,B1i1e,B1i0e,B0t1Div
      real*8 m0,m1,p1sq,musq
      real*8 eps

      parameter (eps=1d-7)
      integer id
      logical NotAVH
      parameter (NotAVH=.true.)

c --------------------------------------
c   B0finG2=B0finG(M1,0,s,musq)  s<=>0                               
c --------------------------------------
c   M1,0,s, 
c   M1,0,0
c   0,M,s                ?     = ?B1i1e(s,mis,mus,id)
c   0,M,0
c   B0finG2(m1,0d0,musq)     =   B1i0e(m02,musq,0)

c --------------------------------------
c   B0finG5=B0finG(M1,M2,s,musq)                                 
c --------------------------------------
c  M1,M2,s
c  M,M,s, (but there is another one for this case)
c
c   B0finG5(m0,m1,p1sq,musq) =   B2i1e(p1sq,m02,m12,musq,0)

c --------------------------------------
c   B0finG6=B0finG(M1,M1,0,musq)          
c --------------------------------------
c M,M, 0

c   B0finG6(m1,musq)         =   B0t1(m1,0d0,musq)     
   
c   M1,M2,S    v
c   M1,M2,0    not yet
c   M,M,S      v
c   M,M,0      v
c   M,0,S      v
c   M,0,0      v
c   0,M,S      v
c   0,M,0      v
C   0,0,0      v
c   0,0,s      v
      if(NotAVH) then
      if (m0.ge.eps) then
         if (m1.ge.eps) then
            if (abs(p1sq).ge.eps) then
               if(abs((m0-m1)).gt.eps) then
c B0finG(m0,m1,p1sq,musq)
               B0finGDiv= B0finG5Div(m0,m1,p1sq,musq,id)
               return
               else
c B0finG(m0,m0,p1sq,musq)                  
               B0finGDiv= B0finG4Div(m0,p1sq,musq,id)
               return
               endif
            else   ! p1s =0
c B0finG(m0,m1,0,musq)'
               if(abs((m0-m1)).gt.eps) then
c B0finG(m0,m1,p1sq,musq)
                Print*, 'No program yet'
                stop
               return
               else
c B0finG(m0,m0,0,musq)                  
               B0finGDiv= B0finG6Div(m0,musq,id)
               return
               endif  
               return
            endif
        else   ! m1s=0
c            if (abs(p1sq).ge.eps) then
c B0finG(m0,0,p1sq,musq)
c              B0finG=B0finG2(m0,p1sq,musq)
c               return
c             else
c B0finG(m0,0,0,musq)
             B0finGDiv= B0finG2Div(m0,p1sq,musq,id)! B0finG2(m1,0d0,musq)= B1i0e(m02,musq,0)
               return
c            endif
        endif  ! m1s=0
c m0=0       
      else  
        if (m1.ge.eps) then
c            if (abs(p1sq).ge.eps) then
c B0finG(0,m1,p1sq,musq)                  
c               B0finG= B0finG2(m1,p1sq,musq)
c               return
c               else   ! p1s =0
c B0finG(0,m1,0,musq)'
               B0finGDiv= B0finG2Div(m1,p1sq,musq,id) ! B0finG2(m1,0d0,musq)= B1i0e(m02,musq,0)
               return
c            endif
        else   ! m1s=0
c            if (abs(p1sq).ge.eps) then
c B0finG(0,0,p1sq,musq)
              B0finGDiv=B0finDiv(p1sq,musq,id)
c           print*,"here",B0finGDiv,p1sq,musq,id
             return
c             else
c B0finG(0,0,0,musq)
c             B0finG= 0d0
c             return
c            endif
        endif  ! m1s=0
       endif   ! m0s =0 
      return
      else
#ifdef WITH_ONELOOP
         call avh_olo_mu_set(sqrt(musq))
         call avh_olo_onshell(eps)
         call avh_olo_b0m(rslt,p1sq,m0**2,m1**2)
         B0finGDiv = rslt(id)
#else
         call avh_error()
#endif
         return
      endif
    
      end
      

      function  C0finG(m0,m1,m2,s1,s2,s3,musq)
      IMPLICIT NONE
      complex*16 C0finG
      real*8 m0,m1,m2,s1,s2,s3,musq
      complex*16 C0fin,C0finG1M,C0finG1MPrime,C0finG21,C0finG2M,I3point
      External C0fin,C0finG1M,C0finG1MPrime,C0finG21,C0finG2M,I3point
      real*8 eps

      parameter (eps=1d-7)
      logical UseAVH, UseFF
      real*8 xpi(6)
      integer ier
      parameter (UseAVH=.false.)
      parameter (UseFF=.true.)
cccccccccccccccccccccccccccccccccccccccccccccccc
c Input
c --------------------------------------
c   C0finG1M=C0finG(0d0,0d0,m2,s1,s2,s3,musq)                             
c --------------------------------------
c --------------------------------------
c   C0finG1MPrime=C0finG(0d0,M1,0d0,s1,s2,s3,musq)                        
c --------------------------------------
c
c  C0finG2(M1,M2,M3,s1,s2,s3,musq) only s3=0
c-------------------------------------------
c
c      C0fin(0,0,0,s1,s2,s3)
c----------------------------------------------
      if(UseAVH) then
#ifdef WITH_ONELOOP
         call avh_olo_mu_set(sqrt(musq))
         call avh_olo_onshell(eps)
         call avh_olo_c0m(rslt,s1,s2,s3,m0**2,m1**2,m2**2)
         C0finG = rslt(0) 
#else
         call avh_error()
#endif
      else if (UseFF) then
        if ( (abs(m0).lt.eps) .and.
     &       (abs(m1).lt.eps) .and.
     &       (abs(m2).lt.eps) ) then
           C0finG=C0fin(s1,s2,s3,musq)
        else
          xpi(1) = m0**2
          xpi(2) = m1**2
          xpi(3) = m2**2
          xpi(4) = s1
          xpi(5) = s2
          xpi(6) = s3
          ier=0
          call vbfffxc0(C0finG,xpi,ier)
        endif
      else
       if(m0.ge.eps) then
         if (m1.ge.eps) then
            if (m2.ge.eps) then
c                    M0,M1,M3
               if((m0.eq.m1).and.(m0.eq.m2)) then 
                    C0finG=I3point(m0,s1,s2,s3)
                    else
                   C0finG=C0finG21(M0,M1,M2,s1,s2,s3,musq)
               endif
               else                   ! m2=0
c                    M0,M1,0
                C0finG=C0finG2M(m0,m1,s3,s1,s2,musq)
            endif      
         else                          ! m1=0    
            if (m2.ge.eps) then
c                   M0,0,M2
                C0finG=C0finG2M(m2,m0,s2,s3,s1,musq)
               else                     ! m2=0
c                   M0,0,0
                C0finG=C0finG1MPrime(m0,s3,s1,s2,musq)
                C0finG=C0finG1M(m0,s2,s3,s1,musq) 
            endif      
         endif
      else                            !m0 = 0
         if (m1.ge.eps) then
            if (m2.ge.eps) then
c                   0,M1,M2
                C0finG=C0finG2M(m1,m2,s1,s2,s3,musq)
             else                   ! m2=0
c                   0,M1,0
                C0finG=C0finG1MPrime(m1,s1,s2,s3,musq)
                C0finG=C0finG1M(m1,s3,s1,s2,musq)
            endif      
         else        ! m1=0    
            if (m2.ge.eps) then
c                   0,0,M2
                C0finG=C0finG1MPrime(m2,s2,s3,s1,musq)
                C0finG=C0finG1M(m2,s1,s2,s3,musq)
               else                     ! m2=0
c                   0, 0, 0
                C0finG=C0fin(s1,s2,s3,musq)
                    return
            endif      
         endif
      endif
         
      endif

      end


      function  C0finGDiv(m0,m1,m2,s1,s2,s3,musq,id)
      IMPLICIT NONE
      complex*16 C0finGDiv
      real*8 m0,m1,m2,s1,s2,s3,musq
      complex*16 C0finDiv,C0finG1MDiv,C0finG1MPrimeDiv
      complex*16 C0finG2Div,C0finG2MDiv
      External C0finDiv,C0finG1MDiv,C0finG1MPrimeDiv
      External C0finG2Div,C0finG2MDiv
      real*8 eps

      parameter (eps=1d-7)
      integer id
      logical NotAVH
      parameter (NotAVH=.true.)
cccccccccccccccccccccccccccccccccccccccccccccccc
c Input
c --------------------------------------
c   C0finG1M=C0finG(0d0,0d0,m2,s1,s2,s3,musq)                             c --------------------------------------
c --------------------------------------
c   C0finG1MPrime=C0finG(0d0,M1,0d0,s1,s2,s3,musq)                        
c --------------------------------------
c
c  C0finG2(M1,M2,M3,s1,s2,s3,musq) only s3=0
c-------------------------------------------
c
c      C0fin(0,0,0,s1,s2,s3)
c----------------------------------------------
      if(NotAVH) then
      if(m0.ge.eps) then
         if (m1.ge.eps) then
            if (m2.ge.eps) then
c                    M0,M1,M3
                 C0finGDiv=C0finG2Div(M0,M1,M2,s1,s2,s3,musq,id)
                 stop
               else                   ! m2=0
c                    M0,M1,0
                C0finGDiv=C0finG2MDiv(m0,m1,s3,s1,s2,musq,id)
            endif      
         else                          ! m1=0    
            if (m2.ge.eps) then
c                   M0,0,M2
                C0finGDiv=C0finG2MDiv(m2,m0,s2,s3,s1,musq,id)
               else                     ! m2=0
c                   M0,0,0
c                C0finGDiv=C0finG1MPrimeDiv(m0,s3,s1,s2,musq,id)
                C0finGDiv=C0finG1MDiv(m0,s2,s3,s1,musq,id) 
            endif      
         endif
      else                            !m0 = 0
         if (m1.ge.eps) then
            if (m2.ge.eps) then
c                   0,M1,M2
                C0finGDiv=C0finG2MDiv(m1,m2,s1,s2,s3,musq,id)
             else                   ! m2=0
c                   0,M1,0
c                C0finGDiv=C0finG1MPrimeDiv(m1,s1,s2,s3,musq,id)
                C0finGDiv=C0finG1MDiv(m1,s3,s1,s2,musq,id)
            endif      
         else        ! m1=0    
            if (m2.ge.eps) then
c                   0,0,M2
c                C0finGDiv=C0finG1MPrimeDiv(m2,s2,s3,s1,musq,id)
                C0finGDiv=C0finG1MDiv(m2,s1,s2,s3,musq,id)
               else                     ! m2=0
c                   0, 0, 0
                C0finGDiv=C0finDiv(s1,s2,s3,musq,id)
                    return
            endif      
         endif
      endif

      else
#ifdef WITH_ONELOOP
         call avh_olo_mu_set(sqrt(musq))
         call avh_olo_onshell(eps)
         call avh_olo_c0m(rslt,s1,s2,s3,m0**2,m1**2,m2**2)
         C0finGDiv = rslt(id) 
#else
         call avh_error()
#endif
         return
      endif

      end


      function D0finG(m0,m1,m2,m3,s12,s23,p1,p2,p3,p4,musq)
      IMPLICIT NONE
      complex*16 D0finG
      real*8 m0,m1,m2,m3,p1,p2,p3,p4,s12,s23,musq
      real*8 m0s,m1s,m2s,m3s
      complex*16 D0fin,D0t1,D0t2,D0t3,D04
      External D0fin,D0t1,D0t2,D0t3,D04

      real*8 eps
      parameter (eps=1d-6)
      logical UseAVH, UseFF
      real*8 xpi(13)
      integer ier
      parameter (UseAVH=.false.)
      parameter (UseFF=.true.)
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

      if(UseAVH) then
#ifdef WITH_ONELOOP
         m0s=m0*m0
         m1s=m1*m1
         m2s=m2*m2
         m3s=m3*m3
         call avh_olo_mu_set(sqrt(musq))
         call avh_olo_onshell(eps)
         call avh_olo_d0m(rslt,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s)
         D0finG = rslt(0) 
#else
         call avh_error()
#endif
      else if (UseFF) then
c take out mass combination which FF cannot yet handle
        if ( (abs(m0).lt.eps) .and.
     &       (abs(m1).lt.eps) .and.
     &       (abs(m2).lt.eps) .and.
     &       (abs(m3).lt.eps) ) then
c           print*, "D0_masses0"
            D0finG=D0fin(s12,s23,p1,p2,p3,p4,musq)
         else if ( (m0.eq.m1) .and.
     &             (m0.eq.m2) .and.
     &             (m0.eq.m3) ) then
           D0finG=D04(p1,p2,p3,p4,s12,s23,m0,m1,m2,m3)
        else
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
         call vbfffxd0(D0finG,xpi,ier)
        endif
      else
      if(m0.ge.eps) then
         if(m1.ge.eps) then
            if(m2.ge.eps) then
               if(m3.ge.eps) then
                  D0finG=D04(p1,p2,p3,p4,s12,s23,m0,m1,m2,m3)
               else
                  stop
               endif
            else
               if(m3.ge.eps) then
                  stop
               else
                  stop
               endif
           endif    
         else   
           if(m2.ge.eps) then
               if(m3.ge.eps) then
                  stop
               else
                  stop
               endif
            else
               if(m3.ge.eps) then
                  stop
               else
                  stop
               endif
           endif
         endif  
      else ! m0=0
         if(m1.ge.eps) then
            if(m2.ge.eps) then
               if(m3.ge.eps) then
                  stop
               else
                  stop
               endif
            else
               if(m3.ge.eps) then
                  stop
               else
                  stop
               endif
           endif    
         else   ! m1=0
           if(m2.ge.eps) then
               if(m3.ge.eps) then
c    0,0,m2,m3
c D0t3 mi4s <>0, mi3s <> 0, s1 = 0, s2 <=> 0, s3 <=> 0, s4 <=> 0, s <>0, t<>0
c Eq. 4.12
c      D02356=D0finG(0d0,0d0,M,M,s234,s345,p2sq,s34,p5sq,s16,musq)
c      D03456=D0finG(0d0,0d0,M,M,s34,s45,p3sq,p4sq,p5sq,s345,musq)
c      D01256=D0finG(0d0,0d0,M,M,s56,s16,p1sq,s234,p5sq,p6sq,musq) 
                 if(abs(p1).le.eps) then
                 m0s=m0*m0
                 m1s=m1*m1
                 m2s=m2*m2
                 m3s=m3*m3
cFC                 print*,'D0 Here D02356,D3456,D01256'
cFC                  print*,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq
                 If(abs(p2).le.eps) then
                 D0finG=D0t3(p1,p4,p3,p2,s23,s12,m1s,m0s,m3s,m2s,musq,0)  
                 else
                 D0finG=D0t3(p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq,0)  
                 endif
                 else
c-----------------------------------------------------------------
c      D01356=D0finG(0d0,0d0,M,M,s56,s345,s12,s34,p5sq,p6sq,musq)     
c      D01456=D0finG(0d0,0d0,M,M,s56,s45,s123,p4sq,p5sq,p6sq,musq)    
c      D02456=D0finG(0d0,0d0,M,M,s234,s45,s23,p4sq,p5sq,s16,musq)    
               D0finG= D04(p1,p2,p3,p4,s12,s23,m0,m1,m2,m3)
               endif
               else
c    0,0,m2,0
                  stop
               endif
            else
               if(m3.ge.eps) then
c   0,0,0,M
c D0t1 mi4s <>0, s1 = 0, s2 = 0, s3 <=> 0, s4 <=> 0, s <>0, t<>0 
c Eq. 4.29
                 m0s=m0*m0
                 m1s=m1*m1
                 m2s=m2*m2
                 m3s=m3*m3
                 if(abs(p1).ge.eps) then
                    if(abs(p2).ge.eps) then
                       if(abs(p3).ge.eps) then
                          if(abs(p4).ge.eps) then
c             p1,p2,p3,p4
                             stop
                          else
c             p1,p2,p3,0
                             stop
                          endif
                       else  ! p3=0  
                          if(abs(p4).ge.eps) then
c            p1,p2,0,p4
                             stop
                          else
c            p1,p2,0,0
                             stop
                          endif
                       endif
                    else   ! p2=0     
                       if(abs(p3).ge.eps) then
                          if(abs(p4).ge.eps) then
c           p1,0,p3,p4
                             stop
                          else
c           p1,0,p3,0 
c-----------------------------------------------------------------
c D0t2 mi4s <>0, s1 = 0, s2 <> 0 s3 <=> 0, s4 <> 0, s <>0, t<>0
c Eq. 4.12
c      D01346=D0finG(0d0,0d0,0d0,M,s123,s345,s12,p3sq,s45,p6sq,musq)
cFC                 print*,'D0 Here D01346'
cFC                  print*,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq
                 D0finG=D0t2(p2,p1,p4,p3,s12,s23,m2s,m1s,m0s,m3s,musq,0) 
                          endif
                       else   ! p3=0  
                          if(abs(p4).ge.eps) then
c           p1,0,0,p4  
c      D01345=D0finG(0d0,0d0,0d0,M,s123,s34,  s12,p3sq,p4sq,s56,musq)
c      D0(p1,0,0,p4,s,t,0,0,0,M)= D0(0,p1,p4,0,s,t,0,0,0M)
cFC                 print*,'D0 Here D01345'
cFC                  print*,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq
                 D0finG=D0t2(p2,p1,p4,p3,s12,s23,m2s,m1s,m0s,m3s,musq,0)   
                          else
c           p1,0,0,0
                             stop
                          endif
                       endif
                     endif
                  else  ! p1=0
                    if(abs(p2).ge.eps) then
                       if(abs(p3).ge.eps) then
                          if(abs(p4).ge.eps) then
c          0,p2,p3,p4
                             stop
                          else
c         0,p2,p3,0
c      D01246=D0finG(0d0,0d0,0d0,M,s123,s16,p1sq,s23,s45,p6sq,musq)
cFC                 print*,'D0 Here D01246'
cFC                  print*,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq
                 D0finG=D0t2(p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq,0)
                          endif
                       else  ! p3=0 
                          if(abs(p4).ge.eps) then
c        0,p2,0,p4 
c      D01245=D0finG(0d0,0d0,0d0,M,s123,s234,p1sq,s23,p4sq,s56,musq)
cFC                 print*,'D0 Here D01245'
cFC                  print*,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq
                D0finG=D0t2(p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq,0)
                          else
c        0,p2,0,0
                             stop
                          endif
                       endif
                    else      !p2=0
                       if(abs(p3).ge.eps) then
                          if(abs(p4).ge.eps) then
c        0,0,p3,p4
c      D01235=D0finG(0d0,0d0,0d0,M,s12,s234,p1sq,p2sq,s34,s56,musq)
c      D02346=D0finG(0d0,0d0,0d0,M,s23,s345,p2sq,p3sq,s45,s16,musq)
cFC                  print*,'D0 Here D01235,D2346'
cFC                  print*,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq
                           D0finG=D0t1(p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq,0)
cFC                Print*,"D0finG",D0finG
                          else
c       0,0,p3,0
c      D01236=D0finG(0d0,0d0,0d0,M,s12,s16,p1sq,p2sq,s345,p6sq,musq)
cFC                 print*,'D0 Here D01236'
cFC                  print*,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq
                           D0finG=D0t1(p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq,0)
cFC                 Print*,"D0finG",D0finG
                          endif
                       else  ! p3=0  
                          if(abs(p4).ge.eps) then
c       0,0,0,p4
c      D02345=D0finG(0d0,0d0,0d0,M,s23,s34,p2sq,p3sq,p4sq,s234,musq)
cFC                 print*,'D0 Here D02345'
cFC                  print*,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq
                           D0finG=D0t1(p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq,0)
cFC                 Print*,"D0finG",D0finG
                          else
c       0,0,0,0
                             stop
                          endif
                       endif
                     endif
                  endif   
                else
c      D01234=D0finG(0d0,0d0,0d0,0d0,s12,s23,p1sq,p2sq,p3sq,s123,musq)
cFC               print*,'D0 Here D01234'
cFC                 print*,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq
                  D0finG=D0fin(s12,s23,p1,p2,p3,p4,musq)
               endif
           endif
         endif
       endif


      return
      endif

      return
      end


      function D0finGDiv(m0,m1,m2,m3,s12,s23,p1,p2,p3,p4,musq,id)
      IMPLICIT NONE
      complex*16 D0finGDiv
      real*8 m0,m1,m2,m3,p1,p2,p3,p4,s12,s23,musq
      real*8 m0s,m1s,m2s,m3s
      complex*16 D0finDiv,D0t1,D0t2,D0t3,D04
      External D0finDiv,D0t1,D0t2,D0t3,D04
      real*8 eps

      parameter (eps=1d-6)
      integer id
      logical NotAVH
      parameter (NotAVH=.true.)

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

      if(NotAVH) then
      if(m0.ge.eps) then
         if(m1.ge.eps) then
            if(m2.ge.eps) then
               if(m3.ge.eps) then
                 D0finGDiv=0d0 
               else
                  stop
               endif
            else
               if(m3.ge.eps) then
                  stop
               else
                  stop
               endif
           endif    
         else   
           if(m2.ge.eps) then
               if(m3.ge.eps) then
                  stop
               else
                  stop
               endif
            else
               if(m3.ge.eps) then
                  stop
               else
                  stop
               endif
           endif
         endif  
      else ! m0=0
         if(m1.ge.eps) then
            if(m2.ge.eps) then
               if(m3.ge.eps) then
                  stop
               else
                  stop
               endif
            else
               if(m3.ge.eps) then
                  stop
               else
                  stop
               endif
           endif    
         else   ! m1=0
           if(m2.ge.eps) then
               if(m3.ge.eps) then
c    0,0,m2,m3
c D0t3 mi4s <>0, mi3s <> 0, s1 = 0, s2 <=> 0, s3 <=> 0, s4 <=> 0, s <>0, t<>0
c Eq. 4.12
c      D02356=D0finG(0d0,0d0,M,M,s234,s345,p2sq,s34,p5sq,s16,musq)
c      D03456=D0finG(0d0,0d0,M,M,s34,s45,p3sq,p4sq,p5sq,s345,musq)
c      D01256=D0finG(0d0,0d0,M,M,s56,s16,p1sq,s234,p5sq,p6sq,musq) 
                 if(abs(p1).le.eps) then
                 m0s=m0*m0
                 m1s=m1*m1
                 m2s=m2*m2
                 m3s=m3*m3
cFC                  print*,'D0 Here D02356,D3456,D01256'
cFC                 print*,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq
                 If(abs(p2).le.eps) then
                 D0finGDiv=D0t3(p1,p4,p3,p2,s23,s12,m1s,m0s,m3s,m2s,musq,id)
                 else
                 D0finGDiv=D0t3(p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq,id)
                 endif
                 else
c-----------------------------------------------------------------
c      D01356=D0finG(0d0,0d0,M,M,s56,s345,s12,s34,p5sq,p6sq,musq)     
c      D01456=D0finG(0d0,0d0,M,M,s56,s45,s123,p4sq,p5sq,p6sq,musq)    
c      D02456=D0finG(0d0,0d0,M,M,s234,s45,s23,p4sq,p5sq,s16,musq)    
               D0finGDiv= 0d0
               endif
               else
c    0,0,m2,0
                  stop
               endif
            else
               if(m3.ge.eps) then
c   0,0,0,M
c D0t1 mi4s <>0, s1 = 0, s2 = 0, s3 <=> 0, s4 <=> 0, s <>0, t<>0 
c Eq. 4.29
                 m0s=m0*m0
                 m1s=m1*m1
                 m2s=m2*m2
                 m3s=m3*m3
                 if(abs(p1).ge.eps) then
                    if(abs(p2).ge.eps) then
                       if(abs(p3).ge.eps) then
                          if(abs(p4).ge.eps) then
c             p1,p2,p3,p4
                             stop
                          else
c             p1,p2,p3,0
                             stop
                          endif
                       else  ! p3=0  
                          if(abs(p4).ge.eps) then
c            p1,p2,0,p4
                             stop
                          else
c            p1,p2,0,0
                             stop
                          endif
                       endif
                    else   ! p2=0     
                       if(abs(p3).ge.eps) then
                          if(abs(p4).ge.eps) then
c           p1,0,p3,p4
                             stop
                          else
c           p1,0,p3,0 
c-----------------------------------------------------------------
c D0t2 mi4s <>0, s1 = 0, s2 <> 0 s3 <=> 0, s4 <> 0, s <>0, t<>0
c Eq. 4.12
c      D01346=D0finG(0d0,0d0,0d0,M,s123,s345,s12,p3sq,s45,p6sq,musq)
cFC                 print*,'D0 Here D01346'
cFC                  print*,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq
                 D0finGDiv=D0t2(p2,p1,p4,p3,s12,s23,m2s,m1s,m0s,m3s,musq,id) 
                          endif
                       else   ! p3=0  
                          if(abs(p4).ge.eps) then
c           p1,0,0,p4  
c      D01345=D0finG(0d0,0d0,0d0,M,s123,s34,  s12,p3sq,p4sq,s56,musq)
c      D0(p1,0,0,p4,s,t,0,0,0,M)= D0(0,p1,p4,0,s,t,0,0,0M)
cFC                 print*,'D0 Here D01345'
cFC                  print*,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq
                 D0finGDiv=D0t2(p2,p1,p4,p3,s12,s23,m2s,m1s,m0s,m3s,musq,id)
                          else
c           p1,0,0,0
                             stop
                          endif
                       endif
                     endif
                  else  ! p1=0
                    if(abs(p2).ge.eps) then
                       if(abs(p3).ge.eps) then
                          if(abs(p4).ge.eps) then
c          0,p2,p3,p4
                             stop
                          else
c         0,p2,p3,0
c      D01246=D0finG(0d0,0d0,0d0,M,s123,s16,p1sq,s23,s45,p6sq,musq)
cFC                 print*,'D0 Here D01246'
cFC                  print*,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq
                 D0finGDiv=D0t2(p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq,id)
                          endif
                       else  ! p3=0 
                          if(abs(p4).ge.eps) then
c        0,p2,0,p4 
c      D01245=D0finG(0d0,0d0,0d0,M,s123,s234,p1sq,s23,p4sq,s56,musq)
cFC                 print*,'D0 Here D01245'
cFC                  print*,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq
                D0finGDiv=D0t2(p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq,id)
                          else
c        0,p2,0,0
                             stop
                          endif
                       endif
                    else      !p2=0
                       if(abs(p3).ge.eps) then
                          if(abs(p4).ge.eps) then
c        0,0,p3,p4
c      D01235=D0finG(0d0,0d0,0d0,M,s12,s234,p1sq,p2sq,s34,s56,musq)
c      D02346=D0finG(0d0,0d0,0d0,M,s23,s345,p2sq,p3sq,s45,s16,musq)
cFC                  print*,'D0 Here D01235,D2346'
cFC                  print*,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq
                           D0finGDiv=D0t1(p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq,id)
cFC                 Print*,"D0finG",D0finGDiv
                          else
c       0,0,p3,0
c      D01236=D0finG(0d0,0d0,0d0,M,s12,s16,p1sq,p2sq,s345,p6sq,musq)
cFC                 print*,'D0 Here D01236'
cFC                  print*,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq
                           D0finGDiv=D0t1(p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq,id)
cFC                 Print*,"D0finG",D0finGDiv
                          endif
                       else  ! p3=0  
                          if(abs(p4).ge.eps) then
c       0,0,0,p4
c      D02345=D0finG(0d0,0d0,0d0,M,s23,s34,p2sq,p3sq,p4sq,s234,musq)
cFC                 print*,'D0 Here D02345'
cFC                  print*,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq
                           D0finGDiv=D0t1(p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq,id)
cFC                 Print*,"D0finG",D0finGDiv
                          else
c       0,0,0,0
                             stop
                          endif
                       endif
                     endif
                  endif   
                else
c      D01234=D0finG(0d0,0d0,0d0,0d0,s12,s23,p1sq,p2sq,p3sq,s123,musq)
c               print*,'D0 Here D01234'
c                 print*,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s,musq
                  D0finGDiv=D0finDiv(s12,s23,p1,p2,p3,p4,musq,id)
               endif
           endif
         endif
       endif  
      return
      
      else
         
#ifdef WITH_ONELOOP
         call avh_olo_mu_set(sqrt(musq))
         call avh_olo_onshell(eps)
         m0s=m0*m0
         m1s=m1*m1
         m2s=m2*m2
         m3s=m3*m3
         call avh_olo_d0m(rslt,p1,p2,p3,p4,s12,s23,m0s,m1s,m2s,m3s)
         D0finGDiv = rslt(id) 
#else
         call avh_error()
#endif
         
      endif


      return
      end


c
c
c
c
c
c
c
c
c
c
c
c
c
c

c --------------------------------------
c   A0finG1=A0finG(0,musq)                                 
c --------------------------------------
      complex*16 function A0finG1(musq)
      implicit none
      double precision musq
      A0finG1 = (0.d0,0.d0)
      end


c --------------------------------------
c   A0finG1=A0finG(0,musq)                                 
c --------------------------------------
      complex*16 function A0finG1Div(musq,id)
      implicit none
      double precision musq
      integer id
      A0finG1Div = (0.d0,0.d0)
      end


c --------------------------------------
c   A0finG2=A0finG(M,musq)                                 
c --------------------------------------
      complex*16 function A0finG2(M,musq)
      implicit none
      double precision M,musq
      A0finG2 = M*M*(-dLog(M*M/musq)+1.d0)
      end

c --------------------------------------
c   A0finG2=A0finG(M,musq)                                 
c --------------------------------------
      complex*16 function A0finG2Div(M,musq,id)
      implicit none
      double precision M,musq
      integer id
      if(id.eq.1) then
      A0finG2Div = M*M
      else
      A0finG2Div=0d0
      endif
      end


c --------------------------------------
c   B0finG5=B0finG(M1,M2,s,musq)                                 
c --------------------------------------
      complex*16 function B0finG5(M1,M2,s,musq)
      implicit none
      real*8 pi
      complex*16 Ipi, Ieps 
      parameter(pi=3.14159265358979324d0,
     & Ipi=(0d0,3.14159265358979324d0),Ieps=(0d0,1d-38))
      double precision M1,M2,s,musq,m1sq,m2sq
      complex*16 b,g1,g2

      m1sq = M1*M1
      m2sq = M2*M2
      b = Sqrt((s-m2sq+m1sq)**2-4.d0*s*(m1sq-Ieps))
      g1 = .5d0*(s-m2sq+m1sq+b)/s
      g2 = .5d0*(s-m2sq+m1sq-b)/s

      B0finG5 = -Log((s-Ieps)/musq) 
     &        + (g1*Log((g1-1.d0)/g1)-Log(g1-1))
     &        + (g2*Log((g2-1.d0)/g2)-Log(g2-1)) + 2.d0 
      end


c --------------------------------------
c   B0finG5=B0finG(M1,M2,s,musq)                                 
c --------------------------------------
      complex*16 function B0finG5Div(M1,M2,s,musq,id)
      implicit none
      integer id
      double precision M1,M2,s,musq
      if(id.eq.1) then
      B0finG5Div = 1.d0
      else
      B0finG5Div = 0d0
      endif
      end


c --------------------------------------
c   B0finG3=B0finG(M1,M1,s,musq)                                 
c --------------------------------------
      complex*16 function B0finG4(M1,s,musq)
      implicit none
      real*8 pi
      complex*16 Ipi, Ieps 
      parameter(pi=3.14159265358979324d0,
     & Ipi=(0d0,3.14159265358979324d0),Ieps=(0d0,1d-38))
      double precision M1,s,musq,msq
      complex*16 b,g1,g2

      msq = M1*M1
      b = Sqrt(s*s-4.d0*s*(msq-Ieps))
      g1 = .5d0*(s+b)/s
      g2 = .5d0*(s-b)/s

      B0finG4 = -Log((s-Ieps)/musq) 
     &        + (g1*Log((g1-1.d0)/g1)-Log(g1-1))
     &        + (g2*Log((g2-1.d0)/g2)-Log(g2-1)) + 2.d0 
      end


c --------------------------------------
c   B0finG3=B0finG(M1,M1,s,musq)                                 
c --------------------------------------
      complex*16 function B0finG4Div(M1,s,musq,id)
      implicit none
      double precision M1,s,musq
      integer id
      if(id.eq.1) then
      B0finG4Div = 1.d0
      else
      B0finG4Div = 0d0
      endif
      end


c --------------------------------------
c   B0finG6=B0finG(M1,M1,0,musq)          
c --------------------------------------
      complex*16 function B0finG6(M1,musq)
      implicit none
      real*8 pi
      complex*16 Ipi, Ieps 
      parameter(pi=3.14159265358979324d0,
     & Ipi=(0d0,3.14159265358979324d0),Ieps=(0d0,1d-38))
      double precision M1,musq,msq
      msq = M1*M1
      B0finG6 = -dLog(msq/musq)
      end


c --------------------------------------
c   B0finG6=B0finG(M1,M1,0,musq)          
c --------------------------------------
      complex*16 function B0finG6Div(M1,musq,id)
      implicit none
      integer id
      double precision M1,musq
      if(id.eq.1) then
      B0finG6Div = 1.d0
      else
      B0finG6Div = 0d0
      endif
      end

c --------------------------------------
c   B0finG2=B0finG(M1,0,s,musq)                                 
c --------------------------------------
      complex*16 function B0finG2Div(M1,s,musq,id)
      implicit none 
      double precision M1,s,musq
      integer id
      if(id.eq.1) then
      B0finG2Div = 1.d0
      else
      B0finG2Div = 0d0
      endif 
      end


c --------------------------------------
c   C0finG2M=C0finG(0d0,m1,m2,s1,s2,s3,musq)                                 
c --------------------------------------
      function  C0finG2M(m1,m2,s1,s2,s3,musq)
      IMPLICIT NONE
      complex*16 C0finG2M
      real*8 m1,m2,s1,s2,s3,musq
      real*8 m1s,m2s
      complex*16 C2i3e,C2i2e,C0finG2
      External   C2i3e,C2i2e,C0finG2
      real*8 eps
      parameter (eps=1d-7)

c	C2i2e(s1,s2,m2s,m3s,mus)	I3d(s1,s2,0;0,mis,mis)	
c	C2i3e(s1,s2,s3,m2s,m3s,mus)	I3d(s1,s2,s3;0,mis,mis) 

c      C0456=C0finG(p4sq=0,p5sq,s45,0d0,M,M,musq)=I3d(s3,s2;0,M,M)
c      C0256=C0finG(s234,p5sq,s16;0d0,M,M;musq)=I3d(s1,s2,s3;0,mis,mis)
c      C0356=C0finG(s34,p5sq,s345;0d0,M,M,musq)=I3d(s1,s2,s3;0,mis,mis)
c      C0156=C0finG(s56,p5sq,0=p6sq;0d0,M,M,musq)=I3d(s1,s2,0;0,mis,mis)

      m1s=m1*m1
      m2s=m2*m2

c   0,M1,M2
               if(abs(s1).ge.eps) then
                  if(abs(s2).ge.eps) then
                     if(abs(s3).ge.eps) then
c                        s1,s2,s3
c                     C0finG2M= C0finG2(0d0,M1,M2,s1,s2,s3,musq)
                     C0finG2M=C2i3e(s1,s2,s3,m1s,m2s,musq,0) 
                   return
                        else
c                        s1,s2,0
c                     C0finG2M= C0finG2(0d0,M1,M2,s1,s2,s3,musq)
                     C0finG2M=C2i2e(s1,s2,m1s,m2s,musq,0)
                     endif
                 else
                     if(abs(s3).ge.eps) then
c                     s1,0,s3
c                     C0finG2M = C0finG2(M2,0d0,M1,s3,s1,s2,musq)
                     C0finG2M = C2i2e(s1,s2,m1s,m2s,musq,0)
                     else
c                     s1,0,0
!                    Print*, "Warning"
                     C0finG2M=C2i3e(s1,s2,s3,m1s,m2s,musq,0) 
                     stop
                     endif
                 endif
               else
                  if(abs(s2).ge.eps) then
                     if(abs(s3).ge.eps) then                  
c                      0,s2,s3
c                     C0finG2M= C0finG2(M1,M2,0d0,s2,s3,s1,musq)
c                     C0finG2M=C2i3e(s1,s2,s3,m1s,m2s,musq,0) 
                     C0finG2M=C2i2e(s3,s2,m2s,m1s,musq,0)
                     else
c                      0,s2,0
!                    Print*, "Warning"
                     C0finG2M=C2i3e(s1,s2,s3,m1s,m2s,musq,0) 
                     stop
                     endif
                 else
                     if(abs(s3).ge.eps) then
c                      0,0,s3
!                    Print*, "Warning"
                     C0finG2M=C2i3e(s1,s2,s3,m1s,m2s,musq,0) 
                     stop
                     else
c                     0,0,0                  
!                    Print*, "Warning"
!                     C0finG=C2i3e(s1,s2,s3,m1s,m2s,musq,0) 
                     stop
                     endif       
                  endif   
               endif 
      
  
          end


c --------------------------------------
c   C0finG2M=C0finG(0d0,m1,m2,s1,s2,s3,musq)                                 
c --------------------------------------
      function  C0finG2MDiv(m1,m2,s1,s2,s3,musq,id)
      IMPLICIT NONE
      complex*16 C0finG2MDiv
      real*8 m1,m2,s1,s2,s3,musq
      real*8 m1s,m2s
      complex*16 C2i3e,C2i2e,C0finG2Div
      External   C2i3e,C2i2e,C0finG2Div
      real*8 eps
      parameter (eps=1d-7)
      integer id

c	C2i2e(s1,s2,m2s,m3s,mus)	I3d(s1,s2,0;0,mis,mis)	
c	C2i3e(s1,s2,s3,m2s,m3s,mus)	I3d(s1,s2,s3;0,mis,mis) 

c      C0456=C0finG(p4sq=0,p5sq,s45,0d0,M,M,musq)=I3d(s3,s2;0,M,M)
c      C0256=C0finG(s234,p5sq,s16;0d0,M,M;musq)=I3d(s1,s2,s3;0,mis,mis)
c      C0356=C0finG(s34,p5sq,s345;0d0,M,M,musq)=I3d(s1,s2,s3;0,mis,mis)
c      C0156=C0finG(s56,p5sq,0=p6sq;0d0,M,M,musq)=I3d(s1,s2,0;0,mis,mis)

      m1s=m1*m1
      m2s=m2*m2

c   0,M1,M2
               if(abs(s1).ge.eps) then
                  if(abs(s2).ge.eps) then
                     if(abs(s3).ge.eps) then
c                        s1,s2,s3
c                     C0finG2MDiv= C0finG2Div(0d0,M1,M2,s1,s2,s3,musq,id)
                     C0finG2MDiv=C2i3e(s1,s2,s3,m1s,m2s,musq,id) 
                   return
                        else
c                        s1,s2,0
c                     C0finG2MDiv= C0finG2Div(0d0,M1,M2,s1,s2,s3,musq,id)
                     C0finG2MDiv=C2i2e(s1,s2,m1s,m2s,musq,id)
                     endif
                 else
                     if(abs(s3).ge.eps) then
c                     s1,0,s3
c                     C0finG2MDiv = C0finG2Div(M2,0d0,M1,s3,s1,s2,musq,id)
                     C0finG2MDiv = C2i2e(s1,s2,m1s,m2s,musq,id)
                     else
c                     s1,0,0
!                    Print*, "Warning"
                     C0finG2MDiv=C2i3e(s1,s2,s3,m1s,m2s,musq,id) 
                     stop
                     endif
                 endif
               else
                  if(abs(s2).ge.eps) then
                     if(abs(s3).ge.eps) then                  
c                      0,s2,s3
c                     C0finG2MDiv= C0finG2Div(M1,M2,0d0,s2,s3,s1,musq,id)
c                     C0finG2MDiv=C2i3e(s1,s2,s3,m1s,m2s,musq,id) 
                     C0finG2MDiv=C2i2e(s3,s2,m2s,m1s,musq,id)
                     else
c                      0,s2,0
!                    Print*, "Warning"
                     C0finG2MDiv=C2i3e(s1,s2,s3,m1s,m2s,musq,id) 
                     stop
                     endif
                 else
                     if(abs(s3).ge.eps) then
c                      0,0,s3
!                    Print*, "Warning"
                     C0finG2MDiv=C2i3e(s1,s2,s3,m1s,m2s,musq,id) 
                     stop
                     else
c                     0,0,0                  
!                    Print*, "Warning"
!                     C0finG=C2i3e(s1,s2,s3,m1s,m2s,musq,0) 
                     stop
                     endif       
                  endif   
               endif 
      
  
          end


c --------------------------------------
c   C0finG1M=C0finG(0d0,0d0,m2,s1,s2,s3,musq)                                 
c --------------------------------------
      function  C0finG1M(m2,s1,s2,s3,musq)
      IMPLICIT NONE
      complex*16 C0finG1M
      real*8 m2,s1,s2,s3,musq
      real*8 m2s
      complex*16 C0fin,C1i1e,C0finG2,C1d2e,C1i2e
      External C0fin,C1i1e,C0finG2,C1d2e,C1i2e
      real*8 eps
      parameter (eps=1d-7)
       m2s=m2*m2

c      C0126=C0finG(0d0,0d0,M,0,s16,0,musq)=C1i1e(s16,Ms2,mus)
c      C0345=C0finG(0d0,0d0,M,0,0,s34,musq)=C1i1e(s34,Ms2,mus)

c      C0125=C0finG(0d0,0d0,M,0,s234,s56,musq)=C1i2e(s234,s56,M,mus)	
c      C0235=C0finG(0d0,0d0,M,0,s34,s234,musq)
c      C0236=C0finG(0d0,0d0,M,0,s345,s16,musq)
c      C0346=C0finG(0d0,0d0,M,0,s45,s345,musq)

c      C0145=C0finG(0d0,0d0,M,s123,0,s56,musq)=C1d2e(s1,s3,m2s,musq
c      C0245=C0finG(0d0,0d0,M,s23,0,s234,musq)

c      C0136=C0finG(0d0,0d0,M,s12,s345,0,musq)=C1d2e(s1,s3,m2s,musq
c      C0146=C0finG(0d0,0d0,M,s123,s45,0,musq)

c      C0135=C0finG(0d0,0d0,M,s12,s34,s56,musq)= f.o.k
c      C0246=C0finG(0d0,0d0,M,s23,s45,s16,musq)

c	C1i1e(s,mis,mus)		I3d(s,0,0;0,mis,0)
c	C1i2e(s1,s2,mis,mus)		I3d(s1,s2,0;0,mis,0)
c	C1d2e(s1,s2,mis,mus)		I3d(s1,s2,0;0,0,mis)

               if(abs(s1).ge.eps) then
                  if(abs(s2).ge.eps) then
                     if(abs(s3).ge.eps) then
c           0,0,M2;s1,s2,s3
                    C0finG1M=C0finG2(0d0,0d0,M2,s1,s2,s3,musq)
                   return
                        else
c             0,0,M2;s1,s2,0
                    C0finG1M=C1d2e(s1,s2,m2s,musq,0)
c                    C0finG1M=C0finG2(0d0,0d0,M2,s1,s2,s3,musq)
                     endif
                 else
                     if(abs(s3).ge.eps) then
c              0,0,M2;s1,0,s3
                    C0finG1M=C1d2e(s1,s3,m2s,musq,0)
                     else
c              0,0,M;s1,0,0
!                    Print*, "Warning"
                     stop
c                    C0finG=C1i1e(s1,m2s,musq)
                     endif
                 endif
               else
                  if(abs(s2).ge.eps) then
                     if(abs(s3).ge.eps) then                  
c              0,0,M; 0,s2,s3
                     C0finG1M=C1i2e(s3,s2,m2s,musq,0)
                     else
c              0,0,M; 0,s2,0
                     C0finG1M=C1i1e(s2,m2s,musq,0)
                     endif
                 else
                     if(abs(s3).ge.eps) then
c               0,0,M; 0,0,s3
                    C0finG1M=C1i1e(s3,m2s,musq,0) 
                    else
c              0,0,M; 0,0,0                  
!                    Print*, "Warning"
                    C0finG1M=0d0
                    stop
                    endif       
                  endif   
               endif 

          end


c --------------------------------------
c   C0finG1M=C0finG(0d0,0d0,m2,s1,s2,s3,musq)                                 
c --------------------------------------
      function  C0finG1MDiv(m2,s1,s2,s3,musq,id)
      IMPLICIT NONE
      complex*16 C0finG1MDiv
      real*8 m2,s1,s2,s3,musq
      real*8 m2s
      complex*16 C1i1e,C0finG2Div,C1d2e,C1i2e
      External C1i1e,C0finG2Div,C1d2e,C1i2e
      real*8 eps
      parameter (eps=1d-7)
      integer id
       m2s=m2*m2

c      C0126=C0finG(0d0,0d0,M,0,s16,0,musq)=C1i1e(s16,Ms2,mus)
c      C0345=C0finG(0d0,0d0,M,0,0,s34,musq)=C1i1e(s34,Ms2,mus)

c      C0125=C0finG(0d0,0d0,M,0,s234,s56,musq)=C1i2e(s234,s56,M,mus)	
c      C0235=C0finG(0d0,0d0,M,0,s34,s234,musq)
c      C0236=C0finG(0d0,0d0,M,0,s345,s16,musq)
c      C0346=C0finG(0d0,0d0,M,0,s45,s345,musq)

c      C0145=C0finG(0d0,0d0,M,s123,0,s56,musq)=C1d2e(s1,s3,m2s,musq
c      C0245=C0finG(0d0,0d0,M,s23,0,s234,musq)

c      C0136=C0finG(0d0,0d0,M,s12,s345,0,musq)=C1d2e(s1,s3,m2s,musq
c      C0146=C0finG(0d0,0d0,M,s123,s45,0,musq)

c      C0135=C0finG(0d0,0d0,M,s12,s34,s56,musq)= f.o.k
c      C0246=C0finG(0d0,0d0,M,s23,s45,s16,musq)

c	C1i1e(s,mis,mus)		I3d(s,0,0;0,mis,0)
c	C1i2e(s1,s2,mis,mus)		I3d(s1,s2,0;0,mis,0)
c	C1d2e(s1,s2,mis,mus)		I3d(s1,s2,0;0,0,mis)

               if(abs(s1).ge.eps) then
                  if(abs(s2).ge.eps) then
                     if(abs(s3).ge.eps) then
c           0,0,M2;s1,s2,s3
                    C0finG1MDiv=C0finG2Div(0d0,0d0,M2,s1,s2,s3,musq,id)
                   return
                        else
c             0,0,M2;s1,s2,0
                    C0finG1MDiv=C1d2e(s1,s2,m2s,musq,id)
c                    C0finG1MDiv=C0finG2Div(0d0,0d0,M2,s1,s2,s3,musq,id)
                     endif
                 else
                     if(abs(s3).ge.eps) then
c              0,0,M2;s1,0,s3
                    C0finG1MDiv=C1d2e(s1,s3,m2s,musq,id)
                     else
c              0,0,M;s1,0,0
!                    Print*, "Warning"
                     stop
c                    C0finG=C1i1e(s1,m2s,musq)
                     endif
                 endif
               else
                  if(abs(s2).ge.eps) then
                     if(abs(s3).ge.eps) then                  
c              0,0,M; 0,s2,s3
                     C0finG1MDiv=C1i2e(s3,s2,m2s,musq,id)
                     else
c              0,0,M; 0,s2,0
                     C0finG1MDiv=C1i1e(s2,m2s,musq,id)
                     endif
                 else
                     if(abs(s3).ge.eps) then
c               0,0,M; 0,0,s3
                    C0finG1MDiv=C1i1e(s3,m2s,musq,id) 
                    else
c              0,0,M; 0,0,0                  
!                    Print*, "Warning"
                    C0finG1MDiv=0d0
                    stop
                    endif       
                  endif   
               endif 

          end


c --------------------------------------
c   C0finG1M=C0finG(0d0,M1,0d0,s1,s2,s3,musq)                                 
c --------------------------------------
      function  C0finG1MPrime(m1,s1,s2,s3,musq)
      IMPLICIT NONE
      complex*16 C0finG1MPrime
      real*8 m1,s1,s2,s3,musq
      complex*16 C0fin,C0finG2,C1i1e,C1i2e,C1d2e
      real*8 m1s
      External C0fin,C0finG2,C1i1e,C1i2e,C1d2e
      real*8 eps
      parameter (eps=1d-7)

      m1s=m1*m1
c	C1i1e(s,mis,mus)		I3d(s,0,0;0,mis,0)
c	C1i2e(s1,s2,mis,mus)		I3d(s1,s2,0;0,mis,0)
c	C1d2e(s1,s2,mis,mus)		I3d(s1,s2,0;0,0,mis)

               if(abs(s1).ge.eps) then
                  if(abs(s2).ge.eps) then
                     if(abs(s3).ge.eps) then
c           0,M1,0;s1,s2,s3
                    C0finG1MPrime=C0finG2(0d0,M1,0d0,s1,s2,s3,musq)
                   return
                        else
c             0,M1,0;s1,s2,0
                    C0finG1MPrime=C1i2e(s1,s2,m1s,musq,0)
c                   C0finG1MPrime=C0finG2(0d0,M1,0d0,s1,s2,s3,musq)
                     endif
                 else
                     if(abs(s3).ge.eps) then
c              0,M1,0;s1,0,s3
                    C0finG1MPrime=C1d2e(s3,s1,m1s,musq,0)
c                   C0finG1MPrime=C0finG2(0d0,0d0,M1,s3,s1,s2,musq)
                     else
c              0,M1,0;s1,0,0
                     C0finG1MPrime=C1i1e(s1,m1s,musq,0)
                     endif
                 endif
               else
                  if(abs(s2).ge.eps) then
                     if(abs(s3).ge.eps) then                  
c              0,M1,0; 0,s2,s3
                     C0finG1MPrime=C1d2e(s3,s2,m1s,musq,0)
c                    C0finG1MPrime=C0finG2(M1,0d0,0d0,s2,s3,s1,musq)
                     else
c              0,M1,0; 0,s2,0
                     C0finG1MPrime=C1i1e(s2,m1s,musq,0)
                     endif
                 else
                     if(abs(s3).ge.eps) then
c               0,M1,0; 0,0,s3
                    Print*, "Warning"
!                    C0finG=C2i3e(s1,s2,s3,m1s,0d0,musq)
!                    C0finG=C2i2e(s3,0d0,m1s,m1s,mus)
                    stop
                    else
c              0,M1,0; 0,0,0                  
!                    Print*, "Warning"
                    stop
                    endif       
                  endif   
               endif 

          end


c --------------------------------------
c   C0finG1M=C0finG(0d0,M1,0d0,s1,s2,s3,musq)                                 
c --------------------------------------
      function  C0finG1MPrimeDiv(m1,s1,s2,s3,musq,id)
      IMPLICIT NONE
      complex*16 C0finG1MPrimeDiv
      real*8 m1,s1,s2,s3,musq
      complex*16 C0finG2Div,C1i1e,C1i2e,C1d2e
      real*8 m1s
      External C0finG2Div,C1i1e,C1i2e,C1d2e
      real*8 eps
      parameter (eps=1d-7)
      integer id

      m1s=m1*m1
c	C1i1e(s,mis,mus)		I3d(s,0,0;0,mis,0)
c	C1i2e(s1,s2,mis,mus)		I3d(s1,s2,0;0,mis,0)
c	C1d2e(s1,s2,mis,mus)		I3d(s1,s2,0;0,0,mis)

               if(abs(s1).ge.eps) then
                  if(abs(s2).ge.eps) then
                     if(abs(s3).ge.eps) then
c           0,M1,0;s1,s2,s3
                    C0finG1MPrimeDiv=C0finG2Div(0d0,M1,0d0,s1,s2,s3,musq,id)
                   return
                        else
c             0,M1,0;s1,s2,0
                    C0finG1MPrimeDiv=C1i2e(s1,s2,m1s,musq,id)
c                   C0finG1MPrimeDiv=C0finG2Div(0d0,M1,0d0,s1,s2,s3,musq,id)
                     endif
                 else
                     if(abs(s3).ge.eps) then
c              0,M1,0;s1,0,s3
                    C0finG1MPrimeDiv=C1d2e(s3,s1,m1s,musq,id)
c                   C0finG1MPrimeDiv=C0finG2Div(0d0,0d0,M1,s3,s1,s2,musq,id)
                     else
c              0,M1,0;s1,0,0
                     C0finG1MPrimeDiv=C1i1e(s1,m1s,musq,id)
                     endif
                 endif
               else
                  if(abs(s2).ge.eps) then
                     if(abs(s3).ge.eps) then                  
c              0,M1,0; 0,s2,s3
                     C0finG1MPrimeDiv=C1d2e(s3,s2,m1s,musq,id)
c                    C0finG1MPrimeDiv=C0finG2Div(M1,0d0,0d0,s2,s3,s1,musq,id)
                     else
c              0,M1,0; 0,s2,0
                     C0finG1MPrimeDiv=C1i1e(s2,m1s,musq,id)
                     endif
                 else
                     if(abs(s3).ge.eps) then
c               0,M1,0; 0,0,s3
                    Print*, "Warning"
!                    C0finG=C2i3e(s1,s2,s3,m1s,0d0,musq)
!                    C0finG=C2i2e(s3,0d0,m1s,m1s,mus)
                    stop
                    else
c              0,M1,0; 0,0,0                  
!                    Print*, "Warning"
                    stop
                    endif       
                  endif   
               endif 

          end


c --------------------------------------
c   C0finG2=C0finG full finite
c   only s3 can be zero
c --------------------------------------
      complex*16 function C0finG2Div(M1,M2,M3,s1,s2,s3,musq,id)
      implicit none
       double precision M1,M2,M3,s1,s2,s3,musq
 
      integer id
 
      C0finG2Div =0d0
      end

c --------------------------------------
c   C0finG2=C0finG full finite
c   only s3 can be zero
c --------------------------------------
      complex*16 function C0finG21(M1,M2,M3,s1,s2,s3,musq)
      implicit none
      real*8 s1,s2,s3,m1,m2,m3,musq,eps,m1sq,m2sq,m3sq
      complex*16 C0finG2,C3i3e
      external C0finG2,C3i3e
      parameter(eps=1d-6)
      
!
               if(abs(s1).ge.eps) then
                  if(abs(s2).ge.eps) then
c           M1,M2,M3;s1,s2,s3
                     m1sq=m1*m1
                     m2sq=m2*m2
                     m3sq=m3*m3
c                    C0finG21=C0finG2(M1,M2,M3,s1,s2,s3,musq)
                    C0finG21=C3i3e(s1,s2,s3,m1sq,m2sq,m3sq,musq,0) 
                   return
                 else
                     if(abs(s3).ge.eps) then
c              M1,M2,M3;s1,0,s3
c                    C0finG21=C0finG2(M2,M1,M3,s1,s3,s2,musq)
                     m1sq=m1*m1
                     m2sq=m2*m2
                     m3sq=m3*m3
                    C0finG21=C3i3e(s1,s3,s2,m2sq,m1sq,m3sq,musq,0)
                     else
c              M1,M2,M3;s1,0,0
                    Print*, "Warning"
                     stop
c                    C0finG=C1i1e(s1,m2s,musq)
                     endif
                 endif
               else
                  if(abs(s2).ge.eps) then
                     if(abs(s3).ge.eps) then                  
c              M1,M2,M3; 0,s2,s3
c                     C0finG21=C0finG2(M2,M3,M1,s2,s3,s1,musq)
                     m1sq=m1*m1
                     m2sq=m2*m2
                     m3sq=m3*m3
                     C0finG21=C3i3e(s2,s3,s1,m2sq,m3sq,m1sq,musq,0)
                     else
c              M1,M2,M3; 0,s2,0
c                     C0finG21=C1i1e(s2,m2s,musq,id)
                    Print*, "Warning"
                     stop
                     endif
                 else
                     if(abs(s3).ge.eps) then
c               M1,M2,M3; 0,0,s3
                    C0finG21=0d0
                    Print*, "Warning"
                   stop
                    else
c              M1,M2,M3;0,0,0                  
                    Print*, "Warning"
                    C0finG21=0d0
                    stop
                    endif       
                  endif   
               endif 

          end     


ccccc COMPLEX MASSES CCCCCCCCCCCCCCCCCCCCCCCCCCCc


      function A0finG_c(m,musq)
      IMPLICIT NONE
      complex*16 A0finG_c
      real*8 eps
      parameter (eps=1d-7)
      real*8 musq
      complex*16 rslt(0:2)
      logical NotAVH
      parameter (NotAVH=.true.)
      complex*16 m ! m is the mass squared

      if(NotAVH) then
        Print*, "Not program"       
      else
#ifdef WITH_ONELOOP
         call avh_olo_mu_set(sqrt(musq))
         call avh_olo_onshell(eps*eps)
         call avh_olo_a0c(rslt,m)
        
#else
         call avh_error()
#endif
         A0finG_c = rslt(0)
      endif

      return
      end


      function A0finGDiv_c(m,musq,id)
      IMPLICIT NONE
      complex*16 A0finGDiv_c
      real*8 eps
      parameter (eps=1d-7)
      real*8 musq
      complex*16 rslt(0:2)
      integer id
      logical NotAVH
      parameter (NotAVH=.true.)
      complex*16 m ! m is the mass squared

      if(NotAVH) then
         Print*, "Not program"       
      else
#ifdef WITH_ONELOOP
         call avh_olo_mu_set(sqrt(musq))
         call avh_olo_onshell(eps*eps)
         call avh_olo_a0c(rslt,m)
#else
         call avh_error()
#endif
         A0finGDiv_c = rslt(id)
         
      endif
         return
      end


      function B0finG_c(m0,m1,p1,musq)
      implicit none
      complex*16 B0finG_c
      complex*16 m0,m1 ! mass squared
      real*8 eps,musq,p1

      complex*16 rslt(0:2)
      parameter (eps=1d-7)
      logical NotAVH
      parameter (NotAVH=.true.)

  
      if(NotAVH) then
         Print*, "Not program"       
      else
#ifdef WITH_ONELOOP
         p1c = dcmplx(p1,0d0)
         call avh_olo_mu_set(sqrt(musq))
         call avh_olo_onshell(eps*eps)
         call avh_olo_b0c(rslt,p1c,m0,m1)
#else
         call avh_error()
#endif
         B0finG_c = rslt(0)
         return
      endif
    
      end


      function B0finGDiv_c(m0,m1,p1,musq,id)
      implicit none
      complex*16 B0finGDiv_c
      complex*16 m0,m1 ! mass squared
      real*8 eps,musq,p1

      complex*16 rslt(0:2)
      parameter (eps=1d-7)
      integer id
      logical NotAVH
      parameter (NotAVH=.true.)


      if(NotAVH) then
        Print*, "Not program"     
      else
#ifdef WITH_ONELOOP
         p1c = dcmplx(p1,0d0)
         
         call avh_olo_mu_set(sqrt(musq))
         call avh_olo_onshell(eps*eps)
         call avh_olo_b0c(rslt,p1c,m0,m1)
#else
         call avh_error()
#endif
         B0finGDiv_c = rslt(id)
         return
      endif
    
      end
      

      function  C0finG_c(m0,m1,m2,s1,s2,s3,musq)
      IMPLICIT NONE
      complex*16 C0finG_c
      complex*16  m0,m1,m2 ! mass squared
      real*8 s1,s2,s3,musq
      complex*16 s1c,s2c,s3c
      real*8 eps
      complex*16 rslt(0:2)
      parameter (eps=1d-7)
      logical NotAVH
      parameter (NotAVH=.true.)

      if(NotAVH) then
         Print*, "Not program"     
      else
         s1c = dcmplx(s1,0d0)
         s2c = dcmplx(s2,0d0)
         s3c = dcmplx(s3,0d0)
 
#ifdef WITH_ONELOOP
         call avh_olo_mu_set(sqrt(musq))
         call avh_olo_onshell(eps)
         call avh_olo_c0c(rslt,s1c,s2c,s3c,m0,m1,m2)
#else
         call avh_error()
#endif
         C0finG_c = rslt(0) 
         return
      endif

      end


      function  C0finGDiv_c(m0,m1,m2,s1,s2,s3,musq,id)
      IMPLICIT NONE
      complex*16 C0finGDiv_c
      complex*16 m0,m1,m2 ! mass squared
      real*8 s1,s2,s3,musq
      complex*16 s1c,s2c,s3c
      real*8 eps
      complex*16 rslt(0:3)
      parameter (eps=1d-7)
      integer id
      logical NotAVH
      parameter (NotAVH=.true.)

      if(NotAVH) then

         Print*, "Not program"     
      else
         s1c = dcmplx(s1,0d0)
         s2c = dcmplx(s2,0d0)
         s3c = dcmplx(s3,0d0)

#ifdef WITH_ONELOOP
         call avh_olo_mu_set(sqrt(musq))
         call avh_olo_onshell(eps)
         call avh_olo_c0c(rslt,s1c,s2c,s3c,m0,m1,m2)
#else
         call avh_error()
#endif
         C0finGDiv_c = rslt(id) 
         return
      endif

      end


      function D0finG_c(m0,m1,m2,m3,s12,s23,p1,p2,p3,p4,musq)
      IMPLICIT NONE
      complex*16 D0finG_c
      complex*16 m0,m1,m2,m3 ! mass squared
      
      real*8 p1,p2,p3,p4,s12,s23,musq
      complex*16 p1c,p2c,p3c,p4c,s12c,s23c
      complex*16 rslt(0:2)
      real*8 eps
      parameter (eps=1d-6)
      logical NotAVH
      parameter (NotAVH=.true.)


      if(NotAVH) then
         Print*, "Not program"        
      else
         
         p1c = dcmplx(p1,0d0)
         p2c = dcmplx(p2,0d0)
         p3c = dcmplx(p3,0d0)
         p4c = dcmplx(p4,0d0)
         s12c = dcmplx(s12,0d0)
         s23c = dcmplx(s23,0d0)

#ifdef WITH_ONELOOP
         call avh_olo_mu_set(sqrt(musq))
         call avh_olo_onshell(eps*eps)
         call avh_olo_d0c(rslt,p1c,p2c,p3c,p4c,s12c,s23c,m0,m1,m2,m3)
#else
         call avh_error()
#endif
         D0finG_c = rslt(0) 
      endif

      return
      end


      function D0finGDiv_c(m0,m1,m2,m3,s12,s23,p1,p2,p3,p4,musq,id)
      IMPLICIT NONE
      complex*16 D0finGDiv_c
      complex*16 m0,m1,m2,m3
     
      real*8 p1,p2,p3,p4,s12,s23,musq
      complex*16 p1c,p2c,p3c,p4c,s12c,s23c
      real*8 eps
      complex*16 rslt(0:2)
      parameter (eps=1d-6)
      integer id
      logical NotAVH
      parameter (NotAVH=.true.)


      if(NotAVH) then

         Print*, "Not program"     
      
      else

         p1c = dcmplx(p1,0d0)
         p2c = dcmplx(p2,0d0)
         p3c = dcmplx(p3,0d0)
         p4c = dcmplx(p4,0d0)
         s12c = dcmplx(s12,0d0)
         s23c = dcmplx(s23,0d0)
         
#ifdef WITH_ONELOOP
         call avh_olo_mu_set(sqrt(musq))
         call avh_olo_onshell(eps*eps)
        
         call avh_olo_d0c(rslt,p1c,p2c,p3c,p4c,s12c,s23c,m0,m1,m2,m3)
#else
         call avh_error()
#endif
         D0finGDiv_c = rslt(id) 
         
      endif


      return
      end

#ifndef WITH_ONELOOP
      subroutine avh_error()
      implicit none
        print *, "Error: OneLoop not linked, but requested"
        print *, "Change type of loop functions in B0C0D01MASG.F"
        stop
      end
#endif

      subroutine vbfff_diverror()
      implicit none
        print *, "Error: divergent part for FF routines not available"
      end

