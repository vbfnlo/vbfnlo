c  A.Denner and S. Dittmair,"One-loop 4-point integrals", arXiv:1005.2076
c
c D0t1 mi4s <>0, s1 = 0, s2 = 0, s3 <=> 0, s4 <=> 0, s <>0, t<>0 
c Eq. 4.29
c D0t2 mi4s <>0, s1 = 0, s2 <> 0 s3 <=> 0, s4 <> 0, s <>0, t<>0
c Eq. 4.12
c D0t3 mi4s <>0, mi3s <> 0, s1 = 0, s2 <=> 0, s3 <=> 0, s4 <=> 0, s <>0, t<>0
c Eq. 4.12
c
c D04 for mi3s <>0, mi4s <> 0, s<>0,t<>0, s1 <>0, s2=0, s3<>0, s4=0 is finite
c
c-----------------------------------------------------------------------
      complex*16 function D0t1(s1,s2,s3,s4,s,t,mi1s,
     $     mi2s,mi3s,mi4s,mus,id)

      implicit none
      integer id
      double precision s1,s2,s3,s4,s,t,mi1s,mi2s,
     $     mi3s,mi4s,delta,mus
      complex*16 ieps,sbar,tbar,s3bar,s4bar
      parameter (ieps=(0d0,1d-16))
      double precision M1,M2,M3,M4,m1s,m2s,m3s,m4s
      double complex tri0m,tri2m,tri3m,D0fin,D0ms ! mass regulated D0
      double complex tri0dreg(-2:0), tri2dreg(-2:0), tri3dreg(-2:0)
      double complex D0dreg(-2:0)
      integer ep
      double complex D04,C1i2e,C0i1e
      external D04
      complex*16 VLI2
      external VLI2
      parameter(delta=1.0d-6)
c
c     delta is the mass regulator 
c
c     1.) mi4s <>0, s3 <=> 0, s4 <=> 0, s <>0, t<>0
c     mi4 must be nonzero 
c      
      m1s = delta**2
      m2s = delta**2
      m3s = delta**2 
      m4s = mi4s


      if(id.eq.0) then
      M1 = sqrt(m1s)
      M2 = sqrt(m2s)
      M3 = sqrt(m3s)
      M4 = sqrt(m4s)
c      qlI3(p1sq,p2sq,p3sq,m1sq,m2sq,m3sq,musq,ep)

c      tri0m = qlI3(s3,t,0.0d0,m2s,m4s,m3s,mus,0)
      
      tbar = t + ieps
      s3bar = s3 + ieps
      s4bar = s4 + ieps
      
      tri0m = 1d0/(s3 - t)*(log((m4s-s3bar)/(delta*M4))**2 - 
     $     log((m4s-tbar)/(delta*M4))**2 + 
     $     VLI2(s3bar/m4s) - VLI2(tbar/m4s))
     
c      tri2m = qlI3(t,s4,0d0,m1s,m4s,m2s,mus,0)

      tri2m = 1d0/(t - s4)*(log((m4s-tbar)/(delta*M4))**2 - 
     $     log((m4s-s4bar)/(delta*M4))**2 + 
     $     VLI2(tbar/m4s) - VLI2(s4bar/m4s))
      
c      tri3m = qlI3(0.0d0,0.0d0,s,m1s,m2s,m3s,mus,0)

      sbar = s+ieps
      tri3m = 1/2d0/s*log(-m1s/sbar)**2
     
c     D0 sigular piece in mass reg.

      D0ms = ((t-s3)*tri0m + (t-s4) * tri2m)/s/(t-m4s) + tri3m/(t-m4s)
      
      if (s.gt.0d0) then
c      D0fin = D04(s1,s2,s3,s4,s,t,M1,M2,M3,M4) - D0ms
      D0fin = D04(s1,s4,s3,s2,t,s,M2,M1,M4,M3) - D0ms
c      D0fin = D04(s2,s1,s4,s3,s,t,M3,M2,M1,M4) - D0ms
      else
      D0fin = D04(s2,s3,s4,s1,t,s,M2,M3,M4,M1) - D0ms
      endif


c FC %      PRINT*, 'fc d0MS', d0MS
c FC %      PRINT*, 'fc d0FIN', d0FIN
      endif

c     we have the regularization ind. finite piece now
c     add divergent piece in the dim. reg. 
c
      ep = id
      tri0dreg(ep) = C1i2e(s3,t,m4s,mus,ep) !qlI3(t,0.0d0,s3,m4s,0.0d0,0.0d0,mus,ep) ! 2 mass ext 1 mass int
      tri2dreg(ep) = C1i2e(t,s4,m4s,mus,ep) !qlI3(0d0,t,s4,0.0d0,0.0d0,m4s,mus,ep) ! 2 mass ex 1 mass int
      tri3dreg(ep) = C0i1e(s,mus,ep) !qlI3(0.0d0,0.0d0,s,0.0d0,0.0d0,0.0d0,mus,ep)! 1 mass
c     D0 singular piece in dim reg.
      D0dreg(ep) = ((t-s3)*tri0dreg(ep) + 
     $     (t-s4) * tri2dreg(ep))/s/(t-m4s) 
     $     + tri3dreg(ep)/(t-m4s)

      if(id.eq.0) then
         D0t1 = D0fin + D0dreg(0)
         return
      elseif(id.eq.-1) then
         D0t1 = D0dreg(-1)
         return
      else
         D0t1 = D0dreg(-2)
         return
      endif
         
      end
c-------------------------
      complex*16 function D0t2(s1,s2,s3,s4,s,t,mi1s,
     $     mi2s,mi3s,mi4s,mus,id)

      implicit none
      integer id
      double precision s1,s2,s3,s4,s,t,mi1s,mi2s,
     $     mi3s,mi4s,delta,mus
      complex*16 ieps,sbar,tbar,s3bar,s4bar,s2bar
      parameter (ieps=(0d0,1d-16))
      double precision M1,M2,M3,M4,m1s,m2s,m3s,m4s
      double complex tri2m,tri3m,D0fin,D0ms ! mass regulated D0
      double complex tri2dreg(-2:0), tri3dreg(-2:0)
      double complex D0dreg(-2:0)
      integer ep
      double complex D04,C1i2e,C0i2e
      external D04
      complex*16 VLI2
      external VLI2
      parameter(delta=1.0d-6)
c
c     delta is the mass regulator 
c
c     2.) mi4s <>0, s1 = 0, s2 <> 0 s3 <=> 0, s4 <> 0, s <>0, t<>0
c     mi4 must be nonzero mi1=mi2=mi3=0
c     mi3 = mi1 = mi2 = 0 
      m1s = delta**2
      m2s = delta**2
      m3s = delta**2 
      m4s = mi4s


      if(id.eq.0) then
      M1 = sqrt(m1s)
      M2 = sqrt(m2s)
      M3 = sqrt(m3s)
      M4 = sqrt(m4s)
c      qlI3(p1sq,p2sq,p3sq,m1sq,m2sq,m3sq,musq,ep)

      tbar = t + ieps
      s3bar = s3 + ieps
      s4bar = s4 + ieps
      s2bar = s2 + ieps
               
c      tri2m = qlI3(t,s4,0d0,m1s,m4s,m2s,mus,0)

      tri2m = 1d0/(t - s4)*(log((m4s-tbar)/(delta*M4))**2 - 
     $     log((m4s-s4bar)/(delta*M4))**2 + 
     $     VLI2(tbar/m4s) - VLI2(s4bar/m4s))
      
c      tri3m = qlI3(0.0d0,s2,s,m1s,m2s,m3s,mus,0)

      sbar = s+ieps

      tri3m = 1d0/(s - s2)*(log((m1s-sbar)/(delta*M1))**2 - 
     $     log((m1s-s2bar)/(delta*M1))**2 + 
     $     VLI2(sbar/m1s) - VLI2(s2bar/m1s))
     
c     D0 sigular piece in mass reg.

      D0ms = ((t-s4)*tri2m + (s-s2)*tri3m)/
     $     ((s)*(t-m4s) - (s2)*(s4-m4s))

      D0fin = D04(s1,s2,s3,s4,s,t,M1,M2,M3,M4) - D0ms
      endif
c     we have the regularization ind. finite piece now
c     add divergent piece in the dim. reg. 
      ep = id
      tri2dreg(ep) = C1i2e(t,s4,m4s,mus,ep) !qlI3(0d0,t,s4,0.0d0,0.0d0,m4s,mus,ep) ! 2 mass ex 1 mass int
      tri3dreg(ep) = C0i2e(s,s2,mus,ep) !qlI3(0.0d0,s2,s,0.0d0,0.0d0,0.0d0,mus,ep)! 1 mass 2 ex
c     D0 singular piece in dim reg.
      D0dreg(ep) = ((t-s4)*tri2dreg(ep) + (s-s2)*tri3dreg(ep))/
     $     ((s)*(t-m4s) - (s2)*(s4-m4s))
      
      if(id.eq.0) then
         D0t2 = D0fin + D0dreg(0)
         return
      elseif(id.eq.-1) then
         D0t2 = D0dreg(-1)
         return
      else
         D0t2 = D0dreg(-2)
         return
      endif
         
      end
cccc
      complex*16 function D0t3(s1,s2,s3,s4,s,t,mi1s,
     $     mi2s,mi3s,mi4s,mus,id)

      implicit none
      integer id
      double precision s1,s2,s3,s4,s,t,mi1s,mi2s,
     $     mi3s,mi4s,delta,mus
      complex*16 ieps,sbar,tbar,s3bar,s4bar,s2bar
      parameter (ieps=(0d0,1d-16))
      double precision M1,M2,M3,M4,m1s,m2s,m3s,m4s
      double complex tri2m,tri3m,D0fin,D0ms ! mass regulated D0
      double complex tri2dreg(-2:0), tri3dreg(-2:0)
      double complex D0dreg(-2:0)
      integer ep
      double complex D04,C1i2e
      external D04
      complex*16 VLI2
      external VLI2
      parameter(delta=1.0d-6)
c
c     3.) mi4s <>0, mi3s <> 0, s1 = 0 s2 <> 0 s3 <> 0, s4 <> 0, s <>0, t<>0
c     mi4 must be nonzero mi1=mi2=0
c     mi3 must be nonzero
c      
      m1s = delta**2
      m2s = delta**2
      m3s = mi3s 
      m4s = mi4s


      if(id.eq.0) then
      M1 = sqrt(m1s)
      M2 = sqrt(m2s)
      M3 = sqrt(m3s)
      M4 = sqrt(m4s)
c      qlI3(p1sq,p2sq,p3sq,m1sq,m2sq,m3sq,musq,ep)

      tbar = t + ieps
      s3bar = s3 + ieps
      s4bar = s4 + ieps
      s2bar = s2 + ieps
               
c      tri2m = qlI3(t,s4,0d0,m1s,m4s,m2s,mus,0)

      tri2m = 1d0/(t - s4)*(log((m4s-tbar)/(delta*M4))**2 - 
     $     log((m4s-s4bar)/(delta*M4))**2 + 
     $     VLI2(tbar/m4s) - VLI2(s4bar/m4s))
      
c      tri3m = qlI3(0.0d0,s2,s,m1s,m2s,m3s,mus,0)

      sbar = s+ieps

      tri3m = 1d0/(s - s2)*(log((m3s-sbar)/(delta*M3))**2 - 
     $     log((m3s-s2bar)/(delta*M3))**2 + 
     $     VLI2(sbar/m3s) - VLI2(s2bar/m3s))
     
c     D0 sigular piece in mass reg.

      D0ms = ((t-s4)*tri2m + (s-s2)*tri3m)/
     $     ((s-m3s)*(t-m4s) - (s2-m3s)*(s4-m4s))

      D0fin = D04(s1,s2,s3,s4,s,t,M1,M2,M3,M4) - D0ms
      endif
c     we have the regularization ind. finite piece now
c     add divergent piece in the dim. reg. 
c
      ep = id
      tri2dreg(ep) = C1i2e(t,s4,m4s,mus,ep) !qlI3(0d0,t,s4,0.0d0,0.0d0,m4s,mus,ep) ! 2 mass ex 1 mass int
      tri3dreg(ep) = C1i2e(s,s2,m3s,mus,ep) !qlI3(0.0d0,s2,s,0.0d0,0.0d0,m3s,mus,ep)! 2 mass ex 1 mass int 
c     D0 singular piece in dim reg.
      D0dreg(ep) = ((t-s4)*tri2dreg(ep) + (s-s2)*tri3dreg(ep))/
     $     ((s-m3s)*(t-m4s) - (s2-m3s)*(s4-m4s))
      
      if(id.eq.0) then
         D0t3 = D0fin + D0dreg(0)
         return
      elseif(id.eq.-1) then
         D0t3 = D0dreg(-1)
         return
      else
         D0t3 = D0dreg(-2)
         return
      endif
         
      end
c---------  scalar two-point function: -------------------------------------
c  
c
c
c     short:               Binoth's notation:
c
c     B1i0e(mis,mus)       I2d(0,0;0,mis)
c     B1i1e(s,mis,mus)     I2d(s,s;0,mis)
c     B2i1e(s,mis,mis,mus) I2d(s,s;mis,mis)
c
      complex*16 function B1i0e(mis,mus,id)
c notation for scalar integrals is as:
c	
c	B0(q1,mi1s,mi2s) = (mu)^(4-d)/(i*pi^(d/2))
c
c	Int d^dq
c                  1
c      --------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s]
c
c
c
c	  	     	      A_k     B_k
c	= Gamma(1+eps) * (  ----- - ----- + C1_k ) +O(eps)
c	          	     eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: no external mass scales
c	      one internal mass	
c     A_k = 0
c     B_k = 1
c     C_k = log(mus/mis)+1
c
c	    
c	Binoth's notation: I3d(0,0;0,mis)
      implicit none
      double precision mus,mis
      integer id

      if(id.eq.0) then
         B1i0e = log(mus/mis)+1.0d0
         return
         
      elseif(id.eq.-1) then
         B1i0e = (1d0,0d0)
         return
      
      elseif(id.eq.-2) then
         B1i0e = (0d0,0d0)
         return
      
      else
         print*,'Wrong id'
         stop
      endif

      end

      complex*16 function B1i1e(s,mis,mus,id)
c notation for scalar integrals is as:
c	
c	B0(q1,mi1s,mi2s) = (mu)^(4-d)/(i*pi^(d/2))
c
c	Int d^dq
c                  1
c      --------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s]
c
c
c
c	  	     	      A_k     B_k
c	= Gamma(1+eps) * (  ----- - ----- + C1_k ) +O(eps)
c	          	     eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: two external mass scales
c	      one internal mass	
c     A_k = 0
c     B_k = 1
c     C_k = 
c
      implicit none 
      double precision mus,mis,s
      double complex C1,sbar
      integer id
      complex*16 ieps,ctwo,czero,cone
      parameter (ieps=(0d0,1d-16),ctwo=(2d0,0d0),czero=(0d0,0d0),
     $     cone=(1d0,0d0))

      sbar = s+ieps
       
      if(id.eq.0) then
         C1 = (mis-s)/s * log((mis-sbar)/mis)
         B1i1e = C1+log(mus/mis)+ctwo
         return
      
      elseif(id.eq.-1) then
         B1i1e = cone
         return

      elseif(id.eq.-2) then
         B1i1e = czero
         
         return
      
      else
         stop
         print*,'wrong id'
      endif

      end

      complex*16 function B2i1e(s,mis1,mis2,mus,id)
c notation for scalar integrals is as:
c	
c	B0(q1,mi1s,mi2s) = (mu)^(4-d)/(i*pi^(d/2))
c
c	Int d^dq
c                  1
c      --------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s]
c
c
c
c	  	     	      A_k     B_k
c	= Gamma(1+eps) * (  ----- - ----- + C1_k ) +O(eps)
c	          	     eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: two external mass scales
c	      two internal mass	
c     A_k = 0
c     B_k = 1
c     C_k = 
c
      implicit none 
      double precision mus,mis1,mis2,s
      double complex mis1bar,gam1,gam2,sbar
      integer id
      complex*16 ieps,ctwo,czero,cone
      parameter (ieps=(0d0,1d-16),ctwo=(2d0,0d0),czero=(0d0,0d0),
     $     cone=(1d0,0d0))

      sbar = s + ieps
      mis1bar = mis1-ieps
      
      
      if(id.eq.0) then
         gam1 = s-mis2+mis1 + sqrt((s-mis2+mis1)**2 - 4.0d0*s*mis1bar)
         gam1 = gam1/(2.0d0*s)
         
         gam2 = s-mis2+mis1 - sqrt((s-mis2+mis1)**2 - 4.0d0*s*mis1bar)
         gam2 = gam2/(2.0d0*s)
         
         B2i1e = ctwo + gam1*log((gam1-cone)/gam1) - log(gam1-cone) +  
     $        gam2*log((gam2-cone)/gam2) - log(gam2-cone) 
     $        - log(sbar/mus)
         return
      
      elseif(id.eq.-1) then
         B2i1e = cone
         return

      elseif(id.eq.-2) then
         B2i1e = czero
         
         return
      
      else
         stop
         print*,'wrong id'
      endif

      end

c      
c this file contains the following scalar triangle functions:
c
c	short:			Binoth's notation:
c
c	C0i1e(s,mus)	 		I3d(s,0,0;0,0,0)
c	C0i2e(s1,s2,mus) 		I3d(s1,s2,0;0,0,0)
c	C0i3e(s1,s2,s3,mus)		I3d(s1,s2,s3;0,0,0)     
c
c	C1i1e(s,mis,mus)		I3d(s,0,0;0,mis,0)
c	C1i2e(s1,s2,mis,mus)		I3d(s1,s2,0;0,mis,0)
c	C1d2e(s1,s2,mis,mus)		I3d(s1,s2,0;0,0,mis)
c	C2i2e(s1,s2,m2s,m3s,mus)	I3d(s1,s2,0;0,mis,mis)	

c	C2i3e(s1,s2,s3,m2s,m3s,mus)	I3d(s1,s2,s3;0,mis,mis) 
c      
c      IR divergent integrals are taken from Dittmaier hep-ph/0308246;
c      other integrals taken from Binoth et al., 0709.3513 [hep-ph]    
c      
c      
c---------  scalar triangle function: massless internal lines  --------------
c
      complex*16 function C0i1e(s,mus,id) 
c
c notation for scalar integrals is as:
c	
c	C0(q1,q2,mi1s,mi2s,mi3s) = (mu)^(4-d)/(i*pi^(d/2))
c
c	Int d^dq
c			       1
c      ----------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s]
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: one external mass scale s
c	      no  internal mass	
c
c		g_1e = 1/s
c
c		A_1e = 1
c		B_1e = log(mus/(-s))
c
c		C_1e = 	1/2*log(mus/(-s))**2 - pi^2/6
c
c	     
c	(see my notes master-integrals.tex)
c
c	Binoth's notation: I3d(s,0,0;0,0,0)
c       checked: TF
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision s,mus

      double complex sbar
      double complex c1e,b1e,a1e
        
      double precision zero
      parameter (zero=0d0)
      complex*16 ieps,VLI2
      parameter (ieps=(0d0,1d-16))
      external VLI2
      double precision pi,pi2o6
      parameter (pi=3.14159 26535 89793d0)
      parameter (pi2o6=pi**2/6d0)
      integer id
          
      sbar = s+ieps


      if(id.eq.0) then
         
         c1e = 0.5d0*(log((-sbar)/mus))**2-pi2o6
         
         C0i1e = c1e/s
         
         return
      elseif(id.eq.-1) then
         
         b1e = -log((-sbar)/mus)
         
         C0i1e = b1e/s
         
         return
      elseif(id.eq.-2) then
         a1e = (1d0,0d0)
         
         C0i1e = a1e/s

         return
      else
         return
      endif
      end      
      
      
c---------  scalar triangle function: massless internal lines  --------------
c
      complex*16 function C0i2e(s1,s2,mus,id) 
c
c notation for scalar integrals is as:
c	
c	C0(q1,q2,mi1s,mi2s,mi3s) = (mu)^(4-d)/(i*pi^(d/2))
c
c	Int d^dq
c			       1
c      ----------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s]
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: two external mass scales s1,s2
c	      no  internal mass	
c
c		g_2e = 1/(s1-s2)
c
c		A_2e = 0
c		B_2e = log((-s2)/(-s1))
c
c		C_2e =	 1/2*log((-s1)/mus)**2 
c			-1/2*log((-s2)/mus)**2 
c
c	     
c	(see my notes master-integrals.tex)
c
c	Binoth's notation: I3d(s1,s2,0;0,0,0)
c       Checked: TF
c ------------------------------------------------------------------------
      implicit none
      
      double precision s1,s2,mus

      double complex s1bar,s2bar
      double complex c2e,b2e,a2e
        
      double precision zero
      parameter (zero=0d0)
      complex*16 ieps,VLI2
      parameter (ieps=(0d0,1d-16))
      external VLI2
      double precision pi,pi2o6
      parameter (pi=3.14159 26535 89793d0)
      parameter (pi2o6=pi**2/6d0)
      integer id
          
      s1bar = s1+ieps
      s2bar = s2+ieps

      if(id.eq.0) then
      
         c2e = (log((-s1bar)/mus))**2-(log((-s2bar)/mus))**2
         
         C0i2e = 0.5d0*c2e/(s1-s2)

         return
      elseif(id.eq.-1) then
         b2e = log(-s2bar/mus) - log(-s1bar/mus)
         C0i2e = b2e/(s1-s2)
         return
      elseif(id.eq.-2) then
         a2e = (0d0,0d0)
         C0i2e = a2e
      else
         print*,'Wrong Id'
         stop
         
      endif
     
      end      
      
            
c---------  scalar triangle function: massless internal lines  --------------
c
      complex*16 function C0i3e(s1,s2,s3,mus,id) 
c
c notation for scalar integrals is as:
c	
c	C0(q1,q2,mi1s,mi2s,mi3s) = (mu)^(4-d)/(i*pi^(d/2))
c
c	Int d^dq
c			       1
c      ----------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s]
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: three external mass scales s1,s2,s3

c		A = 0
c		B = 0
c	     
c	taken from [Binoth et al.], see (A.7)
c
c	Binoth's notation: I3d(s1,s2,s3;0,0,0)
c       Checked: TF
c ------------------------------------------------------------------------
      implicit none
      
      double precision s1,s2,s3,mus
      double complex xp,xm,yp,ym
      double complex gki,c1
      integer id
      
      double precision zero
      parameter (zero=0d0)
      double precision lam
      complex*16 ieps,VLI2,rfunc,eta
      parameter (ieps=(0d0,1d-16))
      external VLI2,rfunc,lam,eta
      double precision pi,pi2o3
      parameter (pi=3.14159 26535 89793d0)
      parameter (pi2o3=pi**2/3.d0)
          
      if(id.eq.0) then
         xp = (s1+s3-s2-sqrt(lam(s1,s2,s3)-ieps*s1))/(2.d0*s1)
         xm = (s1+s3-s2+sqrt(lam(s1,s2,s3)-ieps*s1))/(2.d0*s1)
         
         yp = 1.-xm
         ym = 1.-xp	
         
         c1 =  2.d0*VLI2(-xm/yp)+2.d0*VLI2(-ym/xp)+pi2o3
     &        +0.5d0*(log(xm/yp))**2+0.5d0*(log(ym/xp))**2  
     &        +0.5d0*(log(xp/yp))**2-0.5d0*(log(xm/ym))**2  
         
         gki = sqrt(lam(s1,s2,s3)-ieps*s1)	
         
         C0i3e = -c1/gki
         return
      elseif(id.eq.-1 .or. id.eq.-2) then
         C0i3e = (0d0, 0d0)
         return
      else
         print*,'Wrong Id'
         stop
      endif
      
      end


c---------  scalar triangle function: massive internal lines  --------------
c
      complex*16 function C1i1e(s,mis,mus,id) 
c
c notation for scalar integrals is as:
c	
c	C0(q1,q2,mi0s,mi1s,mi2s) = (mu)^(4-d)/(i*pi^(d/2))
c
c	Int d^dq
c			       1
c      ----------------------------------------------
c	[q^2-mi0s][(q+q1)^2-mi1s][(q+q1+q2)^2-mi2s]
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: one external mass scale s attached to
c	      one internal mass	scale mis
c
c		g_1i1e = 1/s
c
c		A_1i1e = 0
c		B_1i1e = log(mis/(mis-s))
c
c		C_1i1e = log(mus/(mis))*log(mis/(mis-s))
c			+log((mis-s)/mis)**2
c			+VLI2((s)/mis)
c
c
c	     
c	(see my notes master-integrals.tex)
c
c	Binoth's notation: I3d(sq,0,0;0,mis,0)
c       Checked: TF
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision s,mis,mus

      double complex sbar
      double complex c1,b1
      integer id
      
      double precision zero
      parameter (zero=0d0)
      complex*16 ieps,VLI2
      parameter (ieps=(0d0,1d-16))
      external VLI2
      double precision pi,fpi
      parameter (pi=3.14159 26535 89793d0)
      parameter (fpi=4d0*pi)
          
      sbar = s+ieps

      if(id.eq.0) then
         c1 = log(mus/(mis))*log(mis/(mis-sbar))
     &        +log((mis-sbar)/mis)**2
     &        +VLI2(sbar/mis)
         
         C1i1e = c1/s
         return
      elseif(id.eq.-1) then
         b1 = log(mis/mus) - log((mis-sbar)/mus)
         C1i1e = b1/s
         return
      elseif(id.eq.-2) then
         
         C1i1e = (0d0, 0d0)
         return
      else 
         print*,'Wrong Id'
         stop
      endif

      end


c---------  scalar triangle function: massive internal lines  --------------
c
      complex*16 function C1i2e(s1,s2,mis,mus,id) 
c
c notation for scalar integrals is as:
c	
c	C0(q1,q2,mi1s,mi2s,mi3s) = (mu)^(4-d)/(i*pi^(d/2))
c
c	Int d^dq
c			       1
c      ----------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s]
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: two external mass scales s1,s2
c	      one internal mass	scale mis in between
c
c		g_1i2e = 1/(s1-s2)
c
c		A_1i2e = 0
c		B_1i2e = log((mis-s2)/(mis-s1))
c
c		C_1i2e = log(mus/(mis))*log((mis-s2)/(mis-s1))
c			+log((mis-s1)/mis)**2
c			-log((mis-s2)/mis)**2
c			+VLI2((s1)/mis)-VLI2((s2)/mis)
c
c
c	     
c	(see my notes master-integrals.tex)
c
c	Binoth's notation: I3d(s1,s2,0;0,mis,0)
c       Checked: TF
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision s1,s2,mis,mus

      double complex s1bar,s2bar
      double complex c1,a1,b1
      
      double precision zero
      parameter (zero=0d0)
      complex*16 ieps,VLI2
      parameter (ieps=(0d0,1d-16))
      external VLI2
      double precision pi,fpi
      parameter (pi=3.14159 26535 89793d0)
      parameter (fpi=4d0*pi)
      integer id
          
      s1bar = s1+ieps
      s2bar = s2+ieps

      if(id.eq.0) then

      	c1 = log(mus/(mis))*log((mis-s2bar)/(mis-s1bar))
     &	    +(log((mis-s1bar)/mis))**2-(log((mis-s2bar)/mis))**2
     &	    +VLI2(s1bar/mis)-VLI2(s2bar/mis)
     
	C1i2e = c1/(s1-s2)

        return

      elseif(id.eq.-1) then
         b1 = log((mis-s2bar)/mus) - log((mis-s1bar)/mus)
         C1i2e = b1/(s1-s2)

         return
      elseif(id.eq.-2) then
         a1 = (0d0, 0d0)
         C1i2e = a1
         return

      else
         print*,'Wrong id!'
         stop
      endif
		            
      end


c---------  scalar triangle function: massive internal lines  --------------
c
      complex*16 function C1d2e(s1,s2,mis,mus,id) 
c
c notation for scalar integrals is as:
c	
c	C0(q1,q2,mi1s,mi2s,mi3s) = (mu)^(4-d)/(i*pi^(d/2))
c
c	Int d^dq
c			       1
c      ----------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s]
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: two external mass scales s1,s2
c	      one internal mass	scale mis on opposite side
c
c	attention: not symmetric in s1,s2!!
c		
c		g_1d2e = 1/(s1-s2)
c
c		A_1d2e = 0
c		B_1d2e = 0
c
c		C_1d2e =
c
c
c	     
c	taken from [Binoth et al.]
c
c	Binoth's notation: I3d(s1,s2,0;0,0,mis)
c       Checked: TF
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision s1,s2,mis,mus
      double precision x0
      double complex x1bar,s1bar
      double complex c1
      integer id
      
      double precision zero
      parameter (zero=0d0)
      complex*16 ieps,VLI2,rfunc
      parameter (ieps=(0d0,1d-16))
      external VLI2,rfunc
      double precision pi,pi2o6
      parameter (pi=3.14159 26535 89793d0)
      parameter (pi2o6=pi**2/6.d0)
        

      if (id.eq.0) then
         s1bar = s1+ieps
         
         x0 = s1/(s1-s2)
         x1bar = s1bar/(s1-s2+mis)
         
         c1 = rfunc(dcmplx(x0),x1bar)-pi2o6+VLI2(1.-1/(x0-ieps))
         
         C1d2e = c1/(s2-s1)
		            
         return
      elseif(id.eq.-1 .or. id.eq.-2) then

         C1d2e = (0d0,0d0)
         return
      else
         print*,'Wrong id'
         stop
      endif
      end      
         
           
c---------  scalar triangle function: massive internal lines  --------------
c
      complex*16 function C2i2e(s1,s2,m2s,m3s,mus,id) 
c
c notation for scalar integrals is as:
c	
c	C0(q1,q2,mi1s,mi2s,mi3s) = (mu)^(4-d)/(i*pi^(d/2))
c
c	Int d^dq
c			       1
c      ----------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s]
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: two external mass scales s1,s2
c	      two equal internal mass scales m3s=m2s 

c		A = 0
c		B = 0
c	
c	(does not depend on mus up to order eps^0)
c
c	     
c	taken from [Binoth et al.], (A.11)
c
c	Binoth's notation: I3d(s1,s2,0;0,mis,mis)
c       Checked: TF
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision s1,s2,m2s,m3s,mis,mus
      double complex s1bar,s2bar
      double complex xp,xm,eta0
      double complex c1,clog
      double precision gki,x0
      integer id

      double precision zero
      parameter (zero=0d0)
      double complex czero
      parameter (czero=(0d0,0d0))
      complex*16 ieps,VLI2,rfunc,eta
      parameter (ieps=(0d0,1d-16))
      external VLI2,rfunc,eta
          
      if (m3s.ne.m2s) then
      	print*,'C2i2e is implemented for two equal internal mases only!'
	print*,'set C2i2e = 0'
	C2i2e = 0d0
	return
      endif	
      

      if(id.eq.0) then
         mis = m2s
         
         s2bar = s2+ieps
         s1bar = s1+ieps
         
         x0 = 1.d0-s1/s2
         
         xp = 0.5d0*(1+sqrt(1d0-4*mis/s2bar)) 	
         xm = 0.5d0*(1-sqrt(1d0-4*mis/s2bar)) 	
         
         eta0 = eta(1d0-(s1bar/mis)*x0,mis/(mis-s2bar*x0*(1d0-x0))) 
         
         if (eta0.eq.czero) then
            clog = czero
         else
            clog = -eta0*log((1d0-x0)/(-x0))
         endif		
         
         
         c1 = VLI2(s1bar/mis)-VLI2(1d0/xp)-VLI2(1d0/xm)
     &        +rfunc(dcmplx(x0),mis/s1bar)+rfunc(dcmplx(1d0-x0),xm)
     &        -rfunc(dcmplx(x0),xm)
         
         gki = s2-s1
         
         C2i2e = (c1+clog)/gki
         return
      elseif(id.eq.-1 .or. id.eq.-2) then
         C2i2e = (0d0,0d0)
         return
      endif
      end      
     
           
c---------  scalar triangle function: massive internal lines  --------------
c
      complex*16 function C2i3e(s1,s2,s3,m2s,m3s,mus,id) 
c
c notation for scalar integrals is as:
c	
c	C0(q1,q2,mi1s,mi2s,mi3s) = (2*pi*mu)^(4-d)/(i*pi^2)
c
c	Int d^dq
c			       1
c      ----------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s]
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: three external mass scales s1,s2,s3
c	      two equal internal mass scales m2s=m3s 

c		A = 0
c		B = 0
c	     
c	taken from [Binoth et al.]
c
c	Binoth's notation: I3d(s1,s2,s3;0,mis,mis)
c Check this one!
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision s1,s2,s3,m3s,m2s,mus
      double complex x0bar,x1bar,s1bar,s2bar
      double complex xp,xm,etap,etam
      double complex gki,c1
      integer id
      
      double precision zero
      parameter (zero=0d0)
      double precision lam
      complex*16 ieps,VLI2,rfunc,eta
      parameter (ieps=(0d0,1d-16))
      external VLI2,rfunc,lam,eta
      double precision pi,fpi,pi2o6
      parameter (pi=3.14159 26535 89793d0)
      parameter (fpi=4d0*pi,pi2o6=pi**2/6.d0)
          
      if (m3s.ne.m2s) then
      	print*,'C2i3e is implemented for two equal internal mases only!'
	print*,'set C2i3e = 0'
	C2i3e = 0d0
	return
      endif

      if(id.eq.0) then
	
         s2bar = s2+ieps
         s1bar = s1+ieps
         
         xp = (s1+s2-s3-sqrt(lam(s1,s2,s3)-ieps*s2))/(2.d0*s2)	
         xm = (s1+s2-s3+sqrt(lam(s1,s2,s3)-ieps*s2))/(2.d0*s2)	
         
         x0bar = 0.5d0*(1.d0-sqrt(1.-4*m3s/s2bar))
         x1bar = (s1bar-m3s)/(s1-s3)
         
         etap = eta(-s3*xp-s1*(1.-xp)+m3s-ieps,1/(m3s-xp*(1.-xp)*s2bar))
         etam = eta(-s3*xm-s1*(1.-xm)+m3s-ieps,1/(m3s-xm*(1.-xm)*s2bar))
         
         
         c1 = rfunc(xm,x1bar)-rfunc(xp,x1bar)+rfunc(1.-xm,x0bar)
     &        -rfunc(1.-xp,x0bar)-rfunc(xm,x0bar)+rfunc(xp,x0bar)
     &        -etam*log((1.-xm)/(-xm))+etap*log((1.-xp)/(-xp))      
         
         gki = sqrt(lam(s1,s2,s3)-ieps*s2)	
         
         C2i3e = c1/gki
         
      return
      elseif(id.eq.-1 .or. id.eq.-2) then
         C2i3e = (0d0,0d0)
         return
      else
         print*,'Wrong Id'
         stop
      endif
      end  
c---------  scalar triangle function: massive internal lines  --------------
c
      complex*16 function C3i3e(s1,s2,s3,m1s,m2s,m3s,mus,id) 
c
c notation for scalar integrals is as:
c	
c	C0(q1,q2,mi1s,mi2s,mi3s) = (2*pi*mu)^(4-d)/(i*pi^2)
c
c	Int d^dq
c			       1
c      ----------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s]
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: three external mass scales s1,s2,s3
c	      two equal internal mass scales m2s=m3s 

c		A = 0
c		B = 0
c	     
c       General result from A.Denner, Fortschr. Phys. 41 (1993) 307]
c
c	Binoth's notation: I3d(s1,s2,s3;mis1,mis2,mis3)
c
c ------------------------------------------------------------------------  
      implicit none
      double precision s1,s2,s3,m1s,m2s,m3s,mus
      double precision m02,m12,m22,p01,p12,p20
      double precision lam
      external lam
      complex*16 c0_,ieps,alpha,alp(0:2)
      complex*16 y0(0:2),y(0:2,-1:1),x(0:2,-1:1)
      real*8 thp(0:2),thy(0:2)
      integer i,j,id
      double precision eps,pi
      parameter (pi=3.14159 26535 89793d0,eps=1.0d-16)
      
      complex*16 VLI2,eta
      external VLI2,eta
c      
c     T
c
      if((s1*s2.eq.0d0).or. ((m1s+m2s+m3s).eq.0d0)) then
         print*,'Wrong kinematics in C3i3e!'
         stop
      endif

      if(id.eq.0) then
      
c     pi   = 4d0*datan(1d0)                                               
c     eps  = 1d-17
         ieps = dcmplx(0d0,eps)
c     
         m02 = m1s
         m12 = m2s
         m22 = m3s
         p01 = s1
         p12 = s2
         p20 = s3
         
C***  Regular C0 function with p01,p12,p20 =/= 0
         alpha  = sqrt( lam(p01,p12,p20) )
         alp(0) = sqrt( lam(p12,m12,m22)*(1d0+ieps*sign(1d0,p12)) )
         alp(1) = sqrt( lam(p20,m22,m02)*(1d0+ieps*sign(1d0,p20)) )
         alp(2) = sqrt( lam(p01,m02,m12)*(1d0+ieps*sign(1d0,p01)) )
         
         do i=0,2
            if (alp(i).eq.dcmplx(0d0,0d0)) alp(i) = ieps*abs(alpha)
         enddo
         
         y0(0)  = ( p12*(p12-p01-p20+2d0*m02-m12-m22)
     &        - (p20-p01)*(m12-m22)+alpha*(p12-m12+m22) )/2d0/alpha/p12
         if(abs(p20).gt.1d-6) then
           y0(1)  = ( p20*(p20-p12-p01+2d0*m12-m22-m02)
     &          - (p01-p12)*(m22-m02)+alpha*(p20-m22+m02) )/2d0/alpha/p20
         else
           y0(1)=0d0
         endif
         y0(2)  = ( p01*(p01-p20-p12+2d0*m22-m02-m12)
     &        - (p12-p20)*(m02-m12)+alpha*(p01-m02+m12) )/2d0/alpha/p01
         
         do j=-1,1,2
            x(0,j) = (p12-m12+m22+j*alp(0))/2d0/p12
            if(abs(p20).gt.1d-6) then
              x(1,j) = (p20-m22+m02+j*alp(1))/2d0/p20
            else
              x(1,j) = 0d0
            endif
            x(2,j) = (p01-m02+m12+j*alp(2))/2d0/p01
            do i=0,2
               y(i,j) = y0(i)-x(i,j)
            enddo
         enddo
            
         do i=0,2
            thp(i) = 0d0
            thy(i) = 0d0
            if (dimag(y(i,+1)*y(i,-1)).le.0d0) thy(i) = 1d0
         enddo
         if (p12.le.0d0) thp(0) = 1d0
         if (p20.le.0d0) thp(1) = 1d0
         if (p01.le.0d0) thp(2) = 1d0
         
         c0_ = 0d0
         do i=0,2
            if(abs(p20).gt.1d-6.or.i.ne.1) then
            do j=-1,1,2
               c0_ = c0_ + VLI2((y0(i)-1d0)/y(i,j)) - VLI2(y0(i)/y(i,j))
     &              + eta(1d0-x(i,j),1d0/y(i,j))*log((y0(i)-1d0)/y(i,j))
     &              - eta(   -x(i,j),1d0/y(i,j))*log(y0(i)/y(i,j))
            enddo
             c0_ = c0_ - log((y0(i)-1d0)/y0(i))*(
     &            eta(-x(i,+1),-x(i,-1))-eta(y(i,+1),y(i,-1))
     &            - dcmplx(0d0,2d0*pi)*thp(i)*thy(i) )
             endif
          enddo
          
          c0_ = c0_/alpha
          C3i3e = c0_
          return
          
       elseif(id.eq.-1 .or. id.eq.-2) then
          C3i3e = (0d0,0d0)
          return
       else
          print*,'Wrong id'
          stop
       endif
       
       end
c
c-----------------------------------------------------         
c---------  auxiliary functions for computation of loop integrals  --------
c
      complex*16 function rfunc(y0,z) 
c	     
c	taken from [Binoth et al.], see Eq.(A.4)
c
c ------------------------------------------------------------------------
      implicit none

      double complex y0,z
      
      double complex z1,z2,eta1,eta2
      
      double precision zero
      parameter (zero=0d0)
      complex*16 VLI2,eta
      external VLI2,eta     
 
      z1 = y0/(y0-z)
      z2 = (y0-1.d0)/(y0-z)
            
      eta1 = eta(dcmplx(-z),dcmplx(1/(y0-z)))
      eta2 = eta(dcmplx(1.-z),dcmplx(1/(y0-z)))
      
      rfunc = VLI2(dcmplx(z1))-VLI2(dcmplx(z2))+eta1*log(z1)-eta2*log(z2)


      return
      end      
         
c-------------------------------------------------------------------
c
      double precision function lam(x,y,z) 
c	     
c	taken from [Binoth et al.], see Eq.(A.2)
c
c -------------------------------------------------------------------
      implicit none
      
      double precision x,y,z
 
      lam = x**2+y**2+z**2-2*x*y-2*y*z-2*x*z
		            
      return
      end  
                

cc-------------------------------------------------------------------------
c***********************************************************************
c        FUNCTION ETA(C1,C2)                                            
c***********************************************************************
c*       COMPLEX ETA-FUNKTION                                           
c*---------------------------------------------------------------------*
c*       8.06.90    ANSGAR DENNER                                       
c***********************************************************************
cc        IMPLICIT   LOGICAL(A-Z)
c	
c	implicit none
c	                                        
c        COMPLEX*16 ETA,C1,C2                                           
c        REAL*8     PI,IM1,IM2,IM12                                     
c                                                                       
c        PI     = 4D0*DATAN(1D0)                                        
c        IM1    = DIMAG(C1)                                             
c        IM2    = DIMAG(C2)                                             
c        IM12   = DIMAG(C1*C2)                                          
c 
c	if (((IM1.eq.0d0).and.(DREAL(C1).lt.0d0)).or.
c     &	    ((IM2.eq.0d0).and.(DREAL(C2).lt.0d0)).or.
c     &	    ((IM12.eq.0d0).and.(DREAL(C1*C2).lt.0d0))) then
c	  write(*,*) 'eta function on cut !!!'
c	  write(*,*) 'C1    = ',C1
c	  write(*,*) 'C2    = ',C2
c	  write(*,*) 'C1*C2 = ',C1*C2
c	  stop
c	endif
c                                                                      
c        IF(IM1.LT.0D0.AND.IM2.LT.0D0.AND.IM12.GT.0D0) THEN             
c            ETA = DCMPLX(0D0,2D0*PI)                                   
c        ELSE IF (IM1.GT.0D0.AND.IM2.GT.0D0.AND.IM12.LT.0D0) THEN       
c            ETA = DCMPLX(0D0,-2D0*PI)                                  
c        ELSE                                                           
c            ETA = DCMPLX(0D0)                                          
c        END IF                                                         
c        END                                                            
      
c% FC %% 
c% FC %% -----------  spence function  VLI2(z) ------------------------
c% FC %% 
c% FC %%    complex*16 function li2(zin)
c% FC %%    implicit none
c% FC %%    complex*16 zin, z, u, u2, unpo, ans, zext
c% FC %%    double precision r, r2, r2n, fac
c% FC %% 
c% FC %% etermine the value of the dilogarithm 
c% FC %% 
c% FC %%   li2(z) = - int_0^1  log(1-zt)/t dt  with cut along the positive 
c% FC %%                                       real axis, z>1
c% FC %% 
c% FC %%      Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c% FC %%      Initial version:  2000 November 6
c% FC %%      Last modified:    2000 November 12
c% FC %% 
c% FC %%    integer i
c% FC %%    double precision c0,c1,c2,c4,c6,c8,c10,c12,c14,c16,c18,c20,c22
c% FC %%    double precision b0,b1,b2,b4,b6,b8,b10,b12,b14,b16,b18,b20,b22
c% FC %%    double precision d0,d1,d2,d4,d6,d8,d10,d12,d14,d16,d18,d20,d22
c% FC %%    parameter (b0=1d0,            d0 =1d0,      c0= b0/d0)
c% FC %%    parameter (b1=-1d0/2d0,       d1 =d0*2d0,   c1= b1/d1)
c% FC %%    parameter (b2= 1d0/6d0,       d2 =d1*3d0,   c2= b2/d2)
c% FC %%    parameter (b4=-1d0/30d0,      d4 =d2*20d0,  c4= b4/d4)
c% FC %%    parameter (b6=1d0/42d0,       d6 =d4*42d0,  c6= b6/d6)
c% FC %%    parameter (b8=-1d0/30d0,      d8 =d6*72d0,  c8= b8/d8)
c% FC %%    parameter (b10=5d0/66d0,      d10=d8*110d0, c10=b10/d10)
c% FC %%    parameter (b12=-691d0/2730d0, d12=d10*156d0,c12=b12/d12)
c% FC %%    parameter (b14=7d0/6d0,       d14=d12*210d0,c14=b14/d14)
c% FC %%    parameter (b16=-3617d0/510d0, d16=d14*272d0,c16=b16/d16)
c% FC %%    parameter (b18=43867d0/798d0, d18=d16*342d0,c18=b18/d18)
c% FC %%    parameter (b20=-174611d0/330d0,d20=d18*420d0,c20=b20/d20)
c% FC %%    parameter (b22=854513d0/138d0,d22=d20*506d0,c22=b22/d22)
c% FC %%    double precision eps, epst, pi, pi2o6
c% FC %%    parameter (eps=1d-16, epst=1d-3)
c% FC %%    parameter (pi=3.14159 26535 89793238d0, pi2o6=pi**2/6d0)
c% FC %% 
c% FC %% ebug information
c% FC %%    logical ldebug
c% FC %%    parameter (ldebug=.false.)
c% FC %% 
c% FC %%    z = zin
c% FC %%     print*," li2 call with z = ",z
c% FC %%    u = z**2
c% FC %%    r2 = dreal(z)**2+dimag(z)**2 
c% FC %%    if (r2.lt.eps) then
c% FC %%       li2 = z + u/4d0
c% FC %%       return
c% FC %%    elseif (r2.lt.epst) then
c% FC %%       ans = z + u/4d0
c% FC %%       do i = 3,11
c% FC %%          u = u*z
c% FC %%          ans = ans + u/i**2
c% FC %%       enddo
c% FC %%       li2 = ans
c% FC %%       return
c% FC %%    endif
c% FC %%    if (dreal(z).ge.1d0 .and. dimag(z).eq.0 ) then
c% FC %%       z = z + (0d0,1d0)*eps
c% FC %%    endif
c% FC %% 
c% FC %% se z-->1/z and z--> 1-z mappings of the spence function to restrict 
c% FC %% gument to unit circle in the complex plane with Re(z) <= 0.5
c% FC %% 
c% FC %%    zext = (0d0,0d0)
c% FC %%    fac = 1
c% FC %%    if (r2.gt.1d0) then     ! map z ---> 1/z
c% FC %%       fac = -fac
c% FC %%       zext = -pi2o6 - 0.5d0*(log(-z))**2
c% FC %%       z = 1d0/z
c% FC %%    endif
c% FC %%    if (dreal(z).gt.0.5d0) then     ! map new z ---> 1-z
c% FC %%       zext = zext + fac*(pi2o6-log(z)*log(1-z))
c% FC %%       fac = -fac
c% FC %%       z = 1-z
c% FC %%    endif
c% FC %% 
c% FC %% ow use t = 1 - exp(-u) mapping to write Li(z) in terms of Bernoulli 
c% FC %% umbers
c% FC %% 
c% FC %%    u = - log(1-z)
c% FC %%    r2 = abs(u)**2
c% FC %%    u2 = u*u
c% FC %%    ans = u*(c0 + u*(c1+c2*u))
c% FC %%    r2n = r2*r2       !r^4
c% FC %% 
c% FC %%    unpo = u2*u2*u
c% FC %%    ans = ans + c4*unpo
c% FC %% 
c% FC %%    unpo = unpo*u2
c% FC %%    ans = ans + c6*unpo
c% FC %% 
c% FC %%    r = r2n*r2n       !r^8
c% FC %%    unpo = unpo*u2
c% FC %%    ans = ans + c8*unpo
c% FC %% 
c% FC %%    r2n = r*r2        !r^10 
c% FC %%    if ((r2n*c10).gt.eps) then
c% FC %%       unpo = unpo*u2
c% FC %%       ans = ans + c10*unpo
c% FC %%    else
c% FC %%       li2 = fac * ans + zext
c% FC %%       if (ldebug) print*," exit li2s at n=8 "
c% FC %%       return
c% FC %%    endif
c% FC %% 
c% FC %%    unpo = unpo*u2
c% FC %%    ans = ans + c12*unpo
c% FC %% 
c% FC %%    unpo = unpo*u2
c% FC %%    ans = ans + c14*unpo
c% FC %% 
c% FC %%    unpo = unpo*u2
c% FC %%    ans = ans + c16*unpo
c% FC %% 
c% FC %%    r2n = r2n*r
c% FC %%    if ((r2n*c18).gt.eps) then
c% FC %%       unpo = unpo*u2
c% FC %%       ans = ans + c18*unpo
c% FC %%    else
c% FC %%       li2 = fac * ans + zext
c% FC %%       if (ldebug) print*," exit li2s at n=16 "
c% FC %%       return
c% FC %%    endif
c% FC %% 
c% FC %%    unpo = unpo*u2
c% FC %%    ans = ans + c20*unpo
c% FC %% 
c% FC %%    unpo = unpo*u2
c% FC %%    ans = ans + c22*unpo
c% FC %% 
c% FC %%    li2 = fac * ans + zext
c% FC %%    end
c% FC %% ************************************************************************
c% FC %% c        FUNCTION D04(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4)
c% FC %%         FUNCTION D04(P1t,P2t,P3t,P4t,P12t,P23t,M1t,M2t,M3t,M4t)
c% FC %% ************************************************************************
c% FC %% *  SCALAR 4-POINT FUNCTION WITH AT LEAST ONE MASS ZERO                 *
c% FC %% *  P1,P2,P3,P4 = SQUARED EXTERNAL MOMENTA			       *
c% FC %% *  P12 = (p1+p2)**2,  P23 = (p2+p3)**2				       *
c% FC %% *----------------------------------------------------------------------*
c% FC %% *  2.1.92  SD	         					       *
c% FC %% *  rearrangement to exploit massless external momenta   14.3.01  DZ    *
c% FC %% *  Modified: Michael Kubocz                                            *
c% FC %% *  Interception of NANs e.g. caused by log(0) etc. (see below)         *
c% FC %% ************************************************************************
c% FC %%         IMPLICIT REAL*8 (A-Z)
c% FC %% 	REAL*8 M(4),P(4,4),K(4,4)
c% FC %%         real*8 pi,eps,eps1
c% FC %%         real*8 im1,im2
c% FC %%         real*8 m1,m2,m3,m4
c% FC %%         real*8 m1t,m2t,m3t,m4t
c% FC %%         real*8 m02,m12,m22,m32,m42
c% FC %%         real*8 mm0,mm1,mm2,mm3,mm4 
c% FC %%         real*8 p1,p2,p3,p4,p12,p23
c% FC %%         real*8 p1t,p2t,p3t,p4t,p12t,p23t
c% FC %%         real*8 q0,q1,q2,q3,q4,q00,q12,q23
c% FC %% 	COMPLEX*16 A1,A2,A3,A4,SWAP
c% FC %% 	COMPLEX*16 SS(4), XX(2), X(2,4),RS(4,4)
c% FC %% 	COMPLEX*16 S0(4),XX0(2),X0(2,4), R(4,4),G(2)
c% FC %%         COMPLEX*16 D04,LI2,ETA,SQE,ETAS
c% FC %% 	COMPLEX*16 AA,BB,CC,DD,IEPS,H,HH,L1,L2,L3,L4
c% FC %% 	INTEGER I,J, i1,i2,i3,i4
c% FC %%         complex*16 D1,D2,D3,D4
c% FC %% 
c% FC %%         EXTERNAL LI2
c% FC %% 
c% FC %%         D1=DCMPLX(0d0,0d0)
c% FC %%         D2=DCMPLX(0d0,0d0)
c% FC %%         D3=DCMPLX(0d0,0d0)
c% FC %%         D4=DCMPLX(0d0,0d0)
c% FC %%         eps1=1d-7
c% FC %% 
c% FC %%         if(abs(P1t).le.eps1) then
c% FC %%            P1=0d0
c% FC %%         else
c% FC %%            P1=P1t
c% FC %%         endif
c% FC %%         if(abs(P2t).le.eps1) then
c% FC %%            P2=0d0
c% FC %%         else
c% FC %%            P2=P2t
c% FC %%         endif
c% FC %%         if(abs(P3t).le.eps1) then
c% FC %%            P3=0d0
c% FC %%         else
c% FC %%            P3=P3t
c% FC %%         endif
c% FC %%         if(abs(P4t).le.eps1) then
c% FC %%            P4=0d0
c% FC %%         else
c% FC %%            P4=P4t
c% FC %%         endif
c% FC %%         if(abs(P12t).le.eps1) then
c% FC %%            P12=0d0
c% FC %%         else
c% FC %%            P12=P12t
c% FC %%         endif
c% FC %%         if(abs(P23t).le.eps1) then
c% FC %%            P23=0d0
c% FC %%         else
c% FC %%            P23=P23t
c% FC %%         endif
c% FC %%         if(abs(M1t).le.eps1) then
c% FC %%            M1=0d0
c% FC %%         else
c% FC %%            M1=M1t
c% FC %%         endif
c% FC %%         if(abs(M2t).le.eps1) then
c% FC %%            M2=0d0
c% FC %%         else
c% FC %%            M2=M2t
c% FC %%         endif
c% FC %%         if(abs(M3t).le.eps1) then
c% FC %%            M3=0d0
c% FC %%         else
c% FC %%            M3=M3t
c% FC %%         endif
c% FC %%         if(abs(M4t).le.eps1) then
c% FC %%            M4=0d0
c% FC %%         else
c% FC %%            M4=M4t
c% FC %%         endif
c% FC %% 
c% FC %%         MM1=M1
c% FC %%         MM2=M2
c% FC %%         MM3=M3
c% FC %%         MM4=M4
c% FC %%         M12=M1*M1
c% FC %%         M22=M2*M2
c% FC %%         M32=M3*M3
c% FC %%         M42=M4*M4
c% FC %%         Q1=P1
c% FC %%         Q2=P2
c% FC %%         Q3=P3
c% FC %% 	Q4=P4
c% FC %%         Q12=P12
c% FC %%         Q23=P23
c% FC %% 
c% FC %% C	IS AT LEAST ONE MASS ZERO ???
c% FC %% 	IF (MM1*MM2*MM3*MM4.NE.0D0) GOTO 130
c% FC %% 
c% FC %% C	PERMUTATE UNTIL MM3=0D0
c% FC %% 	GOTO 20
c% FC %% 10	CONTINUE
c% FC %% 	MM0=MM1
c% FC %% 	MM1=MM2
c% FC %% 	MM2=MM3
c% FC %% 	MM3=MM4
c% FC %% 	MM4=MM0
c% FC %% 	M02=M12
c% FC %% 	M12=M22
c% FC %% 	M22=M32
c% FC %% 	M32=M42
c% FC %% 	M42=M02
c% FC %% 	Q00=Q12
c% FC %% 	Q12=Q23
c% FC %% 	Q23=Q00
c% FC %% 	Q0=Q1
c% FC %% 	Q1=Q2
c% FC %% 	Q2=Q3
c% FC %% 	Q3=Q4
c% FC %% 	Q4=Q0
c% FC %% 20	IF (MM3.NE.0D0) GOTO 10
c% FC %% C	ONLY MM3 IS ZERO
c% FC %% 	IF (MM1*MM2*MM4.NE.0D0) GOTO 30
c% FC %% C	ONLY MM3 AND MM4 ARE ZERO ==> 3->2, 4->3...
c% FC %% 	IF ((MM1*MM2.NE.0D0).AND.(MM4.EQ.0D0)) GOTO 10
c% FC %% C	ONLY MM2 AND MM3 ARE ZERO
c% FC %% 	IF ((MM1*MM4.NE.0D0).AND.(MM2.EQ.0D0)) GOTO 40
c% FC %% 	WRITE(*,*)'CASE OF THIS SPECIAL D0-FUNCTION NOT IMPLEMENTED!'
c% FC %% 	STOP
c% FC %% 
c% FC %% C	****** NO MASS EQUAL TO ZERO ******
c% FC %% 130	CONTINUE
c% FC %% 	EPS=1D-18
c% FC %% 	IEPS=DCMPLX(0D0,EPS)
c% FC %% c check for massless external momentum: excellent candidate for p13,
c% FC %% c leading to r13 >=1 and real.
c% FC %% c$$$        if (q1*q2*q3*q4.eq.0d0) then !org
c% FC %%         if (abs(q1*q2*q3*q4).le.eps1) then
c% FC %% c$$$           if (q2.eq.0d0) then
c% FC %%            if (abs(q2).le.eps1) then
c% FC %%               I1 = 2
c% FC %%               I2 = 3
c% FC %%               I3 = 1
c% FC %%               I4 = 4
c% FC %% c$$$           elseif (q1.eq.0d0) then
c% FC %%            elseif (abs(q1).le.eps1) then
c% FC %%               I1 = 1
c% FC %%               I2 = 3
c% FC %%               I3 = 2
c% FC %%               I4 = 4
c% FC %% c$$$           elseif (q3.eq.0d0) then
c% FC %%            elseif (abs(q3).le.eps1) then
c% FC %%               I1 = 2
c% FC %%               I2 = 4
c% FC %%               I3 = 1
c% FC %%               I4 = 3
c% FC %%            else
c% FC %%               I1 = 1
c% FC %%               I2 = 4
c% FC %%               I3 = 2
c% FC %%               I4 = 3
c% FC %%            endif
c% FC %%            M(i1)=MM1
c% FC %%            M(i2)=MM2
c% FC %%            M(i3)=MM3
c% FC %%            M(i4)=MM4
c% FC %%            P(i1,i2)=Q1
c% FC %%            P(i3,i2)=Q2
c% FC %%            P(i3,i4)=Q3
c% FC %%            P(i1,i4)=Q4
c% FC %%            P(i1,i3)=Q12
c% FC %%            P(i3,i1)=Q12
c% FC %%            P(i2,i4)=Q23
c% FC %%            P(i4,i2)=Q23
c% FC %% 	ELSEIF( ABS((MM1**2+MM3**2-Q12)/MM1/MM3).LT.2D0 ) THEN
c% FC %% C	R13 WOULD BE NOT REAL. -> PERMUTATION! -> R(2,4) IS NOT REAL.
c% FC %% 	   M(1)=MM2
c% FC %% 	   M(2)=MM3
c% FC %% 	   M(3)=MM4
c% FC %% 	   M(4)=MM1
c% FC %% 	   P(1,2)=Q2
c% FC %% 	   P(1,3)=Q23
c% FC %% 	   P(1,4)=Q1
c% FC %% 	   P(2,3)=Q3
c% FC %% 	   P(2,4)=Q12
c% FC %% 	   P(3,4)=Q4
c% FC %% 	ELSE
c% FC %% C	R(1,3) IS REAL.
c% FC %% 	   M(1)=MM1
c% FC %% 	   M(2)=MM2
c% FC %% 	   M(3)=MM3
c% FC %% 	   M(4)=MM4
c% FC %% 	   P(1,2)=Q1
c% FC %% 	   P(1,3)=Q12
c% FC %% 	   P(1,4)=Q4
c% FC %% 	   P(2,3)=Q2
c% FC %% 	   P(2,4)=Q23
c% FC %% 	   P(3,4)=Q3
c% FC %% 	ENDIF
c% FC %% 
c% FC %% 	DO 11 J=2,4
c% FC %% 	DO 11 I=1,J-1
c% FC %% 	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
c% FC %% 	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
c% FC %%      *	            DCMPLX(1D0,0D0))
c% FC %%         IF( DBLE(K(I,J)).LT.-2D0 ) THEN
c% FC %% c        IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
c% FC %% 	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
c% FC %%      *	               DCMPLX(1D0,0D0))
c% FC %% 	ELSE
c% FC %% 	   RS(I,J)=R(I,J)
c% FC %% 	ENDIF
c% FC %% 11	CONTINUE
c% FC %% 
c% FC %% 	SS(1)=RS(1,2)
c% FC %% 	SS(2)=RS(2,3)
c% FC %% 	SS(3)=RS(3,4)
c% FC %% 	SS(4)=RS(1,4)
c% FC %% 	S0(1)=R(1,2)
c% FC %% 	S0(2)=R(2,3)
c% FC %% 	S0(3)=R(3,4)
c% FC %% 	S0(4)=R(1,4)
c% FC %% 	AA=K(3,4)/R(2,4)+R(1,3)*K(1,2)-K(1,4)*R(1,3)/R(2,4)-K(2,3)
c% FC %% 	BB=(R(2,4)-1D0/R(2,4))*(R(1,3)-1D0/R(1,3))
c% FC %%      *		+K(1,2)*K(3,4)-K(1,4)*K(2,3)
c% FC %% 	CC=K(1,2)/R(1,3)+R(2,4)*K(3,4)-K(1,4)*R(2,4)/R(1,3)-K(2,3)
c% FC %% 	DD=K(2,3)-R(1,3)*K(1,2)-R(2,4)*K(3,4)+R(1,3)*R(2,4)*K(1,4)
c% FC %% 	XX(1)=SQE(AA,BB,CC+IEPS*DD)
c% FC %% 	XX(2)=(CC+IEPS*DD)/AA/XX(1)
c% FC %% 	XX0(1)=SQE(AA,BB,CC)
c% FC %% 	XX0(2)=CC/AA/XX0(1)
c% FC %% c	IF (ABS(DREAL(XX0(1)-XX(2))).LT.ABS(DREAL(XX0(1)-XX(1)))) THEN
c% FC %% 	IF (ABS(XX0(1)-XX(2)).LT.ABS(XX0(1)-XX(1))) THEN
c% FC %% 	  SWAP  =XX0(1)
c% FC %% 	  XX0(1)=XX0(2)
c% FC %% 	  XX0(2)=SWAP
c% FC %% 	ENDIF
c% FC %% 
c% FC %% 	DO 12 I=1,2
c% FC %% 	G(I)  =SIGN( 1D0,DREAL(AA*(XX(I)-XX(3-I))) )
c% FC %% 	 X(I,1)= XX(I)/R(2,4)
c% FC %% 	X0(I,1)=XX0(I)/R(2,4)
c% FC %% 	 X(I,2)= XX(I)/R(2,4)*R(1,3)
c% FC %% 	X0(I,2)=XX0(I)/R(2,4)*R(1,3)
c% FC %% 	 X(I,3)= XX(I)*R(1,3)
c% FC %% 	X0(I,3)=XX0(I)*R(1,3)
c% FC %% 	 X(I,4)= XX(I)
c% FC %% 	X0(I,4)=XX0(I)
c% FC %% 12	CONTINUE
c% FC %% 
c% FC %% 	D04 = DCMPLX(0D0,0D0)
c% FC %% 	DO 13 I=1,2
c% FC %% 	DO 13 J=1,4
c% FC %% 	A1 = 1D0+X0(I,J)*S0(J) + ABS(1D0+X0(I,J)*S0(J))*IEPS*
c% FC %%      *				  SIGN(1D0,DIMAG(X(I,J)*SS(J)))
c% FC %% 	A2 = 1D0+X0(I,J)/S0(J) + ABS(1D0+X0(I,J)/S0(J))*IEPS*
c% FC %%      *				  SIGN(1D0,DIMAG(X(I,J)/SS(J)))
c% FC %% c  org. code:
c% FC %% c$$$           D04 = D04 + (-1D0)**(I+J)*(
c% FC %% c$$$     &          LI2(A1)+ETA(-X(I,J),SS(J))*LOG(A1)
c% FC %% c$$$     &          +LI2(A2)+ETA(-X(I,J),1D0/SS(J))*LOG(A2))
c% FC %% 
c% FC %%         if(abs(ETA(-X(I,J),SS(J))).ne.0d0) then
c% FC %%            D1=ETA(-X(I,J),SS(J))*LOG(A1)
c% FC %%         else
c% FC %%            D1=DCMPLX(0d0,0d0)
c% FC %%         endif
c% FC %%         if(abs(ETA(-X(I,J),1D0/SS(J))).ne.0d0) then
c% FC %%            D2=ETA(-X(I,J),1D0/SS(J))*LOG(A2)
c% FC %%         else
c% FC %%            D2=DCMPLX(0d0,0d0)
c% FC %%         endif
c% FC %%         D04=D04+(-1D0)**(I+J)*(LI2(A1)+LI2(A2)+D1+D2)
c% FC %% c   The enquiry avoids occurrence of NANs causing by LOG(A1) for A1=0 
c% FC %% c   and LOG(A2) for A2=0. At that points also ETA(-X(I,J),1D0/SS(J)) 
c% FC %% c   or ETA(-X(I,J),SS(J) are 0. (Michael Kubocz)
c% FC %% 
c% FC %% 13	CONTINUE
c% FC %% 
c% FC %% c        print*,'DIMAG(R(1,3))',DIMAG(R(1,3))
c% FC %% c	IF( DIMAG(R(1,3)).EQ.0D0 ) THEN !org (makes troubles in squark pentagons)
c% FC %% 	IF( abs(DIMAG(R(1,3))).le.eps1 ) THEN
c% FC %% 	DO 14 I=1,2
c% FC %% 	   A1 = (K(1,3)-2D0*R(1,3))/XX0(I)
c% FC %%      *		      -R(1,3)*K(1,4)+K(3,4)
c% FC %%      	   A2 = ((K(2,4)-2D0*R(2,4))*R(1,3)*XX0(I)
c% FC %%      *		      -R(2,4)*K(3,4)+K(2,3))/DD
c% FC %% 	   A3 = (K(1,3)-2D0*R(1,3))*R(2,4)/XX0(I)
c% FC %%      *		      -R(1,3)*K(1,2)+K(2,3)
c% FC %% 	   A4 = ((K(2,4)-2D0*R(2,4))*XX0(I)
c% FC %%      *		      -R(2,4)*K(1,4)+K(1,2))/DD
c% FC %% 	   L1 = LOG( A1-ABS(A1)*IEPS )
c% FC %%      	   L2 = LOG( A2+ABS(A2)*IEPS*G(I)*SIGN(1D0,DREAL(R(1,3))
c% FC %%      *				        	  *DIMAG(RS(2,4))) )
c% FC %% 	   L3 = LOG( A3-ABS(A3)*IEPS )
c% FC %% 	   L4 = LOG( A4+ABS(A4)*IEPS*G(I)*SIGN(1D0,DIMAG(RS(2,4))) )
c% FC %% 
c% FC %% c org. code:
c% FC %% c$$$	   D04 = D04 
c% FC %% c$$$     &         + (3D0-2D0*I)*(
c% FC %% c$$$     *		       ETAS( -XX(I),R(1,3),RS(1,3) )
c% FC %% c$$$     *		          *( LOG(R(1,3)*XX(I)) + L1 + L2 )
c% FC %% c$$$     *		     + ETAS( -XX(I),1D0/R(2,4),1D0/RS(2,4) )
c% FC %% c$$$     *		          *( LOG(XX(I)/R(2,4)) + L3 + L4 )
c% FC %% c$$$     *		     - ( ETAS( -XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4) )
c% FC %% c$$$     *		       + ETA( RS(1,3),1D0/RS(2,4) )                  )
c% FC %% c$$$     *		        *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 + L2 )
c% FC %% c$$$     *	  	     + ETA( RS(1,3),1D0/RS(2,4) )
c% FC %% c$$$     *		       *ETAS(-XX(I),-R(1,3)/R(2,4),-RS(1,3)/RS(2,4))   )
c% FC %% 
c% FC %%             if(abs(ETAS(-XX(I),R(1,3),RS(1,3))).ne.0d0) then
c% FC %%                D1=ETAS(-XX(I),R(1,3),RS(1,3))*(LOG(R(1,3)*XX(I))+L1+L2)
c% FC %%             else
c% FC %%                D1=DCMPLX(0d0,0d0)
c% FC %%             endif
c% FC %%             if(abs(ETAS(-XX(I),1D0/R(2,4),1D0/RS(2,4))).ne.0d0) then
c% FC %%                D2=ETAS(-XX(I),1D0/R(2,4),1D0/RS(2,4))*(LOG(XX(I)/R(2,4))
c% FC %%      &           +L3+L4)
c% FC %%             else
c% FC %%                D2=DCMPLX(0d0,0d0)
c% FC %%             endif
c% FC %%             if((abs(ETAS(-XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4))).ne.0d0)
c% FC %%      &          .or.(abs(ETA(RS(1,3),1D0/RS(2,4))).ne.0d0)) then 
c% FC %%                D3=-(ETAS(-XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4))
c% FC %%      &              +ETA(RS(1,3),1D0/RS(2,4)))*(LOG(XX(I)*R(1,3)/R(2,4))
c% FC %%      $              +L3+L2)
c% FC %%             else
c% FC %%                D3=DCMPLX(0d0,0d0)
c% FC %%             endif
c% FC %%             D4=ETA(RS(1,3),1D0/RS(2,4))*ETAS(-XX(I),-R(1,3)/R(2,4),
c% FC %%      $           -RS(1,3)/RS(2,4))
c% FC %%             D04=D04+(3D0-2D0*I)*(D1+D2+D3+D4)
c% FC %% c   The enquiry avoids occurrence of NANs causing by LOG(0). At that points 
c% FC %% c   also ETA(...) and ETAS(...) are 0. (Michael Kubocz)
c% FC %% 
c% FC %% 14	CONTINUE
c% FC %% 	ELSE
c% FC %% 	DO 15 I=1,2
c% FC %% 	   L1 = LOG( R(2,4)/XX0(I)+XX0(I)/R(2,4)+K(1,2)
c% FC %%      *		     -XX0(I)/R(2,4)*EPS*BB*G(I) )
c% FC %% 	   L2 = LOG( R(1,3)*XX0(I)+1D0/XX0(I)/R(1,3)+K(3,4)
c% FC %%      *		     -XX0(I)*R(1,3)*EPS*BB*G(I) )
c% FC %% 	   L3 = LOG( R(1,3)/R(2,4)*XX0(I)+R(2,4)/XX0(I)/R(1,3)+K(2,3)
c% FC %%      *		     -XX0(I)*R(1,3)/R(2,4)*EPS*BB*G(I) )
c% FC %% 
c% FC %% 	   D04 = D04 
c% FC %%      &          + (3D0-2D0*I)*(
c% FC %%      *		     ETA(-XX(I),1D0/R(2,4))
c% FC %%      *		      *( LOG(XX(I)/R(2,4)) + L1 )
c% FC %%      *		    +ETA(-XX(I),R(1,3))
c% FC %%      *		      *( LOG(R(1,3)*XX(I)) + L2 )
c% FC %%      *		    -( ETA(-XX(I),R(1,3)/R(2,4))
c% FC %%      *		      +ETA(R(1,3),1D0/R(2,4)) )
c% FC %%      *		      *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 )
c% FC %%      *	  	    +ETA(R(1,3),1D0/R(2,4))
c% FC %%      *		      *ETA(-XX(I),-R(1,3)/R(2,4))
c% FC %%      *		      *(1D0-G(I)*SIGN(1D0,DREAL(BB)))	    )
c% FC %% 15	CONTINUE
c% FC %% 	ENDIF
c% FC %% 
c% FC %% 	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
c% FC %% 	RETURN
c% FC %% 
c% FC %% 
c% FC %% C--->	***************** SPEZIELL ( --> T.SACK-PROMOTION )
c% FC %% C	D1=Q12-M12
c% FC %% C	D2=Q2 -M22
c% FC %% C	D3=Q3 -M42
c% FC %% C	IF ((D1*D2.LE.0D0).OR.(D2*D3.LE.0D0)) THEN
c% FC %% C	   WRITE(*,*) 'THE CASE OF DIFFERENT SIGNS OF THE D1,D2,D3'
c% FC %% C	   WRITE(*,*) 'IN D04(...) IS NOT IMPLEMENTED !!!'
c% FC %% C	   STOP
c% FC %% C	ENDIF
c% FC %% C	NM1=ABS(MM1/D1)
c% FC %% C	NM2=ABS(MM2/D2)
c% FC %% C	NM3=ABS(MM4/D3)
c% FC %% C	NP1=Q2/D2**2+Q12/D1**2+(Q1-Q2-Q12)/D1/D2
c% FC %% C	NP2=Q2/D2**2+ Q3/D3**2+(Q23-Q2-Q3)/D2/D3
c% FC %% C	NP3=Q3/D3**2+Q12/D1**2+(Q4-Q3-Q12)/D1/D3
c% FC %% C	D04=C04(NP1,NP2,NP3,NM1,NM2,NM3)/D1/D2/D3
c% FC %% 
c% FC %% C	*************** ALLGEMEIN
c% FC %% 
c% FC %% 
c% FC %% C	****** ONLY MM3 IS ZERO ******
c% FC %% 30	CONTINUE
c% FC %% 	EPS=1D-17
c% FC %% 	IEPS=DCMPLX(0D0,EPS)
c% FC %% 	M(1)=MM1
c% FC %% 	M(2)=MM2
c% FC %% 	M(3)=10D0
c% FC %% 	M(4)=MM4
c% FC %% 	P(1,2)=Q1
c% FC %% 	P(1,3)=Q12
c% FC %% 	P(1,4)=Q4
c% FC %% 	P(2,3)=Q2
c% FC %% 	P(2,4)=Q23
c% FC %% 	P(3,4)=Q3
c% FC %% 	DO 1 J=2,4
c% FC %% 	DO 1 I=1,J-1
c% FC %% 	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
c% FC %% 	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
c% FC %% 	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
c% FC %% 	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
c% FC %%      *	            DCMPLX(1D0,0D0))
c% FC %% 	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
c% FC %% 	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
c% FC %%      *	               DCMPLX(1D0,0D0))
c% FC %% 	ELSE
c% FC %% 	   RS(I,J)=R(I,J)
c% FC %% 	ENDIF
c% FC %% 1	CONTINUE
c% FC %% 	SS(1)=RS(1,2)
c% FC %% 	SS(2)=RS(2,3)
c% FC %% 	SS(3)=RS(3,4)
c% FC %% 	SS(4)=RS(1,4)
c% FC %% 	AA=K(3,4)/R(2,4)-K(2,3)
c% FC %% 	BB=K(1,3)*(1D0/R(2,4)-R(2,4))+K(1,2)*K(3,4)-K(1,4)*K(2,3)
c% FC %% 	CC=K(1,2)*K(1,3)-K(1,3)*K(1,4)*R(2,4)+R(2,4)*K(3,4)-K(2,3)
c% FC %% 	DD=K(2,3)-R(2,4)*K(3,4)
c% FC %% 	XX(1)=SQE(AA,BB,CC+IEPS*DD)
c% FC %% 	XX(2)=(CC+IEPS*DD)/AA/XX(1)
c% FC %% 	DO 2 I=1,2
c% FC %% 	X(I,1)=XX(I)/R(2,4)
c% FC %% 	X(I,2)=XX(I)/R(2,4)*R(1,3)
c% FC %% 	X(I,3)=XX(I)*R(1,3)
c% FC %% 	X(I,4)=XX(I)
c% FC %% 2	CONTINUE
c% FC %% 	D04 = DCMPLX(0D0,0D0)
c% FC %% 	DO 3 I=1,2
c% FC %% 	D04 = D04 + (2D0*I-3D0)*(
c% FC %%      *		LI2(1D0+SS(4)*X(I,4))
c% FC %%      *	       -LI2(1D0+SS(1)*X(I,1))
c% FC %%      *	       +LI2(1D0+X(I,4)/SS(4))
c% FC %%      *	       -LI2(1D0+X(I,1)/SS(1))
c% FC %%      *	       +ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
c% FC %%      *	       -ETA(-X(I,1),SS(1))*LOG(1D0+SS(1)*X(I,1))
c% FC %%      *	       +ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
c% FC %%      *	       -ETA(-X(I,1),1D0/SS(1))*LOG(1D0+X(I,1)/SS(1))
c% FC %%      *	       -LI2(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
c% FC %%      *	       +LI2(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))
c% FC %%      *	       -ETA(-X(I,4),(K(3,4)-IEPS)/(K(1,3)-IEPS))
c% FC %%      *	           *LOG(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
c% FC %%      *	       +ETA(-X(I,1),(K(2,3)-IEPS)/(K(1,3)-IEPS))
c% FC %%      *	           *LOG(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))   )
c% FC %% 	IF (DIMAG(R(2,4)).NE.0D0) THEN
c% FC %% 	   H=ETA(-1D0/XX(I),R(2,4))
c% FC %% 	ELSE
c% FC %% 	   H=DCMPLX(0D0,0D0)
c% FC %% 	   IF (DREAL(R(2,4)).LT.0D0) THEN
c% FC %% 	      HH=-1D0/XX(I)
c% FC %% 	      IM1=DIMAG(HH)
c% FC %% 	      IM2=DIMAG(RS(2,4))
c% FC %%               pi = 4.D0*datan(1.D0)
c% FC %% 	      IF ((IM1.GT.0D0).AND.(IM2.GT.0D0)) THEN
c% FC %% 	         H=-DCMPLX(0D0,2D0*PI)
c% FC %% 	      ENDIF
c% FC %% 	      IF ((IM1.LT.0D0).AND.(IM2.LT.0D0)) THEN
c% FC %% 	         H=+DCMPLX(0D0,2D0*PI)
c% FC %% 	      ENDIF
c% FC %% 	   ENDIF
c% FC %% 	ENDIF
c% FC %% 	D04 = D04 + (2D0*I-3D0)*
c% FC %%      *	          H*( LOG( (K(1,2)-R(2,4)*K(1,4)
c% FC %%      *			  +XX(I)*(1D0/R(2,4)-R(2,4)))/DD )
c% FC %%      *		     +LOG(K(1,3)-IEPS) )
c% FC %% 3	CONTINUE
c% FC %% 	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
c% FC %% 	RETURN
c% FC %% 
c% FC %% C	****** ONLY MM2 AND MM3 ARE ZERO ******
c% FC %% 40	CONTINUE
c% FC %% 	EPS=1D-17
c% FC %% 	IEPS=DCMPLX(0D0,EPS)
c% FC %% 
c% FC %% 	M(1)=MM1
c% FC %% 	M(2)=10D0
c% FC %% 	M(3)=10D0
c% FC %% 	M(4)=MM4
c% FC %% 	P(1,2)=Q1
c% FC %% 	P(1,3)=Q12
c% FC %% 	P(1,4)=Q4
c% FC %% 	P(2,3)=Q2
c% FC %% 	P(2,4)=Q23
c% FC %% 	P(3,4)=Q3
c% FC %% 	DO 4 J=2,4
c% FC %% 	DO 4 I=1,J-1
c% FC %% 	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
c% FC %% 	IF (I.EQ.2) K(I,J)=K(I,J)-M(I)/M(J)
c% FC %% 	IF (J.EQ.2) K(I,J)=K(I,J)-M(J)/M(I)
c% FC %% 	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
c% FC %% 	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
c% FC %% 	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
c% FC %%      *	            DCMPLX(1D0,0D0))
c% FC %% 	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
c% FC %% 	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
c% FC %%      *	               DCMPLX(1D0,0D0))
c% FC %% 	ELSE
c% FC %% 	   RS(I,J)=R(I,J)
c% FC %% 	ENDIF
c% FC %% 4	CONTINUE
c% FC %% 	SS(1)=RS(1,2)
c% FC %% 	SS(2)=RS(2,3)
c% FC %% 	SS(3)=RS(3,4)
c% FC %% 	SS(4)=RS(1,4)
c% FC %% 	AA=K(2,4)*K(3,4)-K(2,3)
c% FC %% 	BB=K(1,3)*K(2,4)+K(1,2)*K(3,4)-K(1,4)*K(2,3)
c% FC %% 	CC=K(1,2)*K(1,3)-K(2,3)
c% FC %% 	DD=K(2,3)
c% FC %% 	XX(1)=SQE(AA,BB,CC+IEPS*DD)
c% FC %% 	XX(2)=(CC+IEPS*DD)/AA/XX(1)
c% FC %% 	DO 5 I=1,2
c% FC %% 	X(I,1)=XX(I)/R(2,4)
c% FC %% 	X(I,2)=XX(I)/R(2,4)*R(1,3)
c% FC %% 	X(I,3)=XX(I)*R(1,3)
c% FC %% 	X(I,4)=XX(I)
c% FC %% 5	CONTINUE
c% FC %% 	D04 = DCMPLX(0D0,0D0)
c% FC %% 	DO 6 I=1,2
c% FC %% 	D04 = D04 + (2D0*I-3D0)*(
c% FC %%      *		LI2(1D0+SS(4)*X(I,4))
c% FC %%      *	       +LI2(1D0+X(I,4)/SS(4))
c% FC %%      *	       +ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
c% FC %%      *	       +ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
c% FC %%      *	       -LI2(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
c% FC %%      *	       -LI2(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
c% FC %%      *	       -ETA(-XX(I),(K(3,4)-IEPS)/(K(1,3)-IEPS))
c% FC %%      *	           *LOG(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
c% FC %%      *	       -ETA(-XX(I),(K(2,4)-IEPS)/(K(1,2)-IEPS))
c% FC %%      *	           *LOG(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
c% FC %%      *	       +LOG(-XX(I))*( LOG(K(1,2)-IEPS)
c% FC %%      *			     +LOG(K(1,3)-IEPS)-LOG(K(2,3)-IEPS) ) )
c% FC %% 6	CONTINUE
c% FC %% 	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
c% FC %% 
c% FC %% 	RETURN
c% FC %% 
c% FC %% 	END
c% FC %% 
c% FC %% ***********************************************************************
c% FC %%         FUNCTION ETAS(Y,R,RS)                                            
c% FC %% ***********************************************************************
c% FC %% *       MODIFIED ETA-FUNKTION                                           
c% FC %% *---------------------------------------------------------------------*
c% FC %% *       18.1.94   SD                                       
c% FC %% ***********************************************************************
c% FC %%         IMPLICIT   LOGICAL(A-Z)                                        
c% FC %%         COMPLEX*16 ETA,ETAS,Y,R,RS
c% FC %%         REAL*8     PI,IMY,IMRS
c% FC %%                                                                        
c% FC %%         PI     = 4D0*DATAN(1D0)                                        
c% FC %% 
c% FC %% 	IF( DIMAG(R).NE.0D0 ) THEN
c% FC %% 	    ETAS = ETA(Y,R)
c% FC %% 	ELSE	    
c% FC %% 	    IF( DREAL(R).GT.0D0 ) THEN
c% FC %% 		ETAS = DCMPLX(0D0,0D0)
c% FC %% 	    ELSE
c% FC %% 	 	IMY  = DIMAG(Y)
c% FC %% 		IMRS = DIMAG(RS)
c% FC %% 		ETAS = 2D0*DCMPLX(0D0,PI)*(
c% FC %%      *			(1D0+SIGN(1D0,-IMY))*(1D0+SIGN(1D0,-IMRS))-
c% FC %%      *			(1D0+SIGN(1D0, IMY))*(1D0+SIGN(1D0, IMRS))
c% FC %%      *					  )/4D0
c% FC %% 	    ENDIF
c% FC %% 	ENDIF
c% FC %%         END                                                            

c***********************************************************************
c        FUNCTION SQE(A,B,C)                                            
c***********************************************************************
c*       SOLUTION OF QUADRATIC EQUATION				      *
c*---------------------------------------------------------------------*
c*       13.1.92  SD						      *
c***********************************************************************
c        IMPLICIT REAL*8 (A-Z)                                        
c        COMPLEX*16 A,B,C,SQE,X1,X2
c
c	X1=(-B+SQRT(B**2-4D0*A*C))/2D0/A
c	X2=(-B-SQRT(B**2-4D0*A*C))/2D0/A
c
c	IF (ABS(X1).GT.ABS(X2)) THEN
c	   SQE=X1
c	ELSE
c	   SQE=X2
c	ENDIF
c
c        END                                                            
c

