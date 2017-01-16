c
c-------------  spence function  li2(z) ------------------------
c
!      real*16 function li2QUAD(zin)
      function dilog_QUAD(zin)
      implicit none
      real*16 dilog_QUAD
      real*16 zin
      real*16  z, ans, zext,test
      real*16 fac
    

c
c determine the value of the dilogarithm 
c
c    li2(z) = - int_0^1  log(1-zt)/t dt  with cut along the positive 
c                                        real axis, z>1
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2000 November 6
c	Last modified:    2000 November 12
c       Modified: Francisco Campanario to make it Quatruple precision
c       Real argument only! Including c24->c38
c       Relation for negative argument
c !!!ATention only valid for the obtation of real part of Polylog.
c Otherwise, the analytical continuation is wrong!!
c Use instead the complex version: function li2_QUAD(zin)
     
      real*16  c0,c1,c2,c4,c6,c8,c10,c12,c14,c16,c18,c20,c22
      real*16  c24,c26,c28,c30,c32,c34,c36,c38
      real*16  b0,b1,b2,b4,b6,b8,b10,b12,b14,b16,b18,b20,b22
      real*16  b24,b26,b28,b30,b32,b34,b36,b38
      real*16  d0,d1,d2,d4,d6,d8,d10,d12,d14,d16,d18,d20,d22
      real*16  d24,d26,d28,d30,d32,d34,d36,d38
      real*16  I2,I4,I9,I16,I25,I36,I49,I64,I81,I100,I121,I144
      REAL*16  I169
      real*16  t,t1,t2,t4,t8,t16,t32
  
      parameter (b0=1q0,            d0 =1q0,      c0= b0/d0,I2=1q0/2q0, I4=1/4q0 )
      parameter (b1=-1q0/2q0,       d1 =2q0,      c1= b1/d1,I9=1q0/9q0)
      parameter (b2= 1q0/6q0,       d2 =6q0,      c2= b2/d2,I16=1q0/16q0)
      parameter (b4=-1q0/30q0,      d4 =12q1,     c4= b4/d4,I25=1q0/25q0)
      parameter (b6=1q0/42q0,       d6 =504q1,    c6= b6/d6,I36=1q0/36q0)
      parameter (b8=-1q0/30q0,      d8 =36288q1,  c8= b8/d8,I49=1q0/49q0)
      parameter (b10=5q0/66q0,      d10=399168q2,  c10=b10/d10,I64=1q0/64q0)
      parameter (b12=-691q0/2730q0, d12=62270208q2, c12=b12/d12,I81=1q0/81q0)
      parameter (b14=7q0/6q0,       d14=1307674368q3, c14=b14/d14,I100=1q0/100q0)
      parameter (b16=-3617q0/510q0, d16=355687428096q3, c16=b16/d16,I121=1q0/121q0)
      parameter (b18=43867q0/798q0, d18=121645100408832q3, c18=b18/d18,I144=1q0/144q0)
      parameter (b20=-174611q0/330q0,d20=5109094217170944q4, c20=b20/d20)
      parameter (b22=854513q0/138q0,d22=25 852 016 738 884 976 64q4,  c22=b22/d22)

      parameter (b24=-(236364091q0/2730q0),d24=15 511 210 043 330 985 984q6,  c24=b24/d24)
      parameter (b26=8553103q0/6q0,d26=10 888 869 450 418 352 160 768q6, c26=b26/d26)
      parameter (b28=-23749461029q0/870q0
     - ,d28=8 841 761 993 739 701 954 543 616q6, c28=b28/d28)
       parameter (b30=8615841276005q0/14322q0,
     - d30=8 222 838 654 177 922 817 725 562 88q7,  c30=b30/d30)
      parameter (b32=-7709321041217q0/510q0,
     -   d32=8 683 317 618 811 886 495 518 194 401 28q7,  c32=b32/d32)
      parameter (b34=2577687858367q0/6q0,
     -  d34=10 333 147 966 386 144 929 666 651 337 523 2q8, c34=b34/d34)
      parameter (b36=-26315271553053477373q0/1919190q0
     - ,d36=13 763 753 091 226 345 046 315 979 581 580 902 4q8,  
     -  c36=b36/d36)
      parameter (b38=2929993913841559q0/6q0,
     - d38=20 397 882 081 197 443 358 640 281 739 902 897 356 8q8,  
     - c38=b38/d38)
      parameter (I169=1q0/169q0)


      real*16 eps, epst, pi, pi2o6,pi2o3
      parameter (eps=1q-8, epst=1q-3)
      parameter (pi=3.141 592 653 589 793 238 462 643 383 279 502 88q0
     &  , pi2o6=1.644 934 066 848 226 436 472 415 166 646 025 19q0, 
     &    pi2o3=3.289 868 133 696 452 872 944 830 333 292 050 38q0)
c
c debug information
c      logical ldebug
c      parameter (ldebug=.false.)

      z = zin
c      print*,' li2 call with z = ',z
     
      if (z.eq.1) then
      dilog_QUAD=pi2o6
      return
      endif
       

      if (abs(z).lt.eps) then 
         dilog_QUAD =((((I25*z+I16)*z+I9)*z+I4)*z+1q0)*z 
         return
      elseif (abs(z).lE.epst) then
       ans=((((((((((((I169*z+I144)*z+I121)*z+I100)*z+I81)*z+I64)*z+I49)*z+I36)*z+I25)*z+I16)*z+I9)*z+I4)*z+1q0)*z 
      dilog_QUAD = ans
      return
      endif
     
c
c use z-->1/z and z--> 1-z mappings of the spence function to restrict 
c agument to unit circle in the complex plane with Re(z) <= 0.5
c
    
      fac = 1q0
      zext=0q0


      if(z.lt.0q0) then  
        t=Log(-z)
        t1=Log(1Q0-z)
        zext = t1*(t1*I2-t) -pi2o6
         z = 1q0/(1q0-z)
      endif

      if (z.gt.1q0) then     ! map z ---> 1/z
c !!! Here: ATention only valid for the obtation of real part of Polylog.
c Otherwise, the analytical continuation is wrong in this step!!
c Use instead the complex version: function li2_QUAD(zin)
         t=Log(z)
         zext =zext+ (pi2o3 - I2*t*t)
         z = 1q0/z
       fac=-fac
      endif

      if (z.gt.0.5q0) then     ! map new z ---> 1-z
         zext = zext + fac*(pi2o6-Log(z)*Log(1q0-z))
         z = 1q0-z
        fac=-fac
      endif


c
c now use t = 1 - exp(-u) mapping to write Li(z) in terms of Bernoulli 
c numbers
c

      t = -Log(1q0-z)
      t2= t*t
      t4=t2*t2 
      t8=t4*t4
      t16=t8*t8
      t32=t16*t16

      ans = ((c2*t+c1)*t+c0)*t
      ans = (((((c14*t2+c12)*t2+c10)*t2+c8)*t2+c6)*t2+c4)*t4*t+ans 
      test=c16*t16*t
      ans=((((((c30*t2+c28)*t2+c26)*t2+c24)*t2+c22)*t2+c20)*t2+c18)*t2*t*t16+test+ans
      test=c32*t32*t
      ans=((c38*t2+c36)*t2+c34)*t2*t32*t+test+ans
      dilog_QUAD = fac*ans + zext
      Return
      end
