    
c
c-------------  spence function  li2(z) ------------------------
c
      complex*32 function li2_QUAD(zin)
      implicit none
      complex*32 zin, z, unpo, ans, zext,test
      real*16 r, r2, r2n, fac
c
c determine the value of the dilogarithm 
c
c    li2(z) = - int_0^1  log(1-zt)/t dt  with cut along the positive 
c                                        real axis, z>1
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2000 November 6
c	Last modified:    2000 November 12
c  
      integer i
      real*16 c0,c1,c2,c4,c6,c8,c10,c12,c14,c16,c18,c20,c22
      real*16  c24,c26,c28,c30,c32,c34,c36,c38
      real*16 b0,b1,b2,b4,b6,b8,b10,b12,b14,b16,b18,b20,b22
      real*16  b24,b26,b28,b30,b32,b34,b36,b38
      real*16 d0,d1,d2,d4,d6,d8,d10,d12,d14,d16,d18,d20,d22
      real*16  d24,d26,d28,d30,d32,d34,d36,d38
     
      real*16  I2,I4,I9,I16,I25,I36,I49,I64,I81,I100,I121,I144
      real*16 I169
      complex*32  u,u2,u4,u6,u8,u16,u32
     

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
      parameter (b38=2929993913841559/6q0,
     - d38=20 397 882 081 197 443 358 640 281 739 902 897 356 8q8,  
     - c38=b38/d38)
      parameter (I169=1q0/169q0)

      real*16 eps,eps1, epst, pi, pi2o6,pi2o3
      parameter (eps=1q-8, epst=1q-6,eps1=1q-70)
      parameter (pi=3.141 592 653 589 793 238 462 643 383 279 502 88q0
     &  , pi2o6=1.644 934 066 848 226 436 472 415 166 646 025 19q0, 
     &    pi2o3=3.289 868 133 696 452 872 944 830 333 292 050 38q0)

      real*16 zreal,zimag
      complex*32 zlog,zlog1
c
c debug information
      logical ldebug
      parameter (ldebug=.false.)

      z = zin
c      print*,' li2 call with z = ',z
      u = z*z
      zreal=QREAL(z)
      zimag=QIMAG(z)

      r2 = zreal*zreal + zimag*zimag

      if (r2.lT.eps) then
         li2_QUAD = ((((I25*z+I16)*z+I9)*z+I4)*z+1q0)*z 
         return
      elseif (r2.lE.epst) then

         ans=((((((((((((I169*z+I144)*z+I121)*z+I100)*z+I81)*z+I64)*z+I49)*z+I36)
     c     *z+I25)*z+I16)*z+I9)*z+I4)*z+1q0)*z 
         li2_QUAD=ans
         return
      endif
      if (zreal.ge.1q0 .and. zimag.eq.0q0 ) then
         z = z + (0q0,1q0)*eps1
      endif
c
c use z-->1/z and z--> 1-z mappings of the spence function to restrict 
c agument to unit circle in the complex plane with Re(z) <= 0.5
c
      zext = (0q0,0q0)
      fac = 1q0

      if(zREAL.lt.0q0) then  
       ZLOG=Log(-z)
       zlog1=Log(1Q0-z)
       zext = zlog1*(zlog1*I2-zlog) -pi2o6
       z = 1q0/(1q0-z)
       zreal=QREAL(z)
       zimag=QIMAG(z)
      r2 = zreal*zreal + zimag*zimag
      endif




        if (sqrt(r2).gt.1q0) then     ! map z ---> 1/z
         fac = -fac
         ZLOG=Log(-z)
         zext = zext  -(pi2o6 + I2*zlog*zlog)
         z = 1q0/z
        zreal=QREAL(z)
        zimag=QIMAG(z)
      endif



      if (zreal.gt.I2) then     ! map new z ---> 1-z
         zlog=log(z)
         zlog1=log(1q0-z)
         zext = zext + fac*(pi2o6-zlog*zlog1)
         fac = -fac
         z = 1q0-z
      endif
c
c now use t = 1 - exp(-u) mapping to write Li(z) in terms of Bernoulli 
c numbers
c
      u = - log(1q0-z)
      u2 = u*u
      u4=u2*u2
      u8=u4*u4
      u16=u8*u8
      u32=u16*u16

      ans = ((c2*u+c1)*u+c0)*u
      ans = (((((c14*u2+c12)*u2+c10)*u2+c8)*u2+c6)*u2+c4)*u4*u+ans
      test=c16*u16*u
      ans=((((((c30*u2+c28)*u2+c26)*u2+c24)*u2+c22)*u2+c20)*u2+c18)*u2*u*u16+
     c test+ans
      test=c32*u32*u
      ans=((c38*u2+c36)*u2+c34)*u2*u32*u+test+ans
      li2_QUAD = fac*ans + zext
      Return
      end






