      subroutine InitGaugeTest
      implicit none 
#include "GaugeTest.inc"
      real*8 Accuracy
      Accuracy =1d-2
c      print*, "Accuracy",Accuracy
      GaugeAccuracyBoxL=Accuracy
      GaugeAccuracyBoxLnA=Accuracy
      GaugeAccuracyPenL=Accuracy
      GaugeAccuracyPenLnA=Accuracy
      GaugeAccuracyHexL=Accuracy
      GaugeAccuracyHexLnA=Accuracy
c      print*, "GaugeAccuracyBoxL",GaugeAccuracyBoxL      
      end
      

c###################################################################
      subroutine pent_gauge(p1,p2,p3,p4,p5,barpsi_p5,psi_p1,
     -  mup2,mup3,mup4,alpha,musq,result,resultb,Ward)
c 
c Author: Francisco Campanario
C Date: 16/01/2009
C Compute the gauge three gauge test for the pentagons and result .True. if the max of it
c is larger that GaugeAccuracyPenL. In the same subroutines the boxlines called are checked as well.


      IMPLICIT NONE
#include "GaugeTest.inc"
      DOUBLE PRECISION   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p23(0:3),p34(0:3)
      DOUBLE COMPLEX barpsi_p5(2),psi_p1(2),mup2(0:3),mup3(0:3),mup4(0:3) 

      DOUBLE PRECISION musq,test1,test2
      DOUBLE COMPLEX cte
      Integer alpha,i
      DOUBLE COMPLEX result(3),result1,result2,result3,result4
      DOUBLE COMPLEX resultb(3),resultb1,resultb2,resultb3,resultb4
      DOUBLE COMPLEX resultgauge(3),resultgaugeb(3)
      DOUBLE PRECISION pi,pi2,test(3),t
      parameter (pi =3.14159265358979323846264338327950288d0)
      parameter (pi2 =9.86960440108935861883449099987615114d0)
      external dotrr,theta
      DOUBLE PRECISION dotrr,s15,lm,ls,ths,theta
      Logical Ward
      real*8 ratiocheckCan
      common/pengaugeValue/t

      
      t=1d9
      test1=1d9
      test2= 1d9

      Ward=.False.


      do i=0,3
c      emup2(i)=mup2(i)-x2*(-p2(i))
c      emup3(i)=mup3(i)-x3*(-p3(i))
c      emup4(i)=mup4(i)-x4*(-p4(i))
      p34(i)=p3(i)+p4(i)
      p23(i)=p2(i)+p3(i)
      enddo

c           print*, 'zepp', emup2
c           print*, 'xze', x2
c           print*, 'pe', p2(0)
c           print*, 'zmupp' ,emup3
c           print*, 'ataup', emup4
             
c        ctedieter = (3*pi*(0d0,1d0)-pi**2)
c        cte=-(3*pi*(0d0,1d0)-4.0d0*pi**2/3.0d0+7.0d0)
 
        s15=2d0*dotrr(p1,p5)+dotrr(p1,p1)+dotrr(p5,p5)
        lm=Log(musq)
        ls=Log(Abs(s15))
        ths=theta(s15)

c       cte = -cte_dieter but general dependency in musq.
c         The sign is later on corrected in his code in
c          -dotcr(j,p)
        cte= -3d0*(Lm - Ls + (0d0,1d0)*pi*ths) 
     -  - (Lm - Ls + (0d0,1d0)*pi*ths)**2
c   substracted the c_virt
     -  -7d0 + 1d0/3d0*pi2

c box's 2 for contraction 3, 2
      call boxlinemm(p1,p2,p34,p5,barpsi_p5,psi_p1,mup2,mup3
     - ,alpha,musq,3,1,resultgauge,result1,resultgaugeb,resultb1)

c variables for gauge test:test1,test2      
 

      if(abs(resultgaugeb(1)).lt.1d-4) then
      test1=abs(resultgauge(1)- resultgaugeb(1)*(cte))
      else
      test1=abs(resultgauge(1)/(resultgaugeb(1)*(cte))-1d0)
      endif

      if(abs(resultgaugeb(2)).lt.1d-4) then
       test2=abs(resultgauge(2)- (resultgaugeb(2)*(cte)))
      else   
      test2=abs(resultgauge(2)/(resultgaugeb(2)*(cte))-1d0)
      endif

c test NAN

      if(test1+1.eq.test1)then
         test1=1d9
       endif


      if(test2+1.eq.test2)then
         test2=1d9
       endif


      if(test1.ne.test1)then
         test1=1d9
       endif

      if(test2.ne.test2)then
         test2=1d9
       endif


c      test1=abs(resultgauge(1)/(resultgaugeb(1)*(cte))-1d0)
c      test2=abs(resultgauge(2)/(resultgaugeb(2)*(cte))-1d0)


      if(test1.gt.GaugeAccuracyPenL.or.test2.gt.GaugeAccuracyPenL) return ! Ward=false

!      if (test1.lt.1d-2) then
 !     if (test2.lt.1d-2) then
      
      call boxlinemm(p1,p2,p34,p5,barpsi_p5,psi_p1,mup2,mup4
     - ,alpha,musq,9,-1,resultgauge,result2,resultgauge,resultb2)

c box's 1 for contraction 2,1
      
      call boxlinemm(p1,p23,p4,p5,barpsi_p5,psi_p1,mup2,mup4
     - ,alpha,musq,3,1,resultgauge,result3,resultgaugeb,resultb3)

c gauge test: test1,test2
 
      if(abs(resultgaugeb(1)).lt.1d-4) then
      test1=abs(resultgauge(1)- resultgaugeb(1)*(cte))
      else
      test1=abs(resultgauge(1)/(resultgaugeb(1)*(cte))-1d0)
      endif

      if(abs(resultgaugeb(2)).lt.1d-4) then
       test2=abs(resultgauge(2)- (resultgaugeb(2)*(cte)))
      else   
      test2=abs(resultgauge(2)/(resultgaugeb(2)*(cte))-1d0)
      endif

c test NAN

      if(test1+1.eq.test1)then
         test1=1d9
       endif

      if(test2+1.eq.test2)then
         test2=1d9
       endif

      if(test1.ne.test1)then
         test1=1d9
       endif
      if(test2.ne.test2)then
         test2=1d9
       endif


      if(test1.gt.GaugeAccuracyPenL.or.test2.gt.GaugeAccuracyPenL) return ! Ward=false
c      test1=abs(resultgauge(1)/(resultgaugeb(1)*(cte))-1d0)
c      test2=abs(resultgauge(2)/(resultgaugeb(2)*(cte))-1d0)
    
      
!      if (test1.lt.1d-2) then
!      if  (test2.lt.1d-2) then

      call boxlinemm(p1,p23,p4,p5,barpsi_p5,psi_p1,mup3,mup4
     - ,alpha,musq,9,-1,resultgauge,result4,resultgaugeb,resultb4)


      if(abs(resultb4).ne.0d0) then
	ratiocheckCan=Abs(-(result4)/(-cte*(resultb4))-1d0)


  !       print*,'r',ratiocheckCan,cte,resultb4,result4
	if(ratiocheckCan.gt.1d-2) then

	test(1)=abs((result4-cte*resultb4)/(result(1)-cte*resultb(1))+1d0)
	else
	test(1)=abs((result4)/(result(1))+1d0)
	endif
      else
	test(1) = 0d0
      endif


      if(abs(resultb(2)).ne.0d0) then
	ratiocheckCan=Abs(-result(2)/(-cte*resultb(2))-1d0)
  !       print*,'r',ratiocheckCan,cte,resultb(2),result(2)
	if(ratiocheckCan.gt.1d-2) then
	test(2)=abs(((result2-result3)-cte*(resultb2-resultb3))/(result(2)-cte*resultb(2))+1d0)
	
	else
	test(2)=abs((result2-result3)/result(2)+1d0)
	endif
      else
	test(2) = 0d0
      endif
      

      if(abs(resultb1).ne.0d0) then
	ratiocheckCan=Abs(-(-result1)/(-cte*(-resultb1))-1d0)
  !       print*,'r',ratiocheckCan,cte,resultb1,result1
      

	if(ratiocheckCan.gt.1d-2) then
	
	test(3)=abs(((-result1)-cte*(-resultb1))/(result(3)-cte*resultb(3))+1d0)
	else
	test(3)=abs(-result1/result(3)+1d0)
	endif
      else
	test(3) = 0d0
      endif

c test NAN

      do i=1,3
      if(test(i)+1.eq.test(i))then
         test(i)=1d9
      endif

      if(test(i).ne.test(i))then
         test(i)=1d9
      endif

      enddo

      t=Max(test(1),test(2),test(3))
      If(t.lt.GaugeAccuracyPenL) Ward=.True.


      Return
!      endif
!      endif
!      endif
!      endif
      Return
      end
c###################################################################
      subroutine box_gauge(p1,p2,p3,p4,musq,result,resultb,ward)
C Compute the two gauge tests for the boxes and result. .True. if the max of it
c is smaller than GaugeAccuracyPenL. 

      IMPLICIT NONE
#include "GaugeTest.inc"
      DOUBLE PRECISION p1(0:3), p2(0:3), p3(0:3), p4(0:3) 
      DOUBLE COMPLEX result(2), resultb(2)
      DOUBLE PRECISION musq, testb(2)
      DOUBLE COMPLEX cte
      DOUBLE PRECISION dotrr,pi,pi2
      EXTERNAL dotrr
      PARAMETER (pi =3.14159265358979323846264338327950288d0)
      PARAMETER (pi2 =9.86960440108935861883449099987615114d0)
      DOUBLE PRECISION s14,lm,ls,ths,theta
      LOGICAL ward
      
      ward=.False.

      s14=2d0*dotrr(p1,p4)+dotrr(p1,p1)+dotrr(p4,p4)
      lm=Log(musq)
      ls=Log(Abs(s14))
      ths=theta(s14)

c     cte = -cte_dieter but general dependency in musq.
c       The sign is later on corrected in his code in
c        -dotcr(j,p)
      cte= -3d0*(Lm - Ls + (0d0,1d0)*pi*ths) 
     -     - (Lm - Ls + (0d0,1d0)*pi*ths)**2
c   substracted the c_virt
     -     -7d0 + 1d0/3d0*pi2

      IF(ABS(RESULTB(1)).LT.1D-4) THEN
      TESTB(1)=ABS(RESULT(1)-CTE*RESULTB(1))
      ELSE
      testb(1)=abs(result(1)/(cte*resultb(1))-1d0)
      ENDIF
      
      if(ABS(RESULTB(2)).LT.1D-4) THEN
      TESTB(2)=ABS(RESULT(2)-CTE*RESULTB(2))
      ELSE
      testb(2)=abs(result(2)/(cte*resultb(2))-1d0)
      ENDIF
      
      if ( (testb(1).lt.GaugeAccuracyBoxL) .and. (testb(2).lt.GaugeAccuracyBoxL) ) then
      ward = .true.
      else
c      print*,'MAXtestb(1)', MAX(testb(1),testb(2))
c      print*,'testb(1)', testb(1)
c      print*,'testb(2)', testb(2)
c      print*, "cte",cte
c      print*, 'resultb(1)',resultb(1)
c      print*, 'resultb(2)',resultb(2)  
c      print*, 'result(1)',result(1)      
c      print*, 'result(2)',result(2)      
      endif

c      Print*, 'Ward', Ward

      return
      end

********************************************************************************


      subroutine pentNoAbe_gauge(Accuracy,p1,p2,p3,p4,p5,barpsi_p4,psi_p1,
     -  mup2,mup3,mup5,alpha,musq,
     - result,
     - td,Ward)
c
c Author: Francisco Campanario
C Date: 16/02/2009
C Compute the two gauge test for the pentlineNoAbe and returns  .True.
c if the max of it
c is lower than Accuracy. In the same subroutines the boxline(/NoAbe)
called are checked as well with a precision (Accuracy*1d-1)
C result is the result of penlineNoAbemmm for the replacement of mu2_>p2
c and mup3->p3.
c td is the result of largest gauge test which should be lower than
c Accuracy to pass the gauge test.

      IMPLICIT NONE
      double precision   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),pp23(0:3)
      double complex barpsi_p4(2),psi_p1(2),mup2(0:3),mup3(0:3),mup5(0:3)
      double complex mup5test(0:3)

      double precision musq,test1,test2
      double complex cte
      Integer alpha,i
      double complex result(3),result1,result3
      double complex resultb1,resultb3


      double complex resultgauge(3),resultgaugeb(3)
      double complex resultbtemp(5)
      double precision pi,pi2,test(3),t,td,tq,qt
      parameter (pi =3.14159265358979323846264338327950288d0)
      parameter (pi2 =9.86960440108935861883449099987615114d0)
      external dotrr,theta
      double precision dotrr,s14,lm,ls,ths,theta
      Logical Ward,WardB
      double precision Accuracy,Accuracy1
      integer countb1,countb2,countp1,countp2,counth1,counth2
      double complex resultboxNoAbe1,resultboxNoAbe2,resultboxNoAbeG(2)
      double precision btd
      common/gaugetestV/countb1,countb2,countp1,countp2,counth1,counth2  

c      print*, 'HERE BOX'
      Ward=.False.
      WardB=.False. 

      t=1.d8
      td=1.d8
      tq=1.d8
      qt=1.d8
      btd=1.d8

      Accuracy1=Accuracy*1.d-1


      do i=0,3

      pp23(i)=p2(i)+p3(i)
      mup5test(i)=p5(i)
      enddo


       s14=2d0*dotrr(p1,p4)+dotrr(p1,p1)+dotrr(p4,p4)
       lm=Log(musq)
       ls=Log(Abs(s14))
       ths=theta(s14)

       cte= -3d0*(Lm - Ls + (0d0,1d0)*pi*ths)
     -  - (Lm - Ls + (0d0,1d0)*pi*ths)**2
c    substracted the c_virt
     -  -7d0 + 1d0/3d0*pi2


ccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc
c box's 1 for contraction  2,1, (permutation 1)
ccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc

      call boxlineABEmm(p1,pp23,p5,p4,barpsi_p4,psi_p1,mup3,mup5test
     - ,alpha,musq,1,3,3,1,resultgauge,resultbtemp,resultgaugeb,resultb1)


      result1=resultbtemp(2)
c variables for gauge test:test1,test2


      if(abs(resultgaugeb(1)).lt.1d-7) then
      test1=abs(resultgauge(1)- resultgaugeb(1)*(cte))
      else
      test1=abs(resultgauge(1)/(resultgaugeb(1)*(cte))-1d0)
      endif

      if(abs(resultgaugeb(2)).lt.1d-7) then
      test2=abs(resultgauge(2)- (resultgaugeb(2)*(cte)))
      else
      test2=abs(resultgauge(2)/(resultgaugeb(2)*(cte))-1d0)
      endif


c test NAN

      if(test1+1.eq.test1)then
         test1=1d9
       endif

      if(test2+1.eq.test2)then
         test2=1d9
       endif

      if(test1.ne.test1)then
         test1=1d9
       endif

      if(test2.ne.test2)then
         test2=1d9
       endif


      t=max(test1,test2)
c$$$      Print*,"resultgauge(1)",resultgauge(1)
c$$$      Print*,"resultgauge(2)",resultgauge(2)
c$$$      Print*,"resultgaugeb(1)",resultgaugeb(1)
c$$$      Print*,"resultgaugeb(2)",resultgaugeb(2)
c$$$      Print*,"resultgaugeb(1)cte",cte*resultgaugeb(1)
c$$$      Print*,"resultgaugeb(2)cte",cte*resultgaugeb(2)
c$$$      Print*,"test1",test1
c$$$      Print*,"test2",test2
c$$$      
c$$$      Print*,"t1",t


      IF (T.LT.Accuracy1) then
      WardB=.True.
      else
      countb1=countb1+1
      endif

      if (WardB) then
      WardB=.False.
ccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc
c box's 1 for contraction 2,1 (permutation 2)
ccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccc
      call boxlineabemm(p1,p5,pp23,p4,barpsi_p4,psi_p1,mup5test,mup3
     - ,alpha,musq,1,2,3,1,resultgauge,resultbtemp,resultgaugeb,resultb3)

      result3=resultbtemp(2)

c gauge test: test1,test2

      if(abs(resultgaugeb(1)).lt.1d-7) then
      test1=abs(resultgauge(1)- resultgaugeb(1)*(cte))
      else
      test1=abs(resultgauge(1)/(resultgaugeb(1)*(cte))-1d0)
      endif

      if(abs(resultgaugeb(2)).lt.1d-7) then
      test2=abs(resultgauge(2)- (resultgaugeb(2)*(cte)))
      else
      test2=abs(resultgauge(2)/(resultgaugeb(2)*(cte))-1d0)
      endif

c test NAN

      if(test1+1.eq.test1)then
         test1=1d9
       endif

      if(test2+1.eq.test2)then
         test2=1d9
       endif
      if(test1.ne.test1)then
         test1=1d9
       endif

      if(test2.ne.test2)then
         test2=1d9
       endif
      t=max(test1,test2)
c$$$      Print*,"t2",t
      IF (T.LT.Accuracy1) then
      WardB=.True.
      else
      countb1=countb1+1
      endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCC          No ABe box             CCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if(wardb) then
      wardb=.false.

       call boxlineNoABEmm(p1,pp23,p4,p5,barpsi_p4,psi_p1,mup3,mup5,
     - alpha,musq,7,1,resultBoxNoAbeG,resultbtemp(1))

      resultBoxNoAbe1=resultbtemp(1)

      call boxlineNoABEmm(p1,pp23,p4,p5,barpsi_p4,psi_p1,mup3,mup5test,
     - alpha,musq,3,-1,resultbtemp,resultBoxNoAbeG(1))

      resultBoxNoAbeG(1)=resultBoxNoAbeG(1)


      if(abs(resultBoxNoAbeG(1)).lt.1d-7) then
      btd=abs(-(result1+result3)/2d0- resultBoxNoAbeG(1))
      else
      btd=abs(-(result1+result3)/2d0/resultBoxNoAbeG(1)+1d0)
      endif


c test NAN

      if(btd+1.eq.btd)then
         btd=1d9
       endif

       if(btd.ne.btd)then
         btd=1d9
       endif


      if(btd.lt.Accuracy1) wardb=.true.


      if (wardb) then
      call boxlineNoABEmm(p1,pp23,p4,p5,barpsi_p4,psi_p1,mup2,mup5,
     - alpha,musq,7,1,resultBoxNoAbeG,resultbtemp(1))

      resultBoxNoAbe2=resultbtemp(1)
      endif


ccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc
cccccccc    Pent gauge  cccc
ccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc

      if (WardB) then

      if(abs(result(1)).lt.1d-7) then
      test(1)=abs((resultBoxNoAbe1)-(result(1)))
      else
      test(1)=abs((resultBoxNoAbe1)/(result(1))+1d0)
      endif
      if(abs(result(2)).lt.1d-7) then
      test(2)=abs(((-resultBoxNoAbe2))-(result(2)))
      else
      test(2)=abs(((-resultBoxNoAbe2))/(result(2))+1d0)
      endif

c test NAN

      do i=1,2
      if(test(i)+1.eq.test(i))then
         test(i)=1d9
      endif
      if(test(i).ne.test(i))then
         test(i)=1d9
      endif

      enddo

      td=Max(test(1),test(2))


      If(td.lt.Accuracy) then
      Ward=.True.
      else
c$$$      print*, 'HERE'
c$$$      print*, 'td',td
c$$$      print*, result
c$$$      print*, resultQ
c$$$      print*, Qresult
      endif

      endif
      endif
      endif


      end


c################################################################################
      subroutine hexgaugeABE(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,
     -  mup2,mup3,mup4,mup5,alpha,musq,result,resultb,WardH,gauge)
c 
c Author: Francisco Campanario
C Date: 16/01/2009
C Compute the gauge three gauge test for HexagonLine and result .True. if the max of it
c is larger that GaugeAccuracyHexL. In the same subroutines the penlines are checked as well.
      IMPLICIT NONE
#include "GaugeTest.inc"
      Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3),p23(0:3),
     - p34(0:3),p45(0:3)
      Complex*16 barpsi_p6(2),psi_p1(2),mup2(0:3),mup3(0:3),mup4(0:3),
     - mup5(0:3) 
      Complex*16 result1,result2,result3,result4,result5,result6
      Complex*16 resultb1,resultb2,resultb3,resultb4,resultb5,resultb6
      Complex*16 resultgauge(3),resultgaugeb(3)
      Complex*16 result(4),resultb(4)
      real*8 pi,pi2,test(4),t
      parameter (pi =3.14159265358979323846264338327950288d0)
      parameter (pi2 =9.86960440108935861883449099987615114d0)
      external dotrr,theta
      real*8 dotrr,s15,lm,ls,ths,theta,musq
      complex*16 cte
      Integer i,alpha,gauge
      common/GaugeVar/WardP
      Logical WardH,WardP
      logical printGaugePrec
      parameter (printGaugePrec = .false.)
      real*8 comGaugePrec
      common /VVJJGaugePrec/ comGaugePrec
      common/ hexgaugeValue/t
 
      t=1d9
      do i=1,4
            test(i) =1d9
      enddo

      WardH=.False.
      WardP=.False.


      do i=0,3
      p23(i)=p2(i)+p3(i)
      p34(i)=p3(i)+p4(i)
      p45(i)=p4(i)+p5(i)
      enddo

        s15=2d0*dotrr(p1,p5)+dotrr(p1,p1)+dotrr(p5,p5)
        lm=Log(musq)
        ls=Log(Abs(s15))
        ths=theta(s15)

c       cte = -cte_dieter but general dependency in musq.
c         The sign is later on corrected in his code in
c          -dotcr(j,p)
        cte= -3d0*(Lm - Ls + (0d0,1d0)*pi*ths) 
     -  - (Lm - Ls + (0d0,1d0)*pi*ths)**2
c   substracted the c_virt
     -  -7d0 + 1d0/3d0*pi2

        if((gauge.eq.4).or.(gauge.eq.3).or.(gauge.gt.4))then

c pent's 3 for contraction 4, 3
      call penlinemmm(p1,p2,p3,p45,p6,barpsi_p6,psi_p1,mup2,mup3,mup4,
     - alpha,musq,4,1,resultgauge,result1,resultgaugeb,resultb1)


      call pent_gauge(p1,p2,p3,p45,p6,barpsi_p6,psi_p1,
     -  mup2,mup3,mup4,alpha,musq,resultgauge,resultgaugeb,WardP)

 
!      if (WardP) then

      call penlinemmm(p1,p2,p3,p45,p6,barpsi_p6,psi_p1,mup2,mup3,mup5,
     - alpha,musq,9,-1,resultgauge,result2,resultgauge,resultb2)

         endif

       if((gauge.eq.2).or.(gauge.eq.3).or.(gauge.gt.4))then


c pen's 2 for contraction 2,3
      
      call penlinemmm(p1,p2,p34,p5,p6,barpsi_p6,psi_p1,mup2,mup3,mup5,
     - alpha,musq,4,1,resultgauge,result3,resultgaugeb,resultb3)

      call pent_gauge(p1,p2,p34,p5,p6,barpsi_p6,psi_p1,
     -  mup2,mup3,mup5,alpha,musq,resultgauge,resultgaugeb,WardP)

      
!      if (WardP) then

      call penlinemmm(p1,p2,p34,p5,p6,barpsi_p6,psi_p1,mup2,mup4,mup5,
     -  alpha,musq,9,-1,resultgauge,result4,resultgaugeb,resultb4)

        endif

       if((gauge.eq.1).or.(gauge.eq.2).or.(gauge.gt.4))then

c pen's 1 for contraction 1,2
      
      call penlinemmm(p1,p23,p4,p5,p6,barpsi_p6,psi_p1,mup2,mup4,mup5,
     - alpha,musq,4,1,resultgauge,result5,resultgaugeb,resultb5)

      call pent_gauge(p1,p23,p4,p5,p6,barpsi_p6,psi_p1,
     -  mup2,mup4,mup5,alpha,musq,resultgauge,resultgaugeb,WardP)
    
!      if (WardP) then

      call penlinemmm(p1,p23,p4,p5,p6,barpsi_p6,psi_p1,mup3,mup4,mup5,
     - alpha,musq,9,-1,resultgauge,result6,resultgaugeb,resultb6)


         endif


      if(gauge.eq.1) then

      test(1)=abs((result6-cte*resultb6)/(result(1)-cte*resultb(1))-1d0)

      t=test(1)
c NAN
      if(t+1.eq.t) then
          t=1d9
       endif  
      if(t.ne.t) then
          t=1d9
       endif  


      If(t.lt.GaugeAccuracyHexL) WardH=.True.
      if(printGaugePrec.and.t.gt.comGaugePrec) then
	print*,'GaugePrec',t
	comGaugePrec=t
      endif
	
      return

      elseif(gauge.eq.2) then
   
      test(2)=abs(((result4-result5)-cte*(resultb4-resultb5))/(result(2)-cte*resultb(2))-1d0)

      t=test(2)
c NAN
      if(t+1.eq.t) then
          t=1d9
       endif  
      if(t.ne.t) then
          t=1d9
       endif  

      If(t.lt.GaugeAccuracyHexL) WardH=.True.
      if(printGaugePrec.and.t.gt.comGaugePrec) then
	print*,'GaugePrec',t
	comGaugePrec=t
      endif
      return

      elseif(gauge.eq.3) then

      test(3)=abs(((result2-result3)-cte*(resultb2-resultb3))/(result(3)-cte*resultb(3))-1d0)
      t=test(3)
c NAN
      if(t+1.eq.t) then
          t=1d9
       endif  
      if(t.ne.t) then
          t=1d9
       endif  

      If(t.lt.GaugeAccuracyHexL) WardH=.True.
      if(printGaugePrec.and.t.gt.comGaugePrec) then
	print*,'GaugePrec',t
	comGaugePrec=t
      endif
      return

      elseif(gauge.eq.4) then


      test(4)=abs(((-result1)-cte*(-resultb1))/(result(4)-cte*resultb(4))-1d0)
      t=test(4)
c NAN
      if(t+1.eq.t) then
          t=1d9
       endif  
      if(t.ne.t) then
          t=1d9
       endif  

      If(t.lt.GaugeAccuracyHexL) WardH=.True.
      if(printGaugePrec.and.t.gt.comGaugePrec) then
	print*,'GaugePrec',t
	comGaugePrec=t
      endif
      return

      elseif(gauge.gt.4) then


      test(1)=abs((result6-cte*resultb6)/(result(1)-cte*resultb(1))-1d0)
  
      test(2)=abs(((result4-result5)-cte*(resultb4-resultb5))/(result(2)-cte*resultb(2))-1d0)
      test(3)=abs(((result2-result3)-cte*(resultb2-resultb3))/(result(3)-cte*resultb(3))-1d0)
      test(4)=abs(((-result1)-cte*(-resultb1))/(result(4)-cte*resultb(4))-1d0)

c NAN
      do i=1,4
         if(test(i) +1.eq.test(i)) then
            test(i) =1d9
          endif  
      if(test(i).ne.test(i))then
         test(i)=1d9
      endif

       enddo

      t=Max(test(1),test(2),test(3),test(4))

c NAN
      if(t+1.eq.t) then
          t=1d9
       endif  
      if(t.ne.t)then
         t=1d9
       endif
      !Print*, test(1),test(1)
      !Print*, test(2),test(2)
      !Print*, test(3),test(3)
      !Print*, test(4),test(4)
      
!      Print*, 't', t
      
      If(t.lt.GaugeAccuracyHexL) WardH=.True.
      if(printGaugePrec.and.t.gt.comGaugePrec) then
	print*,'GaugePrec',t
	comGaugePrec=t
      endif
     
c      Print*, 'Ward', Ward
      Return
      endif  ! If gauge
!      endif
!      endif
      Return
      end
c###################################################################
      subroutine hexGauge1gNoAbe(Accuracy,p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,
     -  mup2,mup3,mup4,mup6,alpha,musq,
     - resultd,WardH,gauge)
c 
c Author: Francisco Campanario
C Date: 16/01/2009
C Compute the gauge three gauge test for HexagonLine and result .True. if the max of it
c is larger that 1d-2. In the same subroutines the penlines are checked as well.
      IMPLICIT NONE
      Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3),p23(0:3),
     - p34(0:3)
      Complex*16 barpsi_p5(2),psi_p1(2),mup2(0:3),mup3(0:3),mup4(0:3),
     - mup6(0:3) 
      Complex*16 result3,result4,result5,result6

      Complex*16 resultgauge(3)
      Complex*16 resultptemp(8)
      Complex*16 resultd(4)

      real*8 pi,pi2,test(4),t1,t2
      parameter (pi =3.14159265358979323846264338327950288d0)
      parameter (pi2 =9.86960440108935861883449099987615114d0)
      external dotrr,theta
      real*8 dotrr,theta,musq

      Integer i,alpha,gauge
      Logical WardH,WardP,WardP1,WardP2
      real*8 Accuracy,Accuracymy1
      integer countb1,countb2,countp1,countp2,counth1,counth2     
      common/gaugetest/countb1,countb2,countp1,countp2,counth1,counth2
      real*8 ptd
      real*8 td
      real*8 tdmax,tdmin
      logical printGaugePrec
      parameter (printGaugePrec = .false.)
      real*8 comGaugePrec
      common /VVJJGaugePrec/ comGaugePrec
      common/ hexgaugeValue/t
      real*8 t
      t=1d9
CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Accuracymy1=Accuracy!*(1d-7)
CCCCCCCCCCCCCCCCCCCCCCC
      t1=1d8
      t2=1d8
      td=1d8
      tdmin=1d8
      tdmax=1d8
      ptd=1d8
      WardH=.False.
      WardP=.False.
      WardP1=.False.
      WardP2=.False.
CCCCCCCCC
      do i=0,3
      p23(i)=p2(i)+p3(i)
      p34(i)=p3(i)+p4(i)
      enddo

      if((gauge.ge.3).or.(gauge.eq.2))then

ccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc
c pen's 2 for contraction 2,3
ccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc      
      call penlineNoAbemmm(p1,p2,p34,p5,p6,barpsi_p5,psi_p1,mup2,mup3,mup6,
     - alpha,musq,4,1,resultgauge,resultpTEMP(1))
      RESULT3=resultptemp(1)
      call pentNoAbe_gauge(Accuracymy1,p1,p2,p34,p5,p6,barpsi_p5,psi_p1,
     -  mup2,mup3,mup6,alpha,musq,
     - resultgauge,
     - ptd,WardP1)
c       WardP=WardP1
c      if (WardP) then
      call penlineNoAbemmm(p1,p2,p34,p5,p6,barpsi_p5,psi_p1,mup2,mup4,mup6,
     -  alpha,musq,9,-1,resultgauge,resultptemp(1))
      result4=resultptemp(1)

      endif ! If gauge

c      endif
CCCCCCCC
c      wardp=.false.
CCCCCCCCCCCCCCCC
      wardp=.true.
CCCCCCCC
       if (.not.WardP) then
       countp2=countp2+1
       endif

      if (WardP) then
      ptd=1d8


      if((gauge.eq.1).or.(gauge.eq.2).or.(gauge.gt.3))then
ccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc
c pen's 1 for contraction 1,2
ccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc      
      call penlineNoAbemmm(p1,p23,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup4,mup6,
     - alpha,musq,4,1,resultgauge,resultptemp(1))
       result5=resultptemp(1)
      call pentNoAbe_gauge(Accuracymy1,p1,p23,p4,p5,p6,barpsi_p5,psi_p1,
     -  mup2,mup4,mup6,alpha,musq,
     - resultgauge,
     - ptd,WardP1)
      WardP=WardP1
c      if (WardP) then
      call penlineNoAbemmm(p1,p23,p4,p5,p6,barpsi_p5,psi_p1,mup3,mup4,mup6,
     - alpha,musq,9,-1,resultgauge,resultptemp(1))
       result6=resultptemp(1)

       endif ! If gauge

c       endif
CCCCCCCCCCCCCCCCCCCc
c      wardP=.false.
CCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCc
      wardP=.true.
CCCCCCCCCCCCCCCCCCC
       if (.not.WardP) then
       countp2=countp2+1
       endif
       endif
ccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc
cc Gauge test for Hexline cc 
ccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc
      if (WardP) then


      if(gauge.eq.1) then 

      test(1)=abs(-(result6)/(resultd(1))+1d0)
      td=test(1)

c NAN
      if(td+1.eq.td) then
          td=1d9
       endif 
      if(td.ne.td) then
          td=1d9
       endif 


      t=td

      If(td.lt.Accuracy) then
      WardH=.True.
      else
      counth1=counth1+1
      endif
   

      elseif(gauge.eq.2) then 
      test(2)=abs(-((result4-result5))/(resultd(2))+1d0)
      td=test(2)
c NAN
      if(td+1.eq.td) then
          td=1d9
       endif 
      if(td.ne.td) then
          td=1d9
       endif 


      t=td
      If(td.lt.Accuracy) then
      WardH=.True.
      else
      counth1=counth1+1
      endif
   
      elseif(gauge.eq.3) then 

      test(3)=abs(((result3)/(resultd(3))+1d0))
      td=test(3)
c NAN
      if(td+1.eq.td) then
          td=1d9
       endif 
      if(td.ne.td) then
          td=1d9
       endif 

      t=td
      If(td.lt.Accuracy) then
      WardH=.True.
      else
      counth1=counth1+1
      endif
   

      elseif(gauge.gt.3) then 

      test(1)=abs(-(result6)/(resultd(1))+1d0)
      test(2)=abs(-((result4-result5))/(resultd(2))+1d0)
      test(3)=abs(((result3)/(resultd(3))+1d0))
c NAN
      do i=1,3
         if(test(i) +1.eq.test(i)) then
            test(i) =1d9
          endif  
      if(test(i).ne.test(i))then
         test(i)=1d9
      endif

       enddo

C
        IF (test(1) .lt. test(2)) THEN           ! a < b here
            IF (test(1) .lt. test(3)) THEN       ! a < c     : a the smallest
                IF (test(2) .lt. test(3)) THEN   ! b < c  : a < b < c
                td=test(2)
                tdmin=test(1)
                tdmax=test(3)
                ELSE                             ! c <= b : a < c <= b
                td=test(3)
                tdmin=test(1)
                tdmax=test(2)
                END IF
            ELSE                                 ! a >= c    : c <= a < b
            td=test(1)
            tdmin=test(3)
            tdmax=test(2)
            END IF
       ELSE                                      ! b <= a here
            IF (test(2) .lt. test(3)) THEN       ! b < c     : b the smallest
                IF (test(1) .lt. test(3)) THEN   ! a < c   : b <= a < c
                td=test(1)
                tdmin=test(2)
                tdmax=test(3)
                ELSE                             ! a >= c  : b < c <= a
                td=test(3)
                tdmin=test(2)
                tdmax=test(1)
                END IF
             ELSE                                ! c <= b    : c <= b <= a
             td=test(2)
             tdmin=test(3)
             tdmax=test(1)
             END IF
       END IF


       td=(8d0*td+tdmin+tdmax)/(10d0)

c NAN
      if(td+1.eq.td) then
          td=1d9
       endif 
      if(td.ne.td) then
          td=1d9
       endif 


        t=td
      If(td.lt.Accuracy) then
      WardH=.True.
      else
      counth1=counth1+1
      endif
   
      endif ! If gauge
      if(printGaugePrec.and.td.gt.comGaugePrec) then
	print*,'GaugePrec',td
	comGaugePrec=td
      endif
      
      Return
      endif
    
 
      if (.not.WardP) then
      counth1=counth1+1
      counth2=counth2+1
      endif


      Return
      end
c################################################################################


      subroutine box_gaugeMM(Accuracy,p1,p2,p3,p4,musq,
     -  result,resultb,td,WardB)
c
c Author: Francisco Campanario
C Date: 16/01/2009
C Compute the gauge  for the boxes
c When the ward test is not better than Accuraccy then WardB=.false.
c result: result(2) ==resultgauge(1) and resultgauge(2) out of boxline 
cvwith gaugetest=3
C resultb:resultb(2) ==resultgaugeb(1) and resultgaugeb(2) out of
c boxline with gaugetest=3
c p1,...p4, momenta of the boxline
C musq renormalization invariance. Same as used in boxline
c td: result of the gaugetest. For normal points returns ~1d-13 it is
c what the determines how good the point is and drives the result of
c WardB
c Accuracy: Accuracy demanded for the gauge test. With 1d-3 only 1/1d6 
c should be instable.

      IMPLICIT NONE
      Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      Real*8 musq,test1,test2
      Complex*16 cte

      Complex*16 result(2)
      Complex*16 resultb(2)
      real*8 pi,pi2,td
      parameter (pi =3.14159265358979323846264338327950288d0)
      parameter (pi2 =9.86960440108935861883449099987615114d0)
      external dotrr,theta
      real*8 dotrr,s14,lm,ls,ths,theta
      Logical WardB
      double precision Accuracy

      WardB=.False.
      td=1.d08


c        ctedieter = (3*pi*(0d0,1d0)-pi**2)
c        cte=-(3*pi*(0d0,1d0)-4.0d0*pi**2/3.0d0+7.0d0)

      s14=2d0*dotrr(p1,p4)+dotrr(p1,p1)+dotrr(p4,p4)
      lm=Log(musq)
      ls=Log(Abs(s14))
      ths=theta(s14)
      
c       cte = -cte_dieter but general dependency in musq.
c         The sign is later on corrected in his code in
c          -dotcr(j,p)
      cte= -3d0*(Lm - Ls + (0d0,1d0)*pi*ths)
     $     - (Lm - Ls + (0d0,1d0)*pi*ths)**2
c     substracted the c_virt
     $     -7d0 + 1d0/3d0*pi2
c variables for gauge test:test1,test2

      if(abs(resultb(1)).lt.1d-13) then
      test1=abs(result(1)- resultb(1)*(cte))
      else
      test1=abs(result(1)/(resultb(1)*(cte))-1d0)
      endif

      if(abs(resultb(2)).lt.1d-13) then
       test2=abs(result(2)- (resultb(2)*(cte)))
      else
      test2=abs(result(2)/(resultb(2)*(cte))-1d0)
      endif
      td=max(test1,test2)


      IF (Td.LT.Accuracy) then
      WardB=.True.
      ENDIF


      Return
      end

c###################################################################
      subroutine pent_gaugeT(p1,p2,p3,p4,p5,barpsi_p5,psi_p1,
     -     mup2,mup3,mup4,alpha,musq,result,resultb,Ward,Div)
c     
c     Author: Francisco Campanario
C     Date: 16/01/2009
C     Compute the gauge three gauge test for the pentagons and result .True. if the max of it
c     is larger that 1d-2. In the same subroutines the boxlines called are checked as well.


      IMPLICIT NONE
      DOUBLE PRECISION   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p23(0:3),p34(0:3)
      DOUBLE COMPLEX barpsi_p5(2),psi_p1(2),mup2(0:3),mup3(0:3),mup4(0:3) 

      DOUBLE PRECISION musq,test1,test2
      DOUBLE COMPLEX cte,cte1
      Integer alpha,i
      DOUBLE COMPLEX result(3),result1(5),result2(5),result3(5),result4(5)
      DOUBLE COMPLEX resultb(3),resultb1,resultb2,resultb3,resultb4
      DOUBLE COMPLEX resultgauge(2),resultgaugeb(2)
      DOUBLE PRECISION pi,pi2,test(3),t
      parameter (pi =3.14159265358979323846264338327950288d0)
      parameter (pi2 =9.86960440108935861883449099987615114d0)
      external dotrr,theta
      DOUBLE PRECISION dotrr,s15,lm,ls,ths,theta
      Logical Ward
      real*8 ratiocheckCan
      integer Div
#include "GaugeTest.inc"
      Real*8 gaugetestacc
c      parameter(gaugetestacc=GaugeAccuracyPenL)

      gaugetestacc=GaugeAccuracyPenL
      Ward=.False.

c      print*,"gaugetest accuracy",gaugetestacc

      do i=0,3
c     emup2(i)=mup2(i)-x2*(-p2(i))
c     emup3(i)=mup3(i)-x3*(-p3(i))
c     emup4(i)=mup4(i)-x4*(-p4(i))
         p34(i)=p3(i)+p4(i)
         p23(i)=p2(i)+p3(i)
      enddo

c     print*, 'zepp', emup2
c     print*, 'xze', x2
c     print*, 'pe', p2(0)
c     print*, 'zmupp' ,emup3
c     print*, 'ataup', emup4
      

      if(Div.eq.0) then
c     ctedieter = (3*pi*(0d0,1d0)-pi**2)
c     cte=-(3*pi*(0d0,1d0)-4.0d0*pi**2/3.0d0+7.0d0)
         s15=2d0*dotrr(p1,p5)+dotrr(p1,p1)+dotrr(p5,p5)
         lm=Log(musq)
         ls=Log(Abs(s15))
         ths=theta(s15)
c     cte = -cte_dieter but general dependency in musq.
c     The sign is later on corrected in his code in
c     -dotcr(j,p)
         cte1= -3d0*(Lm - Ls + (0d0,1d0)*pi*ths) 
     -        - (Lm - Ls + (0d0,1d0)*pi*ths)**2
c         pause
         cte=cte1
c     substracted the c_virt
     -        -7d0 + 1d0/3d0*pi2
      elseif(Div.eq.2) then
         cte=-2d0
      elseif(Div.eq.1) then
         s15=2d0*dotrr(p1,p5)+dotrr(p1,p1)+dotrr(p5,p5)
         ls=log(Abs(s15))-(0,1)*pi*theta(s15)-log(musq)
c     print*, 'Not programmed yet'
c     stop
         cte=(-3d0+2d0*ls)
      endif

c     box's 2 for contraction 3, 2
c$$$      call boxlinemm(p1,p2,p34,p5,barpsi_p5,psi_p1,mup2,mup3
c$$$     -     ,alpha,musq,3,1,resultgauge,result1,resultgaugeb,resultb1)
      call boxlineABETotal(p1,p2,p34,p5,barpsi_p5,psi_p1,mup2,mup3
     -     ,alpha,musq,0,1,3,1,resultgauge,result1(1),resultgaugeb,resultb1,Div)
c     variables for gauge test:test1,test2      
      

      if(abs(resultgaugeb(1)).lt.1d-4) then
         test1=abs(resultgauge(1)- resultgaugeb(1)*(cte))
      else
         test1=abs(resultgauge(1)/(resultgaugeb(1)*(cte))-1d0)
      endif

      if(abs(resultgaugeb(2)).lt.1d-4) then
         test2=abs(resultgauge(2)- (resultgaugeb(2)*(cte)))
      else   
         test2=abs(resultgauge(2)/(resultgaugeb(2)*(cte))-1d0)
      endif

c$$$      print*,'test1 pen',test1,resultgauge(1),resultgaugeb(1)
c$$$      print*,'test2 pen',test2,resultgauge(2),resultgaugeb(2)

c     test1=abs(resultgauge(1)/(resultgaugeb(1)*(cte))-1d0)
c     test2=abs(resultgauge(2)/(resultgaugeb(2)*(cte))-1d0)


!     if (test1.lt.1d-2) then
!     if (test2.lt.1d-2) then
      
c$$$      call boxlinemm(p1,p2,p34,p5,barpsi_p5,psi_p1,mup2,mup4
c$$$     -     ,alpha,musq,9,-1,resultgauge,result2,resultgauge,resultb2)
      call boxlineABETotal(p1,p2,p34,p5,barpsi_p5,psi_p1,mup2,mup4
     -     ,alpha,musq,0,1,9,-1,resultgauge,result2,resultgauge,resultb2,Div)

c     box's 1 for contraction 2,1
      
c$$$      call boxlinemm(p1,p23,p4,p5,barpsi_p5,psi_p1,mup2,mup4
c$$$     -     ,alpha,musq,3,1,resultgauge,result3,resultgaugeb,resultb3)
      call boxlineABETotal(p1,p23,p4,p5,barpsi_p5,psi_p1,mup2,mup4
     -     ,alpha,musq,0,1,3,1,resultgauge,result3,resultgaugeb,resultb3,Div)
c     gauge test: test1,test2
      
      if(abs(resultgaugeb(1)).lt.1d-4) then
         test1=abs(resultgauge(1)- resultgaugeb(1)*(cte))
      else
         test1=abs(resultgauge(1)/(resultgaugeb(1)*(cte))-1d0)
      endif

      if(abs(resultgaugeb(2)).lt.1d-4) then
         test2=abs(resultgauge(2)- (resultgaugeb(2)*(cte)))
      else   
         test2=abs(resultgauge(2)/(resultgaugeb(2)*(cte))-1d0)
      endif

c     test1=abs(resultgauge(1)/(resultgaugeb(1)*(cte))-1d0)
c     test2=abs(resultgauge(2)/(resultgaugeb(2)*(cte))-1d0)
c$$$      print*,'test1 pen',test1,resultgauge(1),resultgaugeb(1)
c$$$      print*,'test2 pen',test2,resultgauge(2),resultgaugeb(2)
c$$$      pause
      
      
!     if (test1.lt.1d-2) then
!     if  (test2.lt.1d-2) then

c$$$      call boxlinemm(p1,p23,p4,p5,barpsi_p5,psi_p1,mup3,mup4
c$$$     -     ,alpha,musq,9,-1,resultgauge,result4,resultgaugeb,resultb4)

      call boxlineABETotal(p1,p23,p4,p5,barpsi_p5,psi_p1,mup3,mup4
     -     ,alpha,musq,0,1,9,-1,resultgauge,result4,resultgaugeb,resultb4,Div)

      ratiocheckCan=(Abs(-(result4(1))/(-cte*(resultb4)))-1d0)


      if(ratiocheckCan.gt.1d-2) then
         test(1)=abs((result4(1)-cte*resultb4)/(result(1)-cte*resultb(1))+1d0)
      else
         test(1)=abs((result4(1))/(result(1))+1d0)
      endif


      ratiocheckCan=Abs(-result(2)/(-cte*resultb(2)))-1d0


      if(ratiocheckCan.gt.1d-2) then
         test(2)=abs(((result2(1)-result3(1))-cte*(resultb2-resultb3))/(result(2)-cte*resultb(2))+1d0)
         
      else
         test(2)=abs((result2(1)-result3(1))/result(2)+1d0)
      endif


      ratiocheckCan=(Abs(-(-result1(1))/(-cte*(-resultb1)))-1d0)
      

      if(ratiocheckCan.gt.1d-2) then
         
         test(3)=abs(((-result1(1))-cte*(-resultb1))/(result(3)-cte*resultb(3))+1d0)
      else
         test(3)=abs(-result1(1)/result(3)+1d0)
      endif


c Checking for NAN
      If(test(1)+1.eq.test(1))then
        test(1)=1d8
      endif
c Checking for NAN
      If(test(2)+1.eq.test(2))then
        test(2)=1d8
      endif
c Checking for NAN
      If(test(3)+1.eq.test(3))then
        test(3)=1d8
      endif

c Checking for NAN
      If(test(1).ne.test(1))then
        test(1)=1d8
      endif
c Checking for NAN
      If(test(2).ne.test(2))then
        test(2)=1d8
      endif
c Checking for NAN
      If(test(3).ne.test(3))then
        test(3)=1d8
      endif

      
      t=Max(test(1),test(2),test(3))
      If(t.lt.gaugetestacc) then
         Ward=.True.
      else
c         print*, "cte1",cte1
c         print*, "scale",musq
c         print*, "s15", s15
c         print*, "test(1)",test(1),(result4(1)-cte*resultb4),(result(1)-cte*resultb(1))
c         print*, "test(2)",test(2),(result2(1)-result3(1)-cte*(resultb2-resultb3)),(result(2)-cte*resultb(2))
c         print*, "test(3)",test(3),result1(1),result(3)
c         cte=0d0
c         print*, "test(1)",test(1),(result4(1)-cte*resultb4),(result(1)-cte*resultb(1))
c         print*, "test(2)",test(2),(result2(1)-result3(1)),(result(2)-cte*resultb(2))
c         print*, "test(3)",test(3),result1(1),result(3)
       endif


      Return
!     endif
!     endif
!     endif
!     endif
      Return
      end
c#############

c###################################################################
      subroutine box_gaugeT(p1,p2,p3,p4,musq,result,resultb,ward,Div)
C Compute the two gauge tests for the boxes and result. .True. if the max of it
c is smaller than 1d-2. 

      IMPLICIT NONE
      DOUBLE PRECISION p1(0:3), p2(0:3), p3(0:3), p4(0:3) 
      DOUBLE COMPLEX result(2), resultb(2)
      DOUBLE PRECISION musq, testb(2)
      DOUBLE COMPLEX cte
      DOUBLE PRECISION dotrr,pi,pi2
      EXTERNAL dotrr
      PARAMETER (pi =3.14159265358979323846264338327950288d0)
      PARAMETER (pi2 =9.86960440108935861883449099987615114d0)
      DOUBLE PRECISION s14,lm,ls,ths,theta,gaugetestacc
      LOGICAL ward
      Integer Div
#include "GaugeTest.inc"

      gaugetestacc=GaugeAccuracyBoxL

      ward=.False.

      if(Div.eq.0) then
      s14=2d0*dotrr(p1,p4)+dotrr(p1,p1)+dotrr(p4,p4)
      lm=Log(musq)
      ls=Log(Abs(s14))
      ths=theta(s14)


c     cte = -cte_dieter but general dependency in musq.
c       The sign is later on corrected in his code in
c        -dotcr(j,p)
      cte= -3d0*(Lm - Ls + (0d0,1d0)*pi*ths) 
     -     - (Lm - Ls + (0d0,1d0)*pi*ths)**2
c   substracted the c_virt
     -     -7d0 + 1d0/3d0*pi2

      elseif(Div.eq.2) then
         cte=-2d0
      elseif(Div.eq.1) then
        s14=2d0*dotrr(p1,p4)+dotrr(p1,p1)+dotrr(p4,p4)
	ls=log(Abs(s14))-(0,1)*pi*theta(s14)-log(musq)
c         print*, 'Not program yet'
c         stop
         cte=(-3d0+2d0*ls)
       endif


      IF(ABS(RESULTB(1)).LT.1D-4) THEN
      TESTB(1)=ABS(RESULT(1)-CTE*RESULTB(1))
      ELSE
      testb(1)=abs(result(1)/(cte*resultb(1))-1d0)
      ENDIF
      
      if(ABS(RESULTB(2)).LT.1D-4) THEN
      TESTB(2)=ABS(RESULT(2)-CTE*RESULTB(2))
      ELSE
      testb(2)=abs(result(2)/(cte*resultb(2))-1d0)
      ENDIF


c Checking for NAN
      If(testb(1)+1.eq.testb(1))then
        testb(1)=1d8
      endif
c Checking for NAN
      If(testb(2)+1.eq.testb(2))then
        testb(2)=1d8
      endif
      
c Checking for NAN
      If(testb(1).ne.testb(1))then
        testb(1)=1d8
      endif
c Checking for NAN
      If(testb(2).ne.testb(2))then
        testb(2)=1d8
      endif


      if ( (testb(1).lt.gaugetestacc) .and. (testb(2).lt.gaugetestacc)) then
      ward = .true.
      else
cfc      print*, "GaugeAccuracyBoxL",GaugeAccuracyBoxL
cfc      print*, "gaugetestacc",gaugetestacc
cfc      print*,'MAXtestb(1)', MAX(testb(1),testb(2))
cfc      print*,'testb(1)', testb(1)
cfc      print*,'testb(2)', testb(2)
cfc      print*, "cte",cte
cfc      print*, 'resultb(1)',resultb(1)
cfc      print*, 'resultb(2)',resultb(2)  
cfc      print*, 'result(1)',result(1)      
cfc      print*, 'result(2)',result(2)      
c      pause
      endif

c      Print*, 'Ward', Ward

      return
      end


********************************************************************************
