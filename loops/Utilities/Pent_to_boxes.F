      subroutine pent_to_box(p1,p2,p3,p4,p5,barpsi_p5,psi_p1,
     -  mup2,mup3,mup4,x2,x3,x4,alpha,musq,result,resultb)
c 
c Author: Francisco Campanario
C Date: 16/01/2009
C Compute the box contribution when mupj->emupj +xj*pj
C That is the box contribution after the shift of the polarization vectors
C Implement the gauge test here for boxes
C This subroutine is 8 times faster than the equivalent 
c  VVV_J_virtual_tri_box_pent45 with option nlo=5.
c this is due to a factor 4 in boxline itself.And the number 
c of callings done in this subroutine. We can use boxline with option comp=-1
c so that the loop internal dependency is computed only once for each of the
c momentum configurations that appear.
c Note: it seems to be more precise as well: 
c       (True pentagons)+ (reult here)=pentagon with 15 digits of precision.
c If one uses VVV_J_virtual_tri_box_pent45 the precision
c is of 5-6 digits for a good point! Frankly speaking, I do not understand
c the lost of precision in the other subroutine.

      IMPLICIT NONE
      DOUBLE PRECISION p1(0:3), p2(0:3), p3(0:3), p4(0:3) 
      DOUBLE PRECISION p5(0:3), p23(0:3), p34(0:3)
      DOUBLE COMPLEX barpsi_p5(2),psi_p1(2),mup2(0:3),mup3(0:3),mup4(0:3) 
      DOUBLE COMPLEX emup3(0:3),emup4(0:3) 
      DOUBLE PRECISION musq,test1,test2
      DOUBLE COMPLEX x2,x3,x4,cte
      Integer alpha,i
      DOUBLE COMPLEX result,result1,result2,result3,result4
      DOUBLE COMPLEX resultb,resultb1,resultb2,resultb3,resultb4
      DOUBLE COMPLEX resultgauge(2),resultgaugeb(2)
      DOUBLE PRECISION pi,pi2
      parameter (pi =3.14159265358979323846264338327950288d0)
      parameter (pi2 =9.86960440108935861883449099987615114d0)
      external dotrr,theta
      DOUBLE PRECISION dotrr,s15,lm,ls,ths,theta


      do i=0,3
c      emup2(i)=mup2(i)-x2*(-p2(i))
      emup3(i)=mup3(i)-x3*(-p3(i))
      emup4(i)=mup4(i)-x4*(-p4(i))
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
      if(abs(resultgaugeb(1)).lt.1d-6) then
      test1=abs(resultgauge(1)- resultgaugeb(1)*(cte))
      else
      test1=abs(resultgauge(1)/(resultgaugeb(1)*(cte))-1d0)
      endif

      if(abs(resultgaugeb(2)).lt.1d-6) then
       test2=abs(resultgauge(2)- (resultgaugeb(2)*(cte)))
      else   
      test2=abs(resultgauge(2)/(resultgaugeb(2)*(cte))-1d0)
      endif

      if (test1.lt.1d-2) then
      if (test2.lt.1d-2) then
      
      call boxlinemm(p1,p2,p34,p5,barpsi_p5,psi_p1,mup2,emup4
     - ,alpha,musq,9,-1,resultgauge,result2,resultgauge,resultb2)

c box's 1 for contraction 2,1
      
      call boxlinemm(p1,p23,p4,p5,barpsi_p5,psi_p1,mup2,emup4
     - ,alpha,musq,3,1,resultgauge,result3,resultgaugeb,resultb3)

c gauge test: test1,test2

      if(abs(resultgaugeb(1)).lt.1d-6) then
      test1=abs(resultgauge(1)- resultgaugeb(1)*(cte))
      else
      test1=abs(resultgauge(1)/(resultgaugeb(1)*(cte))-1d0)
      endif

      if(abs(resultgaugeb(2)).lt.1d-6) then
       test2=abs(resultgauge(2)- (resultgaugeb(2)*(cte)))
      else   
      test2=abs(resultgauge(2)/(resultgaugeb(2)*(cte))-1d0)
      endif

c      test1=abs(resultgauge(1)/(resultgaugeb(1)*(cte))-1d0)
c      test2=abs(resultgauge(2)/(resultgaugeb(2)*(cte))-1d0)
    
      if (test1.lt.1d-2) then
      if  (test2.lt.1d-2) then

      call boxlinemm(p1,p23,p4,p5,barpsi_p5,psi_p1,emup3,emup4
     - ,alpha,musq,9,-1,resultgauge,result4,resultgaugeb,resultb4)

      result=x4*(-result1)+x3*(result2-result3)+x2*result4
      resultb=x4*(-resultb1)+x3*(resultb2-resultb3)+x2*resultb4

      else
      result=0d0
      resultb=0d0
      endif
      else
      result=0d0
      resultb=0d0
      endif
      else
      result=0d0
      resultb=0d0
      endif
      else
      result=0d0
      resultb=0d0
      endif

      Return
      end
c###################################################################
      subroutine pent_to_box_2hel(p1,p2,p3,p4,p5,barpsi_p5,psi_p1,
     -  mup2,mup3,mup4,x2,x3,x4,musq,result,resultb)
c 
c Author: Francisco Campanario
C Date: 16/01/2009
C Compute the box contribution when mupj->emupj +xj*pj
C That is the box contribution after the shift of the polarization vectors
C Implement the gauge test here for boxes
C This subroutine is 16 times faster than the equivalent 
c  VVV_J_virtual_tri_box_pent45 with option nlo=5.
c this is due to a factor 4 in boxline itself.And the number 
c of callings done in this subroutine. We can use boxline with option comp=-1
c so that the loop internal dependency is compute it only once for each of the
c momenta configuration that appear.
c Note: it seems to be more precise as well: (True pentagons)+ (reult here)=pentagon
c with 15 digits of precision. If one uses VVV_J_virtual_tri_box_pent45 the precision
c is of 5-6 digits for a good point! Frankly speaking, I do not understand the lost
c of precision in the other subroutine.

      IMPLICIT NONE
      DOUBLE PRECISION   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p23(0:3),p34(0:3)
      DOUBLE COMPLEX barpsi_p5(2,-1:1),psi_p1(2,-1:1),mup2(0:3),mup3(0:3),mup4(0:3) 
      DOUBLE COMPLEX emup3(0:3),emup4(0:3) 
      DOUBLE PRECISION musq,test1,test2
      DOUBLE COMPLEX x2,x3,x4,cte
      Integer i
      DOUBLE COMPLEX result(-1:1),result1,result2,result3,result4
      DOUBLE COMPLEX result5,result6,result7,result8
      DOUBLE COMPLEX resultb5,resultb6,resultb7,resultb8
      DOUBLE COMPLEX resultb(-1:1),resultb1,resultb2,resultb3,resultb4
      DOUBLE COMPLEX resultgauge(2),resultgaugeb(2)
      DOUBLE PRECISION pi,pi2
      parameter (pi =3.14159265358979323846264338327950288d0)
      parameter (pi2 =9.86960440108935861883449099987615114d0)
      external dotrr,theta
      DOUBLE PRECISION dotrr,s15,lm,ls,ths,theta


      do i=0,3
c      emup2(i)=mup2(i)-x2*(-p2(i))
      emup3(i)=mup3(i)-x3*(-p3(i))
      emup4(i)=mup4(i)-x4*(-p4(i))
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
c$$$        print*, 'bp5',barpsi_p5
c$$$        print*, 'bp5', barpsi_p5(1,-1)
c$$$        print*, 'bp5', barpsi_p5(2,-1)
c$$$        print*, 'p1', psi_p1
c$$$        print*, 'bp1', psi_p1(1,-1)
c$$$        print*, 'bp1', psi_p1(2,-1)

c box's 2 for contraction 3, 2
      call boxlinemm(p1,p2,p34,p5,barpsi_p5(1,-1),psi_p1(1,-1),mup2,mup3
     - ,-1,musq,3,1,resultgauge,result1,resultgaugeb,resultb1)

c$$$      print*, 'result1', result1
c$$$      print*,'resultb1', resultb1
c variables for gauge test:test1,test2      


      if(abs(resultgaugeb(1)).lt.1d-6) then
      test1=abs(resultgauge(1)- resultgaugeb(1)*(cte))
      else
      test1=abs(resultgauge(1)/(resultgaugeb(1)*(cte))-1d0)
      endif

      if(abs(resultgaugeb(2)).lt.1d-6) then
       test2=abs(resultgauge(2)- (resultgaugeb(2)*(cte)))
      else   
      test2=abs(resultgauge(2)/(resultgaugeb(2)*(cte))-1d0)
      endif

c      test1=abs(resultgauge(1)/(resultgaugeb(1)*(cte))-1d0)
c      test2=abs(resultgauge(2)/(resultgaugeb(2)*(cte))-1d0)

c$$$      PRINT*, 'TEST1', TEST1
c$$$      PRINT*,  'test2', test2
      if (test1.lt.1d-2) then
      if (test2.lt.1d-2) then
      
      call boxlinemm(p1,p2,p34,p5,barpsi_p5(1,1),psi_p1(1,1),mup2,mup3
     - ,1,musq,9,-1,resultgauge,result5,resultgaugeb,resultb5)

      call boxlinemm(p1,p2,p34,p5,barpsi_p5(1,-1),psi_p1(1,-1),mup2,emup4
     - ,-1,musq,9,-1,resultgauge,result2,resultgauge,resultb2)

      call boxlinemm(p1,p2,p34,p5,barpsi_p5(1,1),psi_p1(1,1),mup2,emup4
     - ,1,musq,9,-1,resultgauge,result6,resultgauge,resultb6)

c box's 1 for contraction 2,1
      
      call boxlinemm(p1,p23,p4,p5,barpsi_p5(1,-1),psi_p1(1,-1),mup2,emup4
     - ,-1,musq,3,1,resultgauge,result3,resultgaugeb,resultb3)

c gauge test: test1,test2

      if(abs(resultgaugeb(1)).lt.1d-6) then
      test1=abs(resultgauge(1)- resultgaugeb(1)*(cte))
      else
      test1=abs(resultgauge(1)/(resultgaugeb(1)*(cte))-1d0)
      endif

      if(abs(resultgaugeb(2)).lt.1d-6) then
       test2=abs(resultgauge(2)- (resultgaugeb(2)*(cte)))
      else   
      test2=abs(resultgauge(2)/(resultgaugeb(2)*(cte))-1d0)
      endif

c      test1=abs(resultgauge(1)/(resultgaugeb(1)*(cte))-1d0)
c      test2=abs(resultgauge(2)/(resultgaugeb(2)*(cte))-1d0)
    
      if (test1.lt.1d-2) then
      if  (test2.lt.1d-2) then


      call boxlinemm(p1,p23,p4,p5,barpsi_p5(1,1),psi_p1(1,1),mup2,emup4
     - ,1,musq,9,-1,resultgauge,result7,resultgaugeb,resultb7)

      call boxlinemm(p1,p23,p4,p5,barpsi_p5(1,-1),psi_p1(1,-1),emup3,emup4
     - ,-1,musq,9,-1,resultgauge,result4,resultgaugeb,resultb4)

      call boxlinemm(p1,p23,p4,p5,barpsi_p5(1,1),psi_p1(1,1),emup3,emup4
     - ,1,musq,9,-1,resultgauge,result8,resultgaugeb,resultb8)


      result(-1)=x4*(-result1)+x3*(result2-result3)+x2*result4
      resultb(-1)=x4*(-resultb1)+x3*(resultb2-resultb3)+x2*resultb4

      result(1)=x4*(-result5)+x3*(result6-result7)+x2*result8
      resultb(1)=x4*(-resultb5)+x3*(resultb6-resultb7)+x2*resultb8

      else
      result(-1)=0d0
      resultb(-1)=0d0
      result(1)=0d0
      resultb(1)=0d0
      endif
      else
      result(-1)=0d0
      resultb(-1)=0d0
      result(1)=0d0
      resultb(1)=0d0
      endif
      else
      result(-1)=0d0
      resultb(-1)=0d0
      result(1)=0d0
      resultb(1)=0d0
      endif
      else
      result(-1)=0d0
      resultb(-1)=0d0
      result(1)=0d0
      resultb(1)=0d0
      endif
      Return
      end
