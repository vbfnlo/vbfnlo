
      subroutine pent_to_boxnoab(p1,p2,p3,p4,p5,barpsi_p5,psi_p1,
     -  mup2,mup3,mup4,x2,x3,x4,alpha,musq,glupos,result,resultb)
c 
c Author: Francisco Campanario, 
c modified to non abelian case by Christoph Englert
C Date: 16/01/2009
C Compute the box contribution when mupj->emupj +xj*pj
C That is the box contribution after the shift of the polarization vectors
C Implement the gauge test here for boxes
C This subroutine is 8 times faster than the equivalent 
c  VVV_J_virtual_tri_box_pent45 with option nlo=5.
c this is due to a factor 4 in boxline itself.And the number 
c of callings done in this subroutine. We can use boxline with option comp=-1
c so that the loop internal dependency is compute it only once for each of the
c momenta configuration that appear.
c Note: it seems to be more precise as well: 
c       (True pentagons)+ (reult here)=pentagon
c with 15 digits of precision. If one uses VVV_J_virtual_tri_box_pent45 the 
c precision is of 5-6 digits for a good point! Frankly speaking, I do not 
c understand the loss of precision in the other subroutine.

      IMPLICIT NONE
      Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p23(0:3),p34(0:3)
      Complex*16 barpsi_p5(2),psi_p1(2),mup2(0:3),mup3(0:3),mup4(0:3) 
      Complex*16 emup3(0:3),emup4(0:3) 
      Real*8 musq,test1,test2
      Complex*16 x2,x3,x4,cte
      Integer alpha,i,glupos
      Complex*16 result,resulttmp(5,4),
     &           result1noab,result2noab,result3noab,result4noab
      Complex*16 resultb,resultb1,resultb2,resultb3,resultb4
      Complex*16 resultgauge(3),resultgaugeb(3),cfc,cfcac
      real*8 pi,pi2
      parameter (pi =3.14159265358979323846264338327950288d0)
      parameter (pi2 =9.86960440108935861883449099987615114d0)
      external dotrr,theta
      real*8 dotrr,s15,lm,ls,ths,theta
      integer gluposar(4)
      logical ward

      ward=.true.

      if (glupos.eq.2) then
        gluposar(1)=2
        gluposar(2)=2
        gluposar(3)=2
        gluposar(4)=2
      elseif (glupos.eq.3) then
        gluposar(1)=3
        gluposar(2)=3
        gluposar(3)=3
        gluposar(4)=3
      elseif (glupos.eq.6) then
        gluposar(1)=3
        gluposar(2)=3 !not needed
        gluposar(3)=3 !not needed
        gluposar(4)=2 
      endif


      cfc=dcmplx(4d0/3d0,0d0)
      cfcac=dcmplx(-1d0/6d0,0d0)

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
!       call boxlinemm(p1,p2,p34,p5,barpsi_p5,psi_p1,mup2,mup3
!      - ,alpha,musq,3,1,resultgauge,result1,resultgaugeb,resultb1)

      CALL boxlineABEmm(p1,p2,p34,p5,barpsi_p5,psi_p1,mup2,mup3,
     &      alpha, musq, 1, gluposar(1) ,3,1,resultgauge,resulttmp(1,1),
     &      resultgaugeb,resultb1)
      result1noab=cfc*resulttmp(1,1)+cfcac*resulttmp(2,1)

c variables for gauge test:test1,test2
      if (abs(resultgauge(1)).ge.1e-10) then      
      test1=abs(resultgauge(1)/(resultgaugeb(1)*(cte))-1d0)
      else
      test1=abs(resultgauge(1)-(resultgaugeb(1)*(cte)))
      endif
      if (abs(resultgauge(2)).ge.1e-10) then      
      test2=abs(resultgauge(2)/(resultgaugeb(2)*(cte))-1d0)
      else
      test2=abs(resultgauge(2)-(resultgaugeb(2)*(cte)))
      endif

      if (test1.lt.1d-2) ward=.false.
      if (test2.lt.1d-2) ward=.false.
      
!       call boxlinemm(p1,p2,p34,p5,barpsi_p5,psi_p1,mup2,emup4
!      - ,alpha,musq,9,-1,resultgauge,result2,resultgauge,resultb2)

      call boxlineABEmm(p1,p2,p34,p5,barpsi_p5,psi_p1,mup2,emup4,
     &      alpha, musq, 1, gluposar(2) ,3,1,resultgauge,resulttmp(1,2),
     &      resultgaugeb,resultb2)
      result2noab=cfc*resulttmp(1,2)+cfcac*resulttmp(2,2)

c box's 1 for contraction 2,1
      
!       call boxlinemm(p1,p23,p4,p5,barpsi_p5,psi_p1,mup2,emup4
!      - ,alpha,musq,3,1,resultgauge,result3,resultgaugeb,resultb3)

      call boxlineABEmm(p1,p23,p4,p5,barpsi_p5,psi_p1,mup2,emup4,
     &      alpha, musq, 1, gluposar(3) ,3,1,resultgauge,resulttmp(1,3),
     &      resultgaugeb,resultb3)

      result3noab=cfc*resulttmp(1,3)+cfcac*resulttmp(2,3)

c gauge test: test1,test2
      if (abs(resultgauge(1)).ge.1e-10) then      
      test1=abs(resultgauge(1)/(resultgaugeb(1)*(cte))-1d0)
      else
      test1=abs(resultgauge(1)-(resultgaugeb(1)*(cte)))
      endif
      if (abs(resultgauge(2)).ge.1e-10) then      
      test2=abs(resultgauge(2)/(resultgaugeb(2)*(cte))-1d0)
      else
      test2=abs(resultgauge(2)-(resultgaugeb(2)*(cte)))
      endif
    
      if (test1.lt.1d-2) ward=.false.
      if  (test2.lt.1d-2) ward=.false.

!       call boxlinemm(p1,p23,p4,p5,barpsi_p5,psi_p1,emup3,emup4
!      - ,alpha,musq,9,-1,resultgauge,result4,resultgaugeb,resultb4)

      call boxlineABEmm(p1,p23,p4,p5,barpsi_p5,psi_p1,emup3,emup4,
     &      alpha, musq, 1, gluposar(4) ,3,1,resultgauge,resulttmp(1,4),
     &      resultgaugeb,resultb4)

      result4noab=cfc*resulttmp(1,4)+cfcac*resulttmp(2,4)

!       if (ward) then
        result=x4*(-result1noab)+x3*(result2noab-result3noab)+x2*result4noab
!       else
!         result=dcmplx(0d0,0d0)
!       endif
      resultb=x4*(-resultb1)+x3*(resultb2-resultb3)+x2*resultb4

!       print*, "here"
!       print*, result1/(resulttmp(1,1)+resulttmp(2,1))
!       print*, result2/(resulttmp(1,2)+resulttmp(2,2))
!       print*, result3/(resulttmp(1,3)+resulttmp(2,3))
!       print*, result4/(resulttmp(1,4)+resulttmp(2,4))
!       STOP


!       else
!       result=0d0
!       resultb=0d0
!       endif
!       else
!       result=0d0
!       resultb=0d0
!       endif
!       else
!       result=0d0
!       resultb=0d0
!       endif
!       else
!       result=0d0
!       resultb=0d0
!       endif

      Return
      end


********************************************************************************

