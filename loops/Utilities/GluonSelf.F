c ************************************************************************************
c Author: Matthias Kerner
c Date: 2012/13
c ************************************************************************************

c gluon self energy: gluon/ghost loops, n fermion loops, including the additional propagator
      subroutine gluonSelf(psq,musq,n,res,div)
      implicit none
      integer n,div
      real*8 psq,musq
      complex*16 res,b0fin,b0findiv
      external b0fin,b0findiv
      
      if (div.eq.0) then
	res = B0fin(psq,musq)*(5d0-2d0/3d0*n) + (3d0+2d0*n)/9d0
      else
	res = (-2d0/3d0*n+5d0)*B0finDiv(psq,musq,div)
      endif
      end

c massive fermion loop, including the additional propagator
      subroutine gluonSelfMass(psq,musq,m,res,div)
      implicit none
      integer div
      real*8 psq,musq,m,tiny
      parameter (tiny = 1d-7)
      complex*16 res,B0t1,B0tMDiv,A0fin
      external B0t1,B0tMDiv,A0fin
      
      if (div.eq.0) then
	if(abs(psq).lt.tiny) then
	  res= -2d0/3d0*B0t1(m,0d0,musq) ! do the calculation with psq=0 analytically (propagator 1/0)
	else
	  res= -2d0/9d0*(-1d0+6d0*(M*M-A0fin(M,musq))/psq+B0t1(M,psq,musq)*(3d0+6d0*M*M/psq))
	endif
      else
	res = -2d0/3d0*B0tMDiv(m,psq,div)
      endif
      end
      
c corrections to FFA vertex if ngluon=0; FFG vertex if ngluon=1
c color factors are included
      subroutine vertline(psq,musq,ngluon,res,div)
      implicit none
      integer ngluon, div
      real*8 psq,musq
      complex*16 res,B0fin,B0finDiv,C0fin,C0finDiv
      external B0fin,B0finDiv,C0fin,C0finDiv
      
      if(ngluon.eq.0) then
	if(div.eq.0) then
	  res = -4d0*B0fin(psq,musq)-8d0/3d0*(C0fin(psq,0d0,0d0,musq)*psq+1d0)
	else
	  res = -4d0*B0finDiv(psq,musq,div)-8d0/3d0*C0finDiv(psq,0d0,0d0,musq,div)*psq
	endif
      else if (ngluon.eq.1) then
	if(div.eq.0) then
	  res = -B0fin(psq,musq)+(1d0+psq*C0fin(psq,0d0,0d0,musq))/3d0
	else
	  res = -B0finDiv(psq,musq,div)+psq*C0finDiv(psq,0d0,0d0,musq,div)/3d0
	endif
      else
	print*,"error in vertline"
      endif
      end
      

c 3 Gluon Vertex; all gluons incoming; 1,2 onshell; 3 offshell
c returns born for div < 0
c result has to be multiplied by [t^a,t^b]
      subroutine Vert3G(p1,p2,p3,mup1,mup2,mup3,musq,result,div)
      implicit none
      real*8 p1(0:3),p2(0:3),p3(0:3),musq
      complex*16 mup1(0:3),mup2(0:3),mup3(0:3),result
      integer div
      
      complex*16 mup1mup2,mup1mup3,mup2mup3,p1mup2,p1mup3,p2mup1,p2mup3
      complex*16 p2m1m2m3,p1m2m1m3,pm3m1m2,b0,c0
      real*8 p3s
      integer fin
      
      real*8 dotrr
      complex*16 dotrc, dotcc,B0fin,B0finDiv,C0fin,C0finDiv
      external dotrc, dotcc,dotrr,B0fin,B0finDiv,C0fin,C0finDiv
      
      p3s = dotrr(p3,p3)
      
      mup1mup2 = dotcc(mup1,mup2)
      mup1mup3 = dotcc(mup1,mup3)
      mup2mup3 = dotcc(mup2,mup3)
      
      p1mup2 = dotrc(p1,mup2)
      p1mup3 = dotrc(p1,mup3)
      p2mup1 = dotrc(p2,mup1)
      p2mup3 = dotrc(p2,mup3)
      
      p2m1m2m3 = p2mup1*mup2mup3
      p1m2m1m3 = p1mup2*mup1mup3
      pm3m1m2  = (p1mup3-p2mup3)*mup1mup2
            
      if(div.lt.0) then ! born 
	result = 2d0*(p1m2m1m3-p2m1m2m3)-pm3m1m2
	
      else
	if(div.eq.0) then
	  fin=1
	  b0=B0fin(p3s,musq)
	  c0=C0fin(p3s,0d0,0d0,musq)*p3s
	else
	  fin=0
	  b0=B0finDiv(p3s,musq,div)
	  c0=C0finDiv(p3s,0d0,0d0,musq,div)*p3s
	endif
	result=( (8d0*fin+138d0*b0+36d0*c0)*(p1m2m1m3-p2m1m2m3)
     1	        +(8d0*fin- 78d0*b0-27d0*c0)*pm3m1m2)/(-12d0)
     2         +2d0/p3s*fin*p1mup2*p2mup1*(p1mup3-p2mup3)
      endif
      
      end

