c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c       subroutine dt3(p1s,p2s,s12,ratio)  
c      subroutine dt4(p1s,p2s,p3s,p1p2,p1p3,p2p3,ratio)  
c      subroutine dt4det(p1s,p2s,p3s,p1p2,p1p3,p2p3,det1)  
c       subroutine dt4detp(p1s,p2s,p3s,p1p2,p1p3,p2p3,det1,cdet,ratio)  
c      subroutine dt5(p1,p2,p3,p4,p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,ratio)
c      subroutine dt51m(msq,p1,p2,p3,p4,p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,ratio)
c      subroutine dt6(p1s,p2s,p3s,p4s,p5s,p1p2, p1p3,p1p4,p1p5,p2p3,
c      subroutine Print_Det_up_F(p1,p2,p3,p4,p5,p6)
c      subroutine Calc_Det_up_F(p1,p2,p3,p4,p5,p6,k,countbad,seedbad)
c      subroutine Calc_Det_up_D(p1,p2,p3,p4,k,badF)
c      subroutine Calc_Det_up_E(p1,p2,p3,p4,p5,countbad)
c

c$$$       p1sq = dotrr(p1,p1)
c$$$       p1p2 = dotrr(p1,p2)
c$$$       p1p3 = dotrr(p1,p3)
c$$$       p1p4 = dotrr(p1,p4)
c$$$       p1p5 = dotrr(p1,p5)
c$$$       p2sq = dotrr(p2,p2)
c$$$       p2p3 = dotrr(p2,p3)
c$$$       p2p4 = dotrr(p2,p4)
c$$$       p3sq = dotrr(p3,p3)
c$$$       p3p4 = dotrr(p3,p4)
c$$$       p4sq = dotrr(p4,p4)
c$$$       p4p5 = dotrr(p4,p5)
c$$$       p5sq = dotrr(p5,p5)
c$$$       p2p5= dotrr(p2,p5)
c$$$       p3p5= dotrr(p3,p5)
c$$$       s12 = (p1sq +p2sq+ 2*p1p2) 
c$$$       s15 = (p1sq +p5sq+ 2*p1p5) 
c$$$       s23 = (p2sq +p3sq+ 2*p2p3) 
c$$$       s34 = (p3sq +p4sq+ 2*p3p4) 
c$$$       s45 = (p4sq +p5sq+ 2*p4p5) 
c$$$c C0 functions
c$$$       call  det3(s45,p4sq,p5sq,
c$$$     \p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,i,det33,1,count1)
c$$$      call  det3(p1sq,s15,p5sq,
c$$$     \p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,i,det33,2,count1)
c$$$       call  det3(p1sq,p2sq,s12,
c$$$     \p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,i,det33,3,count1)
c$$$       call  det3(p1sq,s23,s45,
c$$$     \p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,i,det33,4,count1)
c$$$       call  det3(s12,s34,p5sq,
c$$$     \p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,i,det33,5,count1)
c$$$       call  det3(s12,p3sq,s45,
c$$$     \p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,i,det33,6,count1)
c$$$       call  det3(p2sq,p3sq,s23,
c$$$     \p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,i,det33,7,count1)
c$$$       call  det3(p2sq,s34,s15,
c$$$     \p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,i,det33,8,count1)
c$$$       call  det3(s23,p4sq,s15,
c$$$     \p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,i,det33,9,count1)
c$$$       call  det3(p3sq,p4sq,s34,
c$$$     \p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,i,det33,10,count1)
c$$$
c$$$c D0 functions
c$$$
c$$$      call  det(p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,
c$$$     \p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,i,det1,0,det5!,EE0
c$$$     \ ,count) 
c$$$       call  det(s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4,p3p4,
c$$$     \p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,i,det1,1,det5!,EE0
c$$$     \ ,count) 
c$$$       call  det(p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p4+p3p4,
c$$$     \p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,i,det1,2,det5!,EE0
c$$$     \ ,count) 
c$$$       call  det(p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p3+p2p4,
c$$$     \p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,i,det1,3,det5!, EE0
c$$$     \ ,count) 
c$$$       call  det(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,
c$$$     \p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,i,det1,4,det5!, EE0
c$$$     \  ,count) 
c$$$
c$$$
c$$$
c$$$      subroutine det(p1s,p2s,p3s,p1p2,p1p3,p2p3,
c$$$     \ p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,i,det1,j,det5!,EE0
c$$$     \ ,count)  
c$$$      implicit none
c$$$      real * 8 p1s, p2s, p3s, p1p2, p1p3,p2p3,det1,cdet,det5
c$$$      Real*8 p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15
c$$$      Complex*16 EE0
c$$$      Integer i,j,count
c$$$      Character*17 wr2
c$$$      Parameter(wr2='(A4,F28.16)')
c$$$
c$$$      det1=-2.d0*(-2.d0*p1p2*p1p3*p2p3 + p1p3*p1p3*p2s + p1p2*p1p2*p3s 
c$$$     -  + p1s*(p2p3*p2p3 - p2s*p3s))
c$$$
c$$$      cdet =2.d0*(abs(-2.d0*p1p2*p1p3*p2p3) +abs( p1p3*p1p3*p2s) 
c$$$     -  +abs(p1p2*p1p2*p3s)+abs(p1s*p2p3*p2p3)+abs(p1s*p2s*p3s))
c$$$      
c$$$      det5 = 2.d0*s15*(p4sq*s12*(s12*s23 - p2sq*s45) 
c$$$     - +s45*(p3sq*s12*s15 - s12*s23*s34 + p2sq*s34*s45))
c$$$    
c$$$      If ((Abs(det1)/cdet.lt.1.E-02)) then
c$$$c      .and.(Abs(det1)/cdet.gt.1.E-05)) then
c$$$       count=count+1
c$$$c       write(*,*)  '**********************'  
c$$$       write(*,*)  '**********************'  
c$$$       Print*, 'count', count
c$$$       Print* ,'D0 deter=', det1, 
c$$$     - '    cdeter=',cdet
c$$$       Print*, '  deter/cdeter=',det1/cdet
c$$$       write(*,*)     'randon i= ',i, '    box j=',j
c$$$       write(*,*)  'det5', det5 
c$$$c$$$       write(*,*)    'EE0' , EE0
c$$$c$$$       write(*,wr2)  'p1sq', p1sq 
c$$$c$$$       write(*,wr2)  'p2sq' ,p2sq
c$$$c$$$       write(*,wr2)  'p3sq', p3sq
c$$$c$$$       write(*,wr2)  'p4sq', p4sq
c$$$c$$$       write(*,wr2)  'p5sq', p5sq
c$$$c$$$       write(*,wr2)  's12' ,s12
c$$$c$$$       write(*,wr2)  's23' ,s23
c$$$c$$$       write(*,wr2)  's34' ,s34
c$$$c$$$       write(*,wr2)  's45' ,s45   
c$$$c$$$       write(*,wr2)  's15' ,s15
c$$$      endif
c$$$      return
c$$$      end
c$$$
c$$$
c$$$
c$$$       subroutine det3(p1s,p2s,ss12,
c$$$     \ p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,i,det33,j,count)  
c$$$      implicit none
c$$$      real * 8 p1s, p2s, p3s,det33,cdet,p1p2,ss12
c$$$      Real*8 p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15
c$$$      Complex*16 EE0
c$$$      Integer i,j,count
c$$$      Character*17 wr2
c$$$      Parameter(wr2='(A4,F28.16)')
c$$$       
c$$$      p1p2 = (ss12 - p1s - p2s)*0.5d0
c$$$      
c$$$      det33= 2.d0*(p1s*p2s - p1p2*p1p2)
c$$$
c$$$      cdet = 2.d0*(abs(p1s*p2s)+ abs( p1p2*p1p2))
c$$$      
c$$$c      det5 = 2.d0*s15*(p4sq*s12*(s12*s23 - p2sq*s45) 
c$$$c     - +s45*(p3sq*s12*s15 - s12*s23*s34 + p2sq*s34*s45))
c$$$    
c$$$      If ((Abs(det33)/cdet.lt.1.E-02)) then
c$$$c      .and.(Abs(det1)/cdet.gt.1.E-05)) then
c$$$       count=count+1
c$$$c       write(*,*)  '**********************'  
c$$$       write(*,*)  '**********************'  
c$$$       Print*, 'count', count
c$$$       Print* ,'C0 deter=', det33, 
c$$$     - '    cdeter=',cdet
c$$$       Print*, '  deter/cdeter=',det33/cdet
c$$$       write(*,*)     'randon i= ',i, '    box j=',j
c$$$c       write(*,*)  'det5', det5 
c$$$c       write(*,*)    'EE0' , EE0
c$$$c$$$       write(*,wr2)  'p1sq', p1sq 
c$$$c$$$       write(*,wr2)  'p2sq' ,p2sq
c$$$c$$$       write(*,wr2)  'p3sq', p3sq
c$$$c$$$       write(*,wr2)  'p4sq', p4sq
c$$$c$$$       write(*,wr2)  'p5sq', p5sq
c$$$c$$$       write(*,wr2)  's12' ,s12
c$$$c$$$       write(*,wr2)  's23' ,s23
c$$$c$$$       write(*,wr2)  's34' ,s34
c$$$c$$$       write(*,wr2)  's45' ,s45   
c$$$c$$$       write(*,wr2)  's15' ,s15
c$$$      endif
c$$$      return
c$$$      end
c$$$
c$$$
c$$$
c$$$c$$$
c$$$c$$$     subroutine det(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,i,Det)  
c$$$c$$$      implicit none
c$$$c$$$      real * 8 p1sq, p2sq, p3sq, p1p2, p1p3,p2p3,det,cdet
c$$$c$$$      Integer i
c$$$c$$$
c$$$c$$$      det =-2.d0*(-2.d0*p1p2*p1p3*p2p3 + p1p3*p1p3*p2sq + p1p2*p1p2*p3sq 
c$$$c$$$     -  + p1sq*(p2p3*p2p3 - p2sq*p3sq))
c$$$c$$$
c$$$c$$$      cdet =2.d0*(abs(-2.d0*p1p2*p1p3*p2p3) +abs( p1p3*p1p3*p2sq) 
c$$$c$$$     -  +abs(p1p2*p1p2*p3sq)+abs(p1sq*p2p3*p2p3)+abs(p1sq*p2sq*p3sq))
c$$$c$$$
c$$$c$$$      If (Abs(det)/cdet.lt.1.E-08)  Print* ,'D0 deter=', det, 
c$$$c$$$     - '    cdeter',cdet, '      deter/cdeter',det/cdet, 'i',i
c$$$c$$$      
c$$$c$$$      return
c$$$c$$$      end


      subroutine dt3(p1s,p2s,s12,ratio)  
      implicit none
      real * 8 p1s, p2s, s12,ratio,p1p2,det33,cdet
c      Character*17 wr2
c      Parameter(wr2='(A4,F28.16)')
      p1p2 = (s12 - p1s - p2s)*0.5d0
      det33= 2.d0*(p1s*p2s - p1p2*p1p2)
      cdet = 2.d0*(abs(p1s*p2s)+ abs(p1p2*p1p2))
      ratio=abs(det33)/abs(cdet)      
      return
      end


      subroutine dt4(p1s,p2s,p3s,p1p2,p1p3,p2p3,ratio)  
      implicit none
      real * 8 p1s, p2s, p3s, p1p2, p1p3,p2p3,det1,cdet,ratio
      det1=-2.d0*(-2.d0*p1p2*p1p3*p2p3 + p1p3*p1p3*p2s + p1p2*p1p2*p3s 
     -  + p1s*(p2p3*p2p3 - p2s*p3s))
      cdet =2.d0*(abs(-2.d0*p1p2*p1p3*p2p3) +abs( p1p3*p1p3*p2s) 
     -  +abs(p1p2*p1p2*p3s)+abs(p1s*p2p3*p2p3)+abs(p1s*p2s*p3s))
      ratio=abs(det1)/abs(cdet)    
      return
      end

      subroutine dt4det(p1s,p2s,p3s,p1p2,p1p3,p2p3,det1)  
      implicit none
      real * 8 p1s, p2s, p3s, p1p2, p1p3,p2p3,det1,cdet,ratio
      det1=-2.d0*(-2.d0*p1p2*p1p3*p2p3 + p1p3*p1p3*p2s + p1p2*p1p2*p3s 
     -  + p1s*(p2p3*p2p3 - p2s*p3s))
      cdet =2.d0*(abs(-2.d0*p1p2*p1p3*p2p3) +abs( p1p3*p1p3*p2s) 
     -  +abs(p1p2*p1p2*p3s)+abs(p1s*p2p3*p2p3)+abs(p1s*p2s*p3s))
      ratio=abs(det1)/abs(cdet)    
      return
      end

       subroutine dt4detp(p1s,p2s,p3s,p1p2,p1p3,p2p3,det1,cdet,ratio)  
      implicit none
      real * 8 p1s, p2s, p3s, p1p2, p1p3,p2p3,det1,cdet,ratio
      det1=-2.d0*(-2.d0*p1p2*p1p3*p2p3 + p1p3*p1p3*p2s + p1p2*p1p2*p3s 
     -  + p1s*(p2p3*p2p3 - p2s*p3s))
      cdet =2.d0*(abs(-2.d0*p1p2*p1p3*p2p3) +abs( p1p3*p1p3*p2s) 
     -  +abs(p1p2*p1p2*p3s)+abs(p1s*p2p3*p2p3)+abs(p1s*p2s*p3s))
      ratio=abs(det1)/abs(cdet)    
      return
      end


      subroutine dt5(p1,p2,p3,p4,p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,ratio)
      implicit none
      real*8 p1,p2,p3,p4,p5,p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
      real*8 p12,p22,p32,p42,p52,s12s,s23s,s34s,s45s,s15s
      real*8 s12,s23,s34,s45,s15
      real*8  d,cd,ratio

      s12=2.d0*p1p2+p2+p1
      s23=2.d0*p2p3+p2+p3      
      s34=2.d0*p3p4+p3+p4
      s45=p1+p2+p3+2d0*(p1p2+p1p3+p2p3)
      s15=p2+p3+p4+2d0*(p2p3+p2p4+p3p4)
      p5=s15+p1+2d0*(p1p2+p1p3+p1p4)

      p12=p1*p1
      p22=p2*p2
      p32=p3*p3
      p42=p4*p4
      p52=p5*p5

      s12s=s12*s12
      s23s=s23*s23
      s34s=s34*s34
      s45s=s45*s45
      s15s=s15*s15    

       d= 2*(p12*p3*p4*s34+s23*((p3*p5-p4*s12)*(p2*p5-s12*s15)+p5*s12*
     -   s23*s34)+((p2*p4-p3*s15)*(p2*p5-s12*s15)-(p2*p5+s12*s15)*s23
     -   *s34)*s45-p1*(-(p32*p5*s15)+p3*p4*s12*s15+p3*p5*s23*s34+p4*s
     -   12*s23*s34+s34*(p3*s15-s23*s34)*s45+p2*p4*(p3*p5-p4*s12+s34*
     -   s45))+p2*s15*s34*s45s)

       cd=2*Abs(p1*p2*p3*p4*p5) + 2*Abs(p1*p2*p4**2*s12) + 2*Abs(p1*p32*p5*s15) + 
     -  2*Abs(p1*p3*p4*s12*s15) + 2*Abs(p2*p3*p5**2*s23) + 2*Abs(p2*p4*p5*s12*s23) + 
     -  2*Abs(p3*p5*s12*s15*s23) + 2*Abs(p4*s12**2*s15*s23) + 2*Abs(p12*p3*p4*s34) + 
     -  2*Abs(p1*p3*p5*s23*s34) + 2*Abs(p1*p4*s12*s23*s34) + 2*Abs(p5*s12*s23**2*s34) + 
     -  2*Abs(p2**2*p4*p5*s45) + 2*Abs(p2*p3*p5*s15*s45) + 2*Abs(p2*p4*s12*s15*s45) + 
     -  2*Abs(p3*s12*s15**2*s45) + 2*Abs(p1*p2*p4*s34*s45) + 2*Abs(p1*p3*s15*s34*s45) + 
     -  2*Abs(p2*p5*s23*s34*s45) + 2*Abs(s12*s15*s23*s34*s45) + 2*Abs(p1*s23*s34**2*s45) + 
     -  2*Abs(p2*s15*s34*s45s)
      ratio=abs(d)/abs(cd)
      return
      end


      subroutine dt51m(msq,p1,p2,p3,p4,p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,ratio)
c Author: Francisco Campanario
c date: 25.08.2009
c determinte Cayley determinants and the sum of the absolute of their terms  
      implicit none
      real*8 p1,p2,p3,p4,p5,p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
      real*8 p12,p22,p32,p42,p52,s12s,s23s,s34s,s45s,s15s
      real*8 s12,s23,s34,s45,s15
      real*8  d,cd,ratio,msq

      s12=2.d0*p1p2+p2+p1
      s23=2.d0*p2p3+p2+p3      
      s34=2.d0*p3p4+p3+p4
      s45=p1+p2+p3+2d0*(p1p2+p1p3+p2p3)
      s15=p2+p3+p4+2d0*(p2p3+p2p4+p3p4)
      p5=s15+p1+2d0*(p1p2+p1p3+p1p4)

      p12=p1*p1
      p22=p2*p2
      p32=p3*p3
      p42=p4*p4
      p52=p5*p5

      s12s=s12*s12
      s23s=s23*s23
      s34s=s34*s34
      s45s=s45*s45
      s15s=s15*s15    

       d= 2*(16*msq*(p1p3**2*p2p4**2 - p1*p2p4**2*p3 + p1p4**2*(p2p3**2 - p2*p3) + 
     -       p1p2**2*p3p4**2 - p1*p2*p3p4**2 - 
     -       2*(p1p2*p1p3*p2p4*p3p4 - p1*p2p3*p2p4*p3p4 + 
     -          p1p4*(p1p3*p2p3*p2p4 - p1p2*p2p4*p3 - p1p3*p2*p3p4 + p1p2*p2p3*p3p4)) - 
     -       (p1p3**2*p2 - 2*p1p2*p1p3*p2p3 + p1*p2p3**2 + p1p2**2*p3 - p1*p2*p3)*p4) + 
     -    p1**2*p3*p4*s34 + s23*((p3*p5 - p4*s12)*(p2*p5 - s12*s15) + p5*s12*s23*s34) + 
     -    ((p2*p4 - p3*s15)*(p2*p5 - s12*s15) - (p2*p5 + s12*s15)*s23*s34)*s45 - 
     -    p1*(-(p3**2*p5*s15) + p3*p4*s12*s15 + p3*p5*s23*s34 + p4*s12*s23*s34 + 
     -       s34*(p3*s15 - s23*s34)*s45 + p2*p4*(p3*p5 - p4*s12 + s34*s45)) + p2*s15*s34*s45s) 


       cd=   2*(16*Abs(msq*p1p4**2*p2p3**2) + 32*Abs(msq*p1p3*p1p4*p2p3*p2p4) + 
     -    16*Abs(msq*p1p3**2*p2p4**2) + 16*Abs(msq*p1p4**2*p2*p3) + 
     -    32*Abs(msq*p1p2*p1p4*p2p4*p3) + 16*Abs(msq*p1*p2p4**2*p3) + 
     -    32*Abs(msq*p1p3*p1p4*p2*p3p4) + 32*Abs(msq*p1p2*p1p4*p2p3*p3p4) + 
     -    32*Abs(msq*p1p2*p1p3*p2p4*p3p4) + 32*Abs(msq*p1*p2p3*p2p4*p3p4) + 
     -    16*Abs(msq*p1p2**2*p3p4**2) + 16*Abs(msq*p1*p2*p3p4**2) + 16*Abs(msq*p1p3**2*p2*p4) + 
     -    32*Abs(msq*p1p2*p1p3*p2p3*p4) + 16*Abs(msq*p1*p2p3**2*p4) + 
     -    16*Abs(msq*p1p2**2*p3*p4) + 16*Abs(msq*p1*p2*p3*p4) + Abs(p1*p2*p3*p4*p5) + 
     -    Abs(p1*p2*p4**2*s12) + Abs(p1*p3**2*p5*s15) + Abs(p1*p3*p4*s12*s15) + 
     -    Abs(p2*p3*p5**2*s23) + Abs(p2*p4*p5*s12*s23) + Abs(p3*p5*s12*s15*s23) + 
     -    Abs(p4*s12**2*s15*s23) + Abs(p1**2*p3*p4*s34) + Abs(p1*p3*p5*s23*s34) + 
     -    Abs(p1*p4*s12*s23*s34) + Abs(p5*s12*s23**2*s34) + Abs(p2**2*p4*p5*s45) + 
     -    Abs(p2*p3*p5*s15*s45) + Abs(p2*p4*s12*s15*s45) + Abs(p3*s12*s15**2*s45) + 
     -    Abs(p1*p2*p4*s34*s45) + Abs(p1*p3*s15*s34*s45) + Abs(p2*p5*s23*s34*s45) + 
     -    Abs(s12*s15*s23*s34*s45) + Abs(p1*s23*s34**2*s45) + Abs(p2*s15*s34*s45s))


      ratio=abs(d)/abs(cd)
      return
      end


      subroutine dt6(p1s,p2s,p3s,p4s,p5s,p1p2, p1p3,p1p4,p1p5,p2p3,
     #p2p4,p2p5,p3p4,p3p5,p4p5,ratio)
      implicit none 
      real*8 p1s,p2s,p3s,p4s,p5s,p1p2,p1p3,p1p4,p1p5 
      real*8 p2p3,p2p4,p2p5,p3p4,p3p5,p4p5 
      real*8 p6s,p1p6,p5p6,s12,s23,s34,s45,s56,s16,s123,s234,s345 
      real*8 p1s2,p2s2,p3s2,p4s2,p5s2,p6s2,s12s,s23s,s34s,s45s 
      real*8 s56s,s16s,s123s,s234s,s345s 
      real*8 x1,x2,x3,x4,x5,X,Cx,cx1,cx2,cx3,cx4,cx5,ratio
      Integer k 

      p6s=p1s+p2s+p3s+p4s+p5s+2*(p1p2+p1p3+p1p4 
     -+p1p5+p2p3+p2p4+p2p5+p3p4+p3p5+p4p5) 
      p5p6=-(p1p5+p2p5+p3p5+p4p5+p5s)  
      p1p6=-(p1s+p1p2+p1p3+p1p4+p1p5)  
      s12=p1s+p2s+2*p1p2 
      s23=p2s+p3s+2*p2p3 
      s34=p3s+p4s+2*p3p4 
      s45=p4s+p5s+2*p4p5 
      s56=p5s+p6s+2*p5p6 
      s16=p1s+p6s+2*p1p6 
      s123=p1s+p2s+p3s+2*(p1p2+p1p3+p2p3) 
      s234=p2s+p3s+p4s+2*(p2p3+p2p4+p3p4) 
      s345=p3s+p4s+p5s+2*(p3p4+p3p5+p4p5) 
      p1s2=p1s*p1s 
      p2s2=p2s*p2s 
      p3s2=p3s*p3s 
      p4s2=p4s*p4s 
      p5s2=p5s*p5s 
      p6s2=p6s*p6s 
      s12s=s12*s12 
      s23s=s23*s23 
      s34s=s34*s34 
      s45s=s45*s45 
      s56s=s56*s56 
      s16s=s16*s16 
      s123s=s123*s123 
      s234s=s234*s234 
      s345s=s345*s345

      x1=-(p1s*(-4*p5s*(p2s2*p4s+p3s2*s234+s23*(p4s*(s234-s34)+s34*(s
     -   23-s234+s34))-p3s*(p4s*(s234-s34)+s234*(-s234+s34)+s23*(s234
     -   +s34))-p2s*(-p4s2+p3s*(p4s-s23+s234)+(s23-s234)*s34+p4s*(s23
     -   +s234+s34)))-(-p4s-p5s+s45)*(-((-p3s-p4s+s34)*((p2s+p3s-s23)
     -   *(s16-s234+s34-s345)+2*p2s*(p4s-s34+s345-s45)))+(p3s-s23+s23
     -   4-s34)*(-2*p3s*(s16-s234+s34-s345)-(p2s+p3s-s23)*(p4s-s34+s3
     -   45-s45))+(4*p2s*p3s-(p2s+p3s-s23)**2)*(-p4s-p5s+s45))-(s16-s
     -   234+s34-s345)*(-2*p4s*(-2*p3s*(s16-s234+s34-s345)-(p2s+p3s-s
     -   23)*(p4s-s34+s345-s45))-(-p3s2+p2s*(p3s+p4s-s34)+s23*(-p4s+s
     -   34)+p3s*(p4s+s23-2*s234+s34))*(p4s+p5s-s45)-(p3s+p4s-s34)*(-
     -   (s16*s34)+s23*s34+p4s*(s16-s23-s345)-s23*s345+s234*s345+p3s*
     -   (p4s+s16-s234-s45)+s23*s45-s234*s45+s34*s45))+(p4s-s34+s345-
     -   s45)*(-2*p4s*((p2s+p3s-s23)*(s16-s234+s34-s345)+2*p2s*(p4s-s
     -   34+s345-s45))+((p3s-s23)*(p3s-s23+s234-s34)+p2s*(-p3s-2*p4s-
     -   s23+s234+s34))*(-p4s-p5s+s45)+(p3s-s23+s234-s34)*(-(s16*s34)
     -   +s23*s34+p4s*(s16-s23-s345)-s23*s345+s234*s345+p3s*(p4s+s16-
     -   s234-s45)+s23*s45-s234*s45+s34*s45))))-(p6s-s56)*(2*(p2s2*p4
     -   s+p3s2*s234+s23*(p4s*(s234-s34)+s34*(s23-s234+s34))-p3s*(p4s
     -   *(s234-s34)+s234*(-s234+s34)+s23*(s234+s34))-p2s*(-p4s2+p3s*
     -   (p4s-s23+s234)+(s23-s234)*s34+p4s*(s23+s234+s34)))*(p6s-s16+
     -   s234-s56)+(-p4s-p5s+s45)*(-((-p2s2+p1s*(p2s+p3s-s23)+s12*(-p
     -   3s+s23)+p2s*(p3s+s12-2*s123+s23))*(p3s+p4s-s34))-(p2s2-2*p1s
     -   *p3s+p3s*s12+p3s*s123-p3s*s23+s12*s23-s123*s23-p2s*(p3s+s12-
     -   s123+2*s23)+s23s)*(p3s-s23+s234-s34)+(4*p2s*p3s-(p2s+p3s-s23
     -   )**2)*(-s123+s23-s234+s56))+(s16-s234+s34-s345)*(2*p4s*(p2s2
     -   -2*p1s*p3s+p3s*s12+p3s*s123-p3s*s23+s12*s23-s123*s23-p2s*(p3
     -   s+s12-s123+2*s23)+s23s)-(-p3s-p4s+s34)*((p1s+p2s-s12)*(p3s+p
     -   4s-s34)+(p2s-s12+s123-s23)*(-p3s+s23-s234+s34))+(-p3s2+p2s*(
     -   p3s+p4s-s34)+s23*(-p4s+s34)+p3s*(p4s+s23-2*s234+s34))*(-s123
     -   +s23-s234+s56))-(p4s-s34+s345-s45)*(2*p4s*(-p2s2+p1s*(p2s+p3
     -   s-s23)+s12*(-p3s+s23)+p2s*(p3s+s12-2*s123+s23))-(p3s-s23+s23
     -   4-s34)*((p1s+p2s-s12)*(p3s+p4s-s34)+(p2s-s12+s123-s23)*(-p3s
     -   +s23-s234+s34))+((p3s-s23)*(p3s-s23+s234-s34)+p2s*(-p3s-2*p4
     -   s-s23+s234+s34))*(-s123+s23-s234+s56)))+(-s123+s56)*(-((-((-
     -   p3s-p4s+s34)*((p2s+p3s-s23)*(s16-s234+s34-s345)+2*p2s*(p4s-s
     -   34+s345-s45)))+(p3s-s23+s234-s34)*(-2*p3s*(s16-s234+s34-s345
     -   )-(p2s+p3s-s23)*(p4s-s34+s345-s45))+(4*p2s*p3s-(p2s+p3s-s23)
     -   **2)*(-p4s-p5s+s45))*(p6s-s16+s234-s56))+2*p5s*(-((-p2s2+p1s
     -   *(p2s+p3s-s23)+s12*(-p3s+s23)+p2s*(p3s+s12-2*s123+s23))*(p3s
     -   +p4s-s34))-(p2s2-2*p1s*p3s+p3s*s12+p3s*s123-p3s*s23+s12*s23-
     -   s123*s23-p2s*(p3s+s12-s123+2*s23)+s23s)*(p3s-s23+s234-s34)+(
     -   4*p2s*p3s-(p2s+p3s-s23)**2)*(-s123+s23-s234+s56))-(p4s-s34+s
     -   345-s45)*(-((p3s-s23+s234-s34)*(-((p2s-s12+s123-s23)*(s16-s2
     -   34+s34-s345))-(p1s+p2s-s12)*(p4s-s34+s345-s45)))-(-p2s2+p1s*
     -   (p2s+p3s-s23)+s12*(-p3s+s23)+p2s*(p3s+s12-2*s123+s23))*(p4s+
     -   p5s-s45)+((p2s+p3s-s23)*(s16-s234+s34-s345)+2*p2s*(p4s-s34+s
     -   345-s45))*(-s123+s23-s234+s56))+(s16-s234+s34-s345)*(-((-p3s
     -   -p4s+s34)*(-((p2s-s12+s123-s23)*(s16-s234+s34-s345))-(p1s+p2
     -   s-s12)*(p4s-s34+s345-s45)))-(p2s2-2*p1s*p3s+p3s*s12+p3s*s123
     -   -p3s*s23+s12*s23-s123*s23-p2s*(p3s+s12-s123+2*s23)+s23s)*(p4
     -   s+p5s-s45)+(-2*p3s*(s16-s234+s34-s345)-(p2s+p3s-s23)*(p4s-s3
     -   4+s345-s45))*(-s123+s23-s234+s56)))-(-s12+s123)*(-((-2*p4s*(
     -   (p2s+p3s-s23)*(s16-s234+s34-s345)+2*p2s*(p4s-s34+s345-s45))+
     -   ((p3s-s23)*(p3s-s23+s234-s34)+p2s*(-p3s-2*p4s-s23+s234+s34))
     -   *(-p4s-p5s+s45)+(p3s-s23+s234-s34)*(-(s16*s34)+s23*s34+p4s*(
     -   s16-s23-s345)-s23*s345+s234*s345+p3s*(p4s+s16-s234-s45)+s23*
     -   s45-s234*s45+s34*s45))*(p6s-s16+s234-s56))+2*p5s*(2*p4s*(-p2
     -   s2+p1s*(p2s+p3s-s23)+s12*(-p3s+s23)+p2s*(p3s+s12-2*s123+s23)
     -   )-(p3s-s23+s234-s34)*((p1s+p2s-s12)*(p3s+p4s-s34)+(p2s-s12+s
     -   123-s23)*(-p3s+s23-s234+s34))+((p3s-s23)*(p3s-s23+s234-s34)+
     -   p2s*(-p3s-2*p4s-s23+s234+s34))*(-s123+s23-s234+s56))-(-p4s-p
     -   5s+s45)*(-((p3s-s23+s234-s34)*(-((p2s-s12+s123-s23)*(s16-s23
     -   4+s34-s345))-(p1s+p2s-s12)*(p4s-s34+s345-s45)))-(-p2s2+p1s*(
     -   p2s+p3s-s23)+s12*(-p3s+s23)+p2s*(p3s+s12-2*s123+s23))*(p4s+p
     -   5s-s45)+((p2s+p3s-s23)*(s16-s234+s34-s345)+2*p2s*(p4s-s34+s3
     -   45-s45))*(-s123+s23-s234+s56))+(s16-s234+s34-s345)*(-2*p4s*(
     -   -((p2s-s12+s123-s23)*(s16-s234+s34-s345))-(p1s+p2s-s12)*(p4s
     -   -s34+s345-s45))+((p1s+p2s-s12)*(p3s+p4s-s34)+(p2s-s12+s123-s
     -   23)*(-p3s+s23-s234+s34))*(-p4s-p5s+s45)+(-(s16*s34)+s23*s34+
     -   p4s*(s16-s23-s345)-s23*s345+s234*s345+p3s*(p4s+s16-s234-s45)
     -   +s23*s45-s234*s45+s34*s45)*(-s123+s23-s234+s56)))+(-p1s+s12)
     -   *(-((-2*p4s*(-2*p3s*(s16-s234+s34-s345)-(p2s+p3s-s23)*(p4s-s
     -   34+s345-s45))-(-p3s2+p2s*(p3s+p4s-s34)+s23*(-p4s+s34)+p3s*(p
     -   4s+s23-2*s234+s34))*(p4s+p5s-s45)-(p3s+p4s-s34)*(-(s16*s34)+
     -   s23*s34+p4s*(s16-s23-s345)-s23*s345+s234*s345+p3s*(p4s+s16-s
     -   234-s45)+s23*s45-s234*s45+s34*s45))*(p6s-s16+s234-s56))+2*p5
     -   s*(2*p4s*(p2s2-2*p1s*p3s+p3s*s12+p3s*s123-p3s*s23+s12*s23-s1
     -   23*s23-p2s*(p3s+s12-s123+2*s23)+s23s)-(-p3s-p4s+s34)*((p1s+p
     -   2s-s12)*(p3s+p4s-s34)+(p2s-s12+s123-s23)*(-p3s+s23-s234+s34)
     -   )+(-p3s2+p2s*(p3s+p4s-s34)+s23*(-p4s+s34)+p3s*(p4s+s23-2*s23
     -   4+s34))*(-s123+s23-s234+s56))-(-p4s-p5s+s45)*(-((-p3s-p4s+s3
     -   4)*(-((p2s-s12+s123-s23)*(s16-s234+s34-s345))-(p1s+p2s-s12)*
     -   (p4s-s34+s345-s45)))-(p2s2-2*p1s*p3s+p3s*s12+p3s*s123-p3s*s2
     -   3+s12*s23-s123*s23-p2s*(p3s+s12-s123+2*s23)+s23s)*(p4s+p5s-s
     -   45)+(-2*p3s*(s16-s234+s34-s345)-(p2s+p3s-s23)*(p4s-s34+s345-
     -   s45))*(-s123+s23-s234+s56))+(p4s-s34+s345-s45)*(-2*p4s*(-((p
     -   2s-s12+s123-s23)*(s16-s234+s34-s345))-(p1s+p2s-s12)*(p4s-s34
     -   +s345-s45))+((p1s+p2s-s12)*(p3s+p4s-s34)+(p2s-s12+s123-s23)*
     -   (-p3s+s23-s234+s34))*(-p4s-p5s+s45)+(-(s16*s34)+s23*s34+p4s*
     -   (s16-s23-s345)-s23*s345+s234*s345+p3s*(p4s+s16-s234-s45)+s23
     -   *s45-s234*s45+s34*s45)*(-s123+s23-s234+s56)))
       x2=-((-s123+s56)*(-((-((p2s2-2*p1s*p3s+p3s*s12+p3s*s123-p3s*s23
     -   +s12*s23-s123*s23-p2s*(p3s+s12-s123+2*s23)+s23s)*(p4s+p5s-s4
     -   5))+(p3s-s23+s234-s34)*((p2s-s12+s123-s23)*(p4s-s34+s345-s45
     -   )-2*p3s*(p6s-s16+s234-s56))-(-p3s-p4s+s34)*(-((p1s+p2s-s12)*
     -   (p4s-s34+s345-s45))+(p2s+p3s-s23)*(p6s-s16+s234-s56)))*(p6s-
     -   s16+s234-s56))+2*p5s*(-((4*p1s*p3s-(p2s-s12+s123-s23)**2)*(p
     -   3s-s23+s234-s34))+((p2s-s12)*(p2s-s12+s123-s23)+p1s*(-p2s-2*
     -   p3s-s12+s123+s23))*(-p3s-p4s+s34)+(p2s2-2*p1s*p3s+p3s*s12+p3
     -   s*s123-p3s*s23+s12*s23-s123*s23-p2s*(p3s+s12-s123+2*s23)+s23
     -   s)*(-s123+s23-s234+s56))+(s16-s234+s34-s345)*((4*p1s*p3s-(p2
     -   s-s12+s123-s23)**2)*(-p4s-p5s+s45)-(-p3s-p4s+s34)*(2*p1s*(p4
     -   s-s34+s345-s45)-(p2s-s12+s123-s23)*(p6s-s16+s234-s56))+((p2s
     -   -s12+s123-s23)*(p4s-s34+s345-s45)-2*p3s*(p6s-s16+s234-s56))*
     -   (-s123+s23-s234+s56))-(p4s-s34+s345-s45)*(((p2s-s12)*(p2s-s1
     -   2+s123-s23)+p1s*(-p2s-2*p3s-s12+s123+s23))*(-p4s-p5s+s45)-(p
     -   3s-s23+s234-s34)*(2*p1s*(p4s-s34+s345-s45)-(p2s-s12+s123-s23
     -   )*(p6s-s16+s234-s56))+(-((p1s+p2s-s12)*(p4s-s34+s345-s45))+(
     -   p2s+p3s-s23)*(p6s-s16+s234-s56))*(-s123+s23-s234+s56))))+(p6
     -   s-s56)*((-p4s-p5s+s45)*(-((4*p1s*p3s-(p2s-s12+s123-s23)**2)*
     -   (p3s-s23+s234-s34))+((p2s-s12)*(p2s-s12+s123-s23)+p1s*(-p2s-
     -   2*p3s-s12+s123+s23))*(-p3s-p4s+s34)+(p2s2-2*p1s*p3s+p3s*s12+
     -   p3s*s123-p3s*s23+s12*s23-s123*s23-p2s*(p3s+s12-s123+2*s23)+s
     -   23s)*(-s123+s23-s234+s56))+(s16-s234+s34-s345)*(2*p4s*(4*p1s
     -   *p3s-(p2s-s12+s123-s23)**2)-(-p3s-p4s+s34)*(-2*p1s*(p3s+p4s-
     -   s34)+(p2s-s12+s123-s23)*(s123-s23+s234-s56))+(-((p2s-s12+s12
     -   3-s23)*(p3s+p4s-s34))+2*p3s*(s123-s23+s234-s56))*(-s123+s23-
     -   s234+s56))-(p6s-s16+s234-s56)*(2*p4s*(p2s2-2*p1s*p3s+p3s*s12
     -   +p3s*s123-p3s*s23+s12*s23-s123*s23-p2s*(p3s+s12-s123+2*s23)+
     -   s23s)+(p3s-s23+s234-s34)*(-((p2s-s12+s123-s23)*(p3s+p4s-s34)
     -   )+2*p3s*(s123-s23+s234-s56))-(-p3s-p4s+s34)*((p1s+p2s-s12)*(
     -   p3s+p4s-s34)+(p2s+p3s-s23)*(-s123+s23-s234+s56)))-(p4s-s34+s
     -   345-s45)*(2*p4s*((p2s-s12)*(p2s-s12+s123-s23)+p1s*(-p2s-2*p3
     -   s-s12+s123+s23))-(p3s-s23+s234-s34)*(-2*p1s*(p3s+p4s-s34)+(p
     -   2s-s12+s123-s23)*(s123-s23+s234-s56))+(-s123+s23-s234+s56)*(
     -   (p1s+p2s-s12)*(p3s+p4s-s34)+(p2s+p3s-s23)*(-s123+s23-s234+s5
     -   6))))+p1s*(-((-p4s-p5s+s45)*(-((p2s2-2*p1s*p3s+p3s*s12+p3s*s
     -   123-p3s*s23+s12*s23-s123*s23-p2s*(p3s+s12-s123+2*s23)+s23s)*
     -   (p4s+p5s-s45))+(p3s-s23+s234-s34)*((p2s-s12+s123-s23)*(p4s-s
     -   34+s345-s45)-2*p3s*(p6s-s16+s234-s56))-(-p3s-p4s+s34)*(-((p1
     -   s+p2s-s12)*(p4s-s34+s345-s45))+(p2s+p3s-s23)*(p6s-s16+s234-s
     -   56))))+2*p5s*(2*p4s*(p2s2-2*p1s*p3s+p3s*s12+p3s*s123-p3s*s23
     -   +s12*s23-s123*s23-p2s*(p3s+s12-s123+2*s23)+s23s)+(p3s-s23+s2
     -   34-s34)*(-((p2s-s12+s123-s23)*(p3s+p4s-s34))+2*p3s*(s123-s23
     -   +s234-s56))-(-p3s-p4s+s34)*((p1s+p2s-s12)*(p3s+p4s-s34)+(p2s
     -   +p3s-s23)*(-s123+s23-s234+s56)))+(p4s-s34+s345-s45)*(-2*p4s*
     -   (-((p1s+p2s-s12)*(p4s-s34+s345-s45))+(p2s+p3s-s23)*(p6s-s16+
     -   s234-s56))+(-p4s-p5s+s45)*((p1s+p2s-s12)*(p3s+p4s-s34)+(p2s+
     -   p3s-s23)*(-s123+s23-s234+s56))+(p3s-s23+s234-s34)*((p3s+p4s-
     -   s34)*(p6s-s16+s234-s56)+(p4s-s34+s345-s45)*(-s123+s23-s234+s
     -   56)))-(s16-s234+s34-s345)*(-2*p4s*((p2s-s12+s123-s23)*(p4s-s
     -   34+s345-s45)-2*p3s*(p6s-s16+s234-s56))+(-p4s-p5s+s45)*(-((p2
     -   s-s12+s123-s23)*(p3s+p4s-s34))+2*p3s*(s123-s23+s234-s56))+(-
     -   p3s-p4s+s34)*((p3s+p4s-s34)*(p6s-s16+s234-s56)+(p4s-s34+s345
     -   -s45)*(-s123+s23-s234+s56))))+(-s12+s123)*(-((-p4s-p5s+s45)*
     -   (((p2s-s12)*(p2s-s12+s123-s23)+p1s*(-p2s-2*p3s-s12+s123+s23)
     -   )*(-p4s-p5s+s45)-(p3s-s23+s234-s34)*(2*p1s*(p4s-s34+s345-s45
     -   )-(p2s-s12+s123-s23)*(p6s-s16+s234-s56))+(-((p1s+p2s-s12)*(p
     -   4s-s34+s345-s45))+(p2s+p3s-s23)*(p6s-s16+s234-s56))*(-s123+s
     -   23-s234+s56)))+2*p5s*(2*p4s*((p2s-s12)*(p2s-s12+s123-s23)+p1
     -   s*(-p2s-2*p3s-s12+s123+s23))-(p3s-s23+s234-s34)*(-2*p1s*(p3s
     -   +p4s-s34)+(p2s-s12+s123-s23)*(s123-s23+s234-s56))+(-s123+s23
     -   -s234+s56)*((p1s+p2s-s12)*(p3s+p4s-s34)+(p2s+p3s-s23)*(-s123
     -   +s23-s234+s56)))-(p6s-s16+s234-s56)*(-2*p4s*(-((p1s+p2s-s12)
     -   *(p4s-s34+s345-s45))+(p2s+p3s-s23)*(p6s-s16+s234-s56))+(-p4s
     -   -p5s+s45)*((p1s+p2s-s12)*(p3s+p4s-s34)+(p2s+p3s-s23)*(-s123+
     -   s23-s234+s56))+(p3s-s23+s234-s34)*((p3s+p4s-s34)*(p6s-s16+s2
     -   34-s56)+(p4s-s34+s345-s45)*(-s123+s23-s234+s56)))+(s16-s234+
     -   s34-s345)*(-2*p4s*(2*p1s*(p4s-s34+s345-s45)-(p2s-s12+s123-s2
     -   3)*(p6s-s16+s234-s56))+(-p4s-p5s+s45)*(-2*p1s*(p3s+p4s-s34)+
     -   (p2s-s12+s123-s23)*(s123-s23+s234-s56))+(-s123+s23-s234+s56)
     -   *((p3s+p4s-s34)*(p6s-s16+s234-s56)+(p4s-s34+s345-s45)*(-s123
     -   +s23-s234+s56))))-(-p1s+s12)*(-((-p4s-p5s+s45)*((4*p1s*p3s-(
     -   p2s-s12+s123-s23)**2)*(-p4s-p5s+s45)-(-p3s-p4s+s34)*(2*p1s*(
     -   p4s-s34+s345-s45)-(p2s-s12+s123-s23)*(p6s-s16+s234-s56))+((p
     -   2s-s12+s123-s23)*(p4s-s34+s345-s45)-2*p3s*(p6s-s16+s234-s56)
     -   )*(-s123+s23-s234+s56)))+2*p5s*(2*p4s*(4*p1s*p3s-(p2s-s12+s1
     -   23-s23)**2)-(-p3s-p4s+s34)*(-2*p1s*(p3s+p4s-s34)+(p2s-s12+s1
     -   23-s23)*(s123-s23+s234-s56))+(-((p2s-s12+s123-s23)*(p3s+p4s-
     -   s34))+2*p3s*(s123-s23+s234-s56))*(-s123+s23-s234+s56))-(p6s-
     -   s16+s234-s56)*(-2*p4s*((p2s-s12+s123-s23)*(p4s-s34+s345-s45)
     -   -2*p3s*(p6s-s16+s234-s56))+(-p4s-p5s+s45)*(-((p2s-s12+s123-s
     -   23)*(p3s+p4s-s34))+2*p3s*(s123-s23+s234-s56))+(-p3s-p4s+s34)
     -   *((p3s+p4s-s34)*(p6s-s16+s234-s56)+(p4s-s34+s345-s45)*(-s123
     -   +s23-s234+s56)))+(p4s-s34+s345-s45)*(-2*p4s*(2*p1s*(p4s-s34+
     -   s345-s45)-(p2s-s12+s123-s23)*(p6s-s16+s234-s56))+(-p4s-p5s+s
     -   45)*(-2*p1s*(p3s+p4s-s34)+(p2s-s12+s123-s23)*(s123-s23+s234-
     -   s56))+(-s123+s23-s234+s56)*((p3s+p4s-s34)*(p6s-s16+s234-s56)
     -   +(p4s-s34+s345-s45)*(-s123+s23-s234+s56))))
       x3=(-s123+s56)*(-((-((-p2s2+p1s*(p2s+p3s-s23)+s12*(-p3s+s23)+p2
     -   s*(p3s+s12-2*s123+s23))*(p4s+p5s-s45))-(-p3s-p4s+s34)*(-((p1
     -   s+p2s-s12)*(s16-s234+s34-s345))-2*p2s*(p6s-s16+s234-s56))+(p
     -   3s-s23+s234-s34)*((p2s-s12+s123-s23)*(s16-s234+s34-s345)+(p2
     -   s+p3s-s23)*(p6s-s16+s234-s56)))*(p6s-s16+s234-s56))+2*p5s*(-
     -   (((p2s-s12)*(p2s-s12+s123-s23)+p1s*(-p2s-2*p3s-s12+s123+s23)
     -   )*(p3s-s23+s234-s34))+(4*p1s*p2s-(p1s+p2s-s12)**2)*(-p3s-p4s
     -   +s34)+(-p2s2+p1s*(p2s+p3s-s23)+s12*(-p3s+s23)+p2s*(p3s+s12-2
     -   *s123+s23))*(-s123+s23-s234+s56))-(p4s-s34+s345-s45)*((4*p1s
     -   *p2s-(p1s+p2s-s12)**2)*(-p4s-p5s+s45)-(p3s-s23+s234-s34)*(2*
     -   p1s*(s16-s234+s34-s345)+(p1s+p2s-s12)*(p6s-s16+s234-s56))+(-
     -   ((p1s+p2s-s12)*(s16-s234+s34-s345))-2*p2s*(p6s-s16+s234-s56)
     -   )*(-s123+s23-s234+s56))+(s16-s234+s34-s345)*(((p2s-s12)*(p2s
     -   -s12+s123-s23)+p1s*(-p2s-2*p3s-s12+s123+s23))*(-p4s-p5s+s45)
     -   -(-p3s-p4s+s34)*(2*p1s*(s16-s234+s34-s345)+(p1s+p2s-s12)*(p6
     -   s-s16+s234-s56))+((p2s-s12+s123-s23)*(s16-s234+s34-s345)+(p2
     -   s+p3s-s23)*(p6s-s16+s234-s56))*(-s123+s23-s234+s56)))-(-s12+
     -   s123)*(-((-p4s-p5s+s45)*((4*p1s*p2s-(p1s+p2s-s12)**2)*(-p4s-
     -   p5s+s45)-(p3s-s23+s234-s34)*(2*p1s*(s16-s234+s34-s345)+(p1s+
     -   p2s-s12)*(p6s-s16+s234-s56))+(-((p1s+p2s-s12)*(s16-s234+s34-
     -   s345))-2*p2s*(p6s-s16+s234-s56))*(-s123+s23-s234+s56)))+2*p5
     -   s*(2*p4s*(4*p1s*p2s-(p1s+p2s-s12)**2)-(p3s-s23+s234-s34)*(2*
     -   p1s*(p3s-s23+s234-s34)-(p1s+p2s-s12)*(s123-s23+s234-s56))+(-
     -   ((p1s+p2s-s12)*(p3s-s23+s234-s34))+2*p2s*(s123-s23+s234-s56)
     -   )*(-s123+s23-s234+s56))-(p6s-s16+s234-s56)*(-2*p4s*(-((p1s+p
     -   2s-s12)*(s16-s234+s34-s345))-2*p2s*(p6s-s16+s234-s56))+(-p4s
     -   -p5s+s45)*(-((p1s+p2s-s12)*(p3s-s23+s234-s34))+2*p2s*(s123-s
     -   23+s234-s56))+(p3s-s23+s234-s34)*(-((p3s-s23+s234-s34)*(p6s-
     -   s16+s234-s56))+(s16-s234+s34-s345)*(-s123+s23-s234+s56)))+(s
     -   16-s234+s34-s345)*(-2*p4s*(2*p1s*(s16-s234+s34-s345)+(p1s+p2
     -   s-s12)*(p6s-s16+s234-s56))+(-p4s-p5s+s45)*(2*p1s*(p3s-s23+s2
     -   34-s34)-(p1s+p2s-s12)*(s123-s23+s234-s56))+(-s123+s23-s234+s
     -   56)*(-((p3s-s23+s234-s34)*(p6s-s16+s234-s56))+(s16-s234+s34-
     -   s345)*(-s123+s23-s234+s56))))-p1s*(-((-p4s-p5s+s45)*(-((-p2s
     -   2+p1s*(p2s+p3s-s23)+s12*(-p3s+s23)+p2s*(p3s+s12-2*s123+s23))
     -   *(p4s+p5s-s45))-(-p3s-p4s+s34)*(-((p1s+p2s-s12)*(s16-s234+s3
     -   4-s345))-2*p2s*(p6s-s16+s234-s56))+(p3s-s23+s234-s34)*((p2s-
     -   s12+s123-s23)*(s16-s234+s34-s345)+(p2s+p3s-s23)*(p6s-s16+s23
     -   4-s56))))+(p4s-s34+s345-s45)*(-2*p4s*(-((p1s+p2s-s12)*(s16-s
     -   234+s34-s345))-2*p2s*(p6s-s16+s234-s56))+(-p4s-p5s+s45)*(-((
     -   p1s+p2s-s12)*(p3s-s23+s234-s34))+2*p2s*(s123-s23+s234-s56))+
     -   (p3s-s23+s234-s34)*(-((p3s-s23+s234-s34)*(p6s-s16+s234-s56))
     -   +(s16-s234+s34-s345)*(-s123+s23-s234+s56)))+2*p5s*(2*p4s*(-p
     -   2s2+p1s*(p2s+p3s-s23)+s12*(-p3s+s23)+p2s*(p3s+s12-2*s123+s23
     -   ))-(-p3s-p4s+s34)*(-((p1s+p2s-s12)*(p3s-s23+s234-s34))+2*p2s
     -   *(s123-s23+s234-s56))+(p3s-s23+s234-s34)*(s12*s23-s12*s234+s
     -   123*s234+s12*s34-s123*s34+s23*s34-p3s*(s12+s234-s56)-s23*s56
     -   +p2s*(p3s-s123-s34+s56)))-(s16-s234+s34-s345)*(-2*p4s*((p2s-
     -   s12+s123-s23)*(s16-s234+s34-s345)+(p2s+p3s-s23)*(p6s-s16+s23
     -   4-s56))+(-p3s-p4s+s34)*(-((p3s-s23+s234-s34)*(p6s-s16+s234-s
     -   56))+(s16-s234+s34-s345)*(-s123+s23-s234+s56))-(p4s+p5s-s45)
     -   *(s12*s23-s12*s234+s123*s234+s12*s34-s123*s34+s23*s34-p3s*(s
     -   12+s234-s56)-s23*s56+p2s*(p3s-s123-s34+s56))))+(-p1s+s12)*(-
     -   ((-p4s-p5s+s45)*(((p2s-s12)*(p2s-s12+s123-s23)+p1s*(-p2s-2*p
     -   3s-s12+s123+s23))*(-p4s-p5s+s45)-(-p3s-p4s+s34)*(2*p1s*(s16-
     -   s234+s34-s345)+(p1s+p2s-s12)*(p6s-s16+s234-s56))+((p2s-s12+s
     -   123-s23)*(s16-s234+s34-s345)+(p2s+p3s-s23)*(p6s-s16+s234-s56
     -   ))*(-s123+s23-s234+s56)))+(p4s-s34+s345-s45)*(-2*p4s*(2*p1s*
     -   (s16-s234+s34-s345)+(p1s+p2s-s12)*(p6s-s16+s234-s56))+(-p4s-
     -   p5s+s45)*(2*p1s*(p3s-s23+s234-s34)-(p1s+p2s-s12)*(s123-s23+s
     -   234-s56))+(-s123+s23-s234+s56)*(-((p3s-s23+s234-s34)*(p6s-s1
     -   6+s234-s56))+(s16-s234+s34-s345)*(-s123+s23-s234+s56)))-(p6s
     -   -s16+s234-s56)*(-2*p4s*((p2s-s12+s123-s23)*(s16-s234+s34-s34
     -   5)+(p2s+p3s-s23)*(p6s-s16+s234-s56))+(-p3s-p4s+s34)*(-((p3s-
     -   s23+s234-s34)*(p6s-s16+s234-s56))+(s16-s234+s34-s345)*(-s123
     -   +s23-s234+s56))-(p4s+p5s-s45)*(s12*s23-s12*s234+s123*s234+s1
     -   2*s34-s123*s34+s23*s34-p3s*(s12+s234-s56)-s23*s56+p2s*(p3s-s
     -   123-s34+s56)))+2*p5s*(2*p4s*((p2s-s12)*(p2s-s12+s123-s23)+p1
     -   s*(-p2s-2*p3s-s12+s123+s23))-(-p3s-p4s+s34)*(2*p1s*(p3s-s23+
     -   s234-s34)-(p1s+p2s-s12)*(s123-s23+s234-s56))+(-s123+s23-s234
     -   +s56)*(s12*s23-s12*s234+s123*s234+s12*s34-s123*s34+s23*s34-p
     -   3s*(s12+s234-s56)-s23*s56+p2s*(p3s-s123-s34+s56))))-(p6s-s56
     -   )*((-p4s-p5s+s45)*(-(((p2s-s12)*(p2s-s12+s123-s23)+p1s*(-p2s
     -   -2*p3s-s12+s123+s23))*(p3s-s23+s234-s34))+(4*p1s*p2s-(p1s+p2
     -   s-s12)**2)*(-p3s-p4s+s34)+(-p2s2+p1s*(p2s+p3s-s23)+s12*(-p3s
     -   +s23)+p2s*(p3s+s12-2*s123+s23))*(-s123+s23-s234+s56))-(p4s-s
     -   34+s345-s45)*(2*p4s*(4*p1s*p2s-(p1s+p2s-s12)**2)-(p3s-s23+s2
     -   34-s34)*(2*p1s*(p3s-s23+s234-s34)-(p1s+p2s-s12)*(s123-s23+s2
     -   34-s56))+(-((p1s+p2s-s12)*(p3s-s23+s234-s34))+2*p2s*(s123-s2
     -   3+s234-s56))*(-s123+s23-s234+s56))-(p6s-s16+s234-s56)*(2*p4s
     -   *(-p2s2+p1s*(p2s+p3s-s23)+s12*(-p3s+s23)+p2s*(p3s+s12-2*s123
     -   +s23))-(-p3s-p4s+s34)*(-((p1s+p2s-s12)*(p3s-s23+s234-s34))+2
     -   *p2s*(s123-s23+s234-s56))+(p3s-s23+s234-s34)*(s12*s23-s12*s2
     -   34+s123*s234+s12*s34-s123*s34+s23*s34-p3s*(s12+s234-s56)-s23
     -   *s56+p2s*(p3s-s123-s34+s56)))+(s16-s234+s34-s345)*(2*p4s*((p
     -   2s-s12)*(p2s-s12+s123-s23)+p1s*(-p2s-2*p3s-s12+s123+s23))-(-
     -   p3s-p4s+s34)*(2*p1s*(p3s-s23+s234-s34)-(p1s+p2s-s12)*(s123-s
     -   23+s234-s56))+(-s123+s23-s234+s56)*(s12*s23-s12*s234+s123*s2
     -   34+s12*s34-s123*s34+s23*s34-p3s*(s12+s234-s56)-s23*s56+p2s*(
     -   p3s-s123-s34+s56))))
       x4=-(p3s*p4s*p6s*s12*s16)+2*p3s*p5s*s12*s123*s16-p4s2*s12s*s16+
     -   p3s*p4s*s12*s16s-2*p3s*s12*s123*s16s+p4s*s12*s123*s16s-p4s*s
     -   12s*s16s-p3s*p5s*p6s*s12*s23+p4s*p5s*s12s*s23-p3s*p5s*s12*s1
     -   6*s23+2*p3s*p6s*s12*s16*s23-p4s*p6s*s12*s16*s23-p5s*s12*s123
     -   *s16*s23+2*p4s*s12s*s16*s23+p5s*s12s*s16*s23-p3s2*p6s2*s234-
     -   p3s*p4s*p6s*s12*s234+2*p3s*p5s*s12*s123*s234+p3s2*p6s*s16*s2
     -   34-p3s*p4s*s12*s16*s234-p3s*p6s*s12*s16*s234-p3s*p6s*s123*s1
     -   6*s234+4*p3s*s12*s123*s16*s234-p4s*s12*s123*s16*s234+p4s*s12
     -   s*s16*s234+p3s*p6s2*s23*s234-p3s*p5s*s12*s23*s234-2*p3s*p6s*
     -   s12*s23*s234+2*p4s*p6s*s12*s23*s234-p5s*s12*s123*s23*s234-2*
     -   p4s*s12s*s23*s234+p5s*s12s*s23*s234-p3s2*p6s*s234s+p3s*p6s*s
     -   12*s234s+p3s*p6s*s123*s234s-2*p3s*s12*s123*s234s+p5s*p6s*s12
     -   *s23s-2*p5s*s12s*s23s-p3s*p6s*s123*s16*s34+2*p4s*s12*s123*s1
     -   6*s34+p3s*s123*s16s*s34+s12*s123*s16s*s34-s123s*s16s*s34+p3s
     -   *p6s2*s23*s34-p4s*p6s*s12*s23*s34-p5s*s12*s123*s23*s34-p3s*p
     -   6s*s16*s23*s34-p4s*s12*s16*s23*s34-p6s*s12*s16*s23*s34+2*p6s
     -   *s123*s16*s23*s34-2*s12*s123*s16*s23*s34-p3s*p6s*s123*s234*s
     -   34-p3s*s123*s16*s234*s34-s12*s123*s16*s234*s34+s123s*s16*s23
     -   4*s34+2*p3s*p6s*s23*s234*s34-p6s*s12*s23*s234*s34-p6s*s123*s
     -   23*s234*s34+2*s12*s123*s23*s234*s34-p6s2*s23s*s34+p5s*s12*s2
     -   3s*s34+2*p6s*s12*s23s*s34-p4s*s12*s123*s16*s345+2*p4s*p6s*s1
     -   2*s23*s345-p5s*s12*s123*s23*s345-p4s*s12*s16*s23*s345+2*s12*
     -   s123*s16*s23*s345+2*p3s*p6s*s123*s234*s345-p4s*s12*s123*s234
     -   *s345-p3s*s123*s16*s234*s345-s12*s123*s16*s234*s345+s123s*s1
     -   6*s234*s345-p3s*p6s*s23*s234*s345+2*p4s*s12*s23*s234*s345+2*
     -   p6s*s12*s23*s234*s345-p6s*s123*s23*s234*s345-2*s12*s123*s23*
     -   s234*s345+p3s*s123*s234s*s345+s12*s123*s234s*s345-s123s*s234
     -   s*s345+p5s*s12*s23s*s345-2*p6s*s12*s23s*s345+s123s*s16*s34*s
     -   345-p6s*s123*s23*s34*s345-s123*s16*s23*s34*s345+s123s*s234*s
     -   34*s345-s123*s23*s234*s34*s345+p6s*s23s*s34*s345-s123s*s234*
     -   s345s+s123*s23*s234*s345s-s123s*s16*s34s+p6s*s123*s23*s34s+s
     -   123*s16*s23*s34s-p6s*s23s*s34s+p4s*s12s*s16*s45+p5s*s12s*s23
     -   *s45-2*s12s*s16*s23*s45+2*p3s*p6s*s12*s234*s45+p4s*s12s*s234
     -   *s45-p3s*s12*s16*s234*s45-s12*s123*s16*s234*s45+s12s*s16*s23
     -   4*s45-p6s*s12*s23*s234*s45+2*s12s*s23*s234*s45+p3s*s12*s234s
     -   *s45+s12*s123*s234s*s45-s12s*s234s*s45-s12*s123*s16*s34*s45-
     -   p6s*s12*s23*s34*s45+2*s12*s16*s23*s34*s45-s12*s123*s234*s34*
     -   s45-s12*s23*s234*s34*s45+2*s12*s123*s234*s345*s45-s12*s23*s2
     -   34*s345*s45-p1s2*(2*p3s2*p5s+p4s2*s345+s34*s45*(s34-s345+s45
     -   )-p4s*(s345*(-s345+s45)+s34*(s345+s45))-p3s*(p4s*(p5s-2*s34+
     -   2*s345)+2*(s34-s345)*s45+p5s*(s34+s345+s45)))-s12s*s234*s45s
     -   +p3s2*p6s*s16*s56+2*p3s*p4s*s12*s16*s56-p3s2*s16s*s56+p3s*s1
     -   2*s16s*s56+p3s*s123*s16s*s56-p3s*p5s*s12*s23*s56-p3s*p6s*s16
     -   *s23*s56-2*p3s*s12*s16*s23*s56-p4s*s12*s16*s23*s56+p3s2*p6s*
     -   s234*s56+p3s2*s16*s234*s56-p3s*s12*s16*s234*s56-p3s*s123*s16
     -   *s234*s56-p3s*p6s*s23*s234*s56+2*p3s*s12*s23*s234*s56+p5s*s1
     -   2*s23s*s56+2*p3s*s123*s16*s34*s56-p3s*p6s*s23*s34*s56-p3s*s1
     -   6*s23*s34*s56+2*s12*s16*s23*s34*s56-s123*s16*s23*s34*s56+p6s
     -   *s23s*s34*s56-2*s12*s23s*s34*s56-p3s*s123*s16*s345*s56-p3s*p
     -   6s*s23*s345*s56-p4s*s12*s23*s345*s56+2*p3s*s16*s23*s345*s56-
     -   s12*s16*s23*s345*s56-s123*s16*s23*s345*s56-p3s*s123*s234*s34
     -   5*s56-p3s*s23*s234*s345*s56-s12*s23*s234*s345*s56+2*s123*s23
     -   *s234*s345*s56+p6s*s23s*s345*s56+2*s12*s23s*s345*s56-s123*s2
     -   3*s34*s345*s56+s23s*s34*s345*s56+s123*s23*s345s*s56-s23s*s34
     -   5s*s56-p3s*s12*s16*s45*s56+2*s12*s16*s23*s45*s56-p3s*s12*s23
     -   4*s45*s56-s12*s23*s234*s45*s56+2*s12*s23*s34*s45*s56-s12*s23
     -   *s345*s45*s56-p2s2*(p4s2*p6s-p5s*s123*(p6s-2*s123+s45+s56)+s
     -   45*(p6s*(2*s123-s56)+s56*(-2*s123+s45+s56))+p4s*(p6s2-p5s*s1
     -   23+2*s123*s56-s45*s56-p6s*(2*s123+s45+s56)))-p1s*(-(p4s2*s12
     -   *s16)+p4s*p5s*s12*s23+p4s*s12*s16*s34+p4s*s123*s16*s34+p4s*p
     -   6s*s23*s34-2*p4s*s12*s23*s34+p5s*s12*s23*s34-2*p5s*s123*s23*
     -   s34-p4s2*s12*s345-2*p4s*s12*s16*s345+p4s*s123*s16*s345+p4s*p
     -   6s*s23*s345+2*p4s*s12*s23*s345+p5s*s12*s23*s345-2*p5s*s123*s
     -   23*s345+p4s*s12*s234*s345-2*p4s*s123*s234*s345+p4s*s123*s34*
     -   s345+s123*s16*s34*s345+p4s*s23*s34*s345+p6s*s23*s34*s345-4*s
     -   123*s23*s34*s345+s123*s234*s34*s345-p4s*s123*s345s-p4s*s23*s
     -   345s+2*s123*s23*s345s-s123*s234*s345s-s123*s16*s34s-p6s*s23*
     -   s34s+2*s123*s23*s34s+p4s*s12*s16*s45+p5s*s12*s23*s45+p4s*s12
     -   *s234*s45+p4s*s12*s34*s45+s12*s16*s34*s45-2*s123*s16*s34*s45
     -   -2*p6s*s23*s34*s45+2*s12*s23*s34*s45-2*s12*s234*s34*s45+s123
     -   *s234*s34*s45+p4s*s12*s345*s45-2*s12*s23*s345*s45+s12*s234*s
     -   345*s45+s123*s234*s345*s45+s123*s34*s345*s45+s23*s34*s345*s4
     -   5-s123*s34s*s45-s23*s34s*s45-s12*s234*s45s-s12*s34*s45s-2*p4
     -   s*s23*s345*s56+s23*s34*s345*s56-s23*s345s*s56+s23*s34*s45*s5
     -   6+s23*s345*s45*s56-p3s2*(-2*(s16-s234)*(p6s-s56)+p5s*(p6s+s1
     -   6+s234+s56))+p3s*(-2*p6s*s16*s34+2*s123*s16*s34+2*p6s*s23*s3
     -   4+p6s*s234*s34-2*s123*s234*s34-2*s123*s16*s345-2*p6s*s23*s34
     -   5+p6s*s234*s345+2*s123*s234*s345-2*s12*s16*s45+p6s*s234*s45+
     -   2*s12*s234*s45+p6s*s34*s45+s16*s34*s45+s234*s34*s45-2*s234*s
     -   345*s45+s16*s34*s56-2*s23*s34*s56+s16*s345*s56+2*s23*s345*s5
     -   6-2*s234*s345*s56+s16*s45*s56-2*s234*s45*s56+s34*s45*s56-2*s
     -   345*s45*s56+p4s*(p5s*s12+2*s12*s16-2*s12*s234-2*s16*s34+s16*
     -   s345+s234*s345+p6s*(-2*s16+s234-2*s34+s345)+s16*s56+s345*s56
     -   )+p5s*(s123*(s16+s234+s34+s345)+s12*(s16-4*s23+s234+s45)+s23
     -   *(p6s+s34+s345+s56)))+p2s*(p5s*s123*s34+p5s*s123*s345-p4s2*(
     -   p6s-2*s12+s345)-2*p5s*s12*s45+p5s*s123*s45+p6s*s34*s45+2*s12
     -   3*s34*s45-2*s123*s345*s45+2*s12*s45s-s34*s45s-2*s34*s45*s56+
     -   s345*s45*s56-s45s*s56+p4s*(p5s*(-2*s12+s123)-2*s123*s34+2*s1
     -   23*s345-4*s12*s45+s34*s45+s345*s45+p6s*(s34-2*s345+s45)+s345
     -   *s56+s45*s56)+p3s*(p4s*(p5s+2*p6s-2*s56)+2*s45*(-p6s+s56)+p5
     -   s*(p6s-4*s123+s45+s56))))-p3s2*s16*s56s+p3s*s16*s23*s56s+p3s
     -   *s23*s345*s56s-s23s*s345*s56s-p2s*(p5s*s12*s123*s16-p5s*s123
     -   s*s16-p4s2*s12*(p6s+s16)+p5s*p6s*s12*s23+p5s*p6s*s123*s23-4*
     -   p5s*s12*s123*s23+p5s*s12*s123*s234-p5s*s123s*s234-p5s*s123s*
     -   s34+p6s*s123*s16*s34-2*s123s*s16*s34-p6s2*s23*s34+p5s*s123*s
     -   23*s34+2*p6s*s123*s23*s34-2*p6s*s123*s234*s34+2*s123s*s234*s
     -   34-p5s*s123s*s345+2*s123s*s16*s345+p5s*s123*s23*s345-2*p6s*s
     -   123*s23*s345+p6s*s123*s234*s345-2*s123s*s234*s345+p5s*s12*s1
     -   23*s45-2*s12*s123*s16*s45+p5s*s12*s23*s45-2*p6s*s12*s23*s45+
     -   p6s*s12*s234*s45-2*p6s*s123*s234*s45+2*s12*s123*s234*s45-2*p
     -   6s*s123*s34*s45+s123*s16*s34*s45+p6s*s23*s34*s45-2*s123*s234
     -   *s34*s45+s123*s234*s345*s45-s12*s234*s45s+p5s*s12*s23*s56+p5
     -   s*s123*s23*s56+s123*s16*s34*s56+p6s*s23*s34*s56-2*s123*s23*s
     -   34*s56-2*s123*s16*s345*s56+p6s*s23*s345*s56+2*s123*s23*s345*
     -   s56+s123*s234*s345*s56+s12*s16*s45*s56+s123*s16*s45*s56+p6s*
     -   s23*s45*s56+2*s12*s23*s45*s56-2*s12*s234*s45*s56+s123*s234*s
     -   45*s56+s123*s34*s45*s56+s23*s34*s45*s56+s123*s345*s45*s56-2*
     -   s23*s345*s45*s56-s12*s45s*s56+p4s*(2*s12*s123*s16-p6s2*s23+p
     -   5s*s12*(s123+s23)-2*s12*s123*s234+s123*s16*s34-2*s123*s16*s3
     -   45+s123*s234*s345+s12*s16*s45+s12*s234*s45+s12*s16*s56-2*s12
     -   3*s16*s56-2*s12*s23*s56-2*s123*s345*s56+s23*s345*s56+s12*s45
     -   *s56+p6s*(s123*(s16+s234+s34+s345)+s12*(-2*s16+2*s23+s234+s4
     -   5)+s23*(-2*s34+s345+s56)))-s23*s345*s56s-s23*s45*s56s+p3s*(-
     -   2*p6s*s123*s16+2*p6s2*s23-p6s2*s234+2*p6s*s123*s234+p6s*s234
     -   *s45+p6s*s16*s56+2*s123*s16*s56-4*p6s*s23*s56+p6s*s234*s56-2
     -   *s123*s234*s56+p6s*s45*s56-2*s16*s45*s56+s234*s45*s56+p4s*(-
     -   p6s2+s16*s56+p6s*(s16-2*s234+s56))+p5s*(p6s*(s123-2*s23)-2*s
     -   23*s56+s123*(s16+s234+s56))-s16*s56s+2*s23*s56s-s45*s56s))
       x5=-(p4s2*s12s*s16)+p4s*p5s*s12s*s23+p4s*s12s*s16*s23-p3s*p4s*p
     -   6s*s12*s234+2*p3s*p5s*s12*s123*s234-p3s*p4s*s12*s16*s234+2*p
     -   3s*s12*s123*s16*s234-p4s*s12*s123*s16*s234+p4s*s12s*s16*s234
     -   -p3s*p5s*s12*s23*s234-p3s*p6s*s12*s23*s234+2*p4s*p6s*s12*s23
     -   *s234-p5s*s12*s123*s23*s234-2*p4s*s12s*s23*s234+p5s*s12s*s23
     -   *s234-p3s2*p6s*s234s+p3s*p6s*s12*s234s+p3s*p6s*s123*s234s-2*
     -   p3s*s12*s123*s234s-p5s*s12s*s23s+2*p4s*s12*s123*s16*s34-p4s*
     -   p6s*s12*s23*s34-p5s*s12*s123*s23*s34-p4s*s12*s16*s23*s34-s12
     -   *s123*s16*s23*s34-p3s*p6s*s123*s234*s34-p3s*s123*s16*s234*s3
     -   4-s12*s123*s16*s234*s34+s123s*s16*s234*s34+2*p3s*p6s*s23*s23
     -   4*s34-p6s*s12*s23*s234*s34-p6s*s123*s23*s234*s34+2*s12*s123*
     -   s23*s234*s34+p5s*s12*s23s*s34+p6s*s12*s23s*s34-p4s*s12*s123*
     -   s234*s345+2*p4s*s12*s23*s234*s345-s12*s123*s23*s234*s345+p3s
     -   *s123*s234s*s345+s12*s123*s234s*s345-s123s*s234s*s345+s123s*
     -   s234*s34*s345-s123*s23*s234*s34*s345-s123s*s16*s34s+p6s*s123
     -   *s23*s34s+s123*s16*s23*s34s-p6s*s23s*s34s+p4s*s12s*s234*s45+
     -   s12s*s23*s234*s45+p3s*s12*s234s*s45+s12*s123*s234s*s45-s12s*
     -   s234s*s45-s12*s123*s234*s34*s45-s12*s23*s234*s34*s45-p1s2*(p
     -   3s2*p5s+(p4s-s34)*(p4s*s345-s34*s45)-p3s*(p4s*(p5s-2*s34+s34
     -   5)+s34*(p5s+s45)))+2*p3s*p4s*s12*s16*s56-p3s*p5s*s12*s23*s56
     -   -p3s*s12*s16*s23*s56-p4s*s12*s16*s23*s56+p3s2*p6s*s234*s56+p
     -   3s2*s16*s234*s56-p3s*s12*s16*s234*s56-p3s*s123*s16*s234*s56-
     -   p3s*p6s*s23*s234*s56+2*p3s*s12*s23*s234*s56+p5s*s12*s23s*s56
     -   +2*p3s*s123*s16*s34*s56-p3s*p6s*s23*s34*s56-p3s*s16*s23*s34*
     -   s56+2*s12*s16*s23*s34*s56-s123*s16*s23*s34*s56+p6s*s23s*s34*
     -   s56-2*s12*s23s*s34*s56-p4s*s12*s23*s345*s56-p3s*s123*s234*s3
     -   45*s56-p3s*s23*s234*s345*s56-s12*s23*s234*s345*s56+2*s123*s2
     -   3*s234*s345*s56+s12*s23s*s345*s56-s123*s23*s34*s345*s56+s23s
     -   *s34*s345*s56-p3s*s12*s234*s45*s56-s12*s23*s234*s45*s56+2*s1
     -   2*s23*s34*s45*s56-p2s2*(p4s2*p6s+(s123-s56)*(p5s*s123-s45*s5
     -   6)-p4s*(p5s*s123+(-2*s123+s45)*s56+p6s*(s123+s56)))-p1s*(-(p
     -   4s2*s12*s16)+p4s*p5s*s12*s23+p4s*s12*s16*s34+p4s*s123*s16*s3
     -   4+p4s*p6s*s23*s34-2*p4s*s12*s23*s34+p5s*s12*s23*s34-2*p5s*s1
     -   23*s23*s34-p4s2*s12*s345+p4s*s12*s23*s345+p4s*s12*s234*s345-
     -   2*p4s*s123*s234*s345+p4s*s123*s34*s345+p4s*s23*s34*s345-2*s1
     -   23*s23*s34*s345+s123*s234*s34*s345-s123*s16*s34s-p6s*s23*s34
     -   s+2*s123*s23*s34s+p4s*s12*s234*s45+p4s*s12*s34*s45+s12*s23*s
     -   34*s45-2*s12*s234*s34*s45+s123*s234*s34*s45-s123*s34s*s45-s2
     -   3*s34s*s45-2*p4s*s23*s345*s56+s23*s34*s345*s56+s23*s34*s45*s
     -   56-p3s2*(p6s*s234+(s16-2*s234)*s56+p5s*(s234+s56))+p2s*(-(p4
     -   s2*(p6s-2*s12+s345))+s34*(p5s*s123+s123*s45-2*s45*s56)+p3s*(
     -   -2*p5s*s123+p4s*(p5s+p6s-2*s56)+p5s*s56+s45*s56)+p4s*(-2*p5s
     -   *s12+p5s*s123+p6s*s34-2*s123*s34+s123*s345-2*s12*s45+s34*s45
     -   +s345*s56+s45*s56))+p3s*(s123*s16*s34+p6s*s23*s34+p6s*s234*s
     -   34-2*s123*s234*s34+s123*s234*s345+s12*s234*s45+s234*s34*s45+
     -   s16*s34*s56-2*s23*s34*s56+s23*s345*s56-2*s234*s345*s56-2*s23
     -   4*s45*s56+s34*s45*s56+p4s*(p5s*s12+s12*(s16-2*s234)+p6s*s234
     -   -2*p6s*s34-2*s16*s34+s234*s345+s16*s56+s345*s56)+p5s*(s12*(-
     -   2*s23+s234)+s123*(s234+s34)+s23*(s34+s56))))-p3s2*s16*s56s+p
     -   3s*s16*s23*s56s+p3s*s23*s345*s56s-s23s*s345*s56s-p2s*(-(p4s2
     -   *s12*(p6s+s16))+p3s*p6s*s123*s234-s123s*s16*s34+p6s*s123*s23
     -   *s34-2*p6s*s123*s234*s34+2*s123s*s234*s34-s123s*s234*s345+s1
     -   2*s123*s234*s45-2*s123*s234*s34*s45+p3s*s123*s16*s56-2*p3s*p
     -   6s*s23*s56+p3s*p6s*s234*s56-2*p3s*s123*s234*s56+s123*s16*s34
     -   *s56+p6s*s23*s34*s56-2*s123*s23*s34*s56+s123*s23*s345*s56+s1
     -   23*s234*s345*s56+s12*s23*s45*s56+p3s*s234*s45*s56-2*s12*s234
     -   *s45*s56+s123*s234*s45*s56+s123*s34*s45*s56+s23*s34*s45*s56+
     -   p4s*(p5s*s12*(s123+s23)-2*p3s*p6s*s234+p6s*s123*s234+p6s*s12
     -   3*s34+s123*s16*s34-2*p6s*s23*s34+s123*s234*s345+p3s*p6s*s56+
     -   p3s*s16*s56-2*s123*s16*s56+p6s*s23*s56-2*s123*s345*s56+s23*s
     -   345*s56+s12*(s123*(s16-2*s234)+p6s*(s23+s234)+s234*s45+s16*s
     -   56-2*s23*s56+s45*s56))+p5s*(s12*(s123*(-2*s23+s234)+s23*s56)
     -   +p3s*(-2*s23*s56+s123*(s234+s56))+s123*(-(s123*(s234+s34))+s
     -   23*(s34+s56)))-p3s*s16*s56s+2*p3s*s23*s56s-s23*s345*s56s-p3s
     -   *s45*s56s-s23*s45*s56s)


       cx1=Abs(p1s*p2s**2*p4s**2)+Abs(p1s*p2s2*p4s**2)+2*Abs(p1s**2*p3
     -   s**2*p5s)+2*Abs(p1s*p2s**2*p4s*p5s)+2*Abs(p1s*p2s2*p4s*p5s)+
     -   2*Abs(p1s**2*p3s*p4s*p5s)+Abs(p1s*p2s*p3s*p4s*p5s)+4*Abs(p1s
     -   *p2s*p4s**2*p5s)+4*Abs(p1s*p2s*p4s2*p5s)+Abs(p1s*p2s**2*p5s*
     -   *2)+Abs(p1s*p2s2*p5s**2)+2*Abs(p1s**2*p3s*p5s**2)+Abs(p1s*p2
     -   s*p3s*p5s**2)+Abs(p1s*p3s**2*p5s**2)+2*Abs(p1s*p2s*p3s*p4s*p
     -   6s)+Abs(p1s*p3s**2*p4s*p6s)+Abs(p1s*p3s2*p4s*p6s)+2*Abs(p1s*
     -   p2s*p4s**2*p6s)+Abs(p2s2*p4s**2*p6s)+2*Abs(p1s*p2s*p3s*p5s*p
     -   6s)+Abs(p1s*p3s**2*p5s*p6s)+Abs(p1s*p3s2*p5s*p6s)+2*Abs(p1s*
     -   p2s*p4s*p5s*p6s)+Abs(p2s2*p4s*p5s*p6s)+2*Abs(p1s*p3s*p4s*p5s
     -   *p6s)+Abs(p2s*p3s*p4s*p5s*p6s)+2*Abs(p2s2*p4s*p6s**2)+2*Abs(
     -   p2s*p3s*p4s*p6s**2)+2*Abs(p2s*p4s2*p6s**2)+2*Abs(p1s*p2s*p4s
     -   **2*s12)+4*Abs(p1s*p2s*p4s*p5s*s12)+2*Abs(p1s*p3s*p4s*p5s*s1
     -   2)+2*Abs(p1s*p2s*p5s**2*s12)+2*Abs(p1s*p3s*p5s**2*s12)+Abs(p
     -   3s**2*p4s*p6s*s12)+Abs(p3s2*p4s*p6s*s12)+2*Abs(p2s*p4s**2*p6
     -   s*s12)+Abs(p3s**2*p5s*p6s*s12)+Abs(p3s2*p5s*p6s*s12)+2*Abs(p
     -   2s*p4s*p5s*p6s*s12)+2*Abs(p3s*p4s*p5s*p6s*s12)+4*Abs(p1s*p2s
     -   *p3s*p5s*s123)+2*Abs(p1s*p3s**2*p5s*s123)+2*Abs(p1s*p3s2*p5s
     -   *s123)+2*Abs(p1s*p2s*p4s*p5s*s123)+Abs(p2s2*p4s*p5s*s123)+2*
     -   Abs(p1s*p2s*p5s**2*s123)+Abs(p2s2*p5s**2*s123)+2*Abs(p1s*p3s
     -   *p5s**2*s123)+Abs(p2s*p3s*p5s**2*s123)+2*Abs(p2s**2*p4s*p6s*
     -   s123)+2*Abs(p2s**2*p5s*p6s*s123)+2*Abs(p2s*p3s*p5s*p6s*s123)
     -   +4*Abs(p2s*p4s*p5s*p6s*s123)+2*Abs(p3s**2*p5s*s12*s123)+2*Ab
     -   s(p3s2*p5s*s12*s123)+2*Abs(p2s*p4s*p5s*s12*s123)+2*Abs(p2s*p
     -   5s**2*s12*s123)+2*Abs(p3s*p5s**2*s12*s123)+2*Abs(p2s**2*p5s*
     -   s123**2)+2*Abs(p2s*p5s**2*s123**2)+Abs(p1s*p3s**2*p5s*s16)+A
     -   bs(p1s*p3s*p4s*p5s*s16)+2*Abs(p1s*p3s**2*p6s*s16)+4*Abs(p1s*
     -   p3s*p4s*p6s*s16)+Abs(p2s*p3s*p4s*p6s*s16)+2*Abs(p1s*p4s**2*p
     -   6s*s16)+Abs(p2s*p4s**2*p6s*s16)+2*Abs(p2s*p4s2*p6s*s16)+2*Ab
     -   s(p1s*p3s*p4s*s12*s16)+Abs(p3s**2*p4s*s12*s16)+Abs(p3s2*p4s*
     -   s12*s16)+2*Abs(p1s*p4s**2*s12*s16)+Abs(p2s*p4s**2*s12*s16)+2
     -   *Abs(p1s*p3s*p5s*s12*s16)+Abs(p3s**2*p5s*s12*s16)+Abs(p3s2*p
     -   5s*s12*s16)+2*Abs(p1s*p4s*p5s*s12*s16)+Abs(p2s*p4s*p5s*s12*s
     -   16)+Abs(p3s*p4s*p5s*s12*s16)+4*Abs(p2s*p4s*p6s*s12*s16)+2*Ab
     -   s(p3s*p4s*p6s*s12*s16)+2*Abs(p4s**2*p6s*s12*s16)+2*Abs(p4s**
     -   2*s12**2*s16)+2*Abs(p4s*p5s*s12**2*s16)+Abs(p2s**2*p4s*s123*
     -   s16)+Abs(p2s2*p4s*s123*s16)+Abs(p2s**2*p5s*s123*s16)+Abs(p2s
     -   2*p5s*s123*s16)+2*Abs(p1s*p3s*p5s*s123*s16)+Abs(p2s*p3s*p5s*
     -   s123*s16)+2*Abs(p1s*p4s*p5s*s123*s16)+Abs(p2s*p4s*p5s*s123*s
     -   16)+2*Abs(p2s*p3s*p6s*s123*s16)+Abs(p3s**2*p6s*s123*s16)+Abs
     -   (p3s2*p6s*s123*s16)+2*Abs(p2s*p4s*p6s*s123*s16)+2*Abs(p2s*p4
     -   s*s12*s123*s16)+2*Abs(p2s*p5s*s12*s123*s16)+4*Abs(p3s*p5s*s1
     -   2*s123*s16)+2*Abs(p4s*p5s*s12*s123*s16)+2*Abs(p2s*p5s*s123**
     -   2*s16)+Abs(p3s*p4s*s12*s16**2)+Abs(p4s**2*s12*s16**2)+2*Abs(
     -   p4s*s12**2*s16**2)+2*Abs(p3s*s12*s123*s16**2)+2*Abs(p4s*s12*
     -   s123*s16**2)+2*Abs(p1s*p3s**2*p5s*s23)+2*Abs(p1s*p3s2*p5s*s2
     -   3)+Abs(p1s*p3s*p5s**2*s23)+Abs(p2s**2*p4s*p6s*s23)+Abs(p2s2*
     -   p4s*p6s*s23)+Abs(p2s**2*p5s*p6s*s23)+Abs(p2s2*p5s*p6s*s23)+2
     -   *Abs(p1s*p3s*p5s*p6s*s23)+2*Abs(p2s*p3s*p5s*p6s*s23)+2*Abs(p
     -   1s*p4s*p5s*p6s*s23)+Abs(p2s*p4s*p5s*p6s*s23)+2*Abs(p2s*p3s*p
     -   6s**2*s23)+2*Abs(p2s*p4s*p6s**2*s23)+4*Abs(p1s*p3s*p5s*s12*s
     -   23)+2*Abs(p3s**2*p5s*s12*s23)+2*Abs(p3s2*p5s*s12*s23)+2*Abs(
     -   p1s*p4s*p5s*s12*s23)+Abs(p2s*p4s*p5s*s12*s23)+2*Abs(p1s*p5s*
     -   *2*s12*s23)+Abs(p2s*p5s**2*s12*s23)+Abs(p3s*p5s**2*s12*s23)+
     -   2*Abs(p2s*p4s*p6s*s12*s23)+2*Abs(p2s*p5s*p6s*s12*s23)+2*Abs(
     -   p3s*p5s*p6s*s12*s23)+2*Abs(p4s*p5s*p6s*s12*s23)+2*Abs(p4s*p5
     -   s*s12**2*s23)+2*Abs(p5s**2*s12**2*s23)+2*Abs(p2s**2*p5s*s123
     -   *s23)+2*Abs(p2s2*p5s*s123*s23)+2*Abs(p1s*p5s**2*s123*s23)+Ab
     -   s(p2s*p5s**2*s123*s23)+2*Abs(p2s*p5s*p6s*s123*s23)+4*Abs(p2s
     -   *p5s*s12*s123*s23)+2*Abs(p5s**2*s12*s123*s23)+Abs(p3s**2*p6s
     -   *s16*s23)+Abs(p3s2*p6s*s16*s23)+Abs(p3s*p5s*s12*s16*s23)+2*A
     -   bs(p4s*p5s*s12*s16*s23)+2*Abs(p3s*p6s*s12*s16*s23)+2*Abs(p4s
     -   *p6s*s12*s16*s23)+2*Abs(p4s*s12**2*s16*s23)+2*Abs(p5s*s12**2
     -   *s16*s23)+2*Abs(p5s*s12*s123*s16*s23)+Abs(p1s*p4s**2*s23**2)
     -   +2*Abs(p1s*p4s*p5s*s23**2)+Abs(p1s*p5s**2*s23**2)+Abs(p3s*p4
     -   s*p6s*s23**2)+Abs(p3s*p5s*p6s*s23**2)+Abs(p4s**2*s12*s23**2)
     -   +2*Abs(p4s*p5s*s12*s23**2)+2*Abs(p5s*p6s*s12*s23**2)+2*Abs(p
     -   5s*s12**2*s23**2)+2*Abs(p3s*p5s*s123*s23**2)+2*Abs(p4s*p6s*s
     -   16*s23**2)+Abs(p4s*s123*s16*s23**2)+Abs(p5s*s123*s16*s23**2)
     -   +Abs(p4s*p6s*s23**3)+Abs(p5s*p6s*s23**3)+2*Abs(p5s*s123*s23*
     -   *3)+Abs(p1s*p3s**2*p5s*s234)+2*Abs(p1s*p3s2*p5s*s234)+2*Abs(
     -   p1s*p3s**2*p6s*s234)+Abs(p2s**2*p4s*p6s*s234)+Abs(p2s2*p4s*p
     -   6s*s234)+2*Abs(p1s*p3s*p4s*p6s*s234)+2*Abs(p2s*p3s*p4s*p6s*s
     -   234)+2*Abs(p2s*p4s**2*p6s*s234)+2*Abs(p2s*p4s2*p6s*s234)+Abs
     -   (p2s**2*p5s*p6s*s234)+Abs(p2s2*p5s*p6s*s234)+4*Abs(p1s*p3s*p
     -   5s*p6s*s234)+Abs(p2s*p3s*p5s*p6s*s234)+Abs(p3s**2*p5s*p6s*s2
     -   34)+2*Abs(p2s*p3s*p6s**2*s234)+2*Abs(p3s2*p6s**2*s234)+2*Abs
     -   (p2s*p4s*p6s**2*s234)+2*Abs(p3s*p4s*p6s**2*s234)+2*Abs(p1s*p
     -   3s*p4s*s12*s234)+Abs(p3s**2*p4s*s12*s234)+Abs(p3s2*p4s*s12*s
     -   234)+2*Abs(p1s*p3s*p5s*s12*s234)+Abs(p3s**2*p5s*s12*s234)+Ab
     -   s(p3s2*p5s*s12*s234)+2*Abs(p2s*p4s*p6s*s12*s234)+2*Abs(p3s*p
     -   4s*p6s*s12*s234)+2*Abs(p2s*p5s*p6s*s12*s234)+2*Abs(p3s*p5s*p
     -   6s*s12*s234)+Abs(p2s**2*p4s*s123*s234)+Abs(p2s2*p4s*s123*s23
     -   4)+Abs(p2s**2*p5s*s123*s234)+Abs(p2s2*p5s*s123*s234)+2*Abs(p
     -   1s*p3s*p5s*s123*s234)+Abs(p2s*p3s*p5s*s123*s234)+2*Abs(p2s*p
     -   3s*p6s*s123*s234)+Abs(p3s**2*p6s*s123*s234)+Abs(p3s2*p6s*s12
     -   3*s234)+2*Abs(p2s*p4s*p6s*s123*s234)+2*Abs(p2s*p5s*p6s*s123*
     -   s234)+2*Abs(p3s*p5s*p6s*s123*s234)+2*Abs(p2s*p4s*s12*s123*s2
     -   34)+2*Abs(p2s*p5s*s12*s123*s234)+4*Abs(p3s*p5s*s12*s123*s234
     -   )+2*Abs(p2s*p5s*s123**2*s234)+Abs(p3s2*p6s*s16*s234)+Abs(p3s
     -   *p4s*p6s*s16*s234)+Abs(p3s*p4s*s12*s16*s234)+2*Abs(p3s*p5s*s
     -   12*s16*s234)+2*Abs(p3s*p6s*s12*s16*s234)+2*Abs(p4s*p6s*s12*s
     -   16*s234)+2*Abs(p4s*s12**2*s16*s234)+2*Abs(p5s*s12**2*s16*s23
     -   4)+2*Abs(p3s*p5s*s123*s16*s234)+2*Abs(p3s*p6s*s123*s16*s234)
     -   +2*Abs(p4s*p6s*s123*s16*s234)+4*Abs(p3s*s12*s123*s16*s234)+2
     -   *Abs(p4s*s12*s123*s16*s234)+4*Abs(p5s*s12*s123*s16*s234)+2*A
     -   bs(p5s*s123**2*s16*s234)+Abs(p3s**2*p6s*s23*s234)+Abs(p3s2*p
     -   6s*s23*s234)+Abs(p3s*p5s*p6s*s23*s234)+2*Abs(p3s*p6s**2*s23*
     -   s234)+2*Abs(p4s*p6s**2*s23*s234)+Abs(p3s*p5s*s12*s23*s234)+2
     -   *Abs(p3s*p6s*s12*s23*s234)+4*Abs(p4s*p6s*s12*s23*s234)+2*Abs
     -   (p5s*p6s*s12*s23*s234)+2*Abs(p4s*s12**2*s23*s234)+2*Abs(p5s*
     -   s12**2*s23*s234)+2*Abs(p5s*p6s*s123*s23*s234)+2*Abs(p5s*s12*
     -   s123*s23*s234)+Abs(p4s*p6s*s23**2*s234)+Abs(p5s*p6s*s23**2*s
     -   234)+Abs(p4s*s123*s23**2*s234)+Abs(p5s*s123*s23**2*s234)+Abs
     -   (p3s2*p6s*s234**2)+2*Abs(p3s*p6s**2*s234**2)+2*Abs(p3s*p6s*s
     -   12*s234**2)+2*Abs(p3s*p6s*s123*s234**2)+2*Abs(p3s*s12*s123*s
     -   234**2)+Abs(p1s*p4s**2*s23s)+2*Abs(p1s*p4s*p5s*s23s)+Abs(p1s
     -   *p5s**2*s23s)+Abs(p3s*p4s*p6s*s23s)+Abs(p3s*p5s*p6s*s23s)+Ab
     -   s(p4s**2*s12*s23s)+2*Abs(p4s*p5s*s12*s23s)+Abs(p5s**2*s12*s2
     -   3s)+2*Abs(p3s*p5s*s123*s23s)+2*Abs(p4s*p6s*s16*s23s)+Abs(p4s
     -   *s123*s16*s23s)+Abs(p5s*s123*s16*s23s)+Abs(p4s*p6s*s23*s23s)
     -   +Abs(p5s*p6s*s23*s23s)+2*Abs(p5s*s123*s23*s23s)+Abs(p4s*p6s*
     -   s234*s23s)+Abs(p5s*p6s*s234*s23s)+Abs(p4s*s123*s234*s23s)+Ab
     -   s(p5s*s123*s234*s23s)+2*Abs(p1s**2*p3s*p4s*s34)+Abs(p1s*p3s*
     -   *2*p4s*s34)+Abs(p1s*p3s2*p4s*s34)+2*Abs(p1s**2*p3s*p5s*s34)+
     -   Abs(p1s*p3s**2*p5s*s34)+Abs(p1s*p3s2*p5s*s34)+2*Abs(p1s*p2s*
     -   p4s*p6s*s34)+4*Abs(p1s*p3s*p4s*p6s*s34)+2*Abs(p1s*p2s*p5s*p6
     -   s*s34)+2*Abs(p1s*p3s*p5s*p6s*s34)+2*Abs(p2s*p4s*p6s**2*s34)+
     -   2*Abs(p3s*p4s*p6s**2*s34)+2*Abs(p1s*p2s*p4s*s123*s34)+2*Abs(
     -   p1s*p2s*p5s*s123*s34)+2*Abs(p1s*p3s*p5s*s123*s34)+Abs(p3s**2
     -   *p6s*s123*s34)+Abs(p3s2*p6s*s123*s34)+2*Abs(p2s*p4s*p6s*s123
     -   *s34)+2*Abs(p2s*p5s*p6s*s123*s34)+2*Abs(p3s*p5s*p6s*s123*s34
     -   )+2*Abs(p2s*p5s*s123**2*s34)+2*Abs(p1s*p3s*p4s*s16*s34)+Abs(
     -   p1s*p3s*p5s*s16*s34)+4*Abs(p1s*p3s*p6s*s16*s34)+4*Abs(p1s*p4
     -   s*p6s*s16*s34)+Abs(p2s*p4s*p6s*s16*s34)+2*Abs(p3s*p4s*p6s*s1
     -   6*s34)+2*Abs(p1s*p4s*s12*s16*s34)+2*Abs(p1s*p5s*s12*s16*s34)
     -   +2*Abs(p4s*p6s*s12*s16*s34)+2*Abs(p1s*p3s*s123*s16*s34)+2*Ab
     -   s(p1s*p4s*s123*s16*s34)+Abs(p2s*p4s*s123*s16*s34)+2*Abs(p1s*
     -   p5s*s123*s16*s34)+2*Abs(p2s*p5s*s123*s16*s34)+Abs(p3s*p5s*s1
     -   23*s16*s34)+2*Abs(p2s*p6s*s123*s16*s34)+2*Abs(p3s*p6s*s123*s
     -   16*s34)+2*Abs(p4s*p6s*s123*s16*s34)+4*Abs(p4s*s12*s123*s16*s
     -   34)+2*Abs(p5s*s12*s123*s16*s34)+2*Abs(p2s*s123**2*s16*s34)+2
     -   *Abs(p5s*s123**2*s16*s34)+Abs(p4s*s12*s16**2*s34)+Abs(p3s*s1
     -   23*s16**2*s34)+Abs(p4s*s123*s16**2*s34)+2*Abs(s12*s123*s16**
     -   2*s34)+2*Abs(s123**2*s16**2*s34)+Abs(p1s*p3s*p5s*s23*s34)+2*
     -   Abs(p1s*p3s*p6s*s23*s34)+Abs(p3s**2*p6s*s23*s34)+Abs(p3s2*p6
     -   s*s23*s34)+2*Abs(p1s*p4s*p6s*s23*s34)+2*Abs(p2s*p4s*p6s*s23*
     -   s34)+2*Abs(p1s*p5s*p6s*s23*s34)+Abs(p2s*p5s*p6s*s23*s34)+Abs
     -   (p3s*p5s*p6s*s23*s34)+2*Abs(p2s*p6s**2*s23*s34)+2*Abs(p3s*p6
     -   s**2*s23*s34)+2*Abs(p4s*p6s**2*s23*s34)+2*Abs(p1s*p4s*s12*s2
     -   3*s34)+2*Abs(p1s*p5s*s12*s23*s34)+2*Abs(p4s*p6s*s12*s23*s34)
     -   +4*Abs(p5s*p6s*s12*s23*s34)+4*Abs(p1s*p5s*s123*s23*s34)+Abs(
     -   p2s*p5s*s123*s23*s34)+2*Abs(p2s*p6s*s123*s23*s34)+2*Abs(p5s*
     -   p6s*s123*s23*s34)+2*Abs(p5s*s12*s123*s23*s34)+Abs(p3s*p6s*s1
     -   6*s23*s34)+Abs(p4s*p6s*s16*s23*s34)+Abs(p4s*s12*s16*s23*s34)
     -   +Abs(p5s*s12*s16*s23*s34)+2*Abs(p6s*s12*s16*s23*s34)+Abs(p5s
     -   *s123*s16*s23*s34)+4*Abs(p6s*s123*s16*s23*s34)+2*Abs(s12*s12
     -   3*s16*s23*s34)+Abs(p4s*p6s*s23**2*s34)+2*Abs(p6s**2*s23**2*s
     -   34)+Abs(p5s*s12*s23**2*s34)+2*Abs(p6s*s12*s23**2*s34)+Abs(p4
     -   s*s123*s23**2*s34)+Abs(p5s*s123*s23**2*s34)+2*Abs(p1s*p3s*p6
     -   s*s234*s34)+Abs(p3s**2*p6s*s234*s34)+Abs(p3s2*p6s*s234*s34)+
     -   2*Abs(p2s*p6s**2*s234*s34)+2*Abs(p3s*p6s**2*s234*s34)+2*Abs(
     -   p1s*p3s*s123*s234*s34)+4*Abs(p2s*p6s*s123*s234*s34)+2*Abs(p3
     -   s*p6s*s123*s234*s34)+2*Abs(p2s*s123**2*s234*s34)+Abs(p3s*p6s
     -   *s16*s234*s34)+2*Abs(p6s*s12*s16*s234*s34)+Abs(p3s*s123*s16*
     -   s234*s34)+2*Abs(p6s*s123*s16*s234*s34)+2*Abs(s12*s123*s16*s2
     -   34*s34)+2*Abs(s123**2*s16*s234*s34)+2*Abs(p3s*p6s*s23*s234*s
     -   34)+2*Abs(p6s**2*s23*s234*s34)+2*Abs(p6s*s12*s23*s234*s34)+2
     -   *Abs(p6s*s123*s23*s234*s34)+2*Abs(s12*s123*s23*s234*s34)+Abs
     -   (p4s*p6s*s23s*s34)+Abs(p5s*p6s*s23s*s34)+Abs(p4s*s123*s23s*s
     -   34)+Abs(p5s*s123*s23s*s34)+2*Abs(p1s*p6s*s16*s34**2)+2*Abs(p
     -   1s*s123*s16*s34**2)+2*Abs(p6s*s123*s16*s34**2)+2*Abs(s123**2
     -   *s16*s34**2)+Abs(s123*s16**2*s34**2)+2*Abs(p1s*p6s*s23*s34**
     -   2)+2*Abs(p6s**2*s23*s34**2)+2*Abs(p1s*s123*s23*s34**2)+2*Abs
     -   (p6s*s123*s23*s34**2)+Abs(p6s*s16*s23*s34**2)+Abs(s123*s16*s
     -   23*s34**2)+Abs(p6s*s23**2*s34**2)+2*Abs(p1s**2*p3s*p4s*s345)
     -   +Abs(p1s*p3s**2*p4s*s345)+Abs(p1s*p3s2*p4s*s345)+2*Abs(p1s**
     -   2*p4s**2*s345)+Abs(p1s*p2s*p4s**2*s345)+2*Abs(p1s**2*p3s*p5s
     -   *s345)+Abs(p1s*p3s**2*p5s*s345)+Abs(p1s*p3s2*p5s*s345)+2*Abs
     -   (p1s**2*p4s*p5s*s345)+Abs(p1s*p2s*p4s*p5s*s345)+2*Abs(p1s*p3
     -   s*p4s*p5s*s345)+4*Abs(p1s*p2s*p4s*p6s*s345)+2*Abs(p1s*p3s*p4
     -   s*p6s*s345)+2*Abs(p1s*p4s**2*p6s*s345)+Abs(p2s*p4s**2*p6s*s3
     -   45)+2*Abs(p1s*p4s**2*s12*s345)+2*Abs(p1s*p4s*p5s*s12*s345)+2
     -   *Abs(p4s**2*p6s*s12*s345)+2*Abs(p1s*p2s*p4s*s123*s345)+2*Abs
     -   (p1s*p2s*p5s*s123*s345)+2*Abs(p1s*p3s*p5s*s123*s345)+2*Abs(p
     -   1s*p4s*p5s*s123*s345)+Abs(p2s*p4s*p5s*s123*s345)+Abs(p3s**2*
     -   p6s*s123*s345)+Abs(p3s2*p6s*s123*s345)+2*Abs(p2s*p4s*p6s*s12
     -   3*s345)+2*Abs(p4s*p5s*s12*s123*s345)+2*Abs(p2s*p5s*s123**2*s
     -   345)+Abs(p1s*p3s*p4s*s16*s345)+Abs(p1s*p4s**2*s16*s345)+4*Ab
     -   s(p1s*p4s*s12*s16*s345)+Abs(p4s**2*s12*s16*s345)+2*Abs(p1s*p
     -   3s*s123*s16*s345)+2*Abs(p1s*p4s*s123*s16*s345)+2*Abs(p2s*p4s
     -   *s123*s16*s345)+2*Abs(p4s*s12*s123*s16*s345)+2*Abs(p2s*s123*
     -   *2*s16*s345)+Abs(p1s*p3s*p5s*s23*s345)+Abs(p1s*p4s*p5s*s23*s
     -   345)+2*Abs(p1s*p3s*p6s*s23*s345)+Abs(p3s**2*p6s*s23*s345)+Ab
     -   s(p3s2*p6s*s23*s345)+2*Abs(p1s*p4s*p6s*s23*s345)+Abs(p2s*p4s
     -   *p6s*s23*s345)+2*Abs(p1s*p4s*s12*s23*s345)+2*Abs(p1s*p5s*s12
     -   *s23*s345)+Abs(p4s*p5s*s12*s23*s345)+4*Abs(p4s*p6s*s12*s23*s
     -   345)+4*Abs(p1s*p5s*s123*s23*s345)+Abs(p2s*p5s*s123*s23*s345)
     -   +2*Abs(p2s*p6s*s123*s23*s345)+2*Abs(p5s*s12*s123*s23*s345)+A
     -   bs(p4s*s12*s16*s23*s345)+2*Abs(s12*s123*s16*s23*s345)+2*Abs(
     -   p4s*p6s*s23**2*s345)+Abs(p5s*s12*s23**2*s345)+2*Abs(p6s*s12*
     -   s23**2*s345)+Abs(p4s*s123*s23**2*s345)+Abs(p5s*s123*s23**2*s
     -   345)+Abs(p1s*p3s*p4s*s234*s345)+Abs(p1s*p3s*p5s*s234*s345)+2
     -   *Abs(p1s*p3s*p6s*s234*s345)+Abs(p3s**2*p6s*s234*s345)+Abs(p3
     -   s2*p6s*s234*s345)+2*Abs(p1s*p4s*p6s*s234*s345)+Abs(p2s*p4s*p
     -   6s*s234*s345)+Abs(p3s*p4s*p6s*s234*s345)+2*Abs(p1s*p4s*s12*s
     -   234*s345)+2*Abs(p1s*p5s*s12*s234*s345)+4*Abs(p4s*p6s*s12*s23
     -   4*s345)+2*Abs(p1s*p3s*s123*s234*s345)+4*Abs(p1s*p4s*s123*s23
     -   4*s345)+Abs(p2s*p4s*s123*s234*s345)+2*Abs(p1s*p5s*s123*s234*
     -   s345)+2*Abs(p2s*p5s*s123*s234*s345)+Abs(p3s*p5s*s123*s234*s3
     -   45)+2*Abs(p2s*p6s*s123*s234*s345)+4*Abs(p3s*p6s*s123*s234*s3
     -   45)+2*Abs(p4s*p6s*s123*s234*s345)+2*Abs(p4s*s12*s123*s234*s3
     -   45)+2*Abs(p5s*s12*s123*s234*s345)+2*Abs(p2s*s123**2*s234*s34
     -   5)+2*Abs(p5s*s123**2*s234*s345)+Abs(p4s*s12*s16*s234*s345)+A
     -   bs(p3s*s123*s16*s234*s345)+Abs(p4s*s123*s16*s234*s345)+2*Abs
     -   (s12*s123*s16*s234*s345)+2*Abs(s123**2*s16*s234*s345)+Abs(p3
     -   s*p6s*s23*s234*s345)+2*Abs(p4s*p6s*s23*s234*s345)+2*Abs(p4s*
     -   s12*s23*s234*s345)+Abs(p5s*s12*s23*s234*s345)+4*Abs(p6s*s12*
     -   s23*s234*s345)+Abs(p5s*s123*s23*s234*s345)+2*Abs(p6s*s123*s2
     -   3*s234*s345)+2*Abs(s12*s123*s23*s234*s345)+Abs(p3s*p6s*s234*
     -   *2*s345)+2*Abs(p6s*s12*s234**2*s345)+Abs(p3s*s123*s234**2*s3
     -   45)+2*Abs(p6s*s123*s234**2*s345)+2*Abs(s12*s123*s234**2*s345
     -   )+2*Abs(s123**2*s234**2*s345)+2*Abs(p4s*p6s*s23s*s345)+Abs(p
     -   4s*s123*s23s*s345)+Abs(p5s*s123*s23s*s345)+2*Abs(p1s**2*p4s*
     -   s34*s345)+2*Abs(p1s**2*p5s*s34*s345)+2*Abs(p1s*p4s*p6s*s34*s
     -   345)+2*Abs(p1s*p4s*s123*s34*s345)+4*Abs(p1s*p5s*s123*s34*s34
     -   5)+2*Abs(p4s*p6s*s123*s34*s345)+2*Abs(p5s*s123**2*s34*s345)+
     -   Abs(p1s*p4s*s16*s34*s345)+2*Abs(p1s*s123*s16*s34*s345)+Abs(p
     -   4s*s123*s16*s34*s345)+2*Abs(s123**2*s16*s34*s345)+Abs(p1s*p4
     -   s*s23*s34*s345)+2*Abs(p1s*p5s*s23*s34*s345)+2*Abs(p1s*p6s*s2
     -   3*s34*s345)+Abs(p4s*p6s*s23*s34*s345)+4*Abs(p1s*s123*s23*s34
     -   *s345)+2*Abs(p5s*s123*s23*s34*s345)+2*Abs(p6s*s123*s23*s34*s
     -   345)+Abs(s123*s16*s23*s34*s345)+Abs(p6s*s23**2*s34*s345)+2*A
     -   bs(p1s*p6s*s234*s34*s345)+2*Abs(p1s*s123*s234*s34*s345)+2*Ab
     -   s(p6s*s123*s234*s34*s345)+2*Abs(s123**2*s234*s34*s345)+2*Abs
     -   (s123*s16*s234*s34*s345)+Abs(p6s*s23*s234*s34*s345)+Abs(s123
     -   *s23*s234*s34*s345)+2*Abs(p1s**2*p4s*s345**2)+Abs(p1s*p4s**2
     -   *s345**2)+2*Abs(p1s*p4s*s123*s345**2)+Abs(p1s*p4s*s23*s345**
     -   2)+2*Abs(p1s*s123*s23*s345**2)+Abs(p1s*p4s*s234*s345**2)+2*A
     -   bs(p1s*s123*s234*s345**2)+Abs(p4s*s123*s234*s345**2)+2*Abs(s
     -   123**2*s234*s345**2)+Abs(s123*s23*s234*s345**2)+Abs(s123*s23
     -   4**2*s345**2)+2*Abs(p1s*p2s**2*p4s*s45)+2*Abs(p1s*p2s2*p4s*s
     -   45)+2*Abs(p1s*p2s**2*p5s*s45)+2*Abs(p1s*p2s2*p5s*s45)+2*Abs(
     -   p1s**2*p3s*p5s*s45)+Abs(p1s*p2s*p3s*p5s*s45)+2*Abs(p1s**2*p4
     -   s*p5s*s45)+2*Abs(p1s*p2s*p4s*p5s*s45)+2*Abs(p1s*p2s*p3s*p6s*
     -   s45)+Abs(p1s*p3s**2*p6s*s45)+Abs(p1s*p3s2*p6s*s45)+2*Abs(p1s
     -   *p2s*p4s*p6s*s45)+Abs(p2s2*p4s*p6s*s45)+4*Abs(p1s*p2s*p4s*s1
     -   2*s45)+4*Abs(p1s*p2s*p5s*s12*s45)+2*Abs(p1s*p3s*p5s*s12*s45)
     -   +4*Abs(p1s*p4s*p5s*s12*s45)+2*Abs(p2s*p4s*p5s*s12*s45)+Abs(p
     -   3s**2*p6s*s12*s45)+Abs(p3s2*p6s*s12*s45)+2*Abs(p2s*p4s*p6s*s
     -   12*s45)+2*Abs(p4s*p5s*s12**2*s45)+2*Abs(p1s*p2s*p5s*s123*s45
     -   )+Abs(p2s2*p5s*s123*s45)+2*Abs(p2s**2*p6s*s123*s45)+2*Abs(p2
     -   s*p5s*s12*s123*s45)+2*Abs(p1s*p3s*s12*s16*s45)+Abs(p3s**2*s1
     -   2*s16*s45)+Abs(p3s2*s12*s16*s45)+2*Abs(p1s*p4s*s12*s16*s45)+
     -   Abs(p2s*p4s*s12*s16*s45)+2*Abs(p4s*s12**2*s16*s45)+Abs(p2s**
     -   2*s123*s16*s45)+Abs(p2s2*s123*s16*s45)+2*Abs(p2s*s12*s123*s1
     -   6*s45)+Abs(p2s**2*p6s*s23*s45)+Abs(p2s2*p6s*s23*s45)+2*Abs(p
     -   1s*p5s*s12*s23*s45)+Abs(p2s*p5s*s12*s23*s45)+2*Abs(p2s*p6s*s
     -   12*s23*s45)+2*Abs(p5s*s12**2*s23*s45)+2*Abs(s12**2*s16*s23*s
     -   45)+2*Abs(p1s*p4s*s23**2*s45)+2*Abs(p1s*p5s*s23**2*s45)+Abs(
     -   p3s*p6s*s23**2*s45)+2*Abs(p4s*s12*s23**2*s45)+2*Abs(p5s*s12*
     -   s23**2*s45)+Abs(s123*s16*s23**2*s45)+Abs(p6s*s23**3*s45)+Abs
     -   (p1s*p3s*p5s*s234*s45)+Abs(p2s**2*p6s*s234*s45)+Abs(p2s2*p6s
     -   *s234*s45)+2*Abs(p1s*p3s*p6s*s234*s45)+Abs(p2s*p3s*p6s*s234*
     -   s45)+2*Abs(p1s*p4s*p6s*s234*s45)+Abs(p2s*p4s*p6s*s234*s45)+2
     -   *Abs(p1s*p3s*s12*s234*s45)+Abs(p3s**2*s12*s234*s45)+Abs(p3s2
     -   *s12*s234*s45)+2*Abs(p1s*p4s*s12*s234*s45)+Abs(p2s*p4s*s12*s
     -   234*s45)+2*Abs(p1s*p5s*s12*s234*s45)+Abs(p2s*p5s*s12*s234*s4
     -   5)+Abs(p3s*p5s*s12*s234*s45)+2*Abs(p2s*p6s*s12*s234*s45)+4*A
     -   bs(p3s*p6s*s12*s234*s45)+2*Abs(p4s*p6s*s12*s234*s45)+2*Abs(p
     -   4s*s12**2*s234*s45)+2*Abs(p5s*s12**2*s234*s45)+Abs(p2s**2*s1
     -   23*s234*s45)+Abs(p2s2*s123*s234*s45)+2*Abs(p1s*p5s*s123*s234
     -   *s45)+Abs(p2s*p5s*s123*s234*s45)+4*Abs(p2s*p6s*s123*s234*s45
     -   )+2*Abs(p2s*s12*s123*s234*s45)+2*Abs(p5s*s12*s123*s234*s45)+
     -   Abs(p3s*s12*s16*s234*s45)+2*Abs(p4s*s12*s16*s234*s45)+2*Abs(
     -   s12**2*s16*s234*s45)+2*Abs(s12*s123*s16*s234*s45)+2*Abs(p5s*
     -   s12*s23*s234*s45)+2*Abs(p6s*s12*s23*s234*s45)+2*Abs(s12**2*s
     -   23*s234*s45)+Abs(p6s*s23**2*s234*s45)+Abs(s123*s23**2*s234*s
     -   45)+Abs(p3s*p6s*s234**2*s45)+Abs(p3s*s12*s234**2*s45)+2*Abs(
     -   p6s*s12*s234**2*s45)+2*Abs(s12**2*s234**2*s45)+2*Abs(p6s*s12
     -   3*s234**2*s45)+2*Abs(s12*s123*s234**2*s45)+2*Abs(p1s*p4s*s23
     -   s*s45)+2*Abs(p1s*p5s*s23s*s45)+Abs(p3s*p6s*s23s*s45)+2*Abs(p
     -   4s*s12*s23s*s45)+2*Abs(p5s*s12*s23s*s45)+Abs(s123*s16*s23s*s
     -   45)+Abs(p6s*s23*s23s*s45)+Abs(p6s*s234*s23s*s45)+Abs(s123*s2
     -   34*s23s*s45)+2*Abs(p1s**2*p3s*s34*s45)+Abs(p1s*p3s**2*s34*s4
     -   5)+Abs(p1s*p3s2*s34*s45)+2*Abs(p1s**2*p4s*s34*s45)+Abs(p1s*p
     -   2s*p4s*s34*s45)+2*Abs(p1s**2*p5s*s34*s45)+Abs(p1s*p2s*p5s*s3
     -   4*s45)+2*Abs(p1s*p3s*p5s*s34*s45)+2*Abs(p1s*p2s*p6s*s34*s45)
     -   +2*Abs(p1s*p3s*p6s*s34*s45)+2*Abs(p1s*p4s*p6s*s34*s45)+Abs(p
     -   2s*p4s*p6s*s34*s45)+2*Abs(p1s*p4s*s12*s34*s45)+2*Abs(p1s*p5s
     -   *s12*s34*s45)+2*Abs(p4s*p6s*s12*s34*s45)+2*Abs(p1s*p2s*s123*
     -   s34*s45)+2*Abs(p1s*p5s*s123*s34*s45)+Abs(p2s*p5s*s123*s34*s4
     -   5)+4*Abs(p2s*p6s*s123*s34*s45)+2*Abs(p5s*s12*s123*s34*s45)+A
     -   bs(p1s*p3s*s16*s34*s45)+Abs(p1s*p4s*s16*s34*s45)+2*Abs(p1s*s
     -   12*s16*s34*s45)+Abs(p4s*s12*s16*s34*s45)+4*Abs(p1s*s123*s16*
     -   s34*s45)+Abs(p2s*s123*s16*s34*s45)+2*Abs(s12*s123*s16*s34*s4
     -   5)+Abs(p1s*p5s*s23*s34*s45)+4*Abs(p1s*p6s*s23*s34*s45)+Abs(p
     -   2s*p6s*s23*s34*s45)+2*Abs(p1s*s12*s23*s34*s45)+Abs(p5s*s12*s
     -   23*s34*s45)+2*Abs(p6s*s12*s23*s34*s45)+2*Abs(s12*s16*s23*s34
     -   *s45)+Abs(p6s*s23**2*s34*s45)+Abs(s123*s23**2*s34*s45)+Abs(p
     -   1s*p3s*s234*s34*s45)+2*Abs(p1s*p6s*s234*s34*s45)+2*Abs(p2s*p
     -   6s*s234*s34*s45)+Abs(p3s*p6s*s234*s34*s45)+4*Abs(p1s*s12*s23
     -   4*s34*s45)+2*Abs(p6s*s12*s234*s34*s45)+2*Abs(p1s*s123*s234*s
     -   34*s45)+2*Abs(p2s*s123*s234*s34*s45)+4*Abs(p6s*s123*s234*s34
     -   *s45)+2*Abs(s12*s123*s234*s34*s45)+Abs(s12*s16*s234*s34*s45)
     -   +Abs(s123*s16*s234*s34*s45)+Abs(p6s*s23*s234*s34*s45)+Abs(s1
     -   2*s23*s234*s34*s45)+Abs(p6s*s23s*s34*s45)+Abs(s123*s23s*s34*
     -   s45)+2*Abs(p1s**2*s34**2*s45)+2*Abs(p1s*p6s*s34**2*s45)+2*Ab
     -   s(p1s*s123*s34**2*s45)+2*Abs(p6s*s123*s34**2*s45)+Abs(p1s*s1
     -   6*s34**2*s45)+Abs(s123*s16*s34**2*s45)+Abs(p1s*s23*s34**2*s4
     -   5)+Abs(p6s*s23*s34**2*s45)+2*Abs(p1s**2*p3s*s345*s45)+Abs(p1
     -   s*p3s**2*s345*s45)+Abs(p1s*p3s2*s345*s45)+2*Abs(p1s**2*p4s*s
     -   345*s45)+Abs(p1s*p2s*p4s*s345*s45)+2*Abs(p1s*p4s*s12*s345*s4
     -   5)+2*Abs(p1s*p2s*s123*s345*s45)+2*Abs(p1s*s12*s23*s345*s45)+
     -   Abs(s123*s23**2*s345*s45)+2*Abs(p1s*p3s*s234*s345*s45)+Abs(p
     -   1s*p4s*s234*s345*s45)+2*Abs(p1s*s12*s234*s345*s45)+Abs(p4s*s
     -   12*s234*s345*s45)+2*Abs(p1s*s123*s234*s345*s45)+Abs(p2s*s123
     -   *s234*s345*s45)+4*Abs(s12*s123*s234*s345*s45)+Abs(s12*s23*s2
     -   34*s345*s45)+Abs(s12*s234**2*s345*s45)+Abs(s123*s234**2*s345
     -   *s45)+Abs(s123*s23s*s345*s45)+2*Abs(p1s**2*s34*s345*s45)+2*A
     -   bs(p1s*p4s*s34*s345*s45)+2*Abs(p1s*s123*s34*s345*s45)+Abs(p1
     -   s*s23*s34*s345*s45)+Abs(p1s*s234*s34*s345*s45)+Abs(s123*s234
     -   *s34*s345*s45)+Abs(p1s*p2s**2*s45**2)+Abs(p1s*p2s2*s45**2)+2
     -   *Abs(p1s*p2s*s12*s45**2)+Abs(p1s*s23**2*s45**2)+Abs(s12*s23*
     -   *2*s45**2)+2*Abs(p1s*s12*s234*s45**2)+Abs(p2s*s12*s234*s45**
     -   2)+2*Abs(s12**2*s234*s45**2)+Abs(s12*s234**2*s45**2)+Abs(p1s
     -   *s23s*s45**2)+Abs(s12*s23s*s45**2)+2*Abs(p1s**2*s34*s45**2)+
     -   Abs(p1s*p2s*s34*s45**2)+2*Abs(p1s*s12*s34*s45**2)+Abs(p1s*s2
     -   34*s34*s45**2)+Abs(s12*s234*s34*s45**2)+Abs(p1s*s34**2*s45**
     -   2)+2*Abs(p1s*p2s*p3s*p4s*s56)+Abs(p1s*p3s**2*p4s*s56)+Abs(p1
     -   s*p3s2*p4s*s56)+2*Abs(p1s*p2s*p3s*p5s*s56)+Abs(p1s*p3s**2*p5
     -   s*s56)+Abs(p1s*p3s2*p5s*s56)+2*Abs(p2s**2*p4s*p6s*s56)+4*Abs
     -   (p2s2*p4s*p6s*s56)+2*Abs(p2s*p3s*p4s*p6s*s56)+4*Abs(p2s*p4s*
     -   *2*p6s*s56)+4*Abs(p2s*p4s2*p6s*s56)+2*Abs(p2s**2*p5s*p6s*s56
     -   )+4*Abs(p2s*p3s*p5s*p6s*s56)+2*Abs(p3s**2*p5s*p6s*s56)+Abs(p
     -   3s**2*p4s*s12*s56)+Abs(p3s2*p4s*s12*s56)+Abs(p3s**2*p5s*s12*
     -   s56)+Abs(p3s2*p5s*s12*s56)+2*Abs(p2s**2*p4s*s123*s56)+2*Abs(
     -   p2s**2*p5s*s123*s56)+2*Abs(p2s*p3s*p5s*s123*s56)+2*Abs(p1s*p
     -   3s**2*s16*s56)+Abs(p2s**2*p4s*s16*s56)+Abs(p2s2*p4s*s16*s56)
     -   +2*Abs(p1s*p3s*p4s*s16*s56)+Abs(p2s*p3s*p4s*s16*s56)+2*Abs(p
     -   2s*p4s**2*s16*s56)+2*Abs(p2s*p4s2*s16*s56)+Abs(p2s**2*p5s*s1
     -   6*s56)+Abs(p2s2*p5s*s16*s56)+4*Abs(p1s*p3s*p5s*s16*s56)+Abs(
     -   p2s*p3s*p5s*s16*s56)+Abs(p3s**2*p5s*s16*s56)+2*Abs(p2s*p3s*p
     -   6s*s16*s56)+Abs(p3s**2*p6s*s16*s56)+Abs(p3s2*p6s*s16*s56)+2*
     -   Abs(p2s*p4s*p6s*s16*s56)+2*Abs(p3s*p4s*p6s*s16*s56)+2*Abs(p2
     -   s*p4s*s12*s16*s56)+4*Abs(p3s*p4s*s12*s16*s56)+2*Abs(p2s*p5s*
     -   s12*s16*s56)+2*Abs(p3s*p5s*s12*s16*s56)+2*Abs(p2s*p3s*s123*s
     -   16*s56)+Abs(p3s**2*s123*s16*s56)+Abs(p3s2*s123*s16*s56)+4*Ab
     -   s(p2s*p4s*s123*s16*s56)+2*Abs(p2s*p5s*s123*s16*s56)+2*Abs(p3
     -   s*p5s*s123*s16*s56)+Abs(p3s**2*s16**2*s56)+Abs(p3s*p4s*s16**
     -   2*s56)+2*Abs(p3s*s12*s16**2*s56)+2*Abs(p4s*s12*s16**2*s56)+2
     -   *Abs(p3s*s123*s16**2*s56)+2*Abs(p4s*s123*s16**2*s56)+Abs(p2s
     -   **2*p4s*s23*s56)+Abs(p2s2*p4s*s23*s56)+Abs(p2s**2*p5s*s23*s5
     -   6)+Abs(p2s2*p5s*s23*s56)+2*Abs(p1s*p3s*p5s*s23*s56)+2*Abs(p2
     -   s*p3s*p5s*s23*s56)+4*Abs(p2s*p3s*p6s*s23*s56)+2*Abs(p2s*p4s*
     -   p6s*s23*s56)+4*Abs(p2s*p5s*p6s*s23*s56)+4*Abs(p3s*p5s*p6s*s2
     -   3*s56)+2*Abs(p2s*p4s*s12*s23*s56)+2*Abs(p2s*p5s*s12*s23*s56)
     -   +2*Abs(p3s*p5s*s12*s23*s56)+2*Abs(p2s*p5s*s123*s23*s56)+Abs(
     -   p3s**2*s16*s23*s56)+Abs(p3s2*s16*s23*s56)+Abs(p3s*p5s*s16*s2
     -   3*s56)+2*Abs(p3s*p6s*s16*s23*s56)+2*Abs(p4s*p6s*s16*s23*s56)
     -   +2*Abs(p3s*s12*s16*s23*s56)+2*Abs(p4s*s12*s16*s23*s56)+2*Abs
     -   (p5s*s12*s16*s23*s56)+2*Abs(p5s*s123*s16*s23*s56)+Abs(p3s*p4
     -   s*s23**2*s56)+Abs(p3s*p5s*s23**2*s56)+2*Abs(p5s*p6s*s23**2*s
     -   56)+2*Abs(p5s*s12*s23**2*s56)+Abs(p4s*s16*s23**2*s56)+Abs(p5
     -   s*s16*s23**2*s56)+Abs(p4s*s23**3*s56)+Abs(p5s*s23**3*s56)+2*
     -   Abs(p1s*p3s**2*s234*s56)+2*Abs(p2s**2*p4s*s234*s56)+2*Abs(p2
     -   s2*p4s*s234*s56)+2*Abs(p2s*p4s**2*s234*s56)+2*Abs(p2s*p4s2*s
     -   234*s56)+2*Abs(p2s*p3s*p6s*s234*s56)+Abs(p3s**2*p6s*s234*s56
     -   )+3*Abs(p3s2*p6s*s234*s56)+2*Abs(p2s*p3s*s123*s234*s56)+Abs(
     -   p3s**2*s123*s234*s56)+Abs(p3s2*s123*s234*s56)+2*Abs(p3s**2*s
     -   16*s234*s56)+Abs(p3s2*s16*s234*s56)+4*Abs(p3s*p6s*s16*s234*s
     -   56)+2*Abs(p3s*s12*s16*s234*s56)+2*Abs(p3s*s123*s16*s234*s56)
     -   +Abs(p3s**2*s23*s234*s56)+Abs(p3s2*s23*s234*s56)+2*Abs(p3s*p
     -   6s*s23*s234*s56)+2*Abs(p3s*s12*s23*s234*s56)+Abs(p3s**2*s234
     -   **2*s56)+Abs(p3s2*s234**2*s56)+Abs(p3s*p4s*s23s*s56)+Abs(p3s
     -   *p5s*s23s*s56)+Abs(p4s*s16*s23s*s56)+Abs(p5s*s16*s23s*s56)+A
     -   bs(p4s*s23*s23s*s56)+Abs(p5s*s23*s23s*s56)+Abs(p3s**2*p6s*s3
     -   4*s56)+Abs(p3s2*p6s*s34*s56)+Abs(p3s**2*s123*s34*s56)+Abs(p3
     -   s2*s123*s34*s56)+2*Abs(p1s*p3s*s16*s34*s56)+2*Abs(p2s*p6s*s1
     -   6*s34*s56)+2*Abs(p3s*p6s*s16*s34*s56)+2*Abs(p2s*s123*s16*s34
     -   *s56)+4*Abs(p3s*s123*s16*s34*s56)+Abs(p3s*s16**2*s34*s56)+2*
     -   Abs(s12*s16**2*s34*s56)+2*Abs(s123*s16**2*s34*s56)+2*Abs(p1s
     -   *p3s*s23*s34*s56)+Abs(p3s**2*s23*s34*s56)+Abs(p3s2*s23*s34*s
     -   56)+2*Abs(p2s*p6s*s23*s34*s56)+2*Abs(p3s*p6s*s23*s34*s56)+2*
     -   Abs(p2s*s123*s23*s34*s56)+Abs(p3s*s16*s23*s34*s56)+2*Abs(p6s
     -   *s16*s23*s34*s56)+4*Abs(s12*s16*s23*s34*s56)+2*Abs(s123*s16*
     -   s23*s34*s56)+2*Abs(p6s*s23**2*s34*s56)+2*Abs(s12*s23**2*s34*
     -   s56)+Abs(p3s**2*s234*s34*s56)+Abs(p3s2*s234*s34*s56)+2*Abs(p
     -   1s*p2s*p4s*s345*s56)+2*Abs(p1s*p3s*p4s*s345*s56)+2*Abs(p1s*p
     -   2s*p5s*s345*s56)+2*Abs(p1s*p3s*p5s*s345*s56)+Abs(p3s**2*p6s*
     -   s345*s56)+Abs(p3s2*p6s*s345*s56)+2*Abs(p2s*p4s*p6s*s345*s56)
     -   +2*Abs(p3s*p4s*p6s*s345*s56)+Abs(p3s**2*s123*s345*s56)+Abs(p
     -   3s2*s123*s345*s56)+4*Abs(p2s*p4s*s123*s345*s56)+2*Abs(p2s*p5
     -   s*s123*s345*s56)+2*Abs(p3s*p5s*s123*s345*s56)+2*Abs(p1s*p3s*
     -   s16*s345*s56)+2*Abs(p1s*p4s*s16*s345*s56)+2*Abs(p2s*p4s*s16*
     -   s345*s56)+Abs(p3s*p4s*s16*s345*s56)+2*Abs(p4s*s12*s16*s345*s
     -   56)+4*Abs(p2s*s123*s16*s345*s56)+2*Abs(p3s*s123*s16*s345*s56
     -   )+4*Abs(p4s*s123*s16*s345*s56)+2*Abs(p1s*p3s*s23*s345*s56)+A
     -   bs(p3s**2*s23*s345*s56)+Abs(p3s2*s23*s345*s56)+4*Abs(p1s*p4s
     -   *s23*s345*s56)+Abs(p2s*p4s*s23*s345*s56)+2*Abs(p1s*p5s*s23*s
     -   345*s56)+Abs(p2s*p5s*s23*s345*s56)+Abs(p3s*p5s*s23*s345*s56)
     -   +2*Abs(p2s*p6s*s23*s345*s56)+2*Abs(p3s*p6s*s23*s345*s56)+2*A
     -   bs(p4s*p6s*s23*s345*s56)+2*Abs(p4s*s12*s23*s345*s56)+4*Abs(p
     -   5s*s12*s23*s345*s56)+2*Abs(p2s*s123*s23*s345*s56)+2*Abs(p5s*
     -   s123*s23*s345*s56)+2*Abs(p3s*s16*s23*s345*s56)+Abs(p4s*s16*s
     -   23*s345*s56)+2*Abs(s12*s16*s23*s345*s56)+2*Abs(s123*s16*s23*
     -   s345*s56)+Abs(p4s*s23**2*s345*s56)+2*Abs(p6s*s23**2*s345*s56
     -   )+2*Abs(s12*s23**2*s345*s56)+4*Abs(p1s*p3s*s234*s345*s56)+Ab
     -   s(p3s**2*s234*s345*s56)+Abs(p3s2*s234*s345*s56)+2*Abs(p2s*p6
     -   s*s234*s345*s56)+2*Abs(p3s*p6s*s234*s345*s56)+2*Abs(p2s*s123
     -   *s234*s345*s56)+2*Abs(p3s*s123*s234*s345*s56)+Abs(p3s*s16*s2
     -   34*s345*s56)+2*Abs(s12*s16*s234*s345*s56)+2*Abs(s123*s16*s23
     -   4*s345*s56)+Abs(p3s*s23*s234*s345*s56)+2*Abs(p6s*s23*s234*s3
     -   45*s56)+2*Abs(s12*s23*s234*s345*s56)+4*Abs(s123*s23*s234*s34
     -   5*s56)+Abs(p4s*s23s*s345*s56)+Abs(p5s*s23s*s345*s56)+2*Abs(p
     -   1s*s16*s34*s345*s56)+2*Abs(s123*s16*s34*s345*s56)+2*Abs(p1s*
     -   s23*s34*s345*s56)+4*Abs(p6s*s23*s34*s345*s56)+2*Abs(s123*s23
     -   *s34*s345*s56)+Abs(s16*s23*s34*s345*s56)+Abs(s23**2*s34*s345
     -   *s56)+2*Abs(p1s*p4s*s345**2*s56)+2*Abs(p4s*s123*s345**2*s56)
     -   +2*Abs(p1s*s23*s345**2*s56)+Abs(p4s*s23*s345**2*s56)+2*Abs(s
     -   123*s23*s345**2*s56)+Abs(s23**2*s345**2*s56)+2*Abs(p1s*s234*
     -   s345**2*s56)+2*Abs(s123*s234*s345**2*s56)+Abs(s23*s234*s345*
     -   *2*s56)+2*Abs(p1s*p2s*p3s*s45*s56)+Abs(p1s*p3s**2*s45*s56)+A
     -   bs(p1s*p3s2*s45*s56)+2*Abs(p1s*p2s*p4s*s45*s56)+Abs(p2s2*p4s
     -   *s45*s56)+2*Abs(p1s*p2s*p5s*s45*s56)+Abs(p2s2*p5s*s45*s56)+2
     -   *Abs(p1s*p3s*p5s*s45*s56)+Abs(p2s*p3s*p5s*s45*s56)+2*Abs(p2s
     -   **2*p6s*s45*s56)+2*Abs(p2s*p3s*p6s*s45*s56)+4*Abs(p2s*p4s*p6
     -   s*s45*s56)+Abs(p3s**2*s12*s45*s56)+Abs(p3s2*s12*s45*s56)+2*A
     -   bs(p2s*p4s*s12*s45*s56)+2*Abs(p2s*p5s*s12*s45*s56)+2*Abs(p3s
     -   *p5s*s12*s45*s56)+2*Abs(p2s**2*s123*s45*s56)+4*Abs(p2s*p5s*s
     -   123*s45*s56)+Abs(p2s**2*s16*s45*s56)+Abs(p2s2*s16*s45*s56)+2
     -   *Abs(p1s*p3s*s16*s45*s56)+2*Abs(p2s*p3s*s16*s45*s56)+2*Abs(p
     -   1s*p4s*s16*s45*s56)+Abs(p2s*p4s*s16*s45*s56)+2*Abs(p2s*s12*s
     -   16*s45*s56)+2*Abs(p3s*s12*s16*s45*s56)+2*Abs(p4s*s12*s16*s45
     -   *s56)+2*Abs(p2s*s123*s16*s45*s56)+Abs(p2s**2*s23*s45*s56)+Ab
     -   s(p2s2*s23*s45*s56)+2*Abs(p1s*p5s*s23*s45*s56)+Abs(p2s*p5s*s
     -   23*s45*s56)+2*Abs(p2s*p6s*s23*s45*s56)+2*Abs(p2s*s12*s23*s45
     -   *s56)+2*Abs(p5s*s12*s23*s45*s56)+4*Abs(s12*s16*s23*s45*s56)+
     -   Abs(p3s*s23**2*s45*s56)+Abs(s16*s23**2*s45*s56)+Abs(s23**3*s
     -   45*s56)+2*Abs(p2s**2*s234*s45*s56)+2*Abs(p2s2*s234*s45*s56)+
     -   4*Abs(p1s*p3s*s234*s45*s56)+Abs(p2s*p3s*s234*s45*s56)+2*Abs(
     -   p2s*p6s*s234*s45*s56)+2*Abs(p3s*p6s*s234*s45*s56)+4*Abs(p2s*
     -   s12*s234*s45*s56)+2*Abs(p3s*s12*s234*s45*s56)+2*Abs(p2s*s123
     -   *s234*s45*s56)+Abs(p3s*s16*s234*s45*s56)+2*Abs(s12*s16*s234*
     -   s45*s56)+2*Abs(s123*s16*s234*s45*s56)+2*Abs(p6s*s23*s234*s45
     -   *s56)+2*Abs(s12*s23*s234*s45*s56)+2*Abs(s23**2*s234*s45*s56)
     -   +Abs(p3s*s23s*s45*s56)+Abs(s16*s23s*s45*s56)+Abs(s23*s23s*s4
     -   5*s56)+2*Abs(s234*s23s*s45*s56)+4*Abs(p1s*p2s*s34*s45*s56)+2
     -   *Abs(p1s*p3s*s34*s45*s56)+2*Abs(p2s*p6s*s34*s45*s56)+2*Abs(p
     -   3s*p6s*s34*s45*s56)+2*Abs(p2s*s123*s34*s45*s56)+2*Abs(p1s*s1
     -   6*s34*s45*s56)+Abs(p2s*s16*s34*s45*s56)+Abs(p3s*s16*s34*s45*
     -   s56)+4*Abs(s12*s16*s34*s45*s56)+2*Abs(s123*s16*s34*s45*s56)+
     -   2*Abs(p1s*s23*s34*s45*s56)+Abs(p2s*s23*s34*s45*s56)+2*Abs(p6
     -   s*s23*s34*s45*s56)+4*Abs(s12*s23*s34*s45*s56)+2*Abs(s16*s23*
     -   s34*s45*s56)+2*Abs(s23**2*s34*s45*s56)+2*Abs(s23s*s34*s45*s5
     -   6)+2*Abs(p1s*p2s*s345*s45*s56)+4*Abs(p1s*p3s*s345*s45*s56)+2
     -   *Abs(p1s*p4s*s345*s45*s56)+Abs(p2s*p4s*s345*s45*s56)+2*Abs(p
     -   4s*s12*s345*s45*s56)+2*Abs(p2s*s123*s345*s45*s56)+2*Abs(p1s*
     -   s23*s345*s45*s56)+2*Abs(p2s*s23*s345*s45*s56)+2*Abs(s12*s23*
     -   s345*s45*s56)+Abs(s23**2*s345*s45*s56)+4*Abs(p1s*s234*s345*s
     -   45*s56)+Abs(p2s*s234*s345*s45*s56)+2*Abs(p3s*s234*s345*s45*s
     -   56)+2*Abs(s12*s234*s345*s45*s56)+2*Abs(s123*s234*s345*s45*s5
     -   6)+Abs(s23*s234*s345*s45*s56)+Abs(s23s*s345*s45*s56)+2*Abs(p
     -   1s*s34*s345*s45*s56)+2*Abs(s123*s34*s345*s45*s56)+Abs(s23*s3
     -   4*s345*s45*s56)+2*Abs(p1s*p2s*s45**2*s56)+Abs(p2s2*s45**2*s5
     -   6)+2*Abs(p2s*s12*s45**2*s56)+2*Abs(p1s*s234*s45**2*s56)+Abs(
     -   p2s*s234*s45**2*s56)+2*Abs(s12*s234*s45**2*s56)+2*Abs(p1s*s3
     -   4*s45**2*s56)+Abs(p2s*s34*s45**2*s56)+2*Abs(s12*s34*s45**2*s
     -   56)+2*Abs(p2s**2*p4s*s56**2)+2*Abs(p2s2*p4s*s56**2)+2*Abs(p2
     -   s*p4s**2*s56**2)+2*Abs(p2s*p4s2*s56**2)+2*Abs(p2s*p3s*s16*s5
     -   6**2)+Abs(p3s**2*s16*s56**2)+Abs(p3s2*s16*s56**2)+2*Abs(p3s*
     -   s16**2*s56**2)+2*Abs(p2s*p3s*s23*s56**2)+2*Abs(p3s*s16*s23*s
     -   56**2)+Abs(p3s**2*s234*s56**2)+Abs(p3s2*s234*s56**2)+Abs(p3s
     -   **2*s34*s56**2)+Abs(p3s2*s34*s56**2)+Abs(p3s**2*s345*s56**2)
     -   +Abs(p3s2*s345*s56**2)+2*Abs(p2s*s16*s345*s56**2)+2*Abs(p3s*
     -   s16*s345*s56**2)+2*Abs(p2s*s23*s345*s56**2)+2*Abs(p3s*s23*s3
     -   45*s56**2)+2*Abs(s16*s23*s345*s56**2)+2*Abs(s23**2*s345*s56*
     -   *2)+2*Abs(s23*s345**2*s56**2)+2*Abs(p2s**2*s45*s56**2)+2*Abs
     -   (p2s*p3s*s45*s56**2)+2*Abs(p2s*s16*s45*s56**2)+2*Abs(p3s*s16
     -   *s45*s56**2)+2*Abs(p2s*s23*s45*s56**2)+2*Abs(s16*s23*s45*s56
     -   **2)+2*Abs(p2s*s345*s45*s56**2)+2*Abs(p3s*s345*s45*s56**2)+2
     -   *Abs(s23*s345*s45*s56**2)+2*Abs(p2s*s45**2*s56**2)
       cx2=Abs(p1s*p2s**2*p4s**2)+Abs(p1s*p2s2*p4s**2)+2*Abs(p1s**2*p3
     -   s**2*p5s)+2*Abs(p1s*p2s**2*p4s*p5s)+2*Abs(p1s*p2s2*p4s*p5s)+
     -   2*Abs(p1s**2*p3s*p4s*p5s)+Abs(p1s*p2s*p3s*p4s*p5s)+Abs(p1s*p
     -   2s**2*p5s**2)+Abs(p1s*p2s2*p5s**2)+2*Abs(p1s**2*p3s*p5s**2)+
     -   Abs(p1s*p2s*p3s*p5s**2)+2*Abs(p1s*p2s*p3s*p4s*p6s)+2*Abs(p1s
     -   *p2s*p4s**2*p6s)+Abs(p2s**2*p4s**2*p6s)+2*Abs(p1s*p2s*p3s*p5
     -   s*p6s)+Abs(p1s*p3s**2*p5s*p6s)+2*Abs(p1s*p2s*p4s*p5s*p6s)+Ab
     -   s(p2s**2*p4s*p5s*p6s)+Abs(p1s*p3s*p4s*p5s*p6s)+2*Abs(p2s2*p4
     -   s*p6s**2)+Abs(p2s*p3s*p4s*p6s**2)+Abs(p2s*p4s**2*p6s**2)+2*A
     -   bs(p1s*p2s*p4s**2*s12)+4*Abs(p1s*p2s*p4s*p5s*s12)+Abs(p1s*p3
     -   s*p4s*p5s*s12)+2*Abs(p1s*p2s*p5s**2*s12)+Abs(p1s*p3s*p5s**2*
     -   s12)+Abs(p2s*p4s**2*p6s*s12)+Abs(p2s*p4s*p5s*p6s*s12)+4*Abs(
     -   p1s*p2s*p3s*p5s*s123)+2*Abs(p1s*p2s*p4s*p5s*s123)+Abs(p2s**2
     -   *p4s*p5s*s123)+2*Abs(p1s*p2s*p5s**2*s123)+Abs(p2s**2*p5s**2*
     -   s123)+Abs(p1s*p3s*p5s**2*s123)+2*Abs(p2s2*p4s*p6s*s123)+2*Ab
     -   s(p2s2*p5s*p6s*s123)+Abs(p2s*p3s*p5s*p6s*s123)+2*Abs(p2s*p4s
     -   *p5s*p6s*s123)+Abs(p2s*p4s*p5s*s12*s123)+Abs(p2s*p5s**2*s12*
     -   s123)+2*Abs(p2s2*p5s*s123**2)+Abs(p2s*p5s**2*s123**2)+Abs(p1
     -   s*p3s**2*p5s*s16)+Abs(p1s*p3s*p4s*p5s*s16)+2*Abs(p1s*p3s**2*
     -   p6s*s16)+2*Abs(p2s**2*p4s*p6s*s16)+2*Abs(p2s2*p4s*p6s*s16)+4
     -   *Abs(p1s*p3s*p4s*p6s*s16)+Abs(p2s*p3s*p4s*p6s*s16)+2*Abs(p1s
     -   *p4s**2*p6s*s16)+Abs(p2s*p4s**2*p6s*s16)+2*Abs(p1s*p3s*p4s*s
     -   12*s16)+2*Abs(p1s*p4s**2*s12*s16)+Abs(p2s*p4s**2*s12*s16)+2*
     -   Abs(p1s*p3s*p5s*s12*s16)+2*Abs(p1s*p4s*p5s*s12*s16)+Abs(p2s*
     -   p4s*p5s*s12*s16)+4*Abs(p2s*p4s*p6s*s12*s16)+Abs(p3s*p4s*p6s*
     -   s12*s16)+Abs(p4s**2*p6s*s12*s16)+Abs(p4s**2*s12**2*s16)+Abs(
     -   p4s*p5s*s12**2*s16)+Abs(p2s**2*p4s*s123*s16)+Abs(p2s2*p4s*s1
     -   23*s16)+Abs(p2s**2*p5s*s123*s16)+Abs(p2s2*p5s*s123*s16)+2*Ab
     -   s(p1s*p3s*p5s*s123*s16)+Abs(p2s*p3s*p5s*s123*s16)+2*Abs(p1s*
     -   p4s*p5s*s123*s16)+Abs(p2s*p4s*p5s*s123*s16)+2*Abs(p2s*p3s*p6
     -   s*s123*s16)+2*Abs(p2s*p4s*p6s*s123*s16)+2*Abs(p2s*p4s*s12*s1
     -   23*s16)+2*Abs(p2s*p5s*s12*s123*s16)+2*Abs(p3s*p5s*s12*s123*s
     -   16)+Abs(p4s*p5s*s12*s123*s16)+2*Abs(p2s*p5s*s123**2*s16)+Abs
     -   (p3s*p4s*s12*s16**2)+Abs(p4s**2*s12*s16**2)+2*Abs(p4s*s12**2
     -   *s16**2)+2*Abs(p3s*s12*s123*s16**2)+2*Abs(p4s*s12*s123*s16**
     -   2)+Abs(p1s*p3s*p5s**2*s23)+Abs(p2s**2*p4s*p6s*s23)+Abs(p2s2*
     -   p4s*p6s*s23)+Abs(p2s**2*p5s*p6s*s23)+Abs(p2s2*p5s*p6s*s23)+2
     -   *Abs(p1s*p3s*p5s*p6s*s23)+2*Abs(p2s*p3s*p5s*p6s*s23)+2*Abs(p
     -   1s*p4s*p5s*p6s*s23)+Abs(p2s*p4s*p5s*p6s*s23)+2*Abs(p2s*p3s*p
     -   6s**2*s23)+2*Abs(p2s*p4s*p6s**2*s23)+4*Abs(p1s*p3s*p5s*s12*s
     -   23)+2*Abs(p1s*p4s*p5s*s12*s23)+Abs(p2s*p4s*p5s*s12*s23)+2*Ab
     -   s(p1s*p5s**2*s12*s23)+Abs(p2s*p5s**2*s12*s23)+2*Abs(p2s*p4s*
     -   p6s*s12*s23)+2*Abs(p2s*p5s*p6s*s12*s23)+Abs(p3s*p5s*p6s*s12*
     -   s23)+Abs(p4s*p5s*p6s*s12*s23)+Abs(p4s*p5s*s12**2*s23)+Abs(p5
     -   s**2*s12**2*s23)+2*Abs(p2s**2*p5s*s123*s23)+2*Abs(p2s2*p5s*s
     -   123*s23)+2*Abs(p1s*p5s**2*s123*s23)+Abs(p2s*p5s**2*s123*s23)
     -   +2*Abs(p2s*p5s*p6s*s123*s23)+4*Abs(p2s*p5s*s12*s123*s23)+Abs
     -   (p5s**2*s12*s123*s23)+Abs(p3s*p5s*s12*s16*s23)+2*Abs(p4s*p5s
     -   *s12*s16*s23)+2*Abs(p3s*p6s*s12*s16*s23)+2*Abs(p4s*p6s*s12*s
     -   16*s23)+2*Abs(p4s*s12**2*s16*s23)+2*Abs(p5s*s12**2*s16*s23)+
     -   2*Abs(p5s*s12*s123*s16*s23)+Abs(p1s*p4s**2*s23**2)+2*Abs(p1s
     -   *p4s*p5s*s23**2)+Abs(p1s*p5s**2*s23**2)+2*Abs(p4s*p6s**2*s23
     -   **2)+Abs(p5s**2*s12*s23**2)+2*Abs(p5s*p6s*s12*s23**2)+2*Abs(
     -   p5s*s12**2*s23**2)+2*Abs(p4s*p6s*s123*s23**2)+2*Abs(p5s*p6s*
     -   s123*s23**2)+2*Abs(p5s*s123**2*s23**2)+2*Abs(p4s*p6s*s16*s23
     -   **2)+Abs(p4s*s123*s16*s23**2)+Abs(p5s*s123*s16*s23**2)+Abs(p
     -   4s*p6s*s23**3)+Abs(p5s*p6s*s23**3)+2*Abs(p5s*s123*s23**3)+Ab
     -   s(p1s*p3s**2*p5s*s234)+2*Abs(p1s*p3s**2*p6s*s234)+Abs(p2s**2
     -   *p4s*p6s*s234)+Abs(p2s2*p4s*p6s*s234)+2*Abs(p1s*p3s*p4s*p6s*
     -   s234)+2*Abs(p2s*p3s*p4s*p6s*s234)+Abs(p2s**2*p5s*p6s*s234)+A
     -   bs(p2s2*p5s*p6s*s234)+4*Abs(p1s*p3s*p5s*p6s*s234)+Abs(p2s*p3
     -   s*p5s*p6s*s234)+2*Abs(p2s*p3s*p6s**2*s234)+Abs(p3s**2*p6s**2
     -   *s234)+2*Abs(p2s*p4s*p6s**2*s234)+Abs(p3s*p4s*p6s**2*s234)+2
     -   *Abs(p1s*p3s*p4s*s12*s234)+2*Abs(p1s*p3s*p5s*s12*s234)+2*Abs
     -   (p2s*p4s*p6s*s12*s234)+Abs(p3s*p4s*p6s*s12*s234)+2*Abs(p2s*p
     -   5s*p6s*s12*s234)+Abs(p3s*p5s*p6s*s12*s234)+Abs(p2s**2*p4s*s1
     -   23*s234)+Abs(p2s2*p4s*s123*s234)+Abs(p2s**2*p5s*s123*s234)+A
     -   bs(p2s2*p5s*s123*s234)+2*Abs(p1s*p3s*p5s*s123*s234)+Abs(p2s*
     -   p3s*p5s*s123*s234)+2*Abs(p2s*p3s*p6s*s123*s234)+2*Abs(p2s*p4
     -   s*p6s*s123*s234)+2*Abs(p2s*p5s*p6s*s123*s234)+Abs(p3s*p5s*p6
     -   s*s123*s234)+2*Abs(p2s*p4s*s12*s123*s234)+2*Abs(p2s*p5s*s12*
     -   s123*s234)+2*Abs(p3s*p5s*s12*s123*s234)+2*Abs(p2s*p5s*s123**
     -   2*s234)+Abs(p3s**2*p6s*s16*s234)+Abs(p3s*p4s*p6s*s16*s234)+A
     -   bs(p3s*p4s*s12*s16*s234)+2*Abs(p3s*p5s*s12*s16*s234)+2*Abs(p
     -   3s*p6s*s12*s16*s234)+2*Abs(p4s*p6s*s12*s16*s234)+2*Abs(p4s*s
     -   12**2*s16*s234)+2*Abs(p5s*s12**2*s16*s234)+2*Abs(p3s*p5s*s12
     -   3*s16*s234)+2*Abs(p3s*p6s*s123*s16*s234)+2*Abs(p4s*p6s*s123*
     -   s16*s234)+4*Abs(p3s*s12*s123*s16*s234)+2*Abs(p4s*s12*s123*s1
     -   6*s234)+4*Abs(p5s*s12*s123*s16*s234)+2*Abs(p5s*s123**2*s16*s
     -   234)+Abs(p3s*p5s*p6s*s23*s234)+2*Abs(p3s*p6s**2*s23*s234)+2*
     -   Abs(p4s*p6s**2*s23*s234)+Abs(p3s*p5s*s12*s23*s234)+2*Abs(p3s
     -   *p6s*s12*s23*s234)+4*Abs(p4s*p6s*s12*s23*s234)+2*Abs(p5s*p6s
     -   *s12*s23*s234)+2*Abs(p4s*s12**2*s23*s234)+2*Abs(p5s*s12**2*s
     -   23*s234)+2*Abs(p5s*p6s*s123*s23*s234)+2*Abs(p5s*s12*s123*s23
     -   *s234)+Abs(p4s*p6s*s23**2*s234)+Abs(p5s*p6s*s23**2*s234)+Abs
     -   (p4s*s123*s23**2*s234)+Abs(p5s*s123*s23**2*s234)+Abs(p3s**2*
     -   p6s*s234**2)+2*Abs(p3s*p6s**2*s234**2)+2*Abs(p3s*p6s*s12*s23
     -   4**2)+2*Abs(p3s*p6s*s123*s234**2)+2*Abs(p3s*s12*s123*s234**2
     -   )+Abs(p1s*p4s**2*s23s)+2*Abs(p1s*p4s*p5s*s23s)+Abs(p1s*p5s**
     -   2*s23s)+2*Abs(p4s*p6s**2*s23s)+2*Abs(p4s*p6s*s123*s23s)+2*Ab
     -   s(p5s*p6s*s123*s23s)+2*Abs(p5s*s123**2*s23s)+2*Abs(p4s*p6s*s
     -   16*s23s)+Abs(p4s*s123*s16*s23s)+Abs(p5s*s123*s16*s23s)+Abs(p
     -   4s*p6s*s23*s23s)+Abs(p5s*p6s*s23*s23s)+2*Abs(p5s*s123*s23*s2
     -   3s)+Abs(p4s*p6s*s234*s23s)+Abs(p5s*p6s*s234*s23s)+Abs(p4s*s1
     -   23*s234*s23s)+Abs(p5s*s123*s234*s23s)+2*Abs(p1s**2*p3s*p4s*s
     -   34)+2*Abs(p1s**2*p3s*p5s*s34)+2*Abs(p1s*p2s*p4s*p6s*s34)+2*A
     -   bs(p1s*p3s*p4s*p6s*s34)+2*Abs(p1s*p2s*p5s*p6s*s34)+Abs(p1s*p
     -   3s*p5s*p6s*s34)+Abs(p2s*p4s*p6s**2*s34)+2*Abs(p1s*p2s*p4s*s1
     -   23*s34)+2*Abs(p1s*p2s*p5s*s123*s34)+Abs(p1s*p3s*p5s*s123*s34
     -   )+Abs(p2s*p4s*p6s*s123*s34)+Abs(p2s*p5s*p6s*s123*s34)+Abs(p2
     -   s*p5s*s123**2*s34)+2*Abs(p1s*p3s*p4s*s16*s34)+Abs(p1s*p3s*p5
     -   s*s16*s34)+4*Abs(p1s*p3s*p6s*s16*s34)+4*Abs(p1s*p4s*p6s*s16*
     -   s34)+Abs(p2s*p4s*p6s*s16*s34)+2*Abs(p1s*p4s*s12*s16*s34)+2*A
     -   bs(p1s*p5s*s12*s16*s34)+Abs(p4s*p6s*s12*s16*s34)+2*Abs(p1s*p
     -   3s*s123*s16*s34)+2*Abs(p1s*p4s*s123*s16*s34)+Abs(p2s*p4s*s12
     -   3*s16*s34)+2*Abs(p1s*p5s*s123*s16*s34)+2*Abs(p2s*p5s*s123*s1
     -   6*s34)+2*Abs(p2s*p6s*s123*s16*s34)+Abs(p3s*p6s*s123*s16*s34)
     -   +Abs(p4s*p6s*s123*s16*s34)+2*Abs(p4s*s12*s123*s16*s34)+Abs(p
     -   5s*s12*s123*s16*s34)+2*Abs(p2s*s123**2*s16*s34)+Abs(p5s*s123
     -   **2*s16*s34)+Abs(p4s*s12*s16**2*s34)+Abs(p3s*s123*s16**2*s34
     -   )+Abs(p4s*s123*s16**2*s34)+2*Abs(s12*s123*s16**2*s34)+2*Abs(
     -   s123**2*s16**2*s34)+Abs(p1s*p3s*p5s*s23*s34)+2*Abs(p1s*p3s*p
     -   6s*s23*s34)+2*Abs(p1s*p4s*p6s*s23*s34)+2*Abs(p2s*p4s*p6s*s23
     -   *s34)+2*Abs(p1s*p5s*p6s*s23*s34)+Abs(p2s*p5s*p6s*s23*s34)+2*
     -   Abs(p2s*p6s**2*s23*s34)+Abs(p3s*p6s**2*s23*s34)+Abs(p4s*p6s*
     -   *2*s23*s34)+2*Abs(p1s*p4s*s12*s23*s34)+2*Abs(p1s*p5s*s12*s23
     -   *s34)+Abs(p4s*p6s*s12*s23*s34)+2*Abs(p5s*p6s*s12*s23*s34)+4*
     -   Abs(p1s*p5s*s123*s23*s34)+Abs(p2s*p5s*s123*s23*s34)+2*Abs(p2
     -   s*p6s*s123*s23*s34)+Abs(p5s*p6s*s123*s23*s34)+Abs(p5s*s12*s1
     -   23*s23*s34)+Abs(p3s*p6s*s16*s23*s34)+Abs(p4s*p6s*s16*s23*s34
     -   )+Abs(p4s*s12*s16*s23*s34)+Abs(p5s*s12*s16*s23*s34)+2*Abs(p6
     -   s*s12*s16*s23*s34)+Abs(p5s*s123*s16*s23*s34)+4*Abs(p6s*s123*
     -   s16*s23*s34)+2*Abs(s12*s123*s16*s23*s34)+Abs(p5s*p6s*s23**2*
     -   s34)+2*Abs(p6s**2*s23**2*s34)+Abs(p5s*s12*s23**2*s34)+2*Abs(
     -   p6s*s12*s23**2*s34)+2*Abs(p1s*p3s*p6s*s234*s34)+2*Abs(p2s*p6
     -   s**2*s234*s34)+Abs(p3s*p6s**2*s234*s34)+2*Abs(p1s*p3s*s123*s
     -   234*s34)+4*Abs(p2s*p6s*s123*s234*s34)+Abs(p3s*p6s*s123*s234*
     -   s34)+2*Abs(p2s*s123**2*s234*s34)+Abs(p3s*p6s*s16*s234*s34)+2
     -   *Abs(p6s*s12*s16*s234*s34)+Abs(p3s*s123*s16*s234*s34)+2*Abs(
     -   p6s*s123*s16*s234*s34)+2*Abs(s12*s123*s16*s234*s34)+2*Abs(s1
     -   23**2*s16*s234*s34)+2*Abs(p3s*p6s*s23*s234*s34)+2*Abs(p6s**2
     -   *s23*s234*s34)+2*Abs(p6s*s12*s23*s234*s34)+2*Abs(p6s*s123*s2
     -   3*s234*s34)+2*Abs(s12*s123*s23*s234*s34)+2*Abs(p1s*p6s*s16*s
     -   34**2)+2*Abs(p1s*s123*s16*s34**2)+Abs(p6s*s123*s16*s34**2)+A
     -   bs(s123**2*s16*s34**2)+Abs(s123*s16**2*s34**2)+2*Abs(p1s*p6s
     -   *s23*s34**2)+Abs(p6s**2*s23*s34**2)+2*Abs(p1s*s123*s23*s34**
     -   2)+Abs(p6s*s123*s23*s34**2)+Abs(p6s*s16*s23*s34**2)+Abs(s123
     -   *s16*s23*s34**2)+Abs(p6s*s23**2*s34**2)+2*Abs(p1s**2*p3s*p4s
     -   *s345)+2*Abs(p1s**2*p4s**2*s345)+Abs(p1s*p2s*p4s**2*s345)+2*
     -   Abs(p1s**2*p3s*p5s*s345)+2*Abs(p1s**2*p4s*p5s*s345)+Abs(p1s*
     -   p2s*p4s*p5s*s345)+4*Abs(p1s*p2s*p4s*p6s*s345)+Abs(p1s*p3s*p4
     -   s*p6s*s345)+Abs(p1s*p4s**2*p6s*s345)+Abs(p1s*p4s**2*s12*s345
     -   )+Abs(p1s*p4s*p5s*s12*s345)+2*Abs(p1s*p2s*p4s*s123*s345)+2*A
     -   bs(p1s*p2s*p5s*s123*s345)+Abs(p1s*p3s*p5s*s123*s345)+Abs(p1s
     -   *p4s*p5s*s123*s345)+Abs(p2s*p4s*p6s*s123*s345)+Abs(p2s*p5s*s
     -   123**2*s345)+Abs(p1s*p3s*p4s*s16*s345)+Abs(p1s*p4s**2*s16*s3
     -   45)+4*Abs(p1s*p4s*s12*s16*s345)+2*Abs(p1s*p3s*s123*s16*s345)
     -   +2*Abs(p1s*p4s*s123*s16*s345)+2*Abs(p2s*p4s*s123*s16*s345)+A
     -   bs(p4s*s12*s123*s16*s345)+2*Abs(p2s*s123**2*s16*s345)+Abs(p1
     -   s*p3s*p5s*s23*s345)+Abs(p1s*p4s*p5s*s23*s345)+2*Abs(p1s*p3s*
     -   p6s*s23*s345)+2*Abs(p1s*p4s*p6s*s23*s345)+Abs(p2s*p4s*p6s*s2
     -   3*s345)+2*Abs(p1s*p4s*s12*s23*s345)+2*Abs(p1s*p5s*s12*s23*s3
     -   45)+2*Abs(p4s*p6s*s12*s23*s345)+4*Abs(p1s*p5s*s123*s23*s345)
     -   +Abs(p2s*p5s*s123*s23*s345)+2*Abs(p2s*p6s*s123*s23*s345)+Abs
     -   (p5s*s12*s123*s23*s345)+Abs(p4s*s12*s16*s23*s345)+2*Abs(s12*
     -   s123*s16*s23*s345)+Abs(p5s*s12*s23**2*s345)+2*Abs(p6s*s12*s2
     -   3**2*s345)+Abs(p1s*p3s*p4s*s234*s345)+Abs(p1s*p3s*p5s*s234*s
     -   345)+2*Abs(p1s*p3s*p6s*s234*s345)+2*Abs(p1s*p4s*p6s*s234*s34
     -   5)+Abs(p2s*p4s*p6s*s234*s345)+2*Abs(p1s*p4s*s12*s234*s345)+2
     -   *Abs(p1s*p5s*s12*s234*s345)+2*Abs(p4s*p6s*s12*s234*s345)+2*A
     -   bs(p1s*p3s*s123*s234*s345)+4*Abs(p1s*p4s*s123*s234*s345)+Abs
     -   (p2s*p4s*s123*s234*s345)+2*Abs(p1s*p5s*s123*s234*s345)+2*Abs
     -   (p2s*p5s*s123*s234*s345)+2*Abs(p2s*p6s*s123*s234*s345)+2*Abs
     -   (p3s*p6s*s123*s234*s345)+Abs(p4s*p6s*s123*s234*s345)+Abs(p4s
     -   *s12*s123*s234*s345)+Abs(p5s*s12*s123*s234*s345)+2*Abs(p2s*s
     -   123**2*s234*s345)+Abs(p5s*s123**2*s234*s345)+Abs(p4s*s12*s16
     -   *s234*s345)+Abs(p3s*s123*s16*s234*s345)+Abs(p4s*s123*s16*s23
     -   4*s345)+2*Abs(s12*s123*s16*s234*s345)+2*Abs(s123**2*s16*s234
     -   *s345)+Abs(p3s*p6s*s23*s234*s345)+2*Abs(p4s*p6s*s23*s234*s34
     -   5)+2*Abs(p4s*s12*s23*s234*s345)+Abs(p5s*s12*s23*s234*s345)+4
     -   *Abs(p6s*s12*s23*s234*s345)+Abs(p5s*s123*s23*s234*s345)+2*Ab
     -   s(p6s*s123*s23*s234*s345)+2*Abs(s12*s123*s23*s234*s345)+Abs(
     -   p3s*p6s*s234**2*s345)+2*Abs(p6s*s12*s234**2*s345)+Abs(p3s*s1
     -   23*s234**2*s345)+2*Abs(p6s*s123*s234**2*s345)+2*Abs(s12*s123
     -   *s234**2*s345)+2*Abs(s123**2*s234**2*s345)+2*Abs(p1s**2*p4s*
     -   s34*s345)+2*Abs(p1s**2*p5s*s34*s345)+Abs(p1s*p4s*p6s*s34*s34
     -   5)+Abs(p1s*p4s*s123*s34*s345)+2*Abs(p1s*p5s*s123*s34*s345)+A
     -   bs(p1s*p4s*s16*s34*s345)+2*Abs(p1s*s123*s16*s34*s345)+Abs(s1
     -   23**2*s16*s34*s345)+Abs(p1s*p4s*s23*s34*s345)+2*Abs(p1s*p5s*
     -   s23*s34*s345)+2*Abs(p1s*p6s*s23*s34*s345)+4*Abs(p1s*s123*s23
     -   *s34*s345)+Abs(p6s*s123*s23*s34*s345)+Abs(s123*s16*s23*s34*s
     -   345)+Abs(p6s*s23**2*s34*s345)+2*Abs(p1s*p6s*s234*s34*s345)+2
     -   *Abs(p1s*s123*s234*s34*s345)+Abs(p6s*s123*s234*s34*s345)+Abs
     -   (s123**2*s234*s34*s345)+2*Abs(s123*s16*s234*s34*s345)+Abs(p6
     -   s*s23*s234*s34*s345)+Abs(s123*s23*s234*s34*s345)+2*Abs(p1s**
     -   2*p4s*s345**2)+Abs(p1s*p4s*s123*s345**2)+Abs(p1s*p4s*s23*s34
     -   5**2)+2*Abs(p1s*s123*s23*s345**2)+Abs(p1s*p4s*s234*s345**2)+
     -   2*Abs(p1s*s123*s234*s345**2)+Abs(s123**2*s234*s345**2)+Abs(s
     -   123*s23*s234*s345**2)+Abs(s123*s234**2*s345**2)+2*Abs(p1s*p2
     -   s**2*p4s*s45)+2*Abs(p1s*p2s2*p4s*s45)+2*Abs(p1s*p2s**2*p5s*s
     -   45)+2*Abs(p1s*p2s2*p5s*s45)+2*Abs(p1s**2*p3s*p5s*s45)+Abs(p1
     -   s*p2s*p3s*p5s*s45)+2*Abs(p1s**2*p4s*p5s*s45)+2*Abs(p1s*p2s*p
     -   4s*p5s*s45)+2*Abs(p1s*p2s*p3s*p6s*s45)+2*Abs(p1s*p2s*p4s*p6s
     -   *s45)+Abs(p2s**2*p4s*p6s*s45)+4*Abs(p1s*p2s*p4s*s12*s45)+4*A
     -   bs(p1s*p2s*p5s*s12*s45)+Abs(p1s*p3s*p5s*s12*s45)+2*Abs(p1s*p
     -   4s*p5s*s12*s45)+Abs(p2s*p4s*p6s*s12*s45)+2*Abs(p1s*p2s*p5s*s
     -   123*s45)+Abs(p2s**2*p5s*s123*s45)+2*Abs(p2s2*p6s*s123*s45)+A
     -   bs(p2s*p5s*s12*s123*s45)+2*Abs(p1s*p3s*s12*s16*s45)+2*Abs(p1
     -   s*p4s*s12*s16*s45)+Abs(p2s*p4s*s12*s16*s45)+Abs(p4s*s12**2*s
     -   16*s45)+Abs(p2s**2*s123*s16*s45)+Abs(p2s2*s123*s16*s45)+2*Ab
     -   s(p2s*s12*s123*s16*s45)+Abs(p2s**2*p6s*s23*s45)+Abs(p2s2*p6s
     -   *s23*s45)+2*Abs(p1s*p5s*s12*s23*s45)+Abs(p2s*p5s*s12*s23*s45
     -   )+2*Abs(p2s*p6s*s12*s23*s45)+Abs(p5s*s12**2*s23*s45)+2*Abs(s
     -   12**2*s16*s23*s45)+2*Abs(p1s*p4s*s23**2*s45)+2*Abs(p1s*p5s*s
     -   23**2*s45)+2*Abs(p6s*s123*s23**2*s45)+Abs(s123*s16*s23**2*s4
     -   5)+Abs(p6s*s23**3*s45)+Abs(p1s*p3s*p5s*s234*s45)+Abs(p2s**2*
     -   p6s*s234*s45)+Abs(p2s2*p6s*s234*s45)+2*Abs(p1s*p3s*p6s*s234*
     -   s45)+Abs(p2s*p3s*p6s*s234*s45)+2*Abs(p1s*p4s*p6s*s234*s45)+A
     -   bs(p2s*p4s*p6s*s234*s45)+2*Abs(p1s*p3s*s12*s234*s45)+2*Abs(p
     -   1s*p4s*s12*s234*s45)+Abs(p2s*p4s*s12*s234*s45)+2*Abs(p1s*p5s
     -   *s12*s234*s45)+Abs(p2s*p5s*s12*s234*s45)+2*Abs(p2s*p6s*s12*s
     -   234*s45)+2*Abs(p3s*p6s*s12*s234*s45)+Abs(p4s*p6s*s12*s234*s4
     -   5)+Abs(p4s*s12**2*s234*s45)+Abs(p5s*s12**2*s234*s45)+Abs(p2s
     -   **2*s123*s234*s45)+Abs(p2s2*s123*s234*s45)+2*Abs(p1s*p5s*s12
     -   3*s234*s45)+Abs(p2s*p5s*s123*s234*s45)+4*Abs(p2s*p6s*s123*s2
     -   34*s45)+2*Abs(p2s*s12*s123*s234*s45)+Abs(p5s*s12*s123*s234*s
     -   45)+Abs(p3s*s12*s16*s234*s45)+2*Abs(p4s*s12*s16*s234*s45)+2*
     -   Abs(s12**2*s16*s234*s45)+2*Abs(s12*s123*s16*s234*s45)+2*Abs(
     -   p5s*s12*s23*s234*s45)+2*Abs(p6s*s12*s23*s234*s45)+2*Abs(s12*
     -   *2*s23*s234*s45)+Abs(p6s*s23**2*s234*s45)+Abs(s123*s23**2*s2
     -   34*s45)+Abs(p3s*p6s*s234**2*s45)+Abs(p3s*s12*s234**2*s45)+2*
     -   Abs(p6s*s12*s234**2*s45)+2*Abs(s12**2*s234**2*s45)+2*Abs(p6s
     -   *s123*s234**2*s45)+2*Abs(s12*s123*s234**2*s45)+2*Abs(p1s*p4s
     -   *s23s*s45)+2*Abs(p1s*p5s*s23s*s45)+2*Abs(p6s*s123*s23s*s45)+
     -   Abs(s123*s16*s23s*s45)+Abs(p6s*s23*s23s*s45)+Abs(p6s*s234*s2
     -   3s*s45)+Abs(s123*s234*s23s*s45)+2*Abs(p1s**2*p3s*s34*s45)+2*
     -   Abs(p1s**2*p4s*s34*s45)+Abs(p1s*p2s*p4s*s34*s45)+2*Abs(p1s**
     -   2*p5s*s34*s45)+Abs(p1s*p2s*p5s*s34*s45)+2*Abs(p1s*p2s*p6s*s3
     -   4*s45)+Abs(p1s*p3s*p6s*s34*s45)+Abs(p1s*p4s*p6s*s34*s45)+Abs
     -   (p1s*p4s*s12*s34*s45)+Abs(p1s*p5s*s12*s34*s45)+2*Abs(p1s*p2s
     -   *s123*s34*s45)+Abs(p1s*p5s*s123*s34*s45)+2*Abs(p2s*p6s*s123*
     -   s34*s45)+Abs(p1s*p3s*s16*s34*s45)+Abs(p1s*p4s*s16*s34*s45)+2
     -   *Abs(p1s*s12*s16*s34*s45)+4*Abs(p1s*s123*s16*s34*s45)+Abs(p2
     -   s*s123*s16*s34*s45)+Abs(s12*s123*s16*s34*s45)+Abs(p1s*p5s*s2
     -   3*s34*s45)+4*Abs(p1s*p6s*s23*s34*s45)+Abs(p2s*p6s*s23*s34*s4
     -   5)+2*Abs(p1s*s12*s23*s34*s45)+Abs(p6s*s12*s23*s34*s45)+2*Abs
     -   (s12*s16*s23*s34*s45)+Abs(p1s*p3s*s234*s34*s45)+2*Abs(p1s*p6
     -   s*s234*s34*s45)+2*Abs(p2s*p6s*s234*s34*s45)+4*Abs(p1s*s12*s2
     -   34*s34*s45)+Abs(p6s*s12*s234*s34*s45)+2*Abs(p1s*s123*s234*s3
     -   4*s45)+2*Abs(p2s*s123*s234*s34*s45)+2*Abs(p6s*s123*s234*s34*
     -   s45)+Abs(s12*s123*s234*s34*s45)+Abs(s12*s16*s234*s34*s45)+Ab
     -   s(s123*s16*s234*s34*s45)+Abs(p6s*s23*s234*s34*s45)+Abs(s12*s
     -   23*s234*s34*s45)+2*Abs(p1s**2*s34**2*s45)+Abs(p1s*p6s*s34**2
     -   *s45)+Abs(p1s*s123*s34**2*s45)+Abs(p1s*s16*s34**2*s45)+Abs(p
     -   1s*s23*s34**2*s45)+2*Abs(p1s**2*p3s*s345*s45)+2*Abs(p1s**2*p
     -   4s*s345*s45)+Abs(p1s*p2s*p4s*s345*s45)+Abs(p1s*p4s*s12*s345*
     -   s45)+2*Abs(p1s*p2s*s123*s345*s45)+2*Abs(p1s*s12*s23*s345*s45
     -   )+2*Abs(p1s*p3s*s234*s345*s45)+Abs(p1s*p4s*s234*s345*s45)+2*
     -   Abs(p1s*s12*s234*s345*s45)+2*Abs(p1s*s123*s234*s345*s45)+Abs
     -   (p2s*s123*s234*s345*s45)+2*Abs(s12*s123*s234*s345*s45)+Abs(s
     -   12*s23*s234*s345*s45)+Abs(s12*s234**2*s345*s45)+Abs(s123*s23
     -   4**2*s345*s45)+2*Abs(p1s**2*s34*s345*s45)+Abs(p1s*s123*s34*s
     -   345*s45)+Abs(p1s*s23*s34*s345*s45)+Abs(p1s*s234*s34*s345*s45
     -   )+Abs(p1s*p2s**2*s45**2)+Abs(p1s*p2s2*s45**2)+2*Abs(p1s*p2s*
     -   s12*s45**2)+Abs(p1s*s23**2*s45**2)+2*Abs(p1s*s12*s234*s45**2
     -   )+Abs(p2s*s12*s234*s45**2)+Abs(s12**2*s234*s45**2)+Abs(s12*s
     -   234**2*s45**2)+Abs(p1s*s23s*s45**2)+2*Abs(p1s**2*s34*s45**2)
     -   +Abs(p1s*p2s*s34*s45**2)+Abs(p1s*s12*s34*s45**2)+Abs(p1s*s23
     -   4*s34*s45**2)+2*Abs(p1s*p2s*p3s*p4s*s56)+2*Abs(p1s*p2s*p3s*p
     -   5s*s56)+Abs(p1s*p3s**2*p5s*s56)+2*Abs(p2s2*p4s*p6s*s56)+Abs(
     -   p2s*p3s*p4s*p6s*s56)+2*Abs(p2s2*p5s*p6s*s56)+2*Abs(p2s*p3s*p
     -   5s*p6s*s56)+2*Abs(p2s2*p4s*s123*s56)+2*Abs(p2s2*p5s*s123*s56
     -   )+Abs(p2s*p3s*p5s*s123*s56)+2*Abs(p1s*p3s**2*s16*s56)+Abs(p2
     -   s**2*p4s*s16*s56)+Abs(p2s2*p4s*s16*s56)+2*Abs(p1s*p3s*p4s*s1
     -   6*s56)+Abs(p2s*p3s*p4s*s16*s56)+Abs(p2s**2*p5s*s16*s56)+Abs(
     -   p2s2*p5s*s16*s56)+4*Abs(p1s*p3s*p5s*s16*s56)+Abs(p2s*p3s*p5s
     -   *s16*s56)+2*Abs(p2s*p3s*p6s*s16*s56)+Abs(p3s**2*p6s*s16*s56)
     -   +2*Abs(p2s*p4s*p6s*s16*s56)+Abs(p3s*p4s*p6s*s16*s56)+2*Abs(p
     -   2s*p4s*s12*s16*s56)+2*Abs(p3s*p4s*s12*s16*s56)+2*Abs(p2s*p5s
     -   *s12*s16*s56)+Abs(p3s*p5s*s12*s16*s56)+2*Abs(p2s*p3s*s123*s1
     -   6*s56)+4*Abs(p2s*p4s*s123*s16*s56)+2*Abs(p2s*p5s*s123*s16*s5
     -   6)+Abs(p3s*p5s*s123*s16*s56)+Abs(p3s**2*s16**2*s56)+Abs(p3s*
     -   p4s*s16**2*s56)+2*Abs(p3s*s12*s16**2*s56)+2*Abs(p4s*s12*s16*
     -   *2*s56)+2*Abs(p3s*s123*s16**2*s56)+2*Abs(p4s*s123*s16**2*s56
     -   )+Abs(p2s**2*p4s*s23*s56)+Abs(p2s2*p4s*s23*s56)+Abs(p2s**2*p
     -   5s*s23*s56)+Abs(p2s2*p5s*s23*s56)+2*Abs(p1s*p3s*p5s*s23*s56)
     -   +2*Abs(p2s*p3s*p5s*s23*s56)+4*Abs(p2s*p3s*p6s*s23*s56)+2*Abs
     -   (p2s*p4s*p6s*s23*s56)+4*Abs(p2s*p5s*p6s*s23*s56)+2*Abs(p3s*p
     -   5s*p6s*s23*s56)+2*Abs(p2s*p4s*s12*s23*s56)+2*Abs(p2s*p5s*s12
     -   *s23*s56)+Abs(p3s*p5s*s12*s23*s56)+2*Abs(p2s*p5s*s123*s23*s5
     -   6)+Abs(p3s*p5s*s16*s23*s56)+2*Abs(p3s*p6s*s16*s23*s56)+2*Abs
     -   (p4s*p6s*s16*s23*s56)+2*Abs(p3s*s12*s16*s23*s56)+2*Abs(p4s*s
     -   12*s16*s23*s56)+2*Abs(p5s*s12*s16*s23*s56)+2*Abs(p5s*s123*s1
     -   6*s23*s56)+2*Abs(p4s*p6s*s23**2*s56)+2*Abs(p5s*s12*s23**2*s5
     -   6)+2*Abs(p4s*s123*s23**2*s56)+2*Abs(p5s*s123*s23**2*s56)+Abs
     -   (p4s*s16*s23**2*s56)+Abs(p5s*s16*s23**2*s56)+Abs(p4s*s23**3*
     -   s56)+Abs(p5s*s23**3*s56)+2*Abs(p1s*p3s**2*s234*s56)+2*Abs(p2
     -   s*p3s*p6s*s234*s56)+Abs(p3s**2*p6s*s234*s56)+2*Abs(p2s*p3s*s
     -   123*s234*s56)+Abs(p3s**2*s16*s234*s56)+4*Abs(p3s*p6s*s16*s23
     -   4*s56)+2*Abs(p3s*s12*s16*s234*s56)+2*Abs(p3s*s123*s16*s234*s
     -   56)+2*Abs(p3s*p6s*s23*s234*s56)+2*Abs(p3s*s12*s23*s234*s56)+
     -   2*Abs(p4s*p6s*s23s*s56)+2*Abs(p5s*p6s*s23s*s56)+2*Abs(p4s*s1
     -   23*s23s*s56)+2*Abs(p5s*s123*s23s*s56)+Abs(p4s*s16*s23s*s56)+
     -   Abs(p5s*s16*s23s*s56)+Abs(p4s*s23*s23s*s56)+Abs(p5s*s23*s23s
     -   *s56)+2*Abs(p1s*p3s*s16*s34*s56)+2*Abs(p2s*p6s*s16*s34*s56)+
     -   Abs(p3s*p6s*s16*s34*s56)+2*Abs(p2s*s123*s16*s34*s56)+2*Abs(p
     -   3s*s123*s16*s34*s56)+Abs(p3s*s16**2*s34*s56)+2*Abs(s12*s16**
     -   2*s34*s56)+2*Abs(s123*s16**2*s34*s56)+2*Abs(p1s*p3s*s23*s34*
     -   s56)+2*Abs(p2s*p6s*s23*s34*s56)+Abs(p3s*p6s*s23*s34*s56)+2*A
     -   bs(p2s*s123*s23*s34*s56)+Abs(p3s*s16*s23*s34*s56)+2*Abs(p6s*
     -   s16*s23*s34*s56)+4*Abs(s12*s16*s23*s34*s56)+2*Abs(s123*s16*s
     -   23*s34*s56)+2*Abs(p6s*s23**2*s34*s56)+2*Abs(s12*s23**2*s34*s
     -   56)+2*Abs(p1s*p2s*p4s*s345*s56)+Abs(p1s*p3s*p4s*s345*s56)+2*
     -   Abs(p1s*p2s*p5s*s345*s56)+Abs(p1s*p3s*p5s*s345*s56)+Abs(p2s*
     -   p4s*p6s*s345*s56)+2*Abs(p2s*p4s*s123*s345*s56)+Abs(p2s*p5s*s
     -   123*s345*s56)+2*Abs(p1s*p3s*s16*s345*s56)+2*Abs(p1s*p4s*s16*
     -   s345*s56)+2*Abs(p2s*p4s*s16*s345*s56)+Abs(p4s*s12*s16*s345*s
     -   56)+4*Abs(p2s*s123*s16*s345*s56)+Abs(p3s*s123*s16*s345*s56)+
     -   2*Abs(p4s*s123*s16*s345*s56)+2*Abs(p1s*p3s*s23*s345*s56)+4*A
     -   bs(p1s*p4s*s23*s345*s56)+Abs(p2s*p4s*s23*s345*s56)+2*Abs(p1s
     -   *p5s*s23*s345*s56)+Abs(p2s*p5s*s23*s345*s56)+2*Abs(p2s*p6s*s
     -   23*s345*s56)+Abs(p3s*p6s*s23*s345*s56)+Abs(p4s*p6s*s23*s345*
     -   s56)+Abs(p4s*s12*s23*s345*s56)+2*Abs(p5s*s12*s23*s345*s56)+2
     -   *Abs(p2s*s123*s23*s345*s56)+Abs(p5s*s123*s23*s345*s56)+2*Abs
     -   (p3s*s16*s23*s345*s56)+Abs(p4s*s16*s23*s345*s56)+2*Abs(s12*s
     -   16*s23*s345*s56)+2*Abs(s123*s16*s23*s345*s56)+Abs(p5s*s23**2
     -   *s345*s56)+2*Abs(p6s*s23**2*s345*s56)+2*Abs(s12*s23**2*s345*
     -   s56)+4*Abs(p1s*p3s*s234*s345*s56)+2*Abs(p2s*p6s*s234*s345*s5
     -   6)+Abs(p3s*p6s*s234*s345*s56)+2*Abs(p2s*s123*s234*s345*s56)+
     -   Abs(p3s*s123*s234*s345*s56)+Abs(p3s*s16*s234*s345*s56)+2*Abs
     -   (s12*s16*s234*s345*s56)+2*Abs(s123*s16*s234*s345*s56)+Abs(p3
     -   s*s23*s234*s345*s56)+2*Abs(p6s*s23*s234*s345*s56)+2*Abs(s12*
     -   s23*s234*s345*s56)+4*Abs(s123*s23*s234*s345*s56)+2*Abs(p1s*s
     -   16*s34*s345*s56)+Abs(s123*s16*s34*s345*s56)+2*Abs(p1s*s23*s3
     -   4*s345*s56)+2*Abs(p6s*s23*s34*s345*s56)+Abs(s123*s23*s34*s34
     -   5*s56)+Abs(s16*s23*s34*s345*s56)+Abs(s23**2*s34*s345*s56)+Ab
     -   s(p1s*p4s*s345**2*s56)+2*Abs(p1s*s23*s345**2*s56)+Abs(s123*s
     -   23*s345**2*s56)+Abs(s23**2*s345**2*s56)+2*Abs(p1s*s234*s345*
     -   *2*s56)+Abs(s123*s234*s345**2*s56)+Abs(s23*s234*s345**2*s56)
     -   +2*Abs(p1s*p2s*p3s*s45*s56)+2*Abs(p1s*p2s*p4s*s45*s56)+Abs(p
     -   2s**2*p4s*s45*s56)+2*Abs(p1s*p2s*p5s*s45*s56)+Abs(p2s**2*p5s
     -   *s45*s56)+Abs(p1s*p3s*p5s*s45*s56)+2*Abs(p2s2*p6s*s45*s56)+A
     -   bs(p2s*p3s*p6s*s45*s56)+2*Abs(p2s*p4s*p6s*s45*s56)+Abs(p2s*p
     -   4s*s12*s45*s56)+Abs(p2s*p5s*s12*s45*s56)+2*Abs(p2s2*s123*s45
     -   *s56)+2*Abs(p2s*p5s*s123*s45*s56)+Abs(p2s**2*s16*s45*s56)+Ab
     -   s(p2s2*s16*s45*s56)+2*Abs(p1s*p3s*s16*s45*s56)+2*Abs(p2s*p3s
     -   *s16*s45*s56)+2*Abs(p1s*p4s*s16*s45*s56)+Abs(p2s*p4s*s16*s45
     -   *s56)+2*Abs(p2s*s12*s16*s45*s56)+Abs(p3s*s12*s16*s45*s56)+Ab
     -   s(p4s*s12*s16*s45*s56)+2*Abs(p2s*s123*s16*s45*s56)+Abs(p2s**
     -   2*s23*s45*s56)+Abs(p2s2*s23*s45*s56)+2*Abs(p1s*p5s*s23*s45*s
     -   56)+Abs(p2s*p5s*s23*s45*s56)+2*Abs(p2s*p6s*s23*s45*s56)+2*Ab
     -   s(p2s*s12*s23*s45*s56)+Abs(p5s*s12*s23*s45*s56)+4*Abs(s12*s1
     -   6*s23*s45*s56)+2*Abs(p6s*s23**2*s45*s56)+2*Abs(s123*s23**2*s
     -   45*s56)+Abs(s16*s23**2*s45*s56)+Abs(s23**3*s45*s56)+2*Abs(p2
     -   s**2*s234*s45*s56)+2*Abs(p2s2*s234*s45*s56)+4*Abs(p1s*p3s*s2
     -   34*s45*s56)+Abs(p2s*p3s*s234*s45*s56)+2*Abs(p2s*p6s*s234*s45
     -   *s56)+Abs(p3s*p6s*s234*s45*s56)+4*Abs(p2s*s12*s234*s45*s56)+
     -   Abs(p3s*s12*s234*s45*s56)+2*Abs(p2s*s123*s234*s45*s56)+Abs(p
     -   3s*s16*s234*s45*s56)+2*Abs(s12*s16*s234*s45*s56)+2*Abs(s123*
     -   s16*s234*s45*s56)+2*Abs(p6s*s23*s234*s45*s56)+2*Abs(s12*s23*
     -   s234*s45*s56)+2*Abs(s23**2*s234*s45*s56)+2*Abs(p6s*s23s*s45*
     -   s56)+2*Abs(s123*s23s*s45*s56)+Abs(s16*s23s*s45*s56)+Abs(s23*
     -   s23s*s45*s56)+2*Abs(s234*s23s*s45*s56)+4*Abs(p1s*p2s*s34*s45
     -   *s56)+Abs(p1s*p3s*s34*s45*s56)+Abs(p2s*p6s*s34*s45*s56)+Abs(
     -   p2s*s123*s34*s45*s56)+2*Abs(p1s*s16*s34*s45*s56)+Abs(p2s*s16
     -   *s34*s45*s56)+2*Abs(s12*s16*s34*s45*s56)+Abs(s123*s16*s34*s4
     -   5*s56)+2*Abs(p1s*s23*s34*s45*s56)+Abs(p2s*s23*s34*s45*s56)+A
     -   bs(p6s*s23*s34*s45*s56)+2*Abs(s12*s23*s34*s45*s56)+2*Abs(s16
     -   *s23*s34*s45*s56)+2*Abs(p1s*p2s*s345*s45*s56)+2*Abs(p1s*p3s*
     -   s345*s45*s56)+Abs(p1s*p4s*s345*s45*s56)+Abs(p2s*s123*s345*s4
     -   5*s56)+2*Abs(p1s*s23*s345*s45*s56)+2*Abs(p2s*s23*s345*s45*s5
     -   6)+Abs(s12*s23*s345*s45*s56)+4*Abs(p1s*s234*s345*s45*s56)+Ab
     -   s(p2s*s234*s345*s45*s56)+Abs(s12*s234*s345*s45*s56)+Abs(s123
     -   *s234*s345*s45*s56)+Abs(s23*s234*s345*s45*s56)+Abs(p1s*s34*s
     -   345*s45*s56)+2*Abs(p1s*p2s*s45**2*s56)+Abs(p2s**2*s45**2*s56
     -   )+Abs(p2s*s12*s45**2*s56)+2*Abs(p1s*s234*s45**2*s56)+Abs(p2s
     -   *s234*s45**2*s56)+Abs(s12*s234*s45**2*s56)+Abs(p1s*s34*s45**
     -   2*s56)+2*Abs(p2s*p3s*s16*s56**2)+Abs(p3s**2*s16*s56**2)+2*Ab
     -   s(p3s*s16**2*s56**2)+2*Abs(p2s*p3s*s23*s56**2)+2*Abs(p3s*s16
     -   *s23*s56**2)+2*Abs(p2s*s16*s345*s56**2)+Abs(p3s*s16*s345*s56
     -   **2)+2*Abs(p2s*s23*s345*s56**2)+Abs(p3s*s23*s345*s56**2)+2*A
     -   bs(s16*s23*s345*s56**2)+2*Abs(s23**2*s345*s56**2)+Abs(s23*s3
     -   45**2*s56**2)+2*Abs(p2s2*s45*s56**2)+Abs(p2s*p3s*s45*s56**2)
     -   +2*Abs(p2s*s16*s45*s56**2)+Abs(p3s*s16*s45*s56**2)+2*Abs(p2s
     -   *s23*s45*s56**2)+2*Abs(s16*s23*s45*s56**2)+2*Abs(s23**2*s45*
     -   s56**2)+2*Abs(s23s*s45*s56**2)+Abs(p2s*s345*s45*s56**2)+Abs(
     -   s23*s345*s45*s56**2)+Abs(p2s*s45**2*s56**2)
       cx3=Abs(p1s*p2s**2*p4s**2)+Abs(p1s*p2s2*p4s**2)+2*Abs(p1s**2*p3
     -   s**2*p5s)+2*Abs(p1s*p2s**2*p4s*p5s)+2*Abs(p1s*p2s2*p4s*p5s)+
     -   Abs(p1s**2*p3s*p4s*p5s)+Abs(p1s*p2s*p3s*p4s*p5s)+Abs(p1s*p2s
     -   **2*p5s**2)+Abs(p1s*p2s2*p5s**2)+Abs(p1s**2*p3s*p5s**2)+Abs(
     -   p1s*p2s*p3s*p5s**2)+2*Abs(p1s*p2s*p3s*p4s*p6s)+Abs(p1s*p2s*p
     -   4s**2*p6s)+Abs(p2s**2*p4s**2*p6s)+2*Abs(p1s*p2s*p3s*p5s*p6s)
     -   +Abs(p1s*p3s**2*p5s*p6s)+Abs(p1s*p2s*p4s*p5s*p6s)+Abs(p2s**2
     -   *p4s*p5s*p6s)+2*Abs(p2s2*p4s*p6s**2)+Abs(p2s*p3s*p4s*p6s**2)
     -   +2*Abs(p1s*p2s*p4s**2*s12)+4*Abs(p1s*p2s*p4s*p5s*s12)+Abs(p1
     -   s*p3s*p4s*p5s*s12)+2*Abs(p1s*p2s*p5s**2*s12)+Abs(p1s*p3s*p5s
     -   **2*s12)+Abs(p2s*p4s**2*p6s*s12)+Abs(p2s*p4s*p5s*p6s*s12)+4*
     -   Abs(p1s*p2s*p3s*p5s*s123)+Abs(p1s*p2s*p4s*p5s*s123)+Abs(p2s*
     -   *2*p4s*p5s*s123)+Abs(p1s*p2s*p5s**2*s123)+Abs(p2s**2*p5s**2*
     -   s123)+2*Abs(p2s2*p4s*p6s*s123)+2*Abs(p2s2*p5s*p6s*s123)+Abs(
     -   p2s*p3s*p5s*p6s*s123)+Abs(p2s*p4s*p5s*s12*s123)+Abs(p2s*p5s*
     -   *2*s12*s123)+2*Abs(p2s2*p5s*s123**2)+Abs(p1s*p3s**2*p5s*s16)
     -   +2*Abs(p1s*p3s**2*p6s*s16)+2*Abs(p2s**2*p4s*p6s*s16)+2*Abs(p
     -   2s2*p4s*p6s*s16)+2*Abs(p1s*p3s*p4s*p6s*s16)+Abs(p2s*p3s*p4s*
     -   p6s*s16)+2*Abs(p1s*p3s*p4s*s12*s16)+Abs(p1s*p4s**2*s12*s16)+
     -   Abs(p2s*p4s**2*s12*s16)+2*Abs(p1s*p3s*p5s*s12*s16)+Abs(p1s*p
     -   4s*p5s*s12*s16)+Abs(p2s*p4s*p5s*s12*s16)+4*Abs(p2s*p4s*p6s*s
     -   12*s16)+Abs(p3s*p4s*p6s*s12*s16)+Abs(p4s**2*s12**2*s16)+Abs(
     -   p4s*p5s*s12**2*s16)+Abs(p2s**2*p4s*s123*s16)+Abs(p2s2*p4s*s1
     -   23*s16)+Abs(p2s**2*p5s*s123*s16)+Abs(p2s2*p5s*s123*s16)+Abs(
     -   p1s*p3s*p5s*s123*s16)+Abs(p2s*p3s*p5s*s123*s16)+2*Abs(p2s*p3
     -   s*p6s*s123*s16)+Abs(p2s*p4s*p6s*s123*s16)+2*Abs(p2s*p4s*s12*
     -   s123*s16)+2*Abs(p2s*p5s*s12*s123*s16)+2*Abs(p3s*p5s*s12*s123
     -   *s16)+Abs(p2s*p5s*s123**2*s16)+Abs(p3s*p4s*s12*s16**2)+2*Abs
     -   (p4s*s12**2*s16**2)+2*Abs(p3s*s12*s123*s16**2)+Abs(p4s*s12*s
     -   123*s16**2)+Abs(p2s**2*p4s*p6s*s23)+Abs(p2s2*p4s*p6s*s23)+Ab
     -   s(p2s**2*p5s*p6s*s23)+Abs(p2s2*p5s*p6s*s23)+Abs(p1s*p3s*p5s*
     -   p6s*s23)+2*Abs(p2s*p3s*p5s*p6s*s23)+2*Abs(p2s*p3s*p6s**2*s23
     -   )+Abs(p2s*p4s*p6s**2*s23)+4*Abs(p1s*p3s*p5s*s12*s23)+Abs(p1s
     -   *p4s*p5s*s12*s23)+Abs(p2s*p4s*p5s*s12*s23)+Abs(p1s*p5s**2*s1
     -   2*s23)+Abs(p2s*p5s**2*s12*s23)+2*Abs(p2s*p4s*p6s*s12*s23)+2*
     -   Abs(p2s*p5s*p6s*s12*s23)+Abs(p3s*p5s*p6s*s12*s23)+Abs(p4s*p5
     -   s*s12**2*s23)+Abs(p5s**2*s12**2*s23)+2*Abs(p2s**2*p5s*s123*s
     -   23)+2*Abs(p2s2*p5s*s123*s23)+Abs(p2s*p5s*p6s*s123*s23)+4*Abs
     -   (p2s*p5s*s12*s123*s23)+Abs(p3s*p5s*s12*s16*s23)+2*Abs(p3s*p6
     -   s*s12*s16*s23)+Abs(p4s*p6s*s12*s16*s23)+2*Abs(p4s*s12**2*s16
     -   *s23)+2*Abs(p5s*s12**2*s16*s23)+Abs(p5s*s12*s123*s16*s23)+Ab
     -   s(p5s*p6s*s12*s23**2)+2*Abs(p5s*s12**2*s23**2)+Abs(p1s*p3s**
     -   2*p5s*s234)+2*Abs(p1s*p3s**2*p6s*s234)+Abs(p2s**2*p4s*p6s*s2
     -   34)+Abs(p2s2*p4s*p6s*s234)+Abs(p1s*p3s*p4s*p6s*s234)+2*Abs(p
     -   2s*p3s*p4s*p6s*s234)+Abs(p2s**2*p5s*p6s*s234)+Abs(p2s2*p5s*p
     -   6s*s234)+2*Abs(p1s*p3s*p5s*p6s*s234)+Abs(p2s*p3s*p5s*p6s*s23
     -   4)+2*Abs(p2s*p3s*p6s**2*s234)+Abs(p3s**2*p6s**2*s234)+Abs(p2
     -   s*p4s*p6s**2*s234)+2*Abs(p1s*p3s*p4s*s12*s234)+2*Abs(p1s*p3s
     -   *p5s*s12*s234)+2*Abs(p2s*p4s*p6s*s12*s234)+Abs(p3s*p4s*p6s*s
     -   12*s234)+2*Abs(p2s*p5s*p6s*s12*s234)+Abs(p3s*p5s*p6s*s12*s23
     -   4)+Abs(p2s**2*p4s*s123*s234)+Abs(p2s2*p4s*s123*s234)+Abs(p2s
     -   **2*p5s*s123*s234)+Abs(p2s2*p5s*s123*s234)+Abs(p1s*p3s*p5s*s
     -   123*s234)+Abs(p2s*p3s*p5s*s123*s234)+2*Abs(p2s*p3s*p6s*s123*
     -   s234)+Abs(p2s*p4s*p6s*s123*s234)+Abs(p2s*p5s*p6s*s123*s234)+
     -   2*Abs(p2s*p4s*s12*s123*s234)+2*Abs(p2s*p5s*s12*s123*s234)+2*
     -   Abs(p3s*p5s*s12*s123*s234)+Abs(p2s*p5s*s123**2*s234)+Abs(p3s
     -   **2*p6s*s16*s234)+Abs(p3s*p4s*s12*s16*s234)+2*Abs(p3s*p5s*s1
     -   2*s16*s234)+2*Abs(p3s*p6s*s12*s16*s234)+Abs(p4s*p6s*s12*s16*
     -   s234)+2*Abs(p4s*s12**2*s16*s234)+2*Abs(p5s*s12**2*s16*s234)+
     -   Abs(p3s*p6s*s123*s16*s234)+4*Abs(p3s*s12*s123*s16*s234)+Abs(
     -   p4s*s12*s123*s16*s234)+2*Abs(p5s*s12*s123*s16*s234)+Abs(p3s*
     -   p6s**2*s23*s234)+Abs(p3s*p5s*s12*s23*s234)+2*Abs(p3s*p6s*s12
     -   *s23*s234)+2*Abs(p4s*p6s*s12*s23*s234)+Abs(p5s*p6s*s12*s23*s
     -   234)+2*Abs(p4s*s12**2*s23*s234)+2*Abs(p5s*s12**2*s23*s234)+A
     -   bs(p5s*s12*s123*s23*s234)+Abs(p3s**2*p6s*s234**2)+Abs(p3s*p6
     -   s**2*s234**2)+2*Abs(p3s*p6s*s12*s234**2)+Abs(p3s*p6s*s123*s2
     -   34**2)+2*Abs(p3s*s12*s123*s234**2)+2*Abs(p1s**2*p3s*p4s*s34)
     -   +2*Abs(p1s**2*p3s*p5s*s34)+2*Abs(p1s*p2s*p4s*p6s*s34)+2*Abs(
     -   p1s*p3s*p4s*p6s*s34)+2*Abs(p1s*p2s*p5s*p6s*s34)+Abs(p1s*p3s*
     -   p5s*p6s*s34)+Abs(p2s*p4s*p6s**2*s34)+2*Abs(p1s*p2s*p4s*s123*
     -   s34)+2*Abs(p1s*p2s*p5s*s123*s34)+Abs(p1s*p3s*p5s*s123*s34)+A
     -   bs(p2s*p4s*p6s*s123*s34)+Abs(p2s*p5s*p6s*s123*s34)+Abs(p2s*p
     -   5s*s123**2*s34)+2*Abs(p1s*p3s*p4s*s16*s34)+Abs(p1s*p3s*p5s*s
     -   16*s34)+4*Abs(p1s*p3s*p6s*s16*s34)+2*Abs(p1s*p4s*p6s*s16*s34
     -   )+Abs(p2s*p4s*p6s*s16*s34)+2*Abs(p1s*p4s*s12*s16*s34)+2*Abs(
     -   p1s*p5s*s12*s16*s34)+Abs(p4s*p6s*s12*s16*s34)+2*Abs(p1s*p3s*
     -   s123*s16*s34)+Abs(p1s*p4s*s123*s16*s34)+Abs(p2s*p4s*s123*s16
     -   *s34)+Abs(p1s*p5s*s123*s16*s34)+2*Abs(p2s*p5s*s123*s16*s34)+
     -   2*Abs(p2s*p6s*s123*s16*s34)+Abs(p3s*p6s*s123*s16*s34)+2*Abs(
     -   p4s*s12*s123*s16*s34)+Abs(p5s*s12*s123*s16*s34)+2*Abs(p2s*s1
     -   23**2*s16*s34)+Abs(p4s*s12*s16**2*s34)+Abs(p3s*s123*s16**2*s
     -   34)+2*Abs(s12*s123*s16**2*s34)+Abs(s123**2*s16**2*s34)+Abs(p
     -   1s*p3s*p5s*s23*s34)+2*Abs(p1s*p3s*p6s*s23*s34)+Abs(p1s*p4s*p
     -   6s*s23*s34)+2*Abs(p2s*p4s*p6s*s23*s34)+Abs(p1s*p5s*p6s*s23*s
     -   34)+Abs(p2s*p5s*p6s*s23*s34)+2*Abs(p2s*p6s**2*s23*s34)+Abs(p
     -   3s*p6s**2*s23*s34)+2*Abs(p1s*p4s*s12*s23*s34)+2*Abs(p1s*p5s*
     -   s12*s23*s34)+Abs(p4s*p6s*s12*s23*s34)+2*Abs(p5s*p6s*s12*s23*
     -   s34)+2*Abs(p1s*p5s*s123*s23*s34)+Abs(p2s*p5s*s123*s23*s34)+2
     -   *Abs(p2s*p6s*s123*s23*s34)+Abs(p5s*s12*s123*s23*s34)+Abs(p3s
     -   *p6s*s16*s23*s34)+Abs(p4s*s12*s16*s23*s34)+Abs(p5s*s12*s16*s
     -   23*s34)+2*Abs(p6s*s12*s16*s23*s34)+2*Abs(p6s*s123*s16*s23*s3
     -   4)+2*Abs(s12*s123*s16*s23*s34)+Abs(p6s**2*s23**2*s34)+Abs(p5
     -   s*s12*s23**2*s34)+2*Abs(p6s*s12*s23**2*s34)+2*Abs(p1s*p3s*p6
     -   s*s234*s34)+2*Abs(p2s*p6s**2*s234*s34)+Abs(p3s*p6s**2*s234*s
     -   34)+2*Abs(p1s*p3s*s123*s234*s34)+4*Abs(p2s*p6s*s123*s234*s34
     -   )+Abs(p3s*p6s*s123*s234*s34)+2*Abs(p2s*s123**2*s234*s34)+Abs
     -   (p3s*p6s*s16*s234*s34)+2*Abs(p6s*s12*s16*s234*s34)+Abs(p3s*s
     -   123*s16*s234*s34)+Abs(p6s*s123*s16*s234*s34)+2*Abs(s12*s123*
     -   s16*s234*s34)+Abs(s123**2*s16*s234*s34)+2*Abs(p3s*p6s*s23*s2
     -   34*s34)+Abs(p6s**2*s23*s234*s34)+2*Abs(p6s*s12*s23*s234*s34)
     -   +Abs(p6s*s123*s23*s234*s34)+2*Abs(s12*s123*s23*s234*s34)+2*A
     -   bs(p1s*p6s*s16*s34**2)+2*Abs(p1s*s123*s16*s34**2)+Abs(p6s*s1
     -   23*s16*s34**2)+Abs(s123**2*s16*s34**2)+Abs(s123*s16**2*s34**
     -   2)+2*Abs(p1s*p6s*s23*s34**2)+Abs(p6s**2*s23*s34**2)+2*Abs(p1
     -   s*s123*s23*s34**2)+Abs(p6s*s123*s23*s34**2)+Abs(p6s*s16*s23*
     -   s34**2)+Abs(s123*s16*s23*s34**2)+Abs(p6s*s23**2*s34**2)+2*Ab
     -   s(p1s**2*p3s*p4s*s345)+Abs(p1s**2*p4s**2*s345)+Abs(p1s*p2s*p
     -   4s**2*s345)+2*Abs(p1s**2*p3s*p5s*s345)+Abs(p1s**2*p4s*p5s*s3
     -   45)+Abs(p1s*p2s*p4s*p5s*s345)+4*Abs(p1s*p2s*p4s*p6s*s345)+Ab
     -   s(p1s*p3s*p4s*p6s*s345)+Abs(p1s*p4s**2*s12*s345)+Abs(p1s*p4s
     -   *p5s*s12*s345)+2*Abs(p1s*p2s*p4s*s123*s345)+2*Abs(p1s*p2s*p5
     -   s*s123*s345)+Abs(p1s*p3s*p5s*s123*s345)+Abs(p2s*p4s*p6s*s123
     -   *s345)+Abs(p2s*p5s*s123**2*s345)+Abs(p1s*p3s*p4s*s16*s345)+4
     -   *Abs(p1s*p4s*s12*s16*s345)+2*Abs(p1s*p3s*s123*s16*s345)+Abs(
     -   p1s*p4s*s123*s16*s345)+2*Abs(p2s*p4s*s123*s16*s345)+Abs(p4s*
     -   s12*s123*s16*s345)+2*Abs(p2s*s123**2*s16*s345)+Abs(p1s*p3s*p
     -   5s*s23*s345)+2*Abs(p1s*p3s*p6s*s23*s345)+Abs(p1s*p4s*p6s*s23
     -   *s345)+Abs(p2s*p4s*p6s*s23*s345)+2*Abs(p1s*p4s*s12*s23*s345)
     -   +2*Abs(p1s*p5s*s12*s23*s345)+2*Abs(p4s*p6s*s12*s23*s345)+2*A
     -   bs(p1s*p5s*s123*s23*s345)+Abs(p2s*p5s*s123*s23*s345)+2*Abs(p
     -   2s*p6s*s123*s23*s345)+Abs(p5s*s12*s123*s23*s345)+Abs(p4s*s12
     -   *s16*s23*s345)+2*Abs(s12*s123*s16*s23*s345)+Abs(p5s*s12*s23*
     -   *2*s345)+2*Abs(p6s*s12*s23**2*s345)+Abs(p1s*p3s*p4s*s234*s34
     -   5)+Abs(p1s*p3s*p5s*s234*s345)+2*Abs(p1s*p3s*p6s*s234*s345)+A
     -   bs(p1s*p4s*p6s*s234*s345)+Abs(p2s*p4s*p6s*s234*s345)+2*Abs(p
     -   1s*p4s*s12*s234*s345)+2*Abs(p1s*p5s*s12*s234*s345)+2*Abs(p4s
     -   *p6s*s12*s234*s345)+2*Abs(p1s*p3s*s123*s234*s345)+2*Abs(p1s*
     -   p4s*s123*s234*s345)+Abs(p2s*p4s*s123*s234*s345)+Abs(p1s*p5s*
     -   s123*s234*s345)+2*Abs(p2s*p5s*s123*s234*s345)+2*Abs(p2s*p6s*
     -   s123*s234*s345)+2*Abs(p3s*p6s*s123*s234*s345)+Abs(p4s*s12*s1
     -   23*s234*s345)+Abs(p5s*s12*s123*s234*s345)+2*Abs(p2s*s123**2*
     -   s234*s345)+Abs(p4s*s12*s16*s234*s345)+Abs(p3s*s123*s16*s234*
     -   s345)+2*Abs(s12*s123*s16*s234*s345)+Abs(s123**2*s16*s234*s34
     -   5)+Abs(p3s*p6s*s23*s234*s345)+2*Abs(p4s*s12*s23*s234*s345)+A
     -   bs(p5s*s12*s23*s234*s345)+4*Abs(p6s*s12*s23*s234*s345)+Abs(p
     -   6s*s123*s23*s234*s345)+2*Abs(s12*s123*s23*s234*s345)+Abs(p3s
     -   *p6s*s234**2*s345)+2*Abs(p6s*s12*s234**2*s345)+Abs(p3s*s123*
     -   s234**2*s345)+Abs(p6s*s123*s234**2*s345)+2*Abs(s12*s123*s234
     -   **2*s345)+Abs(s123**2*s234**2*s345)+2*Abs(p1s**2*p4s*s34*s34
     -   5)+2*Abs(p1s**2*p5s*s34*s345)+Abs(p1s*p4s*p6s*s34*s345)+Abs(
     -   p1s*p4s*s123*s34*s345)+2*Abs(p1s*p5s*s123*s34*s345)+Abs(p1s*
     -   p4s*s16*s34*s345)+2*Abs(p1s*s123*s16*s34*s345)+Abs(s123**2*s
     -   16*s34*s345)+Abs(p1s*p4s*s23*s34*s345)+2*Abs(p1s*p5s*s23*s34
     -   *s345)+2*Abs(p1s*p6s*s23*s34*s345)+4*Abs(p1s*s123*s23*s34*s3
     -   45)+Abs(p6s*s123*s23*s34*s345)+Abs(s123*s16*s23*s34*s345)+Ab
     -   s(p6s*s23**2*s34*s345)+2*Abs(p1s*p6s*s234*s34*s345)+2*Abs(p1
     -   s*s123*s234*s34*s345)+Abs(p6s*s123*s234*s34*s345)+Abs(s123**
     -   2*s234*s34*s345)+2*Abs(s123*s16*s234*s34*s345)+Abs(p6s*s23*s
     -   234*s34*s345)+Abs(s123*s23*s234*s34*s345)+2*Abs(p1s**2*p4s*s
     -   345**2)+Abs(p1s*p4s*s123*s345**2)+Abs(p1s*p4s*s23*s345**2)+2
     -   *Abs(p1s*s123*s23*s345**2)+Abs(p1s*p4s*s234*s345**2)+2*Abs(p
     -   1s*s123*s234*s345**2)+Abs(s123**2*s234*s345**2)+Abs(s123*s23
     -   *s234*s345**2)+Abs(s123*s234**2*s345**2)+2*Abs(p1s*p2s**2*p4
     -   s*s45)+2*Abs(p1s*p2s2*p4s*s45)+2*Abs(p1s*p2s**2*p5s*s45)+2*A
     -   bs(p1s*p2s2*p5s*s45)+Abs(p1s**2*p3s*p5s*s45)+Abs(p1s*p2s*p3s
     -   *p5s*s45)+2*Abs(p1s*p2s*p3s*p6s*s45)+Abs(p1s*p2s*p4s*p6s*s45
     -   )+Abs(p2s**2*p4s*p6s*s45)+4*Abs(p1s*p2s*p4s*s12*s45)+4*Abs(p
     -   1s*p2s*p5s*s12*s45)+Abs(p1s*p3s*p5s*s12*s45)+Abs(p2s*p4s*p6s
     -   *s12*s45)+Abs(p1s*p2s*p5s*s123*s45)+Abs(p2s**2*p5s*s123*s45)
     -   +2*Abs(p2s2*p6s*s123*s45)+Abs(p2s*p5s*s12*s123*s45)+2*Abs(p1
     -   s*p3s*s12*s16*s45)+Abs(p1s*p4s*s12*s16*s45)+Abs(p2s*p4s*s12*
     -   s16*s45)+Abs(p4s*s12**2*s16*s45)+Abs(p2s**2*s123*s16*s45)+Ab
     -   s(p2s2*s123*s16*s45)+2*Abs(p2s*s12*s123*s16*s45)+Abs(p2s**2*
     -   p6s*s23*s45)+Abs(p2s2*p6s*s23*s45)+Abs(p1s*p5s*s12*s23*s45)+
     -   Abs(p2s*p5s*s12*s23*s45)+2*Abs(p2s*p6s*s12*s23*s45)+Abs(p5s*
     -   s12**2*s23*s45)+2*Abs(s12**2*s16*s23*s45)+Abs(p2s**2*p6s*s23
     -   4*s45)+Abs(p2s2*p6s*s234*s45)+Abs(p1s*p3s*p6s*s234*s45)+Abs(
     -   p2s*p3s*p6s*s234*s45)+2*Abs(p1s*p3s*s12*s234*s45)+Abs(p1s*p4
     -   s*s12*s234*s45)+Abs(p2s*p4s*s12*s234*s45)+Abs(p1s*p5s*s12*s2
     -   34*s45)+Abs(p2s*p5s*s12*s234*s45)+2*Abs(p2s*p6s*s12*s234*s45
     -   )+2*Abs(p3s*p6s*s12*s234*s45)+Abs(p4s*s12**2*s234*s45)+Abs(p
     -   5s*s12**2*s234*s45)+Abs(p2s**2*s123*s234*s45)+Abs(p2s2*s123*
     -   s234*s45)+2*Abs(p2s*p6s*s123*s234*s45)+2*Abs(p2s*s12*s123*s2
     -   34*s45)+Abs(p3s*s12*s16*s234*s45)+2*Abs(s12**2*s16*s234*s45)
     -   +Abs(s12*s123*s16*s234*s45)+Abs(p6s*s12*s23*s234*s45)+2*Abs(
     -   s12**2*s23*s234*s45)+Abs(p3s*s12*s234**2*s45)+Abs(p6s*s12*s2
     -   34**2*s45)+2*Abs(s12**2*s234**2*s45)+Abs(s12*s123*s234**2*s4
     -   5)+2*Abs(p1s**2*p3s*s34*s45)+Abs(p1s**2*p4s*s34*s45)+Abs(p1s
     -   *p2s*p4s*s34*s45)+Abs(p1s**2*p5s*s34*s45)+Abs(p1s*p2s*p5s*s3
     -   4*s45)+2*Abs(p1s*p2s*p6s*s34*s45)+Abs(p1s*p3s*p6s*s34*s45)+A
     -   bs(p1s*p4s*s12*s34*s45)+Abs(p1s*p5s*s12*s34*s45)+2*Abs(p1s*p
     -   2s*s123*s34*s45)+2*Abs(p2s*p6s*s123*s34*s45)+Abs(p1s*p3s*s16
     -   *s34*s45)+2*Abs(p1s*s12*s16*s34*s45)+2*Abs(p1s*s123*s16*s34*
     -   s45)+Abs(p2s*s123*s16*s34*s45)+Abs(s12*s123*s16*s34*s45)+2*A
     -   bs(p1s*p6s*s23*s34*s45)+Abs(p2s*p6s*s23*s34*s45)+2*Abs(p1s*s
     -   12*s23*s34*s45)+Abs(p6s*s12*s23*s34*s45)+2*Abs(s12*s16*s23*s
     -   34*s45)+Abs(p1s*p3s*s234*s34*s45)+Abs(p1s*p6s*s234*s34*s45)+
     -   2*Abs(p2s*p6s*s234*s34*s45)+4*Abs(p1s*s12*s234*s34*s45)+Abs(
     -   p6s*s12*s234*s34*s45)+Abs(p1s*s123*s234*s34*s45)+2*Abs(p2s*s
     -   123*s234*s34*s45)+Abs(s12*s123*s234*s34*s45)+Abs(s12*s16*s23
     -   4*s34*s45)+Abs(s12*s23*s234*s34*s45)+2*Abs(p1s**2*s34**2*s45
     -   )+Abs(p1s*p6s*s34**2*s45)+Abs(p1s*s123*s34**2*s45)+Abs(p1s*s
     -   16*s34**2*s45)+Abs(p1s*s23*s34**2*s45)+2*Abs(p1s**2*p3s*s345
     -   *s45)+Abs(p1s**2*p4s*s345*s45)+Abs(p1s*p2s*p4s*s345*s45)+Abs
     -   (p1s*p4s*s12*s345*s45)+2*Abs(p1s*p2s*s123*s345*s45)+2*Abs(p1
     -   s*s12*s23*s345*s45)+2*Abs(p1s*p3s*s234*s345*s45)+2*Abs(p1s*s
     -   12*s234*s345*s45)+Abs(p1s*s123*s234*s345*s45)+Abs(p2s*s123*s
     -   234*s345*s45)+2*Abs(s12*s123*s234*s345*s45)+Abs(s12*s23*s234
     -   *s345*s45)+Abs(s12*s234**2*s345*s45)+2*Abs(p1s**2*s34*s345*s
     -   45)+Abs(p1s*s123*s34*s345*s45)+Abs(p1s*s23*s34*s345*s45)+Abs
     -   (p1s*s234*s34*s345*s45)+Abs(p1s*p2s**2*s45**2)+Abs(p1s*p2s2*
     -   s45**2)+2*Abs(p1s*p2s*s12*s45**2)+Abs(p1s*s12*s234*s45**2)+A
     -   bs(p2s*s12*s234*s45**2)+Abs(s12**2*s234*s45**2)+Abs(p1s**2*s
     -   34*s45**2)+Abs(p1s*p2s*s34*s45**2)+Abs(p1s*s12*s34*s45**2)+2
     -   *Abs(p1s*p2s*p3s*p4s*s56)+2*Abs(p1s*p2s*p3s*p5s*s56)+Abs(p1s
     -   *p3s**2*p5s*s56)+2*Abs(p2s2*p4s*p6s*s56)+Abs(p2s*p3s*p4s*p6s
     -   *s56)+2*Abs(p2s2*p5s*p6s*s56)+2*Abs(p2s*p3s*p5s*p6s*s56)+2*A
     -   bs(p2s2*p4s*s123*s56)+2*Abs(p2s2*p5s*s123*s56)+Abs(p2s*p3s*p
     -   5s*s123*s56)+2*Abs(p1s*p3s**2*s16*s56)+Abs(p2s**2*p4s*s16*s5
     -   6)+Abs(p2s2*p4s*s16*s56)+Abs(p1s*p3s*p4s*s16*s56)+Abs(p2s*p3
     -   s*p4s*s16*s56)+Abs(p2s**2*p5s*s16*s56)+Abs(p2s2*p5s*s16*s56)
     -   +2*Abs(p1s*p3s*p5s*s16*s56)+Abs(p2s*p3s*p5s*s16*s56)+2*Abs(p
     -   2s*p3s*p6s*s16*s56)+Abs(p3s**2*p6s*s16*s56)+Abs(p2s*p4s*p6s*
     -   s16*s56)+2*Abs(p2s*p4s*s12*s16*s56)+2*Abs(p3s*p4s*s12*s16*s5
     -   6)+2*Abs(p2s*p5s*s12*s16*s56)+Abs(p3s*p5s*s12*s16*s56)+2*Abs
     -   (p2s*p3s*s123*s16*s56)+2*Abs(p2s*p4s*s123*s16*s56)+Abs(p2s*p
     -   5s*s123*s16*s56)+Abs(p3s**2*s16**2*s56)+2*Abs(p3s*s12*s16**2
     -   *s56)+Abs(p4s*s12*s16**2*s56)+Abs(p3s*s123*s16**2*s56)+Abs(p
     -   2s**2*p4s*s23*s56)+Abs(p2s2*p4s*s23*s56)+Abs(p2s**2*p5s*s23*
     -   s56)+Abs(p2s2*p5s*s23*s56)+Abs(p1s*p3s*p5s*s23*s56)+2*Abs(p2
     -   s*p3s*p5s*s23*s56)+4*Abs(p2s*p3s*p6s*s23*s56)+Abs(p2s*p4s*p6
     -   s*s23*s56)+2*Abs(p2s*p5s*p6s*s23*s56)+2*Abs(p2s*p4s*s12*s23*
     -   s56)+2*Abs(p2s*p5s*s12*s23*s56)+Abs(p3s*p5s*s12*s23*s56)+Abs
     -   (p2s*p5s*s123*s23*s56)+Abs(p3s*p6s*s16*s23*s56)+2*Abs(p3s*s1
     -   2*s16*s23*s56)+Abs(p4s*s12*s16*s23*s56)+Abs(p5s*s12*s16*s23*
     -   s56)+Abs(p5s*s12*s23**2*s56)+2*Abs(p1s*p3s**2*s234*s56)+2*Ab
     -   s(p2s*p3s*p6s*s234*s56)+Abs(p3s**2*p6s*s234*s56)+2*Abs(p2s*p
     -   3s*s123*s234*s56)+Abs(p3s**2*s16*s234*s56)+2*Abs(p3s*p6s*s16
     -   *s234*s56)+2*Abs(p3s*s12*s16*s234*s56)+Abs(p3s*s123*s16*s234
     -   *s56)+Abs(p3s*p6s*s23*s234*s56)+2*Abs(p3s*s12*s23*s234*s56)+
     -   2*Abs(p1s*p3s*s16*s34*s56)+2*Abs(p2s*p6s*s16*s34*s56)+Abs(p3
     -   s*p6s*s16*s34*s56)+2*Abs(p2s*s123*s16*s34*s56)+2*Abs(p3s*s12
     -   3*s16*s34*s56)+Abs(p3s*s16**2*s34*s56)+2*Abs(s12*s16**2*s34*
     -   s56)+Abs(s123*s16**2*s34*s56)+2*Abs(p1s*p3s*s23*s34*s56)+2*A
     -   bs(p2s*p6s*s23*s34*s56)+Abs(p3s*p6s*s23*s34*s56)+2*Abs(p2s*s
     -   123*s23*s34*s56)+Abs(p3s*s16*s23*s34*s56)+Abs(p6s*s16*s23*s3
     -   4*s56)+4*Abs(s12*s16*s23*s34*s56)+Abs(s123*s16*s23*s34*s56)+
     -   Abs(p6s*s23**2*s34*s56)+2*Abs(s12*s23**2*s34*s56)+2*Abs(p1s*
     -   p2s*p4s*s345*s56)+Abs(p1s*p3s*p4s*s345*s56)+2*Abs(p1s*p2s*p5
     -   s*s345*s56)+Abs(p1s*p3s*p5s*s345*s56)+Abs(p2s*p4s*p6s*s345*s
     -   56)+2*Abs(p2s*p4s*s123*s345*s56)+Abs(p2s*p5s*s123*s345*s56)+
     -   2*Abs(p1s*p3s*s16*s345*s56)+Abs(p1s*p4s*s16*s345*s56)+2*Abs(
     -   p2s*p4s*s16*s345*s56)+Abs(p4s*s12*s16*s345*s56)+4*Abs(p2s*s1
     -   23*s16*s345*s56)+Abs(p3s*s123*s16*s345*s56)+2*Abs(p1s*p3s*s2
     -   3*s345*s56)+2*Abs(p1s*p4s*s23*s345*s56)+Abs(p2s*p4s*s23*s345
     -   *s56)+Abs(p1s*p5s*s23*s345*s56)+Abs(p2s*p5s*s23*s345*s56)+2*
     -   Abs(p2s*p6s*s23*s345*s56)+Abs(p3s*p6s*s23*s345*s56)+Abs(p4s*
     -   s12*s23*s345*s56)+2*Abs(p5s*s12*s23*s345*s56)+2*Abs(p2s*s123
     -   *s23*s345*s56)+2*Abs(p3s*s16*s23*s345*s56)+2*Abs(s12*s16*s23
     -   *s345*s56)+Abs(s123*s16*s23*s345*s56)+Abs(p6s*s23**2*s345*s5
     -   6)+2*Abs(s12*s23**2*s345*s56)+4*Abs(p1s*p3s*s234*s345*s56)+2
     -   *Abs(p2s*p6s*s234*s345*s56)+Abs(p3s*p6s*s234*s345*s56)+2*Abs
     -   (p2s*s123*s234*s345*s56)+Abs(p3s*s123*s234*s345*s56)+Abs(p3s
     -   *s16*s234*s345*s56)+2*Abs(s12*s16*s234*s345*s56)+Abs(s123*s1
     -   6*s234*s345*s56)+Abs(p3s*s23*s234*s345*s56)+Abs(p6s*s23*s234
     -   *s345*s56)+2*Abs(s12*s23*s234*s345*s56)+2*Abs(s123*s23*s234*
     -   s345*s56)+2*Abs(p1s*s16*s34*s345*s56)+Abs(s123*s16*s34*s345*
     -   s56)+2*Abs(p1s*s23*s34*s345*s56)+2*Abs(p6s*s23*s34*s345*s56)
     -   +Abs(s123*s23*s34*s345*s56)+Abs(s16*s23*s34*s345*s56)+Abs(s2
     -   3**2*s34*s345*s56)+Abs(p1s*p4s*s345**2*s56)+2*Abs(p1s*s23*s3
     -   45**2*s56)+Abs(s123*s23*s345**2*s56)+Abs(s23**2*s345**2*s56)
     -   +2*Abs(p1s*s234*s345**2*s56)+Abs(s123*s234*s345**2*s56)+Abs(
     -   s23*s234*s345**2*s56)+2*Abs(p1s*p2s*p3s*s45*s56)+Abs(p1s*p2s
     -   *p4s*s45*s56)+Abs(p2s**2*p4s*s45*s56)+Abs(p1s*p2s*p5s*s45*s5
     -   6)+Abs(p2s**2*p5s*s45*s56)+2*Abs(p2s2*p6s*s45*s56)+Abs(p2s*p
     -   3s*p6s*s45*s56)+Abs(p2s*p4s*s12*s45*s56)+Abs(p2s*p5s*s12*s45
     -   *s56)+2*Abs(p2s2*s123*s45*s56)+Abs(p2s**2*s16*s45*s56)+Abs(p
     -   2s2*s16*s45*s56)+Abs(p1s*p3s*s16*s45*s56)+2*Abs(p2s*p3s*s16*
     -   s45*s56)+2*Abs(p2s*s12*s16*s45*s56)+Abs(p3s*s12*s16*s45*s56)
     -   +Abs(p2s*s123*s16*s45*s56)+Abs(p2s**2*s23*s45*s56)+Abs(p2s2*
     -   s23*s45*s56)+Abs(p2s*p6s*s23*s45*s56)+2*Abs(p2s*s12*s23*s45*
     -   s56)+2*Abs(s12*s16*s23*s45*s56)+2*Abs(p2s**2*s234*s45*s56)+2
     -   *Abs(p2s2*s234*s45*s56)+2*Abs(p1s*p3s*s234*s45*s56)+Abs(p2s*
     -   p3s*s234*s45*s56)+Abs(p2s*p6s*s234*s45*s56)+4*Abs(p2s*s12*s2
     -   34*s45*s56)+Abs(p3s*s12*s234*s45*s56)+Abs(p2s*s123*s234*s45*
     -   s56)+Abs(s12*s16*s234*s45*s56)+Abs(s12*s23*s234*s45*s56)+4*A
     -   bs(p1s*p2s*s34*s45*s56)+Abs(p1s*p3s*s34*s45*s56)+Abs(p2s*p6s
     -   *s34*s45*s56)+Abs(p2s*s123*s34*s45*s56)+Abs(p1s*s16*s34*s45*
     -   s56)+Abs(p2s*s16*s34*s45*s56)+2*Abs(s12*s16*s34*s45*s56)+Abs
     -   (p1s*s23*s34*s45*s56)+Abs(p2s*s23*s34*s45*s56)+2*Abs(s12*s23
     -   *s34*s45*s56)+2*Abs(p1s*p2s*s345*s45*s56)+2*Abs(p1s*p3s*s345
     -   *s45*s56)+Abs(p2s*s123*s345*s45*s56)+Abs(p1s*s23*s345*s45*s5
     -   6)+2*Abs(p2s*s23*s345*s45*s56)+Abs(s12*s23*s345*s45*s56)+2*A
     -   bs(p1s*s234*s345*s45*s56)+Abs(p2s*s234*s345*s45*s56)+Abs(s12
     -   *s234*s345*s45*s56)+Abs(p1s*s34*s345*s45*s56)+Abs(p1s*p2s*s4
     -   5**2*s56)+Abs(p2s**2*s45**2*s56)+Abs(p2s*s12*s45**2*s56)+2*A
     -   bs(p2s*p3s*s16*s56**2)+Abs(p3s**2*s16*s56**2)+Abs(p3s*s16**2
     -   *s56**2)+2*Abs(p2s*p3s*s23*s56**2)+Abs(p3s*s16*s23*s56**2)+2
     -   *Abs(p2s*s16*s345*s56**2)+Abs(p3s*s16*s345*s56**2)+2*Abs(p2s
     -   *s23*s345*s56**2)+Abs(p3s*s23*s345*s56**2)+Abs(s16*s23*s345*
     -   s56**2)+Abs(s23**2*s345*s56**2)+Abs(s23*s345**2*s56**2)+2*Ab
     -   s(p2s2*s45*s56**2)+Abs(p2s*p3s*s45*s56**2)+Abs(p2s*s16*s45*s
     -   56**2)+Abs(p2s*s23*s45*s56**2)+Abs(p2s*s345*s45*s56**2)
       cx4=2*Abs(p1s2*p3s2*p5s)+Abs(p1s2*p3s*p4s*p5s)+Abs(p1s*p2s*p3s*
     -   p4s*p5s)+2*Abs(p1s*p2s*p3s*p4s*p6s)+Abs(p1s*p2s*p4s2*p6s)+Ab
     -   s(p2s2*p4s2*p6s)+Abs(p1s*p2s*p3s*p5s*p6s)+Abs(p1s*p3s2*p5s*p
     -   6s)+Abs(p2s2*p4s*p6s2)+Abs(p2s*p3s*p4s*p6s2)+2*Abs(p1s*p2s*p
     -   4s2*s12)+2*Abs(p1s*p2s*p4s*p5s*s12)+Abs(p1s*p3s*p4s*p5s*s12)
     -   +Abs(p2s*p4s2*p6s*s12)+4*Abs(p1s*p2s*p3s*p5s*s123)+Abs(p1s*p
     -   2s*p4s*p5s*s123)+Abs(p2s2*p4s*p5s*s123)+2*Abs(p2s2*p4s*p6s*s
     -   123)+Abs(p2s2*p5s*p6s*s123)+Abs(p2s*p3s*p5s*p6s*s123)+Abs(p2
     -   s*p4s*p5s*s12*s123)+2*Abs(p2s2*p5s*s123**2)+Abs(p1s*p3s2*p5s
     -   *s16)+2*Abs(p1s*p3s2*p6s*s16)+2*Abs(p1s*p3s*p4s*p6s*s16)+Abs
     -   (p2s*p3s*p4s*p6s*s16)+2*Abs(p1s*p3s*p4s*s12*s16)+Abs(p1s*p4s
     -   2*s12*s16)+Abs(p2s*p4s2*s12*s16)+Abs(p1s*p3s*p5s*s12*s16)+2*
     -   Abs(p2s*p4s*p6s*s12*s16)+Abs(p3s*p4s*p6s*s12*s16)+Abs(p1s*p3
     -   s*p5s*s123*s16)+Abs(p2s*p3s*p5s*s123*s16)+2*Abs(p2s*p3s*p6s*
     -   s123*s16)+Abs(p2s*p4s*p6s*s123*s16)+2*Abs(p2s*p4s*s12*s123*s
     -   16)+Abs(p2s*p5s*s12*s123*s16)+2*Abs(p3s*p5s*s12*s123*s16)+Ab
     -   s(p2s*p5s*s123s*s16)+Abs(p4s2*s12s*s16)+Abs(p3s*p4s*s12*s16s
     -   )+2*Abs(p3s*s12*s123*s16s)+Abs(p4s*s12*s123*s16s)+Abs(p4s*s1
     -   2s*s16s)+Abs(p1s*p3s*p5s*p6s*s23)+2*Abs(p2s*p3s*p5s*p6s*s23)
     -   +2*Abs(p2s*p3s*p6s2*s23)+Abs(p2s*p4s*p6s2*s23)+4*Abs(p1s*p3s
     -   *p5s*s12*s23)+Abs(p1s*p4s*p5s*s12*s23)+Abs(p2s*p4s*p5s*s12*s
     -   23)+2*Abs(p2s*p4s*p6s*s12*s23)+Abs(p2s*p5s*p6s*s12*s23)+Abs(
     -   p3s*p5s*p6s*s12*s23)+Abs(p2s*p5s*p6s*s123*s23)+4*Abs(p2s*p5s
     -   *s12*s123*s23)+Abs(p4s*p5s*s12s*s23)+Abs(p3s*p5s*s12*s16*s23
     -   )+2*Abs(p3s*p6s*s12*s16*s23)+Abs(p4s*p6s*s12*s16*s23)+Abs(p5
     -   s*s12*s123*s16*s23)+2*Abs(p4s*s12s*s16*s23)+Abs(p5s*s12s*s16
     -   *s23)+Abs(p1s*p3s2*p5s*s234)+2*Abs(p1s*p3s2*p6s*s234)+Abs(p1
     -   s*p3s*p4s*p6s*s234)+2*Abs(p2s*p3s*p4s*p6s*s234)+Abs(p2s*p3s*
     -   p6s2*s234)+Abs(p3s2*p6s2*s234)+2*Abs(p1s*p3s*p4s*s12*s234)+A
     -   bs(p1s*p3s*p5s*s12*s234)+Abs(p2s*p4s*p6s*s12*s234)+Abs(p3s*p
     -   4s*p6s*s12*s234)+Abs(p1s*p3s*p5s*s123*s234)+Abs(p2s*p3s*p5s*
     -   s123*s234)+2*Abs(p2s*p3s*p6s*s123*s234)+Abs(p2s*p4s*p6s*s123
     -   *s234)+2*Abs(p2s*p4s*s12*s123*s234)+Abs(p2s*p5s*s12*s123*s23
     -   4)+2*Abs(p3s*p5s*s12*s123*s234)+Abs(p2s*p5s*s123s*s234)+Abs(
     -   p3s2*p6s*s16*s234)+Abs(p3s*p4s*s12*s16*s234)+Abs(p3s*p6s*s12
     -   *s16*s234)+Abs(p3s*p6s*s123*s16*s234)+4*Abs(p3s*s12*s123*s16
     -   *s234)+Abs(p4s*s12*s123*s16*s234)+Abs(p4s*s12s*s16*s234)+Abs
     -   (p3s*p6s2*s23*s234)+Abs(p3s*p5s*s12*s23*s234)+2*Abs(p3s*p6s*
     -   s12*s23*s234)+2*Abs(p4s*p6s*s12*s23*s234)+Abs(p5s*s12*s123*s
     -   23*s234)+2*Abs(p4s*s12s*s23*s234)+Abs(p5s*s12s*s23*s234)+Abs
     -   (p3s2*p6s*s234s)+Abs(p3s*p6s*s12*s234s)+Abs(p3s*p6s*s123*s23
     -   4s)+2*Abs(p3s*s12*s123*s234s)+Abs(p5s*p6s*s12*s23s)+2*Abs(p5
     -   s*s12s*s23s)+2*Abs(p1s2*p3s*p4s*s34)+Abs(p1s2*p3s*p5s*s34)+A
     -   bs(p1s*p2s*p4s*p6s*s34)+2*Abs(p1s*p3s*p4s*p6s*s34)+2*Abs(p1s
     -   *p2s*p4s*s123*s34)+Abs(p1s*p2s*p5s*s123*s34)+Abs(p1s*p3s*p5s
     -   *s123*s34)+Abs(p2s*p4s*p6s*s123*s34)+Abs(p2s*p5s*s123s*s34)+
     -   2*Abs(p1s*p3s*p4s*s16*s34)+2*Abs(p1s*p3s*p6s*s16*s34)+Abs(p1
     -   s*p4s*s12*s16*s34)+2*Abs(p1s*p3s*s123*s16*s34)+Abs(p1s*p4s*s
     -   123*s16*s34)+Abs(p2s*p4s*s123*s16*s34)+Abs(p2s*p6s*s123*s16*
     -   s34)+Abs(p3s*p6s*s123*s16*s34)+2*Abs(p4s*s12*s123*s16*s34)+2
     -   *Abs(p2s*s123s*s16*s34)+Abs(p3s*s123*s16s*s34)+Abs(s12*s123*
     -   s16s*s34)+Abs(s123s*s16s*s34)+Abs(p1s*p3s*p5s*s23*s34)+2*Abs
     -   (p1s*p3s*p6s*s23*s34)+Abs(p1s*p4s*p6s*s23*s34)+2*Abs(p2s*p4s
     -   *p6s*s23*s34)+Abs(p2s*p6s2*s23*s34)+Abs(p3s*p6s2*s23*s34)+2*
     -   Abs(p1s*p4s*s12*s23*s34)+Abs(p1s*p5s*s12*s23*s34)+Abs(p4s*p6
     -   s*s12*s23*s34)+2*Abs(p1s*p5s*s123*s23*s34)+Abs(p2s*p5s*s123*
     -   s23*s34)+2*Abs(p2s*p6s*s123*s23*s34)+Abs(p5s*s12*s123*s23*s3
     -   4)+Abs(p3s*p6s*s16*s23*s34)+Abs(p4s*s12*s16*s23*s34)+Abs(p6s
     -   *s12*s16*s23*s34)+2*Abs(p6s*s123*s16*s23*s34)+2*Abs(s12*s123
     -   *s16*s23*s34)+Abs(p1s*p3s*p6s*s234*s34)+2*Abs(p1s*p3s*s123*s
     -   234*s34)+2*Abs(p2s*p6s*s123*s234*s34)+Abs(p3s*p6s*s123*s234*
     -   s34)+2*Abs(p2s*s123s*s234*s34)+Abs(p3s*s123*s16*s234*s34)+Ab
     -   s(s12*s123*s16*s234*s34)+Abs(s123s*s16*s234*s34)+2*Abs(p3s*p
     -   6s*s23*s234*s34)+Abs(p6s*s12*s23*s234*s34)+Abs(p6s*s123*s23*
     -   s234*s34)+2*Abs(s12*s123*s23*s234*s34)+Abs(p6s2*s23s*s34)+Ab
     -   s(p5s*s12*s23s*s34)+2*Abs(p6s*s12*s23s*s34)+2*Abs(p1s2*p3s*p
     -   4s*s345)+Abs(p1s2*p4s2*s345)+Abs(p1s*p2s*p4s2*s345)+Abs(p1s2
     -   *p3s*p5s*s345)+2*Abs(p1s*p2s*p4s*p6s*s345)+Abs(p1s*p3s*p4s*p
     -   6s*s345)+Abs(p1s*p4s2*s12*s345)+2*Abs(p1s*p2s*p4s*s123*s345)
     -   +Abs(p1s*p2s*p5s*s123*s345)+Abs(p1s*p3s*p5s*s123*s345)+Abs(p
     -   2s*p4s*p6s*s123*s345)+Abs(p2s*p5s*s123s*s345)+Abs(p1s*p3s*p4
     -   s*s16*s345)+2*Abs(p1s*p4s*s12*s16*s345)+2*Abs(p1s*p3s*s123*s
     -   16*s345)+Abs(p1s*p4s*s123*s16*s345)+2*Abs(p2s*p4s*s123*s16*s
     -   345)+Abs(p4s*s12*s123*s16*s345)+2*Abs(p2s*s123s*s16*s345)+Ab
     -   s(p1s*p3s*p5s*s23*s345)+2*Abs(p1s*p3s*p6s*s23*s345)+Abs(p1s*
     -   p4s*p6s*s23*s345)+Abs(p2s*p4s*p6s*s23*s345)+2*Abs(p1s*p4s*s1
     -   2*s23*s345)+Abs(p1s*p5s*s12*s23*s345)+2*Abs(p4s*p6s*s12*s23*
     -   s345)+2*Abs(p1s*p5s*s123*s23*s345)+Abs(p2s*p5s*s123*s23*s345
     -   )+2*Abs(p2s*p6s*s123*s23*s345)+Abs(p5s*s12*s123*s23*s345)+Ab
     -   s(p4s*s12*s16*s23*s345)+2*Abs(s12*s123*s16*s23*s345)+Abs(p1s
     -   *p3s*p4s*s234*s345)+Abs(p1s*p3s*p6s*s234*s345)+Abs(p1s*p4s*s
     -   12*s234*s345)+2*Abs(p1s*p3s*s123*s234*s345)+2*Abs(p1s*p4s*s1
     -   23*s234*s345)+Abs(p2s*p4s*s123*s234*s345)+Abs(p2s*p6s*s123*s
     -   234*s345)+2*Abs(p3s*p6s*s123*s234*s345)+Abs(p4s*s12*s123*s23
     -   4*s345)+2*Abs(p2s*s123s*s234*s345)+Abs(p3s*s123*s16*s234*s34
     -   5)+Abs(s12*s123*s16*s234*s345)+Abs(s123s*s16*s234*s345)+Abs(
     -   p3s*p6s*s23*s234*s345)+2*Abs(p4s*s12*s23*s234*s345)+2*Abs(p6
     -   s*s12*s23*s234*s345)+Abs(p6s*s123*s23*s234*s345)+2*Abs(s12*s
     -   123*s23*s234*s345)+Abs(p3s*s123*s234s*s345)+Abs(s12*s123*s23
     -   4s*s345)+Abs(s123s*s234s*s345)+Abs(p5s*s12*s23s*s345)+2*Abs(
     -   p6s*s12*s23s*s345)+Abs(p1s2*p4s*s34*s345)+Abs(p1s*p4s*s123*s
     -   34*s345)+Abs(p1s*s123*s16*s34*s345)+Abs(s123s*s16*s34*s345)+
     -   Abs(p1s*p4s*s23*s34*s345)+Abs(p1s*p6s*s23*s34*s345)+4*Abs(p1
     -   s*s123*s23*s34*s345)+Abs(p6s*s123*s23*s34*s345)+Abs(s123*s16
     -   *s23*s34*s345)+Abs(p1s*s123*s234*s34*s345)+Abs(s123s*s234*s3
     -   4*s345)+Abs(s123*s23*s234*s34*s345)+Abs(p6s*s23s*s34*s345)+A
     -   bs(p1s2*p4s*s345**2)+Abs(p1s*p4s*s123*s345s)+Abs(p1s*p4s*s23
     -   *s345s)+2*Abs(p1s*s123*s23*s345s)+Abs(p1s*s123*s234*s345s)+A
     -   bs(s123s*s234*s345s)+Abs(s123*s23*s234*s345s)+Abs(p1s*s123*s
     -   16*s34s)+Abs(s123s*s16*s34s)+Abs(p1s*p6s*s23*s34s)+2*Abs(p1s
     -   *s123*s23*s34s)+Abs(p6s*s123*s23*s34s)+Abs(s123*s16*s23*s34s
     -   )+Abs(p6s*s23s*s34s)+Abs(p1s2*p3s*p5s*s45)+Abs(p1s*p2s*p3s*p
     -   5s*s45)+2*Abs(p1s*p2s*p3s*p6s*s45)+Abs(p1s*p2s*p4s*p6s*s45)+
     -   Abs(p2s2*p4s*p6s*s45)+4*Abs(p1s*p2s*p4s*s12*s45)+2*Abs(p1s*p
     -   2s*p5s*s12*s45)+Abs(p1s*p3s*p5s*s12*s45)+Abs(p2s*p4s*p6s*s12
     -   *s45)+Abs(p1s*p2s*p5s*s123*s45)+Abs(p2s2*p5s*s123*s45)+2*Abs
     -   (p2s2*p6s*s123*s45)+Abs(p2s*p5s*s12*s123*s45)+2*Abs(p1s*p3s*
     -   s12*s16*s45)+Abs(p1s*p4s*s12*s16*s45)+Abs(p2s*p4s*s12*s16*s4
     -   5)+2*Abs(p2s*s12*s123*s16*s45)+Abs(p4s*s12s*s16*s45)+Abs(p1s
     -   *p5s*s12*s23*s45)+Abs(p2s*p5s*s12*s23*s45)+2*Abs(p2s*p6s*s12
     -   *s23*s45)+Abs(p5s*s12s*s23*s45)+2*Abs(s12s*s16*s23*s45)+Abs(
     -   p1s*p3s*p6s*s234*s45)+Abs(p2s*p3s*p6s*s234*s45)+2*Abs(p1s*p3
     -   s*s12*s234*s45)+Abs(p1s*p4s*s12*s234*s45)+Abs(p2s*p4s*s12*s2
     -   34*s45)+Abs(p2s*p6s*s12*s234*s45)+2*Abs(p3s*p6s*s12*s234*s45
     -   )+2*Abs(p2s*p6s*s123*s234*s45)+2*Abs(p2s*s12*s123*s234*s45)+
     -   Abs(p4s*s12s*s234*s45)+Abs(p3s*s12*s16*s234*s45)+Abs(s12*s12
     -   3*s16*s234*s45)+Abs(s12s*s16*s234*s45)+Abs(p6s*s12*s23*s234*
     -   s45)+2*Abs(s12s*s23*s234*s45)+Abs(p3s*s12*s234s*s45)+Abs(s12
     -   *s123*s234s*s45)+Abs(s12s*s234s*s45)+2*Abs(p1s2*p3s*s34*s45)
     -   +Abs(p1s2*p4s*s34*s45)+Abs(p1s*p2s*p4s*s34*s45)+Abs(p1s*p2s*
     -   p6s*s34*s45)+Abs(p1s*p3s*p6s*s34*s45)+Abs(p1s*p4s*s12*s34*s4
     -   5)+2*Abs(p1s*p2s*s123*s34*s45)+2*Abs(p2s*p6s*s123*s34*s45)+A
     -   bs(p1s*p3s*s16*s34*s45)+Abs(p1s*s12*s16*s34*s45)+2*Abs(p1s*s
     -   123*s16*s34*s45)+Abs(p2s*s123*s16*s34*s45)+Abs(s12*s123*s16*
     -   s34*s45)+2*Abs(p1s*p6s*s23*s34*s45)+Abs(p2s*p6s*s23*s34*s45)
     -   +2*Abs(p1s*s12*s23*s34*s45)+Abs(p6s*s12*s23*s34*s45)+2*Abs(s
     -   12*s16*s23*s34*s45)+Abs(p1s*p3s*s234*s34*s45)+2*Abs(p1s*s12*
     -   s234*s34*s45)+Abs(p1s*s123*s234*s34*s45)+2*Abs(p2s*s123*s234
     -   *s34*s45)+Abs(s12*s123*s234*s34*s45)+Abs(s12*s23*s234*s34*s4
     -   5)+Abs(p1s2*s34**2*s45)+2*Abs(p1s2*p3s*s345*s45)+Abs(p1s2*p4
     -   s*s345*s45)+Abs(p1s*p2s*p4s*s345*s45)+Abs(p1s*p4s*s12*s345*s
     -   45)+2*Abs(p1s*p2s*s123*s345*s45)+2*Abs(p1s*s12*s23*s345*s45)
     -   +2*Abs(p1s*p3s*s234*s345*s45)+Abs(p1s*s12*s234*s345*s45)+Abs
     -   (p1s*s123*s234*s345*s45)+Abs(p2s*s123*s234*s345*s45)+2*Abs(s
     -   12*s123*s234*s345*s45)+Abs(s12*s23*s234*s345*s45)+Abs(p1s2*s
     -   34*s345*s45)+Abs(p1s*s123*s34*s345*s45)+Abs(p1s*s23*s34*s345
     -   *s45)+Abs(p1s*s123*s34s*s45)+Abs(p1s*s23*s34s*s45)+Abs(p1s2*
     -   s34*s45**2)+2*Abs(p1s*p2s*s12*s45s)+Abs(p1s*s12*s234*s45s)+A
     -   bs(p2s*s12*s234*s45s)+Abs(s12s*s234*s45s)+Abs(p1s*p2s*s34*s4
     -   5s)+Abs(p1s*s12*s34*s45s)+2*Abs(p1s*p2s*p3s*p4s*s56)+Abs(p1s
     -   *p2s*p3s*p5s*s56)+Abs(p1s*p3s2*p5s*s56)+Abs(p2s2*p4s*p6s*s56
     -   )+Abs(p2s*p3s*p4s*p6s*s56)+2*Abs(p2s2*p4s*s123*s56)+Abs(p2s2
     -   *p5s*s123*s56)+Abs(p2s*p3s*p5s*s123*s56)+2*Abs(p1s*p3s2*s16*
     -   s56)+Abs(p1s*p3s*p4s*s16*s56)+Abs(p2s*p3s*p4s*s16*s56)+Abs(p
     -   2s*p3s*p6s*s16*s56)+Abs(p3s2*p6s*s16*s56)+Abs(p2s*p4s*s12*s1
     -   6*s56)+2*Abs(p3s*p4s*s12*s16*s56)+2*Abs(p2s*p3s*s123*s16*s56
     -   )+2*Abs(p2s*p4s*s123*s16*s56)+Abs(p3s2*s16s*s56)+Abs(p3s*s12
     -   *s16s*s56)+Abs(p3s*s123*s16s*s56)+Abs(p1s*p3s*p5s*s23*s56)+2
     -   *Abs(p2s*p3s*p5s*s23*s56)+4*Abs(p2s*p3s*p6s*s23*s56)+Abs(p2s
     -   *p4s*p6s*s23*s56)+2*Abs(p2s*p4s*s12*s23*s56)+Abs(p2s*p5s*s12
     -   *s23*s56)+Abs(p3s*p5s*s12*s23*s56)+Abs(p2s*p5s*s123*s23*s56)
     -   +Abs(p3s*p6s*s16*s23*s56)+2*Abs(p3s*s12*s16*s23*s56)+Abs(p4s
     -   *s12*s16*s23*s56)+2*Abs(p1s*p3s2*s234*s56)+Abs(p2s*p3s*p6s*s
     -   234*s56)+Abs(p3s2*p6s*s234*s56)+2*Abs(p2s*p3s*s123*s234*s56)
     -   +Abs(p3s2*s16*s234*s56)+Abs(p3s*s12*s16*s234*s56)+Abs(p3s*s1
     -   23*s16*s234*s56)+Abs(p3s*p6s*s23*s234*s56)+2*Abs(p3s*s12*s23
     -   *s234*s56)+Abs(p5s*s12*s23s*s56)+Abs(p1s*p3s*s16*s34*s56)+Ab
     -   s(p2s*s123*s16*s34*s56)+2*Abs(p3s*s123*s16*s34*s56)+2*Abs(p1
     -   s*p3s*s23*s34*s56)+Abs(p2s*p6s*s23*s34*s56)+Abs(p3s*p6s*s23*
     -   s34*s56)+2*Abs(p2s*s123*s23*s34*s56)+Abs(p3s*s16*s23*s34*s56
     -   )+2*Abs(s12*s16*s23*s34*s56)+Abs(s123*s16*s23*s34*s56)+Abs(p
     -   6s*s23s*s34*s56)+2*Abs(s12*s23s*s34*s56)+Abs(p1s*p2s*p4s*s34
     -   5*s56)+Abs(p1s*p3s*p4s*s345*s56)+2*Abs(p2s*p4s*s123*s345*s56
     -   )+Abs(p1s*p3s*s16*s345*s56)+2*Abs(p2s*s123*s16*s345*s56)+Abs
     -   (p3s*s123*s16*s345*s56)+2*Abs(p1s*p3s*s23*s345*s56)+2*Abs(p1
     -   s*p4s*s23*s345*s56)+Abs(p2s*p4s*s23*s345*s56)+Abs(p2s*p6s*s2
     -   3*s345*s56)+Abs(p3s*p6s*s23*s345*s56)+Abs(p4s*s12*s23*s345*s
     -   56)+2*Abs(p2s*s123*s23*s345*s56)+2*Abs(p3s*s16*s23*s345*s56)
     -   +Abs(s12*s16*s23*s345*s56)+Abs(s123*s16*s23*s345*s56)+2*Abs(
     -   p1s*p3s*s234*s345*s56)+Abs(p2s*s123*s234*s345*s56)+Abs(p3s*s
     -   123*s234*s345*s56)+Abs(p3s*s23*s234*s345*s56)+Abs(s12*s23*s2
     -   34*s345*s56)+2*Abs(s123*s23*s234*s345*s56)+Abs(p6s*s23s*s345
     -   *s56)+2*Abs(s12*s23s*s345*s56)+Abs(p1s*s23*s34*s345*s56)+Abs
     -   (s123*s23*s34*s345*s56)+Abs(s23s*s34*s345*s56)+Abs(p1s*s23*s
     -   345s*s56)+Abs(s123*s23*s345s*s56)+Abs(s23s*s345s*s56)+2*Abs(
     -   p1s*p2s*p3s*s45*s56)+Abs(p1s*p2s*p4s*s45*s56)+Abs(p2s2*p4s*s
     -   45*s56)+Abs(p2s2*p6s*s45*s56)+Abs(p2s*p3s*p6s*s45*s56)+Abs(p
     -   2s*p4s*s12*s45*s56)+2*Abs(p2s2*s123*s45*s56)+Abs(p1s*p3s*s16
     -   *s45*s56)+2*Abs(p2s*p3s*s16*s45*s56)+Abs(p2s*s12*s16*s45*s56
     -   )+Abs(p3s*s12*s16*s45*s56)+Abs(p2s*s123*s16*s45*s56)+Abs(p2s
     -   *p6s*s23*s45*s56)+2*Abs(p2s*s12*s23*s45*s56)+2*Abs(s12*s16*s
     -   23*s45*s56)+2*Abs(p1s*p3s*s234*s45*s56)+Abs(p2s*p3s*s234*s45
     -   *s56)+2*Abs(p2s*s12*s234*s45*s56)+Abs(p3s*s12*s234*s45*s56)+
     -   Abs(p2s*s123*s234*s45*s56)+Abs(s12*s23*s234*s45*s56)+2*Abs(p
     -   1s*p2s*s34*s45*s56)+Abs(p1s*p3s*s34*s45*s56)+Abs(p2s*s123*s3
     -   4*s45*s56)+Abs(p1s*s23*s34*s45*s56)+Abs(p2s*s23*s34*s45*s56)
     -   +2*Abs(s12*s23*s34*s45*s56)+Abs(p1s*p2s*s345*s45*s56)+2*Abs(
     -   p1s*p3s*s345*s45*s56)+Abs(p2s*s123*s345*s45*s56)+Abs(p1s*s23
     -   *s345*s45*s56)+2*Abs(p2s*s23*s345*s45*s56)+Abs(s12*s23*s345*
     -   s45*s56)+Abs(p2s2*s45**2*s56)+Abs(p1s*p2s*s45s*s56)+Abs(p2s*
     -   s12*s45s*s56)+Abs(p2s2*s45*s56**2)+Abs(p2s*p3s*s16*s56s)+Abs
     -   (p3s2*s16*s56s)+2*Abs(p2s*p3s*s23*s56s)+Abs(p3s*s16*s23*s56s
     -   )+Abs(p2s*s23*s345*s56s)+Abs(p3s*s23*s345*s56s)+Abs(s23s*s34
     -   5*s56s)+Abs(p2s*p3s*s45*s56s)+Abs(p2s*s23*s45*s56s)
       cx5=Abs(p1s2*p3s2*p5s)+Abs(p1s2*p3s*p4s*p5s)+Abs(p1s*p2s*p3s*p4
     -   s*p5s)+Abs(p1s*p2s*p3s*p4s*p6s)+Abs(p1s*p2s*p4s2*p6s)+Abs(p2
     -   s2*p4s2*p6s)+2*Abs(p1s*p2s*p4s2*s12)+2*Abs(p1s*p2s*p4s*p5s*s
     -   12)+Abs(p1s*p3s*p4s*p5s*s12)+Abs(p2s*p4s2*p6s*s12)+2*Abs(p1s
     -   *p2s*p3s*p5s*s123)+Abs(p1s*p2s*p4s*p5s*s123)+Abs(p2s2*p4s*p5
     -   s*s123)+Abs(p2s2*p4s*p6s*s123)+Abs(p2s*p4s*p5s*s12*s123)+Abs
     -   (p2s2*p5s*s123**2)+Abs(p1s*p3s*p4s*s12*s16)+Abs(p1s*p4s2*s12
     -   *s16)+Abs(p2s*p4s2*s12*s16)+Abs(p2s*p4s*s12*s123*s16)+Abs(p4
     -   s2*s12s*s16)+2*Abs(p1s*p3s*p5s*s12*s23)+Abs(p1s*p4s*p5s*s12*
     -   s23)+Abs(p2s*p4s*p5s*s12*s23)+Abs(p2s*p4s*p6s*s12*s23)+2*Abs
     -   (p2s*p5s*s12*s123*s23)+Abs(p4s*p5s*s12s*s23)+Abs(p4s*s12s*s1
     -   6*s23)+Abs(p1s*p3s2*p5s*s234)+Abs(p1s*p3s2*p6s*s234)+Abs(p1s
     -   *p3s*p4s*p6s*s234)+2*Abs(p2s*p3s*p4s*p6s*s234)+2*Abs(p1s*p3s
     -   *p4s*s12*s234)+Abs(p1s*p3s*p5s*s12*s234)+Abs(p2s*p4s*p6s*s12
     -   *s234)+Abs(p3s*p4s*p6s*s12*s234)+Abs(p1s*p3s*p5s*s123*s234)+
     -   Abs(p2s*p3s*p5s*s123*s234)+Abs(p2s*p3s*p6s*s123*s234)+Abs(p2
     -   s*p4s*p6s*s123*s234)+2*Abs(p2s*p4s*s12*s123*s234)+Abs(p2s*p5
     -   s*s12*s123*s234)+2*Abs(p3s*p5s*s12*s123*s234)+Abs(p2s*p5s*s1
     -   23**2*s234)+Abs(p3s*p4s*s12*s16*s234)+2*Abs(p3s*s12*s123*s16
     -   *s234)+Abs(p4s*s12*s123*s16*s234)+Abs(p4s*s12s*s16*s234)+Abs
     -   (p3s*p5s*s12*s23*s234)+Abs(p3s*p6s*s12*s23*s234)+2*Abs(p4s*p
     -   6s*s12*s23*s234)+Abs(p5s*s12*s123*s23*s234)+2*Abs(p4s*s12s*s
     -   23*s234)+Abs(p5s*s12s*s23*s234)+Abs(p3s2*p6s*s234s)+Abs(p3s*
     -   p6s*s12*s234s)+Abs(p3s*p6s*s123*s234s)+2*Abs(p3s*s12*s123*s2
     -   34s)+Abs(p5s*s12s*s23s)+2*Abs(p1s2*p3s*p4s*s34)+Abs(p1s2*p3s
     -   *p5s*s34)+Abs(p1s*p2s*p4s*p6s*s34)+2*Abs(p1s*p3s*p4s*p6s*s34
     -   )+2*Abs(p1s*p2s*p4s*s123*s34)+Abs(p1s*p2s*p5s*s123*s34)+Abs(
     -   p1s*p3s*p5s*s123*s34)+Abs(p2s*p4s*p6s*s123*s34)+Abs(p2s*p5s*
     -   s123**2*s34)+2*Abs(p1s*p3s*p4s*s16*s34)+Abs(p1s*p4s*s12*s16*
     -   s34)+Abs(p1s*p3s*s123*s16*s34)+Abs(p1s*p4s*s123*s16*s34)+Abs
     -   (p2s*p4s*s123*s16*s34)+2*Abs(p4s*s12*s123*s16*s34)+Abs(p2s*s
     -   123s*s16*s34)+Abs(p1s*p3s*p5s*s23*s34)+Abs(p1s*p3s*p6s*s23*s
     -   34)+Abs(p1s*p4s*p6s*s23*s34)+2*Abs(p2s*p4s*p6s*s23*s34)+2*Ab
     -   s(p1s*p4s*s12*s23*s34)+Abs(p1s*p5s*s12*s23*s34)+Abs(p4s*p6s*
     -   s12*s23*s34)+2*Abs(p1s*p5s*s123*s23*s34)+Abs(p2s*p5s*s123*s2
     -   3*s34)+Abs(p2s*p6s*s123*s23*s34)+Abs(p5s*s12*s123*s23*s34)+A
     -   bs(p4s*s12*s16*s23*s34)+Abs(s12*s123*s16*s23*s34)+Abs(p1s*p3
     -   s*p6s*s234*s34)+2*Abs(p1s*p3s*s123*s234*s34)+2*Abs(p2s*p6s*s
     -   123*s234*s34)+Abs(p3s*p6s*s123*s234*s34)+2*Abs(p2s*s123s*s23
     -   4*s34)+Abs(p3s*s123*s16*s234*s34)+Abs(s12*s123*s16*s234*s34)
     -   +Abs(s123s*s16*s234*s34)+2*Abs(p3s*p6s*s23*s234*s34)+Abs(p6s
     -   *s12*s23*s234*s34)+Abs(p6s*s123*s23*s234*s34)+2*Abs(s12*s123
     -   *s23*s234*s34)+Abs(p5s*s12*s23s*s34)+Abs(p6s*s12*s23s*s34)+A
     -   bs(p1s2*p3s*p4s*s345)+Abs(p1s2*p4s**2*s345)+Abs(p1s*p2s*p4s2
     -   *s345)+Abs(p1s*p4s2*s12*s345)+Abs(p1s*p2s*p4s*s123*s345)+Abs
     -   (p1s*p4s*s12*s23*s345)+Abs(p1s*p3s*p4s*s234*s345)+Abs(p1s*p4
     -   s*s12*s234*s345)+Abs(p1s*p3s*s123*s234*s345)+2*Abs(p1s*p4s*s
     -   123*s234*s345)+Abs(p2s*p4s*s123*s234*s345)+Abs(p4s*s12*s123*
     -   s234*s345)+Abs(p2s*s123s*s234*s345)+2*Abs(p4s*s12*s23*s234*s
     -   345)+Abs(s12*s123*s23*s234*s345)+Abs(p3s*s123*s234s*s345)+Ab
     -   s(s12*s123*s234s*s345)+Abs(s123s*s234s*s345)+Abs(p1s2*p4s*s3
     -   4*s345)+Abs(p1s*p4s*s123*s34*s345)+Abs(p1s*p4s*s23*s34*s345)
     -   +2*Abs(p1s*s123*s23*s34*s345)+Abs(p1s*s123*s234*s34*s345)+Ab
     -   s(s123s*s234*s34*s345)+Abs(s123*s23*s234*s34*s345)+Abs(p1s*s
     -   123*s16*s34s)+Abs(s123s*s16*s34s)+Abs(p1s*p6s*s23*s34s)+2*Ab
     -   s(p1s*s123*s23*s34s)+Abs(p6s*s123*s23*s34s)+Abs(s123*s16*s23
     -   *s34s)+Abs(p6s*s23s*s34s)+2*Abs(p1s*p2s*p4s*s12*s45)+Abs(p1s
     -   *p3s*s12*s234*s45)+Abs(p1s*p4s*s12*s234*s45)+Abs(p2s*p4s*s12
     -   *s234*s45)+Abs(p2s*s12*s123*s234*s45)+Abs(p4s*s12s*s234*s45)
     -   +Abs(s12s*s23*s234*s45)+Abs(p3s*s12*s234s*s45)+Abs(s12*s123*
     -   s234s*s45)+Abs(s12s*s234s*s45)+Abs(p1s2*p3s*s34*s45)+Abs(p1s
     -   2*p4s*s34*s45)+Abs(p1s*p2s*p4s*s34*s45)+Abs(p1s*p4s*s12*s34*
     -   s45)+Abs(p1s*p2s*s123*s34*s45)+Abs(p1s*s12*s23*s34*s45)+Abs(
     -   p1s*p3s*s234*s34*s45)+2*Abs(p1s*s12*s234*s34*s45)+Abs(p1s*s1
     -   23*s234*s34*s45)+2*Abs(p2s*s123*s234*s34*s45)+Abs(s12*s123*s
     -   234*s34*s45)+Abs(s12*s23*s234*s34*s45)+Abs(p1s2*s34**2*s45)+
     -   Abs(p1s*s123*s34s*s45)+Abs(p1s*s23*s34s*s45)+2*Abs(p1s*p2s*p
     -   3s*p4s*s56)+Abs(p1s*p2s*p3s*p5s*s56)+Abs(p1s*p3s2*p5s*s56)+A
     -   bs(p2s2*p4s*p6s*s56)+Abs(p2s*p3s*p4s*p6s*s56)+2*Abs(p2s2*p4s
     -   *s123*s56)+Abs(p2s2*p5s*s123*s56)+Abs(p2s*p3s*p5s*s123*s56)+
     -   Abs(p1s*p3s2*s16*s56)+Abs(p1s*p3s*p4s*s16*s56)+Abs(p2s*p3s*p
     -   4s*s16*s56)+Abs(p2s*p4s*s12*s16*s56)+2*Abs(p3s*p4s*s12*s16*s
     -   56)+Abs(p2s*p3s*s123*s16*s56)+2*Abs(p2s*p4s*s123*s16*s56)+Ab
     -   s(p1s*p3s*p5s*s23*s56)+2*Abs(p2s*p3s*p5s*s23*s56)+2*Abs(p2s*
     -   p3s*p6s*s23*s56)+Abs(p2s*p4s*p6s*s23*s56)+2*Abs(p2s*p4s*s12*
     -   s23*s56)+Abs(p2s*p5s*s12*s23*s56)+Abs(p3s*p5s*s12*s23*s56)+A
     -   bs(p2s*p5s*s123*s23*s56)+Abs(p3s*s12*s16*s23*s56)+Abs(p4s*s1
     -   2*s16*s23*s56)+2*Abs(p1s*p3s2*s234*s56)+Abs(p2s*p3s*p6s*s234
     -   *s56)+Abs(p3s2*p6s*s234*s56)+2*Abs(p2s*p3s*s123*s234*s56)+Ab
     -   s(p3s2*s16*s234*s56)+Abs(p3s*s12*s16*s234*s56)+Abs(p3s*s123*
     -   s16*s234*s56)+Abs(p3s*p6s*s23*s234*s56)+2*Abs(p3s*s12*s23*s2
     -   34*s56)+Abs(p5s*s12*s23s*s56)+Abs(p1s*p3s*s16*s34*s56)+Abs(p
     -   2s*s123*s16*s34*s56)+2*Abs(p3s*s123*s16*s34*s56)+2*Abs(p1s*p
     -   3s*s23*s34*s56)+Abs(p2s*p6s*s23*s34*s56)+Abs(p3s*p6s*s23*s34
     -   *s56)+2*Abs(p2s*s123*s23*s34*s56)+Abs(p3s*s16*s23*s34*s56)+2
     -   *Abs(s12*s16*s23*s34*s56)+Abs(s123*s16*s23*s34*s56)+Abs(p6s*
     -   s23s*s34*s56)+2*Abs(s12*s23s*s34*s56)+Abs(p1s*p2s*p4s*s345*s
     -   56)+Abs(p1s*p3s*p4s*s345*s56)+2*Abs(p2s*p4s*s123*s345*s56)+A
     -   bs(p1s*p3s*s23*s345*s56)+2*Abs(p1s*p4s*s23*s345*s56)+Abs(p2s
     -   *p4s*s23*s345*s56)+Abs(p4s*s12*s23*s345*s56)+Abs(p2s*s123*s2
     -   3*s345*s56)+2*Abs(p1s*p3s*s234*s345*s56)+Abs(p2s*s123*s234*s
     -   345*s56)+Abs(p3s*s123*s234*s345*s56)+Abs(p3s*s23*s234*s345*s
     -   56)+Abs(s12*s23*s234*s345*s56)+2*Abs(s123*s23*s234*s345*s56)
     -   +Abs(s12*s23s*s345*s56)+Abs(p1s*s23*s34*s345*s56)+Abs(s123*s
     -   23*s34*s345*s56)+Abs(s23s*s34*s345*s56)+Abs(p1s*p2s*p3s*s45*
     -   s56)+Abs(p1s*p2s*p4s*s45*s56)+Abs(p2s2*p4s*s45*s56)+Abs(p2s*
     -   p4s*s12*s45*s56)+Abs(p2s2*s123*s45*s56)+Abs(p2s*s12*s23*s45*
     -   s56)+2*Abs(p1s*p3s*s234*s45*s56)+Abs(p2s*p3s*s234*s45*s56)+2
     -   *Abs(p2s*s12*s234*s45*s56)+Abs(p3s*s12*s234*s45*s56)+Abs(p2s
     -   *s123*s234*s45*s56)+Abs(s12*s23*s234*s45*s56)+2*Abs(p1s*p2s*
     -   s34*s45*s56)+Abs(p1s*p3s*s34*s45*s56)+Abs(p2s*s123*s34*s45*s
     -   56)+Abs(p1s*s23*s34*s45*s56)+Abs(p2s*s23*s34*s45*s56)+2*Abs(
     -   s12*s23*s34*s45*s56)+Abs(p2s2*s45*s56**2)+Abs(p2s*p3s*s16*s5
     -   6s)+Abs(p3s2*s16*s56s)+2*Abs(p2s*p3s*s23*s56s)+Abs(p3s*s16*s
     -   23*s56s)+Abs(p2s*s23*s345*s56s)+Abs(p3s*s23*s345*s56s)+Abs(s
     -   23s*s345*s56s)+Abs(p2s*p3s*s45*s56s)+Abs(p2s*s23*s45*s56s)


        If(abs(x1).ge.abs(x2)) then       
         If(abs(x1).ge.abs(x3)) then        
         If(abs(x1).ge.abs(x4)) then        
         If(abs(x1).ge.abs(x5)) then        
         k=1                               
         X=x1  
         cX=cx1                      
         else                               
         k=5                               
         X=x5  
         cX=cx5                             
        endif                              
         else                               
          If(abs(x4).ge.abs(x5)) then        
          k=4                               
         X=x4    
         cX=cx4                              
        else                               
         k=5                                
        X=x5    
        cX=cx5                              
        endif                              
         endif                              
         else                               
         If(abs(x3).ge.abs(x4)) then        
          If(abs(x3).ge.abs(x5)) then        
         k=3                               
           X=x3   
          cX=cx3                            
        else                               
         k=5                               
         X=x5  
         cX=cx5                         
        endif                              
        else                               
         If(abs(x4).ge.abs(x5)) then        
          k=4                               
         X=x4  
         cX=cx4                        
        else                               
         k=5                               
         X=x5       
          cX=cx5                             
        endif                              
        endif                              
        endif                              
        else                               
         If(abs(x2).ge.abs(x3)) then        
          If(abs(x2).ge.abs(x4)) then        
          If(abs(x2).ge.abs(x5)) then        
         k=2                               
         X=x2     
         cX=cx2                              
        else                               
         k=5                               
         X=x5 
         cX=cx5                              
        endif                              
        else                               
         If(abs(x4).ge.abs(x5)) then        
          k=4                               
         X=x4       
        cX=cx4                              
        else                               
        k=5                                
        X=x5 
        cX=cx5                               
        endif                              
        endif                              
        else                               
         If(abs(x3).ge.abs(x4)) then        
         If(abs(x3).ge.abs(x5)) then        
         k=3                               
         X=x3    
         cX=cx3                              
        else                               
         k=5                               
         X=x5       
         cX=cx5                             
        endif                              
        else                               
         If(abs(x4).ge.abs(x5)) then        
          k=4                               
         X=x4       
         cX=cx4                             
        else                               
         k=5                               
         X=x5   
         cX=cx5                             
        endif                              
        endif                              
        endif                              
        endif                              
       ratio=abs(X)/abs(cX)
       Return
       End


      subroutine Print_Det_up_F(p1,p2,p3,p4,p5,p6)
      Implicit none 
      real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3) 
      real*8 p1sq,p2sq,p3sq,p4sq,p5sq,p6sq
      real*8 s12,s23,s34,s45,s56,s16,s123,s234,s345
      real*8 p1p2,p1p3,p1p4,p1p5,p1p6
      real*8 p2p3,p2p4,p2p5
      real*8 p3p4,p3p5
      real*8 p4p5
      real*8 p5p6
      real*8 ratio(20),tempC(15),lossC(15)
      real*8 ratioD(15),lossD(15)


      real*8 dotrr,Cut
      parameter(Cut=1d-21)
c      Common bad
     
       p1sq = dotrr(p1,p1)
       p1p2 = dotrr(p1,p2)
       p1p3 = dotrr(p1,p3)
       p1p4 = dotrr(p1,p4)
       p1p5 = dotrr(p1,p5)
       p1p6 = dotrr(p1,p6)

       p2sq = dotrr(p2,p2)
       p2p3 = dotrr(p2,p3)
       p2p4 = dotrr(p2,p4)
       p2p5 = dotrr(p2,p5)

       p3sq = dotrr(p3,p3)
       p3p4 = dotrr(p3,p4)
       p3p5 = dotrr(p3,p5)
 
       p4sq = dotrr(p4,p4)
       p4p5 = dotrr(p4,p5)

       p5sq = dotrr(p5,p5)
       p5p6 = dotrr(p5,p6)

       p6sq = dotrr(p6,p6)

       s12 = (p1sq +p2sq+ 2*p1p2) 
       s16 = (p1sq +p6sq+ 2*p1p6) 
       s23 = (p2sq +p3sq+ 2*p2p3) 
       s34 = (p3sq +p4sq+ 2*p3p4) 
       s45 = (p4sq +p5sq+ 2*p4p5) 
       s56 = (p5sq +p6sq+ 2*p5p6) 

      s123=p1sq+p2sq+p3sq+2*(p1p2+p1p3+p2p3) 
      s234=p2sq+p3sq+p4sq+2*(p2p3+p2p4+p3p4) 
      s345=p3sq+p4sq+p5sq+2*(p3p4+p3p5+p4p5)

c       C0123
       call dt3(p1sq,p2sq,s12,ratio(1))
c       C0124
       call dt3(p1sq,s23,s123,ratio(2))
c       C0134
       call dt3(s12,p3sq,s123,ratio(3))
c       C0234
       call dt3(p2sq,p3sq,s23,ratio(4))
c       D01234
       call dt4(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,ratioD(1))

       print*, 'det31',ratio(1)
       print*, 'det32',ratio(2)
       print*, 'det32',ratio(3)
       print*, 'det32',ratio(4)
       print*, 'det41',ratioD(1)

c       C123
        tempC(2)=ratio(1)      
c       C0125
       call dt3(p1sq,s234,s56,ratio(5))
c       C0135
       call dt3(s12,s34,s56,ratio(6))
c       C0235
       call dt3(p2sq,s34,s234,ratio(7))
c       D01235
       call dt4(p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p3+p2p4,ratioD(2))


       print*, 'det31',ratio(1)
       print*, 'det32',ratio(5)
       print*, 'det32',ratio(6)
       print*, 'det32',ratio(7)
       print*, 'det42',ratioD(2)


c       C0124
        tempC(3)=ratio(2)      
c       C0125
c        tempC(3)=ratio(5)      
c       C0145
       call dt3(s123,p4sq,s56,ratio(8))
c       C0245
       call dt3(s23,p4sq,s234,ratio(9))
c       D01245
       call dt4(p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p4+p3p4,ratioD(3))


       print*, 'det31',ratio(2)
       print*, 'det32',ratio(5)
       print*, 'det32',ratio(8)
       print*, 'det32',ratio(9)
       print*, 'det43',ratioD(3)


c       C0134
        tempC(4)=ratio(3)
c       C0135
c        tempC(4)=ratio(6)
c       C0145
c        tempC(4)=ratio(8)
c       C0345
        call dt3(p3sq,p4sq,s34,ratio(10))
c       D01345
       call dt4(s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4,p3p4,ratioD(4))
        

       print*, 'det31',ratio(3)
       print*, 'det32',ratio(6)
       print*, 'det32',ratio(8)
       print*, 'det32',ratio(10)
       print*, 'det44',ratioD(4)

c       C0234
        tempC(5)=ratio(4)
c       C0235
c        tempC(5)=ratio(7)
c       C0245
c        tempC(5)=ratio(9)
c       C0345
c        tempC(5)=ratio(10)
c       D02345
       call dt4(p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,ratioD(5))

       print*, 'det31',ratio(4)
       print*, 'det32',ratio(7)
       print*, 'det32',ratio(9)
       print*, 'det32',ratio(10)
       print*, 'det45',ratioD(5)

c       C123
        tempC(6)=ratio(1)    
c       C0126
       call dt3(p1sq,s16,p6sq,ratio(11))
c       C0136
       call dt3(s12,s345,p6sq,ratio(12))
c       C0236
       call dt3(p2sq,s345,s16,ratio(13))

c       D01236
       call dt4(p1sq,p2sq,s345,p1p2,p1p3+p1p4+p1p5,p2p3+p2p4+p2p5,ratioD(6))


       print*, 'det31',ratio(1)
       print*, 'det32',ratio(11)
       print*, 'det32',ratio(12)
       print*, 'det32',ratio(13)
       print*, 'det46',ratioD(6)


c       C0124
        tempC(7)=ratio(2)  
c       C0126
c       call dt3(p1sq,s16,p6sq,ratio(11))
c       C0146
        call dt3(s123,s45,p6sq,ratio(14))
c       C0246
        call dt3(s23,s45,s16,ratio(15))
c       D01246
       call dt4(p1sq,s23,s45,p1p2+p1p3,p1p4+p1p5,p2p4+p2p5+p3p4+p3p5,ratioD(7))
       print*, 'det31',ratio(2)
       print*, 'det32',ratio(11)
       print*, 'det32',ratio(14)
       print*, 'det32',ratio(15)
       print*, 'det47',ratioD(7)

c       C0134
        tempC(8)=ratio(3)
c       C0136
c       call dt3(s12,s345,p6sq,ratio(12))
c       C0146
c        call dt3(s123,s45,p6sq,ratio(14))
c       C0346
       call dt3(p3sq,s45,s345,ratio(16))

c       D01346
       call dt4(s12,p3sq,s45,p1p3+p2p3,p1p4+p1p5+p2p4+p2p5,p3p4+p3p5,ratioD(8))

       print*, 'det31',ratio(3)
       print*, 'det32',ratio(12)
       print*, 'det32',ratio(14)
       print*, 'det32',ratio(16)
       print*, 'det48',ratioD(8)


c       C0234
        tempC(9)=ratio(4)
c       C0236
c       call dt3(p2sq,s345,s16,ratio(13))
c       C0246
c        call dt3(s23,s45,s16,ratio(15))
c       C0346
c       call dt3(p3sq,s45,s345,ratio(16))


c       D02346
       call dt4(p2sq,p3sq,s45,p2p3,p2p4+p2p5,p3p4+p3p5,ratioD(9))

       print*, 'det31',ratio(4)
       print*, 'det32',ratio(13)
       print*, 'det32',ratio(15)
       print*, 'det32',ratio(16)
       print*, 'det49',ratioD(9)


c       C0125
        tempC(10)=ratio(5)
c       C0126
c       call dt3(p1sq,s16,p6sq,ratio(11))
c       C0156
       call dt3(s56,p5sq,p6sq,ratio(17))
c       C0256
       call dt3(s234,p5sq,s16,ratio(18))
c       D01256
       call dt4(p1sq,s234,p5sq,p1p2+p1p3+p1p4,p1p5,p2p5+p3p5+p4p5,ratioD(10))

       print*, 'det31',ratio(5)
       print*, 'det32',ratio(11)
       print*, 'det32',ratio(17)
       print*, 'det32',ratio(18)
       print*, 'det410',ratioD(10)


c       C0135
        tempC(11)=ratio(6)
c       C0136
c       call dt3(s12,s345,p6sq,ratio(12))
c       C0156
c       call dt3(s56,p5sq,p6sq,ratio(17))
c       C0356
       call dt3(s34,p5sq,s345,ratio(19))
c       D01356
       call dt4(s12,s34,p5sq,p1p3+p1p4+p2p3+p2p4,p1p5+p2p5,p3p5+p4p5,ratioD(11))

       print*, 'det31',ratio(6)
       print*, 'det32',ratio(12)
       print*, 'det32',ratio(17)
       print*, 'det32',ratio(19)
       print*, 'det411',ratioD(11)


c       C0235
        tempC(12)=ratio(7)
c       C0236
c       call dt3(p2sq,s345,s16,ratio(13))
c       C0256
       call dt3(s234,p5sq,s16,ratio(18))
c       C0356
       call dt3(s34,p5sq,s345,ratio(19))
c       D02356
       call dt4(p2sq,s34,p5sq,p2p3+p2p4,p2p5,p3p5+p4p5,ratioD(12))

       print*, 'det31',ratio(7)
       print*, 'det32',ratio(13)
       print*, 'det32',ratio(18)
       print*, 'det32',ratio(19)
       print*, 'det412',ratioD(12)


c       C0145
      tempC(13)=ratio(8)
c       C0146
c        call dt3(s123,s45,p6sq,ratio(14))
c       C0156
c       call dt3(s56,p5sq,p6sq,ratio(17))
c       C0456
       call dt3(p4sq,p5sq,s45,ratio(20))
c       D01456
       call dt4(s123,p4sq,p5sq,p1p4+p2p4+p3p4,p1p5+p2p5+p3p5,p4p5,ratioD(13))


       print*, 'det31',ratio(8)
       print*, 'det32',ratio(14)
       print*, 'det32',ratio(17)
       print*, 'det32',ratio(20)
       print*, 'det413',ratioD(13)

c      C0245
       tempC(14)=9
c      C0246
c       tempC(14)=15
c      C0256
c       tempC(14)=18
c      C0456
c       tempC(14)=20
c       D02456
       call dt4(s23,p4sq,p5sq,p2p4+p3p4,p2p5+p3p5,p4p5,ratioD(14))

       print*, 'det31',ratio(9)
       print*, 'det32',ratio(15)
       print*, 'det32',ratio(18)
       print*, 'det32',ratio(20)
       print*, 'det414',ratioD(14)

c      C0345
       tempC(15)=ratio(10)   
c      C0346
c       tempC(15)=ratio(16)
c      C0356
c       tempC(15)=ratio(19)
c      C0456
c     tempC(15)=ratio(20)
c       D03456
       call dt4(p3sq,p4sq,p5sq,p3p4,p3p5,p4p5,ratioD(15))
       lossD(15)=(ratioD(15)*ratioD(15)*ratioD(15))*lossC(15)

       print*, 'det31',ratio(10)
       print*, 'det32',ratio(16)
       print*, 'det32',ratio(19)
       print*, 'det32',ratio(20)
       print*, 'det415',ratioD(15)


      end


      subroutine Calc_Det_up_F(p1,p2,p3,p4,p5,p6,k,countbad,seedbad)
      Implicit none 
      real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3) 
      real*8 p1sq,p2sq,p3sq,p4sq,p5sq,p6sq
      real*8 s12,s23,s34,s45,s56,s16,s123,s234,s345
      real*8 p1p2,p1p3,p1p4,p1p5,p1p6
      real*8 p2p3,p2p4,p2p5
      real*8 p3p4,p3p5
      real*8 p4p5
      real*8 p5p6
      real*8 ratio(20),tempC(15),lossC(15)
      real*8 ratioD(15),lossD(15),tempD(6)
      real*8 lossE(6),tempE(1)

      logical bad
      Integer seedbad(1000000),countbad,k
      real*8 dotrr,Cut
      parameter(Cut=1d-21)
      Common bad
     
       p1sq = dotrr(p1,p1)
       p1p2 = dotrr(p1,p2)
       p1p3 = dotrr(p1,p3)
       p1p4 = dotrr(p1,p4)
       p1p5 = dotrr(p1,p5)
       p1p6 = dotrr(p1,p6)

       p2sq = dotrr(p2,p2)
       p2p3 = dotrr(p2,p3)
       p2p4 = dotrr(p2,p4)
       p2p5 = dotrr(p2,p5)

       p3sq = dotrr(p3,p3)
       p3p4 = dotrr(p3,p4)
       p3p5 = dotrr(p3,p5)
 
       p4sq = dotrr(p4,p4)
       p4p5 = dotrr(p4,p5)

       p5sq = dotrr(p5,p5)
       p5p6 = dotrr(p5,p6)

       p6sq = dotrr(p6,p6)

       s12 = (p1sq +p2sq+ 2*p1p2) 
       s16 = (p1sq +p6sq+ 2*p1p6) 
       s23 = (p2sq +p3sq+ 2*p2p3) 
       s34 = (p3sq +p4sq+ 2*p3p4) 
       s45 = (p4sq +p5sq+ 2*p4p5) 
       s56 = (p5sq +p6sq+ 2*p5p6) 

      s123=p1sq+p2sq+p3sq+2*(p1p2+p1p3+p2p3) 
      s234=p2sq+p3sq+p4sq+2*(p2p3+p2p4+p3p4) 
      s345=p3sq+p4sq+p5sq+2*(p3p4+p3p5+p4p5)

      bad=.False.

c       C0123
       call dt3(p1sq,p2sq,s12,ratio(1))
c       C0124
       call dt3(p1sq,s23,s123,ratio(2))
c       C0134
       call dt3(s12,p3sq,s123,ratio(3))
c       C0234
       call dt3(p2sq,p3sq,s23,ratio(4))

       tempC(1)=ratio(1)
       If(ratio(2).le.tempC(1)) tempC(1)=ratio(2)
       If(ratio(3).le.tempC(1)) tempC(1)=ratio(3)
       If(ratio(4).le.tempC(1)) tempC(1)=ratio(4)
       lossC(1)=tempC(1)*tempC(1)
c       D01234
       call dt4(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,ratioD(1))
       lossD(1)=(ratioD(1)*ratioD(1)*ratioD(1))*lossC(1)

       If(lossD(1).le.Cut)  bad=.True.
       
       If (bad) Then
c       print*, "k", k
c       print*, "lossD(1)", lossD(1)
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else

c       C123
        tempC(2)=ratio(1)      
c       C0125
       call dt3(p1sq,s234,s56,ratio(5))
c       C0135
       call dt3(s12,s34,s56,ratio(6))
c       C0235
       call dt3(p2sq,s34,s234,ratio(7))

       If(ratio(5).le.tempC(2)) tempC(2)=ratio(5)
       If(ratio(6).le.tempC(2)) tempC(2)=ratio(6)
       If(ratio(7).le.tempC(2)) tempC(2)=ratio(7)
       lossC(2)=tempC(2)*tempC(2)
c       D01235
       call dt4(p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p3+p2p4,ratioD(2))
       lossD(2)=(ratioD(2)*ratioD(2)*ratioD(2))*lossC(2)
       If(lossD(2).le.Cut)  bad=.True.

       If (bad) Then
c       print*, "k", k
c       print*, "lossD(2)", lossD(2)
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else

c       C0124
        tempC(3)=ratio(2)      
c       C0125
c        tempC(3)=ratio(5)      
c       C0145
       call dt3(s123,p4sq,s56,ratio(8))
c       C0245
       call dt3(s23,p4sq,s234,ratio(9))
       If(ratio(5).le.tempC(3)) tempC(3)=ratio(5)
       If(ratio(8).le.tempC(3)) tempC(3)=ratio(8)
       If(ratio(9).le.tempC(3)) tempC(3)=ratio(9)
       lossC(3)=tempC(3)*tempC(3)
c       D01245
       call dt4(p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p4+p3p4,ratioD(3))
       lossD(3)=(ratioD(3)*ratioD(3)*ratioD(3))*lossC(3)
       If(lossD(3).le.Cut)  bad=.True.
       If (bad) Then
c       print*, "k", k
c       print*, "lossD(3)", lossD(3)
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else

c       C0134
        tempC(4)=ratio(3)
c       C0135
c        tempC(4)=ratio(6)
c       C0145
c        tempC(4)=ratio(8)
c       C0345
        call dt3(p3sq,p4sq,s34,ratio(10))
       If(ratio(6).le.tempC(4)) tempC(4)=ratio(6)
       If(ratio(8).le.tempC(4)) tempC(4)=ratio(8)
       If(ratio(10).le.tempC(4)) tempC(4)=ratio(10)
       lossC(4)=tempC(4)*tempC(4)
c       D01345
       call dt4(s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4,p3p4,ratioD(4))
       lossD(4)=(ratioD(4)*ratioD(4)*ratioD(4))*lossC(4)
       If(lossD(4).le.Cut)  bad=.True.
       If (bad) Then
c       print*, "k", k
c       print*, "lossD(4)", lossD(4)
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else

c       C0234
        tempC(5)=ratio(4)
c       C0235
c        tempC(5)=ratio(7)
c       C0245
c        tempC(5)=ratio(9)
c       C0345
c        tempC(5)=ratio(10)
       If(ratio(7).le.tempC(5)) tempC(5)=ratio(7)
       If(ratio(9).le.tempC(5)) tempC(5)=ratio(9)
       If(ratio(10).le.tempC(5)) tempC(5)=ratio(10)
       lossC(5)=tempC(5)*tempC(5)
c       D02345
       call dt4(p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,ratioD(5))
       lossD(5)=(ratioD(5)*ratioD(5)*ratioD(5))*lossC(5)
       If(lossD(5).le.Cut)  bad=.True.
       If (bad) Then
c       print*, "k", k
c       print*, "lossD(5)", lossD(5)
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else
c       E012345
       tempD(1)=lossD(1)
       If(lossD(2).le.tempD(1)) tempD(1)=lossD(2)
       If(lossD(3).le.tempD(1)) tempD(1)=lossD(3)
       If(lossD(4).le.tempD(1)) tempD(1)=lossD(4)
       If(lossD(5).le.tempD(1)) tempD(1)=lossD(5)
c       call dt5(p1sq,p2sq,p3sq,p4sq,p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,ratioE(1))
c       lossE(1)=(ratioE(1)*ratioE(1)*ratioE(1)*ratioE(1))*tempD(1)
c       If(lossE(1).le.Cut)  bad=.True.
       If (bad) Then
c       print*, "k", k
c       print*, "lossE(1)", lossE(1)
c       print*,"tempD(1)",tempD(1)
c       print*, "ratioE(1)",ratioE(1)
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else
c       C123
        tempC(6)=ratio(1)    
c       C0126
       call dt3(p1sq,s16,p6sq,ratio(11))
c       C0136
       call dt3(s12,s345,p6sq,ratio(12))
c       C0236
       call dt3(p2sq,s345,s16,ratio(13))
       If(ratio(11).le.tempC(6)) tempC(6)=ratio(11)
       If(ratio(12).le.tempC(6)) tempC(6)=ratio(12)
       If(ratio(13).le.tempC(6)) tempC(6)=ratio(13)
       lossC(6)=tempC(6)*tempC(6)
c       D01234
c       D01236
       call dt4(p1sq,p2sq,s345,p1p2,p1p3+p1p4+p1p5,p2p3+p2p4+p2p5,ratioD(6))
       lossD(6)=(ratioD(6)*ratioD(6)*ratioD(6))*lossC(6)
       If(lossD(6).le.Cut)  bad=.True.
       If (bad) Then
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else

c       C0124
        tempC(7)=ratio(2)  
c       C0126
c       call dt3(p1sq,s16,p6sq,ratio(11))
c       C0146
        call dt3(s123,s45,p6sq,ratio(14))
c       C0246
        call dt3(s23,s45,s16,ratio(15))
       If(ratio(11).le.tempC(7)) tempC(7)=ratio(11)
       If(ratio(14).le.tempC(7)) tempC(7)=ratio(14)
       If(ratio(15).le.tempC(7)) tempC(7)=ratio(15)
       lossC(7)=tempC(7)*tempC(7)
c       D01246
       call dt4(p1sq,s23,s45,p1p2+p1p3,p1p4+p1p5,p2p4+p2p5+p3p4+p3p5,ratioD(7))
       lossD(7)=(ratioD(7)*ratioD(7)*ratioD(7))*lossC(7)
       If(lossD(7).le.Cut)  bad=.True.
       If (bad) Then 
c       print*, "k", k
c       print*, "lossD(7)", lossD(7)
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else
c       C0134
        tempC(8)=ratio(3)
c       C0136
c       call dt3(s12,s345,p6sq,ratio(12))
c       C0146
c        call dt3(s123,s45,p6sq,ratio(14))
c       C0346
       call dt3(p3sq,s45,s345,ratio(16))
       If(ratio(12).le.tempC(8)) tempC(8)=ratio(12)
       If(ratio(14).le.tempC(8)) tempC(8)=ratio(14)
       If(ratio(16).le.tempC(8)) tempC(8)=ratio(16)
       lossC(8)=tempC(8)*tempC(8)
c       D01346
       call dt4(s12,p3sq,s45,p1p3+p2p3,p1p4+p1p5+p2p4+p2p5,p3p4+p3p5,ratioD(8))
       lossD(8)=(ratioD(8)*ratioD(8)*ratioD(8))*lossC(8)
       If(lossD(8).le.Cut)  bad=.True.
       If (bad) Then
c       print*, "k", k
c       print*, "lossD(8)", lossD(8)
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else
c       C0234
        tempC(9)=ratio(4)
c       C0236
c       call dt3(p2sq,s345,s16,ratio(13))
c       C0246
c        call dt3(s23,s45,s16,ratio(15))
c       C0346
c       call dt3(p3sq,s45,s345,ratio(16))
       If(ratio(13).le.tempC(9)) tempC(9)=ratio(13)
       If(ratio(15).le.tempC(9)) tempC(9)=ratio(15)
       If(ratio(16).le.tempC(9)) tempC(9)=ratio(16)
       lossC(9)=tempC(9)*tempC(9)
c       D02346
       call dt4(p2sq,p3sq,s45,p2p3,p2p4+p2p5,p3p4+p3p5,ratioD(9))
       lossD(9)=(ratioD(9)*ratioD(9)*ratioD(9))*lossC(9)
       If(lossD(9).le.Cut)  bad=.True.
       If (bad) Then
c       print*, "k", k
c       print*, "lossD(9)", lossD(9)
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else
c       E012346
       tempD(2)=lossD(1)
       If(lossD(6).le.tempD(2)) tempD(2)=lossD(6)
       If(lossD(7).le.tempD(2)) tempD(2)=lossD(7)
       If(lossD(8).le.tempD(2)) tempD(2)=lossD(8)
       If(lossD(9).le.tempD(2)) tempD(2)=lossD(9)
c       call dt5(p1sq,p2sq,p3sq,s45,p1p2,p1p3,p1p4+p1p5,p2p3,p2p4+p2p5,p3p4+p3p5,ratioE(2))
c       lossE(2)=(ratioE(2)*ratioE(2)*ratioE(2)*ratioE(2))*tempD(2)
c       If(lossE(2).le.Cut)  bad=.True.
       If (bad) Then
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else
c       C0125
        tempC(10)=ratio(5)
c       C0126
c       call dt3(p1sq,s16,p6sq,ratio(11))
c       C0156
       call dt3(s56,p5sq,p6sq,ratio(17))
c       C0256
       call dt3(s234,p5sq,s16,ratio(18))
       If(ratio(11).le.tempC(10)) tempC(10)=ratio(11)
       If(ratio(17).le.tempC(10)) tempC(10)=ratio(17)
       If(ratio(18).le.tempC(10)) tempC(10)=ratio(18)
       lossC(10)=tempC(10)*tempC(10)
c       D01256
       call dt4(p1sq,s234,p5sq,p1p2+p1p3+p1p4,p1p5,p2p5+p3p5+p4p5,ratioD(10))
       lossD(10)=(ratioD(10)*ratioD(10)*ratioD(10))*lossC(10)
       If(lossD(10).le.Cut)  bad=.True.
       If (bad) Then
c       print*, "k", k
c       print*, "lossD(10)", lossD(10)
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else
c       C0135
        tempC(11)=ratio(6)
c       C0136
c       call dt3(s12,s345,p6sq,ratio(12))
c       C0156
c       call dt3(s56,p5sq,p6sq,ratio(17))
c       C0356
       call dt3(s34,p5sq,s345,ratio(19))
       If(ratio(12).le.tempC(11)) tempC(11)=ratio(12)
       If(ratio(17).le.tempC(11)) tempC(11)=ratio(17)
       If(ratio(19).le.tempC(11)) tempC(11)=ratio(19)
       lossC(11)=tempC(11)*tempC(11)
c       D01356
       call dt4(s12,s34,p5sq,p1p3+p1p4+p2p3+p2p4,p1p5+p2p5,p3p5+p4p5,ratioD(11))
       lossD(11)=(ratioD(11)*ratioD(11)*ratioD(11))*lossC(11)
       If(lossD(11).le.Cut)  bad=.True.
       If (bad) Then
c       print*, "k", k
c       print*, "lossD(11)", lossD(11)
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else
c       C0235
        tempC(12)=ratio(7)
c       C0236
c       call dt3(p2sq,s345,s16,ratio(13))
c       C0256
       call dt3(s234,p5sq,s16,ratio(18))
c       C0356
       call dt3(s34,p5sq,s345,ratio(19))
       If(ratio(13).le.tempC(12)) tempC(12)=ratio(13)
       If(ratio(18).le.tempC(12)) tempC(12)=ratio(18)
       If(ratio(19).le.tempC(12)) tempC(12)=ratio(19)
       lossC(12)=tempC(12)*tempC(12)
c       D02356
       call dt4(p2sq,s34,p5sq,p2p3+p2p4,p2p5,p3p5+p4p5,ratioD(12))
       lossD(12)=(ratioD(12)*ratioD(12)*ratioD(12))*lossC(12)
       If(lossD(12).le.Cut)  bad=.True.
       If (bad) Then
c       print*, "k", k
c       print*, "lossD(12)", lossD(12)
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else
c       E012356
       tempD(3)=lossD(2)
       If(lossD(6).le.tempD(3)) tempD(3)=lossD(6)
       If(lossD(10).le.tempD(3)) tempD(3)=lossD(10)
       If(lossD(11).le.tempD(3)) tempD(3)=lossD(11)
       If(lossD(12).le.tempD(3)) tempD(3)=lossD(12)
c       call dt5(p1sq,p2sq,s34,p5sq,p1p2,p1p3+p1p4,p1p5,p2p3+p2p4,p2p5,p3p5+p4p5,ratioE(3))
c       lossE(3)=(ratioE(3)*ratioE(3)*ratioE(3)*ratioE(3))*tempD(3)
c       If(lossE(3).le.Cut)  bad=.True.
       If (bad) Then
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else
c       C0145
      tempC(13)=ratio(8)
c       C0146
c        call dt3(s123,s45,p6sq,ratio(14))
c       C0156
c       call dt3(s56,p5sq,p6sq,ratio(17))
c       C0456
       call dt3(p4sq,p5sq,s45,ratio(20))
       If(ratio(14).le.tempC(13)) tempC(13)=ratio(14)
       If(ratio(17).le.tempC(13)) tempC(13)=ratio(17)
       If(ratio(20).le.tempC(13)) tempC(13)=ratio(20)
       lossC(13)=tempC(13)*tempC(13)
c       D01456
       call dt4(s123,p4sq,p5sq,p1p4+p2p4+p3p4,p1p5+p2p5+p3p5,p4p5,ratioD(13))
       lossD(13)=(ratioD(13)*ratioD(13)*ratioD(13))*lossC(13)
       If(lossD(13).le.Cut)  bad=.True.
       If (bad) Then
c       print*, "k", k
c       print*, "lossD(13)", lossD(13)
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else
c      C0245
       tempC(14)=9
c      C0246
c       tempC(14)=15
c      C0256
c       tempC(14)=18
c      C0456
c       tempC(14)=20
       If(ratio(15).le.tempC(14)) tempC(14)=ratio(15)
       If(ratio(18).le.tempC(14)) tempC(14)=ratio(18)
       If(ratio(20).le.tempC(14)) tempC(14)=ratio(20)
       lossC(14)=tempC(14)*tempC(14)
c       D02456
       call dt4(s23,p4sq,p5sq,p2p4+p3p4,p2p5+p3p5,p4p5,ratioD(14))
       lossD(14)=(ratioD(14)*ratioD(14)*ratioD(14))*lossC(14)
       If(lossD(14).le.Cut)  bad=.True.
       If (bad) Then
c       print*, "k", k
c       print*, "lossD(14)", lossD(14)
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else
c       E012456
       tempD(4)=lossD(3)   
       If(lossD(7).le.tempD(3)) tempD(3)=lossD(7)
       If(lossD(10).le.tempD(3)) tempD(3)=lossD(10)
       If(lossD(13).le.tempD(3)) tempD(3)=lossD(13)
       If(lossD(14).le.tempD(3)) tempD(3)=lossD(14)
c       call dt5(p1sq,s23,p4sq,p5sq,p1p2+p1p3,p1p4,p1p5,p2p4+p3p4,p2p5+p3p5,p4p5,ratioE(4))
c       lossE(4)=(ratioE(4)*ratioE(4)*ratioE(4)*ratioE(4))*tempD(4)
c       If(lossE(4).le.Cut)  bad=.True.
       If (bad) Then
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else
c      C0345
       tempC(15)=ratio(10)   
c      C0346
c       tempC(15)=ratio(16)
c      C0356
c       tempC(15)=ratio(19)
c      C0456
c     tempC(15)=ratio(20)
       If(ratio(16).le.tempC(15)) tempC(15)=ratio(16)
       If(ratio(19).le.tempC(15)) tempC(15)=ratio(19)
       If(ratio(20).le.tempC(15)) tempC(15)=ratio(20)
       lossC(15)=tempC(15)*tempC(15)
c       D03456
       call dt4(p3sq,p4sq,p5sq,p3p4,p3p5,p4p5,ratioD(15))
       lossD(15)=(ratioD(15)*ratioD(15)*ratioD(15))*lossC(15)
       If(lossD(15).le.Cut)  bad=.True.
       If (bad) Then
c       print*, "k", k
c       print*, "lossD(15)", lossD(15)
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else
c       E013456
       tempD(5)=lossD(4)   
       If(lossD(8).le.tempD(5)) tempD(5)=lossD(8)
       If(lossD(11).le.tempD(5)) tempD(5)=lossD(11)
       If(lossD(13).le.tempD(5)) tempD(5)=lossD(13)
       If(lossD(15).le.tempD(5)) tempD(5)=lossD(15)
c       call dt5(s12,p3sq,p4sq,p5sq,p1p3+p2p3,p1p4+p2p4,p1p5+p2p5,p3p4,p3p5,p4p5,ratioE(5))
c       lossE(5)=(ratioE(5)*ratioE(5)*ratioE(5)*ratioE(5))*tempD(5)
c       If(lossE(5).le.Cut)  bad=.True.
       If (bad) Then
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else
c       E023456
       tempD(6)=lossD(5)   
       If(lossD(9).le.tempD(6)) tempD(6)=lossD(9)
       If(lossD(12).le.tempD(6)) tempD(6)=lossD(12)
       If(lossD(14).le.tempD(6)) tempD(6)=lossD(14)
       If(lossD(15).le.tempD(6)) tempD(6)=lossD(15)
c       call dt5(p2sq,p3sq,p4sq,p5sq,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,ratioE(6))
c       lossE(6)=(ratioE(6)*ratioE(6)*ratioE(6)*ratioE(6))*tempD(6)
c       If(lossE(6).le.Cut)  bad=.True.
       If (bad) Then
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       else
c       F0123456
       tempE(1)=lossE(1)
       If(lossE(2).le.tempE(1)) tempE(1)=lossE(2)
       If(lossE(3).le.tempE(1)) tempE(1)=lossE(3)
       If(lossE(4).le.tempE(1)) tempE(1)=lossE(4)
       If(lossE(5).le.tempE(1)) tempE(1)=lossE(5)
       If(lossE(6).le.tempE(1)) tempE(1)=lossE(6)
c       call dt6(p1sq,p2sq,p3sq,p4sq,p5sq,p1p2,p1p3,p1p4,p1p5,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,ratioF(1))
c       lossF(1)=(ratioF(1)*ratioF(1)*ratioF(1)*ratioF(1)*ratioF(1))*tempE(1)
c       If(lossF(1).le.Cut)  bad=.True.
       If (bad) Then
       countbad=countbad+1
       seedbad(countbad)=k
       Return
       Endif
       Endif
       Endif
       Endif
       Endif
       Endif
       Endif
       Endif
       Endif
       Endif
       Endif
       Endif
       Endif
       Endif
       Endif
       Endif
       Endif
       Endif
       Endif
       Endif
       Endif
       Endif
c       print*, 'here',bad
       return
       end


      subroutine Calc_Det_up_D(p1,p2,p3,p4,k,badF)
c Author: Francisco Campanario 
c DAte: 14.08.2009
      Implicit none 
      real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      real*8 p1sq,p2sq,p3sq,p4sq
      real*8 s12,s23
      real*8 p1p2,p1p3,p1p4
      real*8 p2p3,p2p4
      real*8 p3p4
      real*8 ratio(4),tempC(1),lossC(1)
      real*8 ratioD(1),lossD(1)
      logical bad,badF
      Integer k
      real*8 dotrr,Cut
      parameter(Cut=1d-21)
      Common bad
     
       p1sq = dotrr(p1,p1)
       p1p2 = dotrr(p1,p2)
       p1p3 = dotrr(p1,p3)
       p1p4 = dotrr(p1,p4)


       p2sq = dotrr(p2,p2)
       p2p3 = dotrr(p2,p3)
       p2p4 = dotrr(p2,p4)


       p3sq = dotrr(p3,p3)
       p3p4 = dotrr(p3,p4)

       p4sq = dotrr(p4,p4)

       s12 = (p1sq +p2sq+ 2*p1p2) 
       s23 = (p2sq +p3sq+ 2*p2p3) 
 
       bad=.False.

c       C0123
       call dt3(p1sq,p2sq,s12,ratio(1))
c       C0124
       call dt3(p1sq,s23,p4sq,ratio(2))
c       C0134
       call dt3(s12,p3sq,p4sq,ratio(3))
c       C0234
       call dt3(p2sq,p3sq,s23,ratio(4))

       tempC(1)=ratio(1)
       If(ratio(2).le.tempC(1)) tempC(1)=ratio(2)
       If(ratio(3).le.tempC(1)) tempC(1)=ratio(3)
       If(ratio(4).le.tempC(1)) tempC(1)=ratio(4)
       lossC(1)=tempC(1)*tempC(1)
c       D01234
       call dt4(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,ratioD(1))
       lossD(1)=(ratioD(1)*ratioD(1)*ratioD(1))*lossC(1)

       If(lossD(1).le.Cut)  bad=.True.
       
       badF=bad

       end


      subroutine Calc_Det_up_E(p1,p2,p3,p4,p5,countbad)
c Francisco Campanario
c Date:03.03.2010
c Compute gram and Caeley Determinants and write histograms
c to analized their behaviour.
C Need to modified main code to write down the histograms 
c as well as define a common variable. Possible in global.inc
      Implicit none 
      real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
      real*8 p1sq,p2sq,p3sq,p4sq,p5sq
      real*8 s12,s23,s34,s45,s15
      real*8 p1p2,p1p3,p1p4,p1p5
      real*8 p2p3,p2p4,p2p5
      real*8 p3p4,p3p5
      real*8 p4p5
      real*8 ratio(10),tempC(10),lossC(10)
      real*8 ratioD(5),lossD(5),tempD(1)
      real*8 ratioE(1),lossE(1)
      logical bad
      Integer countbad
      real*8 dotrr,Cut
      parameter(Cut=1d-21)
      Common bad
c histograms
      Integer ihC(0:1200),ihCM(0:120),ihCMT(0:1200),ihCL(0:1200),ihCLT(0:1200)
      Integer ihD(0:1200),ihDM(0:1200),ihDL(0:1200),ihDLT(0:1200)
      Integer ihE(0:1200),ihEL(0:1200)
      Integer Small,i
      Real*8 temph
      Integer binsize
      Real*8 offset
      common/histDet/ihC,ihCM,ihCMT,ihCL,ihCLT,
     &               ihD,ihDM,ihDL,ihDLT,
     &               ihE,ihEL,
     &                binsize,offset
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
      binsize=4
      offset=1000d0
       p1sq = dotrr(p1,p1)
       p1p2 = dotrr(p1,p2)
       p1p3 = dotrr(p1,p3)
       p1p4 = dotrr(p1,p4)
       p1p5 = dotrr(p1,p5)
       p2sq = dotrr(p2,p2)
       p2p3 = dotrr(p2,p3)
       p2p4 = dotrr(p2,p4)
       p2p5 = dotrr(p2,p5)
       p3sq = dotrr(p3,p3)
       p3p4 = dotrr(p3,p4)
       p3p5 = dotrr(p3,p5)
       p4sq = dotrr(p4,p4)
       p4p5 = dotrr(p4,p5)
       p5sq = dotrr(p5,p5)
       s12 = (p1sq +p2sq+ 2*p1p2) 
       s23 = (p2sq +p3sq+ 2*p2p3) 
       s34 = (p3sq +p4sq+ 2*p3p4) 
       s45 = (p4sq +p5sq+ 2*p4p5) 
       s15 = (p1sq +p5sq+ 2*p1p5) 
       bad=.False.
c       C0123
       call dt3(p1sq,p2sq,s12,ratio(1))
c       C0124
       call dt3(p1sq,s23,s45,ratio(2))
c       C0134
       call dt3(s12,p3sq,s45,ratio(3))
c       C0234
       call dt3(p2sq,p3sq,s23,ratio(4))
       tempC(1)=ratio(1)
       If(ratio(2).le.tempC(1)) tempC(1)=ratio(2)
       If(ratio(3).le.tempC(1)) tempC(1)=ratio(3)
       If(ratio(4).le.tempC(1)) tempC(1)=ratio(4)
       lossC(1)=tempC(1)*tempC(1)
c       D01234
       call dt4(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,ratioD(1))
       lossD(1)=(ratioD(1)*ratioD(1)*ratioD(1))*lossC(1)
       If(lossD(1).le.Cut)  bad=.True.
       If (bad) Then
       countbad=countbad+1
c       seedbad(countbad)=k
c       Return
c       else
       endif
c       C123
        tempC(2)=ratio(1)      
c       C0125
       call dt3(p1sq,s15,p5sq,ratio(5))
c       C0135
       call dt3(s12,s34,p5sq,ratio(6))
c       C0235
       call dt3(p2sq,s34,s15,ratio(7))
       If(ratio(5).le.tempC(2)) tempC(2)=ratio(5)
       If(ratio(6).le.tempC(2)) tempC(2)=ratio(6)
       If(ratio(7).le.tempC(2)) tempC(2)=ratio(7)
       lossC(2)=tempC(2)*tempC(2)
c       D01235
       call dt4(p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p3+p2p4,ratioD(2))
       lossD(2)=(ratioD(2)*ratioD(2)*ratioD(2))*lossC(2)
       If(lossD(2).le.Cut)  bad=.True.
       If (bad) Then
       countbad=countbad+1
c       seedbad(countbad)=k
c       Return
c       else
       Endif
c       C0124
        tempC(3)=ratio(2)      
c       C0125
c        tempC(3)=ratio(5)      
c       C0145
       call dt3(s45,p4sq,p5sq,ratio(8))
c       C0245
       call dt3(s23,p4sq,s15,ratio(9))
       If(ratio(5).le.tempC(3)) tempC(3)=ratio(5)
       If(ratio(8).le.tempC(3)) tempC(3)=ratio(8)
       If(ratio(9).le.tempC(3)) tempC(3)=ratio(9)
       lossC(3)=tempC(3)*tempC(3)
c       D01245
       call dt4(p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p4+p3p4,ratioD(3))
       lossD(3)=(ratioD(3)*ratioD(3)*ratioD(3))*lossC(3)
       If(lossD(3).le.Cut)  bad=.True.
       If (bad) Then
       countbad=countbad+1
c       seedbad(countbad)=k
c       Return
c       else
       endif
c       C0134
        tempC(4)=ratio(3)
c       C0135
c        tempC(4)=ratio(6)
c       C0145
c        tempC(4)=ratio(8)
c       C0345
        call dt3(p3sq,p4sq,s34,ratio(10))
       If(ratio(6).le.tempC(4)) tempC(4)=ratio(6)
       If(ratio(8).le.tempC(4)) tempC(4)=ratio(8)
       If(ratio(10).le.tempC(4)) tempC(4)=ratio(10)
       lossC(4)=tempC(4)*tempC(4)
c       D01345
       call dt4(s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4,p3p4,ratioD(4))
       lossD(4)=(ratioD(4)*ratioD(4)*ratioD(4))*lossC(4)
       If(lossD(4).le.Cut)  bad=.True.
       If (bad) Then
       countbad=countbad+1
c       seedbad(countbad)=k
c       Return
c       else
       endif
c       C0234
        tempC(5)=ratio(4)
c       C0235
c        tempC(5)=ratio(7)
c       C0245
c        tempC(5)=ratio(9)
c       C0345
c        tempC(5)=ratio(10)
       If(ratio(7).le.tempC(5)) tempC(5)=ratio(7)
       If(ratio(9).le.tempC(5)) tempC(5)=ratio(9)
       If(ratio(10).le.tempC(5)) tempC(5)=ratio(10)
       lossC(5)=tempC(5)*tempC(5)
c       D02345
       call dt4(p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,ratioD(5))
       lossD(5)=(ratioD(5)*ratioD(5)*ratioD(5))*lossC(5)
       If(lossD(5).le.Cut)  bad=.True.
       If (bad) Then
       countbad=countbad+1
c       seedbad(countbad)=k
c       Return
c       else
       endif
c       E012345
       tempD(1)=lossD(1)
       If(lossD(2).le.tempD(1)) tempD(1)=lossD(2)
       If(lossD(3).le.tempD(1)) tempD(1)=lossD(3)
       If(lossD(4).le.tempD(1)) tempD(1)=lossD(4)
       If(lossD(5).le.tempD(1)) tempD(1)=lossD(5)
       call dt5(p1sq,p2sq,p3sq,p4sq,p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,ratioE(1))
       lossE(1)=(ratioE(1)*ratioE(1)*ratioE(1)*ratioE(1))*tempD(1)
       If(lossE(1).le.Cut)  bad=.True.
       If (bad) Then
       countbad=countbad+1
c       seedbad(countbad)=k
c       Return
       ENDIF
CCCCCCCCCCCCCCCCCCCCCc
c   histograms
Cij
       do i=1,10
       temph=Log10(ratio(i))+offset
       Small=temph*binsize
       ihC(Small)=ihC(Small)+1
       enddo
       do i=1,5
        temph=Log10(tempC(i))+offset
        Small=temph*binsize
        ihCM(Small)=ihCM(Small)+1
       enddo
       do i=1,1
        tempC(1)=Min(tempC(1),tempC(2),tempC(3),tempC(4),tempC(5))
        temph=Log10(tempC(i))+offset
        Small=temph*binsize
        ihCMT(Small)=ihCMT(Small)+1
       enddo
       do i=1,5
       temph=Log10(lossC(i))+offset
       Small=temph*binsize
       ihCL(Small)=ihCL(Small)+1
       enddo
       do i=1,1
       lossC(1)=Min(lossC(1),lossC(2),lossC(3),lossC(4),lossC(5))
       temph=Log10(lossC(1))+offset
       Small=temph*binsize
       ihCLT(Small)=ihCLT(Small)+1
       enddo
c   Dij
       do i=1,5
       temph=Log10(ratioD(i))+offset
       Small=temph*binsize
       ihD(Small)=ihD(Small)+1
       enddo
       do i=1,1
        tempD(1)=Min(ratioD(1),ratioD(2),ratioD(3)
     &  ,ratioD(4),ratioD(5))
        temph=Log10(tempD(i))+offset
        Small=temph*binsize
        ihDM(Small)=ihDM(Small)+1
       enddo
       do i=1,5
       temph=Log10(lossD(i))+offset
       Small=temph*binsize
       ihDL(Small)=ihDL(Small)+1
       enddo
       do i=1,1
       lossD(1)=Min(lossD(1),lossD(2),lossD(3),lossD(4),lossD(5))
       temph=Log10(lossD(1))+offset
       Small=temph*binsize
       ihDLT(Small)=ihDLT(Small)+1
       enddo
c Eij
       do i=1,1
       temph=Log10(ratioE(i))+offset
       Small=temph*binsize
       ihE(Small)=ihE(Small)+1
       enddo
       do i=1,1
       temph=Log10(lossE(i))+offset
       Small=temph*binsize
       ihEL(Small)=ihEL(Small)+1
       enddo
c       Endif
c       Endif
c       Endif
c       Endif
c       Endif
c       print*, 'here',bad
       return
       end



