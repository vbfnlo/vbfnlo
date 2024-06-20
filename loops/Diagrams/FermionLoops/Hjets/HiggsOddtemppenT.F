      subroutine gluonpenHiggsOdd_T(m,p1,p2,p3,p4,p5,mup1,mup2
     &                               ,mup3,mup4,result)

       IMPLICIT NONE
      Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
      Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3)
      real*8 m
      complex*16 result
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex*16 mup1temp(0:3),mup2temp(0:3),mup3temp(0:3),mup4temp(0:3)
      complex*16 pen_Ward1,pen_Ward2,pen_Ward3,pen_Ward4
      real*8 pp12(0:3),pp23(0:3),pp34(0:3),pp45(0:3),pp15(0:3)
      real*8 gauge1,gauge2,gauge3,gauge4
      real*8 gauge
      complex*16 box1_for_Ward1, box1_for_Ward2
      complex*16 box2_for_ward2, box2_for_ward3
      complex*16 box3_for_ward3, box3_for_ward4
      complex*16 box4_for_ward4, box4_for_ward1
ccccccccccccccccccccccccccccccccccccccccccccccc
       logical QUAD,Dble1
       integer jtemp
 
        logical HexSingular(2),pensingular(2),boxSingular(2),verSingular(2)
        common/hexSin/HexSingular,pensingular,boxSingular,verSingular
       real*8 gaugelimit 
       common/gauge/gaugelimit
      common/dbleswitch/dble1,quad



      call gluonpenHiggsOdd(m,p1,p2,p3,p4,p5,mup1,mup2,mup3,mup4,
     &   result,1)

c start Ward identity
c The momenta of the Hexline is going to be the same. 
c                       --->  comp=-1 to avoid redundant calculation
c First Gauge test Eps_i -->   p_i
c +  auxiliary momenta
c         print*, "i1--i5", i1,i2,i3,i4,i5
         do jtemp=0,3
            mup1temp(jtemp)= p1(jtemp)
            mup2temp(jtemp)= p2(jtemp)
            mup3temp(jtemp)= p3(jtemp)
            mup4temp(jtemp)= p4(jtemp)
            
            pp12(jtemp)=p1(jtemp)+p2(jtemp)
            pp23(jtemp)=p2(jtemp)+p3(jtemp)
            pp34(jtemp)=p3(jtemp)+p4(jtemp)
            pp45(jtemp)=p4(jtemp)+p5(jtemp)
            pp15(jtemp)=-(p2(jtemp)+p3(jtemp)+p4(jtemp))
         enddo   
c Ward 1
      call gluonpenHiggsOdd(m,p1,p2,p3,p4,p5,
     &   mup1temp,mup2,mup3,mup4,pen_Ward1,-1)
c Ward 2
      call gluonpenHiggsOdd(m,p1,p2,p3,p4,p5,
     &   mup1,mup2temp,mup3,mup4,pen_Ward2,-1)
c Ward 3
      call gluonpenHiggsOdd(m,p1,p2,p3,p4,p5,
     &     mup1,mup2,mup3temp,mup4,pen_Ward3,-1)
c Ward 4
      call gluonpenHiggsOdd(m,p1,p2,p3,p4,p5,
     &   mup1,mup2,mup3,mup4temp,pen_Ward4,-1)
cccccccccccccccccccccccccccccccccccccccccccc 
c           verteces
cccccccccccccccccccccccccccccccccccccccc
c  pp12
      call gluonboxHiggsOdd(m,pp12,p3,p4,p5,
     &      mup2,mup3,mup4,box1_for_Ward1,1)
      call gluonboxHiggsOdd(m,pp12,p3,p4,p5,
     &      mup1,mup3,mup4,box1_for_Ward2,-1)
c  pp23
      call gluonboxHiggsOdd(m,p1,pp23,p4,p5,
     &      mup1,mup3,mup4,box2_for_Ward2,1)
      call gluonboxHiggsOdd(m,p1,pp23,p4,p5,
     &      mup1,mup2,mup4,box2_for_Ward3,-1)
c pp34
      call gluonboxHiggsOdd(m,p1,p2,pp34,p5,
     &     mup1,mup2,mup4,box3_for_Ward3,1)
      call gluonboxHiggsOdd(m,p1,p2,pp34,p5,
     &     mup1,mup2,mup3,box3_for_Ward4,-1)
c pp45
      call gluonboxHiggsOdd(m,p1,p2,p3,pp45,
     &     mup1,mup2,mup3,box4_for_Ward4,1)
c pp15 
      call gluonboxHiggsOdd(m,p2,p3,p4,pp15,
     &     mup2,mup3,mup4,box4_for_Ward1,1)


      
      
c gauge1
         if(abs(Pen_Ward1).gt.1d-2) then
         gauge1=abs(1d0*(box1_for_Ward1-box4_for_Ward1)/(Pen_Ward1)-1d0)
         else
         gauge1=abs(1d0*(box1_for_Ward1-box4_for_Ward1) - (Pen_Ward1))
         endif
c gauge2
         if(abs(Pen_Ward2).gt.1d-2) then
         gauge2=abs(1d0*(box2_for_Ward2-box1_for_Ward2)/(Pen_Ward2)-1d0)
         else
         gauge2=abs(1d0*(box2_for_Ward2-box1_for_Ward2) - (Pen_Ward2))
         endif
c gauge3
         if(abs(Pen_Ward3).gt.1d-2) then
         gauge3=abs(1d0*(box3_for_Ward3-box2_for_Ward3)/(Pen_Ward3)-1d0)
         else
         gauge3=abs(1d0*(box3_for_Ward3-box2_for_Ward3) - (Pen_Ward3))
         endif
c gauge4
         if(abs(Pen_Ward4).gt.1d-2) then
         gauge4=abs(1d0*(box4_for_Ward4-box3_for_Ward4)/(Pen_Ward4)-1d0)
         else
         gauge4=abs(1d0*(box4_for_Ward4-box3_for_Ward4) - (Pen_Ward4))
         endif


         if(gauge1+1.eq.gauge1) then
            gauge1=1d99
         endif
         if(gauge1.ne.gauge1) then
            gauge1=1d99
         endif

         if(gauge2+1.eq.gauge2) then
            gauge2=1d99
         endif
         if(gauge2.ne.gauge2) then
            gauge2=1d99
         endif


         if(gauge3+1.eq.gauge3) then
            gauge3=1d99
         endif
         if(gauge3.ne.gauge3) then
            gauge3=1d99
         endif


         if(gauge4+1.eq.gauge4) then
            gauge4=1d99
         endif
         if(gauge4.ne.gauge4) then
            gauge4=1d99
         endif

 
cccccccccccccccccccccccccccccccccccccccccc
         gauge=Max(gauge1,gauge2,gauge3,gauge4)


         if(gauge+1.eq.gauge) then
            gauge=1d99
         endif

         if(gauge.ne.gauge) then
            gauge=1d99
         endif

c         print*, "gauge1",gauge1
c         print*, "gauge2",gauge2
c         print*, "gauge3",gauge3
         if (gauge.gt.gaugelimit) then
c             PenSingular=.false.
c             print*,"gauge:true:P", result,gauge,boxSingular,gauge1,gauge2,gauge3,gauge4
cfc             else
             penSingular(1)=.true.
cfc             print*,"gauge:false:B",result, gauge,penSingular(1),gauge1,gauge2,gauge3,gauge4
cfc
cfc
cfc             print*, "W1"
cfc             print*, box1_for_Ward1-box4_for_Ward1
cfc             print*, box1_for_Ward1,box4_for_Ward1
cfc             print*, Pen_Ward1
cfc             print*, "W2"
cfc             print*, box2_for_Ward2-box1_for_Ward2
cfc             print*, box2_for_Ward2, box1_for_Ward2
cfc             print*, Pen_Ward2
cfc             print*, "W3"
cfc             print*, box3_for_Ward3-box2_for_Ward3
cfc             print*, box3_for_Ward3, box2_for_Ward3
cfc             print*, Pen_Ward3
cfc             print*, "W4"
cfc             print*, box4_for_Ward4-box3_for_Ward4
cfc             print*, box4_for_Ward4,box3_for_Ward4
cfc             print*, Pen_Ward4,abs(1d0*(box4_for_Ward4-box3_for_Ward4)/(Pen_Ward4)-1d0),1d0*(box4_for_Ward4-box3_for_Ward4) - (Pen_Ward4)



             
             
          endif

         

    
      end subroutine
