       subroutine gluonboxHiggsEven_T(m,p1,p2,p3,p4,mup1,mup2,mup3,resul
     -   t)

       IMPLICIT NONE
      Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      Complex*16 mup1(0:3),mup2(0:3),mup3(0:3)
      real*8 m
      complex*16 result
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex*16 mup1temp(0:3),mup2temp(0:3),mup3temp(0:3)
      complex*16 box_Ward1,box_Ward2,box_Ward3
      real*8 pp12(0:3),pp23(0:3),pp34(0:3),pp14(0:3)
      real*8 gauge1,gauge2,gauge3
      real*8 gauge
      complex*16 ver1_for_Ward1, ver1_for_Ward2
      complex*16 ver2_for_ward2, ver2_for_ward3
      complex*16 ver3_for_ward3, ver3_for_ward1
ccccccccccccccccccccccccccccccccccccccccccccccc
       logical QUAD,Dble1
       integer jtemp
 
        logical HexSingular(2),pensingular(2),boxSingular(2),verSingular(2)
        common/hexSin/HexSingular,pensingular,boxSingular,verSingular
       real*8 gaugelimit 
       common/gauge/gaugelimit
      common/dbleswitch/dble1,quad



      call gluonboxHiggsEven(m,p1,p2,p3,p4,mup1,mup2,mup3,resul
     -   t,1)

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
            pp12(jtemp)=p1(jtemp)+p2(jtemp)
            pp23(jtemp)=p2(jtemp)+p3(jtemp)
            pp34(jtemp)=p3(jtemp)+p4(jtemp)
            pp14(jtemp)=-(p2(jtemp)+p3(jtemp))
         enddo   
c Ward 1
      call gluonboxHiggsEven(m,p1,p2,p3,p4,
     &   mup1temp,mup2,mup3,box_Ward1,-1)
c Ward 2
      call gluonboxHiggsEven(m,p1,p2,p3,p4,
     &   mup1,mup2temp,mup3,box_Ward2,-1)
c Ward 3
      call gluonboxHiggsEven(m,p1,p2,p3,p4,
     &   mup1,mup2,mup3temp,box_Ward3,-1)
cccccccccccccccccccccccccccccccccccccccccccc 
c           verteces
cccccccccccccccccccccccccccccccccccccccc
c  pp12
      call gluonverHiggsEven(m,pp12,p3,p4,
     &      mup2,mup3,1,ver1_for_Ward1)
      call gluonverHiggsEven(m,pp12,p3,p4,
     &      mup1,mup3,-1,ver1_for_Ward2)
c  pp23
      call gluonverHiggsEven(m,p1,pp23,p4,
     &      mup1,mup3,1,ver2_for_Ward2)
      call gluonverHiggsEven(m,p1,pp23,p4,
     &      mup1,mup2,-1,ver2_for_Ward3)
c pp34
      call gluonverHiggsEven(m,p1,p2,pp34,
     &      mup1,mup2,1,ver3_for_Ward3)
c pp14
      call gluonverHiggsEven(m,p2,p3,pp14,
     &      mup2,mup3,1,ver3_for_Ward1)

c gauge1
         if(abs(box_Ward1).gt.1d-2) then
         gauge1=abs(3*m*(ver1_for_Ward1-ver3_for_Ward1)/(Box_Ward1)-1d0)
         else
         gauge1=abs(3*m*(ver1_for_Ward1-ver3_for_Ward1) - (Box_Ward1))
         endif
c gauge2
         if(abs(Box_Ward2).gt.1d-2) then
         gauge2=abs(3*m*(ver2_for_Ward2-ver1_for_Ward2)/(Box_Ward2)-1d0)
         else
         gauge2=abs(3*m*(ver2_for_Ward2-ver1_for_Ward2) - (Box_Ward2))
         endif
c gauge3
         if(abs(Box_Ward3).gt.1d-2) then
         gauge3=abs(3*m*(ver3_for_Ward3-ver2_for_Ward3)/(Box_Ward3)-1d0)
         else
         gauge3=abs(3*m*(ver3_for_Ward3-ver2_for_Ward3) - (Box_Ward3))
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

 
cccccccccccccccccccccccccccccccccccccccccc
         gauge=Max(gauge1,gauge2,gauge3)


         if(gauge+1.eq.gauge) then
            gauge=1d99
         endif

         if(gauge.ne.gauge) then
            gauge=1d99
         endif

c         print*, "gauge1",gauge1
c         print*, "gauge2",gauge2
c         print*, "gauge3",gauge3
         if (gauge.lt.gaugelimit) then
c             PenSingular=.false.
c             print*,"gauge:true:B", result,gauge,boxSingular
             else
             boxSingular(1)=.true.
c            print*,"gauge:false:B",result, gauge,boxSingular(1)
          endif

         

    
      end subroutine
