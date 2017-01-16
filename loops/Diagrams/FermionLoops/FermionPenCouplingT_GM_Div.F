c Author: Francisco Campanario
c  Computes the Box fermion loop with generic couplings and masses
c  cpl(1)-cpl(4)  P_Right couplings
c  cpl(5)-cpl(8)  P_left couplingss

       subroutine FermionPenCoupling_G345_Div(cpl,m1,m2,m3,
     & p1,p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5,
     & musq,gauge,comp,result,resultgauge,Div)
       implicit none
       Real*8 cpl(10),m1,m2,m3
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3),mup5(0:3)
       Complex*16 mup1test(0:3),mup2test(0:3),mup3test(0:3),mup4test(0:3)
       Complex*16 mup5test(0:3)
       Complex*16 result,resultn(1,6)
       Real*8 musq,resultgauge,temp1,temp2,temp3,temp4,temp5
       complex*16 rg1(1,6),rg2(1,6),rg3(1,6),rg4(1,6),rg5(1,6)
       integer gauge,i,comp,Div

        resultgauge=1d99

       call  PenEvenG_Div(cpl,m1,m2,m3,p1,
     & p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5,musq,comp,resultn,Div)

       result=resultn(1,1)

       if (gauge.gt.0) then
       
         do  i=0,3
            mup1test(i)=p1(i)
            mup2test(i)=p2(i)
            mup3test(i)=p3(i)
            mup4test(i)=p4(i)
            mup5test(i)=p5(i)
c            p12(i)= p1(i)+p2(i)
c            p14(i)= p1(i)+p4(i)
c            p23(i)= p2(i)+p3(i)
c            p34(i)= p3(i)+p4(i)
c            p45(i)= p4(i)+p5(i)
         enddo
       
c Call for Gauge

          call PenEvenG_Div(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1test,mup2,mup3,mup4,mup5,
     &     musq,-1,rg1,Div)
          call PenEvenG_Div(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2test,mup3,mup4,mup5,
     &     musq,-1,rg2,Div)
          call PenEvenG_Div(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2,mup3test,mup4,mup5,
     &     musq,-1,rg3,Div)
          call PenEvenG_Div(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2,mup3,mup4test,mup5,
     &     musq,-1,rg4,Div)
          call PenEvenG_Div(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5test,
     &     musq,-1,rg5,Div)


c Gauge 1
          if (abs(rg1(1,1)).gt.1d-7) then
          temp1=abs(rg1(1,1)/(resultn(1,2)) -1d0)
          else
          temp1=abs(rg1(1,1)-resultn(1,2))   
          endif

c Gauge 2
         if (abs(rg2(1,1)).gt.1d-7) then
          temp2=abs(rg2(1,1)/resultn(1,3) -1d0) 
          else
          temp2=abs(rg2(1,1) -resultn(1,3)  ) 
          endif

c Gauge 3
         if (abs(rg3(1,1)).gt.1d-7) then
          temp3=abs(rg3(1,1)/resultn(1,4) -1d0 )
          else
             temp3=abs(rg3(1,1)  -resultn(1,4) )
         endif

c Gauge 4
          if (abs(rg4(1,1)).gt.1d-7) then
          temp4=abs(rg4(1,1)/resultn(1,5) -1d0)
          else
          temp4=abs(rg4(1,1) - resultn(1,5)  )
          endif

c Gauge 5
          if (abs(rg5(1,1)).gt.1d-7) then
          temp5=abs(rg5(1,1)/resultn(1,6) -1d0)
          else
          temp5=abs(rg5(1,1) - resultn(1,6))
          endif
      
cFC          print*, "rg1",rg1
cFC          print*, "rg2",rg2
cFC          print*, "rg3",rg3
cFC          print*, "rg4",rg4
cFC          print*, "temp1",temp1
cFC          print*, "temp2",temp2
cFC          print*, "temp3",temp3
cFC          print*, "temp4",temp4
cFC          print*, "result1",resultn(1)
cFC          print*, "result",resultn(2)
cFC          print*, "result3",resultn(3)
cFC          print*, "result4",resultn(4)
cFC          print*, "result5",resultn(5)
c FC %c	read(*,*)

          resultgauge=max(temp1,temp2,temp3,temp4,temp5)
        
        endif
         
        end

cFC
cFC
cFCcccccccccccccccccccc    ODD PART ccccccccccccccccccccccccccccccccccccccccccccccccccc


       subroutine FermionPenCouplingO_G345_Div(cpl,m1,m2,m3,
     & p1,p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5,
     & musq,gauge,comp,result,resultgauge,Div)
       implicit none
       Real*8 cpl(10),m1,m2,m3
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3),mup5(0:3)
       Complex*16 mup1test(0:3),mup2test(0:3),mup3test(0:3),mup4test(0:3)
       Complex*16 mup5test(0:3)
       Complex*16 result,resultn(1,6)
       Real*8 musq,resultgauge,temp1,temp2,temp3,temp4,temp5
       complex*16 rg1(1,6),rg2(1,6),rg3(1,6),rg4(1,6),rg5(1,6)
       integer gauge,i,comp,Div

        resultgauge=1d99

       call  PenOddG_Div(cpl,m1,m2,m3,p1,
     & p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5,musq,comp,resultn,Div)

       result=resultn(1,1)

       if (gauge.gt.0) then
       
         do  i=0,3
            mup1test(i)=p1(i)
            mup2test(i)=p2(i)
            mup3test(i)=p3(i)
            mup4test(i)=p4(i)
            mup5test(i)=p5(i)
c            p12(i)= p1(i)+p2(i)
c            p14(i)= p1(i)+p4(i)
c            p23(i)= p2(i)+p3(i)
c            p34(i)= p3(i)+p4(i)
         enddo
       
c Call for Gauge

          call PenOddG_Div(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1test,mup2,mup3,mup4,mup5,
     &     musq,-1,rg1,Div)
          call PenOddG_Div(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2test,mup3,mup4,mup5,
     &     musq,-1,rg2,Div)
          call PenOddG_Div(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2,mup3test,mup4,mup5,
     &     musq,-1,rg3,Div)
          call PenOddG_Div(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2,mup3,mup4test,mup5,
     &     musq,-1,rg4,Div)
          call PenOddG_Div(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5test,
     &     musq,-1,rg5,Div)

c Gauge 1
          if (abs(rg1(1,1)).gt.1d-7) then
          temp1=abs(rg1(1,1)/(resultn(1,2)) -1d0)
          else
          temp1=abs(rg1(1,1)-resultn(1,2))   
          endif

c Gauge 2
         if (abs(rg2(1,1)).gt.1d-7) then
          temp2=abs(rg2(1,1)/resultn(1,3) -1d0) 
          else
          temp2=abs(rg2(1,1) -resultn(1,3)  ) 
          endif

c Gauge 3
         if (abs(rg3(1,1)).gt.1d-7) then
          temp3=abs(rg3(1,1)/resultn(1,4) -1d0 )
          else
             temp3=abs(rg3(1,1)  -resultn(1,4) )
         endif

c Gauge 4
          if (abs(rg4(1,1)).gt.1d-7) then
          temp4=abs(rg4(1,1)/resultn(1,5) -1d0)
          else
          temp4=abs(rg4(1,1) - resultn(1,5)  )
          endif

c Gauge 5
          if (abs(rg5(1,1)).gt.1d-7) then
          temp5=abs(rg5(1,1)/resultn(1,6) -1d0)
          else
          temp5=abs(rg5(1,1) - resultn(1,6)  )
          endif

      
c FC %          print*, "rg1O",rg1
c FC %          print*, "rg2O",rg2
c FC %          print*, "rg3O",rg3
c FC %          print*, "rg4O",rg4
c FC %          print*, "temp1",temp1
c FC %          print*, "temp2",temp2
c FC %          print*, "temp3",temp3
c FC %          print*, "temp4",temp4
c FC %c	read(*,*)

          resultgauge=max(temp1,temp2,temp3,temp4,temp5)
        
        endif
         
        end
