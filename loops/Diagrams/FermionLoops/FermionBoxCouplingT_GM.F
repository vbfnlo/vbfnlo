c Author: Francisco Campanario
c  Computes the Box fermion loop with generic couplings and masses
c  cpl(1)-cpl(4)  P_Right couplings
c  cpl(5)-cpl(8)  P_left couplingss

        subroutine FermionBoxCouplingT1_GM(za1,za2,za3,za4,zb1,zb2,zb3,zb4,
     &  m1,m2,m3,m4,p1,p2,p3,p4,
     &  mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge,Div)
       implicit none
       Real*8 za1,za2,za3,za4,zb1,zb2,zb3,zb4
       Real*8  cpl(8),m1,m2,m3,m4
       Real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3)

       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3)
       Complex*16 result
       Real*8 musq,resultgauge
       integer gauge,comp,Div
       cpl(1)=za1
       cpl(2)=za2
       cpl(3)=za3
       cpl(4)=za4
       cpl(5)=zb1
       cpl(6)=zb2
       cpl(7)=zb3
       cpl(8)=zb4
        call FermionBoxCouplingT_GM(cpl,m1,m2,m3,m4,p1,p2,p3,p4,
     &  mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge,Div)

       end


        subroutine FermionBoxCouplingT_GM(cpl,m1,m2,m3,m4,p1,p2,p3,p4,
     &  mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge,Div)
       implicit none
       Real*8  cpl(8),m1,m2,m3,m4
       Real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3)

       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3)
       Complex*16 result,resultE,resultO
       Real*8 musq,resultgauge,resultgaugeE,resultgaugeO
       integer gauge,comp,Div


       call  FermionBoxCouplingTE_G(cpl,m1,m2,m3,m4,p1,p2,p3,p4,
     &  mup1,mup2,mup3,mup4,musq,gauge,comp,resultE,resultgaugeE,Div)


       call  FermionBoxCouplingTO_G(cpl,m1,m2,m3,m4,p1,p2,p3,p4,
     &  mup1,mup2,mup3,mup4,musq,gauge,comp,resultO,resultgaugeO,Div)

c       print*, "resultE", resultE
c       print*, "resultO", resultO

       result=resultE+resultO

       resultgauge=Max(resultgaugeO,resultgaugeE)

       end


cccccccccccccccccccc    Even PART ccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine FermionBoxCouplingTE_G(cpl,m1,m2,m3,m4,p1,p2,p3,p4,
     &  mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge,Div)
       implicit none
       Real*8 cpl(8),m1,m2,m3,m4
       Real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3)

       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3)

       Complex*16 result
       Real*8 musq,resultgauge

       integer gauge,comp,Div

       If (Div.eq.0) then
       call FermionBoxEvenGCoupling(cpl,m1,m2,m3,m4,p1,p2,p3,p4
     & ,mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge)
       else
       call FermionBoxEvenGCoupling_Div(cpl,m1,m2,m3,m4,p1,p2,
     & p3,p4,mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge,Div)
       endif
       end


       subroutine FermionBoxEvenGCoupling(cpl,m1,m2,m3,m4,
     & p1,p2,p3,p4,mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge)
       implicit none
       Real*8 cpl(8),m1,m2,m3,m4
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       Real*8   p12(0:3),p14(0:3),p23(0:3),p34(0:3)
       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3)
       Complex*16 mup1test(0:3),mup2test(0:3),mup3test(0:3),mup4test(0:3)
       Complex*16 result,resultn(5),resultZERO(5)
       Real*8 musq,resultgauge,temp1,temp2,temp3,temp4
       complex*16 rg1,rg2,rg3,rg4
       integer gauge,i,comp

        resultgauge=1d99

       call  FermionBoxEvenGCouplings(cpl,m1,m2,m3,m4,p1,
     & p2,p3,p4,mup1,mup2,mup3,mup4,musq,comp,result,resultn)

       if (gauge.gt.0) then
       
         do  i=0,3
            mup1test(i)=p1(i)
            mup2test(i)=p2(i)
            mup3test(i)=p3(i)
            mup4test(i)=p4(i)
            p12(i)= p1(i)+p2(i)
            p14(i)= p1(i)+p4(i)
            p23(i)= p2(i)+p3(i)
            p34(i)= p3(i)+p4(i)
         enddo
       
c Call for Gauge

          call fermionBoxEvenGCouplings(cpl,m1,m2,m3,m4,
     &     p1,p2,p3,p4,mup1test,mup2,mup3,mup4,musq,-1,rg1,resultZERO)
          call fermionBoxEvenGCouplings(cpl,m1,m2,m3,m4,
     &     p1,p2,p3,p4,mup1,mup2test,mup3,mup4,musq,-1,rg2,resultZERO)
          call fermionBoxEvenGCouplings(cpl,m1,m2,m3,m4,
     &     p1,p2,p3,p4,mup1,mup2,mup3test,mup4,musq,-1,rg3,resultZERO)
          call fermionBoxEvenGCouplings(cpl,m1,m2,m3,m4,
     &     p1,p2,p3,p4,mup1,mup2,mup3,mup4test,musq,-1,rg4,resultZERO)
c Gauge 1
          if (abs(rg1).gt.1d-7) then
          temp1=abs(rg1/(resultn(2)) -1d0)
          else
          temp1=abs(rg1-resultn(2))   
          endif

c Gauge 2
         if (abs(rg2).gt.1d-7) then
          temp2=abs(rg2/resultn(3) -1d0) 
          else
          temp2=abs(rg2 -resultn(3)  ) 
          endif

c Gauge 3
         if (abs(rg3).gt.1d-7) then
          temp3=abs(rg3/resultn(4) -1d0 )
          else
             temp3=abs(rg3  -resultn(4) )
         endif

c Gauge 4
          if (abs(rg4).gt.1d-7) then
          temp4=abs(rg4/resultn(5) -1d0)
          else
          temp4=abs(rg4 - resultn(5)  )
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

          resultgauge=max(temp1,temp2,temp3,temp4)
        
        endif
         
        end

cFC
cFC
cFCcccccccccccccccccccc    ODD PART ccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine FermionBoxCouplingTO_G(cpl,m1,m2,m3,m4,p1,p2,p3,p4,
     &  mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge,Div)
       implicit none
       Real*8 cpl(8),m1,m2,m3,m4
       Real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3)

       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3)

       Complex*16 result
       Real*8 musq,resultgauge

       integer gauge,comp,Div

       If (Div.eq.0) then
       call FermionBoxCouplingO_G(cpl,m1,m2,m3,m4,p1,p2,p3,p4
     & ,mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge)
       else
       call FermionBoxCouplingO_G_Div(cpl,m1,m2,m3,m4,p1,p2,
     & p3,p4,mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge,Div)
       endif
       end


       subroutine FermionBoxCouplingO_G(cpl,m1,m2,m3,m4,
     & p1,p2,p3,p4,mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge)
       implicit none
       Real*8 cpl(8),m1,m2,m3,m4
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       Real*8   p12(0:3),p14(0:3),p23(0:3),p34(0:3)
       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3)
       Complex*16 mup1test(0:3),mup2test(0:3),mup3test(0:3),mup4test(0:3)
       Complex*16 result,resultn(5),resultZERO(5)
       Real*8 musq,resultgauge,temp1,temp2,temp3,temp4
       complex*16 rg1,rg2,rg3,rg4
       integer gauge,i,comp

        resultgauge=1d99

       call  FermionBoxOddCouplings_G(cpl,m1,m2,m3,m4,p1,
     & p2,p3,p4,mup1,mup2,mup3,mup4,musq,comp,result,resultn)

       if (gauge.gt.0) then
       
         do  i=0,3
            mup1test(i)=p1(i)
            mup2test(i)=p2(i)
            mup3test(i)=p3(i)
            mup4test(i)=p4(i)
            p12(i)= p1(i)+p2(i)
            p14(i)= p1(i)+p4(i)
            p23(i)= p2(i)+p3(i)
            p34(i)= p3(i)+p4(i)
         enddo
       
c Call for Gauge

          call fermionBoxOddCouplings_G(cpl,m1,m2,m3,m4,
     &     p1,p2,p3,p4,mup1test,mup2,mup3,mup4,musq,-1,rg1,resultZERO)
          call fermionBoxOddCouplings_G(cpl,m1,m2,m3,m4,
     &     p1,p2,p3,p4,mup1,mup2test,mup3,mup4,musq,-1,rg2,resultZERO)
          call fermionBoxOddCouplings_G(cpl,m1,m2,m3,m4,
     &     p1,p2,p3,p4,mup1,mup2,mup3test,mup4,musq,-1,rg3,resultZERO)
          call fermionBoxOddCouplings_G(cpl,m1,m2,m3,m4,
     &     p1,p2,p3,p4,mup1,mup2,mup3,mup4test,musq,-1,rg4,resultZERO)
c Gauge 1
          if (abs(rg1).gt.1d-7) then
          temp1=abs(rg1/(resultn(2)) -1d0)
          else
          temp1=abs(rg1-resultn(2))   
          endif

c Gauge 2
         if (abs(rg2).gt.1d-7) then
          temp2=abs(rg2/resultn(3) -1d0) 
          else
          temp2=abs(rg2 -resultn(3)  ) 
          endif

c Gauge 3
         if (abs(rg3).gt.1d-7) then
          temp3=abs(rg3/resultn(4) -1d0 )
          else
             temp3=abs(rg3  -resultn(4) )
         endif

c Gauge 4
          if (abs(rg4).gt.1d-7) then
          temp4=abs(rg4/resultn(5) -1d0)
          else
          temp4=abs(rg4 - resultn(5)  )
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

          resultgauge=max(temp1,temp2,temp3,temp4)
        
        endif
         
        end
