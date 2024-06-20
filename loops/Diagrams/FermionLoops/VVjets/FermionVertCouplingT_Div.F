c Author: Francisco Campanario
c  Computes the Box fermion loop with generic couplings for the massless case and 
c the massive with only one mass. Only valid for Z,gamma, gluon, couplings.
c  za1,za2,za3,za4 are the P_Right couplings
c  zb1,zb2,zb3,zb4 are the P_left couplings


ccccccccccccccccccccc    Even PART ccccccccccccccccccccccccccccccccccccccccccccccccccc


       subroutine FermionVertCoupling_Div(za1,za2,za3,zb1,zb2,zb3,
     & m,p1,p2,p3,mup1,mup2,mup3,musq,gauge,comp,result,resultgauge,Div)
       implicit none
       Real*8 za1,za2,za3,zb1,zb2,zb3
       Real*8   p1(0:3),p2(0:3),p3(0:3)
       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3)
       Complex*16 mup1test(0:3),mup2test(0:3),mup3test(0:3)
       Complex*16 result,resultn(4),resultZERO(4)
       Real*8 musq,m,resultgauge,temp1,temp2,temp3
       complex*16 rg1,rg2,rg3
       integer gauge,i,comp,Div

        resultgauge=1d99

       call  FermionVertEvenCouplings_Div(za1,za2,za3,zb1,zb2,zb3,m,p1,
     & p2,p3,mup1,mup2,mup3,musq,comp,result,resultn,Div)

       if (gauge.gt.0) then
       
         do  i=0,3
            mup1test(i)=p1(i)
            mup2test(i)=p2(i)
            mup3test(i)=p3(i)
         enddo
       
c Call for Gauge

          call fermionVertEvenCouplings_Div(za1,za2,za3,zb1,zb2,zb3,m,
     &     p1,p2,p3,mup1test,mup2,mup3,musq,-1,rg1,resultZERO,Div)
          call fermionVertEvenCouplings_Div(za1,za2,za3,zb1,zb2,zb3,m,
     &     p1,p2,p3,mup1,mup2test,mup3,musq,-1,rg2,resultZERO,Div)
          call fermionVertEvenCouplings_Div(za1,za2,za3,zb1,zb2,zb3,m,
     &     p1,p2,p3,mup1,mup2,mup3test,musq,-1,rg3,resultZERO,Div)

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

     
c FC %          print*, "rg1",rg1
c FC %          print*, "rg2",rg2
c FC %          print*, "rg3",rg3
c FC %          print*, "rg4",rg4
c FC %          print*, "temp1",temp1
c FC %          print*, "temp2",temp2
c FC %          print*, "temp3",temp3
c FC %          print*, "temp4",temp4
c FC %c	read(*,*)

          resultgauge=max(temp1,temp2,temp3)
        
        endif
         
        end


cccccccccccccccccccc    ODD PART ccccccccccccccccccccccccccccccccccccccccccccccccccc


       subroutine FermionVertCouplingO_Div(za1,za2,za3,zb1,zb2,zb3,
     & m,p1,p2,p3,mup1,mup2,mup3,musq,gauge,comp,result,resultgauge,Div)
       implicit none
       Real*8 za1,za2,za3,zb1,zb2,zb3
       Real*8   p1(0:3),p2(0:3),p3(0:3)
       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3)
       Complex*16 mup1test(0:3),mup2test(0:3),mup3test(0:3)
       Complex*16 result,resultn(4),resultZERO(4)
       Real*8 musq,m,resultgauge,temp1,temp2,temp3
       complex*16 rg1,rg2,rg3
       integer gauge,i,comp,Div

        resultgauge=1d99

       call  FermionVertOddCouplings_Div(za1,za2,za3,zb1,zb2,zb3,m,p1,
     & p2,p3,mup1,mup2,mup3,musq,comp,result,resultn,Div)

       if (gauge.gt.0) then
       
         do  i=0,3
            mup1test(i)=p1(i)
            mup2test(i)=p2(i)
            mup3test(i)=p3(i)
    
         enddo
       
c Call for Gauge

          call fermionVertOddCouplings_Div(za1,za2,za3,zb1,zb2,zb3,m,
     &     p1,p2,p3,mup1test,mup2,mup3,musq,-1,rg1,resultZERO,Div)
          call fermionVertOddCouplings_Div(za1,za2,za3,zb1,zb2,zb3,m,
     &     p1,p2,p3,mup1,mup2test,mup3,musq,-1,rg2,resultZERO,Div)
          call fermionVertOddCouplings_Div(za1,za2,za3,zb1,zb2,zb3,m,
     &     p1,p2,p3,mup1,mup2,mup3test,musq,-1,rg3,resultZERO,Div)
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


c FC %          print*, "rg1O",rg1
c FC %          print*, "rg2O",rg2
c FC %          print*, "rg3O",rg3
c FC %          print*, "rg4O",rg4
c FC %          print*, "temp1",temp1
c FC %          print*, "temp2",temp2
c FC %          print*, "temp3",temp3
c FC %          print*, "temp4",temp4
c FC %c	read(*,*)

          resultgauge=max(temp1,temp2,temp3)
        
        endif
         
        end



