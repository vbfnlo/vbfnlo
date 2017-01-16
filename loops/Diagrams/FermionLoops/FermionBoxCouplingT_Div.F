c Author: Francisco Campanario
c  Computes the Box fermion loop with generic couplings for the massless case and 
c the massive with only one mass. Only valid for Z,gamma, gluon, couplings.
c  za1,za2,za3,za4 are the P_Right couplings
c  zb1,zb2,zb3,zb4 are the P_left couplings


ccccccccccccccccccccc    Even PART ccccccccccccccccccccccccccccccccccccccccccccccccccc

       subroutine FermionBoxCoupling_Div(za1,za2,za3,za4,zb1,zb2,zb3,zb4,
     & m,p1,p2,p3,p4,mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge,Div)
       implicit none
       Real*8 za1,za2,za3,za4,zb1,zb2,zb3,zb4
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       Real*8   p12(0:3),p14(0:3),p23(0:3),p34(0:3)
       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3)
       Complex*16 mup1test(0:3),mup2test(0:3),mup3test(0:3),mup4test(0:3)
       Complex*16 result,resultn(5),resultZERO(5)
       Real*8 musq,m,resultgauge,temp1,temp2,temp3,temp4
       complex*16 rg1,rg2,rg3,rg4
       integer gauge,i,comp,Div

        resultgauge=1d99

       call  FermionBoxEvenCouplings_Div(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,
     & p2,p3,p4,mup1,mup2,mup3,mup4,musq,comp,result,resultn,Div)

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

          call fermionBoxEvenCouplings_Div(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1test,mup2,mup3,mup4,musq,-1,rg1,resultZERO,Div)
          call fermionBoxEvenCouplings_Div(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1,mup2test,mup3,mup4,musq,-1,rg2,resultZERO,Div)
          call fermionBoxEvenCouplings_Div(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1,mup2,mup3test,mup4,musq,-1,rg3,resultZERO,Div)
          call fermionBoxEvenCouplings_Div(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1,mup2,mup3,mup4test,musq,-1,rg4,resultZERO,Div)
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
      
c	read(*,*)

          resultgauge=max(temp1,temp2,temp3,temp4)
        
        endif
         
        end


ccccccccccccccccccccc    Odd PART ccccccccccccccccccccccccccccccccccccccccccccccccccc

       subroutine FermionBoxCouplingO_Div(za1,za2,za3,za4,zb1,zb2,zb3,zb4,
     & m,p1,p2,p3,p4,mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge,Div)
       implicit none
       Real*8 za1,za2,za3,za4,zb1,zb2,zb3,zb4
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       Real*8   p12(0:3),p14(0:3),p23(0:3),p34(0:3)
       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3)
       Complex*16 mup1test(0:3),mup2test(0:3),mup3test(0:3),mup4test(0:3)
       Complex*16 result,resultn(5),resultZERO(5)
       Real*8 musq,m,resultgauge,temp1,temp2,temp3,temp4
       complex*16 rg1,rg2,rg3,rg4
       integer gauge,i,comp,Div

        resultgauge=1d99

       call  FermionBoxOddCouplings_Div(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,
     & p2,p3,p4,mup1,mup2,mup3,mup4,musq,comp,result,resultn,Div)

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

          call fermionBoxOddCouplings_Div(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1test,mup2,mup3,mup4,musq,-1,rg1,resultZERO,Div)
          call fermionBoxOddCouplings_Div(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1,mup2test,mup3,mup4,musq,-1,rg2,resultZERO,Div)
          call fermionBoxOddCouplings_Div(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1,mup2,mup3test,mup4,musq,-1,rg3,resultZERO,Div)
          call fermionBoxOddCouplings_Div(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1,mup2,mup3,mup4test,musq,-1,rg4,resultZERO,Div)
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
      
c	read(*,*)

          resultgauge=max(temp1,temp2,temp3,temp4)
        
        endif
         
        end
