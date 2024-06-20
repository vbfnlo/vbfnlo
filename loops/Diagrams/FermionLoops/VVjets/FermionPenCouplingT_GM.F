c Author: Francisco Campanario
c  Computes the Box fermion loop with generic couplings and masses
c  cpl(1)-cpl(4)  P_Right couplings
c  cpl(5)-cpl(8)  P_left couplingss

        subroutine FermionPenCouplingTS_GM345(za3,za4,za5,zb3,zb4,zb5,
     &  m1,m2,m3,p1,p2,p3,p4,p5,
     &  mup1,mup2,mup3,mup4,mup5,musq,gauge,comp,result,resultgauge,Div)
       implicit none
       Real*8 za3,za4,za5,zb3,zb4,zb5
       Real*8  cpl(10),m1,m2,m3
       Real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
c       Real*8 p12(0:3),p14(0:3),p23(0:3),p34(0:3)
       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3),mup5(0:3)
       Complex*16 result,resultEO(2)
       Real*8 musq,resultgauge
       integer gauge,comp,Div

       cpl(3)=za3
       cpl(4)=za4
       cpl(5)=za5

       cpl(8)=zb3
       cpl(9)=zb4
       cpl(10)=zb5

       call  FermionPenCouplingT_GM345(cpl,m1,m2,m3,p1,p2,p3,p4,p5,
     &  mup1,mup2,mup3,mup4,mup5,
     &  musq,gauge,comp,resultEO,resultgauge,Div)

       result = resultEO(1)+resultEO(2)

       end


        subroutine FermionPenCouplingT1_GM345(za3,za4,za5,zb3,zb4,zb5,
     &  m1,m2,m3,p1,p2,p3,p4,p5,
     &  mup1,mup2,mup3,mup4,mup5,musq,gauge,comp,resultE,resultO,resultgauge,Div)
       implicit none
       Real*8 za3,za4,za5,zb3,zb4,zb5
       Real*8  cpl(10),m1,m2,m3
       Real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
c       Real*8 p12(0:3),p14(0:3),p23(0:3),p34(0:3)
       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3),mup5(0:3)
       Complex*16 result(2),resultE,resultO
       Real*8 musq,resultgauge
       integer gauge,comp,Div

       cpl(3)=za3
       cpl(4)=za4
       cpl(5)=za5


       cpl(8)=zb3
       cpl(9)=zb4
       cpl(10)=zb5

        call FermionPenCouplingT_GM345(cpl,m1,m2,m3,p1,p2,p3,p4,p5,
     &  mup1,mup2,mup3,mup4,mup5, 
     &  musq,gauge,comp,result,resultgauge,Div)

        resultE=result(1)
        resultO=result(2)

       end


        subroutine FermionPenCouplingT_GM345(cpl,m1,m2,m3,p1,p2,p3,p4,p5,
     &  mup1,mup2,mup3,mup4,mup5,
     &  musq,gauge,comp,result,resultgauge,Div)
       implicit none
       Real*8  cpl(10),m1,m2,m3
       Real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
c       Real*8 p12(0:3),p14(0:3),p23(0:3),p34(0:3)
       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3),mup5(0:3)
       Complex*16 result(2),resultE,resultO
       Real*8 musq,resultgauge,resultgaugeE,resultgaugeO
       integer gauge,comp,Div


       call  FermionPenCouplingTE_G345(cpl,m1,m2,m3,p1,p2,p3,p4,p5,
     &  mup1,mup2,mup3,mup4,mup5,
     &  musq,gauge,comp,resultE,resultgaugeE,Div)


       call  FermionPenCouplingTO_G345(cpl,m1,m2,m3,p1,p2,p3,p4,p5,
     &  mup1,mup2,mup3,mup4,mup5, 
     &  musq,gauge,comp,resultO,resultgaugeO,Div)

c       print*, "resultE", resultE
c       print*, "resultO", resultO

C        result=resultE!resultO

      result(1)=resultE
      result(2)=resultO

       resultgauge=Max(resultgaugeO,resultgaugeE)

       end


cccccccccccccccccccc    Even PART ccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine FermionPenCouplingTE_G345(cpl,m1,m2,m3,p1,p2,p3,p4,p5,
     &  mup1,mup2,mup3,mup4,mup5,musq,gauge,comp,result,resultgauge,Div)
       implicit none
       Real*8 cpl(10),m1,m2,m3
       Real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3),mup5(0:3)
       Complex*16 result
       Real*8 musq,resultgauge
       integer gauge,comp,Div

       If (Div.eq.0) then
       call FermionPenCoupling_G345(cpl,m1,m2,m3,p1,p2,p3,p4,p5,
     & mup1,mup2,mup3,mup4,mup5,musq,gauge,comp,result,resultgauge)
       else
       call FermionPenCoupling_G345_Div(cpl,m1,m2,m3,p1,p2,
     & p3,p4,p5,mup1,mup2,mup3,mup4,mup5,
     & musq,gauge,comp,result,resultgauge,Div)
       endif
       end


       subroutine FermionPenCoupling_G345(cpl,m1,m2,m3,
     & p1,p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5,
     & musq,gauge,comp,result,resultgauge)
       implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/loops/Utilities/GaugeTest.inc"
       Real*8 cpl(10),m1,m2,m3
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3),mup5(0:3)
       Complex*16 mup1test(0:3),mup2test(0:3),mup3test(0:3),mup4test(0:3)
       Complex*16 mup5test(0:3)
       Complex*16 result,resultn1(1,6)
       Real*8 musq,m,resultgauge,temp1,temp2,temp3,temp4
       real*8 resultgauge_cp
       complex*16 result_cp
       complex*16 rg1(1,6),rg2(1,6),rg3(1,6),rg4(1,6)
       integer gauge,i,comp
       common/GramSwitch/GramCS,GramDS
       logical GramCS,GramDS
c Checking Gram massive
       real*8 p1sq,p2sq,p3sq,p1p2,p1p3,p2p3
       real*8 dotrr
       external dotrr
       logical printmy
       common/mprint/printmy
       logical debug
       parameter (debug=.false.)
       real*8 Sensitivity, FinalSensitivity
c For QUAD
       Real*8 resultgauge_q
       complex*16 result_Q
       complex*16 resultn1_Q(1,6)
       complex*16 rg1_Q(1,6),rg2_Q(1,6),rg3_Q(1,6)
       complex*16 rg4_Q(1,6)

        Sensitivity=GaugeAccuracyPenL
        FinalSensitivity=GaugeAccuracyPenL
     

        printmy=.false.
        GramCS=.false.
        GramDS=.false.
        resultgauge=1d99
c        dble1=.false.


c        if(dble1) then
       call  PenEvenGSim(cpl,m1,m2,m3,p1,
     & p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5,musq,comp,resultn1)

       result=resultn1(1,1)    

c       print*, "result",result
c        pause
c       else
c       call  PenEvenGSim_QUAD(cpl,m1,m2,m3,p1,
c     & p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5,musq,comp,resultn1_Q)
c       
c      result_Q=resultn1_Q(1,1)  
c
c      if(abs(result/result_Q -1d0).gt.1d1)then
c        print*, "result_Q",result_Q
c        print*, "ratio", result/result_Q-1d0
cc        pause
c        endif
c       endif
c       stop


      if (gauge.gt.0) then
       
         do  i=0,3
            mup1test(i)=p1(i)
            mup2test(i)=p2(i)
            mup3test(i)=p3(i)
            mup4test(i)=p4(i)
            mup5test(i)=p5(i)
         enddo

c Call for Gauge

          call PenEvenGSim(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1test,mup2,mup3,mup4,mup5,
     &     musq,-1,rg1)
          call PenEvenGSim(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2test,mup3,mup4,mup5,
     &     musq,-1,rg2)
          call PenEvenGSim(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2,mup3test,mup4,mup5,
     &     musq,-1,rg3)
          call PenEvenGSim(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2,mup3,mup4test,mup5,
     &     musq,-1,rg4)

c Gauge 1
          if (abs(rg1(1,1)).gt.1d-7) then
          temp1=abs(rg1(1,1)/(resultn1(1,2)) -1d0)
          else
          temp1=abs(rg1(1,1)-resultn1(1,2))   
          endif

c Gauge 2
         if (abs(rg2(1,1)).gt.1d-7) then
          temp2=abs(rg2(1,1)/resultn1(1,3) -1d0) 
          else
          temp2=abs(rg2(1,1) -resultn1(1,3)  ) 
          endif

c Gauge 3
         if (abs(rg3(1,1)).gt.1d-7) then
          temp3=abs(rg3(1,1)/resultn1(1,4) -1d0 )
          else
             temp3=abs(rg3(1,1)  -resultn1(1,4) )
         endif

c Gauge 4
          if (abs(rg4(1,1)).gt.1d-7) then
          temp4=abs(rg4(1,1)/resultn1(1,5) -1d0)
          else
          temp4=abs(rg4(1,1) - resultn1(1,5)  )
          endif

          resultgauge=max(temp1,temp2,temp3,temp4)
         
c Safe results
          result_cp=result
          resultgauge_cp=resultgauge

cCCCCCCCCCCCc TO RESCUE

        if(resultgauge.gt.Sensitivity) then
cmr        if(.false.) then

        GramCS=.true.
        GramDS=.true.
        resultgauge=1d99

       call  PenEvenGSim(cpl,m1,m2,m3,p1,
     & p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5,musq,comp,resultn1)

       result=resultn1(1,1)    

c Call for Gauge

          call PenEvenGSim(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1test,mup2,mup3,mup4,mup5,
     &     musq,-1,rg1)
          call PenEvenGSim(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2test,mup3,mup4,mup5,
     &     musq,-1,rg2)
          call PenEvenGSim(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2,mup3test,mup4,mup5,
     &     musq,-1,rg3)
          call PenEvenGSim(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2,mup3,mup4test,mup5,
     &     musq,-1,rg4)

c Gauge 1
          if (abs(rg1(1,1)).gt.1d-7) then
          temp1=abs(rg1(1,1)/(resultn1(1,2)) -1d0)
          else
          temp1=abs(rg1(1,1)-resultn1(1,2))   
          endif

c Gauge 2
         if (abs(rg2(1,1)).gt.1d-7) then
          temp2=abs(rg2(1,1)/resultn1(1,3) -1d0) 
          else
          temp2=abs(rg2(1,1) -resultn1(1,3)  ) 
          endif

c Gauge 3
         if (abs(rg3(1,1)).gt.1d-7) then
          temp3=abs(rg3(1,1)/resultn1(1,4) -1d0 )
          else
             temp3=abs(rg3(1,1)  -resultn1(1,4) )
         endif

c Gauge 4
          if (abs(rg4(1,1)).gt.1d-7) then
          temp4=abs(rg4(1,1)/resultn1(1,5) -1d0)
          else
          temp4=abs(rg4(1,1) - resultn1(1,5)  )
          endif

          resultgauge=max(temp1,temp2,temp3,temp4)

CCC  CHECK whether if improves precision

       if(resultgauge.lt. FinalSensitivity
     &    .or.resultgauge_cp.lt.FinalSensitivity) then
           if(resultgauge.lt.resultgauge_cp) then 
cFC             print*, "PE:to rescued works",Real(m1),Real(resultgauge),
cFC     &        Real(resultgauge_cp),Abs(CMPLX(result_cp)/CMPLX(result)-1)
cFC     &        ,CMPLX(result)
           else
cFC             print*
cFC             Print*, "PE:Original Better",Real(m1),Real(resultgauge),
cFC     &       Real(resultgauge_cp),Abs(CMPLX(result_cp)/CMPLX(result)-1),
cFC     &                           CMPLX(result)
             result=result_cp
             resultgauge=resultgauge_cp
           endif   
        else
           if(resultgauge.gt.resultgauge_cp) then 
              result=result_cp
             resultgauge=resultgauge_cp
           endif
cFC           print*
cFC           print*, "PE:to rescued failed",Real(m1),Real(resultgauge),Real(resultgauge_cp),Abs(CMPLX(result_cp)/CMPLX(result)-1),CMPLX(result)
           if(abs(result_cp).eq.0d0) then
cFC             print*
cFC             print*, "PE:failed set to WORKS, amplitude is ZERO", result_cp    
             result=result_cp
            resultgauge=Sensitivity
           endif
c Debugging
             if(debug) then
             printMY=.true.

            print*,"PEN here",m1
             p1sq=dotrr(p1,p1)
             p2sq=dotrr(p2,p2)
             p3sq=dotrr(p3,p3)
c             p4sq=dotrr(p1,p1)
             p1p2=dotrr(p1,p2)
             p1p3=dotrr(p1,p3)
             p2p3=dotrr(p2,p3)

             print*,"inv",p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,p1sq+p2sq+p3sq+ 2d0*(p1p2+p1p3+p2p3)

             print*,"inv2",p1sq,p2sq,p3sq,p1sq+p2sq+p3sq+ 2d0*(p1p2+p1p3+p2p3), 2d0*p1p2+p1sq+p2sq,2d0*p2p3+p2sq+p3sq,m*m


          print*, "rg1",rg1(1,1)
          print*, "rg2",rg2(1,1)
          print*, "rg3",rg3(1,1)
          print*, "rg4",rg4(1,1)
          print*, "result1",resultn1(1,2)
          print*, "result",resultn1(1,3)
          print*, "result3",resultn1(1,4)
          print*, "result4",resultn1(1,5)

          print*, "temp1",temp1
          print*, "temp2",temp2
          print*, "temp3",temp3
          print*, "temp4",temp4

             do i=0,3
            write(*,'(A3,I1,A2,E25.16,E25.16)') "p1(",i,")=",p1(i)
            enddo
             do i=0,3
            write(*,'(A3,I1,A2,E25.16,E25.16)') "p2(",i,")=",p2(i)
            enddo
             do i=0,3
            write(*,'(A3,I1,A2,E25.16,E25.16)') "p3(",i,")=",p3(i)
            enddo
             do i=0,3
            write(*,'(A3,I1,A2,E25.16,E25.16)') "p4(",i,")=",p4(i)
             enddo
             do i=0,3
            write(*,'(A3,I1,A2,E25.16,E25.16)') "p5(",i,")=",p5(i)
             enddo

             PRINT*, "mASS",M1
       call  PenEvenGSim(cpl,m1,m2,m3,p1,
     & p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5,musq,comp,resultn1)
c       pause
c       pause
              endif  ! degug
            endif   ! FinalSensitivity


#ifdef WITH_QUAD
            if(abs(resultgauge).gt.FinalSensitivity) then
cfc            print*, "HERE_QUAD"
cfc
cfc           do i=0,3
cfc            write(*,'(A3,I1,A2,E25.16,E25.16)') "p1(",i,")=",p1(i)
cfc            enddo
cfc             do i=0,3
cfc            write(*,'(A3,I1,A2,E25.16,E25.16)') "p2(",i,")=",p2(i)
cfc            enddo
cfc             do i=0,3
cfc            write(*,'(A3,I1,A2,E25.16,E25.16)') "p3(",i,")=",p3(i)
cfc            enddo
cfc             do i=0,3
cfc            write(*,'(A3,I1,A2,E25.16,E25.16)') "p4(",i,")=",p4(i)
cfc             enddo
cfc             do i=0,3
cfc            write(*,'(A3,I1,A2,E25.16,E25.16)') "p5(",i,")=",p5(i)
cfc             enddo


        resultgauge_Q=1d99

       call  PenEvenGSim_QUAD(cpl,m1,m2,m3,p1,
     & p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5,musq,comp,resultn1_Q)
       

       result_Q=resultn1_Q(1,1)


c Call for Gauge

          call PenEvenGSim_QUAD(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1test,mup2,mup3,mup4,mup5,
     &     musq,-1,rg1_Q)
          call PenEvenGSim_QUAD(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2test,mup3,mup4,mup5,
     &     musq,-1,rg2_Q)
          call PenEvenGSim_QUAD(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2,mup3test,mup4,mup5,
     &     musq,-1,rg3_Q)
          call PenEvenGSim_QUAD(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2,mup3,mup4test,mup5,
     &     musq,-1,rg4_Q)

c Gauge 1
          if (abs(rg1_Q(1,1)).gt.1d-7) then
          temp1=abs(rg1_Q(1,1)/(resultn1_Q(1,2)) -1d0)
          else
          temp1=abs(rg1_Q(1,1)-resultn1_Q(1,2))   
          endif

c Gauge 2
         if (abs(rg2_Q(1,1)).gt.1d-7) then
          temp2=abs(rg2_Q(1,1)/resultn1_Q(1,3) -1d0) 
          else
          temp2=abs(rg2_Q(1,1) -resultn1_Q(1,3)  ) 
          endif

c Gauge 3
         if (abs(rg3_Q(1,1)).gt.1d-7) then
          temp3=abs(rg3_Q(1,1)/resultn1_Q(1,4) -1d0 )
          else
             temp3=abs(rg3_Q(1,1)  -resultn1_Q(1,4) )
         endif

c Gauge 4
          if (abs(rg4_Q(1,1)).gt.1d-7) then
          temp4=abs(rg4_Q(1,1)/resultn1_Q(1,5) -1d0)
          else
          temp4=abs(rg4_Q(1,1) - resultn1_Q(1,5)  )
          endif

          resultgauge_Q=max(temp1,temp2,temp3,temp4)

c       print*, "resultQ",result_Q
c       print*, "result_cp",result_cp
c       print*, "result",result
c
c       print*, "resultgaugeQ",resultgauge_Q
c       print*, "resultgauge_cp",resultgauge_cp
c       print*, "resultgauge",resultgauge
c
c       stop

        if(resultgauge_Q.lt. FinalSensitivity) then
c             print*, "PE:QUAD to rescued works",Real(m1),Real(resultgauge),
c     &        Real(resultgauge_Q),Abs(CMPLX(result)/CMPLX(result_Q)-1)
c     &        ,CMPLX(result_Q)
             result=result_Q
             resultgauge=resultgauge_Q

        else
           if(abs(resultgauge_Q).lt.abs(resultgauge)) then
             result=result_Q
             resultgauge=resultgauge_Q
           endif
c           print*
c           print*, "PE:QUAD to rescued failed",Real(m1),Real(resultgauge),
c     &  Real(resultgauge_cp),Abs(CMPLX(result_cp)/CMPLX(result)-1),CMPLX(result)
           endif
           endif
#endif


       endif  ! Sensitivity
       endif ! Gauge


cFC       print*, "comp",comp
cFC       print*,"PenEvenSim",result
cFC       print*,"Pen Sim Gague 1",resultn1(1,2)
cFC       print*,"Pen Sim Gague 2",resultn1(1,3)
cFC       print*,"Pen Sim Gague 3",resultn1(1,4)
cFC       print*,"Pen Sim Gague 4",resultn1(1,5)
cFC       print*,"Pen Sim Gague 5",resultn1(1,6)
cFC
cFC       call  PenEvenG(cpl,m1,m2,m3,p1,
cFC     & p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5,musq,comp,resultn)
cFC
cFC       result=resultn(1,1)
cFC
cFC       print*,"PenOddTotal",resultn(1,1)/resultn1(1,1)
cFC       print*,"Pen  Gague 1",resultn(1,2)/resultn1(1,2)
cFC       print*,"Pen  Gague 2",resultn(1,3)/resultn1(1,3)
cFC       print*,"Pen  Gague 3",resultn(1,4)/resultn1(1,4)
cFC       print*,"Pen  Gague 4",resultn(1,5)/resultn1(1,5)
cFC       print*,"Pen  Gague 5",resultn(1,6)/resultn1(1,6)
cFCc       pause

cFC          print*, "rg1O",rg1(1,1)
cFC          print*, "rg2O",rg2(1,1)
cFC          print*, "rg3O",rg3(1,1)
cFC          print*, "rg4O",rg4(1,1)
cFC          print*, "rg5O",rg5(1,1)
cFC          print*, "temp1",temp1
cFC          print*, "temp2",temp2
cFC          print*, "temp3",temp3
cFC          print*, "temp4",temp4
cFC          print*, "temp5",temp5
cFC
cFC
cFC       
cFCc Call for Gauge
cFC
cFC          call PenEvenG(cpl,m1,m2,m3,
cFC     &     p1,p2,p3,p4,p5,mup1test,mup2,mup3,mup4,mup5,
cFC     &     musq,-1,rg1)
cFC          call PenEvenG(cpl,m1,m2,m3,
cFC     &     p1,p2,p3,p4,p5,mup1,mup2test,mup3,mup4,mup5,
cFC     &     musq,-1,rg2)
cFC          call PenEvenG(cpl,m1,m2,m3,
cFC     &     p1,p2,p3,p4,p5,mup1,mup2,mup3test,mup4,mup5,
cFC     &     musq,-1,rg3)
cFC          call PenEvenG(cpl,m1,m2,m3,
cFC     &     p1,p2,p3,p4,p5,mup1,mup2,mup3,mup4test,mup5,
cFC     &     musq,-1,rg4)
cFC          call PenEvenG(cpl,m1,m2,m3,
cFC     &     p1,p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5test,
cFC     &     musq,-1,rg5)
cFC
cFC
cFCc Gauge 1
cFC          if (abs(rg1(1,1)).gt.1d-7) then
cFC          temp1=abs(rg1(1,1)/(resultn(1,2)) -1d0)
cFC          else
cFC          temp1=abs(rg1(1,1)-resultn(1,2))   
cFC          endif
cFC
cFCc Gauge 2
cFC         if (abs(rg2(1,1)).gt.1d-7) then
cFC          temp2=abs(rg2(1,1)/resultn(1,3) -1d0) 
cFC          else
cFC          temp2=abs(rg2(1,1) -resultn(1,3)  ) 
cFC          endif
cFC
cFCc Gauge 3
cFC         if (abs(rg3(1,1)).gt.1d-7) then
cFC          temp3=abs(rg3(1,1)/resultn(1,4) -1d0 )
cFC          else
cFC             temp3=abs(rg3(1,1)  -resultn(1,4) )
cFC         endif
cFC
cFCc Gauge 4
cFC          if (abs(rg4(1,1)).gt.1d-7) then
cFC          temp4=abs(rg4(1,1)/resultn(1,5) -1d0)
cFC          else
cFC          temp4=abs(rg4(1,1) - resultn(1,5)  )
cFC          endif
cFC
cFCc Gauge 5
cFC          if (abs(rg5(1,1)).gt.1d-7) then
cFC          temp5=abs(rg5(1,1)/resultn(1,6) -1d0)
cFC          else
cFC          temp5=abs(rg5(1,1) - resultn(1,6))
cFC          endif
cFC
cFC          print*,
cFC          print*, "rg1",rg1(1,1)
cFC          print*, "rg2",rg2(1,1)
cFC          print*, "rg3",rg3(1,1)
cFC          print*, "rg4",rg4(1,1)
cFC          print*, "temp1",temp1
cFC          print*, "temp2",temp2
cFC          print*, "temp3",temp3
cFC          print*, "temp4",temp4
cFC          print*, "temp5",temp5
cFC          pause
cFC          print*, "result1",resultn(1)
cFC          print*, "result",resultn(2)
cFC          print*, "result3",resultn(3)
cFC          print*, "result4",resultn(4)
cFC          print*, "result5",resultn(5)
c FC %c	read(*,*)
        end

cFC
cFC
cFCcccccccccccccccccccc    ODD PART ccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine FermionPenCouplingTO_G345(cpl,m1,m2,m3,p1,p2,p3,p4,p5,
     &  mup1,mup2,mup3,mup4,mup5,musq,gauge,comp,result,resultgauge,Div)
       implicit none
       Real*8 cpl(10),m1,m2,m3
       Real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3),mup5(0:3)


       Complex*16 result
       Real*8 musq,resultgauge

       integer gauge,comp,Div

       If (Div.eq.0) then
       call FermionPenCouplingO_G345(cpl,m1,m2,m3,p1,p2,p3,p4,p5,
     & mup1,mup2,mup3,mup4,mup5,musq,gauge,comp,result,resultgauge)
       else
       call FermionPenCouplingO_G345_Div(cpl,m1,m2,m3,p1,p2,
     & p3,p4,p5,mup1,mup2,mup3,mup4,mup5,musq,gauge,comp,
     & result,resultgauge,Div)
       endif
       end


       subroutine FermionPenCouplingO_G345(cpl,m1,m2,m3,
     & p1,p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5,
     & musq,gauge,comp,result,resultgauge)
       implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/loops/Utilities/GaugeTest.inc"
       Real*8 cpl(10),m1,m2,m3
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3),mup5(0:3)
       Complex*16 mup1test(0:3),mup2test(0:3),mup3test(0:3),mup4test(0:3)
       Complex*16 mup5test(0:3)
       Complex*16 result,resultn1(1,6)
       Real*8 musq,m,resultgauge,temp1,temp2,temp3,temp4,temp5
       real*8 resultgauge_cp
       complex*16 result_cp
       complex*16 rg1(1,6),rg2(1,6),rg3(1,6),rg4(1,6),rg5(1,6)
       integer gauge,i,comp
       complex*16 Fa(4390), Fa1(4390)
       COMMON/FaPenOddGSim/Fa
       COMMON/FaPenOddG/Fa1
       common/GramSwitch/GramCS,GramDS
       logical GramCS,GramDS
c Checking Gram massive
       real*8 p1sq,p2sq,p3sq,p1p2,p1p3,p2p3
       real*8 dotrr
       external dotrr
       logical printmy
       common/mprint/printmy
       logical debug
       parameter (debug=.false.)
       real*8 Sensitivity, FinalSensitivity
c For QUAD
       Real*8 resultgauge_q
       complex*16 result_Q
       complex*16 resultn1_Q(1,6)
       complex*16 rg2_Q(1,6),rg3_Q(1,6)
       complex*16 rg5_Q(1,6)


        Sensitivity=GaugeAccuracyPenL
        FinalSensitivity=GaugeAccuracyPenL

        printmy=.false.
        GramCS=.false.
        GramDS=.false.
        resultgauge=1d99

       call  PenOddGSim(cpl,m1,m2,m3,p1,
     & p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5,musq,comp,resultn1)

       result=resultn1(1,1)    

       if (gauge.gt.0) then
      
         do  i=0,3
            mup1test(i)=p1(i)
            mup2test(i)=p2(i)
            mup3test(i)=p3(i)
            mup4test(i)=p4(i)
            mup5test(i)=p5(i)
         enddo
       
c Call for Gauge

          call PenOddGSim(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2test,mup3,mup4,mup5,
     &     musq,-1,rg2)
          call PenOddGSim(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2,mup3test,mup4,mup5,
     &     musq,-1,rg3)

          call PenOddGSim(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5test,
     &     musq,-1,rg5)

c Gauge 2
         if (abs(rg2(1,1)).gt.1d-7) then
          temp2=abs(rg2(1,1)/resultn1(1,3) -1d0) 
          else
          temp2=abs(rg2(1,1) -resultn1(1,3)  ) 
          endif

c Gauge 3
         if (abs(rg3(1,1)).gt.1d-7) then
          temp3=abs(rg3(1,1)/resultn1(1,4) -1d0 )
          else
             temp3=abs(rg3(1,1)  -resultn1(1,4) )
         endif
c Gauge 5
          if (abs(rg5(1,1)).gt.1d-7) then
          temp5=abs(rg5(1,1)/resultn1(1,6) -1d0)
          else
          temp5=abs(rg5(1,1) - resultn1(1,6)  )
          endif

          resultgauge=max(temp2,temp3,temp5)

c Safe results
          result_cp=result
          resultgauge_cp=resultgauge

cCCCCCCCCCCCc TO RESCUE

cmr        if(.false.) then
        if(resultgauge.gt.Sensitivity) then

        GramCS=.true.
        GramDS=.true.
        resultgauge=1d99


       call  PenOddGSim(cpl,m1,m2,m3,p1,
     & p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5,musq,comp,resultn1)

       result=resultn1(1,1)    
      
c Call for Gauge

          call PenOddGSim(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2test,mup3,mup4,mup5,
     &     musq,-1,rg2)
          call PenOddGSim(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2,mup3test,mup4,mup5,
     &     musq,-1,rg3)

          call PenOddGSim(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5test,
     &     musq,-1,rg5)

c Gauge 2
         if (abs(rg2(1,1)).gt.1d-7) then
          temp2=abs(rg2(1,1)/resultn1(1,3) -1d0) 
          else
          temp2=abs(rg2(1,1) -resultn1(1,3)  ) 
          endif

c Gauge 3
         if (abs(rg3(1,1)).gt.1d-7) then
          temp3=abs(rg3(1,1)/resultn1(1,4) -1d0 )
          else
             temp3=abs(rg3(1,1)  -resultn1(1,4) )
         endif
c Gauge 5
          if (abs(rg5(1,1)).gt.1d-7) then
          temp5=abs(rg5(1,1)/resultn1(1,6) -1d0)
          else
          temp5=abs(rg5(1,1) - resultn1(1,6)  )
          endif

          resultgauge=max(temp2,temp3,temp5)


CCC  CHECK whether if improves precision

       if(resultgauge.lt. FinalSensitivity
     &    .or.resultgauge_cp.lt.FinalSensitivity) then
           if(resultgauge.lt.resultgauge_cp) then 
cFC             print*, "PO:to rescued works",Real(m1),Real(resultgauge),
cFC     &        Real(resultgauge_cp),Abs(CMPLX(result_cp)/CMPLX(result)-1)
cFC     &        ,CMPLX(result)
           else
cFC             print*
cFC             Print*, "PO:Original Better",Real(m1),Real(resultgauge),
cFC     &       Real(resultgauge_cp),Abs(CMPLX(result_cp)/CMPLX(result)-1),
cFC     &                           CMPLX(result)
             result=result_cp
             resultgauge=resultgauge_cp
           endif   
        else
           if(resultgauge.gt.resultgauge_cp) then 
              result=result_cp
             resultgauge=resultgauge_cp
           endif   
cFC           print*
cFC           print*, "PO:to rescued failed",Real(m1),Real(resultgauge)
cFC     &      ,Real(resultgauge_cp),Abs(CMPLX(result_cp)/CMPLX(result)-1)
cFC     &      ,CMPLX(result)
           if(abs(result_cp).eq.0d0) then
cFC             print*
cFC             print*, "PO:failed set to WORKS, amplitude is ZERO", result_cp    
             result=result_cp
            resultgauge=Sensitivity
           endif
c Debugging
             if(debug) then
             printMY=.true.

            print*,"PEN ODD here",m
             p1sq=dotrr(p1,p1)
             p2sq=dotrr(p2,p2)
             p3sq=dotrr(p3,p3)
c             p4sq=dotrr(p1,p1)
             p1p2=dotrr(p1,p2)
             p1p3=dotrr(p1,p3)
             p2p3=dotrr(p2,p3)

             print*,"inv",p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,p1sq+p2sq+p3sq+ 2d0*(p1p2+p1p3+p2p3)

             print*,"inv2",p1sq,p2sq,p3sq,p1sq+p2sq+p3sq+ 2d0*(p1p2+p1p3+p2p3), 2d0*p1p2+p1sq+p2sq,2d0*p2p3+p2sq+p3sq,m*m


          print*, "rg1",rg1(1,1)
          print*, "rg2",rg2(1,1)
          print*, "rg3",rg3(1,1)
          print*, "rg4",rg4(1,1)
          print*, "result1",resultn1(1,2)
          print*, "result",resultn1(1,3)
          print*, "result3",resultn1(1,4)
          print*, "result4",resultn1(1,5)

          print*, "temp1",temp1
          print*, "temp2",temp2
          print*, "temp3",temp3
          print*, "temp4",temp4

             do i=0,3
            write(*,'(A3,I1,A2,E25.16,E25.16)') "p1(",i,")=",p1(i)
            enddo
             do i=0,3
            write(*,'(A3,I1,A2,E25.16,E25.16)') "p2(",i,")=",p2(i)
            enddo
             do i=0,3
            write(*,'(A3,I1,A2,E25.16,E25.16)') "p3(",i,")=",p3(i)
            enddo
             do i=0,3
            write(*,'(A3,I1,A2,E25.16,E25.16)') "p4(",i,")=",p4(i)
             enddo
             do i=0,3
            write(*,'(A3,I1,A2,E25.16,E25.16)') "p5(",i,")=",p5(i)
             enddo

             PRINT*, "mASS",M

       call  PenOddGSim(cpl,m1,m2,m3,p1,
     & p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5,musq,comp,resultn1)
c       pause
              endif  ! degug
            endif   ! FinalSensitivity


#ifdef WITH_QUAD
        if(resultgauge.gt.FinalSensitivity) then

 
        resultgauge_Q=1d99


       call  PenOddGSim_QUAD(cpl,m1,m2,m3,p1,
     & p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5,musq,comp,resultn1_Q)

       result_Q=resultn1_Q(1,1)    
      
c Call for Gauge

          call PenOddGSim_QUAD(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2test,mup3,mup4,mup5,
     &     musq,-1,rg2_Q)
          call PenOddGSim_QUAD(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2,mup3test,mup4,mup5,
     &     musq,-1,rg3_Q)

          call PenOddGSim_QUAD(cpl,m1,m2,m3,
     &     p1,p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5test,
     &     musq,-1,rg5_Q)

c Gauge 2
         if (abs(rg2_Q(1,1)).gt.1d-7) then
          temp2=abs(rg2_Q(1,1)/resultn1_Q(1,3) -1d0) 
          else
          temp2=abs(rg2_Q(1,1) -resultn1_Q(1,3)  ) 
          endif

c Gauge 3
         if (abs(rg3_Q(1,1)).gt.1d-7) then
          temp3=abs(rg3_Q(1,1)/resultn1_Q(1,4) -1d0 )
          else
             temp3=abs(rg3_Q(1,1)  -resultn1_Q(1,4) )
         endif
c Gauge 5
          if (abs(rg5_Q(1,1)).gt.1d-7) then
          temp5=abs(rg5_Q(1,1)/resultn1_Q(1,6) -1d0)
          else
          temp5=abs(rg5_Q(1,1) - resultn1_Q(1,6)  )
          endif

          resultgauge_Q=max(temp2,temp3,temp5)


        if(resultgauge_Q.lt. FinalSensitivity) then
c             print*, "PO:QUAD to rescued works",Real(m1),Real(resultgauge),
c     &        Real(resultgauge_Q),Abs(CMPLX(result)/CMPLX(result_Q)-1)
c     &        ,CMPLX(result_Q)
             result=result_Q
             resultgauge=resultgauge_Q

        else
           if(abs(resultgauge_Q).lt.abs(resultgauge)) then
             result=result_Q
             resultgauge=resultgauge_Q
           endif
c           print*
c           print*, "PO:QUAD to rescued failed",Real(m1),Real(resultgauge),Real(resultgauge_cp),Abs(CMPLX(result_cp)/CMPLX(result)-1),CMPLX(result)
           endif
           endif

#endif

        endif  ! Sensitivity
        endif ! gauge


cFC       print*, "comp",comp
cFC       print*,"PenOddSim",result
cFC       print*,"Pen Sim Gague 1",resultn1(1,2)
cFC       print*,"Pen Sim Gague 2",resultn1(1,3)
cFC       print*,"Pen Sim Gague 3",resultn1(1,4)
cFC       print*,"Pen Sim Gague 4",resultn1(1,5)
cFC       print*,"Pen Sim Gague 5",resultn1(1,6)
cFC
cFC       call  PenOddG(cpl,m1,m2,m3,p1,
cFC     & p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5,musq,comp,resultn)
cFC
cFC       result=resultn(1,1)
cFC
cFC       print*,"PenOddTotal",result
cFC       print*,"Pen  Gague 1",resultn(1,2)/resultn1(1,2)
cFC       print*,"Pen  Gague 2",resultn(1,3)/resultn1(1,3)
cFC       print*,"Pen  Gague 3",resultn(1,4)/resultn1(1,4)
cFC       print*,"Pen  Gague 4",resultn(1,5)/resultn1(1,5)
cFC       print*,"Pen  Gague 5",resultn(1,6)/resultn1(1,6)
cFC          print*, "rg1O",rg1(1,1)
cFC          print*, "rg2O",rg2(1,1)
cFC          print*, "rg3O",rg3(1,1)
cFC          print*, "rg4O",rg4(1,1)
cFC          print*, "rg5O",rg5(1,1)
cFC          print*, "temp1",temp1
cFC          print*, "temp2",temp2
cFC          print*, "temp3",temp3
cFC          print*, "temp4",temp4
cFC          print*, "temp5",temp5
cFCcFC          pause
cFCcFC
cFC
cFC          call PenOddG(cpl,m1,m2,m3,
cFC     &     p1,p2,p3,p4,p5,mup1test,mup2,mup3,mup4,mup5,
cFC     &     musq,-1,rg1)
cFC          call PenOddG(cpl,m1,m2,m3,
cFC     &     p1,p2,p3,p4,p5,mup1,mup2test,mup3,mup4,mup5,
cFC     &     musq,-1,rg2)
cFC          call PenOddG(cpl,m1,m2,m3,
cFC     &     p1,p2,p3,p4,p5,mup1,mup2,mup3test,mup4,mup5,
cFC     &     musq,-1,rg3)
cFC          call PenOddG(cpl,m1,m2,m3,
cFC     &     p1,p2,p3,p4,p5,mup1,mup2,mup3,mup4test,mup5,
cFC     &     musq,-1,rg4)
cFC          call PenOddG(cpl,m1,m2,m3,
cFC     &     p1,p2,p3,p4,p5,mup1,mup2,mup3,mup4,mup5test,
cFC     &     musq,-1,rg5)
cFC
cFC
cFC
cFC
cFCc Gauge 1
cFC          if (abs(rg1(1,1)).gt.1d-7) then
cFC          temp1=abs(rg1(1,1)/(resultn(1,2)) -1d0)
cFC          else
cFC          temp1=abs(rg1(1,1)-resultn(1,2))   
cFC          endif
cFC
cFCc Gauge 2
cFC         if (abs(rg2(1,1)).gt.1d-7) then
cFC          temp2=abs(rg2(1,1)/resultn(1,3) -1d0) 
cFC          else
cFC          temp2=abs(rg2(1,1) -resultn(1,3)  ) 
cFC          endif
cFC
cFCc Gauge 3
cFC         if (abs(rg3(1,1)).gt.1d-7) then
cFC          temp3=abs(rg3(1,1)/resultn(1,4) -1d0 )
cFC          else
cFC             temp3=abs(rg3(1,1)  -resultn(1,4) )
cFC         endif
cFC
cFCc Gauge 4
cFC          if (abs(rg4(1,1)).gt.1d-7) then
cFC          temp4=abs(rg4(1,1)/resultn(1,5) -1d0)
cFC          else
cFC          temp4=abs(rg4(1,1) - resultn(1,5)  )
cFC          endif
cFC
cFCc Gauge 5
cFC          if (abs(rg5(1,1)).gt.1d-7) then
cFC          temp5=abs(rg5(1,1)/resultn(1,6) -1d0)
cFC          else
cFC          temp5=abs(rg5(1,1) - resultn(1,6)  )
cFC          endif
cFC
cFC          print*
cFC          print*, "rg1",rg1(1,1)
cFC          print*, "rg2",rg2(1,1)
cFC          print*, "rg3",rg3(1,1)
cFC          print*, "rg4",rg4(1,1)
cFC          print*, "rg5",rg5(1,1)
cFC          print*, "temp1",temp1
cFC          print*, "temp2",temp2
cFC          print*, "temp3",temp3
cFC          print*, "temp4",temp4
cFC          print*, "temp5",temp5
cFC          pause
c FC %c	read(*,*)
c          resultgauge=max(temp1,temp2,temp3,temp4,temp5)
c          print*, "resultgauge", resultgauge
cFC         if(resultgauge.gt.0.1d0) print*, "ODD",temp2,temp3,temp5

         
        end
