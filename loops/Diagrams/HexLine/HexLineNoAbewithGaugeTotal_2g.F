
       subroutine NoAbe2gaHexLineTotal(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musqIn,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,Div,gauge)
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup5(0:3),mup6(0:3)

       Complex*16 barpsi_p4(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result_graph(6),result_color(2)
       real*8 resultgauge(2)
       Integer Div,Gauge
        musq=musqIn
       if(Div.eq.0) then

        call NoAbe2gaHexLineT(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,gauge)

        else
        call NoAbe2gaHexLineDivT(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,Div,gauge)


        endif


      end


        subroutine  NoAbe2gaHexLineT(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musqIn,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,gauge)

       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup5(0:3),mup6(0:3)
       Complex*16   mup2test(0:3),mup3test(0:3)
       Complex*16 barpsi_p4(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon,i 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result_graph(6),result_color(2),result_color_g(2)
     &            ,resulth(6)
       real*8 resultgauge(2)
       Integer Gauge
        musq=musqIn


        call NoAbe2gaHexLine(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color)
        
        if(gauge.eq.1)then

        do i=0,3
          mup2test(i)=p2(i)
          mup3test(i)=p2(i)+p3(i)
        enddo


        call NoAbe2gaHexLine(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -  ,result_color_g)

        resultgauge(1)=(abs((resulth(1)+resulth(2)+resulth(3)+resulth(4)+resulth(5))/(-resulth(6)) -1d0))
        
        elseif(gauge.eq.2)then

        do i=0,3
          mup2test(i)=p2(i)+p3(i)
          mup3test(i)=p3(i)
        enddo

        call NoAbe2gaHexLine(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -   ,result_color_g)
        
        resultgauge(2)=(abs((resulth(1)+resulth(2)+resulth(3)+resulth(4)+resulth(5))/(-resulth(6)) -1d0))


        elseif(gauge.eq.3)then

        do i=0,3
          mup2test(i)=p2(i)
          mup3test(i)=p2(i)+p3(i)
        enddo


        call NoAbe2gaHexLine(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -  ,result_color_g)

        resultgauge(1)=(abs((resulth(1)+resulth(2)+resulth(3)+resulth(4)+resulth(5))/(-resulth(6)) -1d0))
        
         do i=0,3
          mup2test(i)=p2(i)+p3(i)
          mup3test(i)=p3(i)
        enddo

        call NoAbe2gaHexLine(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -   ,result_color_g)
        
        resultgauge(2)=(abs((resulth(1)+resulth(2)+resulth(3)+resulth(4)+resulth(5))/(-resulth(6)) -1d0))


        endif
        
        
        end


        subroutine  NoAbe2gaHexLineDivT(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musqIn,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,Div,gauge)

       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup5(0:3),mup6(0:3)
       Complex*16   mup2test(0:3),mup3test(0:3)
       Complex*16 barpsi_p4(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon,i 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result_graph(6),result_color(2),result_color_g(2)
     &            ,resulth(6)
       real*8 resultgauge(2)
       Integer Div,Gauge
        musq=musqIn

        call NoAbe2gaHexLineDiv(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color,Div)
        
        if(gauge.eq.1)then

        do i=0,3
          mup2test(i)=p2(i)
          mup3test(i)=p2(i)+p3(i)
        enddo


        call NoAbe2gaHexLineDiv(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -  ,result_color_g,Div)

        resultgauge(1)=(abs((resulth(1)+resulth(2)+resulth(3)+resulth(4)+resulth(5))/(-resulth(6)) -1d0))
        
        elseif(gauge.eq.2)then

        do i=0,3
          mup2test(i)=p2(i)+p3(i)
          mup3test(i)=p3(i)
        enddo

        call NoAbe2gaHexLineDiv(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -   ,result_color_g,Div)
        
        resultgauge(2)=(abs((resulth(1)+resulth(2)+resulth(3)+resulth(4)+resulth(5))/(-resulth(6)) -1d0))


        elseif(gauge.eq.3)then

        do i=0,3
          mup2test(i)=p2(i)
          mup3test(i)=p2(i)+p3(i)
        enddo


        call NoAbe2gaHexLineDiv(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -  ,result_color_g,Div)

        resultgauge(1)=(abs((resulth(1)+resulth(2)+resulth(3)+resulth(4)+resulth(5))/(-resulth(6)) -1d0))
        
         do i=0,3
          mup2test(i)=p2(i)+p3(i)
          mup3test(i)=p3(i)
        enddo

        call NoAbe2gaHexLineDiv(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -   ,result_color_g,Div)
        
        resultgauge(2)=(abs((resulth(1)+resulth(2)+resulth(3)+resulth(4)+resulth(5))/(-resulth(6)) -1d0))


        endif
        
        
        end


c gauge test of Hexagon using Pentagon: result_gaugeHex(1,1) and (1,2) should be the results of the Hexagons with mu2->p2 and mu3->p3; ngluons=2
       subroutine  gauge_NoAbe2gaHexLine(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musqIn,result_gaugeHex,resultgauge,Div,gauge)
       IMPLICIT NONE
#include "VBFNLO/loops/Utilities/GaugeTest.inc"
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup5(0:3),mup6(0:3)
       Complex*16 barpsi_p4(2),psi_p1(2)
       complex*16 result_gaugeHex(2,2)
       real*8 p23(0:3)
       Integer alpha ,mu
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result_graph(3),result_color(2)
       real*8 resultgauge,pengauge(2)
       Integer Div,Gauge
        musq=musqIn

      do mu=0,3
        p23(mu) = p2(mu)+p3(mu)
      enddo
      resultgauge = 0d0

      if(gauge.eq.1.or.gauge.eq.3) then
	call NoAbe2gaPenLineTotal(p1,p23,p4,p5,p6,barpsi_p4,psi_p1,mup3,mup5,mup6,
     &             alpha,musq,2,0,1,result_graph(1),result_color(1),pengauge,div,1)
	if(pengauge(1).gt.GaugeAccuracyHexLnA) then
	  resultgauge = pengauge(1)
	  return
	endif
	resultgauge = abs(-result_color(1)/result_gaugeHex(1,1)-1d0)
      endif


      if(gauge.eq.2.or.gauge.eq.3) then
	call NoAbe2gaPenLineTotal(p1,p23,p4,p5,p6,barpsi_p4,psi_p1,mup2,mup5,mup6,
     &             alpha,musq,2,0,1,result_graph(1),result_color(1),pengauge,div,1)
	if(pengauge(1).gt.GaugeAccuracyHexLnA) then
	  resultgauge = pengauge(1)
	  return
	endif
	resultgauge = max(abs(result_color(1)/result_gaugeHex(1,2)-1d0),resultgauge)
      endif

      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c                                QUAD
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       subroutine NoAbe2gaHexLineTotalQ(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musqIn,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,Div,gauge)
       IMPLICIT NONE
#include "VBFNLO/loops/Utilities/GaugeTest.inc"
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup5(0:3),mup6(0:3)

       Complex*16 barpsi_p4(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result_graph(6),result_color(2)
       real*8 resultgauge(2)
       Integer Div,Gauge
ccccccccccccccccccccccccccccccccccccccccccccccc 
       logical printGaugePrec
       parameter (printGaugePrec = .false.)
       real*8 comGaugePrec
       common /VVJJGaugePrec/ comGaugePrec

       integer jtemp
       common/ hexgaugeValue/t
       real*8 t,tdble,tquad      
       real*8 Accuracy
       Complex*16 result_graph_Q(6),result_color_Q(2)
       real*8 resultgauge_Q(2)
       logical wardpAbe(2),gaugetrue
       common/wardHexNoAbe/wardpAbe

        musq=musqIn
        tdble=1d8
        tquad=1d8
        
        gaugetrue = gauge.gt.0.and.gauge.le.3

      
        resultgauge_Q(1)=1d8
        resultgauge_Q(2)=1d8
        resultgauge(1)=1d8
        resultgauge(2)=1d8
        
       
        Accuracy=GaugeAccuracyHexLnA

        call NoAbe2gaHexLineTotal(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musqIn,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,Div,gauge)

        if(gaugetrue) then
c test NAN

      if(resultgauge(1)+1.eq.resultgauge(1))then
         resultgauge(1)=1d9
       endif

      if(resultgauge(2)+1.eq.resultgauge(2))then
        resultgauge(2)=1d9
       endif

      if(resultgauge(1).ne.resultgauge(1))then
         resultgauge(1)=1d9
       endif

      if(resultgauge(2).ne.resultgauge(2))then
        resultgauge(2)=1d9
       endif

        tdble=Max(resultgauge(1),resultgauge(2))

	if(tdble+1.eq.tdble) then
	  tdble=1d9
	  print*,"NaN in line", __LINE__ , __FILE__
	endif

	if(tdble.ne.tdble) then
	  tdble=1d9
	  print*,"NaN in line", __LINE__ , __FILE__
	endif


        wardpAbe(1) = tdble.lt.Accuracy
        wardpAbe(2) = tdble.lt.Accuracy

       if(tdble+1.eq.tdble) then
          tdble=1d9
       endif   

       if(tdble.ne.tdble) then
          tdble=1d9
       endif   


c       print*,"Inside new"
c       if(tdble.gt.Accuracy) then
c       print*, "gauge Accu Dble:2g", tdble
c       endif
#ifdef WITH_QUAD
       if(tdble.gt.Accuracy) then

        call NoAbe2gaHexLineTotal_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musqIn,ngluon,posgluon,comp,result_gr
     -   aph_Q,result_color_Q,resultgauge_Q,Div,gauge)

c test NAN

      if(resultgauge_Q(1)+1.eq.resultgauge_Q(1))then
         resultgauge_Q(1)=1d9
       endif

      if(resultgauge_Q(2)+1.eq.resultgauge_Q(2))then
         resultgauge_Q(2)=1d9
       endif


      if(resultgauge_Q(1).ne.resultgauge_Q(1))then
         resultgauge_Q(1)=1d9
       endif

      if(resultgauge_Q(2).ne.resultgauge_Q(2))then
         resultgauge_Q(2)=1d9
       endif


        tQUAD=Max(resultgauge_Q(1),resultgauge_Q(2))

	if(tquad+1.eq.tquad) then
	  tquad=1d9
	  print*,"NaN in line", __LINE__ , __FILE__
	endif

	if(tquad.ne.tquad) then
	  tquad=1d9
	  print*,"NaN in line", __LINE__ , __FILE__
	endif

       if(tquad+1.eq.tquad) then
          tquad=1d9
          print*,"QUAD:HexNo2g"
       endif   

       if(tquad.ne.tquad) then
          tquad=1d9
          print*,"QUAD:HexNo2g"
       endif   


c       if(tdble.gt.Accuracy) then
c       print*, "gauge Accu QUAD:2g", tquad
c       do jtemp=1,6
c       print*, "Dble", resultf(jtemp)
c       print*, "QUAD", resultf_Q(jtemp)
c       print*, "ratio:2g", abs(result_graph(jtemp)/result_graph_Q(jtemp) -1)
c       enddo
c       endif ! printing
       wardpAbe(2) = tquad.lt.accuracy

c to rescue Dble precision
       if (tquad.lt.tdble) then
        do jtemp=1,6
           result_graph(jtemp)=result_graph_Q(jtemp)
           if(jtemp.le.2) then
           result_color(jtemp)=result_color_Q(jtemp)
           resultgauge(jtemp)=resultgauge_Q(jtemp)
           endif
        enddo
        else
c       print*, "HexD:2g", tdble             
c       print*, "HexQ:2g", tquad

       endif

       endif  ! tdble


#endif
      endif ! gaugetrue
      end


#ifdef WITH_QUAD


       subroutine NoAbe2gaHexLineTotal_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musqIn,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,Div,gauge)
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup5(0:3),mup6(0:3)

       Complex*16 barpsi_p4(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result_graph(6),result_color(2)
       real*8 resultgauge(2)
       Integer Div,Gauge
        musq=musqIn

       if(Div.eq.0) then

        call NoAbe2gaHexLineT_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,gauge)

        else

        print*,"NO QUAD DIV so far:2g"
        call NoAbe2gaHexLineDivT(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,Div,gauge)


        endif


      end


        subroutine  NoAbe2gaHexLineT_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musqIn,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,gauge)

       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup5(0:3),mup6(0:3)
       Complex*16   mup2test(0:3),mup3test(0:3)
       Complex*16 barpsi_p4(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon,i 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result_graph(6),result_color(2),result_color_g(2)
     &            ,resulth(6)
       real*8 resultgauge(2)
       Integer Gauge
        musq=musqIn


        call NoAbe2gaHexLine_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color)
        
        if(gauge.eq.1)then

        do i=0,3
          mup2test(i)=p2(i)
          mup3test(i)=p2(i)+p3(i)
        enddo


        call NoAbe2gaHexLine_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -  ,result_color_g)

        resultgauge(1)=(abs((resulth(1)+resulth(2)+resulth(3)+resulth(4)+resulth(5))/(-resulth(6)) -1d0))
        
        elseif(gauge.eq.2)then

        do i=0,3
          mup2test(i)=p2(i)+p3(i)
          mup3test(i)=p3(i)
        enddo

        call NoAbe2gaHexLine_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -   ,result_color_g)
        
        resultgauge(2)=(abs((resulth(1)+resulth(2)+resulth(3)+resulth(4)+resulth(5))/(-resulth(6)) -1d0))


        elseif(gauge.eq.3)then

        do i=0,3
          mup2test(i)=p2(i)
          mup3test(i)=p2(i)+p3(i)
        enddo


        call NoAbe2gaHexLine_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -  ,result_color_g)

        resultgauge(1)=(abs((resulth(1)+resulth(2)+resulth(3)+resulth(4)+resulth(5))/(-resulth(6)) -1d0))
        
         do i=0,3
          mup2test(i)=p2(i)+p3(i)
          mup3test(i)=p3(i)
        enddo

        call NoAbe2gaHexLine_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -   ,result_color_g)
        
        resultgauge(2)=(abs((resulth(1)+resulth(2)+resulth(3)+resulth(4)+resulth(5))/(-resulth(6)) -1d0))


        endif
        
        
        end


#endif
