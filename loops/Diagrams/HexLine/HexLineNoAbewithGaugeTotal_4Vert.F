
       subroutine NoAbe4vertHexLineTotal(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musqIn,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,Div,gauge)
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup5(0:3),mup6(0:3)

       Complex*16 barpsi_p4(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result_graph(6,3),result_color(3)
       real*8 resultgauge(2)
       Integer Div,Gauge
        musq=musqIn

       if(Div.eq.0) then

        call NoAbe4vertHexLineT(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,gauge)

        else
        call NoAbe4vertHexLineDivT(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,Div,gauge)


        endif


      end


        subroutine  NoAbe4vertHexLineT(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
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
       Complex*16 result_graph(6,3),result_color(3),result_color_g(3)
     &            ,resulth(6,3)
       real*8 resultgauge(2)
       Integer Gauge
       complex*16 resulttemp
       integer k
       real*8 td(3)
        musq=musqIn


        call NoAbe4vertHexLine(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color)
        
        if(gauge.eq.1)then

        do i=0,3
          mup2test(i)=p2(i)
          mup3test(i)=p2(i)+p3(i)
        enddo


        call NoAbe4vertHexLine(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -  ,result_color_g)


       do k=1,3
       resulttemp=0
       td(k)=1d0

       if (ABS(resulth(1,k)).gt.1d-4)then
       do i=2,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp/(-resulth(1,k)) -1d0)
       
       ELSE
       do i=1,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp)
       ENDIF
       enddo

       resultgauge(1)=max(td(1),td(2),td(3))
        
        elseif(gauge.eq.2)then

        do i=0,3
          mup2test(i)=p2(i)+p3(i)
          mup3test(i)=p3(i)
        enddo

        call NoAbe4vertHexLine(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -   ,result_color_g)


       do k=1,3
       resulttemp=0
       td(k)=1d0

       if (ABS(resulth(1,k)).gt.1d-4)then
       do i=2,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp/(-resulth(1,k)) -1d0)
       
       ELSE
       do i=1,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp)
       ENDIF
       enddo

       resultgauge(2)=max(td(1),td(2),td(3))


        elseif(gauge.eq.3)then

        do i=0,3
          mup2test(i)=p2(i)
          mup3test(i)=p2(i)+p3(i)
        enddo


        call NoAbe4vertHexLine(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -  ,result_color_g)


        do k=1,3
       resulttemp=0
       td(k)=1d0

       if (ABS(resulth(1,k)).gt.1d-4)then
       do i=2,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp/(-resulth(1,k)) -1d0)
       
       ELSE
       do i=1,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp)
       ENDIF
       enddo

       resultgauge(1)=max(td(1),td(2),td(3))

         do i=0,3
          mup2test(i)=p2(i)+p3(i)
          mup3test(i)=p3(i)
        enddo

        call NoAbe4vertHexLine(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -   ,result_color_g)
 

      do k=1,3
       resulttemp=0
       td(k)=1d0

       if (ABS(resulth(1,k)).gt.1d-4)then
       do i=2,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp/(-resulth(1,k)) -1d0)
       
       ELSE
       do i=1,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp)
       ENDIF
       enddo

       resultgauge(2)=max(td(1),td(2),td(3))


        endif
        
        
        end


        subroutine  NoAbe4vertHexLineDivT(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musqIn,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,div,gauge)

       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup5(0:3),mup6(0:3)
       Complex*16   mup2test(0:3),mup3test(0:3)
       Complex*16 barpsi_p4(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon,i 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result_graph(6,3),result_color(3),result_color_g(3)
     &            ,resulth(6,3)
       real*8 resultgauge(2)
       Integer Div,Gauge
       complex*16 resulttemp
       integer k
       real*8 td(3)
        musq=musqIn


        call NoAbe4vertHexLineDiv(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color,Div)
        
        if(gauge.eq.1)then

        do i=0,3
          mup2test(i)=p2(i)
          mup3test(i)=p2(i)+p3(i)
        enddo


        call NoAbe4vertHexLineDiv(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -  ,result_color_g,Div)

       do k=1,3
       resulttemp=0
       td(k)=1d0

       if (ABS(resulth(1,k)).gt.1d-4)then
       do i=2,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp/(-resulth(1,k)) -1d0)
       
       ELSE
       do i=1,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp)
       ENDIF
       enddo


       resultgauge(1)=max(td(1),td(2),td(3))


        elseif(gauge.eq.2)then

        do i=0,3
          mup2test(i)=p2(i)+p3(i)
          mup3test(i)=p3(i)
        enddo

        call NoAbe4vertHexLineDiv(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -   ,result_color_g,Div)


       do k=1,3
       resulttemp=0
       td(k)=1d0

       if (ABS(resulth(1,k)).gt.1d-4)then
       do i=2,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp/(-resulth(1,k)) -1d0)
       
       ELSE
       do i=1,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp)
       ENDIF
       enddo


       resultgauge(2)=max(td(1),td(2),td(3))


        elseif(gauge.eq.3)then

        do i=0,3
          mup2test(i)=p2(i)
          mup3test(i)=p2(i)+p3(i)
        enddo


        call NoAbe4vertHexLineDiv(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -  ,result_color_g,Div)


       do k=1,3
       resulttemp=0
       td(k)=1d0

       if (ABS(resulth(1,k)).gt.1d-4)then
       do i=2,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp/(-resulth(1,k)) -1d0)
       
       ELSE
       do i=1,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp)
       ENDIF
       enddo


       resultgauge(1)=max(td(1),td(2),td(3))
        
         do i=0,3
          mup2test(i)=p2(i)+p3(i)
          mup3test(i)=p3(i)
        enddo

        call NoAbe4vertHexLineDiv(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -   ,result_color_g,Div)


       do k=1,3
       resulttemp=0
       td(k)=1d0

       if (ABS(resulth(1,k)).gt.1d-4)then
       do i=2,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp/(-resulth(1,k)) -1d0)
       
       ELSE
       do i=1,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp)
       ENDIF
       enddo

       resultgauge(2)=max(td(1),td(2),td(3))

        endif
        
        
        end


c gauge test of Hexagon using Pentagon: result_gaugeHex(1,1) and (1,2) should be the results of the Hexagons with mu2->p2 and mu3->p3; ngluons=2
       subroutine  gauge_NoAbe4VertHexLine(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musqIn,result_gaugeColorHex,resultgauge,Div,gauge)
       IMPLICIT NONE
#include "VBFNLO/loops/Utilities/GaugeTest.inc"
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup5(0:3),mup6(0:3)
       Complex*16 barpsi_p4(2),psi_p1(2)
       complex*16 result_gaugeColorHex(3,2)
       real*8 p23(0:3)
       Integer alpha ,mu
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result_graph(3,3),result_color(3)
       real*8 resultgauge,pengauge(2)
       Integer Div,Gauge,i
        musq=musqIn

      do mu=0,3
        p23(mu) = p2(mu)+p3(mu)
      enddo
      resultgauge = 0d0

      if(gauge.eq.1.or.gauge.eq.3) then
	call NoAbe4VertBoxlineTotal(p1,p23,p4,p5,p6,barpsi_p4,psi_p1,mup3,mup5,mup6,
     &             alpha,musq,2,0,1,result_graph(1,1),result_color(1),pengauge,div,1)
      
	if(pengauge(1).gt.GaugeAccuracyHexLnA) then
	  resultgauge = pengauge(1)
	  return
	endif
	do i=1,3
	  resultgauge = max(abs(result_color(i)/result_gaugeColorHex(i,1)+1d0),resultgauge)
	enddo
      endif

      if(gauge.eq.2.or.gauge.eq.3) then
	call NoAbe4VertBoxlineTotal(p1,p23,p4,p5,p6,barpsi_p4,psi_p1,mup2,mup5,mup6,
     &             alpha,musq,2,0,1,result_graph(1,1),result_color(1),pengauge,div,1)
	if(pengauge(1).gt.GaugeAccuracyHexLnA) then
	  resultgauge = pengauge(1)
	  return
	endif
	do i=1,3
	  resultgauge = max(abs(result_color(i)/result_gaugeColorHex(i,2)-1d0),resultgauge)
	enddo
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


       subroutine NoAbe4vertHexLineTotalQ(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
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
       Complex*16 result_graph(6,3),result_color(3)
       real*8 resultgauge(2)
       Integer Div,Gauge
cccccccccccccccccccccccccccccccccccccccccccccccccc
       logical printGaugePrec, gaugetrue
       parameter (printGaugePrec = .false.)


       integer jtemp
       common/ hexgaugeValue/t
       real*8 t,tdble,tquad      
       real*8 Accuracy
       Complex*16 result_graph_Q(6,3),result_color_Q(3)
       real*8 resultgauge_Q(2)
       logical wardpAbe(2)
       common/wardHexNoAbe/wardpAbe

        musq=musqIn
        tdble=1d8
        tquad=1d8
        
        gaugetrue = gauge.gt.0 .and. gauge.le.3
        Accuracy=GaugeAccuracyHexLnA
        resultgauge_Q(1)=1d8
        resultgauge_Q(2)=1d8
        resultgauge(1)=1d8
        resultgauge(2)=1d8


          call NoAbe4vertHexLineTotal(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
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


          call NoAbe4vertHexLineTotal_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
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

        wardpAbe(2) = tquad.lt.Accuracy

       if(tquad+1.eq.tquad) then
          tquad=1d9
       endif   

       if(tquad.ne.tquad) then
          tquad=1d9
       endif

c to rescue Dble precision
       if (tquad.lt.tdble) then

        do jtemp=1,6
           result_graph(jtemp,1)=result_graph_Q(jtemp,1)
           result_graph(jtemp,2)=result_graph_Q(jtemp,2)
           result_graph(jtemp,3)=result_graph_Q(jtemp,3)
        enddo
  
          
           result_color(1)=result_color_Q(1)
           result_color(2)=result_color_Q(2)
           result_color(3)=result_color_Q(3)

           resultgauge(1)=resultgauge_Q(1)
           resultgauge(2)=resultgauge_Q(2)
        else
c       print*, "HexD:4g", tdble             
c       print*, "HexQ:4g", tquad
c       print*, "WardH QUAD:1g",
       endif

       endif  ! tdble


#endif
      endif ! gaugetrue


          end subroutine


#ifdef WITH_QUAD


       subroutine NoAbe4vertHexLineTotal_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musqIn,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,Div,gauge)
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup5(0:3),mup6(0:3)

       Complex*16 barpsi_p4(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result_graph(6,3),result_color(3)
       real*8 resultgauge(2)
       Integer Div,Gauge
        musq=musqIn

       if(Div.eq.0) then

        call NoAbe4vertHexLineT_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,gauge)

        else
        print*,"NO QUAD DIV so far:4g"
        call NoAbe4vertHexLineDivT(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,Div,gauge)


        endif


      end


        subroutine  NoAbe4vertHexLineT_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
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
       Complex*16 result_graph(6,3),result_color(3),result_color_g(3)
     &            ,resulth(6,3)
       real*8 resultgauge(2)
       Integer Gauge
       complex*16 resulttemp
       integer k
       real*8 td(3)
        musq=musqIn


        call NoAbe4vertHexLine_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color)
        
        if(gauge.eq.1)then

        do i=0,3
          mup2test(i)=p2(i)
          mup3test(i)=p2(i)+p3(i)
        enddo


        call NoAbe4vertHexLine_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -  ,result_color_g)


       do k=1,3
       resulttemp=0
       td(k)=1d0

       if (ABS(resulth(1,k)).gt.1d-4)then
       do i=2,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp/(-resulth(1,k)) -1d0)
       
       ELSE
       do i=1,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp)
       ENDIF
       enddo

       resultgauge(1)=max(td(1),td(2),td(3))
        
        elseif(gauge.eq.2)then

        do i=0,3
          mup2test(i)=p2(i)+p3(i)
          mup3test(i)=p3(i)
        enddo

        call NoAbe4vertHexLine_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -   ,result_color_g)


       do k=1,3
       resulttemp=0
       td(k)=1d0

       if (ABS(resulth(1,k)).gt.1d-4)then
       do i=2,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp/(-resulth(1,k)) -1d0)
       
       ELSE
       do i=1,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp)
       ENDIF
       enddo

       resultgauge(2)=max(td(1),td(2),td(3))


        elseif(gauge.eq.3)then

        do i=0,3
          mup2test(i)=p2(i)
          mup3test(i)=p2(i)+p3(i)
        enddo


        call NoAbe4vertHexLine_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -  ,result_color_g)


        do k=1,3
       resulttemp=0
       td(k)=1d0

       if (ABS(resulth(1,k)).gt.1d-4)then
       do i=2,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp/(-resulth(1,k)) -1d0)
       
       ELSE
       do i=1,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp)
       ENDIF
       enddo

       resultgauge(1)=max(td(1),td(2),td(3))

         do i=0,3
          mup2test(i)=p2(i)+p3(i)
          mup3test(i)=p3(i)
        enddo

        call NoAbe4vertHexLine_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
     -   ,result_color_g)
 

      do k=1,3
       resulttemp=0
       td(k)=1d0

       if (ABS(resulth(1,k)).gt.1d-4)then
       do i=2,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp/(-resulth(1,k)) -1d0)
       
       ELSE
       do i=1,6
          resulttemp=resulttemp+resulth(i,k)
       enddo
       td(k)=abs(resulttemp)
       ENDIF
       enddo

       resultgauge(2)=max(td(1),td(2),td(3))


        endif
        
        
        end


cFC        subroutine  NoAbe4vertHexLineDivT_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
cFC     -   up2,mup3,mup5,mup6,alpha,musqIn,ngluon,posgluon,comp,result_gr
cFC     -   aph,result_color,resultgauge,div,gauge)
cFC
cFC       IMPLICIT NONE
cFC       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
cFC       Complex*16   mup2(0:3),mup3(0:3),mup5(0:3),mup6(0:3)
cFC       Complex*16   mup2test(0:3),mup3test(0:3),mup5test(0:3),mup6test(0:3)
cFC       Complex*16 barpsi_p4(2),psi_p1(2)
cFC       Integer comp,alpha, ngluon, posgluon,i 
cFC       Real*8 musqIn,musq
cFC       Common/musqInv/musq
cFC       Complex*16 result_graph(6,3),result_color(3),result_color_g(3)
cFC     &            ,resulth(6,3)
cFC       real*8 resultgauge(2)
cFC       Integer Div,Gauge
cFC       complex*16 resulttemp
cFC       integer k
cFC       real*8 td(3)
cFC        musq=musqIn
cFC
cFC
cFC        call NoAbe4vertHexLineDiv_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
cFC     -   up2,mup3,mup5,mup6,alpha,musq,ngluon,posgluon,comp,result_gr
cFC     -   aph,result_color,Div)
cFC        
cFC        if(gauge.eq.1)then
cFC
cFC        do i=0,3
cFC          mup2test(i)=p2(i)
cFC          mup3test(i)=p2(i)+p3(i)
cFC        enddo
cFC
cFC
cFC        call NoAbe4vertHexLineDiv_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
cFC     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
cFC     -  ,result_color_g,Div)
cFC
cFC       do k=1,3
cFC       resulttemp=0
cFC       td(k)=1d0
cFC
cFC       if (ABS(resulth(1,k)).gt.1d-4)then
cFC       do i=2,6
cFC          resulttemp=resulttemp+resulth(i,k)
cFC       enddo
cFC       td(k)=abs(resulttemp/(-resulth(1,k)) -1d0)
cFC       
cFC       ELSE
cFC       do i=1,6
cFC          resulttemp=resulttemp+resulth(i,k)
cFC       enddo
cFC       td(k)=abs(resulttemp)
cFC       ENDIF
cFC       enddo
cFC
cFC
cFC       resultgauge(1)=max(td(1),td(2),td(3))
cFC
cFC
cFC
cFC
cFC        elseif(gauge.eq.2)then
cFC
cFC        do i=0,3
cFC          mup2test(i)=p2(i)+p3(i)
cFC          mup3test(i)=p3(i)
cFC        enddo
cFC
cFC        call NoAbe4vertHexLineDiv_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
cFC     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
cFC     -   ,result_color_g,Div)
cFC        
cFC
cFC
cFC
cFC       do k=1,3
cFC       resulttemp=0
cFC       td(k)=1d0
cFC
cFC       if (ABS(resulth(1,k)).gt.1d-4)then
cFC       do i=2,6
cFC          resulttemp=resulttemp+resulth(i,k)
cFC       enddo
cFC       td(k)=abs(resulttemp/(-resulth(1,k)) -1d0)
cFC       
cFC       ELSE
cFC       do i=1,6
cFC          resulttemp=resulttemp+resulth(i,k)
cFC       enddo
cFC       td(k)=abs(resulttemp)
cFC       ENDIF
cFC       enddo
cFC
cFC
cFC
cFC       resultgauge(2)=max(td(1),td(2),td(3))
cFC
cFC
cFC
cFC        elseif(gauge.eq.3)then
cFC
cFC        do i=0,3
cFC          mup2test(i)=p2(i)
cFC          mup3test(i)=p2(i)+p3(i)
cFC        enddo
cFC
cFC
cFC        call NoAbe4vertHexLineDiv_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
cFC     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
cFC     -  ,result_color_g,Div)
cFC
cFC  
cFC
cFC       do k=1,3
cFC       resulttemp=0
cFC       td(k)=1d0
cFC
cFC       if (ABS(resulth(1,k)).gt.1d-4)then
cFC       do i=2,6
cFC          resulttemp=resulttemp+resulth(i,k)
cFC       enddo
cFC       td(k)=abs(resulttemp/(-resulth(1,k)) -1d0)
cFC       
cFC       ELSE
cFC       do i=1,6
cFC          resulttemp=resulttemp+resulth(i,k)
cFC       enddo
cFC       td(k)=abs(resulttemp)
cFC       ENDIF
cFC       enddo
cFC
cFC
cFC       resultgauge(1)=max(td(1),td(2),td(3))
cFC        
cFC         do i=0,3
cFC          mup2test(i)=p2(i)+p3(i)
cFC          mup3test(i)=p3(i)
cFC        enddo
cFC
cFC        call NoAbe4vertHexLineDiv_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
cFC     -   up2test,mup3test,mup5,mup6,alpha,musq,ngluon,posgluon,-1,resulth
cFC     -   ,result_color_g,Div)
cFC        
cFC 
cFC
cFC
cFC       do k=1,3
cFC       resulttemp=0
cFC       td(k)=1d0
cFC
cFC       if (ABS(resulth(1,k)).gt.1d-4)then
cFC       do i=2,6
cFC          resulttemp=resulttemp+resulth(i,k)
cFC       enddo
cFC       td(k)=abs(resulttemp/(-resulth(1,k)) -1d0)
cFC       
cFC       ELSE
cFC       do i=1,6
cFC          resulttemp=resulttemp+resulth(i,k)
cFC       enddo
cFC       td(k)=abs(resulttemp)
cFC       ENDIF
cFC       enddo
cFC
cFC       resultgauge(2)=max(td(1),td(2),td(3))
cFC
cFC        endif
cFC        
cFC        
cFC        end


c gauge test of Hexagon using Pentagon: result_gaugeHex(1,1) and (1,2) should be the results of the Hexagons with mu2->p2 and mu3->p3; ngluons=2
       subroutine  gauge_NoAbe4VertHexLine_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musqIn,result_gaugeColorHex,resultgauge,Div,gauge)
       IMPLICIT NONE
#include "VBFNLO/loops/Utilities/GaugeTest.inc"
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup5(0:3),mup6(0:3)
       Complex*16 barpsi_p4(2),psi_p1(2)
       complex*16 result_gaugeColorHex(3,2)
       real*8 p23(0:3)
       Integer alpha ,mu
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result_graph(3,3),result_color(3)
       real*8 resultgauge,pengauge(2)
       Integer Div,Gauge,i
        musq=musqIn

      do mu=0,3
        p23(mu) = p2(mu)+p3(mu)
      enddo
      resultgauge = 0d0

      if(gauge.eq.1.or.gauge.eq.3) then
	call NoAbe4VertBoxlineTotal_QUAD(p1,p23,p4,p5,p6,barpsi_p4,psi_p1,mup3,mup5,mup6,
     &             alpha,musq,2,0,1,result_graph(1,1),result_color(1),pengauge,div,1)
      
	if(pengauge(1).gt.GaugeAccuracyHexLnA) then
	  resultgauge = pengauge(1)
	  return
	endif
	do i=1,3
	  resultgauge = max(abs(result_color(i)/result_gaugeColorHex(i,1)+1d0),resultgauge)
	enddo
      endif

      if(gauge.eq.2.or.gauge.eq.3) then
	call NoAbe4VertBoxlineTotal_QUAD(p1,p23,p4,p5,p6,barpsi_p4,psi_p1,mup2,mup5,mup6,
     &             alpha,musq,2,0,1,result_graph(1,1),result_color(1),pengauge,div,1)
	if(pengauge(1).gt.GaugeAccuracyHexLnA) then
	  resultgauge = pengauge(1)
	  return
	endif
	do i=1,3
	  resultgauge = max(abs(result_color(i)/result_gaugeColorHex(i,2)-1d0),resultgauge)
	enddo
      endif

      end


#endif



