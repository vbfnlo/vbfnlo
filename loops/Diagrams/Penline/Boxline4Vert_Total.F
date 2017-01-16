
       subroutine NoAbe4vertBoxLineTotal(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,m
     -   up2,mup4,mup5,alpha,musqIn,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,Div,gauge)
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16   mup2(0:3),mup4(0:3),mup5(0:3)

       Complex*16 barpsi_p3(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result_graph(3,3),result_color(3)
       real*8 resultgauge(2)
       Integer Div,Gauge
        musq=musqIn


       if(Div.eq.0) then

       call NoAbe4vertBoxLineT(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,m
     -   up2,mup4,mup5,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,gauge)

        else

       call  NoAbe4vertBoxLineDivT(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,m
     -   up2,mup4,mup5,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,Div,gauge)

        endif

      end


       subroutine NoAbe4vertBoxLineT(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,m
     -   up2,mup4,mup5,alpha,musqIn,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,gauge)
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16   mup2(0:3),mup4(0:3),mup5(0:3)
       Complex*16   mup2test(0:3)
       Complex*16 barpsi_p3(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon,i 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result_graph(3,3),result_color(3)
       real*8 resultgauge(2)
       Integer Gauge 

       Complex*16 resulth(3,3),result_color_g(3)
       complex*16 resulttemp
       integer k
       real*8 td(3)
        musq=musqIn


        call NoAbe4vertBoxLine(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,m
     -   up2,mup4,mup5,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color)
        
        if(gauge.eq.1)then

        do i=0,3
          mup2test(i)=p2(i)
        enddo


        call NoAbe4vertBoxLine(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,m
     -   up2test,mup4,mup5,alpha,musq,ngluon,posgluon,-1,resulth
     -  ,result_color_g)


       do k=1,3
       resulttemp=0
       td(k)=1d0
       do i=2,3
          resulttemp=resulttemp+resulth(i,k)
      enddo
	if(abs(resulth(1,k)).gt.1d-5) then
	  td(k)=abs(resulttemp/(-resulth(1,k)) -1d0)
	else
	  td(k)=abs(resulttemp+resulth(1,k))
        endif
       enddo
       resultgauge(1)=max(td(1),td(2),td(3))
    

        endif
        
        
        end


       subroutine NoAbe4vertBoxLineDivT(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,m
     -   up2,mup4,mup5,alpha,musqIn,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,Div,gauge)
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16   mup2(0:3),mup4(0:3),mup5(0:3)
       Complex*16   mup2test(0:3)
       Complex*16 barpsi_p3(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon,i 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result_graph(3,3),result_color(3)
       real*8 resultgauge(2)
       Integer Div,Gauge 

       Complex*16 resulth(3,3),result_color_g(3)
       complex*16 resulttemp
       integer k
       real*8 td(3)
        musq=musqIn


        call NoAbe4vertBoxLineDiv(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,m
     -   up2,mup4,mup5,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color,Div)
        
        if(gauge.eq.1)then

        do i=0,3
          mup2test(i)=p2(i)
        enddo


        call NoAbe4vertBoxLineDiv(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,m
     -   up2test,mup4,mup5,alpha,musq,ngluon,posgluon,-1,resulth
     -  ,result_color_g,Div)


       do k=1,3
       resulttemp=0
       td(k)=1d0
       do i=2,3
          resulttemp=resulttemp+resulth(i,k)
      enddo
	if(abs(resulth(1,k)).gt.1d-5) then
	  td(k)=abs(resulttemp/(-resulth(1,k)) -1d0)
	else
	  td(k)=abs(resulttemp+resulth(1,k))
	endif
       enddo
       resultgauge(1)=max(td(1),td(2),td(3))
    

        endif
        
        
        end

       subroutine NoAbe4vertBoxLineTotalQ(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,m
     -   up2,mup4,mup5,alpha,musqIn,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,Div,gauge)
       IMPLICIT NONE
#include "VBFNLO/loops/Utilities/GaugeTest.inc"
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16   mup2(0:3),mup4(0:3),mup5(0:3)

       Complex*16 barpsi_p3(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result_graph(3,3),result_color(3)
       real*8 resultgauge(2)
       Integer Div,Gauge

       common/pengaugeValue/t
       real*8 t,tdble,tquad

       Complex*16 result_graph_Q(3,3),result_color_Q(3)
       real*8 resultgauge_Q(2)
       integer jtemp,jtemp2
       logical gaugetrue
       logical wardpAbe(2)
       common/wardpenNoAbe/wardpAbe

        musq=musqIn

       tdble=1d9
       tquad=1d9
       resultgauge_Q(1)=1d9
       resultgauge_Q(2)=1d9
       resultgauge(1)=1d9
       resultgauge(2)=1d9

c************************************************************************************

       if(gauge.eq.1) then
       gaugetrue=.true.
       else
       gaugetrue=.false.
       endif

       call  NoAbe4vertBoxLineTotal(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,m
     -   up2,mup4,mup5,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,Div,gauge)

        tdble=resultgauge(1)
       if(tdble+1.eq.tdble) then
	tdble = 1d9
	print*, "NaN in line", __LINE__, __FILE__
       endif

       if(tdble.ne.tdble) then
	tdble = 1d9
	print*, "NaN in line", __LINE__, __FILE__
       endif

       if(tdble+1.eq.tdble) then
          tdble=1d9
       endif   

       if(tdble.ne.tdble) then
          tdble=1d9
       endif   

       if(gaugetrue) then

       wardpAbe(1)=.false.
       wardpAbe(2)=.false.

        if(tdble.gt.GaugeAccuracyPenLnA) then
c        print*, "gauge Accu Dble:P4g", t
c       print*,  "WardP Dble", Ward
       else
       wardpAbe(1)=.true.
       wardpAbe(2)=.true.
       endif

#ifdef WITH_QUAD   

       if(tdble.gt.GaugeAccuracyPenLnA) then
 
       call  NoAbe4vertBoxLineTotal_QUAD(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,m
     -   up2,mup4,mup5,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph_Q,result_color_Q,resultgauge_Q,Div,gauge)

        tquad=resultgauge_Q(1)
       if(tquad+1.eq.tquad) then
	tquad = 1d9
	print*, "NaN in line", __LINE__, __FILE__
       endif

       if(tquad.ne.tquad) then
	tquad = 1d9
	print*, "NaN in line", __LINE__, __FILE__
       endif


       if(tquad+1.eq.tquad) then
          tquad=1d9
          print*,"QUAD:BoxNo4g"
       endif  

       if(tquad.ne.tquad) then
          tquad=1d9
          print*,"QUAD:BoxNo4g"
       endif  


       if(tdble.gt.GaugeAccuracyPenLnA) then
c          print*, "COMP",COMP
c          print*, "gaugetest",gaugetest
c       print*, "gauge Accu QUAD:P4g", tquad
c       do jtemp=1,8
c       print*, "Dble", resultf(jtemp)
c       print*, "QUAD", resultf_Q(jtemp)
c       print*, "ratio:P", abs(result(jtemp)/result_Q(jtemp) -1)
c       enddo
       endif ! printing

c to rescue Dble precision
       if (tquad.lt.tdble) then
         if(tquad.lt.GaugeAccuracyPenLnA) wardpAbe(2)=.true.


        do jtemp=1,3
           result_color(jtemp)=result_color_Q(jtemp)
           if(jtemp.le.2) then
                resultgauge(jtemp)=resultgauge_Q(jtemp)
           endif

           do jtemp2=1,3
             result_graph(jtemp,jtemp2)=result_graph_Q(jtemp,jtemp2)
           enddo 
        enddo

       endif

       endif  ! tdble
#endif


       endif ! gaugetrue


      end


#ifdef WITH_QUAD


      subroutine NoAbe4vertBoxLineTotal_QUAD(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,m
     -   up2,mup4,mup5,alpha,musqIn,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,Div,gauge)
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16   mup2(0:3),mup4(0:3),mup5(0:3)

       Complex*16 barpsi_p3(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result_graph(3,3),result_color(3)
       real*8 resultgauge(2)
       Integer Div,Gauge
        musq=musqIn


       if(Div.eq.0) then

       call NoAbe4vertBoxLineT_QUAD(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,m
     -   up2,mup4,mup5,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,gauge)

        else
       call  NoAbe4vertBoxLineDivT(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,m
     -   up2,mup4,mup5,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,Div,gauge)
           print*, "No 4g pen in QUAD"
!            stop
        endif

      end


       subroutine NoAbe4vertBoxLineT_QUAD(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,m
     -   up2,mup4,mup5,alpha,musqIn,ngluon,posgluon,comp,result_gr
     -   aph,result_color,resultgauge,gauge)
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16   mup2(0:3),mup4(0:3),mup5(0:3)
       Complex*16   mup2test(0:3)
       Complex*16 barpsi_p3(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon,i 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result_graph(3,3),result_color(3)
       real*8 resultgauge(2)
       Integer Gauge 

       Complex*16 resulth(3,3),result_color_g(3)
       complex*16 resulttemp
       integer k
       real*8 td(3)
        musq=musqIn


         call NoAbe4vertBoxLine_QUAD(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,m
     -   up2,mup4,mup5,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color)
        
        if(gauge.eq.1)then

        do i=0,3
          mup2test(i)=p2(i)
        enddo


        call NoAbe4vertBoxLine_QUAD(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,m
     -   up2test,mup4,mup5,alpha,musq,ngluon,posgluon,-1,resulth
     -  ,result_color_g)


       do k=1,3
       resulttemp=0
       td(k)=1d0
       do i=2,3
          resulttemp=resulttemp+resulth(i,k)
      enddo
	if(abs(resulth(1,k)).gt.1d-5) then
	  td(k)=abs(resulttemp/(-resulth(1,k)) -1d0)
	else
	  td(k)=abs(resulttemp+resulth(1,k))
        endif
       enddo
       resultgauge(1)=max(td(1),td(2),td(3))
    

        endif
        
        
        end


#endif
