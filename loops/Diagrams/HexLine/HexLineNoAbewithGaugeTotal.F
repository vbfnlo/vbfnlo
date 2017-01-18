
       subroutine HexLineNoAbe1gTotal(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup3,mup4,
     1  mup6,alpha,musqIn,ngluon,posgluon,comp,
     2 result,resultg,Div,Gauge)

       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup4(0:3),mup6(0:3)

       Complex*16 barpsi_p5(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result(10),resultg(3)

       Integer Div,Gauge
       Complex*16 F(6395)
       Real*8 P(7129)
       COMMON/test_FhlFunctionsNoABE1g/F
       Save/test_FhlFunctionsNoABE1g/
       COMMON/test_PFunctionsNoABE1g/P
       musq=musqIn

       if(Div.eq.0) then

        call HexLineNoAbe1gT(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup3,mup4,
     1  mup6,alpha,musq,ngluon,posgluon,comp,
     2 result,resultg,Gauge)

        else
        call HexLineNoAbe1gDivT(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup3,mup4,
     1  mup6,alpha,musq,ngluon,posgluon,comp,
     2 result,resultg,Div,Gauge)
        endif


      end


       subroutine HexLineNoAbe1gTotalQ(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup3,mup4,
     1  mup6,alpha,musqIn,ngluon,posgluon,comp,
     2 result,resultg,Div,Gauge)

       IMPLICIT NONE
#include "VBFNLO/loops/Utilities/GaugeTest.inc"
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup4(0:3),mup6(0:3)

       Complex*16 barpsi_p5(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result(10),resultg(3)

       Integer Div,Gauge
       Complex*16 F(6395)
       Real*8 P(7129)
       COMMON/test_FhlFunctionsNoABE1g/F
       Save/test_FhlFunctionsNoABE1g/
       COMMON/test_PFunctionsNoABE1g/P

       Logical WardH
       logical printGaugePrec
       parameter (printGaugePrec = .false.)
       real*8 comGaugePrec
       common /VVJJGaugePrec/ comGaugePrec

       integer jtemp
       common/ hexgaugeValue/t
       real*8 t,tdble,tquad

       Complex*16 result_Q(10),resultg_Q(3)
       real*8 Accuracy

       logical gaugetrue
       logical wardpAbe(2)
       common/wardHexNoAbe/wardpAbe

        tdble=1d8
        tquad=1d8
       wardh=.false.
       Accuracy=GaugeAccuracyHexLnA


      if(gauge.gt.0.and.gauge.lt.6) then
       gaugetrue=.true.
       else
       gaugetrue=.false.
       endif
c       print*, "comGaugePrec",comGaugePrec
c       pause
       musq=musqIn


       call HexLineNoAbe1gTotal(p1,p2,p3,p4,p5,p6,
     &            barpsi_p5,psi_p1,mup2,mup3,mup4,mup6,
     &            alpha,musq,ngluon,posgluon,comp,
     2 result,resultg,Div,Gauge)

 
      if(gaugetrue) then

       wardpabe(1)=.false.
       wardpabe(2)=.false.

      call  hexGauge1gNoAbe(Accuracy,p1,p2,p3,p4,p5,p6,
     &            barpsi_p5,psi_p1,mup2,mup3,mup4,mup6,
     &            alpha,musq,
     &            resultg(1),WardH,gauge)

       tdble=t

       if(tdble+1.eq.tdble) then
          tdble=1d9
       endif  

       if(tdble.ne.tdble) then
          tdble=1d9
       endif  


       if(tdble.gt.Accuracy) then
c       print*, "gauge Accu Dble:1g", t
c       print*, "WardH Dble:1g", WardH
       else
       wardpAbe(1)=.true.
       wardpAbe(2)=.true.
       endif
#ifdef WITH_QUAD
       if(tdble.gt.Accuracy) then


       call HexLineNoAbe1gTotal_QUAD(p1,p2,p3,p4,p5,p6,
     &            barpsi_p5,psi_p1,mup2,mup3,mup4,mup6,
     &            alpha,musq,ngluon,posgluon,comp,
     2 result_Q,resultg_Q,Div,Gauge)


       call  hexGauge1gNoAbe_QUAD(Accuracy,p1,p2,p3,p4,p5,p6,
     &            barpsi_p5,psi_p1,mup2,mup3,mup4,mup6,
     &            alpha,musq,
     &            resultg_Q(1),WardH,gauge)

       tquad=t

       if(tquad+1.eq.tquad) then
          tquad=1d9
          print*,"QUAD:Hex1g"
       endif 

       if(tquad.ne.tquad) then
          tquad=1d9
          print*,"QUAD:Hex1g"
       endif 


       if(tdble.gt.Accuracy) then
c       print*, "gauge Accu QUAD:1g", tquad
c       print*, "WardH QUAD:1g", WardH
c       do jtemp=1,10
c       print*, "Dble", resultf(jtemp)
c       print*, "QUAD", resultf_Q(jtemp)
c       print*, "ratio:1g", abs(result(jtemp)/result_Q(jtemp) -1)
c       enddo
       endif ! printing


c to rescue Dble precision
       if (tquad.lt.tdble) then

          if(tquad.lt.Accuracy) wardpAbe(2)=.true.

        do jtemp=1,10
           result(jtemp)=result_Q(jtemp)
           if(jtemp.le.3) then
           resultg(jtemp)=resultg_Q(jtemp)
           endif
        enddo
        else
c       print*, "HexD:1g", tdble             
c       print*, "HexQ:1g", tquad
       endif

       endif  ! tdble

#endif

       endif


      end

#ifdef WITH_QUAD
       subroutine HexLineNoAbe1gTotal_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup3,mup4,
     1  mup6,alpha,musqIn,ngluon,posgluon,comp,
     2 result,resultg,Div,Gauge)

       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup4(0:3),mup6(0:3)

       Complex*16 barpsi_p5(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result(10),resultg(3)

       Integer Div,Gauge
       Complex*16 F(6395)
       Real*8 P(7129)
       COMMON/test_FhlFunctionsNoABE1g/F
       Save/test_FhlFunctionsNoABE1g/
       COMMON/test_PFunctionsNoABE1g/P
       musq=musqIn

       if(Div.eq.0) then

        call HexLineNoAbe1gT_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup3,mup4,
     1  mup6,alpha,musq,ngluon,posgluon,comp,
     2 result,resultg,Gauge)

        else

        print*,"NO QUAD DIV so far:1g"
        call HexLineNoAbe1gDivT(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup3,mup4,
     1  mup6,alpha,musq,ngluon,posgluon,comp,
     2 result,resultg,Div,Gauge)
        endif


      end


#endif
