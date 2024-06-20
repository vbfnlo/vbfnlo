       subroutine hexlineABETotal(p1,p2,p3,p4,p5,p6,
     &     barpsi_p6,psi_p1,
     &     mup2,mup3,mup4,mup5,alpha,musqIn,
     &     ngluon,posgluon,comp,
     &     resultf,resultborn,resultg,resultbg,
     &     Div,Gauge)
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup4(0:3),mup5(0:3)

       Complex*16 barpsi_p6(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 resultf(13),resultborn,resultg(4)
       Complex*16 resultbg(4)
       Complex*16 F(5930)
       Real*8 P(11674)
       COMMON/test_FhlFunctionsABE/F
       Save/test_FhlFunctionsABE/
       COMMON/test_PFunctionsABE/P
       Integer Div,Gauge


       musq=musqIn

       If(Div.eq.0) then

        call hexlineABET(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,
     1  mup2,mup3,mup4,mup5,alpha,musq,ngluon,posgluon,comp,
     2 resultf,resultborn,resultg,resultbg,Gauge)


       else
       
        call hexlineABEDivT(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,
     1  mup2,mup3,mup4,mup5,alpha,musq,ngluon,posgluon,comp,
     2 resultf,resultborn,resultg,resultbg,Div,Gauge)

       endif    
        
       end

cfc
       subroutine hexlineABETotalQ(p1,p2,p3,p4,p5,p6,
     &     barpsi_p6,psi_p1,
     &     mup2,mup3,mup4,mup5,alpha,musqIn,
     &     ngluon,posgluon,comp,
     &     resultf,resultborn,resultg,resultbg,
     &     Div,Gauge)
       IMPLICIT NONE
#include "VBFNLO/loops/Utilities/GaugeTest.inc"
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup4(0:3),mup5(0:3)

       Complex*16 barpsi_p6(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 resultf(13),resultborn,resultg(4)
       Complex*16 resultbg(4)
       Complex*16 F(5930)
       Real*8 P(11674)
       COMMON/test_FhlFunctionsABE/F
       Save/test_FhlFunctionsABE/
       COMMON/test_PFunctionsABE/P
       Integer Div,Gauge
       Logical WardH
       logical printGaugePrec
       parameter (printGaugePrec = .false.)
       real*8 comGaugePrec
       common /VVJJGaugePrec/ comGaugePrec

       Complex*16 resultf_Q(13),resultborn_Q,resultg_Q(4)
       Complex*16 resultbg_Q(4)
       integer jtemp
       common/ hexgaugeValue/t
       real*8 t,tdble,tquad
       logical gaugetrue
       logical wardpAbe(2)
       common/wardHexAbe/wardpAbe       

        tdble=1d8
        tquad=1d8


       if(gauge.gt.0.and.gauge.lt.6) then
       gaugetrue=.true.
       else
       gaugetrue=.false.
       endif
 
c       wardh=.false.

        musq=musqIn

       call hexlineABETotal(p1,p2,p3,p4,p5,p6,
     &     barpsi_p6,psi_p1,
     &     mup2,mup3,mup4,mup5,alpha,musq,
     &     ngluon,posgluon,comp,
     &     resultf,resultborn,resultg,resultbg,
     &     Div,Gauge)

      if(gaugetrue) then

       wardpAbe(1)=.false.
       wardpAbe(2)=.false.  

       call hexGaugeAbe(p1,p2,p3,p4,p5,p6,
     &       barpsi_p6,psi_p1,    
     &       mup2,mup3,mup4,mup5,alpha,musq,
     &       resultg,resultbg,WardH,gauge)


       tdble=t
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


       if(tdble.gt.GaugeAccuracyHexL) then
c       print*, "gauge Accu Dble", t
       else
       wardpAbe(1)=.true.
       wardpAbe(2)=.true.
       endif
#ifdef WITH_QUAD
       if(tdble.gt.GaugeAccuracyHexL) then
       call hexlineABETotal_QUAD(p1,p2,p3,p4,p5,p6,
     &     barpsi_p6,psi_p1,
     &     mup2,mup3,mup4,mup5,alpha,musq,
     &     ngluon,posgluon,comp,
     &     resultf_Q,resultborn_Q,resultg_Q,resultbg_Q,
     &     Div,Gauge)

  
       call hexGaugeAbe_QUAD(p1,p2,p3,p4,p5,p6,
     &       barpsi_p6,psi_p1,    
     &       mup2,mup3,mup4,mup5,alpha,musq,
     &       resultg_Q,resultbg_Q,WardH,gauge)
       tquad=t
       if(tquad+1.eq.tquad) then
          tquad=1d9
          print*,"QUAD:Hex"
       endif 

       if(tquad.ne.tquad) then
          tquad=1d9
          print*,"QUAD:Hex"
       endif 


      if(tdble.gt.GaugeAccuracyHexL) then
c       print*, "gauge Accu QUAD", tquad
c       print*, "WardH QUAD", WardH
c       do jtemp=1,13
c       print*, "Dble", resultf(jtemp)
c       print*, "QUAD", resultf_Q(jtemp)
c       print*, "ratio", abs(resultf(jtemp)/resultf_Q(jtemp) -1)
c       enddo
       endif ! printing

c to rescue Dble precision
       if (tquad.lt.tdble) then

         if(tquad.lt.GaugeAccuracyHexL) wardpAbe(2)=.true.


          resultborn= resultborn_Q
        do jtemp=1,13
           resultf(jtemp)=resultf_Q(jtemp)
           if(jtemp.le.4) then
           resultg(jtemp)=resultg_Q(jtemp)
           resultbg(jtemp)  = resultbg_Q(jtemp)
           endif
        enddo
        else
c       print*, "HexD:0g", tdble             
c       print*, "HexQ:0g", tquad

       endif

       endif  ! tdble
       
#endif

       endif ! gaugetrue
       end


#ifdef WITH_QUAD
       subroutine hexlineABETotal_QUAD(p1,p2,p3,p4,p5,p6,
     &     barpsi_p6,psi_p1,
     &     mup2,mup3,mup4,mup5,alpha,musqIn,
     &     ngluon,posgluon,comp,
     &     resultf,resultborn,resultg,resultbg,
     &     Div,Gauge)
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup4(0:3),mup5(0:3)

       Complex*16 barpsi_p6(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon 
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 resultf(13),resultborn,resultg(4)
       Complex*16 resultbg(4)
       Complex*16 F(5930)
       Real*8 P(11674)
       COMMON/test_FhlFunctionsABE/F
       Save/test_FhlFunctionsABE/
       COMMON/test_PFunctionsABE/P
       Integer Div,Gauge


       musq=musqIn

       If(Div.eq.0) then

        call hexlineABET_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,
     1  mup2,mup3,mup4,mup5,alpha,musq,ngluon,posgluon,comp,
     2 resultf,resultborn,resultg,resultbg,Gauge)


       else
       print*,"NO QUAD DIV so far:HexLine"
        call hexlineABEDivT(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,
     1  mup2,mup3,mup4,mup5,alpha,musq,ngluon,posgluon,comp,
     2 resultf,resultborn,resultg,resultbg,Div,Gauge)

       endif    
        
       end
#endif
