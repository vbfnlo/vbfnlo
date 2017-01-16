c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
       subroutine HexBoxCross1T(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p6,
     &   barpsi_p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,
     &   resultn,resultborn,Div,hel)
       implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16 psi_p1(4),barpsi_p6(4),barpsi_p2(4),psi_p5(4)
       Complex*16   mup3(0:3),mup4(0:3)
       Real*8 musq
       Integer comp,ngluon,posgluon,Div,hel
       Complex*16 result(6),resultn,resultborn
       
       if(hel.eq.77.or.hel.eq.-3) then
	 call HexBoxCross177T(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p6,
     &   barpsi_p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,
     &   resultn,resultborn,Div)
       elseif(hel.eq.67.or.hel.eq.1) then
	 call HexBoxCross167T(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p6,
     &   barpsi_p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,
     &   resultn,resultborn,Div)
       elseif(hel.eq.76.or.hel.eq.-1) then
	 call HexBoxCross176T(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p6,
     &   barpsi_p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,
     &   resultn,resultborn,Div)
       elseif(hel.eq.66.or.hel.eq.3) then
	 call HexBoxCross166T(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p6,
     &   barpsi_p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,
     &   resultn,resultborn,Div)
       else
	 print*,'wrong helicity selected in HexBoxCross1tempT.F'
       endif

       Return
       End


       subroutine HexBoxCross166T(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p6,
     &   barpsi_p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,
     &   resultn,resultborn,Div)
       implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16 psi_p1(4),barpsi_p6(4),barpsi_p2(4),psi_p5(4)
       Complex*16   mup3(0:3),mup4(0:3)
       Real*8 musq
       Integer comp,ngluon,posgluon,Div
       Complex*16 result(6),resultn,resultborn

       if(Div.eq.0) then
          call HexBoxCross166_b(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p6,
     &   barpsi_p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,
     &   resultn,resultborn)


          else 
         call  HexBoxCross166Div_b(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p6,
     &   barpsi_p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,
     &   resultn,resultborn,Div)

       endif

       Return
       End


       subroutine HexBoxCross167T(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p6,
     &   barpsi_p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,
     &   resultn,resultborn,Div)
       implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16 psi_p1(4),barpsi_p6(4),barpsi_p2(4),psi_p5(4)
       Complex*16   mup3(0:3),mup4(0:3)
       Real*8 musq
       Integer comp,ngluon,posgluon,Div
       Complex*16 result(6),resultn,resultborn

       if(Div.eq.0) then
          call HexBoxCross167_b(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p6,
     &   barpsi_p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,
     &   resultn,resultborn)


          else 
         call  HexBoxCross167Div_b(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p6,
     &   barpsi_p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,
     &   resultn,resultborn,Div)

       endif

       Return
       End


       subroutine HexBoxCross176T(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p6,
     &   barpsi_p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,
     &   resultn,resultborn,Div)
       implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16 psi_p1(4),barpsi_p6(4),barpsi_p2(4),psi_p5(4)
       Complex*16   mup3(0:3),mup4(0:3)
       Real*8 musq
       Integer comp,ngluon,posgluon,Div
       Complex*16 result(6),resultn,resultborn

       if(Div.eq.0) then
          call HexBoxCross176_b(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p6,
     &   barpsi_p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,
     &   resultn,resultborn)


          else 
         call  HexBoxCross176Div_b(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p6,
     &   barpsi_p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,
     &   resultn,resultborn,Div)

       endif

       Return
       End


       subroutine HexBoxCross177T(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p6,
     &   barpsi_p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,
     &   resultn,resultborn,Div)
       implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16 psi_p1(4),barpsi_p6(4),barpsi_p2(4),psi_p5(4)
       Complex*16   mup3(0:3),mup4(0:3)
       Real*8 musq
       Integer comp,ngluon,posgluon,Div
       Complex*16 result(6),resultn,resultborn

       if(Div.eq.0) then
          call HexBoxCross177_b(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p6,
     &   barpsi_p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,
     &   resultn,resultborn)


          else 
         call  HexBoxCross177Div_b(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p6,
     &   barpsi_p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,
     &   resultn,resultborn,Div)

       endif

       Return
       End
