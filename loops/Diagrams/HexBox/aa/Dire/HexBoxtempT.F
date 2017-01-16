c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 18/01/2010
c Modified:9/3/2012
c ******************************************************************************
c hel = 66,67,76,77 for rr,rl,...
c or 2*hel1 + hel2 -> 3,1,-1,-3

       subroutine HexBoxT(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn,Div,hel)
       implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16 barpsi_p1(4),psi_p6(4),barpsi_p2(4),psi_p5(4)
       Complex*16   mup3(0:3),mup4(0:3)
       Real*8 musq
       Integer comp,ngluon,posgluon,Div,hel
       Complex*16 result(6),resultn,resultborn

       if (hel.eq.77.or.hel.eq.-3) then
         call HexBox77T(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn,Div)
       elseif (hel.eq.67.or.hel.eq.1) then
         call HexBox67T(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn,Div)
       elseif (hel.eq.76.or.hel.eq.-1) then
         call HexBox76T(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn,Div)
       elseif (hel.eq.66.or.hel.eq.3) then
         call HexBox66T(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn,Div)
       else
	 print*,'wrong helicity selected in HexBoxtempT.F'
       endif
       Return
       End


       subroutine HexBox66T(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn,Div)
       implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16 barpsi_p1(4),psi_p6(4),barpsi_p2(4),psi_p5(4)
       Complex*16   mup3(0:3),mup4(0:3)
       Real*8 musq
       Integer comp,ngluon,posgluon,Div
       Complex*16 result(6),resultn,resultborn

       if(Div.eq.0) then

         call HexBox66(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn)
          else

        call HexBox66Div(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn,Div)

       endif

       Return
       End


       subroutine HexBox77T(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn,Div)
       implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16 barpsi_p1(4),psi_p6(4),barpsi_p2(4),psi_p5(4)
       Complex*16   mup3(0:3),mup4(0:3)
       Real*8 musq
       Integer comp,ngluon,posgluon,Div
       Complex*16 result(6),resultn,resultborn
c************************************************************************************

       if(Div.eq.0) then

         call HexBox77(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn)

          else

        call HexBox77Div(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn,Div)

       endif

       Return
       End


       subroutine HexBox67T(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn,Div)
       implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16 barpsi_p1(4),psi_p6(4),barpsi_p2(4),psi_p5(4)
       Complex*16   mup3(0:3),mup4(0:3)
       Real*8 musq
       Integer comp,ngluon,posgluon,Div
       Complex*16 result(6),resultn,resultborn
c************************************************************************************
c************************************************************************************

       if(Div.eq.0) then

         call HexBox67(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn)

          else

        call HexBox67Div(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn,Div)

       endif

       Return
       End


       subroutine HexBox76T(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn,Div)
       implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16 barpsi_p1(4),psi_p6(4),barpsi_p2(4),psi_p5(4)
       Complex*16   mup3(0:3),mup4(0:3)
       Real*8 musq
       Integer comp,ngluon,posgluon,Div
       Complex*16 result(6),resultn,resultborn
       if(Div.eq.0) then

         call HexBox76(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn)

          else

        call HexBox76Div(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn,Div)

       endif

       Return
       End
