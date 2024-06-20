c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Modified:4/6/2012

c hel = 66,67,76,77 for rr,rl,...
c or 2*hel1 + hel2 -> 3,1,-1,-3

       subroutine PenBoxCrossT(p1,p2,p3,p4,p5,psi_p1,barpsi_p5,barpsi
     &   _p2,psi_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,res
     &   ultborn,Div,hel)
       implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 psi_p1(4),barpsi_p5(4),barpsi_p2(4),psi_p4(4)
       Complex*16   mup3(0:3)
       Real*8 musq
       integer comp,ngluon,posgluon,div,hel
       complex*16 result(3),resultn,resultborn
       
       if(hel.eq.77.or.hel.eq.-3) then
	call PenBoxCross77T(p1,p2,p3,p4,p5,psi_p1,barpsi_p5,barpsi
     &   _p2,psi_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,res
     &   ultborn,Div)
       elseif(hel.eq.67.or.hel.eq.1) then
	call PenBoxCross67T(p1,p2,p3,p4,p5,psi_p1,barpsi_p5,barpsi
     &   _p2,psi_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,res
     &   ultborn,Div)
       elseif(hel.eq.76.or.hel.eq.-1) then
	call PenBoxCross76T(p1,p2,p3,p4,p5,psi_p1,barpsi_p5,barpsi
     &   _p2,psi_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,res
     &   ultborn,Div)
       elseif(hel.eq.66.or.hel.eq.3) then
	call PenBoxCross66T(p1,p2,p3,p4,p5,psi_p1,barpsi_p5,barpsi
     &   _p2,psi_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,res
     &   ultborn,Div)
       else
	 print*,'wrong helicity selected in HexBoxtempT.F'
       endif
       Return
       end


       subroutine PenBoxCross66T(p1,p2,p3,p4,p5,psi_p1,barpsi_p5,barpsi
     &   _p2,psi_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,res
     &   ultborn,Div)
       implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 psi_p1(4),barpsi_p5(4),barpsi_p2(4),psi_p4(4)
       Complex*16   mup3(0:3)
       Real*8 musq
       integer comp,ngluon,posgluon,div
       complex*16 result(3),resultn,resultborn
       
       if (Div.ne.0) then

       call PenBoxCross66Div(p1,p2,p3,p4,p5,psi_p1,barpsi_p5,barpsi
     &   _p2,psi_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,res
     &   ultborn,Div)
       else

       call PenBoxCross66(p1,p2,p3,p4,p5,psi_p1,barpsi_p5,barpsi
     &   _p2,psi_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,res
     &   ultborn)
       endif
       end


       subroutine PenBoxCross67T(p1,p2,p3,p4,p5,psi_p1,barpsi_p5,barpsi
     &   _p2,psi_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,res
     &   ultborn,Div)
       implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 psi_p1(4),barpsi_p5(4),barpsi_p2(4),psi_p4(4)
       Complex*16   mup3(0:3)
       Real*8 musq
       integer comp,ngluon,posgluon,div
       complex*16 result(3),resultn,resultborn
       
       if (Div.ne.0) then

       call PenBoxCross67Div(p1,p2,p3,p4,p5,psi_p1,barpsi_p5,barpsi
     &   _p2,psi_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,res
     &   ultborn,Div)
       else

       call PenBoxCross67(p1,p2,p3,p4,p5,psi_p1,barpsi_p5,barpsi
     &   _p2,psi_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,res
     &   ultborn)
       endif
       end


       subroutine PenBoxCross76T(p1,p2,p3,p4,p5,psi_p1,barpsi_p5,barpsi
     &   _p2,psi_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,res
     &   ultborn,Div)
       implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 psi_p1(4),barpsi_p5(4),barpsi_p2(4),psi_p4(4)
       Complex*16   mup3(0:3)
       Real*8 musq
       integer comp,ngluon,posgluon,div
       complex*16 result(3),resultn,resultborn
       
       if (Div.ne.0) then

       call PenBoxCross76Div(p1,p2,p3,p4,p5,psi_p1,barpsi_p5,barpsi
     &   _p2,psi_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,res
     &   ultborn,Div)
       else

       call PenBoxCross76(p1,p2,p3,p4,p5,psi_p1,barpsi_p5,barpsi
     &   _p2,psi_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,res
     &   ultborn)
       endif
       end


       subroutine PenBoxCross77T(p1,p2,p3,p4,p5,psi_p1,barpsi_p5,barpsi
     &   _p2,psi_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,res
     &   ultborn,Div)
       implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 psi_p1(4),barpsi_p5(4),barpsi_p2(4),psi_p4(4)
       Complex*16   mup3(0:3)
       Real*8 musq
       integer comp,ngluon,posgluon,div
       complex*16 result(3),resultn,resultborn
       
       if (Div.ne.0) then

       call PenBoxCross77Div(p1,p2,p3,p4,p5,psi_p1,barpsi_p5,barpsi
     &   _p2,psi_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,res
     &   ultborn,Div)
       else

       call PenBoxCross77(p1,p2,p3,p4,p5,psi_p1,barpsi_p5,barpsi
     &   _p2,psi_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,res
     &   ultborn)
       endif
       end
