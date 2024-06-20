c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Modified:4/6/2012
c hel = 66,67,76,77 for rr,rl,...
c or 2*hel1 + hel2 -> 3,1,-1,-3
       subroutine PenBoxT(p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2,p
     &   si_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,resultbo
     &   rn,Div,hel)
       Implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p1(4),psi_p5(4),barpsi_p2(4),psi_p4(4)
       Complex*16   mup3(0:3)
       real*8 musq
       integer comp,ngluon,posgluon,Div,hel
       complex*16 result(3),resultn,resultborn
 
       if( hel.eq.77.or.hel.eq.-3) then
	 call PenBox77T(p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2,p
     &   si_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,resultbo
     &   rn,Div)
       elseif( hel.eq.67.or.hel.eq.1) then
	 call PenBox67T(p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2,p
     &   si_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,resultbo
     &   rn,Div)
       elseif( hel.eq.76.or.hel.eq.-1) then
	 call PenBox76T(p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2,p
     &   si_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,resultbo
     &   rn,Div)
       elseif( hel.eq.66.or.hel.eq.3) then
	 call PenBox66T(p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2,p
     &   si_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,resultbo
     &   rn,Div)
       else	 
	print*,'wrong helicity selected in PenBoxtempT.F'
       endif
       
        end


       subroutine PenBox66T(p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2,p
     &   si_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,resultbo
     &   rn,Div)
       Implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p1(4),psi_p5(4),barpsi_p2(4),psi_p4(4)
       Complex*16   mup3(0:3)
       real*8 musq
       integer comp,ngluon,posgluon,Div
       complex*16 result(3),resultn,resultborn
 
 
       if(Div.ne.0) then

       call PenBox66Div(p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2,p
     &   si_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,resultbo
     &   rn,Div)

       else
    
        call PenBox66(p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2,p
     &   si_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,resultbo
     &   rn)

        endif
        end


       subroutine PenBox67T(p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2,p
     &   si_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,resultbo
     &   rn,Div)
       Implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p1(4),psi_p5(4),barpsi_p2(4),psi_p4(4)
       Complex*16   mup3(0:3)
       real*8 musq
       integer comp,ngluon,posgluon,Div
       complex*16 result(3),resultn,resultborn
 
 
       if(Div.ne.0) then

       call PenBox67Div(p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2,p
     &   si_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,resultbo
     &   rn,Div)

       else
    
        call PenBox67(p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2,p
     &   si_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,resultbo
     &   rn)

        endif
        end


       subroutine PenBox76T(p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2,p
     &   si_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,resultbo
     &   rn,Div)
       Implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p1(4),psi_p5(4),barpsi_p2(4),psi_p4(4)
       Complex*16   mup3(0:3)
       real*8 musq
       integer comp,ngluon,posgluon,Div
       complex*16 result(3),resultn,resultborn
 
 
       if(Div.ne.0) then

       call PenBox76Div(p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2,p
     &   si_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,resultbo
     &   rn,Div)

       else
    
        call PenBox76(p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2,p
     &   si_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,resultbo
     &   rn)

        endif
        end


       subroutine PenBox77T(p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2,p
     &   si_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,resultbo
     &   rn,Div)
       Implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p1(4),psi_p5(4),barpsi_p2(4),psi_p4(4)
       Complex*16   mup3(0:3)
       real*8 musq
       integer comp,ngluon,posgluon,Div
       complex*16 result(3),resultn,resultborn
 
 
       if(Div.ne.0) then

       call PenBox77Div(p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2,p
     &   si_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,resultbo
     &   rn,Div)

       else
    
        call PenBox77(p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2,p
     &   si_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,resultbo
     &   rn)

        endif
        end
