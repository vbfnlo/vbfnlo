       subroutine hexlineABET(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,
     1  mup2,mup3,mup4,mup5,alpha,musqIn,ngluon,posgluon,comp,
     2 resultf,resultborn,resultg,resultbg,Gauge)
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup4(0:3),mup5(0:3)
       Complex*16   mup2test(0:3),mup3test(0:3),mup4test(0:3),mup5test(0:3)
       Complex*16 barpsi_p6(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon,i,Gauge 
       Complex*16 resultf(13),resultfg(13),resultborn,resultg(4)
       Complex*16 resultbg(4)
       Complex*16 borntemp,F(5930)
       Real*8 P(11674)
       COMMON/test_FhlFunctionsABE/F
       Save/test_FhlFunctionsABE/
       COMMON/test_PFunctionsABE/P
       Real*8 musqIn,musq
       Common/musqInv/musq
       musq=musqIn
       do i=0,3
          mup2test(i)=p2(i)
          mup3test(i)=p3(i)
          mup4test(i)=p4(i)
          mup5test(i)=p5(i)
      enddo
       call test_hexlineABE(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,mup2,mup3,mup4,
     1 mup5,alpha,musq,ngluon,posgluon,comp,resultf,resultborn)
 
C************************************************8
C**********  DBLE GAUGES *******************************8
C************************************************8

	if(Gauge.eq.1) then

       call test_hexlineABE(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,mup2test,mup3,mup4,
     1 mup5,alpha,musq,0,0,-1,resultfg,borntemp)
 
       resultg(1)=resultfg(1)
       resultbg(1)=borntemp

       elseif(Gauge.eq.2) then

       call test_hexlineABE(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,mup2,mup3test,mup4,
     1 mup5,alpha,musq,0,0,-1,resultfg,borntemp)
 
       resultg(2)=resultfg(1)
       resultbg(2)=borntemp

       elseif(Gauge.eq.3) then

       call test_hexlineABE(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,mup2,mup3,mup4test,
     1 mup5,alpha,musq,0,0,-1,resultfg,borntemp)
 
       resultg(3)=resultfg(1)
       resultbg(3)=borntemp

       elseif(Gauge.eq.4) then

       call test_hexlineABE(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,mup2,mup3,mup4,
     1 mup5test,alpha,musq,0,0,-1,resultfg,borntemp)
 
       resultg(4)=resultfg(1)
       resultbg(4)=borntemp

       
       elseif(Gauge.gt.4) then

       call test_hexlineABE(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,mup2test,mup3,mup4,
     1 mup5,alpha,musq,0,0,-1,resultfg,borntemp)
 
       resultg(1)=resultfg(1)
       resultbg(1)=borntemp


       call test_hexlineABE(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,mup2,mup3test,mup4,
     1 mup5,alpha,musq,0,0,-1,resultfg,borntemp)
 
       resultg(2)=resultfg(1)
       resultbg(2)=borntemp


       call test_hexlineABE(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,mup2,mup3,mup4test,
     1 mup5,alpha,musq,0,0,-1,resultfg,borntemp)
 
       resultg(3)=resultfg(1)
       resultbg(3)=borntemp


       call test_hexlineABE(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,mup2,mup3,mup4,
     1 mup5test,alpha,musq,0,0,-1,resultfg,borntemp)
 
       resultg(4)=resultfg(1)
       resultbg(4)=borntemp

	endif


      return
      end


#if WITH_QUAD


       subroutine hexlineABET_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,
     1  mup2,mup3,mup4,mup5,alpha,musqIn,ngluon,posgluon,comp,
     2 resultf,resultborn,resultg,resultbg,Gauge)
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup4(0:3),mup5(0:3)
       Complex*16   mup2test(0:3),mup3test(0:3),mup4test(0:3),mup5test(0:3)
       Complex*16 barpsi_p6(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon,i,Gauge 
       Complex*16 resultf(13),resultfg(13),resultborn,resultg(4)
       Complex*16 resultbg(4)
       Complex*16 borntemp,F(5930)
       Real*8 P(11674)
       COMMON/test_FhlFunctionsABE/F
       Save/test_FhlFunctionsABE/
       COMMON/test_PFunctionsABE/P
       Real*8 musqIn,musq
       Common/musqInv/musq
       musq=musqIn
       do i=0,3
          mup2test(i)=p2(i)
          mup3test(i)=p3(i)
          mup4test(i)=p4(i)
          mup5test(i)=p5(i)
      enddo
       call test_hexlineABE_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,mup2,mup3,mup4,
     1 mup5,alpha,musq,ngluon,posgluon,comp,resultf,resultborn)
 
C************************************************8
C**********  DBLE GAUGES *******************************8
C************************************************8

	if(Gauge.eq.1) then

       call test_hexlineABE_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,mup2test,mup3,mup4,
     1 mup5,alpha,musq,0,0,-1,resultfg,borntemp)
 
       resultg(1)=resultfg(1)
       resultbg(1)=borntemp

       elseif(Gauge.eq.2) then

       call test_hexlineABE_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,mup2,mup3test,mup4,
     1 mup5,alpha,musq,0,0,-1,resultfg,borntemp)
 
       resultg(2)=resultfg(1)
       resultbg(2)=borntemp

       elseif(Gauge.eq.3) then

       call test_hexlineABE_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,mup2,mup3,mup4test,
     1 mup5,alpha,musq,0,0,-1,resultfg,borntemp)
 
       resultg(3)=resultfg(1)
       resultbg(3)=borntemp

       elseif(Gauge.eq.4) then

       call test_hexlineABE_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,mup2,mup3,mup4,
     1 mup5test,alpha,musq,0,0,-1,resultfg,borntemp)
 
       resultg(4)=resultfg(1)
       resultbg(4)=borntemp

       
       elseif(Gauge.gt.4) then

       call test_hexlineABE_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,mup2test,mup3,mup4,
     1 mup5,alpha,musq,0,0,-1,resultfg,borntemp)
 
       resultg(1)=resultfg(1)
       resultbg(1)=borntemp


       call test_hexlineABE_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,mup2,mup3test,mup4,
     1 mup5,alpha,musq,0,0,-1,resultfg,borntemp)
 
       resultg(2)=resultfg(1)
       resultbg(2)=borntemp


       call test_hexlineABE_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,mup2,mup3,mup4test,
     1 mup5,alpha,musq,0,0,-1,resultfg,borntemp)
 
       resultg(3)=resultfg(1)
       resultbg(3)=borntemp


       call test_hexlineABE_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,mup2,mup3,mup4,
     1 mup5test,alpha,musq,0,0,-1,resultfg,borntemp)
 
       resultg(4)=resultfg(1)
       resultbg(4)=borntemp

	endif


      return
      end

#endif



