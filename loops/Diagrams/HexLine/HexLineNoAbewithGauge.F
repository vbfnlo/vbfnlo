
       subroutine HexLineNoAbe1gT(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup3,mup4,
     1  mup6,alpha,musqIn,ngluon,posgluon,comp,
     2 result,resultg,Gauge)

       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup4(0:3),mup6(0:3)
       Complex*16   mup2test(0:3),mup3test(0:3),mup4test(0:3)
       Complex*16 barpsi_p5(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon,i ,gauge
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result(10),resultgt(10),resultg(3)
       Complex*16 borntemp
       Complex*16 F(6395)
       Real*8 P(7129)
       COMMON/test_FhlFunctionsNoABE1g/F
       Save/test_FhlFunctionsNoABE1g/
       COMMON/test_PFunctionsNoABE1g/P
       musq=musqIn


       do i=0,3
          mup2test(i)=p2(i)
          mup3test(i)=p3(i)
          mup4test(i)=p4(i)
      enddo

      call test_hexlineNoAbe1g(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup3,mup4,
     -mup6,alpha,musq,ngluon,posgluon,comp,result,borntemp)

C************************************************8
C**********  DBLE GAUGES *******************************8
C********************************************************

      if(gauge.eq.1) then

       call test_hexlineNoAbe1g(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2test,mup3,mup4,
     - mup6,alpha,musq,1,9,-1,resultgt,borntemp)

      resultg(1)=resultgt(1)
      
      elseif(gauge.eq.2) then

       call test_hexlineNoAbe1g(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup3test,mup4,
     - mup6,alpha,musq,1,9,-1,resultgt,borntemp)

      resultg(2)=resultgt(1)

      elseif(gauge.eq.3) then

       call test_hexlineNoAbe1g(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup3,mup4test,
     - mup6,alpha,musq,1,9,-1,resultgt,borntemp)

      resultg(3)=resultgt(1)

  
      elseif(gauge.gt.3) then
      
      
       call test_hexlineNoAbe1g(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2test,mup3,mup4,
     - mup6,alpha,musq,1,9,-1,resultgt,borntemp)

      resultg(1)=resultgt(1)


       call test_hexlineNoAbe1g(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup3test,mup4,
     - mup6,alpha,musq,1,9,-1,resultgt,borntemp)

      resultg(2)=resultgt(1)


       call test_hexlineNoAbe1g(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup3,mup4test,
     - mup6,alpha,musq,1,9,-1,resultgt,borntemp)

      resultg(3)=resultgt(1)

      endif
      RETURN


      end


#if WITH_QUAD

       subroutine HexLineNoAbe1gT_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup3,mup4,
     1  mup6,alpha,musqIn,ngluon,posgluon,comp,
     2 result,resultg,Gauge)

       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup4(0:3),mup6(0:3)
       Complex*16   mup2test(0:3),mup3test(0:3),mup4test(0:3)
       Complex*16 barpsi_p5(2),psi_p1(2)
       Integer comp,alpha, ngluon, posgluon,i ,gauge
       Real*8 musqIn,musq
       Common/musqInv/musq
       Complex*16 result(10),resultgt(10),resultg(3)
       Complex*16 borntemp
       Complex*16 F(6395)
       Real*8 P(7129)
       COMMON/test_FhlFunctionsNoABE1g/F
       Save/test_FhlFunctionsNoABE1g/
       COMMON/test_PFunctionsNoABE1g/P
       musq=musqIn


       do i=0,3
          mup2test(i)=p2(i)
          mup3test(i)=p3(i)
          mup4test(i)=p4(i)
      enddo

      call test_hexlineNoAbe1g_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup3,mup4,
     -mup6,alpha,musq,ngluon,posgluon,comp,result,borntemp)

C************************************************8
C**********  DBLE GAUGES *******************************8
C********************************************************

      if(gauge.eq.1) then

       call test_hexlineNoAbe1g_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2test,mup3,mup4,
     - mup6,alpha,musq,1,9,-1,resultgt,borntemp)

      resultg(1)=resultgt(1)
      
      elseif(gauge.eq.2) then

       call test_hexlineNoAbe1g_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup3test,mup4,
     - mup6,alpha,musq,1,9,-1,resultgt,borntemp)

      resultg(2)=resultgt(1)

      elseif(gauge.eq.3) then

       call test_hexlineNoAbe1g_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup3,mup4test,
     - mup6,alpha,musq,1,9,-1,resultgt,borntemp)

      resultg(3)=resultgt(1)

  
      elseif(gauge.gt.3) then
      
      
       call test_hexlineNoAbe1g_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2test,mup3,mup4,
     - mup6,alpha,musq,1,9,-1,resultgt,borntemp)

      resultg(1)=resultgt(1)


       call test_hexlineNoAbe1g_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup3test,mup4,
     - mup6,alpha,musq,1,9,-1,resultgt,borntemp)

      resultg(2)=resultgt(1)


       call test_hexlineNoAbe1g_QUAD(p1,p2,p3,p4,p5,p6,barpsi_p5,psi_p1,mup2,mup3,mup4test,
     - mup6,alpha,musq,1,9,-1,resultgt,borntemp)

      resultg(3)=resultgt(1)

      endif
      RETURN


      end


#endif
