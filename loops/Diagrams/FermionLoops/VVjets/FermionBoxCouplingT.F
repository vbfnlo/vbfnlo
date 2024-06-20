c Author: Francisco Campanario
c  Computes the Box fermion loop with generic couplings for the massless case and 
c the massive with only one mass. Only valid for Z,gamma, gluon, couplings.
c  za1,za2,za3,za4 are the P_Right couplings
c  zb1,zb2,zb3,zb4 are the P_left couplings

cMK when using the Peskin convention (i.e. D_mu=d_mu - igA corresponding to +ig in the vertices) 
cMK this contribution should get a POSITIVE sign


        subroutine FermionBoxCouplingT(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,p2,p3,p4,
     &  mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge,Div)
       implicit none
       Real*8 za1,za2,za3,za4,zb1,zb2,zb3,zb4
       Real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3)

       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3)
       Complex*16 result,resultE,resultO
       Real*8 musq,m,resultgauge,resultgaugeE,resultgaugeO
       integer gauge,comp,Div

       call  FermionBoxCouplingTE(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,p2,p3,p4,
     &  mup1,mup2,mup3,mup4,musq,gauge,comp,resultE,resultgaugeE,Div)
! 
! 
       call  FermionBoxCouplingTO(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,p2,p3,p4,
     &  mup1,mup2,mup3,mup4,musq,gauge,comp,resultO,resultgaugeO,Div)

c      print*, "resultE", resultE
c      print*, "resultO", resultO

cMK the minus sign for the odd part is needed due to the different conventions of eps_0123
       result=resultE-resultO

       resultgauge=Max(resultgaugeO,resultgaugeE)

       end


cccccccccccccccccccc    Even PART ccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine FermionBoxCouplingTE(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,p2,p3,p4,
     &  mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge,Div)
       implicit none
       Real*8 za1,za2,za3,za4,zb1,zb2,zb3,zb4
       Real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3)

       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3)

       Complex*16 result
       Real*8 musq,m,resultgauge

       integer gauge,comp,Div

       If (Div.eq.0) then
       call FermionBoxCoupling(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,p2,p3,p4
     & ,mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge)
       else
       call FermionBoxCoupling_Div(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,p2,
     & p3,p4,mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge,Div)
       endif
       end


       subroutine FermionBoxCoupling(za1,za2,za3,za4,zb1,zb2,zb3,zb4,
     & m,p1,p2,p3,p4,mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge)
       implicit none
#include "VBFNLO/loops/Utilities/GaugeTest.inc"
       Real*8 za1,za2,za3,za4,zb1,zb2,zb3,zb4
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3)
c       Real*8   p12(0:3),p14(0:3),p23(0:3),p34(0:3)
       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3)
       Complex*16 mup1test(0:3),mup2test(0:3),mup3test(0:3),mup4test(0:3)
       Complex*16 result,resultn(5),resultZERO(5)
       Real*8 musq,m,resultgauge,temp1,temp2,temp3,temp4
       real*8 resultgauge_cp
       complex*16 result_cp
       complex*16 rg1,rg2,rg3,rg4
       integer gauge,i,comp
       common/GramSwitch/GramCS,GramDS
       logical GramCS,GramDS
c Checking Gram massive


       real*8 dotrr
       external dotrr
       logical printmy
       common/mprint/printmy
       logical debug
       parameter (debug=.false.)
       real*8 Sensitivity, FinalSensitivity
c For QUAD
       Real*8 resultgauge_q
       complex*16 result_Q
       complex*16 resultn_Q(5)
       complex*16 rg1_Q,rg2_Q,rg3_Q
       complex*16 rg4_Q
       logical wardfb(2)
       common/wardpfermionbox/wardfb
       logical GaugeE,GaugeO
       common/FloopPenggzz/GaugeE,GaugeO
       logical debugodd
       parameter (debugodd=.false.)


       if(debugodd) then
          result=0d0
          resultgauge=0d0
          return
       endif

	wardfb(1)=.false.
        wardfb(2)=.false.

        Sensitivity=GaugeAccuracyBoxL
        FinalSensitivity=GaugeAccuracyBoxL
c        m=4.6d0

        printmy=.false.
        GramCS=.false.
        GramDS=.false.

        resultgauge=1d99

cfc         call  FermionBoxEvenCouplings_QUAD(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,
cfc     & p2,p3,p4,mup1,mup2,mup3,mup4,musq,comp,result_Q,resultn_Q)
                         
       call  FermionBoxEvenCouplings(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,
     & p2,p3,p4,mup1,mup2,mup3,mup4,musq,comp,result,resultn)


       if (gauge.gt.0) then
       
         do  i=0,3
            mup1test(i)=p1(i)
            mup2test(i)=p2(i)
            mup3test(i)=p3(i)
            mup4test(i)=p4(i)
         enddo
       
c Call for Gauge

          call fermionBoxEvenCouplings(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1test,mup2,mup3,mup4,musq,-1,rg1,resultZERO)
          call fermionBoxEvenCouplings(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1,mup2test,mup3,mup4,musq,-1,rg2,resultZERO)
          call fermionBoxEvenCouplings(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1,mup2,mup3test,mup4,musq,-1,rg3,resultZERO)
          call fermionBoxEvenCouplings(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1,mup2,mup3,mup4test,musq,-1,rg4,resultZERO)
c Gauge 1
          if (abs(rg1).gt.1d-7) then
          temp1=abs(rg1/(resultn(2)) -1d0)
          else
          temp1=abs(rg1-resultn(2))   
          endif

c Gauge 2
         if (abs(rg2).gt.1d-7) then
          temp2=abs(rg2/resultn(3) -1d0) 
          else
          temp2=abs(rg2 -resultn(3)  ) 
          endif

c Gauge 3
         if (abs(rg3).gt.1d-7) then
          temp3=abs(rg3/resultn(4) -1d0 )
          else
             temp3=abs(rg3  -resultn(4) )
         endif

c Gauge 4
          if (abs(rg4).gt.1d-7) then
          temp4=abs(rg4/resultn(5) -1d0)
          else
          temp4=abs(rg4 - resultn(5)  )
          endif


       if(temp1+1.eq.temp1) then
          temp1=1d9
       endif   

       if(temp1.ne.temp1) then
          temp1=1d9
       endif

       if(temp2+1.eq.temp2) then
          temp2=1d9
       endif   

       if(temp2.ne.temp2) then
          temp2=1d9
       endif

       if(temp3+1.eq.temp3) then
          temp3=1d9
       endif   

       if(temp3.ne.temp3) then
          temp3=1d9
       endif

       if(temp4+1.eq.temp4) then
          temp4=1d9
       endif   

       if(temp4.ne.temp4) then
          temp4=1d9
       endif


          resultgauge=max(temp1,temp2,temp3,temp4)
        
c Safe results
          result_cp=result
          resultgauge_cp=resultgauge

cfccCCCCCCCCCCCc TO RESCUE
cfc
cfc        if(.false.) then
        if(resultgauge.gt.Sensitivity) then
cfc
cfc        GramCS=.true.
cfc        GramDS=.true.
cfc        resultgauge=1d99
cfc
cfc       call  FermionBoxEvenCouplings(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,
cfc     & p2,p3,p4,mup1,mup2,mup3,mup4,musq,comp,result,resultn)
cfc       
cfcc Call for Gauge
cfc          call fermionBoxEvenCouplings(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
cfc     &     p1,p2,p3,p4,mup1test,mup2,mup3,mup4,musq,-1,rg1,resultZERO)
cfc          call fermionBoxEvenCouplings(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
cfc     &     p1,p2,p3,p4,mup1,mup2test,mup3,mup4,musq,-1,rg2,resultZERO)
cfc          call fermionBoxEvenCouplings(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
cfc     &     p1,p2,p3,p4,mup1,mup2,mup3test,mup4,musq,-1,rg3,resultZERO)
cfc          call fermionBoxEvenCouplings(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
cfc     &     p1,p2,p3,p4,mup1,mup2,mup3,mup4test,musq,-1,rg4,resultZERO)
cfc
cfcc Gauge 1
cfc          if (abs(rg1).gt.1d-7) then
cfc          temp1=abs(rg1/(resultn(2)) -1d0)
cfc          else
cfc          temp1=abs(rg1-resultn(2))   
cfc          endif
cfc
cfcc Gauge 2
cfc         if (abs(rg2).gt.1d-7) then
cfc          temp2=abs(rg2/resultn(3) -1d0) 
cfc          else
cfc          temp2=abs(rg2 -resultn(3)  ) 
cfc          endif
cfc
cfcc Gauge 3
cfc         if (abs(rg3).gt.1d-7) then
cfc          temp3=abs(rg3/resultn(4) -1d0 )
cfc          else
cfc             temp3=abs(rg3  -resultn(4) )
cfc         endif
cfc
cfcc Gauge 4
cfc          if (abs(rg4).gt.1d-7) then
cfc          temp4=abs(rg4/resultn(5) -1d0)
cfc          else
cfc          temp4=abs(rg4 - resultn(5)  )
cfc          endif
cfc
cfc          resultgauge=max(temp1,temp2,temp3,temp4)
cfc
cfcCCCC  CHECK whether if improves precision
cfc
cfc        if(resultgauge.lt. FinalSensitivity
cfc     & .or.resultgauge_cp.lt.FinalSensitivity) then
cfc            if(resultgauge.lt.resultgauge_cp) then 
cfccFc             print*, "BE:to rescued works",Real(m),Real(resultgauge)
cfccFc     &       ,Real(resultgauge_cp),Abs(CMPLX(result_cp)/CMPLX(result)-1)
cfccFC     &        ,CMPLX(result)
cfc
cfcc             print*, "result_new vs result_old", 
cfc           else
cfccFc             print*
cfccFC             Print*, "BE: Original Better",Real(m),Real(resultgauge)
cfccFC     &       ,Real(resultgauge_cp),Abs(CMPLX(result_cp)/CMPLX(result)-1)
cfccFC     &       ,CMPLX(result)
cfc             result=result_cp
cfc             resultgauge=resultgauge_cp
cfc           endif   
cfc          else
cfc            if(resultgauge.gt.resultgauge_cp) then 
cfc              result=result_cp
cfc             resultgauge=resultgauge_cp
cfc           endif   
cfccFC             print*
cfccFC             print*, "BE:to rescued failed",Real(m),Real(resultgauge)
cfccFC     &        ,Real(resultgauge_cp),Abs(CMPLX(result_cp)/CMPLX(result)-1)
cfccFC     &        ,CMPLX(result)
cfc             if(abs(result_cp).eq.0d0) then
cfcc            print*
cfccFC            print*, "BE:failed set to WORKS, amplitude is ZERO", result_cp    
cfc             result=result_cp
cfc             resultgauge=Sensitivity
cfc             endif
cfc             if(debug) then
cfc             printMY=.true.
cfc
cfc            print*,"CXXX here",m
cfc             p1sq=dotrr(p1,p1)
cfc             p2sq=dotrr(p2,p2)
cfc             p3sq=dotrr(p3,p3)
cfcc             p4sq=dotrr(p1,p1)
cfc             p1p2=dotrr(p1,p2)
cfc             p1p3=dotrr(p1,p3)
cfc             p2p3=dotrr(p2,p3)
cfc
cfc             print*,"inv",p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,p1sq+p2sq+p3sq+ 2d0*(p1p2+p1p3+p2p3)
cfc
cfc             print*,"inv2",p1sq,p2sq,p3sq,p1sq+p2sq+p3sq+ 2d0*(p1p2+p1p3+p2p3), 2d0*p1p2+p1sq+p2sq,2d0*p2p3+p2sq+p3sq,m*m
cfc
cfc
cfc          print*, "rg1",rg1
cfc          print*, "rg2",rg2
cfc          print*, "rg3",rg3
cfc          print*, "rg4",rg4
cfc          print*, "result1",resultn(2)
cfc          print*, "result",resultn(3)
cfc          print*, "result3",resultn(4)
cfc          print*, "result4",resultn(5)
cfc
cfc          print*, "temp1",temp1
cfc          print*, "temp2",temp2
cfc          print*, "temp3",temp3
cfc          print*, "temp4",temp4
cfc
cfc             do i=0,3
cfc            write(*,'(A3,I1,A2,E25.16,E25.16)') "p1(",i,")=",p1(i)
cfc            enddo
cfc             do i=0,3
cfc            write(*,'(A3,I1,A2,E25.16,E25.16)') "p2(",i,")=",p2(i)
cfc            enddo
cfc             do i=0,3
cfc            write(*,'(A3,I1,A2,E25.16,E25.16)') "p3(",i,")=",p3(i)
cfc            enddo
cfc             do i=0,3
cfc            write(*,'(A3,I1,A2,E25.16,E25.16)') "p4(",i,")=",p4(i)
cfc             enddo
cfc             PRINT*, "mASS",M
cfc       call  FermionBoxEvenCouplings(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,
cfc     & p2,p3,p4,mup1,mup2,mup3,mup4,musq,comp,result,resultn)
cfcc       pause
cfcc       pause
cfc              endif  ! degug
cfc            endif   ! FinalSensitivity
  
#ifdef WITH_QUAD
c        if(.false.) then
            if(abs(resultgauge).gt.FinalSensitivity) then
c               PRINT*, "fINALsENSITIVITY",FinalSensitivity
c              PAUSE
           GaugeE = .true. 
        resultgauge_Q=1d99


       call  FermionBoxEvenCouplings_QUAD(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,
     & p2,p3,p4,mup1,mup2,mup3,mup4,musq,comp,result_Q,resultn_Q)
       
c Call for Gauge
          call fermionBoxEvenCouplings_QUAD(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1test,mup2,mup3,mup4,musq,-1,rg1_Q,resultZERO)
          call fermionBoxEvenCouplings_QUAD(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1,mup2test,mup3,mup4,musq,-1,rg2_Q,resultZERO)
          call fermionBoxEvenCouplings_QUAD(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1,mup2,mup3test,mup4,musq,-1,rg3_Q,resultZERO)
          call fermionBoxEvenCouplings_QUAD(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1,mup2,mup3,mup4test,musq,-1,rg4_Q,resultZERO)

c Gauge 1
          if (abs(rg1_Q).gt.1d-7) then
          temp1=abs(rg1_Q/(resultn_Q(2)) -1d0)
          else
          temp1=abs(rg1_Q-resultn_Q(2))   
          endif

c Gauge 2
         if (abs(rg2_Q).gt.1d-7) then
          temp2=abs(rg2_Q/resultn_Q(3) -1d0) 
          else
          temp2=abs(rg2_Q -resultn_Q(3)  ) 
          endif

c Gauge 3
         if (abs(rg3_Q).gt.1d-7) then
          temp3=abs(rg3_Q/resultn_Q(4) -1d0 )
          else
             temp3=abs(rg3_Q  -resultn_Q(4) )
         endif

c Gauge 4
          if (abs(rg4_Q).gt.1d-7) then
          temp4=abs(rg4_Q/resultn_Q(5) -1d0)
          else
          temp4=abs(rg4_Q - resultn_Q(5)  )
          endif


       if(temp1+1.eq.temp1) then
          temp1=1d9
       endif   

       if(temp1.ne.temp1) then
          temp1=1d9
       endif

       if(temp2+1.eq.temp2) then
          temp2=1d9
       endif   

       if(temp2.ne.temp2) then
          temp2=1d9
       endif

       if(temp3+1.eq.temp3) then
          temp3=1d9
       endif   

       if(temp3.ne.temp3) then
          temp3=1d9
       endif

       if(temp4+1.eq.temp4) then
          temp4=1d9
       endif   

       if(temp4.ne.temp4) then
          temp4=1d9
       endif


          resultgauge_Q=max(temp1,temp2,temp3,temp4)


        if(resultgauge_Q.lt. FinalSensitivity)  then
              wardfb(2)=.true.

           if(abs(resultgauge_Q).lt.abs(resultgauge)) then
             result=result_Q
             resultgauge=resultgauge_Q

c           print*,"BE:QUAD rescued",Real(m),Real(resultgauge),
c     &  Real(resultgauge_cp),Abs(CMPLX(result_cp)/CMPLX(result)-1),CMPLX(result)

           if(abs(result_Q).gt.1d-6) then
!               print*
!               print*
!               print*
!               print*,"BE:QUAD  rescued faile??"
!               print*
!               print*
!               print*

!        call  FermionBoxEvenCouplings_QUAD_check(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,
!      & p2,p3,p4,mup1,mup2,mup3,mup4,musq,comp,result_Q,resultn_Q)
!           print*,"BE:QUAD  rescued faile??",Real(m),resultgauge,
!      &  resultgauge_Q,Abs(CMPLX(result_Q)/CMPLX(result)-1),CMPLX(result)
!            print*,'p1',p1,dotrr(p1,p1)
!            print*,'p2',p2,dotrr(p2,p2)
!            print*,'p3',p3,dotrr(p3,p3)
!            print*,'p4',p4,dotrr(p4,p4)
!            print*, "result", result
!            print*, "result_Q", result_Q
!            print*, temp1,temp2,temp3,temp4
!            print*, "resultnQ", resultn_Q(1:5)
!            print*, "resultnQ", resultn(1:5)
!            print*
!            pause
           endif
           endif
           else
c           print*,"BE:QUAD to rescued failed",Real(m),resultgauge,
c     &  resultgauge_Q,Abs(CMPLX(result_Q)/CMPLX(result)-1),CMPLX(result)
c           print*, "result", result
c           print*, "result_Q", result_Q
c           print*, temp1,temp2,temp3,temp4
c           print*, "resultnQ", resultn_Q(1:5)
c           print*, "resultnQ", resultn(1:5)
c           print*
cfc           if(abs(result_Q).gt.1d-6) then
cfc       call  FermionBoxEvenCouplings_QUAD_check(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,
cfc     & p2,p3,p4,mup1,mup2,mup3,mup4,musq,comp,result_Q,resultn_Q)
cfc           print*,'p1',p1,dotrr(p1,p1)
cfc           print*,'p2',p2,dotrr(p2,p2)
cfc           print*,'p3',p3,dotrr(p3,p3)
cfc           print*,'p4',p4,dotrr(p4,p4)
c           print*, "result", result
c           print*, "result_Q", result_Q
c           print*, temp1,temp2,temp3,temp4
c           print*, "resultnQ", resultn_Q(1:5)
c           print*, "resultnQ", resultn(1:5)
c           print*

cfc       endif
           

        endif
        endif
#endif
           else  ! sensitivity

	wardfb(1)=.true.
        wardfb(2)=.true.

         endif   ! Sensitivity
        endif ! gauge
         
        end


cccccccccccccccccccc    ODD PART ccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine FermionBoxCouplingTO(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,p2,p3,p4,
     &  mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge,Div)
       implicit none
       Real*8 za1,za2,za3,za4,zb1,zb2,zb3,zb4
       Real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3)

       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3)

       Complex*16 result
       Real*8 musq,m,resultgauge

       integer gauge,comp,Div

       If (Div.eq.0) then
       call FermionBoxCouplingO(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,p2,p3,p4
     & ,mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge)
       else
       call FermionBoxCouplingO_Div(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,p2,
     & p3,p4,mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge,Div)
       endif
       end


       subroutine FermionBoxCouplingO(za1,za2,za3,za4,zb1,zb2,zb3,zb4,
     & m,p1,p2,p3,p4,mup1,mup2,mup3,mup4,musq,gauge,comp,result,resultgauge)
       implicit none
#include "VBFNLO/loops/Utilities/GaugeTest.inc"
       Real*8 za1,za2,za3,za4,zb1,zb2,zb3,zb4
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       Real*8   p12(0:3),p14(0:3),p23(0:3),p34(0:3)
       Complex*16 mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3)
       Complex*16 mup1test(0:3),mup2test(0:3),mup3test(0:3),mup4test(0:3)
       Complex*16 result,resultn(5),resultZERO(5)
       Real*8 musq,m,resultgauge,temp1,temp2,temp3,temp4
       real*8 resultgauge_cp
       complex*16 result_cp
       complex*16 rg1,rg2,rg3,rg4
       integer gauge,i,comp
       real*8 Sensitivity, FinalSensitivity
c For QUAD
       Real*8 resultgauge_q
       complex*16 result_Q
       complex*16 resultn_Q(5)
       complex*16 rg1_Q,rg2_Q,rg3_Q
       complex*16 rg4_Q
       logical wardfb(2)
       common/wardpfermionbox/wardfb
       logical GaugeE,GaugeO
       common/FloopPenggzz/GaugeE,GaugeO
       logical DEBUGodd
       PARAMETER(DEBUGodd=.false.)

       iF(DEBUGodd) THEN
          RESULT=0D0
          RESULTGAUGE=0D0
          return
       ENDIF


	wardfb(1)=.false.
        wardfb(2)=.false.

        Sensitivity=GaugeAccuracyBoxL
        FinalSensitivity=GaugeAccuracyBoxL

c        Sensitivity=1d-9
c        FinalSensitivity=1d-9

        resultgauge=1d99

       call  FermionBoxOddCouplings(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,
     & p2,p3,p4,mup1,mup2,mup3,mup4,musq,comp,result,resultn)

       if (gauge.gt.0) then
       
         do  i=0,3
            mup1test(i)=p1(i)
            mup2test(i)=p2(i)
            mup3test(i)=p3(i)
            mup4test(i)=p4(i)
            p12(i)= p1(i)+p2(i)
            p14(i)= p1(i)+p4(i)
            p23(i)= p2(i)+p3(i)
            p34(i)= p3(i)+p4(i)
         enddo
       
c Call for Gauge

          call fermionBoxOddCouplings(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1test,mup2,mup3,mup4,musq,-1,rg1,resultZERO)
          call fermionBoxOddCouplings(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1,mup2test,mup3,mup4,musq,-1,rg2,resultZERO)
          call fermionBoxOddCouplings(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1,mup2,mup3test,mup4,musq,-1,rg3,resultZERO)
          call fermionBoxOddCouplings(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1,mup2,mup3,mup4test,musq,-1,rg4,resultZERO)
c Gauge 1
          if (abs(rg1).gt.1d-7) then
          temp1=abs(rg1/(resultn(2)) -1d0)
          else
          temp1=abs(rg1-resultn(2))   
          endif

c Gauge 2
         if (abs(rg2).gt.1d-7) then
          temp2=abs(rg2/resultn(3) -1d0) 
          else
          temp2=abs(rg2 -resultn(3)  ) 
          endif

c Gauge 3
         if (abs(rg3).gt.1d-7) then
          temp3=abs(rg3/resultn(4) -1d0 )
          else
             temp3=abs(rg3  -resultn(4) )
         endif

c Gauge 4
          if (abs(rg4).gt.1d-7) then
          temp4=abs(rg4/resultn(5) -1d0)
          else
          temp4=abs(rg4 - resultn(5)  )
          endif
       if(temp1+1.eq.temp1) then
          temp1=1d9
       endif   

       if(temp1.ne.temp1) then
          temp1=1d9
       endif

       if(temp2+1.eq.temp2) then
          temp2=1d9
       endif   

       if(temp2.ne.temp2) then
          temp2=1d9
       endif

       if(temp3+1.eq.temp3) then
          temp3=1d9
       endif   

       if(temp3.ne.temp3) then
          temp3=1d9
       endif

       if(temp4+1.eq.temp4) then
          temp4=1d9
       endif   
      
c FC %          print*, "rg1O",rg1
c FC %          print*, "rg2O",rg2
c FC %          print*, "rg3O",rg3
c FC %          print*, "rg4O",rg4
c FC %          print*, "temp1",temp1
c FC %          print*, "temp2",temp2
c FC %          print*, "temp3",temp3
c FC %          print*, "temp4",temp4
c FC %c	read(*,*)

          resultgauge=max(temp1,temp2,temp3,temp4)


c Safe results
          result_cp=result
          resultgauge_cp=resultgauge

cfccCCCCCCCCCCCc TO RESCUE
cfc
cfc        if(.false.) then
        if(resultgauge.gt.Sensitivity) then


#ifdef WITH_QUAD
c        if(.false.) then
            if(abs(resultgauge).gt.FinalSensitivity) then
c               PRINT*, "fINALsENSITIVITY",FinalSensitivity
c              PAUSE

           GaugeE = .true. 
        resultgauge_Q=1d99

        call  FermionBoxOddCouplings_QUAD(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,p1,
     & p2,p3,p4,mup1,mup2,mup3,mup4,musq,comp,result_Q,resultn_Q)


c Call for Gauge

          call fermionBoxOddCouplings_QUAD(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1test,mup2,mup3,mup4,musq,-1,rg1_Q,resultZERO)
          call fermionBoxOddCouplings_QUAD(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1,mup2test,mup3,mup4,musq,-1,rg2_Q,resultZERO)
          call fermionBoxOddCouplings_QUAD(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1,mup2,mup3test,mup4,musq,-1,rg3_Q,resultZERO)
          call fermionBoxOddCouplings_QUAD(za1,za2,za3,za4,zb1,zb2,zb3,zb4,m,
     &     p1,p2,p3,p4,mup1,mup2,mup3,mup4test,musq,-1,rg4_Q,resultZERO)
c Gauge 1
          if (abs(rg1_Q).gt.1d-7) then
          temp1=abs(rg1_Q/(resultn_Q(2)) -1d0)
          else
          temp1=abs(rg1_Q-resultn_Q(2))   
          endif

c Gauge 2
         if (abs(rg2_Q).gt.1d-7) then
          temp2=abs(rg2_Q/resultn_Q(3) -1d0) 
          else
          temp2=abs(rg2_Q -resultn_Q(3)  ) 
          endif

c Gauge 3
         if (abs(rg3_Q).gt.1d-7) then
          temp3=abs(rg3_Q/resultn_Q(4) -1d0 )
          else
             temp3=abs(rg3_Q  -resultn_Q(4) )
         endif

c Gauge 4
          if (abs(rg4_Q).gt.1d-7) then
          temp4=abs(rg4_Q/resultn_Q(5) -1d0)
          else
          temp4=abs(rg4_Q - resultn_Q(5)  )
          endif


       if(temp1+1.eq.temp1) then
          temp1=1d9
       endif   

       if(temp1.ne.temp1) then
          temp1=1d9
       endif

       if(temp2+1.eq.temp2) then
          temp2=1d9
       endif   

       if(temp2.ne.temp2) then
          temp2=1d9
       endif

       if(temp3+1.eq.temp3) then
          temp3=1d9
       endif   

       if(temp3.ne.temp3) then
          temp3=1d9
       endif

       if(temp4+1.eq.temp4) then
          temp4=1d9
       endif   

       if(temp4.ne.temp4) then
          temp4=1d9
       endif

          resultgauge_Q=max(temp1,temp2,temp3,temp4)


        if(resultgauge_Q.lt. FinalSensitivity) then
c             print*, "BO:QUAD to rescued works",Real(m),Real(resultgauge),
c     &        Real(resultgauge_Q),Abs(CMPLX(result)/CMPLX(result_Q)-1)
c     &        ,CMPLX(result_Q)
             result=result_Q
             resultgauge=resultgauge_Q
           if(abs(resultgauge_Q).lt.abs(resultgauge)) then
             result=result_Q
             resultgauge=resultgauge_Q
           endif
        else
c           print*
c           print*,"BO:QUAD to rescued failed",Real(m),Real(resultgauge_cp),
c     &  Real(resultgauge_Q),Abs(CMPLX(result_cp)/CMPLX(result_Q)-1),CMPLX(result_cp)
           endif
           endif
#endif
           ELSE
	wardfb(1)=.true.
        wardfb(2)=.true.
           endif
        
        endif !if gauge
         
        end
