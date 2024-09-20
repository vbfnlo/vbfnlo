c Francisco Campanario
c Date 2012 05 04
c francam@particle.uni-kalrsruhe.de

      subroutine Init_Instabilities(ps_number)
      implicit none
c include global variables via include-files, e.g. "procID"

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/tau_hel.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/amplitudes/ggf/ggf_process.inc"
      integer ps_number
c     fc Temp for debugging GFHjj
      integer count,countBAD
      Common/ForpentCheck/count,countBAD
      count=0
      countBad=0
cfc Temp for debugging GFHJJ
      

    
      boxcountps(ps_number)=0
      boxcount2ps(ps_number)=0
      box2countps(ps_number)=0
      box2count2ps(ps_number)=0
      pentcountps(ps_number)=0
      pentcount2ps(ps_number)=0
      hexcountps(ps_number)=0
      hexcount2ps(ps_number)=0
                hexinst=0
                peninst=0
                boxinst=0
                verinst=0
      
 
      end


********************************************************************************
********************************************************************************

      subroutine Instabilities(ps_number)

      implicit none

c include global variables via include-files, e.g. "procID"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/tau_hel.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/amplitudes/ggf/ggf_process.inc"
      integer ps_number
      integer gnlo
      common /nlovariable/ gnlo

      Select Case(procID)
      CASE(WPZA,WMZA,WPAA,WMAA,ZAA,AAA,WPAAj,WMAAj,WPAJJ,WMAJJ,ZAJJ,ZAjj_n,
     &     WPhadZA, WPZhadA, WMhadZA, WMZhadA, EWAAjj)

         If(gnlo.eq.2) then

            print *," Box instabilites", boxcountps(ps_number)
            print *,"          total calls", 
     &           boxcountps(ps_number)+boxcount2ps(ps_number)
            if (Dble(boxcountps(ps_number)+boxcount2ps(ps_number)).
     &           ne.0d0) then
               print *,"          instabilities [%]", 
     &              100d0*boxcountps(ps_number)
     &              /Dble(boxcountps(ps_number)+boxcount2ps(ps_number))
            else
               print *,"          instabilities [%]", 0d0
            endif

         elseIf(gnlo.eq.3) then
            print *," Pentagon instabilites", pentcountps(ps_number)
            print *,"          total calls", 
     &           pentcountps(ps_number)+pentcount2ps(ps_number)
            if (Dble(pentcountps(ps_number)+
     &           pentcount2ps(ps_number)).ne.0d0) then
               print *,"          instabilities [%]", 
     &              100d0*pentcountps(ps_number)
     &              /Dble(pentcountps(ps_number)+pentcount2ps(ps_number))
            else
               print *,"          instabilities [%]", 0d0
            endif
            
         elseIf(gnlo.eq.4) then
            print *," Hexagon instabilites", Hexcountps(ps_number)
            print *,"          total calls", 
     &           hexcountps(ps_number)+hexcount2ps(ps_number)
            if (Dble(hexcountps(ps_number)+
     &           hexcount2ps(ps_number)).ne.0d0) then
               print *,"          instabilities [%]", 
     &              100d0*hexcountps(ps_number)
     &              /Dble(hexcountps(ps_number)+hexcount2ps(ps_number))
            else
               print *,"          instabilities [%]", 0d0
            endif
            
         elseIf(gnlo.eq.5) then
                        
            print *," Box2 instabilites", box2countps(ps_number)
            print *,"          total calls", 
     &           box2countps(ps_number)+box2count2ps(ps_number)
            if (Dble(box2countps(ps_number)+
     &           box2count2ps(ps_number)).ne.0d0) then
               print *,"          instabilities [%]", 
     &              100d0*box2countps(ps_number)
     &              /Dble(box2countps(ps_number)+box2count2ps(ps_number))
            else
               print *,"          instabilities [%]", 0d0
            endif
            
         endif

      CASE (GFHjjj) 

         if(Q_loop.ne.0) then
            print*
            print*
           print *," Vertex instabilites:Dble", Verinst(gf_iteration,1,1)
            print *,"          total calls", 
     &            Verinst(gf_iteration,1,1) + Verinst(gf_iteration,2,1)
            if (Dble(Verinst(gf_iteration,1,1) + Verinst(gf_iteration,2,1))
     &                 .ne.0d0) then
               print *,"          instabilities [%]", 
     &              100d0*Verinst(gf_iteration,1,1)
     &              /Dble(Verinst(gf_iteration,1,1) + Verinst(gf_iteration,2,1))
            else
               print *,"          instabilities [%]", 0d0
            endif
#ifdef WITH_QUAD            
            print *," Vertex instabilites:QUAD", Verinst(gf_iteration,1,2)
            print *,"          total calls", 
     &            Verinst(gf_iteration,1,2) + Verinst(gf_iteration,2,2)
            if (Dble(Verinst(gf_iteration,1,2) + Verinst(gf_iteration,2,2))
     &                 .ne.0d0) then
               print *,"          instabilities [%]", 
     &              100d0*Verinst(gf_iteration,1,2)
     &              /Dble(Verinst(gf_iteration,1,2) + Verinst(gf_iteration,2,2))
            else
               print *,"          instabilities [%]", 0d0
            endif
#endif            
          if(boxeva) then
            print*
            print*
           print *," Boxes instabilites:Dble", Boxinst(gf_iteration,1,1)
            print *,"          total calls", 
     &            Boxinst(gf_iteration,1,1) + Boxinst(gf_iteration,2,1)
            if (Dble(Boxinst(gf_iteration,1,1) + Boxinst(gf_iteration,2,1))
     &                 .ne.0d0) then
               print *,"          instabilities [%]", 
     &              100d0*Boxinst(gf_iteration,1,1)
     &              /Dble(Boxinst(gf_iteration,1,1) + Boxinst(gf_iteration,2,1))
            else
               print *,"          instabilities [%]", 0d0
            endif
#ifdef WITH_QUAD           
            print *," Boxes instabilites:QUAD", Boxinst(gf_iteration,1,2)
            print *,"          total calls", 
     &            Boxinst(gf_iteration,1,2) + Boxinst(gf_iteration,2,2)
            if (Dble(Boxinst(gf_iteration,1,2) + Boxinst(gf_iteration,2,2))
     &                 .ne.0d0) then
               print *,"          instabilities [%]", 
     &              100d0*Boxinst(gf_iteration,1,2)
     &              /Dble(Boxinst(gf_iteration,1,2) + Boxinst(gf_iteration,2,2))
            else
               print *,"          instabilities [%]", 0d0
            endif
#endif            
          endif  
          endif
          if (subprQG.or.subprGG) then
             if(peneva) then
                print*
                print*
                print *," Pentagon instabilites:Dble", Peninst(gf_iteration,1,1)
                print *,"          total calls", 
     &               Peninst(gf_iteration,1,1) + Peninst(gf_iteration,2,1)
                if (Dble(Peninst(gf_iteration,1,1) + Peninst(gf_iteration,2,1))
     &               .ne.0d0) then
                   print *,"          instabilities [%]", 
     &                  100d0*Peninst(gf_iteration,1,1)
     &                  /Dble(Peninst(gf_iteration,1,1) + Peninst(gf_iteration,2,1))
                else
               print *,"          instabilities [%]", 0d0
            endif
#ifdef WITH_QUAD            
            print *," Pentagon instabilites:QUAD", Peninst(gf_iteration,1,2)
            print *,"          total calls", 
     &           Peninst(gf_iteration,1,2) + Peninst(gf_iteration,2,2)
            if (Dble(Peninst(gf_iteration,1,2) + Peninst(gf_iteration,2,2))
     &           .ne.0d0) then
               print *,"          instabilities [%]", 
     &              100d0*Peninst(gf_iteration,1,2)
     &              /Dble(Peninst(gf_iteration,1,2) + Peninst(gf_iteration,2,2))
            else
               print *,"          instabilities [%]", 0d0
            endif
#endif            
            print*
            print*
         endif
      endif
            if (subprGG) then
               if(hexeva) then
                  print *," Hexagon instabilites:Dble", Hexinst(gf_iteration,1,1)
               print *,"          total calls", 
     &                 Hexinst(gf_iteration,1,1) + Hexinst(gf_iteration,2,1)
               if (Dble(Hexinst(gf_iteration,1,1) + Hexinst(gf_iteration,2,1))
     &              .ne.0d0) then
                  print *,"          instabilities [%]", 
     &                 100d0*Hexinst(gf_iteration,1,1)
     &                 /Dble(Hexinst(gf_iteration,1,1) + Hexinst(gf_iteration,2,1))
               else
                  print *,"          instabilities [%]", 0d0
               endif
#ifdef WITH_QUAD                        
               print *," Hexagon instabilites:QUAD", Hexinst(gf_iteration,1,2)
               print *,"          total calls", 
     &              Hexinst(gf_iteration,1,2) + Hexinst(gf_iteration,2,2)
               if (Dble(Hexinst(gf_iteration,1,2) + Hexinst(gf_iteration,2,2))
     &              .ne.0d0) then
                  print *,"          instabilities [%]", 
     &                 100d0*Hexinst(gf_iteration,1,2)
     &                 /Dble(Hexinst(gf_iteration,1,2) + Hexinst(gf_iteration,2,2))
               else
                  print *,"          instabilities [%]", 0d0
               endif
#endif               
            endif
         endif

      CASE (GFHjj) 

         if(Q_loop.ne.0) then
            print*
            print*
            print *," Vertex instabilites:Dble", Verinst(gf_iteration,1,1)
            print *,"          total calls", 
     &           Verinst(gf_iteration,1,1) + Verinst(gf_iteration,2,1)
            if (Dble(Verinst(gf_iteration,1,1) + Verinst(gf_iteration,2,1))
     &           .ne.0d0) then
               print *,"          instabilities [%]", 
     &              100d0*Verinst(gf_iteration,1,1)
     &              /Dble(Verinst(gf_iteration,1,1) + Verinst(gf_iteration,2,1))
            else
               print *,"          instabilities [%]", 0d0
            endif
#ifdef WITH_QUAD            
            print *," Vertex instabilites:QUAD", Verinst(gf_iteration,1,2)
            print *,"          total calls", 
     &           Verinst(gf_iteration,1,2) + Verinst(gf_iteration,2,2)
            if (Dble(Verinst(gf_iteration,1,2) + Verinst(gf_iteration,2,2))
     &           .ne.0d0) then
               print *,"          instabilities [%]", 
     &              100d0*Verinst(gf_iteration,1,2)
     &              /Dble(Verinst(gf_iteration,1,2) + Verinst(gf_iteration,2,2))
            else
               print *,"          instabilities [%]", 0d0
            endif
#endif            
            if (subprQG.or.subprGG) then
            
               if(boxeva) then
                  print*
                  print*
                  print *," Boxes instabilites:Dble", Boxinst(gf_iteration,1,1)
                  print *,"          total calls", 
     &                 Boxinst(gf_iteration,1,1) + Boxinst(gf_iteration,2,1)
                  if (Dble(Boxinst(gf_iteration,1,1) + Boxinst(gf_iteration,2,1))
     &                 .ne.0d0) then
                     print *,"          instabilities [%]", 
     &                    100d0*Boxinst(gf_iteration,1,1)
     &                    /Dble(Boxinst(gf_iteration,1,1) + Boxinst(gf_iteration,2,1))
                  else
                     print *,"          instabilities [%]", 0d0
                  endif
#ifdef WITH_QUAD                              
                  print *," Boxes instabilites:QUAD", Boxinst(gf_iteration,1,2)
                  print *,"          total calls", 
     &                 Boxinst(gf_iteration,1,2) + Boxinst(gf_iteration,2,2)
                  if (Dble(Boxinst(gf_iteration,1,2) + Boxinst(gf_iteration,2,2))
     &                 .ne.0d0) then
                     print *,"          instabilities [%]", 
     &                    100d0*Boxinst(gf_iteration,1,2)
     &                    /Dble(Boxinst(gf_iteration,1,2) + Boxinst(gf_iteration,2,2))
                  else
                     print *,"          instabilities [%]", 0d0
                  endif
#endif                  
               endif  
            endif


            if (subprGG) then
               if(peneva) then
                  print*
                  print*
                  print *," Pentagon instabilites:Dble", Peninst(gf_iteration,1,1)
                  print *,"          total calls", 
     &                 Peninst(gf_iteration,1,1) + Peninst(gf_iteration,2,1)
                  if (Dble(Peninst(gf_iteration,1,1) + Peninst(gf_iteration,2,1))
     &                 .ne.0d0) then
                     print *,"          instabilities [%]", 
     &                    100d0*Peninst(gf_iteration,1,1)
     &                    /Dble(Peninst(gf_iteration,1,1) + Peninst(gf_iteration,2,1))
                  else
                     print *,"          instabilities [%]", 0d0
                  endif
#ifdef WITH_QUAD                              
                  print *," Pentagon instabilites:QUAD", Peninst(gf_iteration,1,2)
                  print *,"          total calls", 
     &                 Peninst(gf_iteration,1,2) + Peninst(gf_iteration,2,2)
                  if (Dble(Peninst(gf_iteration,1,2) + Peninst(gf_iteration,2,2))
     &                 .ne.0d0) then
                     print *,"          instabilities [%]", 
     &                    100d0*Peninst(gf_iteration,1,2)
     &                    /Dble(Peninst(gf_iteration,1,2) + Peninst(gf_iteration,2,2))
                  else
                     print *,"          instabilities [%]", 0d0
                  endif
                  print*
                  print*
#endif
              endif

            endif
         endif
c     stop
       
      END SELECT


      END


********************************************************************************
********************************************************************************
      
      subroutine Final_Instabilities

      implicit none
c include global variables via include-files, e.g. "procID"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/tau_hel.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/amplitudes/ggf/ggf_process.inc"
      integer i, Total_Calls
      
      
      Select Case(procID)
      CASE(WPZA,WMZA,WPAA,WMAA,ZAA,AAA,WPAAj,WMAAj,WPAJJ,WMAJJ,ZAJJ,ZAjj_n,
     &     WPhadZA, WPZhadA, WMhadZA, WMZhadA, EWAAjj)

         print *," "
         print *," Box instabilites", boxcountm
         print *,"          total calls", 
     &        boxcountm+boxcount2m
         if (Dble(boxcountm+boxcount2m).ne.0d0) then
            print *,"          instabilities [%]", 
     &           100d0*boxcountm
     &           /Dble(boxcountm+boxcount2m)
         else
            print *,"          instabilities [%]", 0d0
         endif
         

         print *," FBox instabilites", box2countm
         print *,"          total calls", 
     &        box2countm+box2count2m
         if (Dble(box2countm+box2count2m).ne.0d0) then
            print *,"          instabilities [%]", 
     &           100d0*box2countm
     &           /Dble(box2countm+box2count2m)
         else
            print *,"          instabilities [%]", 0d0
         endif
         
         
         print *," Pentagon instabilites", pentcountm
         print *,"          total calls", 
     &        pentcountm+pentcount2m
         if (Dble(pentcountm+pentcount2m).ne.0d0) then
            print *,"          instabilities [%]", 
     &           100d0*pentcountm/Dble(pentcountm+pentcount2m)
         else
            print *,"          instabilities [%]", 0d0
         endif
         
         
         print *," Hexagon instabilites", Hexcountm
         print *,"          total calls", 
     &        hexcountm+hexcount2m
         if (Dble(hexcountm+hexcount2m).ne.0d0) then
            print *,"          instabilities [%]", 
     &           100d0*hexcountm
     &           /Dble(hexcountm+hexcount2m)
         else
            print *,"          instabilities [%]", 0d0
         endif

      CASE(GFWWj,GFZZj)

         print *," "
         print *," Vertex instabilites", box2countm
         print *,"          total calls", 
     &        box2countm+box2count2m
         if (Dble(box2countm+box2count2m).ne.0d0) then
            print *,"          instabilities [%]", 
     &           100d0*box2countm
     &           /Dble(box2countm+box2count2m)
         else
            print *,"          instabilities [%]", 0d0
         endif
         
         
         print *," Box instabilites", boxcountm
         print *,"          total calls", 
     &        boxcountm+boxcount2m
         if (Dble(boxcountm+boxcount2m).ne.0d0) then
            print *,"          instabilities [%]", 
     &           100d0*boxcountm
     &           /Dble(boxcountm+boxcount2m)
         else
            print *,"          instabilities [%]", 0d0
         endif
         

         print *," Pentagon instabilites", pentcountm
         print *,"          total calls", 
     &        pentcountm+pentcount2m
         if (Dble(pentcountm+pentcount2m).ne.0d0) then
            print *,"          instabilities [%]", 
     &           100d0*pentcountm/Dble(pentcountm+pentcount2m)
         else
            print *,"          instabilities [%]", 0d0
         endif
         
      CASE(QCDWPZjj, QCDWMZjj,QCDWPWPjj,QCDWMWMjj,QCDZZjj_ll,QCDZZjj_lnu,QCDZjj_l, QCDZjj_nu, QCDWPjj, QCDWMjj
     &    ,QCDWPAjj, QCDWMAjj,QCDZAjj_l,QCDZAjj_n,QCDAAjj,QCDAjj)

         print *
         print *,  " penbox/ fermionloops instabilites "
         print *,  "     total calls           ", box2count2m
         if (Dble(box2count2m).ne.0d0) then
           print *,"     instabilities dble    ", box2countm, 100d0*box2countm/Dble(box2count2m),"%"
           print *,"     instabilities quad    ", box2countm_Q, 100d0*box2countm_Q/Dble(box2count2m),"%"
         else
           print *,"     instabilities dble    ", box2countm, 0d0, "%"
           print *,"     instabilities quad    ", box2countm, 0d0, "%"
         endif

         print *,  " Box instabilites "
         print *,  "     total calls           ", boxcount2m
         if (Dble(boxcount2m).ne.0d0) then
           print *,"     instabilities dble    ", boxcountm, 100d0*boxcountm/Dble(boxcount2m),"%"
           print *,"     instabilities quad    ", boxcountm_Q, 100d0*boxcountm_Q/Dble(boxcount2m),"%"
         else
           print *,"     instabilities dble    ", boxcountm, 0d0, "%"
           print *,"     instabilities quad    ", boxcountm, 0d0, "%"
         endif

         print *,  " Pentagon instabilites "
         print *,  "     total calls           ", pentcount2m
         if (Dble(pentcount2m).ne.0d0) then
           print *,"     instabilities dble    ", pentcountm, 100d0*pentcountm/Dble(pentcount2m),"%"
           print *,"     instabilities quad    ", pentcountm_Q, 100d0*pentcountm_Q/Dble(pentcount2m),"%"
         else
           print *,"     instabilities dble    ", pentcountm, 0d0, "%"
           print *,"     instabilities quad    ", pentcountm, 0d0, "%"
         endif

         print *,  " Hexagon instabilites "
         print *,  "     total calls           ", hexcount2m
         if (Dble(hexcount2m).ne.0d0) then
           print *,"     instabilities dble    ", hexcountm, 100d0*hexcountm/Dble(hexcount2m),"%"
           print *,"     instabilities quad    ", hexcountm_Q, 100d0*hexcountm_Q/Dble(hexcount2m),"%"
         else
           print *,"     instabilities dble    ", hexcountm, 0d0, "%"
           print *,"     instabilities quad    ", hexcountm, 0d0, "%"
         endif


         instab_error = 0d0
         do i=1,Loops_sub_LO
            if (boxcount2sub(i).ne.0) then
           instab_error = instab_error+ 
     &                        abs(final_xsec(i,0,2)+final_xsec(i,0,4))
     &                        *dble(box2countsub_Q(i)+boxcountsub_Q(i)+pentcountsub_Q(i)+hexcountsub_Q(i))
     &              /boxcount2sub(i)
           endif
         enddo




      CASE (GFHjjj) 

         do i=gf_iteration-1,1,-1
         Verinst(gf_iteration,1,1)=Verinst(i,1,1)+ Verinst(gf_iteration,1,1)
         Verinst(gf_iteration,2,1)=Verinst(i,2,1)+ Verinst(gf_iteration,2,1)
         Verinst(gf_iteration,1,2)=Verinst(i,1,2)+ Verinst(gf_iteration,1,2)
         Verinst(gf_iteration,2,2)=Verinst(i,2,2)+ Verinst(gf_iteration,2,2)
         Boxinst(gf_iteration,1,1)=Boxinst(i,1,1)+ Boxinst(gf_iteration,1,1)
         Boxinst(gf_iteration,2,1)=Boxinst(i,2,1)+ Boxinst(gf_iteration,2,1)
         Boxinst(gf_iteration,1,2)=Boxinst(i,1,2)+ Boxinst(gf_iteration,1,2)
         Boxinst(gf_iteration,2,2)=Boxinst(i,2,2)+ Boxinst(gf_iteration,2,2)
         Peninst(gf_iteration,1,1)=Peninst(i,1,1)+ Peninst(gf_iteration,1,1)
         Peninst(gf_iteration,2,1)=Peninst(i,2,1)+ Peninst(gf_iteration,2,1)
         Peninst(gf_iteration,1,2)=Peninst(i,1,2)+ Peninst(gf_iteration,1,2)
         Peninst(gf_iteration,2,2)=Peninst(i,2,2)+ Peninst(gf_iteration,2,2)
         Hexinst(gf_iteration,1,1)=Hexinst(i,1,1)+ Hexinst(gf_iteration,1,1)
         Hexinst(gf_iteration,2,1)=Hexinst(i,2,1)+ Hexinst(gf_iteration,2,1)
         Hexinst(gf_iteration,1,2)=Hexinst(i,1,2)+ Hexinst(gf_iteration,1,2)
         Hexinst(gf_iteration,2,2)=Hexinst(i,2,2)+ Hexinst(gf_iteration,2,2)
         enddo

         if(Q_loop.ne.0) then         
            write(*,'(A30,I2,A12)') "Total Sum instabilities of", gf_iteration, "iterations"
            print*
            print *," Vertex instabilites:Dble", Verinst(gf_iteration,1,1)
            print *,"          total calls", 
     &            Verinst(gf_iteration,1,1) + Verinst(gf_iteration,2,1)
            if (Dble(Verinst(gf_iteration,1,1) + Verinst(gf_iteration,2,1))
     &                 .ne.0d0) then
               print *,"          instabilities [%]", 
     &              100d0*Verinst(gf_iteration,1,1)
     &              /Dble(Verinst(gf_iteration,1,1) + Verinst(gf_iteration,2,1))
            else
               print *,"          instabilities [%]", 0d0
            endif
#ifdef WITH_QUAD                        
            print *," Vertex instabilites:QUAD", Verinst(gf_iteration,1,2)
            print *,"          total calls", 
     &            Verinst(gf_iteration,1,2) + Verinst(gf_iteration,2,2)
            if (Dble(Verinst(gf_iteration,1,2) + Verinst(gf_iteration,2,2))
     &                 .ne.0d0) then
               print *,"          instabilities [%]", 
     &              100d0*Verinst(gf_iteration,1,2)
     &              /Dble(Verinst(gf_iteration,1,2) + Verinst(gf_iteration,2,2))
            else
               print *,"          instabilities [%]", 0d0
            endif
#endif
         if(boxeva) then
            print*
            print*
           print *," Boxes instabilites:Dble", Boxinst(gf_iteration,1,1)
            print *,"          total calls", 
     &            Boxinst(gf_iteration,1,1) + Boxinst(gf_iteration,2,1)
            if (Dble(Boxinst(gf_iteration,1,1) + Boxinst(gf_iteration,2,1))
     &                 .ne.0d0) then
               print *,"          instabilities [%]", 
     &              100d0*Boxinst(gf_iteration,1,1)
     &              /Dble(Boxinst(gf_iteration,1,1) + Boxinst(gf_iteration,2,1))
            else
               print *,"          instabilities [%]", 0d0
            endif
#ifdef WITH_QUAD                        
            print *," Boxes instabilites:QUAD", Boxinst(gf_iteration,1,2)
            print *,"          total calls", 
     &            Boxinst(gf_iteration,1,2) + Boxinst(gf_iteration,2,2)
            if (Dble(Boxinst(gf_iteration,1,2) + Boxinst(gf_iteration,2,2))
     &                 .ne.0d0) then
               print *,"          instabilities [%]", 
     &              100d0*Boxinst(gf_iteration,1,2)
     &              /Dble(Boxinst(gf_iteration,1,2) + Boxinst(gf_iteration,2,2))
            else
               print *,"          instabilities [%]", 0d0
            endif
#endif      
         endif
         endif
          if(subprQG.or.subprGG) then
             if(peneva) then
                print*
                print*
                print *," Pentagon instabilites:Dble", Peninst(gf_iteration,1,1)
                print *,"          total calls", 
     &               Peninst(gf_iteration,1,1) + Peninst(gf_iteration,2,1)
                if (Dble(Peninst(gf_iteration,1,1) + Peninst(gf_iteration,2,1))
     &               .ne.0d0) then
                   print *,"          instabilities [%]", 
     &                  100d0*Peninst(gf_iteration,1,1)
     &                  /Dble(Peninst(gf_iteration,1,1) + Peninst(gf_iteration,2,1))
                else
                   print *,"          instabilities [%]", 0d0
                endif
#ifdef WITH_QUAD                            
                print *," Pentagon instabilites:QUAD", Peninst(gf_iteration,1,2)
                print *,"          total calls", 
     &               Peninst(gf_iteration,1,2) + Peninst(gf_iteration,2,2)
                if (Dble(Peninst(gf_iteration,1,2) + Peninst(gf_iteration,2,2))
     &               .ne.0d0) then
                   print *,"          instabilities [%]", 
     &                  100d0*Peninst(gf_iteration,1,2)
     &                  /Dble(Peninst(gf_iteration,1,2) + Peninst(gf_iteration,2,2))
                else
                   print *,"          instabilities [%]", 0d0
                endif
#endif          
                print*
                print*
             endif
          endif
          if(subprGG) then
             if(hexeva) then
                print *," Hexagon instabilites:Dble", Hexinst(gf_iteration,1,1)
                print *,"          total calls", 
     &               Hexinst(gf_iteration,1,1) + Hexinst(gf_iteration,2,1)
                if (Dble(Hexinst(gf_iteration,1,1) + Hexinst(gf_iteration,2,1))
     &               .ne.0d0) then
                   print *,"          instabilities [%]", 
     &                  100d0*Hexinst(gf_iteration,1,1)
     &                  /Dble(Hexinst(gf_iteration,1,1) + Hexinst(gf_iteration,2,1))
                else
                   print *,"          instabilities [%]", 0d0
                endif
#ifdef WITH_QUAD                            
                print *," Hexagon instabilites:QUAD", Hexinst(gf_iteration,1,2)
                print *,"          total calls", 
     &            Hexinst(gf_iteration,1,2) + Hexinst(gf_iteration,2,2)
                if (Dble(Hexinst(gf_iteration,1,2) + Hexinst(gf_iteration,2,2))
     &               .ne.0d0) then
                   print *,"          instabilities [%]", 
     &                  100d0*Hexinst(gf_iteration,1,2)
     &                  /Dble(Hexinst(gf_iteration,1,2) + Hexinst(gf_iteration,2,2))
                else
                   print *,"          instabilities [%]", 0d0
                endif
#endif          
         endif
      endif



      CASE (GFHjj) 

         do i=gf_iteration-1,1,-1
         Verinst(gf_iteration,1,1)=Verinst(i,1,1)+ Verinst(gf_iteration,1,1)
         Verinst(gf_iteration,2,1)=Verinst(i,2,1)+ Verinst(gf_iteration,2,1)
         Verinst(gf_iteration,1,2)=Verinst(i,1,2)+ Verinst(gf_iteration,1,2)
         Verinst(gf_iteration,2,2)=Verinst(i,2,2)+ Verinst(gf_iteration,2,2)
         Boxinst(gf_iteration,1,1)=Boxinst(i,1,1)+ Boxinst(gf_iteration,1,1)
         Boxinst(gf_iteration,2,1)=Boxinst(i,2,1)+ Boxinst(gf_iteration,2,1)
         Boxinst(gf_iteration,1,2)=Boxinst(i,1,2)+ Boxinst(gf_iteration,1,2)
         Boxinst(gf_iteration,2,2)=Boxinst(i,2,2)+ Boxinst(gf_iteration,2,2)
         Peninst(gf_iteration,1,1)=Peninst(i,1,1)+ Peninst(gf_iteration,1,1)
         Peninst(gf_iteration,2,1)=Peninst(i,2,1)+ Peninst(gf_iteration,2,1)
         Peninst(gf_iteration,1,2)=Peninst(i,1,2)+ Peninst(gf_iteration,1,2)
         Peninst(gf_iteration,2,2)=Peninst(i,2,2)+ Peninst(gf_iteration,2,2)
         enddo

         if(Q_loop.ne.0) then         
            write(*,'(A30,I2,A12)') "Total Sum instabilities of", gf_iteration, "iterations"
            print*
            print *," Vertex instabilites:Dble", Verinst(gf_iteration,1,1)
            print *,"          total calls", 
     &            Verinst(gf_iteration,1,1) + Verinst(gf_iteration,2,1)
            if (Dble(Verinst(gf_iteration,1,1) + Verinst(gf_iteration,2,1))
     &                 .ne.0d0) then
               print *,"          instabilities [%]", 
     &              100d0*Verinst(gf_iteration,1,1)
     &              /Dble(Verinst(gf_iteration,1,1) + Verinst(gf_iteration,2,1))
            else
               print *,"          instabilities [%]", 0d0
            endif
#ifdef WITH_QUAD                        
            print *," Vertex instabilites:QUAD", Verinst(gf_iteration,1,2)
            print *,"          total calls", 
     &            Verinst(gf_iteration,1,2) + Verinst(gf_iteration,2,2)
            if (Dble(Verinst(gf_iteration,1,2) + Verinst(gf_iteration,2,2))
     &                 .ne.0d0) then
               print *,"          instabilities [%]", 
     &              100d0*Verinst(gf_iteration,1,2)
     &              /Dble(Verinst(gf_iteration,1,2) + Verinst(gf_iteration,2,2))
            else
               print *,"          instabilities [%]", 0d0
            endif
#endif      
            if(subprQG.or.subprGG) then
               if(boxeva) then
                  print*
                  print*
                  print *," Boxes instabilites:Dble", Boxinst(gf_iteration,1,1)
                  print *,"          total calls", 
     &                 Boxinst(gf_iteration,1,1) + Boxinst(gf_iteration,2,1)
                  if (Dble(Boxinst(gf_iteration,1,1) + Boxinst(gf_iteration,2,1))
     &                 .ne.0d0) then
                     print *,"          instabilities [%]", 
     &                    100d0*Boxinst(gf_iteration,1,1)
     &                    /Dble(Boxinst(gf_iteration,1,1) + Boxinst(gf_iteration,2,1))
                  else
                     print *,"          instabilities [%]", 0d0
                  endif
#ifdef WITH_QUAD                              
                  print *," Boxes instabilites:QUAD", Boxinst(gf_iteration,1,2)
                  print *,"          total calls", 
     &                 Boxinst(gf_iteration,1,2) + Boxinst(gf_iteration,2,2)
                  if (Dble(Boxinst(gf_iteration,1,2) + Boxinst(gf_iteration,2,2))
     &                 .ne.0d0) then
                     print *,"          instabilities [%]", 
     &                    100d0*Boxinst(gf_iteration,1,2)
     &                    /Dble(Boxinst(gf_iteration,1,2) + Boxinst(gf_iteration,2,2))
                  else
                     print *,"          instabilities [%]", 0d0
                  endif
#endif            
               endif
            endif
         endif
         if(subprGG) then
            if(peneva) then
               print*
               print*
               print *," Pentagon instabilites:Dble", Peninst(gf_iteration,1,1)
                print *,"          total calls", 
     &               Peninst(gf_iteration,1,1) + Peninst(gf_iteration,2,1)
                if (Dble(Peninst(gf_iteration,1,1) + Peninst(gf_iteration,2,1))
     &               .ne.0d0) then
                   print *,"          instabilities [%]", 
     &                  100d0*Peninst(gf_iteration,1,1)
     &                  /Dble(Peninst(gf_iteration,1,1) + Peninst(gf_iteration,2,1))
                else
                   print *,"          instabilities [%]", 0d0
                endif
#ifdef WITH_QUAD                            
                print *," Pentagon instabilites:QUAD", Peninst(gf_iteration,1,2)
                print *,"          total calls", 
     &               Peninst(gf_iteration,1,2) + Peninst(gf_iteration,2,2)
                if (Dble(Peninst(gf_iteration,1,2) + Peninst(gf_iteration,2,2))
     &               .ne.0d0) then
                   print *,"          instabilities [%]", 
     &                  100d0*Peninst(gf_iteration,1,2)
     &                  /Dble(Peninst(gf_iteration,1,2) + Peninst(gf_iteration,2,2))
                else
                   print *,"          instabilities [%]", 0d0
                endif
#endif          
                print*
                print*
             endif
          endif
       END SELECT



       SELECT CASE(procID)
      CASE(GFHjjj,GFHjj)
c     FC add to divide only if different to zero. Otherwise one gets NAN.
         instab_error = 0d0

         Total_calls=Hexinst(gf_iteration,1,1) + Hexinst(gf_iteration,2,1)
         if(total_calls.ne.0) then
            instab_error = instab_error 
     &        +      (100d0*Hexinst(gf_iteration,1,1)
     &           /Dble(Total_calls))
         endif
         Total_calls=Peninst(gf_iteration,1,1) + Peninst(gf_iteration,2,1)
         if(total_calls.ne.0) then
            instab_error = instab_error 
     &        +      (100d0*Peninst(gf_iteration,1,1)
     &           /Dble(Total_calls))
         endif

         Total_calls=Boxinst(gf_iteration,1,1) + Boxinst(gf_iteration,2,1)
         if(total_calls.ne.0) then
            instab_error = instab_error 
     &        +      (100d0*Boxinst(gf_iteration,1,1)
     &           /Dble(Total_calls))
         endif

        Total_calls=Verinst(gf_iteration,1,1) + Verinst(gf_iteration,2,1)
         if(total_calls.ne.0) then
            instab_error = instab_error 
     &        +      (100d0*Verinst(gf_iteration,1,1)
     &           /Dble(Total_calls))
         endif

c            print*, final_xsec(0,0,-1)
            print*, "Estimated error due to instabilities only Dble:",
     &       final_xsec(0,0,-1)*instab_error/100d0, "[fb]",instab_error,"[%]"
            
            instab_error = 0d0

            Total_calls=Hexinst(gf_iteration,1,2) + Hexinst(gf_iteration,2,2)
         if(total_calls.ne.0) then
            instab_error = instab_error 
     &        +      (100d0*Hexinst(gf_iteration,1,2)
     &           /Dble(Total_calls))
         endif

         Total_calls=Peninst(gf_iteration,1,2) + Peninst(gf_iteration,2,2)
         if(total_calls.ne.0) then
            instab_error = instab_error 
     &        +      (100d0*Peninst(gf_iteration,1,2)
     &           /Dble(Total_calls))
         endif

         Total_calls=Boxinst(gf_iteration,1,2) + Boxinst(gf_iteration,2,2)
         if(total_calls.ne.0) then
            instab_error = instab_error 
     &        +      (100d0*Boxinst(gf_iteration,1,2)
     &           /Dble(Total_calls))
         endif

         Total_calls=Verinst(gf_iteration,1,2) + Verinst(gf_iteration,2,2)
         if(total_calls.ne.0) then
            instab_error = instab_error 
     &        +      (100d0*Verinst(gf_iteration,1,2)
     &           /Dble(Total_calls))
         endif
#ifdef WITH_QUAD                      
            print*, "Estimated error due to instabilities using QUAD:",
     &       final_xsec(0,0,-1)*instab_error/100d0, "[fb]",instab_error,"[%]"
#endif     
            
            END SELECT

             return

      end


********************************************************************************
********************************************************************************
      
      subroutine Final_Printing

      implicit none

c include global variables via include-files, e.g. "procID"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/tau_hel.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/scales.inc"

      integer ps_number,gnlo

      if (PS_Loops .gt. 1) then
         write(*,*) "LO for different PHASESPACE:"
         do ps_number =1,PS_Loops
            write(*,*)  final_xsec(0,ps_number,-1)," +- " ,sqrt(final_sdev2(0,ps_number,-1)) ," fb", 
     &           sqrt(final_sdev2(0,ps_number,-1))/final_xsec(0,ps_number,-1)*100d0
         enddo
      endif
      write(*,*) ""


      if(procID/1000.ne.3) then ! don't show this for QCD induced processes


      if (PS_Loops .gt. 1 .and. doNLO) then
         if (NLO_loops .ge. 1) then
            write(*,*)" "
            write(*,*)" "
            write(*,*)"VIRTUAL CONTRIBUTIONS for different PHASESPACE:"


            write(*,*) " borntype for different PHASESPACE:"
            do ps_number =1,PS_Loops
               write(*,*)  final_xsec(0,ps_number,1)," +- " ,sqrt(final_sdev2(0,ps_number,1)) ,
     &              " fb", 100d0*sqrt(final_sdev2(0,ps_number,1))/final_xsec(0,ps_number,1)
            enddo
         end if
         
         if (NLO_loops .ge. 2) then
            write(*,*) " boxes for different PHASESPACE:"     
            do ps_number =1,PS_Loops
               write(*,*)  final_xsec(0,ps_number,2)," +- " ,sqrt(final_sdev2(0,ps_number,2)) ,
     &              " fb", 100d0*sqrt(final_sdev2(0,ps_number,2))/final_xsec(0,ps_number,2)
            enddo
         end if
         
         if (NLO_loops .ge. 5) then
            write(*,*) " boxes2 for different PHASESPACE:"     
            do ps_number =1,PS_Loops
               write(*,*)  final_xsec(0,ps_number,5)," +- " ,sqrt(final_sdev2(0,ps_number,5)) ,
     &              " fb", 100d0*sqrt(final_sdev2(0,ps_number,5))/final_xsec(0,ps_number,5)
            enddo
         end if
         
         if (NLO_loops .ge. 3) then
            write(*,*) "pentagons for different PHASESPACE:"
            do ps_number =1,PS_Loops
               write(*,*)  final_xsec(0,ps_number,3)," +- " ,sqrt(final_sdev2(0,ps_number,3)) ,
     &              " fb", 100d0*sqrt(final_sdev2(0,ps_number,3))/final_xsec(0,ps_number,3)
            enddo
         end if
         
         if (NLO_loops .ge. 4) then
            write(*,*) "hexagons for different PHASESPACE:"
            do ps_number =1,PS_Loops
               write(*,*)  final_xsec(0,ps_number,4)," +- " ,sqrt(final_sdev2(0,ps_number,4)) ,
     &              " fb", 100d0*sqrt(final_sdev2(0,ps_number,4))/final_xsec(0,ps_number,4)
            enddo
         end if
      end if


c      write(*,*) "K-Factor",  (ps_xsec(2)+ps_xsec(3))/ps_xsec(1)
c      write(*,*) "K-Factor-virtual",  (ps_xsec(2))/ps_xsec(1)
c      write(*,*) ""
      
** NOW THE TOTAL FOR EACH CONTRIBUTION  
      if (NLO_loops .ge. 1 .and. doNLO) then
         write(*,*)"  "
         write(*,*) "TOTAL VIRTUAL CONTRIBUTIONS:"
         write(*,*) "TOTAL borntype:"
         write(*,*)  final_xsec(0,0,1),sqrt(final_sdev2(0,0,1)),
     &        100d0*Sqrt(final_sdev2(0,0,1))/final_xsec(0,0,1)," %"
      end if
         
      if (NLO_loops .ge. 2 .and. doNLO) then
         write(*,*) "TOTAL boxes:"
         write(*,*)  final_xsec(0,0,2),sqrt(final_sdev2(0,0,2)),
     &        100d0*Sqrt(final_sdev2(0,0,2))/final_xsec(0,0,2)," %"
      end if
         
      if (NLO_loops .ge. 5 .and. doNLO) then
         write(*,*) "TOTAL boxes2:"
         write(*,*)  final_xsec(0,0,5),sqrt(final_sdev2(0,0,5)),
     &        100d0*Sqrt(final_sdev2(0,0,5))/final_xsec(0,0,5)," %"
      end if

      if (NLO_loops .ge. 3 .and. doNLO) then
         SELECT CASE(process)
         CASE(WPHJ:WMHJ_ZZ_lnu)
            write(*,*) "TOTAL toploops:"
         CASE DEFAULT
            write(*,*) "TOTAL pentagons:"
         END SELECT
         write(*,*)  final_xsec(0,0,3),sqrt(final_sdev2(0,0,3)),
     &        100d0*Sqrt(final_sdev2(0,0,3))/final_xsec(0,0,3)," %"
      end if
      
      if (NLO_loops .ge. 4 .and. doNLO) then
         write(*,*) "TOTAL hexagons:"
         write(*,*)  final_xsec(0,0,4),sqrt(final_sdev2(0,0,4)),
     &        100d0*Sqrt(final_sdev2(0,0,4))/final_xsec(0,0,4)," %"
      end if


* Real emission:      
      if (PS_loops .gt. 1 .and. doNLO) then
         write(*,*) " "
         write(*,*) " "
         write(*,*) " REAL EMISSION for different PHASESPACE:"
         do ps_number =1,PS_Loops
         write(*,*)  final_xsec(0,ps_number,NLO_loops+1),sqrt(final_sdev2(0,ps_number,NLO_loops+1)),
     &        100d0*Sqrt(final_sdev2(0,ps_number,NLO_loops+1))/final_xsec(0,ps_number,NLO_loops+1)," %"
         enddo
      end if


      else ! QCD processes


      if(doNLO) then
      print*
      print*,'========================================'
      print*,'   results of individual contributions'
      print*,'========================================'
c temporary remove RE contribution from results
      final_xsec(:,:,0) = final_xsec(:,:,0) - final_xsec(:,:,NLO_loops+1)
      final_sdev2(:,:,0) = final_sdev2(:,:,0) - final_sdev2(:,:,NLO_loops+1)
      final_time(:,:,0) = final_time(:,:,0) - final_time(:,:,NLO_loops+1)
      
      print*
      print*,'virtual contributions'

      do gnlo=1,NLO_loops
        print*,'gnlo',gnlo,final_xsec(0,0,gnlo)," +- ",sqrt(final_sdev2(0,0,gnlo))

      enddo
      print*
      print*,'total virtual contributions for different phase spaces'
      do ps_number=1,PS_Loops
        print*,'subproc',ps_number,final_xsec(0,ps_number,0)," +- ",sqrt(final_sdev2(0,ps_number,0))

      enddo
      print*
      print*,'total virtual contributions for different subprocesses'
      do sub_number=1,Loops_sub_LO
        print*,'subproc',sub_number,final_xsec(sub_number,0,0)," +- ",sqrt(final_sdev2(sub_number,0,0))

      enddo
      print*
      print*,'RE subprocesses'
      do ps_number = 1, PS_Loops
      print*,'PS ', ps_number
      do sub_number=1,Loops_sub_NLO
      print*,'subproc',sub_number,final_xsec(sub_number,ps_number,NLO_loops+1)," +- ",
     &                 sqrt(final_sdev2(sub_number,ps_number,NLO_loops+1))

      enddo
      enddo

c add back removed RE contribution form results
      final_xsec(:,:,0) = final_xsec(:,:,0) + final_xsec(:,:,NLO_loops+1)
      final_sdev2(:,:,0) = final_sdev2(:,:,0) + final_sdev2(:,:,NLO_loops+1)
      final_time(:,:,0) = final_time(:,:,0) + final_time(:,:,NLO_loops+1)
      endif ! doNLO


      endif ! QCD process


      return
      end

      logical function IsNaN(dbl,cplx,Ndbl,Ncplx)
      implicit none
      double precision dbl(*)
      double complex cplx(*)
      integer Ndbl,Ncplx,i
      
      IsNaN = .false.
      
      do i=1,Ndbl
        if(dbl(i).ne.dbl(i)) IsNaN=.true.
        if(dbl(i)+1d0.eq.dbl(i)) IsNaN =  .true.
      enddo
      do i=1,Ncplx
        if(cplx(i).ne.cplx(i)) IsNaN =  .true.
        if(cplx(i)+1d0.eq.cplx(i)) IsNaN =  .true.
      enddo
      end
      
      
