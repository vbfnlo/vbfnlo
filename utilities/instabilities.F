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

      integer ps_number

      boxcountps(ps_number)=0
      boxcount2ps(ps_number)=0
      box2countps(ps_number)=0
      box2count2ps(ps_number)=0
      pentcountps(ps_number)=0
      pentcount2ps(ps_number)=0
      hexcountps(ps_number)=0
      hexcount2ps(ps_number)=0
      
 
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

      integer ps_number
      integer gnlo
      common /nlovariable/ gnlo

      Select Case(procID)
      CASE(WPZA,WMZA,WPAA,WMAA,ZAA,AAA,WPAAj,WMAAj,WPAJJ,WMAJJ,ZAJJ,ZAjj_n,
     &     WPhadZA, WPZhadA, WMhadZA, WMZhadA)

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

      integer i
      
      Select Case(procID)
      CASE(WPZA,WMZA,WPAA,WMAA,ZAA,AAA,WPAAj,WMAAj,WPAJJ,WMAJJ,ZAJJ,ZAjj_n,
     &     WPhadZA, WPZhadA, WMhadZA, WMZhadA)

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
     &    ,QCDWPAjj, QCDWMAjj,QCDZAjj_l,QCDZAjj_n,QCDAAjj)

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
           instab_error = instab_error+ 
     &                        abs(final_xsec(i,0,2)+final_xsec(i,0,4))
     &                        *dble(box2countsub_Q(i)+boxcountsub_Q(i)+pentcountsub_Q(i)+hexcountsub_Q(i))
     &                                       /boxcount2sub(i)
         enddo

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

      write(*,*) ""

      if (PS_Loops .gt. 1) then
         write(*,*) ""
         write(*,*) "LO for different PHASESPACE:"
         do ps_number =1,PS_Loops
            write(*,*)  final_xsec(0,ps_number,-1)," +- " ,sqrt(final_sdev2(0,ps_number,-1)) ," fb", 
     &           sqrt(final_sdev2(0,ps_number,-1))/final_xsec(0,ps_number,-1)*100d0
         enddo
      endif


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
c temporary remove RE contribution form results
      final_xsec(:,:,0) = final_xsec(:,:,0) - final_xsec(:,:,NLO_loops+1)
      final_sdev2(:,:,0) = final_sdev2(:,:,0) - final_sdev2(:,:,NLO_loops+1)
      final_time(:,:,0) = final_xsec(:,:,0) - final_time(:,:,NLO_loops+1)
      
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
      
      
