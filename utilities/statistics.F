c Francisco Campanario
c
c subroutine statistics_proc
c subroutine statistics_proc_NLO
c subroutine Instabilities
c subroutine AddToXsec(xsed,sdev2,time,sub_number,ps_number,nlo_contrib)

      subroutine statistics_proc(ps_number)
      implicit none
c include global variables via include-files, e.g. "procID"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/tau_hel.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/scales.inc"

c     variables for the integration loops
      integer iteration
      integer*8 Ncall, two
      integer ps_number
      parameter (two = 2)
      common/Ncall1/Ncall,iteration


      Ncall = N_points(0)-N_iterations1 + iteration 
      
      if(.not.reducesubstat) then
        Ncall = two**Ncall
	return
      endif

      if (.not.lreaddata) then
c            print*,'lreaddata',lreaddata

      Select Case(procID)
      CASE(QCDWPWPjj,QCDWMWMjj)
	if(sub_number.eq.1) then
	  Ncall = two**(Ncall)
	elseif(sub_number.eq.2) then
	  Ncall = two**(Ncall-1)
	else 
	  Ncall = two**(Ncall-3)
	endif
      CASE(QCDWPjj,QCDWMjj)
	Select case(sub_number)
	case (1)
	  Ncall = two**(Ncall)
	case (5)
	  Ncall = two**(Ncall -2)
	case (2)
	  Ncall = two**(Ncall -2)
	case (3)
	  Ncall = two**(Ncall -3)
	case (4)
	  Ncall = two**(Ncall- 4)
	case (8)
	  Ncall = two**(Ncall -5)
	case (6,7)
	  Ncall = two**(Ncall -3)
	end select
      CASE(QCDZjj_l,QCDZjj_nu)
        Select case(sub_number)
        case (1)
          Ncall = two**(Ncall)
        case (2)
          Ncall = two**(Ncall -1)
        case (5,6)
          Ncall = two**(Ncall -2)
        case (3,4)
          Ncall = two**(Ncall -3)
        case (7)
          Ncall = two**(Ncall- 4)
        end select
      CASE(QCDWPZjj,QCDWMZjj)
        Select case(sub_number)
        case (1)
          Ncall = two**(Ncall)
        case (5)
          Ncall = two**(Ncall -2)
        case (2)
          Ncall = two**(Ncall -2)
        case (3)
          Ncall = two**(Ncall -2)
        case (4)
          Ncall = two**(Ncall- 3)
        case (8)
          Ncall = two**(Ncall -4)
        case (6,7)
          Ncall = two**(Ncall -3)
        end select
      CASE(QCDWPAjj,QCDWMAjj)
        if(ps_number.eq.1) Ncall = Ncall -1
	Select case(sub_number)
	case (1)
	  Ncall = two**(Ncall)
	case (5,7)
	  Ncall = two**(Ncall -2)
	case (2)
	  Ncall = two**(Ncall -1)
	case (3)
	  Ncall = two**(Ncall -2)
	case (4)
	  Ncall = two**(Ncall- 3)
	case (6,8)
	  Ncall = two**(Ncall -3)
	end select
      CASE(QCDZAjj_l)
        if(ps_number.eq.1) then
         Ncall = Ncall -1
	Select case(sub_number)
	case (1)
	  Ncall = two**(Ncall)
	case (2)
	  Ncall = two**(Ncall -1)
	case (3)
	  Ncall = two**(Ncall -2)
	case (4)
	  Ncall = two**(Ncall- 3)
	case (5)
	  Ncall = two**(Ncall -2)
	case (6)
	  Ncall = two**(Ncall -2)
	case (7)
	  Ncall = two**(Ncall -4)
	end select
        else
	Select case(sub_number)
	case (1)
	  Ncall = two**(Ncall)
	case (2)
	  Ncall = two**(Ncall -1)
	case (3)
	  Ncall = two**(Ncall -2)
	case (4)
	  Ncall = two**(Ncall- 3)
	case (5)
	  Ncall = two**(Ncall -2)
	case (6,8)
	  Ncall = two**(Ncall -2)
	case (7)
	  Ncall = two**(Ncall -5)
	end select
        endif
      CASE(QCDZAjj_n)
	Select case(sub_number)
	case (1)
	  Ncall = two**(Ncall)
	case (2)
	  Ncall = two**(Ncall -1)
	case (3)
	  Ncall = two**(Ncall -2)
	case (4)
	  Ncall = two**(Ncall- 3)
	case (5)
	  Ncall = two**(Ncall -2)
	case (6,8)
	  Ncall = two**(Ncall -2)
	case (7)
	  Ncall = two**(Ncall -5)
	end select
      CASE(QCDZZjj_ll,QCDZZjj_lnu,QCDAAjj)
        Select case(sub_number)
        case (1,5)
          Ncall = two**(Ncall)
        case (6)
          Ncall = two**(Ncall-1)
        case (2)
          Ncall = two**(Ncall-2)
        case (3,4,7)
          Ncall = two**(Ncall-3)
        end select
      CASE(ZAjj)
        if(ps_number.eq.1) Ncall = Ncall - 2
	Ncall = two**Ncall
      CASE(WPZA, WMZA, WPhadZA, WPZhadA, WMhadZA, WMZhadA)
         if (ps_number.eq.1) then
            Ncall = two**( N_points(0) - N_iterations1 + iteration )
         else
            Ncall = two**( N_points(0) + 1 - N_iterations1 + iteration )
         endif
      CASE(WPAA, WMAA,ZAA)
         if (ps_number.eq.1) then
            Ncall = two**( N_points(0) +1 - N_iterations1 + iteration )
         else
            Ncall = two**( N_points(0)  - N_iterations1 + iteration )
         endif
      CASE(WPAAj,WMAAj)
c         print*,'here',ps_number,sub_number,N_iterations1,iteration,Ncall
         if ((N_iterations1-iteration).eq.0) then
            if (ps_number.eq.1) then
               if (sub_number.gt.4) then
                  Ncall = two**( N_points(0)+1)
               else
                  Ncall = two**( N_points(0))
               endif
            elseif(ps_number.lt.4) then
               if (sub_number.gt.2) then
                  Ncall = two**( N_points(0)-1)
               else
                  Ncall = two**( N_points(0)-2)
               endif
            else
               if (sub_number.gt.2) then
                  Ncall = two**( N_points(0)-2)
               else
                  Ncall = two**( N_points(0)-3)
               endif
            endif
         elseif((N_iterations1-iteration).eq.1) then
            if (ps_number.eq.1) then
               if(sub_number.gt.4) then
                  Ncall = two**( N_points(0))
               else
                  Ncall = two**( N_points(0)-1)
               endif
            elseif(ps_number.lt.4)then
               if(sub_number.gt.2) then
                  Ncall = two**( N_points(0)-2)
               else
                  Ncall = two**( N_points(0)-3)
               endif   
            else
               if(sub_number.gt.2) then
                  Ncall = two**( N_points(0)-3)
               else
                  Ncall = two**( N_points(0)-4)
               endif   
            endif
         elseif((N_iterations1-iteration).gt.1) then
            if (ps_number.eq.1) then
               if(sub_number.gt.4) then
                  Ncall = two**( N_points(0)-1)
               else
                  Ncall = two**( N_points(0)-2)
               endif   
            elseif(ps_number.lt.4)then
               if(sub_number.gt.2) then
                  Ncall = two**( N_points(0)-3)
               else
                  Ncall = two**( N_points(0)-4)
               endif
            else
               if(sub_number.gt.2) then
                  Ncall = two**( N_points(0)-4)
               else
                  Ncall = two**( N_points(0)-5)
               endif
            endif
         endif
      CASE DEFAULT
         Ncall = two**( N_points(0) - N_iterations1 + iteration )
      END SELECT
      endif
c     print*,'here2',ps_number,sub_number,N_iterations1,iteration,Ncall
      return
      end


********************************************************************************
********************************************************************************

      subroutine statistics_proc_NLO(ps_number)

      implicit none

c include global variables via include-files, e.g. "procID"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/tau_hel.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/scales.inc"

c     variables for the integration loops
      integer ps_number
      integer iteration
      integer*8 Ncall, two
      parameter (two = 2)
      common/Ncall1/Ncall,iteration
      integer gnlo
      common /nlovariable/ gnlo
      integer add(7)

      Ncall = N_points(gnlo)!-N_iterations1 + iteration 
      
      if(.not.reducesubstat) then
        Ncall = two**Ncall
	return
      endif

      if (.not.lreaddata) then
c          print*,'here',ps_number,sub_number,Ncall
         Select Case(procID)
      CASE(QCDWPWPjj,QCDWMWMjj)
! 	endif
	if(sub_number.eq.1) then
	  Ncall = two**(Ncall)
	elseif(sub_number.eq.2) then
	  Ncall = two**(Ncall-1)
	else
	  Ncall = two**(Ncall-3)
	endif
      CASE(QCDWPjj,QCDWMjj)
	SELECT CASE(sub_number)
	CASE(1)
	  Ncall = two**(Ncall)
	CASE(5)
	  Ncall = two**(Ncall-2)
	CASE(2)
	  Ncall = two**(Ncall-1)
	CASE(3,6,7)
	  Ncall = two**(Ncall-3)
	CASE(4,8)
	  Ncall = two**(Ncall-4)
	END SELECT
      CASE(QCDWPZjj,QCDWMZjj)
	SELECT CASE(sub_number)
	CASE(1,5)
	  Ncall = two**(Ncall)
	CASE(2,3,6,7)
	  Ncall = two**(Ncall-2)
	CASE(4,8)
	  Ncall = two**(Ncall-3)
	END SELECT
      CASE(QCDWPAjj,QCDWMAjj)
        if(ps_number.eq.1) Ncall = Ncall -1
        Select case(sub_number)
        case (1)
          Ncall = two**(Ncall)
        case (5,7)
          Ncall = two**(Ncall -1)
        case (2)
          Ncall = two**(Ncall -1)
        case (3,6)
          Ncall = two**(Ncall -2)
        case (4)
          Ncall = two**(Ncall- 3)
        case (8)
          Ncall = two**(Ncall -3)
        end select
      CASE(QCDZAjj_l)
        if(ps_number.eq.1) then
        Ncall = Ncall -1
         if(gnlo.eq.1) then
            add(1)=0
            add(2)=0
            add(3)=0
            add(4)=1
            add(5)=-1
            add(6)=0
            add(7)=-2
         elseif(gnlo.eq.2) then
            add(1)=1
            add(2)=1
            add(3)=0
            add(4)=1
            add(5)=1
            add(6)=2
            add(7)=-1
         elseif(gnlo.eq.3) then
            add(1)=0
            add(2)=0
            add(3)=0
            add(4)=1
            add(5)=-1
            add(6)=-1
            add(7)=-2
            else
            add(1)=0
            add(2)=0
            add(3)=0
            add(4)=0
            add(5)=0
            add(6)=0
            add(7)=0
         endif
         else
         if(gnlo.eq.1) then
            add(1)=0
            add(2)=-1
            add(3)=0
            add(4)=0
            add(5)=-1
            add(6)=0
            add(7)=-2
         elseif(gnlo.eq.2) then
            add(1)=0
            add(2)=0
            add(3)=0
            add(4)=-1
            add(5)=0
            add(6)=1
            add(7)=-1
         elseif(gnlo.eq.3) then
            add(1)=0
            add(2)=0
            add(3)=-1
            add(4)=0
            add(5)=-1
            add(6)=-1
            add(7)=-2
            else
            add(1)=0
            add(2)=0
            add(3)=0
            add(4)=0
            add(5)=0
            add(6)=0
            add(7)=0
         endif
       endif
        Select case(sub_number)
        case (1)
          Ncall = two**(Ncall+add(1))
        case (2)
          Ncall = two**(Ncall -1+add(2))
        case (3)
          Ncall = two**(Ncall -2+add(3))
        case (4)
          Ncall = two**(Ncall- 3+add(4))
        case (5)
          Ncall = two**(Ncall -1+add(5))
        case (6)
          Ncall = two**(Ncall -2+add(6))
        case (7)
          Ncall = two**(Ncall -1+add(7))
        end select
      CASE(QCDZAjj_n)
        Select case(sub_number)
        case (1)
          Ncall = two**(Ncall)
        case (2)
          Ncall = two**(Ncall -1)
        case (3,6)
          Ncall = two**(Ncall -2)
        case (4)
          Ncall = two**(Ncall- 3)
        case (5)
          Ncall = two**(Ncall -1)
        case (7)
          Ncall = two**(Ncall -3)
        case (8)
          Ncall = two**(Ncall -3)
        end select
      CASE(QCDZZjj_ll,QCDZZjj_lnu,QCDAAjj)
        if(gnlo.ge.4.and.(sub_number.eq.3.or.sub_number.eq.4)) Ncall=Ncall+1 
        SELECT CASE(sub_number)
        CASE(1,5,6)
          Ncall = two**(Ncall)
        CASE(2)
          Ncall = two**(Ncall-1)
        CASE(3)
          Ncall = two**(Ncall-2)
        CASE(4,7)
          Ncall = two**(Ncall-3)
        END SELECT
      CASE(ZAjj)
        if(ps_number.eq.1) Ncall = Ncall - 2
	Ncall = two**Ncall
      CASE(WPZA, WMZA, WPhadZA, WPZhadA, WMhadZA, WMZhadA)
         if (ps_number.eq.1) then
            Ncall = two**(N_points(gnlo))
         else
            Ncall = two**(N_points(gnlo)+1)
         endif
      CASE(WPAA,WMAA,ZAA,WPAAj,WMAAj)
         if (ps_number.eq.1) then
            if((sub_number.gt.4).or.(sub_number.lt.3))then
               if (sub_number.gt.4) then
                  N_points(3)=N_points(3)+1   
                  Ncall = two**(N_points(gnlo)+1)
                  N_points(3)=N_points(3)-1   
               else
                  Ncall = two**(N_points(gnlo)+1)
               endif
            else
               Ncall = two**(N_points(gnlo))
            endif
         else if (ps_number.eq.4) then
            N_points(1)=N_points(1)-1
            if(sub_number.gt.2)then
               Ncall = two**(N_points(gnlo)-2)
            else
               Ncall = two**(N_points(gnlo)-3)
            endif
            N_points(1)=N_points(1)+1
         else
            if(sub_number.gt.2)then
               Ncall = two**(N_points(gnlo)-1)
            else
               N_points(3)=N_points(3)+1   
               Ncall = two**(N_points(gnlo)-2)
               N_points(3)=N_points(3)-1   
            endif
         endif
      CASE DEFAULT
         Ncall = two**N_points(gnlo)
      END SELECT
      endif
c     print*,'here1',ps_number,sub_number,Ncall
      return
      end
      

********************************************************************************
********************************************************************************

      subroutine statistics_proc_Real(ps_number)

      implicit none

c include global variables via include-files, e.g. "procID"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/tau_hel.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/scales.inc"

      integer ps_number
      integer iteration
      integer*8 Ncall, two
      parameter (two = 2)
      common/Ncall1/Ncall,iteration

      integer gnlo
      common /nlovariable/ gnlo

      Ncall = ( N_points(-1) - N_iterations2 + iteration )

      if(.not.reducesubstat) then
        Ncall = two**Ncall
	return
      endif

      if (.not.lreaddata) then
         Select Case(procID)
      CASE(QCDWPAjj,QCDWMAjj,QCDZAjj_l)
	select case(ps_number)
	case(2)
	  Ncall = two**Ncall
	case(1)
	  Ncall = two**(Ncall-1)
	case default
	  Ncall = two**(Ncall)
	end select
      CASE(ZAjj)
        if(ps_number.eq.1) Ncall = Ncall - 2
	Ncall = two**Ncall
      CASE(QCDZAjj_n)
	  Ncall = two**(Ncall)
      CASE(WPZA, WMZA, WPhadZA, WPZhadA, WMhadZA, WMZhadA)
         if (ps_number.eq.1) then
            Ncall = two**( N_points(-1) - N_iterations2 + iteration )
         else
            Ncall = two**( N_points(-1)+1 - N_iterations2 + iteration )
         endif
      CASE(WPAA, WMAA,ZAA)
         if (ps_number.eq.1) then
            Ncall = two**( N_points(-1)+ 1 - N_iterations2 + iteration )
            Ncall = two**( N_points(-1)+1) ! - N_iterations2 + iteration )
         else
            Ncall = two**( N_points(-1) - N_iterations2 + iteration )
            Ncall = two**( N_points(-1)) ! - N_iterations2 + iteration )
         endif
      CASE(WPAAj,WMAAj)
         if(sub_number.eq.2)then
            N_points(-1)=N_points(-1)-1
         endif                
         if ((N_iterations2-iteration).eq.0) then
            if (ps_number.eq.1) then
               Ncall = two**( N_points(-1)+1)
            else if (ps_number.eq.4) then
               Ncall = two**( N_points(-1)-2)
            else
               Ncall = two**( N_points(-1))
            endif
         elseif((N_iterations2-iteration).eq.1) then
            if (ps_number.eq.1) then
               Ncall = two**( N_points(-1))
            else if (ps_number.eq.4) then
               Ncall = two**( N_points(-1)-3)
            else
               Ncall = two**( N_points(-1)-1)
            endif
         elseif((N_iterations2-iteration).gt.1) then
            if (ps_number.eq.1) then
               Ncall = two**( N_points(-1)-1)
            else if (ps_number.eq.4) then
               Ncall = two**( N_points(-1)-4)
            else
               Ncall = two**( N_points(-1)-2)
            endif
         endif
         if(sub_number.eq.2)then
            N_points(-1)=N_points(-1)+1
         endif
      CASE DEFAULT
         Ncall = two**( N_points(-1) - N_iterations2 + iteration )
      END SELECT
      endif
      
      return
      
      end


********************************************************************************
********************************************************************************

      subroutine statistics_proc_FLoops(ps_number)
      implicit none
c include global variables via include-files, e.g. "procID"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/tau_hel.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/scales.inc"
c     variables for the integration loops
      integer iteration
      integer*8 Ncall, two
      integer ps_number
      parameter (two = 2)
      common/Ncall1/Ncall,iteration

      if (.not.lreaddata) then
c            print*,'lreaddata',lreaddata

         Select Case(procID)
      CASE(WW,WPhadWMlep,WPlepWMhad,AA,WWJ,WPHADWMJ,WPWMHADJ)
         Ncall = two**( N_points(0) -6 - N_iterations1 + iteration )
      CASE(ZZ,ZZhad,ZA,ZZJ)
         Ncall = two**( N_points(0) -7 - N_iterations1 + iteration )
      CASE DEFAULT
         Ncall = two**( N_points(0) - N_iterations1 + iteration )
      END SELECT
      endif
      return
      end

      subroutine AddToXsec(xsec,sdev2,time,sub,ps,nlo)
#include "VBFNLO/utilities/global.inc"
      double precision xsec,sdev2
      real time
      integer sub, ps, nlo
      
      final_xsec(sub,ps,nlo) = xsec
      final_sdev2(sub,ps,nlo) = sdev2
      final_time(sub,ps,nlo) = time

      final_xsec(sub,0,nlo) = final_xsec(sub,0,nlo) + xsec
      final_sdev2(sub,0,nlo)= final_sdev2(sub,0,nlo)+ sdev2
      final_time(sub,0,nlo) = final_time(sub,0,nlo) + time

      final_xsec(0,ps,nlo) = final_xsec(0,ps,nlo) + xsec
      final_sdev2(0,ps,nlo)= final_sdev2(0,ps,nlo)+ sdev2
      final_time(0,ps,nlo) = final_time(0,ps,nlo) + time

      final_xsec(0,0,nlo) = final_xsec(0,0,nlo) + xsec
      final_sdev2(0,0,nlo)= final_sdev2(0,0,nlo)+ sdev2
      final_time(0,0,nlo) = final_time(0,0,nlo) + time
      
      if(nlo.gt.0) then ! don't add LO to NLO result
        final_xsec(sub,ps,0) = final_xsec(sub,ps,0) + xsec
        final_sdev2(sub,ps,0)= final_sdev2(sub,ps,0)+ sdev2
        final_time(sub,ps,0) = final_time(sub,ps,0) + time
        
        final_xsec(sub,0,0) = final_xsec(sub,0,0) + xsec
        final_sdev2(sub,0,0)= final_sdev2(sub,0,0)+ sdev2
        final_time(sub,0,0) = final_time(sub,0,0) + time

        final_xsec(0,ps,0) = final_xsec(0,ps,0) + xsec
        final_sdev2(0,ps,0)= final_sdev2(0,ps,0)+ sdev2
        final_time(0,ps,0) = final_time(0,ps,0) + time

        final_xsec(0,0,0) = final_xsec(0,0,0) + xsec
        final_sdev2(0,0,0)= final_sdev2(0,0,0)+ sdev2
        final_time(0,0,0) = final_time(0,0,0) + time
      endif
      end

