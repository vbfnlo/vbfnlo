c------------------------------------------------------------
c
c     Calculates Kaluza-Klein - STM Gauge Boson Couplings
c
c-------------------------------------------------------------
c 
      subroutine kk_coupl_inp_man(vbfcalcswitch)
          use readinput
      
      implicit none

      logical vbfcalcswitch
c
c  Local Variables
      integer i,j,k,l
c 
c  Global Variables
#include "VBFNLO/utilities/coupl.inc"
#include "kkcoupl.inc"
#include "VBFNLO/utilities/global.inc"
      double precision CLR,XM2,XMG,B,V,A !difference SM<->KK
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

** default and temporary values
      double precision w1mass_def(0:kkw), w1mass_temp(0:kkz)
      double precision a1mass_temp(0:kka), a1mass_def(0:kka)
      double precision z1mass_def(0:kkw), z1mass_temp(0:kkz)
      double precision gwwz_temp(0:kkz), gwwz_def(0:kkz)
      double precision gwwa_temp(0:kka), gwwa_def(0:kka)
      double precision gwwww_temp(0:kkw), gwwww_def(0:kkw)
      double precision gwwzz_temp(0:kkz), gwwzz_def(0:kkz)
      double precision gwwaa_temp(0:kka), gwwaa_def(0:kka)
      double precision gwwaz_temp(0:kkz), gwwaz_def(0:kkz)

** dummy variable used for labelling couplings
      character*3 idstring1, idstring2, idstring3
      integer lenid1, lenid2, lenid3

      w1mass_k(0)=wmass
      z1mass_k(0)=zmass
      h1mass_k(0)=hmass
      a1mass_k(0)=0D0
      w1width_k(0)=wwidth
      z1width_k(0)=zwidth
      h1width_k(0)=hwidth
      gwwz_k(0,0,0)=gwwz
      gwwa_k(0,0,0)=gwwa
      hswitch=.false.

c--------------------------------------------------------------------
c  Read in Parameters from kk_coupl_input.dat
c--------------------------------------------------------------------
** set parameters 

      if (vbfcalcswitch) then
         call loadfileKK("kk_coupl_inp.dat",.true.)  
      else
         call loadfile("kk_coupl_inp.dat",.true.)  
      end if

      hswitchww = .false.
      hswitchwz = .false.
      hswitchzz = .false.
      hswitchwww = .false.
      hswitchwwz = .false.
      hswitchzzw = .false.
      hswitchzzz = .false.

      write(*,*)' '
      call read_int("KKMAXW",kkhw,1)
      call read_int("KKMAXZ",kkhz,1)
      call read_int("KKMAXG",kkha,1)
         
      kkhh = 0
      
      w1mass_def(0) = 500D0
      do i = 1, kkw
         w1mass_def(i) = w1mass_def(i-1) + 50D0
      end do
      call read_realList("KK_MASSW",kkhw,w1mass_temp,w1mass_def,.true.)
      do i=1,kkhw
         w1mass_k(i) = w1mass_temp(i-1)
      enddo
      
      z1mass_def(0) = 500D0
      do i = 1, kkz
         z1mass_def(i) = z1mass_def(i-1) + 50D0
      end do
      call read_realList("KK_MASSZ",kkhz,z1mass_temp,z1mass_def,.true.)
      do i=1,kkhz
         z1mass_k(i) = z1mass_temp(i-1)
      enddo
      
      do i = 0, kka
         a1mass_def(i) = 0D0
      end do
      if (.not. mhlmswitch) then
         call read_realList("KK_MASSG",kkha,a1mass_temp,a1mass_def,.true.)
         do i=1,kkha
            a1mass_k(i) = a1mass_temp(i-1)
         enddo
      else
         do i=1,kkha
            a1mass_k(i) = a1mass_def(i-1)
         enddo
      end if

** reading in the wwz coupling:     
      do i = 0, kkz
         gwwz_def(i) = 0D0
      end do
      call read_realList("CPL_W0-W0-ZX",kkhz+1,gwwz_temp,gwwz_def,
     -     .true.)
      do i = 0, kkhz
         gwwz_k(0,0,i) = gwwz_temp(i)
      end do
      
      do i = 1, kkhw
         do j = 0, i
            write( idstring1, '(i3)' )  i!+1
            do lenid1 = 1, 3
               if (idstring1(lenid1:lenid1).ne.' ') then
                  exit
               endif
            enddo
            write( idstring2, '(i3)' )  j!+1
            do lenid2 = 1, 3
               if (idstring2(lenid2:lenid2).ne.' ') then
                  exit
               endif
            enddo
            call read_realList("CPL_W"//idstring1(lenid1:3)//
     -           "-W"//idstring2(lenid2:3)//"-ZX",kkhz+1,gwwz_temp,
     -           gwwz_def,.true.)
            do k = 0, kkhz
               gwwz_k(i,j,k) = gwwz_temp(k)
            end do
         end do
      end do
      
** read in the wwz' coupling   
      do i = 0, kka
         gwwa_def(i) = 0D0
      end do
      call read_realList("CPL_W0-W0-GX",kkha,gwwa_temp,gwwa_def,
     -     .true.)
      do i = 1, kkha
         gwwa_k(0,0,i) = gwwa_temp(i-1)
      end do
      
      do i = 1, kkhw
         do j = 0, i
            write( idstring1, '(i3)' )  i!+1
            do lenid1 = 1, 3
               if (idstring1(lenid1:lenid1).ne.' ') then
                  exit
               endif
            enddo
            write( idstring2, '(i3)' )  j!+1
            do lenid2 = 1, 3
               if (idstring2(lenid2:lenid2).ne.' ') then
                  exit
               endif
            enddo
            call read_realList("CPL_W"//idstring1(lenid1:3)//
     -           "-W"//idstring2(lenid2:3)//"-GX",kkha+1,gwwa_temp,
     -           gwwa_def,.true.)
            do k = 0, kkha
               gwwa_k(i,j,k) = gwwa_temp(k)
            end do
         end do
      end do
      
** read-in wwww coupling
      do i = 0, kkw
         gwwww_def(i) = 0D0
      end do
      do i=0,kkhw
         do j=0,i
            do k=0,j
               write( idstring1, '(i3)' )  i!+1
               do lenid1 = 1, 3
                  if (idstring1(lenid1:lenid1).ne.' ') then
                     exit
                  endif
               enddo
               write( idstring2, '(i3)' )  j!+1
               do lenid2 = 1, 3
                  if (idstring2(lenid2:lenid2).ne.' ') then
                     exit
                  endif
               enddo
               write( idstring3, '(i3)' )  k!+1
               do lenid3 = 1, 3
                  if (idstring3(lenid3:lenid3).ne.' ') then
                     exit
                  endif
               enddo
               
               call read_realList("CPL_W"//idstring1(lenid1:3)//
     -              "-W"//idstring2(lenid2:3)//"-W"//idstring3(lenid3:3)//
     -              "-WX",k+1,gwwww_temp,gwwww_def,.true.)
               do l=0,k
                  gwwww_k(i,j,k,l) = gwwww_temp(l)
               enddo
            enddo
         enddo
      enddo
         
** read-in the wwzz couplings 
      do i = 0, kkz
         gwwzz_def(i) = 0D0
      end do     
      
      do k=0,kkhz
         write( idstring1, '(i3)' )  k!+1
         do lenid1 = 1, 3
            if (idstring1(lenid1:lenid1).ne.' ') then
               exit
            endif
         enddo
         call read_realList("CPL_W0-W0-Z"//idstring1(lenid1:3)//
     -        "-ZX",k+1,gwwzz_temp,gwwzz_def,.true.)
         do l = 0, k
            gwwzz_k(0,0,k,l) = gwwzz_temp(l)
         end do
      end do
      
      do i = 1, kkhw
         do j = 0, i
            do k = 0, kkhz
               write( idstring1, '(i3)' )  i!+1
               do lenid1 = 1, 3
                  if (idstring1(lenid1:lenid1).ne.' ') then
                     exit
                  endif
               enddo
               write( idstring2, '(i3)' )  j!+1
               do lenid2 = 1, 3
                  if (idstring2(lenid2:lenid2).ne.' ') then
                     exit
                  endif
               enddo
               write( idstring3, '(i3)' )  k!+1
               do lenid3 = 1, 3
                  if (idstring3(lenid3:lenid3).ne.' ') then
                     exit
                  endif
               enddo
               
               call read_realList("CPL_W"//idstring1(lenid1:3)//
     -              "-W"//idstring2(lenid2:3)//"-Z"//idstring3(lenid3:3)//
     -              "-ZX",k+1,gwwzz_temp,gwwzz_def,.true.)
               do l=0,k
                  gwwzz_k(i,j,k,l) = gwwzz_temp(l)
               enddo
            end do
         end do
      end do
         
** read in the wwz'z' couplings
      do i = 0, kka
         gwwaa_def(i) = 0D0
      end do     
      
      do k = 0, kkha
         write( idstring1, '(i3)' )  k
         do lenid1 = 1, 3
            if (idstring1(lenid1:lenid1).ne.' ') then
               exit
            endif
         enddo
         call read_realList("CPL_W0-W0-G"//idstring1(lenid1:3)//
     -        "-GX",k+1,gwwaa_temp,gwwaa_def,.true.)
         do l = 0, k
            gwwaa_k(0,0,k,l) = gwwaa_temp(l)
         end do
      end do
      
      do i = 1, kkhw
         do j = 0, i
            do k = 0, kkha
               write( idstring1, '(i3)' )  i!+1
               do lenid1 = 1, 3
                  if (idstring1(lenid1:lenid1).ne.' ') then
                     exit
                  endif
               enddo
               write( idstring2, '(i3)' )  j!+1
               do lenid2 = 1, 3
                  if (idstring2(lenid2:lenid2).ne.' ') then
                     exit
                  endif
               enddo
               write( idstring3, '(i3)' )  k
               do lenid3 = 1, 3
                  if (idstring3(lenid3:lenid3).ne.' ') then
                     exit
                  endif
               enddo
               
               call read_realList("CPL_W"//idstring1(lenid1:3)//
     -              "-W"//idstring2(lenid2:3)//"-G"//idstring3(lenid3:3)//
     -              "-GX",k+1,gwwaa_temp,gwwaa_def,.true.)
               do l = 0, k
                  gwwaa_k(i,j,k,l) = gwwaa_temp(l)
               end do
            end do
         end do
      end do

** read-in wwz'z couplings
      do i = 0, kkz
         gwwaz_def(i) = 0D0
      end do     
      
      do i = 0, kkhw
         do j = 0, i
            do k = 0, kkha
               write( idstring1, '(i3)' )  i!+1
               do lenid1 = 1, 3
                  if (idstring1(lenid1:lenid1).ne.' ') then
                     exit
                  endif
               enddo
               write( idstring2, '(i3)' )  j!+1
               do lenid2 = 1, 3
                  if (idstring2(lenid2:lenid2).ne.' ') then
                     exit
                  endif
               enddo
               write( idstring3, '(i3)' )  k
               do lenid3 = 1, 3
                  if (idstring3(lenid3:lenid3).ne.' ') then
                     exit
                  endif
               enddo
               call read_realList("CPL_W"//idstring1(lenid1:3)//
     -              "-W"//idstring2(lenid2:3)//"-G"//idstring3(lenid3:3)//
     -              "-ZX",kkhz+1,gwwaz_temp,gwwaz_def,.true.)
               
               do l = 0, kkhz
                  gwwaz_k(i,j,k,l) = gwwaz_temp(l)
               end do
            end do
         end do
      end do
         
      
      call closefile

c$$$  OLD READ-IN ROUTINE
c$$$      open (unit=50,file="kk_coupl_inp.dat")
c$$$      read(50,*) hswitchww
c$$$      read(50,*) hswitchwz
c$$$      read(50,*) hswitchzz
c$$$      read(50,*) hswitchwww
c$$$      read(50,*) hswitchwwz
c$$$      read(50,*) hswitchzzw
c$$$      read(50,*) hswitchzzz
c$$$      read(50,*) kkhw
c$$$      read(50,*) kkhz
c$$$      read(50,*) kkha
c$$$      read(50,*) kkhh
c$$$
c$$$      read (50,*) gtext
c$$$      read (50,*) gtext
c$$$      do i=1,kkhw
c$$$      read (50,*)  w1mass_k(i)
c$$$      enddo
c$$$c--------------------------------------------------------------------
c$$$      read (50,*) gtext
c$$$      do i=1,kkhz
c$$$      read (50,*)  z1mass_k(i)
c$$$      enddo
c$$$c--------------------------------------------------------------------
c$$$      read (50,*) gtext
c$$$      do i=1,kkha
c$$$      read (50,*)  a1mass_k(i)
c$$$      enddo
c$$$c--------------------------------------------------------------------
c$$$      if (hswitch) then
c$$$      read (50,*) gtext
c$$$      do i=0,kkhh
c$$$      read (50,*)  h1mass_k(i)
c$$$      read (50,*)  h1width_k(i)
c$$$      enddo
c$$$      endif
c$$$c--------------------------------------------------------------------
c$$$      read (50,*) gtext
c$$$      do i=0,kkhz
c$$$          read (50,*)  gwwz_k(0,0,i)
c$$$      enddo
c$$$      do i=1,kkhw
c$$$          do j=0,i
c$$$             do k=0,kkhz
c$$$              read (50,*)  gwwz_k(i,j,k) 
c$$$             enddo
c$$$          enddo
c$$$      enddo
c$$$c--------------------------------------------------------------------
c$$$      read (50,*) gtext
c$$$      do i=1,kkha
c$$$          read (50,*)  gwwa_k(0,0,i) 
c$$$      enddo
c$$$      do i=1,kkhw
c$$$          do j=0,i
c$$$             do k=0,kkha
c$$$              read (50,*)  gwwa_k(i,j,k)
c$$$             enddo
c$$$          enddo
c$$$      enddo
c$$$c--------------------------------------------------------------------
c$$$      if (hswitch) then
c$$$      read (50,*) gtext
c$$$      do i=0,kkhw
c$$$          do j=0,i
c$$$             do k=0,kkhh
c$$$              read (50,*)  gwwh_k(i,j,k) 
c$$$             enddo
c$$$          enddo
c$$$      enddo
c$$$      endif
c$$$c--------------------------------------------------------------------
c$$$      if (hswitch) then
c$$$      read (50,*) gtext
c$$$      do i=0,kkhz
c$$$          do j=0,i
c$$$             do k=0,kkhh
c$$$              read (50,*)  gzzh_k(i,j,k) 
c$$$             enddo
c$$$          enddo
c$$$      enddo
c$$$      endif
c$$$c--------------------------------------------------------------------
c$$$      read (50,*) gtext
c$$$      do i=0,kkhw
c$$$          do j=0,i
c$$$             do k=0,j
c$$$                do l=0,k
c$$$                  read (50,*)  gwwww_k(i,j,k,l)
c$$$                 enddo
c$$$             enddo
c$$$          enddo
c$$$      enddo
c$$$c--------------------------------------------------------------------
c$$$      read (50,*) gtext
c$$$      do k=0,kkhz
c$$$         do l=0,k
c$$$             read (50,*)  gwwzz_k(0,0,k,l) 
c$$$         enddo
c$$$      enddo
c$$$      do i=1,kkhw
c$$$          do j=0,i
c$$$             do k=0,kkhz
c$$$                do l=0,k
c$$$                  read (50,*)  gwwzz_k(i,j,k,l) 
c$$$                enddo
c$$$             enddo
c$$$          enddo
c$$$      enddo
c$$$c--------------------------------------------------------------------
c$$$      read (50,*) gtext
c$$$      do k=0,kkha
c$$$         do l=0,k
c$$$             read (50,*)  gwwaa_k(0,0,k,l)
c$$$         enddo
c$$$      enddo
c$$$      do i=1,kkhw
c$$$          do j=0,i
c$$$             do k=0,kkha
c$$$                do l=0,k
c$$$                  read (50,*)  gwwaa_k(i,j,k,l)
c$$$                 enddo
c$$$             enddo
c$$$          enddo
c$$$      enddo
c$$$c--------------------------------------------------------------------
c$$$      read (50,*) gtext
c$$$      do i=0,kkhw
c$$$          do j=0,i
c$$$             do k=0,kkha
c$$$                do l=0,kkhz
c$$$                  read (50,*)  gwwaz_k(i,j,k,l) 
c$$$                enddo
c$$$             enddo
c$$$          enddo
c$$$      enddo
c$$$      close (50)


      gwwz=gwwz_k(0,0,0)
c
c  redefine koppln.f gauge boson coupling gwwz
c      
      B(3,4,2) = -gwwz_k(0,0,0)
      B(4,2,3) = -gwwz_k(0,0,0)
      B(2,3,4) = -gwwz_k(0,0,0)
      B(2,4,3) =  gwwz_k(0,0,0)
      B(4,3,2) =  gwwz_k(0,0,0)
      B(3,2,4) =  gwwz_k(0,0,0)
c
c  redefine Higgs mass for phasespace adatption
c      
      xm2(6)=w1mass_k(1)**2
      Print*,"Higgs mass squared is set to ", xm2(6) , 
     -     " for phasespace generator." 
c
c  Set up boson widths
c 
      call widthkk
      xmg(6)=w1mass_k(1)*w1width_k(1)*10d0
      Print*,"Higgs width*mass information set to ", xmg(6) , 
     -     " for phasespace generator." 
c 


      END


      subroutine kkcouplcheck
      
      implicit none
c
c  Local Variables
      integer i,j,k
      double precision h1,h2,pi
      parameter(pi =3.141592653589793d0)      
c 
c  Global Variables
#include "VBFNLO/utilities/coupl.inc"
#include "kkcoupl.inc"
      open(unit=40, file="kkcheck.dat")
      write(40,*) "Tower height W+/W-   ", kkhw
      write(40,*) "Tower height Z ", kkhz
      write(40,*) "Tower height Z'", kkha
      write(40,*) "Tower height Higgs   ", kkhh
      write(40,*) "----------------------------------------------------"
      do i=0,kkhw
      write (40,*) "W",i+1," -mass  ", w1mass_k(i)
      write (40,*) "W",i+1," -width ", w1width_k(i)
      enddo
      write(40,*) "----------------------------------------------------"
      do i=0,kkhz
      write (40,*)  "Z'",i+1," -mass  ",z1mass_k(i)
      write (40,*)  "Z'",i+1," -width ",z1width_k(i)
      enddo
      write(40,*) "----------------------------------------------------"
      do i=0,kkha
      write (40,*)  "Gamma'",i," -mass ", a1mass_k(i)
      write (40,*)  "Gamma'",i," -width ",a1width_k(i)
      enddo
      write(40,*) "----------------------------------------------------"
      if (hswitch) then
      do i=1,kkhh
      write (40,*)  "H",i," -mass  ", h1mass_k(i)
      write (40,*)  "H",i," -width ", h1width_k(i)
      enddo
      write(40,*) "----------q-----------------------------------------"
      endif
      do i=0,kkhz
      write (40,*)  "W 1 - W 1 - Z",i+1," -Coupling ",gwwz_k(0,0,i) 
      enddo
!       do i=1,kkhw
!           do j=0,kkhw
!              do k=0,kkhz
!               write (40,*) "W",i+1," - W",j+1," - Z",k+1," -Coupling ", gwwz_k(i,j,k) 
!              enddo
!           enddo
!       enddo
      do i=1,kkhw
              write (40,*) "W",i+1," - W 1 - Z 1 -Coupling ", gwwz_k(i,0,0) 
      enddo
      write(40,*) "----------------------------------------------------"
      do i=0,kkha
          write (40,*) "W 1 - W 1 - A",i," -Coupling ",gwwa_k(0,0,i) 
      enddo
!       do i=1,kkhw
!           do j=0,kkhw
!              do k=0,kkha
!               write (40,*) "W",i+1," - W",j+1," - A",k," -Coupling ",gwwa_k(i,j,k)
!              enddo
!           enddo
!       enddo
      write(40,*) "----------------------------------------------------"
      if(hswitch) then
      do i=0,kkhw
          do j=0,i
             do k=0,kkhh
              write (40,*)  "W",i+1," - W",j+1," - H",k+1," -Coupling ",  gwwh_k(i,j,k) 
             enddo
          enddo
      enddo
      write(40,*) "----------------------------------------------------"
      do i=0,kkhz
          do j=0,i
             do k=0,kkhh
              write (40,*) "Z",i+1," - Z",j+1," - H",k+1," -Coupling ", gzzh_k(i,j,k) 
             enddo
          enddo
      enddo
      write(40,*) "----------------------------------------------------"
      endif
!       do i=0,kkhw
!           do j=0,i
!              do k=0,j
!                 do l=0,k
!                   write (40,*) "W",i+1," - W",j+1," - W",k+1," - W",l+1," -Coupling ",  gwwww_k(i,j,k,l) 
!                 enddo
!              enddo
!           enddo
!       enddo
      write (40,*) "W",1," - W",1," - W",1," - W",1," -Coupling ",  gwwww_k(0,0,0,0)
      write(40,*) "----------------------------------------------------"
!       do k=0,kkhz
!          do l=0,k
             write (40,*) "W 1 - W 1 - Z",1," - Z",1," -Coupling ", gwwzz_k(0,0,0,0) 
!          enddo
!       enddo
!       do i=1,kkhw
!           do j=0,i
!              do k=0,kkhz
!                 do l=0,k
!                   write (40,*) "W",i+1," - W",j+1," - W",k+1," - W",l+1," -Coupling ",  gwwzz_k(i,j,k,l) 
!                 enddo
!              enddo
!           enddo
!       enddo
      write(40,*) "----------------------------------------------------"
!       do k=0,kkhz
!          do l=0,k
             write (40,*) "W 1 - W 1 - A",0," - A",0," -Coupling ", gwwaa_k(0,0,0,0) 
!          enddo
!       enddo
!       do i=1,kkhw
!           do j=0,i
!              do k=0,kkhz
!                 do l=0,k
!                  write (40,*) "W",i+1," - W",j+1," - A",k," - A",l," -Coupling ", gwwaa_k(i,j,k,l) 
!                 enddo
!              enddo
!           enddo
!       enddo
      write(40,*) "----------------------------------------------------"
!       do i=0,kkhw
!           do j=0,i
!              do k=0,kkhz
!                 do l=0,kkhz
!                   write (40,*) "W",i+1," - W",j+1," - A",k," - Z",l+1," -Coupling ", gwwaz_k(i,j,k,l) 
!                 enddo
!              enddo
!           enddo
!       enddo

      write (40,*) "W",1," - W",1," - A",0," - Z",1," -Coupling ", gwwaz_k(0,0,0,0) 
      write(40,*) "----------------------------------------------------"
      write(40,*) "----------------Sum - Rules ------------------------"
      write(40,*) "----------------------------------------------------"
      h1=gwwww_k(0,0,0,0)**2
	  h2=4*wmass**2*gwwww_k(0,0,0,0)**2
	  do i=0,kkhz
      h1=h1-gwwz_k(0,0,i)**2
      h2=h2-3*z1mass_k(i)**2*gwwz_k(0,0,i)**2
	  enddo
      do i=0,kkha
      h1=h1-gwwa_k(0,0,i)**2
      h2=h2-3*a1mass_k(i)**2*gwwa_k(0,0,i)**2
	  enddo
      if (hswitchww) then
      do i=0,kkhh
      h2=h2-wmass**2*gwwh_k(0,0,i)**2
      enddo
      endif
      write(40,*) "W1 W1 -> W1 W1: ~E^4", abs(h1/gwwww_k(0,0,0,0)**2)*100, " %", "     ",
     &            " ~E^2", abs(h2/(4*w1mass_k(0)**2*gwwww_k(0,0,0,0)**2))*100,  " %", "     "
      write(40,*) "----------------------------------------------------"
      h1=gwwzz_k(0,0,0,0)
      h2=2*(z1mass_k(0)**2+w1mass_k(0)**2)*gwwzz_k(0,0,0,0)
      do i=0,kkhw
      h1=h1-gwwz_k(i,0,0)**2
      h2=h2-3*(w1mass_k(i)**2)*gwwz_k(i,0,0)**2+((zmass**2-wmass**2)**2/w1mass_k(i)**2)*gwwz_k(i,0,0)**2
      enddo
      if (hswitchww) then
      do i=0,kkhh
      h2=h2-wmass**2*gwwh_k(0,0,i)*gzzh_k(0,0,i)
      enddo
      endif
      write(40,*) "W1 Z1 -> W1 Z1: ~E^4", abs(h1/gwwzz_k(0,0,0,0))*100, " %", "     ",
     &            " ~E^2", abs(h2/((2*z1mass_k(0)**2+2*w1mass_k(0)**2)*gwwzz_k(0,0,0,0)))*100,  " %", "     "
      write(40,*) "----------------------------------------------------"
      END


!    calculate widths of kk-bosons
!    still using standard width for w+-, z and h

      subroutine widthkk
      IMPLICIT NONE

c    in:       l:  decaying boson W/Z
c              i:  KK-Index of decaying boson
c              j:  KK-Index of decay product
c              k:  kk-Index of decay product
c    out:      width
c
#include "kkcoupl.inc"
#include "coupl.inc"

      integer i,j,k
      double precision pi
      Parameter(pi=3.14159265358979)

      do i=1,kkhw
      w1width_k(i)=0d0
      enddo
      do j=1,kkha
      a1width_k(j)=0d0
	  enddo
	  do j=1,kkhz
      z1width_k(j)=0d0
      enddo

      do i=1,kkhw
      do j=0,i-1
      do k=0,kkhz
      if (w1mass_K(i).ge. (w1mass_k(j)+z1mass_k(k))) then

!       w1width_k(i)=w1width_k(i)+(gwwz_k(i,j,k)**2*(((w1mass_k(i) - w1mass_k(j) - z1mass_k(k))*
!      &  (w1mass_k(i) + w1mass_k(j) - z1mass_k(k))*(w1mass_k(i) - w1mass_k(j) + 
!      &   z1mass_k(k))*(w1mass_k(i) + w1mass_k(j) + z1mass_k(k)))/w1mass_k(i)**2)**1.5D0*
!      & -    (w1mass_k(i)**4 + w1mass_k(j)**4 + 10*w1mass_k(j)**2*z1mass_k(k)**2 + 
!      &    z1mass_k(k)**4 + 10*w1mass_k(i)**2*(w1mass_k(j)**2 + z1mass_k(k)**2)))/
!      &    -  (96D0*w1mass_k(i)*w1mass_k(j)**2*z1mass_k(k)**2*(Sqrt((w1mass_k(i)**2 + w1mass_k(j)**2 
!      &  - z1mass_k(k)**2)**2/w1mass_k(i)**2) + Sqrt((w1mass_k(i)**2 - w1mass_k(j)**2 + z1mass_k(k)**2)**2/w1mass_k(i)**2))*Pi)

      w1width_k(i)=w1width_k(i)+(gwwz_k(i,j,k)**2*(((w1mass_k(i) - w1mass_k(j) - z1mass_k(k))
     &    *(w1mass_k(i) + w1mass_k(j) - z1mass_k(k))*
     &    (w1mass_k(i) - w1mass_k(j) + z1mass_k(k))*(w1mass_k(i) + w1mass_k(j) + z1mass_k(k)))/w1mass_k(i)**2)**1.5D0*
     &    (w1mass_k(i)**4 + 4*w1mass_k(j)**4 - 2*w1mass_k(j)**2*z1mass_k(k)**2 
     &    + z1mass_k(k)**4 - 2*w1mass_k(i)**2*(w1mass_k(j)**2 - 5*z1mass_k(k)**2)))/
     &  (  96D0*w1mass_k(i)*w1mass_k(j)**2*z1mass_k(k)**2*(sqrt((w1mass_k(i)**2 
     &  + w1mass_k(j)**2 - z1mass_k(k)**2)**2/w1mass_k(i)**2)+ sqrt((w1mass_k(i)**2 - w1mass_k(j)**2 
     &   + z1mass_k(k)**2)**2/w1mass_k(i)**2))*Pi)

      endif
      enddo
      enddo
      enddo

      do i=1,kkhw
      do j=0,i-1
      do k=1,kkha
      if (w1mass_K(i).ge. (w1mass_k(j)+a1mass_k(k))) then

!       w1width_k(i)=w1width_k(i)+(gwwa_k(i,j,k)**2*(((w1mass_k(i) - w1mass_k(j) - a1mass_k(k))*
!      &  (w1mass_k(i) + w1mass_k(j) - a1mass_k(k))*(w1mass_k(i) - w1mass_k(j) + 
!      &   a1mass_k(k))*(w1mass_k(i) + w1mass_k(j) + a1mass_k(k)))/w1mass_k(i)**2)**1.5D0*
!      & -    (w1mass_k(i)**4 + w1mass_k(j)**4 + 10*w1mass_k(j)**2*a1mass_k(k)**2 + 
!      &    a1mass_k(k)**4 + 10*w1mass_k(i)**2*(w1mass_k(j)**2 + a1mass_k(k)**2)))/
!      &    -  (96D0*w1mass_k(i)*w1mass_k(j)**2*a1mass_k(k)**2*(Sqrt((w1mass_k(i)**2 + w1mass_k(j)**2 
!      &  - a1mass_k(k)**2)**2/w1mass_k(i)**2) + Sqrt((w1mass_k(i)**2 - w1mass_k(j)**2 + a1mass_k(k)**2)**2/w1mass_k(i)**2))*Pi)
    
      w1width_k(i)=w1width_k(i)+(gwwa_k(i,j,k)**2*(((w1mass_k(i) - w1mass_k(j) - a1mass_k(k))
     &    *(w1mass_k(i) + w1mass_k(j) - a1mass_k(k))*
     &    (w1mass_k(i) - w1mass_k(j) + a1mass_k(k))*(w1mass_k(i) + w1mass_k(j) + a1mass_k(k)))/w1mass_k(i)**2)**1.5D0*
     &    (w1mass_k(i)**4 + 4*w1mass_k(j)**4 - 2*w1mass_k(j)**2*a1mass_k(k)**2 
     &    + a1mass_k(k)**4 - 2*w1mass_k(i)**2*(w1mass_k(j)**2 - 5*a1mass_k(k)**2)))/
     &  (  96D0*w1mass_k(i)*w1mass_k(j)**2*a1mass_k(k)**2*(sqrt((w1mass_k(i)**2 
     &  + w1mass_k(j)**2 - a1mass_k(k)**2)**2/w1mass_k(i)**2)+ sqrt((w1mass_k(i)**2 - w1mass_k(j)**2 
     &   + a1mass_k(k)**2)**2/w1mass_k(i)**2))*Pi)


      endif
      enddo
      enddo
      enddo


      do i=1,kkhz
      do k=0,kkhw
      do j=0,k
      if (z1mass_K(i).ge. (w1mass_k(j)+w1mass_k(k))) then

      z1width_k(i)=z1width_k(i)+(gwwz_k(k,j,i)**2*(((w1mass_k(j) - w1mass_k(k) - z1mass_k(i))*(w1mass_k(j) + w1mass_k(k) 
     &   - z1mass_k(i))*(w1mass_k(j) - w1mass_k(k) + z1mass_k(i))*(w1mass_k(j) + w1mass_k(k) + 
     &  z1mass_k(i)))/z1mass_k(i)**2)**1.5D0*
     &  (- (w1mass_k(j)**4 + w1mass_k(k)**4 + 10*w1mass_k(k)**2*z1mass_k(i)**2 + 
     &  z1mass_k(i)**4 + 10*w1mass_k(j)**2*(w1mass_k(k)**2 + z1mass_k(i)**2))))/
     &  (-  (96D0*w1mass_k(j)**2*w1mass_k(k)**2*z1mass_k(i)*(Sqrt((w1mass_k(j)**2 
     &  - w1mass_k(k)**2 + z1mass_k(i)**2)**2/z1mass_k(i)**2) + 
     &  Sqrt((-w1mass_k(j)**2 + w1mass_k(k)**2 + z1mass_k(i)**2)**2/z1mass_k(i)**2))*Pi))

      endif
      enddo
      enddo
      enddo
      do i=1,kkha
      do k=0,kkhw
      do j=0,k
      if (a1mass_K(i).ge. (w1mass_k(j)+w1mass_k(k))) then

      a1width_k(i)=a1width_k(i)+((gwwa_k(k,j,i)**2*(((w1mass_k(j) - w1mass_k(k) - a1mass_k(i))*(w1mass_k(j) + w1mass_k(k) 
     &   - a1mass_k(i))*(w1mass_k(j) - w1mass_k(k) + a1mass_k(i))*(w1mass_k(j) + w1mass_k(k) + 
     &  a1mass_k(i)))/a1mass_k(i)**2)**1.5D0*
     &  (- (w1mass_k(j)**4 + w1mass_k(k)**4 + 10*w1mass_k(k)**2*a1mass_k(i)**2 + 
     &  a1mass_k(i)**4 + 10*w1mass_k(j)**2*(w1mass_k(k)**2 + a1mass_k(i)**2))))/
     &  (-  (96D0*w1mass_k(j)**2*w1mass_k(k)**2*a1mass_k(i)*(Sqrt((w1mass_k(j)**2 
     &  - w1mass_k(k)**2 + a1mass_k(i)**2)**2/a1mass_k(i)**2) + 
     &  Sqrt((-w1mass_k(j)**2 + w1mass_k(k)**2 + a1mass_k(i)**2)**2/a1mass_k(i)**2))*Pi)))

      endif
      enddo
      enddo
      enddo
      END


********************************************************************************
********************************************************************************
     
      subroutine kk_coupl_inp
          use readinput
      
      implicit none
      
c  Local Variables
      double precision kk1,kk2,kk3,kk4,rdwn
      double precision wpmass
      integer*4 k,kkmaxw,kkmaxz,kkmaxg
      logical vbfcalcswitch

c
c  Global Variables
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/kkcoupl.inc"
#include "VBFNLO/utilities/global.inc"

c -----------------------------------------
      call loadfile("kk_input.dat",.true.)

      call read_logical("VBFCALC_SWITCH",vbfcalcswitch,.true.) 
      call read_logical("MHLM_SWITCH",mhlmswitch,.false.)
      if (mhlmswitch) then
         call read_Real("WPMASS",wpmass,500D0)
      else
         call read_Real("RDWN",rdwn,9.75D-9)
         call read_Int("KKMAXW",kkmaxw,1)
         call read_Int("KKMAXZ",kkmaxz,1)
         call read_Int("KKMAXG",kkmaxg,1)
      end if

      call closefile

      kkmaxw=kkmaxw!+1
      kkmaxz=kkmaxz!+1
      kkmaxg=kkmaxg!+1

      if (mhlmswitch) then
         Print*, "**************************************************"
         Print*, "Initiated 3-Site Higgsless Calculation"
         Print*, "**************************************************"
         Print*, "ALL PREVIOUS HIGGS-RELATED INFORMATION IS OBSOLETE"
         Print*, "**************************************************"
         if (vbfcalcswitch) then
            Print*, "Couplings and Z' mass are calculated by VBFNLO"
            Print*, "Heavy Fermions and Couplings of heavy bosons to"
            Print*, "light Fermions are excluded!"
            Print*, "Kaluza-Klein related Information is written to"
            Print*, "kk_coupl_inp.dat:  Parameters included in the "
            Print*, "proceeding run"
            Print*, "kkcheck.dat:  Kaluza-Klein information and "
            Print*, "relative sum rule violation"
            Print*, "**************************************************"
            call MHLM_coupl(wmass,gwwa,wpmass,zmass)
            call mhlm_fermion(wpmass)
         else                   !vbfcalcswitch
            Print*, "MANUAL INPUT OF COUPLINGS MASSES"
            Print*, "Kaluza-Klein related Information is written to"
            Print*, "kk_coupl_inp.dat:  Parameters included in the "
            Print*, "proceeding run"
            Print*, "kkcheck.dat:  Kaluza-Klein information and"
            Print*, "relative sum rule violation"
            Print*, "**************************************************"
            call mhlm_fermion(wpmass)
         endif                  !vbfclacswitch
      else                      !mhlmswitch
         
         Print*, "**************************************************"
         Print*, "Initiated Warped Higgsless Calculation"
         Print*, "**************************************************"
         Print*, "ALL PREVIOUS HIGGS-RELATED INFORMATION IS OBSOLETE"
         Print*, "**************************************************"
         if (vbfcalcswitch) then
            Print*, "Couplings and Kaluza-Klein masses are calculated" 
            Print*, "by VBFNLO"
            Print*, "Kaluza-Klein related Information is written to"
            Print*, "kk_coupl_inp.dat:  Parameters included in the "
            Print*, "proceeding run"
            Print*, "kkcheck.dat:  Kaluza-Klein information and "
            Print*, "relative sum rule violation"
            Print*, "**************************************************"
c     C++ Calculation of Couplings and Masses
            Call calckkcoupl( gwwa , wmass, zmass, rdwn, kkmaxw, kkmaxz, kkmaxg )
         else
            Print*, "MANUAL INPUT OF COUPLINGS MASSES"
            Print*, "Kaluza-Klein related Information is written to"
            Print*, "kk_coupl_inp.dat:  Parameters included in the "
            Print*, "proceeding run"
            Print*, "kkcheck.dat:  Kaluza-Klein information and "
            Print*, "relative sum rule violation"
            Print*, "**************************************************"
        ENDIF
      endif
c
c  Read in Parameters from kk_input.dat
c     
      call kk_coupl_inp_man(vbfcalcswitch)


c calculate Z' couplings if necessary
      if (mgtest) then
         call setmgtest
      endif
c     
c
c  Sum-Rules WW -> WW
c

      if (.not. mhlmswitch) then
         kk1=0d0
         kk2=0d0
         do k=0,kkhz-1
            kk1= gwwz_k(0,0,k)**2+kk1
            kk2=kk2+3*z1mass_k(k)**2*gwwz_k(0,0,k)**2
         enddo
         do k=0,kkha
            kk1= gwwa_k(0,0,k)**2+kk1
            kk2=kk2+3*a1mass_k(k)**2*gwwa_k(0,0,k)**2
         enddo
         
         gwwz_k(0,0,kkhz)=sqrt((4d0*wmass**2*kk1-kk2)/(3d0*z1mass_k(kkhz)**2-4*wmass**2))
         
         gwwww_k(0,0,0,0)=sqrt(kk1+gwwz_k(0,0,kkhz)**2)
	  
c     
c  Sum-Rules WZ -> WZ
c
         kk1=0d0
         do k=0,kkhw-1
            kk1=kk1+2*(wmass**2+zmass**2)*gwwz_k(k,0,0)**2
         enddo
         kk2=0d0
         kk3=zmass**2-wmass**2
         do k=0,kkhw-1
            kk2=kk2+gwwz_k(k,0,0)**2*(3*w1mass_k(k)**2-(kk3**2/w1mass_k(k)**2))
         enddo
         kk4=3*w1mass_k(kkhw)**2 - (kk3**2/w1mass_k(kkhw)**2) -2d0*(wmass**2+zmass**2)
         gwwz_k(kkhw,0,0)=sqrt((kk1-kk2)/kk4)
         
         gwwzz_k(0,0,0,0)=0d0
         do k=0,kkhw
            gwwzz_k(0,0,0,0)=gwwzz_k(0,0,0,0)+gwwz_k(k,0,0)**2
         enddo	  
c
c
c  Sum-Rules AZ -> WW
c      
      gwwaz_k(0,0,0,0)=gwwa*gwwz
c
c  Sum-Rules AA -> WW
c      
      gwwaa_k(0,0,0,0)=gwwa**2

	endif ! mhlmswitch for sum rules check

c  Dump Kaluza-Klein Information
c 
      call kkcouplcheck


      END
