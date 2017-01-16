
      subroutine read_spin2coupl
          use readinput
c reads in initial set of spin-2 parameters
c from file spin2coupl.dat. The parameters
c are stored in common blocks in spin2coupl.inc      
      implicit none

#include "global.inc"   
#include "spin2coupl.inc"   

      print *," "
      print *,"  Information on spin-2 parameters  "
      print *,"------------------------------------------------"
      call loadfile("spin2coupl.dat",.true.)
      call read_real("F1",f1_0,0.0d0)
      call read_real("F2",f2_0,0.0d0)
      call read_real("F3",f3_0,0.0d0)
      call read_real("F4",f4_0,0.0d0)
      call read_real("F5",f5_0,0.0d0)
      call read_real("F9",f9_0,0.0d0)
      call read_real("F6",f6_0,0.0d0)
      call read_real("F7",f7_0,0.0d0)
      call read_real("LAMBDA",lambda_0,0.0d0)
      call read_real("SP2MASS",sp2mass_0,0.0d0)
      call read_real("MSP2TRIPPM",msp2trippm_0,0.0d0)
      call read_real("MSP2TRIPN",msp2tripn_0,0.0d0)
      call read_real("BRRAT",brrat_0,1.0d0)
      call read_real("BRRATTRIPPM",brrattrippm_0,1.0d0)
      call read_real("BRRATTRIPN",brrattripn_0,1.0d0)
      call read_real("NFF",nff_0,0.0d0)
      call read_real("LAMBDAFF",lambdaff_0,1.0d5)
 
      call closefile


** checking branching ratios are between 0 and 1
      if ((brrat_0 .lt. 0d0) .or. (brrat_0 .gt. 1d0)) then
         brrat_0 = 1d0
         write(*,*)" "
         write(*,*)"WARNING! BRRAT is outside the allowed range"
         write(*,*)"Using default BRRAT = 1d0 instead"
      end if
      if ((brrattrippm_0 .lt. 0d0) .or. (brrattrippm_0 .gt. 1d0)) then
         brrattrippm_0 = 1d0
         write(*,*)" "
         write(*,*)"WARNING! BRRATTRIPPM is outside the allowed range"
         write(*,*)"Using default BRRATTRIPPM = 1d0 instead"
      end if
      if ((brrattripn_0 .lt. 0d0) .or. (brrattripn_0 .gt. 1d0)) then
         brrattripn_0 = 1d0
         write(*,*)" "
         write(*,*)"WARNING! BRRATTRIPN is outside the allowed range"
         write(*,*)"Using default BRRATTRIPN = 1d0 instead"
      end if
         

      f1=f1_0
      f2=f2_0
      f3=f3_0
      f4=f4_0
      f5=f5_0
      f9=f9_0
      f6=f6_0
      f7=f7_0
      lambda=lambda_0
      sp2mass=sp2mass_0
      msp2trippm=msp2trippm_0
      msp2tripn=msp2tripn_0
      brrat=brrat_0
      brrattrippm=brrattrippm_0
      brrattripn=brrattripn_0
      nff=nff_0
      lambdaff=lambdaff_0


      end
