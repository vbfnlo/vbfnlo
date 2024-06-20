      subroutine spin2width
       
c this subroutine computes the total decay widths of spin-2 particles:
c sp2width: width of spin-2 singlet particle
c sp2trippmwidth: width of charged spin-2 triplet particle
c sp2tripnwidth: width of neutral spin-2 triplet particle

c by Jessica Frank (frankj@particle.uni-karlsruhe.de)

      implicit none

#include "koppln.inc"
#include "global.inc"
#include "spin2coupl.inc"

c variables for spin2 width:
      double precision sp2widthaa, sp2widthaz, sp2widthzz,
     &     sp2widthww, sp2widthgg, sinw, cosw, gw1, v1, cos2w
      
      double precision sp2tripnwidthww, sp2tripnwidthzz,
     &     sp2tripnwidthaz, sp2tripnwidthaa
      double precision sp2trippmwidthwa, sp2trippmwidthwz
      
      double complex       gal(2), gad(2), gau(2), gwf(2),
     &     gzn(2), gzl(2), gzd(2), gzu(2)
      double precision     gw, gwwa, gwwz
      common /COUPL_GAUGE/ gal   , gad   , gau   , gwf   ,
     &     gzn   , gzl   , gzd   , gzu   ,
     &     gw, gwwa, gwwz
      double complex       gg(2)
      double precision     alpha_test,ee, sin2w_test,
     &     gfermi, alfas_test,g_test
      common /COUPL_BASIC/ gg, alpha_test,ee, sin2w_test, gfermi,
     &     alfas_test,g_test  
      
      
      cos2w = 1.0d0 - sin2w
      cosw = sqrt(cos2w)
      sinw = sqrt(sin2w)
      gw1 = ee / cosw
      v1 = 2.0d0*xmw/gw

c calculation of sp2width

      write(*,*)' '

      if ((sp2mass**2/4.0d0) .gt. (xmw**2)) then
         sp2widthww = ((24.0d0*f2**2*(sp2mass**4-3.0d0*sp2mass**2*xmw**2
     &        +6.0d0*xmw**4)+40.0d0*f2*f5*gw**2*v1**2*(sp2mass**2
     &        -xmw**2))/(12.0d0*LAMBDA**2)+f5**2*gw**4*v1**4*
     &        (sp2mass**4+12.0d0*sp2mass**2*xmw**2+56.0d0*xmw**4)/
     &        (96.0d0*LAMBDA**2*xmw**4))
     &        *sqrt((sp2mass**2/4.0d0-xmw**2))
     &        /(32.0d0*PI**2*sp2mass**2)*1.0d0/5.0d0*4.0d0*PI
      else
         sp2widthww = 0D0
      end if
      

      if ((sp2mass**2/4.0d0) .gt. (xmz**2)) then
         sp2widthzz = ((24.0d0*f2**2*cos2w**2*(sp2mass**4-3.0d0*sp2mass**2
     &        *xmz**2+6.0d0*xmz**4)+8.0d0*cos2w*f2*
     &        (6.0d0*f1*sin2w*(sp2mass**4-3.0d0*sp2mass**2
     &        *xmz**2+6.0d0*xmz**4)+5.0d0*f5*v1**2*(gw**2+gw1**2)
     &        *(sp2mass**2-xmz**2))+24.0d0*f1**2*sin2w**2*
     &        (sp2mass**4-3.0d0*sp2mass**2*xmz**2+6.0d0*xmz**4)+
     &        40.0d0*f1*f5*sin2w*v1**2*(gw**2+gw1**2)*
     &        (sp2mass**2-xmz**2))
     &        /(12.0d0*LAMBDA**2)+f5**2*v1**4*(gw**2+gw1**2)**2*
     &        (sp2mass**4+12.0d0*sp2mass**2*xmz**2+56.0d0*xmz**4)
     &        /(96.0d0*LAMBDA**2*xmz**4))
     &        *1.0d0/10.0d0
     &        *sqrt((sp2mass**2/4.0d0-xmz**2))
     &        /(32.0d0*PI**2*sp2mass**2)*4.0d0*PI
      else 
         sp2widthzz = 0D0
      end if

      sp2widthaa = 2.0d0*sp2mass**4*(cos2w*f1+f2*sin2w)**2/(LAMBDA**2)
     &     /(64.0d0*PI**2*sp2mass)*1.0d0/10.0d0*4.0d0*PI
      
      sp2widthaz = (cos2w*sin2w*(f1-f2)**2*(sp2mass**2-xmz**2)**2
     &     *(6.0d0*sp2mass**4+3.0d0*sp2mass**2*xmz**2+xmz**4
     &     ))/(3.0d0*LAMBDA**2*sp2mass**4)*1.0d0/5.0d0*4.0d0*PI
     &     *(sp2mass**2-xmz**2)/(64.0d0*PI**2*sp2mass**3)

c new: gluon-gluon contribution (coupling: f9)
c symmetry/colorfactor: *1/2 (symm ident gluons) *8 (colour delta)
c multiplication with 1.7: account for higher-order QCD correction 
c factor to partial width 

      sp2widthgg = 2.0d0*sp2mass**4*f9**2/(LAMBDA**2)
     &     /(64.0d0*PI**2*sp2mass)*1.0d0/10.0d0*8.0d0*4.0d0*PI
     &     *1.7d0

      
      if ((f1.EQ.0.0d0) .and. (f2.EQ.0.0d0) .and. (f3.EQ.0.0d0) .and.
     &     (f4.EQ.0.0d0) .and. (f5.EQ.0.0d0) .and. (f9.EQ.0.0d0)) then
       write(*,*) "sp2width is set to 1"
       write(*,*) "(it is 0 because the spin-2 singlet couplings are 0)"
         sp2width = 1.0d0

      else

         if (((sp2mass**2/4.0d0-xmw**2).LT.0.0d0) .and.
     &        ((sp2mass**2/4.0d0-xmz**2).LT.0.0d0) .and.
     &        ((sp2mass**2-xmz**2).LT.0.0d0)) then
            sp2width = (sp2widthaa+sp2widthgg)*1.0d0/brrat
         elseif (((sp2mass**2/4.0d0-xmw**2).LT.0.0d0) .and.
     &           ((sp2mass**2/4.0d0-xmz**2).LT.0.0d0) .and.
     &           ((sp2mass**2-xmz**2).GE.0.0d0)) then
            sp2width = (sp2widthaa+sp2widthaz+sp2widthgg)*1.0d0/brrat
         elseif (((sp2mass**2/4.0d0-xmw**2).GE.0.0d0) .and. 
     &           ((sp2mass**2-xmz**2).GE.0.0d0) .and.
     &           ((sp2mass**2/4.0d0-xmz**2).LT.0.0d0)) then
            sp2width = (sp2widthaa+sp2widthaz+sp2widthww+sp2widthgg)
     &           *1.0d0/brrat
         else
            sp2width = (sp2widthaa+sp2widthaz+sp2widthww+sp2widthzz
     &                  +sp2widthgg)*1.0d0/brrat
         endif
      
      endif


c calculation of sp2tripnwidth

      if ((msp2tripn**2/4.0d0) .gt. (xmw**2)) then
         sp2tripnwidthww = (f6**2*gw**4*v1**4*
     &        (msp2tripn**4+12.0d0*msp2tripn**2*xmw**2+56.0d0*xmw**4)/
     &        (384.0d0*LAMBDA**2*xmw**4))
     &        *sqrt((msp2tripn**2/4.0d0-xmw**2))
     &        /(32.0d0*PI**2*msp2tripn**2)*1.0d0/5.0d0*4.0d0*PI
      else
         sp2tripnwidthww = 0D0
      end if
      

      if ((msp2tripn**2/4.0d0) .gt. (xmz**2)) then
         sp2tripnwidthzz = ((768.0d0*f7**2*cos2w*sin2w*xmz**4*
     &        (msp2tripn**4-3.0d0*msp2tripn**2
     &        *xmz**2+6.0d0*xmz**4)+640.0d0*cosw*f6*f7*xmz**4
     &        *sinw*v1**2*(gw**2+gw1**2)
     &        *(msp2tripn**2-xmz**2)+f6**2*v1**4*(gw**2+gw1**2)**2*
     &        (msp2tripn**4+12.0d0*msp2tripn**2*xmz**2+56.0d0*xmz**4))
     &        /(384.0d0*LAMBDA**2*xmz**4))
     &        *1.0d0/10.0d0
     &        *sqrt((msp2tripn**2/4.0d0-xmz**2))
     &        /(32.0d0*PI**2*msp2tripn**2)*4.0d0*PI
      else
         sp2tripnwidthzz = 0D0
      end if

      sp2tripnwidthaa = 2.0d0*msp2tripn**4*cos2w*sin2w*f7**2
     &     /(LAMBDA**2)/(64.0d0*PI**2*msp2tripn)
     &     *1.0d0/10.0d0*4.0d0*PI
      
      sp2tripnwidthaz = ((cos2w-sin2w)**2*f7**2*
     &     (msp2tripn**2-xmz**2)**2
     &     *(6.0d0*msp2tripn**4+3.0d0*msp2tripn**2*xmz**2+xmz**4
     &     ))/(12.0d0*LAMBDA**2*msp2tripn**4)*1.0d0/5.0d0*4.0d0*PI
     &     *(msp2tripn**2-xmz**2)/(64.0d0*PI**2*msp2tripn**3)
      
      if ((f6.EQ.0.0d0) .and. (f7.EQ.0.0d0)) then
       write(*,*) "sp2tripnwidth is set to 1"
       write(*,*) "(it is 0 because the spin-2 triplet couplings are 0)"
         sp2tripnwidth = 1.0d0

      else
         
         if (((msp2tripn**2/4.0d0-xmw**2).LT.0.0d0) .and.
     &        ((msp2tripn**2/4.0d0-xmz**2).LT.0.0d0) .and.
     &        ((msp2tripn**2-xmz**2).LT.0.0d0)) then
            sp2tripnwidth = sp2tripnwidthaa*1.0d0/brrattripn
         elseif (((msp2tripn**2/4.0d0-xmw**2).LT.0.0d0) .and.
     &           ((msp2tripn**2/4.0d0-xmz**2).LT.0.0d0) .and.
     &           ((msp2tripn**2-xmz**2).GE.0.0d0)) then
            sp2tripnwidth = (sp2tripnwidthaa+sp2tripnwidthaz)
     &           *1.0d0/brrattripn
         elseif (((msp2tripn**2/4.0d0-xmw**2).GE.0.0d0) .and. 
     &           ((msp2tripn**2-xmz**2).GE.0.0d0) .and.
     &           ((msp2tripn**2/4.0d0-xmz**2).LT.0.0d0)) then
            sp2tripnwidth = (sp2tripnwidthaa+sp2tripnwidthaz
     &           +sp2tripnwidthww)
     &           *1.0d0/brrattripn
         else
            sp2tripnwidth = (sp2tripnwidthaa+sp2tripnwidthaz
     &           +sp2tripnwidthww+sp2tripnwidthzz)
     &           *1.0d0/brrattripn
         endif
         
      endif

c calculation of sp2trippmwidth

      sp2trippmwidthwa = (cos2w*f7**2*(msp2trippm**2-xmw**2)**2
     &     *(6.0d0*msp2trippm**4+3.0d0*msp2trippm**2*xmw**2+xmw**4
     &     ))/(12.0d0*LAMBDA**2*msp2trippm**4)*1.0d0/5.0d0*4.0d0*PI
     &     *(msp2trippm**2-xmw**2)/(64.0d0*PI**2*msp2trippm**3)
      
      if (msp2trippm .gt. (xmw+xmz)) then
         sp2trippmwidthwz = ((msp2trippm**2*xmw**2*(msp2trippm**2*xmz**2
     &        *(13.0d0*f6**2*gw**2*v1**4*(gw**2+gw1**2)+256.0d0*f7**2
     &        *xmw**2*xmz**2*sin2w)+1.0d0/4.0d0*(msp2trippm**2-xmw**2
     &        +xmz**2)**2*(7.0d0*f6**2*gw**2*v1**4*(gw**2+gw1**2)-
     &        96.0d0*f7**2*xmw**2*xmz**2*sin2w))+1.0d0/4.0d0*
     &        (msp2trippm**2+xmw**2-xmz**2)**2*(msp2trippm**2*xmz**2
     &        *(7.0d0*f6**2*gw**2*v1**4*(gw**2+gw1**2)-96.0d0*f7**2*
     &        xmw**2*xmz**2*sin2w)+(msp2trippm**2-xmw**2+xmz**2)**2*
     &        (f6**2*gw**2*v1**4*(gw**2+gw1**2)+32.0d0*f7**2*xmw**2*
     &        xmz**2*sin2w))+msp2trippm**2*(msp2trippm**2-xmw**2-
     &        xmz**2)*(1.0d0/4.0d0*(msp2trippm**2+xmw**2-xmz**2)*
     &        (msp2trippm**2-xmw**2+xmz**2)*(128.0d0*f7**2*xmw**2*xmz**2
     &        *sin2w-f6**2*gw**2*v1**4*(gw**2+gw1**2))+40.0d0*f6*f7*gw*
     &        msp2trippm**2*xmw**2*xmz**2*sinw*v1**2*sqrt((gw**2+gw1**2)
     &        ))+1.0d0/4.0d0*msp2trippm**4*(-msp2trippm**2+xmw**2
     &        +xmz**2)**2*(f6**2*gw**2*v1**4*(gw**2+gw1**2)+32.0d0*f7**2
     &        *xmw**2*xmz**2*sin2w)+40.0d0*f6*f7*gw*msp2trippm**2*xmw**2
     &        *xmz**2*sinw*v1**2*sqrt((gw**2+gw1**2))*(msp2trippm**2
     &        +xmw**2-xmz**2)*(msp2trippm**2-xmw**2+xmz**2))
     &        /(96.0d0*LAMBDA**2*msp2trippm**4*xmw**2*xmz**2))
     &        *1.0d0/5.0d0*4.0d0*PI*
     &        sqrt((msp2trippm**2-xmw**2-xmz**2)**2-4.0d0*xmw**2*xmz**2
     &        )/(64.0d0*PI**2*msp2trippm**3)
      else 
         sp2trippmwidthwz = 0D0
      end if


      if ((f6.EQ.0.0d0) .and. (f7.EQ.0.0d0)) then
       write(*,*) "sp2trippmwidth is set to 1"
       write(*,*) "(it is 0 because the spin-2 triplet couplings are 0)"
        sp2trippmwidth = 1.0d0

      else

         if ((msp2trippm**2-xmw**2).LT.0.0d0) then
            write(*,*) "charged spin-2 triplet particle cannot decay"
            write(*,*) "(it is too light)"
            sp2trippmwidth = 0.0d0
         elseif (((msp2trippm**2-xmw**2).GE.0.0d0) .and.
     &           ((msp2trippm-xmw-xmz).LT.0.0d0)) then
            sp2trippmwidth = sp2trippmwidthwa*1.0d0/brrattrippm
         else 
            sp2trippmwidth = (sp2trippmwidthwa+sp2trippmwidthwz)
     &           *1.0d0/brrattrippm
         endif
         
      endif


      write(*,*) "sp2width", sp2width
      write(*,*) "sp2widthww", sp2widthww
      write(*,*) "sp2widthzz", sp2widthzz
      write(*,*) "sp2widthaa", sp2widthaa
      write(*,*) "sp2widthaz", sp2widthaz
      write(*,*) "sp2widthgg", sp2widthgg
      
      write(*,*) "sp2tripnwidth", sp2tripnwidth
      write(*,*) "sp2tripnwidthww", sp2tripnwidthww
      write(*,*) "sp2tripnwidthzz", sp2tripnwidthzz
      write(*,*) "sp2tripnwidthaa", sp2tripnwidthaa
      write(*,*) "sp2tripnwidthaz", sp2tripnwidthaz
      
      write(*,*) "sp2trippmwidth", sp2trippmwidth
      write(*,*) "sp2trippmwidthwa", sp2trippmwidthwa
      write(*,*) "sp2trippmwidthwz", sp2trippmwidthwz
            
      RETURN
      END

