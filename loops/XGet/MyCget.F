      subroutine MyCget(p1sq,p2sq,s12,musq,C0123R,C0123I,Cij)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 6/4/2010
c Returns C0 and Cij     
*****
      implicit none
      Complex*16 B0fin,C0fin
      External B0fin,C0fin
      Real*8 p1sq,p2sq,s12
      Complex*16 B012,B013,B023
      Real*8 B012R,B013R,B023R
      Real*8 B012I,B013I,B023I
      Complex*16 C0123
      Real*8 C0123R
      Real*8 C0123I
      Real*8 Cij123R(4,2)
      Real*8 Cij123I(4,2)
      complex*16 Cij(30,9)
      real*8 musq


       B012=B0fin(p1sq,musq)    
c       print*, 'B012', B012
     
       B013=B0fin(s12,musq)  
c        print*, 'B013', B013

       B023=B0fin(p2sq,musq) 
c       print*, 'B023', B023


       B012R=Dble(B012)          
       B013R=Dble(B013) 
       B023R=Dble(B023 ) 
       B012I=DIMAG(B012)          
       B013I=DIMAG(B013)
       B023I=DIMAG(B023) 


       C0123=C0fin(p1sq,p2sq,s12,musq)


       call tens_red3_new_Re_Com(p1sq,p2sq,s12,B023,B013,B012,C0123,C0
     &   123R,C0123I,Cij123R,Cij123I) 
 
       call ten_red_LUdecomR5(p1sq,p2sq,s12,B023,B013,B012, 
     &                     Cij123r,Cij123I,Cij)

    
       return
       end
