      subroutine MyCget_M(m0,p1sq,p2sq,s12,musq,C0123R,C0123I,Cij)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 6/4/2010
c Returns C0 and Cij     
*****
      implicit none
      Complex*16 B0finG,C0finG
      External B0finG,C0finG
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
      REAL*8 MUSQCP
      common/musqINv/musqcp
      real*8 m0

      MUSQCP=MUSQ

       B012=B0finG(m0,m0,p1sq,musq)    
c       print*, 'B012', B012
     
       B013=B0finG(m0,m0,s12,musq)  
c        print*, 'B013', B013

       B023=B0finG(m0,m0,p2sq,musq) 
c       print*, 'B023', B023


       B012R=Dble(B012)          
       B013R=Dble(B013) 
       B023R=Dble(B023 ) 
       B012I=DIMAG(B012)          
       B013I=DIMAG(B013)
       B023I=DIMAG(B023) 


       C0123=C0finG(m0,m0,m0,p1sq,p2sq,s12,musq)  

c       print*,"C0123_M",C0123
 
   
       call tens_red3_new_Re_Com_1M(m0,p1sq,p2sq,s12,B023,B013,B012,C0123,C0
     &   123R,C0123I,Cij123R,Cij123I) 
 
c       print*,'here'

       call ten_red_LUdecomR5_M(m0,p1sq,p2sq,s12,B023,B013,B012, 
     &                     Cij123r,Cij123I,Cij)

c           print*,'here3'
       return
       end
