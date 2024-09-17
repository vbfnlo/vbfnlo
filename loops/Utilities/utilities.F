! author: Francisco Campanario
C Subroutines:
c       subroutine print_Cij(C0,Cij)
c       subroutine print_BijR(B0,BijR,BijI)
c       subroutine ToZero_BijR(B0,BijR,BijI)
c       subroutine print_CijR(C0,CijR,CijI)
c       subroutine ToZero_CijR(C0,CijR,CijI)
c       subroutine ToDen_CijR(CijR,CijI,CijDeR,CijDeI)
c       subroutine  ToDen_DijR(DijR,DijI,DijDeR,DijDeI)
c       subroutine print_Dij(D0,Dij)
c       subroutine print_DijR(D0,DijR,DijI)
c       subroutine rename_DijR(D0,DijR,DijI,DijR1,DijI1)
c       subroutine ratio_DijR(D0,DijR,DijI,DijR1,DijI1)
c       subroutine ToZero_DijR(D0,DijR,DijI)
c       subroutine print_Eij(E0,Eij)
c       subroutine print_EijR(E0,EijR,EijI)
c       subroutine ToZero_EijR(E0,EijR,EijI)
c       subroutine print_FijR(F0,FijR,FijI)
c 
     

      subroutine print_Cij(C0,Cij)
      implicit none
      complex * 16 C0,Cij(4,2)
      integer indx_C(2),indx_D(3),indx_E(4)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_C(1)=2
      indx_C(2)=4
      write(*,'(A4,E30.22,E30.22)') 'C0 =',C0
      do i1=1,2
         do i2=1,indx_C(i1)
            write(*,*) 'Cij(',i2,i1,') =',Cij(i2,i1)
         enddo
      enddo
      end

      subroutine print_BijR(B0,BijR,BijI)
      implicit none
      complex * 16 B0
      real*8  BijR,BijI
      integer indx_B(2)
      integer i1,i2
      indx_B(1)=1
      write(*,'(A4,E30.22,E30.22)') 'B0 =',B0

         do i2=1,indx_B(1)
            write(*,*) 'Bij(',1,1,') =',DCMPLX(BijR,BijI)
            write(*,*) '-B0/2(',1,1,') =',-B0/2d0
         enddo
         return
      end

      subroutine print_BijRRank2(B0,BijR,BijI)
      implicit none
      complex * 16 B0
      real*8  BijR(2,2),BijI(2,2)
      integer indx_B(2)
      integer i1,i2
      indx_B(1)=1
      indx_B(2)=2

      write(*,'(A4,E30.22,E30.22)') 'B0 =',B0
      do i1=1,2
         do i2=1,indx_B(i1)
            write(*,*) 'Bij(',i2,i1,') =',DCMPLX(BijR(i2,i1),BijI(i2,i1))
         enddo
      enddo
         return
      end




      subroutine ToZero_BijR(B0,BijR,BijI)
      implicit none
      complex * 16 B0
      real*8  BijR,BijI
      integer indx_B(2)
      integer i1,i2
      indx_B(1)=1
      B0=0d0
      BijR=0d0
      BijI=0d0
      return
      end

      subroutine print_CijR(C0,CijR,CijI)
      implicit none
      complex*16 C0
      real*8  CijR(4,2),CijI(4,2)
      integer indx_C(2)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_C(1)=2
      indx_C(2)=4
      write(*,'(A4,E34.16,E34.16,A1)') 'C0=(',C0,')'
      do i1=1,2
         do i2=1,indx_C(i1)
            write(*,'(A4,I1,A1,I1,A3,E34.16,E34.16,A1)') 
     -   'Cij(',i2,',',i1,')=(',CMPLX(CijR(i2,i1),CijI(i2,i1),16),')'
         enddo
      enddo
      end

!      subroutine print_CijR(C0,CijR,CijI)
!      implicit none
!      complex * 16 C0
!      real*8  CijR(4,2),CijI(4,2)
!      integer indx_C(2),indx_D(3),indx_E(4)
!c      common/indexes/indx_C,indx_D,indx_E
!      integer i1,i2
!      indx_C(1)=2
!      indx_C(2)=4
!      write(*,'(A4,E30.22,E30.22)') 'C0 =',C0
!      do i1=1,2
!         do i2=1,indx_C(i1)
!            write(*,*) 'Cij(',i2,i1,') =',DCMPLX(CijR(i2,i1),CijI(i2,i1))
!         enddo
!      enddo
!      return
!      end



      subroutine print_CijR1(Cij,C0,CijR,CijI)
      implicit none
      character*6 Cij
      complex * 16 C0
      real*8  CijR(4,2),CijI(4,2)
      integer indx_C(2),indx_D(3),indx_E(4)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_C(1)=2
      indx_C(2)=4
      write(*,'(A6,E30.22,E30.22)') 'C0 =',C0
      do i1=1,2
         do i2=1,indx_C(i1)
      write(*,'(A7,2I2,A2,E30.22,E30.22)') Cij//'(',i2,i1,')=',DCMPLX(CijR(i2,i1),CijI(i2,i1))
         enddo
      enddo
      return
      end


      subroutine print_CijR1_Q(Cij,C0,CijR,CijI)
      implicit none
      character*6 Cij
      complex * 32 C0
      real*16  CijR(4,2),CijI(4,2)
      integer indx_C(2),indx_D(3),indx_E(4)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_C(1)=2
      indx_C(2)=4
      write(*,'(A6,E30.22,E30.22)') 'C0 =',C0
      do i1=1,2
         do i2=1,indx_C(i1)
      write(*,'(A7,2I2,A2,E30.22,E30.22)') Cij//'(',i2,i1,')=',DCMPLX(CijR(i2,i1),CijI(i2,i1))
         enddo
      enddo
      return
      end




      subroutine print_CijR_R4(C0,CijR,CijI)
      implicit none
      complex * 16 C0
      real*8  CijR(9,4),CijI(9,4)
      integer indx_C(4),indx_D(3),indx_E(4)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_C(1)=2
      indx_C(2)=4
      indx_C(3)=6
      indx_C(4)=9
      write(*,'(A4,E30.22,E30.22)') 'C0 =',C0
      do i1=1,4
         do i2=1,indx_C(i1)
            write(*,*) 'Cij(',i2,i1,') =',DCMPLX(CijR(i2,i1),CijI(i2,i1))
         enddo
      enddo
      return
      end




      subroutine print_CijRRank3(C0,CijR,CijI)
      implicit none
      complex * 16 C0
      real*8  CijR(6,3),CijI(6,3)
      integer indx_C(3),indx_D(3),indx_E(4)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_C(1)=2
      indx_C(2)=4
      indx_C(3)=6
      write(*,'(A4,E30.22,E30.22)') 'C0 =',C0
      do i1=1,3
         do i2=1,indx_C(i1)
            write(*,*) 'Cij(',i2,i1,') =',DCMPLX(CijR(i2,i1),CijI(i2,i1))
         enddo
      enddo
      return
      end










      subroutine ToZero_CijR(C0,CijR,CijI)
      implicit none
      complex * 16 C0
      real*8  CijR(4,2),CijI(4,2)
      integer indx_C(2),indx_D(3),indx_E(4)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_C(1)=2
      indx_C(2)=4
      C0=0d0
      do i1=1,2
         do i2=1,indx_C(i1)
            CijR(i2,i1)=0d0
            CijI(i2,i1)=0d0
         enddo
      enddo
      return
      end





      subroutine ToDen_CijR(CijR,CijI,CijDeR,CijDeI)
      implicit none
      complex * 16 C0
      real*8  CijR(4,2),CijI(4,2)
      real*8  CijDeR(4,2),CijDeI(4,2) 
      integer indx_C(2),indx_D(3),indx_E(4)
c      common/indexes/indx_C,indx_D,indx_E
c      integer i1,i2
   

       CijDeR(1,1)=CijR(1,1)-CijR(2,1)
       CijDeR(2,1)=CijR(2,1)


       CijDeR(1,2)=CijR(1,2)-2d0*CijR(2,2)+CijR(3,2)
       CijDeR(2,2)=CijR(2,2)-CijR(3,2)
       CijDeR(3,2)=CijR(3,2)
       CijDeR(4,2)=CijR(4,2)
      
       CijDeI(1,1)=CijI(1,1)-CijI(2,1)
       CijDeI(2,1)=CijI(2,1)


       CijDeI(1,2)=CijI(1,2)-2d0*CijI(2,2)+CijI(3,2)
       CijDeI(2,2)=CijI(2,2)-CijI(3,2)
       CijDeI(3,2)=CijI(3,2)
       CijDeI(4,2)=CijI(4,2)
      end



      subroutine  ToDen_DijR(DijR,DijI,DijDeR,DijDeI)
      implicit none
      complex * 16 D0
      Real*8 DijR(13,3),DijI(13,3)
      Real*8 DijDeR(13,3),DijDeI(13,3)
      integer indx_C(2),indx_D(3),indx_E(4)
c      common/indexes/indx_C,indx_D,indx_E

c One Index
      DijDeR(1,1)=DijR(1,1)-DijR(2,1)
      DijDeR(2,1)=DijR(2,1)-DijR(3,1)
      DijDeR(3,1)=DijR(3,1)
c Two index
      DijDeR(1,2)=DijR(1,2)-2d0*DijR(4,2)+DijR(2,2)
      DijDeR(2,2)=DijR(2,2)-2d0*DijR(6,2)+DijR(3,2)
      DijDeR(3,2)=DijR(3,2)
      DijDeR(4,2)=DijR(4,2)-DijR(2,2)-DijR(5,2)+DijR(6,2)
      DijDeR(5,2)=DijR(5,2)-DijR(6,2)
      DijDeR(6,2)=DijR(6,2)-DijR(3,2)
      DijDeR(7,2)=DijR(7,2)
c Three indices
c 111  
      DijDeR(1,3)=DijR(1,3)+3d0*(DijR(6,3)-DijR(4,3))-DijR(2,3)
c222
      DijDeR(2,3)=DijR(2,3)+3d0*(DijR(9,3)-DijR(8,3))-DijR(3,3)
c333
      DijDeR(3,3)=DijR(3,3)
c112
      DijDeR(4,3)=DijR(4,3)+2d0*(DijR(10,3)-DijR(6,3))+DijR(2,3)
     \-DijR(5,3)-DijR(8,3)
c113
      DijDeR(5,3)=DijR(5,3)-2d0*DijR(10,3)+DijR(8,3)
c221
      DijDeR(6,3)=DijR(6,3)+2d0*(DijR(8,3)-DijR(10,3))-DijR(2,3)
     \+DijR(7,3)-DijR(9,3)
c331
      DijDeR(7,3)=DijR(7,3)-DijR(9,3)
c223
      DijDeR(8,3)=DijR(8,3)-2d0*DijR(9,3)+DijR(3,3)
c332
      DijDeR(9,3)=DijR(9,3)-DijR(3,3)
c123
       DijDeR(10,3)=DijR(10,3)-DijR(8,3)-DijR(7,3)+DijR(9,3)
c001
       DijDeR(11,3)=DijR(11,3)-DijR(12,3) 
c002
       DijDeR(12,3)=DijR(12,3)-DijR(13,3) 
c003
       DijDeR(13,3)=DijR(13,3) 



c One Index
      DijDeI(1,1)=DijI(1,1)-DijI(2,1)
      DijDeI(2,1)=DijI(2,1)-DijI(3,1)
      DijDeI(3,1)=DijI(3,1)
c Two index
      DijDeI(1,2)=DijI(1,2)-2d0*DijI(4,2)+DijI(2,2)
      DijDeI(2,2)=DijI(2,2)-2d0*DijI(6,2)+DijI(3,2)
      DijDeI(3,2)=DijI(3,2)
      DijDeI(4,2)=DijI(4,2)-DijI(2,2)-DijI(5,2)+DijI(6,2)
      DijDeI(5,2)=DijI(5,2)-DijI(6,2)
      DijDeI(6,2)=DijI(6,2)-DijI(3,2)
      DijDeI(7,2)=DijI(7,2)
c Three indices
c 111  
      DijDeI(1,3)=DijI(1,3)+3d0*(DijI(6,3)-DijI(4,3))-DijI(2,3)
c222
      DijDeI(2,3)=DijI(2,3)+3d0*(DijI(9,3)-DijI(8,3))-DijI(3,3)
c333
      DijDeI(3,3)=DijI(3,3)
c112
      DijDeI(4,3)=DijI(4,3)+2d0*(DijI(10,3)-DijI(6,3))+DijI(2,3)
     \-DijI(5,3)-DijI(8,3)
c113
      DijDeI(5,3)=DijI(5,3)-2d0*DijI(10,3)+DijI(8,3)
c221
      DijDeI(6,3)=DijI(6,3)+2d0*(DijI(8,3)-DijI(10,3))-DijI(2,3)
     \+DijI(7,3)-DijI(9,3)
c331
      DijDeI(7,3)=DijI(7,3)-DijI(9,3)
c223
      DijDeI(8,3)=DijI(8,3)-2d0*DijI(9,3)+DijI(3,3)
c332
      DijDeI(9,3)=DijI(9,3)-DijI(3,3)
c123
       DijDeI(10,3)=DijI(10,3)-DijI(8,3)-DijI(7,3)+DijI(9,3)
c001
       DijDeI(11,3)=DijI(11,3)-DijI(12,3) 
c002
       DijDeI(12,3)=DijI(12,3)-DijI(13,3) 
c003
       DijDeI(13,3)=DijI(13,3) 














      end
    







      subroutine print_Dij(D0,Dij)
      implicit none
      complex * 16 D0,Dij(13,3)
      integer indx_C(2),indx_D(3),indx_E(4)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_D(1)=13
      indx_D(1)=7
      indx_D(1)=13
      write(*,'(A4,E30.22,E30.22)') 'D0 =',D0
      do i1=1,3
         do i2=1,indx_D(i1)
            write(*,*) 'Dij(',i2,i1,') =',Dij(i2,i1)
         enddo
      enddo
      return
      end
     
!FC      subroutine print_DijR(D0,DijR,DijI)
!FC      implicit none
!FC      complex * 16 D0
!FC      Real*8 DijR(13,3),DijI(13,3)
!FC      integer indx_C(2),indx_D(3),indx_E(4)
!FCc      common/indexes/indx_C,indx_D,indx_E
!FC      integer i1,i2
!FC      indx_D(1)=13
!FC      indx_D(2)=7
!FC      indx_D(3)=13
!FC      write(*,'(a4,E30.22,E30.22)') 'D0 =',D0
!FC      do i1=1,3
!FC         Print*, "********rank",i1
!FC         do i2=1,indx_D(i1)
!FC            write(*,*) 'Dij(',i2,i1,') =',DCMPLX(DijR(i2,i1),DijI(i2,i1))
!FC         enddo
!FC         enddo
!FC         return
!FC      end

       subroutine print_DijR(D0,DijR,DijI)
      implicit none
      complex*16 D0
      real*8 DijR(13,3),DijI(13,3)
      integer indx_D(3)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_D(1)=13
      indx_D(2)=7
      indx_D(3)=13
      write(*,'(A4,E34.16,E434.16,A1)') 'D0=(',D0,')'
      do i1=1,3
         Print*, "********rank",i1
         do i2=1,indx_D(i1)
            write(*,'(A4,I2,A1,I1,A3,E34.16,E34.36,A1)')
     -     'Dij(',i2,',',i1,')=(',CMPLX(DijR(i2,i1),DijI(i2,i1),16),')'
         enddo
         enddo
      end




    
      subroutine print_DijRRank4(D0,DijR,DijI)
      implicit none
      complex * 16 D0
      Real*8 DijR(22,4),DijI(22,4)
      integer indx_C(2),indx_D(4),indx_E(4)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_D(1)=13
      indx_D(2)=7
      indx_D(3)=13
      indx_D(4)=22
      write(*,'(a4,E30.22,E30.22)') 'D0 =',D0
      do i1=1,4
         Print*, "********rank",i1
         do i2=1,indx_D(i1)
            write(*,*) 'Dij(',i2,i1,') =',DCMPLX(DijR(i2,i1),DijI(i2,i1))
         enddo
         enddo
         return
      end




      subroutine print_DijR_R5(D0,DijR,DijI)
      implicit none
      complex * 16 D0
      Real*8 DijR(34,5),DijI(34,5)
      integer indx_C(2),indx_D(5),indx_E(4)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_D(1)=13
      indx_D(2)=7
      indx_D(3)=13
      indx_D(4)=22
      indx_D(5)=34
      write(*,'(a4,E30.22,E30.22)') 'D0 =',D0
      do i1=1,5
         Print*, "********rank",i1
         do i2=1,indx_D(i1)
            write(*,*) 'Dij(',i2,i1,') =',DCMPLX(DijR(i2,i1),DijI(i2,i1))
         enddo
         enddo
         return
      end




      subroutine rename_DijR(D0,DijR,DijI,DijR1,DijI1)
      implicit none
      complex * 16 D0
      Real*8 DijR(13,3),DijI(13,3),DijR1(13,3),DijI1(13,3)
      integer indx_C(2),indx_D(3),indx_E(4)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_D(1)=13
      indx_D(2)=7
      indx_D(3)=13
  
      do i1=1,3
  
         do i2=1,indx_D(i1)
            DijR1(i2,i1)=DijR(i2,i1)
            DijI1(i2,i1)=DijI(i2,i1)
         enddo
         enddo
         return
      end

      subroutine ratio_DijR(D0,DijR,DijI,DijR1,DijI1)
      implicit none
      complex * 16 D0
      Real*8 DijR(13,3),DijI(13,3),DijR1(13,3),DijI1(13,3)
      integer indx_C(2),indx_D(3),indx_E(4)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_D(1)=13
      indx_D(2)=7
      indx_D(3)=13
  
      do i1=1,3
         Print*, "********rank",i1
         do i2=1,indx_D(i1)
       if (abs(DijR(i2,i1)).gt.1d-16) then
           if((DijR1(i2,i1)/DijR(i2,i1) -1d0).gt.1d-12) then
           write(*,*) 'DijR(',i2,i1,') =',DijR1(i2,i1)/DijR(i2,i1) -1d0
           endif
       else
           if((DijR1(i2,i1)-DijR(i2,i1)).gt.1d-12) then
           write(*,*) 'DijR(',i2,i1,') =',DijR1(i2,i1)- DijR(i2,i1) 
           endif
       endif

       if (abs(DijI(i2,i1)).gt.1d-16) then

          if((DijI1(i2,i1)/DijI(i2,i1) -1d0).gt.1d-12) then
           write(*,*) 'DijI(',i2,i1,') =',DijI1(i2,i1)/DijI(i2,i1) -1d0
           endif
        else
           if((DijI1(i2,i1)-DijI(i2,i1)).gt.1d-12) then
          write(*,*) 'DijI(',i2,i1,') =',DijI1(i2,i1) - DijI(i2,i1) 
            endif
        endif
         enddo
         enddo
         return
      end











      subroutine ToZero_DijR(D0,DijR,DijI)
      implicit none
      complex * 16 D0
      Real*8 DijR(13,3),DijI(13,3)
      integer indx_C(2),indx_D(3),indx_E(4)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_D(1)=13
      indx_D(2)=7
      indx_D(3)=13
      D0=0d0
      do i1=1,3
         Print*, "********rank",i1
         do i2=1,indx_D(i1)
            DijR(i2,i1)=0d0
            DijI(i2,i1)=0d0
         enddo
         enddo
         return
      end
    












      subroutine print_Eij(E0,Eij)
      implicit none
      complex * 16 E0,Eij(46,4)
      integer indx_C(2),indx_D(3),indx_E(4)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_E(1)=4
      indx_E(2)=11
      indx_E(3)=23
      indx_E(4)=46
      write(*,*) 'E0 =',E0
      do i1=1,4
         do i2=1,indx_E(i1)
            write(*,*) 'Eij(',i2,i1,') =',Eij(i2,i1)
         enddo
      enddo
      end
      

      subroutine print_EijR(E0,EijR,EijI)
      implicit none
      complex*16 E0
      real*8 EijR(46,4),EijI(46,4)
      integer indx_E(4)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_E(1)=4
      indx_E(2)=11
      indx_E(3)=24
      indx_E(4)=46
      write(*,'(A4,E34.16,E34.16,A1)') 'E0=(',E0,')'
      do i1=1,4
        print*, "******rank", i1
         do i2=1,indx_E(i1)
             write(*,'(A4,I2,A1,I1,A3,E34.16,E34.16,A1)')
     -     'Eij(',i2,',',i1,')=(',CMPLX(EijR(i2,i1),EijI(i2,i1),16),')'
         enddo
      enddo
      end

!      subroutine print_EijR(E0,EijR,EijI)
!      implicit none
!      complex * 16 E0
!      real*8 EijR(46,4),EijI(46,4)
!      integer indx_C(2),indx_D(3),indx_E(4)
!c      common/indexes/indx_C,indx_D,indx_E
!      integer i1,i2
!      indx_E(1)=4
!      indx_E(2)=11
!      indx_E(3)=23
!      indx_E(4)=46
!      write(*,*) 'E0 =',E0
!      do i1=1,4
!        print*, "******rank", i1
!         do i2=1,indx_E(i1)
!            write(*,*) 'Eij(',i2,i1,') =',DCMPLX(EijR(i2,i1),EijI(i2,i1))
!         enddo
!      enddo
!      return
!      end




      subroutine ToZero_EijR(E0,EijR,EijI)
      implicit none
      complex * 16 E0
      real*8 EijR(46,4),EijI(46,4)
      integer indx_C(2),indx_D(3),indx_E(4)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_E(1)=4
      indx_E(2)=11
      indx_E(3)=23
      indx_E(4)=46
      E0=0d0
      do i1=1,4
        print*, "******rank", i1
         do i2=1,indx_E(i1)
            EijR(i2,i1)=0d0
            EijI(i2,i1)=0d0
         enddo
      enddo
      return
      end



      subroutine print_FijR(F0,FijR,FijI)
      implicit none
      complex*16 F0
      real*8 FijR(166,5),FijI(166,5)
      integer indx_F(5)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_F(1)=5
      indx_F(2)=16
      indx_F(3)=40
      indx_F(4)=86
      indx_F(5)=166
      write(*,*) 'F0 =',F0
      do i1=1,5
        print*, "******rank", i1
         do i2=1,indx_F(i1)
            write(*,*) 'Fij(',i2,i1,') =',CMPLX(FijR(i2,i1),FijI(i2,i1),16)
         enddo
      enddo
      end

      
!      subroutine print_FijR(F0,FijR,FijI)
!      implicit none
!      complex * 16 F0
!      real*8 FijR(166,5),FijI(166,5)
!      integer indx_F(5)
!c      common/indexes/indx_C,indx_D,indx_E
!      integer i1,i2
!      indx_F(1)=5
!      indx_F(2)=16
!      indx_F(3)=40
!      indx_F(4)=86
!      indx_F(5)=166
!      write(*,*) 'F0 =',F0
!      do i1=1,5
!        print*, "******rank", i1
!         do i2=1,indx_F(i1)
!            write(*,*) 'Fij(',i2,i1,') =',DCMPLX(FijR(i2,i1),FijI(i2,i1))
!         enddo
!      enddo
!      end
      
