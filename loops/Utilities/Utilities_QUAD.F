
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


      subroutine print_Cij_QUAD(C0,Cij)
      implicit none
      complex*32 C0,Cij(4,2)
      integer indx_C(2)
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

      subroutine print_CijR_QUAD(C0,CijR,CijI)
      implicit none
      complex*32 C0
      real*16  CijR(4,2),CijI(4,2)
      integer indx_C(2)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_C(1)=2
      indx_C(2)=4
      write(*,'(A4,E44.36,E44.36,A1)') 'C0=(',C0,')'
      do i1=1,2
         do i2=1,indx_C(i1)
            write(*,'(A4,I1,A1,I1,A3,E44.36,E44.36,A1)') 
     -   'Cij(',i2,',',i1,')=(',CMPLX(CijR(i2,i1),CijI(i2,i1),16),')'
         enddo
      enddo
      end

      subroutine print_CijR_QUAD_D(C0,CijR,CijI,C0_QUAD,CijR_QUAD,CijI_QUAD)
      implicit none
      complex*32 C0_QUAD
      complex*16 C0
      real*16  CijR_QUAD(4,2),CijI_QUAD(4,2)
      real*8  CijR(4,2),CijI(4,2)
      integer indx_C(2)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_C(1)=2
      indx_C(2)=4
      write(*,'(A6,E30.22,E30.22,A1)') 'C0_D=(',C0_QUAD-C0,')'
      do i1=1,2
         do i2=1,indx_C(i1)
            write(*,'(A6,I1,A1,I1,A3,E24.16,E24.16,A1)') 
     -  'Cij_D(',i2,',',i1,')=(',DCMPLX(CijR_QUAD(i2,i1)
     -      -CijR(i2,i1),CijI_QUAD(i2,i1)-CijI(i2,i1)),')'
         enddo
      enddo
      end


      subroutine print_Dij_QUAD(D0,Dij)
      implicit none
      complex*32 D0,Dij(13,3)
      integer indx_D(3)
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
      end
     
      subroutine print_DijR_QUAD(D0,DijR,DijI)
      implicit none
      complex*32 D0
      real*16 DijR(13,3),DijI(13,3)
      integer indx_D(3)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_D(1)=13
      indx_D(2)=7
      indx_D(3)=13
      write(*,'(A4,E44.36,E44.36,A1)') 'D0=(',D0,')'
      do i1=1,3
         Print*, "********rank",i1
         do i2=1,indx_D(i1)
            write(*,'(A4,I2,A1,I1,A3,E44.36,E44.36,A1)')
     -     'Dij(',i2,',',i1,')=(',CMPLX(DijR(i2,i1),DijI(i2,i1),16),')'
         enddo
         enddo
      end
    
      subroutine print_DijR_QUAD_D(D0,DijR,DijI,D0_QUAD,DijR_QUAD,DijI_QUAD)
      implicit none
      complex*32 D0_QUAD
      complex*16 D0
      real*16  DijR_QUAD(13,3),DijI_QUAD(13,3)
      real*8  DijR(13,3),DijI(13,3)
      integer indx_D(3)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_D(1)=13
      indx_D(2)=7
      indx_D(3)=13
      write(*,'(A6,E30.22,E30.22,A1)') 'D0_D=(',D0_QUAD-D0,')'
      do i1=1,3
         do i2=1,indx_D(i1)
            write(*,'(A6,I2,A1,I1,A3,E24.16,E24.16,A1)') 
     -  'Dij_D(',i2,',',i1,')=(',DCMPLX(DijR_QUAD(i2,i1)
     -      -DijR(i2,i1),DijI_QUAD(i2,i1)-DijI(i2,i1)),')'
         enddo
      enddo
      end

      subroutine print_Eij_QUAD(E0,Eij)
      implicit none
      complex*32 E0,Eij(46,4)
      integer indx_E(4)
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
      

      subroutine print_EijR_QUAD_D(E0,EijR,EijI,E0_QUAD,EijR_QUAD,EijI_QUAD)
      implicit none
      complex*32 E0_QUAD
      complex*16 E0
      real*16  EijR_QUAD(46,4),EijI_QUAD(46,4)
      real*8  EijR(46,4),EijI(46,4)
      integer indx_E(4)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_E(1)=4
      indx_E(2)=11
      indx_E(3)=24
      indx_E(4)=46
      write(*,'(A6,E30.22,E30.22,A1)') 'E0_D=(',E0_QUAD-E0,')'
      do i1=1,4
         do i2=1,indx_E(i1)
            write(*,'(A6,I2,A1,I1,A3,E24.16,E24.16,A1)') 
     -  'Eij_D(',i2,',',i1,')=(',DCMPLX(EijR_QUAD(i2,i1)
     -      -EijR(i2,i1),EijI_QUAD(i2,i1)-EijI(i2,i1)),')'
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


      subroutine print_EijR_QUAD(E0,EijR,EijI)
      implicit none
      complex*32 E0
      real*16 EijR(46,4),EijI(46,4)
      integer indx_E(4)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_E(1)=4
      indx_E(2)=11
      indx_E(3)=24
      indx_E(4)=46
      write(*,'(A4,E44.36,E44.36,A1)') 'E0=(',E0,')'
      do i1=1,4
        print*, "******rank", i1
         do i2=1,indx_E(i1)
             write(*,'(A4,I2,A1,I1,A3,E44.36,E44.36,A1)')
     -     'Eij(',i2,',',i1,')=(',CMPLX(EijR(i2,i1),EijI(i2,i1),16),')'
         enddo
      enddo
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
      
      subroutine print_FijR_QUAD(F0,FijR,FijI)
      implicit none
      complex*32 F0
      real*16 FijR(166,5),FijI(166,5)
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
      

      subroutine print_FijR_QUAD_D(F0,FijR,FijI,F0_QUAD,FijR_QUAD,FijI_QUAD)
      implicit none
      complex*32 F0_QUAD
      complex*16 F0
      real*16  FijR_QUAD(166,4),FijI_QUAD(166,4)
      real*8  FijR(46,4),FijI(46,4)
      integer indx_F(5)
c      common/indexes/indx_C,indx_D,indx_E
      integer i1,i2
      indx_F(1)=5
      indx_F(2)=16
      indx_F(3)=40
      indx_F(4)=86
      indx_F(5)=166
      write(*,'(A6,E30.22,E30.22,A1)') 'F0_D=(',F0_QUAD-F0,')'
      do i1=1,5
         do i2=1,indx_F(i1)
            write(*,'(A6,I2,A1,I1,A3,E24.16,E24.16,A1)') 
     -  'Fij_D(',i2,',',i1,')=(',DCMPLX(FijR_QUAD(i2,i1)
     -      -FijR(i2,i1),FijI_QUAD(i2,i1)-FijI(i2,i1)),')'
         enddo
      enddo
      end


      subroutine From_EijR_QUAD_to_Dbl(E0R,E0I,EijR,EijI,
     -                            E0R_QUAD,E0I_QUAD,EijR_QUAD,EijI_QUAD)
      implicit none
      real*16 E0R_QUAD,E0I_QUAD
      real*8 E0R,E0I
      real*16  EijR_QUAD(46,4),EijI_QUAD(46,4)
      real*8  EijR(46,4),EijI(46,4)
      integer indx_E(4)
      integer i1,i2
      indx_E(1)=4
      indx_E(2)=11
      indx_E(3)=24
      indx_E(4)=46
      E0R=Dble(E0R_QUAD)
      E0I=Dble(E0I_QUAD)
      do i1=1,4
         do i2=1,indx_E(i1)
           EijR(i2,i1)=Dble(EijR_QUAD(i2,i1))
           EijI(i2,i1)=Dble(EijI_QUAD(i2,i1))
         enddo
      enddo
      end


      subroutine From_DijR_QUAD_to_Dbl(D0R,D0I,DijR,DijI,
     -                      D0R_QUAD, D0I_QUAD,DijR_QUAD,DijI_QUAD)
      implicit none
      real*16 D0R_QUAD,D0I_QUAD
      real*8 D0R,D0I
      real*16  DijR_QUAD(13,3),DijI_QUAD(13,3)
      real*8  DijR(13,3),DijI(13,3)
      integer indx_D(3)

      integer i1,i2
      indx_D(1)=13
      indx_D(2)=7
      indx_D(3)=13
      D0R=Dble(D0R_QUAD)
      D0I=Dble(D0I_QUAD)
      do i1=1,3
         do i2=1,indx_D(i1)
           DijR(i2,i1)=Dble(DijR_QUAD(i2,i1))
           DijI(i2,i1)=Dble(DijI_QUAD(i2,i1))
         enddo
      enddo
      end


      subroutine From_CijR_QUAD_to_Dbl(C0R,C0I,CijR,CijI,
     -                         C0R_QUAD,C0I_QUAD,CijR_QUAD,CijI_QUAD)
      implicit none
      real*16 C0R_QUAD,C0I_QUAD
      real*8 C0R,C0I
      real*16  CijR_QUAD(4,2),CijI_QUAD(4,2)
      real*8  CijR(4,2),CijI(4,2)
      integer indx_C(2)
      integer i1,i2
      indx_C(1)=2
      indx_C(2)=4
      C0R=Dble(C0R_QUAD)
      C0I=Dble(C0I_QUAD)
      do i1=1,2
         do i2=1,indx_C(i1)
           CijR(i2,i1)=Dble(CijR_QUAD(i2,i1))
           CijI(i2,i1)=Dble(CijI_QUAD(i2,i1))
         enddo
      enddo
      end


      subroutine From_BijR_QUAD_to_Dbl(B0R,B0I,B0R_QUAD,B0I_QUAD)
      implicit none
      real*16 B0R_QUAD,B0I_QUAD
      real*8 B0R,B0I
      B0R=Dble(B0R_QUAD)
      B0I=Dble(B0I_QUAD)
      end


      subroutine From_FijR_QUAD_to_Dbl(F0R,F0I,FijR,FijI,
     -                            F0R_QUAD,F0I_QUAD,FijR_QUAD,FijI_QUAD)
      implicit none
      real*16 F0R_QUAD,F0I_QUAD
      real*8 F0R,F0I
      real*16  FijR_QUAD(166,5),FijI_QUAD(166,5)
      real*8  FijR(166,5),FijI(166,5)
      integer indx_F(5)
      integer i1,i2
      indx_F(1)=5
      indx_F(2)=16
      indx_F(3)=40
      indx_F(4)=86
      indx_F(5)=166

      F0R=Dble(F0R_QUAD)
      F0I=dble(F0I_QUAD)
      do i1=1,5
         do i2=1,indx_F(i1)
           FijR(i2,i1)=Dble(FijR_QUAD(i2,i1))
           FijI(i2,i1)=Dble(FijI_QUAD(i2,i1))
         enddo
      enddo
      end



