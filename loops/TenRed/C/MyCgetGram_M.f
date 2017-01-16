       subroutine MyCgetGram_M(m0,p1sq,p2sq,s12,musq,C0123R,C0123I,Cij)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 6/4/2010
c Returns C0 and Cij     
*****
       implicit none
       Complex*16 B0finG,C0finG,C3i3e
       External B0finG,C0finG,C3i3e
       Real*8 p1sq,p2sq,s12
       Complex*16 B012,B013,B023
       Complex*16 C0123
       Real*8 C0123R
       Real*8 C0123I
       complex*16 Cij(30,9)
       real*8 tempjj,tempkl,IX,Zmax
       integer jjinit
       common/Decide/tempjj,tempkl,IX,Zmax,jjinit
       Save/Decide/

       Real*8 zz11,zz12,zz22
       complex*16 ten2g,test2g,ten2p1p2,test2p1p2
       complex*16 ten2p1p2_Dble,ten2g_Dble
       real*8 r1, r2r1, p1p2
       real*8 musq,C0R,C0I
       integer jj,k,l,order
       real*8 m0
       real*8 accuracyC(0:4,5),AccuracyD(0:5,4)

      Common/Accuracy/AccuracyC,AccuracyD
       Common/musqInv/musqcp
      real*8 musqcp
       musqcp=musq

       tempjj=1d18
       tempkl=1d18
c       jjinit=12

       B012=B0finG(m0,m0,p1sq,musq)    
c       print*, 'B012', B012

       B013=B0finG(m0,m0,s12,musq)  
c        print*, 'B013', B013
       B023=B0finG(m0,m0,p2sq,musq) 
c       print*, 'B023', B023

       C0123=C0finG(m0,m0,m0,p1sq,p2sq,s12,musq)  
c       C0123=C3i3e(p1sq,p2sq,s12,m0*m0,m0*m0,m0*m0,musq,0) 
       C0R=Dble(C0123)
       C0I=DImag(C0123)
c       print*,"C0123",C0123
c FC %%      call tens_red3_new_Re_Com(p1sq,p2sq,s12,B023,B013,B012, 
c FC %%     &                     C0123,C0r,C0I,Cijr,CijI)
  
       k=2
       l=2
       jj=2
       order=0

c FC %      call ten_red3_Gram(p1sq,p2sq,s12,B023,B013,B012,C0123,Cij,k,l,jj,order)

       call ten_red3_Gram2_M(m0,p1sq,p2sq,s12,musq,B023,B013,B012,C0123,Cij)

       C0123R=Dble(C0123)
       C0123I=DImag(C0123)
       
c FC %%c CALL Dble
c FC %%       call MyCget1(p1sq,p2sq,s12,musq,C0R2,C0I2,CijR2,CijI2)
c Checking accuracy  
        p1p2 = (s12 - p1sq - p2sq)*0.5d0
        r1 = p1sq
        r2r1 = s12 - r1
      
	zz11=2d0*p1sq
	zz12=2d0*p1p2
	zz22=2d0*p2sq


      test2g=2d0*(B023+ C0123*m0*m0)

      test2p1p2=-B013*(ZZ12+ZZ22)/2d0 +(r2r1)*B023
     1   -r1*(B012-B013-(r2r1)*C0123)
       
       ten2g=ZZ11*Cij(1,2)+ZZ22*Cij(2,2)+ 
     1  2d0*(ZZ12*Cij(3,2) + 4d0*Cij(4,2))
     2 -1d0
       
       ten2p1p2=ZZ11*(ZZ12*Cij(1,2)+ZZ22*Cij(3,2)) + 
     1  ZZ12*(ZZ22*Cij(2,2) + ZZ12*Cij(3,2)+2*Cij(4,2))
         
cFC %       ten2g_Dble=ZZ11*DCMPLX(CijR2(1,2),CijI2(1,2))+ZZ22*DCMPLX(CijR2(2,2),CijI2(2,2))+ 
cFC %     1  2d0*(ZZ12*DCMPLX(CijR2(3,2),CijI2(3,2)) + 4d0*DCMPLX(CijR2(4,2),CijI2(4,2)))
cFC %     2 -1d0
    
cFC %        ten2p1p2_Dble=ZZ11*(ZZ12*DCMPLX(CijR2(1,2),CijI2(1,2))+ZZ22*DCMPLX(CijR2(3,2),CijI2(3,2))) + 
cFC %     1  ZZ12*(ZZ22*DCMPLX(CijR2(2,2),CijI2(2,2)) + ZZ12*DCMPLX(CijR2(3,2),CijI2(3,2))+2*DCMPLX(CijR2(4,2),CijI2(4,2)))
cFC %
cFC %
c          print*, "ten2g",ten2g,ten2p1p2,test2g,test2p1p2
c          print*,  "Accuracy",accuracyC(0,5),accuracyC(1,5),accuracyC(2,5),accuracyC(3,5),accuracyC(4,5)
c          print*,  "Accuracy",accuracyC(0,4),accuracyC(1,4),accuracyC(2,4),accuracyC(3,4),accuracyC(4,4)

	 If(abs(test2g).gt.1d-6) then
         ten2g     =abs(ten2g/test2g-1d0)
	 ten2g_Dble=abs(ten2g_Dble/test2g-1d0)
	  else
	 ten2g     =abs(ten2g-test2g)
	 ten2g_Dble=abs(ten2g_Dble-test2g)
	 endif
	  
	  
	  If(abs(test2p1p2).gt.1d-6) then
         ten2p1p2     =abs(ten2p1p2/test2p1p2-1d0)
	 ten2p1p2_Dble=abs(ten2p1p2_Dble/test2p1p2-1d0)
	  else
	 ten2p1p2     =abs(ten2p1p2-test2p1p2)
	 ten2p1p2_Dble=abs(ten2p1p2_Dble-test2p1p2)
	 endif


c Check NAN
        if(ten2g+1d0.eq.ten2g)ten2g= 1d99
        if(ten2p1p2+1d0.eq.ten2p1p2)ten2p1p2= 1d99
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
	  
	ten2g_dble=Max(abs(ten2g_dble),abs(ten2p1p2_Dble)) 
	ten2g=Max(abs(ten2g),abs(ten2p1p2))


cFC %        if(1d0*abs(ten2g).gt.abs(ten2g_Dble)) then
cFC %	   call MyCget(p1sq,p2sq,s12,musq,C0123R,C0123I,Cij)
cFC %           C0123=DCMPLX(C0123R,C0123I)
cFC %	endif


       return
       end
