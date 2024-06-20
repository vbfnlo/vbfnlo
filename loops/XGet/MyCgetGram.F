      subroutine MyCgetGram(p1sq,p2sq,s12,musq,C0123R,C0123I,Cij)
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
       Complex*16 C0123
      Real*8 C0123R
      Real*8 C0123I
      complex*16 Cij(30,9)
      real*8 tempjj,tempkl,IX,Zmax
      integer jjinit
      common/Decide/tempjj,tempkl,IX,Zmax,jjinit
      Save/Decide/


       real*8 musq,C0R,C0I
       integer jj,k,l,order
       
c       print*, "musq",musq

      tempjj=1d18
      tempkl=1d18
c      jjinit=12

       B012=B0fin(p1sq,musq)    
c       print*, 'B012', B012

       B013=B0fin(s12,musq)  
c        print*, 'B013', B013
       B023=B0fin(p2sq,musq) 
c       print*, 'B023', B023

      C0123=C0fin(p1sq,p2sq,s12,musq)  
      C0R=Dble(C0123)
      C0I=DImag(C0123)
  
c FC %%      call tens_red3_new_Re_Com(p1sq,p2sq,s12,B023,B013,B012, 
c FC %%     &                     C0123,C0r,C0I,Cijr,CijI)
  
      k=2
      l=2
      jj=2
      order=11

c FC %      call ten_red3_Gram(p1sq,p2sq,s12,B023,B013,B012,C0123,Cij,k,l,jj,order)

       call ten_red3_Gram1(p1sq,p2sq,s12,musq,B023,B013,B012,C0123,Cij)
c        print*, "tempCkl",tempkl,"tempCjj",tempjj,IX,ZMax,jjinit
       C0123R=Dble(C0123)
       C0123I=DImag(C0123)
       
c FC %%c CALL Dble
c FC %%       call MyCget1(p1sq,p2sq,s12,musq,C0R2,C0I2,CijR2,CijI2)
c FC %%c Checking accuracy  
c FC %%        p1p2 = (s12 - p1sq - p2sq)*0.5d0
c FC %%        r1 = p1sq
c FC %%        r2r1 = s12 - r1
c FC %%      
c FC %%	zz11=2d0*p1sq
c FC %%	zz12=2d0*p1p2
c FC %%	zz22=2d0*p2sq
c FC %%
c FC %%
c FC %%      test2g=2d0*B023
c FC %%
c FC %%      test2p1p2=-B013*(ZZ12+ZZ22)/2d0 +(r2r1)*B023
c FC %%     1   -r1*(B012-B013-(r2r1)*C0123)
c FC %%       
c FC %%       ten2g=ZZ11*Cij(1,2)+ZZ22*Cij(2,2)+ 
c FC %%     1  2d0*(ZZ12*Cij(3,2) + 4d0*Cij(4,2))
c FC %%     2 -1d0
c FC %%       
c FC %%       ten2p1p2=ZZ11*(ZZ12*Cij(1,2)+ZZ22*Cij(3,2)) + 
c FC %%     1  ZZ12*(ZZ22*Cij(2,2) + ZZ12*Cij(3,2)+2*Cij(4,2))
c FC %%         
c FC %%       ten2g_Dble=ZZ11*DCMPLX(CijR2(1,2),CijI2(1,2))+ZZ22*DCMPLX(CijR2(2,2),CijI2(2,2))+ 
c FC %%     1  2d0*(ZZ12*DCMPLX(CijR2(3,2),CijI2(3,2)) + 4d0*DCMPLX(CijR2(4,2),CijI2(4,2)))
c FC %%     2 -1d0
c FC %%    
c FC %%        ten2p1p2_Dble=ZZ11*(ZZ12*DCMPLX(CijR2(1,2),CijI2(1,2))+ZZ22*DCMPLX(CijR2(3,2),CijI2(3,2))) + 
c FC %%     1  ZZ12*(ZZ22*DCMPLX(CijR2(2,2),CijI2(2,2)) + ZZ12*DCMPLX(CijR2(3,2),CijI2(3,2))+2*DCMPLX(CijR2(4,2),CijI2(4,2)))
c FC %%
c FC %%
c FC %%	 If(abs(test2g).gt.1d-6) then
c FC %%         ten2g     =abs(ten2g/test2g-1d0)
c FC %%	 ten2g_Dble=abs(ten2g_Dble/test2g-1d0)
c FC %%	  else
c FC %%	 ten2g     =abs(ten2g-test2g)
c FC %%	 ten2g_Dble=abs(ten2g_Dble-test2g)
c FC %%	 endif
c FC %%	  
c FC %%	  
c FC %%	  If(abs(test2p1p2).gt.1d-6) then
c FC %%         ten2p1p2     =abs(ten2p1p2/test2p1p2-1d0)
c FC %%	 ten2p1p2_Dble=abs(ten2p1p2_Dble/test2p1p2-1d0)
c FC %%	  else
c FC %%	 ten2p1p2     =abs(ten2p1p2-test2p1p2)
c FC %%	 ten2p1p2_Dble=abs(ten2p1p2_Dble-test2p1p2)
c FC %%	 endif
c FC %%	  
c FC %%	ten2g_dble=Max(abs(ten2g_dble),abs(ten2p1p2_Dble)) 
c FC %%	ten2g=Max(abs(ten2g),abs(ten2p1p2)) 
c FC %%
c FC %%	 
c FC %%        if(1d0*abs(ten2g).gt.abs(ten2g_Dble)) then
c FC %%	   call MyCget(p1sq,p2sq,s12,musq,C0123R,C0123I,Cij)
c FC %%           C0123=DCMPLX(C0123R,C0123I)
c FC %%	endif


       return
       end
