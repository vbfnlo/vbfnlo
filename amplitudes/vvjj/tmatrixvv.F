ccccccccccccccccccccccccc
cc K-matrix unitarization
      SUBROUTINE TMATRIXVV(Aoff,Aoo,s1,pp1,pp2,pp3,pp4,dim,blockdim, Au)

      IMPLICIT NONE
      INTEGER DIM,BLOCKDIM      
      DOUBLE COMPLEX Aoff(-1:1,-1:1,-1:1,-1:1,1:blockdim,1:blockdim,0:2)
      DOUBLE COMPLEX Aon(-1:1,-1:1,-1:1,-1:1,1:blockdim,1:blockdim,0:2)      
      DOUBLE COMPLEX Aoo(-1:1,-1:1,-1:1,-1:1,1:blockdim,1:blockdim,0:2)
      DOUBLE COMPLEX Au(-1:1,-1:1,-1:1,-1:1,1:blockdim,1:blockdim,0:2)
      

      DOUBLE PRECISION S1, PP1, PP2, PP3, PP4
      DOUBLE PRECISION ctte
      DOUBLE COMPLEX trA(0:2),trA2(0:2),trAfrob(0:2)
      
      DOUBLE PRECISION rZERO, rONE, rPI, rHALF
      PARAMETER(rZERO = 0.0D0, rONE = 1.0D0, rHALF=0.5D0)
      PARAMETER(rPI = 3.14159265358979323846d0)
      DOUBLE COMPLEX cZERO, cONE, cIMAG
      PARAMETER(cZERO=(0.0D0,0.0D0),cONE=(1.0D0,0.0D0),cIMAG=(0.0D0,1.0D0))      
      
      INTEGER I,J,K,L,N,M
      INTEGER LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA, BLOCK1, BLOCK2, BLOCK
      INTEGER ind(-1:1,-1:1,1:blockdim)
      
      INTEGER ier   
      DOUBLE PRECISION acc
      DOUBLE COMPLEX Aeigen(dim,dim)
      DOUBLE PRECISION AArray(dim)
      
      DOUBLE COMPLEX Toff(dim,dim,0:2), Ton(dim,dim,0:2)
      DOUBLE COMPLEX Too(dim,dim,0:2), Tu(dim,dim,0:2)
      DOUBLE COMPLEX A(dim,dim,0:2), Binv(0:2)
      DOUBLE COMPLEX C(dim,dim,0:2), D(dim,dim,0:2)
      DOUBLE COMPLEX Asqr(dim,dim,0:2),Afrob(dim,dim,0:2)
      DOUBLE COMPLEX mONE(dim,dim)   
      DOUBLE PRECISION fullEigen(0:2)
      
      
      DOUBLE PRECISION LMBDon,LMBDoff
      DOUBLE PRECISION Ctt, Cst
      
      DOUBLE PRECISION start, finish
      logical fullMaxEigen
      integer Eigenaprox
      integer fullvalue
      
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/an_couplings.inc"
#include "eigenV.inc"
      fullMaxEigen=.true.
      Eigenaprox=2
      fullvalue=1
      call settoCU(dim,dim,mOne)
      Tu(1:dim,1:dim,0:2) = cZero
      Toff(1:dim,1:dim,0:2) = cZero      
      Ton(1:dim,1:dim,0:2) = cZero
      Too(1:dim,1:dim,0:2) = cZero
      A(1:dim,1:dim,0:2) = cZero
      Binv(0:2) = cZero
      Asqr(1:dim,1:dim,0:2) = cZero   
      Afrob(1:dim,1:dim,0:2)= cZero
      C(1:dim,1:dim,0:2) = cZero
      D(1:dim,1:dim,0:2) = cZero
      ratioEigen(0:2) = rZero
      fullEigen(0:2) = rZERO
      maxEigen(0:2) = rZERO
      AArray(1:dim) = rZERO
      Au(-1:1,-1:1,-1:1,-1:1,1:blockdim,1:blockdim,0:2)=cZERO
****************************************************************   
****************************************************************         
!Phase space couplings, needed for correction.
       call TAmp_normalizationvv(s1,pp3,pp4, LMBDon)
       call TAmp_normalizationvv(s1,pp1,pp2, LMBDoff)
       if(LMBDon.lt.1.0d-16*s1*s1)then
        LMBDon=1.0d-16*s1*s1
       endif
       if(LMBDoff.lt.1.0d-16*s1*s1)then
        LMBDoff=1.0d-16*s1*s1
       endif
       Ctt = sqrt(LMBDon)/s1         !On->On
       Cst =(LMBDon**(0.25)/sqrt(s1))*(LMBDoff**(0.25)/sqrt(s1)) !on->off or off->on

****************************************************************                
* from A(l1,l2,l3,l4,b1,b2,j) to To(i,k,j)  
      ind(-1:1,-1:1,1:blockdim)=0

      select case (dim)
        case (21)
          do block = 1,3,2
            ind(-1,-1,block) = 1+(block-1)*9-(block-1)*3
            ind(-1,0,block)  = 2+(block-1)*9-(block-1)*3
            ind(-1,+1,block) = 3+(block-1)*9-(block-1)*3
            ind(0,-1,block)  = 4+(block-1)*9-(block-1)*3
            ind(0,0,block)   = 5+(block-1)*9-(block-1)*3
            ind(0,+1,block)  = 6+(block-1)*9-(block-1)*3
            ind(+1,-1,block) = 7+(block-1)*9-(block-1)*3
            ind(+1,0,block)  = 8+(block-1)*9-(block-1)*3
            ind(+1,+1,block) = 9+(block-1)*9-(block-1)*3
          enddo
          ind(-1,0,2)=10
          ind( 0,0,2)=11
          ind(+1,0,2)=12   

        case (9)
          ind(-1,-1,1) = 1
          ind(-1,0,1)  = 2
          ind(-1,+1,1) = 3
          ind(0,-1,1)  = 4
          ind(0,0,1)   = 5
          ind(0,+1,1)  = 6
          ind(+1,-1,1) = 7
          ind(+1,0,1)  = 8
          ind(+1,+1,1) = 9
        case(67)
          do lambda1 = -1,1
            do lambda2 = -1,1
              ind(lambda1,lambda2,1) = (lambda2+2)+(lambda1+1)*3
              ind(lambda1,0,2) = 9+(lambda1+2)
              ind(lambda1,lambda2,3) = 12+(lambda2+2)+(lambda1+1)*3
              ind(0,lambda2,4) = 21+(lambda2+2)
              ind(0,0,5) = 25
              ind(0,lambda2,6) = 25+(lambda2+2)
              ind(lambda1,lambda2,7) = 28+(lambda2+2)+(lambda1+1)*3
              ind(lambda1,0,8) = 37+(lambda1+2)
              ind(lambda1,lambda2,9) = 40+(lambda2+2)+(lambda1+1)*3
              ind(lambda1,lambda2,10) = 49+(lambda2+2)+(lambda1+1)*3
              ind(lambda1,lambda2,11) = 58+(lambda2+2)+(lambda1+1)*3
            enddo
          enddo
          case(42)
          do lambda1 = -1,1
            do lambda2 = -1,1
              ind(lambda1,lambda2,1) = (lambda2+2)+(lambda1+1)*3
              ind(0,lambda2,2) = 9+(lambda2+2)
              ind(lambda1,lambda2,3) = 12+(lambda2+2)+(lambda1+1)*3
              ind(lambda1,lambda2,4) = 21+(lambda2+2)+(lambda1+1)*3
              ind(lambda1,0,5) = 30+(lambda1+2)
              ind(lambda1,lambda2,6) = 33+(lambda2+2)+(lambda1+1)*3
            enddo
          enddo
          
      end select


      


    
      do i = 0,2
      do lambda1 = -1,1
       do lambda2 = -1,1
        do lambda3 = -1,1
         do lambda4 = -1,1 
          do block1 = 1,blockdim                 
           do block2 = 1,blockdim                
            Aon(lambda1,lambda2,lambda3,lambda4,block1,block2,i)=conjg(Aoff(lambda3,lambda4,lambda1,lambda2,block2,block1,i)) 
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo     
      enddo
              

      do lambda1 = -1,1
       do lambda2 = -1,1
        do lambda3 = -1,1
         do lambda4 = -1,1 
          do block1 = 1,blockdim                      
           do block2 = 1,blockdim                     
            if((ind(lambda1,lambda2,block1).ne.0).and.(ind(lambda3,lambda4,block2).ne.0)) then       
!!! Voff(1,2) -> Von(3,4)                           
             Toff(ind(lambda1,lambda2,block1),ind(lambda3,lambda4,block2),0)= Cst*Aoff(lambda1,lambda2,lambda3,lambda4,block1,block2,0)
             Toff(ind(lambda1,lambda2,block1),ind(lambda3,lambda4,block2),1)= Cst*Aoff(lambda1,lambda2,lambda3,lambda4,block1,block2,1)
             Toff(ind(lambda1,lambda2,block1),ind(lambda3,lambda4,block2),2)= Cst*Aoff(lambda1,lambda2,lambda3,lambda4,block1,block2,2)  
!!! Von(3,4)  -> Voff(1,2)
             Ton(ind(lambda1,lambda2,block1),ind(lambda3,lambda4, block2),0)= Cst*Aon(lambda1,lambda2,lambda3,lambda4,block1,block2,0)
             Ton(ind(lambda1,lambda2,block1),ind(lambda3,lambda4, block2),1)= Cst*Aon(lambda1,lambda2,lambda3,lambda4,block1,block2,1)
             Ton(ind(lambda1,lambda2,block1),ind(lambda3,lambda4, block2),2)= Cst*Aon(lambda1,lambda2,lambda3,lambda4,block1,block2,2)      
!!! Von(3,4)  -> Von(3,4)           
             Too(ind(lambda1,lambda2,block1),ind(lambda3,lambda4, block2),0)= Ctt*Aoo(lambda1,lambda2,lambda3,lambda4,block1,block2,0)
             Too(ind(lambda1,lambda2,block1),ind(lambda3,lambda4, block2),1)= Ctt*Aoo(lambda1,lambda2,lambda3,lambda4,block1,block2,1)
             Too(ind(lambda1,lambda2,block1),ind(lambda3,lambda4, block2),2)= Ctt*Aoo(lambda1,lambda2,lambda3,lambda4,block1,block2,2) 
            endif           
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
***      

      do j=0,2

       do i = 1,dim
        do k=1,dim
         do l = 1,dim
         A(i,k,j) = A(i,k,j) + Ton(i,l,j)*Toff(l,k,j)
         
         enddo
        enddo
       enddo 


         
         

!!!!!!!!!!!!!!!!!!!!!!! 

       if(.not.fullMaxEigen) then 
**       
         trA(j) = 0.0d0
         do i = 1,dim
           trA(j) = trA(j) + A(i,i,j)
         enddo 
     
**         
         if(eigenaprox.eq.1) then
          do i = 1,dim
           do k = 1,dim
            do l = 1,dim
             Asqr(i,k,j) = Asqr(i,k,j) + A(i,l,j)*A(l,k,j)
            enddo
           enddo
          enddo          
          trA2(j) = 0.0d0
          do i = 1,dim
           trA2(j) = trA2(j) + Asqr(i,i,j)
          enddo                
          ctte = 1.0d0/sqrt(1.0d0*dim)
          maxEigen(j) = ctte*(trA2(j)/trA(j))
**          
         else if(eigenaprox.eq.2) then
          do i = 1,dim
           do k = 1,dim
            do l = 1,dim
             Afrob(i,k,j) = Afrob(i,k,j) + A(i,l,j)*conjg(A(k,l,j))
            enddo
           enddo
          enddo          
          trAfrob(j) = 0.0d0
          do i = 1,dim
           trAfrob(j) = trAfrob(j) + Afrob(i,i,j)
          enddo           
          ctte = dim
          maxEigen(j) =  (1.0d0/ctte)*trA(j) + sqrt((ctte-1.0d0)/ctte)*sqrt(trAfrob(j)-((trA(j)**2)/ctte))  
         endif
**        
       else
        if(fullvalue.eq.1) then
          Aeigen(1:dim,1:dim) = cZero        
          acc = 1.0d-6
          Aeigen(1:dim,1:dim) = A(1:dim,1:dim,j)
          
          call eigenvalues(Aeigen, dim,acc, 10000, ier)
          do i=1,dim
            AArray(i) = Real(Aeigen(i,i))
          enddo
          maxEigen(j) = maxval(AArray)


        else if(fullvalue.eq.2) then
          call MAXEigenvalue(A(1,1,j),dim,maxEigen(j))
        endif
       endif   
       
                
!       

       Binv(j) = (1.0d0/(1.0d0 + maxEigen(j)/4.0d0))

       

       do i = 1,dim
        do k = 1,dim
         do l = 1,dim
          C(i,k,j) = C(i,k,j) + Toff(i,l,j)*Too(l,k,j)
         enddo
        enddo
       enddo
       do i = 1,dim
        do k = 1,dim
        D(i,k,j) = Toff(i,k,j) + (cIMAG/2.0d0)*C(i,k,j)
        enddo
       enddo
       
       do i =1,dim
        do k = 1,dim
         Tu(i,k,j) = Tu(i,k,j) + D(i,k,j)*Binv(j)
        enddo 
       enddo         

      
      enddo

****************************************************************               
****************************************************************            
 
****************************************************************      
****************************************************************      
      do lambda1 = -1,1
       do lambda2 = -1,1
        do lambda3 = -1,1
         do lambda4 = -1,1
          do block1 = 1,blockdim                       
           do block2 = 1,blockdim                      
            if((ind(lambda1,lambda2,block1).ne.0).and.(ind(lambda3,lambda4,block2).ne.0)) then        
            Au(lambda1,lambda2,lambda3,lambda4,block1,block2,0) = TU(ind(lambda1,lambda2,block1),ind(lambda3,lambda4,block2),0)*(rONE/Cst)
            Au(lambda1,lambda2,lambda3,lambda4,block1,block2,1) = TU(ind(lambda1,lambda2,block1),ind(lambda3,lambda4,block2),1)*(rONE/Cst)
            Au(lambda1,lambda2,lambda3,lambda4,block1,block2,2) = TU(ind(lambda1,lambda2,block1),ind(lambda3,lambda4,block2),2)*(rONE/Cst)
           

            endif
           enddo
          enddo
         enddo 
        enddo
       enddo       
      enddo

                      
      
      END
      
ccccccccccccccccccccccccc
cc Normalization of states. Zuber, p. 242

      SUBROUTINE TAmp_normalizationvv(x1,x2,x3, LMBD)
      IMPLICIT NONE

      DOUBLE PRECISION X1, X2, X3, LMBD      
      LMBD = 0.0d0
      
      LMBD = (x1**2.0d0 + x2**2.0d0 + x3**2.0d0) - 2.0d0*x1*x2 - 2.0d0*x2*x3 - 2.0d0*x3*x1
      
      END
      
