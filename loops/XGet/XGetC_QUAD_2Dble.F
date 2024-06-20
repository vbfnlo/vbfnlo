       subroutine PrintIntVert_QUAD(m_d,p1,p2,p3,musq)
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 27/10/2008
c Modified:17/11/2008
c ************************************************************************************
c determine the Hexagon of
c p1(gluon)    ---->------->--------->--------->---   p6(B) p6^2=M_h
c                $         $          $         $     
c                $         $          $         $     
c                $         $          $         $     
c                V         V          V         V     
c                $         $          $         $     
c                $         $          $         $     
c               p2         p3         p4        p5     
c              mu_p2      mu_p3      mu_p4     mu_p5      
c Note: No couplings are assigned
c Note:m the mass of the fermion in the loop  
c Note: To make it shorter in the promgram: mu_p2,...->mup2,... 
c Notation of External momenta: p1+p2+p3+p4+p5+p6=0 
c mu_p1,mu_p2,mu_p3,mu_p4,mu_p5 should be think as external current 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3)
      Real*8 psq(3)
        Real*8   p1sq, p1p2, p1p3
       Real*8   p2sq, p2p3
       Real*8   p3sq
       Real*8   s12
       Complex*16  B012,B013
       Complex*16  B023
       Real*8  B012R,B013R
       Real*8  B023R
       Real*8  B012I,B013I
       Real*8  B023I
       Complex*16 C0123
       Real*8 C0123R
       Real*8 C0123I
       Real*8 C123R(4,2)
       Real*8 C123I(4,2)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  QUAD precision start here
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       Real*16  p1_QUAD(0:3),p2_QUAD(0:3),p3_QUAD(0:3)
       Real*16 ptemp(0:3,3)
       Real*16   p1sq_QUAD, p1p2_QUAD, p1p3_QUAD
       Real*16   p2sq_QUAD, p2p3_QUAD
       Real*16   p3sq_QUAD
       Real*16   s12_QUAD
       Complex*32  B012_QUAD,B013_QUAD
       Complex*32  B023_QUAD
       Real*16  B012R_QUAD,B013R_QUAD
       Real*16  B023R_QUAD
       Real*16  B012I_QUAD,B013I_QUAD
       Real*16  B023I_QUAD
       Complex*32 C0123_QUAD
       Real*16 C0123R_QUAD
       Real*16 C0123I_QUAD
       Real*16 C123R_QUAD(4,2)
       Real*16 C123I_QUAD(4,2)
       Real*16 musq_QUAD, dotrr_QUAD
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      End QUAD precision
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       COMMON/CVALFunctionsC/ C0123R,C123R,C0123I,C123I 

       COMMON/BVALFunctionsC/ B012R,B012I, B013R,B013I,
     &      B023R,B023I
       COMMON/InvariantsC/p1sq,p2sq,p3sq,s12
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       Complex*16 B0fin,C0fin
       EXTERNAL B0fin,C0fin
       Logical PrintB,PrintC
       Integer comp  
       External dotrr_QUAD ,dotrr
      Complex*16 dotrc,dotcc
       Real*8 dotrr
       Real*8 musq
      EXTERNAL   dotrc,dotcc
       Integer i,j,jtemp
       integer j1,j2,j3,j4,j5,j6
       real*8 m_d
       real*16 m
      double precision     hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass
      common /COUPL_MASS/  hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass      
       real*16 c,d,e03,e04
       integer zaxis
       real*8 limit
       real*8 temphiggs(2),tempdiff
       limit=1d-14
c      Integer alpha
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the scalar products. Not inlcueded the contraction of the
c       moments with the external currents  
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c  New recunstraction of momenta

        m=m_d

       do i=0,3 
       ptemp(i,1)=p1(i)
       ptemp(i,2)=p2(i)
       ptemp(i,3)=p3(i)
       enddo
       psq(1)=dotrr(p1,p1)
       psq(2)=dotrr(p2,p2)
       psq(3)=dotrr(p3,p3)

cfc       do jtemp=1,4
cfc          print*, "jtemp",jtemp
cfc       do i=0,3
cfc          print*, ptemp(i,jtemp)
cfc       enddo
cfc       print*
cfc       enddo



cccccccc Identify the Higgs cccccccccccccccccccccccc
       j1=0
       temphiggs(1)=1d99
       temphiggs(2)=1d99
       do i=1,3
         temphiggs(2)=abs(psq(i)/(HMASS*HMASS)-1d0)
       if((j1.eq.0).or.(temphiggs(2).lt.temphiggs(1))) then
       if(temphiggs(2).lt.1d-6) then
        temphiggs(1)=temphiggs(2)
        j1=i
       endif
       endif
       enddo
 
ccccccccccccccccccccccccccccccccccc
c  H
       if(j1.ne.0) then
       jtemp=j1
       if (ptemp(0,jtemp).ge.0d0) then
       ptemp(0,jtemp)=SQRT(HMASS*HMASS
     &                 +ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
       else
       ptemp(0,jtemp)=-SQRT(HMASS*HMASS
     &                 +ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
       endif
       endif
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c  identify gluons in the axis and out of the axis
       j2=0
       j3=0
       j4=0
       j5=0
       j6=0
ccccccccccccccccccccccccccccccccccccccccccccccccccc  
       zaxis=0

      if(j1.eq.0) then

       do jtemp=1,3
c          print*, "j",jtemp
       if(abs(ptemp(1,jtemp)).le.10d-7.and.
     &    abs(ptemp(2,jtemp)).le.10d-7)then
c p_0 >0 axes
       if(j3.eq.0) then
       j3=jtemp
       zaxis=zaxis+1
       else
       j4=jtemp
       zaxis=zaxis + 1
       endif
       else
       if(j6.eq.0) then
       j6=jtemp
       elseif(j5.eq.0) then
       j5=jtemp
       elseif(j2.eq.0) then
       j2=jtemp
       endif
       endif
       enddo

       else !j1 ne. 0

       do jtemp=1,3
       if(jtemp.ne.j1.and.jtemp.ne.j2) then
c          print*, "j",jtemp
       if(abs(ptemp(1,jtemp)).le.10d-7.and.
     &    abs(ptemp(2,jtemp)).le.10d-7)then
c p_0 >0 axes
       if(j3.eq.0) then
       j3=jtemp
       zaxis=zaxis+1
       else
       j4=jtemp
       zaxis=zaxis + 1
       endif
       else
       if(j5.eq.0) then
       j5=jtemp
       elseif(j2.eq.0) then
       j2=jtemp
       endif
       endif
       endif
       enddo
       endif

       if(zaxis.eq.2) then 
cfc       print*, "j1",j1
cfc       print*, "j2",j2
cfc       print*, "j3",j3
cfc       print*, "j4",j4
cfc       print*, "j5",j5
c       print*, "j6",j6       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c On-shell for second gluon out of axis:
       if(j6.ne.0) then
       jtemp=j6

       if(abs(psq(jtemp)).lt.1d-7) then

       if (ptemp(0,jtemp).ge.0d0) then
       ptemp(0,jtemp)=SQRT(ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
       else
      ptemp(0,jtemp)=-SQRT(ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
      endif
      endif
      endif
        

cfcccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cfcc On-shell for second gluon out of axis:
cfc       jtemp=j2
cfc
cfc       if(abs(psq(jtemp)).lt.1d-7) then
cfc
cfc       if (ptemp(0,jtemp).ge.0d0) then
cfc       ptemp(0,jtemp)=SQRT(ptemp(1,jtemp)*ptemp(1,jtemp)+
cfc     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
cfc     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
cfc       else
cfc      ptemp(0,jtemp)=-SQRT(ptemp(1,jtemp)*ptemp(1,jtemp)+
cfc     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
cfc     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
cfc      endif
cfc      endif
cfcccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cfcc On-shell for second gluon out of axis:
cfc       jtemp=j6
cfc       if (ptemp(0,jtemp).ge.0d0) then
cfc       ptemp(0,jtemp)=SQRT(ptemp(1,jtemp)*ptemp(1,jtemp)+
cfc     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
cfc     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
cfc       else
cfc      ptemp(0,jtemp)=-SQRT(ptemp(1,jtemp)*ptemp(1,jtemp)+
cfc     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
cfc     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
cfc      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Impose x, y conservation gluon out of axes with
c H bosons

c          do i=0,3
c            write(*,'(A3,I1,A2,E25.18,E25.18)') "p1(",i,")=",ptemp(i,j5)
c            enddo
       jtemp=j5
c       print*,"j5",j5
c       print*,"x",ptemp(1,j5)+ptemp(1,j1)+ptemp(1,j2)
c       print*,"y",ptemp(2,j5)+ptemp(2,j1)+ptemp(2,j2)
       if(j6.eq.0) then
       ptemp(1,jtemp)=-ptemp(1,j1)!-ptemp(1,j2)!-ptemp(1,j6)
       ptemp(2,jtemp)=-ptemp(2,j1)!-ptemp(2,j2)!-ptemp(2,j6)
       else
       ptemp(1,jtemp)=-ptemp(1,j6)!-ptemp(1,j2)
       ptemp(2,jtemp)=-ptemp(2,j6)!-ptemp(2,j2)    
       endif

c       print*,"x",ptemp(1,j5)+ptemp(1,j1)+ptemp(1,j2)
c       print*,"y",ptemp(2,j5)+ptemp(2,j1)+ptemp(2,j2)
      if(abs(psq(jtemp)).lt.1d-7) then
       if (ptemp(0,jtemp).ge.0d0) then
       ptemp(0,jtemp)=SQRT(
     &                 +ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
       else
       ptemp(0,jtemp)=-SQRT(
     &                 +ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
       endif
       endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Now we have to reconstuct the pgluons of the z axis
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       print*, "j1",j1
c       print*, "j2",j2
c       print*, "j3",j3
c       print*, "j4",j4
c       print*, "j5",j5

      jtemp=j3
      if(j6.eq.0) then
         if(j2.eq.0) then
      c= ptemp(0,j1)+ ptemp(0,j5)!+ ptemp(0,j6)
      d= ptemp(3,j1)+ ptemp(3,j5)!+ ptemp(3,j6)
         else
      c= ptemp(0,j1)+ ptemp(0,j2)+ ptemp(0,j5)!+ ptemp(0,j6)
      d= ptemp(3,j1)+ ptemp(3,j2)+ ptemp(3,j5)!+ ptemp(3,j6)
         endif
      else
         if(j2.eq.0) then
      c=   ptemp(0,j5)+ ptemp(0,j6)
      d=   ptemp(3,j5)+ ptemp(3,j6)
          else
      c=  ptemp(0,j2)+ ptemp(0,j5)+ ptemp(0,j6)
      d=  ptemp(3,j2)+ ptemp(3,j5)+ ptemp(3,j6)
         endif
      endif

      if (ptemp(0,jtemp).ge.0d0) then   
         If(ptemp(3,jtemp).lt.0d0) then
            if (ptemp(0,j4).ge.0d0) then   
               If(ptemp(3,j4).ge.0d0) then
c     +- ++
                   e04=(-c-d)/2q0
                   e03=(d-c)/2q0
c j3
                   ptemp(0,jtemp)=e03
                   ptemp(1,jtemp)=0q0
                   ptemp(2,jtemp)=0q0
                   ptemp(3,jtemp)=-e03
c j4
                   ptemp(0,j4)=e04
                   ptemp(1,j4)=0q0
                   ptemp(2,j4)=0q0
                   ptemp(3,j4)=e04
                   else
c   +-   +-
                    print*,"PenQUAD_queck momenta"
                   endif
            else  
c + -   - +
                  If(ptemp(3,j4).ge.0d0) then
                    print*,"PenQUAD_queck momenta"
                  else
c +-     - -        
                   e04=(c+d)/2q0
                   e03=(d-c)/2q0
                   ptemp(0,jtemp)=e03
                   ptemp(1,jtemp)=0q0
                   ptemp(2,jtemp)=0q0
                   ptemp(3,jtemp)=-e03
c j4
                   ptemp(0,j4)=-e04
                   ptemp(1,j4)=0q0
                   ptemp(2,j4)=0q0
                   ptemp(3,j4)=-e04
                 endif   
            endif
         else  
c   + +
            if (ptemp(0,j4).ge.0d0) then   
               If(ptemp(3,j4).ge.0d0) then
c     ++     ++
                print*,"PenQUAD_queck momenta"
                   else
c   ++    +-
c                      print*, "HERE:recons"
                e04=(d-c)/2q0
                e03=(-c-d)/2q0
                   ptemp(0,jtemp)=e03
                   ptemp(1,jtemp)=0q0
                   ptemp(2,jtemp)=0q0
                   ptemp(3,jtemp)=e03
c j4
                   ptemp(0,j4)=e04
                   ptemp(1,j4)=0q0
                   ptemp(2,j4)=0q0
                   ptemp(3,j4)=-e04
                   endif
            else  
c + +   - +
                  If(ptemp(3,j4).ge.0d0) then
                e04=(-d-c)/2q0
                e03=(d-c)/2q0
                   ptemp(0,jtemp)=e03
                   ptemp(1,jtemp)=0q0
                   ptemp(2,jtemp)=0q0
                   ptemp(3,jtemp)=e03
c j4
                   ptemp(0,j4)=-e04
                   ptemp(1,j4)=0q0
                   ptemp(2,j4)=0q0
                   ptemp(3,j4)=e04
                  else
c ++     - -        
                print*,"PenQUAD_queck momenta"
                 endif   
            endif
        endif
      else
c -
         If(ptemp(3,jtemp).lt.0d0) then
            if (ptemp(0,j4).ge.0d0) then   
               If(ptemp(3,j4).ge.0d0) then
c     -- ++
                    print*,"PenQUAD_queck momenta"
                   else
c   --   +-
                   e04=(d-c)/2q0
                   e03=(c+d)/2q0
                   ptemp(0,jtemp)=-e03
                   ptemp(1,jtemp)=0q0
                   ptemp(2,jtemp)=0q0
                   ptemp(3,jtemp)=-e03
c j4
                   ptemp(0,j4)=e04
                   ptemp(1,j4)=0q0
                   ptemp(2,j4)=0q0
                   ptemp(3,j4)=-e04

                   endif
            else  
c  - -   - +
                  If(ptemp(3,j4).ge.0d0) then
                   e04=(c-d)/2q0
                   e03=(c+d)/2q0
                   ptemp(0,jtemp)=-e03
                   ptemp(1,jtemp)=0q0
                   ptemp(2,jtemp)=0q0
                   ptemp(3,jtemp)=-e03
c j4
                   ptemp(0,j4)=-e04
                   ptemp(1,j4)=0q0
                   ptemp(2,j4)=0q0
                   ptemp(3,j4)=e04
                  else
c --     - -        
                 print*,"PenQUAD_queck momenta"
                 endif   
            endif
         else  
c   - +
            if (ptemp(0,j4).ge.0d0) then   
               If(ptemp(3,j4).ge.0d0) then
c     -+     ++
                   e04=(-c-d)/2q0
                   e03=(c-d)/2q0
                   ptemp(0,jtemp)=-e03
                   ptemp(1,jtemp)=0q0
                   ptemp(2,jtemp)=0q0
                   ptemp(3,jtemp)=+e03
c j4
                   ptemp(0,j4)=e04
                   ptemp(1,j4)=0q0
                   ptemp(2,j4)=0q0
                   ptemp(3,j4)=e04
                   else
c   -+    +-
                print*,"PenQUAD_queck momenta"
                   endif
            else  
c - +   - +
                  If(ptemp(3,j4).ge.0d0) then
                print*,"PenQUAD_queck momenta"
                  else
c -+     - -        
                   e04=(c+d)/2q0
                   e03=(c-d)/2q0
                   ptemp(0,jtemp)=-e03
                   ptemp(1,jtemp)=0q0
                   ptemp(2,jtemp)=0q0
                   ptemp(3,jtemp)=+e03
c j4
                   ptemp(0,j4)=-e04
                   ptemp(1,j4)=0q0
                   ptemp(2,j4)=0q0
                   ptemp(3,j4)=-e04
                 endif   
            endif
        endif
       endif 


       else ! zaxis no equal=2, Then, some gluon is off-shell

 

        do jtemp=1,3
           if(jtemp.ne.j1) then
       if(abs(psq(jtemp)).le.1d-7) THEN 
             J3=JTEMP
           if (ptemp(0,jtemp).ge.0d0) then
           ptemp(0,jtemp)=SQRT(
     &                     +ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                      ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                      ptemp(3,jtemp)*ptemp(3,jtemp))
           else
           ptemp(0,jtemp)=-SQRT(
     &                     +ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                      ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                      ptemp(3,jtemp)*ptemp(3,jtemp))
           endif
   
      else
         j4=jtemp
      ENDIF
      endif
      ENDDO
 

cc off-shell 
       do i=0,3
       ptemp(i,j4)=0q0
       do jtemp=1,3
       if (jtemp.ne.j4) then
       ptemp(i,j4)=ptemp(i,j4)-ptemp(i,jtemp) 
       endif
       enddo
       enddo
       

       endif  ! Zaxis 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      print*
c      do i=0,3
c       print*, "XCGET_QUAD:conv",i,ptemp(i,1)+ptemp(i,2)+ptemp(i,3)
c      enddo
c     i=1

      
      do i=0,3
         if(abs(p1(i)).lt.1d-3) then
            tempdiff=abs(p1(i)-ptemp(i,1))
          else
             tempdiff=abs(p1(i)/ptemp(i,1)-1d0)
          endif
         if(tempdiff.gt.1d-6) then
            write(*,'(A5,I1,A2,E25.16,E25.16)') "Cp1(",i,")=",p1(i)-ptemp(i,1)
            limit=1d0
            GOTO 30
c            stop
         endif
      enddo
      do i=0,3
         if(abs(p2(i)).lt.1d-3) then
            tempdiff=abs(p2(i)-ptemp(i,2))
          else
             tempdiff=abs(p2(i)/ptemp(i,2)-1d0)
          endif
         if(tempdiff.gt.1d-6)then
            write(*,'(A5,I1,A2,E25.16,E25.16)') "Cp2(",i,")=",p2(i)-ptemp(i,2)
         limit=1d0
            GOTO 30
c            stop
         endif
      enddo
      do i=0,3 
          if(abs(p3(i)).lt.1d-3) then
            tempdiff=abs(p3(i)-ptemp(i,3))
          else
             tempdiff=abs(p3(i)/ptemp(i,3)-1d0)
          endif
         if(tempdiff.gt.1d-6)then
            write(*,'(A5,I1,A2,E25.16,E25.16)') "Cp3(",i,")=",p3(i)-ptemp(i,3)
         limit=1d0
            GOTO 30
c            stop
         endif
      enddo
      
       do i=0,3
       if(abs(ptemp(i,1)+ptemp(i,2)+ptemp(i,3)
     &   ).gt.limit) then
 30      jtemp=jtemp
         do jtemp=0,3
        print*, "XCGET_QUAD:conv",jtemp,ptemp(jtemp,1)+ptemp(jtemp,2)+ptemp(jtemp,3)
               print*, "XCGET_Dble:conv",jtemp,p1(jtemp)+p2(jtemp)+p3(jtemp)
         enddo

         print*
         print*, "j1",j1
         print*, "j2",j2
         print*, "j3",j3
         print*, "j4",j4
         print*, "j5",j5
         print*, "j6",j6   
         Print*, "HERE", zaxis
         print*, "j4",j4,psq(j4)
         print*, "j1",j1,psq(j1)
         do jtemp=1,3
            print*, "jtemp",jtemp,psq(jtemp)
            print*
         enddo
         
         do jtemp=0,3
            print*, p1(jtemp)
         enddo
         print*
         
         do jtemp=0,3
            print*, p2(jtemp)
         enddo
         print*
         
         do jtemp=0,3
            print*, p3(jtemp)
         enddo
         print*
         
         do jtemp=1,3
            print*, "jtemp",jtemp
            do j=0,3
               print*, ptemp(j,jtemp)
            enddo
            print*
         enddo
         
c         stop
      endif
      enddo


      do i=0,3 
      p1_QUAD(i)=ptemp(i,1)
      p2_QUAD(i)=ptemp(i,2)
      p3_QUAD(i)=ptemp(i,3)
       enddo

       p1sq_QUAD = dotrr_QUAD(p1_QUAD,p1_QUAD)
       p1p2_QUAD = dotrr_QUAD(p1_QUAD,p2_QUAD)
       p1p3_QUAD = dotrr_QUAD(p1_QUAD,p3_QUAD)
       p2sq_QUAD = dotrr_QUAD(p2_QUAD,p2_QUAD)
       p2p3_QUAD = dotrr_QUAD(p2_QUAD,p3_QUAD)
       p3sq_QUAD = dotrr_QUAD(p3_QUAD,p3_QUAD)
       s12_QUAD = (p1sq_QUAD +p2sq_QUAD+ 2*p1p2_QUAD) 


cfc       Write(*,'(a5,E44.36)')," p1sq_QUAD ", p1sq_QUAD 
cfc       Write(*,'(a5,E44.36)')," p1p2_QUAD ", p1p2_QUAD
cfc       Write(*,'(a5,E44.36)')," p1p3_QUAD ", p1p3_QUAD
cfc       Write(*,'(a5,E44.36)')," p2sq_QUAD ", p2sq_QUAD 
cfc       Write(*,'(a5,E44.36)')," p2p3_QUAD ", p2p3_QUAD
cfc       Write(*,'(a5,E44.36)')," p3sq_QUAD ", p3sq_QUAD 


       p1sq=p1sq_QUAD 
       p2sq=p2sq_QUAD 
       p3sq=p3sq_QUAD 
       s12=s12_QUAD


       do i=0,3
         p1(i)=p1_QUAD(i)
         p2(i)=p2_QUAD(i)
         p3(i)=p3_QUAD(i)
      enddo

      If(abs(musq-s12).lt.10d-6) then
          musq_QUAD=s12_QUAD
          else
          musq_QUAD=musq
          endif
c      Write(*,*) 'musq',musq_QUAD
      PrintB=.false. 
      PrintC=.false. 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c        if (comp.gt.0) then    
c    Calling C_ij,D_ij,E_ij,F_ij Functions    
c************************************************************************************
c************************************************************************************
c*******************************************************************************

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                     B Functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       call GFHiggs_B0finG_QUAD(m,p1sq_QUAD,musq_QUAD,B012_QUAD)
       B012R_QUAD=QREAL(B012_QUAD)
       B012I_QUAD=QIMAG(B012_QUAD)
       call From_BijR_QUAD_to_Dbl(B012R,B012I,B012R_QUAD,B012I_QUAD)
       B012=DCMPLX(B012R,B012I)
       call GFHiggs_B0finG_QUAD(m,p2sq_QUAD,musq_QUAD,B023_QUAD)
       B023R_QUAD=QREAL(B023_QUAD)
       B023I_QUAD=QIMAG(B023_QUAD)
       call From_BijR_QUAD_to_Dbl(B023R,B023I,B023R_QUAD,B023I_QUAD)
       B023=DCMPLX(B023R,B023I)
       call GFHiggs_B0finG_QUAD(m,s12_QUAD,musq_QUAD,B013_QUAD)
       B013R_QUAD=QREAL(B013_QUAD)
       B013I_QUAD=QIMAG(B013_QUAD)
       call From_BijR_QUAD_to_Dbl(B013R,B013I,B013R_QUAD,B013I_QUAD)
       B013=DCMPLX(B013R,B013I)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                     C Functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       call C0finG1_QUAD(m,m,m,p1sq_QUAD,p2sq_QUAD,s12_QUAD,musq_QUAD,C0123_QU
     -   AD)
       call tens_red3_new_Re_Com_1M_QUAD(m,p1sq_QUAD,p2sq_QUAD,s12_QUAD,B02
     -   3_QUAD,B013_QUAD,B012_QUAD,C0123_QUAD,C0123R_QUAD,C0123I_QUA
     -   D,C123R_QUAD,C123I_QUAD)
       C0123R=0.d0
       call  From_CijR_QUAD_to_Dbl(C0123R,C0123I,C123R,C123I,C0123R_QU
     -   AD,C0123I_QUAD,C123R_QUAD,C123I_QUAD)
       C0123=DCMPLX(C0123R,C0123I)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c   Printing C_ij,D_ij,E_ij Functions    
c************************************************************************************
c************************************************************************************
       If(PrintB) then
       write(*,*)"B12"
       write(*,*) "B12=", DCMPLX(B012R,B012I)
       write(*,*)"B23"
       write(*,*) "B23=", DCMPLX(B023R,B023I)
       write(*,*)"B13"
       write(*,*) "B13=", DCMPLX(B013R,B013I)
      Endif    
       If(PrintC) then
       write(*,*)"C123"
       call print_CijR_QUAD(C0123_QUAD,C123R_QUAD,C123I_QUAD)
       call print_CijR_QUAD_D(C0123,C123R,C123I,C0123_QUAD,C123R_QUAD,
     -   C123I_QUAD)
       call print_CijR(C0123,C123R,C123I)
       endif
c*****************

****************************************************************
       Return
       End
