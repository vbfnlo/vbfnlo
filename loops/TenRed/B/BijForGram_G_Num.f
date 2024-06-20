       subroutine ten_red2_forGram_G_Num(m0,m1,musq,p1sq,A2,A1,B012
     &                                  ,B12)
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
        Implicit none
        Real*8 p1sq,m0,m1,m0sq,m1sq,musq
        Complex*16 B012,A2,A1,A202,A402,A602,A802,A1002
        Complex*16 B12(6,11)
        Real*8 Inv2,Inv3,Inv5,Inv6,Inv8,Inv12,Inv18,Inv20,Inv24,Inv36
     &         ,Inv48,Inv60
        Real*8 Inv72,Inv192,Inv300,Inv360,Inv600,Inv1344,Inv1440,Inv1470
     &         ,Inv1920
        Real*8 Inv3600,Inv4536,Inv5880,Inv7920,Inv10890,Inv23040
     &         ,Inv23520,Inv26880
        Real*8 Inv45360,Inv53760,Inv100800,Inv190080,Inv201600,Inv544320
     &         ,Inv609840
        Real*8 Inv1088640,Inv2419200,Inv2439360,Inv7983360,Inv14636160,
     &         Inv58544640
        Real*8 Inv63866880,Inv127733760
        Real*8 p1sq2,p1sq3,p1sq4,p1sq5,m0sq2,m0sq3,m0sq4,m0sq5,m1sq2,
     &         m1sq3,m1sq4,m1sq5
        Real*8 Invp1sq
        Real*8 eps1
        real*8 Inv22680,Inv25200,Inv304920,Inv332640,Inv560,Inv980
c-------the size of eps1 should still be checked!
        Parameter(eps1=1d-7)
        integer I

c 	------ variables for numeric Calculation of B12 -----
        complex*16 x1,x2,f_n,epsilon
        parameter(epsilon=(0,1d-38)) 

c FC %       index(1)=1
c FC %       index(2)=2
c FC %       index(3)=2
c FC %       index(4)=3
c FC %       index(5)=3
c FC %       index(6)=4
c FC %       index(7)=4
c FC %       index(8)=5
c FC %       index(9)=5
c FC %       index(10)=6
c FC %       index(11)=6

!------ If is not needed anymore with Denner/Dittmaier method
!       If(abs(p1sq).gt.1d-8) then
        Inv980=1d0/980d0
       Inv560=1d0/560d0
       Inv25200=1d0/25200d0
       Inv304920  =1D0/304920d0
       Inv22680=1d0/22680d0 
       Inv332640    =    1d0/332640d0


        Inv2=1d0/2d0
	Inv3=1d0/3d0
        Inv5=1d0/5d0
        Inv6=1d0/6d0
        Inv8=1d0/8d0
        Inv12=1d0/12d0
        Inv18=1d0/18d0
        Inv20=1d0/20d0
        Inv24=1d0/24d0
        Inv36=1d0/36d0
        Inv48=1d0/48d0
        Inv60=1d0/60d0
        Inv72=1d0/72d0
        Inv192=1d0/192d0
        Inv300=1d0/300d0
        Inv360=1d0/360d0
        Inv600=1d0/600d0
        Inv1344=1d0/1344d0
        Inv1440=1d0/1440d0
        Inv1470=1d0/1470d0
        Inv1920=1d0/1920d0
        Inv3600=1d0/3600d0
        Inv4536=1d0/4536d0
        Inv5880=1d0/5880d0
        Inv7920=1d0/7920d0
        Inv10890=1d0/10890d0
        Inv23040=1d0/23040d0
        Inv23520=1d0/23520d0
        Inv26880=1d0/26880d0
      
        Inv45360=1d0/45360d0
        Inv53760=1d0/53760d0
        Inv100800=1d0/100800d0
        Inv190080=1d0/190080d0
        Inv201600=1d0/201600d0
        Inv544320=1d0/544320d0
        Inv609840=1d0/609840d0
        Inv1088640=1d0/1088640d0
        Inv2419200=1d0/2419200d0
        Inv2439360=1d0/2439360d0
        Inv7983360=1d0/7983360d0
        Inv14636160=1d0/14636160d0
        Inv58544640=1d0/58544640d0
        Inv63866880=1d0/63866880d0
        Inv127733760=1d0/127733760d0
        p1sq2=p1sq*p1sq
        p1sq3=p1sq2*p1sq
        p1sq4=p1sq3*p1sq
        p1sq5=p1sq4*p1sq
        m1sq=m1*m1
        m1sq2=m1sq*m1sq
        m1sq3=m1sq2*m1sq
        m1sq4=m1sq3*m1sq
        m1sq5=m1sq4*m1sq
        m0sq=m0*m0
        m0sq2=m0sq*m0sq
        m0sq3=m0sq2*m0sq
        m0sq4=m0sq3*m0sq
        m0sq5=m0sq4*m0sq
       if(abs(p1sq).gt.eps1) then
        Invp1sq=1d0/p1sq        
       endif
       A202=Inv8*m1sq*(2*A2+m1sq)
       A402=Inv24*(A2+5*Inv6*m1sq)*m1sq2
       A602=Inv192*(A2+13*Inv12*m1sq)*m1sq3
       A802=Inv1920*(A2+77*Inv60*m1sq)*m1sq4
       A1002=Inv23040*(A2+29*Inv20*m1sq)*m1sq5


       if(abs(p1sq).gt.eps1) then
          if (m0.gt.eps1) then
             if(m1.gt.eps1) then
c               B0(p1sq,m0,m1)                
                if(abs(p1sq-m0sq).lt.eps1) then
                   if(abs(p1sq-m1sq).lt.eps1) then
c                       B0(m0sq,m0sq,m0sq)
                      goto 101
                   else
c                        B0(m0sq,m0sq,m1sq)
                      goto 101
                   endif
                else   
                   if(abs(p1sq-m1sq).lt.eps1) then
c                       B0(m1sq,m0sq,m1sq)
                      goto 101
                   else
                      if(abs(m1-m0).lt.eps1) then
c                        B0(p1sq,m0sq,m0sq)                         
                         goto 101
                      else
c                        B0(p1sq,m0sq,m1sq)
                         goto 101
                      endif   
                   endif
                endif   
c  m1 =0          B(p1sq,m0,0)     
             else 
               if(abs(p1sq-m0sq).lt.eps1) then
c                  B0(m0sq,m0sq,0)
                  goto 101
               else
c                 B0(p1sq,m0,0) 
                  goto 101
               endif
             endif
c  m0=0    
          else  
             if(m1.gt.eps1) then
c          B0(p1sq,0,m1)
               if(abs(p1sq-m1sq).lt.eps1) then
c                  B0(m1sq,0,m1sq)
                  goto 103
               else
c                 B0(p1sq,0,m1sq) 
                  goto 102
               endif
             else
c         B0(p1sq,0,0)
                goto 104
             endif   
          endif
        else
c p1sq =0  
          if (m0.gt.eps1) then
             if(m1.gt.eps1) then
c               B0(0,m0,m1)                
                if(abs(m1sq-m0sq).lt.eps1) then
c                       B0(0,m0sq,m0sq)
                      goto 106
                else   
c                       B0(0,m0sq,m1sq)
                     goto 105
                endif   
c  m1 =0          
             else 
c                 B0(0,m0,0) 
                goto 105
             endif
c  m0=0    
          else  
             if(m1.gt.eps1) then
c          B0(0,0,m1)
              goto 107  
             else
c         B0(0,0,0)
              goto 108
             endif   
          endif
         endif

c  General Case
 101      call solve_xk(p1sq,m0,m1,x1,x2) 
c          write(*,*) 'General case 101',p1sq, m0, m1,musq
           do I=1,11
           B12(1,I)=(log(musq/(m0sq))-f_n(I,x1)-f_n(I,x2))
     &     /(I+1)*(-1)**I        
           end do
           goto 100
c B0(p1sq,0,m1sq)
 102      do I=1,11
            B12(1,I)=(log(musq/(m1sq-p1sq-epsilon))+1d0/(dble(I)+1d0)
     &       -f_n(I,1d0-(m1sq-epsilon)*invp1sq))/(dble(I)+1d0)*(-1)**I
             end do
c             write(*,*) '!!m0 is small 102!!'
           goto 100
c B0(m1sq,0,m1sq)
 103        do I=1,11
             B12(1,I)=(log(musq/m1sq)+2d0/(I+1d0))/(I+1d0)*(-1)**I
             end do
c            write(*,*) '!!m0/m1 and p1sq-m1sq is small 103!!'
           goto 100
c  B0(p1sq,0,0)
 104   B12(1,1)=-(B012*Inv2)
       B12(1,2)=(1d0+6d0*B012)*Inv18
       B12(1,3)=(-1d0-3d0*B012)*Inv12
       B12(1,4)=29d0*Inv300+B012*Inv5
       B12(1,5)=-Inv360*(60d0*B012+37d0)
       B12(1,6)= Inv980*(103d0+140d0*B012)
       B12(1,7)=-Inv560*(59d0+70d0*B012)
       B12(1,8)=Inv22680*(2369d0+2520d0*B012)
       B12(1,9)=(-2593d0 - 2520d0*B012)*Inv25200
       B12(1,10)=(30791d0 + 27720d0*B012)*Inv304920
       B12(1,11)=(-32891d0 - 27720d0*B012)*Inv332640
c            write(*,*) '!!m0 and m1 are small 104!!'
           goto 100
c B0(0,m0,m1) 
 105   do I=1,11
          B12(1,I)=(log(musq/(m0sq))
     &    -f_n(I,(m0sq-epsilon)/(m0sq-m1sq)))/(dble(I)+1d0)*(-1)**I
       end do          
c           write(*,*) '105 Small p1sq! p1sq/m0sq = ',p1sq/m0sq,p1sq,m0sq,musq
           goto 100
c B0(0,m0,m0)       
 106    do I=1,11
c             print*, "I",I
          B12(1,I)=log(musq/(m0sq))/(dble(I)+1d0)*(-1)**I
          end do
c          print*,"106 m0 equal m1,p1sq=0", m0,m1, p1sq ,musq
           goto 100
c B0(0,0,m1)
 107      print*, "NOT IMPLEMENTED YET, SET TO ZERO"
 108  B12(1,1)=0d0
       B12(2,2)=0d0
       B12(1,2)=0d0
       B12(2,3)=0d0
       B12(1,3)=0d0
       B12(3,4)=0d0
       B12(2,4)=0d0
       B12(1,4)=0d0
       B12(3,5)=0d0
       B12(2,5)=0d0
       B12(1,5)=0d0
       B12(4,6)=0d0
       B12(3,6)=0d0
       B12(2,6)=0d0
       B12(1,6)=0d0
       B12(4,7)=0d0
       B12(3,7)=0d0
       B12(2,7)=0d0
       B12(1,7)=0d0
       B12(5,8)=0d0
       B12(4,8)=0d0
       B12(3,8)=0d0
       B12(2,8)=0d0
       B12(1,8)=0d0
       B12(5,9)=0d0
       B12(4,9)=0d0
       B12(3,9)=0d0
       B12(2,9)=0d0
       B12(1,9)=0d0
       B12(6,10)=0d0
       B12(5,10)=0d0
       B12(4,10)=0d0
       B12(3,10)=0d0
       B12(2,10)=0d0
       B12(1,10)=0d0
       B12(6,11)=0d0
       B12(5,11)=0d0
       B12(4,11)=0d0
       B12(3,11)=0d0
       B12(2,11)=0d0
       B12(1,11)=0d0

       return


cFCc------new
cFCc   B0(0,m0,m0)
cFC       if((m0.gt.eps1).and.abs(m0/m1-1d0).lt.eps1.and.p1sq.lt.eps1) then
cFCc     
cFC        print*,"m0 equal m1", m0,m1, p1sq ,musq
cFC          do I=1,11
cFCc             print*, "I",I
cFC          B12_New(1,I)=log(musq/(m0sq))/(dble(I)+1d0)*(-1)**I
cFC          end do
cFCc B0(0,m0,m2)
cFCc------for small values of p1sq:
cFC       elseif (m0.gt.eps1 .AND. (abs(p1sq/m0sq).lt.eps1 .OR. abs(p1sq/m1sq).lt.eps1)) then
cFC          write(*,*) 'Small p1sq! p1sq/m0sq = ',p1sq/m0sq,p1sq,m0sq,musq
cFC          do I=1,11
cFCc             print*, "I",I
cFC          B12_New(1,I)=(log(musq/(m0sq))
cFC     &    -f_n(I,(m0sq-epsilon)/(m0sq-m1sq)))/(dble(I)+1d0)*(-1)**I
cFC          end do
cFC       else
cFCc---------  m0/m1? or different choice?
cFC!          if ( (m0/p1sq).lt.eps1 .OR. (m0/m1).lt.eps1) then
cFCcB0(p1sq?,0,m1)
cFC          if ( m0.lt.eps1) then
cFC            write(*,*) '!!m0 is small!!'
cFCc--------- case: p1sq==m1sq
cFCc B0(m1sq,0,m1sq)
cFC            if (abs(p1sq-m1sq).lt.eps1) then
cFC            write(*,*) '!!m0/m1 and p1sq-m1sq is small!!'
cFC             do I=1,11
cFC             B12_New(1,I)=(log(musq/m1sq)+2d0/(I+1d0))/(I+1d0)*(-1)**I
cFC             end do
cFCc B0(p1sq,0,m1sq)
cFC            else
cFC            write(*,*) '!!m0/m1 is small!!'
cFC             do I=1,11
cFC             B12_New(1,I)=(log(musq/(m1sq-p1sq-epsilon))+1d0/(dble(I)+1d0)
cFC     &       -f_n(I,1d0-(m1sq-epsilon)*invp1sq))/(dble(I)+1d0)*(-1)**I
cFC             end do
cFC            end if
cFC          else
cFC          write(*,*) 'General case'
cFC           call solve_xk(p1sq,m0,m1,x1,x2) 
cFC           do I=1,11
cFC           B12_New(1,I)=(log(musq/(m0sq))-f_n(I,x1)-f_n(I,x2))
cFC     &     /(I+1)*(-1)**I        
cFC           end do
cFC          end if
cFC       end if

 100   B12(2,2)=Inv18*(3*A2+6*B012*m0sq+3*(m0sq+m1sq)-p1sq+3*(m0sq-m1sq+
     &  p1sq)*B12(1,1))
       B12(2,3)=Inv48*(-6*A2-2*m0sq-4*m1sq+p1sq+12*m0sq*B12(1,1)+6*(m0sq
     &  -m1sq+p1sq)*B12(1,2))
       B12(3,4)=Inv600*(10*(m0sq2+m0sq*m1sq+m1sq2)-5*(m0sq+m1sq)*p1sq+p1
     &  sq2+60*(A202+2*m0sq*B12(2,2)+(m0sq-m1sq+p1sq)*B12(2,3)))
       B12(2,4)=Inv300*(5*m0sq+15*m1sq-3*p1sq+30*(A2+2*m0sq*B12(1,2)+(m0
     &  sq-m1sq+p1sq)*B12(1,3)))
       B12(3,5)=Inv1440*(-5*(m0sq2+2*m0sq*m1sq+3*m1sq2)+4*m0sq*p1sq+6*m1
     &  sq*p1sq-p1sq2+120*(-A202+2*m0sq*B12(2,3)+(m0sq-m1sq+p1sq)*B12(2,
     &  4)))
       B12(2,5)=Inv360*(-30*A2-3*(m0sq+4*m1sq)+2*p1sq+60*m0sq*B12(1,3)+3
     &  0*(m0sq-m1sq+p1sq)*B12(1,4))
       B12(4,6)=Inv23520*(1680*A402+35*(m0sq+m1sq)*(m0sq2+m1sq2)-7*(3*m0
     &  sq2+4*m0sq*m1sq+3*m1sq2)*p1sq+7*(m0sq+m1sq)*p1sq2-p1sq3+1680*(2*
     &  m0sq*B12(3,4)+(m0sq-m1sq+p1sq)*B12(3,5)))
       B12(3,6)=Inv5880*(7*(m0sq2+3*m0sq*m1sq+6*m1sq2)-7*(m0sq+2*m1sq)*p
     &  1sq+2*p1sq2+420*(A202+2*m0sq*B12(2,4)+(m0sq-m1sq+p1sq)*B12(2,5))
     &  )
       B12(2,6)=Inv1470*(7*(m0sq+5*m1sq)-5*p1sq+105*(A2+2*m0sq*B12(1,4)+
     &  (m0sq-m1sq+p1sq)*B12(1,5)))
       B12(4,7)=Inv53760*(-14*(m0sq3+2*m0sq2*m1sq+3*m0sq*m1sq2+4*m1sq3)+
     &  14*(m0sq2+2*m0sq*m1sq+2*m1sq2)*p1sq-2*(3*m0sq+4*m1sq)*p1sq2+p1sq
     &  3+3360*(-A402+2*m0sq*B12(3,5)+(m0sq-m1sq+p1sq)*B12(3,6)))
       B12(3,7)=Inv26880*(-14*(m0sq2+4*m0sq*m1sq+10*m1sq2)+8*(2*m0sq+5*m
     &  1sq)*p1sq-5*p1sq2+1680*(-A202+2*m0sq*B12(2,5)+(m0sq-m1sq+p1sq)*B
     &  12(2,6)))
       B12(2,7)=Inv1344*(-4*(m0sq+6*m1sq)+3*p1sq+84*(-A2+2*m0sq*B12(1,5)
     &  +(m0sq-m1sq+p1sq)*B12(1,6)))
       B12(5,8)=Inv1088640*(126*(m0sq4+m0sq3*m1sq+m0sq2*m1sq2+m0sq*m1sq3
     &  +m1sq4)-42*(m0sq+m1sq)*(2*m0sq2+m0sq*m1sq+2*m1sq2)*p1sq+18*(2*m0
     &  sq2+3*m0sq*m1sq+2*m1sq2)*p1sq2-9*(m0sq+m1sq)*p1sq3+p1sq4+60480*(
     &  A602+2*m0sq*B12(4,6)+(m0sq-m1sq+p1sq)*B12(4,7)))
       B12(4,8)=Inv544320*(42*(m0sq3+3*m0sq2*m1sq+6*m0sq*m1sq2+10*m1sq3)
     &  -18*(3*m0sq2+8*m0sq*m1sq+10*m1sq2)*p1sq+9*(3*m0sq+5*m1sq)*p1sq2-
     &  5*p1sq3+30240*(A402+2*m0sq*B12(3,6)+(m0sq-m1sq+p1sq)*B12(3,7)))
       B12(3,8)=Inv45360*(12*(m0sq2+5*m0sq*m1sq+15*m1sq2)-15*(m0sq+3*m1s
     &  q)*p1sq+5*p1sq2+2520*(A202+2*m0sq*B12(2,6)+(m0sq-m1sq+p1sq)*B12(
     &  2,7)))
       B12(2,8)=Inv4536*(9*(m0sq+7*m1sq)-7*p1sq+252*(A2+2*m0sq*B12(1,6)+
     &  (m0sq-m1sq+p1sq)*B12(1,7)))
       B12(5,9)=Inv2419200*(-42*(m0sq4+2*m0sq3*m1sq+3*m0sq2*m1sq2+4*m0sq
     &  *m1sq3+5*m1sq4)+12*(4*m0sq3+9*m0sq2*m1sq+12*m0sq*m1sq2+10*m1sq3)
     &  *p1sq-9*(3*m0sq2+6*m0sq*m1sq+5*m1sq2)*p1sq2+2*(4*m0sq+5*m1sq)*p1
     &  sq3-p1sq4+120960*(-A602+2*m0sq*B12(4,7)+(m0sq-m1sq+p1sq)*B12(4,8
     &  )))
       B12(4,9)=Inv201600*(-6*(m0sq3+4*m0sq2*m1sq+10*m0sq*m1sq2+20*m1sq3
     &  )+3*(3*m0sq2+10*m0sq*m1sq+15*m1sq2)*p1sq-5*(m0sq+2*m1sq)*p1sq2+p
     &  1sq3+10080*(-A402+2*m0sq*B12(3,7)+(m0sq-m1sq+p1sq)*B12(3,8)))
       B12(3,9)=Inv100800*(-15*(m0sq2+6*m0sq*m1sq+21*m1sq2)+10*(2*m0sq+7
     &  *m1sq)*p1sq-7*p1sq2+5040*(-A202+2*m0sq*B12(2,7)+(m0sq-m1sq+p1sq)
     &  *B12(2,8)))
       B12(2,9)=Inv3600*(-5*(m0sq+8*m1sq)+4*p1sq+180*(-A2+2*m0sq*B12(1,7
     &  )+(m0sq-m1sq+p1sq)*B12(1,8)))
       B12(6,10)=Inv58544640*(2661120*A802+462*(m0sq+m1sq)*(m0sq2-m0sq*m
     &  1sq+m1sq2)*(m0sq2+m0sq*m1sq+m1sq2)-66*(5*m0sq4+8*m0sq3*m1sq+9*m0
     &  sq2*m1sq2+8*m0sq*m1sq3+5*m1sq4)*p1sq+33*(m0sq+m1sq)*(5*m0sq2+4*m
     &  0sq*m1sq+5*m1sq2)*p1sq2-11*(5*m0sq2+8*m0sq*m1sq+5*m1sq2)*p1sq3+1
     &  1*(m0sq+m1sq)*p1sq4-p1sq5+2661120*(2*m0sq*B12(5,8)+(m0sq-m1sq+p1
     &  sq)*B12(5,9)))
       B12(5,10)=Inv14636160*(66*(m0sq4+3*m0sq3*m1sq+6*m0sq2*m1sq2+10*m0
     &  sq*m1sq3+15*m1sq4)-99*(m0sq3+3*m0sq2*m1sq+5*m0sq*m1sq2+5*m1sq3)*
     &  p1sq+33*(2*m0sq2+5*m0sq*m1sq+5*m1sq2)*p1sq2-11*(2*m0sq+3*m1sq)*p
     &  1sq3+3*p1sq4+665280*(A602+2*m0sq*B12(4,8)+(m0sq-m1sq+p1sq)*B12(4
     &  ,9)))
       B12(4,10)=Inv2439360*(33*(m0sq3+5*m0sq2*m1sq+15*m0sq*m1sq2+35*m1s
     &  q3)-55*(m0sq2+4*m0sq*m1sq+7*m1sq2)*p1sq+11*(3*m0sq+7*m1sq)*p1sq2
     &  -7*p1sq3+110880*(A402+2*m0sq*B12(3,8)+(m0sq-m1sq+p1sq)*B12(3,9))
     &  )
       B12(3,10)=Inv609840*(55*(m0sq2+7*m0sq*m1sq+28*m1sq2)-77*(m0sq+4*m
     &  1sq)*p1sq+28*p1sq2+27720*(A202+2*m0sq*B12(2,8)+(m0sq-m1sq+p1sq)*
     &  B12(2,9)))
       B12(2,10)=Inv10890*(11*(m0sq+9*m1sq)-9*p1sq+495*(A2+2*m0sq*B12(1,
     &  8)+(m0sq-m1sq+p1sq)*B12(1,9)))
       B12(6,11)=Inv127733760*(-132*(m0sq5+2*m0sq4*m1sq+3*m0sq3*m1sq2+4*
     &  m0sq2*m1sq3+5*m0sq*m1sq4+6*m1sq5)+33*(5*m0sq4+12*m0sq3*m1sq+18*m
     &  0sq2*m1sq2+20*m0sq*m1sq3+15*m1sq4)*p1sq-22*(5*m0sq3+12*m0sq2*m1s
     &  q+15*m0sq*m1sq2+10*m1sq3)*p1sq2+22*(2*m0sq2+4*m0sq*m1sq+3*m1sq2)
     &  *p1sq3-2*(5*m0sq+6*m1sq)*p1sq4+p1sq5+5322240*(-A802+2*m0sq*B12(5
     &  ,9)+(m0sq-m1sq+p1sq)*B12(5,10)))
       B12(5,11)=Inv63866880*(-99*(m0sq4+4*m0sq3*m1sq+10*m0sq2*m1sq2+20*
     &  m0sq*m1sq3+35*m1sq4)+44*(4*m0sq3+15*m0sq2*m1sq+30*m0sq*m1sq2+35*
     &  m1sq3)*p1sq-66*(2*m0sq2+6*m0sq*m1sq+7*m1sq2)*p1sq2+12*(4*m0sq+7*
     &  m1sq)*p1sq3-7*p1sq4+2661120*(-A602+2*m0sq*B12(4,9)+(m0sq-m1sq+p1
     &  sq)*B12(4,10)))
       B12(4,11)=Inv7983360*(-55*(m0sq3+6*m0sq2*m1sq+21*m0sq*m1sq2+56*m1
     &  sq3)+33*(3*m0sq2+14*m0sq*m1sq+28*m1sq2)*p1sq-21*(3*m0sq+8*m1sq)*
     &  p1sq2+14*p1sq3+332640*(-A402+2*m0sq*B12(3,9)+(m0sq-m1sq+p1sq)*B1
     &  2(3,10)))
       B12(3,11)=Inv190080*(-11*(m0sq2+8*m0sq*m1sq+36*m1sq2)+8*(2*m0sq+9
     &  *m1sq)*p1sq-6*p1sq2+7920*(-A202+2*m0sq*B12(2,9)+(m0sq-m1sq+p1sq)
     &  *B12(2,10)))
       B12(2,11)=Inv7920*(-6*(m0sq+10*m1sq)+5*p1sq+330*(-A2+2*m0sq*B12(1
     &  ,9)+(m0sq-m1sq+p1sq)*B12(1,10)))

c FC %       DO I1=1,11
c FC %          print*,''
c FC %          print*, 'Rank=', i1
c FC %              DO I2=1,INDEX(I1)
c FC %                 print*,'B12(I2,I1)',B12(I2,I1)
c FC %              enddo
c FC %       enddo 


!       else
!       print*,'not yet implemented'
!       stop
!       Endif
       return
       End


c     ---------------------------------------------------------------
c     find solutions to the quadratic equation in ln(...) in B_1...1
c     ---------------------------------------------------------------
      subroutine solve_xk(p1sq,m0,m1,x1,x2)
      implicit none

      complex*16  epsilon
      parameter(epsilon=(0,1d-38)) 

      real*8  p1sq,m0,m1,m0sq,m1sq
      complex*16 dsqrt
      complex*16 i1
      parameter (i1=(0.0d0,1.0d0))
      complex*16 x1,x2
      real*8 inv2,invp1,eps1
      parameter(eps1=1d-7)

      inv2=1d0/2d0
      invp1=1d0/p1sq
      m0sq=m0*m0
      m1sq=m1*m1

      dsqrt=sqrt( (p1sq+m0sq-m1sq)*(p1sq+m0sq-m1sq)
     & -4*p1sq*m0sq+4*p1sq*epsilon )

      if (abs(dsqrt).lt.abs(epsilon) .OR. abs(m0).lt.abs(epsilon)) then
         x1=1d0-(m1sq-epsilon)*invp1
         x2=(0d0,0d0)
      else
         x1=( (p1sq+m0sq-m1sq) + dsqrt )*invp1*inv2
         x2=( (p1sq+m0sq-m1sq) - dsqrt )*invp1*inv2
      end if
      return
      end subroutine solve_xk

c     ---------------------------------------------------------------
c     calculate the auxiliary function f_n(x_k) for B12_New
c     ---------------------------------------------------------------
      complex*16 function f_n(n,x)
      implicit none

      integer l,n
      complex*16 sum,test,x

      l=0
      sum=(0d0,0d0)
      test=(0d0,0d0)
c      print*, "n",n,"x",x
c      if (abs(x)<=10d0) then
      if (abs(x).le.10d0) then
            do while (l.le.n)
            sum=sum+(x**(n-l))/dble(l+1)
            l=l+1
            end do
            f_n=(1d0-x**(n+1))*log( (x-1)/x )-sum
c            test2=(1d0-x**(n+1))*log( (x-1)/x )-sum

      else 
            f_n=(0d0,0d0)
            if (abs(x).gt.1d20)  return
            l=n+1
            do
            sum=sum+x**(n-l)/dble(l+1)
            l=l+1
c            if (test==sum) exit
            if (abs(test/sum-1d0).lt.1d-15) exit
            test=sum
            end do
            f_n=log(1-1/x)+sum
      end if
      return
      end function

