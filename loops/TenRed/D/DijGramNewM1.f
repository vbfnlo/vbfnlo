      subroutine ten_red4_Gram1_M(m0,p1sq,p2sq,p3sq,p4sq
     &  ,s12,s23,musq,C0234,C0134,C0124,C0123,D0,Dij,order)
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
      Implicit None
      Real*8 p1sq,p2sq,p3sq,p4sq,s12,s23,musq
      Real*8 p1sqsq,p2sqsq,p3sqsq,p4sqsq,s23sq,s12sq
      Real*8 r10,r21,r32,m0,m0sq,m1sq,m2sq,m3sq
      Complex*16 D0,Dij(34,5)
      Complex*16 C0123,C0124,C0134,C0234
      Real*8 C0123R,C0124R,C0134R,C0234R
      Real*8 C0123I,C0124I,C0134I,C0234I
      Complex*16 Cij123(30,9),Cij124(30,9),Cij134(30,9),Cij234(30,9)
      Complex*16 S31(3)
      Complex*16 S3h001(3)
      Complex*16 S300,S3001(3)
      Complex*16 S30000,S300001(3)
      Complex*16 S3h00001(3)
      Complex*16 S3000000,S30000001(3),S3h0000001(3)
      Complex*16 auxD40,tempD40
      Complex*16 auxD400,tempD400
      Complex*16 auxD40000,tempD40000
      Complex*16 auxD4000000,tempD4000000
      Complex*16 aux1(3),temp1(3)
      Complex*16 aux2(3,3),temp2(3,3)
      Complex*16 aux3(3,3,3),temp3(3,3,3)
      Complex*16 aux001(3),temp001(3)
      Complex*16 aux002(3,3),temp002(3,3)
      Complex*16 aux003(3,3,3),temp003(3,3,3)
      Complex*16 aux00001(3),temp00001(3)
      Complex*16 aux00002(3,3),temp00002(3,3)
      Complex*16 aux00003(3,3,3),temp00003(3,3,3)
      Complex*16 aux0000001(3),temp0000001(3)
      Real*8 Z(3,3),ZZ(3,3,3,3),IX,ZMax,I4Z,I8Z,I12Z,I16Z,I20Z,I24Z,I28Z,I32Z
      Real*8 TX1,TX2,TX3
      Integer k,l,jj
      Real*8 F(7),det4,F1
      Real*8 Inv40,Inv60,Inv80
           include 'DijGramM.inc'
      integer jjtemp,ktemp,ltemp,jjinit,kinit,linit,cont
      real*8 tempjj,tempkl,IXinit,Zmaxinit
      real*8 tempjj1,tempjj2,tempjj3
      real*8 tempkl1,tempkl2,tempkl3,tempkl4,tempkl5,tempkl6
      real*8 IXtemp,Zmaxtemp
      Common/DDecide/tempjj,tempkl,IX,Zmax,jjinit
      Save/DDecide/
      real*8 tempCjj,tempCkl,IXC,ZMaxC,ratio
      Integer jjinitC
      Common/Decide/tempCjj,tempCkl,IXC,ZMaxC,jjinitC
      Save/Decide/
      real*8 det4abs
      integer order
c
      real*8 accuracyC(0:4,5),AccuracyD(0:5,4)
      real*8 accuracyDR(22,0:4,4)
      Common/Accuracy/AccuracyC,AccuracyD
      integer i1,i2,index(0:4),ac
       logical printmy
       common/mprint/printmy
 
      index(0)=1
      index(1)=3
      index(2)=7
      index(3)=13
      index(4)=22
   
      m0sq=m0*m0
      m1sq=m0sq
      m2sq=m0sq
      m3sq=m0sq

       cont=0
c############
      call dt3(p1sq,p2sq,s12,ratio)
c      print*, 'ratio1',ratio
      if(ratio.lt.1d-2) then
        call MyCgetGram_M(m0,p1sq,p2sq,s12,musq,C0123R,C0123I,Cij123)
        C0123=DCMPLX(C0123R,C0123I)
c             if((tempCkl.ge.2d-3).or.(tempCjj.ge.1d-2))then
c                    call MyCget(p1sq,p2sq,s12,musq,C0123R,C0123I,Cij123)
c                  C0123=DCMPLX(C0123R,C0123I)
c              endif
      else
c         call MyCget_QUAD_D1_M(m0,p1sq,p2sq,s12,musq,C0123R,C0123I,Cij123)
         call MyCget_M(m0,p1sq,p2sq,s12,musq,C0123R,C0123I,Cij123)
         C0123=DCMPLX(C0123R,C0123I)
      endif
c#########################
      call dt3(p1sq,s23,p4sq,ratio) 
c      print*, 'ratio2',ratio

      if(ratio.lt.1d-2) then
         call MyCgetGram_M(m0,p1sq,s23,p4sq,musq,C0124R,C0124I,Cij124)
         C0124=DCMPLX(C0124R,C0124I)
c                if((tempCkl.ge.2.3d-3).and.(tempCjj.ge.1d-2))then
             !if(tempCkl.ge.2.3d-3) then
c                 call MyCget(p1sq,s23,p4sq,musq,C0124R,C0124I,Cij124)
c                 C0124=DCMPLX(C0124R,C0124I)
c              endif
      else
c          call MyCget_QUAD_D1_M(m0,p1sq,s23,p4sq,musq,C0124R,C0124I,Cij124)
           call MyCget_M(m0,p1sq,s23,p4sq,musq,C0124R,C0124I,Cij124)
          C0124=DCMPLX(C0124R,C0124I)
      endif
c#########################
      call dt3(s12,p3sq,p4sq,ratio)
c      print*, 'ratio3',ratio
      if(ratio.lt.1d-2) then
            call MyCgetGram_M(m0,s12,p3sq,p4sq,musq,C0134R,C0134I,Cij134)
            C0134=DCMPLX(C0134R,C0134I)
c       if((tempCkl.ge.2.3d-3).and.(tempCjj.ge.1d-2))then
                    !if(tempCkl.ge.2.3d-3) then
c           call MyCget(s12,p3sq,p4sq,musq,C0134R,C0134I,Cij134)
c                C0134=DCMPLX(C0134R,C0134I)
c              endif
       else
c         call MyCget_QUAD_D1_M(m0,s12,p3sq,p4sq,musq,C0134R,C0134I,Cij134)
          call MyCget_M(m0,s12,p3sq,p4sq,musq,C0134R,C0134I,Cij134)
         C0134=DCMPLX(C0134R,C0134I)
       endif
c#########################
      call dt3(p2sq,p3sq,s23,ratio)
c      print*, 'ratio4',ratio
      if(ratio.lt.1d-2) then
         call MyCgetGram_M(m0,p2sq,p3sq,s23,musq,C0234R,C0234I,Cij234)
         C0234=DCMPLX(C0234R,C0234I)
c      if((tempCkl.ge.2.3d-3).and.(tempCjj.ge.1d-2))then
              !if(tempCkl.ge.2.3d-3) then
c                 call MyCget(p2sq,p3sq,s23,musq,C0234R,C0234I,Cij234)
c                 C0234=DCMPLX(C0234R,C0234I)
c              endif
       else
c         call MyCget_QUAD_D1_M(m0,p2sq,p3sq,s23,musq,C0234R,C0234I,Cij234)
          call MyCget_M(m0,p2sq,p3sq,s23,musq,C0234R,C0234I,Cij234)
         C0234=DCMPLX(C0234R,C0234I)
       endif  
c#########################
c       Print*, 'HERE'

       p1sqsq=p1sq*p1sq
       p2sqsq=p2sq*p2sq
       p3sqsq=p3sq*p3sq
       p4sqsq=p4sq*p4sq
       s23sq=s23*s23
       s12sq=s12*s12

       Z(1,1)=-p2sqsq-(p3sq-s23)**2+2*p2sq*(p3sq+s23)
       Z(2,1)=-p2sqsq+2*p1sq*p3sq-p3sq*(p4sq+s12)+(p3sq+p4sq-s12)*s23+p2
     &  sq*(p3sq-p4sq+s12+2*s23)-s23sq
       Z(3,1)=-p2sqsq+p1sq*(p2sq+p3sq-s23)+s12*(-p3sq+s23)+p2sq*(p3sq-2*
     &  p4sq+s12+s23)
       Z(1,2)=-p2sqsq+2*p1sq*p3sq-p3sq*(p4sq+s12)+(p3sq+p4sq-s12)*s23+p2
     &  sq*(p3sq-p4sq+s12+2*s23)-s23sq
       Z(2,2)=4*p1sq*p3sq-(p2sq+p4sq-s12-s23)**2
       Z(3,2)=-((p2sq-s12)*(p2sq+p4sq-s12-s23))+p1sq*(p2sq+2*p3sq-p4sq+s
     &  12-s23)
       Z(1,3)=-p2sqsq+p1sq*(p2sq+p3sq-s23)+s12*(-p3sq+s23)+p2sq*(p3sq-2*
     &  p4sq+s12+s23)
       Z(2,3)=-((p2sq-s12)*(p2sq+p4sq-s12-s23))+p1sq*(p2sq+2*p3sq-p4sq+s
     &  12-s23)
       Z(3,3)=-p1sqsq-(p2sq-s12)**2+2*p1sq*(p2sq+s12)
      If(Abs(Z(1,1)).ge.Abs(Z(1,2))) then
               If(Abs(Z(1,1)).ge.Abs(Z(1,3)))then
                       If(Abs(Z(1,1)).ge.Abs(Z(2,2)))then
                               If(Abs(Z(1,1)).ge.Abs(Z(2,3)))then
                                       If(Abs(Z(1,1)).ge.Abs(Z(3,3)))then
                                               k=1
                                               l=1
                                               Zmax=Z(1,1)
                                       else
                                               k=3
                                               l=3
                                               Zmax=Z(3,3)
                                       endif
                               else
                                       If(Abs(Z(2,3)).ge.Abs(Z(3,3)))then
                                               k=3
                                               l=2
                                               Zmax=Z(2,3)
                                       else
                                               k=3
                                               l=3
                                               Zmax=Z(3,3)
                                       endif
                               endif
                       else
                               If(Abs(Z(2,2)).ge.Abs(Z(2,3)))then
                                       If(Abs(Z(2,2)).ge.Abs(Z(3,3)))then
                                               k=2
                                               l=2
                                               Zmax=Z(2,2)
                                       else
                                               k=3
                                               l=3
                                               Zmax=Z(3,3)
                                       endif
                               else
                                       If(Abs(Z(2,3)).ge.Abs(Z(3,3)))then
                                               k=3
                                               l=2
                                               Zmax=Z(2,3)
                                       else
                                               k=3
                                               l=3
                                               Zmax=Z(3,3)
                                       endif
                               endif
                       endif
               else
                       If(Abs(Z(1,3)).ge.Abs(Z(2,2)))then
                               If(Abs(Z(1,3)).ge.Abs(Z(2,3)))then
                                       If(Abs(Z(1,3)).ge.Abs(Z(3,3)))then
                                               k=3
                                               l=1
                                               Zmax=Z(1,3)
                                       else
                                               k=3
                                               l=3
                                               Zmax=Z(3,3)
                                       endif
                               else
                                       If(Abs(Z(2,3)).ge.Abs(Z(3,3)))then
                                               k=3
                                               l=2
                                               Zmax=Z(2,3)
                                       else
                                               k=3
                                               l=3
                                               Zmax=Z(3,3)
                                       endif
                               endif
                       else
                               If(Abs(Z(2,2)).ge.Abs(Z(2,3)))then
                                       If(Abs(Z(2,2)).ge.Abs(Z(3,3)))then
                                               k=2
                                               l=2
                                               Zmax=Z(2,2)
                                        else
                                               k=3
                                               l=3
                                               Zmax=Z(3,3)
                                        endif
                               else
                                        If(Abs(Z(2,3)).ge.Abs(Z(3,3)))then
                                               k=3
                                               l=2
                                               Zmax=Z(2,3)
                                        else
                                                k=3
                                                l=3
                                                Zmax=Z(3,3)
                                        endif
                               endif
                       endif
               endif
        else
               If(Abs(Z(1,2)).ge.Abs(Z(1,3)))then
                       If(Abs(Z(1,2)).ge.Abs(Z(2,2)))then
                               If(Abs(Z(1,2)).ge.Abs(Z(2,3)))then
                                       If(Abs(Z(1,2)).ge.Abs(Z(3,3)))then
                                               k=2
                                               l=1
                                               Zmax=Z(1,2)
                                       else
                                               k=3
                                               l=3
                                               Zmax=Z(3,3)
                                       endif
                               else
                                       If(Abs(Z(2,3)).ge.Abs(Z(3,3)))then
                                               k=3
                                               l=2
                                               Zmax=Z(2,3)
                                       else
                                               k=3
                                               l=3
                                               Zmax=Z(3,3)
                                       endif
                               endif
                       else
                               If(Abs(Z(2,2)).ge.Abs(Z(2,3)))then
                                       If(Abs(Z(2,2)).ge.Abs(Z(3,3)))then
                                               k=2
                                               l=2
                                               Zmax=Z(2,2)
                                       else
                                               k=3
                                               l=3
                                               Zmax=Z(3,3)
                                       endif
                               else
                                       If(Abs(Z(2,3)).ge.Abs(Z(3,3)))then
                                               k=3
                                               l=2
                                               Zmax=Z(2,3)
                                       else
                                               k=3
                                               l=3
                                               Zmax=Z(3,3)
                                       endif
                               endif
                       endif
               else
                       If(Abs(Z(1,3)).ge.Abs(Z(2,2)))then
                               If(Abs(Z(1,3)).ge.Abs(Z(2,3)))then
                                       If(Abs(Z(1,3)).ge.Abs(Z(3,3)))then
                                               k=3
                                               l=1
                                               Zmax=Z(1,3)
                                        else
                                                k=3
                                                l=3
                                                Zmax=Z(3,3)
                                        endif
                               else
                                       If(Abs(Z(2,3)).ge.Abs(Z(3,3)))then
                                               k=3
                                               l=2
                                               Zmax=Z(2,3)
                                       else
                                                k=3
                                                l=3
                                                Zmax=Z(3,3)
                                       endif
                                 endif
                       else
                               If(Abs(Z(2,2)).ge.Abs(Z(2,3)))then
                                       If(Abs(Z(2,2)).ge.Abs(Z(3,3)))then
                                               k=2
                                               l=2
                                               Zmax=Z(2,2)
                                       else
                                               k=3
                                               l=3
                                               Zmax=Z(3,3)
                                       endif
                               else
                                       If(Abs(Z(2,3)).ge.Abs(Z(3,3)))then
                                               k=3
                                               l=2
                                               Zmax=Z(2,3)
                                       else
                                               k=3
                                               l=3
                                               Zmax=Z(3,3)
                                       endif
                               endif
                       endif
               endif
       endif

      r10=p1sq
      r21=s12-p1sq
      r32=p4sq-s12
       TX1=p3sqsq*r10+p2sqsq*(r10+r21+r32)+s23*(-(p4sq*r21)+p1sq*r32+(r2
     &  1-r32)*s12+(r10+r21)*s23)+p3sq*(p4sq*r21-p1sq*(2*r21+r32)+(r21+r
     &  32)*s12-(2*r10+r21)*s23)-p2sq*(p3sq*(2*r10+r21+r32)-p4sq*(r21+2*
     &  r32)+r21*s12+2*(r10+r21)*s23+r32*(p1sq+s12+s23))
       TX2=p3sq*p4sq*r10+p4sqsq*r21+p2sqsq*(r10+r21+r32)+p3sq*r10*s12-2*
     &  p4sq*r21*s12-p4sq*r32*s12+r21*s12sq+r32*s12sq-(p3sq*r10+p4sq*(r1
     &  0+2*r21)-(r10+2*r21+r32)*s12)*s23+p1sq*(-2*p3sq*(r10+2*r21+r32)+
     &  r32*(p4sq-s12+s23))-p2sq*(p3sq*r10-p4sq*(r10+2*r21+r32)+(r10+2*r
     &  21)*s12+2*(r10+r21)*s23+r32*(p1sq+2*s12+s23))+(r10+r21)*s23sq
       TX3=p1sqsq*r32+p2sqsq*(r10+r21+r32)-p1sq*(-(p4sq*r21)+p3sq*(r10+2
     &  *r21)+p2sq*(r10+r21+2*r32)+r21*s12+2*r32*s12)+p1sq*(r10+r21)*s23
     &  +s12*(p3sq*r10-p4sq*r21+(r21+r32)*s12+(-r10+r21)*s23)-p2sq*(p3sq
     &  *r10-p4sq*(2*r10+r21)+(r10+2*(r21+r32))*s12+(r10+r21)*s23)
      If(abs(TX1).ge.abs(TX2))then
           If(abs(TX1).ge.abs(TX3))then
              jj=1
              IX=1d0/TX1
           else
               jj=3
               IX=1d0/TX3
           endif
       else
           If(abs(TX2).ge.abs(TX3))then
              jj=2
              IX=1d0/TX2
           else
              jj=3
              IX=1d0/TX3
           endif
      endif


       kinit=k
       linit=l
       jjinit=jj

       Zmaxinit=Zmax
       IXinit=IX

c$$$      Print*, 'IX',IX
c$$$      print*, 'IZax4', I4Z
       ZZ(1,1,1,1)=0
       ZZ(2,1,1,1)=0
       ZZ(3,1,1,1)=0
       ZZ(1,2,1,1)=0
       ZZ(2,2,1,1)=0
       ZZ(3,2,1,1)=0
       ZZ(1,3,1,1)=0
       ZZ(2,3,1,1)=0
       ZZ(3,3,1,1)=0
       ZZ(1,1,2,1)=0
       ZZ(2,1,2,1)=-2*p3sq
       ZZ(3,1,2,1)=-p2sq-p3sq+s23
       ZZ(1,2,2,1)=2*p3sq
       ZZ(2,2,2,1)=0
       ZZ(3,2,2,1)=-p2sq-p4sq+s12+s23
       ZZ(1,3,2,1)=p2sq+p3sq-s23
       ZZ(2,3,2,1)=p2sq+p4sq-s12-s23
       ZZ(3,3,2,1)=0
       ZZ(1,1,3,1)=0
       ZZ(2,1,3,1)=-p2sq-p3sq+s23
       ZZ(3,1,3,1)=-2*p2sq
       ZZ(1,2,3,1)=p2sq+p3sq-s23
       ZZ(2,2,3,1)=0
       ZZ(3,2,3,1)=-p1sq-p2sq+s12
       ZZ(1,3,3,1)=2*p2sq
       ZZ(2,3,3,1)=p1sq+p2sq-s12
       ZZ(3,3,3,1)=0
       ZZ(1,1,1,2)=0
       ZZ(2,1,1,2)=2*p3sq
       ZZ(3,1,1,2)=p2sq+p3sq-s23
       ZZ(1,2,1,2)=-2*p3sq
       ZZ(2,2,1,2)=0
       ZZ(3,2,1,2)=p2sq+p4sq-s12-s23
       ZZ(1,3,1,2)=-p2sq-p3sq+s23
       ZZ(2,3,1,2)=-p2sq-p4sq+s12+s23
       ZZ(3,3,1,2)=0
       ZZ(1,1,2,2)=0
       ZZ(2,1,2,2)=0
       ZZ(3,1,2,2)=0
       ZZ(1,2,2,2)=0
       ZZ(2,2,2,2)=0
       ZZ(3,2,2,2)=0
       ZZ(1,3,2,2)=0
       ZZ(2,3,2,2)=0
       ZZ(3,3,2,2)=0
       ZZ(1,1,3,2)=0
       ZZ(2,1,3,2)=-p2sq-p4sq+s12+s23
       ZZ(3,1,3,2)=-p1sq-p2sq+s12
       ZZ(1,2,3,2)=p2sq+p4sq-s12-s23
       ZZ(2,2,3,2)=0
       ZZ(3,2,3,2)=-2*p1sq
       ZZ(1,3,3,2)=p1sq+p2sq-s12
       ZZ(2,3,3,2)=2*p1sq
       ZZ(3,3,3,2)=0
       ZZ(1,1,1,3)=0
       ZZ(2,1,1,3)=p2sq+p3sq-s23
       ZZ(3,1,1,3)=2*p2sq
       ZZ(1,2,1,3)=-p2sq-p3sq+s23
       ZZ(2,2,1,3)=0
       ZZ(3,2,1,3)=p1sq+p2sq-s12
       ZZ(1,3,1,3)=-2*p2sq
       ZZ(2,3,1,3)=-p1sq-p2sq+s12
       ZZ(3,3,1,3)=0
       ZZ(1,1,2,3)=0
       ZZ(2,1,2,3)=p2sq+p4sq-s12-s23
       ZZ(3,1,2,3)=p1sq+p2sq-s12
       ZZ(1,2,2,3)=-p2sq-p4sq+s12+s23
       ZZ(2,2,2,3)=0
       ZZ(3,2,2,3)=2*p1sq
       ZZ(1,3,2,3)=-p1sq-p2sq+s12
       ZZ(2,3,2,3)=-2*p1sq
       ZZ(3,3,2,3)=0
       ZZ(1,1,3,3)=0
       ZZ(2,1,3,3)=0
       ZZ(3,1,3,3)=0
       ZZ(1,2,3,3)=0
       ZZ(2,2,3,3)=0
       ZZ(3,2,3,3)=0
       ZZ(1,3,3,3)=0
       ZZ(2,3,3,3)=0
       ZZ(3,3,3,3)=0
       Inv40=1d0/6d0
       Inv60=(5*(m0sq+m1sq+m2sq+m3sq)-p1sq-p2sq-p3sq-p4sq-s12-s23)/240.d
     &  0
       Inv80=(42*(m0sq**2+m1sq**2+m2sq**2+m2sq*m3sq+m3sq**2)+p1sq*p3sq+p
     &  2sq*p4sq+s12*s23-14*(m2sq*(p2sq+p3sq+s12)+m3sq*(p3sq+p4sq+s23))+
     &  2*(p1sq**2+p2sq**2+p3sq**2+p3sq*p4sq+p4sq**2+p3sq*s12+p4sq*s12+s
     &  12**2+p3sq*s23+p4sq*s23+s23**2+p2sq*(p3sq+s12+s23)+p1sq*(p2sq+p4
     &  sq+s12+s23))-7*(m2sq*p1sq+m3sq*p1sq+m3sq*p2sq+m2sq*p4sq+m3sq*s12
     &  +m2sq*s23+m0sq*(-6*m1sq-6*m2sq-6*m3sq+2*p1sq+p2sq+p3sq+2*p4sq+2*
     &  s12+s23)+m1sq*(-6*m2sq-6*m3sq+2*p1sq+2*p2sq+p3sq+p4sq+s12+2*s23)
     &  ))/20160.d0
       Inv401=-1d0/8d0
       Inv402=-1d0/12d0
       Inv403=-1d0/24d0
       Inv4011=1d0/10d0
       Inv4021=1d0/15d0
       Inv4022=1d0/20d0
       Inv4031=1d0/30d0
       Inv4032=1d0/40d0
       Inv4033=1d0/60d0
       Inv40111=-1d0/12d0
       Inv40211=-1d0/18d0
       Inv40221=-1d0/24d0
       Inv40222=-1d0/30d0
       Inv40311=-1d0/36d0
       Inv40321=-1d0/48d0
       Inv40322=-1d0/60d0
       Inv40331=-1d0/72d0
       Inv40332=-1d0/90d0
       Inv40333=-1d0/120d0
       Inv401111=1d0/14d0
       Inv402111=1d0/21d0
       Inv402211=1d0/28d0
       Inv402221=1d0/35d0
       Inv402222=1d0/42d0
       Inv403111=1d0/42d0
       Inv403211=1d0/56d0
       Inv403221=1d0/70d0
       Inv403222=1d0/84d0
       Inv403311=1d0/84d0
       Inv403321=1d0/105d0
       Inv403322=1d0/126d0
       Inv403331=1d0/140d0
       Inv403332=1d0/168d0
       Inv403333=1d0/210d0
       Inv4011111=-1d0/16d0
       Inv4021111=-1d0/24d0
       Inv4022111=-1d0/32d0
       Inv4022211=-1d0/40d0
       Inv4022221=-1d0/48d0
       Inv4022222=-1d0/56d0
       Inv4031111=-1d0/48d0
       Inv4032111=-1d0/64d0
       Inv4032211=-1d0/80d0
       Inv4032221=-1d0/96d0
       Inv4032222=-1d0/112d0
       Inv4033111=-1d0/96d0
       Inv4033211=-1d0/120d0
       Inv4033221=-1d0/144d0
       Inv4033222=-1d0/168d0
       Inv4033311=-1d0/160d0
       Inv4033321=-1d0/192d0
       Inv4033322=-1d0/224d0
       Inv4033331=-1d0/240d0
       Inv4033332=-1d0/280d0
       Inv4033333=-1d0/336d0
       Inv601=(-18*m0sq-24*(m1sq+m2sq+m3sq)+4*(p1sq+p4sq+s12)+5*(p2sq+p3
     &  sq+s23))/1440.d0
       Inv602=(-12*(m0sq+m1sq)-18*(m2sq+m3sq)+2*p1sq+4*p3sq+3*(p2sq+p4sq
     &  +s12+s23))/1440.d0
       Inv603=(-6*(m0sq+m1sq+m2sq)-12*m3sq+p1sq+p2sq+s12+2*(p3sq+p4sq+s2
     &  3))/1440.d0
       Inv6011=(42*m0sq+5*(14*(m1sq+m2sq+m3sq)-2*(p1sq+p4sq+s12)-3*(p2sq
     &  +p3sq+s23)))/5040.d0
       Inv6021=(56*m0sq+70*m1sq+105*(m2sq+m3sq)-10*p1sq-24*p3sq-15*(p4sq
     &  +s12)-18*(p2sq+s23))/10080.d0
       Inv6022=(21*(m0sq+m1sq)+42*(m2sq+m3sq)-3*p1sq-10*p3sq-6*(p2sq+p4s
     &  q+s12+s23))/5040.d0
       Inv6031=(28*m0sq+35*(m1sq+m2sq)+70*m3sq-6*p2sq-10*p4sq-5*(p1sq+s1
     &  2)-12*(p3sq+s23))/10080.d0
       Inv6032=(21*(m0sq+m1sq)+28*m2sq+56*m3sq-3*p1sq-10*p3sq-4*(p2sq+s1
     &  2)-8*(p4sq+s23))/10080.d0
       Inv6033=(7*(m0sq+m1sq+m2sq)+21*m3sq-p1sq-p2sq-s12-3*(p3sq+p4sq+s2
     &  3))/5040.d0
       Inv60111=(-16*m0sq-32*(m1sq+m2sq+m3sq)+4*(p1sq+p4sq+s12)+7*(p2sq+
     &  p3sq+s23))/2688.d0
       Inv60211=(-160*m0sq-240*m1sq-360*(m2sq+m3sq)+30*p1sq+84*p3sq+45*(
     &  p4sq+s12)+63*(p2sq+s23))/40320.d0
       Inv60221=(-60*m0sq-72*m1sq-144*(m2sq+m3sq)+9*p1sq+35*p3sq+18*(p4s
     &  q+s12)+21*(p2sq+s23))/20160.d0
       Inv60222=(-16*(m0sq+m1sq)-40*(m2sq+m3sq)+2*p1sq+10*p3sq+5*(p2sq+p
     &  4sq+s12+s23))/6720.d0
       Inv60311=(-80*m0sq+3*(-40*(m1sq+m2sq)-80*m3sq+7*p2sq+10*p4sq+5*(p
     &  1sq+s12)+14*(p3sq+s23)))/40320.d0
       Inv60321=(-60*m0sq-72*m1sq-96*m2sq-192*m3sq+9*p1sq+14*p2sq+35*p3s
     &  q+24*p4sq+12*s12+28*s23)/40320.d0
       Inv60322=(-24*(m0sq+m1sq)-40*m2sq-80*m3sq+3*p1sq+15*p3sq+5*(p2sq+
     &  s12)+10*(p4sq+s23))/20160.d0
       Inv60331=(-40*m0sq-48*(m1sq+m2sq)-144*m3sq+7*p2sq+18*p4sq+6*(p1sq
     &  +s12)+21*(p3sq+s23))/40320.d0
       Inv60332=(-32*(m0sq+m1sq)-40*m2sq-120*m3sq+4*p1sq+18*p3sq+5*(p2sq
     &  +s12)+15*(p4sq+s23))/40320.d0
       Inv60333=(-8*(m0sq+m1sq+m2sq)-32*m3sq+p1sq+p2sq+s12+4*(p3sq+p4sq+
     &  s23))/13440.d0
       Inv801=(-84*m0sq**2-140*(m1sq**2+m2sq**2+m2sq*m3sq+m3sq**2)+40*(m
     &  3sq*p4sq+m2sq*s12)+20*(m2sq*(p1sq+p4sq)+m3sq*(p1sq+s12))-5*(p1sq
     &  **2+p4sq**2+p4sq*s12+s12**2+p1sq*(p4sq+s12))+24*(m3sq*p2sq+m2sq*
     &  s23)-3*(p1sq*p3sq+p2sq*p4sq+s12*s23)-6*(p2sq*s12+p3sq*(p4sq+s12)
     &  +p4sq*s23+p1sq*(p2sq+s23))+48*(m2sq*(p2sq+p3sq)+m3sq*(p3sq+s23))
     &  -7*(p2sq**2+p3sq**2+p3sq*s23+s23**2+p2sq*(p3sq+s23))-4*(m1sq*(35
     &  *(m2sq+m3sq)-10*p1sq-6*p3sq-5*(p4sq+s12)-12*(p2sq+s23))+m0sq*(28
     &  *(m1sq+m2sq+m3sq)-8*(p1sq+p4sq+s12)-5*(p2sq+p3sq+s23))))/80640.d
     &  0
       Inv802=(-56*(m0sq**2+m1sq**2)-112*(m2sq**2+m2sq*m3sq+m3sq**2)+12*
     &  (m2sq+m3sq)*p1sq+40*(m2sq+m3sq)*p3sq-6*p3sq**2-3*p1sq*(p2sq+p4sq
     &  +s12+s23)-5*p3sq*(p2sq+p4sq+s12+s23)-2*(p1sq**2+p1sq*p3sq+p2sq*p
     &  4sq+s12*s23)+16*(m3sq*(p2sq+s12)+m2sq*(p4sq+s23))+32*(m2sq*(p2sq
     &  +s12)+m3sq*(p4sq+s23))-4*(p2sq**2+p4sq**2+p2sq*s12+p4sq*s12+s12*
     &  *2+p2sq*s23+p4sq*s23+s23**2+m1sq*(21*(m2sq+m3sq)-4*(p1sq+p3sq)-3
     &  *(p4sq+s12)-6*(p2sq+s23))+m0sq*(14*m1sq+21*(m2sq+m3sq)-4*(p1sq+p
     &  3sq)-6*(p4sq+s12)-3*(p2sq+s23))))/80640.d0
       Inv803=(-28*(m0sq**2+m1sq**2+m2sq**2)-56*m2sq*m3sq-84*m3sq**2-p1s
     &  q**2-p1sq*p2sq-p2sq**2+16*m2sq*p3sq-p1sq*p3sq-p2sq*p4sq-p1sq*s12
     &  -p2sq*s12-s12**2-s12*s23+24*m3sq*(p3sq+p4sq+s23)-2*((p3sq+p4sq)*
     &  s12+p2sq*(p3sq+s23)+p1sq*(p4sq+s23))-3*(p3sq**2+p4sq**2+p4sq*s23
     &  +s23**2+p3sq*(p4sq+s23))+8*(m3sq*(p1sq+p2sq+s12)+m2sq*(p2sq+p4sq
     &  +s12+s23))+4*(m2sq*p1sq-m1sq*(7*m2sq+14*m3sq-2*(p1sq+p2sq+p3sq+p
     &  4sq)-s12-4*s23)-m0sq*(7*(m1sq+m2sq)+14*m3sq-p2sq-4*p4sq-2*(p1sq+
     &  p3sq+s12+s23))))/80640.d0
 100  I4Z=1d0/(ZMax*4d0)
      I8Z=1d0/(ZMax*8d0)
      I12Z=1d0/(ZMax*12d0)
      I16Z=1d0/(ZMax*16d0)
      I20Z=1d0/(ZMax*20d0)
      I24Z=1d0/(Zmax*24d0)
      I28Z=1d0/(Zmax*28d0)
      I32Z=1d0/(Zmax*32d0)
      
       F(1)=+r10*ZZ(k,1,l,1)+r21*ZZ(k,2,l,1)+r32*ZZ(k,3,l,1)
       F(2)=+r10*ZZ(k,1,l,2)+r21*ZZ(k,2,l,2)+r32*ZZ(k,3,l,2)
       F(3)=+r10*ZZ(k,1,l,3)+r21*ZZ(k,2,l,3)+r32*ZZ(k,3,l,3)
       F(4)=2*r10*ZZ(k,1,l,1)+r21*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+
     &       r32*(ZZ(k,1,l,3)+ZZ(k,3,l,1))
       F(5)=r10**2*ZZ(k,1,l,1)+r10*r21*ZZ(k,1,l,2)+r10*r32*ZZ(k,1,l,3)
     &      +r10*r21*ZZ(k,2,l,1)+r21**2*ZZ(k,2,l,2)+r21*r32*ZZ(k,2,l,3)
     &      +r10*r32*ZZ(k,3,l,1)+r21*r32*ZZ(k,3,l,2)+r32**2*ZZ(k,3,l,3)+2*m0sq*Z(k,l)
       F(6)=r10*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+2*r21*ZZ(k,2,l,2)+
     &        r32*(ZZ(k,2,l,3)+ZZ(k,3,l,2))
       F(7)=r10*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+r21*(ZZ(k,2,l,3)+
     &      ZZ(k,3,l,2))+2*r32*ZZ(k,3,l,3)

       det4=-2*(p3sq*p1sqsq-(-p3sqsq+(p4sq+s12+s23)*p3sq
     &      +p2sq*(p3sq+p4sq-s12)+(s12-p4sq)*s23)*p1sq+p2sqsq*p4sq
     &      -p2sq*(-p4sq**2+(s12+s23)*p4sq+p3sq*(p4sq-s23)
     &      +s12*s23)+s12*(p3sq*(p4sq-s23)+s23*(-p4sq+s12+s23)))


       det4Abs= 2d0*(4*Abs(p1sq*p2sq*p3sq) + Abs(p3sq*(p1sq + p2sq - s12)**2) +
     1    Abs(p1sq*(p2sq + p3sq - s23)**2) +
     2    Abs((p1sq + p2sq - s12)*(p2sq + p3sq - s23)*(p2sq + p4sq - s12 - s23)) +
     3    Abs(p2sq*(p2sq + p4sq - s12 - s23)**2))
       IF(PRINTMY) THEN
       print*, 'IZax4', I4Z
       print*,'Z(1,1)',Z(1,1)
       print*,'Z(1,2)',Z(1,2)
       print*,'Z(1,3)',Z(1,3)
       print*,'Z(2,2)',Z(2,2)
       print*,'Z(2,3)',Z(2,3)
       print*,'Z(3,3)',Z(3,3)
       Print*, 'jj',jj
       Print*, 'k',k
       Print*, 'l',l
       print*, "T1",TX1
       print*, "T1ABS",TX1/( Abs(p2sqsq*r10) + 
     -  2*Abs(p2sq*p3sq*r10) + 
     -  Abs(p3sqsq*r10) + Abs(p2sqsq*r21) + 
     -  2*Abs(p1sq*p3sq*r21) + 
     -  Abs(p2sq*p3sq*r21) + 
     -  Abs(p2sq*p4sq*r21) + 
     -  Abs(p3sq*p4sq*r21) + 
     -  Abs(p1sq*p2sq*r32) + Abs(p2sqsq*r32) + 
     -  Abs(p1sq*p3sq*r32) + 
     -  Abs(p2sq*p3sq*r32) + 
     -  2*Abs(p2sq*p4sq*r32) + 
     -  Abs(p2sq*r21*s12) + 
     -  Abs(p3sq*r21*s12) + 
     -  Abs(p2sq*r32*s12) + 
     -  Abs(p3sq*r32*s12) + 
     -  2*Abs(p2sq*r10*s23) + 
     -  2*Abs(p3sq*r10*s23) + 
     -  2*Abs(p2sq*r21*s23) + 
     -  Abs(p3sq*r21*s23) + 
     -  Abs(p4sq*r21*s23) + 
     -  Abs(p1sq*r32*s23) + 
     -  Abs(p2sq*r32*s23) + Abs(r21*s12*s23) + 
     -  Abs(r32*s12*s23) + Abs(r10*s23**2) + 
     -  Abs(r21*s23**2))


       print*, "T2",TX2
       print*, "T2ABS", TX2/(Abs(p2sqsq*r10) + 
     -  2*Abs(p1sq*p3sq*r10) + 
     -  Abs(p2sq*p3sq*r10) + 
     -  Abs(p2sq*p4sq*r10) + 
     -  Abs(p3sq*p4sq*r10) + Abs(p2sqsq*r21) + 
     -  4*Abs(p1sq*p3sq*r21) + 
     -  2*Abs(p2sq*p4sq*r21) + 
     -  Abs(p4sqsq*r21) + Abs(p1sq*p2sq*r32) + 
     -  Abs(p2sqsq*r32) + 
     -  2*Abs(p1sq*p3sq*r32) + 
     -  Abs(p1sq*p4sq*r32) + 
     -  Abs(p2sq*p4sq*r32) + 
     -  Abs(p2sq*r10*s12) + 
     -  Abs(p3sq*r10*s12) + 
     -  2*Abs(p2sq*r21*s12) + 
     -  2*Abs(p4sq*r21*s12) + 
     -  Abs(p1sq*r32*s12) + 
     -  2*Abs(p2sq*r32*s12) + 
     -  Abs(p4sq*r32*s12) + Abs(r21*s12sq) + 
     -  Abs(r32*s12sq) + 2*Abs(p2sq*r10*s23) + 
     -  Abs(p3sq*r10*s23) + 
     -  Abs(p4sq*r10*s23) + 
     -  2*Abs(p2sq*r21*s23) + 
     -  2*Abs(p4sq*r21*s23) + 
     -  Abs(p1sq*r32*s23) + 
     -  Abs(p2sq*r32*s23) + Abs(r10*s12*s23) + 
     -  2*Abs(r21*s12*s23) + 
     -  Abs(r32*s12*s23) + Abs(r10*s23sq) + 
     -  Abs(r21*s23sq))


       print*, "T3",TX3
       print*, "T3ABS",      TX3/( Abs(p1sq*p2sq*r10) + Abs(p2sqsq*r10) + Abs(p1sq*p3sq*r10) + Abs(p2sq*p3sq*r10) + 
     -  2*Abs(p2sq*p4sq*r10) + Abs(p1sq*p2sq*r21) + Abs(p2sqsq*r21) + 2*Abs(p1sq*p3sq*r21) + 
     -  Abs(p1sq*p4sq*r21) + Abs(p2sq*p4sq*r21) + Abs(p1sqsq*r32) + 2*Abs(p1sq*p2sq*r32) + 
     -  Abs(p2sqsq*r32) + Abs(p2sq*r10*s12) + Abs(p3sq*r10*s12) + Abs(p1sq*r21*s12) + 
     -  2*Abs(p2sq*r21*s12) + Abs(p4sq*r21*s12) + 2*Abs(p1sq*r32*s12) + 2*Abs(p2sq*r32*s12) + 
     -  Abs(r21*s12**2) + Abs(r32*s12**2) + Abs(p1sq*r10*s23) + Abs(p2sq*r10*s23) + 
     -  Abs(p1sq*r21*s23) + Abs(p2sq*r21*s23) + Abs(r10*s12*s23) + Abs(r21*s12*s23))
       Print*, 'det4*IX',det4*IX
       Print*, 'det4*IZmax4',det4*I4Z
          ENDIF
c                Iteration0
       S31(1)=C0134-C0234
       S31(2)=C0124-C0134
       S31(3)=C0123-C0124
       auxD40=-(S31(1)*Z(jj,1))-S31(2)*Z(jj,2)-S31(3)*Z(jj,3)
       tempD40=auxD40*IX
       if(printmy) then
       print*, "tempD40 0",tempD40,order,k,l,jj
       endif
       if(order.eq.0) goto 500
c                Iteration1
c                Step1
       S321(1)=C0234+Cij134(1,1)
       S321(2)=Cij134(1,1)-Cij234(1,1)
       S321(3)=Cij134(2,1)-Cij234(2,1)
       S322(1)=Cij124(1,1)-Cij134(1,1)
       S322(2)=Cij124(2,1)-Cij134(1,1)
       S322(3)=Cij124(2,1)-Cij134(2,1)
       S323(1)=Cij123(1,1)-Cij124(1,1)
       S323(2)=Cij123(2,1)-Cij124(2,1)
       S323(3)=-Cij124(2,1)
       S300=2*C0234
       auxD400=-(F(1)*S31(1))-F(2)*S31(2)-F(3)*S31(3)+S321(k)*Z(1,l)+S32
     &  2(k)*Z(2,l)+S323(k)*Z(3,l)+(S300-S321(1)-S322(2)-S323(3))*Z(k,l)
       tempD400=I4Z*(auxD400+tempD40*F(5))
       aux1(1)=-(S321(1)*Z(jj,1))-S322(1)*Z(jj,2)-S323(1)*Z(jj,3)
       aux1(2)=-(S321(2)*Z(jj,1))-S322(2)*Z(jj,2)-S323(2)*Z(jj,3)
       aux1(3)=-(S321(3)*Z(jj,1))-S322(3)*Z(jj,2)-S323(3)*Z(jj,3)
       temp1(1)=IX*(aux1(1)+2*tempD400*Z(jj,1))
       temp1(2)=IX*(aux1(2)+2*tempD400*Z(jj,2))
       temp1(3)=IX*(aux1(3)+2*tempD400*Z(jj,3))
c                Step2
       tempD40=IX*(auxD40+det4*temp1(jj))
       if(printmy) then
       print*, "tempD40 1",tempD40
       endif
       if(order.eq.1) goto 500
c                Iteration2
c                Step1
       S3311(1)=-C0234+Cij134(1,2)
       S3311(2)=Cij134(1,2)+Cij234(1,1)
       S3311(3)=Cij134(3,2)+Cij234(2,1)
       S3312(1)=Cij134(1,2)+Cij234(1,1)
       S3312(2)=Cij134(1,2)-Cij234(1,2)
       S3312(3)=Cij134(3,2)-Cij234(3,2)
       S3313(1)=Cij134(3,2)+Cij234(2,1)
       S3313(2)=Cij134(3,2)-Cij234(3,2)
       S3313(3)=Cij134(2,2)-Cij234(2,2)
       S3321(1)=Cij124(1,2)-Cij134(1,2)
       S3321(2)=Cij124(3,2)-Cij134(1,2)
       S3321(3)=Cij124(3,2)-Cij134(3,2)
       S3322(1)=Cij124(3,2)-Cij134(1,2)
       S3322(2)=Cij124(2,2)-Cij134(1,2)
       S3322(3)=Cij124(2,2)-Cij134(3,2)
       S3323(1)=Cij124(3,2)-Cij134(3,2)
       S3323(2)=Cij124(2,2)-Cij134(3,2)
       S3323(3)=Cij124(2,2)-Cij134(2,2)
       S3331(1)=Cij123(1,2)-Cij124(1,2)
       S3331(2)=Cij123(3,2)-Cij124(3,2)
       S3331(3)=-Cij124(3,2)
       S3332(1)=Cij123(3,2)-Cij124(3,2)
       S3332(2)=Cij123(2,2)-Cij124(2,2)
       S3332(3)=-Cij124(2,2)
       S3333(1)=-Cij124(3,2)
       S3333(2)=-Cij124(2,2)
       S3333(3)=-Cij124(2,2)
       S3h001(1)=Cij134(4,2)-Cij234(4,2)
       S3h001(2)=Cij124(4,2)-Cij134(4,2)
       S3h001(3)=Cij123(4,2)-Cij124(4,2)
       S3001(1)=-2*C0234
       S3001(2)=2*Cij234(1,1)
       S3001(3)=2*Cij234(2,1)
       aux001(1)=-(F(1)*S321(1))-F(2)*S322(1)-F(3)*S323(1)+S3311(k)*Z(1,
     &  l)+S3321(k)*Z(2,l)+S3331(k)*Z(3,l)+S3001(1)*Z(k,l)-S3311(1)*Z(k,
     &  l)-S3322(1)*Z(k,l)-S3333(1)*Z(k,l)-2*S3h001(1)*ZZ(k,1,l,1)-2*S3h
     &  001(2)*ZZ(k,1,l,2)-2*S3h001(3)*ZZ(k,1,l,3)
       aux001(2)=-(F(1)*S321(2))-F(2)*S322(2)-F(3)*S323(2)+S3312(k)*Z(1,
     &  l)+S3322(k)*Z(2,l)+S3332(k)*Z(3,l)+S3001(2)*Z(k,l)-S3312(1)*Z(k,
     &  l)-S3322(2)*Z(k,l)-S3333(2)*Z(k,l)-2*S3h001(1)*ZZ(k,2,l,1)-2*S3h
     &  001(2)*ZZ(k,2,l,2)-2*S3h001(3)*ZZ(k,2,l,3)
       aux001(3)=-(F(1)*S321(3))-F(2)*S322(3)-F(3)*S323(3)+S3313(k)*Z(1,
     &  l)+S3323(k)*Z(2,l)+S3333(k)*Z(3,l)+S3001(3)*Z(k,l)-S3313(1)*Z(k,
     &  l)-S3323(2)*Z(k,l)-S3333(3)*Z(k,l)-2*S3h001(1)*ZZ(k,3,l,1)-2*S3h
     &  001(2)*ZZ(k,3,l,2)-2*S3h001(3)*ZZ(k,3,l,3)
       temp001(1)=I8Z*(aux001(1)+2*tempD400*F(4)+F(5)*temp1(1))
       temp001(2)=I8Z*(aux001(2)+2*tempD400*F(6)+F(5)*temp1(2))
       temp001(3)=I8Z*(aux001(3)+2*tempD400*F(7)+F(5)*temp1(3))
       aux2(1,1)=-(S3311(1)*Z(jj,1))-S3321(1)*Z(jj,2)-S3331(1)*Z(jj,3)
       aux2(2,1)=-(S3312(1)*Z(jj,1))-S3322(1)*Z(jj,2)-S3332(1)*Z(jj,3)
       aux2(2,2)=-(S3312(2)*Z(jj,1))-S3322(2)*Z(jj,2)-S3332(2)*Z(jj,3)
       aux2(3,1)=-(S3313(1)*Z(jj,1))-S3323(1)*Z(jj,2)-S3333(1)*Z(jj,3)
       aux2(3,2)=-(S3313(2)*Z(jj,1))-S3323(2)*Z(jj,2)-S3333(2)*Z(jj,3)
       aux2(3,3)=-(S3313(3)*Z(jj,1))-S3323(3)*Z(jj,2)-S3333(3)*Z(jj,3)
       temp2(1,1)=IX*(aux2(1,1)+4*temp001(1)*Z(jj,1))
       temp2(2,1)=IX*(aux2(2,1)+2*temp001(2)*Z(jj,1)+2*temp001(1)*Z(jj,2
     &  ))
       temp2(2,2)=IX*(aux2(2,2)+4*temp001(2)*Z(jj,2))
       temp2(3,1)=IX*(aux2(3,1)+2*temp001(3)*Z(jj,1)+2*temp001(1)*Z(jj,3
     &  ))
       temp2(3,2)=IX*(aux2(3,2)+2*temp001(3)*Z(jj,2)+2*temp001(2)*Z(jj,3
     &  ))
       temp2(3,3)=IX*(aux2(3,3)+4*temp001(3)*Z(jj,3))
       temp2(1,2)=temp2(2,1)
       temp2(1,3)=temp2(3,1)
       temp2(2,3)=temp2(3,2)
c                Step2
       tempD400=I4Z*(auxD400+tempD40*F(5)-det4*temp2(k,l))
       temp1(1)=IX*(aux1(1)+det4*temp2(1,jj)+2*tempD400*Z(jj,1))
       temp1(2)=IX*(aux1(2)+det4*temp2(2,jj)+2*tempD400*Z(jj,2))
       temp1(3)=IX*(aux1(3)+det4*temp2(3,jj)+2*tempD400*Z(jj,3))
c                Step3
       tempD40=IX*(auxD40+det4*temp1(jj))

c FC %       do jjtemp=1,3
c FC %       tempjj1=abs(det4*temp2(1,jjtemp)+2*tempD400*Z(jjtemp,1))/abs(aux1(1))
c FC %       tempjj2=abs(det4*temp2(2,jjtemp)+2*tempD400*Z(jjtemp,2))/abs(aux1(2))
c FC %       tempjj3=abs(det4*temp2(3,jjtemp)+2*tempD400*Z(jjtemp,3))/abs(aux1(3))
c FC %       Print*,'tempjj1',tempjj1
c FC %       Print*,'tempjj2',tempjj2
c FC %       Print*,'tempjj3',tempjj3
c FC %       enddo


c FC1 %       If(tempjj1.lt.tempjj2)then
c FC1 %          if(tempjj1.lt.tempjj3) then
c FC1 %             jjtemp=1
c FC1 %             IXtemp=1/TX1
c FC1 %             tempjj=tempjj1/abs(auxD40)*abs(TX1)
c FC1 %          else
c FC1 %             jjtemp=3
c FC1 %             IXtemp=1/TX3
c FC1 %             tempjj=tempjj3/abs(auxD40)*abs(TX3)
c FC1 %          endif
c FC1 %       else
c FC1 %          if(tempjj2.lt.tempjj3) then
c FC1 %             jjtemp=2
c FC1 %             IXtemp=1/TX2
c FC1 %             tempjj=tempjj2/abs(auxD40)*abs(TX2)
c FC1 %          else
c FC1 %             jjtemp=3
c FC1 %             IXtemp=1/TX3
c FC1 %             tempjj=tempjj3/abs(auxD40)*abs(TX3)
c FC1 %          endif
c FC1 %       endif
c FC1 %
c FC1 %
c FC1 %c FC %       print*, 'tempjj',tempjj
c FC1 %c FC %       print*, 'jjtemp',jjtemp
c FC1 % 
c FC1 %
       tempjj1=abs(det4*temp1(1)/TX1)
       tempjj2=abs(det4*temp1(2)/TX2)
       tempjj3=abs(det4*temp1(3)/TX3)
       if(printmy) then
       Print*,'tempjj3',tempjj3
       endif
       tempjj3=tempjj3*0.62d0

       If(tempjj1.lt.tempjj2)then
          if(tempjj1.lt.tempjj3) then
             jjtemp=1
             IXtemp=1/TX1
             tempjj=tempjj1/abs(auxD40)*abs(TX1)
          else
             jjtemp=3
             IXtemp=1/TX3
             tempjj=tempjj3/abs(auxD40)*abs(TX3)
          endif
       else
          if(tempjj2.lt.tempjj3) then
             jjtemp=2
             IXtemp=1/TX2
             tempjj=tempjj2/abs(auxD40)*abs(TX2)
          else
             jjtemp=3
             IXtemp=1/TX3
             tempjj=tempjj3/abs(auxD40)*abs(TX3)
          endif
       endif

       if(printmy) then
       Print*,'tempjj1',tempjj1
       Print*,'tempjj2',tempjj2
       Print*,'tempjj3',tempjj3
       print*, 'tempjj',tempjj
       print*, 'jjtemp',jjtemp
       print*, "F(5)",F(5)
       endif


       tempkl1=abs(det4*temp2(1,1)/(Z(1,1)*4d0))
c       tempkl1=tempkl1*1.5d0
       tempkl2=abs(det4*temp2(2,2)/(Z(2,2)*4d0))
       tempkl3=abs(det4*temp2(3,3)/(Z(3,3)*4d0))
       tempkl3=tempkl3*0.62d0
       tempkl4=abs(det4*temp2(2,1)/(Z(2,1)*4d0))
       tempkl5=abs(det4*temp2(3,1)/(Z(3,1)*4d0))
       tempkl6=abs(det4*temp2(3,2)/(Z(3,2)*4d0))

        if(printmy) then
       print*, "tempkl1",tempkl1,temp2(1,1)
       print*, "tempkl2",tempkl2,temp2(2,2)
       print*, "tempkl3",tempkl3,temp2(3,3)
       print*, "tempkl4",tempkl4,temp2(2,1)
       print*, "tempkl5",tempkl5,temp2(3,1)
       print*, "tempkl6",tempkl6,temp2(3,2)
           endif


       if(tempkl1.lt.tempkl2) then
          if(tempkl1*2.0d0.lt.tempkl3) then
             if(tempkl1.lt.tempkl4) then
                if(tempkl1.lt.tempkl5) then
                   if(tempkl1.lt.tempkl6)then
                      ktemp=1
                      ltemp=1
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(1,1)
                      tempkl=tempkl1/abs(auxD400+tempD40*F1)*abs(Z(1,1)*4d0)
                    else
                       ktemp=3
                       ltemp=2
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                       Zmaxtemp=Z(3,2)
                       tempkl=tempkl6/abs(auxD400+tempD40*F1)*abs(Z(3,2)*4d0)
                    endif
                else
                   if(tempkl5.lt.tempkl6)then
                      ktemp=3
                      ltemp=1
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,1)
                      tempkl=tempkl5/abs(auxD400+tempD40*F1)*abs(Z(3,1)*4d0)
                    else
                      ktemp=3
                      ltemp=2
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,2)
                      tempkl=tempkl6/abs(auxD400+tempD40*F1)*abs(Z(3,2)*4d0)
                    endif
                endif    
             else  
                 if(tempkl4.lt.tempkl5) then
                   if(tempkl4.lt.tempkl6)then
                      ktemp=2
                      ltemp=1
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(2,1)
                      tempkl=tempkl4/abs(auxD400+tempD40*F1)*abs(Z(2,1)*4d0)
                    else
                      ktemp=3
                      ltemp=2
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,2)
                      tempkl=tempkl6/abs(auxD400+tempD40*F1)*abs(Z(3,2)*4d0)
                    endif
                else
                   if(tempkl5.lt.tempkl6)then
                      ktemp=3
                      ltemp=1
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,1)
                      tempkl=tempkl5/abs(auxD400+tempD40*F1)*abs(Z(3,1)*4d0)
                    else
                      ktemp=3
                      ltemp=2
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,2)
                      tempkl=tempkl6/abs(auxD400+tempD40*F1)*abs(Z(3,2)*4d0)
                    endif
                endif   
             endif
          else   
             if(tempkl3.lt.tempkl4) then
                if(tempkl3.lt.tempkl5) then
                   if(tempkl3.lt.tempkl6)then
                      ktemp=3
                      ltemp=3
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,3)
                      tempkl=tempkl3/abs(auxD400+tempD40*F1)*abs(Z(3,3)*4d0)
                    else
                      ktemp=3
                      ltemp=2
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,2)
                      tempkl=tempkl6/abs(auxD400+tempD40*F1)*abs(Z(3,2)*4d0)
                    endif
                else
                   if(tempkl5.lt.tempkl6)then
                      ktemp=3
                      ltemp=1
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,1)
                      tempkl=tempkl5/abs(auxD400+tempD40*F1)*abs(Z(3,1)*4d0)
                    else
                      ktemp=3
                      ltemp=2
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,2)
                      tempkl=tempkl6/abs(auxD400+tempD40*F1)*abs(Z(3,2)*4d0)
                    endif
                endif    
             else  
                 if(tempkl4.lt.tempkl5) then
                   if(tempkl4.lt.tempkl6)then
                      ktemp=2
                      ltemp=1
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(2,1)
                      tempkl=tempkl4/abs(auxD400+tempD40*F1)*abs(Z(2,1)*4d0)
                    else
                      ktemp=3
                      ltemp=2
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,2)
                      tempkl=tempkl6/abs(auxD400+tempD40*F1)*abs(Z(3,2)*4d0)
                    endif
                else
                   if(tempkl5.lt.tempkl6)then
                      ktemp=3
                      ltemp=1
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,1)
                      tempkl=tempkl5/abs(auxD400+tempD40*F1)*abs(Z(3,1)*4d0)
                    else
                      ktemp=3
                      ltemp=2
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,2)
                      tempkl=tempkl6/abs(auxD400+tempD40*F1)*abs(Z(3,2)*4d0)
                    endif
                endif   
             endif
          endif
       else   
         if(tempkl2.lt.tempkl3) then
             if(tempkl2.lt.tempkl4) then
                if(tempkl2.lt.tempkl5) then
                   if(tempkl2.lt.tempkl6)then
                      ktemp=2
                      ltemp=2
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(2,2)
                      tempkl=tempkl2/abs(auxD400+tempD40*F1)*abs(Z(2,2)*4d0)
                    else
                      ktemp=3
                      ltemp=2
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,2)
                      tempkl=tempkl6/abs(auxD400+tempD40*F1)*abs(Z(3,2)*4d0)
                    endif
                else
                   if(tempkl5.lt.tempkl6)then
                      ktemp=3
                      ltemp=1
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,1)
                      tempkl=tempkl5/abs(auxD400+tempD40*F1)*abs(Z(3,1)*4d0)
                    else
                      ktemp=3
                      ltemp=2
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,2)
                      tempkl=tempkl5/abs(auxD400+tempD40*F1)*abs(Z(3,2)*4d0)
                    endif
                endif    
             else  
                 if(tempkl4.lt.tempkl5) then
                   if(tempkl4.lt.tempkl6)then
                      ktemp=2
                      ltemp=1
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(2,1)
                      tempkl=tempkl4/abs(auxD400+tempD40*F1)*abs(Z(2,1)*4d0)
                    else
                      ktemp=3
                      ltemp=2
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,2)
                      tempkl=tempkl6/abs(auxD400+tempD40*F1)*abs(Z(3,2)*4d0)
                    endif
                else
                   if(tempkl5.lt.tempkl6)then
                      ktemp=3
                      ltemp=1
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,1)
                      tempkl=tempkl5/abs(auxD400+tempD40*F1)*abs(Z(3,1)*4d0)
                    else
                      ktemp=3
                      ltemp=2
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,2)
                      tempkl=tempkl6/abs(auxD400+tempD40*F1)*abs(Z(3,2)*4d0)
                    endif
                endif   
             endif
          else   
             if(tempkl3.lt.tempkl4) then
                if(tempkl3.lt.tempkl5) then
                   if(tempkl3.lt.tempkl6)then
                      ktemp=3
                      ltemp=3
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,3)
                      tempkl=tempkl3/abs(auxD400+tempD40*F1)*abs(Z(3,3)*4d0)
                    else
                      ktemp=3
                      ltemp=2
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,2)
                      tempkl=tempkl6/abs(auxD400+tempD40*F1)*abs(Z(3,2)*4d0)
                    endif
                else
                   if(tempkl5.lt.tempkl6)then
                      ktemp=3
                      ltemp=1
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,1)
                      tempkl=tempkl5/abs(auxD400+tempD40*F1)*abs(Z(3,1)*4d0)
                    else
                      ktemp=3
                      ltemp=2
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,2)
                      tempkl=tempkl6/abs(auxD400+tempD40*F1)*abs(Z(3,2)*4d0)
                    endif
                endif    
             else  
                 if(tempkl4.lt.tempkl5) then
                   if(tempkl4.lt.tempkl6)then
                      ktemp=2
                      ltemp=1
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(2,1)
                      tempkl=tempkl4/abs(auxD400+tempD40*F1)*abs(Z(2,1)*4d0)
                    else
                      ktemp=3
                      ltemp=2
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,2)
                      tempkl=tempkl6/abs(auxD400+tempD40*F1)*abs(Z(3,2)*4d0)
                    endif
                else
                   if(tempkl5.lt.tempkl6)then
                      ktemp=3
                      ltemp=1
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,1)
                      tempkl=tempkl5/abs(auxD400+tempD40*F1)*abs(Z(3,1)*4d0)
                    else
                      ktemp=3
                      ltemp=2
             F1=r10**2*ZZ(ktemp,1,ltemp,1)+r10*r21*ZZ(ktemp,1,ltemp,2)+r10*r32*ZZ(ktemp,1,ltemp,3)
     &      +r10*r21*ZZ(ktemp,2,ltemp,1)+r21**2*ZZ(ktemp,2,ltemp,2)+r21*r32*ZZ(ktemp,2,ltemp,3)
     &      +r10*r32*ZZ(ktemp,3,ltemp,1)+r21*r32*ZZ(ktemp,3,ltemp,2)+r32**2*ZZ(ktemp,3,ltemp,3)+2*m0sq*Z(ktemp,ltemp)
                      Zmaxtemp=Z(3,2)
                      tempkl=tempkl6/abs(auxD400+tempD40*F1)*abs(Z(3,2)*4d0)
                    endif
                endif   
             endif
          endif
      endif  
    

c FC %      print*, 'ktemp',ktemp
c FC %       print*, 'ltemp',ltemp
cFC       Print*,'tempjj=',tempjj
cFC       Print*,'tempkl=',tempkl
c FC %

      if(cont.eq.4)  goto 200
      if(cont.lt.4) then
         if(ktemp.eq.k) then
             if (ltemp.eq.l) then
                  if (jjtemp.eq.jj) then
                    goto 200
                  else
                   jj=jjtemp
                   IX=IXtemp
                   k=ktemp
                   l=ltemp
                   Zmax=Zmaxtemp
                   cont=cont+1
                   
                   goto 100
                 endif
             else
                   jj=jjtemp
                   IX=IXtemp
                   k=ktemp
                   l=ltemp
                   Zmax=Zmaxtemp
                   cont=cont+1
                   goto 100 
            endif
         else
                    jj=jjtemp
                   IX=IXtemp
                   k=ktemp
                   l=ltemp
                   Zmax=Zmaxtemp
                   cont=cont+1
                   goto 100
        endif
        else
                   jj=jjinit
                   IX=IXinit
                   k=kinit
                   l=linit
                   Zmax=Zmaxinit
                   goto 100
                   cont=4
        endif


 200  If((abs(tempkl).ge.3.1d-2).and.(abs(tempjj).ge.1.5d-2)) then
       !return
       endif
!       IF(abs(tempjj).ge.1.5d-2) then
!       return
!       endif    
       IF(abs(IX*det4).ge.1.0d0) then!.and.abs(det4/det4Abs).gt.1.5d-2) then 
        if(printmy) then
        print*, "CAUTION SMALL CAYLEY,IX*det4 in D", abs(IX*det4),m0
        endif
       !return
       endif
c FC      Print*,"IX*det4", IX*det4
c FC       print*, 'k',k
c FC       print*, 'l',l
c FC       print*, 'jj',jj
c FC       print*, 'cont',cont
        if(printmy) then
       print*, "tempD40 2",tempD40,k,l,jj,IX,1d0/TX3
        endif
       if(order.eq.2) goto 500
c                Iteration3
c                Step1
       S3h0021(1)=Cij134(5,3)+Cij234(4,2)
       S3h0021(2)=Cij134(5,3)-Cij234(5,3)
       S3h0021(3)=Cij134(6,3)-Cij234(6,3)
       S3h0022(1)=Cij124(5,3)-Cij134(5,3)
       S3h0022(2)=Cij124(6,3)-Cij134(5,3)
       S3h0022(3)=Cij124(6,3)-Cij134(6,3)
       S3h0023(1)=Cij123(5,3)-Cij124(5,3)
       S3h0023(2)=Cij123(6,3)-Cij124(6,3)
       S3h0023(3)=-Cij124(6,3)
       S30000=2*Cij234(4,2)
       auxD40000=-(F(1)*S3h001(1))-F(2)*S3h001(2)-F(3)*S3h001(3)+S3h0021
     &  (k)*Z(1,l)+S3h0022(k)*Z(2,l)+S3h0023(k)*Z(3,l)+(Inv40+S30000-S3h
     &  0021(1)-S3h0022(2)-S3h0023(3))*Z(k,l)
       tempD40000=I8Z*(auxD40000+tempD400*F(5))
c FC       print*, 'tempD40000',tempD40000
c FC       print*, 'auxD40000',auxD40000
c FC       print*, 'tempD400',tempD400
c FC       print*, 'F(5)',F(5),F(4),F(3)
c FC       print*, ' S3h0021(1)', S3h0021(1), S3h0021(2), S3h0021(3)
c FC       print*, ' S3h0022(1)', S3h0022(1), S3h0022(2), S3h0022(3)
c FC       print*, ' S3h0023(1)', S3h0023(1), S3h0023(2), S3h0023(3)
c FC       print*, 'inv40', inv40
c FC       print*, ' S3h001(1)', S3h001(1), S3h001(2), S3h001(3)
c FC       print*, 'S30000',S30000

       S34111(1)=C0234+Cij134(1,3)
       S34111(2)=Cij134(1,3)-Cij234(1,1)
       S34111(3)=Cij134(3,3)-Cij234(2,1)
       S34121(1)=Cij134(1,3)-Cij234(1,1)
       S34121(2)=Cij134(1,3)+Cij234(1,2)
       S34121(3)=Cij134(3,3)+Cij234(3,2)
       S34122(1)=Cij134(1,3)+Cij234(1,2)
       S34122(2)=Cij134(1,3)-Cij234(1,3)
       S34122(3)=Cij134(3,3)-Cij234(3,3)
       S34131(1)=Cij134(3,3)-Cij234(2,1)
       S34131(2)=Cij134(3,3)+Cij234(3,2)
       S34131(3)=Cij134(4,3)+Cij234(2,2)
       S34132(1)=Cij134(3,3)+Cij234(3,2)
       S34132(2)=Cij134(3,3)-Cij234(3,3)
       S34132(3)=Cij134(4,3)-Cij234(4,3)
       S34133(1)=Cij134(4,3)+Cij234(2,2)
       S34133(2)=Cij134(4,3)-Cij234(4,3)
       S34133(3)=Cij134(2,3)-Cij234(2,3)
       S34211(1)=Cij124(1,3)-Cij134(1,3)
       S34211(2)=Cij124(3,3)-Cij134(1,3)
       S34211(3)=Cij124(3,3)-Cij134(3,3)
       S34221(1)=Cij124(3,3)-Cij134(1,3)
       S34221(2)=Cij124(4,3)-Cij134(1,3)
       S34221(3)=Cij124(4,3)-Cij134(3,3)
       S34222(1)=Cij124(4,3)-Cij134(1,3)
       S34222(2)=Cij124(2,3)-Cij134(1,3)
       S34222(3)=Cij124(2,3)-Cij134(3,3)
       S34231(1)=Cij124(3,3)-Cij134(3,3)
       S34231(2)=Cij124(4,3)-Cij134(3,3)
       S34231(3)=Cij124(4,3)-Cij134(4,3)
       S34232(1)=Cij124(4,3)-Cij134(3,3)
       S34232(2)=Cij124(2,3)-Cij134(3,3)
       S34232(3)=Cij124(2,3)-Cij134(4,3)
       S34233(1)=Cij124(4,3)-Cij134(4,3)
       S34233(2)=Cij124(2,3)-Cij134(4,3)
       S34233(3)=Cij124(2,3)-Cij134(2,3)
       S34311(1)=Cij123(1,3)-Cij124(1,3)
       S34311(2)=Cij123(3,3)-Cij124(3,3)
       S34311(3)=-Cij124(3,3)
       S34321(1)=Cij123(3,3)-Cij124(3,3)
       S34321(2)=Cij123(4,3)-Cij124(4,3)
       S34321(3)=-Cij124(4,3)
       S34322(1)=Cij123(4,3)-Cij124(4,3)
       S34322(2)=Cij123(2,3)-Cij124(2,3)
       S34322(3)=-Cij124(2,3)
       S34331(1)=-Cij124(3,3)
       S34331(2)=-Cij124(4,3)
       S34331(3)=-Cij124(4,3)
       S34332(1)=-Cij124(4,3)
       S34332(2)=-Cij124(2,3)
       S34332(3)=-Cij124(2,3)
       S34333(1)=-Cij124(4,3)
       S34333(2)=-Cij124(2,3)
       S34333(3)=-Cij124(2,3)
       S30021(1)=2*C0234
       S30021(2)=-2*Cij234(1,1)
       S30021(3)=-2*Cij234(2,1)
       S30022(1)=-2*Cij234(1,1)
       S30022(2)=2*Cij234(1,2)
       S30022(3)=2*Cij234(3,2)
       S30023(1)=-2*Cij234(2,1)
       S30023(2)=2*Cij234(3,2)
       S30023(3)=2*Cij234(2,2)
       aux002(1,1)=-(F(1)*S3311(1))-F(2)*S3321(1)-F(3)*S3331(1)+S34111(k
     &  )*Z(1,l)+S34211(k)*Z(2,l)+S34311(k)*Z(3,l)+(S30021(1)-S34111(1)-
     &  S34221(1)-S34331(1))*Z(k,l)-4*S3h0021(1)*ZZ(k,1,l,1)-4*S3h0022(1
     &  )*ZZ(k,1,l,2)-4*S3h0023(1)*ZZ(k,1,l,3)
       aux002(2,1)=-(F(1)*S3312(1))-F(2)*S3322(1)-F(3)*S3332(1)+S34121(k
     &  )*Z(1,l)+S34221(k)*Z(2,l)+S34321(k)*Z(3,l)+(S30022(1)-S34121(1)-
     &  S34222(1)-S34332(1))*Z(k,l)-2*S3h0021(2)*ZZ(k,1,l,1)-2*S3h0022(2
     &  )*ZZ(k,1,l,2)-2*S3h0023(2)*ZZ(k,1,l,3)-2*S3h0021(1)*ZZ(k,2,l,1)-
     &  2*S3h0022(1)*ZZ(k,2,l,2)-2*S3h0023(1)*ZZ(k,2,l,3)
       aux002(2,2)=-(F(1)*S3312(2))-F(2)*S3322(2)-F(3)*S3332(2)+S34122(k
     &  )*Z(1,l)+S34222(k)*Z(2,l)+S34322(k)*Z(3,l)+(S30022(2)-S34122(1)-
     &  S34222(2)-S34332(2))*Z(k,l)-4*S3h0021(2)*ZZ(k,2,l,1)-4*S3h0022(2
     &  )*ZZ(k,2,l,2)-4*S3h0023(2)*ZZ(k,2,l,3)
       aux002(3,1)=-(F(1)*S3313(1))-F(2)*S3323(1)-F(3)*S3333(1)+S34131(k
     &  )*Z(1,l)+S34231(k)*Z(2,l)+S34331(k)*Z(3,l)+(S30023(1)-S34131(1)-
     &  S34232(1)-S34333(1))*Z(k,l)-2*S3h0021(3)*ZZ(k,1,l,1)-2*S3h0022(3
     &  )*ZZ(k,1,l,2)-2*S3h0023(3)*ZZ(k,1,l,3)-2*S3h0021(1)*ZZ(k,3,l,1)-
     &  2*S3h0022(1)*ZZ(k,3,l,2)-2*S3h0023(1)*ZZ(k,3,l,3)
       aux002(3,2)=-(F(1)*S3313(2))-F(2)*S3323(2)-F(3)*S3333(2)+S34132(k
     &  )*Z(1,l)+S34232(k)*Z(2,l)+S34332(k)*Z(3,l)+(S30023(2)-S34132(1)-
     &  S34232(2)-S34333(2))*Z(k,l)-2*S3h0021(3)*ZZ(k,2,l,1)-2*S3h0022(3
     &  )*ZZ(k,2,l,2)-2*S3h0023(3)*ZZ(k,2,l,3)-2*S3h0021(2)*ZZ(k,3,l,1)-
     &  2*S3h0022(2)*ZZ(k,3,l,2)-2*S3h0023(2)*ZZ(k,3,l,3)
       aux002(3,3)=-(F(1)*S3313(3))-F(2)*S3323(3)-F(3)*S3333(3)+S34133(k
     &  )*Z(1,l)+S34233(k)*Z(2,l)+S34333(k)*Z(3,l)+(S30023(3)-S34133(1)-
     &  S34233(2)-S34333(3))*Z(k,l)-4*S3h0021(3)*ZZ(k,3,l,1)-4*S3h0022(3
     &  )*ZZ(k,3,l,2)-4*S3h0023(3)*ZZ(k,3,l,3)
       temp002(1,1)=I12Z*(aux002(1,1)+4*F(4)*temp001(1)+F(5)*temp2(1,1)+
     &  8*tempD40000*ZZ(k,1,l,1))
       temp002(2,1)=I12Z*(aux002(2,1)+2*(F(6)*temp001(1)+F(4)*temp001(2)
     &  )+F(5)*temp2(2,1)+4*tempD40000*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp002(2,2)=I12Z*(aux002(2,2)+4*F(6)*temp001(2)+F(5)*temp2(2,2)+
     &  8*tempD40000*ZZ(k,2,l,2))
       temp002(3,1)=I12Z*(aux002(3,1)+2*(F(7)*temp001(1)+F(4)*temp001(3)
     &  )+F(5)*temp2(3,1)+4*tempD40000*(ZZ(k,1,l,3)+ZZ(k,3,l,1)))
       temp002(3,2)=I12Z*(aux002(3,2)+2*(F(7)*temp001(2)+F(6)*temp001(3)
     &  )+F(5)*temp2(3,2)+4*tempD40000*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp002(3,3)=I12Z*(aux002(3,3)+4*F(7)*temp001(3)+F(5)*temp2(3,3)+
     &  8*tempD40000*ZZ(k,3,l,3))
       temp002(1,2)=temp002(2,1)
       temp002(1,3)=temp002(3,1)
       temp002(2,3)=temp002(3,2)
       aux3(1,1,1)=-(S34111(1)*Z(jj,1))-S34211(1)*Z(jj,2)-S34311(1)*Z(jj
     &  ,3)
       aux3(2,1,1)=-(S34121(1)*Z(jj,1))-S34221(1)*Z(jj,2)-S34321(1)*Z(jj
     &  ,3)
       aux3(2,2,1)=-(S34122(1)*Z(jj,1))-S34222(1)*Z(jj,2)-S34322(1)*Z(jj
     &  ,3)
       aux3(2,2,2)=-(S34122(2)*Z(jj,1))-S34222(2)*Z(jj,2)-S34322(2)*Z(jj
     &  ,3)
       aux3(3,1,1)=-(S34131(1)*Z(jj,1))-S34231(1)*Z(jj,2)-S34331(1)*Z(jj
     &  ,3)
       aux3(3,2,1)=-(S34132(1)*Z(jj,1))-S34232(1)*Z(jj,2)-S34332(1)*Z(jj
     &  ,3)
       aux3(3,2,2)=-(S34132(2)*Z(jj,1))-S34232(2)*Z(jj,2)-S34332(2)*Z(jj
     &  ,3)
       aux3(3,3,1)=-(S34133(1)*Z(jj,1))-S34233(1)*Z(jj,2)-S34333(1)*Z(jj
     &  ,3)
       aux3(3,3,2)=-(S34133(2)*Z(jj,1))-S34233(2)*Z(jj,2)-S34333(2)*Z(jj
     &  ,3)
       aux3(3,3,3)=-(S34133(3)*Z(jj,1))-S34233(3)*Z(jj,2)-S34333(3)*Z(jj
     &  ,3)
       temp3(1,1,1)=IX*(aux3(1,1,1)+6*temp002(1,1)*Z(jj,1))
       temp3(2,1,1)=IX*(aux3(2,1,1)+4*temp002(2,1)*Z(jj,1)+2*temp002(1,1
     &  )*Z(jj,2))
       temp3(2,2,1)=IX*(aux3(2,2,1)+2*temp002(2,2)*Z(jj,1)+4*temp002(2,1
     &  )*Z(jj,2))
       temp3(2,2,2)=IX*(aux3(2,2,2)+6*temp002(2,2)*Z(jj,2))
       temp3(3,1,1)=IX*(aux3(3,1,1)+4*temp002(3,1)*Z(jj,1)+2*temp002(1,1
     &  )*Z(jj,3))
       temp3(3,2,1)=IX*(aux3(3,2,1)+2*temp002(3,2)*Z(jj,1)+2*temp002(3,1
     &  )*Z(jj,2)+2*temp002(2,1)*Z(jj,3))
       temp3(3,2,2)=IX*(aux3(3,2,2)+4*temp002(3,2)*Z(jj,2)+2*temp002(2,2
     &  )*Z(jj,3))
       temp3(3,3,1)=IX*(aux3(3,3,1)+2*temp002(3,3)*Z(jj,1)+4*temp002(3,1
     &  )*Z(jj,3))
       temp3(3,3,2)=IX*(aux3(3,3,2)+2*temp002(3,3)*Z(jj,2)+4*temp002(3,2
     &  )*Z(jj,3))
       temp3(3,3,3)=IX*(aux3(3,3,3)+6*temp002(3,3)*Z(jj,3))
       temp3(1,1,2)=temp3(2,1,1)
       temp3(1,1,3)=temp3(3,1,1)
       temp3(1,2,1)=temp3(2,1,1)
       temp3(1,2,2)=temp3(2,2,1)
       temp3(1,2,3)=temp3(3,2,1)
       temp3(1,3,1)=temp3(3,1,1)
       temp3(1,3,2)=temp3(3,2,1)
       temp3(1,3,3)=temp3(3,3,1)
       temp3(2,1,2)=temp3(2,2,1)
       temp3(2,1,3)=temp3(3,2,1)
       temp3(2,2,3)=temp3(3,2,2)
       temp3(2,3,1)=temp3(3,2,1)
       temp3(2,3,2)=temp3(3,2,2)
       temp3(2,3,3)=temp3(3,3,2)
       temp3(3,1,2)=temp3(3,2,1)
       temp3(3,1,3)=temp3(3,3,1)
       temp3(3,2,3)=temp3(3,3,2)
c                Step2
       temp001(1)=I8Z*(aux001(1)+2*tempD400*F(4)+F(5)*temp1(1)-det4*temp
     &  3(1,k,l))
       temp001(2)=I8Z*(aux001(2)+2*tempD400*F(6)+F(5)*temp1(2)-det4*temp
     &  3(2,k,l))
       temp001(3)=I8Z*(aux001(3)+2*tempD400*F(7)+F(5)*temp1(3)-det4*temp
     &  3(3,k,l))
       temp2(1,1)=IX*(aux2(1,1)+det4*temp3(1,1,jj)+4*temp001(1)*Z(jj,1))
       temp2(2,1)=IX*(aux2(2,1)+det4*temp3(2,1,jj)+2*temp001(2)*Z(jj,1)+
     &  2*temp001(1)*Z(jj,2))
       temp2(3,1)=IX*(aux2(3,1)+det4*temp3(3,1,jj)+2*temp001(3)*Z(jj,1)+
     &  2*temp001(1)*Z(jj,3))
       temp2(2,2)=IX*(aux2(2,2)+det4*temp3(2,2,jj)+4*temp001(2)*Z(jj,2))
       temp2(3,2)=IX*(aux2(3,2)+det4*temp3(3,2,jj)+2*temp001(3)*Z(jj,2)+
     &  2*temp001(2)*Z(jj,3))
       temp2(3,3)=IX*(aux2(3,3)+det4*temp3(3,3,jj)+4*temp001(3)*Z(jj,3))
       temp2(1,2)=temp2(2,1)
       temp2(1,3)=temp2(3,1)
       temp2(2,3)=temp2(3,2)
c                Step3
       tempD400=I4Z*(auxD400+tempD40*F(5)-det4*temp2(k,l))
       temp1(1)=IX*(aux1(1)+det4*temp2(1,jj)+2*tempD400*Z(jj,1))
       temp1(2)=IX*(aux1(2)+det4*temp2(2,jj)+2*tempD400*Z(jj,2))
       temp1(3)=IX*(aux1(3)+det4*temp2(3,jj)+2*tempD400*Z(jj,3))
c                Step4
       tempD40=IX*(auxD40+det4*temp1(jj))
c FC       print*, 'auxD40',auxD40
c FC       print*, 'tempD40',tempD40
        if(printmy) then
       print*, "tempD40 3",tempD40,k,l,jj,IX,1d0/TX3
       endif
       if(order.eq.3) goto 500
c                Iteration4
c                Step1
       S3h00311(1)=Cij134(6,4)-Cij234(4,2)
       S3h00311(2)=Cij134(6,4)+Cij234(5,3)
       S3h00311(3)=Cij134(8,4)+Cij234(6,3)
       S3h00312(1)=Cij134(6,4)+Cij234(5,3)
       S3h00312(2)=Cij134(6,4)-Cij234(6,4)
       S3h00312(3)=Cij134(8,4)-Cij234(8,4)
       S3h00313(1)=Cij134(8,4)+Cij234(6,3)
       S3h00313(2)=Cij134(8,4)-Cij234(8,4)
       S3h00313(3)=Cij134(7,4)-Cij234(7,4)
       S3h00321(1)=Cij124(6,4)-Cij134(6,4)
       S3h00321(2)=Cij124(8,4)-Cij134(6,4)
       S3h00321(3)=Cij124(8,4)-Cij134(8,4)
       S3h00322(1)=Cij124(8,4)-Cij134(6,4)
       S3h00322(2)=Cij124(7,4)-Cij134(6,4)
       S3h00322(3)=Cij124(7,4)-Cij134(8,4)
       S3h00323(1)=Cij124(8,4)-Cij134(8,4)
       S3h00323(2)=Cij124(7,4)-Cij134(8,4)
       S3h00323(3)=Cij124(7,4)-Cij134(7,4)
       S3h00331(1)=Cij123(6,4)-Cij124(6,4)
       S3h00331(2)=Cij123(8,4)-Cij124(8,4)
       S3h00331(3)=-Cij124(8,4)
       S3h00332(1)=Cij123(8,4)-Cij124(8,4)
       S3h00332(2)=Cij123(7,4)-Cij124(7,4)
       S3h00332(3)=-Cij124(7,4)
       S3h00333(1)=-Cij124(8,4)
       S3h00333(2)=-Cij124(7,4)
       S3h00333(3)=-Cij124(7,4)
       S300311(1)=-2*C0234
       S300311(2)=2*Cij234(1,1)
       S300311(3)=2*Cij234(2,1)
       S300312(1)=2*Cij234(1,1)
       S300312(2)=-2*Cij234(1,2)
       S300312(3)=-2*Cij234(3,2)
       S300313(1)=2*Cij234(2,1)
       S300313(2)=-2*Cij234(3,2)
       S300313(3)=-2*Cij234(2,2)
       S300321(1)=2*Cij234(1,1)
       S300321(2)=-2*Cij234(1,2)
       S300321(3)=-2*Cij234(3,2)
       S300322(1)=-2*Cij234(1,2)
       S300322(2)=2*Cij234(1,3)
       S300322(3)=2*Cij234(3,3)
       S300323(1)=-2*Cij234(3,2)
       S300323(2)=2*Cij234(3,3)
       S300323(3)=2*Cij234(4,3)
       S300331(1)=2*Cij234(2,1)
       S300331(2)=-2*Cij234(3,2)
       S300331(3)=-2*Cij234(2,2)
       S300332(1)=-2*Cij234(3,2)
       S300332(2)=2*Cij234(3,3)
       S300332(3)=2*Cij234(4,3)
       S300333(1)=-2*Cij234(2,2)
       S300333(2)=2*Cij234(4,3)
       S300333(3)=2*Cij234(2,3)
       S300001(1)=-2*Cij234(4,2)
       S300001(2)=2*Cij234(5,3)
       S300001(3)=2*Cij234(6,3)
       S3h00001(1)=Cij134(9,4)-Cij234(9,4)
       S3h00001(2)=Cij124(9,4)-Cij134(9,4)
       S3h00001(3)=Cij123(9,4)-Cij124(9,4)
       aux00001(1)=-(F(1)*S3h0021(1))-F(2)*S3h0022(1)-F(3)*S3h0023(1)+S3
     &  h00311(k)*Z(1,l)+S3h00321(k)*Z(2,l)+S3h00331(k)*Z(3,l)+(Inv401+S
     &  300001(1)-S3h00311(1)-S3h00322(1)-S3h00333(1))*Z(k,l)-2*(S3h0000
     &  1(1)*ZZ(k,1,l,1)+S3h00001(2)*ZZ(k,1,l,2))-2*S3h00001(3)*ZZ(k,1,l
     &  ,3)
       aux00001(2)=-(F(1)*S3h0021(2))-F(2)*S3h0022(2)-F(3)*S3h0023(2)+S3
     &  h00312(k)*Z(1,l)+S3h00322(k)*Z(2,l)+S3h00332(k)*Z(3,l)+(Inv402+S
     &  300001(2)-S3h00312(1)-S3h00322(2)-S3h00333(2))*Z(k,l)-2*(S3h0000
     &  1(1)*ZZ(k,2,l,1)+S3h00001(2)*ZZ(k,2,l,2))-2*S3h00001(3)*ZZ(k,2,l
     &  ,3)
       aux00001(3)=-(F(1)*S3h0021(3))-F(2)*S3h0022(3)-F(3)*S3h0023(3)+S3
     &  h00313(k)*Z(1,l)+S3h00323(k)*Z(2,l)+S3h00333(k)*Z(3,l)+(Inv403+S
     &  300001(3)-S3h00313(1)-S3h00323(2)-S3h00333(3))*Z(k,l)-2*(S3h0000
     &  1(1)*ZZ(k,3,l,1)+S3h00001(2)*ZZ(k,3,l,2))-2*S3h00001(3)*ZZ(k,3,l
     &  ,3)
       temp00001(1)=I12Z*(aux00001(1)+2*tempD40000*F(4)+F(5)*temp001(1))
       temp00001(2)=I12Z*(aux00001(2)+2*tempD40000*F(6)+F(5)*temp001(2))
       temp00001(3)=I12Z*(aux00001(3)+2*tempD40000*F(7)+F(5)*temp001(3))
       S351111(1)=-C0234+Cij134(1,4)
       S351111(2)=Cij134(1,4)+Cij234(1,1)
       S351111(3)=Cij134(3,4)+Cij234(2,1)
       S351211(1)=Cij134(1,4)+Cij234(1,1)
       S351211(2)=Cij134(1,4)-Cij234(1,2)
       S351211(3)=Cij134(3,4)-Cij234(3,2)
       S351221(1)=Cij134(1,4)-Cij234(1,2)
       S351221(2)=Cij134(1,4)+Cij234(1,3)
       S351221(3)=Cij134(3,4)+Cij234(3,3)
       S351222(1)=Cij134(1,4)+Cij234(1,3)
       S351222(2)=Cij134(1,4)-Cij234(1,4)
       S351222(3)=Cij134(3,4)-Cij234(3,4)
       S351311(1)=Cij134(3,4)+Cij234(2,1)
       S351311(2)=Cij134(3,4)-Cij234(3,2)
       S351311(3)=Cij134(4,4)-Cij234(2,2)
       S351321(1)=Cij134(3,4)-Cij234(3,2)
       S351321(2)=Cij134(3,4)+Cij234(3,3)
       S351321(3)=Cij134(4,4)+Cij234(4,3)
       S351322(1)=Cij134(3,4)+Cij234(3,3)
       S351322(2)=Cij134(3,4)-Cij234(3,4)
       S351322(3)=Cij134(4,4)-Cij234(4,4)
       S351331(1)=Cij134(4,4)-Cij234(2,2)
       S351331(2)=Cij134(4,4)+Cij234(4,3)
       S351331(3)=Cij134(5,4)+Cij234(2,3)
       S351332(1)=Cij134(4,4)+Cij234(4,3)
       S351332(2)=Cij134(4,4)-Cij234(4,4)
       S351332(3)=Cij134(5,4)-Cij234(5,4)
       S351333(1)=Cij134(5,4)+Cij234(2,3)
       S351333(2)=Cij134(5,4)-Cij234(5,4)
       S351333(3)=Cij134(2,4)-Cij234(2,4)
       S352111(1)=Cij124(1,4)-Cij134(1,4)
       S352111(2)=Cij124(3,4)-Cij134(1,4)
       S352111(3)=Cij124(3,4)-Cij134(3,4)
       S352211(1)=Cij124(3,4)-Cij134(1,4)
       S352211(2)=Cij124(4,4)-Cij134(1,4)
       S352211(3)=Cij124(4,4)-Cij134(3,4)
       S352221(1)=Cij124(4,4)-Cij134(1,4)
       S352221(2)=Cij124(5,4)-Cij134(1,4)
       S352221(3)=Cij124(5,4)-Cij134(3,4)
       S352222(1)=Cij124(5,4)-Cij134(1,4)
       S352222(2)=Cij124(2,4)-Cij134(1,4)
       S352222(3)=Cij124(2,4)-Cij134(3,4)
       S352311(1)=Cij124(3,4)-Cij134(3,4)
       S352311(2)=Cij124(4,4)-Cij134(3,4)
       S352311(3)=Cij124(4,4)-Cij134(4,4)
       S352321(1)=Cij124(4,4)-Cij134(3,4)
       S352321(2)=Cij124(5,4)-Cij134(3,4)
       S352321(3)=Cij124(5,4)-Cij134(4,4)
       S352322(1)=Cij124(5,4)-Cij134(3,4)
       S352322(2)=Cij124(2,4)-Cij134(3,4)
       S352322(3)=Cij124(2,4)-Cij134(4,4)
       S352331(1)=Cij124(4,4)-Cij134(4,4)
       S352331(2)=Cij124(5,4)-Cij134(4,4)
       S352331(3)=Cij124(5,4)-Cij134(5,4)
       S352332(1)=Cij124(5,4)-Cij134(4,4)
       S352332(2)=Cij124(2,4)-Cij134(4,4)
       S352332(3)=Cij124(2,4)-Cij134(5,4)
       S352333(1)=Cij124(5,4)-Cij134(5,4)
       S352333(2)=Cij124(2,4)-Cij134(5,4)
       S352333(3)=Cij124(2,4)-Cij134(2,4)
       S353111(1)=Cij123(1,4)-Cij124(1,4)
       S353111(2)=Cij123(3,4)-Cij124(3,4)
       S353111(3)=-Cij124(3,4)
       S353211(1)=Cij123(3,4)-Cij124(3,4)
       S353211(2)=Cij123(4,4)-Cij124(4,4)
       S353211(3)=-Cij124(4,4)
       S353221(1)=Cij123(4,4)-Cij124(4,4)
       S353221(2)=Cij123(5,4)-Cij124(5,4)
       S353221(3)=-Cij124(5,4)
       S353222(1)=Cij123(5,4)-Cij124(5,4)
       S353222(2)=Cij123(2,4)-Cij124(2,4)
       S353222(3)=-Cij124(2,4)
       S353311(1)=-Cij124(3,4)
       S353311(2)=-Cij124(4,4)
       S353311(3)=-Cij124(4,4)
       S353321(1)=-Cij124(4,4)
       S353321(2)=-Cij124(5,4)
       S353321(3)=-Cij124(5,4)
       S353322(1)=-Cij124(5,4)
       S353322(2)=-Cij124(2,4)
       S353322(3)=-Cij124(2,4)
       S353331(1)=-Cij124(4,4)
       S353331(2)=-Cij124(5,4)
       S353331(3)=-Cij124(5,4)
       S353332(1)=-Cij124(5,4)
       S353332(2)=-Cij124(2,4)
       S353332(3)=-Cij124(2,4)
       S353333(1)=-Cij124(5,4)
       S353333(2)=-Cij124(2,4)
       S353333(3)=-Cij124(2,4)
       aux003(1,1,1)=-(F(1)*S34111(1))-F(2)*S34211(1)-F(3)*S34311(1)+S35
     &  1111(k)*Z(1,l)+S352111(k)*Z(2,l)+S353111(k)*Z(3,l)+(S300311(1)-S
     &  351111(1)-S352211(1)-S353311(1))*Z(k,l)-6*S3h00311(1)*ZZ(k,1,l,1
     &  )-6*S3h00321(1)*ZZ(k,1,l,2)-6*S3h00331(1)*ZZ(k,1,l,3)
       aux003(2,1,1)=-(F(1)*S34121(1))-F(2)*S34221(1)-F(3)*S34321(1)+S35
     &  1211(k)*Z(1,l)+S352211(k)*Z(2,l)+S353211(k)*Z(3,l)+(S300321(1)-S
     &  351211(1)-S352221(1)-S353321(1))*Z(k,l)-4*S3h00312(1)*ZZ(k,1,l,1
     &  )-4*S3h00322(1)*ZZ(k,1,l,2)-4*S3h00332(1)*ZZ(k,1,l,3)-2*S3h00311
     &  (1)*ZZ(k,2,l,1)-2*S3h00321(1)*ZZ(k,2,l,2)-2*S3h00331(1)*ZZ(k,2,l
     &  ,3)
       aux003(2,2,1)=-(F(1)*S34122(1))-F(2)*S34222(1)-F(3)*S34322(1)+S35
     &  1221(k)*Z(1,l)+S352221(k)*Z(2,l)+S353221(k)*Z(3,l)+(S300322(1)-S
     &  351221(1)-S352222(1)-S353322(1))*Z(k,l)-2*S3h00312(2)*ZZ(k,1,l,1
     &  )-2*S3h00322(2)*ZZ(k,1,l,2)-2*S3h00332(2)*ZZ(k,1,l,3)-4*S3h00312
     &  (1)*ZZ(k,2,l,1)-4*S3h00322(1)*ZZ(k,2,l,2)-4*S3h00332(1)*ZZ(k,2,l
     &  ,3)
       aux003(2,2,2)=-(F(1)*S34122(2))-F(2)*S34222(2)-F(3)*S34322(2)+S35
     &  1222(k)*Z(1,l)+S352222(k)*Z(2,l)+S353222(k)*Z(3,l)+(S300322(2)-S
     &  351222(1)-S352222(2)-S353322(2))*Z(k,l)-6*S3h00312(2)*ZZ(k,2,l,1
     &  )-6*S3h00322(2)*ZZ(k,2,l,2)-6*S3h00332(2)*ZZ(k,2,l,3)
       aux003(3,1,1)=-(F(1)*S34131(1))-F(2)*S34231(1)-F(3)*S34331(1)+S35
     &  1311(k)*Z(1,l)+S352311(k)*Z(2,l)+S353311(k)*Z(3,l)+(S300331(1)-S
     &  351311(1)-S352321(1)-S353331(1))*Z(k,l)-4*S3h00313(1)*ZZ(k,1,l,1
     &  )-4*S3h00323(1)*ZZ(k,1,l,2)-4*S3h00333(1)*ZZ(k,1,l,3)-2*S3h00311
     &  (1)*ZZ(k,3,l,1)-2*S3h00321(1)*ZZ(k,3,l,2)-2*S3h00331(1)*ZZ(k,3,l
     &  ,3)
       aux003(3,2,1)=-(F(1)*S34132(1))-F(2)*S34232(1)-F(3)*S34332(1)+S35
     &  1321(k)*Z(1,l)+S352321(k)*Z(2,l)+S353321(k)*Z(3,l)+(S300332(1)-S
     &  351321(1)-S352322(1)-S353332(1))*Z(k,l)-2*S3h00313(2)*ZZ(k,1,l,1
     &  )-2*S3h00323(2)*ZZ(k,1,l,2)-2*S3h00333(2)*ZZ(k,1,l,3)-2*S3h00313
     &  (1)*ZZ(k,2,l,1)-2*S3h00323(1)*ZZ(k,2,l,2)-2*S3h00333(1)*ZZ(k,2,l
     &  ,3)-2*S3h00312(1)*ZZ(k,3,l,1)-2*S3h00322(1)*ZZ(k,3,l,2)-2*S3h003
     &  32(1)*ZZ(k,3,l,3)
       aux003(3,2,2)=-(F(1)*S34132(2))-F(2)*S34232(2)-F(3)*S34332(2)+S35
     &  1322(k)*Z(1,l)+S352322(k)*Z(2,l)+S353322(k)*Z(3,l)+(S300332(2)-S
     &  351322(1)-S352322(2)-S353332(2))*Z(k,l)-4*S3h00313(2)*ZZ(k,2,l,1
     &  )-4*S3h00323(2)*ZZ(k,2,l,2)-4*S3h00333(2)*ZZ(k,2,l,3)-2*S3h00312
     &  (2)*ZZ(k,3,l,1)-2*S3h00322(2)*ZZ(k,3,l,2)-2*S3h00332(2)*ZZ(k,3,l
     &  ,3)
       aux003(3,3,1)=-(F(1)*S34133(1))-F(2)*S34233(1)-F(3)*S34333(1)+S35
     &  1331(k)*Z(1,l)+S352331(k)*Z(2,l)+S353331(k)*Z(3,l)+(S300333(1)-S
     &  351331(1)-S352332(1)-S353333(1))*Z(k,l)-2*S3h00313(3)*ZZ(k,1,l,1
     &  )-2*S3h00323(3)*ZZ(k,1,l,2)-2*S3h00333(3)*ZZ(k,1,l,3)-4*S3h00313
     &  (1)*ZZ(k,3,l,1)-4*S3h00323(1)*ZZ(k,3,l,2)-4*S3h00333(1)*ZZ(k,3,l
     &  ,3)
       aux003(3,3,2)=-(F(1)*S34133(2))-F(2)*S34233(2)-F(3)*S34333(2)+S35
     &  1332(k)*Z(1,l)+S352332(k)*Z(2,l)+S353332(k)*Z(3,l)+(S300333(2)-S
     &  351332(1)-S352332(2)-S353333(2))*Z(k,l)-2*S3h00313(3)*ZZ(k,2,l,1
     &  )-2*S3h00323(3)*ZZ(k,2,l,2)-2*S3h00333(3)*ZZ(k,2,l,3)-4*S3h00313
     &  (2)*ZZ(k,3,l,1)-4*S3h00323(2)*ZZ(k,3,l,2)-4*S3h00333(2)*ZZ(k,3,l
     &  ,3)
       aux003(3,3,3)=-(F(1)*S34133(3))-F(2)*S34233(3)-F(3)*S34333(3)+S35
     &  1333(k)*Z(1,l)+S352333(k)*Z(2,l)+S353333(k)*Z(3,l)+(S300333(3)-S
     &  351333(1)-S352333(2)-S353333(3))*Z(k,l)-6*S3h00313(3)*ZZ(k,3,l,1
     &  )-6*S3h00323(3)*ZZ(k,3,l,2)-6*S3h00333(3)*ZZ(k,3,l,3)
       temp003(1,1,1)=I16Z*(aux003(1,1,1)+6*F(4)*temp002(1,1)+F(5)*temp3
     &  (1,1,1)+24*temp00001(1)*ZZ(k,1,l,1))
       temp003(2,1,1)=I16Z*(aux003(2,1,1)+2*F(6)*temp002(1,1)+4*F(4)*tem
     &  p002(2,1)+F(5)*temp3(2,1,1)+8*(temp00001(2)*ZZ(k,1,l,1)+temp0000
     &  1(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))))
       temp003(2,2,1)=I16Z*(aux003(2,2,1)+4*F(6)*temp002(2,1)+2*F(4)*tem
     &  p002(2,2)+F(5)*temp3(2,2,1)+8*(temp00001(2)*(ZZ(k,1,l,2)+ZZ(k,2,
     &  l,1))+temp00001(1)*ZZ(k,2,l,2)))
       temp003(2,2,2)=I16Z*(aux003(2,2,2)+6*F(6)*temp002(2,2)+F(5)*temp3
     &  (2,2,2)+24*temp00001(2)*ZZ(k,2,l,2))
       temp003(3,1,1)=I16Z*(aux003(3,1,1)+2*F(7)*temp002(1,1)+4*F(4)*tem
     &  p002(3,1)+F(5)*temp3(3,1,1)+8*(temp00001(3)*ZZ(k,1,l,1)+temp0000
     &  1(1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))))
       temp003(3,2,1)=I16Z*(aux003(3,2,1)+2*F(7)*temp002(2,1)+2*F(6)*tem
     &  p002(3,1)+2*F(4)*temp002(3,2)+F(5)*temp3(3,2,1)+4*(temp00001(3)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00001(2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))
     &  )+4*temp00001(1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp003(3,2,2)=I16Z*(aux003(3,2,2)+2*F(7)*temp002(2,2)+4*F(6)*tem
     &  p002(3,2)+F(5)*temp3(3,2,2)+8*(temp00001(3)*ZZ(k,2,l,2)+temp0000
     &  1(2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))))
       temp003(3,3,1)=I16Z*(aux003(3,3,1)+4*F(7)*temp002(3,1)+2*F(4)*tem
     &  p002(3,3)+F(5)*temp3(3,3,1)+8*(temp00001(3)*(ZZ(k,1,l,3)+ZZ(k,3,
     &  l,1))+temp00001(1)*ZZ(k,3,l,3)))
       temp003(3,3,2)=I16Z*(aux003(3,3,2)+4*F(7)*temp002(3,2)+2*F(6)*tem
     &  p002(3,3)+F(5)*temp3(3,3,2)+8*(temp00001(3)*(ZZ(k,2,l,3)+ZZ(k,3,
     &  l,2))+temp00001(2)*ZZ(k,3,l,3)))
       temp003(3,3,3)=I16Z*(aux003(3,3,3)+6*F(7)*temp002(3,3)+F(5)*temp3
     &  (3,3,3)+24*temp00001(3)*ZZ(k,3,l,3))
       temp003(1,1,2)=temp003(2,1,1)
       temp003(1,1,3)=temp003(3,1,1)
       temp003(1,2,1)=temp003(2,1,1)
       temp003(1,2,2)=temp003(2,2,1)
       temp003(1,2,3)=temp003(3,2,1)
       temp003(1,3,1)=temp003(3,1,1)
       temp003(1,3,2)=temp003(3,2,1)
       temp003(1,3,3)=temp003(3,3,1)
       temp003(2,1,2)=temp003(2,2,1)
       temp003(2,1,3)=temp003(3,2,1)
       temp003(2,2,3)=temp003(3,2,2)
       temp003(2,3,1)=temp003(3,2,1)
       temp003(2,3,2)=temp003(3,2,2)
       temp003(2,3,3)=temp003(3,3,2)
       temp003(3,1,2)=temp003(3,2,1)
       temp003(3,1,3)=temp003(3,3,1)
       temp003(3,2,3)=temp003(3,3,2)
       aux41(1,1,1)=-(S351111(1)*Z(jj,1))-S352111(1)*Z(jj,2)-S353111(1)*
     &  Z(jj,3)
       aux42(1,1,1)=-(S351211(1)*Z(jj,1))-S352211(1)*Z(jj,2)-S353211(1)*
     &  Z(jj,3)
       aux42(2,1,1)=-(S351221(1)*Z(jj,1))-S352221(1)*Z(jj,2)-S353221(1)*
     &  Z(jj,3)
       aux42(2,2,1)=-(S351222(1)*Z(jj,1))-S352222(1)*Z(jj,2)-S353222(1)*
     &  Z(jj,3)
       aux42(2,2,2)=-(S351222(2)*Z(jj,1))-S352222(2)*Z(jj,2)-S353222(2)*
     &  Z(jj,3)
       aux43(1,1,1)=-(S351311(1)*Z(jj,1))-S352311(1)*Z(jj,2)-S353311(1)*
     &  Z(jj,3)
       aux43(2,1,1)=-(S351321(1)*Z(jj,1))-S352321(1)*Z(jj,2)-S353321(1)*
     &  Z(jj,3)
       aux43(2,2,1)=-(S351322(1)*Z(jj,1))-S352322(1)*Z(jj,2)-S353322(1)*
     &  Z(jj,3)
       aux43(2,2,2)=-(S351322(2)*Z(jj,1))-S352322(2)*Z(jj,2)-S353322(2)*
     &  Z(jj,3)
       aux43(3,1,1)=-(S351331(1)*Z(jj,1))-S352331(1)*Z(jj,2)-S353331(1)*
     &  Z(jj,3)
       aux43(3,2,1)=-(S351332(1)*Z(jj,1))-S352332(1)*Z(jj,2)-S353332(1)*
     &  Z(jj,3)
       aux43(3,2,2)=-(S351332(2)*Z(jj,1))-S352332(2)*Z(jj,2)-S353332(2)*
     &  Z(jj,3)
       aux43(3,3,1)=-(S351333(1)*Z(jj,1))-S352333(1)*Z(jj,2)-S353333(1)*
     &  Z(jj,3)
       aux43(3,3,2)=-(S351333(2)*Z(jj,1))-S352333(2)*Z(jj,2)-S353333(2)*
     &  Z(jj,3)
       aux43(3,3,3)=-(S351333(3)*Z(jj,1))-S352333(3)*Z(jj,2)-S353333(3)*
     &  Z(jj,3)
       temp41(1,1,1)=IX*(aux41(1,1,1)+8*temp003(1,1,1)*Z(jj,1))
       temp42(1,1,1)=IX*(aux42(1,1,1)+6*temp003(2,1,1)*Z(jj,1)+2*temp003
     &  (1,1,1)*Z(jj,2))
       temp42(2,1,1)=IX*(aux42(2,1,1)+4*(temp003(2,2,1)*Z(jj,1)+temp003(
     &  2,1,1)*Z(jj,2)))
       temp42(2,2,1)=IX*(aux42(2,2,1)+2*temp003(2,2,2)*Z(jj,1)+6*temp003
     &  (2,2,1)*Z(jj,2))
       temp42(2,2,2)=IX*(aux42(2,2,2)+8*temp003(2,2,2)*Z(jj,2))
       temp43(1,1,1)=IX*(aux43(1,1,1)+6*temp003(3,1,1)*Z(jj,1)+2*temp003
     &  (1,1,1)*Z(jj,3))
       temp43(2,1,1)=IX*(aux43(2,1,1)+4*temp003(3,2,1)*Z(jj,1)+2*(temp00
     &  3(3,1,1)*Z(jj,2)+temp003(2,1,1)*Z(jj,3)))
       temp43(2,2,1)=IX*(aux43(2,2,1)+4*temp003(3,2,1)*Z(jj,2)+2*(temp00
     &  3(3,2,2)*Z(jj,1)+temp003(2,2,1)*Z(jj,3)))
       temp43(2,2,2)=IX*(aux43(2,2,2)+6*temp003(3,2,2)*Z(jj,2)+2*temp003
     &  (2,2,2)*Z(jj,3))
       temp43(3,1,1)=IX*(aux43(3,1,1)+4*(temp003(3,3,1)*Z(jj,1)+temp003(
     &  3,1,1)*Z(jj,3)))
       temp43(3,2,1)=IX*(aux43(3,2,1)+2*(temp003(3,3,2)*Z(jj,1)+temp003(
     &  3,3,1)*Z(jj,2))+4*temp003(3,2,1)*Z(jj,3))
       temp43(3,2,2)=IX*(aux43(3,2,2)+4*(temp003(3,3,2)*Z(jj,2)+temp003(
     &  3,2,2)*Z(jj,3)))
       temp43(3,3,1)=IX*(aux43(3,3,1)+2*temp003(3,3,3)*Z(jj,1)+6*temp003
     &  (3,3,1)*Z(jj,3))
       temp43(3,3,2)=IX*(aux43(3,3,2)+2*temp003(3,3,3)*Z(jj,2)+6*temp003
     &  (3,3,2)*Z(jj,3))
       temp43(3,3,3)=IX*(aux43(3,3,3)+8*temp003(3,3,3)*Z(jj,3))
       temp41(1,1,2)=temp42(1,1,1)
       temp41(1,1,3)=temp43(1,1,1)
       temp41(1,2,1)=temp42(1,1,1)
       temp41(1,2,2)=temp42(2,1,1)
       temp41(1,2,3)=temp43(2,1,1)
       temp41(1,3,1)=temp43(1,1,1)
       temp41(1,3,2)=temp43(2,1,1)
       temp41(1,3,3)=temp43(3,1,1)
       temp42(1,1,2)=temp42(2,1,1)
       temp42(1,1,3)=temp43(2,1,1)
       temp42(1,2,1)=temp42(2,1,1)
       temp42(1,2,2)=temp42(2,2,1)
       temp42(1,2,3)=temp43(2,2,1)
       temp42(1,3,1)=temp43(2,1,1)
       temp42(1,3,2)=temp43(2,2,1)
       temp42(1,3,3)=temp43(3,2,1)
       temp42(2,1,2)=temp42(2,2,1)
       temp42(2,1,3)=temp43(2,2,1)
       temp42(2,2,3)=temp43(2,2,2)
       temp42(2,3,1)=temp43(2,2,1)
       temp42(2,3,2)=temp43(2,2,2)
       temp42(2,3,3)=temp43(3,2,2)
       temp43(1,1,2)=temp43(2,1,1)
       temp43(1,1,3)=temp43(3,1,1)
       temp43(1,2,1)=temp43(2,1,1)
       temp43(1,2,2)=temp43(2,2,1)
       temp43(1,2,3)=temp43(3,2,1)
       temp43(1,3,1)=temp43(3,1,1)
       temp43(1,3,2)=temp43(3,2,1)
       temp43(1,3,3)=temp43(3,3,1)
       temp43(2,1,2)=temp43(2,2,1)
       temp43(2,1,3)=temp43(3,2,1)
       temp43(2,2,3)=temp43(3,2,2)
       temp43(2,3,1)=temp43(3,2,1)
       temp43(2,3,2)=temp43(3,2,2)
       temp43(2,3,3)=temp43(3,3,2)
       temp43(3,1,2)=temp43(3,2,1)
       temp43(3,1,3)=temp43(3,3,1)
       temp43(3,2,3)=temp43(3,3,2)
c                Step2
       tempD40000=I8Z*(auxD40000+tempD400*F(5)-det4*temp002(k,l))
       temp002(1,1)=I12Z*(aux002(1,1)+4*F(4)*temp001(1)+F(5)*temp2(1,1)-
     &  det4*temp41(1,k,l)+8*tempD40000*ZZ(k,1,l,1))
       temp002(2,1)=I12Z*(aux002(2,1)+2*(F(6)*temp001(1)+F(4)*temp001(2)
     &  )+F(5)*temp2(2,1)-det4*temp42(1,k,l)+4*tempD40000*(ZZ(k,1,l,2)+Z
     &  Z(k,2,l,1)))
       temp002(2,2)=I12Z*(aux002(2,2)+4*F(6)*temp001(2)+F(5)*temp2(2,2)-
     &  det4*temp42(2,k,l)+8*tempD40000*ZZ(k,2,l,2))
       temp002(3,1)=I12Z*(aux002(3,1)+2*(F(7)*temp001(1)+F(4)*temp001(3)
     &  )+F(5)*temp2(3,1)-det4*temp43(1,k,l)+4*tempD40000*(ZZ(k,1,l,3)+Z
     &  Z(k,3,l,1)))
       temp002(3,2)=I12Z*(aux002(3,2)+2*(F(7)*temp001(2)+F(6)*temp001(3)
     &  )+F(5)*temp2(3,2)-det4*temp43(2,k,l)+4*tempD40000*(ZZ(k,2,l,3)+Z
     &  Z(k,3,l,2)))
       temp002(3,3)=I12Z*(aux002(3,3)+4*F(7)*temp001(3)+F(5)*temp2(3,3)-
     &  det4*temp43(3,k,l)+8*tempD40000*ZZ(k,3,l,3))
       temp002(1,2)=temp002(2,1)
       temp002(1,3)=temp002(3,1)
       temp002(2,3)=temp002(3,2)
       temp3(1,1,1)=IX*(aux3(1,1,1)+det4*temp41(1,1,jj)+6*temp002(1,1)*Z
     &  (jj,1))
       temp3(2,1,1)=IX*(aux3(2,1,1)+det4*temp42(1,1,jj)+4*temp002(2,1)*Z
     &  (jj,1)+2*temp002(1,1)*Z(jj,2))
       temp3(2,2,1)=IX*(aux3(2,2,1)+det4*temp42(2,1,jj)+2*temp002(2,2)*Z
     &  (jj,1)+4*temp002(2,1)*Z(jj,2))
       temp3(2,2,2)=IX*(aux3(2,2,2)+det4*temp42(2,2,jj)+6*temp002(2,2)*Z
     &  (jj,2))
       temp3(3,1,1)=IX*(aux3(3,1,1)+det4*temp43(1,1,jj)+4*temp002(3,1)*Z
     &  (jj,1)+2*temp002(1,1)*Z(jj,3))
       temp3(3,2,1)=IX*(aux3(3,2,1)+det4*temp43(2,1,jj)+2*(temp002(3,2)*
     &  Z(jj,1)+temp002(3,1)*Z(jj,2))+2*temp002(2,1)*Z(jj,3))
       temp3(3,2,2)=IX*(aux3(3,2,2)+det4*temp43(2,2,jj)+4*temp002(3,2)*Z
     &  (jj,2)+2*temp002(2,2)*Z(jj,3))
       temp3(3,3,1)=IX*(aux3(3,3,1)+det4*temp43(3,1,jj)+2*temp002(3,3)*Z
     &  (jj,1)+4*temp002(3,1)*Z(jj,3))
       temp3(3,3,2)=IX*(aux3(3,3,2)+det4*temp43(3,2,jj)+2*temp002(3,3)*Z
     &  (jj,2)+4*temp002(3,2)*Z(jj,3))
       temp3(3,3,3)=IX*(aux3(3,3,3)+det4*temp43(3,3,jj)+6*temp002(3,3)*Z
     &  (jj,3))
       temp3(1,1,2)=temp3(2,1,1)
       temp3(1,1,3)=temp3(3,1,1)
       temp3(1,2,1)=temp3(2,1,1)
       temp3(1,2,2)=temp3(2,2,1)
       temp3(1,2,3)=temp3(3,2,1)
       temp3(1,3,1)=temp3(3,1,1)
       temp3(1,3,2)=temp3(3,2,1)
       temp3(1,3,3)=temp3(3,3,1)
       temp3(2,1,2)=temp3(2,2,1)
       temp3(2,1,3)=temp3(3,2,1)
       temp3(2,2,3)=temp3(3,2,2)
       temp3(2,3,1)=temp3(3,2,1)
       temp3(2,3,2)=temp3(3,2,2)
       temp3(2,3,3)=temp3(3,3,2)
       temp3(3,1,2)=temp3(3,2,1)
       temp3(3,1,3)=temp3(3,3,1)
       temp3(3,2,3)=temp3(3,3,2)
c                Step3
       temp001(1)=I8Z*(aux001(1)+2*tempD400*F(4)+F(5)*temp1(1)-det4*temp
     &  3(1,k,l))
       temp001(2)=I8Z*(aux001(2)+2*tempD400*F(6)+F(5)*temp1(2)-det4*temp
     &  3(2,k,l))
       temp001(3)=I8Z*(aux001(3)+2*tempD400*F(7)+F(5)*temp1(3)-det4*temp
     &  3(3,k,l))
       temp2(1,1)=IX*(aux2(1,1)+det4*temp3(1,1,jj)+4*temp001(1)*Z(jj,1))
       temp2(2,1)=IX*(aux2(2,1)+det4*temp3(2,1,jj)+2*temp001(2)*Z(jj,1)+
     &  2*temp001(1)*Z(jj,2))
       temp2(2,2)=IX*(aux2(2,2)+det4*temp3(2,2,jj)+4*temp001(2)*Z(jj,2))
       temp2(3,1)=IX*(aux2(3,1)+det4*temp3(3,1,jj)+2*temp001(3)*Z(jj,1)+
     &  2*temp001(1)*Z(jj,3))
       temp2(3,2)=IX*(aux2(3,2)+det4*temp3(3,2,jj)+2*temp001(3)*Z(jj,2)+
     &  2*temp001(2)*Z(jj,3))
       temp2(3,3)=IX*(aux2(3,3)+det4*temp3(3,3,jj)+4*temp001(3)*Z(jj,3))
       temp2(1,2)=temp2(2,1)
       temp2(1,3)=temp2(3,1)
       temp2(2,3)=temp2(3,2)
c                Step4
       tempD400=I4Z*(auxD400+tempD40*F(5)-det4*temp2(k,l))
       temp1(1)=IX*(aux1(1)+det4*temp2(1,jj)+2*tempD400*Z(jj,1))
       temp1(2)=IX*(aux1(2)+det4*temp2(2,jj)+2*tempD400*Z(jj,2))
       temp1(3)=IX*(aux1(3)+det4*temp2(3,jj)+2*tempD400*Z(jj,3))
c                Step5
       tempD40=IX*(auxD40+det4*temp1(jj))
        if(printmy) then
       print*, "tempD40 4",tempD40,k,l,jj,IX,1d0/TX3
        endif
       if(order.eq.4) goto 500
c                Iteration5
c                Step1
       S3000000=2*Cij234(9,4)
       S3h000021(1)=Cij134(11,5)+Cij234(9,4)
       S3h000021(2)=Cij134(11,5)-Cij234(11,5)
       S3h000021(3)=Cij134(12,5)-Cij234(12,5)
       S3h000022(1)=Cij124(11,5)-Cij134(11,5)
       S3h000022(2)=Cij124(12,5)-Cij134(11,5)
       S3h000022(3)=Cij124(12,5)-Cij134(12,5)
       S3h000023(1)=Cij123(11,5)-Cij124(11,5)
       S3h000023(2)=Cij123(12,5)-Cij124(12,5)
       S3h000023(3)=-Cij124(12,5)
       auxD4000000=-(F(1)*S3h00001(1))-F(2)*S3h00001(2)-F(3)*S3h00001(3)
     &  +S3h000021(k)*Z(1,l)+S3h000022(k)*Z(2,l)+S3h000023(k)*Z(3,l)+(In
     &  v60+S3000000-S3h000021(1)-S3h000022(2)-S3h000023(3))*Z(k,l)
       tempD4000000=I12Z*(auxD4000000+tempD40000*F(5))
       S3000021(1)=2*Cij234(4,2)
       S3000021(2)=-2*Cij234(5,3)
       S3000021(3)=-2*Cij234(6,3)
       S3000022(1)=-2*Cij234(5,3)
       S3000022(2)=2*Cij234(6,4)
       S3000022(3)=2*Cij234(8,4)
       S3000023(1)=-2*Cij234(6,3)
       S3000023(2)=2*Cij234(8,4)
       S3000023(3)=2*Cij234(7,4)
       S3h004111(1)=Cij134(7,5)+Cij234(4,2)
       S3h004111(2)=Cij134(7,5)-Cij234(5,3)
       S3h004111(3)=Cij134(9,5)-Cij234(6,3)
       S3h004121(1)=Cij134(7,5)-Cij234(5,3)
       S3h004121(2)=Cij134(7,5)+Cij234(6,4)
       S3h004121(3)=Cij134(9,5)+Cij234(8,4)
       S3h004122(1)=Cij134(7,5)+Cij234(6,4)
       S3h004122(2)=Cij134(7,5)-Cij234(7,5)
       S3h004122(3)=Cij134(9,5)-Cij234(9,5)
       S3h004131(1)=Cij134(9,5)-Cij234(6,3)
       S3h004131(2)=Cij134(9,5)+Cij234(8,4)
       S3h004131(3)=Cij134(10,5)+Cij234(7,4)
       S3h004132(1)=Cij134(9,5)+Cij234(8,4)
       S3h004132(2)=Cij134(9,5)-Cij234(9,5)
       S3h004132(3)=Cij134(10,5)-Cij234(10,5)
       S3h004133(1)=Cij134(10,5)+Cij234(7,4)
       S3h004133(2)=Cij134(10,5)-Cij234(10,5)
       S3h004133(3)=Cij134(8,5)-Cij234(8,5)
       S3h004211(1)=Cij124(7,5)-Cij134(7,5)
       S3h004211(2)=Cij124(9,5)-Cij134(7,5)
       S3h004211(3)=Cij124(9,5)-Cij134(9,5)
       S3h004221(1)=Cij124(9,5)-Cij134(7,5)
       S3h004221(2)=Cij124(10,5)-Cij134(7,5)
       S3h004221(3)=Cij124(10,5)-Cij134(9,5)
       S3h004222(1)=Cij124(10,5)-Cij134(7,5)
       S3h004222(2)=Cij124(8,5)-Cij134(7,5)
       S3h004222(3)=Cij124(8,5)-Cij134(9,5)
       S3h004231(1)=Cij124(9,5)-Cij134(9,5)
       S3h004231(2)=Cij124(10,5)-Cij134(9,5)
       S3h004231(3)=Cij124(10,5)-Cij134(10,5)
       S3h004232(1)=Cij124(10,5)-Cij134(9,5)
       S3h004232(2)=Cij124(8,5)-Cij134(9,5)
       S3h004232(3)=Cij124(8,5)-Cij134(10,5)
       S3h004233(1)=Cij124(10,5)-Cij134(10,5)
       S3h004233(2)=Cij124(8,5)-Cij134(10,5)
       S3h004233(3)=Cij124(8,5)-Cij134(8,5)
       S3h004311(1)=Cij123(7,5)-Cij124(7,5)
       S3h004311(2)=Cij123(9,5)-Cij124(9,5)
       S3h004311(3)=-Cij124(9,5)
       S3h004321(1)=Cij123(9,5)-Cij124(9,5)
       S3h004321(2)=Cij123(10,5)-Cij124(10,5)
       S3h004321(3)=-Cij124(10,5)
       S3h004322(1)=Cij123(10,5)-Cij124(10,5)
       S3h004322(2)=Cij123(8,5)-Cij124(8,5)
       S3h004322(3)=-Cij124(8,5)
       S3h004331(1)=-Cij124(9,5)
       S3h004331(2)=-Cij124(10,5)
       S3h004331(3)=-Cij124(10,5)
       S3h004332(1)=-Cij124(10,5)
       S3h004332(2)=-Cij124(8,5)
       S3h004332(3)=-Cij124(8,5)
       S3h004333(1)=-Cij124(10,5)
       S3h004333(2)=-Cij124(8,5)
       S3h004333(3)=-Cij124(8,5)
       aux00002(1,1)=-(F(1)*S3h00311(1))-F(2)*S3h00321(1)-F(3)*S3h00331(
     &  1)+S3h004111(k)*Z(1,l)+S3h004211(k)*Z(2,l)+S3h004311(k)*Z(3,l)+(
     &  Inv4011+S3000021(1)-S3h004111(1)-S3h004221(1)-S3h004331(1))*Z(k,
     &  l)-4*S3h000021(1)*ZZ(k,1,l,1)-4*S3h000022(1)*ZZ(k,1,l,2)-4*S3h00
     &  0023(1)*ZZ(k,1,l,3)
       aux00002(2,1)=-(F(1)*S3h00312(1))-F(2)*S3h00322(1)-F(3)*S3h00332(
     &  1)+S3h004121(k)*Z(1,l)+S3h004221(k)*Z(2,l)+S3h004321(k)*Z(3,l)+(
     &  Inv4021+S3000022(1)-S3h004121(1)-S3h004222(1)-S3h004332(1))*Z(k,
     &  l)-2*S3h000021(2)*ZZ(k,1,l,1)-2*S3h000022(2)*ZZ(k,1,l,2)-2*S3h00
     &  0023(2)*ZZ(k,1,l,3)-2*S3h000021(1)*ZZ(k,2,l,1)-2*S3h000022(1)*ZZ
     &  (k,2,l,2)-2*S3h000023(1)*ZZ(k,2,l,3)
       aux00002(2,2)=-(F(1)*S3h00312(2))-F(2)*S3h00322(2)-F(3)*S3h00332(
     &  2)+S3h004122(k)*Z(1,l)+S3h004222(k)*Z(2,l)+S3h004322(k)*Z(3,l)+(
     &  Inv4022+S3000022(2)-S3h004122(1)-S3h004222(2)-S3h004332(2))*Z(k,
     &  l)-4*S3h000021(2)*ZZ(k,2,l,1)-4*S3h000022(2)*ZZ(k,2,l,2)-4*S3h00
     &  0023(2)*ZZ(k,2,l,3)
       aux00002(3,1)=-(F(1)*S3h00313(1))-F(2)*S3h00323(1)-F(3)*S3h00333(
     &  1)+S3h004131(k)*Z(1,l)+S3h004231(k)*Z(2,l)+S3h004331(k)*Z(3,l)+(
     &  Inv4031+S3000023(1)-S3h004131(1)-S3h004232(1)-S3h004333(1))*Z(k,
     &  l)-2*S3h000021(3)*ZZ(k,1,l,1)-2*S3h000022(3)*ZZ(k,1,l,2)-2*S3h00
     &  0023(3)*ZZ(k,1,l,3)-2*S3h000021(1)*ZZ(k,3,l,1)-2*S3h000022(1)*ZZ
     &  (k,3,l,2)-2*S3h000023(1)*ZZ(k,3,l,3)
       aux00002(3,2)=-(F(1)*S3h00313(2))-F(2)*S3h00323(2)-F(3)*S3h00333(
     &  2)+S3h004132(k)*Z(1,l)+S3h004232(k)*Z(2,l)+S3h004332(k)*Z(3,l)+(
     &  Inv4032+S3000023(2)-S3h004132(1)-S3h004232(2)-S3h004333(2))*Z(k,
     &  l)-2*S3h000021(3)*ZZ(k,2,l,1)-2*S3h000022(3)*ZZ(k,2,l,2)-2*S3h00
     &  0023(3)*ZZ(k,2,l,3)-2*S3h000021(2)*ZZ(k,3,l,1)-2*S3h000022(2)*ZZ
     &  (k,3,l,2)-2*S3h000023(2)*ZZ(k,3,l,3)
       aux00002(3,3)=-(F(1)*S3h00313(3))-F(2)*S3h00323(3)-F(3)*S3h00333(
     &  3)+S3h004133(k)*Z(1,l)+S3h004233(k)*Z(2,l)+S3h004333(k)*Z(3,l)+(
     &  Inv4033+S3000023(3)-S3h004133(1)-S3h004233(2)-S3h004333(3))*Z(k,
     &  l)-4*S3h000021(3)*ZZ(k,3,l,1)-4*S3h000022(3)*ZZ(k,3,l,2)-4*S3h00
     &  0023(3)*ZZ(k,3,l,3)
       temp00002(1,1)=I16Z*(aux00002(1,1)+4*F(4)*temp00001(1)+F(5)*temp0
     &  02(1,1)+8*tempD4000000*ZZ(k,1,l,1))
       temp00002(2,1)=I16Z*(aux00002(2,1)+2*(F(6)*temp00001(1)+F(4)*temp
     &  00001(2))+F(5)*temp002(2,1)+4*tempD4000000*(ZZ(k,1,l,2)+ZZ(k,2,l
     &  ,1)))
       temp00002(2,2)=I16Z*(aux00002(2,2)+4*F(6)*temp00001(2)+F(5)*temp0
     &  02(2,2)+8*tempD4000000*ZZ(k,2,l,2))
       temp00002(3,1)=I16Z*(aux00002(3,1)+2*(F(7)*temp00001(1)+F(4)*temp
     &  00001(3))+F(5)*temp002(3,1)+4*tempD4000000*(ZZ(k,1,l,3)+ZZ(k,3,l
     &  ,1)))
       temp00002(3,2)=I16Z*(aux00002(3,2)+2*(F(7)*temp00001(2)+F(6)*temp
     &  00001(3))+F(5)*temp002(3,2)+4*tempD4000000*(ZZ(k,2,l,3)+ZZ(k,3,l
     &  ,2)))
       temp00002(3,3)=I16Z*(aux00002(3,3)+4*F(7)*temp00001(3)+F(5)*temp0
     &  02(3,3)+8*tempD4000000*ZZ(k,3,l,3))
       temp00002(1,2)=temp00002(2,1)
       temp00002(1,3)=temp00002(3,1)
       temp00002(2,3)=temp00002(3,2)
       S3004111(1)=2*C0234
       S3004111(2)=-2*Cij234(1,1)
       S3004111(3)=-2*Cij234(2,1)
       S3004121(1)=-2*Cij234(1,1)
       S3004121(2)=2*Cij234(1,2)
       S3004121(3)=2*Cij234(3,2)
       S3004122(1)=2*Cij234(1,2)
       S3004122(2)=-2*Cij234(1,3)
       S3004122(3)=-2*Cij234(3,3)
       S3004131(1)=-2*Cij234(2,1)
       S3004131(2)=2*Cij234(3,2)
       S3004131(3)=2*Cij234(2,2)
       S3004132(1)=2*Cij234(3,2)
       S3004132(2)=-2*Cij234(3,3)
       S3004132(3)=-2*Cij234(4,3)
       S3004133(1)=2*Cij234(2,2)
       S3004133(2)=-2*Cij234(4,3)
       S3004133(3)=-2*Cij234(2,3)
       S3004211(1)=-2*Cij234(1,1)
       S3004211(2)=2*Cij234(1,2)
       S3004211(3)=2*Cij234(3,2)
       S3004221(1)=2*Cij234(1,2)
       S3004221(2)=-2*Cij234(1,3)
       S3004221(3)=-2*Cij234(3,3)
       S3004222(1)=-2*Cij234(1,3)
       S3004222(2)=2*Cij234(1,4)
       S3004222(3)=2*Cij234(3,4)
       S3004231(1)=2*Cij234(3,2)
       S3004231(2)=-2*Cij234(3,3)
       S3004231(3)=-2*Cij234(4,3)
       S3004232(1)=-2*Cij234(3,3)
       S3004232(2)=2*Cij234(3,4)
       S3004232(3)=2*Cij234(4,4)
       S3004233(1)=-2*Cij234(4,3)
       S3004233(2)=2*Cij234(4,4)
       S3004233(3)=2*Cij234(5,4)
       S3004311(1)=-2*Cij234(2,1)
       S3004311(2)=2*Cij234(3,2)
       S3004311(3)=2*Cij234(2,2)
       S3004321(1)=2*Cij234(3,2)
       S3004321(2)=-2*Cij234(3,3)
       S3004321(3)=-2*Cij234(4,3)
       S3004322(1)=-2*Cij234(3,3)
       S3004322(2)=2*Cij234(3,4)
       S3004322(3)=2*Cij234(4,4)
       S3004331(1)=2*Cij234(2,2)
       S3004331(2)=-2*Cij234(4,3)
       S3004331(3)=-2*Cij234(2,3)
       S3004332(1)=-2*Cij234(4,3)
       S3004332(2)=2*Cij234(4,4)
       S3004332(3)=2*Cij234(5,4)
       S3004333(1)=-2*Cij234(2,3)
       S3004333(2)=2*Cij234(5,4)
       S3004333(3)=2*Cij234(2,4)
       S3611111(1)=C0234+Cij134(1,5)
       S3611111(2)=Cij134(1,5)-Cij234(1,1)
       S3611111(3)=Cij134(3,5)-Cij234(2,1)
       S3612111(1)=Cij134(1,5)-Cij234(1,1)
       S3612111(2)=Cij134(1,5)+Cij234(1,2)
       S3612111(3)=Cij134(3,5)+Cij234(3,2)
       S3612211(1)=Cij134(1,5)+Cij234(1,2)
       S3612211(2)=Cij134(1,5)-Cij234(1,3)
       S3612211(3)=Cij134(3,5)-Cij234(3,3)
       S3612221(1)=Cij134(1,5)-Cij234(1,3)
       S3612221(2)=Cij134(1,5)+Cij234(1,4)
       S3612221(3)=Cij134(3,5)+Cij234(3,4)
       S3612222(1)=Cij134(1,5)+Cij234(1,4)
       S3612222(2)=Cij134(1,5)-Cij234(1,5)
       S3612222(3)=Cij134(3,5)-Cij234(3,5)
       S3613111(1)=Cij134(3,5)-Cij234(2,1)
       S3613111(2)=Cij134(3,5)+Cij234(3,2)
       S3613111(3)=Cij134(4,5)+Cij234(2,2)
       S3613211(1)=Cij134(3,5)+Cij234(3,2)
       S3613211(2)=Cij134(3,5)-Cij234(3,3)
       S3613211(3)=Cij134(4,5)-Cij234(4,3)
       S3613221(1)=Cij134(3,5)-Cij234(3,3)
       S3613221(2)=Cij134(3,5)+Cij234(3,4)
       S3613221(3)=Cij134(4,5)+Cij234(4,4)
       S3613222(1)=Cij134(3,5)+Cij234(3,4)
       S3613222(2)=Cij134(3,5)-Cij234(3,5)
       S3613222(3)=Cij134(4,5)-Cij234(4,5)
       S3613311(1)=Cij134(4,5)+Cij234(2,2)
       S3613311(2)=Cij134(4,5)-Cij234(4,3)
       S3613311(3)=Cij134(5,5)-Cij234(2,3)
       S3613321(1)=Cij134(4,5)-Cij234(4,3)
       S3613321(2)=Cij134(4,5)+Cij234(4,4)
       S3613321(3)=Cij134(5,5)+Cij234(5,4)
       S3613322(1)=Cij134(4,5)+Cij234(4,4)
       S3613322(2)=Cij134(4,5)-Cij234(4,5)
       S3613322(3)=Cij134(5,5)-Cij234(5,5)
       S3613331(1)=Cij134(5,5)-Cij234(2,3)
       S3613331(2)=Cij134(5,5)+Cij234(5,4)
       S3613331(3)=Cij134(6,5)+Cij234(2,4)
       S3613332(1)=Cij134(5,5)+Cij234(5,4)
       S3613332(2)=Cij134(5,5)-Cij234(5,5)
       S3613332(3)=Cij134(6,5)-Cij234(6,5)
       S3613333(1)=Cij134(6,5)+Cij234(2,4)
       S3613333(2)=Cij134(6,5)-Cij234(6,5)
       S3613333(3)=Cij134(2,5)-Cij234(2,5)
       S3621111(1)=Cij124(1,5)-Cij134(1,5)
       S3621111(2)=Cij124(3,5)-Cij134(1,5)
       S3621111(3)=Cij124(3,5)-Cij134(3,5)
       S3622111(1)=Cij124(3,5)-Cij134(1,5)
       S3622111(2)=Cij124(4,5)-Cij134(1,5)
       S3622111(3)=Cij124(4,5)-Cij134(3,5)
       S3622211(1)=Cij124(4,5)-Cij134(1,5)
       S3622211(2)=Cij124(5,5)-Cij134(1,5)
       S3622211(3)=Cij124(5,5)-Cij134(3,5)
       S3622221(1)=Cij124(5,5)-Cij134(1,5)
       S3622221(2)=Cij124(6,5)-Cij134(1,5)
       S3622221(3)=Cij124(6,5)-Cij134(3,5)
       S3622222(1)=Cij124(6,5)-Cij134(1,5)
       S3622222(2)=Cij124(2,5)-Cij134(1,5)
       S3622222(3)=Cij124(2,5)-Cij134(3,5)
       S3623111(1)=Cij124(3,5)-Cij134(3,5)
       S3623111(2)=Cij124(4,5)-Cij134(3,5)
       S3623111(3)=Cij124(4,5)-Cij134(4,5)
       S3623211(1)=Cij124(4,5)-Cij134(3,5)
       S3623211(2)=Cij124(5,5)-Cij134(3,5)
       S3623211(3)=Cij124(5,5)-Cij134(4,5)
       S3623221(1)=Cij124(5,5)-Cij134(3,5)
       S3623221(2)=Cij124(6,5)-Cij134(3,5)
       S3623221(3)=Cij124(6,5)-Cij134(4,5)
       S3623222(1)=Cij124(6,5)-Cij134(3,5)
       S3623222(2)=Cij124(2,5)-Cij134(3,5)
       S3623222(3)=Cij124(2,5)-Cij134(4,5)
       S3623311(1)=Cij124(4,5)-Cij134(4,5)
       S3623311(2)=Cij124(5,5)-Cij134(4,5)
       S3623311(3)=Cij124(5,5)-Cij134(5,5)
       S3623321(1)=Cij124(5,5)-Cij134(4,5)
       S3623321(2)=Cij124(6,5)-Cij134(4,5)
       S3623321(3)=Cij124(6,5)-Cij134(5,5)
       S3623322(1)=Cij124(6,5)-Cij134(4,5)
       S3623322(2)=Cij124(2,5)-Cij134(4,5)
       S3623322(3)=Cij124(2,5)-Cij134(5,5)
       S3623331(1)=Cij124(5,5)-Cij134(5,5)
       S3623331(2)=Cij124(6,5)-Cij134(5,5)
       S3623331(3)=Cij124(6,5)-Cij134(6,5)
       S3623332(1)=Cij124(6,5)-Cij134(5,5)
       S3623332(2)=Cij124(2,5)-Cij134(5,5)
       S3623332(3)=Cij124(2,5)-Cij134(6,5)
       S3623333(1)=Cij124(6,5)-Cij134(6,5)
       S3623333(2)=Cij124(2,5)-Cij134(6,5)
       S3623333(3)=Cij124(2,5)-Cij134(2,5)
       S3631111(1)=Cij123(1,5)-Cij124(1,5)
       S3631111(2)=Cij123(3,5)-Cij124(3,5)
       S3631111(3)=-Cij124(3,5)
       S3632111(1)=Cij123(3,5)-Cij124(3,5)
       S3632111(2)=Cij123(4,5)-Cij124(4,5)
       S3632111(3)=-Cij124(4,5)
       S3632211(1)=Cij123(4,5)-Cij124(4,5)
       S3632211(2)=Cij123(5,5)-Cij124(5,5)
       S3632211(3)=-Cij124(5,5)
       S3632221(1)=Cij123(5,5)-Cij124(5,5)
       S3632221(2)=Cij123(6,5)-Cij124(6,5)
       S3632221(3)=-Cij124(6,5)
       S3632222(1)=Cij123(6,5)-Cij124(6,5)
       S3632222(2)=Cij123(2,5)-Cij124(2,5)
       S3632222(3)=-Cij124(2,5)
       S3633111(1)=-Cij124(3,5)
       S3633111(2)=-Cij124(4,5)
       S3633111(3)=-Cij124(4,5)
       S3633211(1)=-Cij124(4,5)
       S3633211(2)=-Cij124(5,5)
       S3633211(3)=-Cij124(5,5)
       S3633221(1)=-Cij124(5,5)
       S3633221(2)=-Cij124(6,5)
       S3633221(3)=-Cij124(6,5)
       S3633222(1)=-Cij124(6,5)
       S3633222(2)=-Cij124(2,5)
       S3633222(3)=-Cij124(2,5)
       S3633311(1)=-Cij124(4,5)
       S3633311(2)=-Cij124(5,5)
       S3633311(3)=-Cij124(5,5)
       S3633321(1)=-Cij124(5,5)
       S3633321(2)=-Cij124(6,5)
       S3633321(3)=-Cij124(6,5)
       S3633322(1)=-Cij124(6,5)
       S3633322(2)=-Cij124(2,5)
       S3633322(3)=-Cij124(2,5)
       S3633331(1)=-Cij124(5,5)
       S3633331(2)=-Cij124(6,5)
       S3633331(3)=-Cij124(6,5)
       S3633332(1)=-Cij124(6,5)
       S3633332(2)=-Cij124(2,5)
       S3633332(3)=-Cij124(2,5)
       S3633333(1)=-Cij124(6,5)
       S3633333(2)=-Cij124(2,5)
       S3633333(3)=-Cij124(2,5)
       aux0041(1,1,1)=-(F(1)*S351111(1))-F(2)*S352111(1)-F(3)*S353111(1)
     &  +S3611111(k)*Z(1,l)+S3621111(k)*Z(2,l)+S3631111(k)*Z(3,l)+(S3004
     &  111(1)-S3611111(1)-S3622111(1)-S3633111(1))*Z(k,l)-8*S3h004111(1
     &  )*ZZ(k,1,l,1)-8*S3h004211(1)*ZZ(k,1,l,2)-8*S3h004311(1)*ZZ(k,1,l
     &  ,3)
       aux0042(1,1,1)=-(F(1)*S351211(1))-F(2)*S352211(1)-F(3)*S353211(1)
     &  +S3612111(k)*Z(1,l)+S3622111(k)*Z(2,l)+S3632111(k)*Z(3,l)+(S3004
     &  211(1)-S3612111(1)-S3622211(1)-S3633211(1))*Z(k,l)-6*S3h004121(1
     &  )*ZZ(k,1,l,1)-6*S3h004221(1)*ZZ(k,1,l,2)-6*S3h004321(1)*ZZ(k,1,l
     &  ,3)-2*S3h004111(1)*ZZ(k,2,l,1)-2*S3h004211(1)*ZZ(k,2,l,2)-2*S3h0
     &  04311(1)*ZZ(k,2,l,3)
       aux0042(2,1,1)=-(F(1)*S351221(1))-F(2)*S352221(1)-F(3)*S353221(1)
     &  +S3612211(k)*Z(1,l)+S3622211(k)*Z(2,l)+S3632211(k)*Z(3,l)+(S3004
     &  221(1)-S3612211(1)-S3622221(1)-S3633221(1))*Z(k,l)-4*S3h004122(1
     &  )*ZZ(k,1,l,1)-4*S3h004222(1)*ZZ(k,1,l,2)-4*S3h004322(1)*ZZ(k,1,l
     &  ,3)-4*S3h004121(1)*ZZ(k,2,l,1)-4*S3h004221(1)*ZZ(k,2,l,2)-4*S3h0
     &  04321(1)*ZZ(k,2,l,3)
       aux0042(2,2,1)=-(F(1)*S351222(1))-F(2)*S352222(1)-F(3)*S353222(1)
     &  +S3612221(k)*Z(1,l)+S3622221(k)*Z(2,l)+S3632221(k)*Z(3,l)+(S3004
     &  222(1)-S3612221(1)-S3622222(1)-S3633222(1))*Z(k,l)-2*S3h004122(2
     &  )*ZZ(k,1,l,1)-2*S3h004222(2)*ZZ(k,1,l,2)-2*S3h004322(2)*ZZ(k,1,l
     &  ,3)-6*S3h004122(1)*ZZ(k,2,l,1)-6*S3h004222(1)*ZZ(k,2,l,2)-6*S3h0
     &  04322(1)*ZZ(k,2,l,3)
       aux0042(2,2,2)=-(F(1)*S351222(2))-F(2)*S352222(2)-F(3)*S353222(2)
     &  +S3612222(k)*Z(1,l)+S3622222(k)*Z(2,l)+S3632222(k)*Z(3,l)+(S3004
     &  222(2)-S3612222(1)-S3622222(2)-S3633222(2))*Z(k,l)-8*S3h004122(2
     &  )*ZZ(k,2,l,1)-8*S3h004222(2)*ZZ(k,2,l,2)-8*S3h004322(2)*ZZ(k,2,l
     &  ,3)
       aux0043(1,1,1)=-(F(1)*S351311(1))-F(2)*S352311(1)-F(3)*S353311(1)
     &  +S3613111(k)*Z(1,l)+S3623111(k)*Z(2,l)+S3633111(k)*Z(3,l)+(S3004
     &  311(1)-S3613111(1)-S3623211(1)-S3633311(1))*Z(k,l)-6*S3h004131(1
     &  )*ZZ(k,1,l,1)-6*S3h004231(1)*ZZ(k,1,l,2)-6*S3h004331(1)*ZZ(k,1,l
     &  ,3)-2*S3h004111(1)*ZZ(k,3,l,1)-2*S3h004211(1)*ZZ(k,3,l,2)-2*S3h0
     &  04311(1)*ZZ(k,3,l,3)
       aux0043(2,1,1)=-(F(1)*S351321(1))-F(2)*S352321(1)-F(3)*S353321(1)
     &  +S3613211(k)*Z(1,l)+S3623211(k)*Z(2,l)+S3633211(k)*Z(3,l)+(S3004
     &  321(1)-S3613211(1)-S3623221(1)-S3633321(1))*Z(k,l)-4*S3h004132(1
     &  )*ZZ(k,1,l,1)-4*S3h004232(1)*ZZ(k,1,l,2)-4*S3h004332(1)*ZZ(k,1,l
     &  ,3)-2*S3h004131(1)*ZZ(k,2,l,1)-2*S3h004231(1)*ZZ(k,2,l,2)-2*S3h0
     &  04331(1)*ZZ(k,2,l,3)-2*S3h004121(1)*ZZ(k,3,l,1)-2*S3h004221(1)*Z
     &  Z(k,3,l,2)-2*S3h004321(1)*ZZ(k,3,l,3)
       aux0043(2,2,1)=-(F(1)*S351322(1))-F(2)*S352322(1)-F(3)*S353322(1)
     &  +S3613221(k)*Z(1,l)+S3623221(k)*Z(2,l)+S3633221(k)*Z(3,l)+(S3004
     &  322(1)-S3613221(1)-S3623222(1)-S3633322(1))*Z(k,l)-2*S3h004132(2
     &  )*ZZ(k,1,l,1)-2*S3h004232(2)*ZZ(k,1,l,2)-2*S3h004332(2)*ZZ(k,1,l
     &  ,3)-4*S3h004132(1)*ZZ(k,2,l,1)-4*S3h004232(1)*ZZ(k,2,l,2)-4*S3h0
     &  04332(1)*ZZ(k,2,l,3)-2*S3h004122(1)*ZZ(k,3,l,1)-2*S3h004222(1)*Z
     &  Z(k,3,l,2)-2*S3h004322(1)*ZZ(k,3,l,3)
       aux0043(2,2,2)=-(F(1)*S351322(2))-F(2)*S352322(2)-F(3)*S353322(2)
     &  +S3613222(k)*Z(1,l)+S3623222(k)*Z(2,l)+S3633222(k)*Z(3,l)+(S3004
     &  322(2)-S3613222(1)-S3623222(2)-S3633322(2))*Z(k,l)-6*S3h004132(2
     &  )*ZZ(k,2,l,1)-6*S3h004232(2)*ZZ(k,2,l,2)-6*S3h004332(2)*ZZ(k,2,l
     &  ,3)-2*S3h004122(2)*ZZ(k,3,l,1)-2*S3h004222(2)*ZZ(k,3,l,2)-2*S3h0
     &  04322(2)*ZZ(k,3,l,3)
       aux0043(3,1,1)=-(F(1)*S351331(1))-F(2)*S352331(1)-F(3)*S353331(1)
     &  +S3613311(k)*Z(1,l)+S3623311(k)*Z(2,l)+S3633311(k)*Z(3,l)+(S3004
     &  331(1)-S3613311(1)-S3623321(1)-S3633331(1))*Z(k,l)-4*S3h004133(1
     &  )*ZZ(k,1,l,1)-4*S3h004233(1)*ZZ(k,1,l,2)-4*S3h004333(1)*ZZ(k,1,l
     &  ,3)-4*S3h004131(1)*ZZ(k,3,l,1)-4*S3h004231(1)*ZZ(k,3,l,2)-4*S3h0
     &  04331(1)*ZZ(k,3,l,3)
       aux0043(3,2,1)=-(F(1)*S351332(1))-F(2)*S352332(1)-F(3)*S353332(1)
     &  +S3613321(k)*Z(1,l)+S3623321(k)*Z(2,l)+S3633321(k)*Z(3,l)+(S3004
     &  332(1)-S3613321(1)-S3623322(1)-S3633332(1))*Z(k,l)-2*S3h004133(2
     &  )*ZZ(k,1,l,1)-2*S3h004233(2)*ZZ(k,1,l,2)-2*S3h004333(2)*ZZ(k,1,l
     &  ,3)-2*S3h004133(1)*ZZ(k,2,l,1)-2*S3h004233(1)*ZZ(k,2,l,2)-2*S3h0
     &  04333(1)*ZZ(k,2,l,3)-4*S3h004132(1)*ZZ(k,3,l,1)-4*S3h004232(1)*Z
     &  Z(k,3,l,2)-4*S3h004332(1)*ZZ(k,3,l,3)
       aux0043(3,2,2)=-(F(1)*S351332(2))-F(2)*S352332(2)-F(3)*S353332(2)
     &  +S3613322(k)*Z(1,l)+S3623322(k)*Z(2,l)+S3633322(k)*Z(3,l)+(S3004
     &  332(2)-S3613322(1)-S3623322(2)-S3633332(2))*Z(k,l)-4*S3h004133(2
     &  )*ZZ(k,2,l,1)-4*S3h004233(2)*ZZ(k,2,l,2)-4*S3h004333(2)*ZZ(k,2,l
     &  ,3)-4*S3h004132(2)*ZZ(k,3,l,1)-4*S3h004232(2)*ZZ(k,3,l,2)-4*S3h0
     &  04332(2)*ZZ(k,3,l,3)
       aux0043(3,3,1)=-(F(1)*S351333(1))-F(2)*S352333(1)-F(3)*S353333(1)
     &  +S3613331(k)*Z(1,l)+S3623331(k)*Z(2,l)+S3633331(k)*Z(3,l)+(S3004
     &  333(1)-S3613331(1)-S3623332(1)-S3633333(1))*Z(k,l)-2*S3h004133(3
     &  )*ZZ(k,1,l,1)-2*S3h004233(3)*ZZ(k,1,l,2)-2*S3h004333(3)*ZZ(k,1,l
     &  ,3)-6*S3h004133(1)*ZZ(k,3,l,1)-6*S3h004233(1)*ZZ(k,3,l,2)-6*S3h0
     &  04333(1)*ZZ(k,3,l,3)
       aux0043(3,3,2)=-(F(1)*S351333(2))-F(2)*S352333(2)-F(3)*S353333(2)
     &  +S3613332(k)*Z(1,l)+S3623332(k)*Z(2,l)+S3633332(k)*Z(3,l)+(S3004
     &  333(2)-S3613332(1)-S3623332(2)-S3633333(2))*Z(k,l)-2*S3h004133(3
     &  )*ZZ(k,2,l,1)-2*S3h004233(3)*ZZ(k,2,l,2)-2*S3h004333(3)*ZZ(k,2,l
     &  ,3)-6*S3h004133(2)*ZZ(k,3,l,1)-6*S3h004233(2)*ZZ(k,3,l,2)-6*S3h0
     &  04333(2)*ZZ(k,3,l,3)
       aux0043(3,3,3)=-(F(1)*S351333(3))-F(2)*S352333(3)-F(3)*S353333(3)
     &  +S3613333(k)*Z(1,l)+S3623333(k)*Z(2,l)+S3633333(k)*Z(3,l)+(S3004
     &  333(3)-S3613333(1)-S3623333(2)-S3633333(3))*Z(k,l)-8*S3h004133(3
     &  )*ZZ(k,3,l,1)-8*S3h004233(3)*ZZ(k,3,l,2)-8*S3h004333(3)*ZZ(k,3,l
     &  ,3)
       temp0041(1,1,1)=I20Z*(aux0041(1,1,1)+8*F(4)*temp003(1,1,1)+F(5)*t
     &  emp41(1,1,1)+48*temp00002(1,1)*ZZ(k,1,l,1))
       temp0042(1,1,1)=I20Z*(aux0042(1,1,1)+2*F(6)*temp003(1,1,1)+6*F(4)
     &  *temp003(2,1,1)+F(5)*temp42(1,1,1)+24*temp00002(2,1)*ZZ(k,1,l,1)
     &  +12*temp00002(1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0042(2,1,1)=I20Z*(aux0042(2,1,1)+4*F(6)*temp003(2,1,1)+4*F(4)
     &  *temp003(2,2,1)+F(5)*temp42(2,1,1)+16*temp00002(2,1)*(ZZ(k,1,l,2
     &  )+ZZ(k,2,l,1))+8*(temp00002(2,2)*ZZ(k,1,l,1)+temp00002(1,1)*ZZ(k
     &  ,2,l,2)))
       temp0042(2,2,1)=I20Z*(aux0042(2,2,1)+6*F(6)*temp003(2,2,1)+2*F(4)
     &  *temp003(2,2,2)+F(5)*temp42(2,2,1)+12*temp00002(2,2)*(ZZ(k,1,l,2
     &  )+ZZ(k,2,l,1))+24*temp00002(2,1)*ZZ(k,2,l,2))
       temp0042(2,2,2)=I20Z*(aux0042(2,2,2)+8*F(6)*temp003(2,2,2)+F(5)*t
     &  emp42(2,2,2)+48*temp00002(2,2)*ZZ(k,2,l,2))
       temp0043(1,1,1)=I20Z*(aux0043(1,1,1)+2*F(7)*temp003(1,1,1)+6*F(4)
     &  *temp003(3,1,1)+F(5)*temp43(1,1,1)+24*temp00002(3,1)*ZZ(k,1,l,1)
     &  +12*temp00002(1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)))
       temp0043(2,1,1)=I20Z*(aux0043(2,1,1)+2*F(7)*temp003(2,1,1)+2*F(6)
     &  *temp003(3,1,1)+4*F(4)*temp003(3,2,1)+F(5)*temp43(2,1,1)+8*(temp
     &  00002(3,2)*ZZ(k,1,l,1)+temp00002(3,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
     &  +8*temp00002(2,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+4*temp00002(1,1)*(ZZ
     &  (k,2,l,3)+ZZ(k,3,l,2)))
       temp0043(2,2,1)=I20Z*(aux0043(2,2,1)+2*F(7)*temp003(2,2,1)+4*F(6)
     &  *temp003(3,2,1)+2*F(4)*temp003(3,2,2)+F(5)*temp43(2,2,1)+8*(temp
     &  00002(3,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00002(3,1)*ZZ(k,2,l,2))
     &  +4*temp00002(2,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*temp00002(2,1)*(ZZ
     &  (k,2,l,3)+ZZ(k,3,l,2)))
       temp0043(2,2,2)=I20Z*(aux0043(2,2,2)+2*F(7)*temp003(2,2,2)+6*F(6)
     &  *temp003(3,2,2)+F(5)*temp43(2,2,2)+24*temp00002(3,2)*ZZ(k,2,l,2)
     &  +12*temp00002(2,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp0043(3,1,1)=I20Z*(aux0043(3,1,1)+4*F(7)*temp003(3,1,1)+4*F(4)
     &  *temp003(3,3,1)+F(5)*temp43(3,1,1)+16*temp00002(3,1)*(ZZ(k,1,l,3
     &  )+ZZ(k,3,l,1))+8*(temp00002(3,3)*ZZ(k,1,l,1)+temp00002(1,1)*ZZ(k
     &  ,3,l,3)))
       temp0043(3,2,1)=I20Z*(aux0043(3,2,1)+2*F(6)*temp003(3,3,1)+2*F(4)
     &  *temp003(3,3,2)+F(5)*temp43(3,2,1)+4*(F(7)*temp003(3,2,1)+temp00
     &  002(3,3)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))+8*temp00002(3,2)*(ZZ(k,1,l,3
     &  )+ZZ(k,3,l,1))+8*temp00002(3,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*temp
     &  00002(2,1)*ZZ(k,3,l,3))
       temp0043(3,2,2)=I20Z*(aux0043(3,2,2)+4*F(7)*temp003(3,2,2)+4*F(6)
     &  *temp003(3,3,2)+F(5)*temp43(3,2,2)+16*temp00002(3,2)*(ZZ(k,2,l,3
     &  )+ZZ(k,3,l,2))+8*(temp00002(3,3)*ZZ(k,2,l,2)+temp00002(2,2)*ZZ(k
     &  ,3,l,3)))
       temp0043(3,3,1)=I20Z*(aux0043(3,3,1)+6*F(7)*temp003(3,3,1)+2*F(4)
     &  *temp003(3,3,3)+F(5)*temp43(3,3,1)+12*temp00002(3,3)*(ZZ(k,1,l,3
     &  )+ZZ(k,3,l,1))+24*temp00002(3,1)*ZZ(k,3,l,3))
       temp0043(3,3,2)=I20Z*(aux0043(3,3,2)+6*F(7)*temp003(3,3,2)+2*F(6)
     &  *temp003(3,3,3)+F(5)*temp43(3,3,2)+12*temp00002(3,3)*(ZZ(k,2,l,3
     &  )+ZZ(k,3,l,2))+24*temp00002(3,2)*ZZ(k,3,l,3))
       temp0043(3,3,3)=I20Z*(aux0043(3,3,3)+8*F(7)*temp003(3,3,3)+F(5)*t
     &  emp43(3,3,3)+48*temp00002(3,3)*ZZ(k,3,l,3))
       temp0041(1,1,2)=temp0042(1,1,1)
       temp0041(1,1,3)=temp0043(1,1,1)
       temp0041(1,2,1)=temp0042(1,1,1)
       temp0041(1,2,2)=temp0042(2,1,1)
       temp0041(1,2,3)=temp0043(2,1,1)
       temp0041(1,3,1)=temp0043(1,1,1)
       temp0041(1,3,2)=temp0043(2,1,1)
       temp0041(1,3,3)=temp0043(3,1,1)
       temp0042(1,1,2)=temp0042(2,1,1)
       temp0042(1,1,3)=temp0043(2,1,1)
       temp0042(1,2,1)=temp0042(2,1,1)
       temp0042(1,2,2)=temp0042(2,2,1)
       temp0042(1,2,3)=temp0043(2,2,1)
       temp0042(1,3,1)=temp0043(2,1,1)
       temp0042(1,3,2)=temp0043(2,2,1)
       temp0042(1,3,3)=temp0043(3,2,1)
       temp0042(2,1,2)=temp0042(2,2,1)
       temp0042(2,1,3)=temp0043(2,2,1)
       temp0042(2,2,3)=temp0043(2,2,2)
       temp0042(2,3,1)=temp0043(2,2,1)
       temp0042(2,3,2)=temp0043(2,2,2)
       temp0042(2,3,3)=temp0043(3,2,2)
       temp0043(1,1,2)=temp0043(2,1,1)
       temp0043(1,1,3)=temp0043(3,1,1)
       temp0043(1,2,1)=temp0043(2,1,1)
       temp0043(1,2,2)=temp0043(2,2,1)
       temp0043(1,2,3)=temp0043(3,2,1)
       temp0043(1,3,1)=temp0043(3,1,1)
       temp0043(1,3,2)=temp0043(3,2,1)
       temp0043(1,3,3)=temp0043(3,3,1)
       temp0043(2,1,2)=temp0043(2,2,1)
       temp0043(2,1,3)=temp0043(3,2,1)
       temp0043(2,2,3)=temp0043(3,2,2)
       temp0043(2,3,1)=temp0043(3,2,1)
       temp0043(2,3,2)=temp0043(3,2,2)
       temp0043(2,3,3)=temp0043(3,3,2)
       temp0043(3,1,2)=temp0043(3,2,1)
       temp0043(3,1,3)=temp0043(3,3,1)
       temp0043(3,2,3)=temp0043(3,3,2)
       aux511(1,1,1)=-(S3611111(1)*Z(jj,1))-S3621111(1)*Z(jj,2)-S3631111
     &  (1)*Z(jj,3)
       aux521(1,1,1)=-(S3612111(1)*Z(jj,1))-S3622111(1)*Z(jj,2)-S3632111
     &  (1)*Z(jj,3)
       aux522(1,1,1)=-(S3612211(1)*Z(jj,1))-S3622211(1)*Z(jj,2)-S3632211
     &  (1)*Z(jj,3)
       aux522(2,1,1)=-(S3612221(1)*Z(jj,1))-S3622221(1)*Z(jj,2)-S3632221
     &  (1)*Z(jj,3)
       aux522(2,2,1)=-(S3612222(1)*Z(jj,1))-S3622222(1)*Z(jj,2)-S3632222
     &  (1)*Z(jj,3)
       aux522(2,2,2)=-(S3612222(2)*Z(jj,1))-S3622222(2)*Z(jj,2)-S3632222
     &  (2)*Z(jj,3)
       aux531(1,1,1)=-(S3613111(1)*Z(jj,1))-S3623111(1)*Z(jj,2)-S3633111
     &  (1)*Z(jj,3)
       aux532(1,1,1)=-(S3613211(1)*Z(jj,1))-S3623211(1)*Z(jj,2)-S3633211
     &  (1)*Z(jj,3)
       aux532(2,1,1)=-(S3613221(1)*Z(jj,1))-S3623221(1)*Z(jj,2)-S3633221
     &  (1)*Z(jj,3)
       aux532(2,2,1)=-(S3613222(1)*Z(jj,1))-S3623222(1)*Z(jj,2)-S3633222
     &  (1)*Z(jj,3)
       aux532(2,2,2)=-(S3613222(2)*Z(jj,1))-S3623222(2)*Z(jj,2)-S3633222
     &  (2)*Z(jj,3)
       aux533(1,1,1)=-(S3613311(1)*Z(jj,1))-S3623311(1)*Z(jj,2)-S3633311
     &  (1)*Z(jj,3)
       aux533(2,1,1)=-(S3613321(1)*Z(jj,1))-S3623321(1)*Z(jj,2)-S3633321
     &  (1)*Z(jj,3)
       aux533(2,2,1)=-(S3613322(1)*Z(jj,1))-S3623322(1)*Z(jj,2)-S3633322
     &  (1)*Z(jj,3)
       aux533(2,2,2)=-(S3613322(2)*Z(jj,1))-S3623322(2)*Z(jj,2)-S3633322
     &  (2)*Z(jj,3)
       aux533(3,1,1)=-(S3613331(1)*Z(jj,1))-S3623331(1)*Z(jj,2)-S3633331
     &  (1)*Z(jj,3)
       aux533(3,2,1)=-(S3613332(1)*Z(jj,1))-S3623332(1)*Z(jj,2)-S3633332
     &  (1)*Z(jj,3)
       aux533(3,2,2)=-(S3613332(2)*Z(jj,1))-S3623332(2)*Z(jj,2)-S3633332
     &  (2)*Z(jj,3)
       aux533(3,3,1)=-(S3613333(1)*Z(jj,1))-S3623333(1)*Z(jj,2)-S3633333
     &  (1)*Z(jj,3)
       aux533(3,3,2)=-(S3613333(2)*Z(jj,1))-S3623333(2)*Z(jj,2)-S3633333
     &  (2)*Z(jj,3)
       aux533(3,3,3)=-(S3613333(3)*Z(jj,1))-S3623333(3)*Z(jj,2)-S3633333
     &  (3)*Z(jj,3)
       temp511(1,1,1)=IX*(aux511(1,1,1)+10*temp0041(1,1,1)*Z(jj,1))
       temp521(1,1,1)=IX*(aux521(1,1,1)+8*temp0042(1,1,1)*Z(jj,1)+2*temp
     &  0041(1,1,1)*Z(jj,2))
       temp522(1,1,1)=IX*(aux522(1,1,1)+6*temp0042(2,1,1)*Z(jj,1)+4*temp
     &  0042(1,1,1)*Z(jj,2))
       temp522(2,1,1)=IX*(aux522(2,1,1)+4*temp0042(2,2,1)*Z(jj,1)+6*temp
     &  0042(2,1,1)*Z(jj,2))
       temp522(2,2,1)=IX*(aux522(2,2,1)+2*temp0042(2,2,2)*Z(jj,1)+8*temp
     &  0042(2,2,1)*Z(jj,2))
       temp522(2,2,2)=IX*(aux522(2,2,2)+10*temp0042(2,2,2)*Z(jj,2))
       temp531(1,1,1)=IX*(aux531(1,1,1)+8*temp0043(1,1,1)*Z(jj,1)+2*temp
     &  0041(1,1,1)*Z(jj,3))
       temp532(1,1,1)=IX*(aux532(1,1,1)+6*temp0043(2,1,1)*Z(jj,1)+2*(tem
     &  p0043(1,1,1)*Z(jj,2)+temp0042(1,1,1)*Z(jj,3)))
       temp532(2,1,1)=IX*(aux532(2,1,1)+4*(temp0043(2,2,1)*Z(jj,1)+temp0
     &  043(2,1,1)*Z(jj,2))+2*temp0042(2,1,1)*Z(jj,3))
       temp532(2,2,1)=IX*(aux532(2,2,1)+6*temp0043(2,2,1)*Z(jj,2)+2*(tem
     &  p0043(2,2,2)*Z(jj,1)+temp0042(2,2,1)*Z(jj,3)))
       temp532(2,2,2)=IX*(aux532(2,2,2)+8*temp0043(2,2,2)*Z(jj,2)+2*temp
     &  0042(2,2,2)*Z(jj,3))
       temp533(1,1,1)=IX*(aux533(1,1,1)+6*temp0043(3,1,1)*Z(jj,1)+4*temp
     &  0043(1,1,1)*Z(jj,3))
       temp533(2,1,1)=IX*(aux533(2,1,1)+2*temp0043(3,1,1)*Z(jj,2)+4*(tem
     &  p0043(3,2,1)*Z(jj,1)+temp0043(2,1,1)*Z(jj,3)))
       temp533(2,2,1)=IX*(aux533(2,2,1)+2*temp0043(3,2,2)*Z(jj,1)+4*(tem
     &  p0043(3,2,1)*Z(jj,2)+temp0043(2,2,1)*Z(jj,3)))
       temp533(2,2,2)=IX*(aux533(2,2,2)+6*temp0043(3,2,2)*Z(jj,2)+4*temp
     &  0043(2,2,2)*Z(jj,3))
       temp533(3,1,1)=IX*(aux533(3,1,1)+4*temp0043(3,3,1)*Z(jj,1)+6*temp
     &  0043(3,1,1)*Z(jj,3))
       temp533(3,2,1)=IX*(aux533(3,2,1)+2*(temp0043(3,3,2)*Z(jj,1)+temp0
     &  043(3,3,1)*Z(jj,2))+6*temp0043(3,2,1)*Z(jj,3))
       temp533(3,2,2)=IX*(aux533(3,2,2)+4*temp0043(3,3,2)*Z(jj,2)+6*temp
     &  0043(3,2,2)*Z(jj,3))
       temp533(3,3,1)=IX*(aux533(3,3,1)+2*temp0043(3,3,3)*Z(jj,1)+8*temp
     &  0043(3,3,1)*Z(jj,3))
       temp533(3,3,2)=IX*(aux533(3,3,2)+2*temp0043(3,3,3)*Z(jj,2)+8*temp
     &  0043(3,3,2)*Z(jj,3))
       temp533(3,3,3)=IX*(aux533(3,3,3)+10*temp0043(3,3,3)*Z(jj,3))
       temp511(1,1,2)=temp521(1,1,1)
       temp511(1,1,3)=temp531(1,1,1)
       temp511(1,2,1)=temp521(1,1,1)
       temp511(1,2,2)=temp522(1,1,1)
       temp511(1,2,3)=temp532(1,1,1)
       temp511(1,3,1)=temp531(1,1,1)
       temp511(1,3,2)=temp532(1,1,1)
       temp511(1,3,3)=temp533(1,1,1)
       temp521(1,1,2)=temp522(1,1,1)
       temp521(1,1,3)=temp532(1,1,1)
       temp521(1,2,1)=temp522(1,1,1)
       temp521(1,2,2)=temp522(2,1,1)
       temp521(1,2,3)=temp532(2,1,1)
       temp521(1,3,1)=temp532(1,1,1)
       temp521(1,3,2)=temp532(2,1,1)
       temp521(1,3,3)=temp533(2,1,1)
       temp522(1,1,2)=temp522(2,1,1)
       temp522(1,1,3)=temp532(2,1,1)
       temp522(1,2,1)=temp522(2,1,1)
       temp522(1,2,2)=temp522(2,2,1)
       temp522(1,2,3)=temp532(2,2,1)
       temp522(1,3,1)=temp532(2,1,1)
       temp522(1,3,2)=temp532(2,2,1)
       temp522(1,3,3)=temp533(2,2,1)
       temp522(2,1,2)=temp522(2,2,1)
       temp522(2,1,3)=temp532(2,2,1)
       temp522(2,2,3)=temp532(2,2,2)
       temp522(2,3,1)=temp532(2,2,1)
       temp522(2,3,2)=temp532(2,2,2)
       temp522(2,3,3)=temp533(2,2,2)
       temp531(1,1,2)=temp532(1,1,1)
       temp531(1,1,3)=temp533(1,1,1)
       temp531(1,2,1)=temp532(1,1,1)
       temp531(1,2,2)=temp532(2,1,1)
       temp531(1,2,3)=temp533(2,1,1)
       temp531(1,3,1)=temp533(1,1,1)
       temp531(1,3,2)=temp533(2,1,1)
       temp531(1,3,3)=temp533(3,1,1)
       temp532(1,1,2)=temp532(2,1,1)
       temp532(1,1,3)=temp533(2,1,1)
       temp532(1,2,1)=temp532(2,1,1)
       temp532(1,2,2)=temp532(2,2,1)
       temp532(1,2,3)=temp533(2,2,1)
       temp532(1,3,1)=temp533(2,1,1)
       temp532(1,3,2)=temp533(2,2,1)
       temp532(1,3,3)=temp533(3,2,1)
       temp532(2,1,2)=temp532(2,2,1)
       temp532(2,1,3)=temp533(2,2,1)
       temp532(2,2,3)=temp533(2,2,2)
       temp532(2,3,1)=temp533(2,2,1)
       temp532(2,3,2)=temp533(2,2,2)
       temp532(2,3,3)=temp533(3,2,2)
       temp533(1,1,2)=temp533(2,1,1)
       temp533(1,1,3)=temp533(3,1,1)
       temp533(1,2,1)=temp533(2,1,1)
       temp533(1,2,2)=temp533(2,2,1)
       temp533(1,2,3)=temp533(3,2,1)
       temp533(1,3,1)=temp533(3,1,1)
       temp533(1,3,2)=temp533(3,2,1)
       temp533(1,3,3)=temp533(3,3,1)
       temp533(2,1,2)=temp533(2,2,1)
       temp533(2,1,3)=temp533(3,2,1)
       temp533(2,2,3)=temp533(3,2,2)
       temp533(2,3,1)=temp533(3,2,1)
       temp533(2,3,2)=temp533(3,2,2)
       temp533(2,3,3)=temp533(3,3,2)
       temp533(3,1,2)=temp533(3,2,1)
       temp533(3,1,3)=temp533(3,3,1)
       temp533(3,2,3)=temp533(3,3,2)
c                Step2
       temp00001(1)=I12Z*(aux00001(1)+2*tempD40000*F(4)+F(5)*temp001(1)-
     &  det4*temp003(1,k,l))
       temp00001(2)=I12Z*(aux00001(2)+2*tempD40000*F(6)+F(5)*temp001(2)-
     &  det4*temp003(2,k,l))
       temp00001(3)=I12Z*(aux00001(3)+2*tempD40000*F(7)+F(5)*temp001(3)-
     &  det4*temp003(3,k,l))
       temp003(1,1,1)=I16Z*(aux003(1,1,1)+6*F(4)*temp002(1,1)+F(5)*temp3
     &  (1,1,1)-det4*temp511(1,k,l)+24*temp00001(1)*ZZ(k,1,l,1))
       temp003(2,1,1)=I16Z*(aux003(2,1,1)+2*F(6)*temp002(1,1)+4*F(4)*tem
     &  p002(2,1)+F(5)*temp3(2,1,1)-det4*temp521(1,k,l)+8*(temp00001(2)*
     &  ZZ(k,1,l,1)+temp00001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))))
       temp003(2,2,1)=I16Z*(aux003(2,2,1)+4*F(6)*temp002(2,1)+2*F(4)*tem
     &  p002(2,2)+F(5)*temp3(2,2,1)-det4*temp522(1,k,l)+8*(temp00001(2)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00001(1)*ZZ(k,2,l,2)))
       temp003(2,2,2)=I16Z*(aux003(2,2,2)+6*F(6)*temp002(2,2)+F(5)*temp3
     &  (2,2,2)-det4*temp522(2,k,l)+24*temp00001(2)*ZZ(k,2,l,2))
       temp003(3,1,1)=I16Z*(aux003(3,1,1)+2*F(7)*temp002(1,1)+4*F(4)*tem
     &  p002(3,1)+F(5)*temp3(3,1,1)-det4*temp531(1,k,l)+8*(temp00001(3)*
     &  ZZ(k,1,l,1)+temp00001(1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))))
       temp003(3,2,1)=I16Z*(aux003(3,2,1)+2*F(7)*temp002(2,1)+2*F(6)*tem
     &  p002(3,1)+2*F(4)*temp002(3,2)+F(5)*temp3(3,2,1)-det4*temp532(1,k
     &  ,l)+4*(temp00001(3)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00001(2)*(ZZ(k
     &  ,1,l,3)+ZZ(k,3,l,1)))+4*temp00001(1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp003(3,2,2)=I16Z*(aux003(3,2,2)+2*F(7)*temp002(2,2)+4*F(6)*tem
     &  p002(3,2)+F(5)*temp3(3,2,2)-det4*temp532(2,k,l)+8*(temp00001(3)*
     &  ZZ(k,2,l,2)+temp00001(2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))))
       temp003(3,3,1)=I16Z*(aux003(3,3,1)+4*F(7)*temp002(3,1)+2*F(4)*tem
     &  p002(3,3)+F(5)*temp3(3,3,1)-det4*temp533(1,k,l)+8*(temp00001(3)*
     &  (ZZ(k,1,l,3)+ZZ(k,3,l,1))+temp00001(1)*ZZ(k,3,l,3)))
       temp003(3,3,2)=I16Z*(aux003(3,3,2)+4*F(7)*temp002(3,2)+2*F(6)*tem
     &  p002(3,3)+F(5)*temp3(3,3,2)-det4*temp533(2,k,l)+8*(temp00001(3)*
     &  (ZZ(k,2,l,3)+ZZ(k,3,l,2))+temp00001(2)*ZZ(k,3,l,3)))
       temp003(3,3,3)=I16Z*(aux003(3,3,3)+6*F(7)*temp002(3,3)+F(5)*temp3
     &  (3,3,3)-det4*temp533(3,k,l)+24*temp00001(3)*ZZ(k,3,l,3))
       temp003(1,1,2)=temp003(2,1,1)
       temp003(1,1,3)=temp003(3,1,1)
       temp003(1,2,1)=temp003(2,1,1)
       temp003(1,2,2)=temp003(2,2,1)
       temp003(1,2,3)=temp003(3,2,1)
       temp003(1,3,1)=temp003(3,1,1)
       temp003(1,3,2)=temp003(3,2,1)
       temp003(1,3,3)=temp003(3,3,1)
       temp003(2,1,2)=temp003(2,2,1)
       temp003(2,1,3)=temp003(3,2,1)
       temp003(2,2,3)=temp003(3,2,2)
       temp003(2,3,1)=temp003(3,2,1)
       temp003(2,3,2)=temp003(3,2,2)
       temp003(2,3,3)=temp003(3,3,2)
       temp003(3,1,2)=temp003(3,2,1)
       temp003(3,1,3)=temp003(3,3,1)
       temp003(3,2,3)=temp003(3,3,2)
       temp41(1,1,1)=IX*(aux41(1,1,1)+det4*temp511(1,1,jj)+8*temp003(1,1
     &  ,1)*Z(jj,1))
       temp42(1,1,1)=IX*(aux42(1,1,1)+det4*temp521(1,1,jj)+6*temp003(2,1
     &  ,1)*Z(jj,1)+2*temp003(1,1,1)*Z(jj,2))
       temp42(2,1,1)=IX*(aux42(2,1,1)+det4*temp522(1,1,jj)+4*(temp003(2,
     &  2,1)*Z(jj,1)+temp003(2,1,1)*Z(jj,2)))
       temp42(2,2,1)=IX*(aux42(2,2,1)+det4*temp522(2,1,jj)+2*temp003(2,2
     &  ,2)*Z(jj,1)+6*temp003(2,2,1)*Z(jj,2))
       temp42(2,2,2)=IX*(aux42(2,2,2)+det4*temp522(2,2,jj)+8*temp003(2,2
     &  ,2)*Z(jj,2))
       temp43(1,1,1)=IX*(aux43(1,1,1)+det4*temp531(1,1,jj)+6*temp003(3,1
     &  ,1)*Z(jj,1)+2*temp003(1,1,1)*Z(jj,3))
       temp43(2,1,1)=IX*(aux43(2,1,1)+det4*temp532(1,1,jj)+4*temp003(3,2
     &  ,1)*Z(jj,1)+2*(temp003(3,1,1)*Z(jj,2)+temp003(2,1,1)*Z(jj,3)))
       temp43(2,2,1)=IX*(aux43(2,2,1)+det4*temp532(2,1,jj)+4*temp003(3,2
     &  ,1)*Z(jj,2)+2*(temp003(3,2,2)*Z(jj,1)+temp003(2,2,1)*Z(jj,3)))
       temp43(2,2,2)=IX*(aux43(2,2,2)+det4*temp532(2,2,jj)+6*temp003(3,2
     &  ,2)*Z(jj,2)+2*temp003(2,2,2)*Z(jj,3))
       temp43(3,1,1)=IX*(aux43(3,1,1)+det4*temp533(1,1,jj)+4*(temp003(3,
     &  3,1)*Z(jj,1)+temp003(3,1,1)*Z(jj,3)))
       temp43(3,2,1)=IX*(aux43(3,2,1)+det4*temp533(2,1,jj)+2*(temp003(3,
     &  3,2)*Z(jj,1)+temp003(3,3,1)*Z(jj,2))+4*temp003(3,2,1)*Z(jj,3))
       temp43(3,2,2)=IX*(aux43(3,2,2)+det4*temp533(2,2,jj)+4*(temp003(3,
     &  3,2)*Z(jj,2)+temp003(3,2,2)*Z(jj,3)))
       temp43(3,3,1)=IX*(aux43(3,3,1)+det4*temp533(3,1,jj)+2*temp003(3,3
     &  ,3)*Z(jj,1)+6*temp003(3,3,1)*Z(jj,3))
       temp43(3,3,2)=IX*(aux43(3,3,2)+det4*temp533(3,2,jj)+2*temp003(3,3
     &  ,3)*Z(jj,2)+6*temp003(3,3,2)*Z(jj,3))
       temp43(3,3,3)=IX*(aux43(3,3,3)+det4*temp533(3,3,jj)+8*temp003(3,3
     &  ,3)*Z(jj,3))
       temp41(1,1,2)=temp42(1,1,1)
       temp41(1,1,3)=temp43(1,1,1)
       temp41(1,2,1)=temp42(1,1,1)
       temp41(1,2,2)=temp42(2,1,1)
       temp41(1,2,3)=temp43(2,1,1)
       temp41(1,3,1)=temp43(1,1,1)
       temp41(1,3,2)=temp43(2,1,1)
       temp41(1,3,3)=temp43(3,1,1)
       temp42(1,1,2)=temp42(2,1,1)
       temp42(1,1,3)=temp43(2,1,1)
       temp42(1,2,1)=temp42(2,1,1)
       temp42(1,2,2)=temp42(2,2,1)
       temp42(1,2,3)=temp43(2,2,1)
       temp42(1,3,1)=temp43(2,1,1)
       temp42(1,3,2)=temp43(2,2,1)
       temp42(1,3,3)=temp43(3,2,1)
       temp42(2,1,2)=temp42(2,2,1)
       temp42(2,1,3)=temp43(2,2,1)
       temp42(2,2,3)=temp43(2,2,2)
       temp42(2,3,1)=temp43(2,2,1)
       temp42(2,3,2)=temp43(2,2,2)
       temp42(2,3,3)=temp43(3,2,2)
       temp43(1,1,2)=temp43(2,1,1)
       temp43(1,1,3)=temp43(3,1,1)
       temp43(1,2,1)=temp43(2,1,1)
       temp43(1,2,2)=temp43(2,2,1)
       temp43(1,2,3)=temp43(3,2,1)
       temp43(1,3,1)=temp43(3,1,1)
       temp43(1,3,2)=temp43(3,2,1)
       temp43(1,3,3)=temp43(3,3,1)
       temp43(2,1,2)=temp43(2,2,1)
       temp43(2,1,3)=temp43(3,2,1)
       temp43(2,2,3)=temp43(3,2,2)
       temp43(2,3,1)=temp43(3,2,1)
       temp43(2,3,2)=temp43(3,2,2)
       temp43(2,3,3)=temp43(3,3,2)
       temp43(3,1,2)=temp43(3,2,1)
       temp43(3,1,3)=temp43(3,3,1)
       temp43(3,2,3)=temp43(3,3,2)
c                Step3
       tempD40000=I8Z*(auxD40000+tempD400*F(5)-det4*temp002(k,l))
       temp002(1,1)=I12Z*(aux002(1,1)+4*F(4)*temp001(1)+F(5)*temp2(1,1)-
     &  det4*temp41(1,k,l)+8*tempD40000*ZZ(k,1,l,1))
       temp002(2,1)=I12Z*(aux002(2,1)+2*(F(6)*temp001(1)+F(4)*temp001(2)
     &  )+F(5)*temp2(2,1)-det4*temp42(1,k,l)+4*tempD40000*(ZZ(k,1,l,2)+Z
     &  Z(k,2,l,1)))
       temp002(2,2)=I12Z*(aux002(2,2)+4*F(6)*temp001(2)+F(5)*temp2(2,2)-
     &  det4*temp42(2,k,l)+8*tempD40000*ZZ(k,2,l,2))
       temp002(3,1)=I12Z*(aux002(3,1)+2*(F(7)*temp001(1)+F(4)*temp001(3)
     &  )+F(5)*temp2(3,1)-det4*temp43(1,k,l)+4*tempD40000*(ZZ(k,1,l,3)+Z
     &  Z(k,3,l,1)))
       temp002(3,2)=I12Z*(aux002(3,2)+2*(F(7)*temp001(2)+F(6)*temp001(3)
     &  )+F(5)*temp2(3,2)-det4*temp43(2,k,l)+4*tempD40000*(ZZ(k,2,l,3)+Z
     &  Z(k,3,l,2)))
       temp002(3,3)=I12Z*(aux002(3,3)+4*F(7)*temp001(3)+F(5)*temp2(3,3)-
     &  det4*temp43(3,k,l)+8*tempD40000*ZZ(k,3,l,3))
       temp002(1,2)=temp002(2,1)
       temp002(1,3)=temp002(3,1)
       temp002(2,3)=temp002(3,2)
       temp3(1,1,1)=IX*(aux3(1,1,1)+det4*temp41(1,1,jj)+6*temp002(1,1)*Z
     &  (jj,1))
       temp3(2,1,1)=IX*(aux3(2,1,1)+det4*temp42(1,1,jj)+4*temp002(2,1)*Z
     &  (jj,1)+2*temp002(1,1)*Z(jj,2))
       temp3(2,2,1)=IX*(aux3(2,2,1)+det4*temp42(2,1,jj)+2*temp002(2,2)*Z
     &  (jj,1)+4*temp002(2,1)*Z(jj,2))
       temp3(2,2,2)=IX*(aux3(2,2,2)+det4*temp42(2,2,jj)+6*temp002(2,2)*Z
     &  (jj,2))
       temp3(3,1,1)=IX*(aux3(3,1,1)+det4*temp43(1,1,jj)+4*temp002(3,1)*Z
     &  (jj,1)+2*temp002(1,1)*Z(jj,3))
       temp3(3,2,1)=IX*(aux3(3,2,1)+det4*temp43(2,1,jj)+2*(temp002(3,2)*
     &  Z(jj,1)+temp002(3,1)*Z(jj,2))+2*temp002(2,1)*Z(jj,3))
       temp3(3,2,2)=IX*(aux3(3,2,2)+det4*temp43(2,2,jj)+4*temp002(3,2)*Z
     &  (jj,2)+2*temp002(2,2)*Z(jj,3))
       temp3(3,3,1)=IX*(aux3(3,3,1)+det4*temp43(3,1,jj)+2*temp002(3,3)*Z
     &  (jj,1)+4*temp002(3,1)*Z(jj,3))
       temp3(3,3,2)=IX*(aux3(3,3,2)+det4*temp43(3,2,jj)+2*temp002(3,3)*Z
     &  (jj,2)+4*temp002(3,2)*Z(jj,3))
       temp3(3,3,3)=IX*(aux3(3,3,3)+det4*temp43(3,3,jj)+6*temp002(3,3)*Z
     &  (jj,3))
       temp3(1,1,2)=temp3(2,1,1)
       temp3(1,1,3)=temp3(3,1,1)
       temp3(1,2,1)=temp3(2,1,1)
       temp3(1,2,2)=temp3(2,2,1)
       temp3(1,2,3)=temp3(3,2,1)
       temp3(1,3,1)=temp3(3,1,1)
       temp3(1,3,2)=temp3(3,2,1)
       temp3(1,3,3)=temp3(3,3,1)
       temp3(2,1,2)=temp3(2,2,1)
       temp3(2,1,3)=temp3(3,2,1)
       temp3(2,2,3)=temp3(3,2,2)
       temp3(2,3,1)=temp3(3,2,1)
       temp3(2,3,2)=temp3(3,2,2)
       temp3(2,3,3)=temp3(3,3,2)
       temp3(3,1,2)=temp3(3,2,1)
       temp3(3,1,3)=temp3(3,3,1)
       temp3(3,2,3)=temp3(3,3,2)
c                Step4
       temp001(1)=I8Z*(aux001(1)+2*tempD400*F(4)+F(5)*temp1(1)-det4*temp
     &  3(1,k,l))
       temp001(2)=I8Z*(aux001(2)+2*tempD400*F(6)+F(5)*temp1(2)-det4*temp
     &  3(2,k,l))
       temp001(3)=I8Z*(aux001(3)+2*tempD400*F(7)+F(5)*temp1(3)-det4*temp
     &  3(3,k,l))
       temp2(1,1)=IX*(aux2(1,1)+det4*temp3(1,1,jj)+4*temp001(1)*Z(jj,1))
       temp2(2,1)=IX*(aux2(2,1)+det4*temp3(2,1,jj)+2*temp001(2)*Z(jj,1)+
     &  2*temp001(1)*Z(jj,2))
       temp2(2,2)=IX*(aux2(2,2)+det4*temp3(2,2,jj)+4*temp001(2)*Z(jj,2))
       temp2(3,1)=IX*(aux2(3,1)+det4*temp3(3,1,jj)+2*temp001(3)*Z(jj,1)+
     &  2*temp001(1)*Z(jj,3))
       temp2(3,2)=IX*(aux2(3,2)+det4*temp3(3,2,jj)+2*temp001(3)*Z(jj,2)+
     &  2*temp001(2)*Z(jj,3))
       temp2(3,3)=IX*(aux2(3,3)+det4*temp3(3,3,jj)+4*temp001(3)*Z(jj,3))
       temp2(1,2)=temp2(2,1)
       temp2(1,3)=temp2(3,1)
       temp2(2,3)=temp2(3,2)
c                Step5
       tempD400=I4Z*(auxD400+tempD40*F(5)-det4*temp2(k,l))
       temp1(1)=IX*(aux1(1)+det4*temp2(1,jj)+2*tempD400*Z(jj,1))
       temp1(2)=IX*(aux1(2)+det4*temp2(2,jj)+2*tempD400*Z(jj,2))
       temp1(3)=IX*(aux1(3)+det4*temp2(3,jj)+2*tempD400*Z(jj,3))
c                Step6
       tempD40=IX*(auxD40+det4*temp1(jj))
        D0=tempD40
        Dij(1,1)=temp1(1)
        Dij(2,1)=temp1(2)
        Dij(3,1)=temp1(3)
        Dij(1,2)=temp2(1,1)
        Dij(2,2)=temp2(2,2)
        Dij(3,2)=temp2(3,3)
        Dij(4,2)=temp2(2,1)
        Dij(5,2)=temp2(3,1)
        Dij(6,2)=temp2(3,2)
        Dij(7,2)=tempD400
        Dij(1,3)=temp3(1,1,1)
        Dij(2,3)=temp3(2,2,2)
        Dij(3,3)=temp3(3,3,3)
        Dij(4,3)=temp3(2,1,1)
        Dij(5,3)=temp3(3,1,1)
        Dij(6,3)=temp3(2,2,1)
        Dij(7,3)=temp3(3,3,1)
        Dij(8,3)=temp3(3,2,2)
        Dij(9,3)=temp3(3,3,2)
        Dij(10,3)=temp3(3,2,1)
        Dij(11,3)=temp001(1)
        Dij(12,3)=temp001(2)
        Dij(13,3)=temp001(3)
        Dij(7,1)=temp002(1,1)
        Dij(8,1)=temp002(2,2)
        Dij(9,1)=temp002(3,3)
        Dij(10,1)=temp002(2,1)
        Dij(11,1)=temp002(3,1)
        Dij(12,1)=temp002(3,2)
        Dij(13,1)=tempD40000
        Dij(1,4)=temp41(1,1,1)
         Dij(2,4)=temp42(2,2,2)
         Dij(3,4)=temp43(3,3,3)
         Dij(4,4)=temp42(1,1,1)
         Dij(5,4)=temp43(1,1,1)
         Dij(6,4)=temp42(2,1,1)
         Dij(7,4)=temp43(2,1,1)
         Dij(8,4)=temp43(3,1,1)
         Dij(9,4)=temp42(2,2,1)
         Dij(10,4)=temp43(2,2,1)
         Dij(11,4)=temp43(3,2,1)
         Dij(12,4)=temp43(3,3,1)
         Dij(13,4)=temp43(2,2,2)
         Dij(14,4)=temp43(3,2,2)
         Dij(15,4)=temp43(3,3,2)
         Dij(16,4)=temp002(1,1)
         Dij(17,4)=temp002(2,2)
         Dij(18,4)=temp002(3,3)
         Dij(19,4)=temp002(2,1)
         Dij(20,4)=temp002(3,1)
         Dij(21,4)=temp002(3,2)
         Dij(22,4)=tempD40000
        if(printmy) then
       print*, "tempD40 5",tempD40
       endif
       if(order.eq.5) goto 500
c                Iteration6
c                Step1
       S30000001(1)=-2*Cij234(9,4)
       S30000001(2)=2*Cij234(11,5)
       S30000001(3)=2*Cij234(12,5)
       S3h0000001(1)=Cij134(16,6)-Cij234(16,6)
       S3h0000001(2)=Cij124(16,6)-Cij134(16,6)
       S3h0000001(3)=Cij123(16,6)-Cij124(16,6)
       S3h0000311(1)=Cij134(13,6)-Cij234(9,4)
       S3h0000311(2)=Cij134(13,6)+Cij234(11,5)
       S3h0000311(3)=Cij134(15,6)+Cij234(12,5)
       S3h0000312(1)=Cij134(13,6)+Cij234(11,5)
       S3h0000312(2)=Cij134(13,6)-Cij234(13,6)
       S3h0000312(3)=Cij134(15,6)-Cij234(15,6)
       S3h0000313(1)=Cij134(15,6)+Cij234(12,5)
       S3h0000313(2)=Cij134(15,6)-Cij234(15,6)
       S3h0000313(3)=Cij134(14,6)-Cij234(14,6)
       S3h0000321(1)=Cij124(13,6)-Cij134(13,6)
       S3h0000321(2)=Cij124(15,6)-Cij134(13,6)
       S3h0000321(3)=Cij124(15,6)-Cij134(15,6)
       S3h0000322(1)=Cij124(15,6)-Cij134(13,6)
       S3h0000322(2)=Cij124(14,6)-Cij134(13,6)
       S3h0000322(3)=Cij124(14,6)-Cij134(15,6)
       S3h0000323(1)=Cij124(15,6)-Cij134(15,6)
       S3h0000323(2)=Cij124(14,6)-Cij134(15,6)
       S3h0000323(3)=Cij124(14,6)-Cij134(14,6)
       S3h0000331(1)=Cij123(13,6)-Cij124(13,6)
       S3h0000331(2)=Cij123(15,6)-Cij124(15,6)
       S3h0000331(3)=-Cij124(15,6)
       S3h0000332(1)=Cij123(15,6)-Cij124(15,6)
       S3h0000332(2)=Cij123(14,6)-Cij124(14,6)
       S3h0000332(3)=-Cij124(14,6)
       S3h0000333(1)=-Cij124(15,6)
       S3h0000333(2)=-Cij124(14,6)
       S3h0000333(3)=-Cij124(14,6)
       aux0000001(1)=-(F(1)*S3h000021(1))-F(2)*S3h000022(1)-F(3)*S3h0000
     &  23(1)+S3h0000311(k)*Z(1,l)+S3h0000321(k)*Z(2,l)+S3h0000331(k)*Z(
     &  3,l)+(Inv601+S30000001(1)-S3h0000311(1)-S3h0000322(1)-S3h0000333
     &  (1))*Z(k,l)-2*S3h0000001(1)*ZZ(k,1,l,1)-2*S3h0000001(2)*ZZ(k,1,l
     &  ,2)-2*S3h0000001(3)*ZZ(k,1,l,3)
       aux0000001(2)=-(F(1)*S3h000021(2))-F(2)*S3h000022(2)-F(3)*S3h0000
     &  23(2)+S3h0000312(k)*Z(1,l)+S3h0000322(k)*Z(2,l)+S3h0000332(k)*Z(
     &  3,l)+(Inv602+S30000001(2)-S3h0000312(1)-S3h0000322(2)-S3h0000333
     &  (2))*Z(k,l)-2*S3h0000001(1)*ZZ(k,2,l,1)-2*S3h0000001(2)*ZZ(k,2,l
     &  ,2)-2*S3h0000001(3)*ZZ(k,2,l,3)
       aux0000001(3)=-(F(1)*S3h000021(3))-F(2)*S3h000022(3)-F(3)*S3h0000
     &  23(3)+S3h0000313(k)*Z(1,l)+S3h0000323(k)*Z(2,l)+S3h0000333(k)*Z(
     &  3,l)+(Inv603+S30000001(3)-S3h0000313(1)-S3h0000323(2)-S3h0000333
     &  (3))*Z(k,l)-2*S3h0000001(1)*ZZ(k,3,l,1)-2*S3h0000001(2)*ZZ(k,3,l
     &  ,2)-2*S3h0000001(3)*ZZ(k,3,l,3)
       temp0000001(1)=I16Z*(aux0000001(1)+2*tempD4000000*F(4)+F(5)*temp0
     &  0001(1))
       temp0000001(2)=I16Z*(aux0000001(2)+2*tempD4000000*F(6)+F(5)*temp0
     &  0001(2))
       temp0000001(3)=I16Z*(aux0000001(3)+2*tempD4000000*F(7)+F(5)*temp0
     &  0001(3))
       S3h0051111(1)=Cij134(8,6)-Cij234(4,2)
       S3h0051111(2)=Cij134(8,6)+Cij234(5,3)
       S3h0051111(3)=Cij134(10,6)+Cij234(6,3)
       S3h0051211(1)=Cij134(8,6)+Cij234(5,3)
       S3h0051211(2)=Cij134(8,6)-Cij234(6,4)
       S3h0051211(3)=Cij134(10,6)-Cij234(8,4)
       S3h0051221(1)=Cij134(8,6)-Cij234(6,4)
       S3h0051221(2)=Cij134(8,6)+Cij234(7,5)
       S3h0051221(3)=Cij134(10,6)+Cij234(9,5)
       S3h0051222(1)=Cij134(8,6)+Cij234(7,5)
       S3h0051222(2)=Cij134(8,6)-Cij234(8,6)
       S3h0051222(3)=Cij134(10,6)-Cij234(10,6)
       S3h0051311(1)=Cij134(10,6)+Cij234(6,3)
       S3h0051311(2)=Cij134(10,6)-Cij234(8,4)
       S3h0051311(3)=Cij134(11,6)-Cij234(7,4)
       S3h0051321(1)=Cij134(10,6)-Cij234(8,4)
       S3h0051321(2)=Cij134(10,6)+Cij234(9,5)
       S3h0051321(3)=Cij134(11,6)+Cij234(10,5)
       S3h0051322(1)=Cij134(10,6)+Cij234(9,5)
       S3h0051322(2)=Cij134(10,6)-Cij234(10,6)
       S3h0051322(3)=Cij134(11,6)-Cij234(11,6)
       S3h0051331(1)=Cij134(11,6)-Cij234(7,4)
       S3h0051331(2)=Cij134(11,6)+Cij234(10,5)
       S3h0051331(3)=Cij134(12,6)+Cij234(8,5)
       S3h0051332(1)=Cij134(11,6)+Cij234(10,5)
       S3h0051332(2)=Cij134(11,6)-Cij234(11,6)
       S3h0051332(3)=Cij134(12,6)-Cij234(12,6)
       S3h0051333(1)=Cij134(12,6)+Cij234(8,5)
       S3h0051333(2)=Cij134(12,6)-Cij234(12,6)
       S3h0051333(3)=Cij134(9,6)-Cij234(9,6)
       S3h0052111(1)=Cij124(8,6)-Cij134(8,6)
       S3h0052111(2)=Cij124(10,6)-Cij134(8,6)
       S3h0052111(3)=Cij124(10,6)-Cij134(10,6)
       S3h0052211(1)=Cij124(10,6)-Cij134(8,6)
       S3h0052211(2)=Cij124(11,6)-Cij134(8,6)
       S3h0052211(3)=Cij124(11,6)-Cij134(10,6)
       S3h0052221(1)=Cij124(11,6)-Cij134(8,6)
       S3h0052221(2)=Cij124(12,6)-Cij134(8,6)
       S3h0052221(3)=Cij124(12,6)-Cij134(10,6)
       S3h0052222(1)=Cij124(12,6)-Cij134(8,6)
       S3h0052222(2)=Cij124(9,6)-Cij134(8,6)
       S3h0052222(3)=Cij124(9,6)-Cij134(10,6)
       S3h0052311(1)=Cij124(10,6)-Cij134(10,6)
       S3h0052311(2)=Cij124(11,6)-Cij134(10,6)
       S3h0052311(3)=Cij124(11,6)-Cij134(11,6)
       S3h0052321(1)=Cij124(11,6)-Cij134(10,6)
       S3h0052321(2)=Cij124(12,6)-Cij134(10,6)
       S3h0052321(3)=Cij124(12,6)-Cij134(11,6)
       S3h0052322(1)=Cij124(12,6)-Cij134(10,6)
       S3h0052322(2)=Cij124(9,6)-Cij134(10,6)
       S3h0052322(3)=Cij124(9,6)-Cij134(11,6)
       S3h0052331(1)=Cij124(11,6)-Cij134(11,6)
       S3h0052331(2)=Cij124(12,6)-Cij134(11,6)
       S3h0052331(3)=Cij124(12,6)-Cij134(12,6)
       S3h0052332(1)=Cij124(12,6)-Cij134(11,6)
       S3h0052332(2)=Cij124(9,6)-Cij134(11,6)
       S3h0052332(3)=Cij124(9,6)-Cij134(12,6)
       S3h0052333(1)=Cij124(12,6)-Cij134(12,6)
       S3h0052333(2)=Cij124(9,6)-Cij134(12,6)
       S3h0052333(3)=Cij124(9,6)-Cij134(9,6)
       S3h0053111(1)=Cij123(8,6)-Cij124(8,6)
       S3h0053111(2)=Cij123(10,6)-Cij124(10,6)
       S3h0053111(3)=-Cij124(10,6)
       S3h0053211(1)=Cij123(10,6)-Cij124(10,6)
       S3h0053211(2)=Cij123(11,6)-Cij124(11,6)
       S3h0053211(3)=-Cij124(11,6)
       S3h0053221(1)=Cij123(11,6)-Cij124(11,6)
       S3h0053221(2)=Cij123(12,6)-Cij124(12,6)
       S3h0053221(3)=-Cij124(12,6)
       S3h0053222(1)=Cij123(12,6)-Cij124(12,6)
       S3h0053222(2)=Cij123(9,6)-Cij124(9,6)
       S3h0053222(3)=-Cij124(9,6)
       S3h0053311(1)=-Cij124(10,6)
       S3h0053311(2)=-Cij124(11,6)
       S3h0053311(3)=-Cij124(11,6)
       S3h0053321(1)=-Cij124(11,6)
       S3h0053321(2)=-Cij124(12,6)
       S3h0053321(3)=-Cij124(12,6)
       S3h0053322(1)=-Cij124(12,6)
       S3h0053322(2)=-Cij124(9,6)
       S3h0053322(3)=-Cij124(9,6)
       S3h0053331(1)=-Cij124(11,6)
       S3h0053331(2)=-Cij124(12,6)
       S3h0053331(3)=-Cij124(12,6)
       S3h0053332(1)=-Cij124(12,6)
       S3h0053332(2)=-Cij124(9,6)
       S3h0053332(3)=-Cij124(9,6)
       S3h0053333(1)=-Cij124(12,6)
       S3h0053333(2)=-Cij124(9,6)
       S3h0053333(3)=-Cij124(9,6)
       S30000311(1)=-2*Cij234(4,2)
       S30000311(2)=2*Cij234(5,3)
       S30000311(3)=2*Cij234(6,3)
       S30000312(1)=2*Cij234(5,3)
       S30000312(2)=-2*Cij234(6,4)
       S30000312(3)=-2*Cij234(8,4)
       S30000313(1)=2*Cij234(6,3)
       S30000313(2)=-2*Cij234(8,4)
       S30000313(3)=-2*Cij234(7,4)
       S30000321(1)=2*Cij234(5,3)
       S30000321(2)=-2*Cij234(6,4)
       S30000321(3)=-2*Cij234(8,4)
       S30000322(1)=-2*Cij234(6,4)
       S30000322(2)=2*Cij234(7,5)
       S30000322(3)=2*Cij234(9,5)
       S30000323(1)=-2*Cij234(8,4)
       S30000323(2)=2*Cij234(9,5)
       S30000323(3)=2*Cij234(10,5)
       S30000331(1)=2*Cij234(6,3)
       S30000331(2)=-2*Cij234(8,4)
       S30000331(3)=-2*Cij234(7,4)
       S30000332(1)=-2*Cij234(8,4)
       S30000332(2)=2*Cij234(9,5)
       S30000332(3)=2*Cij234(10,5)
       S30000333(1)=-2*Cij234(7,4)
       S30000333(2)=2*Cij234(10,5)
       S30000333(3)=2*Cij234(8,5)
       aux00003(1,1,1)=-(F(1)*S3h004111(1))-F(2)*S3h004211(1)-F(3)*S3h00
     &  4311(1)+S3h0051111(k)*Z(1,l)+S3h0052111(k)*Z(2,l)+S3h0053111(k)*
     &  Z(3,l)+(Inv40111+S30000311(1)-S3h0051111(1)-S3h0052211(1)-S3h005
     &  3311(1))*Z(k,l)-6*S3h0000311(1)*ZZ(k,1,l,1)-6*S3h0000321(1)*ZZ(k
     &  ,1,l,2)-6*S3h0000331(1)*ZZ(k,1,l,3)
       aux00003(2,1,1)=-(F(1)*S3h004121(1))-F(2)*S3h004221(1)-F(3)*S3h00
     &  4321(1)+S3h0051211(k)*Z(1,l)+S3h0052211(k)*Z(2,l)+S3h0053211(k)*
     &  Z(3,l)+(Inv40211+S30000321(1)-S3h0051211(1)-S3h0052221(1)-S3h005
     &  3321(1))*Z(k,l)-4*S3h0000312(1)*ZZ(k,1,l,1)-4*S3h0000322(1)*ZZ(k
     &  ,1,l,2)-4*S3h0000332(1)*ZZ(k,1,l,3)-2*S3h0000311(1)*ZZ(k,2,l,1)-
     &  2*S3h0000321(1)*ZZ(k,2,l,2)-2*S3h0000331(1)*ZZ(k,2,l,3)
       aux00003(2,2,1)=-(F(1)*S3h004122(1))-F(2)*S3h004222(1)-F(3)*S3h00
     &  4322(1)+S3h0051221(k)*Z(1,l)+S3h0052221(k)*Z(2,l)+S3h0053221(k)*
     &  Z(3,l)+(Inv40221+S30000322(1)-S3h0051221(1)-S3h0052222(1)-S3h005
     &  3322(1))*Z(k,l)-2*S3h0000312(2)*ZZ(k,1,l,1)-2*S3h0000322(2)*ZZ(k
     &  ,1,l,2)-2*S3h0000332(2)*ZZ(k,1,l,3)-4*S3h0000312(1)*ZZ(k,2,l,1)-
     &  4*S3h0000322(1)*ZZ(k,2,l,2)-4*S3h0000332(1)*ZZ(k,2,l,3)
       aux00003(2,2,2)=-(F(1)*S3h004122(2))-F(2)*S3h004222(2)-F(3)*S3h00
     &  4322(2)+S3h0051222(k)*Z(1,l)+S3h0052222(k)*Z(2,l)+S3h0053222(k)*
     &  Z(3,l)+(Inv40222+S30000322(2)-S3h0051222(1)-S3h0052222(2)-S3h005
     &  3322(2))*Z(k,l)-6*S3h0000312(2)*ZZ(k,2,l,1)-6*S3h0000322(2)*ZZ(k
     &  ,2,l,2)-6*S3h0000332(2)*ZZ(k,2,l,3)
       aux00003(3,1,1)=-(F(1)*S3h004131(1))-F(2)*S3h004231(1)-F(3)*S3h00
     &  4331(1)+S3h0051311(k)*Z(1,l)+S3h0052311(k)*Z(2,l)+S3h0053311(k)*
     &  Z(3,l)+(Inv40311+S30000331(1)-S3h0051311(1)-S3h0052321(1)-S3h005
     &  3331(1))*Z(k,l)-4*S3h0000313(1)*ZZ(k,1,l,1)-4*S3h0000323(1)*ZZ(k
     &  ,1,l,2)-4*S3h0000333(1)*ZZ(k,1,l,3)-2*S3h0000311(1)*ZZ(k,3,l,1)-
     &  2*S3h0000321(1)*ZZ(k,3,l,2)-2*S3h0000331(1)*ZZ(k,3,l,3)
       aux00003(3,2,1)=-(F(1)*S3h004132(1))-F(2)*S3h004232(1)-F(3)*S3h00
     &  4332(1)+S3h0051321(k)*Z(1,l)+S3h0052321(k)*Z(2,l)+S3h0053321(k)*
     &  Z(3,l)+(Inv40321+S30000332(1)-S3h0051321(1)-S3h0052322(1)-S3h005
     &  3332(1))*Z(k,l)-2*(S3h0000313(2)*ZZ(k,1,l,1)+S3h0000323(2)*ZZ(k,
     &  1,l,2)+S3h0000333(2)*ZZ(k,1,l,3)+S3h0000313(1)*ZZ(k,2,l,1)+S3h00
     &  00323(1)*ZZ(k,2,l,2)+S3h0000333(1)*ZZ(k,2,l,3)+S3h0000312(1)*ZZ(
     &  k,3,l,1)+S3h0000322(1)*ZZ(k,3,l,2)+S3h0000332(1)*ZZ(k,3,l,3))
       aux00003(3,2,2)=-(F(1)*S3h004132(2))-F(2)*S3h004232(2)-F(3)*S3h00
     &  4332(2)+S3h0051322(k)*Z(1,l)+S3h0052322(k)*Z(2,l)+S3h0053322(k)*
     &  Z(3,l)+(Inv40322+S30000332(2)-S3h0051322(1)-S3h0052322(2)-S3h005
     &  3332(2))*Z(k,l)-4*S3h0000313(2)*ZZ(k,2,l,1)-4*S3h0000323(2)*ZZ(k
     &  ,2,l,2)-4*S3h0000333(2)*ZZ(k,2,l,3)-2*S3h0000312(2)*ZZ(k,3,l,1)-
     &  2*S3h0000322(2)*ZZ(k,3,l,2)-2*S3h0000332(2)*ZZ(k,3,l,3)
       aux00003(3,3,1)=-(F(1)*S3h004133(1))-F(2)*S3h004233(1)-F(3)*S3h00
     &  4333(1)+S3h0051331(k)*Z(1,l)+S3h0052331(k)*Z(2,l)+S3h0053331(k)*
     &  Z(3,l)+(Inv40331+S30000333(1)-S3h0051331(1)-S3h0052332(1)-S3h005
     &  3333(1))*Z(k,l)-2*S3h0000313(3)*ZZ(k,1,l,1)-2*S3h0000323(3)*ZZ(k
     &  ,1,l,2)-2*S3h0000333(3)*ZZ(k,1,l,3)-4*S3h0000313(1)*ZZ(k,3,l,1)-
     &  4*S3h0000323(1)*ZZ(k,3,l,2)-4*S3h0000333(1)*ZZ(k,3,l,3)
       aux00003(3,3,2)=-(F(1)*S3h004133(2))-F(2)*S3h004233(2)-F(3)*S3h00
     &  4333(2)+S3h0051332(k)*Z(1,l)+S3h0052332(k)*Z(2,l)+S3h0053332(k)*
     &  Z(3,l)+(Inv40332+S30000333(2)-S3h0051332(1)-S3h0052332(2)-S3h005
     &  3333(2))*Z(k,l)-2*S3h0000313(3)*ZZ(k,2,l,1)-2*S3h0000323(3)*ZZ(k
     &  ,2,l,2)-2*S3h0000333(3)*ZZ(k,2,l,3)-4*S3h0000313(2)*ZZ(k,3,l,1)-
     &  4*S3h0000323(2)*ZZ(k,3,l,2)-4*S3h0000333(2)*ZZ(k,3,l,3)
       aux00003(3,3,3)=-(F(1)*S3h004133(3))-F(2)*S3h004233(3)-F(3)*S3h00
     &  4333(3)+S3h0051333(k)*Z(1,l)+S3h0052333(k)*Z(2,l)+S3h0053333(k)*
     &  Z(3,l)+(Inv40333+S30000333(3)-S3h0051333(1)-S3h0052333(2)-S3h005
     &  3333(3))*Z(k,l)-6*S3h0000313(3)*ZZ(k,3,l,1)-6*S3h0000323(3)*ZZ(k
     &  ,3,l,2)-6*S3h0000333(3)*ZZ(k,3,l,3)
       temp00003(1,1,1)=I20Z*(aux00003(1,1,1)+6*F(4)*temp00002(1,1)+F(5)
     &  *temp003(1,1,1)+24*temp0000001(1)*ZZ(k,1,l,1))
       temp00003(2,1,1)=I20Z*(aux00003(2,1,1)+2*F(6)*temp00002(1,1)+4*F(
     &  4)*temp00002(2,1)+F(5)*temp003(2,1,1)+8*(temp0000001(2)*ZZ(k,1,l
     &  ,1)+temp0000001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))))
       temp00003(2,2,1)=I20Z*(aux00003(2,2,1)+4*F(6)*temp00002(2,1)+2*F(
     &  4)*temp00002(2,2)+F(5)*temp003(2,2,1)+8*(temp0000001(2)*(ZZ(k,1,
     &  l,2)+ZZ(k,2,l,1))+temp0000001(1)*ZZ(k,2,l,2)))
       temp00003(2,2,2)=I20Z*(aux00003(2,2,2)+6*F(6)*temp00002(2,2)+F(5)
     &  *temp003(2,2,2)+24*temp0000001(2)*ZZ(k,2,l,2))
       temp00003(3,1,1)=I20Z*(aux00003(3,1,1)+2*F(7)*temp00002(1,1)+4*F(
     &  4)*temp00002(3,1)+F(5)*temp003(3,1,1)+8*(temp0000001(3)*ZZ(k,1,l
     &  ,1)+temp0000001(1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))))
       temp00003(3,2,1)=I20Z*(aux00003(3,2,1)+2*F(7)*temp00002(2,1)+2*F(
     &  6)*temp00002(3,1)+2*F(4)*temp00002(3,2)+F(5)*temp003(3,2,1)+4*(t
     &  emp0000001(3)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0000001(2)*(ZZ(k,1,l
     &  ,3)+ZZ(k,3,l,1)))+4*temp0000001(1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp00003(3,2,2)=I20Z*(aux00003(3,2,2)+2*F(7)*temp00002(2,2)+4*F(
     &  6)*temp00002(3,2)+F(5)*temp003(3,2,2)+8*(temp0000001(3)*ZZ(k,2,l
     &  ,2)+temp0000001(2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))))
       temp00003(3,3,1)=I20Z*(aux00003(3,3,1)+4*F(7)*temp00002(3,1)+2*F(
     &  4)*temp00002(3,3)+F(5)*temp003(3,3,1)+8*(temp0000001(3)*(ZZ(k,1,
     &  l,3)+ZZ(k,3,l,1))+temp0000001(1)*ZZ(k,3,l,3)))
       temp00003(3,3,2)=I20Z*(aux00003(3,3,2)+4*F(7)*temp00002(3,2)+2*F(
     &  6)*temp00002(3,3)+F(5)*temp003(3,3,2)+8*(temp0000001(3)*(ZZ(k,2,
     &  l,3)+ZZ(k,3,l,2))+temp0000001(2)*ZZ(k,3,l,3)))
       temp00003(3,3,3)=I20Z*(aux00003(3,3,3)+6*F(7)*temp00002(3,3)+F(5)
     &  *temp003(3,3,3)+24*temp0000001(3)*ZZ(k,3,l,3))
       temp00003(1,1,2)=temp00003(2,1,1)
       temp00003(1,1,3)=temp00003(3,1,1)
       temp00003(1,2,1)=temp00003(2,1,1)
       temp00003(1,2,2)=temp00003(2,2,1)
       temp00003(1,2,3)=temp00003(3,2,1)
       temp00003(1,3,1)=temp00003(3,1,1)
       temp00003(1,3,2)=temp00003(3,2,1)
       temp00003(1,3,3)=temp00003(3,3,1)
       temp00003(2,1,2)=temp00003(2,2,1)
       temp00003(2,1,3)=temp00003(3,2,1)
       temp00003(2,2,3)=temp00003(3,2,2)
       temp00003(2,3,1)=temp00003(3,2,1)
       temp00003(2,3,2)=temp00003(3,2,2)
       temp00003(2,3,3)=temp00003(3,3,2)
       temp00003(3,1,2)=temp00003(3,2,1)
       temp00003(3,1,3)=temp00003(3,3,1)
       temp00003(3,2,3)=temp00003(3,3,2)
       S30051111(1)=-2*C0234
       S30051111(2)=2*Cij234(1,1)
       S30051111(3)=2*Cij234(2,1)
       S30051211(1)=2*Cij234(1,1)
       S30051211(2)=-2*Cij234(1,2)
       S30051211(3)=-2*Cij234(3,2)
       S30051221(1)=-2*Cij234(1,2)
       S30051221(2)=2*Cij234(1,3)
       S30051221(3)=2*Cij234(3,3)
       S30051222(1)=2*Cij234(1,3)
       S30051222(2)=-2*Cij234(1,4)
       S30051222(3)=-2*Cij234(3,4)
       S30051311(1)=2*Cij234(2,1)
       S30051311(2)=-2*Cij234(3,2)
       S30051311(3)=-2*Cij234(2,2)
       S30051321(1)=-2*Cij234(3,2)
       S30051321(2)=2*Cij234(3,3)
       S30051321(3)=2*Cij234(4,3)
       S30051322(1)=2*Cij234(3,3)
       S30051322(2)=-2*Cij234(3,4)
       S30051322(3)=-2*Cij234(4,4)
       S30051331(1)=-2*Cij234(2,2)
       S30051331(2)=2*Cij234(4,3)
       S30051331(3)=2*Cij234(2,3)
       S30051332(1)=2*Cij234(4,3)
       S30051332(2)=-2*Cij234(4,4)
       S30051332(3)=-2*Cij234(5,4)
       S30051333(1)=2*Cij234(2,3)
       S30051333(2)=-2*Cij234(5,4)
       S30051333(3)=-2*Cij234(2,4)
       S30052111(1)=2*Cij234(1,1)
       S30052111(2)=-2*Cij234(1,2)
       S30052111(3)=-2*Cij234(3,2)
       S30052211(1)=-2*Cij234(1,2)
       S30052211(2)=2*Cij234(1,3)
       S30052211(3)=2*Cij234(3,3)
       S30052221(1)=2*Cij234(1,3)
       S30052221(2)=-2*Cij234(1,4)
       S30052221(3)=-2*Cij234(3,4)
       S30052222(1)=-2*Cij234(1,4)
       S30052222(2)=2*Cij234(1,5)
       S30052222(3)=2*Cij234(3,5)
       S30052311(1)=-2*Cij234(3,2)
       S30052311(2)=2*Cij234(3,3)
       S30052311(3)=2*Cij234(4,3)
       S30052321(1)=2*Cij234(3,3)
       S30052321(2)=-2*Cij234(3,4)
       S30052321(3)=-2*Cij234(4,4)
       S30052322(1)=-2*Cij234(3,4)
       S30052322(2)=2*Cij234(3,5)
       S30052322(3)=2*Cij234(4,5)
       S30052331(1)=2*Cij234(4,3)
       S30052331(2)=-2*Cij234(4,4)
       S30052331(3)=-2*Cij234(5,4)
       S30052332(1)=-2*Cij234(4,4)
       S30052332(2)=2*Cij234(4,5)
       S30052332(3)=2*Cij234(5,5)
       S30052333(1)=-2*Cij234(5,4)
       S30052333(2)=2*Cij234(5,5)
       S30052333(3)=2*Cij234(6,5)
       S30053111(1)=2*Cij234(2,1)
       S30053111(2)=-2*Cij234(3,2)
       S30053111(3)=-2*Cij234(2,2)
       S30053211(1)=-2*Cij234(3,2)
       S30053211(2)=2*Cij234(3,3)
       S30053211(3)=2*Cij234(4,3)
       S30053221(1)=2*Cij234(3,3)
       S30053221(2)=-2*Cij234(3,4)
       S30053221(3)=-2*Cij234(4,4)
       S30053222(1)=-2*Cij234(3,4)
       S30053222(2)=2*Cij234(3,5)
       S30053222(3)=2*Cij234(4,5)
       S30053311(1)=-2*Cij234(2,2)
       S30053311(2)=2*Cij234(4,3)
       S30053311(3)=2*Cij234(2,3)
       S30053321(1)=2*Cij234(4,3)
       S30053321(2)=-2*Cij234(4,4)
       S30053321(3)=-2*Cij234(5,4)
       S30053322(1)=-2*Cij234(4,4)
       S30053322(2)=2*Cij234(4,5)
       S30053322(3)=2*Cij234(5,5)
       S30053331(1)=2*Cij234(2,3)
       S30053331(2)=-2*Cij234(5,4)
       S30053331(3)=-2*Cij234(2,4)
       S30053332(1)=-2*Cij234(5,4)
       S30053332(2)=2*Cij234(5,5)
       S30053332(3)=2*Cij234(6,5)
       S30053333(1)=-2*Cij234(2,4)
       S30053333(2)=2*Cij234(6,5)
       S30053333(3)=2*Cij234(2,5)
       S37111111(1)=-C0234+Cij134(1,6)
       S37111111(2)=Cij134(1,6)+Cij234(1,1)
       S37111111(3)=Cij134(3,6)+Cij234(2,1)
       S37121111(1)=Cij134(1,6)+Cij234(1,1)
       S37121111(2)=Cij134(1,6)-Cij234(1,2)
       S37121111(3)=Cij134(3,6)-Cij234(3,2)
       S37122111(1)=Cij134(1,6)-Cij234(1,2)
       S37122111(2)=Cij134(1,6)+Cij234(1,3)
       S37122111(3)=Cij134(3,6)+Cij234(3,3)
       S37122211(1)=Cij134(1,6)+Cij234(1,3)
       S37122211(2)=Cij134(1,6)-Cij234(1,4)
       S37122211(3)=Cij134(3,6)-Cij234(3,4)
       S37122221(1)=Cij134(1,6)-Cij234(1,4)
       S37122221(2)=Cij134(1,6)+Cij234(1,5)
       S37122221(3)=Cij134(3,6)+Cij234(3,5)
       S37122222(1)=Cij134(1,6)+Cij234(1,5)
       S37122222(2)=Cij134(1,6)-Cij234(1,6)
       S37122222(3)=Cij134(3,6)-Cij234(3,6)
       S37131111(1)=Cij134(3,6)+Cij234(2,1)
       S37131111(2)=Cij134(3,6)-Cij234(3,2)
       S37131111(3)=Cij134(4,6)-Cij234(2,2)
       S37132111(1)=Cij134(3,6)-Cij234(3,2)
       S37132111(2)=Cij134(3,6)+Cij234(3,3)
       S37132111(3)=Cij134(4,6)+Cij234(4,3)
       S37132211(1)=Cij134(3,6)+Cij234(3,3)
       S37132211(2)=Cij134(3,6)-Cij234(3,4)
       S37132211(3)=Cij134(4,6)-Cij234(4,4)
       S37132221(1)=Cij134(3,6)-Cij234(3,4)
       S37132221(2)=Cij134(3,6)+Cij234(3,5)
       S37132221(3)=Cij134(4,6)+Cij234(4,5)
       S37132222(1)=Cij134(3,6)+Cij234(3,5)
       S37132222(2)=Cij134(3,6)-Cij234(3,6)
       S37132222(3)=Cij134(4,6)-Cij234(4,6)
       S37133111(1)=Cij134(4,6)-Cij234(2,2)
       S37133111(2)=Cij134(4,6)+Cij234(4,3)
       S37133111(3)=Cij134(5,6)+Cij234(2,3)
       S37133211(1)=Cij134(4,6)+Cij234(4,3)
       S37133211(2)=Cij134(4,6)-Cij234(4,4)
       S37133211(3)=Cij134(5,6)-Cij234(5,4)
       S37133221(1)=Cij134(4,6)-Cij234(4,4)
       S37133221(2)=Cij134(4,6)+Cij234(4,5)
       S37133221(3)=Cij134(5,6)+Cij234(5,5)
       S37133222(1)=Cij134(4,6)+Cij234(4,5)
       S37133222(2)=Cij134(4,6)-Cij234(4,6)
       S37133222(3)=Cij134(5,6)-Cij234(5,6)
       S37133311(1)=Cij134(5,6)+Cij234(2,3)
       S37133311(2)=Cij134(5,6)-Cij234(5,4)
       S37133311(3)=Cij134(6,6)-Cij234(2,4)
       S37133321(1)=Cij134(5,6)-Cij234(5,4)
       S37133321(2)=Cij134(5,6)+Cij234(5,5)
       S37133321(3)=Cij134(6,6)+Cij234(6,5)
       S37133322(1)=Cij134(5,6)+Cij234(5,5)
       S37133322(2)=Cij134(5,6)-Cij234(5,6)
       S37133322(3)=Cij134(6,6)-Cij234(6,6)
       S37133331(1)=Cij134(6,6)-Cij234(2,4)
       S37133331(2)=Cij134(6,6)+Cij234(6,5)
       S37133331(3)=Cij134(7,6)+Cij234(2,5)
       S37133332(1)=Cij134(6,6)+Cij234(6,5)
       S37133332(2)=Cij134(6,6)-Cij234(6,6)
       S37133332(3)=Cij134(7,6)-Cij234(7,6)
       S37133333(1)=Cij134(7,6)+Cij234(2,5)
       S37133333(2)=Cij134(7,6)-Cij234(7,6)
       S37133333(3)=Cij134(2,6)-Cij234(2,6)
       S37211111(1)=Cij124(1,6)-Cij134(1,6)
       S37211111(2)=Cij124(3,6)-Cij134(1,6)
       S37211111(3)=Cij124(3,6)-Cij134(3,6)
       S37221111(1)=Cij124(3,6)-Cij134(1,6)
       S37221111(2)=Cij124(4,6)-Cij134(1,6)
       S37221111(3)=Cij124(4,6)-Cij134(3,6)
       S37222111(1)=Cij124(4,6)-Cij134(1,6)
       S37222111(2)=Cij124(5,6)-Cij134(1,6)
       S37222111(3)=Cij124(5,6)-Cij134(3,6)
       S37222211(1)=Cij124(5,6)-Cij134(1,6)
       S37222211(2)=Cij124(6,6)-Cij134(1,6)
       S37222211(3)=Cij124(6,6)-Cij134(3,6)
       S37222221(1)=Cij124(6,6)-Cij134(1,6)
       S37222221(2)=Cij124(7,6)-Cij134(1,6)
       S37222221(3)=Cij124(7,6)-Cij134(3,6)
       S37222222(1)=Cij124(7,6)-Cij134(1,6)
       S37222222(2)=Cij124(2,6)-Cij134(1,6)
       S37222222(3)=Cij124(2,6)-Cij134(3,6)
       S37231111(1)=Cij124(3,6)-Cij134(3,6)
       S37231111(2)=Cij124(4,6)-Cij134(3,6)
       S37231111(3)=Cij124(4,6)-Cij134(4,6)
       S37232111(1)=Cij124(4,6)-Cij134(3,6)
       S37232111(2)=Cij124(5,6)-Cij134(3,6)
       S37232111(3)=Cij124(5,6)-Cij134(4,6)
       S37232211(1)=Cij124(5,6)-Cij134(3,6)
       S37232211(2)=Cij124(6,6)-Cij134(3,6)
       S37232211(3)=Cij124(6,6)-Cij134(4,6)
       S37232221(1)=Cij124(6,6)-Cij134(3,6)
       S37232221(2)=Cij124(7,6)-Cij134(3,6)
       S37232221(3)=Cij124(7,6)-Cij134(4,6)
       S37232222(1)=Cij124(7,6)-Cij134(3,6)
       S37232222(2)=Cij124(2,6)-Cij134(3,6)
       S37232222(3)=Cij124(2,6)-Cij134(4,6)
       S37233111(1)=Cij124(4,6)-Cij134(4,6)
       S37233111(2)=Cij124(5,6)-Cij134(4,6)
       S37233111(3)=Cij124(5,6)-Cij134(5,6)
       S37233211(1)=Cij124(5,6)-Cij134(4,6)
       S37233211(2)=Cij124(6,6)-Cij134(4,6)
       S37233211(3)=Cij124(6,6)-Cij134(5,6)
       S37233221(1)=Cij124(6,6)-Cij134(4,6)
       S37233221(2)=Cij124(7,6)-Cij134(4,6)
       S37233221(3)=Cij124(7,6)-Cij134(5,6)
       S37233222(1)=Cij124(7,6)-Cij134(4,6)
       S37233222(2)=Cij124(2,6)-Cij134(4,6)
       S37233222(3)=Cij124(2,6)-Cij134(5,6)
       S37233311(1)=Cij124(5,6)-Cij134(5,6)
       S37233311(2)=Cij124(6,6)-Cij134(5,6)
       S37233311(3)=Cij124(6,6)-Cij134(6,6)
       S37233321(1)=Cij124(6,6)-Cij134(5,6)
       S37233321(2)=Cij124(7,6)-Cij134(5,6)
       S37233321(3)=Cij124(7,6)-Cij134(6,6)
       S37233322(1)=Cij124(7,6)-Cij134(5,6)
       S37233322(2)=Cij124(2,6)-Cij134(5,6)
       S37233322(3)=Cij124(2,6)-Cij134(6,6)
       S37233331(1)=Cij124(6,6)-Cij134(6,6)
       S37233331(2)=Cij124(7,6)-Cij134(6,6)
       S37233331(3)=Cij124(7,6)-Cij134(7,6)
       S37233332(1)=Cij124(7,6)-Cij134(6,6)
       S37233332(2)=Cij124(2,6)-Cij134(6,6)
       S37233332(3)=Cij124(2,6)-Cij134(7,6)
       S37233333(1)=Cij124(7,6)-Cij134(7,6)
       S37233333(2)=Cij124(2,6)-Cij134(7,6)
       S37233333(3)=Cij124(2,6)-Cij134(2,6)
       S37311111(1)=Cij123(1,6)-Cij124(1,6)
       S37311111(2)=Cij123(3,6)-Cij124(3,6)
       S37311111(3)=-Cij124(3,6)
       S37321111(1)=Cij123(3,6)-Cij124(3,6)
       S37321111(2)=Cij123(4,6)-Cij124(4,6)
       S37321111(3)=-Cij124(4,6)
       S37322111(1)=Cij123(4,6)-Cij124(4,6)
       S37322111(2)=Cij123(5,6)-Cij124(5,6)
       S37322111(3)=-Cij124(5,6)
       S37322211(1)=Cij123(5,6)-Cij124(5,6)
       S37322211(2)=Cij123(6,6)-Cij124(6,6)
       S37322211(3)=-Cij124(6,6)
       S37322221(1)=Cij123(6,6)-Cij124(6,6)
       S37322221(2)=Cij123(7,6)-Cij124(7,6)
       S37322221(3)=-Cij124(7,6)
       S37322222(1)=Cij123(7,6)-Cij124(7,6)
       S37322222(2)=Cij123(2,6)-Cij124(2,6)
       S37322222(3)=-Cij124(2,6)
       S37331111(1)=-Cij124(3,6)
       S37331111(2)=-Cij124(4,6)
       S37331111(3)=-Cij124(4,6)
       S37332111(1)=-Cij124(4,6)
       S37332111(2)=-Cij124(5,6)
       S37332111(3)=-Cij124(5,6)
       S37332211(1)=-Cij124(5,6)
       S37332211(2)=-Cij124(6,6)
       S37332211(3)=-Cij124(6,6)
       S37332221(1)=-Cij124(6,6)
       S37332221(2)=-Cij124(7,6)
       S37332221(3)=-Cij124(7,6)
       S37332222(1)=-Cij124(7,6)
       S37332222(2)=-Cij124(2,6)
       S37332222(3)=-Cij124(2,6)
       S37333111(1)=-Cij124(4,6)
       S37333111(2)=-Cij124(5,6)
       S37333111(3)=-Cij124(5,6)
       S37333211(1)=-Cij124(5,6)
       S37333211(2)=-Cij124(6,6)
       S37333211(3)=-Cij124(6,6)
       S37333221(1)=-Cij124(6,6)
       S37333221(2)=-Cij124(7,6)
       S37333221(3)=-Cij124(7,6)
       S37333222(1)=-Cij124(7,6)
       S37333222(2)=-Cij124(2,6)
       S37333222(3)=-Cij124(2,6)
       S37333311(1)=-Cij124(5,6)
       S37333311(2)=-Cij124(6,6)
       S37333311(3)=-Cij124(6,6)
       S37333321(1)=-Cij124(6,6)
       S37333321(2)=-Cij124(7,6)
       S37333321(3)=-Cij124(7,6)
       S37333322(1)=-Cij124(7,6)
       S37333322(2)=-Cij124(2,6)
       S37333322(3)=-Cij124(2,6)
       S37333331(1)=-Cij124(6,6)
       S37333331(2)=-Cij124(7,6)
       S37333331(3)=-Cij124(7,6)
       S37333332(1)=-Cij124(7,6)
       S37333332(2)=-Cij124(2,6)
       S37333332(3)=-Cij124(2,6)
       S37333333(1)=-Cij124(7,6)
       S37333333(2)=-Cij124(2,6)
       S37333333(3)=-Cij124(2,6)
       aux00511(1,1,1)=-(F(1)*S3611111(1))-F(2)*S3621111(1)-F(3)*S363111
     &  1(1)+S37111111(k)*Z(1,l)+S37211111(k)*Z(2,l)+S37311111(k)*Z(3,l)
     &  +(S30051111(1)-S37111111(1)-S37221111(1)-S37331111(1))*Z(k,l)-10
     &  *S3h0051111(1)*ZZ(k,1,l,1)-10*S3h0052111(1)*ZZ(k,1,l,2)-10*S3h00
     &  53111(1)*ZZ(k,1,l,3)
       aux00521(1,1,1)=-(F(1)*S3612111(1))-F(2)*S3622111(1)-F(3)*S363211
     &  1(1)+S37121111(k)*Z(1,l)+S37221111(k)*Z(2,l)+S37321111(k)*Z(3,l)
     &  +(S30052111(1)-S37121111(1)-S37222111(1)-S37332111(1))*Z(k,l)-8*
     &  S3h0051211(1)*ZZ(k,1,l,1)-8*S3h0052211(1)*ZZ(k,1,l,2)-8*S3h00532
     &  11(1)*ZZ(k,1,l,3)-2*S3h0051111(1)*ZZ(k,2,l,1)-2*S3h0052111(1)*ZZ
     &  (k,2,l,2)-2*S3h0053111(1)*ZZ(k,2,l,3)
       aux00522(1,1,1)=-(F(1)*S3612211(1))-F(2)*S3622211(1)-F(3)*S363221
     &  1(1)+S37122111(k)*Z(1,l)+S37222111(k)*Z(2,l)+S37322111(k)*Z(3,l)
     &  +(S30052211(1)-S37122111(1)-S37222211(1)-S37332211(1))*Z(k,l)-6*
     &  S3h0051221(1)*ZZ(k,1,l,1)-6*S3h0052221(1)*ZZ(k,1,l,2)-6*S3h00532
     &  21(1)*ZZ(k,1,l,3)-4*S3h0051211(1)*ZZ(k,2,l,1)-4*S3h0052211(1)*ZZ
     &  (k,2,l,2)-4*S3h0053211(1)*ZZ(k,2,l,3)
       aux00522(2,1,1)=-(F(1)*S3612221(1))-F(2)*S3622221(1)-F(3)*S363222
     &  1(1)+S37122211(k)*Z(1,l)+S37222211(k)*Z(2,l)+S37322211(k)*Z(3,l)
     &  +(S30052221(1)-S37122211(1)-S37222221(1)-S37332221(1))*Z(k,l)-4*
     &  S3h0051222(1)*ZZ(k,1,l,1)-4*S3h0052222(1)*ZZ(k,1,l,2)-4*S3h00532
     &  22(1)*ZZ(k,1,l,3)-6*S3h0051221(1)*ZZ(k,2,l,1)-6*S3h0052221(1)*ZZ
     &  (k,2,l,2)-6*S3h0053221(1)*ZZ(k,2,l,3)
       aux00522(2,2,1)=-(F(1)*S3612222(1))-F(2)*S3622222(1)-F(3)*S363222
     &  2(1)+S37122221(k)*Z(1,l)+S37222221(k)*Z(2,l)+S37322221(k)*Z(3,l)
     &  +(S30052222(1)-S37122221(1)-S37222222(1)-S37332222(1))*Z(k,l)-2*
     &  S3h0051222(2)*ZZ(k,1,l,1)-2*S3h0052222(2)*ZZ(k,1,l,2)-2*S3h00532
     &  22(2)*ZZ(k,1,l,3)-8*S3h0051222(1)*ZZ(k,2,l,1)-8*S3h0052222(1)*ZZ
     &  (k,2,l,2)-8*S3h0053222(1)*ZZ(k,2,l,3)
       aux00522(2,2,2)=-(F(1)*S3612222(2))-F(2)*S3622222(2)-F(3)*S363222
     &  2(2)+S37122222(k)*Z(1,l)+S37222222(k)*Z(2,l)+S37322222(k)*Z(3,l)
     &  +(S30052222(2)-S37122222(1)-S37222222(2)-S37332222(2))*Z(k,l)-10
     &  *S3h0051222(2)*ZZ(k,2,l,1)-10*S3h0052222(2)*ZZ(k,2,l,2)-10*S3h00
     &  53222(2)*ZZ(k,2,l,3)
       aux00531(1,1,1)=-(F(1)*S3613111(1))-F(2)*S3623111(1)-F(3)*S363311
     &  1(1)+S37131111(k)*Z(1,l)+S37231111(k)*Z(2,l)+S37331111(k)*Z(3,l)
     &  +(S30053111(1)-S37131111(1)-S37232111(1)-S37333111(1))*Z(k,l)-8*
     &  S3h0051311(1)*ZZ(k,1,l,1)-8*S3h0052311(1)*ZZ(k,1,l,2)-8*S3h00533
     &  11(1)*ZZ(k,1,l,3)-2*S3h0051111(1)*ZZ(k,3,l,1)-2*S3h0052111(1)*ZZ
     &  (k,3,l,2)-2*S3h0053111(1)*ZZ(k,3,l,3)
       aux00532(1,1,1)=-(F(1)*S3613211(1))-F(2)*S3623211(1)-F(3)*S363321
     &  1(1)+S37132111(k)*Z(1,l)+S37232111(k)*Z(2,l)+S37332111(k)*Z(3,l)
     &  +(S30053211(1)-S37132111(1)-S37232211(1)-S37333211(1))*Z(k,l)-6*
     &  S3h0051321(1)*ZZ(k,1,l,1)-6*S3h0052321(1)*ZZ(k,1,l,2)-6*S3h00533
     &  21(1)*ZZ(k,1,l,3)-2*S3h0051311(1)*ZZ(k,2,l,1)-2*S3h0052311(1)*ZZ
     &  (k,2,l,2)-2*S3h0053311(1)*ZZ(k,2,l,3)-2*S3h0051211(1)*ZZ(k,3,l,1
     &  )-2*S3h0052211(1)*ZZ(k,3,l,2)-2*S3h0053211(1)*ZZ(k,3,l,3)
       aux00532(2,1,1)=-(F(1)*S3613221(1))-F(2)*S3623221(1)-F(3)*S363322
     &  1(1)+S37132211(k)*Z(1,l)+S37232211(k)*Z(2,l)+S37332211(k)*Z(3,l)
     &  +(S30053221(1)-S37132211(1)-S37232221(1)-S37333221(1))*Z(k,l)-4*
     &  S3h0051322(1)*ZZ(k,1,l,1)-4*S3h0052322(1)*ZZ(k,1,l,2)-4*S3h00533
     &  22(1)*ZZ(k,1,l,3)-4*S3h0051321(1)*ZZ(k,2,l,1)-4*S3h0052321(1)*ZZ
     &  (k,2,l,2)-4*S3h0053321(1)*ZZ(k,2,l,3)-2*S3h0051221(1)*ZZ(k,3,l,1
     &  )-2*S3h0052221(1)*ZZ(k,3,l,2)-2*S3h0053221(1)*ZZ(k,3,l,3)
       aux00532(2,2,1)=-(F(1)*S3613222(1))-F(2)*S3623222(1)-F(3)*S363322
     &  2(1)+S37132221(k)*Z(1,l)+S37232221(k)*Z(2,l)+S37332221(k)*Z(3,l)
     &  +(S30053222(1)-S37132221(1)-S37232222(1)-S37333222(1))*Z(k,l)-2*
     &  S3h0051322(2)*ZZ(k,1,l,1)-2*S3h0052322(2)*ZZ(k,1,l,2)-2*S3h00533
     &  22(2)*ZZ(k,1,l,3)-6*S3h0051322(1)*ZZ(k,2,l,1)-6*S3h0052322(1)*ZZ
     &  (k,2,l,2)-6*S3h0053322(1)*ZZ(k,2,l,3)-2*S3h0051222(1)*ZZ(k,3,l,1
     &  )-2*S3h0052222(1)*ZZ(k,3,l,2)-2*S3h0053222(1)*ZZ(k,3,l,3)
       aux00532(2,2,2)=-(F(1)*S3613222(2))-F(2)*S3623222(2)-F(3)*S363322
     &  2(2)+S37132222(k)*Z(1,l)+S37232222(k)*Z(2,l)+S37332222(k)*Z(3,l)
     &  +(S30053222(2)-S37132222(1)-S37232222(2)-S37333222(2))*Z(k,l)-8*
     &  S3h0051322(2)*ZZ(k,2,l,1)-8*S3h0052322(2)*ZZ(k,2,l,2)-8*S3h00533
     &  22(2)*ZZ(k,2,l,3)-2*S3h0051222(2)*ZZ(k,3,l,1)-2*S3h0052222(2)*ZZ
     &  (k,3,l,2)-2*S3h0053222(2)*ZZ(k,3,l,3)
       aux00533(1,1,1)=-(F(1)*S3613311(1))-F(2)*S3623311(1)-F(3)*S363331
     &  1(1)+S37133111(k)*Z(1,l)+S37233111(k)*Z(2,l)+S37333111(k)*Z(3,l)
     &  +(S30053311(1)-S37133111(1)-S37233211(1)-S37333311(1))*Z(k,l)-6*
     &  S3h0051331(1)*ZZ(k,1,l,1)-6*S3h0052331(1)*ZZ(k,1,l,2)-6*S3h00533
     &  31(1)*ZZ(k,1,l,3)-4*S3h0051311(1)*ZZ(k,3,l,1)-4*S3h0052311(1)*ZZ
     &  (k,3,l,2)-4*S3h0053311(1)*ZZ(k,3,l,3)
       aux00533(2,1,1)=-(F(1)*S3613321(1))-F(2)*S3623321(1)-F(3)*S363332
     &  1(1)+S37133211(k)*Z(1,l)+S37233211(k)*Z(2,l)+S37333211(k)*Z(3,l)
     &  +(S30053321(1)-S37133211(1)-S37233221(1)-S37333321(1))*Z(k,l)-4*
     &  S3h0051332(1)*ZZ(k,1,l,1)-4*S3h0052332(1)*ZZ(k,1,l,2)-4*S3h00533
     &  32(1)*ZZ(k,1,l,3)-2*S3h0051331(1)*ZZ(k,2,l,1)-2*S3h0052331(1)*ZZ
     &  (k,2,l,2)-2*S3h0053331(1)*ZZ(k,2,l,3)-4*S3h0051321(1)*ZZ(k,3,l,1
     &  )-4*S3h0052321(1)*ZZ(k,3,l,2)-4*S3h0053321(1)*ZZ(k,3,l,3)
       aux00533(2,2,1)=-(F(1)*S3613322(1))-F(2)*S3623322(1)-F(3)*S363332
     &  2(1)+S37133221(k)*Z(1,l)+S37233221(k)*Z(2,l)+S37333221(k)*Z(3,l)
     &  +(S30053322(1)-S37133221(1)-S37233222(1)-S37333322(1))*Z(k,l)-2*
     &  S3h0051332(2)*ZZ(k,1,l,1)-2*S3h0052332(2)*ZZ(k,1,l,2)-2*S3h00533
     &  32(2)*ZZ(k,1,l,3)-4*S3h0051332(1)*ZZ(k,2,l,1)-4*S3h0052332(1)*ZZ
     &  (k,2,l,2)-4*S3h0053332(1)*ZZ(k,2,l,3)-4*S3h0051322(1)*ZZ(k,3,l,1
     &  )-4*S3h0052322(1)*ZZ(k,3,l,2)-4*S3h0053322(1)*ZZ(k,3,l,3)
       aux00533(2,2,2)=-(F(1)*S3613322(2))-F(2)*S3623322(2)-F(3)*S363332
     &  2(2)+S37133222(k)*Z(1,l)+S37233222(k)*Z(2,l)+S37333222(k)*Z(3,l)
     &  +(S30053322(2)-S37133222(1)-S37233222(2)-S37333322(2))*Z(k,l)-6*
     &  S3h0051332(2)*ZZ(k,2,l,1)-6*S3h0052332(2)*ZZ(k,2,l,2)-6*S3h00533
     &  32(2)*ZZ(k,2,l,3)-4*S3h0051322(2)*ZZ(k,3,l,1)-4*S3h0052322(2)*ZZ
     &  (k,3,l,2)-4*S3h0053322(2)*ZZ(k,3,l,3)
       aux00533(3,1,1)=-(F(1)*S3613331(1))-F(2)*S3623331(1)-F(3)*S363333
     &  1(1)+S37133311(k)*Z(1,l)+S37233311(k)*Z(2,l)+S37333311(k)*Z(3,l)
     &  +(S30053331(1)-S37133311(1)-S37233321(1)-S37333331(1))*Z(k,l)-4*
     &  S3h0051333(1)*ZZ(k,1,l,1)-4*S3h0052333(1)*ZZ(k,1,l,2)-4*S3h00533
     &  33(1)*ZZ(k,1,l,3)-6*S3h0051331(1)*ZZ(k,3,l,1)-6*S3h0052331(1)*ZZ
     &  (k,3,l,2)-6*S3h0053331(1)*ZZ(k,3,l,3)
       aux00533(3,2,1)=-(F(1)*S3613332(1))-F(2)*S3623332(1)-F(3)*S363333
     &  2(1)+S37133321(k)*Z(1,l)+S37233321(k)*Z(2,l)+S37333321(k)*Z(3,l)
     &  +(S30053332(1)-S37133321(1)-S37233322(1)-S37333332(1))*Z(k,l)-2*
     &  S3h0051333(2)*ZZ(k,1,l,1)-2*S3h0052333(2)*ZZ(k,1,l,2)-2*S3h00533
     &  33(2)*ZZ(k,1,l,3)-2*S3h0051333(1)*ZZ(k,2,l,1)-2*S3h0052333(1)*ZZ
     &  (k,2,l,2)-2*S3h0053333(1)*ZZ(k,2,l,3)-6*S3h0051332(1)*ZZ(k,3,l,1
     &  )-6*S3h0052332(1)*ZZ(k,3,l,2)-6*S3h0053332(1)*ZZ(k,3,l,3)
       aux00533(3,2,2)=-(F(1)*S3613332(2))-F(2)*S3623332(2)-F(3)*S363333
     &  2(2)+S37133322(k)*Z(1,l)+S37233322(k)*Z(2,l)+S37333322(k)*Z(3,l)
     &  +(S30053332(2)-S37133322(1)-S37233322(2)-S37333332(2))*Z(k,l)-4*
     &  S3h0051333(2)*ZZ(k,2,l,1)-4*S3h0052333(2)*ZZ(k,2,l,2)-4*S3h00533
     &  33(2)*ZZ(k,2,l,3)-6*S3h0051332(2)*ZZ(k,3,l,1)-6*S3h0052332(2)*ZZ
     &  (k,3,l,2)-6*S3h0053332(2)*ZZ(k,3,l,3)
       aux00533(3,3,1)=-(F(1)*S3613333(1))-F(2)*S3623333(1)-F(3)*S363333
     &  3(1)+S37133331(k)*Z(1,l)+S37233331(k)*Z(2,l)+S37333331(k)*Z(3,l)
     &  +(S30053333(1)-S37133331(1)-S37233332(1)-S37333333(1))*Z(k,l)-2*
     &  S3h0051333(3)*ZZ(k,1,l,1)-2*S3h0052333(3)*ZZ(k,1,l,2)-2*S3h00533
     &  33(3)*ZZ(k,1,l,3)-8*S3h0051333(1)*ZZ(k,3,l,1)-8*S3h0052333(1)*ZZ
     &  (k,3,l,2)-8*S3h0053333(1)*ZZ(k,3,l,3)
       aux00533(3,3,2)=-(F(1)*S3613333(2))-F(2)*S3623333(2)-F(3)*S363333
     &  3(2)+S37133332(k)*Z(1,l)+S37233332(k)*Z(2,l)+S37333332(k)*Z(3,l)
     &  +(S30053333(2)-S37133332(1)-S37233332(2)-S37333333(2))*Z(k,l)-2*
     &  S3h0051333(3)*ZZ(k,2,l,1)-2*S3h0052333(3)*ZZ(k,2,l,2)-2*S3h00533
     &  33(3)*ZZ(k,2,l,3)-8*S3h0051333(2)*ZZ(k,3,l,1)-8*S3h0052333(2)*ZZ
     &  (k,3,l,2)-8*S3h0053333(2)*ZZ(k,3,l,3)
       aux00533(3,3,3)=-(F(1)*S3613333(3))-F(2)*S3623333(3)-F(3)*S363333
     &  3(3)+S37133333(k)*Z(1,l)+S37233333(k)*Z(2,l)+S37333333(k)*Z(3,l)
     &  +(S30053333(3)-S37133333(1)-S37233333(2)-S37333333(3))*Z(k,l)-10
     &  *S3h0051333(3)*ZZ(k,3,l,1)-10*S3h0052333(3)*ZZ(k,3,l,2)-10*S3h00
     &  53333(3)*ZZ(k,3,l,3)
       temp00511(1,1,1)=I24Z*(aux00511(1,1,1)+10*F(4)*temp0041(1,1,1)+F(
     &  5)*temp511(1,1,1)+80*temp00003(1,1,1)*ZZ(k,1,l,1))
       temp00521(1,1,1)=I24Z*(aux00521(1,1,1)+2*F(6)*temp0041(1,1,1)+8*F
     &  (4)*temp0042(1,1,1)+F(5)*temp521(1,1,1)+48*temp00003(2,1,1)*ZZ(k
     &  ,1,l,1)+16*temp00003(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp00522(1,1,1)=I24Z*(aux00522(1,1,1)+4*F(6)*temp0042(1,1,1)+6*F
     &  (4)*temp0042(2,1,1)+F(5)*temp522(1,1,1)+24*(temp00003(2,2,1)*ZZ(
     &  k,1,l,1)+temp00003(2,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))+8*temp00003
     &  (1,1,1)*ZZ(k,2,l,2))
       temp00522(2,1,1)=I24Z*(aux00522(2,1,1)+6*F(6)*temp0042(2,1,1)+4*F
     &  (4)*temp0042(2,2,1)+F(5)*temp522(2,1,1)+8*temp00003(2,2,2)*ZZ(k,
     &  1,l,1)+24*(temp00003(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00003(
     &  2,1,1)*ZZ(k,2,l,2)))
       temp00522(2,2,1)=I24Z*(aux00522(2,2,1)+8*F(6)*temp0042(2,2,1)+2*F
     &  (4)*temp0042(2,2,2)+F(5)*temp522(2,2,1)+16*temp00003(2,2,2)*(ZZ(
     &  k,1,l,2)+ZZ(k,2,l,1))+48*temp00003(2,2,1)*ZZ(k,2,l,2))
       temp00522(2,2,2)=I24Z*(aux00522(2,2,2)+10*F(6)*temp0042(2,2,2)+F(
     &  5)*temp522(2,2,2)+80*temp00003(2,2,2)*ZZ(k,2,l,2))
       temp00531(1,1,1)=I24Z*(aux00531(1,1,1)+2*F(7)*temp0041(1,1,1)+8*F
     &  (4)*temp0043(1,1,1)+F(5)*temp531(1,1,1)+48*temp00003(3,1,1)*ZZ(k
     &  ,1,l,1)+16*temp00003(1,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)))
       temp00532(1,1,1)=I24Z*(aux00532(1,1,1)+2*F(7)*temp0042(1,1,1)+2*F
     &  (6)*temp0043(1,1,1)+6*F(4)*temp0043(2,1,1)+F(5)*temp532(1,1,1)+2
     &  4*temp00003(3,2,1)*ZZ(k,1,l,1)+12*(temp00003(3,1,1)*(ZZ(k,1,l,2)
     &  +ZZ(k,2,l,1))+temp00003(2,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)))+4*temp
     &  00003(1,1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp00532(2,1,1)=I24Z*(aux00532(2,1,1)+2*F(7)*temp0042(2,1,1)+4*F
     &  (6)*temp0043(2,1,1)+4*F(4)*temp0043(2,2,1)+F(5)*temp532(2,1,1)+1
     &  6*temp00003(3,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp00003(3,2,2)
     &  *ZZ(k,1,l,1)+temp00003(3,1,1)*ZZ(k,2,l,2))+8*temp00003(2,2,1)*(Z
     &  Z(k,1,l,3)+ZZ(k,3,l,1))+8*temp00003(2,1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l
     &  ,2)))
       temp00532(2,2,1)=I24Z*(aux00532(2,2,1)+2*F(7)*temp0042(2,2,1)+6*F
     &  (6)*temp0043(2,2,1)+2*F(4)*temp0043(2,2,2)+F(5)*temp532(2,2,1)+2
     &  4*temp00003(3,2,1)*ZZ(k,2,l,2)+4*temp00003(2,2,2)*(ZZ(k,1,l,3)+Z
     &  Z(k,3,l,1))+12*(temp00003(3,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0
     &  0003(2,2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))))
       temp00532(2,2,2)=I24Z*(aux00532(2,2,2)+2*F(7)*temp0042(2,2,2)+8*F
     &  (6)*temp0043(2,2,2)+F(5)*temp532(2,2,2)+48*temp00003(3,2,2)*ZZ(k
     &  ,2,l,2)+16*temp00003(2,2,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp00533(1,1,1)=I24Z*(aux00533(1,1,1)+4*F(7)*temp0043(1,1,1)+6*F
     &  (4)*temp0043(3,1,1)+F(5)*temp533(1,1,1)+24*(temp00003(3,3,1)*ZZ(
     &  k,1,l,1)+temp00003(3,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)))+8*temp00003
     &  (1,1,1)*ZZ(k,3,l,3))
       temp00533(2,1,1)=I24Z*(aux00533(2,1,1)+4*F(7)*temp0043(2,1,1)+2*F
     &  (6)*temp0043(3,1,1)+4*F(4)*temp0043(3,2,1)+F(5)*temp533(2,1,1)+8
     &  *(temp00003(3,3,2)*ZZ(k,1,l,1)+temp00003(3,3,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))+16*temp00003(3,2,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*temp0
     &  0003(3,1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*temp00003(2,1,1)*ZZ(k,3,
     &  l,3))
       temp00533(2,2,1)=I24Z*(aux00533(2,2,1)+4*F(7)*temp0043(2,2,1)+4*F
     &  (6)*temp0043(3,2,1)+2*F(4)*temp0043(3,2,2)+F(5)*temp533(2,2,1)+8
     &  *(temp00003(3,3,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00003(3,3,1)*ZZ
     &  (k,2,l,2))+8*temp00003(3,2,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+16*temp0
     &  0003(3,2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*temp00003(2,2,1)*ZZ(k,3,
     &  l,3))
       temp00533(2,2,2)=I24Z*(aux00533(2,2,2)+4*F(7)*temp0043(2,2,2)+6*F
     &  (6)*temp0043(3,2,2)+F(5)*temp533(2,2,2)+24*(temp00003(3,3,2)*ZZ(
     &  k,2,l,2)+temp00003(3,2,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))+8*temp00003
     &  (2,2,2)*ZZ(k,3,l,3))
       temp00533(3,1,1)=I24Z*(aux00533(3,1,1)+6*F(7)*temp0043(3,1,1)+4*F
     &  (4)*temp0043(3,3,1)+F(5)*temp533(3,1,1)+8*temp00003(3,3,3)*ZZ(k,
     &  1,l,1)+24*(temp00003(3,3,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+temp00003(
     &  3,1,1)*ZZ(k,3,l,3)))
       temp00533(3,2,1)=I24Z*(aux00533(3,2,1)+6*F(7)*temp0043(3,2,1)+2*F
     &  (6)*temp0043(3,3,1)+2*F(4)*temp0043(3,3,2)+F(5)*temp533(3,2,1)+4
     &  *temp00003(3,3,3)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*(temp00003(3,3,2)
     &  *(ZZ(k,1,l,3)+ZZ(k,3,l,1))+temp00003(3,3,1)*(ZZ(k,2,l,3)+ZZ(k,3,
     &  l,2)))+24*temp00003(3,2,1)*ZZ(k,3,l,3))
       temp00533(3,2,2)=I24Z*(aux00533(3,2,2)+6*F(7)*temp0043(3,2,2)+4*F
     &  (6)*temp0043(3,3,2)+F(5)*temp533(3,2,2)+8*temp00003(3,3,3)*ZZ(k,
     &  2,l,2)+24*(temp00003(3,3,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+temp00003(
     &  3,2,2)*ZZ(k,3,l,3)))
       temp00533(3,3,1)=I24Z*(aux00533(3,3,1)+8*F(7)*temp0043(3,3,1)+2*F
     &  (4)*temp0043(3,3,3)+F(5)*temp533(3,3,1)+16*temp00003(3,3,3)*(ZZ(
     &  k,1,l,3)+ZZ(k,3,l,1))+48*temp00003(3,3,1)*ZZ(k,3,l,3))
       temp00533(3,3,2)=I24Z*(aux00533(3,3,2)+8*F(7)*temp0043(3,3,2)+2*F
     &  (6)*temp0043(3,3,3)+F(5)*temp533(3,3,2)+16*temp00003(3,3,3)*(ZZ(
     &  k,2,l,3)+ZZ(k,3,l,2))+48*temp00003(3,3,2)*ZZ(k,3,l,3))
       temp00533(3,3,3)=I24Z*(aux00533(3,3,3)+10*F(7)*temp0043(3,3,3)+F(
     &  5)*temp533(3,3,3)+80*temp00003(3,3,3)*ZZ(k,3,l,3))
       temp00511(1,1,2)=temp00521(1,1,1)
       temp00511(1,1,3)=temp00531(1,1,1)
       temp00511(1,2,1)=temp00521(1,1,1)
       temp00511(1,2,2)=temp00522(1,1,1)
       temp00511(1,2,3)=temp00532(1,1,1)
       temp00511(1,3,1)=temp00531(1,1,1)
       temp00511(1,3,2)=temp00532(1,1,1)
       temp00511(1,3,3)=temp00533(1,1,1)
       temp00521(1,1,2)=temp00522(1,1,1)
       temp00521(1,1,3)=temp00532(1,1,1)
       temp00521(1,2,1)=temp00522(1,1,1)
       temp00521(1,2,2)=temp00522(2,1,1)
       temp00521(1,2,3)=temp00532(2,1,1)
       temp00521(1,3,1)=temp00532(1,1,1)
       temp00521(1,3,2)=temp00532(2,1,1)
       temp00521(1,3,3)=temp00533(2,1,1)
       temp00522(1,1,2)=temp00522(2,1,1)
       temp00522(1,1,3)=temp00532(2,1,1)
       temp00522(1,2,1)=temp00522(2,1,1)
       temp00522(1,2,2)=temp00522(2,2,1)
       temp00522(1,2,3)=temp00532(2,2,1)
       temp00522(1,3,1)=temp00532(2,1,1)
       temp00522(1,3,2)=temp00532(2,2,1)
       temp00522(1,3,3)=temp00533(2,2,1)
       temp00522(2,1,2)=temp00522(2,2,1)
       temp00522(2,1,3)=temp00532(2,2,1)
       temp00522(2,2,3)=temp00532(2,2,2)
       temp00522(2,3,1)=temp00532(2,2,1)
       temp00522(2,3,2)=temp00532(2,2,2)
       temp00522(2,3,3)=temp00533(2,2,2)
       temp00531(1,1,2)=temp00532(1,1,1)
       temp00531(1,1,3)=temp00533(1,1,1)
       temp00531(1,2,1)=temp00532(1,1,1)
       temp00531(1,2,2)=temp00532(2,1,1)
       temp00531(1,2,3)=temp00533(2,1,1)
       temp00531(1,3,1)=temp00533(1,1,1)
       temp00531(1,3,2)=temp00533(2,1,1)
       temp00531(1,3,3)=temp00533(3,1,1)
       temp00532(1,1,2)=temp00532(2,1,1)
       temp00532(1,1,3)=temp00533(2,1,1)
       temp00532(1,2,1)=temp00532(2,1,1)
       temp00532(1,2,2)=temp00532(2,2,1)
       temp00532(1,2,3)=temp00533(2,2,1)
       temp00532(1,3,1)=temp00533(2,1,1)
       temp00532(1,3,2)=temp00533(2,2,1)
       temp00532(1,3,3)=temp00533(3,2,1)
       temp00532(2,1,2)=temp00532(2,2,1)
       temp00532(2,1,3)=temp00533(2,2,1)
       temp00532(2,2,3)=temp00533(2,2,2)
       temp00532(2,3,1)=temp00533(2,2,1)
       temp00532(2,3,2)=temp00533(2,2,2)
       temp00532(2,3,3)=temp00533(3,2,2)
       temp00533(1,1,2)=temp00533(2,1,1)
       temp00533(1,1,3)=temp00533(3,1,1)
       temp00533(1,2,1)=temp00533(2,1,1)
       temp00533(1,2,2)=temp00533(2,2,1)
       temp00533(1,2,3)=temp00533(3,2,1)
       temp00533(1,3,1)=temp00533(3,1,1)
       temp00533(1,3,2)=temp00533(3,2,1)
       temp00533(1,3,3)=temp00533(3,3,1)
       temp00533(2,1,2)=temp00533(2,2,1)
       temp00533(2,1,3)=temp00533(3,2,1)
       temp00533(2,2,3)=temp00533(3,2,2)
       temp00533(2,3,1)=temp00533(3,2,1)
       temp00533(2,3,2)=temp00533(3,2,2)
       temp00533(2,3,3)=temp00533(3,3,2)
       temp00533(3,1,2)=temp00533(3,2,1)
       temp00533(3,1,3)=temp00533(3,3,1)
       temp00533(3,2,3)=temp00533(3,3,2)
       aux6111(1,1,1)=-(S37111111(1)*Z(jj,1))-S37211111(1)*Z(jj,2)-S3731
     &  1111(1)*Z(jj,3)
       aux6211(1,1,1)=-(S37121111(1)*Z(jj,1))-S37221111(1)*Z(jj,2)-S3732
     &  1111(1)*Z(jj,3)
       aux6221(1,1,1)=-(S37122111(1)*Z(jj,1))-S37222111(1)*Z(jj,2)-S3732
     &  2111(1)*Z(jj,3)
       aux6222(1,1,1)=-(S37122211(1)*Z(jj,1))-S37222211(1)*Z(jj,2)-S3732
     &  2211(1)*Z(jj,3)
       aux6222(2,1,1)=-(S37122221(1)*Z(jj,1))-S37222221(1)*Z(jj,2)-S3732
     &  2221(1)*Z(jj,3)
       aux6222(2,2,1)=-(S37122222(1)*Z(jj,1))-S37222222(1)*Z(jj,2)-S3732
     &  2222(1)*Z(jj,3)
       aux6222(2,2,2)=-(S37122222(2)*Z(jj,1))-S37222222(2)*Z(jj,2)-S3732
     &  2222(2)*Z(jj,3)
       aux6311(1,1,1)=-(S37131111(1)*Z(jj,1))-S37231111(1)*Z(jj,2)-S3733
     &  1111(1)*Z(jj,3)
       aux6321(1,1,1)=-(S37132111(1)*Z(jj,1))-S37232111(1)*Z(jj,2)-S3733
     &  2111(1)*Z(jj,3)
       aux6322(1,1,1)=-(S37132211(1)*Z(jj,1))-S37232211(1)*Z(jj,2)-S3733
     &  2211(1)*Z(jj,3)
       aux6322(2,1,1)=-(S37132221(1)*Z(jj,1))-S37232221(1)*Z(jj,2)-S3733
     &  2221(1)*Z(jj,3)
       aux6322(2,2,1)=-(S37132222(1)*Z(jj,1))-S37232222(1)*Z(jj,2)-S3733
     &  2222(1)*Z(jj,3)
       aux6322(2,2,2)=-(S37132222(2)*Z(jj,1))-S37232222(2)*Z(jj,2)-S3733
     &  2222(2)*Z(jj,3)
       aux6331(1,1,1)=-(S37133111(1)*Z(jj,1))-S37233111(1)*Z(jj,2)-S3733
     &  3111(1)*Z(jj,3)
       aux6332(1,1,1)=-(S37133211(1)*Z(jj,1))-S37233211(1)*Z(jj,2)-S3733
     &  3211(1)*Z(jj,3)
       aux6332(2,1,1)=-(S37133221(1)*Z(jj,1))-S37233221(1)*Z(jj,2)-S3733
     &  3221(1)*Z(jj,3)
       aux6332(2,2,1)=-(S37133222(1)*Z(jj,1))-S37233222(1)*Z(jj,2)-S3733
     &  3222(1)*Z(jj,3)
       aux6332(2,2,2)=-(S37133222(2)*Z(jj,1))-S37233222(2)*Z(jj,2)-S3733
     &  3222(2)*Z(jj,3)
       aux6333(1,1,1)=-(S37133311(1)*Z(jj,1))-S37233311(1)*Z(jj,2)-S3733
     &  3311(1)*Z(jj,3)
       aux6333(2,1,1)=-(S37133321(1)*Z(jj,1))-S37233321(1)*Z(jj,2)-S3733
     &  3321(1)*Z(jj,3)
       aux6333(2,2,1)=-(S37133322(1)*Z(jj,1))-S37233322(1)*Z(jj,2)-S3733
     &  3322(1)*Z(jj,3)
       aux6333(2,2,2)=-(S37133322(2)*Z(jj,1))-S37233322(2)*Z(jj,2)-S3733
     &  3322(2)*Z(jj,3)
       aux6333(3,1,1)=-(S37133331(1)*Z(jj,1))-S37233331(1)*Z(jj,2)-S3733
     &  3331(1)*Z(jj,3)
       aux6333(3,2,1)=-(S37133332(1)*Z(jj,1))-S37233332(1)*Z(jj,2)-S3733
     &  3332(1)*Z(jj,3)
       aux6333(3,2,2)=-(S37133332(2)*Z(jj,1))-S37233332(2)*Z(jj,2)-S3733
     &  3332(2)*Z(jj,3)
       aux6333(3,3,1)=-(S37133333(1)*Z(jj,1))-S37233333(1)*Z(jj,2)-S3733
     &  3333(1)*Z(jj,3)
       aux6333(3,3,2)=-(S37133333(2)*Z(jj,1))-S37233333(2)*Z(jj,2)-S3733
     &  3333(2)*Z(jj,3)
       aux6333(3,3,3)=-(S37133333(3)*Z(jj,1))-S37233333(3)*Z(jj,2)-S3733
     &  3333(3)*Z(jj,3)
       temp6111(1,1,1)=IX*(aux6111(1,1,1)+12*temp00511(1,1,1)*Z(jj,1))
       temp6211(1,1,1)=IX*(aux6211(1,1,1)+10*temp00521(1,1,1)*Z(jj,1)+2*
     &  temp00511(1,1,1)*Z(jj,2))
       temp6221(1,1,1)=IX*(aux6221(1,1,1)+8*temp00522(1,1,1)*Z(jj,1)+4*t
     &  emp00521(1,1,1)*Z(jj,2))
       temp6222(1,1,1)=IX*(aux6222(1,1,1)+6*(temp00522(2,1,1)*Z(jj,1)+te
     &  mp00522(1,1,1)*Z(jj,2)))
       temp6222(2,1,1)=IX*(aux6222(2,1,1)+4*temp00522(2,2,1)*Z(jj,1)+8*t
     &  emp00522(2,1,1)*Z(jj,2))
       temp6222(2,2,1)=IX*(aux6222(2,2,1)+2*temp00522(2,2,2)*Z(jj,1)+10*
     &  temp00522(2,2,1)*Z(jj,2))
       temp6222(2,2,2)=IX*(aux6222(2,2,2)+12*temp00522(2,2,2)*Z(jj,2))
       temp6311(1,1,1)=IX*(aux6311(1,1,1)+10*temp00531(1,1,1)*Z(jj,1)+2*
     &  temp00511(1,1,1)*Z(jj,3))
       temp6321(1,1,1)=IX*(aux6321(1,1,1)+8*temp00532(1,1,1)*Z(jj,1)+2*(
     &  temp00531(1,1,1)*Z(jj,2)+temp00521(1,1,1)*Z(jj,3)))
       temp6322(1,1,1)=IX*(aux6322(1,1,1)+6*temp00532(2,1,1)*Z(jj,1)+4*t
     &  emp00532(1,1,1)*Z(jj,2)+2*temp00522(1,1,1)*Z(jj,3))
       temp6322(2,1,1)=IX*(aux6322(2,1,1)+4*temp00532(2,2,1)*Z(jj,1)+6*t
     &  emp00532(2,1,1)*Z(jj,2)+2*temp00522(2,1,1)*Z(jj,3))
       temp6322(2,2,1)=IX*(aux6322(2,2,1)+8*temp00532(2,2,1)*Z(jj,2)+2*(
     &  temp00532(2,2,2)*Z(jj,1)+temp00522(2,2,1)*Z(jj,3)))
       temp6322(2,2,2)=IX*(aux6322(2,2,2)+10*temp00532(2,2,2)*Z(jj,2)+2*
     &  temp00522(2,2,2)*Z(jj,3))
       temp6331(1,1,1)=IX*(aux6331(1,1,1)+8*temp00533(1,1,1)*Z(jj,1)+4*t
     &  emp00531(1,1,1)*Z(jj,3))
       temp6332(1,1,1)=IX*(aux6332(1,1,1)+6*temp00533(2,1,1)*Z(jj,1)+2*t
     &  emp00533(1,1,1)*Z(jj,2)+4*temp00532(1,1,1)*Z(jj,3))
       temp6332(2,1,1)=IX*(aux6332(2,1,1)+4*(temp00533(2,2,1)*Z(jj,1)+te
     &  mp00533(2,1,1)*Z(jj,2))+4*temp00532(2,1,1)*Z(jj,3))
       temp6332(2,2,1)=IX*(aux6332(2,2,1)+2*temp00533(2,2,2)*Z(jj,1)+6*t
     &  emp00533(2,2,1)*Z(jj,2)+4*temp00532(2,2,1)*Z(jj,3))
       temp6332(2,2,2)=IX*(aux6332(2,2,2)+8*temp00533(2,2,2)*Z(jj,2)+4*t
     &  emp00532(2,2,2)*Z(jj,3))
       temp6333(1,1,1)=IX*(aux6333(1,1,1)+6*(temp00533(3,1,1)*Z(jj,1)+te
     &  mp00533(1,1,1)*Z(jj,3)))
       temp6333(2,1,1)=IX*(aux6333(2,1,1)+4*temp00533(3,2,1)*Z(jj,1)+2*t
     &  emp00533(3,1,1)*Z(jj,2)+6*temp00533(2,1,1)*Z(jj,3))
       temp6333(2,2,1)=IX*(aux6333(2,2,1)+2*temp00533(3,2,2)*Z(jj,1)+4*t
     &  emp00533(3,2,1)*Z(jj,2)+6*temp00533(2,2,1)*Z(jj,3))
       temp6333(2,2,2)=IX*(aux6333(2,2,2)+6*(temp00533(3,2,2)*Z(jj,2)+te
     &  mp00533(2,2,2)*Z(jj,3)))
       temp6333(3,1,1)=IX*(aux6333(3,1,1)+4*temp00533(3,3,1)*Z(jj,1)+8*t
     &  emp00533(3,1,1)*Z(jj,3))
       temp6333(3,2,1)=IX*(aux6333(3,2,1)+2*(temp00533(3,3,2)*Z(jj,1)+te
     &  mp00533(3,3,1)*Z(jj,2))+8*temp00533(3,2,1)*Z(jj,3))
       temp6333(3,2,2)=IX*(aux6333(3,2,2)+4*temp00533(3,3,2)*Z(jj,2)+8*t
     &  emp00533(3,2,2)*Z(jj,3))
       temp6333(3,3,1)=IX*(aux6333(3,3,1)+2*temp00533(3,3,3)*Z(jj,1)+10*
     &  temp00533(3,3,1)*Z(jj,3))
       temp6333(3,3,2)=IX*(aux6333(3,3,2)+2*temp00533(3,3,3)*Z(jj,2)+10*
     &  temp00533(3,3,2)*Z(jj,3))
       temp6333(3,3,3)=IX*(aux6333(3,3,3)+12*temp00533(3,3,3)*Z(jj,3))
       temp6111(1,1,2)=temp6211(1,1,1)
       temp6111(1,1,3)=temp6311(1,1,1)
       temp6111(1,2,1)=temp6211(1,1,1)
       temp6111(1,2,2)=temp6221(1,1,1)
       temp6111(1,2,3)=temp6321(1,1,1)
       temp6111(1,3,1)=temp6311(1,1,1)
       temp6111(1,3,2)=temp6321(1,1,1)
       temp6111(1,3,3)=temp6331(1,1,1)
       temp6211(1,1,2)=temp6221(1,1,1)
       temp6211(1,1,3)=temp6321(1,1,1)
       temp6211(1,2,1)=temp6221(1,1,1)
       temp6211(1,2,2)=temp6222(1,1,1)
       temp6211(1,2,3)=temp6322(1,1,1)
       temp6211(1,3,1)=temp6321(1,1,1)
       temp6211(1,3,2)=temp6322(1,1,1)
       temp6211(1,3,3)=temp6332(1,1,1)
       temp6221(1,1,2)=temp6222(1,1,1)
       temp6221(1,1,3)=temp6322(1,1,1)
       temp6221(1,2,1)=temp6222(1,1,1)
       temp6221(1,2,2)=temp6222(2,1,1)
       temp6221(1,2,3)=temp6322(2,1,1)
       temp6221(1,3,1)=temp6322(1,1,1)
       temp6221(1,3,2)=temp6322(2,1,1)
       temp6221(1,3,3)=temp6332(2,1,1)
       temp6222(1,1,2)=temp6222(2,1,1)
       temp6222(1,1,3)=temp6322(2,1,1)
       temp6222(1,2,1)=temp6222(2,1,1)
       temp6222(1,2,2)=temp6222(2,2,1)
       temp6222(1,2,3)=temp6322(2,2,1)
       temp6222(1,3,1)=temp6322(2,1,1)
       temp6222(1,3,2)=temp6322(2,2,1)
       temp6222(1,3,3)=temp6332(2,2,1)
       temp6222(2,1,2)=temp6222(2,2,1)
       temp6222(2,1,3)=temp6322(2,2,1)
       temp6222(2,2,3)=temp6322(2,2,2)
       temp6222(2,3,1)=temp6322(2,2,1)
       temp6222(2,3,2)=temp6322(2,2,2)
       temp6222(2,3,3)=temp6332(2,2,2)
       temp6311(1,1,2)=temp6321(1,1,1)
       temp6311(1,1,3)=temp6331(1,1,1)
       temp6311(1,2,1)=temp6321(1,1,1)
       temp6311(1,2,2)=temp6322(1,1,1)
       temp6311(1,2,3)=temp6332(1,1,1)
       temp6311(1,3,1)=temp6331(1,1,1)
       temp6311(1,3,2)=temp6332(1,1,1)
       temp6311(1,3,3)=temp6333(1,1,1)
       temp6321(1,1,2)=temp6322(1,1,1)
       temp6321(1,1,3)=temp6332(1,1,1)
       temp6321(1,2,1)=temp6322(1,1,1)
       temp6321(1,2,2)=temp6322(2,1,1)
       temp6321(1,2,3)=temp6332(2,1,1)
       temp6321(1,3,1)=temp6332(1,1,1)
       temp6321(1,3,2)=temp6332(2,1,1)
       temp6321(1,3,3)=temp6333(2,1,1)
       temp6322(1,1,2)=temp6322(2,1,1)
       temp6322(1,1,3)=temp6332(2,1,1)
       temp6322(1,2,1)=temp6322(2,1,1)
       temp6322(1,2,2)=temp6322(2,2,1)
       temp6322(1,2,3)=temp6332(2,2,1)
       temp6322(1,3,1)=temp6332(2,1,1)
       temp6322(1,3,2)=temp6332(2,2,1)
       temp6322(1,3,3)=temp6333(2,2,1)
       temp6322(2,1,2)=temp6322(2,2,1)
       temp6322(2,1,3)=temp6332(2,2,1)
       temp6322(2,2,3)=temp6332(2,2,2)
       temp6322(2,3,1)=temp6332(2,2,1)
       temp6322(2,3,2)=temp6332(2,2,2)
       temp6322(2,3,3)=temp6333(2,2,2)
       temp6331(1,1,2)=temp6332(1,1,1)
       temp6331(1,1,3)=temp6333(1,1,1)
       temp6331(1,2,1)=temp6332(1,1,1)
       temp6331(1,2,2)=temp6332(2,1,1)
       temp6331(1,2,3)=temp6333(2,1,1)
       temp6331(1,3,1)=temp6333(1,1,1)
       temp6331(1,3,2)=temp6333(2,1,1)
       temp6331(1,3,3)=temp6333(3,1,1)
       temp6332(1,1,2)=temp6332(2,1,1)
       temp6332(1,1,3)=temp6333(2,1,1)
       temp6332(1,2,1)=temp6332(2,1,1)
       temp6332(1,2,2)=temp6332(2,2,1)
       temp6332(1,2,3)=temp6333(2,2,1)
       temp6332(1,3,1)=temp6333(2,1,1)
       temp6332(1,3,2)=temp6333(2,2,1)
       temp6332(1,3,3)=temp6333(3,2,1)
       temp6332(2,1,2)=temp6332(2,2,1)
       temp6332(2,1,3)=temp6333(2,2,1)
       temp6332(2,2,3)=temp6333(2,2,2)
       temp6332(2,3,1)=temp6333(2,2,1)
       temp6332(2,3,2)=temp6333(2,2,2)
       temp6332(2,3,3)=temp6333(3,2,2)
       temp6333(1,1,2)=temp6333(2,1,1)
       temp6333(1,1,3)=temp6333(3,1,1)
       temp6333(1,2,1)=temp6333(2,1,1)
       temp6333(1,2,2)=temp6333(2,2,1)
       temp6333(1,2,3)=temp6333(3,2,1)
       temp6333(1,3,1)=temp6333(3,1,1)
       temp6333(1,3,2)=temp6333(3,2,1)
       temp6333(1,3,3)=temp6333(3,3,1)
       temp6333(2,1,2)=temp6333(2,2,1)
       temp6333(2,1,3)=temp6333(3,2,1)
       temp6333(2,2,3)=temp6333(3,2,2)
       temp6333(2,3,1)=temp6333(3,2,1)
       temp6333(2,3,2)=temp6333(3,2,2)
       temp6333(2,3,3)=temp6333(3,3,2)
       temp6333(3,1,2)=temp6333(3,2,1)
       temp6333(3,1,3)=temp6333(3,3,1)
       temp6333(3,2,3)=temp6333(3,3,2)
c                Step2
       tempD4000000=I12Z*(auxD4000000+tempD40000*F(5)-det4*temp00002(k,l
     &  ))
       temp00002(1,1)=I16Z*(aux00002(1,1)+4*F(4)*temp00001(1)+F(5)*temp0
     &  02(1,1)-det4*temp0041(1,k,l)+8*tempD4000000*ZZ(k,1,l,1))
       temp00002(2,1)=I16Z*(aux00002(2,1)+2*(F(6)*temp00001(1)+F(4)*temp
     &  00001(2))+F(5)*temp002(2,1)-det4*temp0042(1,k,l)+4*tempD4000000*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp00002(2,2)=I16Z*(aux00002(2,2)+4*F(6)*temp00001(2)+F(5)*temp0
     &  02(2,2)-det4*temp0042(2,k,l)+8*tempD4000000*ZZ(k,2,l,2))
       temp00002(3,1)=I16Z*(aux00002(3,1)+2*(F(7)*temp00001(1)+F(4)*temp
     &  00001(3))+F(5)*temp002(3,1)-det4*temp0043(1,k,l)+4*tempD4000000*
     &  (ZZ(k,1,l,3)+ZZ(k,3,l,1)))
       temp00002(3,2)=I16Z*(aux00002(3,2)+2*(F(7)*temp00001(2)+F(6)*temp
     &  00001(3))+F(5)*temp002(3,2)-det4*temp0043(2,k,l)+4*tempD4000000*
     &  (ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp00002(3,3)=I16Z*(aux00002(3,3)+4*F(7)*temp00001(3)+F(5)*temp0
     &  02(3,3)-det4*temp0043(3,k,l)+8*tempD4000000*ZZ(k,3,l,3))
       temp00002(1,2)=temp00002(2,1)
       temp00002(1,3)=temp00002(3,1)
       temp00002(2,3)=temp00002(3,2)
       temp0041(1,1,1)=I20Z*(aux0041(1,1,1)+8*F(4)*temp003(1,1,1)+F(5)*t
     &  emp41(1,1,1)-det4*temp6111(1,k,l)+48*temp00002(1,1)*ZZ(k,1,l,1))
       temp0042(1,1,1)=I20Z*(aux0042(1,1,1)+2*F(6)*temp003(1,1,1)+6*F(4)
     &  *temp003(2,1,1)+F(5)*temp42(1,1,1)-det4*temp6211(1,k,l)+24*temp0
     &  0002(2,1)*ZZ(k,1,l,1)+12*temp00002(1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)
     &  ))
       temp0042(2,1,1)=I20Z*(aux0042(2,1,1)+4*F(6)*temp003(2,1,1)+4*F(4)
     &  *temp003(2,2,1)+F(5)*temp42(2,1,1)-det4*temp6221(1,k,l)+16*temp0
     &  0002(2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp00002(2,2)*ZZ(k,1,l,1
     &  )+temp00002(1,1)*ZZ(k,2,l,2)))
       temp0042(2,2,1)=I20Z*(aux0042(2,2,1)+6*F(6)*temp003(2,2,1)+2*F(4)
     &  *temp003(2,2,2)+F(5)*temp42(2,2,1)-det4*temp6222(1,k,l)+12*temp0
     &  0002(2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*temp00002(2,1)*ZZ(k,2,l,2
     &  ))
       temp0042(2,2,2)=I20Z*(aux0042(2,2,2)+8*F(6)*temp003(2,2,2)+F(5)*t
     &  emp42(2,2,2)-det4*temp6222(2,k,l)+48*temp00002(2,2)*ZZ(k,2,l,2))
       temp0043(1,1,1)=I20Z*(aux0043(1,1,1)+2*F(7)*temp003(1,1,1)+6*F(4)
     &  *temp003(3,1,1)+F(5)*temp43(1,1,1)-det4*temp6311(1,k,l)+24*temp0
     &  0002(3,1)*ZZ(k,1,l,1)+12*temp00002(1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)
     &  ))
       temp0043(2,1,1)=I20Z*(aux0043(2,1,1)+2*F(7)*temp003(2,1,1)+2*F(6)
     &  *temp003(3,1,1)+4*F(4)*temp003(3,2,1)+F(5)*temp43(2,1,1)-det4*te
     &  mp6321(1,k,l)+8*(temp00002(3,2)*ZZ(k,1,l,1)+temp00002(3,1)*(ZZ(k
     &  ,1,l,2)+ZZ(k,2,l,1)))+8*temp00002(2,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))
     &  +4*temp00002(1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp0043(2,2,1)=I20Z*(aux0043(2,2,1)+2*F(7)*temp003(2,2,1)+4*F(6)
     &  *temp003(3,2,1)+2*F(4)*temp003(3,2,2)+F(5)*temp43(2,2,1)-det4*te
     &  mp6322(1,k,l)+8*(temp00002(3,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00
     &  002(3,1)*ZZ(k,2,l,2))+4*temp00002(2,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))
     &  +8*temp00002(2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp0043(2,2,2)=I20Z*(aux0043(2,2,2)+2*F(7)*temp003(2,2,2)+6*F(6)
     &  *temp003(3,2,2)+F(5)*temp43(2,2,2)-det4*temp6322(2,k,l)+24*temp0
     &  0002(3,2)*ZZ(k,2,l,2)+12*temp00002(2,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)
     &  ))
       temp0043(3,1,1)=I20Z*(aux0043(3,1,1)+4*F(7)*temp003(3,1,1)+4*F(4)
     &  *temp003(3,3,1)+F(5)*temp43(3,1,1)-det4*temp6331(1,k,l)+16*temp0
     &  0002(3,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*(temp00002(3,3)*ZZ(k,1,l,1
     &  )+temp00002(1,1)*ZZ(k,3,l,3)))
       temp0043(3,2,1)=I20Z*(aux0043(3,2,1)+2*F(6)*temp003(3,3,1)+2*F(4)
     &  *temp003(3,3,2)+F(5)*temp43(3,2,1)-det4*temp6332(1,k,l)+4*(F(7)*
     &  temp003(3,2,1)+temp00002(3,3)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))+8*temp0
     &  0002(3,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*temp00002(3,1)*(ZZ(k,2,l,3
     &  )+ZZ(k,3,l,2))+8*temp00002(2,1)*ZZ(k,3,l,3))
       temp0043(3,2,2)=I20Z*(aux0043(3,2,2)+4*F(7)*temp003(3,2,2)+4*F(6)
     &  *temp003(3,3,2)+F(5)*temp43(3,2,2)-det4*temp6332(2,k,l)+16*temp0
     &  0002(3,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*(temp00002(3,3)*ZZ(k,2,l,2
     &  )+temp00002(2,2)*ZZ(k,3,l,3)))
       temp0043(3,3,1)=I20Z*(aux0043(3,3,1)+6*F(7)*temp003(3,3,1)+2*F(4)
     &  *temp003(3,3,3)+F(5)*temp43(3,3,1)-det4*temp6333(1,k,l)+12*temp0
     &  0002(3,3)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+24*temp00002(3,1)*ZZ(k,3,l,3
     &  ))
       temp0043(3,3,2)=I20Z*(aux0043(3,3,2)+6*F(7)*temp003(3,3,2)+2*F(6)
     &  *temp003(3,3,3)+F(5)*temp43(3,3,2)-det4*temp6333(2,k,l)+12*temp0
     &  0002(3,3)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+24*temp00002(3,2)*ZZ(k,3,l,3
     &  ))
       temp0043(3,3,3)=I20Z*(aux0043(3,3,3)+8*F(7)*temp003(3,3,3)+F(5)*t
     &  emp43(3,3,3)-det4*temp6333(3,k,l)+48*temp00002(3,3)*ZZ(k,3,l,3))
       temp0041(1,1,2)=temp0042(1,1,1)
       temp0041(1,1,3)=temp0043(1,1,1)
       temp0041(1,2,1)=temp0042(1,1,1)
       temp0041(1,2,2)=temp0042(2,1,1)
       temp0041(1,2,3)=temp0043(2,1,1)
       temp0041(1,3,1)=temp0043(1,1,1)
       temp0041(1,3,2)=temp0043(2,1,1)
       temp0041(1,3,3)=temp0043(3,1,1)
       temp0042(1,1,2)=temp0042(2,1,1)
       temp0042(1,1,3)=temp0043(2,1,1)
       temp0042(1,2,1)=temp0042(2,1,1)
       temp0042(1,2,2)=temp0042(2,2,1)
       temp0042(1,2,3)=temp0043(2,2,1)
       temp0042(1,3,1)=temp0043(2,1,1)
       temp0042(1,3,2)=temp0043(2,2,1)
       temp0042(1,3,3)=temp0043(3,2,1)
       temp0042(2,1,2)=temp0042(2,2,1)
       temp0042(2,1,3)=temp0043(2,2,1)
       temp0042(2,2,3)=temp0043(2,2,2)
       temp0042(2,3,1)=temp0043(2,2,1)
       temp0042(2,3,2)=temp0043(2,2,2)
       temp0042(2,3,3)=temp0043(3,2,2)
       temp0043(1,1,2)=temp0043(2,1,1)
       temp0043(1,1,3)=temp0043(3,1,1)
       temp0043(1,2,1)=temp0043(2,1,1)
       temp0043(1,2,2)=temp0043(2,2,1)
       temp0043(1,2,3)=temp0043(3,2,1)
       temp0043(1,3,1)=temp0043(3,1,1)
       temp0043(1,3,2)=temp0043(3,2,1)
       temp0043(1,3,3)=temp0043(3,3,1)
       temp0043(2,1,2)=temp0043(2,2,1)
       temp0043(2,1,3)=temp0043(3,2,1)
       temp0043(2,2,3)=temp0043(3,2,2)
       temp0043(2,3,1)=temp0043(3,2,1)
       temp0043(2,3,2)=temp0043(3,2,2)
       temp0043(2,3,3)=temp0043(3,3,2)
       temp0043(3,1,2)=temp0043(3,2,1)
       temp0043(3,1,3)=temp0043(3,3,1)
       temp0043(3,2,3)=temp0043(3,3,2)
       temp511(1,1,1)=IX*(aux511(1,1,1)+det4*temp6111(1,1,jj)+10*temp004
     &  1(1,1,1)*Z(jj,1))
       temp521(1,1,1)=IX*(aux521(1,1,1)+det4*temp6211(1,1,jj)+8*temp0042
     &  (1,1,1)*Z(jj,1)+2*temp0041(1,1,1)*Z(jj,2))
       temp522(1,1,1)=IX*(aux522(1,1,1)+det4*temp6221(1,1,jj)+6*temp0042
     &  (2,1,1)*Z(jj,1)+4*temp0042(1,1,1)*Z(jj,2))
       temp522(2,1,1)=IX*(aux522(2,1,1)+det4*temp6222(1,1,jj)+4*temp0042
     &  (2,2,1)*Z(jj,1)+6*temp0042(2,1,1)*Z(jj,2))
       temp522(2,2,1)=IX*(aux522(2,2,1)+det4*temp6222(2,1,jj)+2*temp0042
     &  (2,2,2)*Z(jj,1)+8*temp0042(2,2,1)*Z(jj,2))
       temp522(2,2,2)=IX*(aux522(2,2,2)+det4*temp6222(2,2,jj)+10*temp004
     &  2(2,2,2)*Z(jj,2))
       temp531(1,1,1)=IX*(aux531(1,1,1)+det4*temp6311(1,1,jj)+8*temp0043
     &  (1,1,1)*Z(jj,1)+2*temp0041(1,1,1)*Z(jj,3))
       temp532(1,1,1)=IX*(aux532(1,1,1)+det4*temp6321(1,1,jj)+6*temp0043
     &  (2,1,1)*Z(jj,1)+2*(temp0043(1,1,1)*Z(jj,2)+temp0042(1,1,1)*Z(jj,
     &  3)))
       temp532(2,1,1)=IX*(aux532(2,1,1)+det4*temp6322(1,1,jj)+4*(temp004
     &  3(2,2,1)*Z(jj,1)+temp0043(2,1,1)*Z(jj,2))+2*temp0042(2,1,1)*Z(jj
     &  ,3))
       temp532(2,2,1)=IX*(aux532(2,2,1)+det4*temp6322(2,1,jj)+6*temp0043
     &  (2,2,1)*Z(jj,2)+2*(temp0043(2,2,2)*Z(jj,1)+temp0042(2,2,1)*Z(jj,
     &  3)))
       temp532(2,2,2)=IX*(aux532(2,2,2)+det4*temp6322(2,2,jj)+8*temp0043
     &  (2,2,2)*Z(jj,2)+2*temp0042(2,2,2)*Z(jj,3))
       temp533(1,1,1)=IX*(aux533(1,1,1)+det4*temp6331(1,1,jj)+6*temp0043
     &  (3,1,1)*Z(jj,1)+4*temp0043(1,1,1)*Z(jj,3))
       temp533(2,1,1)=IX*(aux533(2,1,1)+det4*temp6332(1,1,jj)+2*temp0043
     &  (3,1,1)*Z(jj,2)+4*(temp0043(3,2,1)*Z(jj,1)+temp0043(2,1,1)*Z(jj,
     &  3)))
       temp533(2,2,1)=IX*(aux533(2,2,1)+det4*temp6332(2,1,jj)+2*temp0043
     &  (3,2,2)*Z(jj,1)+4*(temp0043(3,2,1)*Z(jj,2)+temp0043(2,2,1)*Z(jj,
     &  3)))
       temp533(2,2,2)=IX*(aux533(2,2,2)+det4*temp6332(2,2,jj)+6*temp0043
     &  (3,2,2)*Z(jj,2)+4*temp0043(2,2,2)*Z(jj,3))
       temp533(3,1,1)=IX*(aux533(3,1,1)+det4*temp6333(1,1,jj)+4*temp0043
     &  (3,3,1)*Z(jj,1)+6*temp0043(3,1,1)*Z(jj,3))
       temp533(3,2,1)=IX*(aux533(3,2,1)+det4*temp6333(2,1,jj)+2*(temp004
     &  3(3,3,2)*Z(jj,1)+temp0043(3,3,1)*Z(jj,2))+6*temp0043(3,2,1)*Z(jj
     &  ,3))
       temp533(3,2,2)=IX*(aux533(3,2,2)+det4*temp6333(2,2,jj)+4*temp0043
     &  (3,3,2)*Z(jj,2)+6*temp0043(3,2,2)*Z(jj,3))
       temp533(3,3,1)=IX*(aux533(3,3,1)+det4*temp6333(3,1,jj)+2*temp0043
     &  (3,3,3)*Z(jj,1)+8*temp0043(3,3,1)*Z(jj,3))
       temp533(3,3,2)=IX*(aux533(3,3,2)+det4*temp6333(3,2,jj)+2*temp0043
     &  (3,3,3)*Z(jj,2)+8*temp0043(3,3,2)*Z(jj,3))
       temp533(3,3,3)=IX*(aux533(3,3,3)+det4*temp6333(3,3,jj)+10*temp004
     &  3(3,3,3)*Z(jj,3))
       temp511(1,1,2)=temp521(1,1,1)
       temp511(1,1,3)=temp531(1,1,1)
       temp511(1,2,1)=temp521(1,1,1)
       temp511(1,2,2)=temp522(1,1,1)
       temp511(1,2,3)=temp532(1,1,1)
       temp511(1,3,1)=temp531(1,1,1)
       temp511(1,3,2)=temp532(1,1,1)
       temp511(1,3,3)=temp533(1,1,1)
       temp521(1,1,2)=temp522(1,1,1)
       temp521(1,1,3)=temp532(1,1,1)
       temp521(1,2,1)=temp522(1,1,1)
       temp521(1,2,2)=temp522(2,1,1)
       temp521(1,2,3)=temp532(2,1,1)
       temp521(1,3,1)=temp532(1,1,1)
       temp521(1,3,2)=temp532(2,1,1)
       temp521(1,3,3)=temp533(2,1,1)
       temp522(1,1,2)=temp522(2,1,1)
       temp522(1,1,3)=temp532(2,1,1)
       temp522(1,2,1)=temp522(2,1,1)
       temp522(1,2,2)=temp522(2,2,1)
       temp522(1,2,3)=temp532(2,2,1)
       temp522(1,3,1)=temp532(2,1,1)
       temp522(1,3,2)=temp532(2,2,1)
       temp522(1,3,3)=temp533(2,2,1)
       temp522(2,1,2)=temp522(2,2,1)
       temp522(2,1,3)=temp532(2,2,1)
       temp522(2,2,3)=temp532(2,2,2)
       temp522(2,3,1)=temp532(2,2,1)
       temp522(2,3,2)=temp532(2,2,2)
       temp522(2,3,3)=temp533(2,2,2)
       temp531(1,1,2)=temp532(1,1,1)
       temp531(1,1,3)=temp533(1,1,1)
       temp531(1,2,1)=temp532(1,1,1)
       temp531(1,2,2)=temp532(2,1,1)
       temp531(1,2,3)=temp533(2,1,1)
       temp531(1,3,1)=temp533(1,1,1)
       temp531(1,3,2)=temp533(2,1,1)
       temp531(1,3,3)=temp533(3,1,1)
       temp532(1,1,2)=temp532(2,1,1)
       temp532(1,1,3)=temp533(2,1,1)
       temp532(1,2,1)=temp532(2,1,1)
       temp532(1,2,2)=temp532(2,2,1)
       temp532(1,2,3)=temp533(2,2,1)
       temp532(1,3,1)=temp533(2,1,1)
       temp532(1,3,2)=temp533(2,2,1)
       temp532(1,3,3)=temp533(3,2,1)
       temp532(2,1,2)=temp532(2,2,1)
       temp532(2,1,3)=temp533(2,2,1)
       temp532(2,2,3)=temp533(2,2,2)
       temp532(2,3,1)=temp533(2,2,1)
       temp532(2,3,2)=temp533(2,2,2)
       temp532(2,3,3)=temp533(3,2,2)
       temp533(1,1,2)=temp533(2,1,1)
       temp533(1,1,3)=temp533(3,1,1)
       temp533(1,2,1)=temp533(2,1,1)
       temp533(1,2,2)=temp533(2,2,1)
       temp533(1,2,3)=temp533(3,2,1)
       temp533(1,3,1)=temp533(3,1,1)
       temp533(1,3,2)=temp533(3,2,1)
       temp533(1,3,3)=temp533(3,3,1)
       temp533(2,1,2)=temp533(2,2,1)
       temp533(2,1,3)=temp533(3,2,1)
       temp533(2,2,3)=temp533(3,2,2)
       temp533(2,3,1)=temp533(3,2,1)
       temp533(2,3,2)=temp533(3,2,2)
       temp533(2,3,3)=temp533(3,3,2)
       temp533(3,1,2)=temp533(3,2,1)
       temp533(3,1,3)=temp533(3,3,1)
       temp533(3,2,3)=temp533(3,3,2)
c                Step3
       temp00001(1)=I12Z*(aux00001(1)+2*tempD40000*F(4)+F(5)*temp001(1)-
     &  det4*temp003(1,k,l))
       temp00001(2)=I12Z*(aux00001(2)+2*tempD40000*F(6)+F(5)*temp001(2)-
     &  det4*temp003(2,k,l))
       temp00001(3)=I12Z*(aux00001(3)+2*tempD40000*F(7)+F(5)*temp001(3)-
     &  det4*temp003(3,k,l))
       temp003(1,1,1)=I16Z*(aux003(1,1,1)+6*F(4)*temp002(1,1)+F(5)*temp3
     &  (1,1,1)-det4*temp511(1,k,l)+24*temp00001(1)*ZZ(k,1,l,1))
       temp003(2,1,1)=I16Z*(aux003(2,1,1)+2*F(6)*temp002(1,1)+4*F(4)*tem
     &  p002(2,1)+F(5)*temp3(2,1,1)-det4*temp521(1,k,l)+8*(temp00001(2)*
     &  ZZ(k,1,l,1)+temp00001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))))
       temp003(2,2,1)=I16Z*(aux003(2,2,1)+4*F(6)*temp002(2,1)+2*F(4)*tem
     &  p002(2,2)+F(5)*temp3(2,2,1)-det4*temp522(1,k,l)+8*(temp00001(2)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00001(1)*ZZ(k,2,l,2)))
       temp003(2,2,2)=I16Z*(aux003(2,2,2)+6*F(6)*temp002(2,2)+F(5)*temp3
     &  (2,2,2)-det4*temp522(2,k,l)+24*temp00001(2)*ZZ(k,2,l,2))
       temp003(3,1,1)=I16Z*(aux003(3,1,1)+2*F(7)*temp002(1,1)+4*F(4)*tem
     &  p002(3,1)+F(5)*temp3(3,1,1)-det4*temp531(1,k,l)+8*(temp00001(3)*
     &  ZZ(k,1,l,1)+temp00001(1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))))
       temp003(3,2,1)=I16Z*(aux003(3,2,1)+2*F(7)*temp002(2,1)+2*F(6)*tem
     &  p002(3,1)+2*F(4)*temp002(3,2)+F(5)*temp3(3,2,1)-det4*temp532(1,k
     &  ,l)+4*(temp00001(3)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00001(2)*(ZZ(k
     &  ,1,l,3)+ZZ(k,3,l,1)))+4*temp00001(1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp003(3,2,2)=I16Z*(aux003(3,2,2)+2*F(7)*temp002(2,2)+4*F(6)*tem
     &  p002(3,2)+F(5)*temp3(3,2,2)-det4*temp532(2,k,l)+8*(temp00001(3)*
     &  ZZ(k,2,l,2)+temp00001(2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))))
       temp003(3,3,1)=I16Z*(aux003(3,3,1)+4*F(7)*temp002(3,1)+2*F(4)*tem
     &  p002(3,3)+F(5)*temp3(3,3,1)-det4*temp533(1,k,l)+8*(temp00001(3)*
     &  (ZZ(k,1,l,3)+ZZ(k,3,l,1))+temp00001(1)*ZZ(k,3,l,3)))
       temp003(3,3,2)=I16Z*(aux003(3,3,2)+4*F(7)*temp002(3,2)+2*F(6)*tem
     &  p002(3,3)+F(5)*temp3(3,3,2)-det4*temp533(2,k,l)+8*(temp00001(3)*
     &  (ZZ(k,2,l,3)+ZZ(k,3,l,2))+temp00001(2)*ZZ(k,3,l,3)))
       temp003(3,3,3)=I16Z*(aux003(3,3,3)+6*F(7)*temp002(3,3)+F(5)*temp3
     &  (3,3,3)-det4*temp533(3,k,l)+24*temp00001(3)*ZZ(k,3,l,3))
       temp003(1,1,2)=temp003(2,1,1)
       temp003(1,1,3)=temp003(3,1,1)
       temp003(1,2,1)=temp003(2,1,1)
       temp003(1,2,2)=temp003(2,2,1)
       temp003(1,2,3)=temp003(3,2,1)
       temp003(1,3,1)=temp003(3,1,1)
       temp003(1,3,2)=temp003(3,2,1)
       temp003(1,3,3)=temp003(3,3,1)
       temp003(2,1,2)=temp003(2,2,1)
       temp003(2,1,3)=temp003(3,2,1)
       temp003(2,2,3)=temp003(3,2,2)
       temp003(2,3,1)=temp003(3,2,1)
       temp003(2,3,2)=temp003(3,2,2)
       temp003(2,3,3)=temp003(3,3,2)
       temp003(3,1,2)=temp003(3,2,1)
       temp003(3,1,3)=temp003(3,3,1)
       temp003(3,2,3)=temp003(3,3,2)
       temp41(1,1,1)=IX*(aux41(1,1,1)+det4*temp511(1,1,jj)+8*temp003(1,1
     &  ,1)*Z(jj,1))
       temp42(1,1,1)=IX*(aux42(1,1,1)+det4*temp521(1,1,jj)+6*temp003(2,1
     &  ,1)*Z(jj,1)+2*temp003(1,1,1)*Z(jj,2))
       temp42(2,1,1)=IX*(aux42(2,1,1)+det4*temp522(1,1,jj)+4*(temp003(2,
     &  2,1)*Z(jj,1)+temp003(2,1,1)*Z(jj,2)))
       temp42(2,2,1)=IX*(aux42(2,2,1)+det4*temp522(2,1,jj)+2*temp003(2,2
     &  ,2)*Z(jj,1)+6*temp003(2,2,1)*Z(jj,2))
       temp42(2,2,2)=IX*(aux42(2,2,2)+det4*temp522(2,2,jj)+8*temp003(2,2
     &  ,2)*Z(jj,2))
       temp43(1,1,1)=IX*(aux43(1,1,1)+det4*temp531(1,1,jj)+6*temp003(3,1
     &  ,1)*Z(jj,1)+2*temp003(1,1,1)*Z(jj,3))
       temp43(2,1,1)=IX*(aux43(2,1,1)+det4*temp532(1,1,jj)+4*temp003(3,2
     &  ,1)*Z(jj,1)+2*(temp003(3,1,1)*Z(jj,2)+temp003(2,1,1)*Z(jj,3)))
       temp43(2,2,1)=IX*(aux43(2,2,1)+det4*temp532(2,1,jj)+4*temp003(3,2
     &  ,1)*Z(jj,2)+2*(temp003(3,2,2)*Z(jj,1)+temp003(2,2,1)*Z(jj,3)))
       temp43(2,2,2)=IX*(aux43(2,2,2)+det4*temp532(2,2,jj)+6*temp003(3,2
     &  ,2)*Z(jj,2)+2*temp003(2,2,2)*Z(jj,3))
       temp43(3,1,1)=IX*(aux43(3,1,1)+det4*temp533(1,1,jj)+4*(temp003(3,
     &  3,1)*Z(jj,1)+temp003(3,1,1)*Z(jj,3)))
       temp43(3,2,1)=IX*(aux43(3,2,1)+det4*temp533(2,1,jj)+2*(temp003(3,
     &  3,2)*Z(jj,1)+temp003(3,3,1)*Z(jj,2))+4*temp003(3,2,1)*Z(jj,3))
       temp43(3,2,2)=IX*(aux43(3,2,2)+det4*temp533(2,2,jj)+4*(temp003(3,
     &  3,2)*Z(jj,2)+temp003(3,2,2)*Z(jj,3)))
       temp43(3,3,1)=IX*(aux43(3,3,1)+det4*temp533(3,1,jj)+2*temp003(3,3
     &  ,3)*Z(jj,1)+6*temp003(3,3,1)*Z(jj,3))
       temp43(3,3,2)=IX*(aux43(3,3,2)+det4*temp533(3,2,jj)+2*temp003(3,3
     &  ,3)*Z(jj,2)+6*temp003(3,3,2)*Z(jj,3))
       temp43(3,3,3)=IX*(aux43(3,3,3)+det4*temp533(3,3,jj)+8*temp003(3,3
     &  ,3)*Z(jj,3))
       temp41(1,1,2)=temp42(1,1,1)
       temp41(1,1,3)=temp43(1,1,1)
       temp41(1,2,1)=temp42(1,1,1)
       temp41(1,2,2)=temp42(2,1,1)
       temp41(1,2,3)=temp43(2,1,1)
       temp41(1,3,1)=temp43(1,1,1)
       temp41(1,3,2)=temp43(2,1,1)
       temp41(1,3,3)=temp43(3,1,1)
       temp42(1,1,2)=temp42(2,1,1)
       temp42(1,1,3)=temp43(2,1,1)
       temp42(1,2,1)=temp42(2,1,1)
       temp42(1,2,2)=temp42(2,2,1)
       temp42(1,2,3)=temp43(2,2,1)
       temp42(1,3,1)=temp43(2,1,1)
       temp42(1,3,2)=temp43(2,2,1)
       temp42(1,3,3)=temp43(3,2,1)
       temp42(2,1,2)=temp42(2,2,1)
       temp42(2,1,3)=temp43(2,2,1)
       temp42(2,2,3)=temp43(2,2,2)
       temp42(2,3,1)=temp43(2,2,1)
       temp42(2,3,2)=temp43(2,2,2)
       temp42(2,3,3)=temp43(3,2,2)
       temp43(1,1,2)=temp43(2,1,1)
       temp43(1,1,3)=temp43(3,1,1)
       temp43(1,2,1)=temp43(2,1,1)
       temp43(1,2,2)=temp43(2,2,1)
       temp43(1,2,3)=temp43(3,2,1)
       temp43(1,3,1)=temp43(3,1,1)
       temp43(1,3,2)=temp43(3,2,1)
       temp43(1,3,3)=temp43(3,3,1)
       temp43(2,1,2)=temp43(2,2,1)
       temp43(2,1,3)=temp43(3,2,1)
       temp43(2,2,3)=temp43(3,2,2)
       temp43(2,3,1)=temp43(3,2,1)
       temp43(2,3,2)=temp43(3,2,2)
       temp43(2,3,3)=temp43(3,3,2)
       temp43(3,1,2)=temp43(3,2,1)
       temp43(3,1,3)=temp43(3,3,1)
       temp43(3,2,3)=temp43(3,3,2)
c                Step4
       tempD40000=I8Z*(auxD40000+tempD400*F(5)-det4*temp002(k,l))
       temp002(1,1)=I12Z*(aux002(1,1)+4*F(4)*temp001(1)+F(5)*temp2(1,1)-
     &  det4*temp41(1,k,l)+8*tempD40000*ZZ(k,1,l,1))
       temp002(2,1)=I12Z*(aux002(2,1)+2*(F(6)*temp001(1)+F(4)*temp001(2)
     &  )+F(5)*temp2(2,1)-det4*temp42(1,k,l)+4*tempD40000*(ZZ(k,1,l,2)+Z
     &  Z(k,2,l,1)))
       temp002(2,2)=I12Z*(aux002(2,2)+4*F(6)*temp001(2)+F(5)*temp2(2,2)-
     &  det4*temp42(2,k,l)+8*tempD40000*ZZ(k,2,l,2))
       temp002(3,1)=I12Z*(aux002(3,1)+2*(F(7)*temp001(1)+F(4)*temp001(3)
     &  )+F(5)*temp2(3,1)-det4*temp43(1,k,l)+4*tempD40000*(ZZ(k,1,l,3)+Z
     &  Z(k,3,l,1)))
       temp002(3,2)=I12Z*(aux002(3,2)+2*(F(7)*temp001(2)+F(6)*temp001(3)
     &  )+F(5)*temp2(3,2)-det4*temp43(2,k,l)+4*tempD40000*(ZZ(k,2,l,3)+Z
     &  Z(k,3,l,2)))
       temp002(3,3)=I12Z*(aux002(3,3)+4*F(7)*temp001(3)+F(5)*temp2(3,3)-
     &  det4*temp43(3,k,l)+8*tempD40000*ZZ(k,3,l,3))
       temp002(1,2)=temp002(2,1)
       temp002(1,3)=temp002(3,1)
       temp002(2,3)=temp002(3,2)
       temp3(1,1,1)=IX*(aux3(1,1,1)+det4*temp41(1,1,jj)+6*temp002(1,1)*Z
     &  (jj,1))
       temp3(2,1,1)=IX*(aux3(2,1,1)+det4*temp42(1,1,jj)+4*temp002(2,1)*Z
     &  (jj,1)+2*temp002(1,1)*Z(jj,2))
       temp3(2,2,1)=IX*(aux3(2,2,1)+det4*temp42(2,1,jj)+2*temp002(2,2)*Z
     &  (jj,1)+4*temp002(2,1)*Z(jj,2))
       temp3(2,2,2)=IX*(aux3(2,2,2)+det4*temp42(2,2,jj)+6*temp002(2,2)*Z
     &  (jj,2))
       temp3(3,1,1)=IX*(aux3(3,1,1)+det4*temp43(1,1,jj)+4*temp002(3,1)*Z
     &  (jj,1)+2*temp002(1,1)*Z(jj,3))
       temp3(3,2,1)=IX*(aux3(3,2,1)+det4*temp43(2,1,jj)+2*(temp002(3,2)*
     &  Z(jj,1)+temp002(3,1)*Z(jj,2))+2*temp002(2,1)*Z(jj,3))
       temp3(3,2,2)=IX*(aux3(3,2,2)+det4*temp43(2,2,jj)+4*temp002(3,2)*Z
     &  (jj,2)+2*temp002(2,2)*Z(jj,3))
       temp3(3,3,1)=IX*(aux3(3,3,1)+det4*temp43(3,1,jj)+2*temp002(3,3)*Z
     &  (jj,1)+4*temp002(3,1)*Z(jj,3))
       temp3(3,3,2)=IX*(aux3(3,3,2)+det4*temp43(3,2,jj)+2*temp002(3,3)*Z
     &  (jj,2)+4*temp002(3,2)*Z(jj,3))
       temp3(3,3,3)=IX*(aux3(3,3,3)+det4*temp43(3,3,jj)+6*temp002(3,3)*Z
     &  (jj,3))
       temp3(1,1,2)=temp3(2,1,1)
       temp3(1,1,3)=temp3(3,1,1)
       temp3(1,2,1)=temp3(2,1,1)
       temp3(1,2,2)=temp3(2,2,1)
       temp3(1,2,3)=temp3(3,2,1)
       temp3(1,3,1)=temp3(3,1,1)
       temp3(1,3,2)=temp3(3,2,1)
       temp3(1,3,3)=temp3(3,3,1)
       temp3(2,1,2)=temp3(2,2,1)
       temp3(2,1,3)=temp3(3,2,1)
       temp3(2,2,3)=temp3(3,2,2)
       temp3(2,3,1)=temp3(3,2,1)
       temp3(2,3,2)=temp3(3,2,2)
       temp3(2,3,3)=temp3(3,3,2)
       temp3(3,1,2)=temp3(3,2,1)
       temp3(3,1,3)=temp3(3,3,1)
       temp3(3,2,3)=temp3(3,3,2)
c                Step5
       temp001(1)=I8Z*(aux001(1)+2*tempD400*F(4)+F(5)*temp1(1)-det4*temp
     &  3(1,k,l))
       temp001(2)=I8Z*(aux001(2)+2*tempD400*F(6)+F(5)*temp1(2)-det4*temp
     &  3(2,k,l))
       temp001(3)=I8Z*(aux001(3)+2*tempD400*F(7)+F(5)*temp1(3)-det4*temp
     &  3(3,k,l))
       temp2(1,1)=IX*(aux2(1,1)+det4*temp3(1,1,jj)+4*temp001(1)*Z(jj,1))
       temp2(2,1)=IX*(aux2(2,1)+det4*temp3(2,1,jj)+2*temp001(2)*Z(jj,1)+
     &  2*temp001(1)*Z(jj,2))
       temp2(2,2)=IX*(aux2(2,2)+det4*temp3(2,2,jj)+4*temp001(2)*Z(jj,2))
       temp2(3,1)=IX*(aux2(3,1)+det4*temp3(3,1,jj)+2*temp001(3)*Z(jj,1)+
     &  2*temp001(1)*Z(jj,3))
       temp2(3,2)=IX*(aux2(3,2)+det4*temp3(3,2,jj)+2*temp001(3)*Z(jj,2)+
     &  2*temp001(2)*Z(jj,3))
       temp2(3,3)=IX*(aux2(3,3)+det4*temp3(3,3,jj)+4*temp001(3)*Z(jj,3))
       temp2(1,2)=temp2(2,1)
       temp2(1,3)=temp2(3,1)
       temp2(2,3)=temp2(3,2)
c                Step6
       tempD400=I4Z*(auxD400+tempD40*F(5)-det4*temp2(k,l))
       temp1(1)=IX*(aux1(1)+det4*temp2(1,jj)+2*tempD400*Z(jj,1))
       temp1(2)=IX*(aux1(2)+det4*temp2(2,jj)+2*tempD400*Z(jj,2))
       temp1(3)=IX*(aux1(3)+det4*temp2(3,jj)+2*tempD400*Z(jj,3))
c                Step7
       tempD40=IX*(auxD40+det4*temp1(jj))

        ac=1
        accuracyDR(1,0,ac)  = abs(D0        /tempD40         -1d0)
        accuracyDR(1,1,ac)  = abs(Dij(1,1)  /temp1(1)        -1d0)
        accuracyDR(2,1,ac)  = abs(Dij(2,1)  /temp1(2)        -1d0)
        accuracyDR(3,1,ac)  = abs(Dij(3,1)  /temp1(3)        -1d0)
        accuracyDR(1,2,ac)  = abs(Dij(1,2)  /temp2(1,1)      -1d0)
        accuracyDR(2,2,ac)  = abs(Dij(2,2)  /temp2(2,2)      -1d0)
        accuracyDR(3,2,ac)  = abs(Dij(3,2)  /temp2(3,3)      -1d0)
        accuracyDR(4,2,ac)  = abs(Dij(4,2)  /temp2(2,1)      -1d0)
        accuracyDR(5,2,ac)  = abs(Dij(5,2)  /temp2(3,1)      -1d0)
        accuracyDR(6,2,ac)  = abs(Dij(6,2)  /temp2(3,2)      -1d0)
        accuracyDR(7,2,ac)  = abs(Dij(7,2)  /tempD400        -1d0)
        accuracyDR(1,3,ac)  = abs(Dij(1,3)  /temp3(1,1,1)    -1d0)
        accuracyDR(2,3,ac)  = abs(Dij(2,3)  /temp3(2,2,2)    -1d0)
        accuracyDR(3,3,ac)  = abs(Dij(3,3)  /temp3(3,3,3)    -1d0)
        accuracyDR(4,3,ac)  = abs(Dij(4,3)  /temp3(2,1,1)    -1d0)
        accuracyDR(5,3,ac)  = abs(Dij(5,3)  /temp3(3,1,1)    -1d0)
        accuracyDR(6,3,ac)  = abs(Dij(6,3)  /temp3(2,2,1)    -1d0)
        accuracyDR(7,3,ac)  = abs(Dij(7,3)  /temp3(3,3,1)    -1d0)
        accuracyDR(8,3,ac)  = abs(Dij(8,3)  /temp3(3,2,2)    -1d0)
        accuracyDR(9,3,ac)  = abs(Dij(9,3)  /temp3(3,3,2)    -1d0)
        accuracyDR(10,3,ac) = abs(Dij(10,3) /temp3(3,2,1)    -1d0)
        accuracyDR(11,3,ac) = abs(Dij(11,3) /temp001(1)      -1d0)
        accuracyDR(12,3,ac) = abs(Dij(12,3) /temp001(2)      -1d0)
        accuracyDR(13,3,ac) = abs(Dij(13,3) /temp001(3)      -1d0)
        accuracyDR(7,1,ac)  = abs(Dij(7,1)  /temp002(1,1)    -1d0)
        accuracyDR(8,1,ac)  = abs(Dij(8,1)  /temp002(2,2)    -1d0)
        accuracyDR(9,1,ac)  = abs(Dij(9,1)  /temp002(3,3)    -1d0)
        accuracyDR(10,1,ac) = abs(Dij(10,1) /temp002(2,1)    -1d0)
        accuracyDR(11,1,ac) = abs(Dij(11,1) /temp002(3,1)    -1d0)
        accuracyDR(12,1,ac) = abs(Dij(12,1) /temp002(3,2)    -1d0)
        accuracyDR(13,1,ac) = abs(Dij(13,1) /tempD40000      -1d0)
        accuracyDR(1,4,ac)  = abs(Dij(1,4)  /temp41(1,1,1)   -1d0)
        accuracyDR(2,4,ac)  = abs(Dij(2,4)  /temp42(2,2,2)   -1d0)
        accuracyDR(3,4,ac)  = abs(Dij(3,4)  /temp43(3,3,3)   -1d0)
        accuracyDR(4,4,ac)  = abs(Dij(4,4)  /temp42(1,1,1)   -1d0)
        accuracyDR(5,4,ac)  = abs(Dij(5,4)  /temp43(1,1,1)   -1d0)
        accuracyDR(6,4,ac)  = abs(Dij(6,4)  /temp42(2,1,1)   -1d0)
        accuracyDR(7,4,ac)  = abs(Dij(7,4)  /temp43(2,1,1)   -1d0)
        accuracyDR(8,4,ac)  = abs(Dij(8,4)  /temp43(3,1,1)   -1d0)
        accuracyDR(9,4,ac)  = abs(Dij(9,4)  /temp42(2,2,1)   -1d0)
        accuracyDR(10,4,ac) = abs(Dij(10,4) /temp43(2,2,1)   -1d0)
        accuracyDR(11,4,ac) = abs(Dij(11,4) /temp43(3,2,1)   -1d0)
        accuracyDR(12,4,ac) = abs(Dij(12,4) /temp43(3,3,1)   -1d0)
        accuracyDR(13,4,ac) = abs(Dij(13,4) /temp43(2,2,2)   -1d0)
        accuracyDR(14,4,ac) = abs(Dij(14,4) /temp43(3,2,2)   -1d0)
        accuracyDR(15,4,ac) = abs(Dij(15,4) /temp43(3,3,2)   -1d0)
        accuracyDR(16,4,ac) = abs(Dij(16,4) /temp002(1,1)    -1d0)
        accuracyDR(17,4,ac) = abs(Dij(17,4) /temp002(2,2)    -1d0)
        accuracyDR(18,4,ac) = abs(Dij(18,4) /temp002(3,3)    -1d0)
        accuracyDR(19,4,ac) = abs(Dij(19,4) /temp002(2,1)    -1d0)
        accuracyDR(20,4,ac) = abs(Dij(20,4) /temp002(3,1)    -1d0)
        accuracyDR(21,4,ac) = abs(Dij(21,4) /temp002(3,2)    -1d0)
        accuracyDR(22,4,ac) = abs(Dij(22,4) /tempD40000      -1d0)


      DO I1=0,4
           accuracyD(i1,ac)=accuracyDR(1,i1,ac)
            DO I2=2,INDEX(I1)  
       if(accuracyDR(i2,i1,ac).gt.accuracyD(i1,ac)) then
          accuracyD(i1,ac)=accuracyDR(i2,i1,ac)
       endif
          enddo
        enddo

        D0=tempD40
        Dij(1,1)=temp1(1)
        Dij(2,1)=temp1(2)
        Dij(3,1)=temp1(3)
        Dij(1,2)=temp2(1,1)
        Dij(2,2)=temp2(2,2)
        Dij(3,2)=temp2(3,3)
        Dij(4,2)=temp2(2,1)
        Dij(5,2)=temp2(3,1)
        Dij(6,2)=temp2(3,2)
        Dij(7,2)=tempD400
        Dij(1,3)=temp3(1,1,1)
        Dij(2,3)=temp3(2,2,2)
        Dij(3,3)=temp3(3,3,3)
        Dij(4,3)=temp3(2,1,1)
        Dij(5,3)=temp3(3,1,1)
        Dij(6,3)=temp3(2,2,1)
        Dij(7,3)=temp3(3,3,1)
        Dij(8,3)=temp3(3,2,2)
        Dij(9,3)=temp3(3,3,2)
        Dij(10,3)=temp3(3,2,1)
        Dij(11,3)=temp001(1)
        Dij(12,3)=temp001(2)
        Dij(13,3)=temp001(3)
        Dij(7,1)=temp002(1,1)
        Dij(8,1)=temp002(2,2)
        Dij(9,1)=temp002(3,3)
        Dij(10,1)=temp002(2,1)
        Dij(11,1)=temp002(3,1)
        Dij(12,1)=temp002(3,2)
        Dij(13,1)=tempD40000
        Dij(1,4)=temp41(1,1,1)
         Dij(2,4)=temp42(2,2,2)
         Dij(3,4)=temp43(3,3,3)
         Dij(4,4)=temp42(1,1,1)
         Dij(5,4)=temp43(1,1,1)
         Dij(6,4)=temp42(2,1,1)
         Dij(7,4)=temp43(2,1,1)
         Dij(8,4)=temp43(3,1,1)
         Dij(9,4)=temp42(2,2,1)
         Dij(10,4)=temp43(2,2,1)
         Dij(11,4)=temp43(3,2,1)
         Dij(12,4)=temp43(3,3,1)
         Dij(13,4)=temp43(2,2,2)
         Dij(14,4)=temp43(3,2,2)
         Dij(15,4)=temp43(3,3,2)
         Dij(16,4)=temp002(1,1)
         Dij(17,4)=temp002(2,2)
         Dij(18,4)=temp002(3,3)
         Dij(19,4)=temp002(2,1)
         Dij(20,4)=temp002(3,1)
         Dij(21,4)=temp002(3,2)
         Dij(22,4)=tempD40000
        if(printmy) then
       print*, "tempD40 6",tempD40
        endif
       if(order.eq.6) goto 500
c                Iteration7
c                Step1
       S300000000=2*Cij234(16,6)
       S300000021(1)=2*Cij234(9,4)
       S300000021(2)=-2*Cij234(11,5)
       S300000021(3)=-2*Cij234(12,5)
       S300000022(1)=-2*Cij234(11,5)
       S300000022(2)=2*Cij234(13,6)
       S300000022(3)=2*Cij234(15,6)
       S300000023(1)=-2*Cij234(12,5)
       S300000023(2)=2*Cij234(15,6)
       S300000023(3)=2*Cij234(14,6)
       S3h00000021(1)=Cij134(19,7)+Cij234(16,6)
       S3h00000021(2)=Cij134(19,7)-Cij234(19,7)
       S3h00000021(3)=Cij134(20,7)-Cij234(20,7)
       S3h00000022(1)=Cij124(19,7)-Cij134(19,7)
       S3h00000022(2)=Cij124(20,7)-Cij134(19,7)
       S3h00000022(3)=Cij124(20,7)-Cij134(20,7)
       S3h00000023(1)=Cij123(19,7)-Cij124(19,7)
       S3h00000023(2)=Cij123(20,7)-Cij124(20,7)
       S3h00000023(3)=-Cij124(20,7)
       auxD400000000=-(F(1)*S3h0000001(1))-F(2)*S3h0000001(2)-F(3)*S3h00
     &  00001(3)+S3h00000021(k)*Z(1,l)+S3h00000022(k)*Z(2,l)+S3h00000023
     &  (k)*Z(3,l)+(Inv80+S300000000-S3h00000021(1)-S3h00000022(2)-S3h00
     &  000023(3))*Z(k,l)
       tempD400000000=I16Z*(auxD400000000+tempD4000000*F(5))
       S3h00004111(1)=Cij134(15,7)+Cij234(9,4)
       S3h00004111(2)=Cij134(15,7)-Cij234(11,5)
       S3h00004111(3)=Cij134(17,7)-Cij234(12,5)
       S3h00004121(1)=Cij134(15,7)-Cij234(11,5)
       S3h00004121(2)=Cij134(15,7)+Cij234(13,6)
       S3h00004121(3)=Cij134(17,7)+Cij234(15,6)
       S3h00004122(1)=Cij134(15,7)+Cij234(13,6)
       S3h00004122(2)=Cij134(15,7)-Cij234(15,7)
       S3h00004122(3)=Cij134(17,7)-Cij234(17,7)
       S3h00004131(1)=Cij134(17,7)-Cij234(12,5)
       S3h00004131(2)=Cij134(17,7)+Cij234(15,6)
       S3h00004131(3)=Cij134(18,7)+Cij234(14,6)
       S3h00004132(1)=Cij134(17,7)+Cij234(15,6)
       S3h00004132(2)=Cij134(17,7)-Cij234(17,7)
       S3h00004132(3)=Cij134(18,7)-Cij234(18,7)
       S3h00004133(1)=Cij134(18,7)+Cij234(14,6)
       S3h00004133(2)=Cij134(18,7)-Cij234(18,7)
       S3h00004133(3)=Cij134(16,7)-Cij234(16,7)
       S3h00004211(1)=Cij124(15,7)-Cij134(15,7)
       S3h00004211(2)=Cij124(17,7)-Cij134(15,7)
       S3h00004211(3)=Cij124(17,7)-Cij134(17,7)
       S3h00004221(1)=Cij124(17,7)-Cij134(15,7)
       S3h00004221(2)=Cij124(18,7)-Cij134(15,7)
       S3h00004221(3)=Cij124(18,7)-Cij134(17,7)
       S3h00004222(1)=Cij124(18,7)-Cij134(15,7)
       S3h00004222(2)=Cij124(16,7)-Cij134(15,7)
       S3h00004222(3)=Cij124(16,7)-Cij134(17,7)
       S3h00004231(1)=Cij124(17,7)-Cij134(17,7)
       S3h00004231(2)=Cij124(18,7)-Cij134(17,7)
       S3h00004231(3)=Cij124(18,7)-Cij134(18,7)
       S3h00004232(1)=Cij124(18,7)-Cij134(17,7)
       S3h00004232(2)=Cij124(16,7)-Cij134(17,7)
       S3h00004232(3)=Cij124(16,7)-Cij134(18,7)
       S3h00004233(1)=Cij124(18,7)-Cij134(18,7)
       S3h00004233(2)=Cij124(16,7)-Cij134(18,7)
       S3h00004233(3)=Cij124(16,7)-Cij134(16,7)
       S3h00004311(1)=Cij123(15,7)-Cij124(15,7)
       S3h00004311(2)=Cij123(17,7)-Cij124(17,7)
       S3h00004311(3)=-Cij124(17,7)
       S3h00004321(1)=Cij123(17,7)-Cij124(17,7)
       S3h00004321(2)=Cij123(18,7)-Cij124(18,7)
       S3h00004321(3)=-Cij124(18,7)
       S3h00004322(1)=Cij123(18,7)-Cij124(18,7)
       S3h00004322(2)=Cij123(16,7)-Cij124(16,7)
       S3h00004322(3)=-Cij124(16,7)
       S3h00004331(1)=-Cij124(17,7)
       S3h00004331(2)=-Cij124(18,7)
       S3h00004331(3)=-Cij124(18,7)
       S3h00004332(1)=-Cij124(18,7)
       S3h00004332(2)=-Cij124(16,7)
       S3h00004332(3)=-Cij124(16,7)
       S3h00004333(1)=-Cij124(18,7)
       S3h00004333(2)=-Cij124(16,7)
       S3h00004333(3)=-Cij124(16,7)
       aux0000002(1,1)=-(F(1)*S3h0000311(1))-F(2)*S3h0000321(1)-F(3)*S3h
     &  0000331(1)+S3h00004111(k)*Z(1,l)+S3h00004211(k)*Z(2,l)+S3h000043
     &  11(k)*Z(3,l)+(Inv6011+S300000021(1)-S3h00004111(1)-S3h00004221(1
     &  )-S3h00004331(1))*Z(k,l)-4*S3h00000021(1)*ZZ(k,1,l,1)-4*S3h00000
     &  022(1)*ZZ(k,1,l,2)-4*S3h00000023(1)*ZZ(k,1,l,3)
       aux0000002(2,1)=-(F(1)*S3h0000312(1))-F(2)*S3h0000322(1)-F(3)*S3h
     &  0000332(1)+S3h00004121(k)*Z(1,l)+S3h00004221(k)*Z(2,l)+S3h000043
     &  21(k)*Z(3,l)+(Inv6021+S300000022(1)-S3h00004121(1)-S3h00004222(1
     &  )-S3h00004332(1))*Z(k,l)-2*S3h00000021(2)*ZZ(k,1,l,1)-2*S3h00000
     &  022(2)*ZZ(k,1,l,2)-2*S3h00000023(2)*ZZ(k,1,l,3)-2*S3h00000021(1)
     &  *ZZ(k,2,l,1)-2*S3h00000022(1)*ZZ(k,2,l,2)-2*S3h00000023(1)*ZZ(k,
     &  2,l,3)
       aux0000002(2,2)=-(F(1)*S3h0000312(2))-F(2)*S3h0000322(2)-F(3)*S3h
     &  0000332(2)+S3h00004122(k)*Z(1,l)+S3h00004222(k)*Z(2,l)+S3h000043
     &  22(k)*Z(3,l)+(Inv6022+S300000022(2)-S3h00004122(1)-S3h00004222(2
     &  )-S3h00004332(2))*Z(k,l)-4*S3h00000021(2)*ZZ(k,2,l,1)-4*S3h00000
     &  022(2)*ZZ(k,2,l,2)-4*S3h00000023(2)*ZZ(k,2,l,3)
       aux0000002(3,1)=-(F(1)*S3h0000313(1))-F(2)*S3h0000323(1)-F(3)*S3h
     &  0000333(1)+S3h00004131(k)*Z(1,l)+S3h00004231(k)*Z(2,l)+S3h000043
     &  31(k)*Z(3,l)+(Inv6031+S300000023(1)-S3h00004131(1)-S3h00004232(1
     &  )-S3h00004333(1))*Z(k,l)-2*S3h00000021(3)*ZZ(k,1,l,1)-2*S3h00000
     &  022(3)*ZZ(k,1,l,2)-2*S3h00000023(3)*ZZ(k,1,l,3)-2*S3h00000021(1)
     &  *ZZ(k,3,l,1)-2*S3h00000022(1)*ZZ(k,3,l,2)-2*S3h00000023(1)*ZZ(k,
     &  3,l,3)
       aux0000002(3,2)=-(F(1)*S3h0000313(2))-F(2)*S3h0000323(2)-F(3)*S3h
     &  0000333(2)+S3h00004132(k)*Z(1,l)+S3h00004232(k)*Z(2,l)+S3h000043
     &  32(k)*Z(3,l)+(Inv6032+S300000023(2)-S3h00004132(1)-S3h00004232(2
     &  )-S3h00004333(2))*Z(k,l)-2*S3h00000021(3)*ZZ(k,2,l,1)-2*S3h00000
     &  022(3)*ZZ(k,2,l,2)-2*S3h00000023(3)*ZZ(k,2,l,3)-2*S3h00000021(2)
     &  *ZZ(k,3,l,1)-2*S3h00000022(2)*ZZ(k,3,l,2)-2*S3h00000023(2)*ZZ(k,
     &  3,l,3)
       aux0000002(3,3)=-(F(1)*S3h0000313(3))-F(2)*S3h0000323(3)-F(3)*S3h
     &  0000333(3)+S3h00004133(k)*Z(1,l)+S3h00004233(k)*Z(2,l)+S3h000043
     &  33(k)*Z(3,l)+(Inv6033+S300000023(3)-S3h00004133(1)-S3h00004233(2
     &  )-S3h00004333(3))*Z(k,l)-4*S3h00000021(3)*ZZ(k,3,l,1)-4*S3h00000
     &  022(3)*ZZ(k,3,l,2)-4*S3h00000023(3)*ZZ(k,3,l,3)
       temp0000002(1,1)=I20Z*(aux0000002(1,1)+4*F(4)*temp0000001(1)+F(5)
     &  *temp00002(1,1)+8*tempD400000000*ZZ(k,1,l,1))
       temp0000002(2,1)=I20Z*(aux0000002(2,1)+2*(F(6)*temp0000001(1)+F(4
     &  )*temp0000001(2))+F(5)*temp00002(2,1)+4*tempD400000000*(ZZ(k,1,l
     &  ,2)+ZZ(k,2,l,1)))
       temp0000002(2,2)=I20Z*(aux0000002(2,2)+4*F(6)*temp0000001(2)+F(5)
     &  *temp00002(2,2)+8*tempD400000000*ZZ(k,2,l,2))
       temp0000002(3,1)=I20Z*(aux0000002(3,1)+2*(F(7)*temp0000001(1)+F(4
     &  )*temp0000001(3))+F(5)*temp00002(3,1)+4*tempD400000000*(ZZ(k,1,l
     &  ,3)+ZZ(k,3,l,1)))
       temp0000002(3,2)=I20Z*(aux0000002(3,2)+2*(F(7)*temp0000001(2)+F(6
     &  )*temp0000001(3))+F(5)*temp00002(3,2)+4*tempD400000000*(ZZ(k,2,l
     &  ,3)+ZZ(k,3,l,2)))
       temp0000002(3,3)=I20Z*(aux0000002(3,3)+4*F(7)*temp0000001(3)+F(5)
     &  *temp00002(3,3)+8*tempD400000000*ZZ(k,3,l,3))
       temp0000002(1,2)=temp0000002(2,1)
       temp0000002(1,3)=temp0000002(3,1)
       temp0000002(2,3)=temp0000002(3,2)
       S300004111(1)=2*Cij234(4,2)
       S300004111(2)=-2*Cij234(5,3)
       S300004111(3)=-2*Cij234(6,3)
       S300004121(1)=-2*Cij234(5,3)
       S300004121(2)=2*Cij234(6,4)
       S300004121(3)=2*Cij234(8,4)
       S300004122(1)=2*Cij234(6,4)
       S300004122(2)=-2*Cij234(7,5)
       S300004122(3)=-2*Cij234(9,5)
       S300004131(1)=-2*Cij234(6,3)
       S300004131(2)=2*Cij234(8,4)
       S300004131(3)=2*Cij234(7,4)
       S300004132(1)=2*Cij234(8,4)
       S300004132(2)=-2*Cij234(9,5)
       S300004132(3)=-2*Cij234(10,5)
       S300004133(1)=2*Cij234(7,4)
       S300004133(2)=-2*Cij234(10,5)
       S300004133(3)=-2*Cij234(8,5)
       S300004211(1)=-2*Cij234(5,3)
       S300004211(2)=2*Cij234(6,4)
       S300004211(3)=2*Cij234(8,4)
       S300004221(1)=2*Cij234(6,4)
       S300004221(2)=-2*Cij234(7,5)
       S300004221(3)=-2*Cij234(9,5)
       S300004222(1)=-2*Cij234(7,5)
       S300004222(2)=2*Cij234(8,6)
       S300004222(3)=2*Cij234(10,6)
       S300004231(1)=2*Cij234(8,4)
       S300004231(2)=-2*Cij234(9,5)
       S300004231(3)=-2*Cij234(10,5)
       S300004232(1)=-2*Cij234(9,5)
       S300004232(2)=2*Cij234(10,6)
       S300004232(3)=2*Cij234(11,6)
       S300004233(1)=-2*Cij234(10,5)
       S300004233(2)=2*Cij234(11,6)
       S300004233(3)=2*Cij234(12,6)
       S300004311(1)=-2*Cij234(6,3)
       S300004311(2)=2*Cij234(8,4)
       S300004311(3)=2*Cij234(7,4)
       S300004321(1)=2*Cij234(8,4)
       S300004321(2)=-2*Cij234(9,5)
       S300004321(3)=-2*Cij234(10,5)
       S300004322(1)=-2*Cij234(9,5)
       S300004322(2)=2*Cij234(10,6)
       S300004322(3)=2*Cij234(11,6)
       S300004331(1)=2*Cij234(7,4)
       S300004331(2)=-2*Cij234(10,5)
       S300004331(3)=-2*Cij234(8,5)
       S300004332(1)=-2*Cij234(10,5)
       S300004332(2)=2*Cij234(11,6)
       S300004332(3)=2*Cij234(12,6)
       S300004333(1)=-2*Cij234(8,5)
       S300004333(2)=2*Cij234(12,6)
       S300004333(3)=2*Cij234(9,6)
       S300611111(1)=2*C0234
       S300611111(2)=-2*Cij234(1,1)
       S300611111(3)=-2*Cij234(2,1)
       S300612111(1)=-2*Cij234(1,1)
       S300612111(2)=2*Cij234(1,2)
       S300612111(3)=2*Cij234(3,2)
       S300612211(1)=2*Cij234(1,2)
       S300612211(2)=-2*Cij234(1,3)
       S300612211(3)=-2*Cij234(3,3)
       S300612221(1)=-2*Cij234(1,3)
       S300612221(2)=2*Cij234(1,4)
       S300612221(3)=2*Cij234(3,4)
       S300612222(1)=2*Cij234(1,4)
       S300612222(2)=-2*Cij234(1,5)
       S300612222(3)=-2*Cij234(3,5)
       S300613111(1)=-2*Cij234(2,1)
       S300613111(2)=2*Cij234(3,2)
       S300613111(3)=2*Cij234(2,2)
       S300613211(1)=2*Cij234(3,2)
       S300613211(2)=-2*Cij234(3,3)
       S300613211(3)=-2*Cij234(4,3)
       S300613221(1)=-2*Cij234(3,3)
       S300613221(2)=2*Cij234(3,4)
       S300613221(3)=2*Cij234(4,4)
       S300613222(1)=2*Cij234(3,4)
       S300613222(2)=-2*Cij234(3,5)
       S300613222(3)=-2*Cij234(4,5)
       S300613311(1)=2*Cij234(2,2)
       S300613311(2)=-2*Cij234(4,3)
       S300613311(3)=-2*Cij234(2,3)
       S300613321(1)=-2*Cij234(4,3)
       S300613321(2)=2*Cij234(4,4)
       S300613321(3)=2*Cij234(5,4)
       S300613322(1)=2*Cij234(4,4)
       S300613322(2)=-2*Cij234(4,5)
       S300613322(3)=-2*Cij234(5,5)
       S300613331(1)=-2*Cij234(2,3)
       S300613331(2)=2*Cij234(5,4)
       S300613331(3)=2*Cij234(2,4)
       S300613332(1)=2*Cij234(5,4)
       S300613332(2)=-2*Cij234(5,5)
       S300613332(3)=-2*Cij234(6,5)
       S300613333(1)=2*Cij234(2,4)
       S300613333(2)=-2*Cij234(6,5)
       S300613333(3)=-2*Cij234(2,5)
       S300621111(1)=-2*Cij234(1,1)
       S300621111(2)=2*Cij234(1,2)
       S300621111(3)=2*Cij234(3,2)
       S300622111(1)=2*Cij234(1,2)
       S300622111(2)=-2*Cij234(1,3)
       S300622111(3)=-2*Cij234(3,3)
       S300622211(1)=-2*Cij234(1,3)
       S300622211(2)=2*Cij234(1,4)
       S300622211(3)=2*Cij234(3,4)
       S300622221(1)=2*Cij234(1,4)
       S300622221(2)=-2*Cij234(1,5)
       S300622221(3)=-2*Cij234(3,5)
       S300622222(1)=-2*Cij234(1,5)
       S300622222(2)=2*Cij234(1,6)
       S300622222(3)=2*Cij234(3,6)
       S300623111(1)=2*Cij234(3,2)
       S300623111(2)=-2*Cij234(3,3)
       S300623111(3)=-2*Cij234(4,3)
       S300623211(1)=-2*Cij234(3,3)
       S300623211(2)=2*Cij234(3,4)
       S300623211(3)=2*Cij234(4,4)
       S300623221(1)=2*Cij234(3,4)
       S300623221(2)=-2*Cij234(3,5)
       S300623221(3)=-2*Cij234(4,5)
       S300623222(1)=-2*Cij234(3,5)
       S300623222(2)=2*Cij234(3,6)
       S300623222(3)=2*Cij234(4,6)
       S300623311(1)=-2*Cij234(4,3)
       S300623311(2)=2*Cij234(4,4)
       S300623311(3)=2*Cij234(5,4)
       S300623321(1)=2*Cij234(4,4)
       S300623321(2)=-2*Cij234(4,5)
       S300623321(3)=-2*Cij234(5,5)
       S300623322(1)=-2*Cij234(4,5)
       S300623322(2)=2*Cij234(4,6)
       S300623322(3)=2*Cij234(5,6)
       S300623331(1)=2*Cij234(5,4)
       S300623331(2)=-2*Cij234(5,5)
       S300623331(3)=-2*Cij234(6,5)
       S300623332(1)=-2*Cij234(5,5)
       S300623332(2)=2*Cij234(5,6)
       S300623332(3)=2*Cij234(6,6)
       S300623333(1)=-2*Cij234(6,5)
       S300623333(2)=2*Cij234(6,6)
       S300623333(3)=2*Cij234(7,6)
       S300631111(1)=-2*Cij234(2,1)
       S300631111(2)=2*Cij234(3,2)
       S300631111(3)=2*Cij234(2,2)
       S300632111(1)=2*Cij234(3,2)
       S300632111(2)=-2*Cij234(3,3)
       S300632111(3)=-2*Cij234(4,3)
       S300632211(1)=-2*Cij234(3,3)
       S300632211(2)=2*Cij234(3,4)
       S300632211(3)=2*Cij234(4,4)
       S300632221(1)=2*Cij234(3,4)
       S300632221(2)=-2*Cij234(3,5)
       S300632221(3)=-2*Cij234(4,5)
       S300632222(1)=-2*Cij234(3,5)
       S300632222(2)=2*Cij234(3,6)
       S300632222(3)=2*Cij234(4,6)
       S300633111(1)=2*Cij234(2,2)
       S300633111(2)=-2*Cij234(4,3)
       S300633111(3)=-2*Cij234(2,3)
       S300633211(1)=-2*Cij234(4,3)
       S300633211(2)=2*Cij234(4,4)
       S300633211(3)=2*Cij234(5,4)
       S300633221(1)=2*Cij234(4,4)
       S300633221(2)=-2*Cij234(4,5)
       S300633221(3)=-2*Cij234(5,5)
       S300633222(1)=-2*Cij234(4,5)
       S300633222(2)=2*Cij234(4,6)
       S300633222(3)=2*Cij234(5,6)
       S300633311(1)=-2*Cij234(2,3)
       S300633311(2)=2*Cij234(5,4)
       S300633311(3)=2*Cij234(2,4)
       S300633321(1)=2*Cij234(5,4)
       S300633321(2)=-2*Cij234(5,5)
       S300633321(3)=-2*Cij234(6,5)
       S300633322(1)=-2*Cij234(5,5)
       S300633322(2)=2*Cij234(5,6)
       S300633322(3)=2*Cij234(6,6)
       S300633331(1)=2*Cij234(2,4)
       S300633331(2)=-2*Cij234(6,5)
       S300633331(3)=-2*Cij234(2,5)
       S300633332(1)=-2*Cij234(6,5)
       S300633332(2)=2*Cij234(6,6)
       S300633332(3)=2*Cij234(7,6)
       S300633333(1)=-2*Cij234(2,5)
       S300633333(2)=2*Cij234(7,6)
       S300633333(3)=2*Cij234(2,6)
       S3h00611111(1)=Cij134(9,7)+Cij234(4,2)
       S3h00611111(2)=Cij134(9,7)-Cij234(5,3)
       S3h00611111(3)=Cij134(11,7)-Cij234(6,3)
       S3h00612111(1)=Cij134(9,7)-Cij234(5,3)
       S3h00612111(2)=Cij134(9,7)+Cij234(6,4)
       S3h00612111(3)=Cij134(11,7)+Cij234(8,4)
       S3h00612211(1)=Cij134(9,7)+Cij234(6,4)
       S3h00612211(2)=Cij134(9,7)-Cij234(7,5)
       S3h00612211(3)=Cij134(11,7)-Cij234(9,5)
       S3h00612221(1)=Cij134(9,7)-Cij234(7,5)
       S3h00612221(2)=Cij134(9,7)+Cij234(8,6)
       S3h00612221(3)=Cij134(11,7)+Cij234(10,6)
       S3h00612222(1)=Cij134(9,7)+Cij234(8,6)
       S3h00612222(2)=Cij134(9,7)-Cij234(9,7)
       S3h00612222(3)=Cij134(11,7)-Cij234(11,7)
       S3h00613111(1)=Cij134(11,7)-Cij234(6,3)
       S3h00613111(2)=Cij134(11,7)+Cij234(8,4)
       S3h00613111(3)=Cij134(12,7)+Cij234(7,4)
       S3h00613211(1)=Cij134(11,7)+Cij234(8,4)
       S3h00613211(2)=Cij134(11,7)-Cij234(9,5)
       S3h00613211(3)=Cij134(12,7)-Cij234(10,5)
       S3h00613221(1)=Cij134(11,7)-Cij234(9,5)
       S3h00613221(2)=Cij134(11,7)+Cij234(10,6)
       S3h00613221(3)=Cij134(12,7)+Cij234(11,6)
       S3h00613222(1)=Cij134(11,7)+Cij234(10,6)
       S3h00613222(2)=Cij134(11,7)-Cij234(11,7)
       S3h00613222(3)=Cij134(12,7)-Cij234(12,7)
       S3h00613311(1)=Cij134(12,7)+Cij234(7,4)
       S3h00613311(2)=Cij134(12,7)-Cij234(10,5)
       S3h00613311(3)=Cij134(13,7)-Cij234(8,5)
       S3h00613321(1)=Cij134(12,7)-Cij234(10,5)
       S3h00613321(2)=Cij134(12,7)+Cij234(11,6)
       S3h00613321(3)=Cij134(13,7)+Cij234(12,6)
       S3h00613322(1)=Cij134(12,7)+Cij234(11,6)
       S3h00613322(2)=Cij134(12,7)-Cij234(12,7)
       S3h00613322(3)=Cij134(13,7)-Cij234(13,7)
       S3h00613331(1)=Cij134(13,7)-Cij234(8,5)
       S3h00613331(2)=Cij134(13,7)+Cij234(12,6)
       S3h00613331(3)=Cij134(14,7)+Cij234(9,6)
       S3h00613332(1)=Cij134(13,7)+Cij234(12,6)
       S3h00613332(2)=Cij134(13,7)-Cij234(13,7)
       S3h00613332(3)=Cij134(14,7)-Cij234(14,7)
       S3h00613333(1)=Cij134(14,7)+Cij234(9,6)
       S3h00613333(2)=Cij134(14,7)-Cij234(14,7)
       S3h00613333(3)=Cij134(10,7)-Cij234(10,7)
       S3h00621111(1)=Cij124(9,7)-Cij134(9,7)
       S3h00621111(2)=Cij124(11,7)-Cij134(9,7)
       S3h00621111(3)=Cij124(11,7)-Cij134(11,7)
       S3h00622111(1)=Cij124(11,7)-Cij134(9,7)
       S3h00622111(2)=Cij124(12,7)-Cij134(9,7)
       S3h00622111(3)=Cij124(12,7)-Cij134(11,7)
       S3h00622211(1)=Cij124(12,7)-Cij134(9,7)
       S3h00622211(2)=Cij124(13,7)-Cij134(9,7)
       S3h00622211(3)=Cij124(13,7)-Cij134(11,7)
       S3h00622221(1)=Cij124(13,7)-Cij134(9,7)
       S3h00622221(2)=Cij124(14,7)-Cij134(9,7)
       S3h00622221(3)=Cij124(14,7)-Cij134(11,7)
       S3h00622222(1)=Cij124(14,7)-Cij134(9,7)
       S3h00622222(2)=Cij124(10,7)-Cij134(9,7)
       S3h00622222(3)=Cij124(10,7)-Cij134(11,7)
       S3h00623111(1)=Cij124(11,7)-Cij134(11,7)
       S3h00623111(2)=Cij124(12,7)-Cij134(11,7)
       S3h00623111(3)=Cij124(12,7)-Cij134(12,7)
       S3h00623211(1)=Cij124(12,7)-Cij134(11,7)
       S3h00623211(2)=Cij124(13,7)-Cij134(11,7)
       S3h00623211(3)=Cij124(13,7)-Cij134(12,7)
       S3h00623221(1)=Cij124(13,7)-Cij134(11,7)
       S3h00623221(2)=Cij124(14,7)-Cij134(11,7)
       S3h00623221(3)=Cij124(14,7)-Cij134(12,7)
       S3h00623222(1)=Cij124(14,7)-Cij134(11,7)
       S3h00623222(2)=Cij124(10,7)-Cij134(11,7)
       S3h00623222(3)=Cij124(10,7)-Cij134(12,7)
       S3h00623311(1)=Cij124(12,7)-Cij134(12,7)
       S3h00623311(2)=Cij124(13,7)-Cij134(12,7)
       S3h00623311(3)=Cij124(13,7)-Cij134(13,7)
       S3h00623321(1)=Cij124(13,7)-Cij134(12,7)
       S3h00623321(2)=Cij124(14,7)-Cij134(12,7)
       S3h00623321(3)=Cij124(14,7)-Cij134(13,7)
       S3h00623322(1)=Cij124(14,7)-Cij134(12,7)
       S3h00623322(2)=Cij124(10,7)-Cij134(12,7)
       S3h00623322(3)=Cij124(10,7)-Cij134(13,7)
       S3h00623331(1)=Cij124(13,7)-Cij134(13,7)
       S3h00623331(2)=Cij124(14,7)-Cij134(13,7)
       S3h00623331(3)=Cij124(14,7)-Cij134(14,7)
       S3h00623332(1)=Cij124(14,7)-Cij134(13,7)
       S3h00623332(2)=Cij124(10,7)-Cij134(13,7)
       S3h00623332(3)=Cij124(10,7)-Cij134(14,7)
       S3h00623333(1)=Cij124(14,7)-Cij134(14,7)
       S3h00623333(2)=Cij124(10,7)-Cij134(14,7)
       S3h00623333(3)=Cij124(10,7)-Cij134(10,7)
       S3h00631111(1)=Cij123(9,7)-Cij124(9,7)
       S3h00631111(2)=Cij123(11,7)-Cij124(11,7)
       S3h00631111(3)=-Cij124(11,7)
       S3h00632111(1)=Cij123(11,7)-Cij124(11,7)
       S3h00632111(2)=Cij123(12,7)-Cij124(12,7)
       S3h00632111(3)=-Cij124(12,7)
       S3h00632211(1)=Cij123(12,7)-Cij124(12,7)
       S3h00632211(2)=Cij123(13,7)-Cij124(13,7)
       S3h00632211(3)=-Cij124(13,7)
       S3h00632221(1)=Cij123(13,7)-Cij124(13,7)
       S3h00632221(2)=Cij123(14,7)-Cij124(14,7)
       S3h00632221(3)=-Cij124(14,7)
       S3h00632222(1)=Cij123(14,7)-Cij124(14,7)
       S3h00632222(2)=Cij123(10,7)-Cij124(10,7)
       S3h00632222(3)=-Cij124(10,7)
       S3h00633111(1)=-Cij124(11,7)
       S3h00633111(2)=-Cij124(12,7)
       S3h00633111(3)=-Cij124(12,7)
       S3h00633211(1)=-Cij124(12,7)
       S3h00633211(2)=-Cij124(13,7)
       S3h00633211(3)=-Cij124(13,7)
       S3h00633221(1)=-Cij124(13,7)
       S3h00633221(2)=-Cij124(14,7)
       S3h00633221(3)=-Cij124(14,7)
       S3h00633222(1)=-Cij124(14,7)
       S3h00633222(2)=-Cij124(10,7)
       S3h00633222(3)=-Cij124(10,7)
       S3h00633311(1)=-Cij124(12,7)
       S3h00633311(2)=-Cij124(13,7)
       S3h00633311(3)=-Cij124(13,7)
       S3h00633321(1)=-Cij124(13,7)
       S3h00633321(2)=-Cij124(14,7)
       S3h00633321(3)=-Cij124(14,7)
       S3h00633322(1)=-Cij124(14,7)
       S3h00633322(2)=-Cij124(10,7)
       S3h00633322(3)=-Cij124(10,7)
       S3h00633331(1)=-Cij124(13,7)
       S3h00633331(2)=-Cij124(14,7)
       S3h00633331(3)=-Cij124(14,7)
       S3h00633332(1)=-Cij124(14,7)
       S3h00633332(2)=-Cij124(10,7)
       S3h00633332(3)=-Cij124(10,7)
       S3h00633333(1)=-Cij124(14,7)
       S3h00633333(2)=-Cij124(10,7)
       S3h00633333(3)=-Cij124(10,7)
       aux000041(1,1,1)=-(F(1)*S3h0051111(1))-F(2)*S3h0052111(1)-F(3)*S3
     &  h0053111(1)+S3h00611111(k)*Z(1,l)+S3h00621111(k)*Z(2,l)+S3h00631
     &  111(k)*Z(3,l)+(Inv401111+S300004111(1)-S3h00611111(1)-S3h0062211
     &  1(1)-S3h00633111(1))*Z(k,l)-8*S3h00004111(1)*ZZ(k,1,l,1)-8*S3h00
     &  004211(1)*ZZ(k,1,l,2)-8*S3h00004311(1)*ZZ(k,1,l,3)
       aux000042(1,1,1)=-(F(1)*S3h0051211(1))-F(2)*S3h0052211(1)-F(3)*S3
     &  h0053211(1)+S3h00612111(k)*Z(1,l)+S3h00622111(k)*Z(2,l)+S3h00632
     &  111(k)*Z(3,l)+(Inv402111+S300004211(1)-S3h00612111(1)-S3h0062221
     &  1(1)-S3h00633211(1))*Z(k,l)-6*S3h00004121(1)*ZZ(k,1,l,1)-6*S3h00
     &  004221(1)*ZZ(k,1,l,2)-6*S3h00004321(1)*ZZ(k,1,l,3)-2*S3h00004111
     &  (1)*ZZ(k,2,l,1)-2*S3h00004211(1)*ZZ(k,2,l,2)-2*S3h00004311(1)*ZZ
     &  (k,2,l,3)
       aux000042(2,1,1)=-(F(1)*S3h0051221(1))-F(2)*S3h0052221(1)-F(3)*S3
     &  h0053221(1)+S3h00612211(k)*Z(1,l)+S3h00622211(k)*Z(2,l)+S3h00632
     &  211(k)*Z(3,l)+(Inv402211+S300004221(1)-S3h00612211(1)-S3h0062222
     &  1(1)-S3h00633221(1))*Z(k,l)-4*S3h00004122(1)*ZZ(k,1,l,1)-4*S3h00
     &  004222(1)*ZZ(k,1,l,2)-4*S3h00004322(1)*ZZ(k,1,l,3)-4*S3h00004121
     &  (1)*ZZ(k,2,l,1)-4*S3h00004221(1)*ZZ(k,2,l,2)-4*S3h00004321(1)*ZZ
     &  (k,2,l,3)
       aux000042(2,2,1)=-(F(1)*S3h0051222(1))-F(2)*S3h0052222(1)-F(3)*S3
     &  h0053222(1)+S3h00612221(k)*Z(1,l)+S3h00622221(k)*Z(2,l)+S3h00632
     &  221(k)*Z(3,l)+(Inv402221+S300004222(1)-S3h00612221(1)-S3h0062222
     &  2(1)-S3h00633222(1))*Z(k,l)-2*S3h00004122(2)*ZZ(k,1,l,1)-2*S3h00
     &  004222(2)*ZZ(k,1,l,2)-2*S3h00004322(2)*ZZ(k,1,l,3)-6*S3h00004122
     &  (1)*ZZ(k,2,l,1)-6*S3h00004222(1)*ZZ(k,2,l,2)-6*S3h00004322(1)*ZZ
     &  (k,2,l,3)
       aux000042(2,2,2)=-(F(1)*S3h0051222(2))-F(2)*S3h0052222(2)-F(3)*S3
     &  h0053222(2)+S3h00612222(k)*Z(1,l)+S3h00622222(k)*Z(2,l)+S3h00632
     &  222(k)*Z(3,l)+(Inv402222+S300004222(2)-S3h00612222(1)-S3h0062222
     &  2(2)-S3h00633222(2))*Z(k,l)-8*S3h00004122(2)*ZZ(k,2,l,1)-8*S3h00
     &  004222(2)*ZZ(k,2,l,2)-8*S3h00004322(2)*ZZ(k,2,l,3)
       aux000043(1,1,1)=-(F(1)*S3h0051311(1))-F(2)*S3h0052311(1)-F(3)*S3
     &  h0053311(1)+S3h00613111(k)*Z(1,l)+S3h00623111(k)*Z(2,l)+S3h00633
     &  111(k)*Z(3,l)+(Inv403111+S300004311(1)-S3h00613111(1)-S3h0062321
     &  1(1)-S3h00633311(1))*Z(k,l)-6*S3h00004131(1)*ZZ(k,1,l,1)-6*S3h00
     &  004231(1)*ZZ(k,1,l,2)-6*S3h00004331(1)*ZZ(k,1,l,3)-2*S3h00004111
     &  (1)*ZZ(k,3,l,1)-2*S3h00004211(1)*ZZ(k,3,l,2)-2*S3h00004311(1)*ZZ
     &  (k,3,l,3)
       aux000043(2,1,1)=-(F(1)*S3h0051321(1))-F(2)*S3h0052321(1)-F(3)*S3
     &  h0053321(1)+S3h00613211(k)*Z(1,l)+S3h00623211(k)*Z(2,l)+S3h00633
     &  211(k)*Z(3,l)+(Inv403211+S300004321(1)-S3h00613211(1)-S3h0062322
     &  1(1)-S3h00633321(1))*Z(k,l)-4*S3h00004132(1)*ZZ(k,1,l,1)-4*S3h00
     &  004232(1)*ZZ(k,1,l,2)-4*S3h00004332(1)*ZZ(k,1,l,3)-2*S3h00004131
     &  (1)*ZZ(k,2,l,1)-2*S3h00004231(1)*ZZ(k,2,l,2)-2*S3h00004331(1)*ZZ
     &  (k,2,l,3)-2*S3h00004121(1)*ZZ(k,3,l,1)-2*S3h00004221(1)*ZZ(k,3,l
     &  ,2)-2*S3h00004321(1)*ZZ(k,3,l,3)
       aux000043(2,2,1)=-(F(1)*S3h0051322(1))-F(2)*S3h0052322(1)-F(3)*S3
     &  h0053322(1)+S3h00613221(k)*Z(1,l)+S3h00623221(k)*Z(2,l)+S3h00633
     &  221(k)*Z(3,l)+(Inv403221+S300004322(1)-S3h00613221(1)-S3h0062322
     &  2(1)-S3h00633322(1))*Z(k,l)-2*S3h00004132(2)*ZZ(k,1,l,1)-2*S3h00
     &  004232(2)*ZZ(k,1,l,2)-2*S3h00004332(2)*ZZ(k,1,l,3)-4*S3h00004132
     &  (1)*ZZ(k,2,l,1)-4*S3h00004232(1)*ZZ(k,2,l,2)-4*S3h00004332(1)*ZZ
     &  (k,2,l,3)-2*S3h00004122(1)*ZZ(k,3,l,1)-2*S3h00004222(1)*ZZ(k,3,l
     &  ,2)-2*S3h00004322(1)*ZZ(k,3,l,3)
       aux000043(2,2,2)=-(F(1)*S3h0051322(2))-F(2)*S3h0052322(2)-F(3)*S3
     &  h0053322(2)+S3h00613222(k)*Z(1,l)+S3h00623222(k)*Z(2,l)+S3h00633
     &  222(k)*Z(3,l)+(Inv403222+S300004322(2)-S3h00613222(1)-S3h0062322
     &  2(2)-S3h00633322(2))*Z(k,l)-6*S3h00004132(2)*ZZ(k,2,l,1)-6*S3h00
     &  004232(2)*ZZ(k,2,l,2)-6*S3h00004332(2)*ZZ(k,2,l,3)-2*S3h00004122
     &  (2)*ZZ(k,3,l,1)-2*S3h00004222(2)*ZZ(k,3,l,2)-2*S3h00004322(2)*ZZ
     &  (k,3,l,3)
       aux000043(3,1,1)=-(F(1)*S3h0051331(1))-F(2)*S3h0052331(1)-F(3)*S3
     &  h0053331(1)+S3h00613311(k)*Z(1,l)+S3h00623311(k)*Z(2,l)+S3h00633
     &  311(k)*Z(3,l)+(Inv403311+S300004331(1)-S3h00613311(1)-S3h0062332
     &  1(1)-S3h00633331(1))*Z(k,l)-4*S3h00004133(1)*ZZ(k,1,l,1)-4*S3h00
     &  004233(1)*ZZ(k,1,l,2)-4*S3h00004333(1)*ZZ(k,1,l,3)-4*S3h00004131
     &  (1)*ZZ(k,3,l,1)-4*S3h00004231(1)*ZZ(k,3,l,2)-4*S3h00004331(1)*ZZ
     &  (k,3,l,3)
       aux000043(3,2,1)=-(F(1)*S3h0051332(1))-F(2)*S3h0052332(1)-F(3)*S3
     &  h0053332(1)+S3h00613321(k)*Z(1,l)+S3h00623321(k)*Z(2,l)+S3h00633
     &  321(k)*Z(3,l)+(Inv403321+S300004332(1)-S3h00613321(1)-S3h0062332
     &  2(1)-S3h00633332(1))*Z(k,l)-2*S3h00004133(2)*ZZ(k,1,l,1)-2*S3h00
     &  004233(2)*ZZ(k,1,l,2)-2*S3h00004333(2)*ZZ(k,1,l,3)-2*S3h00004133
     &  (1)*ZZ(k,2,l,1)-2*S3h00004233(1)*ZZ(k,2,l,2)-2*S3h00004333(1)*ZZ
     &  (k,2,l,3)-4*S3h00004132(1)*ZZ(k,3,l,1)-4*S3h00004232(1)*ZZ(k,3,l
     &  ,2)-4*S3h00004332(1)*ZZ(k,3,l,3)
       aux000043(3,2,2)=-(F(1)*S3h0051332(2))-F(2)*S3h0052332(2)-F(3)*S3
     &  h0053332(2)+S3h00613322(k)*Z(1,l)+S3h00623322(k)*Z(2,l)+S3h00633
     &  322(k)*Z(3,l)+(Inv403322+S300004332(2)-S3h00613322(1)-S3h0062332
     &  2(2)-S3h00633332(2))*Z(k,l)-4*S3h00004133(2)*ZZ(k,2,l,1)-4*S3h00
     &  004233(2)*ZZ(k,2,l,2)-4*S3h00004333(2)*ZZ(k,2,l,3)-4*S3h00004132
     &  (2)*ZZ(k,3,l,1)-4*S3h00004232(2)*ZZ(k,3,l,2)-4*S3h00004332(2)*ZZ
     &  (k,3,l,3)
       aux000043(3,3,1)=-(F(1)*S3h0051333(1))-F(2)*S3h0052333(1)-F(3)*S3
     &  h0053333(1)+S3h00613331(k)*Z(1,l)+S3h00623331(k)*Z(2,l)+S3h00633
     &  331(k)*Z(3,l)+(Inv403331+S300004333(1)-S3h00613331(1)-S3h0062333
     &  2(1)-S3h00633333(1))*Z(k,l)-2*S3h00004133(3)*ZZ(k,1,l,1)-2*S3h00
     &  004233(3)*ZZ(k,1,l,2)-2*S3h00004333(3)*ZZ(k,1,l,3)-6*S3h00004133
     &  (1)*ZZ(k,3,l,1)-6*S3h00004233(1)*ZZ(k,3,l,2)-6*S3h00004333(1)*ZZ
     &  (k,3,l,3)
       aux000043(3,3,2)=-(F(1)*S3h0051333(2))-F(2)*S3h0052333(2)-F(3)*S3
     &  h0053333(2)+S3h00613332(k)*Z(1,l)+S3h00623332(k)*Z(2,l)+S3h00633
     &  332(k)*Z(3,l)+(Inv403332+S300004333(2)-S3h00613332(1)-S3h0062333
     &  2(2)-S3h00633333(2))*Z(k,l)-2*S3h00004133(3)*ZZ(k,2,l,1)-2*S3h00
     &  004233(3)*ZZ(k,2,l,2)-2*S3h00004333(3)*ZZ(k,2,l,3)-6*S3h00004133
     &  (2)*ZZ(k,3,l,1)-6*S3h00004233(2)*ZZ(k,3,l,2)-6*S3h00004333(2)*ZZ
     &  (k,3,l,3)
       aux000043(3,3,3)=-(F(1)*S3h0051333(3))-F(2)*S3h0052333(3)-F(3)*S3
     &  h0053333(3)+S3h00613333(k)*Z(1,l)+S3h00623333(k)*Z(2,l)+S3h00633
     &  333(k)*Z(3,l)+(Inv403333+S300004333(3)-S3h00613333(1)-S3h0062333
     &  3(2)-S3h00633333(3))*Z(k,l)-8*S3h00004133(3)*ZZ(k,3,l,1)-8*S3h00
     &  004233(3)*ZZ(k,3,l,2)-8*S3h00004333(3)*ZZ(k,3,l,3)
       temp000041(1,1,1)=I24Z*(aux000041(1,1,1)+8*F(4)*temp00003(1,1,1)+
     &  F(5)*temp0041(1,1,1)+48*temp0000002(1,1)*ZZ(k,1,l,1))
       temp000042(1,1,1)=I24Z*(aux000042(1,1,1)+2*F(6)*temp00003(1,1,1)+
     &  6*F(4)*temp00003(2,1,1)+F(5)*temp0042(1,1,1)+24*temp0000002(2,1)
     &  *ZZ(k,1,l,1)+12*temp0000002(1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp000042(2,1,1)=I24Z*(aux000042(2,1,1)+4*F(6)*temp00003(2,1,1)+
     &  4*F(4)*temp00003(2,2,1)+F(5)*temp0042(2,1,1)+16*temp0000002(2,1)
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp0000002(2,2)*ZZ(k,1,l,1)+temp0
     &  000002(1,1)*ZZ(k,2,l,2)))
       temp000042(2,2,1)=I24Z*(aux000042(2,2,1)+6*F(6)*temp00003(2,2,1)+
     &  2*F(4)*temp00003(2,2,2)+F(5)*temp0042(2,2,1)+12*temp0000002(2,2)
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*temp0000002(2,1)*ZZ(k,2,l,2))
       temp000042(2,2,2)=I24Z*(aux000042(2,2,2)+8*F(6)*temp00003(2,2,2)+
     &  F(5)*temp0042(2,2,2)+48*temp0000002(2,2)*ZZ(k,2,l,2))
       temp000043(1,1,1)=I24Z*(aux000043(1,1,1)+2*F(7)*temp00003(1,1,1)+
     &  6*F(4)*temp00003(3,1,1)+F(5)*temp0043(1,1,1)+24*temp0000002(3,1)
     &  *ZZ(k,1,l,1)+12*temp0000002(1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)))
       temp000043(2,1,1)=I24Z*(aux000043(2,1,1)+2*F(7)*temp00003(2,1,1)+
     &  2*F(6)*temp00003(3,1,1)+4*F(4)*temp00003(3,2,1)+F(5)*temp0043(2,
     &  1,1)+8*(temp0000002(3,2)*ZZ(k,1,l,1)+temp0000002(3,1)*(ZZ(k,1,l,
     &  2)+ZZ(k,2,l,1)))+8*temp0000002(2,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+4*
     &  temp0000002(1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp000043(2,2,1)=I24Z*(aux000043(2,2,1)+2*F(7)*temp00003(2,2,1)+
     &  4*F(6)*temp00003(3,2,1)+2*F(4)*temp00003(3,2,2)+F(5)*temp0043(2,
     &  2,1)+8*(temp0000002(3,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0000002(3
     &  ,1)*ZZ(k,2,l,2))+4*temp0000002(2,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*
     &  temp0000002(2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp000043(2,2,2)=I24Z*(aux000043(2,2,2)+2*F(7)*temp00003(2,2,2)+
     &  6*F(6)*temp00003(3,2,2)+F(5)*temp0043(2,2,2)+24*temp0000002(3,2)
     &  *ZZ(k,2,l,2)+12*temp0000002(2,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp000043(3,1,1)=I24Z*(aux000043(3,1,1)+4*F(7)*temp00003(3,1,1)+
     &  4*F(4)*temp00003(3,3,1)+F(5)*temp0043(3,1,1)+16*temp0000002(3,1)
     &  *(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*(temp0000002(3,3)*ZZ(k,1,l,1)+temp0
     &  000002(1,1)*ZZ(k,3,l,3)))
       temp000043(3,2,1)=I24Z*(aux000043(3,2,1)+2*F(6)*temp00003(3,3,1)+
     &  2*F(4)*temp00003(3,3,2)+F(5)*temp0043(3,2,1)+4*(F(7)*temp00003(3
     &  ,2,1)+temp0000002(3,3)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))+8*temp0000002(
     &  3,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*temp0000002(3,1)*(ZZ(k,2,l,3)+Z
     &  Z(k,3,l,2))+8*temp0000002(2,1)*ZZ(k,3,l,3))
       temp000043(3,2,2)=I24Z*(aux000043(3,2,2)+4*F(7)*temp00003(3,2,2)+
     &  4*F(6)*temp00003(3,3,2)+F(5)*temp0043(3,2,2)+16*temp0000002(3,2)
     &  *(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*(temp0000002(3,3)*ZZ(k,2,l,2)+temp0
     &  000002(2,2)*ZZ(k,3,l,3)))
       temp000043(3,3,1)=I24Z*(aux000043(3,3,1)+6*F(7)*temp00003(3,3,1)+
     &  2*F(4)*temp00003(3,3,3)+F(5)*temp0043(3,3,1)+12*temp0000002(3,3)
     &  *(ZZ(k,1,l,3)+ZZ(k,3,l,1))+24*temp0000002(3,1)*ZZ(k,3,l,3))
       temp000043(3,3,2)=I24Z*(aux000043(3,3,2)+6*F(7)*temp00003(3,3,2)+
     &  2*F(6)*temp00003(3,3,3)+F(5)*temp0043(3,3,2)+12*temp0000002(3,3)
     &  *(ZZ(k,2,l,3)+ZZ(k,3,l,2))+24*temp0000002(3,2)*ZZ(k,3,l,3))
       temp000043(3,3,3)=I24Z*(aux000043(3,3,3)+8*F(7)*temp00003(3,3,3)+
     &  F(5)*temp0043(3,3,3)+48*temp0000002(3,3)*ZZ(k,3,l,3))
       temp000041(1,1,2)=temp000042(1,1,1)
       temp000041(1,1,3)=temp000043(1,1,1)
       temp000041(1,2,1)=temp000042(1,1,1)
       temp000041(1,2,2)=temp000042(2,1,1)
       temp000041(1,2,3)=temp000043(2,1,1)
       temp000041(1,3,1)=temp000043(1,1,1)
       temp000041(1,3,2)=temp000043(2,1,1)
       temp000041(1,3,3)=temp000043(3,1,1)
       temp000042(1,1,2)=temp000042(2,1,1)
       temp000042(1,1,3)=temp000043(2,1,1)
       temp000042(1,2,1)=temp000042(2,1,1)
       temp000042(1,2,2)=temp000042(2,2,1)
       temp000042(1,2,3)=temp000043(2,2,1)
       temp000042(1,3,1)=temp000043(2,1,1)
       temp000042(1,3,2)=temp000043(2,2,1)
       temp000042(1,3,3)=temp000043(3,2,1)
       temp000042(2,1,2)=temp000042(2,2,1)
       temp000042(2,1,3)=temp000043(2,2,1)
       temp000042(2,2,3)=temp000043(2,2,2)
       temp000042(2,3,1)=temp000043(2,2,1)
       temp000042(2,3,2)=temp000043(2,2,2)
       temp000042(2,3,3)=temp000043(3,2,2)
       temp000043(1,1,2)=temp000043(2,1,1)
       temp000043(1,1,3)=temp000043(3,1,1)
       temp000043(1,2,1)=temp000043(2,1,1)
       temp000043(1,2,2)=temp000043(2,2,1)
       temp000043(1,2,3)=temp000043(3,2,1)
       temp000043(1,3,1)=temp000043(3,1,1)
       temp000043(1,3,2)=temp000043(3,2,1)
       temp000043(1,3,3)=temp000043(3,3,1)
       temp000043(2,1,2)=temp000043(2,2,1)
       temp000043(2,1,3)=temp000043(3,2,1)
       temp000043(2,2,3)=temp000043(3,2,2)
       temp000043(2,3,1)=temp000043(3,2,1)
       temp000043(2,3,2)=temp000043(3,2,2)
       temp000043(2,3,3)=temp000043(3,3,2)
       temp000043(3,1,2)=temp000043(3,2,1)
       temp000043(3,1,3)=temp000043(3,3,1)
       temp000043(3,2,3)=temp000043(3,3,2)
       S381111111(1)=C0234+Cij134(1,7)
       S381111111(2)=Cij134(1,7)-Cij234(1,1)
       S381111111(3)=Cij134(3,7)-Cij234(2,1)
       S381211111(1)=Cij134(1,7)-Cij234(1,1)
       S381211111(2)=Cij134(1,7)+Cij234(1,2)
       S381211111(3)=Cij134(3,7)+Cij234(3,2)
       S381221111(1)=Cij134(1,7)+Cij234(1,2)
       S381221111(2)=Cij134(1,7)-Cij234(1,3)
       S381221111(3)=Cij134(3,7)-Cij234(3,3)
       S381222111(1)=Cij134(1,7)-Cij234(1,3)
       S381222111(2)=Cij134(1,7)+Cij234(1,4)
       S381222111(3)=Cij134(3,7)+Cij234(3,4)
       S381222211(1)=Cij134(1,7)+Cij234(1,4)
       S381222211(2)=Cij134(1,7)-Cij234(1,5)
       S381222211(3)=Cij134(3,7)-Cij234(3,5)
       S381222221(1)=Cij134(1,7)-Cij234(1,5)
       S381222221(2)=Cij134(1,7)+Cij234(1,6)
       S381222221(3)=Cij134(3,7)+Cij234(3,6)
       S381222222(1)=Cij134(1,7)+Cij234(1,6)
       S381222222(2)=Cij134(1,7)-Cij234(1,7)
       S381222222(3)=Cij134(3,7)-Cij234(3,7)
       S381311111(1)=Cij134(3,7)-Cij234(2,1)
       S381311111(2)=Cij134(3,7)+Cij234(3,2)
       S381311111(3)=Cij134(4,7)+Cij234(2,2)
       S381321111(1)=Cij134(3,7)+Cij234(3,2)
       S381321111(2)=Cij134(3,7)-Cij234(3,3)
       S381321111(3)=Cij134(4,7)-Cij234(4,3)
       S381322111(1)=Cij134(3,7)-Cij234(3,3)
       S381322111(2)=Cij134(3,7)+Cij234(3,4)
       S381322111(3)=Cij134(4,7)+Cij234(4,4)
       S381322211(1)=Cij134(3,7)+Cij234(3,4)
       S381322211(2)=Cij134(3,7)-Cij234(3,5)
       S381322211(3)=Cij134(4,7)-Cij234(4,5)
       S381322221(1)=Cij134(3,7)-Cij234(3,5)
       S381322221(2)=Cij134(3,7)+Cij234(3,6)
       S381322221(3)=Cij134(4,7)+Cij234(4,6)
       S381322222(1)=Cij134(3,7)+Cij234(3,6)
       S381322222(2)=Cij134(3,7)-Cij234(3,7)
       S381322222(3)=Cij134(4,7)-Cij234(4,7)
       S381331111(1)=Cij134(4,7)+Cij234(2,2)
       S381331111(2)=Cij134(4,7)-Cij234(4,3)
       S381331111(3)=Cij134(5,7)-Cij234(2,3)
       S381332111(1)=Cij134(4,7)-Cij234(4,3)
       S381332111(2)=Cij134(4,7)+Cij234(4,4)
       S381332111(3)=Cij134(5,7)+Cij234(5,4)
       S381332211(1)=Cij134(4,7)+Cij234(4,4)
       S381332211(2)=Cij134(4,7)-Cij234(4,5)
       S381332211(3)=Cij134(5,7)-Cij234(5,5)
       S381332221(1)=Cij134(4,7)-Cij234(4,5)
       S381332221(2)=Cij134(4,7)+Cij234(4,6)
       S381332221(3)=Cij134(5,7)+Cij234(5,6)
       S381332222(1)=Cij134(4,7)+Cij234(4,6)
       S381332222(2)=Cij134(4,7)-Cij234(4,7)
       S381332222(3)=Cij134(5,7)-Cij234(5,7)
       S381333111(1)=Cij134(5,7)-Cij234(2,3)
       S381333111(2)=Cij134(5,7)+Cij234(5,4)
       S381333111(3)=Cij134(6,7)+Cij234(2,4)
       S381333211(1)=Cij134(5,7)+Cij234(5,4)
       S381333211(2)=Cij134(5,7)-Cij234(5,5)
       S381333211(3)=Cij134(6,7)-Cij234(6,5)
       S381333221(1)=Cij134(5,7)-Cij234(5,5)
       S381333221(2)=Cij134(5,7)+Cij234(5,6)
       S381333221(3)=Cij134(6,7)+Cij234(6,6)
       S381333222(1)=Cij134(5,7)+Cij234(5,6)
       S381333222(2)=Cij134(5,7)-Cij234(5,7)
       S381333222(3)=Cij134(6,7)-Cij234(6,7)
       S381333311(1)=Cij134(6,7)+Cij234(2,4)
       S381333311(2)=Cij134(6,7)-Cij234(6,5)
       S381333311(3)=Cij134(7,7)-Cij234(2,5)
       S381333321(1)=Cij134(6,7)-Cij234(6,5)
       S381333321(2)=Cij134(6,7)+Cij234(6,6)
       S381333321(3)=Cij134(7,7)+Cij234(7,6)
       S381333322(1)=Cij134(6,7)+Cij234(6,6)
       S381333322(2)=Cij134(6,7)-Cij234(6,7)
       S381333322(3)=Cij134(7,7)-Cij234(7,7)
       S381333331(1)=Cij134(7,7)-Cij234(2,5)
       S381333331(2)=Cij134(7,7)+Cij234(7,6)
       S381333331(3)=Cij134(8,7)+Cij234(2,6)
       S381333332(1)=Cij134(7,7)+Cij234(7,6)
       S381333332(2)=Cij134(7,7)-Cij234(7,7)
       S381333332(3)=Cij134(8,7)-Cij234(8,7)
       S381333333(1)=Cij134(8,7)+Cij234(2,6)
       S381333333(2)=Cij134(8,7)-Cij234(8,7)
       S381333333(3)=Cij134(2,7)-Cij234(2,7)
       S382111111(1)=Cij124(1,7)-Cij134(1,7)
       S382111111(2)=Cij124(3,7)-Cij134(1,7)
       S382111111(3)=Cij124(3,7)-Cij134(3,7)
       S382211111(1)=Cij124(3,7)-Cij134(1,7)
       S382211111(2)=Cij124(4,7)-Cij134(1,7)
       S382211111(3)=Cij124(4,7)-Cij134(3,7)
       S382221111(1)=Cij124(4,7)-Cij134(1,7)
       S382221111(2)=Cij124(5,7)-Cij134(1,7)
       S382221111(3)=Cij124(5,7)-Cij134(3,7)
       S382222111(1)=Cij124(5,7)-Cij134(1,7)
       S382222111(2)=Cij124(6,7)-Cij134(1,7)
       S382222111(3)=Cij124(6,7)-Cij134(3,7)
       S382222211(1)=Cij124(6,7)-Cij134(1,7)
       S382222211(2)=Cij124(7,7)-Cij134(1,7)
       S382222211(3)=Cij124(7,7)-Cij134(3,7)
       S382222221(1)=Cij124(7,7)-Cij134(1,7)
       S382222221(2)=Cij124(8,7)-Cij134(1,7)
       S382222221(3)=Cij124(8,7)-Cij134(3,7)
       S382222222(1)=Cij124(8,7)-Cij134(1,7)
       S382222222(2)=Cij124(2,7)-Cij134(1,7)
       S382222222(3)=Cij124(2,7)-Cij134(3,7)
       S382311111(1)=Cij124(3,7)-Cij134(3,7)
       S382311111(2)=Cij124(4,7)-Cij134(3,7)
       S382311111(3)=Cij124(4,7)-Cij134(4,7)
       S382321111(1)=Cij124(4,7)-Cij134(3,7)
       S382321111(2)=Cij124(5,7)-Cij134(3,7)
       S382321111(3)=Cij124(5,7)-Cij134(4,7)
       S382322111(1)=Cij124(5,7)-Cij134(3,7)
       S382322111(2)=Cij124(6,7)-Cij134(3,7)
       S382322111(3)=Cij124(6,7)-Cij134(4,7)
       S382322211(1)=Cij124(6,7)-Cij134(3,7)
       S382322211(2)=Cij124(7,7)-Cij134(3,7)
       S382322211(3)=Cij124(7,7)-Cij134(4,7)
       S382322221(1)=Cij124(7,7)-Cij134(3,7)
       S382322221(2)=Cij124(8,7)-Cij134(3,7)
       S382322221(3)=Cij124(8,7)-Cij134(4,7)
       S382322222(1)=Cij124(8,7)-Cij134(3,7)
       S382322222(2)=Cij124(2,7)-Cij134(3,7)
       S382322222(3)=Cij124(2,7)-Cij134(4,7)
       S382331111(1)=Cij124(4,7)-Cij134(4,7)
       S382331111(2)=Cij124(5,7)-Cij134(4,7)
       S382331111(3)=Cij124(5,7)-Cij134(5,7)
       S382332111(1)=Cij124(5,7)-Cij134(4,7)
       S382332111(2)=Cij124(6,7)-Cij134(4,7)
       S382332111(3)=Cij124(6,7)-Cij134(5,7)
       S382332211(1)=Cij124(6,7)-Cij134(4,7)
       S382332211(2)=Cij124(7,7)-Cij134(4,7)
       S382332211(3)=Cij124(7,7)-Cij134(5,7)
       S382332221(1)=Cij124(7,7)-Cij134(4,7)
       S382332221(2)=Cij124(8,7)-Cij134(4,7)
       S382332221(3)=Cij124(8,7)-Cij134(5,7)
       S382332222(1)=Cij124(8,7)-Cij134(4,7)
       S382332222(2)=Cij124(2,7)-Cij134(4,7)
       S382332222(3)=Cij124(2,7)-Cij134(5,7)
       S382333111(1)=Cij124(5,7)-Cij134(5,7)
       S382333111(2)=Cij124(6,7)-Cij134(5,7)
       S382333111(3)=Cij124(6,7)-Cij134(6,7)
       S382333211(1)=Cij124(6,7)-Cij134(5,7)
       S382333211(2)=Cij124(7,7)-Cij134(5,7)
       S382333211(3)=Cij124(7,7)-Cij134(6,7)
       S382333221(1)=Cij124(7,7)-Cij134(5,7)
       S382333221(2)=Cij124(8,7)-Cij134(5,7)
       S382333221(3)=Cij124(8,7)-Cij134(6,7)
       S382333222(1)=Cij124(8,7)-Cij134(5,7)
       S382333222(2)=Cij124(2,7)-Cij134(5,7)
       S382333222(3)=Cij124(2,7)-Cij134(6,7)
       S382333311(1)=Cij124(6,7)-Cij134(6,7)
       S382333311(2)=Cij124(7,7)-Cij134(6,7)
       S382333311(3)=Cij124(7,7)-Cij134(7,7)
       S382333321(1)=Cij124(7,7)-Cij134(6,7)
       S382333321(2)=Cij124(8,7)-Cij134(6,7)
       S382333321(3)=Cij124(8,7)-Cij134(7,7)
       S382333322(1)=Cij124(8,7)-Cij134(6,7)
       S382333322(2)=Cij124(2,7)-Cij134(6,7)
       S382333322(3)=Cij124(2,7)-Cij134(7,7)
       S382333331(1)=Cij124(7,7)-Cij134(7,7)
       S382333331(2)=Cij124(8,7)-Cij134(7,7)
       S382333331(3)=Cij124(8,7)-Cij134(8,7)
       S382333332(1)=Cij124(8,7)-Cij134(7,7)
       S382333332(2)=Cij124(2,7)-Cij134(7,7)
       S382333332(3)=Cij124(2,7)-Cij134(8,7)
       S382333333(1)=Cij124(8,7)-Cij134(8,7)
       S382333333(2)=Cij124(2,7)-Cij134(8,7)
       S382333333(3)=Cij124(2,7)-Cij134(2,7)
       S383111111(1)=Cij123(1,7)-Cij124(1,7)
       S383111111(2)=Cij123(3,7)-Cij124(3,7)
       S383111111(3)=-Cij124(3,7)
       S383211111(1)=Cij123(3,7)-Cij124(3,7)
       S383211111(2)=Cij123(4,7)-Cij124(4,7)
       S383211111(3)=-Cij124(4,7)
       S383221111(1)=Cij123(4,7)-Cij124(4,7)
       S383221111(2)=Cij123(5,7)-Cij124(5,7)
       S383221111(3)=-Cij124(5,7)
       S383222111(1)=Cij123(5,7)-Cij124(5,7)
       S383222111(2)=Cij123(6,7)-Cij124(6,7)
       S383222111(3)=-Cij124(6,7)
       S383222211(1)=Cij123(6,7)-Cij124(6,7)
       S383222211(2)=Cij123(7,7)-Cij124(7,7)
       S383222211(3)=-Cij124(7,7)
       S383222221(1)=Cij123(7,7)-Cij124(7,7)
       S383222221(2)=Cij123(8,7)-Cij124(8,7)
       S383222221(3)=-Cij124(8,7)
       S383222222(1)=Cij123(8,7)-Cij124(8,7)
       S383222222(2)=Cij123(2,7)-Cij124(2,7)
       S383222222(3)=-Cij124(2,7)
       S383311111(1)=-Cij124(3,7)
       S383311111(2)=-Cij124(4,7)
       S383311111(3)=-Cij124(4,7)
       S383321111(1)=-Cij124(4,7)
       S383321111(2)=-Cij124(5,7)
       S383321111(3)=-Cij124(5,7)
       S383322111(1)=-Cij124(5,7)
       S383322111(2)=-Cij124(6,7)
       S383322111(3)=-Cij124(6,7)
       S383322211(1)=-Cij124(6,7)
       S383322211(2)=-Cij124(7,7)
       S383322211(3)=-Cij124(7,7)
       S383322221(1)=-Cij124(7,7)
       S383322221(2)=-Cij124(8,7)
       S383322221(3)=-Cij124(8,7)
       S383322222(1)=-Cij124(8,7)
       S383322222(2)=-Cij124(2,7)
       S383322222(3)=-Cij124(2,7)
       S383331111(1)=-Cij124(4,7)
       S383331111(2)=-Cij124(5,7)
       S383331111(3)=-Cij124(5,7)
       S383332111(1)=-Cij124(5,7)
       S383332111(2)=-Cij124(6,7)
       S383332111(3)=-Cij124(6,7)
       S383332211(1)=-Cij124(6,7)
       S383332211(2)=-Cij124(7,7)
       S383332211(3)=-Cij124(7,7)
       S383332221(1)=-Cij124(7,7)
       S383332221(2)=-Cij124(8,7)
       S383332221(3)=-Cij124(8,7)
       S383332222(1)=-Cij124(8,7)
       S383332222(2)=-Cij124(2,7)
       S383332222(3)=-Cij124(2,7)
       S383333111(1)=-Cij124(5,7)
       S383333111(2)=-Cij124(6,7)
       S383333111(3)=-Cij124(6,7)
       S383333211(1)=-Cij124(6,7)
       S383333211(2)=-Cij124(7,7)
       S383333211(3)=-Cij124(7,7)
       S383333221(1)=-Cij124(7,7)
       S383333221(2)=-Cij124(8,7)
       S383333221(3)=-Cij124(8,7)
       S383333222(1)=-Cij124(8,7)
       S383333222(2)=-Cij124(2,7)
       S383333222(3)=-Cij124(2,7)
       S383333311(1)=-Cij124(6,7)
       S383333311(2)=-Cij124(7,7)
       S383333311(3)=-Cij124(7,7)
       S383333321(1)=-Cij124(7,7)
       S383333321(2)=-Cij124(8,7)
       S383333321(3)=-Cij124(8,7)
       S383333322(1)=-Cij124(8,7)
       S383333322(2)=-Cij124(2,7)
       S383333322(3)=-Cij124(2,7)
       S383333331(1)=-Cij124(7,7)
       S383333331(2)=-Cij124(8,7)
       S383333331(3)=-Cij124(8,7)
       S383333332(1)=-Cij124(8,7)
       S383333332(2)=-Cij124(2,7)
       S383333332(3)=-Cij124(2,7)
       S383333333(1)=-Cij124(8,7)
       S383333333(2)=-Cij124(2,7)
       S383333333(3)=-Cij124(2,7)
       aux006111(1,1,1)=-(F(1)*S37111111(1))-F(2)*S37211111(1)-F(3)*S373
     &  11111(1)+S381111111(k)*Z(1,l)+S382111111(k)*Z(2,l)+S383111111(k)
     &  *Z(3,l)+(S300611111(1)-S381111111(1)-S382211111(1)-S383311111(1)
     &  )*Z(k,l)-12*S3h00611111(1)*ZZ(k,1,l,1)-12*S3h00621111(1)*ZZ(k,1,
     &  l,2)-12*S3h00631111(1)*ZZ(k,1,l,3)
       aux006211(1,1,1)=-(F(1)*S37121111(1))-F(2)*S37221111(1)-F(3)*S373
     &  21111(1)+S381211111(k)*Z(1,l)+S382211111(k)*Z(2,l)+S383211111(k)
     &  *Z(3,l)+(S300621111(1)-S381211111(1)-S382221111(1)-S383321111(1)
     &  )*Z(k,l)-10*S3h00612111(1)*ZZ(k,1,l,1)-10*S3h00622111(1)*ZZ(k,1,
     &  l,2)-10*S3h00632111(1)*ZZ(k,1,l,3)-2*S3h00611111(1)*ZZ(k,2,l,1)-
     &  2*S3h00621111(1)*ZZ(k,2,l,2)-2*S3h00631111(1)*ZZ(k,2,l,3)
       aux006221(1,1,1)=-(F(1)*S37122111(1))-F(2)*S37222111(1)-F(3)*S373
     &  22111(1)+S381221111(k)*Z(1,l)+S382221111(k)*Z(2,l)+S383221111(k)
     &  *Z(3,l)+(S300622111(1)-S381221111(1)-S382222111(1)-S383322111(1)
     &  )*Z(k,l)-8*S3h00612211(1)*ZZ(k,1,l,1)-8*S3h00622211(1)*ZZ(k,1,l,
     &  2)-8*S3h00632211(1)*ZZ(k,1,l,3)-4*S3h00612111(1)*ZZ(k,2,l,1)-4*S
     &  3h00622111(1)*ZZ(k,2,l,2)-4*S3h00632111(1)*ZZ(k,2,l,3)
       aux006222(1,1,1)=-(F(1)*S37122211(1))-F(2)*S37222211(1)-F(3)*S373
     &  22211(1)+S381222111(k)*Z(1,l)+S382222111(k)*Z(2,l)+S383222111(k)
     &  *Z(3,l)+(S300622211(1)-S381222111(1)-S382222211(1)-S383322211(1)
     &  )*Z(k,l)-6*S3h00612221(1)*ZZ(k,1,l,1)-6*S3h00622221(1)*ZZ(k,1,l,
     &  2)-6*S3h00632221(1)*ZZ(k,1,l,3)-6*S3h00612211(1)*ZZ(k,2,l,1)-6*S
     &  3h00622211(1)*ZZ(k,2,l,2)-6*S3h00632211(1)*ZZ(k,2,l,3)
       aux006222(2,1,1)=-(F(1)*S37122221(1))-F(2)*S37222221(1)-F(3)*S373
     &  22221(1)+S381222211(k)*Z(1,l)+S382222211(k)*Z(2,l)+S383222211(k)
     &  *Z(3,l)+(S300622221(1)-S381222211(1)-S382222221(1)-S383322221(1)
     &  )*Z(k,l)-4*S3h00612222(1)*ZZ(k,1,l,1)-4*S3h00622222(1)*ZZ(k,1,l,
     &  2)-4*S3h00632222(1)*ZZ(k,1,l,3)-8*S3h00612221(1)*ZZ(k,2,l,1)-8*S
     &  3h00622221(1)*ZZ(k,2,l,2)-8*S3h00632221(1)*ZZ(k,2,l,3)
       aux006222(2,2,1)=-(F(1)*S37122222(1))-F(2)*S37222222(1)-F(3)*S373
     &  22222(1)+S381222221(k)*Z(1,l)+S382222221(k)*Z(2,l)+S383222221(k)
     &  *Z(3,l)+(S300622222(1)-S381222221(1)-S382222222(1)-S383322222(1)
     &  )*Z(k,l)-2*S3h00612222(2)*ZZ(k,1,l,1)-2*S3h00622222(2)*ZZ(k,1,l,
     &  2)-2*S3h00632222(2)*ZZ(k,1,l,3)-10*S3h00612222(1)*ZZ(k,2,l,1)-10
     &  *S3h00622222(1)*ZZ(k,2,l,2)-10*S3h00632222(1)*ZZ(k,2,l,3)
       aux006222(2,2,2)=-(F(1)*S37122222(2))-F(2)*S37222222(2)-F(3)*S373
     &  22222(2)+S381222222(k)*Z(1,l)+S382222222(k)*Z(2,l)+S383222222(k)
     &  *Z(3,l)+(S300622222(2)-S381222222(1)-S382222222(2)-S383322222(2)
     &  )*Z(k,l)-12*S3h00612222(2)*ZZ(k,2,l,1)-12*S3h00622222(2)*ZZ(k,2,
     &  l,2)-12*S3h00632222(2)*ZZ(k,2,l,3)
       aux006311(1,1,1)=-(F(1)*S37131111(1))-F(2)*S37231111(1)-F(3)*S373
     &  31111(1)+S381311111(k)*Z(1,l)+S382311111(k)*Z(2,l)+S383311111(k)
     &  *Z(3,l)+(S300631111(1)-S381311111(1)-S382321111(1)-S383331111(1)
     &  )*Z(k,l)-10*S3h00613111(1)*ZZ(k,1,l,1)-10*S3h00623111(1)*ZZ(k,1,
     &  l,2)-10*S3h00633111(1)*ZZ(k,1,l,3)-2*S3h00611111(1)*ZZ(k,3,l,1)-
     &  2*S3h00621111(1)*ZZ(k,3,l,2)-2*S3h00631111(1)*ZZ(k,3,l,3)
       aux006321(1,1,1)=-(F(1)*S37132111(1))-F(2)*S37232111(1)-F(3)*S373
     &  32111(1)+S381321111(k)*Z(1,l)+S382321111(k)*Z(2,l)+S383321111(k)
     &  *Z(3,l)+(S300632111(1)-S381321111(1)-S382322111(1)-S383332111(1)
     &  )*Z(k,l)-8*S3h00613211(1)*ZZ(k,1,l,1)-8*S3h00623211(1)*ZZ(k,1,l,
     &  2)-8*S3h00633211(1)*ZZ(k,1,l,3)-2*S3h00613111(1)*ZZ(k,2,l,1)-2*S
     &  3h00623111(1)*ZZ(k,2,l,2)-2*S3h00633111(1)*ZZ(k,2,l,3)-2*S3h0061
     &  2111(1)*ZZ(k,3,l,1)-2*S3h00622111(1)*ZZ(k,3,l,2)-2*S3h00632111(1
     &  )*ZZ(k,3,l,3)
       aux006322(1,1,1)=-(F(1)*S37132211(1))-F(2)*S37232211(1)-F(3)*S373
     &  32211(1)+S381322111(k)*Z(1,l)+S382322111(k)*Z(2,l)+S383322111(k)
     &  *Z(3,l)+(S300632211(1)-S381322111(1)-S382322211(1)-S383332211(1)
     &  )*Z(k,l)-6*S3h00613221(1)*ZZ(k,1,l,1)-6*S3h00623221(1)*ZZ(k,1,l,
     &  2)-6*S3h00633221(1)*ZZ(k,1,l,3)-4*S3h00613211(1)*ZZ(k,2,l,1)-4*S
     &  3h00623211(1)*ZZ(k,2,l,2)-4*S3h00633211(1)*ZZ(k,2,l,3)-2*S3h0061
     &  2211(1)*ZZ(k,3,l,1)-2*S3h00622211(1)*ZZ(k,3,l,2)-2*S3h00632211(1
     &  )*ZZ(k,3,l,3)
       aux006322(2,1,1)=-(F(1)*S37132221(1))-F(2)*S37232221(1)-F(3)*S373
     &  32221(1)+S381322211(k)*Z(1,l)+S382322211(k)*Z(2,l)+S383322211(k)
     &  *Z(3,l)+(S300632221(1)-S381322211(1)-S382322221(1)-S383332221(1)
     &  )*Z(k,l)-4*S3h00613222(1)*ZZ(k,1,l,1)-4*S3h00623222(1)*ZZ(k,1,l,
     &  2)-4*S3h00633222(1)*ZZ(k,1,l,3)-6*S3h00613221(1)*ZZ(k,2,l,1)-6*S
     &  3h00623221(1)*ZZ(k,2,l,2)-6*S3h00633221(1)*ZZ(k,2,l,3)-2*S3h0061
     &  2221(1)*ZZ(k,3,l,1)-2*S3h00622221(1)*ZZ(k,3,l,2)-2*S3h00632221(1
     &  )*ZZ(k,3,l,3)
       aux006322(2,2,1)=-(F(1)*S37132222(1))-F(2)*S37232222(1)-F(3)*S373
     &  32222(1)+S381322221(k)*Z(1,l)+S382322221(k)*Z(2,l)+S383322221(k)
     &  *Z(3,l)+(S300632222(1)-S381322221(1)-S382322222(1)-S383332222(1)
     &  )*Z(k,l)-2*S3h00613222(2)*ZZ(k,1,l,1)-2*S3h00623222(2)*ZZ(k,1,l,
     &  2)-2*S3h00633222(2)*ZZ(k,1,l,3)-8*S3h00613222(1)*ZZ(k,2,l,1)-8*S
     &  3h00623222(1)*ZZ(k,2,l,2)-8*S3h00633222(1)*ZZ(k,2,l,3)-2*S3h0061
     &  2222(1)*ZZ(k,3,l,1)-2*S3h00622222(1)*ZZ(k,3,l,2)-2*S3h00632222(1
     &  )*ZZ(k,3,l,3)
       aux006322(2,2,2)=-(F(1)*S37132222(2))-F(2)*S37232222(2)-F(3)*S373
     &  32222(2)+S381322222(k)*Z(1,l)+S382322222(k)*Z(2,l)+S383322222(k)
     &  *Z(3,l)+(S300632222(2)-S381322222(1)-S382322222(2)-S383332222(2)
     &  )*Z(k,l)-10*S3h00613222(2)*ZZ(k,2,l,1)-10*S3h00623222(2)*ZZ(k,2,
     &  l,2)-10*S3h00633222(2)*ZZ(k,2,l,3)-2*S3h00612222(2)*ZZ(k,3,l,1)-
     &  2*S3h00622222(2)*ZZ(k,3,l,2)-2*S3h00632222(2)*ZZ(k,3,l,3)
       aux006331(1,1,1)=-(F(1)*S37133111(1))-F(2)*S37233111(1)-F(3)*S373
     &  33111(1)+S381331111(k)*Z(1,l)+S382331111(k)*Z(2,l)+S383331111(k)
     &  *Z(3,l)+(S300633111(1)-S381331111(1)-S382332111(1)-S383333111(1)
     &  )*Z(k,l)-8*S3h00613311(1)*ZZ(k,1,l,1)-8*S3h00623311(1)*ZZ(k,1,l,
     &  2)-8*S3h00633311(1)*ZZ(k,1,l,3)-4*S3h00613111(1)*ZZ(k,3,l,1)-4*S
     &  3h00623111(1)*ZZ(k,3,l,2)-4*S3h00633111(1)*ZZ(k,3,l,3)
       aux006332(1,1,1)=-(F(1)*S37133211(1))-F(2)*S37233211(1)-F(3)*S373
     &  33211(1)+S381332111(k)*Z(1,l)+S382332111(k)*Z(2,l)+S383332111(k)
     &  *Z(3,l)+(S300633211(1)-S381332111(1)-S382332211(1)-S383333211(1)
     &  )*Z(k,l)-6*S3h00613321(1)*ZZ(k,1,l,1)-6*S3h00623321(1)*ZZ(k,1,l,
     &  2)-6*S3h00633321(1)*ZZ(k,1,l,3)-2*S3h00613311(1)*ZZ(k,2,l,1)-2*S
     &  3h00623311(1)*ZZ(k,2,l,2)-2*S3h00633311(1)*ZZ(k,2,l,3)-4*S3h0061
     &  3211(1)*ZZ(k,3,l,1)-4*S3h00623211(1)*ZZ(k,3,l,2)-4*S3h00633211(1
     &  )*ZZ(k,3,l,3)
       aux006332(2,1,1)=-(F(1)*S37133221(1))-F(2)*S37233221(1)-F(3)*S373
     &  33221(1)+S381332211(k)*Z(1,l)+S382332211(k)*Z(2,l)+S383332211(k)
     &  *Z(3,l)+(S300633221(1)-S381332211(1)-S382332221(1)-S383333221(1)
     &  )*Z(k,l)-4*S3h00613322(1)*ZZ(k,1,l,1)-4*S3h00623322(1)*ZZ(k,1,l,
     &  2)-4*S3h00633322(1)*ZZ(k,1,l,3)-4*S3h00613321(1)*ZZ(k,2,l,1)-4*S
     &  3h00623321(1)*ZZ(k,2,l,2)-4*S3h00633321(1)*ZZ(k,2,l,3)-4*S3h0061
     &  3221(1)*ZZ(k,3,l,1)-4*S3h00623221(1)*ZZ(k,3,l,2)-4*S3h00633221(1
     &  )*ZZ(k,3,l,3)
       aux006332(2,2,1)=-(F(1)*S37133222(1))-F(2)*S37233222(1)-F(3)*S373
     &  33222(1)+S381332221(k)*Z(1,l)+S382332221(k)*Z(2,l)+S383332221(k)
     &  *Z(3,l)+(S300633222(1)-S381332221(1)-S382332222(1)-S383333222(1)
     &  )*Z(k,l)-2*S3h00613322(2)*ZZ(k,1,l,1)-2*S3h00623322(2)*ZZ(k,1,l,
     &  2)-2*S3h00633322(2)*ZZ(k,1,l,3)-6*S3h00613322(1)*ZZ(k,2,l,1)-6*S
     &  3h00623322(1)*ZZ(k,2,l,2)-6*S3h00633322(1)*ZZ(k,2,l,3)-4*S3h0061
     &  3222(1)*ZZ(k,3,l,1)-4*S3h00623222(1)*ZZ(k,3,l,2)-4*S3h00633222(1
     &  )*ZZ(k,3,l,3)
       aux006332(2,2,2)=-(F(1)*S37133222(2))-F(2)*S37233222(2)-F(3)*S373
     &  33222(2)+S381332222(k)*Z(1,l)+S382332222(k)*Z(2,l)+S383332222(k)
     &  *Z(3,l)+(S300633222(2)-S381332222(1)-S382332222(2)-S383333222(2)
     &  )*Z(k,l)-8*S3h00613322(2)*ZZ(k,2,l,1)-8*S3h00623322(2)*ZZ(k,2,l,
     &  2)-8*S3h00633322(2)*ZZ(k,2,l,3)-4*S3h00613222(2)*ZZ(k,3,l,1)-4*S
     &  3h00623222(2)*ZZ(k,3,l,2)-4*S3h00633222(2)*ZZ(k,3,l,3)
       aux006333(1,1,1)=-(F(1)*S37133311(1))-F(2)*S37233311(1)-F(3)*S373
     &  33311(1)+S381333111(k)*Z(1,l)+S382333111(k)*Z(2,l)+S383333111(k)
     &  *Z(3,l)+(S300633311(1)-S381333111(1)-S382333211(1)-S383333311(1)
     &  )*Z(k,l)-6*S3h00613331(1)*ZZ(k,1,l,1)-6*S3h00623331(1)*ZZ(k,1,l,
     &  2)-6*S3h00633331(1)*ZZ(k,1,l,3)-6*S3h00613311(1)*ZZ(k,3,l,1)-6*S
     &  3h00623311(1)*ZZ(k,3,l,2)-6*S3h00633311(1)*ZZ(k,3,l,3)
       aux006333(2,1,1)=-(F(1)*S37133321(1))-F(2)*S37233321(1)-F(3)*S373
     &  33321(1)+S381333211(k)*Z(1,l)+S382333211(k)*Z(2,l)+S383333211(k)
     &  *Z(3,l)+(S300633321(1)-S381333211(1)-S382333221(1)-S383333321(1)
     &  )*Z(k,l)-4*S3h00613332(1)*ZZ(k,1,l,1)-4*S3h00623332(1)*ZZ(k,1,l,
     &  2)-4*S3h00633332(1)*ZZ(k,1,l,3)-2*S3h00613331(1)*ZZ(k,2,l,1)-2*S
     &  3h00623331(1)*ZZ(k,2,l,2)-2*S3h00633331(1)*ZZ(k,2,l,3)-6*S3h0061
     &  3321(1)*ZZ(k,3,l,1)-6*S3h00623321(1)*ZZ(k,3,l,2)-6*S3h00633321(1
     &  )*ZZ(k,3,l,3)
       aux006333(2,2,1)=-(F(1)*S37133322(1))-F(2)*S37233322(1)-F(3)*S373
     &  33322(1)+S381333221(k)*Z(1,l)+S382333221(k)*Z(2,l)+S383333221(k)
     &  *Z(3,l)+(S300633322(1)-S381333221(1)-S382333222(1)-S383333322(1)
     &  )*Z(k,l)-2*S3h00613332(2)*ZZ(k,1,l,1)-2*S3h00623332(2)*ZZ(k,1,l,
     &  2)-2*S3h00633332(2)*ZZ(k,1,l,3)-4*S3h00613332(1)*ZZ(k,2,l,1)-4*S
     &  3h00623332(1)*ZZ(k,2,l,2)-4*S3h00633332(1)*ZZ(k,2,l,3)-6*S3h0061
     &  3322(1)*ZZ(k,3,l,1)-6*S3h00623322(1)*ZZ(k,3,l,2)-6*S3h00633322(1
     &  )*ZZ(k,3,l,3)
       aux006333(2,2,2)=-(F(1)*S37133322(2))-F(2)*S37233322(2)-F(3)*S373
     &  33322(2)+S381333222(k)*Z(1,l)+S382333222(k)*Z(2,l)+S383333222(k)
     &  *Z(3,l)+(S300633322(2)-S381333222(1)-S382333222(2)-S383333322(2)
     &  )*Z(k,l)-6*S3h00613332(2)*ZZ(k,2,l,1)-6*S3h00623332(2)*ZZ(k,2,l,
     &  2)-6*S3h00633332(2)*ZZ(k,2,l,3)-6*S3h00613322(2)*ZZ(k,3,l,1)-6*S
     &  3h00623322(2)*ZZ(k,3,l,2)-6*S3h00633322(2)*ZZ(k,3,l,3)
       aux006333(3,1,1)=-(F(1)*S37133331(1))-F(2)*S37233331(1)-F(3)*S373
     &  33331(1)+S381333311(k)*Z(1,l)+S382333311(k)*Z(2,l)+S383333311(k)
     &  *Z(3,l)+(S300633331(1)-S381333311(1)-S382333321(1)-S383333331(1)
     &  )*Z(k,l)-4*S3h00613333(1)*ZZ(k,1,l,1)-4*S3h00623333(1)*ZZ(k,1,l,
     &  2)-4*S3h00633333(1)*ZZ(k,1,l,3)-8*S3h00613331(1)*ZZ(k,3,l,1)-8*S
     &  3h00623331(1)*ZZ(k,3,l,2)-8*S3h00633331(1)*ZZ(k,3,l,3)
       aux006333(3,2,1)=-(F(1)*S37133332(1))-F(2)*S37233332(1)-F(3)*S373
     &  33332(1)+S381333321(k)*Z(1,l)+S382333321(k)*Z(2,l)+S383333321(k)
     &  *Z(3,l)+(S300633332(1)-S381333321(1)-S382333322(1)-S383333332(1)
     &  )*Z(k,l)-2*S3h00613333(2)*ZZ(k,1,l,1)-2*S3h00623333(2)*ZZ(k,1,l,
     &  2)-2*S3h00633333(2)*ZZ(k,1,l,3)-2*S3h00613333(1)*ZZ(k,2,l,1)-2*S
     &  3h00623333(1)*ZZ(k,2,l,2)-2*S3h00633333(1)*ZZ(k,2,l,3)-8*S3h0061
     &  3332(1)*ZZ(k,3,l,1)-8*S3h00623332(1)*ZZ(k,3,l,2)-8*S3h00633332(1
     &  )*ZZ(k,3,l,3)
       aux006333(3,2,2)=-(F(1)*S37133332(2))-F(2)*S37233332(2)-F(3)*S373
     &  33332(2)+S381333322(k)*Z(1,l)+S382333322(k)*Z(2,l)+S383333322(k)
     &  *Z(3,l)+(S300633332(2)-S381333322(1)-S382333322(2)-S383333332(2)
     &  )*Z(k,l)-4*S3h00613333(2)*ZZ(k,2,l,1)-4*S3h00623333(2)*ZZ(k,2,l,
     &  2)-4*S3h00633333(2)*ZZ(k,2,l,3)-8*S3h00613332(2)*ZZ(k,3,l,1)-8*S
     &  3h00623332(2)*ZZ(k,3,l,2)-8*S3h00633332(2)*ZZ(k,3,l,3)
       aux006333(3,3,1)=-(F(1)*S37133333(1))-F(2)*S37233333(1)-F(3)*S373
     &  33333(1)+S381333331(k)*Z(1,l)+S382333331(k)*Z(2,l)+S383333331(k)
     &  *Z(3,l)+(S300633333(1)-S381333331(1)-S382333332(1)-S383333333(1)
     &  )*Z(k,l)-2*S3h00613333(3)*ZZ(k,1,l,1)-2*S3h00623333(3)*ZZ(k,1,l,
     &  2)-2*S3h00633333(3)*ZZ(k,1,l,3)-10*S3h00613333(1)*ZZ(k,3,l,1)-10
     &  *S3h00623333(1)*ZZ(k,3,l,2)-10*S3h00633333(1)*ZZ(k,3,l,3)
       aux006333(3,3,2)=-(F(1)*S37133333(2))-F(2)*S37233333(2)-F(3)*S373
     &  33333(2)+S381333332(k)*Z(1,l)+S382333332(k)*Z(2,l)+S383333332(k)
     &  *Z(3,l)+(S300633333(2)-S381333332(1)-S382333332(2)-S383333333(2)
     &  )*Z(k,l)-2*S3h00613333(3)*ZZ(k,2,l,1)-2*S3h00623333(3)*ZZ(k,2,l,
     &  2)-2*S3h00633333(3)*ZZ(k,2,l,3)-10*S3h00613333(2)*ZZ(k,3,l,1)-10
     &  *S3h00623333(2)*ZZ(k,3,l,2)-10*S3h00633333(2)*ZZ(k,3,l,3)
       aux006333(3,3,3)=-(F(1)*S37133333(3))-F(2)*S37233333(3)-F(3)*S373
     &  33333(3)+S381333333(k)*Z(1,l)+S382333333(k)*Z(2,l)+S383333333(k)
     &  *Z(3,l)+(S300633333(3)-S381333333(1)-S382333333(2)-S383333333(3)
     &  )*Z(k,l)-12*S3h00613333(3)*ZZ(k,3,l,1)-12*S3h00623333(3)*ZZ(k,3,
     &  l,2)-12*S3h00633333(3)*ZZ(k,3,l,3)
       temp006111(1,1,1)=I28Z*(aux006111(1,1,1)+12*F(4)*temp00511(1,1,1)
     &  +F(5)*temp6111(1,1,1)+120*temp000041(1,1,1)*ZZ(k,1,l,1))
       temp006211(1,1,1)=I28Z*(aux006211(1,1,1)+2*F(6)*temp00511(1,1,1)+
     &  10*F(4)*temp00521(1,1,1)+F(5)*temp6211(1,1,1)+80*temp000042(1,1,
     &  1)*ZZ(k,1,l,1)+20*temp000041(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp006221(1,1,1)=I28Z*(aux006221(1,1,1)+4*F(6)*temp00521(1,1,1)+
     &  F(5)*temp6221(1,1,1)+48*temp000042(2,1,1)*ZZ(k,1,l,1)+32*temp000
     &  042(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(F(4)*temp00522(1,1,1)+te
     &  mp000041(1,1,1)*ZZ(k,2,l,2)))
       temp006222(1,1,1)=I28Z*(aux006222(1,1,1)+6*F(6)*temp00522(1,1,1)+
     &  6*F(4)*temp00522(2,1,1)+F(5)*temp6222(1,1,1)+36*temp000042(2,1,1
     &  )*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*(temp000042(2,2,1)*ZZ(k,1,l,1)+te
     &  mp000042(1,1,1)*ZZ(k,2,l,2)))
       temp006222(2,1,1)=I28Z*(aux006222(2,1,1)+4*F(4)*temp00522(2,2,1)+
     &  F(5)*temp6222(2,1,1)+8*(F(6)*temp00522(2,1,1)+temp000042(2,2,2)*
     &  ZZ(k,1,l,1))+32*temp000042(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*t
     &  emp000042(2,1,1)*ZZ(k,2,l,2))
       temp006222(2,2,1)=I28Z*(aux006222(2,2,1)+10*F(6)*temp00522(2,2,1)
     &  +2*F(4)*temp00522(2,2,2)+F(5)*temp6222(2,2,1)+20*temp000042(2,2,
     &  2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+80*temp000042(2,2,1)*ZZ(k,2,l,2))
       temp006222(2,2,2)=I28Z*(aux006222(2,2,2)+12*F(6)*temp00522(2,2,2)
     &  +F(5)*temp6222(2,2,2)+120*temp000042(2,2,2)*ZZ(k,2,l,2))
       temp006311(1,1,1)=I28Z*(aux006311(1,1,1)+2*F(7)*temp00511(1,1,1)+
     &  10*F(4)*temp00531(1,1,1)+F(5)*temp6311(1,1,1)+80*temp000043(1,1,
     &  1)*ZZ(k,1,l,1)+20*temp000041(1,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)))
       temp006321(1,1,1)=I28Z*(aux006321(1,1,1)+2*F(7)*temp00521(1,1,1)+
     &  2*F(6)*temp00531(1,1,1)+8*F(4)*temp00532(1,1,1)+F(5)*temp6321(1,
     &  1,1)+48*temp000043(2,1,1)*ZZ(k,1,l,1)+16*(temp000043(1,1,1)*(ZZ(
     &  k,1,l,2)+ZZ(k,2,l,1))+temp000042(1,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)
     &  ))+4*temp000041(1,1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp006322(1,1,1)=I28Z*(aux006322(1,1,1)+2*F(7)*temp00522(1,1,1)+
     &  4*F(6)*temp00532(1,1,1)+6*F(4)*temp00532(2,1,1)+F(5)*temp6322(1,
     &  1,1)+24*(temp000043(2,2,1)*ZZ(k,1,l,1)+temp000043(2,1,1)*(ZZ(k,1
     &  ,l,2)+ZZ(k,2,l,1)))+8*temp000043(1,1,1)*ZZ(k,2,l,2)+12*temp00004
     &  2(2,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*temp000042(1,1,1)*(ZZ(k,2,l
     &  ,3)+ZZ(k,3,l,2)))
       temp006322(2,1,1)=I28Z*(aux006322(2,1,1)+2*F(7)*temp00522(2,1,1)+
     &  6*F(6)*temp00532(2,1,1)+4*F(4)*temp00532(2,2,1)+F(5)*temp6322(2,
     &  1,1)+24*temp000043(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*temp00004
     &  3(2,1,1)*ZZ(k,2,l,2)+8*(temp000043(2,2,2)*ZZ(k,1,l,1)+temp000042
     &  (2,2,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)))+12*temp000042(2,1,1)*(ZZ(k,2,
     &  l,3)+ZZ(k,3,l,2)))
       temp006322(2,2,1)=I28Z*(aux006322(2,2,1)+2*F(7)*temp00522(2,2,1)+
     &  8*F(6)*temp00532(2,2,1)+2*F(4)*temp00532(2,2,2)+F(5)*temp6322(2,
     &  2,1)+48*temp000043(2,2,1)*ZZ(k,2,l,2)+4*temp000042(2,2,2)*(ZZ(k,
     &  1,l,3)+ZZ(k,3,l,1))+16*(temp000043(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,
     &  1))+temp000042(2,2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))))
       temp006322(2,2,2)=I28Z*(aux006322(2,2,2)+2*F(7)*temp00522(2,2,2)+
     &  10*F(6)*temp00532(2,2,2)+F(5)*temp6322(2,2,2)+80*temp000043(2,2,
     &  2)*ZZ(k,2,l,2)+20*temp000042(2,2,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp006331(1,1,1)=I28Z*(aux006331(1,1,1)+4*F(7)*temp00531(1,1,1)+
     &  F(5)*temp6331(1,1,1)+48*temp000043(3,1,1)*ZZ(k,1,l,1)+32*temp000
     &  043(1,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*(F(4)*temp00533(1,1,1)+te
     &  mp000041(1,1,1)*ZZ(k,3,l,3)))
       temp006332(1,1,1)=I28Z*(aux006332(1,1,1)+4*F(7)*temp00532(1,1,1)+
     &  2*F(6)*temp00533(1,1,1)+6*F(4)*temp00533(2,1,1)+F(5)*temp6332(1,
     &  1,1)+12*temp000043(3,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*(temp0000
     &  43(3,2,1)*ZZ(k,1,l,1)+temp000043(2,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)
     &  ))+8*temp000043(1,1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*temp000042(1,
     &  1,1)*ZZ(k,3,l,3))
       temp006332(2,1,1)=I28Z*(aux006332(2,1,1)+4*F(7)*temp00532(2,1,1)+
     &  4*F(6)*temp00533(2,1,1)+4*F(4)*temp00533(2,2,1)+F(5)*temp6332(2,
     &  1,1)+16*temp000043(3,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp00004
     &  3(3,2,2)*ZZ(k,1,l,1)+temp000043(3,1,1)*ZZ(k,2,l,2))+16*temp00004
     &  3(2,2,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+16*temp000043(2,1,1)*(ZZ(k,2,
     &  l,3)+ZZ(k,3,l,2))+8*temp000042(2,1,1)*ZZ(k,3,l,3))
       temp006332(2,2,1)=I28Z*(aux006332(2,2,1)+4*F(7)*temp00532(2,2,1)+
     &  6*F(6)*temp00533(2,2,1)+2*F(4)*temp00533(2,2,2)+F(5)*temp6332(2,
     &  2,1)+12*temp000043(3,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*temp000043
     &  (2,2,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+24*(temp000043(3,2,1)*ZZ(k,2,l
     &  ,2)+temp000043(2,2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))+8*temp000042(2,
     &  2,1)*ZZ(k,3,l,3))
       temp006332(2,2,2)=I28Z*(aux006332(2,2,2)+4*F(7)*temp00532(2,2,2)+
     &  F(5)*temp6332(2,2,2)+48*temp000043(3,2,2)*ZZ(k,2,l,2)+32*temp000
     &  043(2,2,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*(F(6)*temp00533(2,2,2)+te
     &  mp000042(2,2,2)*ZZ(k,3,l,3)))
       temp006333(1,1,1)=I28Z*(aux006333(1,1,1)+6*F(7)*temp00533(1,1,1)+
     &  6*F(4)*temp00533(3,1,1)+F(5)*temp6333(1,1,1)+36*temp000043(3,1,1
     &  )*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+24*(temp000043(3,3,1)*ZZ(k,1,l,1)+te
     &  mp000043(1,1,1)*ZZ(k,3,l,3)))
       temp006333(2,1,1)=I28Z*(aux006333(2,1,1)+6*F(7)*temp00533(2,1,1)+
     &  2*F(6)*temp00533(3,1,1)+4*F(4)*temp00533(3,2,1)+F(5)*temp6333(2,
     &  1,1)+8*(temp000043(3,3,2)*ZZ(k,1,l,1)+temp000043(3,3,1)*(ZZ(k,1,
     &  l,2)+ZZ(k,2,l,1)))+24*temp000043(3,2,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)
     &  )+12*temp000043(3,1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+24*temp000043(2
     &  ,1,1)*ZZ(k,3,l,3))
       temp006333(2,2,1)=I28Z*(aux006333(2,2,1)+6*F(7)*temp00533(2,2,1)+
     &  4*F(6)*temp00533(3,2,1)+2*F(4)*temp00533(3,2,2)+F(5)*temp6333(2,
     &  2,1)+8*(temp000043(3,3,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp000043(3
     &  ,3,1)*ZZ(k,2,l,2))+12*temp000043(3,2,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)
     &  )+24*temp000043(3,2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+24*temp000043(2
     &  ,2,1)*ZZ(k,3,l,3))
       temp006333(2,2,2)=I28Z*(aux006333(2,2,2)+6*F(7)*temp00533(2,2,2)+
     &  6*F(6)*temp00533(3,2,2)+F(5)*temp6333(2,2,2)+36*temp000043(3,2,2
     &  )*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+24*(temp000043(3,3,2)*ZZ(k,2,l,2)+te
     &  mp000043(2,2,2)*ZZ(k,3,l,3)))
       temp006333(3,1,1)=I28Z*(aux006333(3,1,1)+4*F(4)*temp00533(3,3,1)+
     &  F(5)*temp6333(3,1,1)+8*(F(7)*temp00533(3,1,1)+temp000043(3,3,3)*
     &  ZZ(k,1,l,1))+32*temp000043(3,3,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+48*t
     &  emp000043(3,1,1)*ZZ(k,3,l,3))
       temp006333(3,2,1)=I28Z*(aux006333(3,2,1)+8*F(7)*temp00533(3,2,1)+
     &  2*F(6)*temp00533(3,3,1)+2*F(4)*temp00533(3,3,2)+F(5)*temp6333(3,
     &  2,1)+4*temp000043(3,3,3)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+16*(temp00004
     &  3(3,3,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+temp000043(3,3,1)*(ZZ(k,2,l,3
     &  )+ZZ(k,3,l,2)))+48*temp000043(3,2,1)*ZZ(k,3,l,3))
       temp006333(3,2,2)=I28Z*(aux006333(3,2,2)+4*F(6)*temp00533(3,3,2)+
     &  F(5)*temp6333(3,2,2)+8*(F(7)*temp00533(3,2,2)+temp000043(3,3,3)*
     &  ZZ(k,2,l,2))+32*temp000043(3,3,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+48*t
     &  emp000043(3,2,2)*ZZ(k,3,l,3))
       temp006333(3,3,1)=I28Z*(aux006333(3,3,1)+10*F(7)*temp00533(3,3,1)
     &  +2*F(4)*temp00533(3,3,3)+F(5)*temp6333(3,3,1)+20*temp000043(3,3,
     &  3)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+80*temp000043(3,3,1)*ZZ(k,3,l,3))
       temp006333(3,3,2)=I28Z*(aux006333(3,3,2)+10*F(7)*temp00533(3,3,2)
     &  +2*F(6)*temp00533(3,3,3)+F(5)*temp6333(3,3,2)+20*temp000043(3,3,
     &  3)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+80*temp000043(3,3,2)*ZZ(k,3,l,3))
       temp006333(3,3,3)=I28Z*(aux006333(3,3,3)+12*F(7)*temp00533(3,3,3)
     &  +F(5)*temp6333(3,3,3)+120*temp000043(3,3,3)*ZZ(k,3,l,3))
       temp006111(1,1,2)=temp006211(1,1,1)
       temp006111(1,1,3)=temp006311(1,1,1)
       temp006111(1,2,1)=temp006211(1,1,1)
       temp006111(1,2,2)=temp006221(1,1,1)
       temp006111(1,2,3)=temp006321(1,1,1)
       temp006111(1,3,1)=temp006311(1,1,1)
       temp006111(1,3,2)=temp006321(1,1,1)
       temp006111(1,3,3)=temp006331(1,1,1)
       temp006211(1,1,2)=temp006221(1,1,1)
       temp006211(1,1,3)=temp006321(1,1,1)
       temp006211(1,2,1)=temp006221(1,1,1)
       temp006211(1,2,2)=temp006222(1,1,1)
       temp006211(1,2,3)=temp006322(1,1,1)
       temp006211(1,3,1)=temp006321(1,1,1)
       temp006211(1,3,2)=temp006322(1,1,1)
       temp006211(1,3,3)=temp006332(1,1,1)
       temp006221(1,1,2)=temp006222(1,1,1)
       temp006221(1,1,3)=temp006322(1,1,1)
       temp006221(1,2,1)=temp006222(1,1,1)
       temp006221(1,2,2)=temp006222(2,1,1)
       temp006221(1,2,3)=temp006322(2,1,1)
       temp006221(1,3,1)=temp006322(1,1,1)
       temp006221(1,3,2)=temp006322(2,1,1)
       temp006221(1,3,3)=temp006332(2,1,1)
       temp006222(1,1,2)=temp006222(2,1,1)
       temp006222(1,1,3)=temp006322(2,1,1)
       temp006222(1,2,1)=temp006222(2,1,1)
       temp006222(1,2,2)=temp006222(2,2,1)
       temp006222(1,2,3)=temp006322(2,2,1)
       temp006222(1,3,1)=temp006322(2,1,1)
       temp006222(1,3,2)=temp006322(2,2,1)
       temp006222(1,3,3)=temp006332(2,2,1)
       temp006222(2,1,2)=temp006222(2,2,1)
       temp006222(2,1,3)=temp006322(2,2,1)
       temp006222(2,2,3)=temp006322(2,2,2)
       temp006222(2,3,1)=temp006322(2,2,1)
       temp006222(2,3,2)=temp006322(2,2,2)
       temp006222(2,3,3)=temp006332(2,2,2)
       temp006311(1,1,2)=temp006321(1,1,1)
       temp006311(1,1,3)=temp006331(1,1,1)
       temp006311(1,2,1)=temp006321(1,1,1)
       temp006311(1,2,2)=temp006322(1,1,1)
       temp006311(1,2,3)=temp006332(1,1,1)
       temp006311(1,3,1)=temp006331(1,1,1)
       temp006311(1,3,2)=temp006332(1,1,1)
       temp006311(1,3,3)=temp006333(1,1,1)
       temp006321(1,1,2)=temp006322(1,1,1)
       temp006321(1,1,3)=temp006332(1,1,1)
       temp006321(1,2,1)=temp006322(1,1,1)
       temp006321(1,2,2)=temp006322(2,1,1)
       temp006321(1,2,3)=temp006332(2,1,1)
       temp006321(1,3,1)=temp006332(1,1,1)
       temp006321(1,3,2)=temp006332(2,1,1)
       temp006321(1,3,3)=temp006333(2,1,1)
       temp006322(1,1,2)=temp006322(2,1,1)
       temp006322(1,1,3)=temp006332(2,1,1)
       temp006322(1,2,1)=temp006322(2,1,1)
       temp006322(1,2,2)=temp006322(2,2,1)
       temp006322(1,2,3)=temp006332(2,2,1)
       temp006322(1,3,1)=temp006332(2,1,1)
       temp006322(1,3,2)=temp006332(2,2,1)
       temp006322(1,3,3)=temp006333(2,2,1)
       temp006322(2,1,2)=temp006322(2,2,1)
       temp006322(2,1,3)=temp006332(2,2,1)
       temp006322(2,2,3)=temp006332(2,2,2)
       temp006322(2,3,1)=temp006332(2,2,1)
       temp006322(2,3,2)=temp006332(2,2,2)
       temp006322(2,3,3)=temp006333(2,2,2)
       temp006331(1,1,2)=temp006332(1,1,1)
       temp006331(1,1,3)=temp006333(1,1,1)
       temp006331(1,2,1)=temp006332(1,1,1)
       temp006331(1,2,2)=temp006332(2,1,1)
       temp006331(1,2,3)=temp006333(2,1,1)
       temp006331(1,3,1)=temp006333(1,1,1)
       temp006331(1,3,2)=temp006333(2,1,1)
       temp006331(1,3,3)=temp006333(3,1,1)
       temp006332(1,1,2)=temp006332(2,1,1)
       temp006332(1,1,3)=temp006333(2,1,1)
       temp006332(1,2,1)=temp006332(2,1,1)
       temp006332(1,2,2)=temp006332(2,2,1)
       temp006332(1,2,3)=temp006333(2,2,1)
       temp006332(1,3,1)=temp006333(2,1,1)
       temp006332(1,3,2)=temp006333(2,2,1)
       temp006332(1,3,3)=temp006333(3,2,1)
       temp006332(2,1,2)=temp006332(2,2,1)
       temp006332(2,1,3)=temp006333(2,2,1)
       temp006332(2,2,3)=temp006333(2,2,2)
       temp006332(2,3,1)=temp006333(2,2,1)
       temp006332(2,3,2)=temp006333(2,2,2)
       temp006332(2,3,3)=temp006333(3,2,2)
       temp006333(1,1,2)=temp006333(2,1,1)
       temp006333(1,1,3)=temp006333(3,1,1)
       temp006333(1,2,1)=temp006333(2,1,1)
       temp006333(1,2,2)=temp006333(2,2,1)
       temp006333(1,2,3)=temp006333(3,2,1)
       temp006333(1,3,1)=temp006333(3,1,1)
       temp006333(1,3,2)=temp006333(3,2,1)
       temp006333(1,3,3)=temp006333(3,3,1)
       temp006333(2,1,2)=temp006333(2,2,1)
       temp006333(2,1,3)=temp006333(3,2,1)
       temp006333(2,2,3)=temp006333(3,2,2)
       temp006333(2,3,1)=temp006333(3,2,1)
       temp006333(2,3,2)=temp006333(3,2,2)
       temp006333(2,3,3)=temp006333(3,3,2)
       temp006333(3,1,2)=temp006333(3,2,1)
       temp006333(3,1,3)=temp006333(3,3,1)
       temp006333(3,2,3)=temp006333(3,3,2)
       aux71111(1,1,1)=-(S381111111(1)*Z(jj,1))-S382111111(1)*Z(jj,2)-S3
     &  83111111(1)*Z(jj,3)
       aux72111(1,1,1)=-(S381211111(1)*Z(jj,1))-S382211111(1)*Z(jj,2)-S3
     &  83211111(1)*Z(jj,3)
       aux72211(1,1,1)=-(S381221111(1)*Z(jj,1))-S382221111(1)*Z(jj,2)-S3
     &  83221111(1)*Z(jj,3)
       aux72221(1,1,1)=-(S381222111(1)*Z(jj,1))-S382222111(1)*Z(jj,2)-S3
     &  83222111(1)*Z(jj,3)
       aux72222(1,1,1)=-(S381222211(1)*Z(jj,1))-S382222211(1)*Z(jj,2)-S3
     &  83222211(1)*Z(jj,3)
       aux72222(2,1,1)=-(S381222221(1)*Z(jj,1))-S382222221(1)*Z(jj,2)-S3
     &  83222221(1)*Z(jj,3)
       aux72222(2,2,1)=-(S381222222(1)*Z(jj,1))-S382222222(1)*Z(jj,2)-S3
     &  83222222(1)*Z(jj,3)
       aux72222(2,2,2)=-(S381222222(2)*Z(jj,1))-S382222222(2)*Z(jj,2)-S3
     &  83222222(2)*Z(jj,3)
       aux73111(1,1,1)=-(S381311111(1)*Z(jj,1))-S382311111(1)*Z(jj,2)-S3
     &  83311111(1)*Z(jj,3)
       aux73211(1,1,1)=-(S381321111(1)*Z(jj,1))-S382321111(1)*Z(jj,2)-S3
     &  83321111(1)*Z(jj,3)
       aux73221(1,1,1)=-(S381322111(1)*Z(jj,1))-S382322111(1)*Z(jj,2)-S3
     &  83322111(1)*Z(jj,3)
       aux73222(1,1,1)=-(S381322211(1)*Z(jj,1))-S382322211(1)*Z(jj,2)-S3
     &  83322211(1)*Z(jj,3)
       aux73222(2,1,1)=-(S381322221(1)*Z(jj,1))-S382322221(1)*Z(jj,2)-S3
     &  83322221(1)*Z(jj,3)
       aux73222(2,2,1)=-(S381322222(1)*Z(jj,1))-S382322222(1)*Z(jj,2)-S3
     &  83322222(1)*Z(jj,3)
       aux73222(2,2,2)=-(S381322222(2)*Z(jj,1))-S382322222(2)*Z(jj,2)-S3
     &  83322222(2)*Z(jj,3)
       aux73311(1,1,1)=-(S381331111(1)*Z(jj,1))-S382331111(1)*Z(jj,2)-S3
     &  83331111(1)*Z(jj,3)
       aux73321(1,1,1)=-(S381332111(1)*Z(jj,1))-S382332111(1)*Z(jj,2)-S3
     &  83332111(1)*Z(jj,3)
       aux73322(1,1,1)=-(S381332211(1)*Z(jj,1))-S382332211(1)*Z(jj,2)-S3
     &  83332211(1)*Z(jj,3)
       aux73322(2,1,1)=-(S381332221(1)*Z(jj,1))-S382332221(1)*Z(jj,2)-S3
     &  83332221(1)*Z(jj,3)
       aux73322(2,2,1)=-(S381332222(1)*Z(jj,1))-S382332222(1)*Z(jj,2)-S3
     &  83332222(1)*Z(jj,3)
       aux73322(2,2,2)=-(S381332222(2)*Z(jj,1))-S382332222(2)*Z(jj,2)-S3
     &  83332222(2)*Z(jj,3)
       aux73331(1,1,1)=-(S381333111(1)*Z(jj,1))-S382333111(1)*Z(jj,2)-S3
     &  83333111(1)*Z(jj,3)
       aux73332(1,1,1)=-(S381333211(1)*Z(jj,1))-S382333211(1)*Z(jj,2)-S3
     &  83333211(1)*Z(jj,3)
       aux73332(2,1,1)=-(S381333221(1)*Z(jj,1))-S382333221(1)*Z(jj,2)-S3
     &  83333221(1)*Z(jj,3)
       aux73332(2,2,1)=-(S381333222(1)*Z(jj,1))-S382333222(1)*Z(jj,2)-S3
     &  83333222(1)*Z(jj,3)
       aux73332(2,2,2)=-(S381333222(2)*Z(jj,1))-S382333222(2)*Z(jj,2)-S3
     &  83333222(2)*Z(jj,3)
       aux73333(1,1,1)=-(S381333311(1)*Z(jj,1))-S382333311(1)*Z(jj,2)-S3
     &  83333311(1)*Z(jj,3)
       aux73333(2,1,1)=-(S381333321(1)*Z(jj,1))-S382333321(1)*Z(jj,2)-S3
     &  83333321(1)*Z(jj,3)
       aux73333(2,2,1)=-(S381333322(1)*Z(jj,1))-S382333322(1)*Z(jj,2)-S3
     &  83333322(1)*Z(jj,3)
       aux73333(2,2,2)=-(S381333322(2)*Z(jj,1))-S382333322(2)*Z(jj,2)-S3
     &  83333322(2)*Z(jj,3)
       aux73333(3,1,1)=-(S381333331(1)*Z(jj,1))-S382333331(1)*Z(jj,2)-S3
     &  83333331(1)*Z(jj,3)
       aux73333(3,2,1)=-(S381333332(1)*Z(jj,1))-S382333332(1)*Z(jj,2)-S3
     &  83333332(1)*Z(jj,3)
       aux73333(3,2,2)=-(S381333332(2)*Z(jj,1))-S382333332(2)*Z(jj,2)-S3
     &  83333332(2)*Z(jj,3)
       aux73333(3,3,1)=-(S381333333(1)*Z(jj,1))-S382333333(1)*Z(jj,2)-S3
     &  83333333(1)*Z(jj,3)
       aux73333(3,3,2)=-(S381333333(2)*Z(jj,1))-S382333333(2)*Z(jj,2)-S3
     &  83333333(2)*Z(jj,3)
       aux73333(3,3,3)=-(S381333333(3)*Z(jj,1))-S382333333(3)*Z(jj,2)-S3
     &  83333333(3)*Z(jj,3)
       temp71111(1,1,1)=IX*(aux71111(1,1,1)+14*temp006111(1,1,1)*Z(jj,1)
     &  )
       temp72111(1,1,1)=IX*(aux72111(1,1,1)+12*temp006211(1,1,1)*Z(jj,1)
     &  +2*temp006111(1,1,1)*Z(jj,2))
       temp72211(1,1,1)=IX*(aux72211(1,1,1)+10*temp006221(1,1,1)*Z(jj,1)
     &  +4*temp006211(1,1,1)*Z(jj,2))
       temp72221(1,1,1)=IX*(aux72221(1,1,1)+8*temp006222(1,1,1)*Z(jj,1)+
     &  6*temp006221(1,1,1)*Z(jj,2))
       temp72222(1,1,1)=IX*(aux72222(1,1,1)+6*temp006222(2,1,1)*Z(jj,1)+
     &  8*temp006222(1,1,1)*Z(jj,2))
       temp72222(2,1,1)=IX*(aux72222(2,1,1)+4*temp006222(2,2,1)*Z(jj,1)+
     &  10*temp006222(2,1,1)*Z(jj,2))
       temp72222(2,2,1)=IX*(aux72222(2,2,1)+2*temp006222(2,2,2)*Z(jj,1)+
     &  12*temp006222(2,2,1)*Z(jj,2))
       temp72222(2,2,2)=IX*(aux72222(2,2,2)+14*temp006222(2,2,2)*Z(jj,2)
     &  )
       temp73111(1,1,1)=IX*(aux73111(1,1,1)+12*temp006311(1,1,1)*Z(jj,1)
     &  +2*temp006111(1,1,1)*Z(jj,3))
       temp73211(1,1,1)=IX*(aux73211(1,1,1)+10*temp006321(1,1,1)*Z(jj,1)
     &  +2*(temp006311(1,1,1)*Z(jj,2)+temp006211(1,1,1)*Z(jj,3)))
       temp73221(1,1,1)=IX*(aux73221(1,1,1)+8*temp006322(1,1,1)*Z(jj,1)+
     &  4*temp006321(1,1,1)*Z(jj,2)+2*temp006221(1,1,1)*Z(jj,3))
       temp73222(1,1,1)=IX*(aux73222(1,1,1)+6*(temp006322(2,1,1)*Z(jj,1)
     &  +temp006322(1,1,1)*Z(jj,2))+2*temp006222(1,1,1)*Z(jj,3))
       temp73222(2,1,1)=IX*(aux73222(2,1,1)+4*temp006322(2,2,1)*Z(jj,1)+
     &  8*temp006322(2,1,1)*Z(jj,2)+2*temp006222(2,1,1)*Z(jj,3))
       temp73222(2,2,1)=IX*(aux73222(2,2,1)+10*temp006322(2,2,1)*Z(jj,2)
     &  +2*(temp006322(2,2,2)*Z(jj,1)+temp006222(2,2,1)*Z(jj,3)))
       temp73222(2,2,2)=IX*(aux73222(2,2,2)+12*temp006322(2,2,2)*Z(jj,2)
     &  +2*temp006222(2,2,2)*Z(jj,3))
       temp73311(1,1,1)=IX*(aux73311(1,1,1)+10*temp006331(1,1,1)*Z(jj,1)
     &  +4*temp006311(1,1,1)*Z(jj,3))
       temp73321(1,1,1)=IX*(aux73321(1,1,1)+8*temp006332(1,1,1)*Z(jj,1)+
     &  2*temp006331(1,1,1)*Z(jj,2)+4*temp006321(1,1,1)*Z(jj,3))
       temp73322(1,1,1)=IX*(aux73322(1,1,1)+6*temp006332(2,1,1)*Z(jj,1)+
     &  4*(temp006332(1,1,1)*Z(jj,2)+temp006322(1,1,1)*Z(jj,3)))
       temp73322(2,1,1)=IX*(aux73322(2,1,1)+6*temp006332(2,1,1)*Z(jj,2)+
     &  4*(temp006332(2,2,1)*Z(jj,1)+temp006322(2,1,1)*Z(jj,3)))
       temp73322(2,2,1)=IX*(aux73322(2,2,1)+2*temp006332(2,2,2)*Z(jj,1)+
     &  8*temp006332(2,2,1)*Z(jj,2)+4*temp006322(2,2,1)*Z(jj,3))
       temp73322(2,2,2)=IX*(aux73322(2,2,2)+10*temp006332(2,2,2)*Z(jj,2)
     &  +4*temp006322(2,2,2)*Z(jj,3))
       temp73331(1,1,1)=IX*(aux73331(1,1,1)+8*temp006333(1,1,1)*Z(jj,1)+
     &  6*temp006331(1,1,1)*Z(jj,3))
       temp73332(1,1,1)=IX*(aux73332(1,1,1)+2*temp006333(1,1,1)*Z(jj,2)+
     &  6*(temp006333(2,1,1)*Z(jj,1)+temp006332(1,1,1)*Z(jj,3)))
       temp73332(2,1,1)=IX*(aux73332(2,1,1)+4*(temp006333(2,2,1)*Z(jj,1)
     &  +temp006333(2,1,1)*Z(jj,2))+6*temp006332(2,1,1)*Z(jj,3))
       temp73332(2,2,1)=IX*(aux73332(2,2,1)+2*temp006333(2,2,2)*Z(jj,1)+
     &  6*(temp006333(2,2,1)*Z(jj,2)+temp006332(2,2,1)*Z(jj,3)))
       temp73332(2,2,2)=IX*(aux73332(2,2,2)+8*temp006333(2,2,2)*Z(jj,2)+
     &  6*temp006332(2,2,2)*Z(jj,3))
       temp73333(1,1,1)=IX*(aux73333(1,1,1)+6*temp006333(3,1,1)*Z(jj,1)+
     &  8*temp006333(1,1,1)*Z(jj,3))
       temp73333(2,1,1)=IX*(aux73333(2,1,1)+4*temp006333(3,2,1)*Z(jj,1)+
     &  2*temp006333(3,1,1)*Z(jj,2)+8*temp006333(2,1,1)*Z(jj,3))
       temp73333(2,2,1)=IX*(aux73333(2,2,1)+2*temp006333(3,2,2)*Z(jj,1)+
     &  4*temp006333(3,2,1)*Z(jj,2)+8*temp006333(2,2,1)*Z(jj,3))
       temp73333(2,2,2)=IX*(aux73333(2,2,2)+6*temp006333(3,2,2)*Z(jj,2)+
     &  8*temp006333(2,2,2)*Z(jj,3))
       temp73333(3,1,1)=IX*(aux73333(3,1,1)+4*temp006333(3,3,1)*Z(jj,1)+
     &  10*temp006333(3,1,1)*Z(jj,3))
       temp73333(3,2,1)=IX*(aux73333(3,2,1)+2*(temp006333(3,3,2)*Z(jj,1)
     &  +temp006333(3,3,1)*Z(jj,2))+10*temp006333(3,2,1)*Z(jj,3))
       temp73333(3,2,2)=IX*(aux73333(3,2,2)+4*temp006333(3,3,2)*Z(jj,2)+
     &  10*temp006333(3,2,2)*Z(jj,3))
       temp73333(3,3,1)=IX*(aux73333(3,3,1)+2*temp006333(3,3,3)*Z(jj,1)+
     &  12*temp006333(3,3,1)*Z(jj,3))
       temp73333(3,3,2)=IX*(aux73333(3,3,2)+2*temp006333(3,3,3)*Z(jj,2)+
     &  12*temp006333(3,3,2)*Z(jj,3))
       temp73333(3,3,3)=IX*(aux73333(3,3,3)+14*temp006333(3,3,3)*Z(jj,3)
     &  )
       temp71111(1,1,2)=temp72111(1,1,1)
       temp71111(1,1,3)=temp73111(1,1,1)
       temp71111(1,2,1)=temp72111(1,1,1)
       temp71111(1,2,2)=temp72211(1,1,1)
       temp71111(1,2,3)=temp73211(1,1,1)
       temp71111(1,3,1)=temp73111(1,1,1)
       temp71111(1,3,2)=temp73211(1,1,1)
       temp71111(1,3,3)=temp73311(1,1,1)
       temp72111(1,1,2)=temp72211(1,1,1)
       temp72111(1,1,3)=temp73211(1,1,1)
       temp72111(1,2,1)=temp72211(1,1,1)
       temp72111(1,2,2)=temp72221(1,1,1)
       temp72111(1,2,3)=temp73221(1,1,1)
       temp72111(1,3,1)=temp73211(1,1,1)
       temp72111(1,3,2)=temp73221(1,1,1)
       temp72111(1,3,3)=temp73321(1,1,1)
       temp72211(1,1,2)=temp72221(1,1,1)
       temp72211(1,1,3)=temp73221(1,1,1)
       temp72211(1,2,1)=temp72221(1,1,1)
       temp72211(1,2,2)=temp72222(1,1,1)
       temp72211(1,2,3)=temp73222(1,1,1)
       temp72211(1,3,1)=temp73221(1,1,1)
       temp72211(1,3,2)=temp73222(1,1,1)
       temp72211(1,3,3)=temp73322(1,1,1)
       temp72221(1,1,2)=temp72222(1,1,1)
       temp72221(1,1,3)=temp73222(1,1,1)
       temp72221(1,2,1)=temp72222(1,1,1)
       temp72221(1,2,2)=temp72222(2,1,1)
       temp72221(1,2,3)=temp73222(2,1,1)
       temp72221(1,3,1)=temp73222(1,1,1)
       temp72221(1,3,2)=temp73222(2,1,1)
       temp72221(1,3,3)=temp73322(2,1,1)
       temp72222(1,1,2)=temp72222(2,1,1)
       temp72222(1,1,3)=temp73222(2,1,1)
       temp72222(1,2,1)=temp72222(2,1,1)
       temp72222(1,2,2)=temp72222(2,2,1)
       temp72222(1,2,3)=temp73222(2,2,1)
       temp72222(1,3,1)=temp73222(2,1,1)
       temp72222(1,3,2)=temp73222(2,2,1)
       temp72222(1,3,3)=temp73322(2,2,1)
       temp72222(2,1,2)=temp72222(2,2,1)
       temp72222(2,1,3)=temp73222(2,2,1)
       temp72222(2,2,3)=temp73222(2,2,2)
       temp72222(2,3,1)=temp73222(2,2,1)
       temp72222(2,3,2)=temp73222(2,2,2)
       temp72222(2,3,3)=temp73322(2,2,2)
       temp73111(1,1,2)=temp73211(1,1,1)
       temp73111(1,1,3)=temp73311(1,1,1)
       temp73111(1,2,1)=temp73211(1,1,1)
       temp73111(1,2,2)=temp73221(1,1,1)
       temp73111(1,2,3)=temp73321(1,1,1)
       temp73111(1,3,1)=temp73311(1,1,1)
       temp73111(1,3,2)=temp73321(1,1,1)
       temp73111(1,3,3)=temp73331(1,1,1)
       temp73211(1,1,2)=temp73221(1,1,1)
       temp73211(1,1,3)=temp73321(1,1,1)
       temp73211(1,2,1)=temp73221(1,1,1)
       temp73211(1,2,2)=temp73222(1,1,1)
       temp73211(1,2,3)=temp73322(1,1,1)
       temp73211(1,3,1)=temp73321(1,1,1)
       temp73211(1,3,2)=temp73322(1,1,1)
       temp73211(1,3,3)=temp73332(1,1,1)
       temp73221(1,1,2)=temp73222(1,1,1)
       temp73221(1,1,3)=temp73322(1,1,1)
       temp73221(1,2,1)=temp73222(1,1,1)
       temp73221(1,2,2)=temp73222(2,1,1)
       temp73221(1,2,3)=temp73322(2,1,1)
       temp73221(1,3,1)=temp73322(1,1,1)
       temp73221(1,3,2)=temp73322(2,1,1)
       temp73221(1,3,3)=temp73332(2,1,1)
       temp73222(1,1,2)=temp73222(2,1,1)
       temp73222(1,1,3)=temp73322(2,1,1)
       temp73222(1,2,1)=temp73222(2,1,1)
       temp73222(1,2,2)=temp73222(2,2,1)
       temp73222(1,2,3)=temp73322(2,2,1)
       temp73222(1,3,1)=temp73322(2,1,1)
       temp73222(1,3,2)=temp73322(2,2,1)
       temp73222(1,3,3)=temp73332(2,2,1)
       temp73222(2,1,2)=temp73222(2,2,1)
       temp73222(2,1,3)=temp73322(2,2,1)
       temp73222(2,2,3)=temp73322(2,2,2)
       temp73222(2,3,1)=temp73322(2,2,1)
       temp73222(2,3,2)=temp73322(2,2,2)
       temp73222(2,3,3)=temp73332(2,2,2)
       temp73311(1,1,2)=temp73321(1,1,1)
       temp73311(1,1,3)=temp73331(1,1,1)
       temp73311(1,2,1)=temp73321(1,1,1)
       temp73311(1,2,2)=temp73322(1,1,1)
       temp73311(1,2,3)=temp73332(1,1,1)
       temp73311(1,3,1)=temp73331(1,1,1)
       temp73311(1,3,2)=temp73332(1,1,1)
       temp73311(1,3,3)=temp73333(1,1,1)
       temp73321(1,1,2)=temp73322(1,1,1)
       temp73321(1,1,3)=temp73332(1,1,1)
       temp73321(1,2,1)=temp73322(1,1,1)
       temp73321(1,2,2)=temp73322(2,1,1)
       temp73321(1,2,3)=temp73332(2,1,1)
       temp73321(1,3,1)=temp73332(1,1,1)
       temp73321(1,3,2)=temp73332(2,1,1)
       temp73321(1,3,3)=temp73333(2,1,1)
       temp73322(1,1,2)=temp73322(2,1,1)
       temp73322(1,1,3)=temp73332(2,1,1)
       temp73322(1,2,1)=temp73322(2,1,1)
       temp73322(1,2,2)=temp73322(2,2,1)
       temp73322(1,2,3)=temp73332(2,2,1)
       temp73322(1,3,1)=temp73332(2,1,1)
       temp73322(1,3,2)=temp73332(2,2,1)
       temp73322(1,3,3)=temp73333(2,2,1)
       temp73322(2,1,2)=temp73322(2,2,1)
       temp73322(2,1,3)=temp73332(2,2,1)
       temp73322(2,2,3)=temp73332(2,2,2)
       temp73322(2,3,1)=temp73332(2,2,1)
       temp73322(2,3,2)=temp73332(2,2,2)
       temp73322(2,3,3)=temp73333(2,2,2)
       temp73331(1,1,2)=temp73332(1,1,1)
       temp73331(1,1,3)=temp73333(1,1,1)
       temp73331(1,2,1)=temp73332(1,1,1)
       temp73331(1,2,2)=temp73332(2,1,1)
       temp73331(1,2,3)=temp73333(2,1,1)
       temp73331(1,3,1)=temp73333(1,1,1)
       temp73331(1,3,2)=temp73333(2,1,1)
       temp73331(1,3,3)=temp73333(3,1,1)
       temp73332(1,1,2)=temp73332(2,1,1)
       temp73332(1,1,3)=temp73333(2,1,1)
       temp73332(1,2,1)=temp73332(2,1,1)
       temp73332(1,2,2)=temp73332(2,2,1)
       temp73332(1,2,3)=temp73333(2,2,1)
       temp73332(1,3,1)=temp73333(2,1,1)
       temp73332(1,3,2)=temp73333(2,2,1)
       temp73332(1,3,3)=temp73333(3,2,1)
       temp73332(2,1,2)=temp73332(2,2,1)
       temp73332(2,1,3)=temp73333(2,2,1)
       temp73332(2,2,3)=temp73333(2,2,2)
       temp73332(2,3,1)=temp73333(2,2,1)
       temp73332(2,3,2)=temp73333(2,2,2)
       temp73332(2,3,3)=temp73333(3,2,2)
       temp73333(1,1,2)=temp73333(2,1,1)
       temp73333(1,1,3)=temp73333(3,1,1)
       temp73333(1,2,1)=temp73333(2,1,1)
       temp73333(1,2,2)=temp73333(2,2,1)
       temp73333(1,2,3)=temp73333(3,2,1)
       temp73333(1,3,1)=temp73333(3,1,1)
       temp73333(1,3,2)=temp73333(3,2,1)
       temp73333(1,3,3)=temp73333(3,3,1)
       temp73333(2,1,2)=temp73333(2,2,1)
       temp73333(2,1,3)=temp73333(3,2,1)
       temp73333(2,2,3)=temp73333(3,2,2)
       temp73333(2,3,1)=temp73333(3,2,1)
       temp73333(2,3,2)=temp73333(3,2,2)
       temp73333(2,3,3)=temp73333(3,3,2)
       temp73333(3,1,2)=temp73333(3,2,1)
       temp73333(3,1,3)=temp73333(3,3,1)
       temp73333(3,2,3)=temp73333(3,3,2)
c                Step2
       temp0000001(1)=I16Z*(aux0000001(1)+2*tempD4000000*F(4)+F(5)*temp0
     &  0001(1)-det4*temp00003(1,k,l))
       temp0000001(2)=I16Z*(aux0000001(2)+2*tempD4000000*F(6)+F(5)*temp0
     &  0001(2)-det4*temp00003(2,k,l))
       temp0000001(3)=I16Z*(aux0000001(3)+2*tempD4000000*F(7)+F(5)*temp0
     &  0001(3)-det4*temp00003(3,k,l))
       temp00003(1,1,1)=I20Z*(aux00003(1,1,1)+6*F(4)*temp00002(1,1)+F(5)
     &  *temp003(1,1,1)-det4*temp00511(1,k,l)+24*temp0000001(1)*ZZ(k,1,l
     &  ,1))
       temp00003(2,1,1)=I20Z*(aux00003(2,1,1)+2*F(6)*temp00002(1,1)+4*F(
     &  4)*temp00002(2,1)+F(5)*temp003(2,1,1)-det4*temp00521(1,k,l)+8*(t
     &  emp0000001(2)*ZZ(k,1,l,1)+temp0000001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1
     &  ))))
       temp00003(2,2,1)=I20Z*(aux00003(2,2,1)+4*F(6)*temp00002(2,1)+2*F(
     &  4)*temp00002(2,2)+F(5)*temp003(2,2,1)-det4*temp00522(1,k,l)+8*(t
     &  emp0000001(2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0000001(1)*ZZ(k,2,l,
     &  2)))
       temp00003(2,2,2)=I20Z*(aux00003(2,2,2)+6*F(6)*temp00002(2,2)+F(5)
     &  *temp003(2,2,2)-det4*temp00522(2,k,l)+24*temp0000001(2)*ZZ(k,2,l
     &  ,2))
       temp00003(3,1,1)=I20Z*(aux00003(3,1,1)+2*F(7)*temp00002(1,1)+4*F(
     &  4)*temp00002(3,1)+F(5)*temp003(3,1,1)-det4*temp00531(1,k,l)+8*(t
     &  emp0000001(3)*ZZ(k,1,l,1)+temp0000001(1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1
     &  ))))
       temp00003(3,2,1)=I20Z*(aux00003(3,2,1)+2*F(7)*temp00002(2,1)+2*F(
     &  6)*temp00002(3,1)+2*F(4)*temp00002(3,2)+F(5)*temp003(3,2,1)-det4
     &  *temp00532(1,k,l)+4*(temp0000001(3)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+te
     &  mp0000001(2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)))+4*temp0000001(1)*(ZZ(k,2
     &  ,l,3)+ZZ(k,3,l,2)))
       temp00003(3,2,2)=I20Z*(aux00003(3,2,2)+2*F(7)*temp00002(2,2)+4*F(
     &  6)*temp00002(3,2)+F(5)*temp003(3,2,2)-det4*temp00532(2,k,l)+8*(t
     &  emp0000001(3)*ZZ(k,2,l,2)+temp0000001(2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2
     &  ))))
       temp00003(3,3,1)=I20Z*(aux00003(3,3,1)+4*F(7)*temp00002(3,1)+2*F(
     &  4)*temp00002(3,3)+F(5)*temp003(3,3,1)-det4*temp00533(1,k,l)+8*(t
     &  emp0000001(3)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+temp0000001(1)*ZZ(k,3,l,
     &  3)))
       temp00003(3,3,2)=I20Z*(aux00003(3,3,2)+4*F(7)*temp00002(3,2)+2*F(
     &  6)*temp00002(3,3)+F(5)*temp003(3,3,2)-det4*temp00533(2,k,l)+8*(t
     &  emp0000001(3)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+temp0000001(2)*ZZ(k,3,l,
     &  3)))
       temp00003(3,3,3)=I20Z*(aux00003(3,3,3)+6*F(7)*temp00002(3,3)+F(5)
     &  *temp003(3,3,3)-det4*temp00533(3,k,l)+24*temp0000001(3)*ZZ(k,3,l
     &  ,3))
       temp00003(1,1,2)=temp00003(2,1,1)
       temp00003(1,1,3)=temp00003(3,1,1)
       temp00003(1,2,1)=temp00003(2,1,1)
       temp00003(1,2,2)=temp00003(2,2,1)
       temp00003(1,2,3)=temp00003(3,2,1)
       temp00003(1,3,1)=temp00003(3,1,1)
       temp00003(1,3,2)=temp00003(3,2,1)
       temp00003(1,3,3)=temp00003(3,3,1)
       temp00003(2,1,2)=temp00003(2,2,1)
       temp00003(2,1,3)=temp00003(3,2,1)
       temp00003(2,2,3)=temp00003(3,2,2)
       temp00003(2,3,1)=temp00003(3,2,1)
       temp00003(2,3,2)=temp00003(3,2,2)
       temp00003(2,3,3)=temp00003(3,3,2)
       temp00003(3,1,2)=temp00003(3,2,1)
       temp00003(3,1,3)=temp00003(3,3,1)
       temp00003(3,2,3)=temp00003(3,3,2)
       temp00511(1,1,1)=I24Z*(aux00511(1,1,1)+10*F(4)*temp0041(1,1,1)+F(
     &  5)*temp511(1,1,1)-det4*temp71111(1,k,l)+80*temp00003(1,1,1)*ZZ(k
     &  ,1,l,1))
       temp00521(1,1,1)=I24Z*(aux00521(1,1,1)+2*F(6)*temp0041(1,1,1)+8*F
     &  (4)*temp0042(1,1,1)+F(5)*temp521(1,1,1)-det4*temp72111(1,k,l)+48
     &  *temp00003(2,1,1)*ZZ(k,1,l,1)+16*temp00003(1,1,1)*(ZZ(k,1,l,2)+Z
     &  Z(k,2,l,1)))
       temp00522(1,1,1)=I24Z*(aux00522(1,1,1)+4*F(6)*temp0042(1,1,1)+6*F
     &  (4)*temp0042(2,1,1)+F(5)*temp522(1,1,1)-det4*temp72211(1,k,l)+24
     &  *(temp00003(2,2,1)*ZZ(k,1,l,1)+temp00003(2,1,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))+8*temp00003(1,1,1)*ZZ(k,2,l,2))
       temp00522(2,1,1)=I24Z*(aux00522(2,1,1)+6*F(6)*temp0042(2,1,1)+4*F
     &  (4)*temp0042(2,2,1)+F(5)*temp522(2,1,1)-det4*temp72221(1,k,l)+8*
     &  temp00003(2,2,2)*ZZ(k,1,l,1)+24*(temp00003(2,2,1)*(ZZ(k,1,l,2)+Z
     &  Z(k,2,l,1))+temp00003(2,1,1)*ZZ(k,2,l,2)))
       temp00522(2,2,1)=I24Z*(aux00522(2,2,1)+8*F(6)*temp0042(2,2,1)+2*F
     &  (4)*temp0042(2,2,2)+F(5)*temp522(2,2,1)-det4*temp72222(1,k,l)+16
     &  *temp00003(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*temp00003(2,2,1)*
     &  ZZ(k,2,l,2))
       temp00522(2,2,2)=I24Z*(aux00522(2,2,2)+10*F(6)*temp0042(2,2,2)+F(
     &  5)*temp522(2,2,2)-det4*temp72222(2,k,l)+80*temp00003(2,2,2)*ZZ(k
     &  ,2,l,2))
       temp00531(1,1,1)=I24Z*(aux00531(1,1,1)+2*F(7)*temp0041(1,1,1)+8*F
     &  (4)*temp0043(1,1,1)+F(5)*temp531(1,1,1)-det4*temp73111(1,k,l)+48
     &  *temp00003(3,1,1)*ZZ(k,1,l,1)+16*temp00003(1,1,1)*(ZZ(k,1,l,3)+Z
     &  Z(k,3,l,1)))
       temp00532(1,1,1)=I24Z*(aux00532(1,1,1)+2*F(7)*temp0042(1,1,1)+2*F
     &  (6)*temp0043(1,1,1)+6*F(4)*temp0043(2,1,1)+F(5)*temp532(1,1,1)-d
     &  et4*temp73211(1,k,l)+24*temp00003(3,2,1)*ZZ(k,1,l,1)+12*(temp000
     &  03(3,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00003(2,1,1)*(ZZ(k,1,l,3
     &  )+ZZ(k,3,l,1)))+4*temp00003(1,1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp00532(2,1,1)=I24Z*(aux00532(2,1,1)+2*F(7)*temp0042(2,1,1)+4*F
     &  (6)*temp0043(2,1,1)+4*F(4)*temp0043(2,2,1)+F(5)*temp532(2,1,1)-d
     &  et4*temp73221(1,k,l)+16*temp00003(3,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1
     &  ))+8*(temp00003(3,2,2)*ZZ(k,1,l,1)+temp00003(3,1,1)*ZZ(k,2,l,2))
     &  +8*temp00003(2,2,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*temp00003(2,1,1)
     &  *(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp00532(2,2,1)=I24Z*(aux00532(2,2,1)+2*F(7)*temp0042(2,2,1)+6*F
     &  (6)*temp0043(2,2,1)+2*F(4)*temp0043(2,2,2)+F(5)*temp532(2,2,1)-d
     &  et4*temp73222(1,k,l)+24*temp00003(3,2,1)*ZZ(k,2,l,2)+4*temp00003
     &  (2,2,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+12*(temp00003(3,2,2)*(ZZ(k,1,l
     &  ,2)+ZZ(k,2,l,1))+temp00003(2,2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))))
       temp00532(2,2,2)=I24Z*(aux00532(2,2,2)+2*F(7)*temp0042(2,2,2)+8*F
     &  (6)*temp0043(2,2,2)+F(5)*temp532(2,2,2)-det4*temp73222(2,k,l)+48
     &  *temp00003(3,2,2)*ZZ(k,2,l,2)+16*temp00003(2,2,2)*(ZZ(k,2,l,3)+Z
     &  Z(k,3,l,2)))
       temp00533(1,1,1)=I24Z*(aux00533(1,1,1)+4*F(7)*temp0043(1,1,1)+6*F
     &  (4)*temp0043(3,1,1)+F(5)*temp533(1,1,1)-det4*temp73311(1,k,l)+24
     &  *(temp00003(3,3,1)*ZZ(k,1,l,1)+temp00003(3,1,1)*(ZZ(k,1,l,3)+ZZ(
     &  k,3,l,1)))+8*temp00003(1,1,1)*ZZ(k,3,l,3))
       temp00533(2,1,1)=I24Z*(aux00533(2,1,1)+4*F(7)*temp0043(2,1,1)+2*F
     &  (6)*temp0043(3,1,1)+4*F(4)*temp0043(3,2,1)+F(5)*temp533(2,1,1)-d
     &  et4*temp73321(1,k,l)+8*(temp00003(3,3,2)*ZZ(k,1,l,1)+temp00003(3
     &  ,3,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))+16*temp00003(3,2,1)*(ZZ(k,1,l,3
     &  )+ZZ(k,3,l,1))+8*temp00003(3,1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*te
     &  mp00003(2,1,1)*ZZ(k,3,l,3))
       temp00533(2,2,1)=I24Z*(aux00533(2,2,1)+4*F(7)*temp0043(2,2,1)+4*F
     &  (6)*temp0043(3,2,1)+2*F(4)*temp0043(3,2,2)+F(5)*temp533(2,2,1)-d
     &  et4*temp73322(1,k,l)+8*(temp00003(3,3,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1
     &  ))+temp00003(3,3,1)*ZZ(k,2,l,2))+8*temp00003(3,2,2)*(ZZ(k,1,l,3)
     &  +ZZ(k,3,l,1))+16*temp00003(3,2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*te
     &  mp00003(2,2,1)*ZZ(k,3,l,3))
       temp00533(2,2,2)=I24Z*(aux00533(2,2,2)+4*F(7)*temp0043(2,2,2)+6*F
     &  (6)*temp0043(3,2,2)+F(5)*temp533(2,2,2)-det4*temp73322(2,k,l)+24
     &  *(temp00003(3,3,2)*ZZ(k,2,l,2)+temp00003(3,2,2)*(ZZ(k,2,l,3)+ZZ(
     &  k,3,l,2)))+8*temp00003(2,2,2)*ZZ(k,3,l,3))
       temp00533(3,1,1)=I24Z*(aux00533(3,1,1)+6*F(7)*temp0043(3,1,1)+4*F
     &  (4)*temp0043(3,3,1)+F(5)*temp533(3,1,1)-det4*temp73331(1,k,l)+8*
     &  temp00003(3,3,3)*ZZ(k,1,l,1)+24*(temp00003(3,3,1)*(ZZ(k,1,l,3)+Z
     &  Z(k,3,l,1))+temp00003(3,1,1)*ZZ(k,3,l,3)))
       temp00533(3,2,1)=I24Z*(aux00533(3,2,1)+6*F(7)*temp0043(3,2,1)+2*F
     &  (6)*temp0043(3,3,1)+2*F(4)*temp0043(3,3,2)+F(5)*temp533(3,2,1)-d
     &  et4*temp73332(1,k,l)+4*temp00003(3,3,3)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)
     &  )+12*(temp00003(3,3,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+temp00003(3,3,1
     &  )*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))+24*temp00003(3,2,1)*ZZ(k,3,l,3))
       temp00533(3,2,2)=I24Z*(aux00533(3,2,2)+6*F(7)*temp0043(3,2,2)+4*F
     &  (6)*temp0043(3,3,2)+F(5)*temp533(3,2,2)-det4*temp73332(2,k,l)+8*
     &  temp00003(3,3,3)*ZZ(k,2,l,2)+24*(temp00003(3,3,2)*(ZZ(k,2,l,3)+Z
     &  Z(k,3,l,2))+temp00003(3,2,2)*ZZ(k,3,l,3)))
       temp00533(3,3,1)=I24Z*(aux00533(3,3,1)+8*F(7)*temp0043(3,3,1)+2*F
     &  (4)*temp0043(3,3,3)+F(5)*temp533(3,3,1)-det4*temp73333(1,k,l)+16
     &  *temp00003(3,3,3)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+48*temp00003(3,3,1)*
     &  ZZ(k,3,l,3))
       temp00533(3,3,2)=I24Z*(aux00533(3,3,2)+8*F(7)*temp0043(3,3,2)+2*F
     &  (6)*temp0043(3,3,3)+F(5)*temp533(3,3,2)-det4*temp73333(2,k,l)+16
     &  *temp00003(3,3,3)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+48*temp00003(3,3,2)*
     &  ZZ(k,3,l,3))
       temp00533(3,3,3)=I24Z*(aux00533(3,3,3)+10*F(7)*temp0043(3,3,3)+F(
     &  5)*temp533(3,3,3)-det4*temp73333(3,k,l)+80*temp00003(3,3,3)*ZZ(k
     &  ,3,l,3))
       temp00511(1,1,2)=temp00521(1,1,1)
       temp00511(1,1,3)=temp00531(1,1,1)
       temp00511(1,2,1)=temp00521(1,1,1)
       temp00511(1,2,2)=temp00522(1,1,1)
       temp00511(1,2,3)=temp00532(1,1,1)
       temp00511(1,3,1)=temp00531(1,1,1)
       temp00511(1,3,2)=temp00532(1,1,1)
       temp00511(1,3,3)=temp00533(1,1,1)
       temp00521(1,1,2)=temp00522(1,1,1)
       temp00521(1,1,3)=temp00532(1,1,1)
       temp00521(1,2,1)=temp00522(1,1,1)
       temp00521(1,2,2)=temp00522(2,1,1)
       temp00521(1,2,3)=temp00532(2,1,1)
       temp00521(1,3,1)=temp00532(1,1,1)
       temp00521(1,3,2)=temp00532(2,1,1)
       temp00521(1,3,3)=temp00533(2,1,1)
       temp00522(1,1,2)=temp00522(2,1,1)
       temp00522(1,1,3)=temp00532(2,1,1)
       temp00522(1,2,1)=temp00522(2,1,1)
       temp00522(1,2,2)=temp00522(2,2,1)
       temp00522(1,2,3)=temp00532(2,2,1)
       temp00522(1,3,1)=temp00532(2,1,1)
       temp00522(1,3,2)=temp00532(2,2,1)
       temp00522(1,3,3)=temp00533(2,2,1)
       temp00522(2,1,2)=temp00522(2,2,1)
       temp00522(2,1,3)=temp00532(2,2,1)
       temp00522(2,2,3)=temp00532(2,2,2)
       temp00522(2,3,1)=temp00532(2,2,1)
       temp00522(2,3,2)=temp00532(2,2,2)
       temp00522(2,3,3)=temp00533(2,2,2)
       temp00531(1,1,2)=temp00532(1,1,1)
       temp00531(1,1,3)=temp00533(1,1,1)
       temp00531(1,2,1)=temp00532(1,1,1)
       temp00531(1,2,2)=temp00532(2,1,1)
       temp00531(1,2,3)=temp00533(2,1,1)
       temp00531(1,3,1)=temp00533(1,1,1)
       temp00531(1,3,2)=temp00533(2,1,1)
       temp00531(1,3,3)=temp00533(3,1,1)
       temp00532(1,1,2)=temp00532(2,1,1)
       temp00532(1,1,3)=temp00533(2,1,1)
       temp00532(1,2,1)=temp00532(2,1,1)
       temp00532(1,2,2)=temp00532(2,2,1)
       temp00532(1,2,3)=temp00533(2,2,1)
       temp00532(1,3,1)=temp00533(2,1,1)
       temp00532(1,3,2)=temp00533(2,2,1)
       temp00532(1,3,3)=temp00533(3,2,1)
       temp00532(2,1,2)=temp00532(2,2,1)
       temp00532(2,1,3)=temp00533(2,2,1)
       temp00532(2,2,3)=temp00533(2,2,2)
       temp00532(2,3,1)=temp00533(2,2,1)
       temp00532(2,3,2)=temp00533(2,2,2)
       temp00532(2,3,3)=temp00533(3,2,2)
       temp00533(1,1,2)=temp00533(2,1,1)
       temp00533(1,1,3)=temp00533(3,1,1)
       temp00533(1,2,1)=temp00533(2,1,1)
       temp00533(1,2,2)=temp00533(2,2,1)
       temp00533(1,2,3)=temp00533(3,2,1)
       temp00533(1,3,1)=temp00533(3,1,1)
       temp00533(1,3,2)=temp00533(3,2,1)
       temp00533(1,3,3)=temp00533(3,3,1)
       temp00533(2,1,2)=temp00533(2,2,1)
       temp00533(2,1,3)=temp00533(3,2,1)
       temp00533(2,2,3)=temp00533(3,2,2)
       temp00533(2,3,1)=temp00533(3,2,1)
       temp00533(2,3,2)=temp00533(3,2,2)
       temp00533(2,3,3)=temp00533(3,3,2)
       temp00533(3,1,2)=temp00533(3,2,1)
       temp00533(3,1,3)=temp00533(3,3,1)
       temp00533(3,2,3)=temp00533(3,3,2)
       temp6111(1,1,1)=IX*(aux6111(1,1,1)+det4*temp71111(1,1,jj)+12*temp
     &  00511(1,1,1)*Z(jj,1))
       temp6211(1,1,1)=IX*(aux6211(1,1,1)+det4*temp72111(1,1,jj)+10*temp
     &  00521(1,1,1)*Z(jj,1)+2*temp00511(1,1,1)*Z(jj,2))
       temp6221(1,1,1)=IX*(aux6221(1,1,1)+det4*temp72211(1,1,jj)+8*temp0
     &  0522(1,1,1)*Z(jj,1)+4*temp00521(1,1,1)*Z(jj,2))
       temp6222(1,1,1)=IX*(aux6222(1,1,1)+det4*temp72221(1,1,jj)+6*(temp
     &  00522(2,1,1)*Z(jj,1)+temp00522(1,1,1)*Z(jj,2)))
       temp6222(2,1,1)=IX*(aux6222(2,1,1)+det4*temp72222(1,1,jj)+4*temp0
     &  0522(2,2,1)*Z(jj,1)+8*temp00522(2,1,1)*Z(jj,2))
       temp6222(2,2,1)=IX*(aux6222(2,2,1)+det4*temp72222(2,1,jj)+2*temp0
     &  0522(2,2,2)*Z(jj,1)+10*temp00522(2,2,1)*Z(jj,2))
       temp6222(2,2,2)=IX*(aux6222(2,2,2)+det4*temp72222(2,2,jj)+12*temp
     &  00522(2,2,2)*Z(jj,2))
       temp6311(1,1,1)=IX*(aux6311(1,1,1)+det4*temp73111(1,1,jj)+10*temp
     &  00531(1,1,1)*Z(jj,1)+2*temp00511(1,1,1)*Z(jj,3))
       temp6321(1,1,1)=IX*(aux6321(1,1,1)+det4*temp73211(1,1,jj)+8*temp0
     &  0532(1,1,1)*Z(jj,1)+2*(temp00531(1,1,1)*Z(jj,2)+temp00521(1,1,1)
     &  *Z(jj,3)))
       temp6322(1,1,1)=IX*(aux6322(1,1,1)+det4*temp73221(1,1,jj)+6*temp0
     &  0532(2,1,1)*Z(jj,1)+4*temp00532(1,1,1)*Z(jj,2)+2*temp00522(1,1,1
     &  )*Z(jj,3))
       temp6322(2,1,1)=IX*(aux6322(2,1,1)+det4*temp73222(1,1,jj)+4*temp0
     &  0532(2,2,1)*Z(jj,1)+6*temp00532(2,1,1)*Z(jj,2)+2*temp00522(2,1,1
     &  )*Z(jj,3))
       temp6322(2,2,1)=IX*(aux6322(2,2,1)+det4*temp73222(2,1,jj)+8*temp0
     &  0532(2,2,1)*Z(jj,2)+2*(temp00532(2,2,2)*Z(jj,1)+temp00522(2,2,1)
     &  *Z(jj,3)))
       temp6322(2,2,2)=IX*(aux6322(2,2,2)+det4*temp73222(2,2,jj)+10*temp
     &  00532(2,2,2)*Z(jj,2)+2*temp00522(2,2,2)*Z(jj,3))
       temp6331(1,1,1)=IX*(aux6331(1,1,1)+det4*temp73311(1,1,jj)+8*temp0
     &  0533(1,1,1)*Z(jj,1)+4*temp00531(1,1,1)*Z(jj,3))
       temp6332(1,1,1)=IX*(aux6332(1,1,1)+det4*temp73321(1,1,jj)+6*temp0
     &  0533(2,1,1)*Z(jj,1)+2*temp00533(1,1,1)*Z(jj,2)+4*temp00532(1,1,1
     &  )*Z(jj,3))
       temp6332(2,1,1)=IX*(aux6332(2,1,1)+det4*temp73322(1,1,jj)+4*(temp
     &  00533(2,2,1)*Z(jj,1)+temp00533(2,1,1)*Z(jj,2))+4*temp00532(2,1,1
     &  )*Z(jj,3))
       temp6332(2,2,1)=IX*(aux6332(2,2,1)+det4*temp73322(2,1,jj)+2*temp0
     &  0533(2,2,2)*Z(jj,1)+6*temp00533(2,2,1)*Z(jj,2)+4*temp00532(2,2,1
     &  )*Z(jj,3))
       temp6332(2,2,2)=IX*(aux6332(2,2,2)+det4*temp73322(2,2,jj)+8*temp0
     &  0533(2,2,2)*Z(jj,2)+4*temp00532(2,2,2)*Z(jj,3))
       temp6333(1,1,1)=IX*(aux6333(1,1,1)+det4*temp73331(1,1,jj)+6*(temp
     &  00533(3,1,1)*Z(jj,1)+temp00533(1,1,1)*Z(jj,3)))
       temp6333(2,1,1)=IX*(aux6333(2,1,1)+det4*temp73332(1,1,jj)+4*temp0
     &  0533(3,2,1)*Z(jj,1)+2*temp00533(3,1,1)*Z(jj,2)+6*temp00533(2,1,1
     &  )*Z(jj,3))
       temp6333(2,2,1)=IX*(aux6333(2,2,1)+det4*temp73332(2,1,jj)+2*temp0
     &  0533(3,2,2)*Z(jj,1)+4*temp00533(3,2,1)*Z(jj,2)+6*temp00533(2,2,1
     &  )*Z(jj,3))
       temp6333(2,2,2)=IX*(aux6333(2,2,2)+det4*temp73332(2,2,jj)+6*(temp
     &  00533(3,2,2)*Z(jj,2)+temp00533(2,2,2)*Z(jj,3)))
       temp6333(3,1,1)=IX*(aux6333(3,1,1)+det4*temp73333(1,1,jj)+4*temp0
     &  0533(3,3,1)*Z(jj,1)+8*temp00533(3,1,1)*Z(jj,3))
       temp6333(3,2,1)=IX*(aux6333(3,2,1)+det4*temp73333(2,1,jj)+2*(temp
     &  00533(3,3,2)*Z(jj,1)+temp00533(3,3,1)*Z(jj,2))+8*temp00533(3,2,1
     &  )*Z(jj,3))
       temp6333(3,2,2)=IX*(aux6333(3,2,2)+det4*temp73333(2,2,jj)+4*temp0
     &  0533(3,3,2)*Z(jj,2)+8*temp00533(3,2,2)*Z(jj,3))
       temp6333(3,3,1)=IX*(aux6333(3,3,1)+det4*temp73333(3,1,jj)+2*temp0
     &  0533(3,3,3)*Z(jj,1)+10*temp00533(3,3,1)*Z(jj,3))
       temp6333(3,3,2)=IX*(aux6333(3,3,2)+det4*temp73333(3,2,jj)+2*temp0
     &  0533(3,3,3)*Z(jj,2)+10*temp00533(3,3,2)*Z(jj,3))
       temp6333(3,3,3)=IX*(aux6333(3,3,3)+det4*temp73333(3,3,jj)+12*temp
     &  00533(3,3,3)*Z(jj,3))
       temp6111(1,1,2)=temp6211(1,1,1)
       temp6111(1,1,3)=temp6311(1,1,1)
       temp6111(1,2,1)=temp6211(1,1,1)
       temp6111(1,2,2)=temp6221(1,1,1)
       temp6111(1,2,3)=temp6321(1,1,1)
       temp6111(1,3,1)=temp6311(1,1,1)
       temp6111(1,3,2)=temp6321(1,1,1)
       temp6111(1,3,3)=temp6331(1,1,1)
       temp6211(1,1,2)=temp6221(1,1,1)
       temp6211(1,1,3)=temp6321(1,1,1)
       temp6211(1,2,1)=temp6221(1,1,1)
       temp6211(1,2,2)=temp6222(1,1,1)
       temp6211(1,2,3)=temp6322(1,1,1)
       temp6211(1,3,1)=temp6321(1,1,1)
       temp6211(1,3,2)=temp6322(1,1,1)
       temp6211(1,3,3)=temp6332(1,1,1)
       temp6221(1,1,2)=temp6222(1,1,1)
       temp6221(1,1,3)=temp6322(1,1,1)
       temp6221(1,2,1)=temp6222(1,1,1)
       temp6221(1,2,2)=temp6222(2,1,1)
       temp6221(1,2,3)=temp6322(2,1,1)
       temp6221(1,3,1)=temp6322(1,1,1)
       temp6221(1,3,2)=temp6322(2,1,1)
       temp6221(1,3,3)=temp6332(2,1,1)
       temp6222(1,1,2)=temp6222(2,1,1)
       temp6222(1,1,3)=temp6322(2,1,1)
       temp6222(1,2,1)=temp6222(2,1,1)
       temp6222(1,2,2)=temp6222(2,2,1)
       temp6222(1,2,3)=temp6322(2,2,1)
       temp6222(1,3,1)=temp6322(2,1,1)
       temp6222(1,3,2)=temp6322(2,2,1)
       temp6222(1,3,3)=temp6332(2,2,1)
       temp6222(2,1,2)=temp6222(2,2,1)
       temp6222(2,1,3)=temp6322(2,2,1)
       temp6222(2,2,3)=temp6322(2,2,2)
       temp6222(2,3,1)=temp6322(2,2,1)
       temp6222(2,3,2)=temp6322(2,2,2)
       temp6222(2,3,3)=temp6332(2,2,2)
       temp6311(1,1,2)=temp6321(1,1,1)
       temp6311(1,1,3)=temp6331(1,1,1)
       temp6311(1,2,1)=temp6321(1,1,1)
       temp6311(1,2,2)=temp6322(1,1,1)
       temp6311(1,2,3)=temp6332(1,1,1)
       temp6311(1,3,1)=temp6331(1,1,1)
       temp6311(1,3,2)=temp6332(1,1,1)
       temp6311(1,3,3)=temp6333(1,1,1)
       temp6321(1,1,2)=temp6322(1,1,1)
       temp6321(1,1,3)=temp6332(1,1,1)
       temp6321(1,2,1)=temp6322(1,1,1)
       temp6321(1,2,2)=temp6322(2,1,1)
       temp6321(1,2,3)=temp6332(2,1,1)
       temp6321(1,3,1)=temp6332(1,1,1)
       temp6321(1,3,2)=temp6332(2,1,1)
       temp6321(1,3,3)=temp6333(2,1,1)
       temp6322(1,1,2)=temp6322(2,1,1)
       temp6322(1,1,3)=temp6332(2,1,1)
       temp6322(1,2,1)=temp6322(2,1,1)
       temp6322(1,2,2)=temp6322(2,2,1)
       temp6322(1,2,3)=temp6332(2,2,1)
       temp6322(1,3,1)=temp6332(2,1,1)
       temp6322(1,3,2)=temp6332(2,2,1)
       temp6322(1,3,3)=temp6333(2,2,1)
       temp6322(2,1,2)=temp6322(2,2,1)
       temp6322(2,1,3)=temp6332(2,2,1)
       temp6322(2,2,3)=temp6332(2,2,2)
       temp6322(2,3,1)=temp6332(2,2,1)
       temp6322(2,3,2)=temp6332(2,2,2)
       temp6322(2,3,3)=temp6333(2,2,2)
       temp6331(1,1,2)=temp6332(1,1,1)
       temp6331(1,1,3)=temp6333(1,1,1)
       temp6331(1,2,1)=temp6332(1,1,1)
       temp6331(1,2,2)=temp6332(2,1,1)
       temp6331(1,2,3)=temp6333(2,1,1)
       temp6331(1,3,1)=temp6333(1,1,1)
       temp6331(1,3,2)=temp6333(2,1,1)
       temp6331(1,3,3)=temp6333(3,1,1)
       temp6332(1,1,2)=temp6332(2,1,1)
       temp6332(1,1,3)=temp6333(2,1,1)
       temp6332(1,2,1)=temp6332(2,1,1)
       temp6332(1,2,2)=temp6332(2,2,1)
       temp6332(1,2,3)=temp6333(2,2,1)
       temp6332(1,3,1)=temp6333(2,1,1)
       temp6332(1,3,2)=temp6333(2,2,1)
       temp6332(1,3,3)=temp6333(3,2,1)
       temp6332(2,1,2)=temp6332(2,2,1)
       temp6332(2,1,3)=temp6333(2,2,1)
       temp6332(2,2,3)=temp6333(2,2,2)
       temp6332(2,3,1)=temp6333(2,2,1)
       temp6332(2,3,2)=temp6333(2,2,2)
       temp6332(2,3,3)=temp6333(3,2,2)
       temp6333(1,1,2)=temp6333(2,1,1)
       temp6333(1,1,3)=temp6333(3,1,1)
       temp6333(1,2,1)=temp6333(2,1,1)
       temp6333(1,2,2)=temp6333(2,2,1)
       temp6333(1,2,3)=temp6333(3,2,1)
       temp6333(1,3,1)=temp6333(3,1,1)
       temp6333(1,3,2)=temp6333(3,2,1)
       temp6333(1,3,3)=temp6333(3,3,1)
       temp6333(2,1,2)=temp6333(2,2,1)
       temp6333(2,1,3)=temp6333(3,2,1)
       temp6333(2,2,3)=temp6333(3,2,2)
       temp6333(2,3,1)=temp6333(3,2,1)
       temp6333(2,3,2)=temp6333(3,2,2)
       temp6333(2,3,3)=temp6333(3,3,2)
       temp6333(3,1,2)=temp6333(3,2,1)
       temp6333(3,1,3)=temp6333(3,3,1)
       temp6333(3,2,3)=temp6333(3,3,2)
c                Step3
       tempD4000000=I12Z*(auxD4000000+tempD40000*F(5)-det4*temp00002(k,l
     &  ))
       temp00002(1,1)=I16Z*(aux00002(1,1)+4*F(4)*temp00001(1)+F(5)*temp0
     &  02(1,1)-det4*temp0041(1,k,l)+8*tempD4000000*ZZ(k,1,l,1))
       temp00002(2,1)=I16Z*(aux00002(2,1)+2*(F(6)*temp00001(1)+F(4)*temp
     &  00001(2))+F(5)*temp002(2,1)-det4*temp0042(1,k,l)+4*tempD4000000*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp00002(2,2)=I16Z*(aux00002(2,2)+4*F(6)*temp00001(2)+F(5)*temp0
     &  02(2,2)-det4*temp0042(2,k,l)+8*tempD4000000*ZZ(k,2,l,2))
       temp00002(3,1)=I16Z*(aux00002(3,1)+2*(F(7)*temp00001(1)+F(4)*temp
     &  00001(3))+F(5)*temp002(3,1)-det4*temp0043(1,k,l)+4*tempD4000000*
     &  (ZZ(k,1,l,3)+ZZ(k,3,l,1)))
       temp00002(3,2)=I16Z*(aux00002(3,2)+2*(F(7)*temp00001(2)+F(6)*temp
     &  00001(3))+F(5)*temp002(3,2)-det4*temp0043(2,k,l)+4*tempD4000000*
     &  (ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp00002(3,3)=I16Z*(aux00002(3,3)+4*F(7)*temp00001(3)+F(5)*temp0
     &  02(3,3)-det4*temp0043(3,k,l)+8*tempD4000000*ZZ(k,3,l,3))
       temp00002(1,2)=temp00002(2,1)
       temp00002(1,3)=temp00002(3,1)
       temp00002(2,3)=temp00002(3,2)
       temp0041(1,1,1)=I20Z*(aux0041(1,1,1)+8*F(4)*temp003(1,1,1)+F(5)*t
     &  emp41(1,1,1)-det4*temp6111(1,k,l)+48*temp00002(1,1)*ZZ(k,1,l,1))
       temp0042(1,1,1)=I20Z*(aux0042(1,1,1)+2*F(6)*temp003(1,1,1)+6*F(4)
     &  *temp003(2,1,1)+F(5)*temp42(1,1,1)-det4*temp6211(1,k,l)+24*temp0
     &  0002(2,1)*ZZ(k,1,l,1)+12*temp00002(1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)
     &  ))
       temp0042(2,1,1)=I20Z*(aux0042(2,1,1)+4*F(6)*temp003(2,1,1)+4*F(4)
     &  *temp003(2,2,1)+F(5)*temp42(2,1,1)-det4*temp6221(1,k,l)+16*temp0
     &  0002(2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp00002(2,2)*ZZ(k,1,l,1
     &  )+temp00002(1,1)*ZZ(k,2,l,2)))
       temp0042(2,2,1)=I20Z*(aux0042(2,2,1)+6*F(6)*temp003(2,2,1)+2*F(4)
     &  *temp003(2,2,2)+F(5)*temp42(2,2,1)-det4*temp6222(1,k,l)+12*temp0
     &  0002(2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*temp00002(2,1)*ZZ(k,2,l,2
     &  ))
       temp0042(2,2,2)=I20Z*(aux0042(2,2,2)+8*F(6)*temp003(2,2,2)+F(5)*t
     &  emp42(2,2,2)-det4*temp6222(2,k,l)+48*temp00002(2,2)*ZZ(k,2,l,2))
       temp0043(1,1,1)=I20Z*(aux0043(1,1,1)+2*F(7)*temp003(1,1,1)+6*F(4)
     &  *temp003(3,1,1)+F(5)*temp43(1,1,1)-det4*temp6311(1,k,l)+24*temp0
     &  0002(3,1)*ZZ(k,1,l,1)+12*temp00002(1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)
     &  ))
       temp0043(2,1,1)=I20Z*(aux0043(2,1,1)+2*F(7)*temp003(2,1,1)+2*F(6)
     &  *temp003(3,1,1)+4*F(4)*temp003(3,2,1)+F(5)*temp43(2,1,1)-det4*te
     &  mp6321(1,k,l)+8*(temp00002(3,2)*ZZ(k,1,l,1)+temp00002(3,1)*(ZZ(k
     &  ,1,l,2)+ZZ(k,2,l,1)))+8*temp00002(2,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))
     &  +4*temp00002(1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp0043(2,2,1)=I20Z*(aux0043(2,2,1)+2*F(7)*temp003(2,2,1)+4*F(6)
     &  *temp003(3,2,1)+2*F(4)*temp003(3,2,2)+F(5)*temp43(2,2,1)-det4*te
     &  mp6322(1,k,l)+8*(temp00002(3,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00
     &  002(3,1)*ZZ(k,2,l,2))+4*temp00002(2,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))
     &  +8*temp00002(2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp0043(2,2,2)=I20Z*(aux0043(2,2,2)+2*F(7)*temp003(2,2,2)+6*F(6)
     &  *temp003(3,2,2)+F(5)*temp43(2,2,2)-det4*temp6322(2,k,l)+24*temp0
     &  0002(3,2)*ZZ(k,2,l,2)+12*temp00002(2,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)
     &  ))
       temp0043(3,1,1)=I20Z*(aux0043(3,1,1)+4*F(7)*temp003(3,1,1)+4*F(4)
     &  *temp003(3,3,1)+F(5)*temp43(3,1,1)-det4*temp6331(1,k,l)+16*temp0
     &  0002(3,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*(temp00002(3,3)*ZZ(k,1,l,1
     &  )+temp00002(1,1)*ZZ(k,3,l,3)))
       temp0043(3,2,1)=I20Z*(aux0043(3,2,1)+2*F(6)*temp003(3,3,1)+2*F(4)
     &  *temp003(3,3,2)+F(5)*temp43(3,2,1)-det4*temp6332(1,k,l)+4*(F(7)*
     &  temp003(3,2,1)+temp00002(3,3)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))+8*temp0
     &  0002(3,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*temp00002(3,1)*(ZZ(k,2,l,3
     &  )+ZZ(k,3,l,2))+8*temp00002(2,1)*ZZ(k,3,l,3))
       temp0043(3,2,2)=I20Z*(aux0043(3,2,2)+4*F(7)*temp003(3,2,2)+4*F(6)
     &  *temp003(3,3,2)+F(5)*temp43(3,2,2)-det4*temp6332(2,k,l)+16*temp0
     &  0002(3,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*(temp00002(3,3)*ZZ(k,2,l,2
     &  )+temp00002(2,2)*ZZ(k,3,l,3)))
       temp0043(3,3,1)=I20Z*(aux0043(3,3,1)+6*F(7)*temp003(3,3,1)+2*F(4)
     &  *temp003(3,3,3)+F(5)*temp43(3,3,1)-det4*temp6333(1,k,l)+12*temp0
     &  0002(3,3)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+24*temp00002(3,1)*ZZ(k,3,l,3
     &  ))
       temp0043(3,3,2)=I20Z*(aux0043(3,3,2)+6*F(7)*temp003(3,3,2)+2*F(6)
     &  *temp003(3,3,3)+F(5)*temp43(3,3,2)-det4*temp6333(2,k,l)+12*temp0
     &  0002(3,3)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+24*temp00002(3,2)*ZZ(k,3,l,3
     &  ))
       temp0043(3,3,3)=I20Z*(aux0043(3,3,3)+8*F(7)*temp003(3,3,3)+F(5)*t
     &  emp43(3,3,3)-det4*temp6333(3,k,l)+48*temp00002(3,3)*ZZ(k,3,l,3))
       temp0041(1,1,2)=temp0042(1,1,1)
       temp0041(1,1,3)=temp0043(1,1,1)
       temp0041(1,2,1)=temp0042(1,1,1)
       temp0041(1,2,2)=temp0042(2,1,1)
       temp0041(1,2,3)=temp0043(2,1,1)
       temp0041(1,3,1)=temp0043(1,1,1)
       temp0041(1,3,2)=temp0043(2,1,1)
       temp0041(1,3,3)=temp0043(3,1,1)
       temp0042(1,1,2)=temp0042(2,1,1)
       temp0042(1,1,3)=temp0043(2,1,1)
       temp0042(1,2,1)=temp0042(2,1,1)
       temp0042(1,2,2)=temp0042(2,2,1)
       temp0042(1,2,3)=temp0043(2,2,1)
       temp0042(1,3,1)=temp0043(2,1,1)
       temp0042(1,3,2)=temp0043(2,2,1)
       temp0042(1,3,3)=temp0043(3,2,1)
       temp0042(2,1,2)=temp0042(2,2,1)
       temp0042(2,1,3)=temp0043(2,2,1)
       temp0042(2,2,3)=temp0043(2,2,2)
       temp0042(2,3,1)=temp0043(2,2,1)
       temp0042(2,3,2)=temp0043(2,2,2)
       temp0042(2,3,3)=temp0043(3,2,2)
       temp0043(1,1,2)=temp0043(2,1,1)
       temp0043(1,1,3)=temp0043(3,1,1)
       temp0043(1,2,1)=temp0043(2,1,1)
       temp0043(1,2,2)=temp0043(2,2,1)
       temp0043(1,2,3)=temp0043(3,2,1)
       temp0043(1,3,1)=temp0043(3,1,1)
       temp0043(1,3,2)=temp0043(3,2,1)
       temp0043(1,3,3)=temp0043(3,3,1)
       temp0043(2,1,2)=temp0043(2,2,1)
       temp0043(2,1,3)=temp0043(3,2,1)
       temp0043(2,2,3)=temp0043(3,2,2)
       temp0043(2,3,1)=temp0043(3,2,1)
       temp0043(2,3,2)=temp0043(3,2,2)
       temp0043(2,3,3)=temp0043(3,3,2)
       temp0043(3,1,2)=temp0043(3,2,1)
       temp0043(3,1,3)=temp0043(3,3,1)
       temp0043(3,2,3)=temp0043(3,3,2)
       temp511(1,1,1)=IX*(aux511(1,1,1)+det4*temp6111(1,1,jj)+10*temp004
     &  1(1,1,1)*Z(jj,1))
       temp521(1,1,1)=IX*(aux521(1,1,1)+det4*temp6211(1,1,jj)+8*temp0042
     &  (1,1,1)*Z(jj,1)+2*temp0041(1,1,1)*Z(jj,2))
       temp522(1,1,1)=IX*(aux522(1,1,1)+det4*temp6221(1,1,jj)+6*temp0042
     &  (2,1,1)*Z(jj,1)+4*temp0042(1,1,1)*Z(jj,2))
       temp522(2,1,1)=IX*(aux522(2,1,1)+det4*temp6222(1,1,jj)+4*temp0042
     &  (2,2,1)*Z(jj,1)+6*temp0042(2,1,1)*Z(jj,2))
       temp522(2,2,1)=IX*(aux522(2,2,1)+det4*temp6222(2,1,jj)+2*temp0042
     &  (2,2,2)*Z(jj,1)+8*temp0042(2,2,1)*Z(jj,2))
       temp522(2,2,2)=IX*(aux522(2,2,2)+det4*temp6222(2,2,jj)+10*temp004
     &  2(2,2,2)*Z(jj,2))
       temp531(1,1,1)=IX*(aux531(1,1,1)+det4*temp6311(1,1,jj)+8*temp0043
     &  (1,1,1)*Z(jj,1)+2*temp0041(1,1,1)*Z(jj,3))
       temp532(1,1,1)=IX*(aux532(1,1,1)+det4*temp6321(1,1,jj)+6*temp0043
     &  (2,1,1)*Z(jj,1)+2*(temp0043(1,1,1)*Z(jj,2)+temp0042(1,1,1)*Z(jj,
     &  3)))
       temp532(2,1,1)=IX*(aux532(2,1,1)+det4*temp6322(1,1,jj)+4*(temp004
     &  3(2,2,1)*Z(jj,1)+temp0043(2,1,1)*Z(jj,2))+2*temp0042(2,1,1)*Z(jj
     &  ,3))
       temp532(2,2,1)=IX*(aux532(2,2,1)+det4*temp6322(2,1,jj)+6*temp0043
     &  (2,2,1)*Z(jj,2)+2*(temp0043(2,2,2)*Z(jj,1)+temp0042(2,2,1)*Z(jj,
     &  3)))
       temp532(2,2,2)=IX*(aux532(2,2,2)+det4*temp6322(2,2,jj)+8*temp0043
     &  (2,2,2)*Z(jj,2)+2*temp0042(2,2,2)*Z(jj,3))
       temp533(1,1,1)=IX*(aux533(1,1,1)+det4*temp6331(1,1,jj)+6*temp0043
     &  (3,1,1)*Z(jj,1)+4*temp0043(1,1,1)*Z(jj,3))
       temp533(2,1,1)=IX*(aux533(2,1,1)+det4*temp6332(1,1,jj)+2*temp0043
     &  (3,1,1)*Z(jj,2)+4*(temp0043(3,2,1)*Z(jj,1)+temp0043(2,1,1)*Z(jj,
     &  3)))
       temp533(2,2,1)=IX*(aux533(2,2,1)+det4*temp6332(2,1,jj)+2*temp0043
     &  (3,2,2)*Z(jj,1)+4*(temp0043(3,2,1)*Z(jj,2)+temp0043(2,2,1)*Z(jj,
     &  3)))
       temp533(2,2,2)=IX*(aux533(2,2,2)+det4*temp6332(2,2,jj)+6*temp0043
     &  (3,2,2)*Z(jj,2)+4*temp0043(2,2,2)*Z(jj,3))
       temp533(3,1,1)=IX*(aux533(3,1,1)+det4*temp6333(1,1,jj)+4*temp0043
     &  (3,3,1)*Z(jj,1)+6*temp0043(3,1,1)*Z(jj,3))
       temp533(3,2,1)=IX*(aux533(3,2,1)+det4*temp6333(2,1,jj)+2*(temp004
     &  3(3,3,2)*Z(jj,1)+temp0043(3,3,1)*Z(jj,2))+6*temp0043(3,2,1)*Z(jj
     &  ,3))
       temp533(3,2,2)=IX*(aux533(3,2,2)+det4*temp6333(2,2,jj)+4*temp0043
     &  (3,3,2)*Z(jj,2)+6*temp0043(3,2,2)*Z(jj,3))
       temp533(3,3,1)=IX*(aux533(3,3,1)+det4*temp6333(3,1,jj)+2*temp0043
     &  (3,3,3)*Z(jj,1)+8*temp0043(3,3,1)*Z(jj,3))
       temp533(3,3,2)=IX*(aux533(3,3,2)+det4*temp6333(3,2,jj)+2*temp0043
     &  (3,3,3)*Z(jj,2)+8*temp0043(3,3,2)*Z(jj,3))
       temp533(3,3,3)=IX*(aux533(3,3,3)+det4*temp6333(3,3,jj)+10*temp004
     &  3(3,3,3)*Z(jj,3))
       temp511(1,1,2)=temp521(1,1,1)
       temp511(1,1,3)=temp531(1,1,1)
       temp511(1,2,1)=temp521(1,1,1)
       temp511(1,2,2)=temp522(1,1,1)
       temp511(1,2,3)=temp532(1,1,1)
       temp511(1,3,1)=temp531(1,1,1)
       temp511(1,3,2)=temp532(1,1,1)
       temp511(1,3,3)=temp533(1,1,1)
       temp521(1,1,2)=temp522(1,1,1)
       temp521(1,1,3)=temp532(1,1,1)
       temp521(1,2,1)=temp522(1,1,1)
       temp521(1,2,2)=temp522(2,1,1)
       temp521(1,2,3)=temp532(2,1,1)
       temp521(1,3,1)=temp532(1,1,1)
       temp521(1,3,2)=temp532(2,1,1)
       temp521(1,3,3)=temp533(2,1,1)
       temp522(1,1,2)=temp522(2,1,1)
       temp522(1,1,3)=temp532(2,1,1)
       temp522(1,2,1)=temp522(2,1,1)
       temp522(1,2,2)=temp522(2,2,1)
       temp522(1,2,3)=temp532(2,2,1)
       temp522(1,3,1)=temp532(2,1,1)
       temp522(1,3,2)=temp532(2,2,1)
       temp522(1,3,3)=temp533(2,2,1)
       temp522(2,1,2)=temp522(2,2,1)
       temp522(2,1,3)=temp532(2,2,1)
       temp522(2,2,3)=temp532(2,2,2)
       temp522(2,3,1)=temp532(2,2,1)
       temp522(2,3,2)=temp532(2,2,2)
       temp522(2,3,3)=temp533(2,2,2)
       temp531(1,1,2)=temp532(1,1,1)
       temp531(1,1,3)=temp533(1,1,1)
       temp531(1,2,1)=temp532(1,1,1)
       temp531(1,2,2)=temp532(2,1,1)
       temp531(1,2,3)=temp533(2,1,1)
       temp531(1,3,1)=temp533(1,1,1)
       temp531(1,3,2)=temp533(2,1,1)
       temp531(1,3,3)=temp533(3,1,1)
       temp532(1,1,2)=temp532(2,1,1)
       temp532(1,1,3)=temp533(2,1,1)
       temp532(1,2,1)=temp532(2,1,1)
       temp532(1,2,2)=temp532(2,2,1)
       temp532(1,2,3)=temp533(2,2,1)
       temp532(1,3,1)=temp533(2,1,1)
       temp532(1,3,2)=temp533(2,2,1)
       temp532(1,3,3)=temp533(3,2,1)
       temp532(2,1,2)=temp532(2,2,1)
       temp532(2,1,3)=temp533(2,2,1)
       temp532(2,2,3)=temp533(2,2,2)
       temp532(2,3,1)=temp533(2,2,1)
       temp532(2,3,2)=temp533(2,2,2)
       temp532(2,3,3)=temp533(3,2,2)
       temp533(1,1,2)=temp533(2,1,1)
       temp533(1,1,3)=temp533(3,1,1)
       temp533(1,2,1)=temp533(2,1,1)
       temp533(1,2,2)=temp533(2,2,1)
       temp533(1,2,3)=temp533(3,2,1)
       temp533(1,3,1)=temp533(3,1,1)
       temp533(1,3,2)=temp533(3,2,1)
       temp533(1,3,3)=temp533(3,3,1)
       temp533(2,1,2)=temp533(2,2,1)
       temp533(2,1,3)=temp533(3,2,1)
       temp533(2,2,3)=temp533(3,2,2)
       temp533(2,3,1)=temp533(3,2,1)
       temp533(2,3,2)=temp533(3,2,2)
       temp533(2,3,3)=temp533(3,3,2)
       temp533(3,1,2)=temp533(3,2,1)
       temp533(3,1,3)=temp533(3,3,1)
       temp533(3,2,3)=temp533(3,3,2)
c                Step4
       temp00001(1)=I12Z*(aux00001(1)+2*tempD40000*F(4)+F(5)*temp001(1)-
     &  det4*temp003(1,k,l))
       temp00001(2)=I12Z*(aux00001(2)+2*tempD40000*F(6)+F(5)*temp001(2)-
     &  det4*temp003(2,k,l))
       temp00001(3)=I12Z*(aux00001(3)+2*tempD40000*F(7)+F(5)*temp001(3)-
     &  det4*temp003(3,k,l))
       temp003(1,1,1)=I16Z*(aux003(1,1,1)+6*F(4)*temp002(1,1)+F(5)*temp3
     &  (1,1,1)-det4*temp511(1,k,l)+24*temp00001(1)*ZZ(k,1,l,1))
       temp003(2,1,1)=I16Z*(aux003(2,1,1)+2*F(6)*temp002(1,1)+4*F(4)*tem
     &  p002(2,1)+F(5)*temp3(2,1,1)-det4*temp521(1,k,l)+8*(temp00001(2)*
     &  ZZ(k,1,l,1)+temp00001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))))
       temp003(2,2,1)=I16Z*(aux003(2,2,1)+4*F(6)*temp002(2,1)+2*F(4)*tem
     &  p002(2,2)+F(5)*temp3(2,2,1)-det4*temp522(1,k,l)+8*(temp00001(2)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00001(1)*ZZ(k,2,l,2)))
       temp003(2,2,2)=I16Z*(aux003(2,2,2)+6*F(6)*temp002(2,2)+F(5)*temp3
     &  (2,2,2)-det4*temp522(2,k,l)+24*temp00001(2)*ZZ(k,2,l,2))
       temp003(3,1,1)=I16Z*(aux003(3,1,1)+2*F(7)*temp002(1,1)+4*F(4)*tem
     &  p002(3,1)+F(5)*temp3(3,1,1)-det4*temp531(1,k,l)+8*(temp00001(3)*
     &  ZZ(k,1,l,1)+temp00001(1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))))
       temp003(3,2,1)=I16Z*(aux003(3,2,1)+2*F(7)*temp002(2,1)+2*F(6)*tem
     &  p002(3,1)+2*F(4)*temp002(3,2)+F(5)*temp3(3,2,1)-det4*temp532(1,k
     &  ,l)+4*(temp00001(3)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00001(2)*(ZZ(k
     &  ,1,l,3)+ZZ(k,3,l,1)))+4*temp00001(1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp003(3,2,2)=I16Z*(aux003(3,2,2)+2*F(7)*temp002(2,2)+4*F(6)*tem
     &  p002(3,2)+F(5)*temp3(3,2,2)-det4*temp532(2,k,l)+8*(temp00001(3)*
     &  ZZ(k,2,l,2)+temp00001(2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))))
       temp003(3,3,1)=I16Z*(aux003(3,3,1)+4*F(7)*temp002(3,1)+2*F(4)*tem
     &  p002(3,3)+F(5)*temp3(3,3,1)-det4*temp533(1,k,l)+8*(temp00001(3)*
     &  (ZZ(k,1,l,3)+ZZ(k,3,l,1))+temp00001(1)*ZZ(k,3,l,3)))
       temp003(3,3,2)=I16Z*(aux003(3,3,2)+4*F(7)*temp002(3,2)+2*F(6)*tem
     &  p002(3,3)+F(5)*temp3(3,3,2)-det4*temp533(2,k,l)+8*(temp00001(3)*
     &  (ZZ(k,2,l,3)+ZZ(k,3,l,2))+temp00001(2)*ZZ(k,3,l,3)))
       temp003(3,3,3)=I16Z*(aux003(3,3,3)+6*F(7)*temp002(3,3)+F(5)*temp3
     &  (3,3,3)-det4*temp533(3,k,l)+24*temp00001(3)*ZZ(k,3,l,3))
       temp003(1,1,2)=temp003(2,1,1)
       temp003(1,1,3)=temp003(3,1,1)
       temp003(1,2,1)=temp003(2,1,1)
       temp003(1,2,2)=temp003(2,2,1)
       temp003(1,2,3)=temp003(3,2,1)
       temp003(1,3,1)=temp003(3,1,1)
       temp003(1,3,2)=temp003(3,2,1)
       temp003(1,3,3)=temp003(3,3,1)
       temp003(2,1,2)=temp003(2,2,1)
       temp003(2,1,3)=temp003(3,2,1)
       temp003(2,2,3)=temp003(3,2,2)
       temp003(2,3,1)=temp003(3,2,1)
       temp003(2,3,2)=temp003(3,2,2)
       temp003(2,3,3)=temp003(3,3,2)
       temp003(3,1,2)=temp003(3,2,1)
       temp003(3,1,3)=temp003(3,3,1)
       temp003(3,2,3)=temp003(3,3,2)
       temp41(1,1,1)=IX*(aux41(1,1,1)+det4*temp511(1,1,jj)+8*temp003(1,1
     &  ,1)*Z(jj,1))
       temp42(1,1,1)=IX*(aux42(1,1,1)+det4*temp521(1,1,jj)+6*temp003(2,1
     &  ,1)*Z(jj,1)+2*temp003(1,1,1)*Z(jj,2))
       temp42(2,1,1)=IX*(aux42(2,1,1)+det4*temp522(1,1,jj)+4*(temp003(2,
     &  2,1)*Z(jj,1)+temp003(2,1,1)*Z(jj,2)))
       temp42(2,2,1)=IX*(aux42(2,2,1)+det4*temp522(2,1,jj)+2*temp003(2,2
     &  ,2)*Z(jj,1)+6*temp003(2,2,1)*Z(jj,2))
       temp42(2,2,2)=IX*(aux42(2,2,2)+det4*temp522(2,2,jj)+8*temp003(2,2
     &  ,2)*Z(jj,2))
       temp43(1,1,1)=IX*(aux43(1,1,1)+det4*temp531(1,1,jj)+6*temp003(3,1
     &  ,1)*Z(jj,1)+2*temp003(1,1,1)*Z(jj,3))
       temp43(2,1,1)=IX*(aux43(2,1,1)+det4*temp532(1,1,jj)+4*temp003(3,2
     &  ,1)*Z(jj,1)+2*(temp003(3,1,1)*Z(jj,2)+temp003(2,1,1)*Z(jj,3)))
       temp43(2,2,1)=IX*(aux43(2,2,1)+det4*temp532(2,1,jj)+4*temp003(3,2
     &  ,1)*Z(jj,2)+2*(temp003(3,2,2)*Z(jj,1)+temp003(2,2,1)*Z(jj,3)))
       temp43(2,2,2)=IX*(aux43(2,2,2)+det4*temp532(2,2,jj)+6*temp003(3,2
     &  ,2)*Z(jj,2)+2*temp003(2,2,2)*Z(jj,3))
       temp43(3,1,1)=IX*(aux43(3,1,1)+det4*temp533(1,1,jj)+4*(temp003(3,
     &  3,1)*Z(jj,1)+temp003(3,1,1)*Z(jj,3)))
       temp43(3,2,1)=IX*(aux43(3,2,1)+det4*temp533(2,1,jj)+2*(temp003(3,
     &  3,2)*Z(jj,1)+temp003(3,3,1)*Z(jj,2))+4*temp003(3,2,1)*Z(jj,3))
       temp43(3,2,2)=IX*(aux43(3,2,2)+det4*temp533(2,2,jj)+4*(temp003(3,
     &  3,2)*Z(jj,2)+temp003(3,2,2)*Z(jj,3)))
       temp43(3,3,1)=IX*(aux43(3,3,1)+det4*temp533(3,1,jj)+2*temp003(3,3
     &  ,3)*Z(jj,1)+6*temp003(3,3,1)*Z(jj,3))
       temp43(3,3,2)=IX*(aux43(3,3,2)+det4*temp533(3,2,jj)+2*temp003(3,3
     &  ,3)*Z(jj,2)+6*temp003(3,3,2)*Z(jj,3))
       temp43(3,3,3)=IX*(aux43(3,3,3)+det4*temp533(3,3,jj)+8*temp003(3,3
     &  ,3)*Z(jj,3))
       temp41(1,1,2)=temp42(1,1,1)
       temp41(1,1,3)=temp43(1,1,1)
       temp41(1,2,1)=temp42(1,1,1)
       temp41(1,2,2)=temp42(2,1,1)
       temp41(1,2,3)=temp43(2,1,1)
       temp41(1,3,1)=temp43(1,1,1)
       temp41(1,3,2)=temp43(2,1,1)
       temp41(1,3,3)=temp43(3,1,1)
       temp42(1,1,2)=temp42(2,1,1)
       temp42(1,1,3)=temp43(2,1,1)
       temp42(1,2,1)=temp42(2,1,1)
       temp42(1,2,2)=temp42(2,2,1)
       temp42(1,2,3)=temp43(2,2,1)
       temp42(1,3,1)=temp43(2,1,1)
       temp42(1,3,2)=temp43(2,2,1)
       temp42(1,3,3)=temp43(3,2,1)
       temp42(2,1,2)=temp42(2,2,1)
       temp42(2,1,3)=temp43(2,2,1)
       temp42(2,2,3)=temp43(2,2,2)
       temp42(2,3,1)=temp43(2,2,1)
       temp42(2,3,2)=temp43(2,2,2)
       temp42(2,3,3)=temp43(3,2,2)
       temp43(1,1,2)=temp43(2,1,1)
       temp43(1,1,3)=temp43(3,1,1)
       temp43(1,2,1)=temp43(2,1,1)
       temp43(1,2,2)=temp43(2,2,1)
       temp43(1,2,3)=temp43(3,2,1)
       temp43(1,3,1)=temp43(3,1,1)
       temp43(1,3,2)=temp43(3,2,1)
       temp43(1,3,3)=temp43(3,3,1)
       temp43(2,1,2)=temp43(2,2,1)
       temp43(2,1,3)=temp43(3,2,1)
       temp43(2,2,3)=temp43(3,2,2)
       temp43(2,3,1)=temp43(3,2,1)
       temp43(2,3,2)=temp43(3,2,2)
       temp43(2,3,3)=temp43(3,3,2)
       temp43(3,1,2)=temp43(3,2,1)
       temp43(3,1,3)=temp43(3,3,1)
       temp43(3,2,3)=temp43(3,3,2)
c                Step5
       tempD40000=I8Z*(auxD40000+tempD400*F(5)-det4*temp002(k,l))
       temp002(1,1)=I12Z*(aux002(1,1)+4*F(4)*temp001(1)+F(5)*temp2(1,1)-
     &  det4*temp41(1,k,l)+8*tempD40000*ZZ(k,1,l,1))
       temp002(2,1)=I12Z*(aux002(2,1)+2*(F(6)*temp001(1)+F(4)*temp001(2)
     &  )+F(5)*temp2(2,1)-det4*temp42(1,k,l)+4*tempD40000*(ZZ(k,1,l,2)+Z
     &  Z(k,2,l,1)))
       temp002(2,2)=I12Z*(aux002(2,2)+4*F(6)*temp001(2)+F(5)*temp2(2,2)-
     &  det4*temp42(2,k,l)+8*tempD40000*ZZ(k,2,l,2))
       temp002(3,1)=I12Z*(aux002(3,1)+2*(F(7)*temp001(1)+F(4)*temp001(3)
     &  )+F(5)*temp2(3,1)-det4*temp43(1,k,l)+4*tempD40000*(ZZ(k,1,l,3)+Z
     &  Z(k,3,l,1)))
       temp002(3,2)=I12Z*(aux002(3,2)+2*(F(7)*temp001(2)+F(6)*temp001(3)
     &  )+F(5)*temp2(3,2)-det4*temp43(2,k,l)+4*tempD40000*(ZZ(k,2,l,3)+Z
     &  Z(k,3,l,2)))
       temp002(3,3)=I12Z*(aux002(3,3)+4*F(7)*temp001(3)+F(5)*temp2(3,3)-
     &  det4*temp43(3,k,l)+8*tempD40000*ZZ(k,3,l,3))
       temp002(1,2)=temp002(2,1)
       temp002(1,3)=temp002(3,1)
       temp002(2,3)=temp002(3,2)
       temp3(1,1,1)=IX*(aux3(1,1,1)+det4*temp41(1,1,jj)+6*temp002(1,1)*Z
     &  (jj,1))
       temp3(2,1,1)=IX*(aux3(2,1,1)+det4*temp42(1,1,jj)+4*temp002(2,1)*Z
     &  (jj,1)+2*temp002(1,1)*Z(jj,2))
       temp3(2,2,1)=IX*(aux3(2,2,1)+det4*temp42(2,1,jj)+2*temp002(2,2)*Z
     &  (jj,1)+4*temp002(2,1)*Z(jj,2))
       temp3(2,2,2)=IX*(aux3(2,2,2)+det4*temp42(2,2,jj)+6*temp002(2,2)*Z
     &  (jj,2))
       temp3(3,1,1)=IX*(aux3(3,1,1)+det4*temp43(1,1,jj)+4*temp002(3,1)*Z
     &  (jj,1)+2*temp002(1,1)*Z(jj,3))
       temp3(3,2,1)=IX*(aux3(3,2,1)+det4*temp43(2,1,jj)+2*(temp002(3,2)*
     &  Z(jj,1)+temp002(3,1)*Z(jj,2))+2*temp002(2,1)*Z(jj,3))
       temp3(3,2,2)=IX*(aux3(3,2,2)+det4*temp43(2,2,jj)+4*temp002(3,2)*Z
     &  (jj,2)+2*temp002(2,2)*Z(jj,3))
       temp3(3,3,1)=IX*(aux3(3,3,1)+det4*temp43(3,1,jj)+2*temp002(3,3)*Z
     &  (jj,1)+4*temp002(3,1)*Z(jj,3))
       temp3(3,3,2)=IX*(aux3(3,3,2)+det4*temp43(3,2,jj)+2*temp002(3,3)*Z
     &  (jj,2)+4*temp002(3,2)*Z(jj,3))
       temp3(3,3,3)=IX*(aux3(3,3,3)+det4*temp43(3,3,jj)+6*temp002(3,3)*Z
     &  (jj,3))
       temp3(1,1,2)=temp3(2,1,1)
       temp3(1,1,3)=temp3(3,1,1)
       temp3(1,2,1)=temp3(2,1,1)
       temp3(1,2,2)=temp3(2,2,1)
       temp3(1,2,3)=temp3(3,2,1)
       temp3(1,3,1)=temp3(3,1,1)
       temp3(1,3,2)=temp3(3,2,1)
       temp3(1,3,3)=temp3(3,3,1)
       temp3(2,1,2)=temp3(2,2,1)
       temp3(2,1,3)=temp3(3,2,1)
       temp3(2,2,3)=temp3(3,2,2)
       temp3(2,3,1)=temp3(3,2,1)
       temp3(2,3,2)=temp3(3,2,2)
       temp3(2,3,3)=temp3(3,3,2)
       temp3(3,1,2)=temp3(3,2,1)
       temp3(3,1,3)=temp3(3,3,1)
       temp3(3,2,3)=temp3(3,3,2)
c                Step6
       temp001(1)=I8Z*(aux001(1)+2*tempD400*F(4)+F(5)*temp1(1)-det4*temp
     &  3(1,k,l))
       temp001(2)=I8Z*(aux001(2)+2*tempD400*F(6)+F(5)*temp1(2)-det4*temp
     &  3(2,k,l))
       temp001(3)=I8Z*(aux001(3)+2*tempD400*F(7)+F(5)*temp1(3)-det4*temp
     &  3(3,k,l))
       temp2(1,1)=IX*(aux2(1,1)+det4*temp3(1,1,jj)+4*temp001(1)*Z(jj,1))
       temp2(2,1)=IX*(aux2(2,1)+det4*temp3(2,1,jj)+2*temp001(2)*Z(jj,1)+
     &  2*temp001(1)*Z(jj,2))
       temp2(2,2)=IX*(aux2(2,2)+det4*temp3(2,2,jj)+4*temp001(2)*Z(jj,2))
       temp2(3,1)=IX*(aux2(3,1)+det4*temp3(3,1,jj)+2*temp001(3)*Z(jj,1)+
     &  2*temp001(1)*Z(jj,3))
       temp2(3,2)=IX*(aux2(3,2)+det4*temp3(3,2,jj)+2*temp001(3)*Z(jj,2)+
     &  2*temp001(2)*Z(jj,3))
       temp2(3,3)=IX*(aux2(3,3)+det4*temp3(3,3,jj)+4*temp001(3)*Z(jj,3))
       temp2(1,2)=temp2(2,1)
       temp2(1,3)=temp2(3,1)
       temp2(2,3)=temp2(3,2)
c                Step7
       tempD400=I4Z*(auxD400+tempD40*F(5)-det4*temp2(k,l))
       temp1(1)=IX*(aux1(1)+det4*temp2(1,jj)+2*tempD400*Z(jj,1))
       temp1(2)=IX*(aux1(2)+det4*temp2(2,jj)+2*tempD400*Z(jj,2))
       temp1(3)=IX*(aux1(3)+det4*temp2(3,jj)+2*tempD400*Z(jj,3))
c                Step8
       tempD40=IX*(auxD40+det4*temp1(jj))

        ac=2
        accuracyDR(1,0,ac)  = abs(D0        /tempD40         -1d0)
        accuracyDR(1,1,ac)  = abs(Dij(1,1)  /temp1(1)        -1d0)
        accuracyDR(2,1,ac)  = abs(Dij(2,1)  /temp1(2)        -1d0)
        accuracyDR(3,1,ac)  = abs(Dij(3,1)  /temp1(3)        -1d0)
        accuracyDR(1,2,ac)  = abs(Dij(1,2)  /temp2(1,1)      -1d0)
        accuracyDR(2,2,ac)  = abs(Dij(2,2)  /temp2(2,2)      -1d0)
        accuracyDR(3,2,ac)  = abs(Dij(3,2)  /temp2(3,3)      -1d0)
        accuracyDR(4,2,ac)  = abs(Dij(4,2)  /temp2(2,1)      -1d0)
        accuracyDR(5,2,ac)  = abs(Dij(5,2)  /temp2(3,1)      -1d0)
        accuracyDR(6,2,ac)  = abs(Dij(6,2)  /temp2(3,2)      -1d0)
        accuracyDR(7,2,ac)  = abs(Dij(7,2)  /tempD400        -1d0)
        accuracyDR(1,3,ac)  = abs(Dij(1,3)  /temp3(1,1,1)    -1d0)
        accuracyDR(2,3,ac)  = abs(Dij(2,3)  /temp3(2,2,2)    -1d0)
        accuracyDR(3,3,ac)  = abs(Dij(3,3)  /temp3(3,3,3)    -1d0)
        accuracyDR(4,3,ac)  = abs(Dij(4,3)  /temp3(2,1,1)    -1d0)
        accuracyDR(5,3,ac)  = abs(Dij(5,3)  /temp3(3,1,1)    -1d0)
        accuracyDR(6,3,ac)  = abs(Dij(6,3)  /temp3(2,2,1)    -1d0)
        accuracyDR(7,3,ac)  = abs(Dij(7,3)  /temp3(3,3,1)    -1d0)
        accuracyDR(8,3,ac)  = abs(Dij(8,3)  /temp3(3,2,2)    -1d0)
        accuracyDR(9,3,ac)  = abs(Dij(9,3)  /temp3(3,3,2)    -1d0)
        accuracyDR(10,3,ac) = abs(Dij(10,3) /temp3(3,2,1)    -1d0)
        accuracyDR(11,3,ac) = abs(Dij(11,3) /temp001(1)      -1d0)
        accuracyDR(12,3,ac) = abs(Dij(12,3) /temp001(2)      -1d0)
        accuracyDR(13,3,ac) = abs(Dij(13,3) /temp001(3)      -1d0)
        accuracyDR(7,1,ac)  = abs(Dij(7,1)  /temp002(1,1)    -1d0)
        accuracyDR(8,1,ac)  = abs(Dij(8,1)  /temp002(2,2)    -1d0)
        accuracyDR(9,1,ac)  = abs(Dij(9,1)  /temp002(3,3)    -1d0)
        accuracyDR(10,1,ac) = abs(Dij(10,1) /temp002(2,1)    -1d0)
        accuracyDR(11,1,ac) = abs(Dij(11,1) /temp002(3,1)    -1d0)
        accuracyDR(12,1,ac) = abs(Dij(12,1) /temp002(3,2)    -1d0)
        accuracyDR(13,1,ac) = abs(Dij(13,1) /tempD40000      -1d0)
        accuracyDR(1,4,ac)  = abs(Dij(1,4)  /temp41(1,1,1)   -1d0)
        accuracyDR(2,4,ac)  = abs(Dij(2,4)  /temp42(2,2,2)   -1d0)
        accuracyDR(3,4,ac)  = abs(Dij(3,4)  /temp43(3,3,3)   -1d0)
        accuracyDR(4,4,ac)  = abs(Dij(4,4)  /temp42(1,1,1)   -1d0)
        accuracyDR(5,4,ac)  = abs(Dij(5,4)  /temp43(1,1,1)   -1d0)
        accuracyDR(6,4,ac)  = abs(Dij(6,4)  /temp42(2,1,1)   -1d0)
        accuracyDR(7,4,ac)  = abs(Dij(7,4)  /temp43(2,1,1)   -1d0)
        accuracyDR(8,4,ac)  = abs(Dij(8,4)  /temp43(3,1,1)   -1d0)
        accuracyDR(9,4,ac)  = abs(Dij(9,4)  /temp42(2,2,1)   -1d0)
        accuracyDR(10,4,ac) = abs(Dij(10,4) /temp43(2,2,1)   -1d0)
        accuracyDR(11,4,ac) = abs(Dij(11,4) /temp43(3,2,1)   -1d0)
        accuracyDR(12,4,ac) = abs(Dij(12,4) /temp43(3,3,1)   -1d0)
        accuracyDR(13,4,ac) = abs(Dij(13,4) /temp43(2,2,2)   -1d0)
        accuracyDR(14,4,ac) = abs(Dij(14,4) /temp43(3,2,2)   -1d0)
        accuracyDR(15,4,ac) = abs(Dij(15,4) /temp43(3,3,2)   -1d0)
        accuracyDR(16,4,ac) = abs(Dij(16,4) /temp002(1,1)    -1d0)
        accuracyDR(17,4,ac) = abs(Dij(17,4) /temp002(2,2)    -1d0)
        accuracyDR(18,4,ac) = abs(Dij(18,4) /temp002(3,3)    -1d0)
        accuracyDR(19,4,ac) = abs(Dij(19,4) /temp002(2,1)    -1d0)
        accuracyDR(20,4,ac) = abs(Dij(20,4) /temp002(3,1)    -1d0)
        accuracyDR(21,4,ac) = abs(Dij(21,4) /temp002(3,2)    -1d0)
        accuracyDR(22,4,ac) = abs(Dij(22,4) /tempD40000      -1d0)


      DO I1=0,4
           accuracyD(i1,ac)=accuracyDR(1,i1,ac)
            DO I2=2,INDEX(I1)  
       if(accuracyDR(i2,i1,ac).gt.accuracyD(i1,ac)) then
          accuracyD(i1,ac)=accuracyDR(i2,i1,ac)
       endif
          enddo
        enddo

        D0=tempD40
        Dij(1,1)=temp1(1)
        Dij(2,1)=temp1(2)
        Dij(3,1)=temp1(3)
        Dij(1,2)=temp2(1,1)
        Dij(2,2)=temp2(2,2)
        Dij(3,2)=temp2(3,3)
        Dij(4,2)=temp2(2,1)
        Dij(5,2)=temp2(3,1)
        Dij(6,2)=temp2(3,2)
        Dij(7,2)=tempD400
        Dij(1,3)=temp3(1,1,1)
        Dij(2,3)=temp3(2,2,2)
        Dij(3,3)=temp3(3,3,3)
        Dij(4,3)=temp3(2,1,1)
        Dij(5,3)=temp3(3,1,1)
        Dij(6,3)=temp3(2,2,1)
        Dij(7,3)=temp3(3,3,1)
        Dij(8,3)=temp3(3,2,2)
        Dij(9,3)=temp3(3,3,2)
        Dij(10,3)=temp3(3,2,1)
        Dij(11,3)=temp001(1)
        Dij(12,3)=temp001(2)
        Dij(13,3)=temp001(3)
        Dij(7,1)=temp002(1,1)
        Dij(8,1)=temp002(2,2)
        Dij(9,1)=temp002(3,3)
        Dij(10,1)=temp002(2,1)
        Dij(11,1)=temp002(3,1)
        Dij(12,1)=temp002(3,2)
        Dij(13,1)=tempD40000
        Dij(1,4)=temp41(1,1,1)
         Dij(2,4)=temp42(2,2,2)
         Dij(3,4)=temp43(3,3,3)
         Dij(4,4)=temp42(1,1,1)
         Dij(5,4)=temp43(1,1,1)
         Dij(6,4)=temp42(2,1,1)
         Dij(7,4)=temp43(2,1,1)
         Dij(8,4)=temp43(3,1,1)
         Dij(9,4)=temp42(2,2,1)
         Dij(10,4)=temp43(2,2,1)
         Dij(11,4)=temp43(3,2,1)
         Dij(12,4)=temp43(3,3,1)
         Dij(13,4)=temp43(2,2,2)
         Dij(14,4)=temp43(3,2,2)
         Dij(15,4)=temp43(3,3,2)
         Dij(16,4)=temp002(1,1)
         Dij(17,4)=temp002(2,2)
         Dij(18,4)=temp002(3,3)
         Dij(19,4)=temp002(2,1)
         Dij(20,4)=temp002(3,1)
         Dij(21,4)=temp002(3,2)
         Dij(22,4)=tempD40000
         Dij(1,5)=temp511(1,1,1)
         Dij(2,5)=temp522(2,2,2)
         Dij(3,5)=temp533(3,3,3)
         Dij(4,5)=temp521(1,1,1)
         Dij(5,5)=temp531(1,1,1)
         Dij(6,5)=temp522(1,1,1)
         Dij(7,5)=temp532(1,1,1)
         Dij(8,5)=temp533(1,1,1)
         Dij(9,5)=temp522(2,1,1)
         Dij(10,5)=temp532(2,1,1)
         Dij(11,5)=temp533(2,1,1)
         Dij(12,5)=temp533(3,1,1)
         Dij(13,5)=temp522(2,2,1)
         Dij(14,5)=temp532(2,2,1)
         Dij(15,5)=temp533(2,2,1)
         Dij(16,5)=temp533(3,2,1)
         Dij(17,5)=temp533(3,3,1)
         Dij(18,5)=temp532(2,2,2)
         Dij(19,5)=temp533(2,2,2)
         Dij(20,5)=temp533(3,2,2)
         Dij(21,5)=temp533(3,3,2)
         Dij(22,5)=temp003(1,1,1)
         Dij(23,5)=temp003(2,2,2)
         Dij(24,5)=temp003(3,3,3)
         Dij(25,5)=temp003(2,1,1)
         Dij(26,5)=temp003(3,1,1)
         Dij(27,5)=temp003(2,2,1)
         Dij(28,5)=temp003(3,2,1)
         Dij(29,5)=temp003(3,3,1)
         Dij(30,5)=temp003(3,2,2)
         Dij(31,5)=temp003(3,3,2)
         Dij(32,5)=temp00001(1)
         Dij(33,5)=temp00001(2)
         Dij(34,5)=temp00001(3)

       if(order.eq.7) goto 500
c                Iteration8
c                Step1
       S3000000001(1)=-2*Cij234(16,6)
       S3000000001(2)=2*Cij234(19,7)
       S3000000001(3)=2*Cij234(20,7)
       S3h000000001(1)=Cij134(25,8)-Cij234(25,8)
       S3h000000001(2)=Cij124(25,8)-Cij134(25,8)
       S3h000000001(3)=Cij123(25,8)-Cij124(25,8)
       S3h000000311(1)=Cij134(22,8)-Cij234(16,6)
       S3h000000311(2)=Cij134(22,8)+Cij234(19,7)
       S3h000000311(3)=Cij134(24,8)+Cij234(20,7)
       S3h000000312(1)=Cij134(22,8)+Cij234(19,7)
       S3h000000312(2)=Cij134(22,8)-Cij234(22,8)
       S3h000000312(3)=Cij134(24,8)-Cij234(24,8)
       S3h000000313(1)=Cij134(24,8)+Cij234(20,7)
       S3h000000313(2)=Cij134(24,8)-Cij234(24,8)
       S3h000000313(3)=Cij134(23,8)-Cij234(23,8)
       S3h000000321(1)=Cij124(22,8)-Cij134(22,8)
       S3h000000321(2)=Cij124(24,8)-Cij134(22,8)
       S3h000000321(3)=Cij124(24,8)-Cij134(24,8)
       S3h000000322(1)=Cij124(24,8)-Cij134(22,8)
       S3h000000322(2)=Cij124(23,8)-Cij134(22,8)
       S3h000000322(3)=Cij124(23,8)-Cij134(24,8)
       S3h000000323(1)=Cij124(24,8)-Cij134(24,8)
       S3h000000323(2)=Cij124(23,8)-Cij134(24,8)
       S3h000000323(3)=Cij124(23,8)-Cij134(23,8)
       S3h000000331(1)=Cij123(22,8)-Cij124(22,8)
       S3h000000331(2)=Cij123(24,8)-Cij124(24,8)
       S3h000000331(3)=-Cij124(24,8)
       S3h000000332(1)=Cij123(24,8)-Cij124(24,8)
       S3h000000332(2)=Cij123(23,8)-Cij124(23,8)
       S3h000000332(3)=-Cij124(23,8)
       S3h000000333(1)=-Cij124(24,8)
       S3h000000333(2)=-Cij124(23,8)
       S3h000000333(3)=-Cij124(23,8)
       aux000000001(1)=-(F(1)*S3h00000021(1))-F(2)*S3h00000022(1)-F(3)*S
     &  3h00000023(1)+S3h000000311(k)*Z(1,l)+S3h000000321(k)*Z(2,l)+S3h0
     &  00000331(k)*Z(3,l)+(Inv801+S3000000001(1)-S3h000000311(1)-S3h000
     &  000322(1)-S3h000000333(1))*Z(k,l)-2*S3h000000001(1)*ZZ(k,1,l,1)-
     &  2*S3h000000001(2)*ZZ(k,1,l,2)-2*S3h000000001(3)*ZZ(k,1,l,3)
       aux000000001(2)=-(F(1)*S3h00000021(2))-F(2)*S3h00000022(2)-F(3)*S
     &  3h00000023(2)+S3h000000312(k)*Z(1,l)+S3h000000322(k)*Z(2,l)+S3h0
     &  00000332(k)*Z(3,l)+(Inv802+S3000000001(2)-S3h000000312(1)-S3h000
     &  000322(2)-S3h000000333(2))*Z(k,l)-2*S3h000000001(1)*ZZ(k,2,l,1)-
     &  2*S3h000000001(2)*ZZ(k,2,l,2)-2*S3h000000001(3)*ZZ(k,2,l,3)
       aux000000001(3)=-(F(1)*S3h00000021(3))-F(2)*S3h00000022(3)-F(3)*S
     &  3h00000023(3)+S3h000000313(k)*Z(1,l)+S3h000000323(k)*Z(2,l)+S3h0
     &  00000333(k)*Z(3,l)+(Inv803+S3000000001(3)-S3h000000313(1)-S3h000
     &  000323(2)-S3h000000333(3))*Z(k,l)-2*S3h000000001(1)*ZZ(k,3,l,1)-
     &  2*S3h000000001(2)*ZZ(k,3,l,2)-2*S3h000000001(3)*ZZ(k,3,l,3)
       temp000000001(1)=I20Z*(aux000000001(1)+2*tempD400000000*F(4)+F(5)
     &  *temp0000001(1))
       temp000000001(2)=I20Z*(aux000000001(2)+2*tempD400000000*F(6)+F(5)
     &  *temp0000001(2))
       temp000000001(3)=I20Z*(aux000000001(3)+2*tempD400000000*F(7)+F(5)
     &  *temp0000001(3))
       S3h000051111(1)=Cij134(17,8)-Cij234(9,4)
       S3h000051111(2)=Cij134(17,8)+Cij234(11,5)
       S3h000051111(3)=Cij134(19,8)+Cij234(12,5)
       S3h000051211(1)=Cij134(17,8)+Cij234(11,5)
       S3h000051211(2)=Cij134(17,8)-Cij234(13,6)
       S3h000051211(3)=Cij134(19,8)-Cij234(15,6)
       S3h000051221(1)=Cij134(17,8)-Cij234(13,6)
       S3h000051221(2)=Cij134(17,8)+Cij234(15,7)
       S3h000051221(3)=Cij134(19,8)+Cij234(17,7)
       S3h000051222(1)=Cij134(17,8)+Cij234(15,7)
       S3h000051222(2)=Cij134(17,8)-Cij234(17,8)
       S3h000051222(3)=Cij134(19,8)-Cij234(19,8)
       S3h000051311(1)=Cij134(19,8)+Cij234(12,5)
       S3h000051311(2)=Cij134(19,8)-Cij234(15,6)
       S3h000051311(3)=Cij134(20,8)-Cij234(14,6)
       S3h000051321(1)=Cij134(19,8)-Cij234(15,6)
       S3h000051321(2)=Cij134(19,8)+Cij234(17,7)
       S3h000051321(3)=Cij134(20,8)+Cij234(18,7)
       S3h000051322(1)=Cij134(19,8)+Cij234(17,7)
       S3h000051322(2)=Cij134(19,8)-Cij234(19,8)
       S3h000051322(3)=Cij134(20,8)-Cij234(20,8)
       S3h000051331(1)=Cij134(20,8)-Cij234(14,6)
       S3h000051331(2)=Cij134(20,8)+Cij234(18,7)
       S3h000051331(3)=Cij134(21,8)+Cij234(16,7)
       S3h000051332(1)=Cij134(20,8)+Cij234(18,7)
       S3h000051332(2)=Cij134(20,8)-Cij234(20,8)
       S3h000051332(3)=Cij134(21,8)-Cij234(21,8)
       S3h000051333(1)=Cij134(21,8)+Cij234(16,7)
       S3h000051333(2)=Cij134(21,8)-Cij234(21,8)
       S3h000051333(3)=Cij134(18,8)-Cij234(18,8)
       S3h000052111(1)=Cij124(17,8)-Cij134(17,8)
       S3h000052111(2)=Cij124(19,8)-Cij134(17,8)
       S3h000052111(3)=Cij124(19,8)-Cij134(19,8)
       S3h000052211(1)=Cij124(19,8)-Cij134(17,8)
       S3h000052211(2)=Cij124(20,8)-Cij134(17,8)
       S3h000052211(3)=Cij124(20,8)-Cij134(19,8)
       S3h000052221(1)=Cij124(20,8)-Cij134(17,8)
       S3h000052221(2)=Cij124(21,8)-Cij134(17,8)
       S3h000052221(3)=Cij124(21,8)-Cij134(19,8)
       S3h000052222(1)=Cij124(21,8)-Cij134(17,8)
       S3h000052222(2)=Cij124(18,8)-Cij134(17,8)
       S3h000052222(3)=Cij124(18,8)-Cij134(19,8)
       S3h000052311(1)=Cij124(19,8)-Cij134(19,8)
       S3h000052311(2)=Cij124(20,8)-Cij134(19,8)
       S3h000052311(3)=Cij124(20,8)-Cij134(20,8)
       S3h000052321(1)=Cij124(20,8)-Cij134(19,8)
       S3h000052321(2)=Cij124(21,8)-Cij134(19,8)
       S3h000052321(3)=Cij124(21,8)-Cij134(20,8)
       S3h000052322(1)=Cij124(21,8)-Cij134(19,8)
       S3h000052322(2)=Cij124(18,8)-Cij134(19,8)
       S3h000052322(3)=Cij124(18,8)-Cij134(20,8)
       S3h000052331(1)=Cij124(20,8)-Cij134(20,8)
       S3h000052331(2)=Cij124(21,8)-Cij134(20,8)
       S3h000052331(3)=Cij124(21,8)-Cij134(21,8)
       S3h000052332(1)=Cij124(21,8)-Cij134(20,8)
       S3h000052332(2)=Cij124(18,8)-Cij134(20,8)
       S3h000052332(3)=Cij124(18,8)-Cij134(21,8)
       S3h000052333(1)=Cij124(21,8)-Cij134(21,8)
       S3h000052333(2)=Cij124(18,8)-Cij134(21,8)
       S3h000052333(3)=Cij124(18,8)-Cij134(18,8)
       S3h000053111(1)=Cij123(17,8)-Cij124(17,8)
       S3h000053111(2)=Cij123(19,8)-Cij124(19,8)
       S3h000053111(3)=-Cij124(19,8)
       S3h000053211(1)=Cij123(19,8)-Cij124(19,8)
       S3h000053211(2)=Cij123(20,8)-Cij124(20,8)
       S3h000053211(3)=-Cij124(20,8)
       S3h000053221(1)=Cij123(20,8)-Cij124(20,8)
       S3h000053221(2)=Cij123(21,8)-Cij124(21,8)
       S3h000053221(3)=-Cij124(21,8)
       S3h000053222(1)=Cij123(21,8)-Cij124(21,8)
       S3h000053222(2)=Cij123(18,8)-Cij124(18,8)
       S3h000053222(3)=-Cij124(18,8)
       S3h000053311(1)=-Cij124(19,8)
       S3h000053311(2)=-Cij124(20,8)
       S3h000053311(3)=-Cij124(20,8)
       S3h000053321(1)=-Cij124(20,8)
       S3h000053321(2)=-Cij124(21,8)
       S3h000053321(3)=-Cij124(21,8)
       S3h000053322(1)=-Cij124(21,8)
       S3h000053322(2)=-Cij124(18,8)
       S3h000053322(3)=-Cij124(18,8)
       S3h000053331(1)=-Cij124(20,8)
       S3h000053331(2)=-Cij124(21,8)
       S3h000053331(3)=-Cij124(21,8)
       S3h000053332(1)=-Cij124(21,8)
       S3h000053332(2)=-Cij124(18,8)
       S3h000053332(3)=-Cij124(18,8)
       S3h000053333(1)=-Cij124(21,8)
       S3h000053333(2)=-Cij124(18,8)
       S3h000053333(3)=-Cij124(18,8)
       S3000000311(1)=-2*Cij234(9,4)
       S3000000311(2)=2*Cij234(11,5)
       S3000000311(3)=2*Cij234(12,5)
       S3000000312(1)=2*Cij234(11,5)
       S3000000312(2)=-2*Cij234(13,6)
       S3000000312(3)=-2*Cij234(15,6)
       S3000000313(1)=2*Cij234(12,5)
       S3000000313(2)=-2*Cij234(15,6)
       S3000000313(3)=-2*Cij234(14,6)
       S3000000321(1)=2*Cij234(11,5)
       S3000000321(2)=-2*Cij234(13,6)
       S3000000321(3)=-2*Cij234(15,6)
       S3000000322(1)=-2*Cij234(13,6)
       S3000000322(2)=2*Cij234(15,7)
       S3000000322(3)=2*Cij234(17,7)
       S3000000323(1)=-2*Cij234(15,6)
       S3000000323(2)=2*Cij234(17,7)
       S3000000323(3)=2*Cij234(18,7)
       S3000000331(1)=2*Cij234(12,5)
       S3000000331(2)=-2*Cij234(15,6)
       S3000000331(3)=-2*Cij234(14,6)
       S3000000332(1)=-2*Cij234(15,6)
       S3000000332(2)=2*Cij234(17,7)
       S3000000332(3)=2*Cij234(18,7)
       S3000000333(1)=-2*Cij234(14,6)
       S3000000333(2)=2*Cij234(18,7)
       S3000000333(3)=2*Cij234(16,7)
       aux0000003(1,1,1)=-(F(1)*S3h00004111(1))-F(2)*S3h00004211(1)-F(3)
     &  *S3h00004311(1)+S3h000051111(k)*Z(1,l)+S3h000052111(k)*Z(2,l)+S3
     &  h000053111(k)*Z(3,l)+(Inv60111+S3000000311(1)-S3h000051111(1)-S3
     &  h000052211(1)-S3h000053311(1))*Z(k,l)-6*S3h000000311(1)*ZZ(k,1,l
     &  ,1)-6*S3h000000321(1)*ZZ(k,1,l,2)-6*S3h000000331(1)*ZZ(k,1,l,3)
       aux0000003(2,1,1)=-(F(1)*S3h00004121(1))-F(2)*S3h00004221(1)-F(3)
     &  *S3h00004321(1)+S3h000051211(k)*Z(1,l)+S3h000052211(k)*Z(2,l)+S3
     &  h000053211(k)*Z(3,l)+(Inv60211+S3000000321(1)-S3h000051211(1)-S3
     &  h000052221(1)-S3h000053321(1))*Z(k,l)-4*S3h000000312(1)*ZZ(k,1,l
     &  ,1)-4*S3h000000322(1)*ZZ(k,1,l,2)-4*S3h000000332(1)*ZZ(k,1,l,3)-
     &  2*S3h000000311(1)*ZZ(k,2,l,1)-2*S3h000000321(1)*ZZ(k,2,l,2)-2*S3
     &  h000000331(1)*ZZ(k,2,l,3)
       aux0000003(2,2,1)=-(F(1)*S3h00004122(1))-F(2)*S3h00004222(1)-F(3)
     &  *S3h00004322(1)+S3h000051221(k)*Z(1,l)+S3h000052221(k)*Z(2,l)+S3
     &  h000053221(k)*Z(3,l)+(Inv60221+S3000000322(1)-S3h000051221(1)-S3
     &  h000052222(1)-S3h000053322(1))*Z(k,l)-2*S3h000000312(2)*ZZ(k,1,l
     &  ,1)-2*S3h000000322(2)*ZZ(k,1,l,2)-2*S3h000000332(2)*ZZ(k,1,l,3)-
     &  4*S3h000000312(1)*ZZ(k,2,l,1)-4*S3h000000322(1)*ZZ(k,2,l,2)-4*S3
     &  h000000332(1)*ZZ(k,2,l,3)
       aux0000003(2,2,2)=-(F(1)*S3h00004122(2))-F(2)*S3h00004222(2)-F(3)
     &  *S3h00004322(2)+S3h000051222(k)*Z(1,l)+S3h000052222(k)*Z(2,l)+S3
     &  h000053222(k)*Z(3,l)+(Inv60222+S3000000322(2)-S3h000051222(1)-S3
     &  h000052222(2)-S3h000053322(2))*Z(k,l)-6*S3h000000312(2)*ZZ(k,2,l
     &  ,1)-6*S3h000000322(2)*ZZ(k,2,l,2)-6*S3h000000332(2)*ZZ(k,2,l,3)
       aux0000003(3,1,1)=-(F(1)*S3h00004131(1))-F(2)*S3h00004231(1)-F(3)
     &  *S3h00004331(1)+S3h000051311(k)*Z(1,l)+S3h000052311(k)*Z(2,l)+S3
     &  h000053311(k)*Z(3,l)+(Inv60311+S3000000331(1)-S3h000051311(1)-S3
     &  h000052321(1)-S3h000053331(1))*Z(k,l)-4*S3h000000313(1)*ZZ(k,1,l
     &  ,1)-4*S3h000000323(1)*ZZ(k,1,l,2)-4*S3h000000333(1)*ZZ(k,1,l,3)-
     &  2*S3h000000311(1)*ZZ(k,3,l,1)-2*S3h000000321(1)*ZZ(k,3,l,2)-2*S3
     &  h000000331(1)*ZZ(k,3,l,3)
       aux0000003(3,2,1)=-(F(1)*S3h00004132(1))-F(2)*S3h00004232(1)-F(3)
     &  *S3h00004332(1)+S3h000051321(k)*Z(1,l)+S3h000052321(k)*Z(2,l)+S3
     &  h000053321(k)*Z(3,l)+(Inv60321+S3000000332(1)-S3h000051321(1)-S3
     &  h000052322(1)-S3h000053332(1))*Z(k,l)-2*S3h000000313(2)*ZZ(k,1,l
     &  ,1)-2*S3h000000323(2)*ZZ(k,1,l,2)-2*S3h000000333(2)*ZZ(k,1,l,3)-
     &  2*S3h000000313(1)*ZZ(k,2,l,1)-2*S3h000000323(1)*ZZ(k,2,l,2)-2*S3
     &  h000000333(1)*ZZ(k,2,l,3)-2*S3h000000312(1)*ZZ(k,3,l,1)-2*S3h000
     &  000322(1)*ZZ(k,3,l,2)-2*S3h000000332(1)*ZZ(k,3,l,3)
       aux0000003(3,2,2)=-(F(1)*S3h00004132(2))-F(2)*S3h00004232(2)-F(3)
     &  *S3h00004332(2)+S3h000051322(k)*Z(1,l)+S3h000052322(k)*Z(2,l)+S3
     &  h000053322(k)*Z(3,l)+(Inv60322+S3000000332(2)-S3h000051322(1)-S3
     &  h000052322(2)-S3h000053332(2))*Z(k,l)-4*S3h000000313(2)*ZZ(k,2,l
     &  ,1)-4*S3h000000323(2)*ZZ(k,2,l,2)-4*S3h000000333(2)*ZZ(k,2,l,3)-
     &  2*S3h000000312(2)*ZZ(k,3,l,1)-2*S3h000000322(2)*ZZ(k,3,l,2)-2*S3
     &  h000000332(2)*ZZ(k,3,l,3)
       aux0000003(3,3,1)=-(F(1)*S3h00004133(1))-F(2)*S3h00004233(1)-F(3)
     &  *S3h00004333(1)+S3h000051331(k)*Z(1,l)+S3h000052331(k)*Z(2,l)+S3
     &  h000053331(k)*Z(3,l)+(Inv60331+S3000000333(1)-S3h000051331(1)-S3
     &  h000052332(1)-S3h000053333(1))*Z(k,l)-2*S3h000000313(3)*ZZ(k,1,l
     &  ,1)-2*S3h000000323(3)*ZZ(k,1,l,2)-2*S3h000000333(3)*ZZ(k,1,l,3)-
     &  4*S3h000000313(1)*ZZ(k,3,l,1)-4*S3h000000323(1)*ZZ(k,3,l,2)-4*S3
     &  h000000333(1)*ZZ(k,3,l,3)
       aux0000003(3,3,2)=-(F(1)*S3h00004133(2))-F(2)*S3h00004233(2)-F(3)
     &  *S3h00004333(2)+S3h000051332(k)*Z(1,l)+S3h000052332(k)*Z(2,l)+S3
     &  h000053332(k)*Z(3,l)+(Inv60332+S3000000333(2)-S3h000051332(1)-S3
     &  h000052332(2)-S3h000053333(2))*Z(k,l)-2*S3h000000313(3)*ZZ(k,2,l
     &  ,1)-2*S3h000000323(3)*ZZ(k,2,l,2)-2*S3h000000333(3)*ZZ(k,2,l,3)-
     &  4*S3h000000313(2)*ZZ(k,3,l,1)-4*S3h000000323(2)*ZZ(k,3,l,2)-4*S3
     &  h000000333(2)*ZZ(k,3,l,3)
       aux0000003(3,3,3)=-(F(1)*S3h00004133(3))-F(2)*S3h00004233(3)-F(3)
     &  *S3h00004333(3)+S3h000051333(k)*Z(1,l)+S3h000052333(k)*Z(2,l)+S3
     &  h000053333(k)*Z(3,l)+(Inv60333+S3000000333(3)-S3h000051333(1)-S3
     &  h000052333(2)-S3h000053333(3))*Z(k,l)-6*S3h000000313(3)*ZZ(k,3,l
     &  ,1)-6*S3h000000323(3)*ZZ(k,3,l,2)-6*S3h000000333(3)*ZZ(k,3,l,3)
       temp0000003(1,1,1)=I24Z*(aux0000003(1,1,1)+6*F(4)*temp0000002(1,1
     &  )+F(5)*temp00003(1,1,1)+24*temp000000001(1)*ZZ(k,1,l,1))
       temp0000003(2,1,1)=I24Z*(aux0000003(2,1,1)+2*F(6)*temp0000002(1,1
     &  )+4*F(4)*temp0000002(2,1)+F(5)*temp00003(2,1,1)+8*(temp000000001
     &  (2)*ZZ(k,1,l,1)+temp000000001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))))
       temp0000003(2,2,1)=I24Z*(aux0000003(2,2,1)+4*F(6)*temp0000002(2,1
     &  )+2*F(4)*temp0000002(2,2)+F(5)*temp00003(2,2,1)+8*(temp000000001
     &  (2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp000000001(1)*ZZ(k,2,l,2)))
       temp0000003(2,2,2)=I24Z*(aux0000003(2,2,2)+6*F(6)*temp0000002(2,2
     &  )+F(5)*temp00003(2,2,2)+24*temp000000001(2)*ZZ(k,2,l,2))
       temp0000003(3,1,1)=I24Z*(aux0000003(3,1,1)+2*F(7)*temp0000002(1,1
     &  )+4*F(4)*temp0000002(3,1)+F(5)*temp00003(3,1,1)+8*(temp000000001
     &  (3)*ZZ(k,1,l,1)+temp000000001(1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))))
       temp0000003(3,2,1)=I24Z*(aux0000003(3,2,1)+2*F(7)*temp0000002(2,1
     &  )+2*F(6)*temp0000002(3,1)+2*F(4)*temp0000002(3,2)+F(5)*temp00003
     &  (3,2,1)+4*(temp000000001(3)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp000000
     &  001(2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)))+4*temp000000001(1)*(ZZ(k,2,l,3
     &  )+ZZ(k,3,l,2)))
       temp0000003(3,2,2)=I24Z*(aux0000003(3,2,2)+2*F(7)*temp0000002(2,2
     &  )+4*F(6)*temp0000002(3,2)+F(5)*temp00003(3,2,2)+8*(temp000000001
     &  (3)*ZZ(k,2,l,2)+temp000000001(2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))))
       temp0000003(3,3,1)=I24Z*(aux0000003(3,3,1)+4*F(7)*temp0000002(3,1
     &  )+2*F(4)*temp0000002(3,3)+F(5)*temp00003(3,3,1)+8*(temp000000001
     &  (3)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+temp000000001(1)*ZZ(k,3,l,3)))
       temp0000003(3,3,2)=I24Z*(aux0000003(3,3,2)+4*F(7)*temp0000002(3,2
     &  )+2*F(6)*temp0000002(3,3)+F(5)*temp00003(3,3,2)+8*(temp000000001
     &  (3)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+temp000000001(2)*ZZ(k,3,l,3)))
       temp0000003(3,3,3)=I24Z*(aux0000003(3,3,3)+6*F(7)*temp0000002(3,3
     &  )+F(5)*temp00003(3,3,3)+24*temp000000001(3)*ZZ(k,3,l,3))
       temp0000003(1,1,2)=temp0000003(2,1,1)
       temp0000003(1,1,3)=temp0000003(3,1,1)
       temp0000003(1,2,1)=temp0000003(2,1,1)
       temp0000003(1,2,2)=temp0000003(2,2,1)
       temp0000003(1,2,3)=temp0000003(3,2,1)
       temp0000003(1,3,1)=temp0000003(3,1,1)
       temp0000003(1,3,2)=temp0000003(3,2,1)
       temp0000003(1,3,3)=temp0000003(3,3,1)
       temp0000003(2,1,2)=temp0000003(2,2,1)
       temp0000003(2,1,3)=temp0000003(3,2,1)
       temp0000003(2,2,3)=temp0000003(3,2,2)
       temp0000003(2,3,1)=temp0000003(3,2,1)
       temp0000003(2,3,2)=temp0000003(3,2,2)
       temp0000003(2,3,3)=temp0000003(3,3,2)
       temp0000003(3,1,2)=temp0000003(3,2,1)
       temp0000003(3,1,3)=temp0000003(3,3,1)
       temp0000003(3,2,3)=temp0000003(3,3,2)
       S3000051111(1)=-2*Cij234(4,2)
       S3000051111(2)=2*Cij234(5,3)
       S3000051111(3)=2*Cij234(6,3)
       S3000051211(1)=2*Cij234(5,3)
       S3000051211(2)=-2*Cij234(6,4)
       S3000051211(3)=-2*Cij234(8,4)
       S3000051221(1)=-2*Cij234(6,4)
       S3000051221(2)=2*Cij234(7,5)
       S3000051221(3)=2*Cij234(9,5)
       S3000051222(1)=2*Cij234(7,5)
       S3000051222(2)=-2*Cij234(8,6)
       S3000051222(3)=-2*Cij234(10,6)
       S3000051311(1)=2*Cij234(6,3)
       S3000051311(2)=-2*Cij234(8,4)
       S3000051311(3)=-2*Cij234(7,4)
       S3000051321(1)=-2*Cij234(8,4)
       S3000051321(2)=2*Cij234(9,5)
       S3000051321(3)=2*Cij234(10,5)
       S3000051322(1)=2*Cij234(9,5)
       S3000051322(2)=-2*Cij234(10,6)
       S3000051322(3)=-2*Cij234(11,6)
       S3000051331(1)=-2*Cij234(7,4)
       S3000051331(2)=2*Cij234(10,5)
       S3000051331(3)=2*Cij234(8,5)
       S3000051332(1)=2*Cij234(10,5)
       S3000051332(2)=-2*Cij234(11,6)
       S3000051332(3)=-2*Cij234(12,6)
       S3000051333(1)=2*Cij234(8,5)
       S3000051333(2)=-2*Cij234(12,6)
       S3000051333(3)=-2*Cij234(9,6)
       S3000052111(1)=2*Cij234(5,3)
       S3000052111(2)=-2*Cij234(6,4)
       S3000052111(3)=-2*Cij234(8,4)
       S3000052211(1)=-2*Cij234(6,4)
       S3000052211(2)=2*Cij234(7,5)
       S3000052211(3)=2*Cij234(9,5)
       S3000052221(1)=2*Cij234(7,5)
       S3000052221(2)=-2*Cij234(8,6)
       S3000052221(3)=-2*Cij234(10,6)
       S3000052222(1)=-2*Cij234(8,6)
       S3000052222(2)=2*Cij234(9,7)
       S3000052222(3)=2*Cij234(11,7)
       S3000052311(1)=-2*Cij234(8,4)
       S3000052311(2)=2*Cij234(9,5)
       S3000052311(3)=2*Cij234(10,5)
       S3000052321(1)=2*Cij234(9,5)
       S3000052321(2)=-2*Cij234(10,6)
       S3000052321(3)=-2*Cij234(11,6)
       S3000052322(1)=-2*Cij234(10,6)
       S3000052322(2)=2*Cij234(11,7)
       S3000052322(3)=2*Cij234(12,7)
       S3000052331(1)=2*Cij234(10,5)
       S3000052331(2)=-2*Cij234(11,6)
       S3000052331(3)=-2*Cij234(12,6)
       S3000052332(1)=-2*Cij234(11,6)
       S3000052332(2)=2*Cij234(12,7)
       S3000052332(3)=2*Cij234(13,7)
       S3000052333(1)=-2*Cij234(12,6)
       S3000052333(2)=2*Cij234(13,7)
       S3000052333(3)=2*Cij234(14,7)
       S3000053111(1)=2*Cij234(6,3)
       S3000053111(2)=-2*Cij234(8,4)
       S3000053111(3)=-2*Cij234(7,4)
       S3000053211(1)=-2*Cij234(8,4)
       S3000053211(2)=2*Cij234(9,5)
       S3000053211(3)=2*Cij234(10,5)
       S3000053221(1)=2*Cij234(9,5)
       S3000053221(2)=-2*Cij234(10,6)
       S3000053221(3)=-2*Cij234(11,6)
       S3000053222(1)=-2*Cij234(10,6)
       S3000053222(2)=2*Cij234(11,7)
       S3000053222(3)=2*Cij234(12,7)
       S3000053311(1)=-2*Cij234(7,4)
       S3000053311(2)=2*Cij234(10,5)
       S3000053311(3)=2*Cij234(8,5)
       S3000053321(1)=2*Cij234(10,5)
       S3000053321(2)=-2*Cij234(11,6)
       S3000053321(3)=-2*Cij234(12,6)
       S3000053322(1)=-2*Cij234(11,6)
       S3000053322(2)=2*Cij234(12,7)
       S3000053322(3)=2*Cij234(13,7)
       S3000053331(1)=2*Cij234(8,5)
       S3000053331(2)=-2*Cij234(12,6)
       S3000053331(3)=-2*Cij234(9,6)
       S3000053332(1)=-2*Cij234(12,6)
       S3000053332(2)=2*Cij234(13,7)
       S3000053332(3)=2*Cij234(14,7)
       S3000053333(1)=-2*Cij234(9,6)
       S3000053333(2)=2*Cij234(14,7)
       S3000053333(3)=2*Cij234(10,7)
       S3007111111(1)=-2*C0234
       S3007111111(2)=2*Cij234(1,1)
       S3007111111(3)=2*Cij234(2,1)
       S3007121111(1)=2*Cij234(1,1)
       S3007121111(2)=-2*Cij234(1,2)
       S3007121111(3)=-2*Cij234(3,2)
       S3007122111(1)=-2*Cij234(1,2)
       S3007122111(2)=2*Cij234(1,3)
       S3007122111(3)=2*Cij234(3,3)
       S3007122211(1)=2*Cij234(1,3)
       S3007122211(2)=-2*Cij234(1,4)
       S3007122211(3)=-2*Cij234(3,4)
       S3007122221(1)=-2*Cij234(1,4)
       S3007122221(2)=2*Cij234(1,5)
       S3007122221(3)=2*Cij234(3,5)
       S3007122222(1)=2*Cij234(1,5)
       S3007122222(2)=-2*Cij234(1,6)
       S3007122222(3)=-2*Cij234(3,6)
       S3007131111(1)=2*Cij234(2,1)
       S3007131111(2)=-2*Cij234(3,2)
       S3007131111(3)=-2*Cij234(2,2)
       S3007132111(1)=-2*Cij234(3,2)
       S3007132111(2)=2*Cij234(3,3)
       S3007132111(3)=2*Cij234(4,3)
       S3007132211(1)=2*Cij234(3,3)
       S3007132211(2)=-2*Cij234(3,4)
       S3007132211(3)=-2*Cij234(4,4)
       S3007132221(1)=-2*Cij234(3,4)
       S3007132221(2)=2*Cij234(3,5)
       S3007132221(3)=2*Cij234(4,5)
       S3007132222(1)=2*Cij234(3,5)
       S3007132222(2)=-2*Cij234(3,6)
       S3007132222(3)=-2*Cij234(4,6)
       S3007133111(1)=-2*Cij234(2,2)
       S3007133111(2)=2*Cij234(4,3)
       S3007133111(3)=2*Cij234(2,3)
       S3007133211(1)=2*Cij234(4,3)
       S3007133211(2)=-2*Cij234(4,4)
       S3007133211(3)=-2*Cij234(5,4)
       S3007133221(1)=-2*Cij234(4,4)
       S3007133221(2)=2*Cij234(4,5)
       S3007133221(3)=2*Cij234(5,5)
       S3007133222(1)=2*Cij234(4,5)
       S3007133222(2)=-2*Cij234(4,6)
       S3007133222(3)=-2*Cij234(5,6)
       S3007133311(1)=2*Cij234(2,3)
       S3007133311(2)=-2*Cij234(5,4)
       S3007133311(3)=-2*Cij234(2,4)
       S3007133321(1)=-2*Cij234(5,4)
       S3007133321(2)=2*Cij234(5,5)
       S3007133321(3)=2*Cij234(6,5)
       S3007133322(1)=2*Cij234(5,5)
       S3007133322(2)=-2*Cij234(5,6)
       S3007133322(3)=-2*Cij234(6,6)
       S3007133331(1)=-2*Cij234(2,4)
       S3007133331(2)=2*Cij234(6,5)
       S3007133331(3)=2*Cij234(2,5)
       S3007133332(1)=2*Cij234(6,5)
       S3007133332(2)=-2*Cij234(6,6)
       S3007133332(3)=-2*Cij234(7,6)
       S3007133333(1)=2*Cij234(2,5)
       S3007133333(2)=-2*Cij234(7,6)
       S3007133333(3)=-2*Cij234(2,6)
       S3007211111(1)=2*Cij234(1,1)
       S3007211111(2)=-2*Cij234(1,2)
       S3007211111(3)=-2*Cij234(3,2)
       S3007221111(1)=-2*Cij234(1,2)
       S3007221111(2)=2*Cij234(1,3)
       S3007221111(3)=2*Cij234(3,3)
       S3007222111(1)=2*Cij234(1,3)
       S3007222111(2)=-2*Cij234(1,4)
       S3007222111(3)=-2*Cij234(3,4)
       S3007222211(1)=-2*Cij234(1,4)
       S3007222211(2)=2*Cij234(1,5)
       S3007222211(3)=2*Cij234(3,5)
       S3007222221(1)=2*Cij234(1,5)
       S3007222221(2)=-2*Cij234(1,6)
       S3007222221(3)=-2*Cij234(3,6)
       S3007222222(1)=-2*Cij234(1,6)
       S3007222222(2)=2*Cij234(1,7)
       S3007222222(3)=2*Cij234(3,7)
       S3007231111(1)=-2*Cij234(3,2)
       S3007231111(2)=2*Cij234(3,3)
       S3007231111(3)=2*Cij234(4,3)
       S3007232111(1)=2*Cij234(3,3)
       S3007232111(2)=-2*Cij234(3,4)
       S3007232111(3)=-2*Cij234(4,4)
       S3007232211(1)=-2*Cij234(3,4)
       S3007232211(2)=2*Cij234(3,5)
       S3007232211(3)=2*Cij234(4,5)
       S3007232221(1)=2*Cij234(3,5)
       S3007232221(2)=-2*Cij234(3,6)
       S3007232221(3)=-2*Cij234(4,6)
       S3007232222(1)=-2*Cij234(3,6)
       S3007232222(2)=2*Cij234(3,7)
       S3007232222(3)=2*Cij234(4,7)
       S3007233111(1)=2*Cij234(4,3)
       S3007233111(2)=-2*Cij234(4,4)
       S3007233111(3)=-2*Cij234(5,4)
       S3007233211(1)=-2*Cij234(4,4)
       S3007233211(2)=2*Cij234(4,5)
       S3007233211(3)=2*Cij234(5,5)
       S3007233221(1)=2*Cij234(4,5)
       S3007233221(2)=-2*Cij234(4,6)
       S3007233221(3)=-2*Cij234(5,6)
       S3007233222(1)=-2*Cij234(4,6)
       S3007233222(2)=2*Cij234(4,7)
       S3007233222(3)=2*Cij234(5,7)
       S3007233311(1)=-2*Cij234(5,4)
       S3007233311(2)=2*Cij234(5,5)
       S3007233311(3)=2*Cij234(6,5)
       S3007233321(1)=2*Cij234(5,5)
       S3007233321(2)=-2*Cij234(5,6)
       S3007233321(3)=-2*Cij234(6,6)
       S3007233322(1)=-2*Cij234(5,6)
       S3007233322(2)=2*Cij234(5,7)
       S3007233322(3)=2*Cij234(6,7)
       S3007233331(1)=2*Cij234(6,5)
       S3007233331(2)=-2*Cij234(6,6)
       S3007233331(3)=-2*Cij234(7,6)
       S3007233332(1)=-2*Cij234(6,6)
       S3007233332(2)=2*Cij234(6,7)
       S3007233332(3)=2*Cij234(7,7)
       S3007233333(1)=-2*Cij234(7,6)
       S3007233333(2)=2*Cij234(7,7)
       S3007233333(3)=2*Cij234(8,7)
       S3007311111(1)=2*Cij234(2,1)
       S3007311111(2)=-2*Cij234(3,2)
       S3007311111(3)=-2*Cij234(2,2)
       S3007321111(1)=-2*Cij234(3,2)
       S3007321111(2)=2*Cij234(3,3)
       S3007321111(3)=2*Cij234(4,3)
       S3007322111(1)=2*Cij234(3,3)
       S3007322111(2)=-2*Cij234(3,4)
       S3007322111(3)=-2*Cij234(4,4)
       S3007322211(1)=-2*Cij234(3,4)
       S3007322211(2)=2*Cij234(3,5)
       S3007322211(3)=2*Cij234(4,5)
       S3007322221(1)=2*Cij234(3,5)
       S3007322221(2)=-2*Cij234(3,6)
       S3007322221(3)=-2*Cij234(4,6)
       S3007322222(1)=-2*Cij234(3,6)
       S3007322222(2)=2*Cij234(3,7)
       S3007322222(3)=2*Cij234(4,7)
       S3007331111(1)=-2*Cij234(2,2)
       S3007331111(2)=2*Cij234(4,3)
       S3007331111(3)=2*Cij234(2,3)
       S3007332111(1)=2*Cij234(4,3)
       S3007332111(2)=-2*Cij234(4,4)
       S3007332111(3)=-2*Cij234(5,4)
       S3007332211(1)=-2*Cij234(4,4)
       S3007332211(2)=2*Cij234(4,5)
       S3007332211(3)=2*Cij234(5,5)
       S3007332221(1)=2*Cij234(4,5)
       S3007332221(2)=-2*Cij234(4,6)
       S3007332221(3)=-2*Cij234(5,6)
       S3007332222(1)=-2*Cij234(4,6)
       S3007332222(2)=2*Cij234(4,7)
       S3007332222(3)=2*Cij234(5,7)
       S3007333111(1)=2*Cij234(2,3)
       S3007333111(2)=-2*Cij234(5,4)
       S3007333111(3)=-2*Cij234(2,4)
       S3007333211(1)=-2*Cij234(5,4)
       S3007333211(2)=2*Cij234(5,5)
       S3007333211(3)=2*Cij234(6,5)
       S3007333221(1)=2*Cij234(5,5)
       S3007333221(2)=-2*Cij234(5,6)
       S3007333221(3)=-2*Cij234(6,6)
       S3007333222(1)=-2*Cij234(5,6)
       S3007333222(2)=2*Cij234(5,7)
       S3007333222(3)=2*Cij234(6,7)
       S3007333311(1)=-2*Cij234(2,4)
       S3007333311(2)=2*Cij234(6,5)
       S3007333311(3)=2*Cij234(2,5)
       S3007333321(1)=2*Cij234(6,5)
       S3007333321(2)=-2*Cij234(6,6)
       S3007333321(3)=-2*Cij234(7,6)
       S3007333322(1)=-2*Cij234(6,6)
       S3007333322(2)=2*Cij234(6,7)
       S3007333322(3)=2*Cij234(7,7)
       S3007333331(1)=2*Cij234(2,5)
       S3007333331(2)=-2*Cij234(7,6)
       S3007333331(3)=-2*Cij234(2,6)
       S3007333332(1)=-2*Cij234(7,6)
       S3007333332(2)=2*Cij234(7,7)
       S3007333332(3)=2*Cij234(8,7)
       S3007333333(1)=-2*Cij234(2,6)
       S3007333333(2)=2*Cij234(8,7)
       S3007333333(3)=2*Cij234(2,7)
       S3h007111111(1)=Cij134(10,8)-Cij234(4,2)
       S3h007111111(2)=Cij134(10,8)+Cij234(5,3)
       S3h007111111(3)=Cij134(12,8)+Cij234(6,3)
       S3h007121111(1)=Cij134(10,8)+Cij234(5,3)
       S3h007121111(2)=Cij134(10,8)-Cij234(6,4)
       S3h007121111(3)=Cij134(12,8)-Cij234(8,4)
       S3h007122111(1)=Cij134(10,8)-Cij234(6,4)
       S3h007122111(2)=Cij134(10,8)+Cij234(7,5)
       S3h007122111(3)=Cij134(12,8)+Cij234(9,5)
       S3h007122211(1)=Cij134(10,8)+Cij234(7,5)
       S3h007122211(2)=Cij134(10,8)-Cij234(8,6)
       S3h007122211(3)=Cij134(12,8)-Cij234(10,6)
       S3h007122221(1)=Cij134(10,8)-Cij234(8,6)
       S3h007122221(2)=Cij134(10,8)+Cij234(9,7)
       S3h007122221(3)=Cij134(12,8)+Cij234(11,7)
       S3h007122222(1)=Cij134(10,8)+Cij234(9,7)
       S3h007122222(2)=Cij134(10,8)-Cij234(10,8)
       S3h007122222(3)=Cij134(12,8)-Cij234(12,8)
       S3h007131111(1)=Cij134(12,8)+Cij234(6,3)
       S3h007131111(2)=Cij134(12,8)-Cij234(8,4)
       S3h007131111(3)=Cij134(13,8)-Cij234(7,4)
       S3h007132111(1)=Cij134(12,8)-Cij234(8,4)
       S3h007132111(2)=Cij134(12,8)+Cij234(9,5)
       S3h007132111(3)=Cij134(13,8)+Cij234(10,5)
       S3h007132211(1)=Cij134(12,8)+Cij234(9,5)
       S3h007132211(2)=Cij134(12,8)-Cij234(10,6)
       S3h007132211(3)=Cij134(13,8)-Cij234(11,6)
       S3h007132221(1)=Cij134(12,8)-Cij234(10,6)
       S3h007132221(2)=Cij134(12,8)+Cij234(11,7)
       S3h007132221(3)=Cij134(13,8)+Cij234(12,7)
       S3h007132222(1)=Cij134(12,8)+Cij234(11,7)
       S3h007132222(2)=Cij134(12,8)-Cij234(12,8)
       S3h007132222(3)=Cij134(13,8)-Cij234(13,8)
       S3h007133111(1)=Cij134(13,8)-Cij234(7,4)
       S3h007133111(2)=Cij134(13,8)+Cij234(10,5)
       S3h007133111(3)=Cij134(14,8)+Cij234(8,5)
       S3h007133211(1)=Cij134(13,8)+Cij234(10,5)
       S3h007133211(2)=Cij134(13,8)-Cij234(11,6)
       S3h007133211(3)=Cij134(14,8)-Cij234(12,6)
       S3h007133221(1)=Cij134(13,8)-Cij234(11,6)
       S3h007133221(2)=Cij134(13,8)+Cij234(12,7)
       S3h007133221(3)=Cij134(14,8)+Cij234(13,7)
       S3h007133222(1)=Cij134(13,8)+Cij234(12,7)
       S3h007133222(2)=Cij134(13,8)-Cij234(13,8)
       S3h007133222(3)=Cij134(14,8)-Cij234(14,8)
       S3h007133311(1)=Cij134(14,8)+Cij234(8,5)
       S3h007133311(2)=Cij134(14,8)-Cij234(12,6)
       S3h007133311(3)=Cij134(15,8)-Cij234(9,6)
       S3h007133321(1)=Cij134(14,8)-Cij234(12,6)
       S3h007133321(2)=Cij134(14,8)+Cij234(13,7)
       S3h007133321(3)=Cij134(15,8)+Cij234(14,7)
       S3h007133322(1)=Cij134(14,8)+Cij234(13,7)
       S3h007133322(2)=Cij134(14,8)-Cij234(14,8)
       S3h007133322(3)=Cij134(15,8)-Cij234(15,8)
       S3h007133331(1)=Cij134(15,8)-Cij234(9,6)
       S3h007133331(2)=Cij134(15,8)+Cij234(14,7)
       S3h007133331(3)=Cij134(16,8)+Cij234(10,7)
       S3h007133332(1)=Cij134(15,8)+Cij234(14,7)
       S3h007133332(2)=Cij134(15,8)-Cij234(15,8)
       S3h007133332(3)=Cij134(16,8)-Cij234(16,8)
       S3h007133333(1)=Cij134(16,8)+Cij234(10,7)
       S3h007133333(2)=Cij134(16,8)-Cij234(16,8)
       S3h007133333(3)=Cij134(11,8)-Cij234(11,8)
       S3h007211111(1)=Cij124(10,8)-Cij134(10,8)
       S3h007211111(2)=Cij124(12,8)-Cij134(10,8)
       S3h007211111(3)=Cij124(12,8)-Cij134(12,8)
       S3h007221111(1)=Cij124(12,8)-Cij134(10,8)
       S3h007221111(2)=Cij124(13,8)-Cij134(10,8)
       S3h007221111(3)=Cij124(13,8)-Cij134(12,8)
       S3h007222111(1)=Cij124(13,8)-Cij134(10,8)
       S3h007222111(2)=Cij124(14,8)-Cij134(10,8)
       S3h007222111(3)=Cij124(14,8)-Cij134(12,8)
       S3h007222211(1)=Cij124(14,8)-Cij134(10,8)
       S3h007222211(2)=Cij124(15,8)-Cij134(10,8)
       S3h007222211(3)=Cij124(15,8)-Cij134(12,8)
       S3h007222221(1)=Cij124(15,8)-Cij134(10,8)
       S3h007222221(2)=Cij124(16,8)-Cij134(10,8)
       S3h007222221(3)=Cij124(16,8)-Cij134(12,8)
       S3h007222222(1)=Cij124(16,8)-Cij134(10,8)
       S3h007222222(2)=Cij124(11,8)-Cij134(10,8)
       S3h007222222(3)=Cij124(11,8)-Cij134(12,8)
       S3h007231111(1)=Cij124(12,8)-Cij134(12,8)
       S3h007231111(2)=Cij124(13,8)-Cij134(12,8)
       S3h007231111(3)=Cij124(13,8)-Cij134(13,8)
       S3h007232111(1)=Cij124(13,8)-Cij134(12,8)
       S3h007232111(2)=Cij124(14,8)-Cij134(12,8)
       S3h007232111(3)=Cij124(14,8)-Cij134(13,8)
       S3h007232211(1)=Cij124(14,8)-Cij134(12,8)
       S3h007232211(2)=Cij124(15,8)-Cij134(12,8)
       S3h007232211(3)=Cij124(15,8)-Cij134(13,8)
       S3h007232221(1)=Cij124(15,8)-Cij134(12,8)
       S3h007232221(2)=Cij124(16,8)-Cij134(12,8)
       S3h007232221(3)=Cij124(16,8)-Cij134(13,8)
       S3h007232222(1)=Cij124(16,8)-Cij134(12,8)
       S3h007232222(2)=Cij124(11,8)-Cij134(12,8)
       S3h007232222(3)=Cij124(11,8)-Cij134(13,8)
       S3h007233111(1)=Cij124(13,8)-Cij134(13,8)
       S3h007233111(2)=Cij124(14,8)-Cij134(13,8)
       S3h007233111(3)=Cij124(14,8)-Cij134(14,8)
       S3h007233211(1)=Cij124(14,8)-Cij134(13,8)
       S3h007233211(2)=Cij124(15,8)-Cij134(13,8)
       S3h007233211(3)=Cij124(15,8)-Cij134(14,8)
       S3h007233221(1)=Cij124(15,8)-Cij134(13,8)
       S3h007233221(2)=Cij124(16,8)-Cij134(13,8)
       S3h007233221(3)=Cij124(16,8)-Cij134(14,8)
       S3h007233222(1)=Cij124(16,8)-Cij134(13,8)
       S3h007233222(2)=Cij124(11,8)-Cij134(13,8)
       S3h007233222(3)=Cij124(11,8)-Cij134(14,8)
       S3h007233311(1)=Cij124(14,8)-Cij134(14,8)
       S3h007233311(2)=Cij124(15,8)-Cij134(14,8)
       S3h007233311(3)=Cij124(15,8)-Cij134(15,8)
       S3h007233321(1)=Cij124(15,8)-Cij134(14,8)
       S3h007233321(2)=Cij124(16,8)-Cij134(14,8)
       S3h007233321(3)=Cij124(16,8)-Cij134(15,8)
       S3h007233322(1)=Cij124(16,8)-Cij134(14,8)
       S3h007233322(2)=Cij124(11,8)-Cij134(14,8)
       S3h007233322(3)=Cij124(11,8)-Cij134(15,8)
       S3h007233331(1)=Cij124(15,8)-Cij134(15,8)
       S3h007233331(2)=Cij124(16,8)-Cij134(15,8)
       S3h007233331(3)=Cij124(16,8)-Cij134(16,8)
       S3h007233332(1)=Cij124(16,8)-Cij134(15,8)
       S3h007233332(2)=Cij124(11,8)-Cij134(15,8)
       S3h007233332(3)=Cij124(11,8)-Cij134(16,8)
       S3h007233333(1)=Cij124(16,8)-Cij134(16,8)
       S3h007233333(2)=Cij124(11,8)-Cij134(16,8)
       S3h007233333(3)=Cij124(11,8)-Cij134(11,8)
       S3h007311111(1)=Cij123(10,8)-Cij124(10,8)
       S3h007311111(2)=Cij123(12,8)-Cij124(12,8)
       S3h007311111(3)=-Cij124(12,8)
       S3h007321111(1)=Cij123(12,8)-Cij124(12,8)
       S3h007321111(2)=Cij123(13,8)-Cij124(13,8)
       S3h007321111(3)=-Cij124(13,8)
       S3h007322111(1)=Cij123(13,8)-Cij124(13,8)
       S3h007322111(2)=Cij123(14,8)-Cij124(14,8)
       S3h007322111(3)=-Cij124(14,8)
       S3h007322211(1)=Cij123(14,8)-Cij124(14,8)
       S3h007322211(2)=Cij123(15,8)-Cij124(15,8)
       S3h007322211(3)=-Cij124(15,8)
       S3h007322221(1)=Cij123(15,8)-Cij124(15,8)
       S3h007322221(2)=Cij123(16,8)-Cij124(16,8)
       S3h007322221(3)=-Cij124(16,8)
       S3h007322222(1)=Cij123(16,8)-Cij124(16,8)
       S3h007322222(2)=Cij123(11,8)-Cij124(11,8)
       S3h007322222(3)=-Cij124(11,8)
       S3h007331111(1)=-Cij124(12,8)
       S3h007331111(2)=-Cij124(13,8)
       S3h007331111(3)=-Cij124(13,8)
       S3h007332111(1)=-Cij124(13,8)
       S3h007332111(2)=-Cij124(14,8)
       S3h007332111(3)=-Cij124(14,8)
       S3h007332211(1)=-Cij124(14,8)
       S3h007332211(2)=-Cij124(15,8)
       S3h007332211(3)=-Cij124(15,8)
       S3h007332221(1)=-Cij124(15,8)
       S3h007332221(2)=-Cij124(16,8)
       S3h007332221(3)=-Cij124(16,8)
       S3h007332222(1)=-Cij124(16,8)
       S3h007332222(2)=-Cij124(11,8)
       S3h007332222(3)=-Cij124(11,8)
       S3h007333111(1)=-Cij124(13,8)
       S3h007333111(2)=-Cij124(14,8)
       S3h007333111(3)=-Cij124(14,8)
       S3h007333211(1)=-Cij124(14,8)
       S3h007333211(2)=-Cij124(15,8)
       S3h007333211(3)=-Cij124(15,8)
       S3h007333221(1)=-Cij124(15,8)
       S3h007333221(2)=-Cij124(16,8)
       S3h007333221(3)=-Cij124(16,8)
       S3h007333222(1)=-Cij124(16,8)
       S3h007333222(2)=-Cij124(11,8)
       S3h007333222(3)=-Cij124(11,8)
       S3h007333311(1)=-Cij124(14,8)
       S3h007333311(2)=-Cij124(15,8)
       S3h007333311(3)=-Cij124(15,8)
       S3h007333321(1)=-Cij124(15,8)
       S3h007333321(2)=-Cij124(16,8)
       S3h007333321(3)=-Cij124(16,8)
       S3h007333322(1)=-Cij124(16,8)
       S3h007333322(2)=-Cij124(11,8)
       S3h007333322(3)=-Cij124(11,8)
       S3h007333331(1)=-Cij124(15,8)
       S3h007333331(2)=-Cij124(16,8)
       S3h007333331(3)=-Cij124(16,8)
       S3h007333332(1)=-Cij124(16,8)
       S3h007333332(2)=-Cij124(11,8)
       S3h007333332(3)=-Cij124(11,8)
       S3h007333333(1)=-Cij124(16,8)
       S3h007333333(2)=-Cij124(11,8)
       S3h007333333(3)=-Cij124(11,8)
       aux0000511(1,1,1)=-(F(1)*S3h00611111(1))-F(2)*S3h00621111(1)-F(3)
     &  *S3h00631111(1)+S3h007111111(k)*Z(1,l)+S3h007211111(k)*Z(2,l)+S3
     &  h007311111(k)*Z(3,l)+(Inv4011111+S3000051111(1)-S3h007111111(1)-
     &  S3h007221111(1)-S3h007331111(1))*Z(k,l)-10*S3h000051111(1)*ZZ(k,
     &  1,l,1)-10*S3h000052111(1)*ZZ(k,1,l,2)-10*S3h000053111(1)*ZZ(k,1,
     &  l,3)
       aux0000521(1,1,1)=-(F(1)*S3h00612111(1))-F(2)*S3h00622111(1)-F(3)
     &  *S3h00632111(1)+S3h007121111(k)*Z(1,l)+S3h007221111(k)*Z(2,l)+S3
     &  h007321111(k)*Z(3,l)+(Inv4021111+S3000052111(1)-S3h007121111(1)-
     &  S3h007222111(1)-S3h007332111(1))*Z(k,l)-8*S3h000051211(1)*ZZ(k,1
     &  ,l,1)-8*S3h000052211(1)*ZZ(k,1,l,2)-8*S3h000053211(1)*ZZ(k,1,l,3
     &  )-2*S3h000051111(1)*ZZ(k,2,l,1)-2*S3h000052111(1)*ZZ(k,2,l,2)-2*
     &  S3h000053111(1)*ZZ(k,2,l,3)
       aux0000522(1,1,1)=-(F(1)*S3h00612211(1))-F(2)*S3h00622211(1)-F(3)
     &  *S3h00632211(1)+S3h007122111(k)*Z(1,l)+S3h007222111(k)*Z(2,l)+S3
     &  h007322111(k)*Z(3,l)+(Inv4022111+S3000052211(1)-S3h007122111(1)-
     &  S3h007222211(1)-S3h007332211(1))*Z(k,l)-6*S3h000051221(1)*ZZ(k,1
     &  ,l,1)-6*S3h000052221(1)*ZZ(k,1,l,2)-6*S3h000053221(1)*ZZ(k,1,l,3
     &  )-4*S3h000051211(1)*ZZ(k,2,l,1)-4*S3h000052211(1)*ZZ(k,2,l,2)-4*
     &  S3h000053211(1)*ZZ(k,2,l,3)
       aux0000522(2,1,1)=-(F(1)*S3h00612221(1))-F(2)*S3h00622221(1)-F(3)
     &  *S3h00632221(1)+S3h007122211(k)*Z(1,l)+S3h007222211(k)*Z(2,l)+S3
     &  h007322211(k)*Z(3,l)+(Inv4022211+S3000052221(1)-S3h007122211(1)-
     &  S3h007222221(1)-S3h007332221(1))*Z(k,l)-4*S3h000051222(1)*ZZ(k,1
     &  ,l,1)-4*S3h000052222(1)*ZZ(k,1,l,2)-4*S3h000053222(1)*ZZ(k,1,l,3
     &  )-6*S3h000051221(1)*ZZ(k,2,l,1)-6*S3h000052221(1)*ZZ(k,2,l,2)-6*
     &  S3h000053221(1)*ZZ(k,2,l,3)
       aux0000522(2,2,1)=-(F(1)*S3h00612222(1))-F(2)*S3h00622222(1)-F(3)
     &  *S3h00632222(1)+S3h007122221(k)*Z(1,l)+S3h007222221(k)*Z(2,l)+S3
     &  h007322221(k)*Z(3,l)+(Inv4022221+S3000052222(1)-S3h007122221(1)-
     &  S3h007222222(1)-S3h007332222(1))*Z(k,l)-2*S3h000051222(2)*ZZ(k,1
     &  ,l,1)-2*S3h000052222(2)*ZZ(k,1,l,2)-2*S3h000053222(2)*ZZ(k,1,l,3
     &  )-8*S3h000051222(1)*ZZ(k,2,l,1)-8*S3h000052222(1)*ZZ(k,2,l,2)-8*
     &  S3h000053222(1)*ZZ(k,2,l,3)
       aux0000522(2,2,2)=-(F(1)*S3h00612222(2))-F(2)*S3h00622222(2)-F(3)
     &  *S3h00632222(2)+S3h007122222(k)*Z(1,l)+S3h007222222(k)*Z(2,l)+S3
     &  h007322222(k)*Z(3,l)+(Inv4022222+S3000052222(2)-S3h007122222(1)-
     &  S3h007222222(2)-S3h007332222(2))*Z(k,l)-10*S3h000051222(2)*ZZ(k,
     &  2,l,1)-10*S3h000052222(2)*ZZ(k,2,l,2)-10*S3h000053222(2)*ZZ(k,2,
     &  l,3)
       aux0000531(1,1,1)=-(F(1)*S3h00613111(1))-F(2)*S3h00623111(1)-F(3)
     &  *S3h00633111(1)+S3h007131111(k)*Z(1,l)+S3h007231111(k)*Z(2,l)+S3
     &  h007331111(k)*Z(3,l)+(Inv4031111+S3000053111(1)-S3h007131111(1)-
     &  S3h007232111(1)-S3h007333111(1))*Z(k,l)-8*S3h000051311(1)*ZZ(k,1
     &  ,l,1)-8*S3h000052311(1)*ZZ(k,1,l,2)-8*S3h000053311(1)*ZZ(k,1,l,3
     &  )-2*S3h000051111(1)*ZZ(k,3,l,1)-2*S3h000052111(1)*ZZ(k,3,l,2)-2*
     &  S3h000053111(1)*ZZ(k,3,l,3)
       aux0000532(1,1,1)=-(F(1)*S3h00613211(1))-F(2)*S3h00623211(1)-F(3)
     &  *S3h00633211(1)+S3h007132111(k)*Z(1,l)+S3h007232111(k)*Z(2,l)+S3
     &  h007332111(k)*Z(3,l)+(Inv4032111+S3000053211(1)-S3h007132111(1)-
     &  S3h007232211(1)-S3h007333211(1))*Z(k,l)-6*S3h000051321(1)*ZZ(k,1
     &  ,l,1)-6*S3h000052321(1)*ZZ(k,1,l,2)-6*S3h000053321(1)*ZZ(k,1,l,3
     &  )-2*S3h000051311(1)*ZZ(k,2,l,1)-2*S3h000052311(1)*ZZ(k,2,l,2)-2*
     &  S3h000053311(1)*ZZ(k,2,l,3)-2*S3h000051211(1)*ZZ(k,3,l,1)-2*S3h0
     &  00052211(1)*ZZ(k,3,l,2)-2*S3h000053211(1)*ZZ(k,3,l,3)
       aux0000532(2,1,1)=-(F(1)*S3h00613221(1))-F(2)*S3h00623221(1)-F(3)
     &  *S3h00633221(1)+S3h007132211(k)*Z(1,l)+S3h007232211(k)*Z(2,l)+S3
     &  h007332211(k)*Z(3,l)+(Inv4032211+S3000053221(1)-S3h007132211(1)-
     &  S3h007232221(1)-S3h007333221(1))*Z(k,l)-4*S3h000051322(1)*ZZ(k,1
     &  ,l,1)-4*S3h000052322(1)*ZZ(k,1,l,2)-4*S3h000053322(1)*ZZ(k,1,l,3
     &  )-4*S3h000051321(1)*ZZ(k,2,l,1)-4*S3h000052321(1)*ZZ(k,2,l,2)-4*
     &  S3h000053321(1)*ZZ(k,2,l,3)-2*S3h000051221(1)*ZZ(k,3,l,1)-2*S3h0
     &  00052221(1)*ZZ(k,3,l,2)-2*S3h000053221(1)*ZZ(k,3,l,3)
       aux0000532(2,2,1)=-(F(1)*S3h00613222(1))-F(2)*S3h00623222(1)-F(3)
     &  *S3h00633222(1)+S3h007132221(k)*Z(1,l)+S3h007232221(k)*Z(2,l)+S3
     &  h007332221(k)*Z(3,l)+(Inv4032221+S3000053222(1)-S3h007132221(1)-
     &  S3h007232222(1)-S3h007333222(1))*Z(k,l)-2*S3h000051322(2)*ZZ(k,1
     &  ,l,1)-2*S3h000052322(2)*ZZ(k,1,l,2)-2*S3h000053322(2)*ZZ(k,1,l,3
     &  )-6*S3h000051322(1)*ZZ(k,2,l,1)-6*S3h000052322(1)*ZZ(k,2,l,2)-6*
     &  S3h000053322(1)*ZZ(k,2,l,3)-2*S3h000051222(1)*ZZ(k,3,l,1)-2*S3h0
     &  00052222(1)*ZZ(k,3,l,2)-2*S3h000053222(1)*ZZ(k,3,l,3)
       aux0000532(2,2,2)=-(F(1)*S3h00613222(2))-F(2)*S3h00623222(2)-F(3)
     &  *S3h00633222(2)+S3h007132222(k)*Z(1,l)+S3h007232222(k)*Z(2,l)+S3
     &  h007332222(k)*Z(3,l)+(Inv4032222+S3000053222(2)-S3h007132222(1)-
     &  S3h007232222(2)-S3h007333222(2))*Z(k,l)-8*S3h000051322(2)*ZZ(k,2
     &  ,l,1)-8*S3h000052322(2)*ZZ(k,2,l,2)-8*S3h000053322(2)*ZZ(k,2,l,3
     &  )-2*S3h000051222(2)*ZZ(k,3,l,1)-2*S3h000052222(2)*ZZ(k,3,l,2)-2*
     &  S3h000053222(2)*ZZ(k,3,l,3)
       aux0000533(1,1,1)=-(F(1)*S3h00613311(1))-F(2)*S3h00623311(1)-F(3)
     &  *S3h00633311(1)+S3h007133111(k)*Z(1,l)+S3h007233111(k)*Z(2,l)+S3
     &  h007333111(k)*Z(3,l)+(Inv4033111+S3000053311(1)-S3h007133111(1)-
     &  S3h007233211(1)-S3h007333311(1))*Z(k,l)-6*S3h000051331(1)*ZZ(k,1
     &  ,l,1)-6*S3h000052331(1)*ZZ(k,1,l,2)-6*S3h000053331(1)*ZZ(k,1,l,3
     &  )-4*S3h000051311(1)*ZZ(k,3,l,1)-4*S3h000052311(1)*ZZ(k,3,l,2)-4*
     &  S3h000053311(1)*ZZ(k,3,l,3)
       aux0000533(2,1,1)=-(F(1)*S3h00613321(1))-F(2)*S3h00623321(1)-F(3)
     &  *S3h00633321(1)+S3h007133211(k)*Z(1,l)+S3h007233211(k)*Z(2,l)+S3
     &  h007333211(k)*Z(3,l)+(Inv4033211+S3000053321(1)-S3h007133211(1)-
     &  S3h007233221(1)-S3h007333321(1))*Z(k,l)-4*S3h000051332(1)*ZZ(k,1
     &  ,l,1)-4*S3h000052332(1)*ZZ(k,1,l,2)-4*S3h000053332(1)*ZZ(k,1,l,3
     &  )-2*S3h000051331(1)*ZZ(k,2,l,1)-2*S3h000052331(1)*ZZ(k,2,l,2)-2*
     &  S3h000053331(1)*ZZ(k,2,l,3)-4*S3h000051321(1)*ZZ(k,3,l,1)-4*S3h0
     &  00052321(1)*ZZ(k,3,l,2)-4*S3h000053321(1)*ZZ(k,3,l,3)
       aux0000533(2,2,1)=-(F(1)*S3h00613322(1))-F(2)*S3h00623322(1)-F(3)
     &  *S3h00633322(1)+S3h007133221(k)*Z(1,l)+S3h007233221(k)*Z(2,l)+S3
     &  h007333221(k)*Z(3,l)+(Inv4033221+S3000053322(1)-S3h007133221(1)-
     &  S3h007233222(1)-S3h007333322(1))*Z(k,l)-2*S3h000051332(2)*ZZ(k,1
     &  ,l,1)-2*S3h000052332(2)*ZZ(k,1,l,2)-2*S3h000053332(2)*ZZ(k,1,l,3
     &  )-4*S3h000051332(1)*ZZ(k,2,l,1)-4*S3h000052332(1)*ZZ(k,2,l,2)-4*
     &  S3h000053332(1)*ZZ(k,2,l,3)-4*S3h000051322(1)*ZZ(k,3,l,1)-4*S3h0
     &  00052322(1)*ZZ(k,3,l,2)-4*S3h000053322(1)*ZZ(k,3,l,3)
       aux0000533(2,2,2)=-(F(1)*S3h00613322(2))-F(2)*S3h00623322(2)-F(3)
     &  *S3h00633322(2)+S3h007133222(k)*Z(1,l)+S3h007233222(k)*Z(2,l)+S3
     &  h007333222(k)*Z(3,l)+(Inv4033222+S3000053322(2)-S3h007133222(1)-
     &  S3h007233222(2)-S3h007333322(2))*Z(k,l)-6*S3h000051332(2)*ZZ(k,2
     &  ,l,1)-6*S3h000052332(2)*ZZ(k,2,l,2)-6*S3h000053332(2)*ZZ(k,2,l,3
     &  )-4*S3h000051322(2)*ZZ(k,3,l,1)-4*S3h000052322(2)*ZZ(k,3,l,2)-4*
     &  S3h000053322(2)*ZZ(k,3,l,3)
       aux0000533(3,1,1)=-(F(1)*S3h00613331(1))-F(2)*S3h00623331(1)-F(3)
     &  *S3h00633331(1)+S3h007133311(k)*Z(1,l)+S3h007233311(k)*Z(2,l)+S3
     &  h007333311(k)*Z(3,l)+(Inv4033311+S3000053331(1)-S3h007133311(1)-
     &  S3h007233321(1)-S3h007333331(1))*Z(k,l)-4*S3h000051333(1)*ZZ(k,1
     &  ,l,1)-4*S3h000052333(1)*ZZ(k,1,l,2)-4*S3h000053333(1)*ZZ(k,1,l,3
     &  )-6*S3h000051331(1)*ZZ(k,3,l,1)-6*S3h000052331(1)*ZZ(k,3,l,2)-6*
     &  S3h000053331(1)*ZZ(k,3,l,3)
       aux0000533(3,2,1)=-(F(1)*S3h00613332(1))-F(2)*S3h00623332(1)-F(3)
     &  *S3h00633332(1)+S3h007133321(k)*Z(1,l)+S3h007233321(k)*Z(2,l)+S3
     &  h007333321(k)*Z(3,l)+(Inv4033321+S3000053332(1)-S3h007133321(1)-
     &  S3h007233322(1)-S3h007333332(1))*Z(k,l)-2*S3h000051333(2)*ZZ(k,1
     &  ,l,1)-2*S3h000052333(2)*ZZ(k,1,l,2)-2*S3h000053333(2)*ZZ(k,1,l,3
     &  )-2*S3h000051333(1)*ZZ(k,2,l,1)-2*S3h000052333(1)*ZZ(k,2,l,2)-2*
     &  S3h000053333(1)*ZZ(k,2,l,3)-6*S3h000051332(1)*ZZ(k,3,l,1)-6*S3h0
     &  00052332(1)*ZZ(k,3,l,2)-6*S3h000053332(1)*ZZ(k,3,l,3)
       aux0000533(3,2,2)=-(F(1)*S3h00613332(2))-F(2)*S3h00623332(2)-F(3)
     &  *S3h00633332(2)+S3h007133322(k)*Z(1,l)+S3h007233322(k)*Z(2,l)+S3
     &  h007333322(k)*Z(3,l)+(Inv4033322+S3000053332(2)-S3h007133322(1)-
     &  S3h007233322(2)-S3h007333332(2))*Z(k,l)-4*S3h000051333(2)*ZZ(k,2
     &  ,l,1)-4*S3h000052333(2)*ZZ(k,2,l,2)-4*S3h000053333(2)*ZZ(k,2,l,3
     &  )-6*S3h000051332(2)*ZZ(k,3,l,1)-6*S3h000052332(2)*ZZ(k,3,l,2)-6*
     &  S3h000053332(2)*ZZ(k,3,l,3)
       aux0000533(3,3,1)=-(F(1)*S3h00613333(1))-F(2)*S3h00623333(1)-F(3)
     &  *S3h00633333(1)+S3h007133331(k)*Z(1,l)+S3h007233331(k)*Z(2,l)+S3
     &  h007333331(k)*Z(3,l)+(Inv4033331+S3000053333(1)-S3h007133331(1)-
     &  S3h007233332(1)-S3h007333333(1))*Z(k,l)-2*S3h000051333(3)*ZZ(k,1
     &  ,l,1)-2*S3h000052333(3)*ZZ(k,1,l,2)-2*S3h000053333(3)*ZZ(k,1,l,3
     &  )-8*S3h000051333(1)*ZZ(k,3,l,1)-8*S3h000052333(1)*ZZ(k,3,l,2)-8*
     &  S3h000053333(1)*ZZ(k,3,l,3)
       aux0000533(3,3,2)=-(F(1)*S3h00613333(2))-F(2)*S3h00623333(2)-F(3)
     &  *S3h00633333(2)+S3h007133332(k)*Z(1,l)+S3h007233332(k)*Z(2,l)+S3
     &  h007333332(k)*Z(3,l)+(Inv4033332+S3000053333(2)-S3h007133332(1)-
     &  S3h007233332(2)-S3h007333333(2))*Z(k,l)-2*S3h000051333(3)*ZZ(k,2
     &  ,l,1)-2*S3h000052333(3)*ZZ(k,2,l,2)-2*S3h000053333(3)*ZZ(k,2,l,3
     &  )-8*S3h000051333(2)*ZZ(k,3,l,1)-8*S3h000052333(2)*ZZ(k,3,l,2)-8*
     &  S3h000053333(2)*ZZ(k,3,l,3)
       aux0000533(3,3,3)=-(F(1)*S3h00613333(3))-F(2)*S3h00623333(3)-F(3)
     &  *S3h00633333(3)+S3h007133333(k)*Z(1,l)+S3h007233333(k)*Z(2,l)+S3
     &  h007333333(k)*Z(3,l)+(Inv4033333+S3000053333(3)-S3h007133333(1)-
     &  S3h007233333(2)-S3h007333333(3))*Z(k,l)-10*S3h000051333(3)*ZZ(k,
     &  3,l,1)-10*S3h000052333(3)*ZZ(k,3,l,2)-10*S3h000053333(3)*ZZ(k,3,
     &  l,3)
       temp0000511(1,1,1)=I28Z*(aux0000511(1,1,1)+10*F(4)*temp000041(1,1
     &  ,1)+F(5)*temp00511(1,1,1)+80*temp0000003(1,1,1)*ZZ(k,1,l,1))
       temp0000521(1,1,1)=I28Z*(aux0000521(1,1,1)+2*F(6)*temp000041(1,1,
     &  1)+8*F(4)*temp000042(1,1,1)+F(5)*temp00521(1,1,1)+48*temp0000003
     &  (2,1,1)*ZZ(k,1,l,1)+16*temp0000003(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,
     &  1)))
       temp0000522(1,1,1)=I28Z*(aux0000522(1,1,1)+4*F(6)*temp000042(1,1,
     &  1)+6*F(4)*temp000042(2,1,1)+F(5)*temp00522(1,1,1)+24*(temp000000
     &  3(2,2,1)*ZZ(k,1,l,1)+temp0000003(2,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)
     &  ))+8*temp0000003(1,1,1)*ZZ(k,2,l,2))
       temp0000522(2,1,1)=I28Z*(aux0000522(2,1,1)+6*F(6)*temp000042(2,1,
     &  1)+4*F(4)*temp000042(2,2,1)+F(5)*temp00522(2,1,1)+8*temp0000003(
     &  2,2,2)*ZZ(k,1,l,1)+24*(temp0000003(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,
     &  1))+temp0000003(2,1,1)*ZZ(k,2,l,2)))
       temp0000522(2,2,1)=I28Z*(aux0000522(2,2,1)+8*F(6)*temp000042(2,2,
     &  1)+2*F(4)*temp000042(2,2,2)+F(5)*temp00522(2,2,1)+16*temp0000003
     &  (2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*temp0000003(2,2,1)*ZZ(k,2,l
     &  ,2))
       temp0000522(2,2,2)=I28Z*(aux0000522(2,2,2)+10*F(6)*temp000042(2,2
     &  ,2)+F(5)*temp00522(2,2,2)+80*temp0000003(2,2,2)*ZZ(k,2,l,2))
       temp0000531(1,1,1)=I28Z*(aux0000531(1,1,1)+2*F(7)*temp000041(1,1,
     &  1)+8*F(4)*temp000043(1,1,1)+F(5)*temp00531(1,1,1)+48*temp0000003
     &  (3,1,1)*ZZ(k,1,l,1)+16*temp0000003(1,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,
     &  1)))
       temp0000532(1,1,1)=I28Z*(aux0000532(1,1,1)+2*F(7)*temp000042(1,1,
     &  1)+2*F(6)*temp000043(1,1,1)+6*F(4)*temp000043(2,1,1)+F(5)*temp00
     &  532(1,1,1)+24*temp0000003(3,2,1)*ZZ(k,1,l,1)+12*(temp0000003(3,1
     &  ,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0000003(2,1,1)*(ZZ(k,1,l,3)+ZZ
     &  (k,3,l,1)))+4*temp0000003(1,1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp0000532(2,1,1)=I28Z*(aux0000532(2,1,1)+2*F(7)*temp000042(2,1,
     &  1)+4*F(6)*temp000043(2,1,1)+4*F(4)*temp000043(2,2,1)+F(5)*temp00
     &  532(2,1,1)+16*temp0000003(3,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(te
     &  mp0000003(3,2,2)*ZZ(k,1,l,1)+temp0000003(3,1,1)*ZZ(k,2,l,2))+8*t
     &  emp0000003(2,2,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*temp0000003(2,1,1)
     &  *(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp0000532(2,2,1)=I28Z*(aux0000532(2,2,1)+2*F(7)*temp000042(2,2,
     &  1)+6*F(6)*temp000043(2,2,1)+2*F(4)*temp000043(2,2,2)+F(5)*temp00
     &  532(2,2,1)+24*temp0000003(3,2,1)*ZZ(k,2,l,2)+4*temp0000003(2,2,2
     &  )*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+12*(temp0000003(3,2,2)*(ZZ(k,1,l,2)+
     &  ZZ(k,2,l,1))+temp0000003(2,2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))))
       temp0000532(2,2,2)=I28Z*(aux0000532(2,2,2)+2*F(7)*temp000042(2,2,
     &  2)+8*F(6)*temp000043(2,2,2)+F(5)*temp00532(2,2,2)+48*temp0000003
     &  (3,2,2)*ZZ(k,2,l,2)+16*temp0000003(2,2,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,
     &  2)))
       temp0000533(1,1,1)=I28Z*(aux0000533(1,1,1)+4*F(7)*temp000043(1,1,
     &  1)+6*F(4)*temp000043(3,1,1)+F(5)*temp00533(1,1,1)+24*(temp000000
     &  3(3,3,1)*ZZ(k,1,l,1)+temp0000003(3,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)
     &  ))+8*temp0000003(1,1,1)*ZZ(k,3,l,3))
       temp0000533(2,1,1)=I28Z*(aux0000533(2,1,1)+4*F(7)*temp000043(2,1,
     &  1)+2*F(6)*temp000043(3,1,1)+4*F(4)*temp000043(3,2,1)+F(5)*temp00
     &  533(2,1,1)+8*(temp0000003(3,3,2)*ZZ(k,1,l,1)+temp0000003(3,3,1)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1)))+16*temp0000003(3,2,1)*(ZZ(k,1,l,3)+ZZ
     &  (k,3,l,1))+8*temp0000003(3,1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*temp
     &  0000003(2,1,1)*ZZ(k,3,l,3))
       temp0000533(2,2,1)=I28Z*(aux0000533(2,2,1)+4*F(7)*temp000043(2,2,
     &  1)+4*F(6)*temp000043(3,2,1)+2*F(4)*temp000043(3,2,2)+F(5)*temp00
     &  533(2,2,1)+8*(temp0000003(3,3,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0
     &  000003(3,3,1)*ZZ(k,2,l,2))+8*temp0000003(3,2,2)*(ZZ(k,1,l,3)+ZZ(
     &  k,3,l,1))+16*temp0000003(3,2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*temp
     &  0000003(2,2,1)*ZZ(k,3,l,3))
       temp0000533(2,2,2)=I28Z*(aux0000533(2,2,2)+4*F(7)*temp000043(2,2,
     &  2)+6*F(6)*temp000043(3,2,2)+F(5)*temp00533(2,2,2)+24*(temp000000
     &  3(3,3,2)*ZZ(k,2,l,2)+temp0000003(3,2,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)
     &  ))+8*temp0000003(2,2,2)*ZZ(k,3,l,3))
       temp0000533(3,1,1)=I28Z*(aux0000533(3,1,1)+6*F(7)*temp000043(3,1,
     &  1)+4*F(4)*temp000043(3,3,1)+F(5)*temp00533(3,1,1)+8*temp0000003(
     &  3,3,3)*ZZ(k,1,l,1)+24*(temp0000003(3,3,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,
     &  1))+temp0000003(3,1,1)*ZZ(k,3,l,3)))
       temp0000533(3,2,1)=I28Z*(aux0000533(3,2,1)+6*F(7)*temp000043(3,2,
     &  1)+2*F(6)*temp000043(3,3,1)+2*F(4)*temp000043(3,3,2)+F(5)*temp00
     &  533(3,2,1)+4*temp0000003(3,3,3)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+12*(te
     &  mp0000003(3,3,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+temp0000003(3,3,1)*(Z
     &  Z(k,2,l,3)+ZZ(k,3,l,2)))+24*temp0000003(3,2,1)*ZZ(k,3,l,3))
       temp0000533(3,2,2)=I28Z*(aux0000533(3,2,2)+6*F(7)*temp000043(3,2,
     &  2)+4*F(6)*temp000043(3,3,2)+F(5)*temp00533(3,2,2)+8*temp0000003(
     &  3,3,3)*ZZ(k,2,l,2)+24*(temp0000003(3,3,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,
     &  2))+temp0000003(3,2,2)*ZZ(k,3,l,3)))
       temp0000533(3,3,1)=I28Z*(aux0000533(3,3,1)+8*F(7)*temp000043(3,3,
     &  1)+2*F(4)*temp000043(3,3,3)+F(5)*temp00533(3,3,1)+16*temp0000003
     &  (3,3,3)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+48*temp0000003(3,3,1)*ZZ(k,3,l
     &  ,3))
       temp0000533(3,3,2)=I28Z*(aux0000533(3,3,2)+8*F(7)*temp000043(3,3,
     &  2)+2*F(6)*temp000043(3,3,3)+F(5)*temp00533(3,3,2)+16*temp0000003
     &  (3,3,3)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+48*temp0000003(3,3,2)*ZZ(k,3,l
     &  ,3))
       temp0000533(3,3,3)=I28Z*(aux0000533(3,3,3)+10*F(7)*temp000043(3,3
     &  ,3)+F(5)*temp00533(3,3,3)+80*temp0000003(3,3,3)*ZZ(k,3,l,3))
       temp0000511(1,1,2)=temp0000521(1,1,1)
       temp0000511(1,1,3)=temp0000531(1,1,1)
       temp0000511(1,2,1)=temp0000521(1,1,1)
       temp0000511(1,2,2)=temp0000522(1,1,1)
       temp0000511(1,2,3)=temp0000532(1,1,1)
       temp0000511(1,3,1)=temp0000531(1,1,1)
       temp0000511(1,3,2)=temp0000532(1,1,1)
       temp0000511(1,3,3)=temp0000533(1,1,1)
       temp0000521(1,1,2)=temp0000522(1,1,1)
       temp0000521(1,1,3)=temp0000532(1,1,1)
       temp0000521(1,2,1)=temp0000522(1,1,1)
       temp0000521(1,2,2)=temp0000522(2,1,1)
       temp0000521(1,2,3)=temp0000532(2,1,1)
       temp0000521(1,3,1)=temp0000532(1,1,1)
       temp0000521(1,3,2)=temp0000532(2,1,1)
       temp0000521(1,3,3)=temp0000533(2,1,1)
       temp0000522(1,1,2)=temp0000522(2,1,1)
       temp0000522(1,1,3)=temp0000532(2,1,1)
       temp0000522(1,2,1)=temp0000522(2,1,1)
       temp0000522(1,2,2)=temp0000522(2,2,1)
       temp0000522(1,2,3)=temp0000532(2,2,1)
       temp0000522(1,3,1)=temp0000532(2,1,1)
       temp0000522(1,3,2)=temp0000532(2,2,1)
       temp0000522(1,3,3)=temp0000533(2,2,1)
       temp0000522(2,1,2)=temp0000522(2,2,1)
       temp0000522(2,1,3)=temp0000532(2,2,1)
       temp0000522(2,2,3)=temp0000532(2,2,2)
       temp0000522(2,3,1)=temp0000532(2,2,1)
       temp0000522(2,3,2)=temp0000532(2,2,2)
       temp0000522(2,3,3)=temp0000533(2,2,2)
       temp0000531(1,1,2)=temp0000532(1,1,1)
       temp0000531(1,1,3)=temp0000533(1,1,1)
       temp0000531(1,2,1)=temp0000532(1,1,1)
       temp0000531(1,2,2)=temp0000532(2,1,1)
       temp0000531(1,2,3)=temp0000533(2,1,1)
       temp0000531(1,3,1)=temp0000533(1,1,1)
       temp0000531(1,3,2)=temp0000533(2,1,1)
       temp0000531(1,3,3)=temp0000533(3,1,1)
       temp0000532(1,1,2)=temp0000532(2,1,1)
       temp0000532(1,1,3)=temp0000533(2,1,1)
       temp0000532(1,2,1)=temp0000532(2,1,1)
       temp0000532(1,2,2)=temp0000532(2,2,1)
       temp0000532(1,2,3)=temp0000533(2,2,1)
       temp0000532(1,3,1)=temp0000533(2,1,1)
       temp0000532(1,3,2)=temp0000533(2,2,1)
       temp0000532(1,3,3)=temp0000533(3,2,1)
       temp0000532(2,1,2)=temp0000532(2,2,1)
       temp0000532(2,1,3)=temp0000533(2,2,1)
       temp0000532(2,2,3)=temp0000533(2,2,2)
       temp0000532(2,3,1)=temp0000533(2,2,1)
       temp0000532(2,3,2)=temp0000533(2,2,2)
       temp0000532(2,3,3)=temp0000533(3,2,2)
       temp0000533(1,1,2)=temp0000533(2,1,1)
       temp0000533(1,1,3)=temp0000533(3,1,1)
       temp0000533(1,2,1)=temp0000533(2,1,1)
       temp0000533(1,2,2)=temp0000533(2,2,1)
       temp0000533(1,2,3)=temp0000533(3,2,1)
       temp0000533(1,3,1)=temp0000533(3,1,1)
       temp0000533(1,3,2)=temp0000533(3,2,1)
       temp0000533(1,3,3)=temp0000533(3,3,1)
       temp0000533(2,1,2)=temp0000533(2,2,1)
       temp0000533(2,1,3)=temp0000533(3,2,1)
       temp0000533(2,2,3)=temp0000533(3,2,2)
       temp0000533(2,3,1)=temp0000533(3,2,1)
       temp0000533(2,3,2)=temp0000533(3,2,2)
       temp0000533(2,3,3)=temp0000533(3,3,2)
       temp0000533(3,1,2)=temp0000533(3,2,1)
       temp0000533(3,1,3)=temp0000533(3,3,1)
       temp0000533(3,2,3)=temp0000533(3,3,2)
       S3911111111(1)=-C0234+Cij134(1,8)
       S3911111111(2)=Cij134(1,8)+Cij234(1,1)
       S3911111111(3)=Cij134(3,8)+Cij234(2,1)
       S3912111111(1)=Cij134(1,8)+Cij234(1,1)
       S3912111111(2)=Cij134(1,8)-Cij234(1,2)
       S3912111111(3)=Cij134(3,8)-Cij234(3,2)
       S3912211111(1)=Cij134(1,8)-Cij234(1,2)
       S3912211111(2)=Cij134(1,8)+Cij234(1,3)
       S3912211111(3)=Cij134(3,8)+Cij234(3,3)
       S3912221111(1)=Cij134(1,8)+Cij234(1,3)
       S3912221111(2)=Cij134(1,8)-Cij234(1,4)
       S3912221111(3)=Cij134(3,8)-Cij234(3,4)
       S3912222111(1)=Cij134(1,8)-Cij234(1,4)
       S3912222111(2)=Cij134(1,8)+Cij234(1,5)
       S3912222111(3)=Cij134(3,8)+Cij234(3,5)
       S3912222211(1)=Cij134(1,8)+Cij234(1,5)
       S3912222211(2)=Cij134(1,8)-Cij234(1,6)
       S3912222211(3)=Cij134(3,8)-Cij234(3,6)
       S3912222221(1)=Cij134(1,8)-Cij234(1,6)
       S3912222221(2)=Cij134(1,8)+Cij234(1,7)
       S3912222221(3)=Cij134(3,8)+Cij234(3,7)
       S3912222222(1)=Cij134(1,8)+Cij234(1,7)
       S3912222222(2)=Cij134(1,8)-Cij234(1,8)
       S3912222222(3)=Cij134(3,8)-Cij234(3,8)
       S3913111111(1)=Cij134(3,8)+Cij234(2,1)
       S3913111111(2)=Cij134(3,8)-Cij234(3,2)
       S3913111111(3)=Cij134(4,8)-Cij234(2,2)
       S3913211111(1)=Cij134(3,8)-Cij234(3,2)
       S3913211111(2)=Cij134(3,8)+Cij234(3,3)
       S3913211111(3)=Cij134(4,8)+Cij234(4,3)
       S3913221111(1)=Cij134(3,8)+Cij234(3,3)
       S3913221111(2)=Cij134(3,8)-Cij234(3,4)
       S3913221111(3)=Cij134(4,8)-Cij234(4,4)
       S3913222111(1)=Cij134(3,8)-Cij234(3,4)
       S3913222111(2)=Cij134(3,8)+Cij234(3,5)
       S3913222111(3)=Cij134(4,8)+Cij234(4,5)
       S3913222211(1)=Cij134(3,8)+Cij234(3,5)
       S3913222211(2)=Cij134(3,8)-Cij234(3,6)
       S3913222211(3)=Cij134(4,8)-Cij234(4,6)
       S3913222221(1)=Cij134(3,8)-Cij234(3,6)
       S3913222221(2)=Cij134(3,8)+Cij234(3,7)
       S3913222221(3)=Cij134(4,8)+Cij234(4,7)
       S3913222222(1)=Cij134(3,8)+Cij234(3,7)
       S3913222222(2)=Cij134(3,8)-Cij234(3,8)
       S3913222222(3)=Cij134(4,8)-Cij234(4,8)
       S3913311111(1)=Cij134(4,8)-Cij234(2,2)
       S3913311111(2)=Cij134(4,8)+Cij234(4,3)
       S3913311111(3)=Cij134(5,8)+Cij234(2,3)
       S3913321111(1)=Cij134(4,8)+Cij234(4,3)
       S3913321111(2)=Cij134(4,8)-Cij234(4,4)
       S3913321111(3)=Cij134(5,8)-Cij234(5,4)
       S3913322111(1)=Cij134(4,8)-Cij234(4,4)
       S3913322111(2)=Cij134(4,8)+Cij234(4,5)
       S3913322111(3)=Cij134(5,8)+Cij234(5,5)
       S3913322211(1)=Cij134(4,8)+Cij234(4,5)
       S3913322211(2)=Cij134(4,8)-Cij234(4,6)
       S3913322211(3)=Cij134(5,8)-Cij234(5,6)
       S3913322221(1)=Cij134(4,8)-Cij234(4,6)
       S3913322221(2)=Cij134(4,8)+Cij234(4,7)
       S3913322221(3)=Cij134(5,8)+Cij234(5,7)
       S3913322222(1)=Cij134(4,8)+Cij234(4,7)
       S3913322222(2)=Cij134(4,8)-Cij234(4,8)
       S3913322222(3)=Cij134(5,8)-Cij234(5,8)
       S3913331111(1)=Cij134(5,8)+Cij234(2,3)
       S3913331111(2)=Cij134(5,8)-Cij234(5,4)
       S3913331111(3)=Cij134(6,8)-Cij234(2,4)
       S3913332111(1)=Cij134(5,8)-Cij234(5,4)
       S3913332111(2)=Cij134(5,8)+Cij234(5,5)
       S3913332111(3)=Cij134(6,8)+Cij234(6,5)
       S3913332211(1)=Cij134(5,8)+Cij234(5,5)
       S3913332211(2)=Cij134(5,8)-Cij234(5,6)
       S3913332211(3)=Cij134(6,8)-Cij234(6,6)
       S3913332221(1)=Cij134(5,8)-Cij234(5,6)
       S3913332221(2)=Cij134(5,8)+Cij234(5,7)
       S3913332221(3)=Cij134(6,8)+Cij234(6,7)
       S3913332222(1)=Cij134(5,8)+Cij234(5,7)
       S3913332222(2)=Cij134(5,8)-Cij234(5,8)
       S3913332222(3)=Cij134(6,8)-Cij234(6,8)
       S3913333111(1)=Cij134(6,8)-Cij234(2,4)
       S3913333111(2)=Cij134(6,8)+Cij234(6,5)
       S3913333111(3)=Cij134(7,8)+Cij234(2,5)
       S3913333211(1)=Cij134(6,8)+Cij234(6,5)
       S3913333211(2)=Cij134(6,8)-Cij234(6,6)
       S3913333211(3)=Cij134(7,8)-Cij234(7,6)
       S3913333221(1)=Cij134(6,8)-Cij234(6,6)
       S3913333221(2)=Cij134(6,8)+Cij234(6,7)
       S3913333221(3)=Cij134(7,8)+Cij234(7,7)
       S3913333222(1)=Cij134(6,8)+Cij234(6,7)
       S3913333222(2)=Cij134(6,8)-Cij234(6,8)
       S3913333222(3)=Cij134(7,8)-Cij234(7,8)
       S3913333311(1)=Cij134(7,8)+Cij234(2,5)
       S3913333311(2)=Cij134(7,8)-Cij234(7,6)
       S3913333311(3)=Cij134(8,8)-Cij234(2,6)
       S3913333321(1)=Cij134(7,8)-Cij234(7,6)
       S3913333321(2)=Cij134(7,8)+Cij234(7,7)
       S3913333321(3)=Cij134(8,8)+Cij234(8,7)
       S3913333322(1)=Cij134(7,8)+Cij234(7,7)
       S3913333322(2)=Cij134(7,8)-Cij234(7,8)
       S3913333322(3)=Cij134(8,8)-Cij234(8,8)
       S3913333331(1)=Cij134(8,8)-Cij234(2,6)
       S3913333331(2)=Cij134(8,8)+Cij234(8,7)
       S3913333331(3)=Cij134(9,8)+Cij234(2,7)
       S3913333332(1)=Cij134(8,8)+Cij234(8,7)
       S3913333332(2)=Cij134(8,8)-Cij234(8,8)
       S3913333332(3)=Cij134(9,8)-Cij234(9,8)
       S3913333333(1)=Cij134(9,8)+Cij234(2,7)
       S3913333333(2)=Cij134(9,8)-Cij234(9,8)
       S3913333333(3)=Cij134(2,8)-Cij234(2,8)
       S3921111111(1)=Cij124(1,8)-Cij134(1,8)
       S3921111111(2)=Cij124(3,8)-Cij134(1,8)
       S3921111111(3)=Cij124(3,8)-Cij134(3,8)
       S3922111111(1)=Cij124(3,8)-Cij134(1,8)
       S3922111111(2)=Cij124(4,8)-Cij134(1,8)
       S3922111111(3)=Cij124(4,8)-Cij134(3,8)
       S3922211111(1)=Cij124(4,8)-Cij134(1,8)
       S3922211111(2)=Cij124(5,8)-Cij134(1,8)
       S3922211111(3)=Cij124(5,8)-Cij134(3,8)
       S3922221111(1)=Cij124(5,8)-Cij134(1,8)
       S3922221111(2)=Cij124(6,8)-Cij134(1,8)
       S3922221111(3)=Cij124(6,8)-Cij134(3,8)
       S3922222111(1)=Cij124(6,8)-Cij134(1,8)
       S3922222111(2)=Cij124(7,8)-Cij134(1,8)
       S3922222111(3)=Cij124(7,8)-Cij134(3,8)
       S3922222211(1)=Cij124(7,8)-Cij134(1,8)
       S3922222211(2)=Cij124(8,8)-Cij134(1,8)
       S3922222211(3)=Cij124(8,8)-Cij134(3,8)
       S3922222221(1)=Cij124(8,8)-Cij134(1,8)
       S3922222221(2)=Cij124(9,8)-Cij134(1,8)
       S3922222221(3)=Cij124(9,8)-Cij134(3,8)
       S3922222222(1)=Cij124(9,8)-Cij134(1,8)
       S3922222222(2)=Cij124(2,8)-Cij134(1,8)
       S3922222222(3)=Cij124(2,8)-Cij134(3,8)
       S3923111111(1)=Cij124(3,8)-Cij134(3,8)
       S3923111111(2)=Cij124(4,8)-Cij134(3,8)
       S3923111111(3)=Cij124(4,8)-Cij134(4,8)
       S3923211111(1)=Cij124(4,8)-Cij134(3,8)
       S3923211111(2)=Cij124(5,8)-Cij134(3,8)
       S3923211111(3)=Cij124(5,8)-Cij134(4,8)
       S3923221111(1)=Cij124(5,8)-Cij134(3,8)
       S3923221111(2)=Cij124(6,8)-Cij134(3,8)
       S3923221111(3)=Cij124(6,8)-Cij134(4,8)
       S3923222111(1)=Cij124(6,8)-Cij134(3,8)
       S3923222111(2)=Cij124(7,8)-Cij134(3,8)
       S3923222111(3)=Cij124(7,8)-Cij134(4,8)
       S3923222211(1)=Cij124(7,8)-Cij134(3,8)
       S3923222211(2)=Cij124(8,8)-Cij134(3,8)
       S3923222211(3)=Cij124(8,8)-Cij134(4,8)
       S3923222221(1)=Cij124(8,8)-Cij134(3,8)
       S3923222221(2)=Cij124(9,8)-Cij134(3,8)
       S3923222221(3)=Cij124(9,8)-Cij134(4,8)
       S3923222222(1)=Cij124(9,8)-Cij134(3,8)
       S3923222222(2)=Cij124(2,8)-Cij134(3,8)
       S3923222222(3)=Cij124(2,8)-Cij134(4,8)
       S3923311111(1)=Cij124(4,8)-Cij134(4,8)
       S3923311111(2)=Cij124(5,8)-Cij134(4,8)
       S3923311111(3)=Cij124(5,8)-Cij134(5,8)
       S3923321111(1)=Cij124(5,8)-Cij134(4,8)
       S3923321111(2)=Cij124(6,8)-Cij134(4,8)
       S3923321111(3)=Cij124(6,8)-Cij134(5,8)
       S3923322111(1)=Cij124(6,8)-Cij134(4,8)
       S3923322111(2)=Cij124(7,8)-Cij134(4,8)
       S3923322111(3)=Cij124(7,8)-Cij134(5,8)
       S3923322211(1)=Cij124(7,8)-Cij134(4,8)
       S3923322211(2)=Cij124(8,8)-Cij134(4,8)
       S3923322211(3)=Cij124(8,8)-Cij134(5,8)
       S3923322221(1)=Cij124(8,8)-Cij134(4,8)
       S3923322221(2)=Cij124(9,8)-Cij134(4,8)
       S3923322221(3)=Cij124(9,8)-Cij134(5,8)
       S3923322222(1)=Cij124(9,8)-Cij134(4,8)
       S3923322222(2)=Cij124(2,8)-Cij134(4,8)
       S3923322222(3)=Cij124(2,8)-Cij134(5,8)
       S3923331111(1)=Cij124(5,8)-Cij134(5,8)
       S3923331111(2)=Cij124(6,8)-Cij134(5,8)
       S3923331111(3)=Cij124(6,8)-Cij134(6,8)
       S3923332111(1)=Cij124(6,8)-Cij134(5,8)
       S3923332111(2)=Cij124(7,8)-Cij134(5,8)
       S3923332111(3)=Cij124(7,8)-Cij134(6,8)
       S3923332211(1)=Cij124(7,8)-Cij134(5,8)
       S3923332211(2)=Cij124(8,8)-Cij134(5,8)
       S3923332211(3)=Cij124(8,8)-Cij134(6,8)
       S3923332221(1)=Cij124(8,8)-Cij134(5,8)
       S3923332221(2)=Cij124(9,8)-Cij134(5,8)
       S3923332221(3)=Cij124(9,8)-Cij134(6,8)
       S3923332222(1)=Cij124(9,8)-Cij134(5,8)
       S3923332222(2)=Cij124(2,8)-Cij134(5,8)
       S3923332222(3)=Cij124(2,8)-Cij134(6,8)
       S3923333111(1)=Cij124(6,8)-Cij134(6,8)
       S3923333111(2)=Cij124(7,8)-Cij134(6,8)
       S3923333111(3)=Cij124(7,8)-Cij134(7,8)
       S3923333211(1)=Cij124(7,8)-Cij134(6,8)
       S3923333211(2)=Cij124(8,8)-Cij134(6,8)
       S3923333211(3)=Cij124(8,8)-Cij134(7,8)
       S3923333221(1)=Cij124(8,8)-Cij134(6,8)
       S3923333221(2)=Cij124(9,8)-Cij134(6,8)
       S3923333221(3)=Cij124(9,8)-Cij134(7,8)
       S3923333222(1)=Cij124(9,8)-Cij134(6,8)
       S3923333222(2)=Cij124(2,8)-Cij134(6,8)
       S3923333222(3)=Cij124(2,8)-Cij134(7,8)
       S3923333311(1)=Cij124(7,8)-Cij134(7,8)
       S3923333311(2)=Cij124(8,8)-Cij134(7,8)
       S3923333311(3)=Cij124(8,8)-Cij134(8,8)
       S3923333321(1)=Cij124(8,8)-Cij134(7,8)
       S3923333321(2)=Cij124(9,8)-Cij134(7,8)
       S3923333321(3)=Cij124(9,8)-Cij134(8,8)
       S3923333322(1)=Cij124(9,8)-Cij134(7,8)
       S3923333322(2)=Cij124(2,8)-Cij134(7,8)
       S3923333322(3)=Cij124(2,8)-Cij134(8,8)
       S3923333331(1)=Cij124(8,8)-Cij134(8,8)
       S3923333331(2)=Cij124(9,8)-Cij134(8,8)
       S3923333331(3)=Cij124(9,8)-Cij134(9,8)
       S3923333332(1)=Cij124(9,8)-Cij134(8,8)
       S3923333332(2)=Cij124(2,8)-Cij134(8,8)
       S3923333332(3)=Cij124(2,8)-Cij134(9,8)
       S3923333333(1)=Cij124(9,8)-Cij134(9,8)
       S3923333333(2)=Cij124(2,8)-Cij134(9,8)
       S3923333333(3)=Cij124(2,8)-Cij134(2,8)
       S3931111111(1)=Cij123(1,8)-Cij124(1,8)
       S3931111111(2)=Cij123(3,8)-Cij124(3,8)
       S3931111111(3)=-Cij124(3,8)
       S3932111111(1)=Cij123(3,8)-Cij124(3,8)
       S3932111111(2)=Cij123(4,8)-Cij124(4,8)
       S3932111111(3)=-Cij124(4,8)
       S3932211111(1)=Cij123(4,8)-Cij124(4,8)
       S3932211111(2)=Cij123(5,8)-Cij124(5,8)
       S3932211111(3)=-Cij124(5,8)
       S3932221111(1)=Cij123(5,8)-Cij124(5,8)
       S3932221111(2)=Cij123(6,8)-Cij124(6,8)
       S3932221111(3)=-Cij124(6,8)
       S3932222111(1)=Cij123(6,8)-Cij124(6,8)
       S3932222111(2)=Cij123(7,8)-Cij124(7,8)
       S3932222111(3)=-Cij124(7,8)
       S3932222211(1)=Cij123(7,8)-Cij124(7,8)
       S3932222211(2)=Cij123(8,8)-Cij124(8,8)
       S3932222211(3)=-Cij124(8,8)
       S3932222221(1)=Cij123(8,8)-Cij124(8,8)
       S3932222221(2)=Cij123(9,8)-Cij124(9,8)
       S3932222221(3)=-Cij124(9,8)
       S3932222222(1)=Cij123(9,8)-Cij124(9,8)
       S3932222222(2)=Cij123(2,8)-Cij124(2,8)
       S3932222222(3)=-Cij124(2,8)
       S3933111111(1)=-Cij124(3,8)
       S3933111111(2)=-Cij124(4,8)
       S3933111111(3)=-Cij124(4,8)
       S3933211111(1)=-Cij124(4,8)
       S3933211111(2)=-Cij124(5,8)
       S3933211111(3)=-Cij124(5,8)
       S3933221111(1)=-Cij124(5,8)
       S3933221111(2)=-Cij124(6,8)
       S3933221111(3)=-Cij124(6,8)
       S3933222111(1)=-Cij124(6,8)
       S3933222111(2)=-Cij124(7,8)
       S3933222111(3)=-Cij124(7,8)
       S3933222211(1)=-Cij124(7,8)
       S3933222211(2)=-Cij124(8,8)
       S3933222211(3)=-Cij124(8,8)
       S3933222221(1)=-Cij124(8,8)
       S3933222221(2)=-Cij124(9,8)
       S3933222221(3)=-Cij124(9,8)
       S3933222222(1)=-Cij124(9,8)
       S3933222222(2)=-Cij124(2,8)
       S3933222222(3)=-Cij124(2,8)
       S3933311111(1)=-Cij124(4,8)
       S3933311111(2)=-Cij124(5,8)
       S3933311111(3)=-Cij124(5,8)
       S3933321111(1)=-Cij124(5,8)
       S3933321111(2)=-Cij124(6,8)
       S3933321111(3)=-Cij124(6,8)
       S3933322111(1)=-Cij124(6,8)
       S3933322111(2)=-Cij124(7,8)
       S3933322111(3)=-Cij124(7,8)
       S3933322211(1)=-Cij124(7,8)
       S3933322211(2)=-Cij124(8,8)
       S3933322211(3)=-Cij124(8,8)
       S3933322221(1)=-Cij124(8,8)
       S3933322221(2)=-Cij124(9,8)
       S3933322221(3)=-Cij124(9,8)
       S3933322222(1)=-Cij124(9,8)
       S3933322222(2)=-Cij124(2,8)
       S3933322222(3)=-Cij124(2,8)
       S3933331111(1)=-Cij124(5,8)
       S3933331111(2)=-Cij124(6,8)
       S3933331111(3)=-Cij124(6,8)
       S3933332111(1)=-Cij124(6,8)
       S3933332111(2)=-Cij124(7,8)
       S3933332111(3)=-Cij124(7,8)
       S3933332211(1)=-Cij124(7,8)
       S3933332211(2)=-Cij124(8,8)
       S3933332211(3)=-Cij124(8,8)
       S3933332221(1)=-Cij124(8,8)
       S3933332221(2)=-Cij124(9,8)
       S3933332221(3)=-Cij124(9,8)
       S3933332222(1)=-Cij124(9,8)
       S3933332222(2)=-Cij124(2,8)
       S3933332222(3)=-Cij124(2,8)
       S3933333111(1)=-Cij124(6,8)
       S3933333111(2)=-Cij124(7,8)
       S3933333111(3)=-Cij124(7,8)
       S3933333211(1)=-Cij124(7,8)
       S3933333211(2)=-Cij124(8,8)
       S3933333211(3)=-Cij124(8,8)
       S3933333221(1)=-Cij124(8,8)
       S3933333221(2)=-Cij124(9,8)
       S3933333221(3)=-Cij124(9,8)
       S3933333222(1)=-Cij124(9,8)
       S3933333222(2)=-Cij124(2,8)
       S3933333222(3)=-Cij124(2,8)
       S3933333311(1)=-Cij124(7,8)
       S3933333311(2)=-Cij124(8,8)
       S3933333311(3)=-Cij124(8,8)
       S3933333321(1)=-Cij124(8,8)
       S3933333321(2)=-Cij124(9,8)
       S3933333321(3)=-Cij124(9,8)
       S3933333322(1)=-Cij124(9,8)
       S3933333322(2)=-Cij124(2,8)
       S3933333322(3)=-Cij124(2,8)
       S3933333331(1)=-Cij124(8,8)
       S3933333331(2)=-Cij124(9,8)
       S3933333331(3)=-Cij124(9,8)
       S3933333332(1)=-Cij124(9,8)
       S3933333332(2)=-Cij124(2,8)
       S3933333332(3)=-Cij124(2,8)
       S3933333333(1)=-Cij124(9,8)
       S3933333333(2)=-Cij124(2,8)
       S3933333333(3)=-Cij124(2,8)
       aux0071111(1,1,1)=-(F(1)*S381111111(1))-F(2)*S382111111(1)-F(3)*S
     &  383111111(1)+S3911111111(k)*Z(1,l)+S3921111111(k)*Z(2,l)+S393111
     &  1111(k)*Z(3,l)+(S3007111111(1)-S3911111111(1)-S3922111111(1)-S39
     &  33111111(1))*Z(k,l)-14*S3h007111111(1)*ZZ(k,1,l,1)-14*S3h0072111
     &  11(1)*ZZ(k,1,l,2)-14*S3h007311111(1)*ZZ(k,1,l,3)
       aux0072111(1,1,1)=-(F(1)*S381211111(1))-F(2)*S382211111(1)-F(3)*S
     &  383211111(1)+S3912111111(k)*Z(1,l)+S3922111111(k)*Z(2,l)+S393211
     &  1111(k)*Z(3,l)+(S3007211111(1)-S3912111111(1)-S3922211111(1)-S39
     &  33211111(1))*Z(k,l)-12*S3h007121111(1)*ZZ(k,1,l,1)-12*S3h0072211
     &  11(1)*ZZ(k,1,l,2)-12*S3h007321111(1)*ZZ(k,1,l,3)-2*S3h007111111(
     &  1)*ZZ(k,2,l,1)-2*S3h007211111(1)*ZZ(k,2,l,2)-2*S3h007311111(1)*Z
     &  Z(k,2,l,3)
       aux0072211(1,1,1)=-(F(1)*S381221111(1))-F(2)*S382221111(1)-F(3)*S
     &  383221111(1)+S3912211111(k)*Z(1,l)+S3922211111(k)*Z(2,l)+S393221
     &  1111(k)*Z(3,l)+(S3007221111(1)-S3912211111(1)-S3922221111(1)-S39
     &  33221111(1))*Z(k,l)-10*S3h007122111(1)*ZZ(k,1,l,1)-10*S3h0072221
     &  11(1)*ZZ(k,1,l,2)-10*S3h007322111(1)*ZZ(k,1,l,3)-4*S3h007121111(
     &  1)*ZZ(k,2,l,1)-4*S3h007221111(1)*ZZ(k,2,l,2)-4*S3h007321111(1)*Z
     &  Z(k,2,l,3)
       aux0072221(1,1,1)=-(F(1)*S381222111(1))-F(2)*S382222111(1)-F(3)*S
     &  383222111(1)+S3912221111(k)*Z(1,l)+S3922221111(k)*Z(2,l)+S393222
     &  1111(k)*Z(3,l)+(S3007222111(1)-S3912221111(1)-S3922222111(1)-S39
     &  33222111(1))*Z(k,l)-8*S3h007122211(1)*ZZ(k,1,l,1)-8*S3h007222211
     &  (1)*ZZ(k,1,l,2)-8*S3h007322211(1)*ZZ(k,1,l,3)-6*S3h007122111(1)*
     &  ZZ(k,2,l,1)-6*S3h007222111(1)*ZZ(k,2,l,2)-6*S3h007322111(1)*ZZ(k
     &  ,2,l,3)
       aux0072222(1,1,1)=-(F(1)*S381222211(1))-F(2)*S382222211(1)-F(3)*S
     &  383222211(1)+S3912222111(k)*Z(1,l)+S3922222111(k)*Z(2,l)+S393222
     &  2111(k)*Z(3,l)+(S3007222211(1)-S3912222111(1)-S3922222211(1)-S39
     &  33222211(1))*Z(k,l)-6*S3h007122221(1)*ZZ(k,1,l,1)-6*S3h007222221
     &  (1)*ZZ(k,1,l,2)-6*S3h007322221(1)*ZZ(k,1,l,3)-8*S3h007122211(1)*
     &  ZZ(k,2,l,1)-8*S3h007222211(1)*ZZ(k,2,l,2)-8*S3h007322211(1)*ZZ(k
     &  ,2,l,3)
       aux0072222(2,1,1)=-(F(1)*S381222221(1))-F(2)*S382222221(1)-F(3)*S
     &  383222221(1)+S3912222211(k)*Z(1,l)+S3922222211(k)*Z(2,l)+S393222
     &  2211(k)*Z(3,l)+(S3007222221(1)-S3912222211(1)-S3922222221(1)-S39
     &  33222221(1))*Z(k,l)-4*S3h007122222(1)*ZZ(k,1,l,1)-4*S3h007222222
     &  (1)*ZZ(k,1,l,2)-4*S3h007322222(1)*ZZ(k,1,l,3)-10*S3h007122221(1)
     &  *ZZ(k,2,l,1)-10*S3h007222221(1)*ZZ(k,2,l,2)-10*S3h007322221(1)*Z
     &  Z(k,2,l,3)
       aux0072222(2,2,1)=-(F(1)*S381222222(1))-F(2)*S382222222(1)-F(3)*S
     &  383222222(1)+S3912222221(k)*Z(1,l)+S3922222221(k)*Z(2,l)+S393222
     &  2221(k)*Z(3,l)+(S3007222222(1)-S3912222221(1)-S3922222222(1)-S39
     &  33222222(1))*Z(k,l)-2*S3h007122222(2)*ZZ(k,1,l,1)-2*S3h007222222
     &  (2)*ZZ(k,1,l,2)-2*S3h007322222(2)*ZZ(k,1,l,3)-12*S3h007122222(1)
     &  *ZZ(k,2,l,1)-12*S3h007222222(1)*ZZ(k,2,l,2)-12*S3h007322222(1)*Z
     &  Z(k,2,l,3)
       aux0072222(2,2,2)=-(F(1)*S381222222(2))-F(2)*S382222222(2)-F(3)*S
     &  383222222(2)+S3912222222(k)*Z(1,l)+S3922222222(k)*Z(2,l)+S393222
     &  2222(k)*Z(3,l)+(S3007222222(2)-S3912222222(1)-S3922222222(2)-S39
     &  33222222(2))*Z(k,l)-14*S3h007122222(2)*ZZ(k,2,l,1)-14*S3h0072222
     &  22(2)*ZZ(k,2,l,2)-14*S3h007322222(2)*ZZ(k,2,l,3)
       aux0073111(1,1,1)=-(F(1)*S381311111(1))-F(2)*S382311111(1)-F(3)*S
     &  383311111(1)+S3913111111(k)*Z(1,l)+S3923111111(k)*Z(2,l)+S393311
     &  1111(k)*Z(3,l)+(S3007311111(1)-S3913111111(1)-S3923211111(1)-S39
     &  33311111(1))*Z(k,l)-12*S3h007131111(1)*ZZ(k,1,l,1)-12*S3h0072311
     &  11(1)*ZZ(k,1,l,2)-12*S3h007331111(1)*ZZ(k,1,l,3)-2*S3h007111111(
     &  1)*ZZ(k,3,l,1)-2*S3h007211111(1)*ZZ(k,3,l,2)-2*S3h007311111(1)*Z
     &  Z(k,3,l,3)
       aux0073211(1,1,1)=-(F(1)*S381321111(1))-F(2)*S382321111(1)-F(3)*S
     &  383321111(1)+S3913211111(k)*Z(1,l)+S3923211111(k)*Z(2,l)+S393321
     &  1111(k)*Z(3,l)+(S3007321111(1)-S3913211111(1)-S3923221111(1)-S39
     &  33321111(1))*Z(k,l)-10*S3h007132111(1)*ZZ(k,1,l,1)-10*S3h0072321
     &  11(1)*ZZ(k,1,l,2)-10*S3h007332111(1)*ZZ(k,1,l,3)-2*S3h007131111(
     &  1)*ZZ(k,2,l,1)-2*S3h007231111(1)*ZZ(k,2,l,2)-2*S3h007331111(1)*Z
     &  Z(k,2,l,3)-2*S3h007121111(1)*ZZ(k,3,l,1)-2*S3h007221111(1)*ZZ(k,
     &  3,l,2)-2*S3h007321111(1)*ZZ(k,3,l,3)
       aux0073221(1,1,1)=-(F(1)*S381322111(1))-F(2)*S382322111(1)-F(3)*S
     &  383322111(1)+S3913221111(k)*Z(1,l)+S3923221111(k)*Z(2,l)+S393322
     &  1111(k)*Z(3,l)+(S3007322111(1)-S3913221111(1)-S3923222111(1)-S39
     &  33322111(1))*Z(k,l)-8*S3h007132211(1)*ZZ(k,1,l,1)-8*S3h007232211
     &  (1)*ZZ(k,1,l,2)-8*S3h007332211(1)*ZZ(k,1,l,3)-4*S3h007132111(1)*
     &  ZZ(k,2,l,1)-4*S3h007232111(1)*ZZ(k,2,l,2)-4*S3h007332111(1)*ZZ(k
     &  ,2,l,3)-2*S3h007122111(1)*ZZ(k,3,l,1)-2*S3h007222111(1)*ZZ(k,3,l
     &  ,2)-2*S3h007322111(1)*ZZ(k,3,l,3)
       aux0073222(1,1,1)=-(F(1)*S381322211(1))-F(2)*S382322211(1)-F(3)*S
     &  383322211(1)+S3913222111(k)*Z(1,l)+S3923222111(k)*Z(2,l)+S393322
     &  2111(k)*Z(3,l)+(S3007322211(1)-S3913222111(1)-S3923222211(1)-S39
     &  33322211(1))*Z(k,l)-6*S3h007132221(1)*ZZ(k,1,l,1)-6*S3h007232221
     &  (1)*ZZ(k,1,l,2)-6*S3h007332221(1)*ZZ(k,1,l,3)-6*S3h007132211(1)*
     &  ZZ(k,2,l,1)-6*S3h007232211(1)*ZZ(k,2,l,2)-6*S3h007332211(1)*ZZ(k
     &  ,2,l,3)-2*S3h007122211(1)*ZZ(k,3,l,1)-2*S3h007222211(1)*ZZ(k,3,l
     &  ,2)-2*S3h007322211(1)*ZZ(k,3,l,3)
       aux0073222(2,1,1)=-(F(1)*S381322221(1))-F(2)*S382322221(1)-F(3)*S
     &  383322221(1)+S3913222211(k)*Z(1,l)+S3923222211(k)*Z(2,l)+S393322
     &  2211(k)*Z(3,l)+(S3007322221(1)-S3913222211(1)-S3923222221(1)-S39
     &  33322221(1))*Z(k,l)-4*S3h007132222(1)*ZZ(k,1,l,1)-4*S3h007232222
     &  (1)*ZZ(k,1,l,2)-4*S3h007332222(1)*ZZ(k,1,l,3)-8*S3h007132221(1)*
     &  ZZ(k,2,l,1)-8*S3h007232221(1)*ZZ(k,2,l,2)-8*S3h007332221(1)*ZZ(k
     &  ,2,l,3)-2*S3h007122221(1)*ZZ(k,3,l,1)-2*S3h007222221(1)*ZZ(k,3,l
     &  ,2)-2*S3h007322221(1)*ZZ(k,3,l,3)
       aux0073222(2,2,1)=-(F(1)*S381322222(1))-F(2)*S382322222(1)-F(3)*S
     &  383322222(1)+S3913222221(k)*Z(1,l)+S3923222221(k)*Z(2,l)+S393322
     &  2221(k)*Z(3,l)+(S3007322222(1)-S3913222221(1)-S3923222222(1)-S39
     &  33322222(1))*Z(k,l)-2*S3h007132222(2)*ZZ(k,1,l,1)-2*S3h007232222
     &  (2)*ZZ(k,1,l,2)-2*S3h007332222(2)*ZZ(k,1,l,3)-10*S3h007132222(1)
     &  *ZZ(k,2,l,1)-10*S3h007232222(1)*ZZ(k,2,l,2)-10*S3h007332222(1)*Z
     &  Z(k,2,l,3)-2*S3h007122222(1)*ZZ(k,3,l,1)-2*S3h007222222(1)*ZZ(k,
     &  3,l,2)-2*S3h007322222(1)*ZZ(k,3,l,3)
       aux0073222(2,2,2)=-(F(1)*S381322222(2))-F(2)*S382322222(2)-F(3)*S
     &  383322222(2)+S3913222222(k)*Z(1,l)+S3923222222(k)*Z(2,l)+S393322
     &  2222(k)*Z(3,l)+(S3007322222(2)-S3913222222(1)-S3923222222(2)-S39
     &  33322222(2))*Z(k,l)-12*S3h007132222(2)*ZZ(k,2,l,1)-12*S3h0072322
     &  22(2)*ZZ(k,2,l,2)-12*S3h007332222(2)*ZZ(k,2,l,3)-2*S3h007122222(
     &  2)*ZZ(k,3,l,1)-2*S3h007222222(2)*ZZ(k,3,l,2)-2*S3h007322222(2)*Z
     &  Z(k,3,l,3)
       aux0073311(1,1,1)=-(F(1)*S381331111(1))-F(2)*S382331111(1)-F(3)*S
     &  383331111(1)+S3913311111(k)*Z(1,l)+S3923311111(k)*Z(2,l)+S393331
     &  1111(k)*Z(3,l)+(S3007331111(1)-S3913311111(1)-S3923321111(1)-S39
     &  33331111(1))*Z(k,l)-10*S3h007133111(1)*ZZ(k,1,l,1)-10*S3h0072331
     &  11(1)*ZZ(k,1,l,2)-10*S3h007333111(1)*ZZ(k,1,l,3)-4*S3h007131111(
     &  1)*ZZ(k,3,l,1)-4*S3h007231111(1)*ZZ(k,3,l,2)-4*S3h007331111(1)*Z
     &  Z(k,3,l,3)
       aux0073321(1,1,1)=-(F(1)*S381332111(1))-F(2)*S382332111(1)-F(3)*S
     &  383332111(1)+S3913321111(k)*Z(1,l)+S3923321111(k)*Z(2,l)+S393332
     &  1111(k)*Z(3,l)+(S3007332111(1)-S3913321111(1)-S3923322111(1)-S39
     &  33332111(1))*Z(k,l)-8*S3h007133211(1)*ZZ(k,1,l,1)-8*S3h007233211
     &  (1)*ZZ(k,1,l,2)-8*S3h007333211(1)*ZZ(k,1,l,3)-2*S3h007133111(1)*
     &  ZZ(k,2,l,1)-2*S3h007233111(1)*ZZ(k,2,l,2)-2*S3h007333111(1)*ZZ(k
     &  ,2,l,3)-4*S3h007132111(1)*ZZ(k,3,l,1)-4*S3h007232111(1)*ZZ(k,3,l
     &  ,2)-4*S3h007332111(1)*ZZ(k,3,l,3)
       aux0073322(1,1,1)=-(F(1)*S381332211(1))-F(2)*S382332211(1)-F(3)*S
     &  383332211(1)+S3913322111(k)*Z(1,l)+S3923322111(k)*Z(2,l)+S393332
     &  2111(k)*Z(3,l)+(S3007332211(1)-S3913322111(1)-S3923322211(1)-S39
     &  33332211(1))*Z(k,l)-6*S3h007133221(1)*ZZ(k,1,l,1)-6*S3h007233221
     &  (1)*ZZ(k,1,l,2)-6*S3h007333221(1)*ZZ(k,1,l,3)-4*S3h007133211(1)*
     &  ZZ(k,2,l,1)-4*S3h007233211(1)*ZZ(k,2,l,2)-4*S3h007333211(1)*ZZ(k
     &  ,2,l,3)-4*S3h007132211(1)*ZZ(k,3,l,1)-4*S3h007232211(1)*ZZ(k,3,l
     &  ,2)-4*S3h007332211(1)*ZZ(k,3,l,3)
       aux0073322(2,1,1)=-(F(1)*S381332221(1))-F(2)*S382332221(1)-F(3)*S
     &  383332221(1)+S3913322211(k)*Z(1,l)+S3923322211(k)*Z(2,l)+S393332
     &  2211(k)*Z(3,l)+(S3007332221(1)-S3913322211(1)-S3923322221(1)-S39
     &  33332221(1))*Z(k,l)-4*S3h007133222(1)*ZZ(k,1,l,1)-4*S3h007233222
     &  (1)*ZZ(k,1,l,2)-4*S3h007333222(1)*ZZ(k,1,l,3)-6*S3h007133221(1)*
     &  ZZ(k,2,l,1)-6*S3h007233221(1)*ZZ(k,2,l,2)-6*S3h007333221(1)*ZZ(k
     &  ,2,l,3)-4*S3h007132221(1)*ZZ(k,3,l,1)-4*S3h007232221(1)*ZZ(k,3,l
     &  ,2)-4*S3h007332221(1)*ZZ(k,3,l,3)
       aux0073322(2,2,1)=-(F(1)*S381332222(1))-F(2)*S382332222(1)-F(3)*S
     &  383332222(1)+S3913322221(k)*Z(1,l)+S3923322221(k)*Z(2,l)+S393332
     &  2221(k)*Z(3,l)+(S3007332222(1)-S3913322221(1)-S3923322222(1)-S39
     &  33332222(1))*Z(k,l)-2*S3h007133222(2)*ZZ(k,1,l,1)-2*S3h007233222
     &  (2)*ZZ(k,1,l,2)-2*S3h007333222(2)*ZZ(k,1,l,3)-8*S3h007133222(1)*
     &  ZZ(k,2,l,1)-8*S3h007233222(1)*ZZ(k,2,l,2)-8*S3h007333222(1)*ZZ(k
     &  ,2,l,3)-4*S3h007132222(1)*ZZ(k,3,l,1)-4*S3h007232222(1)*ZZ(k,3,l
     &  ,2)-4*S3h007332222(1)*ZZ(k,3,l,3)
       aux0073322(2,2,2)=-(F(1)*S381332222(2))-F(2)*S382332222(2)-F(3)*S
     &  383332222(2)+S3913322222(k)*Z(1,l)+S3923322222(k)*Z(2,l)+S393332
     &  2222(k)*Z(3,l)+(S3007332222(2)-S3913322222(1)-S3923322222(2)-S39
     &  33332222(2))*Z(k,l)-10*S3h007133222(2)*ZZ(k,2,l,1)-10*S3h0072332
     &  22(2)*ZZ(k,2,l,2)-10*S3h007333222(2)*ZZ(k,2,l,3)-4*S3h007132222(
     &  2)*ZZ(k,3,l,1)-4*S3h007232222(2)*ZZ(k,3,l,2)-4*S3h007332222(2)*Z
     &  Z(k,3,l,3)
       aux0073331(1,1,1)=-(F(1)*S381333111(1))-F(2)*S382333111(1)-F(3)*S
     &  383333111(1)+S3913331111(k)*Z(1,l)+S3923331111(k)*Z(2,l)+S393333
     &  1111(k)*Z(3,l)+(S3007333111(1)-S3913331111(1)-S3923332111(1)-S39
     &  33333111(1))*Z(k,l)-8*S3h007133311(1)*ZZ(k,1,l,1)-8*S3h007233311
     &  (1)*ZZ(k,1,l,2)-8*S3h007333311(1)*ZZ(k,1,l,3)-6*S3h007133111(1)*
     &  ZZ(k,3,l,1)-6*S3h007233111(1)*ZZ(k,3,l,2)-6*S3h007333111(1)*ZZ(k
     &  ,3,l,3)
       aux0073332(1,1,1)=-(F(1)*S381333211(1))-F(2)*S382333211(1)-F(3)*S
     &  383333211(1)+S3913332111(k)*Z(1,l)+S3923332111(k)*Z(2,l)+S393333
     &  2111(k)*Z(3,l)+(S3007333211(1)-S3913332111(1)-S3923332211(1)-S39
     &  33333211(1))*Z(k,l)-6*S3h007133321(1)*ZZ(k,1,l,1)-6*S3h007233321
     &  (1)*ZZ(k,1,l,2)-6*S3h007333321(1)*ZZ(k,1,l,3)-2*S3h007133311(1)*
     &  ZZ(k,2,l,1)-2*S3h007233311(1)*ZZ(k,2,l,2)-2*S3h007333311(1)*ZZ(k
     &  ,2,l,3)-6*S3h007133211(1)*ZZ(k,3,l,1)-6*S3h007233211(1)*ZZ(k,3,l
     &  ,2)-6*S3h007333211(1)*ZZ(k,3,l,3)
       aux0073332(2,1,1)=-(F(1)*S381333221(1))-F(2)*S382333221(1)-F(3)*S
     &  383333221(1)+S3913332211(k)*Z(1,l)+S3923332211(k)*Z(2,l)+S393333
     &  2211(k)*Z(3,l)+(S3007333221(1)-S3913332211(1)-S3923332221(1)-S39
     &  33333221(1))*Z(k,l)-4*S3h007133322(1)*ZZ(k,1,l,1)-4*S3h007233322
     &  (1)*ZZ(k,1,l,2)-4*S3h007333322(1)*ZZ(k,1,l,3)-4*S3h007133321(1)*
     &  ZZ(k,2,l,1)-4*S3h007233321(1)*ZZ(k,2,l,2)-4*S3h007333321(1)*ZZ(k
     &  ,2,l,3)-6*S3h007133221(1)*ZZ(k,3,l,1)-6*S3h007233221(1)*ZZ(k,3,l
     &  ,2)-6*S3h007333221(1)*ZZ(k,3,l,3)
       aux0073332(2,2,1)=-(F(1)*S381333222(1))-F(2)*S382333222(1)-F(3)*S
     &  383333222(1)+S3913332221(k)*Z(1,l)+S3923332221(k)*Z(2,l)+S393333
     &  2221(k)*Z(3,l)+(S3007333222(1)-S3913332221(1)-S3923332222(1)-S39
     &  33333222(1))*Z(k,l)-2*S3h007133322(2)*ZZ(k,1,l,1)-2*S3h007233322
     &  (2)*ZZ(k,1,l,2)-2*S3h007333322(2)*ZZ(k,1,l,3)-6*S3h007133322(1)*
     &  ZZ(k,2,l,1)-6*S3h007233322(1)*ZZ(k,2,l,2)-6*S3h007333322(1)*ZZ(k
     &  ,2,l,3)-6*S3h007133222(1)*ZZ(k,3,l,1)-6*S3h007233222(1)*ZZ(k,3,l
     &  ,2)-6*S3h007333222(1)*ZZ(k,3,l,3)
       aux0073332(2,2,2)=-(F(1)*S381333222(2))-F(2)*S382333222(2)-F(3)*S
     &  383333222(2)+S3913332222(k)*Z(1,l)+S3923332222(k)*Z(2,l)+S393333
     &  2222(k)*Z(3,l)+(S3007333222(2)-S3913332222(1)-S3923332222(2)-S39
     &  33333222(2))*Z(k,l)-8*S3h007133322(2)*ZZ(k,2,l,1)-8*S3h007233322
     &  (2)*ZZ(k,2,l,2)-8*S3h007333322(2)*ZZ(k,2,l,3)-6*S3h007133222(2)*
     &  ZZ(k,3,l,1)-6*S3h007233222(2)*ZZ(k,3,l,2)-6*S3h007333222(2)*ZZ(k
     &  ,3,l,3)
       aux0073333(1,1,1)=-(F(1)*S381333311(1))-F(2)*S382333311(1)-F(3)*S
     &  383333311(1)+S3913333111(k)*Z(1,l)+S3923333111(k)*Z(2,l)+S393333
     &  3111(k)*Z(3,l)+(S3007333311(1)-S3913333111(1)-S3923333211(1)-S39
     &  33333311(1))*Z(k,l)-6*S3h007133331(1)*ZZ(k,1,l,1)-6*S3h007233331
     &  (1)*ZZ(k,1,l,2)-6*S3h007333331(1)*ZZ(k,1,l,3)-8*S3h007133311(1)*
     &  ZZ(k,3,l,1)-8*S3h007233311(1)*ZZ(k,3,l,2)-8*S3h007333311(1)*ZZ(k
     &  ,3,l,3)
       aux0073333(2,1,1)=-(F(1)*S381333321(1))-F(2)*S382333321(1)-F(3)*S
     &  383333321(1)+S3913333211(k)*Z(1,l)+S3923333211(k)*Z(2,l)+S393333
     &  3211(k)*Z(3,l)+(S3007333321(1)-S3913333211(1)-S3923333221(1)-S39
     &  33333321(1))*Z(k,l)-4*S3h007133332(1)*ZZ(k,1,l,1)-4*S3h007233332
     &  (1)*ZZ(k,1,l,2)-4*S3h007333332(1)*ZZ(k,1,l,3)-2*S3h007133331(1)*
     &  ZZ(k,2,l,1)-2*S3h007233331(1)*ZZ(k,2,l,2)-2*S3h007333331(1)*ZZ(k
     &  ,2,l,3)-8*S3h007133321(1)*ZZ(k,3,l,1)-8*S3h007233321(1)*ZZ(k,3,l
     &  ,2)-8*S3h007333321(1)*ZZ(k,3,l,3)
       aux0073333(2,2,1)=-(F(1)*S381333322(1))-F(2)*S382333322(1)-F(3)*S
     &  383333322(1)+S3913333221(k)*Z(1,l)+S3923333221(k)*Z(2,l)+S393333
     &  3221(k)*Z(3,l)+(S3007333322(1)-S3913333221(1)-S3923333222(1)-S39
     &  33333322(1))*Z(k,l)-2*S3h007133332(2)*ZZ(k,1,l,1)-2*S3h007233332
     &  (2)*ZZ(k,1,l,2)-2*S3h007333332(2)*ZZ(k,1,l,3)-4*S3h007133332(1)*
     &  ZZ(k,2,l,1)-4*S3h007233332(1)*ZZ(k,2,l,2)-4*S3h007333332(1)*ZZ(k
     &  ,2,l,3)-8*S3h007133322(1)*ZZ(k,3,l,1)-8*S3h007233322(1)*ZZ(k,3,l
     &  ,2)-8*S3h007333322(1)*ZZ(k,3,l,3)
       aux0073333(2,2,2)=-(F(1)*S381333322(2))-F(2)*S382333322(2)-F(3)*S
     &  383333322(2)+S3913333222(k)*Z(1,l)+S3923333222(k)*Z(2,l)+S393333
     &  3222(k)*Z(3,l)+(S3007333322(2)-S3913333222(1)-S3923333222(2)-S39
     &  33333322(2))*Z(k,l)-6*S3h007133332(2)*ZZ(k,2,l,1)-6*S3h007233332
     &  (2)*ZZ(k,2,l,2)-6*S3h007333332(2)*ZZ(k,2,l,3)-8*S3h007133322(2)*
     &  ZZ(k,3,l,1)-8*S3h007233322(2)*ZZ(k,3,l,2)-8*S3h007333322(2)*ZZ(k
     &  ,3,l,3)
       aux0073333(3,1,1)=-(F(1)*S381333331(1))-F(2)*S382333331(1)-F(3)*S
     &  383333331(1)+S3913333311(k)*Z(1,l)+S3923333311(k)*Z(2,l)+S393333
     &  3311(k)*Z(3,l)+(S3007333331(1)-S3913333311(1)-S3923333321(1)-S39
     &  33333331(1))*Z(k,l)-4*S3h007133333(1)*ZZ(k,1,l,1)-4*S3h007233333
     &  (1)*ZZ(k,1,l,2)-4*S3h007333333(1)*ZZ(k,1,l,3)-10*S3h007133331(1)
     &  *ZZ(k,3,l,1)-10*S3h007233331(1)*ZZ(k,3,l,2)-10*S3h007333331(1)*Z
     &  Z(k,3,l,3)
       aux0073333(3,2,1)=-(F(1)*S381333332(1))-F(2)*S382333332(1)-F(3)*S
     &  383333332(1)+S3913333321(k)*Z(1,l)+S3923333321(k)*Z(2,l)+S393333
     &  3321(k)*Z(3,l)+(S3007333332(1)-S3913333321(1)-S3923333322(1)-S39
     &  33333332(1))*Z(k,l)-2*S3h007133333(2)*ZZ(k,1,l,1)-2*S3h007233333
     &  (2)*ZZ(k,1,l,2)-2*S3h007333333(2)*ZZ(k,1,l,3)-2*S3h007133333(1)*
     &  ZZ(k,2,l,1)-2*S3h007233333(1)*ZZ(k,2,l,2)-2*S3h007333333(1)*ZZ(k
     &  ,2,l,3)-10*S3h007133332(1)*ZZ(k,3,l,1)-10*S3h007233332(1)*ZZ(k,3
     &  ,l,2)-10*S3h007333332(1)*ZZ(k,3,l,3)
       aux0073333(3,2,2)=-(F(1)*S381333332(2))-F(2)*S382333332(2)-F(3)*S
     &  383333332(2)+S3913333322(k)*Z(1,l)+S3923333322(k)*Z(2,l)+S393333
     &  3322(k)*Z(3,l)+(S3007333332(2)-S3913333322(1)-S3923333322(2)-S39
     &  33333332(2))*Z(k,l)-4*S3h007133333(2)*ZZ(k,2,l,1)-4*S3h007233333
     &  (2)*ZZ(k,2,l,2)-4*S3h007333333(2)*ZZ(k,2,l,3)-10*S3h007133332(2)
     &  *ZZ(k,3,l,1)-10*S3h007233332(2)*ZZ(k,3,l,2)-10*S3h007333332(2)*Z
     &  Z(k,3,l,3)
       aux0073333(3,3,1)=-(F(1)*S381333333(1))-F(2)*S382333333(1)-F(3)*S
     &  383333333(1)+S3913333331(k)*Z(1,l)+S3923333331(k)*Z(2,l)+S393333
     &  3331(k)*Z(3,l)+(S3007333333(1)-S3913333331(1)-S3923333332(1)-S39
     &  33333333(1))*Z(k,l)-2*S3h007133333(3)*ZZ(k,1,l,1)-2*S3h007233333
     &  (3)*ZZ(k,1,l,2)-2*S3h007333333(3)*ZZ(k,1,l,3)-12*S3h007133333(1)
     &  *ZZ(k,3,l,1)-12*S3h007233333(1)*ZZ(k,3,l,2)-12*S3h007333333(1)*Z
     &  Z(k,3,l,3)
       aux0073333(3,3,2)=-(F(1)*S381333333(2))-F(2)*S382333333(2)-F(3)*S
     &  383333333(2)+S3913333332(k)*Z(1,l)+S3923333332(k)*Z(2,l)+S393333
     &  3332(k)*Z(3,l)+(S3007333333(2)-S3913333332(1)-S3923333332(2)-S39
     &  33333333(2))*Z(k,l)-2*S3h007133333(3)*ZZ(k,2,l,1)-2*S3h007233333
     &  (3)*ZZ(k,2,l,2)-2*S3h007333333(3)*ZZ(k,2,l,3)-12*S3h007133333(2)
     &  *ZZ(k,3,l,1)-12*S3h007233333(2)*ZZ(k,3,l,2)-12*S3h007333333(2)*Z
     &  Z(k,3,l,3)
       aux0073333(3,3,3)=-(F(1)*S381333333(3))-F(2)*S382333333(3)-F(3)*S
     &  383333333(3)+S3913333333(k)*Z(1,l)+S3923333333(k)*Z(2,l)+S393333
     &  3333(k)*Z(3,l)+(S3007333333(3)-S3913333333(1)-S3923333333(2)-S39
     &  33333333(3))*Z(k,l)-14*S3h007133333(3)*ZZ(k,3,l,1)-14*S3h0072333
     &  33(3)*ZZ(k,3,l,2)-14*S3h007333333(3)*ZZ(k,3,l,3)
       temp0071111(1,1,1)=I32Z*(aux0071111(1,1,1)+14*F(4)*temp006111(1,1
     &  ,1)+F(5)*temp71111(1,1,1)+168*temp0000511(1,1,1)*ZZ(k,1,l,1))
       temp0072111(1,1,1)=I32Z*(aux0072111(1,1,1)+2*F(6)*temp006111(1,1,
     &  1)+12*F(4)*temp006211(1,1,1)+F(5)*temp72111(1,1,1)+120*temp00005
     &  21(1,1,1)*ZZ(k,1,l,1)+24*temp0000511(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,
     &  l,1)))
       temp0072211(1,1,1)=I32Z*(aux0072211(1,1,1)+4*F(6)*temp006211(1,1,
     &  1)+10*F(4)*temp006221(1,1,1)+F(5)*temp72211(1,1,1)+80*temp000052
     &  2(1,1,1)*ZZ(k,1,l,1)+40*temp0000521(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l
     &  ,1))+8*temp0000511(1,1,1)*ZZ(k,2,l,2))
       temp0072221(1,1,1)=I32Z*(aux0072221(1,1,1)+6*F(6)*temp006221(1,1,
     &  1)+8*F(4)*temp006222(1,1,1)+F(5)*temp72221(1,1,1)+48*(temp000052
     &  2(2,1,1)*ZZ(k,1,l,1)+temp0000522(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)
     &  ))+24*temp0000521(1,1,1)*ZZ(k,2,l,2))
       temp0072222(1,1,1)=I32Z*(aux0072222(1,1,1)+8*F(6)*temp006222(1,1,
     &  1)+6*F(4)*temp006222(2,1,1)+F(5)*temp72222(1,1,1)+24*temp0000522
     &  (2,2,1)*ZZ(k,1,l,1)+48*(temp0000522(2,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l
     &  ,1))+temp0000522(1,1,1)*ZZ(k,2,l,2)))
       temp0072222(2,1,1)=I32Z*(aux0072222(2,1,1)+10*F(6)*temp006222(2,1
     &  ,1)+4*F(4)*temp006222(2,2,1)+F(5)*temp72222(2,1,1)+8*temp0000522
     &  (2,2,2)*ZZ(k,1,l,1)+40*temp0000522(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,
     &  1))+80*temp0000522(2,1,1)*ZZ(k,2,l,2))
       temp0072222(2,2,1)=I32Z*(aux0072222(2,2,1)+12*F(6)*temp006222(2,2
     &  ,1)+2*F(4)*temp006222(2,2,2)+F(5)*temp72222(2,2,1)+24*temp000052
     &  2(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+120*temp0000522(2,2,1)*ZZ(k,2
     &  ,l,2))
       temp0072222(2,2,2)=I32Z*(aux0072222(2,2,2)+14*F(6)*temp006222(2,2
     &  ,2)+F(5)*temp72222(2,2,2)+168*temp0000522(2,2,2)*ZZ(k,2,l,2))
       temp0073111(1,1,1)=I32Z*(aux0073111(1,1,1)+2*F(7)*temp006111(1,1,
     &  1)+12*F(4)*temp006311(1,1,1)+F(5)*temp73111(1,1,1)+120*temp00005
     &  31(1,1,1)*ZZ(k,1,l,1)+24*temp0000511(1,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,
     &  l,1)))
       temp0073211(1,1,1)=I32Z*(aux0073211(1,1,1)+2*F(7)*temp006211(1,1,
     &  1)+2*F(6)*temp006311(1,1,1)+10*F(4)*temp006321(1,1,1)+F(5)*temp7
     &  3211(1,1,1)+80*temp0000532(1,1,1)*ZZ(k,1,l,1)+20*(temp0000531(1,
     &  1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0000521(1,1,1)*(ZZ(k,1,l,3)+Z
     &  Z(k,3,l,1)))+4*temp0000511(1,1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp0073221(1,1,1)=I32Z*(aux0073221(1,1,1)+2*F(7)*temp006221(1,1,
     &  1)+4*F(6)*temp006321(1,1,1)+F(5)*temp73221(1,1,1)+48*temp0000532
     &  (2,1,1)*ZZ(k,1,l,1)+32*temp0000532(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,
     &  1))+8*(F(4)*temp006322(1,1,1)+temp0000531(1,1,1)*ZZ(k,2,l,2))+16
     &  *temp0000522(1,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*temp0000521(1,1,
     &  1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp0073222(1,1,1)=I32Z*(aux0073222(1,1,1)+2*F(7)*temp006222(1,1,
     &  1)+6*F(6)*temp006322(1,1,1)+6*F(4)*temp006322(2,1,1)+F(5)*temp73
     &  222(1,1,1)+36*temp0000532(2,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*(t
     &  emp0000532(2,2,1)*ZZ(k,1,l,1)+temp0000532(1,1,1)*ZZ(k,2,l,2))+12
     &  *temp0000522(2,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+12*temp0000522(1,1
     &  ,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp0073222(2,1,1)=I32Z*(aux0073222(2,1,1)+2*F(7)*temp006222(2,1,
     &  1)+8*F(6)*temp006322(2,1,1)+4*F(4)*temp006322(2,2,1)+F(5)*temp73
     &  222(2,1,1)+32*temp0000532(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*te
     &  mp0000532(2,1,1)*ZZ(k,2,l,2)+8*(temp0000532(2,2,2)*ZZ(k,1,l,1)+t
     &  emp0000522(2,2,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)))+16*temp0000522(2,1,
     &  1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp0073222(2,2,1)=I32Z*(aux0073222(2,2,1)+2*F(7)*temp006222(2,2,
     &  1)+10*F(6)*temp006322(2,2,1)+2*F(4)*temp006322(2,2,2)+F(5)*temp7
     &  3222(2,2,1)+80*temp0000532(2,2,1)*ZZ(k,2,l,2)+4*temp0000522(2,2,
     &  2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+20*(temp0000532(2,2,2)*(ZZ(k,1,l,2)
     &  +ZZ(k,2,l,1))+temp0000522(2,2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))))
       temp0073222(2,2,2)=I32Z*(aux0073222(2,2,2)+2*F(7)*temp006222(2,2,
     &  2)+12*F(6)*temp006322(2,2,2)+F(5)*temp73222(2,2,2)+120*temp00005
     &  32(2,2,2)*ZZ(k,2,l,2)+24*temp0000522(2,2,2)*(ZZ(k,2,l,3)+ZZ(k,3,
     &  l,2)))
       temp0073311(1,1,1)=I32Z*(aux0073311(1,1,1)+4*F(7)*temp006311(1,1,
     &  1)+10*F(4)*temp006331(1,1,1)+F(5)*temp73311(1,1,1)+80*temp000053
     &  3(1,1,1)*ZZ(k,1,l,1)+40*temp0000531(1,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l
     &  ,1))+8*temp0000511(1,1,1)*ZZ(k,3,l,3))
       temp0073321(1,1,1)=I32Z*(aux0073321(1,1,1)+4*F(7)*temp006321(1,1,
     &  1)+2*F(6)*temp006331(1,1,1)+F(5)*temp73321(1,1,1)+48*temp0000533
     &  (2,1,1)*ZZ(k,1,l,1)+16*temp0000533(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,
     &  1))+32*temp0000532(1,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*(F(4)*temp
     &  006332(1,1,1)+temp0000531(1,1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))+8*te
     &  mp0000521(1,1,1)*ZZ(k,3,l,3))
       temp0073322(1,1,1)=I32Z*(aux0073322(1,1,1)+4*F(7)*temp006322(1,1,
     &  1)+4*F(6)*temp006332(1,1,1)+6*F(4)*temp006332(2,1,1)+F(5)*temp73
     &  322(1,1,1)+24*(temp0000533(2,2,1)*ZZ(k,1,l,1)+temp0000533(2,1,1)
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1)))+8*temp0000533(1,1,1)*ZZ(k,2,l,2)+24*
     &  temp0000532(2,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+16*temp0000532(1,1,
     &  1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*temp0000522(1,1,1)*ZZ(k,3,l,3))
       temp0073322(2,1,1)=I32Z*(aux0073322(2,1,1)+4*F(7)*temp006322(2,1,
     &  1)+6*F(6)*temp006332(2,1,1)+4*F(4)*temp006332(2,2,1)+F(5)*temp73
     &  322(2,1,1)+24*temp0000533(2,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*te
     &  mp0000533(2,1,1)*ZZ(k,2,l,2)+16*temp0000532(2,2,1)*(ZZ(k,1,l,3)+
     &  ZZ(k,3,l,1))+24*temp0000532(2,1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*(
     &  temp0000533(2,2,2)*ZZ(k,1,l,1)+temp0000522(2,1,1)*ZZ(k,3,l,3)))
       temp0073322(2,2,1)=I32Z*(aux0073322(2,2,1)+4*F(7)*temp006322(2,2,
     &  1)+2*F(4)*temp006332(2,2,2)+F(5)*temp73322(2,2,1)+16*temp0000533
     &  (2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*temp0000533(2,2,1)*ZZ(k,2,l
     &  ,2)+8*(F(6)*temp006332(2,2,1)+temp0000532(2,2,2)*(ZZ(k,1,l,3)+ZZ
     &  (k,3,l,1)))+32*temp0000532(2,2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*te
     &  mp0000522(2,2,1)*ZZ(k,3,l,3))
       temp0073322(2,2,2)=I32Z*(aux0073322(2,2,2)+4*F(7)*temp006322(2,2,
     &  2)+10*F(6)*temp006332(2,2,2)+F(5)*temp73322(2,2,2)+80*temp000053
     &  3(2,2,2)*ZZ(k,2,l,2)+40*temp0000532(2,2,2)*(ZZ(k,2,l,3)+ZZ(k,3,l
     &  ,2))+8*temp0000522(2,2,2)*ZZ(k,3,l,3))
       temp0073331(1,1,1)=I32Z*(aux0073331(1,1,1)+6*F(7)*temp006331(1,1,
     &  1)+8*F(4)*temp006333(1,1,1)+F(5)*temp73331(1,1,1)+48*(temp000053
     &  3(3,1,1)*ZZ(k,1,l,1)+temp0000533(1,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)
     &  ))+24*temp0000531(1,1,1)*ZZ(k,3,l,3))
       temp0073332(1,1,1)=I32Z*(aux0073332(1,1,1)+6*F(7)*temp006332(1,1,
     &  1)+2*F(6)*temp006333(1,1,1)+6*F(4)*temp006333(2,1,1)+F(5)*temp73
     &  332(1,1,1)+12*temp0000533(3,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+36*te
     &  mp0000533(2,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+12*temp0000533(1,1,1)
     &  *(ZZ(k,2,l,3)+ZZ(k,3,l,2))+24*(temp0000533(3,2,1)*ZZ(k,1,l,1)+te
     &  mp0000532(1,1,1)*ZZ(k,3,l,3)))
       temp0073332(2,1,1)=I32Z*(aux0073332(2,1,1)+6*F(7)*temp006332(2,1,
     &  1)+4*F(6)*temp006333(2,1,1)+4*F(4)*temp006333(2,2,1)+F(5)*temp73
     &  332(2,1,1)+16*temp0000533(3,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(te
     &  mp0000533(3,2,2)*ZZ(k,1,l,1)+temp0000533(3,1,1)*ZZ(k,2,l,2))+24*
     &  temp0000533(2,2,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+24*temp0000533(2,1,
     &  1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+24*temp0000532(2,1,1)*ZZ(k,3,l,3))
       temp0073332(2,2,1)=I32Z*(aux0073332(2,2,1)+6*F(7)*temp006332(2,2,
     &  1)+6*F(6)*temp006333(2,2,1)+2*F(4)*temp006333(2,2,2)+F(5)*temp73
     &  332(2,2,1)+24*temp0000533(3,2,1)*ZZ(k,2,l,2)+12*(temp0000533(3,2
     &  ,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0000533(2,2,2)*(ZZ(k,1,l,3)+ZZ
     &  (k,3,l,1)))+36*temp0000533(2,2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+24*t
     &  emp0000532(2,2,1)*ZZ(k,3,l,3))
       temp0073332(2,2,2)=I32Z*(aux0073332(2,2,2)+6*F(7)*temp006332(2,2,
     &  2)+8*F(6)*temp006333(2,2,2)+F(5)*temp73332(2,2,2)+48*(temp000053
     &  3(3,2,2)*ZZ(k,2,l,2)+temp0000533(2,2,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)
     &  ))+24*temp0000532(2,2,2)*ZZ(k,3,l,3))
       temp0073333(1,1,1)=I32Z*(aux0073333(1,1,1)+8*F(7)*temp006333(1,1,
     &  1)+6*F(4)*temp006333(3,1,1)+F(5)*temp73333(1,1,1)+24*temp0000533
     &  (3,3,1)*ZZ(k,1,l,1)+48*(temp0000533(3,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l
     &  ,1))+temp0000533(1,1,1)*ZZ(k,3,l,3)))
       temp0073333(2,1,1)=I32Z*(aux0073333(2,1,1)+8*F(7)*temp006333(2,1,
     &  1)+2*F(6)*temp006333(3,1,1)+4*F(4)*temp006333(3,2,1)+F(5)*temp73
     &  333(2,1,1)+8*(temp0000533(3,3,2)*ZZ(k,1,l,1)+temp0000533(3,3,1)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1)))+32*temp0000533(3,2,1)*(ZZ(k,1,l,3)+ZZ
     &  (k,3,l,1))+16*temp0000533(3,1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+48*te
     &  mp0000533(2,1,1)*ZZ(k,3,l,3))
       temp0073333(2,2,1)=I32Z*(aux0073333(2,2,1)+8*F(7)*temp006333(2,2,
     &  1)+4*F(6)*temp006333(3,2,1)+2*F(4)*temp006333(3,2,2)+F(5)*temp73
     &  333(2,2,1)+8*(temp0000533(3,3,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0
     &  000533(3,3,1)*ZZ(k,2,l,2))+16*temp0000533(3,2,2)*(ZZ(k,1,l,3)+ZZ
     &  (k,3,l,1))+32*temp0000533(3,2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+48*te
     &  mp0000533(2,2,1)*ZZ(k,3,l,3))
       temp0073333(2,2,2)=I32Z*(aux0073333(2,2,2)+8*F(7)*temp006333(2,2,
     &  2)+6*F(6)*temp006333(3,2,2)+F(5)*temp73333(2,2,2)+24*temp0000533
     &  (3,3,2)*ZZ(k,2,l,2)+48*(temp0000533(3,2,2)*(ZZ(k,2,l,3)+ZZ(k,3,l
     &  ,2))+temp0000533(2,2,2)*ZZ(k,3,l,3)))
       temp0073333(3,1,1)=I32Z*(aux0073333(3,1,1)+10*F(7)*temp006333(3,1
     &  ,1)+4*F(4)*temp006333(3,3,1)+F(5)*temp73333(3,1,1)+8*temp0000533
     &  (3,3,3)*ZZ(k,1,l,1)+40*temp0000533(3,3,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,
     &  1))+80*temp0000533(3,1,1)*ZZ(k,3,l,3))
       temp0073333(3,2,1)=I32Z*(aux0073333(3,2,1)+10*F(7)*temp006333(3,2
     &  ,1)+2*F(6)*temp006333(3,3,1)+2*F(4)*temp006333(3,3,2)+F(5)*temp7
     &  3333(3,2,1)+4*temp0000533(3,3,3)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+20*(t
     &  emp0000533(3,3,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+temp0000533(3,3,1)*(
     &  ZZ(k,2,l,3)+ZZ(k,3,l,2)))+80*temp0000533(3,2,1)*ZZ(k,3,l,3))
       temp0073333(3,2,2)=I32Z*(aux0073333(3,2,2)+10*F(7)*temp006333(3,2
     &  ,2)+4*F(6)*temp006333(3,3,2)+F(5)*temp73333(3,2,2)+8*temp0000533
     &  (3,3,3)*ZZ(k,2,l,2)+40*temp0000533(3,3,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,
     &  2))+80*temp0000533(3,2,2)*ZZ(k,3,l,3))
       temp0073333(3,3,1)=I32Z*(aux0073333(3,3,1)+12*F(7)*temp006333(3,3
     &  ,1)+2*F(4)*temp006333(3,3,3)+F(5)*temp73333(3,3,1)+24*temp000053
     &  3(3,3,3)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+120*temp0000533(3,3,1)*ZZ(k,3
     &  ,l,3))
       temp0073333(3,3,2)=I32Z*(aux0073333(3,3,2)+12*F(7)*temp006333(3,3
     &  ,2)+2*F(6)*temp006333(3,3,3)+F(5)*temp73333(3,3,2)+24*temp000053
     &  3(3,3,3)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+120*temp0000533(3,3,2)*ZZ(k,3
     &  ,l,3))
       temp0073333(3,3,3)=I32Z*(aux0073333(3,3,3)+14*F(7)*temp006333(3,3
     &  ,3)+F(5)*temp73333(3,3,3)+168*temp0000533(3,3,3)*ZZ(k,3,l,3))
       temp0071111(1,1,2)=temp0072111(1,1,1)
       temp0071111(1,1,3)=temp0073111(1,1,1)
       temp0071111(1,2,1)=temp0072111(1,1,1)
       temp0071111(1,2,2)=temp0072211(1,1,1)
       temp0071111(1,2,3)=temp0073211(1,1,1)
       temp0071111(1,3,1)=temp0073111(1,1,1)
       temp0071111(1,3,2)=temp0073211(1,1,1)
       temp0071111(1,3,3)=temp0073311(1,1,1)
       temp0072111(1,1,2)=temp0072211(1,1,1)
       temp0072111(1,1,3)=temp0073211(1,1,1)
       temp0072111(1,2,1)=temp0072211(1,1,1)
       temp0072111(1,2,2)=temp0072221(1,1,1)
       temp0072111(1,2,3)=temp0073221(1,1,1)
       temp0072111(1,3,1)=temp0073211(1,1,1)
       temp0072111(1,3,2)=temp0073221(1,1,1)
       temp0072111(1,3,3)=temp0073321(1,1,1)
       temp0072211(1,1,2)=temp0072221(1,1,1)
       temp0072211(1,1,3)=temp0073221(1,1,1)
       temp0072211(1,2,1)=temp0072221(1,1,1)
       temp0072211(1,2,2)=temp0072222(1,1,1)
       temp0072211(1,2,3)=temp0073222(1,1,1)
       temp0072211(1,3,1)=temp0073221(1,1,1)
       temp0072211(1,3,2)=temp0073222(1,1,1)
       temp0072211(1,3,3)=temp0073322(1,1,1)
       temp0072221(1,1,2)=temp0072222(1,1,1)
       temp0072221(1,1,3)=temp0073222(1,1,1)
       temp0072221(1,2,1)=temp0072222(1,1,1)
       temp0072221(1,2,2)=temp0072222(2,1,1)
       temp0072221(1,2,3)=temp0073222(2,1,1)
       temp0072221(1,3,1)=temp0073222(1,1,1)
       temp0072221(1,3,2)=temp0073222(2,1,1)
       temp0072221(1,3,3)=temp0073322(2,1,1)
       temp0072222(1,1,2)=temp0072222(2,1,1)
       temp0072222(1,1,3)=temp0073222(2,1,1)
       temp0072222(1,2,1)=temp0072222(2,1,1)
       temp0072222(1,2,2)=temp0072222(2,2,1)
       temp0072222(1,2,3)=temp0073222(2,2,1)
       temp0072222(1,3,1)=temp0073222(2,1,1)
       temp0072222(1,3,2)=temp0073222(2,2,1)
       temp0072222(1,3,3)=temp0073322(2,2,1)
       temp0072222(2,1,2)=temp0072222(2,2,1)
       temp0072222(2,1,3)=temp0073222(2,2,1)
       temp0072222(2,2,3)=temp0073222(2,2,2)
       temp0072222(2,3,1)=temp0073222(2,2,1)
       temp0072222(2,3,2)=temp0073222(2,2,2)
       temp0072222(2,3,3)=temp0073322(2,2,2)
       temp0073111(1,1,2)=temp0073211(1,1,1)
       temp0073111(1,1,3)=temp0073311(1,1,1)
       temp0073111(1,2,1)=temp0073211(1,1,1)
       temp0073111(1,2,2)=temp0073221(1,1,1)
       temp0073111(1,2,3)=temp0073321(1,1,1)
       temp0073111(1,3,1)=temp0073311(1,1,1)
       temp0073111(1,3,2)=temp0073321(1,1,1)
       temp0073111(1,3,3)=temp0073331(1,1,1)
       temp0073211(1,1,2)=temp0073221(1,1,1)
       temp0073211(1,1,3)=temp0073321(1,1,1)
       temp0073211(1,2,1)=temp0073221(1,1,1)
       temp0073211(1,2,2)=temp0073222(1,1,1)
       temp0073211(1,2,3)=temp0073322(1,1,1)
       temp0073211(1,3,1)=temp0073321(1,1,1)
       temp0073211(1,3,2)=temp0073322(1,1,1)
       temp0073211(1,3,3)=temp0073332(1,1,1)
       temp0073221(1,1,2)=temp0073222(1,1,1)
       temp0073221(1,1,3)=temp0073322(1,1,1)
       temp0073221(1,2,1)=temp0073222(1,1,1)
       temp0073221(1,2,2)=temp0073222(2,1,1)
       temp0073221(1,2,3)=temp0073322(2,1,1)
       temp0073221(1,3,1)=temp0073322(1,1,1)
       temp0073221(1,3,2)=temp0073322(2,1,1)
       temp0073221(1,3,3)=temp0073332(2,1,1)
       temp0073222(1,1,2)=temp0073222(2,1,1)
       temp0073222(1,1,3)=temp0073322(2,1,1)
       temp0073222(1,2,1)=temp0073222(2,1,1)
       temp0073222(1,2,2)=temp0073222(2,2,1)
       temp0073222(1,2,3)=temp0073322(2,2,1)
       temp0073222(1,3,1)=temp0073322(2,1,1)
       temp0073222(1,3,2)=temp0073322(2,2,1)
       temp0073222(1,3,3)=temp0073332(2,2,1)
       temp0073222(2,1,2)=temp0073222(2,2,1)
       temp0073222(2,1,3)=temp0073322(2,2,1)
       temp0073222(2,2,3)=temp0073322(2,2,2)
       temp0073222(2,3,1)=temp0073322(2,2,1)
       temp0073222(2,3,2)=temp0073322(2,2,2)
       temp0073222(2,3,3)=temp0073332(2,2,2)
       temp0073311(1,1,2)=temp0073321(1,1,1)
       temp0073311(1,1,3)=temp0073331(1,1,1)
       temp0073311(1,2,1)=temp0073321(1,1,1)
       temp0073311(1,2,2)=temp0073322(1,1,1)
       temp0073311(1,2,3)=temp0073332(1,1,1)
       temp0073311(1,3,1)=temp0073331(1,1,1)
       temp0073311(1,3,2)=temp0073332(1,1,1)
       temp0073311(1,3,3)=temp0073333(1,1,1)
       temp0073321(1,1,2)=temp0073322(1,1,1)
       temp0073321(1,1,3)=temp0073332(1,1,1)
       temp0073321(1,2,1)=temp0073322(1,1,1)
       temp0073321(1,2,2)=temp0073322(2,1,1)
       temp0073321(1,2,3)=temp0073332(2,1,1)
       temp0073321(1,3,1)=temp0073332(1,1,1)
       temp0073321(1,3,2)=temp0073332(2,1,1)
       temp0073321(1,3,3)=temp0073333(2,1,1)
       temp0073322(1,1,2)=temp0073322(2,1,1)
       temp0073322(1,1,3)=temp0073332(2,1,1)
       temp0073322(1,2,1)=temp0073322(2,1,1)
       temp0073322(1,2,2)=temp0073322(2,2,1)
       temp0073322(1,2,3)=temp0073332(2,2,1)
       temp0073322(1,3,1)=temp0073332(2,1,1)
       temp0073322(1,3,2)=temp0073332(2,2,1)
       temp0073322(1,3,3)=temp0073333(2,2,1)
       temp0073322(2,1,2)=temp0073322(2,2,1)
       temp0073322(2,1,3)=temp0073332(2,2,1)
       temp0073322(2,2,3)=temp0073332(2,2,2)
       temp0073322(2,3,1)=temp0073332(2,2,1)
       temp0073322(2,3,2)=temp0073332(2,2,2)
       temp0073322(2,3,3)=temp0073333(2,2,2)
       temp0073331(1,1,2)=temp0073332(1,1,1)
       temp0073331(1,1,3)=temp0073333(1,1,1)
       temp0073331(1,2,1)=temp0073332(1,1,1)
       temp0073331(1,2,2)=temp0073332(2,1,1)
       temp0073331(1,2,3)=temp0073333(2,1,1)
       temp0073331(1,3,1)=temp0073333(1,1,1)
       temp0073331(1,3,2)=temp0073333(2,1,1)
       temp0073331(1,3,3)=temp0073333(3,1,1)
       temp0073332(1,1,2)=temp0073332(2,1,1)
       temp0073332(1,1,3)=temp0073333(2,1,1)
       temp0073332(1,2,1)=temp0073332(2,1,1)
       temp0073332(1,2,2)=temp0073332(2,2,1)
       temp0073332(1,2,3)=temp0073333(2,2,1)
       temp0073332(1,3,1)=temp0073333(2,1,1)
       temp0073332(1,3,2)=temp0073333(2,2,1)
       temp0073332(1,3,3)=temp0073333(3,2,1)
       temp0073332(2,1,2)=temp0073332(2,2,1)
       temp0073332(2,1,3)=temp0073333(2,2,1)
       temp0073332(2,2,3)=temp0073333(2,2,2)
       temp0073332(2,3,1)=temp0073333(2,2,1)
       temp0073332(2,3,2)=temp0073333(2,2,2)
       temp0073332(2,3,3)=temp0073333(3,2,2)
       temp0073333(1,1,2)=temp0073333(2,1,1)
       temp0073333(1,1,3)=temp0073333(3,1,1)
       temp0073333(1,2,1)=temp0073333(2,1,1)
       temp0073333(1,2,2)=temp0073333(2,2,1)
       temp0073333(1,2,3)=temp0073333(3,2,1)
       temp0073333(1,3,1)=temp0073333(3,1,1)
       temp0073333(1,3,2)=temp0073333(3,2,1)
       temp0073333(1,3,3)=temp0073333(3,3,1)
       temp0073333(2,1,2)=temp0073333(2,2,1)
       temp0073333(2,1,3)=temp0073333(3,2,1)
       temp0073333(2,2,3)=temp0073333(3,2,2)
       temp0073333(2,3,1)=temp0073333(3,2,1)
       temp0073333(2,3,2)=temp0073333(3,2,2)
       temp0073333(2,3,3)=temp0073333(3,3,2)
       temp0073333(3,1,2)=temp0073333(3,2,1)
       temp0073333(3,1,3)=temp0073333(3,3,1)
       temp0073333(3,2,3)=temp0073333(3,3,2)
       aux811111(1,1,1)=-(S3911111111(1)*Z(jj,1))-S3921111111(1)*Z(jj,2)
     &  -S3931111111(1)*Z(jj,3)
       aux821111(1,1,1)=-(S3912111111(1)*Z(jj,1))-S3922111111(1)*Z(jj,2)
     &  -S3932111111(1)*Z(jj,3)
       aux822111(1,1,1)=-(S3912211111(1)*Z(jj,1))-S3922211111(1)*Z(jj,2)
     &  -S3932211111(1)*Z(jj,3)
       aux822211(1,1,1)=-(S3912221111(1)*Z(jj,1))-S3922221111(1)*Z(jj,2)
     &  -S3932221111(1)*Z(jj,3)
       aux822221(1,1,1)=-(S3912222111(1)*Z(jj,1))-S3922222111(1)*Z(jj,2)
     &  -S3932222111(1)*Z(jj,3)
       aux822221(2,1,1)=-(S3912222211(1)*Z(jj,1))-S3922222211(1)*Z(jj,2)
     &  -S3932222211(1)*Z(jj,3)
       aux822221(2,2,1)=-(S3912222221(1)*Z(jj,1))-S3922222221(1)*Z(jj,2)
     &  -S3932222221(1)*Z(jj,3)
       aux822221(2,2,2)=-(S3912222222(1)*Z(jj,1))-S3922222222(1)*Z(jj,2)
     &  -S3932222222(1)*Z(jj,3)
       aux822222(1,1,1)=-(S3912222211(1)*Z(jj,1))-S3922222211(1)*Z(jj,2)
     &  -S3932222211(1)*Z(jj,3)
       aux822222(2,1,1)=-(S3912222221(1)*Z(jj,1))-S3922222221(1)*Z(jj,2)
     &  -S3932222221(1)*Z(jj,3)
       aux822222(2,2,1)=-(S3912222222(1)*Z(jj,1))-S3922222222(1)*Z(jj,2)
     &  -S3932222222(1)*Z(jj,3)
       aux822222(2,2,2)=-(S3912222222(2)*Z(jj,1))-S3922222222(2)*Z(jj,2)
     &  -S3932222222(2)*Z(jj,3)
       aux831111(1,1,1)=-(S3913111111(1)*Z(jj,1))-S3923111111(1)*Z(jj,2)
     &  -S3933111111(1)*Z(jj,3)
       aux832111(1,1,1)=-(S3913211111(1)*Z(jj,1))-S3923211111(1)*Z(jj,2)
     &  -S3933211111(1)*Z(jj,3)
       aux832211(1,1,1)=-(S3913221111(1)*Z(jj,1))-S3923221111(1)*Z(jj,2)
     &  -S3933221111(1)*Z(jj,3)
       aux832221(1,1,1)=-(S3913222111(1)*Z(jj,1))-S3923222111(1)*Z(jj,2)
     &  -S3933222111(1)*Z(jj,3)
       aux832221(2,1,1)=-(S3913222211(1)*Z(jj,1))-S3923222211(1)*Z(jj,2)
     &  -S3933222211(1)*Z(jj,3)
       aux832221(2,2,1)=-(S3913222221(1)*Z(jj,1))-S3923222221(1)*Z(jj,2)
     &  -S3933222221(1)*Z(jj,3)
       aux832221(2,2,2)=-(S3913222222(1)*Z(jj,1))-S3923222222(1)*Z(jj,2)
     &  -S3933222222(1)*Z(jj,3)
       aux832222(1,1,1)=-(S3913222211(1)*Z(jj,1))-S3923222211(1)*Z(jj,2)
     &  -S3933222211(1)*Z(jj,3)
       aux832222(2,1,1)=-(S3913222221(1)*Z(jj,1))-S3923222221(1)*Z(jj,2)
     &  -S3933222221(1)*Z(jj,3)
       aux832222(2,2,1)=-(S3913222222(1)*Z(jj,1))-S3923222222(1)*Z(jj,2)
     &  -S3933222222(1)*Z(jj,3)
       aux832222(2,2,2)=-(S3913222222(2)*Z(jj,1))-S3923222222(2)*Z(jj,2)
     &  -S3933222222(2)*Z(jj,3)
       aux833111(1,1,1)=-(S3913311111(1)*Z(jj,1))-S3923311111(1)*Z(jj,2)
     &  -S3933311111(1)*Z(jj,3)
       aux833211(1,1,1)=-(S3913321111(1)*Z(jj,1))-S3923321111(1)*Z(jj,2)
     &  -S3933321111(1)*Z(jj,3)
       aux833221(1,1,1)=-(S3913322111(1)*Z(jj,1))-S3923322111(1)*Z(jj,2)
     &  -S3933322111(1)*Z(jj,3)
       aux833221(2,1,1)=-(S3913322211(1)*Z(jj,1))-S3923322211(1)*Z(jj,2)
     &  -S3933322211(1)*Z(jj,3)
       aux833221(2,2,1)=-(S3913322221(1)*Z(jj,1))-S3923322221(1)*Z(jj,2)
     &  -S3933322221(1)*Z(jj,3)
       aux833221(2,2,2)=-(S3913322222(1)*Z(jj,1))-S3923322222(1)*Z(jj,2)
     &  -S3933322222(1)*Z(jj,3)
       aux833222(1,1,1)=-(S3913322211(1)*Z(jj,1))-S3923322211(1)*Z(jj,2)
     &  -S3933322211(1)*Z(jj,3)
       aux833222(2,1,1)=-(S3913322221(1)*Z(jj,1))-S3923322221(1)*Z(jj,2)
     &  -S3933322221(1)*Z(jj,3)
       aux833222(2,2,1)=-(S3913322222(1)*Z(jj,1))-S3923322222(1)*Z(jj,2)
     &  -S3933322222(1)*Z(jj,3)
       aux833222(2,2,2)=-(S3913322222(2)*Z(jj,1))-S3923322222(2)*Z(jj,2)
     &  -S3933322222(2)*Z(jj,3)
       aux833311(1,1,1)=-(S3913331111(1)*Z(jj,1))-S3923331111(1)*Z(jj,2)
     &  -S3933331111(1)*Z(jj,3)
       aux833321(1,1,1)=-(S3913332111(1)*Z(jj,1))-S3923332111(1)*Z(jj,2)
     &  -S3933332111(1)*Z(jj,3)
       aux833321(2,1,1)=-(S3913332211(1)*Z(jj,1))-S3923332211(1)*Z(jj,2)
     &  -S3933332211(1)*Z(jj,3)
       aux833321(2,2,1)=-(S3913332221(1)*Z(jj,1))-S3923332221(1)*Z(jj,2)
     &  -S3933332221(1)*Z(jj,3)
       aux833321(2,2,2)=-(S3913332222(1)*Z(jj,1))-S3923332222(1)*Z(jj,2)
     &  -S3933332222(1)*Z(jj,3)
       aux833322(1,1,1)=-(S3913332211(1)*Z(jj,1))-S3923332211(1)*Z(jj,2)
     &  -S3933332211(1)*Z(jj,3)
       aux833322(2,1,1)=-(S3913332221(1)*Z(jj,1))-S3923332221(1)*Z(jj,2)
     &  -S3933332221(1)*Z(jj,3)
       aux833322(2,2,1)=-(S3913332222(1)*Z(jj,1))-S3923332222(1)*Z(jj,2)
     &  -S3933332222(1)*Z(jj,3)
       aux833322(2,2,2)=-(S3913332222(2)*Z(jj,1))-S3923332222(2)*Z(jj,2)
     &  -S3933332222(2)*Z(jj,3)
       aux833331(1,1,1)=-(S3913333111(1)*Z(jj,1))-S3923333111(1)*Z(jj,2)
     &  -S3933333111(1)*Z(jj,3)
       aux833331(2,1,1)=-(S3913333211(1)*Z(jj,1))-S3923333211(1)*Z(jj,2)
     &  -S3933333211(1)*Z(jj,3)
       aux833331(2,2,1)=-(S3913333221(1)*Z(jj,1))-S3923333221(1)*Z(jj,2)
     &  -S3933333221(1)*Z(jj,3)
       aux833331(2,2,2)=-(S3913333222(1)*Z(jj,1))-S3923333222(1)*Z(jj,2)
     &  -S3933333222(1)*Z(jj,3)
       aux833331(3,1,1)=-(S3913333311(1)*Z(jj,1))-S3923333311(1)*Z(jj,2)
     &  -S3933333311(1)*Z(jj,3)
       aux833331(3,2,1)=-(S3913333321(1)*Z(jj,1))-S3923333321(1)*Z(jj,2)
     &  -S3933333321(1)*Z(jj,3)
       aux833331(3,2,2)=-(S3913333322(1)*Z(jj,1))-S3923333322(1)*Z(jj,2)
     &  -S3933333322(1)*Z(jj,3)
       aux833331(3,3,1)=-(S3913333331(1)*Z(jj,1))-S3923333331(1)*Z(jj,2)
     &  -S3933333331(1)*Z(jj,3)
       aux833331(3,3,2)=-(S3913333332(1)*Z(jj,1))-S3923333332(1)*Z(jj,2)
     &  -S3933333332(1)*Z(jj,3)
       aux833331(3,3,3)=-(S3913333333(1)*Z(jj,1))-S3923333333(1)*Z(jj,2)
     &  -S3933333333(1)*Z(jj,3)
       aux833332(1,1,1)=-(S3913333211(1)*Z(jj,1))-S3923333211(1)*Z(jj,2)
     &  -S3933333211(1)*Z(jj,3)
       aux833332(2,1,1)=-(S3913333221(1)*Z(jj,1))-S3923333221(1)*Z(jj,2)
     &  -S3933333221(1)*Z(jj,3)
       aux833332(2,2,1)=-(S3913333222(1)*Z(jj,1))-S3923333222(1)*Z(jj,2)
     &  -S3933333222(1)*Z(jj,3)
       aux833332(2,2,2)=-(S3913333222(2)*Z(jj,1))-S3923333222(2)*Z(jj,2)
     &  -S3933333222(2)*Z(jj,3)
       aux833332(3,1,1)=-(S3913333321(1)*Z(jj,1))-S3923333321(1)*Z(jj,2)
     &  -S3933333321(1)*Z(jj,3)
       aux833332(3,2,1)=-(S3913333322(1)*Z(jj,1))-S3923333322(1)*Z(jj,2)
     &  -S3933333322(1)*Z(jj,3)
       aux833332(3,2,2)=-(S3913333322(2)*Z(jj,1))-S3923333322(2)*Z(jj,2)
     &  -S3933333322(2)*Z(jj,3)
       aux833332(3,3,1)=-(S3913333332(1)*Z(jj,1))-S3923333332(1)*Z(jj,2)
     &  -S3933333332(1)*Z(jj,3)
       aux833332(3,3,2)=-(S3913333332(2)*Z(jj,1))-S3923333332(2)*Z(jj,2)
     &  -S3933333332(2)*Z(jj,3)
       aux833332(3,3,3)=-(S3913333333(2)*Z(jj,1))-S3923333333(2)*Z(jj,2)
     &  -S3933333333(2)*Z(jj,3)
       aux833333(1,1,1)=-(S3913333311(1)*Z(jj,1))-S3923333311(1)*Z(jj,2)
     &  -S3933333311(1)*Z(jj,3)
       aux833333(2,1,1)=-(S3913333321(1)*Z(jj,1))-S3923333321(1)*Z(jj,2)
     &  -S3933333321(1)*Z(jj,3)
       aux833333(2,2,1)=-(S3913333322(1)*Z(jj,1))-S3923333322(1)*Z(jj,2)
     &  -S3933333322(1)*Z(jj,3)
       aux833333(2,2,2)=-(S3913333322(2)*Z(jj,1))-S3923333322(2)*Z(jj,2)
     &  -S3933333322(2)*Z(jj,3)
       aux833333(3,1,1)=-(S3913333331(1)*Z(jj,1))-S3923333331(1)*Z(jj,2)
     &  -S3933333331(1)*Z(jj,3)
       aux833333(3,2,1)=-(S3913333332(1)*Z(jj,1))-S3923333332(1)*Z(jj,2)
     &  -S3933333332(1)*Z(jj,3)
       aux833333(3,2,2)=-(S3913333332(2)*Z(jj,1))-S3923333332(2)*Z(jj,2)
     &  -S3933333332(2)*Z(jj,3)
       aux833333(3,3,1)=-(S3913333333(1)*Z(jj,1))-S3923333333(1)*Z(jj,2)
     &  -S3933333333(1)*Z(jj,3)
       aux833333(3,3,2)=-(S3913333333(2)*Z(jj,1))-S3923333333(2)*Z(jj,2)
     &  -S3933333333(2)*Z(jj,3)
       aux833333(3,3,3)=-(S3913333333(3)*Z(jj,1))-S3923333333(3)*Z(jj,2)
     &  -S3933333333(3)*Z(jj,3)
       temp811111(1,1,1)=IX*(aux811111(1,1,1)+16*temp0071111(1,1,1)*Z(jj
     &  ,1))
       temp821111(1,1,1)=IX*(aux821111(1,1,1)+14*temp0072111(1,1,1)*Z(jj
     &  ,1)+2*temp0071111(1,1,1)*Z(jj,2))
       temp822111(1,1,1)=IX*(aux822111(1,1,1)+12*temp0072211(1,1,1)*Z(jj
     &  ,1)+4*temp0072111(1,1,1)*Z(jj,2))
       temp822211(1,1,1)=IX*(aux822211(1,1,1)+10*temp0072221(1,1,1)*Z(jj
     &  ,1)+6*temp0072211(1,1,1)*Z(jj,2))
       temp822221(1,1,1)=IX*(aux822221(1,1,1)+8*(temp0072222(1,1,1)*Z(jj
     &  ,1)+temp0072221(1,1,1)*Z(jj,2)))
       temp822221(2,1,1)=IX*(aux822221(2,1,1)+4*temp0072222(1,2,1)*Z(jj,
     &  1)+8*temp0072221(2,1,1)*Z(jj,2)+2*(temp0072222(2,1,1)*Z(jj,1)+te
     &  mp0072222(1,1,1)*Z(jj,2)))
       temp822221(2,2,1)=IX*(aux822221(2,2,1)+2*(temp0072222(1,2,2)*Z(jj
     &  ,1)+temp0072222(2,2,1)*Z(jj,1))+8*temp0072221(2,2,1)*Z(jj,2)+4*t
     &  emp0072222(1,2,1)*Z(jj,2))
       temp822221(2,2,2)=IX*(aux822221(2,2,2)+2*temp0072222(2,2,2)*Z(jj,
     &  1)+8*temp0072221(2,2,2)*Z(jj,2)+6*temp0072222(1,2,2)*Z(jj,2))
       temp822222(1,1,1)=IX*(aux822222(1,1,1)+6*temp0072222(2,1,1)*Z(jj,
     &  1)+10*temp0072222(1,1,1)*Z(jj,2))
       temp822222(2,1,1)=IX*(aux822222(2,1,1)+4*temp0072222(2,2,1)*Z(jj,
     &  1)+12*temp0072222(2,1,1)*Z(jj,2))
       temp822222(2,2,1)=IX*(aux822222(2,2,1)+2*temp0072222(2,2,2)*Z(jj,
     &  1)+14*temp0072222(2,2,1)*Z(jj,2))
       temp822222(2,2,2)=IX*(aux822222(2,2,2)+16*temp0072222(2,2,2)*Z(jj
     &  ,2))
       temp831111(1,1,1)=IX*(aux831111(1,1,1)+14*temp0073111(1,1,1)*Z(jj
     &  ,1)+2*temp0071111(1,1,1)*Z(jj,3))
       temp832111(1,1,1)=IX*(aux832111(1,1,1)+12*temp0073211(1,1,1)*Z(jj
     &  ,1)+2*(temp0073111(1,1,1)*Z(jj,2)+temp0072111(1,1,1)*Z(jj,3)))
       temp832211(1,1,1)=IX*(aux832211(1,1,1)+10*temp0073221(1,1,1)*Z(jj
     &  ,1)+4*temp0073211(1,1,1)*Z(jj,2)+2*temp0072211(1,1,1)*Z(jj,3))
       temp832221(1,1,1)=IX*(aux832221(1,1,1)+8*temp0073222(1,1,1)*Z(jj,
     &  1)+6*temp0073221(1,1,1)*Z(jj,2)+2*temp0072221(1,1,1)*Z(jj,3))
       temp832221(2,1,1)=IX*(aux832221(2,1,1)+4*temp0073222(1,2,1)*Z(jj,
     &  1)+6*temp0073221(2,1,1)*Z(jj,2)+2*(temp0073222(2,1,1)*Z(jj,1)+te
     &  mp0073222(1,1,1)*Z(jj,2))+2*temp0072221(2,1,1)*Z(jj,3))
       temp832221(2,2,1)=IX*(aux832221(2,2,1)+2*(temp0073222(1,2,2)*Z(jj
     &  ,1)+temp0073222(2,2,1)*Z(jj,1))+6*temp0073221(2,2,1)*Z(jj,2)+4*t
     &  emp0073222(1,2,1)*Z(jj,2)+2*temp0072221(2,2,1)*Z(jj,3))
       temp832221(2,2,2)=IX*(aux832221(2,2,2)+6*temp0073221(2,2,2)*Z(jj,
     &  2)+6*temp0073222(1,2,2)*Z(jj,2)+2*(temp0073222(2,2,2)*Z(jj,1)+te
     &  mp0072221(2,2,2)*Z(jj,3)))
       temp832222(1,1,1)=IX*(aux832222(1,1,1)+6*temp0073222(2,1,1)*Z(jj,
     &  1)+8*temp0073222(1,1,1)*Z(jj,2)+2*temp0072222(1,1,1)*Z(jj,3))
       temp832222(2,1,1)=IX*(aux832222(2,1,1)+4*temp0073222(2,2,1)*Z(jj,
     &  1)+10*temp0073222(2,1,1)*Z(jj,2)+2*temp0072222(2,1,1)*Z(jj,3))
       temp832222(2,2,1)=IX*(aux832222(2,2,1)+12*temp0073222(2,2,1)*Z(jj
     &  ,2)+2*(temp0073222(2,2,2)*Z(jj,1)+temp0072222(2,2,1)*Z(jj,3)))
       temp832222(2,2,2)=IX*(aux832222(2,2,2)+14*temp0073222(2,2,2)*Z(jj
     &  ,2)+2*temp0072222(2,2,2)*Z(jj,3))
       temp833111(1,1,1)=IX*(aux833111(1,1,1)+12*temp0073311(1,1,1)*Z(jj
     &  ,1)+4*temp0073111(1,1,1)*Z(jj,3))
       temp833211(1,1,1)=IX*(aux833211(1,1,1)+10*temp0073321(1,1,1)*Z(jj
     &  ,1)+2*temp0073311(1,1,1)*Z(jj,2)+4*temp0073211(1,1,1)*Z(jj,3))
       temp833221(1,1,1)=IX*(aux833221(1,1,1)+8*temp0073322(1,1,1)*Z(jj,
     &  1)+4*(temp0073321(1,1,1)*Z(jj,2)+temp0073221(1,1,1)*Z(jj,3)))
       temp833221(2,1,1)=IX*(aux833221(2,1,1)+2*temp0073322(2,1,1)*Z(jj,
     &  1)+2*temp0073322(1,1,1)*Z(jj,2)+4*(temp0073322(1,2,1)*Z(jj,1)+te
     &  mp0073321(2,1,1)*Z(jj,2))+4*temp0073221(2,1,1)*Z(jj,3))
       temp833221(2,2,1)=IX*(aux833221(2,2,1)+2*(temp0073322(1,2,2)*Z(jj
     &  ,1)+temp0073322(2,2,1)*Z(jj,1))+4*temp0073321(2,2,1)*Z(jj,2)+4*t
     &  emp0073322(1,2,1)*Z(jj,2)+4*temp0073221(2,2,1)*Z(jj,3))
       temp833221(2,2,2)=IX*(aux833221(2,2,2)+2*temp0073322(2,2,2)*Z(jj,
     &  1)+6*temp0073322(1,2,2)*Z(jj,2)+4*(temp0073321(2,2,2)*Z(jj,2)+te
     &  mp0073221(2,2,2)*Z(jj,3)))
       temp833222(1,1,1)=IX*(aux833222(1,1,1)+6*(temp0073322(2,1,1)*Z(jj
     &  ,1)+temp0073322(1,1,1)*Z(jj,2))+4*temp0073222(1,1,1)*Z(jj,3))
       temp833222(2,1,1)=IX*(aux833222(2,1,1)+8*temp0073322(2,1,1)*Z(jj,
     &  2)+4*(temp0073322(2,2,1)*Z(jj,1)+temp0073222(2,1,1)*Z(jj,3)))
       temp833222(2,2,1)=IX*(aux833222(2,2,1)+2*temp0073322(2,2,2)*Z(jj,
     &  1)+10*temp0073322(2,2,1)*Z(jj,2)+4*temp0073222(2,2,1)*Z(jj,3))
       temp833222(2,2,2)=IX*(aux833222(2,2,2)+12*temp0073322(2,2,2)*Z(jj
     &  ,2)+4*temp0073222(2,2,2)*Z(jj,3))
       temp833311(1,1,1)=IX*(aux833311(1,1,1)+10*temp0073331(1,1,1)*Z(jj
     &  ,1)+6*temp0073311(1,1,1)*Z(jj,3))
       temp833321(1,1,1)=IX*(aux833321(1,1,1)+8*temp0073332(1,1,1)*Z(jj,
     &  1)+2*temp0073331(1,1,1)*Z(jj,2)+6*temp0073321(1,1,1)*Z(jj,3))
       temp833321(2,1,1)=IX*(aux833321(2,1,1)+4*temp0073332(1,2,1)*Z(jj,
     &  1)+2*temp0073332(1,1,1)*Z(jj,2)+2*(temp0073332(2,1,1)*Z(jj,1)+te
     &  mp0073331(2,1,1)*Z(jj,2))+6*temp0073321(2,1,1)*Z(jj,3))
       temp833321(2,2,1)=IX*(aux833321(2,2,1)+2*(temp0073332(1,2,2)*Z(jj
     &  ,1)+temp0073332(2,2,1)*Z(jj,1))+2*temp0073331(2,2,1)*Z(jj,2)+4*t
     &  emp0073332(1,2,1)*Z(jj,2)+6*temp0073321(2,2,1)*Z(jj,3))
       temp833321(2,2,2)=IX*(aux833321(2,2,2)+6*temp0073332(1,2,2)*Z(jj,
     &  2)+2*(temp0073332(2,2,2)*Z(jj,1)+temp0073331(2,2,2)*Z(jj,2))+6*t
     &  emp0073321(2,2,2)*Z(jj,3))
       temp833322(1,1,1)=IX*(aux833322(1,1,1)+4*temp0073332(1,1,1)*Z(jj,
     &  2)+6*(temp0073332(2,1,1)*Z(jj,1)+temp0073322(1,1,1)*Z(jj,3)))
       temp833322(2,1,1)=IX*(aux833322(2,1,1)+4*temp0073332(2,2,1)*Z(jj,
     &  1)+6*(temp0073332(2,1,1)*Z(jj,2)+temp0073322(2,1,1)*Z(jj,3)))
       temp833322(2,2,1)=IX*(aux833322(2,2,1)+2*temp0073332(2,2,2)*Z(jj,
     &  1)+8*temp0073332(2,2,1)*Z(jj,2)+6*temp0073322(2,2,1)*Z(jj,3))
       temp833322(2,2,2)=IX*(aux833322(2,2,2)+10*temp0073332(2,2,2)*Z(jj
     &  ,2)+6*temp0073322(2,2,2)*Z(jj,3))
       temp833331(1,1,1)=IX*(aux833331(1,1,1)+8*(temp0073333(1,1,1)*Z(jj
     &  ,1)+temp0073331(1,1,1)*Z(jj,3)))
       temp833331(2,1,1)=IX*(aux833331(2,1,1)+4*temp0073333(1,2,1)*Z(jj,
     &  1)+2*(temp0073333(2,1,1)*Z(jj,1)+temp0073333(1,1,1)*Z(jj,2))+8*t
     &  emp0073331(2,1,1)*Z(jj,3))
       temp833331(2,2,1)=IX*(aux833331(2,2,1)+2*(temp0073333(1,2,2)*Z(jj
     &  ,1)+temp0073333(2,2,1)*Z(jj,1))+4*temp0073333(1,2,1)*Z(jj,2)+8*t
     &  emp0073331(2,2,1)*Z(jj,3))
       temp833331(2,2,2)=IX*(aux833331(2,2,2)+2*temp0073333(2,2,2)*Z(jj,
     &  1)+6*temp0073333(1,2,2)*Z(jj,2)+8*temp0073331(2,2,2)*Z(jj,3))
       temp833331(3,1,1)=IX*(aux833331(3,1,1)+4*temp0073333(1,3,1)*Z(jj,
     &  1)+8*temp0073331(3,1,1)*Z(jj,3)+2*(temp0073333(3,1,1)*Z(jj,1)+te
     &  mp0073333(1,1,1)*Z(jj,3)))
       temp833331(3,2,1)=IX*(aux833331(3,2,1)+2*(temp0073333(1,3,2)*Z(jj
     &  ,1)+temp0073333(3,2,1)*Z(jj,1))+2*temp0073333(1,3,1)*Z(jj,2)+8*t
     &  emp0073331(3,2,1)*Z(jj,3)+2*temp0073333(1,2,1)*Z(jj,3))
       temp833331(3,2,2)=IX*(aux833331(3,2,2)+4*temp0073333(1,3,2)*Z(jj,
     &  2)+8*temp0073331(3,2,2)*Z(jj,3)+2*(temp0073333(3,2,2)*Z(jj,1)+te
     &  mp0073333(1,2,2)*Z(jj,3)))
       temp833331(3,3,1)=IX*(aux833331(3,3,1)+2*(temp0073333(1,3,3)*Z(jj
     &  ,1)+temp0073333(3,3,1)*Z(jj,1))+8*temp0073331(3,3,1)*Z(jj,3)+4*t
     &  emp0073333(1,3,1)*Z(jj,3))
       temp833331(3,3,2)=IX*(aux833331(3,3,2)+2*(temp0073333(3,3,2)*Z(jj
     &  ,1)+temp0073333(1,3,3)*Z(jj,2))+8*temp0073331(3,3,2)*Z(jj,3)+4*t
     &  emp0073333(1,3,2)*Z(jj,3))
       temp833331(3,3,3)=IX*(aux833331(3,3,3)+2*temp0073333(3,3,3)*Z(jj,
     &  1)+8*temp0073331(3,3,3)*Z(jj,3)+6*temp0073333(1,3,3)*Z(jj,3))
       temp833332(1,1,1)=IX*(aux833332(1,1,1)+6*temp0073333(2,1,1)*Z(jj,
     &  1)+2*temp0073333(1,1,1)*Z(jj,2)+8*temp0073332(1,1,1)*Z(jj,3))
       temp833332(2,1,1)=IX*(aux833332(2,1,1)+4*(temp0073333(2,2,1)*Z(jj
     &  ,1)+temp0073333(2,1,1)*Z(jj,2))+8*temp0073332(2,1,1)*Z(jj,3))
       temp833332(2,2,1)=IX*(aux833332(2,2,1)+2*temp0073333(2,2,2)*Z(jj,
     &  1)+6*temp0073333(2,2,1)*Z(jj,2)+8*temp0073332(2,2,1)*Z(jj,3))
       temp833332(2,2,2)=IX*(aux833332(2,2,2)+8*(temp0073333(2,2,2)*Z(jj
     &  ,2)+temp0073332(2,2,2)*Z(jj,3)))
       temp833332(3,1,1)=IX*(aux833332(3,1,1)+4*temp0073333(2,3,1)*Z(jj,
     &  1)+8*temp0073332(3,1,1)*Z(jj,3)+2*(temp0073333(3,1,1)*Z(jj,2)+te
     &  mp0073333(2,1,1)*Z(jj,3)))
       temp833332(3,2,1)=IX*(aux833332(3,2,1)+2*temp0073333(3,2,1)*Z(jj,
     &  2)+2*(temp0073333(2,3,2)*Z(jj,1)+temp0073333(2,3,1)*Z(jj,2))+8*t
     &  emp0073332(3,2,1)*Z(jj,3)+2*temp0073333(2,2,1)*Z(jj,3))
       temp833332(3,2,2)=IX*(aux833332(3,2,2)+4*temp0073333(2,3,2)*Z(jj,
     &  2)+8*temp0073332(3,2,2)*Z(jj,3)+2*(temp0073333(3,2,2)*Z(jj,2)+te
     &  mp0073333(2,2,2)*Z(jj,3)))
       temp833332(3,3,1)=IX*(aux833332(3,3,1)+2*(temp0073333(2,3,3)*Z(jj
     &  ,1)+temp0073333(3,3,1)*Z(jj,2))+8*temp0073332(3,3,1)*Z(jj,3)+4*t
     &  emp0073333(2,3,1)*Z(jj,3))
       temp833332(3,3,2)=IX*(aux833332(3,3,2)+2*(temp0073333(2,3,3)*Z(jj
     &  ,2)+temp0073333(3,3,2)*Z(jj,2))+8*temp0073332(3,3,2)*Z(jj,3)+4*t
     &  emp0073333(2,3,2)*Z(jj,3))
       temp833332(3,3,3)=IX*(aux833332(3,3,3)+2*temp0073333(3,3,3)*Z(jj,
     &  2)+8*temp0073332(3,3,3)*Z(jj,3)+6*temp0073333(2,3,3)*Z(jj,3))
       temp833333(1,1,1)=IX*(aux833333(1,1,1)+6*temp0073333(3,1,1)*Z(jj,
     &  1)+10*temp0073333(1,1,1)*Z(jj,3))
       temp833333(2,1,1)=IX*(aux833333(2,1,1)+4*temp0073333(3,2,1)*Z(jj,
     &  1)+2*temp0073333(3,1,1)*Z(jj,2)+10*temp0073333(2,1,1)*Z(jj,3))
       temp833333(2,2,1)=IX*(aux833333(2,2,1)+2*temp0073333(3,2,2)*Z(jj,
     &  1)+4*temp0073333(3,2,1)*Z(jj,2)+10*temp0073333(2,2,1)*Z(jj,3))
       temp833333(2,2,2)=IX*(aux833333(2,2,2)+6*temp0073333(3,2,2)*Z(jj,
     &  2)+10*temp0073333(2,2,2)*Z(jj,3))
       temp833333(3,1,1)=IX*(aux833333(3,1,1)+4*temp0073333(3,3,1)*Z(jj,
     &  1)+12*temp0073333(3,1,1)*Z(jj,3))
       temp833333(3,2,1)=IX*(aux833333(3,2,1)+2*(temp0073333(3,3,2)*Z(jj
     &  ,1)+temp0073333(3,3,1)*Z(jj,2))+12*temp0073333(3,2,1)*Z(jj,3))
       temp833333(3,2,2)=IX*(aux833333(3,2,2)+4*temp0073333(3,3,2)*Z(jj,
     &  2)+12*temp0073333(3,2,2)*Z(jj,3))
       temp833333(3,3,1)=IX*(aux833333(3,3,1)+2*temp0073333(3,3,3)*Z(jj,
     &  1)+14*temp0073333(3,3,1)*Z(jj,3))
       temp833333(3,3,2)=IX*(aux833333(3,3,2)+2*temp0073333(3,3,3)*Z(jj,
     &  2)+14*temp0073333(3,3,2)*Z(jj,3))
       temp833333(3,3,3)=IX*(aux833333(3,3,3)+16*temp0073333(3,3,3)*Z(jj
     &  ,3))
       temp811111(1,1,2)=temp821111(1,1,1)
       temp811111(1,1,3)=temp831111(1,1,1)
       temp811111(1,2,1)=temp821111(1,1,1)
       temp811111(1,2,2)=temp822111(1,1,1)
       temp811111(1,2,3)=temp832111(1,1,1)
       temp811111(1,3,1)=temp831111(1,1,1)
       temp811111(1,3,2)=temp832111(1,1,1)
       temp811111(1,3,3)=temp833111(1,1,1)
       temp821111(1,1,2)=temp822111(1,1,1)
       temp821111(1,1,3)=temp832111(1,1,1)
       temp821111(1,2,1)=temp822111(1,1,1)
       temp821111(1,2,2)=temp822211(1,1,1)
       temp821111(1,2,3)=temp832211(1,1,1)
       temp821111(1,3,1)=temp832111(1,1,1)
       temp821111(1,3,2)=temp832211(1,1,1)
       temp821111(1,3,3)=temp833211(1,1,1)
       temp822111(1,1,2)=temp822211(1,1,1)
       temp822111(1,1,3)=temp832211(1,1,1)
       temp822111(1,2,1)=temp822211(1,1,1)
       temp822111(1,2,2)=temp822221(1,1,1)
       temp822111(1,2,3)=temp832221(1,1,1)
       temp822111(1,3,1)=temp832211(1,1,1)
       temp822111(1,3,2)=temp832221(1,1,1)
       temp822111(1,3,3)=temp833221(1,1,1)
       temp822211(1,1,2)=temp822221(1,1,1)
       temp822211(1,1,3)=temp832221(1,1,1)
       temp822211(1,2,1)=temp822221(1,1,1)
       temp822211(1,2,2)=temp822222(1,1,1)
       temp822211(1,2,3)=temp832222(1,1,1)
       temp822211(1,3,1)=temp832221(1,1,1)
       temp822211(1,3,2)=temp832222(1,1,1)
       temp822211(1,3,3)=temp833222(1,1,1)
       temp822221(1,1,2)=temp822222(1,1,1)
       temp822221(1,1,3)=temp832222(1,1,1)
       temp822221(1,2,1)=temp822222(1,1,1)
       temp822221(1,2,2)=temp822222(2,1,1)
       temp822221(1,2,3)=temp832222(2,1,1)
       temp822221(1,3,1)=temp832222(1,1,1)
       temp822221(1,3,2)=temp832222(2,1,1)
       temp822221(1,3,3)=temp833222(2,1,1)
       temp822222(1,1,2)=temp822222(2,1,1)
       temp822222(1,1,3)=temp832222(2,1,1)
       temp822222(1,2,1)=temp822222(2,1,1)
       temp822222(1,2,2)=temp822222(2,2,1)
       temp822222(1,2,3)=temp832222(2,2,1)
       temp822222(1,3,1)=temp832222(2,1,1)
       temp822222(1,3,2)=temp832222(2,2,1)
       temp822222(1,3,3)=temp833222(2,2,1)
       temp822222(2,1,2)=temp822222(2,2,1)
       temp822222(2,1,3)=temp832222(2,2,1)
       temp822222(2,2,3)=temp832222(2,2,2)
       temp822222(2,3,1)=temp832222(2,2,1)
       temp822222(2,3,2)=temp832222(2,2,2)
       temp822222(2,3,3)=temp833222(2,2,2)
       temp831111(1,1,2)=temp832111(1,1,1)
       temp831111(1,1,3)=temp833111(1,1,1)
       temp831111(1,2,1)=temp832111(1,1,1)
       temp831111(1,2,2)=temp832211(1,1,1)
       temp831111(1,2,3)=temp833211(1,1,1)
       temp831111(1,3,1)=temp833111(1,1,1)
       temp831111(1,3,2)=temp833211(1,1,1)
       temp831111(1,3,3)=temp833311(1,1,1)
       temp832111(1,1,2)=temp832211(1,1,1)
       temp832111(1,1,3)=temp833211(1,1,1)
       temp832111(1,2,1)=temp832211(1,1,1)
       temp832111(1,2,2)=temp832221(1,1,1)
       temp832111(1,2,3)=temp833221(1,1,1)
       temp832111(1,3,1)=temp833211(1,1,1)
       temp832111(1,3,2)=temp833221(1,1,1)
       temp832111(1,3,3)=temp833321(1,1,1)
       temp832211(1,1,2)=temp832221(1,1,1)
       temp832211(1,1,3)=temp833221(1,1,1)
       temp832211(1,2,1)=temp832221(1,1,1)
       temp832211(1,2,2)=temp832222(1,1,1)
       temp832211(1,2,3)=temp833222(1,1,1)
       temp832211(1,3,1)=temp833221(1,1,1)
       temp832211(1,3,2)=temp833222(1,1,1)
       temp832211(1,3,3)=temp833322(1,1,1)
       temp832221(1,1,2)=temp832222(1,1,1)
       temp832221(1,1,3)=temp833222(1,1,1)
       temp832221(1,2,1)=temp832222(1,1,1)
       temp832221(1,2,2)=temp832222(2,1,1)
       temp832221(1,2,3)=temp833222(2,1,1)
       temp832221(1,3,1)=temp833222(1,1,1)
       temp832221(1,3,2)=temp833222(2,1,1)
       temp832221(1,3,3)=temp833322(2,1,1)
       temp832222(1,1,2)=temp832222(2,1,1)
       temp832222(1,1,3)=temp833222(2,1,1)
       temp832222(1,2,1)=temp832222(2,1,1)
       temp832222(1,2,2)=temp832222(2,2,1)
       temp832222(1,2,3)=temp833222(2,2,1)
       temp832222(1,3,1)=temp833222(2,1,1)
       temp832222(1,3,2)=temp833222(2,2,1)
       temp832222(1,3,3)=temp833322(2,2,1)
       temp832222(2,1,2)=temp832222(2,2,1)
       temp832222(2,1,3)=temp833222(2,2,1)
       temp832222(2,2,3)=temp833222(2,2,2)
       temp832222(2,3,1)=temp833222(2,2,1)
       temp832222(2,3,2)=temp833222(2,2,2)
       temp832222(2,3,3)=temp833322(2,2,2)
       temp833111(1,1,2)=temp833211(1,1,1)
       temp833111(1,1,3)=temp833311(1,1,1)
       temp833111(1,2,1)=temp833211(1,1,1)
       temp833111(1,2,2)=temp833221(1,1,1)
       temp833111(1,2,3)=temp833321(1,1,1)
       temp833111(1,3,1)=temp833311(1,1,1)
       temp833111(1,3,2)=temp833321(1,1,1)
       temp833111(1,3,3)=temp833331(1,1,1)
       temp833211(1,1,2)=temp833221(1,1,1)
       temp833211(1,1,3)=temp833321(1,1,1)
       temp833211(1,2,1)=temp833221(1,1,1)
       temp833211(1,2,2)=temp833222(1,1,1)
       temp833211(1,2,3)=temp833322(1,1,1)
       temp833211(1,3,1)=temp833321(1,1,1)
       temp833211(1,3,2)=temp833322(1,1,1)
       temp833211(1,3,3)=temp833332(1,1,1)
       temp833221(1,1,2)=temp833222(1,1,1)
       temp833221(1,1,3)=temp833322(1,1,1)
       temp833221(1,2,1)=temp833222(1,1,1)
       temp833221(1,2,2)=temp833222(2,1,1)
       temp833221(1,2,3)=temp833322(2,1,1)
       temp833221(1,3,1)=temp833322(1,1,1)
       temp833221(1,3,2)=temp833322(2,1,1)
       temp833221(1,3,3)=temp833332(2,1,1)
       temp833222(1,1,2)=temp833222(2,1,1)
       temp833222(1,1,3)=temp833322(2,1,1)
       temp833222(1,2,1)=temp833222(2,1,1)
       temp833222(1,2,2)=temp833222(2,2,1)
       temp833222(1,2,3)=temp833322(2,2,1)
       temp833222(1,3,1)=temp833322(2,1,1)
       temp833222(1,3,2)=temp833322(2,2,1)
       temp833222(1,3,3)=temp833332(2,2,1)
       temp833222(2,1,2)=temp833222(2,2,1)
       temp833222(2,1,3)=temp833322(2,2,1)
       temp833222(2,2,3)=temp833322(2,2,2)
       temp833222(2,3,1)=temp833322(2,2,1)
       temp833222(2,3,2)=temp833322(2,2,2)
       temp833222(2,3,3)=temp833332(2,2,2)
       temp833311(1,1,2)=temp833321(1,1,1)
       temp833311(1,1,3)=temp833331(1,1,1)
       temp833311(1,2,1)=temp833321(1,1,1)
       temp833311(1,2,2)=temp833322(1,1,1)
       temp833311(1,2,3)=temp833332(1,1,1)
       temp833311(1,3,1)=temp833331(1,1,1)
       temp833311(1,3,2)=temp833332(1,1,1)
       temp833311(1,3,3)=temp833333(1,1,1)
       temp833321(1,1,2)=temp833322(1,1,1)
       temp833321(1,1,3)=temp833332(1,1,1)
       temp833321(1,2,1)=temp833322(1,1,1)
       temp833321(1,2,2)=temp833322(2,1,1)
       temp833321(1,2,3)=temp833332(2,1,1)
       temp833321(1,3,1)=temp833332(1,1,1)
       temp833321(1,3,2)=temp833332(2,1,1)
       temp833321(1,3,3)=temp833333(2,1,1)
       temp833322(1,1,2)=temp833322(2,1,1)
       temp833322(1,1,3)=temp833332(2,1,1)
       temp833322(1,2,1)=temp833322(2,1,1)
       temp833322(1,2,2)=temp833322(2,2,1)
       temp833322(1,2,3)=temp833332(2,2,1)
       temp833322(1,3,1)=temp833332(2,1,1)
       temp833322(1,3,2)=temp833332(2,2,1)
       temp833322(1,3,3)=temp833333(2,2,1)
       temp833322(2,1,2)=temp833322(2,2,1)
       temp833322(2,1,3)=temp833332(2,2,1)
       temp833322(2,2,3)=temp833332(2,2,2)
       temp833322(2,3,1)=temp833332(2,2,1)
       temp833322(2,3,2)=temp833332(2,2,2)
       temp833322(2,3,3)=temp833333(2,2,2)
       temp833331(1,1,2)=temp833332(1,1,1)
       temp833331(1,1,3)=temp833333(1,1,1)
       temp833331(1,2,1)=temp833332(1,1,1)
       temp833331(1,2,2)=temp833332(2,1,1)
       temp833331(1,2,3)=temp833333(2,1,1)
       temp833331(1,3,1)=temp833333(1,1,1)
       temp833331(1,3,2)=temp833333(2,1,1)
       temp833331(1,3,3)=temp833333(3,1,1)
       temp833332(1,1,2)=temp833332(2,1,1)
       temp833332(1,1,3)=temp833333(2,1,1)
       temp833332(1,2,1)=temp833332(2,1,1)
       temp833332(1,2,2)=temp833332(2,2,1)
       temp833332(1,2,3)=temp833333(2,2,1)
       temp833332(1,3,1)=temp833333(2,1,1)
       temp833332(1,3,2)=temp833333(2,2,1)
       temp833332(1,3,3)=temp833333(3,2,1)
       temp833332(2,1,2)=temp833332(2,2,1)
       temp833332(2,1,3)=temp833333(2,2,1)
       temp833332(2,2,3)=temp833333(2,2,2)
       temp833332(2,3,1)=temp833333(2,2,1)
       temp833332(2,3,2)=temp833333(2,2,2)
       temp833332(2,3,3)=temp833333(3,2,2)
       temp833333(1,1,2)=temp833333(2,1,1)
       temp833333(1,1,3)=temp833333(3,1,1)
       temp833333(1,2,1)=temp833333(2,1,1)
       temp833333(1,2,2)=temp833333(2,2,1)
       temp833333(1,2,3)=temp833333(3,2,1)
       temp833333(1,3,1)=temp833333(3,1,1)
       temp833333(1,3,2)=temp833333(3,2,1)
       temp833333(1,3,3)=temp833333(3,3,1)
       temp833333(2,1,2)=temp833333(2,2,1)
       temp833333(2,1,3)=temp833333(3,2,1)
       temp833333(2,2,3)=temp833333(3,2,2)
       temp833333(2,3,1)=temp833333(3,2,1)
       temp833333(2,3,2)=temp833333(3,2,2)
       temp833333(2,3,3)=temp833333(3,3,2)
       temp833333(3,1,2)=temp833333(3,2,1)
       temp833333(3,1,3)=temp833333(3,3,1)
       temp833333(3,2,3)=temp833333(3,3,2)
c                Step2
       tempD400000000=I16Z*(auxD400000000+tempD4000000*F(5)-det4*temp000
     &  0002(k,l))
       temp0000002(1,1)=I20Z*(aux0000002(1,1)+4*F(4)*temp0000001(1)+F(5)
     &  *temp00002(1,1)-det4*temp000041(1,k,l)+8*tempD400000000*ZZ(k,1,l
     &  ,1))
       temp0000002(2,1)=I20Z*(aux0000002(2,1)+2*(F(6)*temp0000001(1)+F(4
     &  )*temp0000001(2))+F(5)*temp00002(2,1)-det4*temp000042(1,k,l)+4*t
     &  empD400000000*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp0000002(2,2)=I20Z*(aux0000002(2,2)+4*F(6)*temp0000001(2)+F(5)
     &  *temp00002(2,2)-det4*temp000042(2,k,l)+8*tempD400000000*ZZ(k,2,l
     &  ,2))
       temp0000002(3,1)=I20Z*(aux0000002(3,1)+2*(F(7)*temp0000001(1)+F(4
     &  )*temp0000001(3))+F(5)*temp00002(3,1)-det4*temp000043(1,k,l)+4*t
     &  empD400000000*(ZZ(k,1,l,3)+ZZ(k,3,l,1)))
       temp0000002(3,2)=I20Z*(aux0000002(3,2)+2*(F(7)*temp0000001(2)+F(6
     &  )*temp0000001(3))+F(5)*temp00002(3,2)-det4*temp000043(2,k,l)+4*t
     &  empD400000000*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp0000002(3,3)=I20Z*(aux0000002(3,3)+4*F(7)*temp0000001(3)+F(5)
     &  *temp00002(3,3)-det4*temp000043(3,k,l)+8*tempD400000000*ZZ(k,3,l
     &  ,3))
       temp0000002(1,2)=temp0000002(2,1)
       temp0000002(1,3)=temp0000002(3,1)
       temp0000002(2,3)=temp0000002(3,2)
       temp000041(1,1,1)=I24Z*(aux000041(1,1,1)+8*F(4)*temp00003(1,1,1)+
     &  F(5)*temp0041(1,1,1)-det4*temp006111(1,k,l)+48*temp0000002(1,1)*
     &  ZZ(k,1,l,1))
       temp000042(1,1,1)=I24Z*(aux000042(1,1,1)+2*F(6)*temp00003(1,1,1)+
     &  6*F(4)*temp00003(2,1,1)+F(5)*temp0042(1,1,1)-det4*temp006211(1,k
     &  ,l)+24*temp0000002(2,1)*ZZ(k,1,l,1)+12*temp0000002(1,1)*(ZZ(k,1,
     &  l,2)+ZZ(k,2,l,1)))
       temp000042(2,1,1)=I24Z*(aux000042(2,1,1)+4*F(6)*temp00003(2,1,1)+
     &  4*F(4)*temp00003(2,2,1)+F(5)*temp0042(2,1,1)-det4*temp006221(1,k
     &  ,l)+16*temp0000002(2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp0000002
     &  (2,2)*ZZ(k,1,l,1)+temp0000002(1,1)*ZZ(k,2,l,2)))
       temp000042(2,2,1)=I24Z*(aux000042(2,2,1)+6*F(6)*temp00003(2,2,1)+
     &  2*F(4)*temp00003(2,2,2)+F(5)*temp0042(2,2,1)-det4*temp006222(1,k
     &  ,l)+12*temp0000002(2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*temp0000002
     &  (2,1)*ZZ(k,2,l,2))
       temp000042(2,2,2)=I24Z*(aux000042(2,2,2)+8*F(6)*temp00003(2,2,2)+
     &  F(5)*temp0042(2,2,2)-det4*temp006222(2,k,l)+48*temp0000002(2,2)*
     &  ZZ(k,2,l,2))
       temp000043(1,1,1)=I24Z*(aux000043(1,1,1)+2*F(7)*temp00003(1,1,1)+
     &  6*F(4)*temp00003(3,1,1)+F(5)*temp0043(1,1,1)-det4*temp006311(1,k
     &  ,l)+24*temp0000002(3,1)*ZZ(k,1,l,1)+12*temp0000002(1,1)*(ZZ(k,1,
     &  l,3)+ZZ(k,3,l,1)))
       temp000043(2,1,1)=I24Z*(aux000043(2,1,1)+2*F(7)*temp00003(2,1,1)+
     &  2*F(6)*temp00003(3,1,1)+4*F(4)*temp00003(3,2,1)+F(5)*temp0043(2,
     &  1,1)-det4*temp006321(1,k,l)+8*(temp0000002(3,2)*ZZ(k,1,l,1)+temp
     &  0000002(3,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))+8*temp0000002(2,1)*(ZZ(k
     &  ,1,l,3)+ZZ(k,3,l,1))+4*temp0000002(1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)
     &  ))
       temp000043(2,2,1)=I24Z*(aux000043(2,2,1)+2*F(7)*temp00003(2,2,1)+
     &  4*F(6)*temp00003(3,2,1)+2*F(4)*temp00003(3,2,2)+F(5)*temp0043(2,
     &  2,1)-det4*temp006322(1,k,l)+8*(temp0000002(3,2)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1))+temp0000002(3,1)*ZZ(k,2,l,2))+4*temp0000002(2,2)*(ZZ(k
     &  ,1,l,3)+ZZ(k,3,l,1))+8*temp0000002(2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)
     &  ))
       temp000043(2,2,2)=I24Z*(aux000043(2,2,2)+2*F(7)*temp00003(2,2,2)+
     &  6*F(6)*temp00003(3,2,2)+F(5)*temp0043(2,2,2)-det4*temp006322(2,k
     &  ,l)+24*temp0000002(3,2)*ZZ(k,2,l,2)+12*temp0000002(2,2)*(ZZ(k,2,
     &  l,3)+ZZ(k,3,l,2)))
       temp000043(3,1,1)=I24Z*(aux000043(3,1,1)+4*F(7)*temp00003(3,1,1)+
     &  4*F(4)*temp00003(3,3,1)+F(5)*temp0043(3,1,1)-det4*temp006331(1,k
     &  ,l)+16*temp0000002(3,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*(temp0000002
     &  (3,3)*ZZ(k,1,l,1)+temp0000002(1,1)*ZZ(k,3,l,3)))
       temp000043(3,2,1)=I24Z*(aux000043(3,2,1)+2*F(6)*temp00003(3,3,1)+
     &  2*F(4)*temp00003(3,3,2)+F(5)*temp0043(3,2,1)-det4*temp006332(1,k
     &  ,l)+4*(F(7)*temp00003(3,2,1)+temp0000002(3,3)*(ZZ(k,1,l,2)+ZZ(k,
     &  2,l,1)))+8*temp0000002(3,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*temp0000
     &  002(3,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*temp0000002(2,1)*ZZ(k,3,l,3
     &  ))
       temp000043(3,2,2)=I24Z*(aux000043(3,2,2)+4*F(7)*temp00003(3,2,2)+
     &  4*F(6)*temp00003(3,3,2)+F(5)*temp0043(3,2,2)-det4*temp006332(2,k
     &  ,l)+16*temp0000002(3,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*(temp0000002
     &  (3,3)*ZZ(k,2,l,2)+temp0000002(2,2)*ZZ(k,3,l,3)))
       temp000043(3,3,1)=I24Z*(aux000043(3,3,1)+6*F(7)*temp00003(3,3,1)+
     &  2*F(4)*temp00003(3,3,3)+F(5)*temp0043(3,3,1)-det4*temp006333(1,k
     &  ,l)+12*temp0000002(3,3)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+24*temp0000002
     &  (3,1)*ZZ(k,3,l,3))
       temp000043(3,3,2)=I24Z*(aux000043(3,3,2)+6*F(7)*temp00003(3,3,2)+
     &  2*F(6)*temp00003(3,3,3)+F(5)*temp0043(3,3,2)-det4*temp006333(2,k
     &  ,l)+12*temp0000002(3,3)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+24*temp0000002
     &  (3,2)*ZZ(k,3,l,3))
       temp000043(3,3,3)=I24Z*(aux000043(3,3,3)+8*F(7)*temp00003(3,3,3)+
     &  F(5)*temp0043(3,3,3)-det4*temp006333(3,k,l)+48*temp0000002(3,3)*
     &  ZZ(k,3,l,3))
       temp000041(1,1,2)=temp000042(1,1,1)
       temp000041(1,1,3)=temp000043(1,1,1)
       temp000041(1,2,1)=temp000042(1,1,1)
       temp000041(1,2,2)=temp000042(2,1,1)
       temp000041(1,2,3)=temp000043(2,1,1)
       temp000041(1,3,1)=temp000043(1,1,1)
       temp000041(1,3,2)=temp000043(2,1,1)
       temp000041(1,3,3)=temp000043(3,1,1)
       temp000042(1,1,2)=temp000042(2,1,1)
       temp000042(1,1,3)=temp000043(2,1,1)
       temp000042(1,2,1)=temp000042(2,1,1)
       temp000042(1,2,2)=temp000042(2,2,1)
       temp000042(1,2,3)=temp000043(2,2,1)
       temp000042(1,3,1)=temp000043(2,1,1)
       temp000042(1,3,2)=temp000043(2,2,1)
       temp000042(1,3,3)=temp000043(3,2,1)
       temp000042(2,1,2)=temp000042(2,2,1)
       temp000042(2,1,3)=temp000043(2,2,1)
       temp000042(2,2,3)=temp000043(2,2,2)
       temp000042(2,3,1)=temp000043(2,2,1)
       temp000042(2,3,2)=temp000043(2,2,2)
       temp000042(2,3,3)=temp000043(3,2,2)
       temp000043(1,1,2)=temp000043(2,1,1)
       temp000043(1,1,3)=temp000043(3,1,1)
       temp000043(1,2,1)=temp000043(2,1,1)
       temp000043(1,2,2)=temp000043(2,2,1)
       temp000043(1,2,3)=temp000043(3,2,1)
       temp000043(1,3,1)=temp000043(3,1,1)
       temp000043(1,3,2)=temp000043(3,2,1)
       temp000043(1,3,3)=temp000043(3,3,1)
       temp000043(2,1,2)=temp000043(2,2,1)
       temp000043(2,1,3)=temp000043(3,2,1)
       temp000043(2,2,3)=temp000043(3,2,2)
       temp000043(2,3,1)=temp000043(3,2,1)
       temp000043(2,3,2)=temp000043(3,2,2)
       temp000043(2,3,3)=temp000043(3,3,2)
       temp000043(3,1,2)=temp000043(3,2,1)
       temp000043(3,1,3)=temp000043(3,3,1)
       temp000043(3,2,3)=temp000043(3,3,2)
       temp006111(1,1,1)=I28Z*(aux006111(1,1,1)+12*F(4)*temp00511(1,1,1)
     &  +F(5)*temp6111(1,1,1)-det4*temp811111(1,k,l)+120*temp000041(1,1,
     &  1)*ZZ(k,1,l,1))
       temp006211(1,1,1)=I28Z*(aux006211(1,1,1)+2*F(6)*temp00511(1,1,1)+
     &  10*F(4)*temp00521(1,1,1)+F(5)*temp6211(1,1,1)-det4*temp821111(1,
     &  k,l)+80*temp000042(1,1,1)*ZZ(k,1,l,1)+20*temp000041(1,1,1)*(ZZ(k
     &  ,1,l,2)+ZZ(k,2,l,1)))
       temp006221(1,1,1)=I28Z*(aux006221(1,1,1)+4*F(6)*temp00521(1,1,1)+
     &  F(5)*temp6221(1,1,1)-det4*temp822111(1,k,l)+48*temp000042(2,1,1)
     &  *ZZ(k,1,l,1)+32*temp000042(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(F
     &  (4)*temp00522(1,1,1)+temp000041(1,1,1)*ZZ(k,2,l,2)))
       temp006222(1,1,1)=I28Z*(aux006222(1,1,1)+6*F(6)*temp00522(1,1,1)+
     &  6*F(4)*temp00522(2,1,1)+F(5)*temp6222(1,1,1)-det4*temp822211(1,k
     &  ,l)+36*temp000042(2,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*(temp00004
     &  2(2,2,1)*ZZ(k,1,l,1)+temp000042(1,1,1)*ZZ(k,2,l,2)))
       temp006222(2,1,1)=I28Z*(aux006222(2,1,1)+4*F(4)*temp00522(2,2,1)+
     &  F(5)*temp6222(2,1,1)-det4*temp822221(1,k,l)+8*(F(6)*temp00522(2,
     &  1,1)+temp000042(2,2,2)*ZZ(k,1,l,1))+32*temp000042(2,2,1)*(ZZ(k,1
     &  ,l,2)+ZZ(k,2,l,1))+48*temp000042(2,1,1)*ZZ(k,2,l,2))
       temp006222(2,2,1)=I28Z*(aux006222(2,2,1)+10*F(6)*temp00522(2,2,1)
     &  +2*F(4)*temp00522(2,2,2)+F(5)*temp6222(2,2,1)-det4*temp822222(1,
     &  k,l)+20*temp000042(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+80*temp00004
     &  2(2,2,1)*ZZ(k,2,l,2))
       temp006222(2,2,2)=I28Z*(aux006222(2,2,2)+12*F(6)*temp00522(2,2,2)
     &  +F(5)*temp6222(2,2,2)-det4*temp822222(2,k,l)+120*temp000042(2,2,
     &  2)*ZZ(k,2,l,2))
       temp006311(1,1,1)=I28Z*(aux006311(1,1,1)+2*F(7)*temp00511(1,1,1)+
     &  10*F(4)*temp00531(1,1,1)+F(5)*temp6311(1,1,1)-det4*temp831111(1,
     &  k,l)+80*temp000043(1,1,1)*ZZ(k,1,l,1)+20*temp000041(1,1,1)*(ZZ(k
     &  ,1,l,3)+ZZ(k,3,l,1)))
       temp006321(1,1,1)=I28Z*(aux006321(1,1,1)+2*F(7)*temp00521(1,1,1)+
     &  2*F(6)*temp00531(1,1,1)+8*F(4)*temp00532(1,1,1)+F(5)*temp6321(1,
     &  1,1)-det4*temp832111(1,k,l)+48*temp000043(2,1,1)*ZZ(k,1,l,1)+16*
     &  (temp000043(1,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp000042(1,1,1)*(
     &  ZZ(k,1,l,3)+ZZ(k,3,l,1)))+4*temp000041(1,1,1)*(ZZ(k,2,l,3)+ZZ(k,
     &  3,l,2)))
       temp006322(1,1,1)=I28Z*(aux006322(1,1,1)+2*F(7)*temp00522(1,1,1)+
     &  4*F(6)*temp00532(1,1,1)+6*F(4)*temp00532(2,1,1)+F(5)*temp6322(1,
     &  1,1)-det4*temp832211(1,k,l)+24*(temp000043(2,2,1)*ZZ(k,1,l,1)+te
     &  mp000043(2,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))+8*temp000043(1,1,1)*Z
     &  Z(k,2,l,2)+12*temp000042(2,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*temp
     &  000042(1,1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp006322(2,1,1)=I28Z*(aux006322(2,1,1)+2*F(7)*temp00522(2,1,1)+
     &  6*F(6)*temp00532(2,1,1)+4*F(4)*temp00532(2,2,1)+F(5)*temp6322(2,
     &  1,1)-det4*temp832221(1,k,l)+24*temp000043(2,2,1)*(ZZ(k,1,l,2)+ZZ
     &  (k,2,l,1))+24*temp000043(2,1,1)*ZZ(k,2,l,2)+8*(temp000043(2,2,2)
     &  *ZZ(k,1,l,1)+temp000042(2,2,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)))+12*tem
     &  p000042(2,1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp006322(2,2,1)=I28Z*(aux006322(2,2,1)+2*F(7)*temp00522(2,2,1)+
     &  8*F(6)*temp00532(2,2,1)+2*F(4)*temp00532(2,2,2)+F(5)*temp6322(2,
     &  2,1)-det4*temp832222(1,k,l)+48*temp000043(2,2,1)*ZZ(k,2,l,2)+4*t
     &  emp000042(2,2,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+16*(temp000043(2,2,2)
     &  *(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp000042(2,2,1)*(ZZ(k,2,l,3)+ZZ(k,3
     &  ,l,2))))
       temp006322(2,2,2)=I28Z*(aux006322(2,2,2)+2*F(7)*temp00522(2,2,2)+
     &  10*F(6)*temp00532(2,2,2)+F(5)*temp6322(2,2,2)-det4*temp832222(2,
     &  k,l)+80*temp000043(2,2,2)*ZZ(k,2,l,2)+20*temp000042(2,2,2)*(ZZ(k
     &  ,2,l,3)+ZZ(k,3,l,2)))
       temp006331(1,1,1)=I28Z*(aux006331(1,1,1)+4*F(7)*temp00531(1,1,1)+
     &  F(5)*temp6331(1,1,1)-det4*temp833111(1,k,l)+48*temp000043(3,1,1)
     &  *ZZ(k,1,l,1)+32*temp000043(1,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*(F
     &  (4)*temp00533(1,1,1)+temp000041(1,1,1)*ZZ(k,3,l,3)))
       temp006332(1,1,1)=I28Z*(aux006332(1,1,1)+4*F(7)*temp00532(1,1,1)+
     &  2*F(6)*temp00533(1,1,1)+6*F(4)*temp00533(2,1,1)+F(5)*temp6332(1,
     &  1,1)-det4*temp833211(1,k,l)+12*temp000043(3,1,1)*(ZZ(k,1,l,2)+ZZ
     &  (k,2,l,1))+24*(temp000043(3,2,1)*ZZ(k,1,l,1)+temp000043(2,1,1)*(
     &  ZZ(k,1,l,3)+ZZ(k,3,l,1)))+8*temp000043(1,1,1)*(ZZ(k,2,l,3)+ZZ(k,
     &  3,l,2))+8*temp000042(1,1,1)*ZZ(k,3,l,3))
       temp006332(2,1,1)=I28Z*(aux006332(2,1,1)+4*F(7)*temp00532(2,1,1)+
     &  4*F(6)*temp00533(2,1,1)+4*F(4)*temp00533(2,2,1)+F(5)*temp6332(2,
     &  1,1)-det4*temp833221(1,k,l)+16*temp000043(3,2,1)*(ZZ(k,1,l,2)+ZZ
     &  (k,2,l,1))+8*(temp000043(3,2,2)*ZZ(k,1,l,1)+temp000043(3,1,1)*ZZ
     &  (k,2,l,2))+16*temp000043(2,2,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+16*tem
     &  p000043(2,1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*temp000042(2,1,1)*ZZ(
     &  k,3,l,3))
       temp006332(2,2,1)=I28Z*(aux006332(2,2,1)+4*F(7)*temp00532(2,2,1)+
     &  6*F(6)*temp00533(2,2,1)+2*F(4)*temp00533(2,2,2)+F(5)*temp6332(2,
     &  2,1)-det4*temp833222(1,k,l)+12*temp000043(3,2,2)*(ZZ(k,1,l,2)+ZZ
     &  (k,2,l,1))+8*temp000043(2,2,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+24*(tem
     &  p000043(3,2,1)*ZZ(k,2,l,2)+temp000043(2,2,1)*(ZZ(k,2,l,3)+ZZ(k,3
     &  ,l,2)))+8*temp000042(2,2,1)*ZZ(k,3,l,3))
       temp006332(2,2,2)=I28Z*(aux006332(2,2,2)+4*F(7)*temp00532(2,2,2)+
     &  F(5)*temp6332(2,2,2)-det4*temp833222(2,k,l)+48*temp000043(3,2,2)
     &  *ZZ(k,2,l,2)+32*temp000043(2,2,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*(F
     &  (6)*temp00533(2,2,2)+temp000042(2,2,2)*ZZ(k,3,l,3)))
       temp006333(1,1,1)=I28Z*(aux006333(1,1,1)+6*F(7)*temp00533(1,1,1)+
     &  6*F(4)*temp00533(3,1,1)+F(5)*temp6333(1,1,1)-det4*temp833311(1,k
     &  ,l)+36*temp000043(3,1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+24*(temp00004
     &  3(3,3,1)*ZZ(k,1,l,1)+temp000043(1,1,1)*ZZ(k,3,l,3)))
       temp006333(2,1,1)=I28Z*(aux006333(2,1,1)+6*F(7)*temp00533(2,1,1)+
     &  2*F(6)*temp00533(3,1,1)+4*F(4)*temp00533(3,2,1)+F(5)*temp6333(2,
     &  1,1)-det4*temp833321(1,k,l)+8*(temp000043(3,3,2)*ZZ(k,1,l,1)+tem
     &  p000043(3,3,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))+24*temp000043(3,2,1)*(
     &  ZZ(k,1,l,3)+ZZ(k,3,l,1))+12*temp000043(3,1,1)*(ZZ(k,2,l,3)+ZZ(k,
     &  3,l,2))+24*temp000043(2,1,1)*ZZ(k,3,l,3))
       temp006333(2,2,1)=I28Z*(aux006333(2,2,1)+6*F(7)*temp00533(2,2,1)+
     &  4*F(6)*temp00533(3,2,1)+2*F(4)*temp00533(3,2,2)+F(5)*temp6333(2,
     &  2,1)-det4*temp833322(1,k,l)+8*(temp000043(3,3,2)*(ZZ(k,1,l,2)+ZZ
     &  (k,2,l,1))+temp000043(3,3,1)*ZZ(k,2,l,2))+12*temp000043(3,2,2)*(
     &  ZZ(k,1,l,3)+ZZ(k,3,l,1))+24*temp000043(3,2,1)*(ZZ(k,2,l,3)+ZZ(k,
     &  3,l,2))+24*temp000043(2,2,1)*ZZ(k,3,l,3))
       temp006333(2,2,2)=I28Z*(aux006333(2,2,2)+6*F(7)*temp00533(2,2,2)+
     &  6*F(6)*temp00533(3,2,2)+F(5)*temp6333(2,2,2)-det4*temp833322(2,k
     &  ,l)+36*temp000043(3,2,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+24*(temp00004
     &  3(3,3,2)*ZZ(k,2,l,2)+temp000043(2,2,2)*ZZ(k,3,l,3)))
       temp006333(3,1,1)=I28Z*(aux006333(3,1,1)+4*F(4)*temp00533(3,3,1)+
     &  F(5)*temp6333(3,1,1)-det4*temp833331(1,k,l)+8*(F(7)*temp00533(3,
     &  1,1)+temp000043(3,3,3)*ZZ(k,1,l,1))+32*temp000043(3,3,1)*(ZZ(k,1
     &  ,l,3)+ZZ(k,3,l,1))+48*temp000043(3,1,1)*ZZ(k,3,l,3))
       temp006333(3,2,1)=I28Z*(aux006333(3,2,1)+8*F(7)*temp00533(3,2,1)+
     &  2*F(6)*temp00533(3,3,1)+2*F(4)*temp00533(3,3,2)+F(5)*temp6333(3,
     &  2,1)-det4*temp833332(1,k,l)+4*temp000043(3,3,3)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1))+16*(temp000043(3,3,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+temp00
     &  0043(3,3,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))+48*temp000043(3,2,1)*ZZ(k
     &  ,3,l,3))
       temp006333(3,2,2)=I28Z*(aux006333(3,2,2)+4*F(6)*temp00533(3,3,2)+
     &  F(5)*temp6333(3,2,2)-det4*temp833332(2,k,l)+8*(F(7)*temp00533(3,
     &  2,2)+temp000043(3,3,3)*ZZ(k,2,l,2))+32*temp000043(3,3,2)*(ZZ(k,2
     &  ,l,3)+ZZ(k,3,l,2))+48*temp000043(3,2,2)*ZZ(k,3,l,3))
       temp006333(3,3,1)=I28Z*(aux006333(3,3,1)+10*F(7)*temp00533(3,3,1)
     &  +2*F(4)*temp00533(3,3,3)+F(5)*temp6333(3,3,1)-det4*temp833333(1,
     &  k,l)+20*temp000043(3,3,3)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+80*temp00004
     &  3(3,3,1)*ZZ(k,3,l,3))
       temp006333(3,3,2)=I28Z*(aux006333(3,3,2)+10*F(7)*temp00533(3,3,2)
     &  +2*F(6)*temp00533(3,3,3)+F(5)*temp6333(3,3,2)-det4*temp833333(2,
     &  k,l)+20*temp000043(3,3,3)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+80*temp00004
     &  3(3,3,2)*ZZ(k,3,l,3))
       temp006333(3,3,3)=I28Z*(aux006333(3,3,3)+12*F(7)*temp00533(3,3,3)
     &  +F(5)*temp6333(3,3,3)-det4*temp833333(3,k,l)+120*temp000043(3,3,
     &  3)*ZZ(k,3,l,3))
       temp006111(1,1,2)=temp006211(1,1,1)
       temp006111(1,1,3)=temp006311(1,1,1)
       temp006111(1,2,1)=temp006211(1,1,1)
       temp006111(1,2,2)=temp006221(1,1,1)
       temp006111(1,2,3)=temp006321(1,1,1)
       temp006111(1,3,1)=temp006311(1,1,1)
       temp006111(1,3,2)=temp006321(1,1,1)
       temp006111(1,3,3)=temp006331(1,1,1)
       temp006211(1,1,2)=temp006221(1,1,1)
       temp006211(1,1,3)=temp006321(1,1,1)
       temp006211(1,2,1)=temp006221(1,1,1)
       temp006211(1,2,2)=temp006222(1,1,1)
       temp006211(1,2,3)=temp006322(1,1,1)
       temp006211(1,3,1)=temp006321(1,1,1)
       temp006211(1,3,2)=temp006322(1,1,1)
       temp006211(1,3,3)=temp006332(1,1,1)
       temp006221(1,1,2)=temp006222(1,1,1)
       temp006221(1,1,3)=temp006322(1,1,1)
       temp006221(1,2,1)=temp006222(1,1,1)
       temp006221(1,2,2)=temp006222(2,1,1)
       temp006221(1,2,3)=temp006322(2,1,1)
       temp006221(1,3,1)=temp006322(1,1,1)
       temp006221(1,3,2)=temp006322(2,1,1)
       temp006221(1,3,3)=temp006332(2,1,1)
       temp006222(1,1,2)=temp006222(2,1,1)
       temp006222(1,1,3)=temp006322(2,1,1)
       temp006222(1,2,1)=temp006222(2,1,1)
       temp006222(1,2,2)=temp006222(2,2,1)
       temp006222(1,2,3)=temp006322(2,2,1)
       temp006222(1,3,1)=temp006322(2,1,1)
       temp006222(1,3,2)=temp006322(2,2,1)
       temp006222(1,3,3)=temp006332(2,2,1)
       temp006222(2,1,2)=temp006222(2,2,1)
       temp006222(2,1,3)=temp006322(2,2,1)
       temp006222(2,2,3)=temp006322(2,2,2)
       temp006222(2,3,1)=temp006322(2,2,1)
       temp006222(2,3,2)=temp006322(2,2,2)
       temp006222(2,3,3)=temp006332(2,2,2)
       temp006311(1,1,2)=temp006321(1,1,1)
       temp006311(1,1,3)=temp006331(1,1,1)
       temp006311(1,2,1)=temp006321(1,1,1)
       temp006311(1,2,2)=temp006322(1,1,1)
       temp006311(1,2,3)=temp006332(1,1,1)
       temp006311(1,3,1)=temp006331(1,1,1)
       temp006311(1,3,2)=temp006332(1,1,1)
       temp006311(1,3,3)=temp006333(1,1,1)
       temp006321(1,1,2)=temp006322(1,1,1)
       temp006321(1,1,3)=temp006332(1,1,1)
       temp006321(1,2,1)=temp006322(1,1,1)
       temp006321(1,2,2)=temp006322(2,1,1)
       temp006321(1,2,3)=temp006332(2,1,1)
       temp006321(1,3,1)=temp006332(1,1,1)
       temp006321(1,3,2)=temp006332(2,1,1)
       temp006321(1,3,3)=temp006333(2,1,1)
       temp006322(1,1,2)=temp006322(2,1,1)
       temp006322(1,1,3)=temp006332(2,1,1)
       temp006322(1,2,1)=temp006322(2,1,1)
       temp006322(1,2,2)=temp006322(2,2,1)
       temp006322(1,2,3)=temp006332(2,2,1)
       temp006322(1,3,1)=temp006332(2,1,1)
       temp006322(1,3,2)=temp006332(2,2,1)
       temp006322(1,3,3)=temp006333(2,2,1)
       temp006322(2,1,2)=temp006322(2,2,1)
       temp006322(2,1,3)=temp006332(2,2,1)
       temp006322(2,2,3)=temp006332(2,2,2)
       temp006322(2,3,1)=temp006332(2,2,1)
       temp006322(2,3,2)=temp006332(2,2,2)
       temp006322(2,3,3)=temp006333(2,2,2)
       temp006331(1,1,2)=temp006332(1,1,1)
       temp006331(1,1,3)=temp006333(1,1,1)
       temp006331(1,2,1)=temp006332(1,1,1)
       temp006331(1,2,2)=temp006332(2,1,1)
       temp006331(1,2,3)=temp006333(2,1,1)
       temp006331(1,3,1)=temp006333(1,1,1)
       temp006331(1,3,2)=temp006333(2,1,1)
       temp006331(1,3,3)=temp006333(3,1,1)
       temp006332(1,1,2)=temp006332(2,1,1)
       temp006332(1,1,3)=temp006333(2,1,1)
       temp006332(1,2,1)=temp006332(2,1,1)
       temp006332(1,2,2)=temp006332(2,2,1)
       temp006332(1,2,3)=temp006333(2,2,1)
       temp006332(1,3,1)=temp006333(2,1,1)
       temp006332(1,3,2)=temp006333(2,2,1)
       temp006332(1,3,3)=temp006333(3,2,1)
       temp006332(2,1,2)=temp006332(2,2,1)
       temp006332(2,1,3)=temp006333(2,2,1)
       temp006332(2,2,3)=temp006333(2,2,2)
       temp006332(2,3,1)=temp006333(2,2,1)
       temp006332(2,3,2)=temp006333(2,2,2)
       temp006332(2,3,3)=temp006333(3,2,2)
       temp006333(1,1,2)=temp006333(2,1,1)
       temp006333(1,1,3)=temp006333(3,1,1)
       temp006333(1,2,1)=temp006333(2,1,1)
       temp006333(1,2,2)=temp006333(2,2,1)
       temp006333(1,2,3)=temp006333(3,2,1)
       temp006333(1,3,1)=temp006333(3,1,1)
       temp006333(1,3,2)=temp006333(3,2,1)
       temp006333(1,3,3)=temp006333(3,3,1)
       temp006333(2,1,2)=temp006333(2,2,1)
       temp006333(2,1,3)=temp006333(3,2,1)
       temp006333(2,2,3)=temp006333(3,2,2)
       temp006333(2,3,1)=temp006333(3,2,1)
       temp006333(2,3,2)=temp006333(3,2,2)
       temp006333(2,3,3)=temp006333(3,3,2)
       temp006333(3,1,2)=temp006333(3,2,1)
       temp006333(3,1,3)=temp006333(3,3,1)
       temp006333(3,2,3)=temp006333(3,3,2)
       temp71111(1,1,1)=IX*(aux71111(1,1,1)+det4*temp811111(1,1,jj)+14*t
     &  emp006111(1,1,1)*Z(jj,1))
       temp72111(1,1,1)=IX*(aux72111(1,1,1)+det4*temp821111(1,1,jj)+12*t
     &  emp006211(1,1,1)*Z(jj,1)+2*temp006111(1,1,1)*Z(jj,2))
       temp72211(1,1,1)=IX*(aux72211(1,1,1)+det4*temp822111(1,1,jj)+10*t
     &  emp006221(1,1,1)*Z(jj,1)+4*temp006211(1,1,1)*Z(jj,2))
       temp72221(1,1,1)=IX*(aux72221(1,1,1)+det4*temp822211(1,1,jj)+8*te
     &  mp006222(1,1,1)*Z(jj,1)+6*temp006221(1,1,1)*Z(jj,2))
       temp72222(1,1,1)=IX*(aux72222(1,1,1)+det4*temp822221(1,1,jj)+6*te
     &  mp006222(2,1,1)*Z(jj,1)+8*temp006222(1,1,1)*Z(jj,2))
       temp72222(2,1,1)=IX*(aux72222(2,1,1)+det4*temp822222(1,1,jj)+4*te
     &  mp006222(2,2,1)*Z(jj,1)+10*temp006222(2,1,1)*Z(jj,2))
       temp72222(2,2,1)=IX*(aux72222(2,2,1)+det4*temp822222(2,1,jj)+2*te
     &  mp006222(2,2,2)*Z(jj,1)+12*temp006222(2,2,1)*Z(jj,2))
       temp72222(2,2,2)=IX*(aux72222(2,2,2)+det4*temp822222(2,2,jj)+14*t
     &  emp006222(2,2,2)*Z(jj,2))
       temp73111(1,1,1)=IX*(aux73111(1,1,1)+det4*temp831111(1,1,jj)+12*t
     &  emp006311(1,1,1)*Z(jj,1)+2*temp006111(1,1,1)*Z(jj,3))
       temp73211(1,1,1)=IX*(aux73211(1,1,1)+det4*temp832111(1,1,jj)+10*t
     &  emp006321(1,1,1)*Z(jj,1)+2*(temp006311(1,1,1)*Z(jj,2)+temp006211
     &  (1,1,1)*Z(jj,3)))
       temp73221(1,1,1)=IX*(aux73221(1,1,1)+det4*temp832211(1,1,jj)+8*te
     &  mp006322(1,1,1)*Z(jj,1)+4*temp006321(1,1,1)*Z(jj,2)+2*temp006221
     &  (1,1,1)*Z(jj,3))
       temp73222(1,1,1)=IX*(aux73222(1,1,1)+det4*temp832221(1,1,jj)+6*(t
     &  emp006322(2,1,1)*Z(jj,1)+temp006322(1,1,1)*Z(jj,2))+2*temp006222
     &  (1,1,1)*Z(jj,3))
       temp73222(2,1,1)=IX*(aux73222(2,1,1)+det4*temp832222(1,1,jj)+4*te
     &  mp006322(2,2,1)*Z(jj,1)+8*temp006322(2,1,1)*Z(jj,2)+2*temp006222
     &  (2,1,1)*Z(jj,3))
       temp73222(2,2,1)=IX*(aux73222(2,2,1)+det4*temp832222(2,1,jj)+10*t
     &  emp006322(2,2,1)*Z(jj,2)+2*(temp006322(2,2,2)*Z(jj,1)+temp006222
     &  (2,2,1)*Z(jj,3)))
       temp73222(2,2,2)=IX*(aux73222(2,2,2)+det4*temp832222(2,2,jj)+12*t
     &  emp006322(2,2,2)*Z(jj,2)+2*temp006222(2,2,2)*Z(jj,3))
       temp73311(1,1,1)=IX*(aux73311(1,1,1)+det4*temp833111(1,1,jj)+10*t
     &  emp006331(1,1,1)*Z(jj,1)+4*temp006311(1,1,1)*Z(jj,3))
       temp73321(1,1,1)=IX*(aux73321(1,1,1)+det4*temp833211(1,1,jj)+8*te
     &  mp006332(1,1,1)*Z(jj,1)+2*temp006331(1,1,1)*Z(jj,2)+4*temp006321
     &  (1,1,1)*Z(jj,3))
       temp73322(1,1,1)=IX*(aux73322(1,1,1)+det4*temp833221(1,1,jj)+6*te
     &  mp006332(2,1,1)*Z(jj,1)+4*(temp006332(1,1,1)*Z(jj,2)+temp006322(
     &  1,1,1)*Z(jj,3)))
       temp73322(2,1,1)=IX*(aux73322(2,1,1)+det4*temp833222(1,1,jj)+6*te
     &  mp006332(2,1,1)*Z(jj,2)+4*(temp006332(2,2,1)*Z(jj,1)+temp006322(
     &  2,1,1)*Z(jj,3)))
       temp73322(2,2,1)=IX*(aux73322(2,2,1)+det4*temp833222(2,1,jj)+2*te
     &  mp006332(2,2,2)*Z(jj,1)+8*temp006332(2,2,1)*Z(jj,2)+4*temp006322
     &  (2,2,1)*Z(jj,3))
       temp73322(2,2,2)=IX*(aux73322(2,2,2)+det4*temp833222(2,2,jj)+10*t
     &  emp006332(2,2,2)*Z(jj,2)+4*temp006322(2,2,2)*Z(jj,3))
       temp73331(1,1,1)=IX*(aux73331(1,1,1)+det4*temp833311(1,1,jj)+8*te
     &  mp006333(1,1,1)*Z(jj,1)+6*temp006331(1,1,1)*Z(jj,3))
       temp73332(1,1,1)=IX*(aux73332(1,1,1)+det4*temp833321(1,1,jj)+2*te
     &  mp006333(1,1,1)*Z(jj,2)+6*(temp006333(2,1,1)*Z(jj,1)+temp006332(
     &  1,1,1)*Z(jj,3)))
       temp73332(2,1,1)=IX*(aux73332(2,1,1)+det4*temp833322(1,1,jj)+4*(t
     &  emp006333(2,2,1)*Z(jj,1)+temp006333(2,1,1)*Z(jj,2))+6*temp006332
     &  (2,1,1)*Z(jj,3))
       temp73332(2,2,1)=IX*(aux73332(2,2,1)+det4*temp833322(2,1,jj)+2*te
     &  mp006333(2,2,2)*Z(jj,1)+6*(temp006333(2,2,1)*Z(jj,2)+temp006332(
     &  2,2,1)*Z(jj,3)))
       temp73332(2,2,2)=IX*(aux73332(2,2,2)+det4*temp833322(2,2,jj)+8*te
     &  mp006333(2,2,2)*Z(jj,2)+6*temp006332(2,2,2)*Z(jj,3))
       temp73333(1,1,1)=IX*(aux73333(1,1,1)+det4*temp833331(1,1,jj)+6*te
     &  mp006333(3,1,1)*Z(jj,1)+8*temp006333(1,1,1)*Z(jj,3))
       temp73333(2,1,1)=IX*(aux73333(2,1,1)+det4*temp833332(1,1,jj)+4*te
     &  mp006333(3,2,1)*Z(jj,1)+2*temp006333(3,1,1)*Z(jj,2)+8*temp006333
     &  (2,1,1)*Z(jj,3))
       temp73333(2,2,1)=IX*(aux73333(2,2,1)+det4*temp833332(2,1,jj)+2*te
     &  mp006333(3,2,2)*Z(jj,1)+4*temp006333(3,2,1)*Z(jj,2)+8*temp006333
     &  (2,2,1)*Z(jj,3))
       temp73333(2,2,2)=IX*(aux73333(2,2,2)+det4*temp833332(2,2,jj)+6*te
     &  mp006333(3,2,2)*Z(jj,2)+8*temp006333(2,2,2)*Z(jj,3))
       temp73333(3,1,1)=IX*(aux73333(3,1,1)+det4*temp833333(1,1,jj)+4*te
     &  mp006333(3,3,1)*Z(jj,1)+10*temp006333(3,1,1)*Z(jj,3))
       temp73333(3,2,1)=IX*(aux73333(3,2,1)+det4*temp833333(2,1,jj)+2*(t
     &  emp006333(3,3,2)*Z(jj,1)+temp006333(3,3,1)*Z(jj,2))+10*temp00633
     &  3(3,2,1)*Z(jj,3))
       temp73333(3,2,2)=IX*(aux73333(3,2,2)+det4*temp833333(2,2,jj)+4*te
     &  mp006333(3,3,2)*Z(jj,2)+10*temp006333(3,2,2)*Z(jj,3))
       temp73333(3,3,1)=IX*(aux73333(3,3,1)+det4*temp833333(3,1,jj)+2*te
     &  mp006333(3,3,3)*Z(jj,1)+12*temp006333(3,3,1)*Z(jj,3))
       temp73333(3,3,2)=IX*(aux73333(3,3,2)+det4*temp833333(3,2,jj)+2*te
     &  mp006333(3,3,3)*Z(jj,2)+12*temp006333(3,3,2)*Z(jj,3))
       temp73333(3,3,3)=IX*(aux73333(3,3,3)+det4*temp833333(3,3,jj)+14*t
     &  emp006333(3,3,3)*Z(jj,3))
       temp71111(1,1,2)=temp72111(1,1,1)
       temp71111(1,1,3)=temp73111(1,1,1)
       temp71111(1,2,1)=temp72111(1,1,1)
       temp71111(1,2,2)=temp72211(1,1,1)
       temp71111(1,2,3)=temp73211(1,1,1)
       temp71111(1,3,1)=temp73111(1,1,1)
       temp71111(1,3,2)=temp73211(1,1,1)
       temp71111(1,3,3)=temp73311(1,1,1)
       temp72111(1,1,2)=temp72211(1,1,1)
       temp72111(1,1,3)=temp73211(1,1,1)
       temp72111(1,2,1)=temp72211(1,1,1)
       temp72111(1,2,2)=temp72221(1,1,1)
       temp72111(1,2,3)=temp73221(1,1,1)
       temp72111(1,3,1)=temp73211(1,1,1)
       temp72111(1,3,2)=temp73221(1,1,1)
       temp72111(1,3,3)=temp73321(1,1,1)
       temp72211(1,1,2)=temp72221(1,1,1)
       temp72211(1,1,3)=temp73221(1,1,1)
       temp72211(1,2,1)=temp72221(1,1,1)
       temp72211(1,2,2)=temp72222(1,1,1)
       temp72211(1,2,3)=temp73222(1,1,1)
       temp72211(1,3,1)=temp73221(1,1,1)
       temp72211(1,3,2)=temp73222(1,1,1)
       temp72211(1,3,3)=temp73322(1,1,1)
       temp72221(1,1,2)=temp72222(1,1,1)
       temp72221(1,1,3)=temp73222(1,1,1)
       temp72221(1,2,1)=temp72222(1,1,1)
       temp72221(1,2,2)=temp72222(2,1,1)
       temp72221(1,2,3)=temp73222(2,1,1)
       temp72221(1,3,1)=temp73222(1,1,1)
       temp72221(1,3,2)=temp73222(2,1,1)
       temp72221(1,3,3)=temp73322(2,1,1)
       temp72222(1,1,2)=temp72222(2,1,1)
       temp72222(1,1,3)=temp73222(2,1,1)
       temp72222(1,2,1)=temp72222(2,1,1)
       temp72222(1,2,2)=temp72222(2,2,1)
       temp72222(1,2,3)=temp73222(2,2,1)
       temp72222(1,3,1)=temp73222(2,1,1)
       temp72222(1,3,2)=temp73222(2,2,1)
       temp72222(1,3,3)=temp73322(2,2,1)
       temp72222(2,1,2)=temp72222(2,2,1)
       temp72222(2,1,3)=temp73222(2,2,1)
       temp72222(2,2,3)=temp73222(2,2,2)
       temp72222(2,3,1)=temp73222(2,2,1)
       temp72222(2,3,2)=temp73222(2,2,2)
       temp72222(2,3,3)=temp73322(2,2,2)
       temp73111(1,1,2)=temp73211(1,1,1)
       temp73111(1,1,3)=temp73311(1,1,1)
       temp73111(1,2,1)=temp73211(1,1,1)
       temp73111(1,2,2)=temp73221(1,1,1)
       temp73111(1,2,3)=temp73321(1,1,1)
       temp73111(1,3,1)=temp73311(1,1,1)
       temp73111(1,3,2)=temp73321(1,1,1)
       temp73111(1,3,3)=temp73331(1,1,1)
       temp73211(1,1,2)=temp73221(1,1,1)
       temp73211(1,1,3)=temp73321(1,1,1)
       temp73211(1,2,1)=temp73221(1,1,1)
       temp73211(1,2,2)=temp73222(1,1,1)
       temp73211(1,2,3)=temp73322(1,1,1)
       temp73211(1,3,1)=temp73321(1,1,1)
       temp73211(1,3,2)=temp73322(1,1,1)
       temp73211(1,3,3)=temp73332(1,1,1)
       temp73221(1,1,2)=temp73222(1,1,1)
       temp73221(1,1,3)=temp73322(1,1,1)
       temp73221(1,2,1)=temp73222(1,1,1)
       temp73221(1,2,2)=temp73222(2,1,1)
       temp73221(1,2,3)=temp73322(2,1,1)
       temp73221(1,3,1)=temp73322(1,1,1)
       temp73221(1,3,2)=temp73322(2,1,1)
       temp73221(1,3,3)=temp73332(2,1,1)
       temp73222(1,1,2)=temp73222(2,1,1)
       temp73222(1,1,3)=temp73322(2,1,1)
       temp73222(1,2,1)=temp73222(2,1,1)
       temp73222(1,2,2)=temp73222(2,2,1)
       temp73222(1,2,3)=temp73322(2,2,1)
       temp73222(1,3,1)=temp73322(2,1,1)
       temp73222(1,3,2)=temp73322(2,2,1)
       temp73222(1,3,3)=temp73332(2,2,1)
       temp73222(2,1,2)=temp73222(2,2,1)
       temp73222(2,1,3)=temp73322(2,2,1)
       temp73222(2,2,3)=temp73322(2,2,2)
       temp73222(2,3,1)=temp73322(2,2,1)
       temp73222(2,3,2)=temp73322(2,2,2)
       temp73222(2,3,3)=temp73332(2,2,2)
       temp73311(1,1,2)=temp73321(1,1,1)
       temp73311(1,1,3)=temp73331(1,1,1)
       temp73311(1,2,1)=temp73321(1,1,1)
       temp73311(1,2,2)=temp73322(1,1,1)
       temp73311(1,2,3)=temp73332(1,1,1)
       temp73311(1,3,1)=temp73331(1,1,1)
       temp73311(1,3,2)=temp73332(1,1,1)
       temp73311(1,3,3)=temp73333(1,1,1)
       temp73321(1,1,2)=temp73322(1,1,1)
       temp73321(1,1,3)=temp73332(1,1,1)
       temp73321(1,2,1)=temp73322(1,1,1)
       temp73321(1,2,2)=temp73322(2,1,1)
       temp73321(1,2,3)=temp73332(2,1,1)
       temp73321(1,3,1)=temp73332(1,1,1)
       temp73321(1,3,2)=temp73332(2,1,1)
       temp73321(1,3,3)=temp73333(2,1,1)
       temp73322(1,1,2)=temp73322(2,1,1)
       temp73322(1,1,3)=temp73332(2,1,1)
       temp73322(1,2,1)=temp73322(2,1,1)
       temp73322(1,2,2)=temp73322(2,2,1)
       temp73322(1,2,3)=temp73332(2,2,1)
       temp73322(1,3,1)=temp73332(2,1,1)
       temp73322(1,3,2)=temp73332(2,2,1)
       temp73322(1,3,3)=temp73333(2,2,1)
       temp73322(2,1,2)=temp73322(2,2,1)
       temp73322(2,1,3)=temp73332(2,2,1)
       temp73322(2,2,3)=temp73332(2,2,2)
       temp73322(2,3,1)=temp73332(2,2,1)
       temp73322(2,3,2)=temp73332(2,2,2)
       temp73322(2,3,3)=temp73333(2,2,2)
       temp73331(1,1,2)=temp73332(1,1,1)
       temp73331(1,1,3)=temp73333(1,1,1)
       temp73331(1,2,1)=temp73332(1,1,1)
       temp73331(1,2,2)=temp73332(2,1,1)
       temp73331(1,2,3)=temp73333(2,1,1)
       temp73331(1,3,1)=temp73333(1,1,1)
       temp73331(1,3,2)=temp73333(2,1,1)
       temp73331(1,3,3)=temp73333(3,1,1)
       temp73332(1,1,2)=temp73332(2,1,1)
       temp73332(1,1,3)=temp73333(2,1,1)
       temp73332(1,2,1)=temp73332(2,1,1)
       temp73332(1,2,2)=temp73332(2,2,1)
       temp73332(1,2,3)=temp73333(2,2,1)
       temp73332(1,3,1)=temp73333(2,1,1)
       temp73332(1,3,2)=temp73333(2,2,1)
       temp73332(1,3,3)=temp73333(3,2,1)
       temp73332(2,1,2)=temp73332(2,2,1)
       temp73332(2,1,3)=temp73333(2,2,1)
       temp73332(2,2,3)=temp73333(2,2,2)
       temp73332(2,3,1)=temp73333(2,2,1)
       temp73332(2,3,2)=temp73333(2,2,2)
       temp73332(2,3,3)=temp73333(3,2,2)
       temp73333(1,1,2)=temp73333(2,1,1)
       temp73333(1,1,3)=temp73333(3,1,1)
       temp73333(1,2,1)=temp73333(2,1,1)
       temp73333(1,2,2)=temp73333(2,2,1)
       temp73333(1,2,3)=temp73333(3,2,1)
       temp73333(1,3,1)=temp73333(3,1,1)
       temp73333(1,3,2)=temp73333(3,2,1)
       temp73333(1,3,3)=temp73333(3,3,1)
       temp73333(2,1,2)=temp73333(2,2,1)
       temp73333(2,1,3)=temp73333(3,2,1)
       temp73333(2,2,3)=temp73333(3,2,2)
       temp73333(2,3,1)=temp73333(3,2,1)
       temp73333(2,3,2)=temp73333(3,2,2)
       temp73333(2,3,3)=temp73333(3,3,2)
       temp73333(3,1,2)=temp73333(3,2,1)
       temp73333(3,1,3)=temp73333(3,3,1)
       temp73333(3,2,3)=temp73333(3,3,2)
c                Step3
       temp0000001(1)=I16Z*(aux0000001(1)+2*tempD4000000*F(4)+F(5)*temp0
     &  0001(1)-det4*temp00003(1,k,l))
       temp0000001(2)=I16Z*(aux0000001(2)+2*tempD4000000*F(6)+F(5)*temp0
     &  0001(2)-det4*temp00003(2,k,l))
       temp0000001(3)=I16Z*(aux0000001(3)+2*tempD4000000*F(7)+F(5)*temp0
     &  0001(3)-det4*temp00003(3,k,l))
       temp00003(1,1,1)=I20Z*(aux00003(1,1,1)+6*F(4)*temp00002(1,1)+F(5)
     &  *temp003(1,1,1)-det4*temp00511(1,k,l)+24*temp0000001(1)*ZZ(k,1,l
     &  ,1))
       temp00003(2,1,1)=I20Z*(aux00003(2,1,1)+2*F(6)*temp00002(1,1)+4*F(
     &  4)*temp00002(2,1)+F(5)*temp003(2,1,1)-det4*temp00521(1,k,l)+8*(t
     &  emp0000001(2)*ZZ(k,1,l,1)+temp0000001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1
     &  ))))
       temp00003(2,2,1)=I20Z*(aux00003(2,2,1)+4*F(6)*temp00002(2,1)+2*F(
     &  4)*temp00002(2,2)+F(5)*temp003(2,2,1)-det4*temp00522(1,k,l)+8*(t
     &  emp0000001(2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp0000001(1)*ZZ(k,2,l,
     &  2)))
       temp00003(2,2,2)=I20Z*(aux00003(2,2,2)+6*F(6)*temp00002(2,2)+F(5)
     &  *temp003(2,2,2)-det4*temp00522(2,k,l)+24*temp0000001(2)*ZZ(k,2,l
     &  ,2))
       temp00003(3,1,1)=I20Z*(aux00003(3,1,1)+2*F(7)*temp00002(1,1)+4*F(
     &  4)*temp00002(3,1)+F(5)*temp003(3,1,1)-det4*temp00531(1,k,l)+8*(t
     &  emp0000001(3)*ZZ(k,1,l,1)+temp0000001(1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1
     &  ))))
       temp00003(3,2,1)=I20Z*(aux00003(3,2,1)+2*F(7)*temp00002(2,1)+2*F(
     &  6)*temp00002(3,1)+2*F(4)*temp00002(3,2)+F(5)*temp003(3,2,1)-det4
     &  *temp00532(1,k,l)+4*(temp0000001(3)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+te
     &  mp0000001(2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)))+4*temp0000001(1)*(ZZ(k,2
     &  ,l,3)+ZZ(k,3,l,2)))
       temp00003(3,2,2)=I20Z*(aux00003(3,2,2)+2*F(7)*temp00002(2,2)+4*F(
     &  6)*temp00002(3,2)+F(5)*temp003(3,2,2)-det4*temp00532(2,k,l)+8*(t
     &  emp0000001(3)*ZZ(k,2,l,2)+temp0000001(2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2
     &  ))))
       temp00003(3,3,1)=I20Z*(aux00003(3,3,1)+4*F(7)*temp00002(3,1)+2*F(
     &  4)*temp00002(3,3)+F(5)*temp003(3,3,1)-det4*temp00533(1,k,l)+8*(t
     &  emp0000001(3)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+temp0000001(1)*ZZ(k,3,l,
     &  3)))
       temp00003(3,3,2)=I20Z*(aux00003(3,3,2)+4*F(7)*temp00002(3,2)+2*F(
     &  6)*temp00002(3,3)+F(5)*temp003(3,3,2)-det4*temp00533(2,k,l)+8*(t
     &  emp0000001(3)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+temp0000001(2)*ZZ(k,3,l,
     &  3)))
       temp00003(3,3,3)=I20Z*(aux00003(3,3,3)+6*F(7)*temp00002(3,3)+F(5)
     &  *temp003(3,3,3)-det4*temp00533(3,k,l)+24*temp0000001(3)*ZZ(k,3,l
     &  ,3))
       temp00003(1,1,2)=temp00003(2,1,1)
       temp00003(1,1,3)=temp00003(3,1,1)
       temp00003(1,2,1)=temp00003(2,1,1)
       temp00003(1,2,2)=temp00003(2,2,1)
       temp00003(1,2,3)=temp00003(3,2,1)
       temp00003(1,3,1)=temp00003(3,1,1)
       temp00003(1,3,2)=temp00003(3,2,1)
       temp00003(1,3,3)=temp00003(3,3,1)
       temp00003(2,1,2)=temp00003(2,2,1)
       temp00003(2,1,3)=temp00003(3,2,1)
       temp00003(2,2,3)=temp00003(3,2,2)
       temp00003(2,3,1)=temp00003(3,2,1)
       temp00003(2,3,2)=temp00003(3,2,2)
       temp00003(2,3,3)=temp00003(3,3,2)
       temp00003(3,1,2)=temp00003(3,2,1)
       temp00003(3,1,3)=temp00003(3,3,1)
       temp00003(3,2,3)=temp00003(3,3,2)
       temp00511(1,1,1)=I24Z*(aux00511(1,1,1)+10*F(4)*temp0041(1,1,1)+F(
     &  5)*temp511(1,1,1)-det4*temp71111(1,k,l)+80*temp00003(1,1,1)*ZZ(k
     &  ,1,l,1))
       temp00521(1,1,1)=I24Z*(aux00521(1,1,1)+2*F(6)*temp0041(1,1,1)+8*F
     &  (4)*temp0042(1,1,1)+F(5)*temp521(1,1,1)-det4*temp72111(1,k,l)+48
     &  *temp00003(2,1,1)*ZZ(k,1,l,1)+16*temp00003(1,1,1)*(ZZ(k,1,l,2)+Z
     &  Z(k,2,l,1)))
       temp00522(1,1,1)=I24Z*(aux00522(1,1,1)+4*F(6)*temp0042(1,1,1)+6*F
     &  (4)*temp0042(2,1,1)+F(5)*temp522(1,1,1)-det4*temp72211(1,k,l)+24
     &  *(temp00003(2,2,1)*ZZ(k,1,l,1)+temp00003(2,1,1)*(ZZ(k,1,l,2)+ZZ(
     &  k,2,l,1)))+8*temp00003(1,1,1)*ZZ(k,2,l,2))
       temp00522(2,1,1)=I24Z*(aux00522(2,1,1)+6*F(6)*temp0042(2,1,1)+4*F
     &  (4)*temp0042(2,2,1)+F(5)*temp522(2,1,1)-det4*temp72221(1,k,l)+8*
     &  temp00003(2,2,2)*ZZ(k,1,l,1)+24*(temp00003(2,2,1)*(ZZ(k,1,l,2)+Z
     &  Z(k,2,l,1))+temp00003(2,1,1)*ZZ(k,2,l,2)))
       temp00522(2,2,1)=I24Z*(aux00522(2,2,1)+8*F(6)*temp0042(2,2,1)+2*F
     &  (4)*temp0042(2,2,2)+F(5)*temp522(2,2,1)-det4*temp72222(1,k,l)+16
     &  *temp00003(2,2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+48*temp00003(2,2,1)*
     &  ZZ(k,2,l,2))
       temp00522(2,2,2)=I24Z*(aux00522(2,2,2)+10*F(6)*temp0042(2,2,2)+F(
     &  5)*temp522(2,2,2)-det4*temp72222(2,k,l)+80*temp00003(2,2,2)*ZZ(k
     &  ,2,l,2))
       temp00531(1,1,1)=I24Z*(aux00531(1,1,1)+2*F(7)*temp0041(1,1,1)+8*F
     &  (4)*temp0043(1,1,1)+F(5)*temp531(1,1,1)-det4*temp73111(1,k,l)+48
     &  *temp00003(3,1,1)*ZZ(k,1,l,1)+16*temp00003(1,1,1)*(ZZ(k,1,l,3)+Z
     &  Z(k,3,l,1)))
       temp00532(1,1,1)=I24Z*(aux00532(1,1,1)+2*F(7)*temp0042(1,1,1)+2*F
     &  (6)*temp0043(1,1,1)+6*F(4)*temp0043(2,1,1)+F(5)*temp532(1,1,1)-d
     &  et4*temp73211(1,k,l)+24*temp00003(3,2,1)*ZZ(k,1,l,1)+12*(temp000
     &  03(3,1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00003(2,1,1)*(ZZ(k,1,l,3
     &  )+ZZ(k,3,l,1)))+4*temp00003(1,1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp00532(2,1,1)=I24Z*(aux00532(2,1,1)+2*F(7)*temp0042(2,1,1)+4*F
     &  (6)*temp0043(2,1,1)+4*F(4)*temp0043(2,2,1)+F(5)*temp532(2,1,1)-d
     &  et4*temp73221(1,k,l)+16*temp00003(3,2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1
     &  ))+8*(temp00003(3,2,2)*ZZ(k,1,l,1)+temp00003(3,1,1)*ZZ(k,2,l,2))
     &  +8*temp00003(2,2,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*temp00003(2,1,1)
     &  *(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp00532(2,2,1)=I24Z*(aux00532(2,2,1)+2*F(7)*temp0042(2,2,1)+6*F
     &  (6)*temp0043(2,2,1)+2*F(4)*temp0043(2,2,2)+F(5)*temp532(2,2,1)-d
     &  et4*temp73222(1,k,l)+24*temp00003(3,2,1)*ZZ(k,2,l,2)+4*temp00003
     &  (2,2,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+12*(temp00003(3,2,2)*(ZZ(k,1,l
     &  ,2)+ZZ(k,2,l,1))+temp00003(2,2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))))
       temp00532(2,2,2)=I24Z*(aux00532(2,2,2)+2*F(7)*temp0042(2,2,2)+8*F
     &  (6)*temp0043(2,2,2)+F(5)*temp532(2,2,2)-det4*temp73222(2,k,l)+48
     &  *temp00003(3,2,2)*ZZ(k,2,l,2)+16*temp00003(2,2,2)*(ZZ(k,2,l,3)+Z
     &  Z(k,3,l,2)))
       temp00533(1,1,1)=I24Z*(aux00533(1,1,1)+4*F(7)*temp0043(1,1,1)+6*F
     &  (4)*temp0043(3,1,1)+F(5)*temp533(1,1,1)-det4*temp73311(1,k,l)+24
     &  *(temp00003(3,3,1)*ZZ(k,1,l,1)+temp00003(3,1,1)*(ZZ(k,1,l,3)+ZZ(
     &  k,3,l,1)))+8*temp00003(1,1,1)*ZZ(k,3,l,3))
       temp00533(2,1,1)=I24Z*(aux00533(2,1,1)+4*F(7)*temp0043(2,1,1)+2*F
     &  (6)*temp0043(3,1,1)+4*F(4)*temp0043(3,2,1)+F(5)*temp533(2,1,1)-d
     &  et4*temp73321(1,k,l)+8*(temp00003(3,3,2)*ZZ(k,1,l,1)+temp00003(3
     &  ,3,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))+16*temp00003(3,2,1)*(ZZ(k,1,l,3
     &  )+ZZ(k,3,l,1))+8*temp00003(3,1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*te
     &  mp00003(2,1,1)*ZZ(k,3,l,3))
       temp00533(2,2,1)=I24Z*(aux00533(2,2,1)+4*F(7)*temp0043(2,2,1)+4*F
     &  (6)*temp0043(3,2,1)+2*F(4)*temp0043(3,2,2)+F(5)*temp533(2,2,1)-d
     &  et4*temp73322(1,k,l)+8*(temp00003(3,3,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1
     &  ))+temp00003(3,3,1)*ZZ(k,2,l,2))+8*temp00003(3,2,2)*(ZZ(k,1,l,3)
     &  +ZZ(k,3,l,1))+16*temp00003(3,2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*te
     &  mp00003(2,2,1)*ZZ(k,3,l,3))
       temp00533(2,2,2)=I24Z*(aux00533(2,2,2)+4*F(7)*temp0043(2,2,2)+6*F
     &  (6)*temp0043(3,2,2)+F(5)*temp533(2,2,2)-det4*temp73322(2,k,l)+24
     &  *(temp00003(3,3,2)*ZZ(k,2,l,2)+temp00003(3,2,2)*(ZZ(k,2,l,3)+ZZ(
     &  k,3,l,2)))+8*temp00003(2,2,2)*ZZ(k,3,l,3))
       temp00533(3,1,1)=I24Z*(aux00533(3,1,1)+6*F(7)*temp0043(3,1,1)+4*F
     &  (4)*temp0043(3,3,1)+F(5)*temp533(3,1,1)-det4*temp73331(1,k,l)+8*
     &  temp00003(3,3,3)*ZZ(k,1,l,1)+24*(temp00003(3,3,1)*(ZZ(k,1,l,3)+Z
     &  Z(k,3,l,1))+temp00003(3,1,1)*ZZ(k,3,l,3)))
       temp00533(3,2,1)=I24Z*(aux00533(3,2,1)+6*F(7)*temp0043(3,2,1)+2*F
     &  (6)*temp0043(3,3,1)+2*F(4)*temp0043(3,3,2)+F(5)*temp533(3,2,1)-d
     &  et4*temp73332(1,k,l)+4*temp00003(3,3,3)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)
     &  )+12*(temp00003(3,3,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+temp00003(3,3,1
     &  )*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))+24*temp00003(3,2,1)*ZZ(k,3,l,3))
       temp00533(3,2,2)=I24Z*(aux00533(3,2,2)+6*F(7)*temp0043(3,2,2)+4*F
     &  (6)*temp0043(3,3,2)+F(5)*temp533(3,2,2)-det4*temp73332(2,k,l)+8*
     &  temp00003(3,3,3)*ZZ(k,2,l,2)+24*(temp00003(3,3,2)*(ZZ(k,2,l,3)+Z
     &  Z(k,3,l,2))+temp00003(3,2,2)*ZZ(k,3,l,3)))
       temp00533(3,3,1)=I24Z*(aux00533(3,3,1)+8*F(7)*temp0043(3,3,1)+2*F
     &  (4)*temp0043(3,3,3)+F(5)*temp533(3,3,1)-det4*temp73333(1,k,l)+16
     &  *temp00003(3,3,3)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+48*temp00003(3,3,1)*
     &  ZZ(k,3,l,3))
       temp00533(3,3,2)=I24Z*(aux00533(3,3,2)+8*F(7)*temp0043(3,3,2)+2*F
     &  (6)*temp0043(3,3,3)+F(5)*temp533(3,3,2)-det4*temp73333(2,k,l)+16
     &  *temp00003(3,3,3)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+48*temp00003(3,3,2)*
     &  ZZ(k,3,l,3))
       temp00533(3,3,3)=I24Z*(aux00533(3,3,3)+10*F(7)*temp0043(3,3,3)+F(
     &  5)*temp533(3,3,3)-det4*temp73333(3,k,l)+80*temp00003(3,3,3)*ZZ(k
     &  ,3,l,3))
       temp00511(1,1,2)=temp00521(1,1,1)
       temp00511(1,1,3)=temp00531(1,1,1)
       temp00511(1,2,1)=temp00521(1,1,1)
       temp00511(1,2,2)=temp00522(1,1,1)
       temp00511(1,2,3)=temp00532(1,1,1)
       temp00511(1,3,1)=temp00531(1,1,1)
       temp00511(1,3,2)=temp00532(1,1,1)
       temp00511(1,3,3)=temp00533(1,1,1)
       temp00521(1,1,2)=temp00522(1,1,1)
       temp00521(1,1,3)=temp00532(1,1,1)
       temp00521(1,2,1)=temp00522(1,1,1)
       temp00521(1,2,2)=temp00522(2,1,1)
       temp00521(1,2,3)=temp00532(2,1,1)
       temp00521(1,3,1)=temp00532(1,1,1)
       temp00521(1,3,2)=temp00532(2,1,1)
       temp00521(1,3,3)=temp00533(2,1,1)
       temp00522(1,1,2)=temp00522(2,1,1)
       temp00522(1,1,3)=temp00532(2,1,1)
       temp00522(1,2,1)=temp00522(2,1,1)
       temp00522(1,2,2)=temp00522(2,2,1)
       temp00522(1,2,3)=temp00532(2,2,1)
       temp00522(1,3,1)=temp00532(2,1,1)
       temp00522(1,3,2)=temp00532(2,2,1)
       temp00522(1,3,3)=temp00533(2,2,1)
       temp00522(2,1,2)=temp00522(2,2,1)
       temp00522(2,1,3)=temp00532(2,2,1)
       temp00522(2,2,3)=temp00532(2,2,2)
       temp00522(2,3,1)=temp00532(2,2,1)
       temp00522(2,3,2)=temp00532(2,2,2)
       temp00522(2,3,3)=temp00533(2,2,2)
       temp00531(1,1,2)=temp00532(1,1,1)
       temp00531(1,1,3)=temp00533(1,1,1)
       temp00531(1,2,1)=temp00532(1,1,1)
       temp00531(1,2,2)=temp00532(2,1,1)
       temp00531(1,2,3)=temp00533(2,1,1)
       temp00531(1,3,1)=temp00533(1,1,1)
       temp00531(1,3,2)=temp00533(2,1,1)
       temp00531(1,3,3)=temp00533(3,1,1)
       temp00532(1,1,2)=temp00532(2,1,1)
       temp00532(1,1,3)=temp00533(2,1,1)
       temp00532(1,2,1)=temp00532(2,1,1)
       temp00532(1,2,2)=temp00532(2,2,1)
       temp00532(1,2,3)=temp00533(2,2,1)
       temp00532(1,3,1)=temp00533(2,1,1)
       temp00532(1,3,2)=temp00533(2,2,1)
       temp00532(1,3,3)=temp00533(3,2,1)
       temp00532(2,1,2)=temp00532(2,2,1)
       temp00532(2,1,3)=temp00533(2,2,1)
       temp00532(2,2,3)=temp00533(2,2,2)
       temp00532(2,3,1)=temp00533(2,2,1)
       temp00532(2,3,2)=temp00533(2,2,2)
       temp00532(2,3,3)=temp00533(3,2,2)
       temp00533(1,1,2)=temp00533(2,1,1)
       temp00533(1,1,3)=temp00533(3,1,1)
       temp00533(1,2,1)=temp00533(2,1,1)
       temp00533(1,2,2)=temp00533(2,2,1)
       temp00533(1,2,3)=temp00533(3,2,1)
       temp00533(1,3,1)=temp00533(3,1,1)
       temp00533(1,3,2)=temp00533(3,2,1)
       temp00533(1,3,3)=temp00533(3,3,1)
       temp00533(2,1,2)=temp00533(2,2,1)
       temp00533(2,1,3)=temp00533(3,2,1)
       temp00533(2,2,3)=temp00533(3,2,2)
       temp00533(2,3,1)=temp00533(3,2,1)
       temp00533(2,3,2)=temp00533(3,2,2)
       temp00533(2,3,3)=temp00533(3,3,2)
       temp00533(3,1,2)=temp00533(3,2,1)
       temp00533(3,1,3)=temp00533(3,3,1)
       temp00533(3,2,3)=temp00533(3,3,2)
       temp6111(1,1,1)=IX*(aux6111(1,1,1)+det4*temp71111(1,1,jj)+12*temp
     &  00511(1,1,1)*Z(jj,1))
       temp6211(1,1,1)=IX*(aux6211(1,1,1)+det4*temp72111(1,1,jj)+10*temp
     &  00521(1,1,1)*Z(jj,1)+2*temp00511(1,1,1)*Z(jj,2))
       temp6221(1,1,1)=IX*(aux6221(1,1,1)+det4*temp72211(1,1,jj)+8*temp0
     &  0522(1,1,1)*Z(jj,1)+4*temp00521(1,1,1)*Z(jj,2))
       temp6222(1,1,1)=IX*(aux6222(1,1,1)+det4*temp72221(1,1,jj)+6*(temp
     &  00522(2,1,1)*Z(jj,1)+temp00522(1,1,1)*Z(jj,2)))
       temp6222(2,1,1)=IX*(aux6222(2,1,1)+det4*temp72222(1,1,jj)+4*temp0
     &  0522(2,2,1)*Z(jj,1)+8*temp00522(2,1,1)*Z(jj,2))
       temp6222(2,2,1)=IX*(aux6222(2,2,1)+det4*temp72222(2,1,jj)+2*temp0
     &  0522(2,2,2)*Z(jj,1)+10*temp00522(2,2,1)*Z(jj,2))
       temp6222(2,2,2)=IX*(aux6222(2,2,2)+det4*temp72222(2,2,jj)+12*temp
     &  00522(2,2,2)*Z(jj,2))
       temp6311(1,1,1)=IX*(aux6311(1,1,1)+det4*temp73111(1,1,jj)+10*temp
     &  00531(1,1,1)*Z(jj,1)+2*temp00511(1,1,1)*Z(jj,3))
       temp6321(1,1,1)=IX*(aux6321(1,1,1)+det4*temp73211(1,1,jj)+8*temp0
     &  0532(1,1,1)*Z(jj,1)+2*(temp00531(1,1,1)*Z(jj,2)+temp00521(1,1,1)
     &  *Z(jj,3)))
       temp6322(1,1,1)=IX*(aux6322(1,1,1)+det4*temp73221(1,1,jj)+6*temp0
     &  0532(2,1,1)*Z(jj,1)+4*temp00532(1,1,1)*Z(jj,2)+2*temp00522(1,1,1
     &  )*Z(jj,3))
       temp6322(2,1,1)=IX*(aux6322(2,1,1)+det4*temp73222(1,1,jj)+4*temp0
     &  0532(2,2,1)*Z(jj,1)+6*temp00532(2,1,1)*Z(jj,2)+2*temp00522(2,1,1
     &  )*Z(jj,3))
       temp6322(2,2,1)=IX*(aux6322(2,2,1)+det4*temp73222(2,1,jj)+8*temp0
     &  0532(2,2,1)*Z(jj,2)+2*(temp00532(2,2,2)*Z(jj,1)+temp00522(2,2,1)
     &  *Z(jj,3)))
       temp6322(2,2,2)=IX*(aux6322(2,2,2)+det4*temp73222(2,2,jj)+10*temp
     &  00532(2,2,2)*Z(jj,2)+2*temp00522(2,2,2)*Z(jj,3))
       temp6331(1,1,1)=IX*(aux6331(1,1,1)+det4*temp73311(1,1,jj)+8*temp0
     &  0533(1,1,1)*Z(jj,1)+4*temp00531(1,1,1)*Z(jj,3))
       temp6332(1,1,1)=IX*(aux6332(1,1,1)+det4*temp73321(1,1,jj)+6*temp0
     &  0533(2,1,1)*Z(jj,1)+2*temp00533(1,1,1)*Z(jj,2)+4*temp00532(1,1,1
     &  )*Z(jj,3))
       temp6332(2,1,1)=IX*(aux6332(2,1,1)+det4*temp73322(1,1,jj)+4*(temp
     &  00533(2,2,1)*Z(jj,1)+temp00533(2,1,1)*Z(jj,2))+4*temp00532(2,1,1
     &  )*Z(jj,3))
       temp6332(2,2,1)=IX*(aux6332(2,2,1)+det4*temp73322(2,1,jj)+2*temp0
     &  0533(2,2,2)*Z(jj,1)+6*temp00533(2,2,1)*Z(jj,2)+4*temp00532(2,2,1
     &  )*Z(jj,3))
       temp6332(2,2,2)=IX*(aux6332(2,2,2)+det4*temp73322(2,2,jj)+8*temp0
     &  0533(2,2,2)*Z(jj,2)+4*temp00532(2,2,2)*Z(jj,3))
       temp6333(1,1,1)=IX*(aux6333(1,1,1)+det4*temp73331(1,1,jj)+6*(temp
     &  00533(3,1,1)*Z(jj,1)+temp00533(1,1,1)*Z(jj,3)))
       temp6333(2,1,1)=IX*(aux6333(2,1,1)+det4*temp73332(1,1,jj)+4*temp0
     &  0533(3,2,1)*Z(jj,1)+2*temp00533(3,1,1)*Z(jj,2)+6*temp00533(2,1,1
     &  )*Z(jj,3))
       temp6333(2,2,1)=IX*(aux6333(2,2,1)+det4*temp73332(2,1,jj)+2*temp0
     &  0533(3,2,2)*Z(jj,1)+4*temp00533(3,2,1)*Z(jj,2)+6*temp00533(2,2,1
     &  )*Z(jj,3))
       temp6333(2,2,2)=IX*(aux6333(2,2,2)+det4*temp73332(2,2,jj)+6*(temp
     &  00533(3,2,2)*Z(jj,2)+temp00533(2,2,2)*Z(jj,3)))
       temp6333(3,1,1)=IX*(aux6333(3,1,1)+det4*temp73333(1,1,jj)+4*temp0
     &  0533(3,3,1)*Z(jj,1)+8*temp00533(3,1,1)*Z(jj,3))
       temp6333(3,2,1)=IX*(aux6333(3,2,1)+det4*temp73333(2,1,jj)+2*(temp
     &  00533(3,3,2)*Z(jj,1)+temp00533(3,3,1)*Z(jj,2))+8*temp00533(3,2,1
     &  )*Z(jj,3))
       temp6333(3,2,2)=IX*(aux6333(3,2,2)+det4*temp73333(2,2,jj)+4*temp0
     &  0533(3,3,2)*Z(jj,2)+8*temp00533(3,2,2)*Z(jj,3))
       temp6333(3,3,1)=IX*(aux6333(3,3,1)+det4*temp73333(3,1,jj)+2*temp0
     &  0533(3,3,3)*Z(jj,1)+10*temp00533(3,3,1)*Z(jj,3))
       temp6333(3,3,2)=IX*(aux6333(3,3,2)+det4*temp73333(3,2,jj)+2*temp0
     &  0533(3,3,3)*Z(jj,2)+10*temp00533(3,3,2)*Z(jj,3))
       temp6333(3,3,3)=IX*(aux6333(3,3,3)+det4*temp73333(3,3,jj)+12*temp
     &  00533(3,3,3)*Z(jj,3))
       temp6111(1,1,2)=temp6211(1,1,1)
       temp6111(1,1,3)=temp6311(1,1,1)
       temp6111(1,2,1)=temp6211(1,1,1)
       temp6111(1,2,2)=temp6221(1,1,1)
       temp6111(1,2,3)=temp6321(1,1,1)
       temp6111(1,3,1)=temp6311(1,1,1)
       temp6111(1,3,2)=temp6321(1,1,1)
       temp6111(1,3,3)=temp6331(1,1,1)
       temp6211(1,1,2)=temp6221(1,1,1)
       temp6211(1,1,3)=temp6321(1,1,1)
       temp6211(1,2,1)=temp6221(1,1,1)
       temp6211(1,2,2)=temp6222(1,1,1)
       temp6211(1,2,3)=temp6322(1,1,1)
       temp6211(1,3,1)=temp6321(1,1,1)
       temp6211(1,3,2)=temp6322(1,1,1)
       temp6211(1,3,3)=temp6332(1,1,1)
       temp6221(1,1,2)=temp6222(1,1,1)
       temp6221(1,1,3)=temp6322(1,1,1)
       temp6221(1,2,1)=temp6222(1,1,1)
       temp6221(1,2,2)=temp6222(2,1,1)
       temp6221(1,2,3)=temp6322(2,1,1)
       temp6221(1,3,1)=temp6322(1,1,1)
       temp6221(1,3,2)=temp6322(2,1,1)
       temp6221(1,3,3)=temp6332(2,1,1)
       temp6222(1,1,2)=temp6222(2,1,1)
       temp6222(1,1,3)=temp6322(2,1,1)
       temp6222(1,2,1)=temp6222(2,1,1)
       temp6222(1,2,2)=temp6222(2,2,1)
       temp6222(1,2,3)=temp6322(2,2,1)
       temp6222(1,3,1)=temp6322(2,1,1)
       temp6222(1,3,2)=temp6322(2,2,1)
       temp6222(1,3,3)=temp6332(2,2,1)
       temp6222(2,1,2)=temp6222(2,2,1)
       temp6222(2,1,3)=temp6322(2,2,1)
       temp6222(2,2,3)=temp6322(2,2,2)
       temp6222(2,3,1)=temp6322(2,2,1)
       temp6222(2,3,2)=temp6322(2,2,2)
       temp6222(2,3,3)=temp6332(2,2,2)
       temp6311(1,1,2)=temp6321(1,1,1)
       temp6311(1,1,3)=temp6331(1,1,1)
       temp6311(1,2,1)=temp6321(1,1,1)
       temp6311(1,2,2)=temp6322(1,1,1)
       temp6311(1,2,3)=temp6332(1,1,1)
       temp6311(1,3,1)=temp6331(1,1,1)
       temp6311(1,3,2)=temp6332(1,1,1)
       temp6311(1,3,3)=temp6333(1,1,1)
       temp6321(1,1,2)=temp6322(1,1,1)
       temp6321(1,1,3)=temp6332(1,1,1)
       temp6321(1,2,1)=temp6322(1,1,1)
       temp6321(1,2,2)=temp6322(2,1,1)
       temp6321(1,2,3)=temp6332(2,1,1)
       temp6321(1,3,1)=temp6332(1,1,1)
       temp6321(1,3,2)=temp6332(2,1,1)
       temp6321(1,3,3)=temp6333(2,1,1)
       temp6322(1,1,2)=temp6322(2,1,1)
       temp6322(1,1,3)=temp6332(2,1,1)
       temp6322(1,2,1)=temp6322(2,1,1)
       temp6322(1,2,2)=temp6322(2,2,1)
       temp6322(1,2,3)=temp6332(2,2,1)
       temp6322(1,3,1)=temp6332(2,1,1)
       temp6322(1,3,2)=temp6332(2,2,1)
       temp6322(1,3,3)=temp6333(2,2,1)
       temp6322(2,1,2)=temp6322(2,2,1)
       temp6322(2,1,3)=temp6332(2,2,1)
       temp6322(2,2,3)=temp6332(2,2,2)
       temp6322(2,3,1)=temp6332(2,2,1)
       temp6322(2,3,2)=temp6332(2,2,2)
       temp6322(2,3,3)=temp6333(2,2,2)
       temp6331(1,1,2)=temp6332(1,1,1)
       temp6331(1,1,3)=temp6333(1,1,1)
       temp6331(1,2,1)=temp6332(1,1,1)
       temp6331(1,2,2)=temp6332(2,1,1)
       temp6331(1,2,3)=temp6333(2,1,1)
       temp6331(1,3,1)=temp6333(1,1,1)
       temp6331(1,3,2)=temp6333(2,1,1)
       temp6331(1,3,3)=temp6333(3,1,1)
       temp6332(1,1,2)=temp6332(2,1,1)
       temp6332(1,1,3)=temp6333(2,1,1)
       temp6332(1,2,1)=temp6332(2,1,1)
       temp6332(1,2,2)=temp6332(2,2,1)
       temp6332(1,2,3)=temp6333(2,2,1)
       temp6332(1,3,1)=temp6333(2,1,1)
       temp6332(1,3,2)=temp6333(2,2,1)
       temp6332(1,3,3)=temp6333(3,2,1)
       temp6332(2,1,2)=temp6332(2,2,1)
       temp6332(2,1,3)=temp6333(2,2,1)
       temp6332(2,2,3)=temp6333(2,2,2)
       temp6332(2,3,1)=temp6333(2,2,1)
       temp6332(2,3,2)=temp6333(2,2,2)
       temp6332(2,3,3)=temp6333(3,2,2)
       temp6333(1,1,2)=temp6333(2,1,1)
       temp6333(1,1,3)=temp6333(3,1,1)
       temp6333(1,2,1)=temp6333(2,1,1)
       temp6333(1,2,2)=temp6333(2,2,1)
       temp6333(1,2,3)=temp6333(3,2,1)
       temp6333(1,3,1)=temp6333(3,1,1)
       temp6333(1,3,2)=temp6333(3,2,1)
       temp6333(1,3,3)=temp6333(3,3,1)
       temp6333(2,1,2)=temp6333(2,2,1)
       temp6333(2,1,3)=temp6333(3,2,1)
       temp6333(2,2,3)=temp6333(3,2,2)
       temp6333(2,3,1)=temp6333(3,2,1)
       temp6333(2,3,2)=temp6333(3,2,2)
       temp6333(2,3,3)=temp6333(3,3,2)
       temp6333(3,1,2)=temp6333(3,2,1)
       temp6333(3,1,3)=temp6333(3,3,1)
       temp6333(3,2,3)=temp6333(3,3,2)
c                Step4
       tempD4000000=I12Z*(auxD4000000+tempD40000*F(5)-det4*temp00002(k,l
     &  ))
       temp00002(1,1)=I16Z*(aux00002(1,1)+4*F(4)*temp00001(1)+F(5)*temp0
     &  02(1,1)-det4*temp0041(1,k,l)+8*tempD4000000*ZZ(k,1,l,1))
       temp00002(2,1)=I16Z*(aux00002(2,1)+2*(F(6)*temp00001(1)+F(4)*temp
     &  00001(2))+F(5)*temp002(2,1)-det4*temp0042(1,k,l)+4*tempD4000000*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1)))
       temp00002(2,2)=I16Z*(aux00002(2,2)+4*F(6)*temp00001(2)+F(5)*temp0
     &  02(2,2)-det4*temp0042(2,k,l)+8*tempD4000000*ZZ(k,2,l,2))
       temp00002(3,1)=I16Z*(aux00002(3,1)+2*(F(7)*temp00001(1)+F(4)*temp
     &  00001(3))+F(5)*temp002(3,1)-det4*temp0043(1,k,l)+4*tempD4000000*
     &  (ZZ(k,1,l,3)+ZZ(k,3,l,1)))
       temp00002(3,2)=I16Z*(aux00002(3,2)+2*(F(7)*temp00001(2)+F(6)*temp
     &  00001(3))+F(5)*temp002(3,2)-det4*temp0043(2,k,l)+4*tempD4000000*
     &  (ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp00002(3,3)=I16Z*(aux00002(3,3)+4*F(7)*temp00001(3)+F(5)*temp0
     &  02(3,3)-det4*temp0043(3,k,l)+8*tempD4000000*ZZ(k,3,l,3))
       temp00002(1,2)=temp00002(2,1)
       temp00002(1,3)=temp00002(3,1)
       temp00002(2,3)=temp00002(3,2)
       temp0041(1,1,1)=I20Z*(aux0041(1,1,1)+8*F(4)*temp003(1,1,1)+F(5)*t
     &  emp41(1,1,1)-det4*temp6111(1,k,l)+48*temp00002(1,1)*ZZ(k,1,l,1))
       temp0042(1,1,1)=I20Z*(aux0042(1,1,1)+2*F(6)*temp003(1,1,1)+6*F(4)
     &  *temp003(2,1,1)+F(5)*temp42(1,1,1)-det4*temp6211(1,k,l)+24*temp0
     &  0002(2,1)*ZZ(k,1,l,1)+12*temp00002(1,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)
     &  ))
       temp0042(2,1,1)=I20Z*(aux0042(2,1,1)+4*F(6)*temp003(2,1,1)+4*F(4)
     &  *temp003(2,2,1)+F(5)*temp42(2,1,1)-det4*temp6221(1,k,l)+16*temp0
     &  0002(2,1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+8*(temp00002(2,2)*ZZ(k,1,l,1
     &  )+temp00002(1,1)*ZZ(k,2,l,2)))
       temp0042(2,2,1)=I20Z*(aux0042(2,2,1)+6*F(6)*temp003(2,2,1)+2*F(4)
     &  *temp003(2,2,2)+F(5)*temp42(2,2,1)-det4*temp6222(1,k,l)+12*temp0
     &  0002(2,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+24*temp00002(2,1)*ZZ(k,2,l,2
     &  ))
       temp0042(2,2,2)=I20Z*(aux0042(2,2,2)+8*F(6)*temp003(2,2,2)+F(5)*t
     &  emp42(2,2,2)-det4*temp6222(2,k,l)+48*temp00002(2,2)*ZZ(k,2,l,2))
       temp0043(1,1,1)=I20Z*(aux0043(1,1,1)+2*F(7)*temp003(1,1,1)+6*F(4)
     &  *temp003(3,1,1)+F(5)*temp43(1,1,1)-det4*temp6311(1,k,l)+24*temp0
     &  0002(3,1)*ZZ(k,1,l,1)+12*temp00002(1,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1)
     &  ))
       temp0043(2,1,1)=I20Z*(aux0043(2,1,1)+2*F(7)*temp003(2,1,1)+2*F(6)
     &  *temp003(3,1,1)+4*F(4)*temp003(3,2,1)+F(5)*temp43(2,1,1)-det4*te
     &  mp6321(1,k,l)+8*(temp00002(3,2)*ZZ(k,1,l,1)+temp00002(3,1)*(ZZ(k
     &  ,1,l,2)+ZZ(k,2,l,1)))+8*temp00002(2,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))
     &  +4*temp00002(1,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp0043(2,2,1)=I20Z*(aux0043(2,2,1)+2*F(7)*temp003(2,2,1)+4*F(6)
     &  *temp003(3,2,1)+2*F(4)*temp003(3,2,2)+F(5)*temp43(2,2,1)-det4*te
     &  mp6322(1,k,l)+8*(temp00002(3,2)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00
     &  002(3,1)*ZZ(k,2,l,2))+4*temp00002(2,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))
     &  +8*temp00002(2,1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp0043(2,2,2)=I20Z*(aux0043(2,2,2)+2*F(7)*temp003(2,2,2)+6*F(6)
     &  *temp003(3,2,2)+F(5)*temp43(2,2,2)-det4*temp6322(2,k,l)+24*temp0
     &  0002(3,2)*ZZ(k,2,l,2)+12*temp00002(2,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)
     &  ))
       temp0043(3,1,1)=I20Z*(aux0043(3,1,1)+4*F(7)*temp003(3,1,1)+4*F(4)
     &  *temp003(3,3,1)+F(5)*temp43(3,1,1)-det4*temp6331(1,k,l)+16*temp0
     &  0002(3,1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*(temp00002(3,3)*ZZ(k,1,l,1
     &  )+temp00002(1,1)*ZZ(k,3,l,3)))
       temp0043(3,2,1)=I20Z*(aux0043(3,2,1)+2*F(6)*temp003(3,3,1)+2*F(4)
     &  *temp003(3,3,2)+F(5)*temp43(3,2,1)-det4*temp6332(1,k,l)+4*(F(7)*
     &  temp003(3,2,1)+temp00002(3,3)*(ZZ(k,1,l,2)+ZZ(k,2,l,1)))+8*temp0
     &  0002(3,2)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+8*temp00002(3,1)*(ZZ(k,2,l,3
     &  )+ZZ(k,3,l,2))+8*temp00002(2,1)*ZZ(k,3,l,3))
       temp0043(3,2,2)=I20Z*(aux0043(3,2,2)+4*F(7)*temp003(3,2,2)+4*F(6)
     &  *temp003(3,3,2)+F(5)*temp43(3,2,2)-det4*temp6332(2,k,l)+16*temp0
     &  0002(3,2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+8*(temp00002(3,3)*ZZ(k,2,l,2
     &  )+temp00002(2,2)*ZZ(k,3,l,3)))
       temp0043(3,3,1)=I20Z*(aux0043(3,3,1)+6*F(7)*temp003(3,3,1)+2*F(4)
     &  *temp003(3,3,3)+F(5)*temp43(3,3,1)-det4*temp6333(1,k,l)+12*temp0
     &  0002(3,3)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))+24*temp00002(3,1)*ZZ(k,3,l,3
     &  ))
       temp0043(3,3,2)=I20Z*(aux0043(3,3,2)+6*F(7)*temp003(3,3,2)+2*F(6)
     &  *temp003(3,3,3)+F(5)*temp43(3,3,2)-det4*temp6333(2,k,l)+12*temp0
     &  0002(3,3)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))+24*temp00002(3,2)*ZZ(k,3,l,3
     &  ))
       temp0043(3,3,3)=I20Z*(aux0043(3,3,3)+8*F(7)*temp003(3,3,3)+F(5)*t
     &  emp43(3,3,3)-det4*temp6333(3,k,l)+48*temp00002(3,3)*ZZ(k,3,l,3))
       temp0041(1,1,2)=temp0042(1,1,1)
       temp0041(1,1,3)=temp0043(1,1,1)
       temp0041(1,2,1)=temp0042(1,1,1)
       temp0041(1,2,2)=temp0042(2,1,1)
       temp0041(1,2,3)=temp0043(2,1,1)
       temp0041(1,3,1)=temp0043(1,1,1)
       temp0041(1,3,2)=temp0043(2,1,1)
       temp0041(1,3,3)=temp0043(3,1,1)
       temp0042(1,1,2)=temp0042(2,1,1)
       temp0042(1,1,3)=temp0043(2,1,1)
       temp0042(1,2,1)=temp0042(2,1,1)
       temp0042(1,2,2)=temp0042(2,2,1)
       temp0042(1,2,3)=temp0043(2,2,1)
       temp0042(1,3,1)=temp0043(2,1,1)
       temp0042(1,3,2)=temp0043(2,2,1)
       temp0042(1,3,3)=temp0043(3,2,1)
       temp0042(2,1,2)=temp0042(2,2,1)
       temp0042(2,1,3)=temp0043(2,2,1)
       temp0042(2,2,3)=temp0043(2,2,2)
       temp0042(2,3,1)=temp0043(2,2,1)
       temp0042(2,3,2)=temp0043(2,2,2)
       temp0042(2,3,3)=temp0043(3,2,2)
       temp0043(1,1,2)=temp0043(2,1,1)
       temp0043(1,1,3)=temp0043(3,1,1)
       temp0043(1,2,1)=temp0043(2,1,1)
       temp0043(1,2,2)=temp0043(2,2,1)
       temp0043(1,2,3)=temp0043(3,2,1)
       temp0043(1,3,1)=temp0043(3,1,1)
       temp0043(1,3,2)=temp0043(3,2,1)
       temp0043(1,3,3)=temp0043(3,3,1)
       temp0043(2,1,2)=temp0043(2,2,1)
       temp0043(2,1,3)=temp0043(3,2,1)
       temp0043(2,2,3)=temp0043(3,2,2)
       temp0043(2,3,1)=temp0043(3,2,1)
       temp0043(2,3,2)=temp0043(3,2,2)
       temp0043(2,3,3)=temp0043(3,3,2)
       temp0043(3,1,2)=temp0043(3,2,1)
       temp0043(3,1,3)=temp0043(3,3,1)
       temp0043(3,2,3)=temp0043(3,3,2)
       temp511(1,1,1)=IX*(aux511(1,1,1)+det4*temp6111(1,1,jj)+10*temp004
     &  1(1,1,1)*Z(jj,1))
       temp521(1,1,1)=IX*(aux521(1,1,1)+det4*temp6211(1,1,jj)+8*temp0042
     &  (1,1,1)*Z(jj,1)+2*temp0041(1,1,1)*Z(jj,2))
       temp522(1,1,1)=IX*(aux522(1,1,1)+det4*temp6221(1,1,jj)+6*temp0042
     &  (2,1,1)*Z(jj,1)+4*temp0042(1,1,1)*Z(jj,2))
       temp522(2,1,1)=IX*(aux522(2,1,1)+det4*temp6222(1,1,jj)+4*temp0042
     &  (2,2,1)*Z(jj,1)+6*temp0042(2,1,1)*Z(jj,2))
       temp522(2,2,1)=IX*(aux522(2,2,1)+det4*temp6222(2,1,jj)+2*temp0042
     &  (2,2,2)*Z(jj,1)+8*temp0042(2,2,1)*Z(jj,2))
       temp522(2,2,2)=IX*(aux522(2,2,2)+det4*temp6222(2,2,jj)+10*temp004
     &  2(2,2,2)*Z(jj,2))
       temp531(1,1,1)=IX*(aux531(1,1,1)+det4*temp6311(1,1,jj)+8*temp0043
     &  (1,1,1)*Z(jj,1)+2*temp0041(1,1,1)*Z(jj,3))
       temp532(1,1,1)=IX*(aux532(1,1,1)+det4*temp6321(1,1,jj)+6*temp0043
     &  (2,1,1)*Z(jj,1)+2*(temp0043(1,1,1)*Z(jj,2)+temp0042(1,1,1)*Z(jj,
     &  3)))
       temp532(2,1,1)=IX*(aux532(2,1,1)+det4*temp6322(1,1,jj)+4*(temp004
     &  3(2,2,1)*Z(jj,1)+temp0043(2,1,1)*Z(jj,2))+2*temp0042(2,1,1)*Z(jj
     &  ,3))
       temp532(2,2,1)=IX*(aux532(2,2,1)+det4*temp6322(2,1,jj)+6*temp0043
     &  (2,2,1)*Z(jj,2)+2*(temp0043(2,2,2)*Z(jj,1)+temp0042(2,2,1)*Z(jj,
     &  3)))
       temp532(2,2,2)=IX*(aux532(2,2,2)+det4*temp6322(2,2,jj)+8*temp0043
     &  (2,2,2)*Z(jj,2)+2*temp0042(2,2,2)*Z(jj,3))
       temp533(1,1,1)=IX*(aux533(1,1,1)+det4*temp6331(1,1,jj)+6*temp0043
     &  (3,1,1)*Z(jj,1)+4*temp0043(1,1,1)*Z(jj,3))
       temp533(2,1,1)=IX*(aux533(2,1,1)+det4*temp6332(1,1,jj)+2*temp0043
     &  (3,1,1)*Z(jj,2)+4*(temp0043(3,2,1)*Z(jj,1)+temp0043(2,1,1)*Z(jj,
     &  3)))
       temp533(2,2,1)=IX*(aux533(2,2,1)+det4*temp6332(2,1,jj)+2*temp0043
     &  (3,2,2)*Z(jj,1)+4*(temp0043(3,2,1)*Z(jj,2)+temp0043(2,2,1)*Z(jj,
     &  3)))
       temp533(2,2,2)=IX*(aux533(2,2,2)+det4*temp6332(2,2,jj)+6*temp0043
     &  (3,2,2)*Z(jj,2)+4*temp0043(2,2,2)*Z(jj,3))
       temp533(3,1,1)=IX*(aux533(3,1,1)+det4*temp6333(1,1,jj)+4*temp0043
     &  (3,3,1)*Z(jj,1)+6*temp0043(3,1,1)*Z(jj,3))
       temp533(3,2,1)=IX*(aux533(3,2,1)+det4*temp6333(2,1,jj)+2*(temp004
     &  3(3,3,2)*Z(jj,1)+temp0043(3,3,1)*Z(jj,2))+6*temp0043(3,2,1)*Z(jj
     &  ,3))
       temp533(3,2,2)=IX*(aux533(3,2,2)+det4*temp6333(2,2,jj)+4*temp0043
     &  (3,3,2)*Z(jj,2)+6*temp0043(3,2,2)*Z(jj,3))
       temp533(3,3,1)=IX*(aux533(3,3,1)+det4*temp6333(3,1,jj)+2*temp0043
     &  (3,3,3)*Z(jj,1)+8*temp0043(3,3,1)*Z(jj,3))
       temp533(3,3,2)=IX*(aux533(3,3,2)+det4*temp6333(3,2,jj)+2*temp0043
     &  (3,3,3)*Z(jj,2)+8*temp0043(3,3,2)*Z(jj,3))
       temp533(3,3,3)=IX*(aux533(3,3,3)+det4*temp6333(3,3,jj)+10*temp004
     &  3(3,3,3)*Z(jj,3))
       temp511(1,1,2)=temp521(1,1,1)
       temp511(1,1,3)=temp531(1,1,1)
       temp511(1,2,1)=temp521(1,1,1)
       temp511(1,2,2)=temp522(1,1,1)
       temp511(1,2,3)=temp532(1,1,1)
       temp511(1,3,1)=temp531(1,1,1)
       temp511(1,3,2)=temp532(1,1,1)
       temp511(1,3,3)=temp533(1,1,1)
       temp521(1,1,2)=temp522(1,1,1)
       temp521(1,1,3)=temp532(1,1,1)
       temp521(1,2,1)=temp522(1,1,1)
       temp521(1,2,2)=temp522(2,1,1)
       temp521(1,2,3)=temp532(2,1,1)
       temp521(1,3,1)=temp532(1,1,1)
       temp521(1,3,2)=temp532(2,1,1)
       temp521(1,3,3)=temp533(2,1,1)
       temp522(1,1,2)=temp522(2,1,1)
       temp522(1,1,3)=temp532(2,1,1)
       temp522(1,2,1)=temp522(2,1,1)
       temp522(1,2,2)=temp522(2,2,1)
       temp522(1,2,3)=temp532(2,2,1)
       temp522(1,3,1)=temp532(2,1,1)
       temp522(1,3,2)=temp532(2,2,1)
       temp522(1,3,3)=temp533(2,2,1)
       temp522(2,1,2)=temp522(2,2,1)
       temp522(2,1,3)=temp532(2,2,1)
       temp522(2,2,3)=temp532(2,2,2)
       temp522(2,3,1)=temp532(2,2,1)
       temp522(2,3,2)=temp532(2,2,2)
       temp522(2,3,3)=temp533(2,2,2)
       temp531(1,1,2)=temp532(1,1,1)
       temp531(1,1,3)=temp533(1,1,1)
       temp531(1,2,1)=temp532(1,1,1)
       temp531(1,2,2)=temp532(2,1,1)
       temp531(1,2,3)=temp533(2,1,1)
       temp531(1,3,1)=temp533(1,1,1)
       temp531(1,3,2)=temp533(2,1,1)
       temp531(1,3,3)=temp533(3,1,1)
       temp532(1,1,2)=temp532(2,1,1)
       temp532(1,1,3)=temp533(2,1,1)
       temp532(1,2,1)=temp532(2,1,1)
       temp532(1,2,2)=temp532(2,2,1)
       temp532(1,2,3)=temp533(2,2,1)
       temp532(1,3,1)=temp533(2,1,1)
       temp532(1,3,2)=temp533(2,2,1)
       temp532(1,3,3)=temp533(3,2,1)
       temp532(2,1,2)=temp532(2,2,1)
       temp532(2,1,3)=temp533(2,2,1)
       temp532(2,2,3)=temp533(2,2,2)
       temp532(2,3,1)=temp533(2,2,1)
       temp532(2,3,2)=temp533(2,2,2)
       temp532(2,3,3)=temp533(3,2,2)
       temp533(1,1,2)=temp533(2,1,1)
       temp533(1,1,3)=temp533(3,1,1)
       temp533(1,2,1)=temp533(2,1,1)
       temp533(1,2,2)=temp533(2,2,1)
       temp533(1,2,3)=temp533(3,2,1)
       temp533(1,3,1)=temp533(3,1,1)
       temp533(1,3,2)=temp533(3,2,1)
       temp533(1,3,3)=temp533(3,3,1)
       temp533(2,1,2)=temp533(2,2,1)
       temp533(2,1,3)=temp533(3,2,1)
       temp533(2,2,3)=temp533(3,2,2)
       temp533(2,3,1)=temp533(3,2,1)
       temp533(2,3,2)=temp533(3,2,2)
       temp533(2,3,3)=temp533(3,3,2)
       temp533(3,1,2)=temp533(3,2,1)
       temp533(3,1,3)=temp533(3,3,1)
       temp533(3,2,3)=temp533(3,3,2)
c                Step5
       temp00001(1)=I12Z*(aux00001(1)+2*tempD40000*F(4)+F(5)*temp001(1)-
     &  det4*temp003(1,k,l))
       temp00001(2)=I12Z*(aux00001(2)+2*tempD40000*F(6)+F(5)*temp001(2)-
     &  det4*temp003(2,k,l))
       temp00001(3)=I12Z*(aux00001(3)+2*tempD40000*F(7)+F(5)*temp001(3)-
     &  det4*temp003(3,k,l))
       temp003(1,1,1)=I16Z*(aux003(1,1,1)+6*F(4)*temp002(1,1)+F(5)*temp3
     &  (1,1,1)-det4*temp511(1,k,l)+24*temp00001(1)*ZZ(k,1,l,1))
       temp003(2,1,1)=I16Z*(aux003(2,1,1)+2*F(6)*temp002(1,1)+4*F(4)*tem
     &  p002(2,1)+F(5)*temp3(2,1,1)-det4*temp521(1,k,l)+8*(temp00001(2)*
     &  ZZ(k,1,l,1)+temp00001(1)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))))
       temp003(2,2,1)=I16Z*(aux003(2,2,1)+4*F(6)*temp002(2,1)+2*F(4)*tem
     &  p002(2,2)+F(5)*temp3(2,2,1)-det4*temp522(1,k,l)+8*(temp00001(2)*
     &  (ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00001(1)*ZZ(k,2,l,2)))
       temp003(2,2,2)=I16Z*(aux003(2,2,2)+6*F(6)*temp002(2,2)+F(5)*temp3
     &  (2,2,2)-det4*temp522(2,k,l)+24*temp00001(2)*ZZ(k,2,l,2))
       temp003(3,1,1)=I16Z*(aux003(3,1,1)+2*F(7)*temp002(1,1)+4*F(4)*tem
     &  p002(3,1)+F(5)*temp3(3,1,1)-det4*temp531(1,k,l)+8*(temp00001(3)*
     &  ZZ(k,1,l,1)+temp00001(1)*(ZZ(k,1,l,3)+ZZ(k,3,l,1))))
       temp003(3,2,1)=I16Z*(aux003(3,2,1)+2*F(7)*temp002(2,1)+2*F(6)*tem
     &  p002(3,1)+2*F(4)*temp002(3,2)+F(5)*temp3(3,2,1)-det4*temp532(1,k
     &  ,l)+4*(temp00001(3)*(ZZ(k,1,l,2)+ZZ(k,2,l,1))+temp00001(2)*(ZZ(k
     &  ,1,l,3)+ZZ(k,3,l,1)))+4*temp00001(1)*(ZZ(k,2,l,3)+ZZ(k,3,l,2)))
       temp003(3,2,2)=I16Z*(aux003(3,2,2)+2*F(7)*temp002(2,2)+4*F(6)*tem
     &  p002(3,2)+F(5)*temp3(3,2,2)-det4*temp532(2,k,l)+8*(temp00001(3)*
     &  ZZ(k,2,l,2)+temp00001(2)*(ZZ(k,2,l,3)+ZZ(k,3,l,2))))
       temp003(3,3,1)=I16Z*(aux003(3,3,1)+4*F(7)*temp002(3,1)+2*F(4)*tem
     &  p002(3,3)+F(5)*temp3(3,3,1)-det4*temp533(1,k,l)+8*(temp00001(3)*
     &  (ZZ(k,1,l,3)+ZZ(k,3,l,1))+temp00001(1)*ZZ(k,3,l,3)))
       temp003(3,3,2)=I16Z*(aux003(3,3,2)+4*F(7)*temp002(3,2)+2*F(6)*tem
     &  p002(3,3)+F(5)*temp3(3,3,2)-det4*temp533(2,k,l)+8*(temp00001(3)*
     &  (ZZ(k,2,l,3)+ZZ(k,3,l,2))+temp00001(2)*ZZ(k,3,l,3)))
       temp003(3,3,3)=I16Z*(aux003(3,3,3)+6*F(7)*temp002(3,3)+F(5)*temp3
     &  (3,3,3)-det4*temp533(3,k,l)+24*temp00001(3)*ZZ(k,3,l,3))
       temp003(1,1,2)=temp003(2,1,1)
       temp003(1,1,3)=temp003(3,1,1)
       temp003(1,2,1)=temp003(2,1,1)
       temp003(1,2,2)=temp003(2,2,1)
       temp003(1,2,3)=temp003(3,2,1)
       temp003(1,3,1)=temp003(3,1,1)
       temp003(1,3,2)=temp003(3,2,1)
       temp003(1,3,3)=temp003(3,3,1)
       temp003(2,1,2)=temp003(2,2,1)
       temp003(2,1,3)=temp003(3,2,1)
       temp003(2,2,3)=temp003(3,2,2)
       temp003(2,3,1)=temp003(3,2,1)
       temp003(2,3,2)=temp003(3,2,2)
       temp003(2,3,3)=temp003(3,3,2)
       temp003(3,1,2)=temp003(3,2,1)
       temp003(3,1,3)=temp003(3,3,1)
       temp003(3,2,3)=temp003(3,3,2)
       temp41(1,1,1)=IX*(aux41(1,1,1)+det4*temp511(1,1,jj)+8*temp003(1,1
     &  ,1)*Z(jj,1))
       temp42(1,1,1)=IX*(aux42(1,1,1)+det4*temp521(1,1,jj)+6*temp003(2,1
     &  ,1)*Z(jj,1)+2*temp003(1,1,1)*Z(jj,2))
       temp42(2,1,1)=IX*(aux42(2,1,1)+det4*temp522(1,1,jj)+4*(temp003(2,
     &  2,1)*Z(jj,1)+temp003(2,1,1)*Z(jj,2)))
       temp42(2,2,1)=IX*(aux42(2,2,1)+det4*temp522(2,1,jj)+2*temp003(2,2
     &  ,2)*Z(jj,1)+6*temp003(2,2,1)*Z(jj,2))
       temp42(2,2,2)=IX*(aux42(2,2,2)+det4*temp522(2,2,jj)+8*temp003(2,2
     &  ,2)*Z(jj,2))
       temp43(1,1,1)=IX*(aux43(1,1,1)+det4*temp531(1,1,jj)+6*temp003(3,1
     &  ,1)*Z(jj,1)+2*temp003(1,1,1)*Z(jj,3))
       temp43(2,1,1)=IX*(aux43(2,1,1)+det4*temp532(1,1,jj)+4*temp003(3,2
     &  ,1)*Z(jj,1)+2*(temp003(3,1,1)*Z(jj,2)+temp003(2,1,1)*Z(jj,3)))
       temp43(2,2,1)=IX*(aux43(2,2,1)+det4*temp532(2,1,jj)+4*temp003(3,2
     &  ,1)*Z(jj,2)+2*(temp003(3,2,2)*Z(jj,1)+temp003(2,2,1)*Z(jj,3)))
       temp43(2,2,2)=IX*(aux43(2,2,2)+det4*temp532(2,2,jj)+6*temp003(3,2
     &  ,2)*Z(jj,2)+2*temp003(2,2,2)*Z(jj,3))
       temp43(3,1,1)=IX*(aux43(3,1,1)+det4*temp533(1,1,jj)+4*(temp003(3,
     &  3,1)*Z(jj,1)+temp003(3,1,1)*Z(jj,3)))
       temp43(3,2,1)=IX*(aux43(3,2,1)+det4*temp533(2,1,jj)+2*(temp003(3,
     &  3,2)*Z(jj,1)+temp003(3,3,1)*Z(jj,2))+4*temp003(3,2,1)*Z(jj,3))
       temp43(3,2,2)=IX*(aux43(3,2,2)+det4*temp533(2,2,jj)+4*(temp003(3,
     &  3,2)*Z(jj,2)+temp003(3,2,2)*Z(jj,3)))
       temp43(3,3,1)=IX*(aux43(3,3,1)+det4*temp533(3,1,jj)+2*temp003(3,3
     &  ,3)*Z(jj,1)+6*temp003(3,3,1)*Z(jj,3))
       temp43(3,3,2)=IX*(aux43(3,3,2)+det4*temp533(3,2,jj)+2*temp003(3,3
     &  ,3)*Z(jj,2)+6*temp003(3,3,2)*Z(jj,3))
       temp43(3,3,3)=IX*(aux43(3,3,3)+det4*temp533(3,3,jj)+8*temp003(3,3
     &  ,3)*Z(jj,3))
       temp41(1,1,2)=temp42(1,1,1)
       temp41(1,1,3)=temp43(1,1,1)
       temp41(1,2,1)=temp42(1,1,1)
       temp41(1,2,2)=temp42(2,1,1)
       temp41(1,2,3)=temp43(2,1,1)
       temp41(1,3,1)=temp43(1,1,1)
       temp41(1,3,2)=temp43(2,1,1)
       temp41(1,3,3)=temp43(3,1,1)
       temp42(1,1,2)=temp42(2,1,1)
       temp42(1,1,3)=temp43(2,1,1)
       temp42(1,2,1)=temp42(2,1,1)
       temp42(1,2,2)=temp42(2,2,1)
       temp42(1,2,3)=temp43(2,2,1)
       temp42(1,3,1)=temp43(2,1,1)
       temp42(1,3,2)=temp43(2,2,1)
       temp42(1,3,3)=temp43(3,2,1)
       temp42(2,1,2)=temp42(2,2,1)
       temp42(2,1,3)=temp43(2,2,1)
       temp42(2,2,3)=temp43(2,2,2)
       temp42(2,3,1)=temp43(2,2,1)
       temp42(2,3,2)=temp43(2,2,2)
       temp42(2,3,3)=temp43(3,2,2)
       temp43(1,1,2)=temp43(2,1,1)
       temp43(1,1,3)=temp43(3,1,1)
       temp43(1,2,1)=temp43(2,1,1)
       temp43(1,2,2)=temp43(2,2,1)
       temp43(1,2,3)=temp43(3,2,1)
       temp43(1,3,1)=temp43(3,1,1)
       temp43(1,3,2)=temp43(3,2,1)
       temp43(1,3,3)=temp43(3,3,1)
       temp43(2,1,2)=temp43(2,2,1)
       temp43(2,1,3)=temp43(3,2,1)
       temp43(2,2,3)=temp43(3,2,2)
       temp43(2,3,1)=temp43(3,2,1)
       temp43(2,3,2)=temp43(3,2,2)
       temp43(2,3,3)=temp43(3,3,2)
       temp43(3,1,2)=temp43(3,2,1)
       temp43(3,1,3)=temp43(3,3,1)
       temp43(3,2,3)=temp43(3,3,2)
c                Step6
       tempD40000=I8Z*(auxD40000+tempD400*F(5)-det4*temp002(k,l))
       temp002(1,1)=I12Z*(aux002(1,1)+4*F(4)*temp001(1)+F(5)*temp2(1,1)-
     &  det4*temp41(1,k,l)+8*tempD40000*ZZ(k,1,l,1))
       temp002(2,1)=I12Z*(aux002(2,1)+2*(F(6)*temp001(1)+F(4)*temp001(2)
     &  )+F(5)*temp2(2,1)-det4*temp42(1,k,l)+4*tempD40000*(ZZ(k,1,l,2)+Z
     &  Z(k,2,l,1)))
       temp002(2,2)=I12Z*(aux002(2,2)+4*F(6)*temp001(2)+F(5)*temp2(2,2)-
     &  det4*temp42(2,k,l)+8*tempD40000*ZZ(k,2,l,2))
       temp002(3,1)=I12Z*(aux002(3,1)+2*(F(7)*temp001(1)+F(4)*temp001(3)
     &  )+F(5)*temp2(3,1)-det4*temp43(1,k,l)+4*tempD40000*(ZZ(k,1,l,3)+Z
     &  Z(k,3,l,1)))
       temp002(3,2)=I12Z*(aux002(3,2)+2*(F(7)*temp001(2)+F(6)*temp001(3)
     &  )+F(5)*temp2(3,2)-det4*temp43(2,k,l)+4*tempD40000*(ZZ(k,2,l,3)+Z
     &  Z(k,3,l,2)))
       temp002(3,3)=I12Z*(aux002(3,3)+4*F(7)*temp001(3)+F(5)*temp2(3,3)-
     &  det4*temp43(3,k,l)+8*tempD40000*ZZ(k,3,l,3))
       temp002(1,2)=temp002(2,1)
       temp002(1,3)=temp002(3,1)
       temp002(2,3)=temp002(3,2)
       temp3(1,1,1)=IX*(aux3(1,1,1)+det4*temp41(1,1,jj)+6*temp002(1,1)*Z
     &  (jj,1))
       temp3(2,1,1)=IX*(aux3(2,1,1)+det4*temp42(1,1,jj)+4*temp002(2,1)*Z
     &  (jj,1)+2*temp002(1,1)*Z(jj,2))
       temp3(2,2,1)=IX*(aux3(2,2,1)+det4*temp42(2,1,jj)+2*temp002(2,2)*Z
     &  (jj,1)+4*temp002(2,1)*Z(jj,2))
       temp3(2,2,2)=IX*(aux3(2,2,2)+det4*temp42(2,2,jj)+6*temp002(2,2)*Z
     &  (jj,2))
       temp3(3,1,1)=IX*(aux3(3,1,1)+det4*temp43(1,1,jj)+4*temp002(3,1)*Z
     &  (jj,1)+2*temp002(1,1)*Z(jj,3))
       temp3(3,2,1)=IX*(aux3(3,2,1)+det4*temp43(2,1,jj)+2*(temp002(3,2)*
     &  Z(jj,1)+temp002(3,1)*Z(jj,2))+2*temp002(2,1)*Z(jj,3))
       temp3(3,2,2)=IX*(aux3(3,2,2)+det4*temp43(2,2,jj)+4*temp002(3,2)*Z
     &  (jj,2)+2*temp002(2,2)*Z(jj,3))
       temp3(3,3,1)=IX*(aux3(3,3,1)+det4*temp43(3,1,jj)+2*temp002(3,3)*Z
     &  (jj,1)+4*temp002(3,1)*Z(jj,3))
       temp3(3,3,2)=IX*(aux3(3,3,2)+det4*temp43(3,2,jj)+2*temp002(3,3)*Z
     &  (jj,2)+4*temp002(3,2)*Z(jj,3))
       temp3(3,3,3)=IX*(aux3(3,3,3)+det4*temp43(3,3,jj)+6*temp002(3,3)*Z
     &  (jj,3))
       temp3(1,1,2)=temp3(2,1,1)
       temp3(1,1,3)=temp3(3,1,1)
       temp3(1,2,1)=temp3(2,1,1)
       temp3(1,2,2)=temp3(2,2,1)
       temp3(1,2,3)=temp3(3,2,1)
       temp3(1,3,1)=temp3(3,1,1)
       temp3(1,3,2)=temp3(3,2,1)
       temp3(1,3,3)=temp3(3,3,1)
       temp3(2,1,2)=temp3(2,2,1)
       temp3(2,1,3)=temp3(3,2,1)
       temp3(2,2,3)=temp3(3,2,2)
       temp3(2,3,1)=temp3(3,2,1)
       temp3(2,3,2)=temp3(3,2,2)
       temp3(2,3,3)=temp3(3,3,2)
       temp3(3,1,2)=temp3(3,2,1)
       temp3(3,1,3)=temp3(3,3,1)
       temp3(3,2,3)=temp3(3,3,2)
c                Step7
       temp001(1)=I8Z*(aux001(1)+2*tempD400*F(4)+F(5)*temp1(1)-det4*temp
     &  3(1,k,l))
       temp001(2)=I8Z*(aux001(2)+2*tempD400*F(6)+F(5)*temp1(2)-det4*temp
     &  3(2,k,l))
       temp001(3)=I8Z*(aux001(3)+2*tempD400*F(7)+F(5)*temp1(3)-det4*temp
     &  3(3,k,l))
       temp2(1,1)=IX*(aux2(1,1)+det4*temp3(1,1,jj)+4*temp001(1)*Z(jj,1))
       temp2(2,1)=IX*(aux2(2,1)+det4*temp3(2,1,jj)+2*temp001(2)*Z(jj,1)+
     &  2*temp001(1)*Z(jj,2))
       temp2(2,2)=IX*(aux2(2,2)+det4*temp3(2,2,jj)+4*temp001(2)*Z(jj,2))
       temp2(3,1)=IX*(aux2(3,1)+det4*temp3(3,1,jj)+2*temp001(3)*Z(jj,1)+
     &  2*temp001(1)*Z(jj,3))
       temp2(3,2)=IX*(aux2(3,2)+det4*temp3(3,2,jj)+2*temp001(3)*Z(jj,2)+
     &  2*temp001(2)*Z(jj,3))
       temp2(3,3)=IX*(aux2(3,3)+det4*temp3(3,3,jj)+4*temp001(3)*Z(jj,3))
       temp2(1,2)=temp2(2,1)
       temp2(1,3)=temp2(3,1)
       temp2(2,3)=temp2(3,2)
c                Step8
       tempD400=I4Z*(auxD400+tempD40*F(5)-det4*temp2(k,l))
       temp1(1)=IX*(aux1(1)+det4*temp2(1,jj)+2*tempD400*Z(jj,1))
       temp1(2)=IX*(aux1(2)+det4*temp2(2,jj)+2*tempD400*Z(jj,2))
       temp1(3)=IX*(aux1(3)+det4*temp2(3,jj)+2*tempD400*Z(jj,3))
c                Step9
       tempD40=IX*(auxD40+det4*temp1(jj))
        ac=3
        accuracyDR(1,0,ac)  = abs(D0        /tempD40         -1d0)
        accuracyDR(1,1,ac)  = abs(Dij(1,1)  /temp1(1)        -1d0)
        accuracyDR(2,1,ac)  = abs(Dij(2,1)  /temp1(2)        -1d0)
        accuracyDR(3,1,ac)  = abs(Dij(3,1)  /temp1(3)        -1d0)
        accuracyDR(1,2,ac)  = abs(Dij(1,2)  /temp2(1,1)      -1d0)
        accuracyDR(2,2,ac)  = abs(Dij(2,2)  /temp2(2,2)      -1d0)
        accuracyDR(3,2,ac)  = abs(Dij(3,2)  /temp2(3,3)      -1d0)
        accuracyDR(4,2,ac)  = abs(Dij(4,2)  /temp2(2,1)      -1d0)
        accuracyDR(5,2,ac)  = abs(Dij(5,2)  /temp2(3,1)      -1d0)
        accuracyDR(6,2,ac)  = abs(Dij(6,2)  /temp2(3,2)      -1d0)
        accuracyDR(7,2,ac)  = abs(Dij(7,2)  /tempD400        -1d0)
        accuracyDR(1,3,ac)  = abs(Dij(1,3)  /temp3(1,1,1)    -1d0)
        accuracyDR(2,3,ac)  = abs(Dij(2,3)  /temp3(2,2,2)    -1d0)
        accuracyDR(3,3,ac)  = abs(Dij(3,3)  /temp3(3,3,3)    -1d0)
        accuracyDR(4,3,ac)  = abs(Dij(4,3)  /temp3(2,1,1)    -1d0)
        accuracyDR(5,3,ac)  = abs(Dij(5,3)  /temp3(3,1,1)    -1d0)
        accuracyDR(6,3,ac)  = abs(Dij(6,3)  /temp3(2,2,1)    -1d0)
        accuracyDR(7,3,ac)  = abs(Dij(7,3)  /temp3(3,3,1)    -1d0)
        accuracyDR(8,3,ac)  = abs(Dij(8,3)  /temp3(3,2,2)    -1d0)
        accuracyDR(9,3,ac)  = abs(Dij(9,3)  /temp3(3,3,2)    -1d0)
        accuracyDR(10,3,ac) = abs(Dij(10,3) /temp3(3,2,1)    -1d0)
        accuracyDR(11,3,ac) = abs(Dij(11,3) /temp001(1)      -1d0)
        accuracyDR(12,3,ac) = abs(Dij(12,3) /temp001(2)      -1d0)
        accuracyDR(13,3,ac) = abs(Dij(13,3) /temp001(3)      -1d0)
        accuracyDR(7,1,ac)  = abs(Dij(7,1)  /temp002(1,1)    -1d0)
        accuracyDR(8,1,ac)  = abs(Dij(8,1)  /temp002(2,2)    -1d0)
        accuracyDR(9,1,ac)  = abs(Dij(9,1)  /temp002(3,3)    -1d0)
        accuracyDR(10,1,ac) = abs(Dij(10,1) /temp002(2,1)    -1d0)
        accuracyDR(11,1,ac) = abs(Dij(11,1) /temp002(3,1)    -1d0)
        accuracyDR(12,1,ac) = abs(Dij(12,1) /temp002(3,2)    -1d0)
        accuracyDR(13,1,ac) = abs(Dij(13,1) /tempD40000      -1d0)
        accuracyDR(1,4,ac)  = abs(Dij(1,4)  /temp41(1,1,1)   -1d0)
        accuracyDR(2,4,ac)  = abs(Dij(2,4)  /temp42(2,2,2)   -1d0)
        accuracyDR(3,4,ac)  = abs(Dij(3,4)  /temp43(3,3,3)   -1d0)
        accuracyDR(4,4,ac)  = abs(Dij(4,4)  /temp42(1,1,1)   -1d0)
        accuracyDR(5,4,ac)  = abs(Dij(5,4)  /temp43(1,1,1)   -1d0)
        accuracyDR(6,4,ac)  = abs(Dij(6,4)  /temp42(2,1,1)   -1d0)
        accuracyDR(7,4,ac)  = abs(Dij(7,4)  /temp43(2,1,1)   -1d0)
        accuracyDR(8,4,ac)  = abs(Dij(8,4)  /temp43(3,1,1)   -1d0)
        accuracyDR(9,4,ac)  = abs(Dij(9,4)  /temp42(2,2,1)   -1d0)
        accuracyDR(10,4,ac) = abs(Dij(10,4) /temp43(2,2,1)   -1d0)
        accuracyDR(11,4,ac) = abs(Dij(11,4) /temp43(3,2,1)   -1d0)
        accuracyDR(12,4,ac) = abs(Dij(12,4) /temp43(3,3,1)   -1d0)
        accuracyDR(13,4,ac) = abs(Dij(13,4) /temp43(2,2,2)   -1d0)
        accuracyDR(14,4,ac) = abs(Dij(14,4) /temp43(3,2,2)   -1d0)
        accuracyDR(15,4,ac) = abs(Dij(15,4) /temp43(3,3,2)   -1d0)
        accuracyDR(16,4,ac) = abs(Dij(16,4) /temp002(1,1)    -1d0)
        accuracyDR(17,4,ac) = abs(Dij(17,4) /temp002(2,2)    -1d0)
        accuracyDR(18,4,ac) = abs(Dij(18,4) /temp002(3,3)    -1d0)
        accuracyDR(19,4,ac) = abs(Dij(19,4) /temp002(2,1)    -1d0)
        accuracyDR(20,4,ac) = abs(Dij(20,4) /temp002(3,1)    -1d0)
        accuracyDR(21,4,ac) = abs(Dij(21,4) /temp002(3,2)    -1d0)
        accuracyDR(22,4,ac) = abs(Dij(22,4) /tempD40000      -1d0)


      DO I1=0,4
           accuracyD(i1,ac)=accuracyDR(1,i1,ac)
            DO I2=2,INDEX(I1)  
       if(accuracyDR(i2,i1,ac).gt.accuracyD(i1,ac)) then
          accuracyD(i1,ac)=accuracyDR(i2,i1,ac)
       endif
          enddo
        enddo

      DO I1=0,4
      if (accuracyD(i1,3).gt.accuracyD(i1,2)) then
         if(printmy) then
         Print*, "Accuracy AC 3 worse than AC2",i1,accuracyD(i1,3),accuracyD(i1,2)
         endif
       return  
      endif
      enddo


cFC        D0=tempD40
cFC        Dij(1,1)=temp1(1)
cFC        Dij(2,1)=temp1(2)
cFC        Dij(3,1)=temp1(3)
cFC        Dij(1,2)=temp2(1,1)
cFC        Dij(2,2)=temp2(2,2)
cFC        Dij(3,2)=temp2(3,3)
cFC        Dij(4,2)=temp2(2,1)
cFC        Dij(5,2)=temp2(3,1)
cFC        Dij(6,2)=temp2(3,2)
cFC        Dij(7,2)=tempD400
cFC        Dij(1,3)=temp3(1,1,1)
cFC        Dij(2,3)=temp3(2,2,2)
cFC        Dij(3,3)=temp3(3,3,3)
cFC        Dij(4,3)=temp3(2,1,1)
cFC        Dij(5,3)=temp3(3,1,1)
cFC        Dij(6,3)=temp3(2,2,1)
cFC        Dij(7,3)=temp3(3,3,1)
cFC        Dij(8,3)=temp3(3,2,2)
cFC        Dij(9,3)=temp3(3,3,2)
cFC        Dij(10,3)=temp3(3,2,1)
cFC        Dij(11,3)=temp001(1)
cFC        Dij(12,3)=temp001(2)
cFC        Dij(13,3)=temp001(3)
cFC        Dij(7,1)=temp002(1,1)
cFC        Dij(8,1)=temp002(2,2)
cFC        Dij(9,1)=temp002(3,3)
cFC        Dij(10,1)=temp002(2,1)
cFC        Dij(11,1)=temp002(3,1)
cFC        Dij(12,1)=temp002(3,2)
cFC        Dij(13,1)=tempD40000
cFC        Dij(1,4)=temp41(1,1,1)
cFC         Dij(2,4)=temp42(2,2,2)
cFC         Dij(3,4)=temp43(3,3,3)
cFC         Dij(4,4)=temp42(1,1,1)
cFC         Dij(5,4)=temp43(1,1,1)
cFC         Dij(6,4)=temp42(2,1,1)
cFC         Dij(7,4)=temp43(2,1,1)
cFC         Dij(8,4)=temp43(3,1,1)
cFC         Dij(9,4)=temp42(2,2,1)
cFC         Dij(10,4)=temp43(2,2,1)
cFC         Dij(11,4)=temp43(3,2,1)
cFC         Dij(12,4)=temp43(3,3,1)
cFC         Dij(13,4)=temp43(2,2,2)
cFC         Dij(14,4)=temp43(3,2,2)
cFC         Dij(15,4)=temp43(3,3,2)
cFC         Dij(16,4)=temp002(1,1)
cFC         Dij(17,4)=temp002(2,2)
cFC         Dij(18,4)=temp002(3,3)
cFC         Dij(19,4)=temp002(2,1)
cFC         Dij(20,4)=temp002(3,1)
cFC         Dij(21,4)=temp002(3,2)
cFC         Dij(22,4)=tempD40000
cFC

 500    D0=tempD40
        Dij(1,1)=temp1(1)
        Dij(2,1)=temp1(2)
        Dij(3,1)=temp1(3)
        Dij(1,2)=temp2(1,1)
        Dij(2,2)=temp2(2,2)
        Dij(3,2)=temp2(3,3)
        Dij(4,2)=temp2(2,1)
        Dij(5,2)=temp2(3,1)
        Dij(6,2)=temp2(3,2)
        Dij(7,2)=tempD400
        Dij(1,3)=temp3(1,1,1)
        Dij(2,3)=temp3(2,2,2)
        Dij(3,3)=temp3(3,3,3)
        Dij(4,3)=temp3(2,1,1)
        Dij(5,3)=temp3(3,1,1)
        Dij(6,3)=temp3(2,2,1)
        Dij(7,3)=temp3(3,3,1)
        Dij(8,3)=temp3(3,2,2)
        Dij(9,3)=temp3(3,3,2)
        Dij(10,3)=temp3(3,2,1)
        Dij(11,3)=temp001(1)
        Dij(12,3)=temp001(2)
        Dij(13,3)=temp001(3)
        Dij(7,1)=temp002(1,1)
        Dij(8,1)=temp002(2,2)
        Dij(9,1)=temp002(3,3)
        Dij(10,1)=temp002(2,1)
        Dij(11,1)=temp002(3,1)
        Dij(12,1)=temp002(3,2)
        Dij(13,1)=tempD40000
        Dij(1,4)=temp41(1,1,1)
         Dij(2,4)=temp42(2,2,2)
         Dij(3,4)=temp43(3,3,3)
         Dij(4,4)=temp42(1,1,1)
         Dij(5,4)=temp43(1,1,1)
         Dij(6,4)=temp42(2,1,1)
         Dij(7,4)=temp43(2,1,1)
         Dij(8,4)=temp43(3,1,1)
         Dij(9,4)=temp42(2,2,1)
         Dij(10,4)=temp43(2,2,1)
         Dij(11,4)=temp43(3,2,1)
         Dij(12,4)=temp43(3,3,1)
         Dij(13,4)=temp43(2,2,2)
         Dij(14,4)=temp43(3,2,2)
         Dij(15,4)=temp43(3,3,2)
         Dij(16,4)=temp002(1,1)
         Dij(17,4)=temp002(2,2)
         Dij(18,4)=temp002(3,3)
         Dij(19,4)=temp002(2,1)
         Dij(20,4)=temp002(3,1)
         Dij(21,4)=temp002(3,2)
         Dij(22,4)=tempD40000
         Dij(1,5)=temp511(1,1,1)
         Dij(2,5)=temp522(2,2,2)
         Dij(3,5)=temp533(3,3,3)
         Dij(4,5)=temp521(1,1,1)
         Dij(5,5)=temp531(1,1,1)
         Dij(6,5)=temp522(1,1,1)
         Dij(7,5)=temp532(1,1,1)
         Dij(8,5)=temp533(1,1,1)
         Dij(9,5)=temp522(2,1,1)
         Dij(10,5)=temp532(2,1,1)
         Dij(11,5)=temp533(2,1,1)
         Dij(12,5)=temp533(3,1,1)
         Dij(13,5)=temp522(2,2,1)
         Dij(14,5)=temp532(2,2,1)
         Dij(15,5)=temp533(2,2,1)
         Dij(16,5)=temp533(3,2,1)
         Dij(17,5)=temp533(3,3,1)
         Dij(18,5)=temp532(2,2,2)
         Dij(19,5)=temp533(2,2,2)
         Dij(20,5)=temp533(3,2,2)
         Dij(21,5)=temp533(3,3,2)
         Dij(22,5)=temp003(1,1,1)
         Dij(23,5)=temp003(2,2,2)
         Dij(24,5)=temp003(3,3,3)
         Dij(25,5)=temp003(2,1,1)
         Dij(26,5)=temp003(3,1,1)
         Dij(27,5)=temp003(2,2,1)
         Dij(28,5)=temp003(3,2,1)
         Dij(29,5)=temp003(3,3,1)
         Dij(30,5)=temp003(3,2,2)
         Dij(31,5)=temp003(3,3,2)
         Dij(32,5)=temp00001(1)
         Dij(33,5)=temp00001(2)
         Dij(34,5)=temp00001(3)
cFC         Dij(1,6)=temp6111(1,1,1)
cFC         Dij(2,6)=temp6222(2,2,2)
cFC         Dij(3,6)=temp6333(3,3,3)
cFC         Dij(4,6)=temp6211(1,1,1)
cFC         Dij(5,6)=temp6311(1,1,1)
cFC         Dij(6,6)=temp6221(1,1,1)
cFC         Dij(7,6)=temp6321(1,1,1)
cFC         Dij(8,6)=temp6331(1,1,1)
cFC         Dij(9,6)=temp6222(1,1,1)
cFC         Dij(10,6)=temp6322(1,1,1)
cFC         Dij(11,6)=temp6332(1,1,1)
cFC         Dij(12,6)=temp6333(1,1,1)
cFC         Dij(13,6)=temp6222(2,1,1)
cFC         Dij(14,6)=temp6322(2,1,1)
cFC         Dij(15,6)=temp6332(2,1,1)
cFC         Dij(16,6)=temp6333(2,1,1)
cFC         Dij(17,6)=temp6333(3,1,1)
cFC         Dij(18,6)=temp6222(2,2,1)
cFC         Dij(19,6)=temp6322(2,2,1)
cFC         Dij(20,6)=temp6332(2,2,1)
cFC         Dij(21,6)=temp6333(2,2,1)
cFC         Dij(22,6)=temp6333(3,2,1)
cFC         Dij(23,6)=temp6333(3,3,1)
cFC         Dij(24,6)=temp6322(2,2,2)
cFC         Dij(25,6)=temp6332(2,2,2)
cFC         Dij(26,6)=temp6333(2,2,2)
cFC         Dij(27,6)=temp6333(3,2,2)
cFC         Dij(28,6)=temp6333(3,3,2)
cFC         Dij(29,6)=temp0041(1,1,1)
cFC         Dij(30,6)=temp0042(2,2,2)
cFC         Dij(31,6)=temp0043(3,3,3)
cFC         Dij(32,6)=temp0042(1,1,1)
cFC         Dij(33,6)=temp0043(1,1,1)
cFC         Dij(34,6)=temp0042(2,1,1)
cFC         Dij(35,6)=temp0043(2,1,1)
cFC         Dij(36,6)=temp0043(3,1,1)
cFC         Dij(37,6)=temp0042(2,2,1)
cFC         Dij(38,6)=temp0043(2,2,1)
cFC         Dij(39,6)=temp0043(3,2,1)
cFC         Dij(40,6)=temp0043(3,3,1)
cFC         Dij(41,6)=temp0043(2,2,2)
cFC         Dij(42,6)=temp0043(3,2,2)
cFC         Dij(43,6)=temp0043(3,3,2)
cFC         Dij(44,6)=temp00002(1,1)
cFC         Dij(45,6)=temp00002(2,2)
cFC         Dij(46,6)=temp00002(3,3)
cFC         Dij(47,6)=temp00002(2,1)
cFC         Dij(48,6)=temp00002(3,1)
cFC         Dij(49,6)=temp00002(3,2)
cFC         Dij(50,6)=tempD4000000
cFC         Dij(1,7) =temp71111(1,1,1)
cFC        Dij(2,7) =temp72222(2,2,2)
cFC        Dij(3,7) =temp73333(3,3,3)
cFC        Dij(4,7) =temp72111(1,1,1)
cFC        Dij(5,7) =temp73111(1,1,1)
cFC        Dij(6,7) =temp72211(1,1,1)
cFC        Dij(7,7) =temp73211(1,1,1)
cFC        Dij(8,7) =temp73311(1,1,1)
cFC        Dij(9,7) =temp72221(1,1,1)
cFC        Dij(10,7)=temp73221(1,1,1)
cFC        Dij(11,7)=temp73321(1,1,1)
cFC        Dij(12,7)=temp73331(1,1,1)
cFC        Dij(13,7)=temp72222(1,1,1)
cFC        Dij(14,7)=temp73222(1,1,1)
cFC        Dij(15,7)=temp73322(1,1,1)
cFC        Dij(16,7)=temp73332(1,1,1)
cFC        Dij(17,7)=temp73333(1,1,1)
cFC        Dij(18,7)=temp72222(2,1,1)
cFC        Dij(19,7)=temp73222(2,1,1)
cFC        Dij(20,7)=temp73322(2,1,1)
cFC        Dij(21,7)=temp73332(2,1,1)
cFC        Dij(22,7)=temp73333(2,1,1)
cFC        Dij(23,7)=temp73333(3,1,1)
cFC        Dij(24,7)=temp72222(2,2,1)
cFC        Dij(25,7)=temp73222(2,2,1)
cFC        Dij(26,7)=temp73322(2,2,1)
cFC        Dij(27,7)=temp73332(2,2,1)
cFC        Dij(28,7)=temp73333(2,2,1)
cFC        Dij(29,7)=temp73333(3,2,1)
cFC        Dij(30,7)=temp73333(3,3,1)
cFC        Dij(31,7)=temp73222(2,2,2)
cFC        Dij(32,7)=temp73322(2,2,2)
cFC        Dij(33,7)=temp73332(2,2,2)
cFC        Dij(34,7)=temp73333(2,2,2)
cFC        Dij(35,7)=temp73333(3,2,2)
cFC        Dij(36,7)=temp73333(3,3,2)
cFC         Dij(37,7)=temp00511(1,1,1)
cFC         Dij(38,7)=temp00522(2,2,2)
cFC         Dij(39,7)=temp00533(3,3,3)
cFC         Dij(40,7)=temp00521(1,1,1)
cFC         Dij(41,7)=temp00531(1,1,1)
cFC         Dij(42,7)=temp00522(1,1,1)
cFC         Dij(43,7)=temp00532(1,1,1)
cFC         Dij(44,7)=temp00533(1,1,1)
cFC         Dij(45,7)=temp00522(2,1,1)
cFC         Dij(46,7)=temp00532(2,1,1)
cFC         Dij(47,7)=temp00533(2,1,1)
cFC         Dij(48,7)=temp00533(3,1,1)
cFC         Dij(49,7)=temp00522(2,2,1)
cFC         Dij(50,7)=temp00532(2,2,1)
cFC         Dij(51,7)=temp00533(2,2,1)
cFC         Dij(52,7)=temp00533(3,2,1)
cFC         Dij(53,7)=temp00533(3,3,1)
cFC         Dij(54,7)=temp00532(2,2,2)
cFC         Dij(55,7)=temp00533(2,2,2)
cFC         Dij(56,7)=temp00533(3,2,2)
cFC         Dij(57,7)=temp00533(3,3,2)
cFC         Dij(58,7)=temp00003(1,1,1)
cFC         Dij(59,7)=temp00003(2,2,2)
cFC         Dij(60,7)=temp00003(3,3,3)
cFC         Dij(61,7)=temp00003(2,1,1)
cFC         Dij(62,7)=temp00003(3,1,1)
cFC         Dij(63,7)=temp00003(2,2,1)
cFC         Dij(64,7)=temp00003(3,2,1)
cFC         Dij(65,7)=temp00003(3,3,1)
cFC         Dij(66,7)=temp00003(3,2,2)
cFC         Dij(67,7)=temp00003(3,3,2)
cFC         Dij(68,7)=temp0000001(1)
cFC         Dij(69,7)=temp0000001(2)
cFC         Dij(70,7)=temp0000001(3)
cFC        Dij(1,8)= temp811111(1,1,1)
cFC        Dij(2,8)= temp822222(2,2,2)
cFC        Dij(3,8)= temp833333(3,3,3)
cFC        Dij(4,8)= temp821111(1,1,1)
cFC        Dij(5,8)= temp831111(1,1,1)
cFC        Dij(6,8)= temp822111(1,1,1)
cFC        Dij(7,8)= temp832111(1,1,1)
cFC        Dij(8,8)= temp833111(1,1,1)
cFC        Dij(9,8)= temp822211(1,1,1)
cFC        Dij(10,8)=temp832211(1,1,1)
cFC        Dij(11,8)=temp833211(1,1,1)
cFC        Dij(12,8)=temp833311(1,1,1)
cFC        Dij(13,8)=temp822221(1,1,1)
cFC        Dij(14,8)=temp832221(1,1,1)
cFC        Dij(15,8)=temp833221(1,1,1)
cFC        Dij(16,8)=temp833321(1,1,1)
cFC        Dij(17,8)=temp833331(1,1,1)
cFC        Dij(18,8)=temp822222(1,1,1)
cFC        Dij(19,8)=temp832222(1,1,1)
cFC        Dij(20,8)=temp833222(1,1,1)
cFC        Dij(21,8)=temp833322(1,1,1)
cFC        Dij(22,8)=temp833332(1,1,1)
cFC        Dij(23,8)=temp833333(1,1,1)
cFC        Dij(24,8)=temp822222(2,1,1)
cFC        Dij(25,8)=temp832222(2,1,1)
cFC        Dij(26,8)=temp833222(2,1,1)
cFC        Dij(27,8)=temp833322(2,1,1)
cFC        Dij(28,8)=temp833332(2,1,1)
cFC        Dij(29,8)=temp833333(2,1,1)
cFC        Dij(30,8)=temp833333(3,1,1)
cFC        Dij(31,8)=temp822222(2,2,1)
cFC        Dij(32,8)=temp832222(2,2,1)
cFC        Dij(33,8)=temp833222(2,2,1)
cFC        Dij(34,8)=temp833322(2,2,1)
cFC        Dij(35,8)=temp833332(2,2,1)
cFC        Dij(36,8)=temp833333(2,2,1)
cFC        Dij(37,8)=temp833333(3,2,1)
cFC        Dij(38,8)=temp833333(3,3,1)
cFC        Dij(39,8)=temp832222(2,2,2)
cFC        Dij(40,8)=temp833222(2,2,2)
cFC        Dij(41,8)=temp833322(2,2,2)
cFC        Dij(42,8)=temp833332(2,2,2)
cFC        Dij(43,8)=temp833333(2,2,2)
cFC        Dij(44,8)=temp833333(3,2,2)
cFC        Dij(47,8)=temp833333(3,3,2)
cFC        Dij(48,8)=temp006111(1,1,1)
cFC        Dij(49,8)=temp006222(2,2,2)
cFC        Dij(50,8)=temp006333(3,3,3)
cFC        Dij(51,8)=temp006211(1,1,1)
cFC        Dij(52,8)=temp006311(1,1,1)
cFC        Dij(53,8)=temp006221(1,1,1)
cFC        Dij(54,8)=temp006321(1,1,1)
cFC        Dij(55,8)=temp006331(1,1,1)
cFC        Dij(56,8)=temp006222(1,1,1)
cFC        Dij(57,8)=temp006322(1,1,1)
cFC        Dij(58,8)=temp006332(1,1,1)
cFC        Dij(59,8)=temp006333(1,1,1)
cFC        Dij(60,8)=temp006222(2,1,1)
cFC        Dij(61,8)=temp006322(2,1,1)
cFC        Dij(62,8)=temp006332(2,1,1)
cFC        Dij(63,8)=temp006333(2,1,1)
cFC        Dij(64,8)=temp006333(3,1,1)
cFC        Dij(65,8)=temp006222(2,2,1)
cFC        Dij(66,8)=temp006322(2,2,1)
cFC        Dij(67,8)=temp006332(2,2,1)
cFC        Dij(68,8)=temp006333(2,2,1)
cFC        Dij(69,8)=temp006333(3,2,1)
cFC        Dij(70,8)=temp006333(3,3,1)
cFC        Dij(71,8)=temp006322(2,2,2)
cFC        Dij(72,8)=temp006332(2,2,2)
cFC        Dij(73,8)=temp006333(2,2,2)
cFC        Dij(74,8)=temp006333(3,2,2)
cFC        Dij(75,8)=temp006333(3,3,2)
cFC        Dij(76,8)=temp000041(1,1,1)
cFC        Dij(77,8)=temp000042(2,2,2)
cFC        Dij(78,8)=temp000043(3,3,3)
cFC        Dij(79,8)=temp000042(1,1,1)
cFC        Dij(80,8)=temp000043(1,1,1)
cFC        Dij(81,8)=temp000042(2,1,1)
cFC        Dij(82,8)=temp000043(2,1,1)
cFC        Dij(83,8)=temp000043(3,1,1)
cFC        Dij(84,8)=temp000042(2,2,1)
cFC        Dij(85,8)=temp000043(2,2,1)
cFC        Dij(86,8)=temp000043(3,2,1)
cFC        Dij(87,8)=temp000043(3,3,1)
cFC        Dij(88,8)=temp000043(2,2,2)
cFC        Dij(89,8)=temp000043(3,2,2)
cFC        Dij(90,8)=temp000043(3,3,2)
cFC        Dij(91,8)=temp0000002(1,1) 
cFC        Dij(92,8)=temp0000002(2,2) 
cFC        Dij(93,8)=temp0000002(3,3) 
cFC        Dij(94,8)=temp0000002(2,1) 
cFC        Dij(95,8)=temp0000002(3,1) 
cFC        Dij(96,8)=temp0000002(3,2) 
cFC        Dij(97,8)=tempD400000000   
cFC        Dij(97,8)=tempD400000000   
cFC        Dij(58,9)=temp0071111(1,1,1)
cFC        Dij(59,9)=temp0072222(2,2,2)
cFC        Dij(60,9)=temp0073333(3,3,3)
cFC        Dij(61,9)=temp0072111(1,1,1)
cFC        Dij(62,9)=temp0073111(1,1,1)
cFC        Dij(63,9)=temp0072211(1,1,1)
cFC        Dij(64,9)=temp0073211(1,1,1)
cFC        Dij(65,9)=temp0073311(1,1,1)
cFC        Dij(66,9)=temp0072221(1,1,1)
cFC        Dij(67,9)=temp0073221(1,1,1)
cFC        Dij(68,9)=temp0073321(1,1,1)
cFC        Dij(69,9)=temp0073331(1,1,1)
cFC        Dij(70,9)=temp0072222(1,1,1)
cFC        Dij(71,9)=temp0073222(1,1,1)
cFC        Dij(72,9)=temp0073322(1,1,1)
cFC        Dij(73,9)=temp0073332(1,1,1)
cFC        Dij(74,9)=temp0073333(1,1,1)
cFC        Dij(75,9)=temp0072222(2,1,1)
cFC        Dij(76,9)=temp0073222(2,1,1)
cFC        Dij(77,9)=temp0073322(2,1,1)
cFC        Dij(78,9)=temp0073332(2,1,1)
cFC        Dij(79,9)=temp0073333(2,1,1)
cFC        Dij(80,9)=temp0073333(3,1,1)
cFC        Dij(81,9)=temp0072222(2,2,1)
cFC        Dij(82,9)=temp0073222(2,2,1)
cFC        Dij(83,9)=temp0073322(2,2,1)
cFC        Dij(84,9)=temp0073332(2,2,1)
cFC        Dij(85,9)=temp0073333(2,2,1)
cFC        Dij(86,9)=temp0073333(3,2,1)
cFC        Dij(87,9)=temp0073333(3,3,1)
cFC        Dij(88,9)=temp0073222(2,2,2)
cFC        Dij(89,9)=temp0073322(2,2,2)
cFC        Dij(90,9)=temp0073332(2,2,2)
cFC        Dij(91,9)=temp0073333(2,2,2)
cFC        Dij(92,9)=temp0073333(3,2,2)
cFC        Dij(93,9)=temp0073333(3,3,2)
cFC        Dij(94,9) =temp0000511(1,1,1)
cFC        Dij(95,9) =temp0000522(2,2,2)
cFC        Dij(96,9) =temp0000533(3,3,3)
cFC        Dij(97,9) =temp0000521(1,1,1)
cFC        Dij(98,9) =temp0000531(1,1,1)
cFC        Dij(99,9) =temp0000522(1,1,1)
cFC        Dij(100,9)=temp0000532(1,1,1)
cFC        Dij(101,9)=temp0000533(1,1,1)
cFC        Dij(102,9)=temp0000522(2,1,1)
cFC        Dij(103,9)=temp0000532(2,1,1)
cFC        Dij(104,9)=temp0000533(2,1,1)
cFC        Dij(105,9)=temp0000533(3,1,1)
cFC        Dij(106,9)=temp0000522(2,2,1)
cFC        Dij(107,9)=temp0000532(2,2,1)
cFC        Dij(108,9)=temp0000533(2,2,1)
cFC        Dij(109,9)=temp0000533(3,2,1)
cFC        Dij(110,9)=temp0000533(3,3,1)
cFC        Dij(111,9)=temp0000532(2,2,2)
cFC        Dij(112,9)=temp0000533(2,2,2)
cFC        Dij(113,9)=temp0000533(3,2,2)
cFC        Dij(114,9)=temp0000533(3,3,2)
cFC        Dij(115,9)= temp0000003(1,1,1)
cFC        Dij(116,9)= temp0000003(2,2,2)
cFC        Dij(117,9)= temp0000003(3,3,3)
cFC        Dij(118,9)= temp0000003(2,1,1)
cFC        Dij(119,9)= temp0000003(3,1,1)
cFC        Dij(120,9)= temp0000003(2,2,1)
cFC        Dij(121,9)= temp0000003(3,2,1)
cFC        Dij(122,9)= temp0000003(3,3,1)
cFC        Dij(123,9)= temp0000003(3,2,2)
cFC        Dij(124,9)= temp0000003(3,3,2)
cFC        Dij(125,9)= temp000000001(1)
cFC        Dij(126,9)= temp000000001(2)
cFC        Dij(127,9)= temp000000001(3)
                    
      return
      End
