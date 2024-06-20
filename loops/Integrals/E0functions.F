
c E0fin
c Author: Francisco Campanario
c Date:2008-04-16
c Generalized the old one by Carlo Oleari
C This is valid for general kinematics pisq<=>0

      function E0fin(p1s,p2s,p3s,p4s,p5s,s12,s23,s34,s45,s15, 
     &     D02345,D01345,D01245,D01235,D01234)
      implicit none
      complex * 16 E0fin,D02345,D01345,D01245,D01235,D01234
      real * 8 p1s,p2s,p3s,p4s,p5s,s12,s23,s34,s45,s15
      real * 8 p1s2,p2s2,p3s2,p4s2,p5s2,s12s,s23s,s34s,s45s,s15s
      real * 8 d,x1,x2,x3,x4,x5

      p1s2=p1s*p1s
      p2s2=p2s*p2s
      p3s2=p3s*p3s
      p4s2=p4s*p4s
      p5s2=p5s*p5s
      s12s=s12*s12
      s23s=s23*s23
      s34s=s34*s34
      s45s=s45*s45
      s15s=s15*s15

       d=2*(p1s2*p3s*p4s*s34+p2s2*p4s*p5s*s45+p2s*(p3s*p5s*(p5s*s23-s1
     -   5*s45)+s34*s45*(-(p5s*s23)+s15*s45)-p4s*s12*(p5s*s23+s15*s45
     -   ))+s12*(p3s*s15*(-(p5s*s23)+s15*s45)+s23*(p4s*s12*s15+p5s*s2
     -   3*s34-s15*s34*s45))-p1s*(-(p3s2*p5s*s15)+s23*s34*(p4s*s12-s3
     -   4*s45)+p2s*p4s*(p3s*p5s-p4s*s12+s34*s45)+p3s*(p4s*s12*s15+p5
     -   s*s23*s34+s15*s34*s45)))
       x1=-(p3s2*p5s*s15)+p3s*p4s*s12*s15+p3s2*s15s-p3s*s12*s15s+p3s*p
     -   5s*s15*s23+p3s*s12*s15*s23-2*p4s*s12*s15*s23+p3s*p5s*s23*s34
     -   +p4s*s12*s23*s34-2*p3s*s15*s23*s34+s12*s15*s23*s34-p5s*s23s*
     -   s34-s12*s23s*s34+p1s*(-(p3s2*s15)+s23*(p4s-s34)*s34+p2s*p4s*
     -   (p3s-p4s+s34)+p3s*(p4s*(s15-2*s34)+(s15+s23)*s34))+s23s*s34s
     -   +p2s2*p4s*(p4s-p5s-s45)-p3s*s15s*s45+p3s*s15*s34*s45+s15*s23
     -   *s34*s45-s23*s34s*s45+p2s*(-(p4s2*s12)+p3s*(p4s*(p5s-2*s15)+
     -   p5s*s15-2*p5s*s23+s15*s45)+s34*(p5s*s23+(-2*s15+s23)*s45)+p4
     -   s*(p5s*s23+s12*(s15+s23)-2*s23*s34+s15*s45+s34*s45))
       x2=p3s2*p5s2-2*p3s*p4s*p5s*s12+p4s2*s12s-p3s2*p5s*s15+p3s*p4s*s
     -   12*s15+p3s*p5s*s12*s15-p4s*s12s*s15-p3s*p5s2*s23+p3s*p5s*s12
     -   *s23+p4s*p5s*s12*s23-p4s*s12s*s23+p3s*p5s*s23*s34+p4s*s12*s2
     -   3*s34-2*p5s*s12*s23*s34+p3s*p5s*s15*s45-2*p3s*s12*s15*s45+p4
     -   s*s12*s15*s45-2*p3s*p5s*s34*s45-2*p4s*s12*s34*s45+p3s*s15*s3
     -   4*s45+s12*s15*s34*s45+p5s*s23*s34*s45+s12*s23*s34*s45-s23*s3
     -   4s*s45+p2s*(-(p4s2*s12)+s34*(p5s-s45)*s45+p3s*p5s*(p4s-p5s+s
     -   45)+p4s*(p5s*(s12-2*s45)+(s12+s34)*s45))+p1s*(-(p3s2*p5s)-(p
     -   4s-s34)*(p4s*s12-s34*s45)+p3s*(p4s*(p5s+s12-2*s34)+s34*(p5s+
     -   s45)))-s15*s34*s45s+s34s*s45s
       x3=-(p3s*p5s2*s23)+p4s*p5s*s12*s23+p3s*p5s*s15*s23-2*p4s*s12*s1
     -   5*s23+p5s*s12*s15*s23+p5s2*s23s-p5s*s12*s23s+p1s2*p4s*(-p3s+
     -   p4s-s34)-p5s*s23s*s34+p3s*p5s*s15*s45+p4s*s12*s15*s45-p3s*s1
     -   5s*s45-s12*s15s*s45-2*p5s*s15*s23*s45+s12*s15*s23*s45+p5s*s2
     -   3*s34*s45+s15*s23*s34*s45+p1s*(-(p4s2*s12)+p4s*s12*s15-2*p4s
     -   *p5s*s23+p4s*s12*s23+p4s*s23*s34+p5s*s23*s34-2*p4s*s15*s45+p
     -   4s*s34*s45+s15*s34*s45-2*s23*s34*s45+p2s*p4s*(-p4s+p5s+s45)+
     -   p3s*(-2*p5s*s15+p4s*(p5s+s15)+p5s*s23+s15*s45))+p2s*(-((p5s-
     -   s45)*(p5s*s23-s15*s45))+p4s*(p5s*(s23-2*s45)+s15*s45))+s15s*
     -   s45s-s15*s34*s45s
       x4=p1s2*s34*(-p3s-p4s+s34)+p2s2*p5s*(-p4s+p5s-s45)+p2s*(-2*p5s*
     -   s12*s15+p3s*p5s*(-p5s+s15)+p4s*s12*(p5s+s15)-p5s2*s23+p5s*s1
     -   2*s23+p5s*s23*s34+p5s*s15*s45+s12*s15*s45+p5s*s34*s45-2*s15*
     -   s34*s45)+s12*(-(p4s*s12*s15)+p3s*(p5s-s15)*s15+s12*s15s+p5s*
     -   s15*s23-s12*s15*s23-2*p5s*s23*s34+s15*s23*s34-s15s*s45+s15*s
     -   34*s45)+p1s*(p4s*s12*s15+p4s*s12*s34-2*s12*s15*s34+p5s*s23*s
     -   34+s12*s23*s34+p3s*(s15*(s12+s34)+p5s*(-2*s15+s34))-s23*s34s
     -   +s15*s34*s45-s34s*s45+p2s*(p3s*p5s+p4s*(p5s-2*s12+s34)+s34*(
     -   -2*p5s+s45)))
       x5=p1s2*p3s*(p3s-p4s-s34)+p2s2*s45*(-p4s-p5s+s45)+s12*(p3s*(p5s
     -   *s23+s15*(s23-2*s45))+s23*(-(p4s*s12)-s12*s15-p5s*s23+s12*s2
     -   3-s23*s34+s15*s45+s34*s45))+p1s*(-(p3s2*(p5s+s15))+s23*(p4s*
     -   s12+s12*s34-2*s34*s45)+p2s*(-2*p4s*s12+p3s*(p4s+p5s-2*s45)+p
     -   4s*s45+s34*s45)+p3s*(p4s*s12+s12*s15+p5s*s23-2*s12*s23+s23*s
     -   34+s15*s45+s34*s45))+p2s*(p5s*s12*s23+s12*s15*s45+p5s*s23*s4
     -   5-2*s12*s23*s45+s23*s34*s45+p4s*s12*(s23+s45)+p3s*(s15*s45+p
     -   5s*(-2*s23+s45))-s15*s45s-s34*s45s)

       E0fin=-(x1*D02345+x2*D01345+x3*D01245+x4*D01235+x5*D01234)/d
     
      end


*----------------------------------------------------------
      function E01M(m,p1,p2,p3,p4,p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,
     -     D0_2345,D0_1345,D0_1245,D0_1235,D0_1234)
cC     Francisco Campanario
c General E0 function for one massive particle in the looop and p_j<=>0
      implicit none
      complex * 16 E01M
      real * 8 p1,p2,p3,p4,p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,m,msq
      real * 8  tempX5,ItempX5,tempX4(0:4)
      complex * 16 D0_2345,D0_1345,D0_1245,D0_1235,D0_1234

       msq=m*m

       tempX5=-8*(p1**2*p2p4**2*p3+4*p1*p1p2*p2p4**2*p3+4*p1p2**2*p2p4*
     &   *2*p3+2*p1*p1p3*p2p4**2*p3+4*p1p2*p1p3*p2p4**2*p3+p1*p2p4**2
     &   *p3**2+p1p4**2*p2*p3*(p2+2*p2p3+p3)+2*p1*p1p3*p2*p2p4*p3p4+4
     &   *p1p2*p1p3*p2*p2p4*p3p4+4*p1p3**2*p2*p2p4*p3p4-2*p1**2*p2p3*
     &   p2p4*p3p4-8*p1*p1p2*p2p3*p2p4*p3p4-8*p1p2**2*p2p3*p2p4*p3p4-
     &   4*p1*p1p3*p2p3*p2p4*p3p4-8*p1p2*p1p3*p2p3*p2p4*p3p4+2*p1*p1p
     &   2*p2p4*p3*p3p4+4*p1p2**2*p2p4*p3*p3p4+4*p1p2*p1p3*p2p4*p3*p3
     &   p4-2*p1*p2*p2p4*p3*p3p4-4*p1*p2p3*p2p4*p3*p3p4+p1**2*p2*p3p4
     &   **2+2*p1*p1p2*p2*p3p4**2+4*p1*p1p3*p2*p3p4**2+4*p1p2*p1p3*p2
     &   *p3p4**2+4*p1p3**2*p2*p3p4**2+p1*p2**2*p3p4**2-4*p1*p1p2*p2p
     &   3*p3p4**2-8*p1p2**2*p2p3*p3p4**2-8*p1p2*p1p3*p2p3*p3p4**2+4*
     &   p1*p2*p2p3*p3p4**2+4*p1*p2p3**2*p3p4**2+(p1**2*(p2p3**2-p2*p
     &   3)+p1p3**2*p2*(p2+2*p2p4+4*p3p4)+p1*(-(p3*(p2**2+2*p2p3*p2p4
     &   +p2*(2*(p2p3+p2p4)+p3)))+2*p1p2*(2*p2p3**2+(-p2+p2p4)*p3+p2p
     &   3*(p3-p3p4))+2*(p2p3*(p2+2*p2p3)-p2*p3)*p3p4)+p1p2**2*(4*p2p
     &   3**2+4*p2p3*(p3-p3p4)+p3*(4*p2p4+p3+2*p3p4))-2*p1p3*(p1*(p2p
     &   3*p2p4+p2*(p2p3+p3-p3p4))+p1p2*(2*p2p3*p2p4-p2p4*p3+p2*(2*p2
     &   p3+p3-p3p4)+4*p2p3*p3p4)))*p4+(p1p3**2*p2-2*p1p2*p1p3*p2p3+p
     &   1p2**2*p3+p1*(p2p3**2-p2*p3))*p4**2+4*msq*(p1p3**2*p2p4**2-p
     &   1*p2p4**2*p3+p1p4**2*(p2p3**2-p2*p3)-2*p1p2*p1p3*p2p4*p3p4+2
     &   *p1*p2p3*p2p4*p3p4+p1p2**2*p3p4**2-p1*p2*p3p4**2-2*p1p4*(p1p
     &   3*p2p3*p2p4-p1p2*p2p4*p3-p1p3*p2*p3p4+p1p2*p2p3*p3p4)-(p1p3*
     &   *2*p2-2*p1p2*p1p3*p2p3+p1*p2p3**2+p1p2**2*p3-p1*p2*p3)*p4)-2
     &   *p1p4*(p1p3*p2*(p2p4*p3+p2*p3p4+(p2p3+p3)*(2*p3p4+p4))+p1*(p
     &   2*(-(p2p3*p3p4)+p3*(p2p4+p3p4+p4))+p2p3*(p2p4*p3-p2p3*(2*p3p
     &   4+p4)))+p1p2*(p2*(-2*p2p3*p3p4+p3*(2*p2p4+p3p4+p4))-(2*p2p3+
     &   p3)*(-(p2p4*p3)+p2p3*(2*p3p4+p4)))))

      ItempX5=1d0/tempX5


       tempX4(0)=-8*(p1+2*p1p2+p1p3)*p2p4**2*p3-8*p2p4*(p1p3*(p2-2*p2p
     &   3)-2*(p1+2*p1p2)*p2p3+p1p2*p3)*p3p4-8*((p1+p1p2+2*p1p3)*p2-2
     &   *p1p2*p2p3)*p3p4**2+8*(p1*(-p2p3**2+p2*p3)+p1p3*(p2p3*p2p4+p
     &   2*(p2p3+p3-p3p4))+p1p2*(-2*p2p3**2+p2*p3-p2p3*p3-p2p4*p3+p2p
     &   3*p3p4))*p4+8*p1p4*(p2*(-(p2p3*p3p4)+p3*(p2p4+p3p4+p4))+p2p3
     &   *(p2p4*p3-p2p3*(2*p3p4+p4)))
       tempX4(1)=-8*(p1p4**2*(p2+p2p3)*p3-p1*p2p4**2*p3-2*p1p2*p2p4**2
     &   *p3-p1p3*p2p4**2*p3+p1*p1p3*p2p4*p3p4+2*p1p2*p1p3*p2p4*p3p4+
     &   2*p1p3**2*p2p4*p3p4-p1p3*p2*p2p4*p3p4+2*p1*p2p3*p2p4*p3p4+4*
     &   p1p2*p2p3*p2p4*p3p4+2*p1p3*p2p3*p2p4*p3p4-p1*p2p4*p3*p3p4-p1
     &   p2*p2p4*p3*p3p4-p1*p1p2*p3p4**2-2*p1p2**2*p3p4**2-2*p1p2*p1p
     &   3*p3p4**2-p1p2*p2*p3p4**2-2*p1p3*p2*p3p4**2+2*p1*p2p3*p3p4**
     &   2+2*p1p2*p2p3*p3p4**2+(p1p3**2*(p2+p2p4)-p1*(p2p3**2+(-p1p2+
     &   p2p4)*p3+p2p3*(p3-p3p4))+p1p3*(-(p1*p2p3)-2*p1p2*p2p3+p2*p2p
     &   3+p2p3*p2p4+p1p2*p3+p2*p3-(p1p2+p2)*p3p4)+p1p2*(-2*p2p3**2+(
     &   2*p1p2+p2-p2p4)*p3+p2p3*(-p3+p3p4)))*p4+p1p4*(-(p1*p2p4*p3)-
     &   2*p1p2*p2p4*p3-p1p3*p2p4*p3+p2*p2p4*p3+p2p3*p2p4*p3-2*p1p3*p
     &   2*p3p4+p1*p2p3*p3p4+2*p1p2*p2p3*p3p4-2*p1p3*p2p3*p3p4-p2*p2p
     &   3*p3p4-2*p2p3**2*p3p4+p1p2*p3*p3p4+p2*p3*p3p4+(-(p2p3*(p1p3+
     &   p2p3))+(p1p2+p2)*p3)*p4))
       tempX4(2)=8*(p1p4**2*p2p3*(p2+2*p2p3+p3)+p1*(p2p4+p3p4)*(p1p3*p
     &   2p4-p2p4*p3+(-p1p2+p2+2*p2p3)*p3p4)-p1*(p2*p2p3+2*p2p3**2+p1
     &   p3*(p2+p2p3)+p2p3*p2p4+p2p3*p3+p2p4*p3-p1p2*(p2p3+p3)-(p2+p2
     &   p3)*p3p4)*p4+(p1p2+p1p3)*(2*p1p3*p2p4*(p2p4+p3p4)-2*p1p2*p3p
     &   4*(p2p4+p3p4)+p1p3*(-p2+p2p4)*p4+p1p2*(2*p2p3+p3-p3p4)*p4)-p
     &   1p4*(p1*p2p4*(p2p3+p3)-p1*(p2+p2p3)*p3p4+p1p2*(2*p2p3*p2p4-(
     &   p2+p3)*p3p4-(p2p3+p3)*p4)+p1p3*(p2p4*p3+p2*(p2p4+p4)+p2p3*(4
     &   *p2p4+2*p3p4+p4))))
       tempX4(3)=-8*(-(p1*p1p2*p2p4*p3)-2*p1p2**2*p2p4*p3+p1*p2*p2p4*p
     &   3+p1*p2p3*p2p4*p3+p1*p2p4**2*p3+p1p4**2*p2*(p2p3+p3)+p1*p1p2
     &   *p2p3*p3p4+2*p1p2**2*p2p3*p3p4-p1*p2*p2p3*p3p4-2*p1*p2p3**2*
     &   p3p4-p1*p1p2*p2p4*p3p4-2*p1p2**2*p2p4*p3p4+p1*p2*p2p4*p3p4-2
     &   *p1*p2p3*p2p4*p3p4-p1p2**2*p3*p3p4+p1*p2*p3*p3p4-2*p1p2**2*p
     &   3p4**2+2*p1*p2*p3p4**2+(-(p1*p2p3*(-p1p2+p2+p2p3+p2p4))+p1p2
     &   **2*(2*p2p3-p3p4)+p1*p2*p3p4)*p4-p1p3**2*p2*(p2p4+2*p3p4+p4)
     &   +p1p4*(p1*(-(p2p3*(p2p3+p2p4))+p2*(p3+p3p4))+p1p2*(-2*p2p3**
     &   2-2*p2p4*p3+p2*(p3+p3p4)+p2p3*(-2*p2p4-p3+2*p3p4+p4)))+p1p3*
     &   (p1p4*p2*(p2p3-p2p4+p3-2*p3p4-p4)+p1*(p2p4*(p2p3+p2p4)-p2*(p
     &   3p4+p4))+p1p2*(-(p2*(p3p4+p4))+2*p2p3*(p2p4+2*p3p4+p4)+p2p4*
     &   (2*p2p4-p3+2*p3p4+p4))))
       tempX4(4)=8*(-(p1p3**2*p2*(p2p4+2*p3p4+p4))+p1*(-(p1p4*p2p3**2)
     &   +p1p4*p2*p3-p1p2*p2p4*p3+p2*p2p4*p3+p2p3*p2p4*p3+p1p2*p2p3*p
     &   3p4-p2*p2p3*p3p4-2*p2p3**2*p3p4+p2*p3*p3p4-p2p3**2*p4+p2*p3*
     &   p4)+p1p2*(p1p4*(p2*p3-p2p3*(2*p2p3+p3))-p1p2*(-2*p2p3*p3p4+p
     &   3*(2*p2p4+p3p4+p4)))+p1p3*(p1p4*p2*(p2p3+p3)+p1*(p2p3*p2p4-p
     &   2*p3p4)+p1p2*(-(p2p4*p3)-p2*p3p4+2*p2p3*(p2p4+2*p3p4+p4))))


       E01M=(D0_2345*tempX4(0)+D0_1345*tempX4(1)+D0_1245*tempX4(2)
     & +D0_1235*tempX4(3)+D0_1234*tempX4(4))*ItempX5    

      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c


      function E02(psq, qsq, lsq, tsq,
     #     pq, pl, pt, ql, qt, lt, 
     #     D0_2345,D0_1345,D0_1245,D0_1235,D0_1234)
      implicit none
      complex * 16 E02
      real * 8 psq, qsq, lsq, tsq, pq, pl, pt, ql, qt, lt
      real * 8 q1s, q2s, ppp, pq1, ppq1, pq2, ppq2, q1q2, pps
      real * 8 det, coeff(1:5)
      complex * 16 D0_2345,D0_1345,D0_1245,D0_1235,D0_1234, J(1:5)
      integer i
      real * 8 tiny
      logical lpr
      common /e0print/lpr

      tiny = 1d-4

      if (abs(psq).gt.tiny) then
         write(*,*) "p^2 =",psq
         write(*,*) "E0 computed only with p^2 and pp^2 = 0"
         stop
      endif

      q1s  = qsq
      q2s  = lsq
      ppp  = psq+pq+pl+pt
      pq1  = pq
      ppq1 = pq + qsq+ql+qt
      pq2  = pl
      ppq2 = pl + ql + lsq + lt
      q1q2 = ql
      pps = lsq+2*lt+2*ql+2*pl+tsq+2*qt+2*pt+qsq+2*pq+psq
      if (abs(pps).gt.tiny) then
         write(*,*) "pp^2 =",pps
         write(*,*) "E0 computed only with p^2 and pp^2 = 0"
         lpr = .true.
c         stop
      endif


      det = -8*ppp*(4*q1s*pq2**2*ppq1+q1s**2*q2s*ppp-4*pq1*pq2*q1s*ppq2-
     #4*pq1*pq2*ppq1*q2s+4*q1s*pq2*q1q2*ppp+4*q1s*pq2*q2s*ppp-4*pq1*q2s*
     #q1q2*ppp+4*q1s*q2s*pq1*ppp-8*pq1*pq2*ppq1*q1q2-8*pq1*q1q2*pq2*ppp+
     #4*q1s*pq1*q1q2*ppq2+2*q1s*pq1*q2s*ppq2-8*pq1*q2s*ppq1*q1q2+2*q2s*p
     #q2*q1s*ppq1+4*q1s*pq2*ppq1*q1q2+2*q1s*q1q2*q2s*ppp-2*q1s**2*pq2*pp
     #q2+4*pq1**2*q2s*ppq2+8*pq1**2*q1q2*ppq2+4*q1s*pq2**2*ppp-8*q1q2**2
     #*pq1*ppp+4*pq1**2*q2s*ppp+q1s*q2s**2*ppp-8*q1q2**2*pq1*ppq1-2*pq1*
     #q2s**2*ppq1)

      
c      write(*,*)  "DETERMINANTE in E0",det

      coeff(1)=(-16*pq1**2*q2s*ppp-8*q1s*pq1*ppq2**2+16*q1q2**2*ppq1*ppp
     #-8*pq1**2*q2s*ppq2-16*pq1**2*q1q2*ppq2-16*q1q2*pq2*ppq1**2-8*q1s*p
     #q2**2*ppq1-16*q1s*pq2**2*ppp+16*q1q2**2*pq1*ppp+16*q1q2**2*ppp**2-
     #8*q2s*pq2*ppq1**2-16*q1s*q2s*ppp**2+8*pq1*pq2*q1s*ppq2+8*pq1*pq2*p
     #pq1*q2s-8*q1s*pq2*q1q2*ppp-8*q1s*pq2*q2s*ppp+8*pq1*q2s*q1q2*ppp-8*
     #q1s*q2s*pq1*ppp+8*q1s*pq2*ppq1*ppq2-8*q1s*q1q2*ppq2*ppp-8*q1s*q2s*
     #ppq1*ppp-8*q1s*q2s*ppq2*ppp+8*pq1*q2s*ppq1*ppq2+8*q1q2*q2s*ppq1*pp
     #p+16*pq1*pq2*ppq1*q1q2+16*q1s*pq2*ppq2*ppp+16*pq1*q2s*ppq1*ppp+32*
     #pq1*q1q2*pq2*ppp-16*pq1*q1q2*ppq2*ppp-16*q1q2*pq2*ppq1*ppp+16*pq1*
     #q1q2*ppq1*ppq2)/det


      coeff(2)=(16*pq1*q2s*ppp**2+8*q2s*pq2*ppq1*ppp+16*pq2**2*ppq1**2-1
     #6*pq1*pq2*ppq2*ppp-16*q1q2*pq2*ppp**2+16*pq1**2*ppq2**2-32*pq1*pq2
     #*ppq1*ppq2+16*pq1*q2s*ppq1*ppp-8*q1q2*q2s*ppp**2+8*q1s*q2s*ppp**2+
     #16*pq2**2*ppq1*ppp+8*pq1*q2s*ppq2*ppp-16*q1q2**2*ppp**2+8*q1s*pq1*
     #ppq2**2-16*q1q2**2*ppq1*ppp+8*q2s*pq2*ppq1**2-8*q1s*pq2*ppq1*ppq2-
     #16*pq1*q1q2*ppq1*ppq2+8*q1s*q1q2*ppq2*ppp+16*q1q2*pq2*ppq1**2+8*q1
     #s*q2s*ppq1*ppp+8*q1s*q2s*ppq2*ppp-8*pq1*q2s*ppq1*ppq2-8*q1q2*q2s*p
     #pq1*ppp)/det


      coeff(3)=8*ppp*(-2*pq1*q1q2*ppp+2*q1s*pq2*ppq2-2*pq1*q2s*ppq1-4*pq
     #1*q1q2*ppq1-2*pq1*q1q2*ppq2-2*pq1*pq2*ppq1+q1s*q1q2*ppp+2*pq1**2*p
     #pq2+2*q1s*pq2*ppp+q1s*pq1*ppq2+q1s*pq2*ppq1-pq1*q2s*ppq2-q2s*pq2*p
     #pq1+2*pq1*pq2*ppq2+2*q1q2**2*ppp+q1q2*q2s*ppp-2*pq2**2*ppq1-2*pq1*
     #q2s*ppp-2*q1q2*pq2*ppq1+2*q1q2*pq2*ppp)/det


      coeff(4)=-8*ppp*(-2*pq1*q1q2*ppp-2*pq1*q2s*ppq1-4*pq1*q1q2*ppq1-2*
     #pq1*pq2*ppq1+q1s*q2s*ppp+q1s*q1q2*ppp+2*pq1**2*ppq2+2*q1s*pq2*ppp+
     #q1s*pq1*ppq2+q1s*pq2*ppq1)/det


      coeff(5)=(16*pq1**2*q1q2*ppq2+8*pq1**2*q2s*ppq2+16*q1s*pq2**2*ppp-
     #8*pq1*pq2*q1s*ppq2-16*pq1*pq2*ppq1*q1q2-8*pq1*pq2*ppq1*q2s+8*q1s*p
     #q2*q1q2*ppp+8*q1s*pq2*q2s*ppp-8*pq1*q2s*q1q2*ppp-32*pq1*q1q2*pq2*p
     #pp-16*q1q2**2*pq1*ppp+8*q1s*pq2**2*ppq1+16*pq1**2*q2s*ppp+8*q1s*q2
     #s*pq1*ppp)/det


c D0(q1,q2,pp-p-q1-q2) 
      J(1) = D0_2345
c D0(p+q1,q2,pp-p-q1-q2)
      J(2) = D0_1345
c D0(p,q1+q2,pp-p-q1-q2)
      J(3) = D0_1245
c D0(p,q1,pp-p-q1)
      J(4) = D0_1235
c D0(p,q1,q2)
      J(5) = D0_1234
      
      E02 = 0d0
      do i=1,5
         E02 = E02 + J(i)*coeff(i)
      enddo
      end

