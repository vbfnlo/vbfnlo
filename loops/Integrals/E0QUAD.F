
c E0fin_QUAD
c Author: Francisco Campanario
c Date:2008-04-16
c Generalized the old one by Carlo Oleari
C This is valid for general kinematics pisq<=>0

c      complex*32 function E0fin_QUAD(p1s,p2s,p3s,p4s,p5s,s12,s23,s34,s45,s15, 
      subroutine E0fin_QUAD(p1s,p2s,p3s,p4s,p5s,s12,s23,s34,s45,s15, 
     #     D02345,D01345,D01245,D01235,D01234,E0_fin_QUAD)
      implicit none
      complex*32 E0_fin_QUAD
      complex*32 D02345,D01345,D01245,D01235,D01234
      real*16  p1s,p2s,p3s,p4s,p5s,s12,s23,s34,s45,s15
      real*16  p1s2,p2s2,p3s2,p4s2,p5s2,s12s,s23s,s34s,s45s,s15s
      real*16  d,x1,x2,x3,x4,x5

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

      E0_fin_QUAD=-(x1*D02345+x2*D01345+x3*D01245+x4*D01235+x5*D01234)/d
      return
      end
