
c.....C0fin(p1sq,p2sq,p3sq,musq)
c
c             
c              /|===== p3^2
c             / | 
c            /  |
c           /   |
c p1^2 =====    | 
c           \   |
c            \  |
c             \ |
c              \|===== p2^2
c             
c
c.....musq = mu^2 is an external scale
c 
c.....int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2 = 
c          N_ep * C0fin(p1sq,p2sq,p3sq,musq)
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)
c
c     ris(1) = finite part
c     ris(2) = coeff of 1/ep 
c     ris(3) = coeff of 1/ep^2      

      function C0fin(p1sq,p2sq,p3sq,musq)
      implicit none
      complex * 16 C0fin
      real * 8 p1sq,p2sq,p3sq,musq
      complex * 16 ris(3)
      real * 8 qsq(3),tmp(3)
      real * 8 dilog
      complex * 16 lr,l1,l2,lr2,lr3
      complex * 16 I
      parameter (I=(0,1))
      double precision pi
      parameter (pi=3.141592653589793238462643383279502884197D0)
      complex * 16 ipi
      parameter (ipi=(0,3.141592653589793238462643383279502884197D0))
      integer j,offshellleg,imax,ii

      real * 8 r3,r2,detsq,det,x,y,max,lomx,lomy
      real * 8 tiny
      parameter (tiny=1d-7)
      external dilog
      
      if (musq.lt.0) then
         write(*,*) "ERROR in C0fin: mu^2 MUST be a positive number"
         stop
      endif
      
c      write(*,*) "C0 with args",p1sq,p2sq,p3sq,musq

      offshellleg = 0
      do j=1,3 
         qsq(j) = 0.d0
      enddo
      if (abs(p1sq).gt.tiny) then
         offshellleg = offshellleg +1
         qsq(offshellleg) = p1sq
      endif
      if (abs(p2sq).gt.tiny) then
         offshellleg = offshellleg +1
         qsq(offshellleg) = p2sq
      endif
      if (abs(p3sq).gt.tiny) then
         offshellleg = offshellleg +1
         qsq(offshellleg) = p3sq
      endif

      if (offshellleg.eq.1) then
c     TRIANGLE WITH ONLY ONE EXTERNAL INVARIANT      
         l1 =  log(abs(qsq(1)/musq))
         if (qsq(1).lt.0) then
c     do nothing
         else
            l1 = l1 -ipi
         endif
         ris(3) = 1.d0
         ris(2) = -l1
         ris(1) = 1.d0/2*l1**2 - pi**2/6         
         do j=1,3
            ris(j) = ris(j)/qsq(1)
         enddo

      elseif (offshellleg.eq.2) then
c TRIANGLE WITH TWO EXTERNAL INVARIANTS  
         lr = log(abs(qsq(2)/qsq(1)))
         l1 = log(abs(qsq(1)/musq))
         l2 = log(abs(qsq(2)/musq))
         if ((qsq(1).lt.0).and.(qsq(2).lt.0)) then
c     do nothing            
         elseif ((qsq(1).gt.0).and.(qsq(2).lt.0)) then
            lr = lr + ipi
            l1 = l1 - ipi
         elseif ((qsq(1).lt.0).and.(qsq(2).gt.0)) then
            lr = lr - ipi
            l2 = l2 - ipi
         elseif ((qsq(1).gt.0).and.(qsq(2).gt.0)) then
            l1 = l1 - ipi
            l2 = l2 - ipi
         endif
         ris(3) = 0.d0
         ris(2) = lr
         ris(1) = 1.d0/2*(l1**2 - l2**2)
         do j=1,3
            ris(j) = ris(j)/(qsq(1)-qsq(2))
         enddo
      elseif (offshellleg.eq.3) then         
c TRIANGLE WITH THREE EXTERNAL INVARIANTS           
         ris(3) = 0.d0
         ris(2) = 0.d0

c order the qsq(i) with absolute max value in the first position
         max = 0d0
         do ii=1,3
            if (abs(qsq(ii)).gt.max) then
               max = abs(qsq(ii))
               imax = ii
            endif
         enddo
         tmp(1) = qsq(imax)

         j = 2
         do ii=1,3
            if (ii.ne.imax) then               
               tmp(j) = qsq(ii)
               j = j+1
            endif
         enddo
         do ii=1,3
            qsq(ii) = tmp(ii)
         enddo

c         write(*,*) "qsq(i)  ============> ", qsq

         r3 = qsq(3)/qsq(1)
         r2 = qsq(2)/qsq(1)

c         if (sqrt(r2)+sqrt(r3).gt.1d0) then 
c            write(*,*) "NOT YET IMPLEMENTED"
c            stop
c         endif
         
         detsq = (1-r2-r3)**2 - 4*r2*r3
         if (detsq.lt.0d0) then
c            write(*,*) qsq,detsq
c            write(*,*) r2,r3
c            write(*,*) "WARNING: this case has NOT yet been implemented",detsq
c            write(*,*) "RETURN 0 from C0fin function"
            C0fin = 0d0
            detsq = 0
            RETURN
         endif
         
         det = sqrt(detsq)
         x = 1.d0/2/r2*(r3+r2-1+det)
         y = 1.d0/2/r3*(r3+r2-1+det)

c         write(*,*) "x, y ==> ",x,y

         if ((x.gt.1d0).or.(y.gt.1d0)) then             
            write(*,*) "ERROR in C0fin: x and/or y have values 
     #           bigger than one", x, y
            stop
         endif
        
         lomx = log(1-x)
         lomy = log(1-y)
         lr2 = log(abs(r2))
         lr3 = log(abs(r3))
         if (r2.lt.0d0) then
            if (qsq(1).lt.0) then
               lr2 = lr2 - ipi
            else
               lr2 = lr2 + ipi
            endif
         endif
         if (r3.lt.0d0) then
            if (qsq(1).lt.0) then
               lr3 = lr3 - ipi
            else
               lr3 = lr3 + ipi
            endif
         endif


         ris(1) = 1/qsq(1)*(1-x)*(1-y)/(1-x*y)*(2*dilog(x)+2*dilog(y)+
     #        (lomx+lomy)**2+2*lr2*lomy+2*lr3*lomx+lr3*lr2+Pi**2/3)
 

      else
c TRIANGLE WITH ZERO EXTERNAL INVARIANTS
         write(*,*) "Warning: C0fin called with external "//
     #        "invariants equal to zero"
         write(*,*) "dot prods",p1sq,p2sq,p3sq,musq
         ris(3) = 0.d0
         ris(2) = 0.d0
         ris(1) = 0.d0
      endif

      C0fin = ris(1)
      end


c------------- C(m,q1^2,q2^2,P^2)  3-point fuction ---------------------------
c
      complex*16 function C0fin1M(m,q1sq,q2sq,Psq)
      implicit none
      double precision m, q1sq, q2sq, Psq,eps 
      double complex I3point
      External I3point
      parameter(eps=1d-07 )

c Author: Francisco Campanario
c Date: 15 08 2008
c(,,+)
      If (Psq.gt.eps) then
c(-,,+)
      If (q1sq.lt.-eps) then
      C0fin1M=I3point(m,Psq,q2sq,q1sq) 
      Return   
c(,-,+)
      else if(q2sq.lt.-eps) then
      C0fin1M=I3point(m,q1sq,Psq,q2sq)
      Return  
      else 
c(+,+,+)
c      Print*,3
C      C0fin1M=-Dble(I3point(m,-q1sq,-Psq,-q2sq))+(0,1)*DIMAG(I3point(m,-q1sq,-Psq,-q2sq))   
      C0fin1M=I3point(m,q1sq,q2sq,Psq)
      Return   
      endif
c(,,-)
      elseif(Psq.lt.-eps) then
      C0fin1M=I3point(m,q1sq,q2sq,Psq)
c      Print*,4
      Return  
      else
c(-,,0)
      If (q1sq.lt.-eps) then
      C0fin1M=I3point(m,Psq,q2sq,q1sq)   
      Return  
c(,-,0)
      else if(q2sq.lt.-eps) then
      C0fin1M=I3point(m,q1sq,Psq,q2sq)
      Return  
      else 
c(+,+,0)
c      C0fin1M=-Dble(I3point(m,-q1sq,-Psq,-q2sq))+(0,1)*DIMAG(I3point(m,-q1sq,-Psq,-q2sq)) 
      C0fin1M=I3point(m,q1sq,q2sq,Psq)
      Return  
      endif
      endif
      end


c --------------------------------------
c   C0finG2=C0finG full finite
c   only s3 can be zero
c --------------------------------------
      complex*16 function C0finG2(M1,M2,M3,s1,s2,s3,musq)
      implicit none
      real*8 pi
      complex*16 Ipi, Ieps 
      parameter(pi=3.14159265358979324d0,
     & Ipi=(0d0,3.14159265358979324d0),Ieps=(0d0,1d-38))
      double precision M1,M2,M3,s1,s2,s3,musq,m1sq,m2sq,m3sq
      complex*16 x(3), ys(3,2), a, b, c, d, e, f,
     &                 alfa, deno, sqr,ce,ad,abc,de
      complex*16  cdilog, x1, x2, C0
      external    cdilog
      integer i,j
      REAL*8 eps
      parameter (eps=1d-9)
     
      m1sq = M1*M1
      m2sq = M2*M2
      m3sq = M3*M3
c     ---
      a = s3
      b = s2
      c = s1-s2-s3
      d = m3sq-m1sq-s3
      e = m2sq-m3sq+s3-s1
      f = m1sq-Ieps
      ce=m2sq-m3sq-s2
      ad=m3sq-m1sq
      abc=s1
      de=m2sq-m1sq-s1      
c     ---
      alfa = (-c + Sqrt(c*c-4.d0*b*a-2.d0*Ieps*(s1+s2+s3)))/(2.d0*b)
      deno = c + 2.d0*alfa*b

      x(1) = - (d+2.d0*a+(ce)*alfa)/deno
      x(2) = - (d+e*alfa)/((1.d0-alfa)*deno)
      x(3) = (d+e*alfa)/(alfa*deno)

      sqr = Sqrt((ce)*(ce)-4.d0*b*(ad+f))
      ys(1,1) = (-(ce)+sqr)/(2.d0*b)
      ys(1,2) = (-(ce)-sqr)/(2.d0*b)


      sqr = Sqrt((de)*(de)-4.d0*f*(abc))
      ys(2,1) = (-(de)+sqr)/(2.d0*(abc))
      ys(2,2) = (-(de)-sqr)/(2.d0*(abc))


      C0 = (0.d0,0.d0)

      if (abs((a*a)).gt.eps) then 
 

      sqr = Sqrt(d*d-4.d0*a*f)
      ys(3,1) = (-d+sqr)/(2.d0*a)
      ys(3,2) = (-d-sqr)/(2.d0*a)
      do i = 1,3
         do j = 1,2
            x1 = x(i)/(x(i)-ys(i,j))
            x2 = (x(i)-1.d0)/(x(i)-ys(i,j))
            C0 = C0 + ((-1.d0)**i)*(cdilog(x1) -  cdilog(x2))
         enddo
      enddo 
      else
      do i = 1,2
         do j = 1,2
           x1 = x(i)/(x(i)-ys(i,j))
            x2 = (x(i)-1.d0)/(x(i)-ys(i,j))
            C0 = C0 + ((-1.d0)**i)*(cdilog(x1) -  cdilog(x2))

         enddo
      enddo 
      endif
      C0finG2 = C0/deno
      end



