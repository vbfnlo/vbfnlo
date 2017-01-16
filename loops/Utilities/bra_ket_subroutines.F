
c p = physical momentum  p[0] = E > 0
c ferm_type = fermion (+1) or antifermion (-1)
c gives ket(1..2,-1:1) for massless spinors
c  -1 : 1 helicity
      subroutine ket(p,ferm_type,psi)
      implicit none
      real * 8 p(0:3)
      integer ferm_type
      complex * 16 psi(1:2,-1:1)
      real * 8 E,px,py,pz
      complex * 16 I
      parameter (I=(0,1))
      real * 8 tiny
      parameter (tiny=1d-13)
      
      E  = p(0)
      px = p(1)
      py = p(2)
      pz = p(3)      
      if (ferm_type.eq.1) then         
c     fermion => u(p)
         if (abs(E+pz).gt.tiny) then
            psi(1,-1) = (-px+I*py)/sqrt(E+pz)
            psi(2,-1) = (E+pz)/sqrt(E+pz)
            psi(1,1)  = (E+pz)/sqrt(E+pz)
            psi(2,1)  = (px+I*py)/sqrt(E+pz)
         else
            psi(1,-1) = -sqrt(2*E)
            psi(2,-1) = 0
            psi(1,1)  = 0
            psi(2,1)  = sqrt(2*E)
         endif
      else
c     antifermion => v(p)
         if (abs(E+pz).gt.tiny) then
            psi(1,1) = -(E+pz)/sqrt(E+pz)
            psi(2,1) = (-px-I*py)/sqrt(E+pz)
            psi(1,-1)  = (px-I*py)/sqrt(E+pz)
            psi(2,-1)  = -(E+pz)/sqrt(E+pz)
         else
            psi(1,1) = 0
            psi(2,1) = -sqrt(2*E)
            psi(1,-1)  = sqrt(2*E)
            psi(2,-1)  = 0
         endif
      endif
      end


c p = physical momentum  p[0] = E > 0
c ferm_type = fermion (+1) or antifermion (-1)
c gives bra(1..2,-1:1) for massless spinors
c  -1 : 1 helicity
      subroutine bra(p,ferm_type,psi)
      implicit none
      real * 8 p(0:3)
      integer ferm_type
      complex * 16 psi(1:2,-1:1)
      real * 8 E,px,py,pz
      complex * 16 I
      parameter (I=(0,1))
      
      real * 8 tiny
      parameter (tiny=1d-13)
      
      E  = p(0)
      px = p(1)
      py = p(2)
      pz = p(3)
      if (ferm_type.eq.1) then
c fermion => u(p)
         if (abs(E+pz).gt.tiny) then
            psi(1,-1) = (-px-I*py)/sqrt(E+pz)
            psi(2,-1) = (E+pz)/sqrt(E+pz)
            psi(1,1)  = (E+pz)/sqrt(E+pz)
            psi(2,1)  = (px-I*py)/sqrt(E+pz)
         else
            psi(1,-1) = -sqrt(2*E)
            psi(2,-1) = 0
            psi(1,1)  = 0
            psi(2,1)  = sqrt(2*E)
         endif
      else
c     antifermion => v(p)
         if (abs(E+pz).gt.tiny) then
            psi(1,1) = -(E+pz)/sqrt(E+pz)
            psi(2,1) = (-px+I*py)/sqrt(E+pz)
            psi(1,-1)  = (px+I*py)/sqrt(E+pz)
            psi(2,-1)  = -(E+pz)/sqrt(E+pz)
         else
            psi(1,1) = 0
            psi(2,1) = -sqrt(2*E)
            psi(1,-1)  = sqrt(2*E)
            psi(2,-1)  = 0
         endif
      endif
      end


c                           ( ris(1))
c gives p_slash |ket,hel> = (       )
c                           ( ris(2))
      subroutine slash_ket(p,ket,hel,ris)
      implicit none
      complex * 16 p(0:3)
      complex * 16 ket(1:2)
      integer hel
      complex * 16 ris(1:2)
      complex * 16 pslash(1:2,1:2)
      complex * 16 I
      parameter (I=(0,1))

      pslash(1,1) = p(0) - hel*p(3)
      pslash(1,2) = -hel*(p(1) - I*p(2))
      pslash(2,1) = -hel*(p(1) + I*p(2))
      pslash(2,2) = p(0) + hel*p(3)
      ris(1) = pslash(1,1)*ket(1)+pslash(1,2)*ket(2)
      ris(2) = pslash(2,1)*ket(1)+pslash(2,2)*ket(2)
      end


c gives <bra,hel | p_slash = ( ris(1), ris(2) )      

      subroutine slash_bra(p,bra,hel,ris)
      implicit none
      complex * 16 p(0:3)
      complex * 16 bra(1:2)
      integer hel
      complex * 16 ris(1:2)
      complex * 16 pslash(1:2,1:2)
      complex * 16 I
      parameter (I=(0,1))

      pslash(1,1) = p(0) - hel*p(3)
      pslash(1,2) = -hel*(p(1) - I*p(2))
      pslash(2,1) = -hel*(p(1) + I*p(2))
      pslash(2,2) = p(0) + hel*p(3)
      ris(1) = pslash(1,1)*bra(1)+pslash(2,1)*bra(2)
      ris(2) = pslash(1,2)*bra(1)+pslash(2,2)*bra(2)
      end


c <bra, helbra| slash_p_1 slash_p_2 ... slash_p_n |ket, helket>
c where list_mom = [p_1, p_2,...., p_n]
c helket = -(-1)^n*helbra, otherwise the result is zero!!!
      subroutine bra_slash_ket(bra,ket,helket,list_mom,n,ris)
      implicit none
      complex * 16 bra(2,-1:1),ket(2,-1:1)
      integer helket
      integer n
      complex * 16 list_mom(0:3,n)
c      dimension list_mom(0:3,n)
      complex * 16 ris

      complex * 16 b(2),k(2),new_ket(2)
      complex * 16 I
      parameter (I=(0,1))
      integer helbra,ii,hel
      
      helbra = -(-1)**n*helket
      hel = helket
      b(1) = bra(1,helbra)
      b(2) = bra(2,helbra)
      k(1) = ket(1,helket)
      k(2) = ket(2,helket)
      
      do ii=n,1,-1
         call slash_ket(list_mom(0,ii),k,hel,new_ket)
         k(1) = new_ket(1)
         k(2) = new_ket(2)
         hel = -hel
      enddo
      ris = b(1)*k(1)+b(2)*k(2)     
      end


C this subroutine adds the emission of a gluon/photon to a fermionic line 
c that has only ONE intermediate propagator!!!  This means that, BEFORE 
c emitting the gluon/phot, the expression of the line is
c <bra, helbra| slash_curr_vertex_1, slash_p1, slash_curr_vertex_2|ket, helket>
c where list_mom = [slash_curr_vertex_1, slash_p1, slash_curr_vertex_2]
c helket = -(-1)^n*helbra, otherwise the result is zero!!!
c NB!!! the momentum of the intermediate propagator is the momentum with 
c NO emission. On the other hand, the external momenta are such that 
c momentum conservation is satisfied!!
c You have to provide the momentum of the emitted gluon/photon (pg) 
c and its pol vector
c This procedure attach the right propagators too, when the gluon/photon 
c is emitted 
      subroutine bra_slash_ket_gluon
     -     (bra,ket,helket,pbra,pket,pg,eps,list_mom,ris)
      implicit none
      complex * 16 bra(2,-1:1),ket(2,-1:1)
      integer helket
      real * 8 pbra(0:3),pket(0:3),pg(0:3),eps(0:3)
      complex * 16 list_mom(0:3,3)
c      dimension list_mom(0:3,n)
      complex * 16 ris

      complex * 16 ris1,ris2,ris3
c      complex * 16 I
c      parameter (I=(0,1))
      complex*16 local_list(0:3,5)
      integer mu
      complex * 16 dotcc
c emissione dalla gamba esterna bra
      do mu=0,3
         local_list(mu,1) = eps(mu)
         local_list(mu,2) = pg(mu)+pbra(mu)
         local_list(mu,3) = list_mom(mu,1)
         local_list(mu,4) = list_mom(mu,2)
         local_list(mu,5) = list_mom(mu,3)
      enddo
      call bra_slash_ket(bra,ket,helket,local_list,5,ris1)

      ris1 = ris1 / dotcc(local_list(0,2),local_list(0,2))/
     #     dotcc(local_list(0,4),local_list(0,4))


c emissione dal propagatore interno
      do mu=0,3
         local_list(mu,1) = list_mom(mu,1)
         local_list(mu,2) = list_mom(mu,2)-pg(mu)
         local_list(mu,3) = eps(mu)
         local_list(mu,4) = list_mom(mu,2)
         local_list(mu,5) = list_mom(mu,3)
      enddo
      call bra_slash_ket(bra,ket,helket,local_list,5,ris2)
      ris2 = ris2 / dotcc(local_list(0,2),local_list(0,2)) /
     #     dotcc(local_list(0,4),local_list(0,4))


c emissione dalla gamba esterna ket
      do mu=0,3
         local_list(mu,1) = list_mom(mu,1)
         local_list(mu,2) = list_mom(mu,2)-pg(mu)
         local_list(mu,3) = list_mom(mu,3)
         local_list(mu,4) = pket(mu)-pg(mu)
         local_list(mu,5) = eps(mu)
      enddo
      call bra_slash_ket(bra,ket,helket,local_list,5,ris3)
      ris3 = ris3 / dotcc(local_list(0,2),local_list(0,2))/
     #     dotcc(local_list(0,4),local_list(0,4))
      
      ris = ris1+ris2+ris3
      end


c <bra, helbra| slash_p_1 slash_p_2 ... slash_p_n |ket, helket>
c where list_mom = [p_1, p_2,...., p_n]
c helket = -(-1)^n*helbra, otherwise the result is zero!!!
      subroutine bra_slash_ket1(bra,ket,helbra,list_mom,n,ris)
      implicit none
      complex * 16 bra(2,-1:1),ket(2,-1:1)
      integer helket
      integer n
      complex * 16 list_mom(0:3,n)
c      dimension list_mom(0:3,n)
      complex * 16 ris

      complex * 16 b(2),k(2),new_bra(2)
      complex * 16 I
      parameter (I=(0,1))
      integer helbra,ii,hel
      
      helket = -(-1)**n*helbra
      hel = helbra
      b(1) = bra(1,helbra)
      b(2) = bra(2,helbra)
      k(1) = ket(1,helket)
      k(2) = ket(2,helket)
      
      do ii=1,n
         call slash_bra(list_mom(0,ii),b,hel,new_bra)
         b(1) = new_bra(1)
         b(2) = new_bra(2)
         hel = -hel
      enddo
      ris = b(1)*k(1)+b(2)*k(2)     
      end


c the helicity of the bra and the ket are the same, otherwise this is zero!
c sig^mu = <bra, helbra| gamma^mu |ket, helket>
      subroutine bra_gamma_ket(bra,ket,hel,sig)
      implicit none
      complex * 16 bra(2,-1:1),ket(2,-1:1)
      complex * 16 b(2),k(2)
      integer hel
      complex * 16 sig(0:3)
      complex * 16 I
      parameter (I=(0,1))

      b(1) = bra(1,hel)
      b(2) = bra(2,hel)
      k(1) = ket(1,hel)
      k(2) = ket(2,hel)

      sig(0) = b(1)*k(1)+k(2)*b(2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c PLEASE NOTE the sign of the next components!!!!
c in fact  slash_p = gamma_0*p_0 - \vec{gamma} dot \vec(p)
c here we have just all the indexes up or down ==> the minus sign on the
c previous equation IS NOT there!!!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      sig(1) = hel*(b(1)*k(2)+k(1)*b(2))
      sig(2) = hel*(-I*b(1)*k(2)+I*k(1)*b(2))
      sig(3) = hel*(b(1)*k(1)-k(2)*b(2))
      end


c this procedure attach a gluon/photon with polarization vector eps along 
c all the fermionic line  (NO COLOR STRUCTURE!!)
c the helicity of the bra and the ket are the same, otherwise this is zero!
c sig^mu = 1/dot0p(pbra+pglu,pbra+pglu) * 
c   <bra, pbra, hel| eps_slash, (pbra+pglu)_slash, gamma^mu |ket, pket, hel>+
c          1/dot0p(pket-pglu,pket-pglu) * 
c   <bra, pbra, hel| gamma^mu, (pket-pglu)_slash, eps_slash |ket, pket, hel>

      subroutine bra_gamma_ket_gluon(bra,ket,hel,pbra,pket,pglu,eps,sig)
      implicit none
      complex * 16 bra(2,-1:1),ket(2,-1:1)
      integer hel
      real * 8 pbra(0:3),pket(0:3),pglu(0:3),eps(0:3)
      complex * 16 sig(0:3)
      complex * 16 I
      parameter (I=(0,1))
      complex*16 list_mom(0:3,3)
      real * 8 gamma(0:3),prop1,prop2,dotrr
      integer mu,nu
      complex * 16 J1(0:3),J2(0:3)


      do mu=0,3
         list_mom(mu,1) = eps(mu)
         list_mom(mu,2) = pglu(mu)+pbra(mu)
      enddo         

      do mu=0,3
         do nu=0,3
            gamma(nu)=0.d0
         enddo
         if (mu.eq.0) then
            gamma(mu) = 1.d0
         else
            gamma(mu) = -1.d0
         endif
         do nu=0,3
            list_mom(nu,3) = gamma(nu)
         enddo 
         call bra_slash_ket(bra,ket,hel,list_mom,3,J1(mu))
      enddo


      do mu=0,3
         list_mom(mu,2) = pket(mu)-pglu(mu)
         list_mom(mu,3) = eps(mu)
      enddo         

      do mu=0,3
         do nu=0,3
            gamma(nu)=0.d0
         enddo
         if (mu.eq.0) then
            gamma(mu) = 1.d0
         else
            gamma(mu) = -1.d0
         endif
         do nu=0,3
            list_mom(nu,1) = gamma(nu)
         enddo 
         call bra_slash_ket(bra,ket,hel,list_mom,3,J2(mu))
      enddo


      prop1 = 2*dotrr(pglu,pbra)
      prop2 = -2*dotrr(pglu,pket)

      do mu=0,3
         sig(mu) = J1(mu)/prop1 + J2(mu)/prop2
      enddo
      end


*-- Author :    F. James, modified by Mike Seymour
C-----------------------------------------------------------------------
      FUNCTION random2(iseed1,iseed2)
C     MAIN RANDOM NUMBER GENERATOR
C     USES METHOD OF l"Ecuyer, (VIA F.JAMES, COMP PHYS COMM 60(1990)329)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION random2
      INTEGER ISEED1,iseed2,K,IZ
      K=ISEED1/53668
      ISEED1=40014*(ISEED1-K*53668)-K*12211
      IF (ISEED1.LT.0) ISEED1=ISEED1+2147483563
      K=ISEED2/52774
      ISEED2=40692*(ISEED2-K*52774)-K*3791
      IF (ISEED2.LT.0) ISEED2=ISEED2+2147483399
      IZ=ISEED1-ISEED2
      IF (IZ.LT.1) IZ=IZ+2147483562
      RANDOM2=DBLE(IZ)/2147483589
c      RANDOM2=DBLE(IZ)*4.656613001013252D-10
C--->                (4.656613001013252D-10 = 1.D0/2147483589)
      end


      subroutine from_p6_to_p(p1,p2,p3,p4,p5,p6,p)
      implicit none
      real * 8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
      real * 8 p(0:3,6)
      integer mu
      do mu=0,3
         p(mu,1) = p1(mu)
         p(mu,2) = p2(mu)
         p(mu,3) = p3(mu)
         p(mu,4) = p4(mu)
         p(mu,5) = p5(mu)
         p(mu,6) = p6(mu)
      enddo
      end

      subroutine from_p7_to_p(p1,p2,p3,p4,p5,p6,p7,p)
      implicit none
      real * 8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3),p7(0:3)
      real * 8 p(0:3,7)
      integer mu
      do mu=0,3
         p(mu,1) = p1(mu)
         p(mu,2) = p2(mu)
         p(mu,3) = p3(mu)
         p(mu,4) = p4(mu)
         p(mu,5) = p5(mu)
         p(mu,6) = p6(mu)
         p(mu,7) = p7(mu)
      enddo
      end

      subroutine from_p_to_p6(p,p1,p2,p3,p4,p5,p6)
      implicit none
      real * 8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
      real * 8 p(0:3,6)
      integer mu
      do mu=0,3
         p1(mu) = p(mu,1)
         p2(mu) = p(mu,2)
         p3(mu) = p(mu,3)
         p4(mu) = p(mu,4)
         p5(mu) = p(mu,5)
         p6(mu) = p(mu,6)
      enddo
      end

      subroutine from_p_to_p7(p,p1,p2,p3,p4,p5,p6,p7)
      implicit none
      real * 8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3),p7(0:3)
      real * 8 p(0:3,7)
      integer mu
      do mu=0,3
         p1(mu) = p(mu,1)
         p2(mu) = p(mu,2)
         p3(mu) = p(mu,3)
         p4(mu) = p(mu,4)
         p5(mu) = p(mu,5)
         p6(mu) = p(mu,6)
         p7(mu) = p(mu,7)
      enddo
      end


      subroutine exchange_mom(p,i,j,dim,pnew)
      implicit none
      real * 8 p(0:3,*),pnew(0:3,*)
      integer i,j,dim
      integer mu,k
      real * 8 tmp(0:3)
      do k=1,dim
         do mu=0,3
            pnew(mu,k) = p(mu,k)
         enddo
      enddo

      do mu=0,3
         tmp(mu) = p(mu,j)
         pnew(mu,j) = p(mu,i)
         pnew(mu,i) = tmp(mu)
      enddo
      end


      subroutine write_p(p,num_part)
      implicit none
      real * 8 p(0:3,7)
      integer num_part
      integer i,mu
      do i=1,num_part
         write(*,*) "mom ",i,(p(mu,i),mu=0,3)
      enddo
      end

      subroutine write_plong(plong,num_part)
      implicit none
      real * 8 plong(0:7,7)
      integer num_part
      integer i,mu
      do i=1,num_part
c         write(*,"(a4,1x,i1,1x,7(g22.14))") "mom ",i,
c     -        (plong(mu,i),mu=0,7)
         write(*,"(a3,1x,i1,1x,7(f15.9,1x))") "mom",i,
     -        (plong(mu,i),mu=0,3)
         write(*,"(6x,7(f15.9,1x))") 
     -        (plong(mu,i),mu=4,7)
      enddo
      end


c Program to read numbers from strings
      subroutine reads(string,nstr,rarr,narr,karr)
      implicit none
      integer nstr, narr, karr, iend, iperiod, istart, isign, istr
      integer j, js, k
      double precision rarr
      dimension rarr(narr),isign(2)
      real * 8 num(2)
      character * (*) string
      character * 1 ch
      karr=0
c get token
      istr=1
 1    continue
c skip blanks
      if(istr.le.nstr.and.string(istr:istr).eq." ") then
         istr=istr+1
         goto 1
      endif
      if(istr.gt.nstr) goto 999
      istart=istr
c find next blank
 2    if(istr.le.nstr.and.string(istr:istr).ne." ") then
         istr=istr+1
         goto 2
      endif
      iend=istr-1 
      iperiod=0
c value
      num(1)=0
c exponent
      num(2)=0
      k=1
      js=istart
 10   if(string(js:js).eq."-") then
         isign(k)=-1
         js=js+1
      elseif(string(js:js).eq."+") then
         isign(k)=1
         js=js+1
      else
         isign(k)=1
      endif
      do j=js,iend
         ch=string(j:j)
         if(ch.le."9".and.ch.ge."0") then
            num(k)=num(k)*10+ichar(ch)-ichar("0")
         elseif(ch.eq.".") then
            if(iperiod.ne.0.or.k.eq.2)goto 998
            iperiod=j-iend
         elseif(ch.eq."e".or.ch.eq."E".or.ch.eq."d".or.ch.eq."D")then
            if(j.eq.1.or.k.eq.2) goto 998
            if(iperiod.ne.0) iperiod=iperiod+iend-j+1
            k=2
            js=j+1
            goto 10
         else
            goto 999
         endif
       enddo
       karr=karr+1
       rarr(karr)=isign(1)*num(1)*(10.d0)**(isign(2)*num(2)+iperiod)
       if(karr.eq.narr) goto 999
       if(iend.lt.nstr) goto 1
       goto 999
 998   continue
       stop
 999   end

c Program to read integers from strings
      subroutine ireads(string,nstr,iarr,narr,karr)
      implicit none
      integer nstr, iarr, narr, karr, iend, isign, istart
      integer istr, num, j, js
      dimension iarr(narr)
      character * (*) string
      character * 1 ch
      karr=0
c get token
      istr=1
 1    continue
c skip blanks
      if(string(istr:istr).eq." ".and.istr.le.nstr) then
         istr=istr+1
         goto 1
      endif
      if(istr.gt.nstr) goto 999
      istart=istr
c find next blank
 2    if(string(istr:istr).ne." ".and.istr.le.nstr) then
         istr=istr+1
         goto 2
      endif
      iend=istr-1 
c value
      num=0
      js=istart
      if(string(js:js).eq."-") then
         isign=-1
         js=js+1
      elseif(string(js:js).eq."+") then
         isign=1
         js=js+1
      else
         isign=1
      endif
      do j=js,iend
         ch=string(j:j)
         if(ch.le."9".and.ch.ge."0") then
            num=num*10+ichar(ch)-ichar("0")
         else
            goto 999
         endif
       enddo
       karr=karr+1
       iarr(karr)=isign*num
       if(karr.eq.narr) goto 999
       if(iend.lt.nstr) goto 1
       goto 999
       continue
       stop
 999   end


c program to read numbers from input line
       subroutine iread(iun,iarr,nel,nread)
       implicit none
       integer iun, iarr, nel, nread
       dimension iarr(40)
       character * 80 str
       read(iun,"(a)") str
       call ireads(str,80,iarr,nel,nread)
       end


      subroutine read_all_results()
      implicit none
      integer event1,event2
      integer k
      real * 8 r1,r2


      open(unit=10,file="all_results.txt",status="old")
      do k=1,1000000
c         read(10,"(i4,g22.14,9x,i4,g22.14)") ! ,err=888
c     #        event1,r1, event2,r2
c         write(*,*) event1,r1,r2
         read(10,*,err=888) 
     #        event1,r1, event2,r2
c         write(*,*) event1,r1,r2
         if (((r1.eq.0d0).and.(r2.ne.0d0)).or.
     #        ((r1.ne.0d0).and.(r2.eq.0d0))) then
            write(*,"(a17,i4,2g22.14)") "ONLY ONE is ZERO ",event1,
     #           r1,r2
         endif
         if (r1*r2.lt.0d0) then
            write(*,"(a17,i4,2g22.14)") "WRONG SIGN ",event1,
     #           r1,r2
         endif         
      enddo
 888  continue
      close(10)
      end


c if quark_type=1 it computes the counterterm with incoming quark in the Born
c if quark_type=2 it computes the counterterm with incoming antiquark 
c     in the Born
c p1 = in gluon
c p2 = in c 
c p3 = out u
c p4 = out s
c p5 = out e+
c p6 = out ve
c p7 = out ubar
      subroutine rescaled_momenta_g_old(p,quark_type,rescaledp,Q2,x,z1)
      implicit none
      real * 8 p(0:3,7)
      integer quark_type

      real * 8 rescaledp(0:3,6)
      real * 8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3),p7(0:3)
      real * 8 q(0:3)
      real * 8 Q2,x,z1
      real * 8 dotrr
      integer mu,i

      call from_p_to_p7(p,p1,p2,p3,p4,p5,p6,p7)
      do i=1,6
         do mu=0,3
            rescaledp(mu,i) = p(mu,i)
         enddo
      enddo
         

      if (quark_type.eq.1) then
c     quark interacting with the rest of the diagram
         do mu=0,3
            q(mu) = p2(mu)-p4(mu)-p5(mu)-p6(mu)
         enddo
         
         Q2 = - dotrr(q,q)
         x = Q2/(2*dotrr(p1,q))
         z1 = dotrr(p3,p1)/dotrr(p1,q)
         
c         write(*,*) "uu line: x, z1 ", x, z1
         
         do mu=0,3
c     rescale incoming momentum
            rescaledp(mu,1) = x*p(mu,1)
c     rescale outgoing momentum
            rescaledp(mu,3) = q(mu) + x*p(mu,1)
         enddo     
         
      elseif (quark_type.eq.2) then
         
c     antiquark intereacting with the rest of the diagram

         do mu=0,3
            q(mu) = p2(mu)-p4(mu)-p5(mu)-p6(mu)
         enddo
         
         Q2 = - dotrr(q,q)
         x = Q2/(2*dotrr(p1,q))
         z1 = dotrr(p7,p1)/dotrr(p1,q)
         
         do mu=0,3
c     rescale incoming momentum
            rescaledp(mu,1) = x*p(mu,1)
c     rescale outgoing momentum
            rescaledp(mu,3) = q(mu) + x*p(mu,1)
         enddo     
         
      else
         write(*,*) "ONLY two quark lines in u_and_ubarc_usepve_count"
         stop
      endif                     
      end
      

c return a random number between 0 and 1
      function rand_num()
      implicit none
      real * 8 rand_num
      real * 8 random2
      integer num,num2
      COMMON/SEED/num,num2
c set random seeds
      data NUM/12345/,NUM2/67890/
      rand_num = random2(num,num2)
      end
      
c return a complex random number with real and imaginary part between 0 and 1
      function crand_num()
      implicit none
      complex * 16 crand_num
      real * 8 rand_num, re, im
      re = rand_num()
      im = rand_num()
      crand_num = dcmplx(re, im)
      end
      
