c---------------------------------------------------------------------------
c	lepton helicities have to be summed over externally 
c		( this routine is for one fixed set of lep. helicities only ) 
c---------------------------------------------------------------------------
      
      subroutine qqsp2aaqq(pbar,sign, nlo, L,bos,
     1                  uucc,uuss,ddcc,ddss,udsc,ducs)
      implicit none

C  qqsp2aaqq calculates the matrix elements**2 for 
C        q1 q3 (-> spin-2) ->   q2 q4 gamma gamma
C for process 240: AAjj

C This subroutine is analogue to qqzzqq.F, but without boxes and pentagons
C since they do not occur in process AAjj 

C by Jessica Frank

#include "VBFNLO/utilities/global.inc"
#include "tensorz.inc"
#include "tensorzn.inc"
c for check of individual graphs :
#include "tampz-lists.inc" 
#include "tensorsp2aa.inc"
c
c electroweak couplings are taken from KOPPLN
c
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
c
c variables for the main part of the program
c
      double precision  pbar(0:3,4+max_v)
      double precision uucc,uuss,ddcc,ddss,udsc,ducs,res(6),resv(6)
      double precision  p(0:3,4+max_v), p21(0:4), p43(0:4)

      integer  sign(4+max_v), nlo, mu, i, k,
     1         isig1, isig3,bos
      integer  ifl(4,6), js1, js3, L, Lold, is1, is3
      double complex prop21(4), prop43(4)
      double complex mat(6,-1:1,-1:1,9), matv(6,-1:1,-1:1,9)
      double complex mm(6,-1:1,-1:1), 
     1               mv12(6,-1:1,-1:1), mv34(6,-1:1,-1:1)

      double complex maa, maz, mza, mzz, mzz5, mzz6

      double complex psi(2,-1:1,4), jqq(0:5,-1:1,2) !, eps(0:5)

      double precision dummy(0:4)
      double complex zm2i(2:4)

      double complex contract_Tjj, dotcc, dotrc, dotqj, s1c
      external contract_Tjj, dotcc, dotrc, dotqj, s1c
      logical ldebug, ldebugv, ldebugvf,linit,lerror,lgc(4)
      data linit /.true./, lerror /.false./, ldebug /.false./
      data lgc /4*.false./
      
      logical vvdebug,vbdebug, vlldebug,bbdebug,pdebug
      parameter (vvdebug = .false.,vbdebug = .false.,vlldebug =.false.,
     #		 bbdebug = .false.,pdebug = .false.)
      
      integer*8 icb1, icount1, icb2, icount2
      data icb1/0/,icount1/0/,icb2/0/,icount2/0/,Lold/1/

      parameter (ldebugv=.false., ldebugvf=.false.)
      save ifl, Lold, zm2i, linit, lgc, icb1, icb2, icount1, icount2
      double complex  zero
      parameter (zero = (0d0,0d0) )


c
cc
c variables for virtual corrections
c
	logical bad_gauge,bad_gauge_sin,bad_gauge_isig
	common /gauge / bad_gauge,bad_gauge_sin,bad_gauge_isig
c
      double precision c2,c2o4pi     !,pi2o3, cvirtc
c      parameter (pi2o3=pi**2/3d0, cvirtc=pi2o3-7d0)
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      logical lnlo, lbox, lpent, lpt, lpq
c      parameter (lnlo=.true.)
      lnlo = NLO.ne.0    ! include some virtual stuff if T
      lbox = NLO.eq.1 .or. NLO.eq.-4 .or. NLO.eq.5    ! call boxline if T
      lpt  = NLO.eq.1 .or. NLO.le.-5                  ! true pentagon contributions on
      lpq  = NLO.eq.5 .or. NLO.eq.-4                  ! q^mu/m_V terms of Pentagon contributions only
      lpent= lpt .or. lpq                             ! include pentagons if T
c
	bad_gauge = .false. ! set F at beginning of each run
c
c
c define flavors of external quarks for the 4 NC and 2 CC subprocesses
c
      if (linit) then
         linit = .false.
         k = 1                  ! uucc
         ifl(1,k) = 3
         ifl(2,k) = 3
         ifl(3,k) = 3
         ifl(4,k) = 3
         k = 2                  ! uuss
         ifl(1,k) = 3
         ifl(2,k) = 3
         ifl(3,k) = 4
         ifl(4,k) = 4
         k = 3                  ! ddcc
         ifl(1,k) = 4
         ifl(2,k) = 4
         ifl(3,k) = 3
         ifl(4,k) = 3
         k = 4                  ! ddss
         ifl(1,k) = 4
         ifl(2,k) = 4
         ifl(3,k) = 4
         ifl(4,k) = 4
         k = 5                  ! udsc
         ifl(1,k) = 3
         ifl(2,k) = 4
         ifl(3,k) = 4
         ifl(4,k) = 3
         k =6                   ! ducs
         ifl(1,k) = 4
         ifl(2,k) = 3
         ifl(3,k) = 3
         ifl(4,k) = 4
         zm2i(2) = 1/dcmplx(xm2(2),-xmg(2))
         zm2i(3) = 1/dcmplx(xm2(3),-xmg(3))
         zm2i(4) = 1/dcmplx(xm2(4),-xmg(4))
      endif

 701  if (ldebug.or.ldebugv) then
         print*," " 
         print*," enter qqsp2aaqq "
      endif
      do k = 1,6
         do isig1 = -1,1,2
            do isig3 = -1,1,2
               do i = 1,9
                  mat(k,isig1,isig3,i) = 0
                  matv(k,isig1,isig3,i) = 0
               enddo
            enddo
         enddo
      enddo
c
c identify fermion line sign factors
c
      is1 = sign(1)
      is3 = sign(3)
      js1 = (3+sign(1))/2       ! 1 for sign1=-1,2 for sign1=+1
      js3 = (7+sign(3))/2       ! 3 for sign3=-1,4 for sign3=+1
      if (ldebugv) then
         print*," fsign = ",(sign(i),i=1,4)
         print*," is1, is3 = ",is1,is3," js1, js3 = ",js1,js3
         print*," lzs = ",lzs," lzs1 = ",lzs(js1), " lzs3 = ",lzs(js3) 
         if (sign(1).ne.sign(2) .or. sign(3).ne.sign(4) ) then
            print*," quark sign factors ",(sign(i),i=1,4)
            print*," are not implemented yet in qqsp2aaqq"
            stop
         endif
      endif

c
c define the internal momenta
c
      do mu = 0,3
         do i = 1,4+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo

            p21(mu) =   p(mu,2) - p(mu,1)
            p43(mu) =   p(mu,4) - p(mu,3)

      enddo
      p21(4) = p21(0)**2 - p21(1)**2 - p21(2)**2 - p21(3)**2
      p43(4) = p43(0)**2 - p43(1)**2 - p43(2)**2 - p43(3)**2
	 

      if (lwarn) then
         do mu=0,3
            dummy(mu) = p(mu,1)+p(mu,3)+p(mu,5)+p(mu,7)-
     1                 (p(mu,2)+p(mu,4)+p(mu,6)+p(mu,8))
         enddo
         if (abs(dummy(0)).gt.1d-8 .or.
     1       abs(dummy(1)).gt.1d-8 .or.
     2       abs(dummy(2)).gt.1d-8 .or.
     3       abs(dummy(3)).gt.1d-8 ) then
            print*," momentum mismatch in qqsp2aaqq "
            print 102," dummy = ",(dummy(mu),mu=0,3)
            print 103," sign  = ",sign
 101        format( " p(", i1, ") = ", 4(f10.3, 2x) )
 102        format(a,4f10.3)
 103        format(a,6i5)
            do i = 1,8
               write(6,101) i, p(0,i), p(1,i), p(2,i), p(3,i)
            end do
         endif 
      endif
c
c get the vector boson propagator factors
c
      prop21(1) = 1/p21(4)
      prop21(2) = 1/dcmplx(p21(4)-xm2(2),xmg(2))
      prop21(3) = 1/dcmplx(p21(4)-xm2(3),xmg(3))
      prop21(4) = prop21(3)

      prop43(1) = 1/p43(4)
      prop43(2) = 1/dcmplx(p43(4)-xm2(2),xmg(2))
      prop43(3) = 1/dcmplx(p43(4)-xm2(3),xmg(3))
      prop43(4) = prop43(3)

c
c get the external quark spinors (including factor sqrt(2E) )
c
      call psi0m(4,pbar(0,1),sign(1),psi)
c
c get the f-fbar currents J21^mu=jqq(mu,*,1), J43^mu=jqq(mu,*,2) 
c
      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1,1))
      call curr6(1,psi(1,-1,4),p(0,4),psi(1,-1,3),p(0,3),jqq(0,-1,2))
      
c
c -------------------------------------------------------------------
c
c contract with vvtoww tensors to get Vertex-Vertex scattering diagrams
c
      do isig1 = -1,1,2
         do isig3 = -1,1,2
            maa = contract_Tjj(sp2aaaa(0,0,L),jqq(0,isig1,1),jqq(0,isig3,2))
            maz = contract_Tjj(sp2azaa(0,0,L),jqq(0,isig1,1),jqq(0,isig3,2))
            mza = contract_Tjj(sp2zaaa(0,0,L),jqq(0,isig3,2),jqq(0,isig1,1))
            mzz = contract_Tjj(sp2zzaa(0,0,L),jqq(0,isig1,1),jqq(0,isig3,2))
            do k = 1,4
               mat(k,isig1,isig3,1) = 
     1              maa*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)
     2            + maz*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)
     3            + mza*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)
     4            + mzz*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)
               matv(k,isig1,isig3,1) = (0d0,0d0)
            enddo
	    
         enddo !isig3
      enddo !isig1
      
c----------------------      
      
      mzz6 = contract_Tjj(sp2wwaa6(0,0,L),jqq(0,-1,2),jqq(0,-1,1))
      mzz5 = contract_Tjj(sp2wwaa5(0,0,L),jqq(0,-1,1),jqq(0,-1,2))
      mat(5,-1,-1,1) = mzz5*clr(3,3,-1)**2
      mat(6,-1,-1,1) = mzz6*clr(3,3,-1)**2
      matv(5,-1,-1,1) = (0d0,0d0)
      matv(6,-1,-1,1) = (0d0,0d0)
     

c for checking individual pieces only:
	if(.false.) then
	do isig1 = -1,1,2
	   do isig3 = -1,1,2
	      do  i = 1,8
	         do k = 1,6
		    mat(k,isig1,isig3,i) = 0d0
		 enddo	
 	      enddo
	    enddo
	  enddo 
	endif
		      

      do k = 1,6
         res(k) = 0
         resv(k) = 0
         do isig1 = -1,1,2
            do isig3 = -1,1,2
               mm(k,isig1,isig3) = 0
               do i = 1,8
                  mm(k,isig1,isig3) = mm(k,isig1,isig3) + 
     1                                mat(k,isig1,isig3,i)
               enddo
	       
               res(k) = res(k) + dreal(mm(k,isig1,isig3))**2
     &                         + dimag(mm(k,isig1,isig3))**2
                if (lnlo) then
                  mv12(k,isig1,isig3) = 
     1         +  matv(k,isig1,isig3,2) + matv(k,isig1,isig3,4) + 
     2            matv(k,isig1,isig3,6) + matv(k,isig1,isig3,7) 
                  mv34(k,isig1,isig3) =
     1          + matv(k,isig1,isig3,3) + matv(k,isig1,isig3,5) +
     2            matv(k,isig1,isig3,9) + matv(k,isig1,isig3,8)


c  add Born type term and multiply by F_q = alphas*C_2/4pi
c  the factor pi^2/3+9/2 for the born term is after adding the subtraction term
c  and the counter term for the renormalization of the pdfs
                  if (nlo.gt.0) then
                     mv12(k,isig1,isig3) = als(1,1)*c2o4pi*
     1                ( mv12(k,isig1,isig3) + mm(k,isig1,isig3)*cvirt )
                     mv34(k,isig1,isig3) = als(2,1)*c2o4pi*
     1                ( mv34(k,isig1,isig3) + mm(k,isig1,isig3)*cvirt )
                  else
                     mv12(k,isig1,isig3) = 
     1                    als(1,1)*c2o4pi*mv12(k,isig1,isig3)
                     mv34(k,isig1,isig3) = 
     1                    als(2,1)*c2o4pi*mv34(k,isig1,isig3)
                  endif
                  resv(k) = resv(k) + 2*dreal(
     1                 mm(k,isig1,isig3)   *
     1              conjg( mv12(k,isig1,isig3)+mv34(k,isig1,isig3) )  )
               endif

           enddo
         enddo  
	       
         if (nlo.ge.0) then
            res(k) = (res(k)+resv(k))*9d0      ! 9 is the color sum factor
         else
            res(k) = resv(k)*9d0               ! 9 is the color sum factor
         endif
      enddo
 
      if (ldebug) then
	print*,"uucc",uucc,res(1)/uucc
	print*,"uuss",uuss,res(2)/uuss
	print*,"ddcc",ddcc,res(3)/ddcc
	print*,"ddss",ddss,res(4)/ddss
	print*,"udsc",udsc,res(5)/udsc
	print*,"ducs",ducs,res(6)/ducs
       print*
      endif
c      if (ldebugv) 
c      read(*,"(a)") text
  
      if (lwarn) then
         resv(1) = uucc
         resv(2) = uuss
         resv(3) = ddcc
         resv(4) = ddss
         resv(5) = udsc
         resv(6) = ducs
         do k = 1,6
            if ((res(k)/resv(k)-1).gt.1d-7) then
               if (.not.ldebug) then
                  print*," k =",k,":",res(k)/resv(k)-1
                  ldebug = .true.
                  lerror = .true.
                  goto 701
               endif
            endif
         enddo
         if (lerror) then
            print*," "
            lerror = .false.
            ldebug = .false.
         endif
      endif


      uucc = res(1)
      uuss = res(2)
      ddcc = res(3)
      ddss = res(4)
      udsc = res(5)
      ducs = res(6)

      Lold = L
      return
      end


