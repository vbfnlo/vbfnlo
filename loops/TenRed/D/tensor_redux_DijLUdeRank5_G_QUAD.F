       subroutine tens_red4_new_Re_Com_R5_G_QUAD(m0,m1,m2,m3,
     &                  p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,
     &                     C0r_234, C0r_134, C0r_124, C0r_123,
     &                     Cijr_234, Cijr_134, Cijr_124, Cijr_123,
     &                     C0I_234, C0I_134, C0I_124, C0I_123,
     &                     CijI_234, CijI_134, CijI_124, CijI_123,
     &                     D0, D0r, D0I,Dijr,DijI)
C                tens_red3 = 3-point tensors
C     Francisco Campanario
C    Email: Francam@particle.uni-karlsruhe.de
C Date: 25/02/2010
c
c
c  determine the Passarino-Veltman tensor decomposition for the four-point
c  tensor integrals
c
c                                          d^4k
c   D0; D_mu; D_mu,nu; D_mu,nu,rho =  Int ------
c                                         (2pi)^4
c
c              1;  k_mu;   k_mu k_nu; k_mu k_nu k_rho
c      -------------------------------------------------------------------
c         [k^2-m^2][(k+p1)^2-m^2][(k+p1+p2)^2-m^2][(k+p1+p2+p3)^2-m^2]
c with
c
c   D_mu = p1_mu D11  +  p2_mu D12  +  p3_mu D13
c
c   D_munu = p1_mu p2_nu D21 + p2_mu p2_nu D22 + ...
c
c  for notation see Passarino&Veltman, NP B160 (1979) 151
c
C INPUT:  psq, pq,...                        kinematics invariants
C         C0_123 = C0(1,2,3) = C0(p1,p2)     scalar three point
C         C0_124 = C0(1,2,4) = C0(p1,p2+p3)  functions in PV notation
C         C0_134 = C0(1,3,4) = C0(p1+p2,p3)
C         C0_234 = C0(2,3,4) = C0(p2,p3)
C         Cij_123(n,m) = C_nm(1,2,3) ....    higher C_nm form factors
C                                            as in tens_red3
c         D0 = D0(p,q,l)                  scalar four point function
c
c OUTPUT: Dij(n,m) = D_nm                    form factors in the tensor
c                                            integrals a la PV
c         nm = 11, 21, 31                    ff"s for D_mu
c         nm = 12, 22, 32, 42, 52, 62, 72    ff"s for D_munu
c         nm = 13, 23, 33, ..., 93, 103, 113, 123  ff"s for D_mu,nu,rho
c
       Implicit none
       real*16 p1sq, p2sq, p3sq, p1p2, p1p3,p2p3
       complex*32 D0
       real*16 r1, r2r1, r3r2, det,r1r0
       real*16 z11,z12,z13,z21,z22,z23,z31,z32,z33
       real*16 iz11,iz22,iz33,piv
       real*16 C0r_234, C0r_134, C0r_124, C0r_123
       real*16 Cijr_234(9,4), Cijr_134(9,4), Cijr_124(9,4), Cijr_123(9,4)
       real*16 C0I_234, C0I_134, C0I_124, C0I_123
       real*16 CijI_234(9,4), CijI_134(9,4), CijI_124(9,4), CijI_123(9,4)
       real*16 D0r, Dijr(34,5)
       real*16 D0I, DijI(34,5)
       real*16 m0,m1,m2,m3,m0sq,m1sq,m2sq,m3sq
       integer Div
      common/DetCountML/Singular
       integer j,indexD(3)
c New for Gram
      logical Singular,Singular1
      common/GramSwitch/GramCS,GramDS
      logical GramCS,GramDS

      real*16 para
cFC      real*16 tempjjC,tempklC,IXC,ZmaxC
cFC      Common/Decide/tempjjC,tempklC,IXC,ZMaxC,jjinitC
cFC      Save/Decide/
cFC      Common/DDecide/tempjj,tempkl,IX,ZMax,jjinit
cFC      Save/DDecide/
cFC      Integer jjinit
cFC      real*16 tempjj,tempkl,IX,Zmax
c New variables to compare
c      real*16 musq
c      Common/musqInv/musq
cFC      complex*32 Dij(34,5)
cFC      real*16 zz11,zz12,zz13,zz22,zz23,zz33
cFC      real*16 cijrtemp(9,4),cijitemp(9,4),C0ITEMP,C0RTEMP
cFC      complex*32 tentest3_Dble,tentest3
cFC      complex*32 tentest2_Dble,tentest2
cFC      complex*32 ten3,ten2
cFC      real*16 tenI3,tenR3,tenR2,tenI2
cFC      complex*32 ten3Dble,ten2Dble
cFCc
cFC      integer j1,j2,index(5)
cFC      real*16 accuracyC(0:4,5),AccuracyD(0:5,4)
cFC      real*16 accuracyDR(22,0:4,5)
cFC      Common/Accuracy/AccuracyC,AccuracyD
      real*16 ratio
cFC      complex*32 C0234,C0134,C0124,C0123
cFC      complex*32 D01234
      real*16 cdet
c      integer i,order
      real*16 s12,s23,p4sq
c      external B0finG
c      complex*32 B0finG
c      complex*32 B023,B024,B034
      Common/ratio/para
       logical printmy
       common/mprint/printmy
       integer i
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       singular1=.false.
      s23= p2sq+p3sq+2q0*p2p3
      s12= p1sq+p2sq+2q0*p1p2
      p4sq = p1sq+p2sq+p3sq+ 2q0*(p1p2+p1p3+p2p3)


      ratio=1q-6 ! use at the end


       indexD(1)=13
       indexD(2)=7
       indexD(3)=13

       det =-2.q0*(-2.q0*p1p2*p1p3*p2p3 + p1p3*p1p3*p2sq + p1p2*p1p2*p3sq
     -  + p1sq*(p2p3*p2p3 - p2sq*p3sq))
       cdet =2.q0*(abs(-2.q0*p1p2*p1p3*p2p3) +abs( p1p3*p1p3*p2sq)
     -  +abs(p1p2*p1p2*p3sq)+abs(p1sq*p2p3*p2p3)+abs(p1sq*p2sq*p3sq))

        para=abs(det/cdet)


cFC      if(printmy) then
cFCc Subdeterminants
cFC      call dt3(p1sq,p2sq,s12,ratio)
cFC      print*, 'ratio1',ratio
cFC      call dt3(p1sq,s23,p4sq,ratio)
cFC      print*, 'ratio2',ratio
cFC      call dt3(s12,p3sq,p4sq,ratio)
cFC      print*, 'ratio3',ratio
cFC      call dt3(p2sq,p3sq,s23,ratio)
cFC      print*, 'ratio4',ratio
cFC
cFC      print*, "DET4", para
cFCc
cFC      endif

       If(para.lt.4.6q-9) then
         Singular=.true.
         singular1=.true.
c         print*, "Here1",abs(det)/abs(cdet),det,Singular, GramDS
         do i=1,5
            do j=1,34
               Dijr(j,i)=0q0
               DijI(j,i)=0q0
               D0r=0q0
               D0I=0q0
            enddo
          enddo
c         stop
c         return
         endif


       Div=0
      m0sq=m0*m0
      m1sq=m1*m1
      m2sq=m2*m2
      m3sq=m3*m3

      r1   = p1sq-m1sq
      r1r0 = r1+m0sq
      r2r1 = p2sq+2.q0*p1p2-m2sq+m1sq
      r3r2 = 2.q0*(p1p3+p2p3)+p3sq-m3sq+m2sq

      r1=r1r0

       D0r=REAL(D0)
       D0I=IMAG(D0)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	if(abs(p1sq).gt.abs(p1p2)) then
		if(abs(p1sq).gt.abs(p1p3)) then
			z11=2q0*p1sq
			iz11=1q0/z11
			z12=2q0*p1p2
			z13=2q0*p1p3
c
			z21=z12*iz11
			z31=z13*iz11
c
			z22=2q0*p2sq-z21*z12
			z32=2q0*p2p3-z31*z12
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
			if(abs(z32).gt.abs(z22)) then
c				iorder(1)=1
c				iorder(2)=3
c                                iorder(3)=2
				piv=z21
				z21=z31
				z31=piv
				piv=z22
				z22=z32
				iz22=1q0/z22
				z32=piv
c
				z23=2q0*p3sq
				z33=2q0*p2p3
c				det=-det
				z32=z32*iz22
	            z23=z23-z21*z13
	            z33=z33-z31*z13-z32*z23
	            iz33=1q0/z33
        	  det=-z11*z22*z33
c$$$             print*, "det1",det


       call ten4_132_R5_G_QUAD(Div,m0sq,m1sq,m2sq,m3sq,p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,r1,r2r1,r3r2,
     &  z11,z12,z13,z21,z22,z23,z31,z32,z33,iz11,iz22,iz33,
     &                     C0r_234, C0r_134, C0r_124, C0r_123,
     &                     Cijr_234, Cijr_134, Cijr_124, Cijr_123,
     &                     C0I_234, C0I_134, C0I_124, C0I_123,
     &                     CijI_234, CijI_134, CijI_124, CijI_123,
     &                     D0R,D0I,DijR,DijI)

c       if(Singular1)goto 300
	RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
			else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c				iorder(1)=1
c				iorder(2)=2
c				iorder(3)=3
c

				iz22=1q0/z22
                                z23=2q0*p2p3
				z33=2q0*p3sq
				z32=z32*iz22
				z23=z23-z21*z13
				z33=z33-z31*z13-z32*z23
				iz33=1q0/z33
	  det=z11*z22*z33
c$$$         print*, "det2",det

          call ten4_123_R5_G_QUAD(Div,m0sq,m1sq,m2sq,m3sq,p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,r1,r2r1,r3r2,
     &	  z11,z12,z13,z21,z22,z23,z31,z32,z33,iz11,iz22,iz33,
     &                     C0r_234, C0r_134, C0r_124, C0r_123,
     &                     Cijr_234, Cijr_134, Cijr_124, Cijr_123,
     &                     C0I_234, C0I_134, C0I_124, C0I_123,
     &                     CijI_234, CijI_134, CijI_124, CijI_123,
     &                     D0R,D0I,DijR,DijI)
c       if(Singular1)goto 300
          RETURN
         		endif
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
              	else   ! p1p3 is the pivot
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

         		z11=2q0*p1p3
         		iz11=1q0/z11
         		z12=2q0*p2p3
         		z13=2q0*p3sq
c
         		z21=2q0*p1p2*iz11
         		z31=2q0*p1sq*iz11

         		z22=2q0*p2sq-z21*z12
         		z32=2q0*p1p2-z31*z12
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
         		if(abs(z32).gt.abs(z22)) then
c        			iorder(1)=2
c        			iorder(2)=3
c                               iorder(3)=1
         			piv=z21
         			z21=z31
         			z31=piv
         			piv=z22
         			z22=z32
         			z32=piv
c
         			iz22=1q0/z22
         			z23=2q0*p1p3
         			z33=2q0*p2p3
         			z32=z32*iz22
         			z23=z23-z21*z13
         			z33=z33-z31*z13-z32*z23
         			iz33=1q0/z33
           det=z11*z22*z33
c$$$               print*, "det3",det

       call ten4_231_R5_G_QUAD(Div,m0sq,m1sq,m2sq,m3sq,p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,r1,r2r1,r3r2,
     &  z11,z12,z13,z21,z22,z23,z31,z32,z33,iz11,iz22,iz33,
     &                     C0r_234, C0r_134, C0r_124, C0r_123,
     &                     Cijr_234, Cijr_134, Cijr_124, Cijr_123,
     &                     C0I_234, C0I_134, C0I_124, C0I_123,
     &                     CijI_234, CijI_134, CijI_124, CijI_123,
     &                     D0R,D0I,DijR,DijI)
c       if(Singular1)goto 300

             RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                       else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c        		iorder(1)=3
c        		iorder(2)=2
c                       iorder(3)=1
         	        iz22=1q0/z22
                        z23=2q0*p2p3
                        z33=2q0*p1p3
                        z32=z32*iz22
         		z23=z23-z21*z13
         		z33=z33-z31*z13-z32*z23
         		iz33=1q0/z33
           det=-z11*z22*z33
c$$$               print*, "det4",det


       call ten4_321_R5_G_QUAD(Div,m0sq,m1sq,m2sq,m3sq,p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,r1,r2r1,r3r2,
     & z11,z12,z13,z21,z22,z23,z31,z32,z33,iz11,iz22,iz33,
     &                     C0r_234, C0r_134, C0r_124, C0r_123,
     &                     Cijr_234, Cijr_134, Cijr_124, Cijr_123,
     &                     C0I_234, C0I_134, C0I_124, C0I_123,
     &                     CijI_234, CijI_134, CijI_124, CijI_123,
     &                     D0R,D0I,DijR,DijI)
c       if(Singular1)goto 300
      RETURN
         		endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         	endif  ! p1sq or p1p3 is the pivot
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         else  ! p1p2 or p1p3 is the pivot
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
             if(abs(p1p2).gt.abs(p1p3)) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         		z11=2q0*p1p2
         		iz11=1q0/z11
         		z12=2q0*p2sq
         		z13=2q0*p2p3

         		z21=2q0*p1sq*iz11
         		z31=2q0*p1p3*iz11

         		z22=z11-z21*z12
         		z32=z13-z31*z12
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        		if(abs(z32).gt.abs(z22)) then
c        			iorder(1)=3
c        			iorder(2)=1
c                               iorder(3)=2
         			piv=z21
         			z21=z31
         			z31=piv
         			piv=z22
         			z22=z32
         			z32=piv
c
         			iz22=1q0/z22
         			z23=2q0*p3sq
         			z33=2q0*p1p3
         			z32=z32*iz22
         			z23=z23-z21*z13
         			z33=z33-z31*z13-z32*z23
           			iz33=1q0/z33
           det=z11*z22*z33
c$$$           print*, "det5",det

       call ten4_312_R5_G_QUAD(Div,m0sq,m1sq,m2sq,m3sq,p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,r1,r2r1,r3r2,
     & z11,z12,z13,z21,z22,z23,z31,z32,z33,iz11,iz22,iz33,
     &                     C0r_234, C0r_134, C0r_124, C0r_123,
     &                     Cijr_234, Cijr_134, Cijr_124, Cijr_123,
     &                     C0I_234, C0I_134, C0I_124, C0I_123,
     &                     CijI_234, CijI_134, CijI_124, CijI_123,
     &                     D0R,D0I,DijR,DijI)
c       if(Singular1)goto 300
             	RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              		else
c        			iorder(1)=2
c        			iorder(2)=1
c        			iorder(3)=3
c
         			z23=2q0*p1p3
         			z33=2q0*p3sq
         			iz22=1q0/z22
c        			det=-det
         			z32=z32*iz22
         			z23=z23-z21*z13
         			z33=z33-z31*z13-z32*z23
         			iz33=1q0/z33
            det=-z11*z22*z33

c$$$           print*, "det6",det

      call ten4_213_R5_G_QUAD(Div,m0sq,m1sq,m2sq,m3sq,p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,r1,r2r1,r3r2,
     &  z11,z12,z13,z21,z22,z23,z31,z32,z33,iz11,iz22,iz33,
     &                     C0r_234, C0r_134, C0r_124, C0r_123,
     &                     Cijr_234, Cijr_134, Cijr_124, Cijr_123,
     &                     C0I_234, C0I_134, C0I_124, C0I_123,
     &                     CijI_234, CijI_134, CijI_124, CijI_123,
     &                     D0R,D0I,DijR,DijI)
c       if(Singular1)goto 300
      RETURN
	      		endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              else
         		z11=2q0*p1p3
         		iz11=1q0/z11
         		z12=2q0*p2p3
         		z13=2q0*p3sq
c
         		z21=2q0*p1p2*iz11
         		z31=2q0*p1sq*iz11

         		z22=2q0*p2sq-z21*z12
         		z32=2q0*p1p2-z31*z12
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         		if(abs(z32).gt.abs(z22)) then
c        			iorder(1)=2
c        			iorder(2)=3
c                               iorder(3)=1
         			piv=z21
         			z21=z31
         			z31=piv
         			piv=z22
         			z22=z32
         			z32=piv
c
         			iz22=1q0/z22
         			z23=2q0*p1p3
         			z33=2q0*p2p3
         			z32=z32*iz22
         			z23=z23-z21*z13
         			z33=z33-z31*z13-z32*z23
         			iz33=1q0/z33
           det=z11*z22*z33
c$$$           print*, "det7",det


       call ten4_231_R5_G_QUAD(Div,m0sq,m1sq,m2sq,m3sq,p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,r1,r2r1,r3r2,
     & z11,z12,z13,z21,z22,z23,z31,z32,z33,iz11,iz22,iz33,
     &                     C0r_234, C0r_134, C0r_124, C0r_123,
     &                     Cijr_234, Cijr_134, Cijr_124, Cijr_123,
     &                     C0I_234, C0I_134, C0I_124, C0I_123,
     &                     CijI_234, CijI_134, CijI_124, CijI_123,
     &                     D0R,D0I,DijR,DijI)
c       if(Singular1)goto 300

              RETURN

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         		else
c        			iorder(1)=3
c        			iorder(2)=2
c        			iorder(3)=1
c
         			z23=2q0*p2p3
         			z33=2q0*p1p3
c        			det=-det

         		  iz22=1q0/z22
           	  	  z32=z32*iz22
           	  	  z23=z23-z21*z13
           	  	  z33=z33-z31*z13-z32*z23
           	  	  iz33=1q0/z33
                   det=-z11*z22*z33
c$$$         	 print*, "det8",det


       call ten4_321_R5_G_QUAD(Div,m0sq,m1sq,m2sq,m3sq,p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,r1,r2r1,r3r2,
     &  z11,z12,z13,z21,z22,z23,z31,z32,z33,iz11,iz22,iz33,
     &                     C0r_234, C0r_134, C0r_124, C0r_123,
     &                     Cijr_234, Cijr_134, Cijr_124, Cijr_123,
     &                     C0I_234, C0I_134, C0I_124, C0I_123,
     &                     CijI_234, CijI_134, CijI_124, CijI_123,
     &                     D0R,D0I,DijR,DijI)
c       if(Singular1)goto 300

       RETURN

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
			endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     	    endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     	endif

cFC 300   if(GramDS.and.singular1) then
cFC           order=8
cFC        if(m0*m1*m2.eq.0d0) then
cFC          C0234=QCMPLX(C0r_234,C0I_234)
cFC          C0134=QCMPLX(C0r_134,C0I_134)
cFC          C0124=QCMPLX(C0r_124,C0I_124)
cFC          C0123=QCMPLX(C0r_123,C0I_123)
cFC          D01234=D0
cFC
cFC          call ten_red4_Gram1(p1sq,p2sq,p3sq,p4sq,s12,s23,musq,C0234,
cFC     &                 C0134,C0124,C0123,D01234,Dij,order)
cFC
cFC        elseif(abs(m0**2/(m1*m2)-1d0).lt.1d-3) then
cFCc
cFC          call ten_red4_Gram1_M(m0,p1sq,p2sq,p3sq,p4sq,s12,s23,musq,C0234,
cFC     &                 C0134,C0124,C0123,D01234,Dij,order)
cFC
cFC        endif ! masses selection
cFC
cFC         if(printmy) then
cFC           print*, "D0",D0
cFC           PRINT*, "d01234gram",d01234
cFC           PRINT*
cFC           PRINT*,  "ix",IX,TEMPJJ,TEMPKL
cFC         print*, "AccuracyD", AccuracyD(0,1),AccuracyD(1,1),AccuracyD(2,1),AccuracyD(3,1),AccuracyD(4,1)
cFC         print*, "AccuracyD", AccuracyD(0,2),AccuracyD(1,2),AccuracyD(2,2),AccuracyD(3,2),AccuracyD(4,2)
cFC         print*, "AccuracyD", AccuracyD(0,3),AccuracyD(1,3),AccuracyD(2,3),AccuracyD(3,3),AccuracyD(4,3)
cFC         endif
cFC
cFC       ZZ11=2d0*p1sq
cFC       ZZ12=2d0*p1p2
cFC       ZZ13=2d0*p1p3
cFC       ZZ22=2d0*p2sq
cFC       ZZ23=2d0*p2p3
cFC       ZZ33=2d0*p3sq
cFC
cFCc	B023=B0fin(p2sq,musq)
cFCc	B024=B0fin(s23,musq)
cFCc     	B034=B0fin(p3sq,musq)
cFCc        tentest3=4d0*(B023-B024-(r3r2)*C0234)
cFCc      tentest2=4d0*(B024-B034-(r2r1)*C0234)
cFC
cFC       B023=B0finG(m1,m2,p2sq,musq)
cFC       B024=B0finG(m1,m3,s23,musq)
cFC       B034=B0finG(m2,m3,p3sq,musq)
cFC
cFC       tentest3_Dble=4d0*(B023-B024-(r3r2)*C0234+m0sq*(C0123-C0124-r3r2*D0))
cFC       tentest2_Dble=4d0*(B024-B034-(r2r1)*C0234 + m0sq*(C0124-C0134-r2r1*D0 ))
cFC       tentest3=4d0*(B023-B024-(r3r2)*C0234+m0sq*(C0123-C0124-r3r2*D01234))
cFC       tentest2=4d0*(B024-B034-(r2r1)*C0234 + m0sq*(C0124-C0134-r2r1*D01234))
cFC
cFC       ten3= 2*(ZZ33**2*Dij(3,3)+ZZ11*(ZZ13*Dij(1,3)+ZZ23*
cFC     -   Dij(4,3)+ZZ33*Dij(5,3))+ZZ22*(ZZ23*Dij(2,3)+ZZ13*
cFC     -   Dij(6,3)+ZZ33*Dij(8,3))+3*(ZZ13*ZZ33*Dij(7,3)+ZZ2
cFC     -   3*ZZ33*Dij(9,3))+4*ZZ13*ZZ23*Dij(10,3)+2*(ZZ12*ZZ13*Dij
cFC     -   (4,3)+ZZ13**2*Dij(5,3)+ZZ12*ZZ23*Dij(6,3)+ZZ23**2
cFC     -   *Dij(8,3)+ZZ12*ZZ33*Dij(10,3))+12*(ZZ13*Dij(11,3
cFC     -   )+ZZ23*Dij(12,3)+ZZ33*Dij(13,3)))
cFC
cFC
cFC       tenI3= 2*(ZZ33**2*DijI(3,3)+ZZ11*(ZZ13*DijI(1,3)+ZZ23*
cFC     -   DijI(4,3)+ZZ33*DijI(5,3))+ZZ22*(ZZ23*DijI(2,3)+ZZ13*
cFC     -   DijI(6,3)+ZZ33*DijI(8,3))+3*(ZZ13*ZZ33*DijI(7,3)+ZZ2
cFC     -   3*ZZ33*DijI(9,3))+4*ZZ13*ZZ23*DijI(10,3)+2*(ZZ12*ZZ13*Dij
cFC     -   I(4,3)+ZZ13**2*DijI(5,3)+ZZ12*ZZ23*DijI(6,3)+ZZ23**2
cFC     -   *DijI(8,3)+ZZ12*ZZ33*DijI(10,3))+12*(ZZ13*DijI(11,3
cFC     -   )+ZZ23*DijI(12,3)+ZZ33*DijI(13,3)))
cFC
cFC       tenR3= 2*(ZZ33**2*DijR(3,3)+ZZ11*(ZZ13*DijR(1,3)+ZZ23*
cFC     -   DijR(4,3)+ZZ33*DijR(5,3))+ZZ22*(ZZ23*DijR(2,3)+ZZ13*
cFC     -   DijR(6,3)+ZZ33*DijR(8,3))+3*(ZZ13*ZZ33*DijR(7,3)+ZZ2
cFC     -   3*ZZ33*DijR(9,3))+4*ZZ13*ZZ23*DijR(10,3)+2*(ZZ12*ZZ13*Dij
cFC     -   R(4,3)+ZZ13**2*DijR(5,3)+ZZ12*ZZ23*DijR(6,3)+ZZ23**2
cFC     -   *DijR(8,3)+ZZ12*ZZ33*DijR(10,3))+12*(ZZ13*DijR(11,3
cFC     -   )+ZZ23*DijR(12,3)+ZZ33*DijR(13,3)))
cFC
cFC
cFC	ten2= 2*(ZZ22**2*Dij(2,3)+ZZ23*ZZ33*Dij(3,3)+ZZ11*(
cFC     -   ZZ12*Dij(1,3)+ZZ22*Dij(4,3)+ZZ23*Dij(5,3))+ZZ12*ZZ
cFC     -   33*Dij(7,3)+2*(ZZ12**2*Dij(4,3)+ZZ12*ZZ13*Dij(5,
cFC     -   3)+ZZ13*ZZ23*Dij(7,3)+ZZ23**2*Dij(9,3))+4*ZZ12*ZZ23*Dij
cFC     -   (10,3)+ZZ22*(3*(ZZ12*Dij(6,3)+ZZ23*Dij(8,3))+ZZ33
cFC     -   *Dij(9,3)+2*ZZ13*Dij(10,3)+12*Dij(12,3))+12*(ZZ
cFC     -   12*Dij(11,3)+ZZ23*Dij(13,3)))
cFC
cFC
cFC       tenR2= 2*(ZZ22**2*DijR(2,3)+ZZ23*ZZ33*DijR(3,3)+ZZ11*(
cFC     -   ZZ12*DijR(1,3)+ZZ22*DijR(4,3)+ZZ23*DijR(5,3))+ZZ12*ZZ
cFC     -   33*DijR(7,3)+2*(ZZ12**2*DijR(4,3)+ZZ12*ZZ13*DijR(5,
cFC     -   3)+ZZ13*ZZ23*DijR(7,3)+ZZ23**2*DijR(9,3))+4*ZZ12*ZZ23*Dij
cFC     -   R(10,3)+ZZ22*(3*(ZZ12*DijR(6,3)+ZZ23*DijR(8,3))+ZZ33
cFC     -   *DijR(9,3)+2*ZZ13*DijR(10,3)+12*DijR(12,3))+12*(ZZ
cFC     -   12*DijR(11,3)+ZZ23*DijR(13,3)))
cFC
cFC
cFC       tenI2= 2*(ZZ22**2*DijI(2,3)+ZZ23*ZZ33*DijI(3,3)+ZZ11*(
cFC     -   ZZ12*DijI(1,3)+ZZ22*DijI(4,3)+ZZ23*DijI(5,3))+ZZ12*ZZ
cFC     -   33*DijI(7,3)+2*(ZZ12**2*DijI(4,3)+ZZ12*ZZ13*DijI(5,
cFC     -   3)+ZZ13*ZZ23*DijI(7,3)+ZZ23**2*DijI(9,3))+4*ZZ12*ZZ23*Dij
cFC     -   I(10,3)+ZZ22*(3*(ZZ12*DijI(6,3)+ZZ23*DijI(8,3))+ZZ33
cFC     -   *DijI(9,3)+2*ZZ13*DijI(10,3)+12*DijI(12,3))+12*(ZZ
cFC     -   12*DijI(11,3)+ZZ23*DijI(13,3)))
cFC
cFC
cFC        if(printmy) then
cFC       print*, "ten3",ten3, tentest3
cFC       print*, "ten2",ten2, tentest2
cFC       print*, "ten3_dble",QCMPLX(tenR3,tenI3), tentest3_dble
cFC       print*, "ten2_dble",QCMPLX(tenR2,tenI2), tentest2_dble
cFC         endif
cFC
cFC	If(abs(tentest3).gt.1d-6) then
cFC	ten3=ten3/tentest3-1d0
cFC	ten3Dble=QCMPLX(tenR3,tenI3)/tentest3_dble-1d0
cFC	else
cFC	ten3=ten3-tentest3
cFC	ten3Dble=QCMPLX(tenR3,tenI3)-tentest3_dble
cFC	endif
cFC
cFC	If(abs(tentest2).gt.1d-6) then
cFC	ten2=ten2/tentest2-1d0
cFC	ten2Dble=QCMPLX(tenR2,tenI2)/tentest2_Dble-1d0
cFC	else
cFC	ten2=ten2-tentest2
cFC	ten2Dble=QCMPLX(tenR2,tenI2)-tentest2_Dble
cFC	endif
cFC
cFC        if(printmy) then
cFC          print*, "TEST",tentest3,tentest2
cFC         print*, "TEST_dble",tentest3_dble,tentest2_dble
cFC         print*, "DiGram",ten3,ten2
cFC         endif
cFC
cFC	ten3=Max(Abs(ten3),Abs(ten2))
cFC	ten3Dble=Max(Abs(ten3Dble),Abs(ten2Dble))
cFC
cFC
cFC        if(printmy) then
cFC          print*, "DiGram",ten3,ten3Dble
cFC        endif
cFC
cFC
cFC
cFCc        print*, "DiGram",ten3,ten3Dble
cFC
cFCcFC         if(abs(ten3).lt.abs(ten3dble)) then
cFCcFC         write(*,'(A12,10E12.3)'),"GOOD RESCUED",abs(ten3),abs(ten3Dble)
cFCcFC     & ,det,abs(det)/abs(4d0*cdet),abs(IX*det),accuracyD(4,1),accuracyD(3,1)
cFCcFC     &   ,accuracyD(2,1),accuracyD(1,1),accuracyD(0,1)
cFCcFC         else
cFCcFC         write(*,'(A12,10E12.3)'),"BAD RESCUED",abs(ten3),abs(ten3Dble)
cFCcFC     & ,det,abs(det)/abs(4d0*cdet),abs(IX*det),accuracyD(4,1),accuracyD(3,1)
cFCcFC     & ,accuracyD(2,1),accuracyD(1,1),accuracyD(0,1)
cFCcFC         endif
cFC
cFCc         pause
cFC
cFC
cFC
cFC       index(1)=3
cFC       index(2)=7
cFC       index(3)=13
cFC       index(4)=22
cFC       index(5)=34
cFC
cFC
cFC
cFC
cFC
cFC         if(abs(ten3).lt.(abs(ten3dble))) then
cFC
cFC         do j1=1,5
cFC          do j2= 1,index(j1)
cFC             DijR(j2,j1)=Dble(Dij(j2,j1))
cFC             DijI(j2,j1)=DImag(Dij(j2,j1))
cFC          enddo
cFC        enddo
cFC
cFC
cFC
cFC
cFC        D0r=Dble(D01234)
cFC        D0I=DImag(D01234)
cFC        if(printmy) then
cFC           print*, "D0",D0
cFC        endif
cFC        D0=D01234
cFC        if(printmy) then
cFC           print*, "D0_cp",D0
cFC        endif
cFC
cFC
cFCc        C0=C0123
cFC
cFC       tentest3_Dble=4d0*(B023-B024-(r3r2)*C0234+m0sq*(C0123-C0124-r3r2*D0))
cFC       tentest2_Dble=4d0*(B024-B034-(r2r1)*C0234 + m0sq*(C0124-C0134-r2r1*D0 ))
cFC
cFC
cFC       tenI3= 2*(ZZ33**2*DijI(3,3)+ZZ11*(ZZ13*DijI(1,3)+ZZ23*
cFC     -   DijI(4,3)+ZZ33*DijI(5,3))+ZZ22*(ZZ23*DijI(2,3)+ZZ13*
cFC     -   DijI(6,3)+ZZ33*DijI(8,3))+3*(ZZ13*ZZ33*DijI(7,3)+ZZ2
cFC     -   3*ZZ33*DijI(9,3))+4*ZZ13*ZZ23*DijI(10,3)+2*(ZZ12*ZZ13*Dij
cFC     -   I(4,3)+ZZ13**2*DijI(5,3)+ZZ12*ZZ23*DijI(6,3)+ZZ23**2
cFC     -   *DijI(8,3)+ZZ12*ZZ33*DijI(10,3))+12*(ZZ13*DijI(11,3
cFC     -   )+ZZ23*DijI(12,3)+ZZ33*DijI(13,3)))
cFC
cFC       tenR3= 2*(ZZ33**2*DijR(3,3)+ZZ11*(ZZ13*DijR(1,3)+ZZ23*
cFC     -   DijR(4,3)+ZZ33*DijR(5,3))+ZZ22*(ZZ23*DijR(2,3)+ZZ13*
cFC     -   DijR(6,3)+ZZ33*DijR(8,3))+3*(ZZ13*ZZ33*DijR(7,3)+ZZ2
cFC     -   3*ZZ33*DijR(9,3))+4*ZZ13*ZZ23*DijR(10,3)+2*(ZZ12*ZZ13*Dij
cFC     -   R(4,3)+ZZ13**2*DijR(5,3)+ZZ12*ZZ23*DijR(6,3)+ZZ23**2
cFC     -   *DijR(8,3)+ZZ12*ZZ33*DijR(10,3))+12*(ZZ13*DijR(11,3
cFC     -   )+ZZ23*DijR(12,3)+ZZ33*DijR(13,3)))
cFC
cFC
cFC       tenR2= 2*(ZZ22**2*DijR(2,3)+ZZ23*ZZ33*DijR(3,3)+ZZ11*(
cFC     -   ZZ12*DijR(1,3)+ZZ22*DijR(4,3)+ZZ23*DijR(5,3))+ZZ12*ZZ
cFC     -   33*DijR(7,3)+2*(ZZ12**2*DijR(4,3)+ZZ12*ZZ13*DijR(5,
cFC     -   3)+ZZ13*ZZ23*DijR(7,3)+ZZ23**2*DijR(9,3))+4*ZZ12*ZZ23*Dij
cFC     -   R(10,3)+ZZ22*(3*(ZZ12*DijR(6,3)+ZZ23*DijR(8,3))+ZZ33
cFC     -   *DijR(9,3)+2*ZZ13*DijR(10,3)+12*DijR(12,3))+12*(ZZ
cFC     -   12*DijR(11,3)+ZZ23*DijR(13,3)))
cFC
cFC
cFC       tenI2= 2*(ZZ22**2*DijI(2,3)+ZZ23*ZZ33*DijI(3,3)+ZZ11*(
cFC     -   ZZ12*DijI(1,3)+ZZ22*DijI(4,3)+ZZ23*DijI(5,3))+ZZ12*ZZ
cFC     -   33*DijI(7,3)+2*(ZZ12**2*DijI(4,3)+ZZ12*ZZ13*DijI(5,
cFC     -   3)+ZZ13*ZZ23*DijI(7,3)+ZZ23**2*DijI(9,3))+4*ZZ12*ZZ23*Dij
cFC     -   I(10,3)+ZZ22*(3*(ZZ12*DijI(6,3)+ZZ23*DijI(8,3))+ZZ33
cFC     -   *DijI(9,3)+2*ZZ13*DijI(10,3)+12*DijI(12,3))+12*(ZZ
cFC     -   12*DijI(11,3)+ZZ23*DijI(13,3)))
cFC
cFC	If(abs(tentest3).gt.1d-6) then
cFC	ten3Dble=QCMPLX(tenR3,tenI3)/tentest3_dble-1d0
cFC	else
cFC	ten3Dble=QCMPLX(tenR3,tenI3)-tentest3_dble
cFC	endif
cFC
cFC	If(abs(tentest2).gt.1d-6) then
cFC	ten2Dble=QCMPLX(tenR2,tenI2)/tentest2_Dble-1d0
cFC	else
cFC	ten2Dble=QCMPLX(tenR2,tenI2)-tentest2_Dble
cFC	endif
cFC
cFC	ten3Dble=Max(Abs(ten3Dble),Abs(ten2Dble))
cFC
cFC        if(printmy) then
cFC        print*, "here",ten3/ten3Dble
cFC        endif
cFC
cFC        Singular=.false.
cFC        return
cFC
cFC        else
cFC
cFC         if(abs(ten3).lt.(abs(20*ten3dble))) then
cFC            if(printmy) then
cFC            PRINT*,"TEN3 LT 20*TEN3DLE Rank 4,5 rewritten"
cFC            endif
cFC          do j1=4,5
cFC          do j2= 1,index(j1)
cFC             DijR(j2,j1)=Dble(Dij(j2,j1))
cFC             DijI(j2,j1)=DImag(Dij(j2,j1))
cFC          enddo
cFC        enddo
cFC        endif
cFC
cFC
cFC         Singular=.false.
cFC         If( para.lt.ratio) then
cFC         Singular=.true.
cFCc         singular1=.true.
cFCc         print*, "Here1HereCCRank4",abs(deter)/abs(detAbs),Singular, GramCS
cFC         do i=1,5
cFC            do j=1,34
cFC               Dijr(j,i)=0d0
cFC               DijI(j,i)=0d0
cFC               D0r=0d0
cFC               D0I=0d0
cFC            enddo
cFC         enddo
cFCc         stop
cFC         return
cFC         endif
cFC         endif
cFC
cFC
cFC
cFC         else ! GramDS singular 1
cFC
cFC         Singular=.false.
cFC         If( para.lt.ratio) then
cFC         Singular=.true.
cFCc         singular1=.true.
cFCc         print*, "Here2HereDDRank4",abs(deter)/abs(detAbs),Singular, GramCS
cFC         do i=1,5
cFC            do j=1,34
cFC               Dijr(j,i)=0d0
cFC               DijI(j,i)=0d0
cFC               D0r=0d0
cFC               D0I=0d0
cFC            enddo
cFC         enddo
cFCc         stop
cFC         return
cFC         endif
cFC
cFC	return
cFC
cFC
cFC       endif ! end GRAMDS singular1

        end

