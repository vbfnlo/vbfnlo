       subroutine tens_red4_new_Re_ComRank4(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,
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
       implicit none
       real * 8 p1sq, p2sq, p3sq, p1p2, p1p3,p2p3
       complex*16 D0
       real*8 r1, r2r1, r3r2, det
       real*8 z11,z12,z13,z21,z22,z23,z31,z32,z33
       real*8 iz11,iz22,iz33,piv


       real*8 C0r_234, C0r_134, C0r_124, C0r_123
       real*8 Cijr_234(6,3), Cijr_134(6,3), Cijr_124(6,3), Cijr_123(6,3)
       real*8 C0I_234, C0I_134, C0I_124, C0I_123
       real*8 CijI_234(6,3), CijI_134(6,3), CijI_124(6,3), CijI_123(6,3)
       real*8 D0r, Dijr(22,4)
       real*8 D0I, DijI(22,4)
!       real*8 i(3)
       integer j, i1
       logical Singular
       common/DetCountML/Singular
       real*8 deter, detAbs


       r1   = p1sq
       r2r1 = p2sq+2.d0*p1p2
       r3r2 = 2.d0*(p1p3+p2p3)+p3sq

      deter =-2.d0*(-2.d0*p1p2*p1p3*p2p3 + p1p3*p1p3*p2sq + p1p2*p1p2*p3sq
     &  + p1sq*(p2p3*p2p3 - p2sq*p3sq))


      detAbs=abs(2.d0*(abs(2.d0*p1p2*p1p3*p2p3) + abs(p1p3*p1p3*p2sq)
     &   +abs(p1p2*p1p2*p3sq)
     &  + abs(p1sq*(abs(p2p3*p2p3) + abs(p2sq*p3sq)))))

ccc$      If( abs(det)/detAbs.le.1d-6)Singular=.true.

       If( abs(deter)/abs(detAbs).lt.1d-6) then
         Singular=.true.
         do j=1,4
            do i1=1,22
               Dijr(i1,j)=0d0
               DijI(i1,j)=0d0
               D0r=0d0
               D0I=0d0
           enddo
         enddo
        Return
        endif


       D0r=Dble(D0)
       D0I=DImag(D0)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	if(abs(p1sq).gt.abs(p1p2)) then
		if(abs(p1sq).gt.abs(p1p3)) then
			z11=2d0*p1sq
			iz11=1d0/z11
			z12=2d0*p1p2
			z13=2d0*p1p3
c
			z21=z12*iz11
			z31=z13*iz11
c
			z22=2d0*p2sq-z21*z12
			z32=2d0*p2p3-z31*z12
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
				iz22=1d0/z22
				z32=piv
c
				z23=2d0*p3sq
				z33=2d0*p2p3
c				det=-det
				z32=z32*iz22
	            z23=z23-z21*z13
	            z33=z33-z31*z13-z32*z23
	            iz33=1d0/z33
        	  det=-z11*z22*z33
c$$$             print*, "det1",det
       call ten4_132rank4(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,r1,r2r1,r3r2,
     &  z11,z12,z13,z21,z22,z23,z31,z32,z33,iz11,iz22,iz33,
     &                     C0r_234, C0r_134, C0r_124, C0r_123,
     &                     Cijr_234, Cijr_134, Cijr_124, Cijr_123,
     &                     C0I_234, C0I_134, C0I_124, C0I_123,
     &                     CijI_234, CijI_134, CijI_124, CijI_123,
     &                     D0R,D0I,DijR,DijI)

	RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
			else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c				iorder(1)=1
c				iorder(2)=2
c				iorder(3)=3
c

				iz22=1d0/z22
                                z23=2d0*p2p3
				z33=2d0*p3sq
				z32=z32*iz22
				z23=z23-z21*z13
				z33=z33-z31*z13-z32*z23
				iz33=1d0/z33
	  det=z11*z22*z33
c$$$         print*, "det2",det

          call ten4_123rank4(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,r1,r2r1,r3r2,
     &	  z11,z12,z13,z21,z22,z23,z31,z32,z33,iz11,iz22,iz33,
     &                     C0r_234, C0r_134, C0r_124, C0r_123,
     &                     Cijr_234, Cijr_134, Cijr_124, Cijr_123,
     &                     C0I_234, C0I_134, C0I_124, C0I_123,
     &                     CijI_234, CijI_134, CijI_124, CijI_123,
     &                     D0R,D0I,DijR,DijI)
          RETURN
         		endif
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
              	else   ! p1p3 is the pivot
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         		z11=2d0*p1p3
         		iz11=1d0/z11
         		z12=2d0*p2p3
         		z13=2d0*p3sq
c
         		z21=2d0*p1p2*iz11
         		z31=2d0*p1sq*iz11

         		z22=2d0*p2sq-z21*z12
         		z32=2d0*p1p2-z31*z12
c
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
         			iz22=1d0/z22
         			z23=2d0*p1p3
         			z33=2d0*p2p3
         			z32=z32*iz22
         			z23=z23-z21*z13
         			z33=z33-z31*z13-z32*z23
         			iz33=1d0/z33
           det=z11*z22*z33
c$$$               print*, "det3",det

       call ten4_231rank4(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,r1,r2r1,r3r2,
     &  z11,z12,z13,z21,z22,z23,z31,z32,z33,iz11,iz22,iz33,
     &                     C0r_234, C0r_134, C0r_124, C0r_123,
     &                     Cijr_234, Cijr_134, Cijr_124, Cijr_123,
     &                     C0I_234, C0I_134, C0I_124, C0I_123,
     &                     CijI_234, CijI_134, CijI_124, CijI_123,
     &                     D0R,D0I,DijR,DijI)
             RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                       else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c        		iorder(1)=3
c        		iorder(2)=2
c                       iorder(3)=1
         	        iz22=1d0/z22
                        z23=2d0*p2p3
                        z33=2d0*p1p3
                        z32=z32*iz22
         		z23=z23-z21*z13
         		z33=z33-z31*z13-z32*z23
         		iz33=1d0/z33
           det=-z11*z22*z33
c$$$               print*, "det4",det


       call ten4_321rank4(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,r1,r2r1,r3r2,
     & z11,z12,z13,z21,z22,z23,z31,z32,z33,iz11,iz22,iz33,
     &                     C0r_234, C0r_134, C0r_124, C0r_123,
     &                     Cijr_234, Cijr_134, Cijr_124, Cijr_123,
     &                     C0I_234, C0I_134, C0I_124, C0I_123,
     &                     CijI_234, CijI_134, CijI_124, CijI_123,
     &                     D0R,D0I,DijR,DijI)

      RETURN
         		endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         	endif  ! p1sq or p1p3 is the pivot
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         else  ! p1p2 or p1p3 is the pivot

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
             if(abs(p1p2).gt.abs(p1p3)) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         		z11=2d0*p1p2
         		iz11=1d0/z11
         		z12=2d0*p2sq
         		z13=2d0*p2p3

         		z21=2d0*p1sq*iz11
         		z31=2d0*p1p3*iz11

         		z22=z11-z21*z12
         		z32=z13-z31*z12
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
         			iz22=1d0/z22
         			z23=2d0*p3sq
         			z33=2d0*p1p3
         			z32=z32*iz22
         			z23=z23-z21*z13
         			z33=z33-z31*z13-z32*z23
           			iz33=1d0/z33
           det=z11*z22*z33
c$$$           print*, "det5",det

       call ten4_312rank4(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,r1,r2r1,r3r2,
     & z11,z12,z13,z21,z22,z23,z31,z32,z33,iz11,iz22,iz33,
     &                     C0r_234, C0r_134, C0r_124, C0r_123,
     &                     Cijr_234, Cijr_134, Cijr_124, Cijr_123,
     &                     C0I_234, C0I_134, C0I_124, C0I_123,
     &                     CijI_234, CijI_134, CijI_124, CijI_123,
     &                     D0R,D0I,DijR,DijI)

             	RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              		else
c        			iorder(1)=2
c        			iorder(2)=1
c        			iorder(3)=3
c
         			z23=2d0*p1p3
         			z33=2d0*p3sq
         			iz22=1d0/z22
c        			det=-det
         			z32=z32*iz22
         			z23=z23-z21*z13
         			z33=z33-z31*z13-z32*z23
         			iz33=1d0/z33
            det=-z11*z22*z33

c$$$           print*, "det6",det

      call ten4_213rank4(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,r1,r2r1,r3r2,
     &  z11,z12,z13,z21,z22,z23,z31,z32,z33,iz11,iz22,iz33,
     &                     C0r_234, C0r_134, C0r_124, C0r_123,
     &                     Cijr_234, Cijr_134, Cijr_124, Cijr_123,
     &                     C0I_234, C0I_134, C0I_124, C0I_123,
     &                     CijI_234, CijI_134, CijI_124, CijI_123,
     &                     D0R,D0I,DijR,DijI)

      RETURN
	      		endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              else
         		z11=2d0*p1p3
         		iz11=1d0/z11
         		z12=2d0*p2p3
         		z13=2d0*p3sq
c
         		z21=2d0*p1p2*iz11
         		z31=2d0*p1sq*iz11

         		z22=2d0*p2sq-z21*z12
         		z32=2d0*p1p2-z31*z12
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
         			iz22=1d0/z22
         			z23=2d0*p1p3
         			z33=2d0*p2p3
         			z32=z32*iz22
         			z23=z23-z21*z13
         			z33=z33-z31*z13-z32*z23
         			iz33=1d0/z33
           det=z11*z22*z33
c$$$           print*, "det7",det


       call ten4_231rank4(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,r1,r2r1,r3r2,
     & z11,z12,z13,z21,z22,z23,z31,z32,z33,iz11,iz22,iz33,
     &                     C0r_234, C0r_134, C0r_124, C0r_123,
     &                     Cijr_234, Cijr_134, Cijr_124, Cijr_123,
     &                     C0I_234, C0I_134, C0I_124, C0I_123,
     &                     CijI_234, CijI_134, CijI_124, CijI_123,
     &                     D0R,D0I,DijR,DijI)
              RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         		else
c        			iorder(1)=3
c        			iorder(2)=2
c        			iorder(3)=1
c
         			z23=2d0*p2p3
         			z33=2d0*p1p3
c        			det=-det

         		  iz22=1d0/z22
           	  	  z32=z32*iz22
           	  	  z23=z23-z21*z13
           	  	  z33=z33-z31*z13-z32*z23
           	  	  iz33=1d0/z33
                   det=-z11*z22*z33
c$$$         	 print*, "det8",det
       call ten4_321rank4(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,r1,r2r1,r3r2,
     &  z11,z12,z13,z21,z22,z23,z31,z32,z33,iz11,iz22,iz33,
     &                     C0r_234, C0r_134, C0r_124, C0r_123,
     &                     Cijr_234, Cijr_134, Cijr_124, Cijr_123,
     &                     C0I_234, C0I_134, C0I_124, C0I_123,
     &                     CijI_234, CijI_134, CijI_124, CijI_123,
     &                     D0R,D0I,DijR,DijI)
       RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
			endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     	    endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     	endif
        end
