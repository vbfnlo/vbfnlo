c this routine is used to call routines which compute pentagons 
c	(depending on NLO "true pentagons" or "boxlike" contributions 
c		are computed)
c	
c	this is a generalization of "J_virtual_tri_box_pent45" to
c	quarkline with arbitrary helicity isig
c
c
c	last modified by Barbara Jager on Sept. 21, 2005
c
c 
c--------------------------------------------------------------------------- 
      subroutine  J_virtual_tri_box_pent_isig45(psi,p,nq1,nq2, sign, musq,
c     input currents (any ORDERED couple among [J1,J2,J3]) in DHELAS format
c     They can be  jmu1,jmu2  or  jmu1,jmu3   or   jmu2,jmu3
     #     jal,jbe,
c     factor of proportionality of qal/qbe in jal = jal"+xal*qal vv.	
     #	   xal,xbe,    
c     uncontracted external current. It can be 1, 2 or 3
     #     num_ret_current,
c     specify piece to be returned: 
c          NLO = +5#include leading q^mu/m_V terms of 
c	                        pentagons which can be related to boxes,
c				set other pentagon contributions zero
c          NLO = -5      finite Pentagon contributions only (no box-type terms)
     #	   NLO,		 		     
c     returned current
     #     Jret)
c     the corresponding Born current
c     #     Jborn)

      implicit none
      double complex  psi(2,-1:1,4)
      double precision p(0:3,4)
      double precision musq
      double complex  jal(0:5),jbe(0:5),xal,xbe
      integer num_ret_current,nlo,nq1,nq2
      double complex jret(0:3,-1:1)
      integer mu,kb, sign(4),isig
      double precision qa(0:4),qb(0:4),fa(0:4,4),fb(0:4,4),
     # 		       fba(0:4,4)
      double complex psia(2,-1:1,4),psib(2,-1:1,4),
     #		     psiba(2,-1:1,4)		      
      double complex  jborn(0:5,-1:1),jret_nlo(0:3,-1:1)	      
   
      double precision  qak(0:4,2),qbk(0:4,2)	      
      double complex cak(0:5,2),cbk(0:5,2),jak(0:5,2),jbk(0:5,2),
     #		     j0ak(0:5,-1:1,2),j0bk(0:5,-1:1,2),
     #		     j0aa(0:5,-1:1,2),j0bb(0:5,-1:1,2),
     #		     jvirt1(0:5,-1:1),jvirt2(0:5,-1:1),
     #		     jvirt3(0:5,-1:1),jvirt4(0:5,-1:1),
     #		     jvabd(0:5,-1:1,3)
  
      double complex psiak(2,-1:1,2,4),psibk(2,-1:1,2,4),
     #  	     psipak(2,-1:1,2,4),psipbk(2,-1:1,2,4)
      double precision  fqak(0:4,2,4),fqbk(0:4,2,4),
     # 			fpak(0:4,2,4),fpbk(0:4,2,4)

c variables for virtual corections
c
      double precision pi,pi2o3,c2,c2o4pi,cvirtc
      parameter (pi = 3.141592653589793d0)
      parameter (pi2o3=pi**2/3d0, cvirtc=pi2o3-7d0)
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)


      ! compute (outgoing) momenta asociated with external currents:
      
      qa(0) = dreal(jal(4))
      qa(1) = dreal(jal(5))
      qa(2) = dimag(jal(5))
      qa(3) = dimag(jal(4))
      qb(0) = dreal(jbe(4))
      qb(1) = dreal(jbe(5))
      qb(2) = dimag(jbe(5))
      qb(3) = dimag(jbe(4))
      
      qa(4) = qa(0)**2-qa(1)**2-qa(2)**2-qa(3)**2
      qb(4) = qb(0)**2-qb(1)**2-qb(2)**2-qb(3)**2
 
      if (nlo.eq.-5) then       ! need "genuine" pentagon contributions
 
  ! compute Born current for pentagon line:
        	
         if (num_ret_current.eq.1) then ! mu,jal,jbe
	
     	   do isig = -1,1,2 ! hel.of quark line
              call bra2c(psi(1,isig,nq2),.true.,p(0,nq2),isig,qb,
     #           jbe,psib(1,isig,nq2),fb(0,nq2))
              call bra2c(psib(1,isig,nq2),.false.,fb(0,nq2),isig,qa,
     #           jal,psiba(1,isig,nq2),fba(0,nq2))
           enddo !isig	 
   
           call curr6(1,psiba(1,-1,nq2),fba(0,nq2),psi(1,-1,nq1),p(0,nq1),
     #           jborn(0,-1))	
     
         elseif (num_ret_current.eq.2) then ! jal,mu,jbe
	
      	   do isig = -1,1,2 ! hel.of quark line
              call bra2c(psi(1,isig,nq2),.true.,p(0,nq2),isig,qb,
     #           jbe,psib(1,isig,nq2),fb(0,nq2))
              call ket2c(psi(1,isig,nq1),.true.,p(0,nq1),isig,qa,
     #           jal,psia(1,isig,nq1),fa(0,nq1))
           enddo !isig	 
           
	   call curr6(1,psib(1,-1,nq2),fb(0,nq2),psia(1,-1,nq1),fa(0,nq1),
     #           jborn(0,-1))	
 
         elseif (num_ret_current.eq.3) then ! jal,jbe,mu
	
      	   do isig = -1,1,2 ! hel.of quark line
              call ket2c(psi(1,isig,nq1),.true.,p(0,nq1),isig,qa,
     #           jal,psia(1,isig,nq1),fa(0,nq1))
              call ket2c(psia(1,isig,nq1),.false.,fa(0,nq1),isig,qb,
     #           jbe,psiba(1,isig,nq1),fba(0,nq1))
           enddo !isig	 
           
	   call curr6(1,psi(1,-1,nq2),p(0,nq2),psiba(1,-1,nq1),fba(0,nq1),
     #           jborn(0,-1))	
 
         else 
      	
            print*,"this configuration for pentagon line is not allowed" 
      
         endif


      !  compute pentagon line including born contribution:
         
	 do isig = -1,1,2
	    
         call J_virtual_tri_box_pent2(psi,p,nq1,nq2,isig,sign, musq,
     1        jal, jbe, num_ret_current, jret_nlo(0,isig) )

	! subtract born from nlo current:

         do mu = 0,3
            jret(mu,isig) = jret_nlo(mu,isig) - cvirtc*jborn(mu,isig) 
         enddo

	enddo !isig

      endif                     ! nlo = -5

c ----------------------------------------

c now compute box-type contributions that come from contraction 
c  of momenta with pentagon :
	
      if (nlo.eq.5) then

c assume that jal stands for term in jalf = jal0 + xal*jal prop. to qa, 
c 	e.g. for jalf = wp have: wp = wpp + xp*qp 
c		-> in this part of the code xp*qp/xm+qm-type contributions 
c		   (and various mixed terms) are computed for specific 
c	            ordering of ext. currents,
c			(xp contained in jal ~ xp*qp)
c 	(analogous for jbe)
 
c
 	! need modified momenta for different contractions:
	! index kb is used to name different sums of momenta:
	!	kb = 1 ... 	qa -> qak(mu,kb=1) = qa
	!	kb = 2 ... 	qa -> qak(mu,kb=2) = qa + qb
	
         do mu = 0,3
            qak(mu,1) = qa(mu)
            qbk(mu,1) = qb(mu)
            qak(mu,2) = qa(mu)+qb(mu)
            qbk(mu,2) = qb(mu)+qa(mu) ! = qak(mu,2)
         enddo 
         do kb = 1,2	
            qak(4,kb) = qak(0,kb)**2-qak(1,kb)**2-qak(2,kb)**2-qak(3,kb)**2
            qbk(4,kb) = qbk(0,kb)**2-qbk(1,kb)**2-qbk(2,kb)**2-qbk(3,kb)**2
         enddo


c-------------------------------------------
          
       ! calculate LO and NLO currents (include gluon loop) for "boxes":
  		
	! need trick to consider internal momentum (e.g. qa + qb) 
	! different from mom. of attached polarization vector (e.g. qa):
	! define wak(mu = 0:5,kb) with mu=4:5 -> internal momentum for 
	! case kb
	
         do mu = 0,3
            do kb = 1,2
               cak(mu,kb)  = dcmplx(qa(mu)) !complex momenta for contraction in "ket"
               cbk(mu,kb)  = dcmplx(qb(mu))   
               jak(mu,kb) = jal(mu)
               jbk(mu,kb) = jbe(mu)
            enddo
         enddo 
         do kb = 1,2
            cak(4,kb)  = dcmplx(qak(0,kb),qak(3,kb)) 
            cak(5,kb)  = dcmplx(qak(1,kb),qak(2,kb))
            cbk(4,kb)  = dcmplx(qbk(0,kb),qbk(3,kb))    
            cbk(5,kb)  = dcmplx(qbk(1,kb),qbk(2,kb))
            jak(4,kb) = cak(4,kb)    
            jak(5,kb) = cak(5,kb)
            jbk(4,kb) = cbk(4,kb)
            jbk(5,kb) = cbk(5,kb)
         enddo	
			
c---------
			
  ! for LO currents ( no gluon radiated):
  
  	! generate psiwppk(1:2,isig=-1:1,kb,i=FID) with momentum fqpp(mu=0:4,kb,i)
	!    for wpp(qp) attached to boson with momentum qppk(mu=0:4,kb) 
	!    next to fermion i
	!
	! corresp. current is jwpk(mu,isig,kb,i)
	
     
c first graph (k=3): order = jal, jbe, open index (LO)

         if (num_ret_current.eq.3) then ! jal,jbe,jexch
				
	    do isig = -1,1,2				
					
            call ket2c(psi(1,isig,nq1),.true.,p(0,nq1),isig,qak(0,1),
     #           jal(0),psiak(1,isig,1,nq1),fqak(0,1,nq1))
     
            call ket2c(psi(1,isig,nq1),.true.,p(0,nq1),isig,qak(0,2),
     #           jal(0),psiak(1,isig,2,nq1),fqak(0,2,nq1))
     
            call ket2c(psi(1,isig,nq1),.true.,p(0,nq1),isig,qak(0,2),
     #           jbe(0),psibk(1,isig,2,nq1),fqbk(0,2,nq1))
    
            call ket2c(psi(1,isig,nq1),.true.,p(0,nq1),isig,qbk(0,2),
     #           cbk(0,2),psipbk(1,isig,2,nq1),fpbk(0,2,nq1))     
     
     	    enddo !isig
     
		
            call curr6(1,psi(1,-1,nq2),p(0,nq2),
     1           psiak(1,-1,1,nq1),fqak(0,1,nq1),j0ak(0,-1,1)  )
		
            call curr6(1,psi(1,-1,nq2),p(0,nq2),
     1           psiak(1,-1,2,nq1),fqak(0,2,nq1),j0ak(0,-1,2)  )

            call curr6(1,psi(1,-1,nq2),p(0,nq2),
     1           psibk(1,-1,2,nq1),fqbk(0,2,nq1),j0bk(0,-1,2)  )

            call curr6(1,psi(1,-1,nq2),p(0,nq2),
     1           psipbk(1,-1,2,nq1),fpbk(0,2,nq1),j0bb(0,-1,2) )
     
     
                ! have psiak(1:2,-1:1,kb=1:2,nq) with fqak(mu,kb,nq) 	
    		! have j0ak(mu=0:5,isig,kb=1:2)    	
    
c
c     
c k=3: order = jal, jbe, open index at NLO
     
c           call boxlinec(isigmax, psi1(2,-1:1),psi2,k1(0:3),k2(0:3),
c     1                    ld0calc, jext(0:5), pos.of open index ivmu = 1:2, 
c     2                    jborn(0:5,-1:1), jvirt(0:5,-1:1) )
 
           call boxlinec(1,psi(1,-1,nq1),psi(1,-1,nq2),p(0,nq1),p(0,nq2),
     #           .true.,jak(0,1),2,j0ak(0,-1,1),jvirt1(0,-1))

           call boxlinec(1,psi(1,-1,nq1),psi(1,-1,nq2),p(0,nq1),p(0,nq2),
     #	   	.true.,jak(0,2),2,j0ak(0,-1,2),jvirt2(0,-1))
     
      	   call boxlinec(1,psi(1,-1,nq1),psi(1,-1,nq2),p(0,nq1),p(0,nq2),
     #	   	.true.,jbk(0,2),2,j0bk(0,-1,2),jvirt3(0,-1))
    
       	   call boxlinec(1,psi(1,-1,nq1),psi(1,-1,nq2),p(0,nq1),p(0,nq2),
     #	   	.true.,cbk(0,2),2,j0bb(0,-1,2),jvirt4(0,-1))
    
                 
    ! emit jal, then jbe from line nq1-nq2: 
    
    	do isig = -1,1,2
            do mu = 0,5
               jvabd(mu,isig,3)= xbe*(jvirt1(mu,isig) - jvirt2(mu,isig)) 
     #                     + xal*jvirt3(mu,isig) + xal*xbe*jvirt4(mu,isig)    
     	! have jvabd(mu,isig = -1:1,k=1:3 for config)
            enddo               !mu
	enddo !isig
	
	
         endif                  ! k=3 configuration	
  		
c -------------

c second graph (k=2): order = jal, open index, jbe (LO)

         if (num_ret_current.eq.2) then ! jal,jexch,jbe
	
	    do isig = -1,1,2
	
            call ket2c(psi(1,isig,nq1),.true.,p(0,nq1),isig,qak(0,1),
     #           jal(0),psiak(1,isig,1,nq1),fqak(0,1,nq1))
		
            call bra2c(psi(1,isig,nq2),.true.,p(0,nq2),isig,qbk(0,1),
     #           jbe(0),psibk(1,isig,1,nq2),fqbk(0,1,nq2))

	! the third contriution gives zero analytically -> don"t compute numerically:    
	!   		call ket2c(psi(1,-1,nq1),.true.,p(0,nq1),-1,qak(0,1),
        !	#		cak(0,1),psipak(1,1,nq1),fpak(0,1,nq1))     
     
            enddo !isig


            call curr6(1,psi(1,-1,nq2),p(0,nq2),
     1           psiak(1,-1,1,nq1),fqak(0,1,nq1),j0ak(0,-1,1)  )
	
            call curr6(1,psibk(1,-1,1,nq2),fqbk(0,1,nq2),
     1           psi(1,-1,nq1),p(0,nq1),j0bk(0,-1,1)   )

c            call curr6(-1,psi(1,-1,nq2),p(0,nq2),
c     1                 psipak(1,1,nq1),fpak(0,1,nq1),j0aa(0,-1,1)   )


c k=2: order = jal, open index, jbe (NLO)
 
     	    call boxlinec(1,psi(1,-1,nq1),psi(1,-1,nq2),p(0,nq1),p(0,nq2),
     #	   	 .true.,jak(0,1),2,j0ak(0,-1,1),jvirt1(0,-1))

            call boxlinec(1,psi(1,-1,nq1),psi(1,-1,nq2),p(0,nq1),p(0,nq2),
     #	   	 .true.,jbk(0,1),1,j0bk(0,-1,1),jvirt2(0,-1))
        
c          call boxlinec(-1,psi(1,-1,nq1),psi(1,-1,nq2),p(0,nq1),p(0,nq2),
c     #	   	.true.,cak(0,1),2,j0aa(0,-1,1),jvirt3(0,-1))


    ! emit jal, then jbe from upper line: 
    	   
	   do isig = -1,1,2
            do mu = 0,5
               jvabd(mu,isig,2)= -xbe*jvirt1(mu,isig) + xal*jvirt2(mu,isig)    
            enddo               !mu
	   enddo !isig


         endif                  ! k=2 configuration	
  		
c -------------

c third graph (k=1): order = open index, jal, jbe (LO)


         if (num_ret_current.eq.1) then ! jexch,jal,jbe

	    do isig = -1,1,2			
				
            call bra2c(psi(1,isig,nq2),.true.,p(0,nq2),isig,qak(0,2),
     #           jal(0),psiak(1,isig,2,nq2),fqak(0,2,nq2))
         
            call bra2c(psi(1,isig,nq2),.true.,p(0,nq2),isig,qbk(0,2),
     #           jbe(0),psibk(1,isig,2,nq2),fqbk(0,2,nq2))
     
            call bra2c(psi(1,isig,nq2),.true.,p(0,nq2),isig,qbk(0,1),
     #           jbe(0),psibk(1,isig,1,nq2),fqbk(0,1,nq2))
    
            call bra2c(psi(1,isig,nq2),.true.,p(0,nq2),isig,qak(0,2),
     #           cak(0,2),psipak(1,isig,2,nq2),fpak(0,2,nq2))     
     
	    enddo !isig	
		
		
            call curr6(1,psiak(1,-1,2,nq2),fqak(0,2,nq2),
     1           psi(1,-1,nq1),p(0,nq1),j0ak(0,-1,2)   )
		
            call curr6(1,psibk(1,-1,2,nq2),fqbk(0,2,nq2),
     1           psi(1,-1,nq1),p(0,nq1),j0bk(0,-1,2)   )

            call curr6(1,psibk(1,-1,1,nq2),fqbk(0,1,nq2),
     1           psi(1,-1,nq1),p(0,nq1),j0bk(0,-1,1)   )

            call curr6(1,psipak(1,-1,2,nq2),fpak(0,2,nq2),
     1           psi(1,-1,nq1),p(0,nq1),j0aa(0,-1,2)   )


c third graph (k=1): NLO for order = open index, jal, jbe 
c 
c
            call boxlinec(1,psi(1,-1,nq1),psi(1,-1,nq2),p(0,nq1),p(0,nq2),
     #           .true.,jak(0,2),1,j0ak(0,-1,2),jvirt1(0,-1))

            call boxlinec(1,psi(1,-1,nq1),psi(1,-1,nq2),p(0,nq1),p(0,nq2),
     #           .true.,jbk(0,2),1,j0bk(0,-1,2),jvirt2(0,-1))
     
            call boxlinec(1,psi(1,-1,nq1),psi(1,-1,nq2),p(0,nq1),p(0,nq2),
     #	   	.true.,jbk(0,1),1,j0bk(0,-1,1),jvirt3(0,-1))
    
            call boxlinec(1,psi(1,-1,nq1),psi(1,-1,nq2),p(0,nq1),p(0,nq2),
     #	   	.true.,cak(0,2),1,j0aa(0,-1,2),jvirt4(0,-1))
    
		 
    ! emit wp, then wm from upper line: 
    
           do isig = -1,1,2
            do mu = 0,5
               jvabd(mu,isig,1)= -xbe*jvirt1(mu,isig) 
     #                    + xal*(jvirt2(mu,isig) - jvirt3(mu,isig))
     #                    - xbe*xal*jvirt4(mu,isig)    
            enddo               !mu
  	   enddo !isig	
		
         endif                  ! k=1 configuration	


c ---------------------------------------------------------------------

	do isig = -1,1,2
         do mu = 0,3
            jret(mu,isig) = -jvabd(mu,isig,num_ret_current) 
  		! ... add (-1) here to have same sign as for "genuine" pentagon 
         enddo
	enddo !isig 
		
      endif                     ! nlo=5	
	
      return
      end

c ---------------------------------------------------------------------
c ---------------------------------------------------------------------


      subroutine  J_virtual_tri_box_pent2(psi,p,nq1,nq2,isig,sign, musq,
c     input currents (any ORDERED couple among [J1,J2,J3]) in DHELAS format
c     They can be  jmu1,jmu2  or  jmu1,jmu3   or   jmu2,jmu3
     #     jal,jbe,
c     uncontracted external current. It can be 1, 2 or 3
     #     num_ret_current,
c     the corresponding Born current
c           Jborn,
c     returned current
     #     Jret)

      implicit none
      complex * 16 psi(2,-1:1,4)
      real * 8 p(0:3,4),q(0:3),l(0:3),t(0:3)
      real * 8 musq
      complex * 16 jal(0:5),jbe(0:5)
      integer num_ret_current, nq1,nq2,nq3,nq4, isig,sign(4)
      complex * 16 jret(0:3)
      logical lpr
      common /e0print/lpr

      integer mu
c map input momenta onto incoming standard and make use of 
c DHELAS information in currents
      if (nq1.eq.1 .and. nq2.eq.2) then
         nq3 = 3
         nq4 = 4
      elseif (nq1.eq.3 .and. nq2.eq.4) then
         nq3 = 1
         nq4 = 2
      else
         print*," check nq1 etc. in J_virtual_tri_box_pent2 call "
         print*,nq1,nq2,nq3,nq4
         stop
      endif
         
      if (num_ret_current.eq.1) then     ! p3=q+p4 for q incoming on 12 line
         do mu = 0,3
            q(mu) = p(mu,nq3)-p(mu,nq4)
         enddo
         l(0) = -dreal(jal(4))
         l(1) = -dreal(jal(5))
         l(2) = -dimag(jal(5))
         l(3) = -dimag(jal(4))
         t(0) = -dreal(jbe(4))
         t(1) = -dreal(jbe(5))
         t(2) = -dimag(jbe(5))
         t(3) = -dimag(jbe(4))
      elseif (num_ret_current.eq.2) then
         do mu = 0,3
            l(mu) = p(mu,nq3)-p(mu,nq4)
         enddo
         q(0) = -dreal(jal(4))
         q(1) = -dreal(jal(5))
         q(2) = -dimag(jal(5))
         q(3) = -dimag(jal(4))
         t(0) = -dreal(jbe(4))
         t(1) = -dreal(jbe(5))
         t(2) = -dimag(jbe(5))
         t(3) = -dimag(jbe(4))
      elseif (num_ret_current.eq.3) then
         do mu = 0,3
            t(mu) = p(mu,nq3)-p(mu,nq4)
         enddo
         q(0) = -dreal(jal(4))
         q(1) = -dreal(jal(5))
         q(2) = -dimag(jal(5))
         q(3) = -dimag(jal(4))
         l(0) = -dreal(jbe(4))
         l(1) = -dreal(jbe(5))
         l(2) = -dimag(jbe(5))
         l(3) = -dimag(jbe(4))
      endif

      do mu = 0,3
         jret(mu) = 0
      enddo
      call  J_virtual_tri_box_pent (
c     bra(p"), ket(p) and helicity 
     #     psi(1,-1,nq2),psi(1,-1,nq1),isig,
c     momenta and squared renormalization scale
     #     p(0,nq1),q,l,t,p(0,nq2),musq,
c     input currents (any ORDERED couple among [J1,J2,J3])
c     They can be  jmu1,jmu2  or  jmu1,jmu3   or   jmu2,jmu3
     #     jal,jbe,
c     uncontracted external current. It can be 1, 2 or 3
     #     num_ret_current,
c     returned current
     #     Jret)

c check momentum conservation and subtract Born contribution
      do mu = 0,3
         t(mu) = q(mu)+l(mu)+t(mu)+p(mu,nq1)-p(mu,nq2) 
         if (abs(t(mu)).gt.1d-8) then
            print*," no momentum conservation in j_virtual_tri_box_pent"
            print*," mom diff(",mu,") = ",t(mu)
         endif
      enddo
      return
      end

