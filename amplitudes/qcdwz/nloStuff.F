c ************************************************************************************
c Author: Matthias Kerner
c Date: 2012/13
c ************************************************************************************

c #define QUAD_ENHANCED
c #define NINH_COMPARE
      subroutine precalcQCDWZjj(p,lokt,lephel,nbos)
      implicit none
#include "brakets.inc"
#include "VBFNLO/utilities/global.inc"
      real*8 p(0:3,max_p,max_kin)
      logical lokt(16)
      integer lephel,lepkin,nbos

      real*8 mjj2,dotrr,dotrrQ
      external mjj2,dotrr,dotrrQ
      real*8 pdia(0:4,5),mpdia(0:4,5)
      integer fsignp(4),fsignm(4),id(4),h,i,j,idb(4),kin,copy(5,16),mu,pol1,pol2,qtype(4)
      data fsignp / 1,1,1,1/, fsignm / -1,-1,-1,-1 /, id /1,3,2,4/,idb/4,2,3,1/
      complex*16 psiTemp(2,-1:1)
      real*8 momTemp(0:4)

      complex*16 epsW(0:3,7),epsWV(0:3,-1:1,7) 
      complex*16 epsV(0:3,-1:1,3:4,-1:1,7) ! ( , mu-hel, u/d , quark-hel), coupl. const to quarks included
      complex*16 epsA(0:3,-1:1),epsZ(0:3,-1:1)
      real*8 momV(0:4,7),momW(0:4,7),momWV(0:4,7)
      common /qcdwznew/ epsW, epsWV, epsV,epsA,epsZ,momV,momW,momWV

c dont' recalculated spinors, that are the same as in kin 1
      data copy /0,0,0,0,0, 0,2,0,0,0, 0,2,0,0,0, 0,2,0,0,0, ! RE kin, ii kin
     $                      1,0,0,0,0, 1,0,0,0,0, 1,0,0,0,0,
     $                      1,2,0,0,0, 1,2,0,0,0, 1,2,0,0,0, ! ff kin
     $                      0,2,0,5,0, 0,2,0,4,0, 0,2,0,3,0, ! if, fi
     $                      1,0,0,5,0, 1,0,0,4,0, 1,0,0,3,0/
      
c generate psi 

      if(nbos.eq.32.or.nbos.eq.3.or.nbos.eq.31) then
	qtype(1) = 3
	qtype(2) = 4
      else
	qtype(1) = 4
	qtype(2) = 3
      endif

      do kin=1,16
      if(lokt(kin)) then 

c use pdia for psi, -pdia for psibar
c -pdia for bosons !!!!!
      do mu=0,3
        pdia(mu,1) = p(mu,1,kin) ! psi in initial state has sign 1
        pdia(mu,2) = p(mu,2,kin)
        pdia(mu,3) =-p(mu,3,kin) ! psi in final state has sign -1
        pdia(mu,4) =-p(mu,4,kin) 
        pdia(mu,5) =-p(mu,5,kin) 
        do i=1,5
          mpdia(mu,i) = -pdia(mu,i)
        enddo
      enddo
      do i=1,5
        pdia(4,i) = 0d0
      enddo
      
      lepkin=kin
      if(kin.gt.7) lepkin = 1

c calc spinors
c kets are psi(:,:, >0 , :)
c bras are psi(:,:, <0 , :)
        call psi0mAll(5,p(0,1,kin),psi(1,-1,-5,kin))
      
      do i=1,5
        call polvec(p(0,i,kin),1,epsG(0,i,1,kin))
        call polvec(p(0,i,kin),2,epsG(0,i,2,kin))
      enddo

      do i=1,4; do j=i+1,5
      do pol1=1,2; do pol2=1,2
      do mu = 0,3
c gluon copplings
c add a minus sign because of pdia convention 
        epsGG(mu,i,j,pol1,pol2,kin)= -(2d0*dotrr(pdia(0,j),epsG(0,i,pol1,kin))*epsG(mu,j,pol2,kin)
     &		    -2d0*dotrr(pdia(0,i),epsG(0,j,pol2,kin))*epsG(mu,i,pol1,kin)
     &        +dotrr(epsG(0,i,pol1,kin),epsG(0,j,pol2,kin))*(pdia(mu,i)-pdia(mu,j)))/mjj2(pdia(0,i),pdia(0,j))
      enddo
      enddo;enddo
      enddo;enddo 

c psiG
      do i=1,5
      if(.not.(kin.ne.1.and.i.eq.5)) then
        do j=1,5
        if(i.ne.j) then
          call bra2r(psi(1,-1,-i,kin),.true.,mpdia(0,i),-1,mpdia(0,j),epsG(0,j,1,kin),psiG(1,-1,-i,j,1,kin),momPsiG(0,-i,j,kin))
          call bra2r(psi(1, 1,-i,kin),.true.,mpdia(0,i), 1,mpdia(0,j),epsG(0,j,1,kin),psiG(1, 1,-i,j,1,kin),momPsiG(0,-i,j,kin))
          call bra2r(psi(1,-1,-i,kin),.true.,mpdia(0,i),-1,mpdia(0,j),epsG(0,j,2,kin),psiG(1,-1,-i,j,2,kin),momPsiG(0,-i,j,kin))
          call bra2r(psi(1, 1,-i,kin),.true.,mpdia(0,i), 1,mpdia(0,j),epsG(0,j,2,kin),psiG(1, 1,-i,j,2,kin),momPsiG(0,-i,j,kin))

          call ket2r(psi(1,-1, i,kin),.true.,pdia(0,i),-1,mpdia(0,j),epsG(0,j,1,kin),psiG(1,-1, i,j,1,kin),momPsiG(0,i,j,kin))
          call ket2r(psi(1, 1, i,kin),.true.,pdia(0,i), 1,mpdia(0,j),epsG(0,j,1,kin),psiG(1, 1, i,j,1,kin),momPsiG(0,i,j,kin))
          call ket2r(psi(1,-1, i,kin),.true.,pdia(0,i),-1,mpdia(0,j),epsG(0,j,2,kin),psiG(1,-1, i,j,2,kin),momPsiG(0,i,j,kin))
          call ket2r(psi(1, 1, i,kin),.true.,pdia(0,i), 1,mpdia(0,j),epsG(0,j,2,kin),psiG(1, 1, i,j,2,kin),momPsiG(0,i,j,kin))
!           print*,psiG(:,-1,-i,j,:,kin)
        endif
        enddo

c psiV, psiW
        call bra2c(psi(1,-1,-i,kin),.true.,mpdia(0,i),-1,momW(0,lepkin),epsW(0,lepkin),
     $             psiW(1,-i,kin),momPsiV(0,-i,1,kin))
        call ket2c(psi(1,-1, i,kin),.true., pdia(0,i),-1,momW(0,lepkin),epsW(0,lepkin),
     $             psiW(1, i,kin),momPsiV(0,i,1,kin))

        do j=3,4
        do h=-1,1,2
          call bra2c(psi(1,h,-i,kin),.true.,mpdia(0,i),h,momV(0,lepkin),epsV(0,lephel,j,h,lepkin),
     $               psiV(1,h,-i,j,kin),momPsiV(0,-i,2,kin))
          call ket2c(psi(1,h, i,kin),.true.,pdia(0,i),h,momV(0,lepkin),epsV(0,lephel,j,h,lepkin),
     $               psiV(1,h, i,j,kin),momPsiV(0,i,2,kin))
        enddo
        enddo
c psiVV
        call bra2c(psi(1,-1,-i,kin),.true.,mpdia(0,i),-1,momWV(0,lepkin),epsWV(0,lephel,lepkin),
     $             psiWV(1,-i,kin),momPsiV(0,-i,3,kin))
        call ket2c(psi(1,-1, i,kin),.true.,pdia(0,i),-1,momWV(0,lepkin),epsWV(0,lephel,lepkin),psiWV(1, i,kin),momPsiV(0,i,3,kin))

	call ket2c(psiV(1,-1,i,qtype(1),kin),.false.,momPsiV(0,i,2,kin),-1,momW(0,lepkin), epsW(0,lepkin),PsiTemp(1,-1),momTemp)
	call ket2c(psiW(1,i,kin),.false.,momPsiV(0,i,1,kin),-1,momV(0,lepkin), epsV(0,lephel,qtype(2),-1,lepkin),
     $             PsiTemp(1,1),momTemp)
	psiWV(1,i,kin) = psitemp(1,-1) + psitemp(1,1) - psiWV(1,i,kin)
	psiWV(2,i,kin) = psitemp(2,-1) + psitemp(2,1) - psiWV(2,i,kin)

	call bra2c(psiV(1,-1,-i,qtype(2),kin),.false.,momPsiV(0,-i,2,kin),-1,momW(0,lepkin), epsW(0,lepkin),
     $             PsiTemp(1,-1),momTemp)
	call bra2c(psiW(1,-i,kin),.false.,momPsiV(0,-i,1,kin),-1,momV(0,lepkin), epsV(0,lephel,qtype(1),-1,lepkin),
     $             PsiTemp(1,1),momTemp)
	psiWV(1,-i,kin) = psitemp(1,-1) + psitemp(1,1) - psiWV(1,-i,kin)
	psiWV(2,-i,kin) = psitemp(2,-1) + psitemp(2,1) - psiWV(2,-i,kin)

c psiWG
        do j=1,5
        if(i.ne.j) then
	  do pol1 = 1,2
	    call ket2c(psiG(1,-1,i,j,pol1,kin),.false.,momPsiG(0,i,j,kin),-1,momW(0,lepkin),epsW(0,lepkin),PsiTemp(1,-1),momtemp)
	    call ket2r(psiW(1,i,kin),.false.,momPsiV(0,i,1,kin),-1,mpdia(0,j),epsG(0,j,pol1,kin),
     $                 PsiTemp(1,1),momPsiWG(0,i,j,kin))
	    PsiWG(1,i,j,pol1,kin) = PsiTemp(1,-1)+PsiTemp(1,1)
            PsiWG(2,i,j,pol1,kin) = PsiTemp(2,-1)+PsiTemp(2,1)

	    call bra2c(psiG(1,-1,-i,j,pol1,kin),.false.,momPsiG(0,-i,j,kin),-1,momW(0,lepkin),epsW(0,lepkin),
     $                 PsiTemp(1,-1),momtemp)
	    call bra2r(psiW(1,-i,kin),.false.,momPsiV(0,-i,1,kin),-1,mpdia(0,j),epsG(0,j,pol1,kin),
     $                 PsiTemp(1,1),momPsiWG(0,-i,j,kin))
	    PsiWG(1,-i,j,pol1,kin) = PsiTemp(1,-1)+PsiTemp(1,1)
            PsiWG(2,-i,j,pol1,kin) = PsiTemp(2,-1)+PsiTemp(2,1)
	  enddo
	endif
        enddo
   
! 
      endif
      enddo
      
      endif ! lokt
      enddo

      end

      subroutine precalcQCDZZjj(p,lokt,lephel)
      implicit none
#include "cacheZZ.inc"
#include "VBFNLO/utilities/global.inc"
      real*8 p(0:3,max_p,max_kin)
      logical lokt(16)
      integer lephel,lepkin

      real*8 mjj2,dotrr,dotrrQ
      external mjj2,dotrr,dotrrQ
      real*8 pdia(0:4,5),mpdia(0:4,5)
      integer fsignp(4),fsignm(4),id(4),i,j,idb(4),kin,copy(5,16),mu,pol1,qh,qt,iV,iV2,nfmax
      data fsignp / 1,1,1,1/, fsignm / -1,-1,-1,-1 /, id /1,3,2,4/,idb/4,2,3,1/
      complex*16 psiTemp(2,-1:1)
      real*8 momTemp(0:4)


c dont' recalculated spinors, that are the same as in kin 1
      data copy /0,0,0,0,0, 0,2,0,0,0, 0,2,0,0,0, 0,2,0,0,0, ! RE kin, ii kin
     $                      1,0,0,0,0, 1,0,0,0,0, 1,0,0,0,0,
     $                      1,2,0,0,0, 1,2,0,0,0, 1,2,0,0,0, ! ff kin
     $                      0,2,0,5,0, 0,2,0,4,0, 0,2,0,3,0, ! if, fi
     $                      1,0,0,5,0, 1,0,0,4,0, 1,0,0,3,0/

c V currents
c V1/V2/V12,
c    for qhel=0:
c       AZ=1,2 : A and Z current wo coupling
c    qhel = -1:1:
c       AZ=3,4 : V coupling to u/d quark with qhel inclouding couplings
      complex*16 epsV,w4l(0:4) ! (mu, qhel, 12=AZ 34=including coupl to ud, V1 V2 V12, hel->0=rnd hel, kin)
      real*8 momV
      common /QCDZZcurrents/ epsV(0:3,-1:1,4,3,0:4,7),w4l,momV(0:4,3,7)
      logical lAZ(0:2,3) ! V.or.A A V,    V1 V2 V12
      common /QCDZZcacheL/ lAZ
      
      
c generate psi 
      do kin=1,16
      if(lokt(kin)) then 
        if (kin.eq.1) then
          nfmax = 5
        else
          nfmax = 4
        endif
c use pdia for psi, -pdia for psibar
c -pdia for bosons !!!!!
      do mu=0,3
        pdia(mu,1) = p(mu,1,kin) ! psi in initial state has sign 1
        pdia(mu,2) = p(mu,2,kin)
        pdia(mu,3) =-p(mu,3,kin) ! psi in final state has sign -1
        pdia(mu,4) =-p(mu,4,kin) 
        pdia(mu,5) =-p(mu,5,kin) 
        do i=1,5
          mpdia(mu,i) = -pdia(mu,i)
        enddo
      enddo
      do i=1,5
        pdia(4,i) = 0d0
      enddo
      
      lepkin=kin
      if(kin.gt.7) lepkin = 1

c calc spinors
c kets are psi(:,:, >0 , :)
c bras are psi(:,:, <0 , :)
        call psi0mAll(nfmax,p(0,1,kin),psi(1,-1,-nfmax,kin))
      
      do i=1,nfmax
        call polvec(p(0,i,kin),1,epsG(0,i,1,kin))
        call polvec(p(0,i,kin),2,epsG(0,i,2,kin))
      enddo
! 
!       do i=1,4; do j=i+1,5
!       do pol1=1,2; do pol2=1,2
!       do mu = 0,3
! c gluon copplings
! c add a minus sign because of pdia convention 
!         epsGG(mu,i,j,pol1,pol2,kin)= -(2d0*dotrrQ(pdia(0,j),epsG(0,i,pol1,kin))*epsG(mu,j,pol2,kin)
!      &              -2d0*dotrrQ(pdia(0,i),epsG(0,j,pol2,kin))*epsG(mu,i,pol1,kin)
!      &        +dotrrQ(epsG(0,i,pol1,kin),epsG(0,j,pol2,kin))*(pdia(mu,i)-pdia(mu,j)))/mjj2Q(pdia(0,i),pdia(0,j))
!       enddo
!       enddo;enddo
!       enddo;enddo 

c psiG
      do i=1,nfmax
        do j=1,nfmax
        if(i.ne.j) then
          call bra2r(psi(1,-1,-i,kin),.true.,mpdia(0,i),-1,mpdia(0,j),epsG(0,j,1,kin),psiG(1,-1,-i,j,1,kin),momPsiG(0,-i,j,kin))
          call bra2r(psi(1, 1,-i,kin),.true.,mpdia(0,i), 1,mpdia(0,j),epsG(0,j,1,kin),psiG(1, 1,-i,j,1,kin),momPsiG(0,-i,j,kin))
          call bra2r(psi(1,-1,-i,kin),.true.,mpdia(0,i),-1,mpdia(0,j),epsG(0,j,2,kin),psiG(1,-1,-i,j,2,kin),momPsiG(0,-i,j,kin))
          call bra2r(psi(1, 1,-i,kin),.true.,mpdia(0,i), 1,mpdia(0,j),epsG(0,j,2,kin),psiG(1, 1,-i,j,2,kin),momPsiG(0,-i,j,kin))

          call ket2r(psi(1,-1, i,kin),.true.,pdia(0,i),-1,mpdia(0,j),epsG(0,j,1,kin),psiG(1,-1, i,j,1,kin),momPsiG(0,i,j,kin))
          call ket2r(psi(1, 1, i,kin),.true.,pdia(0,i), 1,mpdia(0,j),epsG(0,j,1,kin),psiG(1, 1, i,j,1,kin),momPsiG(0,i,j,kin))
          call ket2r(psi(1,-1, i,kin),.true.,pdia(0,i),-1,mpdia(0,j),epsG(0,j,2,kin),psiG(1,-1, i,j,2,kin),momPsiG(0,i,j,kin))
          call ket2r(psi(1, 1, i,kin),.true.,pdia(0,i), 1,mpdia(0,j),epsG(0,j,2,kin),psiG(1, 1, i,j,2,kin),momPsiG(0,i,j,kin))
!           print*,psiG(:,-1,-i,j,:,kin)
        endif
        enddo

c psiV(V)
        do qh=-1,1,2
        do qt=3,4 
          do iV=3,1,-1
          if(lAZ(0,iV)) then
            call bra2c(psi(1,qh,-i,kin),.true.,mpdia(0,i),qh,momV(0,iV,lepkin),epsV(0,qh,qt,iV,lephel,lepkin),
     $             psiV(1,qh,qt,-i,iV,kin),momPsiV(0,-i,iV,kin))
            call ket2c(psi(1,qh, i,kin),.true., pdia(0,i),qh,momV(0,iV,lepkin),epsV(0,qh,qt,iV,lephel,lepkin),
     $             psiV(1,qh,qt, i,iV,kin),momPsiV(0, i,iV,kin))
        
            if(iV.ne.3) then
             iV2=3-iV
             if(lAZ(0,iV2)) then
              call bra2c(psiV(1,qh,qt,-i,iV,kin),.false.,momPsiV(0,-i,iV,kin),qh,momV(0,iV2,lepkin),epsV(0,qh,qt,iV2,lephel,lepkin),
     $             psitemp(1,1),momTemp)
                  psiV(:,qh,qt,-i,3,kin) = psiV(:,qh,qt,-i,3,kin) + psitemp(:,1)
              call ket2c(psiV(1,qh,qt, i,iV,kin),.false.,momPsiV(0, i,iV,kin),qh,momV(0,iV2,lepkin),epsV(0,qh,qt,iV2,lephel,lepkin),
     $                  psitemp(1,1),momTemp)
                    psiV(:,qh,qt, i,3,kin) = psiV(:,qh,qt, i,3,kin) + psitemp(:,1)
              endif
             endif
          endif
          enddo

! c psiVG 
          if(kin.eq.1) then ! only needed for N+1
          do iV=3,1,-1
          if (lAZ(0,iV)) then
          do j=1,5
          if(i.ne.j) then
            do pol1 = 1,2
              call ket2c(psiG(1,qh,i,j,pol1,kin),.false.,momPsiG(0,i,j,kin),qh,momV(0,iV,lepkin),epsV(0,qh,qt,iV,lephel,lepkin),
     $                   PsiTemp(1,-1),momtemp)
              call ket2r(psiV(1,qh,qt, i,iV,kin),.false.,momPsiV(0,i,iV,kin),qh,mpdia(0,j),epsG(0,j,pol1,kin),
     $                 PsiTemp(1,1),momPsiVG(0,i,iV,j,kin))
              PsiVG(1,qh,qt,i,iV,j,pol1) = PsiTemp(1,-1)+PsiTemp(1,1)
              PsiVG(2,qh,qt,i,iV,j,pol1) = PsiTemp(2,-1)+PsiTemp(2,1)

              call bra2c(psiG(1,qh,-i,j,pol1,kin),.false.,momPsiG(0,-i,j,kin),qh,momV(0,iV,lepkin),epsV(0,qh,qt,iV,lephel,lepkin),
     $                 PsiTemp(1,-1),momtemp)
              call bra2r(psiV(1,qh,qt,-i,iV,kin),.false.,momPsiV(0,-i,iV,kin),qh,mpdia(0,j),epsG(0,j,pol1,kin),
     $                 PsiTemp(1,1),momPsiVG(0,-i,iV,j,kin))
              PsiVG(1,qh,qt,-i,iV,j,pol1) = PsiTemp(1,-1)+PsiTemp(1,1)
              PsiVG(2,qh,qt,-i,iV,j,pol1) = PsiTemp(2,-1)+PsiTemp(2,1)

              if(iV.ne.3) then
                iV2=3-iV
                call ket2c(psiVG(1,qh,qt,i,iV,j,pol1),.false.,momPsiVG(0,i,iV,j,kin),qh,momV(0,iV2,lepkin),
     $                     epsV(0,qh,qt,iV2,lephel,lepkin),PsiTemp(1,-1),momtemp)
                PsiVG(1,qh,qt,i,3,j,pol1) = PsiVG(1,qh,qt,i,3,j,pol1) + PsiTemp(1,-1)
                PsiVG(2,qh,qt,i,3,j,pol1) = PsiVG(2,qh,qt,i,3,j,pol1) + PsiTemp(2,-1)

                call bra2c(psiVG(1,qh,qt,-i,iV,j,pol1),.false.,momPsiVG(0,-i,iV,j,kin),qh,momV(0,iV2,lepkin),
     $                     epsV(0,qh,qt,iV2,lephel,lepkin),PsiTemp(1,-1),momtemp)
                PsiVG(1,qh,qt,-i,3,j,pol1) = PsiVG(1,qh,qt,-i,3,j,pol1) + PsiTemp(1,-1)
                PsiVG(2,qh,qt,-i,3,j,pol1) = PsiVG(2,qh,qt,-i,3,j,pol1) + PsiTemp(2,-1)
              endif
                
              
            enddo
          endif
          enddo
          endif
          enddo
          endif
        
        enddo ! qt
        enddo ! qh
! 
      enddo
      
      endif ! lokt
      enddo

      end


c resturn all spinor psi(::,-NF:NF) where a negative sign of NF corresponds to psibar
      SUBROUTINE PSI0MAll( NF,PBAR,PSI )
      IMPLICIT NONE
      INTEGER  IF, NF
      DOUBLE COMPLEX  PSI(2,-1:1,-NF:NF)
      DOUBLE PRECISION  PBAR(0:3,NF),  PBARX, PBARY, PBARZ
      real*8 PABS,PAPZ
      DOUBLE PRECISION  NORMAL, NX, NY, NZ, EPS
      PARAMETER ( EPS=1.d-30 )
CC
      DO IF = 1,NF
        PABS  = PBAR(0,IF)
        PBARX = PBAR(1,IF)
        PBARY = PBAR(2,IF)
        PBARZ = PBAR(3,IF)
        IF ( PBARZ.GT.0.d0 ) THEN
           PAPZ = PABS + PBARZ
        ELSE
           PAPZ = (PBARX**2 + PBARY**2)/(PABS - PBARZ)
        ENDIF
c
c treat the pbar along negative z-axis case first
c
        IF ( PAPZ.LE.EPS*PABS ) THEN
           NORMAL = SQRT(2d0*PABS)
           if(if.ge.3) normal = -normal
           PSI(2,-1,IF) = 0.d0
           PSI(1,1,IF)  = 0.d0
           PSI(1,-1,IF) = -NORMAL 
           PSI(2,1,IF)  = NORMAL

           PSI(2,-1,-IF) = 0.d0
           PSI(1,1,-IF)  = 0.d0
           PSI(1,-1,-IF) = NORMAL 
           PSI(2,1,-IF)  = -NORMAL
        ELSE
c
c and now the general case
c
           NORMAL = 1d0/SQRT(PAPZ)
           if(if.ge.3) normal = -normal
           NX = NORMAL*PBARX
           NY = NORMAL*PBARY
           NZ = NORMAL*PAPZ
           PSI(1,-1,IF) = DCMPLX(-NX,NY)
           PSI(2,-1,IF) = NZ
           PSI(1,1,IF)  = NZ
           PSI(2,1,IF)  = DCMPLX(NX,NY)
           PSI(1,-1,-IF) = DCMPLX(NX,NY)
           PSI(2,-1,-IF) = -NZ
           PSI(1,1,-IF)  = -NZ
           PSI(2,1,-IF)  = DCMPLX(-NX,NY)
        ENDIF
      ENDDO
ccc
      END


c amp(0) pol != pol', amp(i): pol(i)
c buff==nbos
      subroutine CScolAmpGNew(p,kin,idB,fsignB,emitB,spectB,nbos,res)
c get <TiTj> for 2q2g processes, amp = (Born^2, T1T2,T1T3)
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "common.inc"
      real*8 p(0:3,max_p,max_kin),res(-1:2,3:4), amp(0:2,-1:2)
      integer emitB,spectB,mini,maxi, kin, idB(4),fsignB(4),nbos,i,j,k,nres
      real*8 CF,CA
      parameter (CF=4d0/3d0,CA=3d0)
      logical ldebug,usePreCalc
      
      real*8 borncache(0:2,0:2,0:2,3:4,1:4,1:4,16)
      save borncache
      integer helswitch,q
      common/helswitchTemp/ helswitch

      parameter (ldebug = .false.)
      usePreCalc=.false.
      goto 23
      entry CScolAmpGNewUsePrecalc(p,kin,idB,fsignB,emitB,spectB,nbos,res)
      UsePreCalc=.true.

      goto 23
      entry CScolAmpGNewReset()

      do i=1,4; do j=1,4; do k=1,16
        borncache(0,0,0,3,i,j,k) = -1d60
      enddo; enddo; enddo

      return

 23   continue

      if(idB(4).lt.idB(2)) then
        i=idB(4)
        idB(4)=idB(2)
        idB(2)=i
        if(mod(emitB,2).eq.0) emitB=6-emitB
        if(mod(spectB,2).eq.0) spectB=6-spectB
      endif

      if(ldebug) then
	if(idB(4).lt.idB(2)) print*,'id(2) and id(4) should be ordered in PS for borncache'
	if(emitB.eq.spectB) print*,'calculate TiTi, maybe something''s wrong'
      endif
      nres=1
      if(nbos/10.eq.2.or.abs(nbos).eq.2) nres=2 ! ZZ and Z
c look for the amplitudes in borncache
      if(borncache(0,0,0,3,idB(1),idB(3),kin).lt.0d0) then
	if(nbos.eq.32.or.nbos.eq.42.or.nbos.eq.31.or.nbos.eq.41) then !WZ
#ifdef WITH_QCDVV
#ifdef NINH_COMPARE
! 	  lokt=.false.
! 	  call precalcQCDVVjj(p,lokt,-1,nbos)
! 	    ug2dgwzCSNew(p,kin,idB,fsignB,nbos,-1,borncache(0,0,0,idB(1),idB(3),kin))
! 	  call precalcQCDVVjj(p,lokt,-1,nbos)
! 	    ug2dgwzCSNew(p,kin,idB,fsignB,nbos,-1,amptemp)
! 	  borncache(:,:,:,idB(1),idB(3),kin)=borncache(:,:,:,idB(1),idB(3),kin)+amptemp
	  call ug2dgwzCSNew(p,kin,idB,fsignB,nbos,helswitch,.false.,borncache(0,0,0,3,idB(1),idB(3),kin))
#else
          call ug2dgwzCSNew(p,kin,idB,fsignB,nbos,0,usePreCalc,borncache(0,0,0,3,idB(1),idB(3),kin))
          call ug2dgwzCSNew(p,kin,idB,fsignB,nbos,0,usePreCalc,borncache(0,0,0,3,idB(1),idB(3),kin))
#endif
#endif
        elseif(nbos.eq.23.or.nbos.eq.22.or.nbos.eq.21.or.nbos.eq.20.or.abs(nbos).eq.2.or.nbos.eq.24) then
          nres=2
          call amp2q2gZZcs(p,kin,idB,fsignB,0,usePreCalc,borncache(0,0,0,3,idB(1),idB(3),kin))
	elseif(nbos.eq.3.or.nbos.eq.4) then ! W
 	  call ug2dgwCS(p,kin,idB,fsignB,nbos,borncache(0,0,0,3,idB(1),idB(3),kin))
        else
          print*,"undefined Born Matrix element",__FILE__,__LINE__
          stop
	endif
      endif
        
      do q=3,2+nres
      do i=0,2
        if(emitB.eq.2) then ! sum over pol4
          do j=0,2
            amp(i,j) = borncache(i,j,1,q,idB(1),idB(3),kin) + borncache(i,j,2,q,idB(1),idB(3),kin) 
          enddo
          amp(i,-1) = amp(i,0) ! offdiag
          amp(i,0) = borncache(i,0,0,q,idB(1),idB(3),kin) ! summed
        else
          do j=0,2
            amp(i,j) = borncache(i,1,j,q,idB(1),idB(3),kin) + borncache(i,2,j,q,idB(1),idB(3),kin)
          enddo
          amp(i,-1) = amp(i,0)
          amp(i,0) = borncache(i,0,0,q,idB(1),idB(3),kin)
        endif
      enddo

c calculate TiTj
      if(emitB.lt.spectB) then
	mini=emitB
	maxi=spectB
      else
	mini=spectB
	maxi=emitB
      endif
      do i=-1,2
        if (mini.eq.0) then
          res(i,q) = amp(0,i) ! note that T_i^2 isn't included here
        else if(mini.eq.1) then
          if(maxi.eq.2) res(i,q)=amp(1,i)
          if(maxi.eq.3) res(i,q)=amp(2,i)
          if(maxi.eq.4) res(i,q)=-CF*amp(0,i)-amp(1,i)-amp(2,i)
        else if(mini.eq.2) then
          if(maxi.eq.3) res(i,q)=-CF*amp(0,i)-amp(1,i)-amp(2,i)
          if(maxi.eq.4) res(i,q)=(CF-CA)*amp(0,i)+amp(2,i)
        else 
            res(i,q)=amp(1,i)
        endif
      enddo
      enddo ! q
      
      end
! ! 


      subroutine CScolAmpQZZ(p,kin,idB,fsignB,emitB,spectB,nbos,res)
c get <TiTj> for 2q2g processes, amp = (Born^2, T1T2,T1T3)
      implicit none
#include "VBFNLO/utilities/global.inc"
      real*8 p(0:3,max_p,max_kin), amp(0:2),res(-1:2,3:4,3:4,0:2)
      integer h,i,j,ii,jj,mini,maxi, kin, idB(4),idtemp(4),signtemp(4),diffFl(0:2),fsignB(4)
      integer nbos, factorial(4),df,per,emitB,spectB
      integer a,b,q1,q2
      real*8 CF,CA
      parameter (CF=4d0/3d0,CA=3d0)
      logical ldebug
      parameter (ldebug = .false.)
      real*8 borncache(0:2,3:4,3:4,0:2,24,16)

      save borncache

      goto 25
      entry CScolAmpQZZReset()
      do i=1,24; do j=1,16
        borncache(0,3,3,0,i,j) = -1d60
      enddo; enddo
      return
      
      ! convert idb to factorial number

 25   difffl(0) = 0
      difffl(1) = 1
      difffl(2) = 2
      ii=emitB
      jj=spectB
      do i=1,4
        idtemp(i)=idb(i)
        signtemp(i)=fsignb(i)
      enddo
      select case (nbos)! get a unique permutation index, even if id quarks areinterchanged
      case(33,44) 
        if(idb(1).lt.idb(2)) then
          idtemp(1)=idb(1); signtemp(1)=fsignB(1)
          idtemp(2)=idb(2); signtemp(2)=fsignB(2)
        else
          idtemp(1)=idb(2); signtemp(1)=fsignB(2)
          idtemp(2)=idb(1); signtemp(2)=fsignB(1)
          difffl(1) = 2
          difffl(2) = 1
          if(ii.eq.1.or.ii.eq.2) ii=3-ii
          if(jj.eq.1.or.jj.eq.2) jj=3-jj
        endif
        if(idb(3).lt.idb(4)) then
          idtemp(3)=idb(3); signtemp(3)=fsignB(3)
          idtemp(4)=idb(4); signtemp(4)=fsignB(4)
        else
          idtemp(3)=idb(4); signtemp(3)=fsignB(4)
          idtemp(4)=idb(3); signtemp(4)=fsignB(3)
          diffFl(1)=3-diffFl(1)
          diffFl(2)=3-diffFl(2)
          if(ii.eq.3.or.ii.eq.4) ii=7-ii
          if(jj.eq.3.or.jj.eq.4) jj=7-jj
        endif
      end select
        
      do a=1,4
        factorial(a) = idtemp(a)-1
        do b=1,a-1
          if (idtemp(b) .lt. idtemp(a)) factorial(a)=factorial(a)-1
        enddo
      enddo
      per = factorial(1)*6+factorial(2)*2+factorial(3) +1

c look for the amplitudes in borncache
      if(borncache(0,3,3,0,per,kin).lt.0d0) then
        call amp4qZZcs(p,kin,idB,fsignB,0,
     1                  borncache(0,3,3,0,per,kin))
      endif

c calculate TiTj
      if(ii.lt.jj) then
        mini=ii
        maxi=jj
      else
        mini=jj
        maxi=ii
      endif

      do q1=3,4
      do q2=3,4
      do df=0,2
        res(-1,q1,q2,df) = 0d0
        res(1,q1,q2,df) = 0d0
        res(2,q1,q2,df) = 0d0
        do h=0,2
          amp(h) = borncache(h,q1,q2,difffl(df),per,kin)
        enddo
        if(mini.eq.0) then
          res(0,q1,q2,df) = amp(0) ! note that CF isn't included here
        else if(mini.eq.1) then
          if(maxi.eq.2) res(0,q1,q2,df)=amp(1)
          if(maxi.eq.3) res(0,q1,q2,df)=amp(2)
          if(maxi.eq.4) res(0,q1,q2,df)=-CF*amp(0)-amp(1)-amp(2)
        else if(mini.eq.2) then
          if(maxi.eq.3) res(0,q1,q2,df)=-CF*amp(0)-amp(1)-amp(2)
          if(maxi.eq.4) res(0,q1,q2,df)=amp(2)
        else
          res(0,q1,q2,df)=amp(1)
        endif
      enddo
      enddo
      enddo
!       print*,res(0,:,:,:)
      end

      subroutine CScolAmpQNew(p,kin,idB,fsignB,emitB,spectB,nbos,res)
c get <TiTj> for 2q2g processes, amp = (Born^2, T1T2,T1T3)
      implicit none
#include "VBFNLO/utilities/global.inc"
      real*8 p(0:3,max_p,max_kin), amp(0:2),res(-1:2,3:4,0:2)
      integer h,i,j,ii,jj,mini,maxi, kin, idB(4),idtemp(4),signtemp(4),diffFl(0:2),fsignB(4)
      integer nbos, factorial(4),q,df,per,emitB,spectB
      integer a,b
      real*8 CF,CA
      parameter (CF=4d0/3d0,CA=3d0)
      logical ldebug
      parameter (ldebug = .false.)
      real*8 borncache(0:2,3:4,0:2,24,16)
      real*8 tempamp(0:2,0:2)
      save borncache

      goto 24
      entry CScolAmpQNewReset()
      do i=1,24; do j=1,16
        borncache(0,3,0,i,j) = -1d60
      enddo; enddo
      return
      
      ! convert idb to factorial number

 24   difffl(0) = 0
      difffl(1) = 1
      difffl(2) = 2
      ii=emitB
      jj=spectB
      do i=1,4
        idtemp(i)=idb(i)
        signtemp(i)=fsignb(i)
      enddo
      select case (nbos)! get a unique permutation index, even if id quarks areinterchanged
      case(33,44) 
        if(idb(1).lt.idb(2)) then
          idtemp(1)=idb(1); signtemp(1)=fsignB(1)
          idtemp(2)=idb(2); signtemp(2)=fsignB(2)
        else
          idtemp(1)=idb(2); signtemp(1)=fsignB(2)
          idtemp(2)=idb(1); signtemp(2)=fsignB(1)
          difffl(1) = 2
          difffl(2) = 1
          if(ii.eq.1.or.ii.eq.2) ii=3-ii
          if(jj.eq.1.or.jj.eq.2) jj=3-jj
        endif
        if(idb(3).lt.idb(4)) then
          idtemp(3)=idb(3); signtemp(3)=fsignB(3)
          idtemp(4)=idb(4); signtemp(4)=fsignB(4)
        else
          idtemp(3)=idb(4); signtemp(3)=fsignB(4)
          idtemp(4)=idb(3); signtemp(4)=fsignB(3)
          diffFl(1)=3-diffFl(1)
          diffFl(2)=3-diffFl(2)
          if(ii.eq.3.or.ii.eq.4) ii=7-ii
          if(jj.eq.3.or.jj.eq.4) jj=7-jj
        endif
      end select
        
      do a=1,4
        factorial(a) = idtemp(a)-1
        do b=1,a-1
          if (idtemp(b) .lt. idtemp(a)) factorial(a)=factorial(a)-1
        enddo
      enddo
      per = factorial(1)*6+factorial(2)*2+factorial(3) +1

c look for the amplitudes in borncache
      if(borncache(0,3,0,per,kin).lt.0d0) then
        if(nbos.eq.32.or.nbos.eq.42.or.nbos.eq.31.or.nbos.eq.41) then
#ifdef WITH_QCDVV
          call amp4qWZcs(p,kin,idB,fsignB,nbos,0,
     1                  borncache(0,3,0,per,kin))
        elseif (nbos.eq.33.or.nbos.eq.44) then
          call amp4qWWcs(p,kin,idtemp,signtemp,tempamp)
          do i=0,2; do j=0,2
            borncache(i,3,j,per,kin) = tempamp(i,j)
          enddo; enddo
#endif
        elseif( nbos.eq.3.or.nbos.eq.4) then
          call amp4qWcs(p,kin,idB,fsignB,nbos,0,
     1                  borncache(0,3,0,per,kin))
        else
          print*,'bad value for nbos in',__FILE__,__LINE__
          stop
        endif
      endif

c calculate TiTj
      if(ii.lt.jj) then
        mini=ii
        maxi=jj
      else
        mini=jj
        maxi=ii
      endif

      do q=3,4
      do df=0,2
        res(-1,q,df) = 0d0
        res(1,q,df) = 0d0
        res(2,q,df) = 0d0
        do h=0,2
          amp(h) = borncache(h,q,difffl(df),per,kin)
        enddo
        if(mini.eq.0) then
          res(0,q,df) = amp(0) ! note that CF isn't included here
        else if(mini.eq.1) then
          if(maxi.eq.2) res(0,q,df)=amp(1)
          if(maxi.eq.3) res(0,q,df)=amp(2)
          if(maxi.eq.4) res(0,q,df)=-CF*amp(0)-amp(1)-amp(2)
        else if(mini.eq.2) then
          if(maxi.eq.3) res(0,q,df)=-CF*amp(0)-amp(1)-amp(2)
          if(maxi.eq.4) res(0,q,df)=amp(2)
        else
          res(0,q,df)=amp(1)
        endif
      enddo
      enddo
      end


      real*8 function QCDVVjjFinCollBornKin(p,xi,id,fsign,iflav,difffamily,nbos,pdf)
c     calculate finite collinear terms for processes pp > VV jj
c     this routine should be called with p(1:4): quarks, p(5): gluon 
c     initial state splitting g -> q qb, g ->gg will be calulated here as well
      implicit none
#include "common.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"
      real*8 p(0:3,max_p,max_kin),xi(2,2)
      integer id(4),fsign(4),iflav(4),in,nbos, difffamily
      logical ldebug
      parameter (ldebug = .true.)

      double precision  alphas, xmt, alfa, xmz, xmw,s2w, xmh
      common /bkopin/   alphas, xmt, alfa, xmz, xmw, s2w, xmh

      integer beamsign(2),i, dia(4),di,qtype,hel9(9),hel9temp(9),q1,q2
      data hel9 / 0,0,0,0,0,0,0,0,0/
      data hel9temp / -1,-1,-1,-1,-1,-1,-1,-1,-1/
      real*8 x,omx,xa,xxa,omxxa,logx,logomx,logomxxa,T(4,4),logMU2oSij(1,4)
      real*8 pdfqxa,pdfqxxa, pdfgxa,pdfgxxa,pdfb, ampborn(-1:2,3:4,0:2),born,csamp(-1:2,3:4)
      real*8 CF, CA, dotrr, pdf(-6:6,2,2),gamma(3:4),ampbornZZ(-1:2,3:4,3:4,0:2)
      real*8 deltaCReal, theta,logalpha,logax,QCDVVjjFinColltemp,CScolAmpG
      parameter (deltaCReal=0d0)
      parameter (CF=4d0/3d0, CA=3d0)
      external di,CScolAmpgNew,theta,CScolAmpG
!       print*,'ivlav',iflav
!       print*,'id',id

c process specific stuff
      if(nbos.eq.33.or.nbos.eq.44) then
        qtype = 3
        q1=3
      else if (nbos.eq.32.or.nbos.eq.42.or.nbos.eq.31.or.nbos.eq.41.or.nbos.eq.3.or.nbos.eq.4) then
        q1=3
        qtype = abs(iflav(id(2)))
        if(qtype.ne.0) qtype = 3+mod(qtype,2) ! 3 or 4, or 0 for glu
      else if (nbos.eq.23.or.nbos.eq.22.or.nbos.eq.21.or.abs(nbos).eq.2.or.nbos.eq.20.or.nbos.eq.24) then
        qtype = abs(iflav(id(2))) ! only need information if q or g
          q1 = abs(iflav(id(3)))
          q1 = 3+mod(q1,2) ! 3 or 4, or 0 for glu
        if(qtype.ne.0) then
          q2 = abs(iflav(id(4)))
          q2 = 3+mod(q2,2) ! 3 or 4, or 0 for glu
        endif
      else ! don't forget hel factor 
        print*,'wrong input in fincoll'
        stop
      endif
      
      do i=1,4
        dia(i) = di(id,i)
      enddo

      T(3,3)=CF
      T(4,4)=CF 
      if(iflav(3).eq.0) T(3,3) = CA
      if(iflav(4).eq.0) T(4,4) = CA
      
      
      QCDVVjjFinCollBornKin = 0d0
      do in = 1,2
      
      T(1,1) = CF
      if(iflav(in).eq.0) T(1,1) = CA

c get color amp
      if(qtype.eq.0) then ! 2q2g amplitudes
        call CScolAmpGNew(p,1,id,fsign,dia(1),dia(2),nbos,csamp) ! use hel summed amplitudes
        T(1,2) = csamp(0,q1)
        call CScolAmpGNew(p,1,id,fsign,dia(in),dia(3),nbos,csamp)
        T(1,3) = csamp(0,q1)
        call CScolAmpGNew(p,1,id,fsign,dia(in),dia(4),nbos,csamp)
        T(1,4) = csamp(0,q1)
        call CScolAmpGNew(p,1,id,fsign,0,0,nbos,csamp)
        born   = csamp(0,q1)
      else ! 4q amplitudes
        select case (nbos)
        case(33,44,32,42,31,41,3,4)
          call CScolAmpqNew(p,1,id,fsign,dia(1),dia(2),nbos,ampborn)
          T(1,2) = ampborn(0,qtype,difffamily) ! difffamily is her same family
          call CScolAmpqNew(p,1,id,fsign,dia(in),dia(3),nbos,ampborn)
          T(1,3) = ampborn(0,qtype,difffamily)
          call CScolAmpqNew(p,1,id,fsign,dia(in),dia(4),nbos,ampborn)
          T(1,4) = ampborn(0,qtype,difffamily)
          call CScolAmpqNew(p,1,id,fsign,0,0,nbos,ampborn)
          born = ampborn(0,qtype,difffamily)
        case(23,22,21,20,24,2,-2)
!           print*,q1,q2,difffamily
!           print*
          call CScolAmpqZZ(p,1,id,fsign,dia(1),dia(2),nbos,ampbornZZ)
          T(1,2) = ampbornZZ(0,q1,q2,difffamily) ! difffamily is her same family
          call CScolAmpqZZ(p,1,id,fsign,dia(in),dia(3),nbos,ampbornZZ)
          T(1,3) = ampbornZZ(0,q1,q2,difffamily) ! difffamily is her same family
          call CScolAmpqZZ(p,1,id,fsign,dia(in),dia(4),nbos,ampbornZZ)
          T(1,4) = ampbornZZ(0,q1,q2,difffamily) ! difffamily is her same family
          call CScolAmpqZZ(p,1,id,fsign,0,0,nbos,ampbornZZ)
          born = ampbornZZ(0,q1,q2,difffamily) ! difffamily is her same family
! ! 
! !           call amp4qzzBorn(P,id,fsign,0,ampbornZZ(0:2,:,:,:))
! !           print*,born/ampbornZZ(0,q1,q2,difffamily)
! ! !           print*
! ! !           print*,born,T(1,2),T(1,3),T(1,4)
! ! !           pause
        case default
          print*,__FILE__, __LINE__
          stop
        end select
      endif

c      born = 1d0*4*3*8
c      T(1,2:4) = 1d0*4*3*8
c      if(iflav(in).eq.0) then
c      print*,born/4/3/8
c      print*,T(1,2:4)/4/3/8
c      endif

      x=xi(in,1)/xi(in,2)
      omx = 1d0-x
      xa = xi(in,2)
      xxa = xi(in,1)
      omxxa = 1d0 -xxa
      logx = log(x)
      logomx = log(omx)
      logomxxa = log(omxxa)
      logalpha = log(alphaDip)
      logax = log(min(1d0,alphaDip/omx))
      
      logmu2osij(1,2) = log(mufsq(in,1)/(2d0*dotrr(p(0,1,1),p(0,2,1))))
      logmu2osij(1,3) = log(mufsq(in,1)/(2d0*dotrr(p(0,in,1),p(0,3,1))))
      logmu2osij(1,4) = log(mufsq(in,1)/(2d0*dotrr(p(0,in,1),p(0,4,1))))

c      print*,'log',in,logmu2osij(1,2)
c      print*,'log',in,logmu2osij(1,3)
c      print*,'log',in,logmu2osij(1,4)
      if(xi(2,2).lt.xi(2,1).or.xi(2,2).gt.1d0.or.xi(2,1).lt.0d0) print *,xi
      if(xi(2,2).lt.xi(2,1).or.xi(2,2).gt.1d0.or.xi(2,1).lt.0d0) print *,xi
      if(xi(2,2).lt.xi(2,1).or.xi(2,2).gt.1d0.or.xi(2,1).lt.0d0) print *,xi
      
      beamsign(1) = sign1
      beamsign(2) = sign2
       
      gamma(3) = 2d0 ! gamma(q)
      gamma(4) = 2d0
      if(iflav(3).eq.0) gamma(3) = 5.5d0-nfl/3d0 ! gamma(g)
      if(iflav(4).eq.0) gamma(4) = 5.5d0-nfl/3d0
           
      pdfb = pdf(beamsign(3-in)*iflav(3-in),3-in,1)

      if(ldebug) then
        if(  abs((T(1,2)+T(1,3)+T(1,4))/born/T(1,1) + 1d0) .gt. 1d-5) then
            print*,"something's wrong in FinColl"
            print*,'iflav',iflav
        endif
      endif

      
c set pdf and calc fincoll
      if(iflav(in).eq.0) then ! born with g initial state
        pdfgxa = pdf(0,in,2)
        pdfgxxa= pdf(0,in,1)
        pdfqxa = 0d0
        do i=-nfl,nfl
          if(i.ne.0) pdfqxa = pdfqxa + pdf(i,in,2)
        enddo


! ! ! ! ! cccccccccccccccccccccccccccc
! ! ! ! ! ccccccc Ninh test
!       pdfqxa = 0d0
!       pdfgxa = 0d0
!       pdfgxxa = 0d0
! ! ! ! ! ccccccc Ninh test
! ! ! ! ! cccccccccccccccccccccccccccc

        if(alphaDip.eq.1d0) then
          QCDVVjjFinColltemp = (2*Born*(-((12*omxxa*pdfqxa+pdfgxxa*(-150+27*logomxxa**2+8*nfl+1
     &   8*Pi**2))*(-1+x)*x**2)-6*logx*omxxa*(9*pdfgxa*(1-x+x**2)**2-
     &   2*pdfqxa*(-2+4*x-3*x**2+x**3))+6*logomx*omxxa*(-9*pdfgxxa*x*
     &   *2+9*pdfgxa*(1-x+x**2)**2-2*pdfqxa*(-2+4*x-3*x**2+x**3)))+(1
     &   8*omxxa*x*(pdfgxa-pdfgxxa*x)*(gamma(4)*T(1,1)*T(1,4)*T(3,3)+
     &   gamma(3)*T(1,1)*T(1,3)*T(4,4)+6*logMU2oSij(1,2)*T(1,2)*T(3,3
     &   )*T(4,4)-2*logomx*T(1,1)*T(1,2)*T(3,3)*T(4,4)+6*logMU2oSij(1
     &   ,3)*T(1,3)*T(3,3)*T(4,4)+6*logMU2oSij(1,4)*T(1,4)*T(3,3)*T(4
     &   ,4))+(-1+x)*(4*omxxa*(-2*pdfqxa*(2-2*x+x**2)+9*pdfgxa*(-1+2*
     &   x-x**2+x**3))*(3*logMU2oSij(1,2)*T(1,2)-logomx*T(1,1)*T(1,2)
     &   +3*logMU2oSij(1,3)*T(1,3)+3*logMU2oSij(1,4)*T(1,4))*T(3,3)*T
     &   (4,4)-3*pdfgxxa*x**2*(6*gamma(4)*T(1,1)*T(1,4)*T(3,3)+6*gamm
     &   a(3)*T(1,1)*T(1,3)*T(4,4)+(33+36*logomxxa-2*nfl)*logMU2oSij(1
     &   ,2)*T(1,2)*T(3,3)*T(4,4)-6*logomxxa**2*T(1,1)*T(1,2)*T(3,3)*
     &   T(4,4)+2*Pi**2*T(1,1)*T(1,2)*T(3,3)*T(4,4)+33*logMU2oSij(1,3
     &   )*T(1,3)*T(3,3)*T(4,4)-2*nfl*logMU2oSij(1,3)*T(1,3)*T(3,3)*T(
     &   4,4)+33*logMU2oSij(1,4)*T(1,4)*T(3,3)*T(4,4)-2*nfl*logMU2oSij
     &   (1,4)*T(1,4)*T(3,3)*T(4,4)+6*logomxxa*(6*logMU2oSij(1,3)*T(1
     &   ,3)*T(3,3)*T(4,4)+6*logMU2oSij(1,4)*T(1,4)*T(3,3)*T(4,4)+T(1
     &   ,1)*(gamma(4)*T(1,4)*T(3,3)+gamma(3)*T(1,3)*T(4,4))))))/(T(1
     &   ,1)*T(3,3)*T(4,4)))/(18.*omx**2*omxxa*x)
        else
          QCDVVjjFinColltemp = ((4*pdfqxa*(Born*(x**2+(logalpha+logomx-logx)*(2-2*x+x**2))-((-
     &   logalpha+logax+logomx)*(2-2*x+x**2)*T(1,2))/3.+((2-2*x+x**2)
     &   *(logMU2oSij(1,2)*T(1,2)+logMU2oSij(1,3)*T(1,3)+logMU2oSij(1
     &   ,4)*T(1,4)))/T(1,1)))/x+3*pdfgxa*((-6*Born*logx)/omx+(6*(1-2
     &   *x+x**2-x**3)*(logMU2oSij(1,2)*T(1,2)+logMU2oSij(1,3)*T(1,3)
     &   +logMU2oSij(1,4)*T(1,4)))/(x*T(1,1))+(6*Born*((logalpha+logo
     &   mx-logx)*(1-3*x+3*x**2-2*x**3+x**4)+x*(logalpha-Log(alphaDip
     &   +omx)+Log(2-x))+x*(logomx-Log(2-x))*Theta(-alphaDip+omx)))/(
     &   omx*x)+(2*T(1,2)*(logalpha-(logax*(omx+x**2)**2+(-logalpha+l
     &   ogomx)*(1-3*x+3*x**2-2*x**3+x**4))/x-Log(alphaDip+omx)+Log(2
     &   -x)*Theta(-1+alphaDip+x)))/omx+(6*(Born*logomx*T(1,1)+logMU2
     &   oSij(1,2)*T(1,2)+logMU2oSij(1,3)*T(1,3)+logMU2oSij(1,4)*T(1,
     &   4))*T(3,3)*T(4,4)+T(1,1)*(gamma(4)*T(1,4)*T(3,3)+gamma(3)*T(
     &   1,3)*T(4,4)-2*logomx*T(1,2)*T(3,3)*T(4,4))*Theta(-1+alphaDip
     &   +x))/(omx*T(1,1)*T(3,3)*T(4,4)))+3*pdfgxxa*x*(-((6*(Born*log
     &   omx*T(1,1)+logMU2oSij(1,2)*T(1,2)+logMU2oSij(1,3)*T(1,3)+log
     &   MU2oSij(1,4)*T(1,4))*T(3,3)*T(4,4)+T(1,1)*(gamma(4)*T(1,4)*T
     &   (3,3)+gamma(3)*T(1,3)*T(4,4)-2*logomx*T(1,2)*T(3,3)*T(4,4))*
     &   Theta(-1+alphaDip+x))/(omx*T(1,1)*T(3,3)*T(4,4)))+(-2*Born*(
     &   -3*logomxxa**2+Pi**2)+(12*logomxxa*(logMU2oSij(1,2)*T(1,2)+l
     &   ogMU2oSij(1,3)*T(1,3)+logMU2oSij(1,4)*T(1,4)))/T(1,1)+(-201*
     &   Born-99*alphaDip*Born+108*Born*logalpha**2+Born*logalpha*(99
     &   -6*nfl)+10*Born*nfl+6*alphaDip*Born*nfl+54*Born*Pi**2+6*Pi**2*T
     &   (1,2)+(3*(33-2*nfl)*logMU2oSij(1,2)*T(1,2))/T(1,1)+(99*logMU2
     &   oSij(1,3)*T(1,3))/T(1,1)-(6*nfl*logMU2oSij(1,3)*T(1,3))/T(1,1
     &   )+(99*logMU2oSij(1,4)*T(1,4))/T(1,1)-(6*nfl*logMU2oSij(1,4)*T
     &   (1,4))/T(1,1)+(18*alphaDip*gamma(3)*T(1,3))/T(3,3)+(18*alpha
     &   Dip*gamma(4)*T(1,4))/T(4,4))/9.+2*(logalpha**2-logomxxa**2)*
     &   T(1,2)*Theta(-1+alphaDip+xxa)+2*(-logalpha+logomxxa)*((gamma
     &   (3)*T(1,3))/T(3,3)+(gamma(4)*T(1,4))/T(4,4))*Theta(-1+alphaD
     &   ip+xxa))/(2.*omxxa)))/(3.*omx)
        endif
        

! ! ! cccccccccccccccccccccccccccccc
! ! ! cccccccccc Ninh test
! ! ! cccccccc only Operator P
!         QCDVVjjFinColltemp = ((-((33 + 36*logomxxa - 2*nfl)*pdfgxxa*(-1 + x)*x**2) + 
!      -      4*omxxa*(-9*pdfgxxa*x**2 + 9*pdfgxa*(1 - x + x**2)**2 - 2*pdfqxa*(-2 + 4*x - 3*x**2 + x**3)))*
!      -    (logMU2oSij(1,2)*T(1,2) + logMU2oSij(1,3)*T(1,3) + logMU2oSij(1,4)*T(1,4)))/(6.*omx**2*omxxa*x*T(1,1))
! !       QCDVVjjFinColl = sum(T(1,2:4))
! ! ! c only Operator K
!           QCDVVjjFinColltemp = (Born*(-((12*omxxa*pdfqxa + pdfgxxa*(-150 + 27*logomxxa**2 + 8*nfl + 18*Pi**2))*(-1 + x)*x**2) - 
!      -       6*logx*omxxa*(9*pdfgxa*(1 - x + x**2)**2 - 2*pdfqxa*(-2 + 4*x - 3*x**2 + x**3)) + 
!      -       6*logomx*omxxa*(-9*pdfgxxa*x**2 + 9*pdfgxa*(1 - x + x**2)**2 - 2*pdfqxa*(-2 + 4*x - 3*x**2 + x**3)))*T(3,3)*T(4,4) - 
!      -    2*logomx*omxxa*(-9*pdfgxxa*x**2 + 9*pdfgxa*(1 - x + x**2)**2 - 
!      -    2*pdfqxa*(-2 + 4*x - 3*x**2 + x**3))*T(1,2)*T(3,3)*T(4,4) + 
!      -    3*x*(-3*logomxxa**2*omx*pdfgxxa*x*T(1,2)*T(3,3)*T(4,4) + 3*logomxxa*omx*pdfgxxa*x*(gamma(4)*T(1,4)*T(3,3) + 
!      -    gamma(3)*T(1,3)*T(4,4)) + 
!      -       3*omxxa*(pdfgxa - pdfgxxa*x)*(gamma(4)*T(1,4)*T(3,3) + gamma(3)*T(1,3)*T(4,4)) + 
!      -       omx*pdfgxxa*x*(3*gamma(4)*T(1,4)*T(3,3) + 3*gamma(3)*T(1,3)*T(4,4) + Pi**2*T(1,2)*T(3,3)*T(4,4))))/
!      -  (9.*omx**2*omxxa*x*T(3,3)*T(4,4))


! !     QCDVVjjFinColltemp = 0d0

! ! ! cccccccccc Ninh test
! ! ! ccccccccccccccccccccccccccccccc
      else
        pdfqxa = pdf(beamsign(in)*iflav(in),in,2)
        pdfgxa = pdf(0,in,2)
        pdfqxxa= pdf(beamsign(in)*iflav(in),in,1)

! cccccccccccccccccccccccccccc
! ccccccc Ninh test

!       pdfgxa = 0d0
!       pdfqxa = 0d0
!       pdfqxxa= 0d0
! !       
! ccccccc Ninh test
! cccccccccccccccccccccccccccc

        if(alphaDip.eq.1d0) then
          QCDVVjjFinColltemp = (3*pdfgxa*(Born*(omx*x+(logomx-logx)*(0.5-x+x**2))-(3*logomx*(1
     &   -2*x+2*x**2)*T(1,2))/8.+((1-2*x+2*x**2)*(logMU2oSij(1,2)*T(1
     &   ,2)+logMU2oSij(1,3)*T(1,3)+logMU2oSij(1,4)*T(1,4)))/(2.*T(1,
     &   1)))+pdfqxxa*x*(-((8*Born*logomx-6*logomx*T(1,2)+(8*logMU2oS
     &   ij(1,2)*T(1,2))/T(1,1)+(8*logMU2oSij(1,3)*T(1,3))/T(1,1)+(8*
     &   logMU2oSij(1,4)*T(1,4))/T(1,1)+(3*gamma(3)*T(1,3))/T(3,3)+(3
     &   *gamma(4)*T(1,4))/T(4,4))/omx)+(4*Born*(logomxxa**2-Pi**2/3.
     &   )+4*Born*(-5+Pi**2)-3*logomxxa**2*T(1,2)+Pi**2*T(1,2)+(2*(3+
     &   4*logomxxa)*(logMU2oSij(1,2)*T(1,2)+logMU2oSij(1,3)*T(1,3)+l
     &   ogMU2oSij(1,4)*T(1,4)))/T(1,1)+(3*gamma(3)*T(1,3))/T(3,3)+3*
     &   logomxxa*((gamma(3)*T(1,3))/T(3,3)+(gamma(4)*T(1,4))/T(4,4))
     &   +(3*gamma(4)*T(1,4))/T(4,4))/omxxa)+(pdfqxa*(3*gamma(4)*T(1,
     &   1)*T(1,4)*T(3,3)+3*gamma(3)*T(1,1)*T(1,3)*T(4,4)+4*Born*(1+l
     &   ogomx-logx-2*x+(1+logomx-logx)*x**2)*T(1,1)*T(3,3)*T(4,4)+4*
     &   (1+x**2)*logMU2oSij(1,2)*T(1,2)*T(3,3)*T(4,4)-3*logomx*T(1,1
     &   )*T(1,2)*T(3,3)*T(4,4)+4*logMU2oSij(1,3)*T(1,3)*T(3,3)*T(4,4
     &   )+4*logMU2oSij(1,4)*T(1,4)*T(3,3)*T(4,4)+x**2*(-3*logomx*T(1
     &   ,1)*T(1,2)+4*logMU2oSij(1,3)*T(1,3)+4*logMU2oSij(1,4)*T(1,4)
     &   )*T(3,3)*T(4,4)))/(omx*T(1,1)*T(3,3)*T(4,4)))/(3.*omx)
        else
          QCDVVjjFinColltemp = (3*pdfgxa*(Born*(omx*x+(logalpha+logomx-logx)*(0.5-x+x**2))-(3*
     &   (-logalpha+logax+logomx)*(1-2*x+2*x**2)*T(1,2))/8.+((1-2*x+2
     &   *x**2)*(logMU2oSij(1,2)*T(1,2)+logMU2oSij(1,3)*T(1,3)+logMU2
     &   oSij(1,4)*T(1,4)))/(2.*T(1,1)))+pdfqxa*((-8*Born*logx)/omx-(
     &   4*(1+x)*(logMU2oSij(1,2)*T(1,2)+logMU2oSij(1,3)*T(1,3)+logMU
     &   2oSij(1,4)*T(1,4)))/T(1,1)+(8*(Born*logomx*T(1,1)+logMU2oSij
     &   (1,2)*T(1,2)+logMU2oSij(1,3)*T(1,3)+logMU2oSij(1,4)*T(1,4))*
     &   T(3,3)*T(4,4)+3*T(1,1)*(gamma(4)*T(1,4)*T(3,3)+gamma(3)*T(1,
     &   3)*T(4,4)-2*logomx*T(1,2)*T(3,3)*T(4,4))*Theta(-1+alphaDip+x
     &   ))/(omx*T(1,1)*T(3,3)*T(4,4))+(4*Born*(1+logalpha-logomx+log
     &   x-2*x+(1+logalpha+logomx-logx)*x**2+2*Log(2-x)-2*Log(1+alpha
     &   Dip-x)+2*Log((-1+x)/(-2+x))*Theta(1-alphaDip-x))+3*T(1,2)*(l
     &   ogalpha-logax+logomx+(logalpha-logax-logomx)*x**2-2*Log(alph
     &   aDip+omx)+2*Log(2-x)*Theta(-1+alphaDip+x)))/omx)+pdfqxxa*x*(
     &   -((8*(Born*logomx*T(1,1)+logMU2oSij(1,2)*T(1,2)+logMU2oSij(1
     &   ,3)*T(1,3)+logMU2oSij(1,4)*T(1,4))*T(3,3)*T(4,4)+3*T(1,1)*(g
     &   amma(4)*T(1,4)*T(3,3)+gamma(3)*T(1,3)*T(4,4)-2*logomx*T(1,2)
     &   *T(3,3)*T(4,4))*Theta(-1+alphaDip+x))/(omx*T(1,1)*T(3,3)*T(4
     &   ,4)))+(6*Born*logalpha+8*Born*logalpha**2+4*Born*(logomxxa**
     &   2-Pi**2/3.)+2*Born*(-7-3*alphaDip+2*Pi**2)+Pi**2*T(1,2)+(2*(
     &   3+4*logomxxa)*(logMU2oSij(1,2)*T(1,2)+logMU2oSij(1,3)*T(1,3)
     &   +logMU2oSij(1,4)*T(1,4)))/T(1,1)+(3*alphaDip*gamma(3)*T(1,3)
     &   )/T(3,3)+(3*alphaDip*gamma(4)*T(1,4))/T(4,4)+3*(logalpha**2-
     &   logomxxa**2)*T(1,2)*Theta(-1+alphaDip+xxa)+3*(-logalpha+logo
     &   mxxa)*((gamma(3)*T(1,3))/T(3,3)+(gamma(4)*T(1,4))/T(4,4))*Th
     &   eta(-1+alphaDip+xxa))/omxxa))/(3.*omx)
        endif
         
cccccccccccccccccc
ccc  only P 
! !     QCDVVjjFinColltemp = -((4*(3+4*logomxxa)*pdfqxxa*(-1+x)*x+omxxa*(-8*(pdfqxa-2*pdfqxxa*x+pdfqxa*x**2)+3*pdfgxa*(-1+
! !      &   3*x-4*x**2+2*x**3)))*(logMU2oSij(1,2)*T(1,2)+logMU2oSij(1,3)*T(1,3)+logMU2oSij(1,4)*T(1,4)
! !      &   ))/(6.*omx**2*omxxa*T(1,1))
! ! ccc only K 
! !     QCDVVjjFinColltemp = (-4*Born*(-2*(-1+x)*(-4*pdfqxxa*(-15+3*logomxxa**2+2*Pi**2)*x+3*omxxa*(-1+x)*(4*pdfqxa+3*pdfg
! !      &   xa*x))-3*logx*omxxa*(-8*pdfqxa*(1+x**2)+3*pdfgxa*(-1+3*x-4*x**2+2*x**3))+3*logomx*omxxa*(-
! !      &   8*(pdfqxa-2*pdfqxxa*x+pdfqxa*x**2)+3*pdfgxa*(-1+3*x-4*x**2+2*x**3)))*T(3,3)*T(4,4)+3*(3*lo
! !      &   gomx*omxxa*(-8*(pdfqxa-2*pdfqxxa*x+pdfqxa*x**2)+3*pdfgxa*(-1+3*x-4*x**2+2*x**3))*T(1,2)*T(
! !      &   3,3)*T(4,4)+8*(-3*logomxxa**2*omx*pdfqxxa*x*T(1,2)*T(3,3)*T(4,4)+3*logomxxa*omx*pdfqxxa*x*
! !      &   (gamma(4)*T(1,4)*T(3,3)+gamma(3)*T(1,3)*T(4,4))+3*omxxa*(pdfqxa-pdfqxxa*x)*(gamma(4)*T(1,4
! !      &   )*T(3,3)+gamma(3)*T(1,3)*T(4,4))+omx*pdfqxxa*x*(3*gamma(4)*T(1,4)*T(3,3)+3*gamma(3)*T(1,3)
! !      &   *T(4,4)+Pi**2*T(1,2)*T(3,3)*T(4,4)))))/(72.*omx**2*omxxa*T(3,3)*T(4,4))
cccc
cccccccccccccccccc


cccccccccccccccccc
ccc Ninh compare
! !     if (qtype.ne.0) then
! !       count = count+1
! !       print*,count
! !       QCDVVjjFinColltemp = 0d0
! !     endif
ccc Ninh compare
cccccccccccccccccc

      endif


c calc born instead
! !       if(iflav(in).eq.0) then
! !       QCDVVjjFinColltemp = pdfgxxa*Born/2d0/(als(1,1)/2d0/pi)/omx*x/omxxa  ! this is the born, factor 2, because we caount each subpr twice (loop), omx, which is included in mathematica and cancelled with *omx in a few lines
! !       else
! !       QCDVVjjFinColltemp = pdfqxxa*Born/2d0/(als(1,1)/2d0/pi)/omx*x/omxxa  ! this is the born
! !       endif
      
      QCDVVjjFinColltemp = QCDVVjjFinColltemp*als(1,1)/2d0/pi*omx ! 1/omx was included in mathematica file
      QCDVVjjFinCollBornKin = QCDVVjjFinCollBornKin + QCDVVjjFinColltemp*pdfb/x*omxxa ! dx/x, within real emission this 1/x is put to the pa*pb

      enddo
      
      QCDVVjjFinCollBornKin = QCDVVjjFinCollBornKin/4d0/9d0 ! 2nd pdf, col/spin average of born amplitude (assuming q initial state)
      
      if(iflav(1).eq.0) QCDVVjjFinCollBornKin = QCDVVjjFinCollBornKin*3d0/8d0
      if(iflav(2).eq.0) QCDVVjjFinCollBornKin = QCDVVjjFinCollBornKin*3d0/8d0
      if (iflav(3).eq.iflav(4)) QCDVVjjFinCollBornKin = QCDVVjjFinCollBornKin/2d0
c      print*,'finColl',qcdvvjjfincollBornKin
      end


! 
