c General routine to calculate Catani Seymour dipoles
c Author: Matthias Kerner
c 
c xyuz: splitting parameters (1,*)= x or y, (2,*) = u or z.
c splitting: "qg", "qq", or "gg"
c emit1,emit2,spect: phase space index in N+1 kinematics
c hel*: only hel=0 implemented -> helicity summed dipoles
c idB,fsignB: conversion from diagram to ps labels for born process
c emitB,spectB: diagram labels for born ME
c BornME: subroutine to calc Born ME, will be called with
c      call BornME(p,kin,idB,fsignB,emitB,spectB,buff,born)
c    external gluons have to be calculated with polvec
c    born(0,*) |ME|^2
c    born(1:2,*) |ME|^2 with polarization of emitB fixed
c    born(-1,*) 2*Real(<1|TiTj|2>)
c    born(:,1:Nres) allows to calculate the dipoles with different born amplitudes (e.g. for different quark flavours)
c    fsignB: sign factor for the born kinematic;
c             if emitB is a quark, fsignB(emitB) should be the sign of the quark in N+1 kinematics, sign will be changed for initial state splitting g->q qb
c                           gluon, fsignB(emitB) will be assigned here
c    buff(*), additional int array that allows to pass additional information to the BornME, such as process id or helicity flags
c Alpha Parameter can be set with       entry CSGeneralSetAlpha(alpha_dip)
      subroutine CSgeneral(p,kin,xyuz,splitting,emit1,emit2,spect,hel1,hel2,
     &                     BornME,idB,fsignB,emitB,spectB,buff,Nres,res)
      implicit none
#include "global.inc"
      real*8 p(0:3,1:max_p,1:max_kin),res(*),xyuz(2,2:max_kin)
      integer idB(*),fsignb(*),emit1,emit2,spect,hel1,hel2,kin,Nres,emitB,spectB,buff(*)
      character*2 splitting
      
      integer em1,em2,i,j,signBackup,mu
      real*8 p2,dotrr,factor(-1:2),CF,CA,TR,x,omx,z,omz,born(-1:2,Nres),hVec(0:3),pipj,pipk,pjpk,hReal,dotrrQ
      real*8 alpha/1d0/, alpha_dip
      save alpha
      complex*16 dotcc,dotcr,dotrc
      real*8 epsMu(0:3,-1:2)
      external dotrr,dotcc,dotcr,dotrc,dotrrQ
      parameter (CF = 4d0/3d0,CA=3d0,TR=0.5d0)
      
      goto 33
      
      entry CSGeneralSetAlpha(alpha_dip)
      alpha = alpha_dip
      return
      
 33   continue

      do i=-1,2
        factor(i) = 0d0
      enddo
      
      em1 = min(emit1,emit2)
      em2 = max(emit1,emit2)

      if(em2.le.2) then
	print*,"CSgeneral called with em1,em2 in initial state"
	stop
      endif
      
c set conventions of ptile.F
      if(em1.ge.3.and.spect.ge.3) then ! ff kin
	if(idB(emitB).ne.3) then
	  i=idB(emitB)
	  idB(emitB)=idB(spectB)
	  idB(spectB)=i
	endif
      endif
      if(em1.ge.3.and.spect.le.2) then !fi
	if(idB(emitB).ne.3) then
	  do j=1,4
	    if(idB(j).eq.3) i=j
	  enddo
	  idB(i) = 4
	  idB(emitB)=3
	endif
      endif
      if(em1.le.2.and.spect.ge.3) then !if
	if(idB(spectB).ne.3) then
	  do j=1,4
	    if(idB(j).eq.3) i=j
	  enddo
	  idB(i) = 4
	  idB(spectB)=3
	endif
      endif
	  
      
      signBackup = fsignB(emitB)
      if(splitting.eq."qg") then
	if (emit2.le.2) fsignB(emitB) = -fsignB(emitB) ! initial state g -> q qb
      else
	if (em1.le.2) then
	  fsignB(emitB) = -1
	else 
	  fsignB(emitB) = 1
	endif
      endif
      
      p2 = dotrr(p(0,em1,1),p(0,em2,1))
      x= xyuz(1,kin)
      if(x.ge.0d0) then
	omx = 1d0-x
      else
	omx = -x
	x = 1d0-omx
      endif
      z = xyuz(2,kin)
      if(z.ge.0d0) then
	omz = 1d0-z
      else
	omz = -z
	z = 1d0-omz
      endif
c change z for conventions in ptilde
      if (    em1.le.2.and.spect.ge.3.and.em2.gt.spect
     $    .or.em1.ge.3.and.emit1.ge.emit2) then
        hReal=z
        z=omz
        omz=hReal
      endif
        
c check alpha parameter
      if(alpha.lt.1d0) then
        if(em1.le.2) then
         if(  (spect.ge.3.and.z.gt.alpha) .or. (spect.le.2.and.p2/dotrr(p(0,1,1),p(0,2,1)).gt.alpha)) then
            do i=1,nres
              res(i)=0d0
            enddo
            fsignB(emitB)=signBackup
            return
          endif
        else
          if( (spect.le.2.and.omx.gt.alpha) .or. (spect.ge.3.and.x.gt.alpha)) then
            do i=1,nres
              res(i)=0d0
            enddo
            fsignB(emitB)=signBackup
            return
          endif
        endif
      endif

      if(splitting.eq."qg") then
	if(em1.le.2) then
	  if(spect.le.2) then
c ii dipoles
	    if(emit1.le.2) then ! q->q g_i
	      if(hel1.eq.0) then
		factor(0)=8d0*pi*(1d0+x*x)/omx
	      endif
	    else ! g-> q q_i
	      if(hel1.eq.0) then
		factor(0)=8d0*pi*TR/CF*(1d0-2d0*x*omx)*8d0/3d0
	      endif
	    endif
	  else
c if dipoles
	    if(emit1.le.2) then ! q->q g_i
	      if(hel1.eq.0) then
		factor(0)=8d0*pi*(2d0/(omx+z)-x-1d0)
	      endif
	    else ! g-> q q_i
	      if(hel1.eq.0) then
		factor(0)=8d0*pi*TR/CF*(1d0-2d0*x*omx)*8d0/3d0
	      endif
	    endif
	  endif
	else
	  if(spect.le.2) then
c fi dipoles
	    if(hel1.eq.0) then
	      factor(0) = 8d0*pi*(2d0/(omz+omx)-z-1d0)
	    endif
	  else
c ff dipoles
	    if(hel1.eq.0) then
	      factor(0) = 8d0*pi*(2d0/(1d0-z*omx)-z-1d0)
	    endif
	  endif
	endif
      elseif(splitting.eq."gg") then
        call polvec(p(0,idB(emitB),kin),1,epsMu(0,1))
        call polvec(p(0,idB(emitB),kin),2,epsMu(0,2))
        if(em1.le.2) then
          if(spect.le.2) then
c ii dipoles
            pipj=dotrr(p(0,em1,1),p(0,em2,1))
            pipk=dotrr(p(0,1,1),p(0,2,1))
            pjpk=dotrr(p(0,em2,1),p(0,spect,1))
            
            hVec(0) = p(0,em2,1) - pipj/pipk*p(0,spect,1)
            hVec(1) = p(1,em2,1)
            hVec(2) = p(2,em2,1)
            hVec(3) = p(3,em2,1) - pipj/pipk*p(3,spect,1)

            hReal = omx/x*pipk/pipj/pjpk
            factor(0) = x/omx+x*omx
          else
c if dipole
            pjpk=dotrr(p(0,em2,1),p(0,spect,1))

            do mu =0,3
              hVec(mu) = p(mu,em2,1)/z-p(mu,spect,1)/omz
            enddo
            hReal = omx/x*z*omz/pjpk
            factor(0) = 1d0/(omx+z)-1d0+x*omx
          endif
        else
          if (spect.le.2) then
c fi dipole
            do mu=0,3
              hVec(mu) = p(mu,em1,1)*z-p(mu,em2,1)*omz
            enddo
            hReal = 1d0/p2
            factor(0) = 1d0/(omx+z)+1d0/(omx+omz)-2d0
          else
c ff dipole
            do mu=0,3
              hVec(mu) = p(mu,em1,1)*z-p(mu,em2,1)*omz
            enddo
            hReal = 1d0/p2
            factor(0) = 1d0/(1d0-z*omx)+1d0/(1d0-omz*omx)-2d0
          endif
        endif
c common for all gg splittings
        factor(1) = dotrr(epsMu(0,1),hVec)
        factor(2) = dotrr(epsMu(0,2),hVec)
        factor(-1)= factor(1)*factor(2)*hReal
        factor(1) = factor(1)*factor(1)*hReal
        factor(2) = factor(2)*factor(2)*hReal
        do i=-1,2
          factor(i)=factor(i)*16d0*pi
        enddo

c 
c qq 
c
      elseif(splitting.eq."qq") then
        call polvec(p(0,idB(emitB),kin),1,epsMu(0,1))
        call polvec(p(0,idB(emitB),kin),2,epsMu(0,2))
        if(em1.le.2) then
          if(spect.le.2) then
c ii dipoles
            pipj=dotrr(p(0,em1,1),p(0,em2,1))
            pipk=dotrr(p(0,1,1),p(0,2,1))
            pjpk=dotrr(p(0,em2,1),p(0,spect,1))
            
            hVec(0) = p(0,em2,1) - pipj/pipk*p(0,spect,1)
            hVec(1) = p(1,em2,1)
            hVec(2) = p(2,em2,1)
            hVec(3) = p(3,em2,1) - pipj/pipk*p(3,spect,1)

            hReal = 2d0*omx/x*pipk/pipj/pjpk
            factor(0) = x
          else
c if dipole
            pjpk=dotrr(p(0,em2,1),p(0,spect,1))

            do mu=0,3
              hVec(mu) = p(mu,em2,1)/z-p(mu,spect,1)/omz
            enddo
            hReal = 2d0*omx/x*z*omz/pjpk
            factor(0) = x
          endif
        else
          if (spect.le.2) then
c fi dipole
            do mu=0,3
              hVec(mu)= p(mu,em1,1)*z-p(mu,em2,1)*omz
            enddo
            hReal = -2d0/p2
            factor(0) = 1d0
          else
c ff dipole
            do mu=0,3
              hVec(mu) = p(mu,em1,1)*z-p(mu,em2,1)*omz
            enddo
            hReal = -2d0/p2
            factor(0) = 1d0
          endif
        endif
c common for all gg splittings
        factor(1) = dotrr(epsMu(0,1),hVec)
        factor(2) = dotrr(epsMu(0,2),hVec)
        factor(-1)= factor(1)*factor(2)*hReal
        factor(1) = factor(1)*factor(1)*hReal
        factor(2) = factor(2)*factor(2)*hReal
        do i=-1,2
          factor(i)=factor(i)*8d0*pi*TR/CA ! for initial states : CF*3/8=TR
        enddo
      else
	print*,"unknown splitting"
	stop
      endif
      
      if(em1.le.2.or.spect.le.2) then
        do i=-1,2
          factor(i) = factor(i)/x
        enddo
      endif
!       print*,idB(1:4),spectB,emitB,kin,fsignB(1:4),buff(1)
      call BornME(p,kin,idB,fsignB,emitB,spectB,buff,born)

!       print*,born(0,1:Nres)
      do i=1,Nres
        res(i)=0d0
      do j=-1,2 ! sum over born(0), pol states 1,2 and offdiag -1 
        res(i) =res(i) -factor(j)*born(j,i)/2d0/p2
      enddo;enddo
      
      fsignB(emitB) = signBackup

      end
      
      
