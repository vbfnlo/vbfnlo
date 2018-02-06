c******************************************************************
c
c   begin subroutine wbf_h4jq
c
c*****************************************************************
      subroutine wbf_h4jq(xi,p,v,physToDiag,fsign,bos, 
     1     nlo,lok,xuz,
     2     pdf,res,nmin,nmax)
      use globalvars, only: ldoblha
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
#include "color.inc"
      double precision p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), xi(nx), 
     $     xuz(2,2:max_kin)
      double precision pdf(-6:6,2), res(maxnumsubproc,0:max_kin)
      integer physToDiag(6), fsign(6+max_v), gsign, 
     $     bos, nlo, nmin, nmax
      logical lok(max_kin)
c
      double precision uucc_c(3),uuss_c(3),ddcc_c(3),ddss_c(3),
     $     udsc_c(3),ducs_c(3)
      double precision uuccuu(0:2),uussuu(0:2),ddccuu(0:2),ddssuu(0:2)
      double precision uuccdd(0:2),uussdd(0:2),ddccdd(0:2),ddssdd(0:2)
      double precision udsc(0:2),ducs(0:2)

      
c wbf_h4q calls the amplitude square routines 
c             qqh4qnc_mg    for qQ-->qQ qq~H     for bos=6
c for the subtraction terms for the NLO cross section calculation
c
c  INPUT:  p(0:3,np,mxd)      external physical parton momenta
c          v(0:3,nv,mxd)      for Higgs production only
c                          the sum q(mu) = v(mu,1)+...+v(mu,nv) is needed
c          physToDiag(6)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign,gsign     sign factors for amplitude calls; 
c          nlo             nlo = 0: calculate LO only (i.e. no subtraction)
c                          nlo = 1: calculate full NLO subtraction
c          lok(mxd)          lok(ID)=T means momenta set ID passes acceptance
c                          cuts, i.e res(k,ID) needs to be calculated
c
c     need to put in nlo switch
c  OUTPUT:

c  res(k,nd)            
c                            
c  In and Output
c    nmin, nmax            range of process numbers for this call
c
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
c
      integer ip(2:max_kin),kp(2:max_kin),jp(2:max_kin),nd
      double precision pbar(0:3,6+max_v),qbar(0:4)
      double precision kbar(0:3,4+max_v)
      double precision sub(4,2:max_kin),subv(2:max_kin),qvec(0:3,2:max_kin)
      double precision dip,msq_b(0:2,2:max_kin),msq_v(0:2,2:max_kin)
      double precision NCmatrix_r(0:1,0:1,0:1,1:2,0:1)
      double precision CCmatrix_r(0:1,0:1,0:1,1:4,0:1)
      double precision NCmatrix_b(0:1,0:1,2,2:max_kin)
      double precision CCmatrix_b(0:1,2,2:max_kin)
      double precision NCmatrix_v(0:1,0:1,2,2:max_kin)
      double precision CCmatrix_v(0:1,2,2:max_kin)
      double precision polcolqq,polcolgq
      double precision sum,delta,newdelta
      integer deltapair
      integer iflav(6,4),diagToPhys(6),FL_H4j
      double precision dot0p
      external FL_H4j
      external dot0p
      logical ldebug, ldebug1
      parameter(ldebug=.false.,ldebug1=.false.)
      logical ldebugmg
      parameter (ldebugmg=.false.)
      integer i,if1,if2,if3,j,mu,jd,ndmap(16)
      integer k
      integer interference(4)
      integer if1l,if2l,icol
      integer fc(4+max_v)
      double precision symmfac(4)
      integer ndmax
      integer gq,qq,gg,qg
      parameter(qq = 1, qg = 2, gq = 3, gg = 4)
      logical lsymmcontrib, linterference
      parameter (lsymmcontrib = .true.)
      parameter (linterference = .true.)

      polcolqq = 1d0/(NC**2 * 4d0 * xi(1)*xi(2)) 
      polcolgq = 1d0/(NC*(NC**2-1d0) * 4d0 * xi(1)*xi(2))
c     set up ndmap
      if (lsymmcontrib) then
        ! independent of fsign
        ndmap( 1)=4   ! (56)-3
        ndmap( 2)=14  ! (56)-4
        ndmap( 3)=7   ! (56)-1
        ndmap( 4)=17  ! (56)-2
        ndmap( 9)=23  ! (13)-5
        ndmap(10)=24  ! (13)-6
        ndmap(13)=25  ! (24)-5
        ndmap(14)=22  ! (24)-6
        ! upper line
        if (fsign(1).gt.0) then
          ndmap( 5)=3   ! (53)-6
          ndmap( 6)=6   ! (53)-1
          ndmap(11)=8   ! (16)-3
          ndmap(12)=10  ! (16)-5
        else
          ndmap( 5)=2   ! (63)-5
          ndmap( 6)=5   ! (63)-1
          ndmap(11)=9   ! (15)-3
          ndmap(12)=11  ! (15)-6
        endif
        ! lower line
        if (fsign(3).gt.0) then
          ndmap( 7)=13  ! (54)-6
          ndmap( 8)=16  ! (54)-2
          ndmap(15)=18  ! (26)-4
          ndmap(16)=20  ! (26)-5
        else
          ndmap( 7)=12  ! (64)-5
          ndmap( 8)=15  ! (64)-2
          ndmap(15)=19  ! (25)-4
          ndmap(16)=21  ! (25)-6
        endif
        ndmax=16
      else
        ndmap(1)=4
        ndmap(2)=14
        ndmap(3)=7
        ndmap(4)=17
        ndmap(5)=23
        ndmap(6)=24
        ndmap(7)=25
        ndmap(8)=22
        ndmax=8
      endif
c     ave over polarizations and colors
c     symmetry factor for two gluons
      if(ldebug1) then
         print*,"momentum"
         do i = 1,n_p
            print*,"p(mu,i,1)=",(p(mu,i,1),mu=0,3)
         enddo
         print*,"phystodiag",physToDiag(1:6)
         print*,"fsign",fsign(1:6)
      endif
c            
      nmin = nmax+1
c
      do i = 1,6
         diagToPhys(physToDiag(i)) = i
      enddo
         
c get the real emission amplitude squared, 
      do mu = 0,3
         do i = 1,6
            pbar(mu,physToDiag(i))=p(mu,i,1)
         enddo         
      enddo
      
      do mu = 0,3
         pbar(mu,7) = 0         ! dummy momentum
         pbar(mu,8) = 0
         pbar(mu,9) = 0
         pbar(mu,10) = 0
         do i = 1,n_v
            pbar(mu,8) = pbar(mu,8) + v(mu,i,1) ! Higgs momentum
         enddo
      enddo
      fsign(7) = 0
      fsign(8) = 1
      fsign(9) = 0
      fsign(10) = 0

      if (lok(1)) then
        call qqh4q(pbar,fsign,NCmatrix_r,CCmatrix_r,lsymmcontrib)
        if (ldebugmg) then
          call qqh4qnc_mg_c(pbar,fsign,
     $      uuccuu,uussuu,ddccuu,ddssuu,
     $      uuccdd,uussdd,ddccdd,ddssdd)
          write(*,*) "Hjjj 4q NC Comparison:"
          do i=1,2                  !loop over color flow 
            write(*,*) "  colour structure: c = ", i
            write(*,*) "  NC(0,0,0,c)", NCmatrix_r(0,0,0,i,0), uuccuu(i), abs(1-uuccuu(i)/NCmatrix_r(0,0,0,i,0))
            write(*,*) "  NC(0,1,0,c)", NCmatrix_r(0,1,0,i,0), uussuu(i), abs(1-uussuu(i)/NCmatrix_r(0,1,0,i,0))
            write(*,*) "  NC(1,0,0,c)", NCmatrix_r(1,0,0,i,0), ddccuu(i), abs(1-ddccuu(i)/NCmatrix_r(1,0,0,i,0))
            write(*,*) "  NC(1,1,0,c)", NCmatrix_r(1,1,0,i,0), ddssuu(i), abs(1-ddssuu(i)/NCmatrix_r(1,1,0,i,0))
            write(*,*) "  NC(0,0,1,c)", NCmatrix_r(0,0,1,i,0), uuccdd(i), abs(1-uuccdd(i)/NCmatrix_r(0,0,1,i,0))
            write(*,*) "  NC(0,1,1,c)", NCmatrix_r(0,1,1,i,0), uussdd(i), abs(1-uussdd(i)/NCmatrix_r(0,1,1,i,0))
            write(*,*) "  NC(1,0,1,c)", NCmatrix_r(1,0,1,i,0), ddccdd(i), abs(1-ddccdd(i)/NCmatrix_r(1,0,1,i,0))
            write(*,*) "  NC(1,1,1,c)", NCmatrix_r(1,1,1,i,0), ddssdd(i), abs(1-ddssdd(i)/NCmatrix_r(1,1,1,i,0))
          enddo

          call qqh4qcc_mg_c(pbar,fsign,udsc,ducs)
          write(*,*) "Hjjj 4q CC Comparison:"
          do i=1,2                  !loop over color flow 
            write(*,*) "  CC(0,c)", CCmatrix_r(0,1,0,i,0), udsc(i), abs(1-udsc(i)/CCmatrix_r(0,1,0,i,0))
            write(*,*) "  CC(1,c)", CCmatrix_r(1,0,0,i,0), ducs(i), abs(1-ducs(i)/CCmatrix_r(1,0,0,i,0))
          enddo
        endif
      endif

      if (ldoblha) then
        return
      endif

c for the NLO case get the subtraction terms; 

#ifdef WITH_NLO
         if (nlo.eq.1) then
c     create dipole subtraction terms
            call create_ikj(ip,kp,jp)
            do jd = 1,ndmax
               nd = ndmap(jd)
c            do nd = 2,n_kin
               call dips(nd,p,xuz,ip(nd),kp(nd),jp(nd),sub(1,nd),
     $              subv(nd),qvec(0,nd))
            enddo
c          
c     Compute born matrix elements with transformed momenta:
c     
c     
            do jd = 1,ndmax
               kbar = -999
               qbar = -999
               nd = ndmap(jd)   !mapping to dipoles set as data statement
               if (lok(nd).and.(nlo.eq.1)) then
c
                  do mu = 0,3
                     kbar(mu,5) = 0 ! dummy momentum
                     kbar(mu,6) = 0
                     kbar(mu,7) = 0
                     kbar(mu,8) = 0
                     do i = 1,n_v
                        kbar(mu,6) = kbar(mu,6) + v(mu,i,nd) ! Higgs mom.
                     enddo
                  enddo
                  fc(5) = 0
                  fc(6) = 1
                  fc(7) = 0
                  fc(8) = 0
c             
                  if(jd.le.4) then     !qQ -> qQgH born process , g from (56)
                     do mu = 0,3          
                        do i = 1,4
                           kbar(mu,physToDiag(i))=p(mu,i,nd)
                        enddo
                        qbar(mu) = p(mu,5,nd) !gluon momenta
                     enddo
                     qbar(4) = 0.0d0
                     gsign = 1
c

                     do i=1,4
                        fc(i)=fsign(i)
                     enddo
c
                  elseif(lsymmcontrib .and. jd.le.6) then !qQ -> qQgH , g from (53)
                     do mu = 0,3
                        kbar(mu,physToDiag(2))=p(mu,2,nd)
                        kbar(mu,physToDiag(4))=p(mu,4,nd)
                        kbar(mu,physToDiag(1))=p(mu,1,nd)
                        kbar(mu,physToDiag(3))=p(mu,5,nd)
                        qbar(mu) = p(mu,3,nd) !gluon momenta
                     enddo
                     qbar(4) = 0.0d0
                     gsign = 1
c

                     do i=1,4
                        fc(i)=fsign(i)
                     enddo
c
                  elseif(lsymmcontrib .and. jd.le.8) then !qQ -> qQgH , g from (54)
                     do mu = 0,3
                        kbar(mu,physToDiag(1))=p(mu,1,nd)
                        kbar(mu,physToDiag(3))=p(mu,3,nd)
                        kbar(mu,physToDiag(2))=p(mu,2,nd)
                        kbar(mu,physToDiag(4))=p(mu,5,nd)
                        qbar(mu) = p(mu,4,nd) !gluon momenta
                     enddo
                     qbar(4) = 0.0d0
                     gsign = 1
c
                     
                     do i=1,4
                        fc(i)=fsign(i)
                     enddo
c
                  elseif( (lsymmcontrib .and. jd.le.10) .or. jd.le.6 ) then !gQ->Qq"q"~H , g from (13)
                     do mu =0,3
                        qbar(mu) = p(mu,1,nd) !i.s. gluon
                        kbar(mu,phystodiag(2)) = p(mu,2,nd)!i.s. quark
                        kbar(mu,phystodiag(4)) = p(mu,3,nd)!f.s. quark
                        kbar(mu,1) = p(mu,4,nd)
                        kbar(mu,2) = p(mu,5,nd) 
                     enddo
                     qbar(4) = 0.0d0
                     gsign = -1
                     fc(1) = fsign(5)
                     fc(2) = fsign(6)
                     fc(3) = fsign(3)
                     fc(4) = fsign(4)
c
                  elseif(lsymmcontrib .and. jd.le.12) then !gQ->Qq"q"~H , g from (16)
                     do mu =0,3
                        qbar(mu) = p(mu,1,nd) !i.s. gluon
                        kbar(mu,phystodiag(2)) = p(mu,2,nd)!i.s. quark
                        kbar(mu,phystodiag(4)) = p(mu,4,nd)!f.s. quark
                        kbar(mu,phystodiag(1)) = p(mu,5,nd)
                        kbar(mu,phystodiag(3)) = p(mu,3,nd)
                     enddo
                     qbar(4) = 0.0d0
                     gsign = -1
                     fc(1) = -1
                     fc(2) = +1
                     fc(3) = fsign(3)
                     fc(4) = fsign(4)
c
                  elseif( (lsymmcontrib .and. jd.le.14) .or. jd.le.8) then !qg->qq"q"~H , g from (26)
                     do mu =0,3
                        qbar(mu) = p(mu,2,nd) !i.s. gluon
                        kbar(mu,phystodiag(1))=p(mu,1,nd)
                        kbar(mu,phystodiag(3))=p(mu,3,nd)
                        kbar(mu,3) = p(mu,4,nd)!4
                        kbar(mu,4) = p(mu,5,nd)!5
                     enddo
                     qbar(4) = 0.0d0
                     gsign = -1
                     fc(1) = fsign(1)
                     fc(2) = fsign(2)
                     fc(3) = fsign(5)
                     fc(4) = fsign(6)                        
                  elseif(lsymmcontrib .and. jd.le.16) then !qg->qq"q"~H , g from (24)
                     do mu =0,3
                        qbar(mu) = p(mu,2,nd) !i.s. gluon
                        kbar(mu,phystodiag(1))=p(mu,1,nd)
                        kbar(mu,phystodiag(3))=p(mu,3,nd)
                        kbar(mu,phystodiag(4))=p(mu,4,nd)!4
                        kbar(mu,phystodiag(2))=p(mu,5,nd)!5
                     enddo
                     qbar(4) = 0.0d0
                     gsign = -1
                     fc(1) = fsign(1)
                     fc(2) = fsign(2)
                     fc(3) = -1
                     fc(4) = +1
                  else
                     write(*,*) "Error in amplitudes/hjjj/qqh4q_sub.F"
                     write(*,*) "parameters:", lsymmcontrib, jd
                     stop
                  endif
c     
                  if(ldebug1) then
                     print*,"tilde-momentum", jd, ": mapping", nd
                     do i = 1,n_p-1
                        print*,"p(mu,i,nd)=",(p(mu,i,nd),mu=0,3)
                     enddo
                     do i = 1,n_p-2
                        print*,"ktilde(mu,i)=",(kbar(mu,phystodiag(i)),mu=0,3)
                     enddo
                     print*,"ktilde(mu,i)=",(qbar(mu),mu=0,3)
                  endif

                  call qqhqqj_c_virt(0,kbar,fc,qbar,gsign,nd,
     $                 uucc_c,uuss_c,ddcc_c,
     $                 ddss_c,udsc_c,ducs_c)
c                          
c     j = 1 is gluon emission of 21 line of born
c     j = 2 is gluon emission of 43 line of born
c          
                  do j = 2,3
                     NCmatrix_b(0,0,j-1,nd)=uucc_c(j) !|M_b|^2
                     NCmatrix_b(0,1,j-1,nd)=uuss_c(j)
                     NCmatrix_b(1,0,j-1,nd)=ddcc_c(j)
                     NCmatrix_b(1,1,j-1,nd)=ddss_c(j)
                     CCmatrix_b(0,j-1,nd)=udsc_c(j)
                     CCmatrix_b(1,j-1,nd)=ducs_c(j)
                  enddo

c     for spin correlations                  
                  call qqHqqj_spcor(qvec(0,nd),kbar,fc,
     $                 qbar,gsign,nd,uucc_c,uuss_c,
     $                 ddcc_c,ddss_c,udsc_c,ducs_c)
                  do j = 2,3
                     NCmatrix_v(0,0,j-1,nd)=uucc_c(j) !|M_v|^2
                     NCmatrix_v(0,1,j-1,nd)=uuss_c(j)
                     NCmatrix_v(1,0,j-1,nd)=ddcc_c(j)
                     NCmatrix_v(1,1,j-1,nd)=ddss_c(j)
                     CCmatrix_v(0,j-1,nd)=udsc_c(j)
                     CCmatrix_v(1,j-1,nd)=ducs_c(j)
                  enddo

               endif
            enddo               !dipole loop
         endif !nlo
#endif /* WITH_NLO */

! NC case
         do if1=1,nflVBF    !(nfl/2)*2
            iflav(diagToPhys(1),:)=if1*fsign(1)
            iflav(diagToPhys(2),:)=if1*fsign(2)
            do if2=1,nflVBF    !(nfl/2)*2
               iflav(diagToPhys(3),:)=if2*fsign(3)
               iflav(diagToPhys(4),:)=if2*fsign(4)
               do if3=1,nflVBF    !(nfl/2)*2
                  iflav(diagToPhys(5),:) = if3*fsign(5)
                  iflav(diagToPhys(6),:) = if3*fsign(6)

                  ! symmetry factor
                  symmfac = 1
                  interference = 0
                  if (linterference .and. if1.eq.if3) then
                    symmfac(1)=1/2d0
                    interference(1) = 1
                  endif
                  if (linterference .and. if2.eq.if3) then
                    symmfac(2)=1/2d0
                    interference(2) = 1
                  endif
                  if (.not.lsymmcontrib) symmfac=1

                  if(nlo.eq.1) then
#ifdef WITH_NLO
                     k=FL_H4j(iflav,diagToPhys,2)
                     if (lok(1)) then
                       res(k,1)=pdf(sign1*iflav(1,1),1)*
     1                       pdf(sign2*iflav(2,1),2)*
     2                       ( NCmatrix_r(mod(if1,2),mod(if2,2),mod(if3,2),1,interference(1))*symmfac(1)
     2                        +NCmatrix_r(mod(if1,2),mod(if2,2),mod(if3,2),2,interference(2))*symmfac(2))*
     3                       polcolqq
                     else
                       res(k,1) = 0.0d0
                     endif
                     res(k,0) = res(k,1)
#endif /* WITH_NLO */
                  else             !color flow 
                     do icol = 1,2
                       k=FL_H4j(iflav,diagToPhys,icol)
                       if (lok(1)) then
                         res(k,1)=pdf(sign1*iflav(1,1),1)*
     1                       pdf(sign2*iflav(2,1),2)*
     2                       NCmatrix_r(mod(if1,2),mod(if2,2),mod(if3,2),icol,interference(icol))*
     2                       symmfac(icol)*
     3                       polcolqq
                       else
                         res(k,1) = 0.0d0
                       endif
                       res(k,0) = res(k,1)

                     enddo
                  endif !nlo
                  
                  do jd = 1,ndmax
                     nd = ndmap(jd) !dipole map
                     if (lok(nd).and.(nlo.eq.1)) then
#ifdef WITH_NLO
                        if( (lsymmcontrib .and. jd.le.8) .or. jd.le.4) then
                           if1l = if1
                           if2l = if2
                        elseif( (lsymmcontrib .and. jd.le.12) .or. jd.le.6) then
                           if1l = if3
                           if2l = if2
                        elseif( (lsymmcontrib .and. jd.le.16) .or. jd.le.8) then
                           if1l = if1 !1
                           if2l = if3 !3
                        endif
                        msq_b(1,nd) = NCmatrix_b(mod(if1l,2),
     $                       mod(if2l,2),1,nd)
                        msq_b(2,nd) = NCmatrix_b(mod(if1l,2),
     $                       mod(if2l,2),2,nd)
                        msq_v(1,nd) = NCmatrix_v(mod(if1l,2),
     $                       mod(if2l,2),1,nd)
                        msq_v(2,nd) = NCmatrix_v(mod(if1l,2),
     $                       mod(if2l,2),2,nd)
c     
c     compute dipole subtraction terms                
c     
                        dip = 0.0d0
c
c additional splittings only for same-flavour cases
                        if(nd.eq.4)
     $                       dip = TF/2.0d0*(msq_b(1,4)*sub(gq,4) - 
     $                       msq_v(1,4)*subv(4))*polcolqq*
     $                       symmfac(1)
                        if(nd.eq.14)
     $                       dip = TF/2.0d0*(msq_b(2,14)*sub(gq,14) - 
     $                       msq_v(2,14)*subv(14))*polcolqq*
     $                       symmfac(2)
                        if(nd.eq.7)
     $                       dip = TF/2.0d0*(msq_b(1,7)*sub(gq,7) - 
     $                       msq_v(1,7)*subv(7))*polcolqq*
     $                       symmfac(1)
                        if(nd.eq.17)
     $                       dip =  TF/2.0d0*(msq_b(2,17)*sub(gq,17)- 
     $                       msq_v(2,17)*subv(17))*polcolqq*
     $                       symmfac(2)
                        if(nd.eq.23) 
     $                       dip = CF/2.d0*(msq_b(1,23)*sub(gq,23) + 
     $                       msq_v(1,23)*subv(23))*polcolgq*
     $                       symmfac(1)
                        if(nd.eq.24) 
     $                       dip = CF/2.0d0*(msq_b(1,24)*sub(gq,24) + 
     $                       msq_v(1,24)*subv(24))*polcolgq*
     $                       symmfac(1)
                        if(nd.eq.25) 
     $                       dip = CF/2.0d0*(msq_b(2,25)*sub(gq,25) + 
     $                       msq_v(2,25)*subv(25))*polcolgq*
     $                       symmfac(2)
                        if(nd.eq.22) 
     $                       dip = CF/2.0d0*(msq_b(2,22)*sub(gq,22) + 
     $                       msq_v(2,22)*subv(22))*polcolgq*
     $                       symmfac(2)
                        if (if1.eq.if3) then
                          ! fsign(1).gt.0
                          if(nd.eq.3)
     $                         dip = TF/2.0d0*(msq_b(1,3)*sub(gq,3) - 
     $                         msq_v(1,3)*subv(3))*polcolqq*
     $                         symmfac(1)
                          if(nd.eq.6)
     $                         dip = TF/2.0d0*(msq_b(1,6)*sub(gq,6) - 
     $                         msq_v(1,6)*subv(6))*polcolqq*
     $                         symmfac(1)
                          if(nd.eq.8) 
     $                         dip = CF/2.0d0*(msq_b(1,8)*sub(gq,8) + 
     $                         msq_v(1,8)*subv(8))*polcolgq*
     $                         symmfac(1)
                          if(nd.eq.10) 
     $                         dip = CF/2.0d0*(msq_b(1,10)*sub(gq,10) + 
     $                         msq_v(1,10)*subv(10))*polcolgq*
     $                         symmfac(1)
                          ! fsign(1).lt.0
                          if(nd.eq.2)
     $                         dip = TF/2.0d0*(msq_b(1,2)*sub(gq,2) - 
     $                         msq_v(1,2)*subv(2))*polcolqq*
     $                         symmfac(1)
                          if(nd.eq.5)
     $                         dip = TF/2.0d0*(msq_b(1,5)*sub(gq,5) - 
     $                         msq_v(1,5)*subv(5))*polcolqq*
     $                         symmfac(1)
                          if(nd.eq.9) 
     $                         dip = CF/2.0d0*(msq_b(1,9)*sub(gq,9) + 
     $                         msq_v(1,9)*subv(9))*polcolgq*
     $                         symmfac(1)
                          if(nd.eq.11) 
     $                         dip = CF/2.0d0*(msq_b(1,11)*sub(gq,11) + 
     $                         msq_v(1,11)*subv(11))*polcolgq*
     $                         symmfac(1)
                        endif
                        if (if2.eq.if3) then
                          ! fsign(3).gt.0
                          if(nd.eq.13)
     $                         dip = TF/2.0d0*(msq_b(2,13)*sub(gq,13) - 
     $                         msq_v(2,13)*subv(13))*polcolqq*
     $                         symmfac(2)
                          if(nd.eq.16)
     $                         dip =  TF/2.0d0*(msq_b(2,16)*sub(gq,16)- 
     $                         msq_v(2,16)*subv(16))*polcolqq*
     $                         symmfac(2)
                          if(nd.eq.18) 
     $                         dip = CF/2.0d0*(msq_b(2,18)*sub(gq,18) + 
     $                         msq_v(2,18)*subv(18))*polcolgq*
     $                         symmfac(2)
                          if(nd.eq.20) 
     $                         dip = CF/2.0d0*(msq_b(2,20)*sub(gq,20) + 
     $                         msq_v(2,20)*subv(20))*polcolgq*
     $                         symmfac(2)
                          ! fsign(3).lt.0
                          if(nd.eq.12)
     $                         dip = TF/2.0d0*(msq_b(2,12)*sub(gq,12) - 
     $                         msq_v(2,12)*subv(12))*polcolqq*
     $                         symmfac(2)
                          if(nd.eq.15)
     $                         dip =  TF/2.0d0*(msq_b(2,15)*sub(gq,15)- 
     $                         msq_v(2,15)*subv(15))*polcolqq*
     $                         symmfac(2)
                          if(nd.eq.19) 
     $                         dip = CF/2.0d0*(msq_b(2,19)*sub(gq,19) + 
     $                         msq_v(2,19)*subv(19))*polcolgq*
     $                         symmfac(2)
                          if(nd.eq.21) 
     $                         dip = CF/2.0d0*(msq_b(2,21)*sub(gq,21) + 
     $                         msq_v(2,21)*subv(21))*polcolgq*
     $                         symmfac(2)
                        endif
c                        
                        res(k,nd) = -pdf(sign1*iflav(1,1),1)*
     $                       pdf(sign2*iflav(2,1),2)*dip
     
                        res(k,0) = res(k,0)+res(k,nd)
#endif /* WITH_NLO */
                     else
                        res(k,nd) = 0.0d0
                     endif
                  enddo         !loop over dipoles
                  
                  if(ldebug) then
c     take the limit of pi.pj -> 0               
                     delta = abs(dot0p(p(0,5,1),p(0,6,1)))
                     deltapair=56
                     newdelta= abs(dot0p(p(0,1,1),p(0,3,1)))
                     if (newdelta.lt.delta) then
                       delta=newdelta
                       deltapair=13
                     endif
                     newdelta= abs(dot0p(p(0,2,1),p(0,4,1)))
                     if (newdelta.lt.delta) then
                       delta=newdelta
                       deltapair=24
                     endif
                     if ( if1 .eq. if3 ) then
                       if ( fsign(1) .gt. 0 ) then
                         newdelta= abs(dot0p(p(0,1,1),p(0,6,1)))
                         if (newdelta.lt.delta) then
                           delta=newdelta
                           deltapair=16
                         endif
                         newdelta= abs(dot0p(p(0,3,1),p(0,5,1)))
                         if (newdelta.lt.delta) then
                           delta=newdelta
                           deltapair=35
                         endif
                       else
                         newdelta= abs(dot0p(p(0,1,1),p(0,5,1)))
                         if (newdelta.lt.delta) then
                           delta=newdelta
                           deltapair=15
                         endif
                         newdelta= abs(dot0p(p(0,3,1),p(0,6,1)))
                         if (newdelta.lt.delta) then
                           delta=newdelta
                           deltapair=36
                         endif
                       endif
                     endif
                     if ( if2 .eq. if3 ) then
                       if ( fsign(3) .gt. 0 ) then
                         newdelta= abs(dot0p(p(0,2,1),p(0,6,1)))
                         if (newdelta.lt.delta) then
                           delta=newdelta
                           deltapair=26
                         endif
                         newdelta= abs(dot0p(p(0,4,1),p(0,5,1)))
                         if (newdelta.lt.delta) then
                           delta=newdelta
                           deltapair=45
                         endif
                       else
                         newdelta= abs(dot0p(p(0,2,1),p(0,5,1)))
                         if (newdelta.lt.delta) then
                           delta=newdelta
                           deltapair=25
                         endif
                         newdelta= abs(dot0p(p(0,4,1),p(0,6,1)))
                         if (newdelta.lt.delta) then
                           delta=newdelta
                           deltapair=46
                         endif
                       endif
                     endif
                     if(abs(delta).lt.1d0 .and. lok(1))then
                        print*,"NC" 
                        print*,"momentum"
                        do i = 1,n_p
                        print*,"p(mu,",i,")=",(p(mu,i,1),mu=0,3)
                        enddo
                        print*,"p(mu,H)=",(pbar(mu,8),mu=0,3)
                        write(*,*) "fsign:", fsign(1:6)
                        sum=0.0d0
                        do jd = 1,ndmax
                           nd = ndmap(jd)
                           if(lok(nd)) then 
                              sum = sum + res(k,nd)
                              if(res(k,1).ne.0d0) then
                                 print*,"D(",ip(nd),kp(nd),jp(nd),
     -                                ")/|m2s|",nd,
     $                                res(k,nd)/res(k,1),lok(nd)
                              else
                                 print*,"D(",ip(nd),kp(nd),jp(nd),")",nd
     $                                ,res(k,nd),lok(nd)
                              endif
                           endif
                        enddo
                        print*,"f1,f2,f3=",if1,if2,if3
                        print*,"m2real=",res(k,1)
                        print*,"dip_sum=",sum
                        print*,"delta=(",deltapair,")",delta
                        if(res(k,1).ne.0d0) then
                          print*,"dip_sum/res(k,1)",abs(sum/res(k,1))
                        endif
                        print*,"-----------------"     
                     endif  
                  endif         !debug

               enddo            !if3
            enddo               !if2
         enddo                  !if1

         ! Please note:
         ! initial-state bbar might get generated from if.=6 entries.
         ! We use the full list including top here and sort out all
         ! forbidden combinations later.
         do if1=1,6
            do if2=1,6
               do if3=1,6
                  iflav(diagToPhys(1),:)=if1*fsign(1)
                  iflav(diagToPhys(2),:)=if1*fsign(2)
                  iflav(diagToPhys(3),:)=if2*fsign(3)
                  iflav(diagToPhys(4),:)=if2*fsign(4)
                  iflav(diagToPhys(5),:)=if3*fsign(5)
                  iflav(diagToPhys(6),:)=if3*fsign(6)

                  ! get flavours for CC case
                  if (mod(if1,2).eq.mod(if2,2)) then
                    if (fsign(1).eq.-1) then
                      iflav(diagToPhys(1),1) = (if1+2*mod(if1,2)-1)*fsign(1)
                    endif
                    iflav(diagToPhys(2),1) = iflav(diagToPhys(1),1)
                    if (fsign(3).eq.-1) then
                      iflav(diagToPhys(3),2) = (if2+2*mod(if2,2)-1)*fsign(3)
                    endif
                    iflav(diagToPhys(2),1) = iflav(diagToPhys(1),1)
                    iflav(diagToPhys(2),2) = (if1+2*mod(if1,2)-1)*fsign(2)
                    iflav(diagToPhys(4),1) = (if2+2*mod(if2,2)-1)*fsign(4)
                    iflav(diagToPhys(4),2) = iflav(diagToPhys(3),2)
                    iflav(diagToPhys(6),:) = (if3+2*mod(if3,2)-1)*fsign(6)
                  else
                    iflav(diagToPhys(2),:) = (if1+2*mod(if1,2)-1)*fsign(2)
                    iflav(diagToPhys(4),:) = (if2+2*mod(if2,2)-1)*fsign(4)
                    if (fsign(1).eq.1) then
                      iflav(diagToPhys(2),3) = (if3+2*mod(if3,2)-1)*fsign(2)
                      iflav(diagToPhys(6),3) = iflav(diagToPhys(1),3)
                    elseif (fsign(1).eq.-1) then
                      iflav(diagToPhys(1),3) = (if3+2*mod(if3,2)-1)*fsign(1)
                      iflav(diagToPhys(5),3) = iflav(diagToPhys(2),3)
                    endif
                    if (fsign(3).eq.1) then
                      iflav(diagToPhys(4),4) = (if3+2*mod(if3,2)-1)*fsign(4)
                      iflav(diagToPhys(6),4) = iflav(diagToPhys(3),4)
                    elseif (fsign(3).eq.-1) then
                      iflav(diagToPhys(3),4) = (if3+2*mod(if3,2)-1)*fsign(3)
                      iflav(diagToPhys(5),4) = iflav(diagToPhys(4),4)
                    endif
                  endif

                  ! symmetry factor
                  symmfac = 1
                  if (lsymmcontrib .and. linterference) then
                    if ( iflav(3,1).eq.iflav(5,1) .or.
     &                   iflav(3,1).eq.iflav(6,1) ) then
                      symmfac(1)=1/2d0
                    endif
                    if ( iflav(4,2).eq.iflav(5,2) .or.
     &                   iflav(4,2).eq.iflav(6,2) ) then
                      symmfac(2)=1/2d0
                    endif
                    ! icol 3 and 4 cannot have identical particles
                  endif

                  ! interference
                  interference = 0
                  if (linterference .and. (if1-1)/2.eq.(if3-1)/2) interference(1)=1
                  if (linterference .and. (if2-1)/2.eq.(if3-1)/2) interference(2)=1
                  if ( interference(1) .eq. 1 .or.
     &                 fsign(1).eq.+1 .and. .not.
     &                   ( mod(if1,2).ne.mod(if2,2) .and. mod(if1,2).eq.mod(if3,2) ) .or.
     &                 fsign(1).eq.-1 .and. .not.
     &                   ( mod(if1,2).ne.mod(if2,2) .and.  mod(if1,2).ne.mod(if3,2) )
     &               ) symmfac(3) = 0
                  if ( interference(2) .eq. 1 .or.
     &                 fsign(3).eq.+1 .and. .not.
     &                   ( mod(if1,2).ne.mod(if2,2) .and. mod(if2,2).eq.mod(if3,2) ) .or.
     &                 fsign(3).eq.-1 .and. .not.
     &                   ( mod(if1,2).ne.mod(if2,2) .and.  mod(if2,2).ne.mod(if3,2) )
     &               ) symmfac(4) = 0

                  ! invalid combination
                  if ( mod(if1,2).eq.mod(if3,2) .and.
     &                 mod(if2,2).eq.mod(if3,2)
     &               ) symmfac = 0

                  ! check for anything larger than nflvbf
                  do i=1,4
                    if ( any(abs(iflav(:,i)).gt.nflvbf) ) then
                      symmfac(i) = 0
                    endif
                  enddo

                  ! cycle if no contribution
                  if (symmfac(1).eq.0 .and. symmfac(3).eq.0 .and.
     &                symmfac(2).eq.0 .and. symmfac(4).eq.0) cycle

                  if(nlo.eq.1) then
#ifdef WITH_NLO
                     k=FL_H4j(iflav,diagToPhys,2)
                     if (lok(1)) then
                       res(k,1)=pdf(sign1*iflav(1,1),1)*
     1                       pdf(sign2*iflav(2,1),2)*
     2                       ( CCmatrix_r(mod(if1,2),mod(if2,2),mod(if3,2),1,interference(1))*symmfac(1)
     2                        +CCmatrix_r(mod(if1,2),mod(if2,2),mod(if3,2),2,interference(2))*symmfac(2)
     2                        +CCmatrix_r(mod(if1,2),mod(if2,2),mod(if3,2),3,interference(3))*symmfac(3)
     2                        +CCmatrix_r(mod(if1,2),mod(if2,2),mod(if3,2),4,interference(4))*symmfac(4))*
     3                       polcolqq
                     else
                       res(k,1) = 0.0d0
                     endif
                     res(k,0) = res(k,1)
#endif /* WITH_NLO */
                  else             !color flow 
                     do icol = 1,4
                       k=FL_H4j(iflav(1,icol),diagToPhys,icol)
                       if (lok(1)) then
                          res(k,1)=pdf(sign1*iflav(1,icol),1)*
     1                         pdf(sign2*iflav(2,icol),2)*
     2                         CCmatrix_r(mod(if1,2),mod(if2,2),mod(if3,2),icol,interference(icol))*
     2                         symmfac(icol)*
     3                         polcolqq
                       else
                          res(k,1) = 0.0d0
                       endif
                       res(k,0) = res(k,1)
                     enddo
                  endif !nlo
                  
                  do jd = 1,ndmax
                     nd = ndmap(jd) !dipole map
                     if (lok(nd).and.(nlo.eq.1)) then
#ifdef WITH_NLO
! calculate flavour as flavour opposite to other line
                       if1l = mod(abs(iflav(diagtophys(3),2))+1,2)
                       if2l = mod(abs(iflav(diagtophys(1),1))+1,2)
                       msq_b(1,nd) = CCmatrix_b(if1l,1,nd)
                       msq_b(2,nd) = CCmatrix_b(if2l,2,nd)
                       msq_v(1,nd) = CCmatrix_v(if1l,1,nd)
                       msq_v(2,nd) = CCmatrix_v(if2l,2,nd)
c     
c     compute dipole subtraction terms                
c     
                        dip = 0.0d0
c
c check if splitting is an allowed combination
                        if (iflav(5,1) .eq. -iflav(6,1)) then
                          ! iflav(5,2) .eq. -iflav(6,2) is same condition
                          if(nd.eq.4)
     $                         dip = TF/2.0d0*(msq_b(1,4)*sub(gq,4) -
     $                         msq_v(1,4)*subv(4))*polcolqq*
     $                         symmfac(1)
                          if(nd.eq.14)
     $                         dip = TF/2.0d0*(msq_b(2,14)*sub(gq,14) -
     $                         msq_v(2,14)*subv(14))*polcolqq*
     $                         symmfac(2)
                          if(nd.eq.7)
     $                         dip = TF/2.0d0*(msq_b(1,7)*sub(gq,7) -
     $                         msq_v(1,7)*subv(7))*polcolqq*
     $                         symmfac(1)
                          if(nd.eq.17)
     $                         dip = TF/2.0d0*(msq_b(2,17)*sub(gq,17) -
     $                         msq_v(2,17)*subv(17))*polcolqq*
     $                         symmfac(2)
                        endif
                        if (iflav(1,1) .eq. iflav(3,1)) then
                          if(nd.eq.23)
     $                         dip = CF/2.d0*(msq_b(1,23)*sub(gq,23) +
     $                         msq_v(1,23)*subv(23))*polcolgq*
     $                         symmfac(1)
                          if(nd.eq.24)
     $                         dip = CF/2.0d0*(msq_b(1,24)*sub(gq,24) +
     $                         msq_v(1,24)*subv(24))*polcolgq*
     $                         symmfac(1)
                        endif
                        if (iflav(2,2) .eq. iflav(4,2)) then
                          if(nd.eq.25)
     $                         dip = CF/2.0d0*(msq_b(2,25)*sub(gq,25) +
     $                         msq_v(2,25)*subv(25))*polcolgq*
     $                         symmfac(2)
                          if(nd.eq.22)
     $                         dip = CF/2.0d0*(msq_b(2,22)*sub(gq,22) +
     $                         msq_v(2,22)*subv(22))*polcolgq*
     $                         symmfac(2)
                        endif
                        ! fsign(1).gt.0
                        if ( (iflav(3,1) .eq. -iflav(5,1)) .or.
     &                       (symmfac(3) .ne. 0) .and.
     &                         (iflav(3,3) .eq. -iflav(5,3)) ) then
                          if(nd.eq.3)
     $                         dip = TF/2.0d0*(msq_b(1,3)*sub(gq,3) -
     $                         msq_v(1,3)*subv(3))*polcolqq*
     $                         symmfac(1)
                          if(nd.eq.6)
     $                         dip = TF/2.0d0*(msq_b(1,6)*sub(gq,6) -
     $                         msq_v(1,6)*subv(6))*polcolqq*
     $                         symmfac(1)
                        endif
                        if ( (iflav(1,1) .eq. iflav(6,1)) .or.
     &                       (symmfac(3) .ne. 0) .and.
     &                         (iflav(1,3) .eq. iflav(6,3)) ) then
                          if(nd.eq.8)
     $                         dip = CF/2.0d0*(msq_b(1,8)*sub(gq,8) +
     $                         msq_v(1,8)*subv(8))*polcolgq*
     $                         symmfac(1)
                          if(nd.eq.10)
     $                         dip = CF/2.0d0*(msq_b(1,10)*sub(gq,10) +
     $                         msq_v(1,10)*subv(10))*polcolgq*
     $                         symmfac(1)
                        endif
                        ! fsign(1).lt.0
                        if ( (iflav(3,1) .eq. -iflav(6,1)) .or.
     &                       (symmfac(3) .ne. 0) .and.
     &                         (iflav(3,3) .eq. -iflav(6,3)) ) then
                          if(nd.eq.2)
     $                         dip = TF/2.0d0*(msq_b(1,2)*sub(gq,2) -
     $                         msq_v(1,2)*subv(2))*polcolqq*
     $                         symmfac(1)
                          if(nd.eq.5)
     $                         dip = TF/2.0d0*(msq_b(1,5)*sub(gq,5) -
     $                         msq_v(1,5)*subv(5))*polcolqq*
     $                         symmfac(1)
                        endif
                        if ( (iflav(1,1) .eq. iflav(5,1)) .or.
     &                       (symmfac(3) .ne. 0) .and.
     &                         (iflav(1,3) .eq. iflav(5,3)) ) then
                          if(nd.eq.9)
     $                         dip = CF/2.0d0*(msq_b(1,9)*sub(gq,9) +
     $                         msq_v(1,9)*subv(9))*polcolgq*
     $                         symmfac(1)
                          if(nd.eq.11)
     $                         dip = CF/2.0d0*(msq_b(1,11)*sub(gq,11) +
     $                         msq_v(1,11)*subv(11))*polcolgq*
     $                         symmfac(1)
                        endif
                        ! fsign(3).gt.0
                        if ( (iflav(4,2) .eq. -iflav(5,2)) .or.
     &                       (symmfac(4) .ne. 0) .and.
     &                         (iflav(4,4) .eq. -iflav(5,4)) ) then
                          if(nd.eq.13)
     $                         dip = TF/2.0d0*(msq_b(2,13)*sub(gq,13) -
     $                         msq_v(2,13)*subv(13))*polcolqq*
     $                         symmfac(2)
                          if(nd.eq.16)
     $                         dip = TF/2.0d0*(msq_b(2,16)*sub(gq,16) -
     $                         msq_v(2,16)*subv(16))*polcolqq*
     $                         symmfac(2)
                        endif
                        if ( (iflav(2,2) .eq. iflav(6,2)) .or.
     &                       (symmfac(4) .ne. 0) .and.
     &                         (iflav(2,4) .eq. iflav(6,4)) ) then
                          if(nd.eq.18)
     $                         dip = CF/2.0d0*(msq_b(2,18)*sub(gq,18) +
     $                         msq_v(2,18)*subv(18))*polcolgq*
     $                         symmfac(2)
                          if(nd.eq.20)
     $                         dip = CF/2.0d0*(msq_b(2,20)*sub(gq,20) +
     $                         msq_v(2,20)*subv(20))*polcolgq*
     $                         symmfac(2)
                        endif
                        ! fsign(3).lt.0
                        if ( (iflav(4,2) .eq. -iflav(6,2)) .or.
     &                       (symmfac(4) .ne. 0) .and.
     &                         (iflav(4,4) .eq. -iflav(6,4)) ) then
                          if(nd.eq.12)
     $                         dip = TF/2.0d0*(msq_b(2,12)*sub(gq,12) -
     $                         msq_v(2,12)*subv(12))*polcolqq*
     $                         symmfac(2)
                          if(nd.eq.15)
     $                         dip = TF/2.0d0*(msq_b(2,15)*sub(gq,15) -
     $                         msq_v(2,15)*subv(15))*polcolqq*
     $                         symmfac(2)
                        endif
                        if ( (iflav(2,2) .eq. iflav(5,2)) .or.
     &                       (symmfac(4) .ne. 0) .and.
     &                         (iflav(2,4) .eq. iflav(5,4)) ) then
                          if(nd.eq.19)
     $                         dip = CF/2.0d0*(msq_b(2,19)*sub(gq,19) +
     $                         msq_v(2,19)*subv(19))*polcolgq*
     $                         symmfac(2)
                          if(nd.eq.21)
     $                         dip = CF/2.0d0*(msq_b(2,21)*sub(gq,21) +
     $                         msq_v(2,21)*subv(21))*polcolgq*
     $                         symmfac(2)
                        endif
c                        
                        res(k,nd) = -pdf(sign1*iflav(1,1),1)*
     $                       pdf(sign2*iflav(2,1),2)*dip
     
                        res(k,0) = res(k,0)+res(k,nd)
#endif /* WITH_NLO */
                     else
                        res(k,nd) = 0.0d0
                     endif
                  enddo         !loop over dipoles
                  
                  if(ldebug) then
c     take the limit of pi.pj -> 0               
                     delta = 1d99
                     if (iflav(5,1) .eq. -iflav(6,1)) then
                       newdelta = abs(dot0p(p(0,5,1),p(0,6,1)))
                       if (newdelta.lt.delta) then
                         delta=newdelta
                         deltapair=56
                       endif
                     endif
                     if (iflav(1,1) .eq. iflav(3,1)) then
                       newdelta= abs(dot0p(p(0,1,1),p(0,3,1)))
                       if (newdelta.lt.delta) then
                         delta=newdelta
                         deltapair=13
                       endif
                     endif
                     if (iflav(2,2) .eq. iflav(4,2)) then
                       newdelta= abs(dot0p(p(0,2,1),p(0,4,1)))
                       if (newdelta.lt.delta) then
                         delta=newdelta
                         deltapair=24
                       endif
                     endif
                     if ( (iflav(3,1) .eq. -iflav(5,1)) .or.
     &                    (symmfac(3) .ne. 0) .and.
     &                      (iflav(3,3) .eq. -iflav(5,3)) ) then
                       newdelta= abs(dot0p(p(0,3,1),p(0,5,1)))
                       if (newdelta.lt.delta) then
                         delta=newdelta
                         deltapair=35
                       endif
                     endif
                     if ( (iflav(1,1) .eq. iflav(6,1)) .or.
     &                    (symmfac(3) .ne. 0) .and.
     &                      (iflav(1,3) .eq. iflav(6,3)) ) then
                       newdelta= abs(dot0p(p(0,1,1),p(0,6,1)))
                       if (newdelta.lt.delta) then
                         delta=newdelta
                         deltapair=16
                       endif
                     endif
                     if ( (iflav(3,1) .eq. -iflav(6,1)) .or.
     &                    (symmfac(3) .ne. 0) .and.
     &                      (iflav(3,3) .eq. -iflav(6,3)) ) then
                       newdelta= abs(dot0p(p(0,3,1),p(0,6,1)))
                       if (newdelta.lt.delta) then
                         delta=newdelta
                         deltapair=36
                       endif
                     endif
                     if ( (iflav(1,1) .eq. iflav(5,1)) .or.
     &                    (symmfac(3) .ne. 0) .and.
     &                      (iflav(1,3) .eq. iflav(5,3)) ) then
                       newdelta= abs(dot0p(p(0,1,1),p(0,5,1)))
                       if (newdelta.lt.delta) then
                         delta=newdelta
                         deltapair=15
                       endif
                     endif
                     if ( (iflav(4,2) .eq. -iflav(5,2)) .or.
     &                    (symmfac(4) .ne. 0) .and.
     &                      (iflav(4,4) .eq. -iflav(5,4)) ) then
                       newdelta= abs(dot0p(p(0,4,1),p(0,5,1)))
                       if (newdelta.lt.delta) then
                         delta=newdelta
                         deltapair=45
                       endif
                     endif
                     if ( (iflav(2,2) .eq. iflav(6,2)) .or.
     &                    (symmfac(4) .ne. 0) .and.
     &                      (iflav(2,4) .eq. iflav(6,4)) ) then
                       newdelta= abs(dot0p(p(0,2,1),p(0,6,1)))
                       if (newdelta.lt.delta) then
                         delta=newdelta
                         deltapair=26
                       endif
                     endif
                     if ( (iflav(4,2) .eq. -iflav(6,2)) .or.
     &                    (symmfac(4) .ne. 0) .and.
     &                      (iflav(4,4) .eq. -iflav(6,4)) ) then
                       newdelta= abs(dot0p(p(0,4,1),p(0,6,1)))
                       if (newdelta.lt.delta) then
                         delta=newdelta
                         deltapair=46
                       endif
                     endif
                     if ( (iflav(2,2) .eq. iflav(5,2)) .or.
     &                    (symmfac(4) .ne. 0) .and.
     &                      (iflav(2,4) .eq. iflav(5,4)) ) then
                       newdelta= abs(dot0p(p(0,2,1),p(0,5,1)))
                       if (newdelta.lt.delta) then
                         delta=newdelta
                         deltapair=25
                       endif
                     endif

                     if(abs(delta).lt.1d0 .and. lok(1))then
                        print*,"CC" 
                        print*,"momentum"
                        do i = 1,n_p
                        print*,"p(mu,",i,")=",(p(mu,i,1),mu=0,3)
                        enddo
                        print*,"p(mu,H)=",(pbar(mu,8),mu=0,3)
                        write(*,*) "fsign:", fsign(1:6)
                        write(*,*) "iflav 1", symmfac(1), iflav(1:6,1)
                        write(*,*) "iflav 2", symmfac(2), iflav(1:6,2)
                        write(*,*) "iflav 3", symmfac(3), iflav(1:6,3)
                        write(*,*) "iflav 4", symmfac(4), iflav(1:6,4)
                        sum=0.0d0
                        do jd = 1,ndmax
                           nd = ndmap(jd)
                           if(lok(nd)) then 
                              sum = sum + res(k,nd)
                              if(res(k,1).ne.0d0) then
                                 print*,"D(",ip(nd),kp(nd),jp(nd),
     -                                ")/|m2s|",nd,
     $                                res(k,nd)/res(k,1),lok(nd)
                              else
                                 print*,"D(",ip(nd),kp(nd),jp(nd),")",nd
     $                                ,res(k,nd),lok(nd)
                              endif
                           endif
                        enddo
                        print*,"f1,f2,f3=",if1,if2,if3
                        print*,"m2real=",res(k,1)
                        print*,"dip_sum=",sum
                        print*,"delta=(",deltapair,")",delta
                        if(res(k,1).ne.0d0) then
                          print*,"dip_sum/res(k,1)",abs(sum/res(k,1))
                        endif
                        print*,"-----------------"     
                     endif  
                  endif         !debug

               enddo            !if3
            enddo               !if2
         enddo                  !if1

         nmax = FL_H4j(iflav,diagToPhys,-2)
         return
         end
c******************************************************************
c
c   end subroutine wbf_h4j
c
c*****************************************************************
