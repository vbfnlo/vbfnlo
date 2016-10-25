c******************************************************************
c
c   begin subroutine wbf_h4q
c
c*****************************************************************
      subroutine wbf_h4qnc(xi,p,v,physToDiag,fsign,bos, 
     1     nlo,lok,xuz,
     2     pdf,res,nmin,nmax)
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
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
      double precision NCmatrix_r(0:1,0:1,0:1,0:2)
      double precision NCmatrix_b(0:1,0:1,2,2:max_kin)
      double precision NCmatrix_v(0:1,0:1,2,2:max_kin)
      double precision polcolqq,polcolgq
      double precision sum,delta
      integer iflav(6), diagToPhys(6),FL_H4j
      double precision dot0p
      external FL_H4j
      external dot0p
      logical ldebug, ldebug1,ldecay, lres
      parameter(ldebug=.false.,ldebug1=.false.)
      integer i,if1,if2,if3,j,k,mu,jd,ndmap(8)
      integer if1l,if2l,icol
      integer fc(4+max_v),gq,qq,gg,qg



      polcolqq = 1d0/(NC**2 * 4d0 * xi(1)*xi(2)) 
      polcolgq = 1d0/(NC*(NC**2-1d0) * 4d0 * xi(1)*xi(2))
c     set up ndmap
      ndmap(1)=4
      ndmap(2)=14
      ndmap(3)=7
      ndmap(4)=17
      ndmap(5)=23
      ndmap(6)=24
      ndmap(7)=25
      ndmap(8)=22
c     ave over polarizations and colors
c     symmetry factor for two gluons
      if(ldebug1) then
         print*,"momentum"
         do i = 1,n_p
            print*,"p(mu,i,1)=",(p(mu,i,1),mu=0,3)
         enddo
      endif
c     initialize matrix elements squared
      do i = 0,2
      uuccuu(i) = 0.0d0
      uussuu(i) = 0.0d0 
      ddccuu(i) = 0.0d0 
      ddssuu(i) = 0.0d0 
      uuccdd(i) = 0.0d0 
      uussdd(i) = 0.0d0 
      ddccdd(i) = 0.0d0 
      ddssdd(i) = 0.0d0 
      enddo
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
c
c              
      if (lok(1)) then
         call qqh4qnc_mg_c(pbar,fsign,
     $        uuccuu,uussuu,ddccuu,ddssuu,
     $        uuccdd,uussdd,ddccdd,ddssdd)
      endif
c     
      do i = 0,2
         NCmatrix_r(0,0,0,i) = uuccuu(i)
         NCmatrix_r(0,1,0,i) = uussuu(i)
         NCmatrix_r(1,0,0,i) = ddccuu(i)
         NCmatrix_r(1,1,0,i) = ddssuu(i)
         NCmatrix_r(0,0,1,i) = uuccdd(i)
         NCmatrix_r(0,1,1,i) = uussdd(i)
         NCmatrix_r(1,0,1,i) = ddccdd(i)
         NCmatrix_r(1,1,1,i) = ddssdd(i)
      enddo

c for the NLO case get the subtraction terms; 

#ifdef WITH_NLO
         if (nlo.eq.1) then
c     create dipole subtraction terms
            call create_ikj(ip,kp,jp)
            do jd = 1,8
               nd = ndmap(jd)
c            do nd = 2,n_kin
               call dips(nd,p,xuz,ip(nd),kp(nd),jp(nd),sub(1,nd),
     $              subv(nd),qvec(0,nd))
            enddo
c          
c     Compute born matrix elements with transformed momenta:
c     
c     
            do jd = 1,8
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
                  if(jd.le.4) then     !qQ -> qQgH born process
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
                  elseif(jd.le.6) then !gQ->Qq"q"~H
                     do mu =0,3
                        qbar(mu) = p(mu,1,nd) !i.s. gluon
                        kbar(mu,phystodiag(2)) = p(mu,2,nd)!i.s. quark
c                        
                        kbar(mu,phystodiag(4)) = p(mu,3,nd)!f.s. quark
c                       
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
                  elseif(jd.le.8) then !qg->qq"q"~H
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
                  else
                     stop
                  endif
c     
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
                  enddo

               endif
            enddo               !dipole loop
         endif !nlo
#endif /* WITH_NLO */

         do if1=1,nflVBF    !(nfl/2)*2
            iflav(1)=if1*fsign(physToDiag(1))
            iflav(3)=if1*fsign(physToDiag(3))
            do if2=1,nflVBF    !(nfl/2)*2
               iflav(2)=if2*fsign(physToDiag(2))
               iflav(4)=if2*fsign(physToDiag(4))
               do if3=1,nflVBF    !(nfl/2)*2
                  iflav(5) = if3*fsign(physToDiag(5))
                  iflav(6) = if3*fsign(physToDiag(6))
c
                  if(nlo.eq.1) then
#ifdef WITH_NLO
                     k=FL_H4j(iflav,diagToPhys,2)
                     if (lok(1)) then
                        res(k,1)=pdf(sign1*iflav(1),1)*
     1                       pdf(sign2*iflav(2),2)*
     2                       NCmatrix_r(mod(if1,2),mod(if2,2),
     $                       mod(if3,2),0)*
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
                           res(k,1)=pdf(sign1*iflav(1),1)*
     1                          pdf(sign2*iflav(2),2)*
     2                          NCmatrix_r(mod(if1,2),mod(if2,2),
     $                          mod(if3,2),icol)*
     3                          polcolqq
                        else
                           res(k,1) = 0.0d0
                        endif

                        res(k,0) = res(k,1)
                     enddo
                  endif !nlo
                  
                  do jd = 1,8
                     nd = ndmap(jd) !dipole map
                     if (lok(nd).and.(nlo.eq.1)) then
#ifdef WITH_NLO
                        if(jd.le.4) then
                           if1l = if1
                           if2l = if2
                        elseif(jd.le.6) then
                           if1l = if3
                           if2l = if2
                        elseif(jd.le.8) then !bug was is real emmission code
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
                        qq = 1  !labels for splittings 
                        qg = 2
                        gq = 3
                        gg = 4
c
                        dip = 0.0d0
c
                        if(nd.eq.4)
     $                       dip = TF/2.0d0*(msq_b(1,4)*sub(gq,4) - 
     $                       msq_v(1,4)*subv(4))*polcolqq 
                        if(nd.eq.14)
     $                       dip = TF/2.0d0*(msq_b(2,14)*sub(gq,14) - 
     $                       msq_v(2,14)*subv(14))*polcolqq 
                        if(nd.eq.7)
     $                       dip = TF/2.0d0*(msq_b(1,7)*sub(gq,7) - 
     $                       msq_v(1,7)*subv(7))*polcolqq
                        if(nd.eq.17)
     $                       dip =  TF/2.0d0*(msq_b(2,17)*sub(gq,17)- 
     $                       msq_v(2,17)*subv(17))*polcolqq 
                        if(nd.eq.23) 
     $                       dip = CF/2.d0*(msq_b(1,23)*sub(gq,23) + 
     $                       msq_v(1,23)*subv(23))*polcolgq
                        if(nd.eq.24) 
     $                       dip = CF/2.0d0*(msq_b(1,24)*sub(gq,24) + 
     $                       msq_v(1,24)*subv(24))*polcolgq
c
                        if(nd.eq.25) 
     $                       dip = CF/2.0d0*(msq_b(2,25)*sub(gq,25) + 
     $                       msq_v(2,25)*subv(25))*polcolgq
                        if(nd.eq.22) 
     $                       dip = CF/2.0d0*(msq_b(2,22)*sub(gq,22) + 
     $                       msq_v(2,22)*subv(22))*polcolgq

c                        
                        res(k,nd) = -pdf(sign1*iflav(1),1)*
     $                       pdf(sign2*iflav(2),2)*dip
     
                        res(k,0) = res(k,0)+res(k,nd)
#endif /* WITH_NLO */
                     else
                        res(k,nd) = 0.0d0
                     endif
                  enddo         !loop over dipoles
                  
                  if(ldebug) then
c     take the limit of pi.pj -> 0               
                     delta = dot0p(p(0,5,1),p(0,6,1))
                     if(abs(delta).lt.5d0 .and. lok(1))then
                        print*,"NC" 
                        sum=0.0d0
                        do jd = 1,8
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
                        print*,"delta=",delta
                        if(res(k,1).ne.0d0) 
     $                       print*,"dip_sum/res(k,1)",abs(sum/res(k,1))
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
