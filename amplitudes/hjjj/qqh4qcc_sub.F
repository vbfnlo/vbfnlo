c******************************************************************
c
c   begin subroutine wbf_h4q
c
c*****************************************************************
      subroutine wbf_h4qcc(xi,p,v,physToDiag,fsign,bos, 
     1     nlo,lok,xuz,
     2     pdf,res,nmin,nmax)
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "color.inc"
      double precision p(0:3,max_p,max_kin), v(0:3,max_v,max_kin), xi(nx), 
     $     xuz(2,2:max_kin)
      double precision pdf(-6:6,2,5), res(maxnumsubproc,0:max_kin)
      integer physToDiag(6), fsign(6+max_v), gsign, 
     $     bos, nlo, nmin, nmax
      logical lok(max_kin)
c
      double precision uucc_c(3),uuss_c(3),ddcc_c(3),ddss_c(3),
     $     udsc_c(3),ducs_c(3)
      double precision udsc(0:2),ducs(0:2)
      
c wbf_h4q calls the amplitude square routines 
c             qqh4qcc_mg    for qQ-->qQ qq~H     for bos=6
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
      integer jdmin,jdmax
      double precision pbar(0:3,6+max_v),qbar(0:4)
      double precision kbar(0:3,4+max_v)
      double precision sub(4,2:max_kin),subv(2:max_kin),qvec(0:3,2:max_kin)
      double precision dip,msq_b(0:2,2:max_kin),msq_v(0:2,2:max_kin)
      double precision CCmatrix_r(0:1,0:2)
      double precision CCmatrix_b(0:1,2,2:max_kin)
      double precision CCmatrix_v(0:1,2,2:max_kin)
      double precision polcolqq,polcolgq
      double precision sum,delta
      integer iflav(6), diagToPhys(6),FL_H4j
      double precision dot0p
      external FL_H4j
      external dot0p
      logical ldebug, ldebug1,ldecay, lres
      parameter(ldebug=.false.,ldebug1=.false.)
      integer i,if1,if2,if3,j,k,mu,jd,ndmap(8),icc1,icol
      logical chargedcurrent,oneAntiparticle,sametype
      integer fc(4+max_v),gq,qq,gg,qg



      polcolqq = 1.0d0/(NC**2 * 4.0d0 * xi(1)*xi(2)) 
      polcolgq = 1.0d0/(NC*(NC**2-1.0d0) * 4.0d0 * xi(1)*xi(2)) 
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
         udsc(i) = 0.0d0
         ducs(i) = 0.0d0
      enddo
            
      nmin = nmax+1

      do i = 1,6
         diagToPhys(physToDiag(i)) = i
      enddo
         
c get the real emission amplitude squared, 
      do mu = 0,3
         do i = 1,6
            pbar(mu,physToDiag(i))=p(mu,i,1)
         enddo         
      enddo
c
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
              
      if (lok(1)) 
     $    call qqh4qcc_mg_c(pbar,fsign,
     $     udsc,ducs)
 
      do i = 0,2
         CCmatrix_r(0,i) = udsc(i)
         CCmatrix_r(1,i) = ducs(i)
      enddo           

c for the NLO case get the subtraction terms;       
      if(fsign(5).eq.-fsign(6)) then !q5 and q6 are f.s.

#ifdef WITH_NLO
         if (nlo.eq.1) then
c     create dipole subtraction terms
            call create_ikj(ip,kp,jp)
            do jd = 1,8
               nd = ndmap(jd)
               call dips(nd,p,xuz,ip(nd),kp(nd),jp(nd),sub(1,nd),
     $              subv(nd),qvec(0,nd))
            enddo                  
c     
c     Compute born matrix elements with transformed momenta:
c     
c     
            do jd = 1,4
               nd = ndmap(jd) !mapping to dipoles set as data statement
               if (lok(nd).and.(nlo.eq.1)) then
C     qQ-> qQgH born process
                  do mu = 0,3          
                     do i = 1,4
                        kbar(mu,physToDiag(i))=p(mu,i,nd)
                      enddo
                     qbar(mu) = p(mu,5,nd) !gluon momenta
                  enddo
                  qbar(4) = 0.0d0
                  gsign = 1
                  do i=1,4
                     fc(i)=fsign(i)
                  enddo

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

                  call qqhqqj_c_virt(0,kbar,fc,qbar,gsign,nd,
     $                 uucc_c,uuss_c,ddcc_c,
     $                 ddss_c,udsc_c,ducs_c)
c     
c     j = 1 is gluon emission of 21 line of born
c     j = 2 is gluon emission of 43 line of born
c     
                  do j = 2,3
                     CCmatrix_b(0,j-1,nd)=udsc_c(j) !|M_b|^2
                     CCmatrix_b(1,j-1,nd)=ducs_c(j)
                  enddo
                              
c     for spin correlations                  
                  call qqHqqj_spcor(qvec(0,nd),kbar,fc,
     $                 qbar,gsign,nd,uucc_c,uuss_c,
     $                 ddcc_c,ddss_c,udsc_c,ducs_c)
                  do j = 2,3
                     CCmatrix_v(0,j-1,nd)=udsc_c(j) !|M_b|^2
                     CCmatrix_v(1,j-1,nd)=ducs_c(j)            
                  enddo
                     
               endif !lok(nd)
            enddo            !dipole loop
         endif !nlo.eq.1
#endif /* WITH_NLO */

         do if1=1,(nflVBF/2)*2    !(nfl/2)*2
            iflav(1)=if1*fsign(physToDiag(1))
            iflav(3)=if1*fsign(physToDiag(3))
            do if2=1,(nflVBF/2)*2    !(nfl/2)*2
               iflav(2)=if2*fsign(physToDiag(2))
               iflav(4)=if2*fsign(physToDiag(4))
               do if3=1,nflVBF    !(nfl/2)*2
               !particles 56 from gluon splitting, can be bb~
                  iflav(5) = if3*fsign(physToDiag(5))
                  iflav(6) = if3*fsign(physToDiag(6))
c     
C     Now check if there is a CC contribution for this choice of initial state
C     flavors (i,j). First: check whether initial uu or dd, not ud
                  sametype=(mod(if1,2)).eq.(mod(if2,2))
                  oneAntiparticle=
     &                 ((fsign(physToDiag(1))*fsign(physToDiag(2))).
     $                 eq.-1) 
c     true if only one particle is an antiparticle
                  ChargedCurrent=(oneAntiparticle.and.sametype) .or.
     &                 ( (.not.oneAntiparticle).and.(.not.sametype) )
                  if( ChargedCurrent .and. if1.le.4 .and. 
     $                 if2.le.4 ) then   
c     change 1<-->2 and 3<-->4 in outgoing quark flavors
                     iflav(3)=(if1+2*mod(if1,2)-1)*fsign(physToDiag(3))
                     iflav(4)=(if2+2*mod(if2,2)-1)*fsign(physToDiag(4))
                     icc1 = abs(iflav(diagtophys(1)))

                     if(nlo.eq.1) then
#ifdef WITH_NLO
                        k=FL_H4j(iflav,diagToPhys,2)
                        do nd = 0,n_kin
                           res(k,nd) = 0.0d0
                        enddo
                        if (lok(1)) then
                           res(k,1)=pdf(sign1*iflav(1),1,1)*
     $                          pdf(sign2*iflav(2),2,1)*
     $                          CCmatrix_r(mod(icc1,2),0)*
     $                          polcolqq
                        else
                           res(k,1)=0
                        endif
                           res(k,0) = res(k,1) 
#endif     
                     else    !color flow decomposition
                        do icol = 1,2
                           k=FL_H4j(iflav,diagToPhys,icol)
                           do nd = 0,n_kin
                              res(k,nd) = 0.0d0
                           enddo
                           if (lok(1)) then
                              res(k,1)=pdf(sign1*iflav(1),1,1)*
     $                             pdf(sign2*iflav(2),2,1)*
     $                             CCmatrix_r(mod(icc1,2),icol)*
     $                             polcolqq
                           else     
                              res(k,1)=0
                           endif      
                           res(k,0) = res(k,1)
                        enddo
                     endif   !nlo
                        
                     do jd = 1,4
                        nd = ndmap(jd) !dipole map

                        if (lok(nd).and.(nlo.eq.1)) then
                           msq_b(1,nd) = CCmatrix_b(mod(icc1,2),1,nd)
                           msq_b(2,nd) = CCmatrix_b(mod(icc1,2),2,nd)
                           msq_v(1,nd) = CCmatrix_v(mod(icc1,2),1,nd)
                           msq_v(2,nd) = CCmatrix_v(mod(icc1,2),2,nd)

c     compute dipole subtraction terms                     
                           qq = 1 !labels for splittings 
                           qg = 2
                           gq = 3
                           gg = 4

                           dip = 0.0d0
                           if(nd.eq.4)
     $                          dip = TF/2.0d0*(msq_b(1,4)*sub(gq,4) - 
     $                          msq_v(1,4)*subv(4))*polcolqq 
                           if(nd.eq.14)
     $                          dip = TF/2.0d0*(msq_b(2,14)*sub(gq,14) - 
     $                          msq_v(2,14)*subv(14))*polcolqq 
                           if(nd.eq.7)
     $                          dip = TF/2.0d0*(msq_b(1,7)*sub(gq,7) - 
     $                          msq_v(1,7)*subv(7))*polcolqq
                           if(nd.eq.17)
     $                          dip =  TF/2.0d0*(msq_b(2,17)*sub(gq,17)- 
     $                          msq_v(2,17)*subv(17))*polcolqq  
                              
                           res(k,nd) = -pdf(sign1*iflav(1),1,1)*
     $                          pdf(sign2*iflav(2),2,1)*dip
                              
                           res(k,0) = res(k,0)+res(k,nd)             
                        else
                              
c     res(k,0) = res(k,1)
                           res(k,nd) = 0.0d0
                        endif !lok(nd)

                     enddo   !loop over dipoles (jd)
                        
                     if(ldebug) then
c     take the limit of pi.pj -> 0               
                        delta = dot0p(p(0,5,1),p(0,6,1))
                        sum = 0.0d0
                        if(abs(delta).lt.1d0 .and. lok(1)) then
                           print*,"CC"                
                           do jd = 1,4
                              nd = ndmap(jd)
                              if(lok(nd)) then 
                                 sum = sum + res(k,nd)
                                 if(res(k,1).ne.0d0) then
                                    print*,"D(",ip(nd),kp(nd),jp(nd),
     -                                   ")/|m2s|",
     $                                      nd,res(k,nd)/res(k,1),lok(nd)
                                 else
                                    print*,"D(",ip(nd),kp(nd),jp(nd),")",nd
     $                                     ,res(k,nd),lok(nd)
                                 endif

                              endif
                           enddo
                           print*,"f1,f2,f3=",if1,if2,if3
                           print*,"m2real=",res(k,1)
                           print*,"dip_sum=",sum
                           if(res(k,1).ne.0d0) 
     $                          print*,"dip_sum/res(k,1)",abs(sum/res(k,1))
                           print*,"-----------------"

                        endif  
                     endif   !debug

                  endif      !endif(cc)
               enddo         !if3
            enddo            !if2
         enddo               !if1

      elseif(fsign(5).eq.fsign(6)) then !q5 or q6 in i.s.

c     compute born matrix elements with transformed momenta

#ifdef WITH_NLO
         if (nlo.eq.1) then
c     create dipole subtraction terms
            call create_ikj(ip,kp,jp)
            do jd = 1,8
               nd = ndmap(jd)
               call dips(nd,p,xuz,ip(nd),kp(nd),jp(nd),sub(1,nd),
     $              subv(nd),qvec(0,nd))
            enddo                  

c     Compute born matrix elements with transformed momenta:     

            if((fsign(1).eq.-fsign(2)).and.
     $           (fsign(3).eq.fsign(4))) then !q->qg on 21 line
     
               jdmin = 5
               jdmax = 6
            elseif((fsign(1).eq.fsign(2)).and.
     $              (fsign(3).eq.-fsign(4))) then !q->qg on 43 line
               jdmin = 7
               jdmax = 8
            else
               stop
            endif
               
            do jd = jdmin,jdmax
               nd = ndmap(jd) !mapping to dipoles set as data statement

               if (lok(nd).and.(nlo.eq.1)) then
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
            
                  if(jd.le.6) then !gQ->Qq"q"~H
                     do mu =0,3
                        qbar(mu) = p(mu,1,nd) !i.s. gluon
                           
                        kbar(mu,phystodiag(2)) = p(mu,2,nd) !i.s. quark
                        kbar(mu,phystodiag(4)) = p(mu,3,nd) !f.s. quark
                        kbar(mu,phystodiag(5)) = p(mu,4,nd)
                        kbar(mu,phystodiag(6)) = p(mu,5,nd) 
                     enddo
                     qbar(4) = 0.0d0
                     gsign = -1
                     fc(1) = fsign(1)
                     fc(2) = fsign(2)
                     fc(3) = fsign(3)
                     fc(4) = fsign(4)
     
                  elseif(jd.le.8) then !qg->qq"q"~H
                     do mu =0,3
                        qbar(mu) = p(mu,2,nd) !i.s. gluon
                        kbar(mu,phystodiag(3))= p(mu,3,nd)
                        kbar(mu,phystodiag(5))= p(mu,4,nd)
                        kbar(mu,phystodiag(1))= p(mu,1,nd)
                        kbar(mu,phystodiag(6))= p(mu,5,nd)
                     enddo
                     qbar(4) = 0.0d0
                     gsign = -1
                     fc(1) = fsign(1)
                     fc(2) = fsign(2)
                     fc(3) = fsign(3)
                     fc(4) = fsign(4)
                  else
                     stop
                  endif

                  call qqhqqj_c_virt(0,kbar,fc,qbar,gsign,nd,
     $                 uucc_c,uuss_c,ddcc_c,
     $                 ddss_c,udsc_c,ducs_c)
     
c     j = 1 is gluon emission of 21 line of born
c     j = 2 is gluon emission of 43 line of born

                  do j = 2,3
                     CCmatrix_b(0,j-1,nd)=udsc_c(j) !|M_b|^2
                     CCmatrix_b(1,j-1,nd)=ducs_c(j)
                  enddo
                              
c     for spin correlations                  
                  call qqHqqj_spcor(qvec(0,nd),kbar,fc,
     $                 qbar,gsign,nd,uucc_c,uuss_c,
     $                 ddcc_c,ddss_c,udsc_c,ducs_c)
                  do j = 2,3
                     CCmatrix_v(0,j-1,nd)=udsc_c(j) !|M_b|^2
                     CCmatrix_v(1,j-1,nd)=ducs_c(j)            
                  enddo
                     
               endif !lok(nd)
            enddo            !dipole loop
         endif !nlo.eq.1
#endif /* WITH_NLO */

C     Next the CC contributions; no b quarks allowed here
         do if3 = 1,nflVBF      !sum over all flavors 
            !3rd flavour from g->qqbar
            iflav(diagToPhys(5)) = if3*fsign(5)
            iflav(diagToPhys(6)) = if3*fsign(6)
            do if1 = 1,(nflVBF/2)*2 !4     ! sum over all flavors for quark 1 and 2
               iflav(diagToPhys(1)) = fsign(1)*if1
               iflav(diagToPhys(2)) = fsign(2)*(if1-1+2*mod(if1,2))
               do if2 = mod(if1,2)+1,(nflVBF/2)*2,2  !4,2 ! flavor of q3 set by q2 (mod 2)
                  iflav(diagToPhys(3)) = fsign(3)*if2
                  iflav(diagToPhys(4)) = fsign(4)*(if2-1+2*mod(if2,2))

                  if(nlo.eq.1) then
#ifdef WITH_NLO
                     k=FL_H4j(iflav,diagTophys,2)
                     do nd = 0,n_kin
                        res(k,nd) = 0.0d0
                     enddo
         
                     if (lok(1)) then
                        res(k,1)=pdf(sign1*iflav(1),1,1)
     &                       *pdf(sign2*iflav(2),2,1)
     &                       *CCmatrix_r(mod(if1,2),0)*polcolqq
                     else
                        res(k,1) = 0.0d0
                     endif
                     res(k,0) = res(k,1)
#endif /* WITH_NLO */

                  else       !color flow 
                     do icol = 1,2
                        k=FL_H4j(iflav,diagToPhys,icol)

                        do nd = 0,n_kin
                           res(k,nd) = 0.0d0
                        enddo
                           
                        if (lok(1)) then
                           res(k,1)=pdf(sign1*iflav(1),1,1)
     &                             *pdf(sign2*iflav(2),2,1)
     &                             *CCmatrix_r(mod(if1,2),icol)*polcolqq
                        else
                           res(k,1) = 0.0d0
                        endif
                        res(k,0) = res(k,1)
                     enddo
                  endif !nlo

c     compute subtraction terms
c     loop over dipoles
#ifdef WITH_NLO
                  if(nlo.eq.1) then
                     do jd = jdmin,jdmax
                        nd =  ndmap(jd)                        
                        if (lok(nd)) then
                           msq_b(1,nd) = CCmatrix_b(mod(if1,2),1,nd)
                           msq_b(2,nd) = CCmatrix_b(mod(if1,2),2,nd)
                           msq_v(1,nd) = CCmatrix_v(mod(if1,2),1,nd)
                           msq_v(2,nd) = CCmatrix_v(mod(if1,2),2,nd)

c     compute dipole subtraction terms                     
                           qq = 1 !labels for splittings 
                           qg = 2
                           gq = 3
                           gg = 4

                           dip = 0.0d0
                           if(nd.eq.23) 
     $                          dip = CF/2.d0*(msq_b(1,23)*sub(gq,23) + 
     $                          msq_v(1,23)*subv(23))*polcolgq 
                           if(nd.eq.24) 
     $                          dip = CF/2.0d0*(msq_b(1,24)*sub(gq,24) + 
     $                          msq_v(1,24)*subv(24))*polcolgq
                           if(nd.eq.25) 
     $                          dip = CF/2.0d0*(msq_b(2,25)*sub(gq,25) + 
     $                          msq_v(2,25)*subv(25))*polcolgq 
                           if(nd.eq.22) 
     $                          dip = CF/2.0d0*(msq_b(2,22)*sub(gq,22) + 
     $                          msq_v(2,22)*subv(22))*polcolgq 
                              
                           res(k,nd) = -pdf(sign1*iflav(1),1,1)*
     $                             pdf(sign2*iflav(2),2,1)*dip
                              
                           res(k,0) = res(k,0)+res(k,nd)             
                        else
                           res(k,nd) = 0.0d0
                        endif !lok(nd)

                     enddo   !dipoles
                  endif  !nlo.eq.1
#endif

                  if(ldebug) then
c     take the limit of pi.pj -> 0               
                     delta = dot0p(p(0,2,1),p(0,4,1))
                     sum = 0.0d0
                     if(abs(delta).lt.1d0 .and. lok(1)) then
                        print*,"CC" 
                        do jd = jdmin,jdmax
                           nd = ndmap(jd)
                           if(lok(nd)) then 
                              sum = sum + res(k,nd)
                              if(res(k,1).ne.0d0) then
                                 print*,"D(",ip(nd),kp(nd),jp(nd),
     $                                ")/|m2s|",nd,res(k,nd)/res(k,1),
     $                                lok(nd)
                              else
                                 print*,"D(",ip(nd),kp(nd),jp(nd),")",nd,
     $                                      res(k,nd),lok(nd)
                              endif
                           endif
                        enddo
                        print*,"x1,x2",xi(1),xi(2)
                        print*,"f1,f2,f3=",if1,if2,if3
                        print*,"m2real=",res(k,1)
                        print*,"dip_sum=",sum
                        if(res(k,1).ne.0d0) 
     $                       print*,"dip_sum/res(k,1)",abs(sum/res(k,1))
                        print*,"-----------------"
                     endif  
                  endif         !debug

               enddo            !if3
            enddo               !if1
         enddo                  !if2
         
      else
         stop  
      endif !fsign
     
      nmax = FL_H4j(iflav,diagToPhys,-2)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
