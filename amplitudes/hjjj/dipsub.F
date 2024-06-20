c     This subroutine computes the dipole subtractions 
c     for q q -> q q g g H processes.
c
      subroutine dipsub(gsign,fsign,msq_b,msq_v,sub,subv,dip)
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "color.inc"
c
      real*8 dip(2:max_kin),sub(4,2:max_kin),subv(2:max_kin)
      real*8 msq_b(0:2,2:max_kin),msq_v(0:2,2:max_kin)
      integer gsign(2),fsign(4)          !sign factors for gluons
      integer qq,gg,qg,gq,nd
      real*8 polcolqq,polcolgq    
c
c     dipole subtraction terms
c     index labels for spitting functions
      qq = 1 
      qg = 2
      gq = 3
      gg = 4
c
      polcolqq = 1d0/(NC*NC*4d0)
      polcolgq = 1d0/((NC**2-1)*NC*4d0)
      do nd = 2,max_kin
         dip(nd) = 0.0d0
      enddo
c
      if((gsign(1).eq.1).and.(gsign(2).eq.1)) then !two final state gluons
c     Include subtraction terms for gluons splitting into q qbar
c     This needs aspect needs adjusting.
c
         dip(2)= CAc/2.0d0 * msq_b(1,2)*sub(qq,2) 
c     
         dip(3)= CAc/2.0d0 * msq_b(1,3)*sub(qq,3)
c     
         dip(4) = CAc/2.d0 * (msq_b(1,4)*sub(gg,4) + 
     $        msq_v(1,4)*subv(4))     
c     
         dip(5)= sub(qq,5)*(cf*msq_b(2,5) + 
     $        (cf-CAc/2.0d0)*msq_b(1,5))
c
         dip(6)=(cf*msq_b(2,6)+
     $        (cf-CAc/2.0d0)*msq_b(1,6))*sub(qq,6)
c     
         dip(7)=CAc/2.0d0 * (msq_b(1,7)*sub(gg,7)+ !g->gg
     $        msq_v(1,7)*subv(7))    
c
         dip(8)= sub(qq,8)*(cf*msq_b(2,8) + 
     $        (cf-CAc/2.0d0)*msq_b(1,8))
c
         dip(9) = (cf*msq_b(2,9)+
     $        (cf-CAc/2.0d0)*msq_b(1,9))*sub(qq,9)
c
         dip(10) = CAc/2.0d0 * msq_b(1,10) * sub(qq,10)
c     
         dip(11) = CAc/2.0d0 * msq_b(1,11) * sub(qq,11)
ccccccccc
         dip(12)= CAc/2.0d0 * msq_b(2,12)*sub(qq,12)
c
         dip(13)= CAc/2.0d0 * msq_b(2,13)*sub(qq,13)
c
         dip(14) = CAc/2.d0 * (msq_b(2,14)*sub(gg,14)+ !g->gg
     $        msq_v(2,14)*subv(14))    
c
         dip(15)= sub(qq,15)*(cf*msq_b(1,15) + 
     $        (cf-CAc/2.0d0)*msq_b(2,15))
c
         dip(16)=(cf*msq_b(1,16)+
     $        (cf-CAc/2.0d0)*msq_b(2,16))*sub(qq,16)
c
         dip(17)=CAc/2.0d0 * (msq_b(2,17)*sub(gg,17)+ !g->gg
     $        msq_v(2,17)*subv(17))   
c
         dip(18)= sub(qq,18)*(cf*msq_b(1,18) + 
     $        (cf-CAc/2.0d0)*msq_b(2,18))
c
         dip(19) = (cf*msq_b(1,19)+
     $        (cf-CAc/2.d0)*msq_b(2,19))*sub(qq,19)
c
         dip(20) = CAc/2.0d0 * msq_b(2,20) * sub(qq,20)
c     
         dip(21) = CAc/2.0d0 * msq_b(2,21) * sub(qq,21)
c
      elseif((gsign(1).eq.-1).and.(gsign(2).eq.1)) then !g1 initial g2 final
         if(fsign(1).eq.-fsign(2)) then !i.s. gluon on 21 line
c     
            dip(9) = Tf *(msq_b(2,9)+
     $           (cf-CAc/2d0)*msq_b(1,9)/cf)*sub(qg,9)*polcolqq
c     
            dip(11) = Tf *
     $           CAc/2d0*msq_b(1,11)/cf*sub(qg,11)*polcolqq
c     
            dip(23) = Tf *(msq_b(2,23) +
     $           (cf- CAc/2d0)*msq_b(1,23)/cf)*sub(qg,23)*polcolqq
c     
            dip(24) = Tf *
     $           CAc/2d0*msq_b(1,24)/cf*sub(qg,24)*polcolqq
c
            nd = 8
            dip(nd)= CAc/2d0 *(msq_b(1,nd)*sub(gg,nd)+
     $           msq_v(1,nd)*subv(nd))*polcolgq
            nd = 10
            dip(nd) = CAc/2d0 *(msq_b(1,nd)*sub(gg,nd)+
     $           msq_v(1,nd)*subv(nd))*polcolgq
            nd = 18
            dip(nd) = CF*msq_b(1,nd)*sub(qq,nd)*polcolgq
            nd = 15
            dip(nd) = CF*msq_b(1,nd)*sub(qq,nd)*polcolgq
            nd = 5
            dip(nd) = CAc/2d0*msq_b(1,nd)*sub(qq,nd)*polcolgq
            nd = 7
            dip(nd) = CAc/2d0*msq_b(1,nd)*sub(qq,nd)*polcolgq
            nd = 2
            dip(nd) = (CF-CAc/2d0)*msq_b(1,nd)*
     $           sub(qq,nd)*polcolgq
            nd = 4
            dip(nd) = (CF-CAc/2d0)*msq_b(1,nd)*
     $           sub(qq,nd)*polcolgq
         elseif(fsign(3).eq.-fsign(4)) then !i.s. gluon on 43 line
ccc   work these out! these do not work!
            nd = 19
            dip(nd) = Tf *(msq_b(1,nd)+
     $           (cf-CAc/2d0)*msq_b(2,nd)/cf)*sub(qg,nd)*polcolqq
c     
            nd = 22
            dip(nd) = Tf *
     $           CAc/2d0*msq_b(2,nd)/cf*sub(qg,nd)*polcolqq
c     
            nd = 25
            dip(nd) = Tf *(msq_b(1,nd) +
     $           (cf-CAc/2d0)*msq_b(2,nd)/cf)*sub(qg,nd)*polcolqq
c     
            nd = 21
            dip(nd) = Tf *
     $           CAc/2d0*msq_b(2,nd)/cf*sub(qg,nd)*polcolqq
cccc  
            nd = 18
            dip(nd)= CAc/2d0 *(msq_b(2,nd)*sub(gg,nd)+
     $           msq_v(2,nd)*subv(nd))*polcolgq
            nd = 20
            dip(nd) = CAc/2d0 *(msq_b(2,nd)*sub(gg,nd)+
     $        msq_v(2,nd)*subv(nd))*polcolgq
            nd = 8
            dip(nd) = CF*msq_b(2,nd)*sub(qq,nd)*polcolgq
            nd = 5
            dip(nd) = CF*msq_b(2,nd)*sub(qq,nd)*polcolgq
            nd = 15
            dip(nd) = CAc/2d0*msq_b(2,nd)*sub(qq,nd)*polcolgq
            nd = 17
            dip(nd) = CAc/2d0*msq_b(2,nd)*sub(qq,nd)*polcolgq
            nd = 12
            dip(nd) = (CF-CAc/2d0)*msq_b(2,nd)*
     $           sub(qq,nd)*polcolgq
            nd = 14
            dip(nd) = (CF-CAc/2d0)*msq_b(2,nd)*
     $           sub(qq,nd)*polcolgq
c         stop
            
         endif
         
      elseif((gsign(1).eq.-1).and.(gsign(2).eq.-1)) then !two i.s. gluons
         nd = 22                !D^b2,4
         dip(nd) = TF * msq_b(1,nd)*sub(qg,nd)*polcolgq
         nd = 18                !D^b4,2
         dip(nd) = TF * msq_b(1,nd)*sub(qg,nd)*polcolgq
         nd = 23                !D^a1,3
         dip(nd) = TF * msq_b(2,nd)*sub(qg,nd)*polcolgq
         nd = 9                 !D^a3,1
         dip(nd) = TF * msq_b(2,nd)*sub(qg,nd)*polcolgq
      else
         print*,'program will stop'
         stop
      endif
    
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     subroutine which transforms the np parton flavors of the real emission
c     diagram to the (np-1) parton flavors for the born diagram
c
c     transform PhysToDiag map
c
      subroutine transflav(ip,kp,flavp,phystodiag,phystodiagbrn)
      implicit none
#include "VBFNLO/utilities/global.inc"
c     input/output
      integer flavp(max_p),ip,kp
      integer phystodiag(max_p),phystodiagbrn(max_p)
c     flavp(np) = 1 for a quark
c               =-1 for antiquark
c               = 0 for gluon
c               =-2 not a parton 
c     ip is the emitter and kp is the emitted parton
c
c     local variables
      integer lflavp(max_p),ipart,ic,lphystodiag(max_p)

      do ic = 1,n_p
         lflavp(ic)=flavp(ic)
         lphystodiag(ic) = phystodiag(ic)
      enddo
c     set last to -2
c      lflavp(np) = -2           !this is not a parton,dummy variable
c
      if(ip.le.2) then          !ip in initial state
         ipart = 1
         do ic=1,n_p
            if(ic.eq.ip) then   !i emits k 
               flavp(ipart)=lflavp(ip)-lflavp(kp)
c
               if(lflavp(ip).eq.0 .and. lflavp(kp).eq.0) then 
c     gluon -> gluon + gluon
                  physToDiagBrn(ipart) = lphysToDiag(ip)
               elseif(lflavp(ip).eq.0 .and. lflavp(kp).ne.0) then
c     gluon -> quark + anti-quark
                  physToDiagBrn(ipart) = lphysToDiag(kp)
               elseif(lflavp(ip).ne.0 .and. lflavp(kp).eq.0) then
c     quark -> quark + gluon
                  physToDiagBrn(ipart) = lphysToDiag(ip) 
                  
               endif
              
            elseif(ic.eq.kp) then
               goto 20
            else
               flavp(ipart) = lflavp(ic)
               physToDiagBrn(ipart) = lphysToDiag(ic) 
            endif         
            ipart = ipart + 1
 20      continue
      enddo

      elseif(ip.gt.2) then      !ip in final state
         ipart = 1
         do ic=1,n_p
            if(ic.eq.ip) then   !i emits k 
               flavp(ipart)=lflavp(ip)+lflavp(kp)
c     
                if(lflavp(ip).eq.0 .and. lflavp(kp).eq.0) then 
c     gluon -> gluon + gluon
                   physToDiagBrn(ipart) = lphysToDiag(ip)
                elseif(lflavp(ip).ne.0 .and. lflavp(kp).eq.0) then
c     quark -> quark + gluon
                   physToDiagBrn(ipart) = lphysToDiag(ip)
                elseif(lflavp(ip).ne.0 .and. lflavp(kp).ne.0) then
c     gluon -> quark + anti-quark
                   physToDiagBrn(ipart) = lphysToDiag(ip)
               endif
            elseif(ic.eq.kp) then
               goto 30
            else
               flavp(ipart) = lflavp(ic)
               physToDiagBrn(ipart) = lphysToDiag(ic) 
            endif         
            ipart = ipart + 1
 30      continue
      enddo  
      endif
c     set last to -2
      flavp(n_p) = -2    
      phystodiagbrn(n_p) = -2
      return
      end
