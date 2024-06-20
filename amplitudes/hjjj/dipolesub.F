c     Terrance Figy
c     last modified 18 April 2005
c     This subroutine computes the dipole factors,
c     the q_mu for a set of ptilde or transformed momenta
      
      subroutine dips(nd,p,xu,ip,kp,jp,sub,subv,vec)
      implicit none
c     input:p(0:3,np,nd) parton momenta
c           Only the parton momenta are needed!
c     nd = dipole kinematic conf. number
c     xu(2,2:nd) 
c     ip emitter, kp is emitted, and jp is spectator
c     output: dipole factors: sub(4),subv
c     vec(mu) is q_mu for config. nd
c     Will be used in constructing the |mborn|^2.
c
#include "VBFNLO/utilities/global.inc"
c     
      integer ii,nd,mu,ip,jp,kp,qq,qg,gq,gg
      real*8 p(0:3,max_p,max_kin),xu(2,2:max_kin)
      real*8 sub(4),subv
      real*8 vec(0:3) ! this q_mu for ik,j
      real*8 y,omy,z,omz,u,omu,x,omx
      real*8 sik,sjk,sij,vecsq
      real*8 dot0p
      external dot0p
      real*8 gsq,fpi                 !g_strong^2 = 4Pi alpha_s
      parameter(fpi=4d0*pi)
      real*8 one,two,four ! will be put into global
      parameter(one=1.0d0,two=2.0d0,four=4.0d0)
      logical ldebug
      parameter(ldebug = .false.)
c     get g^strong
#include "VBFNLO/utilities/scales.inc"
ccc
c     set gsq 
c     included here is an extra factor of alpha_s 
c     it was in qqhqqgg originally. 
c
      gsq = fpi*als(1,nd)!*als(2,nd) 
c      print*,'gsq=',gsq
c    
c     initialize sub and subv
      do ii=1,4
         sub(ii) = 0.0d0
      enddo
      subv=0.0d0
      do mu =0,3
         vec(mu) = 0.d0
      enddo

c     Label splittings
      qq = 1
      qg = 2
      gq = 3
      gg = 4
c      
c     dot products sik is p_i /cdot p_k
c
c      sik = p(0,ip,1)*p(0,kp,1)-p(1,ip,1)*p(1,kp,1)-
c     $     p(2,ip,1)*p(2,kp,1)-p(3,ip,1)*p(3,kp,1)
c      sjk = p(0,jp,1)*p(0,kp,1)-p(1,jp,1)*p(1,kp,1)-
c     $     p(2,jp,1)*p(2,kp,1)-p(3,jp,1)*p(3,kp,1)
c      sij = p(0,ip,1)*p(0,jp,1)-p(1,ip,1)*p(1,jp,1)-
c     $     p(2,ip,1)*p(2,jp,1)-p(3,ip,1)*p(3,jp,1) 
c
      sik = dot0p(p(0,ip,1),p(0,kp,1))
      sjk = dot0p(p(0,jp,1),p(0,kp,1))
      sij = dot0p(p(0,ip,1),p(0,jp,1))
c
      if(ldebug) print*,'sik,sjk,sij=',sik,sjk,sij
c     final-final
      if((ip.gt.2).and.(jp.gt.2)) then
         if(ldebug) print*,'final-final'
         z = xu(2,nd)
         y = xu(1,nd) 
         omy = one - y
         omz = one - z
c     q_mu            
         do mu = 0,3
            vec(mu) = z*p(mu,ip,1) - omz*p(mu,kp,1)
         enddo 
         if(ldebug) print*,'vec dot pt_i=',dot0p(vec,p(0,ip,nd))
c     dipole factors
         sub(gg) = gsq/sik*(two/(one-z*omy) + two/(one-omz*omy)-four)
         sub(gq) = gsq/sik
         sub(qq) = gsq/sik*(two/(one-z*omy)-one-z)
c     for spin corr. D^ikj_ff,gg = -D^ikj_ff,gq
         subv = gsq/sik * two/sik  !gg

c     Final-initial
c
      elseif((ip.gt.2).and.(jp.le.2)) then
         if(ldebug) print*,'final-initial'
         x = xu(1,nd)  
         z = xu(2,nd) 
         omx = one - x
         omz = one - z
c     q_mu            
         do mu = 0,3
            vec(mu) = z*p(mu,ip,1) - omz*p(mu,kp,1)
         enddo 
         if(ldebug) print*,'vec dot pt_i=',dot0p(vec,p(0,ip,nd))
c     dipole factors
         sub(gg) = gsq/sik/x*(two/(omz+omx)+two/(omx+z)-four)
         sub(gq) = gsq/sik/x
         sub(qq) = gsq/sik/x*(two/(omz+omx)-one-z)
c     for spin corr. gg = -gq 
         subv = gsq/sik/x*two/sik !gg

c     Initial-final
         elseif((ip.le.2).and.(jp.gt.2)) then
            if(ldebug) print*,'initial-final'
c
            x = xu(1,nd)
            u = xu(2,nd)
            omx = one - x
            omu = one - u
c     q_mu            
           do mu = 0,3
              vec(mu) = p(mu,kp,1)/u - p(mu,jp,1)/omu
           enddo 
           if(ldebug) print*,'vec dot p_i=',dot0p(vec,p(0,ip,1))
c     dipole factors
            sub(qq) = gsq/x/sik*(two/(omx+u)-(one+x))
            sub(qg) = gsq/x/sik*(one-two*x*omx)
            sub(gq) = gsq/sik
            sub(gg) = gsq/x/sik*(two/(omx+u)-two+two*x*omx)
c     for spin corr. gg = gq
            subv = gsq/x/sik*two*omx*u*omu/x/sjk
c
c     initial -initial
         elseif((ip.le.2).and.(jp.le.2)) then 
            if(ldebug) print*,'initial-initial'
        
            x = xu(1,nd)
            omx = one - x 
            sub(qq) = gsq/x/sik*(two/omx -(one+x))
            sub(qg) = gsq/x/sik*(one - two*x*omx)
            sub(gq) = gsq/sik
            sub(gg) = gsq/x/sik*(two*x/omx + two*x*omx)
c     for spin corr. gg = gq
c     q_mu = vec(mu)
            do mu = 0,3
               vec(mu) = p(mu,kp,1) - sik*p(mu,jp,1)/sij
            enddo
            if(ldebug) print*,'vec dot p_i=',dot0p(vec,p(0,jp,1))
c
c     q dot q = vecsq
            vecsq = vec(0)**2-vec(1)**2 -vec(2)**2 -vec(3)**2
c            
            subv = -four*gsq/x/sik*omx/x/vecsq 
            
         else
            print*,'program will stop'
            stop
         endif

      end
c
