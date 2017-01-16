c*************************************************************************
      LOGICAL FUNCTION IsJet(ID)
c*************************************************************************
c     This function returns true, if the particle ID is consistent with 
c     considering the particle as a jet, i.e. a quark or gluon
c*************************************************************************
      implicit none
      integer id

      IsJet=(ID.eq.93).or.(ID.eq.94).or.(abs(ID).le.5).or.(ID.eq.21)

      RETURN
      END


c*************************************************************************
      SUBROUTINE defineJets(p, v, Lmax, jets_tag, Njets)
c*************************************************************************
c     This subroutine is the interface to various jet-defining functions.
c     All jets are considered: those from the parton-array and those, which
c     are in the decay-particle array.
c     NOTE: so far the R-separation parameter is defined just here. This
c           has to be modified such, that it can be set in some input file.
c*************************************************************************
*     INPUT 
*     p     : the array containing the partons from the basic 2->2(3)+X process
*     Lmax  : number of momentum configurations
*     v     : the other particles
*
*     OUTPUT
*     jets  : 4-momenta, mass, pt, rapidity and phi-angle of the jets after 
*                jet-definition
*     Njets : the number of combined jets
c*************************************************************************
      implicit none
#include "global.inc"
#include "cuts.inc"

      real*8 partons(0:7,max_jets,max_kin)
      real*8 jets(0:7,max_jets,max_kin),jets_tag(0:7,max_jets,max_kin)
      real*8 p(0:3,max_p,max_kin)
      real*8 v(0:3,max_v,max_kin)
      integer L, Lmax, Njets(max_kin), Npartons(max_kin), NpartonsInCuts
      integer partInSameJet
      
      real*8 pt,y,phi
      integer jrec
      integer i,j,mu,icount

      logical IsJet

      common/partons/partons, Npartons

      logical ldebug
      parameter (ldebug=.false.)


**c common block for transferring information on parton content of jets
**c parton inside jet identified by idef_part(jetID,partonID)
      integer idef_part(max_p+max_v-nx, max_p+max_v-nx)
      common /defpartkt/ idef_part


c-----------------------------------------------------------

      if (ldebug) then
         print*, 'in jetdefinition.F'
         print*, n_p, nx, Lmax
         print*, n_v 
         print*, particle_IDs
         print*, p(0:3,3,1)
      endif


      jrec = -1                 ! recombination scheme

      do L = 1, Lmax ! for all momentum configurations
c     select all jets from decay-particles and partons and fill everything into "partons"
         Npartons(L)=0
         if (L.eq.1) then  
            do i = nx+1, n_p
               do mu=0,3
                  partons(mu,Npartons(L)+1,L) = p(mu,i,L)
               end do
               Npartons(L)=Npartons(L)+1
            end do
         else     ! L > 1 are the ptilde configurations
            do i = nx+1, n_p-1
               do mu=0,3
                  partons(mu,Npartons(L)+1,L) = p(mu,i,L)
               end do
               Npartons(L)=Npartons(L)+1
            end do
         endif

         ! it is important that the partons from p are filled in first for the semileptonic
         ! decays with SINGLE_DECAYJET = 1  !
         do i = 1, n_v
            if ( IsJet(particle_IDs(i)) ) then
               do mu=0,3
                  partons(mu,Npartons(L)+1,L) = v(mu,i,L)
               end do
               Npartons(L)=Npartons(L)+1
            endif
         enddo
c     one has to calculate rapidity, pt, phi-angle for the jkt-algorithms to work

         do i=1, Npartons(L)
            call lego( partons(0,i,L), pt, y, phi )
            partons(4,i,L) = 0.0   ! mass  (you see that all jets are treated as massless)
            partons(5,i,L) = pt    ! transverse momentum  
            partons(6,i,L) = y     ! rapidity 
            partons(7,i,L) = phi   ! azimuthal angle
            if (ldebug) print*, 'pt y phi', pt, y, phi
         end do

         njets(L)=0 
         do i=1,Npartons(L)
           if(abs(partons(6,i,L)).lt.ymax_part) then
              njets(L)=njets(L)+1
              do mu=0,7
                 jets(mu,njets(L),L)=partons(mu,i,L)
              enddo
           endif
         enddo
         NpartonsInCuts = njets(L)

c     now call the jet-defining functions
c     the plan is to switch here between different jet-algorithms
      
         call jkt(njets(L), Rjjmin, jrec, min_jets, jetgenktp, jets(0,1,L)) 

c     for semileptonic decays:
c     check if both decay jets form a single jet
c     partons from the hadronic decay are Npartons(L) and Npartons(L)-1
c     (and NpartonsInCuts / NpartonsInCuts-1 after the eta cut),
c     as they com from v() which is filled last into partons().
         if (semilep_fatjet.eq.1) then
           decayPartonsCombined(L)=.false.
           ! both decay partons must lie inside parton cuts
           if (abs(partons(6,Npartons(L),L)).lt.ymax_part .and. abs(partons(6,Npartons(L)-1,L)).lt.ymax_part) then
             do i = 1, njets(L)
               partInSameJet = 0
               do j = 1, NpartonsInCuts
                 if ( idef_part(i,j).eq.NpartonsInCuts .or. idef_part(i,j).eq.(NpartonsInCuts-1) ) 
     &                  partInSameJet =  partInSameJet + 1
               enddo
               if (partInSameJet.eq.2) decayPartonsCombined(L)=.true.
               if (partInSameJet.gt.2) print*, "Error..."
             enddo
           endif
         endif

c     sort out the jets, which are not inside the detector

         icount=0
         do i=1,njets(L)
            if (abs(jets(6,i,L)).lt.ymax_def.and.jets(5,i,L).gt.ptmin_def(icount+1)) then
               icount=icount+1
               do mu=0,7
                  jets_tag(mu,icount,L)=jets(mu,i,L)
               enddo 
            endif
         enddo
         njets(L)=icount
      enddo
      ! set flag that includes bbar final states only if they are
      ! recombined into one jet
      lincludebbfinaldo = lincludebbfinal .and. (njets(1) .eq. 1)

      if(ldebug) then 
         print*, 'njets', njets
         print*, jets(:,1,1)
         print*, jets(:,2,1)
         print*, jets(:,3,1)
      endif

      RETURN
      END


C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C       k_T-scheme as defined in
C       M.Seymour, Nucl.Phys.B513:269-300,1998 [hep-ph/9707338]
C
c	Initial version:  1998, April 30
C	Last modified  :  2003, April 9
C
C       modified for VBFNLO:  2006, September 16
c          - removed the "defpartkt"-stuff 
c          - adapted array boundaries to the new "global.inc"
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine jkt( njet,       !in and out: number of identified clusters
     &                rjjmin,     !in: size of jet
     &                jrec,       !in: ID for recombination scheme: 1=E
     &                njetmin,    !in: minimal numbers of surviving clusters
     &                jetp,       !in: exponent for generalised kT
     &                p_j)        !in and out: recombined partons

      implicit none

c declare input/output variables and cut parameters
c
#include "global.inc"
c
      integer njet,jrec,njetmin
      real*8 rjjmin, jetp, p_j(0:7,max_jets), p_ja(0:7,max_jets)
c
      real*8 rjj,dij,dib, rsepf
      external rsepf
      integer idef_jetl(max_jets), idpl(max_jets,max_jets)
      integer i,j,ii,jj,iaux,mu,ik,jk
      integer ipoint1,ipoint2, ipoint3,j_complete
c
**c common block for transferring information on parton content of jets
**c parton inside jet identified by idef_part(jetID,partonID)
**      real*8 p_part(0:7,np-nx,3)
**      common /defpartkt/ p_part, idef_part
      integer idef_part(max_p+max_v-nx, max_p+max_v-nx)
      common /defpartkt/ idef_part
c
c define flags
c
      logical ldebug
      parameter (ldebug = .false. )
      integer ipass, icount
c      data icount/0/
c      icount = icount+1
      
      j_complete = 0
c
c prepare information for jet shapes plots: save initial parton momenta
      do ik = 1,max_jets
         idef_jetl(ik) = 0
         do jk = 1,max_jets
            idpl(jk,ik)=0
         enddo
         do mu = 0,7
**            p_part(mu,ik,L) = p_j(mu,ik)
            p_ja(mu,ik) = p_j(mu,ik)
         enddo
      enddo
      do i = 1,njet
         idef_jetl(i) = i
         idpl(i,1) = i
      enddo
     
c      do ii = 1,njet
c         write(6,11) " P_j(",ii,") = ",(p_j(mu,ii),mu=0,7)
c      enddo
      if (ldebug) then
         ipass = 0
         write(6,*) " "
         write(6,*) "------------------------------------------------- "
         write(6,*) " new event in jkt, #",icount
         write(6,*) " input variables "
         do ii = 1,njet
            write(6,11) " P_j(",ii,") = ",(p_j(mu,ii),mu=0,7)
         enddo
         write(6,12) " idef_jetl = ",idef_jetl
 11      format(a,i2,a,8f12.6)
 12      format(a,10i3)
      endif
c
 50   continue       
      ipoint1 = 0
      ipoint2 = 0
      ipoint3 = 0
c
      if (ldebug) then
         ipass = ipass+1
         write(6,*) " cluster momenta at pass # ",ipass
         do ii = 1,njet
            i = idef_jetl(ii)
            write(6,"(2i3,3f12.4)")ii,i,(p_ja(j,i),j=5,7)
         enddo
      endif
c
c  step1: find cluster closest to the beam
      dib = 1d33
      do ii = j_complete+1, njet
         i = idef_jetl(ii)
         rjj = (p_ja(5,i)**jetp)*rjjmin
         if (rjj .lt. dib) then
            dib     = rjj
            ipoint3 = ii
         endif
      enddo
c
      if ((ldebug) .and. (ipoint3.gt.0)) then 
        i = idef_jetl(ipoint3)
        print*," idef_jetl(ipoint3) = ",i
      endif
c
c  step2: search for pair of closest partons
      dij = 1d33
      do ii = j_complete+1,njet-1
         i = idef_jetl(ii)
         do jj = ii+1, njet
            j = idef_jetl(jj)
            rjj = rsepf( p_ja(0,i), p_ja(0,j) )
            rjj = min( (p_ja(5,i)**jetp), (p_ja(5,j)**jetp) )*rjj
            if(ldebug) then
               print*," ii, jj = ",ii,jj," j_complete = ",j_complete
               print*," dij = ",rjj
            endif
            if (rjj .lt. dij) then
               dij     = rjj
               ipoint1 = ii 
               ipoint2 = jj
            endif
         enddo
      enddo
c
c  step3: merge i and j if dij < dib, raise number of complete jets otherwise
      if (ldebug) then
         write (6,*) " ipass = ",ipass," dij = ",dij," dib = ",dib
         write (6,"(3(a,i3))") " ipoint1 = ",ipoint1," ipoint2 = ",
     1               ipoint2, " ipoint3 = ",ipoint3
      endif
c
      if ( dij .lt. dib ) then
         if (ldebug) then 
            write (6,*)" join clusters # ",ipoint1,ipoint2
         endif
         njet = njet-1
         i = idef_jetl(ipoint1)
         j = idef_jetl(ipoint2)
         call rec( jrec,p_ja(0,i),p_ja(0,j) )
         call merge_idef(idpl(1,1),i,j)
c
         do ii = ipoint2,njet
            idef_jetl(ii) = idef_jetl(ii+1)
         enddo
      else
         j_complete = j_complete+1
         if (ldebug) then
            write(6,*)" increase j_complete to ",j_complete,
     1                " njet = ",njet
         endif
C move pointer to ipoint3 to number j_complete, 
         if (ipoint3.gt.j_complete) then
            iaux = idef_jetl(ipoint3)
            if (ldebug) then
               write(6,*) " move jet ",ipoint3,
     1                    " to position ",j_complete
            endif
            do ii = ipoint3, j_complete+1, -1
               idef_jetl(ii) = idef_jetl(ii-1)
            enddo
            idef_jetl(j_complete) = iaux
            ipoint3 = j_complete
         endif
      endif
c
      if (njet .gt. j_complete .and. njet.ge.njetmin) goto 50
c
c   reorder clusters according to E_T
      do i = 1, njet-1
         do j = i+1,njet
            if ( 
     1       p_ja(5,idef_jetl(i)).lt.p_ja(5,idef_jetl(j)) 
     2      ) then
               iaux = idef_jetl(j)
               idef_jetl(j) = idef_jetl(i)
               idef_jetl(i) = iaux
            endif
         enddo
      enddo
c
c map momenta back onto pj
      do ii=1,njet
         i = idef_jetl(ii)
         do mu=0,7
            p_j(mu,ii) = p_ja(mu,i)
         enddo
         do j = 1,n_p+n_v-nx
            idef_part(ii,j) = idpl(i,j)
         enddo
      enddo
      do ii=njet+1,max_jets
         do mu=0,7
            p_j(mu,ii) = 0
         enddo
         do j = 1,n_p+n_v-nx
            idef_part(ii,j) = 0
         enddo
      enddo
c
c         write(6,*) " final cluster pT,y,phi in kT "
c         do i = 1,njet
c            write(6,"(i3,8f12.4,a,10i3)")i,(p_j(j,i),j=0,7),
c     1           " with partons",(idef_part(i,jj,L),jj=1,njet)
c         enddo

**      if (ldebug) then
**         write(6,*) " final cluster pT,y,phi in kT "
**         do i = 1,np-nx
**            write(6,"(i3,3f12.4,a,10i3)")i,(p_j(j,i),j=5,7),
**     1           " with partons",(idef_part(i,jj,L),jj=1,np-nx)
**         enddo
**        if (L.eq.1.and.njet.le.2) read 113,text
** 113     format(a)
**      endif
c
      return
      end


c
c------------  merge_idef  ------------------------------------------
c
      subroutine merge_idef(idef_part,iget,igive)
      implicit none
c move all nonzero entries from idef_part(igive,*) to idef_part(iget,*)
#include "global.inc"
      integer idef_part(max_jets,max_jets),iget,igive
      integer i,j,k,nj,ni
      j = 1
      nj = idef_part(igive,j)
      do while(nj.ne.0)
         do i = 1,max_jets
            ni = idef_part(iget,i)
            if(nj.eq.ni) then
               goto 10              ! already in iget list, move to next
            elseif(nj.lt.ni) then   ! insert nj into iget list
               do k = max_jets,i+1,-1
                  idef_part(iget,k) = idef_part(iget,k-1)
               enddo
               idef_part(iget,i) = nj
               goto 10
            elseif(ni.eq.0) then
               idef_part(iget,i)=nj ! mv nj onto iget list
               goto 10
            endif
         enddo
 10      j = j+1
         nj = idef_part(igive,j)
      enddo
      do i = 1,max_jets
         idef_part(igive,i)=0
      enddo
      end

C
C********************  REC   ***************************************
C
      subroutine rec(jrec,p1,p2)
c     Initial version:  1995 Oct. 1   by Dieter Zeppenfeld
*     Modified       :  1998 Jan.  6: by Jose Repond; 
*                                     implemented Et recombination: jrec=4
c     Last modified  :  1998 Dec. 29    by Dieter Zeppenfeld

*     p1, p2 arrays: 0 = E
*                    1 = px
*                    2 = py
*                    3 = pz
*                    4 = mass
*                    5 = pt if jrec < 0, else ET
*                    6 = eta
*                    7 = phi

c     jrec :   in, recombination flag
c     |jrec| = 1 : E-scheme, p_rec = p1+p2
c              2 : E0 scheme, p_rec(0) = p1(0)+p2(0)
c              3 : p scheme,  p_rec(i) = p1(i)+p2(i)  i = 1,2,3
c      jrec =  4 : Et-weighted recombination to massless parton
c     p1,p2:   in, two 8-momenta to be recombined according to jrec
c     p1   :   out, the recombined momentum is put onto p1
      implicit none
      real*8 p1(0:7), p2(0:7)
      real*8 delta
      real*8 pt12, eta12, phi12, pt1, pt2
      real*8 r_pi, twopi
      parameter( r_pi=3.141592653589793238d0, twopi=2d0*r_pi )
      integer jrec, mu

      if(jrec.le.4) then
        do mu = 0,3
           p1(mu) = p1(mu) + p2(mu)
        enddo
        pt12 = p1(5)+p2(5)
      endif

      if ( abs(jrec).eq.1 ) then
         p1(4) = sqrt( p1(0)**2 - p1(1)**2 - p1(2)**2 - p1(3)**2 )
      elseif ( abs(jrec).eq.2 ) then
         delta = p1(0) / sqrt( p1(1)**2+p1(2)**2+p1(3)**2 )
         do mu = 1,3
            p1(mu) = delta*p1(mu)
         enddo
         p1(4) = 0d0
      elseif  (abs(jrec) .eq. 3) then
         p1(0) = sqrt( p1(1)**2+p1(2)**2+p1(3)**2 )
         p1(4) = 0d0
      elseif  (jrec .eq. 4) then
c recombination to massive cluster
         pt1 = p1(5)
         pt2 = p2(5) 
         pt12=pt1+pt2
         eta12=(pt1*p1(6)+pt2*p2(6))/pt12
         if (abs(p1(7)-p2(7)).le.r_pi) then
            phi12=(pt1*p1(7)+pt2*p2(7))/pt12
         elseif (p1(7).lt.p2(7)) then
            phi12=(pt1*(p1(7)+twopi)+pt2*p2(7))/pt12
         else
            phi12=(pt1*p1(7)+pt2*(p2(7)+twopi))/pt12
         endif
         if (phi12.gt.r_pi) then
            phi12 = phi12 - twopi
         elseif (phi12.lt.-r_pi) then
            phi12 = phi12 + twopi
         endif
         p1(4)=p1(0)**2-p1(1)**2-p1(2)**2-p1(3)**2
         if (p1(4).gt.0) then
            p1(4) = sqrt(p1(4))
         elseif (p1(4).lt.0) then
            p1(4) = -sqrt(-p1(4))
         else
            p1(4) = 0d0
         endif
         p1(5)=pt12
         p1(6)=eta12
         p1(7)=phi12
         return
      elseif  (jrec .eq. 5) then
c recombination to massless cluster
         pt1 = sqrt(p1(1)**2+p1(2)**2)
         pt2 = sqrt(p2(1)**2+p2(2)**2)
         pt12=pt1+pt2
         eta12=(pt1*p1(6)+pt2*p2(6))/pt12
         if (abs(p1(7)-p2(7)).le.r_pi) then
            phi12=(pt1*p1(7)+pt2*p2(7))/pt12
         elseif (p1(7).lt.p2(7)) then
            phi12=(pt1*(p1(7)+twopi)+pt2*p2(7))/pt12
         else
            phi12=(pt1*p1(7)+pt2*(p2(7)+twopi))/pt12
         endif
         if (phi12.gt.r_pi) then
            phi12 = phi12 - twopi
         elseif (phi12.lt.-r_pi) then
            phi12 = phi12 + twopi
         endif
         p1(0)=pt12*cosh(eta12)
         p1(1)=pt12*cos(phi12)
         p1(2)=pt12*sin(phi12)
         p1(3)=pt12*sinh(eta12)
         p1(4)=0.0d0
         p1(5)=pt12
         p1(6)=eta12
         p1(7)=phi12
         return
      elseif  (jrec .eq. 6) then
c recombination to massless cluster, but p1(5) = ET of massive jet
         pt1 = sqrt(p1(1)**2+p1(2)**2)
         pt2 = sqrt(p2(1)**2+p2(2)**2)
         pt12=pt1+pt2
         eta12=(pt1*p1(6)+pt2*p2(6))/pt12
         if (abs(p1(7)-p2(7)).le.r_pi) then
            phi12=(pt1*p1(7)+pt2*p2(7))/pt12
         elseif (p1(7).lt.p2(7)) then
            phi12=(pt1*(p1(7)+twopi)+pt2*p2(7))/pt12
         else
            phi12=(pt1*p1(7)+pt2*(p2(7)+twopi))/pt12
         endif
         do mu = 0,3
            p1(mu) = p1(mu)+p2(mu)
         enddo
         p1(4) = sqrt( p1(0)**2-p1(1)**2-p1(2)**2-p1(3)**2 )
         pt12 = sqrt( p1(0)**2 - p1(3)**2 )
         if (phi12.gt.r_pi) then
            phi12 = phi12 - twopi
         elseif (phi12.lt.-r_pi) then
            phi12 = phi12 + twopi
         endif
         p1(0)=pt12*cosh(eta12)
         p1(1)=pt12*cos(phi12)
         p1(2)=pt12*sin(phi12)
         p1(3)=pt12*sinh(eta12)
         p1(5)=pt12
         p1(6)=eta12
         p1(7)=phi12
         return
      endif
c this is for jrec = +-1, +-2, +-3 only
      call legoy( p1, p1(5), p1(6), p1(7) )
      if (jrec.gt.0) then
         p1(5) = pt12
      endif

      return
      end
