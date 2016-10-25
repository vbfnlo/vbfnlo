      subroutine oxxxxx(p,fmass,nhel,nsf , fo)
c
c This subroutine computes a fermion wavefunction with the flowing-OUT
c fermion number.
c
c input:
c       real    p(0:3)         : four-momentum of fermion
c       real    fmass          : mass          of fermion
c       integer nhel = -1 or 1 : helicity      of fermion
c       integer nsf  = -1 or 1 : +1 for particle, -1 for anti-particle
c
c output:
c       complex fo(6)          : fermion wavefunction               <fo|
c     
      implicit none
      double complex fo(6),chi(2)
      double precision p(0:3),sf(2),sfomeg(2),omega(2),fmass,
     &     pp,pp3,sqp0p3,sqm
      integer nhel,nsf,nh,ip,im

      double precision rZero, rHalf, rTwo
      parameter( rZero = 0.0d0, rHalf = 0.5d0, rTwo = 2.0d0 )

#ifdef HELAS_CHECK
      double precision p2
      double precision epsi
      parameter( epsi = 2.0d-5 )
      integer stdo
      parameter( stdo = 6 )
#endif
c
#ifdef HELAS_CHECK
      pp = sqrt(p(1)**2+p(2)**2+p(3)**2)
      if ( abs(p(0))+pp.eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : p(0:3) in oxxxxx is zero momentum'
      endif
      if ( p(0).le.rZero ) then
         write(stdo,*)
     &        ' helas-error : p(0:3) in oxxxxx has non-positive energy'
         write(stdo,*)
     &        '         : p(0) = ',p(0)
      endif
      p2 = (p(0)-pp)*(p(0)+pp)
      if ( abs(p2-fmass**2).gt.p(0)**2*epsi ) then
         write(stdo,*)
     &        ' helas-error : p(0:3) in oxxxxx has inappropriate mass'
         write(stdo,*)
     &        '             : p**2 = ',p2,' : fmass**2 = ',fmass**2
      endif
      if ( abs(nhel).ne.1 ) then
         write(stdo,*) ' helas-error : nhel in oxxxxx is not -1,1'
         write(stdo,*) '             : nhel = ',nhel
      endif
      if ( abs(nsf).ne.1 ) then
         write(stdo,*) ' helas-error : nsf in oxxxxx is not -1,1'
         write(stdo,*) '             : nsf = ',nsf
      endif
#endif

      fo(5) = dcmplx(p(0),p(3))*nsf
      fo(6) = dcmplx(p(1),p(2))*nsf

      nh = nhel*nsf

      if ( fmass.ne.rZero ) then

         pp = min(p(0),dsqrt(p(1)**2+p(2)**2+p(3)**2))

         if ( pp.eq.rZero ) then
            
            sqm = dsqrt(fmass)
            ip = -((1+nh)/2)
            im =  (1-nh)/2
            
            fo(1) = im     * sqm
            fo(2) = ip*nsf * sqm
            fo(3) = im*nsf * sqm
            fo(4) = ip     * sqm
            
         else
            
            pp = min(p(0),dsqrt(p(1)**2+p(2)**2+p(3)**2))
            sf(1) = dble(1+nsf+(1-nsf)*nh)*rHalf
            sf(2) = dble(1+nsf-(1-nsf)*nh)*rHalf
            omega(1) = dsqrt(p(0)+pp)
            omega(2) = fmass/omega(1)
            ip = (3+nh)/2
            im = (3-nh)/2
            sfomeg(1) = sf(1)*omega(ip)
            sfomeg(2) = sf(2)*omega(im)
            pp3 = max(pp+p(3),rZero)
            chi(1) = dcmplx( dsqrt(pp3*rHalf/pp) )
            if ( pp3.eq.rZero ) then
               chi(2) = dcmplx(-nh )
            else
               chi(2) = dcmplx( nh*p(1) , -p(2) )/dsqrt(rTwo*pp*pp3)
            endif
            
            fo(1) = sfomeg(2)*chi(im)
            fo(2) = sfomeg(2)*chi(ip)
            fo(3) = sfomeg(1)*chi(im)
            fo(4) = sfomeg(1)*chi(ip)

         endif
         
      else
         
         sqp0p3 = dsqrt(max(p(0)+p(3),rZero))*nsf
         chi(1) = dcmplx( sqp0p3 )
         if ( sqp0p3.eq.rZero ) then
            chi(2) = dcmplx(-nhel )*dsqrt(rTwo*p(0))
         else
            chi(2) = dcmplx( nh*p(1), -p(2) )/sqp0p3
         endif
         if ( nh.eq.1 ) then
            fo(1) = chi(1)
            fo(2) = chi(2)
            fo(3) = dcmplx( rZero )
            fo(4) = dcmplx( rZero )
         else
            fo(1) = dcmplx( rZero )
            fo(2) = dcmplx( rZero )
            fo(3) = chi(2)
            fo(4) = chi(1)
         endif
         
      endif
c
      return
      end
