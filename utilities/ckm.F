!
! Includes routines used for treating the ckm correctly
!

! takes existing pdfs and multiplies ckm matrix ontop of them such
! that just by using those pdfs the ckm is taken into account
!
! pdfckm(:,1:2) only have up-type such that double-counting is
! avoided, pdfckm(:,3:4) also has down-types (e.g. for gluon+quark
! cases)
      subroutine pdfckmtransform(pdf,pdfckm)

         double precision pdf(-6:6,2), pdfckm(-6:6,2)
         integer j

         do j=1,2
            ! loop over pdfs of the two protons
            call pdfckmtransformsingle(pdf(-6,j),pdfckm(-6,j))
         enddo

      end subroutine pdfckmtransform


      subroutine pdfckmtransformsingle(pdf,pdfckm)

#include "VBFNLO/utilities/ckm.inc"
      double precision pdf(-6:6), pdfckm(-6:6)
      integer i
      integer u,d,ku,kd,k
      logical ldebug
      parameter (ldebug = .false.)

      ! set gluonpdf zero, such that gluons are not double counted
      ! in general pdfckm should not be used for gluons at all
      pdfckm(0)=0d0

      !disable certain flavors
      !pdf( 1) = 0d0 
      !pdf(-1) = 0d0 
      !pdf( 2) = 0d0 
      !pdf(-2) = 0d0 
      !pdf( 3) = 0d0 
      !pdf(-3) = 0d0 
      !pdf( 4) = 0d0 
      !pdf(-4) = 0d0 

      ! loop over all generations and modify the pdfs
      ! such that pdf u -> V_ckm * pdf(u,c,t)
      ! therefore pdf(u)*pdfckm(d) gives all the couplings ud,us,ub with right ckm factors
      do i=1,3
         ! loop over the 3 generations
         ! here one could use the information of PDFsets/pftopdg.F to save some steps
         u=2*i
            d=2*i-1
            pdfckm(u) = 0d0 
            pdfckm(-u) = 0d0
            pdfckm(d) = 0d0
            pdfckm(-d) = 0d0

            ! this could be simplifies assuming only cabbibo angle
            do k=1,3
               ku=2*k
               kd=2*k-1
               pdfckm(d)  = pdfckm( d)  +pdf(kd)  *vckm(k,i)
               pdfckm(-d) = pdfckm(-d)  +pdf(-kd) *vckm(k,i)
               pdfckm(u)  = pdfckm( u)  +pdf(ku)  *vckm(i,k)
               pdfckm(-u) = pdfckm(-u)  +pdf(-ku) *vckm(i,k)
            enddo
         enddo

         if(ldebug) then
            print*, 'in pdfckmtransform'
            print*, pdf
            print*, 'ckmmod'
            print*, pdfckm
            print*, 'ckm'
            print*, vckm
         endif

      end subroutine pdfckmtransformsingle


      ! routine to return pdf values
      function getpdf(pdf,part1,part2)
         double precision pdf(-6:6,2)
         integer part1, part2
         double precision getpdf
         getpdf = pdf(part1,1)*pdf(part2,2)
         return
      end function getpdf

      ! routines to return pdf values, modified by the ckm matrix
      ! subscript 1 or 2 forces only one of the partons gets the ckm
      ! modification (the one that couples to a w)
      function getpdfckm(pdf,pdfckm,part1,part2)
         double precision pdf(-6:6,2), pdfckm(-6:6,2)
         integer part1, part2
         double precision getpdfckm
         ! factor 0.5 removes double counting
         getpdfckm = 0.5d0*pdf(part1,1)*pdfckm(part2,2) + 0.5d0*pdfckm(part1,1)*pdf(part2,2)
         return
      end function getpdfckm

      function getpdfckm1(pdf,pdfckm,part1,part2)
         double precision pdf(-6:6,2), pdfckm(-6:6,2)
         integer part1, part2
         double precision getpdfckm1
         getpdfckm1 = pdfckm(part1,1)*pdf(part2,2)
         return
      end function getpdfckm1

      function getpdfckm2(pdf,pdfckm,part1,part2)
         double precision pdf(-6:6,2), pdfckm(-6:6,2)
         integer part1, part2
         double precision getpdfckm2
         getpdfckm2 = pdf(part1,1)*pdfckm(part2,2)
         return
      end function getpdfckm2


