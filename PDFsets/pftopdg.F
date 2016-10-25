      subroutine pftopdg_cteq(x,q,pdf)

      implicit none

      real*8 x, q, pdf(-6:6),u,d,Ctq6Pdf, CT10PDF
      external Ctq6Pdf, CT10PDF
      integer i, cteqtype
      COMMON/cteqchoice/ cteqtype
c      print*,' WARNING: mu_f fixed to 91 GeV in pftopdg '
c      q = 91
      do i = -5,2
         if (cteqtype .eq. 6) then 
           pdf(i)=Ctq6Pdf (i, X, Q)
         else if (cteqtype .eq. 10) then 
           pdf(i)=CT10PDF (i, X, Q)
         endif
      enddo
c switch up and down quarks
      u = pdf(1)
      d = pdf(2)
      pdf(1) = d
      pdf(2) = u
      u = pdf(-1)
      d = pdf(-2)
      pdf(-1) = d
      pdf(-2) = u
      do i = 3,5
         pdf(i) = pdf(-i)
      enddo
      pdf(-6)=0
      pdf(6) =0
      do i=-5,5
         pdf(i) = pdf(i)*x
      enddo
      end

