      subroutine finitecollqqREM(xa,x,patilde,pb,pin,mufac,quark,result,pdfchange)

      implicit none

#include "global.inc"

      integer quark
      integer pdfchange

      double precision xa,x,patilde(0:3),pb(0:3),pin(0:3),result,dotrr,
     &   xxa,pdfqxa,pdfgxa,pdfqxxa,pdfgxxa,logomx, logomxxa,logx,omx,omxxa,
     &   Logmupapb2,logmupapi2,NF,mufacsq,mufac,pdftmp(-6:6),one,inv2

      parameter (nf=5d0,one=1d0,inv2=1d0/2d0)

      external dotrr

      mufacsq=mufac**2
      xxa=x*xa
      omxxa=one-xxa
      omx=one-x
!       print*, xa,xxa, " in rems"

      logx=Log(x)
!       if (x.ge.0.999995) logx=0d0

      logomx=Log(omx)
!       if (x.le.0.000005) logomx=0d0

      logomxxa=Log(omxxa)
!       if (xxa.le.0.000005) logomxxa=0d0

      
      logmupapb2=Log(mufacsq*inv2/dotrr(patilde,pb))
      logmupapi2=Log(mufacsq*inv2/dotrr(patilde,pin))

      call pdfproton( xa, mufac, pdftmp(-6) )
      pdfqxa=pdftmp(pdfchange*quark)/xa
      pdfgxa=pdftmp(0)/xa
      call pdfproton( xxa, mufac, pdftmp(-6) )   
      pdfqxxa=pdftmp(pdfchange*quark)/xxa
      pdfgxxa=pdftmp(0)/xxa
  
      result=(logmupapb2*pdfgxa)/(16.*omx) - (9.*logmupapi2*pdfgxa)/(16.*omx) + (7.*Logomx*pdfgxa)/(16.*omx) - 
     -  (Logx*pdfgxa)/(2.*omx) - (11.*pdfqxa)/(4.*omx**2) + (logmupapb2*pdfqxa)/(3.*omx**2) - 
     -  (3.*logmupapi2*pdfqxa)/omx**2 + (7.*Logomx*pdfqxa)/(3.*omx**2) - (4.*Logx*pdfqxa)/omx**2 + 
     -  (Nf*pdfqxa)/(6.*omx**2) + (4.*pdfqxa)/(3.*omx) - (logmupapb2*pdfqxa)/(6.*omx) + 
     -  (3.*logmupapi2*pdfqxa)/(2.*omx) - (7.*Logomx*pdfqxa)/(6.*omx) + (8.*Logx*pdfqxa)/(3.*omx) + (pdfgxa*x)/omx - 
     -  (logmupapb2*pdfgxa*x)/(8.*omx) + (9.*logmupapi2*pdfgxa*x)/(8.*omx) - (7.*Logomx*pdfgxa*x)/(8.*omx) + 
     -  (Logx*pdfgxa*x)/omx - (4.*pdfqxa*x)/(3.*omx) - (logmupapb2*pdfqxa*x)/(6.*omx) + 
     -  (3.*logmupapi2*pdfqxa*x)/(2.*omx) - (7.*Logomx*pdfqxa*x)/(6.*omx) + (8.*Logx*pdfqxa*x)/(3.*omx) + 
     -  (11.*pdfqxxa*x)/(4.*omx**2) - (logmupapb2*pdfqxxa*x)/(3.*omx**2) + (3.*logmupapi2*pdfqxxa*x)/omx**2 - 
     -  (7.*Logomx*pdfqxxa*x)/(3.*omx**2) - (Nf*pdfqxxa*x)/(6.*omx**2) - (113.*pdfqxxa*x)/(12.*omx*omxxa) + 
     -  (logmupapb2*pdfqxxa*x)/(4.*omx*omxxa) - (9.*logmupapi2*pdfqxxa*x)/(4.*omx*omxxa) - 
     -  (11.*Logomxxa*pdfqxxa*x)/(4.*omx*omxxa) + (logmupapb2*Logomxxa*pdfqxxa*x)/(3.*omx*omxxa) - 
     -  (3.*logmupapi2*Logomxxa*pdfqxxa*x)/(omx*omxxa) + (7.*Logomxxa**2*pdfqxxa*x)/(6.*omx*omxxa) + 
     -  (Nf*pdfqxxa*x)/(6.*omx*omxxa) + (Logomxxa*Nf*pdfqxxa*x)/(6.*omx*omxxa) + 
     -  (17.*pdfqxxa*Pi**2*x)/(18.*omx*omxxa) - (pdfgxa*x**2)/omx + (logmupapb2*pdfgxa*x**2)/(8.*omx) - 
     -  (9.*logmupapi2*pdfgxa*x**2)/(8.*omx) + (7.*Logomx*pdfgxa*x**2)/(8.*omx) - (Logx*pdfgxa*x**2)/omx + 
     -  (4.*Logx*pdfqxa*x**2)/(3.*omx**2)


      end


      subroutine finitecollgqREM(xa,x,patilde,pb,pin,mufac,result,pdfchange)

      implicit none

#include "global.inc"

      integer quark
      integer pdfchange

      double precision xa,x,patilde(0:3),pb(0:3),pin(0:3),result,dotrr,
     &   xxa,pdfqxa,pdfgxa,pdfqxxa,pdfgxxa,logomx, logomxxa,logx,omx,omxxa,
     &   Logmupapb2,logmupapi2,NF,mufacsq,mufac,pdftmp(-6:6),one,inv2,zero
      parameter (nf=5d0,one=1d0,inv2=1d0/2d0,zero=0d0)

      external dotrr

      mufacsq=mufac**2
      xxa=x*xa
      omxxa=one-xxa
      omx=one-x

      logx=Log(x)
!       if (x.ge.0.999995) logx=0d0

      logomx=Log(omx)
!       if (x.le.0.000005) logomx=0d0

      logomxxa=Log(omxxa)
!       if (xxa.le.0.000005) logomxxa=0d0

      logmupapb2=Log(mufacsq*inv2/dotrr(patilde,pb))
      logmupapi2=Log(mufacsq*inv2/dotrr(patilde,pin))

      pdfqxa=zero
      pdfqxxa=zero
      
      call pdfproton( xa, mufac, pdftmp(-6) )
      pdfgxa=pdftmp(0)/xa

      do quark=1,4  ! sum over quark gluon splittings
      pdfqxa=pdfqxa+pdftmp(pdfchange*quark)/xa
      enddo
      do quark=-4,-1
      pdfqxa=pdfqxa+pdftmp(pdfchange*quark)/xa
      enddo

      call pdfproton( xxa, mufac, pdftmp(-6) )    
      pdfgxxa=pdftmp(0)/xxa
      do quark=1,4  ! sum over quark gluon splittings
      pdfqxxa=pdfqxxa+pdftmp(pdfchange*quark)/xxa
      enddo
      do quark=-4,-1
      pdfqxxa=pdfqxxa+pdftmp(pdfchange*quark)/xxa
      enddo

      result=(-9.*pdfgxa)/(4.*omx**2) - (3.*logmupapb2*pdfgxa)/omx**2 - (3.*logmupapi2*pdfgxa)/omx**2 + 
     -  (9.*Logomx*pdfgxa)/omx**2 - (6.*Logx*pdfgxa)/omx**2 + (6.*logmupapb2*pdfgxa)/omx + (6.*logmupapi2*pdfgxa)/omx - 
     -  (18.*Logomx*pdfgxa)/omx + (12.*Logx*pdfgxa)/omx + (4.*logmupapb2*pdfqxa)/(3.*omx) + 
     -  (4.*logmupapi2*pdfqxa)/(3.*omx) - (4.*Logomx*pdfqxa)/omx + (8.*Logx*pdfqxa)/(3.*omx) - 
     -  (3.*logmupapb2*pdfgxa)/(omx*x) - (3.*logmupapi2*pdfgxa)/(omx*x) + (9.*Logomx*pdfgxa)/(omx*x) - 
     -  (6.*Logx*pdfgxa)/(omx*x) - (4.*logmupapb2*pdfqxa)/(3.*omx*x) - (4.*logmupapi2*pdfqxa)/(3.*omx*x) + 
     -  (4.*Logomx*pdfqxa)/(omx*x) - (8.*Logx*pdfqxa)/(3.*omx*x) - (3.*logmupapb2*pdfgxa*x)/omx - 
     -  (3.*logmupapi2*pdfgxa*x)/omx + (9.*Logomx*pdfgxa*x)/omx - (6.*Logx*pdfgxa*x)/omx + (9.*pdfgxxa*x)/(4.*omx**2) + 
     -  (3.*logmupapb2*pdfgxxa*x)/omx**2 + (3.*logmupapi2*pdfgxxa*x)/omx**2 - (9.*Logomx*pdfgxxa*x)/omx**2 - 
     -  (227.*pdfgxxa*x)/(12.*omx*omxxa) - (11.*logmupapb2*pdfgxxa*x)/(4.*omx*omxxa) - 
     -  (11.*logmupapi2*pdfgxxa*x)/(4.*omx*omxxa) - (9.*Logomxxa*pdfgxxa*x)/(4.*omx*omxxa) - 
     -  (3.*logmupapb2*Logomxxa*pdfgxxa*x)/(omx*omxxa) - (3.*logmupapi2*Logomxxa*pdfgxxa*x)/(omx*omxxa) + 
     -  (9.*Logomxxa**2*pdfgxxa*x)/(2.*omx*omxxa) + (8.*NF*pdfgxxa*x)/(9.*omx*omxxa) + 
     -  (logmupapb2*NF*pdfgxxa*x)/(6.*omx*omxxa) + (logmupapi2*NF*pdfgxxa*x)/(6.*omx*omxxa) + (4.*pdfqxa*x)/(3.*omx) - 
     -  (2.*logmupapb2*pdfqxa*x)/(3.*omx) - (2.*logmupapi2*pdfqxa*x)/(3.*omx) + (2.*Logomx*pdfqxa*x)/omx - 
     -  (4.*Logx*pdfqxa*x)/(3.*omx) + (3.*pdfgxxa*Pi**2*x)/(2.*omx*omxxa) + (3.*logmupapb2*pdfgxa*x**2)/omx + 
     -  (3.*logmupapi2*pdfgxa*x**2)/omx - (9.*Logomx*pdfgxa*x**2)/omx +
     -  (6.*Logx*pdfgxa*x**2)/omx

      end


      subroutine finitecollqgREM(xa,x,patilde,pb,pin,mufac,quark,result,pdfchange)

      implicit none

#include "global.inc"


      integer quark
      integer pdfchange

      double precision xa,x,patilde(0:3),pb(0:3),pin(0:3),result,dotrr,
     &   xxa,pdfqxa,pdfgxa,pdfqxxa,pdfgxxa,logomx, logomxxa,logx,omx,omxxa,
     &   Logmupapb2,logmupapi2,NF,mufacsq,mufac,pdftmp(-6:6),one,inv2,zero
      parameter (nf=5d0,one=1d0,inv2=1d0/2d0,zero=0d0)

      external dotrr

      mufacsq=mufac**2
      xxa=x*xa

      omxxa=one-xxa
      omx=one-x

      logx=Log(x)
!       if (x.ge.0.999995) logx=0d0

      logomx=Log(omx)
!       if (x.le.0.000005) logomx=0d0

      logomxxa=Log(omxxa)
!       if (xxa.le.0.000005) logomxxa=0d0

      logmupapb2=Log(mufacsq*inv2/dotrr(patilde,pb))
      logmupapi2=Log(mufacsq*inv2/dotrr(patilde,pin))

      pdfqxa=zero
      pdfqxxa=zero

      call pdfproton( xa, mufac, pdftmp(-6) )
      pdfgxa=pdftmp(0)/xa
      pdfqxa=pdftmp(pdfchange*quark)/xa
         
      call pdfproton( xxa, mufac, pdftmp(-6) )    
      pdfgxxa=pdftmp(0)/xxa
      pdfqxxa=pdftmp(pdfchange*quark)/xxa


      result=(-9.*logmupapb2*pdfgxa)/(16.*omx) + (logmupapi2*pdfgxa)/(16.*omx) + (17.*Logomx*pdfgxa)/(16.*omx) - 
     -  (Logx*pdfgxa)/(2.*omx) + (4.*pdfqxa)/3. + pdfqxa/(4.*omx**2) - (3.*logmupapb2*pdfqxa)/(2.*omx**2) + 
     -  (logmupapi2*pdfqxa)/(6.*omx**2) + (17.*Logomx*pdfqxa)/(6.*omx**2) - (4.*Logx*pdfqxa)/(3.*omx**2) + pdfgxa*x + 
     -  (9.*logmupapb2*pdfgxa*x)/8. - (logmupapi2*pdfgxa*x)/8. - (17.*Logomx*pdfgxa*x)/8. + Logx*pdfgxa*x - 
     -  (pdfqxxa*x)/(4.*omx**2) + (3.*logmupapb2*pdfqxxa*x)/omx**2 - (logmupapi2*pdfqxxa*x)/(3.*omx**2) - 
     -  (17.*Logomx*pdfqxxa*x)/(3.*omx**2) - (77.*pdfqxxa*x)/(12.*omx*omxxa) - (9.*logmupapb2*pdfqxxa*x)/(4.*omx*omxxa) + 
     -  (logmupapi2*pdfqxxa*x)/(4.*omx*omxxa) + (Logomxxa*pdfqxxa*x)/(4.*omx*omxxa) - 
     -  (3.*logmupapb2*Logomxxa*pdfqxxa*x)/(omx*omxxa) + (logmupapi2*Logomxxa*pdfqxxa*x)/(3.*omx*omxxa) + 
     -  (17.*Logomxxa**2*pdfqxxa*x)/(6.*omx*omxxa) + (7.*pdfqxxa*Pi**2*x)/(18.*omx*omxxa) - 
     -  (3.*logmupapb2*pdfqxa*x**2)/(2.*omx**2) + (logmupapi2*pdfqxa*x**2)/(6.*omx**2) + 
     -  (17.*Logomx*pdfqxa*x**2)/(6.*omx**2) - (4.*Logx*pdfqxa*x**2)/(3.*omx**2)


      end
