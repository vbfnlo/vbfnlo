c ------------------------------------------------------------------------
c ---- complex dilogarithm -----------------------------------------------
c ------------------------------------------------------------------------
       complex*16 function cdilog(z)                                    
       implicit none                                   
       complex*16 z,zl,coef,dilog1,u,caux                               
       real*8 pi,sign                                                   
       integer n,i 
       pi=3.141592653589793238462643d0                                  
       zl=z                                                             
       dilog1=dcmplx(pi**2/6.d0)                                        
       if(dreal(zl).eq.1.and.dimag(zl).eq.0.) then                      
          cdilog=dilog1                                                    
          return                                                           
       else if (cdabs(zl).lt.1.d-2) then   
          n=-40./dlog(cdabs(zl))                                           
          caux=(0.d0,0.d0)                                                 
          do i=1,n                                                         
             caux=caux+zl**i/dble(i**2)                                       
          enddo                                                           
          cdilog=caux                                                      
          return                                                           
       else if(cdabs(zl).lt.1.) then                                    
          sign=1.d0                                                       
          coef=dcmplx(dble(0.))                                           
       else                                                            
          coef=-cdlog(-zl)**2/2.d0-dilog1                                 
          sign=-1.d0                                                      
          zl=1.d0/zl                                                      
       endif                                                          
       if(dreal(zl).gt.0.5) then                   
          coef=coef+sign*(dilog1-cdlog(zl)*cdlog(1.d0-zl))                
          sign=-sign                                                      
          zl=1.d0-zl                                                      
       else   
       endif  
       u=-cdlog(1.d0-zl)                                               
       cdilog=u-u**2/4.d0+u**3/36.d0-u**5/3600.d0+u**7/211680.d0       
     &  -u**9/10886400.d0+u**11*5.d0/2634508800.d0                     
       cdilog=cdilog-u**13*691.d0/2730.d0/6227020800.d0                
       cdilog=cdilog+u**15*7.d0/6.d0/1.307674368d12                    
       cdilog=cdilog-u**17*3617.d0/510.d0/3.5568742810d14              
       cdilog=cdilog+u**19*43867.d0/798.d0/1.2164510041d17              
       cdilog=cdilog-u**21*174611.d0/330.d0/5.1090942172d19            
       cdilog=sign*cdilog+coef                                         
       return                                                          
       end                                                             
