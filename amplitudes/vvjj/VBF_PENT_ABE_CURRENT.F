      subroutine PENT_VBF_CURRENT(p1,pbos1,pbos2,pjext,p2,
     $      psi2,psi1, 
     $     extbos1,extbos2,jvirt,jborn,
     $     isig1max,
     $     comp, current_number,
     $     scale,divmax,
     $     mat21g,mat21gb,WardTest)
C     the purpose of this routine is to compute the needed pieces for all
C     triangle and box diagrams

      implicit none
      integer divMax
      logical Ward1,WardTest
      integer isig1max
      double precision scale
      integer isig1,isig2,div,comp,i
      double precision p1(0:3),p2(0:3),pjext(0:3),pbos1(0:3),pbos2(0:3)
      double complex PSI1(2,-1:1),PSI2(2,-1:1),extbos1(6),extbos2(6)
      double complex mat21(8,-1:1,0:3,0:2), ! check dime
     $     mat21g(3,-1:1,0:3,0:2), ! check dimen
     $     mat21gb(3,-1:1,0:3,0:2), ! check dimen
     $     mat21b(-1:1,0:3,0:2)  ! check dimen
      real*8 dotrr
      external theta,dotrr
      real*8 theta
cdc      complex*16 ratioFact,lns,cteEps1
c
      real*8 delta
      external delta
      double complex V1(0:3)
      double complex jvirt(0:3,-1:1),jborn(0:3,-1:1)
      integer current_number, comp1
      double precision pi2o3,cvirtc,pi
      parameter (pi=3.141592653589793d0)
      parameter (pi2o3=pi**2/3d0)


      if(comp.eq.1) then
      Ward1=.false.
      endif
      WardTest=.true.
      comp1=comp

cfc      call SetPentMeToZero1( mat21g,
cfc     $     mat21, mat21gb,
cfc     $     mat21b,divmax )

cfc      print*,"vbf_pent"
cfc      print*,"p1",(p1(i),i=0,3)
cfc      print*,"pjext",(pjext(i),i=0,3)
cfc      print*,"pbos1",(pbos1(i),i=0,3)
cfc      print*,"pbos2",(pbos2(i),i=0,3)
cfc      print*,"p2",(p2(i),i=0,3)
cfc      print*,"sum",((p1(i)+pbos1(i)+pbos2(i)+pjext(i)+p2(i)),i=0,3)


c$$$      print*
c$$$
c$$$      do isig1 = -1,1,2
c$$$         print*,(jext(i,isig1),i=0,3)
c$$$      enddo


c     V1 V2 j
      do div = 0,divmax

      if(div.eq.0) then
         cvirtc=pi2o3-7d0
      else
         cvirtc=0d0
         comp1=1
          Ward1=.false.
      endif

      do i=0,3
        v1(0)=delta(i,0)
        v1(1)=-delta(i,1)
        v1(2)=-delta(i,2)
        v1(3)=-delta(i,3)

         do isig1=-1,isig1max,2
            isig2 = isig1       ! kronecker delta_{isig1,isig2}

         if(current_number.eq.3) then

               call penlineABETotal( p1,pbos1,pbos2,pjext,p2,
     $              psi2(1,isig2),psi1(1,isig1),
     $              extbos1,extbos2,V1,isig1,
     $              scale,0,1,
c     FC Here the comp variable, also the gauge check should be computed in the first run
     $              comp1*4,comp1,
     $              mat21g(1,isig1,i,div),
     $              mat21(1,isig1,i,div),
     $              mat21gb(1,isig1,i,div),
     $              mat21b(isig1,i,div),div )

               if(comp1.eq.1) then
c                  print*, 'Div',Div
                  call pent_gaugeT( p1,pbos1,pbos2,pjext,p2,
     $                 psi2(1,isig2),psi1(1,isig1),
     $                 extbos1,extbos2,V1,isig1,
     $                 scale,mat21g(1,isig1,0,div),
     &                 mat21gb(1,isig1,0,div),Ward1,Div )


                  WardTest = WardTest.and.Ward1
                  if(.not.WardTest) then
                     call SetPentMeToZero1( mat21g,
     $                    mat21, mat21gb,
     $                    mat21b,divmax )
                     return
                  endif
               endif     
               comp1=-1

      elseif(current_number.eq.2) then

c     V1 j V2
               call penlineABETotal( p1,pbos1,pjext,pbos2,p2,
     $              psi2(1,isig2),psi1(1,isig1),
     $              extbos1,V1,extbos2,isig1,
     $              scale,0,1,
c     FC Here the comp1 variable, also the gauge check should be comp1uted in the first run
     $              comp1*4,comp1,
     $              mat21g(1,isig1,i,div),
     $              mat21(1,isig1,i,div),
     $              mat21gb(1,isig1,i,div),
     $              mat21b(isig1,i,div),div )

               if(comp1.eq.1) then
c                  print*, 'Div',Div
                  call pent_gaugeT( p1,pbos1,pjext,pbos2,p2,
     $                 psi2(1,isig2),psi1(1,isig1),
     $                 extbos1,V1,extbos2,isig1,
     $                 scale,mat21g(1,isig1,0,div),
     &                 mat21gb(1,isig1,0,div),Ward1,Div )

                  WardTest = WardTest.and.Ward1
                  if(.not.WardTest) then
                     call SetPentMeToZero1( mat21g,
     $                    mat21, mat21gb,
     $                    mat21b,divmax )
                     return
                  endif
               endif     
               comp1=-1

      elseif(current_number.eq.1) then
c     j V1 V2
cfc         print*, "HERE"
cfc
cfc      print*, "p1", p1
cfc      print*, "pjext", pjext
cfc      print*, "pbos1",pbos1
cfc      print*, "pbos2",pbos2
cfc      print*, "p2", p2
cfc      print*, "isig1max", isig1max
cfc      print*, "comp", comp
cfc      print*, "current_number", current_number
cfc      print*, "scale", scale
cfc      print*, "divmax", divmax


               call penlineABETotal(p1,pjext,pbos1,pbos2,p2,
     $              psi2(1,isig2),psi1(1,isig1),
     $              V1,extbos1,extbos2,isig1,
     $              scale,0,1,comp1*4,comp1,
     $              mat21g(1,isig1,i,div),
     $              mat21(1,isig1,i,div),
     $              mat21gb(1,isig1,i,div),
     $              mat21b(isig1,i,div),div )

               if(comp1.eq.1) then
c                  print*, 'Div',Div
                  call pent_gaugeT( p1,pjext,pbos1,pbos2,p2,
     $                 psi2(1,isig2),psi1(1,isig1),
     $                 V1,extbos1,extbos2,isig1,
     $                 scale,mat21g(1,isig1,0,div),
     &                 mat21gb(1,isig1,0,div),Ward1,Div )

c                  print*,"Ward1", Ward1
                  WardTest = WardTest.and.Ward1
                  if(.not.WardTest) then
                     call SetPentMeToZero1( mat21g,
     $                    mat21, mat21gb,
     $                    mat21b,divmax )
                     return
                  endif
               endif     
               comp1=-1
      else
         print*, "Wrong current in VBF_PENT_ABE_CURRENT"
      endif
      enddo
      enddo ! MU of the current
      enddo ! DIV

c NOTE THE CHANGE OF SIGN TO AGREE with CARLO routine
      do isig1=-1,1,2
        do i=0,3
           jvirt(i,isig1)= mat21(1,isig1,i,divmax)
     $            -mat21b(isig1,i,divmax)*(cvirtc)
           jborn(i,isig1)= mat21b(isig1,i,divmax)
           jvirt(i,isig1)= -jvirt(i,isig1)
           jborn(i,isig1)= -jborn(i,isig1)
        enddo
      enddo 

cfcc assign momenta in 4,5 comp1onent
cfc      do isig1=-1,1,2
cfc           jvirt(4,isig1)= -DCMPLX(pjext(0),pjext(3))
cfc           jvirt(5,isig1)= -DCMPLX(pjext(1),pjext(2))
cfc           jborn(4,isig1)= -DCMPLX(pjext(0),pjext(3))
cfc           jborn(5,isig1)= -DCMPLX(pjext(1),pjext(2))
cfc      enddo


      end

c-----------------------------------------------------------------------------

      subroutine SetPentMeToZero1(mat21,mat21g,mat21gb,mat21b,divmax)

      double complex mat21(8,-1:1,0:3,0:2), 
     $     mat21g(3,-1:1,0:3,0:2), 
     $     mat21gb(3,-1:1,0:3,0:2), 
     $     mat21b(-1:1,0:3,0:2)
      integer divmax

      mat21(1:8,-1:1:2,0:3,0:divmax)=0
      mat21g(1:3,-1:1:2,0:3,0:divmax)=0
      mat21gb(1:3,-1:1:2,0:3,0:divmax)=0
      mat21b(-1:1:2,0:3,0:divmax)=0
      end
