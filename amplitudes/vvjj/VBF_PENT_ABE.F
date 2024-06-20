      subroutine PENT_VBF(psi2,psi1,p1,p2,pbos1,pbos2,pjext,
     $     extbos1,extbos2,jext,
     $     isig1max,isig3max,
     $     scale,divmax,
     $     mat21g,mat21,mat21gb,mat21b,WardTest)
C     the purpose of this routine is to compute the needed pieces for all
C     triangle and box diagrams

      implicit none

      integer divMax
      logical Ward1,WardTest
      integer isig1max,isig3max 
      double precision scale
      integer isig1,isig2,isig3,div,comp
      double precision p1(0:3),p2(0:3),pbos1(0:3),pbos2(0:3),pjext(0:3)
      double complex PSI1(2,-1:1),PSI2(2,-1:1),extbos1(6),extbos2(6)
      double complex mat21(8,-1:1,-1:1,6,0:2), 
     $     mat21g(3,-1:1,-1:1,6,0:2), 
     $     mat21gb(3,-1:1,-1:1,6,0:2), 
     $     mat21b(-1:1,-1:1,6,0:2)
      double complex jext(0:3,-1:1)
       real*8 dotrr

      external theta,dotrr
      real*8 theta,pi
      parameter( pi = 3.141592653589793d0)


      Ward1=.false.
      WardTest=.true.

      call SetPentMeToZero( mat21g(1,-1,-1,1,0),
     $     mat21(1,-1,-1,1,0), mat21gb(1,-1,-1,1,0),
     $     mat21b(-1,-1,1,0),divmax )

c$$$      print*,"vbf_pent"
c$$$      print*,"p1",(p1(i),i=0,3)
c$$$      print*,"pbos1",(pbos1(i),i=0,3)
c$$$      print*,"pbos2",(pbos2(i),i=0,3)
c$$$      print*,"pjext",(pjext(i),i=0,3)
c$$$      print*,"p2",(p2(i),i=0,3)
c$$$      print*,"sum",((p1(i)+pbos1(i)+pbos2(i)+pjext(i)+p2(i)),i=0,3)
c$$$      print*
c$$$
c$$$      do isig1 = -1,1,2
c$$$         print*,(jext(i,isig1),i=0,3)
c$$$      enddo


c     V1 V2 j
      do div = 0,divmax
         comp=1
         do isig1=-1,isig1max,2
            isig2 = isig1       ! kronecker delta_{isig1,isig2}
            do isig3 = -1,isig3max,2

               call penlineABETotal( p1,pbos1,pbos2,pjext,p2,
     $              psi2(1,isig2),psi1(1,isig1),
     $              extbos1,extbos2,jext(0,isig3),isig1,
     $              scale,0,1,
c     FC Here the comp variable, also the gauge check should be computed in the first run
     $              comp*4,comp,
     $              mat21g(1,isig1,isig3,1,div),
     $              mat21(1,isig1,isig3,1,div),
     $              mat21gb(1,isig1,isig3,1,div),
     $              mat21b(isig1,isig3,1,div),div )

               if(comp.eq.1) then
c                  print*, 'Div',Div
                  call pent_gaugeT( p1,pbos1,pbos2,pjext,p2,
     $                 psi2(1,isig2),psi1(1,isig1),
     $                 extbos1,extbos2,jext(0,isig3),isig1,
     $                 scale,mat21g(1,isig1,isig3,1,div),
     &                 mat21gb(1,isig1,isig3,1,div),Ward1,Div )


                  WardTest = WardTest.and.Ward1
                  if(.not.WardTest) then
                     call SetPentMeToZero( mat21g(1,-1,-1,1,0),
     $                    mat21(1,-1,-1,1,0), mat21gb(1,-1,-1,1,0),
     $                    mat21b(-1,-1,1,0),divmax )
                     return
                  endif
               endif     
               comp=-1
            enddo
         enddo
      enddo

c     V2 V1 j
      do div = 0,divmax
         comp=1
         do isig1=-1,isig1max,2
            isig2 = isig1       ! kronecker delta_{isig1,isig2}
            do isig3 = -1,isig3max,2

               call penlineABETotal( p1,pbos2,pbos1,pjext,p2,
     $              psi2(1,isig2),psi1(1,isig1),
     $              extbos2,extbos1,jext(0,isig3),isig1,
     $              scale,0,1,comp*4,comp,
     $              mat21g(1,isig1,isig3,2,div),
     $              mat21(1,isig1,isig3,2,div),
     $              mat21gb(1,isig1,isig3,2,div),
     $              mat21b(isig1,isig3,2,div),div )

               if(comp.eq.1) then
c                  print*, 'Div',Div
                  call pent_gaugeT( p1,pbos2,pbos1,pjext,p2,
     $                 psi2(1,isig2),psi1(1,isig1),
     $                 extbos2,extbos1,jext(0,isig3),isig1,
     $                 scale,mat21g(1,isig1,isig3,2,div),
     &                 mat21gb(1,isig1,isig3,2,div),Ward1,Div )

                  WardTest = WardTest.and.Ward1
                  if(.not.WardTest) then
                     call SetPentMeToZero( mat21g(1,-1,-1,1,0),
     $                    mat21(1,-1,-1,1,0), mat21gb(1,-1,-1,1,0),
     $                    mat21b(-1,-1,1,0),divmax )
                     return
                  endif
               endif     
               comp=-1
            enddo
         enddo
      enddo

c     V1 j V2
      do div = 0,divmax
         comp=1
         do isig1=-1,isig1max,2
            isig2 = isig1       ! kronecker delta_{isig1,isig2}
            do isig3 = -1,isig3max,2

               call penlineABETotal( p1,pbos1,pjext,pbos2,p2,
     $              psi2(1,isig2),psi1(1,isig1),
     $              extbos1,jext(0,isig3),extbos2,isig1,
     $              scale,0,1,
c     FC Here the comp variable, also the gauge check should be computed in the first run
     $              comp*4,comp,
     $              mat21g(1,isig1,isig3,3,div),
     $              mat21(1,isig1,isig3,3,div),
     $              mat21gb(1,isig1,isig3,3,div),
     $              mat21b(isig1,isig3,3,div),div )

               if(comp.eq.1) then
c                  print*, 'Div',Div
                  call pent_gaugeT( p1,pbos1,pjext,pbos2,p2,
     $                 psi2(1,isig2),psi1(1,isig1),
     $                 extbos1,jext(0,isig3),extbos2,isig1,
     $                 scale,mat21g(1,isig1,isig3,3,div),
     &                 mat21gb(1,isig1,isig3,3,div),Ward1,Div )

                  WardTest = WardTest.and.Ward1
                  if(.not.WardTest) then
                     call SetPentMeToZero( mat21g(1,-1,-1,1,0),
     $                    mat21(1,-1,-1,1,0), mat21gb(1,-1,-1,1,0),
     $                    mat21b(-1,-1,1,0),divmax )
                     return
                  endif
               endif     
               comp=-1
            enddo
         enddo
      enddo

c     V2 j V1
      do div = 0,divmax
         comp=1
         do isig1=-1,isig1max,2
            isig2 = isig1       ! kronecker delta_{isig1,isig2}
            do isig3 = -1,isig3max,2

               call penlineABETotal( p1,pbos2,pjext,pbos1,p2,
     $              psi2(1,isig2),psi1(1,isig1),
     $              extbos2,jext(0,isig3),extbos1,isig1,
     $              scale,0,1,comp*4,comp,
     $              mat21g(1,isig1,isig3,4,div),
     $              mat21(1,isig1,isig3,4,div),
     $              mat21gb(1,isig1,isig3,4,div),
     $              mat21b(isig1,isig3,4,div),div )

               if(comp.eq.1) then
c                  print*, 'Div',Div
                  call pent_gaugeT( p1,pbos2,pjext,pbos1,p2,
     $                 psi2(1,isig2),psi1(1,isig1),
     $                 extbos2,jext(0,isig3),extbos1,isig1,
     $                 scale,mat21g(1,isig1,isig3,4,div),
     &                 mat21gb(1,isig1,isig3,4,div),Ward1,Div )

                  WardTest = WardTest.and.Ward1
                  if(.not.WardTest) then
                     call SetPentMeToZero( mat21g(1,-1,-1,1,0),
     $                    mat21(1,-1,-1,1,0), mat21gb(1,-1,-1,1,0),
     $                    mat21b(-1,-1,1,0),divmax )
                     return
                  endif
               endif     
               comp=-1
            enddo
         enddo
      enddo

c     j V1 V2
      do div = 0,divmax
         comp=1
         do isig1=-1,isig1max,2
            isig2 = isig1       ! kronecker delta_{isig1,isig2}
            do isig3 = -1,isig3max,2

               call penlineABETotal( p1,pjext,pbos1,pbos2,p2,
     $              psi2(1,isig2),psi1(1,isig1),
     $              jext(0,isig3),extbos1,extbos2,isig1,
     $              scale,0,1,comp*4,comp,
     $              mat21g(1,isig1,isig3,5,div),
     $              mat21(1,isig1,isig3,5,div),
     $              mat21gb(1,isig1,isig3,5,div),
     $              mat21b(isig1,isig3,5,div),div )

               if(comp.eq.1) then
c                  print*, 'Div',Div
                  call pent_gaugeT( p1,pjext,pbos1,pbos2,p2,
     $                 psi2(1,isig2),psi1(1,isig1),
     $                 jext(0,isig3),extbos1,extbos2,isig1,
     $                 scale,mat21g(1,isig1,isig3,5,div),
     &                 mat21gb(1,isig1,isig3,5,div),Ward1,Div )

                  WardTest = WardTest.and.Ward1
                  if(.not.WardTest) then
                     call SetPentMeToZero( mat21g(1,-1,-1,1,0),
     $                    mat21(1,-1,-1,1,0), mat21gb(1,-1,-1,1,0),
     $                    mat21b(-1,-1,1,0),divmax )
                     return
                  endif
               endif     
               comp=-1
            enddo
         enddo
      enddo

c     j V2 V1
      do div = 0,divmax
         comp=1
         do isig1=-1,isig1max,2
            isig2 = isig1       ! kronecker delta_{isig1,isig2}
            do isig3 = -1,isig3max,2

               call penlineABETotal( p1,pjext,pbos2,pbos1,p2,
     $              psi2(1,isig2),psi1(1,isig1),
     $              jext(0,isig3),extbos2,extbos1,isig1,
     $              scale,0,1,comp*4,comp,
     $              mat21g(1,isig1,isig3,6,div),
     $              mat21(1,isig1,isig3,6,div),
     $              mat21gb(1,isig1,isig3,6,div),
     $              mat21b(isig1,isig3,6,div),div )

               if(comp.eq.1) then
c                  print*, 'Div',Div
                  call pent_gaugeT( p1,pjext,pbos2,pbos1,p2,
     $                 psi2(1,isig2),psi1(1,isig1),
     $                 jext(0,isig3),extbos2,extbos1,isig1,
     $                 scale,mat21g(1,isig1,isig3,6,div),
     &                 mat21gb(1,isig1,isig3,6,div),Ward1,Div)

                  WardTest = WardTest.and.Ward1
                  if(.not.WardTest) then
                     call SetPentMeToZero( mat21g(1,-1,-1,1,0),
     $                    mat21(1,-1,-1,1,0), mat21gb(1,-1,-1,1,0),
     $                    mat21b(-1,-1,1,0),divmax )
                     return
                  endif
               endif     
               comp=-1
            enddo
         enddo
      enddo

      end

c-----------------------------------------------------------------------------

      subroutine SetPentMeToZero(mat21,mat21g,mat21gb,mat21b,divmax)
         implicit none

      double complex mat21(8,-1:1,-1:1,2,0:2), 
     $     mat21g(3,-1:1,-1:1,2,0:2), 
     $     mat21gb(3,-1:1,-1:1,2,0:2), 
     $     mat21b(-1:1,-1:1,2,0:2)
      integer isig1, isig3, divmax

      do isig1=-1,1,2
         do isig3 = -1,1,2
            mat21(1:8,isig1,isig3,1:2,0:divmax)=0
            mat21g(1:3,isig1,isig3,1:2,0:divmax)=0
            mat21gb(1:3,isig1,isig3,1:2,0:divmax)=0
            mat21b(isig1,isig3,1:2,0:divmax)=0
         enddo            !isig3
      enddo               !isig1

      end
