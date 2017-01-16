      subroutine BOX_VBF(psi2,psi1,p1,p2,pbos1,pbos2,extbos,jext,
     $     isig1max,isig3max,
     $     scale,divmax,
     $     mat21g,mat21,mat21gb,mat21b,WardTest)
C     this routine computes all the different permutations for a boxline

      implicit none
      integer divMax
      logical Ward1,WardTest
      integer isig1max,isig3max 
      double precision scale
      integer isig1,isig2,isig3,div,comp
      double precision p1(0:3),p2(0:3),pbos1(0:3),pbos2(0:3)
      double complex PSI1(2,-1:1),PSI2(2,-1:1),extbos(1:5)
      double complex mat21(5,-1:1,-1:1,2,0:2), 
     $     mat21g(2,-1:1,-1:1,2,0:2), 
     $     mat21gb(2,-1:1,-1:1,2,0:2), 
     $     mat21b(-1:1,-1:1,2,0:2)
      double complex jext(0:3,-1:1)
      external dotrr
      real*8 dotrr

      Ward1=.false.
      WardTest=.true.

c      print*,"divmax",divmax

      call setmetozero( mat21g,
     $     mat21, mat21gb,
     $     mat21b, divmax )

c$$$      print*,"vbf_box"
c$$$      print*,"p1",(p1(i),i=0,3)
c$$$      print*,"pbos1",(pbos1(i),i=0,3)
c$$$      print*,"pbos2",(pbos2(i),i=0,3)
c$$$      print*,"p2",(p2(i),i=0,3)
c$$$      print*,"sum",((p1(i)+pbos1(i)+pbos2(i)+p2(i)),i=0,3)
c$$$      print*
c$$$      do isig3 = -1,1,2
c$$$      print*,"isig3, jext",isig3,(jext(i,isig3),i=0,3)
c$$$      enddo

c     V j
      do div = 0,divmax,1
c      print*,"div=",div
         comp=1
         do isig1=-1,isig1max,2
            isig2 = isig1       ! kronecker delta_{isig1,isig2}
            do isig3 = -1,isig3max,2

               call boxlineABETotal( p1,pbos1,pbos2,p2,
     $              psi2(1,isig2),psi1(1,isig1),
     $              extbos,jext(0,isig3),isig1,
     $              scale,0,1,
c     FC
c     FC Here the comp variable, also the gauge check should be computed in the first run
     $              comp*3,comp,
c     FC
     $              mat21g(1,isig1,isig3,1,div),
     $              mat21(1,isig1,isig3,1,div),
     $              mat21gb(1,isig1,isig3,1,div),
     $              mat21b(isig1,isig3,1,div),div )

               if(comp.eq.1) then
                  call box_gaugeT(p1,pbos1,pbos2,p2,scale,
     &                 mat21g(1,-1,-1,1,div),mat21gb(1,-1,-1,1,div),
     &                 Ward1,Div)
c                  print*,"ward1",Ward1
                  WardTest = WardTest.and.Ward1
                  if(.not.WardTest) then
                     call setmetozero( mat21g(1,-1,-1,1,0),
     $                    mat21(1,-1,-1,1,0), mat21gb(1,-1,-1,1,0),
     $                    mat21b(-1,-1,1,0),divmax )
                     return

                  endif
               endif   

               comp=-1

            enddo
         enddo
      enddo

c     j V
      do div = 0,divmax,1
c      print*,"div=",div
         comp=1
         do isig1=-1,isig1max,2
            isig2 = isig1       ! kronecker delta_{isig1,isig2}
            do isig3 = -1,isig3max,2

               call boxlineABETotal( p1,pbos2,pbos1,p2,
     $              psi2(1,isig2),psi1(1,isig1),
     $              jext(0,isig3),extbos,isig1,
     $              scale,0,1,comp*3,comp,
     $              mat21g(1,isig1,isig3,2,div),
     $              mat21(1,isig1,isig3,2,div),
     $              mat21gb(1,isig1,isig3,2,div),
     $              mat21b(isig1,isig3,2,div),div )
               
               if(comp.eq.1) then
                  call box_gaugeT(p1,pbos1,pbos2,p2,scale,
     &                 mat21g(1,-1,-1,1,div),mat21gb(1,-1,-1,1,div),
     &                 Ward1,Div)
                  WardTest = WardTest.and.Ward1
c                  print*,"ward2",Ward1
                  if(.not.WardTest) then
                     call setmetozero( mat21g(1,-1,-1,1,0),
     $                    mat21(1,-1,-1,1,0), mat21gb(1,-1,-1,1,0),
     $                    mat21b(-1,-1,1,0),divmax )
                     return
                  endif
               endif     

               comp=-1

            enddo
         enddo
      enddo

c$$$      if(divmax.eq.2) then
c$$$         div = 2
c$$$         j=1
c$$$         print*,"me in vbf_box"
c$$$         do i=1,2
c$$$            do isig1=-1,1,2
c$$$               do isig3 = -1,1,2
c$$$c$$$  do j=1,5
c$$$c$$$  print*,"mat21",mat21(j,isig1,isig3,i,div)
c$$$c$$$  enddo
c$$$c$$$  do j=1,2
c$$$c$$$  print*,"mat21g",mat21g(j,isig1,isig3,i,div)
c$$$c$$$  print*,"mat21gb",mat21gb(j,isig1,isig3,i,div)
c$$$c$$$  enddo
c$$$                  print*,"i,isig1,isig3,mat21",i,isig1,isig3,mat21(j,isig1,isig3,i,div)
c$$$c$$$                  print*,"mat21g",mat21g(j,isig1,isig3,i,div)
c$$$c$$$                  print*,"mat21gb",mat21gb(j,isig1,isig3,i,div)
c$$$c$$$                  print*,"mat21b",mat21b(isig1,isig3,i,div)
c$$$               enddo            !isig3
c$$$            enddo               !isig1
c$$$         enddo                  !i
c$$$      endif                     !divmax

      end

c----------------------------------------------------------------------

      subroutine setmetozero(mat21g,mat21,mat21gb,mat21b,divmax)
         implicit none

      double complex mat21(5,-1:1,-1:1,2,0:2), 
     $     mat21g(2,-1:1,-1:1,2,0:2), 
     $     mat21gb(2,-1:1,-1:1,2,0:2), 
     $     mat21b(-1:1,-1:1,2,0:2)
      integer divmax

c      print*,"setting me to zero"
       mat21(1:5,-1:1:2,-1:1:2,1:2,0:divmax)=0
       mat21g(1:2,-1:1:2,-1:1:2,1:2,0:divmax)=0
       mat21gb(1:2,-1:1:2,-1:1:2,1:2,0:divmax)=0
       mat21b(-1:1:2,-1:1:2,1:2,0:divmax)=0

      end
