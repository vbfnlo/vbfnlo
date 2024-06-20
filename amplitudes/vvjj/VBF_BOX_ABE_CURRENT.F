      subroutine BOX_VBF_CURR(p1,pbos1,pbos2,p2,
     $              psi2,psi1,
     $     extbos,jvirt,jborn,
     $     isig1max,
     $     comp,current_number, 
     $     scale,divmax,
     $     mat21g,mat21gb,WardTest)
C     this routine computes all the different permutations for a boxline

      implicit none
      integer divMax
      logical Ward1,WardTest
      integer isig1max 
      double precision scale
      integer isig1,isig2,div,comp,i,comp1
      double precision p1(0:3),p2(0:3),pbos1(0:3),pbos2(0:3)
      double complex PSI1(2,-1:1),PSI2(2,-1:1),extbos(1:5)
      double complex mat21(5,-1:1,0:3,0:2), 
     $     mat21g(2,-1:1,0:3,0:2), 
     $     mat21gb(2,-1:1,0:3,0:2), 
     $     mat21b(-1:1,0:3,0:2)
      double complex V1(0:3)
      external dotrr
      real*8 dotrr
      Real*8 delta
      External  delta
c external current
      double complex jvirt(0:5,-1:1),jborn(0:5,-1:1)
      integer current_number
      double precision pi2o3,cvirtc,pi
      parameter (pi=3.141592653589793d0)
      parameter (pi2o3=pi**2/3d0)
c old notation from carlo
c                        mu           jext(q2)      for ivmu = 1 
c                     jext(q1)          mu          for ivmu = 2

      if(comp.eq.1) then
      Ward1=.false.
      endif
      WardTest=.true.
      comp1=comp
c      print*,"divmax",divmax

      call setmetozero1( mat21g,
     $     mat21, mat21gb,
     $     mat21b, divmax )

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     V j
      do div = 0,divmax,1
      if(div.eq.0) then
         cvirtc=pi2o3-7d0
      else
         cvirtc=0d0
         comp1=1
          Ward1=.false.
      endif

c      print*,"div=",div
c      print*, "comp1",comp1
c         comp1=1
      do i=0,3
        v1(0)=delta(i,0)
        v1(1)=-delta(i,1)
        v1(2)=-delta(i,2)
        v1(3)=-delta(i,3)

         do isig1=-1,isig1max,2
            isig2 = isig1       ! kronecker delta_{isig1,isig2}
         if(current_number.eq.2) then
cfc               print*, "curr2"
cfc                print*, "k1", p1
cfc                print*, "k2", p2
cfc                print*, "q1", pbos1
cfc                print*, "q2",pbos2
cfc                print*, "psi2(1,isig2)",psi2(:,isig2)
cfc                print*, "psi1(1,isig2)",psi1(:,isig2)
cfc                print*, "extbos", extbos
cfc                print*, "V1", V1
cfc                print*, "isig1",isig1
cfc                print*, "scale",scale
cfc                print*, "comp1",comp1
cfc                print*, "div",div
               call boxlineABETotal(p1,pbos1,pbos2,p2,
     $              psi2(1,isig2),psi1(1,isig1),
     $              extbos,V1,isig1,
     $              scale,0,1,
c     FC
c     FC Here the comp1 variable, also the gauge check should be comp1uted in the first run
     $              comp1*3,comp1,
c     FC
     $              mat21g(1,isig1,i,div),
     $              mat21(1,isig1,i,div),
     $              mat21gb(1,isig1,i,div),
     $              mat21b(isig1,i,div),div )
cfc           print*,"here"
cfc           print*, mat21g(1,isig1,i,div)
cfc           print*,         mat21(1,isig1,i,div)
cfc           print*,         mat21gb(1,isig1,i,div)
cfc           print*,         mat21b(isig1,i,div)
               

             elseif(current_number.eq.1) then
cfc               print*, "curr1"
cfc                print*, "k1", p1
cfc                print*, "k2", p2
cfc                print*, "q1", pbos2
cfc                print*, "q2",pbos1

               call boxlineABETotal(p1,pbos2,pbos1,p2,
     $              psi2(1,isig2),psi1(1,isig1),
     $              V1,extbos,isig1,
     $              scale,0,1,
c     FC
c     FC Here the comp1 variable, also the gauge check should be comp1uted in the first run
     $              comp1*3,comp1,
c     FC
     $              mat21g(1,isig1,i,div),
     $              mat21(1,isig1,i,div),
     $              mat21gb(1,isig1,i,div),
     $              mat21b(isig1,i,div),div )
           else
           print*, "Wrong current number in VBF_Box_abe_current.F"             
           stop
        endif

cfc       print*, "here2"

      if(comp1.eq.1) then
         if(current_number.eq.2) then
                  call box_gaugeT(p1,pbos1,pbos2,p2,scale,
     &                 mat21g(1,-1,0,div),mat21gb(1,-1,0,div),
     &                 Ward1,Div)
cfc                  print*,"ward1",Ward1
cfc                 print*, mat21g(1,isig1,i,div)
cfc                 print*,         mat21(1,isig1,i,div)
cfc                 print*,         mat21gb(1,isig1,i,div)
cfc                  stop
          else
                  call box_gaugeT(p1,pbos2,pbos1,p2,scale,
     &                 mat21g(1,-1,0,div),mat21gb(1,-1,0,div),
     &                 Ward1,Div)
cfc                  print*,"ward1",Ward1
         endif
cfc         print*, "here3"
         WardTest = WardTest.and.Ward1
cfc         print*, "WardTest", WardTest
           if(.not.WardTest) then
               call setmetozero1( mat21g,
     $               mat21, mat21gb,
     $               mat21b,divmax )
                  return
            endif
cfc          print*, "here4"
          endif   
cfc          print*, "here5"
          comp1=-1
cfc          print*, "here6"
            enddo
         enddo
      enddo
cfc       print*, "here7"
cfc TO do Think about poles
      do isig1=-1,1,2
        do i=0,3
           jvirt(i,isig1)= mat21(1,isig1,i,divmax)
     $            -mat21b(isig1,i,divmax)*(cvirtc)
           jborn(i,isig1)= mat21b(isig1,i,divmax)
cfc       IF (DIV.ne.0) then
cfc          print*, "box Div =", Div
cfc          print*, jvirt(i,isig1)/jborn(i,isig1)
cfc       else
cfc          print*, "box Div = 0"
cfc        endif
        enddo
      enddo 
cfc      stop
      IF(current_number.eq.2) then
c assign momenta in 4,5 comp1onent
      do isig1=-1,1,2
           jvirt(4,isig1)= -DCMPLX(pbos2(0),pbos2(3))
           jvirt(5,isig1)= -DCMPLX(pbos2(1),pbos2(2))
           jborn(4,isig1)= -DCMPLX(pbos2(0),pbos2(3))
           jborn(5,isig1)= -DCMPLX(pbos2(1),pbos2(2))
      enddo
      else
      do isig1=-1,1,2
           jvirt(4,isig1)= -DCMPLX(pbos2(0),pbos2(3))
           jvirt(5,isig1)= -DCMPLX(pbos2(1),pbos2(2))
           jborn(4,isig1)= -DCMPLX(pbos2(0),pbos2(3))
           jborn(5,isig1)= -DCMPLX(pbos2(1),pbos2(2))
      enddo
      endif
      end


c----------------------------------------------------------------------

      subroutine setmetozero1(mat21g,mat21,mat21gb,mat21b,divmax)

      double complex mat21(5,-1:1,0:3,0:2), 
     $     mat21g(2,-1:1,0:3,0:2), 
     $     mat21gb(2,-1:1,0:3,0:2), 
     $     mat21b(-1:1,0:3,0:2)
      integer divmax

      mat21(1:5,-1:1:2,0:3,0:divmax)=0
      mat21g(1:2,-1:1:2,0:3,0:divmax)=0
      mat21gb(1:2,-1:1:2,0:3,0:divmax)=0
      mat21b(-1:1:2,0:3,0:divmax)=0
      end
