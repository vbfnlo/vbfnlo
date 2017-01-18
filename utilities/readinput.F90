!*************************************************************************  
!     This file contains subroutines for reading input files. 
!     The input files are in the format: key = value 
!*************************************************************************  

module readinput
      integer, parameter :: maxlines=200
      character*250 lines(maxlines)
      integer n_lines
      integer iunit
      character*250 fname
      logical usedefaults,showinput
      
      character*250 inputpath, pdfpath

contains

      subroutine loadfile(filename,show_input)
!*************************************************************************  
!     Opens a textfile and reads the data in an internal buffer.
!     End-of-line comments, starting with ! or //, are removed.
!     Do only process one single file at a time.
!     DON"T FORGET TO CALL CLOSEFILE AFTERWARDS.
!*************************************************************************  
          use globalvars, only: lglobalprint, ldoblha
      implicit none

      character*(*) filename
     
#include "global.inc"
#include "BLHAhelper.inc"
      
      integer i,le,k,j
      character*250 line
      logical show_input

      character*250 GetInputPath
      character*250 path
      external GetInputPath
      

      iunit=10
      usedefaults=.false.
      showinput=show_input .and. lglobalprint
! add input file path to filename
      path=GetInputPath()
      path = trim(path)//"/"//filename
! Open file
      open (unit=iunit, &
     &      file=trim(path), &
     &      status="old",err=77)

      fname = filename

! read text lines
      n_lines=0
 50   read(iunit,"(a250)",end=99,err=98) lines(n_lines+1)
      n_lines=n_lines+1
      if (n_lines.lt.maxlines) GOTO 50
 98   CONTINUE
      print *," Read error on file",filename
      print *," Using default values."
      usedefaults=.true.
 99   CONTINUE
      n_lines=n_lines+1

! remove comments
      do i=1,n_lines
         line=lines(i)
         if (INDEX(line,"!").ne.0) then
            line=line(1:INDEX(line,"!")-1)
         endif
         if (INDEX(line,"//").ne.0) then
            line=line(1:INDEX(line,"//")-1)
         endif
! upper case everything in front of a "="
         le = index(line,"=")
         if (le.gt.0) then
            do j=1,le-1
               k=ichar(line(j:j))
               if(k.ge.97.and.k.le.122) then
                  k=ichar(line(j:j))-32   
                  line(j:j)=char(k)
               endif
            enddo
         endif
         lines(i)=line
      enddo
!c     if (showinput) write(*,"(A,A,A)"), " reading input from file ",filename," : "
      
      return
 77   CONTINUE
      if (lglobalprint) print *," Open error on file ",filename
      if (.not.ldoblha .and. lglobalprint) then
        print *," Using default values."
      endif
      usedefaults=.true.

      end subroutine

!*************************************************************************  
      subroutine loadfileKK(filename,show_input)
!*************************************************************************  
!     Opens a textfile and reads the data in an internal buffer.
!     End-of-line comments, starting with ! or //, are removed.
!     Do only process one single file at a time.
!     DON"T FORGET TO CALL CLOSEFILE AFTERWARDS.
! IMPORTANT NOTE:  This routine is only used when reading in the KK
!                  data that's just been calculated and stored to file.
!                  It's not a particularly beautiful solution, but it 
!                  works for now.  See wiki.
!*************************************************************************  
      implicit none

      character*(*) filename
     
      integer i,le,k,j
      character*250 line
      logical show_input

      character*250 GetInputPath
      character*250 path
      external GetInputPath
      

      iunit=10
      usedefaults=.false.
      showinput=show_input
! add input file path to filename
      path = filename
!      path=GetInputPath()
!      path = trim(path)//"/"//filename
! Open file
      open (unit=iunit, &
     &      file=trim(path), &
     &      status="old",err=77)

      fname = filename

! read text lines
      n_lines=0
 50   read(iunit,"(a250)",end=99,err=98) lines(n_lines+1)
      n_lines=n_lines+1
      if (n_lines.lt.maxlines) GOTO 50
 98   CONTINUE
      print *," Read error on file",filename
      print *," Using default values."
      usedefaults=.true.
 99   CONTINUE
      n_lines=n_lines+1

! remove comments
      do i=1,n_lines
         line=lines(i)
         if (INDEX(line,"!").ne.0) then
            line=line(1:INDEX(line,"!")-1)
         endif
         if (INDEX(line,"//").ne.0) then
            line=line(1:INDEX(line,"//")-1)
         endif
! upper case everything in front of a "="
         le = index(line,"=")
         if (le.gt.0) then
            do j=1,le-1
               k=ichar(line(j:j))
               if(k.ge.97.and.k.le.122) then
                  k=ichar(line(j:j))-32   
                  line(j:j)=char(k)
               endif
            enddo
         endif
         lines(i)=line
      enddo
!c     if (showinput) write(*,"(A,A,A)"), " reading input from file ",filename," : "
      
      return
 77   CONTINUE
      print *," Open error on file ",filename
      print *," Using default values."
      usedefaults=.true.

      return
      end subroutine

!*************************************************************************  
      subroutine closefile
!*************************************************************************  
!     Closes the current file.
!*************************************************************************  
      implicit none

      close(iunit)
      end subroutine

!*************************************************************************  
      SUBROUTINE read_Real(key, value, default_value, dohide)
!*************************************************************************  
      use globalvars, only: lglobalprint
      implicit none
     
      character*(*) key
      real*8 value, default_value
      character*250 line
      integer i,le,lk,lt
      logical, optional :: dohide
      logical :: hidden


      if (present(dohide)) then
          hidden = dohide
      else
          hidden = .false.
      endif

      value=default_value
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk)) &
     &       .and.verifyWhitespaces(line(1:lk-1)) &
     &       .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9) )) then
               line=" "//line(le+1:len(line))
               read(line,*,ERR=20,END=20) value       
               if (showinput) write(*,"(T4,A,T25,A,G12.6)") key," = ",value
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      if (lglobalprint .and. .not.hidden) &
     &  print *,"Warning: cannot read value for key = ",key, &
     &          ". Using default value : ",default_value
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      END subroutine

!*************************************************************************
      SUBROUTINE read_RealList(key,N, values, default_values, permit_ellipsis)
!*************************************************************************
!     if permit_ellipsis is true, key may only contain a list with less then 
!     N elements
!     missing elements will be filled with the last accepted value
      use globalvars, only: lglobalprint
      implicit none

      integer N,N_part
      logical permit_ellipsis
      character*(*) key
      real*8 values(*), default_values(*)
      character*250 line
      integer i,le,lk,lt,j

      if (N.le.0) return

      do i=1,N
         values(i)=default_values(i)
      enddo
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key)
            if (((le.ne.0).and.(le.gt.lk)) &
     &       .and.verifyWhitespaces(line(1:lk-1)) &
     &       .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9) )) then
               line=" "//line(le+1:len(line))
               if(permit_ellipsis) then
                  N_part = N+1
 30               if(N_part.ge.2) then
                     N_part = N_part-1
                     read(line,*,ERR=20,END=30) (values(j),j=1,N_part)
                     do j=N_part+1,N
                        values(j)=values(N_part)
                     enddo
                  else
                     goto 10
                  endif
               else
                  read(line,*,ERR=20,END=20) (values(j),j=1,N)
               endif
               if (showinput) write(*,"(T4,A,T25,A,50(G12.5))") key, &
     &              " = ",(values(j),j=1,N)
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      do i=1,N
         values(i)=default_values(i)
      enddo
      if (lglobalprint) then
      print *,"Warning: cannot read value for key = ",key, &
     &      ". Using default value : ",(default_values(j),j=1,N)
      endif
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END subroutine

!*************************************************************************  
      SUBROUTINE read_HistList(key, values, default_values)
!*************************************************************************  
      use globalvars, only: lglobalprint
      implicit none
     
      character*(*) key
      integer N
      real*8 values(*), default_values(*)
      character*250 line
      integer i,j,le,lk,lt

      N = 2
      do i=1,N
         values(i)=default_values(i)
      enddo
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk)) &
     &       .and.verifyWhitespaces(line(1:lk-1)).and.(line(lt:lt).eq. &
     &       " ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9) )) then
               line=" "//line(le+1:len(line))
               read(line,*,ERR=20,END=20) (values(j),j=1,N)       
               if (showinput) write(*,"(T4,A,T25,A,2(G12.5))") key, &
     &              " = ",(values(j),j=1,N)
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      END subroutine

!*************************************************************************  
      SUBROUTINE read_Hist2dList(key, values, default_values)
!*************************************************************************  
      implicit none
     
      character*(*) key
      integer N
      real*8 values(*), default_values(*)
      character*250 line
      integer i,j,le,lk,lt

      N = 4
      do i=1,N
         values(i)=default_values(i)
      enddo
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk)) &
     &       .and.verifyWhitespaces(line(1:lk-1)) &
     &       .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or. &
     &           line(lt:lt).eq.achar(9) )) then
               line=" "//line(le+1:len(line))
               read(line,*,ERR=20,END=20) (values(j),j=1,N)       
               if (showinput) write(*,"(T4,A,T25,A,4(G12.5))")  &
     &              key," = ",(values(j),j=1,N)
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END subroutine
      
!*************************************************************************  
      SUBROUTINE read_Int(key, value, default_value, dohide)
!*************************************************************************  
      use globalvars, only: lglobalprint
      implicit none
     
      character*(*) key
      
      integer value, default_value
      character*250 line
      integer i,le,lk,lt
      logical, optional :: dohide
      logical :: hidden


      if (present(dohide)) then
          hidden = dohide
      else
          hidden = .false.
      endif

      value=default_value
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk)) &
     &           .and.verifyWhitespaces(line(1:lk-1)) &
     &           .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9))) then
               line=" "//line(le+1:len(line))
               read(line,*,ERR=20,END=20) value
               if (showinput) write(*,"(T4,A,T25,A,I8)") key," = ",value
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      if (lglobalprint .and. .not. hidden) then
          print *,"Warning: cannot read value for key = ",key, &
     &      ". Using default value : ",default_value
      endif
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END subroutine

!*************************************************************************  
      SUBROUTINE read_IntList(key, N, values, default_values)
!*************************************************************************  
      use globalvars, only: lglobalprint
      implicit none
     
      character*(*) key
      
      integer N, Nnew, Nsafe
      integer values(*), default_values(*)
      character*250 line
      integer i,j,le,lk,lt

      if (N.le.0) return

      Nsafe = N

      do i=1,N
         values(i)=default_values(i)
      enddo
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk)) &
     &           .and.verifyWhitespaces(line(1:lk-1)) &
     &           .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or. &
     &           line(lt:lt).eq.achar(9))) then
               line=" "//line(le+1:len(line))
               if (key .eq. "LEPTONS") then
                  read(line,*,ERR=20,END=20) (values(j),j=1,1)
                  if (values(1) .lt. 98) then
                     read(line,*,ERR=20,END=20) (values(j),j=1,N)
                  end if                  
               elseif (key .eq. "DECAY_QUARKS") then
                  read(line,*,ERR=20,END=20) (values(j),j=1,1)
                  if (values(1) .lt. 93) then
                     read(line,*,ERR=20,END=20) (values(j),j=1,N)
                  end if                  
               else
                  read(line,*,ERR=20,END=20) (values(j),j=1,N)
               end if
               if (showinput) then 
!     maximum of 50 entries hardwired in format statement
                  if (key .eq. "LEPTONS") then
                     Nnew = N
                     do j = 1, N
                        if (values(j) .ge. 98) Nnew = 1
                     end do
                     N = Nnew
                  end if
                  if (key .eq. "DECAY_QUARKS") then
                     Nnew = N
                     do j = 1, N
                        if (values(j) .ge. 93) Nnew = 1
                     end do
                     N = Nnew
                  end if
                  write(*,"(T4,A,T25,A,50(I5))") key," = ", &
     &                 (values(j),j=1,N)

                  N = Nsafe
               endif
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      if (lglobalprint) then
          print *,"Warning: cannot read values for key = ",key,"."
          print *,"Using default values : ",(default_values(j),j=1,N)
      endif
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END subroutine

!*************************************************************************  
      SUBROUTINE read_Logical(key, value, default_value, dohide)
!*************************************************************************  
      use globalvars, only: lglobalprint
      implicit none

      character*(*) key
      Logical value, default_value
      character*250 line
      integer i,le,lk,lt
      logical, optional :: dohide
      logical :: hidden


      if (present(dohide)) then
          hidden = dohide
      else
          hidden = .false.
      endif

      value=default_value
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk)) &
     &           .and.verifyWhitespaces(line(1:lk-1)) &
     &           .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9))) then
               line=" "//line(le+1:len(line))
               read(line,*,ERR=10,END=20) value
               if (showinput) write(*,"(T4,A,T25,A,L8)") key," = ",value
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      if (lglobalprint .and. .not. hidden) then
          print *,"Warning: cannot read value for key = ",key, &
     &      ". Using default value : ",default_value
      endif
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END subroutine

!*************************************************************************  
      SUBROUTINE read_String(key, value, default_value)
!*************************************************************************  
      use globalvars, only: lglobalprint
      implicit none
     
      character*(*) key
      character*(*) value, default_value
      character*250 line
      integer i,le,lk,lt
     
      value=default_value
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk)) &
     &           .and.verifyWhitespaces(line(1:lk-1)) &
     &           .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9))) then           
               line=" "//line(le+1:len(line))
               read(line,*,ERR=20,END=20) value
               if (showinput) write(*,"(T4,A,T25,A,A)") key," = ",value
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      if (lglobalprint) then
          print *,"Warning: cannot read value for key = ",key, &
     &      ". Using default value : ",default_value
      endif
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END subroutine

!*************************************************************************  
      SUBROUTINE read_StringList(key, N, values, default_values)
!*************************************************************************  
      use globalvars, only: lglobalprint
      implicit none
     
      character*(*) key

      integer N
      character*(*) values(*), default_values(*)
      character*250 line
      integer i,j,le,lk,lt

     
      if (N.le.0) return

      do j = 1,N
         values(j)=default_values(j)
      enddo
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk)) &
     &           .and.verifyWhitespaces(line(1:lk-1)) &
     &           .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9))) then           
               line=" "//line(le+1:len(line))
               read(line,*,ERR=20,END=20) (values(j),j=1,N)
               if (showinput) then
!     maximum of 50 entries hardwired in format statement
!                  write(*,"(T4,A,T25,A,50(A))") key," = ",(values(j),j=1,N)
                   write(*,"(T4,A,T25,A,A)") key," = ",values(1)
                   do j = 2,N
                      write(*,"(T28,A)") values(j)
                   enddo
               endif
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      if (lglobalprint) then
          print *,"Warning: cannot read value for key = ",key, &
     &      ". Using default values : ",(default_values(j),j=1,N)
      endif
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      END subroutine


!*************************************************************************  
      logical FUNCTION verifyWhitespaces(s)
!*************************************************************************  
      implicit none
      integer i
      character*(*) s
      do i=1,len_trim(s)
        if ((s(i:i).ne." ").and.(s(i:i).ne.achar(9))) then
           verifyWhitespaces = .false.
           return
        endif
      enddo
      verifyWhitespaces = .true.
      END function
      
!*************************************************************************
      SUBROUTINE read_cplx(key,value, default_value)
!*************************************************************************
!     if permit_ellipsis is true, key may only contain a list with less then 
!     2 elements 
!     missing elements will be filled with 0D0
      use globalvars, only: lglobalprint

      implicit none

      integer N,N_part
      character*(*) key
      double complex value, default_value
      real*8 values_dum(2), default_values_dum(2)
      character*250 line
      integer i,le,lk,lt,j

      
      N = 2

      value = default_value
      do i=1,2
         values_dum(i)=default_values_dum(i)
      enddo
      if (usedefaults) RETURN


      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key)
            if (((le.ne.0).and.(le.gt.lk)) &
     &       .and.verifyWhitespaces(line(1:lk-1)) &
     &       .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or. &
     &           line(lt:lt).eq.achar(9) )) then
               line=" "//line(le+1:len(line))

               N_part = N+1
 30            if(N_part.ge.2) then
                  N_part = N_part-1
                  read(line,*,ERR=20,END=30) (values_dum(j),j=1,N_part)
                  do j=N_part+1,N
                     values_dum(2)=0D0
                  enddo
               else
                  goto 10
               endif

               value = CMPLX(values_dum(1),values_dum(2))
               if (showinput) write(*,"(T4,A,T25,A,50(G12.5))") key, &
     &              " = ",(values_dum(j),j=1,N)
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      do i=1,N
         values_dum(i)=default_values_dum(i)
      enddo
      value = default_value
      if (lglobalprint) then
          print *,"Warning: cannot read value for key = ",key, &
     &      ". Using default value : ",default_value
      endif
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      END subroutine

end module
