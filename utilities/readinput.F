c*************************************************************************  
c     This file contains subroutines for reading input files. 
c     The input files are in the format: key = value 
c*************************************************************************  
c   LIST OF ALL FUNCTIONS AND SUBROUTINES IN THIS FILE:
c
c     SUBROUTINE loadfile(filename)
c     SUBROUTINE closefile
c     SUBROUTINE read_Real(key, value, default_value)
c     SUBROUTINE read_RealList(key, N, values, default_values)
c     SUBROUTINE read_Int(key, value, default_value)
c     SUBROUTINE read_IntList(key, N, values, default_values)
c     SUBROUTINE read_Logical(key, value, default_value)
c     SUBROUTINE read_String(key, value, default_value)
c     SUBROUTINE read_StringList(key, N, values, default_values)
c     FUNCTION   GetInputPath()
c
c*************************************************************************  



c*************************************************************************  
      subroutine loadfile(filename,show_input)
c*************************************************************************  
c     Opens a textfile and reads the data in an internal buffer.
c     End-of-line comments, starting with ! or //, are removed.
c     Do only process one single file at a time.
c     DON"T FORGET TO CALL CLOSEFILE AFTERWARDS.
c*************************************************************************  
      implicit none

      character*(*) filename
     
#include "readinput.inc"
      
      integer i,le,k,j
      character*250 line
      logical show_input

      character*250 GetInputPath
      character*250 path
      external GetInputPath
      

      iunit=10
      usedefaults=.false.
      showinput=show_input
c add input file path to filename
      path=GetInputPath()
      path = path(1:len_trim(path))//"/"//filename
c Open file
      open (unit=iunit,
     &      file=path(1:len_trim(path)),
     &      status="old",err=77)

      fname = filename

c read text lines
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

c remove comments
      do i=1,n_lines
         line=lines(i)
         if (INDEX(line,"!").ne.0) then
            line=line(1:INDEX(line,"!")-1)
         endif
         if (INDEX(line,"//").ne.0) then
            line=line(1:INDEX(line,"//")-1)
         endif
c upper case everything in front of a "="
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
cc     if (showinput) write(*,"(A,A,A)"), " reading input from file ",filename," : "
      
      return
 77   CONTINUE
      print *," Open error on file ",filename
      print *," Using default values."
      usedefaults=.true.

      return
      end

c*************************************************************************  
      subroutine loadfileKK(filename,show_input)
c*************************************************************************  
c     Opens a textfile and reads the data in an internal buffer.
c     End-of-line comments, starting with ! or //, are removed.
c     Do only process one single file at a time.
c     DON"T FORGET TO CALL CLOSEFILE AFTERWARDS.
c IMPORTANT NOTE:  This routine is only used when reading in the KK
c                  data that's just been calculated and stored to file.
c                  It's not a particularly beautiful solution, but it 
c                  works for now.  See wiki.
c*************************************************************************  
      implicit none

      character*(*) filename
     
#include "readinput.inc"
      
      integer i,le,k,j
      character*250 line
      logical show_input

      character*250 GetInputPath
      character*250 path
      external GetInputPath
      

      iunit=10
      usedefaults=.false.
      showinput=show_input
c add input file path to filename
      path = filename
*      path=GetInputPath()
*      path = path(1:len_trim(path))//"/"//filename
c Open file
      open (unit=iunit,
     &      file=path(1:len_trim(path)),
     &      status="old",err=77)

      fname = filename

c read text lines
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

c remove comments
      do i=1,n_lines
         line=lines(i)
         if (INDEX(line,"!").ne.0) then
            line=line(1:INDEX(line,"!")-1)
         endif
         if (INDEX(line,"//").ne.0) then
            line=line(1:INDEX(line,"//")-1)
         endif
c upper case everything in front of a "="
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
cc     if (showinput) write(*,"(A,A,A)"), " reading input from file ",filename," : "
      
      return
 77   CONTINUE
      print *," Open error on file ",filename
      print *," Using default values."
      usedefaults=.true.

      return
      end

c*************************************************************************  
      subroutine closefile
c*************************************************************************  
c     Closes the current file.
c*************************************************************************  
      implicit none

#include "readinput.inc"
      
      close(iunit)
      return
      end

c*************************************************************************  
      SUBROUTINE read_Real(key, value, default_value)
c*************************************************************************  
      implicit none
     
#include "readinput.inc"

      character*(*) key
      real*8 value, default_value
      character*250 line
      integer i,le,lk,lt
      logical verifyWhitespaces
      external verifyWhitespaces
     
      value=default_value
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk))
     &       .and.verifyWhitespaces(line(1:lk-1))
     &       .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9) )) then
               line=" "//line(le+1:len(line))
               read(line,*,ERR=20,END=20) value       
               if (showinput) write(*,"(T4,A,T25,A,G12.6)") key," = ",value
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      print *,"Warning: cannot read value for key = ",key,
     &      ". Using default value : ",default_value
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END

c*************************************************************************
      SUBROUTINE read_RealList(key,N, values, default_values, permit_ellipsis)
c*************************************************************************
c     if permit_ellipsis is true, key may only contain a list with less then 
c     N elements
c     missing elements will be filled with the last accepted value
      implicit none

#include "readinput.inc"

      integer N,N_part
      logical permit_ellipsis
      character*(*) key
      real*8 values(*), default_values(*)
      character*250 line
      integer i,le,lk,lt,j
      logical verifyWhitespaces
      external verifyWhitespaces

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
            if (((le.ne.0).and.(le.gt.lk))
     &       .and.verifyWhitespaces(line(1:lk-1))
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
               if (showinput) write(*,"(T4,A,T25,A,50(G12.5))") key,
     -              " = ",(values(j),j=1,N)
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      do i=1,N
         values(i)=default_values(i)
      enddo
      print *,"Warning: cannot read value for key = ",key,
     &      ". Using default value : ",(default_values(j),j=1,N)
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END

c*************************************************************************  
      SUBROUTINE read_HistList(key, values, default_values)
c*************************************************************************  
      implicit none
     
#include "readinput.inc"
      
      character*(*) key
      integer N
      real*8 values(*), default_values(*)
      character*250 line
      integer i,j,le,lk,lt
      logical verifyWhitespaces
      external verifyWhitespaces

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
            if (((le.ne.0).and.(le.gt.lk))
     &       .and.verifyWhitespaces(line(1:lk-1)).and.(line(lt:lt).eq.
     &       " ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9) )) then
               line=" "//line(le+1:len(line))
               read(line,*,ERR=20,END=20) (values(j),j=1,N)       
               if (showinput) write(*,"(T4,A,T25,A,2(G12.5))") key,
     -              " = ",(values(j),j=1,N)
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
!      print *,"Warning: cannot read value for key = ",key,
!     &      ". Using default value : ",(default_values(j),j=1,N)
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END

c*************************************************************************  
      SUBROUTINE read_Hist2dList(key, values, default_values)
c*************************************************************************  
      implicit none
     
#include "readinput.inc"
      
      character*(*) key
      integer N
      real*8 values(*), default_values(*)
      character*250 line
      integer i,j,le,lk,lt
      logical verifyWhitespaces
      external verifyWhitespaces

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
            if (((le.ne.0).and.(le.gt.lk))
     &       .and.verifyWhitespaces(line(1:lk-1))
     &       .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.
     &           line(lt:lt).eq.achar(9) )) then
               line=" "//line(le+1:len(line))
               read(line,*,ERR=20,END=20) (values(j),j=1,N)       
               if (showinput) write(*,"(T4,A,T25,A,4(G12.5))") 
     &              key," = ",(values(j),j=1,N)
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
!      print *,"Warning: cannot read value for key = ",key,
!     &      ". Using default value : ",(default_values(j),j=1,N)
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END
      
c*************************************************************************  
      SUBROUTINE read_Int(key, value, default_value)
c*************************************************************************  
      implicit none
     
#include "readinput.inc"

      character*(*) key
      
      integer value, default_value
      character*250 line
      integer i,le,lk,lt
      logical verifyWhitespaces
      external verifyWhitespaces
      logical hidden
      
      hidden = .false.
      goto 1
      
      entry read_Int_hidden(key, value, default_value)
      hidden=.true.
      
 1    continue

      value=default_value
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk))
     &           .and.verifyWhitespaces(line(1:lk-1))
     &           .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9))) then
               line=" "//line(le+1:len(line))
               read(line,*,ERR=20,END=20) value
               if (showinput) write(*,"(T4,A,T25,A,I8)") key," = ",value
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      if(.not.hidden) print *,"Warning: cannot read value for key = ",key,
     &      ". Using default value : ",default_value
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END

c*************************************************************************  
      SUBROUTINE read_IntList(key, N, values, default_values)
c*************************************************************************  
      implicit none
     
#include "readinput.inc"

      character*(*) key
      
      integer N, Nnew, Nsafe
      integer values(*), default_values(*)
      character*250 line
      integer i,j,le,lk,lt
      logical verifyWhitespaces
      external verifyWhitespaces

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
            if (((le.ne.0).and.(le.gt.lk))
     &           .and.verifyWhitespaces(line(1:lk-1))
     &           .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.
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
c     maximum of 50 entries hardwired in format statement
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
                  write(*,"(T4,A,T25,A,50(I5))") key," = ",
     &                 (values(j),j=1,N)

                  N = Nsafe
               endif
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      print *,"Warning: cannot read values for key = ",key,"."
      print *,"Using default values : ",(default_values(j),j=1,N)
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END

c*************************************************************************  
      SUBROUTINE read_Logical(key, value, default_value)
c*************************************************************************  
      implicit none

#include "readinput.inc"
      
      character*(*) key
      Logical value, default_value
      character*250 line
      integer i,le,lk,lt
      logical verifyWhitespaces
      external verifyWhitespaces
      logical hidden
      
      hidden=.false.
      goto 1
      
      entry read_Logical_hidden(key, value, default_value)
      hidden=.true.
      
 1    continue

      value=default_value
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk))
     &           .and.verifyWhitespaces(line(1:lk-1))
     &           .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9))) then
               line=" "//line(le+1:len(line))
               read(line,*,ERR=10,END=20) value
               if (showinput) write(*,"(T4,A,T25,A,L8)") key," = ",value
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      if(.not.hidden) print *,"Warning: cannot read value for key = ",key,
     &      ". Using default value : ",default_value
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END

c*************************************************************************  
      SUBROUTINE read_String(key, value, default_value)
c*************************************************************************  
      implicit none
     
#include "readinput.inc"

      character*(*) key
      character*(*) value, default_value
      character*250 line
      integer i,le,lk,lt
      logical verifyWhitespaces
      external verifyWhitespaces
     
      value=default_value
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk))
     &           .and.verifyWhitespaces(line(1:lk-1))
     &           .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9))) then           
               line=" "//line(le+1:len(line))
               read(line,*,ERR=20,END=20) value
               if (showinput) write(*,"(T4,A,T25,A,A)") key," = ",value
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      print *,"Warning: cannot read value for key = ",key,
     &      ". Using default value : ",default_value
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END

c*************************************************************************  
      SUBROUTINE read_StringList(key, N, values, default_values)
c*************************************************************************  
      implicit none
     
#include "readinput.inc"

      character*(*) key

      integer N
      character*(*) values(*), default_values(*)
      character*250 line
      integer i,j,le,lk,lt
      logical verifyWhitespaces
      external verifyWhitespaces

     
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
            if (((le.ne.0).and.(le.gt.lk))
     &           .and.verifyWhitespaces(line(1:lk-1))
     &           .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9))) then           
               line=" "//line(le+1:len(line))
               read(line,*,ERR=20,END=20) (values(j),j=1,N)
               if (showinput) then
c     maximum of 50 entries hardwired in format statement
c                  write(*,"(T4,A,T25,A,50(A))") key," = ",(values(j),j=1,N)
                   write(*,"(T4,A,T25,A,A)") key," = ",values(1)
                   do j = 2,N
                      write(*,"(T28,A)"),values(j)
                   enddo
               endif
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      print *,"Warning: cannot read value for key = ",key,
     &      ". Using default values : ",(default_values(j),j=1,N)
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END


c*************************************************************************  
      logical FUNCTION verifyWhitespaces(s)
c*************************************************************************  
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
      RETURN
      END
     
c*************************************************************************
      SUBROUTINE read_cplx(key,value, default_value)
c*************************************************************************
c     if permit_ellipsis is true, key may only contain a list with less then 
c     2 elements 
c     missing elements will be filled with 0D0

      implicit none

#include "readinput.inc"

      integer N,N_part
      character*(*) key
      double complex value, default_value
      real*8 values_dum(2), default_values_dum(2)
      character*250 line
      integer i,le,lk,lt,j
      logical verifyWhitespaces
      external verifyWhitespaces

      
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
            if (((le.ne.0).and.(le.gt.lk))
     &       .and.verifyWhitespaces(line(1:lk-1))
     &       .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.
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
               if (showinput) write(*,"(T4,A,T25,A,50(G12.5))") key,
     -              " = ",(values_dum(j),j=1,N)
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      do i=1,N
         values_dum(i)=default_values_dum(i)
      enddo
      value = default_value
      print *,"Warning: cannot read value for key = ",key,
     &      ". Using default value : ",default_value
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END

c*************************************************************************  
      subroutine ParseOptions()
c*************************************************************************  
c go through list of command-line options
c handle everything and warn about all unknown options
c*************************************************************************  
      implicit none
      character*250 s, sorig
      integer narg,i,j,k,lk,le,lt
      
#include "global.inc"
#include "readinput.inc"

      inputpath=""
      pdfpath=""

c go through list
      narg = IARGC()
      if (narg.gt.0)then
         do i=1,narg
            call GETARG(i , s)
            sorig=s
c upper case everything in front of a "="
            le = index(s,"=")
c if there is no "=" upper case everything
            if (le.eq.0) le = len_trim(s)+1
            do j=1,le-1
               k=ichar(s(j:j))
               if(k.ge.97.and.k.le.122) then
                  k=ichar(s(j:j))-32   
                  s(j:j)=char(k)
               endif
            enddo
c remove leading dashes
            do while (s(1:1).eq.'-')
              s = s(2:len(s))
            enddo
            le = index(s,"=")
            if (le.eq.0) le=len_trim(s)+1
            if (s(1:le-1).eq."INPUT") then
                inputpath=s(le+1:len(s))
            elseif (s(1:le-1).eq."PDFSETS") then
                pdfpath=s(le+1:len(s))
            elseif (s(1:le-1).eq."VERSION") then
                call PrintVBFNLOVersion(-1)
                print *, "VBFNLO ", 
     &            vbfnloversionstring(1:len_trim(vbfnloversionstring))
                stop
            else
                print *, "Warning: Unknown option ", 
     &                   sorig(1:len_trim(sorig))
            endif
         enddo
      endif

      END

c*************************************************************************  
      function VBFNLOCapabilities(option)
c*************************************************************************  
c test all options VBFNLO can handle
c*************************************************************************  
      implicit none

#include "global.inc"

      integer VBFNLOCapabilities
      character(*) option


c general
      IF (option(1:len_trim(option)).eq."VERSION") then
      VBFNLOCapabilities = vbfnloversionnumber

c processes
#ifdef WITH_VBF
      ELSEIF (option(1:len_trim(option)).eq."VBF") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_HJJJ
      ELSEIF (option(1:len_trim(option)).eq."HJJJ") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_DIBOSON
      ELSEIF (option(1:len_trim(option)).eq."DIBOSON") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_DIBOSONJET
      ELSEIF (option(1:len_trim(option)).eq."DIBOSONJET") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_TRIBOSON
      ELSEIF (option(1:len_trim(option)).eq."TRIBOSON") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_TRIBOSONJET
      ELSEIF (option(1:len_trim(option)).eq."TRIBOSONJET") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_GGF
      ELSEIF (option(1:len_trim(option)).eq."GGF") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_DY
      ELSEIF (option(1:len_trim(option)).eq."DY") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_QCDV
      ELSEIF (option(1:len_trim(option)).eq."QCDV") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_QCDVV
      ELSEIF (option(1:len_trim(option)).eq."QCDVV") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_TOP
      ELSEIF (option(1:len_trim(option)).eq."TOP") then
      VBFNLOCapabilities = 1
#endif

c options
#ifdef WITH_KK
      ELSEIF (option(1:len_trim(option)).eq."KK") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_NLO
      ELSEIF (option(1:len_trim(option)).eq."NLO") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_SPIN2
      ELSEIF (option(1:len_trim(option)).eq."SPIN2") then
      VBFNLOCapabilities = 1
#endif

c programs
#ifdef WITH_FH
      ELSEIF (option(1:len_trim(option)).eq."FH") then
#ifdef WITH_FH26
      VBFNLOCapabilities = 26
#elif  defined WITH_FH27
      VBFNLOCapabilities = 27
#elif defined WITH_FH28
      VBFNLOCapabilities = 28
#else
      VBFNLOCapabilities = 1
#endif
#endif
#ifdef WITH_HEPMC
      ELSEIF (option(1:len_trim(option)).eq."HEPMC") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_LHA
      ELSEIF (option(1:len_trim(option)).eq."LHA") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_LT
      ELSEIF (option(1:len_trim(option)).eq."LT") then
      VBFNLOCapabilities = 1
#endif
#ifdef WITH_ROOT
      ELSEIF (option(1:len_trim(option)).eq."ROOT") then
      VBFNLOCapabilities = 1
#endif

      ELSE
        VBFNLOCapabilities = 0

      ENDIF

      return

      end
