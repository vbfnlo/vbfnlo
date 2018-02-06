#include "VBFNLOConfig.h"

!!! THIS IS A FORTRAN 90 FILE: CONTINUATION LINES AND COMMENTS HAVE A DIFFERENT SYNTAX !!!


!*************************************************************************
      SUBROUTINE InitPDFs(nlo)
!*************************************************************************
!     Initialize parton distribution functions.
!     There should be two pdfsets defined. One for the LO and one for the
!     NLO calculation. 
!*************************************************************************
          use globalvars, only: lglobalprint
          use readinput
      implicit none

      character*250 pdfname1, pdfname2, pdfname
      integer pdflib, value1, value2, valueLHA, value
      integer pdfnumber1, pdfnumber2, pdfnumber
      COMMON/PDFparameters/pdflib, valueLHA
      COMMON/cteqchoice/ value
      character*20 lhaparm(20)
      double precision lhavalue(20)
      character*250 pdfwithpath, shortpdfname, pdfindex
      logical isSingleMemberPDF

      integer nlo, nmem

#include "global.inc"

!      call loadfile("pdfs.dat")
      call loadfile("vbfnlo.dat",.false.)
      call read_int("PDF_SWITCH",pdflib,0)
      SELECT CASE(pdflib)
      CASE(0)
         pdfname1=''
         pdfname2=''
         pdfnumber1=4
         value1 = 6
         pdfnumber2=100
         value2 = 10
#ifdef WITH_LHA
      CASE(1)
         value1 = -1
         value2 = -1
         call read_string("LO_PDFNAME",pdfname1,"LO_PDFSET")
         if (pdfname1 .eq. "LO_PDFSET") then
           pdfname1 = "cteq6ll.LHpdf"
           call read_int("LO_PDFSET",value1,-1)
         endif
         call read_string("NLO_PDFNAME",pdfname2,"NLO_PDFSET")
         if (pdfname2 .eq. "NLO_PDFSET") then
           pdfname2 = "CT10.LHgrid"
           call read_int("NLO_PDFSET",value2,-1)
         endif
         call read_int("LO_PDFMEMBER",pdfnumber1,0)
         call read_int("NLO_PDFMEMBER",pdfnumber2,0)
#else
      CASE(1)
         print *,"Error: LHAPDF library not included."
         print *,"       Add the flag --with-LHAPDF=/path/to/LHAPDF/ to configure."
         STOP
#endif
! Case(2) is hard-wired mrst2004qed
      CASE(2)
         pdfname1=''
         pdfname2=''
         pdfnumber1=4
         pdfnumber2=1
! Case(3) is hard-wired mstw2008
      CASE(3)
         pdfname1=''
         pdfname2=''
         pdfnumber1=4
         pdfnumber2=1
      CASE DEFAULT
         print *,"Error: invalid choice of pdf : PDF_SWITCH = ",pdflib
         STOP
      END SELECT
     
      call closefile

      if (nlo.eq.0) then
         pdfnumber = pdfnumber1
         pdfname = pdfname1
         value = value1
      else
         pdfnumber = pdfnumber2
         pdfname = pdfname2
         value = value2
      endif

      SELECT CASE(pdflib)
      CASE(0)
         if (value .eq. 6) then
             if (lglobalprint) write(*,*)'Initializing PDF set cteq6l1'
           Call SetCtq6(pdfnumber)
         else if (value .eq. 10) then
             if (lglobalprint) write(*,*)'Initializing PDF set CT10'
           Call SetCT10(pdfnumber)
         else 
           print *,"Unknown type of PDFs (value=",value,") in pdfsets.F"
           stop
         endif
#ifdef WITH_LHA
      CASE(1)

! try to figure out where pdf is located
         if (value .ne. -1) then
! use old method of by number

#ifdef WITH_LHAPDF6
           ! routines setlhaparm, SetPDFPath and GetDESC are not available or empty in LHAPDF 6.
           ! since the "SILENT" option is not available, LHAPDF will write something like
           ! "POWHEG using LHAPDFv6", but it will still work correctly...
           lhaparm(1)='DEFAULT'
           lhavalue(1)=value
           valueLHA = value
           print*, ""
           ! PDFSET throws an exception if the lhavalue does not exist instead of exiting normally.
           ! Since fortran 95 does not have a standardized error handling it the fortran library will handle it,
           ! which will look rather ugly. Checking the validity of the lhavalue beforehand on the other hand
           ! is non-trivial (would imply to parse through all PDF folders and check the .info file).
           ! Given that the default way to set the PDFs is "by name" there is no effort done in checking lhavalue.
           call PDFSET(lhaparm,lhavalue)
           print*, ""
#else
           call setlhaparm('SILENT')
           CALL SetPDFPath(LHAPDFPATH)
           lhaparm(1)='DEFAULT'
           lhavalue(1)=value
           valueLHA = value
           call PDFSET(lhaparm,lhavalue)
           write(*,*) ">>>>>> PDF description: <<<<<<"
           call GetDESC()
           write(*,*) ">>>>>>                  <<<<<<"
#endif

         else  ! use "new" method of identifying the pdf by name

           ! check is the given pdf can be found. provide additionally the path to the PDFsets.index file
           ! and a short pdf name, which is without the path in case of LHAPDF <6 and without file 
           call getLHAPDFname(pdfname, shortpdfname, pdfwithpath, pdfindex)

           print*, ""
#ifdef WITH_LHAPDF6
           ! In LHAPDF6 you don't provide the file name of the pdf set, but instead the pdf set name. 
           ! For compatibility, LHAPDF6 strips the file ending, so you can still provide a
           ! pdf name like CT10.LHgrid .
           ! The PDF should be verified beforehand, since LHAPDF throws an error if it's non-existing,
           ! and FORTRAN 95 does not have an exception handling...
           ! The existence of the PDF is checked by getLHAPDFname, which also strips the file ending,
           ! providing shortpdfname .
           if (trim(pdfwithpath).ne."") then
              call InitPDFsetbyName(shortpdfname)
           else
              print*, "Error: The PDF named ", trim(shortpdfname), " could not be found."
              print*, "The PDFs are expected to reside"
              if (LHAPDFPATH.ne."") print*, "in the folder "//LHAPDFPATH//" or"
              print*, "in a folder given by the environmental variable $LHAPDF_DATA_PATH ."
              print*, ""
              stop
           endif
#else
           ! The path to the PDF is checked by getLHAPDFname.
           if (trim(pdfwithpath).ne."") then
              ! found correct PDF file name myself
              call InitPDFset(pdfwithpath)
           else
              ! maybe LHAPDF itself is more lucky ...
              print*, "Warning: VBFNLO could not find the PDF named ", trim(shortpdfname), "."
              print*, "The PDFs are expected to reside"
              if (LHAPDFPATH.ne."") print*, "in the folder "//LHAPDFPATH//" or"
              print*, "in a folder given by the environmental variable $LHAPATH ."
              print*, "LHAPDF itself will now try to find the PDF..."
              print*, ""
             call InitPDFsetbyName(pdfname)
           endif
#endif

           ! get number of members of current pdf set
           call numberPDF(nmem)     ! returns n-1 for n>1 => same result for n=1 and n=2 => leads to problems for pdfs with n=1
           call getPDFnoLHAPDF(shortpdfname, pdfwithpath, pdfindex, valueLHA, isSingleMemberPDF)      ! Finding lhaglue number of pdf for the event output and check for n=1
           if (nmem.eq.1) then      ! could be 1 or 2 members
              if (isSingleMemberPDF) then
                 nmem = 0
                 if (pdfnumber.eq.1) pdfnumber = 0
              endif
           endif
           if (pdfnumber.lt.0 .or. pdfnumber.gt.nmem) then
              print*, ""
              print*, "The PDF member number given in vbfnlo.dat"
              print*, "(LO_PDFMEMBER or NLO_PDFMEMBER =",pdfnumber,")"
              print*, "is not valid, exiting!"
              stop
           endif
           call InitPDF(pdfnumber)        ! initialize requested member of current pdf set
           print*, ""

           ! lhaglue number for event output
           if (valueLHA .ne. 0) then
              valueLHA = valueLHA + pdfnumber
           else                                 ! could not determine base number of PDF set
              valueLHA = -1
           end if
         endif

! sadly lhapdf does not return nf for LHgrid pdfs (most of them)
! this should be fixed in the future
! currently (v.6.0.5) the getNf function of LHAPDF 6 is completely empty
!        call getNf(lhapdfnf)
!        if ( lhapdfnf .ne. nfl) then
!           write(*,*) 'pdf uses ', lhapdfnf, ' flavors, but vbfnlo uses ', nfl
!           write(*,*) 'make sure to use a consistent flavor scheme!'
!        endif

#endif
      CASE(3)
         call init_mstw(nlo)
      END SELECT


      RETURN
      END

!*************************************************************************  
      SUBROUTINE pdfproton(x,q,pdf)
!*************************************************************************
!     This should be the interface for all PDFs, i.e. this subroutine is
!     called inside the matrix elements. 
!     NOTE: GAMpdf is photon entry.  It is set to zero unless we're using
!           MRST2004qed
!*************************************************************************
      implicit none
      real*8 x,q,pdf(-6:6)

      integer pdflib, mu, valueLHA
      COMMON/PDFparameters/pdflib, valueLHA

#include "Apdf.inc"      
    
! Initialize all PDF to zero
      do mu=-6,6
        pdf(mu) = 0d0
      enddo
      GAMpdf = 0D0

      SELECT CASE(pdflib)

      CASE(0)
         call pftopdg_cteq(x,q,pdf)

#ifdef WITH_LHA
      CASE(1)
         call evolvePDF(x,q,pdf)
#endif

! Calling the hard-wired mrst2004qed
      CASE(2)
         call mrst2004qed(1,x,q,pdf,GAMpdf)

! Calling hard-wired mstw2008
      CASE(3)
         call mstw2008(x,q,pdf)

      END SELECT
    

      RETURN
      END

!*************************************************************************  
      SUBROUTINE printnfl(lchg)
!*************************************************************************
      implicit none
      logical lchg
      integer rnfl

      character quark(0:6)
      data quark /'g','d','u','s','c','b','t'/

#include "global.inc"
#include "process.inc"


!     print *," Including pdf contributions up to ",quark(rnfl)," quarks"
      if (nfl.ge.6) then
         print*," Warning: treating the top-quark as massless"
         print*," Error:   This is not a supported option."
         stop
      endif


! Determine how far up we actually go
      SELECT CASE(procID)
      CASE(Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar, &
          Hjj_WW, Hjj_ZZ_ll, Hjj_ZZ_lnu, &
          Hjjj, Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_bbar, &
          Hjjj_WW, Hjjj_ZZ_ll, Hjjj_ZZ_lnu, &
          HAjj, HAjj_AA, HAjj_mu, HAjj_tau, HAjj_bbar, &
          HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu, &
          Zjj_l, Zjj_nu, Ajj, &
          WPWMjj, ZZjj_ll, ZZjj_lnu, ZAjj, ZAjj_n, AAjj, HHjj, &
          HHjj_bbtau, HHjj_bbAA, &
          WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, &
          WMhadZjj, WMZhadjj, ZZhadjj, WPhadWPjj, WMhadWMjj, &
          Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad, Sp2jj_WW, &
          Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu)
! For VBF processes, bottom quarks are included for neutral currents iff
! nfl = 5 and vbfNFLb = .true. and ewcor_switch .eq. false
        if (vbfNFLb .and. (.not. ewcor_switch)) then
           nflVBF = nfl
        else if (nfl .ge. 4) then
           nflVBF = 4
        else if (nfl .ge. 2) then
           nflVBF = 2
        end if
        if (nflVBF .ge. 5) then
           print*," External b-quark contributions are included in neutral currents"
        else
           print*," External b-quark contributions are excluded"
        end if
      CASE default
! If there is a W attached to the quark line (signalled by lchg==.true.),
! the corresponding up-type quark needs to be covered as well
         if (lchg) then
            rnfl = (nfl/2)*2
         else
            rnfl = nfl
         endif

         if (rnfl.ge.5) then
            print*," External b-quark contributions are included"
         else
            print*," External b-quark contributions are excluded"
         endif
      END SELECT

      return
      end


#ifdef WITH_LHA


!*************************************************************************  
      subroutine getLHAPDFname(pdfname_in,pdfname,fullpath,pdfindex)
!*************************************************************************  
!     Tries to find the full path to the PDF file, the PDF name (without
!     path) and the path to the pdfsets.index file
!*************************************************************************  
      implicit none
      
      character*250, intent(in)  :: pdfname_in
      character*250, intent(out) :: pdfname,fullpath,pdfindex

      logical pdfindexFound, pdfFound
      
      integer lp

      ! defaults
      pdfname = ""
      fullpath = ""
      pdfindex = ""


! get the PDF name without path (and for LHAPDF6 without file ending)
      pdfname = trim(pdfname_in)
#ifdef WITH_LHAPDF6
      ! remove file extensions from pdfname if any. full path name not supported here, so no need to remove it
      lp = INDEX(pdfname,".LHgrid")
      if (lp.gt.0) pdfname = pdfname(1:lp-1)
      lp = INDEX(pdfname,".LHpdf")
      if (lp.gt.0) pdfname = pdfname(1:lp-1)
      if (trim(pdfname).eq."cteq6ll") pdfname = "cteq6l1"      ! for this pdf set: filename != pdf name in LHAPDF<6
#else
      ! remove path from pdfname if any
      lp = INDEX(pdfname,"/",.true.)
      if (lp.gt.0) pdfname = pdfname(lp+1:len_trim(pdfname))
#endif


! get the PDF name with full path (for LHAPDF6: path to .info file of the pdf set)
      pdfFound = .false.
#ifdef WITH_LHAPDF6
      call getenv("LHAPDF_DATA_PATH",fullpath)      ! "prefered" variable of LHAPDF 6
      if (trim(fullpath) .ne. "") then
         fullpath = trim(fullpath) //"/"// trim(pdfname) &
                   //"/"// trim(pdfname) // ".info"
         INQUIRE(FILE=fullpath,EXIST=pdfFound)
      endif
      if (.not.pdfFound) then
!         call getenv("LHAPATH",fullpath)            ! "backup" variable, common for LHAPDF <6, does not seem to work any longer
!         if (trim(fullpath) .ne. "") then
!            fullpath = trim(fullpath) //"/"// trim(pdfname) &
!                      //"/"// trim(pdfname) // ".info"
!            INQUIRE(FILE=fullpath,EXIST=pdfFound)
!         endif
         if (.not.pdfFound) then                    ! if not found yet use PDFPATH from configure stage
            fullpath = LHAPDFPATH
            fullpath = trim(fullpath) //"/"// trim(pdfname) &
                      //"/"// trim(pdfname) // ".info"
            INQUIRE(FILE=fullpath,EXIST=pdfFound)
        endif
      endif
#else
      ! try pdfname as it is
      fullpath =  trim(pdfname_in)
      INQUIRE(FILE=fullpath,EXIST=pdfFound)
      if (.not.pdfFound) then
         ! prepend $LHAPATH
         call getenv("LHAPATH",fullpath)
         ! but only if it is set
         if (trim(fullpath) .ne. "") then
            fullpath = trim(fullpath)//'/'//trim(pdfname_in)
            INQUIRE(FILE=fullpath,EXIST=pdfFound)
         endif
         if (.not.pdfFound) then
            ! prepend PDFPATH from configure stage
            fullpath = LHAPDFPATH // "/" // trim(pdfname_in)
            INQUIRE(FILE=fullpath,EXIST=pdfFound)
         endif
      endif
#endif
      if (.not.pdfFound) fullpath = ""


! get the full path to PDFsets.index / pdfsets.index
      pdfindexFound = .false.
#ifdef WITH_LHAPDF6
      call getenv("LHAPDF_DATA_PATH",pdfindex)      ! "prefered" variable of LHAPDF 6
      if (trim(pdfindex) .ne. "") then
         pdfindex = trim(pdfindex)//"/pdfsets.index"
         INQUIRE(FILE=pdfindex,EXIST=pdfindexFound)
      endif
      if (.not.pdfindexFound) then
!         call getenv("LHAPATH",pdfindex)            ! "backup" variable, common for LHAPDF <6, does not seem to work with LHAPDF6
!         if (trim(pdfindex) .ne. "") then
!            pdfindex = trim(pdfindex)//"/pdfsets.index"
!            INQUIRE(FILE=pdfindex,EXIST=pdfindexFound)
!         endif
         if (.not.pdfindexFound) then               ! if not found yet use PDFPATH from configure stage
            pdfindex = LHAPDFPATH
            pdfindex = trim(pdfindex)//"/pdfsets.index"
            INQUIRE(FILE=pdfindex,EXIST=pdfindexFound)
        endif
      endif
#else
      call getenv("LHAPATH",pdfindex)               ! try environment variable first
      if (trim(pdfindex) .ne. "") then
         lp = INDEX(pdfindex,"PDFsets")             ! the PDFsets.index is not in the PDFsets folder, but in the parent folder
         pdfindex = pdfindex(1:lp-1)
         pdfindex = trim(pdfindex)//"/PDFsets.index"
         INQUIRE(FILE=pdfindex,EXIST=pdfindexFound)
      endif
      if (.not.pdfindexFound) then                  ! if not found yet use PDFPATH from configure stage
         pdfindex = LHAPDFPATH
         lp = INDEX(pdfindex,"PDFsets")             ! the PDFsets.index is not in the PDFsets folder, but in the parent folder
         pdfindex = pdfindex(1:lp-1)
         pdfindex = trim(pdfindex)//"/PDFsets.index"
         INQUIRE(FILE=pdfindex,EXIST=pdfindexFound)
      endif
#endif
      if (.not.pdfindexFound) pdfindex = ""


      return
      
      end


!*******************************************************************************************  
      subroutine getPDFnoLHAPDF(shortpdfname, pdfwithpath, pdfindex, value, isSingleMemberPDF)
!*******************************************************************************************
!     Reads the lhapdf file PDFsets.index to find the PDF number for the
!     les houches event output and furthermore checks if the given pdf has only one member
!*******************************************************************************************  
      implicit none

      character*250, intent(in)  :: shortpdfname, pdfwithpath, pdfindex
      integer, intent(out) :: value
      logical, intent(out) :: isSingleMemberPDF
      
      integer iunit, io_error, nummembers, nextvalue
      logical pdfsetsFound

      character*250 line, templine

      integer i, lk, li, lf

      iunit=10

      ! default: ID not found
      value = 0
      nummembers = 0
      isSingleMemberPDF = .false.


#ifdef WITH_LHAPDF6
      ! check lhaglue number from "SetIndex" in .info file, PDFsets.index, and check "NumMembers" for isSingleMemberPDF
      if (trim(pdfwithpath).eq."") return
      ! Open file
      open (unit=iunit, file=trim(pdfwithpath), status="old", iostat=io_error)
      if (io_error.eq.0) then
         do
            read(iunit,"(a100)",iostat=io_error) line
            if (io_error.ne.0) exit
            ! check lhaglue number
            lk=INDEX(line,"SetIndex:")
            if (lk .ne. 0) then
               ! 2nd and last entry of line is PDF number
               li=INDEX(line,":")
               line=line(li+1:len_trim(line))
               read(line,*) value
            end if
            ! check number of members
            lk=INDEX(line,"NumMembers:")
            if (lk .ne. 0) then
               ! 2nd and last entry of line is number of members
               li=INDEX(line,":")
               line=line(li+1:len_trim(line))
               read(line,*) nummembers
               if (nummembers.eq.1) then
                  isSingleMemberPDF = .true.
               endif
            end if
            if (nummembers.ne.0 .and. value.ne.0) exit
         enddo
         close(iunit)
      endif
#else
      ! check lhaglue number from PDFsets.index, afterwards check if lhaglue+1 still corresponds to the same pdf set
      if (trim(pdfindex).eq."") return
      ! Open file
      open (unit=iunit, file=trim(pdfindex), status="old", iostat=io_error)
      if (io_error.eq.0) then
         ! get lhaglue number
         do
            read(iunit,"(a100)",iostat=io_error) line
            if (io_error.ne.0) exit
            lk=INDEX(line,trim(shortpdfname))
            if (lk .ne. 0) then
               ! 1st entry of line is PDF number with trailing whitespace
               li=INDEX(line," ")
               line=line(li+1:len_trim(line))
               lf=INDEX(line," ")
               line=line(1:lf-1)
               read(line,*) value
               exit
            end if
         enddo
         ! is this a single-member pdf set?
         if (value.ne.0) then
            rewind (unit=iunit, iostat=io_error)
            if (io_error.eq.0) then
               isSingleMemberPDF = .true.
               do
                  read(iunit,"(a100)",iostat=io_error) line
                  if (io_error.ne.0) exit
                  li=INDEX(line," ")
                  line=line(li+1:len_trim(line))
                  lf=INDEX(line," ")
                  templine=line(1:lf-1)
                  if (trim(templine).ne."") then
                     read(templine,*) nextvalue
                  else
                     nextvalue=0
                  endif
                  ! check if lhaglue+1 still corresponds to the same pdf set
                  if (nextvalue.eq.(value+1)) then
                     lk=INDEX(line,trim(shortpdfname))
                     if (lk .ne. 0) then
                        ! lhaglue+1 still corresponds to the same pdf set
                        isSingleMemberPDF = .false.
                        exit
                     endif
                  end if
               enddo
            endif
         endif
         close(iunit)
      endif
#endif

      return
      end


#endif

