	integer Maxlines
        parameter (maxlines=200)
        character*250 lines(maxlines)
        integer n_lines
        integer iunit
	character*250 fname
	logical usedefaults,showinput
        COMMON/READINPUT/iunit,n_lines,usedefaults,showinput,lines,fname

        character*250 inputpath, pdfpath
        COMMON/INPUTPATHES/inputpath, pdfpath
