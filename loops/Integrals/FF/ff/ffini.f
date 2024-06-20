*	$Id: vbfffini.f,v 1.1 1996/03/27 08:05:18 gj Exp $
*
*	glue routine for older versions of FF:
*	define vbfffinit, vbfffexit to be equal to vbfffini, vbfffexi
*	
*	when using CERN libs do *not* include this file, because vbfffinit 
*	already exists in packlib.
*
*	All programs written after 17-mar-1996 should work without this file
*	
 	subroutine vbfffinit
	call vbfffini
	end
 	subroutine vbfffexit
	call vbfffexi
	end
