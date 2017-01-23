AC_DEFUN([VBFNLO_DEPENDENCY_TRACKING],[
dnl Add Fortran 77 dependency checking
_AM_IF_OPTION([no-dependencies],,
[AC_PROVIDE_IFELSE([AC_PROG_F77],
                  [_AM_DEPENDENCIES([F77])],
                  [m4_define([AC_PROG_F77],
                             m4_defn([AC_PROG_F77])[_AM_DEPENDENCIES([F77])])])dnl
])
dnl Add Fortran dependency checking
_AM_IF_OPTION([no-dependencies],,
[AC_PROVIDE_IFELSE([AC_PROG_FC],
                  [_AM_DEPENDENCIES([FC])],
                  [m4_define([AC_PROG_FC],
                             m4_defn([AC_PROG_FC])[_AM_DEPENDENCIES([FC])])])dnl
])
])
