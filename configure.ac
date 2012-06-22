dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(epic, 1.0)


AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE()

AC_PROG_CC
#AC_PROG_CXX

m4_include([m4/m4_ax_prog_cxx_mpi.m4])
# If --with-mpi=auto is used, try to find MPI, but use standard C compiler if it is not found.
# If --with-mpi=yes is used, try to find MPI and fail if it isn't found.
# If --with-mpi=no is used, use a standard C compiler instead.
AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
	[compile with MPI (parallelization) support. If none is found,
	MPI is not used. Default: auto])
	],,[with_mpi=auto])
AX_PROG_CXX_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
	use_mpi=no
	if test x"$with_mpi" = xyes; then
		AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
	else
		AC_MSG_WARN([No MPI compiler found, won't use MPI.])
	fi
	])
AM_CONDITIONAL([USE_MPI], [test x"$use_mpi" = xyes])

AC_PROG_RANLIB

#m4_include([m4/m4_ax_mpi.m4])
#AC_LANG_PUSH(C++)
#AX_MPI
#AX_PROG_CXX_MPI
#AC_LANG_POP(C++)
AC_CONFIG_FILES(Makefile src/Makefile fmdb/Makefile eigen/Makefile \
cuba/Makefile cubature/Makefile boost/Makefile odeint/Makefile vtk/Makefile \
cmake/Makefile)
AC_OUTPUT

