AC_INIT([mpMap2])
AC_ARG_WITH([boost],
            AC_HELP_STRING([--with-boost],
                           [use boost]),
            []
            [with_boost=yes])
if test "x$with_boost" = "xyes" ; then
	PKG_CPPFLAGS=-DUSE_BOOST
	PKG_CFLAGS=-DUSE_BOOST
fi
AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_CFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
