# -------------------------------------------------------------
# MetaPhysicL - a C++ header-only numerics utility library
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_METAPHYSICL],
[
  AC_ARG_ENABLE(metaphysicl,
                AS_HELP_STRING([--disable-metaphysicl],
                               [build without MetaPhysicL suppport]),
                [case "${enableval}" in
                  yes)  enablemetaphysicl=yes ;;
                  no)  enablemetaphysicl=no ;;
                  *)  AC_MSG_ERROR(bad value ${enableval} for --enable-metaphysicl) ;;
                esac],
                [enablemetaphysicl=$enablenested])

  # The MetaPhysicL API is distributed with libmesh, so we don't have
  # to guess where it might be installed.  This needs to be replaced
  # someday with an option to include an external version instead.
  if (test $enablemetaphysicl = yes); then
     METAPHYSICL_INCLUDE="-I\$(top_srcdir)/contrib/metaphysicl/0.2.0/src/numerics/include -I\$(top_srcdir)/contrib/metaphysicl/0.2.0/src/core/include -I\$(top_srcdir)/contrib/metaphysicl/0.2.0/src/utilities/include"
     AC_DEFINE(HAVE_METAPHYSICL, 1, [Flag indicating whether the library will be compiled with MetaPhysicL support])
     AC_MSG_RESULT(<<< Configuring library with MetaPhysicL support >>>)
     AC_CONFIG_SUBDIRS([contrib/metaphysicl/0.2.0])
  else
     METAPHYSICL_INCLUDE=""
     enablemetaphysicl=no
  fi

  AC_SUBST(METAPHYSICL_INCLUDE)
])
