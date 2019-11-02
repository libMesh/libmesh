# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([PASSMESS_COMPILER_FEATURES],
[
AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(------- Configuring optional features -------)
AC_MSG_RESULT(---------------------------------------------)


# --------------------------------------------------------------
# library deprecated code - enable by default
# --------------------------------------------------------------
AC_ARG_ENABLE(deprecated,
              [AS_HELP_STRING([--disable-deprecated],[Deprecated code use gives errors rather than warnings])],
              enabledeprecated=$enableval,
              enabledeprecated=yes)

AC_SUBST(enabledeprecated)
AS_IF([test "$enabledeprecated" != yes],
      [
        AC_MSG_RESULT([>>> INFO: Disabling library deprecated code <<<])
        AC_MSG_RESULT([>>> Configuring library without deprecated code support <<<])
      ],
      [
        AC_MSG_RESULT([<<< Configuring library with deprecated code support >>>])
        AC_DEFINE(ENABLE_DEPRECATED, 1, [Flag indicating if the library should support deprecated code])
      ])
# --------------------------------------------------------------


# --------------------------------------------------------------
# C++ exceptions - enabled by default
# --------------------------------------------------------------
AC_ARG_ENABLE(exceptions,
              AS_HELP_STRING([--disable-exceptions],
                             [exit rather than throw exceptions on unexpected errors]),
              enableexceptions=$enableval,
              enableexceptions=yes)

AS_IF([test "$enableexceptions" != no],
      [
        AC_DEFINE(ENABLE_EXCEPTIONS, 1, [Flag indicating if the library should be built to throw C++ exceptions on unexpected errors])
        AC_MSG_RESULT(<<< Configuring library with exception throwing support >>>)
      ])
# --------------------------------------------------------------


# --------------------------------------------------------------
# __TIME__ __DATE__ stamps - enabled by default
# disabling preprocessor timestamps helps compiler caches such
# as ccache to work more effectively.
# --------------------------------------------------------------
AC_ARG_ENABLE(timestamps,
              AS_HELP_STRING([--disable-timestamps],
                             [do not add preprocessor timestamps to the library (helps ccache)]),
              enabletimestamps=$enableval,
              enabletimestamps=yes)

AS_IF([test "$enabletimestamps" != no],
      [
        AC_DEFINE(ENABLE_TIMESTAMPS, 1, [Flag indicating if the library should be built with compile time and date timestamps])
        AC_MSG_RESULT(<<< Configuring library with compile timestamps >>>)
      ])

# --------------------------------------------------------------


# --------------------------------------------------------------
# library warnings - enable by default
# --------------------------------------------------------------
AC_ARG_ENABLE(warnings,
              [AS_HELP_STRING([--disable-warnings],[Do not warn about deprecated, experimental, or questionable code])],
              enablewarnings=$enableval,
              enablewarnings=yes)

AC_SUBST(enablewarnings)
AS_IF([test "$enablewarnings" != yes],
      [
        AC_MSG_RESULT([>>> INFO: Disabling library warnings <<<])
        AC_MSG_RESULT([>>> Configuring library without warnings <<<])
      ],
      [
        AC_MSG_RESULT([<<< Configuring library with warnings >>>])
        AC_DEFINE(ENABLE_WARNINGS, 1, [Flag indicating if the library should have warnings enabled])
      ])
# --------------------------------------------------------------


AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(----- Done configuring optional features ----)
AC_MSG_RESULT(---------------------------------------------)
])
