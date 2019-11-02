# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([PASSMESS_COMPILER_FEATURES],
[
AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(------- Configuring optional features -------)
AC_MSG_RESULT(---------------------------------------------)

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
AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(----- Done configuring optional features ----)
AC_MSG_RESULT(---------------------------------------------)
])
