<?php

function detect() {
  // Note: I originally found this script at: http://us2.php.net/manual/en/function.get-browser.php

  // Temporary Variables

  // The useragent string (lowercase to simplify testing)
  $_nw_ua = strtolower(@$_SERVER["HTTP_USER_AGENT"]);

  // Browser Detection { ======================================================

  // Version checking, each one of these will take a float value describing the
  // version number, or - if the user is not using that browser - zero.

  // Generic code-name "Mozilla" version
  define("NW_MOZ_VERSION", preg_match('/mozilla\/(\d+\.\d+)/',
				      $_nw_ua, $_nw_v) ? (float)$_nw_v[1] : 0);

  // KDE's Konqueror
  define("NW_IS_KONQ", preg_match('/konqueror\/(\d+\.\d+)/',
				  $_nw_ua, $_nw_v) ? (float) $_nw_v[1] : 0);

  // Opera software Opera
  define("NW_IS_OPERA", preg_match('/opera[\s\/](\d+\.\d+)/',
				   $_nw_ua, $_nw_v) ? (float) $_nw_v[1] : 0);

  // Microsoft Internet Explorer
  define("NW_IS_IE", !NW_IS_OPERA && preg_match('/msie (\d+\.\d+)/',
						$_nw_ua, $_nw_v) ? (float) $_nw_v[1] : 0);

  // Gecko-based browsers, such as Mozilla, Netscape 6, DocZilla,
  // K-Meleon, etc.
  define("NW_IS_GECKO", preg_match('/gecko\/(\d+)/',
				   $_nw_ua, $_nw_v) ? (float) $_nw_v[1] : 0);

  // Netscape Navigator (all versions, including Gecko-based browsers)
  define("NW_IS_NN", NW_IS_GECKO ? (preg_match('/netscape6*\/(\d+.\d+)/', $_nw_ua, $_nw_v) ?
				    (float) $_nw_v[1] : 0) : ((!NW_IS_OPERA && !NW_IS_KONQ && !NW_IS_IE) ?
							      NW_MOZ_VERSION : 0));

  // An old 3rd generation web browser
  define("NW_IS_GEN3", NW_IS_NN < 4 && NW_IS_OPERA < 4 && NW_IS_IE < 4 && NW_MOZ_VERSION < 4);
  // } Browser Detection ======================================================

  // Generic Platform Detection { =============================================
  define("NW_IS_LINUX", strstr($_nw_ua, "linux") !== false);
  define("NW_IS_MAC", strstr($_nw_ua, "mac") !== false);
  define("NW_IS_SOLARIS", (strstr($_nw_ua, "solaris") !== false) ||
	 (strstr($_nw_ua, "sunos") !== false));
  define("NW_IS_X11", strstr($_nw_ua, "x11") !== false);
  define("NW_IS_WINDOWS", strstr($_nw_ua, "win") !== false);
  define("NW_IS_OS2", strstr($_nw_ua, "os2") !== false);
  // } Generic Platform Detection =============================================

  unset($_nw_ua, $_nw_v); // clean-up
}

?>