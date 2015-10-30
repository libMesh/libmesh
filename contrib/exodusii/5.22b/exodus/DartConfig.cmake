# Dashboard is opened for submissions for a 24 hour period starting at
# the specified NIGHLY_START_TIME. Time is specified in 24 hour format.
SET (NIGHTLY_START_TIME "21:00:00 EDT")

# Dart server to submit results (used by client)
SET (DROP_SITE "tigre.ca.sandia.gov")
SET (DROP_LOCATION "Drop/")
SET (DROP_SITE_USER "dart")
SET (DROP_METHOD "scp")
SET (TRIGGER_SITE "http://${DROP_SITE}/~dart/cgi-bin/Dart-vtkSNL.pl")

# Project Home Page
SET (PROJECT_URL "http://tigre.ca.sandia.gov/vtkSNL")
SET (CVS_WEB_URL "http://${DROP_SITE}/vtk/")
SET (CVS_WEB_CVSROOT "vtkSNL")
SET (USE_DOXYGEN "On")
SET (DOXYGEN_URL "http://tigre.ca.sandia.gov/vtkSNL/Documentation/Doxygen/" )
