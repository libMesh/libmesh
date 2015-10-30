#include <cstdlib>
#include "TECIO.h"

namespace libMesh {
  namespace totally_obscure_tecplot_namespace {
    
    void legacy_api()
    {
      // Call legacy API functions to ensure they get imported.
      TECINI (NULL,		 
	      NULL,
	      NULL,
	      NULL,
	      NULL,
	      NULL);

      TECZNE (NULL,
	      NULL,
	      NULL,
	      NULL,
	      NULL,
	      NULL);

      TECDAT (NULL,
	      NULL,
	      NULL);

      TECNOD (NULL);
      
      TECEND ();
    }

    void v112_api()
    {
      // Call v112 API functions to ensure they get imported.
      TECINI112 (NULL,		 
		 NULL,
		 NULL,
		 NULL,
		 NULL,
		 NULL,
		 NULL);

      TECZNE112 (NULL,
		 NULL,
		 NULL,
		 NULL,
		 NULL,
		 NULL,
		 NULL,
		 NULL,
		 NULL,
		 NULL,
		 NULL,
		 NULL,
		 NULL,
		 NULL,
		 NULL,
		 NULL,
		 NULL,
		 NULL,
		 NULL,
		 NULL,
		 NULL);

      TECDAT112 (NULL,
		 NULL,
		 NULL);

      TECNOD112 (NULL);
      
      TECEND112 ();
    }
  }
}

