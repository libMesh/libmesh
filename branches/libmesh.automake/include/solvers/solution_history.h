// A SolutionHistory class that enables the storage and retrieval of timesteps 
// and (in the future) adaptive steps

#ifndef __solution_history_h__
#define __solution_history_h__

// Local Includes
#include "libmesh/system.h"

namespace libMesh
{

  // SolutionHistory class declaration and definition
  class SolutionHistory
  {
  public:
    
    // Constructor
    SolutionHistory() {};
    
    // Destructor
    ~SolutionHistory () {};

    // Function to store a solution, pure virtual
    virtual void store() = 0;

    // Function to retrieve a solution, pure virtual
    virtual void retrieve() = 0;

    // Cloning function for an AutoPtr, pure virtual, used in the 
    // setter function in time_solver.C
    virtual AutoPtr<SolutionHistory > clone() const = 0;
    
  }; // end SolutionHistory class definition

} // end namespace libMesh

#endif // __solution_history_h__
