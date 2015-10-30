// 'Save nothing' subclass of Solution History, this is the default

#ifndef __no_solution_history_h__
#define __no_solution_history_h__

// Local includes
#include "libmesh/solution_history.h"

namespace libMesh
{
  
  // NoSolutionHistory class declaration and definition
  class NoSolutionHistory : public SolutionHistory
  {
  public:
    
    // Constructor
  NoSolutionHistory() : SolutionHistory() {};

    // Destructor
    virtual ~NoSolutionHistory() {};

    // Virtual function store which we will be overriding
    virtual void store();

    // Virtual function retrieve which we will be overriding
    virtual void retrieve();

    // Definition of the clone function needed for the setter function 
    virtual AutoPtr<SolutionHistory > clone() const {
    return AutoPtr<SolutionHistory > 
      (new NoSolutionHistory());}
    
  }; // end NoSolutionHistory class definition

} // end namespace libMesh

#endif // __no_solution_history_h__
