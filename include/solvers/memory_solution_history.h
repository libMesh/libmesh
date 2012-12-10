// Subclass of Solution History that stores the solutions 
// and other important vectors in memory

#ifndef __memory_solution_history_h__
#define __memory_solution_history_h__

// Local includes
#include "libmesh/numeric_vector.h"
#include "libmesh/solution_history.h"
#include <list>

namespace libMesh
{
  
  // MemorySolutionHistory class declaration and definition
  class MemorySolutionHistory : public SolutionHistory
  {
  public:    
    
    // Constructor, reference to system to be passed by user, set the  
    // stored_sols iterator to some initial value
  MemorySolutionHistory(System & system_) : stored_sols(stored_solutions.end()), _system(system_) {} ;

    // Destructor
    ~MemorySolutionHistory();

    // Virtual function store which we will be overriding to store timesteps
    virtual void store();

    // Virtual function retrieve which we will be overriding to retrieve timesteps
    virtual void retrieve();

    // Typedef for Stored Solutions iterator, a list of pairs of the current 
    // system time, map of strings and saved vectors
    typedef std::list<std::pair<Real, std::map<std::string, NumericVector<Number>*> > >::iterator stored_solutions_iterator;

    // Definition of the clone function needed for the setter function 
    virtual AutoPtr<SolutionHistory > clone() const {
    return AutoPtr<SolutionHistory > 
      (new MemorySolutionHistory(_system));}
    
  private:
    
    // This list of pairs will hold the current time and stored vectors 
    // from each timestep
    std::list<std::pair<Real, std::map<std::string, NumericVector<Number>*> > > stored_solutions; 
    
    // The stored solutions iterator
    stored_solutions_iterator stored_sols;
    
    // A system reference
    System & _system ;

  };

} // end namespace libMesh

#endif // __memory_solution_history_h__
