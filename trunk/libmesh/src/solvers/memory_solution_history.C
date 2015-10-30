// Function definitions for MemorySolutionHistory

// Local includes
#include "libmesh/memory_solution_history.h"

namespace libMesh
{
  // The Destructor
  MemorySolutionHistory::~MemorySolutionHistory ()
  {
    stored_solutions_iterator stored_sols_it = stored_solutions.begin();
    const stored_solutions_iterator stored_sols_end = stored_solutions.end();

    for(; stored_sols_it != stored_sols_end; ++stored_sols_it)
      {
	// The saved vectors at this timestep
	std::map<std::string, NumericVector<Number> *> saved_vectors = stored_sols_it->second;

	std::map<std::string, NumericVector<Number> *>::iterator vec = saved_vectors.begin();
	std::map<std::string, NumericVector<Number> *>::iterator vec_end = saved_vectors.end();

	// Loop over all the saved vectors
	for (; vec != vec_end; ++vec)
	  {		    	
	    // Delete this saved vector
	    delete vec->second;	    
	  }
      }
  }

  // This functions saves all the 'projection-worthy' system vectors for
  // future use
  void MemorySolutionHistory::store()
  {    
    // Map of to be projected vectors from this solution step
    std::map<std::string, NumericVector<Number> *> saved_vectors;
    
    // Loop over all the system vectors
    for (System::vectors_iterator vec = _system.vectors_begin(); vec != _system.vectors_end(); ++vec)
      {	
	// The name of this vector
	const std::string& vec_name = vec->first;

	// Check if this vector is to be projected
	bool vector_projection_setting = _system.vector_preservation(vec_name);
	
	// If it is important enough to be projected, it is important enough to be saved
	if(vector_projection_setting)
	  {
	    saved_vectors[vec_name] = vec->second->clone().release();
	  }

      }
    
    // Of course, we will always save the actual solution
    std::string _solution("_solution");    
    saved_vectors[_solution] = _system.solution->clone().release();

    // Put all the vectors from this timestep into stored_solutions
    stored_solutions.push_back(std::make_pair(_system.time, saved_vectors));
    
  }

  void MemorySolutionHistory::retrieve()
  {    
    // To find the required entry in the stored_solutions list, we need to first decrement the stored_sols iterator
    --stored_sols;
    
    // Get the time at which we are recovering the solution vectors
    Real recovery_time = stored_sols->first;

    // Print out what time we are recovering vectors at
    std::cout<<"Recovering solution vectors at time: " <<recovery_time << std::endl;

    // Get the saved vectors at this timestep
    std::map<std::string, NumericVector<Number> *> saved_vectors = stored_sols->second;

    std::map<std::string, NumericVector<Number> *>::iterator vec = saved_vectors.begin();
    std::map<std::string, NumericVector<Number> *>::iterator vec_end = saved_vectors.end();

    // Loop over all the saved vectors
    for (; vec != vec_end; ++vec)
      {	
  	// The name of this vector
  	const std::string& vec_name = vec->first;
	
	// Get the vec_name entry in the saved vectors map and set the current system vec[vec_name] entry to it
	if(vec_name != "_solution")
	  _system.get_vector(vec_name) = *(vec->second);	    
      }
    
    // Of course, we will *always* have to get the actual solution
    std::string _solution("_solution");    
    *(_system.solution) = *(saved_vectors[_solution]);
        
  }

}
