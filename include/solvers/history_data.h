#ifndef LIBMESH_HISTORY_DATA_H
#define LIBMESH_HISTORY_DATA_H

#include <cmath>
#include "libmesh/system.h"

// LOCAL INCLUDES

namespace libMesh
{

    /** The History Data classes are companion classes to SolutionHistory and MeshHistory classes.
     * These provide data structures to store different types of history data (timestamps,
     * pointers, filenames) depending on the type of History being used.
     *
     * \author Vikram Garg
     * \date 2021
     * \brief Provides history data structures and I/O, memory operation functions.
     */

    class HistoryData
    {
        public:

        // Constructor
        HistoryData() : time_stamp(std::numeric_limits<unsigned int>::signaling_NaN()),
        deltat_at(std::numeric_limits<double>::signaling_NaN()), previously_stored(false) {};

        // Destructor
        virtual ~HistoryData() {};

        // Accessors for individual history members, common to all types of histories.
        unsigned int get_time_stamp()
        { return time_stamp; }

        Real get_deltat_at()
        { return deltat_at; }

        bool get_previously_stored()
        { return previously_stored; }

        // Setters for common history members
        void set_time_stamp(unsigned int time_stamp_val)
        { time_stamp = time_stamp_val; }

        void set_deltat_at(Real deltat_at_val)
        { deltat_at = deltat_at_val; }

        void set_previously_stored(bool previously_stored_val)
        { previously_stored = previously_stored_val; }

        // Some history storing functions take an iterator to a map as an argument.
        // This iterator is necessary for setting time stamps and delta_ts.
        typedef std::map<Real, std::unique_ptr<HistoryData>> map_type;
        typedef map_type::iterator stored_data_iterator;

        // Operations for SolutionHistory
        virtual void store_initial_solution() = 0;
        virtual void store_primal_solution(stored_data_iterator stored_datum) = 0;
        virtual void store_adjoint_solution() = 0;
        virtual void rewrite_stored_solution() = 0;

        virtual void retrieve_primal_solution() = 0;
        virtual void retrieve_adjoint_solution() = 0;

        // Operations for MeshHistory

        protected:

        // Variable members common to all derived HistoryData types

        // The index of the current time step
        unsigned int time_stamp;

        // The delta_t (timestep taken) at the current timestep.
        Real deltat_at;

        // To help check if we are storing fresh or overwriting.
        bool previously_stored;

    };
} // end namespace libMesh

#endif // LIBMESH_HISTORY_DATA_H