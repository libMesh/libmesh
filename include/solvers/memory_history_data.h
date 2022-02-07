#ifndef MEMORY_HISTORY_DATA_H
#define MEMORY_HISTORY_DATA_H

// Local includes
#include "libmesh/history_data.h"
#include "libmesh/diff_system.h"

#include "libmesh/numeric_vector.h"

namespace libMesh
{
    /** MemoryHistoryData provides a data structure to store memory history data.
     * This is a companion class to MemorySolutionHistory.
     */

    class MemoryHistoryData : public HistoryData
    {
        public:

        // Constructor
        MemoryHistoryData(DifferentiableSystem & system) : HistoryData(), _system(system), stored_vecs{}, stored_vec(stored_vecs.end()) {};

        // Destructor
        ~MemoryHistoryData() {};

        virtual void store_initial_solution() override;
        virtual void store_primal_solution(stored_data_iterator stored_datum) override;
        virtual void store_adjoint_solution() override;
        virtual void rewrite_stored_solution() override;

        virtual void retrieve_primal_solution() override;
        virtual void retrieve_adjoint_solution() override;

        void store_vectors();
        void retrieve_vectors();

        private:

        DifferentiableSystem & _system;

        typedef std::map<std::string, std::unique_ptr<NumericVector<Number>>> map_type;
        typedef map_type::iterator stored_vecs_iterator;

        // Memory History specific Constituents
        //std::unique_ptr<MeshBase> stored_mesh;
        map_type stored_vecs;
        stored_vecs_iterator stored_vec;

    };

}
#endif
