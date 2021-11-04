#ifndef LIBMESH_FILE_HISTORY_DATA_H
#define LIBMESH_FILE_HISTORY_DATA_H

#include "libmesh/history_data.h"
#include "libmesh/diff_system.h"

namespace libMesh
{

    /** HistoryData subclass that provides a struct to store history data
     * such as timestamps, mesh, primal and adjoint filenames and timestep sizes.
     *
     * \author Vikram Garg
     * \date 2021
     * \brief
     */

    class FileHistoryData : public HistoryData
    {
        public:

        FileHistoryData(DifferentiableSystem & system) : HistoryData(), _system(system),mesh_filename(""), primal_filename(""), adjoint_filename("") {};

        ~FileHistoryData() {};

        // Accessors for FileHistory specific variables
        std::string & get_mesh_filename()
        { return mesh_filename; }

        std::string & get_primal_filename()
        { return primal_filename; }

        std::string & get_adjoint_filename()
        { return adjoint_filename; }

        void set_mesh_filename(std::string & mesh_name)
        { mesh_filename = mesh_name; }

        void set_primal_filename(std::string & primal_sol_name)
        { primal_filename = primal_sol_name; }

        void set_adjoint_filename(std::string & adjoint_sol_name)
        { adjoint_filename = adjoint_sol_name; }

        virtual void store_initial_solution() override;
        virtual void store_primal_solution(stored_data_iterator stored_datum) override;
        virtual void store_adjoint_solution() override;
        virtual void rewrite_stored_solution() override;

        virtual void retrieve_primal_solution() override;
        virtual void retrieve_adjoint_solution() override;

        private:

        // Reference to underlying system
        DifferentiableSystem & _system;

        // File History specific variables
        std::string mesh_filename;
        std::string primal_filename;
        std::string adjoint_filename;

    };
}
#endif