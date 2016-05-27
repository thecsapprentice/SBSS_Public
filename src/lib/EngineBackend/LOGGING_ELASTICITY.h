#include "HYBRID_NONLINEAR_ELASTICITY.h"
#include <string>
namespace PhysBAM{

template<class T_ELASTICITY>
class LOGGING_ELASTICITY
{
    typedef T_ELASTICITY T_DISCRETIZATION;   

    public:

    void ExportElasticityData(T_DISCRETIZATION& discretization, const std::string& filename);
    void ImportElasticityData(const std::string& filename,
                              float& h, VECTOR<int,3>& n,
                              VECTOR<float,3>& origin,
                              float& mu,
                              float& alpha,
                              float& kappa,
                              float& cutoff_value,
                              float& stabilization_factor,
                              ARRAY<CELL_TYPE, VECTOR<int, 3> >& cell_type,
                              int& number_of_mesh_cells, 
                              int& number_of_mesh_nodes,
                              ARRAY<CELL_TYPE>& cell_type_mesh,
                              ARRAY<VECTOR<int,3> >& cell_indices_mesh,
                              ARRAY<VECTOR<int,8> >& cells_mesh);

};


};
