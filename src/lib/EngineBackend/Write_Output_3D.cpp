#if 0
#include "OVERRIDES.h"

#ifdef ENABLE_PHYSBAM_IO
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>

#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>

#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>


#include "ELASTIC_LATTICE_DEFORMER.h"
#include "KEYFRAMED_PARAMETER_ANIMATION.h"
//#include "EMBEDDED_DEFORMER.h"
#include "CLElib.h"
#include "CORE_ENGINE_WRAPPER.h"
#include "NONLINEAR_ELASTICITY.h"
#include "HYBRID_NONLINEAR_ELASTICITY.h"
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include "CG_POLICY.h"
#include "Write_Output.h"

namespace PhysBAM{

//#define DISPLAY_FINE_GRID
    namespace {
//#####################################################################
// Function Fill_Mesh_Particles
//#####################################################################

    template<class T,int d,class T_DISCRETIZATION,bool mesh>
    struct Fill_Mesh_Particles{
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        void operator() (const T_DISCRETIZATION& discretization,
                         const typename T_DISCRETIZATION::T_STATE& state,
                         HASHTABLE<int,int>& particle_mesh_hash,
                         GEOMETRY_PARTICLES<TV>& particles){}
    };
    
    
    template<class T,int d,class T_DISCRETIZATION>
    struct Fill_Mesh_Particles<T,d,T_DISCRETIZATION,true>{
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        void operator() (const T_DISCRETIZATION& discretization,
                         const typename T_DISCRETIZATION::T_STATE& state,
                         HASHTABLE<int,int>& particle_mesh_hash,
                         GEOMETRY_PARTICLES<TV>& particles){
            // Initialize particles on mesh
            for(int m=1;m<=discretization.Number_Of_Mesh_Cells();m++){
                const T_INDEX& cell_index = discretization.cell_indices_mesh(m);
                int v=1;
                for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),v++){        
                    if( discretization.cells_mesh(m)(v) > 0 && !particle_mesh_hash.Contains(discretization.cells_mesh(m)(v)) ){
                        int p = particles.array_collection->Add_Element();
                        particle_mesh_hash.Insert(discretization.cells_mesh(m)(v),p);
                        for(int i=1;i<=d;i++){
                            particles.X(p)(i) =state.x_mesh(i)(discretization.cells_mesh(m)(v))+discretization.grid.Node(iterator.Index())(i);
                        }}}}}};                             
        
//#####################################################################
// Function Create Mesh Tetrahedralized Volume
//#####################################################################

    template<class T,int d,class T_DISCRETIZATION,bool mesh>
    struct Create_Mesh_Tetrahedralized_Volume{
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        void operator() (T_DISCRETIZATION& discretization,
                         HASHTABLE<T_INDEX,int>& particle_hash,
                         HASHTABLE<int,int>& particle_mesh_hash,
                         TETRAHEDRALIZED_VOLUME<T>& interior_mesh_tetrahedralized_volume){}
    };
    
    
    template<class T,int d,class T_DISCRETIZATION>
    struct Create_Mesh_Tetrahedralized_Volume<T,d,T_DISCRETIZATION,true>{
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        void operator() (const T_DISCRETIZATION& discretization,
                         HASHTABLE<T_INDEX,int>& particle_hash,
                         HASHTABLE<int,int>& particle_mesh_hash,
                         TETRAHEDRALIZED_VOLUME<T>& interior_mesh_tetrahedralized_volume){
            for(int m=1;m<=discretization.Number_Of_Mesh_Cells();m++){
                const T_INDEX& cell_index = discretization.cell_indices_mesh(m);
                int i,j,ij;cell_index.Get(i,j,ij);    
                int v=1;
                int corners[8];
                for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));
                    iterator.Valid();iterator.Next(),v++){
                    if(discretization.cells_mesh(m)(v) > 0)
                        corners[v-1] = particle_mesh_hash.Get(discretization.cells_mesh(m)(v));
                    else
                        corners[v-1] = particle_hash.Get(iterator.Index());
                }
                if((i+j+ij)%2==0){
                    interior_mesh_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(corners[0],corners[4],
                                                                                            corners[2],corners[1]));
                    interior_mesh_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(corners[4],corners[5],
                                                                                            corners[7],corners[1]));
                    interior_mesh_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(corners[2],corners[6],
                                                                                            corners[7],corners[4]));
                    interior_mesh_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(corners[3],corners[7],
                                                                                            corners[1],corners[2]));
                    interior_mesh_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(corners[4],corners[1],
                                                                                            corners[7],corners[2]));}
                else{
                    interior_mesh_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(corners[0],corners[4],
                                                                                            corners[6],corners[5]));
                    interior_mesh_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(corners[0],corners[2],
                                                                                            corners[3],corners[6]));
                    interior_mesh_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(corners[3],corners[5],
                                                                                            corners[1],corners[0]));
                    interior_mesh_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(corners[3],corners[7],
                                                                                            corners[5],corners[6]));
                    interior_mesh_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(corners[3],corners[0],
                                                                                            corners[6],corners[5]));}        
            }
            LOG::cout << "Created Mesh Volume with " << interior_mesh_tetrahedralized_volume.mesh.elements.m/5 << " elements." << std::endl;
        }};

//#####################################################################
// Function Create_Mesh_Constraint
//#####################################################################

    template<class T,int d,class T_DISCRETIZATION,bool mesh>
    struct Create_Mesh_Constraint{
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        int operator() (const T_DISCRETIZATION& discretization,
                        const typename T_DISCRETIZATION::T_STATE& state,
                        GEOMETRY_PARTICLES<TV>& particles,
                        const CONSTRAINT_NODE<T,d>& endpoint,
                        FREE_PARTICLES<TV>& free_particles)
        {
            PHYSBAM_NOT_IMPLEMENTED();
            return 1;
        }
    };
    
    
    template<class T,int d,class T_DISCRETIZATION>
    struct Create_Mesh_Constraint<T,d,T_DISCRETIZATION,true>{
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        int operator() (const T_DISCRETIZATION& discretization,
                        const typename T_DISCRETIZATION::T_STATE& state,
                        GEOMETRY_PARTICLES<TV>& particles,
                        const CONSTRAINT_NODE<T,d>& endpoint,
                        FREE_PARTICLES<TV>& free_particles)
        {
            int e = particles.array_collection->Add_Element();
            particles.X(e)=discretization.Deformation_Mesh(endpoint.mesh_index(),
                                                           endpoint.multilinear_coordinates(),
                                                           T_DISCRETIZATION::BASE::View_Convert(state.x),
                                                           T_DISCRETIZATION::View_Convert(state.x_mesh));
            free_particles.nodes.Append(e);
            return e;
        }
    };

//#####################################################################
// Function Create Mesh Tetrahedralized Volume
//#####################################################################

    template<class T,int d,class T_DISCRETIZATION,bool mesh>
    struct Create_Mesh_Only_Variables{
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        void operator() (const T_DISCRETIZATION& discretization,
                         HASHTABLE<T_INDEX,int>& particle_hash,
                         HASHTABLE<int,int>& particle_mesh_hash,
                         FREE_PARTICLES<TV>& mesh_only_variables){}
    };

    template<class T,int d,class T_DISCRETIZATION>
    struct Create_Mesh_Only_Variables<T,d,T_DISCRETIZATION,true>{
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        void operator() (const T_DISCRETIZATION& discretization,
                         HASHTABLE<T_INDEX,int>& particle_hash,
                         HASHTABLE<int,int>& particle_mesh_hash,
                         FREE_PARTICLES<TV>& mesh_only_variables){
            for(int m=1;m<=discretization.Number_Of_Mesh_Cells();m++){
                const T_INDEX& cell_index = discretization.cell_indices_mesh(m);
                int i,j,ij;cell_index.Get(i,j,ij);    
                int v=1;
                int corners[8];
                for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));
                    iterator.Valid();iterator.Next(),v++){
                    if(discretization.cells_mesh(m)(v) > 0)
                        mesh_only_variables.nodes.Append(particle_mesh_hash.Get(discretization.cells_mesh(m)(v)));
                }
            }
        }
    };

//#####################################################################
// Function Create Mesh Tetrahedralized Volume
//#####################################################################

    template<class T,int d,class T_DISCRETIZATION,bool mesh>
    struct Create_Mesh_Shared_Variables{
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        void operator() (const T_DISCRETIZATION& discretization,
                         HASHTABLE<T_INDEX,int>& particle_hash,
                         HASHTABLE<int,int>& particle_mesh_hash,
                         FREE_PARTICLES<TV>& mesh_shared_variables){}
    };

    template<class T,int d,class T_DISCRETIZATION>
    struct Create_Mesh_Shared_Variables<T,d,T_DISCRETIZATION,true>{
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        void operator() (const T_DISCRETIZATION& discretization,
                         HASHTABLE<T_INDEX,int>& particle_hash,
                         HASHTABLE<int,int>& particle_mesh_hash,
                         FREE_PARTICLES<TV>& mesh_shared_variables){
            for(int m=1;m<=discretization.Number_Of_Mesh_Cells();m++){
                const T_INDEX& cell_index = discretization.cell_indices_mesh(m);
                int i,j,ij;cell_index.Get(i,j,ij);    
                int v=1;
                int corners[8];
                for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));
                    iterator.Valid();iterator.Next(),v++){
                    if(discretization.cells_mesh(m)(v) == 0)
                        mesh_shared_variables.nodes.Append(particle_hash.Get(iterator.Index()));
                }
            }
        }
    };

#if 0    

//#####################################################################
// Function Create Mesh Tetrahedralized Volume
//#####################################################################

    template<class T,int d,class T_DISCRETIZATION,bool mesh>
    struct Create_Boundary_Variables{
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        void operator() (T_DISCRETIZATION& discretization,
                         HASHTABLE<T_INDEX,int>& particle_hash,
                         HASHTABLE<int,int>& particle_mesh_hash,
                         FREE_PARTICLES<TV>& boundary_variables){}
    };

    template<class T,int d,class T_DISCRETIZATION>
    struct Create_Boundary_Variables<T,d,T_DISCRETIZATION,false>{
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        void operator() (T_DISCRETIZATION& discretization,
                         HASHTABLE<T_INDEX,int>& particle_hash,
                         HASHTABLE<int,int>& particle_mesh_hash,
                         FREE_PARTICLES<TV>& boundary_variables){

            for(RANGE_ITERATOR<d> iterator(discretization.unpadded_node_domain);iterator.Valid();iterator.Next()){
                if( discretization.node_is_boundary(iterator.Index()) == true ){
                    boundary_variables.nodes.Append( particle_hash.Get( iterator.Index() ) );
                }
            }
        }
    };

    template<class T,int d,class T_DISCRETIZATION>
    struct Create_Boundary_Variables<T,d,T_DISCRETIZATION,true>{
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        void operator() (T_DISCRETIZATION& discretization,
                         HASHTABLE<T_INDEX,int>& particle_hash,
                         HASHTABLE<int,int>& particle_mesh_hash,
                         FREE_PARTICLES<TV>& boundary_variables){

            for(RANGE_ITERATOR<d> iterator(discretization.unpadded_node_domain);iterator.Valid();iterator.Next()){
                if( discretization.node_is_boundary(iterator.Index()) == true ){
                    boundary_variables.nodes.Append( particle_hash.Get( iterator.Index() ) );
                }
            }
            for(int m=1;m<=discretization.Number_Of_Mesh_Nodes();m++){
                if( discretization.node_is_boundary_mesh(m) == true ){
                    boundary_variables.nodes.Append( particle_mesh_hash.Get( m ) );
                }
            }
        }
    };

//#####################################################################
// Function Create Subdomain Volume
//#####################################################################

    template<class T,int d,class T_DISCRETIZATION,bool mesh>
    struct Create_Subdomain_Volume{
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        void operator() (T_DISCRETIZATION& discretization,
                         HASHTABLE<T_INDEX,int>& particle_hash,
                         HASHTABLE<int,int>& particle_mesh_hash,
                         int subdomain,
                         TETRAHEDRALIZED_VOLUME<T>& subdomain_volume){}
    };
    
    template<class T,int d,class T_DISCRETIZATION>
    struct Create_Subdomain_Volume<T,d,T_DISCRETIZATION,false>{
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        void operator() (T_DISCRETIZATION& discretization,
                         HASHTABLE<T_INDEX,int>& particle_hash,
                         HASHTABLE<int,int>& particle_mesh_hash,
                         int subdomain,
                         TETRAHEDRALIZED_VOLUME<T>& subdomain_volume){

            for(RANGE_ITERATOR<d> iterator(discretization.unpadded_cell_domain);iterator.Valid();iterator.Next()){
                const T_INDEX& index=iterator.Index();
                if(discretization.cell_color(index)==subdomain && ( discretization.cell_type(index) == INTERIOR_CELL_TYPE || discretization.cell_type(index) == BOUNDARY_CELL_TYPE  || discretization.cell_type(index) == DIRICHLET_CELL_TYPE )) {
                    int i,j,ij;index.Get(i,j,ij);
                    int v=1;
                    int corners[8];
                    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(index,index+1));
                        iterator.Valid();iterator.Next(), v++){
                        corners[v-1] = particle_hash.Get(iterator.Index());
                    }

                    if((i+j+ij)%2==0){
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[0],corners[4], corners[2],corners[1]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[4],corners[5], corners[7],corners[1]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[2],corners[6], corners[7],corners[4]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[3],corners[7], corners[1],corners[2]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[4],corners[1], corners[7],corners[2]));}
                    else{
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[0],corners[4], corners[6],corners[5]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[0],corners[2], corners[3],corners[6]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[3],corners[5], corners[1],corners[0]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[3],corners[7], corners[5],corners[6]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[3],corners[0], corners[6],corners[5]));}   
                }
            }
        }
    };


    
    template<class T,int d,class T_DISCRETIZATION>
    struct Create_Subdomain_Volume<T,d,T_DISCRETIZATION,true>{
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        void operator() (T_DISCRETIZATION& discretization,
                         HASHTABLE<T_INDEX,int>& particle_hash,
                         HASHTABLE<int,int>& particle_mesh_hash,
                         int subdomain,
                         TETRAHEDRALIZED_VOLUME<T>& subdomain_volume){

            for(RANGE_ITERATOR<d> iterator(discretization.unpadded_cell_domain);iterator.Valid();iterator.Next()){
                const T_INDEX& index=iterator.Index();
                if(discretization.cell_color(index)==subdomain && ( discretization.cell_type(index) == INTERIOR_CELL_TYPE || discretization.cell_type(index) == BOUNDARY_CELL_TYPE  || discretization.cell_type(index) == DIRICHLET_CELL_TYPE )){
                    int i,j,ij;index.Get(i,j,ij);
                    int v=1;
                    int corners[8];
                    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(index,index+1));
                        iterator.Valid();iterator.Next(), v++){
                        corners[v-1] = particle_hash.Get(iterator.Index());
                    }

                    if((i+j+ij)%2==0){
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[0],corners[4], corners[2],corners[1]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[4],corners[5], corners[7],corners[1]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[2],corners[6], corners[7],corners[4]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[3],corners[7], corners[1],corners[2]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[4],corners[1], corners[7],corners[2]));}
                    else{
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[0],corners[4], corners[6],corners[5]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[0],corners[2], corners[3],corners[6]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[3],corners[5], corners[1],corners[0]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[3],corners[7], corners[5],corners[6]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[3],corners[0], corners[6],corners[5]));}   
                }
            }
            
            for(int m=1;m<=discretization.Number_Of_Mesh_Cells();m++){
                if( discretization.cell_color_mesh( m ) == subdomain && ( discretization.cell_type_mesh(m) == INTERIOR_CELL_TYPE || discretization.cell_type_mesh(m) == BOUNDARY_CELL_TYPE  || discretization.cell_type_mesh(m) == DIRICHLET_CELL_TYPE )){
                    const T_INDEX& cell_index = discretization.cell_indices_mesh(m);
                    int i,j,ij;cell_index.Get(i,j,ij);    
                    int v=1;
                    int corners[8];
                    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));
                        iterator.Valid();iterator.Next(),v++){
                        if(discretization.cells_mesh(m)(v) > 0)
                            corners[v-1] = particle_mesh_hash.Get(discretization.cells_mesh(m)(v));
                        else
                            corners[v-1] = particle_hash.Get(iterator.Index());
                    }
                    if((i+j+ij)%2==0){
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[0],corners[4], corners[2],corners[1]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[4],corners[5], corners[7],corners[1]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[2],corners[6], corners[7],corners[4]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[3],corners[7], corners[1],corners[2]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[4],corners[1], corners[7],corners[2]));}
                    else{
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[0],corners[4], corners[6],corners[5]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[0],corners[2], corners[3],corners[6]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[3],corners[5], corners[1],corners[0]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[3],corners[7], corners[5],corners[6]));
                        subdomain_volume.mesh.elements.Append(VECTOR<int,4>(corners[3],corners[0], corners[6],corners[5]));}        
                }
            }
        }};

#endif

}

//#####################################################################
// Function Write_Output (3D)
//#####################################################################

template<typename T_ELASTICITY>
void Write_Output(STREAM_TYPE stream_type,CLElib& deformer,const std::string directory, const int frame, const int impl)
{
#if 0
	typedef float T;
    static const int d=3;

    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    typedef ELASTIC_LATTICE_DEFORMER T_IMPLEMENTATION;
    typedef CORE_ENGINE_WRAPPER::T_DISCRETIZATION T_DISCRETIZATION;
    T_IMPLEMENTATION& implementation=deformer.Implementation<T_IMPLEMENTATION>();
    CORE_ENGINE_WRAPPER& discretization=implementation.Discretization();
    FILE_UTILITIES::Create_Directory(directory+"/"+STRING_UTILITIES::Value_To_String(frame));

    deformer.Update_Embedded_Surfaces();


    Initialize_Geometry_Particle();Initialize_Read_Write_Structures();
    GEOMETRY_PARTICLES<TV> particles;
    ARRAY<STRUCTURE<TV>* > display_structures;
    ARRAY<TRIANGULATED_SURFACE<T>*> bone_surfaces;
    ARRAY<TRIANGULATED_SURFACE<T>*> embedded_surfaces;


    int node_grid_m=discretization.Unpadded_Node_Domain().Edge_Lengths()(1)+1;
    int node_grid_n=discretization.Unpadded_Node_Domain().Edge_Lengths()(2)+1;
    int node_grid_mn=discretization.Unpadded_Node_Domain().Edge_Lengths()(3)+1;

    // Initialize particles
    HASHTABLE<T_INDEX,int> particle_hash;
    HASHTABLE<int,int> particle_mesh_hash;

    HASHTABLE<int,T_INDEX> particle_hash_r;
    int last_grid_particle;
    int last_constraint_particle;
    int last_object_particle;

    particles.array_collection->Add_Elements(node_grid_m*node_grid_n*node_grid_mn);
    {int p=0;
    for(RANGE_ITERATOR<d> iterator(discretization.Unpadded_Node_Domain());iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        p++;
        particle_hash.Get_Or_Insert(index,p);
        particle_hash_r.Get_Or_Insert(p,index);
        if(discretization.node_is_active(index) || discretization.node_is_dirichlet(index)) for(int v=1;v<=d;v++) particles.X(p)(v)=implementation.u_internal.x(v)(index)+discretization.grid.Node(index)(v);
        else
            particles.X(p)=TV();}}

    last_grid_particle = particles.X.m;
    Fill_Mesh_Particles<T,d,T_ELASTICITY,CG_POLICY<T_ELASTICITY>::SUPPORTS_MESH>()(discretization,implementation.u_internal,particle_mesh_hash,particles);


    // Interior tetrahedralized volume
    TETRAHEDRALIZED_VOLUME<T>& interior_tetrahedralized_volume=*TETRAHEDRALIZED_VOLUME<T>::Create(particles);
    interior_tetrahedralized_volume.Update_Number_Nodes();
    for(RANGE_ITERATOR<d> iterator(discretization.unpadded_cell_domain);iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        if(discretization.cell_type(index)==INTERIOR_CELL_TYPE || discretization.cell_type(index)==BOUNDARY_CELL_TYPE){
            int i,j,ij;index.Get(i,j,ij);
            const int q=node_grid_n;
            const int pq=node_grid_mn;
            const int LDB=(i-1)*q*pq+(j-1)*pq+(ij-1)+1;
            const int LDF=(i-1)*q*pq+(j-1)*pq+(ij-0)+1;
            const int LUB=(i-1)*q*pq+(j-0)*pq+(ij-1)+1;
            const int LUF=(i-1)*q*pq+(j-0)*pq+(ij-0)+1;
            const int RDB=(i-0)*q*pq+(j-1)*pq+(ij-1)+1;
            const int RDF=(i-0)*q*pq+(j-1)*pq+(ij-0)+1;
            const int RUB=(i-0)*q*pq+(j-0)*pq+(ij-1)+1;
            const int RUF=(i-0)*q*pq+(j-0)*pq+(ij-0)+1;

            PHYSBAM_ASSERT( LDB <= last_grid_particle && LDF <= last_grid_particle && 
                            LUB <= last_grid_particle && LUF <= last_grid_particle && 
                            RDB <= last_grid_particle && RDF <= last_grid_particle && 
                            RUB <= last_grid_particle && RUF <= last_grid_particle ); 

            if((i+j+ij)%2==0){
                interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LDB,RDB,LUB,LDF));
                interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(RDB,RDF,RUF,LDF));
                interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUB,RUB,RUF,RDB));
                interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,RUF,LDF,LUB));
                interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(RDB,LDF,RUF,LUB));}
            else{
                interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LDB,RDB,RUB,RDF));
                interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LDB,LUB,LUF,RUB));
                interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,RDF,LDF,LDB));
                interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,RUF,RDF,RUB));
                interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,LDB,RUB,RDF));}}}
    display_structures.Append(&interior_tetrahedralized_volume);



    // Interior tetrahedralized volume mesh
    TETRAHEDRALIZED_VOLUME<T>& interior_mesh_tetrahedralized_volume=*TETRAHEDRALIZED_VOLUME<T>::Create(particles);
    interior_mesh_tetrahedralized_volume.Update_Number_Nodes();
    Create_Mesh_Tetrahedralized_Volume<T,d,T_DISCRETIZATION,CG_POLICY<T_DISCRETIZATION>::SUPPORTS_MESH>()(discretization,particle_hash,particle_mesh_hash,interior_mesh_tetrahedralized_volume);
    display_structures.Append(&interior_mesh_tetrahedralized_volume);

    // Interior tetrahedralized volume (assigned to odd blocks)
#if 0
    TETRAHEDRALIZED_VOLUME<T>& odd_partition_tetrahedralized_volume=*TETRAHEDRALIZED_VOLUME<T>::Create(particles);
    odd_partition_tetrahedralized_volume.Update_Number_Nodes();

    for(int partition=1;partition<=discretization.cell_block_partition_offsets.m;partition+=2){
        int block_begin=discretization.cell_block_partition_offsets(partition)+1;
        int block_end=(partition<discretization.cell_block_partition_offsets.m)
            ?discretization.cell_block_partition_offsets(partition+1)
            :discretization.cell_block_base_indices.m;
        for(int block=block_begin;block<=block_end;block++){
            const T_INDEX& base_index=discretization.cell_block_base_indices(block);

            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+1));iterator.Valid();iterator.Next()){
                const T_INDEX& index=iterator.Index();
                if(discretization.unpadded_cell_domain.Lazy_Outside(index)) continue;
                if(discretization.cell_type(index)==INTERIOR_CELL_TYPE || discretization.cell_type(index)==BOUNDARY_CELL_TYPE){
                    int i,j,ij;index.Get(i,j,ij);
                    const int q=node_grid_n;
                    const int pq=node_grid_mn;
                    const int LDB=(i-1)*q*pq+(j-1)*pq+(ij-1)+1;
                    const int LDF=(i-1)*q*pq+(j-1)*pq+(ij-0)+1;
                    const int LUB=(i-1)*q*pq+(j-0)*pq+(ij-1)+1;
                    const int LUF=(i-1)*q*pq+(j-0)*pq+(ij-0)+1;
                    const int RDB=(i-0)*q*pq+(j-1)*pq+(ij-1)+1;
                    const int RDF=(i-0)*q*pq+(j-1)*pq+(ij-0)+1;
                    const int RUB=(i-0)*q*pq+(j-0)*pq+(ij-1)+1;
                    const int RUF=(i-0)*q*pq+(j-0)*pq+(ij-0)+1;

                    if((i+j+ij)%2==0){
                        odd_partition_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LDB,RDB,LUB,LDF));
                        odd_partition_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(RDB,RDF,RUF,LDF));
                        odd_partition_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUB,RUB,RUF,RDB));
                        odd_partition_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,RUF,LDF,LUB));
                        odd_partition_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(RDB,LDF,RUF,LUB));}
                    else{
                        odd_partition_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LDB,RDB,RUB,RDF));
                        odd_partition_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LDB,LUB,LUF,RUB));
                        odd_partition_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,RDF,LDF,LDB));
                        odd_partition_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,RUF,RDF,RUB));
                        odd_partition_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,LDB,RUB,RDF));}}}}}
    display_structures.Append(&odd_partition_tetrahedralized_volume);
#endif

    // Interior segmented curve
    interior_tetrahedralized_volume.mesh.Initialize_Segment_Mesh();
    SEGMENTED_CURVE<TV>& interior_segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);
    interior_segmented_curve.mesh.Initialize_Mesh(*interior_tetrahedralized_volume.mesh.segment_mesh);
    for(int s=interior_segmented_curve.mesh.elements.m;s>=1;s--){
        VECTOR<int,2>& segment=interior_segmented_curve.mesh.elements(s);
        int distance=(segment.x>segment.y)?segment.x-segment.y:segment.y-segment.x;
        if(distance!=1 && distance!=node_grid_mn && distance!=node_grid_n*node_grid_mn) interior_segmented_curve.mesh.elements.Remove_Index_Lazy(s);}
    display_structures.Append(&interior_segmented_curve);

    // Dirichlet tetrahedralized volume
    TETRAHEDRALIZED_VOLUME<T>& dirichlet_tetrahedralized_volume=*TETRAHEDRALIZED_VOLUME<T>::Create(particles);
    dirichlet_tetrahedralized_volume.Update_Number_Nodes();
    for(RANGE_ITERATOR<d> iterator(discretization.unpadded_cell_domain);iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        if(discretization.cell_type(index)==DIRICHLET_CELL_TYPE){
            int i,j,ij;index.Get(i,j,ij);
            const int q=node_grid_n;
            const int pq=node_grid_mn;
            const int LDB=(i-1)*q*pq+(j-1)*pq+(ij-1)+1;
            const int LDF=(i-1)*q*pq+(j-1)*pq+(ij-0)+1;
            const int LUB=(i-1)*q*pq+(j-0)*pq+(ij-1)+1;
            const int LUF=(i-1)*q*pq+(j-0)*pq+(ij-0)+1;
            const int RDB=(i-0)*q*pq+(j-1)*pq+(ij-1)+1;
            const int RDF=(i-0)*q*pq+(j-1)*pq+(ij-0)+1;
            const int RUB=(i-0)*q*pq+(j-0)*pq+(ij-1)+1;
            const int RUF=(i-0)*q*pq+(j-0)*pq+(ij-0)+1;

            if((i+j+ij)%2==0){
                dirichlet_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LDB,RDB,LUB,LDF));
                dirichlet_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(RDB,RDF,RUF,LDF));
                dirichlet_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUB,RUB,RUF,RDB));
                dirichlet_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,RUF,LDF,LUB));
                dirichlet_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(RDB,LDF,RUF,LUB));}
            else{
                dirichlet_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LDB,RDB,RUB,RDF));
                dirichlet_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LDB,LUB,LUF,RUB));
                dirichlet_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,RDF,LDF,LDB));
                dirichlet_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,RUF,RDF,RUB));
                dirichlet_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,LDB,RUB,RDF));}}}
    display_structures.Append(&dirichlet_tetrahedralized_volume);

    // Dirichlet segmented curve
    dirichlet_tetrahedralized_volume.mesh.Initialize_Segment_Mesh();
    SEGMENTED_CURVE<TV>& dirichlet_segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);
    dirichlet_segmented_curve.mesh.Initialize_Mesh(*dirichlet_tetrahedralized_volume.mesh.segment_mesh);
    for(int s=dirichlet_segmented_curve.mesh.elements.m;s>=1;s--){
        VECTOR<int,2>& segment=dirichlet_segmented_curve.mesh.elements(s);
        int distance=(segment.x>segment.y)?segment.x-segment.y:segment.y-segment.x;
        if(distance!=1 && distance!=node_grid_mn && distance!=node_grid_n*node_grid_mn) dirichlet_segmented_curve.mesh.elements.Remove_Index_Lazy(s);}
    display_structures.Append(&dirichlet_segmented_curve);

    // Active variables
    FREE_PARTICLES<TV>& active_variables=*FREE_PARTICLES<TV>::Create(particles);
    active_variables.Update_Number_Nodes();
    for(RANGE_ITERATOR<d> iterator(discretization.Unpadded_Node_Domain());iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        if(discretization.node_is_active(index)){
            const int q=node_grid_n;
            const int pq=node_grid_mn;
            active_variables.nodes.Append((index.x-1)*q*pq+(index.y-1)*pq+index.z);}}
    display_structures.Append(&active_variables);

    // Dirichlet variables
    FREE_PARTICLES<TV>& dirichlet_variables=*FREE_PARTICLES<TV>::Create(particles);
    dirichlet_variables.Update_Number_Nodes();
    for(RANGE_ITERATOR<d> iterator(discretization.Unpadded_Node_Domain());iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        if(discretization.node_is_dirichlet(index)){
            const int& q=node_grid_n;
            const int& pq=node_grid_mn;
            dirichlet_variables.nodes.Append((index.x-1)*q*pq+(index.y-1)*pq+index.z);}}
    display_structures.Append(&dirichlet_variables);

    // Mesh-only nodes
    FREE_PARTICLES<TV>& mesh_only_variables=*FREE_PARTICLES<TV>::Create(particles);
    mesh_only_variables.Update_Number_Nodes();
    Create_Mesh_Only_Variables<T,d,T_DISCRETIZATION,CG_POLICY<T_DISCRETIZATION>::SUPPORTS_MESH>()(discretization,particle_hash,particle_mesh_hash,mesh_only_variables);
    display_structures.Append(&mesh_only_variables);

    // Mesh-only nodes
    FREE_PARTICLES<TV>& mesh_shared_variables=*FREE_PARTICLES<TV>::Create(particles);
    mesh_shared_variables.Update_Number_Nodes();
    Create_Mesh_Shared_Variables<T,d,T_DISCRETIZATION,CG_POLICY<T_DISCRETIZATION>::SUPPORTS_MESH>()(discretization,particle_hash,particle_mesh_hash,mesh_shared_variables);
    display_structures.Append(&mesh_shared_variables);

/*
    // Active variables
    FREE_PARTICLES<TV>& mesh_nodes_on_grid=*FREE_PARTICLES<TV>::Create(particles);
    mesh_nodes_on_grid.Update_Number_Nodes();
    for(int m=1;m<=discretization.Number_Of_Mesh_Cells();m++){
        const T_INDEX& cell_index = discretization.cell_indices_mesh(m);
        int i,j,ij;cell_index.Get(i,j,ij);    
        int v=1;
        int corners[8];
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));
            iterator.Valid();iterator.Next(),v++){  
            if( discretization.cells_mesh(m)(v) == 0 )
                mesh_nodes_on_grid.nodes.Append(particle_hash.Get(iterator.Index()));
        }
    }
    display_structures.Append(&mesh_nodes_on_grid);

    // Active variables
    FREE_PARTICLES<TV>& mesh_nodes_on_mesh=*FREE_PARTICLES<TV>::Create(particles);
    mesh_nodes_on_mesh.Update_Number_Nodes();
    for(int m=1;m<=discretization.Number_Of_Mesh_Cells();m++){
        const T_INDEX& cell_index = discretization.cell_indices_mesh(m);
        int i,j,ij;cell_index.Get(i,j,ij);    
        int v=1;
        int corners[8];
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));
            iterator.Valid();iterator.Next(),v++){  
            if( discretization.cells_mesh(m)(v) > 0 )
                mesh_nodes_on_mesh.nodes.Append(particle_hash.Get(iterator.Index()));
        }
    }
    display_structures.Append(&mesh_nodes_on_mesh);
*/

    // Connection points

    // Connection points
   
    FREE_PARTICLES<TV>& dynamic_constraint_endpoints=*FREE_PARTICLES<TV>::Create();
    SEGMENTED_CURVE<TV>& dynamic_constraint_segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);
    for(int c=1;c<=discretization.dynamic_point_constraints.m;c++){
        int p[2];
        for( int e=0;e<2;e++){
            const CONSTRAINT_NODE<T,d>& endpoint = discretization.dynamic_point_constraints(c).endpoints[e];
            switch(endpoint.type){
            case CONSTRAINT_NODE<T,d>::GRID_FIXED:
                {
                    p[e]=particles.array_collection->Add_Element();
                    particles.X(p[e])=discretization.Deformation(endpoint.grid_index(),
                                                                 endpoint.multilinear_coordinates(),
                                                                 implementation.u_internal);
                    dynamic_constraint_endpoints.nodes.Append(p[e]);
                }
                break;
            case CONSTRAINT_NODE<T,d>::MESH_FIXED:
                {
                    p[e]=Create_Mesh_Constraint<T,d,T_DISCRETIZATION,CG_POLICY<T_DISCRETIZATION>::SUPPORTS_MESH>()(discretization,implementation.u_internal,particles,endpoint,dynamic_constraint_endpoints);
                }
                break;
            case CONSTRAINT_NODE<T,d>::KINEMATIC:
                {
                    p[e]=particles.array_collection->Add_Element();
                    particles.X(p[e])=endpoint.spatial_location();
                    dynamic_constraint_endpoints.nodes.Append(p[e]);
                }
                break;
            }
        }
        dynamic_constraint_segmented_curve.mesh.elements.Append(VECTOR<int,2>(p[0],p[1]));      
    }
    
    display_structures.Append(&dynamic_constraint_endpoints);
    display_structures.Append(&dynamic_constraint_segmented_curve);

    FREE_PARTICLES<TV>& static_constraint_endpoints=*FREE_PARTICLES<TV>::Create();
    SEGMENTED_CURVE<TV>& static_constraint_segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);
    for(int c=1;c<=discretization.static_point_constraints.m;c++){
        int p[2];
        for( int e=0;e<2;e++){
            const CONSTRAINT_NODE<T,d>& endpoint = discretization.static_point_constraints(c).endpoints[e];
            switch(endpoint.type){
            case CONSTRAINT_NODE<T,d>::GRID_FIXED:
                {
                    p[e]=particles.array_collection->Add_Element();
                    particles.X(p[e])=discretization.Deformation(endpoint.grid_index(),
                                                                 endpoint.multilinear_coordinates(), implementation.u_internal);
                    static_constraint_endpoints.nodes.Append(p[e]);
                }
                break;
            
            case CONSTRAINT_NODE<T,d>::MESH_FIXED:
                {
                    p[e]=Create_Mesh_Constraint<T,d,T_DISCRETIZATION,CG_POLICY<T_DISCRETIZATION>::SUPPORTS_MESH>()(discretization,implementation.u_internal,particles,endpoint,static_constraint_endpoints);
                }
                break;
            
            case CONSTRAINT_NODE<T,d>::KINEMATIC:
                {
                    p[e]=particles.array_collection->Add_Element();
                    particles.X(p[e])=endpoint.spatial_location();
                    static_constraint_endpoints.nodes.Append(p[e]);
                }
                break;
            }
        }
        static_constraint_segmented_curve.mesh.elements.Append(VECTOR<int,2>(p[0],p[1]));      
    }
    
    display_structures.Append(&static_constraint_endpoints);
    display_structures.Append(&static_constraint_segmented_curve);
    last_constraint_particle = particles.X.m;

    // Add bone surfaces

    // for(int b=1;b<=implementation.Number_Of_Bones();b++){
    //      const TRIANGULATED_SURFACE<T>& bone_surface=implementation.Bone_Surface(b);
    //      TRIANGULATED_SURFACE<T>& display_bone_surface=dynamic_cast<TRIANGULATED_SURFACE<T>&>(*bone_surface.Append_Particles_And_Create_Copy(particles));
    //      	display_structures.Append(&display_bone_surface);
    //      	bone_surfaces.Append(&display_bone_surface);
    //      }
#if 1
    // Add embedded surfaces
    {
        TRIANGULATED_SURFACE<T>* display_embedded_surface=TRIANGULATED_SURFACE<T>::Create(particles);
        
        int particle_base = particles.X.m;
        


        pthread_mutex_lock(&implementation.geometry_buffer_lock);

        const ARRAY<TV>* vertices=implementation.display_vertices;

        for( int v=1;v<=vertices->m;v++){
            int p = particles.array_collection->Add_Element();
            for(int w=1;w<=d;w++)
                particles.X(p)(w)=vertices->operator()(v)(w);
        }
        for( int t=1;t<=implementation.triangles.m;t++){
            display_embedded_surface->mesh.elements.Append(implementation.triangles(t)+particle_base);
        }
        
        pthread_mutex_unlock(&implementation.geometry_buffer_lock);



        display_structures.Append(display_embedded_surface);
        embedded_surfaces.Append(display_embedded_surface);
        
        // Add static models
        
        LOG::cout << "Adding " << implementation.static_model_vertices.m << " static models." << std::endl;
        
        for( int M = 1; M <= implementation.static_model_vertices.m; M++)
            {
                LOG::cout << "Static Model (" << M << ") has " << implementation.static_model_vertices(M).m << " vertices and " << implementation.static_model_triangles(M).m << " triangles." << std::endl;
                
                TRIANGULATED_SURFACE<T>* display_embedded_surface=TRIANGULATED_SURFACE<T>::Create(particles);
                
                int particle_base = particles.X.m;
                
                for( int v=1;v<=implementation.static_model_vertices(M).m;v++){
                    int p = particles.array_collection->Add_Element();
                particles.X(p)=implementation.static_model_vertices(M)(v);
                //LOG::cout << particles.X(p) << std::endl;
                
                }
                for( int t=1;t<=implementation.static_model_triangles(M).m;t++){
                    display_embedded_surface->mesh.elements.Append(implementation.static_model_triangles(M)(t)+particle_base);
                }
                
                display_structures.Append(display_embedded_surface);
                embedded_surfaces.Append(display_embedded_surface);
            }
    }


//    for(int s=1;s<=implementation.Number_Of_Embedded_Surfaces();s++){
//       const TRIANGULATED_SURFACE<T>& embedded_surface=implementation.Embedded_Surface(s);
//        TRIANGULATED_SURFACE<T>& display_embedded_surface=dynamic_cast<TRIANGULATED_SURFACE<T>&>(*embedded_surface.Append_Particles_And_Create_Copy(particles));
//        display_structures.Append(&display_embedded_surface);
//        embedded_surfaces.Append(&display_embedded_surface);}
    last_object_particle = particles.X.m;
#endif    

#ifdef DISPLAY_FINE_GRID
    // Add Cutting fine grid mesh
   
    // Initialize particles
    int fine_grid_m=implementation.fine_grid.Node_Indices().Edge_Lengths()(1)+1;
    int fine_grid_n=implementation.fine_grid.Node_Indices().Edge_Lengths()(2)+1;
    int fine_grid_mn=implementation.fine_grid.Node_Indices().Edge_Lengths()(3)+1;

    int particle_base = particles.X.m;
    particles.array_collection->Add_Elements(fine_grid_m*fine_grid_n*fine_grid_mn);
    {int p=particle_base;
    for(RANGE_ITERATOR<d> iterator(implementation.unpadded_fine_node_domain);iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        p++;
        //if(implementation.voxmap_node(index)==1)
        for(int v=1;v<=d;v++)
            particles.X(p)(v)=implementation.u_fine(v)(index)+implementation.fine_grid.Node(index)(v);
            //else
            //particles.X(p)=TV();
    }}

    // Interior tetrahedralized volume
    TETRAHEDRALIZED_VOLUME<T>& fine_tetrahedralized_volume=*TETRAHEDRALIZED_VOLUME<T>::Create(particles);
    fine_tetrahedralized_volume.Update_Number_Nodes();
    const ARRAY<T_INDEX>& embedding_map=implementation.embedding_map;

    //for( int v=1;v<=implementation.vertices.m;v++){
    //const T_INDEX& index = embedding_map(v);
    for(RANGE_ITERATOR<d> iterator(implementation.fine_grid.Cell_Indices());iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
    
        if(implementation.voxmap_dirichlet(index)==1 && implementation.voxmap( index )){
            //if(){
            int i,j,ij;index.Get(i,j,ij);
            const int q=fine_grid_n;
            const int pq=fine_grid_mn;
            const int LDB=(i-1)*q*pq+(j-1)*pq+(ij-1)+1+particle_base;
            const int LDF=(i-1)*q*pq+(j-1)*pq+(ij-0)+1+particle_base;
            const int LUB=(i-1)*q*pq+(j-0)*pq+(ij-1)+1+particle_base;
            const int LUF=(i-1)*q*pq+(j-0)*pq+(ij-0)+1+particle_base;
            const int RDB=(i-0)*q*pq+(j-1)*pq+(ij-1)+1+particle_base;
            const int RDF=(i-0)*q*pq+(j-1)*pq+(ij-0)+1+particle_base;
            const int RUB=(i-0)*q*pq+(j-0)*pq+(ij-1)+1+particle_base;
            const int RUF=(i-0)*q*pq+(j-0)*pq+(ij-0)+1+particle_base;

            if((i+j+ij)%2==0){
                fine_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LDB,RDB,LUB,LDF));
                fine_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(RDB,RDF,RUF,LDF));
                fine_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUB,RUB,RUF,RDB));
                fine_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,RUF,LDF,LUB));
                fine_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(RDB,LDF,RUF,LUB));}
            else{
                fine_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LDB,RDB,RUB,RDF));
                fine_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LDB,LUB,LUF,RUB));
                fine_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,RDF,LDF,LDB));
                fine_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,RUF,RDF,RUB));
                fine_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,LDB,RUB,RDF));}}}
    display_structures.Append(&fine_tetrahedralized_volume);
#endif



    // Update particle count on all structures
    for(int s=1;s<=display_structures.m;s++)
        display_structures(s)->Update_Number_Nodes();

#if 0
#if 0
    ARRAY<TV> X_debug(particles.X);
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/X_debug"+STRING_UTILITIES::Value_To_String(frame),X_debug);
#else
    ARRAY<TV> X_debug;
    FILE_UTILITIES::Read_From_File(stream_type,directory+"/X_debug"+STRING_UTILITIES::Value_To_String(frame),X_debug);

    LOG::cout << "    last_grid_particle  " <<  last_grid_particle << std::endl;
    LOG::cout << "    last_constraint_particle  " <<  last_constraint_particle << std::endl;
    LOG::cout << "    last_object_particle  " <<  last_object_particle << std::endl;

    for(int i=1; i<= X_debug.m; i++)
        if( X_debug(i) != particles.X(i) ){
            LOG::cout << "i : " << i << std::endl;
            LOG::cout << "X_debug: " << X_debug(i) << std::endl;
            LOG::cout << "particles: " << particles.X(i) << std::endl;
            LOG::cout << "index: " << particle_hash_r.Get_Or_Insert(i,T_INDEX(-1,-1,-1)) << std::endl;

            for(int m = 1; m <= discretization.Number_Of_Mesh_Cells(); m++)
                if( discretization.cell_indices_mesh(m) == particle_hash_r.Get_Or_Insert(i,T_INDEX(-1,-1,-1)) ){
                    LOG::cout << "This is a mesh cell!  ID: " << m << std::endl;
                    
                    
                }

            exit(1);
        }
#endif
#endif

    DEFORMABLE_GEOMETRY_COLLECTION<TV> collection(particles);

    // Select which structures to display
    
    if(interior_tetrahedralized_volume.mesh.elements.m == 0 )
        delete &interior_tetrahedralized_volume;
    else{
        collection.Add_Structure(&interior_tetrahedralized_volume);
        //delete &interior_tetrahedralized_volume;
    }
    
    if(interior_mesh_tetrahedralized_volume.mesh.elements.m == 0)
        delete &interior_mesh_tetrahedralized_volume;
    else{
        collection.Add_Structure(&interior_mesh_tetrahedralized_volume);
        //delete &interior_mesh_tetrahedralized_volume;
    }
#if 0
    for( int i = 1; i <= discretization.number_of_subdomains; i++){
        if(subdomain_volumes(i)->mesh.elements.m == 0){
            delete subdomain_volumes(i);
            subdomain_volumes(i) = 0;            
        }
        else{
            collection.Add_Structure(subdomain_volumes(i));
            //delete subdomain_volumes(i);
            //subdomain_volumes(i) = 0;   
        }
    }
#endif
    
#if 0
    if(odd_partition_tetrahedralized_volume.mesh.elements.m == 0)
        delete &odd_partition_tetrahedralized_volume;
    else{
        //collection.Add_Structure(&odd_partition_tetrahedralized_volume);
        delete &odd_partition_tetrahedralized_volume;    
    }
#endif

#ifdef DISPLAY_FINE_GRID
    collection.Add_Structure(&fine_tetrahedralized_volume);
		//delete &interior_tetrahedralized_volume;
#endif

    //collection.Add_Structure(&interior_segmented_curve);
    delete &interior_segmented_curve;

    if(dirichlet_tetrahedralized_volume.mesh.elements.m == 0 )
        delete &dirichlet_tetrahedralized_volume;
    else{
        collection.Add_Structure(&dirichlet_tetrahedralized_volume);
        //delete &dirichlet_tetrahedralized_volume;
    }

    // collection.Add_Structure(&dirichlet_segmented_curve);
    delete &dirichlet_segmented_curve;

    // collection.Add_Structure(&active_variables);
    delete &active_variables;

    //collection.Add_Structure(&mesh_only_variables);
    delete &mesh_only_variables;

    //collection.Add_Structure(&mesh_shared_variables);
    delete &mesh_shared_variables;

#if 0
    collection.Add_Structure(&boundary_variables);
    //delete &boundary_variables;
#endif
    // collection.Add_Structure(&dirichlet_variables);
    // delete &dirichlet_variables;

    //collection.Add_Structure(&mesh_nodes_on_grid);
    //delete &mesh_nodes_on_grid;

    //collection.Add_Structure(&mesh_nodes_on_mesh);
    //delete &mesh_nodes_on_mesh;

    //collection.Add_Structure(&dynamic_constraint_endpoints);
    delete &dynamic_constraint_endpoints;

    //collection.Add_Structure(&dynamic_constraint_segmented_curve);
    delete &dynamic_constraint_segmented_curve;

    //collection.Add_Structure(&static_constraint_endpoints);
    delete &static_constraint_endpoints;

    //collection.Add_Structure(&static_constraint_segmented_curve);
    delete &static_constraint_segmented_curve;

#if 0
    for(int b=1;b<=bone_surfaces.m;b++){
        //collection.Add_Structure(bone_surfaces(b));
        		delete bone_surfaces(b);
    }
#endif
    for(int s=1;s<=embedded_surfaces.m;s++){
        collection.Add_Structure(embedded_surfaces(s));
        //delete embedded_surfaces(s);
    }


    // Write PhysBAM-style output
    collection.Write(stream_type,directory,frame,frame,true);
	//FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+STRING_UTILITIES::Value_To_String(frame)+"/u_values",discretization.u);
	
#endif
}
//#####################################################################
//	template<float> void Write_Output (STREAM_TYPE stream_type,EMBEDDED_DEFORMER& deformer,const std::string directory,const int frame);
//	template<double> void Write_Output(STREAM_TYPE stream_type,EMBEDDED_DEFORMER& deformer,const std::string directory,const int frame);

    template
    void Write_Output<HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,3,true,false> > >(STREAM_TYPE stream_type,CLElib& deformer,const std::string directory, const int frame, const int impl);
    
#if 1
    template
    void Write_Output<SKINNING_NONLINEAR_ELASTICITY<float,3,true,false> >(STREAM_TYPE stream_type,CLElib& deformer,const std::string directory, const int frame, const int impl);
#endif

}

#endif
#endif
