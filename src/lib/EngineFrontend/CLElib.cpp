//#####################################################################
// Copyright 2010-2012, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CLElib
//#####################################################################
#include "CLElib.h"

#include <fstream>
#include <pthread.h>
#include <errno.h>

#include <PhysBAM_Tools/Data_Structures/QUEUE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>

#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>

#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>

#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGULATED_SURFACE_SUBDIVISION.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>

#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include "Write_Output.h"

#include <EngineInterface/CONSTRAINTS.h>

#include <Thread_Queueing/PTHREAD_QUEUE.h>

#include "ELASTIC_LATTICE_DEFORMER.h"
#include "EMBEDDING_TOOLS.h"
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include "PHYSBAM_LEVELSET_COLLISION.h"
#include <Common/STENCIL.h>

#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>

// New Cutter
#include <Common_Geometry/Topology/REGULAR_HYPERCUBE_MESH.h>
#include <Common_Geometry/Nonmanifold_Topology_Generation/CUTTER.h>
#include <Common_Geometry/Nonmanifold_Topology_Generation/CUTTER_NONMANIFOLD_LEVELSET_STRATEGY.h>
#include <Common_Geometry/Nonmanifold_Topology_Generation/CUTTER_BASIC_STRATEGY.h>
#include <Common_Geometry/Nonmanifold_Topology_Generation/MATERIAL_PREDICATE_VOXELIZED_VOLUME.h>

#ifdef ENABLE_SELF_COLLISIONS
//TetGen
#define TETLIBRARY
#include <tetgen.h>
#endif

using namespace PhysBAM;

#define GRID_IN_GRID
// Set to 1<<30 for full logging.
#define LOG_DEPTH 1<<30
#define LOG_FILE "CLElib.log"

template<class T>
void Write_To_File(const std::string& filename,const std::vector<T>& v)
{
    std::ofstream output(filename.c_str(),std::ios::binary);
    const int size=v.size();
    output.write(reinterpret_cast<const char*>(&size),sizeof(int));
    output.write(reinterpret_cast<const char*>(&v[0]),v.size()*sizeof(T));
    output.close();
}

template<class T>
void Read_From_File(const std::string& filename,std::vector<T>& v)
{
    v.clear();
    std::ifstream input(filename.c_str(),std::ios::binary);
    int size;
    input.read(reinterpret_cast<char*>(&size),sizeof(int));
    v.resize(size);
    input.read(reinterpret_cast<char*>(&v[0]),v.size()*sizeof(T));
    input.close();
}

// Hack!
template<class T, int d>
struct NONLINEAR_ELASTICITY{
    static int Get_Lower_Cell_Padding() {return 1;};
    static int Get_Upper_Cell_Padding() {return 3;};
    static int Get_Lower_Node_Padding() {return 1;};
    static int Get_Upper_Node_Padding() {return 2;};
};
// Hack!

//#define ENABLE_LOG_MESSAGES
using namespace PhysBAM;
//#define UNIT_ACTIVATION_TEST
//#define ZERO_FRAME_TEST
//#####################################################################
// Constructor
//#####################################################################
#define PRIMARY_ELASTICITY ENGINE_INTERFACE
#define PRIMARY_DEFORMER ELASTIC_LATTICE_DEFORMER


CLElib::
CLElib()
{
#ifdef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::CLElib()");
#endif

    
    implementation=new PRIMARY_DEFORMER;
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation< PRIMARY_DEFORMER >();
//    elastic_lattice_deformer.parameter_animation=0;
    //elastic_lattice_deformer.discretization=0;
    LOG::Initialize_Logging(false, false, LOG_DEPTH , false);

    Initialize_Geometry_Particle();Initialize_Read_Write_Structures();

    refinement = 10;

    frame_counter=-1;
    debug_output = 0;

    Set_Youngs_Modulus(3e4);
    Set_Poissons_Ratio(.3);
    Set_Density(1e-6); // Water, in mm-kg-sec units
    Set_Damping(1e-5);
    Set_Time_Step(FLT_MAX);

    std::ofstream output(LOG_FILE);
    output << "Begin." << std::endl;
    output.close();
}
//#####################################################################
// Destructor
//#####################################################################
CLElib::
~CLElib()
{
#ifdef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::~CLElib()");
#endif

    Destroy_Model();
    delete implementation;
    implementation=0;
    
    LOG::Finish_Logging();
}
//#####################################################################
// Function Set_Youngs_Modulus
//#####################################################################
void CLElib::
Set_Youngs_Modulus(const float youngs_modulus)
{
#ifdef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Set_Youngs_Modulus()");
#endif

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Youngs modulus may only be adjusted prior to creating the model");
    elastic_lattice_deformer.youngs_modulus=youngs_modulus;
}
//#####################################################################
// Function Set_Poissons_Ratio
//#####################################################################
void CLElib::
Set_Poissons_Ratio(const float poissons_ratio)
{
#ifdef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Set_Poissons_Ratio()");
#endif

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Poisson's ratio may only be adjusted prior to creating the model");
    elastic_lattice_deformer.poissons_ratio=poissons_ratio;
}
//#####################################################################
// Function Create_Model
//#####################################################################
void CLElib::
//Create_Model(const std::string& cell_info_filename,const std::string& object_name, const std::string& deformation_filename)
Create_Model(const std::vector<float>& vertices_input,const std::vector<int>& triangles_input,const float dx, const int ref)
{
#ifdef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Create_Model()");
#endif

    refinement = ref;
    STREAM_TYPE stream_type((float)0);

    Write_To_File("vertices.dat",vertices_input);
    Write_To_File("triangles.dat",triangles_input);

    std::ofstream output(LOG_FILE, std::ios::app);
    output << "Create  " << dx << std::endl;
    output.close();

    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    Initialize_Geometry_Particle();

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Create_Model() called while model is already created");

    // Translate from std::vector to PhysBAM arrays (converting 0-base to 1-base)

    PHYSBAM_ASSERT(vertices_input.size()%3==0);
    // Allocate Vertex buffers
    ARRAY<TV>& vertices=elastic_lattice_deformer.vertices;
    vertices.Resize(vertices_input.size()/3);
    ARRAY<TV>& vertices_B=elastic_lattice_deformer.vertices_B;
    vertices_B.Resize(vertices_input.size()/3);
    // Allocate Stress Buffers
    ARRAY<TV>& stress=elastic_lattice_deformer.stress;
    stress.Resize(vertices_input.size()/3);
    ARRAY<TV>& stress_B=elastic_lattice_deformer.stress_B;
    stress_B.Resize(vertices_input.size()/3);
    // Allocate Strain Buffers
    ARRAY<TV>& strain=elastic_lattice_deformer.strain;
    strain.Resize(vertices_input.size()/3);
    ARRAY<TV>& strain_B=elastic_lattice_deformer.strain_B;
    strain_B.Resize(vertices_input.size()/3);


    for(int i=0,p=1;p<=vertices.m;p++)
        for(int v=1;v<=3;v++){
            vertices_B(p)(v)=vertices_input[i];
            vertices(p)(v)=vertices_input[i];
            i++;
        }

    PHYSBAM_ASSERT(triangles_input.size()%3==0);
    ARRAY<VECTOR<int,3> >& triangles=elastic_lattice_deformer.triangles;
    triangles.Resize(triangles_input.size()/3);
    for(int i=0,t=1;t<=triangles.m;t++)
        for(int v=1;v<=3;v++)
            triangles(t)(v)=triangles_input[i++]+1;

    // Sanity checks

    PHYSBAM_ASSERT(ARRAYS_COMPUTATIONS::Max(triangles.Flattened())<=vertices.m);
    PHYSBAM_ASSERT(ARRAYS_COMPUTATIONS::Min(triangles.Flattened())>0);
    for(int t=1;t<=triangles.m;t++){
        const VECTOR<int,3>& tri=triangles(t);
        PHYSBAM_ASSERT(tri(1)!=tri(2) && tri(2)!=tri(3) && tri(3)!=tri(1));}

    {
        TRIANGLE_MESH triangle_mesh(vertices.m,triangles);
        triangle_mesh.Initialize_Node_On_Boundary();
        PHYSBAM_ASSERT(!triangle_mesh.node_on_boundary->Contains(true));
    }

#ifdef ENABLE_SELF_COLLISIONS
    // Build Tetrahedralized Object
    // (w/ TetGen)
    GEOMETRY_PARTICLES<TV> tet_particles;
    TETRAHEDRALIZED_VOLUME<T>* tetvolume = TETRAHEDRALIZED_VOLUME<T>::Create(tet_particles);

    {
        LOG::SCOPE scope("Running TetGen");

        GEOMETRY_PARTICLES<TV> triangle_particles;
        TRIANGULATED_SURFACE<T>* triangle_surface=TRIANGULATED_SURFACE<T>::Create(triangle_particles);
        for( int v=1;v<=vertices.m;v++){
            int p = triangle_particles.array_collection->Add_Element();
            for(int w=1;w<=d;w++)
                triangle_particles.X(p)(w)=vertices(v)(w);
        }        
        for( int t=1;t<=triangles.m;t++){
            triangle_surface->mesh.elements.Append(triangles(t));
        }
        triangle_surface->Update_Number_Nodes();
        triangle_surface->Close_Surface(true, 1e-6, false, true );
        triangle_surface->Update_Number_Nodes();


        tetgenio tetgen_input, tetgen_output;
        tetgen_input.numberofpoints = (triangle_particles.X.Size());
        tetgen_input.pointlist = new REAL[ triangle_particles.X.Size() * 3 ];
        for( int i = 0; i<triangle_particles.X.Size(); i++ )
            for( int w = 0; w<d; w++)
                tetgen_input.pointlist[i*d+w] = triangle_particles.X(i+1)(w+1);
      	
        tetgen_input.numberoffacets = triangle_surface->mesh.elements.Size();
        tetgen_input.facetlist = new tetgenio::facet[tetgen_input.numberoffacets];
        tetgen_input.facetmarkerlist = new int[tetgen_input.numberoffacets];
        for( int t = 0; t < triangle_surface->mesh.elements.Size(); t++){
            auto f = &tetgen_input.facetlist[t];
            f->numberofpolygons = 1;
            f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
            f->numberofholes = 0;
            f->holelist = NULL;
            auto p = &f->polygonlist[0];
            p->numberofvertices = 3;
            p->vertexlist = new int[p->numberofvertices];
            p->vertexlist[0] = triangle_surface->mesh.elements(t+1)(1)-1;
            p->vertexlist[1] = triangle_surface->mesh.elements(t+1)(2)-1;
            p->vertexlist[2] = triangle_surface->mesh.elements(t+1)(3)-1;
            tetgen_input.facetmarkerlist[t] = 0;
        }
        
        tetgenbehavior switches;
        switches.plc = 1;
        switches.docheck = 1;
        switches.verbose = 1;
        
        try{
            tetrahedralize( "pq1.41a0.1dVc", &tetgen_input, &tetgen_output );
        }
        catch (int e ){
            LOG::cout << "Tetgen failed with error code : " << e << std::endl;
            LOG::cout.flush();
            PHYSBAM_FATAL_ERROR();
        }
        
        LOG::cout << tetgen_output.numberofpoints << "   " << vertices_input.size()/3 << std::endl;
        
        for( int v=0; v < tetgen_output.numberofpoints; v++ ){
            int p = tet_particles.array_collection->Add_Element();
            PHYSBAM_ASSERT( p == v+1 );
            for( int w=0;w<d;w++)
                tet_particles.X(p)(w+1)=tetgen_output.pointlist[ v*3 + w ];
        }
        for( int t=0; t<tetgen_output.numberoftetrahedra; t++){
            tetvolume->mesh.elements.Append( VECTOR<int,4>(tetgen_output.tetrahedronlist[t*4+0]+1,
                                                           tetgen_output.tetrahedronlist[t*4+1]+1,
                                                           tetgen_output.tetrahedronlist[t*4+2]+1,
                                                           tetgen_output.tetrahedronlist[t*4+3]+1 ));
            PHYSBAM_ASSERT( tetgen_output.tetrahedronlist[t*4+0]+1 >= 1 && 
                            tetgen_output.tetrahedronlist[t*4+0]+1 <= tetgen_output.numberofpoints );
            PHYSBAM_ASSERT( tetgen_output.tetrahedronlist[t*4+1]+1 >= 1 && 
                            tetgen_output.tetrahedronlist[t*4+1]+1 <= tetgen_output.numberofpoints );
            PHYSBAM_ASSERT( tetgen_output.tetrahedronlist[t*4+2]+1 >= 1 && 
                            tetgen_output.tetrahedronlist[t*4+2]+1 <= tetgen_output.numberofpoints );
            PHYSBAM_ASSERT( tetgen_output.tetrahedronlist[t*4+3]+1 >= 1 && 
                            tetgen_output.tetrahedronlist[t*4+3]+1 <= tetgen_output.numberofpoints );            
        }

        
        tetvolume->Update_Number_Nodes();
        
      {
          FILE_UTILITIES::Create_Directory("tetrahedralized_object/0");
          DEFORMABLE_GEOMETRY_COLLECTION<TV> collection(tet_particles);
          if(tetgen_output.numberoftetrahedra != 0 )
              collection.Add_Structure(tetvolume);
          collection.Add_Structure(triangle_surface);
          collection.Write(stream_type,"tetrahedralized_object",0,0,true);
      }

        if( tetgen_output.numberoftetrahedra == 0 )
            PHYSBAM_FATAL_ERROR( "Tetgen failed to mesh the model.");
    }
#endif
    
    // Build High Res Render Surface
    {
        GEOMETRY_PARTICLES<TV> particles;
        TRIANGULATED_SURFACE<T>* render_surface=TRIANGULATED_SURFACE<T>::Create(particles);
        for( int v=1;v<=vertices.m;v++){
            int p = particles.array_collection->Add_Element();
            for(int w=1;w<=d;w++)
                particles.X(p)(w)=vertices(v)(w);
        }        
        for( int t=1;t<=triangles.m;t++){
            render_surface->mesh.elements.Append(triangles(t));
        }
        
        render_surface->Update_Number_Nodes();

        //for( int i = 0; i <= 1; i ++ )
        //    render_surface->Linearly_Subdivide();      
        
        ARRAY<TV>& render_vertices=elastic_lattice_deformer.render_vertices;
        ARRAY<VECTOR<int,3> >& render_triangles=elastic_lattice_deformer.render_triangles;
        render_vertices.Resize(particles.X.m);
        render_triangles.Resize(render_surface->mesh.elements.m);

       for( int v=1;v<=particles.X.m;v++){
            for(int w=1;w<=d;w++)
                render_vertices(v)(w) = particles.X(v)(w);
        }        
        for( int t=1;t<=render_surface->mesh.elements.m;t++){
            render_triangles(t) = render_surface->mesh.elements(t);
        }
        
        delete render_surface;
    }


    // Compute (padded) bounding box

    VECTOR<int,3> cells;
    VECTOR<float,3> min_corner;

    EMBEDDINGTOOLS<T,d>::Get_Grid_Bounds( vertices, triangles, dx, cells, min_corner);

    // Allocate discretization object

    T lambda=elastic_lattice_deformer.youngs_modulus*elastic_lattice_deformer.poissons_ratio/
        ((1+elastic_lattice_deformer.poissons_ratio)*(1-2*elastic_lattice_deformer.poissons_ratio));
    T mu=elastic_lattice_deformer.youngs_modulus/(2*(1+elastic_lattice_deformer.poissons_ratio));
    T alpha = sqrt( mu * lambda) / dx;
    LOG::cout << "Lambda: "<< lambda << std::endl;
    LOG::cout << "Mu: "<< mu << std::endl;
    LOG::cout << "Alpha: "<< alpha << std::endl;
    LOG::cout << "h: " << dx << std::endl;
    LOG::cout << "Mu * h = " << mu * dx << std::endl;
    LOG::cout << "Alpha^2 / Kappa = " << sqr(alpha) / lambda << std::endl;
    
    //elastic_lattice_deformer.Discretization()=new PRIMARY_DEFORMER::T_DISCRETIZATION(cells,dx,min_corner);
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();
    discretization.CreateEngine(cells,dx,min_corner,mu,lambda,alpha,.25,100);
    elastic_lattice_deformer.engine_created = true;
    elastic_lattice_deformer.single_point_constraint_maxid = 0;
    elastic_lattice_deformer.dual_point_constraint_maxid = 0;
    // Rasterize object boundary (and flag respective cells as interior)
    int coarsen_ratio = refinement;
    float fine_dx = dx / coarsen_ratio;
    VECTOR<int,3> fine_cells=cells*coarsen_ratio;
    
    elastic_lattice_deformer.Initialize_Fine_Domain(fine_cells, fine_dx, min_corner, coarsen_ratio,
                                                    NONLINEAR_ELASTICITY<T,d>::Get_Lower_Cell_Padding(),
                                                    NONLINEAR_ELASTICITY<T,d>::Get_Upper_Cell_Padding(),
                                                    NONLINEAR_ELASTICITY<T,d>::Get_Lower_Node_Padding(),
                                                    NONLINEAR_ELASTICITY<T,d>::Get_Upper_Node_Padding());
    elastic_lattice_deformer.coarse_h = dx;
    elastic_lattice_deformer.coarse_grid=discretization.Grid();
    EMBEDDINGTOOLS<T,d>::Rasterize( vertices, triangles, elastic_lattice_deformer.fine_grid,
                               elastic_lattice_deformer.unpadded_fine_domain,
                               elastic_lattice_deformer.padded_fine_domain,
                               elastic_lattice_deformer.voxmap,
                               elastic_lattice_deformer.voxmap_node );
/*
 
    if(!FILE_UTILITIES::Directory_Exists("output"))
        FILE_UTILITIES::Create_Directory("output");
    if(!FILE_UTILITIES::Directory_Exists("output/common"))
        FILE_UTILITIES::Create_Directory("output/common");

    FILE_UTILITIES::Write_To_File(stream_type,"output/common/grid",elastic_lattice_deformer.fine_grid);
*/
#ifdef GRID_IN_GRID
    ARRAY<CELL_TYPE, T_INDEX> cell_type;
    discretization.GetCellType( cell_type );
    EMBEDDINGTOOLS<T,d>::Coarsen(elastic_lattice_deformer.fine_grid, elastic_lattice_deformer.unpadded_fine_domain, elastic_lattice_deformer.padded_fine_domain, elastic_lattice_deformer.voxmap,
                            discretization.Grid(), discretization.Unpadded_Cell_Domain(), discretization.Padded_Cell_Domain(), cell_type);
    discretization.SetCellType( cell_type );

    const GRID<TV>& grid=elastic_lattice_deformer.coarse_grid;
    const GRID<TV>& fine_grid=elastic_lattice_deformer.fine_grid;
    const RANGE<T_INDEX>& unpadded_cell_domain=elastic_lattice_deformer.unpadded_fine_domain;
    // Initialize embedding
    // Note: We are now embedding into the fine grid, NOT the coarser, simulation grid.
    {
        ARRAY<T_INDEX>& embedding_map=elastic_lattice_deformer.embedding_map;
        ARRAY<TV>& embedding_weights=elastic_lattice_deformer.embedding_weights;
        embedding_map.Resize(vertices.m);embedding_weights.Resize(vertices.m);
        for(int v=1;v<=vertices.m;v++){
            const TV& X=vertices(v);
            const T_INDEX cell_index=fine_grid.Cell(X,0);
            PHYSBAM_ASSERT(unpadded_cell_domain.Lazy_Inside(cell_index));
            embedding_map(v)=cell_index;
            const TV weights=(X-fine_grid.Node(cell_index))/fine_dx;
            PHYSBAM_ASSERT(weights.Min()>-1e-4 && weights.Max()<1+1e-4);
            embedding_weights(v)=weights;
        }
    }
    {
        ARRAY<TV>& render_vertices=elastic_lattice_deformer.render_vertices;
        ARRAY<VECTOR<int,3> >& render_triangles=elastic_lattice_deformer.render_triangles;
        ARRAY<T_INDEX>& embedding_map=elastic_lattice_deformer.render_embedding_map;
        ARRAY<TV>& embedding_weights=elastic_lattice_deformer.render_embedding_weights;
        embedding_map.Resize(render_vertices.m);embedding_weights.Resize(render_vertices.m);
        for(int v=1;v<=render_vertices.m;v++){
            const TV& X=render_vertices(v);
            const T_INDEX cell_index=fine_grid.Cell(X,0);
            PHYSBAM_ASSERT(unpadded_cell_domain.Lazy_Inside(cell_index));
            embedding_map(v)=cell_index;
            const TV weights=(X-fine_grid.Node(cell_index))/fine_dx;
            PHYSBAM_ASSERT(weights.Min()>-1e-4 && weights.Max()<1+1e-4);
            embedding_weights(v)=weights;
        }
    }

#else
    discretization.cell_type.Fill(CELL_TYPE(0));
    const GRID<TV>& grid=discretization.grid;
    const RANGE<T_INDEX>& unpadded_cell_domain=discretization.unpadded_cell_domain;
    for(int t=1;t<=triangles.m;t++){
        TRIANGLE_3D<T> triangle(vertices(triangles(t)(1)),vertices(triangles(t)(2)),vertices(triangles(t)(3)));
        RANGE<TV> triangle_box=triangle.Bounding_Box();
        RANGE<T_INDEX> cell_range;
        cell_range.min_corner=grid.Cell(triangle_box.min_corner,0);PHYSBAM_ASSERT(unpadded_cell_domain.Lazy_Inside(cell_range.min_corner));
        cell_range.max_corner=grid.Cell(triangle_box.max_corner,0);PHYSBAM_ASSERT(unpadded_cell_domain.Lazy_Inside(cell_range.max_corner));
        for(RANGE_ITERATOR<d> iterator(cell_range);iterator.Valid();iterator.Next()){
            const T_INDEX& cell_index=iterator.Index();
            RANGE<TV> cell=grid.Cell_Domain(cell_index);
            if(INTERSECTION::Intersects<T>(cell,triangle,0))
                discretization.cell_type(cell_index)=INTERIOR_CELL_TYPE;}}

    // Flood fill the exterior region up to the boundary

    const RANGE<T_INDEX>& padded_cell_domain=discretization.padded_cell_domain;
    QUEUE<T_INDEX> queue(unpadded_cell_domain.Surface_Area());
    for(RANGE_ITERATOR<d> iterator(padded_cell_domain);iterator.Valid();iterator.Next()){
        const T_INDEX& cell_index=iterator.Index();
        if(unpadded_cell_domain.Lazy_Outside(cell_index)){
            discretization.cell_type(cell_index)=EXTERIOR_CELL_TYPE;
            queue.Safe_Enqueue(cell_index);}}
    while(!queue.Empty()){
        const T_INDEX cell_index=queue.Dequeue();
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index-1,cell_index+1));iterator.Valid();iterator.Next()){
            const T_INDEX& neighbor_index=iterator.Index();
            if(padded_cell_domain.Lazy_Inside(neighbor_index) && discretization.cell_type(neighbor_index)==CELL_TYPE(0)){
                discretization.cell_type(neighbor_index)=EXTERIOR_CELL_TYPE;
                queue.Safe_Enqueue(neighbor_index);}}}
    for(RANGE_ITERATOR<d> iterator(padded_cell_domain);iterator.Valid();iterator.Next()){
        const T_INDEX& cell_index=iterator.Index();
        if(discretization.cell_type(cell_index)==CELL_TYPE(0))
            discretization.cell_type(cell_index)=INTERIOR_CELL_TYPE;}
#ifdef ENABLE_LOG_MESSAGES
    {int interior_count=0;
    for(RANGE_ITERATOR<d> iterator(padded_cell_domain);iterator.Valid();iterator.Next())
        if(discretization.cell_type(iterator.Index())==INTERIOR_CELL_TYPE)
            interior_count++;
    LOG::cout<<"Interior cells : "<<interior_count<<std::endl;}
#endif
    
    // Initialize embedding
    ARRAY<T_INDEX>& embedding_map=elastic_lattice_deformer.embedding_map;
    ARRAY<TV>& embedding_weights=elastic_lattice_deformer.embedding_weights;
    embedding_map.Resize(vertices.m);embedding_weights.Resize(vertices.m);
    for(int v=1;v<=vertices.m;v++){
        const TV& X=vertices(v);
        const T_INDEX cell_index=grid.Cell(X,0);
        PHYSBAM_ASSERT(unpadded_cell_domain.Lazy_Inside(cell_index));
        embedding_map(v)=cell_index;
        const TV weights=(X-grid.Node(cell_index))/dx;
        PHYSBAM_ASSERT(weights.Min()>-1e-4 && weights.Max()<1+1e-4);
        embedding_weights(v)=weights;
    }
#endif
    // Initialize and Run cutting algorithm on the fine voxmap
    
#if 1
// New cutter code

    RANGE<TV> t_range;
    t_range.min_corner = min_corner;
    t_range.max_corner = t_range.min_corner + TV(cells)*dx;
    LOG::cout << "GRID Domain: " << grid.domain << std::endl;
    LOG::cout << "GRID dx: " << grid.dX << std::endl;
    LOG::cout << "GRID counts: " << discretization.Unpadded_Cell_Domain() << std::endl;
    REGULAR_HYPERCUBE_MESH<T,d> template_mesh(dx,cells,t_range);
    LOG::cout << "Mesh Domain: " << template_mesh.domain << std::endl;
    LOG::cout << "Mesh dx: " << template_mesh.dx << std::endl;
    //LOG::cout << "Computing Neighbors..." << std::endl;
    template_mesh.Compute_Neighbors();
    CUTTER<T, d> cutter(template_mesh);
    //LOG::cout << "Creating material predicate strategy..." << std::endl;

    ARRAY< ARRAY<bool, T_INDEX> > per_cell_voxmap;

    per_cell_voxmap.Resize( template_mesh.elements.m );
    for( int i = 1; i <= per_cell_voxmap.m; i++){
        per_cell_voxmap(i).Resize( RANGE<T_INDEX>( T_INDEX::All_Ones_Vector(), T_INDEX::All_Ones_Vector() * refinement ) );
        per_cell_voxmap(i).Fill(false);
        for( RANGE_ITERATOR<d> iterator( RANGE<T_INDEX>( T_INDEX::All_Ones_Vector(), T_INDEX::All_Ones_Vector() * refinement )); iterator.Valid(); iterator.Next() ){            
            per_cell_voxmap(i)(iterator.Index()) = elastic_lattice_deformer.voxmap( (template_mesh.elements(i).index - 1) * refinement + iterator.Index());
        }
    }

    MATERIAL_PREDICATE_VOXELIZED_VOLUME<T,d> mp( refinement, per_cell_voxmap );
    //CUTTER_NONMANIFOLD_LEVELSET_STRATEGY<T,d> strategy;
    CUTTER_BASIC_STRATEGY<T,d> strategy;
    cutter.Set_Material_Predicate( &mp );  
    cutter.Set_Strategy( &strategy );
    cutter.Clear_Linkages();
    cutter.Generate();


    { 
        LOG::SCOPE scope("CLElib::Set Mesh Information from Cutter");

        HASHTABLE<T_INDEX, int > cell_location_to_root_cell;
        HASHTABLE<int, int > subcell_id_to_mesh_node_id;
        HASHTABLE<int, int > subcell_id_to_mesh_cell_id;

        ARRAY<int, T_INDEX> node_identifier; 
        ARRAY<bool, T_INDEX> node_location_is_mesh;
        ARRAY<bool, T_INDEX> cell_location_is_mesh;
        node_identifier.Resize( discretization.Unpadded_Node_Domain() ); node_identifier.Fill( -1 );
        node_location_is_mesh.Resize( discretization.Unpadded_Node_Domain() ); node_location_is_mesh.Fill( false );
        cell_location_is_mesh.Resize( discretization.Unpadded_Cell_Domain() ); cell_location_is_mesh.Fill( false );

        int mesh_cell_count = 0;
        int mesh_node_count = 0;
        int mesh_vertex_counter = 1;
        int mesh_cell_counter = 1;

        ARRAY< ARRAY<int> > root_cells_to_subcells = cutter.Root_Cell_To_Subcell();   
        const ARRAY< CUTTER_REGIONS<T,d>::T_CELL >& subcells = cutter.GetSubCells();
        {
            LOG::SCOPE scope("CLElib::Find all obvious mesh nodes (root has been duplicated)" );                        
            for( int i = 1; i <= template_mesh.elements.m; i++){
                const T_INDEX& cell_index = template_mesh.elements(i).index;

                for( int s = 1; s <= root_cells_to_subcells( i ).m; s++ ){
                    int subcell_id = root_cells_to_subcells( i )(s);
                    int flat_index=1;                
                    for( RANGE_ITERATOR<d> iterator( RANGE<T_INDEX>(cell_index,
                                                                    cell_index+1) );
                         iterator.Valid(); iterator.Next(), flat_index++ ){
                        int subcell_vertex = subcells(subcell_id).vertices(flat_index);
                        if( !(node_identifier( iterator.Index() ) == -1 || node_identifier( iterator.Index() ) == subcell_vertex) )
                            node_location_is_mesh( iterator.Index() ) = true;
                        else
                            node_identifier( iterator.Index() ) = subcell_vertex;
                    }
                }
            }
        }

        {
            LOG::SCOPE scope("CLElib::Find all obvious mesh cells (root has been duplicated)" );                        
            for( int i = 1; i <= template_mesh.elements.m; i++){
                cell_location_to_root_cell.Insert( template_mesh.elements(i).index, i );
                if( root_cells_to_subcells(i).m > 1 ){                    
                    // Mark this cell as definately mesh
                    cell_location_is_mesh(template_mesh.elements(i).index) = true;                    
                    
                    for( int s = 1; s <= root_cells_to_subcells( i ).m; s++ ){
                        int subcell_id = root_cells_to_subcells( i )(s);
                        subcell_id_to_mesh_cell_id.Insert( subcell_id, mesh_cell_counter);
                        mesh_cell_counter++;

                        int flat_index=1;
                        for( RANGE_ITERATOR<d> iterator( RANGE<T_INDEX>(template_mesh.elements(i).index,
                                                                        template_mesh.elements(i).index+1) );
                             iterator.Valid(); iterator.Next(), flat_index++ ){
                            // Mark this cells nodes as definately mesh
                            if( node_location_is_mesh( iterator.Index() ) ){
                                int subcell_vertex = subcells(subcell_id).vertices(flat_index);
                                if( !subcell_id_to_mesh_node_id.Contains( subcell_vertex )){
                                    subcell_id_to_mesh_node_id.Insert(subcell_vertex , mesh_vertex_counter );
                                    mesh_vertex_counter++;
                                    //LOG::cout << "Mapping " << subcell_vertex << " to " << mesh_vertex_counter  << std::endl;
                                }
                            }
                        }                        
                    }                    
                }                    
            }
            mesh_node_count = mesh_vertex_counter-1;
            mesh_cell_count = mesh_cell_counter-1;
            LOG::cout << "Number of mesh cells: " << mesh_cell_count << std::endl;
            LOG::cout << "Number of mesh nodes: " << mesh_node_count << std::endl;
        }
       
        {
            LOG::SCOPE scope("CLElib::Mark all cells containing mesh nodes as mesh" );                        
            for( RANGE_ITERATOR<d> iterator( discretization.Unpadded_Cell_Domain() ); iterator.Valid(); iterator.Next() ){
                int root_id = cell_location_to_root_cell.Get( iterator.Index() );

                if( cell_location_is_mesh( iterator.Index() ) == false && root_cells_to_subcells(root_id).m == 1){
                    int subcell_id = root_cells_to_subcells(root_id)(1);
                    int flat_index=1;
                    for( RANGE_ITERATOR<d> node_iterator( RANGE<T_INDEX>(iterator.Index(), iterator.Index()+1)); node_iterator.Valid(); node_iterator.Next(), flat_index++ ){
                        int subcell_vertex = subcells(subcell_id).vertices(flat_index);
                        // If one of our nodes is mesh, we must be mesh                  
                        if( node_location_is_mesh( node_iterator.Index() ) ){
                            cell_location_is_mesh( iterator.Index() ) = true;                            
                            if( !subcell_id_to_mesh_node_id.Contains( subcell_vertex )){
                                subcell_id_to_mesh_node_id.Insert(subcell_vertex , mesh_vertex_counter );
                                mesh_vertex_counter++;
                                //LOG::cout << "Mapping " << subcell_vertex << " to " << mesh_vertex_counter  << std::endl;
                            }
                        }                        
                    }
                    if( cell_location_is_mesh( iterator.Index() ) ){
                        subcell_id_to_mesh_cell_id.Insert( subcell_id, mesh_cell_counter);
                        mesh_cell_counter++;
                    }
                }
            }
            mesh_node_count = mesh_vertex_counter-1;
            mesh_cell_count = mesh_cell_counter-1;
            LOG::cout << "Number of mesh cells: " << mesh_cell_count << std::endl;
            LOG::cout << "Number of mesh nodes: " << mesh_node_count << std::endl;
        }

        {
            LOG::SCOPE scope("CLElib::Initialize mesh data structures.");
            discretization.Initialize_Mesh(mesh_cell_count,mesh_node_count);
        }

        {
            LOG::SCOPE scope("CLElib::Setting up grid cells.");
            for( RANGE_ITERATOR<d> iterator( discretization.Unpadded_Cell_Domain() ); iterator.Valid(); iterator.Next() ){
                if( cell_location_is_mesh( iterator.Index() ) )
                    cell_type(iterator.Index())=EXTERIOR_CELL_TYPE;
            }
        }
        
        

        {
            LOG::SCOPE scope("CLElib::Setting up mesh cells.");
            ARRAY<CELL_TYPE, int> cell_type_mesh;
            ARRAY<T_INDEX, int> cell_indices_mesh;
            ARRAY<VECTOR<int, 8>, int> cells_mesh;
            discretization.GetCellTypeMesh( cell_type_mesh );
            discretization.GetCellIndicesMesh( cell_indices_mesh );
            discretization.GetCellsMesh( cells_mesh );

            for( RANGE_ITERATOR<d> iterator( discretization.Unpadded_Cell_Domain() ); iterator.Valid();iterator.Next()){
                if( cell_location_is_mesh( iterator.Index() ) ){
                    if(!cell_location_to_root_cell.Contains( iterator.Index() ) ){
                        LOG::cout << "Can't find mapping from " << iterator.Index() << " into template mesh." << std::endl;
                        PHYSBAM_FATAL_ERROR();
                    }
                    int root_index = cell_location_to_root_cell.Get( iterator.Index() );
                    
                    //LOG::cout << "Root cell "<< root_index << "("<< iterator.Index() << ") has " << root_cells_to_subcells( root_index ).m << " subcells." << std::endl;

                    for( int s = 1; s <= root_cells_to_subcells( root_index ).m; s++ ){
                        int subcell_id = root_cells_to_subcells( root_index )(s);
                        VECTOR<int,8> new_cell;
                        int flat_index=1;
                        for( RANGE_ITERATOR<d> node_iterator( RANGE<T_INDEX>(iterator.Index(), iterator.Index()+1)); node_iterator.Valid(); node_iterator.Next(),flat_index++ ){
                            if( node_location_is_mesh( node_iterator.Index() ) ){
                                int subcell_vertex = subcells(subcell_id).vertices(flat_index);
                                //LOG::cout << "Loading mapped vertex for " << subcell_vertex << std::endl;
                                new_cell(flat_index) = subcell_id_to_mesh_node_id.Get( subcell_vertex );
                            }
                            else{
                                new_cell(flat_index) = 0;
                            }
                        }
                        //LOG::cout << "\tSubcell " << s << "(" << subcell_id << ") has vertices: " << new_cell << std::endl;
                        cell_type_mesh( subcell_id_to_mesh_cell_id.Get( subcell_id ) ) = INTERIOR_CELL_TYPE;
                        cell_indices_mesh( subcell_id_to_mesh_cell_id.Get( subcell_id ) ) = iterator.Index();
                        cells_mesh( subcell_id_to_mesh_cell_id.Get( subcell_id ) ) = new_cell;                        
                    }
                }
            }
            discretization.SetCellType( cell_type );
            discretization.SetCellTypeMesh( cell_type_mesh );
            discretization.SetCellIndicesMesh( cell_indices_mesh );
            discretization.SetCellsMesh( cells_mesh );
            PHYSBAM_ASSERT( mesh_vertex_counter == (mesh_node_count+1) );
        }

        // Set fine mesh mapping
        {
            LOG::SCOPE scope("CLElib::Initialize Fine Mapping");
            ARRAY< T, T_INDEX >& fine_to_coarsemesh=elastic_lattice_deformer.fine_to_coarsemesh;           

            fine_to_coarsemesh.Fill(0);
            for( RANGE_ITERATOR<d> iterator( discretization.Unpadded_Cell_Domain() ); iterator.Valid(); iterator.Next() ){
                int root_id = cell_location_to_root_cell.Get( iterator.Index() );
                if( cell_location_is_mesh( iterator.Index() ) ){
                    for( int s = 1; s <= root_cells_to_subcells( root_id ).m; s++ ){
                        int subcell_id = root_cells_to_subcells( root_id )(s);
                        const ARRAY<bool, T_INDEX>& fragment_voxmap = mp.SubcellVoxMaterialFragment(subcell_id);
                        for( RANGE_ITERATOR<d> fiterator( RANGE<T_INDEX>( T_INDEX::All_Ones_Vector(), T_INDEX::All_Ones_Vector() * refinement )); fiterator.Valid(); fiterator.Next() ){            
                            if( fragment_voxmap(fiterator.Index()) )
                                fine_to_coarsemesh( (iterator.Index()-1)*refinement + fiterator.Index() ) = subcell_id_to_mesh_cell_id.Get( subcell_id );
                        }
                    }
                    
                }
            }
        }
        


    }



#else
    elastic_lattice_deformer.vrg = new VOXELIZED_REGION_GENERATOR<T,d>( coarsen_ratio,
                                                                        elastic_lattice_deformer.fine_grid,
                                                                        elastic_lattice_deformer.voxmap );
    elastic_lattice_deformer.vrg->Generate();

    // Fill out the mesh stuff here!!!
#ifdef GRID_IN_GRID
    {
        LOG::SCOPE scope("CLElib::Set Mesh Information from Cutter");

        int grid_cell_count=0;
        int grid_node_count=0;
        int mesh_cell_count=0;
        int mesh_node_count=0;
        
        int mesh_cell_estimate=0;
        int mesh_node_estimate=0; 

        ARRAY< T, T_INDEX >& fine_to_coarsemesh=elastic_lattice_deformer.fine_to_coarsemesh;              
        HASHTABLE<int,int> node_index_reducer;

        {
            LOG::SCOPE scope("Counting Grid and Mesh cells..." );
            // Determine total number of mesh cells and nodes needed.
            for(RANGE_ITERATOR<d> iterator(discretization.Unpadded_Cell_Domain());iterator.Valid();iterator.Next()){
                const T_INDEX& index=iterator.Index();
                // This may need to change
                bool is_mesh=false;
                bool is_grid=false;
                if( (cell_type(index)==INTERIOR_CELL_TYPE || cell_type(index)==BOUNDARY_CELL_TYPE) &&
                    elastic_lattice_deformer.vrg->IsMeshMappable( index ) ){
                    is_mesh = true;
                }
                if( (cell_type(index)==INTERIOR_CELL_TYPE || cell_type(index)==BOUNDARY_CELL_TYPE) &&
                    !elastic_lattice_deformer.vrg->IsMeshMappable( index ) ){
                    is_grid = true;
                }

                if(is_mesh) for( int dup=1; dup <= elastic_lattice_deformer.vrg->DuplicatesAtCoarseIndex( index ); dup++){
                        RANGE<T_INDEX> nodes( index, index+1);
                        mesh_cell_count++;         
                        int i=1;
                        VECTOR<int,8> cutter_vertices = elastic_lattice_deformer.vrg->GetCellVertices( index, dup );
                        for( RANGE_ITERATOR<d> node_itr(nodes); node_itr.Valid(); node_itr.Next(), i++){
                            int mesh_node = cutter_vertices(i);
                            if( !node_index_reducer.Contains( mesh_node ) ){
                                mesh_node_count++;
                                node_index_reducer.Insert(mesh_node,mesh_node_count);}
                        }
                    }
                if(is_grid)
                    grid_cell_count++;
            }

            LOG::cout << "Number of grid cells needed: " << grid_cell_count << std::endl;
            LOG::cout << "Number of grid nodes needed: " << grid_cell_count * 8 << std::endl;
        
            LOG::cout << "Number of mesh cells needed: " << mesh_cell_count << std::endl;
            LOG::cout << "Number of mesh nodes needed: " << mesh_node_count << std::endl;
        }

        discretization.Initialize_Mesh(mesh_cell_count,mesh_node_count);

        mesh_cell_estimate = mesh_cell_count;
        mesh_node_estimate = mesh_node_count;
        mesh_cell_count = 0;
        mesh_node_count = 0;

        ARRAY<CELL_TYPE, int> cell_type_mesh;
        ARRAY<T_INDEX, int> cell_indices_mesh;
        ARRAY<VECTOR<int, 8>, int> cells_mesh;
        discretization.GetCellTypeMesh( cell_type_mesh );
        discretization.GetCellIndicesMesh( cell_indices_mesh );
        discretization.GetCellsMesh( cells_mesh );
        int Number_Of_Mesh_Cells = discretization.Number_Of_Mesh_Cells();
        int Number_Of_Mesh_Nodes = discretization.Number_Of_Mesh_Nodes();
        
        LOG::cout << cell_type_mesh.m << std::endl;
        LOG::cout << Number_Of_Mesh_Cells << std::endl;
        PHYSBAM_ASSERT( cell_type_mesh.m == Number_Of_Mesh_Cells );


        {
            LOG::SCOPE scope( "Setting up mesh cells...");
            // Set up mesh cells
            for(RANGE_ITERATOR<d> iterator(discretization.Unpadded_Cell_Domain());iterator.Valid();iterator.Next()){
                const T_INDEX& index=iterator.Index();
                // This may need to change
                bool is_mesh=false;
                if( (cell_type(index)==INTERIOR_CELL_TYPE || cell_type(index)==BOUNDARY_CELL_TYPE) )
                    if(elastic_lattice_deformer.vrg->IsMeshMappable( index ) )
                        is_mesh = true;
            

                if(is_mesh){
                    // Unset grid data
                    cell_type(index)=EXTERIOR_CELL_TYPE;
                
                    for( int dup=1; dup <= elastic_lattice_deformer.vrg->DuplicatesAtCoarseIndex( index ); dup++){
                        RANGE<T_INDEX> nodes( index, index+1);
                        mesh_cell_count++;         
                        PHYSBAM_ASSERT(  mesh_cell_count > 0 && mesh_cell_count <= mesh_cell_estimate );
                        VECTOR<int,8> cutter_vertices = elastic_lattice_deformer.vrg->GetCellVertices( index, dup );
            
                        // Set mesh data
                        cell_type_mesh(mesh_cell_count)=INTERIOR_CELL_TYPE;
                        cell_indices_mesh(mesh_cell_count)=index;
                        cells_mesh(mesh_cell_count)=cutter_vertices;

                        // Set fine mesh mapping
                        T_INDEX fine_base_index = ((index - 1) * refinement) + 1;
                        for(RANGE_ITERATOR<d> fine_iterator(RANGE<T_INDEX>(fine_base_index,fine_base_index+refinement-1));fine_iterator.Valid();fine_iterator.Next()){
                            const T_INDEX& fine_index = fine_iterator.Index();
                            if( elastic_lattice_deformer.voxmap( fine_index ) == 1){
                                T_INDEX vox_cell;
                                for( int v=1; v<=d; v++)
                                    vox_cell(v) = (fine_index(v) % refinement) == 0 ? refinement : (fine_index(v) % refinement);
                                if(elastic_lattice_deformer.vrg->GetSubcellVoxmap(index,dup)(vox_cell)){
                                    //LOG::cout << "Mapping fine cell " << fine_index << " to coarse mesh cell " << mesh_cell_count << std::endl;
                                    fine_to_coarsemesh(fine_index) = mesh_cell_count;
                                }
                            }
                        }
                    }
                }
            }
        }
#if 0
        // Do a sanity check here: For every mesh cell, there should be no more than REF^3 sub cells that match it.
        HASHTABLE<int,int> check_table;
        for(RANGE_ITERATOR<d> iterator(elastic_lattice_deformer.unpadded_fine_domain);iterator.Valid();iterator.Next()){
            if( fine_to_coarsemesh(iterator.Index()) > 0 ){
                int& check_val = check_table.Get_Or_Insert(fine_to_coarsemesh(iterator.Index()),0);
                check_val += 1;
                PHYSBAM_ASSERT( check_val <=  ( refinement * refinement * refinement ) );
            }
        }

        // Which cells map to mesh 1
        for(RANGE_ITERATOR<d> iterator(elastic_lattice_deformer.unpadded_fine_domain);iterator.Valid();iterator.Next()){
            if( fine_to_coarsemesh(iterator.Index()) == 1 ){
                LOG::cout << iterator.Index() << std::endl;
            }
        }
        
        LOG::cout << fine_to_coarsemesh(T_INDEX(99,40,66)) << std::endl;
        LOG::cout << fine_to_coarsemesh(T_INDEX(104,40,131)) << std::endl;
#endif        
        
        // Set up mesh nodes

        {
            LOG::SCOPE scope( "Setting up mesh nodes...");
            OPERATION_HASH<int> mesh_mappable_nodes(Number_Of_Mesh_Nodes);
            mesh_mappable_nodes.Next_Operation();
            OPERATION_HASH<int> mesh_non_mappable_nodes(Number_Of_Mesh_Nodes);
            mesh_non_mappable_nodes.Next_Operation();

            for( int m = 1; m <= Number_Of_Mesh_Cells; m++){
                const T_INDEX& index = cell_indices_mesh(m);

                int flat_index=1;
                for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>(index,index+1));
                    node_iterator.Valid();node_iterator.Next(),flat_index++){
                    const T_INDEX& node_index = node_iterator.Index();

                    bool is_mesh_mappable = true;
                
                    for(RANGE_ITERATOR<d> node_neighbor_iterator(RANGE<T_INDEX>(node_index-1,node_index));
                        node_neighbor_iterator.Valid();node_neighbor_iterator.Next())
                        if( cell_type( node_neighbor_iterator.Index() ) == INTERIOR_CELL_TYPE || 
                            cell_type( node_neighbor_iterator.Index() ) == BOUNDARY_CELL_TYPE )
                            // This 'mesh' node is directly touching an active grid cell. It must be made grid.
                            is_mesh_mappable = false;
                    
                    if(!is_mesh_mappable){
                        mesh_non_mappable_nodes.Mark(node_index_reducer.Get(cells_mesh(m)(flat_index)));
                        cells_mesh(m)(flat_index) = 0;
                    }
                    else{
                        mesh_mappable_nodes.Mark(node_index_reducer.Get(cells_mesh(m)(flat_index)));
                        cells_mesh(m)(flat_index) = node_index_reducer.Get(cells_mesh(m)(flat_index));
                        PHYSBAM_ASSERT( cells_mesh(m)(flat_index) <= Number_Of_Mesh_Nodes && cells_mesh(m)(flat_index) > 0) ;
                    }

                }
            }
            
            {
                int mappable_count = 0;  int non_mappable_count = 0;
                for( int m=1; m <= Number_Of_Mesh_Nodes; m++){
                    PHYSBAM_ASSERT( mesh_mappable_nodes.Is_Marked_Current(m)  || 
                                    mesh_non_mappable_nodes.Is_Marked_Current(m) );
                    if( mesh_mappable_nodes.Is_Marked_Current(m) )
                        mappable_count++;
                    if( mesh_non_mappable_nodes.Is_Marked_Current(m) )
                        non_mappable_count++;
                }         
                
                
                LOG::cout << "Out of a total " << Number_Of_Mesh_Nodes << " mesh nodes: " << std::endl;
                LOG::cout << "\t" << mappable_count << " remain as mesh nodes," << std::endl;
                LOG::cout << "\t" << non_mappable_count <<" must be assigned to the grid." << std::endl; 
            }
        }


        discretization.SetCellType( cell_type );
        discretization.SetCellTypeMesh( cell_type_mesh );
        discretization.SetCellIndicesMesh( cell_indices_mesh );
        discretization.SetCellsMesh( cells_mesh );

    }
#endif
#endif

    // Done with mesh initialization

    // Create Muscle Data structures   
    discretization.Initialize_Muscles();
    //int num_of_muscles = 0;
    //discretization.muscle_fiber_max_stresses.Resize(num_of_muscles);
    //discretization.muscle_activations.Resize(num_of_muscles);


    // Flag domain as yet uninitialized, need to call Set_Fixed_Triangles to complete initialization
    elastic_lattice_deformer.domain_initialized=false;

}
//#####################################################################
// Function SetTextureParameters
//#####################################################################
void CLElib::
SetTextureParameters(const std::vector<int>& dimensions, const std::vector<float>& texcoords, const std::vector<int>& textriangles, const std::vector<int>& coord_map ){
    LOG::SCOPE scope("CLElib::SetTextureParameters()");

    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("GetTextureStress() called while no dynamic model is already created!");

    // Don't call this once the model is deforming! This will mess up the vertex positions.
    ARRAY<VECTOR<int,d> > triangles;
    for( int i=0; i<textriangles.size(); i+=3)
        triangles.Append( VECTOR<int,d>(textriangles.at(i),textriangles.at(i+1), textriangles.at(i+2)));
    ARRAY<TV>& vertices=elastic_lattice_deformer.vertices;
    
    const int TEXTURE_COORD_SIZE = 2;

    // Build Box hierarchy around texture space triangles for faster lookup
    BOX_HIERARCHY< VECTOR<T,2> > texture_bh;
    ARRAY< RANGE<VECTOR<T,2> > > textri_boxes;
    ARRAY< TRIANGLE_2D<T> > texture_triangles;
    for( int t=1; t <= triangles.m; t++ ){
//        LOG::cout << "Building triangle box for T" << t << ":  ";
        RANGE<VECTOR<T,2> > box;
        VECTOR< VECTOR<T, 2>, 3 > points;
        for( int i=1; i<=d; i++){
            points(i) = VECTOR<T,2>( texcoords[ (triangles(t)(i))*TEXTURE_COORD_SIZE + 0 ],
                                     texcoords[ (triangles(t)(i))*TEXTURE_COORD_SIZE + 1 ] );
            box.Enlarge_To_Include_Point( points(i) );
        }
        textri_boxes.Append( box );
        TRIANGLE_2D<T> Tri(points);
        texture_triangles.Append(Tri);
//        LOG::cout << texture_triangles.Last().Signed_Area() << "  " << textri_boxes.Last() << std::endl;
    }
    texture_bh.Set_Leaf_Boxes( textri_boxes, true );

    RANGE< VECTOR<int,2> >& texture_dimensions = elastic_lattice_deformer.texture_dimensions;
    ARRAY< T_INDEX, VECTOR<int, 2> >& pixel_embedding_map = elastic_lattice_deformer.pixel_embedding_map;
    ARRAY< TV, VECTOR<int, 2> >& pixel_embedding_weights = elastic_lattice_deformer.pixel_embedding_weights;
    RANGE< T_INDEX >& unpadded_cell_domain = elastic_lattice_deformer.unpadded_fine_domain;
    const GRID<TV>& fine_grid=elastic_lattice_deformer.fine_grid;
    
    texture_dimensions= RANGE<VECTOR<int,2> >(VECTOR<int,2>(0,0), VECTOR<int,2>(dimensions[0]-1,dimensions[1]-1));
    pixel_embedding_map.Resize( texture_dimensions );
    pixel_embedding_weights.Resize( texture_dimensions );
    
    int mapped_pixels = 0;
    for( RANGE_ITERATOR<2> iter( texture_dimensions ); iter.Valid(); iter.Next() ){        
        const VECTOR<int,2>& pixel_index = iter.Index();
        VECTOR<T,2> pixel_location( pixel_index(1) / float(dimensions[0]), pixel_index(2) / float(dimensions[1]));

        //LOG::cout << "Testing pixel (" << pixel_index << ") at texture location " << pixel_location << std::endl;

        int embedding_triangle = -1;
        TV embedding_coords = TV();

        // Find the triangle we belong too
        ARRAY<int> hits;
        texture_bh.Intersection_List(pixel_location, hits, 0 );
        //LOG::cout << "Potential Source Triangles: " << hits.m << std::endl;
        for( int h = 1; h <= hits.m && embedding_triangle == -1; h++){
            TRIANGLE_2D<T>& ttriangle = texture_triangles(hits(h));
            TRIANGLE_2D<T> clean_ttriangle = ttriangle;
            clean_ttriangle.Fix_Orientation();
            if( clean_ttriangle.Signed_Area() > 0 &&  clean_ttriangle.Inside( pixel_location ) ) {
                //LOG::cout << "Triangle "<<hits(h)<<" is a match!" << std::endl;
                embedding_triangle = hits(h);
                embedding_coords = ttriangle.Barycentric_Coordinates( pixel_location,
                                                                      ttriangle.X(1),
                                                                      ttriangle.X(2),
                                                                      ttriangle.X(3) );                
            }
        }

        if( embedding_triangle == -1 ){
            pixel_embedding_map(pixel_index)=T_INDEX(-1,-1,-1);
            pixel_embedding_weights(pixel_index)=TV(-1,-1,-1);
            continue; // We are not embedded in a triangle, skip this next part...
        }

        // Locate the pixel in 3D space...
        TRIANGLE_3D<T> spatial_triangle;
        spatial_triangle.x1 = vertices(coord_map.at(triangles(embedding_triangle)(1))+1);
        spatial_triangle.x2 = vertices(coord_map.at(triangles(embedding_triangle)(2))+1);
        spatial_triangle.x3 = vertices(coord_map.at(triangles(embedding_triangle)(3))+1);
        TV spatial_point = spatial_triangle.Point_From_Barycentric_Coordinates(embedding_coords);
        
        // And find the embedding relationship to the fine grid...
        {
            const T_INDEX cell_index=fine_grid.Cell(spatial_point,0);
            PHYSBAM_ASSERT(unpadded_cell_domain.Lazy_Inside(cell_index));
            pixel_embedding_map(pixel_index)=cell_index;
            const TV weights=(spatial_point-fine_grid.Node(cell_index))/elastic_lattice_deformer.h;
            PHYSBAM_ASSERT(weights.Min()>-1e-4 && weights.Max()<1+1e-4);
            pixel_embedding_weights(pixel_index)=weights;
            mapped_pixels++;
        }

    }
    
    LOG::cout << "Mapped " << mapped_pixels << " to the fine grid." << std::endl;


    

}
//#####################################################################
// Function GetTextureStress
//#####################################################################
void  CLElib::
GetTextureStress(std::vector<float>& texdata ){
//#ifdef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::GetTextureStress()");
//#endif
    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("GetTextureStress() called while no dynamic model is already created!");

    RANGE< VECTOR<int,2> >& texture_dimensions = elastic_lattice_deformer.texture_dimensions;
    ARRAY< T_INDEX, VECTOR<int, 2> >& pixel_embedding_map = elastic_lattice_deformer.pixel_embedding_map;
    ARRAY< TV, VECTOR<int, 2> >& pixel_embedding_weights = elastic_lattice_deformer.pixel_embedding_weights;
    RANGE< T_INDEX >& unpadded_cell_domain = elastic_lattice_deformer.unpadded_fine_domain;
    const GRID<TV>& fine_grid=elastic_lattice_deformer.fine_grid;
    
    texdata.erase( texdata.begin(), texdata.end() );
    texdata.reserve( texture_dimensions.Size() * 3);
    int colored_pixels = 0;
    for( RANGE_ITERATOR<2> iter( texture_dimensions ); iter.Valid(); iter.Next() ){
        const VECTOR<int,2>& pixel_index = iter.Index();

        if( pixel_embedding_map(pixel_index) == T_INDEX(-1,-1,-1) ){
            // add a black color for unmapped pixels
            texdata.push_back(0.0);
            texdata.push_back(0.0);
            texdata.push_back(0.0);
        }
        else{
            colored_pixels++;
            // add a white color for mapped pixels
            texdata.push_back( elastic_lattice_deformer.u_fine_stress(1)(pixel_embedding_map(pixel_index)) );
            texdata.push_back( elastic_lattice_deformer.u_fine_stress(2)(pixel_embedding_map(pixel_index)) );
            texdata.push_back( elastic_lattice_deformer.u_fine_stress(3)(pixel_embedding_map(pixel_index)) );
        }        
    }          
    LOG::cout << "Created " << colored_pixels << " active data pixels." << std::endl;
    
}
//#####################################################################
// Function Add_Static_Model
//#####################################################################
void CLElib::
//Create_Model(const std::string& cell_info_filename,const std::string& object_name, const std::string& deformation_filename)
Add_Static_Model(const std::vector<float>& vertices_input,const std::vector<int>& triangles_input)
{
#ifdef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Add_Static_Model()");
#endif


    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Add_Static_Model() called while no dynamic model is already created!");

    ARRAY< ARRAY<TV> >& static_model_vertices=elastic_lattice_deformer.static_model_vertices;
    ARRAY< ARRAY<VECTOR<int,3> > >& static_model_triangles=elastic_lattice_deformer.static_model_triangles;

    int last_static_model = static_model_vertices.m;
    int current_static_model = last_static_model+1;
    //  TV average_vertex;

    static_model_vertices.Append( ARRAY<TV>() );
    static_model_triangles.Append( ARRAY<VECTOR<int,3> >() );

    PHYSBAM_ASSERT(vertices_input.size()%3==0);
    static_model_vertices(current_static_model).Resize(vertices_input.size()/3);
    for(int i=0,p=1;p<=static_model_vertices(current_static_model).m;p++)
        {
            for(int v=1;v<=3;v++){
//                average_vertex(v)+=vertices_input[i];
                static_model_vertices(current_static_model)(p)(v)=vertices_input[i++];
            }
            //LOG::cout << static_model_vertices(current_static_model)(p) << std::endl;
        }


    // average_vertex /= static_model_vertices(current_static_model).m;
    // LOG::cout << "Static_Object" <<  current_static_model << " centroid point: " << average_vertex << std::endl;


    PHYSBAM_ASSERT(triangles_input.size()%3==0);
    static_model_triangles(current_static_model).Resize(triangles_input.size()/3);
    for(int i=0,t=1;t<=static_model_triangles(current_static_model).m;t++)
        for(int v=1;v<=3;v++)
            static_model_triangles(current_static_model)(t)(v)=triangles_input[i++]+1;

}
//#####################################################################
// Function Add_Muscle_Layer
//#####################################################################
void CLElib::
Add_Muscle_Layer(const std::vector<float>& vertices,const std::vector<int>& triangles, const std::vector<float>& fiber_direction, const float& maxstress)
{
    LOG::SCOPE scope("CLElib::Add_Muscle_Layer()");

    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef bool T_FLAG;
    typedef ARRAY<T_FLAG, T_INDEX> T_FLAG_ARRAY;

    typedef ARRAY<int, T_INDEX> ID_ARRAY;
    typedef ARRAY<int, int> ID_ARRAY_MESH;
    typedef ARRAY<float, T_INDEX> DENSITY_ARRAY;
    typedef ARRAY<float, int> DENSITY_ARRAY_MESH;
    typedef ARRAY<TV, T_INDEX> FIBER_ARRAY;
    typedef ARRAY<TV, int> FIBER_ARRAY_MESH;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Add_Static_Model() called while no dynamic model is already created!");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();

    ARRAY< TV > muscle_vertices;
    ARRAY< VECTOR<int,3> > muscle_triangles;
    PHYSBAM_ASSERT(vertices.size()%3==0);
    PHYSBAM_ASSERT(triangles.size()%3==0);

    muscle_vertices.Resize( vertices.size()/3 );
    muscle_triangles.Resize( triangles.size()/3 );
    for(int i=0,p=1;p<=muscle_vertices.m;p++)
        for(int v=1;v<=3;v++)
            muscle_vertices(p)(v)=vertices[i++];

    for(int i=0,p=1;p<=muscle_triangles.m;p++)
        for(int v=1;v<=3;v++)
            muscle_triangles(p)(v)=triangles[i++]+1;
        
    TV fiber_dir;
    fiber_dir(1) = fiber_direction[0];
    fiber_dir(2) = fiber_direction[1];
    fiber_dir(3) = fiber_direction[2];

    // Sanity checks

    PHYSBAM_ASSERT(ARRAYS_COMPUTATIONS::Max(muscle_triangles.Flattened())<=muscle_vertices.m);
    PHYSBAM_ASSERT(ARRAYS_COMPUTATIONS::Min(muscle_triangles.Flattened())>0);
    for(int t=1;t<=muscle_triangles.m;t++){
        const VECTOR<int,3>& tri=muscle_triangles(t);
        PHYSBAM_ASSERT(tri(1)!=tri(2) && tri(2)!=tri(3) && tri(3)!=tri(1));}

    {
        TRIANGLE_MESH triangle_mesh(muscle_vertices.m,muscle_triangles);
        triangle_mesh.Initialize_Node_On_Boundary();
        PHYSBAM_ASSERT(!triangle_mesh.node_on_boundary->Contains(true));
    }


    // Rasterize muscle to find extents
    T_FLAG_ARRAY muscle_voxmap;
    T_FLAG_ARRAY muscle_voxmap_node;
    muscle_voxmap.Resize(elastic_lattice_deformer.padded_fine_domain); muscle_voxmap.Fill(false);
    muscle_voxmap_node.Resize(elastic_lattice_deformer.padded_fine_node_domain); muscle_voxmap.Fill(false);
    

    EMBEDDINGTOOLS<T,d>::Rasterize( muscle_vertices, muscle_triangles, elastic_lattice_deformer.fine_grid,
                               elastic_lattice_deformer.unpadded_fine_domain,
                               elastic_lattice_deformer.padded_fine_domain,
                               muscle_voxmap,
                               muscle_voxmap_node );

    // Compute muscle density per coarse cell
    DENSITY_ARRAY muscle_density;
    DENSITY_ARRAY_MESH muscle_density_mesh;
    muscle_density.Resize( discretization.Padded_Cell_Domain() ); muscle_density.Fill(0.0);
    muscle_density_mesh.Resize( discretization.Number_Of_Mesh_Cells() ); muscle_density_mesh.Fill(0.0);


    EMBEDDINGTOOLS<T,d>::Coarsen_Density(elastic_lattice_deformer.fine_grid,
                                    elastic_lattice_deformer.unpadded_fine_domain,
                                    elastic_lattice_deformer.padded_fine_domain,
                                    muscle_voxmap,
                                    discretization.Grid(), discretization.Unpadded_Cell_Domain(),
                                    discretization.Padded_Cell_Domain(), muscle_density);

    // Assign fiber direction per coarse cell
    FIBER_ARRAY muscle_fiber;
    FIBER_ARRAY_MESH muscle_fiber_mesh;
    muscle_fiber.Resize( discretization.Padded_Cell_Domain() ); muscle_fiber.Fill(TV());
    muscle_fiber_mesh.Resize( discretization.Number_Of_Mesh_Cells() ); muscle_fiber_mesh.Fill(TV());

    //TODO: Fix me!
    // This is not ideal. We need to also support setting fiber direction as directed by the 
    // normal of the embedded object. This will allow us to properly handle curved surfaces.
    // For this we need to grab a levelset of the embedded object....
    for( RANGE_ITERATOR<d> iterator(discretization.Unpadded_Cell_Domain());iterator.Valid();iterator.Next()){
        const T_INDEX& index = iterator.Index();
        muscle_fiber(index) = fiber_dir;
    }
    
    // Generate Muscle ID's
    ID_ARRAY muscle_ids; 
    ID_ARRAY_MESH muscle_ids_mesh;
    muscle_ids.Resize( discretization.Padded_Cell_Domain() ); muscle_ids.Fill( -1 );
    muscle_ids_mesh.Resize( discretization.Number_Of_Mesh_Cells() ); muscle_ids_mesh.Fill( -1 );
    for( RANGE_ITERATOR<d> iterator(discretization.Unpadded_Cell_Domain());iterator.Valid();iterator.Next()){
        const T_INDEX& index = iterator.Index();
        if( muscle_density( index ) > 0 )
            muscle_ids(index) = 1;
    }  

    ARRAY<T_INDEX, int> cell_indices_mesh;
    discretization.GetCellIndicesMesh( cell_indices_mesh );
    int Number_Of_Mesh_Cells = discretization.Number_Of_Mesh_Cells();
    for( int cell = 1; cell < Number_Of_Mesh_Cells; cell++){
        const T_INDEX& index = cell_indices_mesh(cell);

        muscle_ids_mesh(cell) = muscle_ids(index);
        // TODO: Fix me!
        // This is not accurate. Each mesh cell needs to have correct distribution of density
        // which depends on which fine voxel is in which mesh cell. It should be okay for now.
        muscle_density_mesh(cell) = muscle_density(index);

        muscle_fiber_mesh(cell) = muscle_fiber(index);
    }

    elastic_lattice_deformer.Update_Muscles( muscle_ids, muscle_ids_mesh, muscle_density, muscle_density_mesh, muscle_fiber, muscle_fiber_mesh, maxstress);
 
}

//#####################################################################
// Function Set_Collision_Model
//#####################################################################
void CLElib::
Set_Collision_Model(const std::vector<float>& vertices_input,const std::vector<int>& triangles_input)
{
    LOG::SCOPE scope("CLElib::Set_Collision_Model()");

    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Set_Collision_Model() called while no dynamic model is already created!");

    PHYSBAM_ASSERT(vertices_input.size()%3==0);
    PHYSBAM_ASSERT(triangles_input.size()%3==0);

    std::vector< std::vector<float> > vertices;
    std::vector< std::vector<int> > triangles;
    
    for( int i=0; i < vertices_input.size()/3; i++){
        std::vector<float> vertex;
        vertex.push_back( vertices_input[i*3+0] );
        vertex.push_back( vertices_input[i*3+1] );
        vertex.push_back( vertices_input[i*3+2] );
        vertices.push_back( vertex );
    }

    for( int i=0; i < triangles_input.size()/3; i++){
        std::vector<int> triangle;
        triangle.push_back( triangles_input[i*3+0] );
        triangle.push_back( triangles_input[i*3+1] );
        triangle.push_back( triangles_input[i*3+2] );
        triangles.push_back( triangle );
    }

    PHYSBAM_LEVELSET_COLLISION<T,d>* levelset_collision;    
    {
        LOG::SCOPE scope("CLElib::Set_Collision_Model()::Initializing Collision Levelset");
        levelset_collision=
            new PHYSBAM_LEVELSET_COLLISION<T,d>(vertices,
                                                triangles,
                                                elastic_lattice_deformer.coarse_h * .1);
        
    }
    
    elastic_lattice_deformer.Set_Collision_Object(levelset_collision);
}

//#####################################################################
// Function Set_Collision_Proxy
//#####################################################################
void CLElib::
Set_Collision_Proxy(const std::vector<float>& vertices_input,const std::vector<int>& triangles_input)
{
    LOG::SCOPE scope("CLElib::Set_Collision_Model()");

    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Set_Collision_Proxy() called while no dynamic model is already created!");

    PHYSBAM_ASSERT(vertices_input.size()%3==0);
    PHYSBAM_ASSERT(triangles_input.size()%3==0);

// Translate from std::vector to PhysBAM arrays (converting 0-base to 1-base)

    PHYSBAM_ASSERT(vertices_input.size()%3==0);
    ARRAY<TV>& collision_proxy_vertices=elastic_lattice_deformer.collision_proxy_vertices;
    collision_proxy_vertices.Resize(vertices_input.size()/3);
    for(int i=0,p=1;p<=collision_proxy_vertices.m;p++)
        for(int v=1;v<=3;v++){
            collision_proxy_vertices(p)(v)=vertices_input[i];
            i++;
        }

    PHYSBAM_ASSERT(triangles_input.size()%3==0);
    ARRAY<VECTOR<int,3> >& collision_proxy_triangles=elastic_lattice_deformer.collision_proxy_triangles;
    collision_proxy_triangles.Resize(triangles_input.size()/3);
    for(int i=0,t=1;t<=collision_proxy_triangles.m;t++)
        for(int v=1;v<=3;v++)
            collision_proxy_triangles(t)(v)=triangles_input[i++]+1;


    T dx = elastic_lattice_deformer.h;
    int coarsen_ratio = refinement;
    float fine_dx = dx / coarsen_ratio;

    const GRID<TV>& fine_grid=elastic_lattice_deformer.fine_grid;
    const RANGE<T_INDEX>& unpadded_cell_domain=elastic_lattice_deformer.unpadded_fine_domain;
    // Initialize embedding
    ARRAY<T_INDEX>& embedding_map=elastic_lattice_deformer.collision_proxy_embedding_map;
    embedding_map.Resize(collision_proxy_vertices.m);
    for(int v=1;v<=collision_proxy_vertices.m;v++){
        const TV& X=collision_proxy_vertices(v);
        const T_INDEX cell_index=fine_grid.Cell(X,0);
        PHYSBAM_ASSERT(unpadded_cell_domain.Lazy_Inside(cell_index));
        embedding_map(v)=cell_index;
    }

    elastic_lattice_deformer.collision_proxy_set = true;

}

//#####################################################################
// Function Destroy_Model
//#####################################################################
void CLElib::
Destroy_Model()
{
#ifdef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Destroy_Model()");
#endif

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    elastic_lattice_deformer.Destroy_Engine();
}
//#####################################################################
// Function Set_Fixed_Geometry
//#####################################################################
void CLElib::
Set_Fixed_Geometry(const std::vector<float>& vertices, const std::vector<int>& points, const std::vector<int>& segments, const std::vector<int>& triangles)
{
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Set_Fixed_Geometry()");
#endif

    PHYSBAM_ASSERT( vertices.size()%3 == 0 );

    std::vector<float> point_vertices;
    for( std::vector<int>::const_iterator iter = points.begin(); iter != points.end(); iter++){
        const int& vertex = *(iter);
        PHYSBAM_ASSERT( ((vertex+1)*3)-1 <= vertices.size() );
        for( int i = vertex; i < vertex+3; i++)
            point_vertices.push_back( vertices[i] );
    }
    Set_Fixed_Points(point_vertices);
    Set_Fixed_Segments(vertices, segments);
    Set_Fixed_Triangles(vertices, triangles);
}

//#####################################################################
// Function Set_Fixed_Points
//#####################################################################
void CLElib::
Set_Fixed_Points(const std::vector<float>& vertices_input)
{
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Set_Fixed_Points()");
#endif

    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Set_Fixed_Points() called without first creating a model");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();
    if(elastic_lattice_deformer.domain_initialized)
        PHYSBAM_FATAL_ERROR("Set_Fixed_Points() called after model final initializtion");
    
    PHYSBAM_ASSERT(vertices_input.size()%3==0);
    std::vector< std::vector<float> > vertices;   
    for( int i=0; i < vertices_input.size()/3; i++){
        std::vector<float> vertex;
        vertex.push_back( vertices_input[i*3+0] );
        vertex.push_back( vertices_input[i*3+1] );
        vertex.push_back( vertices_input[i*3+2] );
        vertices.push_back( vertex );
    }

    
    // Rasterize fixed segments onto dirichlet cells

    const GRID<TV>& grid=elastic_lattice_deformer.fine_grid;
    const RANGE<T_INDEX>& unpadded_cell_domain=elastic_lattice_deformer.unpadded_fine_domain;

    std::vector<std::vector<int> > marked_cells;
    for(unsigned int i=0;i<vertices.size();i++){
        TV a(vertices[i][0], vertices[i][1], vertices[i][2]);
        const T_INDEX& cell_index = grid.Cell(a,0);
        elastic_lattice_deformer.voxmap_dirichlet(cell_index)=true;
        marked_cells.push_back(std::vector<int>(cell_index.begin(), cell_index.end()));
    }

    UpdateDirichletCells(marked_cells);
}
//#####################################################################
// Function Set_Fixed_Segments
//#####################################################################
void CLElib::
Set_Fixed_Segments(const std::vector<float>& vertices_input, const std::vector<int>& segments_input)
{
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Set_Fixed_Segments()");
#endif

    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Set_Fixed_Segments() called without first creating a model");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();
    if(elastic_lattice_deformer.domain_initialized)
        PHYSBAM_FATAL_ERROR("Set_Fixed_Segments() called after model final initializtion");
    
    PHYSBAM_ASSERT(vertices_input.size()%3==0);
    PHYSBAM_ASSERT(segments_input.size()%2==0);

    std::vector< std::vector<float> > vertices;
    std::vector< std::vector<int> > segments;
    
    for( int i=0; i < vertices_input.size()/3; i++){
        std::vector<float> vertex;
        vertex.push_back( vertices_input[i*3+0] );
        vertex.push_back( vertices_input[i*3+1] );
        vertex.push_back( vertices_input[i*3+2] );
        vertices.push_back( vertex );
    }

    for( int i=0; i < segments_input.size()/2; i++){
        std::vector<int> segment;
        segment.push_back( segments_input[i*2+0] );
        segment.push_back( segments_input[i*2+1] );
        segments.push_back( segment );
    }
    // Rasterize fixed segments onto dirichlet cells

    const GRID<TV>& grid=elastic_lattice_deformer.fine_grid;
    const RANGE<T_INDEX>& unpadded_cell_domain=elastic_lattice_deformer.unpadded_fine_domain;

    std::vector<std::vector<int> > marked_cells;
    for(unsigned int i=0;i<segments.size();i++){
        TV a(vertices[segments[i][0]][0], vertices[segments[i][0]][1], vertices[segments[i][0]][2]);
        TV b(vertices[segments[i][1]][0], vertices[segments[i][1]][1], vertices[segments[i][1]][2]);
        SEGMENT_3D<T> s(a,b);
        RAY<TV> segment(s);
        //LOG::cout << "Setting Fixed Segment: " << s.x1 << " " << s.x2 << std::endl;
        //LOG::cout << "Collision Ray: " << segment << std::endl;
        RANGE<TV> segment_box(a,b);
        RANGE<T_INDEX> cell_range;
        cell_range.min_corner=grid.Cell(segment_box.min_corner,0);
        cell_range.max_corner=grid.Cell(segment_box.max_corner,0);
        if(!cell_range.Lazy_Intersection(unpadded_cell_domain))
            {
                LOG::cout << "Segment " << i << " is out of simulation bounds." << std::endl;
                continue;
            }
        cell_range = RANGE<T_INDEX>::Intersect(cell_range, unpadded_cell_domain);

        //LOG::cout << "Segment Cell Range: " << cell_range << std::endl;
        for(RANGE_ITERATOR<d> iterator(cell_range);iterator.Valid();iterator.Next()){
            const T_INDEX& cell_index=iterator.Index();
            RANGE<TV> cell=grid.Cell_Domain(cell_index);
            //LOG::cout << "Does cell " << cell << " intersect segment? " << std::endl;
            RAY<TV> temp(s);
            bool intersects = INTERSECTION::Intersects<T>(temp,cell,0);
            if(intersects){
                //LOG::cout << "YES!" << std::endl;
                elastic_lattice_deformer.voxmap_dirichlet(cell_index)=true;
                marked_cells.push_back(std::vector<int>(cell_index.begin(), cell_index.end()));
            }
        }
    }

    UpdateDirichletCells(marked_cells);
    
}
//#####################################################################
// Function Set_Fixed_Triangles
//#####################################################################
void CLElib::
Set_Fixed_Triangles(const std::vector<float>& vertices_input, const std::vector<int>& triangles_input)
{
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Set_Fixed_Triangles()");
#endif
    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Set_Fixed_Triangles() called without first creating a model");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();
    if(elastic_lattice_deformer.domain_initialized)
        PHYSBAM_FATAL_ERROR("Set_Fixed_Triangles() called after model final initializtion");
    
    PHYSBAM_ASSERT(vertices_input.size()%3==0);
    PHYSBAM_ASSERT(triangles_input.size()%3==0);

    std::vector< std::vector<float> > vertices;
    std::vector< std::vector<int> > triangles;
    
    for( int i=0; i < vertices_input.size()/3; i++){
        std::vector<float> vertex;
        vertex.push_back( vertices_input[i*3+0] );
        vertex.push_back( vertices_input[i*3+1] );
        vertex.push_back( vertices_input[i*3+2] );
        vertices.push_back( vertex );
    }

    for( int i=0; i < triangles_input.size()/3; i++){
        std::vector<int> triangle;
        triangle.push_back( triangles_input[i*3+0] );
        triangle.push_back( triangles_input[i*3+1] );
        triangle.push_back( triangles_input[i*3+2] );
        triangles.push_back( triangle );
    }
    // Rasterize fixed triangles onto dirichlet cells

    const GRID<TV>& grid=elastic_lattice_deformer.fine_grid;
    const RANGE<T_INDEX>& unpadded_cell_domain=elastic_lattice_deformer.unpadded_fine_domain;

    std::vector<std::vector<int> > marked_cells;
    for(unsigned int i=0;i<triangles.size();i++){
        TV a(vertices[triangles[i][0]][0], vertices[triangles[i][0]][1], vertices[triangles[i][0]][2]);
        TV b(vertices[triangles[i][1]][0], vertices[triangles[i][1]][1], vertices[triangles[i][1]][2]);
        TV c(vertices[triangles[i][2]][0], vertices[triangles[i][2]][1], vertices[triangles[i][2]][2]);
        TRIANGLE_3D<T> triangle(a,b,c);
        RANGE<TV> triangle_box=triangle.Bounding_Box();
        RANGE<T_INDEX> cell_range;
        cell_range.min_corner=grid.Cell(triangle_box.min_corner,0);
        cell_range.max_corner=grid.Cell(triangle_box.max_corner,0);
        if(!cell_range.Lazy_Intersection(unpadded_cell_domain))
            {
                LOG::cout << "Triangle " << i << " is out of simulation bounds." << std::endl;
                continue;
            }
        cell_range = RANGE<T_INDEX>::Intersect(cell_range, unpadded_cell_domain);

        for(RANGE_ITERATOR<d> iterator(cell_range);iterator.Valid();iterator.Next()){
            const T_INDEX& cell_index=iterator.Index();
            RANGE<TV> cell=grid.Cell_Domain(cell_index);
            TRIANGLE_3D<T> temp(a,b,c);
            if(INTERSECTION::Intersects<T>(cell,temp,0)){
                elastic_lattice_deformer.voxmap_dirichlet(cell_index)=true;
                marked_cells.push_back(std::vector<int>(cell_index.begin(), cell_index.end()));
            }
        }
    }

    UpdateDirichletCells(marked_cells);
    
    // Prepare domain for simulation
}
//#####################################################################
// Function Set_Fixed_Volume
//#####################################################################
void CLElib::
Set_Fixed_Volume(const std::vector<float>& min_corner_input, const std::vector<float>& max_corner_input)
{
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Set_Fixed_Volume()");
#endif
    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Set_Fixed_Triangles() called without first creating a model");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();
    if(elastic_lattice_deformer.domain_initialized)
        PHYSBAM_FATAL_ERROR("Set_Fixed_Triangles() called after model final initializtion");
    
    PHYSBAM_ASSERT(min_corner_input.size()%3==0);
    PHYSBAM_ASSERT(max_corner_input.size()%3==0);

    RANGE<TV> volume_box;
    volume_box.min_corner = TV(min_corner_input[0], min_corner_input[1], min_corner_input[2]);
    volume_box.max_corner = TV(max_corner_input[0], max_corner_input[1], max_corner_input[2]);

    const GRID<TV>& grid=elastic_lattice_deformer.fine_grid;
    const RANGE<T_INDEX>& unpadded_cell_domain=elastic_lattice_deformer.unpadded_fine_domain;

    std::vector<std::vector<int> > marked_cells;
    RANGE<T_INDEX> cell_range;
    cell_range.min_corner=grid.Cell(volume_box.min_corner,0);
    cell_range.max_corner=grid.Cell(volume_box.max_corner,0);
    if(!cell_range.Lazy_Intersection(unpadded_cell_domain))
        {
            LOG::cout << "Volume is out of simulation bounds." << std::endl;
            return;
        }

    cell_range = RANGE<T_INDEX>::Intersect(cell_range, unpadded_cell_domain);
    
    for(RANGE_ITERATOR<d> iterator(cell_range);iterator.Valid();iterator.Next()){
        const T_INDEX& cell_index=iterator.Index();
        elastic_lattice_deformer.voxmap_dirichlet(cell_index)=true;
        marked_cells.push_back(std::vector<int>(cell_index.begin(), cell_index.end()));        
    }

    UpdateDirichletCells(marked_cells);
    
    // Prepare domain for simulation
}
//#####################################################################
// Function Set_Fixed_Triangles
//#####################################################################
void CLElib::
Finalize_Initialization()
{
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Finalize_Initialization()");
#endif
    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Finalize_Initialization() called without first creating a model");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();
    if(elastic_lattice_deformer.domain_initialized)
        PHYSBAM_FATAL_ERROR("Finalize_Initialization() called twice for the same model instance");
    
	#if 0
    {    
        LOGGING_ELASTICITY<PRIMARY_DEFORMER::T_DISCRETIZATION>  logging;
        logging.ExportElasticityData(discretization, "Export.New.bundle");
    }
	#endif


    elastic_lattice_deformer.Initialize_Collision_Points();
    elastic_lattice_deformer.Initialize_Solver();
    elastic_lattice_deformer.domain_initialized=true;
}
//#####################################################################
// Function WriteDebug
//#####################################################################
void CLElib::
WriteDebug(bool andDie)
{
#if 0
        STREAM_TYPE stream_type((float)0);

        //Write_Output< PRIMARY_DEFORMER::T_DISCRETIZATION >(stream_type, *this, "debug", debug_output++ );

        if(andDie)
            exit(2);
#endif
}
//#####################################################################
// Function Set_Damping
//#####################################################################
void CLElib::
Set_Damping(const float rayleigh_coefficient)
{
#ifdef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Set_Damping()");
#endif

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    elastic_lattice_deformer.rayleigh_coefficient=rayleigh_coefficient;
}
//#####################################################################
// Function Set_Density
//#####################################################################
void CLElib::
Set_Density(const float density)
{
#ifdef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Set_Density()");
#endif

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    elastic_lattice_deformer.density=density;
}
//#####################################################################
// Function Set_Time_Step
//#####################################################################
void CLElib::
Set_Time_Step(const float dt)
{
#ifdef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Set_Time_Step()");
#endif

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    elastic_lattice_deformer.dt=dt;
}

//#####################################################################
// Function Set_Hook_Stiffness
//#####################################################################
void CLElib::
Set_Hook_Stiffness(const float hook_stiffness)
{
#ifdef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Set_Hook_Stiffness()");
#endif

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    elastic_lattice_deformer.hook_stiffness=hook_stiffness;
}
//#####################################################################
// Function Set_Suture_Stiffness
//#####################################################################
void CLElib::
Set_Suture_Stiffness(const float suture_stiffness)
{
#ifdef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Set_Suture_Stiffness()");
#endif

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    elastic_lattice_deformer.suture_stiffness=suture_stiffness;
}

//#####################################################################
// Single Point Constraint Functions
//#####################################################################
#include "CLElib_Single_Point_Constraints.cpp"

//#####################################################################
// Dual Point Constraint Functions
//#####################################################################
#include "CLElib_Dual_Point_Constraints.cpp"


//#####################################################################
// Function Simulate_One_Time_Step
//#####################################################################
void CLElib::
Advance_One_Time_Step()
{
#ifdef ENABLE_LOG_MESSAGES
   LOG::SCOPE scope("CLElib::Simulate_One_Time_Step()");
#endif
    STREAM_TYPE stream_type((float)0);
    
    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
   
    
    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Simulate_One_Time_Step() called without first creating a model");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();
    if(!elastic_lattice_deformer.domain_initialized)
        PHYSBAM_FATAL_ERROR("Simulate_One_Time_Step() called before completing model initialization");
    
    //if(frame_counter < 0) {
    //    Write_Output<PRIMARY_DEFORMER::T_DISCRETIZATION>(stream_type, *this, "output", 0, 0 );
    //    frame_counter++;
    //    std::cout << std::endl << int(frame_counter / 10) << "   " << discretization.Constraint_Count() <<std::endl;
    //}

    const int write_rate = 100;
    const int needed_constraints = 0;


    if( true ){
        //discretization.Print_Force_Diagnostics(true);
        //std::cout << std::endl << int(frame_counter / write_rate) << "   " << discretization.Constraint_Count() <<std::endl;

        //discretization.Exact_Solve(50,1,1e-5,5e-5, true);
        Schedule_Exact_Solve(50,2,1e-2f,5e-2f,frame_counter);

        //Update_Fine_Displacement();

        //if( frame_counter % write_rate == 0 )
        //    {
        //        LOG::SCOPE scope("Saving Physbam scene.");
        //        Write_Output(stream_type, *this, "NE-Fast-output", int(frame_counter / write_rate)+1 );
        //    }
    }
  


}
//#####################################################################
// Function: Schedule Exact Solve
//#####################################################################
namespace {
class Exact_Solve_Task : public PhysBAM::PTHREAD_QUEUE::TASK
{
    CLElib& _clelib;
    typedef float T;
    static const int d=3;
    
    int _ni;
    int _ki;
    float _nt;
    float _kt;
    int _frame;
public:
    Exact_Solve_Task(CLElib& clelib,  int ki, int ni,  float kt, float nt, int frame ) : _clelib(clelib), _ni(ni), _nt(nt), _ki(ki), _kt(kt), _frame(frame){}
    void Run() { 
        STREAM_TYPE stream_type((float)0);
        LOG::SCOPE scope("Exact_Solve_Task::RUN");
        PRIMARY_DEFORMER& elastic_lattice_deformer=_clelib.Implementation< PRIMARY_DEFORMER >();
        _clelib.Update_Collisions();
        elastic_lattice_deformer.Exact_Solve(_ki, _ni, _kt, _nt, true);
        _clelib.Update_Fine_Displacement();
        _clelib.Update_Embedded_Surfaces();        
        elastic_lattice_deformer.SwapGeometryBuffers();
        if( false )
            Write_Output_Embedded( stream_type, elastic_lattice_deformer, "output", _frame );
        
        pthread_mutex_unlock(&elastic_lattice_deformer.simulate_lock);
    }
};
}

void CLElib::Schedule_Exact_Solve(int krylov_iterations, int newton_iterations,  float krylov_tolerance, float newton_tolerance, long int& frame)
{
    typedef float T;
    static const int d=3;

    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation< PRIMARY_DEFORMER >();

    int lock_result = pthread_mutex_trylock(&elastic_lattice_deformer.simulate_lock);
    if( lock_result == 0 ){
        frame++;
        Exact_Solve_Task* task =
            new Exact_Solve_Task(*this, krylov_iterations, newton_iterations, krylov_tolerance, newton_tolerance, frame);
        elastic_lattice_deformer.simulation_queue->Queue(task);
    }
    else if( lock_result == EBUSY )
        return;
    else {
        PHYSBAM_FATAL_ERROR( "Could not aqquire simulate lock, but it was not busy ");
    }  
}

void CLElib::WaitForSolve(){
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation< PRIMARY_DEFORMER >();
    elastic_lattice_deformer.simulation_queue->Wait();
}
                         
//#####################################################################
// Function Generate_UFine_MeshMap
//#####################################################################
void CLElib::
Generate_UFine_MeshMap()
{
#if 0
#ifdef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Generate_UFine_MeshMap");
#endif

    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Update_Fine_Displacement() called while no model has been created");
    const PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();
    const ARRAY<T_INDEX>& embedding_map=elastic_lattice_deformer.render_embedding_map;
    ARRAY< T, T_INDEX >& fine_to_coarsemesh=elastic_lattice_deformer.fine_to_coarsemesh;
    const VOXELIZED_REGION_GENERATOR<T,d>& vrg=*(elastic_lattice_deformer.vrg);

    fine_to_coarsemesh.Resize(elastic_lattice_deformer.unpadded_fine_domain, true, false, -1);

    for( RANGE_ITERATOR<d> iterator(elastic_lattice_deformer.unpadded_fine_domain);
         iterator.Valid(); iterator.Next() )
        {
            const T_INDEX& index = iterator.Index();
            if( elastic_lattice_deformer.voxmap( index ) == 1){
                TV fine_node = elastic_lattice_deformer.fine_grid.Node( index );
                T_INDEX coarse_cell = ((index - 1) / refinement) + 1;
                int mesh_cells_here = elastic_lattice_deformer.vrg->DuplicatesAtCoarseIndex(coarse_cell);
                if( mesh_cells_here > 1 ){
                    T_INDEX vox_cell;
                    for( int v=1; v<=d; v++)
                        vox_cell(v) = (index(v) % refinement) == 0 ? refinement : (index(v) % refinement);
                    
                    for( int dup=1; dup <= discretization.mesh_to_cell_rmap( coarse_cell ).m; dup++){
                        if( discretization.mesh_cell_voxmaps(discretization.mesh_to_cell_rmap( coarse_cell )(dup))(vox_cell) ){
                            fine_to_coarsemesh( index ) = discretization.mesh_to_cell_rmap( coarse_cell )(dup);
                            break;
                        }            
                    }
                    if( fine_to_coarsemesh(index) == -1 )
                        PHYSBAM_FATAL_ERROR( "Oops! Ran out of mesh cells that this location could belong to!.");
                }
            }
        }
#endif
}
                         
//#####################################################################
// Function Update_Fine_Displacement
//#####################################################################
void CLElib::
Update_Fine_Displacement()
{
//#ifdef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Update_Fine_Displacement");
//#endif
    
    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Update_Fine_Displacement() called while no model has been created");
    const PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();
    const ARRAY<T_INDEX>& embedding_map=elastic_lattice_deformer.render_embedding_map;
    const ARRAY< T, T_INDEX >& fine_to_coarsemesh=elastic_lattice_deformer.fine_to_coarsemesh;
    RANGE< VECTOR<int,2> >& texture_dimensions = elastic_lattice_deformer.texture_dimensions;
    ARRAY< T_INDEX, VECTOR<int, 2> >& pixel_embedding_map = elastic_lattice_deformer.pixel_embedding_map;
    ARRAY< TV, VECTOR<int, 2> >& pixel_embedding_weights = elastic_lattice_deformer.pixel_embedding_weights;

    ARRAY< PAIR< T_INDEX, T_INDEX > > update_nodes;

#if 0
    // Update all fine cells
    for( RANGE_ITERATOR<d> fine_grid_iterator( elastic_lattice_deformer.unpadded_fine_node_domain ); fine_grid_iterator.Valid(); fine_grid_iterator.Next() ){
        const T_INDEX& index = fine_grid_iterator.Index();
        if( elastic_lattice_deformer.voxmap( index ) == 1)
            for(RANGE_ITERATOR<d> iter(RANGE<T_INDEX>(index,index+1));iter.Valid();iter.Next())
                update_nodes.Append( PAIR<T_INDEX,T_INDEX>(index,iter.Index()) );
    }
#endif

    // Update all fine cells with embedded vertices
    for( int i=1;i<=elastic_lattice_deformer.render_vertices.m;i++){
        const T_INDEX& index = embedding_map(i);    
        if( elastic_lattice_deformer.voxmap( index ) == 1)
            for(RANGE_ITERATOR<d> iter(RANGE<T_INDEX>(index,index+1));iter.Valid();iter.Next())
                update_nodes.Append( PAIR<T_INDEX,T_INDEX>(index,iter.Index()) );
    }

    // Update all fine cells with embedded pixels
    for( RANGE_ITERATOR<2> iter( texture_dimensions ); iter.Valid(); iter.Next() ){
        const VECTOR<int,2>& pixel_index = iter.Index();
        const T_INDEX& index = pixel_embedding_map(pixel_index);
        if( index != T_INDEX(-1,-1,-1) && elastic_lattice_deformer.voxmap( index ) == 1 ){
            for(RANGE_ITERATOR<d> iter(RANGE<T_INDEX>(index,index+1));iter.Valid();iter.Next())
                update_nodes.Append( PAIR<T_INDEX,T_INDEX>(index,iter.Index()) );
        }
    }
    update_nodes.Prune_Duplicates();
    
    ARRAY< TRIPLE<int, T_INDEX, TV> > update_queue;
    for(int node = 1; node <= update_nodes.m; node++){
        TV node_location;
        const T_INDEX cellindex = update_nodes(node).x;
        const T_INDEX index = update_nodes(node).y;
        TV fine_node = elastic_lattice_deformer.fine_grid.Node( index );
        T_INDEX coarse_cell = ((cellindex - 1) / refinement) + 1;
        TV fine_weights = ( fine_node - elastic_lattice_deformer.coarse_grid.Node(coarse_cell))/elastic_lattice_deformer.coarse_h;
        PHYSBAM_ASSERT(fine_weights.Min()>-1e-4 && fine_weights.Max()<1+1e-4);
        if( fine_to_coarsemesh( cellindex ) == 0 ){
            update_queue.Append( TRIPLE<int, T_INDEX, TV>( 0, coarse_cell, fine_weights) );
        }
        else{
            update_queue.Append( TRIPLE<int, T_INDEX, TV>( fine_to_coarsemesh( cellindex ), T_INDEX(),fine_weights));
        }
    }
    ARRAY<TV> updates;
    
    // Update Displacements
    discretization.Displacement( update_queue, updates );
    for(int node = 1; node <= update_nodes.m; node++){
        const T_INDEX index = update_nodes(node).y;
        for( int v=1;v<=d;v++) elastic_lattice_deformer.u_fine(v)(index) = updates(node)(v);
    }        

    // Update Stresss
    discretization.Stress( update_queue, updates );
    for(int node = 1; node <= update_nodes.m; node++){
        const T_INDEX index = update_nodes(node).y;
        for( int v=1;v<=d;v++) elastic_lattice_deformer.u_fine_stress(v)(index) = updates(node)(v);
    }  

    // Update Strains
    discretization.Strain( update_queue, updates );
    for(int node = 1; node <= update_nodes.m; node++){
        const T_INDEX index = update_nodes(node).y;
        for( int v=1;v<=d;v++) elastic_lattice_deformer.u_fine_strain(v)(index) = updates(node)(v);
    }        


}

//#####################################################################


//#####################################################################
// Function Update_Embedded_Surfaces
//#####################################################################
void CLElib::
Update_Embedded_Surfaces()
{
//#ifdef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Update_Embedded_Surfaces");
//#endif
    
    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Update_Embedded_Surfaces() called while no model has been created");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();

    // Update from embedding relation
    {
        const ARRAY<T_INDEX>& embedding_map=elastic_lattice_deformer.embedding_map;
        const ARRAY<TV>& embedding_weights=elastic_lattice_deformer.embedding_weights;
        
        for( int v=1;v<=elastic_lattice_deformer.working_vertices->m;v++){
            VECTOR<T,d>& X=elastic_lattice_deformer.working_vertices->operator()(v);                
#ifdef GRID_IN_GRID
            X=EMBEDDINGTOOLS<T,d>::Deformation( elastic_lattice_deformer.fine_grid,
                                           elastic_lattice_deformer.u_fine,
                                           embedding_map(v),
                                           embedding_weights(v),
                                           false );
#else            
            X=discretization.Deformation(embedding_map(v),embedding_weights(v));
#endif
        }
    }

    // Update stress values
    {
        const ARRAY<T_INDEX>& embedding_map=elastic_lattice_deformer.embedding_map;
        const ARRAY<TV>& embedding_weights=elastic_lattice_deformer.embedding_weights;
        
        for( int v=1;v<=elastic_lattice_deformer.working_stress->m;v++){
            VECTOR<T,d>& X=elastic_lattice_deformer.working_stress->operator()(v);                
            X=EMBEDDINGTOOLS<T,d>::Deformation( elastic_lattice_deformer.fine_grid,
                                                elastic_lattice_deformer.u_fine_stress,
                                                embedding_map(v),
                                                embedding_weights(v),
                                                true );
        }
    }

    // Update strain values
    {
        const ARRAY<T_INDEX>& embedding_map=elastic_lattice_deformer.embedding_map;
        const ARRAY<TV>& embedding_weights=elastic_lattice_deformer.embedding_weights;
        
        for( int v=1;v<=elastic_lattice_deformer.working_strain->m;v++){
            VECTOR<T,d>& X=elastic_lattice_deformer.working_strain->operator()(v);                
            X=EMBEDDINGTOOLS<T,d>::Deformation( elastic_lattice_deformer.fine_grid,
                                                elastic_lattice_deformer.u_fine_strain,
                                                embedding_map(v),
                                                embedding_weights(v),
                                                true );
        }
    }

    // Update from embedding relation
#if 0
    {
        const ARRAY<T_INDEX>& embedding_map=elastic_lattice_deformer.render_embedding_map;
        const ARRAY<TV>& embedding_weights=elastic_lattice_deformer.render_embedding_weights;
        
        for( int v=1;v<=elastic_lattice_deformer.render_vertices.m;v++){
            VECTOR<T,d>& X=elastic_lattice_deformer.render_vertices(v);                
#ifdef GRID_IN_GRID
            X=EMBEDDINGTOOLS<T,d>::Deformation( elastic_lattice_deformer.fine_grid,
                                           elastic_lattice_deformer.u_fine,
                                           embedding_map(v),
                                           embedding_weights(v),
                                           false );
#else            
            X=discretization.Deformation(embedding_map(v),embedding_weights(v));
#endif
        }
    }
#endif

}
//#####################################################################
// Function Update_Embedded_Surfaces
//#####################################################################
void CLElib::
UpdateDirichletCells(const std::vector<std::vector<int> >& cells)
{   
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Update_Dirichlet_Cells()");
#endif
    typedef float T;
    static const int d=3;
    typedef VECTOR<int,d> T_INDEX;
    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Update_Embedded_Surfaces() called while no model has been created");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();

    int added_grid_cells = 0;
    int added_mesh_cells = 0;

    ARRAY<CELL_TYPE, T_INDEX> cell_type;
    ARRAY<CELL_TYPE, int> cell_type_mesh;
    discretization.GetCellType( cell_type );
    discretization.GetCellTypeMesh( cell_type_mesh );

    LOG::cout << "Existing Dirichlet Grid Cells: " << cell_type.array.Count_Matches(DIRICHLET_CELL_TYPE) << std::endl;
    LOG::cout << "Existing Dirichlet Mesh Cells: " << cell_type_mesh.Count_Matches(DIRICHLET_CELL_TYPE) << std::endl;

    LOG::cout << "Moving " << cells.size() << " fine cells to coarse cells..." << std::endl;

    const ARRAY< T, T_INDEX >& fine_to_coarsemesh=elastic_lattice_deformer.fine_to_coarsemesh;
    typedef std::vector<std::vector<int> >::const_iterator cell_iterator;
    for(cell_iterator iter=cells.begin(); iter != cells.end(); iter++){
        T_INDEX cell_index( (*iter)[0], (*iter)[1], (*iter)[2] );
        if(/*elastic_lattice_deformer.voxmap(cell_index) &&*/ elastic_lattice_deformer.voxmap_dirichlet(cell_index)){
            if( fine_to_coarsemesh( cell_index ) == 0 ){
                T_INDEX coarse_cell = ((cell_index - 1) / refinement) + 1;
                if(cell_type(coarse_cell) == INTERIOR_CELL_TYPE){
                    cell_type(coarse_cell) =DIRICHLET_CELL_TYPE;
                    added_grid_cells++;
                }
            }
            else{
                if(cell_type_mesh(fine_to_coarsemesh(cell_index)) == INTERIOR_CELL_TYPE){
                    cell_type_mesh(fine_to_coarsemesh(cell_index)) =DIRICHLET_CELL_TYPE;
                    added_mesh_cells++;
                }
            }
        }
    }

    LOG::cout << "Added " << added_grid_cells << " dirichlet grid cells." << std::endl;
    LOG::cout << "Added " << added_mesh_cells << " dirichlet mesh cells." << std::endl;

    LOG::cout << "After Dirichlet Grid Cells: " << cell_type.array.Count_Matches(DIRICHLET_CELL_TYPE) << std::endl;
    LOG::cout << "After Dirichlet Mesh Cells: " << cell_type_mesh.Count_Matches(DIRICHLET_CELL_TYPE) << std::endl;
    
    discretization.SetCellType( cell_type );
    discretization.SetCellTypeMesh( cell_type_mesh );

}
//#####################################################################


//#####################################################################
// Function Get_Vertices
//#####################################################################
void CLElib::
Get_Vertices(std::vector<float>& vertices_output) const
{
#ifndef ENABLE_LOG_MESSAGES
    // LOG::SCOPE scope("CLElib::Get_Vertices()");
#endif
    
    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    const PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<const PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Get_Vertices() called while no model has been created");
    const PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.ConstDiscretization();

    // Update from embedding relation

    pthread_mutex_lock(&elastic_lattice_deformer.geometry_buffer_lock);

    const ARRAY<TV>* vertices=elastic_lattice_deformer.display_vertices;
    const ARRAY<T_INDEX>& embedding_map=elastic_lattice_deformer.embedding_map;
    const ARRAY<TV>& embedding_weights=elastic_lattice_deformer.embedding_weights;

    LOG::cout << "Get_Vertices(): Reserving " << (vertices->m *3) << " more floats for vertex data." << std::endl;
    vertices_output.reserve( vertices_output.size() + (vertices->m *3) );
    for(int i=0,v=1;v<=vertices->m;v++){
        const VECTOR<T,d>& X=vertices->operator()(v);
        for(int w=1;w<=d;w++)
            vertices_output.push_back(X(w));
    }

    pthread_mutex_unlock(&elastic_lattice_deformer.geometry_buffer_lock);
}

//#####################################################################
// Function Get_Stress
//#####################################################################
void CLElib::
Get_Stress(std::vector<float>& stress_output) const
{
#ifndef ENABLE_LOG_MESSAGES
    // LOG::SCOPE scope("CLElib::Get_Vertices()");
#endif
    
    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    const PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<const PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Get_Stress() called while no model has been created");
    const PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.ConstDiscretization();

    // Update from embedding relation

    pthread_mutex_lock(&elastic_lattice_deformer.geometry_buffer_lock);

    const ARRAY<TV>* stress=elastic_lattice_deformer.display_stress;
    const ARRAY<T_INDEX>& embedding_map=elastic_lattice_deformer.embedding_map;
    const ARRAY<TV>& embedding_weights=elastic_lattice_deformer.embedding_weights;

    LOG::cout << "Get_Stress(): Reserving " << (stress->m *3) << " more floats for stress data." << std::endl;
    int start = stress_output.size();
    stress_output.reserve( stress_output.size() + (stress->m *3) );
    for(int i=0,v=1;v<=stress->m;v++){
        const VECTOR<T,d>& X=stress->operator()(v);
        for(int w=1;w<=d;w++){
            stress_output.push_back(X(w));
        }
    }

    pthread_mutex_unlock(&elastic_lattice_deformer.geometry_buffer_lock);
}

//#####################################################################
// Function Get_Strain
//#####################################################################
void CLElib::
Get_Strain(std::vector<float>& strain_output) const
{
#ifndef ENABLE_LOG_MESSAGES
    // LOG::SCOPE scope("CLElib::Get_Vertices()");
#endif
    
    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    const PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<const PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Get_Strain() called while no model has been created");
    const PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.ConstDiscretization();

    // Update from embedding relation

    pthread_mutex_lock(&elastic_lattice_deformer.geometry_buffer_lock);

    const ARRAY<TV>* strain=elastic_lattice_deformer.display_strain;
    const ARRAY<T_INDEX>& embedding_map=elastic_lattice_deformer.embedding_map;
    const ARRAY<TV>& embedding_weights=elastic_lattice_deformer.embedding_weights;

    LOG::cout << "Get_Strain(): Reserving " << (strain->m *3) << " more floats for strain data." << std::endl;
    int start = strain_output.size();
    strain_output.reserve( strain_output.size() + (strain->m *3) );
    for(int i=0,v=1;v<=strain->m;v++){
        const VECTOR<T,d>& X=strain->operator()(v);
        for(int w=1;w<=d;w++){
            strain_output.push_back(X(w));
        }
    }

    pthread_mutex_unlock(&elastic_lattice_deformer.geometry_buffer_lock);
}

//#####################################################################
// Function Get_Vertices
//#####################################################################
void CLElib::
Get_Vertices(std::vector<float>& vertices_output, int& sinceFrame) const
{
#ifndef ENABLE_LOG_MESSAGES
    // LOG::SCOPE scope("CLElib::Get_Vertices()");
#endif
    
    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    const PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<const PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Get_Vertices() called while no model has been created");
    const PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.ConstDiscretization();

    //TODO: This might be wrong. Check for risk of race conditions
    vertices_output.erase(vertices_output.begin(), vertices_output.end());
    if( elastic_lattice_deformer.frame > sinceFrame )
        Get_Vertices(vertices_output);   
    
    sinceFrame = elastic_lattice_deformer.frame;
}

//#####################################################################
// Function Get_Vertex_Data
//#####################################################################
void CLElib::
Get_Vertex_Data(std::vector<float>& vdata, int type, int& sinceFrame) const
{
#ifndef ENABLE_LOG_MESSAGES
    // LOG::SCOPE scope("CLElib::Get_Vertices()");
#endif
    
    PHYSBAM_ASSERT( type & ( V_DATA_TYPE::POSITION | V_DATA_TYPE::STRESS | V_DATA_TYPE::STRAIN ) != 0 );

    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    const PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<const PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Get_Vertices() called while no model has been created");
    const PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.ConstDiscretization();

    vdata.erase(vdata.begin(), vdata.end());
    if( elastic_lattice_deformer.frame > sinceFrame ){
        
        if( type & V_DATA_TYPE::POSITION )
            Get_Vertices(vdata);   

        if( type & V_DATA_TYPE::STRESS )
            Get_Stress(vdata);   

        if( type & V_DATA_TYPE::STRAIN )
            Get_Strain(vdata);           
    }   
    
    sinceFrame = elastic_lattice_deformer.frame;

}

//#####################################################################
// Function Apply_Perturbation
//#####################################################################
void CLElib::
Apply_Perturbation(const float perturb_amount)
{
    #if 0
    LOG::SCOPE scope("Generate random perturbation");
    
    
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef ARRAY<T,T_INDEX> T_SCALAR_VARIABLE;
    typedef VECTOR<T_SCALAR_VARIABLE,d> T_VECTOR_VARIABLE;

    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<const PRIMARY_DEFORMER>();
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();

    T_VECTOR_VARIABLE du;
    for(int v=1;v<=d;v++) du(v).Resize(discretization.padded_node_domain,true,false);

    RANDOM_NUMBERS<T> random_numbers;random_numbers.Set_Seed(1);
    for(RANGE_ITERATOR<d> iterator(discretization.Unpadded_Node_Domain());iterator.Valid();iterator.Next()){
		const T_INDEX& node_index=iterator.Index();
		if(discretization.node_is_active(node_index))
			for(int v=1;v<=d;v++)
				du(v)(node_index)=random_numbers.Get_Uniform_Number(-perturb_amount,perturb_amount);}

    elastic_lattice_deformer.u_internal.x += du;
    #endif
}
//#####################################################################
// Function Update_Collisions
//#####################################################################
void CLElib::
Update_Collisions()
{
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef ARRAY<T,T_INDEX> T_SCALAR_VARIABLE;
    typedef VECTOR<T_SCALAR_VARIABLE,d> T_VECTOR_VARIABLE;

    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Get_Vertices() called while no model has been created");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();

    if( !elastic_lattice_deformer.collision_shape )
        return;

    LOG::SCOPE scope("CLElib:: Update_Collision_Constraints()");

    ARRAY<T> collision_constants;
    ARRAY<TV> collision_locations;

    discretization.GetCollisionConstants( collision_constants );
    discretization.GetCollisionLocations( collision_locations );
    int number_of_collisions = discretization.NumberOfConstraints(discretization.COLLISION);

    ARRAY< TRIPLE< int, T_INDEX, TV > > update_queue;
    ARRAY< TV > updates;
    ARRAY< int > collision_references;

    update_queue.Resize( number_of_collisions );
    collision_references.Resize( number_of_collisions );

    CONSTRAINT_SEGMENT<T,d> pc;
    for( int m=1; m<=number_of_collisions; m++){        
        discretization.GetConstraint(discretization.COLLISION, m, pc);
        PHYSBAM_ASSERT(pc.is_reference == true );
        collision_references(m) = pc.spring_coefficient_ptr;
        if(pc.endpoints[1].type == CONSTRAINT_NODE<T,d>::GRID_FIXED)
            update_queue(m) = TRIPLE<int,T_INDEX,TV>( 0, pc.endpoints[1].grid_index(),
                                                      pc.endpoints[1].multilinear_coordinates() );
        else if(pc.endpoints[1].type == CONSTRAINT_NODE<T,d>::MESH_FIXED)
            update_queue(m) = TRIPLE<int,T_INDEX,TV>( pc.endpoints[1].mesh_index(), T_INDEX(),
                                                      pc.endpoints[1].multilinear_coordinates() );           
        else
            PHYSBAM_FATAL_ERROR();
    }
    discretization.Deformation( update_queue, updates );
    for( int m=1; m<=number_of_collisions; m++){        
        T depth = elastic_lattice_deformer.collision_shape->Phi(updates(m));
        if( depth < 0 ) {
            collision_locations(collision_references(m)) = 
                elastic_lattice_deformer.collision_shape->ClosestPoint(updates(m));
            collision_constants(collision_references(m)) = 1e4;
        }
        else{
            collision_locations(collision_references(m)) = TV();
            collision_constants(collision_references(m)) = 0.0;
        }            
    }
    
    discretization.SetCollisionConstants( collision_constants );
    discretization.SetCollisionLocations( collision_locations );
}
//#####################################################################


