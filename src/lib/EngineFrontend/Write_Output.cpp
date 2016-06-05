
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

#include <Common/RANGE_ITERATOR.h>

#include "Write_Output.h"

namespace PhysBAM{

    void Create_Output_Data_Embedded(const ELASTIC_LATTICE_DEFORMER& deformer,
                                     GEOMETRY_PARTICLES< VECTOR<float,3> >& particles,
                                     ARRAY<STRUCTURE<VECTOR<float,3> >* >& collection_structures)
    {


        LOG::SCOPE scope("Output Structures : Embedded");

        typedef float T;
        static const int d=3;

        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
   
        
        // Add embedded surfaces

        {
            TRIANGULATED_SURFACE<T>* display_embedded_surface=TRIANGULATED_SURFACE<T>::Create(particles);
        
            int particle_base = particles.X.m;
        


            pthread_mutex_lock(&deformer.geometry_buffer_lock);

            const ARRAY<TV>* vertices=deformer.display_vertices;

            for( int v=1;v<=vertices->m;v++){
                int p = particles.array_collection->Add_Element();
                for(int w=1;w<=d;w++)
                    particles.X(p)(w)=vertices->operator()(v)(w);
            }
            for( int t=1;t<=deformer.triangles.m;t++){
                display_embedded_surface->mesh.elements.Append(deformer.triangles(t)+particle_base);
            }
        
            pthread_mutex_unlock(&deformer.geometry_buffer_lock);



            collection_structures.Append(display_embedded_surface);
        
            // Add static models
        
            LOG::cout << "Adding " << deformer.static_model_vertices.m << " static models." << std::endl;
        
            for( int M = 1; M <= deformer.static_model_vertices.m; M++)
                {
                    LOG::cout << "Static Model (" << M << ") has " << deformer.static_model_vertices(M).m << " vertices and " << deformer.static_model_triangles(M).m << " triangles." << std::endl;
                
                    TRIANGULATED_SURFACE<T>* display_embedded_surface=TRIANGULATED_SURFACE<T>::Create(particles);
                
                    int particle_base = particles.X.m;
                
                    for( int v=1;v<=deformer.static_model_vertices(M).m;v++){
                        int p = particles.array_collection->Add_Element();
                        particles.X(p)=deformer.static_model_vertices(M)(v);
                        //LOG::cout << particles.X(p) << std::endl;
                
                    }
                    for( int t=1;t<=deformer.static_model_triangles(M).m;t++){
                        display_embedded_surface->mesh.elements.Append(deformer.static_model_triangles(M)(t)+particle_base);
                    }
                
                    collection_structures.Append(display_embedded_surface);
                }
        }

    
        
        
    }


   void Create_Output_Data_FineGrid(const ELASTIC_LATTICE_DEFORMER& deformer,
                                     GEOMETRY_PARTICLES< VECTOR<float,3> >& particles,
                                     ARRAY<STRUCTURE<VECTOR<float,3> >* >& collection_structures)
    {
        typedef float T;
        static const int d=3;

        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;


        // Add Cutting fine grid mesh
   
        // Initialize particles
        int fine_grid_m=deformer.fine_grid.Node_Indices().Edge_Lengths()(1)+1;
        int fine_grid_n=deformer.fine_grid.Node_Indices().Edge_Lengths()(2)+1;
        int fine_grid_mn=deformer.fine_grid.Node_Indices().Edge_Lengths()(3)+1;
        
        int particle_base = particles.X.m;
        particles.array_collection->Add_Elements(fine_grid_m*fine_grid_n*fine_grid_mn);
        {int p=particle_base;
            for(RANGE_ITERATOR<d> iterator(deformer.unpadded_fine_node_domain);iterator.Valid();iterator.Next()){
                const T_INDEX& index=iterator.Index();
                p++;
                //if(deformer.voxmap_node(index)==1)
                for(int v=1;v<=d;v++)
                    particles.X(p)(v)=deformer.u_fine(v)(index)+deformer.fine_grid.Node(index)(v);
                //else
                //particles.X(p)=TV();
            }}
        
        // Interior tetrahedralized volume
        TETRAHEDRALIZED_VOLUME<T>& fine_tetrahedralized_volume=*TETRAHEDRALIZED_VOLUME<T>::Create(particles);
        fine_tetrahedralized_volume.Update_Number_Nodes();
        const ARRAY<T_INDEX>& embedding_map=deformer.embedding_map;
        
        //for( int v=1;v<=deformer.vertices.m;v++){
        //const T_INDEX& index = embedding_map(v);
        for(RANGE_ITERATOR<d> iterator(deformer.fine_grid.Cell_Indices());iterator.Valid();iterator.Next()){
            const T_INDEX& index=iterator.Index();
            
            if(deformer.voxmap( index )){
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
        collection_structures.Append(&fine_tetrahedralized_volume);
        

    }

//#####################################################################
//
//                      Write_Output_Minimal
//
//#####################################################################




    void Write_Output_Embedded(STREAM_TYPE stream_type, const ELASTIC_LATTICE_DEFORMER& deformer, const std::string directory,const int frame)
    {
        LOG::SCOPE scope( "Write Output Embedded");
        typedef float T;
        static const int d=3;
        
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        
        Initialize_Geometry_Particle();Initialize_Read_Write_Structures();    
        GEOMETRY_PARTICLES<TV> particles;
        ARRAY<STRUCTURE<TV>* > collection_structures;
        deformer.ConstDiscretization().Output_Structures(particles, collection_structures);
        Create_Output_Data_Embedded(deformer, particles, collection_structures);
        //Create_Output_Data_FineGrid(deformer, particles, collection_structures);

        for(int s=1;s<=collection_structures.m;s++)
            collection_structures(s)->Update_Number_Nodes();

        Write_Output_Structures(stream_type, particles, collection_structures, directory, frame);
    }

}
