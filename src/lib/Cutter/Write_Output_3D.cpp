#ifdef ENABLE_PHYSBAM_IO

#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>

#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>


#include "RANGE_ITERATOR.h"
#include "VOXELIZED_REGION_GENERATOR.h"


namespace PhysBAM{
//#####################################################################
// Function Write_Output (3D)
//#####################################################################

//template<class T>
    void Write_Output(STREAM_TYPE stream_type, VOXELIZED_REGION_GENERATOR<float,2>& vrg, 
                      const VOXELIZED_REGIONS<float,2>& regions,
                      const std::string directory, const int frame)
{
    PHYSBAM_FATAL_ERROR();
}
    void Write_Output(STREAM_TYPE stream_type, VOXELIZED_REGION_GENERATOR<double,2>& vrg, 
                      const VOXELIZED_REGIONS<double,2>& regions,
                      const std::string directory, const int frame)
{
    PHYSBAM_FATAL_ERROR();
}
    void Write_Output(STREAM_TYPE stream_type, VOXELIZED_REGION_GENERATOR<double,3>& vrg, 
                      const VOXELIZED_REGIONS<double,3>& regions,
                      const std::string directory, const int frame)
{
    PHYSBAM_FATAL_ERROR();
}

    void Write_Output(STREAM_TYPE stream_type, VOXELIZED_REGION_GENERATOR<float,3>& vrg, 
                      const VOXELIZED_REGIONS<float,3>& regions,
                      const std::string directory, const int frame)
{
	typedef float T;

    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;


    Initialize_Geometry_Particle();Initialize_Read_Write_Structures();
    GEOMETRY_PARTICLES<TV> particles;
    ARRAY<STRUCTURE<TV>* > display_structures;
    ARRAY<TRIANGULATED_SURFACE<T>*> bone_surfaces;
    ARRAY<TRIANGULATED_SURFACE<T>*> embedded_surfaces;

    
    RANGE<T_INDEX> coarse_grid_domain( vrg.grid.Cell_Indices() );
    RANGE<T_INDEX> fine_grid_domain( vrg.fine_grid.Cell_Indices() );

    RANGE<T_INDEX> coarse_grid_node_domain( vrg.grid.Node_Indices() );
    RANGE<T_INDEX> fine_grid_node_domain( vrg.fine_grid.Node_Indices() );

    int fine_grid_m=fine_grid_node_domain.Edge_Lengths()(1)+1;
    int fine_grid_n=fine_grid_node_domain.Edge_Lengths()(2)+1;
    int fine_grid_mn=fine_grid_node_domain.Edge_Lengths()(3)+1;

    //LOG::cout << "Coarse Grid Bounds  " << coarse_grid_m << "  " << 
    //    coarse_grid_n << "  " << coarse_grid_mn << std::endl;

    //LOG::cout << "Fine Grid Bounds  " << fine_grid_m << "  " << 
    //    fine_grid_n << "  " << fine_grid_mn << std::endl;

    // Initialize particles
    particles.array_collection->Add_Elements(fine_grid_m*fine_grid_n*fine_grid_mn);
    {int p=0;
    for(RANGE_ITERATOR<d> iterator(fine_grid_node_domain);iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        p++;

        //LOG::cout << index<< " --> " << p << std::endl;
        //if(vrg.voxmap(index)) for(int v=1;v<=d;v++) particles.X(p)(v)=vrg.fine_grid.Node(index)(v);
        //else particles.X(p)=TV();
        for(int v=1;v<=d;v++) particles.X(p)(v)=vrg.fine_grid.Node(index)(v);
    }}

    // Interior tetrahedralized volume
    ARRAY< TETRAHEDRALIZED_VOLUME<T>* > fine_cell_blocks;
    ARRAY< bool > fine_cell_blocks_active;
    for(RANGE_ITERATOR<d> coarse_iterator(coarse_grid_domain);coarse_iterator.Valid();coarse_iterator.Next()){
        const T_INDEX& coarse_index=coarse_iterator.Index();
        fine_cell_blocks.Append( TETRAHEDRALIZED_VOLUME<T>::Create(particles) );
        TETRAHEDRALIZED_VOLUME<T>& interior_tetrahedralized_volume = *fine_cell_blocks.Last();
        interior_tetrahedralized_volume.Update_Number_Nodes();


        
        fine_cell_blocks_active.Append( false );

        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>((coarse_index-1)*vrg.refinement_factor+1, (coarse_index-1)*vrg.refinement_factor+vrg.refinement_factor));iterator.Valid();iterator.Next()){
            
            const T_INDEX& index=iterator.Index();
            
            if(vrg.voxmap(index)){
                fine_cell_blocks_active.Last() = true;
                int i,j,ij;index.Get(i,j,ij);
                const int q=fine_grid_n;
                const int pq=fine_grid_mn;             
                const int LDB=(i-1)*q*pq+(j-1)*pq+(ij-1)+1;
                const int LDF=(i-1)*q*pq+(j-1)*pq+(ij-0)+1;
                const int LUB=(i-1)*q*pq+(j-0)*pq+(ij-1)+1;
                const int LUF=(i-1)*q*pq+(j-0)*pq+(ij-0)+1;
                const int RDB=(i-0)*q*pq+(j-1)*pq+(ij-1)+1;
                const int RDF=(i-0)*q*pq+(j-1)*pq+(ij-0)+1;
                const int RUB=(i-0)*q*pq+(j-0)*pq+(ij-1)+1;
                const int RUF=(i-0)*q*pq+(j-0)*pq+(ij-0)+1;
                

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
    }




    // Interior tetrahedralized volume
    ARRAY< TETRAHEDRALIZED_VOLUME<T>* > coarse_cell_blocks;
    ARRAY< bool > coarse_cell_blocks_active;
    for(RANGE_ITERATOR<d> coarse_iterator(coarse_grid_domain);coarse_iterator.Valid();coarse_iterator.Next()){
        const T_INDEX& coarse_index=coarse_iterator.Index();
        coarse_cell_blocks.Append( TETRAHEDRALIZED_VOLUME<T>::Create(particles) );
        TETRAHEDRALIZED_VOLUME<T>& interior_tetrahedralized_volume = *coarse_cell_blocks.Last();
        interior_tetrahedralized_volume.Update_Number_Nodes();


        coarse_cell_blocks_active.Append( false );

        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(coarse_index, coarse_index));iterator.Valid();iterator.Next()){
            
            const T_INDEX& index=iterator.Index();
            T_INDEX new_index = (index)*vrg.refinement_factor;

            if(fine_cell_blocks_active(coarse_cell_blocks_active.m) ){
                coarse_cell_blocks_active.Last() = true;
                int offset = vrg.refinement_factor;
                int i,j,ij;new_index.Get(i,j,ij);
                const int q=fine_grid_n;
                const int pq=fine_grid_mn;
                const int LDB=(i-offset)*q*pq+(j-offset)*pq+(ij-offset)+1;
                const int LDF=(i-offset)*q*pq+(j-offset)*pq+(ij-0)+1;
                const int LUB=(i-offset)*q*pq+(j-0)*pq+(ij-offset)+1;
                const int LUF=(i-offset)*q*pq+(j-0)*pq+(ij-0)+1;
                const int RDB=(i-0)*q*pq+(j-offset)*pq+(ij-offset)+1;
                const int RDF=(i-0)*q*pq+(j-offset)*pq+(ij-0)+1;
                const int RUB=(i-0)*q*pq+(j-0)*pq+(ij-offset)+1;
                const int RUF=(i-0)*q*pq+(j-0)*pq+(ij-0)+1;
                
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
    }






    // Interior tetrahedralized volume
    ARRAY< ARRAY< TETRAHEDRALIZED_VOLUME<T>*> , T_INDEX > sub_cell_blocks;
    ARRAY< ARRAY< bool > , T_INDEX > sub_cell_blocks_active;
    sub_cell_blocks.Resize( coarse_grid_domain );
    sub_cell_blocks_active.Resize( coarse_grid_domain );
    int max_sub_cell_depth = 0;   

    for(RANGE_ITERATOR<d> coarse_iterator(coarse_grid_domain);coarse_iterator.Valid();coarse_iterator.Next()){
        const T_INDEX& coarse_index=coarse_iterator.Index();

        int root_cell_index = vrg.root_grid_cell_mapping.Find( coarse_index );
        //LOG::cout << "Root cell " << root_cell_index << " has " << vrg.root_sub_mapping.Count_Matches( root_cell_index ) << " subcells" << std::endl;
        max_sub_cell_depth = max( max_sub_cell_depth, vrg.root_sub_mapping.Count_Matches( root_cell_index ) );
        
        int last_sub_cell = 0;
        for( ; vrg.root_sub_mapping.Find( root_cell_index, last_sub_cell+1, last_sub_cell); )
            {
                sub_cell_blocks(coarse_index).Append( TETRAHEDRALIZED_VOLUME<T>::Create(particles) );
                TETRAHEDRALIZED_VOLUME<T>& interior_tetrahedralized_volume=*sub_cell_blocks(coarse_index).Last();
                interior_tetrahedralized_volume.Update_Number_Nodes();
                sub_cell_blocks_active(coarse_index).Append(false);
                
                RANGE<T_INDEX> fine_domain((coarse_index-1)*vrg.refinement_factor+1,
                                           (coarse_index-1)*vrg.refinement_factor+vrg.refinement_factor);
                RANGE<T_INDEX> vox_domain(T_INDEX::All_Ones_Vector(),
                                          T_INDEX::All_Ones_Vector() * vrg.refinement_factor);

                for(RANGE_ITERATOR<d> iterator(fine_domain), viter(vox_domain);
                    iterator.Valid() && viter.Valid();
                    iterator.Next(), viter.Next()){
                    const T_INDEX& index=iterator.Index();
                    
                    if(vrg.subcell_voxmaps(last_sub_cell)(viter.Index())){
                            sub_cell_blocks_active(coarse_index).Last() = true;
                            int i,j,ij;index.Get(i,j,ij);
                            const int q=fine_grid_n;
                            const int pq=fine_grid_mn;             
                            const int LDB=(i-1)*q*pq+(j-1)*pq+(ij-1)+1;
                            const int LDF=(i-1)*q*pq+(j-1)*pq+(ij-0)+1;
                            const int LUB=(i-1)*q*pq+(j-0)*pq+(ij-1)+1;
                            const int LUF=(i-1)*q*pq+(j-0)*pq+(ij-0)+1;
                            const int RDB=(i-0)*q*pq+(j-1)*pq+(ij-1)+1;
                            const int RDF=(i-0)*q*pq+(j-1)*pq+(ij-0)+1;
                            const int RUB=(i-0)*q*pq+(j-0)*pq+(ij-1)+1;
                            const int RUF=(i-0)*q*pq+(j-0)*pq+(ij-0)+1;
                            
                            
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
            }
    }



    // Update particle count on all structures
    for(int s=1;s<=display_structures.m;s++)
        display_structures(s)->Update_Number_Nodes();


    // Select which structures to display
 
    DEFORMABLE_GEOMETRY_COLLECTION<TV> coarse_collection(particles);
    for( int i =1; i <= coarse_cell_blocks.m; i++)
        if(coarse_cell_blocks_active(i)){
            coarse_collection.Add_Structure(coarse_cell_blocks(i));
        }
    FILE_UTILITIES::Create_Directory(directory+"/"+STRING_UTILITIES::Value_To_String(0));
    coarse_collection.Write(stream_type,directory,0,0,true);



    DEFORMABLE_GEOMETRY_COLLECTION<TV> fine_collection(particles);
    for( int i =1; i <= fine_cell_blocks.m; i++)
        if(fine_cell_blocks_active(i)){
            fine_collection.Add_Structure(fine_cell_blocks(i));
        }
    FILE_UTILITIES::Create_Directory(directory+"/"+STRING_UTILITIES::Value_To_String(1));
    fine_collection.Write(stream_type,directory,1,1,true);

 

    for( int i = 1; i <= max_sub_cell_depth; i++)
        {
            int added_structs = 0;
            DEFORMABLE_GEOMETRY_COLLECTION<TV>* subcell_collection =
                new DEFORMABLE_GEOMETRY_COLLECTION<TV>(particles);

             for(RANGE_ITERATOR<d> coarse_iterator(coarse_grid_domain);
                coarse_iterator.Valid();coarse_iterator.Next()){
                const T_INDEX& coarse_index=coarse_iterator.Index();

                if( i <= sub_cell_blocks( coarse_index ).m ){
                    if(sub_cell_blocks_active(coarse_index)(i)){
                        subcell_collection->Add_Structure( sub_cell_blocks( coarse_index )(i) );
                        added_structs++;
                    }
                }
             }
             
             LOG::cout << "For subcell layer " << i << " writing out " << added_structs << " structures." << std::endl;
             FILE_UTILITIES::Create_Directory(directory+"_SubCells/"+STRING_UTILITIES::Value_To_String(i-1));
             subcell_collection->Write(stream_type,directory+"_SubCells",i-1,i-1,true);
             delete subcell_collection;
        }











    ARRAY<STRUCTURE<TV>* > region_display_structures;

    
    // Interior tetrahedralized volume
    ARRAY< TETRAHEDRALIZED_VOLUME<T>* > region_blocks;
    ARRAY< TETRAHEDRALIZED_VOLUME<T>* > fine_region_blocks;
    for( int region = 1; region <= regions.regions.m; region++ ){
        LOG::cout << "Constructing Region " << region << std::endl;
        region_blocks.Append( TETRAHEDRALIZED_VOLUME<T>::Create(particles) );
        TETRAHEDRALIZED_VOLUME<T>& interior_tetrahedralized_volume = *region_blocks.Last();
        interior_tetrahedralized_volume.Update_Number_Nodes();

        fine_region_blocks.Append( TETRAHEDRALIZED_VOLUME<T>::Create(particles) );
        TETRAHEDRALIZED_VOLUME<T>& fine_interior_tetrahedralized_volume = *fine_region_blocks.Last();
        fine_interior_tetrahedralized_volume.Update_Number_Nodes();

        for( int cell = 1; cell <=regions.grid_regions( region ).m; cell++){
            
            const T_INDEX& index = regions.grid_regions( region )(cell);
            T_INDEX new_index = (index)*vrg.refinement_factor;

            int offset = vrg.refinement_factor;
            int i,j,ij;new_index.Get(i,j,ij);
            const int q=fine_grid_n;
            const int pq=fine_grid_mn;
            const int LDB=(i-offset)*q*pq+(j-offset)*pq+(ij-offset)+1;
            const int LDF=(i-offset)*q*pq+(j-offset)*pq+(ij-0)+1;
            const int LUB=(i-offset)*q*pq+(j-0)*pq+(ij-offset)+1;
            const int LUF=(i-offset)*q*pq+(j-0)*pq+(ij-0)+1;
            const int RDB=(i-0)*q*pq+(j-offset)*pq+(ij-offset)+1;
            const int RDF=(i-0)*q*pq+(j-offset)*pq+(ij-0)+1;
            const int RUB=(i-0)*q*pq+(j-0)*pq+(ij-offset)+1;
            const int RUF=(i-0)*q*pq+(j-0)*pq+(ij-0)+1;
            
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
                interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,LDB,RUB,RDF));
            }

            RANGE<T_INDEX> fine_domain((index-1)*vrg.refinement_factor+1,
                                       (index-1)*vrg.refinement_factor+vrg.refinement_factor);
            RANGE<T_INDEX> vox_domain(T_INDEX::All_Ones_Vector(),
                                      T_INDEX::All_Ones_Vector() * vrg.refinement_factor);
            
            for(RANGE_ITERATOR<d> iterator(fine_domain), viter(vox_domain);
                    iterator.Valid() && viter.Valid();
                    iterator.Next(), viter.Next()){
            
            const T_INDEX& fine_index=iterator.Index();
            
            if(regions.voxmap_regions(region)(cell)(viter.Index())){
                int i,j,ij;fine_index.Get(i,j,ij);
                const int q=fine_grid_n;
                const int pq=fine_grid_mn;             
                const int LDB=(i-1)*q*pq+(j-1)*pq+(ij-1)+1;
                const int LDF=(i-1)*q*pq+(j-1)*pq+(ij-0)+1;
                const int LUB=(i-1)*q*pq+(j-0)*pq+(ij-1)+1;
                const int LUF=(i-1)*q*pq+(j-0)*pq+(ij-0)+1;
                const int RDB=(i-0)*q*pq+(j-1)*pq+(ij-1)+1;
                const int RDF=(i-0)*q*pq+(j-1)*pq+(ij-0)+1;
                const int RUB=(i-0)*q*pq+(j-0)*pq+(ij-1)+1;
                const int RUF=(i-0)*q*pq+(j-0)*pq+(ij-0)+1;
                

                if((i+j+ij)%2==0){
                    fine_interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LDB,RDB,LUB,LDF));
                    fine_interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(RDB,RDF,RUF,LDF));
                    fine_interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUB,RUB,RUF,RDB));
                    fine_interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,RUF,LDF,LUB));
                    fine_interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(RDB,LDF,RUF,LUB));}
                else{
                    fine_interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LDB,RDB,RUB,RDF));
                    fine_interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LDB,LUB,LUF,RUB));
                    fine_interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,RDF,LDF,LDB));
                    fine_interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,RUF,RDF,RUB));
                    fine_interior_tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,LDB,RUB,RDF));}}}
        }
       
        region_display_structures.Append(&fine_interior_tetrahedralized_volume);
        region_display_structures.Append(&interior_tetrahedralized_volume);
    }



// Update particle count on all structures
    for(int s=1;s<=region_display_structures.m;s++)
        region_display_structures(s)->Update_Number_Nodes();



    DEFORMABLE_GEOMETRY_COLLECTION<TV> region_collection(particles);
    for( int i =1; i <= region_blocks.m; i++){
        region_collection.Add_Structure(region_blocks(i));
    }

    FILE_UTILITIES::Create_Directory(directory+"_REGIONS/"+STRING_UTILITIES::Value_To_String(0));
    region_collection.Write(stream_type,directory+"_REGIONS",0,0,true);

    DEFORMABLE_GEOMETRY_COLLECTION<TV> fine_region_collection(particles);
    for( int i =1; i <= region_blocks.m; i++){
        fine_region_collection.Add_Structure(fine_region_blocks(i));
    }

    FILE_UTILITIES::Create_Directory(directory+"_REGIONS/"+STRING_UTILITIES::Value_To_String(1));
    fine_region_collection.Write(stream_type,directory+"_REGIONS",1,1,true);


    for( int i =1; i <= region_blocks.m; i++){
        DEFORMABLE_GEOMETRY_COLLECTION<TV> subregion_collection(particles);
        
        TETRAHEDRALIZED_VOLUME<T>* sub_region = new TETRAHEDRALIZED_VOLUME<T>((*region_blocks(i)).mesh, particles);
        TETRAHEDRALIZED_VOLUME<T>* sub_fine_region = new TETRAHEDRALIZED_VOLUME<T>((*fine_region_blocks(i)).mesh, particles);

        subregion_collection.Add_Structure(sub_region);
        subregion_collection.Add_Structure(sub_fine_region);
        FILE_UTILITIES::Create_Directory(directory+"_REGIONS/"+STRING_UTILITIES::Value_To_String(1+i));
        subregion_collection.Write(stream_type,directory+"_REGIONS",1+i,1+i,true);
    }


	
}
//#####################################################################
//	template<float> void Write_Output (STREAM_TYPE stream_type,EMBEDDED_DEFORMER& deformer,const std::string directory,const int frame);
//	template<double> void Write_Output(STREAM_TYPE stream_type,EMBEDDED_DEFORMER& deformer,const std::string directory,const int frame);
}

#endif
