//#####################################################################
// Copyright 2013, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VOXELIZED_REGION_GENERATOR
//#####################################################################

#ifndef __VOXELIZED_REGION_GENERATOR_H__
#define __VOXELIZED_REGION_GENERATOR_H__

#include "REGION_GENERATOR.h"


namespace PhysBAM {

    template< class T, int d >
        class VOXELIZED_REGIONS : public REGIONS<T,d>
        {
        public:
            explicit VOXELIZED_REGIONS( const REGIONS<T,d>& base ) : REGIONS<T,d>(base) {}

            typedef bool T_FLAG;
            typedef ARRAY<T_FLAG, typename REGIONS<T,d>::T_INDEX> T_FLAG_ARRAY;
            typedef ARRAY< T_FLAG_ARRAY > T_VOXMAP_REGION;
            typedef ARRAY< T_VOXMAP_REGION > T_VOXMAP_REGIONS;

            T_VOXMAP_REGIONS voxmap_regions;

        };
    

    template< class T, int d >
        class VOXELIZED_REGION_GENERATOR : public REGION_GENERATOR<T,d>
        {
            //private:
        public:
            typedef REGION_GENERATOR<T,d> BASE;

            typedef typename BASE::T_INDEX T_INDEX;
            typedef typename BASE::TV TV;
            typedef typename BASE::T_CELL T_CELL;
            typedef typename BASE::T_FACE T_FACE;
            typedef typename BASE::T_GRID T_GRID;

            typedef bool T_FLAG;
            typedef ARRAY<T_FLAG,T_INDEX> T_FLAG_ARRAY;

            const int refinement_factor;
            const T_GRID fine_grid;
            const T_FLAG_ARRAY voxmap;
            ARRAY< T_FLAG_ARRAY > subcell_voxmaps;

            using BASE::root_cells;
            using BASE::root_grid_cell_mapping;
            using BASE::sub_cells;
            using BASE::root_sub_mapping;
            using BASE::rootnode_to_subnode;
            using BASE::root_cell_to_subcell;
            using BASE::region_sets;

            VOXELIZED_REGIONS<T,d> vregions;

        public: 
            

            VOXELIZED_REGION_GENERATOR( int refinement_factor_input,
                                        const T_GRID& fine_grid_input,
                                        const T_FLAG_ARRAY& voxmap_input );

            virtual void Generate();
            
            virtual const VOXELIZED_REGIONS<T,d>* GetRegionData() const;
            
            const T_FLAG_ARRAY& GetSubcellVoxmap( const T_INDEX& cell_index, int duplicate ) const;


            virtual VECTOR<bool,T_CELLTYPE<T, d>::vertices_per_cell> CornersInsideOutside( const T_INDEX& cell_index, const int subcell ) const;

        private:
            
            virtual void ResolveCuts();
            virtual bool IsMaterialContinous(const int root_cellA, const int sub_cellA, const AXIS faceA,
                                             const int root_cellB, const int sub_cellB, const AXIS faceB) const; 



        };


}

#endif
