#ifndef __GENERIC_CELL_H__
#define __GENERIC_CELL_H__

using namespace PhysBAM;

template<class T, int d>
struct GENERIC_CELL{
    typedef VECTOR<int,d> T_INDEX;
    bool is_grid;
    union{
        int _grid_index[d];
        int _mesh_index;
    };
    GENERIC_CELL() { is_grid=false; };
    GENERIC_CELL(const T_INDEX& g) { is_grid=true; grid_index() = g; };
    GENERIC_CELL(const int& m) { is_grid=false; mesh_index() = m; };
    GENERIC_CELL(const GENERIC_CELL& other){
        is_grid = other.is_grid;
        grid_index() = other.grid_index();
    };
    GENERIC_CELL& operator= (const GENERIC_CELL& other){
        if(&other != this){
            is_grid = other.is_grid;
            grid_index() = other.grid_index();
        }
        return *this;
    }    
    const T_INDEX& grid_index() const {return *reinterpret_cast<const T_INDEX*>(_grid_index);}
    T_INDEX& grid_index() {return *reinterpret_cast<T_INDEX*>(_grid_index);}
    const int& mesh_index() const {return _mesh_index;}
    int& mesh_index() {return _mesh_index;}
};

#endif
