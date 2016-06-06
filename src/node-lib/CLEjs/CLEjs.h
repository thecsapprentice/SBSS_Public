#ifndef __CLEJS_H__
#define __CLEJS_H__

#include <EngineFrontend/CLElib.h>
#include <nan.h>

//using namespace PhysBAM;
void InitAll(v8::Local<v8::Object> exports);

class CLEjs : public Nan::ObjectWrap, public CLElib {
 public:
  static void Init(v8::Local<v8::Object> exports);
 
 private:

 explicit CLEjs();
  ~CLEjs();
  
  
 // Wrapper methods
     static void New(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static Nan::Persistent<v8::Function> constructor;

/*     static void __SetConstantMu(const Nan::FunctionCallbackInfo<v8::Value>& info); */
    
  
    // These parameters must be set before the call to Create_Model, or after Destroy_Model has been called

     /*!@brief Set the Youngs Modulus parameter
      *
      * @param youngs_modulus New value of the Youngs Modulus
      */     
     static void __Set_Youngs_Modulus(const Nan::FunctionCallbackInfo<v8::Value>& info);

     /*!@brief Set the Poissons Ratio parameter
      *
      * @param youngs_modulus New value of the Poissons Ratio
      */     
     static void __Set_Poissons_Ratio(const Nan::FunctionCallbackInfo<v8::Value>& info);


     // Model setup methods

     /*!@brief Initialize the embedded surface
      *
      * @param vertices Vertices of the emedded surface
      * @param triangles Triangles of the emedded surface
      * @param dx Discretization size of the emedding lattice
      * @param refine Extra refinement for the purposes of building hybrid lattices
      */  
     static void __Create_Model(const Nan::FunctionCallbackInfo<v8::Value>& info);

     static void __Add_Static_Model(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Add_Muscle_Layer(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Set_Collision_Model(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Set_Collision_Proxy(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Destroy_Model(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __SetTextureParameters(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __GetTextureStress(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Set_Fixed_Geometry(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Set_Fixed_Points(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Set_Fixed_Segments(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Set_Fixed_Triangles(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Set_Fixed_Volume(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Finalize_Initialization(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Set_Damping(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Set_Density(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Set_Time_Step(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Set_Hook_Stiffness(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Set_Suture_Stiffness(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Add_Hook(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Move_Hook(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Delete_Hook(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Add_Suture(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Delete_Suture(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Advance_One_Time_Step(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __WaitForSolve(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Update_Fine_Displacement(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Get_Vertices(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Get_Strain(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Get_Stress(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Get_Vertex_Data(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Update_Embedded_Surfaces(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Generate_UFine_MeshMap(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __UpdateDirichletCells(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Update_Collisions(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __WriteDebug(const Nan::FunctionCallbackInfo<v8::Value>& info);
     static void __Apply_Perturbation(const Nan::FunctionCallbackInfo<v8::Value>& info);

     
};

#endif


  
/*     enum{d=3}; */
/*     typedef float T; */
/*     typedef VECTOR<T,d> TV; */
/*     typedef VECTOR<int,d> T_INDEX; */
/*     typedef ELASTICITY_EXAMPLE<TV> BASE; */
/*     typedef typename REGULAR_HYPERCUBE_MESH<d>::ELEMENT T_ELEMENT; */

/*     // Members */

/*     using BASE::output_directory; */
/*     using BASE::stream_type;  */
/*     using BASE::parse_args; */
/*     using BASE::mesh; */
/*     using BASE::X; */
/*     using BASE::V; */
/*     using BASE::X_resting; */
/*     using BASE::constant_mu; */
/*     using BASE::constant_lambda; */
/*     using BASE::dx; */
/*     using BASE::cell_type; */
/*     using BASE::stiffness_matrix; */
/*     using BASE::subdomain_stiffness_matrix; */
/*     using BASE::subdomain_membership; */
/*     using BASE::node_flags; */
/*     using BASE::max_subdomain; */
/*     using BASE::node_hash; */
/*     using BASE::node_to_ijk; */
/*     using BASE::macroblock_size; */
/*     using BASE::macroblock_min_corner; */
/*     using BASE::frame_rate; */

/*     T_INDEX size; */


/*     // Methods */


    
/*     explicit ExampleBaseNode(); */
/*     ~ExampleBaseNode(); */
    
/* // Wrapper methods */
/*     static void New(const Nan::FunctionCallbackInfo<v8::Value>& info); */
/*     static Nan::Persistent<v8::Function> constructor;    */
/*     static void __SetConstantMu(const Nan::FunctionCallbackInfo<v8::Value>& info); */
/*     static void __SetConstantLambda(const Nan::FunctionCallbackInfo<v8::Value>& info); */
/*     static void __SetDx(const Nan::FunctionCallbackInfo<v8::Value>& info); */
/*     static void __SetSize(const Nan::FunctionCallbackInfo<v8::Value>& info); */
/*     static void __SetNewtonIterations(const Nan::FunctionCallbackInfo<v8::Value>& info); */
/*     static void __SetCGIterations(const Nan::FunctionCallbackInfo<v8::Value>& info); */
/*     static void __SetNewtonTolerance(const Nan::FunctionCallbackInfo<v8::Value>& info); */
/*     static void __SetCGTolerance(const Nan::FunctionCallbackInfo<v8::Value>& info); */
/*     static void __SetUseMacroblockSolver(const Nan::FunctionCallbackInfo<v8::Value>& info); */
        
/*     static void __Initialize(const Nan::FunctionCallbackInfo<v8::Value>& info); */

/*     static void __NodeCount( const Nan::FunctionCallbackInfo<v8::Value>& info); */
/*     static void __CellCount( const Nan::FunctionCallbackInfo<v8::Value>& info); */
/*     static void __TriangleCount( const Nan::FunctionCallbackInfo<v8::Value>& info); */

/*     static void __GetNodalData( const Nan::FunctionCallbackInfo<v8::Value>& info); */
/*     static void __GetCellData( const Nan::FunctionCallbackInfo<v8::Value>& info); */
/*     static void __GetTriangleData( const Nan::FunctionCallbackInfo<v8::Value>& info); */

/*     static void __Reset( const Nan::FunctionCallbackInfo<v8::Value>& info); */

    
/* // Implementation methods */

/*     virtual void Initialize(); */
/*     virtual void Set_Kinematic_Positions(const T time, const int frame, ARRAY_VIEW<TV> X) const; */
/*     virtual void Write_Output_Files(const int frame) const; */
/*     virtual void Set_Spring_Constraints(const T time, const int frame); */
/*     virtual void Set_Example_Other_Time_Parameters(const T time, const int frame);     */
/* }; */

/* #endif */
