#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <Thread_Queueing/PTHREAD_QUEUE.h>

#include "CORE_ENGINE_ACCELERATOR.h"
#include <EngineInterfaceMPI/COMMAND_ID.h>

#define COMMAND_TAG 1
#define DATA_TAG 2
#define DIE_TAG 3

using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;

void HandleCommand( int command, CORE_ENGINE_ACCELERATOR& coreEngine);

int main (int argc, char** argv)
{
    //freopen("/dev/null", "w", stdout);
  LOG::Initialize_Logging(true);
  FILE * pFile;
  pFile = fopen ("/tmp/log.phi" , "w");
  if (pFile == NULL){
      perror("Error: Can't Open Log File");
      return 0;
  }
  //LOG::Instance()->log_file = pFile;
  //LOG::Instance()->xml = false;

  int required=MPI_THREAD_SERIALIZED;
  int provided; 

  int rank, size;

  int threads = 1;
  if( argc > 1 ) 
      threads = atoi( argv[1] );

  pthread_queue=new PhysBAM::PTHREAD_QUEUE(threads); 

  MPI_Init_thread (&argc, &argv, required, &provided);	/* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &size);	/* get number of processes */

  if( provided < required ){
      threads = 1;
      omp_set_num_threads(1);
      LOG::cout << "Warning: This MPI version does not support threads. Downgrading to single threaded mode." << std::endl;
  }

/*
  if( rank != 1 ){
      LOG::cout << "Error: Accelerator can only run as Rank 1!" << std::endl;
      MPI_Finalize();
      return 0;
  }
 
*/

  LOG::cout << "Accelerator Starting..." << std::endl;

  CORE_ENGINE_ACCELERATOR coreEngine(threads);   

  LOG::cout << "Accelerator Online with " << threads << " threads."<< std::endl;
                
  typedef int command_t;
  command_t work;
  MPI_Status status;
  bool running = true;
  while(running){
      //LOG::cout << "Accelerator Waiting for Command." << std::endl;           
      fflush (pFile);
      MPI_Recv(&work, 1, MPI_INT, 0, MPI_ANY_TAG,
               MPI_COMM_WORLD, &status);

      switch(status.MPI_TAG){
      case COMMAND_TAG:
          {
              HandleCommand( work, coreEngine );
          }
          break;
      case DIE_TAG:
          {
              running = false;             
          }
          break;
      default:
          {
              LOG::cout << "Got unexpected input." << std::endl;
              LOG::cout << "Source: " << status.MPI_SOURCE;
              LOG::cout << "Tag: " << status.MPI_TAG;
              LOG::cout << "Data: " << work;
          }
          break;          
      }
  }

  fclose (pFile);
  MPI_Finalize();
  return 0;
}


void HandleCommand( int command, CORE_ENGINE_ACCELERATOR& coreEngine){
  switch( command ){
  case CREATE_ENGINE:
      coreEngine.CreateEngine();
      break;
  case INITIALIZE_MUSCLES:
      coreEngine.Initialize_Muscles();
      break;
  case INITIALIZE_MESH:
      coreEngine.Initialize_Mesh();
      break;
  case INITIALIZE_ENGINE:
      coreEngine.InitializeEngine();
      break;
  case EXACT_SOLVE:
      coreEngine.Exact_Solve();
      break;
  case DESTROY_ENGINE:
      coreEngine.DestroyEngine();
      break;
  case DISPLACEMENT_GRID:
      coreEngine.Displacement_Grid();
      break;
  case DISPLACEMENT_MESH:
      coreEngine.Displacement_Mesh();
      break;
  case H:
      coreEngine.h();
      break;
  case CELL:
      coreEngine.Cell();
      break;
  case NODE:
      coreEngine.Node();
      break;
  case PADDED_NODE_DOMAIN:
      coreEngine.Padded_Node_Domain();
      break;
  case UNPADDED_NODE_DOMAIN:
      coreEngine.Unpadded_Node_Domain();
      break;
  case PADDED_CELL_DOMAIN:
      coreEngine.Padded_Cell_Domain();
      break;
  case UNPADDED_CELL_DOMAIN:
      coreEngine.Unpadded_Cell_Domain();
      break;
  case GET_GRID:
      coreEngine.Grid();
      break;
  case GET_CELL_TYPE:
      coreEngine.GetCellType();
      break;
  case SET_CELL_TYPE:
      coreEngine.SetCellType();
      break;
  case GET_CELL_TYPE_MESH:
      coreEngine.GetCellTypeMesh();
      break;
  case SET_CELL_TYPE_MESH:
      coreEngine.SetCellTypeMesh();
      break;
  case GET_CELL_INDICES_MESH:
      coreEngine.GetCellIndicesMesh();
      break;
  case SET_CELL_INDICES_MESH:
      coreEngine.SetCellIndicesMesh();
      break;
  case GET_CELLS_MESH:
      coreEngine.GetCellsMesh();
      break;
  case SET_CELLS_MESH:
      coreEngine.SetCellsMesh();
      break;
  case NUMBER_OF_MESH_CELLS:
      coreEngine.Number_Of_Mesh_Cells();
      break;
  case NUMBER_OF_MESH_NODES:
      coreEngine.Number_Of_Mesh_Nodes();
      break;
  case CREATE_NEW_CONSTRAINT:
      coreEngine.CreateNewConstraint();
      break;
  case GET_CONSTRAINT:
      coreEngine.GetConstraint();
      break;
  case SET_CONSTRAINT:
      coreEngine.SetConstraint();
      break;
  case REMOVE_CONSTRAINT:
      coreEngine.RemoveConstraint();
      break;
  case NUMBER_OF_CONSTRAINTS:
      coreEngine.NumberOfConstraints();
      break;
  case GET_COLLISION_CONSTANTS:
      coreEngine.GetCollisionConstants();
      break;
  case SET_COLLISION_CONSTANTS:
      coreEngine.SetCollisionConstants();
      break;
  case GET_COLLISION_LOCATIONS:
      coreEngine.GetCollisionLocations();
      break;
  case SET_COLLISION_LOCATIONS:
      coreEngine.SetCollisionLocations();
      break;
  case ENGINE_READY:
      coreEngine.EngineReady();
      break;
  case ENGINE_INITIALIZED:
      coreEngine.EngineInitialized();
      break;
  case DISPLACEMENT:
      coreEngine.Displacement();
      break;
  case DEFORMATION:
      coreEngine.Deformation();
      break;
  case STRESS:
      coreEngine.Stress();
      break;
  case STRAIN:
      coreEngine.Strain();
      break;
  case GET_COARSEGRID:
      coreEngine.GetCoarseGrid();
      break;
  case GET_MUSCLE_DATA:
      coreEngine.GetMuscleData();
      break;
  case SET_MUSCLE_DATA:
      coreEngine.SetMuscleData();
      break;
  case GET_MUSCLE_ACTIVATIONS:
      coreEngine.GetMuscleActivations();
      break;
  case SET_MUSCLE_ACTIVATIONS:
      coreEngine.SetMuscleActivations();
      break;
  case GET_MUSCLE_MAX_STRESS:
      coreEngine.GetMuscleMaxStress();
      break;
  case SET_MUSCLE_MAX_STRESS:
      coreEngine.SetMuscleMaxStress();
      break;
  default:
      LOG::cout << "Invalid Command (" << command << ")." << std::endl;
      break;
  }
  MPI_Barrier( MPI_COMM_WORLD );
}
