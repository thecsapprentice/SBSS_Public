#include "config.h"
#include "surgicalActions.h"
#include "GraphicsUtils/wxGraphics.h"

#include <Thread_Queueing/PTHREAD_QUEUE.h>
#include <time.h>
#include <iostream>
#include <stdlib.h> 
#include <vector>

int main( int argc, char** argv ){

    
    using namespace PhysBAM;
    extern PTHREAD_QUEUE* pthread_queue;
    pthread_queue=new PhysBAM::PTHREAD_QUEUE(12);  
    surgicalActions* _surgAct;
    wxGraphics* _wGG;

    time_t time_start;
    double seconds;  
    time_start = time(NULL);
    srand (time(NULL));

    std::vector<std::string> config_files;
    config_files.push_back( "./SBSS.default.conf" );  // Default config
    config_files.push_back( "/etc/SBSS/SBSS.conf" );  // System config
    config_files.push_back( "./SBSS.conf" );          // User config (DON'T COMMIT THIS!!!!!)        

    Configuration* config;
    config = new Configuration();
    bool configured = false;
    try{
        config->parseFiles( config_files );
        configured = true;
    }
    catch( ConfigurationError e ){
        std::cout << e.what() << std::endl;
    }

    if( argc > 1 ){
        try{
            config->Parse( std::string(argv[1]) );
            configured = true;
        }
        catch( ConfigurationError e ){
            std::cout << e.what() << std::endl;           
        }
    }

    if( !configured )
        exit(1);

    while(true){
        _surgAct = new surgicalActions( config );
        _wGG = new wxGraphics();
        _surgAct->setWxGraphics(_wGG);

        int randomExample;
#if 0
        if( argc > 1 ){
            randomExample = atoi( argv[1] );
            switch( randomExample ){
            case 1:	
                _surgAct->pushCommand("{ \"command\":\"loadHistory\", \"data\":{ \"name\":\"./SurgerySim/data/zPlastyOneZ_UVW.hst\" }}");
                break;
            case 2:
                _surgAct->pushCommand("{ \"command\":\"loadHistory\", \"data\":{ \"name\":\"./SurgerySim/data/zPlastyTwoZ_UVW.hst\" }}");
                break;
            case 3:
                _surgAct->pushCommand("{ \"command\":\"loadHistory\", \"data\":{ \"name\":\"./SurgerySim/data/zPlasty3d85_UVW.hst\" }}");
                break;
            case 4:
                _surgAct->pushCommand("{ \"command\":\"loadHistory\", \"data\":{ \"name\":\"./SurgerySim/data/rhomboidFlatUVW.hst\" }}");
                break;
            case 5:
                _surgAct->pushCommand("{ \"command\":\"loadHistory\", \"data\":{ \"name\":\"./SurgerySim/data/sPlastyFlat_UVW.hst\" }}");
                break;
            case 6:
                _surgAct->pushCommand("{ \"command\":\"loadHistory\", \"data\":{ \"name\":\"./SurgerySim/data/scalpDualS_plastiesUVW.hst\" }}");
                break;
            case 7:
                _surgAct->pushCommand("{ \"command\":\"loadHistory\", \"data\":{ \"name\":\"./SurgerySim/data/scalpDufourmentelMoulyUVW.hst\" }}");
                break;        
            }            
        }        
#endif

        int counter = 0;
        bool c = true;
        while(c){
            _surgAct->update();
			// Enable for automatic exmaple playing
            #if 0
               seconds = difftime(time(NULL), time_start);
               if( seconds > 2){
                   _surgAct->pushCommand( "{ \"command\":\"historyNext\", \"data\":{ } }" );
                   time_start = time(NULL);
               }
            #endif
            counter++;
        }
        delete _surgAct;
        delete _wGG;
    }

        
    return 0;
}
