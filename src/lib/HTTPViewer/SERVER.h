#ifndef __HTTPVIEWER_SERVER_H__
#define __HTTPVIEWER_SERVER_H__

#include <vector>
#include <string>
#include <map>
// #include <functional>  // cbc added to pass callback function
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>

struct mg_server;
struct mg_connection;

extern int curl_initialized;

namespace HTTPVIEWER{

    struct CLIENT_DATA;

    struct CONNECTION_DATA{
        std::string ip;
        int port;
        int revision;
        static std::string buildId( std::string _ip, int _port ){
            std::string id = _ip;
            id += std::to_string( (long long) _port );
            return id;
        }
    };

    class SERVER;

    class WORKER{
    private:
        mg_server* server_handle;
        CLIENT_DATA* data;
        SERVER* master;
        
        boost::mutex* server_lock;
        boost::mutex* data_lock;

        std::map<std::string,CONNECTION_DATA> connections;

    private:
        boost::thread m_Thread;
        bool running;

    public:
        WORKER(boost::mutex*, boost::mutex*);
        void run();
        void SetHandle(mg_server*);
        void SetData(CLIENT_DATA*);
        void SetMaster(SERVER*);
        void start();
        void stop();

        friend int handle_connection( struct mg_connection* conn);
        friend int iterate_callback( struct mg_connection* conn);
        friend int lost_connection( struct mg_connection* conn);
    };

    class SERVER{

    public:
        SERVER();
        SERVER(int port, std::string token, std::string url);
        ~SERVER();
        
        void CreateServer();
        void DestroyServer();

        void StartServer();
        void StopServer();

        std::string getClientCommand(); // gets client command from queue
        void clearClientCommand();  // clears queue and ends client and server wait state

        void UpdateTopology( const std::vector<int>& t);
        void UpdateTextureCoords( const std::vector<double>& t);
        void UpdateVertices( const std::vector<double>& t);
        void UpdateVertexData( const std::vector<double>& v, const std::vector<double>& s1, const std::vector<double>& s2, const std::vector<float>& datatex);

        void UpdateNormals( const std::vector<double>& t);
        void UpdateDynamicData( const std::vector<int>& t, 
                                const std::vector<double>& v, 
                                const std::vector<double>& uv);            

        void SetTextureName( const std::string& t);
        void SetNormalName( const std::string& t);

        void AddStaticObject( const std::vector<int>& topo,
                              const std::vector<double>& verts,
                              const std::vector<double>& normals,
                              const std::vector<double>& uvs,
                              const std::string& texturename);  // obsolete - nuke later

        void AddStaticObject( const std::vector<int>& topo,
                              const std::vector<double>& verts,
                              const std::vector<double>& uvs,
                              const std::string& texturename,
                              const std::string& normalname);

        void RegisterTexture( const std::string& name, const char* data, int length);
        void RegisterTexture( const std::string& name, const std::string resource);
        
        // These routines are now depricated
        void AddHook( int triangle, std::vector<double> pos);
        void MoveHook( int hook, std::vector<double> pos);
        void DeleteHook( int hook );       

        void UpdateHooks( const std::vector< std::pair< int, std::vector<double> > >& );
        
        // These routines are now depricated
        void AddSuture( int triangleA, int edgeA, double UVA, int triangleB, int edgeB, double UVB );
        void DeleteSuture( int suture );

        void UpdateSutures( const std::vector< std::pair< std::vector<int>, std::vector<double> > >& );

        void SetSimulatorState( const int state );

        int getCurrentRevision();


        std::string GetData( std::string resource, std::string resource_type );
        void PostStatus( std::string message );
        std::string MessageSend( std::string url_additional, std::string message );

    protected:
        int port;
        std::string token;
        std::string url;

        mg_server* server_handle;
        CLIENT_DATA* data;
        WORKER* worker;   

        boost::mutex server_lock;
        boost::mutex data_lock;
    };

    int handle_connection( struct mg_connection* conn);
    int iterate_callback( struct mg_connection* conn);
    int lost_connection( struct mg_connection* conn);

}

#endif
