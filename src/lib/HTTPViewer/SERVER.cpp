#define NO_MONGOOSE

#ifndef NO_MONGOOSE
#include <mongoose/mongoose.h>
#include <boost/filesystem.hpp>
#include <boost/system/error_code.hpp>
#include <boost/filesystem/fstream.hpp>

namespace fs = boost::filesystem;
namespace sys = boost::system;
#endif

#include "SERVER.h"
#include "DATA.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <curl/curl.h>

int curl_initialized = 0;

using namespace HTTPVIEWER;

int HTTPVIEWER::lost_connection(struct mg_connection *conn) {
#ifndef NO_MONGOOSE
    WORKER* worker = static_cast<WORKER*>(conn->server_param);

    if (conn->is_websocket) {
        std::cout << "Lost websocket connection from " << conn->remote_ip << ":" << conn->remote_port << std::endl;
        std::string id = CONNECTION_DATA::buildId( conn->remote_ip, conn->remote_port );
        std::map<std::string,CONNECTION_DATA>::iterator it;
        it = worker->connections.find( id );
        if( it == worker->connections.end() ){
            std::cout << "This is a new connection from " << id << ". Why don't we know about it?" << std::endl;
        }
        else{
            worker->connections.erase( it );
            std::cout << "We know this connection already. Removing it from list." << std::endl;
        }
    }
#endif
    return 0;
}

int HTTPVIEWER::iterate_callback(struct mg_connection *conn) {
#ifndef NO_MONGOOSE
    WORKER* worker = static_cast<WORKER*>(conn->server_param);

    if (conn->is_websocket) {
        std::string id = CONNECTION_DATA::buildId( conn->remote_ip, conn->remote_port );
        std::map<std::string,CONNECTION_DATA>::iterator it;
        it = worker->connections.find( id );
        if( it != worker->connections.end() ){
            int sinceTime = it->second.revision;
            int timestamp = 0;
            bool skip = false;
            std::string json;
            worker->data_lock->lock();
            //std::cout << "Client is at " << sinceTime << " and current state is :" << worker->data->global_timestamp << std::endl;
            timestamp = worker->data->global_timestamp;
            // Nathan - note added wait condition. Don't send until data fully composed.
            if(worker->data->global_timestamp <= sinceTime || !worker->data->command_queue.empty())
                skip = true;
            else{
                json = json_spirit::write_string(json_spirit::mValue(worker->data->getJson(sinceTime)));            
                sinceTime = worker->data->global_timestamp;
            }
            worker->data_lock->unlock();
            if(!skip){
                int res = mg_websocket_write(conn, 1, json.c_str(), json.length());
                //std::cout << "Writing to websocket: " << res << " bytes, since " << it->second.revision << ", timestamp: "<< timestamp << std::endl;
                it->second.revision = sinceTime;
            }
        }
    }
    return MG_REQUEST_PROCESSED;
#else
    return 0;
#endif
}

int HTTPVIEWER::handle_connection( struct mg_connection* conn){
#ifndef NO_MONGOOSE

    WORKER* worker = static_cast<WORKER*>(conn->server_param);

    std::cout << "Got request for uri: [" << conn->uri << "]"<< std::endl;
    
    if (conn->is_websocket) {
        std::cout << "Handling incoming websocket request." << std::endl;
        std::string id = CONNECTION_DATA::buildId( conn->remote_ip, conn->remote_port );
        if( worker->connections.find( id ) == worker->connections.end() ){
            std::cout << "This is a new connection from " << id << ". Adding it to list." << std::endl;
            CONNECTION_DATA cd;
            cd.ip = conn->remote_ip;
            cd.port = conn->remote_port;
            cd.revision = 0;
            worker->connections.insert( std::pair<std::string, CONNECTION_DATA>(id, cd ) );            
        }
        else{
            std::cout << "We know this connection already. Its revision is: "<< worker->connections[id].revision << std::endl;
        }

        if(conn->content_len == 4 && !memcmp(conn->content, "exit", 4) ){
            std::cout << "Client is going down!" << std::endl;
            worker->connections.erase(worker->connections.find( id ));
            return MG_CLIENT_CLOSE;
        }

        //mg_websocket_write(conn, 1, conn->content, conn->content_len);
        if( conn->content_len > 0 ){
            char* buf = new char[conn->content_len+1];
            memcpy(buf, conn->content, conn->content_len);
            buf[conn->content_len] = '\0';
            worker->data_lock->lock();
            std::cout << "Received message from client: " << std::string(buf) << std::endl;
            worker->data->command_queue.push_back(std::string(buf));
            //worker->data->timestamp_command = worker->data->Timestamp();
            worker->data_lock->unlock();
            delete [] buf;
        }
        int opcode = 0x0f & conn->wsbits;
        std::cout << "Client opcode is: " << opcode << std::endl;
        // This handler is called for each incoming websocket frame, one or more
        // times for connection lifetime.
        // Echo websocket data back to the client.
              
        return MG_CLIENT_CONTINUE;
    }
    else{
        mg_send_header(conn, "Access-Control-Allow-Origin", "*");    
    
        int sinceTime;
        char sinceTimeBuf[10];
        int res = mg_get_var( conn, "since", sinceTimeBuf, 10 );
        if( res > 0 )
            sinceTime = atoi( sinceTimeBuf );
        else
            sinceTime = -1;
        
        std::string uri( conn->uri );
        
        if( uri == "/data" ){
            worker->data_lock->lock();
            std::string json = json_spirit::write_string(json_spirit::mValue(worker->data->getJson(sinceTime)));
            mg_printf_data( conn, "%s", json.c_str() );
            worker->data_lock->unlock();
        }
        else{
            std::string rootdir("./http");
            fs::path filepath(rootdir + std::string(uri) );
            try{
                fs::path abs_root = fs::canonical( fs::path(rootdir) );
                fs::path abs_filepath = fs::canonical( filepath );
                std::cout << uri << std::endl;
                std::cout << filepath << std::endl;
                std::cout << abs_root << std::endl;
                std::cout << abs_filepath << std::endl;
                
                if( abs_filepath.native().find( abs_root.native() ) == 0 ){ //Ensure we belong in the root
                    if( fs::is_regular_file( abs_filepath ) ){
                        fs::ifstream fin(abs_filepath);
                        fin.seekg (0, fin.end);
                        int length = fin.tellg();
                        fin.seekg (0, fin.beg);
                        char* finData = new char[length];
                        fin.read(finData,length);
                        mg_send_data(conn, finData, length);
                        delete [] finData;
                        fin.close();
                    }
                    else
                        throw fs::filesystem_error(std::string("Can't serve directories"), sys::error_code());                
                }
                else
                    throw fs::filesystem_error(std::string("Permission Denied."), sys::error_code());                
                
            }
            catch( const fs::filesystem_error& ex ){
                std::cout << "[Error] " << ex.what() << std::endl;
                mg_send_status(conn, 404 );
                mg_printf_data( conn, "Error!" );
            }
        }
        
        return MG_REQUEST_PROCESSED;
    }
#else
    return 0;
#endif
}


WORKER::WORKER(boost::mutex* s_lock, boost::mutex* d_lock) :
    data(NULL), server_handle(NULL), server_lock(s_lock), data_lock(d_lock){
}

void WORKER::SetHandle( mg_server* __handle){
    server_handle = __handle;
};

void WORKER::SetData( CLIENT_DATA* __data){
    data = __data;
}

void WORKER::SetMaster( SERVER* __master){
    master = __master;
}

void WORKER::start() {
    running = true;
    m_Thread = boost::thread( &WORKER::run, this );
}

void WORKER::stop() {
    running = false;
    m_Thread.join();
}

void WORKER::run() {
    #ifndef NO_MONGOOSE

    unsigned int current_timer = 0, start_timer = 0;
    bool update_completed = false;
    int last_connections = 0;
    while(running){
        server_lock->lock();
        current_timer = mg_poll_server( server_handle, 0.5);        
        if( start_timer==0 )
            start_timer = current_timer;

        mg_iterate_over_connections(server_handle, iterate_callback, &current_timer);

        if(((current_timer - start_timer) % 3 == 0 && !update_completed) || last_connections != connections.size()){
            std::stringstream ss;
            ss << "{\"clients\":" << connections.size() << "}";
            master->PostStatus( ss.str() );
            update_completed = true;
        }

        if((current_timer - start_timer) % 3 != 0){
            update_completed = false;
        }       

        last_connections = connections.size();
        server_lock->unlock();
    }
    #endif
}



SERVER::SERVER() : port(5566), data(NULL), server_handle(NULL), worker(NULL){
    if( !curl_initialized )
        curl_global_init(CURL_GLOBAL_ALL);
    curl_initialized = 1;
}

SERVER::SERVER(int port, std::string token, std::string url) : port(port), token(token), url(url), data(NULL), server_handle(NULL), worker(NULL){
    if( !curl_initialized )
        curl_global_init(CURL_GLOBAL_ALL);
    curl_initialized = 1;
}


SERVER::~SERVER() {
#ifndef NO_MONGOOSE
    if(data){
        delete data;
        data = NULL;
    }
    if(worker){
        worker->stop();
        delete worker;
        worker = NULL;
    }
    if(server_handle)
        mg_destroy_server( &server_handle );
#endif
}

void SERVER::CreateServer(){
#ifndef NO_MONGOOSE
    std::cout << "Creating Server" << std::endl;

    server_lock.lock();
    if(!data)
        data = new CLIENT_DATA;
    if(!worker)
        worker = new WORKER(&server_lock, &data_lock);

    server_handle = mg_create_server(worker);
    mg_set_option( server_handle, "document_root", ".");
    char _port[16];
    sprintf(_port, "%d", port);
    mg_set_option( server_handle, "listening_port", _port);
    mg_set_request_handler(server_handle, handle_connection);
    mg_set_close_handler(server_handle, lost_connection);
    server_lock.unlock();

    worker->SetHandle( server_handle );
    worker->SetData( data );
    worker->SetMaster( this );
#endif
}


void SERVER::DestroyServer(){
#ifndef NO_MONGOOSE

    std::cout << "Destroying Server" << std::endl;

    server_lock.lock();
    if(worker){
        worker->stop();
        delete worker;
        worker = NULL;
    }
    if(data){
        delete data;
        data = NULL;
    }
    if(server_handle){
        mg_destroy_server( &server_handle );
        server_handle = NULL;
    }
    server_lock.unlock();
#endif

}



void SERVER::StartServer(){
    std::cout << "Starting Server"  << std::endl;
    worker->start();
}

void SERVER::StopServer(){
    std::cout << "Stopping Server"  << std::endl;
    worker->stop();
}

std::string SERVER::getClientCommand(){
    std::string c_command;
    data_lock.lock();
    if(!data->command_queue.empty()){
        c_command = data->command_queue.front();
        data->command_queue.pop_front();
    }
    else
        c_command = "";
    data_lock.unlock();
    return c_command;
}
        
void SERVER::clearClientCommand(){
    data_lock.lock();
    while(!data->command_queue.empty())
        data->command_queue.pop_front();
    data_lock.unlock();
}
        
void SERVER::UpdateTopology( const std::vector<int>& t){
    data_lock.lock();
    data->dynamic_object.topology.assign( t.begin(), t.end() );
    data->dynamic_object.timestamp_topology = data->Timestamp();
    data_lock.unlock();
}

void SERVER::UpdateTextureCoords( const std::vector<double>& t){
    data_lock.lock();
    data->dynamic_object.uvs.assign( t.begin(), t.end() );
    data->dynamic_object.timestamp_uvs = data->Timestamp();
    data_lock.unlock();
}

void SERVER::UpdateVertices( const std::vector<double>& t){
    data_lock.lock();
    data->dynamic_object.vertices.assign( t.begin(), t.end() );
    data->dynamic_object.strain = std::vector<double>( t.size(), 0.0 );    
    data->dynamic_object.timestamp_vertices = data->Timestamp();
    data_lock.unlock();
}

void SERVER::UpdateVertexData( const std::vector<double>& v, const std::vector<double>& s1, const std::vector<double>& s2, const std::vector<float>& datatex ){
    data_lock.lock();
    data->dynamic_object.vertices.assign( v.begin(), v.end() );
    data->dynamic_object.stress.assign( s1.begin(), s1.end() );
    data->dynamic_object.strain.assign( s2.begin(), s2.end() );
    data->dynamic_object.data_texture.assign( datatex.begin(), datatex.end() );
    data->dynamic_object.timestamp_vertices = data->Timestamp();
    data_lock.unlock();
}

void SERVER::UpdateDynamicData( const std::vector<int>& t, 
                                const std::vector<double>& v, 
                                const std::vector<double>& uv)
{
    data_lock.lock();
    int timestamp = data->Timestamp();
    data->dynamic_object.topology.assign( t.begin(), t.end() );
    data->dynamic_object.timestamp_topology = timestamp;
    data->dynamic_object.vertices.assign( v.begin(), v.end() );
    data->dynamic_object.strain = std::vector<double>( v.size(), 0.0 );    
    data->dynamic_object.timestamp_vertices = timestamp;
    data->dynamic_object.uvs.assign( uv.begin(), uv.end() );
    data->dynamic_object.timestamp_uvs = timestamp;
    data_lock.unlock();
}

void SERVER::UpdateNormals( const std::vector<double>& t){
    data_lock.lock();
    data->dynamic_object.normals.assign( t.begin(), t.end() );
    data->dynamic_object.timestamp_normals = data->Timestamp();
    data_lock.unlock();
}

void SERVER::SetTextureName( const std::string& t){
    data_lock.lock();
    data->dynamic_object.texturename = t;
    data->dynamic_object.timestamp_texturename = data->Timestamp();
    data_lock.unlock();
}

void SERVER::SetNormalName( const std::string& t){
    data_lock.lock();
    data->dynamic_object.normalname = t;
    data->dynamic_object.timestamp_normalname = data->Timestamp();
    data_lock.unlock();
}

void SERVER::AddStaticObject( const std::vector<int>& topo,
                              const std::vector<double>& verts, 
                              const std::vector<double>& normals,
                              const std::vector<double>& uvs, 
                              const std::string& texturename) {
    data_lock.lock();  // obsolete - nuke later
    int timestamp = data->Timestamp();
    OBJECT o;
    o.topology.assign( topo.begin(), topo.end() );
    o.vertices.assign( verts.begin(), verts.end() );
    o.normals.assign( normals.begin(), normals.end() );
    o.uvs.assign( uvs.begin(), uvs.end() );
    o.texturename = texturename;   
    o.timestamp_topology = timestamp;
    o.timestamp_uvs = timestamp;
    o.timestamp_vertices = timestamp;
    o.timestamp_normals = timestamp;
    o.timestamp_texturename = timestamp;
    data->static_objects.push_back( o );
    data_lock.unlock();
}

void SERVER::AddStaticObject( const std::vector<int>& topo,
                              const std::vector<double>& verts, 
                              const std::vector<double>& uvs, 
                              const std::string& texturename,
                              const std::string& normalname){
    data_lock.lock();
    int timestamp = data->Timestamp();
    OBJECT o;
    o.topology.assign( topo.begin(), topo.end() );
    o.vertices.assign( verts.begin(), verts.end() );
    o.uvs.assign( uvs.begin(), uvs.end() );
    o.texturename = texturename;   
    o.normalname = normalname;   
    o.timestamp_topology = timestamp;
    o.timestamp_uvs = timestamp;
    o.timestamp_vertices = timestamp;
    o.timestamp_texturename = timestamp;
    o.timestamp_normalname = timestamp;
    data->static_objects.push_back( o );
    data_lock.unlock();
}

// BASE64 ENCODING
static void base64_encode(const unsigned char *src, int src_len, char *dst) {
  static const char *b64 =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  int i, j, a, b, c;

  for (i = j = 0; i < src_len; i += 3) {
    a = src[i];
    b = i + 1 >= src_len ? 0 : src[i + 1];
    c = i + 2 >= src_len ? 0 : src[i + 2];

    dst[j++] = b64[a >> 2];
    dst[j++] = b64[((a & 3) << 4) | (b >> 4)];
    if (i + 1 < src_len) {
      dst[j++] = b64[(b & 15) << 2 | (c >> 6)];
    }
    if (i + 2 < src_len) {
      dst[j++] = b64[c & 63];
    }
  }
  while (j % 4 != 0) {
    dst[j++] = '=';
  }
  dst[j++] = '\0';
}

void SERVER::RegisterTexture( const std::string& name, const char* data, int length){
    std::cout << "WARNING!!! SERVER::RegisterTexture is a depricated function. Use RegisterTexture instead." << std::endl;
}

void SERVER::RegisterTexture( const std::string& name, const std::string resource){
    data_lock.lock();
    TEXTURE tex;
    tex.texture_data = resource;
    data->textures.insert( std::pair<std::string, TEXTURE>( name, tex ));
    data->timestamp_texturedata = data->Timestamp();
    data_lock.unlock();
}


void SERVER::AddHook( int triangle, std::vector<double> pos){
    std::cout << "WARNING!!! SERVER::AddHook is a depricated function. Use UpdateHooks instead." << std::endl;
}

void SERVER::MoveHook( int hook, std::vector<double> pos){
    std::cout << "WARNING!!! SERVER::MoveHook is a depricated function. Use UpdateHooks instead." << std::endl;
}

void SERVER::DeleteHook( int hook ){
    std::cout << "WARNING!!! SERVER::DeleteHook is a depricated function. Use UpdateHooks instead." << std::endl;
}


void SERVER::UpdateHooks( const std::vector< std::pair<int, std::vector<double> > >& new_hooks ){
    typedef std::vector< std::pair<int, std::vector<double> > > hook_desc;
    data_lock.lock();
    data->hooks.clear();
    for( hook_desc::const_iterator it=new_hooks.begin(); it != new_hooks.end(); it++){
        HOOK newhook;
        newhook.triangle_num = it->first;
        newhook.pos.assign( it->second.begin(), it->second.end() );
        data->hooks.push_back( newhook );
    }
    data->timestamp_hook = data->Timestamp();
    data_lock.unlock();
}


void SERVER::AddSuture( int triangleA, int edgeA, double UVA, int triangleB, int edgeB, double UVB ){
    std::cout << "WARNING!!! SERVER::AddSuture is a depricated function. Use UpdateSutures instead." << std::endl;
}

void SERVER::DeleteSuture( int suture ){
    std::cout << "WARNING!!! SERVER::DeleteSuture is a depricated function. Use UpdateSutures instead." << std::endl;
}

void SERVER::UpdateSutures( const std::vector< std::pair< std::vector<int>, std::vector<double> > >& new_sutures) {
    typedef std::vector< std::pair<std::vector<int>, std::vector<double> > > suture_desc;
    data_lock.lock();
    data->sutures.clear();
    for( suture_desc::const_iterator it=new_sutures.begin(); it != new_sutures.end(); it++){
        SUTURE newsuture;
        newsuture.triangleA_num = it->first.at(0);
        newsuture.triangleB_num = it->first.at(1);
        newsuture.uvA.assign( it->second.begin(), it->second.begin()+2 );
        newsuture.uvB.assign( it->second.begin()+2, it->second.begin()+4 );
        data->sutures.push_back( newsuture );
    }
    data->timestamp_suture = data->Timestamp();
    data_lock.unlock();
}


void SERVER::SetSimulatorState( const int state ){
    if( state != CLIENT_DATA::NOT_RUNNING && 
        state != CLIENT_DATA::PROCESSING &&
        state != CLIENT_DATA::RUNNING )
        return;
    data_lock.lock();
    data->simulator_state = state;
    data->timestamp_simulatorstate = data->Timestamp();
    data_lock.unlock();
}

int SERVER::getCurrentRevision(){
    return data->global_timestamp;
}

struct DataBuffer {
    char *readptr;
    char *dataptr;
    long sizeleft;
    long size;
};

static size_t read_callback(void *ptr, size_t size, size_t nmemb, void *userp)
{
  struct DataBuffer *buffer = (struct DataBuffer *)userp;

  if(size*nmemb < 1)
    return 0;

  if(buffer->sizeleft) {
    *(char *)ptr = buffer->readptr[0]; /* copy one single byte */
    buffer->readptr++;                 /* advance pointer */
    buffer->sizeleft--;                /* less data left */
    return 1;                        /* we return 1 byte at a time! */
  }

  return 0;                          /* no more data left to deliver */
}


static size_t write_callback(void *contents, size_t size, size_t nmemb, void *userp)
{
  size_t realsize = size * nmemb;
  struct DataBuffer *mem = (struct DataBuffer *)userp;

  mem->dataptr = (char*)(realloc(mem->dataptr, mem->size + realsize + 1));
  if(mem->dataptr == NULL) {
    /* out of memory! */
    printf("not enough memory (realloc returned NULL)\n");
    return 0;
  }

  memcpy(&(mem->dataptr[mem->size]), contents, realsize);
  mem->size += realsize;
  mem->dataptr[mem->size] = 0;

  return realsize;
}

std::string SERVER::GetData( std::string resource, std::string resource_type ){

    std::stringstream ss_url;
    if( resource_type == "SCENE" )
        ss_url << "/data/scene/" << resource;

    if( resource_type == "HISTORY" )
        ss_url << "/data/history/" << resource;

    if( resource_type == "STATIC" )
        ss_url << "/static/data/" << resource;

    return MessageSend( ss_url.str(), "" );
}



void SERVER::PostStatus( std::string message ){
    std::stringstream ss;
    ss << "{ \"token\":\"" << token << "\", \"data\": "  << message << " }";

    MessageSend( "/command/UpdateServerStats", ss.str() );
}


std::string SERVER::MessageSend( std::string url_additional, std::string message ){
    if( token == "" || url == "" )
        return "";

    CURL *curl;
    CURLcode res;
    struct curl_slist *headers = NULL;
    curl = curl_easy_init();

    std::stringstream ss_url;
    ss_url << url;
    std::size_t next_pos = url_additional.find( '/' );
    while( next_pos != std::string::npos ){
        std::size_t pos_after = url_additional.find( '/', next_pos+1 );        
        char* escaped_url = curl_easy_escape( curl, url_additional.c_str()+next_pos+1,
                                              pos_after == std::string::npos ? 0 : pos_after-(next_pos+1));
        ss_url << "/" << escaped_url;
        curl_free( escaped_url );
        next_pos = pos_after;
    }        
    std::string real_url( ss_url.str() );

    DataBuffer chunk;
    chunk.dataptr = (char*)( malloc(1) );  /* will be grown as needed by the realloc above */
    chunk.size = 0;    /* no data at this point */
    std::string data;


    if(curl) {      
        headers = curl_slist_append(headers, "Accept: application/json");
        headers = curl_slist_append(headers, "Content-Type: application/json");
        headers = curl_slist_append(headers, "charsets: utf-8");
        
        /* Set the URL */
                
        std::cout << real_url << std::endl;

        // This option is really important!!!
        // If signals are allowed, DNS timeouts could cause other threads to recieve
        // signals which they can't handle and cause the program to crash randomly.
        // See: http://curl.haxx.se/libcurl/c/libcurl-tutorial.html#Multi-threading
        curl_easy_setopt(curl, CURLOPT_NOSIGNAL, 1L);

        curl_easy_setopt(curl, CURLOPT_URL, real_url.c_str());
        curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
        if( message.size() > 0 )
            curl_easy_setopt(curl, CURLOPT_POST, 1L);
        
        DataBuffer buffer;
        if( message.size() > 0 ){
            buffer.sizeleft = buffer.size = strlen( message.c_str() );
            
            buffer.dataptr = new char[buffer.sizeleft+1];
            strcpy( buffer.dataptr, message.c_str() );
            buffer.dataptr[buffer.sizeleft] = '\0';
            buffer.readptr = buffer.dataptr;

            curl_easy_setopt(curl, CURLOPT_READFUNCTION, read_callback);
            curl_easy_setopt(curl, CURLOPT_READDATA, &buffer);
            curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, buffer.size);
        }

        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void *)&chunk);

        curl_easy_setopt(curl, CURLOPT_USERAGENT, "libcurl-agent/1.0");

        //curl_easy_setopt(curl, CURLOPT_VERBOSE, 1L);
        res = curl_easy_perform(curl);

        if(res != CURLE_OK) {
            std::cout << "curl_easy_perform() failed: " << curl_easy_strerror(res) << std::endl;
            data = "";
        }
        else{
            data = std::string( (char*)(chunk.dataptr) );
        }

        
        curl_easy_cleanup(curl);
        curl_slist_free_all(headers);
        if(message.size() > 0 )
            delete [] buffer.dataptr;
    }
    free(chunk.dataptr);

    return data;
}
