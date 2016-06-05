#ifndef __HTTPVIEWER_DATA_H__
#define __HTTPVIEWER_DATA_H__

#include <vector>
#include <map>
#include <string>
#include <deque>
// #include <functional>  // for callback

#ifndef JSON_SPIRIT_MVALUE_ENABLED
#define JSON_SPIRIT_MVALUE_ENABLED
#endif

//#include <json_spirit/json_spirit_reader_template.h>
//#include <json_spirit/json_spirit_writer_template.h>
#include <json_spirit_reader_template.h>
#include <json_spirit_writer_template.h>



namespace HTTPVIEWER {


struct OBJECT {
    OBJECT(): timestamp_topology(0), timestamp_uvs(0), timestamp_vertices(0), 
              timestamp_normals(0), timestamp_texturename(0), timestamp_normalname(0) {};

    OBJECT(int t): timestamp_topology(t), timestamp_uvs(t), timestamp_vertices(t), 
                   timestamp_normals(t), timestamp_texturename(t), timestamp_normalname(t) {};
    std::vector< int > topology;
    std::vector< double > uvs;
    std::vector< double > vertices;
    std::vector< double > stress;
    std::vector< double > strain;
    std::vector< double > normals;
    std::vector< float > data_texture;
    std::string texturename;
    std::string normalname;
    
    json_spirit::mObject getJson( int since = -1 );
    int max_timestamp() { return std::max(timestamp_topology, std::max(timestamp_uvs, std::max(timestamp_vertices, std::max(timestamp_normals, std::max(timestamp_texturename, timestamp_normalname))))); };
    int timestamp_topology;
    int timestamp_uvs;
    int timestamp_vertices;
    int timestamp_normals;  
    int timestamp_texturename;
    int timestamp_normalname;
    
};

struct TEXTURE {
    std::string texture_data;
    json_spirit::mObject getJson( int since = -1 );
};

struct HOOK {
    int triangle_num;
    std::vector<double> pos;
    json_spirit::mObject getJson( int since = -1 );
};

struct SUTURE {
    int triangleA_num;
    int triangleB_num;
    std::vector<double> uvA, uvB;
    json_spirit::mObject getJson( int since = -1 );
};

struct CLIENT_DATA{

    enum SIMULATOR_STATES { NOT_RUNNING, PROCESSING, RUNNING };

CLIENT_DATA():     dynamic_object(-1), global_timestamp(0), hook_counter(0), simulator_state(NOT_RUNNING), suture_counter(0), 
                   timestamp_command(0), timestamp_texturedata(0), timestamp_hook(0),
                   timestamp_suture(0) {};


    OBJECT dynamic_object;
    std::vector<OBJECT> static_objects;
    std::map<std::string, TEXTURE> textures;
    std::vector<HOOK> hooks;
    std::vector<SUTURE> sutures;
    int simulator_state;
    std::deque<std::string> command_queue;
    
    int timestamp_command;
    int timestamp_texturedata;
    int timestamp_hook;
    int timestamp_suture;
    int timestamp_simulatorstate;

    json_spirit::mObject getJson( int since = -1 );
    int Timestamp();

    int hook_counter;
    int suture_counter;
    int global_timestamp;
};  

}


#endif
