#include "DATA.h"

using namespace HTTPVIEWER;

json_spirit::mObject CLIENT_DATA::getJson( int since ){

    json_spirit::mObject object;

    if( dynamic_object.max_timestamp() >= since )
        object["dynamic"] = dynamic_object.getJson(since);

    object["static"] = json_spirit::mArray();
    for( std::vector<OBJECT>::iterator it = static_objects.begin(); it != static_objects.end(); it++ )
        object["static"].get_array().push_back( it->getJson( since ) );

    if( timestamp_texturedata >= since ){
        object["textures"] = json_spirit::mObject();
        for( std::map<std::string, TEXTURE>::iterator it = textures.begin(); it != textures.end(); it++ )
            object["textures"].get_obj()[it->first] = it->second.getJson(since);
    }
    
    if( timestamp_hook >= since ){
        object["hooks"] = json_spirit::mArray();
        for( std::vector<HOOK>::iterator it = hooks.begin(); it != hooks.end(); it++ )
            object["hooks"].get_array().push_back( it->getJson() );
    }

    if( timestamp_suture >= since ){
        object["sutures"] = json_spirit::mArray();
        for( std::vector<SUTURE>::iterator it = sutures.begin(); it != sutures.end(); it++ )
            object["sutures"].get_array().push_back( it->getJson() );
    }

    // Always send the simulator state, even if its timestamp is old.
    object["simulator_state"] = json_spirit::mValue( simulator_state );

    //if( timestamp_command >= since ){
    // Nathan - clearing command will end client wait state. See SERVER::setClientWait()
    object["command"] = json_spirit::mArray();
    for( std::deque<std::string>::iterator it = command_queue.begin(); it != command_queue.end(); it++ )
        object["command"].get_array().push_back( json_spirit::mValue(*it) );
    //}

    object["timestamp"] = json_spirit::mValue( global_timestamp );

    return object;
}

int CLIENT_DATA::Timestamp(){
    int time = global_timestamp;
    global_timestamp++;
    return time;
}


json_spirit::mObject OBJECT::getJson( int since ){
    json_spirit::mObject object;
    
    if( timestamp_topology >= since ){
        object["topology"] = json_spirit::mArray();
        for( std::vector<int>::iterator it = topology.begin(); it != topology.end(); it++ )
            object["topology"].get_array().push_back( *it );
    }

    if( timestamp_uvs >= since ){
        object["uvs"] = json_spirit::mArray();
        for( std::vector<double>::iterator it = uvs.begin(); it != uvs.end(); it++ )
            object["uvs"].get_array().push_back( *it );
    }

    if( timestamp_vertices >= since ){
        object["vertices"] = json_spirit::mArray();
        for( std::vector<double>::iterator it = vertices.begin(); it != vertices.end(); it++ )
            object["vertices"].get_array().push_back( *it );
        if( stress.size() > 0 ){
            object["stress"] = json_spirit::mArray();
            for( std::vector<double>::iterator it = stress.begin(); it != stress.end(); it++ )
                object["stress"].get_array().push_back( *it );
        }
        if( strain.size() > 0 ){
            object["strain"] = json_spirit::mArray();
            for( std::vector<double>::iterator it = strain.begin(); it != strain.end(); it++ )
                object["strain"].get_array().push_back( *it );
        }
        if( data_texture.size() > 0 ){
            object["data_texture"] = json_spirit::mArray();
            for( std::vector<float>::iterator it = data_texture.begin(); it != data_texture.end(); it++ )
                object["data_texture"].get_array().push_back( (double)(*it) );
        }
    }

    if( timestamp_normals >= since ){
        object["normals"] = json_spirit::mArray();
        for( std::vector<double>::iterator it = normals.begin(); it != normals.end(); it++ )
            object["normals"].get_array().push_back( *it );
    }

    if( timestamp_texturename >= since ){
        object["texturename"] = json_spirit::mValue(texturename);
    }

    if( timestamp_normalname >= since ){
        object["normalname"] = json_spirit::mValue(normalname);
    }

    return object;
}

json_spirit::mObject TEXTURE::getJson( int since ){
    json_spirit::mObject object;
    object["data"] = json_spirit::mValue( texture_data );
    return object;
}

json_spirit::mObject HOOK::getJson( int since ){
    json_spirit::mObject object;
    object["triangle"] = json_spirit::mValue( triangle_num );
    object["pos"] = json_spirit::mArray();
    for( std::vector<double>::iterator it = pos.begin(); it != pos.end(); it++ )
        object["pos"].get_array().push_back( *it );
    return object;
}

json_spirit::mObject SUTURE::getJson( int since ){
    json_spirit::mObject object;
    object["A"] = json_spirit::mObject();
    object["A"].get_obj()["triangle"] = json_spirit::mValue( triangleA_num );
    object["A"].get_obj()["uv"] = json_spirit::mArray();
    for( std::vector<double>::iterator it = uvA.begin(); it != uvA.end(); it++ )
        object["A"].get_obj()["uv"].get_array().push_back( *it );

    object["B"] = json_spirit::mObject();
    object["B"].get_obj()["triangle"] = json_spirit::mValue( triangleB_num );
    object["B"].get_obj()["uv"] = json_spirit::mArray();
    for( std::vector<double>::iterator it = uvB.begin(); it != uvB.end(); it++ )
        object["B"].get_obj()["uv"].get_array().push_back( *it );

    return object;
}

