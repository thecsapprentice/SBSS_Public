#include "config.h"
#include <fstream>
#include <iostream>

#ifndef JSON_SPIRIT_MVALUE_ENABLED
#define JSON_SPIRIT_MVALUE_ENABLED
#endif

Configuration::Configuration() {
}

Configuration::Configuration( std::string config_filename ) {
    parseFile( config_filename );
}

Configuration::Configuration( std::vector<std::string> config_filenames ) {
    parseFiles( config_filenames );
}

Configuration::~Configuration() {

}

void Configuration::Parse( std::string config ){
    json_spirit::mValue value;
    try{
        json_spirit::read_string_or_throw(config,value);
    }
    catch( json_spirit::Error_position err ){
        std::cout << "Parsing Failed: (Line " << err.line_ << ") (Col " << err.column_ << "): " << err.reason_ << std::endl;
        throw ConfigurationError( "Error parsing configuration.", "" );
    }
    catch (...){
        throw ConfigurationError( "Error reading configuration.","");
    }
    json_spirit::mObject temp_config;
    try{
        temp_config = value.get_obj();
    }
    catch (...){
        throw ConfigurationError( "Invalid configuration format." );
    }
    
    digestGroup( config_store, temp_config, 0 ); 
}

void Configuration::parseFile( std::string filename ){

    std::ifstream infile(filename);
    if( !infile.good() ){
        throw ConfigurationError( "Configuration file not available: ", std::string( filename ) );
    }
    json_spirit::mValue value;
    try{
        json_spirit::read_stream_or_throw(infile,value);
    }
    catch( json_spirit::Error_position err ){
        std::cout << "Parsing Failed: (Line " << err.line_ << ") (Col " << err.column_ << "): " << err.reason_ << std::endl;
        throw ConfigurationError( "Error parsing configuration file: ", std::string( filename ) );
    }
    catch (...){
        throw ConfigurationError( "Error reading configuration file: ", std::string( filename ) );
    }
    infile.close();

    json_spirit::mObject temp_config;
    try{
        temp_config = value.get_obj();
    }
    catch (...){
        throw ConfigurationError( "Invalid configuration format." );
    }
    
    digestGroup( config_store, temp_config, 0 ); 
}

void Configuration::parseFiles( std::vector<std::string> filenames ){
    bool good=false;
    for( std::vector<std::string>::iterator it = filenames.begin(); it != filenames.end(); it++){
        try{
            std::cout << "Parsing " << *it << std::endl;
            parseFile( *it );
            good = true;
        }
        catch( ConfigurationError e ){
            std::cout << e.what() << std::endl;
            std::cout << "Failed to parse configuration file " << *it << ", continuing..." << std::endl;
        }
    }
    if( !good )
        throw ConfigurationError( "Unable to parse any configuration files. Remaining unconfigured." );    
}

void Configuration::digestGroup( json_spirit::mObject& store, json_spirit::mObject& group, int level ){

    if( level > 1 )
        throw ConfigurationError( "Only single level groups are supported." );


    typedef json_spirit::mObject::iterator iter;
    for( iter it = group.begin(); it != group.end(); it++){
        if( it->second.type() == json_spirit::obj_type ){           
            iter key_find = store.find( it->first );
            if( key_find == store.end() )
                store.insert( std::pair<std::string, json_spirit::mValue>( it->first, json_spirit::mValue( json_spirit::mObject())));
            try{
                json_spirit::mObject& new_store = store[it->first].get_obj();
                digestGroup( new_store, it->second.get_obj(), level+1 );            
            }
            catch( ... ){
                throw ConfigurationError( "Failure: Type mismatch on configuration key overwrite." );
            }
        }
        
        if( it->second.type() == json_spirit::str_type ){
            iter key_find = store.find( it->first );
            if( key_find == store.end() )
                store.insert( std::pair<std::string, json_spirit::mValue>( it->first, json_spirit::mValue( it->second.get_str())));
            else{
                if( key_find->second.type() != json_spirit::str_type )
                    throw ConfigurationError( "Failure: Type mismatch on configuration key overwrite." );                
                key_find->second = json_spirit::mValue(it->second.get_str());
            }            
        }

        if( it->second.type() == json_spirit::bool_type ){
            iter key_find = store.find( it->first );
            if( key_find == store.end() )
                store.insert( std::pair<std::string, json_spirit::mValue>( it->first, json_spirit::mValue( it->second.get_bool())));
            else{
                if( key_find->second.type() != json_spirit::bool_type )
                    throw ConfigurationError( "Failure: Type mismatch on configuration key overwrite." );                
                key_find->second = json_spirit::mValue(it->second.get_bool());
            }    
        }

        if( it->second.type() == json_spirit::int_type ){
            iter key_find = store.find( it->first );
            if( key_find == store.end() )
                store.insert( std::pair<std::string, json_spirit::mValue>( it->first, json_spirit::mValue( it->second.get_int())));
            else{
                if( key_find->second.type() != json_spirit::int_type )
                    throw ConfigurationError( "Failure: Type mismatch on configuration key overwrite." );                
                key_find->second = json_spirit::mValue(it->second.get_int());
            }    
        }

        if( it->second.type() == json_spirit::real_type ){
            iter key_find = store.find( it->first );
            if( key_find == store.end() )
                store.insert( std::pair<std::string, json_spirit::mValue>( it->first, json_spirit::mValue( it->second.get_real())));
            else{
                if( key_find->second.type() != json_spirit::real_type )
                    throw ConfigurationError( "Failure: Type mismatch on configuration key overwrite." );                
                key_find->second = json_spirit::mValue(it->second.get_real());
            }    
        }

        if( it->second.type() == json_spirit::array_type || 
            it->second.type() == json_spirit::null_type ){
            throw ConfigurationError( "Cannot support configuration values of this type." );
        }
    }
}

const json_spirit::mObject& Configuration::GetGroup( std::string group ) const {
    typedef json_spirit::mObject::const_iterator iter;

    if( group == "" )
        return config_store;
    else{
        iter group_find = config_store.find( group );
        if( group_find != config_store.end() ){
            if( group_find->second.type() != json_spirit::obj_type )
                throw ConfigurationError( std::string("Group '") + group + std::string("' not found." ));
            
            return group_find->second.get_obj();
        }
        else
            throw ConfigurationError( std::string("Group '") + group + std::string("' not found."));
    }
}


json_spirit::mValue Configuration::Get( std::string group, std::string key ) const {
    typedef json_spirit::mObject::const_iterator iter;
        
    const json_spirit::mObject& g = GetGroup( group );
    iter key_find = g.find( key );
    if( key_find != g.end() )
        return key_find->second;
    else
        throw ConfigurationError( std::string("Key '") + key + std::string("' not found."));
}


void Configuration::Get( int& value, const std::string group, const std::string key ) const {
    int _value = 0;
    try{
        _value = Get(group, key).get_int();
    }
    catch( ... ){
        throw ConfigurationError( std::string("Integer key '") + group +"::"+ key + std::string("' not found.") );
    }
    value = _value;
}
void Configuration::GetDefault( int& value, const int default_value, const std::string group, const std::string key ) const {
    try{
        Get( value, group, key );
    }
    catch( ConfigurationError e ){
        value = default_value;
    }
}

void Configuration::Get( bool& value, const std::string group, const std::string key ) const {
    bool _value = false;
    try{
        _value = Get(group, key).get_bool();
    }
    catch( ... ){
        throw ConfigurationError( std::string("Boolean key '") + group +"::"+ key + std::string("' not found.") );
    }
    value = _value;
}
void Configuration::GetDefault( bool& value, const bool default_value, const std::string group, const std::string key ) const {
    try{
        Get( value, group, key );
    }
    catch( ConfigurationError e ){
        value = default_value;
    }
}


void Configuration::Get( std::string& value, const std::string group, const std::string key ) const {
    std::string _value = "";
    try{
        _value = Get(group, key).get_str();
    }
    catch( ... ){
        throw ConfigurationError( std::string("String key '") + group +"::"+ key + std::string("' not found.") );
    }
    value = _value;
}
void Configuration::GetDefault( std::string& value, const std::string default_value, const std::string group, const std::string key ) const {
    try{
        Get( value, group, key );
    }
    catch( ConfigurationError e ){
        value = default_value;
    }
}

void Configuration::Get( float& value, const std::string group, const std::string key ) const {
    float _value = 0.0f;
    try{
        _value = Get(group, key).get_real();
    }
    catch( ... ){
        throw ConfigurationError( std::string("Float key '") + group +"::"+ key + std::string("' not found.") );
    }
    value = _value;
}
void Configuration::GetDefault( float& value, const float default_value, const std::string group, const std::string key ) const {
    try{
        Get( value, group, key );
    }
    catch( ConfigurationError e ){
        value = default_value;
    }
}

std::vector<std::string> Configuration::GetKeys( const std::string group ) const {
    typedef json_spirit::mObject::const_iterator iter;
    std::vector<std::string> keys;

    const json_spirit::mObject& g = GetGroup( group );
    iter key_find = g.begin();
    while( key_find != g.end() ){
        keys.push_back( key_find->first );
        key_find++;
    }
    
    return keys;
}

