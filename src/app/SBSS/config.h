#ifndef __SBSS_CONFIG_H__
#define __SBSS_CONFIG_H__

#include <string>
#include <vector>
#include <exception>
#include "json_spirit_reader_template.h"
#include "json_spirit_writer_template.h"

#ifndef JSON_SPIRIT_MVALUE_ENABLED
#define JSON_SPIRIT_MVALUE_ENABLED
#endif

class ConfigurationError : public std::exception {
public:
    ConfigurationError( std::string msg ) throw() : _msg(std::string("[CONFIGURATION ERROR] ")+msg )  {}
    ConfigurationError( std::string msg, std::string file) throw() : _msg(std::string("[CONFIGURATION ERROR] ") + msg + file ) {}
    virtual ~ConfigurationError() throw() {};
    virtual const char* what() const throw() {
        return _msg.c_str();
    }
private:
    std::string _msg;
};


class Configuration {
    
 public:
    Configuration();
    Configuration( std::string config_filename );
    Configuration( std::vector<std::string> config_filenames );
    ~Configuration();  
   
    void Parse( std::string config );
    void parseFile( std::string filename );
    void parseFiles( std::vector<std::string> filenames );

    void Get( int& value, const std::string group, const std::string key ) const;
    void Get( bool& value, const std::string group, const std::string key ) const;
    void Get( std::string& value, const std::string group, const std::string key ) const;
    void Get( float& value, const std::string group, const std::string key ) const;

    void GetDefault( int& value, const int default_value,const std::string group, const std::string key ) const;
    void GetDefault( bool& value, const bool default_value,const std::string group, const std::string key ) const;
    void GetDefault( std::string& value, const std::string default_value, const std::string group, const std::string key ) const;
    void GetDefault( float& value, const float default_value, const std::string group, const std::string key ) const;

    std::vector<std::string> GetKeys( const std::string group ) const ;

    template< typename T >
    inline void Get( T& value, const std::string key ) const {
        Get( value, "", key );
    }

    template< typename T >
    inline void GetDefault( T& value, const T default_value, const std::string key ) const {
        GetDefault( value, default_value, "",  key );
    }

 private:
    bool fileExists( std::string filename );
    void digestGroup( json_spirit::mObject&, json_spirit::mObject&, int level );  
    json_spirit::mValue Get( std::string group, std::string key ) const;
    const json_spirit::mObject& GetGroup( std::string group ) const;

    json_spirit::mObject config_store;
};

#endif
