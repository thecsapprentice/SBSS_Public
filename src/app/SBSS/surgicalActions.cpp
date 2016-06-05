//#include "wx/wx.h"
#include "GraphicsUtils/wxGraphics.h"
#include "GraphicsUtils/glslTriangle.h"
#include "GraphicsUtils/Vec3f.h"
#include <sstream>
#include <fstream>
#include "surgicalActions.h"
#include <exception>
#include <stdexcept>
#include <vector>
#include "GraphicsUtils/GLmatrices.h"

#include <HTTPViewer/SERVER.h>
#include <HTTPViewer/DATA.h>

surgicalActions::surgicalActions(const Configuration* config) : _config(config), _cle(config)
{
    int port;
    std::string url;
    std::string token;
    config->Get( port, "port" );
    config->GetDefault( token, std::string(""), "idtoken");
    config->GetDefault( url, std::string(""), "master");
    config->GetDefault( use_network_resources, false, "use_net_resources" );
    http_server = new HTTPVIEWER::SERVER( port, token, url );
    getServerHandle().CreateServer();
    http_server->SetSimulatorState(HTTPVIEWER::CLIENT_DATA::NOT_RUNNING);

    int mode;
    config->Get( mode, "", "defaultMode" );
    std::cout << "Starting with default server mode " << mode << std::endl;
    setServerMode( mode );

	_cle.setSurgicalActions(this);
    current_history_position = command_history.begin();
    next_history_position = command_history.end();

    std::vector<std::string> scene_names = config->GetKeys( "scenes" );
    for( std::vector<std::string>::iterator iter = scene_names.begin(); iter != scene_names.end(); iter++){
        std::string scene_file;
        config->Get( scene_file, "scenes", *iter );
        availableScenes[*iter] = scene_file;
    }
}

surgicalActions::~surgicalActions()
{
    if(_serverMode>0)
        getServerHandle().StopServer();

    getServerHandle().DestroyServer();
    delete http_server;
}

void surgicalActions::setServerMode(int mode)
{
// 0=Off, 1=Passive client, 2=Active client
	if(mode>0) {
        getServerHandle().StartServer();
	}
    else{
        getServerHandle().StopServer();
    }
	_serverMode = mode;    
}

void surgicalActions::update() {
    // Push new commands from the network
    checkClientCommands();

    // Process the next command yet to be processed
    if( !raw_command_queue.empty() ){
        std::string commandStr = raw_command_queue.front();
        raw_command_queue.pop_front();
        try{
            processClientCommand( UnpackCommand(commandStr) );
        }
        catch( std::runtime_error err ){
            std::cout << "Command failed: " + std::string(err.what()) << std::endl;
        }
    }

    // Update the Physics
    updatePhysics();
}

void surgicalActions::pushCommand( std::string commandStr ){
    raw_command_queue.push_back( commandStr );
}    

void surgicalActions::updatePhysics() {
    _cle.cleUpdate();
}


//=========================================================================================
//=========================================================================================
//=========================================================================================
//
//
//
//                                   UTILITY METHODS
//
//
//
//=========================================================================================
//=========================================================================================
//=========================================================================================


const json_spirit::mValue surgicalActions::getOrRaise( const json_spirit::mObject& obj, const std::string& key){
    json_spirit::mObject::const_iterator iter = obj.find( key );
    if( iter != obj.end() )
        return iter->second;
    throw std::out_of_range("Key '"+ key +"' not found.");    
}

//=========================================================================================
//=========================================================================================
//=========================================================================================
//
//
//
//                             COMMAND PROCESSING METHODS
//
//
//
//=========================================================================================
//=========================================================================================
//=========================================================================================



void surgicalActions::checkClientCommands() {
    if( getServerMode() > 1 ){
        // wait state set on client and server by client issuing command       
        std::string clientCommand = getServerHandle().getClientCommand();        
        if(clientCommand != ""){
            std::cout << "Incoming client command: " << clientCommand << std::endl;
            raw_command_queue.push_back( clientCommand );
            //getServerHandle().clearClientCommand();
        }
    }
}

json_spirit::mObject surgicalActions::UnpackCommand( std::string command_str ){
    std::cout << "Parsing client command:" << std::endl << command_str << std::endl;

    json_spirit::mValue value;
    try{
        json_spirit::read_string_or_throw(command_str,value);
    }
    catch (json_spirit::Error_position err){
        std::cout << "Parsing Failed: (Line " << err.line_ << ") (Col " << err.column_ << "): " << err.reason_ << std::endl;
        throw std::runtime_error("Parsing error. Ignoring Command");
    }
    catch (...){
        throw std::runtime_error("Unknown exception parsing command string. Ignoring command.");
    }
    try{
        json_spirit::mObject obj = value.get_obj();
        return obj;
    }
    catch(...){        
        throw std::runtime_error("Command is not structured as an object. Ignoring command.");
    }    
}

void surgicalActions::processClientCommand( const json_spirit::mObject commandObj){
    processCommand( commandObj, true );
}

void surgicalActions::processHistoryCommand( const json_spirit::mObject commandObj){
    processCommand( commandObj, false );
}

void surgicalActions::processCommand( const json_spirit::mObject& commandObj, bool can_invalidate_future_history ){
    std::string command;
    json_spirit::mObject data;
    try{
        command = getOrRaise(commandObj,"command").get_str();
    }
    catch (...){
        throw std::runtime_error("No command identifier found. Ignoring command.");
    }
    try{
        data = getOrRaise(commandObj,"data").get_obj();
    }
    catch (...){
        throw std::runtime_error("No command data found. Ignoring command.");
    }    

    //
    // COMMAND ---- LOAD HISTORY
    //

    if( command == "loadHistory" ){
        std::string history_file;
        try{
            history_file = getOrRaise(data,"name").get_str();
        }
        catch (...){
            throw std::runtime_error("History filename not specified. Ignoring command.");
        }

        if( can_invalidate_future_history == false )
            throw std::runtime_error("Cannot load new history without invalidating future history. Ignoring command" );

        std::cout << "Loading history file " << history_file << std::endl;
        loadHistory(history_file.c_str());
        return;
    }

    //
    // COMMAND ---- SAVE HISTORY
    //

    if( command == "saveHistory" ){
        std::string history_file;
        try{
            history_file = getOrRaise(data,"name").get_str();
        }
        catch (...){
            throw std::runtime_error("History filename not specified. Ignoring command.");
        }

        saveHistory(history_file.c_str());
        return;
    }

    //
    // COMMAND ---- LOAD SCENE
    //

    if( command == "loadScene" ){
        std::string scene;
        std::string scene_file;
        try{
            scene = getOrRaise(data,"name").get_str();
        }
        catch (...){
            throw std::runtime_error("Scene not specified. Ignoring command.");
        }
        std::cout << "Loading scene " << scene << std::endl;
        loadScene(scene.c_str());
        if(can_invalidate_future_history) truncateAndAppendHistory( commandObj );
        return;
    }

    //
    // COMMAND ---- ADVANCE HISTORY
    //

    if(command == "historyNext"){
        nextHistoryAction();
        return;
    }

    //
    // COMMAND ---- ADD HOOK
    //

    if(command == "addHook"){
        int triangle;
        std::vector<float> barycentric_coords;

        try{
            triangle = getOrRaise(data,"triangle").get_int();
        }
        catch (...){
            throw std::runtime_error("Triangle id not specified. Ignoring command.");
        }
        try{
            json_spirit::mArray arr = getOrRaise(data,"coords").get_array();
            barycentric_coords.push_back( (float)(arr.at(0).get_real()) );
            barycentric_coords.push_back( (float)(arr.at(1).get_real()) );
            barycentric_coords.push_back( (float)(arr.at(2).get_real()) );
        }
        catch (...){
            throw std::runtime_error("Barycentric coordinates not specified. Ignoring command.");
        }

        Add_Hook( triangle, barycentric_coords );
        if(can_invalidate_future_history) truncateAndAppendHistory( commandObj );
        return;
    }

    //
    // COMMAND ---- MOVE HOOK
    //

    if(command == "moveHook"){
        int hook_id;
        std::vector<float> coords;

        try{
            hook_id = getOrRaise(data,"hook_id").get_int();
        }
        catch (...){
            throw std::runtime_error("Hook id not specified. Ignoring command.");
        }
        try{
            json_spirit::mArray arr = getOrRaise(data,"coords").get_array();
            coords.push_back( (float)(arr.at(0).get_real()) );
            coords.push_back( (float)(arr.at(1).get_real()) );
            coords.push_back( (float)(arr.at(2).get_real()) );
        }
        catch (...){
            throw std::runtime_error("Coordinates not specified. Ignoring command.");
        }

        Move_Hook( hook_id, coords );
        if(can_invalidate_future_history) truncateAndAppendHistory( commandObj );
        return;
    }

    //
    // COMMAND ---- DELETE HOOK
    //

    if(command == "deleteHook"){
        int hook_id;

        try{
            hook_id = getOrRaise(data,"hook_id").get_int();
        }
        catch (...){
            throw std::runtime_error("Hook id not specified. Ignoring command.");
        }

        Remove_Hook( hook_id );
        if(can_invalidate_future_history) truncateAndAppendHistory( commandObj );
        return;
    }

    //
    // COMMAND ---- ADD SUTURE
    //

    if(command == "addSuture"){
        int triangle1;
        std::vector<float> barycentric_coords1;
        int triangle2;
        std::vector<float> barycentric_coords2;

        try{
            triangle1 = getOrRaise(data,"triangleA").get_int();
        }
        catch (...){
            throw std::runtime_error("Triangle 1 id not specified. Ignoring command.");
        }
        try{
            triangle2 = getOrRaise(data,"triangleB").get_int();
        }
        catch (...){
            throw std::runtime_error("Triangle 2 id not specified. Ignoring command.");
        }
        try{
            json_spirit::mArray arr = getOrRaise(data,"uvA").get_array();
            barycentric_coords1.push_back( (float)(arr.at(0).get_real()) );
            barycentric_coords1.push_back( (float)(arr.at(1).get_real()) );
            barycentric_coords1.push_back( (float)(arr.at(2).get_real()) );
        }
        catch (...){
            throw std::runtime_error("Barycentric coordinates 1 not specified. Ignoring command.");
        }
        try{
            json_spirit::mArray arr = getOrRaise(data,"uvB").get_array();
            barycentric_coords2.push_back( (float)(arr.at(0).get_real()) );
            barycentric_coords2.push_back( (float)(arr.at(1).get_real()) );
            barycentric_coords2.push_back( (float)(arr.at(2).get_real()) );
        }
        catch (...){
            throw std::runtime_error("Barycentric coordinates 2 not specified. Ignoring command.");
        }


        Add_Suture( triangle1, barycentric_coords1, triangle2, barycentric_coords2 );
        if(can_invalidate_future_history) truncateAndAppendHistory( commandObj );
        return;
    }

    //
    // COMMAND ---- DELETE SUTURE
    //

    if(command == "deleteSuture"){
        int suture_id;

        try{
            suture_id = getOrRaise(data,"suture_id").get_int();
        }
        catch (...){
            throw std::runtime_error("Suture id not specified. Ignoring command.");
        }

        Remove_Suture( suture_id );
        if(can_invalidate_future_history) truncateAndAppendHistory( commandObj );
        return;
    }

    //
    // COMMAND ---- MAKE INCISION
    //
    if( command == "makeIncision" ){
        std::vector<int> path_triangles;
        std::vector<float> path_uvs;
        std::vector<float> path_positions;
        std::vector<float> path_normals;
        bool edge_start = false;
        bool edge_end = false;

        json_spirit::mArray path_array;
        json_spirit::mArray::iterator it;
        try{
            path_array = getOrRaise( data, "path").get_array();
        }
        catch(...){
            throw std::runtime_error( "Path data not found. Ignoring command.");
        }

        for( it=path_array.begin(); it != path_array.end(); it++ ){
            try {
                int tri = getOrRaise( it->get_obj(), "triangle" ).get_int();
                path_triangles.push_back( tri );
            }
            catch(...){
                throw std::runtime_error( "Path data missing triangle. Ignoring command.");
            }

            try {
                json_spirit::mArray arr = getOrRaise( it->get_obj() ,"coords").get_array();
                path_uvs.push_back( (float)(arr.at(0).get_real()) );
                path_uvs.push_back( (float)(arr.at(1).get_real()) );
            }
            catch(...){
                throw std::runtime_error( "Path data missing UVs. Ignoring command.");
            }

            try {
                json_spirit::mArray arr = getOrRaise( it->get_obj() ,"position").get_array();
                path_positions.push_back( (float)(arr.at(0).get_real()) );
                path_positions.push_back( (float)(arr.at(1).get_real()) );
                path_positions.push_back( (float)(arr.at(2).get_real()) );
            }
            catch(...){
                throw std::runtime_error( "Path data missing Positions. Ignoring command.");
            }


            try {
                json_spirit::mArray arr = getOrRaise( it->get_obj() ,"normal").get_array();
                path_normals.push_back( (float)(arr.at(0).get_real()) );
                path_normals.push_back( (float)(arr.at(1).get_real()) );
                path_normals.push_back( (float)(arr.at(2).get_real()) );
            }
            catch(...){
                throw std::runtime_error( "Path data missing Normals. Ignoring command.");
            }
        }

        try{
            edge_start = getOrRaise( data, "edge_start").get_bool();
        }
        catch(...){
            throw std::runtime_error( "Path start specifier not found. Ignoring command.");
        }

        try{
            edge_end = getOrRaise( data, "edge_end").get_bool();
        }
        catch(...){
            throw std::runtime_error( "Path end specifier not found. Ignoring command.");
        }
        
        CreateIncision( path_triangles, path_uvs, path_positions, path_normals, edge_start, edge_end );
        if(can_invalidate_future_history) truncateAndAppendHistory( commandObj );
        return;
    }

    //
    // COMMAND ---- EXCISE REGION
    //
    if( command == "exciseRegion" ){
        int triangle_id;

        try{
            triangle_id = getOrRaise(data,"triangle").get_int();
        }
        catch (...){
            throw std::runtime_error("Triangle id not specified. Ignoring command.");
        }

        ExciseRegion( triangle_id );
        if(can_invalidate_future_history) truncateAndAppendHistory( commandObj );
        return;
    }


    throw std::runtime_error("Invalid javascript client command:  " + command);
}


//=========================================================================================
//=========================================================================================
//=========================================================================================
//
//
//
//                             SCENE AND ENGINE METHODS
//
//
//
//=========================================================================================
//=========================================================================================
//=========================================================================================


void surgicalActions::loadScene(const char *scene)
{
    std::string scene_data = loadData( scene, "SCENE" );
    _sutures.clear();
    _hooks.clear();
    reloadEngine();
	bool ret = _cle.loadScene(scene_data);
	if( !ret ){
        throw std::runtime_error( "Failed to load scene file." );
    }

    std::stringstream ss;
    ss << "{\"scene\":\"" << scene << "\"}";
    getServerHandle().PostStatus( ss.str() );
}

std::string surgicalActions::loadData( std::string resource, std::string resource_type ){
    std::string data;

    if( use_network_resources ){        
        data = getServerHandle().GetData( resource, resource_type );
    }
    else{
        std::string data_path;
        if( resource_type == "SCENE" || resource_type == "STATIC" ){
            _config->Get( data_path, "", "scenePath" );

            if( resource_type == "SCENE" ){
                std::map<std::string, std::string>::iterator iter;
                iter = availableScenes.find( resource );
                if( iter == availableScenes.end() )
                    throw std::runtime_error( std::string("Scene unknown. Failed to load.") );
                data_path.append( iter->second );
            }
            else
                data_path.append( resource );
        }
        if( resource_type == "HISTORY" ){
            _config->Get( data_path, "", "historyPath" );
            data_path.append( resource );
            data_path.append( ".hst" );
        }
        
        std::ifstream is(data_path);
        if( !is.good() )
            throw std::runtime_error( std::string("Failed to open data at path: ").append( data_path ) );

        std::string line;
        do{
            std::getline( is, line );
            line.push_back( '\n' );
            data.append( line );
        }while( is.good() );              
    }

    return data;
}


void surgicalActions::reloadEngine(){
    _cle.reloadEngine();

    if(getServerMode()>0){
        getServerHandle().StopServer();
        getServerHandle().DestroyServer();
        getServerHandle().CreateServer();
        getServerHandle().StartServer();
        http_server->SetSimulatorState(HTTPVIEWER::CLIENT_DATA::NOT_RUNNING);
    }
}

//=========================================================================================
//=========================================================================================
//=========================================================================================
//
//
//
//                                     HISTORY METHODS
//
//
//
//=========================================================================================
//=========================================================================================
//=========================================================================================

void surgicalActions::saveHistory(const char *fileName)
{
    std::cout << "Saving history file...  " << fileName << std::endl;
    if( use_network_resources ){
        std::string history_dump = dumpHistory();
        std::stringstream ss_url;
        ss_url << "/command/saveHistory/" << fileName;
        getServerHandle().MessageSend( ss_url.str(), history_dump );
    }
    else{
        std::ofstream outfile;
        outfile.exceptions ( std::ofstream::failbit | std::ofstream::badbit );
        try{
            outfile.open( fileName );
            outfile << dumpHistory();
            outfile.close();
        }
        catch( std::ios_base::failure e ){
            throw std::runtime_error( std::string("Failed to open history file for writing: ") + std::string( e.what() ) );
        }
    }
}

void surgicalActions::loadHistory(const char *history)
{

    std::string history_data = loadData( history, "HISTORY" );

    std::list< json_spirit::mObject > temp_history;
        
    json_spirit::mValue value;
    try{
        json_spirit::read_string_or_throw(history_data,value);
    }
    catch(json_spirit::Error_position err){
        std::cout << "Parsing Failed: (Line " << err.line_ << ") (Col " << err.column_ << "): " << err.reason_ << std::endl;
        throw std::runtime_error("Parsing error. Ignoring Command");
    }
    catch (...){
        throw std::runtime_error("Unknown exception parsing history file. Ignoring command.");
    }
    
    json_spirit::mArray commands;
    try{
        commands = value.get_array();
        for( json_spirit::mArray::iterator it = commands.begin(); it != commands.end(); it++){
            json_spirit::mObject obj = it->get_obj();
            temp_history.push_back( obj );
        }
    }
    catch(...){       
        throw std::runtime_error("Invalid history file format. Ignoring command.");
    }
  
     std::stringstream ss;
     ss << "{\"history\":\"" << history << "\"}";
     getServerHandle().PostStatus( ss.str() );

    // replace our command_history with the one just loaded.
    // NOTE: if the above load fails, our current history has been preserved.
    command_history.swap( temp_history );
    
    current_history_position = command_history.begin();
    next_history_position = current_history_position;
    
    // Automatically exectute first command in history
    nextHistoryAction();
}

void surgicalActions::nextHistoryAction(){
    try {
        consumeNextHistory();
    }
    catch( std::runtime_error e ){
        // This would be a good place to send a message to the client
        // For now, pass the error up the chain.
        throw e;        
    }
    catch( ... ){
        throw std::runtime_error("Performing next history action failed for unknown reason.");
    }
}

void surgicalActions::consumeNextHistory(){
    if( next_history_position == command_history.end() )
        throw std::runtime_error( "End of history reached." ); 
    processHistoryCommand( *next_history_position );
    current_history_position = next_history_position;
    next_history_position++;
}

void surgicalActions::truncateAndAppendHistory( const json_spirit::mObject& commandObj ){
    if( next_history_position != command_history.end() )
        command_history.erase( next_history_position, command_history.end() );
    appendHistory( commandObj );
    current_history_position = command_history.end()--;
    next_history_position = command_history.end();
    std::stringstream ss;
    ss << "{\"history\":\"Custom\"}";
    getServerHandle().PostStatus( ss.str() );
}

void surgicalActions::appendHistory( const json_spirit::mObject& commandObj ){
    command_history.push_back( commandObj );
    //std::cout << dumpHistory() << std::endl;
}

std::string surgicalActions::dumpHistory(){
    json_spirit::mArray commands;
    std::string dump;
    typedef std::list< json_spirit::mObject >::const_iterator iter;
    for( iter it = command_history.begin(); it != command_history.end(); it++){
        commands.push_back( *it );
    }
    dump = json_spirit::write_string( json_spirit::mValue( commands ), json_spirit::pretty_print);
    return dump;
}
//=========================================================================================
//=========================================================================================
//=========================================================================================
//
//
//
//                                     HOOK METHODS
//
//
//
//=========================================================================================
//=========================================================================================
//=========================================================================================


void surgicalActions::Add_Hook( const int triangle, const std::vector<float>& barycentric_coords ){
    float uv[2];
    float pos[3];
    uv[0] = barycentric_coords[1]; // We compute the first coordinate implicitly
    uv[1] = barycentric_coords[2];
    hook newhook;
    newhook.uv.assign( uv, uv+2 );
    newhook.triangle = triangle;
    newhook._cle_hookid = _cle.addHook( triangle, uv, pos );
    newhook.position.assign( pos, pos+3 );
    _hooks.push_back( newhook );
    PublishUpdatedHooks();
}

void surgicalActions::Move_Hook( const int hook_id, const std::vector<float>& coords ){
    float pos[3];
    pos[0] = coords.at(0);
    pos[1] = coords.at(1);
    pos[2] = coords.at(2);

    std::list<hook>::iterator it = _hooks.begin();
    for( int i = 0; i < hook_id; i++)
        it++;
    it->position.assign( coords.begin(), coords.end() );  
    _cle.setHookPosition(pos, it->_cle_hookid );

    PublishUpdatedHooks();
}

void surgicalActions::Remove_Hook( const int hook_id ){
    std::list<hook>::iterator it = _hooks.begin();
    for( int i = 0; i < hook_id; i++)
        it++;

    _cle.deleteHook(it->_cle_hookid);	// must delete the cleNode first
    _hooks.erase(it);

    PublishUpdatedHooks();
}

void surgicalActions::PublishUpdatedHooks(){
    if( getServerMode() > 0 ){       
        typedef std::vector< std::pair<int, std::vector<double> > > hook_desc;
        hook_desc hook_updates;
        std::list<hook>::iterator it = _hooks.begin();
        while( it != _hooks.end() ){
            int triangle = it->triangle;
            std::vector<double> dpos;
            dpos.push_back( it->position.at(0) );
            dpos.push_back( it->position.at(1) );
            dpos.push_back( it->position.at(2) );
            hook_updates.push_back( std::pair<int, std::vector<double> >(triangle, dpos) );
            it++;
        }
        getServerHandle().UpdateHooks( hook_updates );        
    }
}

//=========================================================================================
//=========================================================================================
//=========================================================================================
//
//
//
//                                  SUTURE METHODS
//
//
//
//=========================================================================================
//=========================================================================================
//=========================================================================================


void surgicalActions::Add_Suture( const int triangle1, const std::vector<float>& barycentric_coords1,
                                  const int triangle2, const std::vector<float>& barycentric_coords2){
    float uv1[2];
    float uv2[2];

    uv1[0] = barycentric_coords1[1]; // We compute the first coordinate implicitly
    uv1[1] = barycentric_coords1[2];
    uv2[0] = barycentric_coords2[1]; // We compute the first coordinate implicitly
    uv2[1] = barycentric_coords2[2];

    suture newsuture;
    newsuture.uv0.assign( uv1, uv1+2 );
    newsuture.uv1.assign( uv2, uv2+2 );
    newsuture.triangle0 = triangle1;
    newsuture.triangle1 = triangle2;
    newsuture._cle_sutureid = _cle.addSuture( triangle1, uv1, triangle2, uv2  );

    _sutures.push_back( newsuture );
    PublishUpdatedSutures();
}

void surgicalActions::Remove_Suture( const int suture_id ){
    std::list<suture>::iterator it = _sutures.begin();
    for( int i = 0; i < suture_id; i++)
        it++;

    _cle.deleteSuture(it->_cle_sutureid);	// must delete the cleNode first
    _sutures.erase(it);

    PublishUpdatedSutures();
}

void surgicalActions::PublishUpdatedSutures(){
    if( getServerMode() > 0 ){       
        typedef std::vector< std::pair<std::vector<int>, std::vector<double> > > suture_desc;
        suture_desc suture_updates;
        std::list<suture>::iterator it = _sutures.begin();
        while( it != _sutures.end() ){
            std::vector<int> triangles;
            triangles.push_back( it->triangle0 );
            triangles.push_back( it->triangle1 );
            std::vector<double> uvs;
            uvs.push_back( it->uv0.at(0) );
            uvs.push_back( it->uv0.at(1) );
            uvs.push_back( it->uv1.at(0) );
            uvs.push_back( it->uv1.at(1) );
            suture_updates.push_back( std::pair< std::vector<int>, std::vector<double> >(triangles, uvs) );
            it++;
        }
        getServerHandle().UpdateSutures( suture_updates );        
    }
}

//=========================================================================================
//=========================================================================================
//=========================================================================================
//
//
//
//                                  INCISION METHODS
//
//
//
//=========================================================================================
//=========================================================================================
//=========================================================================================

void surgicalActions::CreateIncision( const std::vector<int>& path_triangles, const std::vector<float>& path_uvs, 
                                      const std::vector<float>& path_positions, const std::vector<float>& path_normals, 
                                      const bool edge_start, const bool edge_end ){
    // This is Courts old incision tool glue code. It will probably need to be rewritten as we transition to a new
    // incision / cutting code.
    
    bool Tout=false,nukeThis=false;
    int n = path_triangles.size();
    skinGraphics *sg=_cle.getElasticSkinGraphics();
    trianglesUVW *tri=sg->getTrianglesUVW();
    std::vector<float> hitV,hitP;
    std::vector<int> hitTr;
    int TstartTriangleEdge=0x0003;
    float TstartParam=-1.0f;
    float uv[2],uvw[3],triangleDuv,minParam=1.0e15f;
    for(int j,i=0; i<n; ++i)	{ // nH,
        uv[0] = path_uvs[i<<1];
        uv[1] = path_uvs[(i<<1)+1];
        if(!texturePickCode(path_triangles.at(i),uv,uvw,triangleDuv)) // user picked an edge triangle
            throw std::runtime_error("Incision tool error: User picked an edge triangle");         

        if((i==0 && edge_start) )	{	// ||(i==n-1 && edge_end) removed because in T-out self intersection there is no getOppositeIncisionEdge() yet
            if(path_uvs[(i<<1)+1]==0.0f)	{
                j = 0;
                uv[0] = path_uvs[i<<1];
            }
            else if(path_uvs[i<<1]==0.0f)	{
                j = 2;
                uv[0] = 1.0f - path_uvs[(i<<1)+1];
            }
            else	{
                j = 1;
                uv[0] = path_uvs[(i<<1)+1];
            }
            TstartTriangleEdge = path_triangles[i]<<2;
            TstartTriangleEdge |= j;
            TstartParam = uv[0];
        }
    }

    if(!_cle.makeIncision(tri,n,(float (*)[3])&(path_positions[0]),(float (*)[3])&(path_normals[0]),
                          TstartTriangleEdge,TstartParam,edge_end))
        throw std::runtime_error("Incision tool error. Not recoverable");         
    
}

void surgicalActions::ExciseRegion( const int triangle ){
    _cle.exciseSubobject(triangle);
}

bool surgicalActions::texturePickCode(const int triangle, const float (&uv)[2], float (&uvw)[3], float &triangleDuv){
    // This is required by Court's old incision glue code. We may not need this as we refactor his stuff away.
    

    // texture seam triangles will have a large deltaUV in cylindrical or spherical texture mapping
	// return true=user selected a top or bottom triangle, false if an edge triangle was selected
	float *tx[3],mm[6]={1e30f,-1e30f,1e30f,-1e30f,1e30f,-1e30f},p=1.0f-uv[0]-uv[1];
	trianglesUVW *tri = _cle.getElasticSkinGraphics()->getTrianglesUVW();
	int *tr = tri->triVerts(triangle);
	if(*tr<0)
		return false;
	for(int j,i=0; i<3; ++i) {
		tx[i] = tri->vertexTexture(tr[i]);
		for(j=0; j<3; ++j) {
			if(mm[j<<1]>tx[i][j])
				mm[j<<1]=tx[i][j];
			if(mm[(j<<1)+1]<tx[i][j])
				mm[(j<<1)+1]=tx[i][j];
		}
	}
	for(int i=0; i<3; ++i)
		uvw[i] = p*tx[0][i] + uv[0]*tx[1][i] + uv[1]*tx[2][i];
	triangleDuv = mm[1]-mm[0] + mm[3]-mm[2];
	if(mm[5]-mm[4]>1e-5f)
		return false;
	if(fabs(1.0f-mm[4])<1e-5f || mm[5]<1e-5f)
		return true;
	return false;
}



//=========================================================================================
//=========================================================================================
//=========================================================================================
//
//
//
//                         OLD FUNCTIONS BELOW FOR REFERENCE
//
//
//
//=========================================================================================
//=========================================================================================
//=========================================================================================




//void surgicalActions::sendUserMessage(const char *message, const char *title, bool closeProgram)
//{
	//_frame->showModalMessageDialog(message,title,true,false);
//}

#if 0
bool surgicalActions::rightMouseDown(std::string objectHit, float (&position)[3], int triangle)
{	// returns true if a surgical action is taken, false if this is a simple viewer command
	if((_toolState==0 && objectHit[1]!='_') || (_toolState>0 && (objectHit.substr(0,2)=="H_" || objectHit.substr(0,2)=="S_")))
		return false;
	// staticTriangle objects are only scenery. If user selects one, just ignore it.
	if(dynamic_cast<staticTriangle*>(_wxg->getNodePtr(objectHit))!=NULL)
		return false;
	int hookNum;
	char s[200];
	if(_toolState==0)	//viewer
	{
		if(objectHit.substr(0,2)=="H_")	// user picked a hook
		{
			_selectedSurgObject = objectHit;
			hookNum = atoi(_selectedSurgObject.c_str()+2);
			_hooks.selectHook(hookNum);
			_sutures.selectSuture(-1);
		}
		else if(objectHit.substr(0,2)=="S_")	// user picked a suture
		{
			_selectedSurgObject = objectHit;
			hookNum = atoi(_selectedSurgObject.c_str()+2);
			_sutures.selectSuture(hookNum);
			_hooks.selectHook(-1);
		}
		else
			;
	}
	else if(_toolState==1)	// create hook mode
	{
		skinGraphics *sg = dynamic_cast<skinGraphics*>(_wxg->getNodePtr(objectHit));
		if(sg==NULL)
			return false;
		trianglesUVW *tr = sg->getTrianglesUVW();
		if(_hooks.getNumberOfHooks()<1) {	// initialize hooks
			_hooks.setHookSize(sg->getRadius()*0.03f);
			_hooks.setShapes(_wxg->getShapes());
			_hooks.setGLmatrices(_wxg->getGLmatrices());
		}
		float uv[2];
		tr->getBarycentricProjection(triangle,position[0],position[1],position[2],uv);	// int closeVert = tr->getClosestVertex(position,triangle);
		float uvw[3],triangleDuv;
		if(!texturePickCode(triangle,uv,uvw,triangleDuv)) { // user picked an edge triangle
			int ejunk;
			float pjunk;
			tr->getNearestHardEdge(position,triangle,ejunk,pjunk);
			tr->getBarycentricProjection(triangle,position[0],position[1],position[2],uv);
			if(!texturePickCode(triangle,uv,uvw,triangleDuv))
				assert(false);
		}
		if((hookNum=_hooks.addHook(tr,triangle,uv))>-1)
		{
			_hooks.selectHook(hookNum);
			_sutures.selectSuture(-1);
			_hooks.setHookCLEnode(hookNum,_cle.addHook(triangle,uv));
			sprintf(s,"H_%d",hookNum);
			_selectedSurgObject = s;
			while(_historyIt!=_historyObject.end())
				_historyIt=_historyObject.erase(_historyIt);
			json_spirit::Object hookObj;
			hookObj.push_back(json_spirit::Pair("hookNum",hookNum));
			hookObj.push_back(json_spirit::Pair("u",(double)uvw[0]));
			hookObj.push_back(json_spirit::Pair("v",(double)uvw[1]));
			hookObj.push_back(json_spirit::Pair("w",(double)uvw[2]));
			hookObj.push_back(json_spirit::Pair("triDuv",(double)triangleDuv));
//			hookObj.push_back(json_spirit::Pair("triangle",triangle));
//			hookObj.push_back(json_spirit::Pair("u",(double)uv[0]));
//			hookObj.push_back(json_spirit::Pair("v",(double)uv[1]));
			_historyObject.push_back(json_spirit::Pair("addHook",hookObj));
			_historyIt= _historyObject.end();
			//_frame->setToolState(0);
			_toolState = 0;
		}
	}
	else if(_toolState==2)	// incision mode
	{
		skinGraphics *sg = dynamic_cast<skinGraphics*>(_wxg->getNodePtr(objectHit));
		if(sg==NULL)
			return false;
		trianglesUVW *tr = sg->getTrianglesUVW();
		if(!_fence.isInitialized()) {	// initialize fence
			_fence.setFenceSize(sg->getRadius()*0.02f);
			_fence.setWxGraphics(_wxg);
		}
		bool endConn = false;
		if(_controlKeyDown || _shiftKeyDown /*|| _frame->IsCtrlOrShiftKeyDown()*/)
			endConn = true;
		float pos[3],norm[3],uv[2]={0.0f,0.0f};
		if(endConn && _fence.getPostNumber()<1)	{
			int edg;
			float param;
			tr->getNearestHardEdge(position,triangle,edg,param);
			if(edg<1)
				uv[0]=param;
			else if(edg>1)
				uv[1]=1.0f-param;
			else	{
				uv[0]=1.0f-param;
				uv[1]=param;	}
		}
		else
			tr->getBarycentricProjection(triangle,position[0],position[1],position[2],uv);
		tr->getBarycentricPosition(triangle,uv,pos);
		tr->getTriangleNormal(triangle,norm);
		_fence.addPost(tr,triangle,pos,norm,endConn);
		_hooks.selectHook(-1);
		_sutures.selectSuture(-1);
		if(_fence.getPostNumber()>1 && endConn)	// this must finish an incision
			onKeyDown(13);	// press enter key for user
	}
	else if(_toolState==3)	// create suture mode
	{
		skinGraphics *sg = dynamic_cast<skinGraphics*>(_wxg->getNodePtr(objectHit));
		trianglesUVW *tr = sg->getTrianglesUVW();
		if(sg==NULL)
			return false;
		if(_sutures.getNumberOfSutures()<1) {	// initialize sutures
			_sutures.setSutureSize(sg->getRadius()*0.01f);
			_sutures.setShapes(_wxg->getShapes());
			_sutures.setGLmatrices(_wxg->getGLmatrices());
		}
		_cle.setPhysicsPause(true);
		int i,eTri=triangle,edg;
		float param;
		tr->getNearestHardEdge(position,eTri,edg,param);
		_dragXyz[0]=position[0]; _dragXyz[1]=position[1]; _dragXyz[2]=position[2];
		i =_sutures.addSuture(tr,eTri,edg,param);
		_sutures.selectSuture(i);
		_hooks.selectHook(-1);	// deselect hooks
		_sutures.setSecondVertexPosition(i,position);
		char s[10];
		sprintf(s,"S_%d",i);
		_selectedSurgObject = s;
	}
	else if(_toolState==4)	// excise mode
	{
		while(_historyIt!=_historyObject.end())
			_historyIt=_historyObject.erase(_historyIt);
		json_spirit::Object exciseObj;
		float uv[2]={0.25f,0.25f},uvw[3],triangleDuv;
		trianglesUVW *tuvw= dynamic_cast<skinGraphics*>(_wxg->getNodePtr(objectHit))->getTrianglesUVW();
		if(!texturePickCode(triangle,uv,uvw,triangleDuv)) { // user picked an edge triangle
			int ejunk;
			float pjunk;
			tuvw->getNearestHardEdge(position,triangle,ejunk,pjunk);
			tuvw->getBarycentricProjection(triangle,position[0],position[1],position[2],uv);
			if(!texturePickCode(triangle,uv,uvw,triangleDuv))
				assert(false);
		}
		exciseObj.push_back(json_spirit::Pair("u",(double)uvw[0]));
		exciseObj.push_back(json_spirit::Pair("v",(double)uvw[1]));
		exciseObj.push_back(json_spirit::Pair("w",(double)uvw[2]));
		exciseObj.push_back(json_spirit::Pair("triDuv",(double)triangleDuv));
		_historyObject.push_back(json_spirit::Pair("exciseSubobjectTriangle",exciseObj));
		_historyIt = _historyObject.end();
		_cle.exciseSubobject(triangle);
		_hooks.selectHook(-1);
		_sutures.selectSuture(-1);
		_selectedSurgObject = "";
		//_frame->setToolState(0);
		_toolState = 0;
	}
	else
		;
	return true;
}
#endif

#if 0
bool surgicalActions::rightMouseUp(std::string objectHit, float (&position)[3], int triangle)
{
	std::string hStr;
	if(_toolState==0 && _selectedSurgObject.substr(0,2)=="H_")	// hook selected. Can only drag hooks.
	{
		Vec3f xyz;
		int hookNum = atoi(_selectedSurgObject.c_str()+2);
		_hooks.getHookPosition(hookNum,xyz._v);
		while(_historyIt!=_historyObject.end())
			_historyIt=_historyObject.erase(_historyIt);
		json_spirit::Array hArr;
		hArr.push_back(hookNum);
		hArr.push_back((double)xyz._v[0]);
		hArr.push_back((double)xyz._v[1]);
		hArr.push_back((double)xyz._v[2]);
		_historyObject.push_back(json_spirit::Pair("moveHook",hArr));
		_historyIt = _historyObject.end();
		// EFTYCHIOS - next line
#ifndef NO_PHYSICS
		_cle.setHookPosition(xyz._v,_hooks.getHookCLEnode(hookNum));
#endif

	}
	else if((_toolState==2 || _toolState==0) && _selectedSurgObject.substr(0,2)=="P_")	// fence post selected in viewer or incision mode
		_selectedSurgObject = "";
	else if(_toolState==3)	{	// finish applying a suture
		assert(_selectedSurgObject.substr(0,2)=="S_");
		trianglesUVW *tr = NULL;
		int i = atoi(_selectedSurgObject.c_str()+2);
		skinGraphics *sg;
		if(objectHit!="")
			sg = dynamic_cast<skinGraphics*>(_wxg->getNodePtr(objectHit));
		if(sg!=NULL)
			tr = sg->getTrianglesUVW();
		bool deleteSuture=false;
		int eTri0,edg0,eTri1=triangle,edg1;
		float param0,param1;
		if(tr==NULL)
			deleteSuture=true;
		else
			tr->getNearestHardEdge(position,eTri1,edg1,param1);
		if(deleteSuture)
			_sutures.deleteSuture(i);
		else	{
			if(_sutures.setSecondEdge(i,tr,eTri1,edg1,param1))	{
				_sutures.setSecondVertexPosition(i,position);
				_sutures.getEdgeAttachment(i,tr,eTri0,edg0,param0,true);
				_sutures.setCLEnode(i,_cle.addSuture(eTri0,edg0,param0,eTri1,edg1,param1));
				while(_historyIt!=_historyObject.end())
					_historyIt=_historyObject.erase(_historyIt);
				json_spirit::Object sObj;
				sObj.push_back(json_spirit::Pair("sutureNumber",i));
				json_spirit::Array sArr;
				float uv[2]={0.25f,0.25f},uvw[3],triangleDuv;
				if(!texturePickCode(eTri0,uv,uvw,triangleDuv)) { // user picked an edge triangle
					int ejunk;
					float pjunk;
					tr->getNearestHardEdge(position,triangle,ejunk,pjunk);
					tr->getBarycentricProjection(triangle,position[0],position[1],position[2],uv);
					if(!texturePickCode(triangle,uv,uvw,triangleDuv))
						assert(false);
				}
				sArr.push_back((double)uvw[0]);
				sArr.push_back((double)uvw[1]);
				sArr.push_back((double)uvw[2]);
				sArr.push_back((double)triangleDuv);
				sObj.push_back(json_spirit::Pair("point0",sArr));
				sArr.clear();
				if(!texturePickCode(eTri1,uv,uvw,triangleDuv)) { // user picked an edge triangle
					int ejunk;
					float pjunk;
					tr->getNearestHardEdge(position,triangle,ejunk,pjunk);
					tr->getBarycentricProjection(triangle,position[0],position[1],position[2],uv);
					if(!texturePickCode(triangle,uv,uvw,triangleDuv))
						assert(false);
				}
				sArr.push_back((double)uvw[0]);
				sArr.push_back((double)uvw[1]);
				sArr.push_back((double)uvw[2]);
				sArr.push_back((double)triangleDuv);
				sObj.push_back(json_spirit::Pair("point1",sArr));
				_historyObject.push_back(json_spirit::Pair("addSuture",sObj));
				_historyIt = _historyObject.end();
			}
			else
				_sutures.deleteSuture(i);
		}
		_cle.setPhysicsPause(false);
		_toolState = 0;
		//_frame->setToolState(0);
	}
	else
		return false;
	return true;
}
#endif


#if 0
bool surgicalActions::mouseMotion(float dScreenX, float dScreenY)
{
	if(_selectedSurgObject.substr(0,2)=="H_")	// hook selected.
	{
		Vec3f xyz,dragVector;
		int hookNum = atoi(_selectedSurgObject.c_str()+2);
		_hooks.getHookPosition(hookNum,xyz._v);
		_wxg->getGLmatrices()->getDragVector(dScreenX,dScreenY,xyz._v,dragVector._v);
		xyz += dragVector;
		_hooks.setHookPosition(hookNum,xyz._v);
		//xyz *= _cle.getScaleFactor(); Enabling this does really bad stuff to the new engine -- Nathan
#ifdef NO_PHYSICS	// EFTYCHIOS - this section tests skullCollision.h insideSkull()
//		if(insideSkull(xyz._v))
//			_hooks.setHookInsideSkull(hookNum);
//		else
//			_hooks.selectHook(hookNum);
#endif
	}
	else if((_toolState==2 || _toolState==0) && _selectedSurgObject.substr(0,2)=="P_")
	{
		int postNum = atoi(_selectedSurgObject.c_str()+2);
//			if(!_controlKeyDown)
//				_sog.attractFencePostNormal(hookNum,nearP._v);
//			else	// COURT - have set _dragStart[3]!!! Use it
//				_sog.attractFencePostDepth(hookNum,nearP._v);
	}
	else if(_toolState==3)	{
		assert(_selectedSurgObject.substr(0,2)=="S_");
		int sutNum = atoi(_selectedSurgObject.c_str()+2);
		Vec3f xyz,dragVector,vtx;
		_wxg->getGLmatrices()->getDragVector(dScreenX,dScreenY,_dragXyz,dragVector._v);
		_dragXyz[0]+=dragVector._v[0]; _dragXyz[1]+=dragVector._v[1]; _dragXyz[2]+=dragVector._v[2];
		const float *mm=_wxg->getGLmatrices()->getFrameAndRotationMatrix();
		transformVector3(_dragXyz,mm,xyz._v);
		xyz *= 0.7f;
		xyz._v[0]-=mm[12]; xyz._v[1]-=mm[13]; xyz._v[2]-=mm[14];
		vtx._v[0] = mm[0]*xyz._v[0] + mm[1]*xyz._v[1] + mm[2]*xyz._v[2];
		vtx._v[1] = mm[4]*xyz._v[0] + mm[5]*xyz._v[1] + mm[6]*xyz._v[2];
		vtx._v[2] = mm[8]*xyz._v[0] + mm[9]*xyz._v[1] + mm[10]*xyz._v[2];
		_sutures.setSecondVertexPosition(sutNum,vtx._v);
	}
	else
		;
	return true;
}
#endif



#if 0
void surgicalActions::onKeyDown(int key)
{
	std::string hStr;
	if(key==308 || key==17)	// ctrl key
		_controlKeyDown = true;
	else if(key==306 || key==16)	// shift key
		_shiftKeyDown = true;
	// COURT - next section is to save the current state of the .obj.  Nuke later
/*	else if(key>65 && key<91)	{ // ASCII code
		if(key==83)	{ //lock current position and save
			glslTriangle *tr = dynamic_cast<glslTriangle*>(_wxg->getNodePtr(std::string("../../scalpMaker/scalpInner")));
			int i,n;
			float vOut[3];
			GLfloat *v= tr->getPositionArray(n);
			GLfloat *mm = tr->getModelViewMatrix();
			for(i=0; i<n; ++i)	{
				transformVector3((const float (&)[3])v[i<<2],mm,vOut);
				v[i<<2]=vOut[0]; v[(i<<2)+1]=vOut[1]; v[(i<<2)+2]=vOut[2];
			}
			loadIdentity4x4(mm);
			tr->writeObjFile("skullTemporalis");
		}
		if(key==88)	_x+=0.00157f;
		if(key==89)	_y+=0.00157f;
		if(key==90)	_z+=0.00157f;
		if(key==85)	_x-=0.00157f;
		if(key==86)	_y-=0.00157f;
		if(key==97)	_z-=0.00157f;
		if(key==82)	_r+=0.01f;
		if(key==76)	_r-=0.01f;
		if(key==69)	_u+=0.01f;
		if(key==68)	_u-=0.01f;
		if(key==70)	_f+=0.01f;
		if(key==66)	_f-=0.01f;
		GLfloat *mm = _wxg->getNodePtr(std::string("../../scalpMaker/scalpInner"))->getModelViewMatrix();
		loadIdentity4x4(mm);
		rotateMatrix4x4(mm,'x',_x);
		rotateMatrix4x4(mm,'y',_y);
		rotateMatrix4x4(mm,'z',_z);
		translateMatrix4x4(mm,_r,_u,_f);
	}	*/
	else if(key==127)	// delete key
	{
		if(_toolState==2)	{
			_fence.deleteLastPost(); }
		if(_selectedSurgObject.substr(0,2)=="H_")
		{
			int hookNum = atoi(_selectedSurgObject.c_str()+2);
#ifndef NO_PHYSICS
			_cle.deleteHook(_hooks.getHookCLEnode(hookNum));	// must delete the cleNode first
#endif
			_hooks.deleteHook(hookNum);
			while(_historyIt!=_historyObject.end())
				_historyIt=_historyObject.erase(_historyIt);
			_historyObject.push_back(json_spirit::Pair("deleteHook",hookNum));
			_historyIt = _historyObject.end();
			//_frame->setToolState(0);
			_toolState = 0;
		}
		else if(_selectedSurgObject.substr(0,2)=="S_")
		{
			int sutNum = atoi(_selectedSurgObject.c_str()+2);
#ifndef NO_PHYSICS
			_cle.deleteSuture(_sutures.getCLEnode(sutNum));	// must delete the cleNode first
#endif
			_sutures.deleteSuture(sutNum);
			while(_historyIt!=_historyObject.end())
				_historyIt=_historyObject.erase(_historyIt);
			_historyObject.push_back(json_spirit::Pair("deleteSuture",sutNum));
			_historyIt = _historyObject.end();
			//_frame->setToolState(0);
			_toolState = 0;
		}
		else
			;
		_controlKeyDown=false; _shiftKeyDown=false;
	}
	else if(key==13 || key==73)	// <enter> key
	{
		if(_toolState==2)	//incision mode
		{
			std::vector<float> positions,normals,postUvs;
			std::vector<int> postTriangles;
			bool edgeStart,edgeEnd,Tout=false,nukeThis=false;
			int n=_fence.getPostData(positions,normals,postTriangles,postUvs,edgeStart,edgeEnd); // lowestP,
			while(_historyIt!=_historyObject.end())
				_historyIt=_historyObject.erase(_historyIt);
			// first parameter is which object is being cut.  For now it is a stub as there is only one object.
			json_spirit::Object sObj;
			sObj.push_back(json_spirit::Pair("incisedObject",0));	// for now only one object incisable
			sObj.push_back(json_spirit::Pair("Tin",edgeStart));
			sObj.push_back(json_spirit::Pair("Tout",edgeEnd));
			sObj.push_back(json_spirit::Pair("pointNumber",n));
			skinGraphics *sg=_cle.getElasticSkinGraphics();
			trianglesUVW *tri=sg->getTrianglesUVW();
			std::vector<float> hitV,hitP;
			std::vector<int> hitTr;
			int TstartTriangleEdge=0x0003;
			float TstartParam=-1.0f;
			float uv[2],uvw[3],triangleDuv,minParam=1.0e15f;
			json_spirit::Array iArr;
			for(int j,i=0; i<n; ++i)	{ // nH,
				uv[0] = postUvs[i<<1];
				uv[1] = postUvs[(i<<1)+1];

			if(!texturePickCode(postTriangles[i],uv,uvw,triangleDuv)) // user picked an edge triangle
				assert(false);
				iArr.push_back((double)uvw[0]);
				iArr.push_back((double)uvw[1]);
				iArr.push_back((double)uvw[2]);
				iArr.push_back((double)triangleDuv);
//				iArr.push_back(postTriangles[i]);
//				iArr.push_back((double)postUvs[i<<1]);
//				iArr.push_back((double)postUvs[(i<<1)+1]);
				if((i==0 && edgeStart) )	{	// ||(i==n-1 && edgeEnd) removed because in T-out self intersection there is no getOppositeIncisionEdge() yet
					if(postUvs[(i<<1)+1]==0.0f)	{
						j = 0;
						uv[0] = postUvs[i<<1];
					}
					else if(postUvs[i<<1]==0.0f)	{
						j = 2;
						uv[0] = 1.0f - postUvs[(i<<1)+1];
					}
					else	{
						j = 1;
						uv[0] = postUvs[(i<<1)+1];
					}
					TstartTriangleEdge = postTriangles[i]<<2;
					TstartTriangleEdge |= j;
					TstartParam = uv[0];
				}
				else	{
/*					nrm[0]=-normals[i*3]; nrm[1]=-normals[i*3+1]; nrm[2]=-normals[i*3+2];
					nH=tri->linePick(&(positions[i*3]),nrm,hitV,hitTr,hitP);
					if(nH<2)	{
						nukeThis=true;	break;	}
					minParam=1.0e15f;
					for(j=0; j<nH; ++j)	{
						if(fabs(hitP[j])<minParam)	{
							lowestP=j;
							minParam = fabs(hitP[j]);
						}
					}
					assert(lowestP+1<nH);
					float *tfl= &hitV[(lowestP+1)*3];
					tri->getBarycentricProjection(hitTr[lowestP+1],tfl[0],tfl[1],tfl[2],uv); */
//					iArr.push_back(hitTr[lowestP+1]);
//					iArr.push_back((double)uv[0]);
//					iArr.push_back((double)uv[1]);
				}
				sObj.push_back(json_spirit::Pair("incisionPoint",iArr));
				iArr.clear();
			}
			_historyObject.push_back(json_spirit::Pair("makeIncision",sObj));
			_historyIt = _historyObject.end();
			if(!nukeThis)	{
				if(!_cle.makeIncision(tri,n,(float (*)[3])&(positions[0]),(float (*)[3])&(normals[0]),TstartTriangleEdge,TstartParam,edgeEnd))	{
					sendUserMessage("Save this history file for debugging.", "Incision tool error-");
					nukeThis=true;
				}
			}
//			if(nukeThis)
//				_hst.truncateFromLastCommand();
			_fence.clear();
			//_frame->setToolState(0);
			_toolState = 0;
		}
		_controlKeyDown=false; _shiftKeyDown=false;
	}
	else
		;
}
#endif

#if 0
void surgicalActions::onKeyUp(int key)
{
	if(key==308 || key==17)	// ctrl key
		_controlKeyDown = false;
	else if(key==306 || key==16)	// shift key
		_shiftKeyDown = false;
	else	;
}
#endif

#if 0
bool surgicalActions::nextHistoryAction()
{

	if(_historyObject.empty())
	{
		wxCommandEvent junkCmd;
		//_frame->LoadHistoryFile(junkCmd);
		return false;
	}
	if(_historyIt==_historyObject.end())
	{
	    if(_cle.getServerMode()<2){
            std::cout<<"There are no more actions found in this history file" << std::endl;
        }
	    return false;
	}

    std::cout << "Executing next history action: " << _historyIt->name_ << std::endl;

	if(_historyIt->name_=="loadSceneFile")
	{
		const json_spirit::Array& fArr = _historyIt->value_.get_array();
		if(!loadScene(fArr[0].get_str().c_str(),fArr[1].get_str().c_str()))	{
		    if(_cle.getServerMode()<2)
                //_frame->showModalMessageDialog("The module file in the history file can't be found-","loadModule ERROR",true,false);
		    _historyObject.clear();
		}
	}
	else if(_historyIt->name_=="addHook")
	{
		skinGraphics *sg = _cle.getElasticSkinGraphics();
		trianglesUVW *tr=sg->getTrianglesUVW();
		if(tr==NULL)
			return false;
		if(_hooks.getNumberOfHooks()<1) {	// initialize hooks
			_hooks.setHookSize(sg->getRadius()*0.1f);
			_hooks.setShapes(_wxg->getShapes());
			_hooks.setGLmatrices(_wxg->getGLmatrices());
		}
		json_spirit::Object hookObj=_historyIt->value_.get_obj();
		int hookNum,triangle;
		float uv[2],uvw[3],dUv;
		assert(hookObj[0].name_=="hookNum");
		assert(hookObj[1].name_=="u");
		uvw[0] = hookObj[1].value_.get_real();
		assert(hookObj[2].name_=="v");
		uvw[1] = hookObj[2].value_.get_real();
		assert(hookObj[3].name_=="w");
		uvw[2] = hookObj[3].value_.get_real();
		assert(hookObj[4].name_=="triDuv");
		dUv = hookObj[4].value_.get_real();
		if(!closestTexturePick(uvw,dUv,triangle,uv))
			assert(false);
		if((hookNum=_hooks.addHook(tr,triangle,uv))>-1)
		{
            std::cout << hookNum << "  " << hookObj[0].value_.get_int() << std::endl;
			assert(hookNum == hookObj[0].value_.get_int());
			_hooks.selectHook(hookNum);
			_sutures.selectSuture(-1);
			_hooks.setHookCLEnode(hookNum,_cle.addHook(triangle,uv));
			char s[20];
#ifdef _WINDOWS
			sprintf_s(s,19,"H_%d",hookNum);
#else
			sprintf(s,"H_%d",hookNum);
#endif
			_selectedSurgObject = s;
			//_frame->setToolState(0);
			_toolState = 0;
		}
	}
	else if(_historyIt->name_=="moveHook")
	{
		Vec3f xyz;
		int hookNum;
		const json_spirit::Array& hArr = _historyIt->value_.get_array();
		hookNum = hArr[0].get_int();
		xyz._v[0] = hArr[1].get_real();
		xyz._v[1] = hArr[2].get_real();
		xyz._v[2] = hArr[3].get_real();
		_hooks.setHookPosition(hookNum,xyz._v);
		_cle.setHookPosition(xyz._v,_hooks.getHookCLEnode(hookNum));
		char s[20];
#ifdef _WINDOWS
		sprintf_s(s,19,"H_%d",hookNum);
#else
		sprintf(s,"H_%d",hookNum);
#endif
		_selectedSurgObject = s;
		//_frame->setToolState(0);
		_toolState = 0;
	}
	else if(_historyIt->name_=="addSuture")
	{
		skinGraphics *sg=_cle.getElasticSkinGraphics();
		trianglesUVW *tr=sg->getTrianglesUVW();
		int i,sutureNum,tr1,tr0,edg0,edg1;
		float uv[2],uvw0[3],uvw1[3],param0,param1,dUv0,dUv1;
		const json_spirit::Object& sObj = _historyIt->value_.get_obj();
		assert(sObj[0].name_=="sutureNumber");
		sutureNum = sObj[0].value_.get_int();
		assert(sObj[1].name_=="point0");
		uvw0[0] = sObj[1].value_.get_array()[0].get_real();
		uvw0[1] = sObj[1].value_.get_array()[1].get_real();
		uvw0[2] = sObj[1].value_.get_array()[2].get_real();
		dUv0 = sObj[1].value_.get_array()[3].get_real();
		assert(sObj[2].name_=="point1");
		uvw1[0] = sObj[2].value_.get_array()[0].get_real();
		uvw1[1] = sObj[2].value_.get_array()[1].get_real();
		uvw1[2] = sObj[2].value_.get_array()[2].get_real();
		dUv1 = sObj[2].value_.get_array()[3].get_real();
		float xyz[3];
		if(!closestTexturePick(uvw0,dUv0,tr0,uv))
			assert(false);
		tr->getBarycentricPosition(tr0,uv,xyz);
		tr->getNearestHardEdge(xyz,tr0,edg0,param0,(int)(uvw0[2]+0.1f));
		if(!closestTexturePick(uvw1,dUv1,tr1,uv))
			assert(false);
		tr->getBarycentricPosition(tr1,uv,xyz);
		tr->getNearestHardEdge(xyz,tr1,edg1,param1,(int)(uvw1[2]+0.1f));
		if(_sutures.getNumberOfSutures()<1) {	// initialize sutures
			_sutures.setSutureSize(sg->getRadius()*0.03f);
			_sutures.setShapes(_wxg->getShapes());
			_sutures.setGLmatrices(_wxg->getGLmatrices());
		}
		i =_sutures.addSuture(tr,tr0,edg0,param0);
		assert(i==sutureNum);
		_sutures.selectSuture(i);
		_hooks.selectHook(-1);	// deselect hooks
		if(_sutures.setSecondEdge(i,tr,tr1,edg1,param1))	{
			GLfloat *v,*v2;
			int *tri;
			tri = tr->triangleVertices(tr1);
			v = tr->vertexCoordinate(tri[edg1]);
			v2 = tr->vertexCoordinate(tri[(edg1+1)%3]);
			Vec3f v1;
			for(int k=0; k<3; ++k)
				v1._v[k]=v2[k]*param1 + v[k]*(1.0f-param1);
			_sutures.setSecondVertexPosition(i,v1._v);
			_sutures.setCLEnode(i,_cle.addSuture(tr0,edg0,param0,tr1,edg1,param1));	}
		_toolState = 0;
		//_frame->setToolState(0);
		char s[20];
#ifdef _WINDOWS
		sprintf_s(s,19,"S_%d",sutureNum);
#else
		sprintf(s,"S_%d",sutureNum);
#endif
		_selectedSurgObject = s;
		//_frame->setToolState(0);
		_toolState = 0;
	}
	else if(_historyIt->name_=="deleteHook")
	{
		int hookNum = _historyIt->value_.get_int();
		_cle.deleteHook(_hooks.getHookCLEnode(hookNum));
		_hooks.deleteHook(hookNum);
		//_frame->setToolState(0);
		_toolState = 0;
	}
	else if(_historyIt->name_=="deleteSuture")
	{
		int sutNum = _historyIt->value_.get_int();
		_cle.deleteSuture(_sutures.getCLEnode(sutNum));
		_sutures.deleteSuture(sutNum);
		//_frame->setToolState(0);
		_toolState = 0;
	}
	else if(_historyIt->name_=="makeIncision")
	{
		// for now ignore parsedLine[1] as there is only one incisable object.  May change this later.
		int i,incisPointNum;
		bool startIncis,endIncis;
		int TstartTriangleEdge=0x0003;
		float TstartParam=-1.0f;
		const json_spirit::Object& sObj = _historyIt->value_.get_obj();
		assert(sObj[0].name_=="incisedObject" && sObj[0].value_.get_int()==0);	// for now only one object incisable
		assert(sObj[1].name_=="Tin");
		startIncis = sObj[1].value_.get_bool();
		assert(sObj[2].name_=="Tout");
		endIncis = sObj[2].value_.get_bool();
		assert(sObj[3].name_=="pointNumber");
		incisPointNum = sObj[3].value_.get_int();
		std::vector<float> positions,normals;
		positions.assign(incisPointNum*3,0.0f); normals.assign(incisPointNum*3,0.0f);
		trianglesUVW *tri=_cle.getElasticSkinGraphics()->getTrianglesUVW();
		int triNum;
		float uv[2],uvw[3],dUv;
		for(i=0; i<incisPointNum; ++i)
		{
			assert(sObj[i+4].name_=="incisionPoint");
			const json_spirit::Array& iArr = sObj[i+4].value_.get_array();
//			if(!_hst.getNextHistoryLine(parsedLine) || parsedLine[0]!="incisionPoint")
//			{
//				_frame->showModalMessageDialog("There is an error in this history file.  Truncating from this point forward-","SURGICAL HISTORY INFORMATION",true,false);
//				_hst.truncateFromLastCommand();
//				return;
//			}
			uvw[0] = iArr[0].get_real();
			uvw[1] = iArr[1].get_real();
			uvw[2] = iArr[2].get_real();
			dUv = iArr[3].get_real();
			if(!closestTexturePick(uvw,dUv,triNum,uv))
				assert(false);
			if(i<1 && startIncis)	{
				float xyz[3];
				int edge;
				tri->getBarycentricPosition(triNum,uv,xyz);
				tri->getNearestHardEdge(xyz,triNum,edge,TstartParam,(int)(uvw[2]+0.1f));
				TstartTriangleEdge = (triNum<<2) + edge;
			}
			tri->getTriangleNormal(triNum,(float (&)[3])normals[i*3]);
			tri->getBarycentricPosition(triNum,uv,(float (&)[3])positions[i*3]);
		}
		if(!_cle.makeIncision(_cle.getElasticSkinGraphics()->getTrianglesUVW(),positions.size()/3,(float (*)[3])(&positions[0]),(float (*)[3])(&normals[0]),TstartTriangleEdge,TstartParam,endIncis))	{
		    //if(_cle.getServerMode()<2)
            std::cout <<"Incision in history file failed. Truncating from this point forward." << std::endl;
//			_hst.truncateFromLastCommand();
		    return false;
		}
		//_frame->setToolState(0);
		_toolState = 0;
	}
	else if(_historyIt->name_=="exciseSubobjectTriangle")
	{
		json_spirit::Object excObj=_historyIt->value_.get_obj();
		float uv[2],uvw[3],dUv;
		assert(excObj[0].name_=="u");
		uvw[0] = excObj[0].value_.get_real();
		assert(excObj[1].name_=="v");
		uvw[1] = excObj[1].value_.get_real();
		assert(excObj[2].name_=="w");
		uvw[2] = excObj[2].value_.get_real();
		assert(excObj[3].name_=="triDuv");
		dUv = excObj[3].value_.get_real();
		int triangle;
		if(!closestTexturePick(uvw,dUv,triangle,uv))
				assert(false);
		_cle.exciseSubobject(triangle);
		//_frame->setToolState(0);
		_toolState = 0;
	}
	else
		int i=3;
	++_historyIt;
    return true;

    return true;
}
#endif

/*
bool surgicalActions::texturePickCode(const int triangle, const float (&uv)[2], float (&uvw)[3], float &triangleDuv)
{ // texture seam triangles will have a large deltaUV in cylindrical or spherical texture mapping
	// return true=user selected a top or bottom triangle, false if an edge triangle was selected
	float *tx[3],mm[6]={1e30f,-1e30f,1e30f,-1e30f,1e30f,-1e30f},p=1.0f-uv[0]-uv[1];
	trianglesUVW *tri = _cle.getElasticSkinGraphics()->getTrianglesUVW();
	int *tr = tri->triVerts(triangle);
	if(*tr<0)
		return false;
	for(int j,i=0; i<3; ++i) {
		tx[i] = tri->vertexTexture(tr[i]);
		for(j=0; j<3; ++j) {
			if(mm[j<<1]>tx[i][j])
				mm[j<<1]=tx[i][j];
			if(mm[(j<<1)+1]<tx[i][j])
				mm[(j<<1)+1]=tx[i][j];
		}
	}
	for(int i=0; i<3; ++i)
		uvw[i] = p*tx[0][i] + uv[0]*tx[1][i] + uv[1]*tx[2][i];
	triangleDuv = mm[1]-mm[0] + mm[3]-mm[2];
	if(mm[5]-mm[4]>1e-5f)
		return false;
	if(fabs(1.0f-mm[4])<1e-5f || mm[5]<1e-5f)
		return true;
	return false;
    }*/

 /*
bool surgicalActions::closestTexturePick(const float (&uvw)[3], const float triangleDuv, int &triangle, float (&uv)[2])
{ // at present search limited to top and bottom triangles. May change later.
	triangle = -1;
	trianglesUVW *tri = _cle.getElasticSkinGraphics()->getTrianglesUVW();
	std::vector<int> *tArr=tri->getTriangleArray();
	int *tr,i,j,n=(int)tArr->size()/3;
	float *tx[3],minErr=1e30f,minDuv=1000.0f;
	float minimax[4];
	for(i=0; i<n; ++i)	{
		minimax[0]=1e30f; minimax[1]=-1e30f; minimax[2]=1e30f; minimax[3]=-1e30f;
		tr = tri->triangleVertices(i);

if(i==22239) // || i==3981 || i==14692 || i==42754)
	i=i;

		if(*tr<0)
			continue;
		tx[0] = tri->vertexTexture(*tr);
		if(fabs(tx[0][2]-uvw[2])>1e-5f)
			continue;
		tx[1] = tri->vertexTexture(tr[1]);
		if(fabs(tx[1][2]-uvw[2])>1e-5f)
			continue;
		tx[2] = tri->vertexTexture(tr[2]);
		if(fabs(tx[2][2]-uvw[2])>1e-5f)
			continue;
		for(j=0; j<3; ++j) {
			if(minimax[0]>tx[j][0])
				minimax[0]=tx[j][0];
			if(minimax[1]<tx[j][0])
				minimax[1]=tx[j][0];
			if(minimax[2]>tx[j][1])
				minimax[2]=tx[j][1];
			if(minimax[3]<tx[j][1])
				minimax[3]=tx[j][1];
		}
		if(uvw[0]+1e-5f<minimax[0] || uvw[0]-1e-5f>minimax[1])
			continue;
		if(uvw[1]+1e-5f<minimax[2] || uvw[1]-1e-5f>minimax[3])
			continue;
//		if(fabs(minimax[1]-minimax[0]+minimax[3]-minimax[2]-dUV)>0.1f) // processes texture seam triangles correctly
//			continue;
		float err=0.0f,u,v,det=(tx[1][0]-tx[0][0])*(tx[2][1]-tx[0][1]) - (tx[1][1]-tx[0][1])*(tx[2][0]-tx[0][0]);
		if(fabs(det)<1e-16f)
			continue;
		u = (uvw[0]-tx[0][0])*(tx[2][1]-tx[0][1]) - (uvw[1]-tx[0][1])*(tx[2][0]-tx[0][0]);
		u /= det;
		v =(tx[1][0]-tx[0][0])*(uvw[1]-tx[0][1]) - (tx[1][1]-tx[0][1])*(uvw[0]-tx[0][0]);
		v /= det;
		if(u<-1e-4f)
			err += u*u;
		else if(u>1.0001f)
			err += (u-1.0f)*(u-1.0f);
		else ;
		if(v<-1e-4f)
			err += v*v;
		else if(v>1.0001f)
			err += (v-1.0f)*(v-1.0f);
		else ;
		if(err==0.0f) {
			if(u+v<1.0001f) {
				det = fabs(minimax[1]-minimax[0]+minimax[3]-minimax[2]-triangleDuv);
				if(det<minDuv) {
					minDuv = det;
					triangle = i;
					uv[0] = u;
					uv[1] = v;
				}
			}
			else {
				err = 1.0001f - u - v;
				err *= err;
			}
		}
		if(err<minErr && minDuv>100.0f) {
			minErr = err;
			triangle = i;
			uv[0] = u;
			uv[1] = v;
		}
	}
	if(minDuv>100.0f)
		return false;
	else
		return true;
}
 */
