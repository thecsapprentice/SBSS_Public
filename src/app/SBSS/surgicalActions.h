#ifndef __SURGICALACTIONS_H__
#define __SURGICALACTIONS_H__

#include <string>
#include <vector>
#include <map>
#include "config.h"
#include "hooks.h"
#include "sutures.h"
#include "fences.h"
//#include "json_spirit.h"	// for json_spirit.lib version
#include "json_spirit_reader_template.h"	// these 2 lines for header only version
#include "json_spirit_writer_template.h"
#include "cleScene.h"

#ifndef JSON_SPIRIT_MVALUE_ENABLED
#define JSON_SPIRIT_MVALUE_ENABLED
#endif

// forward declarations
//class MainFrame;
class wxGraphics;

class surgicalActions
{
public:
    //
    // Class Methods 
    //

	surgicalActions(const Configuration* config);
	~surgicalActions();

	void setWxGraphics(wxGraphics *wxGraphicsPtr) {_wxg=wxGraphicsPtr; _cle.setWxGraphics(wxGraphicsPtr); }
    void update();
    void pushCommand( std::string );
	void setServerMode(int mode);  // 0=Off, 1=Passive, 2=Active
    int getServerMode() {return _serverMode;};
    HTTPVIEWER::SERVER& getServerHandle() { return *http_server; };
    cleScene& getCLEScene() { return _cle; };
    std::string loadData( std::string resource, std::string resource_type );
    
private:
    const Configuration* _config;
    bool use_network_resources;
	cleScene _cle;
    HTTPVIEWER::SERVER* http_server;

	wxGraphics *_wxg;
    int _serverMode; // 0=Off, 1=Passive client, 2=Active client
    int _port;
    std::list<hook> _hooks;
    std::list<suture> _sutures;
    std::map<std::string, std::string> availableScenes;

    std::list< std::string > raw_command_queue;
    std::list< json_spirit::mObject > command_history;
    std::list< json_spirit::mObject >::iterator current_history_position;
    std::list< json_spirit::mObject >::iterator next_history_position;

    //
    // Class Methods 
    //
    
    // Utility Methods
    const json_spirit::mValue getOrRaise( const json_spirit::mObject& obj, const std::string& key);

    // Command processing Methods
    json_spirit::mObject UnpackCommand(std::string command_str);
    void checkClientCommands();
    void processClientCommand( const json_spirit::mObject commandObj );
    void processHistoryCommand( const json_spirit::mObject commandObj );
    void processCommand( const json_spirit::mObject& commandObj, bool can_invalidate_future_history );

    // Scene and Engine Methods
	void loadScene(const char *scene);
	void updatePhysics();
    void reloadEngine();    

    // History processing Methods
	void loadHistory(const char *historyFile);
	void saveHistory(const char *fileName);
	void nextHistoryAction();
    void consumeNextHistory();
    void truncateAndAppendHistory( const json_spirit::mObject& );
    void appendHistory( const json_spirit::mObject& );
    std::string dumpHistory();

    // Hook Methods
    void Add_Hook( const int triangle, const std::vector<float>& barycentric_coords );
    void Move_Hook( const int hook_id, const std::vector<float>& coords );
    void Remove_Hook( const int hook_id );
    void PublishUpdatedHooks();

    // Suture Methods
    void Add_Suture( const int triangle1, const std::vector<float>& barycentric_coords1,
                     const int triangle2, const std::vector<float>& barycentric_coords2);
    void Remove_Suture( const int suture_id );
    void PublishUpdatedSutures();

    // Incision Methods
    void CreateIncision( const std::vector<int>& path_triangles, const std::vector<float>& path_uvs, 
                         const std::vector<float>& path_positions, const std::vector<float>& path_normals, 
                         const bool edge_start, const bool edge_end );
    void ExciseRegion( const int triangle );
    bool texturePickCode(const int triangle, const float (&uv)[2], float (&uvw)[3], float &triangleDuv);


};

#endif // __SURGICALACTIONS_H__
