// Dependenices
var request = require('request');
var SBSS_History = require('./sbss-history.js')
// Code

function SBSS_Director(){    
    this.history = new SBSS_History();
}


SBSS_Director.prototype.Initialize = function(options){
    var self = this
    self.history.Clear();
    self.simulation = options.simulation;
    self.server = options.server;
}

SBSS_Director.prototype.ProcessCommandMessage = function( command, data, historyAppend ){
    var self = this
    console.log( "Handling command '%s'...", command );
    console.log( "With data:")
    console.log( data )
    
    switch( command ){
    case 'loadHistory':
        request('http://localhost:8081/history/' + data.name, function (error, response, body) {
            if (!error && response.statusCode == 200) {
                history_cmds = JSON.parse( body );
                self.history.Clear();
                for( command in history_cmds ){
                    self.history.Push( history_cmds[command] )
                }
                if( self.history.CommandsRemaining() > 0 ){
                    var cmd = self.history.Next();
                    self.ProcessCommandMessage(cmd.command, cmd.data, false)
                }
            }});
        break
    case 'historyNext':
        self.history.List();
        if( self.history.CommandsRemaining() > 0 ){
            var cmd = self.history.Next();
            self.ProcessCommandMessage(cmd.command, cmd.data, false)
        }
        else
            console.log( "History Advanced, but no commands remaining.")
        break;
    case 'loadScene':        
        request('http://localhost:8081/scene/' + data.name, function (error, response, body) {
            if (!error && response.statusCode == 200) {
                scene_desc = JSON.parse(body)
                self.LoadScene(scene_desc);
            }});
        break;
    case 'makeIncision':
        if( historyAppend ) self.history.TruncateAppend( {'command':command, 'data':data} );
        self.MakeIncision( data );
        break;
    case 'exciseRegion':
        if( historyAppend ) self.history.TruncateAppend( {'command':command, 'data':data} );
        self.MakeExcision( data );
        break;
    }
}

SBSS_Director.prototype.LoadScene = function(scene){
    var self=this;
    console.log( scene );
    self.server.ResetConnections();
    self.simulation.LoadScene( scene, function(error){
        if( error ){
            console.log( "Load Scene failed!" );
            console.log( error );                    
        }});
}
                         

SBSS_Director.prototype.MakeIncision = function(data){
    var self=this;
    var path_triangles = []
    var path_uv = []
    var path_positions = []
    var path_normals = []
    
    for( pathNode in data.path ){
        path_triangles.push( data.path[pathNode].triangle );
        path_uv.push( data.path[pathNode].coords[0] );
        path_uv.push( data.path[pathNode].coords[1] );
        path_positions.push( data.path[pathNode].position[0] );
        path_positions.push( data.path[pathNode].position[1] );
        path_positions.push( data.path[pathNode].position[2] );
        path_normals.push( data.path[pathNode].normal[0] ); 
        path_normals.push( data.path[pathNode].normal[1] ); 
        path_normals.push( data.path[pathNode].normal[2] ); 
    }

    self.simulation.Incise( path_triangles, path_uv, path_positions, path_normals, data.edge_start, data.edge_end );
}    


SBSS_Director.prototype.MakeExcision = function(data){
    var self=this
    self.simulation.Excise( data.triangle );
}

// Exports

module.exports = SBSS_Director
