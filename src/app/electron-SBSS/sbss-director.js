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

SBSS_Director.prototype.ProcessCommandMessage = function( command, data ){
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
                    self.ProcessCommandMessage(cmd.command, cmd.data)
                }
            }});
        break
    case 'historyNext':
        self.history.List();
        if( self.history.CommandsRemaining() > 0 ){
            var cmd = self.history.Next();
            self.ProcessCommandMessage(cmd.command, cmd.data)
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
                             


// Exports

module.exports = SBSS_Director
