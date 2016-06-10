//Dependencies
ws = require('ws')

// Code

function Object(){
    this.topology = []
    this.vertices = []
    this.normals = []
    this.uvs = []
    this.texturename = 0;
    this.normalname = 0;

    this.topology_timestamp = 0;
    this.vertices_timestamp = 0;
    this.normals_timestamp = 0;
    this.uvs_timestamp = 0;
    this.texturename_timestamp = 0;
    this.normalname_timestamp = 0;

    this.JSONify = function(since){
        object_json_data = {}

        if( this.topology_timestamp >= since )
            object_json_data.topology = this.topology.slice();
        if( this.vertices_timestamp >= since )
            object_json_data.vertices = this.vertices.slice();
        if( this.normals_timestamp >= since )
            object_json_data.normals = this.normals.slice();
        if( this.uvs_timestamp >= since )
            object_json_data.uvs = this.uvs.slice();
        if( this.texturename_timestamp >= since )
            object_json_data.texturename = this.texturename;
        if( this.normalname_timestamp >= since )
            object_json_data.normalname = this.normalname;
        
        return object_json_data;
    }

    this.maxTimestamp = function(){
        var max_timestamp = -1;
        max_timestamp = Math.max( this.topology_timestamp, max_timestamp );
        max_timestamp = Math.max( this.vertices_timestamp, max_timestamp );
        max_timestamp = Math.max( this.uvs_timestamp, max_timestamp );
        max_timestamp = Math.max( this.normals_timestamp, max_timestamp );
        max_timestamp = Math.max( this.texturename_timestamp, max_timestamp );
        max_timestamp = Math.max( this.normalname_timestamp, max_timestamp );
        return max_timestamp;
    }
}

function StateData(){
    this.dynamic = new Object();
    this.staticobjects = []
    this.textures = {}
    
    this.textures_timestamp = 0;
    
    this.maxTimestamp = function(){
        var max_timestamp = -1;
        max_timestamp = Math.max( this.dynamic.maxTimestamp(), max_timestamp );
        for( obj in this.staticobjects )
            max_timestamp = Math.max( this.staticobjects[obj].maxTimestamp(), max_timestamp );
        max_timestamp = Math.max( this.textures_timestamp, max_timestamp );
        return max_timestamp;
    }

    this.JSONify = function(since){
        json_data = {}

        console.log( "Dynamic Timestamp", this.dynamic.maxTimestamp() );
        
        json_data.dynamic = {}
        if( this.dynamic.maxTimestamp() >= since )
            json_data.dynamic = this.dynamic.JSONify( since );

        console.log( json_data );

        json_data["static"] = []
        for( obj in this.staticobjects )
            if( this.staticobjects[obj].maxTimestamp() >= since )
                json_data["static"].push( this.staticobjects[obj].JSONify(since) )
        
        if( this.textures_timestamp >= since ){
            json_data["textures"] = {}
            for( obj in this.textures )
                json_data["textures"][obj] = {'data':this.textures[obj]};
        }
        
        json_data.timestamp = this.maxTimestamp();

        return json_data;
    }



}


function SBSS_Server() {
    this.wss = undefined;
    this.gen_ConnectionID = 0;
    this.connections = {}
    this.state_data = new StateData();
    this.globaltimestamp = 0;
}

SBSS_Server.prototype.UpdateData = function() {
    var self=this;
    console.log( "Current State is at", self.state_data.maxTimestamp() );

    for( var connectionID in self.connections )
        self.UpdateClient( connectionID )
}

SBSS_Server.prototype.UpdateClient = function(connectionID) {
    var self=this;

    //Update Data on each connection if its revision is older than the current dataset
    console.log( "Connection", connectionID, "is at", self.connections[connectionID].revision );
    var connection_update = self.state_data.JSONify( self.connections[connectionID].revision )
    //console.log( connection_update );
    self.connections[connectionID].socket.send( JSON.stringify(connection_update) );
    self.connections[connectionID].revision = connection_update.timestamp;
    console.log( "Connection", connectionID, "is now at", self.connections[connectionID].revision );  
}

SBSS_Server.prototype.UpdateTopology = function(topology){
    var self=this
    self.state_data.dynamic.topology = topology.slice();
    self.state_data.dynamic.topology_timestamp = self.IncTimestamp();
}

SBSS_Server.prototype.UpdateTextureCoords = function(uv){
    var self=this
    self.state_data.dynamic.uvs = uv.slice();
    self.state_data.dynamic.uvs_timestamp = self.IncTimestamp();
}

SBSS_Server.prototype.UpdateVertices = function(vertices){
    var self=this
    self.state_data.dynamic.vertices = vertices.slice();
    self.state_data.dynamic.vertices_timestamp = self.IncTimestamp();
}

SBSS_Server.prototype.UpdateVertexData = function(vertex_data){


}

SBSS_Server.prototype.UpdateNormals = function(normals){


}

SBSS_Server.prototype.UpdateDynamicData = function(topology, vertices, uv ){
    var self=this
    self.state_data.dynamic.topology = topology.slice();
    self.state_data.dynamic.vertices = vertices.slice();
    self.state_data.dynamic.uvs = uv.slice();

    self.state_data.dynamic.topology_timestamp = self.IncTimestamp();
    self.state_data.dynamic.uvs_timestamp = self.state_data.dynamic.topology_timestamp
    self.state_data.dynamic.vertices_timestamp = self.state_data.dynamic.topology_timestamp
}

SBSS_Server.prototype.SetTextureName = function( texname ){


}

SBSS_Server.prototype.SetNormalName = function( normname ){

    
}

SBSS_Server.prototype.AddStaticObject = function( topology, vertices, normals, uv, texname, normname ){


}

SBSS_Server.prototype.RegisterTexture = function( name, data, length ){

}

SBSS_Server.prototype.RegisterTextureResource = function( name, resource ){
    var self = this;
    self.state_data.textures[name] = resource;
    self.state_data.textures_timestamp = self.IncTimestamp();
}

SBSS_Server.prototype.UpdateHooks = function( hook_data ){

}

SBSS_Server.prototype.UpdateSutures = function( suture_data ){

}

SBSS_Server.prototype.SetSimulatorState = function( state ){

}

SBSS_Server.prototype.GetCurrentRevision = function() {


}

SBSS_Server.prototype.IncTimestamp = function() {
    var self = this;
    self.globaltimestamp += 1;
    return self.globaltimestamp;
}

SBSS_Server.prototype.VerifyClient = function(info){
//    console.log( info )
    return true;
}

SBSS_Server.prototype.Initialize = function(options) {
    var self = this;
    self.globaltimestamp = 0;
    self.director = options.director;
    self.wss = ws.Server({
        port:options.port,
        verifyClient:this.VerifyClient.bind(this)
    });
    self.wss.on('connection', function(ws) {
        ws.id = self.gen_ConnectionID++;
        self.connections[ws.id] = {socket:ws, revision:-1};
        console.log( "connection opened: %s", ws.id )

        self.UpdateClient( ws.id );
        
        ws.on('message', function(message) {
            console.log( "Client", ws.id, "sent: ", message );
            try{
                message = JSON.parse( message )
            }
            catch(e){
                console.log(e)
                console.log("Not JSON command. Discarding...");
                return;                
            }
            console.log('(%s) received: %s', ws.id, message);
            if( 'command' in message && 'data' in message )
                self.director.ProcessCommandMessage( message.command, message.data );
        });
        
        ws.on('close', function() {
            console.log( 'connection closed: %s', ws.id );
            delete self.connections[ws.id]
        });
    }); 
}




// Exports
module.exports = SBSS_Server;
