//Dependencies
ws = require('ws')

// Code
function SBSS_Server() {
    this.wss = undefined;
    this.gen_ConnectionID = 0;
    this.connections = {}
}

SBSS_Server.prototype.UpdateData = function() {
    var self=this;
    
    for( var connectionID in self.connections ){
        //Update Data on each connection if its revision is older than the current dataset

        /*
        if( self.connections[connectionID].revision < ??? ){


        }
        */
    }

}

SBSS_Server.prototype.UpdateTopology = function(topology){


}

SBSS_Server.prototype.UpdateTextureCoords = function(uv){


}

SBSS_Server.prototype.UpdateVertices = function(vertices){

    
}

SBSS_Server.prototype.UpdateVertexData = function(vertex_data){


}

SBSS_Server.prototype.UpdateNormals = function(normals){


}

SBSS_Server.prototype.UpdateDynamicData = function(topology, vertices, uv ){


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


}

SBSS_Server.prototype.UpdateHooks = function( hook_data ){

}

SBSS_Server.prototype.UpdateSutures = function( suture_data ){

}

SBSS_Server.prototype.SetSimulatorState = function( state ){

}

SBSS_Server.prototype.GetCurrentRevision = function() {


}

SBSS_Server.prototype.VerifyClient = function(info){
//    console.log( info )
    return true;
}

SBSS_Server.prototype.Initialize = function(options) {
    var self = this;
    self.director = options.director;
    self.wss = ws.Server({
        port:options.port,
        verifyClient:this.VerifyClient.bind(this)
    });
    self.wss.on('connection', function(ws) {
        ws.id = self.gen_ConnectionID++;
        self.connections[ws.id] = {socket:ws, revision:-1};

        console.log( "connection opened: %s", ws.id )
        
        ws.on('message', function(message) {
            message = JSON.parse( message )
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
