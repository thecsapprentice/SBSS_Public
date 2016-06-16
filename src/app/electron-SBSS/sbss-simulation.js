// Dependencies
var request = require('request');
var async = require('async');
var LegacyCutter = require('legacy_cutter');

var TriangulatedSurface = require('./surface_object.js')

// Code 

function SBSS_Simulation(){
    
}

SBSS_Simulation.prototype.Initialize = function(options){
    var self = this
    self.server = options.server;
    self.cutter = new LegacyCutter.Legacy_Cutter()
}

SBSS_Simulation.prototype.LoadScene = function(scene, fail_callback){
    var self = this;
    console.log( "Loading scene...");

    object_fetch_list = []
    for( dynobj in scene.dynamicObjects ){
        object_fetch_list.push( dynobj )
        break;
    }
    if( 'fixedObjects' in scene ){
        for( fixedobj in scene.fixedObjects ){
            object_fetch_list.push( fixedobj )
        }        
    }

    console.log( "Fetching models..." )

    function wrapped_model_request(item, callback){
        request( 'http://localhost:8081/model/' + item, function(error, response, body ){
            if (!error && response.statusCode == 200) {
                callback( null, body );
            }
            else{
                console.log( error )
                callback( "Failed to fetch file!", null );
            }
        });
    }
    console.log( fail_callback )
    async.map( object_fetch_list, wrapped_model_request, function( err, results ){
        if( err ){            
            fail_callback( "Failed to load some models. Loading scene can't continue!" )
        }
        else{
            this._LoadScene_Phase2( scene, results, fail_callback );
        }
    }.bind(this));

}

SBSS_Simulation.prototype._LoadScene_Phase2 = function( scene, model_data, callback ){
    var self = this;

    // Scene paramters
    var incisionWidth = scene.incisionWidth;
    var physicsDx     = scene.physicsNodeSpacing;
    
    // Textures
    for( texture in scene.textureFiles ){
        console.log( "Loading texture", texture, "into slot", scene.textureFiles[texture] );
        self.server.RegisterTextureResource( scene.textureFiles[texture], texture );
    }
          
    
    var dynamic_model = new TriangulatedSurface()
    var ret = dynamic_model.LoadFromBuffer( model_data[0] ); // Assume first model is the dynamic one
    if( ret )
        console.log( "Loaded Dynamic model." )
    else{
        callback( "Failed to parse and load dynamic model." );
        return;
    }
    

    self.cutter.ParseFile(model_data[0]);
    //console.log( new Float32Array(self.cutter.GetJS_Vertex().buffer) )
    
    var dyn_texture;
    var dyn_normals;

    for( dynobj in scene.dynamicObjects ){
        dyn_texture = scene.dynamicObjects[dynobj].textureMap
        dyn_normals = scene.dynamicObjects[dynobj].normalMap
        break;
    }
         
    var static_models = [];
    for( model_d in model_data ){
        if(model_d == 0)
            continue // Skip! We've just done the dynamic stuff...
        var sm = self.cutter.ParseStaticFile( model_data[model_d] );
        console.log( "Loaded Static model." )
        static_models.push( sm );
    }   
    
    // Update Dynamic Model
    
    // This is really inefficient - we should move onto sending binary data directly to the client
    self.server.UpdateDynamicData(Array.prototype.slice.call(new Uint32Array(self.cutter.GetJS_Topology().buffer)),
                                  Array.prototype.slice.call(new Float32Array(self.cutter.GetJS_Vertex().buffer)), 
                                  Array.prototype.slice.call(new Float32Array(self.cutter.GetJS_UV().buffer)))
    self.server.SetTextureName(dyn_texture)
    self.server.SetNormalName(dyn_normals)

    // Update Static Models
    static_models.forEach( function(elem, index){
        self.server.AddStaticObject( Array.prototype.slice.call(new Uint32Array(elem.topology.buffer)),
                                     Array.prototype.slice.call(new Float32Array(elem.vertices.buffer)),
                                     Array.prototype.slice.call(new Float32Array(elem.uv.buffer)),
                                     0,
                                     0 );
    });

    self.server.UpdateData();    
}



// Exports

module.exports = SBSS_Simulation;
