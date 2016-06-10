// Dependencies
var request = require('request');
var async = require('async');
var TriangulatedSurface = require('./surface_object.js')

// Code 

function SBSS_Simulation(){
    
}

SBSS_Simulation.prototype.Initialize = function(options){
    var self = this
    self.server = options.server;

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
    var dyn_texture;
    var dyn_normals;
    
    var static_models = [];
    for( model_d in model_data ){
        if(model_d == 0)
            continue // Skip! We've just done the dynamic stuff...
        var sm = new TriangulatedSurface();
        var ret = sm.LoadFromBuffer( model_data[model_d] );
        if( !ret )
        {  callback( "Failed to parse and load static model: ", model_d );
           return; }
        console.log( "Loaded Static model." )
        static_models.push( sm );
    }
    

    
    self.server.UpdateDynamicData(dynamic_model.Flatten("vertex"), 
                                  dynamic_model.Flatten("triangles"),
                                  dynamic_model.Flatten("uv"))
    self.server.UpdateData();    
}



// Exports

module.exports = SBSS_Simulation;
