// Dependencies
var request = require('request');
var async = require('async');
var LegacyCutter = require('legacy_cutter');
var CLE = require('CLEjs');
var TriangulatedSurface = require('./surface_object.js')
var fs = require('fs');
// Code 

function SBSS_Simulation(){
    
}

SBSS_Simulation.prototype.Initialize = function(options){
    var self = this
    self.server = options.server;
    self.cutter = new LegacyCutter.Legacy_Cutter()
    self.physicsActive = false;
    self.cle = new CLE.CLEjs();
    self.hookStiffness = 1e4;
    self.sutureStiffness = 1e5;
    self.poissons_ratio = .45;
    self.youngs_modulus = 1e3;
    self.dx = 0.02;
    self.refinement = 10;
    self.update_handle = undefined;
}

SBSS_Simulation.prototype.StartUpdates = function( ){
    var self= this;
    if( self.update_handle === undefined )
        self.update_handle = setInterval( self.Update.bind( this ), 33 /* 30fps */ );
}

SBSS_Simulation.prototype.StopUpdates = function( ){
    var self = this;
    if(self.update_handle != undefined)
        clearInterval( self.update_handle )
    self.update_handle = undefined;
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
            this._LoadScene_Phase2( scene, object_fetch_list, results, fail_callback );
        }
    }.bind(this));

}

SBSS_Simulation.prototype._LoadScene_Phase2 = function( scene, model_list, model_data, callback ){
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
        sm.texture = scene["fixedObjects"][model_list[model_d]].textureMap;
        sm.normal = scene["fixedObjects"][model_list[model_d]].normalMap;
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
                                     elem.texture,
                                     elem.normal );
    });

    self.server.UpdateData();    
}

SBSS_Simulation.prototype.Incise = function( path_triangles, path_uv, path_positions, path_normals, T_in, T_out ){
    var self = this;

    self.cutter.Incise( path_triangles, path_uv, path_positions, path_normals, T_in, T_out );
    self.server.UpdateDynamicData(Array.prototype.slice.call(new Uint32Array(self.cutter.GetJS_Topology().buffer)),
                                  Array.prototype.slice.call(new Float32Array(self.cutter.GetJS_Vertex().buffer)), 
                                  Array.prototype.slice.call(new Float32Array(self.cutter.GetJS_UV().buffer)))
    self.server.UpdateData();
}

SBSS_Simulation.prototype.Excise = function( triangle ){
    var self=this
    self.cutter.Excise( triangle );
    self.server.UpdateDynamicData(Array.prototype.slice.call(new Uint32Array(self.cutter.GetJS_Topology().buffer)),
                                  Array.prototype.slice.call(new Float32Array(self.cutter.GetJS_Vertex().buffer)), 
                                  Array.prototype.slice.call(new Float32Array(self.cutter.GetJS_UV().buffer)))
    self.server.UpdateData();
}

SBSS_Simulation.prototype.AddHook = function( triangle, hook_coords ){
    var self=this

    if(!self.physicsActive)
        self.ReinitializeAndStartPhysics()

    self.cle.Add_Hook( triangle, hook_coords.slice() );
    self.UpdateHooks();
}

SBSS_Simulation.prototype.MoveHook = function( hook_id, location ){
    var self=this;

    if(!self.physicsActive)
        return;

    self.cle.Move_Hook( hook_id, location );
    self.UpdateHooks();
}

SBSS_Simulation.prototype.DeleteHook = function( hook_id){
    var self=this;

    if(!self.physicsActive)
        return;

    self.cle.Delete_Hook( hook_id );
    self.UpdateHooks();
}

SBSS_Simulation.prototype.UpdateHooks = function() {
    var self=this;
    var hook_ids = self.cle.Get_Active_Hooks();
    var hook_positions = new Float32Array( self.cle.Get_Hook_Position(hook_ids).buffer )
    var hook_triangles = new Int32Array( self.cle.Get_Hook_Triangles(hook_ids).buffer )
    hook_ids = new Uint32Array(hook_ids.buffer)
   
    hooks = [];   
    for( var i=0; i < hook_ids.length; i++ ){
        if( hook_triangles[i] != -1 )
            hooks.push( {"hook": hook_ids[i], "triangle": hook_triangles[i], "position": Array.prototype.slice.call( hook_positions, i*3, i*3+3 ) } );
    }
        
    self.server.UpdateHooks( hooks );
    self.server.UpdateData();
}

SBSS_Simulation.prototype.ReinitializeAndStartPhysics = function() {
    var self=this;

    if(self.physicsActive)
	self.cle.Destroy_Model();	// no need to destroy a non-existent model

    var vertices = self.cutter.GetRaw_Vertex();
    var topology = self.cutter.GetRaw_Topology();

    
    var verticesOut = fs.createWriteStream('model.vertices.bin');
    verticesOut.write( Buffer.from( vertices.buffer ) );
    verticesOut.end();

    var topologyOut = fs.createWriteStream('model.topology.bin');
    topologyOut.write( Buffer.from( topology.buffer ) );
    topologyOut.end();
    
    self.cle.Set_Hook_Stiffness(self.hookStiffness);
    self.cle.Set_Suture_Stiffness(self.sutureStiffness);
    self.cle.Set_Poissons_Ratio(self.poissons_ratio);
    self.cle.Set_Youngs_Modulus(self.youngs_modulus);
    self.cle.Create_Model(vertices,topology,self.dx,self.refinement);    

    // std::list<staticTriangle*>::iterator tit;
    // for(tit=_triList.begin(); tit!=_triList.end(); ++tit)	{
    //     std::vector<int> stris;
    //     std::vector<float> sverts;
    //     (*tit)->getTriangulatedSurface(sverts,stris);           
    //     self.cle.Add_Static_Model(sverts,stris);
    // }

    // // INDEPENDENT fixed geometry input with scene load. Not necessarily related to an existing
    // // vertex, line or triangle in the scene. Already computed, just output it.
    // self.cle.Set_Fixed_Geometry(_fixedVerts,_fixedPoints,_fixedSegments,_fixedTriangles);

    // for( int i = 0; i < _fixedRegions.size(); i++ ){
    //     self.cle.Set_Fixed_Volume( _fixedRegions[i].first, _fixedRegions[i].second );
    // }

    // if( _fixedmeshVerts.size() > 0 ){
    //     self.cle.Set_Fixed_Triangles(_fixedmeshVerts,_fixedmeshTriangles);  
    // }

    // // This is for collisions
    // if(!_collisionTris.empty()){
    //     self.cle.Set_Collision_Model(_collisionVerts,_collisionTris);
    // }
    
    // for( int i = 0; i<_muscleVerts.size(); i++)
    // {
    //     self.cle.Add_Muscle_Layer( _muscleVerts[i], _muscleTris[i], _muscleFibers[i], _muscleMaxStress[i]  );
    // }

    self.cle.Finalize_Initialization();
    self.physicsActive=true;
    
}

SBSS_Simulation.prototype.Update = function() {
    var self = this;
    
    if( self.physicsActive ){
        //self.cle.Advance_One_Time_Step();


    }

}

// Exports

module.exports = SBSS_Simulation;
