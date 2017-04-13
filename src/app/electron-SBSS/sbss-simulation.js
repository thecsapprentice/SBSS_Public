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
    self.incisionWidth = 0.002;
    self.dx = 0.02;
    self.refinementRatio = 10;
    self.update_handle = undefined;
    self.static_models = [];
    self.fixed_geometry_models = [];
    self.collision_models = [];
    self.fixed_regions = [];
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

SBSS_Simulation.prototype.LoadScene = function(scene, fail_callback, success_callback){
    var self = this;
    console.log( "Loading scene...");

    self.hookStiffness = 1e4;
    self.sutureStiffness = 1e5;
    self.poissons_ratio = .35;
    self.youngs_modulus = 3e4;
    self.incisionWidth = 0.002;
    self.dx = 0.02;
    self.refinementRatio = 3;

    object_fetch_list = []
    for( dynobj in scene.dynamicObjects ){
        object_fetch_list.push( {'name':dynobj, 'type':"dynamic"} )
        break;
    }
    if( 'fixedObjects' in scene ){
        for( fixedobj in scene.fixedObjects ){
            object_fetch_list.push( {'name':fixedobj, 'type':"static"} )
        }        
    }
    if( 'fixedGeometry' in scene )
        if( 'mesh' in scene.fixedGeometry )
            object_fetch_list.push( {'name':scene.fixedGeometry.mesh, 'type' : "fixed" } )
                
    if( 'collisionObjects' in scene ){
        for( item in scene.collisionObjects )
            object_fetch_list.push( {'name':scene.collisionObjects[item], 'type' : "collision" } )
    }

    console.log( "Fetching models..." )

    function wrapped_model_request(item, callback){
        request( 'http://localhost:8081/model/' + item.name, function(error, response, body ){
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
            this._LoadScene_Phase2( scene, object_fetch_list, results, fail_callback, success_callback );
        }
    }.bind(this));

}

SBSS_Simulation.prototype._LoadScene_Phase2 = function( scene, model_list, model_data, fail_callback, success_callback ){
    var self = this;

    // Scene paramters
    if( 'refinementRatio' in scene )
        self.refinementRatio = scene.refinementRatio;
    if( 'incisionWidth' in scene )
        self.incisionWidth = scene.incisionWidth;
    if( 'physicsNodeSpacing' in scene )
        self.dx    = scene.physicsNodeSpacing;
    if( 'hookStiffness' in scene )
        self.hookStiffness = scene.hookStiffness;
    if( 'sutureStiffness' in scene )
        self.sutureStiffness = scene.sutureStiffness;
    
    


    // Textures
    for( texture in scene.textureFiles ){
        console.log( "Loading texture", texture, "into slot", scene.textureFiles[texture] );
        self.server.RegisterTextureResource( scene.textureFiles[texture], texture );
    }
          
  
    var dynamic_model_index = -1;
    for( item in model_list ){
        if( model_list[item].type == 'dynamic' ){
            dynamic_model_index = item;
            break;
        }
    }
    if( dynamic_model_index == -1 )
        fail_callback( "Failed to parse and load dynamic model: No dynamic model found." );

    self.cutter.ParseFile(model_data[dynamic_model_index], self.incisionWidth);
    
    var dyn_texture;
    var dyn_normals;

    for( dynobj in scene.dynamicObjects ){
        dyn_texture = scene.dynamicObjects[dynobj].textureMap
        dyn_normals = scene.dynamicObjects[dynobj].normalMap
        break;
    }
         
    self.static_models.length = 0;
    self.fixed_geometry_models.length = 0;
    self.collision_models.length = 0;
    for( item in model_list ){
        if( model_list[item].type == 'static' ){
            var sm = self.cutter.ParseStaticFile( model_data[item] );
            console.log( "Loaded Static model." )
            sm.texture = scene["fixedObjects"][model_list[item].name].textureMap;
            sm.normal = scene["fixedObjects"][model_list[item].name].normalMap;
            self.static_models.push( sm );
        }
        if( model_list[item].type == 'fixed' ){
            var fm = self.cutter.ParseStaticFile( model_data[item] );
            console.log( "Loaded Fixing model." )
            self.fixed_geometry_models.push( fm );
        }
        if( model_list[item].type == 'collision' ){
            var cm = self.cutter.ParseStaticFile( model_data[item] );
            console.log( "Loaded Collision model." )
            self.collision_models.push( cm );
        }
    }   

    self.fixed_regions.length = 0;
    if( 'fixedGeometry' in scene )
        if( 'regions' in scene.fixedGeometry ) {
            for( region in scene.fixedGeometry.regions ){
                self.fixed_regions.push( [ scene.fixedGeometry.regions[region][0], scene.fixedGeometry.regions[region][1] ]  );
            }
        }
    
    // Update Dynamic Model
    
    // This is really inefficient - we should move onto sending binary data directly to the client
    self.server.UpdateDynamicData(Array.prototype.slice.call(new Uint32Array(self.cutter.GetJS_Topology().buffer)),
                                  Array.prototype.slice.call(new Float32Array(self.cutter.GetJS_Vertex().buffer)), 
                                  Array.prototype.slice.call(new Float32Array(self.cutter.GetJS_UV().buffer)))
    self.server.SetTextureName(dyn_texture)
    self.server.SetNormalName(dyn_normals)

    // Update Static Models
    self.static_models.forEach( function(elem, index){
        self.server.AddStaticObject( Array.prototype.slice.call(new Uint32Array(elem.topology.buffer)),
                                     Array.prototype.slice.call(new Float32Array(elem.vertices.buffer)),
                                     Array.prototype.slice.call(new Float32Array(elem.uv.buffer)),
                                     elem.texture,
                                     elem.normal );
    });

    self.server.UpdateData();
    success_callback();
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

SBSS_Simulation.prototype.AddSuture = function( triangleA, suture_uv_A, triangleB, suture_uv_B ){
    var self=this

    if(!self.physicsActive)
        self.ReinitializeAndStartPhysics()

    self.cle.Add_Suture( triangleA, suture_uv_A, triangleB, suture_uv_B );
    self.UpdateSutures();
}

SBSS_Simulation.prototype.DeleteSuture = function( suture_id ){
    var self=this;

    if(!self.physicsActive)
        return;

    self.cle.Delete_Suture( suture_id );
    self.UpdateSutures();
}

SBSS_Simulation.prototype.UpdateSutures = function() {
    var self=this;
    var suture_ids = self.cle.Get_Active_Sutures();
    var suture_uvs = new Float32Array( self.cle.Get_Suture_UV(suture_ids).buffer )
    var suture_triangles = new Int32Array( self.cle.Get_Suture_Triangles(suture_ids).buffer )
    suture_ids = new Uint32Array(suture_ids.buffer)
    console.log( suture_triangles.length )

    sutures = [];   
    for( var i=0; i < suture_ids.length; i++ ){
        console.log( suture_triangles[i*2+0] )
        console.log( suture_triangles[i*2+1] )
        if( suture_triangles[i*2+0] != -1 && suture_triangles[i*2+1] != -1 )
            sutures.push( {"suture": suture_ids[i],
                           "triangleA": suture_triangles[i*2+0], "uvA": Array.prototype.slice.call( suture_uvs, i*6+0, i*6+3 ),
                           "triangleB": suture_triangles[i*2+1], "uvB": Array.prototype.slice.call( suture_uvs, i*6+3, i*6+6 )
                          });
    }
        
    self.server.UpdateSutures( sutures );
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
    var refinement = Math.ceil( self.dx / (self.incisionWidth / self.refinementRatio) );
    self.cle.Create_Model(vertices,topology,self.dx,refinement);    

    // std::list<staticTriangle*>::iterator tit;
    // for(tit=_triList.begin(); tit!=_triList.end(); ++tit)	{
    //     std::vector<int> stris;
    //     std::vector<float> sverts;
    //     (*tit)->getTriangulatedSurface(sverts,stris);           
    //     self.cle.Add_Static_Model(sverts,stris);
    // }

    // Enable fixed geometry
    //self.cle.Set_Fixed_Geometry(_fixedVerts,_fixedPoints,_fixedSegments,_fixedTriangles);

    for( item in self.fixed_regions ){
        var region = self.fixed_regions[item];
        self.cle.Set_Fixed_Volume( region[0], region[1] );
    }

    for( item in self.fixed_geometry_models ){
        var fm = self.fixed_geometry_models[item];
        self.cle.Set_Fixed_Triangles(fm.vertices, fm.topology);  
    }

    for( item in self.collision_models ){
        var cm = self.collision_models[item];
        self.cle.Set_Collision_Model(cm.vertices, cm.topology);   
        break; // Only ever have one collision model.
    }
    
    // for( int i = 0; i<_muscleVerts.size(); i++)
    // {
    //     self.cle.Add_Muscle_Layer( _muscleVerts[i], _muscleTris[i], _muscleFibers[i], _muscleMaxStress[i]  );
    // }

    self.cle.Finalize_Initialization();
    self.physicsActive=true;
    self.frame = 0;
    
}

SBSS_Simulation.prototype.Update = function() {
    var self = this;
    
    if( self.physicsActive ){
        self.cle.Advance_One_Time_Step();

        var data = self.cle.Get_Vertex_Data(self.cle.V_DATA_POSITION, self.frame );

        if( self.frame != data.frame ){
            //console.log( "New frame data ready. Updating." );
            var splitVertexData = self.cutter.RemapVertexData(data.vdata);
            self.server.UpdateVertices(Array.prototype.slice.call(new Float32Array(splitVertexData.buffer)));
            self.server.UpdateData();
            self.frame = data.frame;
        }
    }



}

SBSS_Simulation.prototype.Teardown = function() {
    console.log( "Simulation Teardown." )
    var self = this;
    self.StopUpdates();
    self.physicsActive = false;
    self.cle.Destroy_Model();	// no need to destroy a non-existent model
    self.frame = 0;
}

// Exports

module.exports = SBSS_Simulation;
