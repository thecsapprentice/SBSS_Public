// Class to handle all the model components of the scene.
//
//
//
//

var async = window.nodeRequire('async');


var ModelScene = function(masterscene) {
    this.masterscene = masterscene;
    this.loader = new THREE.TextureLoader();

    this.textures = {};
    this.static_meshes = [];
    this.vertexShader = null;
    this.fragmentShader = null;
    
    this.grid_helper = null;
    this.dynamic_mesh = null;
    this.dynamic_geometry = null;
   
    this.last_update = 0;

    // Signals
    this.dynamicTopologyUpdated = new signals.Signal();
    this.dynamicVerticesUpdated = new signals.Signal();
    
    // Sharp Edges
    this.sharp_edge_hash = {}
    this.sharp_edge_aabb = null;

    // Debugging
    this.debug = false;
    this.wireframe_helper = null;
    this.sharpedge_helper = null;
    
    // Visualization
    this.data_overlay_percent = 1.0;

    // Shader
    this.attributes = {
        //                "tangent" : { type: 'v3', value: null }
    };
    this.ambVal = new THREE.Vector3(this.masterscene.ambientLight.color.r,
                                   this.masterscene.ambientLight.color.g,
                                   this.masterscene.ambientLight.color.b);

    this.AddTexture( {name:"ColorRamp", path:{data:"colorramp.png"} }, function(err, result){
        console.log( "Loaded ColorRamp." );
    });
    var dummyRGBA = new Float32Array(256 * 256 * 3);
    for(var i=0; i< 256 * 256; i++){
        // RGB from 0 to 255
        dummyRGBA[3*i] = dummyRGBA[3*i + 1] = dummyRGBA[3*i + 2] = 1.0;
    }

    this.data_texture = new THREE.DataTexture( dummyRGBA, 256, 256, THREE.RGBFormat, THREE.FloatType );
    this.data_texture.needsUpdate = true;

    this.uniforms = {
        "dirIntensity" : { type: "f", value: this.masterscene.cameraLight.intensity },
        "diffuse" : { type: "v3", value: new THREE.Vector3(0.93333333, 0.93333333, 0.93333333) },
        "ambient" : { type: "v3", value: this.ambVal },
        "specular" : { type: "v3", value: new THREE.Vector3(0.06666666, 0.06666666, 0.06666666) },
        "shininess" : { type: "f", value: 25.0 },
        "normalScale" : { type: "v2", value: new THREE.Vector2(1.0, 1.0) },
        "dataOverlayPercent" : { type: "f", value: 1.0 },
        "dataOverlayStressRange": {type: "v2", value: new THREE.Vector2(0.0, 1.0) },
        "dataOverlayStrainRange": {type: "v2", value: new THREE.Vector2(1.0, 1.0) },
        "dataOverlayUseMagnitude" : {type: "i", value: 0 },
        "dataOverlaySelector" : {type: "i", value: 0 },
        "dataColorRamp": { type: "t", value: undefined },
        "dataTexture": {type: "t", value: this.data_texture }
    };

    //this.loader.load( this.textures["ColorRamp"].texture_path, function(texture){
    //   this.uniforms["dataColorRamp"].value = texture;
    //}.bind(this));




}

ModelScene.prototype.initialize = function() {
    $.ajax({
        async: false,
        url: 'static/client/js/vertexShader.vshader',
        success: function (data) {
            this.vertexShader = data;
        }.bind(this),
        dataType: 'text'
    });
    
    $.ajax({
        async: false,
        url: 'static/client/js/dynamic_fragmentShader.fshader',
        success: function (data) {
            this.dynamic_fragmentShader = data;
        }.bind(this),
        dataType: 'text'
    });

    $.ajax({
        async: false,
        url: 'static/client/js/fragmentShader.fshader',
        success: function (data) {
            this.fragmentShader = data;
        }.bind(this),
        dataType: 'text'
    });

    this.clearScene();
    this.grid_helper = new THREE.GridHelper( 10, 10, 0x0000ff, 0x808080 );
    this.grid_helper.position.y = - 150;
    this.masterscene.scene.add( this.grid_helper );    
}

ModelScene.prototype.clearScene = function() {
    if(!(this.dynamic_geometry === null)){ 
        this.dynamic_geometry.dispose();
        if( this.debug ){
            if( this.wireframe_helper != null )
                this.masterscene.scene.remove(this.wireframe_helper)
            if( this.sharpedge_helper != null )
                this.masterscene.scene.remove(this.sharpedge_helper)
        }
        this.masterscene.scene.remove(this.dynamic_mesh);
    }
    if( !(this.grid_helper  === null ) ){
        this.masterscene.scene.remove(this.grid_helper);
        this.grid_helper = null;
    }
    this.dynamic_mesh = null;
    this.dynamic_geometry = null;
    this.sharp_edge_aabb = new jsBVH( 3, 6 );

    this.static_meshes.forEach( function( item, index, array ){
        this.masterscene.scene.remove( item );
    });
    this.static_meshes=[];
    this.textures={};
    this.last_update = 0;
}


ModelScene.prototype.processData = function(data, callback) {
    this.last_update = data.timestamp;

    //console.log("ModelScene: Processing Timestamp: " + data.timestamp );
    //console.log( data );
    
    if( "textures" in data && !(data.textures === undefined)){
        var texToLoad = []
        for( var texId in data.textures ){
            texToLoad.push( {name: texId, path: data.textures[texId] })
        }
        async.map( texToLoad, this.AddTexture.bind(this), function( err, results ){
            if( "dynamic" in data && !(data.dynamic === undefined))
                this.ProcessObject( true, data.dynamic );
            
            if( "static" in data && !(data.static === undefined))
                data.static.forEach(this.AddStaticObject.bind(this))   

            callback( true );
        }.bind(this));        
    }
    else{
        if( "dynamic" in data && !(data.dynamic === undefined))
            this.ProcessObject( true, data.dynamic );
        
        if( "static" in data && !(data.static === undefined))
            data.static.forEach(this.AddStaticObject.bind(this))           
        
        callback(true);
    }


}

ModelScene.prototype.AddStaticObject = function( data, index, array ){
    this.ProcessObject( false, data );
}

ModelScene.prototype.ProcessObject = function( isDynamic, data ){
    var newDynObject = false;
    var dynMaterial;
    var geoNow;
    var needGeoCleanup = false;
    var signal_topology = false;
    var signal_vertices = false;
    var deformation_update = false;
    if(isDynamic) {
        if( "vertices" in data && !(
            "uvs" in data && "topology" in data ) ){
            deformation_update = true;
            geoNow = this.dynamic_geometry;
        }
        if( "uvs" in data && "topology" in data ){
            needGeoCleanup = true;
            geoNow = new THREE.BufferGeometry();        
        }
    }
    else
        geoNow = new THREE.BufferGeometry();        
    geoNow.dynamic = isDynamic;
    if( "topology" in data ) {
        var indices = new Uint32Array( data.topology.length );
        for( var i = 0; i < data.topology.length; i++)
            indices[i] = data.topology[i];
        geoNow.setIndex( new THREE.BufferAttribute( indices, 1 ) );
        geoNow.elementsNeedUpdate = true;
    }
    else if(!deformation_update) {
        var indices = new Uint32Array( 0 );
        geoNow.setIndex( new THREE.BufferAttribute( indices, 1 ) );
        geoNow.elementsNeedUpdate = true;      
    }
    if( "uvs" in data ) {
        var uvs = new Float32Array( data.uvs.length );
        for( var i = 0; i <  data.uvs.length; i++)
            uvs[i] = data.uvs[i];
        geoNow.addAttribute( 'uv', new THREE.BufferAttribute( uvs, 2 ) );
        geoNow.elementsNeedUpdate = true;
    }
    else  if(!deformation_update) {
        var uvs = new Float32Array( 0 );
        geoNow.addAttribute( 'uv', new THREE.BufferAttribute( uvs, 2 ) );
        geoNow.elementsNeedUpdate = true;
    }
    if( "vertices" in data ) {
        var positions = new Float32Array( data.vertices.length );
        for( var i = 0; i <  data.vertices.length; i++)
            positions[i] = data.vertices[i];
        geoNow.addAttribute( 'position', new THREE.BufferAttribute( positions, 3 ) );
        geoNow.elementsNeedUpdate = true;
        geoNow.computeVertexNormals();
        geoNow.attributes.position.needsUpdate = true;
        geoNow.attributes.normal.needsUpdate = true;   
        signal_vertices = true;
    }
    else{
        var positions = new Float32Array(0);
        geoNow.addAttribute( 'position', new THREE.BufferAttribute( positions, 3 ) );
        geoNow.elementsNeedUpdate = true;
        geoNow.computeVertexNormals();
        geoNow.attributes.position.needsUpdate = true;
        geoNow.attributes.normal.needsUpdate = true;   
        signal_vertices = true;
    }
    if( "strain" in data && this.uniforms.dataOverlaySelector.value == 0) {        
        var strains = new Float32Array( data.strain.length );        
        for( var i = 0; i <  data.strain.length; i++)
            strains[i] = data.strain[i];
        console.log( data.strain.length );
        geoNow.addAttribute( 'data_overlay', new THREE.BufferAttribute( strains, 3 ) );
        geoNow.elementsNeedUpdate = true;
        geoNow.attributes.data_overlay.needsUpdate = true;
        signal_vertices = true;
    }
    if( "stress" in data && this.uniforms.dataOverlaySelector.value == 1 ) {
        var stresss = new Float32Array( data.stress.length );
        for( var i = 0; i <  data.stress.length; i++)
            stresss[i] = data.stress[i];
        console.log( data.strain.length );
        geoNow.addAttribute( 'data_overlay', new THREE.BufferAttribute( stresss, 3 ) );
        geoNow.elementsNeedUpdate = true;
        geoNow.attributes.data_overlay.needsUpdate = true;
        signal_vertices = true;
    }
    if( "data_texture" in data) {
        this.uniforms.dataTexture.value.image.data = new Float32Array( data.data_texture.length );
        for( var i = 0; i <  data.data_texture.length; i++)
            this.uniforms.dataTexture.value.image.data[i] =  data.data_texture[i];
        this.uniforms.dataTexture.value.image.height = Math.sqrt( (this.data_texture.image.data.length/3) );
        this.uniforms.dataTexture.value.image.width = Math.sqrt( (this.data_texture.image.data.length/3) );
        this.uniforms.dataTexture.value.needsUpdate = true;
    }      

    
    //if(isDynamic && newDynObject===false && "uvs" in data) {
    //    dynMesh = new THREE.Mesh( geoNow, dynMaterial );
    //    scene.add( dynMesh );
    //    dynGeometry = dynMesh.geometry;
    //}
    if( "texturename" in data && "normalname" in data ){
        //if(isDynamic && newDynObject===false)
        //    alert("error in ProcessObject(true)-");    
        if(data.texturename != "" && data.normalname != ""){
            console.log( data.texturename );
            console.log( data.normalname );
            
            if(isDynamic) {               
                
                //this.attributes["data_overlay"] = { type: 'v3', value: undefined }
                //this.attributes["stress"] = { type: 'v3', value: undefined }
                this.uniforms["map"] = { type: "t",
                                         value: this.textures[data.texturename].data };
                this.uniforms["normalMap"] = { type: "t",
                                               value: this.textures[data.normalname].data };

                
                
                this.uniforms.dataTexture.value = this.data_texture;
                if( dynMaterial=== undefined ){
                    dynMaterial = new THREE.ShaderMaterial( {
                        uniforms: this.uniforms,
                        vertexShader: this.vertexShader,
                        fragmentShader: this.dynamic_fragmentShader
                    });
                    dynMaterial.extensions.derivatives = true;
                }
                
                needGeoCleanup = true; 
            }
            else {

                var uniforms = {
                    "dirIntensity" : { type: "f", value: this.masterscene.cameraLight.intensity },
                    "diffuse" : { type: "v3", value: new THREE.Vector3(0.93333333, 0.93333333, 0.93333333) },
                    "ambient" : { type: "v3", value: this.ambVal },
                    "specular" : { type: "v3", value: new THREE.Vector3(0.06666666, 0.06666666, 0.06666666) },
                    "shininess" : { type: "f", value: 25.0 },
                    "normalScale" : { type: "v2", value: new THREE.Vector2(1.0, 1.0) },
                    "map" : { type: "t", value: this.textures[data.texturename].data },
                    "normalMap" : { type: "t", value: this.textures[data.normalname].data }
                };
                
                staticMaterial = new THREE.ShaderMaterial( {
                    uniforms: uniforms,
                    vertexShader: this.vertexShader,
                    fragmentShader: this.fragmentShader
                });
                staticMaterial.extensions.derivatives = true;
                var s = new THREE.Mesh( geoNow, staticMaterial);
                this.masterscene.scene.add( s );
                this.static_meshes.push( s );
            }
        }
    }

    if( needGeoCleanup ){
        if(this.dynamic_geometry){
            this.dynamic_geometry.dispose();
            if(dynMaterial === undefined)
                dynMaterial = this.dynamic_mesh.material;
            if(this.debug)
                this.masterscene.scene.remove(this.wireframe_helper)
            this.masterscene.scene.remove(this.dynamic_mesh);
        }
        else{
            geoNow.computeBoundingSphere();
            if( geoNow.boundingSphere.radius > 0 ){
                var manip_scale = geoNow.boundingSphere.radius;
                this.masterscene.toolscene.setScalingFactor( manip_scale );
                //suture_object.children[0].scale.multiplyScalar( manip_scale*0.015 );
                //suture_object.children[1].scale.set( manip_scale*0.008,  1.0, manip_scale*0.008);
                this.masterscene.camera.position = geoNow.boundingSphere.center;
                this.masterscene.camera.position.z += geoNow.boundingSphere.radius*3.2;
                this.masterscene.cameraLight.position = this.masterscene.camera.position;
                this.masterscene.camera.lookAt( geoNow.boundingSphere.center );
            }
        }
        
        this.dynamic_mesh = new THREE.Mesh( geoNow, dynMaterial );
        this.masterscene.scene.add( this.dynamic_mesh );        

        this.dynamic_geometry = this.dynamic_mesh.geometry;
        signal_topology = true;
        signal_vertices = true;
    }

    if(signal_topology){
        this.dynamicTopologyUpdated.dispatch();
        this.CollectSharpEdges();
    }
    if(signal_vertices){
        this.dynamicVerticesUpdated.dispatch();
        if(this.debug){
            if( this.wireframe_helper != null )
                this.masterscene.scene.remove( this.wireframe_helper );
            this.wireframe_helper = new THREE.WireframeHelper( this.dynamic_mesh );
            this.wireframe_helper.material.depthTest = true;
	    this.wireframe_helper.material.opacity = 1;
	    this.wireframe_helper.material.transparent = false;
            this.masterscene.scene.add( this.wireframe_helper );
            this.createSharpEdgeHelper();
        }
        this.UpdateSharpEdgeAABB();
    }

//    return 0;
}

ModelScene.prototype.AddTexture = function(data,callback){
    if( data.name in this.textures )
        return;
    var scene_texture = {
        texture_path: "http://localhost:8081/texture/"+data.path.data
    }
    
    this.loader.load( scene_texture.texture_path, function(texture){
        scene_texture.data = texture;
        this.textures[data.name] = scene_texture;
        callback( null, scene_texture );    
    }.bind(this));
   
}


ModelScene.prototype.CollectSharpEdges = function(){       
    this.sharp_edge_hash = {};

    var sortFunction = function ( a, b ) { return a < b };
    var uvs = this.dynamic_geometry.attributes[ "uv" ].array;
    var indices = this.dynamic_geometry.getIndex().array;
    var positions = this.dynamic_geometry.attributes[ "position" ].array;      

    var edge = [0,0];

    for( var i = 0; i < indices.length/3; i++){
        tri = [ indices[i*3+0], indices[i*3+1], indices[i*3+2] ];
        for( var e = 0; e < 3; e ++){
            edge[0] = tri[e];
            edge[1] = tri[ (e+1)%3 ];            
            edge.sort( sortFunction );            
            var key = edge.toString();          
            if( this.sharp_edge_hash[key] === undefined )
                this.sharp_edge_hash[key] = { vert1: edge[0], vert2: edge[1], 
                                              seg1: undefined, seg2: undefined,
                                              face1: i, face2: undefined };
            else
                this.sharp_edge_hash[key].face2 = i;
        }
    }

    this.UpdateSharpEdgeAABB();
}

ModelScene.prototype.UpdateSharpEdgeAABB = function() {
    var positions = this.dynamic_geometry.attributes[ "position" ].array;      
    var uvs = this.dynamic_geometry.attributes[ "uv" ].array;
    this.sharp_edge_aabb = new jsBVH( 3, 6 );
    var segment = [null,null];
    var uv = [null, null];

    //console.log( "Filling AABB with detected sharp edges.")
    var collected_edges = 0;
    for ( var key in this.sharp_edge_hash ) {
        var h = this.sharp_edge_hash[ key ];
        if( h.face2 === undefined ){
            segment[0] = new THREE.Vector3().fromArray(positions, h.vert1*3)
            segment[1] = new THREE.Vector3().fromArray(positions, h.vert2*3)
            uv[0] = new THREE.Vector2().fromArray(uvs, h.vert1*2)
            uv[1] = new THREE.Vector2().fromArray(uvs, h.vert2*2)

            // Skip edges which belong to side triangles
            // Identified by the UV values --- REALLY NEED BETTER WAY TO DO THIS!!!!
            if((uv[0].x > 1.05 || uv[1].x > 1.05) && ( uv[0].x < 2.15 || uv[1].x < 2.15 )){
                h.face2 = h.face1; // Cheap hack to remove this key from consideration.
                continue;
            }

            h.seg1 = segment[0];
            h.seg2 = segment[1];
            var segment_box = new THREE.Box3();
            segment_box.setFromPoints( segment );
            var size = segment_box.size();
            var min = segment_box.min;
            this.sharp_edge_aabb.insert( { intervals:
                                           [ {a:min.x, b:size.x},
                                             {a:min.y, b:size.y},
                                             {a:min.z, b:size.z} ],
                                           object: key
                                         }
                                       );
            collected_edges++;
            
        }
    }
    //console.log( "Added " + collected_edges + " sharp edges to heirarchy.");

}


ModelScene.prototype.IsEdgeSharp = function( edge ){
    var sortFunction = function ( a, b ) { return a < b };    
    var sedge = [ edge[0], edge[1] ];
    sedge.sort( sortFunction );
    var key = sedge.toString();
    var h = this.sharp_edge_hash[key];
    if( h === undefined)
        return false;

    if( h.face2 === undefined )
        return true;
    else
        return false;
}

ModelScene.prototype.ClosestPointOnSharpEdge = function( point, radius, origin ){
    var indices = this.dynamic_geometry.getIndex().array;
    var positions = this.dynamic_geometry.attributes[ "position" ].array;

    var searchbox = new THREE.Box3();
    searchbox.expandByPoint( point );
    searchbox.expandByScalar( radius );
    var size = searchbox.size();
    var min = searchbox.min;
    var nearedges = this.sharp_edge_aabb.search({intervals:[{a:min.x, b:size.x},
                                                            {a:min.y, b:size.y},
                                                            {a:min.z, b:size.z}]});
    var result = null;    
    var best_distance = 1e8;
    for( var k = 0; k < nearedges.length; k++ ){
        var key = nearedges[k];
        var h = this.sharp_edge_hash[ key ];
        var closestPoint = new THREE.Vector3();
        if( !(origin===undefined )){
            var direction = new THREE.Vector3().copy( point );
            direction.sub( origin ).normalize();
            var ray = new THREE.Ray( origin, direction );
            ray.distanceSqToSegment(h.seg1, h.seg2, null, closestPoint );
        }
        else{
            var line = new THREE.Line3( h.seg1, h.seg2 );
            line.closestPointToPoint( point, true, closestPoint );
        }
        var dist = point.distanceTo( closestPoint );
        if( dist < best_distance ){
            var tri = new THREE.Triangle( new THREE.Vector3().fromArray( positions, indices[h.face1*3+0]*3 ),
                                          new THREE.Vector3().fromArray( positions, indices[h.face1*3+1]*3 ),
                                          new THREE.Vector3().fromArray( positions, indices[h.face1*3+2]*3 ) );

            var bary = tri.barycoordFromPoint( closestPoint );
            result = {edge_point: closestPoint,
                      tri: tri,
                      triangleID: h.face1,
                      bary_coords: bary,
                      indices: [ indices[h.face1*3+0], indices[h.face1*3+1], indices[h.face1*3+2] ],
                      distance: dist
                     };
            best_distance = dist;
        }
    }
    return result;
}


ModelScene.prototype.toggleDebugMode = function( debugFlag ) {
    console.log( "Toggling debug to " + debugFlag + " state for model display." );
    this.debug = debugFlag;

    if( !this.debug ){
        if( this.wireframe_helper != null )
            this.masterscene.scene.remove( this.wireframe_helper );
        this.wireframe_helper = null;
        if( this.sharpedge_helper != null ){
            this.masterscene.scene.remove( this.sharpedge_helper );
            this.sharpedge_helper.geometry.dispose();
            this.sharpedge_helper = null;
        }
    }
    else{
        if( this.wireframe_helper != null )
            this.masterscene.scene.remove( this.wireframe_helper );
        if( this.dynamic_mesh ){
            this.wireframe_helper = new THREE.WireframeHelper( this.dynamic_mesh );
            this.wireframe_helper.material.depthTest = true;
	    this.wireframe_helper.material.opacity = 1;
	    this.wireframe_helper.material.transparent = false;
            this.masterscene.scene.add( this.wireframe_helper );
            this.createSharpEdgeHelper();
        }
    }
}

ModelScene.prototype.createSharpEdgeHelper = function( ) {
    if( this.sharpedge_helper != null ){
        this.masterscene.scene.remove( this.sharpedge_helper );
        this.sharpedge_helper.geometry.dispose();
        this.sharpedge_helper = null;
    }
    var indices = this.dynamic_geometry.getIndex().array;
    var positions = this.dynamic_geometry.attributes[ "position" ].array;

    var geometry = new THREE.BufferGeometry();
    var coords = new Float32Array( Object.keys(this.sharp_edge_hash).length * 2 * 3 );
    var index = 0;
    for ( var key in this.sharp_edge_hash ) {
        var h = this.sharp_edge_hash[ key ];
        if ( h.face2 === undefined  ) { // hardwired const OK
            coords[ index ++ ] = positions[ h.vert1*3+0 ];
	    coords[ index ++ ] = positions[ h.vert1*3+1 ];
	    coords[ index ++ ] = positions[ h.vert1*3+2 ];          
	    coords[ index ++ ] = positions[ h.vert2*3+0 ];
	    coords[ index ++ ] = positions[ h.vert2*3+1 ];
	    coords[ index ++ ] = positions[ h.vert2*3+2 ];            
	}
    }
    geometry.addAttribute( 'position', new THREE.BufferAttribute( coords, 3 ) );
    this.sharpedge_helper = new THREE.Line( geometry, new THREE.LineBasicMaterial( { color: 0xff00ff, linewidth: 3 } ), THREE.LineSegments );
    this.masterscene.scene.add( this.sharpedge_helper );
}

ModelScene.prototype.DumpModelAsObj = function (){
    if( this.dynamic_mesh == null ) 
        return;

    console.log( "Dumping model as OBJ" );

    $("#model").empty()
    
    var arr_index = this.dynamic_mesh.geometry.getIndex().array;
    var arr_normals = this.dynamic_mesh.geometry.attributes[ "normal" ].array;
    var arr_uv = this.dynamic_mesh.geometry.attributes[ "uv" ].array;
    var arr_position = this.dynamic_mesh.geometry.attributes[ "position" ].array;

    var BlobData = []

    console.log( "Dumping " + String(arr_position.length/3) + " vertices." )
    for( var v = 0; v < arr_position.length; v+=3){
        var obj_text = "v "
        for( var vp = v; vp < v+3; vp++ )
            obj_text += String(arr_position[vp]) + " "
        obj_text += "\n"
        //$("#model").append(obj_text);
        BlobData.push( obj_text );
    }

    console.log( "Dumping " + String(arr_uv.length/2) + " uvs." )
    for( var v = 0; v < arr_uv.length; v+=2){
        var obj_text = "vt "
        for( var vp = v; vp < v+2; vp++ )
            obj_text += String(arr_uv[vp]) + " "
        obj_text += "\n"
        //$("#model").append(obj_text);
        BlobData.push( obj_text );
    }

    console.log( "Dumping " + String(arr_normals.length/3) + " normals." )
    for( var v = 0; v < arr_normals.length; v+=3){
        var obj_text = "vn "
        for( var vp = v; vp < v+3; vp++ )
            obj_text += String(arr_normals[vp]) + " "
        obj_text += "\n"
        //$("#model").append(obj_text);
        BlobData.push( obj_text );
    }

    console.log( "Dumping " + String(arr_index.length/3) + " faces." )
    for( var v = 0; v < arr_index.length; v+=3){
        var obj_text = "f "
        for( var vp = v; vp < v+3; vp++ )
            obj_text += String(arr_index[vp]+1) + "/" + String(arr_index[vp]+1) + "/" + String(arr_index[vp]+1) + " "
        obj_text += "\n"
        //$("#model").append(obj_text);
        BlobData.push( obj_text );
    }

    var oBlob = new Blob(BlobData, {type : 'text/plain;charset=utf-8'}); // the blob
    saveAs( oBlob, "object_dump.obj" );
        
}


ModelScene.prototype.setDataOverlayPercent = function(newPercent) {
    this.uniforms.dataOverlayPercent.value = Math.min( Math.max(newPercent, 0.0 ) , 1.0 );
}

ModelScene.prototype.setDataOverlayStrainRange = function(rLow, rHigh) {
    this.uniforms.dataOverlayStrainRange.value = new THREE.Vector2( rLow, rHigh );
}

ModelScene.prototype.setDataOverlayStressRange = function(rLow, rHigh) {
    this.uniforms.dataOverlayStressRange.value = new THREE.Vector2( rLow, rHigh );
}

ModelScene.prototype.setDataOverlayUseMagnitude = function( toggle ){
    if( toggle )
        this.uniforms.dataOverlayUseMagnitude.value = 1;
    else
        this.uniforms.dataOverlayUseMagnitude.value = 0;
}

ModelScene.prototype.setDataOverlaySource = function( source ) {
    if( source == "strain" )
        this.uniforms.dataOverlaySelector.value = 0;
    else
        this.uniforms.dataOverlaySelector.value = 1;
}
