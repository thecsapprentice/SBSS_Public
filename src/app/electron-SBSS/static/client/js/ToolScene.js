// Class encapsulating scene information for tools
//
//
//

var ToolScene = function(masterscene){
    this.hooks = [];
    this.hook_data = [];
    this.sutures = [];
    this.suture_data = [];
    this.masterscene = masterscene;
    this.last_update = 0;

    this.tool_sub_mode = null;

    this.hook_mesh = null;
    this.hook_geometry = null;
    this.hook_material = null;

    this.suture_object = null;

    this.selected_hook = null;
    this.selected_hook_index = null;
    this.selected_hook_position = null;
    this.selected_suture = null;
    this.selected_suture_index = null;

    this.scaling_factor = 1.0;
    this.mouse = new THREE.Vector2();

    this.keyflags = { 
        shift: false,
        alt: false,
        ctrl: false,
        meta: false
    };

    this.plane = null;
    this.offset = null;
    
    this.draw_suture_group = null;
    this.draw_suture_segment = null;
    this.draw_suture_fixed = null;
    this.draw_suture_free = null;
    this.suture_placement = {
        start: null,
        end: null
    }

    this.draw_incision = null;
    this.incision_placement = [];

    //Debugging
    this.debug = false;
       
};

ToolScene.prototype.initialize = function(){

    var modifier = new THREE.SubdivisionModifier( 2 );
    
    // Build Hook Mesh
    var points = [];
    points.push( new THREE.Vector3(0,0,0) );
    points.push( new THREE.Vector3(.1,0,0) );
    points.push( new THREE.Vector3(.3,0,1) );
    points.push( new THREE.Vector3(0,0,1) );
    this.hook_geometry = new THREE.IcosahedronGeometry(1.0, 1);//new THREE.LatheGeometry( points, 8 );
    this.hook_geometry.mergeVertices();
    this.hook_geometry.computeFaceNormals();
    this.hook_geometry.computeVertexNormals();
    modifier.modify(this.hook_geometry); // Subdivide for smoothness!
    this.hook_geometry.normalsNeedUpdate = true;
    this.hook_material = new THREE.MeshLambertMaterial( { color: 0xffff00, opacity: 0.8 } );
    this.hook_selected_material = new THREE.MeshLambertMaterial( { color: 0x00ff00, opacity: 0.8 } );
    this.hook_mesh = new THREE.Mesh( this.hook_geometry, this.hook_material );

    // Build Suture Mesh
    var suture_sphere = new THREE.IcosahedronGeometry( .4, 1);
    modifier.modify(suture_sphere);
    var suture_cylinder = new THREE.CylinderGeometry(1.0, 1.0, 1.0);
    this.suture_material = new THREE.MeshLambertMaterial( { color: 0xffffee, opacity: 0.8 } );
    this.suture_selected_material = new THREE.MeshLambertMaterial( { color: 0x00ff00, opacity: 0.8 } );
    this.suture_bad_material = new THREE.MeshLambertMaterial( { color: 0xff3333, opacity: 0.8 } );

    this.suture_object = new THREE.Object3D();
    var smesh = new THREE.Mesh( suture_sphere, this.suture_material );
    this.suture_object.add(smesh);
    smesh = new THREE.Mesh( suture_sphere, this.suture_material );
    this.suture_object.add(smesh);
    smesh = new THREE.Mesh( suture_cylinder, this.suture_material );
    this.suture_object.add(smesh);

    // Hooks
    this.masterscene.renderer.domElement.addEventListener('mousemove', this.handleMouseMove.bind(this), false);
    this.masterscene.renderer.domElement.addEventListener('mousedown', this.handleMouseDown.bind(this), false);
    this.masterscene.renderer.domElement.addEventListener('mouseup', this.handleMouseUp.bind(this), false);


    // Move plane
    this.plane = new THREE.Mesh(
	new THREE.PlaneBufferGeometry( 2000, 2000, 8, 8 ),
	new THREE.MeshBasicMaterial( { wireframe:true, color: 0xFFFFFF, opacity: 0.25, transparent: false, side:THREE.DoubleSide } )
    );
    this.plane.material.visible = false;
    this.masterscene.scene.add( this.plane );
    this.offset = new THREE.Vector3(0,0,0);


    // Construct a template suture for placement purposes
    this.draw_suture_group = new THREE.Object3D();
    this.draw_suture_fixed = new THREE.Mesh( suture_sphere, this.suture_material );
    this.draw_suture_group.add( this.draw_suture_fixed );
    this.draw_suture_free = new THREE.Mesh( suture_sphere, this.suture_material );
    this.draw_suture_group.add( this.draw_suture_free );
    this.draw_suture_segment = new THREE.Mesh( suture_cylinder, this.suture_material );
    this.draw_suture_group.add( this.draw_suture_segment );
    
    this.draw_suture_group.visible = false;
    this.masterscene.scene.add( this.draw_suture_group );


    document.addEventListener('keydown', this.handleKeyDown.bind(this), false);
    document.addEventListener('keyup', this.handleKeyUp.bind(this), false);

    this.clearScene();
}

ToolScene.prototype.setScalingFactor = function(factor) {
    this.scaling_factor = factor;
}

ToolScene.prototype.clearScene = function() {
    this.last_update = 0;

    this.hooks.forEach( function( item, index, array ){
        this.masterscene.scene.remove( item );
    });
    this.hooks = [];

    this.sutures.forEach( function( item ){  // , index, array
        this.masterscene.scene.remove( item );
    });
    this.suture_data = [];
    this.sutures = [];

    this.clearMode();
}

ToolScene.prototype.processData = function(data, callback) {
    this.last_update = data.timestamp;

    //console.log("ToolScene: Processing Timestamp: " + data.timestamp );
    //console.log( data );
    
    if( "simulator_state" in data ){       
        if( data.simulator_state == 1 )
            this.lockMode();
        else{
            var mode = this.getMode();
            if( mode.primary == -1 )
                this.clearMode();
        }
    }   

    if( "hooks" in data && !(data.hooks === undefined)){
        this.clearManipulatorSelection();
        this.ProcessHooks( data.hooks );
    }

    if( "sutures" in data && !(data.sutures === undefined)){
        this.clearManipulatorSelection();
        this.ProcessSutures( data.sutures );
    }

    callback( true );
}

ToolScene.prototype.ProcessHooks = function( data ){
    this.hook_data.length = 0;

    data.forEach( function( item ){  // , index, array
        this.hook_data.push( item );
    }.bind(this)); 
    
    this.updateHooks();
}

ToolScene.prototype.ProcessSutures = function( data ){
    this.suture_data.length = 0;

    data.forEach( function( item ){  // , index, array
        this.suture_data.push( item );
    }.bind(this)); 
    
    this.updateSutures();
}


ToolScene.prototype.getIntersection = function( coord, objects ) {

    var Ex = ( coord.x / this.masterscene.renderer.domElement.clientWidth ) * 2 - 1;
    var Ey = - ( coord.y / this.masterscene.renderer.domElement.clientHeight ) * 2 + 1;
    var vector = new THREE.Vector3(Ex, Ey, 0.5).unproject( this.masterscene.camera );
    var cameraPos = this.masterscene.camera.position;
    var rayCaster = new THREE.Raycaster(cameraPos,
                                        vector.sub(cameraPos).normalize());
    var vA, vB, vC;
    var bary;
    var triId;
    var intersections = rayCaster.intersectObjects(objects, true);

    var msg = {
        position: null,
        triangleID: null,
        coord: null
    };

    if (intersections.length > 0) {        
        intersection = intersections[0];
        if( intersection.object.geometry instanceof THREE.BufferGeometry ){
            positionsArray = intersection.object.geometry.attributes.position.array;
            vA = new THREE.Vector3().fromArray(positionsArray, intersection.face.a*3)
            vB = new THREE.Vector3().fromArray(positionsArray, intersection.face.b*3)
            vC = new THREE.Vector3().fromArray(positionsArray, intersection.face.c*3)
            var tri = new THREE.Triangle( vA, vB, vC );
            triId = intersection.faceIndex;
            bary = tri.barycoordFromPoint( intersection.point );
            msg = {
                position: intersection.point,
                triangleID: triId,
                bary_coord: bary,
                object: intersection.object,
                triangle: tri,
                indices: [ intersection.face.a, intersection.face.b, intersection.face.c ]
            };
        }
        else{
            msg = {
                position: intersection.point,
                object: intersection.object
            };
        }
    }

    return msg;
}

ToolScene.prototype.getMode = function () {
    var mode = {
        primary: this.masterscene.toolsgui._toolState,
        secondary: this.tool_sub_mode
    };
    return mode;
}

ToolScene.prototype.incMode = function () {
    if( this.tool_sub_mode == null )
        this.tool_sub_mode = 0;
    this.tool_sub_mode += 1;
}

ToolScene.prototype.clearMode = function () {
    if(this.masterscene.controls)
        this.masterscene.controls.enabled=true;
    this.tool_sub_mode = null;
    this.masterscene.toolsgui.clearToolstate();
}

ToolScene.prototype.lockMode = function () {
    this.masterscene.controls.enabled=false;
    this.tool_sub_mode = null;
    this.masterscene.toolsgui.setToolstate(-1);
}

ToolScene.prototype.updateKeyFlags = function(event){
    this.keyflags.shift = event.shiftKey;
    this.keyflags.ctrl = event.ctrlKey;
    this.keyflags.alt = event.altKey;
    this.keyflags.meta = event.metaKey;    
}

ToolScene.prototype.handleKeyDown = function( event ) {
    var mode = this.getMode();
    console.log( event );
    console.log( "Key Down ( " + event.key + "/" + event.keyIdentifier + " ), Mode: " + JSON.stringify(mode)  );
    this.updateKeyFlags(event);
   
    switch( mode.primary ){
    case -1: // This is the non-user interaction mode
        break;
    case 0: // This is the default (View, Manipulate, Select)
        if( (event.key == "Del" || event.key == "Delete" || event.keyIdentifier == "U+007F") 
            && (this.selected_hook != null || this.selected_suture != null)){
            if( this.selected_hook != null ){
                this.removeHook( this.selected_hook_index );
                this.selected_hook = null;
                this.selected_hook_index = null;
                this.selected_hook_position = null;
            }
            if( this.selected_suture != null ){
                this.removeSuture( this.selected_suture_index );
                this.selected_suture = null;
                this.selected_suture_index = null;
            }         
        }

        break;
    case 1: // This is the place new hook mode
        break;
    case 2: // This is the place new suture mode
        break;
    case 3: // This is the knife tool mode
        if( (event.key == "Enter" || event.keyIdentifier == "Enter") ) {
            this.EndKnife(true);
            this.clearMode();
        }
        break;
    case 4: // This is the excise tool mode
        break;
    default: // This should never happen.
        console.log( "TOOL MODE ERROR: Unknown tool mode encountered." );
    }

}


ToolScene.prototype.handleKeyUp = function( event ) {
    var mode = this.getMode();
    console.log( "Key Up ( " + event.key +" ), Mode: " + JSON.stringify(mode)  );
    this.updateKeyFlags(event);
}


ToolScene.prototype.handleMouseMove = function( event ) {
    var mode = this.getMode();
    console.log( "Mouse Move, Mode: " + JSON.stringify(mode)  );
    this.updateKeyFlags(event);
    this.mouse.x = event.offsetX;
    this.mouse.y = event.offsetY;

    switch( mode.primary ){
    case -1: // This is the non-user interaction mode
        break;

    case 0: // This is the default (View, Manipulate, Select)
        if( mode.secondary != null && this.keyflags.shift){
            if( mode.secondary == 1 )
                this.incMode();
            this.masterscene.controls.enabled=false;           
            intersection = this.getIntersection( this.mouse, [ this.plane ] );
            this.selected_hook.position.copy( intersection.position.sub(this.offset) );
            this.selected_hook_position.copy( this.selected_hook.position );
        }

        // Update move plane
        intersection = this.getIntersection( this.mouse, this.hooks );
        if( intersection.object != null && intersection.object != this.selected_hook ){
	    this.plane.position.copy( intersection.object.position );
	    this.plane.lookAt( this.masterscene.camera.position );
        }

        break;
    case 1: // This is the place new hook mode
        break;
    case 2: // This is the place new suture mode
        if(mode.secondary != null ){
            this.UpdateDrawSuture(!this.keyflags.ctrl);
        }
        break;
    case 3: // This is the knife tool mode
        break;
    case 4: // This is the excise tool mode
        break;
    default: // This should never happen.
        console.log( "TOOL MODE ERROR: Unknown tool mode encountered." );
    } 




}

ToolScene.prototype.handleMouseUp = function( event ) {
    var mode = this.getMode();
    console.log( "Mouse Up, Mode: " + JSON.stringify(mode)  );
    this.updateKeyFlags(event);
    switch( mode.primary ){
    case -1: // This is the non-user interaction mode
        break;

    case 0: // This is the default (View, Manipulate, Select)
        if( mode.secondary > 1 ){
            this.clearMode();
            this.masterscene.controls.enabled=true;
            this.moveHook( this.selected_hook_index, this.selected_hook_position );
        }
        if( mode.secondary == 1 ){
            this.clearMode();
            this.masterscene.controls.enabled=true;
        }
        break;
    case 1: // This is the place new hook mode
        break;
    case 2: // This is the place new suture mode
        if( mode.secondary != null ){
            this.addSuture();
            this.masterscene.controls.enabled=true;
            this.clearMode();
            this.draw_suture_group.visible = false;
        }
        break;
    case 3: // This is the knife tool mode
        break;
    case 4: // This is the excise tool mode
        break;
    default: // This should never happen.
        console.log( "TOOL MODE ERROR: Unknown tool mode encountered." );
    }

}

ToolScene.prototype.handleMouseDown = function( event ) {
    var mode = this.getMode();
    console.log( "Mouse Down, Mode: " + JSON.stringify(mode)  );
    this.updateKeyFlags(event);
    intersection = this.getIntersection( this.mouse, [ this.plane ] );
    if( intersection.position != null )
        this.offset.copy( intersection.position ).sub( this.plane.position );

    switch( mode.primary ){
    case -1: // This is the non-user interaction mode
        break;

    case 0: // This is the default (View, Manipulate, Select)

        // Try to select hooks and sutures
        this.findManipulatorSelection();
        if( this.selected_hook != null ){ // Prepare for potential drag operation
            this.incMode();            
        }
        else{
            this.clearMode();
        }

        break;
    case 1: // This is the place new hook mode
        var results = this.addHook( this.mouse );
        if( results )
            this.clearMode();
        break;
    case 2: // This is the place new suture mode
        if( event.button == 0 )
            if( this.StartSuture(!this.keyflags.ctrl) ){
                this.incMode();
                this.masterscene.controls.enabled=false;
            }
        if( event.button == 2 ){
            this.clearMode();
            this.masterscene.controls.enabled=true;
            this.draw_suture_group.visible = false;
        }            
        break;
    case 3: // This is the knife tool mode
        if( event.button == 0 && (this.keyflags.shift || this.keyflags.ctrl) ){
            if( mode.secondary == null ){
                var ret = false;
                if( this.keyflags.shift)
                    ret = this.StartKnife(false);
                else if( this.keyflags.ctrl)
                    ret = this.StartKnife(true);
                if( !ret ){
                    this.EndKnife(false);
                    this.clearMode();
                }
                else
                    this.incMode();
            }
            else{
                if( this.keyflags.shift)
                    this.AddKnife(false);
                else if( this.keyflags.ctrl)
                    this.AddKnife(true);
            }
        }
        if( event.button == 2 ){
            this.EndKnife(false);
            this.clearMode();
        }            
        break;
    case 4: // This is the excise tool mode
        if( event.button == 0 ){
            this.exciseRegion();
        }
        this.clearMode();
        break;
    default: // This should never happen.
        console.log( "TOOL MODE ERROR: Unknown tool mode encountered." );
    }

}

ToolScene.prototype.clearManipulatorSelection = function(){
    this.plane.material.visible = false;
    if( this.selected_hook_index != null ){
        this.selected_hook.material = this.hook_material;
        this.selected_hook_index = null;       
        this.selected_hook_position = null;
    }
    this.selected_hook = null;
    
    if( this.selected_suture_index != null ){
        this.selected_suture.children[0].material = this.suture_material;
        this.selected_suture.children[1].material = this.suture_material;
        this.selected_suture.children[2].material = this.suture_material;
        this.selected_suture_index = null;
        this.selected_suture = null;
    }
    this.selected_suture = null;

}


ToolScene.prototype.findManipulatorSelection = function(){
    
    // Do just the hooks for now.
    var manipulators = this.hooks.concat( this.sutures );
    //console.log( manipulators );
    var intersection = this.getIntersection( this.mouse, manipulators );

    //console.log( intersection );

    // Clear existing selection if there was one.
    this.clearManipulatorSelection();
    
    if( intersection.position == null ){
        console.log( "No intersection, clearing selection.");
    }
    else{
        //console.log( "Selected this object:" );
        //console.log( intersection.object );
        var hook_index = this.hooks.indexOf( intersection.object );
        var suture_index = this.sutures.indexOf( intersection.object.parent );
        if( hook_index > -1 ){
            if( this.debug )
                this.plane.material.visible = true;
            if( this.selected_hook != null )
                this.selected_hook.material = this.hook_material;
            this.hooks[hook_index].material = this.hook_selected_material;
            this.selected_hook_index = this.hooks[hook_index].hook_id;
            this.selected_hook = this.hooks[hook_index];
            this.selected_hook_position = new THREE.Vector3().copy( this.selected_hook.position );
        }
        if( suture_index > -1 ){
            this.selected_suture_index = suture_index;
            if( this.selected_suture != null ){
                this.selected_suture.children[0].material = this.suture_material;
                this.selected_suture.children[1].material = this.suture_material;
                this.selected_suture.children[2].material = this.suture_material;
            }
            this.sutures[suture_index].material = this.suture_selected_material;
            this.sutures[suture_index].children[0].material = this.suture_selected_material;
            this.sutures[suture_index].children[1].material = this.suture_selected_material;
            this.sutures[suture_index].children[2].material = this.suture_selected_material;

            this.selected_suture = this.sutures[suture_index];
        }        
    }
    
}

ToolScene.prototype.StartSuture = function(snap_edge) {
    var intersection = this.getIntersection( this.mouse, [this.masterscene.modelscene.dynamic_mesh] );
    
    if( intersection.position == null ){
        this.suture_placement.start = null;
        this.suture_placement.end = null;
        return false;
    }

    if( snap_edge ) {
        var closest_sharppoint = this.masterscene.modelscene.ClosestPointOnSharpEdge( intersection.position,.1);
        if(closest_sharppoint===null){
            this.suture_placement.start = null;
            this.suture_placement.end = null;
            return false;
        }
        intersection.bary_coord.copy(closest_sharppoint.bary_coords);
        intersection.position.copy(closest_sharppoint.edge_point);
        intersection.triangle.copy(closest_sharppoint.tri);
        intersection.triangleID = closest_sharppoint.triangleID;
        intersection.indices = closest_sharppoint.indices.slice();
    }
    
    this.draw_suture_group.visible = true;
    
    this.draw_suture_fixed.scale.set( this.scaling_factor*0.01, this.scaling_factor*0.01, this.scaling_factor*0.01 );
    this.draw_suture_free.scale.set( this.scaling_factor*0.01, this.scaling_factor*0.01, this.scaling_factor*0.01 );
    this.draw_suture_group.position.copy( intersection.position );
    this.draw_suture_free.position.set( 0.0, 0.0, 0.0 );
    this.draw_suture_segment.visible = false;
    this.suture_placement.start = intersection
    this.suture_placement.end = intersection;

    this.draw_suture_fixed.material = this.suture_selected_material;
    this.draw_suture_free.material = this.suture_selected_material;
    this.draw_suture_segment.material = this.suture_selected_material;
    this.plane.position.copy( intersection.position );
    this.plane.lookAt( this.masterscene.camera.position );


    return true;
}

ToolScene.prototype.UpdateDrawSuture = function(snap_edge) {

    var Ex = ( this.mouse.x / this.masterscene.renderer.domElement.clientWidth ) * 2 - 1;
    var Ey = - ( this.mouse.y / this.masterscene.renderer.domElement.clientHeight ) * 2 + 1;
    var vector = new THREE.Vector3(Ex, Ey, 0.5).unproject( this.masterscene.camera );

    var intersection = this.getIntersection( this.mouse, [this.masterscene.modelscene.dynamic_mesh] );

    if( intersection.position == null ){
        intersection = this.getIntersection( this.mouse, [ this.plane ] );
        this.draw_suture_free.position.copy( intersection.position );
        this.draw_suture_free.position.sub( this.draw_suture_group.position );
        this.suture_placement.end = null;
        this.draw_suture_free.material = this.suture_bad_material;
        this.draw_suture_segment.material = this.suture_bad_material;
    }
    else{
        if( snap_edge ) {
            var closest_sharppoint = this.masterscene.modelscene.ClosestPointOnSharpEdge( intersection.position,.1);
            if(closest_sharppoint===null){
                intersection = this.getIntersection( this.mouse, [ this.plane ] );
                this.draw_suture_free.position.copy( intersection.position );
                this.draw_suture_free.position.sub( this.draw_suture_group.position );
                this.suture_placement.end = null;
                this.draw_suture_free.material = this.suture_bad_material;
                this.draw_suture_segment.material = this.suture_bad_material;
            }
            else{
                intersection.bary_coord.copy(closest_sharppoint.bary_coords);
                intersection.position.copy(closest_sharppoint.edge_point);
                intersection.triangle.copy(closest_sharppoint.tri);
                intersection.triangleID = closest_sharppoint.triangleID;
                intersection.indices = closest_sharppoint.indices.slice();
                this.draw_suture_free.position.copy( intersection.position );
                this.draw_suture_free.position.sub( this.draw_suture_group.position );
                this.suture_placement.end = intersection;
                this.draw_suture_free.material = this.suture_selected_material;
                this.draw_suture_segment.material = this.suture_selected_material;
            }
        }
        else{
            this.draw_suture_free.position.copy( intersection.position );
            this.draw_suture_free.position.sub( this.draw_suture_group.position );
            this.suture_placement.end = intersection;
            this.draw_suture_free.material = this.suture_selected_material;
            this.draw_suture_segment.material = this.suture_selected_material;
        }
    }

    this.draw_suture_segment.visible = true;
    var center = new THREE.Vector3();
    center.copy( this.draw_suture_free.position );
    center.sub( this.draw_suture_fixed.position );
    var length = center.length();
    center.divideScalar( 2.0 );
    this.draw_suture_segment.position.copy( center );
    this.draw_suture_segment.scale.set( this.scaling_factor*0.008, length, this.scaling_factor*0.008 );
    if( length > 0.000001 ){
        center.normalize();
        this.draw_suture_segment.quaternion.setFromUnitVectors( new THREE.Vector3(0.0, 1.0, 0.0), center);
    }

    console.log( this.suture_placement.start );
    console.log( this.suture_placement.end );

}

ToolScene.prototype.StartKnife = function(snap_edge) {    
    this.draw_incision = new THREE.Object3D();
    this.masterscene.scene.add( this.draw_incision );
    this.incision_placement.length = 0;
    var ret = this.AddKnife(snap_edge);
    return ret;
}

ToolScene.prototype.AddKnife = function(snap_edge) {
    var Ex = ( this.mouse.x / this.masterscene.renderer.domElement.clientWidth ) * 2 - 1;
    var Ey = - ( this.mouse.y / this.masterscene.renderer.domElement.clientHeight ) * 2 + 1;
    var vector = new THREE.Vector3(Ex, Ey, 0.5).unproject( this.masterscene.camera );
    var intersection = this.getIntersection( this.mouse, [this.masterscene.modelscene.dynamic_mesh] );

    if( intersection.position == null )
        return false;

    console.log( intersection.position );
    intersection.edge_snap = false;
    if( snap_edge ) { // snap intersection to the closest "hanging" edge - these are "sharp"
                      // edges, or snap to the clos

        var closest_sharppoint;
        if( this.incision_placement.length > 0 )
            closest_sharppoint = this.masterscene.modelscene.ClosestPointOnSharpEdge( intersection.position,.1,
                                                                                      this.incision_placement[this.incision_placement.length-1].position);
        else
            closest_sharppoint = this.masterscene.modelscene.ClosestPointOnSharpEdge( intersection.position,.1);
        

        var closest_knifepoint = {
            edge_point: new THREE.Vector3(),
            distance: 1e8,
            good: false
        }

        for( var k = 0; k < this.incision_placement.length-3; k++ ){
            var seg_start = this.incision_placement[k];
            var seg_end = this.incision_placement[k+1];
            
            var offset = new THREE.Vector3().copy( intersection.position );            
            var knife_ray = new THREE.Ray( this.incision_placement[this.incision_placement.length-1].position,
                                           offset.sub( this.incision_placement[this.incision_placement.length-1].position).normalize() );
            var segment_point = new THREE.Vector3();
            knife_ray.distanceSqToSegment(seg_start.position, seg_end.position,
                                          null, segment_point );
            var distance = intersection.position.distanceTo( segment_point );

            if( distance < closest_knifepoint.distance ){              
                closest_knifepoint.distance = distance;
                closest_knifepoint.edge_point.copy( segment_point );        
                closest_knifepoint.good = true;
            }
        }

        final_point = {}

        if(!(closest_sharppoint===null)){
            intersection.bary_coord.copy(closest_sharppoint.bary_coords);
            intersection.position.copy(closest_sharppoint.edge_point);
            intersection.triangle.copy(closest_sharppoint.tri);
            intersection.triangleID = closest_sharppoint.triangleID;
            intersection.indices = closest_sharppoint.indices.slice();
            intersection.edge_snap = true;
        }

        if( (closest_knifepoint.good &&
             (closest_sharppoint===null || closest_knifepoint.distance < closest_sharppoint.distance)) ){
            var vector = new THREE.Vector3().copy( closest_knifepoint.edge_point );
            var cameraPos = this.masterscene.camera.position;
            var rayCaster = new THREE.Raycaster(cameraPos,
                                                vector.sub(cameraPos).normalize());
            var intersections = rayCaster.intersectObjects([this.masterscene.modelscene.dynamic_mesh], true);
            if( intersections.length > 0 ){
                knife_intersection = intersections[0];
                if( knife_intersection.object.geometry instanceof THREE.BufferGeometry ){
                    positionsArray = knife_intersection.object.geometry.attributes.position.array;
                    vA = new THREE.Vector3().fromArray(positionsArray, knife_intersection.face.a*3)
                    vB = new THREE.Vector3().fromArray(positionsArray, knife_intersection.face.b*3)
                    vC = new THREE.Vector3().fromArray(positionsArray, knife_intersection.face.c*3)
                    var tri = new THREE.Triangle( vA, vB, vC );
                    triId = knife_intersection.faceIndex;
                    bary = tri.barycoordFromPoint( knife_intersection.point );
                    intersection.position.copy( knife_intersection.point );
                    intersection.triangleID = triId;
                    intersection.bary_coord.copy(bary);
                    intersection.triangle.copy(tri);
                    intersection.indices = [ knife_intersection.face.a, knife_intersection.face.b, knife_intersection.face.c ];
                    intersection.edge_snap= true;
                    
                }
            }           
        }
    }

    this.incision_placement.push( intersection );
    this.UpdateKnifeObject();
    
    return true;
}

ToolScene.prototype.UpdateKnifeObject = function() {
    this.masterscene.scene.remove( this.draw_incision );
    this.draw_incision = new THREE.Object3D();
    var line_material = new THREE.LineBasicMaterial({
	color: 0x0000ff,
        linewidth: this.scaling_factor*25
    });
    var point_material = new THREE.MeshBasicMaterial({
	color: 0x9900ff
    });

    var knife_point = new THREE.IcosahedronGeometry(1.0, 1);
    var geometry = new THREE.Geometry();
    for( var i = 0; i < this.incision_placement.length; i ++){
        geometry.vertices.push( this.incision_placement[i].position );
        var point = new THREE.Mesh( knife_point, point_material );
        point.scale.set( this.scaling_factor*0.01, this.scaling_factor*0.01, this.scaling_factor*0.01 );
        point.position.copy( this.incision_placement[i].position );
        this.draw_incision.add( point );
    }
    var line = new THREE.Line( geometry, line_material );
    this.draw_incision.add( line );
    this.masterscene.scene.add( this.draw_incision );
}

ToolScene.prototype.EndKnife = function(complete) {
    this.masterscene.scene.remove( this.draw_incision );
    this.draw_incision = null;
    if( this.incision_placement.length > 0 && complete ){
        knife_data = []

        var getPosition = function( dynGeometry, tri, uv ){
            var indices = dynGeometry.getIndex().array;
            var positions = dynGeometry.getAttribute( "position" ).array;
            vA = new THREE.Vector3().fromArray(positions, indices[tri*3+0]*3)
            vB = new THREE.Vector3().fromArray(positions, indices[tri*3+1]*3)
            vC = new THREE.Vector3().fromArray(positions, indices[tri*3+2]*3)
            var result = new THREE.Vector3();
            result.copy(vA.multiplyScalar(uv.x));
            result.add( vB.multiplyScalar(uv.y));
            result.add( vC.multiplyScalar(uv.z));
            return result;
        }

        var normals = this.masterscene.modelscene.dynamic_geometry.getAttribute( "normal" ).array;
        var indices = this.masterscene.modelscene.dynamic_geometry.getIndex().array;
        var normal = new THREE.Vector3();
        var position = new THREE.Vector3();
        var i,vNorm = new THREE.Vector3();
        for( var k = 0; k < this.incision_placement.length; k++){
            normal.set( 0,0,0);
            for (i = 0; i < 3; i ++ ) {
                vNorm.set(normals[indices[this.incision_placement[k].triangleID*3+i]*3],
                          normals[indices[this.incision_placement[k].triangleID*3+i]*3+1],
                          normals[indices[this.incision_placement[k].triangleID*3+i]*3+2]);
                normal.add(vNorm); }
            normal.normalize();
            position.copy( getPosition( this.masterscene.modelscene.dynamic_geometry,
                                        this.incision_placement[k].triangleID,
                                        this.incision_placement[k].bary_coord ));

            knife_data.push( { "triangle":  this.incision_placement[k].triangleID, 
                               "coords": [ this.incision_placement[k].bary_coord.x, 
                                           this.incision_placement[k].bary_coord.y ],
                               "position" : [position.x, position.y, position.z],
                               "normal" : [normal.x, normal.y, normal.z] } );                        
        }

        var msg = {
            command: "makeIncision",
            data: {
                path: knife_data,
                edge_start: this.incision_placement[0].edge_snap,
                edge_end: this.incision_placement[this.incision_placement.length-1].edge_snap
            }
        }; 

        this.masterscene.comm.send( JSON.stringify( msg ) );
    }
}


ToolScene.prototype.addHook = function(location) {
    var intersection = this.getIntersection( this.mouse, [this.masterscene.modelscene.dynamic_mesh] );
    
    if( intersection.position == null )
        return false;
    else{
        var msg =  {
            command: "addHook",
            data: {
                triangle: intersection.triangleID,
                coords: [intersection.bary_coord.x,  intersection.bary_coord.y,  intersection.bary_coord.z] 
            }
        }
        
        this.masterscene.comm.send( JSON.stringify( msg ) );
        return true;
    }   
}


ToolScene.prototype.removeHook = function(id) {
    console.log( "Deleting hook " + id );

    msg = {
        command: "deleteHook", 
        data: {
            hook_id: id
        }
    };

    this.masterscene.comm.send( JSON.stringify( msg ) );
    return true;
}

ToolScene.prototype.moveHook = function(id, newlocation) {
    console.log( "Moving hook " + id + " to location " + newlocation.x + " " + newlocation.y + " " + newlocation.z );

    msg = {
        command: "moveHook", 
        data: {
            hook_id: id,
            coords: [newlocation.x, newlocation.y, newlocation.z]
        }
    };

    this.masterscene.comm.send( JSON.stringify( msg ) );
    return true;
}



ToolScene.prototype.addSuture = function() {

    if( this.suture_placement.start == null ||
        this.suture_placement.end == null )
        return false;

    msg = {
        command: "addSuture",
        data : {
            triangleA: this.suture_placement.start.triangleID,
            triangleB: this.suture_placement.end.triangleID,
            uvA: [this.suture_placement.start.bary_coord.x,this.suture_placement.start.bary_coord.y,this.suture_placement.start.bary_coord.z],
            uvB: [this.suture_placement.end.bary_coord.x,this.suture_placement.end.bary_coord.y,this.suture_placement.end.bary_coord.z]
        }
    }

    console.log( msg );

    this.masterscene.comm.send( JSON.stringify( msg ) );
    return true;
}

ToolScene.prototype.removeSuture = function(id) {
    console.log( "Deleting suture " + id );

    msg = {
        command: "deleteSuture", 
        data: {
            suture_id: id
        }
    };

    this.masterscene.comm.send( JSON.stringify( msg ) );
    return true;
}

ToolScene.prototype.exciseRegion = function() {
    console.log( "Excising tissue region." );
    var intersection = this.getIntersection( this.mouse, [this.masterscene.modelscene.dynamic_mesh] );

    msg = {
        command : "exciseRegion",
        data : {
            triangle: intersection.triangleID
        }
    }
    
    this.masterscene.comm.send( JSON.stringify( msg ) );
    return true;

}


ToolScene.prototype.updateHooks = function() {
    this.hooks.forEach( function( item, index, array ){
        this.masterscene.scene.remove( item );
    }.bind(this));
    this.hooks.length=0;

     this.hook_data.forEach( function( item, index, array ){
        var hook = this.hook_mesh.clone();
        hook.position.set( item.pos[0], item.pos[1], item.pos[2] );
        hook.scale.set( this.scaling_factor*0.12,
                        this.scaling_factor*0.12,
                        this.scaling_factor*0.12 );
        var normal = new THREE.Vector3();
        var normals = this.masterscene.modelscene.dynamic_geometry.getAttribute( "normal" ).array;
        var indices = this.masterscene.modelscene.dynamic_geometry.getIndex().array;
        var i,vNorm = new THREE.Vector3();
        for (i = 0; i < 3; i ++ ) {
            vNorm.set(normals[indices[item.triangle*3+i]*3],normals[indices[item.triangle*3+i]*3+1],normals[indices[item.triangle*3+i]*3+2]);
            normal.add(vNorm); }
        normal.normalize();
        var rotation_vec = new THREE.Vector3(0,0,1);
        var angle = rotation_vec.angleTo(normal);
        rotation_vec.cross(normal);
        rotation_vec.normalize();
         hook.rotateOnAxis(rotation_vec, angle );
         hook.hook_id = item.id;
        this.hooks.push( hook );
    }.bind(this));
    this.hooks.forEach( function( item, index, array ){
        if( item.hook_id == this.selected_hook_index ){
            this.selected_hook = item;
            this.selected_hook.position.copy( this.selected_hook_position );
            this.selected_hook.material = this.hook_selected_material;
        }
        this.masterscene.scene.add( item );
    }.bind(this));
}

ToolScene.prototype.updateSutures = function() {
    this.sutures.forEach( function( item ){  // , index, array
        this.masterscene.scene.remove( item );
    }.bind(this));
    this.sutures.length=0; 

    var getSutureMesh = function( suture_object, pntA, pntB, scaling_factor ){
        var sut = suture_object.clone();
        sut.position.copy( pntA );
        sut.children[0].position.set(0.0,0.0,0.0);
        sut.children[1].position.subVectors( pntB, pntA );
        var center = new THREE.Vector3();
        center.copy( sut.children[1].position );
        center.sub(  sut.children[0].position );
        var length = center.length();
        center.divideScalar(2.0);
        sut.children[2].position.copy(center);
        sut.children[0].scale.set( scaling_factor*0.01,  scaling_factor*0.01, scaling_factor*0.01);
        sut.children[1].scale.set( scaling_factor*0.01,  scaling_factor*0.01, scaling_factor*0.01);
        sut.children[2].scale.set( scaling_factor*0.008, length, scaling_factor*0.008);
        if(length > 0.0000001) {
            center.normalize();
            sut.children[2].quaternion.setFromUnitVectors(new THREE.Vector3(0.0, 1.0, 0.0),center); }
        return sut;
    }

    this.suture_data.forEach( function( item, index, array){
        var getPosition = function( dynGeometry, tri, uv ){
            var indices = dynGeometry.getIndex().array;
            var positions = dynGeometry.getAttribute( "position" ).array;
            var uv3 = new THREE.Vector3( 1.0 - (uv[0]+uv[1]), uv[0], uv[1] );
            vA = new THREE.Vector3().fromArray(positions, indices[tri*3+0]*3)
            vB = new THREE.Vector3().fromArray(positions, indices[tri*3+1]*3)
            vC = new THREE.Vector3().fromArray(positions, indices[tri*3+2]*3)
            var result = new THREE.Vector3();
            result.copy(vA.multiplyScalar(uv3.x));
            result.add( vB.multiplyScalar(uv3.y));
            result.add( vC.multiplyScalar(uv3.z));
            return result;
        }
        var ptA = getPosition(this.masterscene.modelscene.dynamic_geometry,
                              item.A.triangle,item.A.uv);
        var ptB = getPosition(this.masterscene.modelscene.dynamic_geometry,
                              item.B.triangle,item.B.uv);
        var sutMesh = getSutureMesh(this.suture_object, ptA, ptB, this.scaling_factor );
     
        this.sutures.push( sutMesh );
    }.bind(this));

    this.sutures.forEach( function( item, index ){  // , index, array
        if( index == this.selected_suture_index ){
            this.selected_suture = item;
            this.selected_suture.children[0].material = this.suture_selected_material;
            this.selected_suture.children[1].material = this.suture_selected_material;
            this.selected_suture.children[2].material = this.suture_selected_material;
        }
        this.masterscene.scene.add( item );
    }.bind(this));
}


ToolScene.prototype.toggleDebugMode = function( debugFlag ) {
    this.debug = debugFlag;

    if( !this.debug ){
        if( this.plane )
            this.plane.material.visible = false;
    }

}

THREE.OBJExporter = function () {};

THREE.OBJExporter.prototype = {

    constructor: THREE.OBJExporter,

    parse: function ( object ) {

	var output = '';

	var indexVertex = 0;
	var indexVertexUvs = 0
	var indexNormals = 0;

	var parseObject = function ( child ) {

	    var nbVertex = 0;
	    var nbVertexUvs = 0;
	    var nbNormals = 0;

	    var geometry = child.geometry;

	    if ( geometry instanceof THREE.Geometry ) {

		output += 'o ' + child.name + '\n';

		for ( var i = 0, l = geometry.vertices.length; i < l; i ++ ) {

		    var vertex = geometry.vertices[ i ].clone();
		    vertex.applyMatrix4( child.matrixWorld );

		    output += 'v ' + vertex.x + ' ' + vertex.y + ' ' + vertex.z + '\n';

		    nbVertex ++;

		}

		// uvs
                /*
		  for ( var i = 0, l = geometry.faceVertexUvs[ 0 ].length; i < l; i ++ ) {

		  var vertexUvs = geometry.faceVertexUvs[ 0 ][ i ];

		  for ( var j = 0; j < vertexUvs.length; j ++ ) {

		  var uv = vertexUvs[ j ];
		  vertex.applyMatrix4( child.matrixWorld );

		  output += 'vt ' + uv.x + ' ' + uv.y + '\n';

		  nbVertexUvs ++;

		  }

		  }

		  // normals

		  for ( var i = 0, l = geometry.faces.length; i < l; i ++ ) {

		  var normals = geometry.faces[ i ].vertexNormals;

		  for ( var j = 0; j < normals.length; j ++ ) {

		  var normal = normals[ j ];
		  output += 'vn ' + normal.x + ' ' + normal.y + ' ' + normal.z + '\n';

		  nbNormals ++;

		  }

		  }
                */

		// faces

		for ( var i = 0, j = 1, l = geometry.faces.length; i < l; i ++, j += 3 ) {

		    var face = geometry.faces[ i ];
                    
		    output += 'f ';
		    output += ( indexVertex + face.a + 1 ) + ' '; // + '/' + ( /*indexVertexUvs*/ + j ) + '/' + ( /*indexNormals*/ + j ) + ' ';
		    output += ( indexVertex + face.b + 1 ) + ' '; // + '/' + ( /*indexVertexUvs*/ + j + 1 ) + '/' + ( /*indexNormals*/ + j + 1 ) + ' ';
		    output += ( indexVertex + face.c + 1 ) + ' '; // + '/' + ( /*indexVertexUvs*/ + j + 2 ) + '/' + ( /*indexNormals*/ + j + 2 ) + '\n';
                    output += '\n';
                    
		}

	    }

	    // update index
	    indexVertex += nbVertex;
	    indexVertexUvs += nbVertexUvs;
	    indexNormals += nbNormals;

	};

        
        if( Array.isArray(object) ){
            for( var i = 0; i < object.length; i++ )
                object[i].traverse( parseObject );
        }
        else
	    object.traverse( parseObject );

	return output;

    }

};





ToolScene.prototype.DumpToolsAsObj = function (){
    console.log( "Dumping Tools as OBJ" );

    $("#model").empty()
    var BlobData = []

    var exporter = new THREE.OBJExporter();
    BlobData.push(exporter.parse(this.sutures));
    var oBlob = new Blob(BlobData, {type : 'text/plain;charset=utf-8'}); // the blob
    saveAs( oBlob, "tools_dump.obj" );
        
}
