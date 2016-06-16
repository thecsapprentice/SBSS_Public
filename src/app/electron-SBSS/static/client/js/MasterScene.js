// Class representing global scene information
//
//
//
//

var MasterScene = function() {
    this.container = null;
    this.comm = null;
    this.modelscene = null;
    this.toolscene = null;
    this.toolsgui = null;
    this.enable_debug = false;
}

MasterScene.prototype.initialize = function() {
    // SCENE
    this.scene = new THREE.Scene();

    // CAMERA
    this.camera = new THREE.PerspectiveCamera(40, window.innerWidth/window.innerHeight, 0.001, 100);
    this.camera.position.z = 0.5;

    // FPS Meter
    //this.stats = new Stats();
    //this.stats.domElement.style.bottom = '0px';
    //this.stats.domElement.style.right = '0px';
    //this.stats.domElement.style.zIndex = 100;
    //this.container.appendChild( this.stats.domElement );

    // RENDERER
    if ( ! Detector.webgl ) Detector.addGetWebGLMessage();
    this.renderer = new THREE.WebGLRenderer( {antialias:true} );
    this.renderer.setSize(this.container.clientWidth, this.container.clientHeight);
    this.container.appendChild( this.renderer.domElement );

    // GUI
    this.toolsgui = new GUI(this);
    this.toolsgui.initialize();

    // TOOL SCENE
    this.toolscene = new ToolScene(this);
    this.toolscene.initialize();

    // CONTROLS
    this.controls = new THREE.TrackballControls( this.camera, this.renderer.domElement );
    this.controls.rotateSpeed = 5.0;
    this.controls.zoomSpeed = 1.2;
    this.controls.panSpeed = 0.8;
    this.controls.noZoom = false;
    this.controls.noPan = false;
    this.controls.staticMoving = true;
    this.controls.dynamicDampingFactor = 0.3;

    // LIGHTS
    this.ambientLight = new THREE.AmbientLight( 0x474747 );
    this.scene.add( this.ambientLight );
    this.cameraLight = new THREE.DirectionalLight( 0xffffff );
    this.cameraLight.intensity = 0.8;
    this.scene.add( this.cameraLight );
    
    // MODEL SCENE
    this.modelscene = new ModelScene(this);
    this.modelscene.initialize();
   

    // SIGNALS
    this.debug_mode_toggled = new signals.Signal();

    // HOOKS
    window.addEventListener( 'resize', this.onWindowResize.bind(this), false );
    window.addEventListener( 'unload', this.onShutdown.bind(this) );
    this.comm.opened.add( this.CommOpened.bind(this) );
    this.comm.messageReceived.add(this.CommMessage.bind(this) );
    this.comm.closed.add( this.CommClosed.bind(this) );
    this.comm.givingup.add( this.CommFailed.bind(this) );
    this.modelscene.dynamicVerticesUpdated.add(this.toolscene.updateSutures.bind(this.toolscene))
    this.modelscene.dynamicVerticesUpdated.add(this.toolscene.updateHooks.bind(this.toolscene))
    this.debug_mode_toggled.add( this.modelscene.toggleDebugMode.bind(this.modelscene) );
    this.debug_mode_toggled.add( this.toolscene.toggleDebugMode.bind(this.toolscene) );

    
}

MasterScene.prototype.onWindowResize = function(event) {
    this.camera.aspect = this.container.clientWidth / this.container.clientHeight;
    this.camera.updateProjectionMatrix();
    this.renderer.setSize(this.container.clientWidth, this.container.clientHeight);
    this.controls.handleResize();
    this.render(0);
}

MasterScene.prototype.onShutdown = function(event) {
    console.log( "Shutting down connection to server.");
    this.comm.teardown();
}

MasterScene.prototype.CommOpened = function() {
    console.log("MasterScene: Communications Established.");
    $("#connection_status").html("SIMULATOR ONLINE");
    this.modelscene.clearScene();
    this.toolscene.clearScene();
    if(this.controls)
        this.controls.reset();
}

MasterScene.prototype.CommMessage = function(data) {
    //console.log("MasterScene: Message Received.");
    //console.log( data );
    var toolData = {}
    var modelData = {}

    toolData.hooks = data.hooks;
    toolData.sutures = data.sutures;
    toolData.simulator_state = data.simulator_state;
    toolData.timestamp = data.timestamp;
    
    modelData.dynamic = data.dynamic;
    modelData.static = data.static;
    modelData.textures = data.textures;
    modelData.timestamp = data.timestamp;

    this.processModels( modelData ); // Process model data first, as tools depend on it.
    this.processTools( toolData );
}

MasterScene.prototype.CommClosed = function() {
    console.log("MasterScene: Communications Terminated.");
    $("#connection_status").html("SIMULATOR OFFLINE");
}

MasterScene.prototype.CommFailed = function() {
    window.location.assign("/");
}


MasterScene.prototype.render = function(timestamp) {
    if(this.start_render_timestamp === null) this.start_render_timestamp = timestamp;
    if(this.last_render_timestamp === null) this.last_render_timestamp = timestamp;

    requestAnimationFrame(this.render.bind(this));
    this.update(timestamp);
    this.renderer.render(this.scene, this.camera);

    this.last_render_timestamp = timestamp;
}

MasterScene.prototype.update = function(timestamp){
    this.controls.update();
    //this.stats.update();
    //if( !connected ){
    //    $("#connection_status").html("SIMULATOR OFFLINE");
    //    connectionNextTime -= timestamp - last_render_timestamp;
    //}
}

MasterScene.prototype.processTools = function(data){
    this.toolscene.processData( data );

}

MasterScene.prototype.processModels = function(data){
    this.modelscene.processData( data );

}
