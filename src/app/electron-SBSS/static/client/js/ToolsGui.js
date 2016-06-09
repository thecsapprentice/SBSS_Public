// Class to handle GUI operations
//
//
//
//
//


var GUI = function(masterscene){
    this._toolState = 0;
    this._exampleFile = 0;
    this._sceneFile = 0;
    this._developersOn = true;
    this.buttons = []
    this.locked = false;
    this.masterscene = masterscene;
}

GUI.prototype.initialize = function() {
    this.button = [1,1,1,1,1,1,1,1,1,1];
    this.buttons[0] = document.getElementById("view");
    this.buttons[0].addEventListener("click", this.clearToolstate.bind(this), false);
    this.buttons[1] = document.getElementById("hook");
    this.buttons[1].addEventListener("click", this.setToolstate.bind(this,1), false);
    this.buttons[2] = document.getElementById("suture");
    this.buttons[2].addEventListener("click", this.setToolstate.bind(this,2), false);
    this.buttons[3] = document.getElementById("knife");
    this.buttons[3].addEventListener("click", this.setToolstate.bind(this,3), false);
    this.buttons[4] = document.getElementById("excise");
    this.buttons[4].addEventListener("click", this.setToolstate.bind(this,4), false);
    this.buttons[5] = document.getElementById("load_scene");
    //this.buttons[5].addEventListener("click", this.loadScene.bind(this), false);
    this.buttons[6] = document.getElementById("load_sequence");
    //this.buttons[6].addEventListener("click", this.loadExampleFile.bind(this), false);
    this.buttons[7] = document.getElementById("next");
    this.buttons[7].addEventListener("click", this.nextButton.bind(this), false);
    this.buttons[8] = document.getElementById("save_sequence");
    this.buttons[8].addEventListener("click", this.saveExampleFile.bind(this), false);
    this.buttons[9] = document.getElementById("help");
    this.buttons[9].addEventListener("click", this.helpFn.bind(this), false);
    this.buttons[10] = document.getElementById("about");
    //this.buttons[10].addEventListener("click", this.aboutFn.bind(this), false);
    this.buttons[11] = document.getElementById("debug");
    this.buttons[11].addEventListener("click", this.debugFn.bind(this), false);
    this.buttons[12] = document.getElementById("disconnect");
    this.buttons[12].addEventListener("click", function() { window.location.assign("/"); }, false );

    $('#exampleBox').on('show.bs.modal', this.initializeExampleBox.bind(this) )
    $('#sceneBox').on('show.bs.modal', this.initializeSceneBox.bind(this) )
    $('#displaysettings').click( function(event) {
        $('#settings-menu').toggleClass('active');
    }.bind(this));

    $("#data-overlay-percent").slider({
	tooltip: 'always'
    });
    $("#data-overlay-percent").on("slide", function(slideEvt) {
        this.masterscene.modelscene.setDataOverlayPercent(slideEvt.value);
    }.bind(this))

    $("#data-overlay-strain_range").slider({
        id: "slider-data-overlay-strain_range",
        min: 0,
        max: 75000,
        range: true,
        value: [0, 100],
	step: 10
    });
    $("#data-overlay-strain_range").on("slide", function(slideEvt) {
        this.masterscene.modelscene.setDataOverlayStrainRange(slideEvt.value[0], slideEvt.value[1]);
    }.bind(this))

    $("#data-overlay-stress_range").slider({
        id: "slider-data-overlay-stress_range",
        min: 1,
        max: 100,
        range: true,
        value: [1, 1],
	step: 1
    });
    $("#data-overlay-stress_range").on("slide", function(slideEvt) {
        this.masterscene.modelscene.setDataOverlayStressRange(slideEvt.value[0], slideEvt.value[1]);
    }.bind(this))

    
    $("#data-overlay-use-magnitude").on("click", function(event) {
        this.masterscene.modelscene.setDataOverlayUseMagnitude( $("#data-overlay-use-magnitude").is(':checked'));
    }.bind(this))

    $("#data-overlay-use-strain").on("click", function(event) {
        this.masterscene.modelscene.setDataOverlaySource( $("input[name=data-overlay-selector]:checked").val() )
    }.bind(this))

    $("#data-overlay-use-stress").on("click", function(event) {
        this.masterscene.modelscene.setDataOverlaySource( $("input[name=data-overlay-selector]:checked").val() )
    }.bind(this))

    
   document.getElementById("dump").addEventListener( "click" , function(){
        this.masterscene.toolscene.DumpToolsAsObj()
        this.masterscene.modelscene.DumpModelAsObj()
    }.bind(this)
                                                      , false );
}

GUI.prototype.setToolstate = function(value) {
    if( this.locked )
        return;

    if( value == -1 ){
        this._toolState = value;
        this.setProcessing(true);
        for(i=0; i<5; i=i+1) {
            $(this.buttons[i]).addClass("disabled");
        }
        this.locked = true;
    }
    else {
        this.setProcessing(false);
        this._toolState = value;
        var i;
        for(i=0; i<5; i=i+1) {
            if(i===value)
                $(this.buttons[i]).removeClass("disabled");
            else
                $(this.buttons[i]).addClass("disabled");
        }
    }
}

GUI.prototype.clearToolstate = function() {
    this.locked = false;
    this.setToolstate( 0 );
}

GUI.prototype.setProcessing = function(toggle) {
    if( toggle ){
        $("#processingMessage").removeClass( "disabled" );   
        var of = $("#container").offset();
        var size = { width: $("#container").innerWidth(), height: $("#container").innerHeight() };
        $("#processingMessage").width( size.width/2 );
        $("#processingMessage").offset( { top: of.top + size.height/2 - size.height/4,
                                   left: of.left + size.width/2 - size.width/4 })
    }
    else
        $("#processingMessage").addClass( "disabled" );
}

GUI.prototype.initializeSceneBox = function(event) {
    var modal = $('#sceneBox')
    
    $("#scenes").html("");
    var requester = new XMLHttpRequest();
    requester.onreadystatechange = function(requester) {
        if( requester.readyState == 4 && requester.status == 200 ){
                console.log( requester.responseText )
            response_data = JSON.parse( requester.responseText );
            response_data.scenes.sort();
            for( var i = 0; i< response_data.scenes.length; i++){
                $("#scenes").append( '<button id="scene_'+(i+1)+'" class="btn btn-block scene_choice" scene_name="'+response_data.scenes[i]+'"scene_number="'+(i+1)+'">'+response_data.scenes[i]+'</button>');
                }
            
            $(".scene_choice").on("click", function( ev ) {
                var scenename = $(ev.target).attr("scene_name")
                var msg = {
                    command: "loadScene",
                    data : {
                            name: scenename
                    }
                };
                
                var jsonStr = JSON.stringify(msg);
                this.masterscene.comm.send(jsonStr);
                    $(this.buttons[5]).click();
            }.bind(this));
        }
    }.bind( this, requester )
    
    requester.open( "GET", "sbss-data://data/scenes", true );
    requester.send();

}

GUI.prototype.initializeExampleBox = function(event) {
    var modal = $('#exampleBox')

    $("#examples").html("");
    var requester = new XMLHttpRequest();
    requester.onreadystatechange = function(requester) {
        if( requester.readyState == 4 && requester.status == 200 ){
            console.log( requester.responseText )
            response_data = JSON.parse( requester.responseText );
            response_data.histories.sort();
            for( var i = 0; i< response_data.histories.length; i++){
                $("#examples").append( '<button id="example_'+(i+1)+'" class="btn btn-block example_choice" history_name="'+response_data.histories[i]+'"history_number="'+(i+1)+'">'+response_data.histories[i]+'</button>');
            }
            
            $(".example_choice").on("click", function( ev ) {
                var examplename = $(ev.target).attr("history_name")
                var msg = {
                    command: "loadHistory",
                    data : {
                        name: examplename
                    }
                };
                
                var jsonStr = JSON.stringify(msg);
                this.masterscene.comm.send(jsonStr);
                $(this.buttons[7]).removeClass("disabled");
                $(this.buttons[6]).click();
            }.bind(this));
        }
    }.bind( this, requester )
    
    requester.open( "GET", "sbss-data://data/histories", true );
    requester.send();        

}

GUI.prototype.nextButton = function() {
    if( $(this.buttons[7]).hasClass("disabled") )
        return
    
    var msg = {
        command: "historyNext",
        data : {}
    };
    var jsonStr = JSON.stringify(msg);
    this.masterscene.comm.send(jsonStr);
}

GUI.prototype.saveExampleFile = function() {
    if( $(this.buttons[8]).hasClass("disabled") )
        return
    
    var seq_Name=prompt("Please enter a name to save your sequence as: ","0");

    var msg = {
        command: "saveHistory",
        data : {
            name: seq_Name
        }
    };
    
    var jsonStr = JSON.stringify(msg);
    this.masterscene.comm.send(jsonStr);  
}


GUI.prototype.doIt = function () {
    alert('Testing, testing-');
}


GUI.prototype.helpFn = function () {
    var hNum=prompt("Mouse drag while pressing:\nleftMouse rotates-\nmiddleMouse zooms-\nrightMouse pans-\n\nWhich item do you need help on:\n\n1. TOOLS\n2. SCENE\n3. EXAMPLES\n","0");
    if(hNum==="1")
        alert("TOOLS not yet implemented. Only radio buttons work.");
    else if(hNum==="2")
        alert("Click LOAD to select a surgical scene file number.");
    else if(hNum==="3")
        alert("First click LOAD to select an example file number.  After scene loads each click on NEXT performs the next action for that example.");
    else ;
}

GUI.prototype.debugFn = function () {
    this.masterscene.enable_debug = !this.masterscene.enable_debug;   
    if( this.masterscene.enable_debug )
        $("#debug").removeClass("disabled");
    else
        $("#debug").addClass("disabled");
    this.masterscene.debug_mode_toggled.dispatch( this.masterscene.enable_debug );
}

