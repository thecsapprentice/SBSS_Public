var CLE = require('CLEjs');
var fs = require('fs');
var async = require('async');

var cle = new CLE.CLEjs();


var verticesIn = fs.createReadStream('model.vertices.bin');
var topologyIn = fs.createReadStream('model.topology.bin');

function ReadFromFile( item, callback ){
    var collected_data = []
    item.on('data', function( chunk ){
        collected_data.push( chunk );
    });
    item.on('end', function(){
        var combined_data_length = 0;
        for( chunk in collected_data )
            combined_data_length += collected_data[chunk].length;
        var combined_data = new Uint8Array(combined_data_length );
        var current_byte = 0;
        for( chunk in collected_data ){
            combined_data.set( collected_data[chunk], current_byte );
            current_byte += collected_data[chunk].length;
        }
        callback( null, combined_data );
    });    
}

async.map( [verticesIn, topologyIn], ReadFromFile, function( err, results ){
    if( err ){            
        console.log("Failure!" );
    }
    else{
        console.log( "Sucess!" );
        var vertices = new Float32Array(results[0].buffer);
        var topology = new Uint32Array(results[1].buffer);

        cle.Set_Hook_Stiffness(1e4);
        cle.Set_Suture_Stiffness(1e5);
        cle.Set_Poissons_Ratio(0.45);
        cle.Set_Youngs_Modulus(1e3);
        cle.Create_Model(new Uint8Array(vertices.buffer),new Uint8Array(topology.buffer),0.05,10);    
        cle.Set_Fixed_Points([0,0,0]);
        cle.Finalize_Initialization();

        var hook_id = cle.Add_Hook( [-0.0353868, 0.025, -0.0229982] );
        cle.Move_Hook( hook_id, [1,1,1] );

        setInterval( cle.Advance_One_Time_Step.bind(cle), 33 /* 30fps */ );
        //cle.Advance_One_Time_Step();
        //cle.WaitForSolve();
    }
})
