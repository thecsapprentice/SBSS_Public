// Dependenices
var request = require('request');

// Code

function SBSS_Director(){


}


SBSS_Director.prototype.Initialize = function(){

}

SBSS_Director.prototype.ProcessCommandMessage = function( command, data ){
    console.log( "Handling command '%s'...", command );
    console.log( "With data:")
    console.log( data )
    
    switch( command ){
    case 'loadHistory':
        request('http://localhost:8081/history/' + data.name, function (error, response, body) {
            if (!error && response.statusCode == 200) {
                console.log(body) // Show the HTML for the Google homepage.
            }});



        break

    }
    
}


// Exports

module.exports = SBSS_Director
