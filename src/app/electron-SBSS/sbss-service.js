//Dependenices

SBSS_Server = require('./sbss-server.js')
SBSS_Director = require('./sbss-director.js')

// Code

function SBSS_Service(){
    this.server = new SBSS_Server();
    this.director = new SBSS_Director();
}


SBSS_Service.prototype.Initialize = function( options ){
    var self = this;
    self.director.Initialize();
    self.server.Initialize({
        port:options.port,
        director: self.director
    });
};


// Exports

module.exports = SBSS_Service;

