//Dependenices

SBSS_Server = require('./sbss-server.js')
SBSS_Director = require('./sbss-director.js')
SBSS_Simulation = require('./sbss-simulation.js')

// Code

function SBSS_Service(){
    this.server = new SBSS_Server();
    this.director = new SBSS_Director();
    this.simulation = new SBSS_Simulation();
}


SBSS_Service.prototype.Initialize = function( options ){
    var self = this;
    self.simulation.Initialize({
        server: self.server
    });
    self.director.Initialize({
        simulation: self.simulation,
        server: self.server
    });
    self.server.Initialize({
        port:options.port,
        director: self.director
    });
};


// Exports

module.exports = SBSS_Service;

