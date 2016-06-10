// Dependencies

// Code 

function SBSS_History() {
    this.command_queue = []
    this.current_command = 0;
    this.verbose = false;
}

SBSS_History.prototype.Next = function() {
    var self = this;
    if( self.current_command >= self.command_queue.length ){
        return {}
        self.current_command = self.command_queue.length;
    }
    else{
        var cmd = self.command_queue[ self.current_command ]
        self.current_command += 1;
        return cmd
    }        
}

SBSS_History.prototype.Clear = function() {
    var self = this;
    self.command_queue.length = 0;
    self.current_command = 0;
}

SBSS_History.prototype.Push = function(new_command) {
    var self = this
    if( 'command' in new_command && 'data' in new_command )
        self.command_queue.push( new_command )
}

SBSS_History.prototype.CommandsRemaining = function(){
    var self = this;
    if( self.current_command >= self.command_queue.length ){
        self.current_command = self.command_queue.length;
        return 0;
    }
    else{
        return self.command_queue.length - self.current_command;
    }
}

SBSS_History.prototype.List = function(){
    var self = this;
    if( self.verbose )
        for( command in self.command_queue ){
            if( self.current_command == command )
                console.log( command, " --> ", self.command_queue[command].command )
            else if( self.current_command > command )
                console.log( command, "  X  ", self.command_queue[command].command )        
            else
                console.log( command, "     ", self.command_queue[command].command )        
        }

}

// Exports

module.exports = SBSS_History
