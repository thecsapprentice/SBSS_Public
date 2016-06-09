// Encapsulates all websocket communications
//
//
//
//
//

var Communications = function(){
    this.websocket = null;
    this.connection_attempts = 1;
    this.reconnect_url = "";
    this.reconnect_time = 0;
    this.max_reconnect_attempts = 10;

    // Signals
    this.opened = new signals.Signal();
    this.messageReceived = new signals.Signal();
    this.closed = new signals.Signal();
    this.givingup = new signals.Signal();
}

Communications.prototype.connect = function( url, reconnect ){
    reconnect = typeof reconnect !== 'undefined' ? reconnect : false;
    this.websocket = new WebSocket( url );
    
    this.websocket.addEventListener("open", this.connectionOpen.bind(this));    
    this.websocket.addEventListener("close", this.connectionClose.bind(this,reconnect));    
    this.websocket.addEventListener("error", this.connectionError.bind(this));    
    this.websocket.addEventListener("message", this.connectionMessage.bind(this));    

    if(reconnect)
        this.reconnect_url = url;
}

Communications.prototype.status = function() {
    switch (socket.readyState) {
    case WebSocket.CONNECTING:
        return "CONNECTING";
        break;
    case WebSocket.OPEN:
        return "OPEN";
        break;
    case WebSocket.CLOSING:
        return "CLOSING";
        break;
    case WebSocket.CLOSED:
        return "CLOSED";
        break;
    default:
        return "UNKNOWN";
        break;
    }
}

Communications.prototype.connectionOpen = function(event) {
    console.log('CONNECTED');
    this.connection_attempts = 1;    
    this.opened.dispatch();
}

Communications.prototype.connectionClose = function(event, reconnect) {
    this.closed.dispatch();
    if(reconnect && this.max_reconnect_attempts > this.connection_attempts){
        this.reconnect_time = this.generateInterval( this.connection_attempts );
        console.log('DISCONNECTED - Reconnecting in ' + this.reconnect_time + ' seconds.');        
        setTimeout( this.reconnect.bind(this),  this.reconnect_time );
    }
    else{
        console.log('DISCONNECTED.');       
        //this.givingup.dispatch();
    }
}

Communications.prototype.connectionError = function(event) {
    console.log( "Websocket Error: " + event.data);
}

Communications.prototype.connectionMessage = function(event) {
    if (!event.data) {
        console.log("PING");
    } else {
        this.messageReceived.dispatch(JSON.parse(event.data));
    }
}

Communications.prototype.reconnect = function() {
    this.connection_attempts ++;
    this.connect(this.reconnect_url, true)
}
    

Communications.prototype.teardown = function() {
    this.websocket.send("exit");
    this.websocket.close();
    this.websocket = null;
}

Communications.prototype.generateInterval = function(k) {
    var max_wait = 2;
    var maxInterval = (Math.pow(2, k) - 1) * 1000;
    if (maxInterval > max_wait*1000) {
        maxInterval = max_wait*1000; 
    }    
    // generate the interval to a random number between 0 and the maxInterval determined from above
    return Math.random() * maxInterval; 
}

Communications.prototype.send = function(message){
    //if(command.length>0){
    //    if(first_rush){
    //        alert("Wait cursor indicates processing occuring. Client commands ignored during this time.");
    //        first_rush = false;
    //    }
    //    return;
    //}
    if(this.websocket.readyState===1 )  {
        this.websocket.send(message);
    }
    else
        console.log('WebSocket not open yet. Message ignored.');
}
