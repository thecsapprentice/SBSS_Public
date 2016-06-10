'use strict';

const electron = require('electron');
// Module to control application life.
const app = electron.app;
// Module to create native browser window.
const BrowserWindow = electron.BrowserWindow;
var path           = require('path')
var fs             = require('fs')
var url            = require('url')
var http           = require('http');
var querystring    = require('querystring')
// Keep a global reference of the window object, if you don't, the window will
// be closed automatically when the JavaScript object is garbage collected.
let mainWindow;
var express = require('express');

var SBSS_Service = require('./sbss-service.js')
let sbss_service;
let http_server;



function LoadAvailableHistories(callback){
    fs.readdir( path.join( path.normalize(__dirname),'data','history'), function(err, history_files){
        console.log( history_files )
        var good_history = []
        for( var i = 0; i < history_files.length; i++ ){
            if(path.extname(history_files[i]) == '.hst')
                good_history.push( path.basename( history_files[i], '.hst' ) )
        }
        callback({mimeType: 'text/json', data: new Buffer(JSON.stringify( { histories: good_history })) });
    });
}

function LoadAvailableScenes(callback){
    fs.readdir( path.join( path.normalize(__dirname),'data','scenes'), function(err, scene_files){
        console.log( scene_files )
        var good_scene = []
        for( var i = 0; i < scene_files.length; i++ ){
            if(path.extname(scene_files[i]) == '.smd')
            good_scene.push( path.basename( scene_files[i], '.smd' ) )
        }
        callback({mimeType: 'text/json', data: new Buffer(JSON.stringify( { scenes: good_scene } ))});
    });
}

function createWindow () {
   
    var protocol       = electron.protocol;

    http_server = express();
    http_server.get( '/history/:name', function( req, res ){
        console.log( "Accessing history file: ", req.params.name );
        var filePath = path.join(__dirname, "data", "history", req.params.name+".hst");
        fs.readFile(filePath, {encoding: 'utf-8'}, function(err,data){
            if (!err){
                var history_cmds = JSON.parse(data);
                res.json( history_cmds);
                res.end();
            }else{
                console.log(err);
                res.status(500).end();
            }            
        });
    });
    http_server.get( '/scene/:name', function( req, res ){
        console.log( "Accessing scene file: ", req.params.name );
        var filePath = path.join(__dirname, "data", "scenes", req.params.name+".smd");
        fs.readFile(filePath, {encoding: 'utf-8'}, function(err,data){
            if (!err){
                var scene_cmds = JSON.parse(data);
                res.json( scene_cmds);
                res.end();
            }else{
                console.log(err);
                res.status(500).end();
            }            
        });
    });
    http_server.get( '/model/:name', function( req, res ){
        console.log( "Accessing model file: ", req.params.name );
        var filePath = path.join(__dirname, "data", "models", req.params.name);
        fs.readFile(filePath, {encoding: 'utf-8'}, function(err,data){
            if (!err){                
                res.send( data );
                res.end();
            }else{
                console.log(err);
                res.status(500).end();
            }            
        });
    });
    http_server.get( '/texture/:name', function( req, res ){
        console.log( "Accessing texture file: ", req.params.name );
        var filePath = path.join(__dirname, "data", "textures", req.params.name);
        fs.readFile(filePath, function(err,data){
            if (!err){
                res.send( data );
                res.end();
            }else{
                console.log(err);
                res.status(500).end();
            }            
        });
    });

    http_server.listen(8081,function(){
        console.log("Server listening on: http://localhost:%s", 8081);
    });
    
    sbss_service = new SBSS_Service();
    sbss_service.Initialize({
        port: 8080
    });
    
    protocol.registerStandardSchemes(['sbss-data']);
    
    protocol.registerBufferProtocol('sbss-data', function(req, callback) {
        
        var uri = url.parse(req.url)
        console.log( uri )
        switch( uri.host ){
        case 'data':
            if( uri.path == '/histories' ){
                console.log( "Loading Histories..." )
                LoadAvailableHistories(callback)
            }
            if( uri.path == '/scenes' ){
                console.log( "Loading Scenes..." )
                LoadAvailableScenes(callback)
            }
            break;
        default:
            console.log( "Unknown Resource Indentifier: %s", uri.host );
            break;
        }
       
                
    })
    
    
  // Create the browser window.
  mainWindow = new BrowserWindow({width: 800, height: 600});

  // and load the index.html of the app.
  mainWindow.loadURL('file://' + __dirname + '/client.html');

  // Open the DevTools.
  mainWindow.webContents.openDevTools();

  // Emitted when the window is closed.
  mainWindow.on('closed', function() {
    // Dereference the window object, usually you would store windows
    // in an array if your app supports multi windows, this is the time
    // when you should delete the corresponding element.
    mainWindow = null;
  });


    
}

// This method will be called when Electron has finished
// initialization and is ready to create browser windows.
app.on('ready', createWindow);

// Quit when all windows are closed.
app.on('window-all-closed', function () {
  // On OS X it is common for applications and their menu bar
  // to stay active until the user quits explicitly with Cmd + Q
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

app.on('activate', function () {
  // On OS X it's common to re-create a window in the app when the
  // dock icon is clicked and there are no other windows open.
  if (mainWindow === null) {
    createWindow();
  }
});
