// Dependencies

// Code 

function TriangulatedSurfaceObject(){
    this.vertices = undefined;
    this.triangles = undefined;
    this.uv_triangles = undefined;
    this.normals = undefined;
    this.uv = undefined;
    this.loaded = false;
}

TriangulatedSurfaceObject.prototype.LoadFromBuffer = function(buffer) {

    vertices = []
    triangles = []
    uv_triangles = []
    normals = []
    uv = []
    this.loaded = false;
    this.vertices = undefined;
    this.triangles = undefined;
    this.uv_triangles = undefined;
    this.normals = undefined;
    this.uv = undefined;

    
    var lines = buffer.split(/\s*\r?\n\s*/);
    while( lines.length > 0 ){
        var line = lines.shift();
        if( line == '' )
            continue;

        var parts = line.split(' ');
        switch( parts[0] ){
        case 'v':
            parts.shift();
            if( parts.length < 3 ){
                console.log( "Failed to parse vertex line:" )
                console.log( parts )
                break;
            }
            vertices.push( parts );
            break;            
        case 'vt':
            parts.shift();
            if( parts.length < 2 ){
                console.log( "Failed to parse uv line:" )
                console.log( parts )
                break;
            }
            uv.push( parts );            
            break;            
        case 'vn':
            parts.shift();
            if( parts.length < 3 ){
                console.log( "Failed to parse normals line:" )
                console.log( parts )
                break;
            }
            normals.push( parts );           
            break;            
        case 'f':
            parts.shift();
            if( parts.length != 3 ){
                console.log( "Failed to parse face line:" )
                console.log( parts )
                break;
            }
            face_vertices = []
            face_uv = []
            for( p in parts ){
                fragments = parts[p].split('/');
                face_vertices.push( fragments.shift()-1 ) // Grab the first identifier as the vertex;
                face_uv.push( fragments.shift()-1 ) // Grab the second identifier as the uv;
            }
            triangles.push( face_vertices );
            uv_triangles.push( face_uv );
            break;
            
        default:
            console.log( "Unknown type flag:", parts[0] );
            break;            
        }
    }

    console.log( "Vertices: ", vertices.length );
    console.log( "Triangles: ", triangles.length );
    console.log( "UV Triangles; ", uv_triangles.length );
    console.log( "VTexture: ", uv.length );
    console.log( "VNormal: ", normals.length );

    if( this._VerifyData( vertices, triangles, uv_triangles, uv, normals ) ){
        this.vertices = vertices;
        this.triangles = triangles;
        this.uv_triangles = uv_triangles;
        this.uv = uv;
        this.normals = normals;
        this.loaded = true;  
    }
        
    return this.loaded;
}

TriangulatedSurfaceObject.prototype._VerifyData = function( v, f, uvf, vt, vn ){
    // There is much to be done here... For now we check some simple stuff...

    if( f.length != uvf.length )
        return false;

    for( triangle in f )
        for(node in f[triangle])
            if( f[triangle][node] >= v.length || f[triangle][node] < 0 )
                return false;

    for( triangle in uvf )
        for(node in uvf[triangle])
            if( uvf[triangle][node] >= vt.length || uvf[triangle][node] < 0 )
                return false;


    return true;
}  


TriangulatedSurfaceObject.prototype.Flatten = function(key){
    if( !this.loaded )
        return [];

    switch( key ){
    case 'vertex':
        data = []
        for( v in this.vertices )
            for( d in [0,1,2] )
                data.push( this.vertices[v][d] );        
        return data;

    case 'topology':
        data = []
        for( t in this.triangles )
            for( d in [0,1,2] )
                data.push( this.triangles[t][d] );        
        return data;

    case 'uv':
        data = []
        for( v in this.uv )
            for( d in [0,1] )
                data.push( this.uv[v][d] );        
        return data;

    default:
        return [];
    }
}

// Exports
module.exports = TriangulatedSurfaceObject
