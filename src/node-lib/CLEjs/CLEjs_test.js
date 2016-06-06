CLEjs = require('CLEjs.node').CLEjs

vertices = new Float32Array( 8 * 3 );
triangles = new Uint32Array( 12 * 3 );

vertices[0] = -0.500000; vertices[1] = -0.500000;  vertices[2] =  0.500000;
vertices[3] =  0.500000; vertices[4] = -0.500000;  vertices[5] =  0.500000;
vertices[6] = -0.500000; vertices[7] =  0.500000;  vertices[8] =  0.500000;
vertices[9] =  0.500000; vertices[10] =  0.500000;  vertices[11] =  0.500000;
vertices[12] = -0.500000; vertices[13] =  0.500000;  vertices[14] = -0.500000;
vertices[15] =  0.500000; vertices[16] =  0.500000;  vertices[17] = -0.500000;
vertices[18] = -0.500000; vertices[19] = -0.500000;  vertices[20] = -0.500000;
vertices[21] =  0.500000; vertices[22] = -0.500000;  vertices[23] = -0.500000;

triangles[0] = 0; triangles[1] = 1; triangles[2] = 2;
triangles[3] = 2; triangles[4] = 1; triangles[5] = 3;

triangles[6] = 2; triangles[7] = 3; triangles[8] = 4;
triangles[9] = 4; triangles[10] = 3; triangles[11] = 5;

triangles[12] = 4; triangles[13] = 5; triangles[14] = 6;
triangles[15] = 6; triangles[16] = 5; triangles[17] = 7;

triangles[18] = 6; triangles[19] = 7; triangles[20] = 0;
triangles[21] = 0; triangles[22] = 7; triangles[23] = 1;

triangles[24] = 1; triangles[25] = 7; triangles[26] = 3;
triangles[27] = 3; triangles[28] = 7; triangles[29] = 5;

triangles[30] = 6; triangles[31] = 0; triangles[32] = 4;
triangles[33] = 4; triangles[34] = 0; triangles[35] = 2;

cle = CLEjs();
cle.Create_Model( new Uint8Array(vertices.buffer), new Uint8Array(triangles.buffer), 0.1, 10);
