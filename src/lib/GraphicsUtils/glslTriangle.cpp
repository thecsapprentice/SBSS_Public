// glslTriangle.cpp
// Author: Court Cutting
// Date: 1/6/2012
// Purpose: Triangle object management class with full topological services which
//       makes maximal use of glsl, speeding graphics and unloading CPU maximally.
//       Copyright 2012 - All rights reserved.

#pragma warning(disable : 4996)
#include <fstream>
#include <set>
#include <assert.h>
#ifdef linux
#include <stdlib.h>
#include <string.h>
#endif
#include "math3d.h"
#include "Vec3f.h"
#include "glslTriangle.h"

//////////////////////// TEMPORARY TEMPORARY TEMPORARY - On SnowLeopard this is suppored, but GLEW doens't hook up properly
//////////////////////// Fixed probably in 10.6.3
#ifdef __APPLE__
#define glGenVertexArrays glGenVertexArraysAPPLE
#define glDeleteVertexArrays  glDeleteVertexArraysAPPLE
#define glBindVertexArray	glBindVertexArrayAPPLE
#endif

/* static const char *GTNormalMakerVertexShader = "#version 140\n"
	"uniform samplerBuffer positionSampler;\n"
	"uniform usamplerBuffer neighborSampler;\n"
	"out vec3 vSurfaceNormal;\n"
	"out vec4 vSurfaceCoord;\n"
	"void main(void) {\n"
	"	uvec4 nei4 = texelFetch(neighborSampler,gl_VertexID);\n"
	"	vec4 tmp,tPosition4 = texelFetch(positionSampler,int(nei4[0]));\n" // gl_VertexID
	"	vSurfaceCoord = tPosition4;\n"
	"	vSurfaceNormal[0] = 0.0f;\n"
	"	vSurfaceNormal[1] = 0.707f;\n"
	"	vSurfaceNormal[2] = 0.707f;\n"
		"	vSurfaceNormal[0] = float(gl_VertexID);\n"
		"	vSurfaceNormal[1] = float(nei4[0]);\n"
		"	vSurfaceNormal[2] = 7.7f;\n"
	"	vec3 norm,vtx,first,prev,now;\n"
	"	vtx = tPosition4.xyz;\n"
	"	norm = vec3(0.0f,0.0f,0.0f);\n"
	"	bool notDone=true,close=true;\n"
	"	if(bool(nei4[1]&0x80000000u)) {\n"
	"		close=false;\n"
	"		nei4[1] &= 0x7fffffffu; }\n"
	"	tmp = texelFetch(positionSampler,int(nei4[1]));\n"
	"	first = prev = tmp.xyz - vtx;\n"
	"	uint pIndx,txIndx;\n"
	"	int i=2;\n"
	"	txIndx = 0x80000000u;\n"
	"	while(notDone) {\n"
	"		if(bool(nei4[i]&0x80000000u)) {\n"
	"			notDone=false;\n"
	"			pIndx = nei4[i] & 0x7fffffffu; }\n"
	"		else {\n"
	"			if(i>2) {\n"
	"				i=0;\n"
	"				if(bool(txIndx&0x80000000u)) {\n"
	"					txIndx = nei4[3];\n"
	"					nei4 = texelFetch(neighborSampler,int(txIndx));\n"
	"					continue; }\n"
	"				else {\n"
	"					pIndx = nei4[i];\n"
	"					++txIndx;\n"
	"					nei4 = texelFetch(neighborSampler,int(txIndx)); } }\n"
	"			else {\n"
	"				pIndx = nei4[i]; ++i; } }\n"
	"		tmp = texelFetch(positionSampler,int(pIndx));\n" // gl_VertexID
	"		now = tmp.xyz - vtx;\n"
	"		norm += cross(now,prev);\n"
	"		prev = now;	}\n"
	"	if(close)\n"
	"		norm += cross(first,prev);\n"
	"	vSurfaceNormal = normalize(norm);\n"
	"}\n"; */

static const char *GTNormalMakerVertexShader = "#version 140\n"
	"uniform samplerBuffer positionSampler;\n"
	"uniform usamplerBuffer neighborSampler;\n"
	"out vec3 vSurfaceNormal;\n"
	"out vec4 vSurfaceCoord;\n"
	"void main(void) {\n"
	"	uvec4 nei4 = texelFetch(neighborSampler,gl_VertexID);\n"
	"	vec4 tmp,tPosition4 = texelFetch(positionSampler,int(nei4[0]));\n" // gl_VertexID
	"	vSurfaceCoord = tPosition4;\n"
	"	vec3 norm,vtx,first,prev,now;\n"
	"	vtx = tPosition4.xyz;\n"
	"	norm = vec3(0.0f,0.0f,0.0f);\n"
	"	bool notDone=true,close=true;\n"
	"	if(bool(nei4[1]&0x80000000u)) {\n"
	"		close=false;\n"
	"		nei4[1] &= 0x7fffffffu; }\n"
	"	tmp = texelFetch(positionSampler,int(nei4[1]));\n"
	"	first = prev = tmp.xyz - vtx;\n"
	"	uint pIndx,txIndx;\n"
	"	int i=2;\n"
	"	txIndx = 0x80000000u;\n"
	"	while(notDone) {\n"
	"		if(bool(nei4[i]&0x80000000u)) {\n"
	"			notDone=false;\n"
	"			pIndx = nei4[i] & 0x7fffffffu; }\n"
	"		else {\n"
	"			if(i>2) {\n"
	"				i=0;\n"
	"				if(bool(txIndx&0x80000000u)) {\n"
	"					txIndx = nei4[3];\n"
	"					nei4 = texelFetch(neighborSampler,int(txIndx));\n"
	"					continue; }\n"
	"				else {\n"
	"					pIndx = nei4[i];\n"
	"					++txIndx;\n"
	"					nei4 = texelFetch(neighborSampler,int(txIndx)); } }\n"
	"			else {\n"
	"				pIndx = nei4[i]; ++i; } }\n"
	"		tmp = texelFetch(positionSampler,int(pIndx));\n" // gl_VertexID
	"		now = tmp.xyz - vtx;\n"
	"		norm += cross(now,prev);\n"
	"		prev = now;	}\n"
	"	if(close)\n"
	"		norm += cross(first,prev);\n"
	"	vSurfaceNormal = normalize(norm);\n"
	"}\n";

void glslTriangle::deleteTriangle(long triangleNumber)
{	// may strand vertices
	for(int i=triangleNumber*3; i<_tris.size()-3; ++i)
		_tris[i] = _tris[i+3];
	_tris.pop_back(); _tris.pop_back(); _tris.pop_back();
	--_triangleNumber;
}

long glslTriangle::addTriangle(unsigned long (&vertices)[3])
{
	_tris.push_back((GLuint)vertices[0]);
	_tris.push_back((GLuint)vertices[1]);
	_tris.push_back((GLuint)vertices[2]);
	_triangleNumber = (unsigned long)(_tris.size()/3);
	if(_adjacenciesComputedP)	{ // done for incision class to make space
		// Really _adjacenciesComputedP should be set false when a new triangle is added.
		_adjsP.push_back(0x00000003);
		_adjsP.push_back(0x00000003);
		_adjsP.push_back(0x00000003); }
	return (long)_triangleNumber-1;
}

long glslTriangle::addCoordNormTexVertices(int numberToAdd)
{
	long retval = _vertexNumber;
	for(unsigned int i=0; i<(unsigned)numberToAdd; ++i)
	{
		_pnI.push_back(positionNormalIndex());
		_pnI.back().posIdx=(long)_Positions.size();
		_pnI.back().posNormIdx=_uniqueNormalNumber+i;
		_Positions.push_back(0.0f);
		_Positions.push_back(0.0f);
		_Positions.push_back(0.0f);
		_Positions.push_back(1.0f);
		_TexCoords.push_back(0.0f);
		_TexCoords.push_back(0.0f);
		if(_adjacenciesComputed)
			_vertexFace.push_back(0x80000000);
		if(_adjacenciesComputedP)
			_vertexFaceP.push_back(0x80000000);
	}
	_vertexNumber += numberToAdd;
	_uniqueNormalNumber += numberToAdd;
	return retval;
}

long glslTriangle::addNormTexVertex(int existingVertex)
{
	if(existingVertex<0 || existingVertex>=(int)_vertexNumber)
		return -1;
	_TexCoords.push_back(0.0f);
	_TexCoords.push_back(0.0f);
	_pnI.push_back(positionNormalIndex());
	_pnI.back().posIdx=_pnI[existingVertex].posIdx;
	_pnI.back().posNormIdx=_uniqueNormalNumber;
	++_uniqueNormalNumber;
	if(_adjacenciesComputed)
		_vertexFace.push_back(_vertexFace[_pnI[existingVertex].posNormIdx]);
	if(_adjacenciesComputedP)
		_vertexFaceP.push_back(_vertexFaceP[_pnI[existingVertex].posIdx>>2]);
	return _vertexNumber++;
}

long glslTriangle::addTexVertex(int existingVertex)
{	// returns the index to the added vertex which shares a coordinate and a normal with input vertex
	if(existingVertex<0 || existingVertex>=(int)_vertexNumber)
		return -1;
	_TexCoords.push_back(0.0f);
	_TexCoords.push_back(0.0f);
	_pnI.push_back(positionNormalIndex());
	_pnI.back().posIdx=_pnI[existingVertex].posIdx;
	_pnI.back().posNormIdx=_pnI[existingVertex].posNormIdx;
	if(_adjacenciesComputed)
		_vertexFace.push_back(_vertexFace[_pnI[existingVertex].posNormIdx]);
	if(_adjacenciesComputedP)
		_vertexFaceP.push_back(_vertexFaceP[_pnI[existingVertex].posIdx>>2]);
	return _vertexNumber++;
}

long glslTriangle::addNormTexVertex2(int existingPositionIndex)
{	// input is an absolute index and should not be divided by 4
	if(existingPositionIndex<0 || existingPositionIndex>=(int)_vertexNumber || (existingPositionIndex&0x00000003))
		return -1;
	_TexCoords.push_back(0.0f);
	_TexCoords.push_back(0.0f);
	_pnI.push_back(positionNormalIndex());
	_pnI.back().posIdx=existingPositionIndex;
	_pnI.back().posNormIdx=_uniqueNormalNumber;
	++_uniqueNormalNumber;
	if(_adjacenciesComputed)	{	// NOT EFFICIENT. Set false before using this routine
		int i=0;
		while(i<_pnI.size())
			if(_pnI[i].posIdx==existingPositionIndex)
				break;
		if(i==_pnI.size())	{	// bad input
			_TexCoords.pop_back();  _TexCoords.pop_back();
			_pnI.pop_back();
			--_uniqueNormalNumber;
			return -1;
		}
		_vertexFace.push_back(_vertexFace[i]);
	}
	if(_adjacenciesComputedP)
		_vertexFaceP.push_back(_vertexFaceP[existingPositionIndex>>2]);
	return _vertexNumber++;
}

void glslTriangle::createNormalMakerProgram()
{
	if(_NormalMakerProgram>0)
		glDeleteProgram(_NormalMakerProgram);
	GLuint normalMakerShader;   
	GLint testVal;
	// Create shader objects
	normalMakerShader = glCreateShader(GL_VERTEX_SHADER);
	GLchar *fsStringPtr[1];
	fsStringPtr[0] = (GLchar *)GTNormalMakerVertexShader;
	glShaderSource(normalMakerShader, 1, (const GLchar **)fsStringPtr, NULL);
	glCompileShader(normalMakerShader);
	glGetShaderiv(normalMakerShader, GL_COMPILE_STATUS, &testVal);
	if(testVal == GL_FALSE){
		GLchar infoLog[2000];
		GLsizei infoLength;
		glGetShaderInfoLog(normalMakerShader,2000,&infoLength,infoLog);
		printf("Normal maker shader failed. Info:\n%s",infoLog);
		glDeleteShader(normalMakerShader);
		return;
	}
	_NormalMakerProgram = glCreateProgram();
	glAttachShader(_NormalMakerProgram, normalMakerShader);
	const char* varying_names[]={"vSurfaceNormal", "vSurfaceCoord"};
	glTransformFeedbackVaryings(_NormalMakerProgram, 2, varying_names, GL_SEPARATE_ATTRIBS); 
	glLinkProgram(_NormalMakerProgram);
	// These are no longer needed
	glDeleteShader(normalMakerShader);
	// Make sure link worked too
	glGetProgramiv(_NormalMakerProgram, GL_LINK_STATUS, &testVal);
	if(testVal == GL_FALSE){
		GLchar infoLog[2000];
		GLsizei infoLength;
		glGetProgramInfoLog(_NormalMakerProgram,2000,&infoLength,infoLog);
		printf("Normal maker shader failed. Info:\n%s",infoLog);
		glDeleteProgram(_NormalMakerProgram);
	}
}

glslTriangle::glslTriangle(bool texturedNotColored, bool staticNotMorphing) {
	sceneNode::_coloredNotTextured=!texturedNotColored;
	setType(TRIANGLES);
	_textured=texturedNotColored;
	_static=staticNotMorphing;
	_deformable=false;
	_Indexes.clear();
	_vertexNumber = 0;
	_uniqueNormalNumber = 0;
	_Verts.clear();
	_Norms.clear();
	_TexCoords.clear();
	_TopoIndex.clear();
	_Positions.clear();
	_glslNeighbors.clear();
	_adjacenciesComputed=false;
	_adjacenciesComputedP=false;
    _nMaxIndexes = 0;
    _nNumIndexes = 0;
    _nNumVerts = 0;
	_vertexArrayBufferObject=0;
	for(int i=0; i<4; ++i)
		_bufferObjects[i]=0;
	_textureBufferObjects[0]=0;
	_textureBufferObjects[1]=0;
	_texBOBuffers[0]=0;
	_texBOBuffers[1]=0;
	_NormalMakerProgram=0;
    _transformFeedbackQuery=0;
}

glslTriangle::~glslTriangle() {
    // Just in case these still are allocated when the object is destroyed
	_Indexes.clear();
	_Verts.clear();
    _Norms.clear();
	_TexCoords.clear();
	_Positions.clear();
	_glslNeighbors.clear();

	if(_NormalMakerProgram)
		glDeleteProgram(_NormalMakerProgram);
	if(_bufferObjects[0]>0)
		glDeleteBuffers(4, _bufferObjects);
	if(_texBOBuffers[0]>0)
		glDeleteBuffers(2,_texBOBuffers);
	// Cleanup textures
	if(_textureBufferObjects[0]>0) {
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_BUFFER_ARB, 0);
		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_BUFFER_ARB, 0);
		glActiveTexture(GL_TEXTURE0);
		glDeleteTextures(2,_textureBufferObjects);
	}
	if(_vertexArrayBufferObject>0)
		glDeleteVertexArrays(1, &_vertexArrayBufferObject);
	if(_transformFeedbackQuery>0)
		glDeleteQueries(1,&_transformFeedbackQuery);
}

bool glslTriangle::getTriangleVertices(unsigned int triangle, int (&vertices)[3])
{
	int i=(triangle<<1)+triangle;
	if(i<_tris.size())	{
		vertices[0] = int(_tris[i]);
		vertices[1] = int(_tris[i+1]);
		vertices[2] = int(_tris[i+2]);
		return true;
	}
	else
		return false;
}

int glslTriangle::getClosestVertex(float (&position)[3], int triangle)
{
	int closeV=-1;
	float *v,d2min=1e30f,d2;
	if(triangle>-1)	{
		GLuint *tr=triVerts(triangle);
		for(int i=0; i<3; ++i) {
			v = vertexCoordinate(long(tr[i]));
			d2 = (*v-position[0])*(*v-position[0]);
			d2 += (v[1]-position[1])*(v[1]-position[1]);
			d2 += (v[2]-position[2])*(v[2]-position[2]);
			if(d2<d2min)	{
				d2min = d2;
				closeV = tr[i];
			}
		}
	}
	else	{
		for(unsigned int i=0; i<_vertexNumber; ++i) {
			v = vertexCoordinate(i);
			d2 = (*v-position[0])*(*v-position[0]);
			d2 += (v[1]-position[1])*(v[1]-position[1]);
			d2 += (v[2]-position[2])*(v[2]-position[2]);
			if(d2<d2min)	{
				d2min = d2;
				closeV = i;
			}
		}
	}
	return closeV;
}

bool glslTriangle::makeObject()	{	// COURT - this goes when debug finished.
	if(_NormalMakerProgram<1)	// do once for the class
		createNormalMakerProgram();
	// Create the buffer objects once
	if(!_vertexArrayBufferObject)
		glGenVertexArrays(1,&_vertexArrayBufferObject);
	if(!_bufferObjects[0])
	    glGenBuffers(4, _bufferObjects);
	// Create 2 new texture buffer objects
	if(!_textureBufferObjects[0])
		glGenTextures(2,_textureBufferObjects);
	if(!_texBOBuffers[0])
		glGenBuffers(2,_texBOBuffers);
    // Just in case this gets called more than once...
	_Indexes.clear();
	_Verts.clear();
	_Norms.clear();
	_TexCoords.clear();
	_TopoIndex.clear();
	_Positions.clear();
	_glslNeighbors.clear();
    
    _nMaxIndexes = 12;
    _nNumIndexes = 12;
    _nNumVerts = 4;
    
    // Allocate new blocks. In reality, the other arrays will be
    // much shorter than the index array
	GLfloat testVerts[] = {0.0f, 0.2f, -0.8f, 1.0f,
		0.0f, 1.8f, -0.8f, 1.0f,
		0.8f, 1.0f, 0.8f, 1.0f,
		-0.8f, 1.0f, 0.8f, 1.0f};
	_Verts.assign(16,0.0f);
//	_Verts.assign(testVerts,testVerts+16);
	_Positions.assign(testVerts,testVerts+16);
	GLuint testIndices[]= {2, 1, 3,
		1, 0, 3,
		0, 1, 2,
		0, 2, 3};
	_Indexes.assign(testIndices,testIndices+12);
//	GLfloat testNorms[] = {-0.2f, -0.2f, -0.2f,
//		0.2f, -0.2f, -0.2f,
//		0.0f, -0.5f, 0.5f,
//		0.0f, 1.0f, 0.0f};
    _Norms.assign(12,0.0f);
//    _Norms.assign(testNorms,testNorms+12);
	GLfloat testTex[] = {0.0f, 0.0f,
		0.0f, 1.0f,
		1.0f, 0.0f,
		1.0f, 1.0f};
    _TexCoords.assign(testTex,testTex+8);
	GLuint testNei[] = {0, 1, 3, 2,
		1, 3,  0, 2,
		2, 1, 0, 3,
		3, 0, 1, 2};
	_glslNeighbors.assign(testNei,testNei+16);
    // Copy data to video memory
    // Vertex data
    glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[0]);	// VERTEX_DATA
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_nNumVerts*4, &(_Verts[0]), GL_DYNAMIC_DRAW);
    // Normal data
    glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[1]);	// NORMAL_DATA
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_nNumVerts*3, &(_Norms[0]), GL_DYNAMIC_DRAW);
    // Texture coordinates
    glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[2]);	// TEXTURE_DATA
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_nNumVerts*2, &(_TexCoords[0]), GL_STATIC_DRAW);
    // Indexes
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bufferObjects[3]);	// INDEX_DATA
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*_Indexes.size(), &(_Indexes[0]), GL_STATIC_DRAW);
	// release for next use
	glBindBuffer( GL_ARRAY_BUFFER, 0);
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0);

	// Create the master vertex array object
	glBindVertexArray(_vertexArrayBufferObject);
    // Vertex data
    glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[0]);	// VERTEX_DATA
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
//	glEnableClientState( GL_VERTEX_ARRAY);
    // Normal data
    glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[1]);	// NORMAL_DATA
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
    // Texture coordinates
    glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[2]);	// TEXTURE_DATA
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, 0);
    // Indexes
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bufferObjects[3]);	// INDEX_DATA
	// release for next use
//    glBindBuffer( GL_ARRAY_BUFFER, 0);
//    glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0);

    // Unbind to anybody
	glBindVertexArray(0);

//	glBindVertexArray(_vertexArrayBufferObject[1]);

	// Texture coordinates
//    glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[2]);	// TEXTURE_DATA
//	glEnableVertexAttribArray(0);	// ENUM 2 GLT_ATTRIBUTE_TEXTURE0
//	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, 0);	// ENUM 2 GLT_ATTRIBUTE_TEXTURE0

	// transform feedback buffers
	glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER,0,_bufferObjects[1]);
//    glBufferData(GL_TRANSFORM_FEEDBACK_BUFFER, sizeof(GLfloat)*_nNumVerts*3, NULL, GL_DYNAMIC_COPY);
	glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER,1,_bufferObjects[0]);
//    glBufferData(GL_TRANSFORM_FEEDBACK_BUFFER, sizeof(GLfloat)*_nNumVerts*4, NULL, GL_DYNAMIC_COPY);

//	glBindVertexArray(0);

	// COURT - probably don't do this next paragraph.
    // Free older, larger arrays
	_Indexes.clear();
	_Verts.clear();
	_Norms.clear();
	_TexCoords.clear();
	_TopoIndex.clear();

	// Load the _Positions first
	glBindBuffer(GL_TEXTURE_BUFFER_ARB, _texBOBuffers[0]);
	glBufferData(GL_TEXTURE_BUFFER_ARB, sizeof(GLfloat)*4*_nNumVerts, &(_Positions[0]), GL_DYNAMIC_DRAW);
//	glBindBuffer(GL_TEXTURE_BUFFER_ARB, 0);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_BUFFER_ARB, _textureBufferObjects[0]);
	glTexBufferARB(GL_TEXTURE_BUFFER_ARB, GL_RGBA32F, _texBOBuffers[0]); 

	// Load thecounterclockwise neighbor info second
	glBindBuffer(GL_TEXTURE_BUFFER_ARB, _texBOBuffers[1]);
	glBufferData(GL_TEXTURE_BUFFER_ARB, sizeof(GLuint)*4*_nNumVerts, &(_glslNeighbors[0]), GL_STATIC_DRAW);
//	glBindBuffer(GL_TEXTURE_BUFFER_ARB, 0);
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_BUFFER_ARB, _textureBufferObjects[1]);
	glTexBufferARB(GL_TEXTURE_BUFFER_ARB, GL_RGBA32UI, _texBOBuffers[1]); 

	glActiveTexture(GL_TEXTURE0);

//    glBindBuffer( GL_TEXTURE_BUFFER_ARB, 0);

	return true;
}

// Draw - make sure you call glEnableClientState for these arrays
void glslTriangle::draw(void) 
{
	// bind the texture buffer objects containing the normals and the vertices
	glBindVertexArray(_vertexArrayBufferObject);
	// assumes texturedShader already in use
	if(_textured) {
		glActiveTexture( GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, _2DtextureBufferNumber);
	}
	//	assumes glUseProgram(_program) has already been called
	glDrawElements(GL_TRIANGLES, (_triangleNumber<<1)+_triangleNumber, GL_UNSIGNED_INT, 0);
    // Unbind to anybody
	glBindVertexArray(0);

GLenum errCode;
if((errCode=glGetError())!=GL_NO_ERROR)
	errCode=errCode;

}    

void glslTriangle::setVerticesComputeNormals()
{
	if(_static)	{	// not frequently morphing
		// Vertex data
	    glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[0]);	// VERTEX_DATA
		GLfloat *vtx,*fptr;
		fptr = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER,GL_WRITE_ONLY);	//|GL_MAP_INVALIDATE_BUFFER_BIT|GL_MAP_FLUSH_EXPLICIT_BIT);
		size_t n = _pnI.size();
		for(unsigned int j,i=0; i<n; ++i)	{
			vtx = &(_Positions[_pnI[i].posIdx]);
			for(j=0; j<4; ++j)
				*(fptr++) = *(vtx++);
		}
		glUnmapBuffer(GL_ARRAY_BUFFER);
/*		glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[0]);
		std::vector<GLfloat> fbCoords;
		fbCoords.assign(_vertexNumber*4,0.0f);
		glGetBufferSubData(GL_ARRAY_BUFFER,0,sizeof(GLfloat)*_vertexNumber*4,&(fbCoords[0])); */
		glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[1]);
		fptr = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER,GL_WRITE_ONLY);	//|GL_MAP_INVALIDATE_BUFFER_BIT|GL_MAP_FLUSH_EXPLICIT_BIT);
		float normal[3];
		for(unsigned int j,i=0; i<n; ++i)	{
			getVertexNormal(i,normal);
			for(j=0; j<3; ++j)
				*(fptr++) = normal[j];
		}
		glUnmapBuffer(GL_ARRAY_BUFFER);
/*		glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[1]);
		fbCoords.assign(_vertexNumber*3,0.0f);
		glGetBufferSubData(GL_ARRAY_BUFFER,0,sizeof(GLfloat)*_vertexNumber*3,&(fbCoords[0])); */
		glBindBuffer(GL_ARRAY_BUFFER,0);
	}
	else	{	// morphing. Load positions and compute normals on graphics card
		// Load the current _Positions. Assumes topological array already loaded
		glBindBuffer(GL_TEXTURE_BUFFER_ARB, _texBOBuffers[0]);
		glBufferData(GL_TEXTURE_BUFFER_ARB, sizeof(GLfloat)*_Positions.size(), &(_Positions[0]), GL_DYNAMIC_DRAW);
		glBindVertexArray(0);
		glUseProgram(_NormalMakerProgram);
		// transform feedback buffers
		glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER,0,_bufferObjects[1]);
	/*	GLint err = glGetError();
		if(err==GL_NO_ERROR)
			int i=0;
		else if(err==GL_INVALID_ENUM)
			int i=0;
		else if(err==GL_INVALID_VALUE)
			int i=0;
		else
			int i=0; */
		glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER,1,_bufferObjects[0]);
	/*	err = glGetError();
		if(err==GL_NO_ERROR)
			int i=0;
		else if(err==GL_INVALID_ENUM)
			int i=0;
		else if(err==GL_INVALID_VALUE)
			int i=0;
		else
			int i=0; */
		glUniform1i( glGetUniformLocation( _NormalMakerProgram, "positionSampler"), 1);
		glUniform1i( glGetUniformLocation( _NormalMakerProgram, "neighborSampler"), 2);
		// bind the texture buffer objects containing the normals and the vertices
		glActiveTexture( GL_TEXTURE1);
		glBindTexture( GL_TEXTURE_BUFFER_EXT, _textureBufferObjects[0]);
		glActiveTexture( GL_TEXTURE2);
		glBindTexture( GL_TEXTURE_BUFFER_EXT, _textureBufferObjects[1]);
		glEnable(GL_RASTERIZER_DISCARD);
		glGenQueries(1,&_transformFeedbackQuery);
		glBeginQuery(GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN, _transformFeedbackQuery);
		glBeginTransformFeedback(GL_POINTS);
		glDrawArrays(GL_POINTS, 0, _vertexNumber);	// _nNumVerts
		glEndTransformFeedback();
		glEndQuery(GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN);
		// Query the number of primitives written in the transform buffer.
		GLuint PrimitivesWritten;
		glGetQueryObjectuiv(_transformFeedbackQuery, GL_QUERY_RESULT, &PrimitivesWritten);
		glDisable(GL_RASTERIZER_DISCARD);
	/*	glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[1]);
		std::vector<GLfloat> fbCoords;
		fbCoords.assign(_vertexNumber*3,0.0f);
		glGetBufferSubData(GL_ARRAY_BUFFER,0,sizeof(GLfloat)*_vertexNumber*3,&(fbCoords[0]));
		glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[0]);
		fbCoords.assign(_vertexNumber*4,0.0f);
		glGetBufferSubData(GL_ARRAY_BUFFER,0,sizeof(GLfloat)*_vertexNumber*4,&(fbCoords[0])); */
	}
}

void glslTriangle::computeLocalBounds()
{
	float minMaxXYZ[6];
	minMaxXYZ[0]=minMaxXYZ[2]=minMaxXYZ[4]=1e30f;
	minMaxXYZ[1]=minMaxXYZ[3]=minMaxXYZ[5]=-1e30f;
	int i,n=(unsigned int)_Positions.size()>>2;
	GLfloat *v=&(_Positions[0]);
	for(i=0; i<n; ++i)	{
		if(*v<minMaxXYZ[0])
			minMaxXYZ[0]=*v;
		if(*v>minMaxXYZ[1])
			minMaxXYZ[1]=*v;
		++v;
		if(*v<minMaxXYZ[2])
			minMaxXYZ[2]=*v;
		if(*v>minMaxXYZ[3])
			minMaxXYZ[3]=*v;
		++v;
		if(*v<minMaxXYZ[4])
			minMaxXYZ[4]=*v;
		if(*v>minMaxXYZ[5])
			minMaxXYZ[5]=*v;
		++v; ++v;
	}
	_localCenter[0] = (minMaxXYZ[0]+minMaxXYZ[1])*0.5f;
	_localCenter[1] = (minMaxXYZ[2]+minMaxXYZ[3])*0.5f;
	_localCenter[2] = (minMaxXYZ[4]+minMaxXYZ[5])*0.5f;
	_radius = (_localCenter[0]-minMaxXYZ[0])*(_localCenter[0]-minMaxXYZ[0]);
	_radius += (_localCenter[1]-minMaxXYZ[2])*(_localCenter[1]-minMaxXYZ[2]);
	_radius += (_localCenter[2]-minMaxXYZ[4])*(_localCenter[2]-minMaxXYZ[4]);
	_radius = sqrt(_radius);
	_boundsComputed=true;
}

//void glslTriangle::getLocalBounds(GLfloat (&localCenter)[3], GLfloat &Radius) {
//	localCenter[0]=_localCenter[0]; localCenter[1]=_localCenter[1]; localCenter[2]=_localCenter[2];
//	Radius = _radius;
//}

bool glslTriangle::parseNextInputFileLine(std::ifstream *infile, std::string &unparsedLine, std::vector<std::string> &parsedLine)
{
	if(infile->eof())
	{
		infile->close();
		return false;
	}
	if(!infile->is_open())
		return false;
	char s[400];
	infile->getline(s,399);
	unparsedLine.assign(s);
	if(unparsedLine=="")
	{
		parsedLine.clear();
		return true;
	}
	parsedLine.clear();
	std::string substr;
	int start=0,next=0;
	while(next<399 && s[next]!= '\0')
	{
		// advance to whitespace
		while(next<399 && s[next]!='\0' && s[next]!=' ' && s[next]!='\t')
			++next;
		substr.assign(s+start,s+next);
		parsedLine.push_back(substr);
		while(s[next]==' ' || s[next]=='\t')
			++next;
		start = next;
	}
	return true;
}

bool glslTriangle::writeObjFile(const char *fileName)
{
	std::string title(fileName);
	if(title.rfind(".obj")>title.size())
		title.append(".obj");
	std::ofstream fin(title.c_str());
    if(!fin.is_open())
        return false;
	if(!_adjacenciesComputed)
		findAdjacentTriangles();
	std::vector<float> normals;
	normals.assign((_uniqueNormalNumber+1)*3,-1000.0f);
	int numNorms=0;
	for(size_t j,i=0; i<_pnI.size(); ++i)	{
		j = _pnI[i].posNormIdx;
		if(j>numNorms)
			numNorms=(int)j;
		assert(j<=_uniqueNormalNumber);
		if(normals[j*3]<-100.0f)
			getVertexNormal((unsigned int)i,(float (&)[3])(normals[j*3]));	}
	char s[400];
	std::string line;
	for(size_t i=0; i<_Positions.size(); i+=4)	{
		sprintf(s,"v %f %f %f\n",_Positions[i],_Positions[i+1],_Positions[i+2]);
		line.assign(s);
		fin.write(line.c_str(),line.size());	}
	for(size_t i=0; i<_TexCoords.size(); i+=2)	{
		sprintf(s,"vt %f %f\n",_TexCoords[i],_TexCoords[i+1]);
		line.assign(s);
		fin.write(line.c_str(),line.size());	}
	for(size_t i=0; i<(numNorms+1)*3; i+=3)	{
		sprintf(s,"vn %f %f %f\n",normals[i],normals[i+1],normals[i+2]);
		line.assign(s);
		fin.write(line.c_str(),line.size());	}
	for(size_t i=0; i<_tris.size(); i+=3)	{
		sprintf(s,"f %d/%d/%d %d/%d/%d %d/%d/%d\n",(_pnI[_tris[i]].posIdx>>2)+1,_tris[i]+1,_pnI[_tris[i]].posNormIdx+1,
			(_pnI[_tris[i+1]].posIdx>>2)+1,_tris[i+1]+1,_pnI[_tris[i+1]].posNormIdx+1,(_pnI[_tris[i+2]].posIdx>>2)+1,_tris[i+2]+1,_pnI[_tris[i+2]].posNormIdx+1);
		line.assign(s);
		fin.write(line.c_str(),line.size());	}
	fin.close();
	return true;
}

int glslTriangle::readObjFile(const char *fileName, bool dataOnlyNoGraphics)
{ // returned error codes: 0=no error, 1=can't open file, 2=non-triangle primitive,
	// 3=exceeded 0x7fff vertex position limit, 4=bad 3D vertex line, 5=bad 2D texture line, 6=nonmanifold surface
	if(!_static && _NormalMakerProgram<1 && !dataOnlyNoGraphics)	// do once for the class
		createNormalMakerProgram();
    std::ifstream fin(fileName);
    if(!fin.is_open())
        return 1;
    std::string unparsedLine;
	std::vector<std::string> parsedLine;
	_vertexNumber=0;
	_triangleNumber=0;
	_uniqueNormalNumber=0;
	cnt cntin;
	CNTMAP cntmap;
	cntmap.clear();
	std::pair<CNTMAP::iterator,bool> cpr;
	std::vector<float> tt,tp;
	tt.clear(); tp.clear();
	int i,j,k,l;
	std::string str;
	char s[100];
	std::vector<GLuint> tempTris;
	while(parseNextInputFileLine(&fin,unparsedLine,parsedLine))
	{
		if(parsedLine.empty())
			continue;
		if(parsedLine[0]=="v")
		{
			if(parsedLine.size()!=4)
				return 4;
			for(i=1; i<4; ++i)
				tp.push_back((float)atof(parsedLine[i].c_str()));
			tp.push_back(1.0f);
		}
		else if(parsedLine[0]=="vt")
		{
			if(parsedLine.size()!=3)
				return 5;
			tt.push_back((float)atof(parsedLine[1].c_str()));
			tt.push_back((float)atof(parsedLine[2].c_str()));
		}
		else if(parsedLine[0]=="vn")
			++_uniqueNormalNumber;
		else if(parsedLine[0]=="f")
		{	// always in vertexPosition/vertexTexture/vertexNormal format. If vP/vT may skip normal. If vP//vN, texture is skipped.
			int numVerts = (int)parsedLine.size();
			if(numVerts<4 || numVerts>5)
				return 2;
			for(i=1; i<numVerts; ++i)
			{
				strncpy(s,parsedLine[i].c_str(),99);
				k=0; l=0;
				for(j=0; j<3; ++j)
				{
					while(s[l]!='/' && s[l]!='\0')
						++l;
					str.clear();
					while(k!=l)
					{
						str.push_back(s[k]);
						++k;
					}
					if(str.empty())
						cntin.v[j] = -1;
					else
						cntin.v[j] = atoi(str.c_str()) - 1;	// remember indexes in obj files start at 1
					++k; ++l;
				}
				cpr = cntmap.insert(std::make_pair(cntin,_vertexNumber));
				if(cpr.second)
					++_vertexNumber;
				tempTris.push_back(cpr.first->second);
			}
			++_triangleNumber;
			if(numVerts>4)	{
				int tts = (int)tempTris.size();
				tts -= 4;
				tempTris.push_back(tempTris[tts]);
				tempTris.push_back(tempTris[tts+2]);
				++_triangleNumber;
			}
		}
		else
			continue;
	}
	fin.close();
	if(	_vertexNumber>0x7fff) {
		return 3;
	}
	_tris.clear();
	_tris.assign(tempTris.begin(),tempTris.end());
	tempTris.clear();
	_Positions.clear();
	_Positions.assign(tp.begin(),tp.end());
	tp.clear();
	_pnI.clear();
	_pnI.assign(cntmap.size(),positionNormalIndex());
	if(_uniqueNormalNumber<1)
		_uniqueNormalNumber = (unsigned long)cntmap.size();
	_TexCoords.clear();
	if(_textured)
		_TexCoords.assign(cntmap.size()<<1,-1.0f);
	GLfloat *tout,*tin;
	CNTMAP::iterator cit=cntmap.begin();
	std::map<edge,int,edgeTest> pnMap;
	std::pair<std::map<edge,int,edgeTest>::iterator,bool> pr;
	int pnCount=0;
	while(cit!=cntmap.end()) {
		_pnI[cit->second].posIdx = (long)cit->first.v[0]<<2;
		if(cit->first.v[2]<0)	{	// no normal specified means if new texture is a hard normal.
			edge edg;
			edg.vtxMin=cit->first.v[0];
			edg.vtxMax=cit->first.v[1];
			pr = pnMap.insert(std::make_pair(edg,pnCount));
			if(pr.second)
				++pnCount;
			_pnI[cit->second].posNormIdx = pr.first->second;
		}
		else
			_pnI[cit->second].posNormIdx = cit->first.v[2];
		if(_textured)	{
			tin = &(tt[(cit->first.v[1])<<1]);
			tout = &(_TexCoords[(cit->second)<<1]);
			*tout = *tin;
			++tout; ++tin;
			*tout = *tin;	}
		++cit;
	}
	cntmap.clear();
	if(!findAdjacentTriangles()) {
		_vertexNumber=0;
		_triangleNumber=0;
		_uniqueNormalNumber=0;
		_tris.clear();
		_Positions.clear();
		_pnI.clear();
		_TexCoords.clear();
		return 6;
	}
	if(!dataOnlyNoGraphics)
		sendToGraphicsCard();
	return 0;
}

bool glslTriangle::findAdjacentTriangles(bool uniqueNormalsNotCoords)
{	// computes all the adjacent triangles from raw triangle input
	// create binary tree to store and find edges efficiently
	// returns false if non-manifold surface is input
	typedef std::set<edge,edgeTest> edgeSet;
	typedef edgeSet::iterator edgeIt;
	std::pair <edgeIt,bool> P;
	edge E;
	edgeIt ei;
	edgeSet M;
	M.clear();
	unsigned long tnow[3],*adjNow;
	unsigned long i,j,tcode,numtris=(unsigned long)_tris.size();
	if(uniqueNormalsNotCoords) {
		_adjs.clear();
		_adjs.assign(numtris,0x00000003);
	}
	else {
		_adjsP.clear();
		_adjsP.assign(numtris,0x00000003);
	}
	for(i=0; i<numtris; i+=3)
	{
		if(_tris[i] == 0xffffffff)	// signals a deleted triangle
			continue;
		for(j=0; j<3; j++) {
			tnow[j] = _tris[i+j];
			if(uniqueNormalsNotCoords)
				tnow[j] = _pnI[tnow[j]].posNormIdx;
			else
				tnow[j] = _pnI[tnow[j]].posIdx>>2;
		}
		if(uniqueNormalsNotCoords)
			adjNow = &(_adjs[i]);
		else
			adjNow = &(_adjsP[i]);
		for(j=0; j<3; j++)
		{
			if(adjNow[j]!=0x00000003)	// adjacency already computed
				continue;
			unsigned int tmp;
			if((tmp=tnow[(j+1)%3])<tnow[j]) {
				E.vtxMin = tmp;
				E.vtxMax = tnow[j];
				E.reversed = 1;	}
			else {
				E.vtxMin = tnow[j];
				E.vtxMax = tmp;
				E.reversed = 0;	}
			E.adjCode = ((i/3)<<2) + j;
			P = M.insert(E);
			// if P.second is true, no match so edge inserted
			if(P.second==false)	// edge match found
			{
				if(P.first->reversed==E.reversed && E.vtxMin!=E.vtxMax)
					return false;
				tcode = P.first->adjCode;
				adjNow[j] = tcode;
				if(uniqueNormalsNotCoords)
					_adjs[(tcode>>2)*3+(tcode&0x00000003)] = E.adjCode;
				else
					_adjsP[(tcode>>2)*3+(tcode&0x00000003)] = E.adjCode;
				M.erase(P.first);
			}
			else
				adjNow[j] = 0x00000003;
		}
	}
//	if(!uniqueNormalsNotCoords)	{	// to check veracity of incision tool
//		checkRing(M);
//		checkRing(M);	}
	M.clear();
	makeVertexToTriangleMap(uniqueNormalsNotCoords);
	if(uniqueNormalsNotCoords)	_adjacenciesComputed = 1;
	else	_adjacenciesComputedP = 1;
	return true;
}

void glslTriangle::checkRing(std::set<edge,edgeTest> &M)
{
	if(M.empty())
		return;
	std::list<int> l1;
	std::set<edge,edgeTest>::const_iterator mit=M.begin();
	l1.push_back(mit->vtxMin);
	l1.push_back(mit->vtxMax);
	M.erase(mit);
	while(l1.front()!=l1.back())	{
		mit=M.begin();
		while( mit!=M.end())	{
			if(mit->vtxMin==l1.front())
				l1.push_front(mit->vtxMax);
			else if(mit->vtxMax==l1.front())
				l1.push_front(mit->vtxMin);
			else if(mit->vtxMin==l1.back())
				l1.push_back(mit->vtxMax);
			else if(mit->vtxMax==l1.back())
				l1.push_back(mit->vtxMin);
			else	{
				++mit;
				continue;	}
			M.erase(mit);
			break;
		}
	}
	return;
}

int glslTriangle::getVertexTriangle(int vertexNumber)	// gets triangle vertex is a member of
{
	if(_adjacenciesComputed)
		return _vertexFace[_pnI[vertexNumber].posNormIdx]&0x3fffffff;
	size_t i,n=_tris.size();
	GLuint glv = vertexNumber;
	for(i=0; i<n; ++i) {
		if(_tris[i]==glv)
			break;
	}
	if(i==n)
		return -1;
	return (int)(i/3);
}

void glslTriangle::makeVertexToTriangleMap(bool uniqueNormalsNotCoords)
{
	int i,j,numtris=(int)_tris.size();
	unsigned int *tnow,vnow;
	// provide each vertex with a face it is a member of
	if(uniqueNormalsNotCoords)	{
		_vertexFace.clear();
		_vertexFace.assign(_vertexNumber,0x80000000);	// initially deleted
	}
	else {
		_vertexFaceP.clear();
		_vertexFaceP.assign(_Positions.size()>>2,0x80000000);	// initially deleted
	}
	for(i=0; i<numtris; i+=3)
	{
		tnow = &(_tris[i]);
		if(tnow==NULL || tnow[0] == 0xffffffff)	// signals a deleted triangle
			continue;
		for(j=0; j<3; j++)
		{
			vnow = tnow[j];
			if(uniqueNormalsNotCoords)	{
				vnow = _pnI[vnow].posNormIdx;
				if(_vertexFace[vnow]&0x40000000)
					continue;	// vertex first on free edge, don't change
				_vertexFace[vnow] = i/3;
				if(_adjs[i+j]==0x00000003)	// vertex first on free edge, lock it for easy neighbor find
					_vertexFace[vnow] |= 0x40000000;
			}
			else	{
				vnow = _pnI[vnow].posIdx>>2;
				if(_vertexFaceP[vnow]&0x40000000)
					continue;	// vertex first on free edge, don't change
				_vertexFaceP[vnow] = i/3;
				if(_adjsP[i+j]==0x00000003)	// vertex first on free edge, lock it for easy neighbor find
					_vertexFaceP[vnow] |= 0x40000000;
			}
		}
	}
}

void glslTriangle::getNeighbors(unsigned int vertex, std::vector<neighborNode> &neighbors, bool uniqueNormalsNotCoords)
{
	unsigned long normVert,trNum,adj,triStart;
	if(uniqueNormalsNotCoords) { // make sure it is a unique normal-vertex
		if(!_adjacenciesComputed)
			findAdjacentTriangles(uniqueNormalsNotCoords);
		vertex = _pnI[vertex].posNormIdx;
		triStart = _vertexFace[vertex];
	}
	else {
		if(!_adjacenciesComputedP)
			findAdjacentTriangles(uniqueNormalsNotCoords);
		vertex = _pnI[vertex].posIdx>>2;
		triStart = _vertexFaceP[vertex];
	}
	neighbors.clear();
	if(triStart&0x80000000)	// unconnected vertex
		return;
	trNum = triStart&0x3fffffff;
	unsigned long *adjs;
	GLuint *tnow = &(_tris[(trNum<<1)+trNum]);
	assert(tnow[0]!=0xffffffff);	// deleted triangle
	int j;
	for(j=0; j<3; ++j)
	{
		normVert = tnow[j];
		if(uniqueNormalsNotCoords)
			normVert = _pnI[normVert].posNormIdx;
		else
			normVert = _pnI[normVert].posIdx>>2;
		if(normVert==vertex)
			break;
	}
	assert(j<3);
	if(uniqueNormalsNotCoords)
		adjs = &(_adjs[(trNum<<1)+trNum]);
	else
		adjs = &(_adjsP[(trNum<<1)+trNum]);
	// set triStart to the end adjacency code for counterclockwise traversal
	neighborNode n;
	if(triStart&0x40000000)	// started on a free edge, will end on one
	{
		triStart = 0x00000003;
		n.triangle = -1;	// code for open ring neighbors
		n.vertex = tnow[(j+1)%3];
		neighbors.push_back(n);
	}
	else	// create adjacency code of starting edge
		triStart = (trNum<<2)+j;
	n.triangle = trNum;
	n.vertex = tnow[(j+2)%3];
	neighbors.push_back(n);
	adj = adjs[(j+2)%3];
	while(adj!=triStart)
	{
		n.triangle = adj>>2;
		tnow = triVerts(n.triangle);
		if(uniqueNormalsNotCoords)
			adjs = &(_adjs[(n.triangle<<1)+n.triangle]);
		else
			adjs = &(_adjsP[(n.triangle<<1)+n.triangle]);
		j = adj&0x00000003;
		n.vertex = tnow[(j+2)%3];
		neighbors.push_back(n);
		adj = adjs[(j+2)%3];
	}
}

bool glslTriangle::sendTopologyChangeToGraphicsCard()
{	// sends current triangle object to glsl capable graphics card
    // Copy data to video memory
	if(!_static)	{
		// Vertex data
		glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, _bufferObjects[0]);	// VERTEX_DATA
		glBufferData(GL_TRANSFORM_FEEDBACK_BUFFER, sizeof(GLfloat)*_pnI.size()*4, NULL, GL_DYNAMIC_COPY);
		// Normal data
		glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, _bufferObjects[1]);	// NORMAL_DATA
		glBufferData(GL_TRANSFORM_FEEDBACK_BUFFER, sizeof(GLfloat)*_pnI.size()*3, NULL, GL_DYNAMIC_COPY);
	}
	else	{
		// Vertex data
		glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[0]);	// VERTEX_DATA
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_pnI.size()*4, NULL, GL_DYNAMIC_DRAW);
		// Normal data
		glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[1]);	// NORMAL_DATA
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_pnI.size()*3, NULL, GL_DYNAMIC_DRAW);
	}
    // Indexes
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bufferObjects[2]);	// INDEX_DATA
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*_tris.size(), &(_tris[0]), GL_STATIC_DRAW);
	if(_textured)	{
	    // Texture coordinates
		glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[3]);	// TEXTURE_DATA
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_pnI.size()*2, &(_TexCoords[0]), GL_STATIC_DRAW);	}
	// release for next use
	glBindBuffer( GL_ARRAY_BUFFER, 0);
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0);
	if(!_static)	{
		// transform feedback buffers
		glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER,0,_bufferObjects[1]);
		glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER,1,_bufferObjects[0]);
		// Load the _Positions first
		glBindBuffer(GL_TEXTURE_BUFFER_ARB, _texBOBuffers[0]);
		glBufferData(GL_TEXTURE_BUFFER_ARB, sizeof(GLfloat)*_Positions.size(), &(_Positions[0]), GL_DYNAMIC_DRAW);
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_BUFFER_ARB, _textureBufferObjects[0]);
		glTexBufferARB(GL_TEXTURE_BUFFER_ARB, GL_RGBA32F, _texBOBuffers[0]); 
		// Load the neighbor info second
		std::vector<GLuint> topoArr;
		createGlslTopoArray(topoArr);
		glBindBuffer(GL_TEXTURE_BUFFER_ARB, _texBOBuffers[1]);
		glBufferData(GL_TEXTURE_BUFFER_ARB, sizeof(GLuint)*topoArr.size(), &(topoArr[0]), GL_STATIC_DRAW);
		topoArr.clear();
		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_BUFFER_ARB, _textureBufferObjects[1]);
		glTexBufferARB(GL_TEXTURE_BUFFER_ARB, GL_RGBA32UI, _texBOBuffers[1]); 
		glActiveTexture(GL_TEXTURE0);
		glBindBuffer( GL_TEXTURE_BUFFER_ARB, 0);
	}
	if(_NormalMakerProgram)	{
		assert(!_static);
		createNormalMakerProgram();	}
	return true;
}

bool glslTriangle::sendToGraphicsCard()
{
	// Create the buffer objects once
	if(!_vertexArrayBufferObject)
		glGenVertexArrays(1,&_vertexArrayBufferObject);
//	else	{
//		glDeleteVertexArrays(1,&_vertexArrayBufferObject);
//		glGenVertexArrays(1,&_vertexArrayBufferObject);	}
//	if(!_vertexArrayBufferObject)
//		glGenVertexArrays(1,&_vertexArrayBufferObject);
	if(_textured)	{
		if(!_bufferObjects[0])
			glGenBuffers(4, _bufferObjects);
		else	{
			glDeleteBuffers(4, _bufferObjects);
			glGenBuffers(4, _bufferObjects);	}
/*		if(!_static)	{
			// Create 2 new texture buffer objects
			if(!_textureBufferObjects[0])
				glGenTextures(2,_textureBufferObjects);
			else	{
				glDeleteTextures(2,_textureBufferObjects);
				glGenTextures(2,_textureBufferObjects);	}
		} */
//		if(!_bufferObjects[0])
//			glGenBuffers(4, _bufferObjects);
	}
	else	{
		if(!_bufferObjects[0])
			glGenBuffers(3, _bufferObjects);	}
/*	if(!_static)	{
		// Create 2 new texture buffer objects
		if(!_textureBufferObjects[0])
			glGenTextures(2,_textureBufferObjects);
		if(!_texBOBuffers[0])
			glGenBuffers(2,_texBOBuffers);	} */
	if(!_static)	{
		// Create 2 new texture buffer objects
		if(!_textureBufferObjects[0])
			glGenTextures(2,_textureBufferObjects);
		else	{
			glDeleteTextures(2,_textureBufferObjects);
			glGenTextures(2,_textureBufferObjects);	}
		if(!_texBOBuffers[0])
			glGenBuffers(2,_texBOBuffers);
		else	{
			glDeleteBuffers(2,_texBOBuffers);
			glGenBuffers(2,_texBOBuffers);
		}
	}

    // Copy data to video memory
	if(!_static)	{
		// Vertex data
		glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, _bufferObjects[0]);	// VERTEX_DATA
		glBufferData(GL_TRANSFORM_FEEDBACK_BUFFER, sizeof(GLfloat)*_pnI.size()*4, NULL, GL_DYNAMIC_COPY);
		// Normal data
		glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, _bufferObjects[1]);	// VERTEX_DATA
		glBufferData(GL_TRANSFORM_FEEDBACK_BUFFER, sizeof(GLfloat)*_pnI.size()*3, NULL, GL_DYNAMIC_COPY);
	}
	else	{
		// Vertex data
		glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[0]);	// VERTEX_DATA
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_pnI.size()*4, NULL, GL_DYNAMIC_DRAW);
		// Normal data
		glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[1]);	// NORMAL_DATA
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_pnI.size()*3, NULL, GL_DYNAMIC_DRAW);
	}
    // Indexes
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bufferObjects[2]);	// INDEX_DATA
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*_tris.size(), &(_tris[0]), GL_STATIC_DRAW);
	if(_textured)	{
	    // Texture coordinates
		glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[3]);	// TEXTURE_DATA
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_pnI.size()*2, &(_TexCoords[0]), GL_STATIC_DRAW);	}
	// release for next use
	glBindBuffer( GL_ARRAY_BUFFER, 0);
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0);

	glBindVertexArray(_vertexArrayBufferObject);
    // Vertex data
    glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[0]);	// VERTEX_DATA
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
//	glEnableClientState( GL_VERTEX_ARRAY);
    // Normal data
    glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[1]);	// NORMAL_DATA
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
    // Indexes
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bufferObjects[2]);	// INDEX_DATA
	if(_textured)	{
	    // Texture coordinates
		glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[3]);	// TEXTURE_DATA
		glEnableVertexAttribArray(2);
		glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, 0);
	}
	// don't release for next use GL_ARRAY_BUFFER or GL_ELEMENT_ARRAY_BUFFER inside a VBO. Unbinds everything.
    // Unbind to anybody
	glBindVertexArray(0);

	if(!_static)	{
		// transform feedback buffers
		glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER,0,_bufferObjects[1]);
	//    glBufferData(GL_TRANSFORM_FEEDBACK_BUFFER, sizeof(GLfloat)*_nNumVerts*3, NULL, GL_DYNAMIC_COPY);
		glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER,1,_bufferObjects[0]);
	//    glBufferData(GL_TRANSFORM_FEEDBACK_BUFFER, sizeof(GLfloat)*_nNumVerts*4, NULL, GL_DYNAMIC_COPY);
	//	glBindVertexArray(0);

		// Load the _Positions first
		glBindBuffer(GL_TEXTURE_BUFFER_ARB, _texBOBuffers[0]);
		glBufferData(GL_TEXTURE_BUFFER_ARB, sizeof(GLfloat)*_Positions.size(), &(_Positions[0]), GL_DYNAMIC_DRAW);
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_BUFFER_ARB, _textureBufferObjects[0]);
		glTexBufferARB(GL_TEXTURE_BUFFER_ARB, GL_RGBA32F, _texBOBuffers[0]); 
		// Load the neighbor info second
		std::vector<GLuint> topoArr;
		createGlslTopoArray(topoArr);
		glBindBuffer(GL_TEXTURE_BUFFER_ARB, _texBOBuffers[1]);
		glBufferData(GL_TEXTURE_BUFFER_ARB, sizeof(GLuint)*topoArr.size(), &(topoArr[0]), GL_STATIC_DRAW);
		topoArr.clear();
	//	glBindBuffer(GL_TEXTURE_BUFFER_ARB, 0);
		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_BUFFER_ARB, _textureBufferObjects[1]);
		glTexBufferARB(GL_TEXTURE_BUFFER_ARB, GL_RGBA32UI, _texBOBuffers[1]); 
		glActiveTexture(GL_TEXTURE0);
	//	glBindBuffer( GL_TEXTURE_BUFFER_ARB, 0);
	}
	return true;
}

void glslTriangle::createGlslTopoArray(std::vector<GLuint> &topoArr)
{
	topoArr.clear();
	std::vector<neighborNode> nei;
	int i,j,k,m,n=(int)_pnI.size();
	topoArr.assign(n<<2,0);
	for(i=0; i<n; ++i)	{
		topoArr[i<<2] = _pnI[i].posIdx>>2;
		getNeighbors(i,nei);
		topoArr[(i<<2)+1] = _pnI[nei.back().vertex].posIdx>>2;
		if(nei[0].triangle<0)
			topoArr[(i<<2)+1] |= 0x80000000;
		m = (int)nei.size();
		for(j=0; j<m-1; ++j) { // order needs to be reversed since vertices listed clockwise
			k = m-j-2;
			if(j<2)
				topoArr[(i<<2)+2+j] = _pnI[nei[k].vertex].posIdx>>2;
			else if(j==2) {
				topoArr.push_back(topoArr[(i<<2)+3]);
				topoArr[(i<<2)+3] = (unsigned long)topoArr.size()>>2;
				topoArr.push_back(_pnI[nei[k].vertex].posIdx>>2);
			}
			else
				topoArr.push_back(_pnI[nei[k].vertex].posIdx>>2);
			if(j+2==m) {
				if(j<2)
					topoArr[(i<<2)+2+j] |= 0x80000000;
				else
					topoArr.back() |= 0x80000000;
			}
		}
		if(j>2) { // pushed onto top of texel array. Pad to even texel
			 m = 4 - ((j-1)&0x0003);
			 for(j=0; j<m; ++j)
				 topoArr.push_back(0);
		}
	}
}

void glslTriangle::getVertexNormal(unsigned int vertex, float (&normal)[3], bool uniqueNormalsNotCoords)
{
	normal[0]=0.0f; normal[1]=0.0f; normal[2]=0.0f;
	if(uniqueNormalsNotCoords && !_adjacenciesComputed)
		findAdjacentTriangles(uniqueNormalsNotCoords);
	if(!uniqueNormalsNotCoords && !_adjacenciesComputedP)
		findAdjacentTriangles(uniqueNormalsNotCoords);
	GLfloat *v,*nv,last[3],now[3];
	v = vertexCoordinate(vertex);
	std::vector<neighborNode> nei;
	getNeighbors(vertex,nei,uniqueNormalsNotCoords);
	int i,n=(int)nei.size();
	if(n<2)	return;
	if(nei[0].triangle>-1) {
		nv = vertexCoordinate(nei.back().vertex);
		i=0; }
	else {
		nv = vertexCoordinate(nei.front().vertex);
		i=1; }
	last[0]=nv[0]-v[0]; last[1]=nv[1]-v[1]; last[2]=nv[2]-v[2];
	while(i<n) {
		nv = vertexCoordinate(nei[i].vertex);
		now[0]=nv[0]-v[0]; now[1]=nv[1]-v[1]; now[2]=nv[2]-v[2];
		normal[0] += last[1]*now[2] - last[2]*now[1];
		normal[1] += last[2]*now[0] - last[0]*now[2];
		normal[2] += last[0]*now[1] - last[1]*now[0];
		last[0]=now[0]; last[1]=now[1]; last[2]=now[2];
		++i;
	}
	float len=sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
	if(len>0.0f) {
		normal[0]/=len; normal[1]/=len; normal[2]/=len; }
}

bool glslTriangle::localPick(const float *lineStart, float *lineDirection, float (&position)[3], int &triangle, float &param)
{ // lineStart and lineDirection are both 3 element vectors
	bool picked=false;
	int j,bigAxis=0;
	if(fabs(lineDirection[1])>fabs(lineDirection[0]))
		bigAxis =1;
	if(fabs(lineDirection[2])>fabs(lineDirection[bigAxis]))
		bigAxis =2;
	unsigned int *tr,i;
	GLfloat *v[3];
	float t,vtx[3],minimax[6],tMin,tMax,rMin,rMax;
	for(i=0; i<_triangleNumber; ++i)	{
		tr = triVerts(i);
		v[0]=vertexCoordinate(tr[0]);
		v[1]=vertexCoordinate(tr[1]);
		v[2]=vertexCoordinate(tr[2]);
		minimax[0]=minimax[1]=v[0][0];
		minimax[2]=minimax[3]=v[0][1];
		minimax[4]=minimax[5]=v[0][2];
		for(j=1; j<3; ++j) {
			if(minimax[0]>v[j][0])
				minimax[0]=v[j][0];
			if(minimax[1]<v[j][0])
				minimax[1]=v[j][0];
			if(minimax[2]>v[j][1])
				minimax[2]=v[j][1];
			if(minimax[3]<v[j][1])
				minimax[3]=v[j][1];
			if(minimax[4]>v[j][2])
				minimax[4]=v[j][2];
			if(minimax[5]<v[j][2])
				minimax[5]=v[j][2];
		}
		tMin = (minimax[bigAxis<<1]-lineStart[bigAxis])/lineDirection[bigAxis];
		tMax = (minimax[(bigAxis<<1)+1]-lineStart[bigAxis])/lineDirection[bigAxis];
		for(j=0; j<3; ++j) {
			if(j==bigAxis)
				continue;
			rMin=lineStart[j]+lineDirection[j]*tMin;
			rMax=lineStart[j]+lineDirection[j]*tMax;
			if(rMin<minimax[j<<1] && rMax<minimax[j<<1])
				break;
			if(rMin>minimax[(j<<1)+1] && rMax>minimax[(j<<1)+1])
				break;
		}
		if(j<3)
			continue;
		if(m3dRayTriangleIntersection(lineStart,lineDirection,v[0],v[1],v[2],t,vtx)) {
			if(t<param){
				param=t;
				picked=true;
				triangle = i;
				position[0]=vtx[0]; position[1]=vtx[1]; position[2]=vtx[2];
			}
		}
	}
	return picked;
}

int glslTriangle::linePick(const float *lineStart, float *lineDirection, std::vector<float> &positions, std::vector<int> &triangles, std::vector<float> &params)
{ // lineStart and lineDirection are both 3 element vectors
	posTri pT;
	std::map<float,posTri> hitMap;
	std::map<float,posTri>::iterator hit;
	int j,bigAxis=0;
	if(fabs(lineDirection[1])>fabs(lineDirection[0]))
		bigAxis =1;
	if(fabs(lineDirection[2])>fabs(lineDirection[bigAxis]))
		bigAxis =2;
	unsigned int *tr,i;
	GLfloat *v[3];
	float t,vtx[3],minimax[6],tMin,tMax,rMin,rMax;
	for(i=0; i<_triangleNumber; ++i)	{
		tr = triVerts(i);
		v[0]=vertexCoordinate(tr[0]);
		v[1]=vertexCoordinate(tr[1]);
		v[2]=vertexCoordinate(tr[2]);
		minimax[0]=minimax[1]=v[0][0];
		minimax[2]=minimax[3]=v[0][1];
		minimax[4]=minimax[5]=v[0][2];
		for(j=1; j<3; ++j) {
			if(minimax[0]>v[j][0])
				minimax[0]=v[j][0];
			if(minimax[1]<v[j][0])
				minimax[1]=v[j][0];
			if(minimax[2]>v[j][1])
				minimax[2]=v[j][1];
			if(minimax[3]<v[j][1])
				minimax[3]=v[j][1];
			if(minimax[4]>v[j][2])
				minimax[4]=v[j][2];
			if(minimax[5]<v[j][2])
				minimax[5]=v[j][2];
		}
		tMin = (minimax[bigAxis<<1]-lineStart[bigAxis])/lineDirection[bigAxis];
		tMax = (minimax[(bigAxis<<1)+1]-lineStart[bigAxis])/lineDirection[bigAxis];
		for(j=0; j<3; ++j) {
			if(j==bigAxis)
				continue;
			rMin=lineStart[j]+lineDirection[j]*tMin;
			rMax=lineStart[j]+lineDirection[j]*tMax;
			if(rMin<minimax[j<<1] && rMax<minimax[j<<1])
				break;
			if(rMin>minimax[(j<<1)+1] && rMax>minimax[(j<<1)+1])
				break;
		}
		if(j<3)
			continue;
		if(m3dRayTriangleIntersection(lineStart,lineDirection,v[0],v[1],v[2],t,vtx)) {
			pT.triangle = i;
			pT.v[0]=vtx[0]; pT.v[1]=vtx[1]; pT.v[2]=vtx[2];
			hitMap.insert(std::make_pair(t,pT));
		}
	}
	positions.clear();	triangles.clear();	params.clear();
	positions.reserve(hitMap.size()*3);
	triangles.reserve(hitMap.size());
	params.reserve(hitMap.size());
	for(hit=hitMap.begin(); hit!=hitMap.end(); ++hit)	{
		params.push_back(hit->first);
		triangles.push_back(hit->second.triangle);
		positions.push_back(hit->second.v[0]);
		positions.push_back(hit->second.v[1]);
		positions.push_back(hit->second.v[2]);
	}
	return (int)hitMap.size();
}

void glslTriangle::getTriangulatedSurface(std::vector<float> &vertices, std::vector<int> &tris)
{
	vertices.clear();
	tris.clear();
	size_t n=_Positions.size()>>2;
	vertices.reserve(n*3);
//	vertices.assign(n*3,0.0f);
	for(int i=0; i<n; ++i)	{
		vertices.push_back(_Positions[i<<2]);
		vertices.push_back(_Positions[(i<<2)+1]);
		vertices.push_back(_Positions[(i<<2)+2]);
//		vertices[i*3]=_Positions[i<<2];
//		vertices[i*3+1]=_Positions[(i<<2)+1];
//		vertices[i*3+2]=_Positions[(i<<2)+2];
	}
	n = _tris.size();
	tris.assign(n,0);
	for(int i=0; i<n; ++i)
		tris[i]=(_pnI[_tris[i]].posIdx)>>2;
}

void glslTriangle::getSurfaceTriangles(std::vector<int> &triangles)
{
    triangles.clear();
    triangles.resize( _tris.size() );
    for( int i = 0; i < _tris.size(); i++){
        triangles[i] = _tris[i];       
    }
}

void glslTriangle::getSurfaceVertices(std::vector<float> &vertices)
{
    vertices.clear();
    vertices.resize( _pnI.size()*3 );
    for( int i = 0; i < _pnI.size(); i++){
        vertices[i*3+0] = _Positions[_pnI[i].posIdx+0];
        vertices[i*3+1] = _Positions[_pnI[i].posIdx+1];
        vertices[i*3+2] = _Positions[_pnI[i].posIdx+2];       
    }
}

void glslTriangle::getSurfaceNormals(std::vector<float> &normals)
{
    //assert( _uniqueNormalNumber == _pnI.size() );
    normals.clear();
    normals.resize(_pnI.size()*3 );
	for(size_t i=0; i<_pnI.size(); ++i)
        getVertexNormal((unsigned int)i,(float (&)[3])(normals[i*3]));
}

void glslTriangle::getSurfaceUVs(std::vector<float> &uv)
{
    uv.clear();
    uv.resize( _pnI.size()*2 );
    for( int i = 0; i < _tris.size(); i++){
        assert( _tris[i] < _TexCoords.size() );
        uv[(_tris[i]*2)+0] = _TexCoords[(_tris[i]*2)+0];
        uv[(_tris[i]*2)+1] = _TexCoords[(_tris[i]*2)+1];
    }
}


void glslTriangle::getClosestBarycentric(int triangle, float (&xyz)[3], float (&uv)[2])
{	// for position xyz return barycentric uv in triangle
	GLfloat *p,*q;
	GLuint *t = triVerts(triangle);
	p = vertexCoordinate(t[0]);
	q = vertexCoordinate(t[1]);
	Vec3f u,v,xmp;
	u.set(q[0]-p[0],q[1]-p[1],q[2]-p[2]);
	q = vertexCoordinate(t[2]);
	v.set(q[0]-p[0],q[1]-p[1],q[2]-p[2]);
	xmp.set(xyz[0]-p[0],xyz[1]-p[1],xyz[2]-p[2]);
	float a,b,c,d;
	a=u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
	b=u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
	c=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
	if(fabs(d=b*b-a*c)<1e-16f)	{	// degenerate triangle
		uv[0]=0.0; uv[1]=0.0f;
		return; }
	uv[1] = ((u*b - v*a)*xmp)/d;
	uv[0] = (xmp*u - uv[1]*b)/a;
}

void glslTriangle::getBarycentricPosition(int triangle, float (&uv)[2], float (&xyz)[3])
{	// for barycentric uv in triangle returns position in xyz
	GLfloat *p,*q;
	GLuint *t = triVerts(triangle);
	p = vertexCoordinate(t[0]);
	q = vertexCoordinate(t[1]);
	Vec3f u,v,r;
	u.set(q[0]-p[0],q[1]-p[1],q[2]-p[2]);
	q = vertexCoordinate(t[2]);
	v.set(q[0]-p[0],q[1]-p[1],q[2]-p[2]);
	r = u*uv[0] + v*uv[1];
	xyz[0]=r.x()+p[0]; xyz[1]=r.y()+p[1]; xyz[2]=r.z()+p[2];
}

void glslTriangle::getBarycentricNormal(int triangle, float (&uv)[2], float (&normal)[3], bool uniqueNormalsNotCoords)
{	// for barycentric uv in triangle returns normal
	GLfloat p[3],q[3];
	GLuint *t = triVerts(triangle);
	getVertexNormal(t[0],p,uniqueNormalsNotCoords);
	getVertexNormal(t[1],q,uniqueNormalsNotCoords);
	Vec3f u,v,r;
	u.set(q[0]-p[0],q[1]-p[1],q[2]-p[2]);
	getVertexNormal(t[2],q,uniqueNormalsNotCoords);
	v.set(q[0]-p[0],q[1]-p[1],q[2]-p[2]);
	r = u*uv[0] + v*uv[1];
	normal[0]=r.x()+p[0]; normal[1]=r.y()+p[1]; normal[2]=r.z()+p[2];
}

void glslTriangle::getTriangleAdjacencies(unsigned int triangle, int (&adjacentTriangles)[3], int (&adjTriEdges)[3], bool uniqueNormalsNotCoords)
{
	if(uniqueNormalsNotCoords && !_adjacenciesComputed)
		findAdjacentTriangles(uniqueNormalsNotCoords);
	if(!uniqueNormalsNotCoords && !_adjacenciesComputedP)
		findAdjacentTriangles(uniqueNormalsNotCoords);
	unsigned long *adjs;
	if(uniqueNormalsNotCoords)
		adjs = &(_adjs[(triangle<<1)+triangle]);
	else
		adjs = &(_adjsP[(triangle<<1)+triangle]);
	for(int i=0; i<3; ++i)	{
		if(adjs[i]==3) {
			adjacentTriangles[i]=-1;
			adjTriEdges[i]=-1;
		}
		else	{
			adjTriEdges[i]=adjs[i]&0x00000003;
			adjacentTriangles[i]=adjs[i]>>2;
		}
	}
}

void glslTriangle::setVertexCoordinate(long vertex, const float (&newCoord)[3])
{	// get base coordinate then update all the copies
	vertex = _pnI[vertex].posIdx;
	GLfloat *v = &(_Positions[vertex]);
	v[0]=newCoord[0]; v[1]=newCoord[1]; v[2]=newCoord[2];
}

void glslTriangle::setVertexTexture(long vertex, const float (&newTex)[2])
{
	GLfloat *tx = &(_TexCoords[vertex<<1]);
	tx[0]=newTex[0]; tx[1]=newTex[1];
}

long glslTriangle::splitTriangleEdge(int triangle, int edge, const float parameter)
{	// Splits a triangle along edge(0-2) by parameter(0-1).
	// Creates 1 or 2 new triangles and one or two new vertices. Returns new vertex number of the input triangle split.
	// If 2 vertices created, user must use position index to find the second.
	// Does fix adjacencyP array so getNeighbors(,,false) works and texture interpolated.
	// Definitely should call getTriangleAdjacencies() ASAP.
	assert(0.0f<=parameter && 1.0f>=parameter);
	if(!_adjacenciesComputedP)	findAdjacentTriangles(false);
	unsigned long t1,newVert= addCoordNormTexVertices();
	GLuint *trVerts = &(_tris[triangle+(triangle<<1)]);
	float gv[3];
	GLfloat *gvp=&_Positions[_pnI[trVerts[edge]].posIdx];
	gv[0]=(1.0f-parameter)*gvp[0]; gv[1]=(1.0f-parameter)*gvp[1]; gv[2]=(1.0f-parameter)*gvp[2];
	gvp=&_Positions[_pnI[trVerts[(edge+1)%3]].posIdx];
	gv[0]+=parameter*gvp[0]; gv[1]+=parameter*gvp[1]; gv[2]+=parameter*gvp[2];
	setVertexCoordinate(newVert,gv);
	interpolateEdgeTextures(triangle,edge,newVert,parameter);
	unsigned long *trAdjs = &(_adjsP[triangle+(triangle<<1)]);
	unsigned long v[3];
	if(trAdjs[edge]==0x00000003)	{
		v[1] = trVerts[(edge+1)%3];
		trVerts[(edge+1)%3] = newVert;
		v[0]=newVert;	v[2]=trVerts[(edge+2)%3];
		t1 = addTriangle(v);	// invalidates old _tris and _adjsP pointers
		_adjsP[t1+(t1<<1)] = 0x00000003;
		_adjsP[t1+(t1<<1)+1] = _adjsP[triangle*3+(edge+1)%3];
		_adjsP[t1+(t1<<1)+2] = (triangle<<2)+((edge+1)%3);
		_adjsP[triangle+(triangle<<1)+(edge+1)%3] = (t1<<2)+2;
		if((_vertexFaceP[_pnI[v[1]].posIdx]&0x3fffffff)==triangle)
			_vertexFaceP[_pnI[v[1]].posIdx]= (t1 | 0x40000000);	// first vertex on a free edge
		return newVert;
	}
	unsigned long t2,v0,v1,va0,va1,newAdjVert,triAdj = trAdjs[edge]>>2;
	int tadjEdge = trAdjs[edge]&0x00000003;
	GLuint *tadjVerts = &(_tris[triAdj+(triAdj<<1)]);
	va0 = tadjVerts[tadjEdge];
	va1 = tadjVerts[(tadjEdge+1)%3];
	v0 = trVerts[edge];
	v1 = trVerts[(edge+1)%3];
	if(va1==v0 && va0==v1)	{
		newAdjVert = newVert;
	}
	else if(_pnI[v0].posNormIdx==_pnI[va1].posNormIdx && _pnI[v1].posNormIdx==_pnI[va0].posNormIdx)	{	// texture seam
		newAdjVert = addTexVertex(newVert);
		interpolateEdgeTextures(triAdj,tadjEdge,newAdjVert,1.0f-parameter);
	}
	else if(_pnI[v0].posIdx==_pnI[va1].posIdx && _pnI[v1].posIdx==_pnI[va0].posIdx)	{	// hard normal seam
		newAdjVert = addNormTexVertex(newVert);
		interpolateEdgeTextures(triAdj,tadjEdge,newAdjVert,1.0f-parameter);
	}
	else
		assert(false);
	trVerts = &_tris[triangle+(triangle<<1)];
	tadjVerts = &_tris[triAdj+(triAdj<<1)];
	unsigned long *tadjAdjs = &_adjsP[triAdj+(triAdj<<1)];
	trVerts[(edge+1)%3] = newVert;
	tadjVerts[(tadjEdge+1)%3] = newAdjVert;
	unsigned long ata,at=_adjsP[triangle+(triangle<<1)+((edge+1)%3)];
	ata=_adjsP[triAdj+(triAdj<<1)+((tadjEdge+1)%3)];
	v[0]=newVert; v[1]=v1; v[2]=trVerts[(edge+2)%3]; 
	t1=addTriangle(v);	// invalidates old _tris and _adjsP pointers
	if(_vertexFaceP[_pnI[v1].posIdx>>2]==triangle)
		_vertexFaceP[_pnI[v1].posIdx>>2] = t1;
	tadjVerts = &_tris[triAdj+(triAdj<<1)];
	v[0]=newAdjVert;	v[1]=va1;	v[2]=tadjVerts[(tadjEdge+2)%3];
	t2=addTriangle(v);	// invalidates old _tris and _adjsP pointers
	if(_vertexFaceP[_pnI[va1].posIdx>>2]==triAdj)
			_vertexFaceP[_pnI[va1].posIdx>>2] = t2;
	_adjsP[t1+(t1<<1)] = (triAdj<<2)+tadjEdge;
	_adjsP[t1+(t1<<1)+1] = at;
	_adjsP[t1+(t1<<1)+2] = (triangle<<2) + ((edge+1)%3);
	_adjsP[t2+(t2<<1)] = (triangle<<2) + edge;
	_adjsP[t2+(t2<<1)+1] = ata;
	_adjsP[t2+(t2<<1)+2] = (triAdj<<2)+((tadjEdge+1)%3);
	_adjsP[triangle+(triangle<<1)+edge] = t2<<2;
	_adjsP[triangle+(triangle<<1)+((edge+1)%3)] = (t1<<2)+2;
	_adjsP[triAdj+(triAdj<<1)+tadjEdge] = t1<<2;
	_adjsP[triAdj+(triAdj<<1)+((tadjEdge+1)%3)] = (t2<<2)+2;
	if(at!=0x00000003)
		_adjsP[(at>>2)*3+(at&0x00000003)] = (t1<<2)+1;
	if(ata!=0x00000003)
		_adjsP[(ata>>2)*3+(ata&0x00000003)] = (t2<<2)+1;
//	if(_vertexFaceP.size()>(unsigned)newVert)	// which one doesn't matter as we are in a closed ring
	assert(_vertexFaceP.size()>(unsigned)_pnI[newVert].posIdx>>2);
	_vertexFaceP[_pnI[newVert].posIdx>>2]=triangle;
	if(newAdjVert>newVert)	// if a second new vertex was added
		_vertexFaceP[_pnI[newAdjVert].posIdx>>2]=t2;
	return newVert;
}

void glslTriangle::interpolateEdgeTextures(int triangle, int edge, int newVert, float param)
{	// assumes triangle hasn't been changed yet by newVert
	GLuint *trVerts = &_tris[triangle+(triangle<<1)];
	GLfloat *txOut=&_TexCoords[newVert<<1],*txIn;
	txIn=&_TexCoords[trVerts[edge]<<1];
	txOut[0]=(1.0f-param)*txIn[0]; txOut[1]=(1.0f-param)*txIn[1];
	txIn=&_TexCoords[trVerts[(edge+1)%3]<<1];
	txOut[0]+=param*txIn[0]; txOut[1]+=param*txIn[1];
}

long glslTriangle::addNewVertexInMidTriangle(int triangle, const float (&uvParameters)[2])
{	// creates 2 new triangles and one new vertex. Returns new vertex number.
	// Input uvParameters[2] are parameters along vectors t1-t0 and t2-t0.
	// Does fix adjacencyP array so getNeighbors(,,false) works and texture&positions interpolated.
	// Definitely should getTriangleAdjacencies() ASAP.
	assert(uvParameters[0]>=0.0f && uvParameters[0]<=1.0f);
	assert(uvParameters[1]>=0.0f && uvParameters[1]<=1.0f);
	assert(uvParameters[0]+uvParameters[1]<=1.0001f);
	GLuint *trVerts = &_tris[triangle+(triangle<<1)];
	if(uvParameters[0]<0.0002f && uvParameters[1]<0.0002f)
		return trVerts[0];
	if(uvParameters[0]>0.9998f)
		return trVerts[1];
	if(uvParameters[1]>0.9998f)
		return trVerts[2];
	if(uvParameters[0]<0.02f)
		return splitTriangleEdge(triangle,2,1.0f-uvParameters[1]);
	if(uvParameters[1]<0.02f)
		return splitTriangleEdge(triangle,0,uvParameters[0]);
	if(uvParameters[0]+uvParameters[1]>0.98f)
		return splitTriangleEdge(triangle,1,1.0f-uvParameters[0]);
	if(!_adjacenciesComputedP)	findAdjacentTriangles(false);
	unsigned long *trAdjs = &_adjsP[triangle+(triangle<<1)];
	unsigned long a1=trAdjs[1],a2=trAdjs[2];
	long t1,t2,oldVert,ret= addCoordNormTexVertices();
	GLfloat p,*pv = &_Positions[_pnI[trVerts[0]].posIdx];
	p = 1.0f - uvParameters[0] - uvParameters[1];
	float vec[3] = {p*pv[0],p*pv[1],p*pv[2]};
	pv =&_TexCoords[trVerts[0]<<1];
	float tx[2]={p*pv[0],p*pv[1]};
	for(int i=0; i<2; ++i)	{
		pv =&_TexCoords[trVerts[i+1]<<1];
		tx[0] += pv[0]*uvParameters[i];
		tx[1] += pv[1]*uvParameters[i];
		pv = &_Positions[_pnI[trVerts[i+1]].posIdx];
		vec[0]+=pv[0]*uvParameters[i]; vec[1]+=pv[1]*uvParameters[i]; vec[2]+=pv[2]*uvParameters[i];
	}
	setVertexCoordinate(ret,vec);
	setVertexTexture(ret,tx);
	oldVert = trVerts[2];
	trVerts[2] = ret;
	unsigned long v[3];
	v[0]=ret; v[1]=trVerts[1]; v[2]=oldVert; 
	t1=addTriangle(v);	// invalidates old tris and adjs pointers
	v[1]=oldVert; v[2]=_tris[triangle+(triangle<<1)];
	t2=addTriangle(v);
	oldVert = _pnI[oldVert].posIdx>>2;
	if((_vertexFaceP[oldVert]&0x3fffffff)==triangle)	{
		_vertexFaceP[oldVert]=t2;
		if(a2==0x00000003)
			_vertexFaceP[oldVert] |= 0x40000000;
	}
	oldVert = _tris[triangle+(triangle<<1)+1];
	oldVert = _pnI[oldVert].posIdx>>2;
	if((_vertexFaceP[oldVert]&0x3fffffff)==triangle)	{
		_vertexFaceP[oldVert]=t1;
		if(a1==0x00000003)
			_vertexFaceP[oldVert] |= 0x40000000;
	}
	_adjsP[t1+(t1<<1)] = (triangle<<2)+1;
	_adjsP[t1+(t1<<1)+1] = a1;
	_adjsP[t1+(t1<<1)+2] = t2<<2;
	_adjsP[t2+(t2<<1)] = (t1<<2)+2;
	_adjsP[t2+(t2<<1)+1] = a2;
	_adjsP[t2+(t2<<1)+2] = (triangle<<2)+2;
	_adjsP[triangle+(triangle<<1)+1] = t1<<2;
	_adjsP[triangle+(triangle<<1)+2] = (t2<<2)+2;
	if(a1!=0x00000003)
		_adjsP[(a1>>2)*3+(a1&0x00000003)] = (t1<<2)+1;
	if(a2!=0x00000003)
		_adjsP[(a2>>2)*3+(a2&0x00000003)] = (t2<<2)+1;
	_vertexFaceP[_pnI[ret].posIdx>>2]=triangle;
	return ret;
}

void glslTriangle::makeLineList(std::vector<GLuint> &lines)
{	// makes a line strip list using 0xffffffff as a primitive restart indicator
	// Indices output are into the _Positions array.
	typedef std::set<edge,edgeTest> edgeSet;
	typedef edgeSet::iterator edgeIt;
	std::pair <edgeIt,bool> P;
	edge E;
	edgeIt ei;
	edgeSet M;
	M.clear();
	unsigned long i,j,tnow[3],numtris=(unsigned long)_tris.size();
	for(i=0; i<numtris; i+=3)
	{
		if(_tris[i] == 0xffffffff)	// signals a deleted triangle
			continue;
		for(j=0; j<3; j++) {
			tnow[j] = _tris[i+j];
			tnow[j] = _pnI[tnow[j]].posIdx>>2;
		}
		for(j=0; j<3; j++)
		{
			unsigned int tmp;
			if((tmp=tnow[(j+1)%3])<tnow[j]) {
				E.vtxMin = tmp;
				E.vtxMax = tnow[j]; }
			else {
				E.vtxMin = tnow[j];
				E.vtxMax = tmp; }
			P = M.insert(E);
		}
	}
	std::list<std::list<GLuint> > strips;
	std::map<GLuint,std::list<std::list<GLuint> >::iterator> heads,tails;
	std::map<GLuint,std::list<std::list<GLuint> >::iterator>::iterator fit;
	std::list<std::list<GLuint> >::iterator h0,h1,t0,t1;
	int tests;
	for(ei=M.begin(); ei!=M.end(); ++ei)	{
		tests=0;
		if((fit=heads.find(ei->vtxMin))!=heads.end())	{
			tests|=0x01;
			h0=fit->second;	}
		if((fit=heads.find(ei->vtxMax))!=heads.end())	{
			tests|=0x04;
			h1=fit->second;	}
		if((fit=tails.find(ei->vtxMin))!=tails.end())	{
			tests|=0x02;
			t0=fit->second;	}
		if((fit=tails.find(ei->vtxMax))!=tails.end())	{
			tests|=0x08;
			t1=fit->second;	}
		if(tests<1)	{
			strips.push_front(std::list<GLuint>());
			strips.front().push_back(ei->vtxMin);
			strips.front().push_back(ei->vtxMax);
			heads.insert(std::make_pair((unsigned long)ei->vtxMin,strips.begin()));
			tails.insert(std::make_pair(ei->vtxMax,strips.begin()));
		}
		else if(tests<4)	{
			if(tests<2)	{
				heads.erase(h0->front());
				h0->push_front(ei->vtxMax);
				heads.insert(std::make_pair(ei->vtxMax,h0));	}
			else	{
				tails.erase(t0->back());
				t0->push_back(ei->vtxMax);
				tails.insert(std::make_pair(ei->vtxMax,t0));	}
		}
		else if((tests&0x03)<1)	{
			if(tests<8)	{
				heads.erase(h1->front());
				h1->push_front(ei->vtxMin);
				heads.insert(std::make_pair((unsigned long)ei->vtxMin,h1));	}
			else	{
				tails.erase(t1->back());
				t1->push_back(ei->vtxMin);
				tails.insert(std::make_pair((unsigned long)ei->vtxMin,t1));	}
		}
		else	{
			if(tests&0x01)	{
				heads.erase(h0->front());
				tails.erase(h0->back());
				t0 = h0;
				t0->reverse();
				heads.insert(std::make_pair(t0->front(),t0));
			}
			else
				tails.erase(t0->back());
			if(tests&0x08)	{
				heads.erase(t1->front());
				tails.erase(t1->back());
				h1 = t1;
				h1->reverse();
			}
			else	{
				heads.erase(h1->front());
				tails.erase(h1->back());	}
			if(t0==h1)	{
				t0->push_back(t0->front());
				heads.insert(std::make_pair(t0->front(),t0));
				tails.insert(std::make_pair(t0->back(),t0));
			}
			else	{
				t0->splice(t0->end(),*h1);
				strips.erase(h1);
				tails.insert(std::make_pair(t0->back(),t0));
			}
		}
	}
	heads.clear();
	tails.clear();
	M.clear();
	lines.clear();
	tests=0;
	for(h0=strips.begin(); h0!=strips.end(); ++h0)
		tests += (int)(h0->size())+1;
	lines.reserve(tests);
	for(h0=strips.begin(); h0!=strips.end(); ++h0)	{
		lines.insert(lines.end(),h0->begin(),h0->end());
		lines.push_back(0xffffffff);
	}
	if(!lines.empty())
		lines.pop_back();
}

void glslTriangle::cleanAndPack(std::vector<int> &newVertexMap, std::vector<int> &newTriangleMap)
{	// warning - invalidates all triangle and vertex indices. Returns the new index numbersof the old vertices and triangles. -1 indicates deletion.
	std::vector<int> pp,pn;
	newVertexMap.clear();
	newVertexMap.assign(_pnI.size(),-1);
	pp.assign(_Positions.size()>>2,-1);
	pn.assign(_uniqueNormalNumber,-1);
	int i,j,k,n=(int)_tris.size();	//,maxPN=0;;
	for(i=0; i<n; i+=3)	{
		if(_tris[i]==0xffffffff)	// deleted triangle
			continue;
		for(j=0; j<3; ++j)	{
			newVertexMap[_tris[i+j]] = 1;
			pp[_pnI[_tris[i+j]].posIdx>>2] = 1;
			pn[_pnI[_tris[i+j]].posNormIdx] = 1;
		}
	}
	n=(int)pp.size();
	j = 0;
	for(i=0; i<n; ++i)	{
		if(pp[i]<0)
			continue;
		if(i>j)	{
			_Positions[j<<2]=_Positions[i<<2];
			_Positions[(j<<2)+1]=_Positions[(i<<2)+1];
			_Positions[(j<<2)+2]=_Positions[(i<<2)+2];
			_Positions[(j<<2)+3]=_Positions[(i<<2)+3];
		}
		pp[i]=j;
		++j;
	}
	_Positions.resize(j<<2);
	n=(int)pn.size();
	j = 0;
	for(i=0; i<n; ++i)	{
		if(pn[i]>-1)	{
			pn[i] = j;
			++j;
		}
	}
	_uniqueNormalNumber=j;
	n=(int)newVertexMap.size();
	j = 0;
	for(i=0; i<n; ++i)	{
		if(newVertexMap[i]<0)
			continue;
		_pnI[j].posIdx=pp[_pnI[i].posIdx>>2]<<2;
		_pnI[j].posNormIdx=pn[_pnI[i].posNormIdx];
		_TexCoords[j<<1]=_TexCoords[i<<1];
		_TexCoords[(j<<1)+1]=_TexCoords[(i<<1)+1];
		newVertexMap[i]=j;
		++j;
	}
	pn.clear();
	_pnI.resize(j);
	_vertexNumber=j;
	_TexCoords.resize(j<<1);
	n=(int)_tris.size();
	newTriangleMap.clear();
	newTriangleMap.assign(n/3,-1);
	k=0;
	for(i=0; i<n; i+=3)	{
		if(_tris[i]==0xffffffff)	// deleted triangle
			continue;
		for(j=0; j<3; ++j)	{
			assert(newVertexMap[_tris[i+j]]>-1);
			_tris[k+j] = newVertexMap[_tris[i+j]];	}
		newTriangleMap[i/3] = k/3;
		k += 3;
	}
	_tris.resize(k);
	_triangleNumber=k/3;
	_adjacenciesComputed=false;
	_adjacenciesComputedP=false;
	_adjs.clear();
	_adjsP.clear();
	_vertexFace.clear();
	_vertexFaceP.clear();
}

void glslTriangle::getNearestHardEdge(float (&xyz)[3], int &triangle, int &edge, float &param)
{	// Input first 2 args, then overwrites all 4 with the point on the nearest hard edge.
	// if triangle<0 search all triangles, else search only side containing input triangle.
	std::set<int> tris;
	std::vector<int> hardEdges;
	if(!_adjacenciesComputed)
			findAdjacentTriangles(true);
	if(triangle<0)	{
		int i,j,n=(int)_tris.size();
		for(i=0; i<n; i+=3)	{
			for(j=0; j<3; ++j)	{
				if(_adjs[i+j]==3)	{
					if(tris.find(i/3)==tris.end())
						recurseHardEdgeSearch(i/3,tris,hardEdges);
					break;
				}
 			}
		}
	}
	else
		recurseHardEdgeSearch(triangle,tris,hardEdges);
	const GLfloat *c;
	float denom,t,minD=1e32f;
	Vec3f *P,Q,R,S,E;
	P = (Vec3f*)&xyz;
	unsigned long i,j,n=(unsigned long)hardEdges.size();
	for(i=0; i<n; ++i)	{
		j = hardEdges[i];
		c = vertexCoordinate(triVerts(j>>2)[j&3]);
		Q.set(c[0],c[1],c[2]);
		c = vertexCoordinate(triVerts(j>>2)[((j&3)+1)%3]);
		R.set(c[0],c[1],c[2]);
		R -= Q;
		denom = R*R;
		if(denom<1e-16)
			continue;
		t = (*P-Q)*R/denom;
		if(t>1.0f)	t=1.0f;
		if(t<0.0f)	t=0.0f;
		S = Q + R*t;
		if((denom=(*P-S).length2())<minD)	{
			minD = denom;
			param = t;
			E = S;
			triangle = j>>2;
			edge = j&3;
		}
	}
	xyz[0]=E.x(); xyz[1]=E.y(); xyz[2]=E.z();
}

void glslTriangle::recurseHardEdgeSearch(int triangle, std::set<int> &visited, std::vector<int> &hardEdges)
{
	if(!visited.insert(triangle).second)
		return;
	int adjT[3],adjE[3];
	getTriangleAdjacencies(triangle,adjT,adjE,true);
	for(int i=0; i<3; ++i)	{
		if(adjT[i]<0)
			hardEdges.push_back((triangle<<2)+i);
		else	{
			recurseHardEdgeSearch(adjT[i],visited,hardEdges);
		}
	}
}

void glslTriangle::getVertexCoordinate(unsigned int vertex, float (&xyz)[3])
{	// type safe version
	GLfloat *v = &_Positions[_pnI[vertex].posIdx];
	xyz[0]=v[0]; xyz[1]=v[1]; xyz[2]=v[2];
}

int glslTriangle::isManifoldConsistent()
{	// topology checker.  Returns # of topological handles or -30000 if inconsistent
	typedef std::set<edge,edgeTest> edgeSet;
	edge E;
	edgeSet M;
	M.clear();
	unsigned long tnow[3];
	unsigned long i,j,numtris=(unsigned long)_tris.size();
	std::vector<int> posVec;
	posVec.assign(_Positions.size()>>2,0);
	for(i=0; i<numtris; i+=3)
	{
		if(_tris[i] == 0xffffffff)	// signals a deleted triangle
			continue;
		for(j=0; j<3; j++)
			tnow[j] = _pnI[_tris[i+j]].posIdx>>2;
		for(j=0; j<3; j++)
		{
			if(tnow[j]>=posVec.size())
				return -30000;
			else
				posVec[tnow[j]] = 1;
			unsigned int tmp;
			if((tmp=tnow[(j+1)%3])<tnow[j]) {
				E.vtxMin = tmp;
				E.vtxMax = tnow[j];	}
			else {
				E.vtxMin = tnow[j];
				E.vtxMax = tmp;	}
			M.insert(E);
		}
	}
	for(i=0; i<posVec.size(); ++i)
		if(!posVec[i])
			return -30000;
	int handles2 = (int)((_tris.size()/3) + posVec.size() - M.size());
	if(handles2&0x0001)
		return -30000;
	return handles2>>1;
}

