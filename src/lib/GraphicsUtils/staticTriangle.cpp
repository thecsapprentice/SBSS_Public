// staticTriangle.cpp
// Author: Court Cutting
// Date: 5/20/2014
// Purpose: Triangle object management class for static triangulated objects.
//    Meant to complement triangleUVW class for management of dynamic objects.
//    Topology and textures should be sent anly once. Vertex positions positions
//    however, may be changed and normals and tangents can be recomputed.
//       Copyright 2014 - All rights reserved.

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
#include "wxGraphics.h"
#include "staticTriangle.h"

//////////////////////// TEMPORARY TEMPORARY TEMPORARY - On SnowLeopard this is suppored, but GLEW doens't hook up properly
//////////////////////// Fixed probably in 10.6.3
#ifdef __APPLE__
#define glGenVertexArrays glGenVertexArraysAPPLE
#define glDeleteVertexArrays  glDeleteVertexArraysAPPLE
#define glBindVertexArray	glBindVertexArrayAPPLE
#endif

const GLchar *staticTriangle::staticVertexShader = "#version 140 \n"
	"in vec4 vVertex;"
	"in vec3 vNormal;"
	"in vec3 vTangent;"
	"in vec2 vTexture;"
	"uniform mat4   mvpMatrix;"
	"uniform mat4   mvMatrix;"
	"uniform mat3   normalMatrix;"
	"uniform vec3   vLightPosition;"
	// next 2 are for bump mapping per Randi Rost
	"smooth out vec3 vLightDir;"
	"smooth out vec3 vEyeDir;"
	"smooth out vec2 vTexCoords;"
	"void main(void) "
	"{"
	"   vEyeDir = vec3(mvMatrix * vVertex);"
	"	vec3 n = normalize(normalMatrix * vNormal);\n"
	"	vec3 t = normalize(normalMatrix * vTangent);\n"
	"	vec3 b = cross(n,t);\n"
	"	vec3 v;\n"
	"	v.x = dot(vLightPosition,t);\n"
	"	v.y = dot(vLightPosition,b);\n"
	"	v.z = dot(vLightPosition,n);\n"
	"	vLightDir = normalize(v);\n"
	"	v.x = dot(vEyeDir,t);\n"
	"	v.y = dot(vEyeDir,b);\n"
	"	v.z = dot(vEyeDir,n);"
	"	vEyeDir = normalize(v);"
	"	vTexCoords = vTexture;"
	"   gl_Position = mvpMatrix * vVertex;"
	"}";

const GLchar *staticTriangle::staticFragmentShader =
	// Adapted from Randi Rost and Bill Licea-Kane
	// OpenGL Shading Language 3rd edition
	"#version 140 \n"
	"out vec4 vFragColor;"
	"uniform vec4 ambientColor;"
	"uniform vec4 diffuseColor;"
	"uniform sampler2D colorMap;"
	"uniform sampler2D normalMap;"
	"smooth in vec3 vLightDir;"
	"smooth in vec3 vEyeDir;"
	"smooth in vec2 vTexCoords;"
	"void main(void)"
	"{"
	// always lighting with white light
	"	const float ambientVal = 0.1;"
	"	const float	fatIncr = 0.5/1024.0;"
	"   vec2 fatD = vec2(fatIncr,0.0);"
	"   vec2 faceUV;"
	"	float lightVal = ambientVal;"
	"	float sn,h,specMult = 0.5;"
	"	const float diffuseVal = 0.9;"
	"	vec3 normDelta = vec3(0.0, 0.0, 1.0);"
	"	vec3 litColor = vec3(1.0, 1.0, 1.0);"
// COURT for debug - remove later
//	"    if(vTexCoords.s>4.0f) "
//	"	    vFragColor = vec4(1.0, 1.0, 0.0, 1.0);"
//	"	else {"	// top of skin
	"	  vec4 tx1 = texture(normalMap,vTexCoords.st);"
	"	  tx1.rgb -= vec3(0.5);"
	"	  normDelta = tx1.rgb*2.0;"
	"	  specMult = 0.11;"
	"	  vFragColor = texture(colorMap, vTexCoords.st);"
	"	lightVal += diffuseVal*max(dot(normDelta,vLightDir), 0.0);"
	"	vFragColor *= lightVal;"
	"	vec3 reflectDir = reflect(vLightDir,normDelta);"
	"	float spec = max(dot(vEyeDir,reflectDir),0.0);"
	"	spec = pow(spec,40.0);"
	"	spec *= specMult;"
	"	litColor = min(vFragColor.rgb + spec, vec3(1.0));"
	"	vFragColor = vec4(litColor, 1.0);"
	"}";

void staticTriangle::deleteTriangle(long triangleNumber)
{	// may strand vertices
	for(int i=triangleNumber*3; i<_tris.size()-3; ++i)
		_tris[i] = _tris[i+3];
	_tris.pop_back(); _tris.pop_back(); _tris.pop_back();
	--_triangleNumber;
}

long staticTriangle::addTriangle(unsigned long (&vertices)[3])
{
	_tris.push_back((GLuint)vertices[0]);
	_tris.push_back((GLuint)vertices[1]);
	_tris.push_back((GLuint)vertices[2]);
	_triangleNumber = (unsigned long)(_tris.size()/3);
	if(_adjacenciesComputed)	{ // done for incision class to make space
		// Really _adjacenciesComputedP should be set false when a new triangle is added.
		_adjs.push_back(0x00000003);
		_adjs.push_back(0x00000003);
		_adjs.push_back(0x00000003); }
	return (long)_triangleNumber-1;
}

long staticTriangle::addCoordNormTexVertices(int numberToAdd)
{
	long retval = _vertexNumber;
	for(unsigned int i=0; i<(unsigned)numberToAdd; ++i)
	{
//		_pnI.push_back(positionNormalIndex());
//		_pnI.back().posIdx=(long)_Positions.size();
//		_pnI.back().posNormIdx=_uniqueNormalNumber+i;
		_xyz1.push_back(0.0f);
		_xyz1.push_back(0.0f);
		_xyz1.push_back(0.0f);
		_xyz1.push_back(1.0f);
		_TexCoords.push_back(0.0f);
		_TexCoords.push_back(0.0f);
		if(_adjacenciesComputed)
			_vertexFace.push_back(0x80000000);
	}
	_vertexNumber += numberToAdd;
	_uniqueNormalNumber += numberToAdd;
	return retval;
}

staticTriangle::staticTriangle(bool texturedNotColored) {
	sceneNode::_coloredNotTextured=!texturedNotColored;
	_textured=texturedNotColored;
	_computeTangents=texturedNotColored;
	setType(TRIANGLES);
	_vertexNumber = 0;
	_uniqueNormalNumber = 0;
//	_Verts.clear();
//	_Norms.clear();
	_TexCoords.clear();
	_xyz1.clear();
	_adjacenciesComputed=false;
    _nMaxIndexes = 0;
    _nNumIndexes = 0;
    _nNumVerts = 0;
	_vertexArrayBufferObject=0;
	for(int i=0; i<5; ++i)
		_bufferObjects[i]=0;
}

void staticTriangle::setTexturedNotColored(bool texturedNotColored) {
	sceneNode::_coloredNotTextured=!texturedNotColored;
	_textured=texturedNotColored;
	_computeTangents=texturedNotColored;
}

staticTriangle::~staticTriangle() {
    // Just in case these still are allocated when the object is destroyed
	_TexCoords.clear();
	_xyz1.clear();
	//if(_bufferObjects[0]>0)
	//	glDeleteBuffers(4, _bufferObjects);
	// Cleanup textures
	//if(_vertexArrayBufferObject>0)
	//	glDeleteVertexArrays(1, &_vertexArrayBufferObject);
}

bool staticTriangle::getTriangleVertices(unsigned int triangle, int (&vertices)[3])
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

int staticTriangle::getClosestVertex(float (&position)[3], int triangle)
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

void staticTriangle::draw(void) 
{
	//glBindVertexArray(_vertexArrayBufferObject);
	//if(_textured) {
	//	glActiveTexture( GL_TEXTURE0);
	//	glBindTexture(GL_TEXTURE_2D, _2DtextureBuffer);
	//	glActiveTexture( GL_TEXTURE1);
	//	glBindTexture(GL_TEXTURE_2D, _2DnormalBuffer);
	//}
	//glDrawElements(GL_TRIANGLES, (GLsizei)_tris.size(), GL_UNSIGNED_INT, 0);
    // Never unbind a GL_ARRAY_BUFFER or GL_ELEMENT_ARRAY_BUFFER inside an active vertex array buffer object
	//glBindVertexArray(0);

//GLenum errCode;
//if((errCode=glGetError())!=GL_NO_ERROR)
//	errCode=errCode;

}    

/* void staticTriangle::setVerticesComputeNormals()
{
	// Vertex data
    glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[0]);	// VERTEX_DATA
	GLfloat *vtx,*fptr;
	fptr = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER,GL_WRITE_ONLY);	//|GL_MAP_INVALIDATE_BUFFER_BIT|GL_MAP_FLUSH_EXPLICIT_BIT);
	size_t n = _pnI.size();
	for(unsigned int j,i=0; i<n; ++i)	{
		vtx = &(_xyz1[_pnI[i].posIdx]);
		for(j=0; j<4; ++j)
			*(fptr++) = *(vtx++);
	}
	glUnmapBuffer(GL_ARRAY_BUFFER);
//		glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[0]);
//		std::vector<GLfloat> fbCoords;
//		fbCoords.assign(_vertexNumber*4,0.0f);
//		glGetBufferSubData(GL_ARRAY_BUFFER,0,sizeof(GLfloat)*_vertexNumber*4,&(fbCoords[0]));
	glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[1]);
	fptr = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER,GL_WRITE_ONLY);	//|GL_MAP_INVALIDATE_BUFFER_BIT|GL_MAP_FLUSH_EXPLICIT_BIT);
	float normal[3];
	for(unsigned int j,i=0; i<n; ++i)	{
		getVertexNormal(i,normal);
		for(j=0; j<3; ++j)
			*(fptr++) = normal[j];
	}
	glUnmapBuffer(GL_ARRAY_BUFFER);
//		glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[1]);
//	fbCoords.assign(_vertexNumber*3,0.0f);
//	glGetBufferSubData(GL_ARRAY_BUFFER,0,sizeof(GLfloat)*_vertexNumber*3,&(fbCoords[0]));
	glBindBuffer(GL_ARRAY_BUFFER,0);
//	glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[1]);
//		std::vector<GLfloat> fbCoords;
//		fbCoords.assign(_vertexNumber*3,0.0f);
//		glGetBufferSubData(GL_ARRAY_BUFFER,0,sizeof(GLfloat)*_vertexNumber*3,&(fbCoords[0]));
//		glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[0]);
//		fbCoords.assign(_vertexNumber*4,0.0f);
//		glGetBufferSubData(GL_ARRAY_BUFFER,0,sizeof(GLfloat)*_vertexNumber*4,&(fbCoords[0]));
} */

void staticTriangle::computeLocalBounds()
{
	float minMaxXYZ[6];
	minMaxXYZ[0]=minMaxXYZ[2]=minMaxXYZ[4]=1e30f;
	minMaxXYZ[1]=minMaxXYZ[3]=minMaxXYZ[5]=-1e30f;
	int i,n=(unsigned int)_xyz1.size()>>2;
	GLfloat *v=&(_xyz1[0]);
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

//void staticTriangle::getLocalBounds(GLfloat (&localCenter)[3], GLfloat &Radius) {
//	localCenter[0]=_localCenter[0]; localCenter[1]=_localCenter[1]; localCenter[2]=_localCenter[2];
//	Radius = _radius;
//}

bool staticTriangle::parseNextInputFileLine(std::stringstream *infile, std::string &unparsedLine, std::vector<std::string> &parsedLine)
{
	if(infile->eof())
	{
		return false;
	}
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

/*bool staticTriangle::writeObjFile(const char *fileName)
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
	for(size_t i=0; i<_xyz1.size(); i+=4)	{
		sprintf(s,"v %f %f %f\n",_xyz1[i],_xyz1[i+1],_xyz1[i+2]);
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
} */

int staticTriangle::readObjFile(const char *data, bool dataOnlyNoGraphics)
{ // returned error codes: 0=no error, 1=can't open file, 2=non-triangle primitive,
	// 3=exceeded 0x7fff vertex position limit, 4=bad 3D vertex line, 5=bad 2D texture line, 6=nonmanifold surface
	_textured=false;
    std::stringstream fin;
    fin << std::string(data);
    fin.seekg( std::ios::beg );
    std::string unparsedLine;
	std::vector<std::string> parsedLine;
	_vertexNumber=0;
	_triangleNumber=0;
	_uniqueNormalNumber=0;
	cnt cntin;
	CNTMAP cntmap;
	cntmap.clear();
	std::pair<CNTMAP::iterator,bool> cpr;
	std::vector<GLfloat> tt,tp;
	tt.clear(); tp.clear();
	int i,j,k,l;
	std::string str;
	char s[100];
	std::vector<GLuint> tempTris;
	while(parseNextInputFileLine(&fin,unparsedLine,parsedLine))
	{
		if(parsedLine.empty() || parsedLine[0]=="#")
			continue;
		if(parsedLine[0]=="v")
		{
			if(parsedLine.size()!=4)
				return 4;
			for(i=1; i<4; ++i)
				tp.push_back((GLfloat)atof(parsedLine[i].c_str()));
		}
		else if(parsedLine[0]=="vt")
		{
			_textured=true;
			if(parsedLine.size()!=3)
				return 5;
			tt.push_back((GLfloat)atof(parsedLine[1].c_str()));
			tt.push_back((GLfloat)atof(parsedLine[2].c_str()));
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
		}
		else
			continue;
	}
	if(	_vertexNumber>0x7fff) {
		return 3;
	}
	_tris.clear();
	_tris.assign(tempTris.begin(),tempTris.end());
	tempTris.clear();
	_xyz1.clear();
	_xyz1.assign(cntmap.size()<<2,1.0f);
	_Normals.clear();
	_Normals.assign(cntmap.size()*3,-1.0f);
	_pnI.clear();
	_pnI.assign(cntmap.size(),-1);
	if(_uniqueNormalNumber<1)
		_uniqueNormalNumber = (unsigned long)cntmap.size();
	_TexCoords.clear();
	_Tangents.clear();
	if(_textured) {
		_TexCoords.assign(cntmap.size()<<1,-1.0f);
	_Tangents.assign(cntmap.size()*3,-1.0f); }
	CNTMAP::iterator cit=cntmap.begin();
	while(cit!=cntmap.end()) {
		_pnI[cit->second] = cit->first.v[0];
		_xyz1[cit->second<<2] = tp[cit->first.v[0]*3];
		_xyz1[(cit->second<<2)+1] = tp[cit->first.v[0]*3+1];
		_xyz1[(cit->second<<2)+2] = tp[cit->first.v[0]*3+2];
		if(_textured)	{
			_TexCoords[cit->second<<1] = tt[(cit->first.v[1])<<1];
			_TexCoords[(cit->second<<1)+1] = tt[((cit->first.v[1])<<1)+1]; }
		++cit;
	}
	cntmap.clear();
	if(!findAdjacentTriangles()) {
		_vertexNumber=0;
		_triangleNumber=0;
		_uniqueNormalNumber=0;
		_tris.clear();
		_xyz1.clear();
		_TexCoords.clear();
		return 6;
	}
//	if(!dataOnlyNoGraphics)
//		sendToGraphicsCard();
	return 0;
}

bool staticTriangle::findAdjacentTriangles()
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
	_adjs.clear();
	_adjs.assign(numtris,0x00000003);
	for(i=0; i<numtris; i+=3)
	{
		if(_tris[i] == 0xffffffff)	// signals a deleted triangle
			continue;
		for(j=0; j<3; j++) {
			tnow[j] = _tris[i+j];
//			tnow[j] = _pnI[tnow[j]].posNormIdx;
		}
		adjNow = &(_adjs[i]);
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
				_adjs[(tcode>>2)*3+(tcode&0x00000003)] = E.adjCode;
				M.erase(P.first);
			}
			else
				adjNow[j] = 0x00000003;
		}
	}
	M.clear();
	makeVertexToTriangleMap();
	_adjacenciesComputed = 1;
	return true;
}

int staticTriangle::getVertexTriangle(int vertexNumber)	// gets triangle vertex is a member of
{
	if(_adjacenciesComputed)
		return _vertexFace[vertexNumber]&0x3fffffff;
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

void staticTriangle::makeVertexToTriangleMap()
{
	int i,j,numtris=(int)_tris.size();
	unsigned int *tnow,vnow;
	// provide each vertex with a face it is a member of
	_vertexFace.clear();
	_vertexFace.assign(_vertexNumber,0x80000000);	// initially deleted
	for(i=0; i<numtris; i+=3)
	{
		tnow = &(_tris[i]);
		if(tnow==NULL || tnow[0] == 0xffffffff)	// signals a deleted triangle
			continue;
		for(j=0; j<3; j++)
		{
			vnow = tnow[j];
//			vnow = _pnI[vnow].posNormIdx;
			if(_vertexFace[vnow]&0x40000000)
				continue;	// vertex first on free edge, don't change
			_vertexFace[vnow] = i/3;
			if(_adjs[i+j]==0x00000003)	// vertex first on free edge, lock it for easy neighbor find
				_vertexFace[vnow] |= 0x40000000;
		}
	}
}

void staticTriangle::getNeighbors(unsigned int vertex, std::vector<neighborNode> &neighbors)
{
	unsigned long normVert,trNum,adj,triStart;
	if(!_adjacenciesComputed)
		findAdjacentTriangles();
//	vertex = _pnI[vertex].posNormIdx;
	triStart = _vertexFace[vertex];
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
//			normVert = _pnI[normVert].posNormIdx;
		if(normVert==vertex)
			break;
	}
	assert(j<3);
	adjs = &(_adjs[(trNum<<1)+trNum]);
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
		adjs = &(_adjs[(n.triangle<<1)+n.triangle]);
		j = adj&0x00000003;
		n.vertex = tnow[(j+2)%3];
		neighbors.push_back(n);
		adj = adjs[(j+2)%3];
	}
}

void staticTriangle::getVertexNormal(unsigned int vertex, float (&normal)[3])
{
	normal[0]=0.0f; normal[1]=0.0f; normal[2]=0.0f;
	if(!_adjacenciesComputed)
		findAdjacentTriangles();
	GLfloat *v,*nv,last[3],now[3];
	v = vertexCoordinate(vertex);
	std::vector<neighborNode> nei;
	getNeighbors(vertex,nei);
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

bool staticTriangle::localPick(const float *lineStart, float *lineDirection, float (&position)[3], int &triangle, float &param)
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

int staticTriangle::linePick(const float *lineStart, float *lineDirection, std::vector<float> &positions, std::vector<int> &triangles, std::vector<float> &params)
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

void staticTriangle::getTriangulatedSurface(std::vector<float> &vertices, std::vector<int> &tris)
{ // gets simple surface with no graphics duplication of vertices
	tris.clear();
	int posMax=-1,i,j,k,n=(int)_tris.size();
	tris.assign(n,-1);
	for(i=0; i<n; i+=3) {
		for(j=0; j<3; ++j) {
			k = _pnI[_tris[i+j]];
			tris[i+j] = k;
			if(posMax<k)
				posMax=k;
		}
	}
	vertices.clear();
	vertices.assign((posMax+1)*3,-1.0f);
	n = (int)_pnI.size();
	for(i=0; i<n; ++i)	{
		j = _pnI[i]*3;
		vertices[j] = _xyz1[i<<2];
		vertices[j+1] = _xyz1[(i<<2)+1];
		vertices[j+2] = _xyz1[(i<<2)+2];
	}
}

void staticTriangle::getSurfaceTriangles(std::vector<int> &triangles)
{
    triangles.clear();
    triangles.resize( _tris.size() );
    for( int i = 0; i < _tris.size(); i++){
        triangles[i] = _tris[i];       
    }
}

void staticTriangle::getSurfaceVertices(std::vector<float> &vertices)
{
    vertices.clear();
    vertices.resize( (_xyz1.size()>>2)*3 );
    for( int i = 0; i < _xyz1.size()>>2; i++){
        vertices[i*3+0] = _xyz1[i<<2];
        vertices[i*3+1] = _xyz1[(i<<2)+1];
        vertices[i*3+2] = _xyz1[(i<<2)+2];       
    }
}

void staticTriangle::getSurfaceNormals(std::vector<float> &normals)
{
    normals.clear();
    normals.resize((_xyz1.size()>>2)*3 );
	for(size_t i=0; i<_xyz1.size()>>2; ++i)
        getVertexNormal((unsigned int)i,(float (&)[3])(normals[i*3]));
}

void staticTriangle::getSurfaceUVs(std::vector<float> &uv)
{
    uv.clear();
    int n = _TexCoords.size();
    uv.reserve(n);
    for( int i = 0; i<n; i++)
        uv.push_back((float)_TexCoords[i]);
}


void staticTriangle::getClosestBarycentric(int triangle, float (&xyz)[3], float (&uv)[2])
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

void staticTriangle::getBarycentricPosition(int triangle, float (&uv)[2], float (&xyz)[3])
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

void staticTriangle::getBarycentricNormal(int triangle, float (&uv)[2], float (&normal)[3])
{	// for barycentric uv in triangle returns normal
	GLfloat p[3],q[3];
	GLuint *t = triVerts(triangle);
	getVertexNormal(t[0],p);
	getVertexNormal(t[1],q);
	Vec3f u,v,r;
	u.set(q[0]-p[0],q[1]-p[1],q[2]-p[2]);
	getVertexNormal(t[2],q);
	v.set(q[0]-p[0],q[1]-p[1],q[2]-p[2]);
	r = u*uv[0] + v*uv[1];
	normal[0]=r.x()+p[0]; normal[1]=r.y()+p[1]; normal[2]=r.z()+p[2];
}

void staticTriangle::setVertexCoordinate(long vertex, const float (&newCoord)[3])
{	// get base coordinate then update all the copies
//	vertex = _pnI[vertex].posIdx;
	GLfloat *v = &(_xyz1[vertex<<2]);
	v[0]=newCoord[0]; v[1]=newCoord[1]; v[2]=newCoord[2];
}

void staticTriangle::setVertexTexture(long vertex, const float (&newTex)[2])
{
	GLfloat *tx = &(_TexCoords[vertex<<1]);
	tx[0]=newTex[0]; tx[1]=newTex[1];
}

void staticTriangle::getVertexCoordinate(unsigned int vertex, float (&xyz)[3])
{	// type safe version
	GLfloat *v = &_xyz1[vertex<<2];
	xyz[0]=v[0]; xyz[1]=v[1]; xyz[2]=v[2];
}

bool staticTriangle::setTextureFileCreateProgram(const char *topTxFile, const char *topNrmFile)
{  // must be set first
	//if(_glslProgram>0)
	//	return false;
	//std::vector<std::string> att;
	//att.assign(4,std::string());
	//att[0] = "vVertex";
	//att[1] = "vNormal";
	//att[2] = "vTangent";
	//att[3] = "vTexture";
	//bool ret = _wxg->addCustomSceneNode(this,topTxFile,topNrmFile,staticVertexShader,staticFragmentShader,att);
	//if(!ret)
	//	return ret;
	//_wxg->getLightsShaders()->useGlslProgram(_glslProgram);  // must be current program. This routine sets other uniforms.
	//if(!_bufferObjects[0])
	//	glGenBuffers(5, _bufferObjects);
//	else	{  // COURT - is this necessary or can I just realloc in same objects?
//			glDeleteBuffers(5, _bufferObjects);
//			glGenBuffers(5, _bufferObjects);	}
	// Vertex data
	//glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[0]);	// VERTEX_DATA
	//glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_xyz1.size(), &(_xyz1[0]), GL_STATIC_DRAW);
	// Normal data
	//glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[1]);	// NORMAL_DATA
	//glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_Normals.size(), &(_Normals[0]), GL_STATIC_DRAW);
	// Tangent data
	//glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[2]);	// TANGENT_DATA
	//glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_Tangents.size(), &(_Tangents[0]), GL_STATIC_DRAW);
    // Texture coordinates
	//glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[3]);	// TEXTURE_DATA
	//glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_TexCoords.size(), &(_TexCoords[0]), GL_STATIC_DRAW);
    //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bufferObjects[4]);	// INDEX_DATA
	//glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*_tris.size(), &(_tris[0]), GL_STATIC_DRAW);
	//if(!_vertexArrayBufferObject)
	//	glGenVertexArrays(1,&_vertexArrayBufferObject);
	// now make vertex array
	//glBindVertexArray(_vertexArrayBufferObject);
    // Position data
    //glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[0]);
	//glEnableVertexAttribArray(0);
	//glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
    // Normal data
    //glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[1]);
	//glEnableVertexAttribArray(1);
	//glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
	// Tangent data
    //glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[2]);
	//glEnableVertexAttribArray(2);
	//glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, 0);
    // Texture coordinates
	//glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[3]);	// TEXTURE_DATA
	//glEnableVertexAttribArray(3);
	//glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, 0, 0);
    // Indexes
    //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bufferObjects[4]);	// INDEX_DATA
	// never unbind a GL_ARRAY_BUFFER or GL_ELEMENT_ARRAY_BUFFER inside a vertexArrayBuffer
	//glBindVertexArray(0);
	return true;
}

void staticTriangle::sendNewColoredTopology()
{
/*
	if(!_bufferObjects[0])
		glGenBuffers(3, _bufferObjects);
	// Vertex data
	glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[0]);	// VERTEX_DATA
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_xyz1.size(), &(_xyz1[0]), GL_STATIC_DRAW);
	// Normal data
	glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[1]);	// NORMAL_DATA
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_Normals.size(), &(_Normals[0]), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bufferObjects[2]);	// INDEX_DATA
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*_tris.size(), &(_tris[0]), GL_STATIC_DRAW);
	if(!_vertexArrayBufferObject)
		glGenVertexArrays(1,&_vertexArrayBufferObject);
	// now make vertex array
	glBindVertexArray(_vertexArrayBufferObject);
    // Position data
    glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[0]);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
    // Normal data
    glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[1]);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
    // Indexes
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bufferObjects[2]);	// INDEX_DATA
	// never unbind a GL_ARRAY_BUFFER or GL_ELEMENT_ARRAY_BUFFER inside a vertexArrayBuffer
	glBindVertexArray(0);
*/
}

void staticTriangle::computeNormalsTangents()
{
	unsigned int i,j,k,n;
	GLfloat *gv[3],*gtx[3];
	_Normals.assign((_TexCoords.size()>>1)*3,0.0f);
	if(_computeTangents)
		_Tangents.assign((_TexCoords.size()>>1)*3,0.0f);
	n = (unsigned int)_tris.size();
	GLuint *tr;
	Vec3f tanV,nrmV,dXyz[2];
	float frct,dUv[2][2],d2;
	for(i=0; i<n; i+=3) {
		tr = &_tris[i];
		for(j=0; j<3; ++j) {
			nrmV[j]=0.0f;
			tanV[j]=0.0f;
			gv[j]=&_xyz1[tr[j]<<2];
			if(_computeTangents)
				gtx[j]=&_TexCoords[tr[j]<<1];
		}
		for(j=0; j<3; ++j) {
			dXyz[0][j] = gv[1][j]-gv[0][j];
			dXyz[1][j] = gv[2][j]-gv[0][j];
			if(_computeTangents && j<2) {
				dUv[0][j] = gtx[1][j]-gtx[0][j];
				dUv[1][j] = gtx[2][j]-gtx[0][j]; }
		}
		nrmV = dXyz[0]^dXyz[1];
		for(j=0; j<3; ++j) {
			_Normals[tr[j]*3] += nrmV[0];
			_Normals[tr[j]*3+1] += nrmV[1];
			_Normals[tr[j]*3+2] += nrmV[2]; }
		if(!_computeTangents)
			continue;
		if(fabs(dUv[0][0])>fabs(dUv[1][0])) {
			j=0; k=1; }
		else {
			j=1; k=0; }
		frct = dUv[k][0]/dUv[j][0];
		dUv[j][1] *= frct;
		dXyz[j] *= frct;
		dUv[k][1] -= dUv[j][1];
		dXyz[k] -= dXyz[j];
		if(dUv[k][1]>0)
			tanV = dXyz[k];
		else
			tanV = -dXyz[k];
		for(j=0; j<3; ++j) {
			_Tangents[tr[j]*3] += tanV[0];
			_Tangents[tr[j]*3+1] += tanV[1];
			_Tangents[tr[j]*3+2] += tanV[2]; }
	}
	n = (int)(_xyz1.size()>>2)*3;
	for(i=0; i<n; i+=3) {
		d2 = _Normals[i]*_Normals[i] + _Normals[i+1]*_Normals[i+1] + _Normals[i+2]*_Normals[i+2];
		if(d2<1e-16f) {
			_Normals[i]=0.0f; _Normals[i+1]=0.0f; _Normals[i+2]=1.0f; }
		else {
			d2 = 1.0f/sqrt(d2);
			_Normals[i]*=d2; _Normals[i+1]*=d2; _Normals[i+2]*=d2;}
		if(!_computeTangents)
			continue;
		Vec3f *np=(Vec3f*)&_Normals[i],*bp=(Vec3f*)&_Tangents[i];
		tanV = *bp^*np;
		d2 = tanV.length2();
		if(d2<1e-16f) {
			_Tangents[i]=0.0f; _Tangents[i+1]=-1.0f; _Tangents[i+2]=0.0f; }
		else {
			d2 = 1.0f/sqrt(d2);
			tanV *= d2;
			_Tangents[i]=tanV[0]; _Tangents[i+1]=tanV[1]; _Tangents[i+2]=tanV[2];}
	}
}

