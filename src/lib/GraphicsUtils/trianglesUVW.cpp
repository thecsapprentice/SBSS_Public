//////////////////////////////////////////////////////////
// File: trianglesUVW.cpp
// Author: Court Cutting, MD
// Date: 4/9/2014
// Purpose: Triangle storage class using only xyz position
//    and uvw texture data for vertices.  Normals are not used
//    and texture seams are not allowed thereby making all
//    vertices unique. An auxilliary vertexNormalUVW class
//    with its included vertex shaders provides vertex normal
//    doubling and tripling for graphics and user purposes.
//    Other auxilliary classes (e.g. skinFragmentUVW) provide
//    a tissue specific fragment shader to procedurally texture
//    the model for graphics purposes.
//////////////////////////////////////////////////////////

#include <assert.h>
#include <fstream>
#include <string.h>
#include <exception>
#include <stdexcept>
#include <functional>
#include <algorithm>
#include <sstream>
#include "trianglesUVW.h"
#include "glslTriangle.h"
#include "math3d.h"
#include "Vec3f.h"

int trianglesUVW::readObjFile(const char *data)
{ // returned error codes: 0=no error, 1=can't open file, 2=non-triangle primitive,
	// 3=bad 3D vertex line, 4=bad 3D texture line, 5=bad uvw face line, 6=exceeds 0x3fffffff vertex limit
    std::stringstream fin;
    fin << std::string( data );
    fin.seekg( std::ios::beg );
	_xyz.clear();
	_uvw.clear();
	_tris.clear();
    std::string unparsedLine;
	std::vector<std::string> parsedLine;
	unsigned int vertexNumber=0;
	unsigned int triangleNumber=0;
	int i,j,k,l;
    long long line_number = 0;
    
	while(parseNextInputFileLine(&fin,unparsedLine,parsedLine))
	{
        line_number++;
		if(parsedLine.empty())
			continue;
		if(parsedLine[0]=="v")
		{
			if(parsedLine.size()!=4){
                for(i=0; i<parsedLine.size(); ++i)
                    std::cout << parsedLine[i] << std::endl;
                throw std::runtime_error( "Cannot read object. Line "+std::to_string(line_number)+
                                          ": Vertex entry malformed, expected 4 tokens but got " + std::to_string((long long)(parsedLine.size())) );
            }
			for(i=1; i<4; ++i)
				_xyz.push_back((float)atof(parsedLine[i].c_str()));
		}
		else if(parsedLine[0]=="vt")
		{
			if(parsedLine.size()!=4)
                throw std::runtime_error( "Cannot read object. Line "+std::to_string(line_number)+
                                          ": Texture Coordinate entry malformed.");
			_uvw.push_back((float)atof(parsedLine[1].c_str()));
			_uvw.push_back((float)atof(parsedLine[2].c_str()));
			_uvw.push_back((float)atof(parsedLine[3].c_str()));
		}
		else if(parsedLine[0]=="f")
		{	// always in vertexPosition/vertexTexture format. If vP/vT may skip normal. If vP//vN, texture is skipped.
			int numVerts = (int)parsedLine.size();
			if(numVerts != 4)
                throw std::runtime_error( "Cannot read object. Line "+std::to_string(line_number)+
                                          ": Non-Triangle Primative.");

			for(i=1; i<numVerts; ++i)
			{
                std::string facevertex_part;
                std::stringstream facevertex_token;
                facevertex_token.str( parsedLine[i] );
                std::vector<int> facevertex_components;
                
                // Determine formating (Count number of '/' characters in token)
                int num_parts = std::count_if( parsedLine[i].begin(), parsedLine[i].end(),
                                               std::bind1st( std::equal_to<char>(),'/'));
                // Process Token
                switch( num_parts ){
                case 0: // Vertex only.
                    throw std::runtime_error( "Cannot read object. Line "+std::to_string(line_number)+
                                              ": No Texture for Face Component.");
                    break;
                case 1: // Vertex and Textures
                    // Vertex Part
                    std::getline(facevertex_token, facevertex_part, '/' );
                    facevertex_components.push_back( std::stoi( facevertex_part )-1 );
                    // Texture Part
                    std::getline(facevertex_token, facevertex_part);
                    facevertex_components.push_back( std::stoi( facevertex_part )-1 );                    
                    break;
                case 2: // Vertex, Textures, Normals
                    // Vertex Part
                    std::getline(facevertex_token, facevertex_part, '/' );
                    facevertex_components.push_back( std::stoi( facevertex_part )-1 );
                    // Texture Part
                    std::getline(facevertex_token, facevertex_part, '/' );
                    facevertex_components.push_back( std::stoi( facevertex_part )-1 );
                    // Normal Part
                    std::getline(facevertex_token, facevertex_part);
                    facevertex_components.push_back( std::stoi( facevertex_part )-1 );
                    break;
                }

                if( facevertex_components[0] != facevertex_components[1] )
                    throw std::runtime_error( "Cannot read object. Line "+std::to_string(line_number)+
                                              ": Bad UVW Face.");
                _tris.push_back( facevertex_components[0] );
               
			}
			++triangleNumber;
		}
		else
			continue;
	}
	if(	vertexNumber>0x3fffffff) {
        throw std::runtime_error( "Cannot read object. Too many vertices.");
	}
	// trim excess capacity
	std::vector<int>(_tris).swap(_tris);
	std::vector<float>(_xyz).swap(_xyz);
	std::vector<float>(_uvw).swap(_uvw);
	_adjacenciesComputed = false;
	return 0;
}

bool trianglesUVW::parseNextInputFileLine(std::stringstream *infile, std::string &unparsedLine, std::vector<std::string> &parsedLine)
{
    std::string token;
    std::stringstream ssline;

	if(infile->eof())
	{
		return false;
	}

    std::getline( *infile, unparsedLine );
    parsedLine.clear();
    if( unparsedLine == "" )
        return true;

    ssline.str( unparsedLine );
    ssline >> token;
    do{
        parsedLine.push_back( token );
        ssline >> token;
    }while( !ssline.eof() );

    return true;
}

bool trianglesUVW::writeObjFile(const char *fileName)
{
	std::string title(fileName);
	if(title.rfind(".obj")>title.size())
		title.append(".obj");
	std::ofstream fin(title.c_str());
    if(!fin.is_open())
        return false;
	char s[400];
	std::string line;
	for(size_t i=0; i<_xyz.size(); i+=3)	{
		sprintf(s,"v %f %f %f\n",_xyz[i],_xyz[i+1],_xyz[i+2]);
		line.assign(s);
		fin.write(line.c_str(),line.size());	}
	for(size_t i=0; i<_uvw.size(); i+=3)	{
		sprintf(s,"vt %f %f %f\n",_uvw[i],_uvw[i+1],_uvw[i+2]);
		line.assign(s);
		fin.write(line.c_str(),line.size());	}
	for(size_t i=0; i<_tris.size(); i+=3)	{
		if(_tris[i]<0)
			continue;
		sprintf(s,"f %d/%d %d/%d %d/%d\n",_tris[i]+1,_tris[i]+1,_tris[i+1]+1,_tris[i+1]+1,_tris[i+2]+1,_tris[i+2]+1);
		line.assign(s);
		fin.write(line.c_str(),line.size());	}
	fin.close();
	return true;
}

void trianglesUVW::getVertexCoordinate(unsigned int vertex, float (&xyz)[3])
{	// type safe version
	float *v = &_xyz[vertex*3];
	xyz[0]=v[0]; xyz[1]=v[1]; xyz[2]=v[2];
}

bool trianglesUVW::getBarycentricProjection(const int triangle, const float x, const float y, const float z, float (&uv)[2])
{	// for position xyz return barycentric uv projection into triangle
	float *p,*q;
	int *t = &_tris[(triangle<<1)+triangle];
	p = vertexCoordinate(t[0]);
	q = vertexCoordinate(t[1]);
	Vec3f u,v,xmp;
	u.set(q[0]-p[0],q[1]-p[1],q[2]-p[2]);
	q = vertexCoordinate(t[2]);
	v.set(q[0]-p[0],q[1]-p[1],q[2]-p[2]);
	xmp.set(x-p[0],y-p[1],z-p[2]);
	float a,b,c,d;
	a=u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
	b=u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
	c=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
	if(fabs(d=b*b-a*c)<1e-16f)	{	// degenerate triangle
		uv[0]=0.0; uv[1]=0.0f;
		return false; }
	uv[1] = ((u*b - v*a)*xmp)/d;
	uv[0] = (xmp*u - uv[1]*b)/a;
	return true;
}

void trianglesUVW::getBarycentricTexture(const int triangle, const float (&uv)[2], float (&uvw)[3])
{
	int *tr = &_tris[(triangle<<1)+triangle];
	float p=1.0f-uv[0]-uv[1],*t0=vertexTexture(tr[0]),*t1=vertexTexture(tr[1]),*t2=vertexTexture(tr[2]);
	for(int i=0; i<3; ++i)
		uvw[i] = t0[i]*p + uv[0]*t1[i] + uv[1]*t2[i];
}

void trianglesUVW::getBarycentricPosition(const int triangle, const float (&uv)[2], float (&xyz)[3])
{	// for barycentric uv in triangle returns position in xyz
	float *p,*q;
	int *t = &_tris[(triangle<<1)+triangle];
	p = vertexCoordinate(t[0]);
	q = vertexCoordinate(t[1]);
	Vec3f u,v,r;
	u.set(q[0]-p[0],q[1]-p[1],q[2]-p[2]);
	q = vertexCoordinate(t[2]);
	v.set(q[0]-p[0],q[1]-p[1],q[2]-p[2]);
	r = u*uv[0] + v*uv[1];
	xyz[0]=r.x()+p[0]; xyz[1]=r.y()+p[1]; xyz[2]=r.z()+p[2];
}

int trianglesUVW::rayIntersect(const float *rayStart, const float *rayDirection, std::vector<int> &triangles, std::vector<float> &params)
{ // lineStart and lineDirection are both 3 element vectors
	std::map<float,int> hitMap;
	std::map<float,int>::iterator hit;
	int j,bigAxis=0;
	if(fabs(rayDirection[1])>fabs(rayDirection[0]))
		bigAxis =1;
	if(fabs(rayDirection[2])>fabs(rayDirection[bigAxis]))
		bigAxis =2;
	int *tr;
	unsigned int i,tNum=(unsigned int)_tris.size()/3;
	float *v[3];
	float minimax[6],tMin,tMax,rMin,rMax;
	for(i=0; i<tNum; ++i)	{
		tr = &_tris[(i<<1)+i];
		if(*tr<0)
			continue;
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
		tMin = (minimax[bigAxis<<1]-rayStart[bigAxis])/rayDirection[bigAxis];
		tMax = (minimax[(bigAxis<<1)+1]-rayStart[bigAxis])/rayDirection[bigAxis];
		for(j=0; j<3; ++j) {
			if(j==bigAxis)
				continue;
			rMin=rayStart[j]+rayDirection[j]*tMin;
			rMax=rayStart[j]+rayDirection[j]*tMax;
			if(rMin<minimax[j<<1] && rMax<minimax[j<<1])
				break;
			if(rMin>minimax[(j<<1)+1] && rMax>minimax[(j<<1)+1])
				break;
		}
		if(j<3)
			continue;
		// now look for triangle intersection
		float b[3],m[9],r[3];
		for(j=0; j<3; ++j) {
			b[j]=rayStart[j]-v[0][j];
			m[(j<<1)+j] = -rayDirection[j];
			m[(j<<1)+j+1] = v[1][j]-v[0][j];
			m[(j<<1)+j+2] = v[2][j]-v[0][j];
		}
		if(!m3dSolveLinearSystem3D(m,b,r))
			continue;
		if(r[1]<0.0f || r[2]<0.0f || r[1]+r[2]>1.0f)
			continue;
		hitMap.insert(std::make_pair(r[0],i));
	}
	triangles.clear();	params.clear();
	triangles.reserve(hitMap.size());
	params.reserve(hitMap.size());
	for(hit=hitMap.begin(); hit!=hitMap.end(); ++hit)	{
		params.push_back(hit->first);
		triangles.push_back(hit->second);
	}
	return (int)hitMap.size();
}

bool trianglesUVW::findAdjacentTriangles(bool fullManifoldTest)
{	// computes all the adjacent triangles from raw triangle input
	// returns false if non-manifold surface is input
	if(_adjacenciesComputed)
		return true;
	typedef std::set<edge,edgeTest> edgeSet;
	typedef edgeSet::iterator edgeIt;
	std::pair <edgeIt,bool> P;
	edge E;
	edgeIt ei;
	edgeSet M;
	M.clear();
	int *tnow;
	unsigned long *adjNow;
	unsigned long i,j,tcode,numtris=(unsigned long)_tris.size();
	_adjs.clear();
	_adjs.assign(numtris,0x00000003);
	for(i=0; i<numtris; i+=3)
	{
		if(_tris[i]<0)	// signals a deleted triangle
			continue;
		tnow = &_tris[i];
		adjNow = &(_adjs[i]);
		for(j=0; j<3; j++) {
			if(adjNow[j]!=0x00000003)	// adjacency already computed
				continue;
			int tmp;
			if((tmp=tnow[(j+1)%3])<tnow[j]) {
				E.vtxMin = tmp;
				E.vtxMax = tnow[j];
				E.reversed = 1;	}
			else {
				E.vtxMin = tnow[j];
				E.vtxMax = tmp;
				E.reversed = 0;	}
			E.matched = 0;
			E.adjCode = ((i/3)<<2) + j;
			P = M.insert(E);
			// if P.second is true, no match so edge inserted
			if(P.second==false)	// edge match found
			{
				if(fullManifoldTest && P.first->matched)
					return false;  // third edge with these coordinates so not manifold
				if(P.first->reversed==E.reversed && E.vtxMin!=E.vtxMax)
					return false;  // triangle ordering error
				tcode = P.first->adjCode;
				adjNow[j] = tcode;
				_adjs[(tcode>>2)*3+(tcode&0x00000003)] = E.adjCode;
				M.erase(P.first);
				if(fullManifoldTest) {  // manifold test for tripled edge
					E.matched = 1;
					M.insert(E); }
			}
			else
				adjNow[j] = 0x00000003;
		}
	}
	M.clear();
	makeVertexToTriangleMap();
	_adjacenciesComputed = true;
	return true;
}

void trianglesUVW::makeVertexToTriangleMap()
{
	int i,j,numtris=(int)_tris.size();
	unsigned int vnow;
	int *tnow;
	// provide each vertex with a face it is a member of
	_vertexFace.clear();
	_vertexFace.assign(_xyz.size()/3,0x80000000);	// initially deleted
	for(i=0; i<numtris; i+=3)
	{
		tnow = &(_tris[i]);
		if(tnow==NULL || tnow[0]<0)	// signals a deleted triangle
			continue;
		for(j=0; j<3; j++)
		{
			vnow = tnow[j];
			if(_vertexFace[vnow]&0x40000000)
				continue;	// vertex first on free edge, don't change
			_vertexFace[vnow] = i/3;
			if(_adjs[i+j]==0x00000003)	// vertex first on free edge, lock it for easy neighbor find
				_vertexFace[vnow] |= 0x40000000;
		}
	}
}

void trianglesUVW::getNeighbors(unsigned int vertex, std::vector<neighborNode> &neighbors)
{
	unsigned long trNum,adj,triStart;
	triStart = _vertexFace[vertex];
	neighbors.clear();
	if(triStart&0x80000000)	// unconnected vertex
		return;
	trNum = triStart&0x3fffffff;
	unsigned long *adjs;
	int *tnow = &(_tris[(trNum<<1)+trNum]);
	assert(tnow[0]>-1);	// deleted triangle
	int j;
	for(j=0; j<3; ++j)
		if(tnow[j]==vertex)
			break;
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
		tnow = &(_tris[(n.triangle<<1)+n.triangle]);
		adjs = &(_adjs[(n.triangle<<1)+n.triangle]);
		j = adj&0x00000003;
		n.vertex = tnow[(j+2)%3];
		neighbors.push_back(n);
		adj = adjs[(j+2)%3];
	}
}

trianglesUVW::trianglesUVW(void) :_adjacenciesComputed(false)
{
}


trianglesUVW::~trianglesUVW(void)
{
}

int* trianglesUVW::getTriangleArray(int &numberOfTriangles)
{
	numberOfTriangles = (int)_tris.size()/3;
	return &_tris[0];
}

float* trianglesUVW::getPositionArray(int &numberOfVertices)
{
	numberOfVertices = (int)_xyz.size()/3;
	return &_xyz[0];
}

float* trianglesUVW::getTextureArray(int &numberOfVertices)
{
	numberOfVertices = (int)_uvw.size()/3;
	return &_uvw[0];
}

bool trianglesUVW::localPick(const float *lineStart, float *lineDirection, float (&position)[3], int &triangle, float &param)
{ // lineStart and lineDirection are both 3 element vectors
	bool picked=false;
	int j,bigAxis=0;
	if(fabs(lineDirection[1])>fabs(lineDirection[0]))
		bigAxis =1;
	if(fabs(lineDirection[2])>fabs(lineDirection[bigAxis]))
		bigAxis =2;
	int *tr,i;
	GLfloat *v[3];
	float t,vtx[3],minimax[6],tMin,tMax,rMin,rMax;
	for(i=0; i<(int)_tris.size(); i+=3)	{
		tr = &_tris[i];
		if(*tr<0)
			continue;
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
	triangle /=3;
	return picked;
}

int trianglesUVW::linePick(const float *lineStart, const float *lineDirection, std::vector<float> &positions, std::vector<int> &triangles, std::vector<float> &params)
{ // lineStart and lineDirection are both 3 element vectors
	struct posTri{
		int triangle;
		float v[3];
	};
	posTri pT;
	std::map<float,posTri> hitMap;
	std::map<float,posTri>::iterator hit;
	int j,bigAxis=0;
	if(fabs(lineDirection[1])>fabs(lineDirection[0]))
		bigAxis =1;
	if(fabs(lineDirection[2])>fabs(lineDirection[bigAxis]))
		bigAxis =2;
	int *tr,i;
	float *v[3];
	float t,vtx[3],minimax[6],tMin,tMax,rMin,rMax;
	for(i=0; i<_tris.size(); i+=3)	{
		tr = &_tris[i];
		if(*tr<0)
			continue;
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
		float pad = (tMax-tMin)*0.1f; // fix roundoff error problem
		for(j=0; j<3; ++j) {
			if(j==bigAxis)
				continue;
			rMin=lineStart[j]+lineDirection[j]*(tMin - pad);
			rMax=lineStart[j]+lineDirection[j]*(tMax + pad);
			if(rMin<minimax[j<<1] && rMax<minimax[j<<1])
				break;
			if(rMin>minimax[(j<<1)+1] && rMax>minimax[(j<<1)+1])
				break;
		}
		if(j<3)
			continue;
		if(m3dRayTriangleIntersection(lineStart,lineDirection,v[0],v[1],v[2],t,vtx)) {
			pT.triangle = i/3;
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

void trianglesUVW::getTriangleNormal(int triangle, float (&normal)[3], bool normalized)
{
	int *tr = &_tris[(triangle<<1)+triangle];
	Vec3f v0,v1,t0((float(&)[3])_xyz[(tr[0]<<1)+tr[0]]),t1((float(&)[3])_xyz[(tr[1]<<1)+tr[1]]),t2((float(&)[3])_xyz[(tr[2]<<1)+tr[2]]);
	v0 = t1-t0;
	v1 = t2-t0;
	t0 = v0^v1;
	if(normalized)
		t0.normalize();
	normal[0]=t0.x(); normal[1]=t0.y(); normal[2]=t0.z();
}

void trianglesUVW::getNearestHardEdge(float (&xyz)[3], int &triangle, int &edge, float &param, int searchLimit)
{	// Input xyz, then overwrites all 4 with the point on the nearest hard edge.
	// If searchLimit==0 only edges where w==0, searchLimit==1 only edges where w==1, and searchLimit==2 allows both.
	if(!_adjacenciesComputed)
		findAdjacentTriangles();
	int *tr,i,j,n=(int)_tris.size();
	Vec3f pt(xyz),e0,e1;
	float w[3],minDsq=1e30f;
	unsigned long *adj;
	for(i=0; i<n; i+=3)	{
		tr = &_tris[i];
		if(*tr<0)
			continue;
		adj = &_adjs[i];
		for(j=0; j<3; ++j)
			w[j] = _uvw[(tr[j]<<1)+tr[j]+2];
		for(j=0; j<3; ++j) {
			if(adj[j]==0x03) { // allow for possible unconnected free edge on top or bottom
				if(searchLimit>1 && (w[j]==0.0f || w[j]==1.0f))
					break;
				else if(w[j]==(float)searchLimit)
					break;
				else ;
			}
			if(w[j]!=w[(j+2)%3] && w[j]==w[(j+1)%3]) {
				if(searchLimit>1 && (w[j]==0.0f || w[j]==1.0f))
					break;
				else if(w[j]==(float)searchLimit)
					break;
				else ;
			}
		}
		if(j>2) continue;
		// candidate edge
		e0.set((float (&)[3])_xyz[tr[j]*3]);
		e1.set((float (&)[3])_xyz[tr[(j+1)%3]*3]);
		Vec3f v0=pt-e0,v1=pt-e1;
		e1 -= e0;
		assert(e1*e1>0.0f);
		float d=v0*v0,t=(v0*e1)/(e1*e1);
		if(t<0.0f ){
			if(minDsq>d) {
				minDsq=d;
				xyz[0] = e0._v[0];
				xyz[1] = e0._v[1];
				xyz[2] = e0._v[2];
				triangle = i;
				edge = j;
				param = 0.0; }
		}
		else if(t>1.0f){
			if(minDsq>(d=v1*v1)) {
				minDsq=d;
				xyz[0] = (e0+e1)._v[0];
				xyz[1] = (e0+e1)._v[1];
				xyz[2] = (e0+e1)._v[2];
				triangle = i;
				edge = j;
				param = 1.0; }
		}
		else {
			v1 = e0 + e1*t;
			v0 = v1 - pt;
			if((d=v0*v0)<minDsq) {
				minDsq=d;
				xyz[0] = v1._v[0];
				xyz[1] = v1._v[1];
				xyz[2] = v1._v[2];
				triangle = i;
				edge = j;
				param = t;
			}
		}
	}
	// always want top or bottom triangle reference
	unsigned long adjV = _adjs[triangle+edge];
	if(adjV!=0x03) {
		triangle = adjV>>2;
		edge = adjV&0x03;
		param = 1.0f - param;
	}
	else
		triangle /= 3;
}

void trianglesUVW::interpolateEdgeTextures(int triangle, int edge, int newVert, float param)
{	// assumes triangle hasn't been changed yet by newVert
	int *trVerts = triVerts(triangle);
	float *txOut = vertexTexture(newVert),*txIn;
	txIn=vertexTexture(trVerts[edge]);
	txOut[0]=(1.0f-param)*txIn[0]; txOut[1]=(1.0f-param)*txIn[1]; txOut[2]=(1.0f-param)*txIn[2];
	txIn=vertexTexture(trVerts[(edge+1)%3]);
	txOut[0]+=param*txIn[0]; txOut[1]+=param*txIn[1]; txOut[2]+=param*txIn[2];
}

int trianglesUVW::splitTriangleEdge(int triangle, int edge, const float parameter)
{	// Splits a triangle along edge(0-2) by parameter(0-1).
	// Creates 1 or 2 new triangles and one new vertex. Returns new vertex number of the input triangle split.
	// Does fix adjacency array so getNeighbors() works and texture interpolated.
	// Definitely should call getTriangleAdjacencies() ASAP.
	assert(0.0f<=parameter && 1.0f>=parameter);
	int *trVerts = triVerts(triangle);
	if(*trVerts<0)
		return -1;
	int tn,newVert= addVertices();
	int ve=trVerts[edge],ve1=trVerts[(edge+1)%3];
	float gv[3];
	float *gvp=vertexCoordinate(ve);
	gv[0]=(1.0f-parameter)*gvp[0]; gv[1]=(1.0f-parameter)*gvp[1]; gv[2]=(1.0f-parameter)*gvp[2];
	gvp=vertexCoordinate(ve1);
	gv[0]+=parameter*gvp[0]; gv[1]+=parameter*gvp[1]; gv[2]+=parameter*gvp[2];
	setVertexCoordinate(newVert,gv);
	interpolateEdgeTextures(triangle,edge,newVert,parameter);
	unsigned long *trAdjs = &(_adjs[triangle+(triangle<<1)]);
	int v[3];
	if(trAdjs[edge]==0x00000003)	{
		v[1] = trVerts[(edge+1)%3];
		trVerts[(edge+1)%3] = newVert;
		v[0]=newVert;	v[2]=trVerts[(edge+2)%3];
		tn = addTriangle(v);	// invalidates old _tris and _adjs pointers
		_adjs[tn+(tn<<1)] = 0x00000003;
		unsigned long adjTE = _adjs[triangle*3+((edge+1)%3)];
		_adjs[tn+(tn<<1)+1] = adjTE;
		_adjs[(adjTE>>2)*3 + (adjTE&0x03)] = (tn<<2)+1;
		_adjs[tn+(tn<<1)+2] = (triangle<<2)+((edge+1)%3);
		_adjs[triangle+(triangle<<1)+((edge+1)%3)] = (tn<<2)+2;
		if((_vertexFace[v[1]]&0x3fffffff)==triangle)
			_vertexFace[v[1]] = (tn | 0x40000000);	// first vertex on a free edge
		_vertexFace[newVert] = (tn | 0x40000000);
		return newVert;
	}
	int ve2=trVerts[(edge+2)%3],va2,tra = trAdjs[edge]>>2;
	int ea = trAdjs[edge]&0x00000003;
	tn=addTriangle(v);	// invalidates old _tris and _adjs pointers
	int tna=addTriangle(v);	// no data yet
	// new vertex assignments
	trVerts = triVerts(triangle);
	ve2 = trVerts[(edge+2)%3];
	trVerts[(edge+1)%3] = newVert;
	trVerts = triVerts(tra);
	assert(ve1==trVerts[ea]);
	assert(ve==trVerts[(ea+1)%3]);
	va2 = trVerts[(ea+2)%3];
	trVerts[(ea+1)%3] = newVert;
	trVerts = triVerts(tn);
	trVerts[0] = ve2;
	trVerts[1] = newVert;
	trVerts[2] = ve1;
	trVerts = triVerts(tna);
	trVerts[0] = va2;
	trVerts[1] = newVert;
	trVerts[2] = ve;
	// new adj assignments
	unsigned long ae1,aa1;
	trAdjs = &(_adjs[triangle+(triangle<<1)]);
	trAdjs[edge] = (tna<<2)+1;
	ae1 = trAdjs[(edge+1)%3];
	trAdjs[(edge+1)%3] = tn<<2;
	trAdjs = &(_adjs[tra+(tra<<1)]);
	assert(trAdjs[ea] == (triangle<<2)+edge);
	trAdjs[ea] = (tn<<2)+1;
	aa1 = trAdjs[(ea+1)%3];
	if(aa1!=0x03)
		assert(_adjs[(aa1>>2)*3+(aa1&0x03)] == (tra<<2)+((ea+1)%3));
	trAdjs[(ea+1)%3] = tna<<2;
	trAdjs = &(_adjs[tn+(tn<<1)]);
	trAdjs[0] = (triangle<<2)+ ((edge+1)%3);
	trAdjs[1] = (tra<<2)+ea; // already done once
	trAdjs[2] = ae1;
	trAdjs = &(_adjs[tna+(tna<<1)]);
	trAdjs[0] = (tra<<2)+ ((ea+1)%3);
	trAdjs[1] = (triangle<<2)+edge;
	trAdjs[2] = aa1;
	if(ae1!=3)
		_adjs[(ae1>>2)*3+(ae1&0x00000003)] = (tn<<2)+2;
	if(aa1!=0x03)
		_adjs[(aa1>>2)*3+(aa1&0x00000003)] = (tna<<2)+2;
	// new vertexFace assignments
	_vertexFace[newVert] = triangle;
	if((_vertexFace[ve1]&0x3fffffff)==triangle) {
		if((_vertexFace[ve1]&0x40000000)>0)
			_vertexFace[ve1] = tn|0x40000000;
		else
			_vertexFace[ve1] = tn;
	}
	if((_vertexFace[ve]&0x3fffffff)==tra) { // va1
		if((_vertexFace[ve]&0x40000000)>0)
			_vertexFace[ve] = tna|0x40000000;
		else
			_vertexFace[ve] = tna;
	}
	return newVert;
}

int trianglesUVW::addNewVertexInMidTriangle(int triangle, const float (&uvParameters)[2])
{	// creates 2 new triangles and one new vertex. Returns new vertex number.
	// Input uvParameters[2] are parameters along vectors t1-t0 and t2-t0.
	// Does fix adjacency array so getNeighbors() works and texture&positions interpolated.
	// Definitely should getTriangleAdjacencies() ASAP.
	assert(uvParameters[0]>=0.0f && uvParameters[0]<=1.0f);
	assert(uvParameters[1]>=0.0f && uvParameters[1]<=1.0f);
	assert(uvParameters[0]+uvParameters[1]<=1.0001f);
	int *trVerts = triVerts(triangle);
	if(*trVerts<0)
		return -1;
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
	// now we know we will add a vertex and two triangles. These operations invalidate pointers due to possible reallocation.
	int v[3],v0=trVerts[0],v1=trVerts[1],oldVert=trVerts[2];
	int t1=addTriangle(v),t2=addTriangle(v),ret= addVertices();
	unsigned long *trAdjs = &_adjs[triangle+(triangle<<1)];
	unsigned long a1=trAdjs[1],a2=trAdjs[2];
	trVerts = triVerts(triangle);
	float p,*pv = vertexCoordinate(v0);
	p = 1.0f - uvParameters[0] - uvParameters[1];
	float vec[3] = {p*pv[0],p*pv[1],p*pv[2]};
	pv = vertexTexture(v0);
	float tx[3]={p*pv[0],p*pv[1],p*pv[2]};
	for(int i=0; i<2; ++i)	{
		pv = vertexTexture(trVerts[i+1]);
		tx[0] += pv[0]*uvParameters[i];
		tx[1] += pv[1]*uvParameters[i];
		tx[2] += pv[2]*uvParameters[i];
		pv = vertexCoordinate(trVerts[i+1]);
		vec[0]+=pv[0]*uvParameters[i]; vec[1]+=pv[1]*uvParameters[i]; vec[2]+=pv[2]*uvParameters[i];
	}
	setVertexCoordinate(ret,vec);
	setVertexTexture(ret,tx);
	// assign vertices
	trVerts[2] = ret;
	trVerts = triVerts(t1);
	trVerts[0]=ret; trVerts[1]=v1; trVerts[2]=oldVert; 
	trVerts = triVerts(t2);
	trVerts[0]=ret; trVerts[1]=oldVert; trVerts[2]=v0;
	// assign adjs
	_adjs[triangle+(triangle<<1)+1] = t1<<2;
	_adjs[triangle+(triangle<<1)+2] = (t2<<2)+2;
	_adjs[t1+(t1<<1)] = (triangle<<2)+1;
	_adjs[t1+(t1<<1)+1] = a1;
	_adjs[t1+(t1<<1)+2] = t2<<2;
	_adjs[t2+(t2<<1)] = (t1<<2)+2;
	_adjs[t2+(t2<<1)+1] = a2;
	_adjs[t2+(t2<<1)+2] = (triangle<<2)+2;
	if(a1!=0x00000003)
		_adjs[(a1>>2)*3+(a1&0x00000003)] = (t1<<2)+1;
	if(a2!=0x00000003)
		_adjs[(a2>>2)*3+(a2&0x00000003)] = (t2<<2)+1;
	// assign vertexFace
	if((_vertexFace[oldVert]&0x3fffffff)==triangle)	{
		_vertexFace[oldVert]=t2;
		if(a2==0x00000003)
			_vertexFace[oldVert] |= 0x40000000;
	}
	if((_vertexFace[v1]&0x3fffffff)==triangle)	{
		_vertexFace[v1]=t1;
		if(a1==0x00000003)
			_vertexFace[v1] |= 0x40000000;
	}
	_vertexFace[ret]=triangle;
	// 	_vertexFace[v0] can stay unchanged
	return ret;
}

int trianglesUVW::addTriangle(int (&vertices)[3])
{
	int retval = (int)_tris.size()/3;
	_tris.push_back(vertices[0]);
	_tris.push_back(vertices[1]);
	_tris.push_back(vertices[2]);
	if(!_adjs.empty()) {
		_adjs.push_back(0x03);
		_adjs.push_back(0x03);
		_adjs.push_back(0x03); }
	_adjacenciesComputed = false;
	return retval;
}

int trianglesUVW::addVertices(int numberToAdd)
{
	int retval = (int)_xyz.size()/3;
	for(unsigned int i=0; i<(unsigned)numberToAdd; ++i)
	{
		_xyz.push_back(0.0f);
		_xyz.push_back(0.0f);
		_xyz.push_back(0.0f);
		_uvw.push_back(0.0f);
		_uvw.push_back(0.0f);
		_uvw.push_back(0.0f);
		if(!_vertexFace.empty())
			_vertexFace.push_back(0x80000000);
	}
	_adjacenciesComputed = false;
	return retval;
}

void trianglesUVW::setVertexCoordinate(int vertex, const float (&newCoord)[3])
{
	float *v = &_xyz[(vertex<<1)+vertex];
	v[0]=newCoord[0];
	v[1]=newCoord[1];
	v[2]=newCoord[2];
}

void trianglesUVW::setVertexTexture(long vertex, const float (&newTex)[3])
{
	float *v = &_uvw[(vertex<<1)+vertex];
	v[0]=newTex[0];
	v[1]=newTex[1];
	v[2]=newTex[2];
}

void trianglesUVW::cleanAndPack()
{	// warning - invalidates all triangle and vertex indices.
	_vertexFace.clear();
	_vertexFace.assign(_xyz.size()/3,0x80000000);
	_adjs.clear();
	_adjs.reserve(_tris.size());
	int i,j,k=0,n=(int)_tris.size();
	for(i=0; i<n; i+=3)	{
		if(_tris[i]<0)	// deleted triangle
			continue;
		for(j=0; j<3; ++j)	{
			_adjs.push_back(_tris[i+j]);
			_vertexFace[_tris[i+j]] = 0;
		}
	}
	_tris.clear();
	_tris.assign(_adjs.begin(),_adjs.end());
	n = (int)_vertexFace.size();
	int bot=0,top,newV=0;
	for(i=0; i<n; ++i)	{
		if(_vertexFace[i])
			continue;
		_vertexFace[i] = newV++;
		top = (i<<1)+i;
		if(bot<top) {
			for(j=0; j<3; ++j) {
				_uvw[bot]=_uvw[top];
				_xyz[bot++]=_xyz[top++];
			}
		}
		else
			bot += 3;
	}
	_xyz.resize(bot);
	_uvw.resize(bot);
	n = (int)_tris.size();
	for(i=0; i<n; ++i)
		_tris[i] = _vertexFace[_tris[i]];
	_adjacenciesComputed=false;
	_adjs.clear();
	_vertexFace.clear();
}

void trianglesUVW::cleanAndPack(std::vector<int> &newVertexMap, std::vector<int> &newTriangleMap)
{	// warning - invalidates all triangle and vertex indices. Returns the new index numbers of the old vertices and triangles. -1 indicates deletion.
	std::vector<int> pp;
	newVertexMap.clear();
	newVertexMap.assign(_xyz.size()/3,-1);
	pp.assign(_xyz.size()/3,-1);
	int i,j,k,n=(int)_tris.size();
	for(i=0; i<n; i+=3)	{
		if(_tris[i]<0)	// deleted triangle
			continue;
		for(j=0; j<3; ++j)	{
			newVertexMap[_tris[i+j]] = 1;
			pp[_tris[i+j]] = 1;
		}
	}
	n=(int)pp.size();
	j = 0;
	for(i=0; i<n; ++i)	{
		if(pp[i]<0)
			continue;
		if(i>j)	{
			_xyz[(j<<1)+j]=_xyz[(i<<1)+i];
			_xyz[(j<<1)+j+1]=_xyz[(i<<1)+i+1];
			_xyz[(j<<1)+j+2]=_xyz[(i<<1)+i+2];
		}
		pp[i]=j;
		++j;
	}
	_xyz.resize((j<<1)+j);
	n=(int)newVertexMap.size();
	j = 0;
	for(i=0; i<n; ++i)	{
		if(newVertexMap[i]<0)
			continue;
		_uvw[(j<<1)+j]=_uvw[(i<<1)+i];
		_uvw[(j<<1)+j+1]=_uvw[(i<<1)+i+1];
		_uvw[(j<<1)+j+2]=_uvw[(i<<1)+i+2];
		newVertexMap[i]=j;
		++j;
	}
	_xyz.resize((j<<1)+j);
	_uvw.resize((j<<1)+j);
	n=(int)_tris.size();
	newTriangleMap.clear();
	newTriangleMap.assign(n/3,-1);
	k=0;
	for(i=0; i<n; i+=3)	{
		if(_tris[i]<0)	// deleted triangle
			continue;
		for(j=0; j<3; ++j)	{
			assert(newVertexMap[_tris[i+j]]>-1);
			_tris[k+j] = newVertexMap[_tris[i+j]];	}
		newTriangleMap[i/3] = k/3;
		k += 3;
	}
	_tris.resize(k);
	_adjacenciesComputed=false;
	_adjs.clear();
	_vertexFace.clear();
}

int trianglesUVW::isManifoldConsistent()
{	// Closed manifold surface topology checker.  Returns # of topological handles or -1 if inconsistent
	typedef std::set<edge,edgeTest> edgeSet;
	edge E;
	edgeSet M;
	M.clear();
	int tnow[3];
	int i,j,numtris=(int)_tris.size();
	std::vector<int> posVec;
	posVec.assign(_xyz.size()/3,0);
	for(i=0; i<numtris; i+=3)
	{
		if(_tris[i]<0)	// signals a deleted triangle
			continue;
		for(j=0; j<3; j++)
			tnow[j] = _tris[i+j];
		for(j=0; j<3; j++)
		{
			if(tnow[j]>=posVec.size())
				return -1;
			else
				posVec[tnow[j]] = 1;
			unsigned int tmp;
			if((tmp=(unsigned int)tnow[(j+1)%3])<(unsigned int)tnow[j]) {
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
			return -1;
	int handles2 = (int)((_tris.size()/3) + posVec.size() - M.size());
	if(handles2&0x0001)
		return -1;
	return handles2>>1;
}

void trianglesUVW::getMeanVertexNormal(int vertex, float (&normal)[3])
{
	if(!_adjacenciesComputed)	findAdjacentTriangles();
	std::vector<neighborNode> nei;
	getNeighbors(vertex,nei);
	if(nei.size()<2) {
		normal[0]=0.0f; normal[1]=0.0f; normal[2]=0.0f;
		return; }
	std::vector<neighborNode>::iterator nit=nei.begin();
	Vec3f last,now,p,mean(0.0f,0.0f,0.0f);
	getVertexCoordinate(vertex,p._v);
	if(nei.front().triangle<0) {
		getVertexCoordinate(nei.front().vertex,last._v);
		++nit; }
	else
		getVertexCoordinate(nei.back().vertex,last._v);
	last -= p;
	while(nit!=nei.end()) {
		getVertexCoordinate(nit->vertex,now._v);
		now -= p;
		mean += now^last;
		last = now;
		++nit;
	}
	mean.normalize();
	normal[0]=mean._v[0]; normal[1]=mean._v[1]; normal[2]=mean._v[2];
}

bool trianglesUVW::deleteEdge(int triangle, int edge, int vertexRemaining)
{
	int *tr=&_tris[triangle*3];
	int ve=tr[edge],ve1=tr[(edge+1)%3];
	int vStay,vGo,atri;
	if(vertexRemaining>0) {
		vStay=ve1; vGo=ve; }
	else {
		vGo=ve1; vStay=ve; }
	unsigned long *vf0=&_vertexFace[ve],*vf1=&_vertexFace[ve1],*avf2=NULL,*vf2=&_vertexFace[tr[(edge+2)%3]];
	unsigned long aae1,aae2,ae,ae1,ae2,*aadj=NULL,*adj = &_adjs[triangle*3];
	ae = adj[edge];
	if(ae!=0x03) {
		atri = ae>>2;
		aadj = &_adjs[(ae>>2)*3];
		aae1 = aadj[((ae&3)+1)%3];
		aae2 = aadj[((ae&3)+2)%3];
		avf2=&_vertexFace[_tris[atri*3+(((ae&3)+2)%3)]];
	}
	ae1 = adj[(edge+1)%3];
	ae2 = adj[(edge+2)%3];
	// next section prevents creation of an illegal surface
	if(ae!=3 && *vf0&0x40000000 && *vf1&0x40000000) { // both on free edge in an isthmus
		// only if one of the triangles on either side of the edge isn't connected elsewhere is this permissable
		if(ae1==3 && ae2==3) {
			_adjs[atri*3+(ae&3)] = 3;
			assert(*vf2==(triangle|0x40000000));
			*vf2 = 0x80000000; // stranded
			adj[0]=3; adj[1]=3; adj[2]=3;
			tr[0] = -1;
			*vf1 = atri|0x40000000;
			// whatever *vf0 was it can stay
			return deleteEdge(atri,ae&3,vertexRemaining==1?0:1);
		}
		else {
			assert(aadj!=NULL && avf2!=NULL);
			if(aae1==3 && aae2==3) {
				_adjs[triangle*3+edge] = 3;
				assert(*avf2==(atri|0x40000000));
				*avf2 = 0x80000000; // stranded
				aadj[0]=3; aadj[1]=3; aadj[2]=3;
				_tris[atri*3] = -1;
				*vf0 = triangle|0x40000000;
				// whatever *vf1 was it can stay
				return deleteEdge(triangle,edge,vertexRemaining);
			}
			else
				return false; // disallow. Would create illegal surface
		}
	}
	std::vector<neighborNode> nei;
	std::vector<neighborNode>::iterator nit;
	getNeighbors(vGo,nei);
	assert(!nei.empty());
	int *goTr;
	for(nit=nei.begin(); nit!=nei.end(); ++nit) {
		if(nit->triangle<0)
			continue;
		goTr = triVerts(nit->triangle);
		if(*goTr<0) assert(false);
		int i;
		for(i=0; i<3; ++i) {
			if(goTr[i]==vGo) {
				goTr[i]=vStay;
				break; }
		}
		assert(i<3);
	}
	if(ae==0x03) {
		assert(*vf0&0x40000000 && *vf1&0x40000000);
		if((*vf1&0x3fffffff)==triangle) { // already known to be on edge
			assert(ae1==0x03);
			if(ae2==0x03) { // single stranded triangle
				*vf0=0x80000000; *vf1=0x80000000; *vf2=0x80000000; // strand vertices
				adj[0]=0x03; adj[1]=0x03; adj[2]=0x03; *tr=-1;
				return true;
			}
			else
				_vertexFace[vStay] = (ae2>>2)|0x40000000;
		}
		else
			_vertexFace[vStay] = *vf1;
		if((*vf2&0x3fffffff)==triangle) {
			if(*vf2&0x40000000) {
				assert(ae1!=0x03); // would have been eliminated by illegal check
				*vf2 = (ae1>>2)|0x40000000; } 
			else {
				assert(ae2!=0x03);
				*vf2 = ae2>>2;
			}
		}
		if(ae1!=0x03)   _adjs[(ae1>>2)*3+(ae1&0x03)] = ae2;
		if(ae2!=0x03)	_adjs[(ae2>>2)*3+(ae2&0x03)] = ae1;
		adj[0]=0x03; adj[1]=0x03; adj[1]=0x03;
		tr[0] = -1;
		_vertexFace[vGo]=0x80000000; // strand vertex. Do last as may have been using its pointer
		return true;
	}
	// triangles on both sides of edge
	tr[0]=-1;
	adj[0]=0x03; adj[1]=0x03; adj[2]=0x03;
	_tris[atri*3]=-1;
	aadj[0]=0x03; aadj[1]=0x03; aadj[2]=0x03;
	if(ae2!=0x03) _adjs[(ae2>>2)*3+(ae2&03)]=ae1;
	if(ae1!=0x03) _adjs[(ae1>>2)*3+(ae1&03)]=ae2;
	if(aae2!=0x03) _adjs[(aae2>>2)*3+(aae2&03)]=aae1;
	if(aae1!=0x03) _adjs[(aae1>>2)*3+(aae1&03)]=aae2;
	// both ve and ve1 not on a free edge
	if(*vf0&0x40000000) {
		_vertexFace[vStay] = *vf0;
		if((*vf0&0x3fffffff)==atri) {
			assert(aae2!=3); // this case should already have been eliminated by illegal check
			_vertexFace[vStay] = (aae2>>2)|0x40000000; } // else leave as is
	}
	else {
		_vertexFace[vStay] = *vf1;
		if(*vf1&0x40000000) {
			if((*vf1&0x3fffffff)==triangle){
				assert(ae1==3 && ae2!=3); // this case should already have been eliminated by illegal check
				_vertexFace[vStay] = (ae2>>2)|0x40000000; } // else leave as is
		}
		else // both completely surrounded
			_vertexFace[vStay] = ae2>>2;
	}
	if((*vf2&0x3fffffff)==triangle) {
		if(*vf2&0x40000000) {
			assert(ae1!=3 && ae2==3); // this case should already have been eliminated by illegal check
			*vf2 = (ae1>>2)|0x40000000; }
		else {
			if(ae1!=3)
				*vf2 = ae1>>2;
			else if(ae2!=0x03)
				*vf2 = ae1>>2;
			else
				assert(false);
		}
	}
	if((*avf2&0x3fffffff)==atri) {
		if(*avf2&0x40000000) {
			assert(aae2==3 && aae1!=3); // this case should already have been eliminated by illegal check
			*avf2 = (aae1>>2)|0x40000000; } // else leave as is
		else {
			if(aae1!=3)
				*avf2 = aae1>>2;
			else if(aae2!=3)
				*avf2 = aae1>>2;
			else
				assert(false);
		}
	}
	_vertexFace[vGo]=0x80000000; // strand vertex. Do last as may have been using its pointer
	return true;
}

bool trianglesUVW::topoCheck()
{
    std::vector<unsigned long> adjs;	// low 2 bits are the edge number of the adjacent triangle.
	adjs.assign(_adjs.begin(),_adjs.end());
	std::vector<unsigned long> vertexFace;
	vertexFace.assign(_vertexFace.begin(),_vertexFace.end());
	_adjacenciesComputed=false; // force computation
	findAdjacentTriangles();
	int i,n;
	n = (int)adjs.size();
	for(i=0; i<n; ++i) {
		if((i%3==0) && _tris[i]<0) {
			i+=2;
			continue; }
		if(adjs[i]!=_adjs[i])
			return false;
	}
	std::vector<trianglesUVW::neighborNode> nei;
	std::vector<trianglesUVW::neighborNode>::iterator nit;
	n = (int)vertexFace.size();
	for(i=0; i<n; ++i) {
		unsigned long vo,vf = vertexFace[i];
		vo = _vertexFace[i];
		if(vo>0xfffffffe || vf>0xfffffffe) {
			if(vo!=vf) // should be stranded vertex
				return false;
			else
				continue;
		}
		if((vo&0x40000000)>0) {
			if(vo!=vf)
				return false;
		}
		else {
			getNeighbors(i,nei);
			for(nit=nei.begin(); nit!=nei.end(); ++nit) {
				if(nit->triangle==vf)
					break;
			}
			if(nit==nei.end()) { // could be stranded vertex
				int j;
				for(j=0; j<n; ++j) {
					if((j%3==0) && _tris[j]<0) {
						j+=2;
						continue; }
					if(_tris[j]==i)
						break;
				}
				if(j<n) // not a stranded vertex
					return false;
			}
		}
	}
	return true;
}

