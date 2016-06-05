// staticTriangle.h
// Author: Court Cutting
// Date: 5/20/2014
// Purpose: Triangle object management class for display of static triangle objects.
//       Copyright 2014 - All rights reserved.

#ifndef __staticTriangle__
#define __staticTriangle__

#include <vector>
#include <string>
#include <list>
#include <map>
#include <set>
#include "sceneNode.h"
#include <sstream>

// forward declarations
class wxGraphics;

class staticTriangle : public sceneNode
{
public:
	void sendNewColoredTopology();
	bool setTextureFileCreateProgram(const char *topTxFile, const char *topNrmFile);  // must be set first before next 2 routines can be called
	void reserveVertices(unsigned int nVertices) {_TexCoords.reserve(nVertices<<1); _xyz1.reserve(nVertices<<2); _Normals.reserve(nVertices*3); _Tangents.reserve(nVertices*3); }
	void reserveTriangles(unsigned int nTriangles) {_tris.reserve(nTriangles*3); }
	long addTriangle(unsigned long (&vertices)[3]);
	void deleteTriangle(long triangleNumber);	// does not remove any stranded vertices
	long addCoordNormTexVertices(int numberToAdd=1);
	struct neighborNode{
		long	vertex;
		long triangle;
	};
	bool findAdjacentTriangles();
	void getNeighbors(unsigned int vertex, std::vector<neighborNode> &neighbors);
	int readObjFile(const char *data, bool dataOnlyNoGraphics=false);
//	bool writeObjFile(const char *fileName);
	void computeLocalBounds();
	void computeNormalsTangents();
	void draw(void);
	inline GLfloat* vertexCoordinate(long vertex) {return (GLfloat*)&_xyz1[vertex<<2];}	// next 2 calls have no error checking for speed. Careful.
	inline GLfloat* vertexTexture(long vertex) {return (GLfloat*)&_TexCoords[vertex<<1];}
	void setVertexCoordinate(long vertex, const float (&newCoord)[3]);
	void setVertexTexture(long vertex, const float (&newTex)[2]);
	void getVertexCoordinate(unsigned int vertex, float (&xyz)[3]);
	void getVertexNormal(unsigned int vertex, float (&normal)[3]);
	int getVertexTriangle(int vertexNumber);	// gets triangle vertex is a member of
	bool getTriangleVertices(unsigned int triangle, int (&vertices)[3]);
	// consider using getTriangleVertices() instead of next routine for type safety
	inline GLuint* triVerts(long triangleNumber) {return &(_tris[(triangleNumber<<1)+triangleNumber]);}	// this routine is dangerous. It returns GLuint*, NOT unsigned int*. Some compilers don't autoconvert this type correctly
	int triangleNumber() {return (int)_tris.size()/3;}
	int vertexNumber() {return (int)_xyz1.size()>>2;}
	// next if true gives first intersection with line at param, returning position and triangle hit
	bool localPick(const float *lineStart, float *lineDirection, float (&position)[3], int &triangle, float &param);
	int linePick(const float *lineStart, float *lineDirection, std::vector<float> &positions, std::vector<int> &triangles, std::vector<float> &params);
	int getClosestVertex(float (&position)[3], int triangle=-1);
	void getTriangulatedSurface(std::vector<float> &vertices, std::vector<int> &tris);	// returns a simple triangulated surface
    void getSurfaceTriangles(std::vector<int> &triangles);
    void getSurfaceVertices(std::vector<float> &vertices);
    void getSurfaceNormals(std::vector<float> &normals);
    void getSurfaceUVs(std::vector<float> &uv);

	GLfloat* getPositionArray(int &numVerts) {numVerts=(int)(_xyz1.size()>>2); return &_xyz1[0]; }	// CAREFUL - direct GL array[4]. vertexCoordinate() safer.
	void getClosestBarycentric(int triangle, float (&xyz)[3], float (&uv)[2]);	// for position xyz return barycentric uv in triangle
	void getBarycentricPosition(int triangle, float (&uv)[2], float (&xyz)[3]);
	void getBarycentricNormal(int triangle, float (&uv)[2], float (&normal)[3]);
	// if above input triangle<0 search all triangles, else search only side containing input triangle.
	inline bool texturedNotColored() {return _textured; }
//	void makeLineList(std::vector<GLuint> &lines);
	const 	std::vector<GLfloat>& getPositionsArray() {return _xyz1;}
	void setWxGraphics(wxGraphics *wxg) {_wxg=wxg;}  // must be set before using this class for graphics output
	void setTexturedNotColored(bool texturedNotColored);
	staticTriangle(bool texturedNotColored=true);
	~staticTriangle();

private:
	static const GLchar *staticVertexShader;
	static const GLchar *staticFragmentShader;
	wxGraphics *_wxg;
	bool 	_textured,_computeTangents;
	std::vector<GLuint> _tris;	// each set of 3 are indices into the vertex arrays that make a triangle 
	std::vector<GLfloat> _xyz1;	// Array of unique positions of objects. For glsl, 4 elements with v[3]=1.0f
	std::vector<GLfloat> _Normals,_Tangents; // 3 compnents each. Not needed for javascript
	std::vector<long> _pnI; // index into the original position array loaded in the .obj file
	std::vector<GLfloat> _TexCoords;    // Array of texture coordinates of dimension 2
	unsigned long _vertexNumber,_triangleNumber,_uniqueNormalNumber;
	bool _adjacenciesComputed;
	std::vector<unsigned long> _adjs;	// low 2 bits are the edge number of the adjacent triangle.
			// If low 2 bits==3 and high order 30 bits==0, there is no adjacent triangle.
			// high 30 bits are the triangle number which must be bit shifted down 2
	struct edge	{
		unsigned long reversed : 1;
		unsigned long vtxMin : 31;
		unsigned long vtxMax;
		unsigned int adjCode;	// adjacency code for this edge
	};
	struct edgeTest : public std::binary_function<edge,edge,bool>
	{		// must be a less than operator
		bool operator()(const edge &e1,const edge &e2) const
		{ 
			if( e1.vtxMin < e2.vtxMin )
				return true;
			else if( e1.vtxMin > e2.vtxMin )
				return false;
			else	// vtxMin values equal
				return (e1.vtxMax < e2.vtxMax);
		}
	};
	std::vector<unsigned long> _vertexFace;	// unique normal version
	struct cnt{
		int v[3];	// 0 is coord index, 1 is texture index, 2 is normal index
//		cnt& operator=(const cnt& lp)	{this->coord=lp.coord; this->norm=lp.norm; this->tex=lp.tex; return *this;}
	};
	struct cntTest : public std::binary_function<cnt,cnt,bool>
	{		// must be a less than operator
		bool operator()(const cnt &c1,const cnt &c2) const
		{ 
			if( c1.v[0] < c2.v[0] )
				return true;
			else if( c1.v[0] > c2.v[0] )
				return false;
			else if( c1.v[1] < c2.v[1] )
				return true;
			else if( c1.v[1] > c2.v[1] )
				return false;
			else
				return (c1.v[2] < c2.v[2]);
		}
	};
	typedef std::map<cnt,int,cntTest> CNTMAP;
	struct posTri{
		int triangle;
		float v[3];
	};
	bool parseNextInputFileLine(std::stringstream *infile, std::string &unparsedLine, std::vector<std::string> &parsedLine);
	void makeVertexToTriangleMap();
	unsigned long _2DtextureBufferNumber;

	GLuint _nMaxIndexes;         // Maximum workspace
	GLuint _nNumIndexes;         // Number of indexes currently used
	GLuint _nNumVerts;           // Number of vertices actually used

	GLuint _bufferObjects[5];
	GLuint _vertexArrayBufferObject;

	friend class uvwConvert;  // COURT - nuke this later
};

#endif    // __staticTriangle__
