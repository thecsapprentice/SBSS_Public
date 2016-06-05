// glslTriangle.h
// Author: Court Cutting
// Date: 1/6/2012
// Purpose: Triangle object management class with full topological services which
//       makes maximal use of glsl, speeding graphics and unloading CPU maximally.
//       Copyright 2012 - All rights reserved.

#ifndef __glslTriangle__
#define __glslTriangle__

#include <vector>
#include <string>
#include <list>
#include <map>
#include <set>
#include "sceneNode.h"

class glslTriangle : public sceneNode
{
public:
	void reserveVertices(unsigned int nVertices) {_pnI.reserve(nVertices); _TexCoords.reserve(nVertices<<1); _Positions.reserve(nVertices<<2); }
	void reserveTriangles(unsigned int nTriangles) {_tris.reserve(nTriangles*3); }
	long addTriangle(unsigned long (&vertices)[3]);
	void deleteTriangle(long triangleNumber);	// does not remove any stranded vertices
	long addCoordNormTexVertices(int numberToAdd=1);
	long addNormTexVertex(int existingVertex);
	long addNormTexVertex2(int existingPositionIndex);	// adds new based only on an existing position array index. Avoid this unless you have to.
	long addTexVertex(int existingVertex);
	struct neighborNode{
		long	vertex;
		long triangle;
	};
	bool findAdjacentTriangles(bool uniqueNormalsNotCoords=true);
	void getNeighbors(unsigned int vertex, std::vector<neighborNode> &neighbors, bool uniqueNormalsNotCoords=true);
	int readObjFile(const char *fileName, bool dataOnlyNoGraphics=false);
	bool writeObjFile(const char *fileName);
	bool sendToGraphicsCard(); // sends current triangle object to glsl capable graphics card
	bool sendTopologyChangeToGraphicsCard(); // sends current triangle object to glsl capable graphics card
	bool makeObject();
	void setVerticesComputeNormals();	// pushes position array to graphics card and forces glsl normal recomputation
	void attach2Dtexture(unsigned long texBufferNumber) {_2DtextureBufferNumber=texBufferNumber;}
	void computeLocalBounds();
	void draw(void);
	void setDeformable(bool deformable=false) { _deformable=deformable;}
	inline GLfloat* vertexCoordinate(long vertex) {return (GLfloat*)&_Positions[_pnI[vertex].posIdx];}	// next 2 calls have no error checking for speed. Careful.
	inline GLfloat* vertexTexture(long vertex) {return (GLfloat*)&_TexCoords[vertex<<1];}
	void setVertexCoordinate(long vertex, const float (&newCoord)[3]);
	void setVertexTexture(long vertex, const float (&newTex)[2]);
	void getVertexCoordinate(unsigned int vertex, float (&xyz)[3]);
	void getVertexNormal(unsigned int vertex, float (&normal)[3], bool uniqueNormalsNotCoords=true);
	int getVertexTriangle(int vertexNumber);	// gets triangle vertex is a member of
	bool getTriangleVertices(unsigned int triangle, int (&vertices)[3]);
	void getTriangleAdjacencies(unsigned int triangle, int (&adjacentTriangles)[3], int (&adjTriEdges)[3], bool uniqueNormalsNotCoords=true);
	// consider using getTriangleVertices() instead of next routine for type safety
	inline GLuint* triVerts(long triangleNumber) {return &(_tris[(triangleNumber<<1)+triangleNumber]);}	// this routine is dangerous. It returns GLuint*, NOT unsigned int*. Some compilers don't autoconvert this type correctly
	int triangleNumber() {return (int)_tris.size()/3;}
	int vertexNumber() {return (int)_pnI.size();}
	// next if true gives first intersection with line at param, returning position and triangle hit
	bool localPick(const float *lineStart, float *lineDirection, float (&position)[3], int &triangle, float &param);
	int linePick(const float *lineStart, float *lineDirection, std::vector<float> &positions, std::vector<int> &triangles, std::vector<float> &params);
	int getClosestVertex(float (&position)[3], int triangle=-1);
	void getTriangulatedSurface(std::vector<float> &vertices, std::vector<int> &tris);	// returns a simple triangulated surface
    void getSurfaceTriangles(std::vector<int> &triangles);
    void getSurfaceVertices(std::vector<float> &vertices);
    void getSurfaceNormals(std::vector<float> &normals);
    void getSurfaceUVs(std::vector<float> &uv);

	GLfloat* getPositionArray(int &numVerts) {numVerts=(int)(_Positions.size()>>2); return &_Positions[0]; }	// CAREFUL - direct GL array[4]. vertexCoordinate() safer.
	int getPositionArrayIndex(int vertexNumber) {return _pnI[vertexNumber].posIdx; }
	int getPositionNormalIndex(int vertexNumber) {return _pnI[vertexNumber].posNormIdx; }
	void getClosestBarycentric(int triangle, float (&xyz)[3], float (&uv)[2]);	// for position xyz return barycentric uv in triangle
	void getBarycentricPosition(int triangle, float (&uv)[2], float (&xyz)[3]);
	void getBarycentricNormal(int triangle, float (&uv)[2], float (&normal)[3], bool uniqueNormalsNotCoords=true);
	void getNearestHardEdge(float (&xyz)[3], int &triangle, int &edge, float &param);	// Input first 2 args, then overwrites all 4 with the point on the nearest hard edge.
	// if above input triangle<0 search all triangles, else search only side containing input triangle.
	long splitTriangleEdge(int triangle, int edge, const float parameter);
	long addNewVertexInMidTriangle(int triangle, const float (&uvParameters)[2]);
	inline bool texturedNotColored() {return _textured; }
	void makeLineList(std::vector<GLuint> &lines);
	// warning - cleanAndPack() invalidates all triangle and vertex indices. Returns the new index numbers of the old vertices and triangles. -1 indicates deletion.
	void cleanAndPack(std::vector<int> &newVertexMap, std::vector<int> &newTriangleMap);
	const 	std::vector<GLfloat>& getPositionsArray() {return _Positions;}
	int isManifoldConsistent();	// topology checker
	glslTriangle(bool texturedNotColored=true, bool staticNotMorphing=false);
	~glslTriangle();

private:
	bool 	_textured,_static;
	std::vector<GLuint> _tris;	// each set of 3 are indices into the vertex arrays that make a triangle 
	std::vector<GLfloat> _Positions;        // Array of unique positions of objects. For glsl, 4 elements with v[3]=1.0f
	struct positionNormalIndex{
		long posIdx;	// index into the _Positions array already*4
		long posNormIdx;	// unique normal index
	};
	std::vector<positionNormalIndex> _pnI;
	std::vector<GLfloat> _TexCoords;    // Array of texture coordinates
	unsigned long _vertexNumber,_triangleNumber,_uniqueNormalNumber;
	void interpolateEdgeTextures(int triangle, int edge,  int newVert, float param);
	bool _adjacenciesComputed;
	bool _adjacenciesComputedP;
	std::vector<unsigned long> _adjs;	// low 2 bits are the edge number of the adjacent triangle.
			// If low 2 bits==3 and high order 30 bits==0, there is no adjacent triangle.
			// high 30 bits are the triangle number which must be bit shifted down 2
	std::vector<unsigned long> _adjsP;	// Unique position version of above
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
	void checkRing(std::set<edge,edgeTest> &M);
	std::vector<unsigned long> _vertexFace;	// unique normal version
	std::vector<unsigned long> _vertexFaceP;	// unique position version
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
	void recurseHardEdgeSearch(int triangle, std::set<int> &visited, std::vector<int> &hardEdges);
	bool parseNextInputFileLine(std::ifstream *infile, std::string &unparsedLine, std::vector<std::string> &parsedLine);
	void makeVertexToTriangleMap(bool uniqueNormalsNotCoords=true);
	unsigned long _2DtextureBufferNumber;
	void createGlslTopoArray(std::vector<GLuint> &topoArr);
	void createNormalMakerProgram();
	GLuint _NormalMakerProgram;	// this contains a vertex shader only
	bool _deformable;
	std::vector<GLuint>  _Indexes;        // Array of indexes
	std::vector<GLfloat> _Verts;        // Array of vertices
	std::vector<GLfloat> _Norms;        // Array of normals
	std::vector<GLfloat> _FeedbackNorms;        // Array of transform feedback normals
	std::vector<GLuint> _TopoIndex;
	std::vector<GLuint> _glslNeighbors;

	GLuint _nMaxIndexes;         // Maximum workspace
	GLuint _nNumIndexes;         // Number of indexes currently used
	GLuint _nNumVerts;           // Number of vertices actually used

	GLuint _bufferObjects[4];
	GLuint _vertexArrayBufferObject;
	GLuint _textureBufferObjects[2];
	GLuint _texBOBuffers[2];
	GLuint _transformFeedbackQuery;

	friend class incision;
	friend class uvwConvert;  // COURT - nuke this later
	friend class MainFrame;  // COURT - nuke this later
};

#endif    // __glslTriangle__
