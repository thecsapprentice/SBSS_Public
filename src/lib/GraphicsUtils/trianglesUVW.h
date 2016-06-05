//////////////////////////////////////////////////////////
// File: trianglesUVW.h
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

#ifndef __TRIANGLES_UVW__
#define __TRIANGLES_UVW__

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

// forward declarations

class trianglesUVW
{
public:
	int readObjFile(const char *data);
	bool writeObjFile(const char *fileName);
	void getVertexCoordinate(unsigned int vertex, float (&xyz)[3]);
	bool getBarycentricProjection(const int triangle, const float x, const float y, const float z, float (&uv)[2]);
	void getBarycentricPosition(const int triangle, const float (&uv)[2], float (&xyz)[3]);
	void getBarycentricTexture(const int triangle, const float (&uv)[2], float (&uvw)[3]);
	int addVertices(int numberToAdd=1);  // Warning these two routines invalidate all pointers and iterators
	int addTriangle(int (&vertices)[3]);
	// ray inputs below are 3 element array pointers. Outputs triangles intersected and parameters along line.
	int rayIntersect(const float *rayStart, const float *rayDirection, std::vector<int> &triangles, std::vector<float> &params);
	struct neighborNode{
		long	vertex;
		long triangle;
	};
	bool findAdjacentTriangles(bool fullManifoldTest=false);    // builds adjacency array for rapid neighbor searches
	void getNeighbors(unsigned int vertex, std::vector<neighborNode> &neighbors);
	inline float* vertexCoordinate(int vertex) {return (float*)&_xyz[(vertex<<1)+vertex];}	// next 3 calls have no error checking for speed. Careful.
	inline float* vertexTexture(int vertex) {return (float*)&_uvw[(vertex<<1)+vertex];}
	int* triangleVertices(unsigned int triangle) {return &_tris[(triangle<<1)+triangle];}
	void setVertexCoordinate(int vertex, const float (&newCoord)[3]);
	void setVertexTexture(long vertex, const float (&newTex)[3]);
	int getVertexTriangle(int vertexNumber){return _vertexFace[vertexNumber];}	// gets triangle vertex is a member of
	void getTriangleNormal(int triangle, float (&normal)[3], bool normalized=true);
	void getMeanVertexNormal(int vertex, float (&normal)[3]);
//	bool getTriangleVertices(unsigned int triangle, int (&vertices)[3]);
	inline int* triVerts(long triangleNumber) {return &(_tris[(triangleNumber<<1)+triangleNumber]);}
	inline int numberOfTriangles() { return (int)_tris.size()/3; }
	inline unsigned long* triAdjs(long triangleNumber) {return &(_adjs[(triangleNumber<<1)+triangleNumber]);}
	inline unsigned long* vertexTriangle(int vertex) {return &(_vertexFace[vertex]);} // careful if you don't know what you're doing

	int* getTriangleArray(int &numberOfTriangles);
	float* getPositionArray(int &numberOfVertices);
	float* getTextureArray(int &numberOfVertices);
	std::vector<int>* getTriangleArray() {return &_tris;}
	std::vector<float>* getPositionArray() {return &_xyz;}
	std::vector<float>* getTextureArray() {return &_uvw;}
	void getNearestHardEdge(float (&xyz)[3], int &triangle, int &edge, float &param, int searchLimit=2);
	bool localPick(const float *lineStart, float *lineDirection, float (&position)[3], int &triangle, float &param);
	int linePick(const float *lineStart, const float *lineDirection, std::vector<float> &positions, std::vector<int> &triangles, std::vector<float> &params);
	int splitTriangleEdge(int triangle, int edge, const float parameter);
	int addNewVertexInMidTriangle(int triangle, const float (&uvParameters)[2]);
	int isManifoldConsistent();  // return -1 if inconsistent and # of topological handles if consistent
	bool deleteEdge(int triangle, int edge, int vertexRemaining); // vertexRemaining==0 leaves first vertex of edge, 1 leaves second
	void cleanAndPack();  // warning - invalidates all triangle and vertex indices.
	void cleanAndPack(std::vector<int> &newVertexMap, std::vector<int> &newTriangleMap); // returns mapping of old indices to new
	bool topoCheck(); // checks current topology versus recomputed from scratch

	trianglesUVW(void);
	~trianglesUVW(void);

private:
    std::vector<int> _tris;	// each set of 3 are indices into the vertex arrays that make a triangle. -1 indicates deletion.
	std::vector<float> _xyz;    // 3 float per vertex position data.
	std::vector<float> _uvw;    // 3 float per vertex texture data
    bool _adjacenciesComputed;
    std::vector<unsigned long> _adjs;	// low 2 bits are the edge number of the adjacent triangle.
        // If low 2 bits==3 and high order 30 bits==0, there is no adjacent triangle.
        // high 30 bits are the triangle number which must be bit shifted down 2
	struct edge	{
		unsigned long reversed : 1;
		unsigned long vtxMin : 31;
		unsigned long matched : 1;  // check for non-manifold surface
		unsigned long vtxMax : 31;
		unsigned int adjCode;	// adjacency code for this edge
	};
//	struct edgeTest : public std::binary_function<edge,edge,bool>
//	{		// must be a less than operator
//		bool operator()(const edge &e1,const edge &e2) const
	struct edgeTest {		// must be a less than operator
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
	void makeVertexToTriangleMap();
	std::vector<unsigned long> _vertexFace;
	bool parseNextInputFileLine(std::stringstream *infile, std::string &unparsedLine, std::vector<std::string> &parsedLine);
	void interpolateEdgeTextures(int triangle, int edge, int newVert, float param);

	friend class uvwConvert;	// this class benefits tremendously from interior access
	friend class incision;	// interior access vital here

};

#endif  // __TRIANGLES_UVW__
