//////////////////////////////////////////////////////////
// File: skinGraphics.h
// Author: Court Cutting, MD
// Date: 4/15/2014
// Purpose: Takes as input a trianglesUVW class and makes it visible
//    on an openGL canvas.  It creates hard normal edges and textures
//    the model given an input texture file for the top skin. Skin side
//    texturing and possibly the bottom fat are done procedurally in
//    the fragment shader.
//////////////////////////////////////////////////////////

#ifndef __SKIN_GRAPHICS__
#define __SKIN_GRAPHICS__

#include <vector>
#include <list>
#include "sceneNode.h" // includes glew.h
#include "trianglesUVW.h"

// forward declarations
class wxGraphics;
class Vec3f;

class skinGraphics : public sceneNode
{
public:
	void draw(void);
	void computeLocalBounds();
	bool setTopTextureFilesCreateProgram(const char *topTextureFile, const char *topNormalFile);  // must be set first before next 2 routines can be called
	void setNewTopology();
	void updatePositionsAndNormals();
	inline trianglesUVW* getTrianglesUVW() {return & _tuvw;}  // gets the uvw textured data class
	void setWxGraphics(wxGraphics *wxg) {_wxg=wxg;}  // must be set before using this class for graphics output
	float getMeanEdgeTriangleLength() {return _meanTriangleEdgeLength; }  // only of value after an initial load of a trianglesUVW object
    void remapJavascriptVertexData(const std::vector<float> &xyz_in, std::vector<float> &xyz_out );
	void getJavascriptData(std::vector<int> &tris, std::vector<float> &xyz, std::vector<float> &uv);
	void getJavascriptPositions(std::vector<float> &xyz);
    void getJavascriptMap(std::vector<int>& map);
	skinGraphics(void);
	~skinGraphics(void);

private:
	// COURT - should there be a static graphics program number for this class?
	static const GLchar *skinVertexShader;
	static const GLchar *skinFragmentShader;
	trianglesUVW _tuvw;
	wxGraphics *_wxg;
	std::vector<GLuint> _tris;  // 0xffffffff signals a deleted triangle
	std::vector<GLfloat> _xyz1;
	std::vector<GLfloat> _uv;  // 3D converted to 2D textures for javascript libraries
	std::vector<GLfloat> _normals;  // these next two won't be transmitted for three.js
	std::vector<GLfloat> _tangents;
	unsigned int _origVertTop;
	std::vector<GLuint> _doubledVerts;
	bool _computeTangents;
	bool _getUvScale;
	float _uvScale;
	float _meanTriangleEdgeLength;
	void computeNormalsTangents();
	float getClosestTopV(int v, std::list<GLuint> &path);
	void recurseMiddle(int v, std::set<GLuint> &m, std::vector<GLuint> &newV);
	bool makeCrease(GLuint topV, bool firstCrease, GLfloat maxDist, GLfloat &maxV);
	GLuint _bufferObjects[5];  // these 2 lines for graphics card buffers
	GLuint _vertexArrayBufferObject;
	bool findAdjacentTriangles();
	void makeVertexToTriangleMap();
	struct neighborNode{
		long	vertex;
		long triangle;
	};
	void getNeighbors(unsigned int vertex, std::vector<neighborNode> &neighbors);
	std::vector<unsigned long> _adjs;	// low 2 bits are the edge number of the adjacent triangle.
			// If low 2 bits==3 and high order 30 bits==0, there is no adjacent triangle.
			// high 30 bits are the triangle number which must be bit shifted down 2
	std::vector<unsigned long> _vertexFace;	// unique normal version
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

};

#endif  // __SKIN_GRAPHICS__
