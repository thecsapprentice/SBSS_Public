// File: fences.h
// Author: Court Cutting, MD
// Date: June 7, 2012
// Purpose: User interface for creating a fence on a glslTriangle object.
//     This will be used to specify a desired incision line.

#ifndef __FENCES_H__
#define __FENCES_H__

#include <list>
#include <set>
#include "GraphicsUtils/wxGraphics.h"

// forward declarations

class fence
{
public:
	void addPost(trianglesUVW *tri, int triangle, float (&xyz)[3], float (&normal)[3], bool connectToNearestEdge);
	int getPostNumber() {return (int)_posts.size();}
	int getPostData(std::vector<float> &positions,std::vector<float> &normals,std::vector<int> &triangles,std::vector<float> &uv,bool &edgeStart,bool &edgeEnd);
	void deleteLastPost();
	inline staticTriangle* getStaticTrianglePtr() {return _slt;}
	void setWxGraphics(wxGraphics *wxg) {_wxg=wxg; _shapes=wxg->getShapes(); _glm=wxg->getGLmatrices(); if(fence::_fenceSize<10000.0f) _initialized=true; }
	void setFenceSize(float size) {fence::_fenceSize=size; if(_wxg!=NULL) _initialized=true;}
	void clear();	// deletes this fence
	bool isInitialized()	{return _initialized;}
	fence();
	~fence();

private:
	struct fencePost{
		bool connectToEdge;
		float xyz[3];
		float nrm[3];
		int triangle;
		float uv[2];
		shape *cylinderShape;
	};
	std::list<fencePost> _posts;
	wxGraphics *_wxg;
	shapes *_shapes;
	GLmatrices *_glm;
	staticTriangle *_slt;
	bool _initialized;
	static float _fenceSize;
	static GLfloat _selectedColor[4],_unselectedColor[4];
};

#endif // #ifndef __FENCES_H__
