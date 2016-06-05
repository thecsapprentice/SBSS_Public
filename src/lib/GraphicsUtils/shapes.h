// File: shapes.h
// Author: Court Cutting
// Date: 3/3/12
// Purpose: Manages untextured cones, cylinders and spheres for graphics.

#ifndef __SHAPES_H__
#define __SHAPES_H__

#include "sceneNode.h"
#include <list>

// forward declaration
class wxGraphics;

class shape: public sceneNode
{
public:
	void computeLocalBounds();
	shape(int type) { setType((sceneNode::nodeType) type);}
	~shape() {}
};

class shapes
{
public:
	enum shapeType{CONE, SPHERE, CYLINDER};
	shape* addShape(shapeType type, char *name);
	bool deleteShape(shape* shapePtr);
	void drawCylinder();
	bool pickCylinder(const float *lineStart, float *lineDirection, float (&position)[3], float &param);
	void drawSphere();
	bool pickSphere(const float *lineStart, float *lineDirection, float (&position)[3], float &param);
	void drawCone();
	bool pickCone(const float *lineStart, float *lineDirection, float (&position)[3], float &param);
	void setWxGraphics(wxGraphics *wgg) {_wgg=wgg;}
	void clear();    // empties graphics card of shape objects
	shapes();
	~shapes();

private:
	wxGraphics *_wgg;
	void getOrCreateConeGraphic(); // creates cone point at origin with base at z=1.0 with 1.0 diameter
	void getOrCreateSphereGraphic(); // creates sphere at origin with diameter 1.0
	void getOrCreateCylinderGraphic(); // creates cylinder at origin with diameter 1.0, from z=-1 to z=1
	std::list<shape> _shapes;
	static GLuint _coneBufferObjects[2];
	static GLuint _coneVertexArrayBufferObject;
	static GLuint _sphereBufferObjects[3];
	static GLuint _sphereVertexArrayBufferObject;
	static GLuint _cylinderBufferObjects[3];
	static GLuint _cylinderVertexArrayBufferObject;

};

#endif	// __SHAPES_H__
