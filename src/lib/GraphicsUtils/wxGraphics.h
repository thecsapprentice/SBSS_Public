// File: wxGraphics.h
// Author: Court Cutting, MD
// Date: January 21, 2012
// Purpose: Interface to wxGraphics library. This library uses wxWidgets as its
//	GUI library and the source of an openGL window.  The wxWidgets library has been
//	modified to use glew.h instead of gl.h and glu.h to create its openGL window.
//	This allows this library to use a more advanced version of openGL which this library
//	requires. This library does all the basic graphics operations for the next level up.
//	Copyright 2012 - Court Cutting - All rights reserved.

#ifndef __WX_GRAPHICS_H__
#define __WX_GRAPHICS_H__

#include "sceneNode.h"
#include "shapes.h"
#include "lines.h"
#include "textures.h"
#include "lightsShaders.h"
#include "trackball.h"
#include "GLmatrices.h"
#include <list>
#include <string>
#include "staticTriangle.h"
#include "skinGraphics.h"

class wxGraphics
{
public:
	bool addCustomSceneNode(sceneNode* sn, const char *texturePath, const char *normalPath, const GLchar *vertexShader, const GLchar *fragmentShader, std::vector<std::string> &attributes);
	void mouseButtonEvent(unsigned short screenX, unsigned short screenY, unsigned short button, bool dragging);
	bool pick(unsigned short x, unsigned short y, std::string &name, float (&position)[3], int &triangle, bool excludeShapes=false);
	bool initializeGraphics(std::string &errorMessage);
	staticTriangle* loadStaticObjFile(const char *filePath, const char *texturePath, const char *normalPath, bool texturedNotColored=true);
	staticTriangle* createEmptyStaticTriangle(bool texturedNotColored=true);
	void drawAll();	// draws everything in the scene
	void getSceneSphere(GLfloat (&center)[3], GLfloat &radius, bool recomputeAll=true);
	void frameScene(bool recomputeAll=true);	// sets view data for currently loaded scene
	void computeAndSetSceneRadius();	// alters scene radius without changing center or zoom
	void setViewport(unsigned short x, unsigned short y, unsigned short xSize, unsigned short ySize);
	inline shapes* getShapes() {return &_shapes;}
	inline lines* getLines() {return &_lines;}
	inline lightsShaders* getLightsShaders() {return &_ls;}
	inline textures* getTextures() {return &_texReader;}
	GLmatrices* getGLmatrices() {return &_glM;}
	void addSceneNode(sceneNode* sn) {_nodes.push_back(sn);}
	void deleteSceneNode(sceneNode* sn);
	sceneNode* getNodePtr(std::string &name);
	void clear();	// empties all graphics
	wxGraphics();
	~wxGraphics();

private:
	trackball _tBall;
	float _rotQuat[4];
	unsigned short _xSize,_ySize,_lastX,_lastY;
	std::list<staticTriangle> _sTris;
	std::list<sceneNode*> _nodes;
	shapes _shapes;
	lines _lines;
    GLfloat _xrot;
    GLfloat _yrot;
	bool _dragging;
	textures _texReader;
	GLmatrices _glM;
	GLuint  texture;
	lightsShaders _ls;
};

#endif	// __WX_GRAPHICS_H__
