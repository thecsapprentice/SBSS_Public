// File: sceneNode.h
// Author: Court Cutting, MD
// Date: February 7, 2012
// Purpose: Basic data and methods every sceneNode must have.

#ifndef __SCENENODE_H__
#define __SCENENODE_H__

// Bring in OpenGL 
// Windows
#ifdef WIN32
#include <windows.h>		// Must have for Windows platform builds
#ifndef GLEW_STATIC
#define GLEW_STATIC
#endif
#include "GL/glew.h"			// OpenGL Extension "autoloader"
//#include <gl\gl.h>			// Microsoft OpenGL headers (version 1.1 by themselves)
#endif

// Mac OS X
#ifdef __APPLE__
#include <TargetConditionals.h>
#if TARGET_OS_IPHONE | TARGET_IPHONE_SIMULATOR
#include <OpenGLES/ES2/gl.h>
#include <OpenGLES/ES2/glext.h>
#define OPENGL_ES
#else
#include <GL/glew.h>
#include <OpenGL/gl.h>		// Apple OpenGL haders (version depends on OS X SDK version)
#endif
#endif

// Linux
#ifdef linux
#define GLEW_STATIC
#include <GL/glew.h>
#endif

#include <GL/glew.h>
#include "GLmatrices.h"
#include <string>

class sceneNode
{
public:
	enum nodeType{CONE, SPHERE, CYLINDER, TRIANGLES, LINES, UVW_TRIANGLES};
	void setType(nodeType type) {_type=type;}
	inline nodeType getType() {return _type;}
	virtual void computeLocalBounds() {}
//	virtual void getLocalBounds(GLfloat (&localCenter)[3], GLfloat &Radius) {Radius=0; localCenter[0]=0;}
	void getLocalBounds(GLfloat (&localCenter)[3], GLfloat &Radius) {
		localCenter[0]=_localCenter[0]; localCenter[1]=_localCenter[1]; localCenter[2]=_localCenter[2];
		Radius = _radius;
	}

	void getBounds(GLfloat (&center)[3], GLfloat &radius, bool recomputeAll)
	{
		if(recomputeAll || !_boundsComputed)
			computeLocalBounds();
		transformVector3(_localCenter,_pat,center);
		GLfloat r = _pat[0]*_radius; r*=r;	// assume no unequal scaling. In some apps this will be a bug. Done for speed.
		r += _pat[4]*_radius*_pat[4]*_radius;
		r += _pat[8]*_radius*_pat[8]*_radius;
		radius = sqrt(r);
	}

	float getRadius() {if(!_boundsComputed) computeLocalBounds(); return _radius;}

	virtual void draw(void) {}
	void setName(const char *name) {_name=name;}
	void getName(std::string &name) {name=_name;}
	void set2DtextureBufferNumber(GLuint texBufNum) {_2DtextureBuffer=texBufNum;}
	void set2DnormalBufferNumber(GLuint nrmBufNum) {_2DnormalBuffer=nrmBufNum;}
	void setGlslProgramNumber(GLuint progNum) {_glslProgram=progNum;}
	inline GLuint getGlslProgramNumber() {return _glslProgram;}
	inline GLfloat* getModelViewMatrix() {return _pat;}
	inline bool coloredNotTextured() {return _coloredNotTextured; }
	GLfloat* getColor() {return _color;}
	inline void setColor(float (&color)[4]) {_color[0]=color[0]; _color[1]=color[1]; _color[2]=color[2]; _color[3]=color[3];}
	sceneNode() {
		_boundsComputed=false;
		loadIdentity4x4(_pat);
		_coloredNotTextured=true;
		_glslProgram = 0;
		_2DtextureBuffer = 0;
		_2DnormalBuffer = 0;
		_color[0]=1.0f; _color[1]=1.0f; _color[2]=1.0f; _color[3]=1.0f;
	}

	~sceneNode() {}

protected:
	bool _coloredNotTextured;
	nodeType _type;
	std::string _name;
	GLfloat _pat[16];	// GL matrix for local position attitude transform
	bool _boundsComputed;
	GLfloat _localCenter[3],_radius,_color[4];
	GLuint _glslProgram;
	GLuint _2DtextureBuffer;
	GLuint _2DnormalBuffer;

};

#endif	// __SCENENODE_H__
