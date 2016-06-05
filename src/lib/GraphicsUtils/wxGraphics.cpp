// File: wxGraphics.cpp
// Author: Court Cutting, MD
// Date: January 21, 2012
// Purpose: Interface to wxGraphics library. This library uses wxWidgets as its
//	GUI library and the source of an openGL window.  The wxWidgets library has been
//	modified to use glew.h instead of gl.h and glu.h to create its openGL window.
//	This allows this library to use a more advanced version of openGL which this library
//	requires. This library does all the basic graphics operations for the next level up.
//	Copyright 2012 - Court Cutting - All rights reserved.

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "math3d.h"
#include "wxGraphics.h"

bool wxGraphics::addCustomSceneNode(sceneNode* sn, const char *texturePath, const char *normalPath, const GLchar *vertexShader, const GLchar *fragmentShader, std::vector<std::string> &attributes)
{
	unsigned long txLoad,nrmLoad;
	if(*texturePath>='0' && *texturePath<='9')
		sscanf(texturePath,"%d",&txLoad);
	else
		txLoad = _texReader.loadTexture(texturePath);
	if(txLoad<1)
		return false;
	sn->set2DtextureBufferNumber(txLoad);
	if(normalPath!=NULL) {
		if(*normalPath>='0' && *normalPath<='9')
			sscanf(normalPath,"%d",&nrmLoad);
		else
			nrmLoad = _texReader.loadTexture(normalPath);
		if(nrmLoad<1)
			return false;
		sn->set2DnormalBufferNumber(nrmLoad);
	}
	GLuint progNum=sn->getGlslProgramNumber();
	if(!_ls.createCustomProgram(progNum,vertexShader,fragmentShader,attributes))
		return false;
	sn->setGlslProgramNumber(progNum);
	addSceneNode(sn);
	char tName[12];
	sprintf(tName,"Custom_%d",(int)progNum);
	sn->setName(tName);
	frameScene();
	return true;
}

staticTriangle* wxGraphics::loadStaticObjFile(const char *filePath, const char *texturePath, const char *normalPath, bool texturedNotColored)
{
	_sTris.push_back(staticTriangle(texturedNotColored));
	addSceneNode(&(_sTris.back()));
	staticTriangle *tr=&_sTris.back();
	tr->setWxGraphics(this);
	char tName[12];
	sprintf(tName,"Tri_%d",(int)_sTris.size());
	tr->setName(tName);
    int err = tr->readObjFile(filePath,true);
	if(err>0)	{
        printf( "Error code: %d\n", err );
		_nodes.pop_front();
		_sTris.pop_back();
		return NULL;
	}
	//tr->computeNormalsTangents();
	//tr->setTextureFileCreateProgram(texturePath,normalPath);
	//frameScene();
	//float c[3],r;
	//tr->getBounds(c,r,true);
	return tr;
}

staticTriangle* wxGraphics::createEmptyStaticTriangle(bool texturedNotColored)
{
	_sTris.push_back(staticTriangle(texturedNotColored));
	staticTriangle *tr=&(_sTris.back());
	addSceneNode(tr);
//	glslTriangle *tr = dynamic_cast<glslTriangle*>(_nodes.front());
	char tName[12];
	sprintf(tName,"Tri_%d",(int)_sTris.size());
	tr->setName(tName);
	if(!texturedNotColored) {
		GLuint pNum = _ls.getOrCreateColorProgram();
		tr->setGlslProgramNumber(pNum); }
	return tr;
}

bool wxGraphics::initializeGraphics(std::string &errorMessage)
{
    const GLubyte *verString = glGetString(GL_VERSION);
	// initialize glew only once
	GLenum err = glewInit();
	if (GLEW_OK != err)
	{
		/* Problem: glewInit failed, something is seriously wrong. */
		errorMessage = std::string((const char *)glewGetErrorString(err));
		return false;
	}
	errorMessage = "Status: Using GLEW ";
	errorMessage += std::string((const char *)glewGetString(GLEW_VERSION));
	_shapes.setWxGraphics(this);
	// Background
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f );
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	_tBall.computeQuat(_rotQuat,0.0f,0.0f,0.0f,0.0f);
	return true;
}

void wxGraphics::setViewport(unsigned short x, unsigned short y, unsigned short xSize, unsigned short ySize)
{
	assert(x==0 && y==0);
	_xSize=xSize;
	_ySize=ySize;
	    // It's up to the application code to update the OpenGL viewport settings.
    // In order to avoid extensive context switching, consider doing this in
    // OnPaint() rather than here, though.
    glViewport(0, 0, (GLint) xSize, (GLint) ySize);
	_glM.setView(0.7f,(float)xSize/ySize);
}

void wxGraphics::mouseButtonEvent(unsigned short screenX, unsigned short screenY, unsigned short button, bool dragging)
{
	if(dragging)	{
		_dragging=true;
		if( button<1)	{	// leftMouse
	        /* drag in progress, simulate trackball */
			float spin_quat[4];
			_tBall.computeQuat(spin_quat,
            (2.0f*_lastX - _xSize) / _xSize,
            (_ySize - 2.0f*_lastY) / _ySize,
            (2.0f*screenX - _xSize)    / _xSize,
            (_ySize - 2.0f*screenY)    / _ySize);
			_tBall.add_quats(spin_quat, _rotQuat, _rotQuat);
		}
		else if(button>1)	{	// rightMouse
			_glM.changeZoom((float)(_lastY-screenY)/_ySize);
		}
		else	{	// middleMouse
			_glM.setPan((float)(_lastX-screenX)/_xSize,(float)(screenY-_lastY)/_ySize);
		}
	}
	_lastX = screenX;
	_lastY = screenY;
}

wxGraphics::wxGraphics()
{
	_dragging=false;
	_lastX=0;
	_lastY=0;
	_ls.setGLmatrices(&_glM);
	_lines.setWxGraphics(this);
	_nodes.clear();
}

wxGraphics::~wxGraphics()
{
}

void wxGraphics::drawAll()
{
	// This is normally only necessary if there is more than one wxGLCanvas
    // or more than one wxGLContext in the application.
    //SetCurrent(*m_glRC);
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    GLfloat m[4][4];
    _tBall.build_rotmatrix( m,_rotQuat);
	// assumes perspective-frame matrix has already been set
	_glM.setFrameAndRotation(&m[0][0]);
	std::list<sceneNode*>::iterator nit;
	int type;
	GLuint currentProgram=0;
//	_ls.useTextureProgram();
	for(nit=_nodes.begin(); nit!=_nodes.end(); ++nit)	{ // textured TRIANGLES will always happen first
		if((*nit)->getGlslProgramNumber()!=currentProgram) {
			currentProgram=(*nit)->getGlslProgramNumber();
			_ls.useGlslProgram(currentProgram);
		}
		_ls.setModelMatrix((*nit)->getModelViewMatrix());	// must reset with new program

//		if(!(*nit)->coloredNotTextured())
//			(*nit)->draw();
//		else
//				_ls.useColorProgram();
//				_ls.setModelMatrix((*nit)->getModelViewMatrix());	// must reset with new program
			type = (*nit)->getType();
			if(type<sceneNode::SPHERE) {	// sceneNode::CONE
 				_ls.setColor((*nit)->getColor());
				_shapes.drawCone(); }
			else if(type<sceneNode::CYLINDER) {	// sceneNode::SPHERE
	 			_ls.setColor((*nit)->getColor());
				_shapes.drawSphere(); }
			else if(type<sceneNode::TRIANGLES) {	// type==sceneNode::CYLINDER
				_ls.setColor((*nit)->getColor());
				_shapes.drawCylinder(); }
			else if(type<sceneNode::LINES)	// type==sceneNode::TRIANGLES
				(*nit)->draw();
			else if(type<sceneNode::UVW_TRIANGLES)	{ // type==sceneNode::LINES
//				_ls.useLineProgram();
				_ls.setColor((*nit)->getColor());
				_ls.setModelMatrix((*nit)->getModelViewMatrix());	// must reset with new program
				_lines.drawLines();
//				_ls.useColorProgram();
			}
			else  // type==sceneNode::UVW_TRIANGLES
				(*nit)->draw();
	}
    glFlush(); // Not really necessary: buffer swapping below implies glFlush()
}

void wxGraphics::computeAndSetSceneRadius()
{ // does not change scene center
	float center[3],radius;
	_glM.getSceneCenter(center);
	radius = _glM.getSceneRadius();
	std::list<sceneNode*>::iterator nit;
	for(nit=_nodes.begin(); nit!=_nodes.end(); ++nit) {
		float c[3],r,d[3];
		(*nit)->getBounds(c,r,false);
		d[0]=center[0]-c[0]; d[1]=center[1]-c[1]; d[2]=center[2]-c[2];
		r += sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
		if(r>radius)
			radius = r;
	}
	_glM.setSceneRadius(radius);
}

void wxGraphics::getSceneSphere(GLfloat (&center)[3], GLfloat &radius, bool recomputeAll)	{
	std::list<sceneNode*>::iterator nit=_nodes.begin();
	(*nit)->getBounds(center,radius,recomputeAll);
	++nit;
	while(nit!=_nodes.end()) {
		float c[3],r,d[3],l;
		(*nit)->getBounds(c,r,recomputeAll);
		d[0]=center[0]-c[0]; d[1]=center[1]-c[1]; d[2]=center[2]-c[2];
		l=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
		if(l<1e-16)	{
			if(r>radius)
				radius=r;
		}
		else if((l+r)>radius) {	// expand scene sphere
			if(r-l>radius) {
				center[0]=c[0]; center[1]=c[1]; center[2]=c[2];
				radius = r;
			}
			else {
				float p = (l+r-radius)*0.5f/l;
				center[0]-=d[0]*p; center[1]-=d[1]*p; center[2]-=d[2]*p;
				radius += l+r;
				radius *= 0.5f;
			}
		}
		else
			;
		++nit;
	}
}

void wxGraphics::frameScene(bool recomputeAll)
{
    /*
	float r,c[3];
	getSceneSphere(c,r,recomputeAll);
	_glM.setSceneCenter(c);
	_glM.setSceneRadius(r);
	_glM.resetPerspective();
    */
}

bool wxGraphics::pick(unsigned short x, unsigned short y, std::string &name, float (&position)[3], int &triangle, bool excludeShapes)
{
	bool picked=false;
	triangle = -1;
	name = "";
	float sc[3],f[3],b[3],minParam=1e30f;
	float zCenter,height,verticalAngle,screenAspect;
	_glM.getCameraData(zCenter,height,verticalAngle,screenAspect);
	height = -zCenter*2.0f*tanf(verticalAngle*0.5f);
	sc[0] = (-0.5f+(float)x/_xSize)*height*screenAspect;
	sc[1] = (0.5f-(float)y/_ySize)*height;
	sc[2] = zCenter;
	GLfloat m[16],invM[16],*om;	// mfp[16],
	const GLfloat *fr=_glM.getFrameAndRotationMatrix();
	std::list<sceneNode*>::iterator nit;
	for(nit=_nodes.begin(); nit!=_nodes.end(); ++nit)	{
		if(excludeShapes && (*nit)->getType()<sceneNode::TRIANGLES)
			continue;
		om = (*nit)->getModelViewMatrix();
		for(int i=0; i<4; ++i)	{
			for(int j=0; j<4; ++j)	// model happens first, then frame-rotation
				m[(i<<2)+j] = fr[j]*om[i<<2] + fr[4+j]*om[(i<<2)+1] + fr[8+j]*om[(i<<2)+2] + fr[12+j]*om[(i<<2)+3];
		}
		invertMatrix4x4(m,invM);
		transformVector3(sc,invM,f);
		// start point is at origin
		b[0]=invM[12]; b[1]=invM[13]; b[2]=invM[14];
		f[0]-=b[0]; f[1]-=b[1]; f[2]-=b[2];
		if((*nit)->getType()==sceneNode::CONE) {
			if(_shapes.pickCone(b,f,position,minParam)) {
				picked=true;
				triangle = -1;
				(*nit)->getName(name);
				float vtx[]={position[0],position[1],position[2]};
				transformVector3(vtx,om,position);
			}
		}
		else if((*nit)->getType()==sceneNode::SPHERE) {
			if(_shapes.pickSphere(b,f,position,minParam)) {
				picked=true;
				triangle = -1;
				(*nit)->getName(name);
				float vtx[]={position[0],position[1],position[2]};
				transformVector3(vtx,om,position);
			}
		}
		else if((*nit)->getType()==sceneNode::CYLINDER) {
			if(_shapes.pickCylinder(b,f,position,minParam)) {
				picked=true;
				triangle = -1;
				(*nit)->getName(name);
				float vtx[]={position[0],position[1],position[2]};
				transformVector3(vtx,om,position);
			}
		}
		else if((*nit)->getType()==sceneNode::LINES)
			continue;
		else if((*nit)->getType()==sceneNode::UVW_TRIANGLES) {
			skinGraphics *sg = static_cast<skinGraphics*>(*nit);
			trianglesUVW *tr = sg->getTrianglesUVW();
			if(tr->localPick(b,f,position,triangle,minParam)) {
				picked=true;
				(*nit)->getName(name);
				float vtx[]={position[0],position[1],position[2]};
				transformVector3(vtx,om,position);
			}
		}
		else {	// (*nit)->getType()==sceneNode::TRIANGLES
			// as of now - do nothing. Static triangles just scenery
/*			staticTriangle *tr = static_cast<staticTriangle*>(*nit);
			if(tr->localPick(b,f,position,triangle,minParam)) {
				picked=true;
				(*nit)->getName(name);
				float vtx[]={position[0],position[1],position[2]};
				transformVector3(vtx,om,position);
			} */
		}
	}
/*	if(picked) {  // this whole section looks like it was trying to find closest vertex, which I no longer want
		if(name.find("_")!=1) {
			sceneNode *sc = getNodePtr(name);
			if(sc->getType()==sceneNode::UVW_TRIANGLES){
			}
			else {  // sc->getType()==sceneNode::TRIANGLES)
				glslTriangle *tr = static_cast<glslTriangle*>(sc);
				unsigned int idx,i,*tri = tr->triVerts(triangle);
				float *v,d,dmin=1e30f;
				for(i=0; i<3; ++i) {
					v = tr->vertexCoordinate(tri[i]);
					d = (v[0]-position[0])*(v[0]-position[0]) + (v[1]-position[1])*(v[1]-position[1]) + (v[2]-position[2])*(v[2]-position[2]);
					if(d<dmin) {
						dmin=d;
						idx=i;
					}
				}
			}
		}
		else
			;
	} */
	return picked;
}

void wxGraphics::deleteSceneNode(sceneNode* sn)
{
	if(sn->getType()<sceneNode::TRIANGLES)
		;
//		_shapes.deleteShape(dynamic_cast<shape*>(sn));	// done in deleteShape() call
	else if(sn->getType()<sceneNode::LINES)	{
		std::list<staticTriangle>::iterator git;
		for(git=_sTris.begin(); git!=_sTris.end(); ++git)	{
			if(&(*git)==sn)	{
				_sTris.erase(git);
				break;
			}
		}
	}
	else	{	// is a LINES sceneNode
		if(!_lines.empty())	{
			_lines.clear();
			return;
		}
	}
	std::list<sceneNode*>::iterator nit;
	for(nit=_nodes.begin(); nit!=_nodes.end(); ++nit)	{
		if(sn==*nit) {
			_nodes.erase(nit);
			return;
		}
	}
}

sceneNode* wxGraphics::getNodePtr(std::string &name)
{
	std::list<sceneNode*>::iterator nit;
	std::string nname;
	for(nit=_nodes.begin(); nit!=_nodes.end(); ++nit)	{
		(*nit)->getName(nname);
		if(nname==name)
			return (*nit);
	}
	return NULL;
}

void wxGraphics::clear() {	// empties all graphics
    _sTris.clear();
    _texReader.clear();
    _lines.clear();
    _ls.clear();
    _shapes.clear();
    _nodes.clear();
}

