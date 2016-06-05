// File: shapes.cpp
// Author: Court Cutting
// Date: 3/3/12
// Purpose: Manages untextured cones, cylinders and spheres for graphics.

#include <vector>
#include "wxGraphics.h"
#include "shapes.h"

GLuint shapes::_coneBufferObjects[]={0,0};
GLuint shapes::_coneVertexArrayBufferObject = 0;
GLuint shapes::_sphereBufferObjects[]={0,0,0};
GLuint shapes::_sphereVertexArrayBufferObject = 0;
GLuint shapes::_cylinderBufferObjects[]={0,0,0};
GLuint shapes::_cylinderVertexArrayBufferObject = 0;

bool shapes::deleteShape(shape* shapePtr)
{
	std::list<shape>::iterator sit;
	for(sit=_shapes.begin(); sit!=_shapes.end(); ++sit) {
		if(shapePtr==&(*sit)) {
			_wgg->deleteSceneNode(shapePtr);
			_shapes.erase(sit);
			return true;
		}
	}
	return false;
}

shape* shapes::addShape(shapeType type, char *name)
{
    return NULL;
/*
	if(type>CYLINDER || type<0)
		return NULL;
	if(type==CONE) {
		getOrCreateConeGraphic();
	}
	else if(type==SPHERE) {
		getOrCreateSphereGraphic();
	}
	else {	// type==CYLINDER
		getOrCreateCylinderGraphic();
	}
	_shapes.push_back(shape(type));
	shape *ret = &_shapes.back();
	_wgg->addSceneNode(ret);
	ret->setGlslProgramNumber(_wgg->getLightsShaders()->getOrCreateColorProgram());
	ret->setName(name);
	return ret;
*/
}

bool shapes::pickSphere(const float *lineStart, float *lineDirection, float (&position)[3], float &param)
{ // assumes param has already been initialized to a large number or to a previous pick
	float t,A,B,C,z;
	bool intrsct=false;
	// intersect with cone
	A = lineDirection[0]*lineDirection[0] + lineDirection[1]*lineDirection[1] + lineDirection[2]*lineDirection[2];
	B = 2.0f*(lineDirection[0]*lineStart[0] + lineDirection[1]*lineStart[1] + lineDirection[2]*lineStart[2]);
	C = lineStart[0]*lineStart[0] + lineStart[1]*lineStart[1] + lineStart[2]*lineStart[2] - 1.0f;
	t = B*B - 4.0f*A*C;
	if(t<0.0f || A==0.0f)
		return false;
	C = sqrt(t);
	for(int i=0; i<2; ++i) {
		if(i<1)
			t = (-B+C)/(2.0f*A);
		else
			t = (-B-C)/(2.0f*A);
		z = lineStart[2]+t*lineDirection[2];
		if(z<1.0f && z>-1.0f) {
			if(param>t) {
				intrsct=true;
				param=t;
				position[0] = lineStart[0]+t*lineDirection[0];
				position[1] = lineStart[1]+t*lineDirection[1];
				position[2] = z;
			}
		}
	}
	return intrsct;
}

bool shapes::pickCylinder(const float *lineStart, float *lineDirection, float (&position)[3], float &param)
{ // assumes param has already been initialized to a large number or to a previous pick
	float t,A,B,C,z;
	bool intrsct=false;
	// intersect with cylinder x^2 +y^2 = 1
	A = lineDirection[0]*lineDirection[0] + lineDirection[1]*lineDirection[1];
	B = 2.0f*(lineDirection[0]*lineStart[0] + lineDirection[1]*lineStart[1]);
	C = lineStart[0]*lineStart[0] + lineStart[1]*lineStart[1] - 1.0f;
	t = B*B - 4.0f*A*C;
	if(t>=0.0f && A!=0.0f) {
		C = sqrt(t);
		for(int i=0; i<2; ++i) {
			if(i<1)
				t = (-B+C)/(2.0f*A);
			else
				t = (-B-C)/(2.0f*A);
			z = lineStart[2]+t*lineDirection[2];
			if(z<1.0f && z>-1.0f) {
				if(param>t) {
					intrsct=true;
					param=t;
					position[0] = lineStart[0]+t*lineDirection[0];
					position[1] = lineStart[1]+t*lineDirection[1];
					position[2] = z;
				}
			}
		}
	}
	// intersect with circles
	if(lineDirection[2]==0.0f)
		return intrsct;
	for(int i=0; i<2; ++i) {
		if(i<1)
			t = (1.0f-lineStart[2])/lineDirection[2];
		else
			t = (-1.0f-lineStart[2])/lineDirection[2];
		if(t>param)
			continue;
		A = t*lineDirection[0]+lineStart[0];
		B = t*lineDirection[1]+lineStart[1];
		if(A*A+B*B<1.0f) {
			if(param>t) {
				intrsct=true;
				param=t;
				position[0] = A;
				position[1] = B;
				if(i<1)
					position[2] = 1.0f;
				else
					position[2] = -1.0f;
			}
		}
	}
	return intrsct;
}

bool shapes::pickCone(const float *lineStart, float *lineDirection, float (&position)[3], float &param)
{ // assumes param has already been initialized to a large number or to a previous pick
	float t,A,B,C,z;
	bool intrsct=false;
	// intersect with cone
	A = lineDirection[0]*lineDirection[0] + lineDirection[1]*lineDirection[1] - lineDirection[2]*lineDirection[2]*0.25f;
	B = 2.0f*(lineDirection[0]*lineStart[0] + lineDirection[1]*lineStart[1] - lineDirection[2]*lineStart[2]*0.25f);
	C = lineStart[0]*lineStart[0] + lineStart[1]*lineStart[1] - lineStart[2]*lineStart[2]*0.25f;
	t = B*B - 4.0f*A*C;
	if(t>=0.0f && A!=0.0f) {
		C = sqrt(t);
		for(int i=0; i<2; ++i) {
			if(i<1)
				t = (-B+C)/(2.0f*A);
			else
				t = (-B-C)/(2.0f*A);
			z = lineStart[2]+t*lineDirection[2];
			if(z<1.0f && z>0.0f) {
				if(param>t) {
					intrsct=true;
					param=t;
					position[0] = lineStart[0]+t*lineDirection[0];
					position[1] = lineStart[1]+t*lineDirection[1];
					position[2] = z;
				}
			}
		}
	}
	// intersect with circle
	if(lineDirection[2]==0.0f)
		return intrsct;
	t = (1.0f-lineStart[2])/lineDirection[2];
	if(t>param)
		return intrsct;
	A = t*lineDirection[0]+lineStart[0];
	B = t*lineDirection[1]+lineStart[1];
	if(A*A+B*B<0.25f) {
		if(param>t) {
			intrsct=true;
			param=t;
			position[0] = A;
			position[1] = B;
			position[2] = 1.0f;
		}
	}
	return intrsct;
}

void shapes::drawCylinder()
{
	glBindVertexArray(_cylinderVertexArrayBufferObject);
	glDrawArrays(GL_TRIANGLE_STRIP,0,42);
	glDrawArrays(GL_TRIANGLE_FAN,42,22);
	glDrawArrays(GL_TRIANGLE_FAN,64,22);
	glBindVertexArray(0);
}

void shapes::drawSphere()
{
	glPrimitiveRestartIndex(0xffffffff);
	glEnable(GL_PRIMITIVE_RESTART);
	glBindVertexArray(_sphereVertexArrayBufferObject);
	//	assumes glUseProgram(_program) has already been called
	glDrawElements(GL_TRIANGLE_STRIP, 419, GL_UNSIGNED_INT, 0);
    // Unbind to anybody
	glBindVertexArray(0);
	glDisable(GL_PRIMITIVE_RESTART);
}

void shapes::drawCone()
{
	glBindVertexArray(_coneVertexArrayBufferObject);
	glDrawArrays(GL_TRIANGLE_FAN,0,22);
	glDrawArrays(GL_TRIANGLE_FAN,22,22);
	glBindVertexArray(0);
}

void shapes::getOrCreateCylinderGraphic()
{ // creates cylinder at origin with diameter 1.0, from z=-1 to z=1
	// Create the buffer objects once
/*
	if(!_cylinderVertexArrayBufferObject)
		glGenVertexArrays(1,&_cylinderVertexArrayBufferObject);
	if(!_cylinderBufferObjects[0])
	    glGenBuffers(3, _cylinderBufferObjects);
	std::vector<GLfloat> vtx,nrm;
	vtx.reserve(344);
	nrm.reserve(258);
	float angle=0.31416f;
	for(int i=0; i<21; ++i) {
		vtx.push_back(cos(i*angle));
		nrm.push_back(vtx.back());
		vtx.push_back(sin(i*angle));
		nrm.push_back(vtx.back());
		vtx.push_back(-1.0f);
		nrm.push_back(0.0f);
		vtx.push_back(1.0f);
		vtx.push_back(cos(i*angle));
		nrm.push_back(vtx.back());
		vtx.push_back(sin(i*angle));
		nrm.push_back(vtx.back());
		vtx.push_back(1.0f);
		nrm.push_back(0.0f);
		vtx.push_back(1.0f);
	}
	vtx.push_back(0.0f);
	vtx.push_back(0.0f);
	vtx.push_back(1.0f);
	vtx.push_back(1.0f);
	nrm.push_back(0.0f);
	nrm.push_back(0.0f);
	nrm.push_back(1.0f);
	float vz=cos(2.0f*angle),vxy=sin(2.0f*angle);
	for(int i=0; i<21; ++i) {
		vtx.push_back(cos(i*angle));
		nrm.push_back(vxy*vtx.back());
		vtx.push_back(sin(i*angle));
		nrm.push_back(vxy*vtx.back());
		vtx.push_back(1.0f);
		nrm.push_back(vz);
		vtx.push_back(1.0f);
	}
	vtx.push_back(0.0f);
	vtx.push_back(0.0f);
	vtx.push_back(-1.0f);
	vtx.push_back(1.0f);
	nrm.push_back(0.0f);
	nrm.push_back(0.0f);
	nrm.push_back(-1.0f);
	for(int i=20; i>-1; --i) {
		vtx.push_back(cos(i*angle));
		nrm.push_back(vxy*vtx.back());
		vtx.push_back(sin(i*angle));
		nrm.push_back(vxy*vtx.back());
		vtx.push_back(-1.0f);
		nrm.push_back(-vz);
		vtx.push_back(1.0f);
	}
    // Vertex and normal data
    glBindBuffer(GL_ARRAY_BUFFER, _cylinderBufferObjects[0]);	// VERTEX_DATA
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*vtx.size(), &(vtx[0]), GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, _cylinderBufferObjects[1]);	// VERTEX_DATA
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*nrm.size(), &(nrm[0]), GL_STATIC_DRAW);
	// Create the master vertex array object
	glBindVertexArray(_cylinderVertexArrayBufferObject);
    // Vertex data
    glBindBuffer(GL_ARRAY_BUFFER, _cylinderBufferObjects[0]);	// VERTEX DATA
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
    // Normal data
    glBindBuffer(GL_ARRAY_BUFFER, _cylinderBufferObjects[1]);	// NORMAL_DATA
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
    // Unbind to anybody
	glBindVertexArray(0);
	// release for next use
    glBindBuffer( GL_ARRAY_BUFFER, 0);
    glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0);
*/
}

void shapes::getOrCreateSphereGraphic()
{ // creates sphere at origin with a radius of 1.0
	// Create the buffer objects once
/*
	if(!_sphereVertexArrayBufferObject)
		glGenVertexArrays(1,&_sphereVertexArrayBufferObject);
	if(!_sphereBufferObjects[0])
	    glGenBuffers(3, _sphereBufferObjects);
	std::vector<GLfloat> vtx,nrm;
	vtx.reserve(728);
	nrm.reserve(546);
	float angle=0.31416f;
	for(int i=1; i<10; ++i) {
		vtx.push_back(sin(i*angle));
		nrm.push_back(vtx.back());
		vtx.push_back(0.0f);
		nrm.push_back(vtx.back());
		vtx.push_back(cos(i*angle));
		nrm.push_back(vtx.back());
		vtx.push_back(1.0f);
	}
	for(int j,i=1; i<20; ++i)	{
		for(j=1; j<10; ++j) {
			vtx.push_back(vtx[(j-1)<<2]*cos(i*angle));
			nrm.push_back(vtx.back());
			vtx.push_back(vtx[(j-1)<<2]*sin(i*angle));
			nrm.push_back(vtx.back());
			vtx.push_back(vtx[((j-1)<<2)+2]);
			nrm.push_back(vtx.back());
			vtx.push_back(1.0f);
		}
	}
	vtx.push_back(0.0f);
	vtx.push_back(0.0f);
	vtx.push_back(1.0f);
	vtx.push_back(1.0f);
	vtx.push_back(0.0f);
	vtx.push_back(0.0f);
	vtx.push_back(-1.0f);
	vtx.push_back(1.0f);
	nrm.push_back(0.0f);
	nrm.push_back(0.0f);
	nrm.push_back(1.0f);
	nrm.push_back(0.0f);
	nrm.push_back(0.0f);
	nrm.push_back(-1.0f);
    // Vertex and normal data
    glBindBuffer(GL_ARRAY_BUFFER, _sphereBufferObjects[0]);	// VERTEX_DATA
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*vtx.size(), &(vtx[0]), GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, _sphereBufferObjects[1]);	// VERTEX_DATA
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*nrm.size(), &(nrm[0]), GL_STATIC_DRAW);
	// Indexes for triangle strips
	std::vector<GLuint> indx;
	int last=19;
	for(int j,i=0; i<20; ++i) {
		indx.push_back(180);
		for(j=0; j<9; ++j) {
			indx.push_back(last*9+j);
			indx.push_back(i*9+j);
		}
		indx.push_back(181);
		indx.push_back(0xffffffff);
		last = i;
	}
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _sphereBufferObjects[2]);	// INDEX_DATA
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*(indx.size()-1), &(indx[0]), GL_STATIC_DRAW);
	// Create the master vertex array object
	glBindVertexArray(_sphereVertexArrayBufferObject);
    // Vertex data
    glBindBuffer(GL_ARRAY_BUFFER, _sphereBufferObjects[0]);	// VERTEX DATA
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
    // Normal data
    glBindBuffer(GL_ARRAY_BUFFER, _sphereBufferObjects[1]);	// NORMAL_DATA
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
    // Indexes
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _sphereBufferObjects[2]);	// INDEX_DATA
    // Unbind to anybody
	glBindVertexArray(0);
	// release for next use
    glBindBuffer( GL_ARRAY_BUFFER, 0);
    glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0);
*/
}

void shapes::getOrCreateConeGraphic()
{
/*
	// Create the buffer objects once
	if(!_coneVertexArrayBufferObject)
		glGenVertexArrays(1,&_coneVertexArrayBufferObject);
	if(!_coneBufferObjects[0])
	    glGenBuffers(2, _coneBufferObjects);
	std::vector<GLfloat> vtx,nrm;
	vtx.assign(44*4,1.0f);
	nrm.assign(44*3,0.0f);
	float angle=0.31416f,mz=-sin(3.1416f/6.0f),mc=cos(3.1416f/6.0f);
	vtx[88]=vtx[89]=nrm[66]=nrm[67]=vtx[0]=vtx[1]=vtx[2]=nrm[0]=nrm[1]=0.0f; // =vtx[164]=vtx[165]
	nrm[2]=-1.0f;
	nrm[68]=1.0f;
	vtx[90]=1.0f;
	float vz=cos(2.0f*angle),vxy=2.0f*sin(2.0f*angle);
	for(int i=0; i<20; ++i)	{
		vtx[(i<<2)+92] = vtx[80-(i<<2)] = 0.5f*cos(angle*i); // vtx[(i<<2)+92] = 
		vtx[(i<<2)+93] = vtx[81-(i<<2)] = 0.5f*sin(angle*i);
		vtx[(i<<2)+94] = vtx[82-(i<<2)] = 1.0f;
		nrm[60-i*3] = mc*cos(angle*i);
		nrm[61-i*3] = mc*sin(angle*i);
		nrm[62-i*3] = mz;
		nrm[i*3+69] = vxy*vtx[(i<<2)+92];
		nrm[i*3+70] = vxy*vtx[(i<<2)+93];
		nrm[i*3+71] = vz;
	}
	vtx[84]=vtx[4];
	vtx[85]=vtx[5];
	vtx[86]=vtx[6];
	nrm[63]=nrm[3];
	nrm[64]=nrm[4];
	nrm[65]=nrm[5];
	vtx[172]=vtx[92];
	vtx[173]=vtx[93];
	vtx[174]=vtx[94];
	nrm[129]=vxy*vtx[92];
	nrm[130]=vxy*vtx[93];
	nrm[131]=vz;
    // Copy data to video memory
    // Vertex data
    glBindBuffer(GL_ARRAY_BUFFER, _coneBufferObjects[0]);	// VERTEX_DATA
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*44*4, &(vtx[0]), GL_DYNAMIC_DRAW);
    // Normal data
    glBindBuffer(GL_ARRAY_BUFFER, _coneBufferObjects[1]);	// NORMAL_DATA
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*44*3, &(nrm[0]), GL_DYNAMIC_DRAW);
	// release for next use
	glBindBuffer( GL_ARRAY_BUFFER, 0);

	// Create the master vertex array object
	glBindVertexArray(_coneVertexArrayBufferObject);
    // Vertex data
    glBindBuffer(GL_ARRAY_BUFFER, _coneBufferObjects[0]);	// VERTEX_DATA
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
//	glEnableClientState( GL_VERTEX_ARRAY);
    // Normal data
    glBindBuffer(GL_ARRAY_BUFFER, _coneBufferObjects[1]);	// NORMAL_DATA
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
    // Unbind to anybody
	glBindVertexArray(0);
	// release for next use
    glBindBuffer( GL_ARRAY_BUFFER, 0);
*/
}

shapes::shapes()
{
}

shapes::~shapes() {
    clear();
}

void shapes::clear()
{
	_shapes.clear();
/*
	if(_coneBufferObjects[0]>0) {
	    glDeleteBuffers(2, _coneBufferObjects);
		_coneBufferObjects[0]=0;
		_coneBufferObjects[1]=0; }
	if(_coneVertexArrayBufferObject>0) {
	    glDeleteVertexArrays(1, &_coneVertexArrayBufferObject);
		_coneVertexArrayBufferObject=0; }
	if(_sphereBufferObjects[0]>0) {
	    glDeleteBuffers(3, _sphereBufferObjects);
		_sphereBufferObjects[0]=0;
		_sphereBufferObjects[1]=0;
		_sphereBufferObjects[2]=0; }
	if(_sphereVertexArrayBufferObject>0) {
	    glDeleteVertexArrays(1, &_sphereVertexArrayBufferObject);
		_sphereVertexArrayBufferObject=0; }
	if(_cylinderBufferObjects[0]>0) {
	    glDeleteBuffers(3, _cylinderBufferObjects);
		_cylinderBufferObjects[0]=0;
		_cylinderBufferObjects[1]=0;
		_cylinderBufferObjects[2]=0; }
	if(_cylinderVertexArrayBufferObject>0) {
	    glDeleteVertexArrays(1, &_cylinderVertexArrayBufferObject);
		_cylinderVertexArrayBufferObject=0; }
*/
}

void shape::computeLocalBounds()
{
    /*
	if(_type<sceneNode::SPHERE) { // _type==sceneNode::CONE
		_localCenter[0]=0.0f; _localCenter[1]=0.0f; _localCenter[2]=0.625f;
		_radius = 0.625f;
	}
	else if(_type<sceneNode::CYLINDER)	{ // _type==sceneNode::SPHERE
		_localCenter[0]=0.0f; _localCenter[1]=0.0f; _localCenter[2]=0.0f;
		_radius = 1.0f;
	}
	else	{ // _type==sceneNode::CYLINDER
		_localCenter[0]=0.0f; _localCenter[1]=0.0f; _localCenter[2]=0.0f;
		_radius = 1.4143f;	// hair bigger than sqrt(2.0f)
	}
    */
	_boundsComputed=true;
}
