// File: fences.h
// Author: Court Cutting, MD
// Date: June 7, 2012
// Purpose: User interface for creating a fence on a glslTriangle object.
//     This will be used to specify a desired incision line.

#include "GraphicsUtils/shapes.h"
#include "GraphicsUtils/GLmatrices.h"
#include "GraphicsUtils/Vec3f.h"
#include <stdio.h>
#include <assert.h>
#include "fences.h"

float fence::_fenceSize=10000.0f;
GLfloat fence::_selectedColor[]={1.0f,1.0f,0.0f,1.0f};
GLfloat fence::_unselectedColor[]={0.0f,1.0f,0.0f,1.0f};

void fence::clear()	// deletes current fence
{
	std::list<fencePost>::iterator pit;
	for(pit=_posts.begin(); pit!=_posts.end(); ++pit)	{
		_wxg->deleteSceneNode(pit->cylinderShape);
	}
	_posts.clear();
	if(_slt!=NULL)	{
		_wxg->deleteSceneNode(_slt);
	}
	_slt = NULL;
}

void fence::deleteLastPost()
{
	assert(false); // fix me
	std::list<fencePost>::reverse_iterator pit=_posts.rbegin();
	_wxg->deleteSceneNode(pit->cylinderShape);
	_posts.pop_back();
	if(_posts.empty())
		return;
	int n = _slt->triangleNumber();
	_slt->deleteTriangle(n-1);
	_slt->deleteTriangle(n-2);
	_slt->findAdjacentTriangles();
	_slt->computeNormalsTangents();
	_slt->sendNewColoredTopology();
}

void fence::addPost(trianglesUVW *tri, int triangle, float (&xyz)[3], float (&normal)[3], bool connectToNearestEdge)
{
	char name[6];
	sprintf(name,"P_%d",(int)_posts.size());
	_posts.push_back(fencePost());
	if(connectToNearestEdge && _posts.empty())	{
		int edge;
		float param;
		tri->getNearestHardEdge(xyz,triangle,edge,param);
		_posts.back().triangle = triangle;
		if(edge<1)	{
			_posts.back().uv[0] = param;
			_posts.back().uv[1] = 0.0f;
		}
		else if(edge<2)	{
			_posts.back().uv[0] = 1.0f - param;
			_posts.back().uv[1] = param;
		}
		else	{
			_posts.back().uv[1] = 1.0f - param;
			_posts.back().uv[0] = 0.0f;
		}
	}
	else	{
		_posts.back().triangle = triangle;
		tri->getBarycentricProjection(triangle,xyz[0],xyz[1],xyz[2],_posts.back().uv);
	}
	_posts.back().xyz[0]=xyz[0]; _posts.back().xyz[1]=xyz[1]; _posts.back().xyz[2]=xyz[2];
	_posts.back().nrm[0]=normal[0]; _posts.back().nrm[1]=normal[1]; _posts.back().nrm[2]=normal[2];
	_posts.back().connectToEdge = connectToNearestEdge;
	shape *sh = _posts.back().cylinderShape = _shapes->addShape(shapes::CYLINDER,name);
	GLfloat *mm = sh->getModelViewMatrix();
	loadIdentity4x4(mm);
	sh->setColor(_selectedColor);
	scaleMatrix4x4(mm,_fenceSize*0.1f,_fenceSize*0.1f,_fenceSize*6);	// was 2
	float angle = atan2(-normal[1],normal[2]);
	rotateMatrix4x4(mm,'x',angle);
	if(angle>-1.5708f && angle<1.5708f)
		rotateMatrix4x4(mm,'y',atan2(normal[0],normal[2]));
	else
		rotateMatrix4x4(mm,'y',atan2(-normal[0],-normal[2]));
	translateMatrix4x4(mm,xyz[0],xyz[1],xyz[2]);
	int n,nv;
	GLfloat* vtx;
	if((n=(int)_posts.size())>1)	{
		if(n<3) {
			_slt = _wxg->createEmptyStaticTriangle(false);
			_slt->addCoordNormTexVertices(4);
			vtx = _slt->vertexCoordinate(0);
			vtx[0]=_posts.front().xyz[0]-_posts.front().nrm[0]*_fenceSize*2; 
			vtx[1]=_posts.front().xyz[1]-_posts.front().nrm[1]*_fenceSize*2; 
			vtx[2]=_posts.front().xyz[2]-_posts.front().nrm[2]*_fenceSize*2;
			vtx = _slt->vertexCoordinate(1);
			vtx[0]=_posts.front().xyz[0]+_posts.front().nrm[0]*_fenceSize*2; 
			vtx[1]=_posts.front().xyz[1]+_posts.front().nrm[1]*_fenceSize*2; 
			vtx[2]=_posts.front().xyz[2]+_posts.front().nrm[2]*_fenceSize*2;
		}
		else {
			_slt->addCoordNormTexVertices(2);
		}
		_slt->setColor(_selectedColor);
		nv = (n-1)<<1;
		vtx = _slt->vertexCoordinate(nv);
		vtx[0]=xyz[0]-normal[0]*_fenceSize*2;
		vtx[1]=xyz[1]-normal[1]*_fenceSize*2;
		vtx[2]=xyz[2]-normal[2]*_fenceSize*2;
		vtx = _slt->vertexCoordinate(nv+1);
		vtx[0]=xyz[0]+normal[0]*_fenceSize*2;
		vtx[1]=xyz[1]+normal[1]*_fenceSize*2;
		vtx[2]=xyz[2]+normal[2]*_fenceSize*2;
		unsigned long tr[3];
		tr[0]=nv-2; tr[1]=nv-1; tr[2]=nv;
		_slt->addTriangle(tr);
		tr[0]=nv; tr[1]=nv-1; tr[2]=nv+1;
		_slt->addTriangle(tr);
		_slt->findAdjacentTriangles();
		_slt->computeNormalsTangents();
		_slt->sendNewColoredTopology();
	}
}

int fence::getPostData(std::vector<float> &positions,std::vector<float> &normals,std::vector<int> &triangles,std::vector<float> &uv,bool &edgeStart,bool &edgeEnd)
{
	positions.clear();
	normals.clear();
	triangles.clear();
	uv.clear();
	int n=(int)_posts.size();
	positions.reserve(n*3);
	normals.reserve(n*3);
	triangles.reserve(n);
	uv.reserve(n<<1);
	edgeStart = _posts.front().connectToEdge;
	edgeEnd = _posts.back().connectToEdge;
	std::list<fencePost>::iterator fit;
	for(fit=_posts.begin(); fit!=_posts.end(); ++fit)	{
		positions.push_back(fit->xyz[0]);
		positions.push_back(fit->xyz[1]);
		positions.push_back(fit->xyz[2]);
		normals.push_back(fit->nrm[0]);
		normals.push_back(fit->nrm[1]);
		normals.push_back(fit->nrm[2]);
		triangles.push_back(fit->triangle);
		uv.push_back(fit->uv[0]);
		uv.push_back(fit->uv[1]);
	}
	return n;
}

	fence::fence():_initialized(false)
	{
		_slt = NULL;
		_wxg=NULL;
	}

	fence::~fence()
	{
	}

