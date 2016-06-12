//////////////////////////////////////////////////////////
// File: skinGraphics.cpp
// Author: Court Cutting, MD
// Date: 4/15/2014
// Purpose: Takes as input a trianglesUVW class and makes it visible
//    on an openGL canvas.  It creates hard normal edges and textures
//    the model given an input texture file for the top skin. Skin side
//    texturing and possibly the bottom fat are done procedurally in
//    the fragment shader.
//////////////////////////////////////////////////////////

#include <algorithm>
#include <set>
#include <assert.h>
#include "Vec3f.h"
#include "lightsShaders.h"
#include "boundingBox.h"
#include "wxGraphics.h"
#include "skinGraphics.h"

const GLchar *skinGraphics::skinVertexShader = "#version 140 \n"
	"in vec4 vVertex;"
	"in vec3 vNormal;"
	"in vec3 vTangent;"
	"in vec2 vTexture;"
	"uniform mat4   mvpMatrix;"
	"uniform mat4   mvMatrix;"
	"uniform mat3   normalMatrix;"
	"uniform vec3   vLightPosition;"
	// next 2 are for bump mapping per Randi Rost
	"smooth out vec3 vLightDir;"
	"smooth out vec3 vEyeDir;"
	"smooth out vec2 vTexCoords;"
	"void main(void) "
	"{"
	"   vEyeDir = vec3(mvMatrix * vVertex);"
	"	vec3 n = normalize(normalMatrix * vNormal);\n"
	"	vec3 t = normalize(normalMatrix * vTangent);\n"
	"	vec3 b = cross(n,t);\n"
	"	vec3 v;\n"
	"	v.x = dot(vLightPosition,t);\n"
	"	v.y = dot(vLightPosition,b);\n"
	"	v.z = dot(vLightPosition,n);\n"
	"	vLightDir = normalize(v);\n"
	"	v.x = dot(vEyeDir,t);\n"
	"	v.y = dot(vEyeDir,b);\n"
	"	v.z = dot(vEyeDir,n);"
	"	vEyeDir = normalize(v);"
	"	vTexCoords = vTexture;"
	"   gl_Position = mvpMatrix * vVertex;"
	"}";

const GLchar *skinGraphics::skinFragmentShader =
	// Adapted from Randi Rost and Bill Licea-Kane
	// OpenGL Shading Language 3rd edition
	"#version 140 \n"
	"out vec4 vFragColor;"
	"uniform vec4 ambientColor;"
	"uniform vec4 diffuseColor;"
	"uniform sampler2D colorMap;"
	"uniform sampler2D normalMap;"
	"smooth in vec3 vLightDir;"
	"smooth in vec3 vEyeDir;"
	"smooth in vec2 vTexCoords;"
	"/* Description : Array and textureless GLSL 2D & 3D simplex noise functions."
	"//      Author : Ian McEwan, Ashima Arts."
	"//  Maintainer : ijm"
	"//     Lastmod : 20110822 (ijm)"
	"//     License : Copyright (C) 2011 Ashima Arts. All rights reserved."
	"//               Distributed under the MIT License. See LICENSE file."
	"//               https://github.com/ashima/webgl-noise */"
	"vec3 mod289(vec3 x) {"
	"  return x - floor(x * (1.0 / 289.0)) * 289.0;}"
	"vec2 mod289(vec2 x) {"
	"  return x - floor(x * (1.0 / 289.0)) * 289.0;}"
	"vec3 permute(vec3 x) {"
	"  return mod289(((x*34.0)+1.0)*x);}"
	"vec4 taylorInvSqrt(vec4 r){"
	"  return 1.79284291400159 - 0.85373472095314 * r;}"
	"float snoise(vec2 v)  {"
	"  const vec4 C = vec4(0.211324865405187,"  // (3.0-sqrt(3.0))/6.0
	"                      0.366025403784439,"  // 0.5*(sqrt(3.0)-1.0)
	"                     -0.577350269189626,"  // -1.0 + 2.0 * C.x
	"                      0.024390243902439);" // 1.0 / 41.0
	// First corner
	"  vec2 i  = floor(v + dot(v, C.yy) );"
	"  vec2 x0 = v -   i + dot(i, C.xx);"
	// Other corners
	"  vec2 i1;"
	  //i1.x = step( x0.y, x0.x ); // x0.x > x0.y ? 1.0 : 0.0
	  //i1.y = 1.0 - i1.x;
	"  i1 = (x0.x > x0.y) ? vec2(1.0, 0.0) : vec2(0.0, 1.0);"
	  // x0 = x0 - 0.0 + 0.0 * C.xx ;
	  // x1 = x0 - i1 + 1.0 * C.xx ;
	  // x2 = x0 - 1.0 + 2.0 * C.xx ;
	"  vec4 x12 = x0.xyxy + C.xxzz;"
	"  x12.xy -= i1;"
	// Permutations
	"  i = mod289(i);" // Avoid truncation effects in permutation
	"  vec3 p = permute( permute( i.y + vec3(0.0, i1.y, 1.0 ))"
	"		+ i.x + vec3(0.0, i1.x, 1.0 ));"
	"  vec3 m = max(0.5 - vec3(dot(x0,x0), dot(x12.xy,x12.xy), dot(x12.zw,x12.zw)), 0.0);"
	"  m = m*m ;"
	"  m = m*m ;"
	// Gradients: 41 points uniformly over a line, mapped onto a diamond.
	// The ring size 17*17 = 289 is close to a multiple of 41 (41*7 = 287)
	"  vec3 x = 2.0 * fract(p * C.www) - 1.0;"
	"  vec3 h = abs(x) - 0.5;"
	"  vec3 ox = floor(x + 0.5);"
	"  vec3 a0 = x - ox;"
	// Normalise gradients implicitly by scaling m
	// Approximation of: m *= inversesqrt( a0*a0 + h*h );
	"  m *= 1.79284291400159 - 0.85373472095314 * ( a0*a0 + h*h );"
	// Compute final noise value at P
	"  vec3 g;"
	"  g.x  = a0.x  * x0.x  + h.x  * x0.y;"
	"  g.yz = a0.yz * x12.xz + h.yz * x12.yw;"
	"  return 130.0 * dot(m, g);}"

	// next routine gets fat from 4 surrounding noise values [-1,1]
	"void getFat(in vec4 nei, out vec4 fatColor, out vec3 normDelta, out float specMult) {"
	"  float h;"
	"  for(int i=0; i<4; ++i){"
	"    h = 1.0 - abs(nei[i]);"
	"    h *= h;"
	"    nei[i] = 1.0 - h; }"
	"  h = 0;"
	"  for(int i=0; i<4; ++i)"
	"    h += nei[i];"
	"  h *= 0.25;"
	"  vec2 p;"
	"  p.x = nei[1]-nei[0];"
	"  p.y = nei[2]-nei[3];"
	"  p *= 130.0;"
	"  p = clamp(p,-1.0,1.0);"
	"  float d,f;"
	"  d = dot(p,p);"
	"  f = inversesqrt(d+1.0);"
	"  p.x = -p.x;"
	"  normDelta = vec3(p,1.0)*f;"
	"  float fatRed, fatGreen, fatBlue;"
	"  if(h<0.04) {"
	"    specMult = 0.2;"
	"    fatRed = (1.0-h)*0.4;"
	"    fatBlue = 0.0;"
	"    fatGreen = (1.0-h)*0.2; }"
	"  else {"
	"    specMult = 1.0;"
	"    fatRed = 0.5 + h*0.8;"
	"    fatBlue = 0.15 + h*0.3;"
	"    fatGreen = 0.35 + h*0.8; }"
	"  fatColor = vec4(fatRed, fatGreen, fatBlue, 1.0); }"

	"void main(void)"
	"{"
	// always lighting with white light
	"	const float ambientVal = 0.1;"
	"	const float	fatIncr = 0.5/1024.0;"
	"   vec2 fatD = vec2(fatIncr,0.0);"
	"   vec2 faceUV;"
	"	float lightVal = ambientVal;"
	"	float sn,h,specMult = 0.5;"
	"	const float diffuseVal = 0.9;"
	"	vec3 normDelta = vec3(0.0, 0.0, 1.0);"
	"	vec3 litColor = vec3(1.0, 1.0, 1.0);"
	"	vec4 nei;"
	"	if(vTexCoords.s>1.05f) {"
	"		if(vTexCoords.s<2.15f) {"  // middle
	"			faceUV = vTexCoords - vec2(1.1,0.0);"
	"			sn = snoise(vec2(faceUV.t*3.0,0.5));"
	"			float val = (sn+1.0)*0.5;"
	"			if(0.40 + 0.2*val > faceUV.s) {"  // fat
	"				faceUV *= vec2(6.4,4.7);"
	"				nei[0] = snoise(faceUV + fatD);"
	"				nei[1] = snoise(faceUV - fatD);"
	"				nei[3] = snoise(faceUV + fatD.yx);"
	"				nei[2] = snoise(faceUV - fatD.yx);"
	"				getFat(nei,vFragColor,normDelta,specMult);"
	"			}"
	"			else if(0.415 + 0.2*val > faceUV.s) {"  // dermal-fat junction
	"				vFragColor = vec4(0.51, 0.44, 0.1412, 1.0);"
	"				specMult = 0.2; }"
	"			else {"  // dermis
	"				specMult = 0.35;"
	"				if(faceUV.s>0.95)"
	"					vFragColor = vec4(0.71, 0.57255, 0.2784, 1.0);"
	"				else if(faceUV.s>0.93)"
	"					vFragColor = vec4(0.843, 0.737, 0.51, 1.0);"
	"				else {"
	"					sn = snoise(vec2((faceUV.t+0.4)*4.2,0.5));"
	"					val = (sn-0.5)*2.0;"
	"					if(0.85 + 0.05*val > faceUV.s)"
	"						vFragColor = vec4(0.9569, 0.8902, 0.71, 1.0);"
	"					else {"
	"						vFragColor = vec4(0.7255, 0.5059, 0.2039, 1.0);"
	"						vFragColor = vFragColor*(3.3-2.8*faceUV.s); }"
	"			} }"
	"		}"
	"		else {"  // bottom
	"			faceUV = vTexCoords - vec2(2.2,0.0);"
	"			nei[0] = snoise(faceUV - fatD);"
	"			nei[1] = snoise(faceUV + fatD);"
	"			nei[2] = snoise(faceUV - fatD.yx);"
	"			nei[3] = snoise(faceUV + fatD.yx);"
	"			getFat(nei,vFragColor,normDelta,specMult);"
	"	} }"
	"	else {"	// top of skin
	"	  vec4 tx1 = texture(normalMap,vTexCoords.st);"
	"	  tx1.rgb -= vec3(0.5);"
	"	  normDelta = tx1.rgb*2.0;"
	"	  specMult = 0.11;"
	"	  vFragColor = texture(colorMap, vTexCoords.st); }"
	"	lightVal += diffuseVal*max(dot(normDelta,vLightDir), 0.0);"
	"	vFragColor *= lightVal;"
	"	vec3 reflectDir = reflect(vLightDir,normDelta);"
	"	float spec = max(dot(vEyeDir,reflectDir),0.0);"
	"	spec = pow(spec,40.0);"
	"	spec *= specMult;"
	"	litColor = min(vFragColor.rgb + spec, vec3(1.0));"
	"	vFragColor = vec4(litColor, 1.0);"
	"}";

bool skinGraphics::setTopTextureFilesCreateProgram(const char *topTextureFile, const char *topNormalFile)
{  // must be set first
	if(_glslProgram>0)
		return false;
	std::vector<std::string> att;
	att.assign(4,std::string());
	att[0] = "vVertex";
	att[1] = "vNormal";
	att[2] = "vTangent";
	att[3] = "vTexture";
	//bool ret = _wxg->addCustomSceneNode(this,topTextureFile,topNormalFile,skinVertexShader,skinFragmentShader,att);
	//if(!ret)
	//	return ret;
	//_wxg->getLightsShaders()->useGlslProgram(_glslProgram);  // must be current program. This routine sets other uniforms.
	if(!_bufferObjects[0])
		//glGenBuffers(5, _bufferObjects);
//	else	{  // COURT - is this necessary or can I just realloc in same objects?
//			glDeleteBuffers(5, _bufferObjects);
//			glGenBuffers(5, _bufferObjects);	}
	// Vertex data
	//glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[0]);	// VERTEX_DATA
	//glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_xyz1.size(), &(_xyz1[0]), GL_DYNAMIC_DRAW);
	// Normal data
	//glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[1]);	// NORMAL_DATA
	//glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_normals.size(), &(_normals[0]), GL_DYNAMIC_DRAW);
	// Tangent data
	//glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[2]);	// TANGENT_DATA
	//glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_tangents.size(), &(_tangents[0]), GL_DYNAMIC_DRAW);
    // Texture coordinates
	//glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[3]);	// TEXTURE_DATA
	//glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_uv.size(), &(_uv[0]), GL_STATIC_DRAW);
    //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bufferObjects[4]);	// INDEX_DATA
	//glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*_tris.size(), &(_tris[0]), GL_STATIC_DRAW);
	if(!_vertexArrayBufferObject)
		//glGenVertexArrays(1,&_vertexArrayBufferObject);
	// now make vertex array
	//glBindVertexArray(_vertexArrayBufferObject);
    // Position data
    //glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[0]);
	//glEnableVertexAttribArray(0);
	//glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
    // Normal data
    //glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[1]);
	//glEnableVertexAttribArray(1);
	//glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
	// Tangent data
    //glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[2]);
	//glEnableVertexAttribArray(2);
	//glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, 0);
    // Texture coordinates
	//glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[3]);	// TEXTURE_DATA
	//glEnableVertexAttribArray(3);
	//glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, 0, 0);
    // Indexes
    //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bufferObjects[4]);	// INDEX_DATA */
	// never unbind a GL_ARRAY_BUFFER or GL_ELEMENT_ARRAY_BUFFER inside a vertexArrayBuffer
	//glBindVertexArray(0);
	return true;
}

void skinGraphics::updatePositionsAndNormals()
{
	int nv;
	const float *vp,*fp = _tuvw.getPositionArray(nv);
	GLfloat *gfp,*dp;
	for(int i=0; i<nv; ++i) {
		vp = &fp[(i<<1)+i];
		gfp = &_xyz1[i<<2];
		gfp[0]=vp[0]; gfp[1]=vp[1]; gfp[2]=vp[2];
	}
	nv = (int)_doubledVerts.size();
	for(int i=0; i<nv; ++i) {
		gfp = &_xyz1[_doubledVerts[i]<<2];
		dp = &_xyz1[(_origVertTop+i)<<2];
		dp[0]=gfp[0]; dp[1]=gfp[1]; dp[2]=gfp[2];
	}
	computeNormalsTangents();
	// COURT - don't do this way. Use mapBuffer() instead.  Interestingly interet speed tests show only ~30% faster
	// I guess comparing graphics card memory realloc vs single vs block data transfer
	// about to switch to javascript graphics so why bother
	// Vertex data
	//glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[0]);	// VERTEX_DATA
	//glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_xyz1.size(), &(_xyz1[0]), GL_DYNAMIC_DRAW);
	// Normal data
	//glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[1]);	// NORMAL_DATA
	//glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_normals.size(), &(_normals[0]), GL_DYNAMIC_DRAW);
	// Tangent data
	//glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[2]);	// TANGENT_DATA
	//glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_tangents.size(), &(_tangents[0]), GL_DYNAMIC_DRAW);

}

void skinGraphics::setNewTopology()
{
	_tris.clear();
	_xyz1.clear();
	_uv.clear();
	_normals.clear();
	_tangents.clear();
	_doubledVerts.clear();
	int n;
	const int *trArr = _tuvw.getTriangleArray(n);
	n *= 3;
	_tris.reserve(n);
	for(int i=0; i<n; ++i) {
		if(trArr[i]<0)  // deleted triangle
			_tris.push_back(0xffffffff);
		else
			_tris.push_back(trArr[i]);
	}
	float *txp;
	const float *pArr = _tuvw.getPositionArray(n);
	_origVertTop=n;
	_xyz1.reserve(n<<2);
	_uv.reserve(n<<1);
	_normals.assign(n*3,-1.0f);
	_tangents.assign(n*3,-1.0f);
	for(int i=0; i<n; ++i) {
		txp =_tuvw.vertexTexture(i);
		if(txp[2]==0.0f)
			_uv.push_back(txp[0]+2.2f);  // mark bottom
		else if(txp[2]<1.0f)
			_uv.push_back(txp[2]+1.1f);  // mark middle, assign u to w
		else
			_uv.push_back(txp[0]);  // top
		_uv.push_back(txp[1]);
	}
	n *= 3;
	for(int j=0,i=0; i<n; ++i) {
		_xyz1.push_back(pArr[i]);
		j++;
		if(j>2) {
			j=0;  _xyz1.push_back(1.0f);  // homogeneous coords
		}
	}
	GLuint *tr;
	float txw[3];
	// find hard edges for duplicating normals
//	std::vector<unsigned int> hardEdges; // COURT - may not need this anymore
//	std::list<std::list<GLuint> > paths;
//	std::list<std::list<GLuint> >::iterator pit,pit2;
	std::set<GLuint> vSet;
	n = (int)_tris.size()/3;
	_tuvw.findAdjacentTriangles();
	for(int i=0; i<n; ++i) {
		tr =  &_tris[i*3];
		if(*tr>0xfffffffe)
			continue;
		txw[0] = _tuvw.vertexTexture(tr[0])[2];
		txw[1] = _tuvw.vertexTexture(tr[1])[2];
		txw[2] = _tuvw.vertexTexture(tr[2])[2];
		for(int j=0; j<3; ++j) {
			if(txw[j]==txw[(j+1)%3] && (txw[j]==1.0f || txw[j]==0.0f) && txw[j]!=txw[(j+2)%3]) {
				vSet.insert(tr[j]);
				vSet.insert(tr[(j+1)%3]);
//				if(txw[j]==1.0f)
//					hardEdges.push_back(_tuvw.triAdjs(i)[j]);
/*				for(pit=paths.begin(); pit!=paths.end(); ++pit) {
					if(tr[j]==pit->front()) {
						assert(false);
						pit->push_front(tr[(j+1)%3]);
						break; }
					else if(tr[j]==pit->back()) {
						pit->push_back(tr[(j+1)%3]);
						break; }
					else if(tr[(j+1)%3]==pit->back()) {
						assert(false);
						pit->push_back(tr[j]);
						break; }
					else if(tr[(j+1)%3]==pit->front()) {
						pit->push_front(tr[j]);
						break; }
					else
						;
				}
				if(pit==paths.end()) {
					paths.push_back(std::list<GLuint>());
					paths.back().push_back(tr[j]);
					paths.back().push_back(tr[(j+1)%3]);
				} */
			}
		}
	}
/*	bool spliceIt=true;
	while(spliceIt) {
		spliceIt=false;
		for(pit=paths.begin(); pit!=paths.end(); ++pit) {
			pit2=pit; ++pit2;
			while(pit2!=paths.end()) {
				if(pit2->back()==pit->front()) {
					pit2->pop_back();
					pit->splice(pit->begin(),*pit2);
					if(pit->front()==pit->back())
						pit->pop_back();
					pit2 = paths.erase(pit2);
					spliceIt=true;
				}
				else if(pit2->front()==pit->back()) {
					pit2->pop_front();
					pit->splice(pit->end(),*pit2);
					if(pit->front()==pit->back())
						pit->pop_back();
					pit2 = paths.erase(pit2);
					spliceIt=true;
				}
				else
					++pit2;
			}
		}
	} */
	// create duplicated vertices
	std::set<GLuint>::iterator vit;
	_doubledVerts.clear();
	_doubledVerts.reserve(vSet.size());
	for(vit=vSet.begin(); vit !=vSet.end(); ++vit) {
		_doubledVerts.push_back(*vit);
		float *fp = _tuvw.vertexCoordinate(*vit);
		_xyz1.push_back(fp[0]); _xyz1.push_back(fp[1]); _xyz1.push_back(fp[2]);
		_xyz1.push_back(1.0f);  // for homogeneous coordinates
		fp = _tuvw.vertexTexture(*vit);
		_uv.push_back(fp[2]+1.1f); _uv.push_back(fp[1]);  // w coord in u for these middle vertices and mark middle
	}
	auto newVert = [&](GLuint &tv) {
		std::vector<GLuint>::iterator dit;
	    dit = std::lower_bound(_doubledVerts.begin(),_doubledVerts.end(),tv);
	    if(dit!=_doubledVerts.end())
			tv = (GLuint)(dit-_doubledVerts.begin())+_origVertTop;
		else
			assert(false);
	};
	std::set<GLuint> topVerts;
// COURT - next section is dumb. Just finds all _doubledVerts with a w==1
	n = (int)_tris.size();
	for(int j,i=0; i<n; i+=3) {
		tr =  &_tris[i];
		if(*tr>0xfffffffe)
			continue;
		txw[0] = _tuvw.vertexTexture(tr[0])[2];
		txw[1] = _tuvw.vertexTexture(tr[1])[2];
		txw[2] = _tuvw.vertexTexture(tr[2])[2];
		if(txw[0]==txw[1] && txw[0]==txw[2]) {
			assert(txw[0]==0.0f || txw[0]==1.0f);
			continue;
		}
		for(j=0; j<3; ++j) {
			if(txw[j]==1.0f || txw[j]==0.0f) {
				newVert(tr[j]);
				if(txw[j]==1.0f)
					topVerts.insert(tr[j]);
			}
		}
	}
	// get all top paths
	findAdjacentTriangles();
	std::list<std::list<GLuint> > topPaths;
	topPaths.push_back(std::list<GLuint>());
	GLuint vNow;
	if(!topVerts.empty())
		vNow = *topVerts.begin();
	std::vector<neighborNode> nei;
	while(!topVerts.empty()) {
		topVerts.erase(vNow);
		topPaths.back().push_back(vNow);
		getNeighbors(vNow,nei);
		assert(nei.front().triangle<0);
		vNow = nei.back().vertex;
		if(topVerts.empty()) {
//			assert(vNow==topPaths.back().front()); // not true in cases of open holes
			break; }
		if(nei.back().vertex==topPaths.back().front()) {
			topPaths.push_back(std::list<GLuint>());
			vNow = *topVerts.begin(); }
	}
	std::list<std::list<GLuint> >::iterator psit;
	std::list<GLuint>::iterator pit;
	GLfloat dist;
	Vec3f *lastV,*nowV;
	for(psit=topPaths.begin(); psit!=topPaths.end(); ++psit) {
		dist=0.0f;
		_uv[(psit->back()<<1)+1] = dist;
		lastV = (Vec3f*)&_xyz1[psit->back()<<2];
		for(pit=psit->begin(); pit!=psit->end(); ++pit) {
			nowV = (Vec3f*)&_xyz1[(*pit)<<2];
			dist += (*nowV-*lastV).length();
			_uv[((*pit)<<1)+1] = dist;
			lastV = nowV;
		}
		if(_meanTriangleEdgeLength<0.0f) { // do first time only on load
			_meanTriangleEdgeLength = dist/psit->size();
		}
	}
	std::vector<neighborNode>::iterator nit;
	std::list<std::vector<GLuint> > pVerts;
	for(psit=topPaths.begin(); psit!=topPaths.end(); ++psit) {
		std::set<GLuint> pSet;
		pVerts.push_back(std::vector<GLuint>());
		pSet.insert(psit->begin(),psit->end());
		getNeighbors(*psit->begin(),nei);

//		assert(nei.size()>2);
//		recurseMiddle(nei[1].vertex,pSet,pVerts.back());

		for(nit=nei.begin(); nit!=nei.end(); ++nit) {
			if(nit->vertex<0)
				continue;
			recurseMiddle(nit->vertex,pSet,pVerts.back());
		}
	}
	std::list<std::vector<GLuint> >::iterator pvit=pVerts.begin();
	for(psit=topPaths.begin(); psit!=topPaths.end(); ++psit) {
		n = (int)pvit->size();
		for(int j,i=0; i<n; ++i) { // mxn - this could be done better
			j = (*pvit)[i];
			_uv[(j<<1)+1] = getClosestTopV(j,*psit); }
		++pvit;
	}
//	_computeTangents = false;
//	computeNormalsTangents();
	// now duplicate vertices in creases in mid zone
	// current duplicated vertices are on top and bottom edges
	findAdjacentTriangles();
	vSet.clear();
	Vec3f v0,v1;
	bool firstCrease;
	GLfloat maxV;
	pvit=pVerts.begin();
	for(psit=topPaths.begin(); psit!=topPaths.end(); ++psit) {
		firstCrease=true;
		GLfloat maxDist = _uv[(psit->back()<<1)+1];
		for(pit=psit->begin(); pit!=psit->end(); ++pit){
			getNeighbors(*pit,nei);
			nit=nei.begin();
			assert(nit->triangle<0);
			++nit;
			_tuvw.getTriangleNormal(nit->triangle,v0._v);
			_tuvw.getTriangleNormal(nei.back().triangle,v1._v);
			// next should insure that creases aren't made at edge of an open hole.  If so rest for it and bypass.
			if(v0*v1<0.3f) { // kink in straight line
				if(makeCrease(*pit,firstCrease,maxDist,maxV) && firstCrease) {
					firstCrease = false;
					pvit->insert(pvit->end(),psit->begin(),psit->end());
					n = (int)pvit->size();
					maxV -= 1e-5f;
					for(int j,i=0; i<n; ++i) {
						j = (*pvit)[i];
						if(_uv[(j<<1)+1]<maxV)
							_uv[(j<<1)+1] += maxDist;
					}
				}
				makeVertexToTriangleMap();
			}
		}
		++pvit;
	}
	std::vector<GLfloat>(_xyz1).swap(_xyz1);
	std::vector<GLfloat>(_uv).swap(_uv);
	std::vector<GLfloat>(_normals).swap(_normals);
	std::vector<GLfloat>(_tangents).swap(_tangents);
	_computeTangents = true;
	_getUvScale=true;
	computeNormalsTangents();
	n = (int)_uv.size();
	for(int i=0; i<n; i+=2) { // scale uv based on object size
		if(_uv[i]>1.05f) {
			_uv[i+1] *= _uvScale;
			if(_uv[i]>2.15f) { // bottom
				_uv[i] -= 2.2f;
				_uv[i] *= _uvScale;
				_uv[i] += 2.2f;
			}
		}
	}
	// Vertex data
	//glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[0]);	// VERTEX_DATA
	//glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_xyz1.size(), &(_xyz1[0]), GL_DYNAMIC_DRAW);
	// Normal data
	//glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[1]);	// NORMAL_DATA
	//glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_normals.size(), &(_normals[0]), GL_DYNAMIC_DRAW);
	// Tangent data
	//glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[2]);	// TANGENT_DATA
	//glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_tangents.size(), &(_tangents[0]), GL_DYNAMIC_DRAW);
    // Texture coordinates
	//glBindBuffer(GL_ARRAY_BUFFER, _bufferObjects[3]);	// TEXTURE_DATA
	//glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*_uv.size(), &(_uv[0]), GL_STATIC_DRAW);
    //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bufferObjects[4]);	// INDEX_DATA
	//glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*_tris.size(), &(_tris[0]), GL_STATIC_DRAW);
}

void skinGraphics::recurseMiddle(int v, std::set<GLuint> &m, std::vector<GLuint> &newV)
{
	if(!m.insert(v).second)
		return;
	newV.push_back(v);
	std::vector<neighborNode> nei;
	std::vector<neighborNode>::iterator nit;
	getNeighbors(v,nei);
	for(nit=nei.begin(); nit!=nei.end(); ++nit) {
		if(nit->vertex<0)
			continue;
		recurseMiddle(nit->vertex,m,newV);
	}
}

float skinGraphics::getClosestTopV(int v, std::list<GLuint> &path)
{
	float minDsq=1e30f,minT=1e30f,*vp=&_xyz1[v<<2];
	Vec3f de,pt(vp[0],vp[1],vp[2]),e0,e1,v0,v1;
	e0.set((float (&)[3])_xyz1[path.back()<<2]);
	v0=pt-e0;
	std::list<GLuint>::iterator lit,clit;
	for(lit=path.begin(); lit!=path.end(); ++lit)	{
		e1.set((float (&)[3])_xyz1[(*lit)<<2]);
		v1=pt-e1;
		de = e1 - e0;
		assert(de*de>0.0f);
		float d,t=(v0*de)/(de*de);
		if(t<=0.0f ){
			d=v0*v0;
			if(minDsq>d) {
				clit = lit;
				minT = 0.0f;
				minDsq=d;
			}
		}
		else if(t>=1.0f){
			d=v1*v1;
			if(minDsq>d) {
				clit = lit;
				minT = 1.0f;
				minDsq=d;
			}
		}
		else {
			v0 = e0 + de*t;
			v0 -= pt;
			d=v0*v0;
			if(d<minDsq) {
				clit = lit;
				minT = t;
				minDsq=d;
			}
		}
		e0 = e1;
		v0 = v1;
	}
	minDsq = 0.0f;
	if(clit==path.begin()) {
		lit = path.end();
		if(minT<0.95f)
			minDsq = _uv[((path.back())<<1)+1]*minT;
	}
	else
		lit = clit;
	--lit;
	minDsq += _uv[((*lit)<<1)+1]*(1.0f-minT);
	minDsq += _uv[((*clit)<<1)+1]*minT;
	return minDsq;
}

bool skinGraphics::makeCrease(GLuint topV, bool firstCrease, GLfloat maxDist, GLfloat &maxV)
{
	std::vector<neighborNode> nei;
	std::vector<neighborNode>::iterator nit,minit,startit;
	GLfloat *v,txTarget;
	GLuint *tr;
	std::vector<GLuint> crease;
	std::set<int> cSideTris;
	cSideTris.insert(-1);
	crease.push_back(topV);
	if(firstCrease)
		maxV = _uv[(topV<<1)+1];
	txTarget = _uv[topV<<1];
	if(txTarget==2.1f)
		txTarget=1.1f;
	else
		txTarget=2.1f;
	Vec3f v0,v1;
	// finished initializing
	do {
		getNeighbors(topV,nei);
		startit = nei.end();
		nit=nei.begin();
		while(true){
			if(cSideTris.find(nit->triangle)!=cSideTris.end())
				startit=nit;
			if(startit!=nei.end() && cSideTris.find(nit->triangle)==cSideTris.end()) {
				startit=nit;
				break;
			}
			++nit;
			if(nit==nei.end())
				nit=nei.begin();
		}
		nit=startit;
		_tuvw.getTriangleNormal(nit->triangle,v0._v);
		++nit;
		if(nit==nei.end())
			nit = nei.begin();
		float d,minDot=1e30f;
		while(cSideTris.find(nit->triangle)==cSideTris.end()) {
			_tuvw.getTriangleNormal(nit->triangle,v1._v);
			d = v0*v1;
			if(minDot>d) {
				minDot=d;
				minit=nit;
			}
			v0 = v1;
			++nit;
			if(nit==nei.end())
				nit=nei.begin();
		}
		nit=startit;
		while(nit!=minit){
			cSideTris.insert(nit->triangle);
			++nit;
			if(nit==nei.end())
				nit=nei.begin();
		}
		if(minit==nei.begin())
			topV = nei.back().vertex;
		else {
			--minit;
			topV = minit->vertex;
		}
		if(_uv[crease.back()<<1]>_uv[topV<<1] && txTarget>1.7f) // no backward movement
			break;
		if(_uv[topV<<1]>_uv[crease.back()<<1] && txTarget<1.7f)
			break;
		crease.push_back(topV);
	}while(txTarget !=_uv[topV<<1]);
	if(txTarget !=_uv[topV<<1])
		return false;
	cSideTris.erase(-1);
	// get last triangles attached to next(last point on other edge)
	getNeighbors(topV,nei);
	if(nei.front().triangle>-1)
		return false;
	startit = nei.end();
	nit=nei.begin();
	++nit;
	while(nit!=nei.end()){
		if(cSideTris.find(nit->triangle)!=cSideTris.end())
			startit=nit;
		if(startit!=nei.end() && cSideTris.find(nit->triangle)==cSideTris.end())
			cSideTris.insert(nit->triangle);
		++nit;
	}
	// now duplicate vertices along this crease
	std::map<GLuint,GLuint> crDups;
	std::vector<GLuint>::iterator crit;
	for(crit=crease.begin(); crit!=crease.end(); ++crit) {
		crDups.insert(std::make_pair(*crit,(GLuint)_doubledVerts.size()+_origVertTop));
		if(*crit<_origVertTop)
			_doubledVerts.push_back(*crit);
		else
			_doubledVerts.push_back(_doubledVerts[*crit-_origVertTop]);
		v = &_xyz1[_doubledVerts.back()<<2];
		GLfloat x=v[0],y=v[1],z=v[2];
		_xyz1.push_back(x); _xyz1.push_back(y); _xyz1.push_back(z);
		_xyz1.push_back(1.0f);  // for homogeneous coordinates
		v = &_uv[(*crit)<<1];
		_uv.push_back(v[0]); _uv.push_back(v[1]);
		if(firstCrease)
			_uv.back() += maxDist;
	}
	std::set<int>::iterator tit;
	std::map<GLuint,GLuint>::iterator dit,ditLast;
	for(tit=cSideTris.begin(); tit!=cSideTris.end(); ++tit) {
		tr = &_tris[*tit*3];
		ditLast = crDups.find(tr[2]);
		if(ditLast!=crDups.end())
			tr[2] = ditLast->second;
		for(int i=0; i<3; ++i) {
			dit = crDups.find(tr[i]);
			if(dit!=crDups.end())
				tr[i] = dit->second;
			if(dit!=crDups.end() && ditLast!=crDups.end()) {  // triangle edge on crease found. Unconnect adjacent triangles.
				unsigned long trAdj = _adjs[*tit*3+((i+2)%3)];
				_adjs[*tit*3+((i+2)%3)] = 0x03;
				_adjs[(trAdj>>2)*3+(trAdj&0x03)] = 0x03;
				break;
			}
			ditLast = dit;
		}
	}
	return true;
}

void skinGraphics::computeNormalsTangents()
{
	unsigned int i,j,k,n;
	GLfloat *gv[3],*gtx[3];
	_normals.assign((_uv.size()>>1)*3,0.0f);
	_tangents.assign((_uv.size()>>1)*3,0.0f);
	n = (unsigned int)_tris.size();
	GLuint *tr;
	Vec3f tanV,nrmV,dXyz[2];
	float frct,dUv[2][2],d2;
	if(_getUvScale)
		_uvScale=0.0f;
	for(i=0; i<n; i+=3) {
		tr = &_tris[i];
		if(*tr>0xfffffffe)
			continue;
		for(j=0; j<3; ++j) {
			nrmV[j]=0.0f;
			tanV[j]=0.0f;
			gv[j]=&_xyz1[tr[j]<<2];
			gtx[j]=&_uv[tr[j]<<1];
		}
		for(j=0; j<3; ++j) {
			dXyz[0][j] = gv[1][j]-gv[0][j];
			dXyz[1][j] = gv[2][j]-gv[0][j];
			if(_computeTangents && j<2) {
				dUv[0][j] = gtx[1][j]-gtx[0][j];
				dUv[1][j] = gtx[2][j]-gtx[0][j]; }
		}
		nrmV = dXyz[0]^dXyz[1];
		for(j=0; j<3; ++j) {
			_normals[tr[j]*3] += nrmV[0];
			_normals[tr[j]*3+1] += nrmV[1];
			_normals[tr[j]*3+2] += nrmV[2]; }
		if(!_computeTangents)
			continue;
		if(fabs(dUv[0][0])>fabs(dUv[1][0])) {
			j=0; k=1; }
		else {
			j=1; k=0; }
		frct = dUv[k][0]/dUv[j][0];
		dUv[j][1] *= frct;
		dXyz[j] *= frct;
		dUv[k][1] -= dUv[j][1];
		dXyz[k] -= dXyz[j];
		if(dUv[k][1]>0)
			tanV = dXyz[k];
		else
			tanV = -dXyz[k];
		if(_getUvScale) {
			if(fabs(dUv[k][1])>1e-16f) {
				dXyz[k] /= dUv[k][1];
				_uvScale += dXyz[k].length();	}
		}
		for(j=0; j<3; ++j) {
			_tangents[tr[j]*3] += tanV[0];
			_tangents[tr[j]*3+1] += tanV[1];
			_tangents[tr[j]*3+2] += tanV[2]; }
	}
	n = (int)(_uv.size()>>1)*3;
	for(i=0; i<n; i+=3) {
		d2 = _normals[i]*_normals[i] + _normals[i+1]*_normals[i+1] + _normals[i+2]*_normals[i+2];
		if(d2<1e-16f) {
			_normals[i]=0.0f; _normals[i+1]=0.0f; _normals[i+2]=1.0f; }
		else {
			d2 = 1.0f/sqrt(d2);
			_normals[i]*=d2; _normals[i+1]*=d2; _normals[i+2]*=d2;}
		if(!_computeTangents)
			continue;
		Vec3f *np=(Vec3f*)&_normals[i],*bp=(Vec3f*)&_tangents[i];
		tanV = *bp^*np;
		d2 = tanV.length2();
		if(d2<1e-16f) {
			_tangents[i]=0.0f; _tangents[i+1]=-1.0f; _tangents[i+2]=0.0f; }
		else {
			d2 = 1.0f/sqrt(d2);
			tanV *= d2;
			_tangents[i]=tanV[0]; _tangents[i+1]=tanV[1]; _tangents[i+2]=tanV[2];}
	}
	if(_getUvScale && _computeTangents) {
		_uvScale /= _uv.size()>>1;
		_uvScale *= 30.0f;
		_getUvScale = false;
	}
}

bool skinGraphics::findAdjacentTriangles()
{	// computes all the adjacent triangles from raw triangle input
	// create binary tree to store and find edges efficiently
	// returns false if non-manifold surface is input
	typedef std::set<edge,edgeTest> edgeSet;
	typedef edgeSet::iterator edgeIt;
	std::pair <edgeIt,bool> P;
	edge E;
	edgeIt ei;
	edgeSet M;
	M.clear();
	unsigned long tnow[3],*adjNow;
	unsigned long i,j,tcode,numtris=(unsigned long)_tris.size();
	_adjs.clear();
	_adjs.assign(numtris,0x00000003);
	for(i=0; i<numtris; i+=3)
	{
		if(_tris[i] == 0xffffffff)	// signals a deleted triangle
			continue;
		for(j=0; j<3; j++)
			tnow[j] = _tris[i+j];
		adjNow = &(_adjs[i]);
		for(j=0; j<3; j++)
		{
			if(adjNow[j]!=0x00000003)	// adjacency already computed
				continue;
			unsigned int tmp;
			if((tmp=tnow[(j+1)%3])<tnow[j]) {
				E.vtxMin = tmp;
				E.vtxMax = tnow[j];
				E.reversed = 1;	}
			else {
				E.vtxMin = tnow[j];
				E.vtxMax = tmp;
				E.reversed = 0;	}
			E.adjCode = ((i/3)<<2) + j;
			P = M.insert(E);
			// if P.second is true, no match so edge inserted
			if(P.second==false)	// edge match found
			{
				if(P.first->reversed==E.reversed && E.vtxMin!=E.vtxMax)
					return false;
				tcode = P.first->adjCode;
				adjNow[j] = tcode;
				_adjs[(tcode>>2)*3+(tcode&0x00000003)] = E.adjCode;
				M.erase(P.first);
			}
			else
				adjNow[j] = 0x00000003;
		}
	}
	M.clear();
	makeVertexToTriangleMap();
	return true;
}

void skinGraphics::makeVertexToTriangleMap()
{
	int i,j,numtris=(int)_tris.size();
	unsigned int *tnow,vnow;
	// provide each vertex with a face it is a member of
	_vertexFace.clear();
	_vertexFace.assign(_xyz1.size()>>2,0x80000000);	// initially deleted
	for(i=0; i<numtris; i+=3)
	{
		tnow = &(_tris[i]);
		if(tnow==NULL || *tnow>0xfffffffe)	// signals a deleted triangle
			continue;
		for(j=0; j<3; j++)
		{
			vnow = tnow[j];
			if(_vertexFace[vnow]&0x40000000)
				continue;	// vertex first on free edge, don't change
			_vertexFace[vnow] = i/3;
			if(_adjs[i+j]==0x00000003)	// vertex first on free edge, lock it for easy neighbor find
				_vertexFace[vnow] |= 0x40000000;
		}
	}
}

void skinGraphics::getNeighbors(unsigned int vertex, std::vector<neighborNode> &neighbors)
{
	unsigned long normVert,trNum,adj,triStart;
	// ASSUMES adjacencies already computed!
	triStart = _vertexFace[vertex];
	neighbors.clear();
	if(triStart&0x80000000)	// unconnected vertex
		return;
	trNum = triStart&0x3fffffff;
	unsigned long *adjs;
	GLuint *tnow = &(_tris[(trNum<<1)+trNum]);
	assert(tnow[0]!=0xffffffff);	// deleted triangle
	int j;
	for(j=0; j<3; ++j)
	{
		normVert = tnow[j];
		if(normVert==vertex)
			break;
	}
	assert(j<3);
	adjs = &(_adjs[(trNum<<1)+trNum]);
	// set triStart to the end adjacency code for counterclockwise traversal
	neighborNode n;
	if(triStart&0x40000000)	// started on a free edge, will end on one
	{
		triStart = 0x00000003;
		n.triangle = -1;	// code for open ring neighbors
		n.vertex = tnow[(j+1)%3];
		neighbors.push_back(n);
	}
	else	// create adjacency code of starting edge
		triStart = (trNum<<2)+j;
	n.triangle = trNum;
	n.vertex = tnow[(j+2)%3];
	neighbors.push_back(n);
	adj = adjs[(j+2)%3];
	while(adj!=triStart)
	{
		n.triangle = adj>>2;
		tnow = &_tris[(n.triangle<<1)+n.triangle];
		adjs = &(_adjs[(n.triangle<<1)+n.triangle]);
		j = adj&0x00000003;
		n.vertex = tnow[(j+2)%3];
		neighbors.push_back(n);
		adj = adjs[(j+2)%3];
	}
}

void skinGraphics::draw(void) 
{
	//glBindVertexArray(_vertexArrayBufferObject);
	//glActiveTexture( GL_TEXTURE0);
	//glBindTexture(GL_TEXTURE_2D, _2DtextureBuffer);
	//glActiveTexture( GL_TEXTURE1);
	//glBindTexture(GL_TEXTURE_2D, _2DnormalBuffer);
	//glDrawElements(GL_TRIANGLES, (GLsizei)_tris.size(), GL_UNSIGNED_INT, 0);
    // Never unbind a GL_ARRAY_BUFFER or GL_ELEMENT_ARRAY_BUFFER inside an active vertex array buffer object
	//glBindVertexArray(0);
}

void skinGraphics::computeLocalBounds()
{
	float minMaxXYZ[6];
	minMaxXYZ[0]=minMaxXYZ[2]=minMaxXYZ[4]=1e30f;
	minMaxXYZ[1]=minMaxXYZ[3]=minMaxXYZ[5]=-1e30f;
	int i,n=(unsigned int)_xyz1.size()>>2;
	GLfloat *v=&(_xyz1[0]);
	for(i=0; i<n; ++i)	{
		if(*v<minMaxXYZ[0])
			minMaxXYZ[0]=*v;
		if(*v>minMaxXYZ[1])
			minMaxXYZ[1]=*v;
		++v;
		if(*v<minMaxXYZ[2])
			minMaxXYZ[2]=*v;
		if(*v>minMaxXYZ[3])
			minMaxXYZ[3]=*v;
		++v;
		if(*v<minMaxXYZ[4])
			minMaxXYZ[4]=*v;
		if(*v>minMaxXYZ[5])
			minMaxXYZ[5]=*v;
		++v; ++v;
	}
	_localCenter[0] = (minMaxXYZ[0]+minMaxXYZ[1])*0.5f;
	_localCenter[1] = (minMaxXYZ[2]+minMaxXYZ[3])*0.5f;
	_localCenter[2] = (minMaxXYZ[4]+minMaxXYZ[5])*0.5f;
	_radius = (_localCenter[0]-minMaxXYZ[0])*(_localCenter[0]-minMaxXYZ[0]);
	_radius += (_localCenter[1]-minMaxXYZ[2])*(_localCenter[1]-minMaxXYZ[2]);
	_radius += (_localCenter[2]-minMaxXYZ[4])*(_localCenter[2]-minMaxXYZ[4]);
	_radius = sqrt(_radius);
	_boundsComputed=true;
}

void skinGraphics::getJavascriptMap(std::vector<int>& map){
    int nv;
    const float *fp = _tuvw.getPositionArray(nv);
    map.reserve( nv + _doubledVerts.size() );
    for(int i=0; i<nv; i++) 
        map.push_back( i );
    for(int j = 0; j < _doubledVerts.size(); j++)
        map.push_back( _doubledVerts.at(j) );
}

void skinGraphics::remapJavascriptVertexData( const std::vector<float> &xyz_in, std::vector<float> &xyz_out)
{
    //std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

    int nv;
    const float *fp = _tuvw.getPositionArray(nv);
    if( xyz_in.size() != nv*3 ){
        std::cout << "FAILURE!!!!" << std::endl;
        return;
    }
    
    xyz_out.clear();
    int n=(int)_xyz1.size();

    //std::cout << n << std::endl;
    //std::cout << nv << std::endl;
    //std::cout << _doubledVerts.size() << std::endl;

//    for(int j = 0; j < _doubledVerts.size(); j++) {
//        std::cout << "J: " << j << ",   DV(J): " << _doubledVerts.at(j) << std::endl;
//    }

    
    xyz_out.reserve( (nv + _doubledVerts.size()) * 3 );
    for(int i=0; i<nv; i++) {
        xyz_out.push_back( xyz_in.at( i*3 + 0) );
        xyz_out.push_back( xyz_in.at( i*3 + 1) );
        xyz_out.push_back( xyz_in.at( i*3 + 2) );
    }

    for(int j = 0; j < _doubledVerts.size(); j++) {
        xyz_out.push_back( xyz_in.at( _doubledVerts.at(j) *3 + 0) );
        xyz_out.push_back( xyz_in.at( _doubledVerts.at(j) *3 + 1) );
        xyz_out.push_back( xyz_in.at( _doubledVerts.at(j) *3 + 2) );
    }

    //std::cout << xyz_out.size() << std::endl;

    //std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;   

}

void skinGraphics::getJavascriptData(std::vector<int> &tris, std::vector<float> &xyz, std::vector<float> &uv)
{
    tris.clear(); xyz.clear(); uv.clear();
	tris.reserve(_tris.size());
	GLuint *tp;
	for(int n=(int)_tris.size(),i=0; i<n; i+=3) {
		tp = &_tris[i];
		if(*tp>0xfffffffe)
			continue;
		tris.push_back(tp[0]); tris.push_back(tp[1]); tris.push_back(tp[2]);
	}
    uv.assign(_uv.begin(),_uv.end());
    int i,n=(int)_xyz1.size();
    xyz.reserve((n>>2)*3);
    for(i=0; i<n; i+=4) {
        xyz.push_back(_xyz1[i]);
        xyz.push_back(_xyz1[i+1]);
        xyz.push_back(_xyz1[i+2]);
    }
}

void skinGraphics::getJavascriptPositions(std::vector<float> &xyz)
{
    xyz.clear();
    int i,n=(int)_xyz1.size();
    xyz.reserve((n>>2)*3);
    for(i=0; i<n; i+=4) {
        xyz.push_back(_xyz1[i]);
        xyz.push_back(_xyz1[i+1]);
        xyz.push_back(_xyz1[i+2]);
    }
}

skinGraphics::skinGraphics(void) : _vertexArrayBufferObject(0), _computeTangents(true), _getUvScale(false), _meanTriangleEdgeLength(-1.0f)
{
	for(int i=0; i<5; ++i)
		_bufferObjects[i]=0;
	_coloredNotTextured=false;
	_type = UVW_TRIANGLES;
	_name.clear();
}

skinGraphics::~skinGraphics(void)
{
}
