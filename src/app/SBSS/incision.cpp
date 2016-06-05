#include "GraphicsUtils/boundingBox.h"
#include "incision.h"
#include <map>
#include <math.h>
#include <assert.h>
#include <deque>
#include <algorithm>
#include "GraphicsUtils/trianglesUVW.h"
#include "GraphicsUtils/lines.h"
#include "GraphicsUtils/wxGraphics.h"
#include "Mat3x3f.h"

incision::incision() : _endHcut(false)
{
	_prefEdgeLength=0.02f; // wild guess. Should set externally
}

incision::~incision()
{
}

void incision::addCornerTrianglePairs()
{ // biggest problem here is possible self Tout. This can invalidate one or more pairs.
	std::vector<trianglesUVW::neighborNode> nei;
	int n=(int)_triPairs.size(),i,j;
	// make sure triUVW::_vertexFace is correct!
	// remove any invalid pairs and fill adjacent triangle data
	for(i=0; i<n; ++i) {
		if((_tri->_vertexFace[_triPairs[i].topV]&0x80000000)>0) {
			_triPairs[i].topV = -1; // mark removed
			continue; }
		if((_tri->_vertexFace[_triPairs[i].botV]&0x80000000)>0) {
			_triPairs[i].topV = -1; // mark removed
			continue; }
		if(_triPairs[i].topCounterclockwise) {
			_tri->getNeighbors(_triPairs[i].topV,nei);
			_triPairs[i].topTri = nei.back().triangle;
			_triPairs[i].botTri = _tri->_vertexFace[_triPairs[i].botV]&0x7fffffff;
		}
		else {
			_triPairs[i].topTri = _tri->_vertexFace[_triPairs[i].topV]&0x7fffffff;
			_tri->getNeighbors(_triPairs[i].botV,nei);
			_triPairs[i].botTri = nei.back().triangle;
		}
	}
	// must trim contiguous pairs and remove duplicates
	for(i=0; i<n-1; ++i) {
		if(_triPairs[i].topV<0)
			continue;
		j=i+1;
		while(j<n) { // bubble sort, but small. Could do better
			if(_triPairs[j].topV<0) {
				++j; continue; }
			if(_triPairs[i].topV==_triPairs[j].topV || _triPairs[i].botV==_triPairs[j].botV)
				break;
			++j;
		}
		if(j<n) {
			_triPairs[i].topV = -1; // mark removed - duplicate
			continue; }
	}
	for(i=0; i<n-1; ++i) {
		if(_triPairs[i].topV<0)
			continue;
		j=i+1;
		while(j<n) { // bubble sort, but small. Could do better
			if(_triPairs[j].topV<0) {
				++j; continue; }
			if(_triPairs[j].topTri==_triPairs[i].topTri || _triPairs[i].botTri==_triPairs[j].botTri) { // adjacent pairs
				_triPairs[i].topLen = -1.0f; // mark for trimming
				_triPairs[j].topLen = -1.0f; } // mark for trimming
			++j;
		}
	}
	auto adjVertTri = [this,&nei](int v, bool cclkWise, bool forceSplit, float lenMax, int &adjV, int &adjTri) {
		Vec3f v0,v1;
		int *tr,edge;
		float param,len;
		_tri->getNeighbors(v,nei);
		assert(nei.front().triangle<0);
		if(cclkWise) {
			adjV = nei.front().vertex;
			adjTri = nei[1].triangle; }
		else {
			adjV = nei.back().vertex;
			adjTri = nei.back().triangle; }
		if(forceSplit) // force adjacent edge split
			param = 0.33f;
		else {
			_tri->getVertexCoordinate(v,v0._v);
			_tri->getVertexCoordinate(adjV,v1._v);
			len = (v0-v1).length();
			if(len>lenMax)
				param = lenMax/len;
			else
				param = -1.0f;
		}
		if(param>=0.0f) { // can't use existing adj must cut edge
			tr = _tri->triVerts(adjTri);
			for(edge=0; edge<3; ++edge)
				if(tr[edge]==v)
					break;
			assert(edge<3);
			if(!cclkWise) {
				edge = (edge+2)%3;
				param = 1.0f - param; }
			adjV = _tri->splitTriangleEdge(adjTri,edge,param);
			float *w = &_tri->_uvw[adjV*3+2]; // fix any interpolation error
			if(*w<1e-5f) *w = 0.0f;
			if(*w>0.99998f) *w = 1.0f;
			if(!cclkWise) {
				_tri->getNeighbors(v,nei);
				adjV = nei.back().vertex;
				adjTri = nei.back().triangle;
			}
		}
	};
	for(i=0; i<n; ++i) {
		if(_triPairs[i].topV<0) // discard
			continue;
		int topAdjV,topAdjTri,botAdjV,botAdjTri;
		int t[2],v[3],top[2],bot[2];
		adjVertTri(_triPairs[i].topV,_triPairs[i].topCounterclockwise,_triPairs[i].topLen<0.0f,_triPairs[i].topLen*0.1f,topAdjV,topAdjTri);
		adjVertTri(_triPairs[i].botV,!_triPairs[i].topCounterclockwise,_triPairs[i].topLen<0.0f,_triPairs[i].botLen*0.1f,botAdjV,botAdjTri);
		if(!_triPairs[i].topCounterclockwise) {
			top[1]=_triPairs[i].topV; top[0]=topAdjV; bot[1]=botAdjV; bot[0]=_triPairs[i].botV; }
		else {
			top[0]=_triPairs[i].topV; top[1]=topAdjV; bot[0]=botAdjV; bot[1]=_triPairs[i].botV; }
		v[0]=top[1]; v[1]=top[0]; v[2]=bot[0];
		t[0] = _tri->addTriangle(v);
		v[0]=bot[1]; v[1]=bot[0]; v[2]=top[0];
		t[1] = _tri->addTriangle(v);
		unsigned long TE;
		int *tr=_tri->triVerts(topAdjTri);
		for(j=0; j<3; ++j) {
			if(tr[j]==top[0])
				break;
		}
		assert(j<3);
		_tri->_adjs[topAdjTri*3+j] = t[0]<<2;
		TE = topAdjTri<<2;
		TE+=j;
		unsigned long *adj=&_tri->_adjs[t[0]*3];
		adj[0] = TE;
		adj[2] = 0x03;
		adj[1] = (t[1]<<2)+1;

		tr=_tri->triVerts(botAdjTri);
		for(j=0; j<3; ++j) {
			if(tr[j]==bot[0])
				break;
		}
		assert(j<3);
		assert(tr[(j+1)%3]==bot[1]);
		_tri->_adjs[botAdjTri*3+j] = t[1]<<2;
		TE = botAdjTri<<2;
		TE+=j;
		adj=&_tri->_adjs[t[1]*3];
		adj[0] = TE;
		adj[2] = 0x03;
		adj[1] = (t[0]<<2)+1;
		_tri->_vertexFace[top[0]] = t[1]|0x40000000;
		_tri->_vertexFace[bot[0]] = t[0]|0x40000000;
	}
}

void incision::cleanCutEdge(std::set<trianglesUVW::edge,trianglesUVW::edgeTest> &borderEdges, float minSize)
{
	std::set<trianglesUVW::edge,trianglesUVW::edgeTest>::iterator bit;
	std::set<unsigned long> done;
	auto cleanVertex = [this,&minSize](unsigned long v) {
		if(_tri->_vertexFace[v]&0x80000000)
			return;
		int triangle = _tri->_vertexFace[v];
		assert(triangle&0x40000000);
		triangle &= 0x3fffffff;
		int *tr=&_tri->_tris[triangle*3];
		if(*tr<0) return;
		Vec3f last,now;
		int j;
		for(j=0; j<3; ++j)
			if(tr[j]==v)
				break;
		assert(j<3);
		_tri->getVertexCoordinate(tr[j],now._v);
		_tri->getVertexCoordinate(tr[(j+1)%3],last._v);
		assert(_tri->_vertexFace[tr[(j+1)%3]]&0x40000000);
		last -=now;
		int vRemaining=2;
		if(fabs(last._v[0])<minSize && fabs(last._v[1])<minSize && fabs(last._v[2])<minSize){ // kill last edge
			if(_corners.find(tr[j])==_corners.end())
				vRemaining=1;
			else if(_corners.find(tr[(j+1)%3])==_corners.end())
				vRemaining=0;
			else ;
			if(vRemaining<2)
				_tri->deleteEdge(triangle,j,vRemaining);
		}
	};
	for(bit=borderEdges.begin(); bit!=borderEdges.end(); ++bit) {
		if(done.insert(bit->vtxMin).second==true)
			cleanVertex(bit->vtxMin);
		if(done.insert(bit->vtxMax).second==true)
			cleanVertex(bit->vtxMax);
	}
}

void incision::midEdgeLocus(const int v0, const int v1, int &midTri, Vec3f &midPoint)
{
	Vec3f p,q,now,last;
	_tri->getVertexCoordinate(v0,p._v);
	_tri->getVertexCoordinate(v1,q._v);
	q -= p;
	last = p;
	float dSq=q*q,t,tLast=0.0f;
	assert(dSq>1e-16f);
	q /= dSq;
	std::vector<trianglesUVW::neighborNode> nei;
	int vLast=v0,vNext=-1;
	while(vNext!=v1) {
		_tri->getNeighbors(vLast,nei);
		if(nei.front().triangle<0) {
			vNext = nei.front().vertex;
			midTri = nei[1].triangle; }
		else {
			std::vector<trianglesUVW::neighborNode>::iterator nit,bit=nei.end();
			float maxN=-1e30f;
			for(nit=nei.begin(); nit!=nei.end(); ++nit) {
				_tri->getVertexCoordinate(nit->vertex,now._v);
				now -= last;
				now.normalize();
				t = now*q;
				if(t>maxN){
					maxN=t;
					bit=nit;
				}
			}
			vNext = bit->vertex;
			++bit;
			if(bit==nei.end())
				midTri = nei.front().triangle;
			else
				midTri = bit->triangle;
		}
		_tri->getVertexCoordinate(vNext,now._v);
		now -= p;
		t = now*q;
		if(t>=0.5f)
			break;
		last = now;
		vLast = vNext;
		tLast = t;
	}
//	t = (0.5f-tLast)/(t-tLast); // 0 div not possible
//	midPoint = now*t + last*(1.0f-t);
	midPoint.set(0.0f,0.0f,0.0f);
	for(int i=0; i<3; ++i) {
		_tri->getVertexCoordinate(_tri->_tris[midTri*3+i],now._v);
		midPoint += now;
	}
	midPoint /= 3.0f;
}

bool incision::Tin(const Vec3f &startPos, const Vec3f &startNorm, const Vec3f &nextPos, const Vec3f &nextNorm, int punchStatus, unsigned long (&Tvertices)[4])
{ // Input punchStatus=0 means no punch done, 1 punch+Tin, 2 punch+Tout
	// condition start and next to assure they are valid

#ifdef _DEBUG
	std::vector<GLuint> lines;
	std::vector<GLfloat> points;
	auto addSegment = [&](int v0, int v1) {
		float *vp = _tri->vertexCoordinate(v0);
		points.push_back(vp[0]);
		points.push_back(vp[1]);
		points.push_back(vp[2]);
		points.push_back(1.0f);
		lines.push_back((int)(points.size()>>2)-1);
		vp = _tri->vertexCoordinate(v1);
		points.push_back(vp[0]);
		points.push_back(vp[1]);
		points.push_back(vp[2]);
		points.push_back(1.0f);
		lines.push_back((int)(points.size()>>2)-1);
		lines.push_back(0xffffffff);
	};
	auto addSegmentPoints = [&](float (&v0)[3], float (&v1)[3]) {
		points.push_back(v0[0]);
		points.push_back(v0[1]);
		points.push_back(v0[2]);
		points.push_back(1.0f);
		lines.push_back(((int)points.size()>>2)-1);
		points.push_back(v1[0]);
		points.push_back(v1[1]);
		points.push_back(v1[2]);
		points.push_back(1.0f);
		lines.push_back(((int)points.size()>>2)-1);
		lines.push_back(0xffffffff);
	};
//	return false;
#endif
	
	int triangle=-1,edge,oppTri,oppEdge,nextTri,nextOppTri;
	float param,oppParam;
	Vec3f startP=startPos,nextP=nextPos,oppP,nextOppP;
	_tri->getNearestHardEdge(startP._v,triangle,edge,param,1);	// Input first 2 args, then overwrites all 4 with the point on the nearest hard edge.
	Vec3f norm=startNorm+nextNorm,negNorm;
	negNorm = -norm;
	if(punchStatus<1) { // no prior punch
		float nextUv[2]={0.0f,0.0f},nextOppUv[2]={0.0f,0.0f};
		if(!getTopBottomHits(nextPos,negNorm,nextTri,nextUv,nextOppTri,nextOppUv))
			return false;
		_tri->getBarycentricPosition(nextTri,nextUv,nextP._v);
		_tri->getBarycentricPosition(nextOppTri,nextOppUv,nextOppP._v);
	}
	else if(punchStatus<2) { // punch+Tin
		midEdgeLocus(_frontPunchVerts[0],_frontPunchVerts[1],nextTri,nextP);
		midEdgeLocus(_frontPunchVerts[3],_frontPunchVerts[2],nextOppTri,nextOppP);
	}
	else { // punch+Tout
		midEdgeLocus(_backPunchVerts[0],_backPunchVerts[1],nextTri,nextP);
		midEdgeLocus(_backPunchVerts[3],_backPunchVerts[2],nextOppTri,nextOppP);
	}
	if(!topoPath2(nextP,nextTri,startP,triangle,norm,edge,param))
		return false;	// COURT - write me.
	int botStart,topStart;
	if(param<0.001f)
		topStart = _tri->_tris[triangle*3+edge];
	else if(param>0.999f)
		topStart = _tri->_tris[triangle*3+((edge+1)%3)];
	else
		topStart=_tri->splitTriangleEdge(triangle,edge,param);
	_tri->getVertexCoordinate(topStart,startP._v);
	// edge normal is in the next triangle beyond top edge
	unsigned long adj = _tri->_adjs[triangle*3+edge];
	if(adj!=0x03)
		_tri->getTriangleNormal(adj>>2,_TEdgeNormal); // will need this to cut closed edge segment
	Vec3f planeN = (startP-nextPos)^norm;
	planeN.normalize();
	float planeD = planeN*startP;
	if(!getOppositeIncisionEdge2(triangle,edge,startP,negNorm,oppTri,oppEdge,oppParam,oppP))
			return false;	// need to debug. Ask user to save history.
	if(!topoPath2(nextOppP,nextOppTri,oppP,oppTri,negNorm,oppEdge,oppParam))
		return false;	// COURT - write me.
	if(oppParam<0.001f)
		botStart = _tri->_tris[oppTri*3+oppEdge];
	else if(param>0.999f)
		botStart = _tri->_tris[oppTri*3+((oppEdge+1)%3)];
	else
		botStart = _tri->splitTriangleEdge(oppTri,oppEdge,oppParam);

	
//	Vec3f junkP = oppP + negNorm;
//	addSegmentPoints(oppP._v,junkP._v);
//	addSegment(topStart,botStart);
//	float color[4]={1.0f,1.0f,0.0f,1.0f};
//	_wgg->getLines()->setColor(color);
//	_wgg->getLines()->addLines(points,lines);
//	return true;



	// now have plane eq. and top and bottom start points
	unsigned long TE[4];
	float prm[4];
	float leftD=planeD-_incisionWidthOver2,rightD=planeD+_incisionWidthOver2,outD;
	outD=leftD;
	bool leftMod=false,rightMod=false;
	if(!edgeBound(topStart,planeN,true,rightD,outD,TE[1],prm[1])) {
		if(prm[1]>1.5f) return false;
		leftMod = true;
		leftD = outD; }
	if(!edgeBound(botStart,planeN,false,rightD,outD,TE[3],prm[3])) {
		if(prm[3]>1.5f) return false;
		edgeBound(topStart,planeN,true,rightD,outD,TE[1],prm[1]);
		if(prm[1]>1.5f) return false;
		leftMod = true;
		leftD = outD; }
	if(leftMod)
		rightD = leftD + 2.0f*_incisionWidthOver2;
	outD = rightD;
	if(!edgeBound(topStart,planeN,false,leftD,outD,TE[0],prm[0])) {
		if(prm[0]>1.5f) return false;
		rightMod = true;
		rightD = outD; }
	if(!edgeBound(botStart,planeN,true,leftD,outD,TE[2],prm[2])) {
		if(prm[2]>1.5f) return false;
		edgeBound(topStart,planeN,false,leftD,outD,TE[0],prm[0]);
		if(prm[0]>1.5f) return false;
		rightMod = true;
		rightD = outD; }
	if(rightMod && !leftMod) { // look for new left side otherwise accept impossible situation
		leftD = rightD - 2.0f*_incisionWidthOver2;
		outD = leftD;
		// repeat left search with new limit
		edgeBound(topStart,planeN,true,rightD,outD,TE[1],prm[1]);
		if(prm[1]>1.5f) return false;
		if(!edgeBound(botStart,planeN,false,rightD,outD,TE[3],prm[3])) {
			if(prm[3]>1.5f) return false;
			edgeBound(topStart,planeN,true,rightD,outD,TE[1],prm[1]);
		}
	}
	// incision edge answers in TE and prm
	for(int i=0; i<4; ++i) {
		if(prm[i]<-0.5f)
			Tvertices[i] = TE[i];
		else {
			if(!(i&0x01) && TE[i]==TE[i+1]) {
				if(prm[i]>prm[i+1])
					prm[i+1] /= prm[i]; // zero divisor not possible
				else {
					prm[i] /= prm[i+1]; // zero divisor not possible
					if(prm[i+1]>0.98f)
						Tvertices[i+1] = _tri->_tris[(TE[i+1]>>2)*3+(((TE[i+1]&3)+1)%3)];
					else if(prm[i+1]>0.02f)
						Tvertices[i] = _tri->splitTriangleEdge(TE[i+1]>>2,TE[i+1]&0x03,prm[i+1]);
					else  // prm[i+1]<0.02f not possible if prm[i] smaller with same TE
						assert(false);
					if(prm[i]<0.02f)
						Tvertices[i] = _tri->_tris[(TE[i]>>2)*3+(TE[i]&3)];
					else if(prm[i]<0.98f)
						Tvertices[i] = _tri->splitTriangleEdge(TE[i]>>2,TE[i]&0x03,prm[i]);
					else  // not possible if prm[i+1] > prm[i]
						Tvertices[i] = _tri->_tris[(TE[i]>>2)*3+(((TE[i]&3)+1)%3)];
					++i; continue;
				}
			}
			if(prm[i]<0.02f) {
				Tvertices[i] = _tri->_tris[(TE[i]>>2)*3+(TE[i]&3)];
				continue; }
			else if(prm[i]>0.98f) {
				Tvertices[i] = _tri->_tris[(TE[i]>>2)*3+(((TE[i]&3)+1)%3)];
				continue; }
			else
				Tvertices[i] = _tri->splitTriangleEdge(TE[i]>>2,TE[i]&0x03,prm[i]);
		}
	}

/*	std::vector<GLfloat> points;
	std::vector<GLuint> lines;
	float xyz[3];
	points.reserve(16);
	_tri->getVertexCoordinate(Tvertices[0],xyz);
	points.push_back(xyz[0]);	points.push_back(xyz[1]);	points.push_back(xyz[2]);	points.push_back(1.0f);
	_tri->getVertexCoordinate(Tvertices[1],xyz);
	points.push_back(xyz[0]);	points.push_back(xyz[1]);	points.push_back(xyz[2]);	points.push_back(1.0f);
	_tri->getVertexCoordinate(Tvertices[3],xyz);
	points.push_back(xyz[0]);	points.push_back(xyz[1]);	points.push_back(xyz[2]);	points.push_back(1.0f);
	_tri->getVertexCoordinate(Tvertices[2],xyz);
	points.push_back(xyz[0]);	points.push_back(xyz[1]);	points.push_back(xyz[2]);	points.push_back(1.0f);
	lines.reserve(5);
	lines.push_back(0);	lines.push_back(1);	lines.push_back(2);	lines.push_back(3);	lines.push_back(0);	
	_wgg->getLines()->addLines(points,lines);
	return false; */

	return true;
}

bool incision::edgeBound(int &vStart, const Vec3f &planeN, bool neiFront, float limitD, float &targetD, unsigned long &triEdge, float &param)
{ // inputs a start vertex v which must be on an edge. In getting  hardEdgeNeighbors of v, if neiFront search in that direction, otherwise
	// search in neiBack direction. planeN is normalized plane normal of T cut. Exceeding limitD triggers false setting param to -1
	// and triEdge to the VERTEX closest to targetD with target D returning its D value. If returns true, targetD was found with
	// result returned in triEdge and param.
	param = -1; // set for false return
	triEdge = vStart;
	int vLast,*tv;
	Vec3f pos;
	float dLast,dNext,dMax,beyondLimit;
	std::vector<trianglesUVW::neighborNode> nei;
	_tri->getVertexCoordinate(vStart,pos._v);
	dNext = planeN*pos;
	bool movingDown;
	auto slideStart = [&](bool doTarget) ->bool{
		int i;
		for(i=0; i<5; ++i) {
			if(!getHardEdgeNeighbors(vStart,nei))
				return false;
			if(neiFront==doTarget)
				vStart = nei.back().vertex;
			else
				vStart = nei.front().vertex;
			_tri->getVertexCoordinate(vStart,pos._v);
			dNext = planeN*pos;
			if(doTarget && ((movingDown && targetD<dNext) || (!movingDown && targetD>dNext)))
				break;
			if(!doTarget && ((movingDown && limitD>dNext) || (!movingDown && limitD<dNext)))
				break;
		}
		if(i>4) {
			param = 2.0f;
			return false; }
		dMax = dNext;
		return true;
	};
	int v;
	auto edgeStraddle = [&]() {
		if(neiFront)
			triEdge = nei[1].triangle;
		else
			triEdge = nei.back().triangle;
		tv = _tri->triVerts(triEdge);
		int i;
		for(i=0; i<3; ++i) {
			if(neiFront) {
				param = targetD-dLast;
				if(tv[i]==vLast) break;
			}
			else {
				param = dNext - targetD;
				if(tv[i]==v) break;
			}
		}
		assert((neiFront && tv[(i+1)%3]==v) || (!neiFront && tv[(i+1)%3]==vLast));
		triEdge <<= 2;   triEdge += i;
		param /= (dNext-dLast); // denom never 0
	};
	dMax = dNext;
	if(limitD<targetD) {
		movingDown = false;
		if(targetD<dNext)
			if(!slideStart(true))
				return false;
		if(limitD>dNext)
			beyondLimit = dNext;
//			if(!slideStart(false))
	//			return false;
		else
			beyondLimit = 1e30f;
	}
	else if(limitD>targetD) {
		movingDown = true;
		if(targetD>dNext)
			if(!slideStart(true))
				return false;
		if(limitD<dNext)
			beyondLimit = dNext;
//			if(!slideStart(false))
//				return false;
		else
			beyondLimit = -1e30f;
	}
	else {
		param = 2.0f;
		return false;
	}
	v = vStart;
	while(true)	{
		if(!getHardEdgeNeighbors(v,nei))
			return false;
		vLast = v;
		dLast = dNext;
		if(neiFront)
			v = nei.front().vertex;
		else
			v = nei.back().vertex;
		_tri->getVertexCoordinate(v,pos._v);
		dNext = planeN*pos;
		if(movingDown) {
			if(limitD<dNext) {
				if(dNext<beyondLimit)
					beyondLimit = dNext;
				else {
					targetD = dMax;
					return false; }
			}
			if(dNext<dMax) {
				dMax=dNext;
				triEdge = v; }
			if(dNext<targetD) {
				edgeStraddle();
				break;
			}
		}
		else {
			if(limitD>dNext) {
				if(dNext>beyondLimit)
					beyondLimit = dNext;
				else {
					targetD = dMax;
					return false; }
			}
			if(dNext>dMax) {
				dMax=dNext;
				triEdge = v; }
			if(dNext>targetD) {
				edgeStraddle();
				break;
			}
		}
	}
	return true;
}

void incision::subdivideLongTriangles(float maxSize, int startTriangle)
{
	assert(maxSize>1e-5f);
	Vec3f last,now;
	float wLast,wNow;
	float len,msSq=maxSize*maxSize,msSmall=0.98f*maxSize;
	startTriangle *= 3;
	for(int j,i=startTriangle; i<(int)_tri->_tris.size(); i+=3) { // routine makes _tris bigger while it runs
		_tri->getVertexCoordinate(_tri->_tris[i+2],last._v);
		wLast = _tri->vertexTexture(_tri->_tris[i+2])[2];
		for(j=0; j<3; ++j) {
			_tri->getVertexCoordinate(_tri->_tris[i+j],now._v);
			wNow = _tri->vertexTexture(_tri->_tris[i+j])[2];
			if(fabs(wLast-wNow)<1e-5f) { // only subdivide middle spanning edges
				last = now;
				wLast = wNow;
				continue; }
			last -= now;
			if((len=last.length2())>msSq) {
				len = sqrt(len);
				len = msSmall/len;
				_tri->splitTriangleEdge(i/3,(j+2)%3,len); // makes new triangle(s) so _tris.size() increases
				i-=3; // redo tri as other edges may be too long
				break;
			}
			last = now;
			wLast = wNow;
		}
	}

	int junk =(int)_tri->_tris.size();
	junk = junk;
}

void incision::fillHoles()
{ // assumes topology arrays are correct
	std::vector<trianglesUVW::neighborNode> nei;
	vert_triEdge vte;
	std::list<vert_triEdge> hole;
	int subdivStart=(int)_tri->_tris.size();
	auto fillHole = [this](std::list<vert_triEdge> &hole) {
		std::list<vert_triEdge>::iterator hit0,hit1,hit2=hole.begin();
		float w0=-2.0f,w1=-2.0f,w2;
		while(true) {
			w2 = _tri->_uvw[hit2->vertex*3+2];
			if(w0>0.9999f && w1>0.9999f && w2<0.998f)
				break;
			w0=w1; w1=w2;
			hit0=hit1; hit1=hit2;
			++hit2;
			if(hit2==hole.end())
				hit2=hole.begin();
		}
		auto addBackTriangle = [&](trianglesUVW *tUVW) {
			int tv[3];
			tv[0]=hit1->vertex; tv[1]=hit0->vertex; tv[2]=hit2->vertex;
			int newT = tUVW->addTriangle(tv);
			tUVW->triAdjs(hit0->triEdge>>2)[hit0->triEdge&0x03] = newT<<2;
			tUVW->triAdjs(hit1->triEdge>>2)[hit1->triEdge&0x03] = (newT<<2)+2;
			unsigned long *adjs=tUVW->triAdjs(newT);
			adjs[0] = hit0->triEdge;
			adjs[1] = 0x03;
			adjs[2] = hit1->triEdge;
			*tUVW->vertexTriangle(hit0->vertex) = newT|0x40000000;
			*tUVW->vertexTriangle(hit1->vertex) = newT; // hit2 vf doesn't change
			hit0->triEdge = (newT<<2)+1;
			hole.erase(hit1); // unnecessary - just for diagnostics
		};
		addBackTriangle(_tri);
		while(w2>0.0001f) {
			hit1=hit2; ++hit2;
			if(hit2==hole.end())
				hit2=hole.begin();
			w2 = _tri->vertexTexture(hit2->vertex)[2];
			addBackTriangle(_tri);
		}
		float d0=-1.0f,d2=-1.0f;
		std::list<vert_triEdge>::iterator next0,next2,next00,next22;
		Vec3f p,q;
		char run=0x03;
		while(run) {
			if(d0<0.0f) {
				next0=hit0;
				if(next0==hole.begin())
					next0=hole.end();
				--next0;
				next00=next0;
				if(next00==hole.begin())
					next00=hole.end();
				--next00;
				if(_tri->vertexTexture(next00->vertex)[2]<0.9999f) {
					d0 = 4e30f;
					run &= 0x02; }
				if(!run)
					break;
			}
			if(d2<0.0f) {
				next2=hit2; ++next2;
				if(next2==hole.end()) next2=hole.begin();
				next22=next2; ++next22;
				if(next22==hole.end()) next22=hole.begin();
				if(_tri->vertexTexture(next22->vertex)[2]>0.0001f) {
					d2 = 4e30f;
					run &= 0x01; }
				if(!run)
					break;
			}
			if(d0<1e30f) {
				p.set((const float (&)[3])*_tri->vertexCoordinate(hit2->vertex));
				q.set((const float (&)[3])*_tri->vertexCoordinate(next0->vertex));
				d0 = (p-q).length2(); }
			if(d2<1e30f) {
				p.set((const float (&)[3])*_tri->vertexCoordinate(hit0->vertex));
				q.set((const float (&)[3])*_tri->vertexCoordinate(next2->vertex));
				d2 = (p-q).length2(); }
			if(d0<d2) {
				hit1 = hit0;
				hit0 = next0;
				d0 = -1.0f;
			}
			else {
				hit1 = hit2;
				hit2 = next2;
				d2 = -1.0f;
			}
			addBackTriangle(_tri);
		}
		do {
			hit1=hit2; ++hit2;
			if(hit2==hole.end())
				hit2=hole.begin();
			addBackTriangle(_tri);
		}while(hit2!=next0);
		_tri->_adjs[(hit0->triEdge>>2)*3+(hit0->triEdge&0x03)] = hit2->triEdge;
		_tri->_adjs[(hit2->triEdge>>2)*3+(hit2->triEdge&0x03)] = hit0->triEdge;
		_tri->_vertexFace[hit0->vertex] = hit0->triEdge>>2;
		_tri->_vertexFace[hit2->vertex] = hit2->triEdge>>2;
	};
	int n=(int)_tri->_adjs.size();
	unsigned long *adjs;
	int *tri;
	for(int k,j,i=0; i<n; i+=3) {
		if(_tri->_tris[i]<0)
			continue;
		adjs = &_tri->_adjs[i];
		for(j=0; j<3; ++j) {
			if(adjs[j]==0x03) { // construct hole and fill it
				int startV = _tri->_tris[i+j];
				vte.vertex = startV;
				vte.triEdge = ((i/3)<<2)+j;
				hole.push_back(vte);
				vte.vertex = _tri->_tris[i+((j+1)%3)];
				while(vte.vertex!=startV) {
					_tri->getNeighbors(vte.vertex,nei);
					assert(nei.front().triangle<0);
					vte.triEdge = nei[1].triangle;
					tri = &_tri->_tris[vte.triEdge*3];
					for(k=0; k<3; ++k)
						if(tri[k]==vte.vertex)
							break;
					vte.triEdge <<= 2;
					vte.triEdge += k;
					hole.push_back(vte);
					vte.vertex = nei.front().vertex;
				}
				fillHole(hole);
				hole.clear();
			}
		}
	}
}

bool incision::cutTborders(unsigned long (&p)[4], unsigned long (&q)[4], bool downQ, std::set<trianglesUVW::edge,trianglesUVW::edgeTest> &borderEdges)
{	// p and q are 4 vertex sets that define the two ends of a T segment to be removed.
	// Looking from center of T segment to either end, first 2 have w==1.0, second 2 w==0.0, left to right.
	// Collects borderEdges for subsequent recursive delete of interior. downQ determines whether cut extends down q end
	std::vector<unsigned int> path;
	Vec3f nrm;
	auto vNorm = [this](unsigned int v, Vec3f &normal, int w) { // w==2 mean a side edge normal is desired
		std::vector<trianglesUVW::neighborNode> nei;
		std::vector<trianglesUVW::neighborNode>::iterator nit;
		_tri->getNeighbors(v,nei);
		normal.set(0.0f,0.0f,0.0f);
		Vec3f center,last,now;
		float *fp=&_tri->_xyz[v*3];
		center.set(fp[0],fp[1],fp[2]);
		assert(!nei.empty());
		float wNow,wLast;
		nit=nei.begin();
		if(nit->triangle<0) {
			fp = &_tri->_xyz[nit->vertex*3];
			last.set(fp[0],fp[1],fp[2]);
			wLast=_tri->_uvw[nit->vertex*3+2];
			++nit; }
		else {
			wLast=_tri->_uvw[nei.back().vertex*3+2];
			fp = &_tri->_xyz[nei.back().vertex*3];
			last.set(fp[0],fp[1],fp[2]); }
		last -= center;
		bool addIt;
		for(; nit!=nei.end(); ++nit) {
			addIt=true;
			wNow = _tri->_uvw[nit->vertex*3+2];
			fp = &_tri->_xyz[nit->vertex*3];
			now.set(fp[0],fp[1],fp[2]);
			now -= center;
			if(w>1 && fabs(wNow-wLast)<1e-5f)
				addIt=false;
			else if(w<1 && (wNow>1e-5f || wLast>1e-5f))
				addIt=false;
			if(w==1 && (wNow<0.99998f || wLast<0.99998f))
				addIt=false;
			if(addIt)
				normal += last^now;
			wLast = wNow;
			last = now;
		}
		normal.normalize();
	};
	Vec3f n0,n1;

/* #ifdef _DEBUG
	std::vector<GLuint> lines;
	std::vector<GLfloat> points;
	auto addSegment = [&](int v0, int v1) {
		float *vp = _tri->vertexCoordinate(v0);
		points.push_back(vp[0]);
		points.push_back(vp[1]);
		points.push_back(vp[2]);
		points.push_back(1.0f);
		lines.push_back((int)(points.size()>>2)-1);
		vp = _tri->vertexCoordinate(v1);
		points.push_back(vp[0]);
		points.push_back(vp[1]);
		points.push_back(vp[2]);
		points.push_back(1.0f);
		lines.push_back((int)(points.size()>>2)-1);
		lines.push_back(0xffffffff);
	};
#endif */

	for(int j=1,i=0; i<4; ++i) {
		vNorm(p[i],n0,i<2?1:0);
		vNorm(q[j],n1,i<2?1:0);
		path.push_back(p[i]);
		_corners.insert(p[i]);
		_corners.insert(q[i]);

//		addSegment(p[i],q[j]);

		if(!segmentPath(p[i],q[j],n0+n1,path))
			return false;
		addBorderEdges(path,borderEdges);
		for(int num=(int)path.size(),k=0; k<num; ++k) { // remove interpolation artifact from w texture
			if(i<2)
				_tri->_uvw[path[k]*3+2] = 1.0f;
			else
				_tri->_uvw[path[k]*3+2] = 0.0f; }
		path.clear();
		if(i==1)
			j=3;
		else 
			--j;
	}
//	auto addDownBorders = [addSegment,this,&vNorm](unsigned int vTop, unsigned int vBottom, std::set<trianglesUVW::edge,trianglesUVW::edgeTest> &borderEdges) ->bool {
	auto addDownBorders = [this,&vNorm](unsigned int vTop, unsigned int vBottom, std::set<trianglesUVW::edge,trianglesUVW::edgeTest> &borderEdges) ->bool {
		std::vector<trianglesUVW::neighborNode> nei;
		_tri->getNeighbors(vTop,nei);
		if(nei.front().triangle<0)
			return true;
#ifdef _DEBUG
		_tri->getNeighbors(vBottom,nei);
		assert(nei.front().triangle>-1);
#endif
		std::vector<unsigned int> path;
		Vec3f n0,n1;
		vNorm(vTop,n0,2);
		vNorm(vBottom,n1,2);
		path.clear();
		path.push_back(vTop);

//		addSegment(vTop,vBottom);

		if(!segmentPath(vTop,vBottom,n0+n1,path))
			return false;
		addBorderEdges(path,borderEdges);
		return true;
	};
	for(int i=0; i<2; ++i) {
		if(!addDownBorders(p[i],p[i+2],borderEdges))
			return false;
		if(downQ && !addDownBorders(q[i],q[i+2],borderEdges))
			return false;
	}

//	_wgg->getLines()->addLines(points,lines);

	return true;
}

bool incision::addBorderEdges(std::vector<unsigned int> &path, std::set<trianglesUVW::edge,trianglesUVW::edgeTest> &borderEdges)
{
	if(path.empty())
		return true;
	unsigned int lastV=path.front();
	trianglesUVW::edge edg;
	for(int n=(int)path.size(),i=1; i<n; ++i) {
		edg.vtxMin = path[i];
		if(lastV<(int)edg.vtxMin)	{
			edg.vtxMax=edg.vtxMin;
			edg.vtxMin=lastV;
			lastV = edg.vtxMax;	}
		else	{
			edg.vtxMax=lastV;
			lastV = edg.vtxMin;	}
		if(!borderEdges.insert(edg).second)
			return false;
	}
	return true;
}

bool incision::punch(std::vector<Vec3f> &positions, std::vector<Vec3f> &normals, bool startT, bool endT)
{
	std::vector<int> topVerts,bottomVerts;
	if(!getCookieCutter2(positions,normals,startT,endT,topVerts,bottomVerts) )
		return false;
	_triPairs.clear();
	_corners.insert(topVerts.begin(),topVerts.end());
	_corners.insert(bottomVerts.begin(),bottomVerts.end());
	_frontPunchVerts[0]=topVerts[0];
	_frontPunchVerts[1]=topVerts[1];
	_frontPunchVerts[2]=bottomVerts[0];
	_frontPunchVerts[3]=bottomVerts[1];
	int lastTop,lastBot,n=(int)topVerts.size(),i;
	_backPunchVerts[0]=topVerts[n-1];
	_backPunchVerts[1]=topVerts[n-2];
	_backPunchVerts[2]=bottomVerts[n-1];
	_backPunchVerts[3]=bottomVerts[n-2];
	std::vector<unsigned int> topPath,botPath;
	if(!segmentPath(topVerts[1],topVerts[0],normals[0],topPath))
		return false;
	topPath.back() |= 0x80000000;
	if(!segmentPath(bottomVerts[1],bottomVerts[0],normals[0],botPath))
		return false;
	botPath.back() |= 0x80000000;
	lastTop = topVerts[0];
	lastBot = bottomVerts[0];
	Vec3f nrm;
	n >>=1;
	for(i=1; i<n; ++i) {
		nrm = normals[i-1] + normals[i];
		if(!segmentPath(lastTop,topVerts[i<<1],nrm,topPath))
			return false;
		topPath.back() |= 0x80000000;
		lastTop = topVerts[i<<1];
		if(!segmentPath(lastBot,bottomVerts[i<<1],nrm,botPath))
			return false;
		botPath.back() |= 0x80000000;
		lastBot = bottomVerts[i<<1];
	}
	if(!segmentPath(lastTop,topVerts[(n<<1)-1],normals[n-1],topPath))
		return false;
	topPath.back() |= 0x80000000;
	lastTop = topVerts[(n<<1)-1];
	if(!segmentPath(lastBot,bottomVerts[(n<<1)-1],normals[n-1],botPath))
		return false;
	botPath.back() |= 0x80000000;
	lastBot = bottomVerts[(n<<1)-1];
	for(i=n-2; i>-1; --i) {
		nrm = normals[i+1] + normals[i];
		if(!segmentPath(lastTop,topVerts[(i<<1)+1],nrm,topPath))
			return false;
		topPath.back() |= 0x80000000;
		lastTop = topVerts[(i<<1)+1];
		if(!segmentPath(lastBot,bottomVerts[(i<<1)+1],nrm,botPath))
			return false;
		botPath.back() |= 0x80000000;
		lastBot = bottomVerts[(i<<1)+1];
	}
	// topPath and botPath should now be topologically continuous
/* #ifdef _DEBUG
	auto topoCorral = [](trianglesUVW *tr, std::vector<unsigned int> &v) {
		std::vector<trianglesUVW::neighborNode> nei;
		std::vector<trianglesUVW::neighborNode>::iterator nit;
		unsigned int lastV=v.back()&0x7fffffff;
		for(int n=(int)v.size(),i=0; i<n; ++i) {
			tr->getNeighbors(v[i]&0x7fffffff,nei);
			for(nit=nei.begin(); nit!=nei.end(); ++nit) {
				if(nit->vertex==lastV)
					break;
			}
			if(nit==nei.end())
				i=i;
			lastV = v[i]&0x7fffffff;
		}
		lastV=0;
	};
	topoCorral(_tri,topPath);
	topoCorral(_tri,botPath);
#endif */
	auto triToGo = [](trianglesUVW *tr, std::vector<unsigned int> &v) -> int{
		int triGo=v[v.size()-2]&0x7fffffff;
		std::vector<trianglesUVW::neighborNode> nei;
		std::vector<trianglesUVW::neighborNode>::iterator nit,startit;
		tr->getNeighbors(v.back()&0x7fffffff,nei);
		if(nei.front().triangle<0)
			return -1;
		for(startit=nei.begin(); startit!=nei.end(); ++startit) {
			if(startit->vertex==triGo)
				break;
		}
		assert(startit!=nei.end());
		nit=startit; ++nit;
		float angle=0.0f;
		Vec3f center,lastV,nowV;
		center.set((const float (&)[3])tr->_xyz[(v.back()&0x7fffffff)*3]);
		lastV.set((const float (&)[3])tr->_xyz[triGo*3]);
		lastV -= center;
		lastV.normalize();
		while(nit!=startit){
			if(nit==nei.end())
				nit=nei.begin();
			nowV.set((const float (&)[3])tr->_xyz[nit->vertex*3]);
			nowV -= center;
			nowV.normalize();
			angle += acos(nowV*lastV);
			if(nit->vertex==(v.front()&0x7fffffff))
				break;
			lastV=nowV;
			++nit;
		}
		if(nit==startit)
			return -1;
		if(angle<3.1416f)
			return nit->triangle;
		else
			return startit->triangle;
	};
	// fix any interpolation errors in top and bottom w textures
	n = (int)topPath.size();
	for(i=0; i<n; ++i)
		_tri->_uvw[(topPath[i]&0x7fffffff)*3+2] = 1.0f;
	n = (int)botPath.size();
	for(i=0; i<n; ++i)
		_tri->_uvw[(botPath[i]&0x7fffffff)*3+2] = 0.0f;
	int triGo=triToGo(_tri,topPath);
	if(triGo<0)
		return false;
	deleteInteriorTriangles(triGo,topPath);
	triGo=triToGo(_tri,botPath);
	if(triGo<0)
		return false;
	deleteInteriorTriangles(triGo,botPath);
	_tri->makeVertexToTriangleMap();
	// add triangle pairs from top to bottom for later hole filler
	int nCorners=0,j=0;
	i=0;
	lastTop = (int)topVerts.size();
	triPairConn tpc;
	while(nCorners<lastTop) {
		while((topPath[i]&0x80000000)<1)
			++i;
		while((botPath[j]&0x80000000)<1)
			++j;
		++nCorners;
		if(nCorners==(lastTop>>1))
			tpc.topCounterclockwise = true;
		else if(nCorners<(lastTop>>1)+2)
			tpc.topCounterclockwise = false;
		else
			tpc.topCounterclockwise = true;
		tpc.topV = topPath[i]&0x7fffffff;
		tpc.botV = botPath[j]&0x7fffffff;
		_triPairs.push_back(tpc);
		++i; ++j;
	}
	Vec3f v0,v1;
	for(i=0; i<lastTop; ++i) {
		if(_triPairs[i].topCounterclockwise)
			j = i-1;
		else
			j = i+1;
		_tri->getVertexCoordinate(_triPairs[j].topV,v0._v);
		_tri->getVertexCoordinate(_triPairs[i].topV,v1._v);
		v1 -= v0;
		_triPairs[i].topLen = v1.length();
		_tri->getVertexCoordinate(_triPairs[j].botV,v0._v);
		_tri->getVertexCoordinate(_triPairs[i].botV,v1._v);
		v1 -= v0;
		_triPairs[i].botLen = v1.length();
	}
	return true;
}

void incision::trianglePairConnect(const unsigned int *topVerts, const unsigned int *botVerts)
{ // Inputs two top and bot verts. Connects top and bot with adjacent triangle pair returned in newTris.
	assert(false);

	std::vector<trianglesUVW::neighborNode> neiTop,neiBot;
	_tri->getNeighbors(topVerts[0]&0x7fffffff,neiTop);
	assert(neiTop.front().triangle<0);
	int top[2],bot[2];
	unsigned long topTE,botTE;
	if((topVerts[1]&0x7fffffff)==neiTop.front().vertex) {
		top[0] = topVerts[0]&0x7fffffff;
		top[1] = topVerts[1]&0x7fffffff;
		bot[1] = botVerts[0]&0x7fffffff;
		bot[0] = botVerts[1]&0x7fffffff;
		topTE = neiTop[1].triangle;
	}
	else {
		assert((topVerts[1]&0x7fffffff)==neiTop.back().vertex);
		top[0] = topVerts[1]&0x7fffffff;
		top[1] = topVerts[0]&0x7fffffff;
		bot[1] = botVerts[1]&0x7fffffff;
		bot[0] = botVerts[0]&0x7fffffff;
		topTE = neiTop.back().triangle;
	}
	_tri->getNeighbors(bot[0]&0x7fffffff,neiBot);
	botTE = neiBot[1].triangle;
	assert(neiBot.front().vertex==bot[1]);
	// don't bother comparing diagonal lengths.  Shouldn't be much skew.
	int t[2],v[3]={top[1],top[0],bot[0]};
	t[0] = _tri->addTriangle(v);
	v[0]=bot[1]; v[1]=bot[0]; v[2]=top[0];
	t[1] = _tri->addTriangle(v);
	int i,*tr=_tri->triVerts(topTE);
	for(i=0; i<3; ++i) {
		if(tr[i]==top[0])
			break;
	}
	assert(i<3);
	_tri->_adjs[topTE*3+i] = t[0]<<2;
	topTE<<=2; topTE+=i;
	unsigned long *adj=&_tri->_adjs[t[0]*3];
	adj[0] = topTE;
	adj[2] = 0x03;
	adj[1] = (t[1]<<2)+1;
	tr=_tri->triVerts(botTE);
	for(i=0; i<3; ++i) {
		if(tr[i]==bot[0])
			break;
	}
	assert(i<3);
	assert(tr[(i+1)%3]==bot[1]);
	_tri->_adjs[botTE*3+i] = t[1]<<2;
	botTE<<=2; botTE+=i;
	adj=&_tri->_adjs[t[1]*3];
	adj[0] = botTE;
	adj[2] = 0x03;
	adj[1] = (t[0]<<2)+1;
	_tri->_vertexFace[top[0]] = t[1]|0x40000000;
	_tri->_vertexFace[bot[0]] = t[0]|0x40000000;
}

void incision::deleteInteriorTriangles(long triangleToGo, std::vector<unsigned int> path)
{ // path normally closed. Bit 0x80000000 set means is a corner point. If path open, open span must be unconnected.
	std::set<trianglesUVW::edge,trianglesUVW::edgeTest> borderEdges;
	trianglesUVW::edge edg;
	unsigned int vNow,lastV=path.back()&0x7fffffff;
	for(int n=(int)path.size(),i=0; i<n; ++i) {
		vNow = path[i]&0x7fffffff;
		if(vNow<lastV) {
			edg.vtxMin=vNow; edg.vtxMax=lastV; }
		else {
			edg.vtxMax=vNow; edg.vtxMin=lastV; }
		borderEdges.insert(edg);
		lastV = vNow;
	}
	deleteInteriorTriangles(triangleToGo,borderEdges);
	_tri->makeVertexToTriangleMap();
	cleanCutEdge(borderEdges,_prefEdgeLength*0.2f);
}

bool incision::makeIncision(trianglesUVW* tri,int nPoints,float (*positions)[3],float (*normals)[3],int TstartTriangleEdge,float TstartParam,bool Tout)
{	// first 4 arguments are simple user inputs. TstartTriangleEdge if 0x03 means there is not a T in, otherwise there is a T in at the top triangle<<2|triangleEdge
	// identifying the top incision edge and TstartParam is the parameter (0-1) along that edge where T in commences. If TstartTriangleEdge<0 means is not a T in.
	// Tout is a simple bool since it is possible to self intersect.  In this case this routine makes the first part of the incision, then this routine finds the point
	// of T out from the last line segment recursively with an H cut.
	// Version 2.0 of the new incision code
	_tri = tri;
	_corners.clear();
	_triPairs.clear();
	int nPunch=nPoints;
	if(Tout) --nPunch;
	if(TstartTriangleEdge!=0x03) --nPunch;
	std::vector<Vec3f> posV,nrmV;
	posV.assign(nPoints,Vec3f());
	nrmV.assign(nPoints,Vec3f());
	for(int i=0; i<nPoints; ++i) {
		posV[i].set(positions[i]);
		nrmV[i].set(normals[i]);
	}
	if(nPunch>1) {
		if(!punch(posV,nrmV,TstartTriangleEdge!=0x03,Tout))
			return false;
	}
	if(nPunch==nPoints){
		int subdivStart = _tri->numberOfTriangles();
		// next two calls add triangles above subdivStart
		addCornerTrianglePairs();
		fillHoles();
//		subdivideLongTriangles(_prefEdgeLength*3.0f,subdivStart);
		_tri->cleanAndPack();
		_tri->_adjacenciesComputed=false;
		_tri->findAdjacentTriangles();
		return true;
	}
	// some sort of T occurs
	auto getDeleteTriangle = [this](unsigned int corner0, unsigned int corner1) ->int {
		std::vector<trianglesUVW::neighborNode> nei;
		if(!getHardEdgeNeighbors(corner0,nei))
			return -1;
		assert(!nei.empty());
		Vec3f c0,c1,v;
		_tri->getVertexCoordinate(corner0,c0._v);
		_tri->getVertexCoordinate(corner1,c1._v);
		c1 -= c0;
		_tri->getVertexCoordinate(nei.front().vertex,v._v);
		float dot = c1*(v-c0);
		_tri->getVertexCoordinate(nei.back().vertex,v._v);
		if(dot>c1*(v-c0))
			return nei[1].triangle;
		else
			return nei.back().triangle;
	};

/* #ifdef _DEBUG
	std::vector<GLuint> lines;
	std::vector<GLfloat> points;
	auto addSegment = [&](int v0, int v1) {
		float *vp = _tri->vertexCoordinate(v0);
		points.push_back(vp[0]);
		points.push_back(vp[1]);
		points.push_back(vp[2]);
		points.push_back(1.0f);
		lines.push_back(((int)points.size()>>2)-1);
		vp = _tri->vertexCoordinate(v1);
		points.push_back(vp[0]);
		points.push_back(vp[1]);
		points.push_back(vp[2]);
		points.push_back(1.0f);
		lines.push_back(((int)points.size()>>2)-1);
		lines.push_back(0xffffffff);
	};
	auto addSegmentPoints = [&](float (&v0)[3], float (&v1)[3]) {
		points.push_back(v0[0]);
		points.push_back(v0[1]);
		points.push_back(v0[2]);
		points.push_back(1.0f);
		lines.push_back(((int)points.size()>>2)-1);
		points.push_back(v1[0]);
		points.push_back(v1[1]);
		points.push_back(v1[2]);
		points.push_back(1.0f);
		lines.push_back(((int)points.size()>>2)-1);
		lines.push_back(0xffffffff);
	};
#endif */

	triPairConn tpc;
	std::set<trianglesUVW::edge,trianglesUVW::edgeTest> borderEdges;
	if(nPunch<1) { // simple H
		assert(nPoints>1);
		unsigned long Tstart[4],Tend[4];
		if(!Tin(posV[0],nrmV[0],posV[1],nrmV[1],0,Tstart))
			return false;
		_backPunchVerts[0]=Tstart[1];
		_backPunchVerts[1]=Tstart[0];
		_backPunchVerts[2]=Tstart[3];
		_backPunchVerts[3]=Tstart[2];
		if(!Tin(posV[1],nrmV[1],posV[0],nrmV[0],2,Tend))
			return false;
		if(!cutTborders(Tstart,Tend,true,borderEdges))
			return false;
		deleteInteriorTriangles(getDeleteTriangle(Tstart[0],Tstart[1]),borderEdges);
		// will never do doubleDelete if fillHoles() called previously
		_tri->makeVertexToTriangleMap();
		cleanCutEdge(borderEdges,_prefEdgeLength*0.2f);
		// no triangle pairs necessary as both ends were filled
		fillHoles();
		_tri->cleanAndPack();
		_tri->_adjacenciesComputed=false;
		_tri->findAdjacentTriangles();
		return true;
	}
	else if(nPunch<2) { // double H, just Tin or just Tout
		std::vector<int> topVerts,bottomVerts;
		if(!getCookieCutter2(posV,nrmV,TstartTriangleEdge!=0x03,Tout,topVerts,bottomVerts) )
			return false;
		_frontPunchVerts[0]=topVerts[0];
		_frontPunchVerts[1]=topVerts[1];
		_frontPunchVerts[2]=bottomVerts[0];
		_frontPunchVerts[3]=bottomVerts[1];
		_backPunchVerts[0]=topVerts[1];
		_backPunchVerts[1]=topVerts[0];
		_backPunchVerts[2]=bottomVerts[1];
		_backPunchVerts[3]=bottomVerts[0];
		// midpoint loaded on both ends
		if(nPoints<3) { // either a Tin or a Tout, not both
			std::vector<unsigned int> path;
			Vec3f n0,n1;
			_tri->getMeanVertexNormal(topVerts[0],n0._v);
			_tri->getMeanVertexNormal(topVerts[1],n1._v);
			n0 += n1; // each already normalized
			path.push_back(topVerts[1]);

//			addSegment(topVerts[1],topVerts[0]);

			if(!segmentPath(topVerts[1],topVerts[0],n0,path))
				return false;
			if(!addBorderEdges(path,borderEdges))
				return false;
			_tri->getMeanVertexNormal(bottomVerts[0],n0._v);
			_tri->getMeanVertexNormal(bottomVerts[1],n1._v);
			n0 += n1; // each already normalized
			path .clear();
			path.push_back(bottomVerts[1]);

//			addSegment(bottomVerts[1],bottomVerts[0]);

			if(!segmentPath(bottomVerts[1],bottomVerts[0],n0,path))
				return false;
			if(!addBorderEdges(path,borderEdges))
				return false;

		}
	}
	else
		;
	// punch with Tin, Tout or both
	if(TstartTriangleEdge!=0x03){ // Tin
		unsigned long Tstart[4];
		if(nPunch>0) {
			if(!Tin(posV[0],nrmV[0],posV[1],nrmV[1],1,Tstart))
				return false; }
		else {
			if(!Tin(posV[0],nrmV[0],posV[1],nrmV[1],0,Tstart))
				return false; }
		if(!cutTborders(Tstart,_frontPunchVerts,false,borderEdges))
			return false;
		// don't do deletes in a double H
		if(nPunch!=1 || nPoints!=3) {
			deleteInteriorTriangles(getDeleteTriangle(Tstart[0],Tstart[1]),borderEdges);
			// never do doubleDeletes in a pure Tin if fillHoles() called previously
			_tri->makeVertexToTriangleMap();
			cleanCutEdge(borderEdges,_prefEdgeLength*0.2f);
			if(nPoints==2) {
				Vec3f v0,v1;
				_tri->getVertexCoordinate(_frontPunchVerts[0],v0._v);
				_tri->getVertexCoordinate(Tstart[1],v1._v);
				v1 -= v0;
				tpc.topLen = v1.length();
				_tri->getVertexCoordinate(_frontPunchVerts[2],v0._v);
				_tri->getVertexCoordinate(Tstart[3],v1._v);
				v1 -= v0;
				tpc.botLen = v1.length();
				tpc.topV = _frontPunchVerts[0];
				tpc.botV = _frontPunchVerts[2];
				tpc.topCounterclockwise = true;
				_triPairs.push_back(tpc);

				_tri->getVertexCoordinate(_frontPunchVerts[1],v0._v);
				_tri->getVertexCoordinate(Tstart[0],v1._v);
				v1 -= v0;
				tpc.topLen = v1.length();
				_tri->getVertexCoordinate(_frontPunchVerts[3],v0._v);
				_tri->getVertexCoordinate(Tstart[2],v1._v);
				v1 -= v0;
				tpc.botLen = v1.length();
				tpc.topV = _frontPunchVerts[1];
				tpc.botV = _frontPunchVerts[3];
				tpc.topCounterclockwise = false;
				_triPairs.push_back(tpc);
			}
		}
	}
	if(nPunch!=1)
		borderEdges.clear();
	if(Tout){ // Tout
		unsigned long Tstart[4];
		int n = (int)posV.size();
		if(nPunch>0) {
			if(!Tin(posV[n-1],nrmV[n-1],posV[n-2],nrmV[n-2],2,Tstart))
				return false; }
		else {
			if(!Tin(posV[n-1],nrmV[n-1],posV[n-2],nrmV[n-2],0,Tstart))
				return false; }
		if(!cutTborders(Tstart,_backPunchVerts,false,borderEdges))
			return false;
		long killTri1,killTri0=getDeleteTriangle(Tstart[0],Tstart[1]);
		deleteInteriorTriangles(killTri0,borderEdges);
		killTri1=getDeleteTriangle(Tstart[2],Tstart[3]);
		if(killTri1<0) {
			_tri->makeVertexToTriangleMap();
			cleanCutEdge(borderEdges,_prefEdgeLength*0.2f); }
		// could be doing Tout to an unfilled punch
		if(killTri1>-1) { // not already deleted so not a filled punch
			tpc.topV = Tstart[0]; // problem is on Tout to beginning of closed loop, these bump into beginning triPairs
			tpc.botV = Tstart[2];
			Vec3f v0,v1;
			_tri->getVertexCoordinate(_backPunchVerts[1],v0._v);
			_tri->getVertexCoordinate(Tstart[0],v1._v);
			v1 -= v0;
			tpc.topLen = v1.length();
			_tri->getVertexCoordinate(_backPunchVerts[3],v0._v);
			_tri->getVertexCoordinate(Tstart[2],v1._v);
			v1 -= v0;
			tpc.botLen = v1.length();
			tpc.topCounterclockwise = true;
			_triPairs.push_back(tpc);

			tpc.topV = Tstart[1];
			tpc.botV = Tstart[3];
			tpc.topCounterclockwise = false;
			_tri->getVertexCoordinate(_backPunchVerts[0],v0._v);
			_tri->getVertexCoordinate(Tstart[1],v1._v);
			v1 -= v0;
			tpc.topLen = v1.length();
			_tri->getVertexCoordinate(_backPunchVerts[2],v0._v);
			_tri->getVertexCoordinate(Tstart[3],v1._v);
			v1 -= v0;
			tpc.botLen = v1.length();
			_triPairs.push_back(tpc);
			deleteInteriorTriangles(killTri1,borderEdges);
			_tri->makeVertexToTriangleMap();
			cleanCutEdge(borderEdges,_prefEdgeLength*0.2f);
		}
		if((nPunch==1 && nPoints==3) || nPoints==2) { // double H or 2 point Tout
			tpc.topV = _backPunchVerts[0];
			tpc.botV = _backPunchVerts[2];
			tpc.topCounterclockwise = true;
			Vec3f v0,v1;
			_tri->getVertexCoordinate(_backPunchVerts[0],v0._v);
			_tri->getVertexCoordinate(Tstart[1],v1._v);
			v1 -= v0;
			tpc.topLen = v1.length();
			_tri->getVertexCoordinate(_backPunchVerts[2],v0._v);
			_tri->getVertexCoordinate(Tstart[3],v1._v);
			v1 -= v0;
			tpc.botLen = v1.length();
			_triPairs.push_back(tpc);

			tpc.topV = _backPunchVerts[1];
			tpc.botV = _backPunchVerts[3];
			tpc.topCounterclockwise = false;
			_tri->getVertexCoordinate(_backPunchVerts[1],v0._v);
			_tri->getVertexCoordinate(Tstart[0],v1._v);
			v1 -= v0;
			tpc.topLen = v1.length();
			_tri->getVertexCoordinate(_backPunchVerts[3],v0._v);
			_tri->getVertexCoordinate(Tstart[2],v1._v);
			v1 -= v0;
			tpc.botLen = v1.length();
			_triPairs.push_back(tpc);
		}
	}
	_tri->makeVertexToTriangleMap();
	addCornerTrianglePairs();
	fillHoles();
	_tri->cleanAndPack();
	_tri->_adjacenciesComputed=false;
	_tri->findAdjacentTriangles();

//	_wgg->getLines()->addLines(points,lines);

	return true;
}

void incision::deleteInteriorTriangles(long triangle, std::set<trianglesUVW::edge,trianglesUVW::edgeTest> &borderEdges)
{ // only fixes _vertexFace[] of touched border vertices, not stranded vertices in middle
	if(*_tri->triVerts(triangle)<0)	// signals a deleted triangle
			return;
	unsigned long at[3];
	at[2] = _tri->_tris[triangle*3+2];
	for(int i=0; i<3; ++i)	{
		trianglesUVW::edge edg;
		edg.vtxMin = _tri->_tris[triangle*3+i];
		if(at[2]<edg.vtxMin) {
			edg.vtxMax = edg.vtxMin;
			edg.vtxMin = at[2];
			at[2] = edg.vtxMax;}
		else {
			edg.vtxMax = at[2];
			at[2] = edg.vtxMin; }
		at[i] = _tri->_adjs[triangle*3+((i+2)%3)];
		_tri->_adjs[triangle*3+((i+2)%3)] = 0x03;
		if(at[i]==0x00000003)
			at[i] = 0xffffffff;
		if(borderEdges.find(edg)!=borderEdges.end()) {
			if(at[i]<0xffffffff) {
				_tri->_vertexFace[_tri->_tris[(at[i]>>2)*3+(at[i]&0x03)]] = (at[i]>>2)|0x40000000;
				_tri->_adjs[(at[i]>>2)*3+(at[i]&0x03)] = 0x03;
			}
			at[i] = 0xffffffff;
		}
	}
	_tri->_tris[triangle*3] = -1;
	for(int i=0; i<3; ++i)
		if(at[i]<0xffffffff) {
			_tri->_adjs[(at[i]>>2)*3+(at[i]&0x03)] = 0x03;
			deleteInteriorTriangles(at[i]>>2,borderEdges);
		}
}

bool incision::getCookieCutter2(const std::vector<Vec3f> &positions, const std::vector<Vec3f> &normals, bool startT, bool endT, std::vector<int> &topVerts,std::vector<int> &bottomVerts)
{	// Inputs user positions and normals of n size. Outputs 2*nPoints newPositions having same normals specifying cookie cutter to make incision.
	int nOut,nPnts=(int)positions.size();
	nOut=nPnts;
	if(startT) --nOut;
	if(endT) --nOut;
	topVerts.clear();
	bottomVerts.clear();
	topVerts.reserve(nOut*2);
	bottomVerts.reserve(nOut*2);
	Vec3f q,plane,planeLast,vP,vM,norm;
	float pD1,pD2,pLD1,pLD2,nD;
	for(int i=0; i<nPnts; ++i)	{
		if(i<nPnts-1)	{
			norm = normals[i] + normals[i+1];
			q = positions[i+1] - positions[i];
		}
		else {
			norm = normals[i] + normals[i-1];
			q = positions[i] - positions[i-1];
		}
		plane = norm^q;
		plane.normalize();
		pD1=pD2=plane*positions[i];
		pD1 -= _incisionWidthOver2;
		pD2 += _incisionWidthOver2;
		if(i<1 || i>nPnts-2)	{
			vM = positions[i] - plane*_incisionWidthOver2;
			vP = positions[i] + plane*_incisionWidthOver2;
			norm = -norm;
			norm.normalize();
		}
		else	{
			if(planeLast*plane<0.0f)
				norm = planeLast-plane;
			else
				norm = planeLast+plane;
			q = positions[i+1] - positions[i-1];
			norm = norm^q;
			if(fabs(norm._v[0])<1.0e-5f && fabs(norm._v[1])<1.0e-5f && fabs(norm._v[2])<1.0e-5f)
				norm = -normals[i];
			else if(norm*normals[i]>=0.0f)
				norm = -norm;
			else
				;
			norm.normalize();
			nD = positions[i]*normals[i];
			Mat3x3f M(planeLast._v[0],plane._v[0],normals[i][0],planeLast._v[1],plane._v[1],normals[i][1],planeLast._v[2],plane._v[2],normals[i][2]);
			q.set(pLD1,pD1,nD);
			vM = M.Robust_Solve_Linear_System(q);
			q.set(pLD2,pD2,nD);
			vP = M.Robust_Solve_Linear_System(q);
		}

norm = -normals[i];

		if((i<1 && startT) || (i>nPnts-2 && endT))
				;
		else {
			int topV,bottomV;
			if(!hitTopBottom(vM,norm,topV,bottomV))
				return false;
			topVerts.push_back(topV);
			bottomVerts.push_back(bottomV);
			if(!hitTopBottom(vP,norm,topV,bottomV))
				return false;
			topVerts.push_back(topV);
			bottomVerts.push_back(bottomV);
		}
		planeLast = plane;
		pLD1=pD1;	pLD2=pD2;
	}
	return true;
}

bool incision::hitTopBottom(Vec3f &pt, Vec3f &nrm, int &topV, int &bottomV)
{
	int ttr,btr;
	float tuv[2],buv[2];
	if(!getTopBottomHits(pt,nrm,ttr,tuv,btr,buv))
		return false;
	auto getMakeV = [this](int tr, float (&uv)[2]) ->int{
	if(uv[0]<0.01f) {
		if(uv[1]<0.01f)
			return _tri->triVerts(tr)[0];
		else if(uv[1]>0.99f)
			return _tri->triVerts(tr)[2];
		else
			return _tri->splitTriangleEdge(tr,2,1.0f-uv[1]);
	}
	else if(uv[1]<0.01f) {
		if(uv[0]>0.99f)
			return _tri->triVerts(tr)[1];
		else
			return _tri->splitTriangleEdge(tr,0,uv[0]);
	}
	else if(uv[0]+uv[1]>0.99f)
		return _tri->splitTriangleEdge(tr,1,uv[1]);
	else
		return _tri->addNewVertexInMidTriangle(tr,uv);
	};
	topV = getMakeV(ttr,tuv);
	bottomV = getMakeV(btr,buv);
	return true;
}

bool incision::getTopBottomHits(const Vec3f &pt, const Vec3f &nrm, int &topTri, float (&topUV)[2], int &botTri, float (&botUV)[2])
{
	std::vector<float> lpos,params;
	std::vector<int> triangles;
	float nD = 1.0e30f;
	int k=30000,i=0,j=_tri->linePick(pt._v,nrm._v,lpos,triangles,params);
	if(j<1)
		return false;
	assert(j>1);
	float nextW;
	while(i<j)	{
		nextW = _tri->_uvw[_tri->_tris[triangles[i]*3]*3+2];
		if(fabs(params[i])<nD && nextW>0.99f) {
			nextW = _tri->_uvw[_tri->_tris[triangles[i+1]*3]*3+2];
			if(i+1<j && nextW<0.01f) {
				nD = fabs(params[i]);
				k=i; }
		}
		++i;
	}
	if(k>j-2)
		return false;
	float *pv=&lpos[k*3];
	_tri->getBarycentricProjection(triangles[k],pv[0],pv[1],pv[2],topUV);
	if(topUV[0]>1.001f || topUV[1]>1.001f || topUV[0]<-0.001f || topUV[1]<-0.001f || topUV[0]+topUV[1]>1.001f)
		return false;
	topTri = triangles[k];
	++k;
	pv=&lpos[k*3];
	_tri->getBarycentricProjection(triangles[k],pv[0],pv[1],pv[2],botUV);
	if(botUV[0]>1.001f || botUV[1]>1.001f || botUV[0]<-0.001f || botUV[1]<-0.001f || botUV[0]+botUV[1]>1.001f)
		return false;
	botTri = triangles[k];
	return true;
}

bool incision::getTriangleIntersect(const float (&position)[3],const float (&normal)[3], const float targetW, int &triangle,float (&UV)[2])
{
	boundingBox<float> rayBox,triBox;
	rayBox.Empty_Box();
	float v[3],rayParam,triUV[2];
	for(int i=0; i<3; ++i)	v[i]=position[i]+normal[i]*10.0f; // look deep along normal, but not superficial
	rayBox.Enlarge_To_Include_Point(v);
	for(int i=0; i<3; ++i)	v[i]=position[i]-normal[i]*0.03f;
	rayBox.Enlarge_To_Include_Point(v);
	std::map<float,rayIntersect> rayIntersects;
	int nTris = (int)_tri->_tris.size()/3;
	for(int j,i=0; i<nTris; ++i)	{
		int *trV = _tri->triVerts(i);
		if(*trV<0)
			continue;
		triBox.Empty_Box();
		for(j=0; j<3; ++j)	{
			GLfloat *gvp = _tri->vertexCoordinate(trV[j]);
			v[0]=gvp[0]; v[1]=gvp[1]; v[2]=gvp[2];
			triBox.Enlarge_To_Include_Point(v);
		}
		if(rayBox.Intersection(triBox))	{
			if(rayTestTriangle(_tri,i,position,normal,rayParam,triUV))
				rayIntersects.insert(std::make_pair(rayParam,rayIntersect(_tri,i,triUV)));
		}
	}
	std::map<float,rayIntersect>::iterator rit=rayIntersects.begin();
	while(rit!=rayIntersects.end()) {
		float junkW = _tri->_uvw[(_tri->_tris[rit->second.triNum*3])*3+2];
		if(targetW==_tri->_uvw[(_tri->_tris[rit->second.triNum*3])*3+2]) // take first one from the top
			break;
		++rit;
	}
	if(rit==rayIntersects.end())
		return false;
	triangle= rit->second.triNum;
	UV[0]=rit->second.triUV[0]; UV[1]=rit->second.triUV[1];
	return true;
}

bool incision::rayTestTriangle(trianglesUVW* tr, int triangleNumber, const float (&position)[3], const float (&normal)[3], float &rayParameter, float (&triUV)[2])
{	// assumes normals have been normalized.
	int *tv=tr->triVerts(triangleNumber);
	const float *V0=tr->vertexCoordinate(tv[0]),
		*V1=tr->vertexCoordinate(tv[1]),
		*V2=tr->vertexCoordinate(tv[2]);
	Vec3f P(V1[0]-V0[0],V1[1]-V0[1],V1[2]-V0[2]),Q(V2[0]-V0[0],V2[1]-V0[1],V2[2]-V0[2]);
	Mat3x3f m((const Vec3f&)normal,P,Q);
	P.set(position[0]-V0[0],position[1]-V0[1],position[2]-V0[2]);
	Vec3f res = m.Robust_Solve_Linear_System(P);
	if(res[0]<-1e-5f || res[1]<-1e-5f || res[2]<-1e-5f || res[1]+res[2]>1.0001f)
		return false;
	rayParameter = res[0];	// res[1]=u, res[2]=v of P(u,v) of input triangle
	triUV[0] = res[1];
	triUV[1] = res[2];
	return true;
}

bool incision::getTriEdgePath(trianglesUVW* tr, int vFrom, int vTo, Vec3f &plane, float d, std::vector<unsigned long> &triEdgePath, std::vector<float> &params)
{
	triEdgePath.clear();
	params.clear();
	const float *now;
	Vec3f P,T,N,dir,from;
	now = tr->vertexCoordinate(vFrom);
	from.set(now[0],now[1],now[2]);
	now = tr->vertexCoordinate(vTo);
	P.set(now[0],now[1],now[2]);
	dir = P-from;
	std::vector<trianglesUVW::neighborNode> nei;
	std::vector<trianglesUVW::neighborNode>::iterator nit,max_nit;
	tr->getNeighbors(vFrom,nei);
	float a,s,t,max_s,max_t=-1.0f,dotMax=-1.0f,dirSq=dir*dir;
	assert(dirSq>1e-16f);
	nit=nei.begin();
	if(nit->triangle<0)	{
		now = tr->vertexCoordinate(nei.front().vertex);
		++nit;	}
	else
		now = tr->vertexCoordinate(nei.back().vertex);
	P.set(now[0],now[1],now[2]);
    max_nit = nei.end();
	for(; nit!=nei.end(); ++nit)	{
		if(nit->vertex==vTo) {
			return true;  } // the points must be connected with no triEdgePath
		now = tr->vertexCoordinate(nit->vertex);
		T.set(now[0],now[1],now[2]);
		T -= P;
		a = plane*T;
		if(fabs(a)<1e-16)	{
			P.set(now[0],now[1],now[2]);
			continue;	}
		s = (d - plane*P)/a;
		if(s<-1e-16f || s>1.0001f)	{
			P.set(now[0],now[1],now[2]);
			continue;	}
		t = ((P + T*s - from)*dir)/dirSq;
		if(t<0.0f || t>1.0001f)	{
			P.set(now[0],now[1],now[2]);
			continue;	}
		N = P + T*s;
		N -= from;
		N.normalize();
		a = N*dir;
		if(a>dotMax)	{
			dotMax = a;
			max_t=t;
			max_s=s;
			max_nit=nit;
		}
		P.set(now[0],now[1],now[2]);
	}
	assert(max_nit!=nei.end());
	int i,*tri=tr->triVerts(max_nit->triangle);
	for(i=0; i<3; ++i)
		if(tri[i]==vFrom)	break;
	assert(i<3);
	now = tr->vertexCoordinate(tri[(i+1)%3]);
	T.set(now[0],now[1],now[2]);
	s = plane*T-d;
	now = tr->vertexCoordinate(tri[(i+2)%3]);
	T.set(now[0],now[1],now[2]);
	t = plane*T-d;
	assert(s<0.0f != t<0.0f);
	params.push_back(1.0f-max_s);
	unsigned long nextTE = tr->_adjs[max_nit->triangle*3 + ((i+1)%3)];
	triEdgePath.push_back(nextTE);
	while(true) {
		if(nextTE==triEdgePath.front() && triEdgePath.size()>1)
			return false; // no path to vTo
		tri=tr->triVerts(nextTE>>2);
		i = tri[((nextTE&0x03)+2)%3];
		if(i==vTo)
			break;
		now = tr->vertexCoordinate(i);
		T.set(now[0],now[1],now[2]);
		a = plane*T-d;
		if(s<0.0f != a<0.0f) {
			t = a;
			i = ((nextTE&0x03)+1)%3; }
		else {
			s = a;
			i = ((nextTE&0x03)+2)%3; }
		nextTE = tr->_adjs[(nextTE>>2)*3 + i];
		triEdgePath.push_back(nextTE);
		if(fabs(t-s)<1e-16f)
			params.push_back(0.5f);
		else
			params.push_back(t/(t-s));
	}
	return true;
}

bool incision::segmentPath(int vertexFrom, int vertexTo, const Vec3f &surfaceNormal, std::vector<unsigned int> &vertexPath)
{
	Vec3f from,dir,target;
	double d,plane[3];
	const float *vp=_tri->vertexCoordinate(vertexFrom);
	from.set(vp[0],vp[1],vp[2]);
	vp =_tri->vertexCoordinate(vertexTo);
	target.set(vp[0],vp[1],vp[2]);
	dir = from - target;
	from = dir^surfaceNormal;
	from.normalize();
	plane[0]=from._v[0]; plane[1]=from._v[1]; plane[2]=from._v[2];
	d = plane[0]*target[0] + plane[1]*target[1] + plane[2]*target[2];
	float param=0.0,dirParam=1.0f,lastDP=1.0f;
	std::vector<unsigned long> triEdgePath;
	std::vector<float> params;
	unsigned long TE,TEorig,unused;
	TE = _tri->_vertexFace[vertexFrom]&0x3fffffff;
	int j,*tri=&_tri->_tris[TE*3];
	for(j=0; j<3; ++j)
		if(tri[j]==vertexFrom)
			break;
	assert(j<3);
	TE <<= 2;
	TE += j;
	TEorig = TE;
	bool pathDone=true,restart=false;
	if(getFirstStep(false,plane,d,target,dir,dirParam,vertexTo,TE,param)==1) {
		lastDP = dirParam;
		triEdgePath.push_back(TE);
		params.push_back(param);
		pathDone=false; }
	while(!pathDone) {
		while(topoStep(plane,d,target,dir,dirParam,vertexTo,TE,param,unused)) {
			if(dirParam>lastDP || dirParam<0.0f) { // possible wrong path picked. Try the other one.
				if(restart) // already tried second path
					return false;
				param=0.0; dirParam=1.0f; lastDP=1.0f;
				triEdgePath.clear();
				params.clear();
				int secondStep;
				secondStep=getFirstStep(true,plane,d,target,dir,dirParam,vertexTo,TE,param);
				if(secondStep==0) // unlikely
					break;
				else if(secondStep==2) // there is no acceptable second path
					return false;
				else ;
				restart=true;
			}
			lastDP = dirParam;
			triEdgePath.push_back(TE);
			params.push_back(param);
		}
		pathDone=true;
	}
	for(int n=(int)triEdgePath.size(),i=0; i<n; ++i) {
		TE = triEdgePath[i];
		if(params[i]==0.0f) {
			j = _tri->_tris[(TE>>2)*3+(TE&0x03)];
			if(!vertexPath.empty() && j==vertexPath.back())
				continue;
			vertexPath.push_back(j);
		}
		else {
			vertexPath.push_back(_tri->splitTriangleEdge(TE>>2,TE&0x03,params[i]));
			param = _tri->_uvw[vertexPath.back()*3+2];
		}
	}

//	topoCheck();  // GOOD - addNewVertexInMidTriangle() leaves topological arrays adjusted appropriately

/*	std::vector<unsigned long> triEdgePath2;
	std::vector<float> params2;
	if(!getTriEdgePath2(_tri,vertexFrom,vertexTo,surfaceNormal,triEdgePath2,params2))
		return false;
	for(int n=(int)triEdgePath.size(),i=0; i<n; ++i) {
		unsigned long te = triEdgePath[i];
		vertexPath.push_back(tr->splitTriangleEdge(te>>2,te&0x03,params[i]));
	} */

	vertexPath.push_back(vertexTo);
	return true;
}


bool incision::getTriEdgePath2(trianglesUVW* tr, int vFrom, int vTo, const Vec3f &surfaceNormal, std::vector<unsigned long> &triEdgePath, std::vector<float> &params)
{
	triEdgePath.clear();
	params.clear();
	const float *now;
	Vec3f P,T,N,dir,from,target,plane;
	now = tr->vertexCoordinate(vFrom);
	from.set(now[0],now[1],now[2]);
	now = tr->vertexCoordinate(vTo);
	target.set(now[0],now[1],now[2]);
	dir = target-from;
	float a,s,t,d,max_s,max_t=-1.0f,dotMax=-1.0f,dirSq=dir*dir;
	plane = dir^surfaceNormal;
	d = plane*target;
	std::vector<trianglesUVW::neighborNode> nei;
	std::vector<trianglesUVW::neighborNode>::iterator nit,max_nit;
	tr->getNeighbors(vFrom,nei);
	if(dirSq<1e-16f)
		return false;
	nit=nei.begin();
	if(nit->triangle<0)	{
		now = tr->vertexCoordinate(nei.front().vertex);
		++nit;	}
	else
		now = tr->vertexCoordinate(nei.back().vertex);
	P.set(now[0],now[1],now[2]);
    max_nit = nei.end();
	for(; nit!=nei.end(); ++nit)	{
		if(nit->vertex==vTo) {
			return true;  } // the points must be connected with no triEdgePath
		now = tr->vertexCoordinate(nit->vertex);
		T.set(now[0],now[1],now[2]);
		T -= P;
		a = plane*T;
		if(fabs(a)<1e-16)	{
			P.set(now[0],now[1],now[2]);
			continue;	}
		s = (d - plane*P)/a;
		if(s<-1e-16f || s>1.0001f)	{
			P.set(now[0],now[1],now[2]);
			continue;	}
		t = ((P + T*s - from)*dir)/dirSq;
		if(t<0.0f || t>1.0001f)	{
			P.set(now[0],now[1],now[2]);
			continue;	}
		N = P + T*s;
		N -= from;
		N.normalize();
		a = N*dir;
		if(a>dotMax)	{
			dotMax = a;
			max_t=t;
			max_s=s;
			max_nit=nit;
		}
		P.set(now[0],now[1],now[2]);
	}
	if(max_nit==nei.end())
		return false;
	int i,*tri=tr->triVerts(max_nit->triangle);
	for(i=0; i<3; ++i)
		if(tri[i]==vFrom)	break;
	assert(i<3);
	unsigned long nextTE = tr->_adjs[max_nit->triangle*3+((i+1)%3)];
	triEdgePath.push_back(nextTE);
	params.push_back(1.0f-max_s);
	return getTriEdgePathToTarget(tr,vTo,surfaceNormal,triEdgePath,params);
}

bool incision::getTriEdgePathToTarget(trianglesUVW* tr, int vTo, const Vec3f &surfaceNormal, std::vector<unsigned long> &triEdgePath, std::vector<float> &params)
{ // assumes triEdgePath and params already seeded with a single starting element
	unsigned long nextTE;
	int vNext,triangle,edge=triEdgePath.front();
	triangle = edge>>2;
	edge &= 0x03;
	Vec3f sV,tV,T,plane,target;
	const float *now;
	now = tr->vertexCoordinate(vTo);
	target.set(now[0],now[1],now[2]);
	float d,s,t,a,dist2,paramLast;
	int *tri = tr->triVerts(triangle);
	now = tr->vertexCoordinate(tri[edge]);
	sV.set(now[0],now[1],now[2]);
	now = tr->vertexCoordinate(tri[(edge+1)%3]);
	tV.set(now[0],now[1],now[2]);
	T = sV*(1.0f-params.front()) + tV*params.front();
	T -= target;
	dist2 = T*T;
	plane = surfaceNormal^T;
	d = plane*target;
	s = plane*sV-d;
	t = plane*tV-d;
	if(s<0.0f == t<0.0f)
		return false;
	while(true) {
		vNext = tri[(edge+2)%3];
		if(vNext==vTo)
			break;
		now = tr->vertexCoordinate(vNext);
		T.set(now[0],now[1],now[2]);
		a = plane*T-d;
		if(s<0.0f != a<0.0f) {
			t = a;
			tV = T;
			edge = (edge+2)%3; }
		else {
			s = a;
			sV = T;
			edge = (edge+1)%3; }
		nextTE = tr->_adjs[triangle*3 + edge];

		if(nextTE==7488)
			a=a;

		if(nextTE==triEdgePath.front() && triEdgePath.size()>1)
			return false; // no path to vTo
		triEdgePath.push_back(nextTE);
		if(fabs(t-s)<1e-16f)
			paramLast = 0.5f;
		else
			paramLast = t/(t-s);
		params.push_back(1.0f-paramLast);
		triangle = nextTE>>2;
		if(triangle==1872)
			a=a;
		edge = nextTE&0x03;
		tri=tr->triVerts(triangle);
		T = sV*paramLast + tV*(1.0f-paramLast);
		T -= target;
		a = T*T;
		if(a>dist2)
			a=a;
//			return false;
		dist2 = a;
		plane = surfaceNormal^T; // correct plane if necessary as target is approached
		d = plane*target;
		s = plane*sV-d;
		t = plane*tV-d;
//		if(fabs(t-s)<1e-16f)
//			params.push_back(0.5f);
//		else
//			params.push_back(-s/(t-s));
		if(s<0.0f == t<0.0f)
			s=s;
	}
	return true;
}

bool incision::topoPath2(Vec3f &fromV, int fromTri, Vec3f &toV, int &toTri, Vec3f &normal, int &stopEdge, float &stopParam)
{	// Tries to find a topologically connected path from position fromV on triangle fromTri in direction of position
	// toV on triangle toTri on a line connecting those positions in a plane containing mean surface normal. Returns true
	// if a hard incision edge is encountered returning the top triangle in toTri, position in toV, the stopping triangle
	// edge in stopEdge and the parameter along that edge in stopParam. Returns false if surface bends out of path plane,
	// so no hard edge T possible and returns data from last triangle edge cut.
	Vec3f from,dir,vLast,vNow;
	float txW=_tri->_uvw[_tri->_tris[toTri*3]*3+2];
	double d,plane[3];
	dir = fromV - toV;
	from = dir^normal;
	from.normalize();
	plane[0]=from._v[0]; plane[1]=from._v[1]; plane[2]=from._v[2];
	d = plane[0]*toV[0] + plane[1]*toV[1] + plane[2]*toV[2];
	double dNow,dLast;
	int edge=-1,i,*tr=&_tri->_tris[fromTri*3];
	_tri->getVertexCoordinate(tr[2],vLast._v);
	dLast = plane[0]*vLast[0] + plane[1]*vLast[1] + plane[2]*vLast[2] - d;
	float t,tBest=1e30f,dirSq=dir*dir;
	for(i=0; i<3; ++i) {
		_tri->getVertexCoordinate(tr[i],vNow._v);
		dNow = plane[0]*vNow[0] + plane[1]*vNow[1] + plane[2]*vNow[2] - d;
		if(dNow>0.0 != dLast>0.0) {
			t = (float)(-dLast/(dNow-dLast)); // 0 div not possible
			vLast = vLast*(1.0f-t) + vNow*t;
			t = dir*(vLast-toV)/dirSq;
			if(t<tBest) {
				tBest=t;
				edge = (i+2)%3;
			}
		}
		vLast = vNow;
		dLast = dNow;
	}
	assert(edge>-1);
	unsigned long prevTE,TE=_tri->_adjs[fromTri*3+edge];
	t=tBest;
	prevTE = TE;
	while(topoStep(plane,d,toV,dir,tBest,-1,TE,stopParam,prevTE)) {
		if(TE==3) { // free edge
			stopEdge = prevTE&3;
			toTri = prevTE>>2; // works OK for filled edge
			stopParam = 1.0f-stopParam;
			return true; }
		float nextW=_tri->_uvw[_tri->_tris[(TE>>2)*3+(((TE&3)+2)%3)]*3+2]; // careful TE=3 gives an answer here. Filter by above.
		if(txW!=nextW) { // filled hard edge
			unsigned long adj = _tri->_adjs[(TE>>2)*3+(TE&3)];
			stopEdge = adj&3;
			toTri = adj>>2; // works OK for filled edge
			stopParam = 1.0f-stopParam;
			return true; }
		if(fabs(tBest)<1e-16f && (fabs(stopParam)<1e-16f || fabs(stopParam-1.0f)<1e-16f)) { // exact vertex hit
			int hitV = _tri->_tris[(TE>>2)*3+(((TE&3)+(stopParam<1e-16?0:1))%3)];
			std::vector<trianglesUVW::neighborNode> nei;
			if(!getHardEdgeNeighbors(hitV,nei))
				return false;
			toTri = nei.back().triangle; // first or last doesn't matter
			for(i=0; i<3; ++i)
				if(_tri->_tris[toTri*3+i]==hitV)
					break;
			assert(i<3);
			stopEdge = i;
			stopParam = 0.0f;
			return true;
		}
		if(tBest>t || tBest<-0.001f) 
			return false;
		t=tBest;
	}
	return true;
}

bool incision::topoPath(Vec3f &fromV, int fromTri, Vec3f &toV, int &toTri, Vec3f &normal, int &stopEdge, float &stopParam)
{	// Tries to find a topologically connected path from position fromV on triangle fromTri in direction of position
	// toV on triangle toTri on a line connecting those positions in a plane containing mean surface normal. Returns true
	// if a hard incision edge is encountered returning the top triangle in toTri, position in toV, the stopping triangle
	// edge in stopEdge and the parameter along that edge in stopParam. Returns false if surface bends out of path plane,
	// so no hard edge T possible and returns data from last triangle edge cut.
	Vec3f plane=toV-fromV;
	plane = plane^normal;
	plane.normalize();
	float D,dLast,dNext,distSq,dMin=9.0e30f;
	D = plane*fromV;
	int i,lastEdge=-1,nextEdge=-1,originalFromTri=fromTri;
	if(!_tri->_adjacenciesComputed)
		_tri->findAdjacentTriangles();
	int *tr;
	bool startup=true;
	Vec3f vLast,vNext,planeV;
	while(fromTri!=toTri)	{
		tr = _tri->triVerts(fromTri);
		_tri->getVertexCoordinate(tr[2],vLast._v);
		dLast = (plane*vLast) -D;
		nextEdge = -1;
		for(i=0; i<3; ++i)	{
			_tri->getVertexCoordinate(tr[i],vNext._v);
			dNext = (plane*vNext) -D;
			if(dNext<0.0f != dLast<0.0f)	{	// edge cut
				if(i==lastEdge)
					;	// already done, skip it
				else if(startup)	{
					planeV = ((vNext*dLast) - (vLast*dNext))/(dLast-dNext);
					planeV -= toV;
					distSq = planeV*planeV;
					lastEdge = i;
					if(distSq<dMin)	{
						if(nextEdge>-1)
							lastEdge = nextEdge;
						dMin = distSq;
						nextEdge = i;	}
				}
				else	{
					unsigned long adj = _tri->_adjs[fromTri*3+((i+2)%3)];
					if(adj==0x003)	{	// incision edge hit
						stopEdge=(i+2)%3;
						stopParam = -dLast/(dNext-dLast);
						assert(stopParam>=0 && stopParam<=1.0f);
						toV = vLast*stopParam + vNext*(1.0f-stopParam);
						toTri = fromTri;
						return true;
					}
					lastEdge = ((adj&0x03)+1)%3;
					fromTri = adj>>2;
					break;
				}
			}
			dLast = dNext;
			vLast = vNext;
		}
		if(!startup && i>2)	{	// surface bends out of path plane
			// WARNING - COURT has not written this. Shouldn't happen.
			// Should do topological circle forever running through start point
			assert(false);
			return false;
		}
		if(!startup && originalFromTri==fromTri)	{	// surface bends out of path plane
			// WARNING - COURT has not tested this
			// This is topological circle forever running through start point
			assert(false);
			return false;
		}
		startup=false;
	}
	// uninterrupted path to original hard edge input. Don't change input data.
	return true;
}

int incision::incisionEdgeT_Band(int &triEdge, Vec3f &edgeV, const Vec3f &nextV, const Vec3f &surfaceNormal, int &posVert, float &posDist, int &negVert, float &negDist)
{	// inputs first 4 args where edgeV is a point on triEdge which is on an incision edge. nextV is the next interior point of a T in. surfaceNormal is normal at
	// edgeV or mean norm of edge and next V. Routine tries to find a band _incisionWidth wide around edgeV and puts result in last 4 args where posVert is next
	// outside vertex on incision edge on positive side of plane which is positiveDistance in along edge from desired plane intersection.  Neg arguments are the
	// same for the negative side of the plane.  This annoying return is to allow for possible later edge splits of same triangle edge for both points.  If not
	// enough edge room on positive side, return bit 0 is set and edgeV is adjusted and returned.  Not enough room on negative side sets bit 1 and returns new edgeV.
	// Both sides sets both bits and edgeV returned as average of the two limits. Assumes topoPath() already done.

	assert(false); // this goes in version 2

	posVert=-1; posDist=1e30f; negVert=-1; negDist=1e30f;
//	assert(_tri->_adjs[(triEdge>>2)*3+(triEdge&0x03)]==0x03);	// must be hard edge start
	int *tv,ret=3,posMaxV=-1,negMaxV=-1;	// neither side found yet
	Vec3f posLast,posNext,planeN;
	float dLast,dNext,posDlimit,negDlimit,planeD,posMax=0.0f,negMax=0.0f;
	planeN = _tPlaneNormal;
	planeD = _tPlaneD;
	posDlimit = _incisionWidthOver2;
	negDlimit = -_incisionWidthOver2;
	tv = _tri->triVerts(triEdge>>2);
	int vNext=tv[triEdge&0x0003],vLast;
	_tri->getVertexCoordinate(vNext,posNext._v);
	dNext = planeN*posNext - planeD;
	bool start=true,neiTop=false,limitCorrectionDone=false;
	std::vector<trianglesUVW::neighborNode> nei;
	while(true)	{	// fix me
		if(!getHardEdgeNeighbors(vNext,nei))
			vNext=vNext;
		vLast = vNext;
		posLast = posNext;
		dLast = dNext;
		if(neiTop)
			vNext = nei.back().vertex;
		else
			vNext = nei.front().vertex;
		_tri->getVertexCoordinate(vNext,posNext._v);
		dNext = planeN*posNext - planeD;
		if(start)	{	// if first point was in middle of input edge, this should never happen
			assert(!(dNext>posDlimit && dLast>posDlimit) && !(dNext<=negDlimit && dLast<=negDlimit));
			start = false;
		}
		if(dNext>posDlimit || dNext<=negDlimit)	{	// this direction found
			float dN,len = (posNext-posLast).length();
			dN = dNext-dLast;
			if(fabs(dN)>1e-16f)	{
				if(dNext>posDlimit)
					len *= (dNext - posDlimit)/dN;
				else
					len *= (dNext - negDlimit)/dN;
			}
			else	len=0;
			assert(len>=0.0f);
			if(dNext>posDlimit)	{
 				if((ret&0x01)<1)	{	// second positive limit hit
					negVert = negMaxV;
					negDist = 0.0f;
					// turn plane to reflect negative limit
					_tri->getVertexCoordinate(negMaxV,posNext._v);
					planeN=posNext-(nextV-(planeN*_incisionWidthOver2*1.02f));
					planeN = planeN^surfaceNormal;
					planeN.normalize();
					planeD = planeN*nextV;
					_tPlaneNormal = planeN;
					_tPlaneD = planeD;
					// restart search
					vNext = negMaxV;
					_tri->getVertexCoordinate(vNext,posNext._v);
					dNext = planeN*posNext - planeD;
					int stopV=vNext;
					do {
						vLast = vNext;
						posLast = posNext;
						dLast = dNext;
						getHardEdgeNeighbors(vNext,nei);
						assert(nei.front().triangle<0);	// moving along hard edge
						vNext = nei.back().vertex;
//						vNext = nei.front().vertex;
						_tri->getVertexCoordinate(vNext,posNext._v);
						dNext = planeN*posNext - planeD;
					}while(dNext<0.0f && vNext!=stopV);
					edgeV = (posLast*dNext - posNext*dLast)/(dNext-dLast);	// can't be 0
					triEdge = nei.back().triangle;
//					triEdge = nei[1].triangle;
					tv = _tri->triVerts(triEdge);
					for(stopV=0; stopV<3; ++stopV)
						if(tv[stopV]==vNext)	break;
//						if(tv[stopV]==vLast)	break;
					triEdge = stopV + (triEdge<<2);
					int retLim;	// second positive limit hit
					while(retLim=incisionEdgeT_Band(triEdge,edgeV,nextV,surfaceNormal,posVert,posDist,negVert,negDist))
						if(retLim==2)	return 4;	// shouldn't happen unless tiny hole in scene or change in incision width
					return 1;
				}
				ret &= 0x02;
				if(neiTop)	{	// this is the one of possibly two
					posVert = vNext;
					posDist = len;	}
			}
			else	{
				if((ret&0x02)<1)	{	// second negative limit hit
					posVert = posMaxV;
					posDist = 0.0f;
					// turn plane to reflect positive limit
					_tri->getVertexCoordinate(posMaxV,posNext._v);
					planeN=posNext-(nextV+(planeN*_incisionWidthOver2*1.02f));
					planeN = planeN^surfaceNormal;
					planeN.normalize();
					planeD = planeN*nextV;
					_tPlaneNormal = planeN;
					_tPlaneD = planeD;
					// restart search
					vNext = posMaxV;
					_tri->getVertexCoordinate(vNext,posNext._v);
					dNext = planeN*posNext - planeD;
					int stopV=vNext;
					do {
						vLast = vNext;
						posLast = posNext;
						dLast = dNext;
						getHardEdgeNeighbors(vNext,nei);
						assert(nei.front().triangle<0);	// moving along hard edge
						vNext = nei.front().vertex;
						_tri->getVertexCoordinate(vNext,posNext._v);
						dNext = planeN*posNext - planeD;
					}while(dNext>=0.0f && vNext!=stopV);
					edgeV = (posLast*dNext - posNext*dLast)/(dNext-dLast);	// can't be 0
					triEdge = nei[1].triangle;
					tv = _tri->triVerts(triEdge);
					for(stopV=0; stopV<3; ++stopV)
						if(tv[stopV]==vLast)	break;
					triEdge = stopV + (triEdge<<2);
					int retLim;	// second negative limit hit
					while(retLim=incisionEdgeT_Band(triEdge,edgeV,nextV,surfaceNormal,posVert,posDist,negVert,negDist))
						if(retLim==1)	return 8;	// shouldn't happen unless tiny hole in scene or change in incision width
					return 2;
				}
				ret &= 0x01;
				if(!neiTop)	{	// this is the one of possibly two
					negVert = vNext;
					negDist = len;	}
			}
			if(ret<1)	// done
				break;
			neiTop = !neiTop;
		}
		else	{ // Need to keep pos and neg maximums in case of failure to exceed limit.
			if(dNext>posMax)	{
				posMax=dNext;
				posMaxV = vNext;	}
			if(dNext<negMax)	{
				negMax=dNext;
				negMaxV = vNext;	}
		}
	}
	return ret;
}

bool incision::getOppositeIncisionEdge2(int topTri, int topEdge, Vec3f &topPos, Vec3f &dir, int &bottomTri, int &bottomEdge, float &bottomParam, Vec3f &bottomV)
{	// uses lastPlane data to find opposite incision edge from top.  If stays in plane returns true, else false returning closest answer in bottom.
	int *tr,i,j,n=(int)_tri->_tris.size();
	tr = &_tri->_tris[topTri*3];
	float targetW = _tri->_uvw[tr[topEdge]*3+2];
	if(targetW>0.999f)
		targetW=0.0f;
	else
		targetW=1.0f;
	Vec3f pt=topPos,e0,e1,topEdgeDir;
	topEdgeDir.set((float (&)[3])_tri->_xyz[tr[(topEdge+1)%3]*3]);
	e0.set((float (&)[3])_tri->_xyz[tr[topEdge]*3]);
	topEdgeDir -= e0;
	float w[3],minDsq=1e30f;
	unsigned long *adj;
	for(i=0; i<n; i+=3)	{
		tr = &_tri->_tris[i];
		if(*tr<0)
			continue;
		adj = &_tri->_adjs[i];
		for(j=0; j<3; ++j)
			w[j] = _tri->_uvw[(tr[j]<<1)+tr[j]+2];
		for(j=0; j<3; ++j) {
			if(adj[j]==0x03) { // allow for possible unconnected free edge on top or bottom
				if(w[j]==targetW)
					break;
			}
			if(w[j]!=w[(j+2)%3] && w[j]==w[(j+1)%3] && w[j]==targetW)
					break;
		}
		if(j>2) continue;
		// candidate edge
		e0.set((float (&)[3])_tri->_xyz[tr[j]*3]);
		e1.set((float (&)[3])_tri->_xyz[tr[(j+1)%3]*3]);
		Vec3f v0=pt-e0,v1=pt-e1;
		e1 -= e0;
		if(e1*e1<1e-15f)  continue;
		if(adj[j]==0x03) {
//			if(e1*topEdgeDir>0.0f) // bottom edge points in wrong direction
//				continue;
		}
		else { // a connected side edge
			if(e1*topEdgeDir<0.0f) // side edge points in wrong direction
				continue; }
		float d,t=(v0*e1)/(e1*e1);
		if(t<0.0f ){
			d=v0*v0;
			if(dir*(e0-topPos)<0.0f) // dir*dir always positive
				continue;
			if(minDsq>d) {
				minDsq=d;
				bottomV[0] = e0._v[0];
				bottomV[1] = e0._v[1];
				bottomV[2] = e0._v[2];
				bottomTri = i;
				bottomEdge = j;
				bottomParam = 0.0; }
		}
		else if(t>1.0f){
			e1 += e0;
			if(dir*(e1-topPos)<0.0f)
				continue;
			if(minDsq>(d=v1*v1)) {
				minDsq=d;
				bottomV[0] = e1._v[0];
				bottomV[1] = e1._v[1];
				bottomV[2] = e1._v[2];
				bottomTri = i;
				bottomEdge = j;
				bottomParam = 1.0; }
		}
		else {
			v1 = e0 + e1*t;
			v0 = v1 - pt;
			if(dir*(v1-topPos)<0.0f)
				continue;
			if((d=v0*v0)<minDsq) {
				minDsq=d;
				bottomV[0] = v1._v[0];
				bottomV[1] = v1._v[1];
				bottomV[2] = v1._v[2];
				bottomTri = i;
				bottomEdge = j;
				bottomParam = t;
			}
		}
	}
	if(minDsq>1e29f)
		return false;
	// always want top or bottom triangle reference
	unsigned long adjV = _tri->_adjs[bottomTri+bottomEdge];
	if(adjV!=0x03) {
		bottomTri = adjV>>2;
		bottomEdge = adjV&0x03;
		bottomParam = 1.0f - bottomParam;
	}
	else
		bottomTri /= 3;
	return true;
}

bool incision::getOppositeIncisionEdge(int topTri, int topEdge, Vec3f &topPos, Vec3f &dir, int &bottomTri, int &bottomEdge, Vec3f &bottomV)
{	// uses lastPlane data to find opposite incision edge from top.  If stays in plane returns true, else false returning closest answer in bottom.
	unsigned long lastEdge = _tri->_adjs[topTri*3+topEdge];	// botEdge,
	bottomTri = lastEdge>>2;
	lastEdge &= 0x0003;
	dir.normalize();
	Vec3f trV[3];	//,p,q;
	float botP,Dvals[3],targetW;
	std::set<int> visitedTris;
	int *tr = _tri->triVerts(bottomTri);
	do	{
		if(visitedTris.empty())
			targetW = fabs(1.0f - _tri->_uvw[tr[lastEdge]*3+2]);
		for(int i=0; i<3; ++i)	{
			_tri->getVertexCoordinate(tr[i],trV[i]._v);
			Dvals[i] = (_tPlaneNormal*trV[i]) - _tPlaneD;
		}
		if((Dvals[lastEdge]<0.0f != Dvals[(lastEdge+2)%3]<0.0f) &&	// len[(lastEdge+2)%3]>maxLen &&
			visitedTris.find((_tri->_adjs[bottomTri*3+((lastEdge+2)%3)])>>2)==visitedTris.end())	{
			botP = Dvals[lastEdge]/(Dvals[lastEdge]-Dvals[(lastEdge+2)%3]);
			lastEdge = (lastEdge+2)%3;
		}
		else if((Dvals[(lastEdge+1)%3]<0.0f != Dvals[(lastEdge+2)%3]<0.0f) &&	// len[(lastEdge+2)%3]>maxLen &&
			visitedTris.find((_tri->_adjs[bottomTri*3+((lastEdge+1)%3)])>>2)==visitedTris.end())	{
			botP = Dvals[(lastEdge+2)%3]/(Dvals[(lastEdge+2)%3]-Dvals[(lastEdge+1)%3]);
			lastEdge = (lastEdge+1)%3;
		}
		else	{	// surface bending out of plane.
			int e1=_tri->_adjs[bottomTri*3+((lastEdge+1)%3)],e2=_tri->_adjs[bottomTri*3+((lastEdge+2)%3)];
			tr = _tri->triVerts(e1>>2);
			e1 = tr[((e1&0x03)+2)%3];
			tr = _tri->triVerts(e2>>2);
			e2 = tr[((e2&0x03)+2)%3];
			Vec3f v1,v2;
			_tri->getVertexCoordinate(e1,v1._v);
			_tri->getVertexCoordinate(e2,v2._v);
			float i1 = (v1-topPos)*dir;
			if(i1>1e-16f)
				i1 =((_tPlaneNormal*v1) - _tPlaneD)/i1;
			float i2 = (v2-topPos)*dir;
			if(21>1e-16f)
				i2 =((_tPlaneNormal*v2) - _tPlaneD)/i2;
			if(i1<i2 && visitedTris.find((_tri->_adjs[bottomTri*3+((lastEdge+1)%3)])>>2)!=visitedTris.end())	{
				lastEdge = (lastEdge+1)%3;
			}
			else if(visitedTris.find((_tri->_adjs[bottomTri*3+((lastEdge+2)%3)])>>2)!=visitedTris.end())	{
				lastEdge = (lastEdge+2)%3;
			}
			else
				return false;	// shouldn't happen
			if(Dvals[lastEdge]<Dvals[(lastEdge+1)%3])
				botP = 1.0f;
			else
				botP = 0.0f;
		}
		lastEdge = _tri->_adjs[bottomTri*3+lastEdge];
		visitedTris.insert(bottomTri);
		bottomTri = lastEdge>>2;
		lastEdge &= 0x0003;
		tr = _tri->triVerts(bottomTri);
	}while(_tri->_uvw[tr[lastEdge]*3+2]!=targetW || _tri->_uvw[tr[(lastEdge+1)%3]*3+2]!=targetW);
	bottomEdge = lastEdge;
	float *v0=_tri->vertexCoordinate(tr[(lastEdge+1)%3]),*v1=_tri->vertexCoordinate(tr[lastEdge]);
	for(int i=0; i<3; ++i)
		bottomV._v[i] = v0[i]*(1.0f-botP) + v1[i]*botP;
	return true;
}

bool incision::getHardEdgeNeighbors(unsigned int vertex, std::vector<trianglesUVW::neighborNode> &neighbors)
{ // assumes vertex already on hard edge
	float targetW=_tri->_uvw[vertex*3+2];
	if(targetW>0.9999f) targetW=1.0f;
	if(targetW<0.0001f) targetW=0.0f;
	if(targetW!=0.0f && targetW!=1.0f)
		return false;
	neighbors.clear();
	std::vector<trianglesUVW::neighborNode> nei;
	std::vector<trianglesUVW::neighborNode>::iterator nit,nStop;
	_tri->getNeighbors(vertex,nei);
	nit = nei.begin();
	if(nit->triangle<0) {
		if(_tri->_uvw[nit->vertex*3+2]!=_tri->_uvw[nei.back().vertex*3+2])
			return false;
		neighbors.assign(nei.begin(),nei.end());
		return true;
	}
	while(nit!=nei.end() && fabs(targetW-_tri->_uvw[nit->vertex*3+2])<1e-8f)
		++nit;
	if(nit==nei.end())
		return false;
	do{
		++nit;
		if(nit==nei.end())
			nit = nei.begin();
	}while(fabs(targetW-_tri->_uvw[nit->vertex*3+2])>0.001f);
	nStop = nit;
	nit->triangle = -1;
	do{
		neighbors.push_back(*nit);
		++nit;
		if(nit==nei.end())
			nit = nei.begin();
	}while(nit!=nStop && fabs(targetW-_tri->_uvw[nit->vertex*3+2])<1e-8f);
	if(nit==nStop)
		return false;
	return true;
}

int incision::getFirstStep(bool rejectFirst, const double (&plane)[3], const double D, const Vec3f &target, const Vec3f &dir, float &dirParam, int vertexTo, unsigned long &TE, float &param)
{ // Return 0==found end vertex. Return 1==got first step but not finished. Return 2==rejectFirst but there wasn't a second.
	assert(param==0.0f); // vertex start point
	int triangle=TE>>2,edge=TE&0x03;
	int *tri=&_tri->_tris[triangle*3];
	double s,t;
	Vec3f sV,tV,aV;
	float dSq=dir*dir;
	auto vAssign = [](trianglesUVW *tr, int v, Vec3f &vOut) {const float *vp=tr->vertexCoordinate(v); vOut.set(vp[0],vp[1],vp[2]);};
	auto dotdf = [](const double (&d)[3], const float (&f)[3]) -> double { return d[0]*f[0]+d[1]*f[1]+d[2]*f[2]; };
	int j,v=tri[edge];
	std::vector<trianglesUVW::neighborNode> nei;
	std::vector<trianglesUVW::neighborNode>::iterator nit,minNit;
	_tri->getNeighbors(v,nei);
	nit=nei.begin();
	minNit = nei.end();
	if(nit->triangle<0) {
		j = nit->vertex;
		if(nit->vertex==vertexTo)
			return 0;
		++nit; }
	else
		j = nei.back().vertex;
	vAssign(_tri,j,sV);
	s = dotdf(plane,sV._v)-D;
	float p,prm,minPrm,nowDir,secondBestDir=-1e30f,bestDir=-1e30f;
	while(nit!=nei.end()) {
		if(nit->vertex==vertexTo)
			return false;
		vAssign(_tri,nit->vertex,tV);
		t = dotdf(plane,tV._v)-D;
		if(t<0.0 != s<0.0) { // edge straddles plane, usually 2 of these
			p = float(t/(t-s)); // can't be 0 divisor
			aV = sV*p + tV*(1.0f-p);
			aV -= target;
			prm = (aV*dir)/dSq;
			aV.normalize();
			if(prm<dirParam) {
				if((nowDir=aV*dir)>bestDir) {
					bestDir = nowDir;
					if(!rejectFirst) {
						minPrm = prm;
						minNit = nit;
						param = p; }
				}
				if(rejectFirst && (secondBestDir<-1e29f || nowDir<bestDir)) {
					secondBestDir = nowDir;
					minPrm = prm;
					minNit = nit;
					param = p;
				}
			}
		}
		s=t;
		sV=tV;
		++nit;
	}
	if(rejectFirst && secondBestDir==bestDir)
		return 2;
	assert(minNit!=nei.end());
	triangle = minNit->triangle;
	tri=&_tri->_tris[triangle*3];
	for(j=0; j<3; ++j)
		if(tri[j]==v)
			break;
	assert(j<3);
	edge = (j+1)%3;
	TE = _tri->_adjs[triangle*3+edge];
	dirParam = minPrm;
	return 1;
}

bool incision::topoStep(const double (&plane)[3], const double D, const Vec3f &target, const Vec3f &dir, float &dirParam, int vertexTo, unsigned long &TE, float &param, unsigned long &prevTE)
{
	int triangle=TE>>2,edge=TE&0x03;
	int *tri=&_tri->_tris[triangle*3];
	double s,t;
	Vec3f sV,tV,aV;
	float dSq=dir*dir;
	auto vAssign = [](trianglesUVW *tr, int v, Vec3f &vOut) {const float *vp=tr->vertexCoordinate(v); vOut.set(vp[0],vp[1],vp[2]);};
	auto dotdf = [](const double (&d)[3], const float (&f)[3]) -> double { return d[0]*f[0]+d[1]*f[1]+d[2]*f[2]; };
	// edge start point
	int oldEdge=edge,nextV = tri[(edge+2)%3];
	if(nextV==vertexTo)
		return false;
	vAssign(_tri,tri[edge],sV);
	s = dotdf(plane,sV._v)-D;
	vAssign(_tri,tri[(edge+1)%3],tV);
	t = dotdf(plane,tV._v)-D;
	double a;
	vAssign(_tri,nextV,aV);
	a = dotdf(plane,aV._v)-D;
	if(s<0.0 == a<0.0) {
		s = a;
		sV = aV;
		edge = (edge+1)%3;	}
	else {
		t = a;
		tV = aV;
		edge = (edge+2)%3;	}
	assert(fabs(t-s)>1e-18); // no straddling of this edge would have occurred
	param = float(-s/(t-s)); // parameter for next edge, not this one
	aV = sV*(1.0f-param) + tV*param;
	TE = _tri->_adjs[triangle*3+edge];
	prevTE = (triangle<<2)+edge;
	dirParam = ((aV-target)*dir)/dSq;
	return true;
}

void incision::vertexNeighbors(int centerV, int vertex, float rSq, std::set<int> &nei)
{ // diagnostic routine. Not for production use.
	Vec3f center,vNow;
	_tri->getVertexCoordinate(centerV,center._v);
	_tri->getVertexCoordinate(vertex,vNow._v);
	vNow -= center;
	if(vNow*vNow>rSq || nei.find(vertex)!=nei.end())
		return;
	nei.insert(vertex);
	std::vector<trianglesUVW::neighborNode> ne;
	std::vector<trianglesUVW::neighborNode>::iterator nit;
	_tri->getNeighbors(vertex,ne);
	for(nit=ne.begin(); nit!=ne.end(); ++nit)
		vertexNeighbors(centerV,nit->vertex,rSq,nei);
}

bool incision::topoCheck()
{
    std::vector<unsigned long> adjs;	// low 2 bits are the edge number of the adjacent triangle.
	adjs.assign(_tri->_adjs.begin(),_tri->_adjs.end());
	std::vector<unsigned long> vertexFace;
	vertexFace.assign(_tri->_vertexFace.begin(),_tri->_vertexFace.end());
	_tri->_adjacenciesComputed=false; // force computation
	_tri->findAdjacentTriangles();
	int i,n;
	n = (int)adjs.size();
	for(i=0; i<n; ++i) {
		if((i%3==0) && _tri->_tris[i]<0) {
			i+=2;
			continue; }
		if(adjs[i]!=_tri->_adjs[i])
			return false;
	}
	std::vector<trianglesUVW::neighborNode> nei;
	std::vector<trianglesUVW::neighborNode>::iterator nit;
	n = (int)vertexFace.size();
	for(i=0; i<n; ++i) {
		unsigned long vo,vf = vertexFace[i];
		vo = _tri->_vertexFace[i];
		if((vo&0x40000000)>0) {
			if(vo!=vf)
				return false;
		}
		else {
			_tri->getNeighbors(i,nei);
			for(nit=nei.begin(); nit!=nei.end(); ++nit) {
				if(nit->triangle==vf)
					break;
			}
			if(nit==nei.end()) { // could be stranded vertex
				int j;
				for(j=0; j<n; ++j) {
					if((j%3==0) && _tri->_tris[j]<0) {
						j+=2;
						continue; }
					if(_tri->_tris[j]==i)
						break;
				}
				if(j<n) // not a stranded vertex
					return false;
			}
		}
	}
	return true;
}
