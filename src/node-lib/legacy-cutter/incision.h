#ifndef __INCISION_H__
#define __INCISION_H__

#include <vector>
#include <set>
#include <list>
#include "GraphicsUtils/trianglesUVW.h"
#include "GraphicsUtils/Vec3f.h"

// forward declarations
class wxGraphics;

class incision
{
public:
	bool makeIncision(trianglesUVW* tri,int nPoints,float (*positions)[3],float (*normals)[3],int TstartTriangleEdge,float TstartParam,bool Tout);
	void setIncisionWidth(float width)	{_incisionWidthOver2=width/2.0f; }
	void setPreferredEdgeLength(float prefLen) {_prefEdgeLength=prefLen; } // sets preferred length of new triangles created
	void setWxGraphics(wxGraphics *wgg) {_wgg=wgg;}
	incision();
	~incision();

private:
	wxGraphics *_wgg;
	trianglesUVW* _tri;
	bool punch(std::vector<Vec3f> &positions, std::vector<Vec3f> &normals, bool startT, bool endT); // center punch not connected to or intersecting open space. Leaves holes.
	bool getCookieCutter2(const std::vector<Vec3f> &positions, const std::vector<Vec3f> &normals, bool startT, bool endT, std::vector<int> &topVerts,std::vector<int> &bottomVerts);
	bool hitTopBottom(Vec3f &pt, Vec3f &nrm, int &topV, int &bottomV);
	bool getTopBottomHits(const Vec3f &pt, const Vec3f &nrm, int &topTri, float (&topUV)[2], int &botTri, float (&botUV)[2]);
	unsigned long _frontPunchVerts[4],_backPunchVerts[4];
	void deleteInteriorTriangles(long triangleToGo, std::vector<unsigned int> path);
	void trianglePairConnect(const unsigned int *topVerts, const unsigned int *botVerts);
	bool getOppositeIncisionEdge2(int topTri, int topEdge, Vec3f &topPos, Vec3f &dir, int &bottomTri, int &bottomEdge, float &bottomParam, Vec3f &bottomV);
	bool cutTborders(unsigned long (&p)[4], unsigned long (&q)[4], bool downQ, std::set<trianglesUVW::edge,trianglesUVW::edgeTest> &borderEdges);
	bool addBorderEdges(std::vector<unsigned int> &path, std::set<trianglesUVW::edge,trianglesUVW::edgeTest> &borderEdges);
	bool Tin(const Vec3f &startPos, const Vec3f &startNorm, const Vec3f &nextPos, const Vec3f &nextNorm, int punchStatus, unsigned long (&Tvertices)[4]);
	void fillHoles();
	void subdivideLongTriangles(float maxSize, int startTriangle);
	float _prefEdgeLength;
	bool edgeBound(int &vStart, const Vec3f &planeN, bool neiFront, float limitD, float &targetD, unsigned long &triEdge, float &param);
	void cleanCutEdge(std::set<trianglesUVW::edge,trianglesUVW::edgeTest> &borderEdges, float minSize);
	void addCornerTrianglePairs();
	bool topoPath2(Vec3f &fromV, int fromTri, Vec3f &toV, int &toTri, Vec3f &normal, int &stopEdge, float &stopParam);
	void midEdgeLocus(const int v0, const int v1, int &midTri, Vec3f &midPoint);
	std::set<unsigned int> _corners;
	struct vert_triEdge{
		int vertex;
		unsigned long triEdge;
	};
	struct triPairConn{
		int topV;
		int topTri;
		int botV;
		int botTri;
		float topLen;
		float botLen;
		bool topCounterclockwise;
	};
	std::vector<triPairConn> _triPairs;


	bool getHardEdgeNeighbors(unsigned int vertex, std::vector<trianglesUVW::neighborNode> &neighbors);
	bool getOppositeIncisionEdge(int topTri, int topEdge, Vec3f &topPos, Vec3f &dir, int &bottomTri, int &bottomEdge, Vec3f &bottomV);
	bool topoPath(Vec3f &fromV, int fromTri, Vec3f &toV, int &toTri, Vec3f &normal, int &stopEdge, float &stopParam);
	int incisionEdgeT_Band(int &triEdge, Vec3f &edgeV, const Vec3f &nextV, const Vec3f &surfaceNormal, int &posVert, float &posDist, int &negVert, float &negDist);
	void deleteInteriorTriangles(long triangle, std::set<trianglesUVW::edge,trianglesUVW::edgeTest> &borderEdges);
	void recurseTriangleSide(trianglesUVW* tri, int triangle, std::set<int> &triSide, std::vector<int> &edgeTriangles);
	bool getTriangleIntersect(const float (&position)[3],const float (&normal)[3], const float targetW, int &triangle,float (&UV)[2]);
	bool rayTestTriangle(trianglesUVW* tr, int triangleNumber, const float (&position)[3], const float (&normal)[3], float &rayParameter, float (&triUV)[2]);
	bool getTriEdgePath(trianglesUVW* tr, int vFrom, int vTo, Vec3f &plane, float d, std::vector<unsigned long> &triEdgePath, std::vector<float> &params);

	bool getTriEdgePath2(trianglesUVW* tr, int vFrom, int vTo, const Vec3f &surfaceNormal, std::vector<unsigned long> &triEdgePath, std::vector<float> &params);
	bool getTriEdgePathToTarget(trianglesUVW* tr, int vTo, const Vec3f &surfaceNormal, std::vector<unsigned long> &triEdgePath, std::vector<float> &params);
	bool topoCheck();

	void vertexNeighbors(int centerV, int vertex, float rSq, std::set<int> &nei);
	int getFirstStep(bool rejectFirst, const double (&plane)[3], const double D, const Vec3f &target, const Vec3f &dir, float &dirParam, int vertexTo, unsigned long &TE, float &param);
	bool topoStep(const double (&plane)[3], const double D, const Vec3f &target, const Vec3f &dir, float &dirParam, int vertexTo, unsigned long &TE, float &param, unsigned long &prevTE);
	bool segmentPath(int vertexFrom, int vertexTo, const Vec3f &surfaceNormal, std::vector<unsigned int> &vertexPath);
	float _incisionWidthOver2,_unmovedEdgeVertex[3],_topTx,_botTx,_TEdgeNormal[3];
	int _txCoord,_TinFrontVertices[2],_TinBackVertices[2],_ToutEndVertices[4];
	bool _endHcut;
	struct rayIntersect{
		trianglesUVW*	tr;
		int triNum;
		float triUV[2];
		rayIntersect(trianglesUVW*	tri,int triNumber,float (&triUv)[2])	{tr=tri; triNum=triNumber; triUV[0]=triUv[0]; triUV[1]=triUv[1];}
	};
	Vec3f _tPlaneNormal;
	float _tPlaneD;
};

#endif	// ifndef __INCISION_H__
