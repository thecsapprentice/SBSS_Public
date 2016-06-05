// File: sutures.h
// Author: Court Cutting
// Date: 2/20/12
// Purpose: Class for handling sutures and associated graphics

#ifndef __SUTURES_H__
#define __SUTURES_H__

#include <vector>
#include <map>
#include "GraphicsUtils/shapes.h"

// forward declarations
class trianglesUVW;
class GLmatrices;
class lightsShaders;

struct suture
{
    suture() {}
	~suture() {}
    
	int triangle0, triangle1;
    std::vector<float> uv0, uv1;
	int _cle_sutureid;		
};


#endif	// __SUTURES_H__
