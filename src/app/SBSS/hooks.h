// File: hooks.h
// Author: Court Cutting
// Date: 2/18/12
// Purpose: Class for handling tissue hooks and associated graphics

#ifndef __HOOKS_H__
#define __HOOKS_H__

#include <map>
#include <vector>
#include "GraphicsUtils/shapes.h"

// forward declarations
class trianglesUVW;
class GLmatrices;
class shapes;

struct hook
{
public:
	hook() 
    {
        position.resize(3);
        uv.resize(2);
    }
	~hook() {}
    
    std::vector<float> position;
    std::vector<float> uv;
    
    int triangle;
    int _cle_hookid;
};

#endif	// __HOOKS_H__
