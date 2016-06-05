#ifndef __LINES_H__
#define __LINES_H__

#include <vector>
#include "sceneNode.h"

// forward declaration
class wxGraphics;

class lines : public sceneNode
{
public:
	void drawLines();
	// addLines() follows. Each element of points array has 4 floats as xyz1. lines uses 0xffffffff for
	// line strip primitive restart index.
	void addLines(const std::vector<GLfloat> &points, const std::vector<GLuint> &lines);
	void clear();
	bool empty() {return _linesVertexArrayBufferObject<1;	}
	void setWxGraphics(wxGraphics *wgg) {_wgg=wgg;}
	lines(void);
	~lines(void);

private:
	wxGraphics *_wgg;
	GLuint _linesBufferObjects[2];
	GLuint _linesVertexArrayBufferObject;
	GLsizei _nPoints;

};

#endif	// __LINES_H__