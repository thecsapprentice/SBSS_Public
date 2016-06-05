#include "wxGraphics.h"
#include "lines.h"

void lines::drawLines()
{
	glPrimitiveRestartIndex(0xffffffff);
	glEnable(GL_PRIMITIVE_RESTART);
	glBindVertexArray(_linesVertexArrayBufferObject);
	//	assumes glUseProgram(_program) has already been called
	glDrawElements(GL_LINE_STRIP, _nPoints, GL_UNSIGNED_INT, 0);
    // Unbind to anybody
	glBindVertexArray(0);
	glDisable(GL_PRIMITIVE_RESTART);
}

void lines::clear()
{
	if(!empty())	{
		if(_linesVertexArrayBufferObject)
			glDeleteVertexArrays(1,&_linesVertexArrayBufferObject);
		_linesVertexArrayBufferObject = 0;
		if(_linesBufferObjects[0])
		    glDeleteBuffers(2, _linesBufferObjects);
		_linesBufferObjects[0] = 0;
		_linesBufferObjects[1] = 0;
		_wgg->deleteSceneNode(this);
	}
}

void lines::addLines(const std::vector<GLfloat> &points, const std::vector<GLuint> &lines)
{	// each element of points array has 4 floats as xyz1. lines uses 0xffffffff for
	// line strip primitive restart index.
	if(lines.empty())
		return;
	if(!_linesVertexArrayBufferObject)
		glGenVertexArrays(1,&_linesVertexArrayBufferObject);
	if(!_linesBufferObjects[0])
	    glGenBuffers(2, _linesBufferObjects);
	setType(LINES);
	_nPoints = (GLsizei)lines.size();
    // Vertex and normal data
    glBindBuffer(GL_ARRAY_BUFFER, _linesBufferObjects[0]);	// VERTEX_DATA
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*points.size(), &(points[0]), GL_DYNAMIC_DRAW);
    // Indexes
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _linesBufferObjects[1]);	// INDEX_DATA
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*lines.size(), &(lines[0]), GL_STATIC_DRAW);

	// Create the master vertex array object
	glBindVertexArray(_linesVertexArrayBufferObject);
    // Vertex data
    glBindBuffer(GL_ARRAY_BUFFER, _linesBufferObjects[0]);	// VERTEX DATA
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
    // Indexes
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _linesBufferObjects[1]);	// INDEX_DATA
    // Unbind to anybody
	glBindVertexArray(0);
	// release for next use
    glBindBuffer( GL_ARRAY_BUFFER, 0);
    glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0);
	_wgg->addSceneNode(this);
}

lines::lines(void) : _linesVertexArrayBufferObject(0)
{
	_linesBufferObjects[0]=0;
	_linesBufferObjects[1]=0;
	_linesVertexArrayBufferObject = 0;

}


lines::~lines(void)
{
	if(_linesBufferObjects[0]>0) {
	    glDeleteBuffers(2, _linesBufferObjects);
		_linesBufferObjects[0]=0;
		_linesBufferObjects[1]=0; }
	if(_linesVertexArrayBufferObject>0) {
	    glDeleteVertexArrays(1, &_linesVertexArrayBufferObject);
		_linesVertexArrayBufferObject=0; }
}
