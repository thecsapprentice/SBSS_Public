// textures.h
// Author: Court Cutting with help from Richard Wright and others
// Date: January 22, 2012
// Purpose: Texture and bump map loader
#ifndef __textures_h__
#define __textures_h__

// Bring in OpenGL 
// Windows
#ifdef WIN32
#include <windows.h>		// Must have for Windows platform builds
#ifndef GLEW_STATIC
#define GLEW_STATIC
#endif
#include "GL/glew.h"			// OpenGL Extension "autoloader"
//#include <gl\gl.h>			// Microsoft OpenGL headers (version 1.1 by themselves)
#endif

// Linux
#ifdef linux
#define GLEW_STATIC
#include <GL/glew.h>
#endif

#include <list>

class textures
{
public:
	unsigned long loadTexture(const char *fileName);	// 0 means load failure, otherwise is 2D texture buffer #
	bool textureExists(unsigned long textureNumber);
	unsigned long textureExists(std::string &textureName); // 0 return is no texture found
	unsigned long firstTexture();
	void clear();	// clears textures from graphics card
	textures();
	~textures();

private:
	struct tex{
		std::string name;
		GLuint texture;
	};
	std::list<tex> _textures;
	bool LoadTGATexture(const char *szFileName, GLenum minFilter, GLenum magFilter, GLenum wrapMode);
	bool loadBMPTexture(const char *fileName, GLenum minFilter, GLenum magFilter, GLenum wrapMode);
	GLbyte* gltReadTGABits(const char *szFileName, GLint *iWidth, GLint *iHeight, GLint *iComponents, GLenum *eFormat);
	// Define targa header. This is only used locally.
	#pragma pack(1)
	typedef struct
	{
		GLbyte	identsize;              // Size of ID field that follows header (0)
		GLbyte	colorMapType;           // 0 = None, 1 = paletted
		GLbyte	imageType;              // 0 = none, 1 = indexed, 2 = rgb, 3 = grey, +8=rle
		unsigned short	colorMapStart;          // First colour map entry
		unsigned short	colorMapLength;         // Number of colors
		unsigned char 	colorMapBits;   // bits per palette entry
		unsigned short	xstart;                 // image x origin
		unsigned short	ystart;                 // image y origin
		unsigned short	width;                  // width in pixels
		unsigned short	height;                 // height in pixels
		GLbyte	bits;                   // bits per pixel (8 16, 24, 32)
		GLbyte	descriptor;             // image descriptor
	} TGAHEADER;
	#pragma pack(8)
	struct RGB { 
		GLbyte blue;
		GLbyte green;
		GLbyte red;
		GLbyte alpha;
	};

};
#endif	// __textures_h__
