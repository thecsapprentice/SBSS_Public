// textures.cpp
// Author: Court Cutting with help from Richard Wright and others
// Date: January 22, 2012
// Purpose: Texture and bump map loader

#include <fstream>
#ifdef linux
#include <stdlib.h>
#endif
#include "Bitmap.h"
#include "textures.h"

textures::textures() {
}

textures::~textures() {
    clear();
}

void textures::clear() {	// clears textures from graphics card
    std::list<tex>::iterator tit;
    for(tit=_textures.begin(); tit!=_textures.end(); ++tit)	{
        //glDeleteTextures(1, &(tit->texture));
    }
    _textures.clear();
}

bool textures::textureExists(unsigned long textureNumber)
{
	std::list<tex>::iterator tit;
	for(tit=_textures.begin(); tit!=_textures.end(); ++tit)	{
		if((unsigned long)tit->texture==textureNumber)
			return true;
	}
	return false;
}

unsigned long textures::textureExists(std::string &textureName)
{	// 0 return is no texture found
	std::list<tex>::iterator tit;
	for(tit=_textures.begin(); tit!=_textures.end(); ++tit)	{
		if(tit->name==textureName)
			return (unsigned long)tit->texture;
	}
	return 0L;
}

unsigned long textures::firstTexture()
{
	if(_textures.empty())
		return 0L;
	return (unsigned long)_textures.front().texture;
}

unsigned long textures::loadTexture(const char *fileName)	{
	_textures.push_back(tex());
	_textures.back().name = fileName;
	//glGenTextures(1, &(_textures.back().texture));
	//glBindTexture(GL_TEXTURE_2D, _textures.back().texture);
	bool ret;
	if(_textures.back().name.size()-_textures.back().name.rfind(".bmp")==4)
		ret = loadBMPTexture(_textures.back().name.c_str(), GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR, GL_TEXTURE_WRAP_S);
	else if(_textures.back().name.size()-_textures.back().name.rfind(".tga")==4)
		ret = LoadTGATexture(_textures.back().name.c_str(), GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR, GL_TEXTURE_WRAP_S);
	else
		ret = false;
	if(!ret)	{
		//glDeleteTextures(1, &(_textures.back().texture));
		_textures.pop_back();
		return 0;
	}
	else
		return _textures.back().texture;
}

bool textures::loadBMPTexture(const char *fileName, GLenum minFilter, GLenum magFilter, GLenum wrapMode)
{
	Bitmap bmp;
	// Read the texture bits
	bmp.loadBMP(fileName);
	if(bmp.data == NULL) 
		return false;
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapMode);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapMode);
	
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minFilter);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magFilter);
    
	//glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, bmp.width, bmp.height, 0, GL_RGB, GL_UNSIGNED_BYTE, bmp.data);
    delete[] bmp.data ;
	bmp.data = NULL;
    if(minFilter == GL_LINEAR_MIPMAP_LINEAR || 
       minFilter == GL_LINEAR_MIPMAP_NEAREST ||
       minFilter == GL_NEAREST_MIPMAP_LINEAR ||
       minFilter == GL_NEAREST_MIPMAP_NEAREST)
        //glGenerateMipmap(GL_TEXTURE_2D);
	return true;
}

// Load a TGA as a 2D Texture. Completely initialize the state
bool textures::LoadTGATexture(const char *szFileName, GLenum minFilter, GLenum magFilter, GLenum wrapMode)
{
	GLbyte *pBits;
	int nWidth, nHeight, nComponents;
	GLenum eFormat;
	
	// Read the texture bits
	pBits = gltReadTGABits(szFileName, &nWidth, &nHeight, &nComponents, &eFormat);
	if(pBits == NULL) 
		return false;
	
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapMode);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapMode);
	
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minFilter);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magFilter);
    
	//glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	//glTexImage2D(GL_TEXTURE_2D, 0, nComponents, nWidth, nHeight, 0, eFormat, GL_UNSIGNED_BYTE, pBits);
	
    free(pBits);
    
    if(minFilter == GL_LINEAR_MIPMAP_LINEAR || 
       minFilter == GL_LINEAR_MIPMAP_NEAREST ||
       minFilter == GL_NEAREST_MIPMAP_LINEAR ||
       minFilter == GL_NEAREST_MIPMAP_NEAREST)
        //glGenerateMipmap(GL_TEXTURE_2D);
    
	return true;
}

////////////////////////////////////////////////////////////////////
// Allocate memory and load targa bits. Returns pointer to new buffer,
// height, and width of texture, and the OpenGL format of data.
// Call free() on buffer when finished!
// This only works on pretty vanilla targas... 8, 24, or 32 bit color
// only, no palettes, no RLE encoding.
GLbyte* textures::gltReadTGABits(const char *szFileName, GLint *iWidth, GLint *iHeight, GLint *iComponents, GLenum *eFormat)
{
    FILE *pFile;			// File pointer
    TGAHEADER tgaHeader;		// TGA file header
    unsigned long lImageSize;		// Size in bytes of image
    short sDepth;			// Pixel depth;
    GLbyte	*pBits = NULL;          // Pointer to bits
    // Default/Failed values
    *iWidth = 0;
    *iHeight = 0;
    *eFormat = GL_RGB;
    *iComponents = GL_RGB;
    // Attempt to open the file
    pFile = fopen(szFileName, "rb");
    if(pFile == NULL)
        return NULL;
    // Read in header (binary)
    fread(&tgaHeader, 18/* sizeof(TGAHEADER)*/, 1, pFile);
    // Do byte swap for big vs little endian
#ifdef __APPLE__
    LITTLE_ENDIAN_WORD(&tgaHeader.colorMapStart);
    LITTLE_ENDIAN_WORD(&tgaHeader.colorMapLength);
    LITTLE_ENDIAN_WORD(&tgaHeader.xstart);
    LITTLE_ENDIAN_WORD(&tgaHeader.ystart);
    LITTLE_ENDIAN_WORD(&tgaHeader.width);
    LITTLE_ENDIAN_WORD(&tgaHeader.height);
#endif
    // Get width, height, and depth of texture
    *iWidth = tgaHeader.width;
    *iHeight = tgaHeader.height;
    sDepth = tgaHeader.bits / 8;
    // Put some validity checks here. Very simply, I only understand
    // or care about 8, 24, or 32 bit targa's.
    if(tgaHeader.bits != 8 && tgaHeader.bits != 24 && tgaHeader.bits != 32)
        return NULL;
    // Calculate size of image buffer
    lImageSize = tgaHeader.width * tgaHeader.height * sDepth;
    // Allocate memory and check for success
    pBits = (GLbyte*)malloc(lImageSize * sizeof(GLbyte));
    if(pBits == NULL)
        return NULL;
    // Read in the bits
    // Check for read error. This should catch RLE or other 
    // weird formats that I don't want to recognize
    if(fread(pBits, lImageSize, 1, pFile) != 1)
		{
        free(pBits);
        return NULL;
		}
    // Set OpenGL format expected
    switch(sDepth)
		{
#ifndef OPENGL_ES
        case 3:     // Most likely case
            *eFormat = GL_BGR;
            *iComponents = GL_RGB;
            break;
#endif
        case 4:
            *eFormat = GL_BGRA;
            *iComponents = GL_RGBA;
            break;
        case 1:
            *eFormat = GL_LUMINANCE;
            *iComponents = GL_LUMINANCE;
            break;
        default:        // RGB
            // If on the iPhone, TGA's are BGR, and the iPhone does not 
            // support BGR without alpha, but it does support RGB,
            // so a simple swizzle of the red and blue bytes will suffice.
            // For faster iPhone loads however, save your TGA's with an Alpha!
        break;
		}
    // Done with File
    fclose(pFile);
    // Return pointer to image data
    return pBits;
}

