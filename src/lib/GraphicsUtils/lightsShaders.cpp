// file: lightsShaders.cpp
// Author: CourtCutting, MD
// Date: January 31,2012
// Purpose: This class does basic lighting and coordinated shaders.  It expects
//	vertex attributes vVertex(vec4), vNormal(vec3), and vTexture(vec2) as inputs.
//	It also expects uniforms mvpMatrix(mat4), mvMatrix(mat4), and mNormal(mat3)
//	which are the modelview-projection, modelview, and inverse modelview rotation
//	matrices respectively. It is expected in the future that this class will be
//	jazzed up significantly by a shader guru, but it at least provides a basic level.

#include <stdio.h>
#include <fstream>
#include <iostream>

#include "GLmatrices.h"
#include "lightsShaders.h"

GLuint lightsShaders::_textureProgram=0;
GLuint lightsShaders::_colorProgram=0;
GLuint lightsShaders::_lineProgram=0;

static const char *GTVertexShaderColoredLine = "#version 140\n"
	"in vec4 vVertex;\n"
	"uniform mat4   mvpMatrix;\n"
	"uniform mat4   mvMatrix;\n"
	"void main(void)\n"
	"{\n"
		// Get vertex position in eye coordinates
	"   vec4 vPosition4 = mvMatrix * vVertex;"
	"   gl_Position = mvpMatrix * vVertex;"
	"}";

static const char *GTFragmentShaderColoredLine = "#version 140\n"
	"out vec4 vFragColor;\n"
	"uniform vec4 objectColor;\n"
	"void main(void)\n"
	"{\n"
	"	vFragColor = objectColor;\n"
	"}";

static const char *GTVertexShaderColoredPhong = "#version 140\n"
	"in vec4 vVertex;\n"
	"in vec3 vNormal;\n"
	"uniform mat4   mvpMatrix;\n"
	"uniform mat4   mvMatrix;\n"
	"uniform mat3   normalMatrix;\n"
	"uniform vec3   vLightPosition;\n"
//	"uniform vec3   halfVector;\n"
	// half vector H=L+E which is the sum of the light vector and the eye vector
	// Blinn-Phong specular component S=(H dot N)^shininess_coef * light_spec_coef * material_spec_coef
	"smooth out vec3 normal,lightDir;\n"
//	"smooth out float dist;\n"
	"void main(void)\n"
	"{\n"
		// Get surface normal in eye coordinates
	"	normal = normalMatrix * vNormal;\n"
		// Get vertex position in eye coordinates
		/* these are the new lines of code to compute the light's direction */
	"   vec4 vPosition4 = mvMatrix * vVertex;"
	"	vec3 vPosition3 = vPosition4.xyz / vPosition4.w;\n"
	"	vec3 aux = vec3(vLightPosition-vPosition3);\n"
	"	lightDir = normalize(aux);\n"
//	"	dist = length(aux);\n"
	"   gl_Position = mvpMatrix * vVertex;"
	"}";

static const char *GTFragmentShaderColoredPhong = "#version 140\n"
	"smooth in vec3 normal,lightDir;\n"
	"out vec4 vFragColor;\n"
	"uniform vec4 objectColor;\n"
	// half vector H=L+E which is the sum of the light vector and the eye vector
	// Blinn-Phong specular component S=(H dot N)^shininess_coef * light_spec_coef * material_spec_coef
//	"uniform vec3   halfVector;\n"	// make sure is already normalized
	"void main(void)\n"
	"{\n"
	"	vec3 n;\n"
	"	float NdotL,NdotHV;\n"
//	"	float att;\n"
	"	vFragColor = objectColor;\n"
	"	vFragColor.rgb *= 0.3f;\n"
	"	n = normalize(normal);\n"
	// compute the dot product between normal and normalized lightdir
	"	NdotL = max(dot(n,normalize(lightDir)),0.0);\n"
	"	if (NdotL > 0.0) {\n"
//			att = 1.0 / (gl_LightSource[0].constantAttenuation +
//					gl_LightSource[0].linearAttenuation * dist +
//					gl_LightSource[0].quadraticAttenuation * dist * dist);
//			vFragColor += att * (diffuse * NdotL + ambient);
"		vFragColor.rgb += NdotL * objectColor.rgb*0.7f;\n"
//	"		NdotHV = max(dot(n,halfVector),0.0);\n"
//			color += att * gl_FrontMaterial.specular * gl_LightSource[0].specular * 
//							pow(NdotHV,gl_FrontMaterial.shininess);
	"		float fSpec = pow(NdotL, 1.0);\n"
	"		vFragColor.rgb += vec3(.5, .5, .5) * fSpec; }\n"
	"}";

static const char *GTVertexShaderDefault = "#version 140\n"
	// Incoming per vertex
	"in vec4 vVertex;\n"
	"in vec3 vNormal;\n"
	"in vec2 vTexture;\n"
	"uniform mat4   mvpMatrix;\n"
	"uniform mat4   mvMatrix;\n"
	"uniform mat3   normalMatrix;\n"
	"uniform vec3   vLightPosition;\n"
	// Color to fragment program
	"smooth out vec3 vVaryingNormal;\n"
	"smooth out vec3 vVaryingLightDir;\n"
	"smooth out vec2 vTexCoords;\n"
	"void main(void)\n"
	"{\n"
		// Get surface normal in eye coordinates
	"	vVaryingNormal = normalMatrix * vNormal;\n"
		// Get vertex position in eye coordinates
	"   vec4 vPosition4 = mvMatrix * vVertex;"
	"	vec3 vPosition3 = vPosition4.xyz / vPosition4.w;\n"
		// Get vector to light source"
	"	vVaryingLightDir = normalize(vLightPosition-vPosition3);\n"
		// Pass along the texture coordinates"
	"	vTexCoords = vTexture.st;\n"
	"   gl_Position = mvpMatrix * vVertex;"
	"}";

static const char *GTFragmentShaderDefault = // ADS Point lighting Shader
	// Adapted from Richard S. Wright Jr.
	// OpenGL SuperBible
	"#version 140\n"
	"out vec4 vFragColor;\n"
	"uniform vec4 ambientColor;\n"
	"uniform vec4 diffuseColor;\n"
//	"uniform vec4 specularColor;\n"
	"uniform sampler2D colorMap;\n"
	"smooth in vec3 vVaryingNormal;\n"
	"smooth in vec3 vVaryingLightDir;\n"
	"smooth in vec2 vTexCoords;\n"
	"void main(void)\n"
	"{\n"
		// Dot product gives us diffuse intensity
	"	float diff = max(0.0, dot(normalize(vVaryingNormal), normalize(vVaryingLightDir)));\n"
		// Multiply intensity by diffuse color, force alpha to 1.0
		// Add in ambient light"
	"	vFragColor = ambientColor;\n"
	"	vFragColor += diff * diffuseColor;\n"
		// Modulate in the texture
	"	vFragColor *= texture(colorMap, vTexCoords);\n"
		// Specular Light
//	"	vec3 vReflection = normalize(reflect(-normalize(vVaryingLightDir), normalize(vVaryingNormal)));\n"
//	"	float spec = max(0.0, dot(normalize(vVaryingNormal), vVaryingLightDir));\n"
	"	if(diff > 0) {\n"
	"		float fSpec = pow(diff, 128.0);\n"
	"		vFragColor.rgb += vec3(.1, .1, .1) * fSpec;\n"
	"	}\n"
	"	if(vTexCoords.s > 1.5)"
	"		vFragColor.rgb = vec3(0.0, 1.0, 0.0);"
	"}";

void lightsShaders::setModelMatrix(GLfloat *model)
{
/*
	// remember setProjectionMatrix() and setFrameAndRotation() must be called first in GLmatrices or will get undefined results
	const GLfloat *fr=_glM->getFrameAndRotationMatrix();
	const GLfloat *pm=_glM->getProjectionMatrix();
	for(int i=0; i<4; ++i)	{
		for(int j=0; j<4; ++j)	{	// model happens first, then frame-rotation
			_modelMat[(i<<2)+j] = fr[j]*model[i<<2] + fr[4+j]*model[(i<<2)+1] + fr[8+j]*model[(i<<2)+2] + fr[12+j]*model[(i<<2)+3];
			if(i<3 && j<3)
				_normMat[(i<<1)+i+j]=_modelMat[(i<<2)+j];
		}
	}
	for(int i=0; i<4; ++i)	{ //columns
		for(int j=0; j<4; ++j)	{ // rows. Model first then passed into projection matrix
			_MVP[(i<<2)+j] = pm[j]*_modelMat[i<<2] + pm[4+j]*_modelMat[(i<<2)+1] + pm[8+j]*_modelMat[(i<<2)+2] + pm[12+j]*_modelMat[(i<<2)+3];
		}
	}
	if(_textureProgram>0) {
		glUniformMatrix4fv(_texUni.locMVP, 1, GL_FALSE, (GLfloat *)_MVP);
		glUniformMatrix4fv(_texUni.locMV, 1, GL_FALSE, (GLfloat *)_modelMat);
		glUniformMatrix3fv(_texUni.locNM, 1, GL_FALSE, (GLfloat *)_normMat); }
	if(_colorProgram>0) {
		glUniformMatrix4fv(_colorUni.locMVP, 1, GL_FALSE, (GLfloat *)_MVP);
		glUniformMatrix4fv(_colorUni.locMV, 1, GL_FALSE, (GLfloat *)_modelMat);
		glUniformMatrix3fv(_colorUni.locNM, 1, GL_FALSE, (GLfloat *)_normMat); }
	if(_lineProgram>0) {
		glUniformMatrix4fv(_lineUni.locMVP, 1, GL_FALSE, (GLfloat *)_MVP);
		glUniformMatrix4fv(_lineUni.locMV, 1, GL_FALSE, (GLfloat *)_modelMat);} 
	progUniforms &pu = _programUniforms[_currentProgram];
    
	if(pu.locMVP>-1)
		glUniformMatrix4fv(pu.locMVP, 1, GL_FALSE, (GLfloat *)_MVP);
	if(pu.locMV>-1)
		glUniformMatrix4fv(pu.locMV, 1, GL_FALSE, (GLfloat *)_modelMat);
	if(pu.locNM>-1)
		glUniformMatrix3fv(pu.locNM, 1, GL_FALSE, (GLfloat *)_normMat);
    
    */
}

/* bool lightsShaders::useLineProgram()
{
	if(_lineProgram<1)	{
		if(!createLineProgram())
			return false;
	}
	glUseProgram(_lineProgram);
	return true;
}

bool lightsShaders::useColorProgram()
{
	if(_colorProgram<1)
		return false;
	glUseProgram(_colorProgram);
	glUniform3fv(_colorUni.locLight, 1, _vEyeLight);
	return true;
} */

void lightsShaders::setColor(GLfloat *color)
{
	// COURT - add back later
/*
	if(_colorProgram>0)
		glUniform4fv(_colorUni.locObjColor, 1, (GLfloat *)color);
	if(_lineProgram>0)
		glUniform4fv(_lineUni.color, 1, (GLfloat *)color);
	glUseProgram(_colorProgram);
	glUniform4fv(glGetUniformLocation(_colorProgram, "objectColor"), 1, (GLfloat *)color);
*/
}

void lightsShaders::useGlslProgram(GLuint programNumber)  // careful - no error checking for validity
{
/*
	glUseProgram(programNumber);
	_currentProgram = programNumber;
	setProgramUniforms(programNumber);
*/
}

/* bool lightsShaders::useTextureProgram()
{
	// first do textured objects
	if(_textureProgram<1)
		return false;
//		glBindTexture(GL_TEXTURE_2D, texture);
	glUseProgram(_textureProgram);
	glUniform4fv(_texUni.locAmbient, 1, _vAmbientColor);
	glUniform4fv(_texUni.locDiffuse, 1, _vDiffuseColor);
	glUniform4fv(_texUni.locSpecular, 1, _vSpecularColor);
	glUniform3fv(_texUni.locLight, 1, _vEyeLight);
	// 2D texture
	glUniform1i(_texUni.locTexture, 0);
	return true;
} */

bool lightsShaders::createLineProgram()
{
/*
	if(_lineProgram>0)
		return true;
	std::vector<std::string> att;
	att.assign(1,std::string());
	att[0] = "vVertex";
	if(createProgramWithAttributes(_lineProgram,GTVertexShaderColoredLine,GTFragmentShaderColoredLine,att))	{
		_lineUni.locMVP = glGetUniformLocation(_lineProgram, "mvpMatrix");
		_lineUni.locMV  = glGetUniformLocation(_lineProgram, "mvMatrix");
		_lineUni.color  = glGetUniformLocation(_lineProgram, "objectColor");
		return true;
	}
	else
		return false;
*/
    return true;
}

GLuint lightsShaders::getOrCreateColorProgram()
{
/*
	if(_colorProgram>0)
		return _colorProgram;
	std::vector<std::string> att;
	att.assign(2,std::string());
	att[0] = "vVertex";
	att[1] = "vNormal";
	createProgramWithAttributes(_colorProgram,GTVertexShaderColoredPhong,GTFragmentShaderColoredPhong,att);
	_colorUni.locMVP = glGetUniformLocation(_colorProgram, "mvpMatrix");
	_colorUni.locMV  = glGetUniformLocation(_colorProgram, "mvMatrix");
	_colorUni.locNM  = glGetUniformLocation(_colorProgram, "normalMatrix");
	_colorUni.locLight = glGetUniformLocation(_colorProgram, "vLightPosition");
	_colorUni.locObjColor  = glGetUniformLocation(_colorProgram, "objectColor"); 
	if(!_programUniforms.insert(std::make_pair(_colorProgram,progUniforms())).second) {
		glDeleteProgram(_colorProgram);
		return 0;
	}
	glUseProgram(_colorProgram);
	progUniforms *pu = &_programUniforms[_colorProgram];
	pu->notDoneOnce = true;
	pu->locAmbient = glGetUniformLocation(_colorProgram, "ambientColor");
	pu->locDiffuse = glGetUniformLocation(_colorProgram, "diffuseColor");
	pu->locSpecular = glGetUniformLocation(_colorProgram, "specularColor");
	pu->locMVP = glGetUniformLocation(_colorProgram, "mvpMatrix");
	pu->locMV  = glGetUniformLocation(_colorProgram, "mvMatrix");
	pu->locNM  = glGetUniformLocation(_colorProgram, "normalMatrix");
	pu->locLight = glGetUniformLocation(_colorProgram, "vLightPosition");
	pu->locObjColor = glGetUniformLocation(_colorProgram, "objectColor");
	pu->locTexture0 = glGetUniformLocation(_colorProgram, "colorMap");
	pu->locTexture1 = glGetUniformLocation(_colorProgram, "texture1");
	pu->locTexture2 = glGetUniformLocation(_colorProgram, "texture2");
	return _colorProgram;
*/
    return 0;
}

void lightsShaders::createTextureProgram()
{
/*
	if(_textureProgram>0)
		return;
	std::vector<std::string> att;
	att.assign(3,std::string());
	att[0] = "vVertex";
	att[1] = "vNormal";
	att[2] = "vTexture";
	createProgramWithAttributes(_textureProgram,GTVertexShaderDefault,GTFragmentShaderDefault,att);
	_texUni.locAmbient = glGetUniformLocation(_textureProgram, "ambientColor");
	_texUni.locDiffuse = glGetUniformLocation(_textureProgram, "diffuseColor");
	_texUni.locSpecular = glGetUniformLocation(_textureProgram, "specularColor");
	_texUni.locMVP = glGetUniformLocation(_textureProgram, "mvpMatrix");
	_texUni.locMV  = glGetUniformLocation(_textureProgram, "mvMatrix");
	_texUni.locNM  = glGetUniformLocation(_textureProgram, "normalMatrix");
	_texUni.locLight = glGetUniformLocation(_textureProgram, "vLightPosition");
	_texUni.locTexture = glGetUniformLocation(_textureProgram, "colorMap"); 
	if(!_programUniforms.insert(std::make_pair(_textureProgram,progUniforms())).second) {
		glDeleteProgram(_textureProgram);
		return;
	}
	progUniforms *pu = &_programUniforms[_textureProgram];
	pu->notDoneOnce = true;
	pu->locAmbient = glGetUniformLocation(_textureProgram, "ambientColor");
	pu->locDiffuse = glGetUniformLocation(_textureProgram, "diffuseColor");
	pu->locSpecular = glGetUniformLocation(_textureProgram, "specularColor");
	pu->locMVP = glGetUniformLocation(_textureProgram, "mvpMatrix");
	pu->locMV  = glGetUniformLocation(_textureProgram, "mvMatrix");
	pu->locNM  = glGetUniformLocation(_textureProgram, "normalMatrix");
	pu->locLight = glGetUniformLocation(_textureProgram, "vLightPosition");
	pu->locObjColor = glGetUniformLocation(_textureProgram, "objectColor");
	pu->locTexture0 = glGetUniformLocation(_textureProgram, "colorMap");
	pu->locTexture1 = glGetUniformLocation(_textureProgram, "texture1");
	pu->locTexture2 = glGetUniformLocation(_textureProgram, "texture2");
*/
}

bool lightsShaders::createCustomProgram(GLuint &program, const char *vertexShader, const char *fragmentShader, std::vector<std::string> &attributes)
{
/*
	if(!createProgramWithAttributes(program,vertexShader,fragmentShader,attributes))
		return false;
	if(!_programUniforms.insert(std::make_pair(program,progUniforms())).second) {
		glDeleteProgram(program);
		return false;
	}
	progUniforms *pu = &_programUniforms[program];
	pu->notDoneOnce = true;
	pu->locAmbient = glGetUniformLocation(program, "ambientColor");
	pu->locDiffuse = glGetUniformLocation(program, "diffuseColor");
	pu->locSpecular = glGetUniformLocation(program, "specularColor");
	pu->locMVP = glGetUniformLocation(program, "mvpMatrix");
	pu->locMV  = glGetUniformLocation(program, "mvMatrix");
	pu->locNM  = glGetUniformLocation(program, "normalMatrix");
	pu->locLight = glGetUniformLocation(program, "vLightPosition");
	pu->locObjColor = glGetUniformLocation(program, "objectColor");
	pu->locTexture0 = glGetUniformLocation(program, "colorMap");
	pu->locTexture1 = glGetUniformLocation(program, "normalMap");
	pu->locTexture2 = glGetUniformLocation(program, "texture2");
*/
	return true;

}

void lightsShaders::setProgramUniforms(GLuint program) {
/*
	progUniforms *pu = &_programUniforms[program];
	if(pu->notDoneOnce) {
		if(pu->locAmbient>-1)
			glUniform4fv(pu->locAmbient, 1, _vAmbientColor);
		if(pu->locDiffuse>-1)
			glUniform4fv(pu->locDiffuse, 1, _vDiffuseColor);
		if(pu->locSpecular>-1)
			glUniform4fv(pu->locSpecular, 1, _vSpecularColor);
		if(pu->locLight>-1)
			glUniform3fv(pu->locLight, 1, _vEyeLight);
		// 2D texture
		if(pu->locTexture0>-1)
			glUniform1i(pu->locTexture0, 0);
		if(pu->locTexture1>-1)
			glUniform1i(pu->locTexture1, 1);
		if(pu->locTexture2>-1)
			glUniform1i(pu->locTexture2, 2);
		pu->notDoneOnce = false;
	}
	if(pu->locMVP>0)
		glUniformMatrix4fv(pu->locMVP, 1, GL_FALSE, (GLfloat *)_MVP);
	if(pu->locMV>0)
		glUniformMatrix4fv(pu->locMV, 1, GL_FALSE, (GLfloat *)_modelMat);
	if(pu->locNM>0)
		glUniformMatrix3fv(pu->locNM, 1, GL_FALSE, (GLfloat *)_normMat);
//	if(pu->locObjColor>-1)
//		glUniform4fv(pu->locObjColor, 1, _vDiffuseColor);	// COURT - change me
*/
}

bool lightsShaders::createProgramWithAttributes(GLuint &program, const char *vertexShader, const char *fragmentShader, std::vector<std::string> &attributes)
{
/*
	if(program>0)
		return true;
	// Temporary Shader objects
	GLuint hVertexShader;
	GLuint hFragmentShader;
	GLint testVal;
	// Create shader objects
	hVertexShader = glCreateShader(GL_VERTEX_SHADER);
	hFragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	GLchar *fsStringPtr[1];
	fsStringPtr[0] = (GLchar *)vertexShader;
	glShaderSource(hVertexShader, 1, (const GLchar **)fsStringPtr, NULL);
	fsStringPtr[0] = (GLchar *)fragmentShader;
	glShaderSource(hFragmentShader, 1, (const GLchar **)fsStringPtr, NULL);
	// Compile them
	glCompileShader(hVertexShader);
	glCompileShader(hFragmentShader);
	// Check for errors
	glGetShaderiv(hVertexShader, GL_COMPILE_STATUS, &testVal);
	if(testVal == GL_FALSE){
		GLchar infoLog[2000];
		GLsizei infoLength;
		glGetShaderInfoLog(hVertexShader,2000,&infoLength,infoLog);
		printf("%s",infoLog);
		glDeleteShader(hVertexShader);
		glDeleteShader(hFragmentShader);
		return false;
	}
	glGetShaderiv(hFragmentShader, GL_COMPILE_STATUS, &testVal);
	if(testVal == GL_FALSE){
		GLchar infoLog[2000];
		GLsizei infoLength;
		glGetShaderInfoLog(hFragmentShader,2000,&infoLength,infoLog);
		printf("%s",infoLog);
		glDeleteShader(hVertexShader);
		glDeleteShader(hFragmentShader);
		return false;
	}
	// Link them - assuming it works...
	program = glCreateProgram();
	glAttachShader(program, hVertexShader);
	glAttachShader(program, hFragmentShader);
	// List of attributes
	int i,n=(int)attributes.size();
	for(i=0; i<n; ++i)
		glBindAttribLocation(program, i, attributes[i].c_str());
//	glBindAttribLocation(_textureProgram, 1, "vNormal");
//	glBindAttribLocation(_textureProgram, 2, "vTexture0");
	glLinkProgram(program);
	// These are no longer needed
	glDeleteShader(hVertexShader);
	glDeleteShader(hFragmentShader);  
	// Make sure link worked too
	glGetProgramiv(program, GL_LINK_STATUS, &testVal);
	if(testVal == GL_FALSE){
		GLchar infoLog[2000];
		GLsizei infoLength;
		glGetProgramInfoLog(program,2000,&infoLength,infoLog);
		printf("%s",infoLog);
		glDeleteProgram(program);
		return false;
	}
*/
	return true;
}

int loadShader(const char** dest, const char* filename, const char* fallback ){
    std::ifstream shader_reader;
    shader_reader.exceptions ( std::ifstream::failbit | std::ifstream::badbit );

    std::cerr << "Reading shader " << filename << std::endl;
    int success = 0;

    char* new_buff = NULL;
    try{
        shader_reader.open( filename );
        if( shader_reader.fail() )
            std::cerr << "Uncaught failure detected!!!!" << std::endl;

        shader_reader.seekg (0, shader_reader.end);
        int length = (int)shader_reader.tellg();
        shader_reader.seekg (0, shader_reader.beg);

        new_buff = new char[length+1];
        new_buff[length] = '\0';

        shader_reader.read( new_buff, length);

        (*dest) = new_buff;        
        success = 1;
        shader_reader.close();
    }
    catch (std::ifstream::failure e) {
        std::cerr << "Could not read " << filename << ", using fallback shader.\n";
        if(new_buff!=NULL)
            delete [] new_buff;
        (*dest) = fallback;
    } 
    return success;
}


lightsShaders::lightsShaders()
{
/*
    _programUniforms.clear();
    if(loadShader( &__GTVertexShaderColoredLine, "shaders/GTVertexShaderColoredLine.shader", GTVertexShaderColoredLine))
        cleanup_shaders.push_back( __GTVertexShaderColoredLine );
    if(loadShader( &__GTVertexShaderColoredPhong, "shaders/GTVertexShaderColoredPhong.shader", GTVertexShaderColoredPhong))
        cleanup_shaders.push_back( __GTVertexShaderColoredPhong );
    if(loadShader( &__GTVertexShaderDefault, "shaders/GTVertexShaderDefault.shader", GTVertexShaderDefault))
        cleanup_shaders.push_back( __GTVertexShaderDefault );
    if(loadShader( &__GTFragmentShaderColoredLine, "shaders/GTFragmentShaderColoredLine.shader", GTFragmentShaderColoredLine))
        cleanup_shaders.push_back( __GTFragmentShaderColoredLine );
    if(loadShader( &__GTFragmentShaderColoredPhong, "shaders/GTFragmentShaderColoredPhong.shader", GTFragmentShaderColoredPhong))
        cleanup_shaders.push_back( __GTFragmentShaderColoredPhong );
    if(loadShader( &__GTFragmentShaderDefault, "shaders/GTFragmentShaderDefault.shader", GTFragmentShaderDefault))
        cleanup_shaders.push_back( __GTFragmentShaderDefault ); 
    
  
	_textureProgram = 0;
	_colorProgram = 0;
	_vEyeLight[0]=0.0f; _vEyeLight[1]=0.0f; _vEyeLight[2]=400.0f;
	_vAmbientColor[0]=0.2f; _vAmbientColor[1]=0.2f; _vAmbientColor[2]=0.2f; _vAmbientColor[3]=1.0f;
	_vDiffuseColor[0]=0.8f; _vDiffuseColor[1]=0.8f; _vDiffuseColor[2]=0.8f; _vDiffuseColor[3]=1.0f;
	_vSpecularColor[0]=1.0f; _vSpecularColor[1]=1.0f; _vSpecularColor[2]=1.0f; _vSpecularColor[3]=1.0f;
*/
}

lightsShaders::~lightsShaders()
{
//    for( int i = 0; i < cleanup_shaders.size(); i++)
//        delete [] cleanup_shaders[i];
//    clear();
}

void lightsShaders::clear() {
//	if(_textureProgram>0)
//		glDeleteProgram(_textureProgram);
//	if(_colorProgram>0)
//		glDeleteProgram(_colorProgram);
//	if(_lineProgram>0)
//		glDeleteProgram(_lineProgram);
//	std::map<GLuint,progUniforms>::iterator pit;
//	for(pit=_programUniforms.begin(); pit!=_programUniforms.end(); ++pit)
//		glDeleteProgram(pit->first);
}

