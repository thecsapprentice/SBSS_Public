//#extension GL_OES_standard_derivatives : enable
#define USE_MAP
#define USE_NORMALMAP

// Description : Array and textureless GLSL 2D & 3D simplex noise functions.
//      Author : Ian McEwan, Ashima Arts.
//  Maintainer : ijm
//     Lastmod : 20110822 (ijm)
//     License : Copyright (C) 2011 Ashima Arts. All rights reserved.
//               Distributed under the MIT License. See LICENSE file.
//               https://github.com/ashima/webgl-noise


vec3 mod289(vec3 x) {
    return x - floor(x * (1.0 / 289.0)) * 289.0;}
vec2 mod289(vec2 x) {
    return x - floor(x * (1.0 / 289.0)) * 289.0;}
vec3 permute(vec3 x) {
    return mod289(((x*34.0)+1.0)*x);}
vec4 taylorInvSqrt(vec4 r){
    return 1.79284291400159 - 0.85373472095314 * r;}
float snoise(vec2 v)  {
    const vec4 C = vec4(0.211324865405187,  // (3.0-sqrt(3.0))/6.0
                        0.366025403784439,  // 0.5*(sqrt(3.0)-1.0)
                        -0.577350269189626,  // -1.0 + 2.0 * C.x
                        0.024390243902439); // 1.0 / 41.0
    // First corner
    vec2 i  = floor(v + dot(v, C.yy) );
    vec2 x0 = v -   i + dot(i, C.xx);
    // Other corners
    vec2 i1;
    //i1.x = step( x0.y, x0.x ); // x0.x > x0.y ? 1.0 : 0.0
    //i1.y = 1.0 - i1.x;
    i1 = (x0.x > x0.y) ? vec2(1.0, 0.0) : vec2(0.0, 1.0);
    // x0 = x0 - 0.0 + 0.0 * C.xx ;
    // x1 = x0 - i1 + 1.0 * C.xx ;
    // x2 = x0 - 1.0 + 2.0 * C.xx ;
    vec4 x12 = x0.xyxy + C.xxzz;
    x12.xy -= i1;
    // Permutations
    i = mod289(i); // Avoid truncation effects in permutation
    vec3 p = permute( permute( i.y + vec3(0.0, i1.y, 1.0 ))
    	              + i.x + vec3(0.0, i1.x, 1.0 ));
    vec3 m = max(0.5 - vec3(dot(x0,x0), dot(x12.xy,x12.xy), dot(x12.zw,x12.zw)), 0.0);
    m = m*m ;
    m = m*m ;
    // Gradients: 41 points uniformly over a line, mapped onto a diamond.
    // The ring size 17*17 = 289 is close to a multiple of 41 (41*7 = 287)
    vec3 x = 2.0 * fract(p * C.www) - 1.0;
    vec3 h = abs(x) - 0.5;
    vec3 ox = floor(x + 0.5);
    vec3 a0 = x - ox;
    // Normalise gradients implicitly by scaling m
    // Approximation of: m *= inversesqrt( a0*a0 + h*h );
    m *= 1.79284291400159 - 0.85373472095314 * ( a0*a0 + h*h );
    // Compute final noise value at P
    vec3 g;
    g.x  = a0.x  * x0.x  + h.x  * x0.y;
    g.yz = a0.yz * x12.xz + h.yz * x12.yw;
    return 130.0 * dot(m, g);}

// Description : Custom skin fragment shader
//      Author : Court Cutting, MD
//        Date : 20140604
//     License : Copyright (C) 2014 Court Cutting. All rights reserved.
//               Distributed under the GNU General Public License. See LICENSE file.
//               https://gnu.org/licenses/gpl.html
// next routine gets fat from 4 surrounding noise values [-1,1]
void getFat(in vec4 nei, out vec4 fatColor, out vec3 normDelta, out float specMult) {
    float h;
    for(int i=0; i<4; ++i){
        h = 1.0 - abs(nei[i]);
        h *= h;
        nei[i] = 1.0 - h; }
    h = 0.0;
    for(int i=0; i<4; ++i)
        h += nei[i];
    h *= 0.25;
    vec2 p;
    p.x = nei[1]-nei[0];
    p.y = nei[2]-nei[3];
    p *= 130.0;
    p = clamp(p,-1.0,1.0);
    float d,f;
    d = dot(p,p);
    f = inversesqrt(d+1.0);
    p.x = -p.x;
    normDelta = vec3(p,1.0)*f;
    float fatRed, fatGreen, fatBlue;
    if(h<0.04) {
        specMult = 0.2;
        fatRed = (1.0-h)*0.4;
        fatBlue = 0.0;
        fatGreen = (1.0-h)*0.2; }
    else {
        specMult = 1.0;
        fatRed = 0.5 + h*0.8;
        fatBlue = 0.15 + h*0.3;
        fatGreen = 0.35 + h*0.8; }
    fatColor = vec4(fatRed, fatGreen, fatBlue, 1.0); }


uniform float dirIntensity;
uniform vec3 diffuse;
uniform vec3 ambient;
uniform vec3 specular;
uniform float shininess;

varying vec2 vUv;
uniform sampler2D map;
//uniform vec3 ambientLightColor;
varying vec3 vViewPosition;
varying vec3 vNormal;
//varying vec3 vData;
uniform float dataOverlayPercent;
uniform vec2 dataOverlayStrainRange;
uniform vec2 dataOverlayStressRange;
uniform int dataOverlayUseMagnitude;
uniform int dataOverlaySelector;
uniform sampler2D dataColorRamp;
uniform sampler2D dataTexture;

#ifdef USE_NORMALMAP
uniform sampler2D normalMap;
uniform vec2 normalScale;

vec3 perturbNormal2Arb( vec3 eye_pos, vec3 surf_norm, vec3 mapN ) {
    vec3 q0 = dFdx( eye_pos.xyz );
    vec3 q1 = dFdy( eye_pos.xyz );
    vec2 st0 = dFdx( vUv.st );
    vec2 st1 = dFdy( vUv.st );
    vec3 S = normalize(  q0 * st1.t - q1 * st0.t );
    vec3 T = normalize( -q0 * st1.s + q1 * st0.s );
    vec3 N = normalize( surf_norm );
    //		vec3 mapN = texture2D( normalMap, vUv ).xyz * 2.0 - 1.0;
    mapN.xy = normalScale * mapN.xy;
    mat3 tsn = mat3( S, T, N );
    return normalize( tsn * mapN );
}
#endif
void main() {
    bool useDoobsPhong = false;
    const float fatIncr = 0.5/1024.0;
    vec2 fatD = vec2(fatIncr,0.0);
    vec2 faceUV;
    float lightVal;
    vec3 normDelta = vec3(0.0, 0.0, 1.0);
    float sn,h,specMult = 1.0;
    float pShininess = shininess;
    vec4 nei,vFragColor;
    if(vUv.s>1.05) {
        if(vUv.s<2.15) {  // middle
            faceUV = vUv - vec2(1.1,0.0);
            sn = snoise(vec2(faceUV.t*3.0,0.5));
            float val = (sn+1.0)*0.5;
            if(0.40 + 0.2*val > faceUV.s) {  // fat
                faceUV *= vec2(6.4,4.7);
                nei[0] = snoise(faceUV + fatD);
                nei[1] = snoise(faceUV - fatD);
                nei[3] = snoise(faceUV - fatD.yx);
                nei[2] = snoise(faceUV + fatD.yx);
                getFat(nei,vFragColor,normDelta,specMult);
                pShininess *= 4.0;
                useDoobsPhong = true;
            }
            else if(0.415 + 0.2*val > faceUV.s) {  // dermal-fat junction
                vFragColor = vec4(0.51, 0.44, 0.1412, 1.0);
                specMult = 0.2; }
            else {  // dermis
                specMult = 0.35;
                if(faceUV.s>0.95)
                    vFragColor = vec4(0.71, 0.57255, 0.2784, 1.0);
                else if(faceUV.s>0.93)
                    vFragColor = vec4(0.843, 0.737, 0.51, 1.0);
                else {
                    sn = snoise(vec2((faceUV.t+0.4)*4.2,0.5));
                    val = (sn-0.5)*2.0;
                    if(0.85 + 0.05*val > faceUV.s)
                        vFragColor = vec4(0.9569, 0.8902, 0.71, 1.0);
                    else {
                        vFragColor = vec4(0.7255, 0.5059, 0.2039, 1.0);
                        vFragColor = vFragColor*(3.3-2.8*faceUV.s); }
                } }
        }
        else {  // bottom
            faceUV = vUv - vec2(2.2,0.0);
            nei[0] = snoise(faceUV - fatD);
            nei[1] = snoise(faceUV + fatD);
            nei[2] = snoise(faceUV - fatD.yx);
            nei[3] = snoise(faceUV + fatD.yx);
            getFat(nei,vFragColor,normDelta,specMult);
            useDoobsPhong = true;
        }
        pShininess *= 4.0;
        gl_FragColor = vFragColor;
    }
    else {
        //	vec4 texelColor = texture2D( map, vUv );
        gl_FragColor = texture2D( map, vUv );        
	normDelta = texture2D( normalMap, vUv ).xyz * 2.0 - 1.0;
        useDoobsPhong = true;
    } 

    if(useDoobsPhong) {

	vec3 normal = normalize( vNormal );
	vec3 viewPosition = normalize( vViewPosition );
#ifdef USE_NORMALMAP
	normal = perturbNormal2Arb( -vViewPosition, normal, normDelta );
#endif
	vec3 dirDiffuse  = vec3( 0.0 );
	vec3 dirSpecular = vec3( 0.0 );
        //		vec4 lDirection = viewMatrix * vec4( directionalLightDirection[ 0 ], 0.0 );
	vec4 lDirection = viewMatrix * vec4( cameraPosition, 0.0 ); // vec4( 0.0, 0.0, 1.0, 0.0 );
	vec3 dirVector = normalize( lDirection.xyz );
	float dotProduct = dot( normal, dirVector );
	float dirDiffuseWeight = max( dotProduct, 0.0 );
	dirDiffuse  += diffuse * dirDiffuseWeight * dirIntensity;
	vec3 dirHalfVector = normalize( dirVector + viewPosition );
	float dirDotNormalHalf = max( dot( normal, dirHalfVector ), 0.0 );
	float dirSpecularWeight = specMult * max( pow( dirDotNormalHalf, pShininess ), 0.0 );
	float specularNormalization = ( pShininess + 2.0001 ) / 8.0;
	vec3 schlick = specular + vec3( 1.0 - specular ) * pow( max( 1.0 - dot( dirVector, dirHalfVector ), 0.0 ), 5.0 );
	dirSpecular += schlick * dirSpecularWeight * dirDiffuseWeight * specularNormalization * dirIntensity;
	gl_FragColor.xyz = gl_FragColor.xyz * ( dirDiffuse + ambient ) + dirSpecular;
    }


    vec3 source;
    if( vUv.y > 1.0 || vUv.y < 0.0 || vUv.x > 1.0 || vUv.x < 0.0 )
        source = vec3(0.0,0.0,0.0);
    else{
        source = texture2D(dataTexture, vec2(vUv.y,1.0-vUv.x) ).rgb;
    }
    
    vec2 range;
    //if( dataOverlaySelector == 0 )
    //    range = dataOverlayStrainRange;
    //else
    //    range = dataOverlayStressRange;

    range = dataOverlayStrainRange;
    float pow_factor = dataOverlayStressRange.y;

    // Mix with stress information
    if( dataOverlayUseMagnitude > 0 ){
        vec2 blend_range;

        vec3 abs_data = abs( source );
        float max_data = max(source.x, max(source.y, source.z));
        float magnitude = length(source);
        float value = max_data;

        float scaled_value = pow(smoothstep( range.x, range.y, value), pow_factor );;
        float start = 0.0;
        float end = 0.0; 
        float factor = 0.0;
        vec3 ramp_color = texture2D( dataColorRamp, vec2( dataOverlaySelector, scaled_value ) ).rgb;
        gl_FragColor.xyz = mix(ramp_color, gl_FragColor.xyz, dataOverlayPercent);
    }   
    else{
        vec3 smoothed_data = smoothstep( range.x, range.y, source.xyz );
        gl_FragColor.xyz = mix(smoothed_data, gl_FragColor.xyz, dataOverlayPercent);
    }


}
