varying vec3 vViewPosition;
varying vec3 vNormal;
varying vec2 vUv;
//varying vec3 vData;
//attribute vec3 data_overlay;

void main() {
    vUv = uv;
    vec3 transformedNormal = normalMatrix * normal;
    vNormal = normalize( transformedNormal ); // probably unnecessary
    vec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );
    //vData = data_overlay;
    gl_Position = projectionMatrix * mvPosition;
    vViewPosition = -mvPosition.xyz;
}
