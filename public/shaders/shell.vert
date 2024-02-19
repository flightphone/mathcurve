uniform float u_time;
uniform vec3 camera;

varying vec3 norm;
varying vec3 rd;
varying vec2 uUv;

#define PI  3.14159265359
#define TAU 6.28318530718
#define rot(a) mat2(cos(a), -sin(a), sin(a), cos(a))

vec3 shell (float u, float v)
{
    //u = -u;
    float a = 3., b = 2.5, m = -0.1, k = 2.5;
    float x = exp(m*u)*cos(u)*(a + b*cos(v));
    float y = exp(m*u)*sin(u)*(a + b*cos(v));
    float z = exp(m*u)*(k*a + b*sin(v)), h = k*a+b;
    z -= h/2.;
    //z = -z;
    vec3 val = vec3(x, y, z);
    val.xz *= rot(PI/2.);
    return val;
}

vec3 sine (float u, float v)
{
    float a = 2.;
    return vec3(a*sin(u), a*sin(v), a*sin(u+v));
}

vec3 sphere (float u, float v)
{
    float a = 2.;
    return vec3(a*sin(v)*cos(u), a*sin(v)*sin(u), a*cos(v));
}


vec3 helicoid (float u, float v)
{
    float a = 2., h = 1.;
    return vec3(a*cos(v) - u*sin(v), a*sin(v)+ u*cos(v), h*v);
}

vec3 gorn(float u, float v)
{
    float h = 1./0.2;
    vec3 val = vec3(u*cos(v), u*sin(v), 1./u - h/2.);
    val.xz *= rot(PI/2.);
    return val;
}

float glz()
{
    float t = u_time/5.;
    float st = mod(floor(t), 4.);
    float res;
    if (st == 0.)
        res = 1.;
    if (st == 1.)    
        res = cos(fract(t)*PI/2.);//(1.- fract(t))*(1.- fract(t));
    if (st == 2.)    
        res = 0.;
    if (st == 3.)
        res = sin(fract(t)*PI/2.); //fract(t)*fract(t);   
    return res;    
}
vec3 surf (float u, float v)
{
    
    float pst = glz();
    vec3 v1 =  shell(u, v);
    vec3 v2 =  gorn(u/12.+0.2, v);
    vec3 res = mix(v1, v2, pst);
    return res;
    
    //return sine(u, v);
    //return helicoid (v, u);
    //return gorn(u, v);
    //return shell(u, v);
    //return sphere(u, PI - v);
}

vec3 surf_normal(float u, float v)
{
    float h = 0.05;
    vec3 du = surf(u+h, v) - surf(u-h, v);
    vec3 dv = surf(u, v+h) - surf(u, v-h);
    return  normalize(cross(du, dv));
}
varying float sides;
void main() {
    sides = 1.;
    uUv = uv;
    vec3 p = position;
    float u = p.x*TAU*7., v = p.y*TAU;
    //float u = uv.x*TAU*7., v = uv.y*TAU;
    //float u = uv.x*TAU + 0.2, v = uv.y*TAU;
    //float u = uv.x*TAU, v = uv.y*PI;
    
    //surface mesh
    vec3 val = surf(u, v);
    norm =  surf_normal(u, v);
    
    
    
    //returned
    norm = (modelViewMatrix * vec4(norm, 0)).xyz;
    vec4 wPosition = modelViewMatrix * vec4(val, 1);
    rd = normalize(wPosition.xyz - camera);
    gl_Position = projectionMatrix * wPosition;
}