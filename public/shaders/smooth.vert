attribute vec2 size;
attribute vec2 size2;

uniform float u_time;
uniform float u_nx;
uniform float u_ny;
uniform vec3 camera;
uniform float u_shape;

varying vec3 norm;
varying vec3 rd;
varying vec2 uUv;
varying float sides;

#define PI  3.14159265359
#define TAU 6.28318530718

#define rot(a) mat2(cos(a), -sin(a), sin(a), cos(a))

float heir_len = 4.;


struct HIT
{
    vec3 val;
    mat3 xyz;
};


vec3 sphere (float u, float v, float a, vec3 c)
{
    
    //v = PI - v;
    norm = vec3(sin(v)*cos(u), sin(v)*sin(u), cos(v));
    return c + a*norm;
}





vec3 curve(float t, int n)
{
    return vec3(0.);
}

//curve mesh
HIT curve_norm(float t, int n)
{
    float h = 0.05;
    vec3 vt = curve(t, n);
    vec3 vth = curve(t+h, n);
    vec3 r1 = normalize(vth - vt);
    vec3 r2 = normalize(curve(t+2.*h, n) - vth*2. + vt);
    
    vec3 x = normalize(cross(r1, r2));
    vec3 y = normalize(cross(x, r1));
    return HIT(vt, mat3(x, y, r1));
}

mat3 getAxis(vec3 a)
{
    a = normalize(a); 
    vec3 z = a, x  = vec3(1, 0, 0);
    if (a.x != 0.0) 
        x = normalize(vec3(a.y, -a.x, 0));
    vec3 y = normalize(cross(x, a));   
    return mat3(x, y, z);
}  

HIT curve_norm2(float t, int n)
{
    float h = 0.05;
    vec3 vt = curve(t, n);
    vec3 vth = curve(t+h, n);
    vec3 r1 = normalize(vth - vt);
    mat3 mt = getAxis(r1);
    return HIT(vt, mt);
}

vec2 P0, P1, P2, P3;
vec2 bezier(float t)
{
    float u = 1.-t;
    return u*u*u*P0 + 3.*t*u*u*P1 + 3.*t*t*u*P2 + t*t*t*P3;
}

vec3 bezier_cyrc(float t, float u)
{
    vec2 b = bezier(t);
    float z = b.x, x = b.y*cos(u), y = b.y*sin(u);
    return vec3(x, y, z);
}


vec3 surf(float u, float v, int n)
{
    if (n == 0)
        return bezier_cyrc(u, v);
    return vec3(0.);
    
}

vec3 surf_normal(float u, float v, int n)
{
    float h = 0.01;
    vec3 du = surf(u+h, v, n) - surf(u-h, v, n);
    vec3 dv = surf(u, v+h, n) - surf(u, v-h, n);
    return  normalize(cross(du, dv));
}

vec3 sphere_map()
{
    vec3 val = vec3(0.);
    //float a = 2.0,
    float h = position.x*2.,
    u = position.y*TAU, v = 0.;
    if (h >= 1.)
    {
        v = (h - 1.)*PI;
        val = sphere(u, v, 1.5, vec3(0., 0., -3.));
    } 
    else
    {
        v = h*PI;
        val = sphere(u, v, 2., vec3(0., 0., 1.5));
        
    }
    float h1 = 0.8, h2 = 1.6;
    if (h >= h1 && h <= h2)
    {
        v = h1*PI;
        val = sphere(0., v, 2., vec3(0., 0., 1.5));
        P0 = vec2(val.z, val.x);
        v += PI/2.;
        P1 = P0 + vec2(cos(v), sin(v));

        v = (h2-1.)*PI;
        val = sphere(0., v, 1.5, vec3(0., 0., -3.));
        P3 = vec2(val.z, val.x);
        v -= PI/2.;
        P2 = P3 + vec2(cos(v), sin(v));
        float t = (h-h1)/(h2 - h1);

        val = surf(t, u, 0);
        norm = surf_normal(t, u, 0);
    }        
    
    
    return val;
}


vec3 curve_map(int n, float k)
{
    
    float t = position.x*TAU*k, f = position.y*TAU, r = 0.1;
    vec3 val = curve(t, n);
    HIT res = curve_norm(t, n);
    norm = res.xyz*vec3(cos(f) , sin(f), 0.);
    val += norm*r;
    return val;
}

vec3 surf_map(int n, float nu, float nv)
{
    float u = position.x*TAU*nu, v = position.y*TAU*nv;
    vec3 val = surf(u, v, n);
    norm = surf_normal(u, v, n);
    return val;
}




void main() {
    sides = 1.;
    vec3 val = vec3(0.);
    val = sphere_map();
    
    
    norm = (modelViewMatrix * vec4(norm, 0)).xyz;
    vec4 worldPosition = modelViewMatrix * vec4(val, 1);
    rd = normalize(worldPosition.xyz - camera);
    gl_Position = projectionMatrix * worldPosition;
}