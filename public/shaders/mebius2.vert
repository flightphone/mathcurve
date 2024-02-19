attribute vec2 size;
attribute vec2 size2;

uniform float u_time;
uniform float u_nx;
uniform float u_ny;
uniform vec3 camera;
uniform float u_shape;
uniform vec3 u_light;

varying vec3 norm;
varying vec3 rd;
varying vec2 uUv;
varying float sides;
varying vec3 light2;

#define PI  3.14159265359
#define TAU 6.28318530718

#define rot(a) mat2(cos(a), -sin(a), sin(a), cos(a))

//========================================surface==============================
vec3 rect(float t)
{
    float w  = 5., h = 1.5;
    float fi = PI/2., r0 = 0.25;
    vec3 d = vec3(r0*cos(t), r0*sin(t), 0);
    t = mod(t,TAU);
    
    float k = floor(t/fi);
    vec3 shift = vec3(0.);
    if (k == 0.)
        shift = vec3(w/2., h/2., 0.);
    if (k == 1.)    
        shift = vec3(-w/2., h/2., 0.);
    if (k == 2.)    
        shift = vec3(-w/2., -h/2., 0.);
    if (k == 3.)    
        shift = vec3(w/2., -h/2., 0.);    
    d += shift;
    return d;   
}

vec3 elli(float u)
{
    float a = 18., b = 10., r = 0., n = 1.;
    
    return vec3((a + r * cos(n * u))*cos(u), (b + r * cos(n * u))*sin(u), r * sin(n * u));
}

vec3 curve(float t, int n)
{
    if (n == 0)
        return elli(t);
    
}
mat3 getAxis(vec3 a) {
    a = normalize(a);
    vec3 z = a, x = vec3(1, 0, 0);
    if(a.x != 0.0)
        x = normalize(vec3(a.y, -a.x, 0));
    vec3 y = normalize(cross(x, a));
    return mat3(x, y, z);
}

struct HIT {
    vec3 val;
    mat3 xyz;
};



//return normal from curve by index
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

//return normal from curve by index
HIT curve_norm2(float t, int n)
{
    float h = 0.05;
    vec3 vt = curve(t, n);
    vec3 vth = curve(t+h, n);
    vec3 r1 = normalize(vth - vt);
    mat3 mt = getAxis(r1);
    return HIT(vt, mt);
}

vec3 mebius(float u, float v)
{
    vec3 rc = rect(v);
    rc.yz *= rot(PI/2.);
    rc.xz *= rot(.5*u);
    rc.xy *=rot(u);
    rc += elli(u);
    
    return rc;
    
}

vec3 sphere (float u, float v)
{
    float a = 2.;
    vec3 res = vec3(a*sin(v)*cos(u), a*sin(v)*sin(u), a*cos(v));
    return res;
}

//return surface by index
vec3 surf(float u, float v, int n)
{
    if (n==0)
        return mebius(u, v);
    if (n==1)    
        return sphere(u, v);
}

//return normal from surface by index
vec3 surf_normal(float u, float v, int n)
{
    float h = 0.01;
    vec3 du = surf(u+h, v, n) - surf(u-h, v, n);
    vec3 dv = surf(u, v+h, n) - surf(u, v-h, n);
    return  normalize(cross(du, dv));
}

vec3 surf_map(int n, float nu, float nv)
{
    float u = position.x*TAU*nu, v = position.y*TAU*nv;
    vec3 val = surf(u, v, n);
    norm = surf_normal(u, v, n);
    return val;
}






void main() {
    sides = size2.x;
    vec3 val = vec3(0.0);
    
    if (size2.x == 0.)
        val = surf_map(0, 1., 1.);
    else
    {
        
        float u = position.x*TAU, v = (1.-position.y)*PI, t = u_time; 
        val = sphere(u, v);
        
        
        float a = length(val);
        vec3 h = surf_normal(t, PI/2., 0)*(a+2.),
        m = elli(t);

        norm = normalize(val);
        val = val + m + h;
        light2 = normalize(m);
        
    }
    

    
    
    //returned
    norm = (modelViewMatrix * vec4(norm, 0)).xyz;
    vec4 worldPosition = modelViewMatrix * vec4(val, 1);
    rd = normalize(worldPosition.xyz - camera);
    gl_Position = projectionMatrix * worldPosition;
}