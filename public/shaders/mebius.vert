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


struct HIT
{
    vec3 val;
    mat3 xyz;
};
//========================================surface==============================
vec3 mebius(float v, float u)
{
    //https://mathcurve.com/surfaces.gb/mobiussurface/mobiussurface.shtml
    float a = 1.5, 
    x = (a + u* cos(v/2.))*cos(v),
    y = (a + u* cos(v/2.))*sin(v),
    z = u*sin(v/2.);
    return vec3(x, y, z);
}
//https://mathcurve.com/surfaces.gb/cycliddedupin/cyclidededupin.shtml
vec3 dedupin(float u, float v)  
{
    float c = .5, b = 1.1, d = .65, a = 1.3,
    dd = (a - c*cos(u)*cos(v)),
    x = (d*(c - a*cos(u)*cos(v)) + b*b*cos(u))/dd,
    y = b*sin(u)*(a - d*cos(v))/dd,
    z = b*sin(v)*(c*cos(u) - d)/dd;
    return vec3(x, y, z);
}


vec3 klein(float u, float v)
{
    sides = 0.;
    float a = 3., b = 4., c = 2., r, x, y, z;
    if (u < PI)
    {
        r = c*(1. - cos(u)/2.);
        x = cos(u)*(a*(1.+sin(u)) + r*cos(v));
        y = sin(u)*(b + r*cos(v));
        z = r*sin(v);
    }
    
    else
    {
        r = c*(1. - cos(u)/2.);
        x = cos(u)*a*(1.+sin(u)) - r*cos(v);
        y = sin(u)*b;
        z = r*sin(v);
    }
    
    vec3 res =  vec3(x, y, z);
    res.yz*=rot(PI);
    return res;

}

vec3 boys(float u, float v)
{
    sides = 0.;
    float k = 1.,
    ka = cos(u)/(sqrt(2.)- k*sin(2.*u)*sin(3.*v)),
    z = ka*3.0*cos(u),
    x = ka*(cos(u)*cos(2.*v) + sqrt(2.)*sin(u)*cos(v)),
    y = ka*(cos(u)*sin(2.*v) - sqrt(2.)*sin(u)*sin(v));
    vec3 res =  vec3(x, y, z);
    res.z -= sqrt(2.);
    return res;
}

vec3 sine (float u, float v)
{
    sides = 0.;
    float a = 2.;
    return vec3(a*sin(u), a*sin(v), a*sin(u+v));
}


vec3 quart(float t)
{
    float r = 2., n  = 4., dn = TAU/n, c = t/dn, a0 = floor(c)*dn, 
    a1 = a0 + dn, d = fract(c);
    vec3 v1 = vec3(cos(a0), sin(a0), 0.);
    vec3 v2 = vec3(cos(a1), sin(a1), 0.);
    vec3 res = r*(v1 + (v2-v1)*d);
    return res;

}

vec3 twist(float u, float v)
{
        vec3 res = quart(TAU - v);
        res.xy *= rot(u);
        float a = 1.;
        res.z = a*u;
        res.z -= TAU;
        res.xz *= rot(PI/2.);
        return res;
}

vec3 Breather(float v, float u)
{
    //sides = 0.;
    u -= 1.*TAU;
    //v-= 12.*TAU;
    float a = 2.5/5., b = 1.-a*a,
    d = a*(b*cosh(a*u)*cosh(a*u) + a*a*sin(sqrt(b)*v)*sin(sqrt(b)*v)),
    x = -u + 2.*b*cosh(a*u)*sinh(a*u)/d,
    y = 2.*b*cosh(a*u)*(-sqrt(b)*cos(v)*cos(sqrt(b)*v) - sin(v)*sin(sqrt(b)*v))/d,
    z = 2.*b*cosh(a*u)*(-sqrt(b)*sin(v)*cos(sqrt(b)*v) + cos(v)*sin(sqrt(b)*v))/d;
    return vec3(x, y, z);
}

vec3 shell (float u, float v)
{
    float a = 3., b = 2.5, m = -0.1, k = 2.5;
    float x = exp(m*u)*cos(u)*(a + b*cos(v));
    float y = exp(m*u)*sin(u)*(a + b*cos(v));
    float z = exp(m*u)*(k*a + b*sin(v)), h = k*a+b;
    z -= h/2.;
    vec3 val = vec3(x, y, z);
    val.xy *= rot(PI/8.);
    val.xz *= rot(-PI/2.5);
    return val;
}
//https://mathcurve.com/surfaces.gb/kuen/kuen.shtml
float ukuen = 1.4305;
vec3 kuen(float u, float v)
{
    //u -= PI;
    u -= ukuen*TAU/2.;
    float a = 5.,
    d = 1. + u*u*sin(v)*sin(v),
    x = a*(cos(u) + u*sin(u))*sin(v)/d,
    y = a*(sin(u)- u*cos(u))*sin(v)/d,
    z = a*(0.5*log(tan(v/2.)) + cos(v)/d);
    return vec3(x, y, z);
}


//return surface by index
vec3 surf(float u, float v, int n)
{
    if (n == 0)
        return mebius(u, v);
    if (n == 1) 
        return dedupin(u, v);  
    if (n == 2) 
        return klein(u, v); 
    if (n == 3)    
        return boys(u, v);
    if (n == 4)    
        return sine(u, v);
    if (n == 5)
       return twist(u, v);
    if (n==6)        
        return Breather(u, v);
    if (n==7)        
        return shell(u, v);    
    if (n==8)        
        return kuen(u, v);        
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




//========================================curves==============================

vec3 liss(float t)
{
    float a = 3.,
    b = 3.,
    c = 1.5,
    n = 1.5,
    m = 2.5,
    f = PI/2.,
    e = .0;
    return vec3(a*sin(t), b*sin(n*t + f), c*sin(m*t + e));
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

//return curve by index
vec3 curve(float t, int n)
{
    if (n == 0)
        return quart(t);
    if (n == 1)
        return liss(t);    
}

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


vec3 curve_map(int n, float k)
{
    
    float t = position.x*TAU*k, f = position.y*TAU, r = 0.3;
    vec3 val = curve(t, n);
    HIT res = curve_norm(t, n);
    norm = res.xyz*vec3(cos(f) , sin(f), 0.);
    val += norm*r;
    return val;
}

vec3 cube_map(int n, float k)
{
    float hh = 1.5, r = 0.2;
    float t = position.x*TAU*k, y = position.y, f = 0.;
    vec3 shift = vec3(0.);
    HIT res =  curve_norm2(t, n);
    vec3 val = res.val;
    float eps = 0.05;
    if (y < 0.5)
    {
        f = y*PI + 1.5*PI;
        norm = res.xyz*vec3(cos(f) , sin(f), 0.);
        val += (norm*r + shift);
    }
    if (y < eps)
    {
        val.z = -r;
        norm = vec3(0., 0., -1.);
    }
    if (y == 0.)
    {
        val = vec3(0., 0., -r);
    }
    if (y >= 0.5)
    {
        f = (y-0.5)*PI;
        norm = res.xyz*vec3(cos(f) , sin(f), 0.);
        shift.z = hh;
        val += (norm*r + shift);
    }
    
    if (y > 1.-eps)
    {
        val.z = hh + r;
        norm = vec3(0., 0., 1.);
    }
    if (y == 1.)
    {
        val = vec3(0., 0., hh+r);
    }
    val.z -= (hh/2.);
    return val;
}



vec3 cube_map2(int n, float k)
{
    float hh = 1.5, r = 0.2;
    float t = position.x*TAU*k, y = position.y, f = 0.;
    vec3 shift = vec3(0.);
    HIT res =  curve_norm2(t, n);
    vec3 val = res.val;
    if (y < 0.5)
    {
        f = y*PI*2.0 + PI;
        norm = res.xyz*vec3(cos(f) , sin(f), 0.);
        val += (norm*r + shift);
    }
    
    if (y >= 0.5)
    {
        f = (y-0.5)*PI*2.0;
        norm = res.xyz*vec3(cos(f) , sin(f), 0.);
        shift.z = hh;
        val += (norm*r + shift);
    }
    if (y == 1.)
        val.z = 0.;
    
    
    val.z -= (hh/2.);
    return val;
}


void main() {
    sides = 1.;
    vec3 val = vec3(0.0);
    if (u_shape == 0.)
        //val = surf_map(8, ukuen, 1.);  //Breather  22.01.2024
        val = surf_map(6, 7.45, 2.5);  //Breather  22.01.2024
    if (u_shape == 1.)
         val = surf_map(0, 2., 1.); //mebius
    if (u_shape == 2.)
        val = surf_map(1, 1., 1.);  //dupin
    if (u_shape == 3.)
        val = surf_map(2, 1., 1.);  //klein
    if (u_shape == 4.)
        val = surf_map(3, .5,.5);  //boys 22.01.2024
    if (u_shape == 5.)
        val = surf_map(4, 1.,1.);  //sine 22.01.2024
    if (u_shape == 6.)
        val = surf_map(5, 2.,1.);  //twist 22.01.2024
    
    if (u_shape == 7.)
        val = surf_map(7, 7.,1.);  //shell 24.01.2024
    if (u_shape == 8.)    
        val = curve_map(1, 5.);
    if (u_shape == 9.)
        val = surf_map(8, ukuen, 1.);  //kuen  01.02.2024    
    
    
    //returned
    norm = (modelViewMatrix * vec4(norm, 0)).xyz;
    vec4 worldPosition = modelViewMatrix * vec4(val, 1);
    rd = normalize(worldPosition.xyz - camera);
    gl_Position = projectionMatrix * worldPosition;
}