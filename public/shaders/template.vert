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

vec3 heir(float t)
{
    float fi = u_time*4., z = t, h = heir_len, 
    x = sin(z * TAU  + fi + 0.8) * h * 0.3 * z * (h - z) / h / h, 
    y = sin(z * TAU  - fi + 0.8) * h * 0.3 * z * (h - z) / h / h;
    return vec3(x, y, z);
}  

vec3 larme(float u, float v)
{
    float a = 2.,
    zoom = 1.5,
    l = a*sin(v)*sin(v/2.)*sin(v/2.)*zoom,
    z = a*cos(v),
    x = l*cos(u),
    y = l*sin(u);
    return vec3(x, y, z);
}

vec3 sphere (float u, float v)
{
    float a = 2.;
    norm = vec3(sin(v)*cos(u), sin(v)*sin(u), cos(v));
    return a*norm;
}

vec3 tube(float t)
{
    float n = 3., fi = TAU/n, w = 4., r0 = 0.25;
    vec3 d = vec3(r0*cos(t), r0*sin(t), 0);
    t = mod(t,TAU);
    float k = floor(t/fi)*fi + fi/2.;
    vec3 shift = vec3(w*cos(k), w*sin(k), 0);
    d += shift;
    return d;   
}

vec3 hexogram(float t)
{
    float r0 = .25, w = 3., sg = 1., angl = -PI/2.;
    vec3 c = vec3(1.5*w+r0+0.2, -w*sqrt(3.)/2., 0.);
    //t = mod(t, (4.*TAU));
    for (int i = 0; i < 13; i++)
    {
        float f = (sg == 1.)? PI*2./3. : PI/3.;
        if (t <= f)
        {
            vec3 d = vec3(r0*cos(t*sg + angl), r0*sin(t*sg + angl), 0) + c;    
            return d;
        }
        else
        {
            angl += (f*sg);
            c += 2.* r0* vec3(cos(angl), sin(angl), 0.);
            c += w*vec3(cos(angl + sg*PI/2.), sin(angl+ sg*PI/2.), 0.);
            angl -= PI;
            t -= f;
            sg *= -1.;
            
        }
    }
    return vec3(-4.);
}

vec3 rose(float t)
{
    float a = 3., n = 2.2, r = a*cos(n*t);
    return vec3(r*cos(t), r*sin(t), 0.0);
}    

vec3 curve(float t, int n)
{
    if (n == 0)
        return heir(t);
    if (n == 1)
        return tube(t);  
    if (n == 2) 
        return rose(t);     
    if (n == 3) 
        return hexogram(t);         
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



vec3 surf(float u, float v, int n)
{
    if (n == 0)
        return larme(u, v);
    if (n == 1)
        return sphere(u, v);
}

vec3 surf_normal(float u, float v, int n)
{
    float h = 0.01;
    vec3 du = surf(u+h, v, n) - surf(u-h, v, n);
    vec3 dv = surf(u, v+h, n) - surf(u, v-h, n);
    return  normalize(cross(du, dv));
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
vec3 curve_map(int n, float k)
{
    
    float t = position.x*TAU*k, f = position.y*TAU, r = 0.1;
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


vec3 heir_map() {
    float h = heir_len;
    float t = -1. + position.x*(h+2.), f = position.y*TAU;
    float r1 = 0.2, r0 = 0.1;
    vec3 val = vec3(0.);
    if (t <= h && t >= 0.)
    {
        HIT res =  curve_norm(t, 0);   
        float r = r1 * (h - t) / h + r0;
        norm = res.xyz*vec3(cos(f) , sin(f), 0.);
        val = res.val + norm*r;
    }    
    else
    {
        float t0 = clamp(t, 0., h),
        fi = clamp((t-t0)*PI/2., -PI/2., PI/2.);
        HIT res =  curve_norm(t0, 0);   
        float r = r1 * (h - t0) / h + r0;
        //vec3 x = res.x, y = res.y, z = res.z;
        norm = res.xyz*vec3(cos(f)*cos(fi), sin(f)*cos(fi), sin(fi));
        val = res.val + norm*r;
    }
    val.z -= h/2.;
    return val;
}

vec3 larme_map()
{
    float h = PI/3.;
    float u = (1. - position.x)*TAU, v = position.y*0.9*PI;
    float ub = (1. - size.x)*TAU, vb = size.y*0.9*PI;
    vec3 val = vec3(0.);
    if (v > h)
    {
        val = larme(u, v);
        norm = surf_normal(u, v, 0);
    }    
    else
    {
        float a = 2., zoom = 1.5, l = a*sin(h)*sin(h/2.)*sin(h/2.)*zoom,
        z = a*cos(h);
        l = l*v/h;
        float x = l*cos(u),
        y = l*sin(u);
        val = vec3(x, y, z);
        norm = vec3(0, 0, -1.);
    }

    return val;
    //block
    
    vec3 valb = vec3(0.);
    vec3 normb = vec3(0.);

    if (vb > h)
    {
        valb = larme(ub, vb);
        normb = surf_normal(ub, vb, 0);
    }    
    else
    {
        float a = 2., zoom = 1.5, l = a*sin(h)*sin(h/2.)*sin(h/2.)*zoom,
        z = a*cos(h);
        l = l*vb/h;
        float x = l*cos(ub),
        y = l*sin(ub);
        valb = vec3(x, y, z);
        normb = vec3(0, 0, -1.);
    }
    
    float pst = 1. - glz();
    val -= valb;
    val.xz *= rot(pst*10.*TAU);
    //val.yz *= rot(pst*3.*TAU);
    norm.xz *= rot(pst*10.*TAU);
    //norm.yz *= rot(pst*3.*TAU);
    val += valb;
    normb = normalize(valb);
    val += normb*pst*50.;
    
    return val;

}
vec3 hexogram2()
{
    vec3 val = cube_map2(1, 1.);
    if (size2.x == 0.)
        val.xy*=rot(PI/6.);
    else    
        val.xy*=rot(-PI/6.);
    return val;
}

vec3 bio()
{
    vec3 val = vec3(0.);
    if (size2.x == 0.)
    {
        float u = position.x*TAU, v=(1.-position.y)*PI;
        val = sphere(u, v);
    }
    else   
    { 
        val = heir_map();
        val.z += (heir_len/2. + 2.);
    }
    return val;    
}
void main() {
    sides = 1.;
    vec3 val = vec3(0.);
    if (u_shape == 0.)
        val = larme_map();
    if (u_shape == 1.)
        val = heir_map();   
    if (u_shape == 2.)
        val = curve_map(2, 5.);
    if (u_shape == 3.)
        val = cube_map(1, 1.);
    if (u_shape == 4.)
        val = cube_map2(2, 5.);
    if (u_shape == 5.)
        val = cube_map2(1, 1.);
    if (u_shape == 6.)
        val = cube_map2(3, 4.);

    //vec3 val = curve_map(2, 5.);
    //vec3 val = cube_map(1, 1.);
    //vec3 val = cube_map2(2, 5.);
    //vec3 val = cube_map2(1, 1.);
    //vec3 val = cube_map2(3, 4.);
    //vec3 val = hexogram2();
    //vec3 val = bio();
    //vec3 val = larme_map();
    //vec3 val = heir_map();

    //vec3 val = cube_map2(3, 4.);
    
    //returned
    norm = (modelViewMatrix * vec4(norm, 0)).xyz;
    vec4 worldPosition = modelViewMatrix * vec4(val, 1);
    rd = normalize(worldPosition.xyz - camera);
    gl_Position = projectionMatrix * worldPosition;
}