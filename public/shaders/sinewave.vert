precision highp float;
uniform float u_time;
uniform vec3 camera;
uniform float u_shape;

varying vec3 norm;
varying vec3 rd;
varying vec2 uUv;

#define PI  3.14159265359
#define TAU 6.28318530718

float numcurve = 0.;

vec3 tennis(float t) {
    float a = 3., b = 1., c = 2. * sqrt(a * b), x = a * cos(t) + b * cos(3. * t), y = a * sin(t) - b * sin(3. * t), z = c * sin(2. * t);
    return vec3(x, y, z);
}

vec3 sinewave(float t) {
    float a = 5.;
    float b = 2.5;
    float m = 3.5;
    float n = 9.;
    t = t;
    return vec3(a * cos(t), a * sin(t), b * sin(n / m * t));
}

vec3 trefoil(float t) {
    float a = 2.;
    return a * vec3(sin(t) + 2. * sin(2. * t), cos(t) - 2. * cos(2. * t), -0.5 * sin(3. * t));
}

vec3 circ(float t) {
    return vec3(2. * cos(t), 2. * sin(t), 0.);
}

vec3 eight_knot(float t) {
    return vec3((2. + cos(2. * t)) * cos(3. * t), (2. + cos(2. * t)) * sin(3. * t), 0.5 * sin(4. * t));
}

vec3 liss(float t) {
    float a = 3., b = 3., c = 2.5, n = 1.2, m = 2., f = PI / 2., e = .0;
    return vec3(a * sin(t), b * sin(n * t + f), c * sin(m * t + e));
}

vec3 solenoid(float t) {
    float n = 1.25, R = 6., r = 1.;
    return vec3((R + r * cos(n * t)) * cos(t), (R + r * cos(n * t)) * sin(t), r * sin(n * t));
}

vec3 rose(float t)
{
    float a = 3.,
    n = 2.2,
    b = 0.,
    r = a*cos(n*t);
    return vec3(r*cos(t), r*sin(t), b*cos(n*t)*cos(n*t));
}

float glz() {
    float t = u_time / 5.;
    float st = mod(floor(t), 4.);
    float res;
    if(st == 0.)
        res = 1.;
    if(st == 1.)
        res = cos(fract(t) * PI / 2.);//(1.- fract(t))*(1.- fract(t));
    if(st == 2.)
        res = 0.;
    if(st == 3.)
        res = sin(fract(t) * PI / 2.); //fract(t)*fract(t);   
    return res;
}
vec3 curve(float t) {
    if(numcurve == 1.)
        return solenoid(t);
    if(numcurve == 2.)
        return liss(t);
    if(numcurve == 3.)
        return sinewave(t);

    if(numcurve == 4.)
        return trefoil(t);    
    if(numcurve == 5.)
        return tennis(t);     

    if(numcurve == 6.)
        return eight_knot(t);           
    
    if(numcurve == 7.)
        return rose(t);
    if(numcurve == 0.) {
        float pst = glz();
        vec3 v1 = liss(t * 5. / 7.);//solenoid(t*4./7.);circ(t);trefoil(t);
        vec3 v2 = sinewave(t);//solenoid(t*4./7.);//eight_knot(t);
        vec3 v = mix(v1, v2, pst);
        return v;
    }

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

HIT curve_norm(float t) {
    float h = 0.1;
    vec3 vt = curve(t);
    vec3 vth = curve(t + h);
    vec3 r1 = normalize(vth - vt);
    //h = 0.2;
    vec3 r2 = normalize(curve(t + 2. * h) - vth * 2. + vt);

    vec3 x = normalize(cross(r1, r2));
    vec3 y = normalize(cross(x, r1));
    return HIT(vt, mat3(x, y, r1));
}
varying float sides;

void main() {
    sides = 1.;
    uUv = uv;
    vec3 p = position;
    float t = p.x * TAU, f = p.y * TAU;
    float r = .2;
    numcurve = u_shape;
    if(numcurve == 1.)
    {
        r = 0.4;
        t*= 4.;
    }
    if(numcurve == 0.)
        t *= 7.;
    if(numcurve == 2.)
        t*=5.;
    if(numcurve == 3.)
    {
        r = .3;
        t*=7.;
    }

    if(numcurve == 4.)
    {
        r = .4;
        t*=2.;  
    }
    if(numcurve == 5.)
    {
        r = .4;
        t*=2.;
    }

    if(numcurve == 6.)
    {
        r = .3;
        t*=2.;
    }
    if(numcurve == 7.)
    {
        r = .2;
        t*=5.;
    }

    HIT mn = curve_norm(t);
    mat3 mt = mn.xyz;
    vec3 vt = mn.val;

    
    norm = mt * vec3(cos(f), sin(f), 0.);
    vec3 val = vt + norm * r;

    //returned
    norm = (modelViewMatrix * vec4(norm, 0)).xyz;
    vec4 worldPosition = modelViewMatrix * vec4(val, 1);
    rd = normalize(worldPosition.xyz - camera);
    gl_Position = projectionMatrix * worldPosition;
}