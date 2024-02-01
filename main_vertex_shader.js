import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { GUI } from 'three/addons/libs/lil-gui.module.min.js';

let nx = 2400;  //widthSegments
let ny = 200;   //heightSegments
let nbx = 1;    //width distinct areas for one object
let nby = 1;    //height distinct areas for one object
let grid = false; //if grid = false then all areas represen xy coordinates from 0 to 1, width = 1, height = 1
//else (grid = true) area represent tail width = 1/nbx, heigth = 1/nby


//init scene
let renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(window.devicePixelRatio);
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);
const collections = ["surfaces", "curves", "transform", "others"]
const shapes = ["breather", "mebius", "dupin", "klein", "boys", "sine", "twist", "shell", "liss", "kuen"];
const curves = ["mix", "solenoid", "liss", "sinewave", "trefoil", "tennis", "eight_knot", "rose"];
const others = ["larme", "heir", "2", "3", "4", "5", "6"];
let manager = { 
    collection: "surfaces",
    surfaces: "breather", 
    curves: "mix",
    others: "larme"
}
const gui = new GUI();
gui.title("Math Collection");
gui.add(manager, "collection", collections);
gui.add(manager, "surfaces", shapes);
gui.add(manager, "curves", curves);
gui.add(manager, "others", others);


const camera = new THREE.PerspectiveCamera(75, 2, 0.1, 500);
const scene = new THREE.Scene();
//look
camera.position.set(0, 0, -10);
camera.up.set(0, 1, 0);
camera.lookAt(0, 0, 0);
//mouse rotate
let controls = new OrbitControls(camera, renderer.domElement);
controls.update();
window.addEventListener('resize', onWindowResize);
scene.background = new THREE.Color(0x888888);

const uniforms = {
    u_time: { value: 0.0 },
    camera: { value: new THREE.Vector3() },
    u_light: { value: new THREE.Vector3(0, 0, 1) },
    u_nx: { value: nx },
    u_ny: { value: ny },
    u_shape: { value: 0.0 }
};

let vertexShader = "";
let fragmentFront = `
uniform vec3 u_light;
varying vec3 norm;
varying vec3 rd;
varying vec2 uUv;
void main() {
    vec3 nor = norm;
    vec3 light = normalize(u_light);
    vec3 col = vec3(.5, 0.4,  1.0);
    vec3 R = reflect(light, nor);
    float specular  = pow(max(abs(dot(R, rd)), 0.), 10.);
    float difu = abs(dot(nor, light));
    col = col*(col*clamp(difu, 0., 1.0) + 0.5) + vec3(.5)*specular*specular;
    gl_FragColor = vec4(col , 1. );
}
`;
let fragmentBack = `
uniform vec3 u_light;
varying vec3 norm;
varying vec3 rd;
varying vec2 uUv;
varying float sides;
void main() {
    vec3 nor = norm;
    vec3 light = normalize(u_light);
    vec3 col = vec3(.2, 1.,  .3);
    if (sides == 0.)
        col = vec3(.5, 0.4,  1.0);
    vec3 R = reflect (light, nor);
    float specular  = pow(max(abs(dot(R, rd)), 0.), 10.);
    float difu = abs(dot(nor, light));
    col = col*(col*clamp(difu, 0., 1.0) + 0.5) + vec3(.5)*specular*specular;
    gl_FragColor = vec4(col , 1. );
}
`;

//init simple plane geometry
const plane = new THREE.Object3D();
scene.add(plane);
const geometry = new THREE.BufferGeometry();
const indices = [];
const vertices = [];
const sizes = [];
const sizes2 = [];
let nvert = 0;
let x0 = 0.;
let x1 = 1.;
let y0 = 0.;
let y1 = 1.;
let startvertex = 0;
for (let ib = 0; ib < nbx; ib++) {
    for (let jb = 0; jb < nby; jb++) {
        if (grid) {
            x0 = ib / nbx;
            x1 = x0 + 1. / nbx;
            y0 = jb / nby;
            y1 = y0 + 1. / nby;
        }
        //let vert = [];
        let blx = (x0 + x1) / 2.;
        let bly = (y0 + y1) / 2.;
        startvertex = nvert;
        for (let i = 0; i < nx; i++) {
            for (let j = 0; j < ny; j++) {
                let x = x0 + (x1 - x0) * i / (nx - 1);
                let y = y0 + (y1 - y0) * j / (ny - 1);
                vertices.push(x);
                vertices.push(y);
                vertices.push(0.);
                sizes.push(blx);
                sizes.push(bly);

                sizes2.push(ib);
                sizes2.push(jb);
                
                nvert += 1;
            }
        }
        for (let i = 0; i < nx - 1; i++) {
            for (let j = 0; j < ny - 1; j++) {
                let a = i * ny + j;
                let b = a + 1;
                let c = a + ny;
                let d = b + ny;
                
                a += startvertex;
                b += startvertex;
                c += startvertex;
                d += startvertex;

                indices.push(a);
                indices.push(c);
                indices.push(b);

                indices.push(b);
                indices.push(c);
                indices.push(d);
            }
        }
    }
}
geometry.setIndex(indices);
geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
geometry.setAttribute('size', new THREE.Float32BufferAttribute(sizes, 2).setUsage(THREE.DynamicDrawUsage));
geometry.setAttribute('size2', new THREE.Float32BufferAttribute(sizes2, 2).setUsage(THREE.DynamicDrawUsage));
const planeFront = new THREE.Mesh(geometry);
const planeBack = new THREE.Mesh(geometry);
plane.add(planeFront)
plane.add(planeBack);


//list shaders for different shapes and curves
const verts = ['shaders/mebius.vert', 'shaders/sinewave.vert', 'shaders/shell.vert', 'shaders/template.vert'];
const allshapes = [shapes, curves, [], others];
let vert = 'shaders/mebius.vert';

loadVertexShader(vert);
requestAnimationFrame(render);




function loadVertexShader(vertUrl)
{
    fetch(vertUrl).then(response => response.text())
    .then(txtShader => {
        vertexShader = txtShader;
        const objShaderFront = {
            uniforms: uniforms,
            side: THREE.FrontSide,
            vertexShader: vertexShader,
            fragmentShader: fragmentFront,
        }

        const objShaderBack = {
            uniforms: uniforms,
            side: THREE.BackSide,
            vertexShader: vertexShader,
            fragmentShader: fragmentBack,
        }
        //shader
        const materialFront = new THREE.ShaderMaterial(objShaderFront);
        const materialBack = new THREE.ShaderMaterial(objShaderBack);
        planeFront.material = materialFront;
        planeBack.material = materialBack;

       
    });
}


function onWindowResize() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
}

function render(time) {
    time *= 0.001; // convert to seconds;
    uniforms.u_time.value = time;
    uniforms.camera.value = camera.position;
    let clt = collections.indexOf(manager.collection);
    let vsh = verts[clt];
    if (vsh!=vert)
    {
        vert = vsh;
        loadVertexShader(vert);
    }
    let cshape = manager.surfaces;
    if (clt == 1)
        cshape = manager.curves;
    if (clt == 3)
        cshape = manager.others;

    uniforms.u_shape.value = allshapes[clt].indexOf(cshape);
    //plane.rotation.x += 0.01;
    //plane.rotation.y += 0.01;
    controls.update();
    renderer.render(scene, camera);
    requestAnimationFrame(render);
}


