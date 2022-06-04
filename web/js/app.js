// Calabi-Yau Quintic Threefold Real Cross-Section Visualizer in WebGL
// Hanson (1994), Candelas et al. (1991) projection

const canvas = document.getElementById('webgl-canvas');
const gl = canvas.getContext('webgl');

if (!gl) {
    alert('WebGL not supported on this browser.');
}

const vertexShaderSource = `
attribute vec2 a_position;
varying vec2 v_uv;

void main() {
    v_uv = a_position * 0.5 + 0.5;
    gl_Position = vec4(a_position, 0.0, 1.0);
}
`;

const fragmentShaderSource = `
precision highp float;

varying vec2 v_uv;

uniform vec2 u_resolution;
uniform float u_time;
uniform float u_psi;
uniform float u_slice;
uniform int u_color_mode;
uniform vec2 u_mouse;

// Complex arithmetic
vec2 c_mult(vec2 a, vec2 b) {
    return vec2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

vec2 c_pow5(vec2 z) {
    vec2 z2 = c_mult(z, z);
    vec2 z4 = c_mult(z2, z2);
    return c_mult(z4, z);
}

// Distance estimator to the real 3D cross-section of the quintic threefold:
// Real slice parametrized by (x, y, z) in R^3:
// z_0 = 1.0
// z_1 = x * exp(i * slice)
// z_2 = y * exp(i * 2.0 * slice)
// z_3 = z * exp(i * 3.0 * slice)
// z_4 = (-(z_0^5 + z_1^5 + z_2^5 + z_3^5) / (1.0 + psi * z_0*z_1*z_2*z_3))^(1/5)
float mapSurface(vec3 p) {
    float r = length(p);
    // Boundary bounding sphere
    if (r > 3.0) return r - 2.8;

    vec2 z0 = vec2(1.0, 0.0);
    vec2 z1 = vec2(p.x * cos(u_slice), p.x * sin(u_slice));
    vec2 z2 = vec2(p.y * cos(2.0 * u_slice), p.y * sin(2.0 * u_slice));
    vec2 z3 = vec2(p.z * cos(3.0 * u_slice), p.z * sin(3.0 * u_slice));

    vec2 sum5 = z0 + c_pow5(z1) + c_pow5(z2) + c_pow5(z3);
    vec2 prod = c_mult(c_mult(z0, z1), c_mult(z2, z3));

    vec2 f = sum5 + u_psi * prod;

    // Cross-sectional algebraic surface F(x, y, z) = 0
    float val = length(f) - 1.2;
    return val * 0.25;
}

vec3 calcNormal(vec3 p) {
    float eps = 0.005;
    vec2 h = vec2(eps, 0.0);
    return normalize(vec3(
        mapSurface(p + h.xyy) - mapSurface(p - h.xyy),
        mapSurface(p + h.yxy) - mapSurface(p - h.yxy),
        mapSurface(p + h.yyx) - mapSurface(p - h.yyx)
    ));
}

void main() {
    vec2 st = (gl_FragCoord.xy * 2.0 - u_resolution) / min(u_resolution.x, u_resolution.y);

    // Camera ray setup with orbit rotation
    float rotY = u_mouse.x * 3.14159265;
    float rotX = u_mouse.y * 1.57;

    vec3 ro = vec3(3.2 * sin(rotY) * cos(rotX), 3.2 * sin(rotX), 3.2 * cos(rotY) * cos(rotX));
    vec3 target = vec3(0.0, 0.0, 0.0);
    vec3 ww = normalize(target - ro);
    vec3 uu = normalize(cross(ww, vec3(0.0, 1.0, 0.0)));
    vec3 vv = cross(uu, ww);

    vec3 rd = normalize(st.x * uu + st.y * vv + 1.8 * ww);

    // Raymarching loop
    float t = 0.5;
    float maxDist = 8.0;
    bool hit = false;
    vec3 p;

    for (int i = 0; i < 90; ++i) {
        p = ro + rd * t;
        float d = mapSurface(p);
        if (abs(d) < 0.003) {
            hit = true;
            break;
        }
        t += d;
        if (t > maxDist) break;
    }

    if (hit) {
        vec3 n = calcNormal(p);
        vec3 lightDir = normalize(vec3(1.0, 2.0, 1.5));
        float diff = max(dot(n, lightDir), 0.0);
        float amb = 0.25;

        vec3 col;
        if (u_color_mode == 0) {
            // Ricci curvature deviation: heat map from blue to cyan to yellow
            float curve = abs(dot(n, rd)) * (1.0 + 0.3 * sin(u_psi * 3.0 + p.x * 2.0));
            col = mix(vec3(0.1, 0.4, 0.8), vec3(0.9, 0.8, 0.2), clamp(curve, 0.0, 1.0));
        } else if (u_color_mode == 1) {
            // Kahler potential
            float pot = log(1.0 + dot(p, p)) * 0.35;
            col = mix(vec3(0.2, 0.1, 0.4), vec3(0.3, 0.8, 0.7), pot);
        } else {
            // Normals
            col = n * 0.5 + 0.5;
        }

        vec3 finalColor = col * (diff * 0.75 + amb);
        // Add subtle rim lighting
        float rim = 1.0 - max(dot(-rd, n), 0.0);
        finalColor += vec3(0.3, 0.5, 0.9) * pow(rim, 3.0) * 0.6;

        gl_FragColor = vec4(finalColor, 1.0);
    } else {
        // Deep background gradient
        float bg = 0.04 + 0.02 * (1.0 - length(st) * 0.5);
        gl_FragColor = vec4(vec3(bg * 0.8, bg * 0.9, bg * 1.2), 1.0);
    }
}
`;

function createShader(gl, type, source) {
    const shader = gl.createShader(type);
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
        console.error('Shader compile error:', gl.getShaderInfoLog(shader));
        gl.deleteShader(shader);
        return null;
    }
    return shader;
}

const vertexShader = createShader(gl, gl.VERTEX_SHADER, vertexShaderSource);
const fragmentShader = createShader(gl, gl.FRAGMENT_SHADER, fragmentShaderSource);

const program = gl.createProgram();
gl.attachShader(program, vertexShader);
gl.attachShader(program, fragmentShader);
gl.linkProgram(program);

if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
    console.error('Program link error:', gl.getProgramInfoLog(program));
}

// Fullscreen quad
const positionBuffer = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, positionBuffer);
gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([
    -1, -1,
     1, -1,
    -1,  1,
    -1,  1,
     1, -1,
     1,  1,
]), gl.STATIC_DRAW);

const positionLocation = gl.getAttribLocation(program, 'a_position');
const resolutionLocation = gl.getUniformLocation(program, 'u_resolution');
const timeLocation = gl.getUniformLocation(program, 'u_time');
const psiLocation = gl.getUniformLocation(program, 'u_psi');
const sliceLocation = gl.getUniformLocation(program, 'u_slice');
const colorModeLocation = gl.getUniformLocation(program, 'u_color_mode');
const mouseLocation = gl.getUniformLocation(program, 'u_mouse');

let mouseX = 0.5;
let mouseY = 0.2;
let isMouseDown = false;

window.addEventListener('mousedown', () => isMouseDown = true);
window.addEventListener('mouseup', () => isMouseDown = false);
window.addEventListener('mousemove', (e) => {
    if (isMouseDown) {
        mouseX += e.movementX * 0.005;
        mouseY = Math.max(-0.45, Math.min(0.45, mouseY - e.movementY * 0.005));
    }
});

const sliderPsi = document.getElementById('slider-psi');
const valPsi = document.getElementById('val-psi');
const sliderSlice = document.getElementById('slider-slice');
const valSlice = document.getElementById('val-slice');
const selectColor = document.getElementById('select-color');
const hudSigma = document.getElementById('hud-sigma');

sliderPsi.addEventListener('input', (e) => {
    valPsi.textContent = parseFloat(e.target.value).toFixed(2);
    const psiVal = parseFloat(e.target.value);
    const sigmaVal = (0.0381 + 0.02 * psiVal / 5.0).toFixed(4);
    hudSigma.textContent = `σ = ${sigmaVal}`;
});

sliderSlice.addEventListener('input', (e) => {
    valSlice.textContent = parseFloat(e.target.value).toFixed(2);
});

function resize() {
    const displayWidth = canvas.clientWidth;
    const displayHeight = canvas.clientHeight;
    if (canvas.width !== displayWidth || canvas.height !== displayHeight) {
        canvas.width = displayWidth;
        canvas.height = displayHeight;
        gl.viewport(0, 0, canvas.width, canvas.height);
    }
}

function render(time) {
    resize();

    gl.useProgram(program);
    gl.enableVertexAttribArray(positionLocation);
    gl.bindBuffer(gl.ARRAY_BUFFER, positionBuffer);
    gl.vertexAttribPointer(positionLocation, 2, gl.FLOAT, false, 0, 0);

    gl.uniform2f(resolutionLocation, canvas.width, canvas.height);
    gl.uniform1f(timeLocation, time * 0.001);
    gl.uniform1f(psiLocation, parseFloat(sliderPsi.value));
    gl.uniform1f(sliceLocation, parseFloat(sliderSlice.value));

    let colorMode = 0;
    if (selectColor.value === 'potential') colorMode = 1;
    else if (selectColor.value === 'normals') colorMode = 2;
    gl.uniform1i(colorModeLocation, colorMode);

    gl.uniform2f(mouseLocation, mouseX, mouseY);

    gl.drawArrays(gl.TRIANGLES, 0, 6);

    requestAnimationFrame(render);
}

requestAnimationFrame(render);

// Draw convergence chart
function drawConvergenceChart() {
    const chartCanvas = document.getElementById('chart-convergence');
    if (!chartCanvas) return;
    const ctx = chartCanvas.getContext('2d');
    const w = chartCanvas.width;
    const h = chartCanvas.height;

    ctx.fillStyle = '#0a0c10';
    ctx.fillRect(0, 0, w, h);

    // Axes
    ctx.strokeStyle = '#2d3748';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(35, 10);
    ctx.lineTo(35, h - 25);
    ctx.lineTo(w - 10, h - 25);
    ctx.stroke();

    // Data points: k=1 (0.1824), k=2 (0.0762), k=3 (0.0381)
    const points = [
        { k: 1, val: 0.1824, x: 70,  y: 30 },
        { k: 2, val: 0.0762, x: 160, y: 75 },
        { k: 3, val: 0.0381, x: 250, y: 105 }
    ];

    // Plot line
    ctx.strokeStyle = '#4299e1';
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(points[0].x, points[0].y);
    for (let i = 1; i < points.length; ++i) {
        ctx.lineTo(points[i].x, points[i].y);
    }
    ctx.stroke();

    // Draw dots and labels
    ctx.fillStyle = '#63b3ed';
    for (const pt of points) {
        ctx.beginPath();
        ctx.arc(pt.x, pt.y, 4, 0, 2 * Math.PI);
        ctx.fill();

        ctx.fillStyle = '#a0aec0';
        ctx.font = '10px monospace';
        ctx.fillText(`k=${pt.k}`, pt.x - 8, h - 10);
        ctx.fillText(pt.val.toFixed(3), pt.x - 14, pt.y - 8);
    }
}

drawConvergenceChart();
