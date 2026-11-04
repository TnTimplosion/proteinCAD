/* WebGL2 renderer and camera.

   Draws instanced spheres (atoms), instanced cylinders (bonds) and triangle meshes (cartoon).
   Everything else that is drawn (sketch plane, gizmo, strokes, labels, chain breaks) goes on the
   2D overlay canvas, using the same camera through project() / ray(). */
(function () {
  'use strict';

  /* ---------------- camera ---------------- */
  function Camera() {
    this.target = [0, 0, 0];
    this.R = V.I3();          // columns: camera right, up, back (world)
    this.dist = 60;
    this.fov = 30;            // degrees, vertical
    this.ortho = false;
    this.w = 1; this.h = 1;
    this.sceneRadius = 30;
  }
  Camera.prototype = {
    right() { return V.col(this.R, 0); },
    up() { return V.col(this.R, 1); },
    back() { return V.col(this.R, 2); },
    eye() { return V.addScaled(this.target, this.back(), this.dist); },
    halfH() { return this.dist * Math.tan(this.fov * V.D2R / 2); },
    nearFar() {
      const near = Math.max(0.5, this.dist - this.sceneRadius * 1.6);
      const far = this.dist + this.sceneRadius * 1.6 + 10;
      return [near, far];
    },
    /* World -> view space (camera looks down -z). */
    toView(p) {
      const d = V.sub(p, this.eye());
      return [V.dot(d, this.right()), V.dot(d, this.up()), V.dot(d, this.back())];
    },
    viewMatrix() {
      const r = this.right(), u = this.up(), b = this.back(), e = this.eye();
      return new Float32Array([
        r[0], u[0], b[0], 0,
        r[1], u[1], b[1], 0,
        r[2], u[2], b[2], 0,
        -V.dot(r, e), -V.dot(u, e), -V.dot(b, e), 1]);
    },
    projMatrix() {
      const [n, f] = this.nearFar();
      const a = this.w / this.h;
      if (this.ortho) {
        const hh = this.halfH(), hw = hh * a;
        return new Float32Array([1 / hw, 0, 0, 0, 0, 1 / hh, 0, 0, 0, 0, -2 / (f - n), 0, 0, 0, -(f + n) / (f - n), 1]);
      }
      const t = 1 / Math.tan(this.fov * V.D2R / 2);
      return new Float32Array([t / a, 0, 0, 0, 0, t, 0, 0, 0, 0, (f + n) / (n - f), -1, 0, 0, 2 * f * n / (n - f), 0]);
    },
    /* World -> screen (CSS px). Returns [x, y, depth] or null when behind the camera. */
    project(p) {
      const v = this.toView(p);
      const [n] = this.nearFar();
      if (!this.ortho && v[2] > -n * 0.5) return null;
      const hh = this.halfH(), a = this.w / this.h;
      let nx, ny;
      if (this.ortho) { nx = v[0] / (hh * a); ny = v[1] / hh; }
      else { const s = this.dist / -v[2]; nx = v[0] * s / (hh * a); ny = v[1] * s / hh; }
      return [(nx + 1) * 0.5 * this.w, (1 - ny) * 0.5 * this.h, -v[2]];
    },
    /* Screen point -> world ray {o, d}. */
    ray(sx, sy) {
      const a = this.w / this.h, hh = this.halfH();
      const nx = (sx / this.w) * 2 - 1, ny = 1 - (sy / this.h) * 2;
      const r = this.right(), u = this.up(), b = this.back();
      if (this.ortho) {
        const o = V.add(V.addScaled(V.addScaled(this.eye(), r, nx * hh * a), u, ny * hh), V.scale(b, 0));
        return { o, d: V.scale(b, -1) };
      }
      const dirV = [nx * hh * a / this.dist, ny * hh / this.dist, -1];
      const d = V.norm(V.add(V.add(V.scale(r, dirV[0]), V.scale(u, dirV[1])), V.scale(b, dirV[2])));
      return { o: this.eye(), d };
    },
    /* World size of one CSS pixel at point p. */
    pixelSize(p) {
      const hh = this.halfH();
      if (this.ortho) return 2 * hh / this.h;
      const z = -this.toView(p)[2];
      return 2 * hh * (z / this.dist) / this.h;
    },
    orbit(dx, dy) {
      const k = 0.008;
      const Ry = V.rotAxis(this.up(), -dx * k * 57.2958);
      const Rx = V.rotAxis(this.right(), -dy * k * 57.2958);
      this.R = orthonormalize(V.m3mul(Rx, V.m3mul(Ry, this.R)));
    },
    pan(dx, dy) {
      const s = this.pixelSize(this.target);
      this.target = V.add(this.target, V.add(V.scale(this.right(), -dx * s), V.scale(this.up(), dy * s)));
    },
    zoom(f) { this.dist = Math.max(3, Math.min(5000, this.dist * f)); },
    /* Look along -dir (dir points from target to the eye), keeping `up` as close as possible. */
    lookAlong(dir, up) {
      const b = V.norm(dir);
      let u = V.perp(up || [0, 1, 0], b);
      if (V.len(u) < 1e-3) u = V.perp([0, 0, 1], b);
      u = V.norm(u);
      this.R = V.m3cols(V.cross(u, b), u, b);
    },
    fit(center, radius) {
      this.target = center.slice();
      this.sceneRadius = Math.max(this.sceneRadius, radius);
      this.dist = Math.max(10, radius / Math.tan(this.fov * V.D2R / 2) * 1.15);
    },
  };

  function orthonormalize(m) {
    let x = V.col(m, 0), y = V.col(m, 1);
    x = V.norm(x); y = V.norm(V.perp(y, x));
    return V.m3cols(x, y, V.cross(x, y));
  }

  /* ---------------- geometry primitives ---------------- */
  function icosphere(sub) {
    const t = (1 + Math.sqrt(5)) / 2;
    let verts = [[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t], [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]].map(V.norm);
    let faces = [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8], [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]];
    for (let s = 0; s < sub; s++) {
      const cache = new Map(), nf = [];
      const mid = (a, b) => {
        const k = a < b ? a + '_' + b : b + '_' + a;
        if (cache.has(k)) return cache.get(k);
        verts.push(V.norm(V.scale(V.add(verts[a], verts[b]), 0.5)));
        cache.set(k, verts.length - 1);
        return verts.length - 1;
      };
      for (const [a, b, c] of faces) {
        const ab = mid(a, b), bc = mid(b, c), ca = mid(c, a);
        nf.push([a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]);
      }
      faces = nf;
    }
    return { pos: new Float32Array(verts.flat()), idx: new Uint16Array(faces.flat()) };
  }

  function cylinder(seg) {
    const pos = [], nrm = [], idx = [];
    for (let i = 0; i <= seg; i++) {
      const a = i / seg * Math.PI * 2, c = Math.cos(a), s = Math.sin(a);
      pos.push(c, 0, s, c, 1, s); nrm.push(c, 0, s, c, 0, s);
    }
    for (let i = 0; i < seg; i++) { const k = i * 2; idx.push(k, k + 1, k + 2, k + 1, k + 3, k + 2); }
    return { pos: new Float32Array(pos), nrm: new Float32Array(nrm), idx: new Uint16Array(idx) };
  }

  /* ---------------- shaders ---------------- */
  const FRAG = `#version 300 es
  precision highp float;
  in vec3 vN; in vec3 vV; in vec4 vC;
  uniform vec3 uFog; uniform vec2 uFogRange;
  out vec4 o;
  void main(){
    vec3 n = normalize(vN);
    if (!gl_FrontFacing) n = -n;
    vec3 L = normalize(vec3(0.35, 0.45, 1.0));
    vec3 E = normalize(-vV);
    float d = max(dot(n, L), 0.0);
    float s = pow(max(dot(n, normalize(L + E)), 0.0), 48.0);
    vec3 base = vC.rgb;
    float st = vC.a;
    vec3 emit = vec3(0.0);
    if (st > 0.5 && st < 1.5) { base = mix(base, vec3(1.0, 0.82, 0.18), 0.6); emit = vec3(0.18, 0.13, 0.0); }
    else if (st > 1.5 && st < 2.5) { base = mix(base, vec3(1.0), 0.35); }
    else if (st > 2.5) { base = mix(base, vec3(0.35, 0.9, 1.0), 0.7); emit = vec3(0.0, 0.08, 0.1); }
    vec3 col = base * (0.32 + 0.68 * d) + vec3(0.28) * s + emit;
    float z = -vV.z;
    float f = clamp((z - uFogRange.x) / max(uFogRange.y - uFogRange.x, 1e-3), 0.0, 1.0);
    col = mix(col, uFog, f * 0.55);
    o = vec4(col, 1.0);
  }`;

  const VS_SPHERE = `#version 300 es
  layout(location=0) in vec3 aP;
  layout(location=1) in vec4 iCR;
  layout(location=2) in vec4 iCol;
  uniform mat4 uV, uP;
  out vec3 vN; out vec3 vV; out vec4 vC;
  void main(){
    vec4 v = uV * vec4(iCR.xyz + aP * iCR.w, 1.0);
    vV = v.xyz; vN = mat3(uV) * aP; vC = iCol;
    gl_Position = uP * v;
  }`;

  const VS_CYL = `#version 300 es
  layout(location=0) in vec3 aP;
  layout(location=1) in vec3 aN;
  layout(location=2) in vec4 iA;
  layout(location=3) in vec3 iB;
  layout(location=4) in vec4 iCol;
  uniform mat4 uV, uP;
  out vec3 vN; out vec3 vV; out vec4 vC;
  void main(){
    vec3 ax = iB - iA.xyz; float L = length(ax);
    vec3 y = ax / max(L, 1e-6);
    vec3 t = abs(y.x) < 0.9 ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 1.0, 0.0);
    vec3 x = normalize(cross(t, y)); vec3 z = cross(x, y);
    vec3 wp = iA.xyz + x * (aP.x * iA.w) + y * (aP.y * L) + z * (aP.z * iA.w);
    vec4 v = uV * vec4(wp, 1.0);
    vV = v.xyz; vN = mat3(uV) * (x * aN.x + z * aN.z); vC = iCol;
    gl_Position = uP * v;
  }`;

  const VS_MESH = `#version 300 es
  layout(location=0) in vec3 aP;
  layout(location=1) in vec3 aN;
  layout(location=2) in vec4 aCol;
  uniform mat4 uV, uP;
  out vec3 vN; out vec3 vV; out vec4 vC;
  void main(){
    vec4 v = uV * vec4(aP, 1.0);
    vV = v.xyz; vN = mat3(uV) * aN; vC = aCol;
    gl_Position = uP * v;
  }`;

  function compile(gl, vs, fs) {
    const mk = (type, src) => {
      const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s);
      if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
      return s;
    };
    const p = gl.createProgram();
    gl.attachShader(p, mk(gl.VERTEX_SHADER, vs)); gl.attachShader(p, mk(gl.FRAGMENT_SHADER, fs));
    gl.linkProgram(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(p));
    return { p, uV: gl.getUniformLocation(p, 'uV'), uP: gl.getUniformLocation(p, 'uP'), uFog: gl.getUniformLocation(p, 'uFog'), uFogRange: gl.getUniformLocation(p, 'uFogRange') };
  }

  /* ---------------- renderer ---------------- */
  function Renderer(canvas) {
    const gl = canvas.getContext('webgl2', { antialias: true, alpha: false });
    if (!gl) throw new Error('WebGL2 is not available in this browser.');
    this.gl = gl; this.canvas = canvas;
    this.cam = new Camera();
    this.bg = [0.105, 0.118, 0.141];
    this.progS = compile(gl, VS_SPHERE, FRAG);
    this.progC = compile(gl, VS_CYL, FRAG);
    this.progM = compile(gl, VS_MESH, FRAG);

    const sph = icosphere(2);
    this.sph = { vao: gl.createVertexArray(), n: sph.idx.length, inst: gl.createBuffer(), count: 0 };
    gl.bindVertexArray(this.sph.vao);
    this._buf(gl.ARRAY_BUFFER, sph.pos); gl.enableVertexAttribArray(0); gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
    this._buf(gl.ELEMENT_ARRAY_BUFFER, sph.idx);
    gl.bindBuffer(gl.ARRAY_BUFFER, this.sph.inst);
    gl.enableVertexAttribArray(1); gl.vertexAttribPointer(1, 4, gl.FLOAT, false, 32, 0); gl.vertexAttribDivisor(1, 1);
    gl.enableVertexAttribArray(2); gl.vertexAttribPointer(2, 4, gl.FLOAT, false, 32, 16); gl.vertexAttribDivisor(2, 1);

    const cyl = cylinder(14);
    this.cyl = { vao: gl.createVertexArray(), n: cyl.idx.length, inst: gl.createBuffer(), count: 0 };
    gl.bindVertexArray(this.cyl.vao);
    this._buf(gl.ARRAY_BUFFER, cyl.pos); gl.enableVertexAttribArray(0); gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
    this._buf(gl.ARRAY_BUFFER, cyl.nrm); gl.enableVertexAttribArray(1); gl.vertexAttribPointer(1, 3, gl.FLOAT, false, 0, 0);
    this._buf(gl.ELEMENT_ARRAY_BUFFER, cyl.idx);
    gl.bindBuffer(gl.ARRAY_BUFFER, this.cyl.inst);
    gl.enableVertexAttribArray(2); gl.vertexAttribPointer(2, 4, gl.FLOAT, false, 44, 0); gl.vertexAttribDivisor(2, 1);
    gl.enableVertexAttribArray(3); gl.vertexAttribPointer(3, 3, gl.FLOAT, false, 44, 16); gl.vertexAttribDivisor(3, 1);
    gl.enableVertexAttribArray(4); gl.vertexAttribPointer(4, 4, gl.FLOAT, false, 44, 28); gl.vertexAttribDivisor(4, 1);

    this.mesh = { vao: gl.createVertexArray(), vbo: gl.createBuffer(), ibo: gl.createBuffer(), n: 0 };
    gl.bindVertexArray(this.mesh.vao);
    gl.bindBuffer(gl.ARRAY_BUFFER, this.mesh.vbo);
    gl.enableVertexAttribArray(0); gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 40, 0);
    gl.enableVertexAttribArray(1); gl.vertexAttribPointer(1, 3, gl.FLOAT, false, 40, 12);
    gl.enableVertexAttribArray(2); gl.vertexAttribPointer(2, 4, gl.FLOAT, false, 40, 24);
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.mesh.ibo);
    gl.bindVertexArray(null);
  }

  Renderer.prototype = {
    _buf(target, data) {
      const gl = this.gl, b = gl.createBuffer();
      gl.bindBuffer(target, b); gl.bufferData(target, data, gl.STATIC_DRAW);
      return b;
    },
    resize(w, h, dpr) {
      this.canvas.width = Math.max(1, Math.round(w * dpr));
      this.canvas.height = Math.max(1, Math.round(h * dpr));
      this.cam.w = w; this.cam.h = h;
    },
    /* scene: { spheres: Float32Array (8 per), cylinders: Float32Array (11 per),
                mesh: {verts: Float32Array (10 per), idx: Uint32Array} } */
    upload(scene) {
      const gl = this.gl;
      gl.bindBuffer(gl.ARRAY_BUFFER, this.sph.inst);
      gl.bufferData(gl.ARRAY_BUFFER, scene.spheres, gl.DYNAMIC_DRAW);
      this.sph.count = scene.spheres.length / 8;
      gl.bindBuffer(gl.ARRAY_BUFFER, this.cyl.inst);
      gl.bufferData(gl.ARRAY_BUFFER, scene.cylinders, gl.DYNAMIC_DRAW);
      this.cyl.count = scene.cylinders.length / 11;
      gl.bindVertexArray(this.mesh.vao);
      gl.bindBuffer(gl.ARRAY_BUFFER, this.mesh.vbo);
      gl.bufferData(gl.ARRAY_BUFFER, scene.mesh.verts, gl.DYNAMIC_DRAW);
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.mesh.ibo);
      gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, scene.mesh.idx, gl.DYNAMIC_DRAW);
      this.mesh.n = scene.mesh.idx.length;
      gl.bindVertexArray(null);
    },
    draw() {
      const gl = this.gl, cam = this.cam;
      gl.viewport(0, 0, this.canvas.width, this.canvas.height);
      gl.clearColor(this.bg[0], this.bg[1], this.bg[2], 1);
      gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
      gl.enable(gl.DEPTH_TEST);
      gl.disable(gl.CULL_FACE);
      const Vm = cam.viewMatrix(), Pm = cam.projMatrix();
      const fogNear = Math.max(0, cam.dist), fogFar = cam.dist + cam.sceneRadius * 1.5;
      const use = (pr) => {
        gl.useProgram(pr.p);
        gl.uniformMatrix4fv(pr.uV, false, Vm); gl.uniformMatrix4fv(pr.uP, false, Pm);
        gl.uniform3fv(pr.uFog, this.bg); gl.uniform2f(pr.uFogRange, fogNear, fogFar);
      };
      if (this.sph.count) { use(this.progS); gl.bindVertexArray(this.sph.vao); gl.drawElementsInstanced(gl.TRIANGLES, this.sph.n, gl.UNSIGNED_SHORT, 0, this.sph.count); }
      if (this.cyl.count) { use(this.progC); gl.bindVertexArray(this.cyl.vao); gl.drawElementsInstanced(gl.TRIANGLES, this.cyl.n, gl.UNSIGNED_SHORT, 0, this.cyl.count); }
      if (this.mesh.n) { use(this.progM); gl.bindVertexArray(this.mesh.vao); gl.drawElements(gl.TRIANGLES, this.mesh.n, gl.UNSIGNED_INT, 0); }
      gl.bindVertexArray(null);
    },
  };

  /* ---------------- scene building ---------------- */
  const ELEMENT_COLOR = { C: null, N: [0.27, 0.45, 0.95], O: [0.93, 0.25, 0.22], S: [0.95, 0.82, 0.25], P: [1.0, 0.55, 0.1], Se: [0.95, 0.6, 0.1] };
  const SS_COLOR = { H: [0.91, 0.33, 0.42], E: [0.98, 0.78, 0.25], L: [0.62, 0.68, 0.72] };
  const CHAIN_PALETTE = [
    [0.33, 0.62, 0.95], [0.40, 0.80, 0.50], [0.95, 0.60, 0.30], [0.75, 0.50, 0.92], [0.35, 0.82, 0.82],
    [0.93, 0.45, 0.62], [0.78, 0.78, 0.40], [0.55, 0.65, 0.98], [0.95, 0.50, 0.45], [0.55, 0.85, 0.65]];
  function chainColor(id) {
    const k = Model.CHAIN_IDS.indexOf(id);
    let h = k;
    if (k < 0) { h = 0; for (const c of id) h = (h * 31 + c.charCodeAt(0)) >>> 0; }
    return CHAIN_PALETTE[h % CHAIN_PALETTE.length];
  }
  const LIGAND_C = [0.55, 0.85, 0.45];

  /* Bonds inside a residue by distance, plus the peptide bond to the next residue. */
  function residueBonds(res) {
    const out = [];
    const A = res.atoms;
    for (let i = 0; i < A.length; i++) for (let j = i + 1; j < A.length; j++) {
      const cut = (A[i].el === 'S' || A[j].el === 'S' || A[i].el === 'Se' || A[j].el === 'Se') ? 2.2 : 1.9;
      if (V.dist2(A[i].p, A[j].p) < cut * cut) out.push([A[i], A[j]]);
    }
    return out;
  }

  /* Build GPU arrays for the model. opts: { atoms: 'none'|'backbone'|'all', style: 'sticks'|'ball',
     cartoon, trace, color: 'chain'|'ss'|'element', sel: Set, hover: uid, ghost: [[res]], hidden: Set(chainId) } */
  function buildScene(model, opts) {
    const sph = [], cyl = [];
    const verts = [], idx = [];
    const BBSET = new Set(['N', 'CA', 'C', 'O', 'CB', 'OXT']);
    const atomR = opts.style === 'ball' ? 0.38 : 0.22, bondR = opts.style === 'ball' ? 0.12 : 0.2;

    const resColor = (res, chainId) => {
      if (res.het) return LIGAND_C;
      if (opts.color === 'ss') return SS_COLOR[res.ss] || SS_COLOR.L;
      return chainColor(chainId);
    };
    const stateOf = (res, ghost) => ghost ? 3 : opts.sel.has(res.uid) ? 1 : opts.hover === res.uid ? 2 : 0;
    const atomColor = (a, base) => (opts.color === 'element' || opts.color === 'chain') && a.el !== 'C' && ELEMENT_COLOR[a.el] ? ELEMENT_COLOR[a.el] : base;

    const drawResidueAtoms = (res, base, st, mode) => {
      const keep = a => mode === 'all' || BBSET.has(a.name) || res.het;
      for (const a of res.atoms) {
        if (!keep(a)) continue;
        const c = atomColor(a, base);
        const r = (a.name === 'CA' && !res.het) ? atomR * 1.25 : atomR;
        sph.push(a.p[0], a.p[1], a.p[2], opts.style === 'sticks' ? Math.max(r, bondR) : r, c[0], c[1], c[2], st);
      }
      for (const [a, b] of residueBonds(res)) {
        if (!keep(a) || !keep(b)) continue;
        const ca = atomColor(a, base), cb = atomColor(b, base);
        const m = V.lerp(a.p, b.p, 0.5);
        cyl.push(a.p[0], a.p[1], a.p[2], bondR, m[0], m[1], m[2], ca[0], ca[1], ca[2], st);
        cyl.push(m[0], m[1], m[2], bondR, b.p[0], b.p[1], b.p[2], cb[0], cb[1], cb[2], st);
      }
    };

    const chains = model.chains.filter(c => !opts.hidden.has(c.id)).map(c => ({ chain: c, ghost: false }));
    for (const g of (opts.ghost || [])) chains.push({ chain: { id: '_', residues: g }, ghost: true });

    for (const { chain, ghost } of chains) {
      const rs = chain.residues;
      const atomsMode = ghost ? (opts.atoms === 'none' ? 'backbone' : opts.atoms) : opts.atoms;
      for (let i = 0; i < rs.length; i++) {
        const res = rs[i];
        const base = resColor(res, chain.id), st = stateOf(res, ghost);
        if (atomsMode !== 'none') {
          drawResidueAtoms(res, base, st, atomsMode);
          const nx = rs[i + 1];
          if (nx && Model.bonded(res, nx)) {
            const C = Geom.pos(res, 'C'), N = Geom.pos(nx, 'N');
            const m = V.lerp(C, N, 0.5);
            const c1 = base, c2 = ELEMENT_COLOR.N && opts.color !== 'ss' ? ELEMENT_COLOR.N : resColor(nx, chain.id);
            cyl.push(C[0], C[1], C[2], bondR, m[0], m[1], m[2], c1[0], c1[1], c1[2], st);
            cyl.push(m[0], m[1], m[2], bondR, N[0], N[1], N[2], c2[0], c2[1], c2[2], stateOf(nx, ghost));
          }
        }
        if (opts.trace || (atomsMode === 'none' && !opts.cartoon)) {
          const CA = Geom.pos(res, 'CA');
          if (CA) {
            sph.push(CA[0], CA[1], CA[2], 0.55, base[0], base[1], base[2], st);
            const nx = rs[i + 1], CA2 = nx && Geom.pos(nx, 'CA');
            if (CA2 && V.dist(CA, CA2) < 4.3) {
              const m = V.lerp(CA, CA2, 0.5), b2 = resColor(nx, chain.id);
              cyl.push(CA[0], CA[1], CA[2], 0.28, m[0], m[1], m[2], base[0], base[1], base[2], st);
              cyl.push(m[0], m[1], m[2], 0.28, CA2[0], CA2[1], CA2[2], b2[0], b2[1], b2[2], stateOf(nx, ghost));
            }
          }
        }
      }
      if (opts.cartoon) cartoon(rs, res => resColor(res, chain.id), res => stateOf(res, ghost), verts, idx);
    }
    for (const res of model.ligands) {
      if (res.chainId && opts.hidden.has('lig')) continue;
      drawResidueAtoms(res, LIGAND_C, stateOf(res, false), 'all');
    }
    return { spheres: new Float32Array(sph), cylinders: new Float32Array(cyl), mesh: { verts: new Float32Array(verts), idx: new Uint32Array(idx) } };
  }

  /* ---------------- cartoon ---------------- */
  function catmull(p0, p1, p2, p3, t) {
    const t2 = t * t, t3 = t2 * t;
    const f = (a, b, c, d) => 0.5 * ((2 * b) + (-a + c) * t + (2 * a - 5 * b + 4 * c - d) * t2 + (-a + 3 * b - 3 * c + d) * t3);
    return [f(p0[0], p1[0], p2[0], p3[0]), f(p0[1], p1[1], p2[1], p3[1]), f(p0[2], p1[2], p2[2], p3[2])];
  }

  /* Ribbon through CA atoms: flat for helix and strand (width along the peptide plane, from the
     C=O direction), round tube for loops, an arrowhead at the end of each strand. */
  function cartoon(rs, colorOf, stateOf, verts, idx) {
    const SEG = 7, RING = 10;
    // Split into continuous segments.
    const segs = []; let cur = [];
    for (let i = 0; i < rs.length; i++) {
      const r = rs[i];
      const CA = Geom.pos(r, 'CA'), O = Geom.pos(r, 'O'), C = Geom.pos(r, 'C');
      if (!CA) { if (cur.length) segs.push(cur); cur = []; continue; }
      if (cur.length) {
        const prev = cur[cur.length - 1];
        if (V.dist(prev.CA, CA) > 4.3) { segs.push(cur); cur = []; }
      }
      cur.push({ r, CA, O: O && C ? V.norm(V.sub(O, C)) : null });
    }
    if (cur.length) segs.push(cur);

    for (const seg of segs) {
      const n = seg.length;
      if (n < 2) continue;
      // Orientation vectors, made consistent along the chain.
      let prevO = null;
      for (let i = 0; i < n; i++) {
        let o = seg[i].O;
        if (!o) {
          const a = seg[Math.max(0, i - 1)].CA, b = seg[Math.min(n - 1, i + 1)].CA;
          o = V.anyPerp(V.norm(V.sub(b, a)));
        }
        if (prevO && V.dot(o, prevO) < 0) o = V.scale(o, -1);
        seg[i].o = o; prevO = o;
      }
      const P = i => seg[Math.max(0, Math.min(n - 1, i))].CA;
      const dims = i => {
        const ss = seg[i].r.ss;
        if (ss === 'H') return [1.25, 0.28];
        if (ss === 'E') return [1.15, 0.28];
        return [0.32, 0.32];
      };
      const rings = [];
      for (let i = 0; i < n - 1; i++) {
        const lastStrand = seg[i].r.ss === 'E' && (i + 1 >= n || seg[i + 1].r.ss !== 'E');
        const arrow = seg[i].r.ss === 'E' && seg[i + 1].r.ss !== 'E';
        for (let s = 0; s < SEG; s++) {
          const t = s / SEG;
          const p = catmull(P(i - 1), P(i), P(i + 1), P(i + 2), t);
          const p2 = catmull(P(i - 1), P(i), P(i + 1), P(i + 2), Math.min(1, t + 0.01));
          let T = V.norm(V.sub(p2, p));
          if (V.len(T) < 1e-6) T = V.norm(V.sub(P(i + 1), P(i)));
          let Nn = V.norm(V.perp(V.lerp(seg[i].o, seg[i + 1].o, t), T));
          if (V.len(Nn) < 1e-6) Nn = V.anyPerp(T);
          const B = V.cross(T, Nn);
          const d0 = dims(i), d1 = dims(i + 1);
          let w, h;
          if (arrow) { w = 1.9 * (1 - t) + 0.32 * t; h = 0.28 * (1 - t) + 0.32 * t; }
          else { const k = t * t * (3 - 2 * t); w = d0[0] + (d1[0] - d0[0]) * k; h = d0[1] + (d1[1] - d0[1]) * k; }
          const r = t < 0.5 ? seg[i].r : seg[i + 1].r;
          rings.push({ p, Nn, B, w, h, col: colorOf(r), st: stateOf(r) });
          void lastStrand;
        }
      }
      // Final ring at the last CA.
      {
        const i = n - 1;
        const T = V.norm(V.sub(P(i), P(i - 1)));
        const Nn = V.norm(V.perp(seg[i].o, T));
        const d = dims(i);
        rings.push({ p: P(i), Nn, B: V.cross(T, Nn), w: d[0], h: d[1], col: colorOf(seg[i].r), st: stateOf(seg[i].r) });
      }
      const base0 = verts.length / 10;
      for (const rg of rings) {
        for (let k = 0; k < RING; k++) {
          const a = k / RING * Math.PI * 2, c = Math.cos(a), s = Math.sin(a);
          const p = V.add(rg.p, V.add(V.scale(rg.Nn, c * rg.w), V.scale(rg.B, s * rg.h)));
          const nn = V.norm(V.add(V.scale(rg.Nn, c * rg.h), V.scale(rg.B, s * rg.w)));
          verts.push(p[0], p[1], p[2], nn[0], nn[1], nn[2], rg.col[0], rg.col[1], rg.col[2], rg.st);
        }
      }
      for (let r = 0; r + 1 < rings.length; r++) {
        for (let k = 0; k < RING; k++) {
          const a = base0 + r * RING + k, b = base0 + r * RING + (k + 1) % RING;
          const c = a + RING, d = b + RING;
          idx.push(a, c, b, b, c, d);
        }
      }
    }
  }

  window.Render = { Camera, Renderer, buildScene, chainColor, SS_COLOR, ELEMENT_COLOR, residueBonds };
})();
