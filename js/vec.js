/* Small vector / matrix helpers. Vectors are plain [x, y, z] arrays; 3x3 matrices are row-major
   arrays of 9 numbers. Everything returns new arrays, nothing mutates its arguments. */
(function () {
  'use strict';

  const D2R = Math.PI / 180;

  const V = {
    D2R,
    add: (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]],
    sub: (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]],
    scale: (a, s) => [a[0] * s, a[1] * s, a[2] * s],
    addScaled: (a, b, s) => [a[0] + b[0] * s, a[1] + b[1] * s, a[2] + b[2] * s],
    dot: (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2],
    cross: (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]],
    len: (a) => Math.hypot(a[0], a[1], a[2]),
    dist: (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]),
    dist2: (a, b) => { const x = a[0] - b[0], y = a[1] - b[1], z = a[2] - b[2]; return x * x + y * y + z * z; },
    lerp: (a, b, t) => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t],
    norm(a) {
      const l = Math.hypot(a[0], a[1], a[2]);
      return l > 1e-12 ? [a[0] / l, a[1] / l, a[2] / l] : [0, 0, 0];
    },
    copy: (a) => [a[0], a[1], a[2]],

    /* Component of a perpendicular to unit vector u. */
    perp: (a, u) => { const d = a[0] * u[0] + a[1] * u[1] + a[2] * u[2]; return [a[0] - u[0] * d, a[1] - u[1] * d, a[2] - u[2] * d]; },

    /* Any unit vector perpendicular to unit vector u. */
    anyPerp(u) {
      const t = Math.abs(u[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
      return V.norm(V.cross(u, t));
    },

    centroid(pts) {
      const c = [0, 0, 0];
      for (const p of pts) { c[0] += p[0]; c[1] += p[1]; c[2] += p[2]; }
      const n = Math.max(pts.length, 1);
      return [c[0] / n, c[1] / n, c[2] / n];
    },

    /* Dihedral p0-p1-p2-p3 in degrees (IUPAC sign: a right-handed rotation of p3 about p1->p2
       increases it). */
    dihedral(p0, p1, p2, p3) {
      const b1 = V.sub(p1, p0), b2 = V.sub(p2, p1), b3 = V.sub(p3, p2);
      const n1 = V.cross(b1, b2), n2 = V.cross(b2, b3);
      const y = V.dot(V.cross(n1, n2), V.norm(b2));
      return Math.atan2(y, V.dot(n1, n2)) / D2R;
    },

    angle(a, b, c) {
      const u = V.norm(V.sub(a, b)), v = V.norm(V.sub(c, b));
      return Math.acos(Math.max(-1, Math.min(1, V.dot(u, v)))) / D2R;
    },

    /* Wrap degrees into (-180, 180]. */
    wrap(d) { d = ((d + 180) % 360 + 360) % 360 - 180; return d === -180 ? 180 : d; },

    /* ---- 3x3 matrices ---- */
    I3: () => [1, 0, 0, 0, 1, 0, 0, 0, 1],
    m3vec: (m, v) => [
      m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
      m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
      m[6] * v[0] + m[7] * v[1] + m[8] * v[2]],
    m3mul(a, b) {
      const r = new Array(9);
      for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++)
        r[i * 3 + j] = a[i * 3] * b[j] + a[i * 3 + 1] * b[3 + j] + a[i * 3 + 2] * b[6 + j];
      return r;
    },
    m3T: (m) => [m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]],
    /* Matrix whose columns are a, b, c. */
    m3cols: (a, b, c) => [a[0], b[0], c[0], a[1], b[1], c[1], a[2], b[2], c[2]],
    col: (m, j) => [m[j], m[3 + j], m[6 + j]],

    /* Right-handed rotation by deg about unit axis u (Rodrigues). */
    rotAxis(u, deg) {
      u = V.norm(u);
      const t = deg * D2R, c = Math.cos(t), s = Math.sin(t), k = 1 - c;
      const [x, y, z] = u;
      return [
        c + x * x * k, x * y * k - z * s, x * z * k + y * s,
        y * x * k + z * s, c + y * y * k, y * z * k - x * s,
        z * x * k - y * s, z * y * k + x * s, c + z * z * k];
    },

    /* Rotate point p about the axis through o along u by deg. */
    rotateAbout(p, o, u, deg) { return V.add(o, V.m3vec(V.rotAxis(u, deg), V.sub(p, o))); },

    /* Rotation that takes orthonormal frame (a1, a2) onto (b1, b2): R a1 = b1, R a2 = b2. */
    frameToFrame(a1, a2, b1, b2) {
      a1 = V.norm(a1); a2 = V.norm(V.perp(a2, a1)); const a3 = V.cross(a1, a2);
      b1 = V.norm(b1); b2 = V.norm(V.perp(b2, b1)); const b3 = V.cross(b1, b2);
      return V.m3mul(V.m3cols(b1, b2, b3), V.m3T(V.m3cols(a1, a2, a3)));
    },

    /* Symmetric 3x3 eigen-decomposition (Jacobi). Returns {values, vectors} sorted by
       descending value; vectors[k] is the unit eigenvector of values[k]. */
    eigSym3(m) {
      const a = [[m[0], m[1], m[2]], [m[3], m[4], m[5]], [m[6], m[7], m[8]]];
      const v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
      for (let sweep = 0; sweep < 50; sweep++) {
        const off = Math.abs(a[0][1]) + Math.abs(a[0][2]) + Math.abs(a[1][2]);
        if (off < 1e-14) break;
        for (const [p, q] of [[0, 1], [0, 2], [1, 2]]) {
          if (Math.abs(a[p][q]) < 1e-18) continue;
          const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
          const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
          const c = 1 / Math.sqrt(t * t + 1), s = t * c;
          for (let k = 0; k < 3; k++) {
            const akp = a[k][p], akq = a[k][q];
            a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq;
          }
          for (let k = 0; k < 3; k++) {
            const apk = a[p][k], aqk = a[q][k];
            a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk;
          }
          for (let k = 0; k < 3; k++) {
            const vkp = v[k][p], vkq = v[k][q];
            v[k][p] = c * vkp - s * vkq; v[k][q] = s * vkp + c * vkq;
          }
        }
      }
      const idx = [0, 1, 2].sort((i, j) => a[j][j] - a[i][i]);
      return { values: idx.map(i => a[i][i]), vectors: idx.map(i => [v[0][i], v[1][i], v[2][i]]) };
    },

    /* Principal axes of a point cloud: {center, axes:[major, middle, minor], values}. */
    pca(pts) {
      const c = V.centroid(pts);
      const m = [0, 0, 0, 0, 0, 0, 0, 0, 0];
      for (const p of pts) {
        const d = V.sub(p, c);
        for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) m[i * 3 + j] += d[i] * d[j];
      }
      const e = V.eigSym3(m);
      return { center: c, axes: e.vectors, values: e.values };
    },

    /* Weighted Kabsch: rotation R and translation t minimising sum w |R x + t - y|^2.
       Returns {R, t}; apply as y ~ R x + t. */
    kabsch(xs, ys, ws) {
      let W = 0; const cx = [0, 0, 0], cy = [0, 0, 0];
      for (let i = 0; i < xs.length; i++) {
        const w = ws ? ws[i] : 1; W += w;
        for (let k = 0; k < 3; k++) { cx[k] += w * xs[i][k]; cy[k] += w * ys[i][k]; }
      }
      for (let k = 0; k < 3; k++) { cx[k] /= W; cy[k] /= W; }
      // Covariance H = sum w (x - cx)(y - cy)^T; use the quaternion (Horn) method for robustness.
      const H = [0, 0, 0, 0, 0, 0, 0, 0, 0];
      for (let i = 0; i < xs.length; i++) {
        const w = ws ? ws[i] : 1;
        const a = V.sub(xs[i], cx), b = V.sub(ys[i], cy);
        for (let p = 0; p < 3; p++) for (let q = 0; q < 3; q++) H[p * 3 + q] += w * a[p] * b[q];
      }
      const [Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz] = H;
      const N = [
        [Sxx + Syy + Szz, Syz - Szy, Szx - Sxz, Sxy - Syx],
        [Syz - Szy, Sxx - Syy - Szz, Sxy + Syx, Szx + Sxz],
        [Szx - Sxz, Sxy + Syx, -Sxx + Syy - Szz, Syz + Szy],
        [Sxy - Syx, Szx + Sxz, Syz + Szy, -Sxx - Syy + Szz]];
      const q = V._maxEig4(N);
      const [w0, x, y, z] = q;
      const R = [
        w0 * w0 + x * x - y * y - z * z, 2 * (x * y - w0 * z), 2 * (x * z + w0 * y),
        2 * (x * y + w0 * z), w0 * w0 - x * x + y * y - z * z, 2 * (y * z - w0 * x),
        2 * (x * z - w0 * y), 2 * (y * z + w0 * x), w0 * w0 - x * x - y * y + z * z];
      const t = V.sub(cy, V.m3vec(R, cx));
      return { R, t };
    },

    /* Eigenvector of the largest eigenvalue of a symmetric 4x4 matrix (Jacobi). */
    _maxEig4(A) {
      const a = A.map(r => r.slice());
      const v = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
      for (let sweep = 0; sweep < 60; sweep++) {
        let off = 0;
        for (let p = 0; p < 4; p++) for (let q = p + 1; q < 4; q++) off += Math.abs(a[p][q]);
        if (off < 1e-14) break;
        for (let p = 0; p < 4; p++) for (let q = p + 1; q < 4; q++) {
          if (Math.abs(a[p][q]) < 1e-18) continue;
          const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
          const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
          const c = 1 / Math.sqrt(t * t + 1), s = t * c;
          for (let k = 0; k < 4; k++) {
            const akp = a[k][p], akq = a[k][q];
            a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq;
          }
          for (let k = 0; k < 4; k++) {
            const apk = a[p][k], aqk = a[q][k];
            a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk;
          }
          for (let k = 0; k < 4; k++) {
            const vkp = v[k][p], vkq = v[k][q];
            v[k][p] = c * vkp - s * vkq; v[k][q] = s * vkp + c * vkq;
          }
        }
      }
      let best = 0;
      for (let i = 1; i < 4; i++) if (a[i][i] > a[best][best]) best = i;
      const q = [v[0][best], v[1][best], v[2][best], v[3][best]];
      const l = Math.hypot(...q);
      return q.map(x => x / l);
    },

    /* Closest-approach parameter s along line (p + s u) to the ray (o + t d). Returns null when
       the two are nearly parallel. */
    lineRayParam(p, u, o, d) {
      const w0 = V.sub(p, o);
      const a = V.dot(u, u), b = V.dot(u, d), c = V.dot(d, d);
      const dd = V.dot(u, w0), e = V.dot(d, w0);
      const den = a * c - b * b;
      if (Math.abs(den) < 1e-6) return null;
      return (b * e - c * dd) / den;
    },

    /* Ray / plane intersection; returns the point or null. */
    rayPlane(o, d, p0, n) {
      const den = V.dot(d, n);
      if (Math.abs(den) < 1e-9) return null;
      const t = V.dot(V.sub(p0, o), n) / den;
      return t > 0 ? V.addScaled(o, d, t) : null;
    },
  };

  window.V = V;
})();
