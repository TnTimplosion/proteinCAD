/* Builders: secondary-structure elements placed relative to a sketch plane, loop closure between
   two chain ends, and fitting a loop to a freehand curve.

   Loops are solved with cyclic coordinate descent (CCD, Canutescu & Dunbrack 2003) generalised to
   many weighted targets: each rotatable backbone bond in turn is rotated to the angle that best
   brings every downstream target onto its goal (closed form), with phi/psi kept inside broad
   allowed Ramachandran regions. Many random starts, scored on closure, fit and clashes. */
(function () {
  'use strict';
  const { BB, CONF } = Geom;

  /* ---------------- helix extrusion ---------------- */

  /* opts: { n, conf, base (point on plane), dir (unit; extrusion direction), inPlane (unit vector
     in the plane, phase reference), phase (deg), endAtPlane (C-terminus sits on the plane), name } */
  function helix(opts) {
    const n = Math.max(1, opts.n | 0);
    const { bb, axis } = Geom.idealSegment(n, opts.conf);
    const iRef = opts.endAtPlane ? n - 1 : 0;
    const a0 = Geom.projectOnAxis(bb[iRef].CA, axis);
    const w = V.sub(bb[iRef].CA, a0);
    const travel = opts.endAtPlane ? V.scale(opts.dir, -1) : opts.dir;
    const ref = V.m3vec(V.rotAxis(travel, opts.phase || 0), opts.inPlane);
    const R = V.frameToFrame(axis.u, w, travel, ref);
    const t = V.sub(opts.base, V.m3vec(R, a0));
    const moved = bb.map(r => ({ N: V.add(V.m3vec(R, r.N), t), CA: V.add(V.m3vec(R, r.CA), t), C: V.add(V.m3vec(R, r.C), t) }));
    return Geom.residuesFromBackbone(moved, opts.name || 'ALA', opts.conf.ss, opts.conf.psi);
  }

  function riseOf(conf) {
    return Geom.idealSegment(3, conf).axis.rise;
  }

  /* ---------------- beta sheet ---------------- */

  /* opts: { p0, p1 (first strand drawn on the plane), normal, nStrands, parallel, spacing (A between
     strand axes), side (+1/-1: which side of the first strand the others go), twist (deg between
     neighbouring strands, about the inter-strand axis), shift (register shift, residues per strand),
     name }. Returns one residue array per strand, in order across the sheet. */
  function sheet(opts) {
    const conf = opts.parallel ? CONF.parallel : CONF.antiparallel;
    const nrm = V.norm(opts.normal);
    const dRaw = V.perp(V.sub(opts.p1, opts.p0), nrm);
    const L = V.len(dRaw);
    const d = L > 1e-6 ? V.norm(dRaw) : V.anyPerp(nrm);
    const rise = riseOf(conf);
    const nres = opts.nres || Math.max(2, Math.round(L / rise) + 1);
    const s = V.scale(V.norm(V.cross(nrm, d)), opts.side < 0 ? -1 : 1);
    const c0 = V.addScaled(opts.p0, d, L / 2);
    const m = Math.floor((nres - 1) / 2);
    const out = [];
    for (let k = 0; k < Math.max(1, opts.nStrands | 0); k++) {
      const fwd = opts.parallel || k % 2 === 0;
      const { bb, axis } = Geom.idealSegment(nres, conf);
      const mk = fwd ? m : nres - 1 - m;
      const aM = Geom.projectOnAxis(bb[mk].CA, axis);
      const cb = Geom.idealCB(bb[mk].N, bb[mk].CA, bb[mk].C);
      const w = V.sub(cb, aM);
      const shiftK = (opts.shift || 0) * k;
      const up = (((shiftK % 2) + 2) % 2 === 0) ? 1 : -1;
      const Rt = V.rotAxis(s, (opts.twist || 0) * k);
      const X = V.m3vec(Rt, fwd ? d : V.scale(d, -1));
      const Y = V.m3vec(Rt, V.scale(nrm, up));
      const R = V.frameToFrame(axis.u, w, X, Y);
      const center = V.addScaled(V.addScaled(c0, s, (opts.spacing || 4.8) * k), d, rise * shiftK);
      const t = V.sub(center, V.m3vec(R, aM));
      const moved = bb.map(r => ({ N: V.add(V.m3vec(R, r.N), t), CA: V.add(V.m3vec(R, r.CA), t), C: V.add(V.m3vec(R, r.C), t) }));
      out.push(Geom.residuesFromBackbone(moved, opts.name || 'ALA', 'E', conf.psi));
    }
    return out;
  }

  /* ---------------- CCD engine ---------------- */

  /* Allowed (phi, psi) for non-glycine loop residues: the broad alpha-R / beta / PPII basins plus
     alpha-L (phi +45..+90), which tight turns and 2-residue hairpins need. In designed proteins
     those positions become Gly/Asn/Asp, which an inverse-folding model will choose itself; the
     score charges a small penalty per alpha-L residue so it is used only when it helps.
     Glycine is unrestricted. */
  const ALLOWED = {
    phi: [[-170, -45], [45, 90]],
    psi: [[-80, 180], [-180, -170]],
  };
  const ALPHA_L_PENALTY = 0.8;
  /* Glycine: anything but phi near zero. */
  const ALLOWED_GLY = { phi: [[-180, -45], [45, 180]], psi: [[-180, 180]] };

  function clampTo(v, intervals) {
    for (const [lo, hi] of intervals) if (v >= lo && v <= hi) return v;
    let best = v, bd = Infinity;
    for (const [lo, hi] of intervals) for (const b of [lo, hi]) {
      const dd = Math.abs(V.wrap(v - b));
      if (dd < bd) { bd = dd; best = b; }
    }
    return best;
  }

  function randomTorsion(rng, gly) {
    const g = () => { let s = 0; for (let i = 0; i < 4; i++) s += rng(); return (s - 2) * 1.7; };
    const r = rng();
    let phi, psi;
    if (r < (gly ? 0.15 : 0.06)) { phi = 70 + 10 * g(); psi = 25 + 15 * g(); }
    else if (r < 0.38) { phi = -65 + 12 * g(); psi = -40 + 12 * g(); }
    else if (r < 0.72) { phi = -120 + 20 * g(); psi = 130 + 20 * g(); }
    else { phi = -70 + 10 * g(); psi = 145 + 15 * g(); }
    const lim = gly ? ALLOWED_GLY : ALLOWED;
    phi = clampTo(V.wrap(phi), lim.phi); psi = clampTo(V.wrap(psi), lim.psi);
    return { phi, psi };
  }

  function mulberry(seed) {
    return function () {
      seed |= 0; seed = seed + 0x6D2B79F5 | 0;
      let t = Math.imul(seed ^ seed >>> 15, 1 | seed);
      t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
      return ((t ^ t >>> 14) >>> 0) / 4294967296;
    };
  }

  /* The atom list is a flat array P of positions; consecutive triples are residues (N CA C in the
     forward direction, C CA N in reverse). bonds: [{a, kind:'phi'|'psi', free}] rotate P[a+2..]
     about P[a] -> P[a+1]. targets: [{i, p, w}]. */
  function ccd(P, bonds, targets, opts) {
    const maxSweeps = opts.sweeps || 300, tol2 = (opts.tol || 0.02) ** 2;
    const freeRoot = !!opts.freeRoot;
    let prev = Infinity;
    for (let sweep = 0; sweep < maxSweeps; sweep++) {
      const order = sweep % 2 === 0 ? bonds : bonds.slice().reverse();
      for (const bd of order) {
        const a = bd.a, b = a + 1;
        const o = P[b], u = V.norm(V.sub(P[b], P[a]));
        let A = 0, B = 0;
        for (const t of targets) {
          if (t.i <= b) continue;
          const r = V.perp(V.sub(P[t.i], o), u);
          const f = V.sub(t.p, o);
          A += t.w * V.dot(f, r);
          B += t.w * V.dot(f, V.cross(u, r));
        }
        if (A === 0 && B === 0) continue;
        let th = Math.atan2(B, A) / V.D2R;
        const p0 = a > 0 ? P[a - 1] : opts.pre;
        if (p0) {
          const cur = V.dihedral(p0, P[a], P[b], P[b + 1]);
          const want = clampTo(V.wrap(cur + th), (bd.gly ? ALLOWED_GLY : ALLOWED)[bd.kind]);
          th = V.wrap(want - cur);
        }
        if (Math.abs(th) < 1e-6) continue;
        const R = V.rotAxis(u, th);
        for (let j = b + 1; j < P.length; j++) P[j] = V.add(o, V.m3vec(R, V.sub(P[j], o)));
      }
      if (freeRoot) {
        const { R, t } = V.kabsch(targets.map(x => P[x.i]), targets.map(x => x.p), targets.map(x => x.w));
        for (let j = 0; j < P.length; j++) P[j] = V.add(V.m3vec(R, P[j]), t);
      }
      const e = wmsd(P, targets);
      if (e < tol2 || Math.abs(prev - e) < 1e-9) break;
      prev = e;
    }
    return P;
  }

  function wmsd(P, targets) {
    let s = 0, w = 0;
    for (const t of targets) { s += t.w * V.dist2(P[t.i], t.p); w += t.w; }
    return w ? s / w : 0;
  }
  function rmsdOf(P, targets) {
    if (!targets.length) return 0;
    let s = 0; for (const t of targets) s += V.dist2(P[t.i], t.p);
    return Math.sqrt(s / targets.length);
  }

  /* Environment atoms for clash scoring: backbone + CB of everything except `skip` residue uids. */
  function envGrid(model, skip) {
    const cell = 4.0, grid = new Map();
    const add = p => {
      const k = Math.floor(p[0] / cell) + ',' + Math.floor(p[1] / cell) + ',' + Math.floor(p[2] / cell);
      if (!grid.has(k)) grid.set(k, []);
      grid.get(k).push(p);
    };
    for (const r of Model.allResidues(model)) {
      if (skip.has(r.uid)) continue;
      for (const a of r.atoms) if (r.het || a.name === 'N' || a.name === 'CA' || a.name === 'C' || a.name === 'O' || a.name === 'CB') add(a.p);
    }
    return {
      count(p, cut) {
        const gx = Math.floor(p[0] / cell), gy = Math.floor(p[1] / cell), gz = Math.floor(p[2] / cell);
        let n = 0; const c2 = cut * cut;
        for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
          const l = grid.get((gx + dx) + ',' + (gy + dy) + ',' + (gz + dz));
          if (l) for (const q of l) if (V.dist2(p, q) < c2) n++;
        }
        return n;
      },
    };
  }

  /* Clash count of loop backbone [{N, CA, C}] (with ideal CB unless glycine) against the
     environment and against itself. The cut-offs sit below backbone H-bond distances (N...O 2.9 A)
     so that tight turns are not rejected. ctx: { start: {N, CA, C} anchor, endN } adds the carbonyl
     O atoms, including the start anchor's, whose position the loop decides. */
  function loopClashes(bb, env, gly, ctx) {
    let n = 0;
    for (let k = 0; k < bb.length; k++) {
      const r = bb[k];
      for (const p of [r.N, r.CA, r.C]) n += env.count(p, 2.7);
      if (!gly) n += env.count(Geom.idealCB(r.N, r.CA, r.C), 2.9);
      for (let j = k + 2; j < bb.length; j++) if (V.dist(r.CA, bb[j].CA) < 4.0) n++;
    }
    const Os = [];
    if (ctx && ctx.start && bb.length) Os.push({ k: -1, p: Geom.carbonylO(ctx.start.N, ctx.start.CA, ctx.start.C, bb[0].N) });
    for (let k = 0; k < bb.length; k++) {
      const nN = k + 1 < bb.length ? bb[k + 1].N : ctx && ctx.endN;
      if (nN) Os.push({ k, p: Geom.carbonylO(bb[k].N, bb[k].CA, bb[k].C, nN) });
    }
    for (const o of Os) {
      n += env.count(o.p, 2.6);
      for (let j = 0; j < bb.length; j++) {
        if (Math.abs(j - o.k) < 2) continue;
        for (const p of [bb[j].N, bb[j].CA, bb[j].C]) if (V.dist2(o.p, p) < 2.6 * 2.6) n++;
      }
    }
    return n;
  }

  /* Core solver. spec:
       start: {N, CA, C, gly, prevC?} anchor residue the loop grows from (C-terminal end of a chain),
              or null. With prevC (the C of the residue before it) the anchor's own phi may change too;
              its new C is returned as startC.
       end:   {N, CA, C, gly} anchor residue the loop must join (N-terminal end), or null
       n:     loop residues
       caTargets: optional array of n points (or nulls) for the loop CAs; caWeights optional
       wClose: weight of the closure atoms
       gly:   loop residues are glycine
     Returns { bb: [{N, CA, C}] (loop residues, N->C), startC, closeRmsd, fitRmsd, clashes, score }. */
  function solveLoop(spec, env, restarts, seed) {
    const n = spec.n;
    const rng = mulberry(seed || 12345);
    const reverse = !spec.start && !!spec.end;
    const free = !spec.start && !spec.end;
    let best = null;
    for (let trial = 0; trial < restarts; trial++) {
      // Initial coordinates in residue order.
      let loop;
      if (spec.start) {
        loop = [];
        let { N, CA, C } = spec.start;
        let psiPrev = randomTorsion(rng, spec.start.gly).psi;
        for (let k = 0; k < n; k++) {
          const tk = randomTorsion(rng, spec.gly);
          const nx = Geom.extendForward(N, CA, C, psiPrev, tk.phi);
          loop.push(nx); ({ N, CA, C } = nx); psiPrev = tk.psi;
        }
        loop.psiLast = psiPrev;
      } else if (spec.end) {
        loop = [];
        let { N, CA, C } = spec.end;
        let phiCur = randomTorsion(rng, spec.end.gly).phi;
        for (let k = 0; k < n; k++) {
          const tk = randomTorsion(rng, spec.gly);
          const pv = Geom.extendBackward(N, CA, C, phiCur, tk.psi);
          loop.unshift(pv); ({ N, CA, C } = pv); phiCur = tk.phi;
        }
      } else {
        const tors = [];
        for (let k = 0; k < n; k++) { const t = randomTorsion(rng, spec.gly); tors.push({ phi: t.phi, psi: t.psi, omega: 180 }); }
        loop = Geom.backboneFromTorsions(tors);
      }

      // Flatten into the CCD atom list.
      const P = [], map = [];
      const bonds = [];
      const targets = [];
      if (!reverse) {
        let base = 0;
        if (spec.start) { P.push(spec.start.N, spec.start.CA, spec.start.C); base = 3; }
        loop.forEach((r, k) => { map[k] = { N: base + 3 * k, CA: base + 3 * k + 1, C: base + 3 * k + 2 }; P.push(r.N, r.CA, r.C); });
        if (spec.end) {
          // Virtual copy of the end anchor, continued from the last loop residue.
          const last = loop[n - 1];
          const ve = Geom.extendForward(last.N, last.CA, last.C, randomTorsion(rng, spec.gly).psi, randomTorsion(rng, spec.end.gly).phi);
          const i0 = P.length; P.push(ve.N, ve.CA, ve.C);
          targets.push({ i: i0, p: spec.end.N, w: spec.wClose }, { i: i0 + 1, p: spec.end.CA, w: spec.wClose }, { i: i0 + 2, p: spec.end.C, w: spec.wClose });
        }
        for (let a = 1; a + 2 < P.length; a++) {
          if (a % 3 === 2) continue;
          const resIdx = Math.floor(a / 3) - (spec.start ? 1 : 0);
          const kind = a % 3 === 0 ? 'phi' : 'psi';
          let gly;
          if (resIdx < 0) gly = spec.start.gly; else if (resIdx >= n) gly = spec.end.gly; else gly = spec.gly;
          bonds.push({ a, kind, gly });
        }
        if (spec.start && spec.start.prevC) bonds.unshift({ a: 0, kind: 'phi', gly: spec.start.gly });
      } else {
        P.push(spec.end.C, spec.end.CA, spec.end.N);
        for (let k = n - 1; k >= 0; k--) {
          const i0 = P.length; P.push(loop[k].C, loop[k].CA, loop[k].N);
          map[k] = { C: i0, CA: i0 + 1, N: i0 + 2 };
        }
        for (let a = 1; a + 2 < P.length; a++) {
          if (a % 3 === 2) continue;
          const resIdx = Math.floor(a / 3) - 1; // 0 = loop residue n-1
          const kind = a % 3 === 0 ? 'psi' : 'phi';
          bonds.push({ a, kind, gly: resIdx < 0 ? spec.end.gly : spec.gly });
        }
      }
      const caT = [];
      if (spec.caTargets) spec.caTargets.forEach((p, k) => { if (p) caT.push({ i: map[k].CA, p, w: spec.caWeights ? spec.caWeights[k] : 1 }); });
      const all = targets.concat(caT);

      if (free) {
        const { R, t } = V.kabsch(caT.map(x => P[x.i]), caT.map(x => x.p));
        for (let j = 0; j < P.length; j++) P[j] = V.add(V.m3vec(R, P[j]), t);
      }
      const pre = spec.start && !reverse ? spec.start.prevC : undefined;
      // With both a drawing and a closure, follow the drawing first and close afterwards: closing
      // first from a random start tends to lock in a conformation that ignores the drawing.
      if (targets.length && caT.length) ccd(P, bonds, caT, { sweeps: 150, pre });
      ccd(P, bonds, all, { freeRoot: free, sweeps: spec.sweeps || 250, pre });
      if (targets.length && caT.length) {
        ccd(P, bonds, targets.concat(caT.map(x => ({ i: x.i, p: x.p, w: 0.05 }))), { sweeps: 120, tol: 0.01, pre });
        if (rmsdOf(P, targets) > 0.05) ccd(P, bonds, targets, { sweeps: 400, tol: 0.01, pre });
      }

      const bb = loop.map((r, k) => ({ N: P[map[k].N], CA: P[map[k].CA], C: P[map[k].C] }));
      const closeRmsd = rmsdOf(P, targets);
      const fitRmsd = rmsdOf(P, caT);
      const clashes = env ? loopClashes(bb, env, spec.gly, { start: spec.start && !reverse ? { N: P[0], CA: P[1], C: P[2] } : null, endN: spec.end ? spec.end.N : null }) : 0;
      const alphaL = countAlphaL(bb, spec, spec.start && !reverse ? P[2] : null);
      const score = closeRmsd * 25 + fitRmsd * 2 + clashes * 1.5 + alphaL * ALPHA_L_PENALTY;
      if (!best || score < best.score) best = { bb, startC: spec.start && !reverse ? P[2] : null, closeRmsd, fitRmsd, clashes, alphaL, score };
      if (!caT.length && closeRmsd < 0.05 && clashes === 0 && trial >= 4) break;
    }
    return best;
  }

  /* Non-glycine residues with positive phi: loop residues, plus the end anchor (whose phi the loop
     decides). startC: the start anchor's C, which precedes the first loop residue. */
  function countAlphaL(bb, spec, startC) {
    let n = 0;
    if (!spec.gly) {
      for (let k = 0; k < bb.length; k++) {
        const prevC = k > 0 ? bb[k - 1].C : startC;
        if (prevC && V.dihedral(prevC, bb[k].N, bb[k].CA, bb[k].C) > 0) n++;
      }
    }
    if (spec.end && !spec.end.gly && bb.length && V.dihedral(bb[bb.length - 1].C, spec.end.N, spec.end.CA, spec.end.C) > 0) n++;
    return n;
  }

  function anchorOf(res) {
    return { N: Geom.pos(res, 'N'), CA: Geom.pos(res, 'CA'), C: Geom.pos(res, 'C'), gly: res.name === 'GLY' };
  }

  /* Close a gap between residue `from` (loop starts after it) and `to` (loop ends before it).
     n = loop length, or 0 for automatic. skip: uids excluded from clash scoring. */
  function closeLoop(model, from, to, n, opts) {
    opts = opts || {};
    const start = anchorOf(from), end = anchorOf(to);
    if (!start.N || !start.CA || !start.C || !end.N || !end.CA || !end.C) return { error: 'Both anchor residues need N, CA and C.' };
    const skip = new Set([from.uid, to.uid, ...(opts.skip || [])]);
    const env = envGrid(model, skip);
    const gly = opts.name === 'GLY';
    const d = V.dist(start.CA, end.CA);
    const nMin = Math.max(1, Math.ceil(d / 3.6) - 1);
    const lengths = n > 0 ? [n] : [nMin, nMin + 1, nMin + 2, nMin + 3, nMin + 4, nMin + 6];
    let best = null;
    for (const len of lengths) {
      const r = solveLoop({ start, end, n: len, wClose: 1, gly }, env, opts.restarts || 30, 1000 + len);
      r.n = len;
      if (!best || r.score + 0.15 * len < best.score + 0.15 * best.n) best = r;
      if (n <= 0 && r.closeRmsd < 0.15 && r.clashes === 0) { best = r; break; }
    }
    best.residues = residuesOf(best.bb, opts.name || 'ALA', end.N);
    best.gapCA = d;
    return best;
  }

  /* Fit a loop to a 3D polyline. spec: { points, start (residue or null), end (residue or null),
     spacing (A along the curve per residue), name }. */
  function fitCurve(model, spec) {
    let pts = spec.points.slice();
    const start = spec.start ? anchorOf(spec.start) : null;
    const end = spec.end ? anchorOf(spec.end) : null;
    if (start) pts[0] = start.CA;
    if (end) pts[pts.length - 1] = end.CA;
    pts = smoothCurve(resample(pts, 0.5), 3);
    if (start) pts[0] = start.CA;
    if (end) pts[pts.length - 1] = end.CA;
    const L = arcLength(pts);
    const spacing = spec.spacing || 3.2;
    // Residue count and arc-length position of each loop CA for a given count.
    const layout = n0 => {
      let n, sAt;
      if (start && end) { n = Math.max(1, n0 - 1); sAt = k => (k + 1) * L / (n + 1); }
      else if (start) { n = Math.max(1, n0); sAt = k => (k + 1) * L / n; }
      else if (end) { n = Math.max(1, n0); sAt = k => k * L / n; }
      else { n = Math.max(2, n0 + 1); sAt = k => k * L / (n - 1); }
      return { n, sAt };
    };
    const n0 = Math.round(L / spacing);
    if (n0 > 200) return { error: 'The drawing is too long (' + n0 + ' residues). Draw a shorter path.' };
    // The first steps out of an anchored end cannot turn sharply (the next CA sits within ~20 deg of
    // the anchor's CA->C bond), so targets right next to an anchor get less weight.
    const taper = (n, k) => {
      let w = 1;
      if (start) w = Math.min(w, 0.35 + 0.3 * k);
      if (end) w = Math.min(w, 0.35 + 0.3 * (n - 1 - k));
      return w;
    };
    const skip = new Set([spec.start && spec.start.uid, spec.end && spec.end.uid].filter(Boolean));
    const env = envGrid(model, skip);
    const gly = spec.name === 'GLY';
    if (start && spec.prevC) start.prevC = spec.prevC;
    let best = null;
    for (const cand of [n0, n0 - 1, n0 + 1]) {
      const { n, sAt } = layout(cand);
      if (best && best.n === n) continue;
      const caTargets = [], caWeights = [];
      for (let k = 0; k < n; k++) { caTargets.push(pointAt(pts, sAt(k))); caWeights.push(taper(n, k)); }
      const r = start && end && n >= 6
        ? fitBothEnds(start, end, caTargets, caWeights, env, gly, spec.restarts || 10, 777 + n)
        : solveLoop({ start, end, n, caTargets, caWeights, wClose: 10, gly, sweeps: 200 }, env, spec.restarts || 10, 777 + n);
      r.n = n;
      if (!best || r.score < best.score - 0.05) best = r;
    }
    best.residues = residuesOf(best.bb, spec.name || 'ALA', end ? end.N : null);
    best.length = L;
    return best;
  }

  /* A drawn loop joining two ends: fit the front part from the start anchor and the back part from
     the end anchor (each following the drawing), then close a short bridge of 3-4 residues between
     them. Closing the whole loop at once drags it away from the drawing. */
  function fitBothEnds(start, end, T, W, env, gly, restarts, seed) {
    let best = solveLoop({ start, end, n: T.length, caTargets: T, caWeights: W, wClose: 10, gly, sweeps: 200 }, env, restarts, seed);
    for (const m of [3, 4, 5]) {
      if (T.length - m < 2) break;
      const r = bridgeFit(start, end, T, W, env, gly, restarts, seed + 10 * m, m);
      if (r.score < best.score) best = r;
    }
    return best;
  }

  function bridgeFit(start, end, T, W, env, gly, restarts, seed, m) {
    const n = T.length;
    const f = Math.floor((n - m) / 2), b = f + m;
    const front = f > 0 ? solveLoop({ start, n: f, caTargets: T.slice(0, f), caWeights: W.slice(0, f), gly }, env, restarts, seed) : null;
    const back = solveLoop({ end, n: n - b, caTargets: T.slice(b), caWeights: W.slice(b), gly }, env, restarts, seed + 1);
    let midStart;
    if (front) {
      const last = front.bb[f - 1];
      const prevC = f > 1 ? front.bb[f - 2].C : (front.startC || start.C);
      midStart = { N: last.N, CA: last.CA, C: last.C, gly, prevC };
    } else midStart = start;
    const midEnd = { N: back.bb[0].N, CA: back.bb[0].CA, C: back.bb[0].C, gly };
    const mid = solveLoop({ start: midStart, end: midEnd, n: m, caTargets: T.slice(f, b), caWeights: W.slice(f, b).map(w => w * 0.3), wClose: 10, gly }, env, restarts * 2, seed + 2);
    const bb = [];
    if (front) {
      front.bb.forEach((r, k) => bb.push(k === f - 1 && mid.startC ? { N: r.N, CA: r.CA, C: mid.startC } : r));
    }
    bb.push(...mid.bb, ...back.bb);
    let s = 0;
    bb.forEach((r, k) => { s += V.dist2(r.CA, T[k]); });
    const fitRmsd = Math.sqrt(s / n);
    const clashes = env ? loopClashes(bb, env, gly, { start: { N: start.N, CA: start.CA, C: (front ? front.startC : mid.startC) || start.C }, endN: end.N }) : 0;
    const closeRmsd = mid.closeRmsd;
    const alphaL = (front ? front.alphaL : 0) + mid.alphaL + back.alphaL;
    return { bb, startC: front ? front.startC : mid.startC, closeRmsd, fitRmsd, clashes, alphaL, score: closeRmsd * 25 + fitRmsd * 2 + clashes * 1.5 + alphaL * ALPHA_L_PENALTY };
  }

  function residuesOf(bb, name, nNext) {
    const ext = nNext ? bb.concat([{ N: nNext, CA: null, C: null }]) : bb;
    const out = [];
    for (let i = 0; i < bb.length; i++) {
      const { N, CA, C } = bb[i];
      const O = Geom.carbonylO(N, CA, C, ext[i + 1] ? ext[i + 1].N : null, 135);
      const r = { name, num: i + 1, icode: '', ss: 'L', atoms: [{ name: 'N', el: 'N', p: N }, { name: 'CA', el: 'C', p: CA }, { name: 'C', el: 'C', p: C }, { name: 'O', el: 'O', p: O }] };
      Geom.buildSideChain(r, name);
      out.push(r);
    }
    return out;
  }

  /* ---------------- curve utilities ---------------- */
  function arcLength(pts) { let s = 0; for (let i = 1; i < pts.length; i++) s += V.dist(pts[i - 1], pts[i]); return s; }

  function pointAt(pts, s) {
    let acc = 0;
    for (let i = 1; i < pts.length; i++) {
      const l = V.dist(pts[i - 1], pts[i]);
      if (acc + l >= s) return V.lerp(pts[i - 1], pts[i], l > 0 ? (s - acc) / l : 0);
      acc += l;
    }
    return pts[pts.length - 1].slice();
  }

  function resample(pts, step) {
    const L = arcLength(pts);
    if (L < 1e-9) return [pts[0], pts[0]];
    const n = Math.max(2, Math.ceil(L / step) + 1);
    const out = [];
    for (let k = 0; k < n; k++) out.push(pointAt(pts, k * L / (n - 1)));
    return out;
  }

  /* Moving-average smoothing that keeps the end points. */
  function smoothCurve(pts, passes) {
    let p = pts;
    for (let it = 0; it < passes; it++) {
      const q = p.map(x => x.slice());
      for (let i = 1; i < p.length - 1; i++) q[i] = V.scale(V.add(V.add(p[i - 1], p[i + 1]), V.scale(p[i], 2)), 0.25);
      p = q;
    }
    return p;
  }

  window.Builders = {
    helix, sheet, riseOf, closeLoop, fitCurve, solveLoop, ccd, clampTo, ALLOWED, arcLength, resample, smoothCurve, pointAt,
  };
})();
