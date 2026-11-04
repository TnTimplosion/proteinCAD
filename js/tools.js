/* Viewport interaction: picking, camera, the transform gizmo, the CAD tools (helix, sheet, draw
   loop), the connect action, and the 2D overlay (plane grid, gizmo, strokes, labels, warnings). */
(function () {
  'use strict';
  const App = window.App, $ = App.$;
  let ov, ctx, dpr = 1;
  let drag = null;
  const cam = () => App.renderer.cam;

  const AXIS_COL = ['#e5534b', '#57ab5a', '#539bf5'];
  const BBSET = new Set(['N', 'CA', 'C', 'O', 'CB', 'OXT']);

  App.initViewport = function () {
    ov = $('overlay'); ctx = ov.getContext('2d');
    App.renderer = new Render.Renderer($('gl'));
    const resize = () => {
      const vp = $('viewport');
      dpr = window.devicePixelRatio || 1;
      const w = vp.clientWidth, h = vp.clientHeight;
      App.renderer.resize(w, h, dpr);
      ov.width = Math.round(w * dpr); ov.height = Math.round(h * dpr);
      App.cameraChanged();
    };
    resize();
    App.resizeViewport = resize;
    new ResizeObserver(resize).observe($('viewport'));
    window.addEventListener('resize', resize);
    ov.addEventListener('pointerdown', onDown);
    window.addEventListener('pointermove', onMove);
    window.addEventListener('pointerup', onUp);
    ov.addEventListener('wheel', onWheel, { passive: false });
    ov.addEventListener('dblclick', onDbl);
    ov.addEventListener('contextmenu', e => e.preventDefault());
    ov.addEventListener('pointerleave', () => { if (!drag && App.hover !== null) { App.hover = null; App.sceneDirty = true; App.refreshHud(); } App.mouse = null; App.overlayDirty = true; });
  };

  function mpos(e) { const r = ov.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; }

  /* ---------------- picking ---------------- */
  function projCache() {
    if (App.projCache) return App.projCache;
    const c = cam(), items = [];
    const mode = App.view.atoms;
    const add = (p, uid, worldR) => {
      const s = c.project(p); if (!s) return;
      items.push({ x: s[0], y: s[1], z: s[2], uid, r: worldR / c.pixelSize(p) });
    };
    for (const chain of App.model.chains) {
      if (App.hidden.has(chain.id)) continue;
      for (const res of chain.residues) {
        if (mode === 'none') { const ca = Geom.pos(res, 'CA'); if (ca) add(ca, res.uid, 1.6); continue; }
        for (const a of res.atoms) if (mode === 'all' || BBSET.has(a.name)) add(a.p, res.uid, a.name === 'CA' && App.view.cartoon ? 1.0 : 0.6);
      }
    }
    if (!App.hidden.has('lig')) for (const res of App.model.ligands) for (const a of res.atoms) add(a.p, res.uid, 0.6);
    App.projCache = items;
    return items;
  }

  function pick(mx, my) {
    let best = null;
    for (const it of projCache()) {
      const dx = it.x - mx, dy = it.y - my, rr = Math.max(it.r, 5);
      if (dx * dx + dy * dy < rr * rr && (!best || it.z < best.z)) best = it;
    }
    return best ? best.uid : null;
  }

  /* Chain terminus near the pointer: which = 'C' (last residue) or 'N' (first residue). */
  function terminalNear(mx, my, which) {
    let best = null, bd = 16;
    for (const chain of App.model.chains) {
      if (App.hidden.has(chain.id) || !chain.residues.length) continue;
      const res = which === 'C' ? chain.residues[chain.residues.length - 1] : chain.residues[0];
      const ca = Geom.pos(res, 'CA'); if (!ca) continue;
      const s = cam().project(ca); if (!s) continue;
      const d = Math.hypot(s[0] - mx, s[1] - my);
      if (d < bd) { bd = d; best = res; }
    }
    return best;
  }

  App.planePoint = function (mx, my, noSnap) {
    const r = cam().ray(mx, my), o = App.planeOrigin(), pl = App.plane;
    const p = V.rayPlane(r.o, r.d, o, pl.n);
    if (!p) return null;
    const s = App.opts.snap;
    if (noSnap || !(s > 0)) return p;
    const d = V.sub(p, o);
    const a = Math.round(V.dot(d, pl.u) / s) * s, b = Math.round(V.dot(d, pl.v) / s) * s;
    return V.addScaled(V.addScaled(o, pl.u, a), pl.v, b);
  };

  /* Point under the pointer at a given view depth (distance in front of the eye). */
  function screenPoint(mx, my, depth) {
    const c = cam(), r = c.ray(mx, my);
    if (c.ortho) return V.addScaled(r.o, r.d, depth);
    return V.addScaled(r.o, r.d, depth / V.dot(r.d, V.scale(c.back(), -1)));
  }

  /* ---------------- gizmo ---------------- */
  function gizmoFrame() {
    if (!App.gizmo.show || App.tool !== 'select' || !App.sel.size) return null;
    const pts = App.selectionPoints();
    if (!pts.length) return null;
    const c = V.centroid(pts);
    let axes = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    if (App.gizmo.space === 'local') {
      if (pts.length >= 3) {
        const p = V.pca(pts);
        const seg = Geom.segmentAxis(App.selectedResidues().filter(r => !r.het));
        const x = seg ? seg.u : p.axes[0];
        const y = V.norm(V.perp(p.axes[1], x));
        axes = [x, y, V.cross(x, y)];
      } else {
        const r = App.selectedResidues()[0];
        const N = Geom.pos(r, 'N'), CA = Geom.pos(r, 'CA'), C = Geom.pos(r, 'C');
        if (N && CA && C) { const f = Geom.residueFrame(N, CA, C); axes = [V.col(f.R, 0), V.col(f.R, 1), V.col(f.R, 2)]; }
      }
    }
    return { c, axes, L: 80 * cam().pixelSize(c) };
  }

  function ringPoints(c, ax, rad) {
    const a = V.anyPerp(ax), b = V.cross(ax, a), out = [];
    for (let i = 0; i <= 64; i++) {
      const t = i / 64 * Math.PI * 2;
      out.push(V.add(c, V.add(V.scale(a, Math.cos(t) * rad), V.scale(b, Math.sin(t) * rad))));
    }
    return out;
  }

  function distToSeg(px, py, ax, ay, bx, by) {
    const dx = bx - ax, dy = by - ay, l2 = dx * dx + dy * dy;
    const t = l2 ? Math.max(0, Math.min(1, ((px - ax) * dx + (py - ay) * dy) / l2)) : 0;
    return [Math.hypot(px - ax - t * dx, py - ay - t * dy), t];
  }

  function gizmoHit(mx, my) {
    const g = gizmoFrame(); if (!g) return null;
    const pc = cam().project(g.c); if (!pc) return null;
    if (Math.hypot(mx - pc[0], my - pc[1]) < 9) return { type: 'free', g };
    for (let k = 0; k < 3; k++) {
      const tp = cam().project(V.addScaled(g.c, g.axes[k], g.L)); if (!tp) continue;
      const [d, t] = distToSeg(mx, my, pc[0], pc[1], tp[0], tp[1]);
      if (d < 7 && t > 0.2) return { type: 'move', k, g };
    }
    let best = null, bd = 6;
    for (let k = 0; k < 3; k++) {
      const pts = ringPoints(g.c, g.axes[k], g.L * 0.72).map(p => cam().project(p));
      for (let i = 1; i < pts.length; i++) {
        if (!pts[i] || !pts[i - 1]) continue;
        const [d] = distToSeg(mx, my, pts[i - 1][0], pts[i - 1][1], pts[i][0], pts[i][1]);
        if (d < bd) { bd = d; best = { type: 'rot', k, g }; }
      }
    }
    return best;
  }

  function startGizmo(hit, mx, my) {
    const atoms = [];
    for (const r of App.selectedResidues()) for (const a of r.atoms) atoms.push({ a, p0: a.p.slice() });
    const g = { c: hit.g.c.slice(), axes: hit.g.axes, L: hit.g.L };
    const r0 = cam().ray(mx, my);
    drag = { kind: 'gizmo', hit, g, atoms, snap: App.snapshot(), x0: mx, y0: my, moved: false, delta: [0, 0, 0] };
    if (hit.type === 'move') drag.s0 = V.lineRayParam(g.c, g.axes[hit.k], r0.o, r0.d);
    if (hit.type === 'free') drag.h0 = V.rayPlane(r0.o, r0.d, g.c, cam().back());
    if (hit.type === 'rot') {
      const ax = g.axes[hit.k];
      drag.h0 = Math.abs(V.dot(r0.d, ax)) > 0.2 ? V.rayPlane(r0.o, r0.d, g.c, ax) : null;
    }
  }

  function moveGizmo(mx, my, e) {
    const d = drag, r = cam().ray(mx, my), snap = e.ctrlKey || e.metaKey;
    const c = d.g.c;
    let T = null;
    if (d.hit.type === 'move') {
      const ax = d.g.axes[d.hit.k];
      const s = V.lineRayParam(c, ax, r.o, r.d);
      if (s === null || d.s0 === null) return;
      let delta = s - d.s0;
      if (snap) delta = Math.round(delta / App.gizmo.snapMove) * App.gizmo.snapMove;
      d.delta = V.scale(ax, delta);
      d.readout = (delta >= 0 ? '+' : '') + delta.toFixed(2) + ' Å';
      T = p => V.add(p, d.delta);
    } else if (d.hit.type === 'free') {
      const h = V.rayPlane(r.o, r.d, c, cam().back());
      if (!h || !d.h0) return;
      let delta = V.sub(h, d.h0);
      if (snap) delta = delta.map(x => Math.round(x / App.gizmo.snapMove) * App.gizmo.snapMove);
      d.delta = delta;
      d.readout = V.len(delta).toFixed(2) + ' Å';
      T = p => V.add(p, delta);
    } else {
      const ax = d.g.axes[d.hit.k];
      let ang;
      if (d.h0) {
        const h = V.rayPlane(r.o, r.d, c, ax);
        if (!h) return;
        const v0 = V.perp(V.sub(d.h0, c), ax), v1 = V.perp(V.sub(h, c), ax);
        ang = Math.atan2(V.dot(V.cross(v0, v1), ax), V.dot(v0, v1)) / V.D2R;
      } else {
        ang = ((mx - d.x0) - (my - d.y0)) * 0.5;
      }
      if (snap) ang = Math.round(ang / App.gizmo.snapRot) * App.gizmo.snapRot;
      d.readout = (ang >= 0 ? '+' : '') + ang.toFixed(1) + '°';
      const R = V.rotAxis(ax, ang);
      T = p => V.add(c, V.m3vec(R, V.sub(p, c)));
    }
    for (const it of d.atoms) it.a.p = T(it.p0);
    d.moved = true;
    App.moved();
  }

  /* ---------------- pointer handlers ---------------- */
  function onDown(e) {
    if (!cam().w || !cam().h) App.resizeViewport();
    const [mx, my] = mpos(e);
    try { ov.setPointerCapture(e.pointerId); } catch (_) { /* ignore */ }
    if (e.button === 1 || (e.button === 2 && e.shiftKey)) { drag = { kind: 'pan', x: mx, y: my }; e.preventDefault(); return; }
    if (e.button === 2) { drag = { kind: 'orbit', x: mx, y: my }; return; }
    if (e.button !== 0) return;
    const gh = gizmoHit(mx, my);
    if (gh) { startGizmo(gh, mx, my); return; }
    const t = App.tool;
    if (t === 'select') {
      const uid = pick(mx, my);
      if (e.shiftKey && uid === null) { drag = { kind: 'box', x0: mx, y0: my, x: mx, y: my, remove: e.ctrlKey || e.metaKey }; return; }
      drag = { kind: 'click', x0: mx, y0: my, uid, shift: e.shiftKey, ctrl: e.ctrlKey || e.metaKey };
      return;
    }
    if (t === 'helix' || t === 'sheet') {
      const p = App.planePoint(mx, my);
      if (!p) { drag = { kind: 'orbit', x: mx, y: my }; App.status('The sketch plane is edge-on or behind the camera here. Orbit or pick another plane.', 'err'); return; }
      drag = t === 'helix' ? { kind: 'helix', base: p, x0: mx, y0: my, n: 0, moved: false } : { kind: 'sheet', p0: p, p1: p, x0: mx, y0: my, moved: false };
      return;
    }
    if (t === 'draw') {
      drag = { kind: 'draw', pts: [[mx, my]], startC: terminalNear(mx, my, 'C'), startN: terminalNear(mx, my, 'N') };
    }
  }

  function onMove(e) {
    const [mx, my] = mpos(e);
    if (!drag) {
      if (e.target !== ov) return;
      App.mouse = [mx, my];
      const uid = App.tool === 'select' || App.tool === 'draw' ? pick(mx, my) : null;
      if (uid !== App.hover) { App.hover = uid; App.sceneDirty = true; App.refreshHud(); App.highlightSeq && App.highlightSeq(); }
      App.overlayDirty = true;
      return;
    }
    App.mouse = [mx, my];
    switch (drag.kind) {
      case 'orbit': cam().orbit(mx - drag.x, my - drag.y); drag.x = mx; drag.y = my; App.cameraChanged(); break;
      case 'pan': cam().pan(mx - drag.x, my - drag.y); drag.x = mx; drag.y = my; App.cameraChanged(); break;
      case 'click':
        if (Math.hypot(mx - drag.x0, my - drag.y0) > 4) { drag = { kind: 'orbit', x: mx, y: my }; }
        break;
      case 'box': drag.x = mx; drag.y = my; App.overlayDirty = true; break;
      case 'gizmo': moveGizmo(mx, my, e); break;
      case 'helix': updateHelix(mx, my); break;
      case 'sheet': updateSheet(mx, my); break;
      case 'draw': {
        const last = drag.pts[drag.pts.length - 1];
        if (Math.hypot(mx - last[0], my - last[1]) > 2.5) { drag.pts.push([mx, my]); App.overlayDirty = true; }
        break;
      }
    }
  }

  function onUp(e) {
    if (!drag) return;
    const d = drag; drag = null;
    const [mx, my] = mpos(e);
    switch (d.kind) {
      case 'click': finishClick(d); break;
      case 'box': finishBox(d); break;
      case 'gizmo':
        if (d.moved) {
          App.pushUndo(d.snap, d.hit.type === 'rot' ? 'Rotate' : 'Move');
          App.changed();
          App.status((d.hit.type === 'rot' ? 'Rotated ' : 'Moved ') + App.sel.size + ' residue(s) ' + (d.readout || '') + '.');
        }
        break;
      case 'helix': commitHelix(d); break;
      case 'sheet': commitSheet(d); break;
      case 'draw': d.pts.push([mx, my]); commitDraw(d, mx, my); break;
    }
    App.ghost = null; App.sceneDirty = true; App.overlayDirty = true;
  }

  function onWheel(e) {
    e.preventDefault();
    cam().zoom(Math.exp(e.deltaY * (e.deltaMode === 1 ? 0.05 : 0.0015)));
    App.cameraChanged();
  }

  function onDbl(e) {
    if (App.tool !== 'select') return;
    const [mx, my] = mpos(e);
    const uid = pick(mx, my);
    if (uid !== null) App.setSelection(App.expand(uid, 'element'), 'set');
  }

  function finishClick(d) {
    if (d.uid === null) { if (!d.shift && !d.ctrl) App.setSelection([], 'set'); return; }
    const units = App.expand(d.uid);
    if (d.ctrl) {
      const all = units.every(u => App.sel.has(u));
      App.setSelection(units, all ? 'remove' : 'add');
    } else if (d.shift && App.lastClick !== undefined) {
      const a = App.index.get(App.lastClick), b = App.index.get(d.uid);
      if (a && b && a.chain && a.chain === b.chain) {
        const lo = Math.min(a.idx, b.idx), hi = Math.max(a.idx, b.idx);
        App.setSelection(a.chain.residues.slice(lo, hi + 1).map(r => r.uid), 'add');
      } else App.setSelection(units, 'add');
    } else {
      App.setSelection(units, 'set');
    }
    App.lastClick = d.uid;
  }

  function finishBox(d) {
    const x0 = Math.min(d.x0, d.x), x1 = Math.max(d.x0, d.x), y0 = Math.min(d.y0, d.y), y1 = Math.max(d.y0, d.y);
    if (x1 - x0 < 3 && y1 - y0 < 3) return;
    const uids = [];
    const test = (res) => {
      const p = (!res.het && Geom.pos(res, 'CA')) || (res.atoms[0] && res.atoms[0].p);
      if (!p) return;
      const s = cam().project(p);
      if (s && s[0] >= x0 && s[0] <= x1 && s[1] >= y0 && s[1] <= y1) uids.push(res.uid);
    };
    for (const c of App.model.chains) if (!App.hidden.has(c.id)) c.residues.forEach(test);
    if (!App.hidden.has('lig')) App.model.ligands.forEach(test);
    App.setSelection(uids, d.remove ? 'remove' : 'add');
    App.status((d.remove ? 'Removed ' : 'Added ') + uids.length + ' residue(s) ' + (d.remove ? 'from' : 'to') + ' the selection.');
  }

  /* ---------------- helix tool ---------------- */
  function helixConf() { return Geom.CONF[App.opts.helixConf] || Geom.CONF.alpha; }

  function helixFrom(base, nres, sign) {
    return Builders.helix({
      n: nres, conf: helixConf(), base, dir: V.scale(App.plane.n, sign), inPlane: App.plane.u,
      phase: App.opts.helixPhase, endAtPlane: App.opts.helixEnd === 'end', name: App.opts.newRes,
    });
  }

  function updateHelix(mx, my) {
    const d = drag;
    if (Math.hypot(mx - d.x0, my - d.y0) > 4) d.moved = true;
    if (!d.moved) return;
    const r = cam().ray(mx, my);
    let h = V.lineRayParam(d.base, App.plane.n, r.o, r.d);
    if (h === null) h = (d.y0 - my) * cam().pixelSize(d.base);
    const rise = Builders.riseOf(helixConf());
    d.n = Math.max(1, Math.round(Math.abs(h) / rise) + 1);
    d.sign = h >= 0 ? 1 : -1;
    d.label = helixConf().label + ' · ' + d.n + ' res · ' + ((d.n - 1) * rise).toFixed(1) + ' Å' + (d.sign < 0 ? ' (below plane)' : '');
    App.ghost = d.n >= 2 ? [helixFrom(d.base, d.n, d.sign)] : null;
    App.sceneDirty = true; App.overlayDirty = true;
  }

  function commitHelix(d) {
    let n = d.n, sign = d.sign || 1;
    if (!d.moved) { n = App.opts.helixLen; sign = 1; }
    if (n < 2) { App.status('Drag further to make a helix of at least 2 residues.'); return; }
    let chain;
    App.edit('Extrude helix', () => { chain = Model.addChain(App.model, helixFrom(d.base, n, sign)); });
    App.setSelection(chain.residues.map(r => r.uid), 'set');
    App.status('Helix ' + chain.id + ': ' + n + ' residues (' + helixConf().label + ').', 'ok');
  }

  /* ---------------- sheet tool ---------------- */
  function sheetFrom(p0, p1) {
    const o = App.opts;
    return Builders.sheet({
      p0, p1, normal: App.plane.n, nStrands: o.sheetStrands, parallel: o.sheetParallel, spacing: o.sheetSpacing,
      side: o.sheetSide, twist: o.sheetTwist, shift: o.sheetShift, name: App.opts.newRes,
    });
  }
  function updateSheet(mx, my) {
    const d = drag;
    const p = App.planePoint(mx, my);
    if (!p) return;
    d.p1 = p;
    const L = V.dist(d.p0, d.p1);
    if (L > 2) d.moved = true;
    if (!d.moved) return;
    const conf = App.opts.sheetParallel ? Geom.CONF.parallel : Geom.CONF.antiparallel;
    const rise = Builders.riseOf(conf);
    const n = Math.max(2, Math.round(L / rise) + 1);
    d.label = App.opts.sheetStrands + ' strands × ' + n + ' res · ' + (App.opts.sheetParallel ? 'parallel' : 'antiparallel');
    App.ghost = sheetFrom(d.p0, d.p1);
    App.sceneDirty = true; App.overlayDirty = true;
  }
  function commitSheet(d) {
    let p1 = d.p1;
    if (!d.moved) {
      const conf = App.opts.sheetParallel ? Geom.CONF.parallel : Geom.CONF.antiparallel;
      p1 = V.addScaled(d.p0, App.plane.u, (App.opts.sheetLen - 1) * Builders.riseOf(conf));
    }
    const created = [];
    App.edit('Sheet', () => { for (const strand of sheetFrom(d.p0, p1)) created.push(Model.addChain(App.model, strand)); });
    App.setSelection(created.flatMap(c => c.residues.map(r => r.uid)), 'set');
    App.status('Sheet of ' + created.length + ' strands (chains ' + created.map(c => c.id).join(', ') + '). Connect neighbouring ends with ⌒ Connect to make hairpins.', 'ok');
  }

  /* ---------------- draw-loop tool ---------------- */
  function commitDraw(d, mx, my) {
    if (d.pts.length < 3) return;
    let pts2 = d.pts.slice();
    let start = d.startC, end = terminalNear(mx, my, 'N');
    // A stroke that starts on an N-terminus runs backwards into it.
    if (!start && d.startN && !end) {
      pts2.reverse();
      end = d.startN;
      start = terminalNear(mx, my, 'C');
    }
    if (start && end && App.index.get(start.uid).chain === App.index.get(end.uid).chain) {
      App.status('Both ends of the drawing are on the same chain; a loop cannot join a chain to itself.', 'err');
      return;
    }
    const surface = App.opts.drawSurface === 'auto' ? (start || end ? 'screen' : 'plane') : App.opts.drawSurface;
    let pts3;
    if (surface === 'plane') {
      pts3 = pts2.map(([x, y]) => App.planePoint(x, y, true)).filter(Boolean);
      if (pts3.length < 2) { App.status('The drawing missed the sketch plane. Orbit so the plane faces you, or set "Draw on" to Screen.', 'err'); return; }
    } else {
      const c = cam();
      const depthOf = p => -c.toView(p)[2];
      const z0 = start ? depthOf(Geom.pos(start, 'CA')) : end ? depthOf(Geom.pos(end, 'CA')) : depthOf(c.target);
      const z1 = end ? depthOf(Geom.pos(end, 'CA')) : z0;
      let total = 0; const acc = [0];
      for (let i = 1; i < pts2.length; i++) { total += Math.hypot(pts2[i][0] - pts2[i - 1][0], pts2[i][1] - pts2[i - 1][1]); acc.push(total); }
      pts3 = pts2.map(([x, y], i) => screenPoint(x, y, z0 + (z1 - z0) * (total ? acc[i] / total : 0)));
    }
    // The C of the residue before the start anchor lets the anchor's own phi bend toward the drawing.
    const prevC = (() => {
      if (!start) return null;
      const i = App.index.get(start.uid), p = i.chain.residues[i.idx - 1];
      return p && Model.bonded(p, start) ? Geom.pos(p, 'C') : null;
    })();
    App.busy(true);
    setTimeout(() => {
      let res;
      try { res = Builders.fitCurve(App.model, { points: pts3, start, end, prevC, spacing: App.opts.drawSpacing, name: App.opts.loopRes }); }
      finally { App.busy(false); }
      if (res.error) { App.status(res.error, 'err'); return; }
      let target;
      App.edit('Draw loop', () => {
        const loop = res.residues.map(r => Model.stamp(App.model, r));
        if (start && prevC && res.startC) {
          const A = App.index.get(start.uid).chain;
          Model.setTorsion(A, A.residues.length - 1, 'phi', V.dihedral(prevC, Geom.pos(start, 'N'), Geom.pos(start, 'CA'), res.startC), 'C');
        }
        if (start && end) {
          const A = App.index.get(start.uid).chain, B = App.index.get(end.uid).chain;
          const ia = A.residues.length - 1;
          loop.forEach((r, k) => { r.num = start.num + 1 + k; });
          A.residues.push(...loop);
          Model.mergeChains(App.model, A, B);
          Model.rebuildO(A, ia);
          if (App.opts.renumber) Model.renumber(A.residues, A.residues[0].num);
          target = A;
        } else if (start) {
          const A = App.index.get(start.uid).chain;
          const ia = A.residues.length - 1;
          loop.forEach((r, k) => { r.num = start.num + 1 + k; });
          A.residues.push(...loop);
          Model.rebuildO(A, ia);
          target = A;
        } else if (end) {
          const B = App.index.get(end.uid).chain;
          loop.forEach((r, k) => { r.num = end.num - loop.length + k; });
          B.residues.unshift(...loop);
          if (App.opts.renumber && B.residues[0].num < 1) Model.renumber(B.residues, 1);
          target = B;
        } else {
          target = Model.addChain(App.model, loop);
        }
        App._newLoop = loop.map(r => r.uid);
      });
      App.setSelection(App._newLoop, 'set');
      const parts = [res.n + '-residue loop in chain ' + target.id, 'fit to drawing ' + res.fitRmsd.toFixed(2) + ' Å CA RMSD'];
      if (end) parts.push('closure ' + res.closeRmsd.toFixed(2) + ' Å');
      if (res.clashes) parts.push(res.clashes + ' close contact(s)');
      App.status(parts.join(', ') + '.', res.closeRmsd > 0.5 ? 'err' : 'ok');
    }, 20);
  }

  App.busy = function (on) { $('busy').hidden = !on; };

  /* ---------------- connect ---------------- */
  /* Which pair (from, to) of residues can be joined: the loop goes from `from`'s C-terminus to
     `to`'s N-terminus. */
  function connectPair(r1, r2) {
    const a = App.index.get(r1.uid), b = App.index.get(r2.uid);
    if (!a || !b || !a.chain || !b.chain) return null;
    if (a.chain === b.chain) {
      const [lo, hi] = a.idx < b.idx ? [a, b] : [b, a];
      if (hi.idx === lo.idx) return null;
      return { from: lo.res, to: hi.res, chain: lo.chain, same: true, lo: lo.idx, hi: hi.idx };
    }
    const isLast = x => x.idx === x.chain.residues.length - 1, isFirst = x => x.idx === 0;
    if (isLast(a) && isFirst(b)) return { from: a.res, to: b.res };
    if (isLast(b) && isFirst(a)) return { from: b.res, to: a.res };
    return null;
  }

  App.connect = function () {
    const sel = App.selectedResidues().filter(r => !r.het);
    const chains = App.selectedChains(true);
    const fullChains = chains.filter(c => c.residues.every(r => App.sel.has(r.uid)));
    let pairs = null;
    if (fullChains.length >= 2 && fullChains.length === chains.length) {
      pairs = [];
      for (let i = 0; i + 1 < fullChains.length; i++) pairs.push({ chainA: fullChains[i], chainB: fullChains[i + 1] });
    } else if (sel.length === 2) {
      const p = connectPair(sel[0], sel[1]);
      if (!p) {
        App.status('Connect needs the last residue of one chain and the first residue of another, two residues of the same chain (to rebuild the stretch between them), or two or more whole chains (joined in the order you selected them).', 'err');
        return;
      }
      pairs = [p];
    } else {
      App.status('Select two chain ends (C-terminus then N-terminus) or whole chains in the order to join them, then Connect.', 'err');
      return;
    }
    App.busy(true);
    setTimeout(() => {
      const snap = App.snapshot();
      const reports = [];
      let lastChain = null;
      try {
        for (const p of pairs) {
          let from, to, same = false, lo, hi;
          if (p.chainA) {
            // An earlier pair may have merged chainA into another chain; its residues are still
            // the same objects, so its last residue is still the end to grow from.
            from = p.chainA.residues[p.chainA.residues.length - 1]; to = p.chainB.residues[0];
          } else ({ from, to, same, lo, hi } = p);
          const skip = same ? App.index.get(from.uid).chain.residues.slice(lo + 1, hi).map(r => r.uid) : [];
          const res = Builders.closeLoop(App.model, from, to, App.opts.loopLen, { name: App.opts.loopRes, skip, restarts: 30 });
          if (res.error) { reports.push(res.error); continue; }
          const loop = res.residues.map(r => Model.stamp(App.model, r));
          const A = App.index.get(from.uid).chain;
          const ia = A.residues.indexOf(from);
          if (same) {
            A.residues.splice(ia + 1, A.residues.indexOf(to) - ia - 1, ...loop);
          } else {
            const B = App.index.get(to.uid).chain;
            A.residues.push(...loop);
            Model.mergeChains(App.model, A, B);
          }
          loop.forEach((r, k) => { r.num = from.num + 1 + k; });
          Model.rebuildO(A, ia);
          if (App.opts.renumber) Model.renumber(A.residues, A.residues[0].num);
          App.index = Model.buildIndex(App.model);
          lastChain = A;
          reports.push(loop.length + ' res, closure ' + res.closeRmsd.toFixed(2) + ' Å' + (res.clashes ? ', ' + res.clashes + ' close contact(s)' : '') + (res.closeRmsd > 0.5 ? ' (NOT closed: try more residues)' : ''));
        }
      } catch (e) {
        App.model = snap; App.changed(); App.busy(false);
        App.status('Connect failed: ' + e.message, 'err'); console.error(e);
        return;
      }
      App.busy(false);
      App.pushUndo(snap, 'Connect');
      App.changed();
      if (lastChain) App.setSelection(lastChain.residues.map(r => r.uid), 'set');
      const bad = reports.some(r => /NOT|need/.test(r));
      App.status('Connected' + (lastChain ? ' into chain ' + lastChain.id : '') + ': ' + reports.join('; ') + '.', bad ? 'err' : 'ok');
    }, 20);
  };

  /* ---------------- refine ---------------- */
  App.refineOpts = function () {
    const o = App.opts;
    return { hold: o.refineHold, clashScale: o.refineClash, maxGap: o.refineGap, iterations: o.refineIters, sideMode: o.refineSide };
  };

  /* Residues to regularise: the selection, or whole chains when `chains` is given. */
  App.refine = function (chains) {
    const residues = chains ? chains.flatMap(c => c.residues) : App.selectedResidues().filter(r => !r.het);
    if (!residues.length) {
      App.status('Select the residues to refine first (or use Refine in Chain tools for whole chains).', 'err');
      return;
    }
    App.busy(true);
    setTimeout(() => {
      const snap = App.snapshot();
      let rep;
      try { rep = Refine.run(App.model, residues, App.refineOpts()); }
      catch (e) {
        App.model = snap; App.changed(); App.busy(false);
        App.status('Refine failed: ' + e.message, 'err'); console.error(e); return;
      }
      App.busy(false);
      if (!rep) { App.status('Nothing to refine in the selection.', 'err'); return; }
      App.pushUndo(snap, 'Refine');
      App.changed();
      App.lastRefine = rep;
      App.status(refineReport(rep), rep.skippedBreaks.length ? 'err' : 'ok');
    }, 20);
  };

  function refineReport(rep) {
    const b = rep.before, a = rep.after;
    const parts = ['Refined ' + rep.residues + ' residue(s)'];
    parts.push('bonds ' + b.bondRms.toFixed(3) + ' → ' + a.bondRms.toFixed(3) + ' Å RMS');
    parts.push('angles ' + b.angleRms.toFixed(3) + ' → ' + a.angleRms.toFixed(3) + ' Å');
    if (b.clashes || a.clashes) parts.push('clashes ' + b.clashes + ' → ' + a.clashes);
    if (b.rama.length || a.rama.length) parts.push('φ outliers ' + b.rama.length + ' → ' + a.rama.length);
    if (rep.pack) {
      parts.push('repacked ' + rep.pack.changed + ' of ' + rep.pack.sites + ' side chain(s)' +
        (rep.pack.committed ? ' (' + rep.pack.committed + ' held by a disulfide or ligand)' : ''));
    }
    if (b.wrongChirality || a.wrongChirality) parts.push('D-centres ' + b.wrongChirality + ' → ' + a.wrongChirality);
    for (const c of rep.closed) parts.push('break after ' + c.from.num + ' closed (' + c.was.toFixed(2) + ' → ' + c.now.toFixed(2) + ' Å)');
    parts.push('moved ' + rep.rmsd.toFixed(2) + ' Å RMS, up to ' + rep.maxShift.toFixed(2) + ' Å');
    let s = parts.join(', ') + '.';
    if (rep.after.rama.length) {
      s += ' φ still outside the allowed range at ' +
        rep.after.rama.slice(0, 6).map(x => x.res.name + x.res.num + ' (' + x.phi.toFixed(0) + '°)').join(', ') + '.';
    }
    if (rep.skippedBreaks.length) {
      s += ' ' + rep.skippedBreaks.length + ' gap(s) too wide to close (' +
        rep.skippedBreaks.map(x => x.gap.toFixed(1)).join(', ') + ' Å) — use ⌒ Connect to build a loop there.';
    }
    return s;
  }
  App.refineReport = refineReport;

  /* ---------------- other actions ---------------- */
  App.deleteSelection = function () {
    if (!App.sel.size) return;
    const n = App.sel.size, uids = new Set(App.sel);
    App.edit('Delete', () => Model.deleteResidues(App.model, uids));
    App.status('Deleted ' + n + ' residue(s).');
  };
  App.duplicateSelection = function () {
    const rs = App.selectedResidues().filter(r => !r.het);
    if (!rs.length) return;
    let c;
    App.edit('Duplicate', () => { c = Model.duplicate(App.model, rs); });
    App.setSelection(c.residues.map(r => r.uid), 'set');
    App.status('Copied ' + rs.length + ' residue(s) into new chain ' + c.id + '. Move it with the gizmo.', 'ok');
  };

  /* ---------------- overlay ---------------- */
  function seg3(a, b) {
    const c = cam();
    const zc = -c.nearFar()[0] * 0.6;
    if (!c.ortho) {
      const za = c.toView(a)[2], zb = c.toView(b)[2];
      if (za > zc && zb > zc) return null;
      if (za > zc) a = V.lerp(a, b, (za - zc) / (za - zb));
      else if (zb > zc) b = V.lerp(b, a, (zb - zc) / (zb - za));
    }
    const pa = c.project(a), pb = c.project(b);
    return pa && pb ? [pa[0], pa[1], pb[0], pb[1]] : null;
  }
  function line3(a, b) { const s = seg3(a, b); if (s) { ctx.moveTo(s[0], s[1]); ctx.lineTo(s[2], s[3]); } return s; }
  function poly3(pts) { for (let i = 1; i < pts.length; i++) line3(pts[i - 1], pts[i]); }
  function text(t, x, y, col, align) {
    ctx.font = '12px system-ui, sans-serif';
    ctx.textAlign = align || 'left';
    ctx.lineWidth = 3; ctx.strokeStyle = 'rgba(0,0,0,0.75)'; ctx.strokeText(t, x, y);
    ctx.fillStyle = col || '#fff'; ctx.fillText(t, x, y);
  }

  function drawPlane() {
    const pl = App.plane, o = App.planeOrigin();
    const b = Model.bounds(App.model);
    let H = 30;
    if (b) H = Math.min(150, Math.max(30, Math.ceil((V.dist(b.center, o) + b.radius) / 10) * 10));
    const step = H > 80 ? 10 : 5;
    ctx.lineWidth = 1;
    ctx.strokeStyle = 'rgba(120,170,230,0.16)';
    ctx.beginPath();
    for (let k = -H; k <= H; k += step) {
      if (k === 0) continue;
      line3(V.addScaled(V.addScaled(o, pl.u, k), pl.v, -H), V.addScaled(V.addScaled(o, pl.u, k), pl.v, H));
      line3(V.addScaled(V.addScaled(o, pl.v, k), pl.u, -H), V.addScaled(V.addScaled(o, pl.v, k), pl.u, H));
    }
    ctx.stroke();
    ctx.strokeStyle = 'rgba(120,170,230,0.45)';
    ctx.beginPath();
    poly3([V.addScaled(V.addScaled(o, pl.u, -H), pl.v, -H), V.addScaled(V.addScaled(o, pl.u, H), pl.v, -H), V.addScaled(V.addScaled(o, pl.u, H), pl.v, H), V.addScaled(V.addScaled(o, pl.u, -H), pl.v, H), V.addScaled(V.addScaled(o, pl.u, -H), pl.v, -H)]);
    ctx.stroke();
    ctx.lineWidth = 1.5;
    ctx.strokeStyle = 'rgba(229,83,75,0.7)'; ctx.beginPath(); line3(o, V.addScaled(o, pl.u, H)); ctx.stroke();
    ctx.strokeStyle = 'rgba(87,171,90,0.7)'; ctx.beginPath(); line3(o, V.addScaled(o, pl.v, H)); ctx.stroke();
    ctx.lineWidth = 2;
    ctx.strokeStyle = 'rgba(83,155,245,0.9)'; ctx.beginPath(); const s = line3(o, V.addScaled(o, pl.n, 8)); ctx.stroke();
    if (s) text('n', s[2] + 4, s[3] - 4, '#8fbaf8');
  }

  function drawFootprint(p) {
    const pl = App.plane;
    for (const [rad, col] of [[5.0, 'rgba(255,255,255,0.55)'], [2.3, 'rgba(255,255,255,0.3)']]) {
      const pts = [];
      for (let i = 0; i <= 48; i++) { const t = i / 48 * Math.PI * 2; pts.push(V.add(p, V.add(V.scale(pl.u, Math.cos(t) * rad), V.scale(pl.v, Math.sin(t) * rad)))); }
      ctx.strokeStyle = col; ctx.lineWidth = 1.5; ctx.beginPath(); poly3(pts); ctx.stroke();
    }
    const s = cam().project(p);
    if (s) { ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(s[0], s[1], 2.5, 0, 7); ctx.fill(); }
  }

  function drawGizmo() {
    const live = drag && drag.kind === 'gizmo';
    const g = live ? { c: drag.hit.type === 'rot' ? drag.g.c : V.add(drag.g.c, drag.delta), axes: drag.g.axes, L: drag.g.L } : gizmoFrame();
    if (!g) return;
    const pc = cam().project(g.c); if (!pc) return;
    const hot = !live && App.mouse ? gizmoHit(App.mouse[0], App.mouse[1]) : (live ? drag.hit : null);
    const isHot = (type, k) => hot && hot.type === type && (type === 'free' || hot.k === k);
    for (let k = 0; k < 3; k++) {
      const pts = ringPoints(g.c, g.axes[k], g.L * 0.72);
      ctx.strokeStyle = isHot('rot', k) ? '#ffd02e' : AXIS_COL[k];
      ctx.globalAlpha = isHot('rot', k) ? 1 : 0.55;
      ctx.lineWidth = isHot('rot', k) ? 3 : 1.6;
      ctx.beginPath(); poly3(pts); ctx.stroke();
    }
    ctx.globalAlpha = 1;
    for (let k = 0; k < 3; k++) {
      const tip = cam().project(V.addScaled(g.c, g.axes[k], g.L)); if (!tip) continue;
      const col = isHot('move', k) ? '#ffd02e' : AXIS_COL[k];
      ctx.strokeStyle = col; ctx.fillStyle = col; ctx.lineWidth = isHot('move', k) ? 3.5 : 2.5;
      ctx.beginPath(); ctx.moveTo(pc[0], pc[1]); ctx.lineTo(tip[0], tip[1]); ctx.stroke();
      const ang = Math.atan2(tip[1] - pc[1], tip[0] - pc[0]);
      ctx.beginPath(); ctx.moveTo(tip[0] + Math.cos(ang) * 8, tip[1] + Math.sin(ang) * 8);
      ctx.lineTo(tip[0] + Math.cos(ang + 2.5) * 7, tip[1] + Math.sin(ang + 2.5) * 7);
      ctx.lineTo(tip[0] + Math.cos(ang - 2.5) * 7, tip[1] + Math.sin(ang - 2.5) * 7); ctx.closePath(); ctx.fill();
      text('XYZ'[k] + (App.gizmo.space === 'local' ? '′' : ''), tip[0] + Math.cos(ang) * 16 - 4, tip[1] + Math.sin(ang) * 16 + 4, col);
    }
    ctx.fillStyle = isHot('free', 0) ? '#ffd02e' : 'rgba(255,255,255,0.85)';
    ctx.fillRect(pc[0] - 5, pc[1] - 5, 10, 10);
    if (live && drag.readout) text(drag.readout, pc[0] + 14, pc[1] - 14, '#ffd02e');
  }

  function drawAxesTriad() {
    const c = cam(), ox = 42, oy = ov.height / dpr - 42;
    const ax = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    ctx.lineWidth = 2;
    ax.forEach((a, k) => {
      const x = V.dot(a, c.right()), y = V.dot(a, c.up());
      ctx.strokeStyle = AXIS_COL[k]; ctx.beginPath(); ctx.moveTo(ox, oy); ctx.lineTo(ox + x * 26, oy - y * 26); ctx.stroke();
      text('XYZ'[k], ox + x * 32 - 3, oy - y * 32 + 4, AXIS_COL[k]);
    });
  }

  App.drawOverlay = function () {
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, ov.width, ov.height);
    const t = App.tool;
    const planeTool = t === 'helix' || t === 'sheet' || (t === 'draw' && App.opts.drawSurface !== 'screen') || App.showPlane;
    if (planeTool) drawPlane();

    // Chain breaks.
    if (App.view.breaks) {
      ctx.setLineDash([4, 4]); ctx.strokeStyle = '#ef5b5b'; ctx.lineWidth = 1.5;
      for (const b of App.derived.breaks) {
        if (App.hidden.has(b.chain.id)) continue;
        ctx.beginPath(); const s = line3(b.a, b.b); ctx.stroke();
        if (s && App.derived.breaks.length < 40) text(b.d.toFixed(1) + ' Å', (s[0] + s[2]) / 2 + 4, (s[1] + s[3]) / 2 - 4, '#ff8a8a');
      }
      ctx.setLineDash([]);
    }
    if (App.view.clashes) {
      ctx.strokeStyle = '#f0a04b'; ctx.lineWidth = 2;
      for (const c of App.derived.clashes.slice(0, 400)) {
        const s = cam().project(V.lerp(c.a, c.b, 0.5)); if (!s) continue;
        ctx.beginPath(); ctx.moveTo(s[0] - 4, s[1] - 4); ctx.lineTo(s[0] + 4, s[1] + 4); ctx.moveTo(s[0] + 4, s[1] - 4); ctx.lineTo(s[0] - 4, s[1] + 4); ctx.stroke();
      }
    }
    // Labels on a small selection.
    if (App.view.labels && App.sel.size && App.sel.size <= 16) {
      for (const r of App.selectedResidues()) {
        const info = App.index.get(r.uid);
        const p = (!r.het && Geom.pos(r, 'CA')) || r.atoms[0].p;
        const s = cam().project(p); if (!s) continue;
        text((info.chain ? info.chain.id : (r.chainId || '')) + r.num + (r.icode || '') + ' ' + r.name, s[0] + 8, s[1] - 8, '#ffe27a');
      }
    }
    // Draw tool: chain termini handles.
    if (t === 'draw') {
      for (const chain of App.model.chains) {
        if (App.hidden.has(chain.id) || !chain.residues.length) continue;
        const ends = [[chain.residues[0], 'N'], [chain.residues[chain.residues.length - 1], 'C']];
        for (const [res, w] of ends) {
          const ca = Geom.pos(res, 'CA'); if (!ca) continue;
          const s = cam().project(ca); if (!s) continue;
          ctx.fillStyle = w === 'C' ? '#ffb347' : '#7fd3ff'; ctx.strokeStyle = '#000';
          ctx.beginPath();
          if (w === 'C') { ctx.moveTo(s[0], s[1] - 7); ctx.lineTo(s[0] + 6, s[1] + 5); ctx.lineTo(s[0] - 6, s[1] + 5); ctx.closePath(); }
          else ctx.arc(s[0], s[1], 5.5, 0, 7);
          ctx.fill(); ctx.stroke();
          text(chain.id + (w === 'C' ? ' C-term' : ' N-term'), s[0] + 9, s[1] + 4, w === 'C' ? '#ffcf8a' : '#b4e6ff');
        }
      }
    }
    // Live drag feedback.
    if (drag && drag.kind === 'box') {
      ctx.strokeStyle = drag.remove ? '#ef5b5b' : '#ffd02e'; ctx.setLineDash([5, 3]); ctx.lineWidth = 1;
      ctx.strokeRect(Math.min(drag.x0, drag.x), Math.min(drag.y0, drag.y), Math.abs(drag.x - drag.x0), Math.abs(drag.y - drag.y0));
      ctx.setLineDash([]);
    }
    if (drag && drag.kind === 'draw') {
      ctx.strokeStyle = '#7fe0a0'; ctx.lineWidth = 2.5; ctx.lineJoin = 'round';
      ctx.beginPath(); drag.pts.forEach(([x, y], i) => i ? ctx.lineTo(x, y) : ctx.moveTo(x, y)); ctx.stroke();
    }
    if (drag && drag.kind === 'helix') {
      drawFootprint(drag.base);
      ctx.strokeStyle = 'rgba(83,155,245,0.8)'; ctx.setLineDash([3, 3]); ctx.beginPath();
      line3(V.addScaled(drag.base, App.plane.n, -40), V.addScaled(drag.base, App.plane.n, 40)); ctx.stroke(); ctx.setLineDash([]);
    }
    if (drag && drag.kind === 'sheet') {
      ctx.strokeStyle = '#ffd02e'; ctx.lineWidth = 2; ctx.beginPath(); line3(drag.p0, drag.p1); ctx.stroke();
    }
    if (drag && drag.label && App.mouse) text(drag.label, App.mouse[0] + 16, App.mouse[1] + 22, '#cfe6ff');
    if (!drag && App.mouse && (t === 'helix' || t === 'sheet')) {
      const p = App.planePoint(App.mouse[0], App.mouse[1]);
      if (p) {
        if (t === 'helix') drawFootprint(p);
        else { const s = cam().project(p); if (s) { ctx.fillStyle = '#ffd02e'; ctx.beginPath(); ctx.arc(s[0], s[1], 3.5, 0, 7); ctx.fill(); } }
        const o = App.planeOrigin(), dd = V.sub(p, o);
        text('(' + V.dot(dd, App.plane.u).toFixed(1) + ', ' + V.dot(dd, App.plane.v).toFixed(1) + ')', App.mouse[0] + 14, App.mouse[1] - 10, '#9fb4cc');
      }
    }
    drawGizmo();
    if (App.view.axes) drawAxesTriad();
  };

  App.refreshHud = function () {
    const h = $('hud');
    const uid = App.hover;
    if (uid === null || uid === undefined) { h.textContent = ''; return; }
    const info = App.index.get(uid); if (!info) { h.textContent = ''; return; }
    const r = info.res;
    let s = (info.chain ? 'Chain ' + info.chain.id : 'Ligand') + '  ' + r.name + ' ' + r.num + (r.icode || '');
    if (info.chain) {
      s += '  [' + ({ H: 'helix', E: 'strand', L: 'loop' }[r.ss] || r.ss) + ']';
      const f = Model.phi(info.chain, info.idx), p = Model.psi(info.chain, info.idx);
      s += '\nφ ' + (f === null ? '—' : f.toFixed(0)) + '°   ψ ' + (p === null ? '—' : p.toFixed(0)) + '°';
    }
    h.textContent = s;
  };

  App.setTool = function (t) {
    App.tool = t;
    drag = null; App.ghost = null;
    document.querySelectorAll('#toolButtons button').forEach(b => b.classList.toggle('active', b.dataset.tool === t));
    const hints = {
      select: 'Click to select · Shift-drag: box · gizmo: move / rotate (Ctrl snaps) · right-drag: orbit · wheel: zoom',
      helix: 'Click the plane for a ' + App.opts.helixLen + '-residue helix, or press and drag along the normal to set its length · right-drag: orbit',
      sheet: 'Drag across the plane to draw the first strand; ' + (App.opts.sheetStrands - 1) + ' more are laid beside it · click for a default sheet',
      draw: 'Draw a path. Start on a C-terminus (▲) and/or end on an N-terminus (●) to attach the loop · right-drag: orbit',
    };
    $('hint').textContent = hints[t] || '';
    App.sceneDirty = true; App.overlayDirty = true;
    App.refreshPanels && App.refreshPanels();
  };
})();
