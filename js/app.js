/* Application state, undo, selection, rendering schedule, files. Tools (tools.js) and panels
   (panels.js) hang their functions on the same App object; main.js wires it all up. */
(function () {
  'use strict';
  const $ = id => document.getElementById(id);

  const App = {
    $,
    model: Model.create(),
    index: new Map(),
    sel: new Set(),
    selOrder: [],           // uids in the order they were selected
    hover: null,
    tool: 'select',
    selMode: 'residue',     // residue | element | chain
    plane: { base: [0, 0, 0], n: [0, 0, 1], u: [1, 0, 0], v: [0, 1, 0], offset: 0 },
    view: { cartoon: true, atoms: 'all', style: 'sticks', trace: false, color: 'chain', breaks: true, clashes: true, labels: true, axes: true },
    gizmo: { show: true, space: 'world', snapMove: 1.0, snapRot: 15 },
    hidden: new Set(),
    opts: {
      newRes: 'ALA', snap: 1.0,
      helixConf: 'alpha', helixLen: 14, helixEnd: 'start', helixPhase: 0,
      sheetStrands: 4, sheetParallel: false, sheetSpacing: 4.8, sheetSide: 1, sheetTwist: 0, sheetShift: 0, sheetLen: 7,
      loopLen: 0, loopRes: 'ALA', renumber: true,
      drawSurface: 'auto', drawSpacing: 3.2,
      extendN: 4, extendConf: 'alpha', torsionSide: 'C', applyConf: 'alpha',
      refineHold: 0.15, refineClash: 0.88, refineGap: 4.0, refineIters: 400, refineSide: 'ideal',
    },
    undoStack: [], redoStack: [],
    ghost: null,
    derived: { breaks: [], clashes: [] },
    renderer: null,
    sceneDirty: true, overlayDirty: true,
  };

  App.planeOrigin = function () { return V.addScaled(App.plane.base, App.plane.n, App.plane.offset); };

  /* ---------------- status ---------------- */
  App.status = function (text, kind) {
    const m = $('msg');
    m.textContent = text;
    m.className = kind || '';
  };

  /* ---------------- undo ---------------- */
  App.snapshot = () => structuredClone(App.model);
  App.pushUndo = function (snap, label) {
    App.undoStack.push({ model: snap, label });
    if (App.undoStack.length > 40) App.undoStack.shift();
    App.redoStack.length = 0;
  };
  /* Run fn as one undoable edit. fn may return false to signal "nothing changed". */
  App.edit = function (label, fn, opts) {
    const snap = App.snapshot();
    let r;
    try { r = fn(); } catch (e) { App.model = snap; App.changed(); App.status(label + ' failed: ' + e.message, 'err'); console.error(e); return false; }
    if (r === false) return false;
    App.pushUndo(snap, label);
    App.changed(opts);
    return true;
  };
  App.undo = function () {
    const u = App.undoStack.pop(); if (!u) return;
    App.redoStack.push({ model: App.snapshot(), label: u.label });
    App.model = u.model; App.pruneSelection(); App.changed();
    App.status('Undid: ' + u.label);
  };
  App.redo = function () {
    const r = App.redoStack.pop(); if (!r) return;
    App.undoStack.push({ model: App.snapshot(), label: r.label });
    App.model = r.model; App.pruneSelection(); App.changed();
    App.status('Redid: ' + r.label);
  };

  /* ---------------- change propagation ---------------- */
  App.modelVersion = 0;
  App.changed = function (opts) {
    App.modelVersion++;
    App.index = Model.buildIndex(App.model);
    App.pruneSelection();
    App.derived.breaks = Model.chainBreaks(App.model);
    App.derived.clashes = App.view.clashes ? Model.clashes(App.model) : [];
    App.sceneDirty = true; App.overlayDirty = true; App.projCache = null;
    App.refineCheckText = '';
    if (!(opts && opts.quiet)) App.refreshPanels && App.refreshPanels();
    App.autosave();
  };
  /* Light refresh while dragging: geometry only. */
  App.moved = function () { App.sceneDirty = true; App.overlayDirty = true; App.projCache = null; };

  App.selectionChanged = function () {
    App.sceneDirty = true; App.overlayDirty = true;
    App.refineCheckText = '';   // the geometry report belongs to the selection it was run on
    App.refreshPanels && App.refreshPanels();
  };

  /* ---------------- selection ---------------- */
  App.pruneSelection = function () {
    for (const u of [...App.sel]) if (!App.index.has(u)) App.sel.delete(u);
    App.selOrder = App.selOrder.filter(u => App.sel.has(u));
    if (App.hover !== null && !App.index.has(App.hover)) App.hover = null;
  };
  App.setSelection = function (uids, mode) {
    // mode: 'set' | 'add' | 'remove' | 'toggle'
    if (mode === 'set' || !mode) { App.sel.clear(); App.selOrder = []; }
    for (const u of uids) {
      if (mode === 'remove' || (mode === 'toggle' && App.sel.has(u))) App.sel.delete(u);
      else if (!App.sel.has(u)) { App.sel.add(u); App.selOrder.push(u); }
    }
    App.selOrder = App.selOrder.filter(u => App.sel.has(u));
    App.selectionChanged();
  };
  /* Expand a residue uid to the current selection unit. */
  App.expand = function (uid, mode) {
    const info = App.index.get(uid);
    if (!info || !info.chain) return [uid];
    mode = mode || App.selMode;
    const rs = info.chain.residues;
    if (mode === 'chain') return rs.map(r => r.uid);
    if (mode === 'element') {
      const ss = info.res.ss;
      let a = info.idx, b = info.idx;
      while (a > 0 && rs[a - 1].ss === ss && Model.bonded(rs[a - 1], rs[a])) a--;
      while (b + 1 < rs.length && rs[b + 1].ss === ss && Model.bonded(rs[b], rs[b + 1])) b++;
      return rs.slice(a, b + 1).map(r => r.uid);
    }
    return [uid];
  };
  App.selectedResidues = function () {
    // In chain order, then ligands.
    const out = [];
    for (const c of App.model.chains) for (const r of c.residues) if (App.sel.has(r.uid)) out.push(r);
    for (const r of App.model.ligands) if (App.sel.has(r.uid)) out.push(r);
    return out;
  };
  /* Chains with at least one selected residue, in model order (or by selection order if asked). */
  App.selectedChains = function (bySelectionOrder) {
    if (bySelectionOrder) {
      const seen = []; for (const u of App.selOrder) { const i = App.index.get(u); if (i && i.chain && !seen.includes(i.chain)) seen.push(i.chain); }
      return seen;
    }
    return App.model.chains.filter(c => c.residues.some(r => App.sel.has(r.uid)));
  };
  App.selectionPoints = function () {
    const pts = [];
    for (const r of App.selectedResidues()) {
      const ca = !r.het && Geom.pos(r, 'CA');
      if (ca) pts.push(ca); else for (const a of r.atoms) pts.push(a.p);
    }
    return pts;
  };

  /* ---------------- rendering ---------------- */
  App.frame = function () {
    const r = App.renderer;
    if (App.sceneDirty) {
      r.upload(Render.buildScene(App.model, {
        atoms: App.view.atoms, style: App.view.style, cartoon: App.view.cartoon, trace: App.view.trace, color: App.view.color,
        sel: App.sel, hover: App.hover, ghost: App.ghost, hidden: App.hidden,
      }));
      App.sceneDirty = false; App.glDirty = true;
    }
    if (App.glDirty) { r.draw(); App.glDirty = false; }
    if (App.overlayDirty) { App.drawOverlay && App.drawOverlay(); App.overlayDirty = false; }
    requestAnimationFrame(App.frame);
  };
  App.cameraChanged = function () { App.glDirty = true; App.overlayDirty = true; App.projCache = null; };

  App.updateSceneRadius = function () {
    const b = Model.bounds(App.model);
    const c = App.renderer.cam;
    const pr = V.dist(App.planeOrigin(), c.target) + 40;
    c.sceneRadius = Math.max(30, b ? b.radius + V.dist(b.center, c.target) : 0, pr);
  };

  App.fitView = function (onlySel) {
    const pts = onlySel && App.sel.size ? App.selectionPoints() : [];
    const cam = App.renderer.cam;
    if (pts.length) {
      const c = V.centroid(pts); let r = 6;
      for (const p of pts) r = Math.max(r, V.dist(p, c) + 3);
      cam.fit(c, r);
    } else {
      const b = Model.bounds(App.model);
      if (b) cam.fit(b.center, b.radius + 2); else cam.fit(App.planeOrigin(), 25);
    }
    App.updateSceneRadius();
    App.cameraChanged();
  };

  App.setView = function (name) {
    const cam = App.renderer.cam;
    // Z is up (the default sketch plane is XY, so helices extrude along Z).
    if (name === 'front') cam.lookAlong([0, -1, 0], [0, 0, 1]);
    else if (name === 'top') cam.lookAlong([0, 0, 1], [0, 1, 0]);
    else if (name === 'right') cam.lookAlong([1, 0, 0], [0, 0, 1]);
    else cam.lookAlong([0.62, -0.78, 0.55], [0, 0, 1]);
    App.cameraChanged();
  };

  /* ---------------- files ---------------- */
  App.loadText = function (text, filename, replace) {
    let m;
    try { m = IO.parse(text, filename); } catch (e) { App.status('Could not read ' + filename + ': ' + e.message, 'err'); return; }
    const nres = m.chains.reduce((s, c) => s + c.residues.length, 0);
    if (!nres && !m.ligands.length) { App.status('No atoms found in ' + filename + '.', 'err'); return; }
    App.edit(replace ? 'Open ' + filename : 'Import ' + filename, () => {
      if (replace) { App.model = Model.create(); App.sel.clear(); App.selOrder = []; App.hidden.clear(); }
      const used = new Set(App.model.chains.map(c => c.id));
      const renamed = [];
      for (const c of m.chains) {
        const res = c.residues.map(r => { delete r.uid; return r; });
        let id = c.id;
        if (used.has(id)) { const nid = Model.newChainId(App.model, used); renamed.push(id + '→' + nid); id = nid; }
        used.add(id);
        Model.addChain(App.model, res, id);
      }
      for (const l of m.ligands) { delete l.uid; App.model.ligands.push(Model.stamp(App.model, l)); }
      App._lastRenamed = renamed;
    });
    App.fitView();
    const rn = App._lastRenamed && App._lastRenamed.length ? ' Renamed chains ' + App._lastRenamed.join(', ') + ' to avoid clashes.' : '';
    App.status((replace ? 'Opened ' : 'Imported ') + filename + ': ' + m.chains.length + ' chain(s), ' + nres + ' residues, ' + m.ligands.length + ' ligand group(s).' + rn, 'ok');
  };

  App.download = function (name, text) {
    const blob = new Blob([text], { type: 'text/plain' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob); a.download = name;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 2000);
  };

  App.exportModel = function (fmt) {
    if (!App.model.chains.length && !App.model.ligands.length) { App.status('Nothing to export.', 'err'); return; }
    const opts = { backboneOnly: $('exBackbone').checked, includeLigands: $('exLigands').checked };
    const long = App.model.chains.filter(c => c.id.length > 1).map(c => c.id);
    if (fmt === 'pdb') {
      App.download('design.pdb', IO.toPDB(App.model, opts));
      App.status('Exported design.pdb' + (long.length ? ' (chain IDs ' + long.join(', ') + ' cut to one character; use mmCIF to keep them)' : '') + '.', long.length ? 'err' : 'ok');
    } else {
      App.download('design.cif', IO.toCIF(App.model, opts));
      App.status('Exported design.cif.', 'ok');
    }
  };

  /* ---------------- autosave (per-browser convenience) ---------------- */
  let saveTimer = null;
  App.autosave = function () {
    clearTimeout(saveTimer);
    saveTimer = setTimeout(() => {
      try { localStorage.setItem('protein_cad_autosave', JSON.stringify(App.model)); } catch (e) { /* storage full or blocked */ }
    }, 800);
  };
  App.restore = function () {
    try {
      const s = localStorage.getItem('protein_cad_autosave');
      if (!s) return false;
      const m = JSON.parse(s);
      if (!m || !Array.isArray(m.chains)) return false;
      if (!m.chains.length && !(m.ligands || []).length) return false;
      m.ligands = m.ligands || [];
      App.model = m;
      return true;
    } catch (e) { return false; }
  };

  window.App = App;
})();
