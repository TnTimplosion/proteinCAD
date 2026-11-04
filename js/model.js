/* The editable model and the operations on it.

   model    = { chains: [chain], ligands: [residue], nextUid }
   chain    = { id, residues: [residue] }               (residues in N->C order)
   residue  = { uid, name, num, icode, ss: 'H'|'E'|'L', atoms: [{ name, el, p:[x,y,z] }], het? }

   Ligands are hetero groups carried along from an import: shown, selectable and movable, written
   back out as HETATM, but not part of any chain. Every operation here mutates the model in place;
   the app takes an undo snapshot first. */
(function () {
  'use strict';

  const CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
  const BOND_MAX = 2.0;   // C(i)-N(i+1) longer than this is a chain break

  function create() { return { chains: [], ligands: [], nextUid: 1 }; }

  function newChainId(model, avoid) {
    const used = new Set(model.chains.map(c => c.id));
    if (avoid) for (const a of avoid) used.add(a);
    for (const c of CHAIN_IDS) if (!used.has(c)) return c;
    let k = 1; while (used.has('X' + k)) k++;
    return 'X' + k;
  }

  function stamp(model, res) { res.uid = model.nextUid++; return res; }

  function addChain(model, residues, id) {
    const chain = { id: id || newChainId(model), residues: [] };
    for (const r of residues) chain.residues.push(r.uid ? r : stamp(model, r));
    model.chains.push(chain);
    return chain;
  }

  /* uid -> { res, chain, idx } (chain null for ligands). */
  function buildIndex(model) {
    const m = new Map();
    for (const chain of model.chains) chain.residues.forEach((res, idx) => m.set(res.uid, { res, chain, idx }));
    for (const res of model.ligands) m.set(res.uid, { res, chain: null, idx: -1 });
    return m;
  }

  function allResidues(model) {
    const out = [];
    for (const c of model.chains) for (const r of c.residues) out.push(r);
    for (const r of model.ligands) out.push(r);
    return out;
  }

  const pos = Geom.pos;

  /* ---- connectivity ---- */
  function bonded(a, b) {
    const C = a && pos(a, 'C'), N = b && pos(b, 'N');
    return !!(C && N && V.dist(C, N) < BOND_MAX);
  }

  function phi(chain, i) {
    const r = chain.residues[i], p = chain.residues[i - 1];
    if (!p || !bonded(p, r)) return null;
    const a = pos(p, 'C'), b = pos(r, 'N'), c = pos(r, 'CA'), d = pos(r, 'C');
    return a && b && c && d ? V.dihedral(a, b, c, d) : null;
  }
  function psi(chain, i) {
    const r = chain.residues[i], n = chain.residues[i + 1];
    if (!n || !bonded(r, n)) return null;
    const a = pos(r, 'N'), b = pos(r, 'CA'), c = pos(r, 'C'), d = pos(n, 'N');
    return a && b && c && d ? V.dihedral(a, b, c, d) : null;
  }
  function omega(chain, i) {
    const r = chain.residues[i], n = chain.residues[i + 1];
    if (!n || !bonded(r, n)) return null;
    const a = pos(r, 'CA'), b = pos(r, 'C'), c = pos(n, 'N'), d = pos(n, 'CA');
    return a && b && c && d ? V.dihedral(a, b, c, d) : null;
  }

  function chainBreaks(model) {
    const out = [];
    for (const chain of model.chains) {
      const rs = chain.residues;
      for (let i = 0; i + 1 < rs.length; i++) {
        const C = pos(rs[i], 'C'), N = pos(rs[i + 1], 'N');
        if (!C || !N) continue;
        const d = V.dist(C, N);
        if (d > BOND_MAX) out.push({ chain, i, a: C, b: N, d });
      }
    }
    return out;
  }

  /* Backbone-level clashes (N CA C O CB) closer than `cut` between residues that are not chain
     neighbours. Uses a spatial hash. */
  function clashes(model, cut) {
    cut = cut || 2.5;
    const cut2 = cut * cut, cell = 3.0;
    const pts = [];
    model.chains.forEach((chain, ci) => chain.residues.forEach((res, ri) => {
      for (const a of res.atoms) if (a.name === 'N' || a.name === 'CA' || a.name === 'C' || a.name === 'O' || a.name === 'CB')
        pts.push({ p: a.p, ci, ri });
    }));
    model.ligands.forEach((res, li) => { for (const a of res.atoms) pts.push({ p: a.p, ci: -1 - li, ri: 0 }); });
    const grid = new Map();
    const key = (x, y, z) => x + ',' + y + ',' + z;
    pts.forEach((q, k) => {
      const g = key(Math.floor(q.p[0] / cell), Math.floor(q.p[1] / cell), Math.floor(q.p[2] / cell));
      if (!grid.has(g)) grid.set(g, []);
      grid.get(g).push(k);
    });
    const out = [];
    pts.forEach((q, k) => {
      const gx = Math.floor(q.p[0] / cell), gy = Math.floor(q.p[1] / cell), gz = Math.floor(q.p[2] / cell);
      for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
        const list = grid.get(key(gx + dx, gy + dy, gz + dz));
        if (!list) continue;
        for (const j of list) {
          if (j <= k) continue;
          const o = pts[j];
          if (o.ci === q.ci && Math.abs(o.ri - q.ri) < 2) continue;
          if (V.dist2(q.p, o.p) < cut2) out.push({ a: q.p, b: o.p, d: V.dist(q.p, o.p) });
        }
      }
    });
    return out;
  }

  /* ---- secondary structure from backbone torsions ---- */
  function autoSS(chain) {
    const n = chain.residues.length;
    const kind = new Array(n).fill('L');
    for (let i = 0; i < n; i++) {
      const f = phi(chain, i), s = psi(chain, i);
      if (f === null || s === null) continue;
      if (f > -100 && f < -30 && s > -80 && s < -5) kind[i] = 'h';
      else if (f < -90 && (s > 90 || s < -160)) kind[i] = 'e';
    }
    const ss = new Array(n).fill('L');
    const runs = (ch, min, out) => {
      for (let i = 0; i < n;) {
        if (kind[i] !== ch) { i++; continue; }
        let j = i; while (j < n && kind[j] === ch) j++;
        if (j - i >= min) for (let k = i; k < j; k++) ss[k] = out;
        i = j;
      }
    };
    runs('h', 4, 'H'); runs('e', 3, 'E');
    chain.residues.forEach((r, i) => { r.ss = ss[i]; });
  }

  /* ---- coordinate edits ---- */
  function moveAtoms(atoms, fn) { for (const a of atoms) a.p = fn(a.p); }

  /* Rebuild the carbonyl O of residue i from its neighbour (keeps O consistent after edits). */
  function rebuildO(chain, i) {
    const r = chain.residues[i]; if (!r) return;
    const N = pos(r, 'N'), CA = pos(r, 'CA'), C = pos(r, 'C');
    if (!N || !CA || !C) return;
    const nx = chain.residues[i + 1];
    const Nn = nx && bonded(r, nx) ? pos(nx, 'N') : null;
    const O = Geom.getAtom(r, 'O');
    if (!Nn && O) return; // terminal: keep whatever O it has
    if (Nn) r.atoms = r.atoms.filter(a => a.name !== 'OXT');
    const p = Geom.carbonylO(N, CA, C, Nn, psi(chain, i) ?? 135);
    if (O) O.p = p; else r.atoms.push({ name: 'O', el: 'O', p });
  }

  /* Set a backbone torsion of residue i ('phi', 'psi' or 'omega' to the next residue) by rotating
     the C-terminal ('C') or N-terminal ('N') part of the chain. */
  function setTorsion(chain, i, kind, value, side) {
    const rs = chain.residues, r = rs[i];
    let cur, a, b, moving = [];
    const others = (res, keep) => res.atoms.filter(x => !keep.includes(x.name));
    if (kind === 'phi') {
      cur = phi(chain, i); a = pos(r, 'N'); b = pos(r, 'CA');
      if (side === 'C') { moving = others(r, ['N', 'H', 'CA']); for (let k = i + 1; k < rs.length; k++) moving.push(...rs[k].atoms); }
      else { moving = r.atoms.filter(x => x.name === 'H'); for (let k = 0; k < i; k++) moving.push(...rs[k].atoms); }
    } else if (kind === 'psi') {
      cur = psi(chain, i); a = pos(r, 'CA'); b = pos(r, 'C');
      if (side === 'C') { moving = r.atoms.filter(x => x.name === 'O' || x.name === 'OXT'); for (let k = i + 1; k < rs.length; k++) moving.push(...rs[k].atoms); }
      else { moving = others(r, ['CA', 'C', 'O', 'OXT']); for (let k = 0; k < i; k++) moving.push(...rs[k].atoms); }
    } else {
      cur = omega(chain, i); a = pos(r, 'C'); b = rs[i + 1] && pos(rs[i + 1], 'N');
      if (side === 'C') { moving = others(rs[i + 1], ['N']); for (let k = i + 2; k < rs.length; k++) moving.push(...rs[k].atoms); }
      else { moving = others(r, ['C', 'O', 'OXT']); for (let k = 0; k < i; k++) moving.push(...rs[k].atoms); }
    }
    if (cur === null || !a || !b) return false;
    let delta = V.wrap(value - cur);
    if (side !== 'C') delta = -delta;
    const R = V.rotAxis(V.sub(b, a), delta);
    moveAtoms(moving, p => V.add(b, V.m3vec(R, V.sub(p, b))));
    return true;
  }

  /* Add n residues at the C ('C') or N ('N') terminus with conformation conf {phi, psi, ss}. */
  function extendChain(model, chain, end, n, conf, name) {
    const rs = chain.residues;
    if (!rs.length) return 0;
    if (end === 'C') {
      const last = rs[rs.length - 1];
      let N = pos(last, 'N'), CA = pos(last, 'CA'), C = pos(last, 'C');
      if (!N || !CA || !C) return 0;
      const added = [];
      for (let k = 0; k < n; k++) {
        const nx = Geom.extendForward(N, CA, C, conf.psi, conf.phi);
        added.push(nx); ({ N, CA, C } = nx);
      }
      const res = Geom.residuesFromBackbone(added, name, conf.ss, conf.psi);
      res.forEach((r, k) => { r.num = last.num + 1 + k; stamp(model, r); rs.push(r); });
      rebuildO(chain, rs.length - n - 1);
      return n;
    }
    const first = rs[0];
    let N = pos(first, 'N'), CA = pos(first, 'CA'), C = pos(first, 'C');
    if (!N || !CA || !C) return 0;
    const added = [];
    for (let k = 0; k < n; k++) {
      const pv = Geom.extendBackward(N, CA, C, conf.phi, conf.psi);
      added.unshift(pv); ({ N, CA, C } = pv);
    }
    const res = Geom.residuesFromBackbone(added.concat([{ N: pos(first, 'N'), CA: pos(first, 'CA'), C: pos(first, 'C') }]), name, conf.ss).slice(0, n);
    res.forEach((r, k) => { r.num = first.num - n + k; stamp(model, r); });
    rs.unshift(...res);
    return n;
  }

  /* ---- residue-type and chain edits ---- */
  function mutate(res, name) {
    if (res.name === name && Geom.SC[name]) {
      return Geom.buildSideChain(res, name, Geom.chiAngles(res));
    }
    return Geom.buildSideChain(res, name);
  }

  /* 'backbone': keep N CA C O (rename to GLY if renameGly). 'ala': N CA C O CB as ALA. 'gly'. */
  function convertResidues(residues, mode, renameGly) {
    let n = 0;
    for (const r of residues) {
      if (r.het) continue;
      if (mode === 'backbone') {
        r.atoms = r.atoms.filter(a => a.name === 'N' || a.name === 'CA' || a.name === 'C' || a.name === 'O' || a.name === 'OXT');
        if (renameGly) r.name = 'GLY';
        n++;
      } else if (mode === 'ala' || mode === 'gly') {
        if (Geom.buildSideChain(r, mode === 'ala' ? 'ALA' : 'GLY')) n++;
      }
    }
    return n;
  }

  function renumber(residues, start) {
    residues.forEach((r, k) => { r.num = start + k; r.icode = ''; });
  }

  /* Remove residues by uid; drops chains that become empty. */
  function deleteResidues(model, uids) {
    for (const c of model.chains) c.residues = c.residues.filter(r => !uids.has(r.uid));
    model.chains = model.chains.filter(c => c.residues.length);
    model.ligands = model.ligands.filter(r => !uids.has(r.uid));
  }

  /* Split chain after residue index i: residues i+1.. become a new chain. */
  function splitChain(model, chain, i) {
    if (i < 0 || i >= chain.residues.length - 1) return null;
    const tail = chain.residues.splice(i + 1);
    const nc = { id: newChainId(model), residues: tail };
    model.chains.splice(model.chains.indexOf(chain) + 1, 0, nc);
    return nc;
  }

  /* Append chain b to chain a (b disappears). */
  function mergeChains(model, a, b) {
    a.residues.push(...b.residues);
    model.chains = model.chains.filter(c => c !== b);
  }

  /* Copy residues into a new chain (new uids). */
  function duplicate(model, residues) {
    const copy = residues.map(r => stamp(model, { name: r.name, num: r.num, icode: r.icode, ss: r.ss, atoms: r.atoms.map(a => ({ name: a.name, el: a.el, p: a.p.slice() })) }));
    return addChain(model, copy);
  }

  function bounds(model) {
    const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    let n = 0;
    for (const r of allResidues(model)) for (const a of r.atoms) {
      n++;
      for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], a.p[k]); hi[k] = Math.max(hi[k], a.p[k]); }
    }
    if (!n) return null;
    return { lo, hi, center: V.scale(V.add(lo, hi), 0.5), radius: Math.max(5, V.dist(lo, hi) / 2) };
  }

  function sequence(chain) { return chain.residues.map(r => Geom.oneLetter(r.name)).join(''); }

  window.Model = {
    CHAIN_IDS, BOND_MAX, create, newChainId, stamp, addChain, buildIndex, allResidues,
    bonded, phi, psi, omega, chainBreaks, clashes, autoSS, rebuildO, setTorsion, extendChain,
    mutate, convertResidues, renumber, deleteResidues, splitChain, mergeChains, duplicate, bounds,
    sequence, moveAtoms,
  };
})();
