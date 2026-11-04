/* Side-chain repacking: a discrete search over chi angles that clears clashes a rigid side chain
   cannot.

   Rotamers come from a staggered set — the three positions about each sp3 bond, plus the ring and
   amide orientations that are actually populated for each residue type — expanded combinatorially
   and pruned of combinations that clash with their own backbone. This is not a backbone-dependent
   library with real frequencies; chi1 preferences enter only as a small tie-break, and the score is
   dominated by sterics. For a design template, where the sequence is about to be replaced by an
   inverse-folding model anyway, that is the right level of detail.

   Only the selected residues move. Everything else, including the side chains of residues outside
   the selection, is fixed environment — so a clash with a fixed neighbour is resolved by moving the
   selected side chain, which is the only thing the user asked to change. CB does not move either:
   its position follows from N, CA and C. */
(function () {
  'use strict';

  /* Candidate values for each chi, by residue type. Aromatic rings and carboxylates are two-fold
     symmetric, so their last chi only needs to span 180 degrees. */
  const CHI = {
    SER: [[62, -177, -65]],
    CYS: [[62, -177, -65]],
    THR: [[62, -175, -61]],
    VAL: [[63, 175, -60]],
    ILE: [[62, -177, -65], [100, 170, -60]],
    LEU: [[62, -177, -65], [65, 175, -60]],
    MET: [[62, -177, -65], [65, 180, -65], [75, 180, -75]],
    PHE: [[62, -177, -65], [90, -85, -30]],
    TYR: [[62, -177, -65], [90, -85, -30]],
    TRP: [[62, -177, -65], [-105, -90, -5, 90]],
    HIS: [[62, -177, -65], [-75, 80, 170]],
    ASN: [[62, -177, -65], [-80, -20, 30, 120]],
    ASP: [[62, -177, -65], [-70, -15, 30]],
    GLN: [[62, -177, -65], [65, 180, -65], [-60, 0, 60]],
    GLU: [[62, -177, -65], [65, 180, -65], [-60, 0, 60]],
    LYS: [[62, -177, -65], [68, 180, -68], [68, 180, -68], [65, 180, -65]],
    ARG: [[62, -177, -65], [65, 180, -65], [65, 180, -65], [85, 180, -85]],
  };

  /* Small preference for the chi1 wells, as a tie-break only: m is the most common, p the least,
     except for the beta-branched residues where t dominates. */
  const CHI1_COST = { p: 0.8, t: 0.3, m: 0.0 };
  const CHI1_COST_BRANCHED = { p: 0.5, t: 0.0, m: 0.5 };
  const BRANCHED = new Set(['VAL', 'ILE', 'THR']);
  const MIN_GAIN = 0.6;       // a rotamer must beat the current one by this much to replace it
  const SELF_CUTOFF = 0.8;    // a rotamer that clashes with its own backbone this badly is dropped
  const MAX_ROTAMERS = 24;    // kept per residue after scoring against the fixed environment
  const PAIR_CUT = 13.0;      // CB-CB distance within which two side chains can reach each other

  const BACKBONE = ['N', 'CA', 'C', 'O', 'OXT', 'CB'];

  function chi1Well(v) {
    const w = V.wrap(v);
    if (w > 0 && w < 120) return 'p';
    if (w >= 120 || w < -120) return 't';
    return 'm';
  }

  /* Which atom pairs inside a residue are too close in the bond graph to be scored against each
     other. Built once per residue type from a reference geometry; bond lengths do not depend on
     chi, so the graph does not either. */
  const topoCache = new Map();
  function topology(name) {
    if (topoCache.has(name)) return topoCache.get(name);
    const bb = Geom.backboneFromTorsions([{ phi: -60, psi: -45, omega: 180 }])[0];
    const res = { name, atoms: [['N', 'N'], ['CA', 'C'], ['C', 'C']].map(([n, el], i) => ({ name: n, el, p: [bb.N, bb.CA, bb.C][i] })) };
    Geom.buildSideChain(res, name, (Geom.SC[name] || { chi: [] }).chi);
    res.atoms.push({ name: 'O', el: 'O', p: Geom.carbonylO(bb.N, bb.CA, bb.C, null, -45) });
    const names = res.atoms.map(a => a.name);
    const adj = new Map(names.map(n => [n, []]));
    for (let i = 0; i < res.atoms.length; i++) {
      for (let j = i + 1; j < res.atoms.length; j++) {
        const cut = (res.atoms[i].el === 'S' || res.atoms[j].el === 'S') ? 2.1 : 1.9;
        if (V.dist(res.atoms[i].p, res.atoms[j].p) < cut) {
          adj.get(names[i]).push(names[j]); adj.get(names[j]).push(names[i]);
        }
      }
    }
    const close = new Set();
    for (const start of names) {
      const seen = new Map([[start, 0]]);
      const queue = [start];
      while (queue.length) {
        const cur = queue.shift(), d = seen.get(cur);
        if (d >= 3) continue;
        for (const nx of adj.get(cur)) if (!seen.has(nx)) { seen.set(nx, d + 1); queue.push(nx); }
      }
      for (const [n, d] of seen) if (d <= 3) close.add(start + '|' + n);
    }
    topoCache.set(name, close);
    return close;
  }

  function packable(res) {
    return !res.het && CHI[res.name] && Geom.pos(res, 'N') && Geom.pos(res, 'CA') && Geom.pos(res, 'C');
  }

  /* A side chain whose chi angles are not free: one end of a disulfide, or covalently attached to
     (or coordinating) a hetero group. Turning these would pull the bridge apart — the packer has no
     notion of a bond, so an S-S at 2.03 A reads to it as a terrible clash and it would happily
     rotate out of it. As in refinement, the test has to be chemical rather than by distance, since
     a bond and an overlap look alike. */
  function committed(res, others) {
    for (const a of res.atoms) {
      if (BACKBONE.includes(a.name)) continue;
      for (const o of others) {
        if (o.res === res) continue;
        const ss = a.el === 'S' && o.el === 'S';
        if (!ss && !o.res.het) continue;
        if (V.dist2(a.p, o.p) < (ss ? 2.6 : 2.4) ** 2) return true;
      }
    }
    return false;
  }

  /* Side-chain atoms (beyond CB) for this residue at these chi angles. */
  function placeRotamer(res, chi) {
    const tmp = { name: res.name, atoms: ['N', 'CA', 'C'].map(n => ({ name: n, el: n[0], p: Geom.pos(res, n) })) };
    if (!Geom.buildSideChain(tmp, res.name, chi)) return null;
    return tmp.atoms.filter(a => !BACKBONE.includes(a.name)).map(a => ({ name: a.name, el: a.el, p: a.p }));
  }

  function cartesian(lists) {
    let out = [[]];
    for (const list of lists) {
      const next = [];
      for (const prefix of out) for (const v of list) next.push(prefix.concat([v]));
      out = next;
    }
    return out;
  }

  function grid(atoms, cell) {
    const g = new Map();
    const key = (x, y, z) => x + ',' + y + ',' + z;
    for (const a of atoms) {
      const k = key(Math.floor(a.p[0] / cell), Math.floor(a.p[1] / cell), Math.floor(a.p[2] / cell));
      if (!g.has(k)) g.set(k, []);
      g.get(k).push(a);
    }
    return {
      near(p) {
        const out = [];
        const gx = Math.floor(p[0] / cell), gy = Math.floor(p[1] / cell), gz = Math.floor(p[2] / cell);
        for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
          const l = g.get(key(gx + dx, gy + dy, gz + dz));
          if (l) for (const a of l) out.push(a);
        }
        return out;
      },
    };
  }

  const overlap = (a, b, scale) => {
    const d0 = Refine.clashTarget(a, b, scale);
    const d2 = V.dist2(a.p, b.p);
    if (d2 >= d0 * d0) return 0;
    const over = d0 - Math.sqrt(d2);
    return over * over;
  };

  /* Clash energy of a set of atoms against a list of environment atoms. */
  function energyAgainst(atoms, env, scale) {
    let e = 0;
    for (const a of atoms) for (const b of env) e += overlap(a, b, scale);
    return e;
  }

  /* Repack the side chains of `residues`. Returns a report, or null when none can be packed. */
  function run(model, residues, opts) {
    opts = opts || {};
    const scale = opts.clashScale !== undefined ? opts.clashScale : Refine.DEFAULTS.clashScale;
    const sweeps = opts.sweeps || 6;
    const all = [];
    for (const res of Model.allResidues(model)) for (const a of res.atoms) all.push({ p: a.p, el: a.el, res, name: a.name });
    const candidates = residues.filter(packable);
    const held = [];
    const targets = candidates.filter(r => {
      if (!committed(r, all)) return true;
      held.push(r); return false;
    });
    if (!targets.length) return { sites: 0, changed: 0, energyBefore: 0, energyAfter: 0, rotamers: 0, committed: held.length };
    const moving = new Set(targets.map(r => r.uid));

    // Fixed environment: every atom in the model that this run will not move.
    const envAtoms = all.filter(a => !moving.has(a.res.uid) || BACKBONE.includes(a.name));
    const env = grid(envAtoms, 6.0);

    // Rotamers per residue, scored against the fixed environment and their own backbone.
    const sites = [];
    for (const res of targets) {
      const close = topology(res.name);
      const own = res.atoms.filter(a => BACKBONE.includes(a.name)).map(a => ({ p: a.p, el: a.el, name: a.name }));
      const current = Geom.chiAngles(res) || [];
      // The rotamer it already has goes in first, so "no change" is always available.
      const combos = cartesian(CHI[res.name]);
      const branched = BRANCHED.has(res.name);
      const rots = [];
      [current.slice()].concat(combos).forEach((chi, idx) => {
        const atoms = placeRotamer(res, chi);
        if (!atoms) return;
        let self = 0;
        for (let i = 0; i < atoms.length; i++) {
          for (const b of own) if (!close.has(atoms[i].name + '|' + b.name)) self += overlap(atoms[i], b, scale);
          for (let j = i + 1; j < atoms.length; j++) if (!close.has(atoms[i].name + '|' + atoms[j].name)) self += overlap(atoms[i], atoms[j], scale);
        }
        if (self > SELF_CUTOFF && idx > 0) return;   // impossible combination of chi angles
        const cost = (branched ? CHI1_COST_BRANCHED : CHI1_COST)[chi1Well(chi[0])] || 0;
        rots.push({ chi, atoms, self, prior: cost, isCurrent: idx === 0 });
      });
      if (rots.length < 2) continue;
      // Score against the fixed environment once; this never changes during the search.
      for (const r of rots) {
        let e = 0;
        for (const a of r.atoms) for (const b of env.near(a.p)) if (b.res !== res) e += overlap(a, b, scale);
        r.fixed = e + r.self + r.prior;
      }
      rots.sort((a, b) => a.fixed - b.fixed);
      const kept = rots.slice(0, MAX_ROTAMERS);
      if (!kept.some(r => r.isCurrent)) kept.push(rots.find(r => r.isCurrent));
      sites.push({ res, rots: kept, cur: kept.findIndex(r => r.isCurrent), cb: Geom.pos(res, 'CB') || Geom.pos(res, 'CA') });
    }
    if (!sites.length) return null;

    // Which sites can reach each other.
    const nbr = sites.map(() => []);
    for (let i = 0; i < sites.length; i++) {
      for (let j = i + 1; j < sites.length; j++) {
        if (V.dist(sites[i].cb, sites[j].cb) < PAIR_CUT) { nbr[i].push(j); nbr[j].push(i); }
      }
    }

    const pairEnergy = (i, ri, j, rj) => energyAgainst(sites[i].rots[ri].atoms, sites[j].rots[rj].atoms, scale);
    const total = pick => {
      let e = 0;
      for (let i = 0; i < sites.length; i++) {
        e += sites[i].rots[pick[i]].fixed;
        for (const j of nbr[i]) if (j > i) e += pairEnergy(i, pick[i], j, pick[j]);
      }
      return e;
    };

    /* Go round the sites choosing the best rotamer for each given its neighbours, until nothing
       changes. Cheap, and with a few starting points it does not get stuck anywhere interesting. */
    const optimise = (pick, rng) => {
      for (let sweep = 0; sweep < sweeps; sweep++) {
        let changed = 0;
        const order = sites.map((_, i) => i);
        for (let i = order.length - 1; i > 0; i--) { const j = Math.floor(rng() * (i + 1)); [order[i], order[j]] = [order[j], order[i]]; }
        for (const i of order) {
          let best = pick[i], bestE = Infinity;
          for (let r = 0; r < sites[i].rots.length; r++) {
            let e = sites[i].rots[r].fixed;
            for (const j of nbr[i]) e += pairEnergy(i, r, j, pick[j]);
            if (e < bestE) { bestE = e; best = r; }
          }
          if (best !== pick[i]) { pick[i] = best; changed++; }
        }
        if (!changed) break;
      }
      return pick;
    };

    let seed = 20261030;
    const rng = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; };
    const starts = [
      sites.map(s => s.cur),                     // what it has now
      sites.map(() => 0),                        // best against the fixed environment
      sites.map(s => Math.floor(rng() * s.rots.length)),
      sites.map(s => Math.floor(rng() * s.rots.length)),
    ];
    const before = total(sites.map(s => s.cur));
    let best = null, bestE = Infinity;
    for (const start of starts) {
      const pick = optimise(start.slice(), rng);
      const e = total(pick);
      if (e < bestE) { bestE = e; best = pick; }
    }

    /* Put back every rotamer the search moved for no real gain. The user placed these side chains;
       a search over a coarse staggered library will happily swap a perfectly good one for a
       marginally better score, which shows up as the model changing more than it needed to. */
    for (let round = 0; round < 2; round++) {
      let reverted = 0;
      for (let i = 0; i < sites.length; i++) {
        if (best[i] === sites[i].cur) continue;
        const at = r => {
          let e = sites[i].rots[r].fixed;
          for (const j of nbr[i]) e += pairEnergy(i, r, j, best[j]);
          return e;
        };
        if (at(best[i]) > at(sites[i].cur) - MIN_GAIN) { best[i] = sites[i].cur; reverted++; }
      }
      if (!reverted) break;
    }
    bestE = total(best);

    // Write the chosen rotamers back.
    let changed = 0;
    sites.forEach((s, i) => {
      if (best[i] === s.cur) return;
      if (Geom.buildSideChain(s.res, s.res.name, s.rots[best[i]].chi)) changed++;
    });
    return { sites: sites.length, changed, energyBefore: before, energyAfter: bestE,
      rotamers: sites.reduce((n, s) => n + s.rots.length, 0), committed: held.length };
  }

  window.Pack = { run, CHI, topology, packable, placeRotamer };
})();
