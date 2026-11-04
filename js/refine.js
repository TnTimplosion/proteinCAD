/* Restrained regularisation of a selected region: pulls backbone bond lengths, angles and peptide
   units back to ideal values, keeps the CA centres L-chiral, pushes apart atoms that overlap, and
   closes small chain breaks — while the rest of the model stays fixed and acts as context.

   Everything is expressed as harmonic restraints on interatomic distances, plus one signed-volume
   term per CA for chirality. A bond angle is restrained as the 1-3 distance implied by the ideal
   angle and its two ideal bond lengths, and the peptide torsion as the 1-4 CA(i)-CA(i+1) distance;
   both have trivial, well-conditioned gradients and no singularities, unlike angle and torsion
   derivatives. Side chains are held rigid by restraining their internal distances to the values
   they have now, so they rotate and translate with the backbone but keep their rotamer.

   Minimised with L-BFGS and an Armijo backtracking line search. */
(function () {
  'use strict';

  /* Ideal geometry (Engh & Huber bond lengths and angles, as in Geom.BB), with the 1-3 and 1-4
     distances they imply. Derived numerically from the same constants; see README. */
  const T = {
    N_CA: 1.458, CA_C: 1.525, C_N: 1.329, C_O: 1.231, CA_CB: 1.530,
    N_C: 2.4616,      // N-CA-C   111.2
    N_CB: 2.4524,     // N-CA-CB  110.3
    C_CB: 2.5047,     // CB-CA-C  110.1
    CA_O: 2.3972,     // CA-C-O   120.5
    CA_N1: 2.4252,    // CA-C-N   116.2
    O_N1: 2.2534,     // O-C-N    123.3
    C_CA1: 2.4348,    // C-N-CA   121.7
    CA_CA1_trans: 3.8040, CA_CA1_cis: 2.7714,
    O_CA1_trans: 2.7791, O_CA1_cis: 3.5817,
    O_OXT: 2.1670,    // O-C-OXT 123.3, the carboxy terminus
    CHIRAL_V: 2.5088, // (N-CA) . ((C-CA) x (CB-CA)) for an L-residue
  };

  /* Restraint weights. Relative values only; bonds are the stiffest thing in the system. */
  const W = {
    bond: 100, angle: 45, peptide: 12, chiral: 60, rigid: 50, nb: 12, conf: 4, rama: 30,
  };

  /* phi of a non-glycine residue is confined to |phi| >= 40 and phi <= 100: the band around zero
     is blocked by its own carbonyl, and the far positive side by the side chain. */
  const PHI_MIN = 40, PHI_MAX_POS = 100;

  const VDW = { C: 1.70, N: 1.55, O: 1.52, S: 1.80, P: 1.80, Se: 1.90 };
  const POLAR = new Set(['N', 'O']);
  const BB_ORDER = ['N', 'CA', 'C', 'O', 'CB'];

  const DEFAULTS = {
    iterations: 300,
    hold: 0.15,        // weight of the restraint to each atom's starting position
    // Soft-sphere target = this x the sum of the van der Waals radii. 0.88 puts C...C at 2.99 A and
    // N...O at 2.69: close enough to leave a good structure and its hydrogen bonds alone, far
    // enough to clear anything that reads as a clash.
    clashScale: 0.88,
    maxGap: 4.0,       // close a chain break only when the C...N gap is at most this
    rigidSide: true,   // keep side chains rigid (otherwise only CB is restrained)
    // 'keep'   leave side chains exactly as they are, restrained to their present shape
    // 'ideal'  rebuild each one from its own current chi angles (regularises it, same rotamer)
    // 'repack' search chi for the rotamer that fits, between minimisation passes
    sideMode: 'ideal',
    moveLigands: false,
    contextCut: 8.0,   // fixed atoms this far from the mobile set still push on it
  };

  /* A residue of this type at these chi angles, built on an ideal backbone. Side-chain restraint
     targets are measured here rather than on the residue being refined: distances taken from a
     distorted backbone disagree with the ideal bond lengths and angles, and the two restraint sets
     then fight each other to a standstill with bonds still 0.03 A out. */
  const refCache = new Map();
  function referenceResidue(name, chi) {
    const key = name + '|' + chi.map(c => Math.round(c)).join(',');
    if (refCache.has(key)) return refCache.get(key);
    const bb = Geom.backboneFromTorsions([{ phi: -60, psi: -45, omega: 180 }])[0];
    const res = { name, atoms: ['N', 'CA', 'C'].map((n, i) => ({ name: n, el: n[0], p: [bb.N, bb.CA, bb.C][i] })) };
    const map = Geom.buildSideChain(res, name, chi) ? Object.fromEntries(res.atoms.map(a => [a.name, a.p])) : null;
    refCache.set(key, map);
    return map;
  }

  /* Has this residue every atom its type's side chain should have? */
  function sideChainComplete(res) {
    const lib = Geom.SC[res.name];
    if (!lib) return false;
    if (res.name === 'GLY') return true;
    return ['CB'].concat(lib.atoms.map(a => a[0])).every(n => Geom.getAtom(res, n));
  }

  function vdw(el) { return VDW[el] || 1.70; }
  function clashTarget(a, b, scale) {
    const d = (vdw(a.el) + vdw(b.el)) * scale;
    return POLAR.has(a.el) && POLAR.has(b.el) ? Math.min(d, 2.80) : d;
  }

  /* ---------------- system assembly ---------------- */

  /* Collect the atoms, the restraints between them and the exclusion sets.
     `residues` are the ones that may move; everything else in the model is context. */
  function build(model, residues, opts) {
    opts = Object.assign({}, DEFAULTS, opts);
    const mobile = new Set(residues.filter(r => opts.moveLigands || !r.het).map(r => r.uid));

    const atoms = [];           // { atom, res, rid, el, free }
    const X = [];               // flat coordinates, 3 per atom
    const slot = [];            // index into the free vector, or -1
    const ridOf = new Map();    // residue uid -> rid
    const atomIdx = [];         // rid -> { atomName: index into atoms }
    let nFree = 0;

    const addAtom = (atom, res, rid, free) => {
      const k = atoms.length;
      atoms.push({ atom, res, rid, el: atom.el, free });
      X.push(atom.p[0], atom.p[1], atom.p[2]);
      slot.push(free ? nFree++ : -1);
      if (!atomIdx[rid]) atomIdx[rid] = {};
      atomIdx[rid][atom.name] = k;
      return k;
    };

    // Residues in chain order, so that neighbours are adjacent rids.
    const chainOf = new Map();
    const order = [];
    for (const chain of model.chains) for (const res of chain.residues) { order.push(res); chainOf.set(res.uid, chain); }
    for (const res of model.ligands) order.push(res);

    // Which residues take part at all: the mobile ones, their chain neighbours (as anchors), and
    // any residue with an atom close to the mobile set.
    const mobilePts = [];
    for (const res of order) if (mobile.has(res.uid)) for (const a of res.atoms) mobilePts.push(a.p);
    if (!mobilePts.length) return null;
    const near = nearbyFilter(mobilePts, opts.contextCut);

    const taking = new Set();
    for (const chain of model.chains) {
      chain.residues.forEach((res, i) => {
        if (!mobile.has(res.uid)) return;
        taking.add(res.uid);
        for (const j of [i - 1, i + 1]) {
          const nb = chain.residues[j];
          if (nb) taking.add(nb.uid);
        }
      });
    }
    for (const res of order) {
      if (taking.has(res.uid) || mobile.has(res.uid)) continue;
      if (res.atoms.some(a => near(a.p))) taking.add(res.uid);
    }

    let rid = 0;
    const rids = [];
    for (const res of order) {
      if (!taking.has(res.uid) && !mobile.has(res.uid)) continue;
      const free = mobile.has(res.uid);
      const myRid = rid++;
      ridOf.set(res.uid, myRid);
      rids.push({ res, rid: myRid, free, chain: chainOf.get(res.uid) || null });
      for (const a of res.atoms) addAtom(a, res, myRid, free);
    }

    const dist = [];            // { a, b, d0, w }   harmonic both ways
    const push1 = [];           // { a, b, d0, w }   one-sided: only penalised when d < d0
    const chiral = [];          // { n, ca, c, cb, v0, w }
    const excl = new Set();     // "a,b" pairs skipped by the soft-sphere term
    const pairKey = (a, b) => a < b ? a + ',' + b : b + ',' + a;
    const addDist = (a, b, d0, w) => { if (a !== undefined && b !== undefined) dist.push({ a, b, d0, w }); };
    const exclude = (a, b) => { if (a !== undefined && b !== undefined) excl.add(pairKey(a, b)); };

    const P = i => [X[3 * i], X[3 * i + 1], X[3 * i + 2]];
    const d2 = (a, b) => { const p = P(a), q = P(b); return Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]); };

    const breaks = [];
    for (const { res, rid: ri, chain } of rids) {
      const A = atomIdx[ri];
      const has = n => A[n] !== undefined;
      const isAA = !res.het && has('N') && has('CA') && has('C');
      const local = new Set();
      const seenPair = (a, b) => { const k = pairKey(a, b); if (local.has(k)) return true; local.add(k); return false; };

      // Everything inside a residue is excluded from the soft-sphere term.
      const all = Object.values(A);
      for (let i = 0; i < all.length; i++) for (let j = i + 1; j < all.length; j++) exclude(all[i], all[j]);

      if (!isAA) {
        // A ligand or other non-standard group: hold its shape, nothing else.
        for (let i = 0; i < all.length; i++) for (let j = i + 1; j < all.length; j++) {
          if (!seenPair(all[i], all[j])) addDist(all[i], all[j], d2(all[i], all[j]), W.rigid);
        }
        continue;
      }

      // Intra-residue backbone.
      seenPair(A.N, A.CA); seenPair(A.CA, A.C); seenPair(A.N, A.C);
      addDist(A.N, A.CA, T.N_CA, W.bond);
      addDist(A.CA, A.C, T.CA_C, W.bond);
      addDist(A.C, A.O, T.C_O, W.bond);
      addDist(A.N, A.C, T.N_C, W.angle);
      addDist(A.CA, A.O, T.CA_O, W.angle);
      if (has('CB')) {
        seenPair(A.CA, A.CB); seenPair(A.N, A.CB); seenPair(A.C, A.CB);
        addDist(A.CA, A.CB, T.CA_CB, W.bond);
        addDist(A.N, A.CB, T.N_CB, W.angle);
        addDist(A.C, A.CB, T.C_CB, W.angle);
        chiral.push({ n: A.N, ca: A.CA, c: A.C, cb: A.CB, v0: T.CHIRAL_V, w: W.chiral });
      }
      if (has('O')) { seenPair(A.C, A.O); seenPair(A.CA, A.O); }
      if (has('OXT')) {
        seenPair(A.C, A.OXT); seenPair(A.CA, A.OXT); seenPair(A.O, A.OXT);
        addDist(A.C, A.OXT, T.C_O, W.bond); addDist(A.CA, A.OXT, T.CA_O, W.angle);
        // Without this the two carboxylate oxygens are each only tied to C and CA, which leaves
        // them free to swing around that axis and even sit on top of each other.
        addDist(A.O, A.OXT, T.O_OXT, W.angle);
      }

      // Side chains: hold their shape and their rotamer, anchored on the backbone so they travel
      // with it. Targets come from a reference residue at the same chi angles when the side chain
      // is complete, and from the present coordinates otherwise (a truncated or non-standard side
      // chain, which is then anchored on CA and CB only).
      const side = res.atoms.map(a => a.name).filter(n => !BB_ORDER.includes(n) && n !== 'OXT');
      if (opts.rigidSide && side.length) {
        const ref = (opts.sideMode !== 'keep' && sideChainComplete(res))
          ? referenceResidue(res.name, Geom.chiAngles(res) || []) : null;
        const anchors = (ref ? ['N', 'CA', 'C', 'CB'] : ['CA', 'CB']).filter(has);
        const names = side.concat(anchors);
        for (let i = 0; i < names.length; i++) for (let j = i + 1; j < names.length; j++) {
          const na = names[i], nb2 = names[j], a = A[na], b = A[nb2];
          if (seenPair(a, b)) continue;
          const target = ref && ref[na] && ref[nb2] ? V.dist(ref[na], ref[nb2]) : d2(a, b);
          addDist(a, b, target, W.rigid);
        }
      }

      // The peptide unit to the next residue.
      if (!chain) continue;
      const idx = chain.residues.indexOf(res);
      const nx = chain.residues[idx + 1];
      if (!nx || ridOf.get(nx.uid) === undefined) continue;
      const B = atomIdx[ridOf.get(nx.uid)];
      if (!has('C') || B.N === undefined) continue;
      const gap = d2(A.C, B.N);
      const broken = gap > Model.BOND_MAX;
      if (broken) {
        breaks.push({ res, next: nx, gap, closing: gap <= opts.maxGap });
        if (gap > opts.maxGap) continue;   // too far: leave it to Connect
      }
      addDist(A.C, B.N, T.C_N, W.bond);
      addDist(A.CA, B.N, T.CA_N1, W.angle);
      addDist(A.O, B.N, T.O_N1, W.angle);
      addDist(A.C, B.CA, T.C_CA1, W.angle);
      // Keep whichever peptide isomer it already has (cis prolines stay cis).
      const om = (A.CA !== undefined && B.CA !== undefined && !broken)
        ? Math.abs(V.dihedral(P(A.CA), P(A.C), P(B.N), P(B.CA))) : 180;
      const cis = om < 90;
      addDist(A.CA, B.CA, cis ? T.CA_CA1_cis : T.CA_CA1_trans, W.peptide);
      addDist(A.O, B.CA, cis ? T.O_CA1_cis : T.O_CA1_trans, W.peptide);

      // Keep the conformation. Nothing above restrains phi or psi — deliberately, since the design
      // is the user's — but that leaves them free to drift under the soft-sphere term, and a
      // terminal residue, whose phi almost nothing else touches, can swing right round. These two
      // 1-4 distances pin psi(i) and phi(i+1) weakly. Their targets come from an ideal two-residue
      // fragment built at the torsions this pair has now, so they agree with the bond and angle
      // restraints instead of carrying over whatever distortion is in the present coordinates.
      if (!broken && A.N !== undefined && B.C !== undefined && A.CA !== undefined && B.CA !== undefined) {
        const psiI = V.dihedral(P(A.N), P(A.CA), P(A.C), P(B.N));
        const omega = V.dihedral(P(A.CA), P(A.C), P(B.N), P(B.CA));
        const phiI1 = V.wrap(V.dihedral(P(A.C), P(B.N), P(B.CA), P(B.C)));
        const frag = phi => Geom.backboneFromTorsions([{ phi: -60, psi: psiI, omega }, { phi, psi: -45, omega: 180 }]);
        const f = frag(phiI1);
        addDist(A.N, B.N, V.dist(f[0].N, f[1].N), W.conf);
        addDist(A.C, B.C, V.dist(f[0].C, f[1].C), W.conf);
        // Push phi out of the band around zero that no non-glycine residue can occupy. The
        // C(i)...C(i+1) distance grows with |phi| and carries no sign, so one target serves both
        // +40 and -40 and the minimiser leaves by whichever side is nearer. One-sided, so it only
        // ever pushes |phi| up and never drags a good torsion back down.
        if (nx.name !== 'GLY' && Math.abs(phiI1) < PHI_MIN) {
          const g2 = frag(PHI_MIN);
          push1.push({ a: A.C, b: B.C, d0: V.dist(g2[0].C, g2[1].C), w: W.rama });
        }
      }
      // Neighbouring residues never push on each other: their contacts are set by phi/psi, which
      // these restraints and the starting positions already decide.
      for (const a of Object.values(A)) for (const b of Object.values(B)) exclude(a, b);
    }

    // Links between residues that are not the peptide bond — disulfides above all, but also a
    // covalently attached ligand or a coordinated metal. Without this the soft-sphere term treats
    // an S-S bond at 2.03 A as a bad clash and pulls the bridge apart.
    //
    // Distance alone cannot tell a bond from an overlap: two badly placed residues can easily be
    // 1.8 A apart. So the test is chemical. Two standard residues are only ever joined by the
    // peptide bond (handled above) or a disulfide; anything else counts as a link only when a
    // hetero group — a ligand or a metal — is one of the partners.
    for (let i = 0; i < atoms.length; i++) {
      for (let j = i + 1; j < atoms.length; j++) {
        if (atoms[i].rid === atoms[j].rid) continue;
        if (!atoms[i].free && !atoms[j].free) continue;
        const d = d2(i, j);
        const ss = atoms[i].el === 'S' && atoms[j].el === 'S';
        const het = atoms[i].res.het || atoms[j].res.het;
        if (ss ? d > 2.6 : !(het && d <= 2.4)) continue;
        if (excl.has(pairKey(i, j))) continue;   // the peptide bond, already restrained
        addDist(i, j, ss ? 2.03 : d, W.bond);
        exclude(i, j);
        // and whatever is bonded to each end, so the term does not fight the link's own geometry
        for (const [p, q] of [[i, j], [j, i]]) {
          for (let k = 0; k < atoms.length; k++) {
            if (atoms[k].rid !== atoms[p].rid || k === p) continue;
            if (d2(k, p) < 1.95) exclude(k, q);
          }
        }
      }
    }

    const x0 = new Float64Array(X);
    return {
      opts, atoms, X: new Float64Array(X), x0, slot, nFree, dist, push1, chiral, excl, breaks, atomIdx, rids,
      nMobile: rids.filter(r => r.free).length,
    };
  }

  /* A test for "within cut of any of these points", on a grid. */
  function nearbyFilter(pts, cut) {
    const cell = Math.max(cut, 2), grid = new Map();
    const key = (x, y, z) => x + ',' + y + ',' + z;
    for (const p of pts) {
      const k = key(Math.floor(p[0] / cell), Math.floor(p[1] / cell), Math.floor(p[2] / cell));
      if (!grid.has(k)) grid.set(k, []);
      grid.get(k).push(p);
    }
    const c2 = cut * cut;
    return p => {
      const gx = Math.floor(p[0] / cell), gy = Math.floor(p[1] / cell), gz = Math.floor(p[2] / cell);
      for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
        const l = grid.get(key(gx + dx, gy + dy, gz + dz));
        if (l) for (const q of l) if (V.dist2(p, q) < c2) return true;
      }
      return false;
    };
  }

  /* Soft-sphere pairs, rebuilt from the current coordinates. Only pairs with at least one mobile
     atom matter; fixed-fixed pairs cannot change. */
  function nbPairs(sys) {
    const { atoms, X, excl, opts } = sys;
    const scale = opts.clashScale;
    const maxD = 2 * 1.9 * scale + 0.8;
    const cell = maxD, grid = new Map();
    const key = (x, y, z) => x + ',' + y + ',' + z;
    for (let i = 0; i < atoms.length; i++) {
      const k = key(Math.floor(X[3 * i] / cell), Math.floor(X[3 * i + 1] / cell), Math.floor(X[3 * i + 2] / cell));
      if (!grid.has(k)) grid.set(k, []);
      grid.get(k).push(i);
    }
    const out = [];
    const seen = new Set();
    for (let i = 0; i < atoms.length; i++) {
      const gx = Math.floor(X[3 * i] / cell), gy = Math.floor(X[3 * i + 1] / cell), gz = Math.floor(X[3 * i + 2] / cell);
      for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
        const l = grid.get(key(gx + dx, gy + dy, gz + dz));
        if (!l) continue;
        for (const j of l) {
          if (j <= i) continue;
          if (!atoms[i].free && !atoms[j].free) continue;
          const k = i + ',' + j;
          if (seen.has(k) || excl.has(k)) continue;
          const d0 = clashTarget(atoms[i], atoms[j], scale);
          const dx2 = X[3 * i] - X[3 * j], dy2 = X[3 * i + 1] - X[3 * j + 1], dz2 = X[3 * i + 2] - X[3 * j + 2];
          if (dx2 * dx2 + dy2 * dy2 + dz2 * dz2 > (d0 + 0.6) * (d0 + 0.6)) continue;
          seen.add(k);
          out.push({ a: i, b: j, d0, w: W.nb });
        }
      }
    }
    return out;
  }

  /* ---------------- energy and gradient ---------------- */

  function energy(sys, nb, g) {
    const { X, x0, dist, chiral, slot, opts } = sys;
    let E = 0;
    if (g) g.fill(0);
    const push = (i, gx, gy, gz) => {
      const s = slot[i];
      if (s < 0 || !g) return;
      g[3 * s] += gx; g[3 * s + 1] += gy; g[3 * s + 2] += gz;
    };

    const harmonic = (list, oneSided) => {
      for (const r of list) {
        const ax = 3 * r.a, bx = 3 * r.b;
        const dx = X[ax] - X[bx], dy = X[ax + 1] - X[bx + 1], dz = X[ax + 2] - X[bx + 2];
        const d = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-9;
        const diff = d - r.d0;
        if (oneSided && diff > 0) continue;
        E += r.w * diff * diff;
        const f = 2 * r.w * diff / d;
        push(r.a, f * dx, f * dy, f * dz);
        push(r.b, -f * dx, -f * dy, -f * dz);
      }
    };
    harmonic(dist, false);
    harmonic(sys.push1, true);
    harmonic(nb, true);

    for (const c of chiral) {
      const ca = [X[3 * c.ca], X[3 * c.ca + 1], X[3 * c.ca + 2]];
      const a = [X[3 * c.n] - ca[0], X[3 * c.n + 1] - ca[1], X[3 * c.n + 2] - ca[2]];
      const b = [X[3 * c.c] - ca[0], X[3 * c.c + 1] - ca[1], X[3 * c.c + 2] - ca[2]];
      const e = [X[3 * c.cb] - ca[0], X[3 * c.cb + 1] - ca[1], X[3 * c.cb + 2] - ca[2]];
      const bxe = V.cross(b, e), exa = V.cross(e, a), axb = V.cross(a, b);
      const vol = a[0] * bxe[0] + a[1] * bxe[1] + a[2] * bxe[2];
      const diff = vol - c.v0;
      E += c.w * diff * diff;
      const f = 2 * c.w * diff;
      push(c.n, f * bxe[0], f * bxe[1], f * bxe[2]);
      push(c.c, f * exa[0], f * exa[1], f * exa[2]);
      push(c.cb, f * axb[0], f * axb[1], f * axb[2]);
      push(c.ca, -f * (bxe[0] + exa[0] + axb[0]), -f * (bxe[1] + exa[1] + axb[1]), -f * (bxe[2] + exa[2] + axb[2]));
    }

    const wp = opts.hold;
    if (wp > 0) {
      for (let i = 0; i < slot.length; i++) {
        if (slot[i] < 0) continue;
        const dx = X[3 * i] - x0[3 * i], dy = X[3 * i + 1] - x0[3 * i + 1], dz = X[3 * i + 2] - x0[3 * i + 2];
        E += wp * (dx * dx + dy * dy + dz * dz);
        push(i, 2 * wp * dx, 2 * wp * dy, 2 * wp * dz);
      }
    }
    return E;
  }

  /* Copy the free coordinates out of / into the full coordinate array. */
  function gather(sys, x) { for (let i = 0; i < sys.slot.length; i++) { const s = sys.slot[i]; if (s >= 0) { x[3 * s] = sys.X[3 * i]; x[3 * s + 1] = sys.X[3 * i + 1]; x[3 * s + 2] = sys.X[3 * i + 2]; } } }
  function scatter(sys, x) { for (let i = 0; i < sys.slot.length; i++) { const s = sys.slot[i]; if (s >= 0) { sys.X[3 * i] = x[3 * s]; sys.X[3 * i + 1] = x[3 * s + 1]; sys.X[3 * i + 2] = x[3 * s + 2]; } } }

  /* ---------------- L-BFGS ---------------- */
  function minimize(sys, iterations) {
    const n = 3 * sys.nFree;
    if (!n) return { iterations: 0, E0: 0, E: 0 };
    const x = new Float64Array(n), g = new Float64Array(n), xNew = new Float64Array(n), gNew = new Float64Array(n);
    gather(sys, x);
    let nb = nbPairs(sys);
    const evaluate = (xv, gv) => { scatter(sys, xv); return energy(sys, nb, gv); };
    let E = evaluate(x, g);
    const E0 = E;
    const m = 8, S = [], Y = [], RHO = [];
    const d = new Float64Array(n);
    let xRef = new Float64Array(x);   // coordinates the pair list was built from
    let step = 1 / Math.max(1e-6, Math.sqrt(dot(g, g) / n));
    let it = 0;
    for (; it < iterations; it++) {
      const gn = Math.sqrt(dot(g, g));
      if (gn < 1e-6) break;
      // Two-loop recursion for the search direction.
      for (let i = 0; i < n; i++) d[i] = -g[i];
      const alpha = new Float64Array(S.length);
      for (let k = S.length - 1; k >= 0; k--) {
        alpha[k] = RHO[k] * dot(S[k], d);
        for (let i = 0; i < n; i++) d[i] -= alpha[k] * Y[k][i];
      }
      if (S.length) {
        const k = S.length - 1;
        const yy = dot(Y[k], Y[k]);
        const scale = yy > 1e-12 ? dot(S[k], Y[k]) / yy : 1;
        for (let i = 0; i < n; i++) d[i] *= scale;
      }
      for (let k = 0; k < S.length; k++) {
        const beta = RHO[k] * dot(Y[k], d);
        for (let i = 0; i < n; i++) d[i] += (alpha[k] - beta) * S[k][i];
      }
      let slope = dot(g, d);
      if (slope >= 0) { for (let i = 0; i < n; i++) d[i] = -g[i]; slope = -dot(g, g); S.length = Y.length = RHO.length = 0; }

      // Armijo backtracking, with the step capped so no atom moves more than 0.3 A at once.
      let dmax = 0;
      for (let i = 0; i < n; i += 3) dmax = Math.max(dmax, Math.hypot(d[i], d[i + 1], d[i + 2]));
      let t = Math.min(S.length ? 1 : step, 0.3 / Math.max(dmax, 1e-9));
      let ok = false, ENew = E;
      for (let ls = 0; ls < 20; ls++) {
        for (let i = 0; i < n; i++) xNew[i] = x[i] + t * d[i];
        ENew = evaluate(xNew, gNew);
        if (ENew <= E + 1e-4 * t * slope) { ok = true; break; }
        t *= 0.5;
      }
      if (!ok) break;
      const s = new Float64Array(n), y = new Float64Array(n);
      for (let i = 0; i < n; i++) { s[i] = xNew[i] - x[i]; y[i] = gNew[i] - g[i]; }
      const sy = dot(s, y);
      if (sy > 1e-12) {
        S.push(s); Y.push(y); RHO.push(1 / sy);
        if (S.length > m) { S.shift(); Y.shift(); RHO.shift(); }
      }
      x.set(xNew); g.set(gNew);
      const dE = E - ENew;
      E = ENew; step = t;
      // The pair list is built with a 0.6 A buffer, so it only goes stale once an atom has moved
      // half of that. Rebuilding changes the energy surface, so the curvature history goes with it
      // — which is why it is not done on a fixed schedule.
      let moved = 0;
      for (let i = 0; i < n; i += 3) moved = Math.max(moved, Math.hypot(x[i] - xRef[i], x[i + 1] - xRef[i + 1], x[i + 2] - xRef[i + 2]));
      if (moved > 0.3) {
        scatter(sys, x); nb = nbPairs(sys); E = evaluate(x, g);
        xRef = new Float64Array(x);
        S.length = Y.length = RHO.length = 0;
      }
      if (dE < 1e-7 * Math.max(1, Math.abs(E))) break;
    }
    scatter(sys, x);
    return { iterations: it, E0, E };
  }

  function dot(a, b) { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; }

  /* ---------------- measurement ---------------- */

  /* Deviations of the current coordinates from the restraint targets, plus the overlaps left. */
  function measure(sys) {
    const { X, dist, chiral, opts } = sys;
    const d = (a, b) => Math.hypot(X[3 * a] - X[3 * b], X[3 * a + 1] - X[3 * b + 1], X[3 * a + 2] - X[3 * b + 2]);
    let bs = 0, bn = 0, bmax = 0, as = 0, an = 0, amax = 0, chi = 0;
    for (const r of dist) {
      if (r.w === W.rigid || r.w === W.conf) continue;
      const e = Math.abs(d(r.a, r.b) - r.d0);
      if (r.w === W.bond) { bs += e * e; bn++; bmax = Math.max(bmax, e); }
      else { as += e * e; an++; amax = Math.max(amax, e); }
    }
    for (const c of chiral) {
      const ca = [X[3 * c.ca], X[3 * c.ca + 1], X[3 * c.ca + 2]];
      const a = V.sub([X[3 * c.n], X[3 * c.n + 1], X[3 * c.n + 2]], ca);
      const b = V.sub([X[3 * c.c], X[3 * c.c + 1], X[3 * c.c + 2]], ca);
      const e = V.sub([X[3 * c.cb], X[3 * c.cb + 1], X[3 * c.cb + 2]], ca);
      if (V.dot(a, V.cross(b, e)) < 0) chi++;
    }
    let clashes = 0, worst = 0;
    for (const p of nbPairs(sys)) {
      const dd = d(p.a, p.b);
      const over = clashTarget(sys.atoms[p.a], sys.atoms[p.b], opts.clashScale) - dd;
      if (over > 0.25) { clashes++; worst = Math.max(worst, over); }
    }
    return {
      bondRms: bn ? Math.sqrt(bs / bn) : 0, bondMax: bmax,
      angleRms: an ? Math.sqrt(as / an) : 0, angleMax: amax,
      wrongChirality: chi, clashes, worstOverlap: worst, rama: ramaOutliers(sys),
    };
  }

  /* Mobile non-glycine residues whose phi is where no residue can go. Read from the working
     coordinates, not the model, so it can be called mid-refinement. */
  function ramaOutliers(sys) {
    const { X, rids, atomIdx } = sys;
    const P = i => [X[3 * i], X[3 * i + 1], X[3 * i + 2]];
    const out = [];
    for (let k = 1; k < rids.length; k++) {
      const cur = rids[k], prev = rids[k - 1];
      if (!cur.free || cur.res.name === 'GLY' || cur.res.het) continue;
      if (!cur.chain || cur.chain !== prev.chain) continue;
      if (cur.chain.residues.indexOf(cur.res) !== cur.chain.residues.indexOf(prev.res) + 1) continue;
      const A = atomIdx[prev.rid], B = atomIdx[cur.rid];
      if (A.C === undefined || B.N === undefined || B.CA === undefined || B.C === undefined) continue;
      if (V.dist(P(A.C), P(B.N)) > Model.BOND_MAX) continue;
      const phi = V.wrap(V.dihedral(P(A.C), P(B.N), P(B.CA), P(B.C)));
      if (Math.abs(phi) < PHI_MIN - 5 || phi > PHI_MAX_POS + 5) out.push({ res: cur.res, phi });
    }
    return out;
  }

  function shifts(sys) {
    let sum = 0, n = 0, max = 0;
    for (let i = 0; i < sys.slot.length; i++) {
      if (sys.slot[i] < 0) continue;
      const dx = sys.X[3 * i] - sys.x0[3 * i], dy = sys.X[3 * i + 1] - sys.x0[3 * i + 1], dz = sys.X[3 * i + 2] - sys.x0[3 * i + 2];
      const d2v = dx * dx + dy * dy + dz * dz;
      sum += d2v; n++; max = Math.max(max, d2v);
    }
    return { rmsd: n ? Math.sqrt(sum / n) : 0, maxShift: Math.sqrt(max) };
  }

  /* ---------------- entry points ---------------- */

  /* What is wrong with these residues, without changing anything. */
  function analyse(model, residues, opts) {
    const sys = build(model, residues, opts);
    if (!sys) return null;
    const m = measure(sys);
    m.breaks = sys.breaks;
    m.residues = sys.nMobile;
    return m;
  }

  /* Regularise `residues` in place. Returns a report, or null when there is nothing to do. */
  function run(model, residues, opts) {
    opts = Object.assign({}, DEFAULTS, opts);
    const movable = residues.filter(r => opts.moveLigands || !r.het);
    if (!movable.length) return null;

    // Rebuilding a complete side chain from its current chi angles clears any distortion in it
    // before the backbone is regularised; incomplete ones (backbone-only models) are left alone.
    let rebuilt = 0;
    if (opts.sideMode !== 'keep') {
      for (const r of movable) {
        if (r.het || r.name === 'GLY' || !sideChainComplete(r) || !Geom.pos(r, 'CA')) continue;
        if (Geom.buildSideChain(r, r.name, Geom.chiAngles(r))) rebuilt++;
      }
    }

    const start = new Map();
    for (const r of movable) for (const a of r.atoms) start.set(a, a.p.slice());

    // Several passes. The conformational and phi restraints are derived from the coordinates the
    // pass starts with, so rebuilding them on the result of the last pass lets strain that was
    // locked in by the starting geometry — a phi pushed into the forbidden band by closing a gap,
    // say — work its way out. Each pass after the first moves very little, so this is cheap.
    let before = null, after = null, sys = null, first = null, iterations = 0;
    let pack = null;
    // Repacking happens before the first minimisation, so the "before" numbers have to be taken
    // here or they would describe the already-repacked model rather than what the user had.
    if (opts.sideMode === 'repack') {
      const s0 = build(model, movable, opts);
      if (!s0) return null;
      before = measure(s0);
      first = s0;
    }
    for (let pass = 0; pass < (opts.passes || 3); pass++) {
      // Repack before each pass: the search needs a backbone to work against, and the backbone
      // the last pass produced is a better one than the backbone this started with.
      if (opts.sideMode === 'repack' && window.Pack) {
        const pr = Pack.run(model, movable, { clashScale: opts.clashScale });
        if (pr) pack = pack
          ? Object.assign({}, pr, { changed: pack.changed + pr.changed, energyBefore: pack.energyBefore })
          : pr;
      }
      sys = build(model, movable, opts);
      if (!sys) return null;
      const m0 = measure(sys);
      if (!before) { before = m0; }
      if (!first) first = sys;
      const res = minimize(sys, opts.iterations);
      iterations += res.iterations;
      after = measure(sys);
      for (let i = 0; i < sys.atoms.length; i++) {
        if (sys.slot[i] < 0) continue;
        sys.atoms[i].atom.p = [sys.X[3 * i], sys.X[3 * i + 1], sys.X[3 * i + 2]];
      }
      const moved = shifts(sys).rmsd;
      if (pass && moved < 0.01 && after.rama.length >= m0.rama.length) break;
    }

    const closed = [];
    for (const b of first.breaks) {
      if (!b.closing) continue;
      const C = Geom.pos(b.res, 'C'), N = Geom.pos(b.next, 'N');
      closed.push({ from: b.res, to: b.next, was: b.gap, now: V.dist(C, N) });
    }
    // Displacement measured against where things were before the first pass.
    let sum = 0, n = 0, max = 0;
    for (const [atom, p0] of start) {
      const d = V.dist2(atom.p, p0);
      sum += d; n++; max = Math.max(max, d);
    }
    return { before, after, rebuilt, pack, residues: sys.nMobile, iterations,
      closed, skippedBreaks: first.breaks.filter(b => !b.closing),
      rmsd: n ? Math.sqrt(sum / n) : 0, maxShift: Math.sqrt(max) };
  }

  window.Refine = { run, analyse, build, measure, minimize, nbPairs, DEFAULTS, T, W, clashTarget,
    /* exposed for the test page's finite-difference gradient check */
    _energy: energy, _gather: gather, _scatter: scatter };
})();
