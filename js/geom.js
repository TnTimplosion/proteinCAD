/* Protein geometry: ideal bond geometry, NeRF placement, backbone from torsions, side chains from
   internal coordinates, screw axes of regular segments.

   A "residue" here is the model's plain object { uid, name, num, icode, ss, atoms:[{name, el, p}] }.
   This module never touches the model as a whole; model.js and builders.js do. */
(function () {
  'use strict';
  const { D2R } = V;

  /* Engh & Huber style ideal backbone. */
  const BB = {
    N_CA: 1.458, CA_C: 1.525, C_N: 1.329, C_O: 1.231, CA_CB: 1.530,
    N_CA_C: 111.2, CA_C_N: 116.2, C_N_CA: 121.7, CA_C_O: 120.5, N_CA_CB: 110.3,
    // Torsion C-N-CA-CB for an L-amino acid (checked against the ProteinMPNN virtual-CB formula).
    C_N_CA_CB: -122.5,
  };

  /* Named (phi, psi) presets. */
  const CONF = {
    alpha: { label: 'α-helix', phi: -57.8, psi: -47.0, ss: 'H' },
    helix310: { label: '3₁₀-helix', phi: -49.0, psi: -26.0, ss: 'H' },
    pi: { label: 'π-helix', phi: -57.0, psi: -70.0, ss: 'H' },
    antiparallel: { label: 'β (antiparallel)', phi: -139.0, psi: 135.0, ss: 'E' },
    parallel: { label: 'β (parallel)', phi: -119.0, psi: 113.0, ss: 'E' },
    ppii: { label: 'Polyproline II', phi: -75.0, psi: 145.0, ss: 'L' },
  };

  /* NeRF: place d so that |cd| = bond, angle b-c-d = ang, dihedral a-b-c-d = tor (degrees). */
  function place(a, b, c, bond, ang, tor) {
    ang *= D2R; tor *= D2R;
    const bc = V.norm(V.sub(c, b));
    const n = V.norm(V.cross(V.sub(b, a), bc));
    const m = V.cross(n, bc);
    const x = -bond * Math.cos(ang), y = bond * Math.sin(ang) * Math.cos(tor), z = bond * Math.sin(ang) * Math.sin(tor);
    return [c[0] + x * bc[0] + y * m[0] + z * n[0],
      c[1] + x * bc[1] + y * m[1] + z * n[1],
      c[2] + x * bc[2] + y * m[2] + z * n[2]];
  }

  /* ---- Residue tables ---- */
  const THREE_TO_ONE = {
    ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G', HIS: 'H', ILE: 'I',
    LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S', THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V',
  };
  const ONE_TO_THREE = Object.fromEntries(Object.entries(THREE_TO_ONE).map(([k, v]) => [v, k]));
  const AA_NAMES = Object.keys(THREE_TO_ONE).sort();
  /* Residue names folded into a standard one on import. */
  const ALIASES = {
    MSE: 'MET', HID: 'HIS', HIE: 'HIS', HIP: 'HIS', HSD: 'HIS', HSE: 'HIS', HSP: 'HIS', CYX: 'CYS',
    CYM: 'CYS', ASH: 'ASP', GLH: 'GLU', LYN: 'LYS', ARN: 'ARG',
  };
  const BACKBONE = new Set(['N', 'CA', 'C', 'O', 'OXT']);

  /* Side-chain internal coordinates: [atom, element, ref1, ref2, ref3, bond, angle, torsion]
     where torsion is a number or 'chiK' / 'chiK+off'. Default chi values are common rotamers
     (Lovell et al. 2000 modal values, rounded). */
  const SC = {
    ALA: { chi: [], atoms: [] },
    GLY: { chi: [], atoms: [] },
    SER: { chi: [64], atoms: [['OG', 'O', 'N', 'CA', 'CB', 1.417, 110.8, 'chi1']] },
    CYS: { chi: [-65], atoms: [['SG', 'S', 'N', 'CA', 'CB', 1.808, 113.8, 'chi1']] },
    THR: {
      chi: [62], atoms: [
        ['OG1', 'O', 'N', 'CA', 'CB', 1.433, 109.2, 'chi1'],
        ['CG2', 'C', 'N', 'CA', 'CB', 1.521, 111.5, 'chi1-120']],
    },
    VAL: {
      chi: [175], atoms: [
        ['CG1', 'C', 'N', 'CA', 'CB', 1.527, 110.7, 'chi1'],
        ['CG2', 'C', 'N', 'CA', 'CB', 1.527, 110.4, 'chi1+120']],
    },
    ILE: {
      chi: [-65, 170], atoms: [
        ['CG1', 'C', 'N', 'CA', 'CB', 1.527, 110.7, 'chi1'],
        ['CG2', 'C', 'N', 'CA', 'CB', 1.527, 110.4, 'chi1-120'],
        ['CD1', 'C', 'CA', 'CB', 'CG1', 1.520, 113.97, 'chi2']],
    },
    LEU: {
      chi: [-65, 175], atoms: [
        ['CG', 'C', 'N', 'CA', 'CB', 1.530, 116.1, 'chi1'],
        ['CD1', 'C', 'CA', 'CB', 'CG', 1.524, 110.3, 'chi2'],
        ['CD2', 'C', 'CA', 'CB', 'CG', 1.525, 110.6, 'chi2+120']],
    },
    MET: {
      chi: [-65, 180, -70], atoms: [
        ['CG', 'C', 'N', 'CA', 'CB', 1.520, 113.7, 'chi1'],
        ['SD', 'S', 'CA', 'CB', 'CG', 1.810, 112.7, 'chi2'],
        ['CE', 'C', 'CB', 'CG', 'SD', 1.790, 100.6, 'chi3']],
    },
    PHE: {
      chi: [-65, 90], atoms: [
        ['CG', 'C', 'N', 'CA', 'CB', 1.502, 113.9, 'chi1'],
        ['CD1', 'C', 'CA', 'CB', 'CG', 1.389, 120.0, 'chi2'],
        ['CD2', 'C', 'CA', 'CB', 'CG', 1.389, 120.0, 'chi2+180'],
        ['CE1', 'C', 'CB', 'CG', 'CD1', 1.389, 120.0, 180],
        ['CE2', 'C', 'CB', 'CG', 'CD2', 1.389, 120.0, 180],
        ['CZ', 'C', 'CG', 'CD1', 'CE1', 1.389, 120.0, 0]],
    },
    TYR: {
      chi: [-65, 90], atoms: [
        ['CG', 'C', 'N', 'CA', 'CB', 1.512, 113.9, 'chi1'],
        ['CD1', 'C', 'CA', 'CB', 'CG', 1.389, 120.8, 'chi2'],
        ['CD2', 'C', 'CA', 'CB', 'CG', 1.389, 120.8, 'chi2+180'],
        ['CE1', 'C', 'CB', 'CG', 'CD1', 1.389, 121.2, 180],
        ['CE2', 'C', 'CB', 'CG', 'CD2', 1.389, 121.2, 180],
        ['CZ', 'C', 'CG', 'CD1', 'CE1', 1.378, 119.6, 0],
        ['OH', 'O', 'CD1', 'CE1', 'CZ', 1.376, 119.9, 180]],
    },
    TRP: {
      chi: [-65, 95], atoms: [
        ['CG', 'C', 'N', 'CA', 'CB', 1.498, 114.1, 'chi1'],
        ['CD1', 'C', 'CA', 'CB', 'CG', 1.365, 127.1, 'chi2'],
        ['CD2', 'C', 'CA', 'CB', 'CG', 1.433, 126.7, 'chi2+180'],
        ['NE1', 'N', 'CB', 'CG', 'CD1', 1.374, 110.2, 180],
        ['CE2', 'C', 'CB', 'CG', 'CD2', 1.409, 107.2, 180],
        ['CE3', 'C', 'CB', 'CG', 'CD2', 1.398, 133.9, 0],
        ['CZ2', 'C', 'CG', 'CD2', 'CE2', 1.394, 122.4, 180],
        ['CZ3', 'C', 'CG', 'CD2', 'CE3', 1.382, 118.7, 180],
        ['CH2', 'C', 'CD2', 'CE2', 'CZ2', 1.368, 117.5, 0]],
    },
    HIS: {
      chi: [-65, -70], atoms: [
        ['CG', 'C', 'N', 'CA', 'CB', 1.497, 113.7, 'chi1'],
        ['ND1', 'N', 'CA', 'CB', 'CG', 1.378, 122.7, 'chi2'],
        ['CD2', 'C', 'CA', 'CB', 'CG', 1.354, 131.0, 'chi2+180'],
        ['CE1', 'C', 'CB', 'CG', 'ND1', 1.321, 109.0, 180],
        ['NE2', 'N', 'CB', 'CG', 'CD2', 1.374, 107.0, 180]],
    },
    ASP: {
      chi: [-70, -15], atoms: [
        ['CG', 'C', 'N', 'CA', 'CB', 1.516, 112.6, 'chi1'],
        ['OD1', 'O', 'CA', 'CB', 'CG', 1.249, 118.4, 'chi2'],
        ['OD2', 'O', 'CA', 'CB', 'CG', 1.249, 118.4, 'chi2+180']],
    },
    ASN: {
      chi: [-65, -20], atoms: [
        ['CG', 'C', 'N', 'CA', 'CB', 1.516, 112.6, 'chi1'],
        ['OD1', 'O', 'CA', 'CB', 'CG', 1.231, 120.8, 'chi2'],
        ['ND2', 'N', 'CA', 'CB', 'CG', 1.328, 116.4, 'chi2+180']],
    },
    GLU: {
      chi: [-65, 180, -10], atoms: [
        ['CG', 'C', 'N', 'CA', 'CB', 1.520, 113.8, 'chi1'],
        ['CD', 'C', 'CA', 'CB', 'CG', 1.516, 112.6, 'chi2'],
        ['OE1', 'O', 'CB', 'CG', 'CD', 1.249, 118.4, 'chi3'],
        ['OE2', 'O', 'CB', 'CG', 'CD', 1.249, 118.4, 'chi3+180']],
    },
    GLN: {
      chi: [-65, 180, -25], atoms: [
        ['CG', 'C', 'N', 'CA', 'CB', 1.520, 113.8, 'chi1'],
        ['CD', 'C', 'CA', 'CB', 'CG', 1.516, 112.6, 'chi2'],
        ['OE1', 'O', 'CB', 'CG', 'CD', 1.231, 120.8, 'chi3'],
        ['NE2', 'N', 'CB', 'CG', 'CD', 1.328, 116.4, 'chi3+180']],
    },
    LYS: {
      chi: [-65, 180, 180, 180], atoms: [
        ['CG', 'C', 'N', 'CA', 'CB', 1.520, 113.8, 'chi1'],
        ['CD', 'C', 'CA', 'CB', 'CG', 1.520, 111.3, 'chi2'],
        ['CE', 'C', 'CB', 'CG', 'CD', 1.520, 111.3, 'chi3'],
        ['NZ', 'N', 'CG', 'CD', 'CE', 1.489, 111.9, 'chi4']],
    },
    ARG: {
      chi: [-65, 180, 180, 180], atoms: [
        ['CG', 'C', 'N', 'CA', 'CB', 1.520, 113.8, 'chi1'],
        ['CD', 'C', 'CA', 'CB', 'CG', 1.520, 111.3, 'chi2'],
        ['NE', 'N', 'CB', 'CG', 'CD', 1.460, 112.0, 'chi3'],
        ['CZ', 'C', 'CG', 'CD', 'NE', 1.329, 124.2, 'chi4'],
        ['NH1', 'N', 'CD', 'NE', 'CZ', 1.326, 120.0, 0],
        ['NH2', 'N', 'CD', 'NE', 'CZ', 1.326, 120.0, 180]],
    },
    PRO: {
      // Exo pucker; with the narrower N-CA-CB angle below, CD lands 1.47 A from N.
      chi: [-25, 34], atoms: [
        ['CG', 'C', 'N', 'CA', 'CB', 1.495, 104.5, 'chi1'],
        ['CD', 'C', 'CA', 'CB', 'CG', 1.502, 105.5, 'chi2']],
    },
  };

  function getAtom(res, name) {
    for (const a of res.atoms) if (a.name === name) return a;
    return null;
  }
  function pos(res, name) { const a = getAtom(res, name); return a ? a.p : null; }

  function torsionValue(spec, chi) {
    if (typeof spec === 'number') return spec;
    const m = /^chi(\d)([+-]\d+)?$/.exec(spec);
    return chi[+m[1] - 1] + (m[2] ? +m[2] : 0);
  }

  /* Ideal CB from the residue's own N, CA, C. Proline's ring pulls N-CA-CB in to about 103.5 deg. */
  function idealCB(N, CA, C, pro) {
    return pro ? place(C, N, CA, BB.CA_CB, 103.5, -120.0) : place(C, N, CA, BB.CA_CB, BB.N_CA_CB, BB.C_N_CA_CB);
  }

  /* Rebuild the side chain of res as type `name` (three-letter). Keeps N, CA, C, O (and OXT).
     chi: optional array overriding the default rotamer. Returns false if the backbone is incomplete. */
  function buildSideChain(res, name, chi) {
    const N = pos(res, 'N'), CA = pos(res, 'CA'), C = pos(res, 'C');
    if (!N || !CA || !C) return false;
    const lib = SC[name];
    res.atoms = res.atoms.filter(a => BACKBONE.has(a.name));
    res.name = name;
    if (!lib || name === 'GLY') return true;
    chi = chi || lib.chi;
    const placed = { N, CA, C, CB: idealCB(N, CA, C, name === 'PRO') };
    res.atoms.push({ name: 'CB', el: 'C', p: placed.CB });
    for (const [an, el, r1, r2, r3, b, ang, tor] of lib.atoms) {
      placed[an] = place(placed[r1], placed[r2], placed[r3], b, ang, torsionValue(tor, chi));
      res.atoms.push({ name: an, el, p: placed[an] });
    }
    return true;
  }

  /* Current chi angles of a residue (for keeping a rotamer across edits). */
  function chiAngles(res) {
    const lib = SC[res.name];
    if (!lib) return null;
    const out = [];
    const defs = {};
    for (const [an, , r1, r2, r3, , , tor] of lib.atoms) {
      if (typeof tor === 'string' && /^chi\d$/.test(tor)) defs[tor] = [r1, r2, r3, an];
    }
    for (let k = 1; k <= lib.chi.length; k++) {
      const d = defs['chi' + k]; if (!d) { out.push(lib.chi[k - 1]); continue; }
      const ps = d.map(n => n === 'CB' ? pos(res, 'CB') : pos(res, n));
      out.push(ps.every(Boolean) ? V.dihedral(...ps) : lib.chi[k - 1]);
    }
    return out;
  }

  /* Carbonyl O of residue i from its own N, CA, C and the next residue's N (in the peptide plane,
     anti to the bisector), or from psi when there is no next residue. */
  function carbonylO(N, CA, C, Nnext, psi) {
    if (Nnext && V.dist(C, Nnext) < 2.0) {
      const b = V.norm(V.add(V.norm(V.sub(C, CA)), V.norm(V.sub(C, Nnext))));
      return V.addScaled(C, b, BB.C_O);
    }
    return place(N, CA, C, BB.C_O, BB.CA_C_O, (psi === undefined ? 135 : psi) + 180);
  }

  /* Build an ideal backbone from per-residue torsions [{phi, psi, omega}], in a canonical frame
     (first N at the origin). Returns [{N, CA, C}] (no O). */
  function backboneFromTorsions(tors) {
    const N0 = [0, 0, 0], CA0 = [BB.N_CA, 0, 0];
    const a = BB.N_CA_C * D2R;
    const C0 = [CA0[0] - BB.CA_C * Math.cos(a), BB.CA_C * Math.sin(a), 0];
    const out = [{ N: N0, CA: CA0, C: C0 }];
    for (let i = 1; i < tors.length; i++) {
      const p = out[i - 1];
      const N = place(p.N, p.CA, p.C, BB.C_N, BB.CA_C_N, tors[i - 1].psi);
      const CA = place(p.CA, p.C, N, BB.N_CA, BB.C_N_CA, tors[i - 1].omega ?? 180);
      const C = place(p.C, N, CA, BB.CA_C, BB.N_CA_C, tors[i].phi);
      out.push({ N, CA, C });
    }
    return out;
  }

  /* Next residue's backbone after (N, CA, C) given this residue's psi, the omega and the next phi. */
  function extendForward(N, CA, C, psi, phiNext, omega) {
    const N1 = place(N, CA, C, BB.C_N, BB.CA_C_N, psi);
    const CA1 = place(CA, C, N1, BB.N_CA, BB.C_N_CA, omega ?? 180);
    const C1 = place(C, N1, CA1, BB.CA_C, BB.N_CA_C, phiNext);
    return { N: N1, CA: CA1, C: C1 };
  }

  /* Previous residue's backbone before (N, CA, C) given this residue's phi, omega and the previous psi. */
  function extendBackward(N, CA, C, phi, psiPrev, omega) {
    const Cp = place(C, CA, N, BB.C_N, BB.C_N_CA, phi);
    const CAp = place(CA, N, Cp, BB.CA_C, BB.CA_C_N, omega ?? 180);
    const Np = place(N, Cp, CAp, BB.N_CA, BB.N_CA_C, psiPrev);
    return { N: Np, CA: CAp, C: Cp };
  }

  /* Turn [{N, CA, C}] into model residues (atoms N, CA, C, O, side chain of `name`). */
  function residuesFromBackbone(bb, name, ss, psiLast) {
    const out = [];
    for (let i = 0; i < bb.length; i++) {
      const { N, CA, C } = bb[i];
      const O = carbonylO(N, CA, C, i + 1 < bb.length ? bb[i + 1].N : null, psiLast);
      const res = {
        name, num: i + 1, icode: '', ss: ss || 'L',
        atoms: [{ name: 'N', el: 'N', p: N }, { name: 'CA', el: 'C', p: CA }, { name: 'C', el: 'C', p: C }, { name: 'O', el: 'O', p: O }],
      };
      buildSideChain(res, name);
      out.push(res);
    }
    return out;
  }

  /* Rigid frame of a residue: origin CA, e1 along CA->C, e2 in the N-CA-C plane. */
  function residueFrame(N, CA, C) {
    const e1 = V.norm(V.sub(C, CA));
    const e2 = V.norm(V.perp(V.sub(N, CA), e1));
    return { R: V.m3cols(e1, e2, V.cross(e1, e2)), t: CA };
  }

  /* Screw axis of the rigid motion taking residue i's frame to residue i+1's: the axis of a regular
     helix or strand. Returns {u (unit, N->C direction of travel), point, rise (A/residue),
     turn (deg/residue)}. */
  function screwAxis(r0, r1) {
    const F0 = residueFrame(r0.N, r0.CA, r0.C), F1 = residueFrame(r1.N, r1.CA, r1.C);
    const R = V.m3mul(F1.R, V.m3T(F0.R));
    const t = V.sub(F1.t, V.m3vec(R, F0.t));
    const cosT = Math.max(-1, Math.min(1, (R[0] + R[4] + R[8] - 1) / 2));
    const theta = Math.acos(cosT);
    let u = [R[7] - R[5], R[2] - R[6], R[3] - R[1]];
    if (V.len(u) < 1e-9) u = V.norm(t); else u = V.norm(u);
    let rise = V.dot(t, u);
    const tp = V.sub(t, V.scale(u, rise));
    const cot = Math.abs(Math.sin(theta / 2)) > 1e-9 ? Math.cos(theta / 2) / Math.sin(theta / 2) : 0;
    const point = V.scale(V.add(tp, V.scale(V.cross(u, tp), cot)), 0.5);
    let turn = theta / D2R;
    if (rise < 0) { u = V.scale(u, -1); rise = -rise; turn = -turn; }
    return { u, point, rise, turn };
  }

  /* Axis of a regular stretch (helix or strand): the mean of the screw axes between consecutive
     residues. Exact for ideal segments, where PCA of the CA atoms tilts by a few degrees for
     non-integer numbers of turns. Returns {u, point} or null when the residues are not a
     consistent regular segment. */
  function segmentAxis(residues) {
    const bb = residues.map(r => ({ N: pos(r, 'N'), CA: pos(r, 'CA'), C: pos(r, 'C') }));
    if (bb.length < 3 || bb.some(b => !b.N || !b.CA || !b.C)) return null;
    let su = [0, 0, 0]; const pts = [];
    for (let i = 0; i + 1 < bb.length; i++) {
      if (V.dist(bb[i].C, bb[i + 1].N) > 2.0) return null;
      const ax = screwAxis(bb[i], bb[i + 1]);
      su = V.add(su, ax.u);
      pts.push(projectOnAxis(bb[i].CA, ax));
    }
    if (V.len(su) / (bb.length - 1) < 0.95) return null;
    return { u: V.norm(su), point: V.centroid(pts) };
  }

  /* Project p onto the line (point, u). */
  function projectOnAxis(p, ax) { return V.addScaled(ax.point, ax.u, V.dot(V.sub(p, ax.point), ax.u)); }

  /* Ideal regular segment of n residues with the given conformation, returned together with its
     screw axis (in the same canonical frame). */
  function idealSegment(n, conf) {
    const tors = [];
    for (let i = 0; i < n; i++) tors.push({ phi: conf.phi, psi: conf.psi, omega: 180 });
    const bb = backboneFromTorsions(tors.length > 1 ? tors : [tors[0], tors[0]]).slice(0, n);
    const probe = backboneFromTorsions([tors[0], tors[0], tors[0]]);
    const ax0 = screwAxis(probe[0], probe[1]);
    // The canonical frame is the same for both builds, so the axis carries over.
    return { bb, axis: ax0 };
  }

  function oneLetter(name) { return THREE_TO_ONE[name] || 'X'; }

  window.Geom = {
    BB, CONF, SC, AA_NAMES, THREE_TO_ONE, ONE_TO_THREE, ALIASES, BACKBONE,
    place, getAtom, pos, idealCB, buildSideChain, chiAngles, carbonylO,
    backboneFromTorsions, extendForward, extendBackward, residuesFromBackbone,
    residueFrame, screwAxis, segmentAxis, projectOnAxis, idealSegment, oneLetter,
  };
})();
