/* Side panels: tool options, sketch plane, residue inspector, chain list and chain tools, and the
   sequence strip. Panels are re-rendered from state with innerHTML; events are delegated. */
(function () {
  'use strict';
  const App = window.App, $ = App.$;
  const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const fx = (x, d) => (x === null || x === undefined || isNaN(x)) ? '' : (+x).toFixed(d === undefined ? 1 : d);

  const aaOptions = (cur, extra) => {
    const names = Geom.AA_NAMES.slice();
    if (cur && !names.includes(cur)) names.unshift(cur);
    return names.map(n => `<option value="${n}"${n === cur ? ' selected' : ''}>${n} (${Geom.oneLetter(n)})</option>`).join('') + (extra || '');
  };
  const confOptions = (cur, keys) => (keys || Object.keys(Geom.CONF)).map(k => `<option value="${k}"${k === cur ? ' selected' : ''}>${Geom.CONF[k].label}</option>`).join('');
  const sel = (opt, pairs) => `<select data-opt="${opt}">` + pairs.map(([v, l]) => `<option value="${v}"${String(App.opts[opt]) === String(v) ? ' selected' : ''}>${l}</option>`).join('') + '</select>';
  const num = (opt, step, w) => `<input type="number" data-opt="${opt}" value="${App.opts[opt]}" step="${step || 1}" style="width:${w || 64}px">`;

  /* ---------------- tool panel ---------------- */
  function loopBlock() {
    return `<h3>Loops (Connect / Draw)</h3><div class="kv">
      <label>Length</label><span>${num('loopLen', 1, 54)} <span class="muted small">0 = automatic</span></span>
      <label>Residue</label>${sel('loopRes', [['ALA', 'ALA'], ['GLY', 'GLY']])}
      <label>Renumber</label><label class="check" style="padding:0"><input type="checkbox" data-opt="renumber"${App.opts.renumber ? ' checked' : ''}> merged chains from their first number</label>
    </div>
    <p class="note">Loop residues stay in the αR, β, PPII and αL regions of the Ramachandran plot. αL (φ ≈ +60°) is allowed but penalised, because tight turns and 2-residue β-hairpins need it; in a design those positions usually become Gly, Asn or Asp.</p>`;
  }
  function refineBlock() {
    const o = App.opts;
    return `<h3>Refine (R)</h3><div class="kv">
      <label>Hold shape</label>${sel('refineHold', [['0.03', 'Loose'], ['0.15', 'Normal'], ['0.6', 'Tight']])}
      <label>Clash target</label>${sel('refineClash', [['0.80', 'Permissive'], ['0.88', 'Normal'], ['0.95', 'Strict']])}
      <label>Close gaps ≤</label><span>${num('refineGap', 0.5, 54)} Å</span>
      <label>Steps</label>${num('refineIters', 50, 54)}
      <label>Side chains</label>${sel('refineSide', [['ideal', 'Rebuild at same χ'], ['repack', 'Repack (search χ)'], ['keep', 'Leave alone']])}
    </div>
    <div class="btns"><button data-act="refine" ${App.sel.size ? '' : 'disabled'}>✦ Refine selection</button><button data-act="check" ${App.sel.size ? '' : 'disabled'}>Check geometry</button></div>
    <div class="result" id="refineResult">${App.refineCheckText || ''}</div>
    <p class="note">Pulls bond lengths, angles and peptide units back to ideal values, pushes apart
    atoms that overlap, and closes small chain breaks. Residues outside the selection stay fixed and
    act as context, so select a margin around whatever you are fixing. φ/ψ are kept close to what
    they are — this regularises a model, it does not refold it.
    <b>Repack</b> also searches each selected side chain's χ angles for a rotamer that fits, which
    clears clashes that moving the backbone alone cannot.</p>`;
  }

  function renderTool() {
    const t = App.tool, o = App.opts;
    let h = '';
    if (t === 'select') {
      h = `<h2>Select &amp; transform</h2><div class="kv">
        <label>Click picks</label><select id="selMode"><option value="residue"${App.selMode === 'residue' ? ' selected' : ''}>Residue</option><option value="element"${App.selMode === 'element' ? ' selected' : ''}>SS element</option><option value="chain"${App.selMode === 'chain' ? ' selected' : ''}>Chain</option></select>
        <label>Gizmo</label><span class="seg"><button data-gspace="world" class="${App.gizmo.space === 'world' ? 'active' : ''}">World</button><button data-gspace="local" class="${App.gizmo.space === 'local' ? 'active' : ''}">Local</button></span>
        <label>Ctrl snap</label><span><input type="number" id="snapMove" value="${App.gizmo.snapMove}" step="0.5" style="width:50px"> Å <input type="number" id="snapRot" value="${App.gizmo.snapRot}" step="5" style="width:50px"> °</span>
        <label>Plane</label><label class="check" style="padding:0"><input type="checkbox" id="showPlane"${App.showPlane ? ' checked' : ''}> always show sketch plane</label>
      </div>
      <p class="note">Local axes follow the selection's principal axes (for a helix, X′ is its long axis), or the residue frame for one residue.</p>` + refineBlock() + loopBlock();
    } else if (t === 'helix') {
      h = `<h2>Helix — extrude from plane</h2><div class="kv">
        <label>Type</label><select data-opt="helixConf">${confOptions(o.helixConf, ['alpha', 'helix310', 'pi'])}</select>
        <label>Click length</label><span>${num('helixLen', 1, 54)} residues</span>
        <label>On plane</label>${sel('helixEnd', [['start', 'N-terminus'], ['end', 'C-terminus']])}
        <label>Phase</label><span>${num('helixPhase', 10, 54)} °</span>
        <label>Residue</label><select data-opt="newRes">${aaOptions(o.newRes)}</select>
      </div>
      <p class="note">Click the plane for a helix of the click length, or press and drag along the dashed normal: up or down sets length and direction. ${(Builders.riseOf(Geom.CONF[o.helixConf])).toFixed(2)} Å rise per residue. The footprint circle (10 Å) is a typical helix-helix packing distance.</p>`;
    } else if (t === 'sheet') {
      h = `<h2>β-sheet — lay on plane</h2><div class="kv">
        <label>Strands</label>${num('sheetStrands', 1, 54)}
        <label>Pairing</label>${sel('sheetParallel', [['false', 'Antiparallel'], ['true', 'Parallel']])}
        <label>Spacing</label><span>${num('sheetSpacing', 0.1, 54)} Å</span>
        <label>Side</label>${sel('sheetSide', [['1', 'Left of stroke'], ['-1', 'Right of stroke']])}
        <label>Twist</label><span>${num('sheetTwist', 5, 54)} ° / strand</span>
        <label>Register</label><span>${num('sheetShift', 1, 54)} res / strand</span>
        <label>Click length</label><span>${num('sheetLen', 1, 54)} residues</span>
        <label>Residue</label><select data-opt="newRes">${aaOptions(o.newRes)}</select>
      </div>
      <p class="note">Drag across the plane to draw the first strand (N→C along the stroke). Strand axes are 4.8 Å apart in a typical sheet. Natural sheets twist left-handed about the inter-strand axis; try −15° to −20°.</p>`;
    } else if (t === 'draw') {
      h = `<h2>Draw loop</h2><div class="kv">
        <label>Draw on</label>${sel('drawSurface', [['auto', 'Auto'], ['plane', 'Sketch plane'], ['screen', 'Screen']])}
        <label>Spacing</label><span>${num('drawSpacing', 0.1, 54)} Å / residue</span>
      </div>
      <p class="note">The path is resampled every “spacing” Å along its length and a backbone is fitted to it (CA positions, phi/psi kept in allowed regions). Start on a C-terminus ▲ to grow from it, end on an N-terminus ● to join it. <b>Auto</b>: on the screen at the depth of the attached ends, otherwise on the sketch plane. 3.8 Å is fully extended; 3.0–3.3 Å follows a curve more loosely.</p>` + loopBlock();
    }
    $('toolPanel').innerHTML = h;
  }

  /* ---------------- plane panel ---------------- */
  function renderPlane() {
    const show = App.tool !== 'select' || App.showPlane;
    if (!show) { $('planePanel').innerHTML = ''; return; }
    const o = App.planeOrigin(), n = App.plane.n;
    $('planePanel').innerHTML = `<h2>Sketch plane</h2>
      <div class="btns">
        <button data-plane="xy" title="Normal +Z">XY</button><button data-plane="yz" title="Normal +X">YZ</button><button data-plane="zx" title="Normal +Y">ZX</button>
        <button data-plane="view" title="Facing the camera through the view centre">View</button>
        <button data-plane="fit" title="Best-fit plane through the selected residues (e.g. a sheet)">Fit sel</button>
        <button data-plane="axis" title="Perpendicular to the selection's long axis, at its C-terminal end (to stack or continue a helix)">⟂ Sel axis</button>
      </div>
      <div class="btns">
        <button data-plane="flip">Flip</button><button data-plane="toSel" title="Move the plane origin to the selection centre">Origin → sel</button>
        <button data-plane="zero">Origin → 0</button><button data-plane="look" title="Look straight at the plane">Look at</button>
      </div>
      <div class="kv">
        <label>Offset</label><span><input type="number" id="planeOffset" value="${fx(App.plane.offset, 1)}" step="1" style="width:64px"> Å along n</span>
        <label>Snap</label><span>${num('snap', 0.5, 54)} Å grid</span>
        <label>Origin</label><span class="muted small">${fx(o[0])}, ${fx(o[1])}, ${fx(o[2])}</span>
        <label>Normal</label><span class="muted small">${fx(n[0], 2)}, ${fx(n[1], 2)}, ${fx(n[2], 2)}</span>
      </div>`;
  }

  App.setPlane = function (kind) {
    const P = App.plane, cam = App.renderer.cam;
    const pts = App.selectionPoints();
    const setFrame = (n, u) => { P.n = V.norm(n); P.u = V.norm(V.perp(u, P.n)); if (V.len(P.u) < 1e-6) P.u = V.anyPerp(P.n); P.v = V.cross(P.n, P.u); };
    const keepOrigin = () => { const o = App.planeOrigin(); P.base = o; P.offset = 0; };
    switch (kind) {
      case 'xy': keepOrigin(); setFrame([0, 0, 1], [1, 0, 0]); break;
      case 'yz': keepOrigin(); setFrame([1, 0, 0], [0, 1, 0]); break;
      case 'zx': keepOrigin(); setFrame([0, 1, 0], [0, 0, 1]); break;
      case 'view': setFrame(cam.back(), cam.right()); P.base = cam.target.slice(); P.offset = 0; break;
      case 'fit': {
        if (pts.length < 3) { App.status('Select at least 3 residues to fit a plane.', 'err'); return; }
        const p = V.pca(pts);
        let n = p.axes[2]; if (V.dot(n, cam.back()) < 0) n = V.scale(n, -1);
        setFrame(n, p.axes[0]); P.base = p.center; P.offset = 0; break;
      }
      case 'axis': {
        if (pts.length < 2) { App.status('Select at least 2 residues (e.g. a helix) to use its axis.', 'err'); return; }
        const p = V.pca(pts);
        const seg = Geom.segmentAxis(App.selectedResidues().filter(r => !r.het));
        let n = seg ? seg.u : p.axes[0];
        const c = seg ? seg.point : p.center;
        if (V.dot(n, V.sub(pts[pts.length - 1], pts[0])) < 0) n = V.scale(n, -1);
        setFrame(n, p.axes[1]);
        P.base = V.addScaled(c, n, V.dot(V.sub(pts[pts.length - 1], c), n));
        P.offset = 0; break;
      }
      case 'flip': P.n = V.scale(P.n, -1); P.v = V.scale(P.v, -1); P.offset = -P.offset; break;
      case 'toSel': if (!pts.length) { App.status('Nothing selected.', 'err'); return; } P.base = V.centroid(pts); P.offset = 0; break;
      case 'zero': P.base = [0, 0, 0]; P.offset = 0; break;
      case 'look': cam.lookAlong(P.n, P.v); cam.target = App.planeOrigin(); App.updateSceneRadius(); App.cameraChanged(); return;
    }
    App.updateSceneRadius();
    App.overlayDirty = true;
    renderPlane();
  };

  /* ---------------- inspector ---------------- */
  function renderInspector() {
    const rs = App.selectedResidues();
    const el = $('inspector');
    if (!rs.length) {
      el.innerHTML = `<h2>Selection</h2><p class="note">Nothing selected.</p>
        <h3>Getting started</h3>
        <ol class="note" style="padding-left:16px;margin:0">
          <li><b>Helix</b> (H): click the grid to extrude an α-helix; drag along the normal to set its length.</li>
          <li><b>Sheet</b> (B): drag across the grid to lay a β-sheet.</li>
          <li><b>Select</b> (S): click residues, then move or rotate them with the gizmo.</li>
          <li>Select a C-terminus, then an N-terminus (Ctrl-click), and press <b>Connect</b> (C) to build a loop. Whole chains selected in order are joined in that order.</li>
          <li><b>Draw loop</b> (D): sketch a path from a ▲ end to a ● end.</li>
          <li><b>Export</b> a PDB or mmCIF as the template for ProteinMPNN or other inverse-folding tools.</li>
        </ol>`;
      return;
    }
    let h = '';
    if (rs.length === 1) h += singleBlock(rs[0]);
    else {
      const chains = App.selectedChains();
      h += `<h2>Selection</h2><p class="note">${rs.length} residues${chains.length ? ' in chain' + (chains.length > 1 ? 's ' : ' ') + chains.map(c => c.id).join(', ') : ''}${rs.some(r => r.het) ? ' (incl. ligands)' : ''}.</p>`;
      h += `<h3>Residue type</h3><div class="row"><select id="mMut">${aaOptions('ALA')}</select><button data-act="mutateAll">Mutate all</button></div>`;
    }
    const anyChain = rs.some(r => !r.het);
    if (anyChain) {
      h += `<h3>Convert</h3><div class="btns">
        <button data-act="conv" data-mode="backbone" title="Keep N, CA, C, O; residue names unchanged">Backbone only</button>
        <button data-act="conv" data-mode="backboneG" title="Keep N, CA, C, O and rename to GLY">Backbone → GLY</button>
        <button data-act="conv" data-mode="ala" title="N, CA, C, O + ideal CB">Poly-ALA</button>
        <button data-act="conv" data-mode="gly">Poly-GLY</button></div>`;
      if (rs.length > 1) {
        h += `<h3>Secondary structure label</h3><div class="btns"><button data-act="ss" data-ss="H">Helix</button><button data-act="ss" data-ss="E">Strand</button><button data-act="ss" data-ss="L">Loop</button><button data-act="autoSS" title="Assign from phi/psi">From φ/ψ</button></div>
        <h3>Conformation</h3><div class="row"><select data-opt="applyConf">${confOptions(App.opts.applyConf)}</select><button data-act="applyConf" title="Set phi/psi of every selected residue, rotating the C-terminal side">Apply φ/ψ</button></div>`;
        const one = App.selectedChains().length === 1;
        h += `<h3>Numbering</h3><div class="row">${one ? `<label>From</label><input type="number" id="renumStart" value="${rs[0].num}" style="width:56px"><button data-act="renumSel">Renumber</button>` : ''}
          <label>Shift</label><input type="number" id="shiftBy" value="0" style="width:50px"><button data-act="shiftNum">Apply</button></div>`;
      }
    }
    h += `<h3>Transform</h3>
      <div class="row tight"><label>Move</label><span class="xyz"><input type="number" id="tx" value="0" step="0.5"><input type="number" id="ty" value="0" step="0.5"><input type="number" id="tz" value="0" step="0.5"></span><button data-act="translate">Apply</button></div>
      <div class="row tight"><label>Rotate</label><input type="number" id="rAng" value="90" step="15" style="width:50px">° about
        <select id="rAxis"><option value="x">X</option><option value="y">Y</option><option value="z">Z</option><option value="major">long axis</option><option value="minor">short axis</option><option value="normal">plane normal</option></select>
        <button data-act="rotate">Apply</button></div>
      <p class="note">Rotation is about the selection centre. Or drag the gizmo in the view.</p>`;
    h += `<h3>Structure</h3><div class="btns">
      <button data-act="refine" title="Regularise bond lengths, angles and clashes (R)">✦ Refine</button>
      <button data-act="connect" title="Build a loop between two ends (C)">⌒ Connect</button>
      <button data-act="duplicate" title="Copy into a new chain (Ctrl+D)">Duplicate</button>
      ${App.selectedChains().length >= 2 ? '<button data-act="join" title="Concatenate the selected chains in selection order without building a loop">Join chains</button>' : ''}
      <button data-act="delete">Delete</button>
      <button data-act="fitSel">Zoom to</button></div>
      <div class="result" id="inspResult"></div>`;
    el.innerHTML = h;
  }

  function singleBlock(r) {
    const info = App.index.get(r.uid);
    if (!info.chain) {
      const c = V.centroid(r.atoms.map(a => a.p));
      return `<h2>Ligand</h2><div class="kv"><label>Name</label><span>${esc(r.name)} ${r.num}${esc(r.icode || '')}</span>
        <label>Chain</label><span>${esc(r.chainId || '')}</span><label>Atoms</label><span>${r.atoms.length}</span>
        <label>Centre</label><span class="muted small">${fx(c[0])}, ${fx(c[1])}, ${fx(c[2])}</span></div>`;
    }
    const ch = info.chain, i = info.idx;
    const ca = Geom.pos(r, 'CA') || [0, 0, 0];
    const f = Model.phi(ch, i), p = Model.psi(ch, i), w = Model.omega(ch, i);
    const term = i === 0 ? ' · N-terminus' : i === ch.residues.length - 1 ? ' · C-terminus' : '';
    return `<h2>Residue</h2><div class="kv">
      <label>Chain</label><span>${esc(ch.id)} <span class="muted small">(${i + 1} of ${ch.residues.length}${term})</span></span>
      <label>Type</label><select id="iType">${aaOptions(r.name)}</select>
      <label>Number</label><span><input type="number" id="iNum" value="${r.num}" style="width:64px"> <input type="text" id="iIcode" value="${esc(r.icode || '')}" maxlength="1" style="width:26px" title="Insertion code"></span>
      <label>SS</label><select id="iSS">${[['H', 'Helix'], ['E', 'Strand'], ['L', 'Loop']].map(([k, l]) => `<option value="${k}"${r.ss === k ? ' selected' : ''}>${l}</option>`).join('')}</select>
      <label>CA</label><span class="xyz"><input type="number" id="cax" value="${fx(ca[0], 2)}" step="0.5"><input type="number" id="cay" value="${fx(ca[1], 2)}" step="0.5"><input type="number" id="caz" value="${fx(ca[2], 2)}" step="0.5"></span>
    </div>
    <h3>Backbone torsions</h3><div class="kv">
      <label>φ</label><span><input type="number" data-tor="phi" value="${fx(f)}" step="5" ${f === null ? 'disabled' : ''}> °</span>
      <label>ψ</label><span><input type="number" data-tor="psi" value="${fx(p)}" step="5" ${p === null ? 'disabled' : ''}> °</span>
      <label>ω</label><span><input type="number" data-tor="omega" value="${fx(w)}" step="5" ${w === null ? 'disabled' : ''}> ° <span class="muted small">to next</span></span>
      <label>Moves</label>${sel('torsionSide', [['C', 'C-terminal side'], ['N', 'N-terminal side']])}
    </div>
    <div class="btns">${i < ch.residues.length - 1 ? '<button data-act="split" title="Make the residues after this one a new chain">Split chain after</button>' : ''}</div>`;
  }

  /* ---------------- chains ---------------- */
  function renderChains() {
    const el = $('chainList');
    const m = App.model;
    if (!m.chains.length && !m.ligands.length) {
      el.innerHTML = '<div class="empty">No chains yet. Use the Helix or Sheet tool, or Open / Import a PDB or mmCIF file.</div>';
    } else {
      el.innerHTML = m.chains.map(c => {
        const col = Render.chainColor(c.id).map(x => Math.round(x * 255)).join(',');
        const n = c.residues.length, nsel = c.residues.filter(r => App.sel.has(r.uid)).length;
        const cls = nsel === n ? 'sel' : nsel ? 'part' : '';
        const rng = n ? c.residues[0].num + '–' + c.residues[n - 1].num : '';
        const brk = App.derived.breaks.filter(b => b.chain === c).length;
        return `<div class="chain ${cls}" data-chain="${esc(c.id)}"><span class="dot" style="background:rgb(${col})"></span><span class="cid">${esc(c.id)}</span>
          <span class="info">${n} res · ${rng}${brk ? ` · <span style="color:var(--bad)">${brk} break${brk > 1 ? 's' : ''}</span>` : ''}</span>
          <button class="eye ${App.hidden.has(c.id) ? 'off' : ''}" data-eye="${esc(c.id)}" title="Show / hide">👁</button></div>`;
      }).join('') + (m.ligands.length ? `<div class="chain" data-chain="__lig"><span class="dot" style="background:rgb(140,217,115)"></span><span class="cid">⬡</span><span class="info">${m.ligands.length} ligand group(s)</span>
        <button class="eye ${App.hidden.has('lig') ? 'off' : ''}" data-eye="lig" title="Show / hide">👁</button></div>` : '');
    }
    // Chain tools
    const target = App.selectedChains();
    const ops = $('chainOps');
    if (!m.chains.length) { ops.innerHTML = ''; return; }
    const names = target.length ? target.map(c => c.id).join(', ') : 'all chains';
    const single = target.length === 1 ? target[0] : (m.chains.length === 1 ? m.chains[0] : null);
    let h = `<h2>Chain tools</h2><p class="note" style="margin-top:0">Apply to: <b>${esc(names)}</b>${target.length ? '' : ' (select residues to narrow this)'}</p>
      <div class="btns"><button data-cop="backbone" title="Keep N, CA, C, O">Backbone</button><button data-cop="backboneG" title="Backbone only, all residues named GLY">Backbone→GLY</button>
      <button data-cop="ala" title="N, CA, C, O + ideal CB, named ALA">Poly-ALA</button><button data-cop="gly">Poly-GLY</button></div>
      <div class="row"><label>Renumber from</label><input type="number" id="cRenum" value="1" style="width:56px"><button data-cop="renum">Apply</button></div>
      <div class="btns"><button data-cop="refine" title="Regularise bond lengths, angles and clashes over the whole chain">✦ Refine</button><button data-cop="autoss" title="Re-assign helix/strand/loop labels from phi/psi">SS from φ/ψ</button><button data-cop="select">Select</button><button data-cop="delete">Delete</button></div>`;
    if (single) {
      h += `<h3>Chain ${esc(single.id)}</h3>
        <div class="row"><label>ID</label><input type="text" id="cId" value="${esc(single.id)}" maxlength="4" style="width:44px"><button data-cop="rename">Rename</button></div>
        <div class="row tight"><label>Extend</label><select id="cEnd"><option value="C">C-term</option><option value="N">N-term</option></select> by
          ${num('extendN', 1, 44)} as <select data-opt="extendConf">${confOptions(App.opts.extendConf)}</select><button data-cop="extend">Add</button></div>`;
    }
    ops.innerHTML = h;
  }

  /* ---------------- sequence ---------------- */
  let seqVersion = -1;
  function renderSeq() {
    const el = $('seqInner');
    if (seqVersion !== App.modelVersion) {
      seqVersion = App.modelVersion;
      const m = App.model;
      if (!m.chains.length) { el.innerHTML = '<div class="empty">Sequence appears here. Click or drag across residues to select them.</div>'; }
      else {
        el.innerHTML = m.chains.map(c => {
          const rs = c.residues; let groups = '';
          for (let i = 0; i < rs.length; i += 10) {
            let aa = '';
            for (let k = i; k < Math.min(i + 10, rs.length); k++) {
              const r = rs[k];
              const brk = k > 0 && !Model.bonded(rs[k - 1], r) ? ' brk' : '';
              aa += `<span class="r ${r.ss}${brk}" data-uid="${r.uid}" title="${esc(c.id)} ${esc(r.name)} ${r.num}${esc(r.icode || '')}">${Geom.oneLetter(r.name)}</span>`;
            }
            groups += `<div class="seqgrp"><span class="num">${rs[i].num}</span><span class="aa">${aa}</span></div>`;
          }
          return `<div class="seqrow"><span class="sid">${esc(c.id)}</span><div class="sres">${groups}</div></div>`;
        }).join('');
      }
    }
    el.querySelectorAll('span.r').forEach(s => s.classList.toggle('sel', App.sel.has(+s.dataset.uid)));
    App.highlightSeq();
  }
  let lastHov = null;
  App.highlightSeq = function () {
    if (lastHov) lastHov.classList.remove('hov');
    lastHov = App.hover !== null ? $('seqInner').querySelector(`span.r[data-uid="${App.hover}"]`) : null;
    if (lastHov) lastHov.classList.add('hov');
  };

  function renderCounts() {
    const m = App.model;
    const nres = m.chains.reduce((s, c) => s + c.residues.length, 0);
    const nb = App.derived.breaks.length, nc = App.derived.clashes.length;
    $('counts').innerHTML = `${m.chains.length} chain(s) · ${nres} residues` + (m.ligands.length ? ` · ${m.ligands.length} ligand(s)` : '') +
      (App.sel.size ? ` · <span style="color:var(--sel)">${App.sel.size} selected</span>` : '') +
      (nb ? ` · <span class="bad">${nb} chain break(s)</span>` : '') + (nc ? ` · <span class="warn">${nc} clash(es)</span>` : '');
  }

  App.refreshPanels = function () {
    renderTool(); renderPlane(); renderInspector(); renderChains(); renderSeq(); renderCounts();
    $('bUndo').disabled = !App.undoStack.length; $('bRedo').disabled = !App.redoStack.length;
    $('bConnect').disabled = App.sel.size < 2;
    $('bRefine').disabled = !App.sel.size;
    $('bDelete').disabled = !App.sel.size;
  };

  /* ---------------- actions ---------------- */
  function chainResidues() { return App.selectedResidues().filter(r => !r.het); }
  function opChains() { const t = App.selectedChains(); return t.length ? t : App.model.chains.slice(); }

  function convert(rs, mode) {
    const m = mode === 'backboneG' ? 'backbone' : mode;
    let n = 0;
    App.edit('Convert', () => { n = Model.convertResidues(rs, m, mode === 'backboneG'); });
    App.status({ backbone: 'Stripped ' + n + ' residue(s) to backbone N, CA, C, O.', backboneG: 'Stripped ' + n + ' residue(s) to backbone and renamed them GLY.', ala: 'Converted ' + n + ' residue(s) to alanine.', gly: 'Converted ' + n + ' residue(s) to glycine.' }[mode], 'ok');
  }

  function rotationAxis(kind, pts) {
    if (kind === 'x') return [1, 0, 0];
    if (kind === 'y') return [0, 1, 0];
    if (kind === 'z') return [0, 0, 1];
    if (kind === 'normal') return App.plane.n;
    if (pts.length < 2) return [0, 0, 1];
    const p = V.pca(pts);
    const seg = kind === 'major' && Geom.segmentAxis(App.selectedResidues().filter(r => !r.het));
    return seg ? seg.u : kind === 'major' ? p.axes[0] : p.axes[2];
  }

  function onInspectorClick(e) {
    const b = e.target.closest('button'); if (!b) return;
    const act = b.dataset.act; if (!act) return;
    const rs = App.selectedResidues();
    switch (act) {
      case 'mutateAll': {
        const name = $('mMut').value; let n = 0;
        App.edit('Mutate', () => { for (const r of rs) if (!r.het && Model.mutate(r, name)) n++; });
        App.status('Mutated ' + n + ' residue(s) to ' + name + '.', 'ok');
        break;
      }
      case 'conv': convert(chainResidues(), b.dataset.mode); break;
      case 'ss': App.edit('Set SS', () => { for (const r of chainResidues()) r.ss = b.dataset.ss; }); break;
      case 'autoSS': {
        // Assign from phi/psi, but only change the selected residues.
        App.edit('SS from φ/ψ', () => {
          for (const c of App.selectedChains()) {
            const keep = c.residues.map(r => r.ss);
            Model.autoSS(c);
            c.residues.forEach((r, k) => { if (!App.sel.has(r.uid)) r.ss = keep[k]; });
          }
        });
        break;
      }
      case 'applyConf': {
        const conf = Geom.CONF[App.opts.applyConf];
        let n = 0;
        App.edit('Apply φ/ψ', () => {
          for (const c of App.selectedChains()) {
            c.residues.forEach((r, i) => {
              if (!App.sel.has(r.uid)) return;
              if (Model.setTorsion(c, i, 'phi', conf.phi, 'C')) n++;
              Model.setTorsion(c, i, 'psi', conf.psi, 'C');
              r.ss = conf.ss;
            });
          }
        });
        App.status('Set ' + conf.label + ' φ/ψ on the selection (' + n + ' residues had a φ to set; the C-terminal side moved).', 'ok');
        break;
      }
      case 'renumSel': {
        const s = parseInt($('renumStart').value);
        if (isNaN(s)) return;
        App.edit('Renumber', () => Model.renumber(chainResidues(), s));
        warnDuplicates();
        break;
      }
      case 'shiftNum': {
        const k = parseInt($('shiftBy').value);
        if (!k) return;
        App.edit('Shift numbering', () => { for (const r of rs) r.num += k; });
        warnDuplicates();
        break;
      }
      case 'translate': {
        const d = [+$('tx').value || 0, +$('ty').value || 0, +$('tz').value || 0];
        App.edit('Move', () => { for (const r of rs) Model.moveAtoms(r.atoms, p => V.add(p, d)); });
        break;
      }
      case 'rotate': {
        const pts = App.selectionPoints();
        const c = V.centroid(pts), ax = rotationAxis($('rAxis').value, pts), ang = +$('rAng').value || 0;
        const R = V.rotAxis(ax, ang);
        App.edit('Rotate', () => { for (const r of rs) Model.moveAtoms(r.atoms, p => V.add(c, V.m3vec(R, V.sub(p, c)))); });
        break;
      }
      case 'refine': App.refine(); break;
      case 'check': checkGeometry(); break;
      case 'connect': App.connect(); break;
      case 'duplicate': App.duplicateSelection(); break;
      case 'delete': App.deleteSelection(); break;
      case 'fitSel': App.fitView(true); break;
      case 'join': {
        const cs = App.selectedChains(true);
        App.edit('Join chains', () => { for (let i = 1; i < cs.length; i++) Model.mergeChains(App.model, cs[0], cs[i]); });
        App.status('Joined chains into ' + cs[0].id + ' without a loop (the gaps show as chain breaks).');
        break;
      }
      case 'split': {
        const info = App.index.get(rs[0].uid);
        let nc;
        App.edit('Split chain', () => { nc = Model.splitChain(App.model, info.chain, info.idx); return nc ? undefined : false; });
        if (nc) App.status('Residues after ' + rs[0].num + ' are now chain ' + nc.id + '.', 'ok');
        break;
      }
    }
  }

  /* Report what is wrong with the selected residues without changing anything. */
  function checkGeometry() {
    const rs = App.selectedResidues().filter(r => !r.het);
    if (!rs.length) return;
    const a = Refine.analyse(App.model, rs, App.refineOpts());
    if (!a) { App.refineCheckText = ''; return; }
    const bits = [
      `${a.residues} residue(s)`,
      `bonds ${a.bondRms.toFixed(3)} Å RMS (worst ${a.bondMax.toFixed(2)})`,
      `angles ${a.angleRms.toFixed(3)} Å (worst ${a.angleMax.toFixed(2)})`,
      `${a.clashes} clash(es)` + (a.clashes ? `, worst overlap ${a.worstOverlap.toFixed(2)} Å` : ''),
    ];
    if (a.wrongChirality) bits.push(`<b style="color:var(--bad)">${a.wrongChirality} D-amino acid centre(s)</b>`);
    if (a.rama.length) bits.push(`${a.rama.length} residue(s) with φ outside the allowed range (${a.rama.slice(0, 4).map(x => x.res.name + x.res.num).join(', ')})`);
    for (const b of a.breaks) bits.push(`break after ${b.res.num}: ${b.gap.toFixed(2)} Å` + (b.closing ? ' (Refine will close it)' : ' (too wide — use Connect)'));
    const clean = a.bondRms < 0.02 && a.angleRms < 0.05 && !a.clashes && !a.breaks.length && !a.wrongChirality && !a.rama.length;
    App.refineCheckText = (clean ? '<span style="color:var(--ok)">Geometry is clean.</span> ' : '') + bits.join(' · ');
    const el = $('refineResult');
    if (el) el.innerHTML = App.refineCheckText;
  }

  function warnDuplicates() {
    for (const c of App.model.chains) {
      const seen = new Set();
      for (const r of c.residues) {
        const k = r.num + (r.icode || '');
        if (seen.has(k)) { App.status('Chain ' + c.id + ' now has duplicate residue number ' + k + '. Most tools need unique numbers per chain.', 'err'); return; }
        seen.add(k);
      }
    }
  }

  function onInspectorChange(e) {
    const t = e.target;
    const rs = App.selectedResidues();
    if (t.dataset.opt) return; // handled by the generic option handler
    if (rs.length !== 1) return;
    const r = rs[0], info = App.index.get(r.uid);
    if (t.id === 'iType') {
      const name = t.value;
      App.edit('Mutate', () => Model.mutate(r, name) ? undefined : false);
      App.status(name === r.name ? 'Residue ' + r.num + ' is now ' + name + '.' : 'Could not rebuild the side chain (backbone incomplete).', name === r.name ? 'ok' : 'err');
    } else if (t.id === 'iNum' || t.id === 'iIcode') {
      const n = parseInt($('iNum').value);
      if (isNaN(n)) return;
      App.edit('Renumber residue', () => { r.num = n; r.icode = $('iIcode').value.trim(); });
      App.status('Residue is now number ' + n + (r.icode || '') + '.', 'ok');
      warnDuplicates();
    } else if (t.id === 'iSS') {
      App.edit('Set SS', () => { r.ss = t.value; });
    } else if (t.id === 'cax' || t.id === 'cay' || t.id === 'caz') {
      const ca = Geom.pos(r, 'CA'); if (!ca) return;
      const target = [+$('cax').value, +$('cay').value, +$('caz').value];
      if (target.some(isNaN)) return;
      const d = V.sub(target, ca);
      App.edit('Move residue', () => Model.moveAtoms(r.atoms, p => V.add(p, d)));
    } else if (t.dataset.tor) {
      const v = parseFloat(t.value); if (isNaN(v)) return;
      const ok = App.edit('Set ' + t.dataset.tor, () => Model.setTorsion(info.chain, info.idx, t.dataset.tor, v, App.opts.torsionSide) ? undefined : false);
      if (!ok) App.status('That torsion is not defined here (chain end or break).', 'err');
    }
  }

  function onChainsClick(e) {
    const eye = e.target.closest('[data-eye]');
    if (eye) {
      const id = eye.dataset.eye;
      if (App.hidden.has(id)) App.hidden.delete(id); else App.hidden.add(id);
      App.projCache = null; App.sceneDirty = true; App.overlayDirty = true; renderChains();
      return;
    }
    const row = e.target.closest('.chain'); if (!row) return;
    const id = row.dataset.chain;
    const uids = id === '__lig' ? App.model.ligands.map(r => r.uid) : (App.model.chains.find(c => c.id === id) || { residues: [] }).residues.map(r => r.uid);
    const add = e.ctrlKey || e.metaKey || e.shiftKey;
    if (add) {
      const all = uids.every(u => App.sel.has(u));
      App.setSelection(uids, all ? 'remove' : 'add');
    } else App.setSelection(uids, 'set');
  }

  function onChainOps(e) {
    const b = e.target.closest('button'); if (!b || !b.dataset.cop) return;
    const chains = opChains();
    const all = chains.flatMap(c => c.residues);
    switch (b.dataset.cop) {
      case 'backbone': case 'backboneG': case 'ala': case 'gly': convert(all, b.dataset.cop); break;
      case 'renum': {
        const s = parseInt($('cRenum').value); if (isNaN(s)) return;
        App.edit('Renumber chains', () => { for (const c of chains) Model.renumber(c.residues, s); });
        App.status('Renumbered ' + chains.map(c => c.id).join(', ') + ' from ' + s + '.', 'ok');
        break;
      }
      case 'refine': App.refine(chains); break;
      case 'autoss': App.edit('SS from φ/ψ', () => { for (const c of chains) Model.autoSS(c); }); break;
      case 'select': App.setSelection(all.map(r => r.uid), 'set'); break;
      case 'delete': {
        const ids = chains.map(c => c.id);
        if (!confirm('Delete chain(s) ' + ids.join(', ') + '?')) return;
        App.edit('Delete chains', () => { App.model.chains = App.model.chains.filter(c => !chains.includes(c)); });
        App.status('Deleted chain(s) ' + ids.join(', ') + '.');
        break;
      }
      case 'rename': {
        const single = App.selectedChains().length === 1 ? App.selectedChains()[0] : App.model.chains[0];
        const id = $('cId').value.trim();
        if (!id || /\s/.test(id)) { App.status('Chain IDs cannot be empty or contain spaces.', 'err'); return; }
        if (App.model.chains.some(c => c !== single && c.id === id)) { App.status('Chain ' + id + ' already exists.', 'err'); return; }
        const old = single.id;
        App.edit('Rename chain', () => { single.id = id; });
        if (App.hidden.has(old)) { App.hidden.delete(old); App.hidden.add(id); }
        App.status('Chain ' + old + ' is now ' + id + (id.length > 1 ? ' (PDB export keeps only the first character; mmCIF keeps it all)' : '') + '.', id.length > 1 ? 'err' : 'ok');
        break;
      }
      case 'extend': {
        const single = App.selectedChains().length === 1 ? App.selectedChains()[0] : App.model.chains[0];
        const end = $('cEnd').value, n = Math.max(1, App.opts.extendN | 0);
        const conf = Geom.CONF[App.opts.extendConf];
        let added = 0;
        App.edit('Extend chain', () => { added = Model.extendChain(App.model, single, end, n, conf, App.opts.newRes); return added ? undefined : false; });
        if (added) {
          const rs = end === 'C' ? single.residues.slice(-added) : single.residues.slice(0, added);
          App.setSelection(rs.map(r => r.uid), 'set');
          App.status('Added ' + added + ' ' + App.opts.newRes + ' (' + conf.label + ') at the ' + end + '-terminus of chain ' + single.id + '.', 'ok');
        } else App.status('Could not extend: the terminal residue needs N, CA and C.', 'err');
        break;
      }
    }
  }

  /* Generic option inputs (data-opt) anywhere in the side panels. */
  function onOption(e) {
    const t = e.target, k = t.dataset.opt;
    if (!k) return;
    let v = t.type === 'checkbox' ? t.checked : t.value;
    if (typeof App.opts[k] === 'number') { v = parseFloat(v); if (isNaN(v)) return; }
    if (typeof App.opts[k] === 'boolean' && typeof v === 'string') v = v === 'true';
    App.opts[k] = v;
    if (k === 'helixConf' || k === 'helixLen' || k === 'sheetStrands') { App.setTool(App.tool); }
    if (k === 'snap' || k === 'drawSurface') App.overlayDirty = true;
  }

  App.initPanels = function () {
    const right = $('right'), left = $('left');
    right.addEventListener('change', onOption);
    left.addEventListener('change', onOption);
    $('inspector').addEventListener('click', onInspectorClick);
    $('inspector').addEventListener('change', onInspectorChange);
    $('chainList').addEventListener('click', onChainsClick);
    $('chainOps').addEventListener('click', onChainOps);
    $('toolPanel').addEventListener('change', e => {
      if (e.target.id === 'selMode') App.selMode = e.target.value;
      if (e.target.id === 'snapMove') App.gizmo.snapMove = Math.max(0.01, +e.target.value || 1);
      if (e.target.id === 'snapRot') App.gizmo.snapRot = Math.max(0.1, +e.target.value || 15);
      if (e.target.id === 'showPlane') { App.showPlane = e.target.checked; renderPlane(); App.overlayDirty = true; }
    });
    $('toolPanel').addEventListener('click', e => {
      const b = e.target.closest('[data-gspace]'); if (!b) return;
      App.gizmo.space = b.dataset.gspace; App.overlayDirty = true; renderTool();
    });
    // The tool panel carries Refine and Check geometry, which are the same actions as in the
    // inspector.
    $('toolPanel').addEventListener('click', onInspectorClick);
    $('planePanel').addEventListener('click', e => { const b = e.target.closest('[data-plane]'); if (b) App.setPlane(b.dataset.plane); });
    $('planePanel').addEventListener('change', e => {
      if (e.target.id === 'planeOffset') { App.plane.offset = parseFloat(e.target.value) || 0; App.updateSceneRadius(); App.overlayDirty = true; renderPlane(); }
    });
    // Sequence strip: click / drag to select ranges.
    let seqDrag = null;
    $('seqInner').addEventListener('mousedown', e => {
      const s = e.target.closest('span.r'); if (!s) return;
      const uid = +s.dataset.uid;
      const mode = e.ctrlKey || e.metaKey ? 'toggle' : e.shiftKey ? 'extend' : 'set';
      if (mode === 'extend' && App.lastClick !== undefined) {
        const a = App.index.get(App.lastClick), b = App.index.get(uid);
        if (a && b && a.chain === b.chain) {
          const lo = Math.min(a.idx, b.idx), hi = Math.max(a.idx, b.idx);
          App.setSelection(a.chain.residues.slice(lo, hi + 1).map(r => r.uid), 'add');
          return;
        }
      }
      seqDrag = { start: uid, base: mode === 'set' ? new Set() : new Set(App.sel), remove: mode === 'toggle' && App.sel.has(uid) };
      App.lastClick = uid;
      applySeqDrag(uid);
      e.preventDefault();
    });
    $('seqInner').addEventListener('mouseover', e => {
      const s = e.target.closest('span.r'); if (!s) return;
      const uid = +s.dataset.uid;
      if (seqDrag) applySeqDrag(uid);
      if (App.hover !== uid) { App.hover = uid; App.sceneDirty = true; App.refreshHud(); App.highlightSeq(); }
    });
    $('seqInner').addEventListener('mouseleave', () => { if (App.hover !== null) { App.hover = null; App.sceneDirty = true; App.refreshHud(); App.highlightSeq(); } });
    window.addEventListener('mouseup', () => { seqDrag = null; });
    function applySeqDrag(uid) {
      const a = App.index.get(seqDrag.start), b = App.index.get(uid);
      let range = [seqDrag.start];
      if (a && b && a.chain === b.chain) {
        const lo = Math.min(a.idx, b.idx), hi = Math.max(a.idx, b.idx);
        range = a.chain.residues.slice(lo, hi + 1).map(r => r.uid);
      }
      const s = new Set(seqDrag.base);
      for (const u of range) { if (seqDrag.remove) s.delete(u); else s.add(u); }
      App.sel.clear(); App.selOrder = App.selOrder.filter(u => s.has(u));
      for (const u of s) { App.sel.add(u); if (!App.selOrder.includes(u)) App.selOrder.push(u); }
      App.selectionChanged();
    }
  };
})();
