/* PDB / mmCIF import and export.

   Import keeps the first model and the first alternate location, drops hydrogens and waters,
   folds common residue-name variants (MSE, HID, CYX, ...) into the standard 20, and puts
   amino-acid residues (anything with N, CA and C) into chains and everything else into ligands. */
(function () {
  'use strict';

  const WATER = new Set(['HOH', 'WAT', 'DOD', 'H2O', 'SOL', 'TIP', 'TIP3']);

  function guessElement(atomName, resName) {
    const n = atomName.replace(/[^A-Za-z]/g, '').toUpperCase();
    if (!n) return 'C';
    if (Geom.THREE_TO_ONE[resName] || Geom.ALIASES[resName]) return n[0];
    const two = n.slice(0, 2);
    if (['CL', 'BR', 'FE', 'ZN', 'MG', 'MN', 'NA', 'CA', 'CU', 'CO', 'NI', 'SE', 'CD', 'HG'].includes(two) && atomName.length >= 2 && atomName[0] !== ' ') return two[0] + two[1].toLowerCase();
    return n[0];
  }

  /* Collect flat atom records into the model. rec: {chain, resSeq, icode, resName, name, el, alt, p, het} */
  function assemble(records, ssRanges) {
    const model = Model.create();
    const chains = new Map();
    let cur = null, curKey = null;
    const altOf = new Map();
    for (const a of records) {
      if (a.el === 'H' || a.el === 'D') continue;
      if (WATER.has(a.resName)) continue;
      const key = a.chain + '|' + a.resSeq + '|' + a.icode + '|' + a.resName;
      if (a.alt) {
        const k2 = key + '|' + a.name;
        if (altOf.has(key) && altOf.get(key) !== a.alt) continue;
        altOf.set(key, a.alt);
        if (altOf.has(k2)) continue;
        altOf.set(k2, true);
      }
      if (key !== curKey) {
        cur = { chain: a.chain, name: a.resName, num: a.resSeq, icode: a.icode, ss: 'L', atoms: [], het: a.het };
        curKey = key;
        if (!chains.has(a.chain)) chains.set(a.chain, []);
        chains.get(a.chain).push(cur);
      }
      if (cur.atoms.some(x => x.name === a.name)) continue;
      cur.atoms.push({ name: a.name, el: a.el, p: a.p });
    }
    for (const [cid, list] of chains) {
      const residues = [];
      for (const r of list) {
        let name = r.name;
        if (Geom.ALIASES[name]) {
          if (name === 'MSE') for (const at of r.atoms) if (at.name === 'SE') { at.name = 'SD'; at.el = 'S'; }
          name = Geom.ALIASES[name];
        }
        const isAA = ['N', 'CA', 'C'].every(n => r.atoms.some(x => x.name === n));
        if (isAA) {
          residues.push({ name, num: r.num, icode: r.icode, ss: 'L', atoms: r.atoms });
        } else {
          model.ligands.push(Model.stamp(model, { name: r.name, num: r.num, icode: r.icode, chainId: cid, ss: 'L', atoms: r.atoms, het: true }));
        }
      }
      if (residues.length) Model.addChain(model, residues, cid.trim() || Model.newChainId(model));
    }
    // Secondary structure: from the file if it has any, else from backbone torsions.
    if (ssRanges && ssRanges.length) {
      for (const chain of model.chains) chain.residues.forEach(r => { r.ss = 'L'; });
      for (const s of ssRanges) {
        const chain = model.chains.find(c => c.id === s.chain);
        if (!chain) continue;
        let on = false;
        for (const r of chain.residues) {
          if (r.num === s.from && (r.icode || '') === (s.fromI || '')) on = true;
          if (on) r.ss = s.type;
          if (r.num === s.to && (r.icode || '') === (s.toI || '')) on = false;
          if (r.num > s.to) on = false;
        }
      }
    } else {
      for (const chain of model.chains) Model.autoSS(chain);
    }
    return model;
  }

  /* ---------------- PDB ---------------- */
  function parsePDB(text) {
    const records = [], ss = [];
    let modelSeen = false;
    for (const line of text.split(/\r?\n/)) {
      const rec = line.slice(0, 6);
      if (rec === 'MODEL ') { if (modelSeen) break; modelSeen = true; continue; }
      if (rec === 'ENDMDL') break;
      if (rec === 'HELIX ') {
        ss.push({ type: 'H', chain: line[19], from: parseInt(line.slice(21, 25)), fromI: line[25].trim(), to: parseInt(line.slice(33, 37)), toI: (line[37] || '').trim() });
        continue;
      }
      if (rec === 'SHEET ') {
        ss.push({ type: 'E', chain: line[21], from: parseInt(line.slice(22, 26)), fromI: line[26].trim(), to: parseInt(line.slice(33, 37)), toI: (line[37] || '').trim() });
        continue;
      }
      if (rec !== 'ATOM  ' && rec !== 'HETATM') continue;
      const name = line.slice(12, 16).trim();
      const resName = line.slice(17, 20).trim();
      let el = line.slice(76, 78).trim();
      if (!el) el = guessElement(line.slice(12, 16), resName);
      el = el.length > 1 ? el[0].toUpperCase() + el.slice(1).toLowerCase() : el.toUpperCase();
      records.push({
        het: rec === 'HETATM',
        name, resName,
        alt: line[16] && line[16] !== ' ' ? line[16] : '',
        chain: line[21] && line[21] !== ' ' ? line[21] : 'A',
        resSeq: parseInt(line.slice(22, 26)),
        icode: (line[26] || ' ').trim(),
        p: [parseFloat(line.slice(30, 38)), parseFloat(line.slice(38, 46)), parseFloat(line.slice(46, 54))],
        el,
      });
    }
    return assemble(records, ss);
  }

  /* ---------------- mmCIF ---------------- */
  function tokenize(line) {
    const out = [];
    let i = 0;
    while (i < line.length) {
      const c = line[i];
      if (c === ' ' || c === '\t') { i++; continue; }
      if (c === '#') break;
      if (c === "'" || c === '"') {
        let j = i + 1;
        while (j < line.length && !(line[j] === c && (j + 1 === line.length || line[j + 1] === ' ' || line[j + 1] === '\t'))) j++;
        out.push(line.slice(i + 1, j)); i = j + 1;
      } else {
        let j = i;
        while (j < line.length && line[j] !== ' ' && line[j] !== '\t') j++;
        out.push(line.slice(i, j)); i = j;
      }
    }
    return out;
  }

  /* Tables by category: { '_atom_site': { cols: [...], rows: [[...]] } }. Handles loop_ and
     single key-value items, quoted strings and ;-delimited text fields. */
  function parseCifTables(text) {
    const lines = text.split(/\r?\n/);
    const tables = {};
    let i = 0;
    const getTable = cat => tables[cat] || (tables[cat] = { cols: [], rows: [] });
    while (i < lines.length) {
      const line = lines[i];
      const t = line.trim();
      if (t === 'loop_') {
        i++;
        const cols = []; let cat = null;
        while (i < lines.length && lines[i].trim().startsWith('_')) {
          const full = lines[i].trim().split(/\s+/)[0];
          const dot = full.indexOf('.');
          cat = full.slice(0, dot); cols.push(full.slice(dot + 1)); i++;
        }
        const tab = getTable(cat); tab.cols = cols;
        let buf = [];
        while (i < lines.length) {
          const l = lines[i];
          const lt = l.trim();
          if (lt === 'loop_' || lt.startsWith('_') || lt.startsWith('data_') || (lt === '#' && buf.length === 0 && tab.rows.length)) break;
          if (l.startsWith(';')) {
            let txt = l.slice(1); i++;
            while (i < lines.length && !lines[i].startsWith(';')) { txt += '\n' + lines[i]; i++; }
            buf.push(txt.trim()); i++;
          } else { buf.push(...tokenize(l)); i++; }
          while (buf.length >= cols.length) { tab.rows.push(buf.slice(0, cols.length)); buf = buf.slice(cols.length); }
        }
        continue;
      }
      if (t.startsWith('_')) {
        const toks = tokenize(t);
        const full = toks[0]; const dot = full.indexOf('.');
        if (dot > 0) {
          const cat = full.slice(0, dot), col = full.slice(dot + 1);
          let val = toks[1];
          i++;
          if (val === undefined) {
            if (i < lines.length && lines[i].startsWith(';')) {
              let txt = lines[i].slice(1); i++;
              while (i < lines.length && !lines[i].startsWith(';')) { txt += '\n' + lines[i]; i++; }
              val = txt.trim(); i++;
            } else if (i < lines.length) { val = tokenize(lines[i])[0]; i++; }
          }
          const tab = getTable(cat);
          if (!tab.rows.length) tab.rows.push([]);
          tab.cols.push(col); tab.rows[0].push(val);
          continue;
        }
      }
      i++;
    }
    return tables;
  }

  function parseCIF(text) {
    const tables = parseCifTables(text);
    const at = tables._atom_site;
    if (!at) throw new Error('No _atom_site table found in the mmCIF file.');
    const ci = name => at.cols.indexOf(name);
    const pick = (row, ...names) => { for (const n of names) { const k = ci(n); if (k >= 0 && row[k] !== '?' && row[k] !== '.') return row[k]; } return null; };
    const records = [];
    let firstModel = null;
    for (const row of at.rows) {
      const mdl = pick(row, 'pdbx_PDB_model_num');
      if (firstModel === null) firstModel = mdl;
      if (mdl !== firstModel) continue;
      const resName = pick(row, 'auth_comp_id', 'label_comp_id') || 'UNK';
      const name = pick(row, 'auth_atom_id', 'label_atom_id') || 'X';
      let el = pick(row, 'type_symbol') || guessElement(name, resName);
      el = el.length > 1 ? el[0].toUpperCase() + el.slice(1).toLowerCase() : el.toUpperCase();
      records.push({
        het: pick(row, 'group_PDB') === 'HETATM',
        name, resName,
        alt: pick(row, 'label_alt_id') || '',
        chain: pick(row, 'auth_asym_id', 'label_asym_id') || 'A',
        resSeq: parseInt(pick(row, 'auth_seq_id', 'label_seq_id') || '0'),
        icode: pick(row, 'pdbx_PDB_ins_code') || '',
        p: [parseFloat(pick(row, 'Cartn_x')), parseFloat(pick(row, 'Cartn_y')), parseFloat(pick(row, 'Cartn_z'))],
        el,
      });
    }
    const ss = [];
    const rangeRows = (tab, type, filter) => {
      if (!tab) return;
      const c = n => tab.cols.indexOf(n);
      for (const row of tab.rows) {
        if (filter && !filter(row, c)) continue;
        const g = (a, b) => { const k = c(a) >= 0 ? c(a) : c(b); return k >= 0 ? row[k] : null; };
        const clean = v => (v === '?' || v === '.' || v === null) ? '' : v;
        ss.push({
          type, chain: g('beg_auth_asym_id', 'beg_label_asym_id'),
          from: parseInt(g('beg_auth_seq_id', 'beg_label_seq_id')), fromI: clean(g('pdbx_beg_PDB_ins_code')),
          to: parseInt(g('end_auth_seq_id', 'end_label_seq_id')), toI: clean(g('pdbx_end_PDB_ins_code')),
        });
      }
    };
    rangeRows(tables._struct_conf, 'H', (row, c) => c('conf_type_id') < 0 || /^HELX/.test(row[c('conf_type_id')]));
    rangeRows(tables._struct_sheet_range, 'E');
    return assemble(records, ss);
  }

  function parse(text, filename) {
    const isCif = /\.(cif|mmcif)$/i.test(filename || '') || /^data_/m.test(text.slice(0, 2000)) && /_atom_site\./.test(text);
    return isCif ? parseCIF(text) : parsePDB(text);
  }

  /* ---------------- export ---------------- */
  function fmtAtomName(name, el) {
    if (name.length >= 4) return name.slice(0, 4);
    return (el.length === 1 ? ' ' + name : name).padEnd(4);
  }

  /* opts: { backboneOnly, includeLigands } */
  function toPDB(model, opts) {
    opts = opts || {};
    const L = [];
    const f = (x, w, d) => x.toFixed(d).padStart(w);
    L.push('REMARK   1 WRITTEN BY PROTEIN CAD');
    // Secondary structure records.
    let hn = 0, sn = 0;
    for (const chain of model.chains) {
      const rs = chain.residues;
      for (let i = 0; i < rs.length;) {
        const t = rs[i].ss;
        let j = i; while (j + 1 < rs.length && rs[j + 1].ss === t) j++;
        const a = rs[i], b = rs[j], cid = chain.id.slice(0, 1);
        if (t === 'H' && j > i) {
          hn++;
          L.push('HELIX  ' + String(hn).padStart(3) + ' ' + String(hn).padStart(3) + ' ' + a.name.padEnd(3) + ' ' + cid + ' ' + String(a.num).padStart(4) + (a.icode || ' ') + ' ' +
            b.name.padEnd(3) + ' ' + cid + ' ' + String(b.num).padStart(4) + (b.icode || ' ') + ' 1'.padStart(2) + ' '.repeat(30) + String(j - i + 1).padStart(6));
        } else if (t === 'E' && j > i) {
          sn++;
          L.push('SHEET  ' + String(sn).padStart(3) + ' ' + ('S' + sn).padStart(3) + ' 1 ' + a.name.padEnd(3) + ' ' + cid + String(a.num).padStart(4) + (a.icode || ' ') + ' ' +
            b.name.padEnd(3) + ' ' + cid + String(b.num).padStart(4) + (b.icode || ' ') + ' 0');
        }
        i = j + 1;
      }
    }
    let serial = 1;
    const BBN = new Set(['N', 'CA', 'C', 'O']);
    const atomLine = (rec, a, res, cid) => rec + String(serial++ % 100000).padStart(5) + ' ' + fmtAtomName(a.name, a.el) + ' ' + res.name.slice(0, 3).padStart(3) + ' ' + cid +
      String(res.num).padStart(4) + (res.icode || ' ') + '   ' + f(a.p[0], 8, 3) + f(a.p[1], 8, 3) + f(a.p[2], 8, 3) + '  1.00  0.00          ' + a.el.toUpperCase().padStart(2);
    for (const chain of model.chains) {
      const cid = chain.id.slice(0, 1);
      let last = null;
      for (const res of chain.residues) {
        const std = !!Geom.THREE_TO_ONE[res.name];
        for (const a of orderAtoms(res)) {
          if (opts.backboneOnly && !BBN.has(a.name)) continue;
          L.push(atomLine(std ? 'ATOM  ' : 'HETATM', a, res, cid));
        }
        last = res;
      }
      if (last) L.push('TER   ' + String(serial++ % 100000).padStart(5) + '      ' + last.name.padStart(3) + ' ' + cid + String(last.num).padStart(4) + (last.icode || ' '));
    }
    if (opts.includeLigands !== false) {
      for (const res of model.ligands) for (const a of res.atoms) L.push(atomLine('HETATM', a, res, (res.chainId || 'Z').slice(0, 1)));
    }
    L.push('END');
    return L.join('\n') + '\n';
  }

  const ORDER = ['N', 'CA', 'C', 'O', 'CB'];
  function orderAtoms(res) {
    return res.atoms.slice().sort((a, b) => {
      const ia = ORDER.indexOf(a.name), ib = ORDER.indexOf(b.name);
      return (ia < 0 ? 99 : ia) - (ib < 0 ? 99 : ib);
    });
  }

  function cifStr(v) {
    v = String(v);
    if (v === '') return '?';
    return /[\s'"]/.test(v) || /^[_#$;]/.test(v) ? '"' + v + '"' : v;
  }

  function toCIF(model, opts) {
    opts = opts || {};
    const L = ['data_protein_cad', '#', 'loop_'];
    const cols = ['group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id', 'label_comp_id', 'label_asym_id', 'label_entity_id', 'label_seq_id',
      'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv', 'auth_seq_id', 'auth_comp_id', 'auth_asym_id', 'auth_atom_id', 'pdbx_PDB_model_num'];
    for (const c of cols) L.push('_atom_site.' + c);
    let serial = 1, entity = 0;
    const BBN = new Set(['N', 'CA', 'C', 'O']);
    for (const chain of model.chains) {
      entity++;
      chain.residues.forEach((res, k) => {
        const std = !!Geom.THREE_TO_ONE[res.name];
        for (const a of orderAtoms(res)) {
          if (opts.backboneOnly && !BBN.has(a.name)) continue;
          L.push([std ? 'ATOM' : 'HETATM', serial++, a.el.toUpperCase(), cifStr(a.name), '.', res.name, chain.id, entity, k + 1, res.icode || '?',
            a.p[0].toFixed(3), a.p[1].toFixed(3), a.p[2].toFixed(3), '1.00', '0.00', res.num, res.name, chain.id, cifStr(a.name), 1].join(' '));
        }
      });
    }
    if (opts.includeLigands !== false) {
      for (const res of model.ligands) {
        entity++;
        const cid = res.chainId || 'Z';
        for (const a of res.atoms) {
          L.push(['HETATM', serial++, a.el.toUpperCase(), cifStr(a.name), '.', res.name, cid + '_' + res.uid, entity, '.', res.icode || '?',
            a.p[0].toFixed(3), a.p[1].toFixed(3), a.p[2].toFixed(3), '1.00', '0.00', res.num, res.name, cid, cifStr(a.name), 1].join(' '));
        }
      }
    }
    L.push('#');
    return L.join('\n') + '\n';
  }

  window.IO = { parse, parsePDB, parseCIF, parseCifTables, toPDB, toCIF, tokenize };
})();
