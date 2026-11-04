/* Wiring: toolbar, keyboard, files, drag-and-drop, start-up. */
(function () {
  'use strict';
  const App = window.App, $ = App.$;

  function init() {
    try { App.initViewport(); } catch (e) {
      document.getElementById('viewport').innerHTML = '<div style="padding:30px;color:#ef5b5b">Could not start WebGL2: ' + e.message + '</div>';
      return;
    }
    App.initPanels();

    // Toolbar
    let openMode = 'replace';
    $('bNew').onclick = () => {
      if ((App.model.chains.length || App.model.ligands.length) && !confirm('Start a new empty model? (Undo can bring the current one back.)')) return;
      App.edit('New model', () => { App.model = Model.create(); App.sel.clear(); App.selOrder = []; });
      App.status('New model.');
    };
    $('bOpen').onclick = () => { openMode = 'replace'; $('fileInput').value = ''; $('fileInput').click(); };
    $('bImport').onclick = () => { openMode = 'add'; $('fileInput').value = ''; $('fileInput').click(); };
    $('fileInput').onchange = e => {
      const f = e.target.files[0]; if (!f) return;
      f.text().then(t => App.loadText(t, f.name, openMode === 'replace'));
    };
    const menu = $('exportMenu');
    $('bExport').onclick = e => {
      const r = e.target.getBoundingClientRect();
      menu.style.left = r.left + 'px'; menu.style.top = (r.bottom + 4) + 'px';
      menu.classList.toggle('open');
      e.stopPropagation();
    };
    menu.addEventListener('click', e => {
      e.stopPropagation();
      const b = e.target.closest('[data-export]');
      if (b) { App.exportModel(b.dataset.export); menu.classList.remove('open'); }
    });
    document.addEventListener('click', () => menu.classList.remove('open'));
    $('bUndo').onclick = App.undo;
    $('bRedo').onclick = App.redo;
    document.querySelectorAll('#toolButtons [data-tool]').forEach(b => { b.onclick = () => App.setTool(b.dataset.tool); });
    $('bConnect').onclick = () => App.connect();
    $('bRefine').onclick = () => App.refine();
    $('bDelete').onclick = () => App.deleteSelection();
    $('bFit').onclick = () => App.fitView(App.sel.size > 0);
    document.querySelectorAll('[data-view]').forEach(b => { b.onclick = () => App.setView(b.dataset.view); });
    $('bOrtho').onclick = () => {
      const c = App.renderer.cam; c.ortho = !c.ortho;
      $('bOrtho').textContent = c.ortho ? 'Ortho' : 'Persp';
      App.cameraChanged();
    };
    $('bHelp').onclick = () => { $('help').hidden = false; };
    $('helpClose').onclick = () => { $('help').hidden = true; };
    $('help').onclick = e => { if (e.target.id === 'help') $('help').hidden = true; };

    // Display options
    const disp = (id, key, prop) => {
      const el = $(id);
      el.addEventListener('change', () => {
        App.view[key] = el.type === 'checkbox' ? el.checked : el.value;
        if (key === 'clashes') App.derived.clashes = App.view.clashes ? Model.clashes(App.model) : [];
        App.projCache = null; App.sceneDirty = true; App.overlayDirty = true;
        App.refreshPanels();
      });
      void prop;
    };
    disp('dCartoon', 'cartoon'); disp('dAtoms', 'atoms'); disp('dStyle', 'style'); disp('dTrace', 'trace'); disp('dColor', 'color');
    disp('dBreaks', 'breaks'); disp('dClashes', 'clashes'); disp('dLabels', 'labels'); disp('dAxes', 'axes');

    // Keyboard
    window.addEventListener('keydown', e => {
      const tag = (e.target.tagName || '').toLowerCase();
      if (tag === 'input' || tag === 'select' || tag === 'textarea') return;
      const k = e.key.toLowerCase(), ctrl = e.ctrlKey || e.metaKey;
      if (ctrl && k === 'z' && !e.shiftKey) { App.undo(); e.preventDefault(); return; }
      if (ctrl && (k === 'y' || (k === 'z' && e.shiftKey))) { App.redo(); e.preventDefault(); return; }
      if (ctrl && k === 'a') { App.setSelection(App.model.chains.flatMap(c => c.residues.map(r => r.uid)), 'set'); e.preventDefault(); return; }
      if (ctrl && k === 'd') { App.duplicateSelection(); e.preventDefault(); return; }
      if (ctrl) return;
      if (k === 'delete' || k === 'backspace') { App.deleteSelection(); e.preventDefault(); return; }
      if (k === 'escape') {
        if (App.tool !== 'select') App.setTool('select'); else App.setSelection([], 'set');
        $('help').hidden = true; return;
      }
      const tools = { s: 'select', h: 'helix', b: 'sheet', d: 'draw' };
      if (tools[k]) { App.setTool(tools[k]); return; }
      if (k === 'c') { App.connect(); return; }
      if (k === 'r') { App.refine(); return; }
      if (k === 'f') { App.fitView(App.sel.size > 0); return; }
      if (k === 'l') { App.gizmo.space = App.gizmo.space === 'world' ? 'local' : 'world'; App.overlayDirty = true; App.refreshPanels(); App.status('Gizmo axes: ' + App.gizmo.space + '.'); return; }
      if (k === '?') { $('help').hidden = false; }
    });

    // Drag and drop files: replaces the model unless Shift is held (then adds).
    window.addEventListener('dragover', e => { e.preventDefault(); });
    window.addEventListener('drop', e => {
      e.preventDefault();
      const f = e.dataTransfer.files[0]; if (!f) return;
      const add = e.shiftKey || App.model.chains.length > 0 && confirm('Add ' + f.name + ' to the current model?\n\nOK = add its chains, Cancel = replace the model.');
      f.text().then(t => App.loadText(t, f.name, !add));
    });

    // Start-up
    const restored = App.restore();
    App.changed();
    App.setView('iso');
    App.fitView();
    App.setTool('select');
    App.status(restored ? 'Restored the model from your last session (stored in this browser). New starts over.' : 'Ready. Pick the Helix tool (H) and click the grid to start, or Open a PDB / mmCIF file.');
    requestAnimationFrame(App.frame);
  }

  window.addEventListener('DOMContentLoaded', init);
})();
