/* Copyright (c) 2026 Neil Yuanting Li. MIT License; see LICENSE. */
(function () {
  'use strict';
  const P = window.TransportPhysics;
  const $ = id => document.getElementById(id);
  const blue = '#145b91', red = '#9b392b', gray = '#666';
  const records = {};
  const downloads = {};
  const defaults = {
    level: { epsilon: 0, gammaL: .5, gammaR: .5, probe: 0 },
    chain: { sites: 5, hopping: 1, contact: 1, onsite: 0, probe: .3 },
    fano: { sideEnergy: 0, coupling: .5, probe: 0 },
    bias: { epsilon: 0, bias: 1, kT: .025, fraction: .5 },
    eta: { gamma: .02, eta: .001, limit: 4, intervals: 1600, alignment: 0 }
  };
  const state = Object.fromEntries(Object.entries(defaults).map(([key, value]) => [key, {...value}]));
  const presets = {
    level: {
      matched: {...defaults.level}, narrow: {epsilon: 0, gammaL: .05, gammaR: .05, probe: 0},
      unequal: {epsilon: 0, gammaL: .1, gammaR: .9, probe: 0}
    },
    chain: {
      matched: {...defaults.chain}, weak: {sites: 5, hopping: 1, contact: .35, onsite: 0, probe: .3},
      barrier: {sites: 5, hopping: 1, contact: 1, onsite: 1, probe: -1.5}
    },
    fano: {
      zero: {...defaults.fano}, shift: {sideEnergy: .7, coupling: .35, probe: .7},
      off: {sideEnergy: 0, coupling: 0, probe: 0}
    },
    bias: {
      equilibrium: {...defaults.bias, bias: 0}, flow: {...defaults.bias, bias: 2},
      asymmetric: {...defaults.bias, bias: 3, fraction: .1}
    },
    eta: {
      resolved: {...defaults.eta},
      absorber: {gamma: .02, eta: .05, limit: 4, intervals: 1600, alignment: 0},
      missed: {gamma: .002, eta: .0001, limit: 4, intervals: 400, alignment: .5},
      isolated: {gamma: 0, eta: .001, limit: 4, intervals: 1600, alignment: 0}
    }
  };
  const fmt = (x, digits = 4) => {
    if (x === null || x === undefined || !Number.isFinite(x)) return '—';
    if (x === 0) return '0';
    return Math.abs(x) < .0001 || Math.abs(x) >= 1e5 ? x.toExponential(2) : Number(x.toFixed(digits)).toString();
  };
  const text = (x,y,s,extra='') => '<text x="'+x+'" y="'+y+'" '+extra+'>'+s+'</text>';
  const line = (x1,y1,x2,y2,extra='') => '<line x1="'+x1+'" y1="'+y1+'" x2="'+x2+'" y2="'+y2+'" '+extra+'/>';
  function sample(fn, lo, hi, count = 700, extra = []) {
    const energies = [...Array.from({length: count + 1}, (_,i) => lo + (hi-lo)*i/count), ...extra.filter(x => x>=lo && x<=hi)];
    return [...new Set(energies)].sort((a,b)=>a-b).map(E => [E, fn(E)]);
  }
  function plot(id, curves, options = {}) {
    const svg = $(id), W = 680, H = 340, left = 62, right = 18, top = 20, bottom = 50;
    const {xmin=-3,xmax=3,ymin=0,ymax=1.05,ylabel='Probability',probe=null,verticals=[],dots=[]} = options;
    const x = E => left + (E-xmin)/(xmax-xmin)*(W-left-right);
    const y = value => H-bottom-(value-ymin)/(ymax-ymin)*(H-top-bottom);
    const clip = id+'-clip';
    let html = '<title>'+svg.getAttribute('aria-label')+'</title><rect width="680" height="340" fill="white"/><defs><clipPath id="'+clip+'"><rect x="'+left+'" y="'+top+'" width="'+(W-left-right)+'" height="'+(H-top-bottom)+'"/></clipPath></defs><g font-family="Times New Roman,serif" font-size="16" fill="#111">';
    for (let i=0;i<=4;i++) {
      const value=ymin+(ymax-ymin)*i/4, py=y(value);
      html+=line(left,py,W-right,py,'stroke="#ddd"')+text(left-8,py+5,fmt(value,2),'text-anchor="end"');
    }
    for (let i=0;i<=6;i++) {
      const E=xmin+(xmax-xmin)*i/6, px=x(E);
      html+=line(px,H-bottom,px,H-bottom+5,'stroke="#111"')+text(px,H-bottom+24,fmt(E,2),'text-anchor="middle"');
    }
    html+=line(left,top,left,H-bottom,'stroke="#111"')+line(left,H-bottom,W-right,H-bottom,'stroke="#111"');
    html+=text((left+W-right)/2,H-7,'Energy E (eV)','text-anchor="middle"');
    html+=text(17,(top+H-bottom)/2,ylabel,'text-anchor="middle" transform="rotate(-90 17 '+((top+H-bottom)/2)+')"');
    html+='<g clip-path="url(#'+clip+')">';
    for (const E of verticals) html+=line(x(E),top,x(E),H-bottom,'stroke="#777" stroke-dasharray="4 5"');
    if (probe!==null) html+=line(x(probe),top,x(probe),H-bottom,'stroke="#888"');
    for (const curve of curves) {
      let d='', inPath=false;
      for (const [E,value] of curve.data) {
        if (!Number.isFinite(value)) {inPath=false;continue;}
        d+=(inPath?'L':'M')+x(E).toFixed(3)+','+y(value).toFixed(3); inPath=true;
      }
      html+='<path d="'+d+'" fill="none" stroke="'+curve.color+'" stroke-width="2"'+(curve.dash?' stroke-dasharray="6 4"':'')+'/>';
    }
    for (const [E,value] of dots) if (E>=xmin && E<=xmax && Number.isFinite(value)) html+='<circle cx="'+x(E)+'" cy="'+y(value)+'" r="2.5" fill="#111"/>';
    html+='</g></g>'; svg.innerHTML=html;
  }
  function schematic(id, mode, p) {
    let s='<title>'+$(id).getAttribute('aria-label')+'</title><g font-family="Times New Roman,serif" font-size="18" fill="#111" stroke="#111" stroke-width="1.5">';
    s+='<rect x="14" y="28" width="93" height="43" fill="white"/><rect x="573" y="28" width="93" height="43" fill="white"/>';
    s+=text(60,55,'L','text-anchor="middle" stroke="none"')+text(620,55,'R','text-anchor="middle" stroke="none"');
    s+=line(107,49,573,49);
    if(mode==='chain') {
      for(let i=0;i<p.sites;i++) {
        const x=160+i*360/Math.max(1,p.sites-1);
        s+='<circle cx="'+(p.sites===1?340:x)+'" cy="49" r="8" fill="white"/>';
      }
      s+=text(340,91,'N = '+p.sites+'; t = '+fmt(p.hopping)+' eV; τ = '+fmt(p.contact)+' eV','text-anchor="middle" stroke="none"');
    } else {
      s+='<circle cx="340" cy="49" r="13" fill="white"/>';
      if(mode==='fano') {
        s+=line(340,62,340,108)+'<circle cx="340" cy="113" r="10" fill="white"/>';
        s+=text(363,119,'εs = '+fmt(p.sideEnergy)+' eV','stroke="none"')+text(360,84,'v = '+fmt(p.coupling)+' eV','stroke="none"');
      } else {
        s+=text(340,91,'ε₀ = '+fmt(p.epsilon)+' eV','text-anchor="middle" stroke="none"');
        s+=text(209,32,'ΓL = '+fmt(p.gammaL)+' eV','text-anchor="middle" stroke="none"');
        s+=text(477,32,'ΓR = '+fmt(p.gammaR)+' eV','text-anchor="middle" stroke="none"');
      }
    }
    $(id).innerHTML=s+'</g>';
  }
  function slider(group,key,label,min,max,step,unit='') {
    const id=group+'-'+key;
    return '<div class="control"><label for="'+id+'">'+label+'<output for="'+id+'" id="'+id+'-value"></output></label><input id="'+id+'" type="range" min="'+min+'" max="'+max+'" step="'+step+'" data-group="'+group+'" data-key="'+key+'"><span class="small">'+min+' to '+max+(unit?' '+unit:'')+'</span></div>';
  }
  function select(group,key,label,items) {
    const id=group+'-'+key;
    return '<div class="control"><label for="'+id+'">'+label+'</label><select id="'+id+'" data-group="'+group+'" data-key="'+key+'">'+items.map(([value,title])=>'<option value="'+value+'">'+title+'</option>').join('')+'</select></div>';
  }
  function sync(group) {
    for (const [key,value] of Object.entries(state[group])) {
      const control=$(group+'-'+key); if(!control)continue;
      control.value=value;
      const output=$(group+'-'+key+'-value'); if(output)output.value=fmt(value);
    }
  }
  function readout(group,paragraphs) { $(group+'-readout').innerHTML=paragraphs.map(s=>'<p>'+s+'</p>').join(''); }
  function record(group,parameters,results,curves) {
    records[group]={format:'quantum-transport-notes/1',version:'1.0.0',experiment:group,
      units:{energy:'eV',current:'microampere',time:'fs',spectralDensity:'1/eV'},
      conventions:{spinSpecies:1,retardedSign:'E + i eta',leadEta:0,leadCurrent:'positive entering device',model:'noninteracting, coherent, elastic; finite eta interpreted as empty absorber'},
      parameters,results,curves};
  }
  function renderLevel() {
    const p=state.level, values=sample(E=>P.resonantLevel(E,p),-3,3,900,[p.epsilon,p.probe]);
    const curves=[{color:blue,data:values.map(([E,r])=>[E,r.T])},{color:red,dash:true,data:values.map(([E,r])=>[E,r.R])}];
    plot('level-plot',curves,{probe:p.probe}); schematic('level-schematic','level',p);
    const r=P.resonantLevel(p.probe,p), peak=P.resonantLevel(p.epsilon,p), gamma=p.gammaL+p.gammaR;
    readout('level',[
      '<strong>At E = '+fmt(p.probe)+' eV: T = '+fmt(r.T)+', R = '+fmt(r.R)+'.</strong> T + R = '+fmt(r.T+r.R,10)+'.',
      'Peak T = '+fmt(peak.T)+'; full width at half maximum Γ = '+fmt(gamma)+' eV. Probability lifetime ℏ/Γ = '+fmt(P.constants.hbarEVfs/gamma)+' fs; amplitude e-folding time = '+fmt(2*P.constants.hbarEVfs/gamma)+' fs.',
      'If the probe energy is the Fermi energy, the zero-temperature small-bias conductance is '+fmt(P.constants.currentMicroampPerEV*r.T)+' μS per spin species.'
    ]); record('level',{...p,eta:0},r,curves);
  }
  function renderChain() {
    const p=state.chain, values=window.TransportSampling.chainCurve({...p,eta:0});
    const curves=[{color:blue,data:values}];
    plot('chain-plot',curves,{xmin:-2.4,xmax:2.4,probe:p.probe,verticals:[-2,2]});
    schematic('chain-schematic','chain',p);
    const r=P.chain(p.probe,p), s=r.selfEnergy;
    readout('chain',[
      '<strong>At E = '+fmt(p.probe)+' eV: T = '+fmt(r.T)+'.</strong> '+(r.channels?'R = '+fmt(r.R)+'; T + R = '+fmt(r.T+r.R,10)+'.':'No propagating lead channel at this energy; an incident-flux reflection probability is not defined.'),
      'Each contact: Σʳ = '+fmt(s.re)+' − i '+fmt(-s.im)+' eV; Γ = −2 Im Σʳ = '+fmt(s.gamma)+' eV.',
      'The finite device has '+p.sites+' site'+(p.sites>1?'s':'')+'. Leads remain semi-infinite, with hopping 1 eV and on-site energy 0.'
    ]); record('chain',{...p,eta:0,leadHopping:1,leadOnsite:0},r,curves);
  }
  function renderFano() {
    const p={...state.fano,epsilon:0,gammaL:.5,gammaR:.5,eta:0};
    // The antiresonance narrows quadratically with v. Resolve its neighborhood
    // explicitly instead of drawing a falsely broad notch between coarse nodes.
    const localWidth=p.coupling*p.coupling*.5/(p.sideEnergy*p.sideEnergy+.25);
    const local=p.coupling>0?Array.from({length:801},(_,i)=>p.sideEnergy+(i-400)*localWidth/20):[];
    const values=sample(E=>P.fano(E,p).T,-3,3,900,[p.sideEnergy,p.probe,...local]);
    const curves=[{color:blue,data:values},{color:gray,dash:true,data:sample(E=>P.resonantLevel(E,p).T,-3,3)}];
    plot('fano-plot',curves,{probe:p.probe}); schematic('fano-schematic','fano',p);
    const r=P.fano(p.probe,p), atSide=P.fano(p.sideEnergy,p);
    readout('fano',[
      '<strong>At E = '+fmt(p.probe)+' eV: T = '+fmt(r.T)+', R = '+fmt(r.R)+'.</strong>',
      p.coupling>0 ? 'At the side-orbital energy E = εs = '+fmt(p.sideEnergy)+' eV, T = '+fmt(atSide.T)+'. The full incident flux is reflected; no probability is absorbed.' : 'The side orbital is disconnected. Its bound state does not enter the transmission. Its occupation requires an initial condition, and its delta-function spectral weight is not drawn.',
      'Move εs while keeping the conducting orbital fixed: the zero follows εs. Changing v controls the width and hybridized resonances.'
    ]); record('fano',p,r,curves);
  }
  function renderBias() {
    const s=state.bias, p={epsilon:s.epsilon,gammaL:.4*s.fraction,gammaR:.4*(1-s.fraction),bias:s.bias,kT:s.kT,eta:0,limit:40,intervals:16000};
    const muL=p.bias/2,muR=-p.bias/2;
    const energies=sample(E=>E,-3,3,900,[muL,muR,p.epsilon]).map(([E])=>E);
    const curves=[
      {color:blue,data:energies.map(E=>[E,P.resonantLevel(E,p).T])},
      {color:red,data:energies.map(E=>[E,P.fermi(E,muL,p.kT)])},
      {color:'#555',dash:true,data:energies.map(E=>[E,P.fermi(E,muR,p.kT)])}
    ];
    plot('bias-plot',curves,{ylabel:'T or occupation',verticals:[muL,muR]});
    const r=P.bias(p), fine=P.bias({...p,intervals:32000});
    const half=.2, exactArea=(Math.atan((40-p.epsilon)/half)-Math.atan((-40-p.epsilon)/half))/Math.PI;
    readout('bias',[
      '<strong>I = '+fmt(fine.I,5)+' μA; n = '+fmt(fine.n,5)+' electron per spin orbital.</strong> μL = '+fmt(muL)+' eV; μR = '+fmt(muR)+' eV.',
      'ΓL = '+fmt(p.gammaL)+' eV, ΓR = '+fmt(p.gammaR)+' eV. IL + IR = '+fmt(fine.IL+fine.IR)+' μA. Temperature = '+fmt(p.kT/P.constants.kBoltzmannEV,1)+' K.',
      'Integration: −40 to +40 eV, 32,000 base intervals (ΔE = 0.0025 eV); at zero temperature integration also splits at the chemical potentials. Halving ΔE changes I by '+fmt(Math.abs(fine.I-r.I),7)+' μA and n by '+fmt(Math.abs(fine.n-r.n),7)+'. The finite window omits '+fmt(1-exactArea,5)+' of the total spectral weight; this bounds the omitted occupation, independently of grid error. The stated n is the finite-window result.'
    ]); record('bias',{...p,intervals:32000}, {...fine,coarser:r,exactWindowSpectralArea:exactArea},curves);
  }
  // Independent closed-form zero-temperature integrals make the sink-current
  // diagnostic insensitive to the intentionally coarse spectral quadrature.
  function exactSinkCurrents(gamma,eta) {
    const width=gamma+2*eta, half=width/2, factor=P.constants.currentMicroampPerEV;
    const occupiedIntegral=mu=>(Math.atan(mu/half)+Math.PI/2)/half;
    const left=occupiedIntegral(.5),right=occupiedIntegral(-.5);
    const pair=factor*gamma*gamma/4*(left-right);
    const sinkL=factor*2*eta*gamma/2*left,sinkR=factor*2*eta*gamma/2*right;
    return {I:pair,IL:pair+sinkL,IR:-pair+sinkR,Isink:sinkL+sinkR,bias:1,kT:0,integration:'analytic infinite window'};
  }
  function renderEta() {
    const s=state.eta, p={epsilon:0,gammaL:s.gamma/2,gammaR:s.gamma/2,eta:s.eta,limit:s.limit,intervals:s.intervals,offset:s.alignment*2*s.limit/s.intervals};
    if(s.gamma+2*s.eta===0) {
      readout('eta',['<strong>An isolated level at η = 0 has a delta-function spectrum.</strong> Choose a positive η to sample it. No ordinary finite curve or uniform-grid integral is defined at its pole.']);
      $('eta-plot').innerHTML='<text x="30" y="160" font-family="serif" font-size="20">δ(E): choose η &gt; 0 to draw a regularized spectrum.</text>';
      $('eta-table').innerHTML=''; record('eta',p,{distribution:'delta(E)'},[]); return;
    }
    const d=P.etaDiagnostics(p), half=d.width/2;
    const extra=Array.from({length:401},(_,i)=>(i-200)/20*half).concat([0]);
    const values=sample(E=>P.resonantLevel(E,p),-.15,.15,800,extra);
    const curves=[{color:blue,data:values.map(([E,r])=>[E,r.A/(2*Math.PI)])},
      {color:red,dash:true,data:values.map(([E,r])=>[E,(r.AL+r.AR)/(2*Math.PI)])}];
    const dots=[];
    for(let i=0;i<=s.intervals;i++) {
      const E=d.lower+i*d.gridStep;
      if(Math.abs(E)<=.15) dots.push([E,P.resonantLevel(E,p).A/(2*Math.PI)]);
    }
    plot('eta-plot',curves,{xmin:-.15,xmax:.15,ymax:d.peak.A/(2*Math.PI)*1.06,ylabel:'Spectral density (1/eV)',dots});
    const currents=exactSinkCurrents(s.gamma,s.eta), ratio=s.gamma/d.width;
    readout('eta',[
      '<strong>Sampled full area = '+fmt(d.area,7)+'; exact same-window area = '+fmt(d.analyticArea,7)+'.</strong>',
      'Quadrature error = '+fmt(d.quadratureError,7)+'; omitted spectral tail = '+fmt(d.omittedTail,7)+'. ΔE = '+fmt(d.gridStep,6)+' eV; full width Γ + 2η = '+fmt(d.width,6)+' eV. '+(d.gridStep>half?'The integration spacing exceeds the half-width; inspect the alignment and refinement table.':'Refine and shift the grid to test resolution.'),
      'At resonance for left incidence: T = '+fmt(d.peakT)+', R = '+fmt(d.peakR)+', absorbed probability = '+fmt(d.peakLoss)+'. Their sum is '+fmt(d.peakT+d.peakR+d.peakLoss,10)+'. The contact-only spectral area on the infinite axis is '+fmt(ratio)+'.',
      'At V = 1 V and temperature zero (analytic integrals): IL = '+fmt(currents.IL)+' μA; IR = '+fmt(currents.IR)+' μA; IL + IR = Isink = '+fmt(currents.Isink)+' μA. Both lead currents are defined entering the device. '+(s.gamma===0?'There is no physical-lead injection.':'The pair current I = '+fmt(currents.I)+' μA alone hides this sink.'),
      s.eta>0 ? 'If η is deliberately interpreted as absorption, its additional probability e-folding time is ℏ/(2η) = '+fmt(P.constants.hbarEVfs/(2*s.eta))+' fs. This is not a lifetime measurement.' : 'η = 0 adds no artificial absorption; any linewidth here comes from physical contacts.'
    ]);
    const table=[];
    for(const eta of [...new Set([s.eta,s.eta/10,0])]) for(const multiplier of [1,2,4]) {
      if(s.gamma+2*eta===0)continue;
      const n=s.intervals*multiplier;
      const result=P.etaDiagnostics({...p,eta,intervals:n,offset:s.alignment*2*s.limit/n});
      table.push({eta,intervals:n,...result});
    }
    $('eta-table').innerHTML='<table><caption>Recomputed convergence, same fractional grid alignment</caption><thead><tr><th>η (eV)</th><th>Intervals</th><th>ΔE (eV)</th><th>Sampled area</th><th>Exact window</th><th>Grid error</th></tr></thead><tbody>'+table.map(r=>'<tr><td>'+fmt(r.eta)+'</td><td>'+r.intervals+'</td><td>'+fmt(r.gridStep,6)+'</td><td>'+fmt(r.area,6)+'</td><td>'+fmt(r.analyticArea,6)+'</td><td>'+fmt(r.quadratureError,6)+'</td></tr>').join('')+'</tbody></table>';
    record('eta',p,{...d,currents,contactOnlyInfiniteArea:ratio,convergence:table,integrationSamplesInPlot:dots},curves);
  }
  const renderers={level:renderLevel,chain:renderChain,fano:renderFano,bias:renderBias,eta:renderEta};
  function render(group) {
    try {sync(group);renderers[group]();$('app-error').hidden=true;}
    catch(error) {$('app-error').textContent='Calculation could not finish: '+error.message;$('app-error').hidden=false;console.error(error);}
  }
  function prepareDownload(group,kind) {
    const r=records[group]; if(!r)return;
    let content;
    if(kind==='json') content=JSON.stringify(r,null,2);
    else {
      const clone=$(group+'-plot').cloneNode(true), ns='http://www.w3.org/2000/svg';
      clone.setAttribute('viewBox','0 0 680 440');clone.setAttribute('width','680');clone.setAttribute('height','440');
      const metadata=document.createElementNS(ns,'metadata');
      metadata.textContent=JSON.stringify({experiment:group,parameters:r.parameters,units:r.units,conventions:r.conventions});
      clone.prepend(metadata);
      const legends={level:['Transmission T','Reflection R'],chain:['Transmission'],
        fano:['With side orbital','Without side orbital'],bias:['Transmission T','fL','fR'],
        eta:['Full A / 2π','Contact part / 2π']};
      const p=r.parameters;
      const captions={
        level:'ε₀ = '+fmt(p.epsilon)+' eV; ΓL = '+fmt(p.gammaL)+' eV; ΓR = '+fmt(p.gammaR)+' eV; η = 0',
        chain:'N = '+p.sites+'; t = '+fmt(p.hopping)+' eV; τ = '+fmt(p.contact)+' eV; εd = '+fmt(p.onsite)+' eV; η = 0',
        fano:'εs = '+fmt(p.sideEnergy)+' eV; v = '+fmt(p.coupling)+' eV; ΓL = ΓR = 0.5 eV; η = 0',
        bias:'ε₀ = '+fmt(p.epsilon)+' eV; V = '+fmt(p.bias)+' V; kBT = '+fmt(p.kT)+' eV; ΓL / Γ = '+fmt(p.gammaL/.4),
        eta:'Γ = '+fmt(p.gammaL+p.gammaR)+' eV; η = '+fmt(p.eta)+' eV; window half-width = '+p.limit+' eV'
      };
      const footer=document.createElementNS(ns,'g');footer.setAttribute('font-family','Times New Roman,serif');footer.setAttribute('font-size','15');footer.setAttribute('fill','#111');
      let s='<rect x="0" y="340" width="680" height="100" fill="white"/>';
      (legends[group]||[]).forEach((label,i)=>{
        const x=34+i*215,curve=r.curves[i]; if(!curve)return;
        s+=line(x,361,x+25,361,'stroke="'+curve.color+'" stroke-width="2"'+(curve.dash?' stroke-dasharray="5 3"':''));
        s+=text(x+31,366,label);
      });
      if(group==='eta')s+='<circle cx="477" cy="361" r="2.5"/>'+text(489,366,'Integration samples');
      s+=text(34,397,captions[group]);
      s+=text(34,422,'One spin species. Copyright (c) 2026 Neil Yuanting Li. MIT License.');
      footer.innerHTML=s;clone.append(footer);
      content=new XMLSerializer().serializeToString(clone);
    }
    const mime=kind==='json'?'application/json':'image/svg+xml';
    const key=group+kind;if(downloads[key])URL.revokeObjectURL(downloads[key]);
    downloads[key]=URL.createObjectURL(new Blob([content],{type:mime}));
    const anchor=$(group+'-download'); anchor.href=downloads[key];anchor.download=group+'.'+kind;anchor.textContent='Save '+group+'.'+kind;anchor.hidden=false;anchor.click();
    const source=$(group+'-export-source');source.hidden=false;
    source.querySelector('summary').textContent='Inspect this '+kind.toUpperCase()+' export (snapshot of the settings above)';
    source.querySelector('pre').textContent=content;
  }
  try {
    if(!P)throw new Error('physics.js did not load');
    $('level-controls').innerHTML='<h3>Change the orbital</h3>'+slider('level','epsilon','Level ε₀',-1.5,1.5,.05,'eV')+slider('level','gammaL','Left width ΓL',.05,1,.05,'eV')+slider('level','gammaR','Right width ΓR',.05,1,.05,'eV')+slider('level','probe','Probe energy',-2.5,2.5,.05,'eV');
    $('chain-controls').innerHTML='<h3>Change the chain</h3>'+slider('chain','sites','Sites N',1,9,1)+slider('chain','hopping','Internal hopping t',.4,1.4,.05,'eV')+slider('chain','contact','Contact hopping τ',.15,1.4,.05,'eV')+slider('chain','onsite','Device on-site εd',-1.5,1.5,.1,'eV')+slider('chain','probe','Probe energy',-2.4,2.4,.05,'eV');
    $('fano-controls').innerHTML='<h3>Change the detour</h3>'+slider('fano','sideEnergy','Side energy εs',-1.5,1.5,.05,'eV')+slider('fano','coupling','Coupling v',0,1,.025,'eV')+slider('fano','probe','Probe energy',-2.5,2.5,.05,'eV');
    $('bias-controls').innerHTML='<h3>Fill the reservoirs</h3>'+slider('bias','epsilon','Level ε₀',-2,2,.05,'eV')+slider('bias','bias','Bias V',-4,4,.1,'V')+slider('bias','kT','Thermal energy kBT',0,.2,.005,'eV')+slider('bias','fraction','ΓL / (ΓL + ΓR)',.05,.95,.05)+'<p class="small">Total physical width ΓL + ΓR = 0.4 eV. Numerical η = 0.</p>';
    $('eta-controls').innerHTML='<h3>Audit the numerics</h3>'+
      select('eta','gamma','Total physical width Γ (eV)',[[0,'0 (isolated)'],[.002,'0.002'],[.02,'0.02'],[.2,'0.2'],[.6,'0.6']])+
      select('eta','eta','Numerical η (eV)',[[0,'0'],[.00001,'0.00001'],[.0001,'0.0001'],[.001,'0.001'],[.01,'0.01'],[.05,'0.05'],[.1,'0.1']])+
      select('eta','limit','Window half-width (eV)',[[1,'1'],[4,'4'],[8,'8'],[20,'20']])+
      select('eta','intervals','Uniform grid intervals',[[200,'200'],[400,'400'],[800,'800'],[1600,'1,600'],[6400,'6,400'],[25600,'25,600']])+
      select('eta','alignment','Grid alignment',[[0,'Peak on a grid point'],[.5,'Shift by half a grid step']])+
      '<p class="small">ε₀ = 0; ΓL = ΓR = Γ/2. Shifting the grid also shifts the two window endpoints; the exact area uses those actual endpoints.</p>';
    for (const group of Object.keys(state)) {
      $(group+'-controls').insertAdjacentHTML('beforeend','<button data-reset="'+group+'">Reset this experiment</button>');
      document.querySelector('[data-export="'+group+'"]').innerHTML='<button data-save="'+group+':svg">Save figure (SVG)</button><button data-save="'+group+':json">Save calculation (JSON)</button><a id="'+group+'-download" hidden></a>';
      document.querySelector('[data-export="'+group+'"]').insertAdjacentHTML('afterend','<details id="'+group+'-export-source" hidden><summary>Inspect export</summary><pre tabindex="0" aria-label="Export source"></pre></details>');
    }
    document.addEventListener('input',event=>{
      const {group,key}=event.target.dataset;
      if(!group || !key)return;
      state[group][key]=Number(event.target.value);render(group);
    });
    document.addEventListener('click',event=>{
      const target=event.target.closest('button');if(!target)return;
      if(target.dataset.preset) {
        const [group,name]=target.dataset.preset.split(':');state[group]={...presets[group][name]};render(group);
      } else if(target.dataset.reset) {
        const group=target.dataset.reset;state[group]={...defaults[group]};render(group);
      } else if(target.dataset.save) {
        prepareDownload(...target.dataset.save.split(':'));
      }
    });
    for(const group of Object.keys(state))render(group);
  } catch(error) {
    $('app-error').textContent='The experiments could not start: '+error.message;$('app-error').hidden=false;console.error(error);
  }
})();
