/* Copyright (c) 2026 Neil Yuanting Li. MIT License; see LICENSE. */
(async function () {
  'use strict';
  const P = window.TransportPhysics;
  const $ = id => document.getElementById(id);
  const blue = '#145b91', red = '#9b392b', gray = '#666';
  const records = {};
  const downloads = {};
  const errors = {};
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
  // All displayed mathematical notation is authored as TeX. KaTeX handles
  // page text; MathJax supplies self-contained vector paths for SVG figures.
  const tex = String.raw;
  const numberTex = (x, digits=4) => {
    if (!Number.isFinite(x)) return String.raw`\text{not defined}`;
    const value=fmt(x,digits), parts=value.split('e');
    return parts.length===2 ? parts[0]+String.raw`\times10^{`+Number(parts[1])+'}' : value;
  };
  const math = (strings,...values) => window.katex.renderToString(String.raw(strings,...values),{throwOnError:true,strict:'error',trust:false});
  const mathNumber = (x,digits=4) => math`${numberTex(x,digits)}`;
  const svgMath = (x,y,source,options={}) => window.TransportSvgMath.render(source,x,y,options);
  const text = (x,y,s,extra='') => '<text x="'+x+'" y="'+y+'" '+extra+'>'+s+'</text>';
  const line = (x1,y1,x2,y2,extra='') => '<line x1="'+x1+'" y1="'+y1+'" x2="'+x2+'" y2="'+y2+'" '+extra+'/>';
  function sample(fn, lo, hi, count = 700, extra = []) {
    const energies = [...Array.from({length: count + 1}, (_,i) => lo + (hi-lo)*i/count), ...extra.filter(x => x>=lo && x<=hi)];
    return [...new Set(energies)].sort((a,b)=>a-b).map(E => [E, fn(E)]);
  }
  function plot(id, curves, options = {}) {
    const spectralPlot = options.density === true;
    const svg = $(id), W = 680, H = 340, left = spectralPlot ? 90 : 62, right = 18, top = 20, bottom = 50;
    const {xmin=-3,xmax=3,ymin=0,ymax=1.05,ylabelTex=tex`\text{Probability}`,probe=null,verticals=[],dots=[]} = options;
    const ticks = options.yticks || (ymin===0 && ymax===1.05 ? [0,.25,.5,.75,1] : Array.from({length:5},(_,i)=>ymin+(ymax-ymin)*i/4));
    const x = E => left + (E-xmin)/(xmax-xmin)*(W-left-right);
    const y = value => H-bottom-(value-ymin)/(ymax-ymin)*(H-top-bottom);
    const clip = id+'-clip';
    let html = '<title>'+svg.getAttribute('aria-label')+'</title><rect width="680" height="340" fill="white"/><defs><clipPath id="'+clip+'"><rect x="'+left+'" y="'+top+'" width="'+(W-left-right)+'" height="'+(H-top-bottom)+'"/></clipPath></defs><g font-family="Times New Roman,serif" font-size="16" fill="#111">';
    for (const value of ticks) {
      const py=y(value);
      html+=line(left,py,W-right,py,'stroke="#ddd"')+svgMath(left-8,py+5,numberTex(spectralPlot?Number(value.toPrecision(3)):value,2),{anchor:'end'});
    }
    for (let i=0;i<=6;i++) {
      const E=xmin+(xmax-xmin)*i/6, px=x(E);
      html+=line(px,H-bottom,px,H-bottom+5,'stroke="#111"')+svgMath(px,H-bottom+24,numberTex(E,2),{anchor:'middle'});
    }
    html+=line(left,top,left,H-bottom,'stroke="#111"')+line(left,H-bottom,W-right,H-bottom,'stroke="#111"');
    html+=svgMath((left+W-right)/2,H-7,tex`\text{Energy }E\;(\mathrm{eV})`,{anchor:'middle'});
    html+=svgMath(17,(top+H-bottom)/2,ylabelTex,{anchor:'middle',rotate:-90});
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
    if(mode==='fano') $(id).setAttribute('aria-label',p.coupling>0?'A side orbital attached to the conducting level':'A side orbital disconnected from the conducting level');
    let s='<title>'+$(id).getAttribute('aria-label')+'</title><g font-family="Times New Roman,serif" font-size="18" fill="#111" stroke="#111" stroke-width="1.5">';
    s+='<rect x="14" y="28" width="93" height="43" fill="white"/><rect x="573" y="28" width="93" height="43" fill="white"/>';
    s+=svgMath(60,55,tex`\mathrm{L}`,{anchor:'middle',fontSize:18})+svgMath(620,55,tex`\mathrm{R}`,{anchor:'middle',fontSize:18});
    s+=line(107,49,573,49);
    if(mode==='chain') {
      for(let i=0;i<p.sites;i++) {
        const x=160+i*360/Math.max(1,p.sites-1);
        s+='<circle cx="'+(p.sites===1?340:x)+'" cy="49" r="8" fill="white"/>';
      }
      s+=svgMath(340,91,tex`N=${p.sites};\quad t=${numberTex(p.hopping)}\,\mathrm{eV};\quad \tau=${numberTex(p.contact)}\,\mathrm{eV}`,{anchor:'middle',fontSize:18});
    } else {
      s+='<circle cx="340" cy="49" r="13" fill="white"/>';
      if(mode==='fano') {
        s+=(p.coupling>0?line(340,62,340,108):'')+'<circle cx="340" cy="113" r="10" fill="white"/>';
        s+=svgMath(363,119,tex`\varepsilon_s=${numberTex(p.sideEnergy)}\,\mathrm{eV}`,{fontSize:18})+svgMath(360,84,tex`v=${numberTex(p.coupling)}\,\mathrm{eV}`,{fontSize:18});
      } else {
        s+=svgMath(340,91,tex`\varepsilon_0=${numberTex(p.epsilon)}\,\mathrm{eV}`,{anchor:'middle',fontSize:18});
        s+=svgMath(209,32,tex`\Gamma_L=${numberTex(p.gammaL)}\,\mathrm{eV}`,{anchor:'middle',fontSize:18});
        s+=svgMath(477,32,tex`\Gamma_R=${numberTex(p.gammaR)}\,\mathrm{eV}`,{anchor:'middle',fontSize:18});
      }
    }
    $(id).innerHTML=s+'</g>';
  }
  const controlNames = {
    epsilon:'Level energy in electronvolts', gammaL:'Left contact escape width in electronvolts',
    gammaR:'Right contact escape width in electronvolts', probe:'Probe energy in electronvolts',
    sites:'Number of device sites', hopping:'Internal hopping in electronvolts',
    contact:'Contact hopping in electronvolts', onsite:'Device on-site energy in electronvolts',
    sideEnergy:'Side orbital energy in electronvolts', coupling:'Side orbital coupling in electronvolts',
    bias:'Bias in volts', kT:'Thermal energy in electronvolts',
    fraction:'Left contact fraction of the total physical width',
    gamma:'Total physical width in electronvolts', eta:'Numerical broadening in electronvolts',
    limit:'Integration window half-width in electronvolts', intervals:'Uniform grid intervals', alignment:'Grid alignment'
  };
  function slider(group,key,label,min,max,step,unit='') {
    const id=group+'-'+key;
    return '<div class="control"><label for="'+id+'"><span>'+label+'</span><output for="'+id+'" id="'+id+'-value"></output></label><input id="'+id+'" type="range" aria-label="'+controlNames[key]+'" min="'+min+'" max="'+max+'" step="'+step+'" data-group="'+group+'" data-key="'+key+'"><span class="small">'+min+' to '+max+(unit?' '+unit:'')+'</span></div>';
  }
  function select(group,key,label,items) {
    const id=group+'-'+key;
    return '<div class="control"><label for="'+id+'"><span>'+label+'</span></label><select id="'+id+'" aria-label="'+controlNames[key]+'" data-group="'+group+'" data-key="'+key+'">'+items.map(([value,title])=>'<option value="'+value+'">'+title+'</option>').join('')+'</select></div>';
  }
  function sync(group) {
    for (const [key,value] of Object.entries(state[group])) {
      const control=$(group+'-'+key); if(!control)continue;
      control.value=value;
      const output=$(group+'-'+key+'-value'); if(output)output.innerHTML=mathNumber(value);
    }
  }
  function readout(group,paragraphs) { $(group+'-readout').innerHTML=paragraphs.map(s=>'<p>'+s+'</p>').join(''); }
  function record(group,parameters,results,curves) {
    records[group]={format:'quantum-transport-notes/1',version:'1.2.0',experiment:group,
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
      `<strong>At ${math`E=${numberTex(p.probe)}\,\mathrm{eV}`}: ${math`T=${numberTex(r.T)},\ R=${numberTex(r.R)}`}.</strong> ${math`T+R=${numberTex(r.T+r.R,10)}`}.`,
      `Peak ${math`T=${numberTex(peak.T)}`}; full width at half maximum ${math`\Gamma=${numberTex(gamma)}\,\mathrm{eV}`}. Probability lifetime ${math`\hbar/\Gamma=${numberTex(P.constants.hbarEVfs/gamma)}\,\mathrm{fs}`}; amplitude e-folding time ${math`2\hbar/\Gamma=${numberTex(2*P.constants.hbarEVfs/gamma)}\,\mathrm{fs}`}.`,
      `If the probe energy is the Fermi energy, the zero-temperature small-bias conductance is ${math`${numberTex(P.constants.currentMicroampPerEV*r.T)}\,\mu\mathrm{S}`} per spin species.`
    ]); record('level',{...p,eta:0},r,curves);
  }
  function renderChain() {
    const p=state.chain, values=window.TransportSampling.chainCurve({...p,eta:0});
    const curves=[{color:blue,data:values}];
    plot('chain-plot',curves,{xmin:-2.4,xmax:2.4,probe:p.probe,verticals:[-2,2]});
    schematic('chain-schematic','chain',p);
    const r=P.chain(p.probe,p), s=r.selfEnergy;
    readout('chain',[
      `<strong>At ${math`E=${numberTex(p.probe)}\,\mathrm{eV}`}: ${math`T=${numberTex(r.T)}`}.</strong> `+(r.channels?`${math`R=${numberTex(r.R)}`}; ${math`T+R=${numberTex(r.T+r.R,10)}`}.`:'No propagating lead channel at this energy; an incident-flux reflection probability is not defined.'),
      `Each contact: ${math`\Sigma^r=(${numberTex(s.re)}-i\,${numberTex(-s.im)})\,\mathrm{eV}`}; ${math`\Gamma=-2\operatorname{Im}\Sigma^r=${numberTex(s.gamma)}\,\mathrm{eV}`}.`,
      `The finite device has ${p.sites} site${p.sites>1?'s':''}. Leads remain semi-infinite, with hopping ${math`1\,\mathrm{eV}`} and on-site energy ${math`0`}.`
    ]); record('chain',{...p,eta:0,leadHopping:1,leadOnsite:0},r,curves);
  }
  function renderFano() {
    const p={...state.fano,epsilon:0,gammaL:.5,gammaR:.5,eta:0};
    const values=window.TransportSampling.fanoCurve(p);
    const curves=[{color:blue,data:values},{color:gray,dash:true,data:sample(E=>P.resonantLevel(E,p).T,-3,3)}];
    plot('fano-plot',curves,{probe:p.probe}); schematic('fano-schematic','fano',p);
    const r=P.fano(p.probe,p), atSide=P.fano(p.sideEnergy,p);
    readout('fano',[
      `<strong>At ${math`E=${numberTex(p.probe)}\,\mathrm{eV}`}: ${math`T=${numberTex(r.T)},\ R=${numberTex(r.R)}`}.</strong>`,
      p.coupling>0 ? `At the side-orbital energy ${math`E=\varepsilon_s=${numberTex(p.sideEnergy)}\,\mathrm{eV}`}, ${math`T=${numberTex(atSide.T)}`}. The full incident flux is reflected; no probability is absorbed.` : 'The side orbital is disconnected. Its bound state does not enter the transmission. Its occupation requires an initial condition, and its delta-function spectral weight is not drawn.',
      p.coupling>0 ? `Move ${math`\varepsilon_s`} while keeping the conducting orbital fixed: the zero follows ${math`\varepsilon_s`}. Changing ${math`v`} controls the width and hybridized resonances.` : `At ${math`v=0`}, moving ${math`\varepsilon_s`} changes only the isolated bound-state energy. Transmission through the conducting orbital is unchanged.`
    ]); record('fano',p,r,curves);
  }
  function renderBias() {
    const s=state.bias, p={epsilon:s.epsilon,gammaL:.4*s.fraction,gammaR:.4*(1-s.fraction),bias:s.bias,kT:s.kT,eta:0,limit:40,intervals:16000};
    const muL=p.bias/2,muR=-p.bias/2;
    const energies=sample(E=>E,-3,3,900,[muL,muR,p.epsilon]).map(([E])=>E);
    // At zero temperature, the occupation has a discontinuity. A vertical
    // segment displays its one-sided limits without inventing thermal width.
    const occupationCurve = mu => p.kT===0 ? [[-3,1],[mu,1],[mu,0],[3,0]] : energies.map(E=>[E,P.fermi(E,mu,p.kT)]);
    const curves=[
      {color:blue,data:energies.map(E=>[E,P.resonantLevel(E,p).T])},
      {color:red,data:occupationCurve(muL)},
      {color:'#555',dash:true,data:occupationCurve(muR)}
    ];
    plot('bias-plot',curves,{ylabelTex:tex`T\text{ or occupation}`,verticals:[muL,muR]});
    const r=P.bias(p), fine=P.bias({...p,intervals:32000});
    const half=.2, exactArea=(Math.atan((40-p.epsilon)/half)-Math.atan((-40-p.epsilon)/half))/Math.PI;
    readout('bias',[
      `<strong>${math`I=${numberTex(fine.I,5)}\,\mu\mathrm{A}`}; ${math`n=${numberTex(fine.n,5)}`} electron per spin orbital.</strong> ${math`\mu_L=${numberTex(muL)}\,\mathrm{eV}`}; ${math`\mu_R=${numberTex(muR)}\,\mathrm{eV}`}.`,
      `${math`\Gamma_L=${numberTex(p.gammaL)}\,\mathrm{eV}`}, ${math`\Gamma_R=${numberTex(p.gammaR)}\,\mathrm{eV}`}. ${math`I_L+I_R=${numberTex(fine.IL+fine.IR)}\,\mu\mathrm{A}`}. Temperature = ${math`${numberTex(p.kT/P.constants.kBoltzmannEV,1)}\,\mathrm{K}`}.`,
      `Integration: ${math`-40`} to ${math`+40\,\mathrm{eV}`}, 32,000 base intervals (${math`\Delta E=0.0025\,\mathrm{eV}`}); at zero temperature integration also splits at the chemical potentials. Halving ${math`\Delta E`} changes ${math`I`} by ${math`${numberTex(Math.abs(fine.I-r.I),7)}\,\mu\mathrm{A}`} and ${math`n`} by ${math`${numberTex(Math.abs(fine.n-r.n),7)}`}. The finite window omits ${mathNumber(1-exactArea,5)} of the total spectral weight; this bounds the omitted occupation, independently of grid error. The stated ${math`n`} is the finite-window result.`+(p.epsilon===0 && p.bias===0?` The exact infinite-window equilibrium occupation is ${math`n=1/2`} by symmetry.`:'')
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
      readout('eta',[`<strong>An isolated level at ${math`\eta=0`} has a delta-function spectrum.</strong> Choose a positive ${math`\eta`} to sample it. No ordinary finite curve or uniform-grid integral is defined at its pole.`]);
      $('eta-plot').innerHTML=svgMath(340,150,tex`\delta(E):\quad\text{choose }\eta>0`,{anchor:'middle',fontSize:24})+svgMath(340,185,tex`\text{to draw a regularized spectrum.}`,{anchor:'middle',fontSize:18});
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
    plot('eta-plot',curves,{xmin:-.15,xmax:.15,ymax:d.peak.A/(2*Math.PI)*1.06,ylabelTex:tex`\text{Spectral density }(\mathrm{eV}^{-1})`,density:true,dots});
    const currents=exactSinkCurrents(s.gamma,s.eta), ratio=s.gamma/d.width;
    readout('eta',[
      `<strong>Sampled full area = ${mathNumber(d.area,7)}; exact same-window area = ${mathNumber(d.analyticArea,7)}.</strong>`,
      `Quadrature error = ${mathNumber(d.quadratureError,7)}; omitted spectral tail = ${mathNumber(d.omittedTail,7)}. ${math`\Delta E=${numberTex(d.gridStep,6)}\,\mathrm{eV}`}; full width ${math`\Gamma+2\eta=${numberTex(d.width,6)}\,\mathrm{eV}`}. `+(d.gridStep>half?'The integration spacing exceeds the half-width; inspect the alignment and refinement table.':'Refine and shift the grid to test resolution.'),
      `At resonance for left incidence: ${math`T=${numberTex(d.peakT)}`}, ${math`R=${numberTex(d.peakR)}`}, absorbed probability = ${mathNumber(d.peakLoss)}. Their sum is ${mathNumber(d.peakT+d.peakR+d.peakLoss,10)}. The contact-only spectral area on the infinite axis is ${mathNumber(ratio)}.`,
      `At ${math`V=1\,\mathrm{V}`} and temperature zero (analytic integrals): ${math`I_L=${numberTex(currents.IL)}\,\mu\mathrm{A}`}; ${math`I_R=${numberTex(currents.IR)}\,\mu\mathrm{A}`}; ${math`I_L+I_R=I_{\mathrm{sink}}=${numberTex(currents.Isink)}\,\mu\mathrm{A}`}. Both lead currents are defined entering the device. `+(s.gamma===0?'There is no physical-lead injection.':`The pair current ${math`I=${numberTex(currents.I)}\,\mu\mathrm{A}`} alone hides this sink.`),
      s.eta>0 ? `If ${math`\eta`} is deliberately interpreted as absorption, its additional probability e-folding time is ${math`\hbar/(2\eta)=${numberTex(P.constants.hbarEVfs/(2*s.eta))}\,\mathrm{fs}`}. This is not a lifetime measurement.` : `${math`\eta=0`} adds no artificial absorption; any linewidth here comes from physical contacts.`
    ]);
    const table=[];
    for(const eta of [...new Set([s.eta,s.eta/10,0])]) for(const multiplier of [1,2,4]) {
      if(s.gamma+2*eta===0)continue;
      const n=s.intervals*multiplier;
      const result=P.etaDiagnostics({...p,eta,intervals:n,offset:s.alignment*2*s.limit/n});
      table.push({eta,intervals:n,...result});
    }
    $('eta-table').innerHTML=`<table><caption>Recomputed convergence, same fractional grid alignment</caption><thead><tr><th>${math`\eta\;(\mathrm{eV})`}</th><th>Intervals</th><th>${math`\Delta E\;(\mathrm{eV})`}</th><th>Sampled area</th><th>Exact window</th><th>Grid error</th></tr></thead><tbody>`+table.map(r=>'<tr><td>'+mathNumber(r.eta)+'</td><td>'+r.intervals+'</td><td>'+mathNumber(r.gridStep,6)+'</td><td>'+mathNumber(r.area,6)+'</td><td>'+mathNumber(r.analyticArea,6)+'</td><td>'+mathNumber(r.quadratureError,6)+'</td></tr>').join('')+'</tbody></table>';
    record('eta',p,{...d,currents,contactOnlyInfiniteArea:ratio,convergence:table,integrationSamplesInPlot:dots},curves);
  }
  const renderers={level:renderLevel,chain:renderChain,fano:renderFano,bias:renderBias,eta:renderEta};
  function render(group) {
    // Exports describe one parameter snapshot. Hide it when those parameters
    // change so a previous download cannot be mistaken for the current plot.
    const source=$(group+'-export-source'), anchor=$(group+'-download');
    if(source) {source.hidden=true;source.open=false;source.querySelector('pre').textContent='';}
    if(anchor) {anchor.hidden=true;anchor.removeAttribute('href');}
    for(const kind of ['json','svg']) {
      const key=group+kind; if(downloads[key]) {URL.revokeObjectURL(downloads[key]);delete downloads[key];}
    }
    try {sync(group);renderers[group]();delete errors[group];}
    catch(error) {delete records[group];errors[group]=error.message;console.error(error);}
    $('app-error').textContent=Object.entries(errors).map(([name,message])=>'Calculation '+name+' could not finish: '+message).join(' ');
    $('app-error').hidden=Object.keys(errors).length===0;
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
      const legends={level:[tex`\text{Transmission }T`,tex`\text{Reflection }R`],chain:[tex`\text{Transmission}`],
        fano:[tex`\text{With side orbital}`,tex`\text{Without side orbital}`],bias:[tex`\text{Transmission }T`,tex`f_L`,tex`f_R`],
        eta:[tex`\text{Full }A/(2\pi)`,tex`(A_L+A_R)/(2\pi)`]};
      const p=r.parameters;
      const captions={
        level:tex`\varepsilon_0=${numberTex(p.epsilon)}\,\mathrm{eV};\quad\Gamma_L=${numberTex(p.gammaL)}\,\mathrm{eV};\quad\Gamma_R=${numberTex(p.gammaR)}\,\mathrm{eV};\quad\eta=0`,
        chain:tex`N=${p.sites};\quad t=${numberTex(p.hopping)}\,\mathrm{eV};\quad\tau=${numberTex(p.contact)}\,\mathrm{eV};\quad\varepsilon_d=${numberTex(p.onsite)}\,\mathrm{eV};\quad\eta=0`,
        fano:tex`\varepsilon_s=${numberTex(p.sideEnergy)}\,\mathrm{eV};\quad v=${numberTex(p.coupling)}\,\mathrm{eV};\quad\Gamma_L=\Gamma_R=0.5\,\mathrm{eV};\quad\eta=0`,
        bias:tex`\varepsilon_0=${numberTex(p.epsilon)}\,\mathrm{eV};\quad V=${numberTex(p.bias)}\,\mathrm{V};\quad k_BT=${numberTex(p.kT)}\,\mathrm{eV};\quad\Gamma_L/\Gamma=${numberTex(p.gammaL/.4)}`,
        eta:tex`\Gamma=${numberTex(p.gammaL+p.gammaR)}\,\mathrm{eV};\quad\eta=${numberTex(p.eta)}\,\mathrm{eV};\quad\text{window half-width }${numberTex(p.limit)}\,\mathrm{eV}`
      };
      const footer=document.createElementNS(ns,'g');footer.setAttribute('font-family','Times New Roman,serif');footer.setAttribute('font-size','15');footer.setAttribute('fill','#111');
      let s='<rect x="0" y="340" width="680" height="100" fill="white"/>';
      (legends[group]||[]).forEach((label,i)=>{
        const x=34+i*215,curve=r.curves[i]; if(!curve)return;
        s+=line(x,361,x+25,361,'stroke="'+curve.color+'" stroke-width="2"'+(curve.dash?' stroke-dasharray="5 3"':''));
        s+=svgMath(x+31,366,label,{fontSize:15});
      });
      if(group==='eta' && r.curves.length)s+='<circle cx="477" cy="361" r="2.5"/>'+text(489,366,'Integration samples');
      s+=svgMath(34,397,captions[group],{fontSize:14});
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
    await window.TransportSvgMath.ready;
    if(!P)throw new Error('physics.js did not load');
    $('level-controls').innerHTML='<h3>Change the orbital</h3>'+slider('level','epsilon','Level '+math`\varepsilon_0`,-1.5,1.5,.05,'eV')+slider('level','gammaL','Left width '+math`\Gamma_L`,.05,1,.05,'eV')+slider('level','gammaR','Right width '+math`\Gamma_R`,.05,1,.05,'eV')+slider('level','probe','Probe energy',-2.5,2.5,.05,'eV');
    $('chain-controls').innerHTML='<h3>Change the chain</h3>'+slider('chain','sites','Sites '+math`N`,1,9,1)+slider('chain','hopping','Internal hopping '+math`t`,.4,1.4,.05,'eV')+slider('chain','contact','Contact hopping '+math`\tau`,.15,1.4,.05,'eV')+slider('chain','onsite','Device on-site '+math`\varepsilon_d`,-1.5,1.5,.1,'eV')+slider('chain','probe','Probe energy',-2.4,2.4,.05,'eV');
    $('fano-controls').innerHTML='<h3>Change the detour</h3>'+slider('fano','sideEnergy','Side energy '+math`\varepsilon_s`,-1.5,1.5,.05,'eV')+slider('fano','coupling','Coupling '+math`v`,0,1,.025,'eV')+slider('fano','probe','Probe energy',-2.5,2.5,.05,'eV');
    $('bias-controls').innerHTML='<h3>Fill the reservoirs</h3>'+slider('bias','epsilon','Level '+math`\varepsilon_0`,-2,2,.05,'eV')+slider('bias','bias','Bias '+math`V`,-4,4,.1,'V')+slider('bias','kT','Thermal energy '+math`k_BT`,0,.2,.005,'eV')+slider('bias','fraction',math`\Gamma_L/(\Gamma_L+\Gamma_R)`,.05,.95,.05)+`<p class="small">Total physical width ${math`\Gamma_L+\Gamma_R=0.4\,\mathrm{eV}`}. Numerical ${math`\eta=0`}.</p>`;
    $('eta-controls').innerHTML='<h3>Audit the numerics</h3>'+
      select('eta','gamma','Total physical width '+math`\Gamma\;(\mathrm{eV})`,[[0,'0 (isolated)'],[.002,'0.002'],[.02,'0.02'],[.2,'0.2'],[.6,'0.6']])+
      select('eta','eta','Numerical '+math`\eta\;(\mathrm{eV})`,[[0,'0'],[.00001,'0.00001'],[.0001,'0.0001'],[.001,'0.001'],[.01,'0.01'],[.05,'0.05'],[.1,'0.1']])+
      select('eta','limit','Window half-width (eV)',[[1,'1'],[4,'4'],[8,'8'],[20,'20']])+
      select('eta','intervals','Uniform grid intervals',[[200,'200'],[400,'400'],[800,'800'],[1600,'1,600'],[6400,'6,400'],[25600,'25,600']])+
      select('eta','alignment','Grid alignment',[[0,'Peak on a grid point'],[.5,'Shift by half a grid step']])+
      `<p class="small">${math`\varepsilon_0=0`}; ${math`\Gamma_L=\Gamma_R=\Gamma/2`}. Shifting the grid also shifts the two window endpoints; the exact area uses those actual endpoints.</p>`;
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
