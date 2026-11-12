function s = fdtd_spectrum(indices, thicknesses, dz, cfg)
% Identical source/monitor geometry in reference and structure simulations.
% s = fdtd_spectrum() computes the default F2-b single-layer spectrum only.
% Run main.m to generate all analyses, figures and output files.
if nargin==0
    cfg=config_f2b();
    indices=cfg.nfilm; thicknesses=cfg.d0; dz=cfg.dz;
    fprintf(['Using default F2-b single-layer parameters.\n',...
        'For all analyses and figures, run main.m.\n']);
elseif nargin~=4
    error('F2B:InputCount', ['fdtd_spectrum requires either no inputs or four inputs. ',...
        'Example: s = fdtd_spectrum(cfg.nfilm,cfg.d0,cfg.dz,cfg). ',...
        'Run main.m for the complete project.']);
end
validateattributes(dz,{'numeric'},{'real','finite','scalar','positive'},mfilename,'dz');
validateattributes(indices,{'numeric'},{'real','finite','positive'},mfilename,'indices');
validateattributes(thicknesses,{'numeric'},{'real','finite','positive'},mfilename,'thicknesses');
if numel(indices)~=numel(thicknesses)
    error('F2B:LayerCount','Each layer must have one refractive index and one thickness.');
end
indices=indices(:).'; thicknesses=thicknesses(:).';
ref=fdtd_run(indices,thicknesses,dz,cfg,true);
obj=fdtd_run(indices,thicknesses,dz,cfg,false);
% Dense FFT sampling prevents interpolation from shifting the reflection notch.
Nfft=2^nextpow2(numel(obj.time)*32);
freq=(0:Nfft/2)/(Nfft*obj.dt);
inc=fft(ref.er,Nfft); incT=fft(ref.et,Nfft);
refl=fft(obj.er-ref.er,Nfft); trans=fft(obj.et,Nfft);
sel=find(freq<=max(cfg.f)+2/(Nfft*obj.dt) & ...
    abs(inc(1:numel(freq)))>1e-8*max(abs(inc)) & ...
    abs(incT(1:numel(freq)))>1e-8*max(abs(incT)));
freq=freq(sel);
% Interpolate power ratios rather than complex phases at remote monitors.
rp=abs(refl(sel)./inc(sel)).^2;
tp=abs(trans(sel)./incT(sel)).^2;
s.R=interp1(freq,rp,cfg.f,'pchip');
s.t2=interp1(freq,tp,cfg.f,'pchip');
s.T=cfg.ns/cfg.n0*s.t2;
s.sourceLevel=interp1(freq,abs(inc(sel))/max(abs(inc(sel))),cfg.f,'linear');
s.valid=s.sourceLevel>1e-3 & cfg.f>0 & isfinite(s.R) & isfinite(s.T);
s.f=cfg.f; s.reference=ref; s.run=obj;
end
