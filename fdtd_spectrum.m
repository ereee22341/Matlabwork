function s = fdtd_spectrum(indices, thicknesses, dz, cfg)
% Identical source/monitor geometry in reference and structure simulations.
ref=fdtd_run(indices,thicknesses,dz,cfg,true);
obj=fdtd_run(indices,thicknesses,dz,cfg,false);
Nfft=2^nextpow2(numel(obj.time)*4);
freq=(0:Nfft/2)/(Nfft*obj.dt);
inc=fft(ref.er,Nfft); incT=fft(ref.et,Nfft);
refl=fft(obj.er-ref.er,Nfft); trans=fft(obj.et,Nfft);
sel=1:Nfft/2+1;
% Interpolate power ratios rather than complex phases at remote monitors.
rp=abs(refl(sel)./inc(sel)).^2;
tp=abs(trans(sel)./incT(sel)).^2;
s.R=interp1(freq,rp,cfg.f,'pchip');
s.t2=interp1(freq,tp,cfg.f,'pchip');
s.T=cfg.ns/cfg.n0*s.t2;
s.sourceLevel=interp1(freq,abs(inc(sel))/max(abs(inc(sel))),cfg.f,'linear');
s.valid=s.sourceLevel>1e-3 & cfg.f>0;
s.f=cfg.f; s.reference=ref; s.run=obj;
end
