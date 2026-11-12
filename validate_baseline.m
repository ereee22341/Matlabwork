function metrics = validate_baseline(cfg)
% Physical checks, independent closed-form spectrum and reference TMM.
bare=tmm_spectrum(cfg.f,[],[],cfg);
theory=tmm_spectrum(cfg.f,cfg.nfilm,cfg.d0,cfg);
closed=cos(pi*cfg.f/(2*cfg.f0)).^2;
closed=closed./(3+closed);
assert(max(abs(theory.R-closed))<1e-12,'TMM closed-form mismatch');
assert(max(abs(bare.R-0.25))<1e-12,'Bare-interface mismatch');
assert(max(abs(theory.R+theory.T-1))<1e-12,'TMM energy mismatch');
num=fdtd_spectrum(cfg.nfilm,cfg.d0,cfg.dz,cfg);
bareNum=fdtd_spectrum([],[],cfg.dz,cfg);
mask=num.valid;
metrics.Rrms=sqrt(mean((num.R(mask)-theory.R(mask)).^2));
metrics.Trms=sqrt(mean((num.T(mask)-theory.T(mask)).^2));
metrics.energyMax=max(abs(num.R(mask)+num.T(mask)-1));
metrics.bareRrms=sqrt(mean((bareNum.R(mask)-bare.R(mask)).^2));
metrics.sourceMin=min(num.sourceLevel);
assert(metrics.Rrms<0.01,'FDTD reflection RMS exceeds 0.01');
assert(metrics.Trms<0.02,'FDTD transmission RMS exceeds 0.02');
assert(metrics.energyMax<0.03,'FDTD energy residual exceeds 0.03');
fprintf('R RMS %.4g; T RMS %.4g; energy max %.4g; bare R RMS %.4g\n',...
    metrics.Rrms,metrics.Trms,metrics.energyMax,metrics.bareRrms);
save(fullfile(cfg.outdir,'baseline.mat'),'cfg','num','theory','metrics','bareNum');
end
