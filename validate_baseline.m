function metrics = validate_baseline(cfg)
% 物理校验：独立闭式解、无膜界面及 TMM 参考谱。
bare=tmm_spectrum(cfg.f,[],[],cfg);
theory=tmm_spectrum(cfg.f,cfg.nfilm,cfg.d0,cfg);
closed=cos(pi*cfg.f/(2*cfg.f0)).^2;
closed=closed./(3+closed);
assert(max(abs(theory.R-closed))<1e-12,'TMM 与单层膜闭式解不一致。');
assert(max(abs(bare.R-0.25))<1e-12,'无膜界面反射率与理论值 0.25 不一致。');
assert(max(abs(theory.R+theory.T-1))<1e-12,'TMM 结果未满足能量守恒。');
num=fdtd_spectrum(cfg.nfilm,cfg.d0,cfg.dz,cfg);
bareNum=fdtd_spectrum([],[],cfg.dz,cfg);
mask=num.valid & bareNum.valid;
assert(all(mask(cfg.f>0)),'入射频谱未完整覆盖所需的正频率范围。');
metrics.Rrms=sqrt(mean((num.R(mask)-theory.R(mask)).^2));
metrics.Trms=sqrt(mean((num.T(mask)-theory.T(mask)).^2));
metrics.energyMax=max(abs(num.R(mask)+num.T(mask)-1));
metrics.bareRrms=sqrt(mean((bareNum.R(mask)-bare.R(mask)).^2));
metrics.sourceMin=min(num.sourceLevel);
assert(metrics.Rrms<0.01,'FDTD 反射率均方根误差超过 0.01。');
assert(metrics.Trms<0.02,'FDTD 透射率均方根误差超过 0.02。');
assert(metrics.energyMax<0.03,'FDTD 能量守恒最大偏差超过 0.03。');
assert(metrics.bareRrms<0.01,'无膜界面 FDTD 反射率均方根误差超过 0.01。');
fprintf('基础检查：反射率 RMS %.4g；透射率 RMS %.4g；最大能量偏差 %.4g；无膜反射率 RMS %.4g。\n',...
    metrics.Rrms,metrics.Trms,metrics.energyMax,metrics.bareRrms);
save(fullfile(cfg.outdir,'baseline.mat'),'cfg','num','theory','metrics','bareNum');
end
