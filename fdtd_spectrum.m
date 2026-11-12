function s = fdtd_spectrum(indices, thicknesses, dz, cfg)
% 参考仿真与结构仿真使用完全相同的源和观察点布局。
% s = fdtd_spectrum() 仅计算默认 F2-b 单层膜频谱。
% 生成全部分析、图表和结果文件，请运行 main.m。
if nargin==0
    cfg=config_f2b();
    indices=cfg.nfilm; thicknesses=cfg.d0; dz=cfg.dz;
    fprintf(['已使用默认 F2-b 参数计算单层膜频谱。\n',...
        '如需全部分析和图表，请运行 main.m。\n']);
elseif nargin~=4
    error('F2B:InputCount', ['fdtd_spectrum 需要零个或四个输入参数。',...
        '示例：s = fdtd_spectrum(cfg.nfilm,cfg.d0,cfg.dz,cfg)。',...
        '运行完整项目请使用 main.m。']);
end
validateattributes(dz,{'numeric'},{'real','finite','scalar','positive'},mfilename,'dz');
validateattributes(indices,{'numeric'},{'real','finite','positive'},mfilename,'indices');
validateattributes(thicknesses,{'numeric'},{'real','finite','positive'},mfilename,'thicknesses');
if numel(indices)~=numel(thicknesses)
    error('F2B:LayerCount','每一层都必须同时指定折射率与厚度，两个数组的元素数量应一致。');
end
indices=indices(:).'; thicknesses=thicknesses(:).';
ref=fdtd_run(indices,thicknesses,dz,cfg,true);
obj=fdtd_run(indices,thicknesses,dz,cfg,false);
% 采用较密的 FFT 采样，避免插值导致反射谷位置偏移。
Nfft=2^nextpow2(numel(obj.time)*32);
freq=(0:Nfft/2)/(Nfft*obj.dt);
inc=fft(ref.er,Nfft); incT=fft(ref.et,Nfft);
refl=fft(obj.er-ref.er,Nfft); trans=fft(obj.et,Nfft);
sel=find(freq<=max(cfg.f)+2/(Nfft*obj.dt) & ...
    abs(inc(1:numel(freq)))>1e-8*max(abs(inc)) & ...
    abs(incT(1:numel(freq)))>1e-8*max(abs(incT)));
freq=freq(sel);
% 对功率比例插值，避免远处观察点的快速相位变化影响插值。
rp=abs(refl(sel)./inc(sel)).^2;
tp=abs(trans(sel)./incT(sel)).^2;
s.R=interp1(freq,rp,cfg.f,'pchip');
s.t2=interp1(freq,tp,cfg.f,'pchip');
s.T=cfg.ns/cfg.n0*s.t2;
s.sourceLevel=interp1(freq,abs(inc(sel))/max(abs(inc(sel))),cfg.f,'linear');
s.valid=s.sourceLevel>1e-3 & cfg.f>0 & isfinite(s.R) & isfinite(s.T);
s.f=cfg.f; s.reference=ref; s.run=obj;
end
