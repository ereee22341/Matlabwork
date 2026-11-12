%% F2-b 主入口：重新计算全部结果并生成图表。
clear; clc; close all;
root=fileparts(mfilename('fullpath')); addpath(root);
cfg=config_f2b();
if ~exist(cfg.outdir,'dir'), mkdir(cfg.outdir); end
timer=tic;
metrics=validate_baseline(cfg);
base=load(fullfile(cfg.outdir,'baseline.mat'));
num=base.num; theory=base.theory;
%% 膜厚偏差分析。
scales=[.9,1,1.1]; thickness=cell(1,3); thickTheory=thickness;
for k=1:3
    thickness{k}=fdtd_spectrum(cfg.nfilm,cfg.d0*scales(k),cfg.dz,cfg);
    thickTheory{k}=tmm_spectrum(cfg.f,cfg.nfilm,cfg.d0*scales(k),cfg);
end
%% 固定物理膜厚，用整数网格划分研究收敛性。
cells=[10,20,40]; grids=cell(1,3); errors=zeros(3,3);
for k=1:3
    grids{k}=fdtd_spectrum(cfg.nfilm,cfg.d0,cfg.d0/cells(k),cfg);
    valid=grids{k}.valid;
    errors(k,:)=[sqrt(mean((grids{k}.R(valid)-theory.R(valid)).^2)),...
        sqrt(mean((grids{k}.T(valid)-theory.T(valid)).^2)),...
        max(abs(grids{k}.R(valid)+grids{k}.T(valid)-1))];
end
assert(all(all(diff(errors(:,1:2),1,1)<0)),'网格加密后，均方根误差未单调减小。');
gridTable=table(cells',cfg.d0./cells'*1e3,errors(:,1),errors(:,2),errors(:,3),...
    'VariableNames',{'FilmCells','dz_mm','R_RMS','T_RMS','EnergyMax'});
writetable(gridTable,fullfile(cfg.outdir,'grid_errors.csv'));
% 物理几何保持不变，因此可用相邻网格的误差比估计收敛阶。
orderR=log2(errors(1:end-1,1)./errors(2:end,1));
orderT=log2(errors(1:end-1,2)./errors(2:end,2));
%% 膜厚取整：分别与设计厚度及实际厚度的 TMM 结果比较。
roundDz=cfg.c0/(2*cfg.f0*cfg.ns*20);
roundD=round(cfg.d0/roundDz)*roundDz;
rounded=fdtd_spectrum(cfg.nfilm,roundD,roundDz,cfg);
roundTheory=tmm_spectrum(cfg.f,cfg.nfilm,roundD,cfg);
%% 必做扩展：双层增透膜，层顺序从空气侧指向基底侧。
dualN=[3^(1/4),3^(3/4)]; dualD=cfg.c0./(4*cfg.f0*dualN);
dualTheory=tmm_spectrum(cfg.f,dualN,dualD,cfg);
dual=fdtd_spectrum(dualN,dualD,cfg.d0/40,cfg);
fineF=linspace(0,2*cfg.f0,12001);
fineSingle=tmm_spectrum(fineF,cfg.nfilm,cfg.d0,cfg);
fineDual=tmm_spectrum(fineF,dualN,dualD,cfg);
bands=[reflection_band(fineF,fineSingle.R,cfg.f0,.01);...
    reflection_band(fineF,fineDual.R,cfg.f0,.01);...
    reflection_band(num.f,num.R,cfg.f0,.01);...
    reflection_band(dual.f,dual.R,cfg.f0,.01)];
bandTable=table({'Single TMM';'Double TMM';'Single FDTD';'Double FDTD'},...
    bands(:,1)/1e9,bands(:,2)/1e9,bands(:,3)*100,...
    'VariableNames',{'Method','Lower_GHz','Upper_GHz','FractionalBW_percent'});
writetable(bandTable,fullfile(cfg.outdir,'bandwidths.csv'));
%% 仿真时长敏感性：记录窗口仍在主要边界往返波到达之前。
longCfg=cfg; longCfg.tend=2.15e-9;
longRun=fdtd_spectrum(cfg.nfilm,cfg.d0,cfg.dz,longCfg);
metrics.timeSensitivityR=max(abs(longRun.R-num.R));
metrics.timeSensitivityT=max(abs(longRun.T-num.T));
metrics.dualRrms=sqrt(mean((dual.R(dual.valid)-dualTheory.R(dual.valid)).^2));
assert(metrics.timeSensitivityR<.005 && metrics.timeSensitivityT<.005,'改变仿真时长后，频谱变化超过允许范围。');
assert(metrics.dualRrms<.01,'双层膜反射谱误差超过允许范围。');
% 在最低点附近进行三点抛物线插值，估计反射谷的频率。
minima=zeros(3,5);
for k=1:3
    [~,ix]=min(thickness{k}.R);
    xx=cfg.f(ix-1:ix+1)/1e9; yy=thickness{k}.R(ix-1:ix+1);
    p=polyfit(xx-xx(2),yy,2); fmin=xx(2)-p(2)/(2*p(1));
    minima(k,:)=[scales(k),cfg.d0*scales(k)*1e3,cfg.f0/1e9/scales(k),fmin,...
        sqrt(mean((thickness{k}.R(thickness{k}.valid)-thickTheory{k}.R(thickness{k}.valid)).^2))];
end
assert(all(abs(minima(:,4)-minima(:,3))<.01),'反射谷频率与理论值的偏差超过 10 MHz。');
thicknessTable=array2table(minima,'VariableNames',...
    {'ThicknessScale','d_mm','TheoryMinimum_GHz','FDTDMinimum_GHz','R_RMS'});
writetable(thicknessTable,fullfile(cfg.outdir,'thickness_minima.csv'));
make_figures(cfg,num,theory,thickness,thickTheory,grids,gridTable,...
    rounded,roundTheory,dual,dualTheory,longRun);
elapsed=toc(timer);
save(fullfile(cfg.outdir,'results.mat'),'cfg','metrics','num','theory',...
    'thickness','thickTheory','grids','gridTable','rounded','roundTheory',...
    'roundD','roundDz','dual','dualN','dualD','dualTheory','bandTable','elapsed',...
    'orderR','orderT','thicknessTable');
spectrumTable=table(cfg.f'/1e9,num.R',theory.R',num.t2',theory.t2',num.T',theory.T',...
    'VariableNames',{'f_GHz','R_FDTD','R_TMM','t2_FDTD','t2_TMM','T_FDTD','T_TMM'});
writetable(spectrumTable,fullfile(cfg.outdir,'baseline_spectra.csv'));
names=fieldnames(metrics); values=struct2cell(metrics);
writetable(table(names,cell2mat(values),'VariableNames',{'Metric','Value'}),...
    fullfile(cfg.outdir,'validation_metrics.csv'));
print_summary(cfg,gridTable,bandTable,thicknessTable,metrics,orderR,orderT,elapsed);
assert(elapsed<300,'完整运行时间超过任务书规定的 5 分钟。');

