%% F2-b: run this entry point to regenerate all analyses and figures.
clear; clc; close all;
root=fileparts(mfilename('fullpath')); addpath(root);
cfg=config_f2b();
if ~exist(cfg.outdir,'dir'), mkdir(cfg.outdir); end
timer=tic;
metrics=validate_baseline(cfg);
base=load(fullfile(cfg.outdir,'baseline.mat'));
num=base.num; theory=base.theory;
%% Thickness sensitivity.
scales=[.9,1,1.1]; thickness=cell(1,3); thickTheory=thickness;
for k=1:3
    thickness{k}=fdtd_spectrum(cfg.nfilm,cfg.d0*scales(k),cfg.dz,cfg);
    thickTheory{k}=tmm_spectrum(cfg.f,cfg.nfilm,cfg.d0*scales(k),cfg);
end
%% Convergence at fixed physical thickness, integer film-cell counts.
cells=[10,20,40]; grids=cell(1,3); errors=zeros(3,3);
for k=1:3
    grids{k}=fdtd_spectrum(cfg.nfilm,cfg.d0,cfg.d0/cells(k),cfg);
    valid=grids{k}.valid;
    errors(k,:)=[sqrt(mean((grids{k}.R(valid)-theory.R(valid)).^2)),...
        sqrt(mean((grids{k}.T(valid)-theory.T(valid)).^2)),...
        max(abs(grids{k}.R(valid)+grids{k}.T(valid)-1))];
end
assert(all(errors(end,1:2)<errors(1,1:2)),'Refinement failed to improve RMS');
gridTable=table(cells',cfg.d0./cells'*1e3,errors(:,1),errors(:,2),errors(:,3),...
    'VariableNames',{'FilmCells','dz_mm','R_RMS','T_RMS','EnergyMax'});
writetable(gridTable,fullfile(cfg.outdir,'grid_errors.csv'));
% Observed orders are meaningful here because physical geometry stays fixed.
orderR=log2(errors(1:end-1,1)./errors(2:end,1));
orderT=log2(errors(1:end-1,2)./errors(2:end,2));
%% Rounding: compare both intended and actual geometry against TMM.
roundDz=cfg.c0/(2*cfg.f0*cfg.ns*20);
roundD=round(cfg.d0/roundDz)*roundDz;
rounded=fdtd_spectrum(cfg.nfilm,roundD,roundDz,cfg);
roundTheory=tmm_spectrum(cfg.f,cfg.nfilm,roundD,cfg);
%% Required double-layer extension: order is air -> substrate.
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
%% Time truncation sensitivity, still before principal Mur round trips.
longCfg=cfg; longCfg.tend=2.15e-9;
longRun=fdtd_spectrum(cfg.nfilm,cfg.d0,cfg.dz,longCfg);
metrics.timeSensitivityR=max(abs(longRun.R-num.R));
metrics.timeSensitivityT=max(abs(longRun.T-num.T));
metrics.dualRrms=sqrt(mean((dual.R(dual.valid)-dualTheory.R(dual.valid)).^2));
assert(metrics.timeSensitivityR<.005 && metrics.timeSensitivityT<.005,'Time sensitivity too large');
assert(metrics.dualRrms<.01,'Double-layer error too large');
% Tabulate minima using parabolic interpolation around the sampled minimum.
minima=zeros(3,5);
for k=1:3
    [~,ix]=min(thickness{k}.R);
    xx=cfg.f(ix-1:ix+1)/1e9; yy=thickness{k}.R(ix-1:ix+1);
    p=polyfit(xx-xx(2),yy,2); fmin=xx(2)-p(2)/(2*p(1));
    minima(k,:)=[scales(k),cfg.d0*scales(k)*1e3,6/scales(k),fmin,...
        sqrt(mean((thickness{k}.R(thickness{k}.valid)-thickTheory{k}.R(thickness{k}.valid)).^2))];
end
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
disp(gridTable); disp(bandTable); disp(thicknessTable); disp(metrics);
fprintf('Reflectance convergence orders: %.3f, %.3f\n',orderR);
fprintf('Completed in %.2f s. Outputs: %s\n',elapsed,cfg.outdir);
assert(elapsed<300,'Runtime exceeds 5 minutes');

