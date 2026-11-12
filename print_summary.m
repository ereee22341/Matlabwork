function print_summary(cfg,gridTable,bandTable,thicknessTable,metrics,orderR,orderT,elapsed)
% 同时输出中文命令行摘要与 UTF-8 文本，保留数值表的原有字段名。
summaryPath=fullfile(cfg.outdir,'results_summary.txt');
fid=fopen(summaryPath,'w','n','UTF-8');
assert(fid~=-1,'无法创建中文结果汇总文件。');
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
say('\nF2-b 计算结果汇总\n========================================\n');
say('定义：R=|r|^2 为功率反射率，T=3|t|^2 为功率透射率。\n');
say('均方根误差为功率比例的绝对误差，未采用零反射点的相对误差。\n\n');
say('网格收敛分析\n');
say('膜厚内网格数 | 网格间距(mm) | 反射率均方根误差 | 透射率均方根误差 | 最大能量偏差\n');
for k=1:height(gridTable)
    say('%12d | %12.6f | %18.6g | %18.6g | %12.6g\n',...
        gridTable.FilmCells(k),gridTable.dz_mm(k),gridTable.R_RMS(k),...
        gridTable.T_RMS(k),gridTable.EnergyMax(k));
end
say('反射率收敛阶：%.3f、%.3f；透射率收敛阶：%.3f、%.3f。\n\n',...
    orderR(1),orderR(2),orderT(1),orderT(2));
say('反射率小于 1%% 的连续频带\n');
labels={'单层膜 TMM','双层膜 TMM','单层膜 FDTD','双层膜 FDTD'};
for k=1:height(bandTable)
    say('%s：%.6f～%.6f GHz，相对带宽 %.3f%%。\n',labels{k},...
        bandTable.Lower_GHz(k),bandTable.Upper_GHz(k),bandTable.FractionalBW_percent(k));
end
say('\n膜厚偏差与反射谷频率\n');
for k=1:height(thicknessTable)
    say('膜厚 %.1f d0（%.6f mm）：理论 %.6f GHz，FDTD %.6f GHz。\n',...
        thicknessTable.ThicknessScale(k),thicknessTable.d_mm(k),...
        thicknessTable.TheoryMinimum_GHz(k),thicknessTable.FDTDMinimum_GHz(k));
end
say('\n基础验证\n');
say('反射率均方根误差：%.6g；透射率均方根误差：%.6g。\n',metrics.Rrms,metrics.Trms);
say('能量守恒最大偏差：%.4f%%；无膜界面反射率均方根误差：%.6g。\n',...
    100*metrics.energyMax,metrics.bareRrms);
say('入射频谱最低归一化幅度：%.6g。\n',metrics.sourceMin);
say('时长变化引起的最大差异：反射率 %.6g，透射率 %.6g。\n',...
    metrics.timeSensitivityR,metrics.timeSensitivityT);
say('双层膜反射率均方根误差：%.6g。\n',metrics.dualRrms);
say('\n计算与绘图用时：%.2f 秒。输出目录：%s\n',elapsed,cfg.outdir);
if cfg.showFigures
    say('已显示并保留 13 个图窗，同时保存 PNG 图片和可编辑 FIG 文件。\n');
else
    say('当前为批量运行模式，图表已保存到 output 文件夹，可打开 PNG 或 FIG 查看。\n');
end
say('中文汇总文件：%s\n',summaryPath);
    function say(varargin)
        fprintf(1,varargin{:});
        fprintf(fid,varargin{:});
    end
end
