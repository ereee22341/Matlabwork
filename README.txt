F2-b 一维增透膜课程设计计算程序
================================

运行环境
已在 Windows / MATLAB R2024a 实际运行验证；代码按 R2020b 及以上版本编写，
不依赖额外工具箱。旧版本兼容性未单独实机验证。

运行方法
在 MATLAB 中打开 main.m 并运行。入口会自动定位工程目录，不需要修改绝对路径。
也可以在命令窗口直接输入 main。
在 MATLAB 桌面中运行时，会保留并显示全部 13 张图，默认停靠在图窗区域，
可以切换图窗标签查看。PNG 和可编辑 FIG 文件仍同步保存在 output 文件夹。
图题、图例、坐标轴和图窗名称使用中文，优先选择微软雅黑字体。
物理符号、单位、函数名及数据字段保持标准写法，便于计算和数据处理。
命令行批量运行时默认只导出文件，不打开图窗。
如需在桌面中也只导出图片，可将 config_f2b.m 中的 cfg.showFigures 改为 false。
查看已有结果无需重新计算，例如：
openfig(fullfile('output','05_reflection_comparison.fig'),'visible');
fdtd_spectrum 是计算单次频谱的函数，不是完整项目入口。现在直接输入
s = fdtd_spectrum; 会使用默认 F2-b 参数计算单层频谱，结果保存在变量 s 中；
这种调用不会生成全套图表。生成所有图表请运行 main。
自定义调用示例：
cfg = config_f2b();
s = fdtd_spectrum(cfg.nfilm,cfg.d0,cfg.dz,cfg);
所有结果写入 output/；重复运行会覆盖同名结果。完整计算与绘图在本机约 17 秒。
运行中自动验证解析解、无膜界面、能量守恒、网格收敛、运行时长敏感性和双层膜结果。
报告和 PPT 尚未制作。

文件作用
main.m              一键运行，组织全部计算、检查、数据表和绘图。
config_f2b.m        F2-b 物理参数、网格、时间步、频率点和输出目录。
tmm_spectrum.m      多层介质特征矩阵参考解；采用 exp(+j*w*t) 约定。
fdtd_run.m          自编 Yee 电磁场更新、软源、局部波速 Mur 和介质分配。
fdtd_spectrum.m     参考/结构两次仿真，反射场相减与 FFT 归一化；支持无参数默认调用。
validate_baseline.m 无膜与单层模型的独立解析校验及 FDTD 基础检查。
reflection_band.m   提取包含设计频率的连续低反射频带。
make_figures.m      自动生成所有 PNG 和可编辑 FIG。
print_summary.m     输出中文运行摘要并生成中文结果汇总文件。
ALGORITHM.md        建模、离散公式、边界和误差定义的计算说明。

必交图与输出对应关系（每张同时生成 .png 和 .fig）
结构示意图           output/01_structure
时域波形             output/03_time_waveforms
反射谱对比           output/05_reflection_comparison
透射谱对比           output/06_transmission_comparison
膜厚 ±10% 反射谱     output/08_thickness_sensitivity
网格误差分析         output/09_grid_convergence、10_thickness_rounding
双层膜必做扩展       output/11_double_layer
FDTD 组共同要求       output/02_yee_grid
补充验证图           output/04_field_snapshots、07_energy_balance、
                     12_source_spectrum、13_time_sensitivity

数据文件
output/results.mat            全部主要计算数据、参数、误差和分析表。
output/baseline.mat           无膜和单层基础验证数据。
output/baseline_spectra.csv    单层的 FDTD/TMM 反射谱、电场透射平方、功率透射谱。
output/grid_errors.csv        三档网格的 RMS 误差及最大功率守恒偏差。
output/bandwidths.csv         单层/双层的 FDTD/TMM 低于 1% 反射率带宽。
output/thickness_minima.csv   厚度变化对应的理论/FDTD 反射谷频率。
output/validation_metrics.csv 基础验证、时长敏感性和扩展误差。
output/results_summary.txt    可直接阅读的中文结果汇总；CSV 原有字段名继续保留。

定义与使用注意
程序中 r、t 是电场系数，R=|r|^2，t2=|t|^2，T=3|t|^2 是功率透射率。
无耗时 R+T=1。数值 RMS 为功率比例的绝对误差，不能称为相对误差。
0 GHz 作为频谱极限展示，不参与 RMS；其余 0.02~12 GHz 全部通过激励幅度检查。
膜层按电场节点控制体积的材料占比平均介电常数，基底延伸到右边界。
主单层使用 20 个膜厚网格；收敛采用 10/20/40；双层使用 d0/40 的空间步长。
output/ 为可重建结果，未纳入 Git；源码、参数、公式说明和本说明纳入 Git。
如需恢复：先用 git log --oneline 查看提交，再选择所需版本；保留当前修改后再操作。
