function make_figures(cfg,num,theory,thickness,thickTheory,grids,gt,...
    rounded,roundTheory,dual,dualTheory,longRun)
% 生成中文图表，并导出高清 PNG 与可编辑 FIG 文件。
showFigures=usejava('desktop');
if isfield(cfg,'showFigures'), showFigures=cfg.showFigures; end
% 对图内所有文字统一设置中文字体，避免导出图片时出现方框。
fontName='Microsoft YaHei';
availableFonts=listfonts;
for candidate={'Microsoft YaHei','SimHei','SimSun','Arial Unicode MS'}
    if any(strcmpi(availableFonts,candidate{1})), fontName=candidate{1}; break; end
end
figureTitles={'结构示意图','Yee 网格示意图','时域波形','电场空间快照',...
    '反射谱对比','透射谱对比','能量守恒检验','膜厚偏差分析',...
    '网格收敛分析','膜厚取整误差分析','双层增透膜扩展',...
    '激励频谱覆盖范围','仿真时长敏感性分析'};
fig=[];
f=cfg.f/1e9; colors=lines(3);
newfig();
subplot(2,1,1); hold on;
patch([0 1 1 0],[0 0 1 1],[.88 .94 1]);
patch([1 2 2 1],[0 0 1 1],[.9 .8 .5]);
patch([2 3 3 2],[0 0 1 1],[.75 .85 .9]);
text(.15,.5,'空气：n = 1'); text(1.08,.65,'增透膜：n = \surd3');
text(1.08,.35,sprintf('d = %.4f mm',cfg.d0*1e3));
text(2.08,.5,'半无限基底：n = 3'); text(.1,1.15,'入射方向 →');
axis([0 3 -.05 1.3]); axis off; title('单层四分之一波长增透膜结构（非等比例示意）');
subplot(2,1,2); plot(num.run.z*1e3,num.run.epsr,'LineWidth',1.5);
xline(num.run.zsrc*1e3,'--','激励源'); xline(num.run.zref*1e3,'--','反射观察点');
xline(num.run.ztrans*1e3,'--','透射观察点');
xlabel('位置 z (mm)'); ylabel('相对介电常数 \epsilon_r'); grid on;
title('计算区域：两端采用当地介质波速的 Mur 吸收边界'); finish('01_structure');
newfig();
subplot(2,1,1); plot(0:4,zeros(1,5),'bo',.5:1:3.5,zeros(1,4),'rs','MarkerFaceColor','auto');
legend('电场 E：整数节点 i','磁场 H：半整数节点 i+1/2'); ylim([-.5 .5]); yticks([]); xlabel('归一化位置 z / \Delta z'); title('Yee 网格的空间交错分布');
subplot(2,1,2); plot(0:3,ones(1,4),'bo-',.5:1:2.5,zeros(1,3),'rs-');
yticks([0 1]); yticklabels({'磁场 H','电场 E'}); ylim([-.5 1.5]); xlabel('归一化时间 t / \Delta t'); grid on; title('电场与磁场的蛙跳式时间推进'); finish('02_yee_grid');
newfig(); plot(num.run.time*1e9,[num.reference.er;num.run.er-num.reference.er;num.run.et]','LineWidth',1.3);
xlabel('时间 t (ns)'); ylabel('电场 E_x (V/m)'); legend('入射场（左观察点）','反射场（左观察点）','透射场（右观察点）'); grid on;
title('宽带脉冲的入射、反射与透射时域波形'); finish('03_time_waveforms');
newfig();
for k=1:3
    subplot(3,1,k); plot(num.run.z*1e3,num.run.fields(k,:),'LineWidth',1.2);
    xline(num.run.zfront*1e3,'--'); xline(num.run.zback*1e3,'--');
    xlabel('位置 z (mm)'); ylabel('电场 E_x (V/m)'); grid on;
    title(sprintf('电场空间分布：t = %.3f ns',num.run.snapTimes(k)*1e9));
end
finish('04_field_snapshots');
newfig(); plot(f,theory.R,'k-',f,num.R,'r--','LineWidth',1.4);
legend('TMM 参考解','FDTD 数值解'); xlabel('频率 f (GHz)'); ylabel('功率反射率 |r|^2'); grid on; title('单层增透膜反射谱对比'); finish('05_reflection_comparison');
newfig();
subplot(2,1,1); plot(f,theory.t2,'k-',f,num.t2,'r--','LineWidth',1.3);
legend('TMM 参考解','FDTD 数值解'); ylabel('电场透射系数平方 |t|^2'); xlabel('频率 f (GHz)'); grid on; title('电场透射系数平方');
subplot(2,1,2); plot(f,theory.T,'k-',f,num.T,'r--','LineWidth',1.3);
legend('TMM 参考解','FDTD 数值解'); ylabel('功率透射率 3|t|^2'); xlabel('频率 f (GHz)'); grid on; title('按介质波阻抗归一化的功率透射率'); finish('06_transmission_comparison');
newfig();
subplot(2,1,1); plot(f,num.R+num.T,'LineWidth',1.3); yline(1,'k--');
xlabel('频率 f (GHz)'); ylabel('反射率与透射率之和 R + T'); grid on; title('无耗介质的能量守恒检验');
subplot(2,1,2); plot(f,100*(num.R+num.T-1),'LineWidth',1.3);
xlabel('频率 f (GHz)'); ylabel('能量守恒有符号偏差 (%)'); grid on; finish('07_energy_balance');
newfig(); hold on;
for k=1:3
    plot(f,thickTheory{k}.R,'-','Color',colors(k,:),'LineWidth',1.4,'DisplayName',sprintf('TMM 参考解：d/d_0=%.1f',.8+.1*k));
    plot(f,thickness{k}.R,'--','Color',colors(k,:),'LineWidth',1,'DisplayName',sprintf('FDTD 数值解：d/d_0=%.1f',.8+.1*k));
end
xlabel('频率 f (GHz)'); ylabel('功率反射率 |r|^2'); legend('Location','best'); grid on; title('膜厚偏差为 -10%、0、+10% 时的反射谱'); finish('08_thickness_sensitivity');
newfig(); subplot(2,1,1); hold on;
for k=1:3, plot(f,grids{k}.R-theory.R,'DisplayName',sprintf('膜厚内 %d 个网格',gt.FilmCells(k))); end
xlabel('频率 f (GHz)'); ylabel('反射率绝对误差'); grid on; legend; title('不同网格密度下的反射率误差');
subplot(2,1,2); loglog(gt.dz_mm,gt.R_RMS,'o-',gt.dz_mm,gt.T_RMS,'s-','LineWidth',1.3);
xlabel('网格间距 \Delta z (mm)'); ylabel('均方根误差（功率比例）'); legend('功率反射率','功率透射率'); grid on; title('网格加密的收敛性'); finish('09_grid_convergence');
newfig(); plot(f,theory.R,'k-',f,roundTheory.R,'b-',f,rounded.R,'r--','LineWidth',1.3);
legend('TMM：设计膜厚','TMM：取整后的膜厚','FDTD：取整后的膜厚');
xlabel('频率 f (GHz)'); ylabel('功率反射率 |r|^2'); grid on; title('区分膜厚取整误差与数值离散误差'); finish('10_thickness_rounding');
newfig(); subplot(2,1,1); plot(f,theory.R,'b-',f,num.R,'b--',f,dualTheory.R,'r-',f,dual.R,'r--','LineWidth',1.3); hold on;
yline(.01,'k:','1% 反射率阈值'); legend('单层膜 TMM','单层膜 FDTD','双层膜 TMM','双层膜 FDTD','Location','best');
xlabel('频率 f (GHz)'); ylabel('功率反射率 |r|^2'); grid on; title('单层与双层四分之一波长增透膜对比');
subplot(2,1,2); plot(f,theory.R,'b-',f,num.R,'b--',f,dualTheory.R,'r-',f,dual.R,'r--','LineWidth',1.3); hold on;
yline(.01,'k:'); bs=reflection_band(cfg.f,num.R,cfg.f0,.01); bd=reflection_band(cfg.f,dual.R,cfg.f0,.01);
plot(bs(1:2)/1e9,[.01 .01],'bo','MarkerFaceColor','b');
plot(bd(1:2)/1e9,[.01 .01],'ro','MarkerFaceColor','r');
xlim([4 8]); ylim([0 .025]); grid on;
xlabel('频率 f (GHz)'); ylabel('功率反射率 |r|^2');
title(sprintf('FDTD 低反射相对带宽：单层膜 %.2f%%，双层膜 %.2f%%',100*bs(3),100*bd(3)));
finish('11_double_layer');
newfig(); plot(f,num.sourceLevel,'LineWidth',1.3); yline(1e-3,'k--');
xlabel('频率 f (GHz)'); ylabel('归一化入射频谱幅度'); grid on; title('宽带激励的频谱覆盖范围'); finish('12_source_spectrum');
newfig(); plot(f,longRun.R-num.R,f,longRun.T-num.T,'LineWidth',1.3);
legend('功率反射率变化','功率透射率变化'); xlabel('频率 f (GHz)'); ylabel('变化量（功率比例）'); grid on;
title(sprintf('仿真时长敏感性：%.2f ns 与 %.2f ns',cfg.tend*1e9,longRun.run.time(end)*1e9)); finish('13_time_sensitivity');
    function newfig()
        fig=figure('Color','w','Position',[100 100 1000 700],'Visible','off');
    end
    function finish(name)
        index=sscanf(name,'%d',1);
        set(fig,'Name',sprintf('%02d %s',index,figureTitles{index}),'NumberTitle','off');
        set(findall(fig,'-property','FontName'),'FontName',fontName);
        set(findall(fig,'Type','axes'),'FontSize',11);
        exportgraphics(fig,fullfile(cfg.outdir,[name,'.png']),'Resolution',180);
        savefig(fig,fullfile(cfg.outdir,[name,'.fig']));
        if showFigures
            if usejava('desktop'), set(fig,'WindowStyle','docked'); end
            set(fig,'Visible','on');
            drawnow;
        else
            close(fig);
        end
    end
end
