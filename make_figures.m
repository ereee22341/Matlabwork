function make_figures(cfg,num,theory,thickness,thickTheory,grids,gt,...
    rounded,roundTheory,dual,dualTheory,longRun)
% Export high-resolution PNG previews and editable MATLAB FIG files.
showFigures=usejava('desktop');
if isfield(cfg,'showFigures'), showFigures=cfg.showFigures; end
fig=[];
f=cfg.f/1e9; colors=lines(3);
newfig();
subplot(2,1,1); hold on;
patch([0 1 1 0],[0 0 1 1],[.88 .94 1]);
patch([1 2 2 1],[0 0 1 1],[.9 .8 .5]);
patch([2 3 3 2],[0 0 1 1],[.75 .85 .9]);
text(.15,.5,'Air: n = 1'); text(1.08,.65,'Film: n = sqrt(3)');
text(1.08,.35,sprintf('d = %.4f mm',cfg.d0*1e3));
text(2.08,.5,'Substrate: n = 3'); text(.1,1.15,'Incident wave -->');
axis([0 3 -.05 1.3]); axis off; title('Single-layer model (not to scale)');
subplot(2,1,2); plot(num.run.z*1e3,num.run.epsr,'LineWidth',1.5);
xline(num.run.zsrc*1e3,'--','Source'); xline(num.run.zref*1e3,'--','Reflection');
xline(num.run.ztrans*1e3,'--','Transmission');
xlabel('z (mm)'); ylabel('Relative permittivity'); grid on;
title('Computational domain: local-speed Mur at both ends'); finish('01_structure');
newfig();
subplot(2,1,1); plot(0:4,zeros(1,5),'bo',.5:1:3.5,zeros(1,4),'rs','MarkerFaceColor','auto');
legend('E at integer i','H at i+1/2'); ylim([-.5 .5]); yticks([]); xlabel('z / Delta z'); title('Yee spatial staggering');
subplot(2,1,2); plot(0:3,ones(1,4),'bo-',.5:1:2.5,zeros(1,3),'rs-');
yticks([0 1]); yticklabels({'H','E'}); ylim([-.5 1.5]); xlabel('t / Delta t'); grid on; title('Leapfrog time stepping'); finish('02_yee_grid');
newfig(); plot(num.run.time*1e9,[num.reference.er;num.run.er-num.reference.er;num.run.et]','LineWidth',1.3);
xlabel('Time (ns)'); ylabel('Electric field (V/m)'); legend('Incident (left)','Reflected (left)','Transmitted (right)'); grid on;
title('Broadband time-domain waveforms'); finish('03_time_waveforms');
newfig();
for k=1:3
    subplot(3,1,k); plot(num.run.z*1e3,num.run.fields(k,:),'LineWidth',1.2);
    xline(num.run.zfront*1e3,'--'); xline(num.run.zback*1e3,'--');
    xlabel('z (mm)'); ylabel('E_x (V/m)'); grid on;
    title(sprintf('t = %.3f ns',num.run.snapTimes(k)*1e9));
end
finish('04_field_snapshots');
newfig(); plot(f,theory.R,'k-',f,num.R,'r--','LineWidth',1.4);
legend('TMM','FDTD'); xlabel('Frequency (GHz)'); ylabel('Reflectance |r|^2'); grid on; title('Single-layer reflectance'); finish('05_reflection_comparison');
newfig();
subplot(2,1,1); plot(f,theory.t2,'k-',f,num.t2,'r--','LineWidth',1.3);
legend('TMM','FDTD'); ylabel('Field ratio |t|^2'); xlabel('Frequency (GHz)'); grid on;
subplot(2,1,2); plot(f,theory.T,'k-',f,num.T,'r--','LineWidth',1.3);
legend('TMM','FDTD'); ylabel('Power transmittance 3|t|^2'); xlabel('Frequency (GHz)'); grid on; finish('06_transmission_comparison');
newfig();
subplot(2,1,1); plot(f,num.R+num.T,'LineWidth',1.3); yline(1,'k--');
xlabel('Frequency (GHz)'); ylabel('R + T'); grid on; title('Lossless energy balance');
subplot(2,1,2); plot(f,100*(num.R+num.T-1),'LineWidth',1.3);
xlabel('Frequency (GHz)'); ylabel('Signed residual (%)'); grid on; finish('07_energy_balance');
newfig(); hold on;
for k=1:3
    plot(f,thickTheory{k}.R,'-','Color',colors(k,:),'LineWidth',1.4,'DisplayName',sprintf('TMM d/d0=%.1f',.8+.1*k));
    plot(f,thickness{k}.R,'--','Color',colors(k,:),'LineWidth',1,'DisplayName',sprintf('FDTD d/d0=%.1f',.8+.1*k));
end
xlabel('Frequency (GHz)'); ylabel('Reflectance |r|^2'); legend('Location','best'); grid on; title('Thickness tolerance'); finish('08_thickness_sensitivity');
newfig(); subplot(2,1,1); hold on;
for k=1:3, plot(f,grids{k}.R-theory.R,'DisplayName',sprintf('%d film cells',gt.FilmCells(k))); end
xlabel('Frequency (GHz)'); ylabel('Reflectance error'); grid on; legend;
subplot(2,1,2); loglog(gt.dz_mm,gt.R_RMS,'o-',gt.dz_mm,gt.T_RMS,'s-','LineWidth',1.3);
xlabel('Grid spacing (mm)'); ylabel('RMS error (power fraction)'); legend('Reflectance','Transmittance'); grid on; finish('09_grid_convergence');
newfig(); plot(f,theory.R,'k-',f,roundTheory.R,'b-',f,rounded.R,'r--','LineWidth',1.3);
legend('TMM: design thickness','TMM: rounded thickness','FDTD: rounded thickness');
xlabel('Frequency (GHz)'); ylabel('Reflectance |r|^2'); grid on; title('Geometry rounding vs discretization'); finish('10_thickness_rounding');
newfig(); subplot(2,1,1); plot(f,theory.R,'b-',f,num.R,'b--',f,dualTheory.R,'r-',f,dual.R,'r--','LineWidth',1.3); hold on;
yline(.01,'k:','1% threshold'); legend('Single TMM','Single FDTD','Double TMM','Double FDTD','Location','best');
xlabel('Frequency (GHz)'); ylabel('Reflectance |r|^2'); grid on; title('Single vs double quarter-wave coating');
subplot(2,1,2); plot(f,theory.R,'b-',f,num.R,'b--',f,dualTheory.R,'r-',f,dual.R,'r--','LineWidth',1.3); hold on;
yline(.01,'k:'); bs=reflection_band(cfg.f,num.R,cfg.f0,.01); bd=reflection_band(cfg.f,dual.R,cfg.f0,.01);
plot(bs(1:2)/1e9,[.01 .01],'bo','MarkerFaceColor','b');
plot(bd(1:2)/1e9,[.01 .01],'ro','MarkerFaceColor','r');
xlim([4 8]); ylim([0 .025]); grid on;
xlabel('Frequency (GHz)'); ylabel('Reflectance |r|^2');
title(sprintf('FDTD fractional bandwidth: single %.2f%%; double %.2f%%',100*bs(3),100*bd(3)));
finish('11_double_layer');
newfig(); plot(f,num.sourceLevel,'LineWidth',1.3); yline(1e-3,'k--');
xlabel('Frequency (GHz)'); ylabel('Normalized spectral amplitude'); grid on; title('Incident source bandwidth'); finish('12_source_spectrum');
newfig(); plot(f,longRun.R-num.R,f,longRun.T-num.T,'LineWidth',1.3);
legend('Reflectance change','Transmittance change'); xlabel('Frequency (GHz)'); ylabel('Change (power fraction)'); grid on;
title('Duration sensitivity: 1.85 ns vs 2.15 ns'); finish('13_time_sensitivity');
    function newfig()
        fig=figure('Color','w','Position',[100 100 1000 700],'Visible','off');
    end
    function finish(name)
        set(fig,'Name',name,'NumberTitle','off');
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
