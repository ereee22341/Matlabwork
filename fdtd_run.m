function run = fdtd_run(indices, thicknesses, dz, cfg, reference)
% 一维 Ex/Hy Yee 求解器：电场位于整数时空节点，磁场位于半整数节点。
% 使用电场软源、当地介质波速的一阶 Mur 边界及节点控制体积平均。
if nargin<5, reference=false; end
validateattributes(cfg.S,{'numeric'},{'real','finite','scalar','positive','<=',1},mfilename,'cfg.S');
dt=cfg.S*dz/cfg.c0;
snap=@(x) round(x/dz)*dz;
zsrc=snap(0.35); zref=snap(0.40); zfront=snap(0.45);
zback=zfront+sum(thicknesses);
ztrans=snap(zback+0.04);
z=0:dz:snap(ztrans+0.12);
N=numel(z); Nt=ceil(cfg.tend/dt)+1; time=(0:Nt-1)*dt;
isrc=round(zsrc/dz)+1; ir=round(zref/dz)+1; it=round(ztrans/dz)+1;
epsr=ones(1,N);
if ~reference
    faces=zfront+[0,cumsum(thicknesses)];
    for k=1:numel(indices)
        fraction=max(0,min(z+dz/2,faces(k+1))-max(z-dz/2,faces(k)))/dz;
        epsr=epsr+(indices(k)^2-1)*fraction;
    end
    fraction=max(0,min(1,(z+dz/2-zback)/dz));
    epsr=epsr+(cfg.ns^2-1)*fraction;
end
ce=dt./(cfg.eps0*epsr*dz); ch=dt/(cfg.mu0*dz);
vleft=cfg.c0/sqrt(epsr(1)); vright=cfg.c0/sqrt(epsr(end));
ql=(vleft*dt-dz)/(vleft*dt+dz); qr=(vright*dt-dz)/(vright*dt+dz);
E=zeros(1,N); H=zeros(1,N-1); er=zeros(1,Nt); et=er;
snapTimes=[0.40,0.65,0.95]*1e-9;
snapSteps=round(snapTimes/dt)+1; fields=zeros(numel(snapSteps),N);
source=2*cfg.S*exp(-((time-cfg.t0)/cfg.tau).^2);
% 首个记录为 E(t=0)=0，后续场更新与记录时刻一一对应。
for j=2:Nt
    oldLeft=E(1:2); oldRight=E(end-1:end);
    H=H-ch*diff(E);
    E(2:end-1)=E(2:end-1)-ce(2:end-1).*diff(H);
    E(isrc)=E(isrc)+source(j);
    E(1)=oldLeft(2)+ql*(E(2)-oldLeft(1));
    E(end)=oldRight(1)+qr*(E(end-1)-oldRight(2));
    er(j)=E(ir); et(j)=E(it);
    k=find(snapSteps==j,1);
    if ~isempty(k), fields(k,:)=E; end
end
run.time=time; run.er=er; run.et=et; run.z=z; run.epsr=epsr;
run.fields=fields; run.snapTimes=(snapSteps-1)*dt;
run.dt=dt; run.dz=dz; run.zsrc=zsrc; run.zref=zref;
run.zfront=zfront; run.zback=zback; run.ztrans=ztrans;
end
