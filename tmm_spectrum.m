function s = tmm_spectrum(f, indices, thicknesses, cfg)
% Characteristic matrices for exp(+j*w*t); order is air -> substrate.
assert(numel(indices)==numel(thicknesses));
A=ones(size(f)); B=zeros(size(f)); C=B; D=A;
for k=1:numel(indices)
    n=indices(k); delta=2*pi*f*n*thicknesses(k)/cfg.c0;
    a=cos(delta); b=1i*sin(delta)/n; c=1i*n*sin(delta);
    [A,B,C,D]=deal(A.*a+B.*c,A.*b+B.*a,C.*a+D.*c,C.*b+D.*a);
end
den=cfg.n0*(A+B*cfg.ns)+C+D*cfg.ns;
s.r=(cfg.n0*(A+B*cfg.ns)-C-D*cfg.ns)./den;
s.t=2*cfg.n0./den;
s.R=abs(s.r).^2;
s.t2=abs(s.t).^2;
s.T=cfg.ns/cfg.n0*s.t2;
s.f=f;
end
