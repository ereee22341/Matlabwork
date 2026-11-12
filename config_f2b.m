function cfg = config_f2b()
% F2-b: lossless, nonmagnetic, nondispersive media; SI units internally.
cfg.c0 = 2.998e8;
cfg.mu0 = 4*pi*1e-7;
cfg.eps0 = 1/(cfg.mu0*cfg.c0^2); % Consistent constants (8.854e-12 rounded).
cfg.f0 = 6e9;
cfg.n0 = 1;
cfg.ns = 3;
cfg.nfilm = sqrt(3);
cfg.d0 = cfg.c0/(4*cfg.f0*cfg.nfilm);
cfg.f = linspace(0,2*cfg.f0,601);
cfg.S = 0.9;
cfg.dz = cfg.d0/20;
cfg.tau = 35e-12;
cfg.t0 = 6*cfg.tau;
cfg.tend = 1.85e-9;
cfg.outdir = fullfile(fileparts(mfilename('fullpath')),'output');
end
