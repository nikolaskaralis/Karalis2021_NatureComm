function logZ = phase_modulation(phases,uniform)

if nargin<2; uniform = 0; end;

if uniform; phases = make_uniform(phases); end;
logZ = log((abs(mean(exp(1i*phases)))^2)*length(phases));
