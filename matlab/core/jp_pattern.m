function [GdB, AF] = jp_pattern(p, phi, tau, theta_grid, f_vec)
%JP_PATTERN  Beamforming gain (dB) of a JPTA weight set on a fine grid.
%
%   GdB(na, nf) is the total BS gain for a UE at azimuth theta_grid(na) and
%   frequency f_vec(nf):
%
%       GdB = 10*log10(AF_h) + 10*log10(Nv) + Gelem
%
%   with AF_h = |sum_m p_m(f) a_m(theta,f)^*|^2 / Nh in [0, Nh].  The vertical
%   dimension is fixed at broadside (the paper's scenario is 2-D), so it only
%   contributes the constant 10*log10(16); Gelem is set in jp_params so
%   that the 16x16 peak gain is exactly the 28 dB of Table I.
%
%   The phasor sum is accumulated element by element to keep memory at
%   O(na*nf) instead of O(na*nf*Nh).

M   = p.Nh;
k0  = (0:M-1).';
f   = f_vec(:).';
na  = numel(theta_grid);
nf  = numel(f_vec);

P   = exp(1j*(phi(:).' + 2*pi*(f(:)*tau(:).')));      % nf x M
cst = -1j*2*pi*(p.d/p.c);                             % phase per (m, f*cos)
cf  = cosd(theta_grid(:)) * f;                        % na x nf

AF  = zeros(na, nf);
for m = 1:M
    AF = AF + exp(cst*m*cf) .* P(:,m).';              % element m contribution
end
AF  = abs(AF).^2 / M;                                 % peak = M = 16

GdB = 10*log10(max(AF, eps)) + 10*log10(p.Nv) + p.Gelem;
end
