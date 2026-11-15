function [phi, tau, info] = jp_c17_design(p, theta_deg, bandIdx, mode, opt)
%JP_C17_DESIGN  Type-1 JPTA codebook, three candidate design rules.
%
%   [phi,tau] = jp_c17_design(p, theta, bandIdx, mode)
%
%   mode = 'ls'        global optimum of C17 eq. (4), no approximation.
%   mode = 'sawtooth'  delay profile structure of C17 Fig. 8:
%                      tau_m = tauStep * mod(m-1, N_UE), i.e. a sawtooth of
%                      slope one TTD step per element that wraps every N_UE
%                      elements.  Max delay = tauStep*(N_UE-1).
%   mode = 'rainbow'   closed form of eqs. (5)-(6) (Type-2), used as a check.
%
%   In all three cases the phase shifter is set to the optimal value for the
%   given delay vector, phi_m = angle( sum_k conj(b_{k,m}) e^{j2*pi*f_k*tau_m} ),
%   which is the closed-form solution of eq. (4) in phi.
%
%   'ls' implementation note: because the subcarriers are equally spaced,
%   A_m(tau) = sum_k conj(b_{k,m}) e^{j 2 pi f_k tau} is a DFT of the target
%   sequence evaluated on the delay axis, so one zero-padded FFT per element
%   searches the *entire* delay period 1/SCS (8.33 us) exactly, with no grid
%   bias.  nfft is chosen so that the delay resolution is opt.dtau.

if nargin < 5, opt = struct(); end
opt = jp_filldef(opt, struct('dtau', 0.02e-9, 'nPS', p.nPS, 'verbose', false));
if nargin < 4 || isempty(mode), mode = 'ls'; end

M   = p.Nh;
k0  = (0:M-1).';
Nu  = numel(theta_deg);
Ct  = reshape(cosd(theta_deg(:).'), [], 1);
Ct  = reshape(Ct(bandIdx), [], 1);

% ideal steering vector, eq. (3)
B = exp(1j*2*pi*(p.d/p.c) * ((p.f .* Ct) * k0.'));   % Nsc x M

switch lower(mode)
    case 'sawtooth'
        tau = p.tauStep * mod(k0, Nu);

    case 'rainbow'
        dTheta = max(theta_deg) - min(theta_deg);
        thc    = mean(theta_deg);
        tau    = (k0/p.BW) * sind(dTheta/2);
        tau    = p.tauStep * round(tau/p.tauStep);

    case 'ls'
        N   = 2^nextpow2(round(1/(opt.dtau*p.SCS)));
        tauAxis = (0:N-1).'/(N*p.SCS);              % covers [0, 1/SCS)
        tau = zeros(M,1);
        phi = zeros(M,1);
        for m = 1:M
            Y = conj(fft(conj(B(:,m)), N));         % |Y| = |A_m(tau_n)|
            [~, i] = max(abs(Y));
            tau(m) = tauAxis(i);
        end

    otherwise
        error('jp_c17_design: unknown mode %s', mode);
end

% optimal phase for the chosen delays (closed form of eq. 4 in phi)
phi = zeros(M,1);
for m = 1:M
    A = sum(B(:,m) .* exp(1j*2*pi*p.f*tau(m)));
    phi(m) = angle(A);
end
if isfinite(opt.nPS)
    q = 2*pi/opt.nPS;
    phi = q*round(phi/q);
end

% residual of eq. (4) and the achieved in-band gain per UE
P    = exp(1j*(phi.' + 2*pi*(p.f*tau.')));          % Nsc x M
info.obj     = sum(abs(P(:)-B(:)).^2);
info.objNorm = info.obj/(2*p.Nsc);
info.mode    = mode;
info.tauMax  = max(tau);
info.gIn     = zeros(1,Nu);
for u = 1:Nu
    gu = jp_c17_pattern(p, phi, tau, theta_deg(u), p.f(bandIdx==u));
    info.gIn(u) = min(gu(:));
end
if opt.verbose
    fprintf('  %-9s tau_max=%6.2f ns  in-band min gain %s dB\n', ...
            mode, max(tau)*1e9, mat2str(round(info.gIn,1)));
end
end
