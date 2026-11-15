function [phi, tau, info] = jp_c17_type1(p, theta_deg, bandIdx, opt)
%JP_C17_TYPE1  Type-1 JPTA codebook by solving C17 eq. (4) exactly.
%
%   The paper's problem is
%       min_{tau_m,phi_m}  sum_k || p_k - b_k ||^2 ,                       (4)
%   with p_k the JPTA weight vector and b_k the ideal steering vector of
%   eq. (3).  Because both vectors are per-antenna, the objective decouples
%   completely over m:
%
%       sum_k sum_m | p_{k,m} - b_{k,m} |^2
%         = const - 2 Re{ sum_m e^{j phi_m} * sum_k e^{j 2 pi f_k tau_m} conj(b_{k,m}) } .
%
%   Hence for a given tau_m the optimal phase is closed form,
%       phi_m = angle( A_m(tau_m) ),   A_m(tau) = sum_k e^{j 2 pi f_k tau} conj(b_{k,m}),
%   and the remaining scalar problem is a one-dimensional search
%       max_tau |A_m(tau)| .
%   This is the *global* optimum of eq. (4), not an alternating heuristic:
%   no iteration is needed and the achieved normalised objective can be
%   reported as the exact optimum.  (The paper says "solved by an iterative
%   algorithm proposed in [6]"; the decoupling makes iteration unnecessary.)

if nargin < 4, opt = struct(); end
opt = jp_filldef(opt, struct('tauMax', p.tauMax, 'tauStep', p.tauStep, ...
                             'nPS', p.nPS, 'verbose', false));

M   = p.Nh;
k0  = (0:M-1).';                                  % element index m
cosT = cosd(theta_deg(:).');
Ct   = reshape(cosT(bandIdx), [], 1);              % Nsc x 1 (one angle per subcarrier)

% ideal steering vector at (theta_{i(k)}, f_k):  eq. (3)
TH = 2*pi*(p.d/p.c) * ((p.f .* Ct) * k0.');        % Nsc x M
B  = exp(1j*TH);

% delay grid (hardware resolution 2.5 ns)
taus = (0:opt.tauStep:opt.tauMax).';
Nt   = numel(taus);

% A(tau, m) = sum_k exp(j2*pi*f_k*tau) * conj(B(k,m))
E = exp(1j*2*pi*(p.f*taus.'));                    % Nsc x Nt
A = E.' * conj(B);                                % Nt x M

[~, ib] = max(abs(A), [], 1);
tau = reshape(taus(ib), M, 1);                         % M x 1
phi = reshape(angle(A(sub2ind([Nt M], ib, 1:M))), M, 1);

if isfinite(opt.nPS)
    q = 2*pi/opt.nPS;
    phi = q*round(phi/q);
end

% achieved objective and the ideal (b == p) value, for reporting
P    = exp(1j*(phi.' + 2*pi*(p.f*tau.')));        % Nsc x M
obj  = sum(abs(P(:)-B(:)).^2);
obj0 = 0;                                          % sum_k ||b_k||^2 = Nsc
info.obj     = obj;
info.objNorm = obj/(2*p.Nsc);                      % 0 => perfect fit, 1 => orthogonal
info.tauMax  = max(tau);
info.tauMin  = min(tau);
info.fit     = 10*log10(1 - info.objNorm/1);       % rough dB view of the residual
if opt.verbose
    fprintf('  Nu=%d: obj/Nsc = %.4f, tau in [%.2f, %.2f] ns\n', ...
            numel(theta_deg), obj/p.Nsc, min(tau)*1e9, max(tau)*1e9);
end
end
