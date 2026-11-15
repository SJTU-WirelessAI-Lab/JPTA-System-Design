% =======================================================================
%  run_b1d.m  --  how fine must the TTD step be for this paper's "<=3 dB" claim?
%
%  Fix the hardware constraint (one PS + one TTD per element), take the
%  Type-1 delay profile of eq. (5) wrapped every N_UE elements,
%      tau_m = dt * mod(m-1, N_UE),
%  optimise only the phase (closed form of eq. 4), and sweep dt continuously.
%  For each dt record the worst in-band gain over all UEs.
%
%  dt is the only knob: the TTD grid of Table I allows dt = 2.5*k ns only.
% =======================================================================
addpath(fullfile(fileparts(mfilename('fullpath')),'core'));
addpath(fullfile(fileparts(mfilename('fullpath')),'alg'));
p   = jp_params();
out = fullfile(fileparts(mfilename('fullpath')), 'out');
if ~exist(out,'dir'), mkdir(out); end
fid = fopen(fullfile(out,'b1d_summary.txt'),'w');

% ---- the eq.(5) slope is never a multiple of the TTD grid ---------------
fprintf(fid,'eq.(5) slope vs the 2.5 ns TTD grid\n');
fprintf(fid,'  %-8s %12s %12s %12s\n','case','sin(d/2)/W','/2.5 ns','on grid?');
for c = 1:numel(p.caseAngles)
    th = p.caseAngles{c};
    d5 = sind((max(th)-min(th))/2)/p.BW;
    on = 'NO';
    if abs(d5/p.tauStep - round(d5/p.tauStep)) < 1e-9, on = 'yes'; end
    fprintf(fid,'  %-8s %12.3f %12.3f %12s\n', p.caseName{c}, d5*1e9, ...
            d5/p.tauStep, on);
end
fprintf(fid,'\nno case sits on the grid, so every design is quantised away from eq.(5).\n\n');

dts = logspace(log10(0.02), log10(3), 60)*1e-9;   % 0.02 ... 3 ns
G   = zeros(numel(p.caseAngles), numel(dts));
fig = figure('Visible','off','Position',[80 80 900 560]); hold on;
mk  = {'-o','-s','-^','-d'}; cl = lines(4);
for c = 1:numel(p.caseAngles)
    th = p.caseAngles{c}; Nu = numel(th); b = jp_bands(p, Nu);
    Ws = (b.edges(2)-b.edges(1))*p.SCS;
    for i = 1:numel(dts)
        tau = dts(i)*mod((0:p.Nh-1).', Nu);
        phi = opt_phi(p, th, b.idx, tau);
        gW  = zeros(1,Nu);
        for u = 1:Nu
            gu = jp_pattern(p, phi, tau, th(u), p.f(b.idx==u));
            gW(u) = min(gu(:));
        end
        G(c,i) = min(gW);
    end
    plot(dts*1e9, G(c,:), mk{c}, 'Color',cl(c,:), 'LineWidth',1.7, ...
         'MarkerFaceColor',cl(c,:), 'MarkerSize',4);
    dmax = dt_cross(dts, G(c,:));
    fprintf(fid,'%-8s Ws=%5.0f MHz  best sampled dt keeping 25 dB = %6.3f ns | at eq.(5) slope = %6.2f dB\n', ...
            p.caseName{c}, Ws/1e6, dmax, worst_gain(p, th, b.idx, sind((max(th)-min(th))/2)/p.BW));
end
plot([2.5 2.5], [0 30], 'k--', 'LineWidth',1.4);
text(2.6, 27, 'Table I: TTD step = 2.5 ns', 'FontSize',9);
yline(p.Gbs-3,'r:','LineWidth',1.3);
text(0.03, p.Gbs-2.2, '28 - 3 dB', 'Color','r', 'FontSize',9);
set(gca,'XScale','log'); xlim([0.02 3]); ylim([-10 30]);
xlabel('TTD step dt (ns)   [Table I allows 2.5\cdotk only]');
ylabel('worst in-band gain over all UEs (dB)');
title('Fig.7 "<=3 dB below 28 dB" vs TTD resolution');
legend(p.caseName,'Location','southwest'); grid on; box on;
print(fig, fullfile(out,'b1d_step_sweep.png'), '-dpng','-r170'); close(fig);

% ---- the 2.5 ns grid vs the exact eq.(5) step, per case -----------------
fprintf(fid,'\nwhat the two candidate steps give:\n');
fprintf(fid,'%-8s %10s %10s %10s %10s\n','case','eq5 step','eq5q step','g@eq5','g@eq5q');
for c = 1:numel(p.caseAngles)
    th = p.caseAngles{c}; Nu = numel(th); b = jp_bands(p, Nu);
    d5  = sind((max(th)-min(th))/2)/p.BW;
    d5q = ceil(d5/p.tauStep)*p.tauStep;
    g5  = worst_gain(p, th, b.idx, d5);
    g5q = worst_gain(p, th, b.idx, d5q);
    fprintf(fid,'%-8s %10.3f %10.3f %10.2f %10.2f\n', p.caseName{c}, d5*1e9, d5q*1e9, g5, g5q);
end
fclose(fid);
type(fullfile(out,'b1d_summary.txt'));

% =======================================================================
function phi = opt_phi(p, th, bandIdx, tau)
M  = p.Nh; k0 = (0:M-1).';
Ct = reshape(cosd(th(:).'),[],1); Ct = reshape(Ct(bandIdx),[],1);
B  = exp(1j*2*pi*(p.d/p.c)*((p.f.*Ct)*k0.'));
phi = zeros(M,1);
for m = 1:M
    phi(m) = angle(sum(B(:,m).*exp(1j*2*pi*p.f*tau(m))));
end
end

function g = worst_gain(p, th, bandIdx, dt)
Nu  = numel(th); tau = dt*mod((0:p.Nh-1).', Nu);
phi = opt_phi(p, th, bandIdx, tau);
gW  = zeros(1,Nu);
for u = 1:Nu
    gu = jp_pattern(p, phi, tau, th(u), p.f(bandIdx==u));
    gW(u) = min(gu(:));
end
g = min(gW);
end

function d = dt_cross(dts, g)
% largest dt at which the worst gain is still >= 25 dB
ok = find(g >= 25, 1, 'last');
if isempty(ok), d = 0; else, d = dts(ok)*1e9; end
end
