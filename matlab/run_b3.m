% =======================================================================
%  run_b3.m  --  Fig. 8 (delay law) + the two structural checks
%    1) delay law  tau_m = (1/W)*mod(m-1, N_UE)  -> max delay (N_UE-1)/W
%    2) symmetry   |AF(theta)| == |AF(-theta)| for any JPTA weight set
%    3) delay needed to steer a *symmetric* user set (analytic, eq. 5)
% =======================================================================
addpath(fullfile(fileparts(mfilename('fullpath')),'core'));
addpath(fullfile(fileparts(mfilename('fullpath')),'alg'));
p   = jp_params();
out = fullfile(fileparts(mfilename('fullpath')), 'out');
if ~exist(out,'dir'), mkdir(out); end
fid = fopen(fullfile(out,'b3_summary.txt'),'w');

% ---------------- 1) delay law ------------------------------------------
nUE  = [2 4 8 16];
tauS = zeros(numel(nUE), p.Nh);
for i = 1:numel(nUE)
    tauS(i,:) = (1/p.BW)*mod(0:p.Nh-1, nUE(i));
end
fprintf(fid,'[1] delay law  tau_m = (1/W) * mod(m-1, N_UE),  1/W = %.2f ns\n', 1e9/p.BW);
fprintf(fid,'  %-6s %10s %10s %8s\n','N_UE','mine (ns)','paper (ns)','err');
paper = [2.5 7.5 17.5 35];
for i = 1:numel(nUE)
    fprintf(fid,'  %-6d %10.2f %10.2f %8.2f\n', nUE(i), max(tauS(i,:))*1e9, ...
            paper(i), max(tauS(i,:))*1e9 - paper(i));
end
fprintf(fid,'  teeth per array = 16/N_UE: %s\n\n', mat2str(16./nUE));

fg = figure('Visible','off','Position',[100 100 620 460]); hold on;
mk = {'-o','-v','-s','-d'}; cl = lines(4);
for i = 1:numel(nUE)
    plot(1:p.Nh, tauS(i,:)*1e9, mk{i}, 'Color',cl(i,:), 'LineWidth',1.8, ...
         'MarkerFaceColor',cl(i,:), 'MarkerSize',6);
end
xlabel('Delay element index'); ylabel('Delay (ns)');
title('Fig.8 reproduction: delay per element');
legend(arrayfun(@(n) sprintf('%d UEs',n), nUE, 'UniformOutput',false), 'Location','northwest');
grid on; box on; ylim([0 38]);
print(fg, fullfile(out,'b3_fig8_delay.png'), '-dpng', '-r180'); close(fg);

% ---------------- 2) symmetry check -------------------------------------
% any weight set (phi, tau) gives |AF(theta,f)| = |AF(-theta,f)|
rng(1);
th = p.caseAngles{2};
[phi, tau, ~] = jp_design(p, th, jp_bands(p,numel(th)).idx, 'sawtooth');
ang = -60:0.25:60;
G1 = jp_pattern(p, phi, tau, ang,  p.f(1:20:end));
G2 = jp_pattern(p, phi, tau, -ang, p.f(1:20:end));
fprintf(fid,'[2] symmetry of the array factor (any phi,tau, any frequency)\n');
fprintf(fid,'  max |G(theta)-G(-theta)| = %.3e dB over %d angles x %d freqs\n', ...
        max(abs(G1(:)-G2(:))), numel(ang), size(G1,2));

% deliberate test with a random weight set as well
phi2 = 2*pi*rand(p.Nh,1); tau2 = 40e-9*rand(p.Nh,1);
Ga = jp_pattern(p, phi2, tau2, ang, p.f(1:40:end));
Gb = jp_pattern(p, phi2, tau2,-ang, p.f(1:40:end));
fprintf(fid,'  same for random weights: %.3e dB\n\n', max(abs(Ga(:)-Gb(:))));

% ---------------- 3) delay a symmetric user set actually needs ------------
fprintf(fid,'[3] delay needed to sweep between two symmetric angles (eq. 5)\n');
fprintf(fid,'  %-8s %10s %12s %12s\n','case','dTheta/2','2sin (cos span)','tau_span (ns)');
for i = 1:numel(p.caseAngles)
    th2 = p.caseAngles{i};
    dth = max(th2) - min(th2);
    cspan = 2*sind(dth/2);          % cos range covered by a monotone sweep
    fprintf(fid,'  %-8s %10.1f %12.3f %12.2f\n', p.caseName{i}, dth/2, cspan, ...
            15*cspan/(2*p.BW)*1e9);
end
fprintf(fid,'  (the paper uses 2.5/7.5/17.5/35 ns, i.e. the sawtooth law, not eq. 5)\n');
fclose(fid);
type(fullfile(out,'b3_summary.txt'));
