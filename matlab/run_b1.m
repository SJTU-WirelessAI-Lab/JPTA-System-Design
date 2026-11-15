% =======================================================================
%  run_b1.m  --  C17 Fig. 7 (Type-1 beam patterns) and Fig. 8 (delay profile)
%  Run on nanyangnorth:  matlab -batch run_b1
% =======================================================================
addpath(fullfile(fileparts(mfilename('fullpath')),'core'));
addpath(fullfile(fileparts(mfilename('fullpath')),'alg'));

p   = jp_params_c17();
out = fullfile(fileparts(mfilename('fullpath')), 'out');
if ~exist(out,'dir'), mkdir(out); end

nCase = numel(p.caseAngles);
angV  = -60:0.25:60;                       % fine azimuth grid for the heatmap
fo   = p.foff;                             % subcarrier offsets

Dmax  = zeros(1,nCase);
fid = fopen(fullfile(out,'b1_summary.txt'),'w');

for c = 1:nCase
    th  = p.caseAngles{c};
    Nu  = numel(th);
    b   = jp_c17_bands(p, Nu);

    t1 = tic;
    [phi, tau, info] = jp_c17_type1(p, th, b.idx, struct('verbose',false));
    dt = toc(t1);

    % ---- beam pattern over the whole band -------------------------------
    G = jp_c17_pattern(p, phi, tau, angV, p.f);

    % gain inside each UE's own subband, at its own angle
    gIn = zeros(1,Nu);
    for u = 1:Nu
        sel = (b.idx == u);
        gu  = jp_c17_pattern(p, phi, tau, th(u), p.f(sel));
        gIn(u) = min(gu(:));
    end
    Dmax(c) = max(tau)*1e9;

    fprintf(fid,'=== %s: angles %s\n', p.caseName{c}, mat2str(th));
    fprintf(fid,'  design time %.2f s, eq.(4) obj/Nsc = %.4f\n', dt, info.obj);
    fprintf(fid,'  delay (ns), element 1..16: %s\n', mat2str(round(tau.'*1e9,2)));
    fprintf(fid,'  max delay %.2f ns, min delay %.2f ns\n', Dmax(c), min(tau)*1e9);
    fprintf(fid,'  worst in-band gain per UE (dB): %s   (peak = 28 dB)\n', ...
            mat2str(round(gIn,2)));
    fprintf(fid,'  global max gain %.2f dB\n\n', max(G(:)));

    % ---- Fig.7 panel ----------------------------------------------------
    fg = figure('Visible','off','Position',[100 100 720 480]);
    imagesc(fo/1e8, angV, G); axis xy; hold on;
    for i = 1:Nu
        plot([fo(b.edges(i)+1) fo(b.edges(i+1))]/1e8, [th(i) th(i)], 'w--', 'LineWidth',1.2);
    end
    cb = colorbar; cb.Label.String = 'Beamforming gain (dB)';
    caxis([12 28]);
    xlabel('Frequency offset (10^8 Hz)'); ylabel('Angle (degree)');
    title(sprintf('C17 Fig.7(%s): Type-1 JPTA, %d UEs, peak 28 dB', ...
          char('a'+c-1), Nu));
    print(fg, fullfile(out, sprintf('b1_fig7_%due.png', Nu)), '-dpng', '-r170');
    close(fg);

    if c == 1
        S.phi = phi; S.tau = tau; S.G = G; S.angV = angV;
    else
        S(c).phi = phi; S(c).tau = tau; S(c).G = G; S(c).angV = angV; %#ok<SAGROW>
    end
end

% ---- Fig.8 : delay vs element index ------------------------------------
fid2 = fopen(fullfile(out,'b1_fig8.csv'),'w');
fprintf(fid2,'element,2UE,4UE,8UE,16UE\n');
for m = 1:p.Nh
    fprintf(fid2,'%d', m-1);
    for c = 1:nCase, fprintf(fid2,',%.2f', S(c).tau(m)*1e9); end
    fprintf(fid2,'\n');
end
fclose(fid2);

fg = figure('Visible','off','Position',[100 100 620 460]); hold on;
mk = {'-o','-s','-^','-d'};
cl = lines(nCase);
for c = 1:nCase
    plot(0:p.Nh-1, S(c).tau*1e9, mk{c}, 'Color', cl(c,:), 'LineWidth',1.6, ...
         'MarkerFaceColor', cl(c,:), 'MarkerSize',6);
end
xlabel('Delay element index'); ylabel('Delay (ns)');
title('C17 Fig.8 reproduction: delay per element');
legend(p.caseName, 'Location','northwest');
grid on; box on;
print(fg, fullfile(out,'b1_fig8_delay.png'), '-dpng', '-r180');
close(fg);

fprintf(fid,'max delay (ns) vs paper:\n');
fprintf(fid,'  %-8s %10s %10s\n','case','mine','paper');
paper = [2.5 7.5 17.5 35];
for c = 1:nCase
    fprintf(fid,'  %-8s %10.2f %10.2f\n', p.caseName{c}, Dmax(c), paper(c));
end
fclose(fid);

fprintf('B1 done: max delays %s ns\n', mat2str(round(Dmax,2)));
save(fullfile(out,'b1_designs.mat'),'S','Dmax');
