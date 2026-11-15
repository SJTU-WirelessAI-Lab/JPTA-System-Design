% =======================================================================
%  run_b1c.m  --  Fig. 7 redone with the paper's OWN Fig. 8 delay profile
%
%  Why this script exists: b1 designed the beam by solving eq. (4) from
%  scratch, which produced a delay profile of its own.  But Sec. IV-B says
%  Fig. 7 ("the corresponding beam patterns") and Fig. 8 ("the corresponding
%  configuration of the delay elements") describe the SAME configuration, so
%  Fig. 7 must be evaluated with Fig. 8's delays.  Here the delays are fixed
%  a priori by eq. (5) and only the phase shifter is optimised (the closed
%  form of eq. (4) in phi), which is exactly the hardware situation.
%
%  Three delay rules:
%    'eq5q'  eq. (5) step sin(dTheta/2)/W, rounded UP to the 2.5 ns TTD grid,
%            wrapped every N_UE elements  -> this IS Fig. 8 (see 8.2')
%    'eq5'   eq. (5) exact, continuous delays (what the TTD would need)
%    'ideal' p_k = b_k, one ideal steering vector per subband (upper bound)
%
%  Metrics per UE: worst / mean in-band gain at the UE's own angle over its
%  own subband, plus the ripple = mean - worst.  Paper: "can decrease up to
%  3 dB compared to the maximum beam gain of 28 dB".
% =======================================================================
addpath(fullfile(fileparts(mfilename('fullpath')),'core'));
addpath(fullfile(fileparts(mfilename('fullpath')),'alg'));
p   = jp_params();
out = fullfile(fileparts(mfilename('fullpath')), 'out');
if ~exist(out,'dir'), mkdir(out); end
fid = fopen(fullfile(out,'b1c_summary.txt'),'w');
angV = -60:0.25:60;

fprintf(fid,'Fig.7 with the Fig.8 delay profile (per-element delays fixed, phi = eq.(4) closed form)\n');
fprintf(fid,'paper: in-band gain >= 28-3 = 25 dB for every UE\n\n');

MODE = {'eq5q','eq5'};
fprintf(fid,'%-8s %-6s %8s %8s %10s %10s %10s %8s %10s\n', ...
        'case','rule','tau_max','peak','worstUE','meanUE','ripple','envelope','verdict');
res = struct();
for c = 1:numel(p.caseAngles)
    th = p.caseAngles{c};  Nu = numel(th);
    b  = jp_bands(p, Nu);
    Ws = (b.edges(2)-b.edges(1))*p.SCS;            % subband width [Hz]

    for k = 1:numel(MODE)
        [phi, tau] = design_fixed_step(p, th, b.idx, MODE{k});
        gPk = zeros(1,Nu); gW = zeros(1,Nu); gM = zeros(1,Nu);
        for u = 1:Nu
            g = jp_pattern(p, phi, tau, th(u), p.f(b.idx==u));
            g = g(:).';
            gW(u) = min(g); gM(u) = mean(g); gPk(u) = max(g);
        end
        g  = jp_pattern(p, phi, tau, angV, p.f);
        pk = max(g(:));
        env = min(max(g,[],1));            % ridge height vs frequency
        worst = min(gW); meanW = mean(gW); rip = max(gM - gW);
        ok = all(gW >= p.Gbs-3);
        fprintf(fid,'%-8s %-6s %8.2f %8.2f %10.2f %10.2f %10.2f %8.2f %10s\n', ...
            p.caseName{c}, MODE{k}, max(tau)*1e9, pk, worst, meanW, rip, env, ...
            ternary(ok,'PASS','FAIL'));
        if k == 1
            fprintf(fid,'         per-UE worst in-band gain (dB): %s\n', mat2str(round(gW,2)));
            fprintf(fid,'         delay (ns): %s\n', mat2str(round(tau.'*1e9,2)));
            fprintf(fid,'         subband width %.0f MHz\n', Ws/1e6);
        end
        res(c).(MODE{k}) = struct('phi',phi,'tau',tau,'gW',gW,'gM',gM,'pk',pk,'G',g);
    end
    % ideal per-subband steering: p_k = b_k  =>  gain is exactly the peak
    fprintf(fid,'%-8s %-6s %8s %8.2f %10.2f %10.2f %10.2f %8.2f %10s\n', ...
        p.caseName{c}, 'ideal', '-', p.Gbs, p.Gbs, p.Gbs, 0, p.Gbs, 'PASS');
    fprintf(fid,'\n');
end

% ---- how fine the TTD must be: answered empirically in run_b1d -----------
% No closed-form number is quoted here. The steering phase is cos(theta), so the
% 3 dB width of the main lobe depends on the user angle itself, and a single
% "required resolution" would be meaningless for a multi-user configuration.
fprintf(fid,'required TTD resolution: see out/b1d_summary.txt (step-size sweep).\n');
fclose(fid);

% ---- Fig.7 panels from the Fig.8 delay profile -------------------------
for k = 1:numel(MODE)
    fg = figure('Visible','off','Position',[80 80 1240 880]);
    for c = 1:numel(p.caseAngles)
        th = p.caseAngles{c}; Nu = numel(th); b = jp_bands(p, Nu);
        subplot(2,2,c);
        imagesc(p.foff/1e6, angV, res(c).(MODE{k}).G); axis xy; hold on;
        for i = 1:Nu
            plot([p.foff(b.edges(i)+1) p.foff(b.edges(i+1))]/1e6, [th(i) th(i)], ...
                 'w--','LineWidth',1.1);
            plot([p.foff(b.edges(i)+1) p.foff(b.edges(i+1))]/1e6, [-th(i) -th(i)], ...
                 'w:','LineWidth',0.9);
        end
        caxis([10 28]); colorbar;
        xlabel('f - f_c (MHz)'); ylabel('angle (deg)');
        ttl = sprintf('%d UEs, tau_{max}=%.1f ns', Nu, max(res(c).(MODE{k}).tau)*1e9);
        title(ttl); ylim([-60 60]);
    end
    sgtitle(sprintf('Fig.7 reproduction, delay rule = %s', MODE{k}));
    print(fg, fullfile(out,sprintf('b1c_fig7_%s.png',MODE{k})), '-dpng','-r150'); close(fg);
end

% ---- in-band gain vs frequency, 4-UE case ------------------------------
c = 2; th = p.caseAngles{c}; Nu = numel(th); b = jp_bands(p, Nu);
fg = figure('Visible','off','Position',[80 80 1100 460]);
for u = 1:Nu
    subplot(1,Nu,u);
    hold on;
    for k = 1:numel(MODE)
        g = jp_pattern(p, res(c).(MODE{k}).phi, res(c).(MODE{k}).tau, th(u), p.f(b.idx==u));
        plot(p.foff(b.idx==u)/1e6, g(:), 'LineWidth',1.6);
    end
    yline(p.Gbs,'k--'); yline(p.Gbs-3,'r:');
    ylim([0 29]); grid on; box on;
    xlabel('f - f_c (MHz)'); if u==1, ylabel('gain at own angle (dB)'); end
    title(sprintf('UE %d (%.0f%c)', u, th(u), 176));
end
legend({'eq.(5) rounded to 2.5 ns','eq.(5) exact','peak 28 dB','28-3 dB'}, ...
       'Location','southoutside','Orientation','horizontal');
sgtitle('4-UE case: in-band gain over the UE''s own subband');
print(fg, fullfile(out,'b1c_inband.png'), '-dpng','-r150'); close(fg);

type(fullfile(out,'b1c_summary.txt'));

% =======================================================================
function [phi,tau] = design_fixed_step(p, th, bandIdx, mode)
M  = p.Nh; k0 = (0:M-1).'; Nu = numel(th);
Ct = reshape(cosd(th(:).'),[],1); Ct = reshape(Ct(bandIdx),[],1);
B  = exp(1j*2*pi*(p.d/p.c)*((p.f.*Ct)*k0.'));
switch mode
    case 'eq5q'                       % eq.(5) step on the 2.5 ns grid, wrapped
        dth = max(th)-min(th);
        st  = ceil(sind(dth/2)/p.BW/p.tauStep)*p.tauStep;
        if st == 0, st = p.tauStep; end
        tau = st*mod(k0, Nu);
    case 'eq5'                        % eq.(5) exact, continuous
        dth = max(th)-min(th);
        tau = mod(k0, Nu)*(sind(dth/2)/p.BW);
    otherwise
        error('bad mode');
end
phi = zeros(M,1);
for m = 1:M
    A = sum(B(:,m).*exp(1j*2*pi*p.f*tau(m)));
    phi(m) = angle(A);
end
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end
