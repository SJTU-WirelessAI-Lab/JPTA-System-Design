% =======================================================================
%  run_b1e.m  --  why Fig. 8 is a sawtooth and not eq. (5)'s straight ramp
%
%  eq. (5) gives a linear ramp over the whole array,
%      tau_m = (m-1) * sin(dTheta/2)/W,      m = 1..M          ("ramp")
%  Fig. 8 wraps the same slope every N_UE elements,
%      tau_m = ((m-1) mod N_UE) * sin(dTheta/2)/W            ("saw")
%  The slopes are identical; only the wrap differs.  A linear ramp is a
%  single frequency-independent beam (that is the whole point of TTD), so it
%  can NOT put N_UE beams in N_UE subbands.  The wrap tiles the aperture into
%  M/N_UE teeth, which is what turns one beam into a comb of N_UE beams.
% =======================================================================
addpath(fullfile(fileparts(mfilename('fullpath')),'core'));
addpath(fullfile(fileparts(mfilename('fullpath')),'alg'));
p   = jp_params_c17();
out = fullfile(fileparts(mfilename('fullpath')), 'out');
if ~exist(out,'dir'), mkdir(out); end
fid = fopen(fullfile(out,'b1e_summary.txt'),'w');
angV = -60:0.25:60;

fprintf(fid,'ramp vs sawtooth, same slope (eq. 5), M=%d elements, W=%.0f MHz\n\n', p.Nh, p.BW/1e6);
fprintf(fid,'%-8s %-5s %8s %10s %12s %10s\n','case','rule','tau_max','#beams<-3dB','worstUE(dB)','spread(dB)');
for c = 1:numel(p.caseAngles)
    th = p.caseAngles{c}; Nu = numel(th); b = jp_c17_bands(p, Nu);
    dth = max(th)-min(th); dt = sind(dth/2)/p.BW;
    rules = {'ramp','saw'};
    for k = 1:2
        switch rules{k}
            case 'ramp', tau = dt*(0:p.Nh-1).';
            case 'saw',  tau = dt*mod((0:p.Nh-1).', Nu);
        end
        phi = opt_phi(p, th, b.idx, tau);
        gW  = zeros(1,Nu);
        for u = 1:Nu
            gu = jp_c17_pattern(p, phi, tau, th(u), p.f(b.idx==u));
            gW(u) = min(gu(:));
        end
        % beams at the centre of the band: local maxima within 3 dB of the peak
        g  = jp_c17_pattern(p, phi, tau, angV, p.fc); g = g(:).';
        lm = find(g(2:end-1) > g(1:end-2) & g(2:end-1) > g(3:end)) + 1;
        nb = sum(g(lm) >= max(g)-3);
        fprintf(fid,'%-8s %-5s %8.2f %10d %12.2f %10.2f\n', p.caseName{c}, rules{k}, ...
                max(tau)*1e9, nb, min(gW), max(gW)-min(gW));

        if c == 2
            E.(rules{k}) = struct('phi',phi,'tau',tau, ...
                'G', jp_c17_pattern(p, phi, tau, angV, p.f));
        end
    end
    fprintf(fid,'\n');
end
fprintf(fid,'note: "beams" counts local maxima at f_c within 3 dB of the peak.\n');
fprintf(fid,'      the ramp forms one beam; only the wrap puts a beam in each subband.\n');
fclose(fid);

% ---- picture: ramp vs saw for the 4-UE case ----------------------------
fg = figure('Visible','off','Position',[80 80 1240 470]);
rules = {'ramp','saw'};
for k = 1:2
    subplot(1,2,k);
    imagesc(p.foff/1e6, angV, E.(rules{k}).G); axis xy; hold on;
    Nu = numel(p.caseAngles{2}); b = jp_c17_bands(p, Nu);
    for i = 1:Nu
        plot([p.foff(b.edges(i)+1) p.foff(b.edges(i+1))]/1e6, [p.caseAngles{2}(i) p.caseAngles{2}(i)], ...
             'w--','LineWidth',1.1);
    end
    caxis([0 28]); colorbar; xlabel('f-f_c (MHz)'); ylabel('angle (deg)');
    ylim([-60 60]);
    title(sprintf('%s: tau_{max}=%.2f ns, %d beams', rules{k}, ...
          max(E.(rules{k}).tau)*1e9, numel(find_peaks(E.(rules{k}).G, angV))));
end
sgtitle('4 UEs: eq.(5) as a straight ramp (left) vs wrapped every N_{UE} (right)');
print(fg, fullfile(out,'b1e_ramp_vs_saw.png'), '-dpng','-r150'); close(fg);
type(fullfile(out,'b1e_summary.txt'));

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

function pk = find_peaks(G, angV)
g  = G(:,round(end/2)); g = g(:).';
lm = find(g(2:end-1) > g(1:end-2) & g(2:end-1) > g(3:end)) + 1;
pk = angV(lm(g(lm) >= max(g)-3));
end
