% =======================================================================
%  run_b2.m  --  Fig. 9: mean UE throughput vs distance (PAA vs JPTA)
%  System-level simulation per Sec. IV-A with the the uplink-coverage paper link abstraction.
% =======================================================================
addpath(fullfile(fileparts(mfilename('fullpath')),'core'));
addpath(fullfile(fileparts(mfilename('fullpath')),'alg'));
p   = jp_params();
T   = jp_mcs();
out = fullfile(fileparts(mfilename('fullpath')), 'out');
if ~exist(out,'dir'), mkdir(out); end

D    = (30:30:1500).';
nUE  = [2 4 8 16];
nA   = p.Npaa;
bAng = linspace(-55,55,nA);            % PAA codebook [assumed uniform]
gEdge = 3;                             % JPTA <=3 dB below peak (paper)

Gh = 20*log10(p.c/(4*pi*p.fc)) - 10*p.beta*log10(D);

R_paa = zeros(numel(D),numel(nUE));
R_jpt = zeros(numel(D),numel(nUE));
Dcov_p = zeros(1,numel(nUE));  Dcov_j = zeros(1,numel(nUE));
gPAA = cell(1,numel(nUE));  rbJU = zeros(1,numel(nUE));

for iu = 1:numel(nUE)
    Nu = nUE(iu);
    th = p.caseAngles{iu};
    bw = p.BW/Nu;
    rbJ = floor(p.Nrb/Nu);
    rbJU(iu) = rbJ;

    % ---- PAA: best of the 16 frequency-flat beams ----------------------
    gp = zeros(1,Nu);
    for u = 1:Nu
        a = exp(-1j*2*pi*(p.d/p.c)*((p.f*cosd(th(u)))*(0:p.Nh-1)));
        gbest = -inf;
        for b = 1:nA
            wm = exp(1j*2*pi*(p.d/p.c)*p.fc*cosd(bAng(b))*(0:p.Nh-1).');
            af = abs(a*wm).^2/p.Nh;
            gbest = max(gbest, max(10*log10(max(af,eps)) + 10*log10(p.Nv) + p.Gelem));
        end
        gp(u) = gbest;
    end
    gPAA{iu} = gp;

    % ---- link tables (EESM is shift invariant in dB) -------------------
    tabP = cell(1,Nu);  tabJ = cell(1,Nu);
    for u = 1:Nu
        tabP{u} = jp_link_table(p, gp(u)*ones(p.Nsc,1), p.Nrb, T);
        sel = ((1:p.Nsc).' > (u-1)*p.Nsc/Nu) & ((1:p.Nsc).' <= u*p.Nsc/Nu);
        fo  = p.foff(sel);
        fo  = fo - mean(fo);                 % offset from the subband centre
        gj  = p.Gbs - gEdge*(2*fo/bw).^2;    % 28 dB centre, >=25 dB at edges
        tabJ{u} = jp_link_table(p, gj, rbJ, T);
    end

    for id = 1:numel(D)
        ap = 0;  aj = 0;
        for u = 1:Nu
            ap = ap + jp_link_eval(p, tabP{u}, T, Gh(id), 1/Nu);
            aj = aj + jp_link_eval(p, tabJ{u}, T, Gh(id), 1);
        end
        R_paa(id,iu) = ap/Nu;
        R_jpt(id,iu) = aj/Nu;
    end
    fprintf('  case %s done: PAA@30m=%.1f JPTA@30m=%.1f\n', p.caseName{iu}, ...
            R_paa(1,iu), R_jpt(1,iu));
    Dcov_p(iu) = coverDist(D, R_paa(:,iu), 1.0);
    Dcov_j(iu) = coverDist(D, R_jpt(:,iu), 1.0);
end

% ---- report -------------------------------------------------------------
fid = fopen(fullfile(out,'b2_summary.txt'),'w');
fprintf(fid,'Fig.9 reproduction | Table-I link budget | beta = %g\n', p.beta);
fprintf(fid,'PAA: all RBs, 1-of-N slots | JPTA: 1/N of RBs, every slot\n');
fprintf(fid,'TDD uplink duty = %.2f, BLER target = %.0f%%\n\n', p.ulDuty, 100*p.blerTgt);
for iu = 1:numel(nUE)
    fprintf(fid,'=== %s  (PAA gain %s dB, JPTA RBs/UE = %d)\n', ...
            p.caseName{iu}, mat2str(round(gPAA{iu},1)), rbJU(iu));
    fprintf(fid,'   dist     PAA     JPTA   ratio\n');
    idx = unique([1:5:numel(D), numel(D)]);
    for id = idx
        fprintf(fid,'  %5.0f  %7.4f %7.4f  %5.2fx\n', D(id), ...
                R_paa(id,iu), R_jpt(id,iu), R_jpt(id,iu)/max(R_paa(id,iu),eps));
    end
    fprintf(fid,'  1 Mbps coverage: PAA %.0f m -> JPTA %.0f m (x%.2f)\n\n', ...
            Dcov_p(iu), Dcov_j(iu), Dcov_j(iu)/max(Dcov_p(iu),eps));
end
fprintf(fid,'ratio at 1500 m: ');
for iu = 1:numel(nUE)
    fprintf(fid,'%s %.1fx (+%.0f%%)  ', p.caseName{iu}, ...
            R_jpt(end,iu)/max(R_paa(end,iu),eps), ...
            100*(R_jpt(end,iu)/max(R_paa(end,iu),eps)-1));
end
fprintf(fid,'\n');
fclose(fid);

csv = fopen(fullfile(out,'b2_throughput.csv'),'w');
fprintf(csv,'dist');
for iu = 1:numel(nUE), fprintf(csv,',paa%d,jpta%d',nUE(iu),nUE(iu)); end
fprintf(csv,'\n');
for id = 1:numel(D)
    fprintf(csv,'%.0f',D(id));
    for iu = 1:numel(nUE), fprintf(csv,',%.5f,%.5f',R_paa(id,iu),R_jpt(id,iu)); end
    fprintf(csv,'\n');
end
fclose(csv);

fg = figure('Visible','off','Position',[100 100 660 480]); hold on;
cl = lines(numel(nUE));
for iu = 1:numel(nUE)
    plot(D, R_paa(:,iu), '-',  'Color',cl(iu,:), 'LineWidth',2.2);
    plot(D, R_jpt(:,iu), '--', 'Color',cl(iu,:), 'LineWidth',2.2);
end
set(gca,'YScale','log'); grid on; box on; ylim([1e-2 1e2]);
xlabel('2D Distance (m)'); ylabel('Mean Throughput (Mbps)');
lg = cell(2*numel(nUE),1);
for iu = 1:numel(nUE)
    lg{2*iu-1} = sprintf('PAA %d UEs', nUE(iu));
    lg{2*iu}   = sprintf('JPTA %d UEs', nUE(iu));
end
legend(lg,'Location','northeast'); title('Fig.9 reproduction');
print(fg, fullfile(out,'b2_fig9_throughput.png'), '-dpng', '-r180');
close(fg);

type(fullfile(out,'b2_summary.txt'));
save(fullfile(out,'b2_sls.mat'),'D','R_paa','R_jpt','Dcov_p','Dcov_j','nUE','gPAA');

% ------------------------------------------------------------------------
function dc = coverDist(D, R, thr)
ok = find(R >= thr);
if isempty(ok), dc = 0; return; end
i = ok(end);
if i == numel(D), dc = D(end); return; end
x = [D(i) D(i+1)];
y = log10(max(R(i:i+1),1e-6));
dc = x(1) + (log10(thr)-y(1))/(y(2)-y(1))*(x(2)-x(1));
end
