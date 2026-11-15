% =======================================================================
%  run_b1b.m  --  compare the three design rules for C17 Type-1 codebooks
%  Diagnostic step: which rule reproduces Fig. 7 (pattern) and Fig. 8 (delay)?
% =======================================================================
addpath(fullfile(fileparts(mfilename('fullpath')),'core'));
addpath(fullfile(fileparts(mfilename('fullpath')),'alg'));
p   = jp_params_c17();
out = fullfile(fileparts(mfilename('fullpath')), 'out');
if ~exist(out,'dir'), mkdir(out); end

modes = {'ls','sawtooth','rainbow'};
angF  = -60:0.5:60;
fid = fopen(fullfile(out,'b1b_design_compare.txt'),'w');

for c = 1:numel(p.caseAngles)
    th = p.caseAngles{c};  Nu = numel(th);
    b  = jp_c17_bands(p, Nu);
    fprintf(fid,'\n===== %s, angles %s\n', p.caseName{c}, mat2str(round(th,1)));
    for k = 1:numel(modes)
        [phi, tau, info] = jp_c17_design(p, th, b.idx, modes{k});

        % peak angle at the centre frequency of each UE's subband
        pk = zeros(1,Nu);  gn = zeros(1,Nu);
        for u = 1:Nu
            sel = find(b.idx==u);
            fc  = p.f(sel(round(end/2)));
            g   = jp_c17_pattern(p, phi, tau, angF, fc);
            [gn(u), i] = max(g);
            pk(u) = angF(i);
        end
        fprintf(fid,'  %-9s tau_max=%6.2f ns | peak angle %s (want %s)\n', ...
                modes{k}, max(tau)*1e9, mat2str(pk), mat2str(round(th)));
        fprintf(fid,'            peak gain  %s dB | worst in-band %s dB\n', ...
                mat2str(round(gn,1)), mat2str(round(info.gIn,1)));
        fprintf(fid,'            tau = %s\n', mat2str(round(tau.'*1e9,1)));
    end
end
fclose(fid);
type(fullfile(out,'b1b_design_compare.txt'));
