function tab = jp_c17_link_table(p, gainDb, nrbMax, T)
%JP_C17_LINK_TABLE  Pre-compute the effective SNR of every (RB, MCS) pair.
%
%   EESM is shift invariant in dB: adding a common dB offset to all
%   subcarrier SNRs shifts the effective SNR by exactly the same amount.
%   The path gain can therefore be applied afterwards, which turns the
%   distance sweep into a table lookup.

cand = [4 8 16 32 64 128 277];
cand = unique([cand(cand <= nrbMax), nrbMax]);

tab.cand = cand(:);
tab.eff0 = zeros(numel(cand), numel(T.mcs));   % effective SNR at gPath = 0
tab.k    = zeros(numel(cand), 1);              % subcarriers used

for ic = 1:numel(cand)
    m = cand(ic);
    k = min(12*m, numel(gainDb));
    g = gainDb(1:k);
    psd   = p.Ptx - 10*log10(k);               % fixed total UE power
    snr0  = psd + p.Gue + g(:).' - p.N0 - p.NF;
    lin   = 10.^(snr0/10);
    tab.k(ic) = k;
    for i = 1:numel(T.mcs)
        tab.eff0(ic,i) = jp_eesm(lin, T.eta(i));
    end
end
end
