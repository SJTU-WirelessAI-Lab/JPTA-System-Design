function b = jp_bands(p, N_UE)
%JP_BANDS  Subband -> UE mapping used in Sec. IV.
%
%   The RB set is split into N_UE equal contiguous blocks; block i is served
%   by the i-th UE (non-overlapping RBs => no inter-user interference, as the
%   paper states).  Block edges are rounded to whole subcarriers.

edges = round(linspace(0, p.Nsc, N_UE+1));
b.idx = zeros(p.Nsc,1);
for i = 1:N_UE
    b.idx(edges(i)+1 : edges(i+1)) = i;
end
b.edges = edges;
b.rbOf  = @(i) floor((edges(i+1)-edges(i))/12);
end
