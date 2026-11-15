function [thr, out] = jp_c17_link_eval(p, tab, T, gPath, duty)
%JP_C17_LINK_EVAL  Best (RB, MCS) throughput for one UE at a given path gain.
%
%   Implements the scheduler of C17 Sec. IV-A:
%       max_{m, MCS} m*12*SCS*SE*(1-BLER)   s.t. BLER <= 10%, m >= 4,
%   with the UE's fixed total power spread over the m RBs it is granted.

effDb = tab.eff0 + gPath + 0;                      % shift by the path gain
bler  = 0.5*erfc((effDb - T.snrTh.' + T.blerOff*T.blerSigma) / ...
                 (sqrt(2)*T.blerSigma));
ok    = bler <= p.blerTgt;
gain  = (tab.cand*12*p.SCS/1e6) * (T.se.') .* (1-bler) .* ok;
[thr, idx] = max(gain(:));
[ic, im]   = ind2sub(size(gain), idx);

thr = thr*duty*p.ulDuty;
out.nrb = tab.cand(ic);  out.mcs = T.mcs(im);  out.se = T.se(im);
out.bler = bler(ic,im);  out.effDb = effDb(ic,im);
end
