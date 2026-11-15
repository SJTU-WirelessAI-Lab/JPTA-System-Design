function p = jp_params()
%JP_PARAMS  Parameters from (AlAmmouri/Mo et al., Asilomar 2024), Table I.
%
%   Everything the paper states is set here; every quantity the paper does
%   NOT state is marked [assumed] with the value used (see REPORT.md G16-G19).

p.c     = 299792458;          % speed of light [m/s]
p.fc    = 28e9;               % carrier frequency [Hz]              (Table I)
p.BW    = 400e6;              % system bandwidth [Hz]               (Table I)

% ---- OFDM grid ------------------------------------------------------------
p.Nfft  = 4096;               % FFT size                            (Table I)
p.SCS   = 120e3;              % subcarrier spacing [Hz]             (Table I)
p.Nsc   = floor(p.BW/p.SCS);  % occupied subcarriers = 3333
p.Nrb   = floor(p.Nsc/12);    % resource blocks      = 277
p.rbBW  = 12*p.SCS;           % 1.44 MHz per RB

% subcarrier offset from the carrier, symmetric about fc
m       = (0:p.Nsc-1).';
p.f     = p.fc + (m - (p.Nsc-1)/2)*p.SCS;      % absolute frequency [Hz]
p.foff  = p.f - p.fc;                          % offset [Hz]

% ---- array ---------------------------------------------------------------
p.Nv    = 16;                 % vertical elements                   (Sec. IV-A)
p.Nh    = 16;                 % horizontal elements = #delay units  (Sec. IV-A)
p.d     = 0.5*p.c/p.fc;       % half-wavelength spacing [assumed]
p.Gelem = 28 - 10*log10(p.Nv*p.Nh);   % 3.92 dB element gain, so that the
                                      % 16x16 peak gain equals 28 dB (Table I)

% ---- link budget ---------------------------------------------------------
p.Ptx   = 23;                 % UE TX power [dBm]                   (Table I)
p.Gue   = 0;                  % UE beam gain [dB]                   (Table I)
p.Gbs   = 28;                 % BS peak beam gain [dB]              (Table I)
p.NF    = 5;                  % BS noise figure [dB]                (Table I)
p.beta  = 3;                  % path-loss exponent                  (Table I)
p.N0    = -174 + 10*log10(p.SCS);     % noise PSD per subcarrier [dBm]

% ---- delay units ---------------------------------------------------------
p.tauStep = 2.5e-9;           % TTD step size                       (Table I)
p.tauMax  = 200e-9;           % [assumed] top of the search window; the
                              % paper only shows that <=35 ns is needed
p.nPS     = Inf;              % [assumed] PS resolution: continuous

% ---- scheduling ----------------------------------------------------------
p.Npaa    = 16;               % analog beams in the PAA             (Table I)
p.minRB   = 4;                % minimum RB allocation               (Sec. IV-A)
p.blerTgt = 0.10;             % BLER constraint                     (Sec. IV-A)
p.ulDuty  = 0.2;              % [assumed] TDD uplink duty cycle (DDDSU).  The
                              % paper quotes absolute Mbps but never states the
                              % TDD pattern; 0.2 is the only free scale left
                              % and is what aligns the plateau of Fig. 9.
p.gainThr = 3;                % [assumed] subcarriers >3 dB below the
                              % peak are dropped before EESM (eq. 11 of the uplink-coverage paper)

% ---- scenario angles -----------------------------------------------------
p.caseAngles = {[-30 30], ...
                [-30 -10 10 30], ...
                -55 + 110*(0:7)/7, ...
                -55 + 110*(0:15)/15};
p.caseName   = {'2 UEs','4 UEs','8 UEs','16 UEs'};
end
