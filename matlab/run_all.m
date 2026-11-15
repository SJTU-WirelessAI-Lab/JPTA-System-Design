function run_all(what)
%RUN_ALL  Single entry point for every experiment in this repository.
%
%   run_all              run every figure in the paper order
%   run_all('fig7raw')   Fig. 7 from the unconstrained optimum of eq. (4)   ~10 s
%   run_all('design')    Fig. 7 for the three design rules of Sec. IV-B      ~1 s
%   run_all('delaylaw')  Fig. 8: delay law, quantization, wrap-around       ~10 s
%   run_all('fig7')      Fig. 7 with the Fig. 8 delays, phase-only design   ~40 s
%   run_all('step')      TTD step-size sweep (how fine must the TTD be)
%   run_all('ramp')      straight ramp versus sawtooth
%   run_all('fig9')      Fig. 9: throughput and coverage
%
%   Every job writes its figures and tables into matlab/out/.

if nargin < 1 || isempty(what), what = 'all'; end
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here,'core'));
addpath(fullfile(here,'alg'));

% key, script, kind ('f' = function file, 's' = script file), description
jobs = {
    'fig7raw', 'run_b1',  's', 'Fig. 7 from the unconstrained global optimum of eq. (4)'
    'design',  'run_b1b', 's', 'Fig. 7 for the three candidate design rules of Sec. IV-B'
    'delaylaw','run_b3',  's', 'Fig. 8: delay law of eq. (5), 2.5 ns quantization, wrap-around'
    'fig7',    'run_b1c', 's', 'Fig. 7 with the Fig. 8 delays and the closed-form phase'
    'step',    'run_b1d', 's', 'TTD step-size sweep: how fine must the delay resolution be'
    'ramp',    'run_b1e', 's', 'straight ramp versus sawtooth delay profile'
    'fig9',    'run_b2',  's', 'Fig. 9: throughput and coverage, EESM + MCS/RB scheduling'
};

if strcmpi(what,'all')
    sel = true(size(jobs,1),1);
else
    sel = strcmpi(jobs(:,1), what);
    if ~any(sel)
        fprintf('known jobs: all, %s\n', strjoin(jobs(:,1).', ', '));
        return
    end
end

for k = find(sel).'
    fprintf('\n########## %s : %s ##########\n', jobs{k,1}, jobs{k,4});
    t0 = tic;
    switch jobs{k,3}
        case 'f', feval(jobs{k,2});
        case 's', evalin('base', sprintf('run(''%s'')', fullfile(here,[jobs{k,2} '.m'])));
        %  A script file runs in the workspace of its caller, so it would overwrite the
        %  loop index and the timer of this function.  Evaluating it in the base workspace
        %  keeps run_all itself intact; all paths are global, so nothing else is needed.
    end
    fprintf('########## %s done in %.1f s ##########\n', jobs{k,1}, toc(t0));
end
