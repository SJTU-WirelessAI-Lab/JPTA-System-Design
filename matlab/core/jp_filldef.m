function opt = jp_filldef(opt, def)
% Merge a user-supplied option struct with defaults.
if isempty(opt), opt = struct(); end
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opt, fn{k})
        opt.(fn{k}) = def.(fn{k});
    end
end
end
