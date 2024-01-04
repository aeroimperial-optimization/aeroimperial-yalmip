function shift = computePSDShift(K)

% Compute shift for indexing in the PSD cone

shift = 0;
if isfield(K,'f') && ~isempty(K.f)
    shift = shift + K.f;
end
if isfield(K,'l') && ~isempty(K.l)
    shift = shift + K.l;
end
if isfield(K,'q') && ~isempty(K.l)
    shift = shift + K.l;
end
if isfield(K,'r') && ~isempty(K.l)
    shift = shift + K.l;
end
end