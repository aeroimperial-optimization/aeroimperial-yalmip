function [sol, opts] = postprocess(sol, opts)

% Postprocess solution: rescaling
if opts.rescale
   sol.y = opts.scaleFactors.y .* sol.y;
   for i = 1:opts.noCliques
      sol.s{i} = opts.scaleFactors.s{i} .* sol.s{i};
      sol.z{i} = opts.scaleFactors.eta{i} .* sol.z{i};
      sol.eta{i} = opts.scaleFactors.eta{i} .* sol.eta{i};
      sol.xi{i} = opts.scaleFactors.xi{i} .* sol.xi{i};
   end
end

end