function [At, b, c, K,  exponents] = makeConicSplittingVariable(AA, b, c, K, exponents, shift, whichClique, clique_all_exponents)

% Build the constraint matrix At. Keep clique moments individual, and use a
% splitting variable to enforce consistency. This helps with first-order
% splitting methods (parallelization)?
% Use sparse indexing for speed. Order of the variables:
% [y_matching; y_1; y_2; ...; y_numCliques]
% input "shift" gives last entry of "exponents" corresponding to moments in
% the objective, so can easily keep track of local vs global variables when
% looping through tows of "exponents"

error('Not implemented yet. Need to think about it -- need better code structure?')

% End function
end
