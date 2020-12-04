function [At,c,K,momentID] = buildConstraintMatrix(num_vars,var_id,omega,numx,num_moments,clique_moments,inclique,g,At,c,K,mass)

% Construct the constraint matrix for an equality constraint with moments up to
% degree 2*omega in the variables with identifier specified by var_id. These are
% a subset of a larger number of variables (numx variables in total).
% We use a probabilistic search method for speed

% Construct exponents in the moment matrix for the given variables
degM = 2*omega-g.degree;
M = monolistcoeff(num_vars, degM, degM);

% Vectorize the exponents in a random way for faster assembly (it works with
% probability 1). Then, add constant moment (0)
hash = randn(numx,1);
moments_hash = clique_moments*hash;
moments_hash = [0; moments_hash];
inclique = [1; inclique+1];
Mhash = M*hash(var_id);

% Initialize empty row/column indices and values
rows = [];
cols = [];
vals = [];

% Loop over the coefficients of the multiplier g
num_eq = size(M,1);
for kk = 1:length(g.coef)
    % Multiply the monomial vector by the current monomial in g
    Qhash = Mhash + g.pows(kk,:)*hash(var_id);
    % Unwind local -> global mapping of variables
    [~,TEMP] = ismembertol(Qhash, moments_hash);
    % Indices and value for LMI in sedumi format
    rows = [rows; (1:num_eq).'];
    cols = [cols; inclique(TEMP)];
    vals = [vals; -g.coef(kk).*ones(size(TEMP,1),1)];
end

% Assemble the LMI in sedumi format
K.f = K.f + num_eq;
isy0 = (cols==1);
cols(~isy0) = cols(~isy0) - 1; % we remove the constant moment!
At.f{end+1} = sparse(rows(~isy0), cols(~isy0), vals(~isy0), num_eq, num_moments);
c.f{end+1} = sparse(rows(isy0), 1, -mass.*vals(isy0), num_eq, 1);

% Find moment indices in terms of ALL variables -- not just those involved
% in the moment matrix. NOTE: This list does NOT include the moment y0,
% whose value is just the mass of the underlying measures (1 by default)
if nargout > 3
    momentID = unique(cols(~isy0));
end