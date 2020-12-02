function [At,c,K,momentID,gramMonomials] = buildMomentMatrix(num_vars,var_id,omega,numx,all_moments,g,At,c,K,mass)

% Construct the moment matrix with moments up to degree 2*omega in the variables
% with identifier specified by var_id. These are a subset of a larger number of
% variables (numx variables in total).
% Use functions from YALMIP to compute the outer product of the monomial vector
% We also use a probabilistic search method for speed

% Is weight g for localizing matrix empty? If so, set it to 1 (in
% decomposed form)
if isempty(g)
    g.coef = 1;
    g.pows = zeros(1,num_vars);
    g.degree = 0;
end

% Construct exponents in the moment matrix for the given variables
degM = omega - ceil(0.5*g.degree);
MMt = monolistcoeff(num_vars, degM, degM);
if nargout > 4
    gramMonomials = MMt;
end
nsdp = size(MMt,1);
[MMt,unique_moments] = monomialproducts({MMt});
MMt = MMt{1};

% Vectorize the exponents in a random way for faster assembly (it works with
% probability 1)
hash = rand(numx,1);
moments_hash = all_moments*hash;
num_moments = length(moments_hash);

% Add constant moment
moments_hash = [0; moments_hash];

% Initialize empty row/column indices and values
rows = [];
cols = [];
vals = [];

% Loop over the coefficients of the multiplier g
for kk = 1:length(g.coef)
    % Multiply the PSD matrix obtained by outer product of monomials by
    % the current monomial in g
    Q = MMt + repmat([0, 0, g.pows(kk,:)], size(MMt,1), 1);
    P = unique_moments + repmat([0, 0, g.pows(kk,:)], size(unique_moments,1), 1);
    Phash = P(:,3:end)*hash(var_id);
    Qhash = Q(:,3:end)*hash(var_id);
    % Unwind local -> global mapping of variables
    [~,LOCB1] = ismembertol(Qhash, Phash);
    [~,LOCB2] = ismembertol(Phash, moments_hash);
    TEMP = LOCB2(LOCB1);
    % Indices and value for LMI in sedumi format
    rows = [rows; Q(:,1)+(Q(:,2)-1).*nsdp];
    cols = [cols; TEMP];
    vals = [vals; -g.coef(kk).*ones(size(TEMP,1),1)];
end

% Assemble the LMI in sedumi format
K.s(end+1) = nsdp;
isy0 = (cols==1);
cols(~isy0) = cols(~isy0) - 1; % we remove the constant moment!
At.s{end+1} = sparse(rows(~isy0), cols(~isy0), vals(~isy0), nsdp^2, num_moments);
c.s{end+1} = sparse(rows(isy0), 1, -mass.*vals(isy0), nsdp^2, 1);

% Find moment indices in terms of ALL variables -- not just those involved
% in the moment matrix. NOTE: This list does NOT include the moment y0,
% whose value is just the mass of the underlying measures (1 by default)
if nargout > 3
    momentID = unique(cols(~isy0));
end