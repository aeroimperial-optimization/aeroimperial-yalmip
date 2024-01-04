function options = admmOpts

% ADMMOPTS
%
% Default options for ADMM solver
%
% Generic solver options
% ----------------------
% options.relTol     = 1e-4;      % tolerance
% options.verbose    = 1;         % print or silent
% options.dispIter   = 50;        % print every dispIter iterations
% options.maxIter    = 1000;      % max # iterations
%
% ADMM options
% --------------------
% options.lambda   = 0.5;       % splitting parameter for the cost function
% options.rho      = 1;         % penalty parameter
% options.adaptive = true;      % adaptive penalty factor?
% options.tau      = 2;         % increase factor for adaptive penalty scheme (must be > 1)
% options.mu       = 10;        % ratio of residuals for adaptive penalty scheme
% options.rhoMax   = 1e6;       % maximum penalty parameter
% options.rhoMin   = 1e-6;      % minimum penalty parameter
% options.rhoIt    = 10;        % update options.rho every options.rhoIt iterations
%
% Advanced options
% ----------------
% options.KKTfact    = 'blk';     % Options for KKT systems
%                                 %  a) 'blk': block elimination,
%                                 %  b) 'ldl': ldl factor,
%                                 %  c) 'inv': invert

% Create options structure
options.adaptive = true;     % adaptive penalty factor?
options.dispIter = 50;        % print every dispIter iterations
options.KKTfact  = 'blk';     % Options for KKT systems
options.lambda   = 0.5;       % add penalty term for Y block to cost
options.maxIter  = 1e4;       % max # iterations
options.mu       = 2;         % increase/decrease factor for adaptive penalty scheme (must be > 1)
options.nu       = 10;        % ratio of residuals for adaptive penalty scheme
options.relTol   = 1e-3;      % tolerance
options.rescale  = 0;         % try to rescale  data (does not work well...)
options.rho      = 1e3;        % penalty parameter
options.rhoIt    = 100;       % if pres/dres>mu (<mu) mu for rhoIt iterations, adapt rho
options.rhoMax   = 1e6;       % maximum penalty parameter
options.rhoMin   = 1e-6;      % minimum penalty parameter
options.verbose  = 1;         % print or silent

