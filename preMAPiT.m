function [pre,tau,itau] = preMAPiT(s,pdf_s,x,pdf_x)

%% MAPiT measure preserving transformation of pseudotime into true time
%
%
%% INPUTS:
% s     - 1xn, pseudotime scale
% pdf_s - 1xn, marginal density on pseudotime scale
% x     - 1xn, true scale
% p_x   - 1xn, marginal density on true scale
%
%% OUTPUTS:
%
% pre   - 1xn, prefactor for transformations of probability distributions
% tau   - 1xn, transformation: pseudotime -> true scale
% itau  - 1xn, transformation: true scale -> pseudotime
%
%% Reference:
% 
% Reconstructing temporal and spatial dynamics from snap-shot data
% of heterogeneous cell populations
% Karsten Kuritz, Daniela Stöhr, Daniela Maichl, Nadine Pollak, Markus Rehm,
% Frank Allgöwer (2019)
%
%% ------------------------------------------------------


%% Parse inputs
% Check for positivity of pdf

if any(pdf_s<0) || any(pdf_x<0)
	error('pdf_s and pdf_x must be positive');
end

%% prepare transformation

% cumulative distribution
cdf_s	= cumtrapz(s,pdf_s);	
cdf_x	= cumtrapz(x,pdf_x);

% forced to cdf(end) = 1
cdf_s = cdf_s./cdf_s(end);
cdf_x = cdf_x./cdf_x(end);



% inverse cdf
tau		= interp1(cdf_x,x,cdf_s);		% x(ss)
itau	= interp1(cdf_s,s,cdf_x);		% s(xx)

% prefactor
pre		= pdf_x ./ interp1(s,pdf_s,itau);


% test_pdf_x = pre .* interp1(s,pdf_s,itau);
% plot(x,pdf_x,'r'); hold on
% plot(x,test_pdf_x,'b--')





