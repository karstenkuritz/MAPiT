function [p_xy, pdf_x] = MAPiT(s, y, p_sy, itau, pre)

%% Apply MAPiT to joint density
% 
% Applying the transformation obtained by preMAPiT to joint density of
% pseudotime and markers values.
%
%% INPUTS
%
% s     - 1xn pseudotime scale
% y     - 1xn marker scale
% p_sy  - nxn joint density p(PT,Y)
% itau  - 1xn, marginal density on pseudotime scale
% pre   - prefactor 
%
%% OUTPUTS
%
% p_xy  - nxn joint density with true scale p(X,Y)
% p_x   - 1xn marginal density on true scale (for validation)
%
%% Reference:
% 
% Reconstructing temporal and spatial dynamics from snap-shot data
% of heterogeneous cell populations
% Karsten Kuritz, Daniela Stöhr, Daniela Maichl, Nadine Pollak, Markus Rehm,
% Frank Allgöwer (2019)
%
%% ------------------------------------------------------

% apply transformation to joint density
p_xy = bsxfun(@times,interp2(s,y,p_sy,itau,y),pre); 

% marginal density on true scale
pdf_x = trapz(y,p_xy,1);

