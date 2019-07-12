function [s,y,pdf_s,p_sy,bw] = jointDensityPseudotimeY(PT,Y,t_scale)

%% Joint and marginal density in Pseudotime and marker
% This function provides the joint and marginal density of marker and
% pseudotime values using kernel density estimation with diffusion kde2d 
% (https://de.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation)
%
%% INPUTS
% PT    - nx1 Pseudotime values
% Y     - nx1 marker intensities 
% t_scale - 1x1 factor for adjusting the bandwidth of KDE (default = 1)
%
%% OUTPUTS
%
% s     - linspace(0,1,256), pseudotime scale
% y     - 1x256 marker scale
% p_sy  - 256x256 joint density p(PT,Y)
% pdf_s - 1x256, marginal density on pseudotime scale
%
%% Reference:
% 
% Reconstructing temporal and spatial dynamics from snap-shot data
% of heterogeneous cell populations
% Karsten Kuritz, Daniela Stöhr, Daniela Maichl, Nadine Pollak, Markus Rehm,
% Frank Allgöwer (2019)
%
%% ------------------------------------------------------

% parse input
if nargin < 3
	t_scale = 1;
end

% Joint density of pseudotime and marker -> pseudotime trajectory
[bw,p_sy,ss,yy] = kde2d_MAPiT([PT,Y],2^8,[0,min(Y)-range(Y)/4],[1,max(Y)+range(Y)/4],t_scale);


% marginal in s <=> pseudotime density
s		= ss(1,:);
y		= yy(:,1);
pdf_s		= trapz(y,p_sy,1);
