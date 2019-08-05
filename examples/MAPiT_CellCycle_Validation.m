%% Comparing geminin trajectories obtained by MAPiT vs. life-cell imaging
% Reproducing Fig 2 in: 
% Reconstructing temporal and spatial dynamics in single-cell experiments
% Karsten Kuritz, Daniela Stöhr, Daniela Maichl, Nadine Pollak, Markus Rehm, Frank Allgöwer
% bioRxiv 697151; doi: https://doi.org/10.1101/697151 
%
% 0) load single cell trajectories 
% 1) Rescale time such that birth is at t=0
% 2) Smoothen curves
% 3) Remove background and scale to fit MAPiT result
%
%% ----------------------------------------------------------------------------
% add toolbox 
addpath(genpath('../'))

% some variables to load single cell trajectories
cell_trackno = [2,3];	% cells with a full cycle 
N = length(cell_trackno); 
K = 6;					% # Detection fields
C = K*N;				% Total cell number
c = 1;

% Preallocate single cell data struct
Cell(C).file = 'test';

%% 0) load single cell trajectories 
for k = 1:K
	
	% Load data from detection field
	file = sprintf('180213KK1_p0001-00%gKK_QTFyDetectionCh1SegMethod1.csv',k);
	T = readtable(file);
	
	for n = 1:N
		
		% get single cell 
		this_cell = T.TrackNumber == cell_trackno(n);
		
		% Save to single cell data struct
		Cell(c).file = file;
		Cell(c).y = T.MeanNoBgCorrectedCh01(this_cell);
		Cell(c).t = T.TimePoint(this_cell);
		
		c = c+1;
	end
end

%% 1) Rescale time such that birth is at t=0
for c = 1:C
	Cell(c).t = (Cell(c).t - Cell(c).t(1)) / 4;
end

% Common time vector
allt = vertcat(Cell.t);
newt = min(allt):0.25:max(allt);

% interpolate single cell trajectories to common time vector
tn = 3; % cell division in last 3 timepoints (45 min)  
yall = zeros(length(newt),C);
for ic = 1:C
	yall(:,ic) = interp1(Cell(ic).t(1:end-tn),Cell(ic).y(1:end-tn),newt,'linear',nan);
end
yall(yall==0) = nan;

%% 2) Smoothen curves
binsize = 0.5;
yma = nanmoving_average(newt,yall',newt,binsize)';
% remove bad values
yma(isinf(yma)|yma<=0.01) = nan;


%% 3) Remove background and scale to fit MAPiT result
% load MAPiT cell cycle data
example_cellcycle

% median from MAPiT result
MAPiT_median = cdflevelsets(p_xy,y,0.5);

% bring median to common time vectore
MAPiT_median_ct = interp1(x,MAPiT_median,newt)';

% find optimal scaling p(1) and background p(2) parameter
normaly = @(y,p) (y - p(2)).*p(1);
obj = @(p) nansum(nansum((exp(MAPiT_median_ct) - (normaly(yma,p))).^2));
popt = fminsearch(obj,[0.05,0.1]);

% scale and take the logarithm
yn = normaly(yma,popt);
yn(yn<=0) = nan;
yn = log(yn);

%% plot results
figure;
pcolor(x,y,p_xy); shading interp
xlim([0, T]);
ylim([-8.5, -2]);
hold on
plot(newt,MAPiT_median_ct,'k--')
plot(newt,yn,'r')

legend('MAPiT','MAPiT median','microscopy','Location','northwest')
xlabel('real-time [h]')
ylabel('log(geminin)')
title('Trajectories in real-time')




