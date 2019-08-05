%% Example Spheroid data
%
% 0) preprocess data (load, normalize, remove outlier, ...)
% 1) Generate pseudotemporal ordering of cells with your favorite algorithm
% 2) Define true-scale distribution 
% 3) Get joint distribution of pseudotime and markers with "jointDensityPseudotimeYpre.m"
% 4) Get transformation with "preMAPiT.m"
% 5) Transform pseudotime trajectories to new scale with "MAPiT.m"
%
%% ----------------------------------------------------------------------------

% add toolbox 
addpath(genpath('../'))

% add wanderlust
wanderlustpath = fullfile(userpath,'/cyt3');
addpath(genpath(wanderlustpath))												% add wanderlust path

%% 0) preprocess data (load, normalize, remove outlier, ...)

% load spheroid data (HCT116, 11 days)
A = importdata('das2015-07-10_HCT_d11_Ki67_p27_HPy_001_P1.csv',';',21);
measured_states = strsplit(A.textdata{4},';');

% Normalize, extract, etc.

YoI			= {'V1 VioBlue-A','B2 PE-A','R1 APC-A','B1 FITC-A','FSC-A','SSC-A'};
y_names		= {'DNA', 'RNA', 'Ki67', 'p27', 'FSC'};
	
y_index = zeros(1,length(YoI));
for i = 1:length(YoI)
	pos = strcmp(YoI{i},measured_states);
    y_index(i) = find(pos);
end
data = A.data(:,y_index);

% remove negative values
neg = any(data<0,2);
data(neg,:) = [];

% remove outlier
out_prctile	= 1;
dataout_percentiles = false(size(data,1),1);
for i = 1:size(data,2)
	lowout = prctile(data(:,i), out_prctile, 1);
	upout  = prctile(data(:,i),100-out_prctile, 1);
	dataout_percentiles = dataout_percentiles | (data(:,i) < lowout) | (data(:,i) > upout);	
end
data(dataout_percentiles,:) = [];

% downsample to k datapoints
[N,n] = size(data);
k = 10000;
rng(2);
rand2use = randperm(N,k);
data = data(rand2use,:);

% RNA, Ki-67 and p27 to log scale
data(:,[2, 3, 4]) = log(data(:,[2, 3, 4]));


%% 1) Generate pseudotemporal ordering of cells with your favorite algorithm

% data for pseudotime algorithm (p27, RNA, Ki-67)
data_pt = data(:,[2, 3, 4]);
% normalize to [0, 1]
data_pt = (data_pt - min(data_pt)) ./ range(data_pt);

% wanderlust parameter
w_param.l				= 15;
w_param.k				= 8;
w_param.num_landmarks	= 80;
w_param.num_graphs		= 40;
w_param.metric			= 'euclidean';
w_param.normalize		= true;
w_param.band_sample		= true;
w_param.voting_scheme	= 'exponential';
w_param.flock_landmarks	= 2;
w_param.verbose			= false;
w_param.sig_factor		= 5;				% factor for WA kernel density estimation bandwidth
w_param.s = find((data_pt(:,1) > 0.9) & (data_pt(:,2) > 0.9) );   % wanderlust root cells 
G = wanderlust(data_pt,w_param);

% take the mean and normalize
PT	= mean(G.T)';
PT	= (PT - min(PT));
PT  = PT./max(PT);

%% 2) Define true-scale distribution 
% <=> real-time density (steady state age structure)

day		= 11;
radius	= 19.2 + 22.4 * day;
d_N		= 270;
r_N		= max([radius-d_N,0]);
x		= linspace(0,radius,55);
pdf_x	= 3 * (radius - x).^2 ./ (radius^3 - r_N^3); 

% prepare figure
a = length(y_names);
b = 2;
rect = [1, 1, 18, 20];
fh1 = figure('Color','w','Units','centimeters','Position',rect);


for i = 1:length(y_names)
	
%% 3) Get joint distribution of pseudotime and markers with "jointDensityPseudotimeYpre.m"

	Y = data(:,i);
	t_scale = 10^2;			% increase bandwith (by default too small)
	[s,y,pdf_s,p_sy] = jointDensityPseudotimeY(PT,Y,t_scale);
	
%% 4) Get transformation with "preMAPiT.m"

	[pre,tau,itau] = preMAPiT(s,pdf_s,x,pdf_x);
	
	
%% 5) Transform pseudotime trajectories to new scale with "MAPiT.m"
	
	[p_xy, test_pdf_x] = MAPiT(s, y, p_sy, itau, pre);
	
	% plot results
	subplot(a,b,1+(i-1)*b)
	pcolor(s,y,p_sy); shading interp
	xlabel('pseudotime')
	ylabel(y_names{i})
	title('Trajectory in pseudotime')
	
	subplot(a,b,2+(i-1)*b)
	pcolor(x,y,p_xy); shading interp
	xlim([0, radius]);
	xlabel('distance from surface')
	ylabel(y_names{i})
	title('Trajectory on spatial scale')
end


rect = [1, 1, 18, 10];
fh2 = figure('Color','w','Units','centimeters','Position',rect);

subplot(1,3,1)
plot(s,pdf_s)
xlabel('pseudotime')
ylabel('cell density')
title('Cell density in pseudotime')

subplot(1,3,2)
plot(x,pdf_x)
xlabel('distance from surface')
xlim([0, radius]);
ylabel('density')
title('Cell density on spatial scale')

subplot(1,3,3)
scatter(interp1(s,tau,PT),PT,'.')
xlabel('distance from surface')
xlim([0, radius]);
ylabel('pseudotime')
title('Transformation \tau^{-1}(x)')








