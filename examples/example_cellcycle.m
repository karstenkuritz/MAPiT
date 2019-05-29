%% Example Cell cycle data
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
% load cell cycle data
A = importdata('das2015-02-04_H460 gem Cyclin B.001_P1.csv',';',21);

% DNA and geminin are in columns 7 and 8
data = A.data(:,[7,8]);


% remove negative values
neg = any(data<0);
data(neg,:) = [];

% remove outlier
out_prctile	= 0.1;
dataout_percentiles = false(size(data,1),1);
for i = 1:size(data,2)
	lowout = prctile(data(:,i), out_prctile, 1);
	upout  = prctile(data(:,i),100-out_prctile, 1);
	dataout_percentiles = dataout_percentiles | (data(:,i) < lowout) | (data(:,i) > upout);	
end
data(dataout_percentiles,:) = [];


% Y2 to log scale
data(:,2) = log(data(:,2));

% normalize
ndata = (data - min(data)) ./ range(data);

%% 1) Generate pseudotemporal ordering of cells with your favorite algorithmw_param					= [];
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
w_param.s = find(ndata(:,2) < 0.1);   % wanderlust root cells 
G = wanderlust(ndata,w_param);

% take the mean and normalize
PT	= mean(G.T)';
PT	= (PT - min(PT));
PT  = PT./max(PT);

Y1 = data(:,1);
Y2 = data(:,2);

%% 2) Define true-scale distribution 
% <=> real-time density (steady state age structure)

T = 16;
x = linspace(0,T,55);
pdf_x = 2*log(2)/T*2.^(-x./T);


%% 3) Get joint distribution of pseudotime and markers with "jointDensityPseudotimeYpre.m"

t_scale = 10^4;			% increase bandwith (by default too small)
[s,y,pdf_s,p_sy] = jointDensityPseudotimeY(PT,Y2,t_scale);


%% 4) Get transformation with "preMAPiT.m"

[pre,tau,itau] = preMAPiT(s,pdf_s,x,pdf_x);


%% 5) Transform pseudotime trajectories to new scale with "MAPiT.m"

[p_xy, test_pdf_x] = MAPiT(s, y, p_sy, itau, pre);


%% plot results

subplot(2,3,1)
scatter(Y1,Y2,[],PT)
h = colorbar;
ylabel(h, 'pseudotime')
xlabel('DNA')
ylabel('log(geminin)')
title('cell cycle data')

subplot(2,3,5)
pcolor(s,y,p_sy); shading interp
xlabel('pseudotime')
ylabel('log(geminin)')
title('Trajectory in pseudotime')

subplot(2,3,6)
pcolor(x,y,p_xy); shading interp
xlim([0, T]);
xlabel('real-time [h]')
ylabel('log(geminin)')
title('Trajectory in real-time')

subplot(2,3,2)
plot(s,pdf_s)
xlabel('pseudotime')
ylabel('cell density')
title('Cell density in pseudotime')

subplot(2,3,3)
plot(x,pdf_x)
hold on
% plot(xx,test_pdf_x,'--')
xlabel('real-time [h]')
xlim([0, T]);
ylabel('density')
title('Cell density in real-time')
% legend('original','transformed')

subplot(2,3,4)
scatter(interp1(s,tau,PT),PT,'.')
xlabel('real-time [h]')
xlim([0, T]);
ylabel('pseudotime')
title('Transformation \tau^{-1}(x)')

