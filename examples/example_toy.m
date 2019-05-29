%% Example Toy Data
%
% 1) Generate toy data with pseudotemporal ordering 
% 2) Define true-scale distribution 
% 3) Get joint distribution of pseudotime and markers with "jointDensityPseudotimeYpre.m"
% 4) Get transformation with "preMAPiT.m"
% 5) Transform pseudotime trajectories to new scale with "MAPiT.m"
%
%% ----------------------------------------------------------------------------

% add toolbox 
addpath(genpath('../'))

%% 1) Generate toy data with pseudotemporal ordering 
N = 10000;
% Pseudotime values
a = 3;
b = 5;
PT = betarnd(a,b,[N,1]);

% marker values 
Y =  randn(N,1)*0.1 + PT;

%% 2) Define true-scale distribution
% <=> real-time density (for example, uniform, steady state age structure, ...)
x = linspace(0,1,55);
pdf_x = ones(size(x)); % uniform

%% 3) Get joint distribution of pseudotime and markers with "jointDensityPseudotimeYpre.m"
tic
[s,y,pdf_s,p_sy] = jointDensityPseudotimeY(PT,Y);

%% 4) Get transformation with "preMAPiT.m"

[pre,tau,itau] = preMAPiT(s,pdf_s,x,pdf_x);

%% 5) Transform pseudotime trajectories to new scale with "MAPiT.m"

[p_xy, test_pdf_x] = MAPiT(s, y, p_sy, itau, pre);
toc

%% plot results
rect = [1, 1, 18, 12];
fh = figure('Color','w','Units','centimeters','Position',rect);

subplot(2,3,1)
scatter(PT,Y)
xlabel('pseudotime')
ylabel('marker')
title('toy data')

subplot(2,3,5)
pcolor(s,y,p_sy); shading interp
xlabel('pseudotime')
ylabel('marker')


subplot(2,3,6)
pcolor(x,y,p_xy); shading interp
xlabel('real-time')
ylabel('marker')

subplot(2,3,2)
plot(s,pdf_s)
xlabel('pseudotime')
ylabel('density')

subplot(2,3,3)
plot(x,pdf_x)
hold on
plot(x,test_pdf_x,'--')
ylim([0,1.1])
xlabel('real-time')
ylabel('density')
legend('original','transformed')

subplot(2,3,4)
scatter(interp1(s,tau,PT),PT,'.')
hold on
plot(x,betainv(x,a,b))
legend('derived by MAPiT','original')
xlabel('real-time')
ylabel('pseudotime')
title('Transformation \tau^{-1}(x)')








