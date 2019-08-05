%% Comparing Ki-67 intensities across spheroid obtained by MAPiT vs. imaging of slices
% Reproducing Fig 4 in: 
% Reconstructing temporal and spatial dynamics in single-cell experiments
% Karsten Kuritz, Daniela Stöhr, Daniela Maichl, Nadine Pollak, Markus Rehm, Frank Allgöwer
% bioRxiv 697151; doi: https://doi.org/10.1101/697151 
%
% 1) Load intensity profiles in spheroid slices
% 2) Smooth raw data
% 3) Load MAPiT results from 'example_spheroid.m'
% 4) Scale to match MAPiT results
%
%% ----------------------------------------------------------------------------

% add toolbox 
addpath(genpath('../'))


%% 1) Load intensity profiles in spheroid slices
slice_files = {'Profil_1_T0_Z0.csv','Profil_2_T0_Z0.csv','Profil_3_T0_Z0.csv',...
	'Profil_4_T0_Z0.csv','Profil_5_T0_Z0.csv','Profil_6_T0_Z0.csv',...
	'Profil_7_T0_Z0.csv','Profil_8_T0_Z0.csv','Profil_9_T0_Z0.csv'};

n = length(slice_files);
for k = 1:n
	A = importdata(slice_files{k},',',8);
	S(k).x = A.data(:,1)/1000;
	S(k).y = A.data(:,4);
end

% plot raw data
fh1 = figure('Color','w');
for k=1:n
	plot(S(k).x,S(k).y)
	hold on
end
legend(slice_files, 'Interpreter', 'none')
title('Transversal quantified Ki-67 intensities')
ylabel('Ki-67')
xlabel('Distance [µm]')


%% 2) extract data from outer rim 200 µm into the spheroid and smoothen it

% x-scale 0-200 µm
xx = linspace(0,200,201);

% preallocate y-matrix 
y_init = zeros(length(xx),2*n);

% 1) start to core ->
for k=1:n
	thisdata = [S(k).x, S(k).y];
	
	% find most outer rim 
	small = find(smooth(S(k).y,5) > 30);
	inside = (S(k).x > 100);
	% take only these 
	thisdata_take = thisdata(small(1):end,:);
	% shift distance such that it starts with 0
	thisdata_take(:,1) = thisdata_take(:,1) - thisdata_take(1,1);
	
	% extract 0-200 µm 
	ywant = interp1(thisdata_take(:,1),thisdata_take(:,2),xx);
	% smooth
	y_init(:,k) = smooth(ywant,50);
end

% 2) end to core <-
for k=1:n
	thisdata = [S(k).x, S(k).y];
	% flip and start from the other side
	thisdata = flipud(thisdata);
	thisdata(:,1) = -(thisdata(:,1)-thisdata(1,1));
	
	% find most outer rim 
	small = find(smooth(thisdata(:,2),5) > 30);
	inside = (thisdata(:,1) > 50);
	% take only these 
	thisdata_take = thisdata(small(1):end,:);
	% shift distance such that it starts with 0
	thisdata_take(:,1) = thisdata_take(:,1) - thisdata_take(1,1);
	
	% extract 0-200 µm 
	ywant = interp1(thisdata_take(:,1),thisdata_take(:,2),xx);
	% smooth
	y_init(:,n+k) = smooth(ywant,50);
end

% plot profiles w.r.t. distance from the surface
fh2 = figure('Color','w');
for k=1:n
	plot(xx,y_init)
	hold on
end
title('Ki-67 profiles from outer rim to the core')
ylabel('Ki-67')
xlabel('Distance from surface [µm]')


%% 3) Load MAPiT Ki67 intensity distribution of Ki67, HCT-116 day 11 #1
% use example spheroid script - for comparison we need x,y,p_xy of MAPiT results

example_spheroid
t_scale = 10^2;			% increase bandwith (by default too small)
Ki67_column = 3;
[s,y,pdf_s,p_sy] = jointDensityPseudotimeY(PT, data(:,Ki67_column), t_scale);
[pre,tau,itau] = preMAPiT(s, pdf_s, x, pdf_x);
[p_xy, test_pdf_x] = MAPiT(s, y, p_sy, itau, pre);

%% 4) find scaling factor by fitting median against density from MAPiT 
% median of individual profiles
ymedian	= nanmedian(log(y_init'));
s0 = -5; % Starting guess
objfun = @(s) -sum(interp2(x,y,p_xy,xx,ymedian+s,'linear',eps));
[opt_s,fval,exitflag,output] = fminunc(@(x) objfun(x),s0);

% scale intensity profiles
yy = log(y_init')+opt_s; % log(y)+s


%% Plot results
fh3 = figure('Color','w')
% conditional probabilty
p_xy_conditional = p_xy./trapz(y,p_xy);
p0 = pcolor(x,y,p_xy_conditional); shading interp
hold on

p1 = plot(xx,yy,'Color',[1 0 0 0.5]);
p2 = plot(xx,nanmedian(yy),'y');
p3 = plot(xx,prctile(yy,[25,75],1),'y--');
ylabel('Ki-67')
xlabel('Distance from surface')
title('Ki-67 profiles: MAPiT vs. microscopy')
legend([p0,p1(1),p2,p3(1)],'MAPiT','microscopy','microscopy median','microscopy quartiles')
xlim([0,200])
ylim([-10,-3])



	
