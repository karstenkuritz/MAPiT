function cdflevelsets = cdflevelsets(pdf,y,levelsets)



[n,m] = size(pdf);
% make sure that not all are 0 in onerow
rowswithallzeros =  all(pdf==0);
pdf(:,rowswithallzeros) = eps;


% Cumulative Probability Distribution of the readout at each position

cdf  = bsxfun(@times,cumtrapz(y,pdf),trapz(y,pdf).^-1);

% levelsets of the cumulative distribution are the important ones and are
% obtained by evaluating the inverse of the cdf. For the inverse with
% interp1 the cumDistribution must be strictly monotonic. Thus we remove
% values that appear more than once.

cdflevelsets = zeros(length(levelsets),m);
for j = 1:m
	[~,notdoubles,~] = unique(cdf(:,j));
	cdflevelsets(:,j) = interp1(cdf(notdoubles,j),y(notdoubles),levelsets);
end

%% Test Plot
% [X,Y] = meshgrid(1:m,y);

% pcolor(X,Y,cdf); shading interp
% hold on
% plot(1:m,cdflevelsets,'r')
