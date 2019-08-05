function ywant = nanmoving_average(x,y,xwant,binsize)

[~,nx] = size(x);
[d,n] = size(y);
xi = xwant(:);
[~,N] = size(xwant);

if nx ~= n
	error('length in the second dimension must be same for x and y')
end

% moving horizont
inbin = (bsxfun(@ge,x,-binsize/2+xi) & bsxfun(@lt,x,binsize/2+xi));

% nans in y
nanas = isnan(y);

% preallocate output
ywant = zeros(d,N);
for i=1:N
	ywant(:,i) = nansum(y(:,inbin(i,:)),2) ./ (sum(inbin(i,:))- sum(nanas(:,inbin(i,:)),2));
end





































