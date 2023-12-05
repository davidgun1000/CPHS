function [newx,neww]=multidimsort(x,w)

% Exaample: x=randn(D,N) where D=dimensionality, N=sample size

N=size(x,2); % number of particles
ndss=zeros(1,N); % preallocate
indxx=1:N; % list of indexes
[~,ind]=min(x(1,:)); % select index of a particle having least value along the x dimension
ndss(1)=ind; % store
indxx(ind)=[]; % delete the selected index from index list
for p=2:N-1,
    svalues=sqrt(sum(bsxfun(@minus,x(:,ind),x(:,indxx)).^2,1)); % distance of selected particle to every other particle
    [~,ix]=min(svalues); % find a particle closest to the selected particle
    ind=indxx(ix); % select its corresponding index
    ndss(p)=ind; % store
    indxx(ix)=[]; % delete the selected index from the index list 
end
ndss(N)=indxx; % the last index is the remaining index

newx=x(:,ndss); % sort particles
neww=w(:,ndss); % sort particles