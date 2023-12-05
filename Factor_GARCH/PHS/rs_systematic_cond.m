function indx_new=rs_systematic_cond(w,conditional,usamp)

%if nargin<3,
%    usamp=rand; % random number
%end

%keyboard
N=length(w); % number of particles/output indices
Ntimesw=N*w(1,1);
if Ntimesw<=1
   usamp=unifrnd(0,Ntimesw);
else
   r1=Ntimesw-floor(Ntimesw);
   prob=(r1*(floor(Ntimesw)+1))/Ntimesw;
   A1=rand();
   if A1<=prob
      usamp=unifrnd(0,r1);
   else
      usamp=unifrnd(r1,1);
   end
end

indx=zeros(1,N); % preallocate index variable
Q=cumsum(w); % cumulative sum
u=((0:N-1)+usamp)/N; % set strata

j=1;
for i=1:N
    while (Q(j)<u(i))
        j=j+1;
    end
    indx(i)=j;
end

cycle = randsample(length(find(indx==1)),1); % select one index containing one
indx_new = [indx(cycle:N),indx(1:cycle-1)]; % cycle the indices



% indx1=find(indx==1);
% indx1_l=length(indx1);
% cycle=randi([0,indx1_l-1]);
% 
% for i=1:N
%     if cycle+i<=N
%         cycle_ind(1,i)=cycle+i;
%     else
%         cycle_ind(1,i)=cycle+i-N;
%     end
% end
% 
% indx_new=indx(1,cycle_ind);


% if indx1_l==0
%    indx_new(1)=1; 
% end




%if conditional
%    indx(1)=1; % static conditional index
%end