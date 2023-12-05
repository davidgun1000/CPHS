function [indx,u_res]=rs_multinomial_cond_corr2(particles,w,u)
%correlated conditional multinomial resampling

N=length(w); % number of particles
reference_particle=particles(1,1);
orig_index=(1:1:N);
col=[particles',w',orig_index'];
col_sort=sortrows(col,1);
particles_sort=col_sort(:,1);
weight_sort=col_sort(:,2);
orig_ind_sort=col_sort(:,3);

%indx=zeros(1,N); % preallocate 
Q=cumsum(weight_sort); % cumulative sum
Q(end)=1;
%u=sort(rand(1,N-1)); % random numbers
j=1;
for i=1:N-1
    while (Q(j)<u(i))
        j=j+1; % climb the ladder
    end
    indx_sort(i)=j; % assign index
end
tol=eps(reference_particle);
id_l=ismembertol(particles_sort, reference_particle*ones(N,1), tol);
ind=find(id_l>0);

%ind=find(orig_ind_sort==1);
if ind==1
   u_ref=unifrnd(0,Q(ind,1));
else
   u_ref=unifrnd(Q(ind-1,1),Q(ind,1));
end
%ll=length(u_ref)
ll=length(u_ref);
u_res=[u_ref(1,1),u]';
indx_sort=[ind(1,1),indx_sort];
%particles_sort=particles_sort';
indx=orig_ind_sort(indx_sort');
indx=indx';


end

