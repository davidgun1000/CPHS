function [indx]=rs_multinomial_corr2(particles,w,u)
%correlated multinomial resampling

N=length(w); % number of particles
orig_index=(1:1:N);
col=[particles',w',orig_index'];
col_sort=sortrows(col,1);
particles_sort=col_sort(:,1);
weight_sort=col_sort(:,2);
orig_ind_sort=col_sort(:,3);
%indx_sort=zeros(1,N); % preallocate 
Q=cumsum(weight_sort); % cumulative sum
Q(end)=1;

i=1;
j=1;
while (Q(j)<u(i)),
   j=j+1;
end;
indx_sort(i)=j;

j=1;
for i=2:N
    while (Q(j)<u(i))
        j=j+1; % climb the ladder
    end
    indx_sort(i)=j; % assign index
end
%particles_sort=particles_sort';
indx=orig_ind_sort(indx_sort');
indx=indx';


end



