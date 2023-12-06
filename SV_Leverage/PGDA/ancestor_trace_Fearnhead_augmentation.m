function [sirhat]=ancestor_trace_Fearnhead_augmentation(particles,w,indx)
%keyboard
% ancestor tracing algorithm for PGDA
T=size(w,1);
N=size(w,2);
outndx=NaN(1,T);
outndx(T)=randsample(N,1,true,w(T,:));
sirhat(:,T)=particles(:,outndx(T),T);
for t=T-1:-1:1,
    outndx(t)=indx(t+1,outndx(t+1));
    sirhat(:,t)=particles(:,outndx(t),t);
end