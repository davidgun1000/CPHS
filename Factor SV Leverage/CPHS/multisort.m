function [newp,neww,nds]=multisort(vals,particles,w)
nds=hilbertsort(vals); % sorting indices
newp=particles(:,nds); % sort particles
neww=w(:,nds); % sort weights