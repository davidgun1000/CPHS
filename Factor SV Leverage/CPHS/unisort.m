function [newp,neww,nds]=unisort(vals,particles,w)
[~,nds]=sort(vals); % sorting indices
newp=particles(:,nds); % sort particles
neww=w(:,nds); % sort weights