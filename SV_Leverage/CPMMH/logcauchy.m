function [logpdf]=logcauchy(tau)
%logcauchy prior
logpdf=-0.5*log(tau)-log(pi*(1+tau));