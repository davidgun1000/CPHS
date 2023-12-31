function f = IACT_used(chain,K0)
% compute IACT using autocorrelations with adaptive cutoff
[N,dim] =size(chain);
if (nargin<2 || isempty(K0)) K0 = 1000; end

for i = 1:dim
    rho = autocorr(chain(:,i),K0);
    rho = rho(2:length(rho));
    %aux = N-(1:K0)';
    aux = N;
    aux = 2./sqrt(aux);
    flag = abs(rho)<aux;
    ind = 1;
    while (ind<=K0)&&(~flag(ind)) ind = ind+1; end
    L = min(K0,ind);
    f(i) = 1+2*sum(rho(1:L));
end
end

