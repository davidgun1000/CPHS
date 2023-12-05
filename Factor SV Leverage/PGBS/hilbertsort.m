function ndss = hilbertsort(x)

% Hilbert sorting

D=size(x,1); % dimensionality
N=size(x,2); % number of particles
if D==1, % univariate sorting
    [~,ndss]=sort(x,2);
elseif D>1, % multivariate hilbert sorting
    rprecision=8; % precision
    hn=zeros(1,N); % preallocate
    xn=floor(2^rprecision*(normcdf(x'))); % translate to uniforms
    xnn=uint32(xn); % translate to 32 bit integers
    for j=1:N,
        h1=HilbertIndexTransposed(rprecision,D,xnn(j,:));
        hn(j)=get_hilbert_general(h1,rprecision);
    end
    [~,ndss]=sort(hn,2);
end