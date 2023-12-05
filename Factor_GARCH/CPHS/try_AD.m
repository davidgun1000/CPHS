z=unifrnd(0,1,100000,1);
x=2*z-1;
sigma=10^-10;
epsi=(10^-10);
y=(x-epsi).*normcdf(x./sigma)+(x+epsi).*(1-normcdf(x./sigma));
yy=(y+1)./2;
