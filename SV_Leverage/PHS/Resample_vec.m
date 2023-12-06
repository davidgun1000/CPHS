function Y=Resample(X,weights,N);

[PH,bin]=histc(rand(N,1),cumsum([0;weights]));

Y=X(bin,:);