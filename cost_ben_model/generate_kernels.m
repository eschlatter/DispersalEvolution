%% Nonparametric
clear

nbins = 6; % how many dispersal bins to use (2<=nbins)
ndiv = 6; % increments by which to shift when generating kernels

bars = nchoosek(1:(nbins+ndiv-1),(nbins-1)); % locations of bars in "bars and stars" proof
bars = [zeros(length(bars),1) bars (nbins+ndiv)*ones(length(bars),1)];

kernels = ones(length(bars),nbins);

for j=1:length(bars)
    for i=1:nbins
        kernels(j,i) = bars(j,i+1)-bars(j,i)-1;
    end
end

kernels = (1/ndiv)*kernels;

%% Weibull

clear

nbins = 6; % how many dispersal bins to use (2<=nbins)

% values of lambda and k s.t. cdf at x=6 is greater than 0.99
ks = 1:0.5:5;
lambdas = (-6*log(0.01)).^(1./ks);
plot(ks,lambdas)
xlabel('k')
ylabel('lambda')


lambda = 5;
k = 1.5;
plot(0:0.1:10,fn_weibull_pdf(0:0.1:10,lambda,k))
fn_weibull_cdf(60,lambda,k)


hold on
plot(0:0.1:2.5,fn_weibull_pdf(0:0.1:2.5,1,1.5))
plot(0:0.1:2.5,fn_weibull_pdf(0:0.1:2.5,2,1.5))
hold off

a = Weibull(1:6,lambda,k);