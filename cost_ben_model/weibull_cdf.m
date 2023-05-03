function F = weibull_cdf(x, lambda, k)
% calculate the CDF of the Weibull distribution
F = 1 - exp(-((x./lambda).^k));