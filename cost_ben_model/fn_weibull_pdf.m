function f = fn_weibull_pdf(x,lambda,k)
% calculate the pdf of the Weibull distribution
f = (k./lambda).*((x./lambda).^(k-1)).*exp(-((x./lambda).^k));