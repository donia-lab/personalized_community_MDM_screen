function power = welch_power(mu1,mu2,sd1,sd2,n,alpha)

%This script estimates the power of a welch t test with unequal variances
%and equal sample sizes using methods outlined in Harrison and Brady (2004)
%The test is one-sided. 

s1 = sd1^2;
s2 = sd2^2; 

nu = (s1/n + s2/n)^2/(1/(n-1)*(s1/n)^2 + 1/(n-1)*(s2/n)^2);

delta = (mu2 - mu1)/sqrt(s1/n + s2/n);

if abs(mu1-mu2) == 0
    power = 0;
else
    t = tinv(1-alpha,nu);
    power = nctcdf(t,nu,delta,'upper');
end

end

