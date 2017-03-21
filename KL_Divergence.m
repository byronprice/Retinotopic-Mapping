function [D_hat] = KL_Divergence(P_data,Q_data)
%KL_Divergence.m
%   Calculate an estimator for the Kullback-Leibler divergence between two
%    empirical distributions.
%
%INPUT: P_data- 1D dataset (the target distribution)
%       Q_data- a second 1D dataset (the comparison
%         distribution)
%OUTPUT: D_hat- estimator for the KL divergence
%
%Created: 2017/03/19, 24 Cummington Mall
%  Byron Price
%Updated: 2017/03/19
% By: Byron Price

% see Fernando Perez-Cruz ... Kullback-Leibler Divergence Estimation of
% Continuous Distributions
% 
% for comparison ... the analytical solution for two univariate Guassian
% distributions with p ~ normal(mu1,sigma1) and q ~ normal(mu2,sigma2)
%  D(P,q) = log(sigma2/sigma1)+(sigma1^2+(mu1-mu2)^2)/(2*sigma2^2)-0.5

minimumP = min(P_data);maximumP = max(P_data);
minimumQ = min(Q_data);maximumQ = max(Q_data);

N = 1e5;
dataRange = linspace(min(minimumP,minimumQ),max(maximumP,maximumQ),N);

% get empirical CDF's
[P,xp] = GetECDF(P_data);
[Q,xq] = GetECDF(Q_data);

% convert to piecewise linear continuous distributions
P_continuous = zeros(N,1);
Q_continuous = zeros(N,1);

P_continuous(dataRange<minimumP) = 0;
P_continuous(dataRange>=maximumP) = 1;
for ii=1:length(xp)-1
    inds = find(dataRange>=xp(ii) & dataRange<xp(ii+1));
    m = (P(ii+1)-P(ii))/(xp(ii+1)-xp(ii));
    b = P(ii)-m*xp(ii);
    P_continuous(inds) = m.*dataRange(inds)+b;
end

Q_continuous(dataRange<minimumQ) = 0;
Q_continuous(dataRange>=maximumQ) = 1;
for ii=1:length(xq)-1
    inds = find(dataRange>=xq(ii) & dataRange<xq(ii+1));
    m = (Q(ii+1)-Q(ii))/(xq(ii+1)-xq(ii));
    b = Q(ii)-m*xq(ii);
    Q_continuous(inds) = m.*dataRange(inds)+b;
end

deltaX = mean(diff(dataRange));

D_hat = 0;
count = 1;
for ii=2:N
    deltaP = P_continuous(ii)-P_continuous(ii-1);
    deltaQ = Q_continuous(ii)-Q_continuous(ii-1);
    if deltaP ~= 0 && deltaQ ~= 0
        D_hat = D_hat+(log2(deltaP/deltaX)-log2(deltaQ/deltaX));
        count = count+1;
    end
end
D_hat = D_hat/count;
display(count/(N-1));
end

function [f,x] = GetECDF(data)
x = sort(data);
N = length(data);
f = zeros(N,1);
for ii=1:N
    f(ii) = sum(
end

end
