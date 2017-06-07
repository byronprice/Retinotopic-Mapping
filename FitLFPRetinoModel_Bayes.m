function [posteriorMean,posteriorInterval,posteriorSamples,posteriorMode] = FitLFPRetinoModel_Bayes(Response,xaxis,yaxis)
%FitLFPRetinoModel_Loglog.m
%   Use data from LFP retinotopic mapping experiment to fit a non-linear
%    model of that retinotopy (data is maximum LFP magnitude in window 
%    from 180 to 270msec minus minimum magnitude in window from 100 to
%    160 msec after stimulus presentation, assumes a
%    log-logistic likelihood and priors based on data previously collected
%

%Created: 2017/05/17, 24 Cummington Mall, Boston
% Byron Price
%Updated: 2017/05/17
% By: Byron Price

%  model has 7 parameters, defined by vector p
%  data ~ log-logistic(log mean, log scale), 
%    where the log mean of the data is as follows
%     log mean = p(1)*exp(-(xpos-p(2)).^2./(2*p(4)*p(4))-
%        (ypos-p(3)).^2./(2*p(5)*p(5)))+p(6);
%
%  and p(7) = log scale

numChans = size(Response,1);

numParameters = 8;

posteriorMode = zeros(numChans,numParameters);
posteriorMean = zeros(numChans,numParameters);
posteriorInterval = zeros(numChans,numParameters,2);
posteriorSamples = zeros(numChans,numParameters,(numIter-burnIn)/skipRate);

for zz=1:numChans
%     display(sprintf('Running Data for Channel %d...',zz));
    priorParams = zeros(numParameters,2);

    priorParams(1,:) = [0,0];% [0.125,0] to left for exp pdf, for inverse gamma [25.8169,12.6828];
    priorParams(2,:) = [max(xaxis)/2,700];% to left for normal, for gamma [50.788,26.021];
    priorParams(3,:) = [3,453/3];% gamma with slightly higher variance than mle [6.339,72.07];
    priorParams(4,:) = [25.508,11.71]; % gamma
    priorParams(5,:) = [19.786,12.543]; % gamma
    priorParams(6,:) = [950.71,0.006]; % gamma
    priorParams(7,:) = [31.9816,0.1048];% gamma for 1/parameter,
    %                for inverse gamma of parameter 7 [30.5995,9.1263];
    priorParams(8,:) = [1e-3,1e3]; % gamma, for precision on first parameter

    lBound = zeros(numParameters,1);lBound(1) = -Inf;
    lBound(2) = min(xaxis)-50;lBound(3) = min(yaxis)-50;lBound(8) = -Inf;
    uBound = [Inf,max(xaxis)+50,max(yaxis)+50,Inf,Inf,Inf,Inf,Inf]';

    numIter = 5.5e5;
    burnIn = 5e4;
    skipRate = 500;

    Data = Response{zz};
    reps = size(Data,1);
    
    flashPoints = Data(:,1:2);
    vepMagnitude = Data(:,3);
    
    posteriorProb = zeros(numIter,1);
    params = zeros(numParameters,numIter);
    
    for ii=[3,4,5,6,7]
        params(ii,1) = gamrnd(priorParams(ii,1),priorParams(ii,2));
    end
    params(1,1) = 1./gamrnd(26.6198,0.0764); % from inverse gamma above
    params(2,1) = normrnd(priorParams(2,1),priorParams(2,2));
    params(8,1) = log(1.5);
    
    inputParams = params(:,1);inputParams(7) = 1/inputParams(7);
    [loglikelihood] = GetLikelihood(reps,inputParams,vepMagnitude,flashPoints);

    logPrior = sum(log(gampdf(params(3:7,1),priorParams(3:7,1),priorParams(3:7,2)))); 
    logPrior = logPrior+log(normpdf(params(2,1),priorParams(2,1),priorParams(2,2)));
    alpha = exp(params(8,1));
    logPrior = logPrior+0.5*log(alpha)-0.5*alpha*params(1,1)^2;
    logPrior = logPrior+log(gampdf(alpha,priorParams(8,1),priorParams(8,2)));
    posteriorProb(1) = loglikelihood+logPrior;
    
    % FOR AUTOMATIC CREATION OF UPDATE MATRIX
    updateParam = logspace(-0.3,-3,burnIn);
    mu = zeros(numParameters,1);sigma = eye(numParameters);
    identity = eye(numParameters);
    for ii=1:numParameters
       sigma(ii,ii) = abs(sigma(ii,ii)*params(ii,1)/2); 
    end
    halfSigma = cholcov(sigma);
    loglambda = log(2.38^2/numParameters);
    updateMu = zeros(numParameters,1);
    
    for ii=[3,4,5,6,7]
        updateMu(ii,1) = gamrnd(priorParams(ii,1),priorParams(ii,2));
    end
    updateMu(1,1) = 1./gamrnd(26.6198,0.0764);
    updateMu(2,1) = normrnd(priorParams(2,1),priorParams(2,2));
    updateMu(8,1) = log(1);
    optimalAccept = 0.44;%0.234;
    
    for ii=2:burnIn
        Z = mvnrnd(mu,exp(loglambda).*sigma)';
        pStar = params(:,ii-1)+Z;
        
        if sum(pStar<=lBound) == 0 & sum(pStar>=uBound) == 0
            inputParams = pStar;inputParams(7) = 1/inputParams(7);
            pStarLogLikelihood = GetLikelihood(reps,inputParams,vepMagnitude,flashPoints);
            pStarLogPrior = sum(log(gampdf(pStar(3:7),priorParams(3:7,1),priorParams(3:7,2))));
            alpha = exp(pStar(8));
            pStarLogPrior = pStarLogPrior+0.5*log(alpha)-0.5*alpha*pStar(1)^2;
            pStarLogPrior = pStarLogPrior+log(normpdf(pStar(2),priorParams(2,1),priorParams(2,2)));
            pStarLogPrior = pStarLogPrior+log(gampdf(alpha,priorParams(8,1),priorParams(8,2)));
%             
%             if pStar(2) < priorParams(2,1) || pStar(2) > priorParams(2,2) || ...
%                     pStar(3) < priorParams(3,1) || pStar(3) > priorParams(3,2)
%                 pStarLogPrior = pStarLogPrior-1e20;  % for uniform prior
%                         % on x and y centroid position
%             end
            
            logA = (pStarLogLikelihood+pStarLogPrior)-posteriorProb(ii-1);
            
            if log(rand) < logA
                params(:,ii) = pStar;
                posteriorProb(ii) = pStarLogLikelihood+pStarLogPrior;
            else
                params(:,ii) = params(:,ii-1);
                posteriorProb(ii) = posteriorProb(ii-1);
            end
            
            meanSubtract = params(:,ii)-updateMu;
            updateMu = updateMu+updateParam(ii).*meanSubtract;
            halfSigma = halfSigma+updateParam(ii).*(triu((halfSigma^-1)*(halfSigma'*halfSigma+meanSubtract*...
            meanSubtract')*((halfSigma^-1)')-identity)-halfSigma);
            sigma = halfSigma'*halfSigma;
            loglambda = loglambda+updateParam(ii).*(exp(min(0,logA))-optimalAccept);
            
        else
            params(:,ii) = params(:,ii-1);
            posteriorProb(ii) = posteriorProb(ii-1);
            loglambda = loglambda+updateParam(ii).*(-optimalAccept);
        end
%         error = abs(1394-updateMu(2))+abs(419-updateMu(3));
%         scatter(ii,posteriorProb(ii));%scatter(ii,error);
%         hold on;pause(0.01);
    end
    [V,D] = eig(cov(params(:,1e4:100:burnIn)'));
    W = V*sqrtm(D);
    eigenvals = diag(W'*W);
    
    tempW = [];
    tempEigs = [];
    for jj=1:numParameters
        if eigenvals(jj) > 1e-5
            tempW = [tempW,W(:,jj)];
            tempEigs = [tempEigs,eigenvals(jj)];
        end
    end
    
    W = fliplr(tempW);
    eigenvals = fliplr(tempEigs);
    q = length(eigenvals);
    p = ones(q,1)./q;
    loglambda = loglambda.*ones(q,1);
    updateParam = 1e-2;
    
    acceptRate = 0;
    for ii=burnIn+1:numIter
        index = find(mnrnd(1,p)==1);
        lambda = loglambda(index);
        stdev = sqrt(exp(lambda).*eigenvals(index));
        pStar = params(:,ii-1)+W(:,index)*normrnd(0,stdev);
        
        if sum((pStar-lBound)<=0) == 0 & sum((pStar-uBound)>0) == 0
            inputParams = pStar;inputParams(7) = 1/inputParams(7);
            pStarLogLikelihood = GetLikelihood(reps,inputParams,vepMagnitude,flashPoints);
            pStarLogPrior = sum(log(gampdf(pStar(3:7),priorParams(3:7,1),priorParams(3:7,2))));
            alpha = exp(pStar(8));
            pStarLogPrior = pStarLogPrior+0.5*log(alpha)-0.5*alpha*pStar(1)^2;
            pStarLogPrior = pStarLogPrior+log(normpdf(pStar(2),priorParams(2,1),priorParams(2,2)));
            pStarLogPrior = pStarLogPrior+log(gampdf(alpha,priorParams(8,1),priorParams(8,2)));
            
%             if pStar(2) < priorParams(2,1) || pStar(2) > priorParams(2,2) || ...
%                     pStar(3) < priorParams(3,1) || pStar(3) > priorParams(3,2)
%                 pStarLogPrior = pStarLogPrior-1e20;
%             end
            
            logA = (pStarLogLikelihood+pStarLogPrior)-posteriorProb(ii-1);
            
            if log(rand) < logA
                params(:,ii) = pStar;
                posteriorProb(ii) = pStarLogLikelihood+pStarLogPrior;
                acceptRate = acceptRate+1;
            else
                params(:,ii) = params(:,ii-1);
                posteriorProb(ii) = posteriorProb(ii-1);
            end
            lambda = lambda+updateParam.*(exp(min(0,logA))-optimalAccept);
        else
            params(:,ii) = params(:,ii-1);
            posteriorProb(ii) = posteriorProb(ii-1);
            lambda = lambda+updateParam.*(-optimalAccept);
        end
        loglambda(index) = lambda;
    end
    params(7,:) = 1./params(7,:);
    params(8,:) = 1./exp(params(8,:));
    posteriorSamples(zz,:,:) = params(:,burnIn+1:skipRate:end);
    
    figure();
    numRows = ceil(numParameters/2);
    for ii=1:numParameters
       data = squeeze(posteriorSamples(zz,ii,:));
       subplot(numRows,2,ii);histogram(data);
    end
    
    [~,ind] = max(posteriorProb);
    posteriorMode(zz,:) = params(:,ind)';
    posteriorMean(zz,:) = mean(squeeze(posteriorSamples(zz,:,:)),2)';
    
    display(posteriorMean(zz,:));
    
    alpha = 0.05;
    posteriorInterval(zz,:,:) = quantile(squeeze(posteriorSamples(zz,:,:)),[alpha/2,1-alpha/2],2);
end
end

% log-logistic likelihood for non-linear 2D Gaussian function
function [loglikelihood] = GetLikelihood(reps,parameterVec,vepMagnitude,flashPoints)
loglikelihood = zeros(reps,1);
for kk=1:reps
    mu = parameterVec(1)*exp(-((flashPoints(kk,1)-parameterVec(2)).^2)./(2*parameterVec(4).^2)-...
        ((flashPoints(kk,2)-parameterVec(3)).^2)./(2*parameterVec(5).^2))+parameterVec(6);
    loglikelihood(kk) = -log(parameterVec(7))-log(vepMagnitude(kk))+(log(vepMagnitude(kk))-mu)/parameterVec(7)-...
        2*log(1+exp((log(vepMagnitude(kk))-mu)/parameterVec(7)));
end
loglikelihood = sum(loglikelihood);
end

% gamma likelihood
% function [loglikelihood] = GetLikelihood(reps,parameterVec,vepMagnitude,flashPoints)
% loglikelihood = zeros(reps,1);
% for kk=1:reps
%     mu = parameterVec(1)*exp(-((flashPoints(kk,1)-parameterVec(2)).^2)./(2*parameterVec(4).^2)-...
%         ((flashPoints(kk,2)-parameterVec(3)).^2)./(2*parameterVec(5).^2)-...
%         parameterVec(7)*(flashPoints(kk,1)-parameterVec(2))*(flashPoints(kk,2)-parameterVec(3))/(2*parameterVec(4)*parameterVec(5)))+parameterVec(6);
%     loglikelihood(kk) = (-vepMagnitude(kk)/mu-log(mu))/parameterVec(7)-log(parameterVec(7))/parameterVec(7)+...
%            (1/parameterVec(7)-1)*log(vepMagnitude(kk))-log(gamma(1/parameterVec(7)));
% end
% end

