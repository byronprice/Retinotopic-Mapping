function [posteriorMedian,posteriorInterval,posteriorSamples,posteriorMode] = FitLFPRetinoModel_Bayes(Response,xaxis,yaxis)
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

numParameters = 7;
priorParams = zeros(numParameters,2);

priorParams(1,:) = [0.125,0];% to left for exp pdf, for inverse gamma [25.8169,12.6828];
priorParams(2,:) = [max(xaxis)/2,500];% to left for normal, for gamma [50.788,26.021]; 
priorParams(3,:) = [4,453/4];% gamma with slightly higher variance than mle [6.339,72.07];
priorParams(4,:) = [25.508,11.71]; % gamma
priorParams(5,:) = [19.786,12.543]; % gamma
priorParams(6,:) = [950.71,0.006]; % gamma
priorParams(7,:) = [31.9816,0.1048];% gamma for 1/parameter,  
%                for inverse gamma of parameter 7 [30.5995,9.1263];

lBound = zeros(numParameters,1);
uBound = [Inf,max(xaxis)+50,max(yaxis)+50,Inf,Inf,Inf,Inf]';

numIter = 7e5;
burnIn = 2e5;
skipRate = 100;

posteriorMode = zeros(numChans,numParameters);
posteriorMedian = zeros(numChans,numParameters);
posteriorInterval = zeros(numChans,numParameters,2);
posteriorSamples = zeros(numChans,numParameters,(numIter-burnIn)/skipRate);

for zz=1:numChans
%     display(sprintf('Running Data for Channel %d...',zz));

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
    
    inputParams = params(:,1);inputParams(7) = 1/inputParams(7);
    [loglikelihood] = GetLikelihood(reps,inputParams,vepMagnitude,flashPoints);
    logPrior = 0;
    for ii=[3,4,5,6,7]
       logPrior = logPrior+log(gampdf(params(ii,1),priorParams(ii,1),priorParams(ii,2))); 
    end
    logPrior = logPrior+log(normpdf(params(2,1),priorParams(2,1),priorParams(2,2)));
    logPrior = logPrior+log(exppdf(params(1,1),priorParams(1,1)));
    posteriorProb(1) = loglikelihood+logPrior;
    
    % FOR AUTOMATIC CREATION OF UPDATE MATRIX
    updateParam = 1e-3;
    mu = zeros(numParameters,1);sigma = diag(ones(numParameters,1));
    for ii=1:numParameters
       sigma(ii,ii) = abs(sigma(ii,ii)*params(ii,1)/2); 
    end
    loglambda = log(2.38^2/numParameters);
    updateMu = zeros(numParameters,1);
    
    for ii=[3,4,5,6,7]
        updateMu(ii,1) = gamrnd(priorParams(ii,1),priorParams(ii,2));
    end
    updateMu(1,1) = 1./gamrnd(26.6198,0.0764);
    updateMu(2,1) = normrnd(priorParams(2,1),priorParams(2,2));
    optimalAccept = 0.44;%0.234;
    
    for ii=2:burnIn
        Z = mvnrnd(mu,exp(loglambda).*sigma)';
        pStar = params(:,ii-1)+Z;
        
        if sum(pStar<=lBound) == 0 & sum(pStar>=uBound) == 0
            inputParams = pStar;inputParams(7) = 1/inputParams(7);
            pStarLogLikelihood = GetLikelihood(reps,inputParams,vepMagnitude,flashPoints);
            pStarLogPrior = 0;
            for jj=[3,4,5,6,7]
                pStarLogPrior = pStarLogPrior+log(gampdf(pStar(jj),priorParams(jj,1),priorParams(jj,2)));
            end
            pStarLogPrior = pStarLogPrior+log(exppdf(pStar(1),priorParams(1,1)));
            pStarLogPrior = pStarLogPrior+log(normpdf(pStar(2),priorParams(2,1),priorParams(2,2)));
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
            
            sigma = sigma+updateParam.*((params(:,ii)-updateMu)*...
                (params(:,ii)-updateMu)'-sigma);
            updateMu = updateMu+updateParam.*(params(:,ii)-updateMu);
            loglambda = loglambda+updateParam.*(exp(min(0,logA))-optimalAccept);
            
        else
            params(:,ii) = params(:,ii-1);
            posteriorProb(ii) = posteriorProb(ii-1);
            sigma = sigma+updateParam.*((params(:,ii)-updateMu)*...
                (params(:,ii)-updateMu)'-sigma);
            updateMu = updateMu+updateParam.*(params(:,ii)-updateMu);
            loglambda = loglambda+updateParam.*(-optimalAccept);
        end
%         error = abs(1394-updateMu(2))+abs(419-updateMu(3));
%         scatter(ii,error);%scatter(ii,error);
%         hold on;pause(0.01);
    end
    
    acceptRate = 0;sigma = exp(loglambda).*sigma;
    for ii=burnIn+1:numIter
        Z = mvnrnd(mu,sigma)';
        pStar = params(:,ii-1)+Z;
        
        if sum((pStar-lBound)<=0) == 0 & sum((pStar-uBound)>0) == 0
            inputParams = pStar;inputParams(7) = 1/inputParams(7);
            pStarLogLikelihood = GetLikelihood(reps,inputParams,vepMagnitude,flashPoints);
            pStarLogPrior = 0;
            for jj=[3,4,5,6,7]
                pStarLogPrior = pStarLogPrior+log(gampdf(pStar(jj),priorParams(jj,1),priorParams(jj,2)));
            end
            pStarLogPrior = pStarLogPrior+log(exppdf(pStar(1),priorParams(1,1)));
            pStarLogPrior = pStarLogPrior+log(normpdf(pStar(2),priorParams(2,1),priorParams(2,2)));
            
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
        else
            params(:,ii) = params(:,ii-1);
            posteriorProb(ii) = posteriorProb(ii-1);
        end
    end
    params(7,:) = 1./params(7,:);
    posteriorSamples(zz,:,:) = params(:,burnIn+1:skipRate:end);
    
    figure();
    numRows = ceil(numParameters/2);
    for ii=1:numParameters
       data = squeeze(posteriorSamples(zz,ii,:));
       subplot(numRows,2,ii);histogram(data);
    end
    
    [~,ind] = max(posteriorProb);
    posteriorMode(zz,:) = params(:,ind)';
    posteriorMedian(zz,:) = median(squeeze(posteriorSamples(zz,:,:)),2)';
    
    display(posteriorMode(zz,:));
    
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

