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
%Updated: 2017/06/14
% By: Byron Price

%  model has 7 parameters, defined by vector p
%  data ~ log-logistic(log mean, log scale), 
%    where the log mean of the data is as follows
%     log mean = p(1)*exp(-(xpos-p(2)).^2./(2*p(4)*p(4))-
%        (ypos-p(3)).^2./(2*p(5)*p(5)))+p(6);
%
%  and p(7) = log scale

numChans = size(Response,1);

maxNumCompThreads(4);
gcp = parpool(min(numChans,4));

numParameters = 7;
numIter = 6e5;
burnIn = 1e5;
skipRate = 500;

posteriorMode = zeros(numChans,numParameters);
posteriorMean = zeros(numChans,numParameters);
posteriorInterval = zeros(numChans,numParameters,2);
posteriorSamples = zeros(numChans,numParameters,(numIter-burnIn)/skipRate);

parfor zz=1:numChans
    fprintf('Running Data for Channel %d...',zz);
    priorParams = zeros(numParameters,2);

    priorParams(1,:) = [0,50];% [0.125,0] for exp pdf, [9.099,25.089] for gamma with log
                             % for inverse gamma without log [25.8169,12.6828];
    priorParams(2,:) = [max(xaxis)/2,1000];% to left for normal, for gamma [50.788,26.021];
    priorParams(3,:) = [500,1000];% normal ... mle [6.339,72.07];
    priorParams(4,:) = [25.508,11.71]; % gamma
    priorParams(5,:) = [19.786,12.543]; % gamma
    priorParams(6,:) = [38.25,8.74]; % gamma with log conversion
    priorParams(7,:) = [30,0.112];% gamma for 1/parameter,
    %                for increased variance [30,0.112]  ... original
    %                   [31.98,0.105]
%     priorParams(8,:) = [1e-3,1e3]; % gamma, for precision on first parameter

    lBound = zeros(numParameters,1);lBound(1) = 0;
    lBound(2) = min(xaxis)-50;lBound(3) = min(yaxis)-50;
    uBound = [Inf,max(xaxis)+50,max(yaxis)+50,Inf,Inf,Inf,Inf]';


    Data = Response{zz};
    reps = size(Data,1);
    
    flashPoints = Data(:,1:2);
    vepMagnitude = Data(:,3);
    
    posteriorProb = zeros(numIter,1);
    params = zeros(numParameters,numIter);
    
    numStarts = 5000;
    maxPost = zeros(numStarts,1);
    paramSet = zeros(numParameters,numStarts);
    for jj=1:numStarts
        
        for ii=4:7
            paramSet(ii,jj) = gamrnd(priorParams(ii,1),priorParams(ii,2));
        end
        paramSet(1,jj) = gamrnd(9.1,25.09); 
        paramSet(2,jj) = normrnd(priorParams(2,1),priorParams(2,2));
        paramSet(3,jj) = normrnd(priorParams(3,1),priorParams(3,2));
%         params(8,1) = log(1.5);
        paramSet(:,jj) = min(max(paramSet(:,jj),lBound),uBound);
        
        inputParams = paramSet(:,jj);inputParams(7) = 1/inputParams(7);
        [loglikelihood] = GetLikelihood(reps,inputParams,vepMagnitude,flashPoints);
        
        logPrior = sum(log(gampdf(paramSet(4:7,jj),priorParams(4:7,1),priorParams(4:7,2))))+...
               sum(log(normpdf(paramSet(2:3,jj),priorParams(2:3,1),priorParams(2:3,2))))+...
               log(exppdf(paramSet(1,jj),priorParams(1,2)));
        maxPost(jj) = loglikelihood+logPrior;
        
        
        updateParam = 0.1;
        mu = zeros(numParameters,1);sigma = eye(numParameters);
        for ii=1:numParameters
            if ii<=3
                sigma(ii,ii) = priorParams(ii,2)^2/2;
            else
                sigma(ii,ii) = priorParams(ii,1)*(priorParams(ii,2)^2)/2;
            end
        end
%         halfSigma = cholcov(sigma);identity = eye(numParameters);
        loglambda = log(2.38^2/numParameters);
        updateMu = zeros(numParameters,1);
        
        for ii=4:7
            updateMu(ii,1) = gamrnd(priorParams(ii,1),priorParams(ii,2));
        end
        updateMu(1,1) = gamrnd(9.1,25.09);
        updateMu(2,1) = normrnd(priorParams(2,1),priorParams(2,2));
        updateMu(3,1) = normrnd(priorParams(3,1),priorParams(3,2));
        updateMu = min(max(updateMu,lBound),uBound);
        optimalAccept = 0.234;%0.44
        %     scatter(1,posteriorProb(1));hold on;pause(0.5);
        
        for ii=1:100
            pStar = paramSet(:,jj)+mvnrnd(mu,exp(loglambda)*sigma)';
            
            if sum(pStar<lBound) == 0 && sum(pStar>=uBound) == 0
                inputParams = pStar;inputParams(7) = 1/inputParams(7);
                pStarLogLikelihood = GetLikelihood(reps,inputParams,vepMagnitude,flashPoints);
                pStarLogPrior = sum(log(gampdf(pStar(4:7),priorParams(4:7,1),priorParams(4:7,2))))+...
                    sum(log(normpdf(pStar(2:3),priorParams(2:3,1),priorParams(2:3,2))))+...
                    log(exppdf(pStar(1),priorParams(1,2)));
                
                logA = (pStarLogLikelihood+pStarLogPrior)-maxPost(jj);
                
                if logA>0
                    paramSet(:,jj) = pStar;
                    maxPost(jj) = pStarLogLikelihood+pStarLogPrior;
                end
                
                meanSubtract = paramSet(:,jj)-updateMu;
                updateMu = updateMu+updateParam.*meanSubtract;
                sigma = sigma+updateParam.*(meanSubtract*meanSubtract'-sigma);
                loglambda = loglambda+updateParam.*(exp(min(0,logA))-optimalAccept);
                
            else
                loglambda = loglambda+updateParam.*(-optimalAccept);
            end
            %         error = abs(1394-updateMu(2))+abs(419-updateMu(3));
%                     scatter(ii,maxPost(jj));%scatter(ii,error);
%                     hold on;pause(0.01);
        end
    end
    [~,ind] = max(maxPost);
    posteriorProb(1) = maxPost(ind);
    params(:,1) = paramSet(:,ind);

    % FOR AUTOMATIC CREATION OF PROPOSAL DISTRIBUTION
    updateParam = 1e-2;
    mu = zeros(numParameters,1);sigma = eye(numParameters);
    for ii=1:numParameters
       if ii<=3
           sigma(ii,ii) = priorParams(ii,2)^2/2;
       else
           sigma(ii,ii) = priorParams(ii,1)*(priorParams(ii,2)^2)/2;
       end
    end
    halfSigma = cholcov(sigma);identity = eye(numParameters);
    loglambda = log(2.38^2/numParameters);
    updateMu = params(:,1)+mvnrnd(mu,exp(loglambda)*sigma)';
    updateMu = min(max(updateMu,lBound),uBound);
    optimalAccept = 0.234;%0.44
%     scatter(1,posteriorProb(1));hold on;pause(0.5);

    for ii=2:burnIn
        pStar = params(:,ii-1)+mvnrnd(mu,exp(loglambda)*sigma)';
        
        if sum(pStar<lBound) == 0 && sum(pStar>=uBound) == 0
            inputParams = pStar;inputParams(7) = 1/inputParams(7);
            pStarLogLikelihood = GetLikelihood(reps,inputParams,vepMagnitude,flashPoints);
            pStarLogPrior = sum(log(gampdf(pStar(4:7),priorParams(4:7,1),priorParams(4:7,2))))+...
                sum(log(normpdf(pStar(2:3),priorParams(2:3,1),priorParams(2:3,2))))+...
                log(exppdf(pStar(1),priorParams(1,2)));
            
            logA = (pStarLogLikelihood+pStarLogPrior)-posteriorProb(ii-1);
            
            if log(rand) < logA
                params(:,ii) = pStar;
                posteriorProb(ii) = pStarLogLikelihood+pStarLogPrior;
            else
                params(:,ii) = params(:,ii-1);
                posteriorProb(ii) = posteriorProb(ii-1);
            end
            
            if mod(ii,100)==0
                meanSubtract = params(:,ii)-updateMu;
                updateMu = updateMu+updateParam.*meanSubtract;
                halfSigma = halfSigma+updateParam.*(triu((inv(halfSigma))*(halfSigma'*halfSigma+meanSubtract*...
                    meanSubtract')*((inv(halfSigma))')-identity)-halfSigma);
                
                sigma = halfSigma'*halfSigma;
                %                 sigma = sigma+updateParam(ii).*(meanSubtract*meanSubtract'-sigma);
                loglambda = loglambda+updateParam.*(exp(min(0,logA))-optimalAccept);
            end
        else
            params(:,ii) = params(:,ii-1);
            posteriorProb(ii) = posteriorProb(ii-1);
            if mod(ii,100)==0
                loglambda = loglambda+updateParam.*(-optimalAccept);
            end
        end
%         error = abs(1394-updateMu(2))+abs(419-updateMu(3));
%         scatter(ii,posteriorProb(ii));%scatter(ii,error);
%         hold on;pause(0.01);
    end
    
%     [V,D] = eig(cov(params(:,1e4:100:burnIn)'));
%     W = V*sqrtm(D);
%     eigenvals = diag(W'*W);
%     
%     
%     tempW = [];
%     tempEigs = [];
%     for jj=1:numParameters
%         if eigenvals(jj) > 1e-5
%             tempW = [tempW,W(:,jj)];
%             tempEigs = [tempEigs,eigenvals(jj)];
%         end
%     end
%     
%     q = length(tempEigs);
%     if q <= 2
%         [V,D] = eig(sigma);
%         W = fliplr(V*sqrtm(D));
%         eigenvals = fliplr(diag(W'*W));
%         q = length(eigenvals);
%         fprintf('Error\n');
%     else
%         W = fliplr(tempW);
%         eigenvals = fliplr(tempEigs);
%         q = length(eigenvals);
%     end
% %     p = ones(q,1)./q;
%     loglambda = loglambda.*ones(q,1);
%     updateParam = 1e-3;
    sigma = cov(params(:,1e4:100:burnIn)');
    lambda = 2.38^2/numParameters;
    acceptRate = 0;
    for ii=burnIn+1:numIter
%         index = random('Discrete Uniform',q);
%         lambda = loglambda(index);
%         stdev = sqrt(exp(lambda).*eigenvals(index));
%         pStar = params(:,ii-1)+W(:,index)*normrnd(0,stdev);
        pStar = params(:,ii-1)+mvnrnd(mu,lambda*sigma)';
        
        if sum(pStar<lBound) == 0 && sum(pStar>=uBound) == 0
            inputParams = pStar;inputParams(7) = 1/inputParams(7);
            pStarLogLikelihood = GetLikelihood(reps,inputParams,vepMagnitude,flashPoints);
            pStarLogPrior = sum(log(gampdf(pStar(4:7),priorParams(4:7,1),priorParams(4:7,2))))+...
                sum(log(normpdf(pStar(2:3),priorParams(2:3,1),priorParams(2:3,2))))+...
                log(exppdf(pStar(1),priorParams(1,2)));
            
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
    
%     test = autocorr(squeeze(posteriorSamples(zz,1,:)));
%     if test(2) > 0.2
%         fprintf('Auto-correlated Samples\n');
%     end
    
%     figure();
%     numRows = ceil(numParameters/2);
%     for ii=1:numParameters
%        data = squeeze(posteriorSamples(zz,ii,:));
%        subplot(numRows,2,ii);histogram(data);
%     end
%     
    [~,ind] = max(posteriorProb);
    posteriorMode(zz,:) = params(:,ind)';
    posteriorMean(zz,:) = mean(squeeze(posteriorSamples(zz,:,:)),2)';
    
%     display(posteriorMean(zz,:));
    
    alpha = 0.05;
    posteriorInterval(zz,:,:) = quantile(squeeze(posteriorSamples(zz,:,:)),[alpha/2,1-alpha/2],2);
end
delete(gcp);
end

% log-logistic likelihood for non-linear 2D Gaussian function
function [loglikelihood] = GetLikelihood(reps,parameterVec,vepMagnitude,flashPoints)
loglikelihood = zeros(reps,1);
for kk=1:reps
    mu = log(parameterVec(1)*exp(-((flashPoints(kk,1)-parameterVec(2)).^2)./(2*parameterVec(4).^2)-...
        ((flashPoints(kk,2)-parameterVec(3)).^2)./(2*parameterVec(5).^2))+parameterVec(6));
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

