function [finalParameters,covariance] = BayesianFitLFPModel(Response,xaxis,yaxis,centerVals)
% BayesianFitLFPModel.m

numChans = size(Response,1);
reps = size(Response,2);

numParameters = 7;
   
finalParameters = zeros(numChans,numParameters);
covariance = zeros(numChans,numParameters,numParameters);

priorParams = [1.75,100;4,300;4,250;1.75,200;1.75,200;1.5,150;1.75,150];
% proposal = [150,100;1000,500;700,500;300,200;250,300;200,150;200,150];

sigma = 1;

for zz=1:numChans
    flashPoints = centerVals(squeeze(Response(zz,:,1)),:);
    peakNegativity = log(abs(squeeze(Response(zz,:,2)));
    N = 3e6; burnIn = 2e6;
    x = zeros(N,numParameters);
    
    [finalParams,finalLikelihood] = GradientAscentInit(flashPoints,peakNegativity,reps);
    x(1,:) = finalParams;
    likelihoodXprev = finalLikelihood;
    
    for ii=2:N
%         priorXstar = 1;priorXprev = 1;
        xStar = x(ii-1,:)+normrnd(0,sigma,[1,numParameters]);

%         for jj=1:numParameters
% %             priorXstar = priorXstar*gampdf(xStar(jj),priorParams(jj,1),priorParams(jj,2));
% %             priorXprev = priorXprev*gampdf(x(ii-1,jj),priorParams(jj,1),priorParams(jj,2));
%         end
        if sum(xStar<0) > 0
            x(ii,:) = x(ii-1,:);
        else
            likelihoodXstar = GetLikelihood(reps,xStar,peakNegativity,flashPoints);
            %         A = min(1,(likelihoodXstar*priorXstar)./(likelihoodXprev*priorXprev));
%             A = min(1,likelihoodXprev/likelihoodXstar);
            logA = likelihoodXstar-likelihoodXprev;
            if log(rand) < logA
                x(ii,:) = xStar;
                likelihoodXprev = likelihoodXstar;
            else
                x(ii,:) = x(ii-1,:);
            end
        end
    end
    figure();
    for ii=1:numParameters
       subplot(4,2,ii);histogram(x(burnIn:end,ii)); 
    end
    finalParameters(zz,:) = median(x(burnIn:end,:),1);
    mode(x(burnIn:end,:))
    covariance(zz,:,:) = cov(x(burnIn:end,:));
end

end

function [loglikelihood] = GetLikelihood(reps,parameterVec,peakNegativity,flashPoints)
loglikelihood = 0;
for kk=1:reps
    distX = flashPoints(kk,1)-parameterVec(2);
    distY = flashPoints(kk,2)-parameterVec(3);
    b = [parameterVec(1),parameterVec(4),parameterVec(5),parameterVec(7)];
    mu = b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-(distY.^2)./(2*b(3)*b(3)))+b(4);
    stdev = parameterVec(6);
    loglikelihood = loglikelihood-(1/2)*log(2*pi*stdev*stdev)-(1/(2*stdev*stdev))*(peakNegativity(kk)-mu).^2;
end
end

function [finalParams,finalLikelihood] = GradientAscentInit(flashPoints,peakNegativity,reps)
%GradientAscentInit

numParameters = 7;
maxITER = 1000;
tolerance = 1e-2;

h = ones(numParameters,1)./100;

% repeat gradient ascent from a number of different starting
%  positions
gradientVec = zeros(1,numParameters);
parameterVec = zeros(maxITER,numParameters);
logLikelihood = zeros(maxITER,1);

priorParams = [1.75,100;4,300;4,250;1.75,200;1.75,200;1.5,150;1.75,150];
for jj=1:numParameters
    parameterVec(1,jj) = gamrnd(priorParams(jj,1),priorParams(jj,2));
end


logLikelihood(1) = GetLikelihood(reps,parameterVec(1,:),peakNegativity,flashPoints);
check = 1;
iter = 1;

while check > tolerance && iter < maxITER
    iter = iter+1;
    % calculate likelihood at the current position in parameter
    %  space
    %
    %             % check if a random jump along one dimension improves the
    %             %  likelihood
    %             temp = parameterVec(iter-1,:);
    %             randInd = random('Discrete Uniform',numParameters,1);
    %             temp(randInd) = Bounds(randInd,1)+(Bounds(randInd,2)-Bounds(randInd,1)).*rand;
    %
    %             [likely] = GetLikelihood(reps,temp,peakNegativity,flashPoints);
    %             if likely > logLikelihood(iter-1) %|| mod(iter,200) == 0
    %                 parameterVec(iter-1,:) = temp';
    %                 logLikelihood(iter-1) = likely;
    %             end
    
    % calculate the gradient by calculating the likelihood
    %  after moving over a small step h along each
    %  dimension
    for jj=1:numParameters
        tempParameterVec = parameterVec(iter-1,:);
        tempParameterVec(jj) = tempParameterVec(jj)+h(jj);
        [gradLikelihoodplus] = GetLikelihood(reps,tempParameterVec,peakNegativity,flashPoints);
        
        tempParameterVec = parameterVec(iter-1,:);
        tempParameterVec(jj) = tempParameterVec(jj)-h(jj);
        [gradLikelihoodminus] = GetLikelihood(reps,tempParameterVec,peakNegativity,flashPoints);
        gradientVec(jj) = (gradLikelihoodplus-gradLikelihoodminus)./(2*h(jj));
    end
    
    % line search to get distance to move along gradient
    alpha = [0,1e-6,1e-4,1e-2,1e-1];
    lineSearchLikelihoods = zeros(length(alpha),1);
    lineSearchLikelihoods(1) = logLikelihood(iter-1);
    
    for ii=2:length(alpha)
        tempParameterVec = parameterVec(iter-1,:)'+gradientVec.*alpha(ii);
        [lineSearchLikelihoods(ii)] = GetLikelihood(reps,tempParameterVec,peakNegativity,flashPoints);
    end
    
    newInds = find(imag(lineSearchLikelihoods)==0);
    newValues = zeros(length(newInds),1);
    for ii=1:length(newInds)
        newValues(ii) = lineSearchLikelihoods(newInds(ii));
    end
    [~,ind] = max(newValues);
    ind = newInds(ind);
    
    logLikelihood(iter) = lineSearchLikelihoods(ind);
    check = diff(logLikelihood(iter-1:iter));
    parameterVec(iter,:) = parameterVec(iter-1,:)+alpha(ind).*gradientVec;
end

finalParams = parameterVec(iter,:);
finalLikelihood = logLikelihood(iter);
end