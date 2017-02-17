function [finalParameters,fisherInfo,ninetyfiveErrors] = FitLFPRetinoModel_SGA(Response,xaxis,yaxis,centerVals)
%FitLFPGaussRetinoModel.m
%   Use data from LFP retinotopic mapping experiment to fit a non-linear
%    model of that retinotopy (using peak negativity in window from 50 to
%    120 msec after stimulus presentation and Gaussian likelihood)
%     Use stoachastic gradient ascent ... extremely slow and not as
%     reliable as Levenberg-Marquardt

%Created: 2017/01/18, 24 Cummington Mall, Boston
% Byron Price
%Updated: 2017/02/15
% By: Byron Price

%  data (in response) ~ N(mu,sigma), 
%    where mu = (b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-(distY.^2)./(2*b(3)*b(3)))+b(4));

Bounds = [-10,1000;min(xaxis),max(xaxis);min(yaxis),max(yaxis);1,1000;1,1000;1,1000;0,1000];

numChans = size(Response,1);
reps = size(Response,2);

numParameters = 7;
   
finalParameters = zeros(numChans,numParameters);
fisherInfo = zeros(numChans,numParameters,numParameters);
ninetyfiveErrors = zeros(numChans,numParameters);
numRepeats = 1;numBatches = 2000;
numRedo = 4;
maxITER = 100;
tolerance = 1e-3;
h = ones(numParameters,1)./100;

% parameters are 
%  1) b(1) - rise of map at center
%  2)+3) the position (x and y) of the center of mass
%  4) b(2) - standard deviation or spread of the map in x
%  5) b(3) - standard deviation or spread of map in y
%  6) rise for sigma, from the Gaussian likelihood, at map center
%  7) b(4) - peak negativity at edges of retinotopic region
%    b(1)+b(4) = peak negativity at retinotopic center of mass


% display('Steepest Ascent ...');
for zz=1:numChans
%     display(sprintf('Running Data for Channel %d...',zz));
    flashPoints = centerVals(squeeze(Response(zz,:,1)),:);
    peakNegativity = squeeze(Response(zz,:,2));
    
    finalParams = zeros(numBatches,numParameters);finalLikely = zeros(numBatches,1);
    
    parfor batch = 1:numBatches
        nu = logspace(2,-1,numRedo);
        bigParams = zeros(numRedo,numParameters);bigLikely = zeros(numRedo,1);
        proposal = [1.75,100;4,250;4,250;1.75,150;1.75,150;1.5,150;1.75,200];
        initialParams = zeros(1,numParameters);
        for ii=1:numParameters
            initialParams(1,ii) = gamrnd(proposal(ii,1),proposal(ii,2));
        end
        
        [bigParams(1,:),bigLikely(1)] = GradientAscent(initialParams,numRepeats,...
            numParameters,flashPoints,peakNegativity,maxITER,tolerance,reps,nu(1),1/numRedo,Bounds);

        for ii=2:numRedo
            [~,index] = max(bigLikely(1:ii-1));
            [bigParams(ii,:),bigLikely(ii)] = GradientAscent(bigParams(index,:),...
                numRepeats,numParameters,flashPoints,peakNegativity,maxITER,tolerance,reps,nu(ii),ii/numRedo,Bounds);
        end
        [finalLikely(batch),index] = max(bigLikely);
        finalParams(batch,:) = bigParams(index,:);
    end
    [~,index] = max(finalLikely);plot(1:numBatches,finalLikely);
    finalParameters(zz,:) = finalParams(index,:);
    [fisherInfo(zz,:,:),ninetyfiveErrors(zz,:)] = getFisherInfo(finalParameters(zz,:),...
        numParameters,h,reps,peakNegativity,flashPoints);
end

end

function [params,maxLikely] = GradientAscent(initialParams,numRepeats,numParameters,flashPoints,peakNegativity,maxITER,tolerance,reps,nu,redoRatio,Bounds)
bigParameterVec = zeros(numRepeats,numParameters);
bigLogLikelihoods = zeros(numRepeats,1);
% repeat gradient ascent from a number of different starting
%  positions
initSigma = 100;

for repeats = 1:numRepeats
    h = ones(numParameters,1)./100;
    gradientVec = zeros(1,numParameters);
    parameterVec = zeros(maxITER,numParameters);
    logLikelihood = zeros(maxITER,1);
    
    parameterVec(1,:) = initialParams+normrnd(0,initSigma*(1-redoRatio),[1,numParameters]);
    
    logLikelihood(1) = GetLikelihood(reps,parameterVec(1,:),peakNegativity,flashPoints);
    check = 1;
    iter = 1;
    % for each starting position, do maxITER iterations
    while abs(check) > tolerance && iter < maxITER
        iter = iter+1;
        
        for jj=1:numParameters
            tempParameterVec = parameterVec(iter-1,:);
            tempParameterVec(jj) = tempParameterVec(jj)+h(jj);
            [gradLikelihoodplus] = GetLikelihood(reps,tempParameterVec,peakNegativity,flashPoints);
            
            tempParameterVec = parameterVec(iter-1,:);
            tempParameterVec(jj) = tempParameterVec(jj)-h(jj);
            [gradLikelihoodminus] = GetLikelihood(reps,tempParameterVec,peakNegativity,flashPoints);
            gradientVec(jj) = (gradLikelihoodplus-gradLikelihoodminus)./(2*h(jj));
        end
        
        parameterVec(iter,:) = max(Bounds(:,1)',min(parameterVec(iter-1,:)+nu*gradientVec,Bounds(:,2)'));
        logLikelihood(iter) = GetLikelihood(reps,parameterVec(iter,:),peakNegativity,flashPoints);
        check = diff(logLikelihood(iter-1:iter));
    end
    [~,index] = max(logLikelihood(1:iter));
    bigParameterVec(repeats,:) = parameterVec(index,:);
    bigLogLikelihoods(repeats) = logLikelihood(index);
end
[maxLikely,index] = max(bigLogLikelihoods);
params = bigParameterVec(index,:);

end

function [loglikelihood] = GetLikelihood(reps,parameterVec,peakNegativity,flashPoints)
loglikelihood = 0;
for kk=1:reps
    mu = parameterVec(1)*exp(-((flashPoints(kk,1)-parameterVec(2)).^2)./(2*parameterVec(4).^2)-...
        ((flashPoints(kk,2)-parameterVec(3)).^2)./(2*parameterVec(5).^2))+parameterVec(7);
%*exp(-(distX.^2)./(2*b(2)*b(2))-b(5)*distX*distY/(2*b(2)*b(3))-(distY.^2)./(2*b(3)*b(3)))+parameterVec(9);
    loglikelihood = loglikelihood-(1/2)*log(2*pi*parameterVec(6)*parameterVec(6))-...
        (1/(2*parameterVec(6)*parameterVec(6)))*(peakNegativity(kk)-mu).^2;

end

end

function [fisherInfo,errors] = getFisherInfo(parameters,numParameters,h,reps,peakNegativity,flashPoints)
fisherInfo = zeros(numParameters,numParameters);
errors = zeros(1,numParameters);

% h = h./100;
% the observed fisher information matrix 
for jj=1:numParameters
    for kk=jj:numParameters
        firstParam = jj;secondParam = kk;
        deltaX = h(firstParam);deltaY = h(secondParam);
        
        if firstParam ~= secondParam
            parameterVec = parameters;parameterVec(firstParam) = parameterVec(firstParam)+deltaX;
            parameterVec(secondParam) = parameterVec(secondParam)+deltaY;
            likelyplusplus = GetLikelihood(reps,parameterVec,peakNegativity,flashPoints);
            
            parameterVec = parameters;parameterVec(firstParam) = parameterVec(firstParam)+deltaX;
            parameterVec(secondParam) = parameterVec(secondParam)-deltaY;
            likelyplusminus = GetLikelihood(reps,parameterVec,peakNegativity,flashPoints);
            
            parameterVec = parameters;parameterVec(firstParam) = parameterVec(firstParam)-deltaX;
            parameterVec(secondParam) = parameterVec(secondParam)+deltaY;
            likelyminusplus = GetLikelihood(reps,parameterVec,peakNegativity,flashPoints);
            
            parameterVec = parameters;parameterVec(firstParam) = parameterVec(firstParam)-deltaX;
            parameterVec(secondParam) = parameterVec(secondParam)-deltaY;
            likelyminusminus = GetLikelihood(reps,parameterVec,peakNegativity,flashPoints);
            
            fisherInfo(jj,kk) = -(likelyplusplus-likelyplusminus-likelyminusplus+likelyminusminus)./(4*deltaX*deltaY);
        else
            likely = GetLikelihood(reps,parameters,peakNegativity,flashPoints);
            
            parameterVec = parameters;parameterVec(firstParam) = parameterVec(firstParam)-deltaX;
            likelyminus = GetLikelihood(reps,parameterVec,peakNegativity,flashPoints);
            
            parameterVec = parameters;parameterVec(firstParam) = parameterVec(firstParam)+deltaX;
            likelyplus = GetLikelihood(reps,parameterVec,peakNegativity,flashPoints);
            
            fisherInfo(jj,kk) = -(likelyminus-2*likely+likelyplus)./(deltaX*deltaX);
        end
    end
end

transpose = fisherInfo';
for ii=1:numParameters
    for jj=1:ii-1
        fisherInfo(ii,jj) = transpose(ii,jj);
    end
end

inverseFisherInfo = inv(fisherInfo);
for ii=1:numParameters
    errors(ii) = sqrt(inverseFisherInfo(ii,ii));
end

if isreal(errors) == 0
    temp = sqrt(errors.*conj(errors));
    errors = 1.96.*temp;
elseif isreal(errors) == 1
    errors = 1.96.*errors;
end
end


