function [finalParameters,fisherInfo,ninetyfiveErrors] = FitLFPRetinoModel_SGA(Response,xaxis,yaxis,centerVals)
%FitLFPRetinoModel_LM.m
%   Use data from LFP retinotopic mapping experiment to fit a non-linear
%    model of that retinotopy (using peak negativity in window from 50 to
%    120 msec after stimulus presentation and Gaussian likelihood)
%      Stochastic gradient ascent with the Jacobian matrix ... not as
%      effective as Levenberg-Marquardt

%Created: 2017/02/16, 24 Cummington Mall, Boston
% Byron Price
%Updated: 2017/02/21
% By: Byron Price

%  data (in response) ~ N(mu,sigma), 
%    where mu = (b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-(distY.^2)./(2*b(3)*b(3)))+b(4));

Bounds = [0,1000;min(xaxis),max(xaxis);min(yaxis),max(yaxis);1,1000;1,1000;1,1000;0,1000];
numChans = size(Response,1);
reps = size(Response,2);

numParameters = 7;
   
finalParameters = zeros(numChans,numParameters);
fisherInfo = zeros(numChans,numParameters,numParameters);
ninetyfiveErrors = zeros(numChans,numParameters);
numRepeats = 1000;
maxITER = 50;
tolerance = 1e-3;


for zz=1:numChans
%     display(sprintf('Running Data for Channel %d...',zz));
    flashPoints = centerVals(squeeze(Response(zz,:,1)),:);
    peakNegativity = squeeze(Response(zz,:,2))';
    h = ones(numParameters,1)./100;
    bigParameterVec = zeros(numParameters,numRepeats);
    bigLikelihoods = zeros(numRepeats,1);
    % repeat gradient ascent from a number of different starting
    %  positions

    for repeats = 1:numRepeats
        parameterVec = zeros(numParameters,maxITER);
        maxLikelihood = zeros(maxITER,1);
        logLikelihoods = zeros(reps,maxITER);

        proposal = [1.75,100;4,250;4,250;1.75,150;1.75,150;1.5,150;1.75,200];
        for ii=1:numParameters 
              parameterVec(ii,1) = gamrnd(proposal(ii,1),proposal(ii,2));
        end
        
        [logLikelihoods(:,1)] = GetLikelihood(reps,parameterVec(:,1),peakNegativity,flashPoints);
        maxLikelihood(1) = sum(logLikelihoods(:,1));

        check = 1;
        iter = 1;
        lambda = 100;
        % for each starting position, do maxITER iterations
        while abs(check) > tolerance && iter < maxITER
            [Jacobian,tempLikelihood] = GetJacobian(reps,parameterVec(:,iter),peakNegativity,flashPoints,numParameters,h,logLikelihoods(:,iter));
%             H = Jacobian'*Jacobian;
%             update = pinv(H+lambda.*diag(diag(H)))*Jacobian'*(logLikelihoods(iter)-tempLikelihood);
            update = Jacobian'*(tempLikelihood-logLikelihoods(:,iter));
            tempParams = parameterVec(:,iter)+lambda.*update;
            
            tempParams = max(Bounds(:,1),min(tempParams,Bounds(:,2)));
            
            [logLikelihoods(:,iter+1)] = GetLikelihood(reps,tempParams,peakNegativity,flashPoints);
            maxLikelihood(iter+1) = sum(logLikelihoods(:,iter+1));
            
            check = diff(maxLikelihood(iter:iter+1));
            if check <= 0
                parameterVec(:,iter+1) = parameterVec(:,iter);
                lambda = max(lambda/10,1e-12);
                check = 1;
                logLikelihoods(:,iter+1) = logLikelihoods(:,iter);
                maxLikelihood(iter+1) = maxLikelihood(iter);
            else
                parameterVec(:,iter+1) = tempParams;
                lambda = min(lambda*10,1e12);
            end
            iter = iter+1;
        end
        [bigLikelihoods(repeats),index] = max(maxLikelihood(1:iter));
        bigParameterVec(:,repeats) = parameterVec(:,index);
    end
    [~,index] = max(bigLikelihoods);
    finalParameters(zz,:) = bigParameterVec(:,index)'
    
    [fisherInfo(zz,:,:),ninetyfiveErrors(zz,:)] = getFisherInfo(finalParameters(zz,:),numParameters,h,reps,peakNegativity,flashPoints);
    display(zz);
    display(finalParameters(zz,:));
    display(ninetyfiveErrors(zz,:));
end
end

function [Jacobian,likelihood] = GetJacobian(reps,parameterVec,peakNegativity,flashPoints,numParameters,h,prevLikelihood)
Jacobian = zeros(reps,numParameters);
likelihood = zeros(reps,1);
for kk=1:reps
    for jj=1:numParameters
       tempParams = parameterVec;tempParams(jj) = tempParams(jj)+h(jj);
       mu = tempParams(1)*exp(-((flashPoints(kk,1)-tempParams(2)).^2)./(2*tempParams(4).^2)-...
        ((flashPoints(kk,2)-tempParams(3)).^2)./(2*tempParams(5).^2))+tempParams(7);
       likelihood(kk) = -(1/2)*log(2*pi*tempParams(6)*tempParams(6))-...
        (1/(2*tempParams(6)*tempParams(6)))*(peakNegativity(kk)-mu).^2;
       Jacobian(kk,jj) = (likelihood(kk)-prevLikelihood(kk))/h(jj);
    end
end
end


function [loglikelihood] = GetLikelihood(reps,parameterVec,peakNegativity,flashPoints)
loglikelihood = zeros(reps,1);
for kk=1:reps
    mu = parameterVec(1)*exp(-((flashPoints(kk,1)-parameterVec(2)).^2)./(2*parameterVec(4).^2)-...
        ((flashPoints(kk,2)-parameterVec(3)).^2)./(2*parameterVec(5).^2))+parameterVec(7);
      %*exp(-(distX.^2)./(2*b(2)*b(2))-b(5)*distX*distY/(2*b(2)*b(3))-(distY.^2)./(2*b(3)*b(3)))+parameterVec(9);
    loglikelihood(kk) = -(1/2)*log(2*pi*parameterVec(6)*parameterVec(6))-...
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


