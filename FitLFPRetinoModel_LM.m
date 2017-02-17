function [finalParameters,fisherInfo,ninetyfiveErrors] = FitLFPRetinoModel_LM(Response,xaxis,yaxis,centerVals)
%FitLFPRetinoModel_LM.m
%   Use data from LFP retinotopic mapping experiment to fit a non-linear
%    model of that retinotopy (using peak negativity in window from 50 to
%    120 msec after stimulus presentation and Gaussian likelihood)
%   Levenberg-Marquardt algorithm ... far superior to gradient ascent for
%    this data (faster and more reliable)

%Created: 2017/02/16, 24 Cummington Mall, Boston
% Byron Price
%Updated: 2017/02/16
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
numRepeats = 500;
maxITER = 50;
tolerance = 1e-3;


for zz=1:numChans
%     display(sprintf('Running Data for Channel %d...',zz));
    flashPoints = centerVals(squeeze(Response(zz,:,1)),:);
    peakNegativity = squeeze(Response(zz,:,2))';
    h = ones(numParameters,1)./100;
    bigParameterVec = zeros(numParameters,numRepeats);
    bigResiduals = zeros(numRepeats,1);
    % repeat gradient ascent from a number of different starting
    %  positions

    parfor repeats = 1:numRepeats
        parameterVec = zeros(numParameters,maxITER);
        squaredResiduals = zeros(maxITER,1);

        proposal = [1.75,100;4,250;4,250;1.75,150;1.75,150;1.5,150;1.75,200];
        for ii=1:numParameters 
              parameterVec(ii,1) = gamrnd(proposal(ii,1),proposal(ii,2));
        end
        
        [yhat] = Getyhat(reps,parameterVec(:,1),flashPoints);
        squaredResiduals(1) = sum((peakNegativity-yhat).^2);
        
        check = 1;
        iter = 1;
        lambda = 100;
        % for each starting position, do maxITER iterations
        while abs(check) > tolerance && iter < maxITER
            [Jacobian] = GetJacobian(reps,parameterVec(:,iter),flashPoints,numParameters,h,yhat);
            H = Jacobian'*Jacobian;
            update = pinv(H+lambda.*diag(diag(H)))*Jacobian'*(peakNegativity-yhat);
            tempParams = parameterVec(:,iter)+update;
            
            tempParams = max(Bounds(:,1),min(tempParams,Bounds(:,2)));
            
            [tempYhat] = Getyhat(reps,tempParams,flashPoints);
            squaredResiduals(iter+1) = sum((peakNegativity-tempYhat).^2);
            check = diff(squaredResiduals(iter:iter+1));
            if check >= 0
                parameterVec(:,iter+1) = parameterVec(:,iter);
                lambda = lambda*10;
                check = 1;
                squaredResiduals(iter+1) = squaredResiduals(iter);
            else
                parameterVec(:,iter+1) = tempParams;
                yhat = tempYhat;
                lambda = lambda/10;
            end
            iter = iter+1;
        end
        [bigResiduals(repeats),index] = min(squaredResiduals(1:iter));
        bigParameterVec(:,repeats) = parameterVec(:,index);
    end
    [~,index] = min(bigResiduals);
    finalParameters(zz,:) = bigParameterVec(:,index)';
    [yhat] = Getyhat(reps,finalParameters(zz,:),flashPoints);
    finalParameters(zz,6) = std(peakNegativity-yhat);

    parameterVec = zeros(maxITER,numParameters);gradientVec = zeros(1,numParameters);
    logLikelihood = zeros(maxITER,1);
    
    parameterVec(1,:) = finalParameters(zz,:);
    [logLikelihood(1)] = GetLikelihood(reps,parameterVec(1,:),peakNegativity,flashPoints);
    check = 1;iter = 1;lambda = 10;
    while abs(check) > tolerance && iter < maxITER
        for jj=1:numParameters
            tempParameterVec = parameterVec(iter,:);
            tempParameterVec(jj) = tempParameterVec(jj)+h(jj);
            [gradLikelihoodplus] = GetLikelihood(reps,tempParameterVec,peakNegativity,flashPoints);
            
            gradientVec(jj) = (gradLikelihoodplus-logLikelihood(iter))./h(jj);
        end
        tempParams = parameterVec(iter,:)+lambda*gradientVec;
        tempLikelihood = GetLikelihood(reps,tempParams,peakNegativity,flashPoints);
        check = logLikelihood(iter)-tempLikelihood;
        
        if check <= 0
            lambda = lambda/10;
            parameterVec(iter+1,:) = parameterVec(iter,:);
            logLikelihood(iter+1) = logLikelihood(iter);
        else
            parameterVec(iter+1,:) = tempParams;
            logLikelihood(iter+1,:) = tempLikelihood;
        end
        iter = iter+1;
    end
    [~,index] = max(logLikelihood(1:iter));
    finalParameters(zz,:) = parameterVec(index,:);
    [fisherInfo(zz,:,:),ninetyfiveErrors(zz,:)] = getFisherInfo(finalParameters(zz,:),numParameters,h,reps,peakNegativity,flashPoints);
    display(zz);
    display(finalParameters(zz,:));
    display(ninetyfiveErrors(zz,:));
end
end

function [Jacobian] = GetJacobian(reps,parameterVec,flashPoints,numParameters,h,yhat)
Jacobian = zeros(reps,numParameters);
for kk=1:reps
    for jj=1:numParameters
       tempParams = parameterVec;tempParams(jj) = tempParams(jj)+h(jj);
       mu = tempParams(1)*exp(-((flashPoints(kk,1)-tempParams(2)).^2)./(2*tempParams(4).^2)-...
        ((flashPoints(kk,2)-tempParams(3)).^2)./(2*tempParams(5).^2))+tempParams(7);
       Jacobian(kk,jj) = (mu-yhat(kk))/h(jj);
    end
end
end

function [yhat] = Getyhat(reps,parameterVec,flashPoints)
yhat = zeros(reps,1);
for kk=1:reps
    yhat(kk) = parameterVec(1)*exp(-((flashPoints(kk,1)-parameterVec(2)).^2)./(2*parameterVec(4).^2)-...
        ((flashPoints(kk,2)-parameterVec(3)).^2)./(2*parameterVec(5).^2))+parameterVec(7);

end
end

function [loglikelihood] = GetLikelihood(reps,parameterVec,peakNegativity,flashPoints)
loglikelihood = 0;
for kk=1:reps
    mu = parameterVec(1)*exp(-((flashPoints(kk,1)-parameterVec(2)).^2)./(2*parameterVec(4).^2)-...
        ((flashPoints(kk,2)-parameterVec(3)).^2)./(2*parameterVec(5).^2))+parameterVec(7);
    stdev = parameterVec(6);%*exp(-(distX.^2)./(2*b(2)*b(2))-b(5)*distX*distY/(2*b(2)*b(3))-(distY.^2)./(2*b(3)*b(3)))+parameterVec(9);
    loglikelihood = loglikelihood-(1/2)*log(2*pi*stdev*stdev)-(1/(2*stdev*stdev))*(peakNegativity(kk)-mu).^2;
%     summation = summation+(peakNegativity(kk)-mu).^2;
end

% loglikelihood = (-reps/2)*log(2*pi*parameterVec(6)*parameterVec(6))-...
%     (1/(2*parameterVec(6)*parameterVec(6)))*summation;
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


