function [finalParameters,fisherInfo,ninetyfiveErrors] = FitLFPRetinoModel_LM(Response,xaxis,yaxis,centerVals)
%FitLFPRetinoModel_LM.m
%   Use data from LFP retinotopic mapping experiment to fit a non-linear
%    model of that retinotopy (using peak negativity in window from 50 to
%    120 msec after stimulus presentation and Gaussian likelihood)
%   Levenberg-Marquardt algorithm

%Created: 2017/02/16, 24 Cummington Mall, Boston
% Byron Price
%Updated: 2017/02/16
% By: Byron Price

%  data (in response) ~ N(mu,sigma), 
%    where mu = (b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-(distY.^2)./(2*b(3)*b(3)))+b(4));

numChans = size(Response,1);
reps = size(Response,2);

numParameters = 7;
   
finalParameters = zeros(numChans,numParameters);
fisherInfo = zeros(numChans,numParameters,numParameters);
ninetyfiveErrors = zeros(numChans,numParameters);
numRepeats = 50;
maxITER = 2000;
tolerance = 1e-2;


for zz=1:numChans
%     display(sprintf('Running Data for Channel %d...',zz));
    flashPoints = centerVals(squeeze(Response(zz,:,1)),:);
    peakNegativity = squeeze(Response(zz,:,2));
    h = ones(numParameters,1)./10;
    bigParameterVec = zeros(numRepeats,numParameters);
    bigLogLikelihoods = zeros(numRepeats,1);
    % repeat gradient ascent from a number of different starting
    %  positions

    parfor repeats = 1:numRepeats
        gradientVec = zeros(numParameters,1);
        parameterVec = zeros(numParameters,maxITER);
        logLikelihood = zeros(maxITER,1);

        proposal = [1.75,100;4,250;4,250;1.75,150;1.75,150;1.5,150;1.75,200];
        for ii=1:numParameters 
              parameterVec(ii,1) = gamrnd(proposal(ii,1),proposal(ii,2));
        end
        logLikelihood(1) = GetLikelihood(reps,parameterVec(:,1),peakNegativity,flashPoints);
        check = 1;
        iter = 1;
        lambda = 100;

        % for each starting position, do maxITER iterations
        while check > tolerance && iter < maxITER
            % calculate the gradient by calculating the likelihood
            %  after moving over a small step h along each
            %  dimension
            for jj=1:numParameters
                tempParameterVec = parameterVec(:,iter);
                tempParameterVec(jj) = tempParameterVec(jj)+h(jj);
                [gradLikelihoodplus] = GetLikelihood(reps,tempParameterVec,peakNegativity,flashPoints);
                
                tempParameterVec = parameterVec(:,iter);
                tempParameterVec(jj) = tempParameterVec(jj)-h(jj);
                [gradLikelihoodminus] = GetLikelihood(reps,tempParameterVec,peakNegativity,flashPoints);
                gradientVec(jj) = (gradLikelihoodplus-gradLikelihoodminus)./(2*h(jj));
            end
            
            H = gradientVec*gradientVec';
            diagH = diag(diag(H));
%             tempParams = parameterVec(:,iter)+inv(H+lambda.*diagH)*gradientVec;
%             tempParams = parameterVec(:,iter)+inv(H+lambda.*eye(numParameters))*gradientVec;
            tempParams = parameterVec(:,iter)+(H+lambda.*diagH)\gradientVec;
            
            if sum(isnan(tempParams)) == 0
                tempLikelihood = GetLikelihood(reps,tempParams,peakNegativity,flashPoints);
                
                % if we do better
                tempCheck = tempLikelihood-logLikelihood(iter);
                if (tempCheck) > 0
                    parameterVec(:,iter+1) = tempParams;
                    lambda = lambda/10;
                    check = tempCheck;
                    logLikelihood(iter+1) = tempLikelihood;
                    % if we do worse
                else
                    check = 1;
                    lambda = lambda*10;
                    parameterVec(:,iter+1) = parameterVec(:,iter);
                    logLikelihood(iter+1) = logLikelihood(iter);
                end
                iter = iter+1;
            else
                break;
            end
        end
        bigParameterVec(repeats,:) = parameterVec(:,iter)';
        bigLogLikelihoods(repeats) = logLikelihood(iter);
    end
    [~,index] = max(bigLogLikelihoods);
    finalParameters(zz,:) = bigParameterVec(index,:);
    [fisherInfo(zz,:,:),ninetyfiveErrors(zz,:)] = getFisherInfo(finalParameters(zz,:),numParameters,h,reps,peakNegativity,flashPoints);
%     display(finalParameters(zz,:));
end

end

function [loglikelihood] = GetLikelihood(reps,parameterVec,peakNegativity,flashPoints)
loglikelihood = 0;
for kk=1:reps
    distX = flashPoints(kk,1)-parameterVec(2);
    distY = flashPoints(kk,2)-parameterVec(3);
    b = [parameterVec(1),parameterVec(4),parameterVec(5),parameterVec(7)];
    mu = b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-(distY.^2)./(2*b(3)*b(3)))+b(4);
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


