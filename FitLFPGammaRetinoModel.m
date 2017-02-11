function [finalParameters,fisherInfo,ninetyfiveErrors] = FitLFPGammaRetinoModel(Response,xaxis,yaxis,centerVals)
%FitLFPGaussRetinoModel.m
%   Use data from LFP retinotopic mapping experiment to fit a non-linear
%    model of that retinotopy (using absolute value of the peak negativity 
%    in window from 50 to 120 msec after stimulus presentation and Gamma likelihood)

%Created: 2017/02/11, 10 Lawrence Street, Cambridge
% Byron Price
%Updated: 2017/02/11
% By: Byron Price

numChans = size(Response,1);
reps = size(Response,2);

numParameters = 7;
   
finalParameters = zeros(numChans,numParameters);
fisherInfo = zeros(numChans,numParameters,numParameters);
ninetyfiveErrors = zeros(numChans,numParameters);
numRepeats = 1000;
maxITER = 200;
tolerance = 1e-3;
h = ones(numParameters,1)./100;

% parameters are 
%  1) b(1) - rise of map at center
%  2)+3) the position (x and y) of the center of mass
%  4) b(2) - standard deviation or spread of the map in x
%  5) b(3) - standard deviation or spread of map in y
%  6) nu parameter from gamma likelihood
%  7) b(4) - peak negativity at edges of retinotopic region
%    b(1)+b(4) = peak negativity at retinotopic center of mass
%  
Bounds = [0,1000;min(xaxis),max(xaxis);min(yaxis),max(yaxis);1,2000;1,2000;0.1,2000;0,1000];

% display('Steepest Ascent ...');
for zz=1:numChans
%     display(sprintf('Running Data for Channel %d...',zz));
    flashPoints = centerVals(squeeze(Response(zz,:,1)),:);
    peakNegativity = abs(squeeze(Response(zz,:,2)));
    
    bigParameterVec = zeros(numRepeats,numParameters);
    bigLogLikelihoods = zeros(numRepeats,1);
    % repeat gradient ascent from a number of different starting
    %  positions
    parfor repeats = 1:numRepeats
        gradientVec = zeros(numParameters,1);
        parameterVec = zeros(maxITER,numParameters);
        logLikelihood = zeros(maxITER,1);
        %parameterVec(1,:) = squeeze(bigParameterVec(repeats,:));
        for ii=1:numParameters
            parameterVec(1,ii) = Bounds(ii,1)+(Bounds(ii,2)-Bounds(ii,1)).*rand;
        end
        
        % for each starting position, do maxITER iterations
        for iter=2:maxITER
            % calculate likelihood at the current position in parameter
            %  space
            [logLikelihood(iter-1)] = GetLikelihood(reps,parameterVec(iter-1,:),peakNegativity,flashPoints);
            
            if iter > 3
%                 check = sum(diff(logLikelihood(iter-3:iter-1)));
                check = diff(logLikelihood(iter-2:iter-1));
                if check < tolerance
                    break;
                end
            end            
            
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
            alpha = [0,1e-8,1e-6,1e-4,1e-2,1e0,1e1];
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
            
            for jj=1:numParameters
                parameterVec(iter,jj) = max(Bounds(jj,1),min(parameterVec(iter-1,jj)+alpha(ind)*gradientVec(jj),Bounds(jj,2)));
            end
        end
        bigParameterVec(repeats,:) = parameterVec(iter-1,:);
        bigLogLikelihoods(repeats) = logLikelihood(iter-1);
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
    nu = parameterVec(6);
    loglikelihood = loglikelihood-log(gamma(nu))+nu*...
       log(nu*peakNegativity(kk)/mu)-nu*peakNegativity(kk)/mu-...
       log(peakNegativity(kk));
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


