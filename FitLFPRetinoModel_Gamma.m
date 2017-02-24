function [finalParameters,fisherInfo,ninetyfiveErrors] = FitLFPRetinoModel_Gamma(Response,xaxis,yaxis,centerVals)
%FitLFPRetinoModel_Gamma.m
%   Use data from LFP retinotopic mapping experiment to fit a non-linear
%    model of that retinotopy (data is maximum LFP magnitude in window 
%    from 150 to 250msec minus minimum magnitude in window from 50 to
%    120 msec after stimulus presentation, assumes a
%    Gamma likelihood)
%

%Created: 2017/02/22, 24 Cummington Mall, Boston
% Byron Price
%Updated: 2017/02/22
% By: Byron Price

%  model has 7 parameters, defined by vector p
%  data ~ N(mu,sigma), 
%    where mu = (p(1)*exp(-(xpos-p(2)).^2./(2*p(4)*p(4))-(ypos-p(3)).^2./(2*p(5)*p(5)))+p(7));
%    and sigma = p(6)

% parameter estimates are constrained to a reasonable range of values
Bounds = [0,600;min(xaxis),max(xaxis);min(yaxis),max(yaxis);1,1000;1,1000;1,2000;0,1000];
numChans = size(Response,1);

numParameters = 7;
   
finalParameters = zeros(numChans,numParameters);
fisherInfo = zeros(numChans,numParameters,numParameters);
ninetyfiveErrors = zeros(numChans,numParameters);
numRepeats = 100;
maxITER = 500;
tolerance = 1e-3;


for zz=1:numChans
%     display(sprintf('Running Data for Channel %d...',zz));

    Data = Response{zz};
    reps = size(Data,1);
    
    flashPoints = centerVals(squeeze(Data(:,1)),:);
    peakNegativity = abs(Data(:,2));
    h = ones(numParameters,1)./100;
    bigParameterVec = zeros(numParameters,numRepeats);
    bigLikelihood = zeros(numRepeats,1);
    % repeat gradient ascent from a number of different starting
    %  positions

    parfor repeats = 1:numRepeats
        parameterVec = zeros(numParameters,maxITER);
        logLikelihood = zeros(reps,maxITER);

        proposal = [1.75,100;4,250;4,250;1.75,150;1.75,150;1.5,5;1.75,200];
        for ii=1:numParameters 
              parameterVec(ii,1) = gamrnd(proposal(ii,1),proposal(ii,2));
        end
        
        [logLikelihood(:,1)] = GetLikelihood(reps,parameterVec(:,1),peakNegativity,flashPoints);
        
        check = 1;
        iter = 1;
        lambda = 100;
        % for each starting position, do maxITER iterations
        while abs(check) > tolerance && iter < maxITER
            [Jacobian,tempLikely] = GetJacobian(reps,parameterVec(:,iter),peakNegativity,flashPoints,numParameters,h,logLikelihood(:,iter));
            H = Jacobian'*Jacobian;
            update = pinv(H+lambda.*diag(diag(H)))*Jacobian'*((tempLikely-logLikelihood(:,iter))./h(1));
            tempParams = parameterVec(:,iter)+update;
            
            tempParams = max(Bounds(:,1),min(tempParams,Bounds(:,2)));
            
            [tempLikely] = GetLikelihood(reps,tempParams,peakNegativity,flashPoints);
            check = sum(tempLikely)-sum(logLikelihood(:,iter));
            if check <= 0
                parameterVec(:,iter+1) = parameterVec(:,iter);
                logLikelihood(:,iter+1) = logLikelihood(:,iter);
                lambda = min(lambda*10,1e15);
                check = 1;
            else
                parameterVec(:,iter+1) = tempParams;
                logLikelihood(:,iter+1) = tempLikely;
                lambda = max(lambda/10,1e-15);
            end
            iter = iter+1;
        end
        maxLikelies = sum(logLikelihood(:,1:iter),1);
        [bigLikelihood(repeats),index] = max(maxLikelies);
        bigParameterVec(:,repeats) = parameterVec(:,index);
    end
    [~,index] = max(bigLikelihood);
    finalParameters(zz,:) = bigParameterVec(:,index)';
    
    [fisherInfo(zz,:,:),ninetyfiveErrors(zz,:)] = getFisherInfo(finalParameters(zz,:),numParameters,h,reps,peakNegativity,flashPoints);
    display(zz);
    display(finalParameters(zz,:));
    display(ninetyfiveErrors(zz,:));
end
end

function [Jacobian,tempLikely] = GetJacobian(reps,parameterVec,peakNegativity,flashPoints,numParameters,h,prevLikely)
Jacobian = zeros(reps,numParameters);
tempLikely = zeros(reps,1);
for kk=1:reps
    for jj=1:numParameters
       tempParams = parameterVec;tempParams(jj) = tempParams(jj)+h(jj);
       mu = tempParams(1)*exp(-((flashPoints(kk,1)-tempParams(2)).^2)./(2*tempParams(4).^2)-...
        ((flashPoints(kk,2)-tempParams(3)).^2)./(2*tempParams(5).^2))+tempParams(7);
       likelihood = (peakNegativity(kk)*(-1/mu)-log(mu))/(1/tempParams(6))+tempParams(6)*log(tempParams(6))+...
           (tempParams(6)-1)*log(peakNegativity(kk))-log(gamma(tempParams(6)));
       Jacobian(kk,jj) = (likelihood-prevLikely(kk))/h(jj);
    end
    tempParams = parameterVec+h;
    mu = tempParams(1)*exp(-((flashPoints(kk,1)-tempParams(2)).^2)./(2*tempParams(4).^2)-...
        ((flashPoints(kk,2)-tempParams(3)).^2)./(2*tempParams(5).^2))+tempParams(7);
    tempLikely(kk) = (peakNegativity(kk)*(-1/mu)-log(mu))/(1/tempParams(6))+tempParams(6)*log(tempParams(6))+...
           (tempParams(6)-1)*log(peakNegativity(kk))-log(gamma(tempParams(6)));
end
end


function [loglikelihood] = GetLikelihood(reps,parameterVec,peakNegativity,flashPoints)
loglikelihood = zeros(reps,1);
for kk=1:reps
    mu = parameterVec(1)*exp(-((flashPoints(kk,1)-parameterVec(2)).^2)./(2*parameterVec(4).^2)-...
        ((flashPoints(kk,2)-parameterVec(3)).^2)./(2*parameterVec(5).^2))+parameterVec(7);
    loglikelihood(kk) = (peakNegativity(kk)*(-1/mu)-log(mu))/(1/parameterVec(6))+parameterVec(6)*log(parameterVec(6))+...
           (parameterVec(6)-1)*log(peakNegativity(kk))-log(gamma(parameterVec(6)));
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
            
            fisherInfo(jj,kk) = -(sum(likelyplusplus)-sum(likelyplusminus)-sum(likelyminusplus)+sum(likelyminusminus))./(4*deltaX*deltaY);
        else
            likely = GetLikelihood(reps,parameters,peakNegativity,flashPoints);
            
            parameterVec = parameters;parameterVec(firstParam) = parameterVec(firstParam)-deltaX;
            likelyminus = GetLikelihood(reps,parameterVec,peakNegativity,flashPoints);
            
            parameterVec = parameters;parameterVec(firstParam) = parameterVec(firstParam)+deltaX;
            likelyplus = GetLikelihood(reps,parameterVec,peakNegativity,flashPoints);
            
            fisherInfo(jj,kk) = -(sum(likelyminus)-2*sum(likely)+sum(likelyplus))./(deltaX*deltaX);
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


