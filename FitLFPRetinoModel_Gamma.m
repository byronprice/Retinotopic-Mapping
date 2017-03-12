function [finalParameters,fisherInfo,ninetyfiveErrors] = FitLFPRetinoModel_Gamma(Response,xaxis,yaxis)
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
%  data ~ gamma(k,theta), 
%    where the mode of the data follows
%     mode = (p(1)*exp(-(xpos-p(2)).^2./(2*p(4)*p(4))-(ypos-p(3)).^2./(2*p(5)*p(5)))+p(6));
%    and mode = (k-1)*theta

% parameter estimates are constrained to a reasonable range of values
Bounds = [0,600;min(xaxis)-50,max(xaxis)+50;min(yaxis)-50,max(yaxis)+50;1,1000;1,1000;0,1000;1,2000];
numChans = size(Response,1);

numParameters = 7;
   
finalParameters = zeros(numChans,numParameters);
fisherInfo = zeros(numChans,numParameters,numParameters);
ninetyfiveErrors = zeros(numChans,numParameters);
numRepeats = 100;
maxITER = 100;
tolerance = 1e-3;


for zz=1:numChans
%     display(sprintf('Running Data for Channel %d...',zz));

    Data = Response{zz};
    reps = size(Data,1);
    
    flashPoints = Data(:,1:2);
    vepMagnitude = abs(Data(:,3));
    
    h = ones(numParameters,1)./100;
    bigParameterVec = zeros(numParameters,numRepeats);
    bigLikelihood = zeros(numRepeats,1);
    % repeat gradient ascent from a number of different starting
    %  positions

    parfor repeats = 1:numRepeats
        parameterVec = zeros(numParameters,maxITER);
        logLikelihood = zeros(reps,maxITER);

        proposal = [1.75,100;4,250;4,250;1.75,150;1.75,150;1.5,200;1.75,100];
        for ii=1:numParameters 
              parameterVec(ii,1) = gamrnd(proposal(ii,1),proposal(ii,2));
        end
        
        [logLikelihood(:,1)] = GetLikelihood(reps,parameterVec(:,1),vepMagnitude,flashPoints);

        check = 1;
        iter = 1;
        lambda = 100;
        % for each starting position, do maxITER iterations
        while abs(check) > tolerance && iter < maxITER
            [Jacobian,tempLikely] = GetJacobian(reps,parameterVec(:,iter),vepMagnitude,flashPoints,numParameters,h,logLikelihood(:,iter));
            H = Jacobian'*Jacobian;
            update = pinv(H+lambda.*diag(diag(H)))*Jacobian'*((tempLikely-logLikelihood(:,iter))); % or /h(1)
            tempParams = parameterVec(:,iter)+update;
            
            tempParams = max(Bounds(:,1),min(tempParams,Bounds(:,2)));
            
            [tempLikely] = GetLikelihood(reps,tempParams,vepMagnitude,flashPoints);
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
    
    [fisherInfo(zz,:,:),ninetyfiveErrors(zz,:)] = getFisherInfo(finalParameters(zz,:),numParameters,h,reps,vepMagnitude,flashPoints);
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
       mode = tempParams(1)*exp(-((flashPoints(kk,1)-tempParams(2)).^2)./(2*tempParams(4).^2)-...
        ((flashPoints(kk,2)-tempParams(3)).^2)./(2*tempParams(5).^2))+tempParams(6);
       variance = tempParams(7)^2;
       theta = max(1e-10,(2*variance)/(mode+sqrt(mode*mode+4*variance)));
       k = max(1e-10,1+mode/theta);
       likelihood = -log(gamma(k))-k*log(theta)+(k-1)*log(peakNegativity(kk))-peakNegativity(kk)/theta;
       Jacobian(kk,jj) = (likelihood-prevLikely(kk))/h(jj);
    end
    tempParams = parameterVec+h;
    variance = tempParams(7)^2;
    mode = tempParams(1)*exp(-((flashPoints(kk,1)-tempParams(2)).^2)./(2*tempParams(4).^2)-...
        ((flashPoints(kk,2)-tempParams(3)).^2)./(2*tempParams(5).^2))+tempParams(6);
    theta = (2*variance)/(mode+sqrt(mode*mode+4*variance));
    k = 1+mode/theta;
    tempLikely(kk) = -log(gamma(k))-k*log(theta)+(k-1)*log(peakNegativity(kk))-peakNegativity(kk)/theta;
end
end


function [loglikelihood] = GetLikelihood(reps,parameterVec,peakNegativity,flashPoints)
loglikelihood = zeros(reps,1);
variance = parameterVec(7)^2;
for kk=1:reps
    mode = parameterVec(1)*exp(-((flashPoints(kk,1)-parameterVec(2)).^2)./(2*parameterVec(4).^2)-...
        ((flashPoints(kk,2)-parameterVec(3)).^2)./(2*parameterVec(5).^2))+parameterVec(6);
    theta = (2*variance)/(mode+sqrt(mode*mode+4*variance));
    k = 1+mode/theta;
    loglikelihood(kk) = -log(gamma(k))-k*log(theta)+(k-1)*log(peakNegativity(kk))-peakNegativity(kk)/theta;
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


