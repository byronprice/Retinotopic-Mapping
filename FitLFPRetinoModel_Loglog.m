function [finalParameters,fisherInfo,ninetyfiveErrors,result,Deviance,chi2p] = FitLFPRetinoModel_Loglog(Response,xaxis,yaxis,numRepeats)
%FitLFPRetinoModel_Loglog.m
%   Use data from LFP retinotopic mapping experiment to fit a non-linear
%    model of that retinotopy (data is maximum LFP magnitude in window 
%    from 180 to 270msec minus minimum magnitude in window from 100 to
%    160 msec after stimulus presentation, assumes a
%    log-logistic likelihood)
%

%Created: 2017/03/17, 24 Cummington Mall, Boston
% Byron Price
%Updated: 2017/04/10
% By: Byron Price

%  model has 7 parameters, defined by vector p
%  data ~ log-logistic(log mean, log scale), 
%    where the log mean of the data is as follows
%     log mean = p(1)*exp(-(xpos-p(2)).^2./(2*p(4)*p(4))-
%        (ypos-p(3)).^2./(2*p(5)*p(5)))+p(6);
%
%  and p(7) = log scale


% parameter estimates are constrained to a reasonable range of values
Bounds = [1e-2,2;min(xaxis)-50,max(xaxis)+50;min(yaxis)-50,max(yaxis)+50;50,1000;50,1000;1e-3,10;1e-3,0.9];
numChans = size(Response,1);

numParameters = 7;
   
finalParameters = zeros(numChans,numParameters);
fisherInfo = zeros(numChans,numParameters,numParameters);
ninetyfiveErrors = zeros(numChans,numParameters);
result = zeros(numChans,1);
chi2p = struct('Saturated_Full',zeros(numChans,1),'Full_Null',zeros(numChans,1));
Deviance = struct('Full',cell(numChans,1),'Null',cell(numChans,1));

% numRepeats = 1e4;
maxITER = 1000;
likelyTolerance = 1e-6;
gradientTolerance = 1e-6;

for zz=1:numChans
%     display(sprintf('Running Data for Channel %d...',zz));

    Data = Response{zz};
    reps = size(Data,1);
    
    flashPoints = Data(:,1:2);
    vepMagnitude = Data(:,3);
    
    h = ones(numParameters,1)./1000;
    bigParameterVec = zeros(numParameters,numRepeats);
    bigLikelihood = zeros(numRepeats,1);
    % repeat gradient ascent from a number of different starting
    %  positions
    
    phat = mle(vepMagnitude,'distribution','loglogistic');
    parfor repeats = 1:numRepeats
        parameterVec = zeros(numParameters,maxITER);
        totalLikelihood = zeros(maxITER,1);
        
        proposal = [0.5,0.5;5,200;2.6,152;6.7,38.6;5.9,40.3;8.6,10;1.5,0.1];
        for ii=1:numParameters-2
              parameterVec(ii,1) = gamrnd(proposal(ii,1),proposal(ii,2));
        end
        parameterVec(6,1) = phat(1)+normrnd(0,1);
        parameterVec(7,1) = phat(2)+normrnd(0,proposal(7,2));

        parameterVec(:,1) = max(Bounds(1:numParameters,1),min(parameterVec(:,1),Bounds(1:numParameters,2)));
        
        logLikelihood = GetLikelihood(reps,parameterVec(:,1),vepMagnitude,flashPoints);
        
        totalLikelihood(1) = sum(logLikelihood);

        check = 1;
        iter = 1;
        lambda = 1;
        update = ones(numParameters,1);
%         figure();scatter(iter,totalLikelihood(1));pause(1);hold on;
            % for each starting position, do maxITER iterations
        try
            while abs(check) > likelyTolerance && iter < maxITER && sum(abs(update)) > gradientTolerance
                [Jacobian] = GetJacobian(reps,parameterVec(:,iter),flashPoints,numParameters,h,logLikelihood,vepMagnitude);
                likelihoodForGradient = GetLikelihood(reps,parameterVec(:,iter)+h,vepMagnitude,flashPoints);
                H = Jacobian'*Jacobian;
                update = pinv(H+lambda.*diag(diag(H)))*Jacobian'*(likelihoodForGradient-logLikelihood);
                tempParams = parameterVec(:,iter)+update;
            
                tempParams = max(Bounds(1:numParameters,1),min(tempParams,Bounds(1:numParameters,2)));
            
                tempLikelihood = GetLikelihood(reps,tempParams,vepMagnitude,flashPoints);
                
                totalLikelihood(iter+1) = sum(tempLikelihood);
                check = diff(totalLikelihood(iter:iter+1));
                
                % lambda updating from Marco Giordan, Federico Vaggi, Ron
                %    Wehrens arXiv  ... maximum likelihood exponential
                %    family
                gain = (sum(tempLikelihood)-sum(logLikelihood))/...
                    (-0.5*(tempParams-parameterVec(:,iter))'*H*(tempParams-parameterVec(:,iter)));
                lambda = lambda*max(1/3,1-(2*gain-1).^2)*(gain>0)+lambda*2*(gain<=0);
                if check <= 0
                    parameterVec(:,iter+1) = parameterVec(:,iter);
                    check = 1;
                    totalLikelihood(iter+1) = totalLikelihood(iter);
                else
                    parameterVec(:,iter+1) = tempParams;
                    logLikelihood = tempLikelihood;
                end
                iter = iter+1;%scatter(iter,totalLikelihood(iter));pause(0.01);
            end
            
            [bigLikelihood(repeats),index] = max(totalLikelihood(1:iter));
            bigParameterVec(:,repeats) = parameterVec(:,index);
        catch
            
        end
    end
    logicalInds = bigLikelihood~=0;%figure();histogram(bigDeviance(logicalInds),1:50:1000);
    [maxLikely,index] = max(bigLikelihood(logicalInds));
    
    tempBigParams = bigParameterVec(:,logicalInds);
    finalParameters(zz,1:numParameters) = tempBigParams(:,index)';

    
    [fisherInfo(zz,:,:),ninetyfiveErrors(zz,:)] = getFisherInfo(finalParameters(zz,:),numParameters,h,reps,vepMagnitude,flashPoints);
    
    [deviance] = GetDeviance(reps,finalParameters(zz,:),vepMagnitude,flashPoints);
    %Deviance.Full(zz) = sum(deviance)./exp(finalParameters(zz,end));
    % need to divide by constant 1.228
    Deviance(zz).Full = deviance./1.228;
    chi2p.Saturated_Full(zz) = 1-chi2cdf(sum(Deviance(zz).Full),reps-numParameters);
    
    
    [nullDeviance] = GetNullDeviance(reps,vepMagnitude,phat);
    Deviance(zz).Null = nullDeviance./1.228;
    chi2p.Full_Null(zz) = 1-chi2cdf(sum(Deviance(zz).Null)-sum(Deviance(zz).Full),numParameters-length(phat));
    
end
end

function [Jacobian] = GetJacobian(reps,parameterVec,flashPoints,numParameters,h,oldLikely,vepMagnitude)
Jacobian = zeros(reps,numParameters);
for kk=1:reps
    for jj=1:numParameters
       tempParams = parameterVec;tempParams(jj) = tempParams(jj)+h(jj);
       tempMu = tempParams(1)*exp(-((flashPoints(kk,1)-tempParams(2)).^2)./(2*tempParams(4).^2)-...
        ((flashPoints(kk,2)-tempParams(3)).^2)./(2*tempParams(5).^2))+tempParams(6);
       tempLikely = -log(tempParams(7))-log(vepMagnitude(kk))+(log(vepMagnitude(kk))-tempMu)/tempParams(7)-...
        2*log(1+exp((log(vepMagnitude(kk))-tempMu)/tempParams(7)));
       Jacobian(kk,jj) = (tempLikely-oldLikely(kk))/h(jj);
    end
end
end

function [loglikelihood] = GetLikelihood(reps,parameterVec,vepMagnitude,flashPoints)
loglikelihood = zeros(reps,1);
for kk=1:reps
    mu = parameterVec(1)*exp(-((flashPoints(kk,1)-parameterVec(2)).^2)./(2*parameterVec(4).^2)-...
        ((flashPoints(kk,2)-parameterVec(3)).^2)./(2*parameterVec(5).^2))+parameterVec(6);
    loglikelihood(kk) = -log(parameterVec(7))-log(vepMagnitude(kk))+(log(vepMagnitude(kk))-mu)/parameterVec(7)-...
        2*log(1+exp((log(vepMagnitude(kk))-mu)/parameterVec(7)));
end
end

function [deviance] = GetDeviance(reps,parameterVec,vepMagnitude,flashPoints)
deviance = zeros(reps,1);

for kk=1:reps
    mu = parameterVec(1)*exp(-((flashPoints(kk,1)-parameterVec(2)).^2)./(2*parameterVec(4).^2)-...
        ((flashPoints(kk,2)-parameterVec(3)).^2)./(2*parameterVec(5).^2))+parameterVec(6);
    deviance(kk) =  2*(-2*log(2)-(log(vepMagnitude(kk))-mu)/parameterVec(7)+...
        2*log(1+exp((log(vepMagnitude(kk))-mu)/parameterVec(7))));
end
end

function [nullDeviance] = GetNullDeviance(reps,vepMagnitude,phat)
nullDeviance = zeros(reps,1);

for kk=1:reps
    nullDeviance(kk) =  2*(-2*log(2)-(log(vepMagnitude(kk))-phat(1))/phat(2)+...
        2*log(1+exp((log(vepMagnitude(kk))-phat(1))/phat(2)))); 
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

inverseFisherInfo = pinv(fisherInfo);
for ii=1:numParameters
    errors(ii) = sqrt(inverseFisherInfo(ii,ii));
end

if isreal(errors) == 0
    temp = sqrt(errors.*conj(errors));
    errors = 1.96.*temp;
%     fprintf('Warning: Complex Errors in Model Fit\n\n');
elseif isreal(errors) == 1
    errors = 1.96.*errors;
end
end