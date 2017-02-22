% RetinotopyAnalysis.m
%  Use data from 2016 summer retinotopic mapping experiments with new
%   retinotopic region models (Gaussian distribution and binomial
%   distribution) ... try to figure which most accurately portrays the 
%   data
%  fileStart = 'RetinoMap*.mat';
function [] = RetinotopyAnalysis(fileStart)
xaxis = 1:2560; yaxis = 1:1440;

files = dir(fileStart);
numFiles = length(files);

AIC = zeros(3,1);numParameters = [6,7,8];
for zz=1:numFiles
   load(files(zz).name);
   
   try
       numChans = MapParams.numChans;
       centerVals = MapParams.centerVals;
       tempResponse = MapParams.Response;
       numStimuli = size(tempResponse,2);
       numReps = size(tempResponse,3);
       temp = MapParams.centerMass;
       
       centerMass = zeros(numChans,2);
       centerMass(:,1) = temp.x;
       centerMass(:,2) = temp.y;
       
       display(strcat('Running file-',files(zz).name));
%        filename = files(zz).name;filename = strcat(filename(1:end-4),'.fig');
%        openfig(filename);
       for ii=1:numChans
           if isnan(centerMass(ii,1)) == 0
%                dataFile = [dataFile;files(zz).name];
               gaussResponse = zeros(1,numStimuli*numReps,2);
               binomResponse = zeros(1,numStimuli*numReps,2);
               
               xPos = unique(centerVals(:,1));yPos = unique(centerVals(:,2));
               VEPs = zeros(length(xPos),length(yPos),numReps);
               count = 1;
               for jj=1:numStimuli
                   for kk=1:numReps
                       gaussResponse(1,count,1) = jj;
                       gaussResponse(1,count,2) = max(tempResponse(ii,jj,kk,150:250))-min(tempResponse(ii,jj,kk,50:120));
%                        gaussResponse(1,count,2) = -min(tempResponse(ii,jj,kk,50:120));
                       
                       x = centerVals(jj,1);y = centerVals(jj,2);
                       xInd = find(xPos == x);yInd = find(yPos == y);
                       VEPs(xInd,yInd,kk) = gaussResponse(1,count,2);
                       binomResponse(1,count,1) = jj;
                       if gaussResponse(1,count,2) < -150
                           binomResponse(1,count,2) = 1;
                       end
                       count = count+1;
                   end
               end
               
%                figure();count = 1;edges = 0:25:600;
%                for kk=length(yPos):-1:1
%                    for jj=1:length(xPos)
%                       temp = squeeze(VEPs(jj,kk,:));
%                       subplot(length(yPos),length(xPos),count);histogram(temp,edges);
%                       stdev = std(temp);meanVal = mean(temp);
%                       title(sprintf('%3.0f  %3.0f',meanVal,stdev));
%                       count = count+1;
%                    end
%                end
                
               [littleAIC] = FitModel1(gaussResponse,xaxis,yaxis,centerVals);
               AIC(1) = AIC(1)+littleAIC;
%                
               [littleAIC] = FitModel2(gaussResponse,xaxis,yaxis,centerVals);
               AIC(2) = AIC(2)+littleAIC;

               [littleAIC] = FitModel3(gaussResponse,xaxis,yaxis,centerVals);
               AIC(3) = AIC(3)+littleAIC;
%                [bernoulliParameters] = FitLFPretinoModel(binomResponse,xaxis,yaxis,0.2,centerVals);
%                allBinomParams = [allBinomParams;bernoulliParameters];
%                tempX = bernoulliParameters(1,2);
%                tempY = bernoulliParameters(1,3);
%                
%                errorX = abs(tempX-centerMass(ii,1))/centerMass(ii,1);
%                errorY = abs(tempY-centerMass(ii,2))/centerMass(ii,2);
%                binomError = (errorX+errorY)/2;display(binomError);
%                allBinomErrors = [allBinomErrors,binomError];
           end
       end
   catch
       display(strcat('Skipping File-'),files(zz).name);
       continue;
   end
end
figure();plot(numParameters,AIC,'b','LineWidth',2);title('Retinotopic-Map Model Comparison: AIC');
xlabel('Number of Parameters');ylabel('AIC');
save('LFPRetinotopy_ModelComparison.mat','AIC','numParameters');
end

% finalIm = zeros(length(xaxis),length(yaxis));
% for ii=1:length(xaxis)
%     for jj=1:length(yaxis)
%         distance = DistFun(finalParameters(1,2:3),[ii,jj]);
%         finalIm(ii,jj) = hyperParameterFun([finalParameters(1,1),finalParameters(1,4),-120],distance);
%     end
% end
% figure();imagesc(finalIm');set(gca,'YDir','normal');

function [AIC] = FitModel1(Response,xaxis,yaxis,centerVals)
%FitLFPRetinoModel_LM.m

% parameter estimates are constrained to a reasonable range of values
Bounds = [0,1000;min(xaxis),max(xaxis);min(yaxis),max(yaxis);1,1000;1,1000;0,1000];
numChans = size(Response,1);
reps = size(Response,2);

numParameters = 6;
   
finalParameters = zeros(numChans,numParameters);
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

        proposal = [1.75,100;4,250;4,250;1.75,150;1.5,150;1.75,200];
        for ii=1:numParameters 
              parameterVec(ii,1) = gamrnd(proposal(ii,1),proposal(ii,2));
        end
        
        [yhat] = Getyhat1(reps,parameterVec(:,1),flashPoints);
        squaredResiduals(1) = sum((peakNegativity-yhat).^2);
        
        check = 1;
        iter = 1;
        lambda = 100;
        % for each starting position, do maxITER iterations
        while abs(check) > tolerance && iter < maxITER
            [Jacobian] = GetJacobian1(reps,parameterVec(:,iter),flashPoints,numParameters,h,yhat);
            H = Jacobian'*Jacobian;
            update = pinv(H+lambda.*diag(diag(H)))*Jacobian'*(peakNegativity-yhat);
            tempParams = parameterVec(:,iter)+update;
            
            tempParams = max(Bounds(:,1),min(tempParams,Bounds(:,2)));
            
            [tempYhat] = Getyhat1(reps,tempParams,flashPoints);
            squaredResiduals(iter+1) = sum((peakNegativity-tempYhat).^2);
            check = diff(squaredResiduals(iter:iter+1));
            if check >= 0
                parameterVec(:,iter+1) = parameterVec(:,iter);
                lambda = min(lambda*10,1e10);
                check = 1;
                squaredResiduals(iter+1) = squaredResiduals(iter);
            else
                parameterVec(:,iter+1) = tempParams;
                yhat = tempYhat;
                lambda = max(lambda/10,1e-10);
            end
            iter = iter+1;
        end
        [bigResiduals(repeats),index] = min(squaredResiduals(1:iter));
        bigParameterVec(:,repeats) = parameterVec(:,index);
    end
    [~,index] = min(bigResiduals);
    finalParameters(zz,:) = bigParameterVec(:,index)';
    [yhat] = Getyhat1(reps,finalParameters(zz,:),flashPoints);
    finalParameters(zz,5) = std(peakNegativity-yhat);

    parameterVec = zeros(maxITER,numParameters);gradientVec = zeros(1,numParameters);
    logLikelihood = zeros(maxITER,1);
    
    parameterVec(1,:) = finalParameters(zz,:);
    [logLikelihood(1)] = GetLikelihood1(reps,parameterVec(1,:),peakNegativity,flashPoints);
    check = 1;iter = 1;lambda = 100;
    while abs(check) > tolerance && iter < maxITER
        for jj=1:numParameters
            tempParameterVec = parameterVec(iter,:);
            tempParameterVec(jj) = tempParameterVec(jj)+h(jj);
            [gradLikelihoodplus] = GetLikelihood1(reps,tempParameterVec,peakNegativity,flashPoints);
            
            gradientVec(jj) = (gradLikelihoodplus-logLikelihood(iter))./h(jj);
        end
        tempParams = parameterVec(iter,:)+lambda*gradientVec;
        tempLikelihood = GetLikelihood1(reps,tempParams,peakNegativity,flashPoints);
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
    [maxLikely] = max(logLikelihood(1:iter));
    AIC = 2*numParameters-2*maxLikely;
end
end

function [Jacobian] = GetJacobian1(reps,parameterVec,flashPoints,numParameters,h,yhat)
Jacobian = zeros(reps,numParameters);
for kk=1:reps
    for jj=1:numParameters
       tempParams = parameterVec;tempParams(jj) = tempParams(jj)+h(jj);
       mu = tempParams(1)*exp(-((flashPoints(kk,1)-tempParams(2)).^2)./(2*tempParams(4).^2)-...
        ((flashPoints(kk,2)-tempParams(3)).^2)./(2*tempParams(4).^2))+tempParams(6);
       Jacobian(kk,jj) = (mu-yhat(kk))/h(jj);
    end
end
end

function [yhat] = Getyhat1(reps,parameterVec,flashPoints)
yhat = zeros(reps,1);
for kk=1:reps
    yhat(kk) = parameterVec(1)*exp(-((flashPoints(kk,1)-parameterVec(2)).^2)./(2*parameterVec(4).^2)-...
        ((flashPoints(kk,2)-parameterVec(3)).^2)./(2*parameterVec(4).^2))+parameterVec(6);

end
end

function [loglikelihood] = GetLikelihood1(reps,parameterVec,peakNegativity,flashPoints)
loglikelihood = 0;
for kk=1:reps
    mu = parameterVec(1)*exp(-((flashPoints(kk,1)-parameterVec(2)).^2)./(2*parameterVec(4).^2)-...
        ((flashPoints(kk,2)-parameterVec(3)).^2)./(2*parameterVec(4).^2))+parameterVec(6);
    stdev = parameterVec(5);%*exp(-(distX.^2)./(2*b(2)*b(2))-b(5)*distX*distY/(2*b(2)*b(3))-(distY.^2)./(2*b(3)*b(3)))+parameterVec(9);
    loglikelihood = loglikelihood-(1/2)*log(2*pi*stdev*stdev)-(1/(2*stdev*stdev))*(peakNegativity(kk)-mu).^2;
end

end

function [AIC] = FitModel2(Response,xaxis,yaxis,centerVals)
%FitLFPRetinoModel_LM.m

% parameter estimates are constrained to a reasonable range of values
Bounds = [0,1000;min(xaxis),max(xaxis);min(yaxis),max(yaxis);1,1000;1,1000;1,1000;0,1000];
numChans = size(Response,1);
reps = size(Response,2);

numParameters = 7;
   
finalParameters = zeros(numChans,numParameters);
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
        
        [yhat] = Getyhat2(reps,parameterVec(:,1),flashPoints);
        squaredResiduals(1) = sum((peakNegativity-yhat).^2);
        
        check = 1;
        iter = 1;
        lambda = 100;
        % for each starting position, do maxITER iterations
        while abs(check) > tolerance && iter < maxITER
            [Jacobian] = GetJacobian2(reps,parameterVec(:,iter),flashPoints,numParameters,h,yhat);
            H = Jacobian'*Jacobian;
            update = pinv(H+lambda.*diag(diag(H)))*Jacobian'*(peakNegativity-yhat);
            tempParams = parameterVec(:,iter)+update;
            
            tempParams = max(Bounds(:,1),min(tempParams,Bounds(:,2)));
            
            [tempYhat] = Getyhat2(reps,tempParams,flashPoints);
            squaredResiduals(iter+1) = sum((peakNegativity-tempYhat).^2);
            check = diff(squaredResiduals(iter:iter+1));
            if check >= 0
                parameterVec(:,iter+1) = parameterVec(:,iter);
                lambda = min(lambda*10,1e10);
                check = 1;
                squaredResiduals(iter+1) = squaredResiduals(iter);
            else
                parameterVec(:,iter+1) = tempParams;
                yhat = tempYhat;
                lambda = max(lambda/10,1e-10);
            end
            iter = iter+1;
        end
        [bigResiduals(repeats),index] = min(squaredResiduals(1:iter));
        bigParameterVec(:,repeats) = parameterVec(:,index);
    end
    [~,index] = min(bigResiduals);
    finalParameters(zz,:) = bigParameterVec(:,index)';
    [yhat] = Getyhat2(reps,finalParameters(zz,:),flashPoints);
    finalParameters(zz,6) = std(peakNegativity-yhat);

    parameterVec = zeros(maxITER,numParameters);gradientVec = zeros(1,numParameters);
    logLikelihood = zeros(maxITER,1);
    
    parameterVec(1,:) = finalParameters(zz,:);
    [logLikelihood(1)] = GetLikelihood2(reps,parameterVec(1,:),peakNegativity,flashPoints);
    check = 1;iter = 1;lambda = 100;
    while abs(check) > tolerance && iter < maxITER
        for jj=1:numParameters
            tempParameterVec = parameterVec(iter,:);
            tempParameterVec(jj) = tempParameterVec(jj)+h(jj);
            [gradLikelihoodplus] = GetLikelihood2(reps,tempParameterVec,peakNegativity,flashPoints);
            
            gradientVec(jj) = (gradLikelihoodplus-logLikelihood(iter))./h(jj);
        end
        tempParams = parameterVec(iter,:)+lambda*gradientVec;
        tempLikelihood = GetLikelihood2(reps,tempParams,peakNegativity,flashPoints);
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
    [maxLikely] = max(logLikelihood(1:iter));
    AIC = 2*numParameters-2*maxLikely;
end
end

function [Jacobian] = GetJacobian2(reps,parameterVec,flashPoints,numParameters,h,yhat)
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

function [yhat] = Getyhat2(reps,parameterVec,flashPoints)
yhat = zeros(reps,1);
for kk=1:reps
    yhat(kk) = parameterVec(1)*exp(-((flashPoints(kk,1)-parameterVec(2)).^2)./(2*parameterVec(4).^2)-...
        ((flashPoints(kk,2)-parameterVec(3)).^2)./(2*parameterVec(5).^2))+parameterVec(7);

end
end

function [loglikelihood] = GetLikelihood2(reps,parameterVec,peakNegativity,flashPoints)
loglikelihood = 0;
for kk=1:reps
    mu = parameterVec(1)*exp(-((flashPoints(kk,1)-parameterVec(2)).^2)./(2*parameterVec(4).^2)-...
        ((flashPoints(kk,2)-parameterVec(3)).^2)./(2*parameterVec(5).^2))+parameterVec(7);
    stdev = parameterVec(6);%*exp(-(distX.^2)./(2*b(2)*b(2))-b(5)*distX*distY/(2*b(2)*b(3))-(distY.^2)./(2*b(3)*b(3)))+parameterVec(9);
    loglikelihood = loglikelihood-(1/2)*log(2*pi*stdev*stdev)-(1/(2*stdev*stdev))*(peakNegativity(kk)-mu).^2;
end

end

function [AIC] = FitModel3(Response,xaxis,yaxis,centerVals)
%FitLFPRetinoModel_LM.m

% parameter estimates are constrained to a reasonable range of values
Bounds = [0,1000;min(xaxis),max(xaxis);min(yaxis),max(yaxis);1,1000;1,1000;1,1000;0,1000;-1,1];
numChans = size(Response,1);
reps = size(Response,2);

numParameters = 8;
   
finalParameters = zeros(numChans,numParameters);
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

        proposal = [1.75,100;4,250;4,250;1.75,150;1.75,150;1.5,150;1.75,200;0.5,0.4];
        for ii=1:numParameters 
              parameterVec(ii,1) = gamrnd(proposal(ii,1),proposal(ii,2));
        end
        
        [yhat] = Getyhat3(reps,parameterVec(:,1),flashPoints);
        squaredResiduals(1) = sum((peakNegativity-yhat).^2);
        
        check = 1;
        iter = 1;
        lambda = 100;
        % for each starting position, do maxITER iterations
        while abs(check) > tolerance && iter < maxITER
            [Jacobian] = GetJacobian3(reps,parameterVec(:,iter),flashPoints,numParameters,h,yhat);
            H = Jacobian'*Jacobian;
            update = pinv(H+lambda.*diag(diag(H)))*Jacobian'*(peakNegativity-yhat);
            tempParams = parameterVec(:,iter)+update;
            
            tempParams = max(Bounds(:,1),min(tempParams,Bounds(:,2)));
            
            [tempYhat] = Getyhat3(reps,tempParams,flashPoints);
            squaredResiduals(iter+1) = sum((peakNegativity-tempYhat).^2);
            check = diff(squaredResiduals(iter:iter+1));
            if check >= 0
                parameterVec(:,iter+1) = parameterVec(:,iter);
                lambda = min(lambda*10,1e10);
                check = 1;
                squaredResiduals(iter+1) = squaredResiduals(iter);
            else
                parameterVec(:,iter+1) = tempParams;
                yhat = tempYhat;
                lambda = max(lambda/10,1e-10);
            end
            iter = iter+1;
        end
        [bigResiduals(repeats),index] = min(squaredResiduals(1:iter));
        bigParameterVec(:,repeats) = parameterVec(:,index);
    end
    [~,index] = min(bigResiduals);
    finalParameters(zz,:) = bigParameterVec(:,index)';
    [yhat] = Getyhat3(reps,finalParameters(zz,:),flashPoints);
    finalParameters(zz,6) = std(peakNegativity-yhat);

    parameterVec = zeros(maxITER,numParameters);gradientVec = zeros(1,numParameters);
    logLikelihood = zeros(maxITER,1);
    
    parameterVec(1,:) = finalParameters(zz,:);
    [logLikelihood(1)] = GetLikelihood3(reps,parameterVec(1,:),peakNegativity,flashPoints);
    check = 1;iter = 1;lambda = 100;
    while abs(check) > tolerance && iter < maxITER
        for jj=1:numParameters
            tempParameterVec = parameterVec(iter,:);
            tempParameterVec(jj) = tempParameterVec(jj)+h(jj);
            [gradLikelihoodplus] = GetLikelihood3(reps,tempParameterVec,peakNegativity,flashPoints);
            
            gradientVec(jj) = (gradLikelihoodplus-logLikelihood(iter))./h(jj);
        end
        tempParams = parameterVec(iter,:)+lambda*gradientVec;
        tempLikelihood = GetLikelihood3(reps,tempParams,peakNegativity,flashPoints);
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
    [maxLikely] = max(logLikelihood(1:iter));
    AIC = 2*numParameters-2*maxLikely;
end
end

function [Jacobian] = GetJacobian3(reps,parameterVec,flashPoints,numParameters,h,yhat)
Jacobian = zeros(reps,numParameters);
for kk=1:reps
    for jj=1:numParameters
       tempParams = parameterVec;tempParams(jj) = tempParams(jj)+h(jj);
       mu = tempParams(1)*exp(-((flashPoints(kk,1)-tempParams(2)).^2)./(2*tempParams(4).^2)-...
        ((flashPoints(kk,2)-tempParams(3)).^2)./(2*tempParams(5).^2)-...
        tempParams(8).*(flashPoints(kk,1)-tempParams(2))*(flashPoints(kk,2)-tempParams(3))./(2*tempParams(4)*tempParams(5)))+tempParams(7);
       Jacobian(kk,jj) = (mu-yhat(kk))/h(jj);
    end
end
end

function [yhat] = Getyhat3(reps,parameterVec,flashPoints)
yhat = zeros(reps,1);
for kk=1:reps
    yhat(kk) = parameterVec(1)*exp(-((flashPoints(kk,1)-parameterVec(2)).^2)./(2*parameterVec(4).^2)-...
        ((flashPoints(kk,2)-parameterVec(3)).^2)./(2*parameterVec(5).^2)-...
        parameterVec(8)*(flashPoints(kk,1)-parameterVec(2))*(flashPoints(kk,2)-parameterVec(3))./(2*parameterVec(4)*parameterVec(5)))+parameterVec(7);

end
end

function [loglikelihood] = GetLikelihood3(reps,parameterVec,peakNegativity,flashPoints)
loglikelihood = 0;
for kk=1:reps
    mu = parameterVec(1)*exp(-((flashPoints(kk,1)-parameterVec(2)).^2)./(2*parameterVec(4).^2)-...
        ((flashPoints(kk,2)-parameterVec(3)).^2)./(2*parameterVec(5).^2)-...
        parameterVec(8)*(flashPoints(kk,1)-parameterVec(2))*(flashPoints(kk,2)-parameterVec(3))./(2*parameterVec(4)*parameterVec(5)))+parameterVec(7);
    stdev = parameterVec(6);%*exp(-(distX.^2)./(2*b(2)*b(2))-b(5)*distX*distY/(2*b(2)*b(3))-(distY.^2)./(2*b(3)*b(3)))+parameterVec(9);
    loglikelihood = loglikelihood-(1/2)*log(2*pi*stdev*stdev)-(1/(2*stdev*stdev))*(peakNegativity(kk)-mu).^2;
end

end