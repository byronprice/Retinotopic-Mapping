function [finalParameters] = BayesianFitLFPModel(Response,xaxis,yaxis,centerVals)
% BayesianFitLFPModel.m

numChans = size(Response,1);
reps = size(Response,2);

numParameters = 7;
   
finalParameters = zeros(numChans,numParameters);

priorParams = [1.75,100;4,300;4,250;1.75,200;1.75,200;1.75,150;1.75,150];
% proposal = [150,100;1000,500;700,500;300,200;250,300;200,150;200,150];
        
N = 1000000;
x = zeros(N,numParameters);

for jj=1:numParameters
   x(1,jj) = gamrnd(priorParams(jj,1),priorParams(jj,2));
end


for zz=1:numChans
    randVals = rand([N,1]);
    flashPoints = centerVals(squeeze(Response(zz,:,1)),:);
    peakNegativity = squeeze(Response(zz,:,2));
    for ii=2:N
%         priorXstar = 1;priorXprev = 1;
        xStar = x(ii-1,:)+normrnd(0,1,[1,numParameters]);
%         check = 0;
%         for jj=1:numParameters
%             if xStar(jj) < 0
%                 check = 1;
%             end
% %             priorXstar = priorXstar*gampdf(xStar(jj),priorParams(jj,1),priorParams(jj,2));
% %             priorXprev = priorXprev*gampdf(x(ii-1,jj),priorParams(jj,1),priorParams(jj,2));
%         end

        likelihoodXstar = GetLikelihood(reps,xStar,peakNegativity,flashPoints);
        likelihoodXprev = GetLikelihood(reps,x(ii-1,:),peakNegativity,flashPoints);
%         A = min(1,(likelihoodXstar*priorXstar)./(likelihoodXprev*priorXprev));
        A = min(1,likelihoodXstar/likelihoodXprev);
        if randVals(ii) > A || sum(xStar<0) > 0
            x(ii,:) = x(ii-1,:);
        else
            x(ii,:) = xStar;
        end
    end
    figure();
    for ii=1:numParameters
       subplot(4,2,ii);histogram(x(100000:end,ii)); 
    end
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