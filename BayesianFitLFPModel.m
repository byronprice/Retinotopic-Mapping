function [finalParameters] = BayesianFitLFPModel(Response,xaxis,yaxis,centerVals)
% BayesianFitLFPModel.m

numChans = size(Response,1);
reps = size(Response,2);

numParameters = 7;
   
finalParameters = zeros(numChans,numParameters);

% priorParams = [1.75,100;1000,500;700,400;1.75,300;1.75,300;1.75,200;1.75,150];
proposal = [-150,100;1000,500;700,500;250,200;250,200;200,150;-200,150];

Bounds = [-500,0;min(xaxis),max(xaxis);min(yaxis),max(yaxis);1,1000;1,1000;1,1000;-500,0];
        
N = 50000;
x = zeros(N,numParameters);

for ii=1:numParameters
   x(1,ii) = Bounds(ii,1)+(Bounds(ii,2)-Bounds(ii,1)).*rand;
end

% for jj=[2,3]
%    x(1,jj) = normrnd(priorParams(jj,1),priorParams(jj,2)); 
% end
% for jj=[4,5,6]
%    x(1,jj) = gamrnd(priorParams(jj,1),priorParams(jj,2));
% end
% for jj=[1,7]
%     x(1,jj) = -gamrnd(priorParams(jj,1),priorParams(jj,2));
% end

for zz=1:numChans
    randVals = rand([N,1]);
    flashPoints = centerVals(squeeze(Response(zz,:,1)),:);
    peakNegativity = squeeze(Response(zz,:,2));
    for ii=2:N
%         priorXstar = 1;priorXprev = 1;qXstar = 1;qXprev = 1;
%         for jj=[2,3]
%             x(ii,jj) = normrnd(x(ii-1,jj),priorParams(jj,2));
%             priorXstar = priorXstar*normpdf(x(ii,jj),priorParams(jj,1),priorParams(jj,2));
%             priorXprev = priorXprev*normpdf(x(ii-1,jj),priorParams(jj,1),priorParams(jj,2));
%             qXstar = qXstar*normpdf(x(ii-1,jj),x(ii,jj),priorParams(jj,2));
%             qXprev = qXprev*normpdf(x(ii,jj),x(ii-1,jj),priorParams(jj,2));
%         end
%         for jj=[4,5,6]
%             x(ii,jj) = gamrnd(priorParams(jj,1),x(ii-1,jj)/priorParams(jj,1));
%             priorXstar = priorXstar*gampdf(x(ii,jj),priorParams(jj,1),priorParams(jj,2));
%             priorXprev = priorXprev*gampdf(x(ii-1,jj),priorParams(jj,1),priorParams(jj,2));
%             qXstar = qXstar*gampdf(x(ii-1,jj),priorParams(jj,1),x(ii,jj)/priorParams(jj,1));
%             qXprev = qXprev*gampdf(x(ii,jj),priorParams(jj,1),x(ii-1,jj)/priorParams(jj,1));
%         end
%         for jj=[1,7]
%             x(ii,jj) = -gamrnd(priorParams(jj,1),-x(ii-1,jj)/priorParams(jj,1));
%             priorXstar = priorXstar*gampdf(-x(ii,jj),priorParams(jj,1),priorParams(jj,2));
%             priorXprev = priorXprev*gampdf(-x(ii-1,jj),priorParams(jj,1),priorParams(jj,2));
%             qXstar = qXstar*gampdf(-x(ii-1,jj),priorParams(jj,1),-x(ii,jj)/priorParams(jj,1));
%             qXprev = qXprev*gampdf(-x(ii,jj),priorParams(jj,1),-x(ii-1,jj)/priorParams(jj,1));
%         end
        priorXstar = 1;priorXprev = 1;qXstar = 1;qXprev = 1;
        for jj=1:numParameters
            x(ii,jj) = normrnd(x(ii-1,jj),proposal(jj,2));
            if x(ii,jj) < Bounds(jj,1) || x(ii,jj) > Bounds(jj,2)
                priorXstar = 0;
            end
            if x(ii-1,jj) < Bounds(jj,1) || x(ii-1,jj) > Bounds(jj,2)
               priorXprev = 0; 
            end
            qXstar = qXstar*normpdf(x(ii-1,jj),x(ii,jj),proposal(jj,2));
            qXprev = qXprev*normpdf(x(ii,jj),x(ii-1,jj),proposal(jj,2));
        end
        
        likelihoodXstar = GetLikelihood(reps,x(ii,:),peakNegativity,flashPoints);
        likelihoodXprev = GetLikelihood(reps,x(ii-1,:),peakNegativity,flashPoints);
        
        A = min(1,(likelihoodXstar*priorXstar*qXstar)./(likelihoodXprev*priorXprev*qXprev));
        
        if randVals(ii) > A
            x(ii,:) = x(ii-1,:);
        end
    end
    figure();
    for ii=1:numParameters
       subplot(4,2,ii);histogram(x(5000:end,ii)); 
       mode(x(5000:end,ii))
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