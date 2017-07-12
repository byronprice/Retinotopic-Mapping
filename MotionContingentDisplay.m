% MotionContingentDisplay.m

s = PL_InitClient(0);
if s == 0
   return
end

% setup server, send signal to stim computer that this computer is ready to
%  begin
display('Opening TCP/IP Server ...');
tcpipServer = tcpip('128.197.59.166',30000,'NetworkRole','server');
bufferSize = 50000;
set(tcpipServer,'OutputBufferSize',bufferSize);
fopen(tcpipServer);

p = PL_GetPars(s);
sampleFreq = p(14);

startEXP = 254;
endEXP = 255;

% approximate sigma = k*MAD 
k = 1.253;% 1.4826; for median absolute deviation

chanInd = 49;

% collect some data to get a baseline of the noise
%  wait for the signal from the stim computer that the experiment is ready
%  to begin
totalHeld = 30*sampleFreq;
D = zeros(30*sampleFreq,1);
display('Obtaining estimate of noise ...');

index = 1;
check = 0;

while check == 0
    pause(0.01);
    [~,tEvs] = PL_GetTS(s);
    [n,~,d] = PL_GetADV(s);

    D(index:index+n-1) = d(:,chanInd);
    index = index+n;
    check = sum(tEvs(:,3) == startEXP);
    
    if index >= totalHeld-n
        index = 1;
    end
end

threshold = (1*k).*mad(abs(D));

if threshold < 1e-7 || threshold > 1e-3
    threshold = 1e-6;
end
fprintf('Threshold: %3.2e\n',threshold);

% this is the actual experiment
%  perform simple data analysis on the VEPs to determine whether or not a
%  significant VEP has occurred, relay the information to the stim computer



display('Beginning experiment ...');
check = 0;
while check == 0
    
    [~, tEvs] = PL_GetTS(s);
    % tEvs contains the event timeStamps ... if tEvs(x,1) = 4 and tEvs(x,2)
    % = 257 , then tEvs(x,3) = strobed event word and tEvs(x,4) = timeStamp
    % this can be compared against t from GetADV and t(x,4) for all of the
    % events in which tEvs(x,2) = 6 (your desired channel)
    [~, ~, d] = PL_GetADV(s);
    temp1 = abs(d(:,chanInd));
    temp = sum(temp1>=threshold)/length(temp1);
    fwrite(tcpipServer,temp,'double');
    
    check = sum(tEvs(:,3) == endEXP);
    pause(1);
end

display('Experiment over.');
% you need to call PL_Close(s) to close the connection
% with the Plexon server
PL_Close(s);
s = 0;
fclose(tcpipServer);
