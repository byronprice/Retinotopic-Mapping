function [] = MapRetWrapper(fileStart)
%MapRetWrapper.m 
%  Wrapper for MapRetinotopy.m ... will use that script across multiple
%  animals and multiple days
%INPUT: Will run through all of the files in a given folder with the name
%         according to fileStart, e.g. fileStart = 'RetinoData*.plx';
%         to get all of the files only for a single animal, write
%         fileStart = 'RetinoData*26881*.plx';
%OUTPUT: stimVals
%        x -
%        y -
%        centerMass - 
%
% Created: 2016/07/19, 24 Cummington, Boston
%  Byron Price
% Updated: 2016/07/19
%  By: Byron Price 

fileList = dir(fileStart);
numFiles = size(fileList,1);

fileStart = fileList(1).name(1:end-18);
AnimalNames = zeros(numFiles,1);
for ii=1:numFiles
    AnimalNames(ii) = str2double(fileList(ii).name(end-8:end-4));
end
AnimalNames = unique(AnimalNames);
numAnimals = length(AnimalNames);

data = struct('Dates',{cell(numAnimals,1)},'Images',{cell(numAnimals,1)},...
    'Stats',{cell(numAnimals,1)},'numChans',zeros(numAnimals,1));
for ii=1:numAnimals
    list = dir(strcat(fileStart,'*',num2str(AnimalNames(ii)),'*.plx'));
    numDates = length(list);
    data.Dates{ii} = zeros(numDates,1);
    data.Images{ii} = cell(numDates,1);
    data.Stats{ii} = cell(numDates,1);
    for jj=1:numDates
        Date = str2double(list(jj).name(end-17:end-10));
        data.Dates{ii}(jj,1) = Date;
        [stimVals,centerMass,numChans] = MapRetinotopy(AnimalNames(ii),Date,0);
        data.Images{ii}{jj} = stimVals;
        data.Stats{ii}{jj} = centerMass;
        data.numChans(ii) = numChans;
    end    
end

for ii=1:numAnimals
    numDates = length(data.Dates{ii});numChans = data.numChans(ii);
    figure();
    stats = zeros(numDates,numChans,4);
    count = 1;
    for jj=1:numDates
        for kk=1:numChans
            subplot(numDates,numChans,count);imagesc(squeeze(data.Images{ii}{jj}(kk,:,:))');set(gca,'YDir','normal');h=colorbar;
            title(sprintf('Retinotopic Heat Map: Day %d, Channel %d',data.Dates{ii}(jj,1),kk));
            ylabel(h,'VEP Magnitude (\muV)');xlabel('Horizontal Screen Position (pixels)');
            ylabel('Vertical Screen Position (pixels)');colormap('jet');
            count = count+1;
        end
        stats(jj,:,:) = data.Stats{ii}{jj};
    end
    figure();count = 1;
    for kk=1:numChans
        for ll=1:4
            subplot(numChans,4,count);histogram(squeeze(stats(:,kk,ll)));
            if ll == 1
                title(sprintf('Histogram Center of Mass X-axis, Channel %d',kk));
                xlabel('Horizontal Position (pixels)');ylabel('Count');
            elseif ll == 2
                title(sprintf('Histogram Center of Mass Y-axis, Channel %d',kk));
                xlabel('Vertical Position (pixels)');ylabel('Count');
            elseif ll == 3
                title(sprintf('Histogram Spatial Spread X-axis, Channel %d',kk));
                xlabel('Standard Deviation (pixels)');ylabel('Count');
            elseif ll == 4
                title(sprintf('Histogram Spatial Spread Y-axis, Channel %d',kk));
                xlabel('Standard Deviation (pixels)');ylabel('Count');
            end
            count = count+1;
        end
    end
end
end

