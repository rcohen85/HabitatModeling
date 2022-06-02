% Script to organize eddy trajectory data from Aviso
% Eddy trajectory atlas handbook: https://www.aviso.altimetry.fr/fileadmin/documents/data/tools/hdbk_eddytrajectory_META3.2_DT.pdf

% there's a lot of eddy data, trim to study period & region to make it more
% manageable
stDt = datenum(2016,5,1,0,0,0); % start date
edDt = datenum(2019,4,30,0,0,0); % end date
latMin = 25;
latMax = 45;
lonMax = -65;
lonMin = -82;
% directory of .nc files from Aviso
inDir = 'J:\Chpt_3\Eddies'; 

fileList = dir(fullfile(inDir,'*.nc'));

eddyTable = table;
for i=1:length(fileList) % for each eddy trajectory file
    
    % load file
    myFile = fullfile(inDir,fileList(i).name); 
    
    % extract variables of interest
    centerLat = ncread(myFile,'latitude');
    centerLon=ncread(myFile,'longitude');
    amp=ncread(myFile,'amplitude');
    area=ncread(myFile,'speed_area');
    speed=ncread(myFile,'speed_average');
    date=ncread(myFile,'time');
    
    % create variable coding for cyclonic/anticyclonic
    if ~isempty(strfind(fileList(i).name,'Anti'))
        polarity = repmat({'A'},length(date),1);
    else
        polarity = repmat({'C'},length(date),1);
    end
    
    % create variable coding for long-lived vs short-lived eddies
    if ~isempty(strfind(fileList(i).name,'long'))
        dur = repmat({'long'},length(date),1);
    elseif ~isempty(strfind(fileList(i).name,'short'))
        dur = repmat({'short'},length(date),1);
    end
    
    % convert timestamp to Matlab datenum
    time_origin = datenum(1950,1,1,0,0,0);
    date = datenum(date) + time_origin;
    
    eddyTable = vertcat(eddyTable,table(date,centerLat,centerLon,polarity,area,amp,speed,dur));
    
end

% Find timestamps within study period
goodInd = find(eddyTable.date>=stDt & eddyTable.date<=edDt);
eddyTable = eddyTable(goodInd,:);

% Find eddies within study region
goodLat = find(eddyTable.centerLat>=latMin & eddyTable.centerLat<=latMax);
goodLon = find(eddyTable.centerLon<=lonMax+360 & eddyTable.centerLon>=lonMin+360);
goodCoords = intersect(goodLat,goodLon);
eddyTable = eddyTable(goodCoords,:);

% Sort table chronologically
eddyTable = sortrows(eddyTable,'date');

% Convert timestamps to strings
eddyTable.date = datestr(eddyTable.date);

% Save table as .csv 
writetable(eddyTable,fullfile(inDir,'/AllEddies.csv'));

