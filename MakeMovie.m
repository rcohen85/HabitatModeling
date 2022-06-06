tic
GSvid = VideoWriter('AllEddies_Temp.avi','Motion JPEG AVI');
GSvid.FrameRate = 5
open(GSvid);
figDir = 'E:\ModelingCovarData\EddyMoviePlots_Temp';
fileList = dir(fullfile(figDir,'*.jpeg'));

for i=1:length(fileList)
    
    I = imread(fullfile(figDir,fileList(i).name));
    writeVideo(GSvid,I)
end

close(GSvid);
toc