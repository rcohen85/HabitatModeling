tic
GSvid = VideoWriter('GSPosition.avi','Motion JPEG AVI');
GSvid.FrameRate = 5
open(GSvid);
figDir = 'E:\hycom\0.08deg\fsle_movie';
fileList = dir(fullfile(figDir,'*.jpeg'));

for i=1:length(fileList)
    
    I = imread(fullfile(figDir,fileList(i).name));
    writeVideo(GSvid,I)
end

close(GSvid);
toc