clear

[filename,path] = uigetfile('C:\Users\Janis\Documents\Munoz Lab\Data\indy\*.xlsx');

task = 'prior04';

[~,~,ext] = fileparts(filename);
if strcmp(ext,'.xlsx')
    [~,~,raw] = xlsread([path,filename],'Sheet1');
    relCell = cellfun(@findstr,repmat({task},[size(raw,1),1]),raw(:,2),'UniformOutput',false);
    relIndex = ~cellfun(@isempty,relCell);
    relevant = raw(relIndex,:);
    filename = relevant(:,1);
    filename = strcat(filename,'-r.mat');
else
    filename = {filename};
end


for n = 1:size(filename,1)
    findTrials([path,filename{n}])
end