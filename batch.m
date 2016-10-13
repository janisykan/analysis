% 1) goes into each file and plot pupil size, rasters, and spk density
% function
% 2) plot average pupil size and spk density for all files.
% 3) compares spk density of each condition with no pop-out condition using
% signed rank test without holm-bonferroni correction
% 4) THIS TASK ONLY FOR excel sheet with info on visual/motor activity. Use runles_working
% to get this info.

clear

[filename,path] = uigetfile('C:\Users\Janis\Documents\Munoz Lab\Data\uly\','multiselect','on');

task = 'multipopout4v2';
xmin = -200;
xmax = 800; % time window to display spikes

if ~iscell(filename)
    [~,~,ext] = fileparts(filename);
    if strcmp(ext,'.xlsx')
        [~,~,raw] = xlsread([path,filename],1);
        relCell = cellfun(@findstr,repmat({task},[size(raw,1),1]),raw(:,2),'UniformOutput',false);
        relIndex = ~cellfun(@isempty,relCell);
%         relIndex = find(strcmp(raw(:,2),task));
        relevant = raw(relIndex,:);
        filename = relevant(:,1);
        filename = strcat(filename,'-r.mat');
        rew1 = mat2cell(cell2mat(relevant(:,8))+cell2mat(relevant(:,9))/2,ones(1,size(relevant,1)));
        rew2 = mat2cell(cell2mat(relevant(:,10))+cell2mat(relevant(:,11))/2,ones(1,size(relevant,1)));
        depth = relevant(:,5);
        unit = relevant(:,3);
        isolation = relevant(:,4);
        actualtask = relevant(:,2);
        layer = relevant(:,14);
    else
        filename = {filename};
        depth = {input('depth of electrode rel. SC surface: ')};
        rew1 = {input('reward 1: ')};
        rew2 = {input('reward 2: ')};
    end
end

file = cell2struct([filename,actualtask,layer,unit,isolation,rew1,rew2,depth],{'name','task','layer','unit','isolation','rew1','rew2','depth'},2);
multi4v3Index = [];

if strcmpi(task,'prior04')
    condNames = {'1 high rew IN','1 high rew OUT', '1 low rew IN', '1 low rew OUT', '2 high rew', '2 low rew', 'high IN low OUT', 'low IN high OUT'};
    colors = {'m-','m:','g-','g:','m--','g--','m-.','g-.'};
    [file.condNames] = deal(condNames);
    [file.colors] = deal(colors);
else
    condNames = {'0 targets','1 target IN','1 target OUT', '2 targets IN', '2 targets OUT', '4 targets'};
    colors = {'k','g','m','c','r','b'};
    [file.condNames] = deal(condNames);
    [file.colors] = deal(colors);
    multi4v3Index = find(cellfun(@strcmpi,repmat({'multipopout4v3'},[length(actualtask),1]),actualtask));
    [file(multi4v3Index).condNames] = deal([condNames, condNames(2:end)]);
    [file(multi4v3Index).colors] = deal({'k','g','m','c','r','b','g--','m--','c--','r--','b--'});
%     if strcmp(task,'multipopout')
%         multi4v3Index = find(cellfun(@strcmpi,repmat({'multipopout4v3'},[length(actualtask),1]),actualtask));
%         [file(multi4v3Index).condNames] = deal([condNames, condNames(2:end)]);
%         [file(multi4v3Index).colors] = deal({'k','g','m','c','r','b','g--','m--','c--','r--','b--'});
%     end
end

[file.path] = deal(path);

n = 0;
SCsN = 0;
SCiN = 0;
SUs = 0;
SUi = 0;
useless = 0;
if ~isempty(multi4v3Index)
    numCond = size(file(multi4v3Index(1)).condNames, 2);
else
    numCond = size(condNames,2);
end

for y = 1:numCond
    eval(['SCs' num2str(y-1) '=[];']);
    eval(['SCi' num2str(y-1) '=[];']);
    eval(['lfpS' num2str(y-1) '=[];']);
    eval(['lfpI' num2str(y-1) '=[];']);
end
for x = 1:size(filename,1) % for each file
    if strcmp(file(x).layer,'s')
        [convspks, lfp] = neural3(file(x));
        if ischar(convspks)
            useless = useless + 1;
            continue
        end
        SCsN = SCsN+1;
        if strcmpi(file(x).isolation,'SU')
            SUs = SUs+1;
        end
        maxSpk = ceil(max(max(convspks,[],2)));
        baseline = ceil(min(min(convspks,[],2)));
        convspks_norm = (convspks-baseline)/(maxSpk-baseline);
        maxSpk = ceil(max(max(abs(lfp),[],2)));
        baseline = ceil(min(min(abs(lfp),[],2)));
        lfp_norm = (lfp-baseline)/(maxSpk-baseline);
        for y = 1:size(convspks,1) % for each condition
            eval(['SCs' num2str(y-1) '= [SCs' num2str(y-1) '; convspks_norm(y,:)];']);
            eval(['lfpS' num2str(y-1) '= [lfpS' num2str(y-1) '; lfp_norm(y,:)];']);
        end
    elseif strcmp(file(x).layer,'i')
        [convspks, lfp] = neural3(file(x));
        if ischar(convspks)
            useless = useless + 1;
            continue
        end
        SCiN = SCiN + 1;
        if strcmpi(file(x).isolation,'SU')
            SUi = SUi+1;
        end
        maxSpk = ceil(max(max(convspks,[],2)));
        baseline = ceil(min(min(convspks,[],2)));
        convspks_norm = (convspks-baseline)/(maxSpk-baseline);
        maxSpk = ceil(max(max(abs(lfp),[],2)));
        baseline = ceil(min(min(abs(lfp),[],2)));
        lfp_norm = (lfp-baseline)/(maxSpk-baseline);
        for y = 1:size(convspks,1)
            eval(['SCi' num2str(y-1) '= [SCi' num2str(y-1) '; convspks_norm(y,:)];']);
            eval(['lfpI' num2str(y-1) '= [lfpI' num2str(y-1) '; lfp_norm(y,:)];']);
        end
    end
    [pupil] = behavior2(file(x));
    if ischar(pupil)
%         display(pupil)
%         file(x).name
        continue
    end
    for y = 1:size(pupil,1)
        eval(['c' num2str(y-1) '(x, :)=pupil(y,:);']);
    end
    n = n+1;
end

close all

%% plot neural data for all files
if SCsN ~=0 || SCiN ~=0
       
    
    figure % SCs spk density
    hold on
    sigma=5; nstd=sigma*5;
    tt=xmin-nstd:xmax+nstd-1;
    for x = 1:size(file(1).colors,2)
        eval(['avg_spksS(x,:)=nanmean(SCs' num2str(x-1) ',1);']);
        hSpk(x) = plot(tt,avg_spksS(x,:),file(1).colors{x});
    end
    
    ymax=ceil(max(max(avg_spksS,[],2)));
    ymin=floor(min(min(avg_spksS,[],2)));
    xlabel('time from target (ms) ')
    ylabel('spikes/s ')
    title(['all cells SCs n = ' num2str(SCsN)])
    legend(hSpk,file(1).condNames)
    set(gca,'XTick', xmin:100:xmax);
    axis([0 200 ymin ymax])
    
    figure % SCi spk density
    hold on
    sigma=5; nstd=sigma*5;
    tt=xmin-nstd:xmax+nstd-1;
    for x = 1:size(file(1).colors,2)
        eval(['avg_spksI(x,:)=nanmean(SCi' num2str(x-1) ',1);']);
        hSpk(x) = plot(tt,avg_spksI(x,:),file(1).colors{x});
    end
    
    ymax=ceil(max(max(avg_spksI,[],2)));
    ymin=floor(min(min(avg_spksI,[],2)));
    xlabel('time from target (ms) ')
    ylabel('spikes/s ')
    title(['all cells SCi n = ' num2str(SCiN)])
    legend(hSpk,file(1).condNames)
    set(gca,'XTick', xmin:100:xmax);
    axis([0 200 ymin ymax])
    
    figure % SCs LFP
    hold on
    tt=xmin:xmax-1;
    if ~isempty(multi4v3Index)
        col = file(multi4v3Index(1)).colors;
    else
        col = colors;
    end
    for x = 1:size(col,2)
        eval(['avg_lfpS(x,:)=nanmean(lfpS' num2str(x-1) ',1);']);
        hSpk(x) = plot(tt,avg_lfpS(x,:),col{x});
    end
    
    ymax=ceil(max(max(avg_lfpS,[],2)));
    ymin=floor(min(min(avg_lfpS,[],2)));
    xlabel('time from target (ms) ')
    ylabel('spikes/s ')
    title(['all cells lfp SCs n = ' num2str(SCsN)])
    legend(hSpk,file(1).condNames)
    set(gca,'XTick', xmin:100:xmax);
%     axis([0 200 ymin ymax])
    
    figure % SCi LFP
    hold on
    tt=xmin:xmax-1;
    for x = 1:size(col,2)
        eval(['avg_lfpI(x,:)=nanmean(lfpI' num2str(x-1) ',1);']);
        hSpk(x) = plot(tt,avg_lfpI(x,:),col{x});
    end
    
    ymax=ceil(max(max(avg_lfpS,[],2)));
    ymin=floor(min(min(avg_lfpS,[],2)));
    xlabel('time from target (ms) ')
    ylabel('spikes/s ')
    title(['all cells lfp SCi n = ' num2str(SCiN)])
    legend(hSpk,file(1).condNames)
    set(gca,'XTick', xmin:100:xmax);
%     axis([0 200 ymin ymax])

%% plotting lfp - spk density
    figure % single condition, all cells
    hold on
    sigma=5; nstd=sigma*5;
    clear hSpk
    for x = 1:size(SCs0,1)
       hSpk(x) =  plot(tt,SCs0(x,nstd+1:end-nstd)-lfpS0(x,:));
    end
    title('lfp - spk. Single condition, all cells')
    
    figure % single cell, all condition
    hold on
    clear hSpk
    for x = 1:size(file(1).colors,2)
        eval(['hSpk(x) = plot(tt,SCs' num2str(x-1) '(1,nstd+1:end-nstd)-lfpS0(x,:),file(1).colors{x});']);
    end
    legend(hSpk,file(1).condNames)
    title('lfp - spk. Single cell, all conditions')
    
    figure % avg SCs cell, all conditions
    hold on
    clear hSpk
    for x = 1:size(file(1).colors,2)
        hSpk(x) = plot(tt,avg_spksS(x,nstd+1:end-nstd)-avg_lfpS(x,:),file(1).colors{x});
    end
    legend(hSpk,file(1).condNames)
    title('lfp - spk. avg SCs cells, all conditions')
    
    figure % avg SCi cell, all conditions
    hold on
    clear hSpk
    for x = 1:size(file(1).colors,2)
        hSpk(x) = plot(tt,avg_spksI(x,nstd+1:end-nstd)-avg_lfpI(x,:),file(1).colors{x});
    end
    legend(hSpk,file(1).condNames) 
    title('lfp - spk. avg SCi cells, all conditions')
    
    keyboard

%% plotting spk density vs lfp

    figure % single condition, all cells
    hold on
    sigma=5; nstd=sigma*5;
    clear hSpk
    for x = 1:size(SCs0,1)
        hSpk(x) = plot(SCs0(x,nstd+1:end-nstd),lfpS0(x,:));
    end
    legend(hSpk,file(1).condNames)
    title('spk vs lfp. single condition, all cells')

    figure % single cell, all condition
    hold on
    clear hSpk
    for x = 1:size(file(1).colors,2)
        eval(['hSpk(x) = plot(SCs' num2str(x-1) '(1,nstd+1:end-nstd),lfpS0(x,:),file(1).colors{x});']);
    end
    legend(hSpk,file(1).condNames)
    title('spk vs lfp. single cell, all conditions')
    
    figure % avg SCs cell, all conditions
    hold on
    clear hSpk
    for x = 1:size(file(1).colors,2)
        hSpk(x) = plot(avg_spksS(x,nstd+1:end-nstd),avg_lfpS(x,:),file(1).colors{x});
    end
    legend(hSpk,file(1).condNames)
    title('spk vs lfp. avg SCs cells, all conditions')
    
    figure % avg SCi cell, all conditions
    hold on
    clear hSpk
    for x = 1:size(file(1).colors,2)
        hSpk(x) = plot(avg_spksI(x,nstd+1:end-nstd),avg_lfpI(x,:),file(1).colors{x});
    end
    legend(hSpk,file(1).condNames)
    title('spk vs lfp. avg SCi cells, all conditions')
    
    keyboard
    
    tt=xmin-nstd:xmax+nstd-1;
    %% Running wilcoxon signed rank test on neural data
    for y = 1:size(convspks,1)-1 % for each condition
        figure
        hold on
        eval(['curSpks = SCs' num2str(y) ';']);
        sigTimes = [];
        for x=3:size(curSpks,2)-2 % running average
            condA = nanmean(curSpks(:,x-2:x+2),2);
            condB = nanmean(SCs0(:,x-2:x+2),2);
            p = signrank(condA,condB);
            if p<=0.05
                plot([x+(xmin-nstd) x+(xmin-nstd)], [-0.04 -0.02],'k') % tick for significant difference
            end
        end
        
        % plot standard error
        errorUpper = nanmean(curSpks,1)+(std(curSpks,0,1)/sqrt(size(curSpks,1)-1));
        errorLower = nanmean(curSpks,1)-(std(curSpks,0,1)/sqrt(size(curSpks,1)-1));
        ymax = max(errorUpper);
        h = fill([tt,fliplr(tt)],[errorUpper fliplr(errorLower)],[1 0.6 0.6]);
        set(h,'EdgeColor','none')
        errorUpper = nanmean(SCs0,1)+(std(SCs0,0,1)/sqrt(size(SCs0,1)-1));
        errorLower = nanmean(SCs0,1)-(std(SCs0,0,1)/sqrt(size(SCs0,1)-1));
        h = fill([tt,fliplr(tt)],[errorUpper fliplr(errorLower)],[0.8 0.8 0.8]);
        set(h,'EdgeColor','none')
        
        xlabel('time from target (ms) ')
        ylabel('spikes/s ')
        title(['SCs signed rank test condition = ' file(1).condNames{y+1}])
        set(gca,'XTick', xmin:100:xmax);
        axis([0 200 -0.04 ymax])
        
        % plot the average curves
        eval(['avg_spks=nanmean(SCs' num2str(y) ',1);']);
        plot(tt,avg_spks(1,:),file(1).colors{y+1});
        avg_spks=nanmean(SCs0,1);
        plot(tt,avg_spks(1,:),'k');
        
        
    end
    
  
    for y = 1:size(convspks,1)-1 % for each condition
        figure
        hold on
        eval(['curSpks = SCi' num2str(y) ';']);
        for x=6:5:size(curSpks,2)-5 % running average
            condA = nanmean(curSpks(:,x-2:x+2),2);
            condB = nanmean(SCi0(:,x-2:x+2),2);
            p = signrank(condA,condB);
            if p<=0.05
                plot([x+(xmin-nstd) x+(xmin-nstd)], [-0.04 -0.02],'k')
            end
        end
        % plot standard error
        errorUpper = nanmean(curSpks,1)+(std(curSpks,0,1)/sqrt(size(curSpks,1)-1));
        errorLower = nanmean(curSpks,1)-(std(curSpks,0,1)/sqrt(size(curSpks,1)-1));
        h = fill([tt,fliplr(tt)],[errorUpper fliplr(errorLower)],[1 0.6 0.6]);
        set(h,'EdgeColor','none')
        errorUpper = nanmean(SCi0,1)+(std(SCi0,0,1)/sqrt(size(SCi0,1)-1));
        errorLower = nanmean(SCi0,1)-(std(SCi0,0,1)/sqrt(size(SCi0,1)-1));
        h = fill([tt,fliplr(tt)],[errorUpper fliplr(errorLower)],[0.8 0.8 0.8]);
        set(h,'EdgeColor','none')
        
        % plot the average curves
        eval(['avg_spks=nanmean(SCi' num2str(y) ',1);']);
        plot(tt,avg_spks(1,:),file(1).colors{y+1});
        avg_spks=nanmean(SCi0,1);
        plot(tt,avg_spks(1,:),'k');
        
        xlabel('time from target (ms) ')
        ylabel('spikes/s ')
        title(['SCi signed rank test condition = ' file(1).condNames{y+1}])
        set(gca,'XTick', xmin:100:xmax);
        axis([0 200 -0.04 ymax])
    end
    
    %% getting LFP data into format for repeated measures ANOVA
    % first, check significant window using all 3 conditions with target in
    % warning, lots of hard-coding in here. Will only work when used with
    % multipopout4v2 or multipopout
    figure
    hold on
    inSpks = cat(3,SCs1,SCs3,SCs5); % combine all conditions with target in
    inSpks = mean(inSpks,3);
    sigTimes = [];
    for x=3:size(curSpks,2)-2 % running average
        condA = nanmean(inSpks(:,x-2:x+2),2);
        condB = nanmean(SCs0(:,x-2:x+2),2);
        p = signrank(condA,condB);
        if p<=0.05
            plot([x+(xmin-nstd) x+(xmin-nstd)], [-0.04 -0.02],'k') % tick for significant difference
            if x+(xmin-nstd)>0 && x+(xmin-nstd)<200
                sigTimes = [sigTimes,x];
            end
        end
    end
    
    % plot standard error
    errorUpper = nanmean(inSpks,1)+(std(inSpks,0,1)/sqrt(size(inSpks,1)-1));
    errorLower = nanmean(inSpks,1)-(std(inSpks,0,1)/sqrt(size(inSpks,1)-1));
    h = fill([tt,fliplr(tt)],[errorUpper fliplr(errorLower)],[1 0.6 0.6]);
    set(h,'EdgeColor','none')
    errorUpper = nanmean(SCs0,1)+(std(SCs0,0,1)/sqrt(size(SCs0,1)-1));
    errorLower = nanmean(SCs0,1)-(std(SCs0,0,1)/sqrt(size(SCs0,1)-1));
    h = fill([tt,fliplr(tt)],[errorUpper fliplr(errorLower)],[0.8 0.8 0.8]);
    set(h,'EdgeColor','none')
    
    xlabel('time from target (ms) ')
    ylabel('spikes/s ')
    title('SCs signed rank test all conditions')
    set(gca,'XTick', xmin:100:xmax);
    axis([-100 300 -0.04 ymax])
    
    avg_spks=nanmean(inSpks,1);
    plot(tt,avg_spks(1,:),file(1).colors{2});
    avg_spks=nanmean(SCs0,1);
    plot(tt,avg_spks(1,:),'k');
    
    ANOVAwin = [sigTimes(1) sigTimes(end)];
    keyboard
    ANOVAs = [];
    for y = 1:size(convspks,1)
        eval(['ANOVAs' num2str(y-1) '= nanmean(SCs' num2str(y-1) '(:,ANOVAwin(1):ANOVAwin(2)),2);']);
        eval(['errorS' num2str(y-1) '= (std(ANOVAs' num2str(y-1) ',0,1)/sqrt(size(ANOVAs' num2str(y-1) ',1)-1));']);
    end
    
    figure
    hold on
    inSpki = cat(3,SCi1,SCi3,SCi5); % combine all conditions with target in
    inSpki = mean(inSpki,3);
    count = 0;
    for x=3:size(curSpks,2)-2 % running average
        condA = nanmean(inSpki(:,x-2:x+2),2);
        condB = nanmean(SCi0(:,x-2:x+2),2);
        p = signrank(condA,condB);
        if p<=0.05
            plot([x+(xmin-nstd) x+(xmin-nstd)], [-0.04 -0.02],'k') % tick for significant difference
            if x+(xmin-nstd)>0 && x+(xmin-nstd)<120
                sigTimes = [sigTimes,x];
            end
        end
    end
    
    % plot standard error
    errorUpper = nanmean(inSpki,1)+(std(inSpki,0,1)/sqrt(size(inSpki,1)-1));
    errorLower = nanmean(inSpki,1)-(std(inSpki,0,1)/sqrt(size(inSpki,1)-1));
    h = fill([tt,fliplr(tt)],[errorUpper fliplr(errorLower)],[1 0.6 0.6]);
    set(h,'EdgeColor','none')
    errorUpper = nanmean(SCi0,1)+(std(SCi0,0,1)/sqrt(size(SCi0,1)-1));
    errorLower = nanmean(SCi0,1)-(std(SCi0,0,1)/sqrt(size(SCi0,1)-1));
    h = fill([tt,fliplr(tt)],[errorUpper fliplr(errorLower)],[0.8 0.8 0.8]);
    set(h,'EdgeColor','none')
    
    xlabel('time from target (ms) ')
    ylabel('spikes/s ')
    title('SCi signed rank test conditions')
    set(gca,'XTick', xmin:100:xmax);
    axis([-100 300 -0.04 ymax])
    
    avg_spks=nanmean(inSpki,1);
    plot(tt,avg_spks(1,:),file(1).colors{2});
    avg_spks=nanmean(SCi0,1);
    plot(tt,avg_spks(1,:),'k');
    
    ANOVAwin = [sigTimes(1) sigTimes(end)];
    ANOVAi = [];
    for y = 1:size(convspks,1)
        eval(['ANOVAi' num2str(y-1) '= nanmean(SCi' num2str(y-1) '(:,ANOVAwin(1):ANOVAwin(2)),2);']);
        eval(['errorI' num2str(y-1) '= (std(ANOVAi' num2str(y-1) ',0,1)/sqrt(size(ANOVAi' num2str(y-1) ',1)-1));']);
    end
    
    ANOVAiIN = [mean(ANOVAi0),mean(ANOVAi1),mean(ANOVAi3),mean(ANOVAi5)];
    ANOVAiOUT = [mean(ANOVAi0),mean(ANOVAi2),mean(ANOVAi4)];
    ANOVAsIN = [mean(ANOVAs0),mean(ANOVAs1),mean(ANOVAs3),mean(ANOVAs5)];
    ANOVAsOUT = [mean(ANOVAs0),mean(ANOVAs2),mean(ANOVAs4)];
    error = [errorS0, errorS1, errorS3, errorS5];
    figure
    bar(ANOVAsIN)
    hold on
    errorbar(ANOVAsIN, error, 'rx')
    ylim([0 0.5])
    figure
    error = [errorS0, errorS2, errorS4];
    bar(ANOVAsOUT)
    hold on
    errorbar(ANOVAsOUT, error, 'rx')
    ylim([0 0.5])
    figure
    error = [errorI0, errorI1, errorI3, errorI5];
    bar(ANOVAiIN)
    hold on
    errorbar(ANOVAiIN, error, 'rx')
    ylim([0 0.5])
    figure
    error = [errorI0, errorI2, errorI4];
    bar(ANOVAiOUT)
    hold on
    errorbar(ANOVAiOUT, error, 'rx')
    ylim([0 0.5])
    keyboard
    
    %% getting data into format for repeated measures ANOVA
    % first, check significant window using all 3 conditions with target in
    % warning, lots of hard-coding in here. Will only work when used with
    % multipopout4v2 or multipopout
    figure
    hold on
    inSpks = cat(3,SCs1,SCs3,SCs5); % combine all conditions with target in
    inSpks = mean(inSpks,3);
    sigTimes = [];
    for x=3:size(curSpks,2)-2 % running average
        condA = nanmean(inSpks(:,x-2:x+2),2);
        condB = nanmean(SCs0(:,x-2:x+2),2);
        p = signrank(condA,condB);
        if p<=0.05
            plot([x+(xmin-nstd) x+(xmin-nstd)], [-0.04 -0.02],'k') % tick for significant difference
            if x+(xmin-nstd)>0 && x+(xmin-nstd)<200
                sigTimes = [sigTimes,x];
            end
        end
    end
    
    % plot standard error
    errorUpper = nanmean(inSpks,1)+(std(inSpks,0,1)/sqrt(size(inSpks,1)-1));
    errorLower = nanmean(inSpks,1)-(std(inSpks,0,1)/sqrt(size(inSpks,1)-1));
    h = fill([tt,fliplr(tt)],[errorUpper fliplr(errorLower)],[1 0.6 0.6]);
    set(h,'EdgeColor','none')
    errorUpper = nanmean(SCs0,1)+(std(SCs0,0,1)/sqrt(size(SCs0,1)-1));
    errorLower = nanmean(SCs0,1)-(std(SCs0,0,1)/sqrt(size(SCs0,1)-1));
    h = fill([tt,fliplr(tt)],[errorUpper fliplr(errorLower)],[0.8 0.8 0.8]);
    set(h,'EdgeColor','none')
    
    xlabel('time from target (ms) ')
    ylabel('spikes/s ')
    title('SCs signed rank test all conditions')
    set(gca,'XTick', xmin:100:xmax);
    axis([-100 300 -0.04 ymax])
    
    avg_spks=nanmean(inSpks,1);
    plot(tt,avg_spks(1,:),file(1).colors{2});
    avg_spks=nanmean(SCs0,1);
    plot(tt,avg_spks(1,:),'k');
    
    ANOVAwin = [sigTimes(1) sigTimes(end)];
    keyboard
    ANOVAs = [];
    for y = 1:size(convspks,1)
        eval(['ANOVAs' num2str(y-1) '= nanmean(SCs' num2str(y-1) '(:,ANOVAwin(1):ANOVAwin(2)),2);']);
        eval(['errorS' num2str(y-1) '= (std(ANOVAs' num2str(y-1) ',0,1)/sqrt(size(ANOVAs' num2str(y-1) ',1)-1));']);
    end
    
    figure
    hold on
    inSpki = cat(3,SCi1,SCi3,SCi5); % combine all conditions with target in
    inSpki = mean(inSpki,3);
    count = 0;
    for x=3:size(curSpks,2)-2 % running average
        condA = nanmean(inSpki(:,x-2:x+2),2);
        condB = nanmean(SCi0(:,x-2:x+2),2);
        p = signrank(condA,condB);
        if p<=0.05
            plot([x+(xmin-nstd) x+(xmin-nstd)], [-0.04 -0.02],'k') % tick for significant difference
            if x+(xmin-nstd)>0 && x+(xmin-nstd)<120
                sigTimes = [sigTimes,x];
            end
        end
    end
    
    % plot standard error
    errorUpper = nanmean(inSpki,1)+(std(inSpki,0,1)/sqrt(size(inSpki,1)-1));
    errorLower = nanmean(inSpki,1)-(std(inSpki,0,1)/sqrt(size(inSpki,1)-1));
    h = fill([tt,fliplr(tt)],[errorUpper fliplr(errorLower)],[1 0.6 0.6]);
    set(h,'EdgeColor','none')
    errorUpper = nanmean(SCi0,1)+(std(SCi0,0,1)/sqrt(size(SCi0,1)-1));
    errorLower = nanmean(SCi0,1)-(std(SCi0,0,1)/sqrt(size(SCi0,1)-1));
    h = fill([tt,fliplr(tt)],[errorUpper fliplr(errorLower)],[0.8 0.8 0.8]);
    set(h,'EdgeColor','none')
    
    xlabel('time from target (ms) ')
    ylabel('spikes/s ')
    title('SCi signed rank test conditions')
    set(gca,'XTick', xmin:100:xmax);
    axis([-100 300 -0.04 ymax])
    
    avg_spks=nanmean(inSpki,1);
    plot(tt,avg_spks(1,:),file(1).colors{2});
    avg_spks=nanmean(SCi0,1);
    plot(tt,avg_spks(1,:),'k');
    
    ANOVAwin = [sigTimes(1) sigTimes(end)];
    ANOVAi = [];
    for y = 1:size(convspks,1)
        eval(['ANOVAi' num2str(y-1) '= nanmean(SCi' num2str(y-1) '(:,ANOVAwin(1):ANOVAwin(2)),2);']);
        eval(['errorI' num2str(y-1) '= (std(ANOVAi' num2str(y-1) ',0,1)/sqrt(size(ANOVAi' num2str(y-1) ',1)-1));']);
    end
    
    ANOVAiIN = [mean(ANOVAi0),mean(ANOVAi1),mean(ANOVAi3),mean(ANOVAi5)];
    ANOVAiOUT = [mean(ANOVAi0),mean(ANOVAi2),mean(ANOVAi4)];
    ANOVAsIN = [mean(ANOVAs0),mean(ANOVAs1),mean(ANOVAs3),mean(ANOVAs5)];
    ANOVAsOUT = [mean(ANOVAs0),mean(ANOVAs2),mean(ANOVAs4)];
    error = [errorS0, errorS1, errorS3, errorS5];
    figure
    bar(ANOVAsIN)
    hold on
    errorbar(ANOVAsIN, error, 'rx')
    ylim([0 0.5])
    figure
    error = [errorS0, errorS2, errorS4];
    bar(ANOVAsOUT)
    hold on
    errorbar(ANOVAsOUT, error, 'rx')
    ylim([0 0.5])
    figure
    error = [errorI0, errorI1, errorI3, errorI5];
    bar(ANOVAiIN)
    hold on
    errorbar(ANOVAiIN, error, 'rx')
    ylim([0 0.5])
    figure
    error = [errorI0, errorI2, errorI4];
    bar(ANOVAiOUT)
    hold on
    errorbar(ANOVAiOUT, error, 'rx')
    ylim([0 0.5])
    keyboard
end % end of neural stuff



%% plot pupil data for all files


if strcmpi(task,'prior04')
    cmb0 = [c0;c1]; % since pupil don't have RF, can combine all trials w/ same # of targets
    cmb1 = [c2;c3];
    cmb2 = c4;
    cmb3 = c5;
    cmb4 = [c6;c7];
    combo = 5;
    legName = {'1 high','1 low','2 high', '2 low', 'high low'};
else
    cmb0 = c0;
    cmb1 = [c1;c2];
    cmb2 = [c3;c4];
    cmb3 = c5;
    combo = 4;
    legName = {'0 pop','1 pop','2 pop','4 pop'};
    if strcmpi(task,'multipopout4v3')
        cmb4 = [c6;c7];
        cmb5 = [c8;c9];
        cmb6 = c10;
        combo = 7;
        legName = {'0 pop','1 pop','2 pop','4 pop','1 tar','2 tar','4 tar'};
    end
end


figure
hold on
xmin = -100;
xmax = 800;
tt=xmin:xmax-1;
for x = 1:combo
    eval(['avg_size=nanmean(cmb' num2str(x-1) ',1);']);
    h(x)=plot(tt,avg_size(1,:),file(1).colors{x});
end

ymax=ceil(max(max(avg_size,[],2)));
ymin=floor(min(min(avg_size,[],2)));

xlabel('time from target (ms) ')
ylabel('pupil size ')
title(['all files n=' num2str(n)]) 
legend(h,legName)
set(gca,'XTick', xmin:100:xmax);
axis([xmin xmax ymin ymax])