function [pupil_m,tt] = behavior(file)

%% define all my event codes

% these codes can be used as is. Presented in approx. the order they appear.
START = 1001;
PARADIGM = 32006;
PAUSECD = 1003; % I think?
NSTIM = 196;
SENDSTIM = 1060;
CDEDIVIDER = 1062;
SENDSTIMEND = 1061;
SENDORI = 1070;
SENDORIEND = 1071;
SENDCOL = 1080;
SENDCOLEND = 1081;
WINDOPN = 800;
BCKLIT = 1027;

% these codes are added to the condition number 6xyy, where x are the
% following codes, yy is the stimulus number.
FPONSET = 200;
EYEINON = 400;
ARRAYON = 600;
FPOFF = 800;
SACON = 1000;

% back to standalone codes
REWON = 1030;
REWOFF = 1040;

% add this to whatever basecode was there (ie. add -2000 to 6xyy to
% get 4xyy)
ABORT = -2000;

numConds = length(file.condNames);

%% open file
try
    load([file.path,file.name],'r_strobed','EL_P','EL_H','EL_V')
catch
    try
        readall([file.path,file.name(1:end-6),'.plx'])
        load([file.path,file.name],'r_strobed','EL_P','EL_H','EL_V')
    catch
        pupil_m = 'corrupted';
        tt = 'corrupted';
        return
    end
end

%% dividing trial by events
AllEvts=r_strobed(r_strobed(:,2)==START,:);% 1005

j=1;
for i=1:length(AllEvts)-1,    
    trlcd_tmp=r_strobed(r_strobed(:,1)>AllEvts(i)-1 & r_strobed(:,1)<AllEvts(i+1),:); % relevant time stamps for current trial
    trlcd_tmp(trlcd_tmp(:,2)==PAUSECD,:)=[];
    if ~isempty(trlcd_tmp(trlcd_tmp(:,2)==REWON,1)), %was there a reward  
        % START
        try
            TrlCodes{j}=trlcd_tmp;
            TrlStart_ts(j,:)=trlcd_tmp(1,1); %trial start
            
            RFxy(j,:) = [trlcd_tmp(find(trlcd_tmp(:,2)==SENDSTIM)+1,2)-10000, trlcd_tmp(find(trlcd_tmp(:,2)==SENDSTIM)+2,2)-10000] / 10; % convert to degrees
            Cond(j,:) = trlcd_tmp(find(trlcd_tmp(:,2)==PARADIGM,1,'last')+1,2);
            if Cond(j,:) < 6000 || Cond(j,:) > 7000 % don't know why sometimes conditions don't start with 6xxx.
                Cond(j,:) = [];
            end
            
            % EVENTS
%             Fixon_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)==(Cond(j)+FPONSET), 1);
%             Fixon_codes(j,:) = trlcd_tmp(trlcd_tmp(:,2)==(Cond(j)+FPONSET), 2);
%             Eyein_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)==(Cond(j)+EYEINON), 1);
%             Eyein_codes(j,:) = trlcd_tmp(trlcd_tmp(:,2)==(Cond(j)+EYEINON), 2);
            Ton_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)==(Cond(j)+ARRAYON), 1);
            Ton_codes(j,:) = trlcd_tmp(trlcd_tmp(:,2)==(Cond(j)+ARRAYON), 2);
%             Fixoff_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)==(Cond(j)+FPOFF),1);
%             Fixoff_codes(j,:) = trlcd_tmp(trlcd_tmp(:,2)==(Cond(j)+FPOFF),2);
%             Sacon_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)==(Cond(j)+SACON),1);
%             Sacon_codes(j,:) = trlcd_tmp(trlcd_tmp(:,2)==(Cond(j)+SACON),2);
%             Rewon_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)==REWON,1);
%             Rewon_codes(j,:) = trlcd_tmp(trlcd_tmp(:,2)==REWON,2);
%             Rewoff_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)==REWOFF,1);
%             Rewoff_codes(j,:) = trlcd_tmp(trlcd_tmp(:,2)==REWOFF,2);
            j=j+1;
        catch
        end
    end     
    clear trlcd_tmp; 
end

%% figure out the condition of each file
try
    condition = rem(Cond,numConds);
    condition(~condition)=numConds; % make condition 5/6 instead of "condition 0"
catch
    pupil_m = 'wrong format?';
    tt = 'wrong format?';
    return
end

for x = 1:numConds
    if sum(condition==x)<1
        pupil_m = 'not enough';
        tt = 'not enough';
        return
    end
end

rewCond = condition;
if strcmpi(file.task,'prior04')
    if file.rew1 < file.rew2 % make rewCond relative to reward size, not color
        rewCond(condition==1) = 3;
        rewCond(condition==2) = 4;
        rewCond(condition==3) = 1;
        rewCond(condition==4) = 2;
        rewCond(condition==5) = 6;
        rewCond(condition==6) = 5;
        rewCond(condition==7) = 8;
        rewCond(condition==8) = 7;
    end
end

%% correct for vertical refresh (LCD 60 Hz, 16.67ms per refresh)
vertical_screen_deg = 50;
refresh_ms=16.67;
vertical_refresh_correction=((vertical_screen_deg - ((vertical_screen_deg / 2) - (RFxy(1,2) * .1))) / vertical_screen_deg) * refresh_ms; 
vertical_refresh_correction=ceil(vertical_refresh_correction);

%% plot align on target
xmin = -100;
xmax = 800; % time window to display spikes
tt=xmin:xmax-1;
pupil_m=[];
b=Ton_ts - vertical_refresh_correction;

figure
hold on
for t=1:length(b),     
    pupildata = EL_P(b(t)+xmin:b(t)+xmax-1)';
    baseline = mean(pupildata(abs(xmin)-20:abs(xmin)+150),2);
    eval(['c' num2str(rewCond(t)-1) '(t, :)=pupildata-baseline;']);   
end

for x=1:numConds
    eval(['pupil_m(x,:)=mean(c' num2str(x-1) ',1);']);

    h(x) = plot(tt,pupil_m(x,:),file.colors{x});
%     set(h(x),'color',colors(x))
end
ymax1=ceil(max(max(pupil_m,[],2)))+1;
ymin1=floor(min(min(pupil_m,[],2)))-1;

xlabel('time from target (ms) ')
ylabel('spikes/s ')
title([file.name(1:end-6) ', x = ' num2str(RFxy(1,1)) ' , y = ' num2str(RFxy(1,2)) ', n = ' num2str(t)]) 
legend(h,file.condNames)
set(gca,'XTick', xmin:100:xmax);
axis([xmin xmax ymin1 ymax1])


vector = sqrt(RFxy(1,1)^2 + RFxy(1,2)^2);
if vector<10
    pupil_m = 'small RF';
    tt = 'small RF';
end
