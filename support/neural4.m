% in this version i'm trying to include analysis for priority task in here
% as well. might get too messy and have to split into 2 separate programs
%
% nevermind, i think i'm gonna make this file for priority task only.
% they're too different from the multipopout task

function [convspks_delay_m, lfpOUT] = neural(file, useOrig)

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
SACEND = 1400;

% back to standalone codes
REWON = 1030;
REWOFF = 1040;

% add this to whatever basecode was there (ie. add -2000 to 6xyy to
% get 4xyy)
ABORT = -2000;

%% open file
channel = ['sig001',file.unit];
try
    load([file.path,file.name],'r_strobed',channel,'AD_heye','AD_veye','AD01')
catch
    try
        readall([file.path,file.name(1:end-6),'.plx'])
        load([file.path,file.name],'r_strobed',channel,'AD_heye','AD_veye','AD01')
    catch
        convspks_delay_m = 'corrupted';
        lfpOUT = 'corrupted';
        return
    end
end


if useOrig
    try
        load([file.path,file.name(1:end-6),'.mat'],'org_result')
    catch
        try
            findTrials([file.path,file.name])
            load([file.path,file.name(1:end-6),'.mat'],'org_result')
        catch
            convspks_delay_m = 'corrupted';
            lfpOUT = 'corrupted';
            return
        end
    end
    trialInfo = org_result;
else
    try
        load([file.path,file.name(1:end-6),'.mat'],'result')
    catch
        disp([file.name(1:end-6) ' has no eyetrace check file.'])
        try
            findTrials([file.path,file.name])
            load([file.path,file.name(1:end-6),'.mat'],'result')
        catch
            convspks_delay_m = 'corrupted';
            lfpOUT = 'corrupted';
            return
        end
    end
    trialInfo = result;
end

% %% dividing trial by events
% AllEvts=r_strobed(r_strobed(:,2)==1001,:);% 1005
% 
% j=1;
% for i=1:length(AllEvts)-1, % for each trial
%     trlcd_tmp=r_strobed(r_strobed(:,1)>AllEvts(i)-1 & r_strobed(:,1)<AllEvts(i+1),:); % relevant time stamps for current trial
%     trlcd_tmp(trlcd_tmp(:,2)==PAUSECD,:)=[];
%     if ~isempty(trlcd_tmp(trlcd_tmp(:,2)==REWON,1)), %was there a reward  
%         % START
%         try
%             TrlCodes{j}=trlcd_tmp;
%             TrlStart_ts(j,:)=trlcd_tmp(1,1); %trial start
%             
%             RFxy(j,:) = [trlcd_tmp(find(trlcd_tmp(:,2)==SENDSTIM)+1,2)-10000, trlcd_tmp(find(trlcd_tmp(:,2)==SENDSTIM)+2,2)-10000] / 10; % convert to degrees
%             Cond(j,:) = trlcd_tmp(find(trlcd_tmp(:,2)==PARADIGM,1,'last')+1,2);
%             if Cond(j,:) < 6000 || Cond(j,:) > 7000
%                 Cond(j,:) = [];
%                 continue
%             else
%                 condCorr(j) = 6000;
%             end
%             
%             % EVENTS
%             Ton_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)==(Cond(j)+ARRAYON), 1);
%             Rewon_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)==REWON,1);
%             
%             j=j+1;
%         catch
%         end
%     end     
%     clear trlcd_tmp; 
% end
% 
% %% figure out the condition of each file
% try
%     numConds = length(file.condNames);
%     condition = rem(Cond-condCorr',numConds);
%     condition(~condition)=numConds; % make condition 5/6 instead of "condition 0"
% catch
%     convspks_delay_m = 'wrong format?';
%     lfpOUT = 'wrong format?';
%     return
% end

%% deleting trials that was manually counted as invalid
delIndex = cellfun(@isempty,{trialInfo.SacON});
trialInfo(delIndex) = [];
SacOFF = [trialInfo.SacOFF]';
Ton_ts = [trialInfo.TarON]';
FPoff = [trialInfo.FPoff];
RFxy = [trialInfo(1).RFx,trialInfo(1).RFy];
condition = [trialInfo.condition]';
numConds = length(file.condNames);

for x = 1:numConds
    if sum(condition==x)<1
        convspks_delay_m = 'not enough';
        lfpOUT = 'not enough';
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
        
        rewCond(condition==9) = 11;
        rewCond(condition==10) = 12;
        rewCond(condition==11) = 9;
        rewCond(condition==12) = 10;
        rewCond(condition==13) = 14;
        rewCond(condition==14) = 13;
        rewCond(condition==15) = 16;
        rewCond(condition==16) = 15;
    end
end

%% correct for vertical refresh (LCD 60 Hz, 16.67ms per refresh)
vertical_screen_deg = 50;
refresh_ms=16.67;
vertical_refresh_correction=((vertical_screen_deg - ((vertical_screen_deg / 2) - (RFxy(1,2) * .1))) / vertical_screen_deg) * refresh_ms; 
vertical_refresh_correction=ceil(vertical_refresh_correction);

%% plot align on target
xmin = -200;ymin=0;
xmax = 800; % time window to display spikes
sigma=5; nstd=sigma*5;
gauswin=normpdf(-nstd:nstd,0,sigma);
tt=xmin-nstd:xmax+nstd-1;
winsize = xmax - xmin;
spks_m=zeros(1,winsize);convspks_delay_m=[];
b=Ton_ts - vertical_refresh_correction;
for x=1:numConds
    eval(['c' num2str(x-1) ' =[];']);
    eval(['lfp' num2str(x-1) ' =[];']);
end
tick_st=0;

% nice function to get rid of HF noise
heye=smooth(AD_heye,50);
veye=smooth(AD_veye,50);

% filting lfp
d = designfilt('lowpassiir',...
    'PassbandFrequency',50,...
    'StopbandFrequency',55,...
    'DesignMethod','butter', ...         % Design method
    'MatchExactly','passband',...
    'SampleRate',1000);
lfpFilt = filtfilt(d,AD01);

for t=1:length(b),     
    
%     % determining where saccade landed to determine which reward was given
%     
%     % find where saccade landed
%     saccX = heye(SacOFF(t));
%     saccY = veye(SacOFF(t));
%     fromRF = sqrt((saccX-RFxy(1,1))^2 + (saccY-RFxy(1,2))^2);
%     
%     figure
%     plot(-100:500,heye(FPoff(t)-100:FPoff(t)+500));
%     hold on
%     plot(-100:500,veye(FPoff(t)-100:FPoff(t)+500),'m');
%     plot([-100 500],[RFxy(1,1) RFxy(1,1)],'b')
%     plot([-100 500],[RFxy(1,2) RFxy(1,2)],'m')
%     plot(SacOFF(t)-FPoff(t),saccX,'bo')
%     plot(SacOFF(t)-FPoff(t),saccY,'mo')
%     title(['condition: ' num2str(condition(t))])
%     hold off
%     
%     figure
%     ezplot(@(x,y)myfun(x,y,RFxy(1,1),RFxy(1,2)),[RFxy(1,1)-8 RFxy(1,1)+8 RFxy(1,2)-8 RFxy(1,2)+8])
%     hold on
%     plot(heye(FPoff(t)-100:FPoff(t)+500),veye(FPoff(t)-100:FPoff(t)+500),'k')
%     plot(saccX,saccY,'rx');
%     hold off
%     axis([-30 30 -30 30])
%     axis square
%     
%     keyboard
%     
%     if fromRF<5
%         inRF = true;
%     else
%         inRF = false;
%     end
    
    eval(['spks = ' channel '(' channel '< b(t)+xmax & ' channel '> b(t)+xmin) - b(t)-xmin;']);
    tempspk =zeros(1,winsize);
    tempspk(1,spks)=1;
    eval(['c' num2str(rewCond(t)-1) '=[c' num2str(rewCond(t)-1) '; tempspk];']);
    
%     baseLFP = mean(lfpFilt(b(t)+(xmin-10):b(t)+(xmin+10)));
    baseLFP = mean(lfpFilt(b(t)-50:b(t)));
    lfp = lfpFilt(b(t)+xmin:b(t)+xmax-1)';
    lfp = lfp-baseLFP;
    
    eval(['lfp' num2str(rewCond(t)-1) '=[lfp' num2str(rewCond(t)-1) '; lfp];']);
end



%% plotting rasters
figure(3)
subplot(2,1,1)
hold on
for x=1:numConds
    eval(['curCond = c' num2str(x-1) ';']);
    for y=1:size(curCond,1)
        tick_st=tick_st+4;
        spks = find(curCond(y,:));
        if ~isempty(spks),
            for m=1:size(spks,2)
                h=plot([spks(m)+xmin spks(m)+xmin], [tick_st tick_st+2],file.colors{x});
                set(h,'linewidth',1.5)
            end
            clear h ;
        end
    end
end


%% plotting spike density function
subplot(2,1,2)
hold on
clear spks_m
for x=1:numConds
    eval(['spks_m(1,:)=mean(c' num2str(x-1) ',1);']);
    % convolve with gaussian
    convspks_delay_m(x,:)=conv(spks_m(1,:)*1000, gauswin);
    
    % spike density function
    h(x)=plot(tt,convspks_delay_m(x,:),file.colors{x});
    clear spks_m
    
    eval(['lfpOUT(x,:)=mean(lfp' num2str(x-1) ',1);']);
end
ymax1=ceil(max(max(convspks_delay_m,[],2)))+10;

if numConds == 5
    extraCond = nan(1,length(convspks_delay_m));
    convspks_delay_m = [extraCond;convspks_delay_m];
end

xlabel('time from target (ms) ')
ylabel('spikes/s ')
title([file.name(1:end-6) ' depth' num2str(file.depth) ', x = ' num2str(RFxy(1,1)) ' , y = ' num2str(RFxy(1,2)) ', isolation: ' file.isolation]) 
legend(h,file.condNames)
set(gca,'XTick', xmin:100:xmax);
axis([xmin xmax 0 30])
close all

end

function z = myfun(x,y,RFx,RFy)
z = (x-RFx)^2/5^2+(y-RFy)^2/5^2-1;
end