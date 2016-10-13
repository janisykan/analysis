function varargout = Manual_MS_correction(varargin)
% MANUAL_MS_CORRECTION MATLAB code for Manual_MS_correction.fig
%
%      THIS IS FOR MANUALLY CHECKING SACCADES IN PRIORITY TASK. Adapted
%      from MingPo Yang's (Dorris Lab @ SIBS) code for checking
%      microsaccade. (Updated 2016 by Janis Kan)
%
%      MANUAL_MS_CORRECTION, by itself, creates a new MANUAL_MS_CORRECTION or raises the existing
%      singleton*.
%
%      H = MANUAL_MS_CORRECTION returns the handle to a new MANUAL_MS_CORRECTION or the handle to
%      the existing singleton*.
%
%      MANUAL_MS_CORRECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MANUAL_MS_CORRECTION.M with the given input arguments.
%
%      MANUAL_MS_CORRECTION('Property','Value',...) creates a new MANUAL_MS_CORRECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Manual_MS_correction_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Manual_MS_correction_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Manual_MS_correction

% Last Modified by GUIDE v2.5 15-Jun-2016 15:26:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Manual_MS_correction_OpeningFcn, ...
                   'gui_OutputFcn',  @Manual_MS_correction_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
               
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Manual_MS_correction is made visible.
function Manual_MS_correction_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Manual_MS_correction (see VARARGIN)

% Choose default command line output for Manual_MS_correction
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Manual_MS_correction wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = Manual_MS_correction_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Load.
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename=uigetfile('*.mat');
if filename==0
    return;
end
load(filename);
color = 'krg';
setappdata(gcf,'filename',filename);
setappdata(gcf,'org_result',org_result);
setappdata(gcf,'result',result);
setappdata(gcf,'trialind',trialind);
cla(handles.Htrace);
cla(handles.Vtrace);
plot(handles.Htrace,result(trialind).time-result(trialind).FPoff,result(trialind).xpos);
plot(handles.Vtrace,result(trialind).time-result(trialind).FPoff,result(trialind).ypos);
plot(handles.Htrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[result(trialind).RFx, result(trialind).RFx],color(result(trialind).RFid+1));
plot(handles.Htrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[-result(trialind).RFx, -result(trialind).RFx],color(result(trialind).antiRFid+1));
plot(handles.Vtrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[result(trialind).RFy, result(trialind).RFy],color(result(trialind).RFid+1));
plot(handles.Vtrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[-result(trialind).RFy, -result(trialind).RFy],color(result(trialind).antiRFid+1));
if ~isempty(result(trialind).SacON)
    index=result(trialind).SacON:result(trialind).SacOFF;
    index = index - result(trialind).time(1)+1;
    plot(handles.Htrace,result(trialind).time(index)-result(trialind).FPoff,result(trialind).xpos(index),'color','g','linewidth',2);
    plot(handles.Vtrace,result(trialind).time(index)-result(trialind).FPoff,result(trialind).ypos(index),'color','g','linewidth',2);
end
set(handles.hortitle,'string',[filename(1:end-4) ' HorTrace' num2str(trialind)]);
set(handles.vertitle,'string',[filename(1:end-4) ' VerTrace' num2str(trialind)]);

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename=getappdata(gcf,'filename');
org_result=getappdata(gcf,'org_result');
result=getappdata(gcf,'result');
trialind=getappdata(gcf,'trialind');
if isempty(result)
    return;
end
save(filename,'org_result','result','trialind');

% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename=getappdata(gcf,'filename');
org_result=getappdata(gcf,'org_result');
result=getappdata(gcf,'result');
trialind=getappdata(gcf,'trialind');

color = 'krg';
if isempty(result)
    return;
end
if mod(trialind,20)==0 % autosave after every 20 trials
    save(filename,'org_result','result','trialind');
end
if trialind<length(result)
    trialind=trialind+1;
    cla(handles.Htrace);
    cla(handles.Vtrace);
    plot(handles.Htrace,result(trialind).time-result(trialind).FPoff,result(trialind).xpos);
    plot(handles.Vtrace,result(trialind).time-result(trialind).FPoff,result(trialind).ypos);
    plot(handles.Htrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[result(trialind).RFx, result(trialind).RFx],color(result(trialind).RFid+1));
    plot(handles.Htrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[-result(trialind).RFx, -result(trialind).RFx],color(result(trialind).antiRFid+1));
    plot(handles.Vtrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[result(trialind).RFy, result(trialind).RFy],color(result(trialind).RFid+1));
    plot(handles.Vtrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[-result(trialind).RFy, -result(trialind).RFy],color(result(trialind).antiRFid+1));
    if ~isempty(result(trialind).SacON)
        index=result(trialind).SacON:result(trialind).SacOFF;
        index = index - result(trialind).time(1)+1;
        plot(handles.Htrace,result(trialind).time(index)-result(trialind).FPoff,result(trialind).xpos(index),'color','g','linewidth',2);
        plot(handles.Vtrace,result(trialind).time(index)-result(trialind).FPoff,result(trialind).ypos(index),'color','g','linewidth',2);
    end
    set(handles.hortitle,'string',[filename(1:end-4) ' HorTrace' num2str(trialind)]);
    set(handles.vertitle,'string',[filename(1:end-4) ' VerTrace' num2str(trialind)]);
    setappdata(gcf,'trialind',trialind);
else
    save(filename,'org_result','result','trialind');
end

% --- Executes on button press in del.
function del_Callback(hObject, eventdata, handles)
% hObject    handle to del (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
color = 'krg';
result=getappdata(gcf,'result');
trialind=getappdata(gcf,'trialind');
if isempty(result)
    return;
end

result(trialind).SacON=[];
result(trialind).SacOFF=[];

cla(handles.Htrace);
cla(handles.Vtrace);
plot(handles.Htrace,result(trialind).time-result(trialind).FPoff,result(trialind).xpos);
plot(handles.Vtrace,result(trialind).time-result(trialind).FPoff,result(trialind).ypos);
plot(handles.Htrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[result(trialind).RFx, result(trialind).RFx],color(result(trialind).RFid+1));
plot(handles.Htrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[-result(trialind).RFx, -result(trialind).RFx],color(result(trialind).antiRFid+1));
plot(handles.Vtrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[result(trialind).RFy, result(trialind).RFy],color(result(trialind).RFid+1));
plot(handles.Vtrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[-result(trialind).RFy, -result(trialind).RFy],color(result(trialind).antiRFid+1));
if ~isempty(result(trialind).SacON)
    index=result(trialind).SacON:result(trialind).SacOFF;
    index = index - result(trialind).time(1)+1;
    plot(handles.Htrace,result(trialind).time(index)-result(trialind).FPoff,result(trialind).xpos(index),'color','g','linewidth',2);
    plot(handles.Vtrace,result(trialind).time(index)-result(trialind).FPoff,result(trialind).ypos(index),'color','g','linewidth',2);
end
% if ~isempty(result(trialind).param)
%     param=result(trialind).param;
%     for k=1:size(param,2)
%         index=[];
%         index=param(1,k):param(2,k);
%         plot(handles.Htrace,result(trialind).time(index),result(trialind).xpos(index),'color','g','linewidth',2);
%         plot(handles.Vtrace,result(trialind).time(index),result(trialind).ypos(index),'color','g','linewidth',2);
%     end
% end

setappdata(gcf,'result',result);
setappdata(gcf,'trialind',trialind);

% --- Executes on button press in add.
function add_Callback(hObject, eventdata, handles)
% hObject    handle to add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result=getappdata(gcf,'result');
trialind=getappdata(gcf,'trialind');
color = 'krg';
if isempty(result)
    return;
end
dcm_obj=datacursormode(gcf);
c_info = getCursorInfo(dcm_obj);
if isempty(c_info)
    return;
end
for k=1:length(c_info)
    cursorpos(k)=c_info(k).Position(1);
end
cursorpos=sort(cursorpos);
% dt=result(1).dt*1000;%change the units to ms, because the time of abscissa is ms
% cursorpos=round(cursorpos/dt);%change the cursorpos to index,because msst&mset is index
if mod(length(cursorpos),2)==1
    cursorpos=[cursorpos length(result(1).time)];
end
for k=1:length(cursorpos)/2
    msst(k)=cursorpos(2*k-1);
    mset(k)=cursorpos(2*k);
end

result(trialind).SacON=msst+result(trialind).FPoff;
result(trialind).SacOFF=mset+result(trialind).FPoff;

cla(handles.Htrace);
cla(handles.Vtrace);
plot(handles.Htrace,result(trialind).time-result(trialind).FPoff,result(trialind).xpos);
plot(handles.Vtrace,result(trialind).time-result(trialind).FPoff,result(trialind).ypos);
plot(handles.Htrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[result(trialind).RFx, result(trialind).RFx],color(result(trialind).RFid+1));
plot(handles.Htrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[-result(trialind).RFx, -result(trialind).RFx],color(result(trialind).antiRFid+1));
plot(handles.Vtrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[result(trialind).RFy, result(trialind).RFy],color(result(trialind).RFid+1));
plot(handles.Vtrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[-result(trialind).RFy, -result(trialind).RFy],color(result(trialind).antiRFid+1));
if ~isempty(result(trialind).SacON)
    index=result(trialind).SacON:result(trialind).SacOFF;
    index = index - result(trialind).time(1)+1;
    plot(handles.Htrace,result(trialind).time(index)-result(trialind).FPoff,result(trialind).xpos(index),'color','g','linewidth',2);
    plot(handles.Vtrace,result(trialind).time(index)-result(trialind).FPoff,result(trialind).ypos(index),'color','g','linewidth',2);
end

setappdata(gcf,'result',result);
setappdata(gcf,'trialind',trialind);

% --- Executes on button press in compare.
function compare_Callback(hObject, eventdata, handles)
% hObject    handle to compare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename=getappdata(gcf,'filename');
org_result=getappdata(gcf,'org_result');
result=getappdata(gcf,'result');
trialind=getappdata(gcf,'trialind');
color = 'krg';
if isempty(result)
    return;
end
cla(handles.Htrace);
cla(handles.Vtrace);
plot(handles.Htrace,result(trialind).time-result(trialind).FPoff,result(trialind).xpos);
plot(handles.Vtrace,result(trialind).time-result(trialind).FPoff,result(trialind).ypos);
plot(handles.Htrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[result(trialind).RFx, result(trialind).RFx],color(result(trialind).RFid+1));
plot(handles.Htrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[-result(trialind).RFx, -result(trialind).RFx],color(result(trialind).antiRFid+1));
plot(handles.Vtrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[result(trialind).RFy, result(trialind).RFy],color(result(trialind).RFid+1));
plot(handles.Vtrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[-result(trialind).RFy, -result(trialind).RFy],color(result(trialind).antiRFid+1));
if ~isempty(result(trialind).SacON)
    index=result(trialind).SacON:result(trialind).SacOFF;
    index = index - result(trialind).time(1)+1;
    plot(handles.Htrace,result(trialind).time(index)-result(trialind).FPoff,result(trialind).xpos(index),'color','g','linewidth',2);
    plot(handles.Vtrace,result(trialind).time(index)-result(trialind).FPoff,result(trialind).ypos(index),'color','g','linewidth',2);
end

if ~isempty(org_result(trialind).SacON)
    index=org_result(trialind).SacON:org_result(trialind).SacOFF;
    index = index - org_result(trialind).time(1)+1;
    plot(handles.Htrace,org_result(trialind).time(index)-org_result(trialind).FPoff,org_result(trialind).xpos(index),'color','r','linestyle','--','linewidth',4);
    plot(handles.Vtrace,org_result(trialind).time(index)-org_result(trialind).FPoff,org_result(trialind).ypos(index),'color','r','linestyle','--','linewidth',4);
end

set(handles.hortitle,'string',[filename(1:end-4) ' HorTrace' num2str(trialind)]);
set(handles.vertitle,'string',[filename(1:end-4) ' VerTrace' num2str(trialind)]);


function Jump_Callback(hObject, eventdata, handles)
% hObject    handle to Jump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Jump as text
%        str2double(get(hObject,'String')) returns contents of Jump as a double
filename=getappdata(gcf,'filename');
result=getappdata(gcf,'result');
color = 'krg';
if isempty(result)
    return;
end
trialind=round(str2double(get(hObject,'String')));
if ~isnumeric(trialind)||isnan(trialind)||trialind<1||trialind>length(result)
    set(hObject,'String','Need numbers!');
    trialind=getappdata(gcf,'trialind');
end
cla(handles.Htrace);
cla(handles.Vtrace);
plot(handles.Htrace,result(trialind).time-result(trialind).FPoff,result(trialind).xpos);
plot(handles.Vtrace,result(trialind).time-result(trialind).FPoff,result(trialind).ypos);
plot(handles.Htrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[result(trialind).RFx, result(trialind).RFx],color(result(trialind).RFid+1));
plot(handles.Htrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[-result(trialind).RFx, -result(trialind).RFx],color(result(trialind).antiRFid+1));
plot(handles.Vtrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[result(trialind).RFy, result(trialind).RFy],color(result(trialind).RFid+1));
plot(handles.Vtrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[-result(trialind).RFy, -result(trialind).RFy],color(result(trialind).antiRFid+1));
if ~isempty(result(trialind).SacON)
    index=result(trialind).SacON:result(trialind).SacOFF;
    index = index - result(trialind).time(1)+1;
    plot(handles.Htrace,result(trialind).time(index)-result(trialind).FPoff,result(trialind).xpos(index),'color','g','linewidth',2);
    plot(handles.Vtrace,result(trialind).time(index)-result(trialind).FPoff,result(trialind).ypos(index),'color','g','linewidth',2);
end
set(handles.hortitle,'string',[filename(1:end-4) ' HorTrace' num2str(trialind)]);
set(handles.vertitle,'string',[filename(1:end-4) ' VerTrace' num2str(trialind)]);
setappdata(gcf,'result',result);
setappdata(gcf,'trialind',trialind);

% --- Executes during object creation, after setting all properties.
function Jump_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Jump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in previous.
function previous_Callback(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename=getappdata(gcf,'filename');
result=getappdata(gcf,'result');
trialind=getappdata(gcf,'trialind');
color = 'krg';
if isempty(result)
    return;
end
trialind=trialind-1;
if trialind<=0
    trialind=1;
end
cla(handles.Htrace);
cla(handles.Vtrace);
plot(handles.Htrace,result(trialind).time-result(trialind).FPoff,result(trialind).xpos);
plot(handles.Vtrace,result(trialind).time-result(trialind).FPoff,result(trialind).ypos);
plot(handles.Htrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[result(trialind).RFx, result(trialind).RFx],color(result(trialind).RFid+1));
plot(handles.Htrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[-result(trialind).RFx, -result(trialind).RFx],color(result(trialind).antiRFid+1));
plot(handles.Vtrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[result(trialind).RFy, result(trialind).RFy],color(result(trialind).RFid+1));
plot(handles.Vtrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[-result(trialind).RFy, -result(trialind).RFy],color(result(trialind).antiRFid+1));
if ~isempty(result(trialind).SacON)
    index=result(trialind).SacON:result(trialind).SacOFF;
    index = index - result(trialind).time(1)+1;
    plot(handles.Htrace,result(trialind).time(index)-result(trialind).FPoff,result(trialind).xpos(index),'color','g','linewidth',2);
    plot(handles.Vtrace,result(trialind).time(index)-result(trialind).FPoff,result(trialind).ypos(index),'color','g','linewidth',2);
end
set(handles.hortitle,'string',[filename(1:end-4) ' HorTrace' num2str(trialind)]);
set(handles.vertitle,'string',[filename(1:end-4) ' VerTrace' num2str(trialind)]);
setappdata(gcf,'result',result);
setappdata(gcf,'trialind',trialind);


% --- Executes on button press in Revert.
function Revert_Callback(hObject, eventdata, handles)
% hObject    handle to Revert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename=getappdata(gcf,'filename');
org_result=getappdata(gcf,'org_result');
result=getappdata(gcf,'result');
trialind=getappdata(gcf,'trialind');
color = 'krg';
if isempty(result)
    return;
end

result(trialind).SacON=org_result(trialind).SacON;
result(trialind).SacOFF=org_result(trialind).SacOFF;

cla(handles.Htrace);
cla(handles.Vtrace);
plot(handles.Htrace,result(trialind).time-result(trialind).FPoff,result(trialind).xpos);
plot(handles.Vtrace,result(trialind).time-result(trialind).FPoff,result(trialind).ypos);
plot(handles.Htrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[result(trialind).RFx, result(trialind).RFx],color(result(trialind).RFid+1));
plot(handles.Htrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[-result(trialind).RFx, -result(trialind).RFx],color(result(trialind).antiRFid+1));
plot(handles.Vtrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[result(trialind).RFy, result(trialind).RFy],color(result(trialind).RFid+1));
plot(handles.Vtrace,[result(trialind).time(1)-result(trialind).FPoff,result(trialind).time(end)-result(trialind).FPoff],[-result(trialind).RFy, -result(trialind).RFy],color(result(trialind).antiRFid+1));
if ~isempty(result(trialind).SacON)
    index=result(trialind).SacON:result(trialind).SacOFF;
    index = index - result(trialind).time(1)+1;
    plot(handles.Htrace,result(trialind).time(index)-result(trialind).FPoff,result(trialind).xpos(index),'color','g','linewidth',2);
    plot(handles.Vtrace,result(trialind).time(index)-result(trialind).FPoff,result(trialind).ypos(index),'color','g','linewidth',2);
end

set(handles.hortitle,'string',[filename(1:end-4) ' HorTrace' num2str(trialind)]);
set(handles.vertitle,'string',[filename(1:end-4) ' VerTrace' num2str(trialind)]);

setappdata(gcf,'result',result);
setappdata(gcf,'trialind',trialind);