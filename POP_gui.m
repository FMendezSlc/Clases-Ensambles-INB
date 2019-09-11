function varargout = POP_gui(varargin)
% POP_GUI MATLAB code for POP_gui.fig
%      POP_GUI, by itself, creates a new POP_GUI or raises the existing
%      singleton*.
%
%      H = POP_GUI returns the handle to a new POP_GUI or the handle to
%      the existing singleton*.
%
%      POP_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POP_GUI.M with the given input arguments.
%
%      POP_GUI('Property','Value',...) creates a new POP_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before POP_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to POP_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help POP_gui

% Last Modified by GUIDE v2.5 11-Sep-2019 17:44:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @POP_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @POP_gui_OutputFcn, ...
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


% --- Executes just before POP_gui is made visible.
function POP_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to POP_gui (see VARARGIN)

% Choose default command line output for POP_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes POP_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = POP_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in kmeans.
function kmeans_Callback(hObject, eventdata, handles)
% hObject    handle to kmeans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eval_clust = evalclusters(handles.score(:,:),'kmeans','CalinskiHarabasz','KList',[2:20]);
idx = kmeans(handles.score, eval_clust.OptimalK);
axes(handles.MainAxes)
gscatter(handles.score(:,1), handles.score(:,2), idx)
hold off
handles.idx = idx;
guidata(hObject, handles)

% --- Executes on button press in gmdistfit.
function gmdistfit_Callback(hObject, eventdata, handles)
% hObject    handle to gmdistfit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eval_clust = evalclusters(handles.score(:,1:2),'gmdistribution','CalinskiHarabasz','KList',[2:20]);
options = statset('MaxIter',1000);
GMModels = fitgmdist(handles.score(:,1:2),eval_clust.OptimalK,'Options',options,'RegularizationValue',0.1);
axes(handles.MainAxes)
gscatter(handles.score(:,1), handles.score(:,2)), hold on
h = gca;
gmPDF = @(x1,x2)reshape(pdf(GMModels,[x1(:) x2(:)]),size(x1));
fcontour(gmPDF,[h.XLim h.YLim],'MeshDensity',100)
hold off
idx = cluster(GMModels,handles.score(:,1:2));
handles.idx = idx;
guidata(hObject, handles)

% --- Executes on button press in pcaDimRed.
function pcaDimRed_Callback(hObject, eventdata, handles)
% hObject    handle to pcaDimRed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[coeff,score,latent,tsquared,explained] = pca(handles.raw_data');
axes(handles.MainAxes)
gscatter(score(:,1), score(:,2))
legend('off');
hold off
handles.score = score(:,1:3);
guidata(hObject, handles)


% --- Executes on button press in tsneDimRed.
function tsneDimRed_Callback(hObject, eventdata, handles)
% hObject    handle to tsneDimRed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
score = tsne(handles.raw_data');
axes(handles.MainAxes)
gscatter(score(:,1), score(:,2))
legend('off');
handles.score = score;
guidata(hObject, handles)
hold off

% --- Executes on button press in loadFile.
function loadFile_Callback(hObject, eventdata, handles)
% hObject    handle to loadFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, path] = uigetfile('*.mat',...
               'Select an icon file','Spikes_significant.mat');
load([path file]);
handles.raw_data = Spikes;
guidata(hObject, handles)
cla(handles.MainAxes, 'reset')
axes(handles.MainAxes)
imagesc(Spikes)
hold off


% --- Executes on button press in showClusters.
function showClusters_Callback(hObject, eventdata, handles)
% hObject    handle to showClusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure()
for ii = 1:length(unique(handles.idx))
subplot(1, length(unique(handles.idx)), ii)
imagesc(handles.raw_data(:,handles.idx==ii))
title(sprintf('Cluster %i', ii))
end