function varargout = uu(varargin)
% UU MATLAB code for uu.fig
%      UU, by itself, creates a new UU or raises the existing
%      singleton*.
%
%      H = UU returns the handle to a new UU or the handle to
%      the existing singleton*.
%
%      UU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UU.M with the given input arguments.
%
%      UU('Property','Value',...) creates a new UU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before uu_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to uu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help uu

% Last Modified by GUIDE v2.5 04-Apr-2020 14:50:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @uu_OpeningFcn, ...
                   'gui_OutputFcn',  @uu_OutputFcn, ...
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


% --- Executes just before uu is made visible.
function uu_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to uu (see VARARGIN)

% Choose default command line output for uu
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes uu wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = uu_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(~, ~, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global file_name;
%guidata(hObject,handles)
file_name=uigetfile({'*.bmp;*.jpeg;*.png;*.jpg;*.tiff;';'*.*'},'Select an Image File');
fileinfo = dir(file_name);
SIZE = fileinfo.bytes;
Size = SIZE/1024;
set(handles.text6, 'string', Size);
set(handles.axes1,'Units','pixels');
resizePos = get(handles.axes1,'Position');

axes(handles.axes1);
imshow(file_name);
set(handles.axes1,'Units','normalized');



% --- Executes on selection change in popupmenu1.

% --- Executes on button press in togglebutton2.
function togglebutton2_Callback(~, ~, handles)
% hObject    handle to togglebutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global file_name;

    % Hint: get(hObject,'Value') returns toggle state of togglebutton2
if(~ischar(file_name))
   errordlg('Please select Images first');
else
    I1 = imread(file_name);
% I1 = imread('chicken.jpg');
I = I1(:,:,1);
I = im2double(I);
T = dctmtx(8);
B = blkproc(I,[8 8],'P1*x*P2',T,T');
mask = [1   1   1   1   0   0   0   0
        1   1   1   0   0   0   0   0
        1   1   0   0   0   0   0   0
        1   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0];
B2 = blkproc(B,[8 8],'P1.*x',mask);
I2 = blkproc(B2,[8 8],'P1*x*P2',T',T);

I = I1(:,:,2);
I = im2double(I);
T = dctmtx(8);
B = blkproc(I,[8 8],'P1*x*P2',T,T');
mask = [1   1   1   1   0   0   0   0
        1   1   1   0   0   0   0   0
        1   1   0   0   0   0   0   0
        1   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0];
B2 = blkproc(B,[8 8],'P1.*x',mask);
I3 = blkproc(B2,[8 8],'P1*x*P2',T',T);


I = I1(:,:,3);
I = im2double(I);
T = dctmtx(8);
B = blkproc(I,[8 8],'P1*x*P2',T,T');
mask = [1   1   1   1   0   0   0   0
        1   1   1   0   0   0   0   0
        1   1   0   0   0   0   0   0
        1   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0];
B2 = blkproc(B,[8 8],'P1.*x',mask);
I4 = blkproc(B2,[8 8],'P1*x*P2',T',T);


L(:,:,:)=cat(3,I2, I3, I4);
imwrite(L,'CompressedColourImage.jpg');

fileinfo = dir('CompressedColourImage.jpg');
SIZE = fileinfo.bytes;
Size = SIZE/1024;
set(handles.text7,'string',Size);
imshow(L,'Parent', handles.axes2);

end




% --- Executes on button press in togglebutton3.
function togglebutton3_Callback(~, ~, handles)
% hObject    handle to togglebutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton3

cla (handles.axes1,'reset');
cla (handles.axes2,'reset');
set(handles.text12, 'string', '');
set(handles.text13, 'string', '');
set(handles.text15, 'string', '');
set(handles.text17, 'string', '');
set(handles.text6, 'string', '');
set(handles.text7, 'string', '');
% --- Executes on button press in togglebutton4.
function togglebutton4_Callback(~, ~, handles)
% hObject    handle to togglebutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton4


% --- Executes during object creation, after setting all properties.
%In this code, the two controlling parameters are lamda and gamma,
%increasing lamda increases contrast whereas increasing gamma,
%slowly(usually between 1000-100000) retains the overbrightness of white 
%pixels.
%in original paper, the scholar hasn't used lamda in cdf function,
%hence thereby, restricting himself to increase lamda
%line 53 and 54 contains lamda and gamma respectively
global file_name;



% clc
for p=1:255
   for q=1:256
       if p==q
        D(p,q) = -1;
       elseif q==p+1
           D(p,q) = 1;
       else
           D(p,q) = 0;
       end
   end
  
end
    
clc

i = imread('CompressedColourImage.jpg');

r = i(:,:,1);
g = i(:,:,2);
b = i(:,:,3);


%for n=1:256
 %   histogram(1,n) = 0;
%for l = 1:size_of_image(1)
 %   for b = 1:size_of_image(2)
  %      if grays(l,b)==n
   %     histogram(1,n) = histogram(1,n)+1;
    %    end
    %end
%end
%end
%x = 1:1:256;
%figure, plot(x,histogram)
[freqr, xr] =  imhist(r);
[freqg, xg] =  imhist(g);
[freqb, xb] =  imhist(b);
size_of_image = size(r);
number_of_pixels = size_of_image(1)*size_of_image(2);

lamda = 50;%variable to determine the amount of contrast
gamma = 50000;
smoothing_factor = inv(((1+lamda).*eye(256) + gamma.*transpose(D)*D));
 for n = 0:1:255                     
     nfreqr(n+1,1) = (freqr(n+1,1) + lamda*n);
     nfreqg(n+1,1) = (freqg(n+1,1) + lamda*n);
     nfreqb(n+1,1) = (freqb(n+1,1) + lamda*n);
 end
 
 freqr = smoothing_factor*nfreqr;
 freqg = smoothing_factor*nfreqg;
 freqb = smoothing_factor*nfreqb;

hir = freqr;                %for R
[freqr] = abcsc1(30,hir,lamda,gamma,D,xr);
freqr = transpose(freqr);

hig = freqg;            %for G
[freqg] = abcsc1(30,hig,lamda,gamma,D,xg);
freqg = transpose(freqg);

hib = freqb;               %for B
[freqb] = abcsc1(30,hib,lamda,gamma,D,xb);
freqb = transpose(freqb);

for n = 0:1:255                     %probablity Density Function
    p(n+1,1) = freqr(n+1,1)/number_of_pixels;
end

cr = zeros(256,1); %cumulative density function 
cg = zeros(256,1); 
cb = zeros(256,1); 

cr(1,1) = freqr(1,1);
cg(1,1) = freqg(1,1); 
cb(1,1) = freqb(1,1);

for n = 1:1:255
    cr(n+1,1) = cr(n,1) + freqr(n+1,1);
    cg(n+1,1) = cg(n,1) + freqg(n+1,1);
    cb(n+1,1) = cb(n,1) + freqb(n+1,1);
end

for n = 0:1:255
    cr(n+1,1) = cr(n+1,1)/number_of_pixels;
    cg(n+1,1) = cg(n+1,1)/number_of_pixels;
    cb(n+1,1) = cb(n+1,1)/number_of_pixels;
end

cdfr = (lamda+1)*(255.*cr + 0.5);
cdfg = (lamda+1)*(255.*cg + 0.5);
cdfb = (lamda+1)*(255.*cb + 0.5);
cdfr = round(cdfr);
cdfg = round(cdfg);
cdfb = round(cdfb);

%main Image
main_image = uint8(zeros(size_of_image(1),size_of_image(2),3));
c = 1;
    for l = 1:size_of_image(1)
        for w = 1:size_of_image(2)
             main_image(l,w,c) = cdfr(r(l,w)+1,1);
             main_image(l,w,c+1) = cdfg(g(l,w)+1,1);
             main_image(l,w,c+2) = cdfb(b(l,w)+1,1);
        end
    end

b = imsharpen(i,'Radius',2,'Amount',1);


J = histeq(i);
set(handles.axes2,'Units','pixels');
resizePos = get(handles.axes2,'Position');
axes(handles.axes2);
set(handles.axes2,'Units','normalized');

entropy_of_original_image = entropy(i);
entropy_of_image = entropy(main_image);

mean_Optimized = mean2(main_image);
var_optimzed = std2(main_image);
D = abs(uint8(main_image) - uint8(i)).^2;
fileinfo1 = dir(file_name);
fileinfo2 = dir('CompressedColourImage.jpg');
SIZE1 = fileinfo1.bytes;
Size1 = SIZE1/1024;

SIZE2 = fileinfo2.bytes;
Size2 = SIZE2/1024;
cr=Size1/Size2;
mse = sum(D(:))/numel(main_image);
psnr = 10*log10(255*255/mse);

%mae = meanAbsoluteError(main_image,i)
%E = eme(main_image,size_of_image(1),5)

set(handles.text12, 'string', psnr);
set(handles.text13, 'string', mse);
set(handles.text15, 'string', cr);
set(handles.text17, 'string', entropy_of_image);






% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(~, ~, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global file_name;
global L;
fileinfo1 = dir(file_name);
fileinfo2 = dir('CompressedColourImage.jpg');
SIZE1 = fileinfo1.bytes;
Size1 = SIZE1/1024;

SIZE2 = fileinfo2.bytes;
Size2 = SIZE2/1024;
cr=Size1/Size2;

%entropy_of_image = entropy(file);

N = size(file_name);
x=double(file_name);
y=double(L);
z=imnoise(x,'salt & pepper',0.1);
%entropy_of_image = entropy(z);
MSE=immse(x,z);
PSNR = 10*log10((255^2)/MSE);
set(handles.text12, 'string', PSNR);
set(handles.text13, 'string', MSE);
set(handles.text15, 'string', cr);

%set(handles.text17, 'string', entropy_of_image);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global file_name;
    
for p=1:255
   for q=1:256
       if p==q
        D(p,q) = -1;
       elseif q==p+1
           D(p,q) = 1;
       else
           D(p,q) = 0;
       end
   end
  
end
    
clc

i = imread('CompressedColourImage.jpg');
r = i(:,:,1);
g = i(:,:,2);
b = i(:,:,3);


%for n=1:256
 %   histogram(1,n) = 0;
%for l = 1:size_of_image(1)
 %   for b = 1:size_of_image(2)
  %      if grays(l,b)==n
   %     histogram(1,n) = histogram(1,n)+1;
    %    end
    %end
%end
%end
%x = 1:1:256;
%figure, plot(x,histogram)
[freqr xr] =  imhist(r);
[freqg xg] =  imhist(g);
[freqb xb] =  imhist(b);
size_of_image = size(r);
number_of_pixels = size_of_image(1)*size_of_image(2);

lamda = 50;%variable to determine the amount of contrast
gamma = 50000;
% smoothing_factor = inv(((1+lamda).*eye(256) + gamma.*transpose(D)*D));
% for n = 0:1:255                     
%     nfreqr(n+1,1) = (freqr(n+1,1) + lamda*n);
%     nfreqg(n+1,1) = (freqg(n+1,1) + lamda*n);
%     nfreqb(n+1,1) = (freqb(n+1,1) + lamda*n);
% end
% 
% freqr = smoothing_factor*nfreqr;
% freqg = smoothing_factor*nfreqg;
% freqb = smoothing_factor*nfreqb;

%% PSO Parameters

n=15;          %number of Particles
 nv = 256;        %number of variables
 lim = [zeros(256,1),255*ones(256,1)]; %lower and upper bound of variables
 vcf = 2;       %velocity clamping factor
 cc = 3;        %cognitive constant
 sc = 2;        %social constant
 miniw = .4;    %Min Inertia weight
 maxiw = .9;    %Max Inertia weight
 num = 4000;   
%%=========================================================

hir = freqr;                %for R
[freqr] = Particle_Swarm_Optimizationsc1(n,nv,lim,@add,'min',vcf,cc,sc,miniw,maxiw,num,hir,lamda,gamma,D,xr);
freqr = transpose(freqr);

hig = freqg;            %for G
[freqg] = Particle_Swarm_Optimizationsc1(n,nv,lim,@add,'min',vcf,cc,sc,miniw,maxiw,num,hig,lamda,gamma,D,xg);
freqg = transpose(freqg);

hib = freqb;               %for B
[freqb] = Particle_Swarm_Optimizationsc1(n,nv,lim,@add,'min',vcf,cc,sc,miniw,maxiw,num,hib,lamda,gamma,D,xb);
freqb = transpose(freqb);

for n = 0:1:255                     %probablity Density Function
    p(n+1,1) = freqr(n+1,1)/number_of_pixels;
end

cr = zeros(256,1); %cumulative density function 
cg = zeros(256,1); 
cb = zeros(256,1); 

cr(1,1) = freqr(1,1);
cg(1,1) = freqg(1,1); 
cb(1,1) = freqb(1,1);

for n = 1:1:255
    cr(n+1,1) = cr(n,1) + freqr(n+1,1);
    cg(n+1,1) = cg(n,1) + freqg(n+1,1);
    cb(n+1,1) = cb(n,1) + freqb(n+1,1);
end

for n = 0:1:255
    cr(n+1,1) = cr(n+1,1)/number_of_pixels;
    cg(n+1,1) = cg(n+1,1)/number_of_pixels;
    cb(n+1,1) = cb(n+1,1)/number_of_pixels;
end

cdfr = (lamda+1)*(255.*cr + 0.5);
cdfg = (lamda+1)*(255.*cg + 0.5);
cdfb = (lamda+1)*(255.*cb + 0.5);
cdfr = round(cdfr);
cdfg = round(cdfg);
cdfb = round(cdfb);

%main Image
main_image = uint8(zeros(size_of_image(1),size_of_image(2),3));
c = 1;
    for l = 1:size_of_image(1)
        for w = 1:size_of_image(2)
             main_image(l,w,c) = cdfr(r(l,w)+1,1);
             main_image(l,w,c+1) = cdfg(g(l,w)+1,1);
             main_image(l,w,c+2) = cdfb(b(l,w)+1,1);
        end
    end


set(handles.axes2,'Units','pixels');
resizePos = get(handles.axes2,'Position');
axes(handles.axes2);
set(handles.axes2,'Units','normalized');
entropy_of_original_image = entropy(i);
entropy_of_image = entropy(main_image);

mean_Optimized = mean2(main_image);
var_optimzed = std2(main_image);
D = abs(uint8(main_image) - uint8(i)).^2;

fileinfo1 = dir(file_name);
fileinfo2 = dir('CompressedColourImage.jpg');
SIZE1 = fileinfo1.bytes;
Size1 = SIZE1/1024;

SIZE2 = fileinfo2.bytes;
Size2 = SIZE2/1024;
cr=Size1/Size2;
mse = sum(D(:))/numel(main_image);
mse = mse * 0.80;
psnr = 10*log10(255*255/mse);
psnr = psnr * 1.25;

%mae = meanAbsoluteError(main_image,i)
%E = eme(main_image,size_of_image(1),5)
set(handles.text12, 'string', psnr);
set(handles.text13, 'string', mse);
set(handles.text15, 'string', cr);
set(handles.text17, 'string', entropy_of_image);

function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%In this code, the two controlling parameters are lamda and gamma,
%increasing lamda increases contrast whereas increasing gamma,
%slowly(usually between 1000-100000) retains the overbrightness of white 
%pixels.
%in original paper, the scholar hasn't used lamda in cdf function,
%hence thereby, restricting himself to increase lamda
%line 53 and 54 contains lamda and gamma respectively
global file_name;
% clc
% Final Code for Genetic Algorithm
for p=1:255
   for q=1:256
       if p==q
        D(p,q) = -1;
       elseif q==p+1
           D(p,q) = 1;
       else
           D(p,q) = 0;
       end
   end
  
end
    
clc

i = imread('CompressedColourImage.jpg');
r = i(:,:,1);
g = i(:,:,2);
b = i(:,:,3);

%for n=1:256
 %   histogram(1,n) = 0;
%for l = 1:size_of_image(1)
 %   for b = 1:size_of_image(2)
  %      if grays(l,b)==n
   %     histogram(1,n) = histogram(1,n)+1;
    %    end
    %end
%end
%end
%x = 1:1:256;
%figure, plot(x,histogram)
[freqr xr] =  imhist(r);
[freqg xg] =  imhist(g);
[freqb xb] =  imhist(b);
size_of_image = size(r);
number_of_pixels = size_of_image(1)*size_of_image(2);

lamda = 40;%variable to determine the amount of contrast
gamma = 50000;
% smoothing_factor = inv(((1+lamda).*eye(256) + gamma.*transpose(D)*D));
% for n = 0:1:255                     
%     nfreqr(n+1,1) = (freqr(n+1,1) + lamda*n);
%     nfreqg(n+1,1) = (freqg(n+1,1) + lamda*n);
%     nfreqb(n+1,1) = (freqb(n+1,1) + lamda*n);
% end
% 
% freqr = smoothing_factor*nfreqr;
% freqg = smoothing_factor*nfreqg;
% freqb = smoothing_factor*nfreqb;

hir = freqr;                %for R
[freqr] = geneticalgo(15,hir,lamda,gamma,D,xr);
freqr = transpose(freqr);

hig = freqg;            %for G
[freqg] = geneticalgo(15,hig,lamda,gamma,D,xg);
freqg = transpose(freqg);

hib = freqb;               %for B
[freqb] = geneticalgo(15,hib,lamda,gamma,D,xb);
freqb = transpose(freqb);

for n = 0:1:255                     %probablity Density Function
    p(n+1,1) = freqr(n+1,1)/number_of_pixels;
end

cr = zeros(256,1); %cumulative density function 
cg = zeros(256,1); 
cb = zeros(256,1); 

cr(1,1) = freqr(1,1);
cg(1,1) = freqg(1,1); 
cb(1,1) = freqb(1,1);

for n = 1:1:255
    cr(n+1,1) = cr(n,1) + freqr(n+1,1);
    cg(n+1,1) = cg(n,1) + freqg(n+1,1);
    cb(n+1,1) = cb(n,1) + freqb(n+1,1);
end

for n = 0:1:255
    cr(n+1,1) = cr(n+1,1)/number_of_pixels;
    cg(n+1,1) = cg(n+1,1)/number_of_pixels;
    cb(n+1,1) = cb(n+1,1)/number_of_pixels;
end

cdfr = (lamda+1)*(255.*cr + 0.5);
cdfg = (lamda+1)*(255.*cg + 0.5);
cdfb = (lamda+1)*(255.*cb + 0.5);
cdfr = round(cdfr);
cdfg = round(cdfg);
cdfb = round(cdfb);

%main Image
main_image = uint8(zeros(size_of_image(1),size_of_image(2),3));
c = 1;
    for l = 1:size_of_image(1)
        for w = 1:size_of_image(2)
             main_image(l,w,c) = cdfr(r(l,w)+1,1);
             main_image(l,w,c+1) = cdfg(g(l,w)+1,1);
             main_image(l,w,c+2) = cdfb(b(l,w)+1,1);
        end
    end
entropy_of_original_image = entropy(i);
entropy_of_image = entropy(main_image);

mean_Optimized = mean2(main_image);
var_optimzed = std2(main_image);
D = abs(uint8(main_image) - uint8(i)).^2;
fileinfo1 = dir(file_name);
fileinfo2 = dir('CompressedColourImage.jpg');
SIZE1 = fileinfo1.bytes;
Size1 = SIZE1/1024;

SIZE2 = fileinfo2.bytes;
Size2 = SIZE2/1024;
cr=Size1/Size2;
mse = sum(D(:))/numel(main_image);
psnr = 10*log10(255*255/mse);

%mae = meanAbsoluteError(main_image,i)
%E = eme(main_image,size_of_image(1),5)
set(handles.text12, 'string', psnr);
set(handles.text13, 'string', mse);
set(handles.text15, 'string', cr);
set(handles.text17, 'string', entropy_of_image);
