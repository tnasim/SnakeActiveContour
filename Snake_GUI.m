function varargout = Snake_GUI(varargin)
% SNAKE_GUI MATLAB code for Snake_GUI.fig
%      SNAKE_GUI, by itself, creates a new SNAKE_GUI or raises the existing
%      singleton*.
%
%      H = SNAKE_GUI returns the handle to a new SNAKE_GUI or the handle to
%      the existing singleton*.
%
%      SNAKE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SNAKE_GUI.M with the given input arguments.
%
%      SNAKE_GUI('Property','Value',...) creates a new SNAKE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Snake_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Snake_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Snake_GUI

% Last Modified by GUIDE v2.5 03-Dec-2020 19:06:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Snake_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Snake_GUI_OutputFcn, ...
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
end 

% --- Executes just before Snake_GUI is made visible.
function Snake_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Snake_GUI (see VARARGIN)

% Choose default command line output for Snake_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
end

% UIWAIT makes Snake_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Snake_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Executes on button press in Image.
function Image_Callback(hObject, eventdata, handles)
% hObject    handle to Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [filename, pathname] = uigetfile(...
        {'*.*'},'File Selector');
    handles.filename = strcat(pathname,'\',filename);
    guidata(hObject, handles);
    handles.filename;

    axes(handles.axes1);
    x = imread(handles.filename);
    handles.image = x;
    img1 = handles.image;
    if size(size(img1), 2) == 3
        img = rgb2gray(img1);
    else
        img = img1;
    end
    imshow(img)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize values on GUI after image is loaded. %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a = 1; % alpha
    b = 1; % beta
    g = 2; % gamma
    
    lineFunction = -1;
    edgeFunction = 0;
    terminalFunction  = 0;
    
    set(handles.Alpha,'string',num2str(a));
    set(handles.Beta,'string',num2str(b));
    set(handles.Gamma,'string',num2str(g));
    
    set(handles.linefunction,'string',num2str(lineFunction));
    set(handles.EdgeFunction,'string',num2str(edgeFunction));
    set(handles.terminalFunction,'string',num2str(terminalFunction));
    guidata(hObject, handles); 
end

function Alpha_Callback(hObject, eventdata, handles)
% hObject    handle to Alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% Hints: get(hObject,'String') returns contents of Alpha as text
%        str2double(get(hObject,'String')) returns contents of Alpha as a double


% --- Executes during object creation, after setting all properties.
function Alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function Beta_Callback(hObject, eventdata, handles)
% hObject    handle to Beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% Hints: get(hObject,'String') returns contents of Beta as text
%        str2double(get(hObject,'String')) returns contents of Beta as a double


% --- Executes during object creation, after setting all properties.
function Beta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function Gamma_Callback(hObject, eventdata, handles)
% hObject    handle to Gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end
% Hints: get(hObject,'String') returns contents of Gamma as text
%        str2double(get(hObject,'String')) returns contents of Gamma as a double


% --- Executes during object creation, after setting all properties.
function Gamma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function linefunction_Callback(hObject, eventdata, handles)
% hObject    handle to linefunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of linefunction as text
%        str2double(get(hObject,'String')) returns contents of linefunction as a double
end

% --- Executes during object creation, after setting all properties.
function linefunction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to linefunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function EdgeFunction_Callback(hObject, eventdata, handles)
% hObject    handle to EdgeFunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EdgeFunction as text
%        str2double(get(hObject,'String')) returns contents of EdgeFunction as a double
end

% --- Executes during object creation, after setting all properties.
function EdgeFunction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgeFunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');


end
end

function terminalFunction_Callback(hObject, eventdata, handles)
% hObject    handle to terminalFunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of terminalFunction as text
%        str2double(get(hObject,'String')) returns contents of terminalFunction as a double

end
% --- Executes during object creation, after setting all properties.
function terminalFunction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to terminalFunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');

end
end
% --- Executes on button press in Run.

function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% obtaining various parameters from the text boxes

global currentSnake
global affectedIndices
global updateSnake
global fixedPoint
global mouseX
global mouseY
global nearestIdx
global isDragging
global line
global done

OPEN_SNAKE = 1;
RECORD_VIDEO = 0;

img1 = handles.image;
if size(size(img1), 2) == 3
    img = rgb2gray(img1);
else
    img = img1;
end

imshow(img)
hold on;

% Control Alpha, Beta and Gamma.
a = str2double(get(handles.Alpha,'String'));
b = str2double(get(handles.Beta,'String'));
g = str2double(get(handles.Gamma,'String'));

% Control Images forces
lineFunction = str2double(get(handles.linefunction,'String'));
edgeFunction = str2double(get(handles.EdgeFunction,'String'));
terminalFunction  = str2double(get(handles.terminalFunction,'String'));

% housekeeping
set(handles.Image,'Enable','off');
set(handles.Run,'Enable','off');

% Mouse info
x_cord = [];
y_cord = [];

button = 1;
count = 1;
hold on

h = plot( [0 0], [0 0] , 'r');
 while sum(button) <=1   % read ginputs until a mouse right-button occurs
   [x_temp,y_temp,button] = ginput(1);
   x_cord(count) = x_temp;
   y_cord(count) = y_temp;
   count = count+1;
   
   set(h, 'XData', x_cord);
   set(h, 'YData', y_cord);
 end
 
 set(h, 'XData', [0 0]);
 set(h, 'YData', [0 0]);
 
%% X&Y points to ellipse
% Reparameterize and interp.
dist = sqrt(diff(x_cord).^2+diff(y_cord).^2);

% Reparameterize
tm = cumsum([0 dist]);
t = linspace(0,tm(end));
xInt = interp1(tm,x_cord,t,"spline");
yInt = interp1(tm,y_cord,t,"spline");
hold on

x = xInt';
y = yInt';

hold on

disp(size(img));
h = plot(x, y, 'g-');


nearestIdx = 1;
mouseX = x(nearestIdx);
mouseY = y(nearestIdx);
line = plot( [0 0], [0 0] , 'r');

currentSnake = [x(:) y(:)];
set (gcf, 'Pointer', 'crosshair');
% set (gcf, 'WindowButtonDownFcn', @onClick);
set (gcf, 'WindowButtonMotionFcn', @move);
set (gcf, 'WindowButtonUpFcn',@drop);
set (gcf, 'WindowButtonDownFcn', @drag);
set (gcf, 'WindowKeyPressFcn', @keyPress);

N = length(x);

% creating tri-diagonal branded matrix:
r = [2 -1 zeros(1,N-2)];
alpha = toeplitz(r); % <-- creates a diagonal branded matrix

% creating penta-diagonal branded matrix:
r2 = [6 -4 1 zeros(1,N-3)];
beta = toeplitz(r2); % <-- creates a diagonal branded matrix

if OPEN_SNAKE == 1
    % update the corner values
    alpha(1, 1) =  1;
    alpha(N, N) =  1;
    
    % update the corner values
    beta(1, 1) =  1;
    beta(N, N) =  1;
    beta(2, 1) = -2;
    beta(1, 2) = -2;
    beta(2, 2) = 5;
    beta(N, N-1) = -2;
    beta(N-1, N) = -2;
    beta(N-1, N-1) = 5;
else
    % update the corner values
    alpha(1, 1) =  2;
    alpha(1, N) = -1;
    alpha(N, 1) = -1;
    alpha(N, N) =  2;
    
    % update the corner values
    beta(1, 1) =  6;
    beta(1, N) = -4;
    beta(N, 1) = -4;
    beta(N, N) =  6;
    beta(1, N-1) = 1;
    beta(2, N) = 1;
    beta(N-1, 1) = 1;
    beta(N, 2) = 1;
    
end

A = a*alpha + b*beta;

% Calculate the first term from equation (19) and (20)
% first_term=inv(A + g.* eye(N));
first_term=(A + g.* eye(N));

% hard constraint init
updateSnake = currentSnake;
affectedIndices = zeros(size(currentSnake));
if RECORD_VIDEO == 1
    writerObj = VideoWriter('snake.avi');
    writerObj.FrameRate = 25;
    open(writerObj);
end

wConHard = 1;

% Image gradient for external energy
lineForce = double(img);
lineForce = lineForce / max(lineForce(:));

[magnitude, direction] = imgradient(img);
magnitude = magnitude / max(magnitude(:));

fixedPoint = currentSnake;
affectedIndices = zeros(size( currentSnake ));

%terminal energy calculation
I=double(img);

Ix=differentiation(I,0.1,'x');  
Iy=differentiation(I,0.1,'y');
Ixx=differentiation(I,0.1,'xx');
Ixy=differentiation(I,0.1,'xy');
Iyy=differentiation(I,0.1,'yy');

terminalEnergy = (Iyy.*Ix.^2 -2*Ixy.*Ix.*Iy + Ixx.*Iy.^2)./((1+Ix.^2 + Iy.^2).^(3/2));

% total image force
Ext = lineFunction * lineForce - edgeFunction * magnitude - terminalFunction * terminalEnergy;
[Fx,Fy] = imgradientxy(Ext);

isDragging = 0;
i = 1;
done = 0;

while 1
   x = currentSnake(:, 1);
   y = currentSnake(:, 2);
   
   fx = interp2(Fx,x,y);
   fy = interp2(Fy,x,y);
   
   % difference among current and updated snake points.
   distance = currentSnake - updateSnake;
   % hard constraint
   constraint = distance + wConHard * (currentSnake - fixedPoint);
   
   % Instead of taking inverse of 'first_term' using A\b format which is
   % faster than INV(A)*b (suggested by matlab docs)   
   x = first_term\(g*x - (fx + constraint(:,1)));
   y = first_term\(g*y - (fy + constraint(:,2)));
   
   if OPEN_SNAKE == 0
       % fill up the gap between last and first point.
       x(1) = ( x(1) + x(100) )/2.0 - 0.1;
       y(1) = ( y(1) + y(100) )/2.0 - 0.1;
       x(100) = ( x(1) + x(100) )/2.0 + 0.1;
       y(100) = ( y(1) + y(100) )/2.0 + 0.1;
   end
   
   currentSnake = [x(:) y(:)];
   updateSnake = currentSnake;
   
   % Use the following for Hard Constraint:
   fixedPoint = affectedIndices .* fixedPoint + (1 - affectedIndices) .* currentSnake;
   
   set(h, 'YData', y);
   set(h, 'XData', x);
   
   if done
       set(h, 'Color', 'r');
       break;
   end
   
%    line = plot( [x(nearestIdx) mouseX], [y(nearestIdx) mouseY] );

   if isDragging == 1
       set(line, 'YData', [y(nearestIdx) mouseY]);
       set(line, 'XData', [x(nearestIdx) mouseX]);
   else
       set(line, 'YData', [0 0]);
       set(line, 'XData', [0 0]);
   end
   
   if RECORD_VIDEO == 1
       writeVideo(writerObj,getframe(gcf));
       i = i + 1;
       disp(i);
       if i > 1500
           break;
       end
   end
   
   drawnow;
end

hold off

% -------------------------------------------

% close the writer object
if RECORD_VIDEO == 1
    close(writerObj);
end
% -------------------------------------------


end


%% Handles mouse click event
function move (object, eventdata)
    global currentSnake
    global affectedIndices
    global updateSnake
    global isDragging
    global mouseX
    global mouseY
    global nearestIdx
    global fixedPoint

    if isDragging == 1
        click_point = get (gca, 'CurrentPoint');
        button = get(object,'SelectionType');
        click_x = click_point(1, 1); click_y = click_point(1, 2);
        p_new = [click_x click_y];
        mouseX = click_x;
        mouseY = click_y;
%         distances = sqrt(sum(( currentSnake - p_new ).^2,2));
        row = nearestIdx
%         nearestIdx = row;

        % Capture both left and right mouse clicks
        if strcmpi(button,'normal')
            updateSnake(row, 1) = click_x;
            updateSnake(row, 2) = click_y;
        elseif strcmpi(button,'alt')
            fixedPoint(row, 1) = click_x;
            fixedPoint(row, 2) = click_y;
            affectedIndices(row, :) = 1;
        end
    end
end

%% Handles mouse click event
function drag (object, eventdata)
    global currentSnake
    global isDragging
    global mouseX
    global mouseY
    global nearestIdx

    click_point = get (gca, 'CurrentPoint');
    button = get(object,'SelectionType');
    click_x = click_point(1, 1); click_y = click_point(1, 2);
    p_new = [click_x click_y];
    mouseX = click_x;
    mouseY = click_y;
    distances = sqrt(sum(( currentSnake - p_new ).^2,2));
    row = find(distances==min(distances));
    nearestIdx = row;
        
    isDragging = 1;
end

%% Handles mouse click event
function drop (object, eventdata)
    global currentSnake
    global isDragging
    global mouseX
    global mouseY
    global nearestIdx
    global line

    if isDragging == 1
        isDragging = 0;
        nearestIdx = -1;
%         mouseX = currentSnake(nearestIdx, 1);
%         mouseY = currentSnake(nearestIdx, 2);
%         line = plot( [mouseX mouseX], [mouseY mouseY] );
    end
end

%% Handles keyboard entry
function keyPress (object, eventdata)
    global done
    
    % determine the key that was pressed
    keyPressed = eventdata.Key;
    
    % if spacebar is pressed, mark it done.
    if strcmpi(keyPressed,'space')
        done = 1;
    end
end



