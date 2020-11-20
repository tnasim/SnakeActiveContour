clear
clc
global currentSnake
global affectedIndices
global updateSnake
global mouseX
global mouseY
global nearestIdx
global isDragging
global line
global done

RECORD_VIDEO = 0
% img1 = imread('images/coins.tif');
img1 = imread('k24s-raw.gif');

if size(size(img1), 2) == 3
    img = rgb2gray(img1);
else
    img = img1;
end

figure()
imshow(img)
axis on

% Mouse info
x_cord = [];
y_cord = [];
my_vertices = [];

button = 1;
count = 1;
hold on

h = plot( [0 0], [0 0] , 'g');
 while sum(button) <=1   % read ginputs until a mouse right-button occurs
   [x_temp,y_temp,button] = ginput(1);
   x_cord(count) = x_temp;
   y_cord(count) = y_temp;
   count = count+1;

   my_vertices = [x_cord; y_cord]';
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

% Uncomment below if you want to see the outline of the points for
% initalization
% my_vertices = [x_cord; y_cord]';
% h = drawpolygon('Position',my_vertices);


%%
% axis([0 500 0 500])
%[x, y] = ellipse(120, 15, 150, 150, .1);

x = xInt';
y = yInt';

hold on
disp(size(img));
% im = plot(img);
% h = plot(x, y, 'b-o','MarkerIndices',1:5:length(y));
h = plot(x, y, 'b-');

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

% parameters
g = 5;
N = length(x);

% creating tri-diagonal branded matrix:
r = [2 -1 zeros(1,N-2)];
alpha = toeplitz(r); % <-- creates a diagonal branded matrix
% update the corner values
alpha(1, 1) =  2;
alpha(1, N) = -1;
alpha(N, 1) = -1;
alpha(N, N) =  2;

% creating penta-diagonal branded matrix:
r2 = [6 -4 1 zeros(1,N-3)];
beta = toeplitz(r2); % <-- creates a diagonal branded matrix
% update the corner values
beta(1, 1) =  6;
beta(1, N) = -4;
beta(N, 1) = -4;
beta(N, N) =  6;
beta(1, N-1) = 1;
beta(2, N) = 1;
beta(N-1, 1) = 1;
beta(N, 2) = 1;
disp(alpha);
disp(beta);
A = alpha + beta;

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


% Adding Images forces
lineFunction = 1;
edgeFunction = 1;

% Image gradient for external energy
lineForce = double(img);
lineForce = lineForce / max(lineForce(:));

[magnitude, direction] = imgradient(img);
magnitude = magnitude / max(magnitude(:));

% total image force
Ext = lineFunction * lineForce - edgeFunction * magnitude;
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
   
   % Instead of taking inverse of 'first_term' using A\b format which is
   % faster than INV(A)*b (suggested by matlab docs)   
   x = first_term\(g*x - (fx + distance(:,1)));
   y = first_term\(g*y - (fy + distance(:,2)));
   
   % fill up the gap between last and first point.
   x(1) = ( x(1) + x(100) )/2.0 - 0.1;
   y(1) = ( y(1) + y(100) )/2.0 - 0.1;
   x(100) = ( x(1) + x(100) )/2.0 + 0.1;
   y(100) = ( y(1) + y(100) )/2.0 + 0.1;
   
   currentSnake = [x(:) y(:)];
   updateSnake = affectedIndices .* updateSnake + (1 - affectedIndices) .* currentSnake;
   
   set(h, 'YData', y);
   set(h, 'XData', x);
   
   if done
       set(h, 'Color', 'r');
%        set(h, 'LineStyle', '-');
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

%% Handles mouse click event
function onClick (object, eventdata)
    global currentSnake
    global affectedIndices
    global updateSnake

    click_point = get (gca, 'CurrentPoint');
    button = get(object,'SelectionType');
    click_x = click_point(1, 1); click_y = click_point(1, 2);
    p_new = [click_x click_y];
    distances = sqrt(sum(( currentSnake - p_new ).^2,2));
    row = find(distances==min(distances));     

    % Capture both left and right mouse clicks
    if strcmpi(button,'normal') | strcmpi(button,'alt')
        updateSnake(row, 1) = click_x;
        updateSnake(row, 2) = click_y;
        affectedIndices = zeros(size( currentSnake ));
        affectedIndices(row, :) = 1;
    end
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
        if strcmpi(button,'normal') | strcmpi(button,'alt')
            updateSnake(row, 1) = click_x;
            updateSnake(row, 2) = click_y;
            affectedIndices = zeros(size( currentSnake ));
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
    global done

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
    
    if strcmpi(button,'alt')
        done = 1;
    end
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

%% Returns an ellipse based on the input parameters.
function [x, y] = ellipse(x1, y1, x2, y2, e)
	a = 1/2*sqrt((x2-x1)^2+(y2-y1)^2);
	b = a*sqrt(1-e^2);
    t = linspace(0,2*pi);
    X = a*cos(t);
    Y = b*sin(t);
    w = atan2(y2-y1,x2-x1);
    x = (x1+x2)/2 + X*cos(w) - Y*sin(w);
    y = (y1+y2)/2 + X*sin(w) + Y*cos(w);
end