
function showmaze(maze,start,goal,seq,sol)
% function showmaze(maze,start,goal,seq,sol)
%
% Visualization tool for maze search algorithms.
%
% INPUTS
%   maze - nr x nc x 4, navigability at every point in the maze along N,E,S,W directions
%   start - 1 x 1, starting position index into the maze
%   goal - 1 x 1, goal position index into the maze
%   seq - len x 1, visiting sequence from start to goal by depth first search
%   sol - solution path
% OUTPUTS
%   display the maze and the visiting sequence

% TJ Keemon <keemon@bc.edu>
% January 27, 2009

close all;
[h w d] = size(maze);

if nargin < 4
    seq = [];
end
if nargin < 3 || isempty(goal)
    goal = h*w;
end
if nargin < 2 || isempty(start)
    start = 1;
end
    

figure; hold on; 

for x = 0:w,
    plot([x x],[h 0],'k','LineWidth',1);
end
for x = 0:h,
    plot([0 w],[x x],'k','LineWidth',1);
end

for i = 1:4
    switch i
        case 1
            off = [0 .005]; dir = [1 0];
        case 2
            off = [.995 0]; dir = [0 1];
        case 3
            off = [0 .995]; dir = [1 0];
        case 4
            off = [.005 0]; dir = [0 1];
    end
    in = find(maze(:,:,i) == 0);
    [x y] = getIndex(in,h);
    
    %plots based on the direction and offset given in the above
    %   switch statement
    for val = 1:numel(x),
        xV = x(val) + off(1);
        yV = y(val) + off(2);
        plot([xV xV+dir(1)],[yV yV+dir(2)],'k-','LineWidth',2);
    end
end


%label start and end
[xs ys] = getIndex(start,h,[.5 .5]);
[xg yg] = getIndex(goal,h,[.5 .5]);
text(xs,ys,'S','FontSize',12,'Color',[1 0 0],'HorizontalAlignment','Center');
text(xg,yg,'G','FontSize',12,'Color',[1 0 0],'HorizontalAlignment','Center');

%plots the sequence of moves
if ~isempty(seq)
    seq = seq(:);
    [Y I] = sort(seq);
    numZeros = sum(seq == 0);
    I(1:numZeros+1) = [];  %removes untraveled cells

    for i = 1:numel(I)
        [x y] = getIndex(I(i),h,[.5 .5]);
        if I(i) ~= goal
            text(x,y,num2str(i),'FontSize',12,'Color',[0 0 0],'HorizontalAlignment','Center');
        end
    end
end

if ~isempty(sol)
    [x y] = getIndex(sol,h,[.5 .5]);
%     x = [x; x(1)]; y = [y; y(1)];
    len = numel(x);
    for i = 1:len-1
        plot([x(i) x(i+1)],[y(i) y(i+1)],'r-','LineWidth',2);
    end
end

axis ij;
axis off;
axis equal;
axis([-.1 w+.1 -.1 h+.1]);

%%
%I hate using subfunctions, but this one didn't warrent its own m file
function [x y] = getIndex(index,h,offset)
% getIndex - returns x,y coordinates from the given linear index
%           coordinates are of the bottom left corner for a given cell
% INPUTS
%   index - the linear index (can be a vector)
%   h - the height of the maze
%   offset - scaler indicating a desired offset for each coordinate
%            [0 0] - default
% OUTPUTS
%   [x y] - x and y coordinates for the plot (each offset by +.5)
if nargin < 3 
    offset = [0 0];
end
x = floor((index-1)./h) + offset(1);
y = mod(index-1,h)  + offset(2);
