
function [len,sol,seq] = dfs(maze,start,goal)
% function [len,sol,seq] = dfs(maze,start,goal)
%
% Maze search using the depth first search algorithm.
%
% INPUTS
%   maze - nr x nc x 4, navigability at every point in the maze along N,E,S,W directions
%   start - 1 x 1, starting position index into the maze
%   goal - 1 x 1, goal position index into the maze
% OUTPUTS
%   len - 1 x 1, length of the path connecting start to goal, inf if not navigable
%   sol - (len+1) x 1, solution path from start to goal by depth first search
%   seq - nr x nc, maze cell visit sequence id during the search, 0 if
%       never visited

% TJ Keemon <keemon@bc.edu>
% February 4, 2009

[h w d] = size(maze);

seq = zeros(h,w);
dir = ones(h*w,1);

stack = start;
seq(start) = -1;

move = 1;
while 1
    %if stack empty, break;
    if isempty(stack)
        break;
    end
    
    current = stack(end); %take current value
    
    [i j] = ind2sub([h w],current); %slowest call in my function
                                    %   but it's not worth rewriting
    
    switch dir(current)
        case 1
            D = maze(i,j,1);
            off = -1;
        case 2
            D = maze(i,j,2);
            off = h;
        case 3
            D = maze(i,j,3);
            off = 1;
        case 4
            D = maze(i,j,4);
            off = -h;
        case 5
            stack(end) = [];
    end

    if D == 1 && seq(current+off) == 0
        seq(current+off) = move;
        stack = [stack current+off];
        move = move + 1;
        if current+off == goal
            break;
        end
    else
        dir(current) = dir(current) + 1;
    end
end

if ~isempty(stack)  %if a solution was found
    sol = stack;
    len = numel(sol)-1;
else
    sol = [];
    len = 0;
end

if nargout == 0
    showmaze(maze,start,goal,seq,sol);
end