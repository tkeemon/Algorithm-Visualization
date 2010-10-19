
function [len,sol,seq] = bfs(maze,start,goal)
% function [len,sol,seq] = bfs(maze,start,goal)
%
% Maze search using the breadth first search algorithm.
%
% INPUTS
%   maze - nr x nc x 4, navigability at every point in the maze along N,E,S,W directions
%   start - 1 x 1, starting position index into the maze
%   goal - 1 x 1, goal position index into the maze
% OUTPUTS
%   len - 1 x 1, length of the path connecting start to goal, inf if not navigable
%   sol - (len+1) x 1, solution path from start to goal by breadth first search
%   seq - nr x nc, maze cell visit sequence id during the search, 0 if
%       never visited

% TJ Keemon <keemon@bc.edu>
% February 12, 2009

[h w d] = size(maze);

seq = zeros(h,w);

queue = start;
tree(1) = start;
seq(start) = -1;
parents(1) = 0;

index = 1;
move = 1;
while 1
    %if queue empty, break;
    if isempty(queue)
        break;
    end
    
    current = queue(1); %take current value
    queue(1) = [];      %pop frontmost value from queue
    
    [i j] = ind2sub([h w],current); %slowest call in my function
                                    %   but it's not worth rewriting
    D = maze(i,j,:);

    %search N
    if D(1) == 1 && seq(current-1) == 0
        seq(current-1) = move;
        queue = [queue current-1];
        move = move + 1;
        tree(index+1) = current-1;
        parents(index+1) = current;
        if current-1 == goal
            break;
        end
    end
    %search E
    if D(2) == 1 && seq(current+h) == 0
        seq(current+h) = move;
        queue = [queue current+h];
        move = move + 1;
        tree(index+2) = current+h;
        parents(index+2) = current;
        if current+h == goal
            break;
        end
    end
    %search S
    if D(3) == 1 && seq(current+1) == 0
        seq(current+1) = move;
        queue = [queue current+1];
        move = move + 1;
        tree(index+3) = current+1;
        parents(index+3) = current;
        if current+1 == goal
            break;
        end
    end
    %search W
    if D(4) == 1 && seq(current-h) == 0
        seq(current-h) = move;
        queue = [queue current-h];
        move = move + 1;
        tree(index+4) = current-h;
        parents(index+4) = current;
        if current-h == goal
            break;
        end
    end
    %showmaze(maze,start,goal,seq);
    index = index + 4;
end

if ~isempty(queue)  %if a solution was found
    current = goal;
    sol = current;
    while 1
        I = find(tree == current);

        current = parents(I);

        sol = [current sol];

        if current == start
            break;
        end
    end

    len = numel(sol)-1;
else
    sol = [];
    len = 0;
end

if nargout == 0
    showmaze(maze,start,goal,seq,sol);
end