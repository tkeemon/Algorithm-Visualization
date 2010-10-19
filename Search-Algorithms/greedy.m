
function [len,sol,seq] = greedy(maze,start,goal)
% function [len,sol,seq] = greedy(maze,start,goal)
%
% Maze search using a greedy algorithm with a city-block metric.
%
% INPUTS
%   maze - nr x nc x 4, navigability at every point in the maze along N,E,S,W directions
%   start - 1 x 1, starting position index into the maze
%   goal - 1 x 1, goal position index into the maze
% OUTPUTS
%   len - 1 x 1, length of the path connecting start to goal, inf if not navigable
%   sol - (len+1) x 1, solution path from start to goal by greedy search on city-block distance
%   seq - nr x nc, maze cell visit sequence id during the search, 0 if
%       never visited

% TJ Keemon <keemon@bc.edu>
% February 23, 2009

if nargin < 1
    help greedy
    return;
end

[h w d] = size(maze);

dmap = zeros(h,w); dmap(goal) = 1;
dmap = bwdist(dmap,'cityblock');
seq = zeros(h,w); seq(start) = -1;

dVals = [-1 0; 0 1; 1 0; 0 -1] * [1;h];
dInd = (0:d-1)'*(h*w);

index = 0;
stack = [start; zeros(h*w,1)]; ptr = 1;
while ptr > 0
    
    current = stack(ptr); %take current value
    
    %if soluton found, break
    if current == goal
        break;
    end
    
    t = find(maze(current+dInd));
    t = t(seq(current+dVals(t))==0);
    
    if ~isempty(t)        
        [v i] = min(dmap(current+dVals(t)));
        
        current = current + dVals(t(i));
        ptr = ptr + 1;
        stack(ptr) = current;
        
        index = index + 1;
        seq(current) = index;
    else
        ptr = ptr - 1;
    end
%     showmaze(maze,start,goal,seq);
end

if ptr ~= 0
    len = ptr - 1;
    sol = stack(1:ptr);
else
    len = inf;
    sol = [];
end

if nargout == 0
    showmaze(maze,start,goal,seq,sol);
end