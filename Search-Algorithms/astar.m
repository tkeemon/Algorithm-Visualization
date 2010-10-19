function [len,sol,seq] = astar(maze,start,goal)
% function [len,sol,seq] = astar(maze,start,goal)
%
% Maze search using the a* algorithm.
%
% INPUTS
%   maze - nr x nc x 4, navigability at every point in the maze along N,E,S,W directions
%   start - 1 x 1, starting position index into the maze
%   goal - 1 x 1, goal position index into the maze
% OUTPUTS
%   len - 1 x 1, length of the path connecting start to goal, inf if not navigable
%   sol - (len+1) x 1, solution path from start to goal by A* search on city-block distance
%   seq - nr x nc, maze cell visit sequence id during the search, 0 if
%       never visited

% TJ Keemon <keemon@bc.edu>
% February 23, 2009

[h w d] = size(maze);

dmap = zeros(h,w); dmap(goal) = 1;
dmap = bwdist(dmap,'cityblock');
gn = inf(h,w); gn(start) = -1;

seq = zeros(h,w); seq(start) = -1;
parents = zeros(h,w); 

dVals = [-1 0; 0 1; 1 0; 0 -1] * [1;h];
dInd = (0:d-1)'*(h*w);

if start == goal
    len = 0;
    sol = [start; goal];
    return
end

fringe = [start dmap(start) 0; zeros(h*w,1) inf(h*w,1) inf(h*w,1)];

index = 1;
while 1
    
    minDist = min(fringe(:,2));
    
    ind = find(fringe(:,2)==min(fringe(:,2)),1);
    current = fringe(ind,1);  %takes min value
    step = 1+fringe(ind,3);
    fringe(ind,2) = inf;
    
    %if no solution found
    if minDist == inf
        break;
    %if soluton found, break
    elseif current == goal && min(fringe(:,2)) > gn(goal)  %%fix min > goal
        break;
    end
    
    t = find(maze(current+dInd));         %finds all available moves
    nodes = current+dVals(t);
    
    %update nodes if closer path found
    update = find(dmap(nodes)+step < gn(nodes));
    gn(nodes) = step;
    parents(nodes(update)) = current;
    
    t = t(seq(current+dVals(t))==0);
    nodes = current+dVals(t);
    
    np = numel(t);
    seq(nodes) = index:index+np-1;
    
    fringe(index:index+np-1,:) = [nodes dmap(nodes)+gn(nodes) repmat(step,np,1)];
    index = index + np;    
end

if minDist ~= inf
    sol = goal;
    len = 1;
    cur = goal;
    
    while cur ~= start
        len = len+1;
        sol = [parents(cur); sol];
        cur = parents(cur);
    end
    len = len - 1;
else
    len = inf;
    sol = [];
end

if nargout == 0
    showmaze(maze,start,goal,seq,sol);
end