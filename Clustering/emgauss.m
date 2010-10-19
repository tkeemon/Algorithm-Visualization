
function [W,M,S,A] = emgauss(X,nk,niter)
% function [W,M,S,A] = emgauss(X,nk,niter)
%
% Implementation of expectation-maximization algorithm
%   for 2D gaussian distributions
%
% INPUTS
%   X - np x nd, number of nd-dimensional points
%   nk - (2), number of clusters
%   niter = (50), number of iterations
% OUTPUTS
%   W - np x nk, soft label for each point
%   M - nk x nd, mean
%   S - nd x nd x nk, variance
%   A - 1 x nk, mixture

% TJ Keemon <keemon@bc.edu> 
% 9 April 2009

if nargin<1 | isempty(X),
    help(mfilename);
    return;
end

if nargin<2 | isempty(nk),
    nk = 2;
end

if nargin<3 | isempty(niter),
    niter = 50;
end

[np,nd] = size(X);

% random initialization
L = [1:nk,zeros(1,np-nk)];
L = find(L(randperm(np)));
M = X(L,:);
S = repmat(eye(nd),[1,1,nk]);
A = repmat(1/nk,1,nk);
W = repmat(1/nk,np,nk);

for iter=1:niter,
    
    % visualization
    showptgauss(X,M,S,W); title(iter,'FontSize',20); pause(0.5);
    
    % E-step: update W
    P = zeros(np,nk);
    for i = 1:nk
        sig = S(:,:,i);
        mu = ones(size(X)) * diag(M(i,:));
        
        base = A(i)/((2*pi)^nd*det(sig))^.5; 
        vals = base*exp(-.5*(X-mu)*inv(sig)*(X-mu)');
        W(:,i) = diag(vals);
    end
    W = W ./ repmat(sum(W')',[1 nk]);
    
    % M-step: update M, S
    [in in] = max(W');
    for i = 1:nk
        Xs = X(in==i,:);
        Ws = W(in==i,i);
        
        M0(i,:) = sum(repmat(Ws,[1 nd]) .* Xs)/numel(Ws);
        S0(:,:,i) = cov(Xs);
        A(i) = numel(in==i)/numel(in);
    end
    % a very simple convergence test
    err = max(abs([M(:)-M0(:); S(:)-S0(:)]));
    M = M0
    S = S0;
    if err<1e-3,
        return;
    end
end
