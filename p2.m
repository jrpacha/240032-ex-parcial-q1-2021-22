clearvars
close all

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P2
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART A
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the meshed domain meshHole.m from the course library of meshed
% domains, Count how many nodes in the mesh lie enterely above the line
% y = 1.9 x - 0.4. That is, the nodes that satisfy y > 1.9 x - 0.4

% Hint: The mean of the nodes strictly above the line is 
% (2.7694e-01, 7.0548e-01)

alpha = 1.9; 
beta = -0.4;

eval('meshHole');
numNod = size(nodes,1);
numElem = size(elem,1);

idx = find(nodes(:,2) > 1.9*nodes(:,1) - 0.4);
n = length(idx);
xMean = sum(nodes(idx,1))/n; yMean = sum(nodes(idx,2))/n;

fprintf('\nPART A:\n')
fprintf('Number of nodes enterely above y = %.4f x %+.4f, n = %d\n',...
    alpha,beta,n)
fprintf('* Hint. Mean of the nodes strictly above the line: (%.4e,%.4e)\n',...
    xMean,yMean)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part B
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We have measured the following values (x_{1}, y_{1}),...,(x_{6}, y_{6})
% for a function y(x)
%
% x = [ 1, 3, 5, 7, 9, 12];
% y = [-1.1, 0.4, 0.5, 1.2, 2.1, 4.4];
%
% We wish to approximate y(x) with an approximation or regression
% polynomial p(x) on the interval [x_1, x_6] of the measured values. To
% decide which degree to choose, we will compute the error of each degree
% d (ranging from 1 to 5) as:
%
% E(d) = (d+1)(Delta_y)^2 + sum_{i = 1}^{6} (y(x_{i})-p(x_{i}))^2
%
% where the Delta_y term depends on the data and is the mean of the
% differences |y_{i+1} - y_{i}|, for i=1,...,5

x = [ 1, 3, 5, 7, 9, 12];
y = [-1.1, 0.4, 0.5, 1.2, 2.1, 4.4];

numNod = size(x,2);

DeltaY = sum(abs(y(2:end)-y(1:end-1)))/(numNod-1);
E = zeros(1,numNod);

for d = 1:numNod
    p = polyfit(x,y,d-1);
    px = polyval(p,x);
    E(d) = sum((y-px).^2) + d*DeltaY^2;
end

[minE,indexMinE]=min(E);

fprintf('\nPART B:\n')
fprintf('Value of Delta_Y = %.4e\n', DeltaY)
fprintf('Degree d that makes minumum E(d), d = %d\n',indexMinE-1)
fprintf('* Hint. For degree d = 5, E(5) = %.4e\n',E(6))