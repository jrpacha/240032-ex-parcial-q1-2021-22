clearvars
close all

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P3
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider the boundary value problem (BVP) in a 1-dim domain [0, L] with 
% L = 1
% 
% -((x+4)u')' + u = x,       0 < x < L,
%           u'(0) = u'_{0},
%           u'(L) = u'_{L}.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART A
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fix u'_{0} = 2.0, u'_{L} = 3.5 in the BVP above, and consider its FEM
% solution with N = 50 linear elements (notice that no fixed nodes are
% present). If {u_{i}}_{i = 1,...,N+1} is the corresponding nodal solution,
% then max_{i=1,...,N+1} | u_{i} |. Hint: The solution for node 42 is 
% u_{42} = 1.09034e-01
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART B
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The interpolated value of u at x = 0.77 is
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART C
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now, consider the approximation of u using a cubic spline built from the
% nodal solution of the BVP computed in part A. The absolute value of the
% difference between the interpolated u at point x = 0.77 found in part B
% and the one obtained using the spline is:

L = 1.0;     %length of the interval
N = 50;      %number of divisions

dudx0 = 2.0;  %du/dx at x = 0 
dudxL = 3.5; %du/dx at x = L

p = 0.77;    %point to compute u

alpha = 1;
beta = 4.0;  % a1(x) = alpha*x + beta
gamma = 1.0; % a0(x) = gamma

a1 = @(x) alpha*x + beta;

h = L/N;
nodes = (0:h:L)';
elem = [(1:N)',(2:N+1)'];

numNod = size(nodes,1);
numElem = size(elem,1);

%Assembly
K = zeros(numNod);
F = zeros(numNod,1);
Q = zeros(numNod,1);

for e = 1:numElem
    rows = [elem(e,1); elem(e,2)];
    cols = rows;
    x1 = nodes(rows(1)); %Position of 1st node of element e
    x2 = nodes(rows(2)); %Position of 2nd node of element e
    Ke = 0.5*alpha*(x1+x2)/h*[1,-1;-1,1] + ...
        beta/h*[1,-1;-1,1] + gamma*h/6*[2,1;1,2]; %local stiffness matrix
    Fe = h*[2*x1+x2; x1+2*x2]/6;                %local load vector                
    K(rows,cols) = K(rows,cols) + Ke;
    F(rows) = F(rows) + Fe;
end

%Boundary conditons
%Natural B.C.:
Q(1) = -a1(0)*dudx0;
Q(numNod) = a1(L)*dudxL;
%Essential B.C.: none

F = Q + F;
u = K\F;

maxU = norm(u,inf);
uP = interp1(nodes,u,p);
upSpline = spline(nodes,u,p);
err = abs(uP-upSpline);

% Output 
fprintf('(a) max |u| = %.5e\n',maxU)
fprintf('    * Hint u(42) = %.5e\n',u(42))
fprintf('(b) u_interp(0.77) = %.5e\n',uP)
fprintf('(c) |u_interp(0.77) - u_spline(0.77)| = %.5e\n',err)
