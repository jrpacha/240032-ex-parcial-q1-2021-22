clearvars 
close all

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A certain 1D process is simulated by the following coupled problems
%
% -(2 exp(x) u')' = 10, 0 < x < 0.5,
% u(0) = 0, u(0.5) = v(0.2)
%
% -(3 v')' = 20,        0 < x < 0.2) 
% v(0) = 0, v'(0.2) = u'(0.5).
%
% which are related by the following compatibility conditions for u(x)
% and v(x)
%
% u(0) = v(0) = 0,
% u(0.5) = v(0.2) = alpha,
% u'(0.5) = v'(0.2) = beta
%
% with alpha, beta parameters to be determined.
%
% We apply the FEM to both equations using the only one linear element 
% in each case and let us assume that K u = F + Q and
%  L v = G + R  are the respective element associated systems
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART A
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The value of K(2,2) is,

syms x Psi1 Psi2 a1;
Psi2 =  @(x) 2*x;
a1 = @(x) 2*exp(x);
K22 = int(a1(x)*diff(Psi2(x),x)*diff(Psi2(x),x),0,0.5);

fprintf('(a) K(2,2) = %.4e\n',K22)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART B
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The value of F_{2} is

F2 = int(10*Psi2(x),0,0.5);
fprintf('(b) F(2) = %.4e\n',F2)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART C
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Taking into account the compatibility conditions (value and derivative
% of u and v) at the end nodes of each element, the value α we obtain for 
% the common node is,
% 
% Sol. Writng down the element equations and setting the corresponding 
% B.C. we have, from the 1st system 
%
% K12 * alpha = Q1 + F1,
% K22 * alpha = 2*exp(0.5) * beta + F2
%
% and, from the second
%
% L12 * alpha = R1 + G1,
% L22 * alpha = 3 * beta + G2
%
% So, taking the second equation of both systems and re-arranging the terms
% we get a linear system in the unknowns alpha, beta. Indeed,
%
% K22 * alpha - 2*exp(0.5) * beta = F2,
% L22 * alpha - 3 * beta          = G2,
% 
% K22 and F2 have been found in parts (a) and (b) whereas L22 and G2 are,
% respectively, the components of the stiffness matrix and the load vector
% of a BVP with constant coefficients, so we can find them using the 
% formulas in the course notes. 

syms Phi2
Phi2 = @(x) 5*x;
L22 = int(3*diff(Phi2(x),x)*diff(Phi2(x),x),0,0.2);
G2 = int(20*Phi2(x),0,0.2);

A = [K22, -a1(0.5); L22, -3];
b = [F2; G2];
y = eval(A\b);

fprintf('alpha = %.4e\n',y(1))
