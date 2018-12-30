function varargout = projectKappaTau(B,m,n,kappa,tau, projTol, projData)
%% Computes the projection of an m x n matrix B onto ||X||_* <= kappa, ||X||_inf <= tau
%
% Usage:   [p,projData] = projectKappaTau(B,m,n,kappa,tau,projTol,projData)
% Inputs:  B - vectorized version of the matrix
%          m - number of rows
%          n - number of columnds
%          kappa - radius of nuclear-norm ball
%          tau - radius of infinity-ball
%          projTol - optional tolerance parameter 
%          projData - optional projection data 
% Outputs: p - projection of x
%          projData - optional projection data 
%
% Most recent change - 10/3/2013
%
% Copyright 2013, M. Davenport, Y. Plan, E. van den Berg, M. Wootters
%
% This file is part of 1BMC Toolbox version 1.2.
%
%    The 1BMC Toolbox is free software: you can redistribute it and/or 
%    modify it under the terms of the GNU General Public License as 
%    published by the Free Software Foundation, either version 3 of the 
%    License, or (at your option) any later version.
%
%    The 1BMC Toolbox is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with the 1BMC Toolbox. If not, see <http://www.gnu.org/licenses/>.

if nargin < 6, projTol = 1e-8; end
if nargin < 7, projData = []; end

% ------------------------------------------------------------------------
% Original problem
% ------------------------------------------------------------------------
% Minimize     1/2 * ||X - B||_F^2
% Subject to   ||X||_inf <= tau, ||X||_* <=kappa
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% Reformulation
% ------------------------------------------------------------------------
% Minimize     1/2 * ||X1 - B||_F^2
% Subject to   ||X1||_inf <= tau, ||X2||_* <= kappa, X1=X2
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% Augmented Lagrangian
% ------------------------------------------------------------------------
% L(X1,X2,Y) = 1/2 * ||X1 - B||_F^2 - <Y,X1-X2> + mu/2*||X1-X2||_F^2
%
% Subject to: ||X1||_inf <= tau, ||X2||_* <= kappa
% ------------------------------------------------------------------------


% Initialization
if (~isempty(projData))
   Y  = projData.Y;
   X1 = projData.X1;
   X2 = projData.X2;
   mu= projData.mu;
else
   Y  = zeros(m*n,1);
   mu = 0.998;
   X1 = zeros(m*n,1); % max(min(B,tau),-tau);
   X2 = zeros(m*n,1);
end

mnMin = min(m,n);
rank  = min(m,n);

% Check feasibility
if (max(abs(B(:))) <= tau)
   [U,S,V] = svd(reshape(B,m,n),'econ');
   s1 = diag(S);
   if (sum(s1) <= kappa)
      X = B;
      
      varargout{1} = X;
      if (nargout > 1), varargout{2} = []; end;
      return;
   end
end

for i=1:1000

   % ---------------------------------------------------------------------
   % Infinity norm
   % ---------------------------------------------------------------------
   % Minimize    1/2 * ||X1 - B||_F^2 - <Y,X1> + mu/2*||X1-X2||_F^2
   %
   % Minimize    1/2 * ||sqrt(1+mu)*X1 - B/sqrt(1+mu) - Y/sqrt(1+mu) - mu*X2/sqrt(1+mu)||_F^2
   %
   % Minimize    1/2 * ||X1 - (B + Y + mu*X2) / (1+mu))||_F^2
   %    X1
   % Subject to  ||X1||_inf <= tau
   % ---------------------------------------------------------------------
   X1 = max(min( (B+Y+mu*X2)/(1+mu) ,tau),-tau);

   % ---------------------------------------------------------------------
   % Nuclear norm
   % ---------------------------------------------------------------------
   % Minimize     1/mu <Y,X2> + 1/2*||X1-X2||_F^2
   % Minimize     1/2*||X1 - 1/mu*Y - X2||_F^2
   % Minimize     1/2*||X2 - (X1 - 1/mu*Y)||_F^2
   %    X2
   % Subject to   ||X2||_* <= kappa
   % ---------------------------------------------------------------------

   % Projection -- Method #1
   [U,S,V] = svd(reshape(X1 - Y/mu,m,n),'econ');
   s1 = diag(S);
   s2 = oneProjectorMex(s1,kappa);
   X2 = U*diag(s2)*V';
   X2 = reshape(X2,m*n,1);
   
   % Projection -- Method #2
   %{
   A = X1-Y/mu; rank = mnMin;
   while (true)
      singVals = min(rank+3,mnMin);
      [U,S,V] = propack.lansvd(A,singVals,'L');
      s1 = diag(S);
      s2 = oneProjectorMex(s1,kappa);
      
      if ((singVals == mnMin) || (s2(end) < 1e-8))
         break;
      else
         rank = min(rank + 2*(singVals - rank),mnMin);
      end
   end
   rank = sum(s2 < 1e-8);
   X2 = U*diag(s2)*V';
   %} 
   

   % Compute difference between X1 and X2
   xFeas  = max(abs(X2));
   xDiff  = norm(X1-X2,'fro');
   xNorm  = norm(X1,'fro');
   infeas = max(xFeas,tau)-tau;
   relInf = infeas / max(1,tau); % Relative infeasibility

   %fprintf('%3d   %9.3e    %9.3e\n',i,xDiff/max(1,xNorm),relInf);
   
   if ((i > 10) && (xDiff / max(1,xNorm) < projTol) && (relInf < projTol))
      break;
   end;

   % Update Lagrange multiplier and increase mu
   Y = Y - mu * (X1-X2);
   
   %mu = mu * 1.1;
   mu = mu * 1.05;
end

X = (X1+X2)/2;

% Set ouput
varargout{1} = X;

if (nargout > 1)
   projData    = struct();
   projData.X1 = X1;
   projData.X2 = X2;
   projData.Y  = Y;
   projData.mu = mu;
   
   varargout{2} = projData;
end
