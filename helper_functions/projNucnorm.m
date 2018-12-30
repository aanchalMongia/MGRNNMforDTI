function varargout = projNucnorm(B,m,n,radius,projTol,projData)
%% Computes the projection of an m x n matrix B onto ||X||_* <= radius
%
% Usage:   [p,projData] = projNucnorm(B,m,n,radius,projTol,projData)
% Inputs:  B - vectorized version of the matrix
%          m - number of rows
%          n - number of columnds
%          radius - radius of nuclear-norm ball
%          projTol - optional tolerance parameter (not used)
%          projData - optional projection data (not used)
% Outputs: p - projection of B
%          projData - optional projection data (not used)
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
try
[U,S,V] = svd(reshape(B,m,n),'econ');

s1 = diag(S);
s2 = oneProjectorMex(s1,radius);

p = U*diag(s2)*V';
p = p(:);

% Set output
varargout{1} = p;

catch
            fprintf('SVD didnt converge!')
end
        
if (nargout > 1)
   projData    = []; % Dummy solution
   varargout{2} = projData;
end
