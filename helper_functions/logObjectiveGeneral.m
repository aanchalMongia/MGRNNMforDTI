function varargout = logObjectiveGeneral(x,y,idx,f,fprime)
%% Computes the negative log-likelihood function and its gradient
%
% Usage:   [F,G] = logObjectiveGeneral(x,y,idx,f,fprime)
% Inputs:  x - vectorized version of the matrix
%          y - observations
%          idx - vector of indices corresponding to indices of observations
%          f - link function (f(x[i]) gives the probability that y[i]=1)
%          fprime - gradient of f with respect to x
% Outputs: F - the negative log-likelihood of x
%          G - the gradient of F (optional)
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

if nargout > 0
   F = -sum(log(y(idx).*f(x(idx)) - (y(idx)-1)/2));
   varargout{1} = F;
end

if nargout > 1
   G = zeros(size(x));
   v = (f(x(idx))+(y(idx)-1)/2);
   w = -fprime(x(idx));
   G(idx) = w./v;
   varargout{2} = G;
end
