function p = gausscdf(x,mu,sigma)

%% Computes the cdf of a Gaussian distribution at specified values 
%
% Usage:   p = gausscdf(x,mu,sigma)
% Inputs:  x - array of values
%          mu - mean of Gaussian
%          sigma - standard deviation of Gaussian
% Outputs: p - values of the CDF at x
%
% Most recent change - 05/16/2014
%
% Copyright 2014, M. Davenport, Y. Plan, E. van den Berg, M. Wootters
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

%% Standardize variables
z = (x-mu)/sigma;

%% Compute values of CDF
p = 0.5 * erfc(-z/sqrt(2));
