function opt = setOptions(opt,varargin)
%% Set field values in options structure
%
%  opt = SetOptions(opt,varargin)
%
%  Varargin contains alternately a string representing the field name, and
%  the default value to assign to it. The updates for all fields that
%  already appear in the options structure are skipped.

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

for i=1:2:length(varargin)
   if ~isfield(opt,varargin{i})
      opt = setfield(opt,varargin{i},varargin{i+1});
   end
end
