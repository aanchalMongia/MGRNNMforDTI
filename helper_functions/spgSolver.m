function [x,info] = spgSolver(funObj, funProj, x0, options)
%% Spectral projected-gradient solver for convex constrained problems
%
% [x,info] = spgSolver(funObj, funProj, x0, options)
%
% solves the convex constrained problem
%
%   minimize  funObj(x)  subj to  x in C
%
% where C is a convex set implicitly defined by funProj. Parameter
% FUNOBJ must be a function handle that takes as input an X and outputs
% either the function value, or the function value and gradient, depending
% on the number of output arguments. FUNPROJ takes as input an X and
% returns the Euclidean projection onto the feasible set C, which is
% implicitly defined by this projection.

% options
% .verbosity    0 = No output
%               1 = Major output
%               2 = All output

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

% ----------------------------------------------------------------------
% This file is derived from SPGL1 <http://www.cs.ubc.ca/~mpf/spgl1/>
% Copyright 2007, Ewout van den Berg and Michael P. Friedlander,
%
% [1] E. van den Berg and M. P. Friedlander, Probing the Pareto frontier
%     for basis pursuit solutions, SIAM J. on Scientific Computing,
%     31(2):890-912, November 2008
% ----------------------------------------------------------------------


% Tic-safe toc for compatibility with older Matlab versions
try, t0 = toc; catch, tic, t0 = toc; end;

%----------------------------------------------------------------------
% Check arguments. 
%----------------------------------------------------------------------
if nargin < 3, error('Too few input arguments'); end
if nargin < 4 || ~exist('options','var'), options = struct(); end;

n = length(x0);

% Process options
options = setOptions(options, ...
  'fid'       ,     1 , ... % File ID for output
  'prefix'    ,    '' , ... % Output prefix
  'verbosity' ,     2 , ... % Verbosity level
  'iterations',  10*n , ... % Max number of iterations
  'nPrevVals' ,     3 , ... % Number previous func values for linesearch
  'optTol'    ,  1e-5 , ... % Optimality tolerance in x
  'dxTol'     ,     0 , ... % Optimality tolerance in dx
  'maxFunEval',   Inf , ... % Maximum number of function evaluations
  'maxRuntime',   Inf , ... % Maximum runtime allowed (in seconds)
  'callback'  ,    [] , ... % Callback function for optimality test
  'stepMin'   , 1e-16 , ... % Minimum spectral step
  'stepMax'   , 1e+05   ... % Maximum spectral step
);
fid           = options.fid;
logLevel      = options.verbosity;
maxIts        = options.iterations;
maxRuntime    = options.maxRuntime;
nPrevVals     = options.nPrevVals;
optTol        = options.optTol;
dxTol         = options.dxTol;
stepMin       = options.stepMin;
stepMax       = options.stepMax;
prefix        = options.prefix;

maxLineErrors = 20;   % Maximum number of line-search failures.

%----------------------------------------------------------------------
% Initialize local variables.
%----------------------------------------------------------------------
iter          = 0; % Total iterations.
stat          = 0;
nProjections  = 0; % Number of projections
nFunctionEval = 0; % Number of function evaluations
nGradientEval = 0; % Number of gradient evaluations
nLineTot      = 0; % Total number of linesearch steps
timeProject   = 0; % Projection time
timeFunEval   = 0; % Function/gradient evaluation time
lastFv        = -inf(nPrevVals,1);  % Last m function values
usecallback   = ~isempty(options.callback);

% Exit conditions (constants).
EXIT_OPTIMAL1         = 1;
EXIT_OPTIMAL2         = 2;
EXIT_ITERATIONS       = 3;
EXIT_LINE_ERROR       = 4;
EXIT_FUNEVAL_LIMIT    = 5;
EXIT_CALLBACK         = 6;
EXIT_RUNTIME_LIMIT    = 7;

%----------------------------------------------------------------------
% Log header.
%----------------------------------------------------------------------
logB = '%s %5i  %13.7e  %9.2e  %9.2e  %9.2e  %5.1e  %2d  %s\n';
logH = '%s %5s  %13s  %9s  %9s  %9s  %6s  %2s  %s\n';
if logLevel > 0
   %%fprintf(fid,'%s\n', prefix);
   %%fprintf(fid,'%s %s\n', prefix,repmat('=',1,80));
   %%fprintf(fid,'%s SPG-Solver\n', prefix);
   %%fprintf(fid,'%s %s\n', prefix, repmat('=',1,80));
   %%fprintf(fid,'%s %-22s: %8.2e %4s %-22s: %8i\n', prefix, 'Optimality tol',optTol,'','Maximum iterations',maxIts);
end

if logLevel > 1
   %%%fprintf(fid,'%s \n', prefix);
   %%%fprintf(fid, logH, prefix,'Iter','Objective','xNorm','dxNorm','gNorm','stepG','LS','Comment');
end

%----------------------------------------------------------------------
% Final preparations.
%----------------------------------------------------------------------

% Project the starting point and evaluate function and gradient.
x     = project(x0);
[f,g] = objective(x);

% Required for nonmonotone strategy.
lastFv(1) = f;
fBest     = f;
xBest     = x;
fOld      = f;

% Compute projected gradient direction and initial steplength.
dx        = project(x - g) - x;
dxNorm    = norm(dx,inf);
dx0Norm   = dxNorm;
gStep     = 1;
lnErr     = 0;
stepG     = 1; % Line search length

history   = zeros(maxIts,5);

projData  = [];
projTol   = 1e-8;  % Projection tolerance must start relatively small;
reproject = false; % otherwise we may find an (infeasible) point whose
                   % objective we cannot improve on.

bbMult = 1;                   
                   
%----------------------------------------------------------------------
% MAIN LOOP.
%----------------------------------------------------------------------
while 1
    commentStr = '';
    
    %------------------------------------------------------------------
    % Test iterate-independent exit conditions.
    %------------------------------------------------------------------
    if usecallback
       [cbstat,data] = options.callback(x,data);
       if cbstat ~= 0
          stat = EXIT_CALLBACK;
       end
    end
    
    if (toc - t0) > maxRuntime
       stat = EXIT_RUNTIME_LIMIT;
    end

    if iter >= maxIts
       stat = EXIT_ITERATIONS;
    end


    %------------------------------------------------------------------
    % Test for unconstrained minimum.
    %------------------------------------------------------------------
    gNorm = norm(g,2);
    xNorm = norm(x,2);
    if (gNorm / max(xNorm,1) < 1e-8) % TODO: Value 1e-8
       
       % Do not exit if gradient is tiny but change in x is still large
       if (projTol > 1e-10)
          projTol = 1e-10;
          commentStr = 'Optimal 1 - tightening proj. tol';
       elseif ((iter > 0) && (norm(s)/max(1,xNorm) < 1e-5))
          if (stat == 0)
             stat = EXIT_OPTIMAL1;
          end
       end
    end
    
    % -----------------------------------------------------------------
    % Compute scaled projected gradient; make sure step is large enough
    % and check optimality conditions
    % -----------------------------------------------------------------
    if (stat == 0)
       if (reproject == false)
          projData = [];
       end
       
       reproject = false;
       [nx,projData] = project(x - max(gStep,1)*g, projTol, projData);
       dx        = nx - x;
       gtd       = real(g'*dx); % Compute inner product
       dxNorm    = norm(dx,2);
       stepG     = 1;
       
   
       % --------------------------------------------------------------
       % Check whether a maximum step would get us anywhere. There are
       % a number of situations to consider:
       % 1. ||g|| is small and we are near optimal (checked above)
       % 2. ||dx|| is small; this will be the case for (1) but could
       %    also happen when x is at the boundary of the feasible set
       %    and -g is (nearly) perpendicular to the supporting
       %    hyperplane of the set.
       % 3. When <g,d> >= 0 we should have that ||dx|| is small,
       %    unless of course the projection is so inaccurate that
       %    <g,d> >> 0.
       % -----------------------------------------------------------------
       if ((gtd >= 0) || (dxNorm < 1e-9 * max(xNorm,1)))
          if (projTol <= 1.1e-10)
             % Projection is presumed to be sufficiently accurate
             bbMult = bbMult * 10; gStep = gStep * 10;
             if (bbMult >= 1e6)
                stat = EXIT_OPTIMAL2;
             else
                commentStr = sprintf('BBMult = %6.2e',bbMult);
                reproject = true;
             end
          else
             % Increate projection accuracy
             projTol = projTol * 0.1;
    
             commentStr = sprintf('Reprojecting (%6.2e)',projTol);
             reproject = true;
          end
       end
    end
    
    %------------------------------------------------------------------
    % Print log, update history and act on exit conditions.
    %------------------------------------------------------------------
    if logLevel > 1 || (stat ~= 0 && logLevel > 0)
       %%fprintf(fid, logB, prefix,iter,f,max(xNorm,1),dxNorm,norm(g,2),gStep,lnErr,commentStr);
    end
    history(iter+1,:) = [f,xNorm,dxNorm,norm(g,2),gStep];
    
    if (reproject), continue; end % Project with increased accuracy

    if stat ~= 0, break; end % Act on exit conditions.


    %==================================================================
    % Iterations begin here.
    %==================================================================
    iter = iter + 1;
    xOld = x;  fOld = f;  gOld = g;

    try
       %---------------------------------------------------------------
       % Projected gradient step and linesearch.
       %---------------------------------------------------------------
       [f,x,nLine,lnErr] = spgLine(x,dx,gtd,fOld,max(lastFv),@objective);
       nLineTot = nLineTot + nLine;

       if lnErr
          % Retry with curvilinear linesearch.
          x = xOld;
          f = fOld;

          % Curvilinear line search
          [f,x,nLine,stepG,lnErr] = ...
              spgLineCurvy(x,gStep*g,max(lastFv),projTol,@objective,@project);
       
          nLineTot = nLineTot + nLine;
       end
       
       if lnErr
           x = xOld;
           f = fOld;

           % Compute projection and display statistics
           dx  = project(x - gStep*g) - x;
           gtd = real(g'*dx);
           
           stat = EXIT_LINE_ERROR;
       end
       
       %---------------------------------------------------------------
       % Update gradient and compute new Barzilai-Borwein scaling.
       %---------------------------------------------------------------
       [f,g] = objective(x);
       s    = x - xOld;
       y    = g - gOld;
       sts  = s'*s;
       sty  = s'*y;
       if   sty <= 0, gStep = stepMax;
       else           gStep = min( stepMax, max(stepMin, sts/sty) );
       end

       % Update step length with BB multiplication factor.
       gStep = gStep * bbMult;

       
    catch % Detect function evaluation limit error
       err = lasterror;
       if strcmp(err.identifier,'SPG:MaximumFunEval')
         stat = EXIT_FUNEVAL_LIMIT;
         iter = iter - 1;
         
         % Restore previous iterate
         x = xOld;  f = fOld;  g = gOld;
         break;
       else
         rethrow(err);
       end
    end
    
    %%%fprintf(fid,'>> %9.2e (stat = %d)\n',norm(s,2),stat);
    
    %------------------------------------------------------------------
    % Update function history.
    %------------------------------------------------------------------
    lastFv(mod(iter,nPrevVals)+1) = f;
    if fBest > f
       fBest = f;
       xBest = x;
    end
    
end % while 1

% Restore best iterate
if f > fBest
   if logLevel > 0
      %%fprintf(fid,'%s\n', prefix);
      %%fprintf(fid,'%s Restoring best iterate to objective %13.7e\n', prefix,fBest);
   end   
   x = xBest;
   f = fBest;
end

% Truncate history information
history(iter+1:end,:) = [];

% Information structure
info = struct();
info.iter          = iter;
info.stat          = stat;
info.timeTotal     = toc - t0;
info.timeProject   = timeProject;
info.timeFunEval   = timeFunEval;
info.nFunctionEval = nFunctionEval;
info.nGradientEval = nGradientEval;
info.nProjections  = nProjections;
info.history       = history;

% Print final output.
switch (stat)
   case EXIT_OPTIMAL1
      info.statMsg = 'Optimal solution found - 1';
   case EXIT_OPTIMAL2
      info.statMsg = 'Optimal solution found - 2';
   case EXIT_ITERATIONS
      info.statMsg = 'Too many iterations';
   case EXIT_FUNEVAL_LIMIT
      info.statMsg = 'Maximum function evaluations reached';
   case EXIT_LINE_ERROR
      info.statMsg = 'Maximum line search iterations reached';
   case EXIT_CALLBACK
      info.statMsg = 'Callback function exit';
   case EXIT_RUNTIME_LIMIT
      info.statMsg = 'Maximum runtime reached';
   otherwise
      info.statMsg = 'Unknown termination condition';
end

if logLevel > 0
   %%fprintf(fid,'%s\n', prefix);
   %%fprintf(fid,'%s EXIT -- %s\n', prefix, info.statMsg)
   %%fprintf(fid,'%s\n', prefix);
   %%fprintf(fid,'%s %-20s:  %6i %6s %-20s:  %6.1f\n', prefix, 'Function evals.',nFunctionEval,'',  'Total time   (secs)',info.timeTotal);
   %%fprintf(fid,'%s %-20s:  %6i %6s %-20s:  %6.1f\n', prefix, 'Gradient evals.',nGradientEval,'', 'Project time (secs)', timeProject);
   %%fprintf(fid,'%s %-20s:  %6i %6s %-20s:  %6.1f\n', prefix, 'Iterations',info.iter,'', 'Fun eval time (secs)',timeFunEval);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED FUNCTIONS.  These share some vars with workspace above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ----------------------------------------------------------------------
function varargout = objective(varargin)
% ----------------------------------------------------------------------
   if (nFunctionEval >= options.maxFunEval)
     error('SPG:MaximumFunEval','');
   end
   tStart = toc;

   if nargout == 1
      varargout{1} = funObj(varargin{:});
   elseif nargout == 2
      [varargout{1},varargout{2}] = funObj(varargin{:});
      nGradientEval = nGradientEval + 1;
   end
   
   nFunctionEval = nFunctionEval + 1;
   
   timeFunEval = timeFunEval + (toc - tStart);
end % function objective

% ----------------------------------------------------------------------
function varargout = project(x,projTol,projData)
% ----------------------------------------------------------------------
   if nargin < 2, projTol  = 1e-9; end;
   if nargin < 3, projData = []; end

   tStart       = toc;
   if (nargout == 1)
      x = funProj(x,projTol,projData);
      varargout{1} = x;
   else
      [x,projData] = funProj(x,projTol,projData);
      varargout{1} = x;
      varargout{2} = projData;
   end
   timeProject  = timeProject + (toc - tStart);
   nProjections = nProjections + 1;
end % function project


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of nested functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % spgGeneral


% ----------------------------------------------------------------------
function [fNew,xNew,iter,stat] = spgLine(x,d,gtd,f,fMax,funObj)
% ----------------------------------------------------------------------
% Nonmonotone linesearch.

EXIT_CONVERGED  = 0;
EXIT_ITERATIONS = 1;

maxIts = 20;
step   =  1;
iter   =  0;
gamma  = 1e-4;

while 1
    % Evaluate trial point and function value.
    xNew = x + step*d;
    fNew = funObj(xNew);

    % Check exit conditions.
    if fNew < fMax + gamma*step*gtd  % Sufficient descent condition.
       stat = EXIT_CONVERGED;
       break
    elseif iter >= maxIts           % Too many linesearch iterations.
       stat = EXIT_ITERATIONS;
       break
    end

    % New linesearch iteration.
    iter = iter + 1;
    
    % Safeguarded quadratic interpolation.
    if step <= 0.1
       step  = step / 2;
    else
       tmp = (-gtd*step^2) / (2*(fNew-f-step*gtd));
       if tmp < 0.1 || tmp > 0.9*step || isnan(tmp)
          tmp = step / 2;
       end
       step = tmp;
    end
end % while 1

end % function spgLine


% ----------------------------------------------------------------------
function [fNew,xNew,iter,step,stat] = spgLineCurvy(x,g,f,projTol,funObj,funProj)
% ----------------------------------------------------------------------
% Projected backtracking linesearch.
% On entry d is the (possibly scaled) steepest descent direction.

EXIT_CONVERGED  = 0;
EXIT_ITERATIONS = 1;

maxIts = 20;
step   =  1;
sNorm  =  1;
scale  =  1; % Safeguard scaling (see below).
iter   =  0;
gamma  =  1e-4;

while 1
    % Evaluate trial point and function value.
    xNew = funProj(x - step*scale*g,projTol);
    
    fNew = funObj(xNew);
    s    = xNew - x;
    gts  = scale * real(g' * s); % This should be negative

    gtsBar = gts / (norm(g)*norm(s));
    %%%fprintf('Curvilinear %9.4e\n', gtsBar);
    
    %if gts < 0
    if (gtsBar < -1e-6) % The angle must be sufficiently large (TODO)
       if fNew < f + gamma*step*gts  % EvdB Jan 31, 2011: Changed from - to +
          stat = EXIT_CONVERGED;
          break;
       end
    end

    if iter >= maxIts
       stat = EXIT_ITERATIONS;
       break;
    end      

    % New linesearch iteration.
    iter = iter + 1;
    step = step / 2;

    % Safeguard: If stepMax is huge, then even damped search directions
    % can give exactly the same point after projection.  If we observe
    % this in adjacent iterations, we drastically damp the next search
    % direction.
    sNorm = norm(s,Inf);
    if sNorm <= 1e-6 * max(norm(x,2),1);
       scale = scale / 10;
    end
end % while 1

end % spgLineCurvy
