function [x,fval,flag] = NewtonMethod(f,df,x0,tol,maxIter)
%NewtonMethod finds the root of a one-dimensional function.
%   NewtownMethod attempts to find a single root of a function with an
%   input function and function derivative. Uses the Newton-Raphson method
%   to itereate from an initial guess search for a root until terminiation
%   criteria are met.
%
%   [X,FVAL] = NewtonMethod(F,DF,X0) inputs an anonymous function F and an
%   anonymous function for the derivative DF and an intial guess for x X0 
%   and returns the location of the root and the value of the function 
%   evaluated at the root.
%
%   [X,FVAL] = NewtonMethod(F,DF,X0,TOL) allows the user to specify the
%   termination criteria tolerance TOL where at when abs(fval) < TOL,
%   iteration is terminated and the current x value returned.
%
%   [X,FVAL] = NewtonMethod(F,DF,X0,TOL,MAXITER) allows the user to
%   specify MAXITER, the maximum number of iteration the algorithm can run 
%   before terminating.
%
%   [X,FVAL,FLAG] = NewtonMethod(F,DF,X0,TOL,MAXITER) returns FLAG
%   indicating if the algorithm diverged. If the algorithm divereged, 
%   FLAG returns 1, if the divergence was not detected, FLAG returns 0. If
%   the algorithm diverges, the initial guess will be returned.
%
%   Ian Elliott 09/13/2017  - Init
%
if (nargin < 4 || isempty(tol)), tol = 1e-15; end
if (nargin < 5 || isempty(maxIter)), maxIter = 100; end

x = x0;
for ii = 1:maxIter
    % Find next point
    x = x - f(x)/df(x);
    fval = f(x);
    
    % Check tolerance
    if abs(fval) < tol
        break
    end
end

% Check for divergence
flag = 0;
if abs(fval) > abs(f(x0))
    warning('Divergence detected. Returning initial guess.')
    x = x0;
    fval = f(x0);
    flag = 1;
end

end
