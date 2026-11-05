% root finding function via Newton's method for multiple dimensions
%INPUTS:
%   fun: the function we are computing the root of
%   Note that fun(x) should output [f,dfdx], where dfdx is the derivative of f
%   (see test_func01 below for example)
%   x0: initial guess for Newton's method
%   dxmin: termination threshold (stop when interval abs(x_{i+1}-x_i) < dxtol)
%   ftol: termination threshold (stop when abs(f(x_{i}))<ftol)
%   max_iter: maximum iteration limit
%   dxmax: threshold checking that the jacobian isn't singular, which is bad
%   terminate when abs(x_{i+1}-x_i)>dxmax, where dxmax is a very large number
%   numerical_diff: boolean (0 or 1)
%           1-> numerically differentiate, 0-> use analytical derivative
%OUTPUTS
%   x: estimate for root of fun
%   exit_flag: an integer indicating that the solver succeeded (1) or
%   failed (0)
function [X, exit_flag, fun_count] = newton_solver(fun,X0,numerical_diff,dxmin,ftol,max_iter,dxmax)
    arguments
        fun (1,:) function_handle
        X0 (:,1) double
        numerical_diff (1,:) double = 1
        dxmin (1,:) double = 10e-14;
        ftol (1,:) double = 10e-14;
        max_iter (1,:) double = 200;
        dxmax (1,:) double = 10e5;
    end

    iter = 0;                  % set iteration variable
    X = X0;

    iteration_flag = true;     % set flags True to start the while loop
    interval_flag = true;
    value_flag = true;
    fun_count = 0;
    % loop through newton's method until the root is found, or until
    % the iteration maximum is hit
    while interval_flag && value_flag && iteration_flag

        if numerical_diff == 0
            [F,J] = fun(X);
        else 
            [J, num_evals] = approximate_jacobian(fun,X);
            fun_count = fun_count + num_evals;
            F = fun(X);
        end

        fun_count = fun_count + 1;
        X0 = X;   
        dx = J\F;     
        X = X0 - dx;               % otherwise continue computing

        % check for a denominator of zero or zero determinant of J
        if any(dxmax < abs(dx)) || det(J*J') < dxmin       
            % disp("Zero denominator error, or oversized update step size.")
            exit_flag = 0;              % if true: mark failure and exit
            return                      % the program
        end

        % add one to the iteration
        iter = iter+1;

        % break while loop if...
        % iterations are too close together (guesses too far left/right)
        % final value guess is 'close enough' to zero
        % maximum iteration values reached
        % or iterations are too far apart
        interval_flag = any(dxmin < abs(X - X0));
        value_flag = any(ftol < abs(F));
        iteration_flag = max_iter > iter;

        % else: continue while loop to try and find roots

    end

    % success is based on whether the final value is 'close enough' to
    % zero. set exit_flag correspondingly
    
    if ~value_flag
        exit_flag = 1;
        return
    elseif ~interval_flag
        %disp("Iterations tending towards a false root.")
        exit_flag = 0;
        return
    elseif ~iteration_flag
        %disp("Maximum iterations reached.")
        exit_flag = 0;
        return
    else
        disp("now you've gone and done it")
    end

   
end