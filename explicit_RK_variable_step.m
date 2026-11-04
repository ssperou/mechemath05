% Function to computes the value of X at the next time step
% for any arbitrary embedded RK method
% Also computes the next step size to use, and whether or not to accept/
% reject the projected value of X(t+h) for variable time step methods
%INPUTS:
% rate_func_in: the function used to compute dXdt. rate_func_in will
%               have the form: dXdt = rate_func_in(t,X) (t is before X)
% t: the value of time at the current step
% XA: the value of X(t)
% h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
% BT_struct: a struct that contains the Butcher tableau
%            BT_struct.A: matrix of a_{ij} values
%            BT_struct.B: vector of b_i values
%            BT_struct.C: vector of c_i values
% p: how error scales with step size (error = k*hˆp)
% error_desired: the desired local truncation error at each step
%OUTPUTS:
% XB: the approximate value for X(t+h)
% num_evals: A count of the number of times that you called
%            rate_func_in when computing the next step
% h_next: the time-step size at the next iteration
% redo: False if the estimated error was less than error_desired
%       True if the estimated error was larger than error_desired
function [XB, fun_count, h_next, redo] = explicit_RK_variable_step...
(rate_func_in,t,XA,h,BT_struct,p,error_desired,growth_limit)
    % run the method and calculate the output error
    [XB1, XB2, fun_count] = embedded_RK_step(rate_func_in,t,XA,h,BT_struct);
    error = norm(XB1 - XB2);

    % set h_next based on the error ratio and error scaling p, or the
    % growth limit
    h_next = min(0.9*(error_desired/error)^(1/p),growth_limit)*h;

    XB = XB1;
    redo = error>error_desired;
end