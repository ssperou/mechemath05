% Runs numerical integration arbitrary RK method using variable time steps
%INPUTS:
% step_failure_rate: the amount of failed integrations/all integration
% attempts
% rate_func_in: the function used to compute dXdt. rate_func_in will
%              have the form: dXdt = rate_func_in(t,X) (t is before X)
% tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
% X0: the vector describing the initial conditions, X(t_start)
% h_ref: the desired value of the average step size (not the actual value)
% BT_struct: a struct that contains the Butcher tableau
%            BT_struct.A: matrix of a_{ij} values
%            BT_struct.B: vector of b_i values
%            BT_struct.C: vector of c_i values
% p: how error scales with step size (error = k*hˆp)
% error_desired: the desired local truncation error at each step
%OUTPUTS:
% t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is 
%         approximated at
% X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
% h_avg: the average step size
% num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg,fun_count,failure_rate] = explicit_RK_variable_step_integration ...
(func,tspan,X0,h_ref,BT_struct,p,error_desired,growth_limit)
    % initialize parameters
    t0 = tspan(1); tf = tspan(2);
    t = t0;
    h = h_ref;
    fun_count = 0;

    % create empty lists for saving timesteps and X values
    % including the first value
    t_list = [t];
    X_list = [X0];
    h_list = [];
    fail_count = 0;

    t_next = t+h;

    if t_next>tf
        t_next = tf;
        h = tf-t;
    end
    
    % loop until we reach the final time value
    while t < tf
        % append the new x-value at every loop
        XA = X_list(:,end);
        % solve the step function
        [XB, evals, h_next, redo] = explicit_RK_variable_step(func,t,XA,h, ...
            BT_struct,p,error_desired,growth_limit);
    
        % if we don't have to redo, update the x-value and current time;
        % otherwise, continue attempting to solve for the same time value
        if redo
            fail_count = fail_count+1;
        else
            X_list(:,end+1) = XB';            
            t = t_next;
            t_list(end+1) = t;
            h_list(end+1) = h;
        end

        % update the function call count, timestep size, and timestep
        fun_count = fun_count + evals;
        h = h_next;
        t_next = t+h;

        % if updating the time would put us over the final time value,
        % set the next time step to tf-t instead of h_next
        if t_next > tf
            t_next = tf;
            h = tf-t;
        end    
    end

    % calculate the average timestep size
    h_avg = mean(h_list);

    % calculate the failure rate
    failure_rate = fail_count/length(t_list);

end