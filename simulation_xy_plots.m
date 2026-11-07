% compares the linear and nonlinear system oscillation
%INPUTS: None
%OUTPUTS: None
function simulation_xy_plots()
    % define system parameters
    box_params = init_box_params();
    x0 = 0;
    y0 = 0;
    theta0 = pi/6;
    vx0 = 0;
    vy0 = 0;
    vtheta0 = pi/4;

    % create a wrapper function
    func_box = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);
    func_newton = @(V_in) box_rate_func(0,V_in,box_params);

    % calculate V_eq
    V0 = [x0;y0;theta0;vx0;vy0;vtheta0];
    V_eq = newton_solver(func_newton,V0);
    
    % set up integration parameters:
    %small number used to scale initial perturbation
    epsilon = 0.1;
    V0 = V_eq + epsilon*[x0;y0;theta0;vx0;vy0;vtheta0];
    % time, reference step size, RK method, desired error, growth limit
    tspan = [0 10];
    href = 0.01;
    [BT_struct,p] = init_BT("Fehlberg");
    error_desired = 1e-7;
    growth_limit = 1;

    % run the integration for the nonlinear system
    [t_nonlinear,V_nonlinear] = explicit_RK_variable_step_integration(func_box,tspan,V0, ...
                    href,BT_struct,p,error_desired,growth_limit);

    % run the integration for the linear system
    J_approx = approximate_jacobian(func_newton, V_eq);
    func_jacobian = @(t_in, V_in) J_approx*(V_in-V_eq);
    [t_linear,V_linear] = explicit_RK_variable_step_integration(func_jacobian,tspan,V0, ...
                        href,BT_struct,p,error_desired,growth_limit);
    
    % plot x and y values
    figure();
    hold on;
    grid on;
    plot(t_nonlinear, V_nonlinear(1, :),DisplayName="Nonlinear X values");
    plot(t_nonlinear, V_nonlinear(2, :),DisplayName="Nonlinear Y values");
    plot(t_linear, V_linear(1, :),DisplayName="Linear X values");
    plot(t_linear, V_linear(2, :),DisplayName="Linear Y values");
    xlabel("$X$ Position (-)", Interpreter="latex");
    ylabel("$Y$ Position (-)", Interpreter="latex");
    title("System Position over Time", Interpreter="latex")
    set(gca,'TickLabelInterpreter','latex')
    legend("Interpreter", "Latex");
    axis([tspan(1) tspan(2) -0.6 1])


    V_eq
    norm(func_newton(V_eq))

end