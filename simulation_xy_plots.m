% compares the linear and nonlinear system oscillation
%INPUTS: None
%OUTPUTS: None
function simulation_xy_plots(modal)
    % optional argument that defaults to off
    % if true or 1, will run modal analysis
    arguments
        modal (1,:) double = 0;
    end
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
    
    % set up different analyses
    if modal
        % loop through eigenvalues of part of the jacobian
        A = approximate_jacobian(func_newton, V_eq);
        Q = A(4:6, 1:3);
        [U_mode, omega_n] = eigs(-Q)
        omega_n = sqrt(omega_n)
        value_list = U_mode;
        epsilon = 1e-6;
    else
        % loop through magnification value
        epsilon = [1e-6, 1e-3, 1e-2];
        value_list = epsilon;
    end
    
    for i = 1:length(value_list)
        % time, reference step size, RK method, desired error, growth limit
        tspan = [0 10];
        href = 0.01;
        [BT_struct,p] = init_BT('Fehlberg');
        error_desired = 1e-7;
        growth_limit = 1;
        
        % if modal analysis
        if modal
            % set initial condition based on eigenvalues
            tspan = [0 5];
            V0 = V_eq + epsilon.*[U_mode(:,i); 0;0;0];

            [t_nonlinear,V_nonlinear] = explicit_RK_variable_step_integration(func_box,tspan,V0, ...
                                href,BT_struct,p,error_desired,growth_limit);

            x_modal = V_eq(1) + epsilon*U_mode(1,i)*cos(omega_n(i,i)*t_nonlinear);
            y_modal = V_eq(2) + epsilon*U_mode(2,i)*cos(omega_n(i,i)*t_nonlinear);
            theta_modal = V_eq(3) + epsilon*U_mode(3,i)*cos(omega_n(i,i)*t_nonlinear);
            V_linear = [x_modal; y_modal; theta_modal];
            t_linear = t_nonlinear;

            subt = sprintf('Mode-shape %g',i);
            label = 'Modal Approximation';
            
        % elseif regular analysis
        else
            % set up integration parameters:
            % small number used to scale initial perturbation
            V0 = V_eq + epsilon(i)*[x0;y0;theta0;vx0;vy0;vtheta0];
    
            % run the integration for the nonlinear system
            [t_nonlinear,V_nonlinear] = explicit_RK_variable_step_integration(func_box,tspan,V0, ...
                            href,BT_struct,p,error_desired,growth_limit);
        
            % run the integration for the linear system
            J_approx = approximate_jacobian(func_newton, V_eq);
            func_jacobian = @(t_in, V_in) J_approx*(V_in-V_eq);
            [t_linear,V_linear] = explicit_RK_variable_step_integration(func_jacobian,tspan,V0, ...
                            href,BT_struct,p,error_desired,growth_limit);

            subt = sprintf('Distance from Equilibrium: Factor of %s',epsilon(i));
            label = 'Linear Approximation';
        end
        
        % initialize the subplots
        figure(i);
        t = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
        title(t,'Position of a Vibrating Box','Interpreter','Latex','FontSize',24);
        subtitle(t,subt,'Interpreter','Latex','FontSize',20)
    
        % plot the x-values
        nexttile;
        hold on;
        grid on;
        plot(t_linear, V_linear(1, :), '-', 'LineWidth', 2);
        plot(t_nonlinear, V_nonlinear(1, :), ':', 'LineWidth', 2);
        set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 14);
        title('X Position','Interpreter','Latex', 'FontSize',18);
        xlabel('Time (s)','Interpreter','Latex', 'FontSize',15);
        ylabel('X (m)','Interpreter','Latex', 'FontSize',15)
        
        % plot the y-values
        nexttile;
        hold on;
        grid on;
        plot(t_linear, V_linear(2, :), '-', 'LineWidth', 2);
        plot(t_nonlinear, V_nonlinear(2, :), ':', 'LineWidth', 2);
        set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 14);
        title('Y Position','Interpreter','Latex', 'FontSize',18);
        xlabel('Time (s)','Interpreter','Latex', 'FontSize',15);
        ylabel('Y (m)','Interpreter','L', 'FontSize',15)
        
        % plot the orientation (theta-values)
        nexttile;
        hold on;
        grid on;
        plot(t_linear, V_linear(3, :), '-', 'LineWidth', 2, 'DisplayName', label);
        plot(t_nonlinear, V_nonlinear(3, :), ':', 'LineWidth', 2, 'DisplayName', 'Non-linear Approximation');
        set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 14);
        title('Orientation $\theta$', 'Interpreter','latex', 'FontSize',18);
        xlabel('Time (s)',  'Interpreter', 'Latex', 'FontSize', 15);
        ylabel('$\theta$ (rad)', 'Interpreter', 'Latex', 'FontSize', 15)

        % set legend
        my_legend = legend();
        my_legend.Layout.Tile = 'south';
        my_legend.ItemTokenSize(1) = 10;
        set(my_legend,'Interpreter','Latex','FontSize',16,'Orientation', 'Horizontal');


    end

    % % if Modal Analysis
    % A = J_approx;
    % Q = A(4:6, 1:3);
    % [U_mode, omega_n] = eigs(Q);
    % 
    % epsilon = 1;
    % V0 = V_eq + epsilon.*[U_mode; zeros(3,3)]
    % 
    % tspan = [0 2];
    % 
    % [t_nonlinear,V_nonlinear] = explicit_RK_variable_step_integration(func_box,tspan,V0, ...
    %                     href,BT_struct,p,error_desired,growth_limit);
    % 
    % x_modal = V_eq(1) + epsilon * U_mode(1) * cos(omega_n(:,1)*t_nonlinear);
    % y_modal = V_eq(2) + epsilon * U_mode(2) *cos(omega_n(:,2)*t_nonlinear);
    % theta_modal = V_eq(3) + epsilon * U_mode(3) *cos(omega_n(:,3)*t_nonlinear);
    % 
    %  % initialize the subplots for each mode shape
    % figure(i);
    % t = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    % title(t,'Position of a Vibrating Box','Interpreter','Latex','FontSize',24);
    % %subtitle(t,sprintf('Distance from Equilibrium: Factor of %s',epsilon(i)), ...
    %     %'Interpreter','Latex','FontSize',20)
    % 
    % % plot the x-values
    % nexttile;
    % hold on;
    % grid on;
    % plot(t_nonlinear, x_modal, '-', 'LineWidth', 2, 'DisplayName', 'Linearized');
    % plot(t_nonlinear, V_nonlinear(1, :), '.', 'LineWidth', 2, 'DisplayName', 'Non-linear');
    % set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 14);
    % title('X Position','Interpreter','Latex', 'FontSize',18);
    % xlabel('Time (s)','Interpreter','Latex', 'FontSize',15);
    % ylabel('X (m)','Interpreter','Latex', 'FontSize',15)
    % 
    % % plot the y-values
    % nexttile;
    % hold on;
    % grid on;
    % plot(t_linear, y_modal, '-', 'LineWidth', 2, 'DisplayName', 'Linearized');
    % plot(t_nonlinear, V_nonlinear(2, :), '.', 'LineWidth', 2, 'DisplayName', 'Non-linear');
    % set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 14);
    % title('Y Position','Interpreter','Latex', 'FontSize',18);
    % xlabel('Time (s)','Interpreter','Latex', 'FontSize',15);
    % ylabel('Y (m)','Interpreter','L', 'FontSize',15)
    % 
    % % plot the orientation (theta-values)
    % nexttile;
    % hold on;
    % grid on;
    % plot(t_linear, theta_modal, '-', 'LineWidth', 2, 'DisplayName', 'Linear Appriximation');
    % plot(t_nonlinear, V_nonlinear(3, :), '.', 'LineWidth', 2, 'DisplayName', 'Non-linear Approximation');
    % set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 14);
    % title('Orientation $\theta$', 'Interpreter','latex', 'FontSize',18);
    % xlabel('Time (s)',  'Interpreter', 'Latex', 'FontSize', 15);
    % ylabel('$\theta$ (rad)', 'Interpreter', 'Latex', 'FontSize', 15)
    % 
    % % set legend
    % my_legend = legend();
    % my_legend.Layout.Tile = 'south';
    % my_legend.ItemTokenSize(1) = 10;
    % set(my_legend,'Interpreter','Latex','FontSize',16,'Orientation', 'Horizontal');


end