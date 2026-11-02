function simulate_box()
    %define system parameters
    box_params = struct();
    box_params.m = %your code here
    box_params.I = %your code here
    box_params.g = %your code here
    box_params.k_list = %your code here
    box_params.l0_list = %your code here
    box_params.P_world = %your code here
    box_params.P_box = %your code here
    %load the system parameters into the rate function
    %via an anonymous function
    my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);
    x0 = %your code here
    y0 = %your code here
    theta0 = %your code here
    vx0 = %your code here
    vy0 = %your code here
    vtheta0 = %your code here
    
    V0 = [x0;y0;theta0;vx0;vy0;vtheta0];
    tspan = %your code here
    %run the integration
    % [tlist,Vlist] = your_integrator(my_rate_func,tspan,V0,...)
end