% simulates a vibrating box
%INPUTS: 
% vidflag: boolean describing whether or not a video should be saved
%OUTPUTS: None
function simulate_box(system, vidflag)
    % optional argument that defaults to off
    % if true or 1, will save video
    arguments
        system (1,:) string = "Nonlinear"
        vidflag (1,:) double = 0
    end
    % define system parameters
    box_params = init_box_params();
    p_box = box_params.P_box;
    p_world = box_params.P_world;
    x0 = 6;
    y0 = 6;
    theta0 = pi/6;
    vx0 = 4;
    vy0 = 4;
    vtheta0 = pi;

    % create a wrapper function
    func_box = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);
    func_newton = @(V_in) box_rate_func(0,V_in,box_params);

    % calculate V_eq
    V0 = [x0;y0;theta0;vx0;vy0;vtheta0];
    V_eq = newton_solver(func_newton,V0);
    
    % set up integration parameters:
    % set up integration parameters:
    %small number used to scale initial perturbation
    epsilon = 0.1;
    V0 = V_eq + epsilon*[x0;y0;theta0;vx0;vy0;vtheta0];
    % time, reference step size, RK method, desired error, growth limit
    tspan = [0 30];
    href = 0.01;
    [BT_struct,p] = init_BT("Fehlberg");
    error_desired = 1e-7;
    growth_limit = 1;

    % run the integration
    if system == "Nonlinear" || system == "Both"
        [tlist,V_nonlinear] = explicit_RK_variable_step_integration(func_box,tspan,V0, ...
                        href,BT_struct,p,error_desired,growth_limit);
    elseif system == "Linear" || system == "Both"
        J_approx = approximate_jacobian(func_newton, V_eq);
        func_jacobian = @(t_in, V_in) J_approx*(V_in-V_eq);
        [tlist,V_linear] = explicit_RK_variable_step_integration(func_jacobian,tspan,V0, ...
                        href,BT_struct,p,error_desired,growth_limit);
    end

    % simulate the nonlinear vibrating box
    % save as video if vidflag true
    if vidflag
        mypath1 = 'C:\Users\ccirone\Downloads\';
        fname='square_animation.avi';
        input_fname = [mypath1, fname];
    
        %create a videowriter, which will write frames to the animation file
        writerObj = VideoWriter(input_fname);
        open(writerObj); %must call open before writing any frames
    end

    % initialize the current figure and save as object
    fig1 = figure(1);
    clf

    %set up plot styling
    hold on;
    axis equal; axis square;
    axis([-10,10,-10,10]);
    xlabel("$X$ Position (-)", Interpreter="latex");
    ylabel("$Y$ Position (-)", Interpreter="latex");
    title("Box Simulation", Interpreter="latex")
    set(gca,'TickLabelInterpreter','latex')
    % xticks([-100 -80 -60 -40 -20 0 20 40 60])
    fontsize(scale=2)
    
    % forces the video to be high-quality
    set(fig1,'units','pixels','position',[0 0 1440 1080])

    %%%%%%%% this section taken from orion's spring code
    % initialize springs
    sp1 = initialize_spring_plot();
    sp2 = initialize_spring_plot();
    sp3 = initialize_spring_plot();
    sp4 = initialize_spring_plot();
    origin = plot(0, 0, 'black');

    for j = 1:length(tlist)
        % initialize parameters
        if system == "Nonlinear"
            x = V_nonlinear(1,j);
            y = V_nonlinear(2,j);
            theta = V_nonlinear(3, j);
        elseif system == "Linear"
            x = V_linear(1,j);
            y = V_linear(2,j);
            theta = V_linear(3, j);
        end
        p_box_reframed = compute_rbt(x, y, theta, p_box);

        % draw springs
        update_spring_plot(sp1,p_box_reframed(:, 1), p_world(:, 1));
        update_spring_plot(sp2,p_box_reframed(:, 2), p_world(:, 2));
        update_spring_plot(sp3,p_box_reframed(:, 3), p_world(:, 3));
        update_spring_plot(sp4,p_box_reframed(:, 4), p_world(:, 4));

        % draw box
        bound = box_params.boundary_pts;
        box_bound = compute_rbt(x, y, theta, bound);
        set(origin, 'xdata', box_bound(1, :), 'ydata', box_bound(2, :));

        drawnow;

        if vidflag
            %capture a frame (what is currently plotted)
            current_frame = getframe(fig1);
            %write the frame to the video
            writeVideo(writerObj,current_frame);
            if length(V_nonlinear) == j
                %must call close after all frames are written to save the video
                close(writerObj);
                length(xfoot)
            end
        end
    end
end


% updates spring plotting object so that spring is plotted
% with ends located at points P1 and P2
function update_spring_plot(spring_plot_struct,P1,P2)
    dP = P2-P1;
    R = [dP(1),-dP(2)/norm(dP);dP(2),dP(1)/norm(dP)];
    plot_pts = R*spring_plot_struct.zig_zag;
    set(spring_plot_struct.line_plot,...
    'xdata',plot_pts(1,:)+P1(1),...
    'ydata',plot_pts(2,:)+P1(2));
    set(spring_plot_struct.point_plot,...
    'xdata',[P1(1),P2(1)],...
    'ydata',[P1(2),P2(2)]);
end

%create a struct containing plotting info for a single spring
%INPUTS:
% num_zigs: number of zig zags in spring drawing
% w: width of the spring drawing
function spring_plot_struct = initialize_spring_plot(num_zigs,w)
    arguments
        num_zigs (1,:) double = 5
        w (1,:) double = 0.1
    end
    spring_plot_struct = struct();
    zig_ending = [.25,.75,1; ...
    -1,1,0];
    zig_zag = zeros(2,3+3*num_zigs);
    zig_zag(:,1) = [-.5;0];
    zig_zag(:,end) = [num_zigs+.5;0];
    for n = 0:(num_zigs-1)
        zig_zag(:,(3+3*n):2+3*(n+1)) = zig_ending + [n,n,n;0,0,0];
    end
    zig_zag(1,:)=(zig_zag(1,:)-zig_zag(1,1))/(zig_zag(1,end)-zig_zag(1,1));
    zig_zag(2,:)=zig_zag(2,:)*w;
    spring_plot_struct.zig_zag = zig_zag;
    spring_plot_struct.line_plot = plot(0,0,'k','linewidth',2);
    spring_plot_struct.point_plot = plot(0,0,'ro','markerfacecolor','r','markersize',7);
end

