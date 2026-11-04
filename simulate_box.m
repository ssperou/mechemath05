% simulates a vibrating box
%INPUTS: None
%OUTPUTS: None
function simulate_box(vidflag)
    %define system parameters
    box_params = init_box_params();

    % create a wrapper function
    func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);

    % set up box parameters
    x0 = 0;
    y0 = 0;
    theta0 = 0;
    vx0 = 1;
    vy0 = 3;
    vtheta0 = 0.01;
    
    % set up integration parameters
    tspan = [0 10];
    V0 = [x0;y0;theta0;vx0;vy0;vtheta0];
    href = 0.01;
    [BT_struct,p] = init_BT("Fehlberg");
    error_desired = 1e-3;
    growth_limit = 3;

    % run the integration
    [tlist,Vlist] = explicit_RK_variable_step_integration(func,tspan,V0, ...
                        href,BT_struct,p,error_desired,growth_limit);

    % might move this to a different piece of code
    % taken from kai-xing's & my strandbeest animation code
    if vidflag
        % define location and filename where video will be stored
        mypath1 = 'C:\Users\jshen\OneDrive - Olin College of Engineering\2026-27\MechEMath\mechemath05\';
        fname='washing_machine';
        input_fname = [mypath1,fname];
        
        %create a videowriter, which will write frames to the animation file
        writerObj = VideoWriter(input_fname);
        
        %must call open before writing any frames
        open(writerObj);
    end
    
    % initialize the current figure and save as object
    clf
    fig1 = figure(1);

    %set up the plotting axis
    hold on;
    % axis([-115,65,-115,65])
    % axis equal;
    xlabel("$X$ Position (-)", Interpreter="latex");
    ylabel("$Y$ Position (-)", Interpreter="latex");
    title("Box Simulation", Interpreter="latex")
    set(gca,'TickLabelInterpreter','latex')
    % xticks([-100 -80 -60 -40 -20 0 20 40 60])
    fontsize(scale=2)
    % forces the video to be high-quality
    set(fig1,'units','pixels','position',[0 0 1440 1080])

    %%%%%%%% this section taken from orion's spring code
    % initialize spring 
    num_zigs = 5;
    w = .1;
    hold on;
    spring_plot_struct = initialize_spring_plot(num_zigs,w);
    axis equal; axis square;
    axis([-3,3,-3,3]);
    for frames = 1:length(V_list)
        P1 = V_list(1,frames);
        P2 = V_list(2,frames);
        update_spring_plot(spring_plot_struct,P1,P2)
        drawnow;
    end
    %%%%%%%%% end

        if vidflag
            %capture a frame (what is currently plotted)
            current_frame = getframe(fig1);
            %write the frame to the video
            writeVideo(writerObj,current_frame);
            if t == frames
                %must call close after all frames are written to save the video
                close(writerObj);
                length(xfoot)
            end
        end
end


%updates spring plotting object so that spring is plotted
%with ends located at points P1 and P2
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
%num_zigs: number of zig zags in spring drawing
%w: width of the spring drawing
function spring_plot_struct = initialize_spring_plot(num_zigs,w)
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

