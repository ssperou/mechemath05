% Computes the linear and angular acceleration of the box given its current 
% position and orientation
%INPUTS:
% x: current x position of the box
% y: current y position of the box
% theta: current orientation of the box
% box_params: a struct containing the parameters that describe the system
%   Fields:
%   box_params.m: mass of the box
%   box_params.I: moment of inertia w/respect to centroid
%   box_params.g: acceleration due to gravity
%   box_params.k_list: list of spring stiffnesses
%   box_params.l0_list: list of spring natural lengths
%   box_params.P_world: 2 x n list of static mounting points for the spring
%                       (in the world frame)
%   box_params.P_box: 2 x n list of mounting points for the spring (in the 
%                     box frame)
%OUTPUTS
% ax: x acceleration of the box
% ay: y acceleration of the box
% atheta: angular acceleration of the box
function [ax,ay,atheta] = compute_accel(x,y,theta,box_params)
   % initialize structs
    m = box_params.m;
    g = box_params.g;
    I = box_params.I;
    k_list = box_params.k_list;
    l0_list = box_params.l0_list;
    P_box = box_params.P_box;
    P_world = box_params.P_world;
    
    % rotate box to our world!!!!
    P_list_world = compute_rbt(x,y,theta,P_box);

    % initialize sum of springs, torque
    F_ext_sum = [0;0]; % ! vertical or horizontal?
    T_ext_sum = 0;


    % loop thru each spring, calc sum of forces
    for i = 1:length(l0_list)
        k_i = k_list(i);
        l0_i = l0_list(i);

        PA = P_world(:,i); % fixed mount
        PB = P_list_world(:,i); % box world mount
        
        % calc force of each spring, add to sum
        F_i = compute_spring_force(k_i, l0_i, PA, PB);
        
        % add force to sumnation
        F_ext_sum = F_ext_sum + F_i;

        % calc torque
        
        T_ext_sum = T_ext_sum + cross((PB - [x;y]), F_i) % #crossproductsbabyyyy # i think this format is wrong
    
    end
    % add gravity to net force
    F_net =  [0; -m * g] + F_ext_sum;

    % calc x,y,theta accelerations
    ax = F_net(1)/m;
    ay = F_net(2)/m;
    atheta = T_ext_sum/I;
end