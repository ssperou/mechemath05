% computes the linear and angular acceleration of the box given its current 
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
%   box_params.p_world: 2 x n list of static mounting points for the spring
%                       (in the world frame)
%   box_params.p_box: 2 x n list of mounting points for the spring (in the 
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
    p_box = box_params.P_box;
    p_world = box_params.P_world;
    
    % rotate box to our world!!!!
    p_list_world = compute_rbt(x,y,theta,p_box);

    % initialize sum of external spring force and torque
    F_ext = [0;0]; 
    T_ext = 0;


    % loop thru each spring, calc sum of forces
    for i = 1:length(l0_list)
        k = k_list(i);
        l0 = l0_list(i);

        PA = p_world(:,i); % fixed mount location
        PB = p_list_world(:,i); % box world mount location
        
        % calc force of this spring and add to sum
        F = compute_spring_force(k, l0, PA, PB);
        F_ext = F_ext + F;

        % calc torque and add to sum
        p_sub = PB - [x;y];
        T = p_sub(1)*F(2) - p_sub(2)*F(1);
        T_ext = T_ext + T;
    
    end
    % add gravity to net force
    F_net =  [0; -m * g] + F_ext;

    % calc x,y,theta accelerations
    ax = F_net(1)/m;
    ay = F_net(2)/m;
    atheta = T_ext/I;
end