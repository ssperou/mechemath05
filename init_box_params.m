% creates a struct describing a box in a spring system
%INPUTS: None
%OUTPUTS:
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
function box_params = init_box_params()
    % define box characteristics
    % i don't know what LG is. arbitrary spring constant thing?
    % this code has not been touched compared to Orion's example except
    % that i made it a square
    LW = 2; LH = 2; LG = 5;

    % define system parameters based on box characteristics
    m = 1; Ic = (1/12)*(LH^2+LW^2);
    g = 1; k = 20; k_list = [.5*k,.5*k,2*k,5*k];
    l0 = 1.5*LG;

    % define box mounting points
    Pbl_box = [-LW;-LH];
    Pbr_box = [LW;-LH];
    Ptl_box = [-LW;LH];
    Ptr_box = [LW;LH];
    boundary_pts = [Pbl_box,Pbr_box,Ptr_box,Ptl_box,Pbl_box];

    % define world mounting points
    % for low speed but sometimes outside of the world-box:
    Pbl_world = Pbl_box + [-LG;-LG];
    Ptl_world = Pbl_box + [LG;-LG];
    Pbr_world = Pbr_box + [LG;LG];
    Ptr_world = Pbr_box + [-LG;LG];

    % for high speed but always inside the world-box:
    Pbl_world = Pbl_box + [-LG;-LG];
    Pbr_world = Pbr_box + [LG;-LG];
    Ptl_world = Ptl_box + [-LG;LG];
    Ptr_world = Ptr_box + [LG;LG];
    P_world = [Pbl_world,Ptl_world,Pbr_world,Ptr_world];
    P_box = [Pbl_box,Ptl_box,Pbr_box,Ptr_box];
    
    % save to a struct
    box_params = struct();
    box_params.m = m;
    box_params.I = Ic;
    box_params.g = g;
    box_params.k_list = k_list;
    box_params.l0_list = l0*ones(size(P_world,2));
    box_params.P_world = P_world;
    box_params.P_box = P_box;
    box_params.boundary_pts = boundary_pts;
end