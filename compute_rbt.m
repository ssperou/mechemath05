% Computes the rigid body transformation that maps a set of points in the 
% box-frame to their corresponding world-frame coordinates
%INPUTS:
% x: the x position of the centroid of the box
% y: the y position of the centroid of the box
% theta: the orientation of the box
% p_list_box: a 2 x n matrix of points in the box frame
%OUTPUTS:
% p_list_world: a 2 x n matrix of points describing the world-frame 
%              coordinates of the points in plist_box
function p_list_world = compute_rbt(x,y,theta,p_list_box)
    rotation = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    p_list_world = rotation * p_list_box' + [x; y];
end