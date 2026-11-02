%This function computes the value of X at the next time step
%for any arbitrary embedded RK method
%INPUTS:
%name: a string representing the name of the Butcher Tableau or associated 
%      method
%OUTPUTS:
%tableau: a struct containing the A, B, and C values for the requested
%         Butcher Tableau
function [tableau, p] = init_BT(name)
    tableau = struct();

    % fixed-step (non-embedded/non adaptive) BTs:
    if name == "Fourth Order"
        tableau.A = [0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0];
        tableau.B = [1/6 1/3 1/3 1/6];
        tableau.C = [0; 0.5; 0.5; 1];
        p = 4;
    elseif name == "Forward Euler"
        tableau.A = 0;
        tableau.B = 1;
        tableau.C = 0;
        p = 1;
    elseif name == "Explicit Midpoint"
        tableau.A = [0 0; 0.5 0];
        tableau.B = [0 1];
        tableau.C = [0; 0.5];
        p = 2;
    elseif name == "Ralston's"
        tableau.A = [0 0; 2/3 0];
        tableau.B = [1/4 3/4];
        tableau.C = [0; 2/3];
        p = 2;
    elseif name == "Kutta's Third-Order"
        tableau.A = [0 0 0; 1/3 0 0; 0 2/3 0];
        tableau.B = [1/4 0 3/4];
        tableau.C = [0; 1/3; 2/3];
        p = 3;
    elseif name == "Van der Houwen's Third-Order"
        tableau.A = [0 0 0; 8/15 0 0; 1/4 5/12 0];
        tableau.B = [1/4 0 3/4];
        tableau.C = [0; 8/15; 2/3];
        p = 3;        
    % embedded methods:
    elseif name == "Dormand Prince"
        tableau.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
        tableau.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
        5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
        tableau.A = [0,0,0,0,0,0,0;
                1/5, 0, 0, 0,0,0,0;...
                3/40, 9/40, 0, 0, 0, 0,0;...
                44/45, -56/15, 32/9, 0, 0, 0,0;...
                19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
                9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
                35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];
        p = 4;
    elseif name == "Fehlberg"
        tableau.C = [0, 1/4, 3/8, 12/13, 1, 1/2];
        tableau.B = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55;...
            25/216, 0, 1408/2565, 2197/4104, -1/5, 0];
        tableau.A = [0,0,0,0,0,0;...
                    1/4, 0,0,0,0,0;...
                    3/32, 9/32, 0,0,0,0;...
                    1932/2197, -7200/2197, 7296/2197, 0,0,0;...
                    439/216, -8, 3680/513, -845/4104, 0,0;...
                    -8/27, 2, -3544/2565, 1859/4104, -11/40, 0];
        p = 4;
    elseif name == "Heun Euler"
        tableau.C = [0,1];
        tableau.B = [1/2,1/2;1,0];
        tableau.A = [0,0;1,0];
        p = 1;
    elseif name == "Fehlberg RK1"
        tableau.C = [0,1/2,1];
        tableau.B = [1/512, 255/256, 1/512; 1/256, 255/256, 0];
        tableau.A = [0,0,0;1/2,0,0;1/256,255/256,0];
        p = 1;
    elseif name == "Bogacki"
        tableau.C = [0,1/2, 3/4, 1];
        tableau.B = [2/9, 1/3, 4/9, 0; 7/24, 1/4, 1/3, 1/8];
        tableau.A = [0,0,0,0; 1/2,0,0,0; 0,3/4,0,0; 2/9,1/3, 4/9, 0];
        p = 2;
    else
        disp("Not a valid method name.")
    end
end