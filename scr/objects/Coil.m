classdef Coil
    properties
        % Geometry
        R                 % Core radius [m]
        wire_d            % Wire diameter [m]
        l                 % Length of winding [m]
        N                 % Turns per layer
        layers            % Number of winding layers
        points_per_turn   % Points per turn for wire resolution

        % Positioning
        position          % [x, y, z] - coil center in space
        rotation          % [ux, uy, uz, angle] - rotation axis + angle (degrees)

        % Current
        I                 % Current through the coil [A]

        % Computed geometry
        points       % [Nx3] coordinates of the wire path
    end

    methods
        function obj = Coil(R, wire_d, l, N, layers, points_per_turn, position, rotation, I)
            obj.R = R;
            obj.wire_d = wire_d;
            obj.l = l;
            obj.N = N;
            obj.layers = layers;
            obj.points_per_turn = points_per_turn;
            obj.position = position;
            obj.rotation = rotation;
            obj.I = I;
        end
        
        % Generate points
        function obj = generate(obj)
            total_points = obj.points_per_turn * obj.N * obj.layers;
            obj.points = zeros(total_points, 3);
            idx = 1;
        
            for k = 1:obj.layers
                % Set layer radius
                Rk = obj.R + (k-1)*obj.wire_d;
                
                % Defining the position of the finite element
                t = linspace(0, 2*pi*obj.N, obj.points_per_turn*obj.N);
                x = Rk * cos(t);
                y = Rk * sin(t);
                
                % Determine Z position
                if mod(k,2) == 1
                    z = linspace(0, obj.l, length(t)) - obj.l/2;
                else
                    z = linspace(obj.l, 0, length(t)) - obj.l/2;
                end
                

                num_points = length(t);
                obj.points(idx:idx+num_points-1, :) = [x', y', z'];
                idx = idx + num_points;
            end

            axis = obj.rotation(1:3);
            angle = deg2rad(obj.rotation(4)); % convert degrees to radians
        
            if norm(axis) > 0
                axis = axis / norm(axis); % normalize axis
                % Rodrigues' rotation formula
                K = [    0       -axis(3)  axis(2);
                      axis(3)     0      -axis(1);
                     -axis(2)  axis(1)     0    ];
                R = eye(3) + sin(angle)*K + (1 - cos(angle))*(K*K);
            else
                R = eye(3); % no rotation
            end
        
            obj.points = (R * obj.points')'; % apply rotation
        
            obj.points = obj.points + obj.position;
        end

        function plot(obj)
            plot3(obj.points(:,1),...
                  obj.points(:,2),...
                  obj.points(:,3), 'r', 'LineWidth', 0.001);
        end

        function obj = setPosition(obj, position, rotation)
            obj.position = position
            obj.rotation = rotation
            axis = obj.rotation(1:3);
            angle = deg2rad(obj.rotation(4)); % convert degrees to radians
        
            if norm(axis) > 0
                axis = axis / norm(axis); % normalize axis
                % Rodrigues' rotation formula
                K = [    0       -axis(3)  axis(2);
                      axis(3)     0      -axis(1);
                     -axis(2)  axis(1)     0    ];
                R = eye(3) + sin(angle)*K + (1 - cos(angle))*(K*K);
            else
                R = eye(3); % no rotation
            end
        
            obj.points = (R * obj.points')'; % apply rotation
        
            obj.points = obj.points + obj.position;
        end

        function [Bx, By, Bz] = calculateField(obj, obs_points)
            mu0 = 4*pi*1e-7;
            Bx = zeros(size(obs_points,1),1);
            By = zeros(size(obs_points,1),1);
            Bz = zeros(size(obs_points,1),1);
        
            R_min = obj.R - obj.wire_d/2;
            R_max = obj.R + (obj.layers-1)*obj.wire_d + obj.wire_d/2;
        
            axis = obj.rotation(1:3);
            angle = deg2rad(obj.rotation(4));
            if norm(axis) > 0
                axis = axis / norm(axis);
                K = [0, -axis(3), axis(2);
                     axis(3), 0, -axis(1);
                    -axis(2), axis(1), 0];
                R = eye(3) + sin(-angle)*K + (1 - cos(-angle))*(K*K);
            else
                R = eye(3);
            end
        
            hWait = waitbar(0, 'Calculating field from coil...');
            nPoints = size(obs_points,1);
        
            for p = 1:nPoints
                obs_global = obs_points(p,:);
        
                obs_local = (R * (obs_global - obj.position)')';
        
                r_xy = norm(obs_local(1:2));
                z_val = obs_local(3);
        
                inside_radial = (r_xy >= R_min) && (r_xy <= R_max);
                inside_height = (z_val >= -obj.l/2) && (z_val <= obj.l/2);
        
                if inside_radial && inside_height
                    Bx(p) = NaN;
                    By(p) = NaN;
                    Bz(p) = NaN;
                    continue;
                end
        
                B = [0, 0, 0];
                for i = 1:(size(obj.points,1)-1)
                    r1 = obj.points(i,:);
                    r2 = obj.points(i+1,:);
                    dl = r2 - r1;
                    r = obs_global - (r1 + r2)/2;
                    r_norm = norm(r);
                    if r_norm < 1e-9
                        continue;
                    end
                    dB = (mu0*obj.I)/(4*pi) * cross(dl, r) / (r_norm^3);
                    B = B + dB;
                end
        
                Bx(p) = B(1);
                By(p) = B(2);
                Bz(p) = B(3);
        
                if mod(p, 100) == 0 || p == nPoints
                    waitbar(p / nPoints, hWait);
                end
            end    
            close(hWait);
        end
    end
end
