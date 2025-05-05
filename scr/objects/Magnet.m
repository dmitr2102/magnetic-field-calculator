classdef Magnet
    properties
        % Geometry
        R          % Radius of the cylindrical magnet [m]
        l          % Length of the magnet [m]
        Ms         % Saturation magnetization [A/m]
        n_radial   % Number of radial divisions
        n_angular  % Number of angular divisions
        n_axial    % Number of axial divisions

        % Positioning
        position   % [x, y, z] position of the magnet center
        rotation   % [ux, uy, uz, angle] — rotation axis and angle [degrees]

        % Internal representation
        points     % [Nx3] coordinates of magnet elements
        dm         % [Nx3] magnetic moments of each element
    end

    methods
        function obj = Magnet(R, l, Ms, n_radial, n_angular, n_axial, position, rotation)
            obj.R = R;
            obj.l = l;
            obj.Ms = Ms;
            obj.n_radial = n_radial;
            obj.n_angular = n_angular;
            obj.n_axial = n_axial;
            obj.position = position;
            obj.rotation = rotation;
        end

        function obj = generate(obj)
            % Discretization step sizes
            dr = obj.R / obj.n_radial;
            dtheta = 2 * pi / obj.n_angular;
            dz = obj.l / obj.n_axial;

            % Generate cylindrical grid
            r = linspace(dr/2, obj.R - dr/2, obj.n_radial);
            theta = linspace(0, 2*pi, obj.n_angular);
            z = linspace(-obj.l/2 + dz/2, obj.l/2 - dz/2, obj.n_axial);
            [Rgrid, Theta, Z] = ndgrid(r, theta, z);

            % Convert to Cartesian coordinates
            X = Rgrid .* cos(Theta);
            Y = Rgrid .* sin(Theta);
            obj.points = [X(:), Y(:), Z(:)];

            % Element volume and magnetization
            dV = Rgrid .* dr .* dtheta .* dz;
            dm = zeros(size(obj.points));
            dm(:,3) = obj.Ms * dV(:); % Initial magnetization along Z

            % Apply rotation
            axis = obj.rotation(1:3);
            angle = deg2rad(obj.rotation(4));
            if norm(axis) > 0
                axis = axis / norm(axis);
                K = [0 -axis(3) axis(2); 
                    axis(3) 0 -axis(1); 
                    -axis(2) axis(1) 0];
                Rmat = eye(3) + sin(angle)*K + (1 - cos(angle))*(K*K);
                obj.points = (Rmat * obj.points')';
                dm = (Rmat * dm')';
            end

            % Apply translation
            obj.points = obj.points + obj.position;
            obj.dm = dm;
        end

        function plot(obj)
            plot3(obj.points(:,1),...
                  obj.points(:,2),...
                  obj.points(:,3), 'r', 'LineWidth', 0.001);
        end

        function [Bx, By, Bz] = calculateField(obj, obs_points)
            mu0 = 4*pi*1e-7;
            N = size(obs_points,1);
            Bx = zeros(N,1);
            By = zeros(N,1);
            Bz = zeros(N,1);
        
            R_max = obj.R;
            Z_half = obj.l / 2;
        
            axis = obj.rotation(1:3);
            angle = deg2rad(obj.rotation(4));
            if norm(axis) > 0
                axis = axis / norm(axis);
                K = [0 -axis(3) axis(2);
                     axis(3) 0 -axis(1);
                    -axis(2) axis(1) 0];
                Rinv = eye(3) + sin(-angle)*K + (1 - cos(-angle))*(K*K);
            else
                Rinv = eye(3);
            end
        
            hWait = waitbar(0, 'Calculating field from magnet...');
        
            for i = 1:N
                r_obs_global = obs_points(i,:);
        
                r_obs_local = (Rinv * (r_obs_global - obj.position)')';
        
                r_xy = norm(r_obs_local(1:2));
                z_val = r_obs_local(3);
        
                inside_radial = (r_xy <= R_max);
                inside_height = (z_val >= -Z_half) && (z_val <= Z_half);
        
                if inside_radial && inside_height
                    Bx(i) = NaN;
                    By(i) = NaN;
                    Bz(i) = NaN;
                    continue;
                end
        
                % Вычисление поля
                B = [0, 0, 0];
        
                for j = 1:size(obj.points,1)
                    r0 = obj.points(j,:);
                    m = obj.dm(j,:);
                    r_vec = r_obs_global - r0;
                    r_norm = norm(r_vec);
        
                    if r_norm < 1e-9
                        continue;
                    end
        
                    term1 = 3 * r_vec * dot(m, r_vec) / (r_norm^5);
                    term2 = m / (r_norm^3);
                    dB = (mu0 / (4*pi)) * (term1 - term2);
                    B = B + dB;
                end
        
                Bx(i) = B(1);
                By(i) = B(2);
                Bz(i) = B(3);
        
                if mod(i,100) == 0 || i == N
                    waitbar(i/N, hWait);
                end
            end
        
            close(hWait);
        end

        function M_total = calculateTorque(obj, xq, yq, zq, Bx, By, Bz)

            M_total = [0, 0, 0];
            for i = 1:size(obj.points,1)
                pos = obj.points(i,:);
        
                Bx_i = interp3(xq, yq, zq, Bx, pos(1), pos(2), pos(3), 'linear', 0);
                By_i = interp3(xq, yq, zq, By, pos(1), pos(2), pos(3), 'linear', 0);
                Bz_i = interp3(xq, yq, zq, Bz, pos(1), pos(2), pos(3), 'linear', 0);
                B_i = [Bx_i, By_i, Bz_i];
        
                dm_i = obj.dm(i,:);
        
                dM = cross(dm_i, B_i);
                
                M_total = M_total + dM;
            end            
         end
    end
end
