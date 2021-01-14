classdef SurfMat
    %SURFMAT Uses Molecular Dynamics output files of a glass slab
    %subjected to a corrosive environment on its top and bottom surfaces
    %to count the topological constraints
    %Author: Karan Doss (2020)
    
    methods (Static)
        function [label, coordinates, scale] = generate_coordinates(filename)
            %GENERATE_COORDINATES Generate coordinates from .xyz MD Files
            %   Generate coordinates from .xyz Dynamics files that have
            %   normalized (x,y,z) coordinates.
            %   # Argument: filename (string) - file_path/file_name.xyz
            %   # Output: label (nx1 double) - labels for atoms
            %   # Output: coordinates (nx3 double) - x, y, z position
            %   column vectors
            %   # Output: scale (3x1 double) - scaling factor used to scale
            %   normalized coordinates
            %   scale = [ub_x - lb_x, ub_y - lb_y, ub_z - lb_z]
            %   where ub_i and lb_i are respectively the upper bound (or
            %   highest value) and lower bound (or lowest value) values
            %   along coordinate i (= x,y,z)
            
            readtext = textread(filename,'%s');
            lower_bounds = str2double(readtext(15:2:20));   % Lower bounds for x,y,z
            upper_bounds = str2double(readtext(16:2:20));   % Upper bounds for x,y,z
            scale = upper_bounds - lower_bounds;    % Scaling factor for normalized coordinates x,y,z
            
            raw = str2double(readtext(28:end)); % raw data read from text
            pd = 5; % period = number of colums for data
            
            % Format into output matrix
            for i = 1:5;
                counts = 1;
                for j = i:pd:length(raw)
                    res(counts,i) = raw(j);
                    counts = counts + 1;
                end
            end
            
            label = res(:,2);   % Atom labels
            coordinates = [res(:,3)*scale(1), res(:,4)*scale(2), res(:,5)*scale(3)]; %x,y,z corrdinates in Angstrom
        end
        
        function [new_labels, new_coordinates] = delete_modifiers(labels, coordinates, modifier_labels)
            %DELETE_MODIFIERS Delete modifier atom coordinates
            %   Delete coordinates of atoms that have modifier labels since
            %   they do not contribute to the rigidity of the network.
            %   # Argument: labels (nx1 double) - labels for atoms
            %   # Argument: coordinates (nx3 double) - x, y, z position
            %   column vectors
            %   # Argument: modifier_labels (mx1 double) - Numbers that
            %   correspond to labels of modifiers.
            %   # Output: new_labels (n'x1 double) - labels vector with
            %   modifier labels deleted
            %   # Output: new_coordinates (n'x3 double) - coordinates of
            %   non-modifier atoms
            
            counts = 1;
            for i = 1:length(labels)
                if ismember(labels(i),modifier_labels) == false
                    new_labels(counts,1) = labels(i);
                    new_coordinates(counts,:) = coordinates(i,:);
                    counts = counts +1;
                end
            end
        end
        
        function [new_labels, new_coordinates, coordination] = ...
                get_coordination(labels,coordinates,cutoff_bond_radius,scale, network_former_labels)
            %GET_COORDINATION Delete NBOs & calculate coordination of atoms
            %   Use a cutoff bond radius to delete all non-bridging oxygens
            %   (NBOs) and evaluate the coordination # of each network
            %   forming atom.
            %   # Argument: labels (n'x1 double) - labels for atoms
            %   # Argument: coordinates (n'x3 double) - x, y, z position
            %   column vectors
            %   # Argument: cutoff_bond_radius (1x1 double) - cutoff bond
            %   radius for an M-O bond in Angstroms.
            %   # Argument: scale (1x1 double) - scaling factor used to
            %   scale normalized coordinates
            %   # Argument: network_former_labels (kx1 double) - Numbers
            %   that correspond to labels of network formers.
            %   # Output: new_labels (n"x1 double) - labels vector with
            %   non-bridging oxygen labels deleted
            %   # Output: new_coordinates (n"x3 double) - coordinates of
            %   bridging oxygen atoms
            %   # Output: coordination (n"x3 double) - Coordination of
            %   network forming atoms.
            
            noc = size(coordinates,1);
            
            % Coordinates for applying Periodic Boundary Conditions for x-y coorinates
            m1 = coordinates - scale(1)*[ones(noc,1), zeros(noc,2)];    % C(x,y,z) - (Lx,0,0)
            m2 = coordinates - scale(2)*[zeros(noc,1), ones(noc,1), zeros(noc,1)];  % C(x,y,z) - (0,Ly,0)
            m3 = coordinates + scale(1)*[ones(noc,1), zeros(noc,2)];    % C(x,y,z) + (Lx,0,0)
            m4 = coordinates + scale(2)*[zeros(noc,1), ones(noc,1), zeros(noc,1)];  % C(x,y,z) + (0,Ly,0)
            pseudo_coordinates = [coordinates; m1; m2; m3; m4];
            
            j = 1;  % Initialize index
            
            for row = 1 : length(coordinates)
                distances = sqrt((coordinates(row,1) - pseudo_coordinates(:, 1)).^2 + ...
                    (coordinates(row, 2) - pseudo_coordinates(:, 2)).^2 + ...
                    (coordinates(row, 3) - pseudo_coordinates(:, 3)).^2);
                
                % Count the number of points that are greater than 0 (to
                % exclude this point itself) but less than or equal to
                % the bond radius.
                
                COUNTS(row) = sum(distances > 0 & distances <= cutoff_bond_radius);
                
                if COUNTS(row) >= 0 && ismember(labels(row), network_former_labels(2:end)) == true
                    res(j,1) = labels(row);
                    res(j,2:4) = coordinates(row,:);
                    res(j,5) = COUNTS(row);
                    j = j+1;
                elseif COUNTS(row) == 2 && ismember(labels(row), network_former_labels(1)) == true
                    res(j,1) = labels(row);
                    res(j,2:4) = coordinates(row,:);
                    res(j,5) = COUNTS(row);
                    j = j+1;
                end
            end
            new_labels = res(:,1);
            new_coordinates = res(:,2:4);
            coordination = res(:,5);
        end
        
        function [top_surface_constraints, bottom_surface_constraints, cs_x, cs_y] = ...
                count_surface_constraints(coordinates, coordination, scale, surface_depth, ngrid)
            %COUNT_SURFACE_CONSTRAINTS Count surface constraints
            %   Partition the top and bottom surfaces into their respective
            %   grids and calculate the average number of constraints per
            %   pixel within the grid.
            %   # Argument: coordinates (n"x3 double) - x, y, z position
            %   column vectors
            %   # Argument: coordination (n"x3 double) - Coordination of
            %   network forming atoms.
            %   # Argument: scale (3x1 double)
            %   scale = [ub_x - lb_x, ub_y - lb_y, ub_z - lb_z]
            %   where ub_i and lb_i are respectively the upper bound (or
            %   highest value) and lower bound (or lowest value) values
            %   along coordinate i (= x,y,z)
            %   # Argument: surface_depth (1x1 double) Depth (in angstroms)
            %   below the surface upto which constraints are counted.
            %   Default value is 6 Angstroms.
            %   # Argument: surface_depth (1x1 double) Depth (in angstroms)
            %   below the surface upto which constraints are counted.
            %   # Output: MAT1
            %   # Output: MAT1_1
            %   # Output: cs_x - Equispaced grid from 0 to scale(1),
            %   default number of points = 10. (along x)
            %   # Output: cs_y - Equispaced grid from 0 to scale(1),
            %   default number of points = 10. (along y)
            
            arguments
                coordinates double
                coordination double
                scale double
                surface_depth (1,1) double = 6  % Angstroms
                ngrid (1,1) double = 10
            end
            
            grid_size_x = scale(1)/ngrid;   % x - Grid Size
            grid_size_y = scale(2)/ngrid;   % y - Grid Size
            grid_size = [grid_size_x, grid_size_y];
            cs_x = linspace(0,scale(1),ngrid);  % Scaled x-axis for contour plot
            cs_y = linspace(0,scale(2),ngrid);  % Scaled y-axis for contour plot
            mat = zeros(ngrid);
            avg = zeros(ngrid);
            mat_1 = zeros(ngrid);
            avg_1 = zeros(ngrid);
            
            % Bounds for cutting the surfaces
            top_lower_limit = max(coordinates(:,3));
            top_upper_limit = max(coordinates(:,3)) - surface_depth;
            bottom_lower_limit = min(coordinates(:,3)) + surface_depth;
            bottom_upper_limit = min(coordinates(:,3));
            
            for x = 1:ngrid
                for y = 1:ngrid
                    for i = 1:length(coordinates)
                        if coordinates(i,3) < top_lower_limit && coordinates(i,3) > top_upper_limit %Top Half
                            if coordinates(i,1:2) < grid_size.*[x y] & coordinates(i,1:2) > (grid_size).*[x-1 y-1]
                                mat(x,y) = mat(x,y) + coordination(i);
                                avg(x,y) = avg(x,y) + 1;
                            end
                        elseif coordinates(i,3) < bottom_lower_limit && coordinates(i,3) > bottom_upper_limit %Bottom Half
                            if coordinates(i,1:2) < grid_size.*[x y] & coordinates(i,1:2) > (grid_size).*[x-1 y-1]
                                mat_1(x,y) = mat_1(x,y) + coordination(i);
                                avg_1(x,y) = avg_1(x,y) + 1;
                            end
                        end
                    end
                end
            end
            
            top_surface_coord(1:ngrid, 1:ngrid) = mat./avg; % Average coordination per atom for each voxel
            top_surface_coord(isnan(top_surface_coord))=0;
            top_surface_coord(isinf(top_surface_coord))=0;
            top_surface_coord(top_surface_coord<0)=0;
            
            top_surface_constraints(1:ngrid, 1:ngrid)  = 0.5*top_surface_coord(1:ngrid, 1:ngrid) + ...
                (2*top_surface_coord(1:ngrid, 1:ngrid) - 3);    % Average number of constraints per atom for each pixel
            top_surface_constraints(isnan(top_surface_constraints))=0;
            top_surface_constraints(isinf(top_surface_constraints))=0;
            top_surface_constraints(top_surface_constraints<0)=0;
            
            bottom_surface_coord(1:ngrid, 1:ngrid) = mat_1./avg_1;
            bottom_surface_coord(isnan(bottom_surface_coord))=0;
            bottom_surface_coord(isinf(bottom_surface_coord))=0;
            bottom_surface_coord(bottom_surface_coord<0)=0;
            
            bottom_surface_constraints(1:ngrid, 1:ngrid)  = 0.5*bottom_surface_coord(1:ngrid, 1:ngrid) + ...
                (2*bottom_surface_coord(1:ngrid, 1:ngrid) - 3);
            bottom_surface_constraints(isnan(bottom_surface_constraints))=0;
            bottom_surface_constraints(isinf(bottom_surface_constraints))=0;
            bottom_surface_constraints(bottom_surface_constraints<0)=0;
        end
        
        function [counts, phase_fraction, avg_phase_fraction] = analyze_surface_constraints(top_surface, bottom_surface, lb, ub)
            %ANALYZE_SURFACE_CONSTRAINTS Analyze fraction of constraints
            %   Analyze the fraction of surface constraints that lie in 
            %   between lower bound (lb, default - 2.8) and upper bound
            %   (ub, default = 3.5)
            %   SThis functionality is still in beta testing!
            arguments
                top_surface double
                bottom_surface double
                lb (1,1) double = 2.8
                ub (1,1) double = 3.5
            end
            
            n = size(top_surface);
            surf = [reshape(top_surface, [n(1)*n(2),n(3)]); reshape(bottom_surface, [n(1)*n(2),n(3)])];
            for i = 1:n(3)
                N(i,:) = histcounts(surf(:,i), [0:0.5:7], 'normalization', 'pdf');
            end
            counts = [[0.5:0.5:7]', N'];
            phase_fraction = sum(surf >= lb & surf <= ub)./length(surf >= lb & surf <= ub);
            avg_phase_fraction = mean(phase_fraction)
        end
    end
end