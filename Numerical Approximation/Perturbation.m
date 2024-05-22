function V = Perturbation(perturbations, size)
    % Input: Perturbations specified as {{type, Gamma_0, fun, site, p},...}
    % Output: A perturbation function V in form of a matrix

    % Initialize a final matrix (the sum of all perturbations)
    V = @(t) zeros(size);
    % Loop through all perturbations given
    for i = 1:length(perturbations)
        % Get the general information about this perturbation
        type = perturbations{i}{1};
        Gamma_0 = perturbations{i}{2};
        fun = perturbations{i}{3};

        % Initialize a matrix for this specific perturbation
        F_i = zeros(size);
        
        % Check which type of perturbation we're dealing with
        if type == "Metal"
            % If metal, find the affected site and p
            site = perturbations{i}{4};
            p = perturbations{i}{5};

            % Define Gamma as a 2x2 matrix for the specified p
            Gamma = -1i*Gamma_0/2*(eye(2,2)+[p(3) p(1)-1i*p(2); p(1)+1i*p(2) -p(3)]);
            
            % Find the corresponding matrix position for the site
            matrix_pos = 2*site-1;
            
            % Add elements of Gamma to the matrix
            F_i(matrix_pos, matrix_pos) = Gamma(1, 1); 
            F_i(matrix_pos, matrix_pos+1) = Gamma(1, 2); 
            F_i(matrix_pos+1, matrix_pos) = Gamma(2, 1); 
            F_i(matrix_pos+1, matrix_pos+1) = Gamma(2, 2);

        elseif type == "E-field"
            % Loop through the molecule and add values on diagonal
            for j = 1:size/2
                % Find the corresponding matrix position for the site
                matrix_pos = 2*j-1;

                % Add diagonal elements to the matrix
                F_i(matrix_pos, matrix_pos) =  (j-1)/(size/2-1);
                F_i(matrix_pos+1, matrix_pos+1) = (j-1)/(size/2-1);
                
                % Add off-diagonal elements to the matrix
                if matrix_pos < size-1
                    F_i(matrix_pos+2, matrix_pos) = 1/6;
                    F_i(matrix_pos+3, matrix_pos+1) = 1/6;
                    F_i(matrix_pos, matrix_pos+2) = -1/6;
                    F_i(matrix_pos+1,matrix_pos+3) = -1/6;
                end
            end
            % Scale with Gamma_0
            F_i = Gamma_0*F_i;
        end
        
        % Add the perturbation matrix to the total matrix
        V = @(t) V(t) + F_i*fun(t);
    end
end