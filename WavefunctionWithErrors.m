function [wavefunction_approximation, total_errors] = WavefunctionWithErrors(n, t, H_0, V, start_guess)
    % Get the size of the matrices
    matrix_size = length(H_0);
    
    % Get the number of time steps
    N = length(t);
    
    
    %% Error estimation
    % Create an array to store integration error for each estimation
    integration_errors = cell(1,n+1);
    integration_errors{1} = cell(1,N);
    
    
    %% Determine psi_0 using our start guess
    % Create array for storing each wavevector for each time step
    wavefunctions = cell(1,n+1);
    wavefunctions{1} = cell(1,N);
    
    % Normalize start guess
    normalized_start_guess = start_guess/norm(start_guess);
    
    % Loop through time and add start guess for each time step
    for k = 1:N
        % Initialize the vector of psi_0
        wavefunctions{1}{k} = normalized_start_guess;
        
        % No integral calculation for psi_0, therefore no integration error
        integration_errors{1}{k} = zeros(matrix_size,1);
    end
    
    
    %% Define function for calculating integrand vector for each time
    function [integrand_values, error_integrand_values] = integrand(n)
        % Define integrand functions as exp(i*H0*t)V(t)exp(-i*H0*t)*f(t)
        integrand_function = @(t,k) (-1i)*expmv(1i*H_0,V(t(k))*expmv(-1i*H_0,wavefunctions{n}{k},t(k)-t(1)),t(k)-t(1));
        error_integrand_function = @(t,k) (-1i)*expmv(1i*H_0,V(t(k))*expmv(-1i*H_0,integration_errors{n}{k},t(k)-t(1)),t(k)-t(1));

        % Initialize arrays to store all integrands for each element in vector
        integrand_values = cell(1, matrix_size);
        error_integrand_values = cell(1, matrix_size);
        
        % Loop through time and calculate integrand vector
        for k = 1:N
            % Evaluate the integrand vectors at the current time
            integrand_result = integrand_function(t,k);
            error_integrand_result = error_integrand_function(t,k);

            % Loop through the integrand vectors
            for element = 1:matrix_size
                % Store integrand in its position and time k in the matrix
                integrand_values{1,element}(k) = integrand_result(element, 1);
                error_integrand_values{1,element}(k) = error_integrand_result(element, 1);
            end
        end
    end


    %% Define function for calculating the wavefunction for n > 0
    function [newWavefunction,error] = wavefunction(previousWavefunction)
        % Initialize arrays to store vectors for different times
        newWavefunction = cell(1, N);
        error = cell(1, N);
        
        % The integral of f(t) from t_0 to t_0 is zero, so we initialize as
        newWavefunction{1} = zeros(matrix_size,1);
        error{1} = zeros(matrix_size,1);
        
        % Loop through time and determine the vectors
        for k = 2:N
            % Loop through the molecule and do the integration for each site
            for element = 1:matrix_size
                % Integrate using the Simpson's rule with error estimation
                [newWavefunction{k}(element,1), error{k}(element,1)] = SimpsonWithErrors(t(1:k),previousWavefunction{1,element}(1:k));
            end
        end
    end


    %% Depending on the order n, return an approximation of a wavevector for each time
    % Loop through values of n and determine the wavefunction vector and errors
    for i = 1:n
        % Determine integrand vectors (for wavefunction and error)
        [integrandVector,errorVector] = integrand(i);

        % Calculate the new wavefunction and all 3 parts of integration error
        [wavefunctions{i+1},error_wavefunction] = wavefunction(integrandVector);
        [error_integrated,error_error] = wavefunction(errorVector);

        % Sum up the total integration error
        for k = 1:N
            integration_errors{i+1}{k} = abs(error_wavefunction{k})+abs(error_integrated{k})+abs(error_error{k});
        end
    end
    
    % Create an array to store the final wavevector and error for each time
    wavefunction_approximation = cell(1,N);
    total_errors = cell(1,N);
    
    % Loop through time and sum all order contributions to the wavevector and errors
    for k = 1:N
        wavefunction_approximation{k} = zeros(matrix_size,1);
        total_errors{k} = zeros(matrix_size,1);
        for i = 1:n+1
            wavefunction_approximation{k} = wavefunction_approximation{k}+wavefunctions{i}{k};
            total_errors{k} = total_errors{k}+integration_errors{i}{k};
        end
        wavefunction_approximation{k} = expmv(-1i*H_0,wavefunction_approximation{k},t(k)-t(1));
    end
end


function [I, error] = SimpsonWithErrors(t, f)
    % Computes the integral (I) and estimation of maximal error (error) 
    % Input: Time vector (t) and function value vector (f)
    
    %% Integral calculation
    % Determine the number of time steps
    N = length(t);
    
    % Determine the length of the time step
    h = t(2)-t(1);
    
    % Reshape the function value-vector
    f = reshape(f,1,N);
    
    % Calculate the Simpson specific weights for calculation of the integral
    Weights = 2*ones(N,1); Weights(2:2:end-1) = 4;
    Weights(1) = 1; Weights(N) = 1;
    
    % Calculate the integral
    I = h/3 * f * Weights;
    
    %% Error estimation
    % Define function for approximating 4th derivative (for error estimation)
    function d4 = diff4(f, h)
        % Define function for approximating 1st derivative
        function d1 = diff1(y, h)
            n=length(y);
            d1(1)=(-3*y(1)+4*y(2)-y(3))/(2*h);
            d1(2:n-1)=(y(3:n)-y(1:n-2))/(2*h);
            d1(n)=(3*y(n)-4*y(n-1)+y(n-2))/(2*h);
        end
        % In each differentiation smoothen the data to eliminate noise
        d1 = diff1(f,h);
        d1_f = smoothdata(d1);
        d2 = diff1(d1_f,h);
        d2_f = smoothdata(d2);
        d3 = diff1(d2_f,h);
        d3_f = smoothdata(d3);
        d4 = diff1(d3_f,h);
    end
    
    % Needs at least 3 function points to approximate error
    if N > 2
        % Calculate the 4th derivate of the given function values
        d4 = diff4(f,h);
        % Max estimation of discretization error (see formula sheet)
        error = -(t(end)-t(1))/180*h^4*max(abs(d4));
    else
        error = 0;
    end
end
