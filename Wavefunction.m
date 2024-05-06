function wavefunction_approximation = Wavefunction(n, t, H_0, V, start_guess)
    % Get the size of the matrices
    matrix_size = length(H_0);
    
    % Get the number of time steps
    N = length(t);
    
    
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
    end
    
    
    %% Define function for calculating integrand vector for each time
    function integrand_values = integrand(n)
        % Define integrand functions as exp(i*H0*t)V(t)exp(-i*H0*t)*psi(t)
        integrand_function = @(t,k) (-1i)*expmv(1i*H_0,V(t(k))*expmv(-1i*H_0,wavefunctions{n}{k},t(k)-t(1)),t(k)-t(1));
            
        % Initialize an array to store all integrands for each element in vector
        integrand_values = cell(1, matrix_size);
        
        % Loop through time and calculate integrand vector
        for k = 1:N
            % Evaluate the integrand vector at the current time
            integrand_result = integrand_function(t,k);

            % Loop through the integrand vector
            for element = 1:matrix_size
                % Store integrand in its position and time k in the matrix
                integrand_values{1,element}(k) = integrand_result(element, 1);
            end
        end
    end


    %% Define function for calculating the wavefunction for n > 0
    function newWavefunction = wavefunction(previousWavefunction)
        % Initialize an array to store wavefunction vectors for different times
        newWavefunction = cell(1, N);
        
        % The integral of f(t) from t_0 to t_0 is zero, so we initialize as
        newWavefunction{1} = zeros(matrix_size,1);
        
        % Loop through time and determine the wavefunction vector
        for k = 2:N
            % Loop through the molecule and do the integration for each site
            for element = 1:matrix_size
                % Integrate using Simpson's rule for integration
                newWavefunction{k}(element,1) = Simpson(t(1:k),previousWavefunction{1,element}(1:k));
            end
        end
    end


    %% Depending on the order n, return an approximation of a wavevector for each time
    % Loop through values of n and determine the wavefunction vectors
    for i = 1:n
        i
        integrandVector = integrand(i);
        wavefunctions{i+1} = wavefunction(integrandVector);
    end
    
    % Create an array to store the final wavevector for each time
    wavefunction_approximation = cell(1,N);
    
    % Create a vector to store norms of wavefunction for each time
    norms = zeros(N,1);

    % Loop through time and sum all order contributions to the wavevector
    for k = 1:N
        wavefunction_approximation{k} = zeros(matrix_size,1);
        for i = 1:n+1
            wavefunction_approximation{k} = wavefunction_approximation{k}+wavefunctions{i}{k};
        end
        wavefunction_approximation{k} = expmv(-1i*H_0,wavefunction_approximation{k},t(k)-t(1));
        
        % Calculate the norm of the current wavefunction
        norms(k) = norm(wavefunction_approximation{k})^2;
    end

    % Determine the time-average of the norms
    avg_norm = (1/(t(end)-t(1)))*Simpson(t,norms);

    % Loop through time and "normalize" the wavefunctions
    for k = 1:N
        wavefunction_approximation{k} = wavefunction_approximation{k}/avg_norm;
    end
end


function result = Simpson(t, f)
    % Computes the integral (I)
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
    result = h/3 * f * Weights;
end

