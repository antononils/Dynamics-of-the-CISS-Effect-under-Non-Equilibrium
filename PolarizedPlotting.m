min_value = 1;
step = 1;
max_value = 10;

%values = min_value:step:max_value;
values = logspace(-1,4,5);

% Create time vector
t_0 = 0; T = 20; dt = 0.2;
t = t_0:dt:T;

spinUp = zeros(length(t),length(values));
spinDown = zeros(length(t),length(values));


length(values)
for i = 1:length(values)
    i
    x = values(i);

    % Parameter values
    M = 5;               % Number of laps
    N = 4;                % Number of atoms/sites per lap
    a = 5*10^(-9);          % Radius of helix
    c = M*30*10^(-9);      % Length of helix
    
    epsilon_0 = 3;        % 1st energy term
    gamma = 1;            % 2nd energy term
    lambda = 10^(-2);           % 3rd energy term
    
    Gamma = 1.8;   % Perturbation term

    % Construction of Hamiltonian
    H_0 = Hamiltonian(N, M, a, c, epsilon_0, gamma, lambda, '+');

    % Construction of perturbation
    fun = @(t) cos(x*t);
    perturbations = {{'E-field', Gamma_0, fun2}};
    V = Perturbation(perturbations, 2*N*M);

    % Make a starting guess
    sites = 1:N*M;
    start_guess = ones(2*N*M,1);
    
    % Determine wavefunctions and relevant outputs
    wavefunctions = Wavefunction(100,t,H_0,V,start_guess);
    [n, m] = Distributions(wavefunctions, t);

    % Polarize the data
    spinUp(:,i) = Matrix2Vector(n{1});
    spinDown(:,i) = Matrix2Vector(n{2});
end

% Normalize the data
spinUp = spinUp./max(abs(spinUp),[],'all');
spinDown = spinDown./max(abs(spinDown),[],'all');

% Plot semi-3D colorplots
ColorPlot(values,t,spinUp,spinDown,'Polarized Probability Density', 'Spin', 'Frequency', 'Time', 'Polarized', 'log')



function vector = Matrix2Vector(matrix)
    weights = linspace(-length(matrix(1,:)),length(matrix(1,:)),length(matrix(1,:)));
    vector = zeros(length(matrix(:,1)),1);
    for k = 1:length(matrix(:,1))
        vector(k) = sum(weights.*matrix(k,:));
    end
end