%% Parameter values
M = 10;                % Number of laps
N = 4;                % Number of atoms/sites per lap
a = 5*10^(-10);       % Radius of helix
c = M*30*10^(-10);    % Length of helix

epsilon_0 = 3;        % 1st energy term
gamma = 1;            % 2nd energy term
lambda = 10^(-1);     % 3rd energy term

Gamma_01 = 1.8;        % Perturbation term
Gamma_02 = 0.03;        % Perturbation term


%% Construction of Hamiltonian
H_01 = Hamiltonian(N, M, a, c, epsilon_0, gamma, lambda, '+');
H_02 = Hamiltonian(N, M, a, c, epsilon_0, gamma, lambda, '-');


%% Construction of perturbation
fun1 = @(t) 1;
fun2 = @(t) -1;
fun3 = @(t) sin(2*pi/20*t);
fun4 = @(t) -exp(-0.5*t);
fun5 = @(t) -2*(heaviside(t-10)-0.5);

perturbations = {{'Metal', Gamma_01, fun2, 20, [0 0 0]}};
V = Perturbation(perturbations, 2*N*M);


%% Create time vector
t_0 = 0; T = 250; dt = 0.2;
t = t_0:dt:T;


%% Make a starting guess
sites = 1:N*M;
start_guess = ones(2*N*M,1);


%% Determine wavefunctions and relevant outputs
wavefunctions1 = Wavefunction(0,t,H_01,V,start_guess);
wavefunctions2 = Wavefunction(0,t,H_02,V,start_guess);
[n1, m1] = Distributions(wavefunctions1, t);
[n2, m2] = Distributions(wavefunctions2, t);


%% Plot semi-3D colorplots
%ColorPlot(sites, t, n1{1}, n1{2}, '\textbf{Probability Density}', 'Spin', 'Site Index', 'Time', 'Alternative', 'linear')
%ColorPlot(sites, t, m1{3}, m2{3}, '\textbf{Spin Polarization}', 'Helicity', 'Site Index', 'Time', 'Polarized', 'linear')


%% Test for convergence
%ConvergenceTest(1000,t,H_0,V,start_guess)
