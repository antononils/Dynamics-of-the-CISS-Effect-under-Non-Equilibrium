%% Parameter values
M = 5;               % Number of laps
N = 4;                % Number of atoms/sites per lap
a = 5*10^(-10);       % Radius of helix
c = M*30*10^(-10);    % Length of helix

epsilon_0 = 3;        % 1st energy term
gamma = 1;            % 2nd energy term
lambda = 10^(-2);     % 3rd energy term

Gamma_0 = 1.8;        % Perturbation term


%% Construction of Hamiltonian
H_01 = Hamiltonian(N, M, a, c, epsilon_0, gamma, lambda, '+');
H_02 = Hamiltonian(N, M, a, c, epsilon_0, gamma, lambda, '-');


%% Construction of perturbation
% fun = @(t) 1;
% V = Perturbation(true,"1",true,'2', Gamma_0, fun, 2*N*M);
V = @(t) Gamma_0*cos(t)*F;


%% Create time vector
t_0 = 0; T = 20; dt = 0.4;
t = t_0:dt:T;


%% Make a starting guess
sites = 1:N*M;
start_guess = ones(2*N*M,1);


%% Determine wavefunctions and relevant outputs
wavefunctions1 = Wavefunction(100,t,H_01,V,start_guess);
wavefunctions2 = Wavefunction(100,t,H_02,V,start_guess);
[n1, m1] = Distributions(wavefunctions1, t);
[n2, m2] = Distributions(wavefunctions2, t);


%% Plot semi-3D colorplots
%ColorPlot(sites, t, n1{1}, n1{2}, '\textbf{Probability Density}', 'Spin', 'Site Index', 'Time', 'Alternative', 'linear')
ColorPlot(sites, t, m1{3}, m2{3}, '\textbf{Spin Polarization}', 'Helicity', 'Site Index', 'Time', 'Polarized', 'linear')


%% Test for convergence
%ConvergenceTest(1000,t,H_0,V,start_guess)






































%% Plot 2D graphs with time development for specific sites for choosen output
% for i = [4,3,2,1]
%     y1 = charge_distribution{1}(:,i);
%     y2 = charge_distribution{2}(:,i);
%     figure(i)
%     plot(t,y1,t,y2)
%     legend('Spin up','Spin down')
% end





%% Examine the integration errors and plot as function of time step
% time_steps = 0.01:0.005:1;
% normsOfErrors = zeros(length(time_steps),1);
% for i = 1:length(time_steps)
%     i
%     t = t_0:time_steps(i):T;
%     [wavefunctions,error] = WavefunctionWithErrorsMIT(1,t,H_0,V,start_guess);
%     for k = 1:length(t)
%         normsOfErrors(i) = normsOfErrors(i)+norm(error{k});
%     end
%     normsOfErrors(i) = normsOfErrors(i)/length(t);
% end
% plot(time_steps,normsOfErrors)





















% f = fun(t);
% ft = fft(f);
% fs = 1/dt;
% f = (0:length(ft)-1)*fs/length(ft);
% plot(f,abs(ft))

function d1 = diff1(y, h)
    n = length(y);
    d1(1) = (-3*y(1)+4*y(2)-y(3))/(2*h);
    d1(2:n-1) = (y(3:n)-y(1:n-2))/(2*h);
    d1(n) = (3*y(n)-4*y(n-1)+y(n-2))/(2*h);
end