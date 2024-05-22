syms t 


%% Defining the parameters
y = 1;
eps = 3;
a = 5e-10;
c = 30e-10;
lambda = 1e-2;
Gamma = 1.8/2;

%% Setting the perturbation
p = [-1 0 0];
chrility = true;
perturbation = "metal_l";
func = 1-exp(-t);

%%  Setting the Hamiltonian for the positiv resp. negativ helicity
H_pos = Perturbation_analytic("pos",a,c,lambda, Gamma, chrility, perturbation,p,func);
H_neg = Perturbation_analytic("neg", a, c, lambda, Gamma, chrility, perturbation, p, func);

%% Calculating the probability denisity for the positive resp. negative helicity
T = 10;
time_step = 100;
[J_pos_up, J_pos_down, t_vec] = Probability_density(H_pos, eps, y, T, time_step);
[J_neg_up, J_neg_down, t_vec]= Probability_density(H_neg, eps, y, T, time_step);

%% Plotting the probability density
Sites = [1 2 3 4];
ColorPlot(Sites, t_vec, J_pos_up-J_pos_down, J_neg_up-J_neg_down, "\textbf{Spin Polarization}", "Helicity", "Site Index", "Time","Polarized", "linear")





