syms t 

y = 1;
eps = 3;
a = 5e-10;
c = 30e-10;
lambda = 1e-3;
Gamma = 50e-3;

p = [0 0 0];
chrility = true;
perturbation = "metal_l";
func = 0;


H_pos = Perturbation_analytic("pos",a,c,lambda, Gamma, chrility, perturbation,p,func);
H_neg = Perturbation_analytic("neg", a, c, lambda, Gamma, chrility, perturbation, p, func);

T = 10;
time_step = 100;
[J_pos_up, J_pos_down, t_vec] = Probablity_density(H_pos, eps, y, T, time_step);
[J_neg_up, J_neg_down, t_vec]= Probablity_density(H_neg, eps, y, T, time_step);



Sites = [1 2 3 4];
ColorPlot(Sites, t_vec, J_pos_up-J_pos_down, J_neg_up-J_neg_down, "\textbf{Spin Polarization}", "Helicity", "Site Index", "Time","Polarized", "linear")





