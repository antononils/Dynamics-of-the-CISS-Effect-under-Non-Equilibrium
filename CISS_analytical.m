syms t 

y = 1;
eps = 3;
a = 5e-10;
c = 30e-10;
lambda = 1e-3;
Gamma = 50e-3;
T = 10;
time_step = 50;
p = [0 0 0];
chrility = true;
perturbation = [];



% f = -1i*50*cos(10*t); 
func = 0;

H_pos = Perturbation_analytic("pos",a,c,lambda, Gamma, chrility, perturbation,p,func);
H_neg = Perturbation_analytic("neg", a, c, lambda, Gamma, chrility, perturbation, p, func);
H_neg
[J_pos_up, J_pos_down, t_vec] = Probablity_density(H_pos, eps, y, T, time_step);
[J_neg_up, J_neg_down, t_vec]= Probablity_density(H_neg, eps, y, T, time_step);


Xdata = [1 2 3 4];
ColorPlot(Xdata, t_vec, J_pos_up-J_pos_down, J_neg_up-J_neg_down, "Probability density", "Spin", "Site", "Time","Alternative")





