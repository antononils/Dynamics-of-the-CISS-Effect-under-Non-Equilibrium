syms t 

y = 1;
eps = 3;
a = 5e-10;
c = 30e-10;
lambda = 1e-3;
Gamma = 50e-3;
T = 10;
time_step = 5;
p = [0 0 0];
chrility = true;
perturbation = None;



% f = -1i*50*cos(10*t); 
func = 0;

H_pos = Perturbation_analytic("pos",a,c,lambda,chrility,perturbation,p,func);
H_neg = Perturbation_analytic("neg",a,c,lambda,chrility,perturbation,p,func);

J_pos = Probablity_density(H_pos, T, time_step);
J_neg = Probablity_density(H_neg, T, time_step);

Xdata = [1 2 3 4];
ColorPlot(Xdata, J_neg(3), J_pos(1)-J_pos(2), J_neg(1)-J_neg(2), "Probability density", "Spin", "Site", "Time","Alternative")





