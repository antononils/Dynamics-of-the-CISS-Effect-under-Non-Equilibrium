syms t x J_up_1(x) J_up_2(x) J_up_3(x) J_up_4(x) J_down_1(x) J_down_2(x) J_down_3(x) J_down_4(x)

y = 1;
A = [0 y 0 0;
    y 0 y 0;
    0 y 0 y;
    0 0 y 0];


[v,l]  = eig(A);
h=1;

l_1 = diag(exp(diag(1i*t*l/h)));
l_2 = diag(exp(-1*diag(1i*t*l/h)));

V_1 = kron(simplify(v*l_1*inv(v)),eye(2));
V_2 = kron(simplify(v*l_2*inv(v)), eye(2));

% f = -1i*50*cos(10*t); 
f = 0;


eps = 3;
a = 5e-10;
c = 30e-10;
lambda = 1e-3;


Initial_value = transpose([1 1 1 1 1 1 1 1]);
Initial_value = Initial_value/norm(Psi);

H = Perturbation_analytic("pos");



Integrand = simplify(1/(1i*h)*V_1*H_pos*V_2*Psi);
I = int(Integrand, t, 0, x);
I = simplify(I);
U = exp(-1i*x*eps/h)*(kron(v*diag(exp(diag(-1i*l*x/h)))*inv(v),eye(2,2)));

Psi_tot = U*(Initial_value + I);


up_1 = [1 0 0 0 0 0 0 0];
Psi_up_1 = (up_1'.*Psi_tot);

up_2 = [0 0 1 0 0 0 0 0];
Psi_up_2 = (up_2'.*Psi_tot);

up_3 = [0 0 0 0 1 0 0 0];
Psi_up_3 = (up_3'.*Psi_tot);

up_4 = [0 0 0 0 0 0 1 0];
Psi_up_4 = (up_4'.*Psi_tot);

down_1 = [0 1 0 0 0 0 0 0];
Psi_down_1 = (down_1'.*Psi_tot);

down_2 = [0 0 0 1 0 0 0 0];
Psi_down_2 = (down_2'.*Psi_tot);

down_3 = [0 0 0 0 0 1 0 0];
Psi_down_3 = (down_3'.*Psi_tot);

down_4 = [0 0 0 0 0 0 0 1];
Psi_down_4 = (down_4'.*Psi_tot);


J_up_1(x) = real(Psi_up_1'*Psi_up_1);
J_up_2(x) = real(Psi_up_2'*Psi_up_2);
J_up_3(x) = real(Psi_up_3'*Psi_up_3);
J_up_4(x) = real(Psi_up_4'*Psi_up_4);



J_down_1(x) = real(Psi_down_1'*Psi_down_1);
J_down_2(x) = real(Psi_down_2'*Psi_down_2);
J_down_3(x) = real(Psi_down_3'*Psi_down_3);
J_down_4(x) = real(Psi_down_4'*Psi_down_4);

T = 10;
time_step = 100;

t_vec = linspace(0,T,time_step);

J_up = [vpa(J_up_1(t_vec))' vpa(J_up_2(t_vec))' vpa(J_up_3(t_vec))' vpa(J_up_4(t_vec))'];
J_down = [vpa(J_down_1(t_vec))' vpa(J_down_2(t_vec))' vpa(J_down_3(t_vec))' vpa(J_down_4(t_vec))'];

Xdata = [1 2 3 4];

ColorPlot(Xdata, t_vec, J_up-J_down, J_down, "Probability density", "Spin", "Site", "Time","Alternative")





