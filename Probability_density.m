function [J_up, J_down, t_vec] = Probability_density(H, eps, y, T, time_steps)
    %% Defining the variables 
    syms t x J_up_1(x) J_up_2(x) J_up_3(x) J_up_4(x) J_down_1(x) J_down_2(x) J_down_3(x) J_down_4(x)
    %% Setting the initial condition
    Psi = transpose([1 1 1 1 1 1 1 1]);
    Psi = Psi/norm(Psi);

    t_vec = linspace(0,T,time_steps);
    
    
    %% Solving for A eigenvaluse and eigenvectors
    A = [0 y 0 0;
        y 0 y 0;
        0 y 0 y;
        0 0 y 0];


    [v,l]  = eig(A);
    l_1 = diag(exp(diag(1i*t*l)));
    l_2 = diag(exp(-1*diag(1i*t*l)));

    V_1 = kron(simplify(v*l_1*inv(v)),eye(2));
    V_2 = kron(simplify(v*l_2*inv(v)), eye(2));

    %% Solving for the first perturbation term    
    Integrand = simplify(1/(1i)*V_1*H*V_2*Psi);
    I = int(Integrand, t, 0, x);
    I = simplify(I);
    U = exp(-1i*x*eps)*(kron(v*diag(exp(diag(-1i*l*x)))*inv(v),eye(2,2)));
    
    
    Psi_tot = U*(Psi + I);

    %% The spin-wave funcitons for each site
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

    %% Calculatging the probability density
    J_up_1(x) = real(Psi_up_1'*Psi_up_1);
    J_up_2(x) = real(Psi_up_2'*Psi_up_2);
    J_up_3(x) = real(Psi_up_3'*Psi_up_3);
    J_up_4(x) = real(Psi_up_4'*Psi_up_4);

    J_down_1(x) = real(Psi_down_1'*Psi_down_1);
    J_down_2(x) = real(Psi_down_2'*Psi_down_2);
    J_down_3(x) = real(Psi_down_3'*Psi_down_3);
    J_down_4(x) = real(Psi_down_4'*Psi_down_4);

    J_up = [vpa(J_up_1(t_vec))' vpa(J_up_2(t_vec))' vpa(J_up_3(t_vec))' vpa(J_up_4(t_vec))'];
    J_down = [vpa(J_down_1(t_vec))' vpa(J_down_2(t_vec))' vpa(J_down_3(t_vec))' vpa(J_down_4(t_vec))'];
    