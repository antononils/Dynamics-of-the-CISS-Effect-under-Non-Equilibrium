% Function to check for convergence
function ConvergenceTest(n_max,t,H_0,V,start_guess)
    % Initialize psi_n with the unperturbed wavefunction
    psi_1 = Wavefunction(0,t,H_0,V,start_guess);

    % Loop through all n's
    for i = 1:n_max
        % Create psi_n+1
        psi_2 = Wavefunction(i,t,H_0,V,start_guess);
        disp(i)

        % Display the difference between the two orders of perturbation
        disp(norm(psi_2{end}-psi_1{end}))

        % Update psi_n
        psi_1 = Wavefunction(i,t,H_0,V,start_guess);
    end
end