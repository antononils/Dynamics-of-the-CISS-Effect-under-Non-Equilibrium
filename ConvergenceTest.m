function ConvergenceTest(n_max,t,H_0,V,start_guess)
    psi_1 = Wavefunction(0,t,H_0,V,start_guess);
    for i = 1:n_max
        psi_2 = Wavefunction(i,t,H_0,V,start_guess);
        disp(i)
        disp(norm(psi_2{end}-psi_1{end}))
        psi_1 = Wavefunction(i,t,H_0,V,start_guess);
    end
end