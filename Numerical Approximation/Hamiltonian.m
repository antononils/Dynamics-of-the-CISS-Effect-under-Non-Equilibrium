function H_0 = Hamiltonian(N, M, a, c, epsilon_0, gamma, lambda, helicity)
    % Computes the (sparse) CISS Hamiltonian (H_0)
    % Input: N (Number of atoms/sites per lap), M (Number of laps), a 
    % (Radius of helix), c (Length of helix), epsilon_0 (1st energy term), 
    % gamma (2nd energy term), lambda (3rd energy term) and helicity
    
    % Fill main diagonal with the 1st energy terms
    H_0 = spdiags(epsilon_0*ones(2*N*M),0,2*N*M,2*N*M);
    
    % Fill the two (+-2) diagonals with the 2nd energy terms
    H_0 = spdiags(gamma*ones(2*N*M,1),-2,H_0);
    H_0 = spdiags(gamma*ones(2*N*M,1),2,H_0);
    
    % Define vector v in SOC part of Hamiltonian to describe helicity
    function v = SpinOrbitCoupling(m,sign,a,c,N,M)
        % Define the unit vector for every site pointing from m+s to m
        function d_hat = d_ms(m,s,a,c,N,M)
            % Define the components of the cylindrical coordinates
            function r = r_m(m,a,c,N,M)
                % Define the angle phi in the cylindrical coordinates
                function phi_m = phi_m(m, N)
                    % Check for helicity (orientation of molecule)
                    if helicity == '+'
                        phi_m = (m-1)*2*pi/N;
                    else
                        phi_m = -(m-1)*2*pi/N;
                    end
                end
                r(1) = a*cos(phi_m(m,N));
                r(2) = a*sin(phi_m(m,N));
                r(3) = (m-1)*c/(M*N-1);
            end
            d = (r_m(m,a,c,N,M)-r_m(m+s,a,c,N,M));
            d_hat = d/norm(d);
        end
        
        % Check if s=1 or s=-1 for this part of Hamiltonian
        % Return the vector v in SOC part of Hamiltonian
        if sign == '+'
            v = cross(d_ms(m,1,a,c,N,M), d_ms(m+1,1,a,c,N,M));
        else
            v = cross(d_ms(m,-1,a,c,N,M), d_ms(m-1,-1,a,c,N,M));
        end
    end

    % Initialize the diagonals of the SOC part (see separate notes)
    v_z_plus = zeros(2*N*M,1);
    v_z_minus = zeros(2*N*M,1);
    v_xy_plus = zeros(2*N*M,1);
    v_xy_plus_T = zeros(2*N*M,1);
    v_xy_minus = zeros(2*N*M,1);
    v_xy_minus_T = zeros(2*N*M,1);

    % Initialize indexes to keep track of position on diagonal
    i_plus = 1;
    i_minus = 1;
    
    % Loop through all sites in the molecule and add SOC part of H_0
    % (See separate notes for formulas for each diagonal)
    for m = 1:N*M
        % Add upper diagonal elements to the Hamiltonian
        if m < N*M-1
            v_plus = SpinOrbitCoupling(m,'+',a,c,N,M);
            v_z_plus(i_plus+4) = 1i*lambda*v_plus(3);
            v_z_plus(i_plus+5) = -1i*lambda*v_plus(3);
            v_xy_plus(i_plus+5) = 1i*lambda*(v_plus(1)-1i*v_plus(2));
            v_xy_plus_T(i_plus+4) = 1i*lambda*(v_plus(1)+1i*v_plus(2));
            i_plus = i_plus+2;
        end
        % Add lower diagonal elements to the Hamiltonian
        if m > 2
            v_minus = SpinOrbitCoupling(m,'-',a,c,N,M);
            v_z_minus(i_minus) = 1i*lambda*v_minus(3);
            v_z_minus(i_minus+1) = -1i*lambda*v_minus(3);
            v_xy_minus(i_minus+1) = 1i*lambda*(v_minus(1)-1i*v_minus(2));
            v_xy_minus_T(i_minus) = 1i*lambda*(v_minus(1)+1i*v_minus(2));
            i_minus = i_minus+2;
        end
    end
    
    % Add the remaining diagonals corresponding to the SOC part of H_0
    H_0 = spdiags(v_z_minus,-4,H_0);
    H_0 = spdiags(v_z_plus,4,H_0);
    H_0 = spdiags(v_xy_minus,-3,H_0);
    H_0 = spdiags(v_xy_minus_T,-5,H_0);
    H_0 = spdiags(v_xy_plus,5,H_0);
    H_0 = spdiags(v_xy_plus_T,3,H_0);
    
end

