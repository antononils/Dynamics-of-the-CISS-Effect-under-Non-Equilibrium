function [charge_distribution, magnetic_distribution] = Distributions(wavefunctions, t)
    % Computes the charge distribution in terms of spin up and spin down
    % and magnetic distribution in components of x, y and z
    % Input: Array of wavevectors for each time (wavefunctions)
    % and time vector (t)
    
    % Determine the number of time steps
    N = length(t);
    
    % Determine the length of the time step
    h = t(2)-t(1);
    
    % Get the size of the matrices
    matrix_size = length(wavefunctions{1});
    
    % Initialize arrays for magnetic and charge distribution
    magnetic_distribution = cell(1,3);
    charge_distribution = cell(1,2);
    
    % Loop through time and split up spin up and spin down parts
    for k = 1:N
        % Loop through each site of the molecule
        for site = 1:matrix_size/2
            % Split up the wavevectors in spin up and down parts
            psi_up = wavefunctions{k}(2*site-1);
            psi_down = wavefunctions{k}(2*site);
            
            % Calculate the magnetic distribution
            magnetic_distribution{1}(k,site) = psi_up'*psi_down + psi_down'*psi_up;
            magnetic_distribution{2}(k,site) = (-1i)*psi_up'*psi_down + (1i)*psi_down'*psi_up;
            magnetic_distribution{3}(k,site) = psi_up'*psi_up - psi_down'*psi_down;
            
            % Calculate the charge distribution
            charge_distribution{1}(k,site) = psi_up'*psi_up;
            charge_distribution{2}(k,site) = psi_down'*psi_down;
        end
    end
end


