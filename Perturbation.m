function V = Perturbation(type, Gamma_0, fun, size)
    if type == "1"
        p = [0 0 1];
    elseif type == "2"
        p = [0 0 -1];
    elseif type == "3"
        p = [0 1 0];
    elseif type == "4"
        p = [0 -1 0];
    elseif type == "5"
        p = [1 0 0];
    elseif type == "6"
        p = [-1 0 0];
    elseif type == "Normal"
        p = [0 0 0];
    end
    Gamma = -1i*Gamma_0/2*(eye(2,2)+[p(3) p(1)-1i*p(2); p(1)+1i*p(2) -p(3)]);
    F = zeros(size); F(1,1) = Gamma(1,1); F(1,2) = Gamma(1,2); F(2,1) = Gamma(2,1); F(2,2) = Gamma(2,2);
    V = @(t) F*fun(t);
end