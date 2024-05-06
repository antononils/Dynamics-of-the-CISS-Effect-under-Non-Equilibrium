function V = Perturbation(type, fun, size)
    if type == "metal"
        F = zeros(size); F(1,1) = 1; F(2,2) = 1;
    elseif type == "field"
        F = eye(size);
    end
    V = @(t) F*fun(t);
end