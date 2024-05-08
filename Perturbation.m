function V = Perturbation(type, left, pL, right, pR, Gamma_0, fun, size)
    if type == "metal"
        GammaL = -1i*Gamma_0/2*(eye(2,2)+[pL(3) pL(1)-1i*pL(2); pL(1)+1i*pL(2) -pL(3)]);
        GammaR = -1i*Gamma_0/2*(eye(2,2)+[pR(3) pR(1)-1i*pR(2); pR(1)+1i*pR(2) -pR(3)]);
        F = zeros(size); 
        if left
            F(1,1) = GammaL(1,1); F(1,2) = GammaL(1,2); F(2,1) = GammaL(2,1); F(2,2) = GammaL(2,2);
        end
        if right
            F(end-1,end-1) = GammaR(1,1); F(end-1,end) = GammaR(1,2); F(end,end-1) = GammaR(2,1); F(end,end) = GammaR(2,2);
        end
    elseif type == "E-field"
        F = eye(2*M*N);
        for j = 1:M*N
            i = 2*j-1;
            F(i,i) = i-1;
            F(i+1,i+1) = i-1;
        end
        F = F./max(F,[],'all');
    end
    V = @(t) F*fun(t);
end