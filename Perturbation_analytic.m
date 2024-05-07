function H = Perturbation_analytic(Helicity, radius, length,lambda, chirility, perturbation,p,func)
    a = radius;
    c = length;
    d= 2*a*c*(1/3);
    b = 2*a^2;
    
    if chirility == false
        g = 0;
    elseif chirility == true
        g = 1/(b+(c^2)/9);
    end
    
    if Helicity == "pos"
        h_13 = g*[1i*b 1i*d;
           1i*d -1i*b];
        h_24 = g*[1i*b d;
          -d  -1i*b];
        h_31 = g*[-1i*b -1i*d;
          -1i*d  1i*b];
        h_42 = g*[-1i*b -d;
            d  1i*b];
            
        H = [h_11 zeros(2,2) h_13 zeros(2,2);
            zeros(2,2) zeros(2,2) zeros(2,2) h_24;
            h_31 zeros(2,2) zeros(2,2) zeros(2,2);
            zeros(2,2) h_42 zeros(2,2) zeros(2,2)];
    elseif Helicity == "neg"
        h_13 = g*[-1i*b -1i*d;
           -1i*d 1i*b];
        h_24 = g*[-1i*b d;
          -d  1i*b];
        h_31 = g*[1i*b 1i*d;
          1i*d  -1i*b];
        h_42 = g*[1i*b -d;
            d  -1i*b];
            
        H = lambda*[h_11 zeros(2,2) h_13 zeros(2,2);
            zeros(2,2) zeros(2,2) zeros(2,2) h_24;
            h_31 zeros(2,2) zeros(2,2) zeros(2,2);
            zeros(2,2) h_42 zeros(2,2) zeros(2,2)];
    end

    V = zeros(8,8);

    if  ismember("metal_r", perturbation)
        v = Gammma*(eye(2,2)+[p(3) p(1)-1i*p(2);
            p(1)+1i*p(2) -p(3)]);
        V = V + [zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2);
            zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2);
            zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2);
            zeros(2,2) zeros(2,2) zeros(2,2) v];
    end

    if  ismember("metal_l", perturbation)
        v = Gammma*(eye(2,2)+[p(3) p(1)-1i*p(2);
            p(1)+1i*p(2) -p(3)]);
        V = V + [v zeros(2,2) zeros(2,2) zeros(2,2);
            zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2);
            zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2);
            zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2)];
    end

    H = H + func*V;
end
