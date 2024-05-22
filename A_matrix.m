syms s f(x) g(x)
y = 1; 

A = [0 y 0 0 0 0 0 0 ;
     y 0 y 0 0 0 0 0 ;
     0 y 0 y 0 0 0 0 ;
     0 0 y 0 y 0 0 0 ;
     0 0 0 y 0 y 0 0 ;
     0 0 0 0 y 0 y 0 ;
     0 0 0 0 0 y 0 y ;
     0 0 0 0 0 0 y 0 ;
    
     ];

det(A-s*eye(8,8))

f(x) = x^4 - 7*x^3 + 15*x^2 - 10*x +1;

solve(f(x) == 0,x)
g(x) = x^3-6*x^2+9*x-1;
solve(g(x)==0,x)