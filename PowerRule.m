t = input('Enter the value of t (this value is the number of terms in your problem): \n');

function [productEquation] = Powerrule(t)
    if t <= 0
        productEquation = sym(0);
        return ;
    end 
    a = input('Enter the value of a: ');
    n = input('Enter the value of n: '); 
    if n == 1
        productEquation = a + Powerrule(t-1);
        return;
    end
    if n == 0
        productEquation = sym(0) + Powerrule(t-1);
        return;
    end
    if n < 0
        coefficient = a * n; 
        n = abs(n);
        syms x;
        productEquation = coefficient / x^(n+1) + Powerrule(t-1);
        return;
    end
    syms x; 
    coefficient = a * n;
    derivative = coefficient * x^(n-1); 
    productEquation = derivative + Powerrule(t-1);
end
disp(Powerrule(t))