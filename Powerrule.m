function productEquation = Powerrule(t,count)
    if t <= 0
        productEquation = sym(0);
        return;
    else
    rval9 = 2;
    while rval9 == 2
         disp(' ')
         %fprintf('You selected Algebraic function.\n');
         fprintf('Now entering variables for term %d.\n',count);
         fprintf('Example formatting: \n')
         fprintf('Ex: ax^n \n')
         fprintf('Ex: g(x) = x^3, a = 1, n = 3 \n')
         a = input('Input your coefficient (a): \n');

         %Error Check
         if isnumeric(a)
         else
             fprintf('Invalid input. Please enter a numeric value for the coefficient.\n');
             a = input('Input your coefficient (a): \n');
         end

         n = input('Input your exponent (n): \n');

         %Error Check
         if isnumeric(n)
         else
            fprintf('Invalid input. Please enter a numeric value for the exponent.\n');
            n = input('Input your exponent (n): \n');
         end

                    
         %Error Check
         corinput = input('Do these values match what you input (1 for yes 2 for no) \n');
         if corinput == 1
             rval9 = 1;
         elseif corinput == 2
         end
    end
    if n == 1
        productEquation = a + Powerrule(t-1, count+1);
        return;
    end
    if n == 0
        productEquation = sym(0) + Powerrule(t-1,count+1);
        return;
    end
    if n < 0
        coefficient = a * n; 
        n = abs(n);
        syms x;
        productEquation = coefficient / x^(n+1) + Powerrule(t-1,count+1);
        return;
    end
    syms x; 
    coefficient = a * n;
    derivative = coefficient * x^(n-1); 
    productEquation = derivative + Powerrule(t-1,count+1);
    end
end

