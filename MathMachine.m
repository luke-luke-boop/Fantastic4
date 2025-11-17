clear
clc
syms x y

% ==============================
%     The Math Machine
%   Derivative Solver (with Log)
% ==============================

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
logfile = sprintf('MathMachine_Log_%s.txt', timestamp);
diary(logfile);  % Start recording session output
disp('==============================')
disp('      The Math Machine')
disp('     Derivative Solver')
disp('==============================')
disp(['Session Log: ', logfile])
disp('')

while true
    % Main Menu
    disp('Select a differentiation method:');
    disp('1. Power Rule');
    disp('2. Product Rule');
    disp('3. Chain Rule');
    disp('4. Quotient Rule');
    disp('5. Implicit Differentiation');
    disp('6. Logarithmic Differentiation');
    disp('7. L''Hopital''s Rule');
    disp('8. Polynomial Differentiation');
    disp('9. Exit');
    choice = input('Enter your choice (1-9): ');

    switch choice
        case 1
            PowerruleInteractive();
        case 2
            productRule();
        case 3
            chainRule();
        case 4
            quotientRule();
        case 5
            implicitDiff();
        case 6
            LogRuleInteractive();
        case 7
            lHopital();
        case 8
            polyDiff();
        case 9
            disp('Exiting The Math Machine. Goodbye!');
            break;
        otherwise
            disp('Invalid selection. Please enter a number between 1 and 9.');
    end

    disp('------------------------------------------');
    next = input('Would you like to return to the main menu? (y/n): ', 's');
    if lower(next) ~= 'y'
        disp('Exiting The Math Machine. Goodbye!');
        break;
    end
    clc
end

diary off;  % Stop recording session

% -----------------------------
% FUNCTION DEFINITIONS BELOW
% -----------------------------

function chainRule()
    clc;
    syms x

    disp('You selected the Chain Rule.');
    disp('------------------------------------------');
    disp('This rule applies when one function is inside another: y = f(g(x))');
    disp('Formula: dy/dx = f''(g(x)) * g''(x)');
    disp('Example: y = sin(x^2) → dy/dx = cos(x^2) * 2*x');
    disp('------------------------------------------');

    func = input('Enter your composite function (e.g., sin(x^2), exp(3*x), (x^2+1)^5): ', 's');

    % Determine outer function
    if contains(func, 'sin(')
        outer = 'sin(u)';
        outer_deriv = 'cos(u)';
    elseif contains(func, 'cos(')
        outer = 'cos(u)';
        outer_deriv = '-sin(u)';
    elseif contains(func, 'tan(')
        outer = 'tan(u)';
        outer_deriv = 'sec(u)^2';
    elseif contains(func, 'exp(')
        outer = 'e^u';
        outer_deriv = 'e^u';
    elseif contains(func, 'log(')
        outer = 'ln(u)';
        outer_deriv = '1/u';
    elseif contains(func, '^')
        outer = 'u^n';
        outer_deriv = 'n*u^(n-1)';
    else
        disp('Function type not recognized. Try using sin, cos, exp, log, or ^.');
        return;
    end

    % Extract inner function inside first parentheses
    startIdx = strfind(func, '(');
    endIdx = strfind(func, ')');

    if ~isempty(startIdx) && ~isempty(endIdx)
        firstClose = endIdx(find(endIdx > startIdx(1), 1));
        inner = func(startIdx(1)+1:firstClose-1);
    else
        disp('Could not detect an inner function automatically. Please enter it manually.');
        inner = input('Enter the inner function g(x): ', 's');
    end

    disp(' ');
    disp('--- Step-by-Step Solution ---');
    fprintf('Step 1: Original function: y = %s\n', func);
    fprintf('Step 2: Outer function f(u) = %s\n', outer);
    fprintf('Step 3: Inner function g(x) = %s\n', inner);
    fprintf('Step 4: Derivative of outer function f''(u) = %s\n', outer_deriv);

    % Calculate derivative of inner function manually
    dydx_inner = PowerruleInteractiveForTerm(inner);
    fprintf('Step 5: Derivative of inner function g''(x) = %s\n', dydx_inner);

    % Apply chain rule
    dy_dx = strrep(outer_deriv, 'u', ['(' inner ')']);
    fprintf('Step 6: Apply chain rule: dy/dx = f''(g(x)) * g''(x) = %s * %s\n', dy_dx, dydx_inner);
    disp('------------------------------------------');
end

function result = PowerruleInteractiveForTerm(inner)
    % Calls the full interactive Power Rule for one term
    disp('Finding derivative of inner function using Power Rule...');
    result = Powerrule(1,1); % simplified: one term
end

function PowerruleInteractive()
    t = input('How many terms are in your polynomial? ');
    derivative = Powerrule(t,1);
    fprintf('Derivative of your polynomial: %s\n', char(derivative));
end

function productRule()
    clc;
    disp('You selected Product Rule.');
    disp('------------------------------------------');
    disp('This section will be implemented later.');
    disp('------------------------------------------');
end

function quotientRule()
    clc;
    disp('You selected Quotient Rule.');
    disp('------------------------------------------');
    disp('This section will be implemented later.');
    disp('------------------------------------------');
end

function implicitDiff()
    clc;
    disp('You selected Implicit Differentiation.');
    disp('------------------------------------------');
    disp('This section will be implemented later.');
    disp('------------------------------------------');
end

function LogRuleInteractive()
    clc;
    disp('You selected Logarithmic Differentiation.');
    disp('------------------------------------------');
    terms = input('How many logarithmic terms? ');
    derivative = LogRule(terms);
    fprintf('Derivative of your logarithmic function: %s\n', char(derivative));
end

function lHopital()
    clc;
    syms x

    disp('You selected L''Hopital''s Rule.');
    disp('------------------------------------------');
    disp('L''Hopital''s Rule is used for limits of the form 0/0 or ∞/∞.');
    disp('If lim f(x) and lim g(x) give 0/0 or ∞/∞, then');
    disp('    lim f(x)/g(x) = lim f''(x)/g''(x),  (if the latter limit exists)');
    disp('------------------------------------------');
    disp('Enter your limit:  lim_{x -> a}  f(x) / g(x)');
    disp(' ');

    % Get user input as strings
    num_str = input('Enter the numerator f(x) (e.g., sin(x)-x): ', 's');
    den_str = input('Enter the denominator g(x) (e.g., x^3): ', 's');
    a_str   = input('Enter the point a for x -> a (e.g., 0, 2, pi, inf, -inf): ', 's');

    % Convert to symbolic expressions
    try
        f = str2sym(num_str);
        g = str2sym(den_str);
        a = str2sym(a_str);
    catch
        disp('Error converting input to symbolic expressions. Please check your syntax.');
        return;
    end

    disp(' ');
    disp('--- Step-by-Step L''Hopital Evaluation ---');
    fprintf('We are evaluating:  lim_{x -> %s} (%s) / (%s)\n\n', char(a), char(f), char(g));

    % Compute the original limit
    originalLimit = limit(f/g, x, a);
    fprintf('Step 1: Try the limit directly (without L''Hopital):\n');
    fprintf('    lim f(x)/g(x) = %s\n', char(originalLimit));
    disp(' ');

    % Evaluate numerator and denominator limits separately
    num_val = limit(f, x, a);
    den_val = limit(g, x, a);

    fprintf('Step 2: Check the form of the limit:\n');
    fprintf('    lim f(x) -> %s\n', char(num_val));
    fprintf('    lim g(x) -> %s\n', char(den_val));

    isZeroNum = isequal(simplify(num_val), sym(0));
    isZeroDen = isequal(simplify(den_val), sym(0));

    % For infinities, convert to double and test isinf (inside try/catch)
    isInfNum = false;
    isInfDen = false;
    try
        isInfNum = isinf(double(num_val));
        isInfDen = isinf(double(den_val));
    catch
        % If conversion fails, leave as false
    end

    if ~( (isZeroNum && isZeroDen) || (isInfNum && isInfDen) )
        disp(' ');
        disp('The limit is NOT of the form 0/0 or ∞/∞, so L''Hopital''s Rule does not apply.');
        fprintf('Final answer: lim_{x -> %s} f(x)/g(x) = %s\n', char(a), char(originalLimit));
        disp('------------------------------------------');
        return;
    end

    disp(' ');
    disp('The limit is an indeterminate form (0/0 or ∞/∞).');
    disp('We can apply L''Hopital''s Rule: differentiate numerator and denominator.');
    disp(' ');

    % Iterative L'Hopital
    maxIter = 10;
    f_curr = f;
    g_curr = g;

    for k = 1:maxIter
        fprintf('--- Application %d of L''Hopital''s Rule ---\n', k);

        % Differentiate numerator and denominator
        f_curr = diff(f_curr, x);
        g_curr = diff(g_curr, x);

        fprintf('New numerator f_%d(x) = %s\n', k, char(f_curr));
        fprintf('New denominator g_%d(x) = %s\n', k, char(g_curr));

        % Compute limits of new numerator and denominator
        num_val = limit(f_curr, x, a);
        den_val = limit(g_curr, x, a);

        fprintf('lim_{x -> %s} f_%d(x) = %s\n', char(a), k, char(num_val));
        fprintf('lim_{x -> %s} g_%d(x) = %s\n', char(a), k, char(den_val));

        isZeroNum = isequal(simplify(num_val), sym(0));
        isZeroDen = isequal(simplify(den_val), sym(0));

        isInfNum = false;
        isInfDen = false;
        try
            isInfNum = isinf(double(num_val));
            isInfDen = isinf(double(den_val));
        catch
        end

        % If no longer indeterminate, we can get the limit
        if ~( (isZeroNum && isZeroDen) || (isInfNum && isInfDen) )
            disp(' ');
            disp('The resulting limit is no longer indeterminate.');
            newLimit = num_val / den_val;
            newLimit = simplify(newLimit);

            fprintf('So, lim_{x -> %s} f(x)/g(x) = lim_{x -> %s} f_%d(x)/g_%d(x) = %s\n', ...
                char(a), char(a), k, k, char(newLimit));
            disp('------------------------------------------');
            return;
        end

        disp('Still an indeterminate form; applying L''Hopital''s Rule again...');
        disp(' ');
    end

    % If we reach here, we never resolved the form
    disp('Reached the maximum number of L''Hopital applications without resolving the limit.');
    fprintf('Symbolic direct limit (if it exists) was: %s\n', char(originalLimit));
    disp('------------------------------------------');
end

function polyDiff()
    clc;
    disp('You selected Polynomial Differentiation.');
    disp('------------------------------------------');
    disp('This section will be implemented later.');
    disp('------------------------------------------');
end

% -----------------------------
% Power Rule Function (your groupmate)
% -----------------------------
function productEquation = Powerrule(t,count)
    if t <= 0
        productEquation = sym(0);
        return;
    else
    rval9 = 2;
    while rval9 == 2
         disp(' ')
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


function productEquation = LogRule(term)
    if term <= 0
        productEquation = sym(0);
        return;
    end
    numerator = sym(0);
    denominator = sym(0);
    syms x
    answer = input('Is your term natural log (1) or log base number (2): \n');
    if answer == 1 % ln scenario
        c = input('How many terms in f(x)? ');
        numerator = 0;
        denominator = 0;
        while c > 0
            b = input('Format of this term? (powerRule, base):','s');
            if strcmp(b,'powerRule')
                result2 = Powerrule(1,1);
                numerator = numerator + result2;
                denominator = denominator + 1; % simple placeholder
            elseif strcmp(b,'base')
                numerator = 1;
                denominator = 1;
            end
            c = c-1;
        end
        productEquation = numerator/denominator + L
    end
end
