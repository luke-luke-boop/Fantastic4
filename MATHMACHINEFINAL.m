% MathMachine.m
clear; clc
syms x y

% ==============================
%     The Math Machine
%   Derivative Solver
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
    disp('6. L''Hopital''s Rule');
    disp('7. Exit');
    choice = input('Enter your choice (1-7): ');

    switch choice
        case 1
            powerRule();
        case 2
            productRule();
        case 3
            chainRule();
        case 4
            quotientRule();
        case 5
            implicitDiff();
        case 6
            lHopital();
        case 7
            disp('Exiting The Math Machine. Goodbye!');
            break;
        otherwise
            disp('Invalid selection. Please enter a number between 1 and 7.');
    end
    clc
end

diary off;  % Stop recording session
%% ==============================
% CHAIN RULE FUNCTION WITH STEPS
%% ==============================
function chainRule()
    clc;
    syms x
    disp('You selected the Chain Rule.');
    disp('------------------------------------------');
    disp('Enter a composite function y = f(g(x))');
    disp('Supported: sin, cos, tan, exp, log, powers, sums, products');
    disp('Example: sin(x^2), exp(3*x), (x^2+1)^5, tan(exp(x^3))');
    disp('------------------------------------------');

    while true
        % ----- Get valid function input -----
        while true
            func_input = input('Enter your function in terms of x: ', 's');
            try
                % standardize ln to log
                func_input = strrep(func_input,'ln','log');
                f = str2sym(func_input);
                break;
            catch
                disp('ERROR: Invalid function input. Please check your syntax and try again.');
            end
        end

        % ----- Display original function -----
        disp(' ');
        disp('Original function:');
        fprintf('y = %s\n', char(f));

        % ----- Compute derivative with steps -----
        try
            [dy_dx, steps] = computeChainDerivativeSteps(f, x);
        catch ME
            disp('ERROR while computing the derivative:');
            disp(ME.message);
            continue;  % ask user to reinput
        end

        % ----- Display steps -----
        disp(' ');
        disp('--- Step-by-Step Derivative ---');
        for k = 1:length(steps)
            fprintf('%s\n', steps{k});
        end

        % ----- Display final derivative -----
        disp(' ');
        disp('--- Final Simplified Derivative ---');
        fprintf('dy/dx = %s\n', char(simplify(dy_dx)));
        disp('------------------------------------------');

        % ----- Ask if retry -----
        retry = '';
        while ~ismember(retry, {'y','n'})
            retry = lower(strtrim(input('Would you like to try another Chain Rule problem? (y/n): ', 's')));
            if ~ismember(retry, {'y','n'})
                disp('Invalid input. Please enter y or n.');
            end
        end
        if strcmp(retry,'y')
            clc;
            continue;
        end

        % ----- Ask if return to main menu or exit -----
        menuReturn = '';
        while ~ismember(menuReturn, {'y','n'})
            menuReturn = lower(strtrim(input('Would you like to return to the main menu? (y/n): ', 's')));
            if ~ismember(menuReturn, {'y','n'})
                disp('Invalid input. Please enter y or n.');
            end
        end

        if strcmp(menuReturn,'y')
            clc;
            return;
        else
            disp('Exiting The Math Machine. Goodbye!');
            error('ExitMathMachine');  % terminate function
        end
    end
end

%% ============================================
% Recursive derivative solver WITH STEPS
%% ============================================
function [dydx, steps] = computeChainDerivativeSteps(f, x)
    steps = {};
    
    if isempty(symvar(f))  % constant
        dydx = sym(0);
        steps{end+1} = sprintf('Derivative of constant %s is 0', char(f));
        return;
    end
    if isequal(f, x)  % variable
        dydx = sym(1);
        steps{end+1} = sprintf('Derivative of variable %s is 1', char(f));
        return;
    end

    op = char(feval(symengine,'op',f,0));
    ch = children(f);

    % Standardize operators
    if strcmp(op,'plus'), op = '_plus'; end
    if strcmp(op,'minus'), op = '_minus'; end
    if strcmp(op,'times'), op = '_times'; end
    if strcmp(op,'power'), op = '_power'; end
    if strcmp(op,'_mult'), op = '_times'; end

    switch op
        case '_plus'
            dydx = 0;
            for k=1:length(ch)
                [dterm, s] = computeChainDerivativeSteps(ch{k}, x);
                dydx = dydx + dterm;
                steps = [steps, s];
            end
        case '_minus'
            if length(ch)==1
                [dterm, s] = computeChainDerivativeSteps(ch{1}, x);
                dydx = -dterm;
                steps = [steps, strcat('Derivative of -(', s{:}, ') = -(', char(dterm), ')')];
            else
                [d1, s1] = computeChainDerivativeSteps(ch{1}, x);
                [d2, s2] = computeChainDerivativeSteps(ch{2}, x);
                dydx = d1 - d2;
                steps = [steps, s1, s2, sprintf('Subtract derivatives: %s - %s', char(d1), char(d2))];
            end
        case '_times'
            dydx = 0;
            terms = length(ch);
            for i = 1:terms
                term_prod = 1;
                term_steps = {};
                for j = 1:terms
                    if i==j
                        [dterm, s] = computeChainDerivativeSteps(ch{j}, x);
                        term_prod = term_prod * dterm;
                        term_steps = [term_steps, s];
                    else
                        term_prod = term_prod * ch{j};
                    end
                end
                dydx = dydx + term_prod;
                steps = [steps, term_steps, sprintf('Term %d contribution to product rule: %s', i, char(term_prod))];
            end
        case '_power'
            base = ch{1};
            exponent = ch{2};
            if isempty(symvar(exponent))
                [dBase, sBase] = computeChainDerivativeSteps(base, x);
                dydx = exponent * base^(exponent-1) * dBase;
                steps = [steps, sBase, sprintf('Power rule: d(%s^%s)/dx = %s', char(base), char(exponent), char(dydx))];
            else
                [dBase, sBase] = computeChainDerivativeSteps(base, x);
                [dExp, sExp] = computeChainDerivativeSteps(exponent, x);
                dydx = f*(dExp*log(base) + (exponent/base)*dBase);
                steps = [steps, sBase, sExp, sprintf('Generalized power rule: d(%s)/dx = %s', char(f), char(dydx))];
            end
        case 'sin'
            [u, sU] = computeChainDerivativeSteps(ch{1}, x);
            dydx = cos(ch{1}) * u;
            steps = [steps, sU, sprintf('d(sin(%s))/dx = cos(%s) * d(%s)/dx = %s', char(ch{1}), char(ch{1}), char(ch{1}), char(dydx))];
        case 'cos'
            [u, sU] = computeChainDerivativeSteps(ch{1}, x);
            dydx = -sin(ch{1}) * u;
            steps = [steps, sU, sprintf('d(cos(%s))/dx = -sin(%s) * d(%s)/dx = %s', char(ch{1}), char(ch{1}), char(ch{1}), char(dydx))];
        case 'tan'
            [u, sU] = computeChainDerivativeSteps(ch{1}, x);
            dydx = sec(ch{1})^2 * u;
            steps = [steps, sU, sprintf('d(tan(%s))/dx = sec(%s)^2 * d(%s)/dx = %s', char(ch{1}), char(ch{1}), char(ch{1}), char(dydx))];
        case 'exp'
            [u, sU] = computeChainDerivativeSteps(ch{1}, x);
            dydx = exp(ch{1}) * u;
            steps = [steps, sU, sprintf('d(exp(%s))/dx = exp(%s) * d(%s)/dx = %s', char(ch{1}), char(ch{1}), char(ch{1}), char(dydx))];
        case {'log','ln'}
            [u, sU] = computeChainDerivativeSteps(ch{1}, x);
            dydx = (1/ch{1}) * u;  % chain rule: derivative of ln(f(x)) = f'/f
            steps = [steps, sU, sprintf('d(log(%s))/dx = 1/(%s) * d(%s)/dx = %s', char(ch{1}), char(ch{1}), char(ch{1}), char(dydx))];
        otherwise
            error("Unsupported operation encountered: " + op);
    end
end


%% -----------------------------
% Power Rule (kept from your earlier version)
%% -----------------------------
function powerRule()
    t = input('How many terms are in your polynomial? ');
    derivative = Powerrule(t,1);
    fprintf('Derivative of your polynomial: %s\n', char(derivative));
end

function productEquation = Powerrule(t, count)
    syms x;

    % Base case
    if t <= 0
        productEquation = sym(0);
        return;
    end

    % Interactive input for each term
    fprintf('Entering term %d\n', count);
    a = input('Input coefficient a: ');
    n = input('Input exponent n: ');

    % Error check
    while ~isnumeric(a)
        fprintf('Invalid coefficient.\n');
        a = input('Input coefficient a: ');
    end
    while ~isnumeric(n)
        fprintf('Invalid exponent.\n');
        n = input('Input exponent n: ');
    end

    % Compute derivative of this term
    if n == 0
        term = sym(0);
    elseif n < 0
        term = a * n / x^(abs(n)+1);
        fprintf("Now displaying the steps for term %d\n", count);
        fprintf("Step 1: Mulitply a * n: %d\n", a*n);
        fprintf("Step 2: Subtract 1 from n: %d\n", n-1);
        fprintf("Step 3: Move the variable x^(-n) under the coeffcient and do the absolute value of n: %d\n", abs(n));
        fprintf('Term %d derivative: %s\n', count, char(term));
    else
        term = a * n * x^(n-1);
        fprintf("Now displaying the steps for term %d\n", count);
        fprintf("Step 1: Mulitply a * n: %d\n", a*n);
        fprintf("Step 2: Subtract 1 from n: %d\n", n-1);
        fprintf('Term %d derivative: %s\n', count, char(term));
    end

    % Recursive call for remaining terms
    productEquation = term + Powerrule(t-1, count+1);
end

%% -----------------------------
% L'Hopital's Rule (kept)
%% -----------------------------
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

%% -----------------------------
% Placeholders for other rules
%% -----------------------------
function productRule()
    clc; syms x
    disp('You selected Product Rule.');
    disp('------------------------------------------');
    disp('Formula: d(uv)/dx = u''*v + u*v''');
    disp('You will enter the two functions u(x) and v(x).');
    disp('Supported types: polynomial, log, trig, exp');
    disp('------------------------------------------');

    % Input functions
    u_str = input('Enter the first function u(x): ','s');
    v_str = input('Enter the second function v(x): ','s');

    % Types
    u_type = input('Type of u(x) (poly/log/trig/exp): ','s');
    v_type = input('Type of v(x) (poly/log/trig/exp): ','s');

    disp(' ');
    disp('--- Step-by-Step Derivative ---');

    % Derivative of u
    fprintf('Step 1: Compute derivative of u(x) = %s\n', u_str);
    du = derivativeStepProduct(u_str, u_type);
    fprintf('Step 2: u''(x) = %s\n', char(du));

    % Derivative of v
    fprintf('Step 3: Compute derivative of v(x) = %s\n', v_str);
    dv = derivativeStepProduct(v_str, v_type);
    fprintf('Step 4: v''(x) = %s\n', char(dv));

    % Apply product rule
    fprintf('Step 5: Apply Product Rule: u''*v + u*v''\n');
    y_prime = simplify(du*str2sym(v_str) + str2sym(u_str)*dv);
    fprintf('Step 6: Final derivative: dy/dx = %s\n', char(y_prime));
    disp('------------------------------------------');

    % ASK NEXT MOVE
    again = input('Would you like to solve another product rule problem? y/n: ','s');
    if lower(again) == 'y'
        productRule();
        return;
    end

    back = input('Would you like to go back to the main menu? y/n: ','s');
    if lower(back) == 'y'
        MathMachine();  % your main machine function
    else
        disp('Exiting Math Machine...');
    end
end


function df = derivativeStepProduct(f_str, f_type)
    syms x

    % Fix log → ln automatically
    f_str = regexprep(f_str, 'log', 'ln');

    switch lower(f_type)

        case 'poly'
            f_sym = str2sym(f_str);
            df = diff(f_sym, x);

        case 'log'
            df = LogRule(str2sym(f_str));

        case 'trig'
            f_sym = str2sym(f_str);
            df = diff(f_sym, x);

        case 'exp'
            f_sym = str2sym(f_str);
            df = diff(f_sym, x);

        otherwise
            error('Unsupported function type. Use poly, log, trig, or exp.');
    end
end


function quotientRule()
    clc; syms x
    disp('You selected Quotient Rule.');
    disp('------------------------------------------');
    disp('Formula: d(u/v)/dx = (u''*v - u*v'') / v^2');
    disp('You will enter the numerator u(x) and denominator v(x).');
    disp('Supported types: polynomial, log, trig, exp, powers, sums, products');
    disp('------------------------------------------');

    while true
        % Input functions
        u_str = input('Enter the numerator u(x) (ex. x^2): ','s');
        v_str = input('Enter the denominator v(x) (ex. exp(x^2)): ','s');

        % Convert to symbolic expressions
        try
            u = str2sym(u_str);
            v = str2sym(v_str);
        catch
            disp('ERROR: Invalid symbolic expression. Please check your input.');
            continue; % ask again
        end

        disp(' ');
        disp('--- Step-by-Step Derivative ---');

        % Derivative of numerator
        fprintf('Step 1: Compute derivative of numerator u(x) = %s\n', u_str);
        [du, steps_u] = computeChainDerivativeSteps(u, x);
        fprintf('Step 2: u''(x) = %s\n', char(du));

        % Derivative of denominator
        fprintf('Step 3: Compute derivative of denominator v(x) = %s\n', v_str);
        [dv, steps_v] = computeChainDerivativeSteps(v, x);
        fprintf('Step 4: v''(x) = %s\n', char(dv));

        % Apply quotient rule
        fprintf('Step 5: Apply Quotient Rule: (u''*v - u*v'') / v^2\n');
        y_prime = simplify((du*v - u*dv)/v^2);
        fprintf('Step 6: Final derivative: dy/dx = %s\n', char(y_prime));
        disp('------------------------------------------');

        % Optional: display detailed steps
        disp('--- Detailed Steps for u(x) ---');
        for k = 1:length(steps_u)
            fprintf('%s\n', steps_u{k});
        end
        disp('--- Detailed Steps for v(x) ---');
        for k = 1:length(steps_v)
            fprintf('%s\n', steps_v{k});
        end

        % Ask if user wants to solve another Quotient Rule problem
        solveAnother = '';
        while ~ismember(solveAnother, {'y','n'})
            solveAnother = lower(strtrim(input('Would you like to solve another Quotient Rule problem? (y/n): ','s')));
        end
        if strcmp(solveAnother, 'y')
            clc;
            continue; % restart loop
        end

        % Ask if user wants to return to main menu
        returnMenu = '';
        while ~ismember(returnMenu, {'y','n'})
            returnMenu = lower(strtrim(input('Would you like to return to the main menu? (y/n): ','s')));
        end
        if strcmp(returnMenu, 'y')
            clc;
            return; % back to main menu
        else
            disp('Exiting The Math Machine. Goodbye!');
            error('ExitMathMachine'); % terminate program
        end
    end
end



function implicitDiff()
    clc;
    syms x y
    disp('You selected Implicit Differentiation.');
    disp('------------------------------------------');
    disp('Enter your equation F(x, y) = 0');
    disp('Example: x^2 + y^2 - 1 = 0');
    disp('------------------------------------------');

    while true
        eqn_input = input('Enter your equation F(x, y) = 0: ', 's');

        try
            F = str2sym(eqn_input);  % convert string to symbolic
        catch
            disp('ERROR: Invalid equation. Please enter a valid symbolic equation.');
            continue;
        end

        fprintf('\nOriginal equation: F(x, y) = %s = 0\n\n', char(F));

        try
            % Compute dF/dx
            disp('Step 1: Compute dF/dx (derivative with respect to x)');
            [dFdx, dx_steps] = termWiseDerivativeClean(F, x);
            for k=1:length(dx_steps)
                disp(dx_steps{k});
            end
            fprintf('dF/dx = %s\n\n', char(dFdx));

            % Compute dF/dy
            disp('Step 2: Compute dF/dy (derivative with respect to y)');
            [dFdy, dy_steps] = termWiseDerivativeClean(F, y);
            for k=1:length(dy_steps)
                disp(dy_steps{k});
            end
            fprintf('dF/dy = %s\n\n', char(dFdy));

            % Solve for dy/dx
            disp('Step 3: Solve for dy/dx = - (dF/dx) / (dF/dy)');
            dydx = -dFdx/dFdy;
            fprintf('Implicit derivative dy/dx = %s\n', char(simplify(dydx)));
            disp('------------------------------------------');

        catch ME
            disp('ERROR while computing implicit derivative:');
            disp(ME.message);
        end

        % Ask user if they want to try another
        retry = '';
        while ~ismember(retry, {'y','n'})
            retry = lower(strtrim(input('Would you like to try another implicit problem? (y/n): ', 's')));
            if ~ismember(retry, {'y','n'})
                disp('Invalid input. Please enter y or n.');
            end
        end
        if strcmp(retry, 'n')
            clc;
            return;
        else
            clc;
            disp('You selected Implicit Differentiation.');
            disp('------------------------------------------');
            disp('Enter your equation F(x, y) = 0');
            disp('Example: x^2 + y^2 - 1 = 0');
            disp('------------------------------------------');
        end
    end
end

%% -----------------------------
%% Term-wise derivative with clean steps
%% -----------------------------
function [deriv, steps] = termWiseDerivativeClean(expr, var)
    syms x y
    steps = {};

    if isempty(find(symvar(expr) == var, 1))
        deriv = sym(0);
        return;
    end

    if isempty(symvar(expr))
        deriv = sym(0);
        return;
    end

    if isequal(expr, var)
        deriv = sym(1);
        steps{end+1} = sprintf('Derivative of %s with respect to %s = 1', char(expr), char(var));
        return;
    end

    ch = children(expr);
    op = char(feval(symengine,'op', expr, 0));

    % Standardize operator names
    if strcmp(op,'plus')
        op = '_plus';
    elseif strcmp(op,'times') || strcmp(op,'_mult')
        op = '_times';
    elseif strcmp(op,'power')
        op = '_power';
    elseif strcmp(op,'minus')
        op = '_minus';
    end

    % Sum
    if strcmp(op,'_plus')
        deriv = 0;
        for k=1:length(ch)
            [dterm, s] = termWiseDerivativeClean(ch{k}, var);
            deriv = deriv + dterm;
            steps = [steps, s]; %#ok<AGROW>
        end
        return;
    end

    % Difference
    if strcmp(op,'_minus')
        if length(ch)==1
            [dterm, s] = termWiseDerivativeClean(ch{1}, var);
            deriv = -dterm;
            steps = [steps, strcat('Derivative of -', char(ch{1}), ' with respect to ', char(var), ' = -', char(dterm))];
        else
            [d1, s1] = termWiseDerivativeClean(ch{1}, var);
            [d2, s2] = termWiseDerivativeClean(ch{2}, var);
            deriv = d1 - d2;
            steps = [steps, s1, s2];
        end
        return;
    end

    % Product
    if strcmp(op,'_times')
        if length(ch)==2
            [u_der, s1] = termWiseDerivativeClean(ch{1}, var);
            [v_der, s2] = termWiseDerivativeClean(ch{2}, var);
            deriv = u_der*ch{2} + ch{1}*v_der;
            steps = [steps, s1, s2];
        else
            deriv = 0;
            for i=1:length(ch)
                term = 1;
                for j=1:length(ch)
                    if i==j
                        [dterm, s] = termWiseDerivativeClean(ch{j}, var);
                        term = term*dterm;
                        steps = [steps, s];
                    else
                        term = term*ch{j};
                    end
                end
                deriv = deriv + term;
            end
        end
        return;
    end

    % Power
    if strcmp(op,'_power')
        base = ch{1};
        exponent = ch{2};
        if isempty(symvar(exponent))
            [base_der, s] = termWiseDerivativeClean(base,var);
            deriv = exponent*base^(exponent-1)*base_der;
            steps{end+1} = sprintf('Derivative of %s with respect to %s = %s', char(expr), char(var), char(deriv));
        else
            deriv = expr*(termWiseDerivativeClean(exponent,var)*log(base) + (exponent/base)*termWiseDerivativeClean(base,var));
            steps{end+1} = sprintf('Derivative of %s with respect to %s (general power rule)', char(expr), char(var));
        end
        return;
    end

    % Trig / exp / log
    if strcmp(op,'sin')
        u = ch{1};
        [u_der, s] = termWiseDerivativeClean(u,var);
        deriv = cos(u)*u_der;
        steps = [steps, s, {sprintf('Derivative of sin(%s) with respect to %s = cos(%s)*d(%s)/d%s', char(u), char(var), char(u), char(u), char(var))}];
        return;
    end
    if strcmp(op,'cos')
        u = ch{1};
        [u_der, s] = termWiseDerivativeClean(u,var);
        deriv = -sin(u)*u_der;
        steps = [steps, s, {sprintf('Derivative of cos(%s) with respect to %s = -sin(%s)*d(%s)/d%s', char(u), char(var), char(u), char(u), char(var))}];
        return;
    end
    if strcmp(op,'tan')
        u = ch{1};
        [u_der, s] = termWiseDerivativeClean(u,var);
        deriv = sec(u)^2*u_der;
        steps = [steps, s, {sprintf('Derivative of tan(%s) with respect to %s = sec(%s)^2*d(%s)/d%s', char(u), char(var), char(u), char(u), char(var))}];
        return;
    end
    if strcmp(op,'exp')
        u = ch{1};
        [u_der, s] = termWiseDerivativeClean(u,var);
        deriv = exp(u)*u_der;
        steps = [steps, s, {sprintf('Derivative of exp(%s) with respect to %s = exp(%s)*d(%s)/d%s', char(u), char(var), char(u), char(u), char(var))}];
        return;
    end
    if strcmp(op,'log')
        u = ch{1};
        [u_der, s] = termWiseDerivativeClean(u,var);
        deriv = (1/u)*u_der;
        steps = [steps, s, {sprintf('Derivative of log(%s) with respect to %s = (1/%s)*d(%s)/d%s', char(u), char(var), char(u), char(u), char(var))}];
        return;
    end

    % Unknown
    error('Unsupported operation: %s', op);
end
