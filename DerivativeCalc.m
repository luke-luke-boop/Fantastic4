clear
clc
%this is a trial for the commit

%Initial menu display
disp('Math Machine: ')
disp('')
mathex = ('Chose which math course your problem relates to: \n 1. Algebra II \n 2. Pre-Calculus \n 3. Calculus I \n 4. Calculus II \n');
userexample = input(mathex);

userexample = round(userexample);

%Error check
if userexample > 4 || userexample < 1
    fprintf('Invalid selection. Please choose a number between 1 and 4. \n');
    userexample = input(mathex);
end

%Switch case to determine which scripts to call
switch userexample
    case 1
        disp('You selected Algebra II.');
        disp('')
        mathprob = ('Chose your type of problem: \n');
        userprob = input(mathprob);

    case 2
        disp('You selected Pre-Calculus.');
        disp('')
        mathprob = ('Chose your type of problem: \n');
        userprob = input(mathprob);

    case 3
        disp('You selected Calculus I.');
        disp('')
        mathprob = ('Chose your type of problem: \n');
        userprob = input(mathprob);

    case 4
        disp('You selected Calculus II.');
        disp('')
        mathprob = ('Chose your type of problem: \n 1. Derivation \n 2. Intergration \n 3. Trig-sub \n 4. Series \n');
        userprob = input(mathprob);
end
