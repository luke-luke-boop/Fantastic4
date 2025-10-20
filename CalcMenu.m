disp('you have selected the calc menu');
disp('here are your options: 1. Limits and Continuity\n 2. Derivatives\n 3. Application of Derivatives\n 4. Integration');


userChoice = ('please choose a number between 1 and 4: \n');
choice = input(userChoice);

choice = round(choice);

if choice > 4 || choice < 1
    fprintf('Invalid selection. Please choose a number between 1 and 4. \n');
    choice = input(userChoice);
end


switch choice
    case 1
        disp('You selected Limits and Continuity. Please enter the function:');
        disp(' ');

    case 2
        disp('You selected Derivatives. Please enter the function:');
        disp(' ');

    case 3
        disp("You selected application of Derivatives. Please enter the function");
        disp(' ');
       
    case 4
        disp('You selected Integration. Please enter the function:');
        disp(' ');
end