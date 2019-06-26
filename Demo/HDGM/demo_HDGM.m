clc
clearvars
close all

while (1)
    disp(' ');
    disp('****************************************************************');
    pause(0.05);
    disp('*                          HDGM demo                           *');
    pause(0.05);
    disp('****************************************************************');
    pause(0.05);
    disp(' ');
    disp('Choose which script to run:');
    pause(0.05);
    disp(' ');
    disp('[1] - Univariate model (< 3 minutes)');
    pause(0.05);
    disp('[2] - Bivariate  model (< 3 minutes)');
    pause(0.05);
    disp(' ');
    disp('[3] - Exit');
    disp(' ');
    choice  =  input('Insert a number and press Enter: ','s');
    switch choice
        case '1'
            run('Scripts/demo_univariate');
        case '2'
            run('Scripts/demo_bivariate');
        case '3'
            break;
        otherwise
            disp(' ');
            disp('Please insert a number between 1 and 3');
            pause(2.5);
    end
end
