clc
clearvars
close all

while(1)
    disp(' ');
    disp('****************************************************************');
    pause(0.05);
    disp('*                          f-HDGM Demo                         *');
    pause(0.05)
    disp('****************************************************************');
    pause(0.05);
    disp(' ');
    disp('Choose which script to run:');
    pause(0.05);
    disp(' ');
    disp('[1] - Ozone data: model estimation (< 10 minutes)');
    pause(0.05);
    disp('[2] - Ozone data: cross-validation (< 3 minutes)');
    pause(0.05);
    disp('[3] - Ozone data: dynamic kriging (> 10 minutes)');
    pause(0.05)
    disp('[4] - ROAB data: model estimation and kriging (< 5 minutes)');
    pause(0.05)
    disp(' ');
    disp('[5] - Exit');
    disp(' ');
    choice  =  input('Insert a number and press Enter: ','s');
    switch choice
        case '1'
            run('O3/Scripts/demo_section4_model_estimate');
        case '2'
            run('O3/Scripts/demo_section4_crossvalidation');
        case '3'
            run('O3/Scripts/demo_section4_kriging');
        case '4'
            run('Radiosonde/demo_section5')
        case '5'
            break;
        otherwise
            disp(' ');
            disp('Please insert a number between 1 and 5');
            pause(2.5);
    end
end    
