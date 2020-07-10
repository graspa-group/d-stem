%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% D-STEM - Distributed Space Time Expecation Maximization              %
%%%                                                                      %
%%% Author: Francesco Finazzi                                            %
%%% E-mail: francesco.finazzi@unibg.it                                   %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Information and                    %
%%%              Production Engineering                                  %
%%% Author website: http://www.unibg.it/pers/?francesco.finazzi          %
%%%                                                                      %
%%% Author: Yaqiong Wang                                                 %
%%% E-mail: yaqiongwang@pku.edu.cn                                       %
%%% Affiliation: Peking University,                                      %
%%%              Guanghua school of management,                          %
%%%              Business Statistics and Econometrics                    %
%%%                                                                      %
%%% Author: Alessandro Fassò                                             %
%%% E-mail: alessandro.fasso@unibg.it                                    %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Information and                    %
%%%              Production Engineering                                  %
%%% Author website: http://www.unibg.it/pers/?alessandro.fasso           %
%%%                                                                      %
%%% Code website: https://github.com/graspa-group/d-stem                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

In order to run a demo:

1) Open MATLAB;
2) Select a demo folder between DCM, HDGM or f-HDGM as Current Folder;
3) Open and run the demo_*.m script;
4) Choose an option from the menu.

In order to estimate a model on your own data set:

1) Add the D-STEM Src folder and sub-folders to the Matlab seach path list
2) Adapt any script under DCM/Scripts, HDGM/Scripts or f-HDGM/Scripts to
   your data set

Please note that:

a) The following MATLAB toolboxes are required:
   - Statistics and Machine Learning Toolbox
   - Optimization Toolbox
   - Mapping Toolbox
   - Parallel Computing Toolbox 

b) The Mapping Toolbox is only required for plotting maps while it is not requires for model estimation

c) The Parallel Computing Toolbox is only needed to exploit the parallel capabilities of D-STEM

d) The code has been developed and tested on MATLAB R2019a