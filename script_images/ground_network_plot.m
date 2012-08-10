clc
clear all
s=shaperead('../Maps/CNTR_BN_03M_2010');
load ../Data/st_model_20120603_135334.mat


mapshow(s,'Color','black');
hold on
