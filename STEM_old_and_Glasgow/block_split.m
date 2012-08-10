%input: full X covariates matrix and mask matrix based on elevation data (water=missing)
%output: a set of .mat files with block_size rows each with only the rows
%of X corresponding to land pixels

clc
load ../Data/krig_X_scotland2009
load ../Data/krig_mask_scotland2009
mask=mask(:);
L=isnotnan(mask);
slv=slv(L,:,:);
d=size(slv,1);
block_size=500;
blocks=0:block_size:d;
if not(blocks(end)==d)
    blocks=[blocks d];
end
for i=1:length(blocks)-1
    X_krig_block=slv(blocks(i)+1:blocks(i+1),:,:);
    code=num2str(i);
    while not(length(code)==5)
        code=['0',code];
    end
    save(['../Data/X_krig_blocks/block',code],'X_krig_block','-v7.3');
    i
end