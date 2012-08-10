classdef stem_residual < stem_data
    properties
        kriging_Var_W_bar_hat=[];
        %leave one out properties
        rmse=[];
        logL=[];
    end
    
    methods
        function obj = stem_residual(stem_data,stem_varset)
            %residual must be a cell array with the same structure of the orignal Y
            obj=obj@stem_data(stem_varset,stem_data.stem_gridlist,stem_data.stem_datestamp,[],stem_data.shape,0); %call the superclass constructor
        end
    end
end