function [colocated_data] = build_colocate_data(raw_data)
% X è il temporaneo di raw_data
% Y è il temporaneo di colocated_data

disp('Build colocated data')
X=raw_data;

Y.Y=X.Y;

for j=1:length(X.Y)
    grid_Y=X.Y{j}.grid;
    Y.Grid{j}=grid_Y;
    Y.Grid_size{j}=X.Y{j}.grid_size;
    s=size(X.Y{j}.data);
    for i=1:length(X.X)
        if 1 %eseguo colocazione su tutte le combinazioni 
            cdisp('Colocating (Y,X) = ', [j,i])
            grid_X=X.X{i}.grid;
            Y.X{j,i}.var_name=X.Y{j}.name;
            Y.X{j,i}.cov_name=X.X{i}.name;
            lat_X=reshape(grid_X(:,1),X.X{i}.grid_size(1,1),X.X{i}.grid_size(1,2));
            long_X=reshape(grid_X(:,2),X.X{i}.grid_size(1,1),X.X{i}.grid_size(1,2));
            for t=1:s(2)
                data_X=reshape(X.X{i}.data(:,t),X.X{i}.grid_size(1,1),X.X{i}.grid_size(1,2));
                if strcmp(X.Y{j}.grid_type,'regular')
                    lat_Y=reshape(grid_Y(:,1),X.Y{j}.grid_size(1,1),X.Y{j}.grid_size(1,2));
                    long_Y=reshape(grid_Y(:,2),X.Y{j}.grid_size(1,1),X.Y{j}.grid_size(1,2));
                    n_Y=length(X.Y{j}.grid);
                    Y.X{j,i}.data(:,t)=reshape(interp2(long_X',lat_X',data_X',long_Y',lat_Y','linear')',n_Y,1);
                else
                    for k=1:s(1)
                        Y.X{j,i}.data(k,t)=interp2(long_X',lat_X',data_X',X.Y{j}.grid(k,2),X.Y{j}.grid(k,1),'linear');
                    end
                end
            end
        end
    end
end

for j=1:length(X.Y)
    Y.Y_name{j}=X.Y{j}.name;
end
for j=1:length(X.X)
    Y.X_name{j}=X.X{j}.name;
end

colocated_data=Y;
