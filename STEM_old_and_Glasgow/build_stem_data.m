function [stem_data_record, colocated_data] = build_stem_data(colocated_data,data_settings)
% X è il temporaneo di colocated_data
% Y è il temporaneo di stem_data

disp('Building stem data...')

n_var=length(colocated_data.Y);
n_cov=size(colocated_data.X,2);
T=size(colocated_data.Y{1,1}.data,2);

for i=1:n_var
    s=size(colocated_data.Y{i}.grid);
    stem_data_record.var_dims(i)=s(1);
end

N=sum(stem_data_record.var_dims);

if sum(strcmp('latitude',data_settings.covariates))
    n_cov=n_cov+1;
    for i=1:n_var
        colocated_data.X{i,n_cov}.var_name=colocated_data.Y_name{i};
        colocated_data.X{i,n_cov}.cov_name='latitude';
        colocated_data.X{i,n_cov}.data=repmat(colocated_data.Grid{1,i}(:,1),1,T);
    end
end

if sum(strcmp('longitude',data_settings.covariates))
    n_cov=n_cov+1;
    for i=1:n_var
        colocated_data.X{i,n_cov}.var_name=colocated_data.Y_name{i};
        colocated_data.X{i,n_cov}.cov_name='longitude';
        colocated_data.X{i,n_cov}.data=repmat(colocated_data.Grid{1,i}(:,2),1,T);
    end
end

if data_settings.add_constant
    n_cov=n_cov+1;
    for i=1:n_var
        colocated_data.X{i,n_cov}.var_name=colocated_data.Y_name{i};
        colocated_data.X{i,n_cov}.cov_name='constant';
        colocated_data.X{i,n_cov}.data=ones(stem_data_record.var_dims(i),T);
    end
end

stem_data_record.Grid=colocated_data.Grid;
stem_data_record.Grid_size=colocated_data.Grid_size;

variables_dim=sum(stem_data_record.var_dims);
beta_dim=n_cov*n_var;

%inizializzazione matrice rettangolare Xbig con valori covariate per ogni t
stem_data_record.X=zeros(variables_dim,beta_dim,T);

%trovare eventualmente modo più furbo per creare la matrice senza dover
%utilizzare n_col e n_row
n_col=1;
n_row=1;
for j=1:length(colocated_data.Y)
    for i=1:length(colocated_data.X)
        for t=1:T
            stem_data_record.X(n_row:n_row+stem_data_record.var_dims(j)-1,n_col,t)=colocated_data.X{j,i}.data(:,t);
        end
        n_col=n_col+1;
    end
    n_row=n_row+stem_data_record.var_dims(j);
end

%costruzione matrice Y di dimenzione NxT
stem_data_record.Y=[];
for i=1:length(colocated_data.Y)
    stem_data_record.Y=[stem_data_record.Y; colocated_data.Y{i}.data];
end

%creazione della matrice delle distanze e salvataggio in Y
if 1
    stem_data_record.DistMat=stem_data.get_distance_matrix(stem_data_record.Grid);
else
    elevation_data=build_elevation_data();
    stem_data_record.DistMat=stem_data.get_orographic_distance_matrix(stem_data_record.Grid,elevation_data);
end
stem_data_record.variables_name=data_settings.variables;
stem_data_record.covariates_name=data_settings.covariates;
stem_data_record.data_settings=data_settings;
end

