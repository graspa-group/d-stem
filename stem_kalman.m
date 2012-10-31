%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef stem_kalman < handle
    
    properties
        stem_model=[];
    end
    
    methods
        function obj = stem_kalman(stem_model)
            if strcmp(class(stem_model),'stem_model')
                obj.stem_model=stem_model;
            else
                error('The input argument must be of class stem_model');
            end
        end
        
        function [st_kalmanfilter_result,sigma_eps,sigma_W_r,sigma_W_g,sigma_Z,aj_rg,aj_g,M,sigma_geo] = filter(obj,compute_logL,time_steps,pathparallel)
            if nargin<2
                compute_logL=0;
            end
            if nargin<3
                pathparallel=[];
                time_steps=[];
            end
            if nargin==3
                error('The pathparallel input argument must be provided');
            end
            disp('    Kalman filter started...');
            ct1=clock;
            
            z0=zeros(obj.stem_model.stem_par.p,1);
            P0=eye(obj.stem_model.stem_par.p);
            time_diagonal=obj.stem_model.stem_par.time_diagonal;
            
            data=obj.stem_model.stem_data;
            par=obj.stem_model.stem_par;            
            
            [sigma_eps,sigma_W_r,sigma_W_g,sigma_geo,sigma_Z,aj_rg,aj_g,M] = obj.stem_model.get_sigma();
            
            tapering=obj.stem_model.tapering;            
            if isempty(pathparallel)
                [zk_f,zk_u,Pk_f,Pk_u,J,logL] = stem_kalman.Kfilter(data.Y,data.X_rg,data.X_beta,data.X_time,data.X_g,par.beta,par.G,par.sigma_eta,sigma_W_r,sigma_W_g,sigma_eps,sigma_geo,aj_rg,aj_g,M,z0,P0,time_diagonal,tapering,compute_logL);
            else
                [zk_f,zk_u,Pk_f,Pk_u,J,logL] = stem_kalman.Kfilter_parallel(data.Y,data.X_rg,data.X_beta,data.X_time,data.X_g,par.beta,par.G,par.sigma_eta,sigma_W_r,sigma_W_g,sigma_eps,sigma_geo,aj_rg,aj_g,M,z0,P0,time_diagonal,time_steps,pathparallel,tapering,compute_logL);
            end
            st_kalmanfilter_result = stem_kalmanfilter_result(zk_f,zk_u,Pk_f,Pk_u,J,logL);
            
            ct2=clock;
            disp(['    Kalman filter ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
        end
        
        function [st_kalmansmoother_result,sigma_eps,sigma_W_r,sigma_W_g,sigma_Z,aj_rg,aj_g,M,sigma_geo] = smoother(obj,compute_logL,time_steps,pathparallel)
            if nargin<2
              compute_logL=0;
            end
            if nargin<3
                pathparallel=[];
                time_steps=[];
            end
            if nargin==3
                error('The pathparallel input argument must be provided');
            end
            disp('    Kalman smoother started...');
            ct1=clock;
            z0=zeros(obj.stem_model.stem_par.p,1);
            P0=eye(obj.stem_model.stem_par.p);
            time_diagonal=obj.stem_model.stem_par.time_diagonal;
            
            data=obj.stem_model.stem_data;
            par=obj.stem_model.stem_par;
            
            [sigma_eps,sigma_W_r,sigma_W_g,sigma_geo,sigma_Z,aj_rg,aj_g,M] = obj.stem_model.get_sigma();
            
            tapering=obj.stem_model.tapering;
            [zk_s,Pk_s,PPk_s,logL] = obj.Ksmoother(data.Y,data.X_rg,data.X_beta,data.X_time,...
                data.X_g,par.beta,par.G,par.sigma_eta,sigma_W_r,...
                sigma_W_g,sigma_eps,sigma_geo,aj_rg,aj_g,M,z0,P0,...
                time_diagonal,time_steps,pathparallel,tapering,compute_logL);
            st_kalmansmoother_result = stem_kalmansmoother_result(zk_s,Pk_s,PPk_s,logL);
            ct2=clock;
            disp(['    Kalman smoother ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
        end
    end
    
    methods (Static)
        
        function [zk_f,zk_u,Pk_f,Pk_u,J,logL] = Kfilter(Y,X_rg,X_beta,X_time,X_g,beta,G,sigma_eta,sigma_W_r,sigma_W_g,sigma_eps,sigma_geo,aj_rg,aj_g,M,z0,P0,time_diagonal,tapering,compute_logL)
            if nargin<13
                error('You have to provide all the input arguments');
            end
                        
            if size(X_beta,2)~=length(beta)
                error('X_beta and beta are not conformable');
            end
           
            if size(G,1)~=size(G,2)
                error('G must be square');
            end
            
            if size(sigma_eta,1)~=size(sigma_eta,2)
                error('sigma_eta must be square');
            end      
            
            if size(G,1)~=size(sigma_eta,1)
                error('G and sigma_eta must have the same dimensions');
            end
            
            if size(Y,1)~=size(sigma_eps,1)
                error('The dimensions of sigma_eps must be equal to the number of rows of Y');
            end
            
            if size(z0,1)~=size(G,1)
                error('The length of z0 must be equal to the dimensions of G');
            end
            
            if size(P0,1)~=size(P0,2)
                error('P0 must be square');
            end
            
            if size(P0,1)~=size(G,1)
                error('The dimensions of P0 must be equal to the dimensions of G and sigma_eps');
            end
            
            if isempty(sigma_geo)
                compute_sigma_geo=1;
            else
                compute_sigma_geo=0;
            end
           
            p=size(G,1);
            N=size(Y,1);
            T=size(Y,2);
            zk_f=zeros(p,T+1);
            zk_u=zeros(p,T+1);
            Pk_f=zeros(p,p,T+1);
            Pk_u=zeros(p,p,T+1);
            J=zeros(p,size(Y,1),T+1);
            innovation=zeros(size(Y,1),T);
            
            zk_u(:,1)=z0;
            Pk_u(:,:,1)=P0;
            
            logL=0;
            for t=2:T+1
                if size(X_time,3)==1
                    tK=2;
                else
                    tK=t; %time variant
                end
                
                if compute_sigma_geo
                    if not(isempty(X_rg))
                        sigma_geo=zeros(N);
                        if size(X_rg,3)>1
                            sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_r,M,'b'),X_rg(:,1,t-1),'b'),aj_rg,'b');
                        else
                            sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_r,M,'b'),X_rg(:,1,1),'b'),aj_rg,'b');
                        end
                    end

                    if not(isempty(X_g))
                        if isempty(X_rg)
                            if tapering
                                sigma_geo=spalloc(size(sigma_W_g{1},1),size(sigma_W_g{1},1),nnz(sigma_W_g{1}));
                            else
                                sigma_geo=zeros(N);
                            end
                        end                        
                        for k=1:size(X_g,4)
                            if size(X_g,3)>1
                               sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_g{k},X_g(:,1,t-1,k),'b'),aj_g(:,k),'b');
                            else
                               sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_g{k},X_g(:,1,1,k),'b'),aj_g(:,k),'b');
                            end
                        end
                    end
                    if isempty(X_g)&&isempty(X_rg)
                        sigma_geo=sigma_eps;
                    else
                        sigma_geo=sigma_geo+sigma_eps;
                    end
                end
                
                if size(X_beta,3)==1
                    tX=2;
                else
                    tX=t; %time variant
                end
                Lt=not(isnan(Y(:,t-1))); %note the t-1
                temp=sigma_geo(Lt,Lt);
                if tapering
                    if not(stem_misc.isdiagonal(temp))
                        r = symamd(temp);
                        c=chol(temp(r,r));
                        temp2=speye(sum(Lt));
                        temp3=full(stem_misc.chol_solve(c,temp2(r,:)));
                        sigma_geo_inv=zeros(size(temp3));
                        sigma_geo_inv(r,:)=temp3;
                        clear temp2
                        clear temp3
                    else
                        d=1./diag(temp);
                        sigma_geo_inv=sparse(1:length(d),1:length(d),d);
                    end
                else
                    if not(stem_misc.isdiagonal(temp))
                        c=chol(sigma_geo(Lt,Lt));
                        sigma_geo_inv=stem_misc.chol_solve(c,eye(sum(Lt)));
                    else
                        sigma_geo_inv=diag(1./diag(temp));
                    end
                end
                
                if not(time_diagonal)
                    %filter
                    zk_f(:,t)=G*zk_u(:,t-1); %(6.19) Stoffer
                    Pk_f(:,:,t)=G*Pk_u(:,:,t-1)*G'+sigma_eta; %(6.20) Stoffer
                    
                    %update
                    %original formula
                    %J(i,Lt,t)=Pk_f(i,i,t)*X_time(Lt,i,tK-1)'/(X_time(Lt,i,tK-1)*Pk_f(i,i,t)*X_time(Lt,i,tK-1)'+sigma_geo(Lt,Lt)); %(6.23) Stoffer
                    %Sherman-Morrison-Woodbury formula: (B*P*B+D)^-1=D^-1-D^-1*B(P^-1+B*D^-1*B)^-1*B*D^-1
                    %J(i,Lt,t)=Pk_f(i,i,t)*X_time(Lt,i,tK-1)'*(sigma_geo_inv-sigma_geo_inv*X_time(Lt,i,tK-1)/(1/Pk_f(i,i,t)+X_time(Lt,i,tK-1)'*sigma_geo_inv*X_time(Lt,i,tK-1))*(X_time(Lt,i,tK-1)'*sigma_geo_inv));
                    
                    temp=X_time(Lt,:,tK-1)'*sigma_geo_inv; %note that temp can be computed in a distributed way before the KF is started
                    temp2=temp*X_time(Lt,:,tK-1);
                    if compute_logL
                        sigma_t_inv=sigma_geo_inv-(temp'/((Pk_f(:,:,t)\eye(size(temp2)))+temp2))*temp;
                    end
                    
                    temp3=sparse(Pk_f(:,:,t)*X_time(Lt,:,tK-1)');
                    J(:,Lt,t)=Pk_f(:,:,t)*temp-temp3*(temp'/(Pk_f(:,:,t)\eye(size(temp2))+temp2)*temp);
                         
                    if not(isempty(X_beta))
                        innovation(Lt,t-1)=Y(Lt,t-1)-X_beta(Lt,:,tX-1)*beta-X_time(Lt,:,tK-1)*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                    else
                        innovation(Lt,t-1)=Y(Lt,t-1)-X_time(Lt,:,tK-1)*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                    end
                    
                    zk_u(:,t)=zk_f(:,t)+J(:,Lt,t)*innovation(Lt,t-1); 
                    Pk_u(:,:,t)=(eye(p)-J(:,Lt,t)*X_time(Lt,:,tK-1))*Pk_f(:,:,t);  %(6.22) Stoffer
                else
                    %filter
                    zk_f(:,t)=diag(G).*zk_u(:,t-1); %(6.19) Stoffer
                    Pk_f(:,:,t)=diag(diag(G).^2.*diag(Pk_u(:,:,t-1))+diag(sigma_eta)); %(6.20) Stoffer
                    
                    %update

                    %original formula
                    %J(i,Lt,t)=Pk_f(i,i,t)*X_time(Lt,i,tK-1)'/(X_time(Lt,i,tK-1)*Pk_f(i,i,t)*X_time(Lt,i,tK-1)'+sigma_geo(Lt,Lt)); %(6.23) Stoffer
                    %Sherman-Morrison-Woodbury formula: (B*P*B+D)^-1=D^-1-D^-1*B(P^-1+B*D^-1*B)^-1*B*D^-1
                    %J(i,Lt,t)=Pk_f(i,i,t)*X_time(Lt,i,tK-1)'*(sigma_geo_inv-sigma_geo_inv*X_time(Lt,i,tK-1)/(1/Pk_f(i,i,t)+X_time(Lt,i,tK-1)'*sigma_geo_inv*X_time(Lt,i,tK-1))*(X_time(Lt,i,tK-1)'*sigma_geo_inv));
                    
                    temp=X_time(Lt,:,tK-1)'*sigma_geo_inv; %note that temp can be computed in a distributed way before the KF is started
                    temp2=temp*X_time(Lt,:,tK-1);
                    P=diag(1./diag(Pk_f(:,:,t)));
                    if compute_logL
                        sigma_t_inv=sigma_geo_inv-(temp'/(P+temp2))*temp;
                    end
                    temp3=Pk_f(:,:,t)*X_time(Lt,:,tK-1)';
                    J(:,Lt,t)=Pk_f(:,:,t)*temp-temp3*(temp'/(P+temp2)*temp);
                    
                    if not(isempty(X_beta))
                        innovation(Lt,t-1)=Y(Lt,t-1)-X_beta(Lt,:,tX-1)*beta-X_time(Lt,:,tK-1)*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                    else
                        innovation(Lt,t-1)=Y(Lt,t-1)-X_time(Lt,:,tK-1)*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                    end
                    
                    zk_u(:,t)=zk_f(:,t)+J(:,Lt,t)*innovation(Lt,t-1);
                    Pk_u(:,:,t)=diag(diag((eye(p)-J(:,Lt,t)*X_time(Lt,:,tK-1))).*diag(Pk_f(:,:,t))); %(6.22) Stoffer
                end
                if compute_logL
                    r = symamd(sigma_t_inv);
                    c=chol(sigma_t_inv(r,r));
                    logL=logL+(-2*sum(log(diag(c))));
                    logL=logL+innovation(Lt,t-1)'*sigma_t_inv*innovation(Lt,t-1);
                end
            end
            logL=-logL/2;
        end
        
        function [zk_f,zk_u,Pk_f,Pk_u,J,logL] = Kfilter_parallel(Y,X_rg,X_beta,X_time,X_g,beta,G,sigma_eta,sigma_W_r,sigma_W_g,sigma_eps,sigma_geo,aj_rg,aj_g,M,z0,P0,time_diagonal,time_steps,pathparallel,tapering,compute_logL)
            if nargin<20
                error('You have to provide all the input arguments');
            end
                        
            if size(X_beta,2)~=length(beta)
                error('X_beta and beta are not conformable');
            end
           
            if size(G,1)~=size(G,2)
                error('G must be square');
            end
            
            if size(sigma_eta,1)~=size(sigma_eta,2)
                error('sigma_eta must be square');
            end      
            
            if size(G,1)~=size(sigma_eta,1)
                error('G and sigma_eta must have the same dimensions');
            end
            
            if size(Y,1)~=size(sigma_eps,1)
                error('The dimensions of sigma_eps must be equal to the number of rows of Y');
            end
            
            if size(z0,1)~=size(G,1)
                error('The length of z0 must be equal to the dimensions of G');
            end
            
            if size(P0,1)~=size(P0,2)
                error('P0 must be square');
            end
            
            if size(P0,1)~=size(G,1)
                error('The dimensions of P0 must be equal to the dimensions of G and sigma_eps');
            end
            
            if isempty(sigma_geo)
                compute_sigma_geo=1;
            else
                compute_sigma_geo=0;
            end
            
            time_steps=time_steps+1; %!!!
            
            min_ts=min(time_steps);
            max_ts=max(time_steps);
            if min_ts==2 %note the 2 due to +1
                server=1;
                %if min_ts==2 it means the call of the function is on the server
            else
                server=0;
            end
            
            p=size(G,1);
            N=size(Y,1);
            T=size(Y,2);
            zk_f=zeros(p,T+1);
            zk_u=zeros(p,T+1);
            Pk_f=zeros(p,p,T+1);
            Pk_u=zeros(p,p,T+1);
            J=zeros(p,size(Y,1),T+1);
            innovation=zeros(size(Y,1),T);
            
            zk_u(:,1)=z0;
            Pk_u(:,:,1)=P0;
            logL=0;
         
            if server
                for t=2:T+1
                    if size(X_time,3)==1
                        tK=2;
                    else
                        tK=t; %time variant
                    end
                    
                    if compute_sigma_geo
                        if not(isempty(X_rg))
                            sigma_geo=zeros(N);
                            if size(X_rg,3)>1
                                sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_r,M,'b'),X_rg(:,1,t-1),'b'),aj_rg,'b');
                            else
                                sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_r,M,'b'),X_rg(:,1,1),'b'),aj_rg,'b');
                            end
                        end
                        
                        if not(isempty(X_g))
                            if isempty(X_rg)
                                if tapering
                                    sigma_geo=spalloc(size(sigma_W_g{1},1),size(sigma_W_g{1},1),nnz(sigma_W_g{1}));
                                else
                                    sigma_geo=zeros(N);
                                end
                            end
                            for k=1:size(X_g,4)
                                if size(X_g,3)>1
                                    sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_g{k},X_g(:,1,t-1,k),'b'),aj_g(:,k),'b');
                                else
                                    sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_g{k},X_g(:,1,1,k),'b'),aj_g(:,k),'b');
                                end
                            end
                        end
                        if isempty(X_g)&&isempty(X_rg)
                            sigma_geo=sigma_eps;
                        else
                            sigma_geo=sigma_geo+sigma_eps;
                        end
                    end
                    
                    if size(X_beta,3)==1
                        tX=2;
                    else
                        tX=t; %time variant
                    end
                    Lt=not(isnan(Y(:,t-1))); %note the t-1

                    temp=sigma_geo(Lt,Lt);
                    if tapering
                        if not(stem_misc.isdiagonal(temp))
                            r = symamd(temp);
                            c=chol(temp(r,r));
                            temp2=speye(sum(Lt));
                            temp3=full(stem_misc.chol_solve(c,temp2(r,:)));
                            sigma_geo_inv=zeros(size(temp3));
                            sigma_geo_inv(r,:)=temp3;
                            clear temp2
                            clear temp3
                        else
                            d=1./diag(temp);
                            sigma_geo_inv=sparse(1:length(d),1:length(d),d);
                        end
                    else
                        if not(stem_misc.isdiagonal(temp))
                            c=chol(sigma_geo(Lt,Lt));
                            sigma_geo_inv=stem_misc.chol_solve(c,eye(sum(Lt)));
                        else
                            sigma_geo_inv=diag(1./diag(temp));
                        end
                    end
                        
                    if t>max_ts
                        %wait for the proper file from the clients
                        exit=0;
                        disp(['        Waiting for kalman_output_',num2str(t)]);
                        while not(exit)
                            exit=exist([pathparallel,'kalman_ouput_',num2str(t),'.mat'],'file');
                        end
                        read=0;
                        while not(read)
                            try
                                load([pathparallel,'kalman_ouput_',num2str(t),'.mat']);
                                read=1;
                                disp(['        kalman_ouput_',num2str(t),' readed']);
                            catch
                            end
                            pause(0.05);
                        end
                        deleted=0;
                        while not(deleted)
                            try
                                delete([pathparallel,'kalman_ouput_',num2str(t),'.mat']);
                                deleted=1;
                                disp(['        kalman_ouput_',num2str(t),' deleted']);
                            catch
                            end
                            pause(0.05);
                        end
                    end
                    
                    if not(time_diagonal)
                        %filter
                        zk_f(:,t)=G*zk_u(:,t-1); %(6.19) Stoffer
                        Pk_f(:,:,t)=G*Pk_u(:,:,t-1)*G'+sigma_eta; %(6.20) Stoffer
                        
                        %update
                        if t<=max_ts %the time steps up to max_ts are computed locally
                            temp=X_time(Lt,:,tK-1)'*sigma_geo_inv;
                            temp2=temp*X_time(Lt,:,tK-1);
                            temp3=Pk_f(:,:,t)*X_time(Lt,:,tK-1)';
                            J(:,Lt,t)=Pk_f(:,:,t)*temp-temp3*(temp'/(Pk_f(:,:,t)\eye(size(temp2))+temp2)*temp);
                        else
                            %temp and temp2 has been already reader from the file
                            temp3=Pk_f(:,:,t)*X_time(Lt,:,tK-1)';
                            J(:,Lt,t)=Pk_f(:,:,t)*temp-temp3*(temp'/(Pk_f(:,:,t)\eye(size(temp2))+temp2)*temp);
                        end
                        if compute_logL
                            sigma_t_inv=sigma_geo_inv-(temp'/((Pk_f(:,:,t)\eye(size(temp2)))+temp2))*temp;
                        end

                        if not(isempty(X_beta))
                            innovation(Lt,t-1)=Y(Lt,t-1)-X_beta(Lt,:,tX-1)*beta-X_time(Lt,:,tK-1)*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                        else
                            innovation(Lt,t-1)=Y(Lt,t-1)-X_time(Lt,:,tK-1)*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                        end
                        
                        zk_u(:,t)=zk_f(:,t)+J(:,Lt,t)*innovation(Lt,t-1);
                        Pk_u(:,:,t)=(eye(p)-J(:,Lt,t)*X_time(Lt,:,tK-1))*Pk_f(:,:,t);  %(6.22) Stoffer
                    else
                        %filter
                        zk_f(:,t)=diag(G).*zk_u(:,t-1); %(6.19) Stoffer
                        Pk_f(:,:,t)=diag(diag(G).^2.*diag(Pk_u(:,:,t-1))+diag(sigma_eta)); %(6.20) Stoffer
                        
                        %update
                        
                        %original formula
                        %J(i,Lt,t)=Pk_f(i,i,t)*X_time(Lt,i,tK-1)'/(X_time(Lt,i,tK-1)*Pk_f(i,i,t)*X_time(Lt,i,tK-1)'+sigma_geo(Lt,Lt)); %(6.23) Stoffer
                        %Sherman-Morrison-Woodbury formula: (B*P*B+D)^-1=D^-1-D^-1*B(P^-1+B*D^-1*B)^-1*B*D^-1
                        %J(i,Lt,t)=Pk_f(i,i,t)*X_time(Lt,i,tK-1)'*(sigma_geo_inv-sigma_geo_inv*X_time(Lt,i,tK-1)/(1/Pk_f(i,i,t)+X_time(Lt,i,tK-1)'*sigma_geo_inv*X_time(Lt,i,tK-1))*(X_time(Lt,i,tK-1)'*sigma_geo_inv));
                        
                        if t<=max_ts %the time steps up to max_ts are computed locally
                            temp=X_time(Lt,:,tK-1)'*sigma_geo_inv;
                            temp2=temp*X_time(Lt,:,tK-1);
                        else
                            %temp and temp2 has been already reader from the file
                        end
                        P=diag(1./diag(Pk_f(:,:,t)));
                        if compute_logL
                            sigma_t_inv=sigma_geo_inv-(temp'/(P+temp2))*temp;
                        end
                        temp3=Pk_f(:,:,t)*X_time(Lt,:,tK-1)';
                        J(:,Lt,t)=temp3*sigma_geo_inv-temp3*(temp'/(P+temp2)*temp);      
                        
                        if not(isempty(X_beta))
                            innovation(Lt,t-1)=Y(Lt,t-1)-X_beta(Lt,:,tX-1)*beta-X_time(Lt,:,tK-1)*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                        else
                            innovation(Lt,t-1)=Y(Lt,t-1)-X_time(Lt,:,tK-1)*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                        end
                        
                        zk_u(:,t)=zk_f(:,t)+J(:,Lt,t)*innovation(Lt,t-1);
                        Pk_u(:,:,t)=diag(diag((eye(p)-J(:,Lt,t)*X_time(Lt,:,tK-1))).*diag(Pk_f(:,:,t))); %(6.22) Stoffer
                    end
                    if compute_logL
                        r = symamd(sigma_t_inv);
                        c=chol(sigma_t_inv(r,r));
                        logL=logL+1/(2*sum(log(diag(c))));
                        logL=logL+innovation(Lt,t-1)'*sigma_t_inv*innovation(Lt,t-1);
                    end
                end
                logL=-logL/2;
            else
                %client computation
                for t=time_steps
                    if size(X_time,3)==1
                        tK=2;
                    else
                        tK=t; %time variant
                    end
                    if compute_sigma_geo
                        if not(isempty(X_rg))
                            sigma_geo=zeros(N);
                            if size(X_rg,3)>1
                                sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_r,M,'b'),X_rg(:,1,t-1),'b'),aj_rg,'b');
                            else
                                sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_r,M,'b'),X_rg(:,1,1),'b'),aj_rg,'b');
                            end
                        end
                        
                        if not(isempty(X_g))
                            if isempty(X_rg)
                                if tapering
                                    sigma_geo=spalloc(size(sigma_W_g{1},1),size(sigma_W_g{1},1),nnz(sigma_W_g{1}));
                                else
                                    sigma_geo=zeros(N);
                                end
                            end
                            for k=1:size(X_g,4)
                                if size(X_g,3)>1
                                    sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_g{k},X_g(:,1,t-1,k),'b'),aj_g(:,k),'b');
                                else
                                    sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_g{k},X_g(:,1,1,k),'b'),aj_g(:,k),'b');
                                end
                            end
                        end
                        if isempty(X_g)&&isempty(X_rg)
                            sigma_geo=sigma_eps;
                        else
                            sigma_geo=sigma_geo+sigma_eps;
                        end
                    end

                    Lt=not(isnan(Y(:,t-1))); %note the t-1
                    
                    temp=sigma_geo(Lt,Lt);
                    if tapering
                        if not(stem_misc.isdiagonal(temp))
                            r = symamd(temp);
                            c=chol(temp(r,r));
                            temp2=speye(sum(Lt));
                            temp3=full(stem_misc.chol_solve(c,temp2(r,:)));
                            sigma_geo_inv=zeros(size(temp3));
                            sigma_geo_inv(r,:)=temp3;
                            clear temp2
                            clear temp3
                        else
                            d=1./diag(temp);
                            sigma_geo_inv=sparse(1:length(d),1:length(d),d);
                        end
                    else
                        if not(stem_misc.isdiagonal(temp))
                            c=chol(sigma_geo(Lt,Lt));
                            sigma_geo_inv=stem_misc.chol_solve(c,eye(sum(Lt)));
                        else
                            sigma_geo_inv=diag(1./diag(temp));
                        end
                    end
                    
                    temp=X_time(Lt,:,tK-1)'*sigma_geo_inv;
                    temp2=temp*X_time(Lt,:,tK-1);
                    save([pathparallel,'temp/kalman_ouput_',num2str(t),'.mat'],'temp','temp2');
                    movefile([pathparallel,'temp/kalman_ouput_',num2str(t),'.mat'],[pathparallel,'kalman_ouput_',num2str(t),'.mat']);
                    disp(['Saved kalman_ouput_',num2str(t),'.mat']);
                end
            end
        end
        
        function [zk_s,Pk_s,PPk_s,logL] = Ksmoother(Y,X_rg,X_beta,X_time,X_g,beta,G,sigma_eta,sigma_W_r,sigma_W_g,sigma_eps,sigma_geo,aj_rg,aj_g,M,z0,P0,time_diagonal,time_steps,pathparallel,tapering,compute_logL)
            if isempty(pathparallel)
                [zk_f,zk_u,Pk_f,Pk_u,J,logL] = stem_kalman.Kfilter(Y,X_rg,X_beta,X_time,X_g,beta,G,sigma_eta,sigma_W_r,sigma_W_g,sigma_eps,sigma_geo,aj_rg,aj_g,M,z0,P0,time_diagonal,tapering,compute_logL);
            else
                [zk_f,zk_u,Pk_f,Pk_u,J,logL] = stem_kalman.Kfilter_parallel(Y,X_rg,X_beta,X_time,X_g,beta,G,sigma_eta,sigma_W_r,sigma_W_g,sigma_eps,sigma_geo,aj_rg,aj_g,M,z0,P0,time_diagonal,time_steps,pathparallel,tapering,compute_logL);
            end
            
            p=size(G,1);
            T=size(Y,2);
            
            H=zeros(p,p,T+1);
            Pk_s=zeros(p,p,T+1);
            zk_s=zeros(p,T+1);
            zk_s(:,end)=zk_u(:,end); %inizializzazione (6.47) Stoffer
            Pk_s(:,:,end)=Pk_u(:,:,end); %inizializzazione (6.48) Stoffer
            PPk_s=zeros(p,p,T+1);
            
            for t=T+1:-1:2
                H(:,:,t-1)=Pk_u(:,:,t-1)*G'/(Pk_f(:,:,t)); %(6.49) Stoffer
                zk_s(:,t-1)=zk_u(:,t-1)+H(:,:,t-1)*(zk_s(:,t)-zk_f(:,t)); %(6.47) Stoffer
                Pk_s(:,:,t-1)=Pk_u(:,:,t-1)+H(:,:,t-1)*(Pk_s(:,:,t)-Pk_f(:,:,t))*H(:,:,t-1)'; %(6.48) Stoffer
            end
            
            Lt=not(isnan(Y(:,end)));
            PPk_s(:,:,end)=(eye(p)-J(:,Lt,end)*X_time(Lt,:,end))*G*Pk_u(:,:,end-1); %(6.55) Stoffer
            for t=T+1:-1:3
                PPk_s(:,:,t-1)=Pk_u(:,:,t-1)*H(:,:,t-2)'+H(:,:,t-1)*(PPk_s(:,:,t)-G*Pk_u(:,:,t-1))*H(:,:,t-2)'; %(6.56) Stoffer
            end
        end
    end
end

