%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef stem_EM < EM
   
    properties
        stem_model=[];               %stem_model object
        stem_EM_options=[];
    end
    
    methods
        function obj = stem_EM(stem_model,stem_EM_options)
            % constructor
            % see properties description
            if nargin<2
                error('All the arguments must be provided');
            end

            if isa(stem_model,'stem_model')
                obj.stem_model=stem_model;
            else
                error('The first argument must be of class stem_model');
            end
            
            if isa(stem_EM_options,'stem_EM_options');
                obj.stem_EM_options=stem_EM_options;
            else
                error('The second argument must be of class stem_EM_options');
            end

            if isempty(obj.stem_model.stem_par_initial)
                error('Initial value estimation for model parameters must be provided first');
            end            
        end
        
        function st_EM_result = estimate(obj)
            % EM estimation
            t1=clock;
            if isempty(obj.stem_model)&&(nargin==0)
                error('You have to set the stem_model property first');
            end
            if isempty(obj.stem_model.stem_par_initial)
                error('Initial value estimation for model parameters must be provided first');
            end
            delta=9999;
            delta_logL=9999;
            last_logL=0;
            last_stem_par=obj.stem_model.stem_par;
            iteration=0;
            st_EM_result=stem_EM_result(); 
            st_EM_result.max_iterations=obj.stem_EM_options.max_iterations;
            st_EM_result.exit_toll=obj.stem_EM_options.exit_toll;
            st_EM_result.machine=computer;
            st_EM_result.date_start=datestr(now);
            while (delta>obj.stem_EM_options.exit_toll)&&(delta_logL>obj.stem_EM_options.exit_toll)&&(iteration<obj.stem_EM_options.max_iterations)
                ct1=clock;
                iteration=iteration+1;
                disp('************************');
                disp(['Iteration ',num2str(iteration),' started...']);
                disp('************************');
                
                clear E_wr_y1
                clear sum_Var_wr_y1
                clear diag_Var_wr_y1
                clear cov_wr_z_y1
                clear E_wg_y1 
                clear sum_Var_wg_y1
                clear diag_Var_wg_y1
                clear cov_wg_z_y1
                clear M_cov_wr_wg_y1
                clear cov_wgk_wgh_y1
                clear diag_Var_e_y1
                clear E_e_y1
                clear sigma_eps
                clear Xbeta

                [E_wr_y1,sum_Var_wr_y1,diag_Var_wr_y1,cov_wr_z_y1,E_wg_y1,sum_Var_wg_y1,diag_Var_wg_y1,cov_wg_z_y1,M_cov_wr_wg_y1,cov_wgk_wgh_y1,diag_Var_e_y1,E_e_y1,sigma_eps,sigma_W_r,sigma_W_g,Xbeta,st_kalmansmoother_result] = obj.E_step();
                obj.M_step(E_wr_y1,sum_Var_wr_y1,diag_Var_wr_y1,cov_wr_z_y1,E_wg_y1,sum_Var_wg_y1,diag_Var_wg_y1,cov_wg_z_y1,M_cov_wr_wg_y1,cov_wgk_wgh_y1,diag_Var_e_y1,E_e_y1,sigma_eps,sigma_W_r,sigma_W_g,st_kalmansmoother_result);

                st_EM_result.stem_par_all(:,iteration)=obj.stem_model.stem_par.vec;
                if not(isempty(st_kalmansmoother_result))
                    if not(st_kalmansmoother_result.logL==0)
                        logL=st_kalmansmoother_result.logL;
                        st_EM_result.logL_all(iteration)=logL;
                        delta_logL=abs(logL-last_logL)/abs(logL);
                        last_logL=logL;
                        disp('****************');
                        disp( ['logL: ',num2str(logL)]);
                        disp(['relative delta logL: ',num2str(delta_logL)]);
                    else
                       delta_logL=9999; 
                    end
                else
                    delta_logL=9999;
                end
                delta=norm(obj.stem_model.stem_par.vec()-last_stem_par.vec())/norm(last_stem_par.vec());
                last_stem_par=obj.stem_model.stem_par;
                disp(['Norm: ',num2str(delta)]);
                obj.stem_model.stem_par.print;
                ct2=clock;
                disp('**********************************************');
                disp(['Iteration ',num2str(iteration),' ended in ',stem_time(etime(ct2,ct1))]);
                disp('**********************************************');
            end
            t2=clock;
            st_EM_result.stem_par=obj.stem_model.stem_par;
            st_EM_result.stem_kalmansmoother_result=st_kalmansmoother_result;
            st_EM_result.E_wg_y1=E_wg_y1;
            st_EM_result.Var_wg_y1=diag_Var_wg_y1;
            st_EM_result.y_hat=obj.stem_model.stem_data.Y;
            st_EM_result.y_hat(isnan(st_EM_result.y_hat))=0;
            st_EM_result.y_hat=st_EM_result.y_hat-E_e_y1;
            st_EM_result.iterations=iteration;
            st_EM_result.computation_time=etime(t2,t1);
        end
        
        function st_EM_result = estimate_parallel(obj,pathparallel)
            % EM estimation
            t1=clock;
            if isempty(obj.stem_model)&&(nargin==0)
                error('You have to set the stem_model property first');
            end
            if isempty(obj.stem_model.stem_par_initial)
                error('Initial value estimation for model parameters must be provided first');
            end
            
            T=obj.stem_model.T;
            K=obj.stem_model.stem_par.k;
            local_efficiency=1;
            delta=9999;
            delta_logL=9999;
            last_logL=0;
            last_stem_par=obj.stem_model.stem_par;
            iteration=0;
            st_EM_result=stem_EM_result(); 
            st_EM_result.max_iterations=obj.stem_EM_options.max_iterations;
            st_EM_result.exit_toll=obj.stem_EM_options.exit_toll;
            st_EM_result.machine=computer;
            st_EM_result.date_start=datestr(now);
            while (delta>obj.stem_EM_options.exit_toll)&&(delta_logL>obj.stem_EM_options.exit_toll)&&(iteration<obj.stem_EM_options.max_iterations)
                ct1_iteration=clock;
                iteration=iteration+1;
                disp('************************');
                disp(['Iteration ',num2str(iteration),' started...']);
                disp('************************');
                
                %repeat the E-step until no timeout occurs
                timeout=1;
                while timeout
                    %set the timeout to zero
                    timeout=0;
                    clear E_wr_y1
                    clear sum_Var_wr_y1
                    clear diag_Var_wr_y1
                    clear cov_wr_z_y1
                    clear E_wg_y1
                    clear sum_Var_wg_y1
                    clear diag_Var_wg_y1
                    clear cov_wg_z_y1
                    clear M_cov_wr_wg_y1
                    clear cov_wgk_wgh_y1
                    clear diag_Var_e_y1
                    clear E_e_y1
                    clear sigma_eps
                    
                    %delete all the file in the exchange directory
                    files=dir([pathparallel,'*.mat']);
                    for i=1:length(files)
                        delete([pathparallel,files(i).name]);
                    end
                    
                    %create the file for the whoishere request
                    disp('    Looking for distributed clients...');
                    whoishere.IDrequest=unifrnd(0,100000,1,1);

                    save([pathparallel,'temp/whoishere.mat'],'whoishere');
                    movefile([pathparallel,'temp/whoishere.mat'],[pathparallel,'whoishere.mat']);
                   
                    if iteration==1
                        hosts=[];
                    end
                    nhosts=length(hosts);
                    
                    %wait for the replies from the clients
                    wait1=clock;
                    exit=0;
                    while not(exit)
                        files=dir([pathparallel,'machine_*.*']);
                        for i=1:length(files)
                            try
                                load([pathparallel,files(i).name])
                                if machine.IDrequest==whoishere.IDrequest
                                    %check if the client is already in the hosts list
                                    idx=[];
                                    for j=1:nhosts
                                        if hosts(j).IPaddress==machine.IPaddress
                                            hosts(j).active=1;
                                            hosts(j).require_stemmodel=machine.require_stemmodel;
                                            idx=j;
                                        end
                                    end
                                    %if not, add the client
                                    if isempty(idx)
                                        nhosts=nhosts+1;
                                        hosts(nhosts).IPaddress=machine.IPaddress;
                                        hosts(nhosts).data_received=0;
                                        %the first time a client is added it has the efficiency of the server
                                        hosts(nhosts).efficiency=local_efficiency;
                                        hosts(nhosts).require_stemmodel=machine.require_stemmodel;
                                        hosts(nhosts).active=1;
                                    end
                                end
                            catch
                            end
                        end
                        wait2=clock;
                        if etime(wait2,wait1)>15 %20 seconds timeout
                            exit=1;
                        end
                        pause(0.1);
                    end
                    delete([pathparallel,'whoishere.mat']);
                    %check for inactive clients
                    idx=[];
                    for i=1:nhosts
                        if hosts(i).active==0
                            idx=[idx i];
                        end
                    end
                    if not(isempty(idx))
                        hosts(idx)=[];
                        nhosts=length(hosts);
                    end
                    
                    %if there is at least one client then distribute the st_model
                    if nhosts>=1
                        disp(['    ',num2str(nhosts),' parallel client(s) found']);
                        disp('    Saving st_model to distribute');
                        st_model=obj.stem_model;
                        for i=1:nhosts
                            if hosts(i).require_stemmodel
                                save([pathparallel,'temp/st_model_parallel_',num2str(hosts(i).IPaddress),'.mat'],'st_model');
                                movefile([pathparallel,'temp/st_model_parallel_',num2str(hosts(i).IPaddress),'.mat'],[pathparallel,'st_model_parallel_',num2str(hosts(i).IPaddress),'.mat']);
                            end
                        end
                    else
                        disp('    No clients found. Only the server is used');
                    end
                   
                    if nhosts>=1
                        %the st_par to be distributed is the same for all the clients
                        disp('    Saving st_par to distribute')
                        st_par=obj.stem_model.stem_par;
                        for i=1:nhosts
                            save([pathparallel,'temp/st_par_parallel_',num2str(hosts(i).IPaddress),'.mat'],'st_par');
                            movefile([pathparallel,'temp/st_par_parallel_',num2str(hosts(i).IPaddress),'.mat'],[pathparallel,'st_par_parallel_',num2str(hosts(i).IPaddress),'.mat']);
                        end
                        clear st_par
                    end
                    
                    veff=local_efficiency;
                    for i=1:nhosts
                        veff=[veff hosts(i).efficiency];
                    end
                    veff=veff/sum(veff);
                    veff=[0 cumsum(veff)];
                    %compute the time_steps for the server
                    time_steps=1:round(veff(2)*T);
                    disp(['    ',num2str(length(time_steps)),' time will be assigned to the server machine']);                    
                    
                    %Kalman smoother
                    if obj.stem_model.stem_par.p>0
                        %distribute the st_par and the data needed to the clients
                        if nhosts>=1
                            disp('    Saving the Kalman data structure to distribute')
                            %send the information for the computation of the parallel kalman
                            data.iteration=iteration;
                            for i=1:nhosts
                                %compute the time_steps for the clients
                                data.time_steps=round(veff(i+1)*T)+1:round(veff(i+2)*T);
                                disp(['    ',num2str(length(data.time_steps)),' time steps assigned to client ',num2str(hosts(i).IPaddress)]);
                                save([pathparallel,'temp/kalman_parallel_',num2str(hosts(i).IPaddress),'.mat'],'data');
                                movefile([pathparallel,'temp/kalman_parallel_',num2str(hosts(i).IPaddress),'.mat'],[pathparallel,'kalman_parallel_',num2str(hosts(i).IPaddress),'.mat']);
                            end
                            %local Kalman Smoother computation
                            disp('    Kalman smoother started...');
                            ct1=clock;
                            st_kalman=stem_kalman(obj.stem_model);
                            [st_kalmansmoother_result,sigma_eps,~,~,~,~,~,~,~] = st_kalman.smoother(obj.stem_EM_options.compute_logL_at_all_steps,time_steps,pathparallel);
                            ct2=clock;
                            disp(['    Kalman smoother ended in ',stem_time(etime(ct2,ct1))]);
                        else
                            %The computation is only local. The standard Kalman smoother is considered
                            disp('    Kalman smoother started...');
                            ct1=clock;
                            st_kalman=stem_kalman(obj.stem_model);
                            [st_kalmansmoother_result,sigma_eps,~,~,~,~,~,~,~] = st_kalman.smoother(obj.stem_EM_options.compute_logL_at_all_steps);
                            ct2=clock;
                            disp(['    Kalman smoother ended in ',stem_time(etime(ct2,ct1))]);
                            time_steps=1:T;
                        end
                    else
                        st_kalmansmoother_result=[];
                        %sigma_eps
                        d=[];
                        dim=obj.stem_model.stem_data.dim;
                        for i=1:obj.stem_model.stem_data.nvar
                            d=[d;repmat(obj.stem_model.stem_par.sigma_eps(i,i),dim(i),1)];
                        end
                        sigma_eps=diag(d);
                    end
                    
                    ct1_distributed=clock;
                    disp('    Saving the E-step data structure to distribute')
                    data.st_kalmansmoother_result=st_kalmansmoother_result;
                    data.iteration=iteration;
                    for i=1:nhosts
                        %compute the time_steps for the clients
                        data.time_steps=round(veff(i+1)*T)+1:round(veff(i+2)*T);
                        disp(['    ',num2str(length(data.time_steps)),' time steps assigned to client ',num2str(hosts(i).IPaddress)]);
                        save([pathparallel,'temp/data_parallel_',num2str(hosts(i).IPaddress),'.mat'],'data');
                        movefile([pathparallel,'temp/data_parallel_',num2str(hosts(i).IPaddress),'.mat'],[pathparallel,'data_parallel_',num2str(hosts(i).IPaddress),'.mat']);
                    end
                    clear data

                    %local E-step computation
                    ct1_local=clock;
                    [E_wr_y1,sum_Var_wr_y1,diag_Var_wr_y1,cov_wr_z_y1,E_wg_y1,sum_Var_wg_y1,diag_Var_wg_y1,cov_wg_z_y1,M_cov_wr_wg_y1,cov_wgk_wgh_y1,diag_Var_e_y1,E_e_y1,cb] = obj.E_step_parallel(time_steps,st_kalmansmoother_result);
                    ct2_local=clock;
                    %local_efficiency=cb/etime(ct2_local,ct1_local);
                    local_efficiency=length(time_steps)/etime(ct2_local,ct1_local);
                    disp(['    Local computation burden: ',num2str(cb)]);
                    disp(['    Local computation time: ',num2str(etime(ct2_local,ct1_local))]);
                    disp(['    Local efficiency: ',num2str(local_efficiency)]);
                    
                    if nhosts>=1
                        disp('    Waiting for the results from the client(s)...');
                        exit=0;
                        wait1=clock;
                        while not(exit)
                            files=dir([pathparallel,'output_*.*']);
                            for i=1:length(files)
                                ct2_distributed=clock;
                                load([pathparallel,files(i).name]);
                                disp(['    Received output file from client ',num2str(output.IPaddress)]);
                                if iteration==output.iteration
                                    idx=[];
                                    for j=1:nhosts
                                        if (hosts(j).IPaddress==output.IPaddress)&&(hosts(j).data_received==0)
                                            idx=j;
                                        end
                                    end
                                    if not(isempty(idx))
                                        disp('    The data from the client was expected');
                                        %hosts(idx).efficiency=output.cb/output.ct;
                                        hosts(idx).efficiency=length(output.time_steps)/output.ct;
                                        disp(['    Computational burden of client ',num2str(hosts(idx).IPaddress),': ',num2str(output.cb)]);
                                        disp(['    Computational time of client ',num2str(hosts(idx).IPaddress),': ',num2str(output.ct)]);
                                        disp(['    Efficiency of client ',num2str(hosts(idx).IPaddress),': ',num2str(hosts(idx).efficiency)]);
                                        tsteps=output.time_steps;
                                        if not(isempty(E_wr_y1))
                                            E_wr_y1(:,tsteps)=output.E_wr_y1;
                                        end
                                        if not(isempty(sum_Var_wr_y1))
                                            %the matrix is recomposed since only the upper triangular part is received
                                            sum_Var_wr_y1=sum_Var_wr_y1+output.sum_Var_wr_y1+triu(output.sum_Var_wr_y1,1)';
                                        end
                                        if not(isempty(diag_Var_wr_y1))
                                            diag_Var_wr_y1(:,tsteps)=output.diag_Var_wr_y1;
                                        end
                                        if not(isempty(cov_wr_z_y1))
                                            cov_wr_z_y1(:,:,tsteps)=output.cov_wr_z_y1;
                                        end
                                        if not(isempty(E_wg_y1))
                                            for k=1:K
                                                E_wg_y1(:,tsteps,k)=output.E_wg_y1(:,:,k);
                                            end
                                        end
                                        if not(isempty(sum_Var_wg_y1))
                                            for k=1:K
                                                %the matrix is recomposed since only the upper triangular part is received
                                                sum_Var_wg_y1{k}=sum_Var_wg_y1{k}+output.sum_Var_wg_y1{k}+triu(output.sum_Var_wg_y1{k},1)';
                                            end
                                        end
                                        if not(isempty(diag_Var_wg_y1))
                                            for k=1:K
                                                diag_Var_wg_y1(:,tsteps,k)=output.diag_Var_wg_y1(:,:,k);
                                            end
                                        end
                                        if not(isempty(cov_wg_z_y1))
                                            for k=1:K
                                                cov_wg_z_y1(:,:,tsteps,k)=output.cov_wg_z_y1(:,:,:,k);
                                            end
                                        end
                                        if not(isempty(M_cov_wr_wg_y1))
                                            for k=1:K
                                                M_cov_wr_wg_y1(:,tsteps,k)=output.M_cov_wr_wg_y1(:,:,k);
                                            end
                                        end
                                        if iscell(cov_wgk_wgh_y1)
                                            for h=1:K
                                                for k=h+1:K
                                                    cov_wgk_wgh_y1{k,h}(:,tsteps)=output.cov_wgk_wgh_y1{k,h};
                                                end
                                            end
                                        end
                                        diag_Var_e_y1(:,tsteps)=output.diag_Var_e_y1;
                                        E_e_y1(:,tsteps)=output.E_e_y1;
                                        
                                        hosts(idx).data_received=1;
                                        clear output
                                    else
                                        disp('    Something is wrong');
                                    end
                                    exit=1;
                                    for j=1:nhosts
                                        if hosts(j).data_received==0
                                            exit=0;
                                        end
                                    end
                                    if exit==1
                                        disp('    All the data from the client(s) have been collected');
                                    end
                                else
                                    disp('    The iteration within the output file does not match. The file is deleted');
                                end
                                deleted=0;
                                while not(deleted)
                                    try
                                        delete([pathparallel,files(i).name]);
                                        deleted=1;
                                    catch
                                    end
                                end
                            end
                            wait2=clock;
                            if etime(wait2,wait1)>72000 %two hours
                                disp('    Timeout');
                                timeout=1;
                                exit=1;
                            end
                            pause(0.02);
                        end
                    end
                    
                    for i=1:nhosts
                        hosts(i).active=0;
                        hosts(i).data_received=0;
                    end
                end
                
                clear data
                if (K<=1)
                    %run the non parallel version of the M-step
                    logL = obj.M_step(E_wr_y1,sum_Var_wr_y1,diag_Var_wr_y1,cov_wr_z_y1,E_wg_y1,sum_Var_wg_y1,diag_Var_wg_y1,cov_wg_z_y1,M_cov_wr_wg_y1,cov_wgk_wgh_y1,diag_Var_e_y1,E_e_y1,sigma_eps,st_kalmansmoother_result);
                    %send the message to the other machine that they don't have to run the M-step
                    for i=1:nhosts
                        data.iteration=iteration;
                        data.index=[];
                        save([pathparallel,'temp/data_parallel_mstep',num2str(hosts(i).IPaddress),'.mat'],'data');
                        movefile([pathparallel,'temp/data_parallel_mstep',num2str(hosts(i).IPaddress),'.mat'],[pathparallel,'data_parallel_mstep',num2str(hosts(i).IPaddress),'.mat']);
                    end
                else
                    step=ceil(K/(nhosts+1));
                    counter=1;
                    for i=1:nhosts+1
                        if (i==1)
                            index_local=counter:step+counter-1;
                        else
                            index{i-1}=counter:step+counter-1;
                            index{i-1}(index{i-1}>K)=[];
                        end
                        counter=counter+step;
                    end
                    %send the messages to the host
                    clear data
                    if issparse(sum_Var_wg_y1{1})
                        r=symamd(sum_Var_wg_y1{1});
                    else
                        r=[];
                    end
                    for i=1:nhosts
                        data.iteration=iteration;
                        data.index=index{i};
                        disp(['     Preparing M-step data for client ',num2str(hosts(i).IPaddress)]);
                        data.sum_Var_wg_y1=sum_Var_wg_y1(index{i});
                        data.E_wg_y1=E_wg_y1(:,:,index{i});
                        data.r=r;
                        disp(['     Sending M-step data to client ',num2str(hosts(i).IPaddress)]);
                        save([pathparallel,'temp/data_parallel_mstep',num2str(hosts(i).IPaddress),'.mat'],'data');
                        movefile([pathparallel,'temp/data_parallel_mstep',num2str(hosts(i).IPaddress),'.mat'],[pathparallel,'data_parallel_mstep',num2str(hosts(i).IPaddress),'.mat']);
                        disp(['     M-Step data sent.']);
                    end
                    %M-step locale
                    logL = obj.M_step_parallel(E_wr_y1,sum_Var_wr_y1,diag_Var_wr_y1,cov_wr_z_y1,E_wg_y1,sum_Var_wg_y1,diag_Var_wg_y1,cov_wg_z_y1,M_cov_wr_wg_y1,cov_wgk_wgh_y1,diag_Var_e_y1,E_e_y1,sigma_eps,st_kalmansmoother_result,index_local);
                end
                
                %Attende la ricezione dagli altri nodi
                if nhosts>0
                    disp(['     Wait for output_mstep from the client(s)']);
                    exit=0;
                    while not(exit)
                        files=dir([pathparallel,'output_mstep_*.*']);
                        for i=1:length(files)
                            % try
                            ct2_distributed=clock;
                            load([pathparallel,files(i).name]);
                            disp(['    Received output_mstep file from client ',num2str(output.IPaddress)]);
                            if iteration==output.iteration
                                idx=[];
                                for j=1:nhosts
                                    if (hosts(j).IPaddress==output.IPaddress)&&(hosts(j).data_received==0)
                                        idx=j;
                                    end
                                end
                                if not(isempty(idx))
                                    disp('    The output_mstep from the client was expected');
                                    if not(isempty(output.index))
                                        for z=1:length(output.index)
                                            obj.stem_model.stem_par.v_g(:,:,output.index(z))=output.mstep_par.v_g(:,:,output.index(z));
                                            obj.stem_model.stem_par.theta_g(output.index(z))=output.mstep_par.theta_g(output.index(z));
                                            disp([num2str(output.index(z)),'th component of vg and theta_g updated']);
                                        end
                                    else
                                        disp('     The output_mstep data from the client is empty');
                                    end
                                    hosts(idx).data_received=1;
                                    clear output
                                else
                                    disp('    Something is wrong');
                                end
                                exit=1;
                                for j=1:nhosts
                                    if hosts(j).data_received==0
                                        exit=0;
                                    end
                                end
                                if exit==1
                                    disp('    All the M-step data from the client(s) have been collected');
                                end
                            else
                                disp('    The iteration within the output file does not match. The file is deleted');
                            end
                            deleted=0;
                            while not(deleted)
                                try
                                    delete([pathparallel,files(i).name]);
                                    deleted=1;
                                catch
                                end
                            end
                            %catch
                            %end
                        end
                        wait2=clock;
                        if etime(wait2,wait1)>72000 %two hours
                            disp('    Timeout');
                            timeout=1;
                            exit=1;
                        end
                        pause(0.05);
                    end
                    for i=1:nhosts
                        hosts(i).active=0;
                        hosts(i).data_received=0;
                    end
                end

                if not(isempty(st_kalmansmoother_result))
                    if not(st_kalmansmoother_result.logL==0)
                        logL=st_kalmansmoother_result.logL;
                        st_EM_result.logL_all(iteration)=logL;
                        delta_logL=abs(logL-last_logL)/abs(logL);
                        last_logL=logL;
                        disp('****************');
                        disp( ['logL: ',num2str(logL)]);
                        disp(['relative delta logL: ',num2str(delta_logL)]);
                    else
                       delta_logL=9999; 
                    end
                else
                    delta_logL=9999;
                end
                delta=norm(obj.stem_model.stem_par.vec()-last_stem_par.vec())/norm(last_stem_par.vec());
                last_stem_par=obj.stem_model.stem_par;
                disp(['Norm: ',num2str(delta)]);
                obj.stem_model.stem_par.print;
                ct2_iteration=clock;
                disp('**********************************************');
                disp(['Iteration ',num2str(iteration),' ended in ',stem_time(etime(ct2_iteration,ct1_iteration))]);
                disp('**********************************************');
            end
            t2=clock;
            st_EM_result.stem_par=obj.stem_model.stem_par;
            st_EM_result.stem_kalmansmoother_result=st_kalmansmoother_result;
            st_EM_result.E_wg_y1=E_wg_y1;
            st_EM_result.Var_wg_y1=diag_Var_wg_y1;
            st_EM_result.iterations=iteration;
            st_EM_result.computation_time=etime(t2,t1);
        end        
        
        function [E_wr_y1,sum_Var_wr_y1,diag_Var_wr_y1,cov_wr_z_y1,E_wg_y1,sum_Var_wg_y1,diag_Var_wg_y1,cov_wg_z_y1,M_cov_wr_wg_y1,cov_wgk_wgh_y1,diag_Var_e_y1,E_e_y1,sigma_eps,sigma_W_r,sigma_W_g,Xbeta,st_kalmansmoother_result] = E_step(obj,T)
            %INPUT
            %T: e' fornito solo nel caso di stima parallela e vengono
            %elaborati i time steps da 1 a T

            %Nota: dalla funzione E-step e' ricavata la funzione di kriging.
            %Le modifiche di questa funzione devono riperquotersi (quando
            %necessario) sulla funzione di kriging

            %E-Step
            if nargin==1
                T=obj.stem_model.stem_data.T;
            end
            N=obj.stem_model.stem_data.N;
            if not(isempty(obj.stem_model.stem_data.stem_varset_r))
                Nr=obj.stem_model.stem_data.stem_varset_r.N;
            else
                Nr=0;
            end
            Ng=obj.stem_model.stem_data.stem_varset_g.N;
            
            K=obj.stem_model.stem_par.k;
            p=obj.stem_model.stem_par.p;
            data=obj.stem_model.stem_data;
            par=obj.stem_model.stem_par;

            disp('  E step started...');
            ct1_estep=clock;
            
            if p>0
                %Kalman smoother
                disp('    Kalman smoother started...');
                ct1=clock;
                st_kalman=stem_kalman(obj.stem_model);
                [st_kalmansmoother_result,sigma_eps,sigma_W_r,sigma_W_g,sigma_Z,aj_rg,aj_g,M,sigma_geo] = st_kalman.smoother(obj.stem_EM_options.compute_logL_at_all_steps);
                ct2=clock;
                disp(['    Kalman smoother ended in ',stem_time(etime(ct2,ct1))]);
                if not(data.X_time_tv)
                    if obj.stem_model.tapering
                        %migliorare la creazione della matrice sparsa!!!
                        var_Zt=sparse(data.X_time(:,:,1))*sparse(sigma_Z)*sparse(data.X_time(:,:,1)'); 
                    else
                        var_Zt=data.X_time(:,:,1)*sigma_Z*data.X_time(:,:,1)';
                    end
                end
                if not(isempty(sigma_geo))
                    var_Yt=sigma_geo+var_Zt;
                end
            else
                disp('    sigma matrices evaluation started...');
                ct1=clock;
                [sigma_eps,sigma_W_r,sigma_W_g,sigma_geo,~,aj_rg,aj_g,M] = obj.stem_model.get_sigma();
                st_kalmansmoother_result=stem_kalmansmoother_result([],[],[]);
                var_Zt=[];
                %variance of Y
                if not(isempty(sigma_geo))
                    var_Yt=sigma_geo; %sigma_geo includes sigma_eps
                end
                ct2=clock;
                disp(['    sigma matrices evaluation ended in ',stem_time(etime(ct2,ct1))]);
            end
            
            res=data.Y;
            res(isnan(res))=0;
            if not(isempty(data.X_beta))
                disp('    Xbeta evaluation started...');
                ct1=clock;
                Xbeta=zeros(N,T);
                if data.X_beta_tv
                    for t=1:T
                        Xbeta(:,t)=data.X_beta(:,:,t)*par.beta;
                    end
                else
                    Xbeta=repmat(data.X_beta(:,:,1)*par.beta,1,T);
                end
                ct2=clock;
                disp(['    Xbeta evaluation ended in ',stem_time(etime(ct2,ct1))]);
                res=res-Xbeta;
            else
                Xbeta=[];
            end
            E_e_y1=res;
            diag_Var_e_y1=zeros(N,T);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Conditional expectation, conditional variance and conditional covariance evaluation  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %sigma_Z=Var(Zt)
            %var_Zt=Var(X_time*Zt*X_time')
            
            disp('    Conditional E, Var, Cov evaluation started...');
            ct1=clock;
            %cov_wr_yz time invariant case
            if not(isempty(data.X_rg))
                if (obj.stem_model.tapering)
                    Lr=find(sigma_W_r);
                    [Ir,Jr]=ind2sub(size(sigma_W_r),Lr);
                    nnz_r=length(Ir);
                end
                
                if not(data.X_rg_tv)
                    cov_wr_y=D_apply(D_apply(M_apply(sigma_W_r,M,'r'),data.X_rg(:,1,1),'r'),aj_rg,'r');
                end
                E_wr_y1=zeros(Nr,T);
                if (obj.stem_model.tapering)
                    sum_Var_wr_y1=spalloc(size(sigma_W_r,1),size(sigma_W_r,2),nnz_r);
                else
                    sum_Var_wr_y1=zeros(Nr);
                end
                diag_Var_wr_y1=zeros(Nr,T);
                cov_wr_z_y1=zeros(Nr,p,T);
            end
            %cov_wg_yz time invariant case
            if not(isempty(data.X_g))
                if obj.stem_model.tapering
                    Lg=find(sigma_W_g{1});
                    [Ig,Jg]=ind2sub(size(sigma_W_g{1}),Lg);
                    nnz_g=length(Ig);
                end
                if not(data.X_g_tv)
                    for k=1:K
                        cov_wg_y{k}=D_apply(D_apply(sigma_W_g{k},data.X_g(:,1,1,k),'r'),aj_g(:,k),'r');
                    end
                end
                for h=1:K
                    for k=h+1:K
                        %VERIFICARE SE PU� ESSERE SPARSA!!!
                        cov_wgk_wgh_y1{k,h}=zeros(Ng,T);
                    end
                end
                E_wg_y1=zeros(Ng,T,K);
                for k=1:K
                    if obj.stem_model.tapering
                        sum_Var_wg_y1{k}=spalloc(size(sigma_W_g{k},1),size(sigma_W_g{k},2),nnz_g);
                    else
                        sum_Var_wg_y1{k}=zeros(Ng,Ng);
                    end
                end
                diag_Var_wg_y1=zeros(Ng,T,K);
                cov_wg_z_y1=zeros(Ng,p,T,K);
            end
            
            if not(isempty(data.X_rg)) && not(isempty(data.X_g))
                M_cov_wr_wg_y1=zeros(N,T,K);
            else
                M_cov_wr_wg_y1=[];
            end
            
            for t=1:T
                %missing at time t
                Lt=not(isnan(data.Y(:,t)));
                
                if data.X_rg_tv
                    tRG=t;
                else
                    tRG=1;
                end
                if data.X_time_tv
                    tT=t;
                else
                    tT=1;
                end
                if data.X_g_tv
                    tG=t;
                else
                    tG=1;
                end
                
                %evaluate var_yt in the time variant case
                if data.X_tv
                    if not(isempty(data.X_rg))
                        sigma_geo=D_apply(D_apply(M_apply(sigma_W_r,M,'b'),data.X_rg(:,1,tRG),'b'),aj_rg,'b');
                    end
                    
                    if not(isempty(data.X_g))
                        if isempty(data.X_rg)
                            if obj.stem_model.tapering
                                sigma_geo=spalloc(size(sigma_W_g{1},1),size(sigma_W_g{1},1),nnz(sigma_W_g{1}));
                            else
                                sigma_geo=zeros(N);
                            end
                        end
                        for k=1:size(data.X_g,4)
                            sigma_geo=sigma_geo+D_apply(D_apply(sigma_W_g{k},data.X_g(:,1,tG,k),'b'),aj_g(:,k),'b');
                        end
                    end
                    if isempty(data.X_g)&&isempty(data.X_rg)
                        sigma_geo=sigma_eps;
                    else
                        sigma_geo=sigma_geo+sigma_eps;
                    end
                    
                    if not(isempty(data.X_time))
                        if not(isempty(data.X_rg))||not(isempty(data.X_g))
                            if data.X_time_tv
                                if obj.stem_model.tapering
                                    var_Zt=sparse(data.X_time(:,:,tT))*sparse(sigma_Z)*sparse(data.X_time(:,:,tT)');
                                else
                                    var_Zt=data.X_time(:,:,tT)*sigma_Z*data.X_time(:,:,tT)';
                                end
                            end
                            var_Yt=sigma_geo+var_Zt;
                        end
                    else
                        if not(isempty(data.X_rg))||not(isempty(data.X_g))
                            var_Yt=sigma_geo;
                        end
                    end
                end
                
                %check if the temporal loadings are time variant
                if not(isempty(data.X_time))
                    temp=data.X_time(:,:,tT)*st_kalmansmoother_result.Pk_s(:,:,t+1);
                    if N>obj.stem_model.system_size
                        blocks=0:80:size(diag_Var_e_y1,1);
                        if not(blocks(end)==size(diag_Var_e_y1,1))
                            blocks=[blocks size(diag_Var_e_y1,1)];
                        end
                        for i=1:length(blocks)-1
                            diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag(temp(blocks(i)+1:blocks(i+1),:)*data.X_time(blocks(i)+1:blocks(i+1),:,tT)');
                        end
                    else
                        diag_Var_e_y1(:,t)=diag(temp*data.X_time(:,:,tT)');
                    end
                    %update E(e|y1)
                    temp=st_kalmansmoother_result.zk_s(:,t+1);
                    E_e_y1(:,t)=E_e_y1(:,t)-data.X_time(:,:,tT)*temp;
                end
                
                if not(isempty(data.X_rg))||not(isempty(data.X_g))
                    %build the Ht matrix
                    if not(isempty(var_Zt))
                        H1t=[var_Yt(Lt,Lt), data.X_time(Lt,:,tT)*sigma_Z; sigma_Z*data.X_time(Lt,:,tT)', sigma_Z];
                    else
                        H1t=var_Yt(Lt,Lt);
                        temp=[];
                    end
                    
                    if obj.stem_model.tapering
                        cs=[];
                        r = symamd(H1t);
                        chol_H1t=chol(H1t(r,r));
                        temp2=[res(Lt,t);temp];
                        cs(r,1)=chol_solve(chol_H1t,temp2(r));
                    else
                        chol_H1t=chol(H1t);
                        cs=chol_solve(chol_H1t,[res(Lt,t);temp]);
                    end
                end
                
                if not(isempty(data.X_rg))
                    %check if the remote loadings are time variant
                    if data.X_rg_tv
                        %cov_wr_yz time variant case
                        cov_wr_y=D_apply(D_apply(M_apply(sigma_W_r,M,'r'),data.X_rg(:,1,tRG),'r'),aj_rg,'r');
                    end
                    cov_wr_y1z=[cov_wr_y(:,Lt),zeros(size(cov_wr_y,1),p)];
                    %compute E(w_r|y1);
                    E_wr_y1(:,t)=cov_wr_y1z*cs;
                    %compute Var(w_r|y1)
                    if obj.stem_model.tapering
                        temp_r(r,:)=chol_solve(chol_H1t,cov_wr_y1z(:,r)',1);
                        temp_r2=cov_wr_y1z*temp_r;
                        temp_r3=temp_r2(Lr);
                        temp_r3=sparse(Ir,Jr,temp_r3);
                        clear temp_r2
                        Var_wr_y1=sigma_W_r-temp_r3;
                        clear temp_r3
                    else
                        temp_r=chol_solve(chol_H1t,cov_wr_y1z');
                        Var_wr_y1=sigma_W_r-cov_wr_y1z*temp_r;
                    end
                    
                    if p>0
                        %compute cov(w_r,z|y1)
                        cov_wr_z_y1(:,:,t)=temp_r(end-p+1:end,:)'*st_kalmansmoother_result.Pk_s(:,:,t+1);
                        Var_wr_y1=Var_wr_y1+cov_wr_z_y1(:,:,t)*temp_r(end-p+1:end,:);
                        %update diag(Var(e|y1))
                        temp=D_apply(D_apply(M_apply(cov_wr_z_y1(:,:,t),M,'l'),data.X_rg(:,1,tRG),'l'),aj_rg,'l');
                        if N>obj.stem_model.system_size
                            blocks=0:80:size(diag_Var_e_y1,1);
                            if not(blocks(end)==size(diag_Var_e_y1,1))
                                blocks=[blocks size(diag_Var_e_y1,1)];
                            end
                            for i=1:length(blocks)-1
                                diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)+2*diag(temp(blocks(i)+1:blocks(i+1),:)*data.X_time(blocks(i)+1:blocks(i+1),:,tT)'); %notare 2*
                            end
                        else
                            %faster for N small
                            diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*diag(temp*data.X_time(:,:,tT)');
                        end
                    else
                        cov_wr_z_y1=[];
                    end
                    %compute diag(Var(w_r|y1))
                    diag_Var_wr_y1(:,t)=diag(Var_wr_y1);
                    %compute sum(Var(w_r|y1))
                    sum_Var_wr_y1=sum_Var_wr_y1+Var_wr_y1;
                    %update E(e|y1)
                    E_e_y1(:,t)=E_e_y1(:,t)-D_apply(D_apply(M_apply(E_wr_y1(:,t),M,'l'),data.X_rg(:,1,tRG),'l'),aj_rg,'l');
                    %update diag(Var(e|y1))
                    diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+D_apply(D_apply(M_apply(diag_Var_wr_y1(:,t),M,'l'),data.X_rg(:,1,tRG),'b'),aj_rg,'b'); %tested
                else
                    E_wr_y1=[];
                    diag_Var_wr_y1=[];
                    sum_Var_wr_y1=[];
                    cov_wr_z_y1=[];
                end
                clear temp_r
                if not(isempty(data.X_g))
                    %check if the ground loadings are time variant
                    if data.X_g_tv
                        %cov_wg_yz time invariant case
                        for k=1:K
                            cov_wg_y{k}=D_apply(D_apply(sigma_W_g{k},data.X_g(:,1,tG,k),'r'),aj_g(:,k),'r');
                        end
                    end
                    for k=1:K
                        cov_wg_y1z=[cov_wg_y{k}(:,Lt) zeros(size(cov_wg_y{k},1),p)];
                        %compute E(w_g_k|y1);
                        E_wg_y1(:,t,k)=cov_wg_y1z*cs;
                        %compute Var(w_g_k|y1)
                        if obj.stem_model.tapering
                            temp_g{k}(r,:)=chol_solve(chol_H1t,cov_wg_y1z(:,r)',1);
                            temp_g2=cov_wg_y1z*temp_g{k};
                            temp_g3=temp_g2(Lg);
                            temp_g3=sparse(Ig,Jg,temp_g3);
                            clear temp_g2
                            Var_wg_y1=sigma_W_g{k}-temp_g3;
                            clear temp_g3
                        else
                            temp_g{k}=chol_solve(chol_H1t,cov_wg_y1z');
                            Var_wg_y1=sigma_W_g{k}-cov_wg_y1z*temp_g{k};
                        end
                        
                        %VALUTARE SE ANCHE TEMP_G{K} VA TRIMMATA E APPORTARE LA MODIFICA AGLI ALTRI METODI!!!!
                        
                        if p>0
                            %compute cov(w_g,z|y1)
                            cov_wg_z_y1(:,:,t,k)=temp_g{k}(end-p+1:end,:)'*st_kalmansmoother_result.Pk_s(:,:,t+1);
                            Var_wg_y1=Var_wg_y1+cov_wg_z_y1(:,:,t,k)*temp_g{k}(end-p+1:end,:);
                            %update diag(Var(e|y1))
                            temp=D_apply(D_apply(cov_wg_z_y1(:,:,t,k),data.X_g(:,1,tG,k),'l'),aj_g(:,k),'l');
                            if N>obj.stem_model.system_size
                                blocks=0:80:size(diag_Var_e_y1,1);
                                if not(blocks(end)==size(diag_Var_e_y1,1))
                                    blocks=[blocks size(diag_Var_e_y1,1)];
                                end
                                for i=1:length(blocks)-1
                                    diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)+2*diag(temp(blocks(i)+1:blocks(i+1),:)*data.X_time(blocks(i)+1:blocks(i+1),:,tT)'); %notare 2*
                                end
                            else
                                diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*diag(temp*data.X_time(:,:,tT)');
                            end
                        else
                            cov_wg_z_y1=[];
                        end
                        diag_Var_wg_y1(:,t,k)=diag(Var_wg_y1);
                        sum_Var_wg_y1{k}=sum_Var_wg_y1{k}+Var_wg_y1;
                        %update E(e|y1)
                        E_e_y1(:,t)=E_e_y1(:,t)-D_apply(D_apply(E_wg_y1(:,t,k),data.X_g(:,1,tG,k),'l'),aj_g(:,k),'l');
                        %update diag(Var(e|y1))
                        diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+D_apply(D_apply(diag_Var_wg_y1(:,t,k),data.X_g(:,:,tG,k),'b'),aj_g(:,k),'b'); %K varianze
                        
                        if not(isempty(data.X_rg))
                            %compute M_cov(w_r,w_g|y1); cio� M*cov(w_r,w_g|y1) da tenere in considerazione nelle forme chiuse!
                            if length(M)>obj.stem_model.system_size
                                for i=1:length(M)
                                    %tested
                                    if p>0
                                        M_cov_wr_wg_y1(i,t,k)=-cov_wr_y1z(M(i),:)*temp_g{k}(:,i)+cov_wr_z_y1(M(i),:,t)*temp_g{k}(end-p+1:end,i); %ha gi� l'M_apply su left!!
                                    else
                                        M_cov_wr_wg_y1(i,t,k)=-cov_wr_y1z(M(i),:)*temp_g{k}(:,i);
                                    end
                                end
                            else
                                if p>0
                                    M_cov_wr_wg_y1(1:length(M),t,k)=diag(-cov_wr_y1z(M,:)*temp_g{k}(:,1:length(M))+cov_wr_z_y1(M,:,t)*temp_g{k}(end-p+1:end,1:length(M))); %ha gi� l'M_apply su left!!
                                else
                                    M_cov_wr_wg_y1(1:length(M),t,k)=diag(-cov_wr_y1z(M,:)*temp_g{k}(:,1:length(M)));
                                end
                            end
                            %update diag(Var(e|y1)) - tested
                            temp=D_apply(D_apply(M_cov_wr_wg_y1(:,t,k),data.X_rg(:,1,tRG),'l'),aj_rg,'l');
                            temp=D_apply(D_apply(temp,[data.X_g(:,1,tG,k);zeros(Nr,1)],'l'),aj_g(:,k),'l');
                            diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*temp;
                        end
                    end
                    
                    if K>1
                        %compute cov(w_gk,w_gh|y1);
                        for h=1:K
                            for k=h+1:K
                                cov_wgk_y1z=[cov_wg_y{k}(:,Lt) zeros(size(cov_wg_y{k},1),p)];
                                if N>obj.stem_model.system_size
                                    blocks=0:80:size(cov_wgk_y1z,1);
                                    if not(blocks(end)==size(cov_wgk_y1z,1))
                                        blocks=[blocks size(cov_wgk_y1z,1)];
                                    end
                                    for i=1:length(blocks)-1
                                        if not(isempty(cov_wg_z_y1))
                                            cov_wgk_wgh_y1{k,h}(blocks(i)+1:blocks(i+1),t)=diag(-cov_wgk_y1z(blocks(i)+1:blocks(i+1),:)*temp_g{h}(:,blocks(i)+1:blocks(i+1))+cov_wg_z_y1(blocks(i)+1:blocks(i+1),:,t,k)*temp_g{h}(end-p+1:end,blocks(i)+1:blocks(i+1)));
                                        else
                                            cov_wgk_wgh_y1{k,h}(blocks(i)+1:blocks(i+1),t)=diag(-cov_wgk_y1z(blocks(i)+1:blocks(i+1),:)*temp_g{h}(:,blocks(i)+1:blocks(i+1)));
                                        end
                                    end
                                else
                                    if not(isempty(cov_wg_z_y1))
                                        cov_wgk_wgh_y1{k,h}(:,t)=diag(-cov_wgk_y1z*temp_g{h}+cov_wg_z_y1(:,:,t,k)*temp_g{h}(end-p+1:end,:));
                                    else
                                        cov_wgk_wgh_y1{k,h}(:,t)=diag(-cov_wgk_y1z*temp_g{h});
                                    end
                                end
                                temp=D_apply(D_apply(cov_wgk_wgh_y1{k,h}(:,t),data.X_g(:,1,tG,k),'l'),aj_g(:,k),'l');
                                temp=D_apply(D_apply(temp,[data.X_g(:,1,tG,h);zeros(Nr,1)],'l'),aj_g(:,h),'l');
                                %update diag(Var(e|y1))
                                diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*temp;
                            end
                        end
                    else
                        cov_wgk_wgh_y1=[];
                    end
                else
                    E_wg_y1=[];
                    diag_Var_wg_y1=[];
                    sum_Var_wg_y1=[];
                    M_cov_wr_wg_y1=[];
                    cov_wg_z_y1=[];
                    cov_wgk_wgh_y1=[];
                end
                %delete the variables the dimension of which changes every t
                clear temp_g
            end
            ct2=clock;
            disp(['    Conditional E, Var, Cov evaluation ended in ',stem_time(etime(ct2,ct1))]);
            ct2_estep=clock;
            disp(['  E step ended in ',stem_time(etime(ct2_estep,ct1_estep))]);
            disp('');
        end
        
        function M_step(obj,E_wr_y1,sum_Var_wr_y1,diag_Var_wr_y1,cov_wr_z_y1,E_wg_y1,sum_Var_wg_y1,diag_Var_wg_y1,cov_wg_z_y1,M_cov_wr_wg_y1,cov_wgk_wgh_y1,diag_Var_e_y1,E_e_y1,sigma_eps,sigma_W_r,sigma_W_g,st_kalmansmoother_result)
            disp('  M step started...');
            ct1_mstep=clock;
            if not(isempty(obj.stem_model.stem_data.stem_varset_r))
                Nr=obj.stem_model.stem_data.stem_varset_r.N;
            else
                Nr=0;
            end
            Ng=obj.stem_model.stem_data.stem_varset_g.N;
            T=obj.stem_model.stem_data.T;
            N=obj.stem_model.stem_data.N;
            K=obj.stem_model.stem_par.k;
            M=obj.stem_model.stem_data.M;
            dim=obj.stem_model.stem_data.dim;
            
            data=obj.stem_model.stem_data;
            par=obj.stem_model.stem_par;
            st_par_em_step=par;
            
            inv_sigma_eps=diag(1./diag(sigma_eps));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             beta update                %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if not(isempty(data.X_beta))
                ct1=clock;
                disp('    beta update started...');
                temp1=zeros(size(data.X_beta,2));
                temp2=zeros(size(data.X_beta,2),1);
                d=diag(inv_sigma_eps);
                for t=1:T
                    Lt=not(isnan(data.Y(:,t)));
                    if data.X_beta_tv
                        tB=t;
                    else
                        tB=1;
                    end
                    temp1=temp1+data.X_beta(Lt,:,tB)'*D_apply(data.X_beta(Lt,:,tB),d(Lt),'l');
                    temp2=temp2+data.X_beta(Lt,:,tB)'*D_apply(E_e_y1(Lt,tB)+data.X_beta(Lt,:,tB)*par.beta,d(Lt),'l');
                end
                st_par_em_step.beta=temp1\temp2;
                ct2=clock;
                disp(['    beta update ended in ',stem_time(etime(ct2,ct1))]);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %              sigma_eps                 %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('    sigma_eps update started...');
            ct1=clock;
            temp=zeros(N,1);
            temp1=zeros(N,1);
            d=diag(sigma_eps);
            for t=1:T
                Lt=not(isnan(data.Y(:,t)));
                %the next two lines are ok only for sigma_eps diagonal
                temp1(Lt)=E_e_y1(Lt,t).^2+diag_Var_e_y1(Lt,t);
                temp1(~Lt)=d(~Lt);
                temp=temp+temp1;
            end
            temp=temp/T;
            blocks=[0 cumsum(dim)];
            for i=1:length(dim)
               st_par_em_step.sigma_eps(i,i)=mean(temp(blocks(i)+1:blocks(i+1)));
            end
            ct2=clock;
            disp(['    sigma_eps update ended in ',stem_time(etime(ct2,ct1))]);
           

            %%%%%%%%%%%%%%%%%%%%%%%%%
            %    G and sigma_eta    %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            if par.p>0
                disp('    G and sigma_eta update started...');
                ct1=clock;
                if not(obj.stem_model.stem_par.time_diagonal)
                    S11=st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,2:end)'+sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3);
                    S00=st_kalmansmoother_result.zk_s(:,1:end-1)*st_kalmansmoother_result.zk_s(:,1:end-1)'+sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3);
                    S10=st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,1:end-1)'+sum(st_kalmansmoother_result.PPk_s(:,:,2:end),3);
                else
                    S11=diag(diag(st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,2:end)'))+diag(diag(sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3)));
                    S00=diag(diag(st_kalmansmoother_result.zk_s(:,1:end-1)*st_kalmansmoother_result.zk_s(:,1:end-1)'))+diag(diag(sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3)));
                    S10=diag(diag(st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,1:end-1)'))+diag(diag(sum(st_kalmansmoother_result.PPk_s(:,:,2:end),3)));
                end
                
                temp=S10/S00;
                if max(eig(temp))<1
                    st_par_em_step.G=temp;
                else
                    warning('G is not stable. The last G is retained.');    
                end

                temp=(S11-S10*par.G'-par.G*S10'+par.G*S00*par.G')/T;
                %st_par_em_step.sigma_eta=(S11-S10*par.G'-par.G*S10'+par.G*S00*par.G')/T;
                %st_par_em_step.sigma_eta=(S11-st_par_em_step.G*S10')/T;
                if min(eig(temp))>0
                    st_par_em_step.sigma_eta=temp;
                else
                    warning('Sigma eta is not s.d.p. The last s.d.p. solution is retained');
                end
                ct2=clock;
                disp(['    G and sigma_eta ended in ',stem_time(etime(ct2,ct1))]);
            end
            
            [aj_rg,aj_g]=obj.stem_model.get_aj();
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %          alpha_rg, theta_r and v_r            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if not(isempty(data.X_rg))
                disp('    alpha_rg update started...');
                ct1=clock;
                for r=1:obj.stem_model.stem_data.nvar
                    [aj_rg_r,j_r] = obj.stem_model.get_jrg(r);
                    sum_num=0;
                    sum_den=0;
                    for t=1:T
                        if data.X_rg_tv
                            tRG=t;
                        else
                            tRG=1;
                        end
                        if data.X_time_tv
                            tT=t;
                        else
                            tT=1;
                        end
                        Lt=not(isnan(data.Y(:,t)));
                        temp1=E_e_y1(:,t)+D_apply(D_apply(M_apply(E_wr_y1(:,t),M,'l'),data.X_rg(:,1,tRG),'l'),aj_rg_r,'l');
                        temp2=D_apply(D_apply(M_apply(E_wr_y1(:,t)',M,'r'),data.X_rg(:,1,tRG),'r'),j_r,'r');
                        sum_num=sum_num+sum(temp1(Lt).*temp2(Lt)');
                        
                        if par.p>0
                            temp1=D_apply(D_apply(M_apply(cov_wr_z_y1(:,:,t),M,'l'),data.X_rg(:,1,tRG),'l'),j_r,'l');
                            temp2=zeros(size(temp1,1),1);
                            if N>obj.stem_model.system_size
                                blocks=0:80:size(temp1,1);
                                if not(blocks(end)==size(temp1,1))
                                    blocks=[blocks size(temp1,1)];
                                end
                                for i=1:length(blocks)-1
                                    temp2(blocks(i)+1:blocks(i+1))=temp1(blocks(i)+1:blocks(i+1),:)*data.X_time(blocks(i)+1:blocks(i+1),:,tT)';
                                end
                            else
                                temp2=temp1*data.X_time(:,:,tT)';
                            end
                            sum_num=sum_num-sum(temp2(Lt));
                        end
                        
                        if par.k>0
                            if data.X_g_tv
                                tG=t;
                            else
                                tG=1;
                            end
                            for k=1:K
                                temp1=D_apply(D_apply(M_cov_wr_wg_y1(:,t,k),data.X_rg(:,1,tRG),'l'),j_r,'l');
                                temp2=[data.X_g(:,1,tG,k);zeros(size(temp1,1)-size(data.X_g(:,1,tG,k),1),1)];
                                temp1=D_apply(D_apply(temp1',temp2,'r'),aj_g(:,k),'r');
                                sum_num=sum_num-sum(temp1(Lt));
                            end
                        end
                        
                        temp1=E_wr_y1(:,t).^2+diag_Var_wr_y1(:,t);
                        temp1=D_apply(D_apply(M_apply(temp1,M,'l'),data.X_rg(:,1,tRG),'b'),j_r,'b');
                        sum_den=sum_den+sum(temp1(Lt));
                    end
                    alpha_rg(r)=sum_num/sum_den;
                end
                st_par_em_step.alpha_rg=alpha_rg';
                ct2=clock;
                disp(['    alpha_rg update ended in ',stem_time(etime(ct2,ct1))]);

                disp('    v_r update started...');
                ct1=clock;
                
                if Nr<=obj.stem_EM_options.mstep_system_size
                    %TEMP FULL DEVE ESSERE CALCOLATO PER LE MATRICI V IN OGNI CASO SE NON SI RISOLVE COREGIONALIZZ!!!
                    temp=zeros(size(sum_Var_wr_y1));
                    for t=1:T
                        temp=temp+E_wr_y1(:,t)*E_wr_y1(:,t)';
                    end
                    temp=temp+sum_Var_wr_y1;
                end
                
                if par.remote_correlated
                    %indices are permutated in order to avoid deadlock
                    r = symamd(sum_Var_wr_y1); %note that sum_Var_wr_y1 and sigma_W_r have the same sparseness structure
                    kindex=randperm(size(par.v_r,1));
                    for k=kindex
                        hindex=randperm(size(par.v_r,1)-k)+k;
                        for h=hindex
                            initial=par.v_r(k,h);
                            if Nr<=obj.stem_EM_options.mstep_system_size
                                min_result = fminsearch(@(x) stem_geo_coreg_function_velement(x,k,h,par.v_r,par.theta_r,par.correlation_type,data.DistMat_r,...
                                    data.stem_varset_r.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_r.tap,r),initial,optimset('MaxIter',50,'TolX',1e-3,'UseParallel','always'));
                            else
                                disp('WARNING: this operation will take a long time');
                                min_result = fminsearch(@(x) stem_geo_coreg_function_velement(x,k,h,par.v_r,par.theta_r,par.correlation_type,data.DistMat_r,...
                                    data.stem_varset_r.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_r.tap,r),initial,optimset('MaxIter',50,'TolX',1e-3,'UseParallel','always'));
                            end
                            st_par_em_step.v_r(k,h)=min_result;
                            st_par_em_step.v_r(h,k)=min_result;
                        end
                    end
                    ct2=clock;
                    disp(['    v_r update ended in ',stem_time(etime(ct2,ct1))]);
                else
                    %da implementare
                end

                disp('    theta_r updating started...');
                ct1=clock;
                initial=par.theta_r;
                if Nr<=obj.stem_EM_options.mstep_system_size
                    min_result = fminsearch(@(x) stem_geo_coreg_function_theta(x,par.v_r,par.correlation_type,data.DistMat_r,...
                        data.stem_varset_r.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_r.tap,r),log(initial),optimset('MaxIter',50,'TolX',1e-3,'UseParallel','always'));
                    st_par_em_step.theta_r=exp(min_result);
                else
                    if data.stem_varset_r.nvar>1
                        disp('WARNING: this operation will take a long time');
                        min_result = fminsearch(@(x) stem_geo_coreg_function_theta(x,par.v_r,par.correlation_type,data.DistMat_r,...
                            data.stem_varset_r.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_r.tap,r),log(initial),optimset('MaxIter',50,'TolX',1e-3,'UseParallel','always'));
                        st_par_em_step.theta_r=exp(min_result);
                    else
                        s=ceil(Nr/obj.stem_EM_options.mstep_system_size);
                        step=ceil(Nr/s);
                        blocks=0:step:Nr;
                        if not(blocks(end)==Nr)
                            blocks=[blocks Nr];
                        end
                        for j=1:length(blocks)-1
                            block_size=blocks(j+1)-blocks(j);
                            idx=blocks(j)+1:blocks(j+1);
                            temp=zeros(block_size);
                            for t=1:T
                                temp=temp+E_wr_y1(idx,t)*E_wr_y1(idx,t)';
                            end
                            temp=temp+sum_Var_wr_y1(idx,idx);
                            r_partial=symamd(sum_Var_wr_y1(idx,idx));
                            min_result(j,:) = fminsearch(@(x) stem_geo_coreg_function_theta(x,par.v_r,par.correlation_type,data.DistMat_r(idx,idx),...
                                length(idx),temp,t,obj.stem_model.stem_data.stem_gridlist_r.tap,r_partial),log(initial),optimset('maxiter',50,'tolx',1e-3,'useparallel','always'));
                        end
                        st_par_em_step.theta_r=exp(mean(min_result));
                        %min_result = fminsearch(@(x) stem_geo_coreg_function_theta2(x,par.v_r,par.correlation_type,data.DistMat_r,...
                        %data.stem_varset_r.dim,E_wr_y1,sum_Var_wr_y1,obj.stem_model.stem_data.stem_gridlist_r.tap,r),initial,optimset('MaxIter',25,'TolX',1e-3,'UseParallel','always'));
                    end
                    ct2=clock;
                    disp(['    theta_r update ended in ',stem_time(etime(ct2,ct1))]);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %          alpha_g               %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if not(isempty(data.X_g))
                disp('    alpha_g update started...');
                for s=1:K
                    for r=1:obj.stem_model.stem_data.stem_varset_g.nvar
                        [aj_g_rs,j_r] = obj.stem_model.get_jg(r,s);
                        sum_num=0;
                        sum_den=0;
                        for t=1:T
                            if data.X_rg_tv
                                tRG=t;
                            else
                                tRG=1;
                            end
                            if data.X_time_tv
                                tT=t;
                            else
                                tT=1;
                            end
                            if data.X_g_tv
                                tG=t;
                            else
                                tG=1;
                            end
                            Lt=not(isnan(data.Y(:,t)));
                            
                            temp1=E_e_y1(:,t)+D_apply(D_apply(E_wg_y1(:,t,s),data.X_g(:,1,tG,s),'l'),aj_g_rs,'l');
                            temp2=D_apply(D_apply(E_wg_y1(:,t,s)',data.X_g(:,1,tG,s),'r'),j_r,'r');
                            sum_num=sum_num+sum(temp1(Lt).*temp2(Lt)');
                            
                            if par.p>0
                                temp1=D_apply(D_apply(cov_wg_z_y1(:,:,t,s),data.X_g(:,1,tG,s),'l'),j_r,'l');
                                temp2=zeros(size(temp1,1),1);
                                if N>obj.stem_model.system_size
                                    blocks=0:80:size(temp1,1);
                                    if not(blocks(end)==size(temp1,1))
                                        blocks=[blocks size(temp1,1)];
                                    end
                                    for i=1:length(blocks)-1
                                        temp2(blocks(i)+1:blocks(i+1))=temp1(blocks(i)+1:blocks(i+1),:)*data.X_time(blocks(i)+1:blocks(i+1),:,tT)';
                                    end
                                else
                                    temp2=temp1*data.X_time(:,:,tT)';
                                end
                                sum_num=sum_num-sum(temp2(Lt));
                            end
                            
                            if K>1
                                for k=1:K
                                    if not(k==s)
                                        if k<s
                                            kk=s;
                                            ss=k;
                                        else
                                            kk=k;
                                            ss=s;
                                        end
                                        temp1=D_apply(D_apply(cov_wgk_wgh_y1{kk,ss}(:,t),data.X_g(:,1,tG,k),'l'),aj_g(:,k),'l');
                                        temp1=D_apply(D_apply(temp1',[data.X_g(:,1,tG,s);zeros(Nr,1)],'r'),j_r,'r');
                                        sum_num=sum_num-sum(temp1(Lt));
                                    end
                                end
                            end
                            
                            if not(isempty(data.X_rg))
                                temp1=D_apply(D_apply(M_cov_wr_wg_y1(:,t,s),data.X_rg(:,1,tRG),'l'),aj_rg,'l');
                                temp2=[data.X_g(:,1,tG,s);zeros(size(temp1,1)-size(data.X_g(:,1,tG,s),1),1)];
                                temp1=D_apply(D_apply(temp1',temp2,'r'),j_r,'r');
                                sum_num=sum_num-sum(temp1(Lt));
                            end
                            
                            temp1=E_wg_y1(:,t,s).^2+diag_Var_wg_y1(:,t,s);
                            temp1=D_apply(D_apply(temp1,data.X_g(:,1,tG,s),'b'),j_r,'b');
                            sum_den=sum_den+sum(temp1(Lt));
                        end
                        alpha_g(r,s)=sum_num/sum_den;
                    end
                end
                st_par_em_step.alpha_g=alpha_g;
                ct2=clock;
                disp(['    alpha_g update ended in ',stem_time(etime(ct2,ct1))]);
                
                %indices are permutated in order to avoid deadlock
                disp('    v_g and theta_g update started...');
                ct1=clock;
                r = symamd(sum_Var_wg_y1{1}); %note that sum_Var_wg_y1{1} and sigma_W_g{k} have the same sparseness structure
                for z=1:K
                    if Ng<=obj.stem_EM_options.mstep_system_size
                        temp=zeros(size(sum_Var_wg_y1{z}));
                        for t=1:T
                            temp=temp+E_wg_y1(:,t,z)*E_wg_y1(:,t,z)';
                        end
                        temp=temp+sum_Var_wg_y1{z};
                    end
                    
                    kindex=randperm(size(par.v_g(:,:,z),1));
                    for k=kindex
                        hindex=randperm(size(par.v_g(:,:,z),1)-k)+k;
                        for h=hindex
                            initial=par.v_g(k,h,z);
                            if Ng<=obj.stem_EM_options.mstep_system_size
                                min_result = fminsearch(@(x) stem_geo_coreg_function_velement(x,k,h,par.v_g(:,:,z),par.theta_g(z),par.correlation_type,data.DistMat_g,...
                                    data.stem_varset_g.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_g.tap,r),initial,optimset('MaxIter',50,'TolX',1e-3));
                            else
                                disp('WARNING: this operation will take a long time');
                                min_result = fminsearch(@(x) stem_geo_coreg_function_velement(x,k,h,par.v_g(:,:,z),par.theta_g(z),par.correlation_type,data.DistMat_g,...
                                    data.stem_varset_g.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_g.tap,r),initial,optimset('MaxIter',50,'TolX',1e-3));
                            end
                            st_par_em_step.v_g(k,h,z)=min_result;
                            st_par_em_step.v_g(h,k,z)=min_result;
                        end
                    end
                    
                    initial=par.theta_g(z);
                    if Ng<=obj.stem_EM_options.mstep_system_size
                        min_result = fminsearch(@(x) stem_geo_coreg_function_theta(x,par.v_g(:,:,z),par.correlation_type,data.DistMat_g,...
                            data.stem_varset_g.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_g.tap,r),log(initial),optimset('MaxIter',50,'TolX',1e-3));
                        st_par_em_step.theta_g(z)=exp(min_result);
                    else
                        if data.stem_varset_g.nvar>1
                            disp('WARNING: this operation will take a long time');
                            min_result = fminsearch(@(x) stem_geo_coreg_function_theta(x,par.v_g(:,:,z),par.correlation_type,data.DistMat_g,...
                                data.stem_varset_g.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_g.tap,r),log(initial),optimset('MaxIter',50,'TolX',1e-3));
                            st_par_em_step.theta_g(z)=exp(min_result);
                        else
                            s=ceil(Ng/obj.stem_EM_options.mstep_system_size);
                            step=ceil(Ng/s);
                            blocks=0:step:Ng;
                            if not(blocks(end)==Ng)
                                blocks=[blocks Ng];
                            end
                            for j=1:length(blocks)-1
                                block_size=blocks(j+1)-blocks(j);
                                idx=blocks(j)+1:blocks(j+1);
                                temp=zeros(block_size);
                                for t=1:T
                                    temp=temp+E_wg_y1(idx,t,z)*E_wg_y1(idx,t,z)';
                                end
                                temp=temp+sum_Var_wg_y1{z}(idx,idx);
                                r_partial=symamd(sum_Var_wg_y1{z}(idx,idx));
                                min_result(j,:) = fminsearch(@(x) stem_geo_coreg_function_theta(x,par.v_g(:,:,z),par.correlation_type,data.DistMat_g(idx,idx),...
                                    length(idx),temp,t,obj.stem_model.stem_data.stem_gridlist_g.tap,r_partial),log(initial),optimset('maxiter',50,'tolx',1e-3,'useparallel','always'));
                            end
                            st_par_em_step.theta_g(z)=exp(mean(min_result));
                        end
                    end
                    ct2=clock;
                    disp(['    v_g and theta_g update ended in ',stem_time(etime(ct2,ct1))]);
                end
            end
            obj.stem_model.stem_par=st_par_em_step;
            ct2_mstep=clock;
            disp(['  M step ended in ',stem_time(etime(ct2_mstep,ct1_mstep))]);
        end
            
        function [E_wr_y1,sum_Var_wr_y1,diag_Var_wr_y1,cov_wr_z_y1,E_wg_y1,sum_Var_wg_y1,diag_Var_wg_y1,cov_wg_z_y1,M_cov_wr_wg_y1,cov_wgk_wgh_y1,diag_Var_e_y1,E_e_y1,cb] = E_step_parallel(obj,time_steps,st_kalmansmoother_result)
            %E-Step
            N=obj.stem_model.stem_data.N;
            if not(isempty(obj.stem_model.stem_data.stem_varset_r))
                Nr=obj.stem_model.stem_data.stem_varset_r.N;
            else
                Nr=0;
            end
            Ng=obj.stem_model.stem_data.stem_varset_g.N;
            T=length(time_steps);
            K=obj.stem_model.stem_par.k;
            p=obj.stem_model.stem_par.p;
            data=obj.stem_model.stem_data;
            par=obj.stem_model.stem_par;
            
            fts=time_steps(1);

            disp('  E step started...');
            ct1_estep=clock;
            
            [sigma_eps,sigma_W_r,sigma_W_g,sigma_geo,sigma_Z,aj_rg,aj_g,M] = obj.stem_model.get_sigma();
            if p>0
                if not(data.X_time_tv)
                    if obj.stem_model.tapering
                        var_Zt=sparse(data.X_time(:,:,1))*sparse(sigma_Z)*sparse(data.X_time(:,:,1)'); 
                    else
                        var_Zt=data.X_time(:,:,1)*sigma_Z*data.X_time(:,:,1)'; 
                    end
                end
                if not(isempty(sigma_geo))
                    var_Yt=sigma_geo+var_Zt;
                end                
            else
                st_kalmansmoother_result=stem_kalmansmoother_result([],[],[]);    
                var_Zt=[];
                %variance of Y
                if not(isempty(sigma_geo))
                    var_Yt=sigma_geo; %sigma_geo includes sigma_eps
                end                
            end            
            
            res=data.Y(:,time_steps);
            res(isnan(res))=0;
            if not(isempty(data.X_beta))
                disp('    Xbeta evaluation started...');
                ct1=clock;
                Xbeta=zeros(N,T);
                if data.X_beta_tv
                    for t=1:T
                        Xbeta(:,t)=data.X_beta(:,:,t+fts-1)*par.beta;
                    end
                else
                    Xbeta=repmat(data.X_beta(:,:,1)*par.beta,1,T);
                end
                ct2=clock;
                disp(['    Xbeta evaluation ended in ',stem_time(etime(ct2,ct1))]);
                res=res-Xbeta;
            else
                Xbeta=[];
            end
            E_e_y1=res;
            diag_Var_e_y1=zeros(N,T);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Conditional expectation, conditional variance and conditional covariance evaluation  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %sigma_Z=Var(Zt)
            %var_Zt=Var(X_time*Zt*X_time')
            
            disp('    Conditional E, Var, Cov evaluation started...');
            ct1=clock;
            %cov_wr_yz time invariant case
            if not(isempty(data.X_rg))
                if obj.stem_model.tapering
                    Lr=find(sigma_W_r);
                    [Ir,Jr]=ind2sub(size(sigma_W_r),Lr);
                    nnz_r=length(Ir);
                end
                if not(data.X_rg_tv)
                    cov_wr_y=D_apply(D_apply(M_apply(sigma_W_r,M,'r'),data.X_rg(:,1,1),'r'),aj_rg,'r');
                end
                E_wr_y1=zeros(Nr,T);
                
                if not(obj.stem_model.tapering)
                    sum_Var_wr_y1=zeros(Nr);
                else
                    sum_Var_wr_y1=spalloc(size(sigma_W_r,1),size(sigma_W_r,2),nnz_r);
                end
                
                diag_Var_wr_y1=zeros(Nr,T);
                cov_wr_z_y1=zeros(Nr,p,T);
            end
            %cov_wg_yz time invariant case
            if not(isempty(data.X_g))
                if obj.stem_model.tapering
                    Lg=find(sigma_W_g{1});
                    [Ig,Jg]=ind2sub(size(sigma_W_g{1}),Lg);
                    nnz_g=length(Ig);
                end
                if not(data.X_g_tv)
                    for k=1:K
                        cov_wg_y{k}=D_apply(D_apply(sigma_W_g{k},data.X_g(:,1,1,k),'r'),aj_g(:,k),'r');
                    end
                end
                for h=1:K
                    for k=h+1:K
                        cov_wgk_wgh_y1{k,h}=zeros(Ng,T);
                    end
                end
                E_wg_y1=zeros(Ng,T,K);
                for k=1:K
                    if not(obj.stem_model.tapering)
                        sum_Var_wg_y1{k}=zeros(Ng,Ng);
                    else
                        sum_Var_wg_y1{k}=spalloc(size(sigma_W_g{k},1),size(sigma_W_g{k},2),nnz_g);
                    end
                end
                diag_Var_wg_y1=zeros(Ng,T,K);
                cov_wg_z_y1=zeros(Ng,p,T,K);
            end
            
            if not(isempty(data.X_rg)) && not(isempty(data.X_g))
                M_cov_wr_wg_y1=zeros(N,T,K);
            else
                M_cov_wr_wg_y1=[];
            end
            
            cb=0; %computational burden
            for t=1:T
                t_partial1=clock;
                %missing at time t
                Lt=not(isnan(data.Y(:,t+fts-1)));
                
                if data.X_rg_tv
                    tRG=t+fts-1;
                else
                    tRG=1;
                end
                if data.X_time_tv
                    tT=t+fts-1;
                else
                    tT=1;
                end
                if data.X_g_tv
                    tG=t+fts-1;
                else
                    tG=1;
                end
                
                %evaluate var_yt in the time variant case
                if data.X_tv
                    if not(isempty(data.X_rg))
                        sigma_geo=D_apply(D_apply(M_apply(sigma_W_r,M,'b'),data.X_rg(:,1,tRG),'b'),aj_rg,'b');
                    end
                    
                    if not(isempty(data.X_g))
                        if isempty(data.X_rg)
                            if obj.stem_model.tapering
                                sigma_geo=spalloc(size(sigma_W_g{1},1),size(sigma_W_g{1},1),nnz(sigma_W_g{1}));
                            else
                                sigma_geo=zeros(N);
                            end
                        end
                        for k=1:size(data.X_g,4)
                            sigma_geo=sigma_geo+D_apply(D_apply(sigma_W_g{k},data.X_g(:,1,tG,k),'b'),aj_g(:,k),'b');
                        end
                    end
                    if isempty(data.X_g)&&isempty(data.X_rg)
                        sigma_geo=sigma_eps;
                    else
                        sigma_geo=sigma_geo+sigma_eps;
                    end
                    
                    if not(isempty(data.X_time))
                        if not(isempty(data.X_rg))||not(isempty(data.X_g))
                            if data.X_time_tv
                                if obj.stem_model.tapering
                                    var_Zt=sparse(data.X_time(:,:,tT))*sparse(sigma_Z)*sparse(data.X_time(:,:,tT)');
                                else
                                    var_Zt=data.X_time(:,:,tT)*sigma_Z*data.X_time(:,:,tT)';
                                end
                            end
                            var_Yt=sigma_geo+var_Zt;
                        end
                    else
                        if not(isempty(data.X_rg))||not(isempty(data.X_g))
                            var_Yt=sigma_geo;
                        end
                    end
                end
                
                %check if the temporal loadings are time variant
                if not(isempty(data.X_time))
                    temp=data.X_time(:,:,tT)*st_kalmansmoother_result.Pk_s(:,:,t+fts-1+1);
                    if N>obj.stem_model.system_size
                        blocks=0:80:size(diag_Var_e_y1,1);
                        if not(blocks(end)==size(diag_Var_e_y1,1))
                            blocks=[blocks size(diag_Var_e_y1,1)];
                        end
                        for i=1:length(blocks)-1
                            %update diag(Var(e|y1))
                            diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag(temp(blocks(i)+1:blocks(i+1),:)*data.X_time(blocks(i)+1:blocks(i+1),:,tT)');
                        end
                    else
                        diag_Var_e_y1(:,t)=diag(temp*data.X_time(:,:,tT)');
                    end
                end
                
                if not(isempty(data.X_rg))||not(isempty(data.X_g))
                    %build the Ht matrix
                    if not(isempty(var_Zt))
                        temp=st_kalmansmoother_result.zk_s(:,t+fts-1+1);
                        H1t=[var_Yt(Lt,Lt), data.X_time(Lt,:,tT)*sigma_Z; sigma_Z*data.X_time(Lt,:,tT)', sigma_Z];
                        %update E(e|y1)
                        E_e_y1(:,t)=E_e_y1(:,t)-data.X_time(:,:,tT)*temp;
                    else
                        H1t=var_Yt(Lt,Lt);
                        temp=[];
                    end
                    cb=cb+size(H1t,1)^3;
                    
                    if obj.stem_model.tapering
                        cs=[];
                        r = symamd(H1t);
                        chol_H1t=chol(H1t(r,r));
                        temp2=[res(Lt,t);temp];
                        cs(r,1)=chol_solve(chol_H1t,temp2(r));
                    else
                        chol_H1t=chol(H1t);
                        cs=chol_solve(chol_H1t,[res(Lt,t);temp]);
                    end
                end
                
                if not(isempty(data.X_rg))
                    %check if the remote loadings are time variant
                    if data.X_rg_tv
                        %cov_wr_yz time variant case
                        cov_wr_y=D_apply(D_apply(M_apply(sigma_W_r,M,'r'),data.X_rg(:,1,tRG),'r'),aj_rg,'r');
                    end
                    cov_wr_y1z=[cov_wr_y(:,Lt),zeros(size(cov_wr_y,1),p)];
                    %compute E(w_r|y1);
                    E_wr_y1(:,t)=cov_wr_y1z*cs;
                    %compute Var(w_r|y1)
                    if obj.stem_model.tapering
                        temp_r(r,:)=chol_solve(chol_H1t,cov_wr_y1z(:,r)',1);
                        temp_r2=cov_wr_y1z*temp_r;
                        temp_r3=temp_r2(Lr);
                        temp_r3=sparse(Ir,Jr,temp_r3);
                        clear temp_r2
                        Var_wr_y1=sigma_W_r-temp_r3;
                        clear temp_r3
                    else
                        temp_r=chol_solve(chol_H1t,cov_wr_y1z');
                        Var_wr_y1=sigma_W_r-cov_wr_y1z*temp_r;
                    end
                    
                    if p>0
                        %compute cov(w_r,z|y1)
                        cov_wr_z_y1(:,:,t)=temp_r(end-p+1:end,:)'*st_kalmansmoother_result.Pk_s(:,:,t+fts-1+1);
                        Var_wr_y1=Var_wr_y1+cov_wr_z_y1(:,:,t)*temp_r(end-p+1:end,:);
                        %update diag(Var(e|y1))
                        temp=D_apply(D_apply(M_apply(cov_wr_z_y1(:,:,t),M,'l'),data.X_rg(:,1,tRG),'l'),aj_rg,'l');
                        if N>obj.stem_model.system_size
                            blocks=0:80:size(diag_Var_e_y1,1);
                            if not(blocks(end)==size(diag_Var_e_y1,1))
                                blocks=[blocks size(diag_Var_e_y1,1)];
                            end
                            for i=1:length(blocks)-1
                                diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)+2*diag(temp(blocks(i)+1:blocks(i+1),:)*data.X_time(blocks(i)+1:blocks(i+1),:,tT)'); %notare 2*
                            end
                        else
                            %faster for N small
                            diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*diag(temp*data.X_time(:,:,tT)');
                        end
                    else
                        cov_wr_z_y1=[];
                    end
                    %compute diag(Var(w_r|y1))
                    diag_Var_wr_y1(:,t)=diag(Var_wr_y1);
                    %compute sum(Var(w_r|y1))
                    sum_Var_wr_y1=sum_Var_wr_y1+Var_wr_y1;
                    %update E(e|y1)
                    E_e_y1(:,t)=E_e_y1(:,t)-D_apply(D_apply(M_apply(E_wr_y1(:,t),M,'l'),data.X_rg(:,1,tRG),'l'),aj_rg,'l');
                    %update diag(Var(e|y1))
                    diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+D_apply(D_apply(M_apply(diag_Var_wr_y1(:,t),M,'l'),data.X_rg(:,1,tRG),'b'),aj_rg,'b'); %tested
                else
                    E_wr_y1=[];
                    diag_Var_wr_y1=[];
                    sum_Var_wr_y1=[];
                    cov_wr_z_y1=[];
                end
                clear temp_r
                if not(isempty(data.X_g))
                    %check if the ground loadings are time variant
                    if data.X_g_tv
                        %cov_wg_yz time invariant case
                        for k=1:K
                            cov_wg_y{k}=D_apply(D_apply(sigma_W_g{k},data.X_g(:,1,tG,k),'r'),aj_g(:,k),'r');
                        end
                    end
                    for k=1:K
                        cov_wg_y1z=[cov_wg_y{k}(:,Lt) zeros(size(cov_wg_y{k},1),p)];
                        %compute E(w_g_k|y1);
                        E_wg_y1(:,t,k)=cov_wg_y1z*cs;
                        %compute Var(w_g_k|y1)
                        if obj.stem_model.tapering
                            temp_g{k}(r,:)=chol_solve(chol_H1t,cov_wg_y1z(:,r)',1);
                            temp_g2=cov_wg_y1z*temp_g{k};
                            temp_g3=temp_g2(Lg);
                            temp_g3=sparse(Ig,Jg,temp_g3);
                            clear temp_g2
                            Var_wg_y1=sigma_W_g{k}-temp_g3;
                            clear temp_g3
                        else
                            temp_g{k}=chol_solve(chol_H1t,cov_wg_y1z');
                            Var_wg_y1=sigma_W_g{k}-cov_wg_y1z*temp_g{k};
                        end
                        
                        if p>0
                            %compute cov(w_g,z|y1)
                            cov_wg_z_y1(:,:,t,k)=temp_g{k}(end-p+1:end,:)'*st_kalmansmoother_result.Pk_s(:,:,t+fts-1+1);
                            Var_wg_y1=Var_wg_y1+cov_wg_z_y1(:,:,t,k)*temp_g{k}(end-p+1:end,:);
                            %update diag(Var(e|y1))
                            temp=D_apply(D_apply(cov_wg_z_y1(:,:,t,k),data.X_g(:,1,tG,k),'l'),aj_g(:,k),'l');
                            if N>obj.stem_model.system_size
                                blocks=0:80:size(diag_Var_e_y1,1);
                                if not(blocks(end)==size(diag_Var_e_y1,1))
                                    blocks=[blocks size(diag_Var_e_y1,1)];
                                end
                                for i=1:length(blocks)-1
                                    diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)+diag(2*temp(blocks(i)+1:blocks(i+1),:)*data.X_time(blocks(i)+1:blocks(i+1),:,tT)'); %notare 2*
                                end
                            else
                                diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*diag(temp*data.X_time(:,:,tT)');
                            end
                        else
                            cov_wg_z_y1=[];
                        end
                        diag_Var_wg_y1(:,t,k)=diag(Var_wg_y1);
                        sum_Var_wg_y1{k}=sum_Var_wg_y1{k}+Var_wg_y1;
                        %update E(e|y1)
                        E_e_y1(:,t)=E_e_y1(:,t)-D_apply(D_apply(E_wg_y1(:,t,k),data.X_g(:,1,tG,k),'l'),aj_g(:,k),'l');
                        %update diag(Var(e|y1))
                        diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+D_apply(D_apply(diag_Var_wg_y1(:,t,k),data.X_g(:,:,tG,k),'b'),aj_g(:,k),'b'); %K varianze
                        
                        if not(isempty(data.X_rg))
                            %compute M_cov(w_r,w_g|y1); cio� M*cov(w_r,w_g|y1) da tenere in considerazione nelle forme chiuse!
                            if length(M)>obj.stem_model.system_size
                                for i=1:length(M)
                                    %tested
                                    if p>0
                                        M_cov_wr_wg_y1(i,t,k)=-cov_wr_y1z(M(i),:)*temp_g{k}(:,i)+cov_wr_z_y1(M(i),:,t)*temp_g{k}(end-p+1:end,i); %ha gi� l'M_apply su left!!
                                    else
                                        M_cov_wr_wg_y1(i,t,k)=-cov_wr_y1z(M(i),:)*temp_g{k}(:,i);
                                    end
                                end
                            else
                                if p>0
                                    M_cov_wr_wg_y1(1:length(M),t,k)=diag(-cov_wr_y1z(M,:)*temp_g{k}(:,1:length(M))+cov_wr_z_y1(M,:,t)*temp_g{k}(end-p+1:end,1:length(M))); %ha gi� l'M_apply su left!!
                                else
                                    M_cov_wr_wg_y1(1:length(M),t,k)=diag(-cov_wr_y1z(M,:)*temp_g{k}(:,1:length(M)));
                                end
                            end
                            %update diag(Var(e|y1)) - tested
                            temp=D_apply(D_apply(M_cov_wr_wg_y1(:,t,k),data.X_rg(:,1,tRG),'l'),aj_rg,'l');
                            temp=D_apply(D_apply(temp,[data.X_g(:,1,tG,k);zeros(Nr,1)],'l'),aj_g(:,k),'l');
                            diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*temp;
                        end
                    end
                    
                    if K>1
                        %compute cov(w_gk,w_gh|y1);
                        for h=1:K
                            for k=h+1:K
                                cov_wgk_y1z=[cov_wg_y{k}(:,Lt) zeros(size(cov_wg_y{k},1),p)];
                                if N>obj.stem_model.system_size
                                    blocks=0:80:size(cov_wgk_y1z,1);
                                    if not(blocks(end)==size(cov_wgk_y1z,1))
                                        blocks=[blocks size(cov_wgk_y1z,1)];
                                    end
                                    for i=1:length(blocks)-1
                                        if not(isempty(cov_wg_z_y1))
                                            cov_wgk_wgh_y1{k,h}(blocks(i)+1:blocks(i+1),t)=diag(-cov_wgk_y1z(blocks(i)+1:blocks(i+1),:)*temp_g{h}(:,blocks(i)+1:blocks(i+1))+cov_wg_z_y1(blocks(i)+1:blocks(i+1),:,t,k)*temp_g{h}(end-p+1:end,blocks(i)+1:blocks(i+1)));
                                        else
                                            cov_wgk_wgh_y1{k,h}(blocks(i)+1:blocks(i+1),t)=diag(-cov_wgk_y1z(blocks(i)+1:blocks(i+1),:)*temp_g{h}(:,blocks(i)+1:blocks(i+1)));
                                        end
                                    end
                                else
                                    if not(isempty(cov_wg_z_y1))
                                        cov_wgk_wgh_y1{k,h}(:,t)=diag(-cov_wgk_y1z*temp_g{h}+cov_wg_z_y1(:,:,t,k)*temp_g{h}(end-p+1:end,:));
                                    else
                                        cov_wgk_wgh_y1{k,h}(:,t)=diag(-cov_wgk_y1z*temp_g{h});
                                    end
                                end
                                temp=D_apply(D_apply(cov_wgk_wgh_y1{k,h}(:,t),data.X_g(:,1,tG,k),'l'),aj_g(:,k),'l');
                                temp=D_apply(D_apply(temp,[data.X_g(:,1,tG,h);zeros(Nr,1)],'l'),aj_g(:,h),'l');
                                %update diag(Var(e|y1))
                                diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*temp;
                            end
                        end
                    else
                        cov_wgk_wgh_y1=[];
                    end
                else
                    E_wg_y1=[];
                    diag_Var_wg_y1=[];
                    sum_Var_wg_y1=[];
                    M_cov_wr_wg_y1=[];
                    cov_wg_z_y1=[];
                    cov_wgk_wgh_y1=[];
                end
                clear temp_g
                t_partial2=clock;
                disp(['      Time step ',num2str(t),' evaluated in ',stem_time(etime(t_partial2,t_partial1))]);
            end
            
            ct2=clock;
            disp(['    Conditional E, Var, Cov evaluation ended in ',stem_time(etime(ct2,ct1))]);
            ct2_estep=clock;
            disp(['  E step ended in ',stem_time(etime(ct2_estep,ct1_estep))]);
            disp('');
        end
        
        function M_step_parallel(obj,E_wr_y1,sum_Var_wr_y1,diag_Var_wr_y1,cov_wr_z_y1,E_wg_y1,sum_Var_wg_y1,diag_Var_wg_y1,cov_wg_z_y1,M_cov_wr_wg_y1,cov_wgk_wgh_y1,diag_Var_e_y1,E_e_y1,sigma_eps,st_kalmansmoother_result,index)
            disp('  M step started...');
            ct1_mstep=clock;
            if not(isempty(obj.stem_model.stem_data.stem_varset_r))
                Nr=obj.stem_model.stem_data.stem_varset_r.N;
            else
                Nr=0;
            end
            Ng=obj.stem_model.stem_data.stem_varset_g.N;
            T=obj.stem_model.stem_data.T;
            N=obj.stem_model.stem_data.N;
            K=obj.stem_model.stem_par.k;
            M=obj.stem_model.stem_data.M;
            dim=obj.stem_model.stem_data.dim;
            
            data=obj.stem_model.stem_data;
            par=obj.stem_model.stem_par;
            st_par_em_step=par;
            
            inv_sigma_eps=diag(1./diag(sigma_eps));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             beta update                %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if not(isempty(data.X_beta))
                ct1=clock;
                disp('    beta update started...');
                temp1=zeros(size(data.X_beta,2));
                temp2=zeros(size(data.X_beta,2),1);
                d=diag(inv_sigma_eps);
                for t=1:T
                    Lt=not(isnan(data.Y(:,t)));
                    if data.X_beta_tv
                        tB=t;
                    else
                        tB=1;
                    end
                    temp1=temp1+data.X_beta(Lt,:,tB)'*D_apply(data.X_beta(Lt,:,tB),d(Lt),'l');
                    temp2=temp2+data.X_beta(Lt,:,tB)'*D_apply(E_e_y1(Lt,tB)+data.X_beta(Lt,:,tB)*par.beta,d(Lt),'l');
                end
                st_par_em_step.beta=temp1\temp2;
                ct2=clock;
                disp(['    beta update ended in ',stem_time(etime(ct2,ct1))]);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %              sigma_eps                 %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('    sigma_eps update started...');
            ct1=clock;
            temp=zeros(N,1);
            temp1=zeros(N,1);
            d=diag(sigma_eps);
            for t=1:T
                Lt=not(isnan(data.Y(:,t)));
                %the next two lines are ok only for sigma_eps diagonal
                temp1(Lt)=E_e_y1(Lt,t).^2+diag_Var_e_y1(Lt,t);
                temp1(~Lt)=d(~Lt);
                temp=temp+temp1;
            end
            temp=temp/T;
            blocks=[0 cumsum(dim)];
            for i=1:length(dim)
               st_par_em_step.sigma_eps(i,i)=mean(temp(blocks(i)+1:blocks(i+1)));
            end
            ct2=clock;
            disp(['    sigma_eps update ended in ',stem_time(etime(ct2,ct1))]);
           

            %%%%%%%%%%%%%%%%%%%%%%%%%
            %    G and sigma_eta    %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            if par.p>0
                disp('    G and sigma_eta update started...');
                ct1=clock;
                if not(obj.stem_model.stem_par.time_diagonal)
                    S11=st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,2:end)'+sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3);
                    S00=st_kalmansmoother_result.zk_s(:,1:end-1)*st_kalmansmoother_result.zk_s(:,1:end-1)'+sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3);
                    S10=st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,1:end-1)'+sum(st_kalmansmoother_result.PPk_s(:,:,2:end),3);
                else
                    S11=diag(diag(st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,2:end)'))+diag(diag(sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3)));
                    S00=diag(diag(st_kalmansmoother_result.zk_s(:,1:end-1)*st_kalmansmoother_result.zk_s(:,1:end-1)'))+diag(diag(sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3)));
                    S10=diag(diag(st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,1:end-1)'))+diag(diag(sum(st_kalmansmoother_result.PPk_s(:,:,2:end),3)));
                end
                
                temp=S10/S00;
                if max(eig(temp))<1
                    st_par_em_step.G=temp;
                else
                    warning('G is not stable. The last G is retained.');    
                end

                temp=(S11-S10*par.G'-par.G*S10'+par.G*S00*par.G')/T;
                %st_par_em_step.sigma_eta=(S11-S10*par.G'-par.G*S10'+par.G*S00*par.G')/T;
                %st_par_em_step.sigma_eta=(S11-st_par_em_step.G*S10')/T;
                if min(eig(temp))>0
                    st_par_em_step.sigma_eta=temp;
                else
                    warning('Sigma eta is not s.d.p. The last s.d.p. solution is retained');
                end
                ct2=clock;
                disp(['    G and sigma_eta ended in ',stem_time(etime(ct2,ct1))]);
            end
            
            [aj_rg,aj_g]=obj.stem_model.get_aj();
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %          alpha_rg, theta_r and v_r            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if not(isempty(data.X_rg))
                disp('    alpha_rg update started...');
                ct1=clock;
                for r=1:obj.stem_model.stem_data.nvar
                    [aj_rg_r,j_r] = obj.stem_model.get_jrg(r);
                    sum_num=0;
                    sum_den=0;
                    for t=1:T
                        if data.X_rg_tv
                            tRG=t;
                        else
                            tRG=1;
                        end
                        if data.X_time_tv
                            tT=t;
                        else
                            tT=1;
                        end
                        Lt=not(isnan(data.Y(:,t)));
                        temp1=E_e_y1(:,t)+D_apply(D_apply(M_apply(E_wr_y1(:,t),M,'l'),data.X_rg(:,1,tRG),'l'),aj_rg_r,'l');
                        temp2=D_apply(D_apply(M_apply(E_wr_y1(:,t)',M,'r'),data.X_rg(:,1,tRG),'r'),j_r,'r');
                        sum_num=sum_num+sum(temp1(Lt).*temp2(Lt)');
                        
                        if par.p>0
                            temp1=D_apply(D_apply(M_apply(cov_wr_z_y1(:,:,t),M,'l'),data.X_rg(:,1,tRG),'l'),j_r,'l');
                            temp2=zeros(size(temp1,1),1);
                            if N>obj.stem_model.system_size
                                blocks=0:80:size(temp1,1);
                                if not(blocks(end)==size(temp1,1))
                                    blocks=[blocks size(temp1,1)];
                                end
                                for i=1:length(blocks)-1
                                    temp2(i)=temp1(blocks(i)+1:blocks(i+1),:)*data.X_time(blocks(i)+1:blocks(i+1),:,tT)';
                                end
                            else
                                temp2=temp1*data.X_time(:,:,tT)';
                            end
                            sum_num=sum_num-sum(temp2(Lt));
                        end
                        
                        if par.k>0
                            if data.X_g_tv
                                tG=t;
                            else
                                tG=1;
                            end
                            for k=1:K
                                temp1=D_apply(D_apply(M_cov_wr_wg_y1(:,t,k),data.X_rg(:,1,tRG),'l'),j_r,'l');
                                temp2=[data.X_g(:,1,tG,k);zeros(size(temp1,1)-size(data.X_g(:,1,tG,k),1),1)];
                                temp1=D_apply(D_apply(temp1',temp2,'r'),aj_g(:,k),'r');
                                sum_num=sum_num-sum(temp1(Lt));
                            end
                        end
                        
                        temp1=E_wr_y1(:,t).^2+diag_Var_wr_y1(:,t);
                        temp1=D_apply(D_apply(M_apply(temp1,M,'l'),data.X_rg(:,1,tRG),'b'),j_r,'b');
                        sum_den=sum_den+sum(temp1(Lt));
                    end
                    alpha_rg(r)=sum_num/sum_den;
                end
                st_par_em_step.alpha_rg=alpha_rg';
                ct2=clock;
                disp(['    alpha_rg update ended in ',stem_time(etime(ct2,ct1))]);

                disp('    v_r update started...');
                ct1=clock;
                if Nr<=obj.stem_EM_options.mstep_system_size
                    temp=zeros(size(sum_Var_wr_y1));
                    for t=1:T
                        temp=temp+E_wr_y1(:,t)*E_wr_y1(:,t)';
                    end
                    temp=temp+sum_Var_wr_y1;
                end
                r = symamd(sum_Var_wr_y1); %note that sum_Var_wr_y1 and sigma_W_r have the same sparseness structure
                if par.remote_correlated
                    %indices are permutated in order to avoid deadlock
                    kindex=randperm(size(par.v_r,1));
                    for k=kindex
                        hindex=randperm(size(par.v_r,1)-k)+k;
                        for h=hindex
                            initial=par.v_r(k,h);
                            if Nr<=obj.stem_EM_options.mstep_system_size
                                min_result = fminsearch(@(x) stem_geo_coreg_function_velement(x,k,h,par.v_r,par.theta_r,par.correlation_type,data.DistMat_r,...
                                    data.stem_varset_r.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_r.tap,r),initial,optimset('MaxIter',50,'TolX',1e-3,'UseParallel','always'));
                            else
                                disp('WARNING: this operation will take a long time');
                                min_result = fminsearch(@(x) stem_geo_coreg_function_velement(x,k,h,par.v_r,par.theta_r,par.correlation_type,data.DistMat_r,...
                                    data.stem_varset_r.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_r.tap,r),initial,optimset('MaxIter',50,'TolX',1e-3,'UseParallel','always'));
                            end
                            st_par_em_step.v_r(k,h)=min_result;
                            st_par_em_step.v_r(h,k)=min_result;
                        end
                    end
                    ct2=clock;
                    disp(['    v_r update ended in ',stem_time(etime(ct2,ct1))]);
                end

                disp('    theta_r updating started...');
                ct1=clock;
                initial=par.theta_r;
                
                if Nr<=obj.stem_EM_options.mstep_system_size
                    min_result = fminsearch(@(x) stem_geo_coreg_function_theta(x,par.v_r,par.correlation_type,data.DistMat_r,...
                        data.stem_varset_r.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_r.tap,r),log(initial),optimset('MaxIter',50,'TolX',1e-3,'UseParallel','always'));
                    st_par_em_step.theta_r=exp(min_result);
                else
                    if data.stem_varset_r.nvar>1
                        disp('WARNING: this operation will take a long time');
                        min_result = fminsearch(@(x) stem_geo_coreg_function_theta(x,par.v_r,par.correlation_type,data.DistMat_r,...
                            data.stem_varset_r.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_r.tap,r),log(initial),optimset('MaxIter',50,'TolX',1e-3,'UseParallel','always'));
                        st_par_em_step.theta_r=exp(min_result);
                    else
                        s=ceil(Nr/obj.stem_EM_options.mstep_system_size);
                        step=ceil(Nr/s);
                        blocks=0:step:Nr;
                        if not(blocks(end)==Nr)
                            blocks=[blocks Nr];
                        end
                        for j=1:length(blocks)-1
                            block_size=blocks(j+1)-blocks(j);
                            idx=blocks(j)+1:blocks(j+1);
                            temp=zeros(block_size);
                            for t=1:T
                                temp=temp+E_wr_y1(idx,t)*E_wr_y1(idx,t)';
                            end
                            temp=temp+sum_Var_wr_y1(idx,idx);
                            r_partial=symamd(sum_Var_wr_y1(idx,idx));
                            min_result(j,:) = fminsearch(@(x) stem_geo_coreg_function_theta(x,par.v_r,par.correlation_type,data.DistMat_r(idx,idx),...
                                length(idx),temp,t,obj.stem_model.stem_data.stem_gridlist_r.tap,r_partial),log(initial),optimset('maxiter',50,'tolx',1e-3,'useparallel','always'));
                        end
                        st_par_em_step.theta_r=exp(mean(min_result));
                    end
                end
                ct2=clock;
                disp(['    theta_r update ended in ',stem_time(etime(ct2,ct1))]);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %          alpha_g               %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if not(isempty(data.X_g))
                disp('    alpha_g update started...');
                for s=1:K
                    for r=1:obj.stem_model.stem_data.stem_varset_g.nvar
                        [aj_g_rs,j_r] = obj.stem_model.get_jg(r,s);
                        sum_num=0;
                        sum_den=0;
                        for t=1:T
                            if data.X_rg_tv
                                tRG=t;
                            else
                                tRG=1;
                            end
                            if data.X_time_tv
                                tT=t;
                            else
                                tT=1;
                            end
                            if data.X_g_tv
                                tG=t;
                            else
                                tG=1;
                            end
                            Lt=not(isnan(data.Y(:,t)));
                            
                            temp1=E_e_y1(:,t)+D_apply(D_apply(E_wg_y1(:,t,s),data.X_g(:,1,tG,s),'l'),aj_g_rs,'l');
                            temp2=D_apply(D_apply(E_wg_y1(:,t,s)',data.X_g(:,1,tG,s),'r'),j_r,'r');
                            sum_num=sum_num+sum(temp1(Lt).*temp2(Lt)');
                            
                            if par.p>0
                                temp1=D_apply(D_apply(cov_wg_z_y1(:,:,t,s),data.X_g(:,1,tG,s),'l'),j_r,'l');
                                temp2=zeros(size(temp1,1),1); 
                                if N>obj.stem_model.system_size
                                    blocks=0:80:size(temp1,1);
                                    if not(blocks(end)==size(temp1,1))
                                        blocks=[blocks size(temp1,1)];
                                    end
                                    for i=1:length(blocks)-1
                                        temp2(blocks(i)+1:blocks(i+1))=temp1(blocks(i)+1:blocks(i+1),:)*data.X_time(blocks(i)+1:blocks(i+1),:,tT)';
                                    end
                                else
                                    temp2=temp1*data.X_time(:,:,tT)';    
                                end
                                sum_num=sum_num-sum(temp2(Lt));
                            end
                            
                            if K>1
                                for k=1:K
                                    if not(k==s)
                                        if k<s
                                            kk=s;
                                            ss=k;
                                        else
                                            kk=k;
                                            ss=s;
                                        end
                                        temp1=D_apply(D_apply(cov_wgk_wgh_y1{kk,ss}(:,t),data.X_g(:,1,tG,k),'l'),aj_g(:,k),'l');
                                        temp1=D_apply(D_apply(temp1',[data.X_g(:,1,tG,s);zeros(Nr,1)],'r'),j_r,'r');
                                        sum_num=sum_num-sum(temp1(Lt));
                                    end
                                end
                            end

                            if not(isempty(data.X_rg))
                                temp1=D_apply(D_apply(M_cov_wr_wg_y1(:,t,s),data.X_rg(:,1,tRG),'l'),aj_rg,'l');
                                temp2=[data.X_g(:,1,tG,s);zeros(size(temp1,1)-size(data.X_g(:,1,tG,s),1),1)];
                                temp1=D_apply(D_apply(temp1',temp2,'r'),j_r,'r');
                                sum_num=sum_num-sum(temp1(Lt));
                            end
                            
                            temp1=E_wg_y1(:,t,s).^2+diag_Var_wg_y1(:,t,s);
                            temp1=D_apply(D_apply(temp1,data.X_g(:,1,tG,s),'b'),j_r,'b');
                            sum_den=sum_den+sum(temp1(Lt));
                        end
                        alpha_g(r,s)=sum_num/sum_den;
                    end
                end
                st_par_em_step.alpha_g=alpha_g;
                ct2=clock;
                disp(['    alpha_g update ended in ',stem_time(etime(ct2,ct1))]);
                %indices are permutated in order to avoid deadlock
                disp('    v_g and theta_g update started...');
                ct1=clock;
                
                %LA DIFFERENZA TRA M-STEP E M-STEP_PARALLEL E' QUESTO CICLO FOR!!! CHE NON VA 1:K MA SU INDEX
                r = symamd(sum_Var_wg_y1{1}); %note that sum_Var_wr_y1 and sigma_W_g{k} have the same sparseness structure
                for z=index
                    if Ng<=obj.stem_EM_options.mstep_system_size
                        temp=zeros(size(sum_Var_wg_y1{z}));
                        for t=1:T
                            temp=temp+E_wg_y1(:,t,z)*E_wg_y1(:,t,z)';
                        end
                        temp=temp+sum_Var_wg_y1{z};
                    end
                    
                    kindex=randperm(size(par.v_g(:,:,z),1));
                    for k=kindex
                        hindex=randperm(size(par.v_g(:,:,z),1)-k)+k;
                        for h=hindex
                            initial=par.v_g(k,h,z);
                            if Ng<=obj.stem_EM_options.mstep_system_size
                                min_result = fminsearch(@(x) stem_geo_coreg_function_velement(x,k,h,par.v_g(:,:,z),par.theta_g(z),par.correlation_type,data.DistMat_g,...
                                    data.stem_varset_g.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_g.tap,r),initial,optimset('MaxIter',50,'TolX',1e-3));
                            else
                                disp('WARNING: this operation will take a long time');
                                min_result = fminsearch(@(x) stem_geo_coreg_function_velement(x,k,h,par.v_g(:,:,z),par.theta_g(z),par.correlation_type,data.DistMat_g,...
                                    data.stem_varset_g.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_g.tap,r),initial,optimset('MaxIter',50,'TolX',1e-3));
                            end                            
                            st_par_em_step.v_g(k,h,z)=min_result;
                            st_par_em_step.v_g(h,k,z)=min_result;
                        end
                    end
                    
                    initial=par.theta_g(z);
                    if Ng<=obj.stem_EM_options.mstep_system_size
                        min_result = fminsearch(@(x) stem_geo_coreg_function_theta(x,par.v_g(:,:,z),par.correlation_type,data.DistMat_g,...
                            data.stem_varset_g.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_g.tap,r),log(initial),optimset('MaxIter',50,'TolX',1e-3));
                        st_par_em_step.theta_g(z)=exp(min_result);
                    else
                        if data.stem_varset_g.nvar>1
                            disp('WARNING: this operation will take a long time');
                            min_result = fminsearch(@(x) stem_geo_coreg_function_theta(x,par.v_g(:,:,z),par.correlation_type,data.DistMat_g,...
                                data.stem_varset_g.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_g.tap,r),log(initial),optimset('MaxIter',50,'TolX',1e-3));
                            st_par_em_step.theta_g(z)=exp(min_result);
                        else
                            s=ceil(Ng/obj.stem_EM_options.mstep_system_size);
                            step=ceil(Ng/s);
                            blocks=0:step:Ng;
                            if not(blocks(end)==Ng)
                                blocks=[blocks Ng];
                            end
                            for j=1:length(blocks)-1
                                block_size=blocks(j+1)-blocks(j);
                                idx=blocks(j)+1:blocks(j+1);
                                temp=zeros(block_size);
                                for t=1:T
                                    temp=temp+E_wg_y1(idx,t,z)*E_wg_y1(idx,t,z)';
                                end
                                temp=temp+sum_Var_wg_y1{z}(idx,idx);
                                r_partial=symamd(sum_Var_wg_y1{z}(idx,idx));
                                min_result(j,:) = fminsearch(@(x) stem_geo_coreg_function_theta(x,par.v_g(:,:,z),par.correlation_type,data.DistMat_g(idx,idx),...
                                    length(idx),temp,t,obj.stem_model.stem_data.stem_gridlist_g.tap,r_partial),log(initial),optimset('maxiter',50,'tolx',1e-3,'useparallel','always'));
                            end
                            st_par_em_step.theta_g(z)=exp(mean(min_result));
                        end
                    end
                end
                ct2=clock;
                disp(['    v_g and theta_g update ended in ',stem_time(etime(ct2,ct1))]);
            end
            obj.stem_model.stem_par=st_par_em_step;
            ct2_mstep=clock;
            disp(['  M step ended in ',stem_time(etime(ct2_mstep,ct1_mstep))]);
        end 
        
        function st_par_em_step = M_step_vg_and_theta(obj,E_wg_y1,sum_Var_wg_y1,index,r)
            st_par_em_step=obj.stem_model.stem_par;
            Ng=obj.stem_model.stem_data.stem_varset_g.N;
            for z=index
                
                if Ng<=obj.stem_EM_options.mstep_system_size
                    temp=zeros(size(sum_Var_wg_y1{z-index(1)+1}));
                    for t=1:size(E_wg_y1,2)
                        temp=temp+E_wg_y1(:,t,z-index(1)+1)*E_wg_y1(:,t,z-index(1)+1)';
                    end
                    temp=temp+sum_Var_wg_y1{z-index(1)+1};
                end
                
                kindex=randperm(size(st_par_em_step.v_g(:,:,z),1));
                for k=kindex
                    hindex=randperm(size(st_par_em_step.v_g(:,:,z),1)-k)+k;
                    for h=hindex
                        initial=st_par_em_step.v_g(k,h,z);
                        if Ng<=obj.stem_EM_options.mstep_system_size
                            min_result = fminsearch(@(x) stem_geo_coreg_function_velement(x,k,h,st_par_em_step.v_g(:,:,z),st_par_em_step.theta_g(z),st_par_em_step.correlation_type,data.DistMat_g,...
                                data.stem_varset_g.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_g.tap,r),initial,optimset('MaxIter',50,'TolX',1e-3));
                        else
                            disp('WARNING: this operation will take a long time');
                            min_result = fminsearch(@(x) stem_geo_coreg_function_velement(x,k,h,st_par_em_step.v_g(:,:,z),st_par_em_step.theta_g(z),st_par_em_step.correlation_type,data.DistMat_g,...
                                data.stem_varset_g.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_g.tap,r),initial,optimset('MaxIter',50,'TolX',1e-3));
                        end
                        st_par_em_step.v_g(k,h,z)=min_result;
                        st_par_em_step.v_g(h,k,z)=min_result;
                    end
                end
                
                initial=st_par_em_step.theta_g(z);
                if Ng<=obj.stem_EM_options.mstep_system_size
                    min_result = fminsearch(@(x) stem_geo_coreg_function_theta(x,st_par_em_step.v_g(:,:,z),st_par_em_step.correlation_type,obj.stem_model.stem_data.DistMat_g,...
                        obj.stem_model.stem_data.stem_varset_g.dim,temp,obj.stem_model.stem_data.T,obj.stem_model.stem_data.stem_gridlist_g.tap,r),log(initial),optimset('MaxIter',25,'TolX',1e-3));
                    st_par_em_step.theta_g(z)=exp(min_result);
                else
                    if obj.stem_model.stem_data.stem_varset_g.nvar>1
                        disp('WARNING: this operation will take a long time');
                        min_result = fminsearch(@(x) stem_geo_coreg_function_theta(x,st_par_em_step.v_g(:,:,z),st_par_em_step.correlation_type,obj.stem_model.stem_data.DistMat_g,...
                            obj.stem_model.stem_data.stem_varset_g.dim,temp,obj.stem_model.stem_data.T,obj.stem_model.stem_data.stem_gridlist_g.tap,r),log(initial),optimset('MaxIter',25,'TolX',1e-3));
                        st_par_em_step.theta_g(z)=exp(min_result);
                    else
                        s=ceil(Ng/obj.stem_EM_options.mstep_system_size);
                        step=ceil(Ng/s);
                        blocks=0:step:Ng;
                        if not(blocks(end)==Ng)
                            blocks=[blocks Ng];
                        end
                        for j=1:length(blocks)-1
                            block_size=blocks(j+1)-blocks(j);
                            idx=blocks(j)+1:blocks(j+1);
                            temp=zeros(block_size);
                            for t=1:size(E_wg_y1,2)
                                temp=temp+E_wg_y1(idx,t,z-index(1)+1)*E_wg_y1(idx,t,z-index(1)+1)';
                            end
                            temp=temp+sum_Var_wg_y1{z-index(1)+1}(idx,idx);
                            r_partial=symamd(sum_Var_wg_y1{z-index(1)+1}(idx,idx));
                            min_result(j,:) = fminsearch(@(x) stem_geo_coreg_function_theta(x,st_par_em_step.v_g(:,:,z),st_par_em_step.correlation_type,obj.stem_model.stem_data.DistMat_g(idx,idx),...
                                length(idx),temp,obj.stem_model.stem_data.T,obj.stem_model.stem_data.stem_gridlist_g.tap,r_partial),log(initial),optimset('MaxIter',25,'TolX',1e-3));
                        end
                        st_par_em_step.theta_g(z)=exp(mean(min_result));
                    end
                end
            end
        end
        
        function set.stem_model(obj,stem_model)
            if strcmp(class(stem_model),'stem_model')
                obj.stem_model=stem_model;
            else
                error('You have to provide an object of class stem_model');
            end            
        end
        
    end
end

