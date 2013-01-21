%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Finazzi                                    %
% e-mail: francesco.finazzi@unibg.it                           %
% Affiliation: University of Bergamo                           %
% Department: Information Technology and Mathematical Methods  %
%                                                              %
% Version: beta                                                %
% Release date: 15/05/2012                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
clear all

pathparallel='/opt/matNfs/';

%address = java.net.InetAddress.getLocalHost;
%split=regexp(char(address.getHostAddress),'\.','split');
h=now;
IPaddress = round((h-floor(h))*1000000);
disp([datestr(now),' - Machine code: ',num2str(IPaddress)]);
timeout=72000;
while(1)
    exit=0;
    disp([datestr(now),' - Waiting whois from server...']);
    while not(exit)
        exit=exist([pathparallel,'whoishere.mat'],'file');
        pause(0.1);
    end
    read=0;
    while not(read)
        try
            load([pathparallel,'whoishere.mat']);
            read=1;
            disp([datestr(now),' - whois received.']);
        catch
        end
        pause(0.1);
    end
    machine.IPaddress=IPaddress;
    machine.IDrequest=whoishere.IDrequest;
    if exist('st_model','var')
        machine.require_stemmodel=0;
    else
        machine.require_stemmodel=1;
    end

    save([pathparallel,'temp/machine_',num2str(IPaddress),'.mat'],'machine');
    pause(0.5);
    movefile([pathparallel,'temp/machine_',num2str(IPaddress),'.mat'],[pathparallel,'machine_',num2str(IPaddress),'.mat']);
    disp([datestr(now),' - Replied to whois.']);
    
    if machine.require_stemmodel
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %           st_model         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp([datestr(now),' - Waiting for st_model...']);
        exit=0;
        ct1=clock;
        while not(exit)
            exit=exist([pathparallel,'st_model_parallel_',num2str(IPaddress),'.mat'],'file');
            pause(0.1);
            ct2=clock;
            if etime(ct2,ct1)>timeout 
                disp('Timeout in waiting for st_model_parallel');
                exit=1;
            end
        end
        read=0;
        ct1=clock;
        while not(read)
            try
                load([pathparallel,'st_model_parallel_',num2str(IPaddress),'.mat']);
                disp([datestr(now),' - st_model received.']);
                read=1;
            catch
            end
            pause(0.1);
            ct2=clock;
            if etime(ct2,ct1)>timeout
                disp('Timeout in reading st_model_parallel');
                read=1;
            end
        end
        deleted=0;
        ct1=clock;
        while not(deleted)
            try
                delete([pathparallel,'st_model_parallel_',num2str(IPaddress),'.mat']);
                disp([datestr(now),' - st_model deleted.']);
                deleted=1;
            catch
            end
            pause(0.1);
            ct2=clock;
            if etime(ct2,ct1)>timeout
                disp('Timeout in deleting st_model_parallel');
                deleted=1;
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           st_par           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp([datestr(now),' - Waiting for st_par']);
    exit=0;
    ct1=clock;
    while not(exit)
        exit=exist([pathparallel,'st_par_parallel_',num2str(IPaddress),'.mat'],'file');
        pause(0.1);
        ct2=clock;
        if etime(ct2,ct1)>timeout 
            disp('Timeout in waiting for st_par_parallel');
            exit=1;
        end
    end
    read=0;
    ct1=clock;
    while not(read)
        try
            load([pathparallel,'st_par_parallel_',num2str(IPaddress),'.mat']);
            disp([datestr(now),' - st_par received.']);
            read=1;
        catch
        end
        pause(0.05);
        ct2=clock;
        if etime(ct2,ct1)>timeout 
            disp('Timeout in reading st_par_parallel');
            read=1;
        end
    end
    deleted=0;
    ct1=clock;
    while not(deleted)
        try
            delete([pathparallel,'st_par_parallel_',num2str(IPaddress),'.mat']);
            disp([datestr(now),' - st_par deleted.']);
            deleted=1;
        catch
        end
        pause(0.1);
        ct2=clock;
        if etime(ct2,ct1)>timeout 
            disp('Timeout in deleting st_par_parallel');
            deleted=1;
        end
    end
    st_model.stem_par=st_par;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           Kalman           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if st_model.stem_par.p>0
        disp([datestr(now),' - Waiting for kalman_parallel...']);
        exit=0;
        ct1=clock;
        while not(exit)
            exit=exist([pathparallel,'kalman_parallel_',num2str(IPaddress),'.mat'],'file');
            pause(0.1);
            ct2=clock;
            if etime(ct2,ct1)>timeout 
                disp('Timeout in waiting for kalman_parallel');
                exit=1;
            end
        end
        read=0;
        ct1=clock;
        while not(read)
            try
                load([pathparallel,'kalman_parallel_',num2str(IPaddress),'.mat']);
                disp([datestr(now),' - kalman_parallel received']);
                read=1;
            catch
            end
            pause(0.1);
            ct2=clock;
            if etime(ct2,ct1)>timeout 
                disp('Timeout in loading kalman_parallel');
                read=1;
            end
        end
        deleted=0;
        ct1=clock;
        while not(deleted)
            try
                delete([pathparallel,'kalman_parallel_',num2str(IPaddress),'.mat']);
                disp([datestr(now),' - kalman_parallel deleted']);
                deleted=1;
            catch
            end
            pause(0.1);
            ct2=clock;
            if etime(ct2,ct1)>timeout 
                disp('Timeout in deleting kalman_parallel');
                deleted=1;
            end
        end
        
        st_kalman=stem_kalman(st_model);
        st_kalman.filter(0,data.time_steps,pathparallel);
        clear data
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         EM -  E step       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp([datestr(now),' - Waiting for data_parallel...']);
    exit=0;
    ct1=clock;
    while not(exit)
        exit=exist([pathparallel,'data_parallel_',num2str(IPaddress),'.mat'],'file');
        pause(0.1);
        ct2=clock;
        if etime(ct2,ct1)>timeout 
            disp('Timeout in waiting for data_parallel');
            exit=1;
        end
    end
    read=0;
    ct1=clock;
    while not(read)
        try
            load([pathparallel,'data_parallel_',num2str(IPaddress),'.mat']);
            disp([datestr(now),' - data_parallel received']);
            read=1;
        catch
        end
        pause(0.1);
        ct2=clock;
        if etime(ct2,ct1)>timeout 
            disp('Timeout in loading data_parallel');
            read=1;
        end
    end
    deleted=0;
    ct1=clock;
    while not(deleted)
        try
            delete([pathparallel,'data_parallel_',num2str(IPaddress),'.mat']);
            disp([datestr(now),' - data_parallel deleted']);
            deleted=1;
        catch
        end
        pause(0.1);
        ct2=clock;
        if etime(ct2,ct1)>timeout 
            disp('Timeout in deleting data_parallel');
            deleted=1;
        end
    end
    
    
    st_EM_options=stem_EM_options(0.001,1,'single',[],0,[]);
    st_EM=stem_EM(st_model,st_EM_options);
    ct1=clock;
    [output.E_wr_y1,output.sum_Var_wr_y1,output.diag_Var_wr_y1,output.cov_wr_z_y1,output.E_wg_y1,...
        output.sum_Var_wg_y1,output.diag_Var_wg_y1,output.cov_wg_z_y1,output.M_cov_wr_wg_y1,...
        output.cov_wgk_wgh_y1,output.diag_Var_e_y1,output.E_e_y1] = st_EM.E_step_parallel(data.time_steps,data.st_kalmansmoother_result);
    output.sum_Var_wr_y1=triu(output.sum_Var_wr_y1);
    for k=1:length(output.sum_Var_wg_y1)
        output.sum_Var_wg_y1{k}=triu(output.sum_Var_wg_y1{k});
    end
    ct2=clock;
    output.ct=etime(ct2,ct1);
    output.cb=data.cb;
    output.iteration=data.iteration;
    output.time_steps=data.time_steps;
    output.IPaddress=IPaddress;
    disp([datestr(now),' - saving output...']);
    save([pathparallel,'temp/output_',num2str(IPaddress),'.mat'],'output');
    pause(0.5);
    movefile([pathparallel,'temp/output_',num2str(IPaddress),'.mat'],[pathparallel,'output_',num2str(IPaddress),'.mat']);
    disp([datestr(now),' - output saved.']);
    clear data
    clear output
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         EM - M-Step        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp([datestr(now),' - Waiting for data_parallel_mstep...']);
    exit=0;
    ct1=clock;
    while not(exit)
        exit=exist([pathparallel,'data_parallel_mstep',num2str(IPaddress),'.mat'],'file');
        pause(0.1);
        ct2=clock;
        if etime(ct2,ct1)>timeout 
            disp('Timeout in waiting for data_parallelmstep');
            exit=1;
        end
    end
    read=0;
    ct1=clock;
    while not(read)
        try
            load([pathparallel,'data_parallel_mstep',num2str(IPaddress),'.mat']);
            disp([datestr(now),' - data_parallelmstep received']);
            read=1;
        catch
        end
        pause(0.1);
        ct2=clock;
        if etime(ct2,ct1)>timeout 
            disp('Timeout in loading data_parallelmstep');
            read=1;
        end
    end
    deleted=0;
    ct1=clock;
    while not(deleted)
        try
            delete([pathparallel,'data_parallel_mstep',num2str(IPaddress),'.mat']);
            disp([datestr(now),' - data_parallelmstep deleted']);
            deleted=1;
        catch
        end
        pause(0.1);
        ct2=clock;
        if etime(ct2,ct1)>timeout 
            disp('Timeout in deleting data_parallelmstep');
            deleted=1;
        end
    end    
    
    output.iteration=data.iteration;
    output.IPaddress=IPaddress;
    if not(isempty(data.index))
        disp([datestr(now),' - m-step running...']);
        ct1=clock;
        output.mstep_par=st_EM.M_step_vg_and_theta(data.E_wg_y1,data.sum_Var_wg_y1,data.index,data.r);
        ct2=clock;
        output.ct=etime(ct2,ct1);
        output.index=data.index;
        disp([datestr(now),' - saving output_mstep...']);
    else
        output.index=[];
        disp([datestr(now),' - saving output_mstep empty...']);
    end
    save([pathparallel,'temp/output_mstep_',num2str(IPaddress),'.mat'],'output');
    pause(0.5);
    movefile([pathparallel,'temp/output_mstep_',num2str(IPaddress),'.mat'],[pathparallel,'output_mstep_',num2str(IPaddress),'.mat']);
    disp([datestr(now),' - output_mstep saved.']);
    clear data
    clear output
end
