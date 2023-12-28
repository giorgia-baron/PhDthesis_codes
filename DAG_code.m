clear all
clc
close all

%% setting params
%only cortical ROIs, otherwise N_nodes=1:74 to include subcortical nodes
N_nodes=1:62;
n=max(N_nodes);

node_labels=node_labels(N_nodes);
N_subs=20;
% EC_all=EC_all(N_nodes,N_nodes,:);
group_num=group_num(N_nodes);
group=group(N_nodes);


files=dir('...'); %path to data
files=files(3:end);


%% find cycles and sort DAG

for ss=1:N_subs
    ss
    
    EC_sub=load(fullfile(files(ss).folder,files(ss).name)).A; 
    EC_sub=EC_sub(N_nodes,N_nodes);

    % estimate matrix S
    NoiseVar2=load(fullfile(files(ss).folder,files(ss).name)).NoiseVar;
    Q_sub=NoiseVar2.*eye(n);
    Sigma_sub=lyap(EC_sub,Q_sub); %zero-lag covariance matrix
    S_sub=(EC_sub*Sigma_sub-Sigma_sub*EC_sub'); %time-lag covariance matrix S
    EC_sub_check=0.5*(-Q_sub+S_sub)*inv(Sigma_sub); %test if EC_sub_check=EC_sub
    S_sub=S_sub'; %matlab convention: outgoing over rows and incoming over columns, need to transpose
    S_sub(S_sub<0)=0; % set negative weights to zero
    MATRIX_S(:,:,ss)=S_sub;
    all_edges_S=sum(length(nonzeros(S_sub)));
    S_sub=digraph(abs(S_sub),node_labels); %create directed graph
    

    for iter=1:100
        iter
        %set seed
        rng(iter,'Threefry')

        all_edges=all_edges_S;
        A_sub=S_sub;
        edges_removed=0;
        counter=1;
        
        while(hascycles(A_sub))
                
               [cycles,edgecycles] = allcycles(A_sub,'MaxNumCycles',100);
               idx_cycle_rand=randperm(length(edgecycles),1); %take random cycle and remove the lowest link
               w_idx=cell2mat(edgecycles(idx_cycle_rand,:));
               min_idx=1;
               w_min=A_sub.Edges.Weight(w_idx(min_idx));
               for ll=2:length(w_idx)
                      w_test=A_sub.Edges.Weight(w_idx(ll));
                       if abs(w_test)<abs(w_min)
                           w_min=w_test;
                           min_idx=ll;
                       end
               end
              
               cycles_rmd_subj(counter)=w_min;
               nodes_cycles_rmd_subj(counter,:)=A_sub.Edges(w_idx(min_idx),:).EndNodes;
               edges_removed=edges_removed+1;
               counter=counter+1;
               A_sub=rmedge(A_sub,w_idx(min_idx));
               clear edgecycles
        end
    
        cycles_rmd(iter,ss).weights=cycles_rmd_subj;
        cycles_rmd(iter,ss).nodes=nodes_cycles_rmd_subj;
        
        %save DAG matrix
        A_sub_mat=full(adjacency(A_sub,"weighted"));
        DAG(:,:,iter,ss)=A_sub_mat; %transposed wrt EC
       
        %save removed edges
        ER_perc(iter,ss)=edges_removed/all_edges*100;
        ER(iter,ss)=edges_removed;
    
        clear cycles_rmd_subj nodes_cycles_rmd_subj A_sub_mat w_idx w_min A_sub all_edges edges_removed counter

    end
  
    
    
    

end


save('...','DAG','MATRIX_EC_ANTISYM','MATRIX_S',...
    'ER','ER_perc','cycles_rmd')
