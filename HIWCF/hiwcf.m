% scalable neighborhood-based model to predict unknown observations
% 
clearvars;close all;clc;

datapath = './GDSC_data/';
cellline_drug_response = load([datapath, 'cell-line_drug_response.txt']);
Cellsim_probe = load([datapath, 'cellline_similarity.txt']);
Drugsim_fig_mt = load([datapath, 'drug_similarity.txt']);

cell_num = size(Cellsim_probe,1);
drug_num = size(Drugsim_fig_mt,1);

% data extraction, cellline id, drug id and the response score of each
% cellline-drug pair
cellline_id = cellline_drug_response(:,1);
drug_id = cellline_drug_response(:,2);
num = cellline_drug_response(:,3);

% baseline estimates
% the overall average rating
mu = mean(num);

% For each item i
lambda2 = 5;
for i = 1:drug_num
    id = find(cellline_drug_response(:,2) == i);
    b_i(i) = sum(num(id) - mu)/(lambda2 + length(id));
end

% For each user u
lambda3 = 2;
for u = 1:cell_num
    id = find(cellline_drug_response(:,1) == u);
    b_u(u) = sum(num(id) - b_i(cellline_drug_response(id,2))' - mu)/(lambda3 + length(id));
end

cellline_drug_IC50_0 = NaN(cell_num, drug_num);
for id = 1:length(cellline_drug_response)
    cellline_drug_IC50(cellline_drug_response(id,1),cellline_drug_response(id,2)) = num(id);
end

N = 15;
lambda4 = 50;
drugwisecorr = NaN(size(Drugsim_fig_mt,1),1);
drugwise_qt = NaN(size(Drugsim_fig_mt,1),1);
drugwiseerr = NaN(size(Drugsim_fig_mt,1),1);
drugwiseerr_qt = NaN(size(Drugsim_fig_mt,1),1);
drugwiserepn = NaN(size(Drugsim_fig_mt,1),1);

% Unknown ones
% cellline_id = cellline_drug_response(:,1);
% drug_id = cellline_drug_response(:,2);
% num = cellline_drug_response(:,3);
response_unknown = NaN(cell_num, drug_num);
for id = 1:length(num)
    response_unknown(cellline_drug_response(id,1),cellline_drug_response(id,2)) = cellline_drug_response(id,3);
end
[row, col] = find(isnan(response_unknown));
unknown_cell_drug_id = [row col];
cellline_response_pred = NaN(length(unknown_cell_drug_id),3);

Response_train = num;
cell_train_id = cellline_id;
drug_train_id = drug_id;

for id = 1:length(unknown_cell_drug_id)
   % ÕÒ³ö´ýÔ¤²âµÄuser idºÍdrug id
   cell_test_id = unknown_cell_drug_id(id,1);
   drug_test_id = unknown_cell_drug_id(id,2);

   tmp_id = find(~isnan(cellline_drug_IC50(cell_test_id,:)));
   r_u = mean(cellline_drug_IC50(cell_test_id,tmp_id));
   tmp_id = find(~isnan(cellline_drug_IC50(:,drug_test_id)));
   [Y,I] = sort(Cellsim_probe(tmp_id,cell_test_id),1,'descend');
   N_u_i = tmp_id (I(2:1+N)); % nearest celllines that have the same response in drug d_i 

   sim = NaN(length(N_u_i),1);
   ri = NaN(length(N_u_i),1);
   sim1 = NaN(length(N_u_i),1);
   for v = 1:length(N_u_i)
        % similarity between cellline u and v from response matrix
        tmp_id = find(~isnan(cellline_drug_IC50(cell_test_id,:)) & ~isnan(cellline_drug_IC50(N_u_i(v),:)));
        if(length(tmp_id) > 1) 
           sim1(v) = corr(cellline_drug_IC50(cell_test_id,tmp_id)',cellline_drug_IC50(N_u_i(v),tmp_id)');
        else
           sim1(v) = 1;
        end
        % disp(length(tmp_id))
        sim1(v) = length(tmp_id)/(length(tmp_id) + lambda4) * sim1(v);  % 

        % similarity between cellline u and v from similarity matrix
        sim(v) = Cellsim_probe(cell_test_id,N_u_i(v)) * sim1(v);%;
        ri(v) = cellline_drug_IC50(N_u_i(v),drug_test_id) - mean(cellline_drug_IC50(N_u_i(v),~isnan(cellline_drug_IC50(N_u_i(v),:))));
    end

    sim = sim .^ 4;
    numpred = r_u + ri' * (sim) /sum(sim);
    if (isnan(numpred))
         disp(r_u);
         disp([ri sim sim1]);
         %return;
         numpred = r_u;
    end

    tmp_id = find(~isnan(cellline_drug_IC50(:,drug_test_id)));
    r_u = mean(cellline_drug_IC50(tmp_id,drug_test_id));
    tmp_id = find(~isnan(cellline_drug_IC50(cell_test_id,:)));
    [Y,I] = sort(Drugsim_fig_mt(tmp_id,drug_test_id),1,'descend');
    if (length(tmp_id) > 1)
        if (length(tmp_id) > N)
            N_u_i = tmp_id (I(2:1+N)); % nearest drugs that have the same response in drug d_i 
        else
            N_u_i = tmp_id (I(2:end)); % nearest drugs that have the same response in drug d_i 
        end

         sim = NaN(length(N_u_i),1);
         ri = NaN(length(N_u_i),1);
         sim1 = NaN(length(N_u_i),1);
         for v = 1:length(N_u_i)
             % similarity between drug d_i and v from response matrix
             tmp_id = find(~isnan(cellline_drug_IC50(:,drug_test_id)) & ~isnan(cellline_drug_IC50(:,N_u_i(v))));
             if(length(tmp_id) > 1) 
                 sim1(v) = corr(cellline_drug_IC50(tmp_id,drug_test_id),cellline_drug_IC50(tmp_id,N_u_i(v)));
             else
                 sim1(v) = 1;
             end
             sim1(v) = length(tmp_id)/(length(tmp_id) + lambda4) * sim1(v);

             % similarity between drug d_i and v from similarity matrix
             sim(v) = Drugsim_fig_mt(N_u_i(v),drug_test_id) * sim1(v);%
             ri(v) = cellline_drug_IC50(cell_test_id,N_u_i(v)) - mean(cellline_drug_IC50(~isnan(cellline_drug_IC50(:,N_u_i(v))),N_u_i(v)));
         end
         
         sim = sim .^ 4;
         numpred = r_u + ri' * (sim) /sum(sim) + numpred;
         if (isnan(numpred))
            disp(r_u);
            disp([ri sim sim1]);
            %return;
            numpred = numpred + r_u;
         end
         numpred = numpred/2;
    else
        numpred = (numpred + r_u)/2;
    end
    cellline_response_pred(id, :) = [cell_test_id drug_test_id numpred];
end

