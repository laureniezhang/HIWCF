% investigate the similarity of gene expression profiles and response
% matrix
clearvars;close all;clc;

cellline_drug_response = load('./CCLE_data/cell-line_drug_response.txt');
Cellsim_probe = load('./CCLE_data/cellline_similarity.txt');
Drugsim_fig_mt = load('./CCLE_data/drug_similarity.txt');
[T,drug_list,RAW]  =xlsread('./CCLE_data/drug_id_name.xlsx');
[T,cellline_list,RAW]  =xlsread('./CCLE_data/cellline_id_name.xlsx');

cell_num = size(Cellsim_probe,1);
drug_num = size(Drugsim_fig_mt,1);

% 提取数据：cellline的id，drug的id和每一对cellline-drug的response score
cellline_id = cellline_drug_response(:,1);
drug_id = cellline_drug_response(:,2);

scale1 = cellline_drug_response(:,3);
num = cellline_drug_response(:,3);%/max(max(scale1),abs(min(scale1)));

%% removing global effect
% removing GE0, overall mean
GE0 = mean(num);
num = num - GE0;

% removing GE1
alpha1 = 3; % constanct when removing the GE effect
% calculate GE user main effect, explanatory variable x_u_i = 1
theta1 = NaN(cell_num,1); 
for i = 1:cell_num
    id_responsed = find(cellline_id == i);
    theta1(i) = sum(num(id_responsed))/(length(id_responsed)+alpha1);
    num(id_responsed) = num(id_responsed) - theta1(i);
end

% removing GE2
alpha2 = 10; % constanct when removing the GE effect
% calculate GE item main effect, explanatory variable x_u_i = 1
theta2 = NaN(drug_num,1); 
for j = 1:drug_num
    id_responsed = find(drug_id == j);
    theta2(j) = sum(num(id_responsed))/(length(id_responsed)+alpha2);
    num(id_responsed) = num(id_responsed) - theta2(j);
end
% num stores the residues after global effect removal

cellline_drug_IC50 = NaN(size(Cellsim_probe,1), size(Drugsim_fig_mt,1));
for i = 1:size(Cellsim_probe,1)
    for j = 1: size(Drugsim_fig_mt,1)
         id = find(cellline_drug_response(:,1) == i & cellline_drug_response(:,2) == j );
         if (~isnan(id))
             cellline_drug_IC50(i,j) = num(id);
         end
    end
end

% 从response matrix计算细胞系类似性
cell_sim = zeros(cell_num,cell_num);
for i = 1:cell_num
    for j = 1:cell_num
        if i <= j
            if i == j
                cell_sim(i,j) = 1;
            else
                nan_id = find(~isnan(cellline_drug_IC50(i,:)) & ~isnan(cellline_drug_IC50(j,:)));
                % disp(nan_id);
                if length(nan_id) > 1
                    cell_sim(i,j) = corr(cellline_drug_IC50(i,nan_id)',cellline_drug_IC50(j,nan_id)');
                else
                    cell_sim(i,j) = 0;
                end
                cell_sim(j,i) = cell_sim(i,j);
            end
        end
    end
end

triu_cell_sim = triu(cell_sim,1);
id = find(triu_cell_sim > 0.8);

% 从response matrix计算药物类似性
drug_sim = zeros(drug_num,drug_num);
for i = 1:drug_num
    for j = 1:drug_num
        if i <= j
            if i == j
                drug_sim(i,j) = 1;
            else
                nan_id = find(~isnan(cellline_drug_IC50(:,i)) & ~isnan(cellline_drug_IC50(:,j)));
                % disp(nan_id);
                if length(nan_id) > 1
                    drug_sim(i,j) = corr(cellline_drug_IC50(nan_id,i),cellline_drug_IC50(nan_id,j));
                else
                    drug_sim(i,j) = 0;
                end
                drug_sim(j,i) = drug_sim(i,j);
            end
        end
    end
end

clustergram(drug_sim, 'Standardize','none', 'RowLabels',drug_list, 'ColumnLabels',drug_list);
Drugsim_fig_mt = Drugsim_fig_mt - 0.6;
clustergram(Drugsim_fig_mt, 'Standardize','none', 'DisplayRange',0.4, 'RowLabels',drug_list, 'ColumnLabels',drug_list);

clustergram(cell_sim, 'Standardize',3, 'RowLabels',cellline_list, 'ColumnLabels',cellline_list);
Cellsim_probe = Cellsim_probe - 0.9;
clustergram(Cellsim_probe, 'Standardize',3, 'DisplayRange',0.2,'RowLabels',cellline_list, 'ColumnLabels',cellline_list);