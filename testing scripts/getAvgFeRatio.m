% EL
% Feb 2023

% Get the average r_eq between Deng20 and Armstrong19 results in
% db\r_eq.xlsx


clear;

data_D20 = readmatrix('\db\r_eq.xlsx', 'Sheet', 'Deng20');
data_A19 = readmatrix('\db\r_eq.xlsx', 'Sheet', 'Armstrong19');

avg_req = zeros(size(data_D20)) + data_D20(:,1);

for i = 2:size(data_D20,2)         %columns, ignore P col
    for j = 1:size(data_D20,1)     %rows
        avg = (data_D20(j,i) + data_A19(j,i))/2;
        avg_req(j,i) = avg;
    end
end


writematrix(avg_req, 'D:\School Stuff\Yale\MOredox_2023\db\r_eq.xlsx', 'Sheet', 'Avg')