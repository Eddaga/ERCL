test_target = readmatrix('dataframes.xlsx','Sheet','test_target');
test_output = readmatrix('dataframes.xlsx','Sheet','test_outputs');

test_data = [test_target(:,2) test_output(:,2)];

for i = 1 : 1 : length(test_data)
    test_data(i,3) = (test_data(i,1) - test_data(i,2)) / test_data(i,1);
end

xAxis = 1: 1: length(test_data);

figure(1)
plot(xAxis,test_data(:,3));

mvdTarget = movmean(test_data(:,1),100);
mvdOutput = movmean(test_data(:,2),100);

figure(2)
mvAveragedGraph = nexttile;
plot(mvAveragedGraph,xAxis,mvdTarget,xAxis,mvdOutput);
title(mvAveragedGraph,'LearnResult');
xlabel(mvAveragedGraph,'tick');
ylabel(mvAveragedGraph,'P(V*I)');
legend(mvAveragedGraph,'target','output');
