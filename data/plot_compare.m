% generate comparison of hilbert, non-blocking, blocking in a bar chart
clear
clc
close all

maxdim = 6;
blocksizes = [1 2 4 8 16];

for d=1:maxdim
    fname_d = strcat('bar', num2str(d), 'data.dat')
    data_d = dlmread(fname_d);
    
    % get the number of blocked loops, block size
    fname_d = strcat('bar', num2str(d), 'data2.dat')
    auxdata_d = dlmread(fname_d);
    
    figure
    lab2 = strcat(num2str(auxdata_d(1)), '-block-', num2str(auxdata_d(2)));
    lab3 = strcat('hilbert block-', num2str(auxdata_d(3)))
    xlabels = {'non-blocked', lab2, lab3}
    

    bar(data_d, 'b')
    set(gca, 'xticklabel', xlabels)
    
    fontsize = 12
    set(gca,'FontSize',fontsize)
%    xlabel(xlab, 'FontSize', fontsize)
    ylabel('execution time (s)')
    graphtitle = strcat(num2str(d), 'd Comparison')
    title(graphtitle)
    grid
    
    sname = strcat('comparison', num2str(d), 'data.png')
    print(sname, '-dpng', '-r600')
end
    
    