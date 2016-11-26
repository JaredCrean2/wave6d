% generate surface plots of nblocked loops vs blocksize
clear
clc
close all

maxdim = 6;
blocksizes = [2 4 8 16];

for d=1:maxdim
    fname_d = strcat('d', num2str(d), 'data.dat')
    data_d = dlmread(fname_d);
    
    figure
    nblocked_loops = 0:d
    
    surf(nblocked_loops, blocksizes, data_d.')
    shading faceted
    xlabel('number of blocked loops')
    ylabel('block size')
    zlabel('execution time (s)')
    graphtitle = strcat(num2str(d), 'd loop blocking')
    title(graphtitle)
    xlim([0 d])
    ylims = [min(blocksizes) max(blocksizes)]
    ylim(ylims)
    
    sname = strcat('d', num2str(d), 'data.png')
    print(sname, '-dpng', '-r600')
    
    view(-180, 0)
    sname = strcat('d', num2str(d), 'data16.png')
    print(sname, '-dpng', '-r600')
end