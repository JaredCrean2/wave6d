% generate plots of hilbert curve block size vs execution time
clear
clc
close all

maxdim = 6;
blocksizes_outer = [1 4 8 16 32];

for d=2:maxdim
    fname_d = strcat('hd', num2str(d), 'data.dat')
    data_d = dlmread(fname_d);
    
    if d == 6
        blocksizes = blocksizes_outer(1:end-1)
    else
        blocksizes = blocksizes_outer
    end
    figure;
    linestyle='-';   % solid line (default)
    marker='x';     c% marker type (blank by default)
    linecolor='b';  % line color blue (default)
    styleString=strcat(linestyle,marker,linecolor)
    linewidth=2;
    
    plot(blocksizes, data_d, styleString, 'LineWidth', linewidth)
    xlabel('block size')
    ylabel('execution time (s)')
    graphtitle = strcat(num2str(d), 'd Hilbert blocking')
    title(graphtitle)
    grid
    
    sname = strcat('hd', num2str(d), 'data.png')
    print(sname, '-dpng', '-r600')
end