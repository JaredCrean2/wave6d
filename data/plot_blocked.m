% generate surface plots of nblocked loops vs blocksize
clear
clc
close all

maxdim = 6;
blocksizes_outer = [4 8 16 32];

for d=1:maxdim
    fname_d = strcat('d', num2str(d), 'data.dat')
    data_d = dlmread(fname_d);
    
    if d == 6
        blocksizes = blocksizes_outer(1:end-1)
    else
        blocksizes = blocksizes_outer
    end
    
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
    
%    view(-180, 0)
%    sname = strcat('d', num2str(d), 'data16.png')
%    print(sname, '-dpng', '-r600')

    % extract the maximum blocksize data and do a line plot
    figure
    linestyle='-';   % solid line (default)
    marker='x';     % marker type (blank by default)
    linecolor='b';  % line color blue (default)
    styleString=strcat(linestyle,marker,linecolor)
    linewidth=2;
    
    maxplot_data = data_d(:, end)
    plot(nblocked_loops, maxplot_data, styleString, 'LineWidth', linewidth)
    xlabel('number of blocked loops')
    ylabel('execution time (s)')
    graphtitle = strcat(num2str(d), 'd loop blocking, blocksize = ', num2str(blocksizes(end)))
    title(graphtitle)
    grid
    
    sname = strcat('d', num2str(d), 'data16.png')
    print(sname, '-dpng', '-r600')
    
end