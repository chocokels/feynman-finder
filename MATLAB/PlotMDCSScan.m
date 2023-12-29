%% Sept 6, 2023 - Kelsey Bates
% Recursive function to find Feynman diagrams of a system.
%
% Input:
%   ts:         
%   dat:        
%   four:  0 for time domain, 1 for frequency domain
%   compx: 0 for absolute value, 1 for real, 2 for imaginary
%
% Output:
%   res:  Outputs figure number used, or 0 on failure

function [dat,ax1,ax2,fignum] = PlotMDCSScan(ts,dat,four,compx)

axs = find(size(dat) > 1);

if length(axs) ~= 2
    disp('dat should have two dimensions of length greater than 1');
    fignum = 0;
    return
end

ax1 = ts{axs(1)};
ax2 = ts{axs(2)};

dat = squeeze(dat);

if four
    dat = ifftshift(ifftshift(ifft(ifft(dat,length(ax1),1),length(ax2),2),1),2);
    
    ax1 = 2*pi/(ax1(2)-ax1(1)) * ((0:length(ax1)-1)/length(ax1)-1/2);
    ax2 = 2*pi/(ax2(2)-ax2(1)) * ((0:length(ax2)-1)/length(ax2)-1/2);
end

if compx == 0
    dat = abs(dat);
    lims = [0 max(dat,[],'all')];
elseif compx == 1
    dat = real(dat);
    lims = [-1 1] * max(abs(dat),[],'all');
elseif compx == 2
    dat = imag(dat);
    lims = [-1 1] * max(abs(dat),[],'all');
elseif compx == 3
    dat = angle(dat);
    lims = [-pi pi];
end

imagesc(ax2,ax1,dat,lims);
set(gca,'YDir','normal');

if length(ts) == 3
    if ~four
        labs = {'\tau','T','t'};
    else
        labs = {'\omega_\tau','\omega_T','\omega_t'};
    end
    xlabel(labs{axs(2)});
    ylabel(labs{axs(1)});
else
    if ~four
        xlabel(['t_' num2str(ax(2))]);
        ylabel(['t_' num2str(ax(1))]);
    else
        xlabel(['\omega_{t' num2str(axs(2)) '}']);
        ylabel(['\omega_{t' num2str(axs(1)) '}']);
    end
end

fignum = get(gcf,'Number');

end