function [ctl_corrs, hypo_corrs, ctlVShypo, mini_winCorr] = PBC_bursts(filename1, filename2)
load(filename1) % Control_5_rafagas.mat
load(filename2) % Post_5_rafagas.mat

ctl_corrs = zeros(5,5);
hypo_corrs = zeros(5,5);
ctlVShypo = zeros(5,5);
all_bursts = [control; post];

for burst = 1:(size(control, 1))
auto_corr = max(xcorr(control(burst,:)));
for next_burst = 1:(size(control, 1))
cross_burst = xcorr(control(burst,:), control(next_burst, :));
% plot(cross_burst)
% shg
% pause;
normal_cross = max(cross_burst)/max(auto_corr);
row = burst;
col  = next_burst;
ctl_corrs(row, col) = normal_cross;
ctl_corrs(col, row) = normal_cross;
end
end

for burst = 1:(size(post, 1))
auto_corr = max(xcorr(post(burst,:)));
for next_burst = 1:(size(post, 1))
cross_burst = xcorr(post(burst,:), post(next_burst, :));
% plot(cross_burst)
% shg
% pause;
normal_cross = max(cross_burst)/max(auto_corr);
row = burst;
col  = next_burst;
hypo_corrs(row, col) = normal_cross;
hypo_corrs(col, row) = normal_cross;
end
end

for burst = 1:(size(control, 1))
for next_burst = 1:(size(post, 1))
cross_burst = xcorr(control(burst,:), post(next_burst, :));
% plot(cross_burst)
% shg
% pause;
%normal_cross = max(cross_burst)/max(auto_corr);
row = burst;
col  = next_burst;
ctlVShypo(row, col) = max(cross_burst);
ctlVShypo(col, row) = max(cross_burst);
end
end
%% Descomponer rafagas en mini ventanas (25ms)
mini_bursts = zeros(10*20, 625);
row_count = 1;
for ii = 1:size(all_bursts, 1)
mini_start = 1;
mini_stop = 625;
for jj = 1:20
    mini_bursts(row_count, :) = all_bursts(ii, (mini_start:mini_stop));
    row_count = row_count + 1;
    mini_start = mini_start + 625;
    mini_stop = mini_stop + 625;
end
end
mini_winCorr = squareform(1-pdist(mini_bursts, 'correlation'));

%%
figure(1)
subplot(2,2,1)
imagesc(ctl_corrs)
colormap 'jet'
colorbar
set(gca, 'YTick', [1:5]);
title('Control');
axis 'square'

subplot(2,2,2)
imagesc(hypo_corrs)
colormap 'jet'
colorbar
set(gca, 'YTick', [1:5]);
title('Post');
axis 'square'

subplot(2,2, [3 4])
imagesc(ctlVShypo)
colormap 'jet'
colorbar
axis 'square'
set(gca, 'YTick', [1:5]);
title('Control vs Post');
xlabel('Control')
ylabel('Post')

figure(2)
imagesc(mini_winCorr)
colormap 'jet'
colorbar
title('Mini Bursts');
axis 'square'

return

