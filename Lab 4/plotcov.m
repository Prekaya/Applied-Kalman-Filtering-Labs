tcl = tiledlayout(2,1);

nexttile;
hold on
plot(tHistory,sqrt(pHistoryb(:,1)),'.-');
plot(tHistory,sqrt(pHistoryc(:,1)),'.-');
hold off
ylim('padded')
ylabel('\surdP(1,1)')
xlabel('Times (s)')

nexttile;
hold on
plot(tHistory,sqrt(pHistoryb(:,2)),'.-');
plot(tHistory,sqrt(pHistoryc(:,2)),'.-');
hold off
ylim('padded')
ylabel('\surdP(2,2)')
xlabel('Times (s)')
title(tcl, 'RMSE Plots')
leg = legend('Position Meas', 'Velocity Meas','Orientation', 'horizontal');
leg.Layout.Tile = 'north';
