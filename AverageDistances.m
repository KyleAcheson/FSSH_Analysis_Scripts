% plot mean distances and std around mean. Can also plot PDFs directly from trajs.

clear all
close all

load('Filtered_Trajectories_Fixed.mat')
Q = combined_geoms;

[natom, ~, ntraj, nts] = size(Q);
dt = 0.5;
time = 0:dt:1000;
pulse = [150 0];
D = zeros(size(Q));
for traj=1:ntraj
        for ts=1:nts
            for a=1:3
                for b=a+1:3
                dd = norm(Q(a,1:3,traj,ts) - Q(b,1:3,traj,ts));
                D(a,b,traj,ts) = dd;
                D(b,a,traj,ts) = dd;
                end
            end
        end
end

DD = zeros(3,3,197,2001);
a = 0;
for traj=1:ntraj
    if ismember(traj,crash_index) ~= 1
       a = a + 1;
       DD(:,:,a,:) = D(:,:,traj,:);
    end
end
clear D
D = DD;
ntraj = 197;

angles = zeros(ntraj,nts);
ssdist = zeros(ntraj,nts);
csdist1 = zeros(ntraj,nts);
csdist2 = zeros(ntraj,nts);

for traj=1:ntraj
    d_end = squeeze(D(:,:,traj,end));
    cs1_end = d_end(2,1);
    cs2_end = d_end(3,1);
    for ts=1:nts
        d = squeeze(D(:,:,traj,ts));
        cs1 = d(2,1);
        cs2 = d(3,1);
        ss = d(3,2);
        costheta = (cs1^2+cs2^2-ss^2)/(2*cs1*cs2);
        ang = acos(costheta)*(180/pi);
        angles(traj,ts) = ang;
        ssdist(traj,ts) = ss;
        if cs1_end > cs2_end
        csdist1(traj,ts) = cs1;
        csdist2(traj,ts) = cs2;
        else
        csdist1(traj,ts) = cs2;
        csdist2(traj,ts) = cs1;
        end
    end
    end

angles(angles==0) = NaN;
avg_angle = sum(angles,1)./ntraj;
std_angle = std(angles,1,1);
x = time;
x2 = [x, fliplr(x)];

for ts=1:nts
    c1(ts) = avg_angle(ts) + std_angle(ts);
    c2(ts) = avg_angle(ts) - std_angle(ts);
end

InBetween = [c2, fliplr(c1)];

figure
fill(x2, InBetween, 'b', 'FaceAlpha', 0.2);
hold on
plot(time, avg_angle, 'r', 'LineWidth', 2);
xlim([0, 1000])
ylabel(['SCS Angle (Degrees)'], 'interpreter','latex')
xlabel('Time (fs)')
legend('Standard Deviation', 'Mean')

csdist1(csdist1==0) = NaN;
avg_cs1 = sum(csdist1,1)./ntraj;
std_cs1 = std(csdist1,1,1);
min_cs1 = min(csdist1,[],1);
x = time;
x2 = [x, fliplr(x)];

for ts=1:nts
    c1(ts) = avg_cs1(ts) + std_cs1(ts);
    c2(ts) = min_cs1(ts); %avg_cs1(ts) + std_cs1(ts);
end

InBetween = [c2, fliplr(c1)];

figure
fill(x2, InBetween, 'b', 'FaceAlpha', 0.2);
hold on;
plot(time, avg_cs1, 'r', 'LineWidth', 2);
xlim([0, 1000])
ylabel(['C-S Bond Length (\AA)'], 'interpreter','latex')
xlabel('Time (fs)')
legend('Standard Deviation', 'Mean')

csdist2(csdist2==0) = NaN;
avg_cs2 = sum(csdist2,1)./ntraj;
std_cs2 = std(csdist2,1,1);
x = time;
x2 = [x, fliplr(x)];

for ts=1:nts
    c1(ts) = avg_cs2(ts) + std_cs2(ts);
    c2(ts) = avg_cs2(ts) - std_cs2(ts);
end

InBetween = [c2, fliplr(c1)];

figure
fill(x2, InBetween, 'b', 'FaceAlpha', 0.2);
hold on;
plot(time, avg_cs2, 'r', 'LineWidth', 2);
xlim([0, 1000])
ylabel(['C-S Bond Length (\AA)'], 'interpreter','latex')
xlabel('Time (fs)')
legend('Standard Deviation', 'Mean')
%%
avg_angle = sum(angles,1)./ntraj;
avg_cs1 = sum(csdist1,1)./ntraj;
avg_cs2 = sum(csdist2,1)./ntraj;

min_angle = min(angles,[],1);
max_angle = max(angles,[],1);
min_cs1 = min(csdist1,[],1);
max_cs1 = max(csdist1,[],1);
min_cs2 = min(csdist2,[],1);
max_cs2 = max(csdist2,[],1);

figure
plot(time, avg_angle, 'r')
hold on
plot(time, min_angle, time, max_angle)

figure
plot(time, avg_cs1, 'r')
hold on
plot(time, min_cs1, time, max_cs1);

figure
plot(time, avg_cs2, 'r')
hold on
plot(time, min_cs2, time, max_cs2);





%%
pulse = [150 0];
bin_edge = linspace(0,20,2000);
lbm = length(bin_edge)-1;
lb = length(bin_edge);
NC1 = zeros(lbm, nts);
NC2 = zeros(lbm, nts);
NSS = zeros(lbm, nts);
C1_edges = zeros(lb, nts);
C2_edges = zeros(lb, nts);
SS_edges = zeros(lb, nts);

for ts=1:nts
[NC1(:,ts), C1_edges(:,ts)] = histcounts(csdist1(:,ts), bin_edge);
[NC2(:,ts), C2_edges(:,ts)] = histcounts(csdist2(:,ts), bin_edge);
[NSS(:,ts), SS_edges(:,ts)] = histcounts(ssdist(:,ts), bin_edge);
end

[DG, TG] = meshgrid(bin_edge(1:end-1), time);

figure
mesh(DG, TG, NC1')
hold on
mesh(DG, TG, NC2')
mesh(DG, TG, NSS')
colormap(parula)
caxis([0 10])
view(0,90)

%%

for ts=1:nts
    dNC1(:,ts) = (NC1(:,ts) - NC1(:,1));
    dNC2(:,ts) = (NC2(:,ts) - NC2(:,1));
    dNSS(:,ts) = (NSS(:,ts) - NSS(:,1));
end

dNC1(isnan(dNC1)) = 0;
dNC2(isnan(dNC2)) = 0;
dNSS(isnan(NSS)) = 0;

[DG, TG] = meshgrid(bin_edge(1:end-1), time);

figure
mesh(DG, TG, dNC1')
hold on
mesh(DG, TG, dNC2')
mesh(DG, TG, dNSS')
colormap(parula)
xlim([0 8])
caxis([0 10])
view(0,90)

%%

[tc,conv_NC1,ntc,tind] = convolution_signal(time,NC1,pulse,dt,nts,lbm);
[tc,conv_NC2,ntc,tind] = convolution_signal(time,NC2,pulse,dt,nts,lbm);
[tc,conv_NSS,ntc,tind] = convolution_signal(time,NSS,pulse,dt,nts,lbm);

%%

[DG, TG] = meshgrid(bin_edge(1:end-1), time);

figure
mesh(DG, TG, NC1')
hold on
mesh(DG, TG, NC2')
mesh(DG, TG, NSS')
colormap(parula)
caxis([0 10])
xlim([0 8])
view(0,90)
xlabel('$\mathit{r}$ (\r{A})', 'interpreter','latex')
ylabel('time (fs)')

[DCG, TCG] = meshgrid(bin_edge(1:end-1), tc);

figure
mesh(DCG, TCG, conv_NC1')
hold on
mesh(DCG, TCG, conv_NC2')
mesh(DCG, TCG, conv_NSS')
colormap(parula)
caxis([0 10])
xlim([0 8])
view(0,90)
xlabel('$\mathit{r}$ (\r{A})', 'interpreter','latex')
ylabel('time (fs)')

fconv_NC1 = conv_NC1(:, tind(1):tind(2));
fconv_NC2 = conv_NC2(:, tind(1):tind(2));
fconv_NSS = conv_NSS(:, tind(1):tind(2));

[DCG, TCG] = meshgrid(bin_edge(1:end-1), time);

figure
mesh(DCG, TCG, fconv_NC1')
hold on
mesh(DCG, TCG, fconv_NC2')
mesh(DCG, TCG, fconv_NSS')
colormap(parula)
caxis([0 10])
xlim([0 8])
view(0,90)
ylabel('time (fs)')
xlabel('$\mathit{r}$ (\r{A})', 'interpreter','latex')


