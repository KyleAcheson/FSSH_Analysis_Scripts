clear all
close all

load('Filtered_Trajectories_Fixed.mat')
load('cslabels_2.mat')

geometries = combined_geoms;
[~, ~, ntraj, nts] = size(geometries);
tt = 0:0.5:1000;

Q = zeros(size(geometries));
for i=1:ntraj
    for t=1:nts
        if t < final_t(i)*2+2
            Q(:,:,i,t) = geometries(:,:,i,t);
        end
    end   
end

D = zeros(size(geometries));
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

angles = zeros(ntraj,nts);
ssdist = zeros(ntraj,nts);
csdist1 = zeros(ntraj,nts);
csdist2 = zeros(ntraj,nts);

for traj=1:ntraj
    for ts=1:nts
        d = squeeze(D(:,:,traj,ts));
        cs1 = d(2,1);
        cs2 = d(3,1);
        ss = d(3,2);
        costheta = (cs1^2+cs2^2-ss^2)/(2*cs1*cs2);
        ang = acos(costheta)*(180/pi);
        angles(traj,ts) = ang;
        ssdist(traj,ts) = ss;
        csdist1(traj,ts) = cs1;
        csdist2(traj,ts) = cs2;
    end
end
angles(angles==0) = NaN;

MinAnglesPT = min(angles, [], 2);
MinSSDist = zeros(213,1);
MinCSDist1 = zeros(213,1);
MinCSDist2 = zeros(213,2);
for i=1:213
    x = find(angles(i,:) == MinAnglesPT(i));
    MinSSDist(i) = ssdist(i,x);
    MinCSDist1(i) = csdist1(i,x);
    MinCSDist2(i) = csdist2(i,x);
end

count = 0;
for i=1:213
    if cslabels(i) == 0
        count = count + 1;
        MinAngF(count, 1) = MinAnglesPT(i,1);
        MinSSDistF(count, 1) = MinSSDist(i,1);
        MinCSDist1F(count, 1) = MinCSDist1(i,1);
        MinCSDist2F(count, 1) = MinCSDist2(i,1);
    end
end


%%
figure
subplot(2,2,1)
scatter(MinAngF, MinSSDistF)
xlabel('Minimum $\theta_{\textrm{scs}}$ (Degrees)', 'interpreter','latex')
ylabel('$\mathit{r}_{\textrm{ss}}$ (\r{A})', 'interpreter','latex')
subplot(2,2,2)
scatter(MinAngF, MinCSDist1F)
xlabel('Minimum $\theta_{\textrm{scs}}$ (Degrees)', 'interpreter','latex')
ylabel('$\mathit{r}_{\textrm{cs}}$ (\r{A})', 'interpreter','latex')
subplot(2,2,3)
scatter(MinAngF, MinSSDistF)
xlabel('Minimum $\theta_{\textrm{scs}}$ (Degrees)', 'interpreter','latex')
ylabel('$\mathit{r}_{\textrm{ss}}$ (\r{A})', 'interpreter','latex')
subplot(2,2,4)
scatter(MinAngF, MinCSDist2F)
xlabel('Minimum $\theta_{\textrm{scs}}$ (Degrees)', 'interpreter','latex')
ylabel('$\mathit{r}_{\textrm{cs}}$ (\r{A})', 'interpreter','latex')

%%
figure
for i=1:ntraj
    if cslabels(i) == 1
        subplot(1,2,1)
        scatter(angles(i,:), ssdist(i,:));
        xlabel('$\theta_{\textrm{scs}}$ (Degrees)', 'interpreter','latex')
        ylabel('$\mathit{r}_{\textrm{ss}}$ (\r{A})', 'interpreter','latex')
        hold on
        %xlim([0 140])
        subplot(1,2,2)
        scatter(angles(i,:), csdist1(i,:))
        hold on
        %xlim([0 140])
        xlabel('$\theta_{\textrm{scs}}$ (Degrees)', 'interpreter','latex')
        ylabel('$\mathit{r}_{\textrm{cs1}}$ (\r{A})', 'interpreter','latex')
    end
end

%%
figure
for i=1:ntraj
    if cslabels(i) == 3
        subplot(1,2,1)
        scatter(angles(i,:), ssdist(i,:));
        xlabel('$\theta_{\textrm{scs}}$ (Degrees)', 'interpreter','latex')
        ylabel('$\mathit{r}_{\textrm{ss}}$ (\r{A})', 'interpreter','latex')
        hold on
        %xlim([0 140])
        ylim([0 10])
        subplot(1,2,2)
        scatter(angles(i,:), csdist2(i,:))
        hold on
        %xlim([0 140])
        ylim([0 10])
        xlabel('$\theta_{\textrm{scs}}$ (Degrees)', 'interpreter','latex')
        ylabel('$\mathit{r}_{\textrm{cs2}}$ (\r{A})', 'interpreter','latex')
    end
end



%%
% figure
% for i=1:ntraj
%     hold on
%     plot(tt,angles(i,:))
% end
% 
% figure
% mean_angles = nanmean(angles,1);
% plot(tt,mean_angles, 'LineWidth', 4, 'color', 'r')
% xlim([0 500])

% count = 0
% for i=1:ntraj
%     for j=1:nts
%         if angles(i,j) < 120
%             count = count + 1;
%         end
%     end
% end
% 
% countb = 0;
% for i=1:ntraj
%     for j=1:nts
%         if angles(i,j) < 100
%             countb = countb + 1;
%         end
%     end
% end
% 
% countc = 0
% for i=1:ntraj
%     for j=1:nts
%         if angles(i,j) < 90
%             countc = countc + 1;
%         end
%     end
% end
% 
% countd = 0
% for i=1:ntraj
%     for j=1:nts
%         if angles(i,j) < 80
%             countd = countd + 1;
%         end
%     end
% end
% 
% under120 = (count/(ntraj*nts))*100;
% under100 = (countb/(ntraj*nts))*100;
% under90 = (countc/(ntraj*nts))*100;
% under80 = (countd/(ntraj*nts))*100;




    
    