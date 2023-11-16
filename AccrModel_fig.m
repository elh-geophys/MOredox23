% EL 
% 11/9/2023
% Get pics of accretion models


model1 = 4;            %accretion models, 1 = W90(high), 2 = W90(low), 3 = H00, 4 = H04, 5 = N21
model2 = 5;            

[t, Accr_model1] = getAccrModel(model1);
[~, Accr_model2] = getAccrModel(model2);

figure(1);
subplot(1,2,1)
box on
plot(t, Accr_model1, 'k', "LineWidth", 1.5)
xticks([0 20 40 60 80 100])
yticks([0 0.2 0.4 0.6 0.8 1.0])
ax = gca;
ax.TickLength = [0.02 0];
ax.XTickLabelRotation = 0;
ax.Layer = 'top';
xlabel("Time (Myr)")
ylabel("Mass fraction")

subplot(1,2,2)
box on
plot(t, Accr_model2, 'k', "LineWidth", 1.5)
xticks([0 20 40 60 80 100])
yticks([0 0.2 0.4 0.6 0.8 1.0])
ax = gca;
ax.TickLength = [0.02 0];
ax.XTickLabelRotation = 0;
ax.Layer = 'top';
xlabel("Time (Myr)")
ylabel("Mass fraction")