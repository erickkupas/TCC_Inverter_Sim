data = load('PR4.txt');
%data = table(uk,vok,rk);
data2 = load('PRAv4.txt');

close all;

%% Plots

subplot(3,1,1)
plot(data(:,3)*1000,data(:,2),'-k','linewidth',1.2); grid on; hold on;
plot(data2(:,3)*1000,data2(:,2),'-b','linewidth',1.2); grid on; hold on;
plot(data(:,3)*1000,data(:,4),'--r','linewidth',1.5);
ylim([-300,300]);
%xlim([0.053,0.115]);
xlabel('Tempo [ms]')
legend('v_{o PR}(k)','v_{o PR Av}(k)','r(k)','location','best')
%legend('v_{o PR Av}(k)','r(k)','location','best')
ylabel('Tensão [V]')
set(gca,'fontname','times') 
set(gca,'fontsize',10)

subplot(3,1,2)
plot(data(:,3)*1000,data(:,1),'-k','linewidth',1.2); grid on; hold on;
plot(data2(:,3)*1000,data2(:,1),'-b','linewidth',1.5); grid on; hold on;
%xlim([0,0.06]);
ylim([-0.7,0.7]);
xlabel('Tempo [ms]')
ylabel('Ação de controle')
legend('u_{PR}(k)','u_{PR Av}(k)','location','best')
%legend('u_{PR Av}(k)','location','best')
set(gca,'fontname','times') 
set(gca,'fontsize',10)

e1 = abs(data(:,4)-data(:,2));

% stop = 1;
% 
% for i=1:length(e1)
%     if e1(i)>1000 & stop ==1
%         aux = i;
%         stop = 0;
%     end
% end
% 
% e1(aux:end) = 200;

e2 = abs(data2(:,4)-data2(:,2));

subplot(3,1,3)
plot(data(:,3)*1000,e1,'-k','linewidth',1.2); grid on; hold on;
plot(data2(:,3)*1000,e2,'-b','linewidth',1); grid on; hold on;
%xlim([0.053,0.115]);
xlabel('Tempo [ms]')
legend('e_{PR}(k)','e_{PR Av}(k)','location','best')
ylabel('Erro [V]')
ylim([0 50])
set(gca,'fontname','times') 
set(gca,'fontsize',10)
