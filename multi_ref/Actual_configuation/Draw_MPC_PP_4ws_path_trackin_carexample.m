function Draw_MPC_PP_4ws_path_trackin_carexample (t,xx,xx1,u_cl,xs,N, step_time, average_mpc_time, cost_f, dim_error, x_p, y_p, beta, arc_length, x_r, a)


set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)

line_width = 1.5;
fontsize_labels = 14;
% %% compute path
% 
% %--------------------------------------------------------------------------
% %-----------------------Simulate robots -----------------------------------
% %--------------------------------------------------------------------------
% 
% x_r_1 = [];
% y_r_1 = [];
% % compute coordinates for s
% x_s = interp1(arc_length,x_p,xx1(:,4,:),'pchip');
% y_s = interp1(arc_length, y_p,xx1(:,4,:),'pchip');
% 
% figure(500)
% 
% % Animate the robot motion
% %figure;
% set(gcf,'Position',[100 100 1280 1280]);
% %set(gcf,'PaperPositionMode','auto')
% set(gcf, 'Color', 'w');
% %set(gcf,'Units','normalized','OuterPosition',[10 0 0.55 1]);
%     plot(x_p, y_p,'-k','linewidth',1);
%     hold on
%     % plot(xx(1,:), xx(2,:),'linewidth',1.5)
% 
% for k = 1:size(xx,2)
%     h_t = 0.2; w_t=0.1; % triangle parameters
%     %plot path
%     plot(x_p, y_p,'-k','linewidth',1);
%     hold on;
%     x1 = xx(1,k); y1 = xx(2,k); th1 = xx(3,k); %state
%     x_r_1 = [x_r_1 x1];
%     y_r_1 = [y_r_1 y1];
%     x1_tri = [ x1+h_t*cos(th1), x1+(w_t/2)*cos((pi/2)-th1), x1-(w_t/2)*cos((pi/2)-th1)];%,x1+(h_t/3)*cos(th1)];
%     y1_tri = [ y1+h_t*sin(th1), y1-(w_t/2)*sin((pi/2)-th1), y1+(w_t/2)*sin((pi/2)-th1)];%,y1+(h_t/3)*sin(th1)];
% 
%     plot(x_r_1,y_r_1,'-r','linewidth',line_width);hold on % plot exhibited trajectory
%     if k < size(xx,2) % plot prediction
%         plot(xx1(1:N,1,k),xx1(1:N,2,k),'r--*')
%         plot(x_r(:,1,k),x_r(:,2,k), 'g--*')
%     end
% 
%     hold on
%     fill(x1_tri, y1_tri, 'blue'); % plot robot position
% 
% 
% 
%     hold off
%     %figure(500)
%     ylabel('$y$-position (m)','interpreter','latex','FontSize',fontsize_labels)
%     xlabel('$x$-position (m)','interpreter','latex','FontSize',fontsize_labels)
%     legend('Path','Executed','Predicted', 'reference')
%     axis([-6 1 -1 4]) 
%     pause(0.001)
%     box on;
%     grid on
%     %aviobj = addframe(aviobj,gcf);
%     drawnow
%     % for video generation
%     F(k) = getframe(gcf); % to get the current frame
% end
% %viobj = close(aviobj)
% % close(gcf)
% % video = VideoWriter('carexample_path_tracking.avi','Uncompressed AVI');
% % video.FrameRate = 10;  % (frames per second) 
% % open(video)
% % writeVideo(video,F)
% % close (video)


figure
subplot(321)
stairs(t,u_cl(:,1),'k','linewidth',1.5); axis([0 t(end) 0 1.2])
ylabel('v_f (m/s)')
grid on
subplot(322)
stairs(t,u_cl(:,2),'k','linewidth',1.5); axis([0 t(end) 0 1.2])
xlabel('time (seconds)')
ylabel('v_r (rad)')
grid on
subplot(323)
stairs(t,u_cl(:,3),'r','linewidth',1.5); axis([0 t(end) -1.2 1.2])
xlabel('time (seconds)')
ylabel('\delta_f (rad)')
grid on
subplot(324)
stairs(t,u_cl(:,4),'r','linewidth',1.5); axis([0 t(end) -1.2 1.2])
xlabel('time (seconds)')
ylabel('\delta_r (rad)')
grid on
subplot(325)
stairs(t,u_cl(:,5),'r','linewidth',1.5); axis([0 t(end) -1.2 1.2])
xlabel('time (seconds)')
ylabel('\alpha (rad)')
grid on
subplot(326)
stairs(t,u_cl(:,6),'r','linewidth',1.5); axis([0 t(end) 0 1.2])
xlabel('time (seconds)')
ylabel('v_{virtual} (m/s)')
grid on


figure
evolution = cumtrapz(t,( u_cl(:,5))) + arc_length(1);
plot(t, evolution, "LineWidth",1.5)
xlabel('time (seconds)')
ylabel('path evolution (rad/s)')
grid on

figure
plot(1:length(step_time),step_time*1000,'linewidth',1);
hold on 
yline(average_mpc_time*1000,'r','linewidth',line_width);
xlabel('mpc step')
ylabel('solving time (ms)')
legend('','average')
xlim([0 length(t)])
grid on

figure
plot(cost_f,'linewidth',line_width)
xlabel('mpc step')
ylabel('Objective function')
xlim([0 length(t)])
grid on

figure
plot(t, dim_error,'linewidth',1.5)
xlabel('time (seconds)')
ylabel('Path-error (meters)')
xlim([0 t(end)])
grid on

figure
plot(t, beta,'linewidth',1.5)
xlabel('time (seconds)')
ylabel('\beta (rad)')
xlim([0 t(end)])
grid on


