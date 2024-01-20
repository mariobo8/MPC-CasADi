function Draw_MPC_path_trackin_carexample (t,xx,xx1,u_cl,xs,N, step_time, average_mpc_time, f2, cost_f)


set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)

line_width = 1.5;
fontsize_labels = 14;
%% compute path
theta_range = linspace(-30, 0, 300);
rho_val = full(f2(theta_range));
%--------------------------------------------------------------------------
%-----------------------Simulate robots -----------------------------------
%--------------------------------------------------------------------------
x_r_1 = [];
y_r_1 = [];


figure(500)

% Animate the robot motion
%figure;%('Position',[200 200 1280 720]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
set(gcf,'Units','normalized','OuterPosition',[0 0 0.55 1]);

for k = 1:size(xx,2)
    h_t = 0.5; w_t=0.3; % triangle parameters
    %plot path
    plot(theta_range, rho_val,'-k','linewidth',line_width);
    hold on
    x1 = xs(1); y1 = xs(2); th1 = xs(3); %final goal
    x1_tri = [ x1+h_t*cos(th1), x1+(w_t/2)*cos((pi/2)-th1), x1-(w_t/2)*cos((pi/2)-th1)];%,x1+(h_t/3)*cos(th1)];
    y1_tri = [ y1+h_t*sin(th1), y1-(w_t/2)*sin((pi/2)-th1), y1+(w_t/2)*sin((pi/2)-th1)];%,y1+(h_t/3)*sin(th1)];
    fill(x1_tri, y1_tri, 'g'); % plot final position
    hold on;
    x1 = xx(1,k,1); y1 = xx(2,k,1); th1 = xx(3,k,1); %state
    x_r_1 = [x_r_1 x1];
    y_r_1 = [y_r_1 y1];
    x1_tri = [ x1+h_t*cos(th1), x1+(w_t/2)*cos((pi/2)-th1), x1-(w_t/2)*cos((pi/2)-th1)];%,x1+(h_t/3)*cos(th1)];
    y1_tri = [ y1+h_t*sin(th1), y1-(w_t/2)*sin((pi/2)-th1), y1+(w_t/2)*sin((pi/2)-th1)];%,y1+(h_t/3)*sin(th1)];

    plot(x_r_1,y_r_1,'-r','linewidth',line_width);hold on % plot exhibited trajectory
    if k < size(xx,2) % plot prediction
        plot(xx1(1:N,1,k),xx1(1:N,2,k),'r--*')
    end
    
    fill(x1_tri, y1_tri, 'r'); % plot robot position

    
   
    hold off
    %figure(500)
    ylabel('$y$-position (m)','interpreter','latex','FontSize',fontsize_labels)
    xlabel('$x$-position (m)','interpreter','latex','FontSize',fontsize_labels)
    axis([-35 5 -5 7.5]) 
    pause(0.1)
    box on;
    grid on
    %aviobj = addframe(aviobj,gcf);
    drawnow
    % for video generation
    F(k) = getframe(gcf); % to get the current frame
end
close(gcf)
% viobj = close(aviobj)
% video = VideoWriter('exp.mp4','MPEG-4');
% video.FrameRate = 5;  % (frames per second) this number depends on the sampling time and the number of frames you have
% open(video)
% writeVideo(video,F)
% close (video)

figure
subplot(311)
stairs(t,u_cl(:,1),'k','linewidth',1.5); xlim([0 t(end)])
ylabel('v (rad/s)')
grid on
subplot(312)
stairs(t,u_cl(:,2),'k','linewidth',1.5); xlim([0 t(end)])
xlabel('time (seconds)')
ylabel('delta (rad/s)')
grid on
subplot(313)
stairs(t,u_cl(:,3),'r','linewidth',1.5); xlim([0 t(end)])
xlabel('time (seconds)')
ylabel('virtual_v (rad)')
grid on

figure
plot(t,step_time*1000,'linewidth',line_width);
hold on 
yline(average_mpc_time*1000,'linewidth',line_width);
xlabel('mpc step')
ylabel('solving time (ms)')
grid on

figure
plot(t, cost_f,'linewidth',line_width)
xlabel('mpc step')
ylabel('Objective function')
grid on


