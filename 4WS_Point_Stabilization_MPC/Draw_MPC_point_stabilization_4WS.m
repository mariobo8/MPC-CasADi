function Draw_MPC_point_stabilization_4WS (t,xx,xx1,u_cl,xs,N,rob_diam,beta)


set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)

line_width = 1.5;
fontsize_labels = 14;

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
    h_t = 0.14; w_t=0.09; % triangle parameters
    
    x1 = xs(1); y1 = xs(2); th1 = xs(3);
    x1_tri = [ x1+h_t*cos(th1), x1+(w_t/2)*cos((pi/2)-th1), x1-(w_t/2)*cos((pi/2)-th1)];%,x1+(h_t/3)*cos(th1)];
    y1_tri = [ y1+h_t*sin(th1), y1-(w_t/2)*sin((pi/2)-th1), y1+(w_t/2)*sin((pi/2)-th1)];%,y1+(h_t/3)*sin(th1)];
    fill(x1_tri, y1_tri, 'g'); % plot reference state
    hold on;
    x1 = xx(1,k,1); y1 = xx(2,k,1); th1 = xx(3,k,1);
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
    axis([-0.2 2.5 -0.2 2.5])
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
video = VideoWriter('4ws_point_stab.mp4','MPEG-4');
video.FrameRate = 10;  % (frames per second) this number depends on the sampling time and the number of frames you have
open(video)
writeVideo(video,F)
close (video)

figure
subplot(221)
stairs(t,u_cl(:,1),'k','linewidth',1.5); axis([0 t(end) -0.35 0.75])
ylabel('v_f (rad/s)')
grid on
subplot(222)
stairs(t,u_cl(:,2),'k','linewidth',1.5); axis([0 t(end) -0.85 0.85])
xlabel('time (seconds)')
ylabel('v_r (rad/s)')
grid on
subplot(223)
stairs(t,u_cl(:,3),'r','linewidth',1.5); axis([0 t(end) -0.85 0.85])
xlabel('time (seconds)')
ylabel('\delta_f (rad)')
grid on
subplot(224)
stairs(t,u_cl(:,4),'r','linewidth',1.5); axis([0 t(end) -0.85 0.85])
xlabel('time (seconds)')
ylabel('\delta_r (rad)')
grid on
% subplot(325)
% stairs(t,u_cl(:,5),'r','linewidth',1.5); axis([0 t(end) -0.85 0.85])
% xlabel('time (seconds)')
% ylabel('alpha (rad)')
% grid on

figure
stairs(t, beta,'b','linewidth',1.5)
xlabel('time (seconds)')
ylabel('\beta (rad)')
grid on
