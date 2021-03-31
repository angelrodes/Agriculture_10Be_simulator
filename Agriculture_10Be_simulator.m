clear
close all hidden

P=[4,4*0.01,4*0.008];
L=[160,1000,500];
l=log(2)/1.38e6;

dz=1; % cm
dt=200; % a
erosion_rate=(dz)/(dt);
disp(['erosion_rate=' num2str(erosion_rate*1e4) 'mm/ka'])
soil_production_rate=erosion_rate;



t_max=201e3;
t=[0:dt:t_max];

rho_bedrock=2.7;
rho_saprolite=2.2;
rho_soil=1.3;
z=(0:dz:500+t_max*dz/dt)';
rho_profile=(rho_bedrock*(z>300)+rho_saprolite*(z<=300&z>150)+rho_soil*(z<=150)).*(z>=0);
effective_depth=cumsum(dz*rho_profile);

lowering=0;

% start figure
set(0,'defaultfigurecolor',[1 1 1])
figure('units','normalized','outerposition',[0 0 0.3 1])

  % create the video writer with 1 fps
  writerObj = VideoWriter('myVideo.avi');
  writerObj.FrameRate = 10;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video


max_x=sum(P./(l+rho_soil*erosion_rate./L))*1.1;
C=0.*z;

% generate quartz grains
numparticles=100000;
particles_z=rand(1,numparticles)*max(z);
particles_x=rand(1,numparticles)*max_x;
particles_value=rand(1,numparticles)*max(rho_profile);
sel_part=particles_value<interp1(z,rho_profile,particles_z,'nearest');
plot(particles_x(sel_part),particles_z(sel_part),'.','Color',[0.7 0.7 0.7],'MarkerSize',10)

for n=1:numel(t)
    
    % erase previus
    plot(100000,0,'.k')
    hold on
    
    
    
    if t(n)<50e3
        rho_profile=rho_bedrock.*(z>=0);
        mixing_depth=0;
        title_string=['Bedrock'];
    else
        rho_profile=(rho_bedrock*(z>300)+rho_saprolite*(z<=300&z>150)+rho_soil*(z<=150)).*(z>=0);
        if t(n)<100e3
            mixing_depth=0;
            title_string=['Soil'];
        elseif t(n)<150e3
            mixing_depth=20;
            title_string=['Soil + bioturbation'];
        elseif t(n)<200e3
            mixing_depth=50;
            title_string=['Sustainable agriculture'];
        else
            lowering=lowering+50*dt/1000;
            mixing_depth=50+lowering;
            erosion_rate=(dz+lowering)/(dt);
            rho_profile=(rho_bedrock*(z>300)+rho_saprolite*(z<=300&z>150)+rho_soil*(z<=150)).*(z>=lowering);
            title_string=['Unustainable agriculture'];
        end
    end
    
    z=z-dz;
    
    % accumualtion
   effective_depth=cumsum(dz*rho_profile);
   C(effective_depth<=0)=0;
    
    % accumulation
    C=C*exp(-l*dt)+(...
        P(1)/l.*exp(-effective_depth./L(1)).*(1-exp(-l.*dt))+...
        P(2)/l.*exp(-effective_depth./L(2)).*(1-exp(-l.*dt))+...
        P(3)/l.*exp(-effective_depth./L(3)).*(1-exp(-l.*dt))...
        ).*(effective_depth>0);
    C(z<mixing_depth & effective_depth>0)=mean(C(z<mixing_depth & effective_depth>0)); % mixing
    
    
    
    % particles
    particles_z=particles_z-dz;
    % mixing
    particles_z(particles_z>0 & particles_z<mixing_depth)=mixing_depth*rand(size(find(particles_z>0 & particles_z<mixing_depth)));
    % apply density
    sel_part=particles_value<interp1(z,rho_profile,particles_z,'nearest');
    
    plot(particles_x(sel_part),particles_z(sel_part),'.','Color',[0.7 0.7 0.7],'MarkerSize',5)
    
    % plot Be-10 profile
    plot(C(effective_depth>0),z(effective_depth>0),'-b','LineWidth',3)
    
    % texts
%     if erosion_rate*1e4<1000
        text(max_x,lowering,['\bf \epsilon=' num2str(erosion_rate*1e4) 'mm/ka'],'HorizontalAlignment','Right','VerticalAlignment','top')
%     else
%         text(max_x,lowering,['\bf \epsilon=' num2str(erosion_rate*1e4/1000) 'mm/a'],'HorizontalAlignment','Right','VerticalAlignment','top')
%     end
    
    if interp1(z,rho_profile,225)~=interp1(z,rho_profile,375)
        text(max_x,95,['\bf soil \rho=' num2str(interp1(z,rho_profile,95)) 'g/cm^{2}'],'HorizontalAlignment','Right')
        text(max_x,225,['\bf saprolite \rho=' num2str(interp1(z,rho_profile,225)) 'g/cm^{2}'],'HorizontalAlignment','Right')
        text(max_x,300,['\bf SPR=' num2str(soil_production_rate*1e4) 'mm/ka'],'HorizontalAlignment','Right')
        text(max_x,375,['\bf bedrock \rho=' num2str(interp1(z,rho_profile,375)) 'g/cm^{2}'],'HorizontalAlignment','Right')
    else
        text(max_x,225,['\bf bedrock \rho=' num2str(interp1(z,rho_profile,225)) 'g/cm^{2}'],'HorizontalAlignment','Right')
    end
    text(max_x/2,480,['\bf t=' num2str((t(n))/1000,'%15.1f') ' ka'],'HorizontalAlignment','center')
    
    % beautify the figure
    ylim([0 500])
    set(gca, 'Ydir', 'reverse')
    set(gca, 'XAxisLocation', 'Top')
    xlabel('[^{10}Be]','Color','b')
    xlim([0 max_x])
    ylabel('cm')
    title(title_string)
    
    hold off
    
        
        F(n) = getframe(gcf);
        frame = F(n) ;    
        drawnow
        writeVideo(writerObj, frame);
    
    if erosion_rate>soil_production_rate
    pause(.1)
    for nv=1:10
        writeVideo(writerObj, frame);
    end
    end
end


pause(4)
for nv=1:25000/200
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
end

% generate apparent profile

app_er=110/1e4; % cm/a

% apparent density
app_rho=0.*effective_depth;
for k=1:numel(z)
    if rho_profile(k)<=0
        app_rho(k)=0;
        surface=k;
    else
        app_rho(k)=mean(rho_profile(surface+1:k));
    end
end
app_C=(...
    P(1)./(l+app_rho.*app_er./L(1)).*exp(-effective_depth./L(1)).*(1-exp(-(l+app_rho.*app_er./L(1)).*t(n)))+...
    P(2)./(l+app_rho.*app_er./L(2)).*exp(-effective_depth./L(2)).*(1-exp(-(l+app_rho.*app_er./L(2)).*t(n)))+...
    P(3)./(l+app_rho.*app_er./L(3)).*exp(-effective_depth./L(3)).*(1-exp(-(l+app_rho.*app_er./L(3)).*t(n)))...
    ).*(effective_depth>0);
% app_C(z<mixing_depth & effective_depth>0)=mean(app_C(z<mixing_depth & effective_depth>0)); % mixing

hold on
plot(app_C(effective_depth>0),z(effective_depth>0),'-r','LineWidth',1.5)

sample_idx=find(z>200,1);
plot(app_C(sample_idx),z(sample_idx),'pr','MarkerFaceColor','r',...
    'MarkerSize',15)
text(app_C(sample_idx),z(sample_idx),['\bf  Apparent SPR=' num2str(app_er*1e4) 'mm/ka'],'Color','r', 'FontSize', 12)

% beautify the figure
ylim([0 500])
set(gca, 'Ydir', 'reverse')
set(gca, 'XAxisLocation', 'Top')
xlabel('[^{10}Be]','Color','b')
xlim([0 max_x])
ylabel('cm')
title(title_string)

pause(4)
for nv=1:25000/200
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
end

% credits
 text(0,525,'\bf More info: www.angelrodes.com','Color','b', 'FontSize', 14)
 % beautify the figure
ylim([0 500])
set(gca, 'Ydir', 'reverse')
set(gca, 'XAxisLocation', 'Top')
xlabel('[^{10}Be]','Color','b')
xlim([0 max_x])
ylabel('cm')
title(title_string)

pause(8)
for nv=1:50000/200
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
end

 
close(writerObj);



