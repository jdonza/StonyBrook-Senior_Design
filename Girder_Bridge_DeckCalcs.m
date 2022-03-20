close all
clear
clc

%%% This code calculates the tensile and compressive stresses that
%%% are experienced within a bamboo culm when loaded in a 3-point
%%% bending scenario.

%Calculates Maximum Moment on an individually Bamboo Deck Member

L = input('Width of Deck in feet =');
l = L*0.3048; %Converts feet into meters
F = input('Point Load Acting in Center of Beam in lbf =');
f = F*4.44822;
moment = (f*(l/2))*1000;

%Calculates Tensile and Compressive Stresses in Bamboo Deck Member

answer = questdlg('Are the Deck Members Split in Half?','User Question','Yes','No','No');

switch answer
    case 'Yes'
        
        %Draws a Bamboo Deck Member Cross-Section
        
        R = 41.98055556;%input('Outer Radius of Bamboo in mm =');
        r = 33.16111111;%input('Inner Radius of Bamboo in mm =');
        th = linspace(0,pi,100);
        x = R*cos(th);
        y = R*sin(th);
        figure
        plot(x,y,'b');
        axis equal;
        ylabel('Millimeters');
        xlabel('Millimeters');
        title('Bamboo Deck Cross-Section');
        
        hold on
        x1 = r*cos(th);
        y1 = r*sin(th);
        plot(x1,y1,'b');
        
        hold on
        x2 = linspace(R,r,100);
        y2 = 0*x2;
        plot(x2,y2,'b');
        
        hold on
        x3 = linspace(-R,-r,100);
        y3 = 0*x3;
        plot(x3,y3,'b');
        hold off
        
        %Calculates Stress in Bamboo Deck Member
        
        Inertia = (0.1098*(R^4-r^4)) - ((0.283*(R^2*r^2)*(R-r))/(R+r));
        y_bar = (4*(R^3-r^3))/(3*pi*(R^2-r^2));
        Cc = R-y_bar;
        Ct = y_bar;
        sigma_T = (Ct*moment)/Inertia;
        sigma_C = (Cc*moment)/Inertia;
        disp('Tensile Stress in N/mm^2 =');
        disp(sigma_T);
        disp('Compressive Stress in N/mm^2 =');
        disp(sigma_C);
        
    case 'No'
        
        % Draws the Cross-Section of the Bamboo Deck Member
        
        R = 41.98055556;%input('Outer Radius of Bamboo in mm =');
        r = 33.16111111;%input('Inner Radius of Bamboo in mm =');
        center = [0,0];
        figure
        viscircles(center,R);
        ylabel('Millimeters');
        xlabel('Millimeters');
        title('Bamboo Deck Cross-Section');
        
        hold on
        viscircles(center,r);
        hold off
        
        %Calculates the Stresses in Bamboo Member
        
        Inertia = (pi/4)*((R^4)-(r^4));
        Ct = R;
        Ct = R;
        sigma_T = (Ct*moment)/Inertia;
        sigma_C = -1*(Ct*moment)/Inertia;
        disp('Tensile Stress in N/mm^2=');
        disp(sigma_T);
        disp('Compressive Stress in N/mm^2=');
        disp(sigma_C);
end