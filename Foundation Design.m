close all
clear
clc
%%%%%%%%%% Determines Weight and dimensions of Foundation
%Inputs for Area 1

Y = input ('Unit weight of Foundation ='); %lbs/ft^3
W1 = input ('Width of Foundation = '); %ft
L1 = input ('Length of Foundation ='); %ft
H1 = input ('Height of Foundation ='); %ft
V1 = W1 * L1 * H1; %Volume of cube foundation in ft^3

%Inputs for Area 2

W2 = input ('Width of Step 1 ='); %ft
L2 = input ('Length of Step 1 ='); %ft
H2 = input ('Height of Step 1 ='); %ft
V2 = W2 * L2 * H2; %ft^3

%Calculate total weight of foundation

weight = (V1 + V2) * Y; %lbs

%%%%%%%%%% Coulomb's Theory
% Determines horizontal earth pressures for active
% conditions using Coulomb's Theory on left side of foundation which will consist of clay. 
% However, sand and gravel is prefered Coulomb's Theory is accompanied by a variety of assumptions. 
%An additional assumption that we are making is that the soil behind the abutment is not sloped and actually horiontal
% to the x-axis. This assumption will actually incorporate in itself a
% factor of safety because we are now assuming that there is more soil
% behind the abutment then there actually is.
% To begin, we will first figure out what configuration of theta and phi produces
% the largest tangent value in the Pa equation. All other variables are held constant.
% Pa(i,n) = (w + Q) / (sind(90-ii+delta) + cosd(90-ii+delta) * tand(90-theta_f+phi));

%Define Variables

    %theta_f = failure angle of soil
    %delta = internal friction angle between wall and soil
    %phi = angle of internal friction
    %Pa = Resultant Earth Force
Q = 0; %surcharge load
ii = 90; %angle of foundation wall from horizontal
i = 0; %Iteration Number
n = 0; %Iteration Number
k = 0; %Iteration Number 
lambda_sat = 71.520; %Unit weight of gravel in lb/ft^3 obtained from measuring on sight http://www.geotechnicalinfo.com/soil_unit_weight.html
%Define Anonymous Fuctions

cosine = @(delta)cosd(90-ii+delta);
sine = @(delta)sind(90-ii+delta);
soil_weight = @(theta_f)((0.5*(H1^2))/tand(theta_f))*lambda_sat;
%We start at 36 degrees because we cannot have an angle in tangent function
%over 90 degrees.

%Calculates maximum Pa values using delta = 14 degrees. Phi and theta_f are the variables.
for delta = 29:1:31
    k = k+1;
    n = 0;
for phi = 32:1:44 %values obtained at http://www.geotechdata.info/parameter/angle-of-friction.html
    i = 0;
    n = n+1;
    for theta_f = 36:1:90
        i = i+1;
        c(i,n) = cosine(delta);
        s(i,n) = sine(delta);
        sw(i,n) = soil_weight(theta_f);
        tangent(i,n) = tand(90-theta_f+phi);
    end
end
Pa= sw ./ (s + (c.*tangent)); %lbs/ft
[Max_Pa(:,k),I] = max(Pa(:));
disp('Angle of Friction between Wall and Soil (delta) =');
disp(delta);
disp('Maximum Resultant Earth Pressure (Pa) ='); 
disp(Max_Pa(:,k));
[I_row,I_column] = ind2sub(size(Pa),I);
disp('I_row =');
disp(I_row);
disp('I_column =');
disp(I_column);
printmat(Pa,'Maximum Pa Values','(Row1)36 (Row2)37 (Row3)38 (Row4)39 (Row5)40 (Row6)41 (Row7)42 (Row8)43 (Row9)44 (Row10)45 (Row11)46 (Row12)47 (Row13)48 (Row14)49 (Row15)50 (Row16)51 (Row17)52 (Row18)53 (Row19)54 (Row20)55 (Row21)56 (Row22)57 (Row23)58 (Row24)59 (Row25)60 (Row26)61 (Row27)62 (Row28)63 (Row29)64 (Row30)65 (Row31)66 (Row32)67 (Row33)68 (Row34)69 (Row35)70 (Row36)71 (Row37)72 (Row38)73 (Row39)74 (Row40)75 (Row41)76 (Row42)77 (Row43)78 (Row44)79 (Row45)80 (Row46)81 (Row47)82 (Row48)83 (Row49)84 (Row50)85 (Row51)86 (Row52)87 (Row53)88 (Row54)89 (Row55)90','(Colum1)32 (Colum2)33 (Colum3)34 (Colum4)35 (Colum5)36 (Colum6)37 (Colum7)38 (Colum8)39 (Colum9)40 (Colum10)41 (Colum11)42 (Colum12)43 (Colum13)44');
end 
maximum_Pa = max(Max_Pa(:)); 
disp('Design Pa Value in lbs/ft');
disp(maximum_Pa);
delta_final = input('Value of delta that corresponds to Design Pa Value =');
PaH = maximum_Pa*cosd(90 - ii + delta_final);%lbs/ft
disp('Horizontal Earth Pressure lbs/ft =');
disp(PaH);

%%%%%%%%%% Determine Xn, e, and N

%Step 1: Split the L-shaped foundation into 2 seperate sections and
%calculate the weight, lever arm, and moment for each section.

    %Step 1.1: Calculates for area 1
    
    weight1 = (W1-W2)*(H1)*(Y); %lbs/ft
    lever_arm1 = W2 + ((W1-W2)/2); %ft
    Moment1 = weight1*lever_arm1; %lbs*ft/ft
    
    %Step 1.2: Calculates for area 2
    
    weight2 = W2*(H1-H2)*Y; %lbs/ft
    lever_arm2 = W2/2; %ft
    Moment2 = weight2*lever_arm2; %lbs*ft/ft
    
    %Step 1.3: Determine the normal force and sum of the moments
    
    N = weight1 + weight2; %lbs/ft
    M = Moment1 + Moment2; %lbs/ft/ft
    
%Step 2: Calculate Xn, set up sum of moments equal to 0. Counter-clockwise
%positive

    Xn = (M - (PaH*(H1/3)))/N; %ft
    
%Step 3: Calculate eccentricty 

    e = (W1/2) - Xn; %ft
    
%Step 4: Check for bearing capacity failure. Because we are not sure whether or not our
%foundation lies on a rigid or plastic soil, we will check the bearing
%capacity failure for both conditions. 
disp('Assess the stability of the foundation against bearing capacity failure');


    %Step 4.1: Assume Rigid Foundation Soil
    
    if e < W1/6
        qmax = (N/W1) + ((6*e*N) / (W1^2)); %psf
        disp('Rigid Foundation Condition'); 
        disp('Case 1');
        disp('Entire footing is under compressive stress with qmin > 0')
    elseif e == W1/6
        qmax = (2*N)/(W1); %psf
        disp('Rigid Foundation Condition'); 
        disp('Case 2');
        disp('Entire Footing is under compression stress with qmin = 0');
    else 
        Wo = 3*Xn;
        qmax = (2*N)/(Wo); %psf
        disp('Rigid Foundation Condition'); 
        disp('Case 3, e > W/6');
        disp('qmin < 0, Ignore tension and assume Wo = 3Xn')
    end 
    
    %Step 4.2: Plastic Foundation 
    disp('Plastic Foundation Condition');
    qmax1 = (N)/(2*Xn); %psf
    
    %Step 4.3: Determine whether the plastic foundation or the rigid foundation
    %assumption yields a larger qmax value. 
    
    if qmax > qmax1
        qmax_final = qmax; %psf
    else 
        qmax_final = qmax1; %psf
    end 
    disp('Maximum pressure applied on the soil from the foundation (qmax) =');
    disp(qmax_final);
   
    %Step 4.4: Determine if foundation is considered shallow or deep
    
    Df = H1;
    if Df <= (8*W1*L1)/(W1+L1)
        disp('Foundation is Shallow');
    else 
        disp('Foundation is Deep and the following calculations do not apply');
    end 
    
    %Step 4.5: Look at Lesson 18 Study Guide to determine the N_gamma and Nq values
    
    disp('Look at Lesson 18 Study Guide. Use the table on the last page, along with the angle of internal friction, to determine N and N_gamma');
    N_gamma = input('N_gamma =');
    Nq = input('Nq =');
    
    %Step 4.6: Use clayey sand cohesion value from http://www.geotechdata.info/parameter/cohesion.html
    
    C = 105; %Units in lbs/ft^2 (Cohesion value) 
    
    %Step 4.7: Determine "Nc" value from using the axial load from bridge
    
    axial_force = input('Axial Force in lbf ='); %lbf
    angle = input('Angle of Axial Force Measured from X-Axis ='); %degrees
    V1 = axial_force*(sind(angle)); %lbf
    H1 = axial_force*(cosd(angle)); %lbf
    
    if Df/W1 <= 2.5
        Nc = 5 * (1 - (1.3 * (H1/V1))) * (1 + (0.2 * (W1/L1))) * (1 + ((0.2 * Df)/(W1)));
    else 
        Nc = 7 * (1 - (1.3 * (H1/V1))) * (1 + (0.2 * (W1/L1)));
    end 
    
    %Step 4.8: Calculate the ultimate bearing capacity "qult"
    
    qult = (C*Nc) + (0.5*W1*lambda_sat*N_gamma) + (lambda_sat*Df*(Nq-1)); 
    
    %Step 4.9: Check the factor of safety against bearing capacity failure
    
    FS_B = qult/qmax_final;
    if FS_B >=3
        disp('Bearing capacity failure will not occur')
    else 
        disp('Bearing capcaity failure could occur')
    end 
    
%Step 5: Assess the stability of foundation against overturning
disp('Assess the stability of foundation against overturning');
if Xn/W1 >= 1/3
    disp('Foundation will not overturn')
else 
    disp('Foundation may overturn')
end 

%Step 6: Assess the stability of foundation against sliding
disp('Assess the stability of foundation against sliding');
S = N * tand(delta_final);
T = PaH;
FS_S = S/T;

if FS_S >= 1.5
    disp('Foundation will not slide');
else
    disp('Foundation may slide');
end 


    
    
    
    
    
    
    
    
    
    
    
