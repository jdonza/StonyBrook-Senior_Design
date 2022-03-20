close all
clear
clc
%This code operates under the assumption that the bamboo members are not
%hollow. 
%%%%%%%%%%%% Live Load Determination

%Step 1: Define the length of the 2 members running in the positive y-axis
%direction. Also, define the width of the deck and the live load distribution.

L = input('Length of Longitudinal Members Across Bridge Span in meters = '); %7 meters + 2(1ft)
W = input('Width of the Decking CL to CL in meters = ');
live_load = input('Live Load Distribution in N/m^2 = ');

%Step 2: Determine the loading distribution (wL) across one of the longitudinal
%members.

wL = (live_load*W)/2;% N/m

%%%%%%%%%%%% Dead Load from Transverse Deck Members

%Step 1: Input Transverse Deck Member Dimensions. Assume that the bamboo
%members are non-hollow cylinders. This actually incorporates a factor of
%safety because now we are assuming that there is dead weight that there actually is.  

d = input('Outer Diameter of Transverse Deck Member Dimension in meters =');
rho = input('Density of Transverse Bamboo in kg/m^3 ='); 
w = input('Width of Cut Transverse Members in meters =');

%Step 2: Determine the loading distribution (wd) subjected to one half of
%the bridge

n = L/d; %Calculates total number of transverse members needed to span deck
v = pi*((d/2)^2)*w; %Calculates the volume of each member
V = v*n; %Calculates the total volume of all bamboo members
mass = rho*V; %Calculates the weight of all the transverse members in kg
T_weight = mass*9.80665; %Calculates the weight of the transverse members in newtons using conversion factor
p = (T_weight)/(L*w); %Calculates the Pressure Distribution due to the transverse members in N/m^2
wd = (p*w)/2; %Calculates the loading distribution across one side of the bridge in N/m

%%%%%%%%%%%% Dead Load from Longitudinal Bamboo Deck Members

%Step 1: Input Material Dimensions

D = input('Outer Diameter of Longitundal Bamboo Deck Members in meters =');%Assume non-hollow cylindrical section
N = input('Number of Longitudinal Members =');%If more than 2, must rewrite code
rho1 = input('Density of Longitudinal Bamboo in kg/m^3 =');
%Step 2: Determine the loading distribution (wd1) subjected to one half of
%the bridge

v1 = pi*((D/2)^2)*L; %Calculates volume of each member
V1 = v1*N; %Calculates total volume of longitudinal members
MASS = rho1*V1; %Calculates the total mass of the bamboo in kg
T1_weight = MASS*9.80665; %Calculates the total weight of the members
wd1 = (T1_weight/N)/L; %Calculates the load distribution

%%%%%%%%%%%% Calculate Factored Load Distribution

WD = wd+wd1; %Adds together the two dead loading conditions
Factored_Loading = (1.2*WD) + (1.6*wL); % N/m
disp('The factored loading condition ='); %N/m
disp(Factored_Loading);


