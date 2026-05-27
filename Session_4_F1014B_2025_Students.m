clear
clc
clf

%% =====================================================================
%% PART 1: FIELD GRID DESIGN & SYSTEM CONFIGURATION
%% =====================================================================
nl=1;                                           % Number of loops in the solenoid
ds=0.1;                                         % Spatial step size for grids
x=(-5:ds:5); y=(-5:ds:5); z=(-3:ds:5);          % Coordinate axes vectors
Lx=length(x); Ly=length(y); Lz=length(z);       % Dimensions of spatial vectors
rw=0.1;                                         % Smoothing/softening factor to avoid singularity divisions
Icurr=30;                                       % Current flowing through the ring (Amperes)
mo=4*pi*1e-7;                                   % Permeability of free space
km=mo*Icurr/(4*pi);                             % Biot-Savart proportionality constant

% Find index intersections near zero crossings
cp=find(x<=0); cx=cp(1);
cp=find(y<=0); cy=cp(1);                        
cp=find(z<=0); cz=cp(1);

N=12;                                           % Number of discrete elements approximating the circle
R=1;                                            % Radius of the loop
sz=1;                                           % Vertical spacing parameter
cp=find(z>=-nl*sz*0.5); cz=cp(1);      
s=1;    
dtheta=2*pi/N;                                  % Angular step size per element
dl=R*dtheta;                                    % Length of each arc element
ang=(0:dtheta:2*pi-dtheta);                     % Angular vector around the unit circle

% Discretize the ring into spatial current element segments
for I=1:nl                           
    Px(s:s+N-1)=R*cos(ang);                     % X positions of current elements
    Py(s:s+N-1)=R*sin(ang);                     % Y positions of current elements
    Pz(s:s+N-1)=-nl/2*sz+(I-1)*sz;              % Z positions of current elements
    dx(s:s+N-1)=-Py(s:s+N-1)*dtheta;            % Differential length vector dx component
    dy(s:s+N-1)=Px(s:s+N-1)*dtheta;             % Differential length vector dy component
    s=s+N;                            
end
dz(1:N*nl)=0;                                   % Element length dz component (flat horizontal ring)

%% =====================================================================
%% DIAGNOSTIC VISUALIZATION: FIGURE 1
%% =====================================================================
figure(1)
hold on
quiver3(Px,Py,Pz,dx,dy,dz,0.5,'-r','LineWidth',3) % Render current flow directions via vectors
plot3(Px,Py,Pz,'*k','LineWidth',5)               % Emphasize coordinate discretization nodes
plot3(Px,Py,Pz,'-k','LineWidth',.3)              % Draw underlying path wire outline
xlabel 'x'; ylabel 'y'; zlabel 'z'; title 'Current elements'
legend('Discrete current elements','Current element location','Loop/Solenoid')
axis([-3 3 -3 3 -3 3]); grid on; view(-34,33)

%% =====================================================================
%% PART 2: BIOT-SAVART INTEGRATION NUMERICAL SOLVER
%% =====================================================================
dBx(1:Lx,1:Ly,1:Lz)=0; dBy=dBx; dBz=dBx;        % Preallocate magnetic field matrices
tic
for I=1:Lx
    for J=1:Ly                  
        for K=1:Lz
            for L=1:nl*N        
                rx=x(I)-Px(L);                  % Distance vector X component
                ry=y(J)-Py(L);                  % Distance vector Y component
                rz=z(K)-Pz(L);                  % Distance vector Z component
                r=sqrt(rx^2+ry^2+rz^2+rw^2);    % Magnitude with smoothing factor
                r3=r^3;                         % Cubed denominator cross-product term
                % Summing magnetic field contributions pointwise
                dBx(I,J,K)=dBx(I,J,K)+km*dy(L)*rz/r3;
                dBy(I,J,K)=dBy(I,J,K)+km*dx(L)*rz/r3;
                dBz(I,J,K)=dBz(I,J,K)+km*(dx(L)*ry-dy(L)*rx)/r3;
            end
        end
    end
end
toc

%% =====================================================================
%% PLOT MAGNETIC STREAMLINES: FIGURE 3
%% =====================================================================
Bmag=sqrt(dBx.^2+dBy.^2+dBz.^2);    
centery=round(Ly/2);                            % Find plane slice intersection center
Bx_xz=squeeze(dBx(:,centery,:));   
Bz_xz=squeeze(dBz(:,centery,:));    
Bxz=squeeze(Bmag(:,centery,:));    

figure(3)
hold on
pcolor(x,z,(Bxz').^(1/3)); shading interp; colormap jet; colorbar
h1=streamslice(x,z,Bx_xz',Bz_xz',1);            % Overlay field direction flow paths
set(h1,'Color', [0.8 1 0.9]);
xlabel 'x'; ylabel 'z'; title 'Magnetic field of a circular current'

%% =====================================================================
%% PART 3: DYNAMIC DIPOLE SIMULATION RUNTIME (ANIMATED)
%% =====================================================================
mag=10000;              % Magnetic moment of the falling magnet
m=0.004;                % Magnet's mass in kg
w=m*-9.81;              % Gravitational force (Weight vector in Newtons)
zo=5;                   % Initial elevation height in meters
dt=0.05;                % Computational sampling discrete time step
zm(1)=zo;               % Tracking vector array for actual magnetic fall
zmfree(1)=zo;           % Tracking vector array for reference uninhibited fall
tt(1)=0;                % Base clock initialization
vz(1)=0;                % Real velocity tracker
vzfree(1)=0;            % Free fall velocity tracker
cc=1;                   % General execution counter

figure(2)
path=animatedline('linewidth', 2, 'linestyle', ':', 'color', 'r');

while zm(cc)>0.0162
    addpoints(path,0,zm(cc));                   % Update path tracker graphics dynamically
    drawnow
    head=scatter(0,zm(cc), 100,'filled');       % Render physical marker tracking coordinates
   
    % Analytical calculations for the Magnetic force interaction profile
    Fm(cc)=(6*mo*Icurr*R^2*mag*(zm(cc)+0.00001))/(4*((zm(cc)+0.00001)^2+R^2)^(5/2));    
    F(cc)=Fm(cc)+w;                             % Sum forces acting on magnet
    a=F(cc)/m;                                  % Dynamic acceleration solver via Newton's 2nd Law
    pause(0.001)
    
    % Kinematic state integration calculations
    zm(cc+1)=zm(cc)+vz(cc)*dt+0.5*a*dt*dt;      % Calculate real position displacement
    zmfree(cc+1)=zmfree(cc)+vzfree(cc)*dt-0.5*9.81*dt*dt; % Calculate reference position displacement
    vz(cc+1)=(zm(cc+1)-zm(cc))/dt;              % Track physical velocity derivative step
    vzfree(cc+1)=(zmfree(cc+1)-zmfree(cc))/dt;  % Track reference free fall derivative step
    cc=cc+1;                                    % Increment loop steps array counters
    delete(head)                                % Wipe temporary locator dot graphic marker
end

%% =====================================================================
%% PART 4: SYSTEM KINETICS EVALUATIONS: FIGURE 4
%% =====================================================================
figure(4) % Opens a new figure window designated as Figure 4
subplot(1,2,1) % Sets up a 1-row, 2-column grid of plots and selects the 1st subplot (left side)
hold on % Retains existing plots on the current axes so new plots don't overwrite them
plot(zm(1:length(Fm)),1000*Fm, '-b', 'LineWidth', 2) % Scaled to milliNewtons (mN) for visualization clarity
plot([0,0],[-0.15,0.15],'-.k','LineWidth', 2) % Draws a vertical black dash-dot line at z=0 to mark the position of the current loop
grid on % Turns on the background grid lines for precise data reading
xlabel 'z position (m)' % Labels the horizontal axis as the magnet's position in meters
ylabel 'Magnetic force (mN)' % Labels the vertical axis as the Magnetic force experienced in milliNewtons
title 'Magnetic force of a Current ring over a falling magnet' % Assigns a descriptive title to the first subplot
legend('Magnetic force in the Z direction','Current loop location','Location','southeast') % Places a legend in the bottom-right corner
subplot(1,2,2) % Selects the 2nd subplot within the 1-row, 2-column layout (right side)
hold on % Retains existing plots on this second subplot axis
tt=0:dt:(cc-1)*dt; % Generates a time vector 'tt' spanning from 0 seconds to the final simulated time step
plot(tt,zm,'-r', 'LineWidth', 2) % Plots the real trajectory (Z-position vs time) of the magnet under magnetic influence as a solid red line
plot(tt,zmfree,'--b', 'LineWidth', 2)   % Plots the theoretical free-fall trajectory (no magnetic force) as a dashed blue line for comparison
plot([0,1.8],[0,0],'-.k','LineWidth', 2) % Draws a horizontal black dash-dot line at z=0 to visually indicate where the current ring sits over time
grid on % Turns on the grid lines for the second subplot
xlabel 'time (s)' % Labels the horizontal axis as elapsed time in seconds
ylabel 'z position (m)' % Labels the vertical axis as the vertical position in meters
title 'Position vs time of a Magnetic dipole falling throug a current ring' % Assigns a descriptive title to the second subplot
legend('Fall over a current ring','Free fall (no Magnetic force)', 'Current loop location','Location','southwest') % Places an explanatory legend in the bottom-left corner
axis([0 1.8 -6 6]) % Restricts the view window of the plot: time from 0 to 1.8s, and position from -6m to 6m

%-------------------------------%
% Challenge Session No. 4
%-------------------------------%
%START%
clear % Clears all variables from the workspace to start fresh
clc   % Clears the command window for clean readability
clf   % Clears the current figure window to remove previous drawings

mag = 500; % Define a "mag" variable to define the magnetic moment equal to 500
Rring = 0.5; % Define a "Rring" variable to define the ring radius in m (e.g., 0.5m)

% Define initial position of magnet (zo = 0.1).
zo = 0.5; 
% Define ring position = 0 (use "zring" variable).
zring = 0; 

dt = 0.0¿1; % Define time step ("dt").

% Define a time vector, and initialize to zero its first entry.
t(1) = 0; 

zm(1) = zo; % Define a Z-position vector, and initialize to zo its first entry.
cc = 1;     % Define, as in previous class, a "cc" counter.
vz(1) = 0;    % Define a vz vector; initialize its first entry to zero.

figure(1) % Call for figure(1).

%--------Part 2---------%
% Use [open and close] a while loop, using the zm variable vector, such that
% you can warranty that all z values for the magnet while falling and
% passing through the coil are covered.
while zm(cc) > (-zo)
    % Inside the while loop, use the following (uncomment), and
    % explain what it does or will do.
    pause(0.001) % Pauses execution for 1 ms to slow down the animation for human eyes
    clf          % Clears the figure window so the next frame updates without overlapping old ones
    
    [x, y, phiB1, Bz] = B_due_M(zm(cc), mag, Rring);
    % What is this? Explain in full detail, take your time!
    % It calls a helper function that calculates the spatial grid arrays (x, y),
    % the vertical magnetic field components (Bz) across the area of the ring, 
    % and the total initial magnetic flux (phiB1) passing through the ring 
    % at the magnet's current position zm(cc).
    
    % Using kinematics, the same as in last week's code actually, write an
    % expression for the zm position of the magnet, while free falling.
    zm(cc+1) = zm(cc) + vz(cc)*dt - 0.5*9.81*dt^2;
    
    % Using kinematics, the same as in last week's code actually, write an
    % expression for the vz position of the magnet, while free falling.
    vz(cc+1) = (zm(cc+1)-zm(cc))/dt;
    
    [x, y, phiB2, Bz] = B_due_M(zm(cc+1), mag, Rring);
    % What is this? Explain in full detail, take your time!
    % This recalculates the magnetic field profile and the new magnetic flux (phiB2) 
    % through the ring area at the magnet's NEXT predicted position zm(cc+1). 
    % We need both phiB1 and phiB2 to evaluate how fast the flux is changing over dt.
    
    %--------Part 4---------% fem CALCULATION | DECLARATION
    % Taking into account that the B_due_M function has been called
    % twice (why)? Calculate a vector function of fem(cc) equal to the
    % "rate change or difference of the magnetic flux, with respecto to
    % time", associated with the area of the ring (coil or metal loop).
    fem(cc) = -(phiB2 - phiB1) / dt;
    
    % Uncomment all of the following, and explain in details what each line is doing.
    subplot(2,2,1) % Divides the figure window into a 2x2 grid and selects the 1st plot (top-left)
    hold on        % Retains the current plot when new data points are added
    grid on        % Turns on the grid lines for easier measurement reading
    xlabel 'time, s'
    ylabel 'fem, mV'
    plot(t(1:cc), 100*fem(1:cc), '-k', 'LineWidth', 1) % Plots the continuous history of electromotive force (EMF) in black
    plot(t(1:cc), 100*fem(1:cc), '*r', 'LineWidth', 2) % Places red asterisks at the discrete calculated data points
    
    subplot(2,2,2) % Selects the 2nd plot window in the 2x2 grid (top-right)
    hold on        % Retains the current plot
    axis([0 0.3 -10 10]) % Constrains the axes limits (Time: 0 to 0.3s, Height: -10 to 10cm)
    grid on        % Turns on grid lines
    xlabel 'time, s'
    ylabel 'magnet heigth, cm'
    plot(t(1:cc), 100*zm(1:cc), 'ob', 'LineWidth', 2) % Plots the falling magnet's height over time as blue circles
    
    subplot(2,2,3) % Selects the 3rd plot window in the 2x2 grid (bottom-left)
    hold on        % Retains the current plot
    % Generates a flat 2D color map representing the density/strength of the magnetic field (Bz) passing through the ring.
    pcolor(x, y, zm(cc)/abs(zm(cc)) * abs(abs(0.005^2 * Bz)).^(1/3)); 
    shading interp; % Smooths the color transitions between calculation nodes
    colormap hot;   % Sets the visual theme to a fire-like hot scale (yellow/red/black)
    colorbar        % Displays a color scale bar showing value mappings
    view(-45, -45)  % Shifts the spatial observation camera perspective angle
    
    subplot(2,2,4) % Selects the 4th plot window in the 2x2 grid (bottom-right)
    hold on        % Retains the current plot
    % Constructs a 3D wireframe mesh representing the magnetic field intensity over the 2D surface area of the loop.
    mesh(x, y, zm(cc)/abs(zm(cc)) * abs(abs(10^2 * Bz)).^(1/3));
    view(-30, -3)   % Adjusts the 3D rotation view parameters of the mesh plot
    axis([-Rring Rring -Rring Rring -5 15]) % Sets strict spatial borders around the ring size bounds
    
    cc = cc + 1;       % Increments the index counter by 1 for the next loop iteration
    t(cc) = t(cc-1) + dt; % Calculates and appends the next time step to our time vector array
end % close the while loop

%% =====================================================================
%% NESTED FUNCTION: MAGNETIC FIELD DUE TO DIPOLE MOMENT
%% =====================================================================
function [x,y,phiB,Bz]=B_due_M(z,mag,Rring)
    mo=4*pi*1e-7;
    ds=0.005;
    x=-Rring:ds:Rring; y=-Rring:ds:Rring;
    Lx=length(x); Ly=length(y);
    Bz(1:Lx,1:Ly)=0;
    phiB=0;
    for i=1:Lx
        for j=1:Ly
            r=sqrt(x(i)^2+y(j)^2);
            if r<Rring
                Bz(i,j)=mo/(4*pi)*(3*z*(mag*z)-mag*(x(i)^2+y(j)^2+z^2))/((x(i)^2+y(j)^2+z^2+(ds/10)^2)^(5/2));
                phiB=phiB+ds^2*Bz(i,j);
            end
        end
    end
end