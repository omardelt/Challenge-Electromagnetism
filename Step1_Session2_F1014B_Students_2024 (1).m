%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Session 2 | Step 1              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rodrigo Gamboa & Francisco Montes | May 2024


%-OBJECTIVE-%

%In this activity you will learn to calculate the force between two sources of magnetic field.
%Specifically, you will calculate the force that a falling magnet (magnetic dipole) will experiment, 
%due to its interaction your former circular coil in the ground (with current I). 
% You will use and extend your code from the previous session.

%Recall commenting and completing all necessary code from now on,
%as this will impact the quality and grade of the associate Deliverable.


%-----------------------------------------------------------------------
%Code starts here (paste this work below your previous code viceversa)
%Complete the code in the left, using the comments in the right.
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------%%
%-----------------------------------------------------------------------%%
%-------------------------Dipole Falling Part---------------------------%%
%-----------------------------------------------------------------------%%
%-----------------------------------------------------------------------%%

%Reseach and explain (in your own words), what is magnetic moment in general, 
%and for a permanent magnet (simple equations are allowed and expected):
%__________________________________________________________________________
%__________________________________________________________________________
%__________________________________________________________________________
%__________________________________________________________________________
%__________________________________________________________________________
%__________________________________________________________________________


%New code starts here!


mag = 2000;               %Declare a variable called "mag" (Magnetic moment of the magnet),
                        %and initialize to 2000. What are its units %_____________?

mass = 0.004;                %Declare the Magnet´s mass in kg, use 0.004 to begin.
                        %How many grams is this_____________?

w = mass*-9.81;              %Calculate Magnet´s weight in N, and store it in a var. called "w". Don´t 
                        %forget the sign.


zo = 5;                   %Declare Magnet's initial position using a var. called "zo".

dt = 0.05;                %Declare a time step "dt" equal to 0.05.

zm(1) = zo;               %Declare a vector called "zm" (magnet position), and store in its first 
                        %value the Magnet's initial position var.

zmfree(1) = zo;           %Same as before, but using now a new vector called "zmfree", which stands
                        %for magnet´s position for the free fall case. Why
                        %do you think we introduced this
                        %_________________________________________________?


tt(1) = 0;                %Declare a vector called "tt" (time), and initialize it to zero.

vz(1) = 0;                %Declare a vector called "vz" (Z component for velocity), and initialize it to zero.

vzfree(1) = 0;            %Same as above but using now a vector labeled as "vzfree". Again, why _____________?

cc = 1;                   %Introduce a var., called "cc" (an index counter), and set it equal to one.



path = animatedline('Color','r','LineWidth',2);
%Use the animatedline function (use Matlab´s upper-right search function bar) to create an 
                    %animated line that has no data and adds it to the current axes. Later (a few lines below), add points to the line 
                    %in a loop using the addpoints function, add data.


 %<start  while>    %Start a while loop running until the position of the magnet (dynamically) in z 
                    %is bigger than -5 (i.e. the bottom of the canvas). 

                    %For diagnostics purposes, print the position of the magnet in z. 
while zm(cc) > -5

    disp(['Magnet z-position: ', num2str(zm(cc))]);




    %Using the addpoints function mentioned above, add points (the zm(cc) ones) 
                              %to the path, centered in x=0.

    addpoints(path,0,zm(cc));
                       %Use drawnow to update and/or modify graphics objects and 
                       %want to see the updates on your canvas immediately.
    drawnow
    
    head=scatter(0,zm(cc), 100,'filled'); %Uncomment this line and explain what its doing. 
   
    %This line creates a filled marker representing the magnet position.

    
    
   %delete(head)                %Leave this commented, uncomment 
                                %in next step.


   %%<end  while>             Close the while loop
end

%---------------------------------------------------------------------
%--Final Part---------------------------------------------------------
%---------------------------------------------------------------------

%Uncomment all of the following, and explain in 
%detail (in next step of the code, not now) what it does or will do. 

figure(4)
subplot(1,2,1)

hold on
plot(zm(1:length(Fm)),1000*Fm, '-b', 'LineWidth', 2) %Why do we have 1000Fm?
plot([0,0],[-150,150],'-.k','LineWidth', 2)
grid on
xlabel 'z position (m)'
ylabel 'Magnetic force (mN)'
title 'Magnetic force of a Current ring over a falling magnet'
legend('Magnetic force in the Z direction','Current loop location','Location','southeast')

subplot(1,2,2)

hold on
tt=0:dt:(cc-1)*dt;
plot(tt,zm,'-r', 'LineWidth', 2)
plot(tt,zmfree,'--b', 'LineWidth', 2)   
plot([0,1.8],[0,0],'-.k','LineWidth', 2)
grid on
xlabel 'time (s)'
ylabel 'z position (m)'
title 'Position vs time of a Magnetic dipole falling throug a current ring'
legend('Fall over a current ring','Free fall (no Magnetic force)', 'Current loop location','Location','southwest')
axis([0 1.8 -6 6])

