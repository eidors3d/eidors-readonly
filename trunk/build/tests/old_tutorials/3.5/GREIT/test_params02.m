r = 0.05; % target radius
Npos = 20; % number of positions
Xpos =  linspace(0,0.9,Npos); % positions to simulated along x-axis
Ypos = zeros(1,Npos); 
Zpos = ones(1,Npos);  %% for off-plane, adjust the level (*1.5)
xyzr = [Xpos; Ypos; Zpos; r*ones(1,Npos)];

[vh,vi] = simulate_movement(imgs, xyzr);
