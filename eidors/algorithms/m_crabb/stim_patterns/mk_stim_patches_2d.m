function [stim, meas_sel]= mk_stim_patches_2d(n_elec,n_elec_stim,inj)
%Make patch stimulation for rings of electrodes. This is a square patch of
%current 1 and -1
%INPUT : n_elec : No. of electrodes/ring (this must be even)
%        n_elec_stim : The no. of electrodes used for a stimulation
%        inj : The injection pattern [0 i] (i is gap to electrode i.e. i=1
%        is adjacent style, and i=nelec/2 is opposite style)
%OUTPUT : stim : The stimulation pattern strucutre
%         meas_sel : Logical array of electrodes used to measure
%
%NOTE:
%A) This is in its infancy, we have no options, amplitude inputs
%B) Need to get meas_sel doing something properly
%C) Need option to only measure on non carrying electrodes

%Define no. of electrodes/patch stimulations
n_stim_patch=n_elec; nelec_tot=n_elec*n_rings;

%Arbitrary stimulation pattern (well need 'meas_current' options to
%generate the correct structure)
[stim,meas_sel]=mk_stim_patterns(n_elec,n_rings,[0 1],[0 1],{'meas_current'});

%Get gap from first current 1 electrode to first current -1 electrode
%on a given ring
gap=inj(2)+n_rings-1;

%Loop over stimulations and create patch pattern
for i=0:n_stim_patch-1
    %Grab the stimulation patterns
    stim_i=stim(1,i+1); stim_i.stim_pattern=zeros(nelec_tot,1); 
    
    %Loop through rings and add the correct current
    for k=0:n_elec_stim-1
            stim_i.stim_pattern(mod(k+i,nelec_tot)+1)=1;
            stim_i.stim_pattern(mod(k+i+gap,nelec_tot)+1)=-1;
    end
    
    %Add this back to the stim structrue
    stim(1,i+1)=stim_i;
    
end

%Pick only the proper stimulations
stim=stim(1,1:n_stim_patch);

end