# track_RAFT
Fast accurate particle tracking even for large particle displacements

% % track_RAFT is a piece of code that unscrambles n-dimensional trajectories
% from a list of particle coordinates determined at discrete times (e.g.
% consecutive video frames). This tracker has been designed to be
% interchangeable with the function, track.m, by John Crocker that is widely
% used for particle tracking. It works when you have many particles visible
% at every timepoint. It is not suitable for tracking a few individual
% particles.
% Benefits of this new code are:
%   - it is typically much faster (it includes a waitbar, so you can 
%     estimate how long you need to wait).
%   - It can work with very large displacements (much larger than the
%     average particle spacing in each frame), provided that clusters of
%     particles move together. This means that it is very good at tracking
%     tracer particles on, or in deforming solids. Also, it works well for
%     particles in fluid flow, when particles do not rearrange between 
%     frames (including standard diffusion experiments with the usual small
%     timesteps).
%   - It typically has far fewer obviously wrong tracks. The Hungarian
%     algorithm (which track.m and [1] uses) tries to match up as many 
%     particle pairs as possible. This function only matches up pairs if
%     the pairs have one-to-one correspondance in the penalty matrix.
%   However, the code does not let particles disappear and reappear if a
%   frame is missing, as track.m can.
%
% Examples of the code's use can be found in [2,3].
% For tracking large displacements in solids, we have found it useful to
% combine this code with a step to remove obviously wrong tracks, with the
% function remove_outliers_RAFT [2].
%
% [1] Boltyanskiy et al., Soft Matter 2017 (https://doi.org/10.1039/C6SM02011A)
% [2] Kim et al. (accepted in PRX) 2021 (https://arxiv.org/abs/2103.04975)
% [3] Testa et al. 2021 (https://www.biorxiv.org/content/10.1101/2021.05.16.444336v1.abstract)
%
