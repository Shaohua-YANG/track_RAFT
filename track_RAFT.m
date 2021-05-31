function trks = track_RAFT(xyzt,maxdisp,varargin)

% % This is a piece of code that unscrambles n-dimensional trajectories
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
% I hope you find this useful. If you use this code, please do cite the
% doi reference to the github page
% It helps us to keep doing science.
% Any feedback/bug reports etc would be very welcome. I will post up to
% date versions on the arxiv and on the Matlab file exchange
% 
%
% e.g. trks=track_RAFT(xyzt,2,'dim',3,'y',1)
%

% Inputs: xyzt: an array with the last column being the time (must be
% discrete... i.e. 1 2 3 4, not real times (i.e. 0.1 0.2 0.3). The first
% dim columns should be the particle points.
% Note - there must be at least n_consider + 1 particles at each time
% point, or the code will creash.

% maxdisp: an estimate of the maximum displacement that a particle will
% possible take between frames

% OPTIONS
% 'min_track_length': the minimum number of timepoints that individual
% particles should contain. Must be 2 or greater. Default is 2
%
% 'dim': the dimensions of the data to be unscrambled. If there are more than
% dim+1 columns, the tracker assumes that the first dim columns are the
% locations, and the last is the frame number. It will return the contents
% of the middle columns in the final result. Default is one less than the
% number of columns
%
% 'n_consider': the number of neighbouring points to consider when tracking
% individual points from frame to frame (the tracker looks for clusters of
% particles moving together, which helps remove drift). Default is 10
%
% 'n_use': the number of particles from n_consider to keep when calculating
% the penalty matrix. Essentially the tracker tries to match n_consider
% points around a particle at time 1,p1 ,to n_consider points around a particle
% at time 2, p2, to see if those points move together relative to the
% displacement from p1 to p2. If n_use<n_consider, the penalty function for
% p1 to p2 takes the smallest n_use out of n_consider displacements of the particle cloud,
% and adds their magnitude together. The idea is that if you know you have
% lost particles due to poor imaging, or particles going in and out of
% focus, you can reduce n_use relative to n_consider, and it will account
% for that. In general, the more you keep track of your particles in all
% your frames, the higher n_use should be relative to n_consider.
% n_use must be <= n_consider. If n_use>n_consider, it will be
% automatically set as n_consider. Default is 8.
%
% 'x': Tell the program that you want tracked displacements to be < this
% in the z direction (1st column). Default is maxdisp
%
% 'y': Tell the program that you want tracked displacements to be < this
% in the y direction (2nd column). Default is maxdisp
%
% 'z': Tell the program that you want tracked displacements to be < this
% in the z direction (3rd column). Default is maxdisp

% Outputs: trks: unscrambled data, where last column is particle number,
% previous columns are the same as in xyzt...

% Created by Rob Style in lockdown, completed 30/05/2021
%
% Notes by Rob. This runs pretty well for large data sets, though if
% needed, step 2 can be optimised to avoid id_tot containing lots of zeros.
% In particular, I would be tempted to remove rows from id_tot when a track
% terminates, and immediately add them to the trks structure.

inputExist = find(cellfun(@(x) strcmpi(x, 'min_track_length') , varargin));
if inputExist
    min_track_length = varargin{inputExist+1};
else
    min_track_length=2;
end
inputExist = find(cellfun(@(x) strcmpi(x, 'dim') , varargin));
if inputExist
    dim = varargin{inputExist+1};
else
    dim = size(xyzt,2)-1;
end
inputExist = find(cellfun(@(x) strcmpi(x, 'n_consider') , varargin));
if inputExist
    n_consider = varargin{inputExist+1};
else
    n_consider = 10;
end
inputExist = find(cellfun(@(x) strcmpi(x, 'n_use') , varargin));
if inputExist
    n_use = min(varargin{inputExist+1},n_consider);
else
    n_use = 8;
end
inputExist = find(cellfun(@(x) strcmpi(x, 'x') , varargin));
if inputExist
    maxdisp_x = varargin{inputExist+1};
else
    maxdisp_x=maxdisp;
end
inputExist = find(cellfun(@(x) strcmpi(x, 'y') , varargin));
if inputExist
    maxdisp_y = varargin{inputExist+1};
else
    maxdisp_y=maxdisp;
end
inputExist = find(cellfun(@(x) strcmpi(x, 'z') , varargin));
if inputExist
    maxdisp_z = varargin{inputExist+1};
else
    maxdisp_z=maxdisp;
end


% Step one, find the the number of time points that are 1 apart. If there
% are none of these (e.g. the last column is [1;3;5]. Then quit and tell
% the user to sort. Also have a few other error catchers

sz=size(xyzt);
if sum(xyzt(2:end,sz(2))-xyzt(1:end-1,sz(2)))==0
    disp('All positions are at the same time!')
    return
elseif min(xyzt(2:end,sz(2))-xyzt(1:end-1,sz(2)))>1
    disp('All time points more than 1 apart!')
    return
elseif sum(xyzt(:,end)-round(xyzt(:,end)))~=0
    disp('Time points must all be integers')
    return
elseif or(dim~=round(dim),dim<1)
    disp('dim must be an integer > 1')
    return
elseif or(min_track_length~=round(min_track_length),min_track_length<2)
    disp('min_track_length must be an integer > 2')
    return
end

% This is a list of the indices of rows that have times that are 1 time
% point apart.
trackable_time_inds=find(xyzt(2:end,sz(2))-xyzt(1:end-1,sz(2))==1);
trackable_time_inds=[trackable_time_inds;trackable_time_inds(end)+1];
trackable_times=xyzt(trackable_time_inds,sz(2));

% Now we need to work through one timepoint at a time, and track between
% the times

f = waitbar(0,'% of the way through time steps');

for k=1:length(trackable_times)
    % Get the information for the right time points
    pts1=xyzt(xyzt(:,sz(2))==trackable_times(k),1:dim);
    pts2=xyzt(xyzt(:,sz(2))==trackable_times(k)+1,1:dim);
    % Get the indices of the matched points
    res(k).idx=flow_tracker(pts1,pts2,maxdisp,n_consider,n_use,'x',maxdisp_x,'y',maxdisp_y,'z',maxdisp_z);
    waitbar(k/(length(trackable_times)-1),f,'% of the way through time steps (step 1/2)');       
end

% Now need to combine all the indices into a big mega index with all the
% particle positions..

t=unique([trackable_times;trackable_times+1]); % These are the times that we have trackable data at
id_tot=sortrows(res(1).idx);
if length(res)>1
    for k=2:length(res)
        if trackable_times(k-1)+1==trackable_times(k)
            no_rows=size(id_tot,1);
            no_cols=size(id_tot,2);
            id_tot=[id_tot,zeros(no_rows,1)];
            % Add the indices that match up to the previous time point
            for i=1:no_rows
                if sum(res(k).idx(:,1)==id_tot(i,no_cols))==1
                    id_tot(i,no_cols+1)=res(k).idx(res(k).idx(:,1)==id_tot(i,no_cols),2);
                    res(k).idx(res(k).idx(:,1)==id_tot(i,no_cols),:)=[];
                elseif sum(id_tot(i,:)~=0)<min_track_length
                    id_tot(i,:)=[]; % delete rows that have too short tracks to save space
                end
            end
        % Now add remaining indices below the index list
        sz_res=size(res(k).idx);
        id_tot=[id_tot ; zeros(sz_res(1),no_cols-1) res(k).idx];
        else
            % If there is a gap in the data so there are no particles at
            % the previous timestep, don't do the matching step..
            no_rows=size(id_tot,1);
            no_cols=size(id_tot,2);
            id_tot=[id_tot,zeros(no_rows,2)];
            sz_res=size(res(k).idx);
            id_tot=[id_tot ; zeros(sz_res(1),no_cols) res(k).idx];
        end
        waitbar(k/(length(res)+1),f,'% of the way through aligning indices (step 2/2)'); 
    end
end

% Throw away too short tracks
[rowIdcs, ~] = find(id_tot~=0);
[counts, ~] = hist(rowIdcs,1:size(id_tot,1)); % This has the number of nonzero elements in each row
id_tot(counts<min_track_length,:)=[]; % Throw away tracks that are shorter than min_track_length

trks=[];
% Convert into the usual tracking structure...
for i=1:length(t)
    % Get the information for the right time points
    pts1=xyzt(xyzt(:,end)==t(i),:);
    [rowIdcs, ~] = find(id_tot(:,i)~=0);
    trks=[trks; pts1(id_tot(id_tot(:,i)~=0,i),:) rowIdcs];
end
% Finally, sort by the particle tracking number...
trks=sortrows(trks,[sz(2)+1,sz(2)]);
close(f)
end



function idx=flow_tracker(pts1,pts2,maxdisp,n_consider,n_use,varargin)
% Rob's flow tracker.. idea is to track particles that move in clusters

% Inputs 
% pts1 List of particles in first frame, [X,Y], N x column vector 
% pts2 List of particles in second frame, [X,Y], N x column vector 
% n_use Number of nearest neighbours to use
% n_consider Number of nearest neighbours to consider
% maxdisp Maximum distance particles can move
% varargin If you want to put limits on how far a particle can move in x,y
% or z, just add 'x',2 or 'y',10 etc

% Outputs
% idx is a two column vector, with the first column being the number of the
% particle in pts1, and the second being the number that it corresponds to
% in pts2


% First set up the bounds on the displacement. See varargin above
k=0;
maxdisp_x=maxdisp;
maxdisp_y=maxdisp;
maxdisp_z=maxdisp;
if ~isempty(find(strcmp(varargin,'x')))
        ind=find(strcmp(varargin,'x'));
        maxdisp_x=min(varargin{ind(1)+1},maxdisp);
end
if ~isempty(find(strcmp(varargin,'y')))
        ind=find(strcmp(varargin,'y'));
        maxdisp_y=min(varargin{ind(1)+1},maxdisp);
end
if ~isempty(find(strcmp(varargin,'z')))
        ind=find(strcmp(varargin,'z'));
        maxdisp_z=min(varargin{ind(1)+1},maxdisp);
end

% For each point in pts1, find the nearest n_consider points in pts1, and do the same
% for each point in pts2.

near_neighb_inds_pts1=zeros(length(pts1),n_consider+1);
near_neighb_inds_pts2=zeros(length(pts2),n_consider+1);

% First get a list of indices of the nearest n_consider neighbours of each
% point in pts1
for i=1:size(pts1,1)
    [~,I]=mink(sum((pts1-pts1(i,:)).^2,2),n_consider+1);
    near_neighb_inds_pts1(i,:)=I';
end

% Next get a list of indices of the nearest n_consider neighbours of each
% point in pts1
for i=1:size(pts2,1)
    [~,I]=mink(sum((pts2-pts2(i,:)).^2,2),n_consider+1);
    near_neighb_inds_pts2(i,:)=I';
end

% Remove the first column, because these are the points themselves - not
% nearest neighbours
near_neighb_inds_pts1(:,1)=[];
near_neighb_inds_pts2(:,1)=[];

siz=size(pts1);
N=siz(2);

nn=zeros(length(pts1(:,1)),3);

for i=1:length(pts1(:,1))
    % Now for each point in pts1, find the nearest points in pts2 that
    % satisfy any constraints you put in on x, y and z displacements
    if N==1
       inds_near=sum((pts2-pts1(i,:)).^2,2)<maxdisp^2 & (pts2(:,1)-pts1(i,1)).^2<maxdisp_x^2;
    elseif N==2
       inds_near=sum((pts2-pts1(i,:)).^2,2)<maxdisp^2 & (pts2(:,1)-pts1(i,1)).^2<maxdisp_x^2 & (pts2(:,2)-pts1(i,2)).^2<maxdisp_y^2; 
    else
       inds_near=sum((pts2-pts1(i,:)).^2,2)<maxdisp^2 & (pts2(:,1)-pts1(i,1)).^2<maxdisp_x^2 & (pts2(:,2)-pts1(i,2)).^2<maxdisp_y^2 & (pts2(:,3)-pts1(i,3)).^2<maxdisp_z^2; 
    end
    % for each particle in pts1, inds_near gives a list of the indices of
    % particles in pts2 that could be matches
    inds_near=find(inds_near);
    
    if ~isempty(inds_near)
        for j=1:length(inds_near)
            % Work out the relative positions of the nearest neighbours to
            % the point in pts1
            ri=pts1(near_neighb_inds_pts1(i,:),:)-pts1(i,:);
            % do the same for pts2
            rj=pts2(near_neighb_inds_pts2(inds_near(j),:),:)-pts2(inds_near(j),:);
            % Calculate the squared distance matrix for each of the
            % relative particle points ri, rj
            dij = bsxfun(@plus,sum(ri.*ri,2),sum(rj.*rj,2)') - 2*ri*rj';
            
            % The cost in the cost matrix is the sum of the distances
            % between n_use points. Note this cheats slightly, as it
            % doesn't use the minimum in both column and row, but I've
            % found this works well in practice.
            pm(j)=sum(sqrt(mink(min(dij,[],2),n_use)));
        end
        % Find the minimum value of the penalty function for all the
        % potential particle pairs, and the index in pts2 of the
        % corresponding particle
        [penalty,ind_nn2]=min(pm);
        % Update a matrix that collects the index of the pts1 particle, the
        % index of its best match particle in pts2 and the penalty
        nn(i,1:3)=[i,inds_near(ind_nn2),penalty];
        clear pm
    else
        nn(i,1:3)=[NaN NaN NaN];
    end
end
% remove the NaN rows from nn
nn(isnan(nn(:,1)),:)=[];
% sort the rows to collect all the rows with the same pts2 particle
% together, and order so the one with the lowest penalty comes first
nn=sortrows(nn,[2,3]);
% throw away the repeat points in column 2 where the penalty isn't the
% lowest value
[~,inds]=unique(nn(:,2));
nn=nn(inds,:);
% sort the rows so that the first column is in order. This step can
% probably be deleted
nn=sortrows(nn,1);

% Spit out the indices of particle pairs from pts1 and pts2
idx=nn(:,1:2);

end
