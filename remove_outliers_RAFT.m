function tks2=remove_outliers_RAFT(tks,varargin)
%
% e.g. tks2=remove_outliers_RAFT(tks,'n',20,'abs_diff',5);
%
% code to identify outliers by comparing the movement of points to that of
% nearby points. If these significantly differ - either by 3 standard
% deviations, or by a fixed amount - then the track is split at the point
% where the deviation happens.
% If you find this useful please consider citing DOI: 10.5281/zenodo.4884065
% INPUT
% tks is the output of track_RAFT
% OPTIONS:
% 'min_track_length': the minimum number of timepoints that individual
% particles should contain. Split tracks that are shorter than this will be
% thrown away. Must be 2 or greater. Default is 2
%
% 'dim' is the number of dimensions the data has. Default is the number of
% columns in the tracking structure - 2
%
% 'n' is the number of neighbours to compare against. Default is 20
%
% 'abs_diff' is the maximum difference that a particles displacement can
% have from the average displacement of its nearest neighbours in any
% cartesian direction. Set to 0 if you don't wish to use this option.
% Default is 0
%
% Warning! If using this on data where you are measuring something like
% diffusivity, make sure it is working like it is meant to (i.e. it is not
% breaking up real tracks), otherwise it will bias your diffusivity data,
% by throwing away bigger jumps.
%
% Outputs
% trks2: the trks data with the particle tracks split at points where bad
% tracking is detected
%
% Created by Rob Style in lockdown, completed 31/05/2021

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
    dim = size(tks,2)-2;
end

inputExist = find(cellfun(@(x) strcmpi(x, 'n') , varargin));
if inputExist
    n_consider = varargin{inputExist+1};
else
    n_consider=20;
end
inputExist = find(cellfun(@(x) strcmpi(x, 'abs_diff') , varargin));
if inputExist
    abs_diff = varargin{inputExist+1};
else
    abs_diff=0; % if zero, it is not applied
end

ts=unique(tks(:,end-1));

f = waitbar(0,'% of the way through time steps');

for i=1:length(ts)
    if ~isempty(tks(tks(:,end-1)==ts(i)+1,:))
        pno_diff=tks(:,end) - circshift(tks(:,end),-1);
        t_diff=tks(:,end-1) - circshift(tks(:,end-1),-1);
        inds_pair_ts1=pno_diff==0 & t_diff==-1 & tks(:,end-1)==ts(i);
        pno_diff=tks(:,end) - circshift(tks(:,end),1);
        t_diff=tks(:,end-1) - circshift(tks(:,end-1),1);
        inds_pair_ts2=pno_diff==0 & t_diff==1 & tks(:,end-1)==ts(i)+1;
        inds_bad=find_outliers_RAFT(tks(inds_pair_ts1,1:dim),tks(inds_pair_ts2,1:dim),'n',n_consider,'abs_diff',abs_diff);
        inds=inds_pair_ts1;
        inds(inds_pair_ts1)=inds_bad; % This has a list of rows which correspond to bad pairs
        % Now take bad tracks and split the tracks there and allocate the new split track a new number
        bad_pnos=tks(inds,end);
        for j=1:length(bad_pnos)
            if sum(tks(:,end)==bad_pnos(j) & tks(:,end-1)>ts(i),'all')>=min_track_length % if there are less than min_track_length particles left in the track after the split, delete. Otherwise assign the rest of the track a new number
                tks(tks(:,end)==bad_pnos(j) & tks(:,end-1)>ts(i),end)=max(tks(:,end))+1;
            else
                tks(tks(:,end)==bad_pnos(j) & tks(:,end-1)>ts(i),:)=[];
            end
            if sum(tks(:,end)==bad_pnos(j) & tks(:,end-1)<=ts(i),'all')<min_track_length % if there are less than min_track_length particles left in the track before the split, delete
                tks(tks(:,end)==bad_pnos(j) & tks(:,end-1)<=ts(i),:)=[];
            end
        end
        %tks(ismember(tks(:,end),tks(inds,end)),:)=[];     
        waitbar(i/(length(ts)),f,'% of the way through time steps');
    end
end
close(f)
tks2=tks;
end

function inds_bad=find_outliers_RAFT(pks1,pks2,varargin)

% code to identify outliers. Idea is that for every point, find the
% n nearest points, and work out what the average displacement vector is for
% them. Then compare and see if they're within 3 standard deviations of
% each other
% INPUTS:
% pks1,pks2 are two arrays that are the same size, and have matched data
% sets. I.e. pks1 contains the positions of points at one time, and pks2
% contains points at the second time, where points in the same row
% correspond to each other
% OPTIONS:
% 'n' is the number of neighbours to compare against. Default is 20
% 'abs_diff' is the maximum difference that a particles displacement can
% have from the average displacement of its nearest neighbours in any
% cartesian direction
inputExist = find(cellfun(@(x) strcmpi(x, 'n') , varargin));
if inputExist
    n_consider = varargin{inputExist+1};
else
    n_consider=20;
end
inputExist = find(cellfun(@(x) strcmpi(x, 'abs_diff') , varargin));
if inputExist
    abs_diff = varargin{inputExist+1};
else
    abs_diff=0;
end

% Start by working out the average distance between particles
near_neighb_inds=zeros(length(pks1),n_consider+1);
for i=1:size(pks1,1)
    [~,I]=mink(sum((pks1-pks1(i,:)).^2,2),n_consider+1);
    near_neighb_inds(i,:)=I';
end

% Get rid of first column, as this is just the particle itself.
near_neighb_inds(:,1)=[];
%[~,ids]=mink(D11,n+1,2);
%ids=near_neighb_inds;
%ids(:,1)=[];
disp=pks2-pks1; % displacements

disp_mean=zeros(size(pks1,1),size(pks1,2));
disp_std=zeros(size(pks1,1),size(pks1,2));

for i=1:size(pks1,1)
    disp_temp=pks2(near_neighb_inds(i,:),:)-pks1(near_neighb_inds(i,:),:);
    disp_mean(i,:)=mean(disp_temp,1);
    disp_std(i,:)=std(disp_temp,1);
end

if abs_diff==0
    bads=abs(disp_mean-disp)>3*disp_std;
else
    bads=[abs(disp_mean-disp)>3*disp_std] | [abs(disp_mean-disp)>abs_diff];
end

inds_bad=any(bads,2); %bads(:,1) | bads(:,2) | bads(:,3);

end
