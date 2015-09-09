function [features, thumbnails] =...
    extract_feature_tracks(video_name, gt_name, detector)
% Extracts sequence of feature vectors for a particular
% individual in consecutive video frames.
%
% INPUTS
%  video_name         - path of video from which to extract feature tracks.
%  gt_name            - path of annotations file (a .vbb file)
%  category_detector  - pre-trained category detector
%
% OUTPUTS
%  features        - A cell array of length num_tracks. Each cell contains the corresponding features of the image patches in thumbs.
%  thumbnails      - A cell array of the image patches extracted from the annotated video.
%
% Created by David Hall Apr-2013 (dhall@caltech.edu)

vid = vid_io(video_name, 'reader'); %load video
A = vbb('vbbLoad', gt_name); %load annotations
stats = vbb('getStats', A); %get annotation statistics

% Load detector parameters
shrink = detector.opts.pPyramid.pChns.shrink;
model_dims = detector.opts.modelDs;
model_dims_pad = detector.opts.modelDsPad;
channel_parameters = detector.opts.pPyramid.pChns;
smooth = detector.opts.pPyramid.smooth;

% Window size parameters. Used to extract image patch from image
pSmp.dims = model_dims([2, 1]);
cropDs = max(8 * shrink, model_dims_pad);
pad = cropDs ./ model_dims - 1;
pSmp.pad = pad([2 1]);

% For each individual track, extract features and thumbnails.
track_num = unique(stats.objnums);
num_tracks = length(track_num);

features = cell(1, num_tracks);
thumbnails = cell(1, num_tracks);

ticId=ticStatus('Begin...', [], [], 1);
for k=1:num_tracks
    current_track = stats.objnums == track_num(k);
    
    frames = stats.frame(current_track);
    bb = stats.pos(current_track, :);
    occluded = stats.occl(current_track);
    
    length_curr_track = sum(~occluded);
    Is = zeros(cropDs(1), cropDs(2), 3, length_curr_track, 'uint8');
    
    % Extract thumbnails
    itr = 1;
    for j=1:length(frames)  
        if(~occluded(j))
            vid.seek(frames(j));
            I = vid.getframe();
            I1 = sampleWins(I, bb(j,:), pSmp);
            Is(:,:,:,itr) = I1{:};
            itr = itr + 1;
        end
    end
    
    % Compute Features
    chns = chnsCompute1(Is, model_dims_pad, channel_parameters, smooth);

    features{k} = reshape(chns, size(chns,4), []);
    thumbnails{k} = Is;
    
    tocStatus(ticId, k./num_tracks);
end

end

function chns = chnsCompute1( Is, modelDsPad, pChns, smooth )
% From Piotr Dollar's Toolbox
% Compute single scale channels of dimensions modelDsPad.
ds=[size(Is,1) size(Is,2) size(Is,3)]; n=size(Is,4); shrink=pChns.shrink;
cr=rem(ds(1:2),shrink); s=floor(cr/2)+1; e=ceil(cr/2);
Is=Is(s(1):end-e(1),s(2):end-e(2),:,:); ds(1:2)=[size(Is,1) size(Is,2)];
if(any(ds(1:2)<modelDsPad)), error('Windows too small.'); end
Is=squeeze(mat2cell2(Is,[1 1 1 n])); chns=cell(1,n);
parfor i=1:n
    C=chnsCompute(Is{i},pChns); C=convTri(cat(3,C.data{:}),smooth);
    ds1=size(C); cr=ds1(1:2)-modelDsPad/shrink; s=floor(cr/2)+1;
    e=ceil(cr/2); C=C(s(1):end-e(1),s(2):end-e(2),:); chns{i} = C;
    assert(numel(C)==prod(modelDsPad/shrink)*size(C,3));
end; chns=cat(4,chns{:}); assert(isa(chns,'single'));
end

function Is = sampleWins( I, bbs, varargin )
% Sample pos or neg windows from an annotated image.
%
% An annotated image can contain both pos and neg examples of a given class
% (such as a pedestrian). This function allows for sampling of only pos
% windows without sampling any negs, or vice-versa. For example, this can
% be quite useful during bootstrapping, to sample high scoring false pos
% without actually sampling any windows that contain true pos.
%
% bbs should contain the candidate bounding boxes, and ibbs should contain
% the bounding boxes that are to be ignored. During sampling, only bbs that
% do not match any ibbs are kept (two bbs match if their area of overlap is
% above the given thr, see bbGt>compOas). Use gtBbs=bbLoad(...) to obtain a
% list of ground truth bbs containing the positive windows (with ignore
% flags set as desired). Let dtBbs contain detected bbs. Then:
%  to sample true-positives, use:   bbs=gtBbs and ibbs=[]
%  to sample false-negatives, use:  bbs=gtBbs and ibbs=dtBbs
%  to sample false-positives, use:  bbs=dtBbs and ibbs=gtBbs
% To sample regular negatives without bootstrapping generate bbs
% systematically or randomly (see for example bbApply>random).
%
% dims determines the dimension of the sampled bbs. If dims has two
% elements [w h], then the aspect ratio (ar) of each bb is set to ar=w/h
% using bbApply>squarify, and the extracted patches are resized to the
% target w and h. If dims has 1 element then ar=dims, but the bbs are not
% resized to a fixed size. If dims==[], the bbs are not altered.
%
% USAGE
%  Is = bbGt( 'sampleWins', I, pSmp )
%
% INPUTS
%  I        - input image from which to sample
%  pSmp     - parameters (struct or name/value pairs)
%   .n        - [inf] max number of bbs to sample
%   .bbs      - [REQ] candidate bbs from which to sample [x y w h ign]
%   .ibbs     - [] bbs that should not be sampled [x y w h ign]
%   .thr      - [.5] overlap threshold between bbs and ibbs
%   .dims     - [] target bb aspect ratio [ar] or dims [w h]
%   .squarify - [1] if squarify expand bb to ar else stretch patch to ar
%   .pad      - [0] frac extra padding for each patch (or [padx pady])
%   .padEl    - ['replicate'] how to pad at boundaries (see bbApply>crop)
%
% OUTPUTS
%  Is       - [nx1] cell of cropped image regions
%
% EXAMPLE
%
% See also bbGt, bbGt>bbLoad, bbApply>crop, bbApply>resize,
% bbApply>squarify bbApply>random, bbGt>compOas
%
% From Piotr Dollar's Toolbox

% get parameters
dfs={'n',inf, 'ibbs',[], 'thr',.5, 'dims',[], ...
    'squarify',1, 'pad',0, 'padEl','replicate' };
[n,ibbs,thr,dims,squarify,pad,padEl] = getPrmDflt(varargin,dfs,1);
% discard any candidate bbs that match the ignore bbs, sample to at most n
nd=size(bbs,2); if(nd==5), bbs=bbs(bbs(:,5)==0,:); end; m=size(bbs,1);
if(isempty(ibbs)), if(m>n), bbs=bbs(randSample(m,n),:); end; else
    if(m>n), bbs=bbs(randperm(m),:); end; K=false(1,m); i=1;
    keep=@(i) all(compOas(bbs(i,:),ibbs,ibbs(:,5))<thr);
    while(sum(K)<n && i<=m), K(i)=keep(i); i=i+1; end; bbs=bbs(K,:);
end
% standardize aspect ratios (by growing bbs), pad bbs and finally crop
if(numel(dims)==2), ar=dims(1)/dims(2); else ar=dims; dims=[]; end
if(numel(pad)==1), pad=[pad pad]; end; if(dims), dims=dims.*(1+pad); end
if(~isempty(ar) && squarify), bbs=bbApply('squarify',bbs,0,ar); end
if(any(pad~=0)), bbs=bbApply('resize',bbs,1+pad(2),1+pad(1)); end
Is=bbApply('crop',I,bbs,padEl,round(dims));

end