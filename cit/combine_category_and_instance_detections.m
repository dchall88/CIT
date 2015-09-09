function [Bboxes, Patches] = combine_category_and_instance_detections(bboxes, patches, minDetectionOverlap, minCategoryDetectorScore)

% minDetectionOverlap - Category detection and Individual detection must overlap by at least this amount
% minCategoryDetectorScore - Category detection minimum confidence score


%% Split Detections into Category and Instance (category detections have id = 1)
bboxes_cat = bboxes(bboxes(:, 6) == 1, :); 
nCategoryBoxes = size(bboxes_cat,1); 
patches1 = patches(1:nCategoryBoxes, :);

bboxes_instance = bboxes(bboxes(:, 6) ~= 1, :); 
nInstanceBoxes = size(bboxes_instance, 1);

%% Filter Detections
% Suppress multiple firings of individual detector (should only be at
% most 1 individual in the frame - choose one with maximum score
if nInstanceBoxes>0
    ind = bboxes_instance(:,6); 
	id = unique(ind); 
	ind = logical(ind);
    for k = 1:length(id)
        ids = find(bboxes_instance(:,6) == id(k));
        [~, J] = max(bboxes_instance(ids, 5)); 
		ind(ids(J)) = 0;
    end
    bboxes_instance(ind, :)=[];
end

%% Determine Assignment:
% 3 options - a) Unique category detection (corresponds to new individual)
%             b) Unique individual detection (corresponds to existing
%                     individual that was missed by category detector)
%             c) Category and individual detections overlap (combine result)

oa1 = bbGt('compOas', bboxes_cat(:,1:4), bboxes_instance(:,1:4));
oa = oa1 > minDetectionOverlap;

% a) Unique Category Detections
ind_cat_unq = all(oa == 0, 2);
bboxes_cat_unq = bboxes_cat(ind_cat_unq, :); 
nCategoryBoxes_UNQ = size(bboxes_cat_unq, 1);
patches1_unq = patches1(ind_cat_unq, :);
% Remove unique category detections that have low confidence score
if(nCategoryBoxes_UNQ > 0)
    ind=bboxes_cat_unq(:, 5) > minCategoryDetectorScore;
    bboxes_cat_unq = bboxes_cat_unq(ind, :);
    patches1_unq = patches1_unq(ind, :);
end
%Create overlapping category detections
bboxes_cat = bboxes_cat(~ind_cat_unq, :);
patches1 = patches1(~ind_cat_unq, :);
oa1 = oa1(~ind_cat_unq, :);

% b) Unique Individual Detections
ind_instance_unq = all(oa == 0, 1);
bboxes_instance_unq = bboxes_instance(ind_instance_unq, :);
%Create overlapping individual detections
bboxes_instance = bboxes_instance(~ind_instance_unq, :);
oa1 = oa1(:, ~ind_instance_unq);

% c) Overlapping Detections
costMat = 1./oa1;
assignMat = lapjv(costMat);
if size(costMat, 1) > size(costMat, 2),
    ind_instance_ol = 1:size(costMat, 2); 
	ind_cat_ol = assignMat;
    %Unassigned category detections get added to unique cat. dets
    ind = setdiff(1:size(costMat, 1), assignMat);
    if ~isempty(ind),
        b = bboxes_cat(ind, :); 
		p = patches1(ind, :);
        ind = b(:, 5) > minCategoryDetectorScore;
        bboxes_cat_unq = [bboxes_cat_unq; b(ind, :)];
        patches1_unq = [patches1_unq; p(ind, :)];
    end
    patches1 = patches1(ind_cat_ol, :);
else
    ind_cat_ol = 1:size(costMat, 1); 
	ind_instance_ol = assignMat;
    %Unassigned individual detections get added to unique instance. dets
    ind = setdiff(1:size(costMat, 2), assignMat);
    if ~isempty(ind) 
		bboxes_instance_unq = [bboxes_instance_unq; bboxes_instance(ind, :)]; 
    end
end

%Category bboxes and instance bboxes that overlap sufficiently get converted into
%a single bb and the scores are added together
bboxes_ol = bbApply('union', bboxes_instance(ind_instance_ol, 1:4), bboxes_cat(ind_cat_ol, 1:4));
bboxes_ol(:,5) = bboxes_instance(ind_instance_ol, 5) + bboxes_cat(ind_cat_ol, 5);
bboxes_ol(:,6) = bboxes_instance(ind_instance_ol, 6);
bboxes_ol(:, 7:11) = bboxes_cat(ind_cat_ol, 1:5);
bboxes_ol(:, 12:16) = bboxes_instance(ind_instance_ol, 1:5);

Bboxes.cat = bboxes_cat_unq; 
Patches.cat = patches1_unq;

Bboxes.both = bboxes_ol; 
Patches.both = patches1;

Bboxes.instance = bboxes_instance_unq;

end