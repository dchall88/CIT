function [individual_classifier] = ...
    create_individual_clf(feature_vector, category_classifier, spread, beta)
% Creates an individual classifier from the feature vector F, the category
% classifier clf, the average spread, and the parameter beta.
%
% INPUTS
%  feature_vector       - the feature vector F (or array of row vectors)
%                         of an instance of the individual to model.
%
%  category_classifier  - the pre-trained, boosted category classifer 
%                         with stumps as weak classifiers.
%                         Refer acfTrain and adaBoost Train in 
%                         Piotr Dollar's toolbox
%                         http://vision.ucsd.edu/~pdollar/toolbox/doc/
%
%  spread               - the spread of each feature for an individual
%                         sequence averaged across many 
%                         sequences. Learnt offline from a training set.
%                         Refer to compute_feature_spread 
%
%  beta                 - tuning parameter. thresholds = beta * spread
%
% OUTPUTS
%  individual_classifier - the individual classifier
%
% Created by David Hall Apr-2013 (dhall@caltech.edu)

%% Only modify the thresholds of the boosted category detector
thresholds = spread .* beta;
individual_classifier =...
    set_thresholds(category_classifier, feature_vector, thresholds);

end

function classifier =...
    set_thresholds(classifier, feature_vector, thresholds)

N = size(classifier.fids, 2); % Number of weak classifiers
depth = classifier.treeDepth; % Depth of weak classifiers
if(depth ~=1)
    error('The category classifier must use stumps for weak classifiers');
end
W = max(classifier.hs);    % Weak classifier weights 

% Get zero-indexed feature indices from the internal nodes of the weak 
% classifiers and convert to one-indexed indices.
ftr_indices = classifier.fids(1, :) + 1;

% The individual classifiers use 'intervals' as weak classifiers. These 
% intervals can be implemented as decision trees. These new decision trees
% have the following parameters.
num_nodes = 5; 
num_leaves = 3; 
child = [2; 0; 4; 0; 0]; % Child nodes of each internal node in the tree
node_index = logical(child); % Indices of internal nodes (non-leaf nodes)
depths = [0; 1; 1; 2; 2]; %Depths of each node in the tree

% Modify classifer with new parameters
classifier.treeDepth = 0;
[classifier.hs, classifier.thrs]=deal(zeros(num_nodes, N, 'single'));
classifier.fids = zeros(num_nodes, N, 'uint32');

% Set thresholds for the weak classifiers
ftr_indices = ftr_indices([1 1], :);
classifier.thrs(node_index,:) = feature_vector(ftr_indices)... 
                              + [-ones(1, N); ones(1, N)]... 
                              .* thresholds(ftr_indices);
classifier.fids(node_index,:) = ftr_indices-1;

% Set weights for the weak classfiers
hs = -ones(num_leaves,N); 
hs(2, :) = 1;
classifier.hs(~node_index, :) = repmat(W, num_leaves, 1) .* hs;

% Update child and depth fields
classifier.child = uint32(repmat(child, 1, N));
classifier.depth = uint32(repmat(depths, 1, N));

end

