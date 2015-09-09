function [spread] = compute_feature_spread(features)
% Computes the spread of a feature vector by calculating the standard 
% deviation of a feature for a particular sequence of an individual. 
% An average is then taken across many different sequences using different 
% individuals.
%
% INPUTS
%  features       - A {1 x num_sequences} cell array. Each cell corresponds 
%                   to a sequence of feature vectors for a particular 
%                   individual in consecutive video frames. It contains a 
%                   [sequence_length x ftr_dimensions] array. 
%
% OUTPUTS
%  spread - an estimate of the spread of a feature
%
% Created by David Hall Apr-2013 (dhall@caltech.edu)

spread = cellfun(@std, features, 'UniformOutput', 0);
spread = cat(1, spread{:});
spread = single(median(spread));
spread = spread';


end