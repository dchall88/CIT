function results_player(vid_name, detected_bbs)
% Simple GUI to play videos with detected bounding boxes overlayed
%
% Uses the vid_player to display a video.
%
% USAGE
%  resultsPlayer(vid_name, detected_bbs)
%
% INPUTS
%  vid_name          - video file name
%  detected_bbs      - detected bounding boxes
%
% Copyright Dec-2013 David Hall.  [dhall-at-caltech.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Simplified BSD License
dispFunc = @(f) drawToFrame(detected_bbs{f+1});
vid_player(vid_name, dispFunc);
end

function hs = drawToFrame(bbs)
cols = uniqueColors(3,8);
hs = bbApply('draw', bbs, cols, 2, '-',[], mod(bbs(:,6) - 1, size(cols,1)) + 1 );
end