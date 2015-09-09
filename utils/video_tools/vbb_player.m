function vbb_player(vdir, anndir )
% Simple GUI to play annotated videos
%
% Uses the vid_player to display a video. See help for the
% vid_player for navigation options (playing forward/backward at variable
% speeds). Note that loading a new video via the file menu will not work
% (as the associated annotation is not loaded). To
% actually alter the annotations, use the vbb_labeler.
%
% USAGE
%  vbb_player(vidNm, annNm )
%
% INPUTS
%  vidNm      - video file name
%  annNm      - annotation file name
%
% See also vid_player, vbb, vbb_labeler
%
% Copyright Dec-2013 David Hall.  [dhall-at-caltech.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Simplified BSD License

if(nargin<1 || isempty(vdir)), vdir=[]; end
if(nargin<2 || isempty(anndir)), 
    [d,fn]=fileparts(vdir);
    if(isempty(d)), d='.'; end;
    d=[d '/'];
    anndir=[d 'Annotations/' fn '.vbb']; 
end

A = vbb('vbbLoad', anndir );
vbb( 'vbbPlayer', A, vdir );

end
