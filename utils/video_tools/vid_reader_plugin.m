function varargout = vid_reader_plugin( cmd, h, varargin )
% Plugin for vid_io to allow reading of video files.
%
% Do not call directly, use as plugin for vid_io instead.
% The following is a list of commands available (srp=vid_reader_plugin):
%  h = srp('open',h,fName)        % Open a video file for reading (h ignored).
%  h = srp('close',h);            % Close video file (output h is -1).
%  I = srp('getframe',h)          % Get current frame (returns [] if invalid).
%  I = srp('getframes',h,frame)   % Get specified frame/frames.
%  I = srp('getnext',h)           % Shortcut for 'next' followed by 'getframe'.
%  info = srp('getinfo',h)        % Return struct with info about video.
%  out = srp('next',h)            % Go to next frame (out=0 on fail).
%  out = srp('seek',h,frame)      % Go to specified frame (out=0 on fail).
%  out = srp('step',h,delta)      % Go to current frame+delta (out=0 on fail).
%
% USAGE
%  varargout = vid_reader_plugin( cmd, h, varargin )
%
% INPUTS
%  cmd        - string indicating operation to perform
%  h          - unique identifier for open video file
%  varargin   - additional options (vary according to cmd)
%
% OUTPUTS
%  varargout  - output (varies according to cmd)
%
% EXAMPLE
%
% See also vid_io, vid_writer_plugin
%
% Copyright Dec-2013 David Hall.  [dhall-at-caltech.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Simplified BSD License

% persistent variables to keep track of all loaded video files
persistent h1 hs cs vids infos
if(isempty(h1)), h1=int32(now); hs=int32([]); infos={}; end
nIn=nargin-2; in=varargin; o2=[]; cmd=lower(cmd);

% open video file
if(strcmp(cmd,'open'))
    chk(nIn,1,2); h=length(hs)+1; hs(h)=h1; varargout={h1}; h1=h1+1;
    cs(h)=0; fName=in{1};
    [infos{h},vids{h}]=open(fName); return;
end

% Get the handle for this instance
[v,h]=ismember(h,hs); if(~v), error('Invalid load plugin handle'); end
c=cs(h); vid=vids{h}; info=infos{h};

% close video file
if(strcmp(cmd,'close'))
    chk(nIn,0); varargout={-1};
    kp=[1:h-1 h+1:length(hs)];
    hs=hs(kp); cs=cs(kp); vids=vids(kp); infos=infos(kp);
    return;
end

% perform appropriate operation
switch( cmd )
    case 'getframe',      chk(nIn,0); o1=getFrame(c,vid,info);
    case 'getframes',     chk(nIn,1); o1=getFrame(in{1},vid,info);
    case 'getinfo',       chk(nIn,0); o1=info; o1.curFrame=c;
    case 'getnext',       chk(nIn,0); c=c+1; o1=getFrame(c,vid,info);
    case 'next',          chk(nIn,0); [o1,c]=valid(c+1,info);
    case 'seek',          chk(nIn,1); [o1,c]=valid(in{1},info);
    case 'step',          chk(nIn,1); [o1,c]=valid(c+in{1},info);
    otherwise,            error(['Unrecognized command: "' cmd '"']);
end
cs(h)=c; varargout={o1,o2};

end

function chk(nIn,nMin,nMax)
if(nargin<3), nMax=nMin; end
if(nIn>0 && nMin==0 && nMax==0), error(['"' cmd '" takes no args.']); end
if(nIn<nMin||nIn>nMax), error(['Incorrect num args for "' cmd '".']); end
end

function [info, vid] = open( fName )
% open video for reading, get header
if(exist(fName,'file')==0), error('video file not found: %s%s',fName); end
fName = GetFullPath(fName);
V = mmread(fName,1); 
info.NumberOfFrames = abs(V.nrFramesTotal);
info.height=V.height; 
info.width=V.width; 
info.fps=V.rate;
vid = fName;
end

function [v,frame] = valid( frame, info )
fst=min(frame); lst=max(frame);
nFrms=abs(info.NumberOfFrames);
v=(fst>=1 && fst<=nFrms) && (lst>=1 && lst<=nFrms);
end

function [I] = getFrame( frame, vid, info)
% get frame image (I) and timestamp (ts) at which frame was recorded
if(~valid(frame,info)), I=[]; return; end
V = mmread(vid, frame);
I = cat(4,V.frames(:).cdata);
end