function varargout = vid_writer_plugin( cmd, h, varargin )
% Plugin for vid_io to allow writing of video files.
%
% Do not call directly, use as plugin for vid_io instead.
% The following is a list of commands available (swp=vidWriterPlugin):
%  h=swp('open',h,fName,info) % Open a video file for writing (h ignored).
%  h=swp('close',h)           % Close video file (output h is -1).
%  swp('addframe',h,I,[ts])   % Writes video frame (and timestamp).
%  info = swp('getinfo',h)    % Return struct with info about video.
%
% The following params may be specified in struct 'info' upon opening:
%  FrameRate      - [30] frames per second
%  codec          - [MPEG-4] string representing codec, refer to VideoWriter
%
% USAGE
%  varargout = vidWriterPlugin( cmd, h, varargin )
%
% INPUTS
%  cmd        - string indicating operation to perform
%  h          - unique identifier for open video file
%  varargin   - additional options (vary according to cmd)
%
% OUTPUTS
%  varargout  - output (varies according to cmd)
%
% Copyright Dec-2013 David Hall.  [dhall-at-caltech.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Simplified BSD License

% persistent variables to keep track of all loaded .seq files
persistent h1 hs vids infos;
if(isempty(h1)), h1=int32(now); hs=int32([]); infos={}; end
nIn=nargin-2; in=varargin; o1=[]; cmd=lower(cmd);

% open video file
if(strcmp(cmd,'open'))
  chk(nIn,2); h=length(hs)+1; hs(h)=h1; varargout={h1}; h1=h1+1;
  fName=in{1};
  [infos{h},vids{h}]=openVideo(fName,in{2}); return;
end

% Get the handle for this instance
[v,h]=ismember(h,hs); if(~v), error('Invalid load plugin handle'); end
vid=vids{h}; info=infos{h}; 

% close video file
if(strcmp(cmd,'close'))
  chk(nIn,0); varargout={-1}; 
  close(vid); kp=[1:h-1 h+1:length(hs)];
  hs=hs(kp); vids=vids(kp); infos=infos(kp);
  return;
end

% perform appropriate operation
switch( cmd )
  case 'addframe',  chk(nIn,1,2); info=addFrame(vid,info,in{:});
  case 'getinfo',   chk(nIn,0); o1=info;
  otherwise,        error(['Unrecognized command: "' cmd '"']);
end
infos{h}=info; varargout={o1};

end

function chk(nIn,nMin,nMax)
if(nargin<3), nMax=nMin; end
if(nIn>0 && nMin==0 && nMax==0), error(['"' cmd '" takes no args.']); end
if(nIn<nMin||nIn>nMax), error(['Incorrect num args for "' cmd '".']); end
end

function [info, vid] = openVideo( fName, info )
% open video for writing, create space for header
if(~isfield(info,'codec')), info.codec='MPEG-4'; end
if(~isfield(info,'FrameRate')), info.FrameRate=30; end
vid=VideoWriter(fName,info.codec);
set(vid,'FrameRate',info.FrameRate);
open(vid);
end

function info = addFrame( vid, info, I)
% write frame
writeVideo(vid,I);
end
