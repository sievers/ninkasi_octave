function[value]=getdata_double_channel_piece(ff,chan,varargin)

first_frame=get_keyval_default('first_frame',0,varargin{:});
first_sample=get_keyval_default('first_sample',0,varargin{:});
num_frames=get_keyval_default('num_frames',0,varargin{:});
num_samples=get_keyval_default('num_samples',0,varargin{:});
%disp([first_frame first_sample num_frames num_samples])
value=getdata_double_channel_piece_c(ff,chan,first_frame,first_sample,num_frames,num_samples);
