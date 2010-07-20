function[varargout]=cexec(to_eval)
mpi_bcast_string(to_eval);

evalin('caller',to_eval);


%eval(['[varargout]={' to_eval '};']);

%eval(varargout=to_eval);