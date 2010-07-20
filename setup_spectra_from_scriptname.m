function[value]=setup_spectra_from_scriptname(fname,bigset,myset,fac,opts)
if ~exist('opts')
  opts.asdfasdf='';
end

opts=set_struct_mem_default(opts,'tag','');
opts=set_struct_mem_default(opts,'kmax',6000);
opts=set_struct_mem_default(opts,'iters',[50 100 200 500 1000]);




[outdir,outname]=get_outname_from_script (fname,bigset,myset,fac);
outdir=[pwd '/' outdir];
outname=turn_outname_to_runSpectra (outname);



specdir=[outdir '/spectra_' num2str(opts.kmax)];
if ~isempty(opts.tag)
  specdir=[specdir '_' opts.tag];
end


system(['mkdir ' specdir]);
pbsname=[specdir '/runSpectra.pbs'];

fid=fopen(pbsname,'w');
if (fid=='-1')
  error(['unable to open pbs script ' pbsname 'for writing.']);
end



fprintf(fid,'#!/bin/bash\n');
fprintf(fid,['#PBS -N spectra_k' num2str(opts.kmax) '_' fname '\n']);
fprintf(fid,'#PBS -q largemem\n');
fprintf(fid,'#PBS -l nodes=1:ppn=16\n');
fprintf(fid,'#PBS -l walltime=48:00:00\n');
fprintf(fid,'#PBS -e /scratch/sievers/${PBS_JOBID}_stderr.txt\n');
fprintf(fid,'#PBS -o /scratch/sievers/${PBS_JOBID}_stdout.txt\n');
fprintf(fid,'#PBS -r n\n');
fprintf(fid,'\n\n');

fprintf(fid,['cd ' specdir '\n']);
fprintf(fid,['runSpectra runSpectra.dict\n']);
fclose(fid);


template_dir='/home/sievers/act/mapping2/spectra/templates/';
system(['cp ' template_dir '/specPatchesNWayMasterIterTemplate.dict ' specdir]);
system(['cp ' template_dir '/processSpectraNWayIterTemplate.dict ' specdir]);
system(['cp ' template_dir '/createPatchesNWayIterTemplate.dict ' specdir]);


infid=fopen([template_dir '/specPatchesNWayMasterIterTemplate.dict'],'r');
fid=fopen([specdir '/specPatchesNWayMasterIterTemplate.dict'],'w');
done=false
while ~done
  mystr=fgetl(infid);
  if ischar(mystr)
    [a,b]=strtok(mystr);
    if (strcmp(a,'lmaxForMaster'))
      kkk=num2str(opts.kmax);
      if sum(kkk=='.')==0
        kkk=[kkk '.'];
      end
      mystr=['lmaxForMaster = ' kkk];
    end
    fprintf(fid,'%s\n',mystr);
  else
    fclose(infid);
    fclose(fid);
    done=true;
  end
  
end


infid=fopen([template_dir '/runSpectra.dict'],'r');
fid=fopen([specdir '/runSpectra.dict'],'w');
done=false
while ~done
  mystr=fgetl(infid);
  if ischar(mystr)
    [a,b]=strtok(mystr);
    if (strcmp(a,'mapDir'))
      mystr=['mapDir = ''' outdir ''''];
    end
    if (strcmp(a,'mapNameTemplate'))
      mystr=['mapNameTemplate = ''' outname ''''];
    end
    if (strcmp(a,'iterations'))
      iters=opts.iters;
      mystr=['iterations = [' num2str(iters(1))];
      for j=2:length(opts.iters)
        mystr=[mystr ',' num2str(iters(j))];
      end
      mystr=[mystr ']'];
    end
    fprintf(fid,'%s\n',mystr);
  else
    fclose(infid);
    fclose(fid);
    done=true;
  end
end


system(['qsub ' pbsname]);    


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[myopts]=set_struct_mem_default(myopts,mystr,myval)
if ~isfield(myopts,mystr)
  eval(['myopts.' mystr ' = myval;']);
end
