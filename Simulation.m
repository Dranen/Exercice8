function Simulation(name, n, alpha, Rnucleus, V0, x0, sigmanorm, dt, tfinal, ndx, question, mode, nbsim, ndx_max, dt_max, alpha_max)

workingfolder = './';
binfilename = 'Exercice8';
input_file = [name '_inp.dat'];
output_file = [name '_out.dat'];

%create the input data file
fid = fopen( [ workingfolder, input_file], 'wt' ); %create or overwrite (empty file, text mode)

%fill the file
fprintf( fid, name );
fprintf( fid, '\n');
fprintf( fid, num2str(n));
fprintf( fid, '\n');
fprintf( fid, num2str(alpha));
fprintf( fid, '\n');
fprintf( fid, num2str(Rnucleus));
fprintf( fid, '\n');
fprintf( fid, num2str(V0));
fprintf( fid, '\n');
fprintf( fid, num2str(x0));
fprintf( fid, '\n');
fprintf( fid, num2str(sigmanorm));
fprintf( fid, '\n');
fprintf( fid, num2str(dt));
fprintf( fid, '\n');
fprintf( fid, num2str(tfinal));
fprintf( fid, '\n');
fprintf( fid, num2str(ndx));
fprintf( fid, '\n');
fprintf( fid, num2str(question));
fprintf( fid, '\n');
fprintf( fid, num2str(mode));
fprintf( fid, '\n');

if(mode == 2)
   fprintf( fid, num2str(nbsim));
   fprintf( fid, '\n'); 
   fprintf( fid, num2str(ndx_max));
   fprintf( fid, '\n'); 
end
if(mode == 3)
   fprintf( fid, num2str(nbsim));
   fprintf( fid, '\n'); 
   fprintf( fid, num2str(dt_max));
   fprintf( fid, '\n'); 
end
if(mode == 4)
   fprintf( fid, num2str(nbsim));
   fprintf( fid, '\n'); 
   fprintf( fid, num2str(alpha_max));
   fprintf( fid, '\n'); 
end

fclose( fid );

eval( [ '!', workingfolder, binfilename, ' < ', input_file, ' > ', output_file] );

fprintf('\n');
delete(input_file);

end

