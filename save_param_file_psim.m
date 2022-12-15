function save_param_file_psim(S_vars,filename)

%Salva dados txt
format long
mycell = {  'Vin' S_vars.Vin;
            'fo' S_vars.fo;
            'ma' S_vars.ma;
            
            'Vopk' S_vars.Vopk;
            'RoLin1' S_vars.RoLin1;
            'RoLin2' S_vars.RoLin2;
            'CLoad' S_vars.CLoad;
            'RoNL1' S_vars.RoNL1;
            'RoNL2' S_vars.RoNL2;
            
            'fs' S_vars.fs;
            'fa' S_vars.fa;
            
            'Lf' S_vars.Lf;
            'Cf' S_vars.Cf;
            
            'rL' S_vars.rL;
            'rC' S_vars.rC;
            
            'Rdson' S_vars.Rdson;
            'Vsdf' S_vars.Vsdf;
            'Cin' S_vars.Cin;
            
            'V3fin' S_vars.V3fin;
            'enRect' S_vars.enRect;
            'dVin' S_vars.dVin;
            
            't_dvin_up' S_vars.t_dvin_up;
            't_dvin_down' S_vars.t_dvin_down;
            
            't_step' S_vars.t_step;
            'sim_time' S_vars.sim_time;
            
            };
[nrows,ncols]= size(mycell);
%filename = 'Buck_psim_V01.txt';
fid = fopen(filename, 'w');
for row=1:nrows
    fprintf(fid, '%s %.12f\r\n', mycell{row,:});
end
fclose(fid);
% type(filename)
fprintf('Parameter file writen correctly\n')