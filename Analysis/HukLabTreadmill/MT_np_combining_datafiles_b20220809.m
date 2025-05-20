% Combining datafiles

file1='/home/declan/data/BrieMTsession/20220809/ephys/2022-08-09_17-31-35/continuous/Neuropix-PXI-101.0/continuous-001.dat';
file2='/home/declan/data/BrieMTsession/20220809/ephys/2022-08-09_12-44-58/experiment1/continuous/Neuropix-PXI-101.0/continuous-002.dat';

fout='/home/declan/data/BrieMTsession/20220809/ephys/Combined/continuousfull.dat';

fid_write = fopen(fout, 'w');

fid_read1 = fopen(file1);
A = fread(fid_read1, '*int16');
fwrite(fid_write, A, 'int16')
fclose(fid_read1);

fid_read2 = fopen(file2);
A = fread(fid_read2, '*int16');
fwrite(fid_write, A, 'int16')
fclose(fid_read2);

fclose(fid_write)
%%
% fid_write = fopen(my_binary_file_path, 'w');
% for j = 1:length(flist)
%     fid_read = fopen(flist{j});
%     A = fread(fid_read, '*int16');
%     fwrite(fid_write, A, 'int16')
%     fclose(fid_read)
% end
% fclose(fid_write)