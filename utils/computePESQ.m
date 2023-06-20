function pesqScore=computePESQ(s,vsr,par)

% COMPUTEPESQ A wrapper m-file function to run a binary PESQ function on
% the speech signals s and vsr.
%
%  Inputs: s                    Clean speech signal
%          vsr                  Reconstructed speeh signal
%          par                  A parameter strcuture containing:
%            .fs                Sampling frequency of the speech signals
%            .PESQpath          A path to the binary PESQ function
%
% Outputs: pesqScore            The PESQ score

v_writewav(s,par.fs,'tmp1.wav');
v_writewav(vsr,par.fs,'tmp2.wav');

pesq_path = par.PESQpath;
command = " +16000 tmp1.wav tmp2.wav | tail -1 | grep -oE '[^ ]+$'";

input = strcat(pesq_path,command);

[~,output] = system(input);

pesqScore = str2double(output);

system('rm tmp1.wav tmp2.wav _pesq_results.txt _pesq_itu_results.txt');
