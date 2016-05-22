function Y_out = admittance_transfer(TR_matrix, Y_in)
% The 'admittance_transfer' function finds the "upstream" admittance,
% Y_out, given the "downstream" admittance, Y_in, and the connecting
% transfer matrix, TR_matrix.
%
% USAGE: Y_out = admittance_transfer(TR_matrix, Y_in);
%
partial = TR_matrix*[1; Y_in];
Y_out = partial(2)/partial(1);
