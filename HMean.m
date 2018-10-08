function HMout = HMean(A,i,j,in,jn)
% A = Array to be averaged
% i,j  = location one
%in,ij = location two
%
%~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Developed by: 
%Ricardo Medina (ricarm3@uci.edu)        University of California, Irvine, Irvine, CA
%Russ Detwiler  (Detwiler@uci.edu)        University of California, Irvine, Irvine, CA
%
% Citation: Medina et al. (201#) FIXME. BEFORE FINAL SUBMISSION 
%
%~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

HMout = 2./( 1./A(i,j) + 1./A(in,jn) );

end

