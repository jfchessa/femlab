function svm = mises_val( s )

% function svm = mises_val( s )
%
% Computes the Mises effective value from a voigt vector
% 

sp=principal_val(s);

svm=abs( sqrt( 0.5*( (sp(2)-sp(1))^2 + (sp(3)-sp(2))^2 + (sp(1)-sp(3))^2 ) ));

