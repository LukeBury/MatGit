function sendEmailFromMatlab(receiver, subject, message)
%%% Description
%       Send an email to signify the end of a script
%       
% ------------------------------------------------------------------------
%%% Inputs
%       receiver - [string] Destination email address
%       subject  - [string] Subject of email
%       message  - [string] Body of email
% ------------------------------------------------------------------------
% Created: 11/14/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% -------------------------------------------------
%%% Setup
% -------------------------------------------------
senderEmail = 'SpacewhalesWorldwide@gmail.com'; 
password = 'SpaceWhales!1';  
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','E_mail',senderEmail);
setpref('Internet','SMTP_Username',senderEmail);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

% -------------------------------------------------
%%% Send email
% -------------------------------------------------
sendmail(receiver, subject, message);

end % function