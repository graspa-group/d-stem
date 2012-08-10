myaddress = 'finazzif@gmail.com';
mypassword = 'ncc1701d';

setpref('Internet','E_mail',myaddress');
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',myaddress);
setpref('Internet','SMTP_Password',mypassword);

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class','javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

%sendmail(myaddress,'STEM4 - Kriging result','Kriging result','st_krig_result.mat');
sendmail(myaddress,'STEM4 - Model','Model','st_model.mat');