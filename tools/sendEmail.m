function sendEmail(subject,content,filePath)
%调用格式为：1.sendEmail(subject,content)
%???????????2.sendEmail(subject,content,filePath)
%???Exp1:sendEmail('你的MATLAB代码执行完毕啦！',strcat('代码执行所花的时间为(s)：',num2str(totalTime)));
%???Exp2:sendEmail('你的MATLAB代码执行完毕啦！',strcat('代码执行所花的时间为(s)：',num2str(totalTime)),filePath);
%???subject:为邮件的主题
%???content：为邮件的内容
%???filePath：附件的路径(要包含具体的文件名)
MailAddress='zx4612@mail.ustc.edu.cn';%此处填写ustc邮箱账号
password='zx19961012';%此处填写密码
setpref('Internet','E_mail',MailAddress);
setpref('Internet','SMTP_Server','mail.ustc.edu.cn');
setpref('Internet','SMTP_Username',MailAddress);
setpref('Internet','SMTP_Password',password);
props=java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
if nargin == 2
    sendmail(MailAddress,subject,content);
elseif nargin==3
    sendmail(MailAddress,subject,content,filePath);
elseif nargin > 3 
    error('Too?many?input?arguments');
elseif nargin <2
    error('Too?less?input?arguments');
end
end
