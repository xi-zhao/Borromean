function sendEmail(subject,content,filePath)
%���ø�ʽΪ��1.sendEmail(subject,content)
%???????????2.sendEmail(subject,content,filePath)
%???Exp1:sendEmail('���MATLAB����ִ���������',strcat('����ִ��������ʱ��Ϊ(s)��',num2str(totalTime)));
%???Exp2:sendEmail('���MATLAB����ִ���������',strcat('����ִ��������ʱ��Ϊ(s)��',num2str(totalTime)),filePath);
%???subject:Ϊ�ʼ�������
%???content��Ϊ�ʼ�������
%???filePath��������·��(Ҫ����������ļ���)
MailAddress='zx4612@mail.ustc.edu.cn';%�˴���дustc�����˺�
password='zx19961012';%�˴���д����
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
