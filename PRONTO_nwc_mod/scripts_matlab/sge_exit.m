function sge_exit(code,str)

if nargin==1
    str = 'ERROR';
end
    
if (usejava('desktop') || usejava('jvm'))
    error(str);
else
    display(str);
    exit(code);
end