function QC_wrapper(step,inputfile, newmaskname, Npcs)

if ischar(step)
    step = str2num(step);
end
if ischar(Npcs)
    Npcs = str2num(Npcs);
end

if step==0 || step==1
    Pipeline_QC1( inputfile )
end
if step==0 || step==2
    Pipeline_QC2( inputfile , newmaskname, Npcs)
end

