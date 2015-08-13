# -*- coding: cp1252 -*-
from Tkinter import Tk, Button,OptionMenu, Toplevel, Message, Checkbutton, Frame, Label, LabelFrame, Entry, Menu, END, LEFT, RIGHT, BOTH, Y, INSERT,IntVar, StringVar, W, E, S, N, Listbox, Scrollbar, VERTICAL, HORIZONTAL, SUNKEN
from tkFileDialog import askopenfilename, askdirectory, asksaveasfilename, askopenfilenames
import tkMessageBox
import ScrolledText
import os,sys
import numpy
from  scipy.io import loadmat, savemat

class CreateToolTip(object):
    '''
    create a tooltip for a given widget
    '''
    def __init__(self, widget, text='widget info'):
        self.widget = widget
        self.text = text
        self.widget.bind("<Enter>", self.enter)
        self.widget.bind("<Leave>", self.close)
    def enter(self, event=None):
        x = y = 0
        x, y, cx, cy = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 20
        # creates a toplevel window
        self.tw = Toplevel(self.widget)
        # Leaves only the label and removes the app window
        self.tw.wm_overrideredirect(True)
        self.tw.wm_geometry("+%d+%d" % (x, y))
        label = Label(self.tw, text=self.text, justify='left',
                       background='yellow', relief='solid', borderwidth=1,
                       font=("times", "10", "normal"))
        label.pack(ipadx=1)
    def close(self, event=None):
        if self.tw:
            self.tw.destroy()
def isfloat(s):
    o = s.split(".")
    if len(o)>2: return False
    if len(o)==0: return True
    if not o[0].isdigit(): return False
    if len(o)==1:return True
    if not o[1]: return True
    if not o[1].isdigit(): return False
    return True
def removeempty(sin):
    o = []
    s = sin.replace(" ",",")
    s = s.split(",")
    for i in s:
        if i: o.append(i)
    so = ",".join(o)
    return so
    
def convert_list(l,k):
    b = []
    for i in l:
        for j in i:
            b.append(j+k)
    o = ",".join([str(x) for x in b])
    return o
class InputLines(object):
    def __init__(self):
        self.sin        = []
        self.out        = []
        self.prefix     = []
        self.physio     = []
        self.customreg  = []
        self.structrual = []
        self.taskinfo   = []
        self.drop       = []
        self.tr_msec = []
        self.xtype = []
        self.block = []
        self.seed = []
        self.unit = []
        self.tasktext = []
        self.conditions = []

    def append(self,sin='',out='',prefix='',physio='',customreg='',structrual='',taskinfo='',drop='',tr_msec='',xtype='',block='',seed='',tasktext='',unit='', conditions={}):
        self.sin.append(sin)
        self.out.append(out)
        self.prefix.append(prefix)
        self.physio.append(physio)
        self.customreg.append(customreg)
        self.structrual.append(structrual)
        self.taskinfo.append(taskinfo)
        self.drop.append(drop)
        self.tr_msec.append(tr_msec)
        self.xtype.append(xtype)
        self.block.append(block)
        self.seed.append(seed)
        self.tasktext.append(tasktext)
        self.unit.append(unit)
        self.conditions.append(conditions.copy())


    def edit(self,number,sin='',out='',prefix='',physio='',customreg='',structrual='',taskinfo='',drop='',tr_msec='',xtype='',block='',seed='',tasktext='',unit='',conditions={}):
        if sin: self.sin[number] = sin
        if out: self.out[number] = out
        if prefix: self.prefix[number] = prefix
        if physio: self.physio[number] = physio
        if customreg: self.customreg[number]   = customreg
        if structrual: self.structrual[number] = structrual
        if taskinfo: self.taskinfo[number] = taskinfo
        if drop: self.drop[number]        = drop
        if tr_msec: self.tr_msec[number]  = tr_msec
        if xtype: self.xtype[number]    = xtype
        if block: self.block[number]    = block
        if seed: self.seed[number]      = seed
        if unit: self.unit[number]      = unit
        if tasktext: self.tasktext[number] = tasktext
        if conditions: self.conditions[number]      = conditions.copy()
        

    def newitem(self):
        self.append()

    def delete(self,number):
        del self.sin[number]
        del self.out[number]
        del self.prefix[number]
        del self.physio[number]
        del self.customreg[number]
        del self.structrual[number]
        del self.taskinfo[number]
        del self.drop[number]
        del self.tr_msec[number]
        del self.xtype[number]
        del self.block[number]
        del self.seed[number]
        del self.tasktext[number]
        del self.unit[number]
        del self.conditions[number]
                
    def delete_last(self):
        self.delete(len(self.sin))
        
    def make_copy(self,number):
        self.sin.append(self.sin[number])
        self.out.append(self.out[number])
        self.prefix.append(self.prefix[number])
        self.physio.append(self.physio[number])
        self.customreg.append(self.customreg[number])
        self.structrual.append(self.structrual[number])
        self.drop.append(self.drop[number])
        self.tr_msec.append(self.tr_msec[number])
        self.xtype.append(self.xtype[number])
        self.block.append(self.block[number])
        self.seed.append(self.seed[number])
        self.tasktext.append(self.tasktext[number])       
        self.unit.append(self.unit[number])       
        self.taskinfo.append(self.taskinfo[number])       
        self.conditions.append(self.conditions[number].copy())       
      

    def clear(self):
        self.sin        = []
        self.out        = []
        self.prefix     = []
        self.physio     = []
        self.customreg  = []
        self.structrual = []
        self.taskinfo   = []
        self.drop       = []
        self.tr_msec = []
        self.xtype = []
        self.block = []
        self.seed = []
        self.tasktext = []
        self.unit      = []
        self.conditions = []
        
class Files(Frame):
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.parent = master
        #self.pack(side=LEFT)
        self.Frame1 = Frame(master,bd=2, relief=SUNKEN)
        self.Frame1.place(relx=0.17, rely=0.015, relheight=0.5, relwidth=0.605)
        self.Menus()

        self.g=0

        # Scroll bar and list box

        self.Frame2 = Frame(master)
        self.Frame2.place(relx=0.01, rely=0.15, relheight=0.45, relwidth=0.151)
        self.scrollbar = Scrollbar(self.Frame2, orient=VERTICAL)
        self.listbox = Listbox(self.Frame2,yscrollcommand=self.scrollbar.set,exportselection=0)
        self.scrollbar.config(command=self.listbox.yview)
        self.listbox.place(relx=0.01, rely=0.05, relheight=0.70, relwidth=0.8)
        self.scrollbar.place(relx=0.81, rely=0.05, relheight=0.70, relwidth=0.1)

        
        self.listbox.bind('<<ListboxSelect>>', self.listbox_callback)
        self.menu = Menu(self.Frame2, tearoff=0)
        self.menu.add_command(label="New subject/run", command=self.Add_subject)
        self.menu.add_command(label="Add PLS SessionFile(s)",  command=self.LoadPLSSessionFiles)
        self.menu.add_command(label="Add a Input TextFile",  command=self.AddInputFile)
        self.menu.add_command(label="Make a copy", command=self.Make_copy_subject)
        self.menu.add_command(label="Delete subject/run",  command=self.Delete_subject)

        self.listbox.bind("<Button-3>", self.popup)
        self.listbox.select_set(0)


        self.lin  = Label(self.Frame1, text="INPUT", anchor=E)
        self.lin.grid(row=0,column=1)
        self.lout = Label(self.Frame1, text="OUTPUT directory",anchor=E)
        self.lout.grid(row=1,column=1)
        self.lprefix = Label(self.Frame1, text="OUTPUT prefix",anchor=E)
        self.lprefix.grid(row=2,column=1)
        self.lphysio = Label(self.Frame1, text="PHYSIO",anchor=E)
        self.lphysio.grid(row=3,column=1)
        self.lcustomreg = Label(self.Frame1, text="CUSTOMREG", anchor=E)
        self.lcustomreg.grid(row=5,column=1)
        self.lstructrual = Label(self.Frame1, text="STRUCT",anchor=E)
        self.lstructrual.grid(row=4,column=1)
        self.ltaskinfo   = Label(self.master, text="Task Info", anchor=E)
        self.ltaskinfo.place(relx=0.470,rely=0.045,relheight=0.05, relwidth=0.05)
        self.ldrop = Label(self.Frame1, text="DROP",anchor=E)
        self.ldrop.grid(row=6,column=1)
        self.itisme = True
        self.str_in  = StringVar()
        self.str_out = StringVar()
        self.str_prefix = StringVar()
        self.str_physio = StringVar()
        self.str_customreg = StringVar()
        self.str_structrual = StringVar()
        self.str_taskinfo   = StringVar()
        self.str_cdetrend   = StringVar()
        self.str_csmooth   = StringVar()
        self.str_modelmenu = StringVar()
        self.str_parametermenu = StringVar()
        self.str_contrast  = StringVar()
        self.str_checkbox_keepmean = StringVar()
        self.str_voxelsize                =StringVar()
        self.str_checkbox_deoblique = StringVar()
        self.str_checkbox_convolve  = StringVar()
        self.str_tpatternmenu       = StringVar()
        self.str_cunit              = StringVar()
        self.str_ctrmsec            = StringVar()
        self.str_ctype              = StringVar()
        self.str_ccensor            =  StringVar()
        self.str_consets        = StringVar()
        self.str_cblocklength        = StringVar()
        self.str_cseed             = StringVar()
        self.str_cname             = StringVar()

        self.str_addcondition      = StringVar()

        self.str_numcores         = StringVar()
        self.str_numcores.set('1')
        self.str_queue         = StringVar()
        self.str_queue.set('bigmem_16.q')
        self.str_environment         = StringVar()
        self.str_environment.set('matlab')
        self.int_autodetect    = IntVar()
        self.int_autodetect.set(0)
        self.str_nresample    = StringVar()
        self.str_nresample.set('10')
        self.str_submissionnode = StringVar()
        self.str_submissionnode.set('headnode')
        
        
        self.str_parameter_decision_model = StringVar()
        self.str_parameter_drf            = StringVar()
        self.str_parameter_Nblock         = StringVar()
        self.str_parameter_WIND           = StringVar()
        self.str_parameter_subspace       = StringVar()
        self.str_parameter_spm            = StringVar()
        self.str_parameter_null            = StringVar()
        self.str_metricmenu               =StringVar()
        self.str_partmenu               = StringVar()
        self.str_referene               = StringVar()
        #self.str_drop                   = StringVar()
        self.str_parameter_decision_model.set("linear")
        self.str_parameter_drf.set("0.5")
        self.str_parameter_Nblock.set("4")
        self.str_parameter_WIND.set("10")
        self.str_parameter_subspace.set("onecomp")
        self.str_parameter_spm.set("zcore")
        self.str_parameter_null.set("")
        self.str_partmenu.set("Run all parts of full optimization pipeline")
        self.str_cunit.set("TR")
        
        self.str_metricmenu.set("dPR")
        self.str_parametermenu.set("drf")
        self.str_modelmenu.set("LDA")
        self.str_ctype.set("block")
        
        self.parameter_dict={'None':[],'LDA':['drf'],'GNB':['decision_model'],'erCVA':['Nblock','WIND','drf','subspace'],'erGNB':['Nblock','WIND'],'SCONN':['spm'],'gPCA':[],'GLM':[],'erGLM':[]}
        self.design_dict={'block':['GNB','LDA','GLM','None'],'event':['erGNB','erCVA','erGLM','None'],'nocontrast':['SCONN','gPCA','None'],'':['None']}
        self.parameters_rf={'drf':self.str_parameter_drf,'decision_model':self.str_parameter_decision_model,'Nblock':self.str_parameter_Nblock,'WIND':self.str_parameter_WIND,'subspace':self.str_parameter_subspace,'spm':self.str_parameter_spm}
        
        self.str_parametermenu.trace("w",self.parametermenu_callback)
        self.str_modelmenu.trace("w",self.modelmenu_callback)
        self.str_ctype.trace("w",self.design_type_call)
        self.str_cname.trace("w",self.condition_callback)

                               
        self.previously_selected_item = 0
        self.current_item             = 0
        self.list_entries = InputLines();
        self.list_entries.clear()


        '''
        for i in range(0,1):

            self.list_entries.append(sin='in subject %03d' % i,out='out subject %03d' % i,prefix='prefix subject %03d' % i, \
                                     physio='physio subject %03d' % i,customreg='customreg subject %03d' % i, \
                                     structrual='structrual subject %03d' % i, taskinfo= '', \
                                     drop = "[0,0]",conditions={'1':['23 18 19','54 96 27'],'2':['23 0 19','54 0 27']}, \
                                     tr_msec='2000',unit='TR')
        '''
        self.ein  = Entry(self.Frame1,textvariable=self.str_in)
        self.ein.grid(row=0, column=2)
        self.eout = Entry(self.Frame1, textvariable=self.str_out)
        self.eout.grid(row=1, column=2)
        self.eprefix =Entry(self.Frame1, textvariable=self.str_prefix)
        self.eprefix.grid(row=2, column=2)
        self.ephysio = Entry(self.Frame1, textvariable=self.str_physio)
        self.ephysio.grid(row=3, column=2)
        self.customreg = Entry(self.Frame1, textvariable=self.str_customreg)
        self.customreg.grid(row=5, column=2)
        self.structural = Entry(self.Frame1, textvariable=self.str_structrual)
        self.structural.grid(row=4, column=2)
        self.taskinfo   = Entry(self.master, textvariable=self.str_taskinfo)
        self.taskinfo.place(relx=0.52,rely=0.04,relheight=0.06, relwidth=0.19)
        vdrop = self.register(self.validate_drop_entry)

        self.drop = Entry(self.Frame1,text="0,0", validate="key",validatecommand=(vdrop,"%P"))
        self.drop.grid(row=6, column=2)
        self.drop.insert(0,"0,0")
        #self.str_drop.set("0,0")

        self.bin = Button(self.Frame1,
                        text = "Browse", command=self.BrowseNiftiFile, justify=RIGHT)
        self.bin.grid(row=0,column=3)

        self.bout = Button(self.Frame1,
                        text = "Browse", command=self.BrowseDir,justify=RIGHT)
        self.bout.grid(row=1,column=3)

        self.bphysio = Button(self.Frame1,
                        text = "Browse", command=self.BrowsePhysio,justify=RIGHT)
        self.bphysio.grid(row=3,column=3)

        self.bcustomreg  = Button(self.Frame1,
                        text = "Browse", command=self.BrowseCustomReg,justify=RIGHT)
        self.bcustomreg.grid(row=5,column=3)
        
        self.bstructrual = Button(self.Frame1,
                        text = "Browse", command=self.BrowseStructural,justify=RIGHT)
        self.bstructrual.grid(row=4,column=3)

        self.btaskinfo = Button(self.master,
                        text = "Load", command=self.BrowseTaskInfo,justify=RIGHT)
        self.btaskinfo.place(relx=0.72,rely=0.04, relheight=0.06, relwidth=0.05)
        
        # 1 STEP 1
        self.lstep1 = Label(self.master, text="1:Input Files",  anchor=W)
        self.lstep1.place(relx=0.01, rely=0.01, relheight=0.06, relwidth=0.15)
       

        # text editor
        self.textpad = ScrolledText.ScrolledText(master)
        #self.textpad.place(relx=0.47,rely=0.22, relheight=0.25, relwidth=0.3)                                      

        self.lunit = Label(self.master,text="UNIT:", anchor=E)
        self.lunit.place(relx=0.455,rely=0.145, relheight=0.05, relwidth=0.05)
        self.cunits = OptionMenu(self.master,self.str_cunit,"TR","msec","sec")
        self.cunits.place(relx=0.503,rely=0.135, relheight=0.07, relwidth=0.065)

        self.ltr_msec = Label(self.master,text="TR_MSEC:", anchor=E)
        self.ltr_msec.place(relx=0.57,rely=0.145, relheight=0.05, relwidth=0.055)
        vtrmsec = self.register(self.validate_trmsec_entry)
        self.ctrmsec = Entry(self.master,validate="key", validatecommand=(vtrmsec,'%P'))
        self.ctrmsec.place(relx=0.623,rely=0.14, relheight=0.06, relwidth=0.05)

        self.ltype = Label(self.master,text="TYPE", anchor=E)
        self.ltype.place(relx=0.01, rely=0.105, relheight=0.03, relwidth=0.03)
        self.ctype = OptionMenu(self.master,self.str_ctype,"block","event","nocontrast")
        self.ctype.place(relx=0.04, rely=0.09, relheight=0.06, relwidth=0.1)


        self.lconditions = Label(self.master,text="Condition(s):", anchor=E)
        self.lconditions.place(relx=0.475,rely=0.235, relheight=0.05, relwidth=0.072)
        self.cconditions = OptionMenu(self.master,self.str_cname,"",command=self.condition_callback)
        self.cconditions.place(relx=0.545,rely=0.2250, relheight=0.07, relwidth=0.12)


        self.optionmenu_contex = Menu(self.master, tearoff=0)
        self.optionmenu_contex.add_command(label="Delete",  command=self.DeleteCondition)
        self.optionmenu_contex.add_command(label="Add", command=self.AddCondition)
        self.optionmenu_contex.add_command(label="Rename",  command=self.RenameCondition)
        #self.optionmenu_contex.add_separator()
        #self.optionmenu_contex.add_command(label="Copy",  command=self.CopyCondition)
        #self.optionmenu_contex.add_command(label="Paste",  command=self.PasteCondition)
        self.cconditions.bind("<Button-3>", self.popup_optionmenu_contex)
        self.cconditions.unbind("<Button-3>")
        
        self.bdeleteconditions = Button(self.master,text = "Delete", command=self.DeleteCondition,justify=RIGHT)
        self.bdeleteconditions.place(relx=0.667,rely=0.2260, relheight=0.06, relwidth=0.05)
        self.baddconditions=Button(self.master,text = "Add", command=self.AddCondition,justify=RIGHT)
        self.baddconditions.place(relx=0.72,rely=0.2260, relheight=0.06, relwidth=0.05)

        self.lonset = Label(self.master,text="Onsets:", anchor=E)
        self.lonset.place(relx=0.475,rely=0.315, relheight=0.05, relwidth=0.072)
        vonsets = self.register(self.validate_onsets_entry)
        self.consets = Entry(self.master,validate="key",validatecommand=(vonsets,'%P','%V'))
        self.consets.place(relx=0.547,rely=0.32, relheight=0.05, relwidth=0.22)

        self.lduration=Label(self.master,text="Duration:", anchor=E)
        self.lduration.place(relx=0.475,rely=0.375, relheight=0.05, relwidth=0.072)
        vblocklength = self.register(self.validate_blocklength_entry)
        self.cblocklength = Entry(self.master,validate="key",validatecommand=(vblocklength,'%P','%V'))
        self.cblocklength.place(relx=0.546,rely=0.375, relheight=0.05, relwidth=0.22)

        self.lseed = Label(self.master,text="Seed:", anchor=E)
        self.lseed.place(relx=0.475,rely=0.4350, relheight=0.05, relwidth=0.072)
        self.cseed = Entry(self.master,textvariable=self.str_cseed)
        self.cseed.place(relx=0.547,rely=0.435, relheight=0.05, relwidth=0.17)
        self.bseed = Button(self.master,text = "Browse", command=self.BrowseSeed,justify=RIGHT)
        self.bseed.place(relx=0.72,rely=0.430, relheight=0.06, relwidth=0.047)


        
        self.drop = Entry(self.Frame1,text="0,0", validate="key",validatecommand=(vdrop,"%P"))
        self.drop.grid(row=6, column=2)
        self.drop.insert(0,"0,0")

        self.listbox_callback([])

        # Pipeline choices
        self.frame_pipelinelist = Frame(master,bd=2, relief=SUNKEN)
        self.frame_pipelinelist.place(relx=0.8,rely=0.007,relheight=0.95, relwidth=0.19)
        
        self.lstep2  = Label(self.frame_pipelinelist, text="2: Pipeline list", anchor=E)
        self.lstep2.place(relx=0.25,rely=0.01,relheight=0.08, relwidth=0.43)
        
        self.lmotcor   = Label(self.frame_pipelinelist, text="MOTCOR", anchor=E)
        self.lmotcor.place(relx=0.13,rely=0.1,relheight=0.08, relwidth=0.4)
        self.cmotcor   = Button(self.frame_pipelinelist,text = "ON/OFF", command=lambda: self.MotcorState(0))
        self.cmotcor.place(relx=0.54,rely=0.1,relheight=0.07, relwidth=0.35)
        
        self.lcensor   = Label(self.frame_pipelinelist, text="CENSOR", anchor=E)
        self.lcensor.place(relx=0.13,rely=0.17,relheight=0.08, relwidth=0.4)
        vcensor = self.register(self.validate_censor_entry)
        self.ccensor   = Entry(self.frame_pipelinelist,text="0,1",validate="key",validatecommand=(vcensor,"%P"))
        self.ccensor.place(relx=0.54,rely=0.17,relheight=0.07, relwidth=0.35)
        self.str_ccensor.set("0,1")
        self.ccensor.insert(0,"0,1")

        self.lretroicor   = Label(self.frame_pipelinelist, text="RETROICOR", anchor=E)
        self.lretroicor.place(relx=0.13,rely=0.24,relheight=0.08, relwidth=0.4)
        self.cretroicor   = Button(self.frame_pipelinelist,text = "ON/OFF", command=lambda: self.MotcorState(2))
        self.cretroicor.place(relx=0.54,rely=0.24,relheight=0.07, relwidth=0.35)
        
        self.ltimecor   = Label(self.frame_pipelinelist, text="TIMECOR", anchor=E)
        self.ltimecor.place(relx=0.13,rely=0.31,relheight=0.08, relwidth=0.4)
        self.ctimecor   = Button(self.frame_pipelinelist,text = "ON/OFF", command=lambda: self.MotcorState(3))
        self.ctimecor.place(relx=0.54,rely=0.31,relheight=0.07, relwidth=0.35)
        
        self.lsmooth   = Label(self.frame_pipelinelist, text="SMOOTH", anchor=E)
        self.lsmooth.place(relx=0.13,rely=0.38,relheight=0.08, relwidth=0.4)
        vsmooth = self.register(self.validate_smooth_entry)
        self.csmooth   = Entry(self.frame_pipelinelist,text="6",validate="key",validatecommand=(vsmooth,"%P"))
        self.csmooth.place(relx=0.54,rely=0.38,relheight=0.07, relwidth=0.35)
        self.str_csmooth.set("6")
        self.csmooth.insert(0,"6")

        
        self.ldetrend   = Label(self.frame_pipelinelist, text="DETREND", anchor=E)
        self.ldetrend.place(relx=0.13,rely=0.45,relheight=0.08, relwidth=0.4)
        vdetrend = self.register(self.validate_detrend_entry)
        self.cdetrend   = Entry(self.frame_pipelinelist,text="0,1,2,3,4,5",validate="key",validatecommand=(vdetrend,"%P"))
        self.cdetrend.place(relx=0.54,rely=0.45,relheight=0.07, relwidth=0.35)
        self.str_cdetrend.set("0,1,2,3,4,5")
        self.cdetrend.insert(0,"0,1,2,3,4,5")


        self.ltask   = Label(self.frame_pipelinelist, text="TASK", anchor=E)
        self.ltask.place(relx=0.13,rely=0.52,relheight=0.08, relwidth=0.4)
        self.ctask   = Button(self.frame_pipelinelist,text = "ON/OFF", command=lambda: self.MotcorState(6))
        self.ctask.place(relx=0.54,rely=0.52,relheight=0.07, relwidth=0.35)
        
        self.lmotreg   = Label(self.frame_pipelinelist, text="MOTREG", anchor=E)
        self.lmotreg.place(relx=0.13,rely=0.59,relheight=0.08, relwidth=0.4)
        self.cmotreg   = Button(self.frame_pipelinelist,text = "ON/OFF", command=lambda: self.MotcorState(7))
        self.cmotreg.place(relx=0.54,rely=0.59,relheight=0.07, relwidth=0.35)
        

        self.lgspc1   = Label(self.frame_pipelinelist, text="GSPC1", anchor=E)
        self.lgspc1.place(relx=0.13,rely=0.66,relheight=0.08, relwidth=0.4)
        self.cgspc1   = Button(self.frame_pipelinelist,text = "ON/OFF", command=lambda: self.MotcorState(8))
        self.cgspc1.place(relx=0.54,rely=0.66,relheight=0.07, relwidth=0.35)

        self.lphyplus   = Label(self.frame_pipelinelist, text="PHYPLUS", anchor=E)
        self.lphyplus.place(relx=0.13,rely=0.73,relheight=0.08, relwidth=0.4)
        self.cphyplus   = Button(self.frame_pipelinelist,text = "ON/OFF", command=lambda: self.MotcorState(9))
        self.cphyplus.place(relx=0.54,rely=0.73,relheight=0.07, relwidth=0.35)

        self.llowpass   = Label(self.frame_pipelinelist, text="LOWPASS", anchor=E)
        self.llowpass.place(relx=0.13,rely=0.80,relheight=0.08, relwidth=0.4)
        self.clowpass   = Button(self.frame_pipelinelist,text = "ON/OFF", command=lambda: self.MotcorState(10))
        self.clowpass.place(relx=0.54,rely=0.80,relheight=0.07, relwidth=0.35)

        self.lcustomreg2   = Label(self.frame_pipelinelist, text="CUSTOMREG", anchor=E)
        self.lcustomreg2.place(relx=0.13,rely=0.87,relheight=0.08, relwidth=0.4)
        self.ccustomreg   = Button(self.frame_pipelinelist,text = "ON/OFF", command=lambda: self.MotcorState(11))
        self.ccustomreg.place(relx=0.54,rely=0.87,relheight=0.07, relwidth=0.35)
        
    # run frame
        
            
        self.Frame_run = Frame(master,bd=2, relief=SUNKEN)
        self.Frame_run.place(relx=0.01, rely=0.57, relheight=0.41, relwidth=0.46)
        
        self.lanalysismodel   = Label(self.Frame_run, text="Analysis Model:", anchor=E)
        self.lanalysismodel.place(relx=0.00,rely=0.03,relheight=0.15, relwidth=0.23)
        self.parameters   = Entry(self.Frame_run,textvariable=self.str_parameter_drf)
        self.parameters.place(relx=0.49,rely=0.23,relheight=0.14, relwidth=0.13)
  
    
        self.lmodelparameter   = Label(self.Frame_run, text="Model Parameters:", anchor=E)
        self.lmodelparameter.place(relx=0.00,rely=0.22,relheight=0.15, relwidth=0.23)
        
        self.parametermenu = OptionMenu(self.Frame_run, self.str_parametermenu, "drf")
        self.parametermenu.place(relx=0.225, rely=0.21, relheight=0.2, relwidth=0.26)
        
        
        self.lmetric   = Label(self.Frame_run, text="Metric:", anchor=E)
        self.lmetric.place(relx=0.00,rely=0.41,relheight=0.15, relwidth=0.23)
        

        self.lcontrast   = Label(self.Frame_run, text="Contrast:", anchor=E)
        self.lcontrast.place(relx=0.00,rely=0.63,relheight=0.15, relwidth=0.23)
        self.contrast   = Entry(self.Frame_run,textvariable=self.str_contrast)
        self.contrast.place(relx=0.23,rely=0.63,relheight=0.15, relwidth=0.4)

        self.lvoxelsize   = Label(self.Frame_run, text="Voxel size:", anchor=E)
        self.lvoxelsize.place(relx=0.00,rely=0.83,relheight=0.15, relwidth=0.23)
        vvoxelsize = self.register(self.validate_voxelsize_entry)
        self.voxelsize   = Entry(self.Frame_run,validate="key",validatecommand=(vvoxelsize,"%P"))
        self.voxelsize.place(relx=0.23,rely=0.83,relheight=0.15, relwidth=0.4)
        
        self.keepmean_checkbox = Checkbutton(self.Frame_run, text="Keep mean", variable=self.str_checkbox_keepmean,justify=LEFT)
        self.keepmean_checkbox.place(relx=0.68,rely=0.22,relheight=0.15, relwidth=0.18)
        self.str_checkbox_keepmean.set(0)

        self.deoblique_checkbox = Checkbutton(self.Frame_run, text="DEOBLIQUE", variable=self.str_checkbox_deoblique,justify=LEFT)
        self.deoblique_checkbox.place(relx=0.685,rely=0.37,relheight=0.15, relwidth=0.18)
        self.str_checkbox_deoblique.set(0)

        self.convolv_checkbox = Checkbutton(self.Frame_run, text="Convolve", variable=self.str_checkbox_convolve,justify=LEFT)
        self.convolv_checkbox.place(relx=0.670,rely=0.52,relheight=0.15, relwidth=0.18)
        self.str_checkbox_convolve.set(0)

        
        self.str_tpatternmenu.set("None")
        self.ltpattern   = Label(self.Frame_run, text="TPATTERN:", anchor=E)
        self.ltpattern.place(relx=0.67,rely=0.04,relheight=0.15, relwidth=0.15)
        self.tpatternmenu = OptionMenu(self.Frame_run, self.str_tpatternmenu, "None","altplus","altimus","seqplus","seqminus")
        self.tpatternmenu.place(relx=0.814, rely=0.01, relheight=0.2, relwidth=0.184)

        self.modelmenu = OptionMenu(self.Frame_run, self.str_modelmenu, *self.parameter_dict.keys())  
        self.modelmenu.place(relx=0.225, rely=0.01, relheight=0.2, relwidth=0.26)
        
        self.metricmenu = OptionMenu(self.Frame_run, self.str_metricmenu, "dPR", "R", "P")
        self.metricmenu.place(relx=0.225, rely=0.41, relheight=0.2, relwidth=0.26)

    # submission Frame
        self.Frame_submission = Frame(master,bd=2, relief=SUNKEN)
        self.Frame_submission.place(relx=0.475, rely=0.57, relheight=0.41, relwidth=0.31)

        self.lselectpipelinepart   = Label(self.Frame_submission, text="Select Pipeline Optimization Part:", anchor=E)
        self.lselectpipelinepart.place(relx=0.00,rely=0.0,relheight=0.15, relwidth=0.75)

        self.partmenu = OptionMenu(self.Frame_submission, self.str_partmenu, "Run all parts of full optimization pipeline", "Part 1: run all pipelines and produce metrics", "Part 2: select optimal pipelines based on metrics","Part 3: do spatial normalization of optimized results")
        self.partmenu.place(relx=0.0, rely=0.15, relheight=0.2, relwidth=1)


        self.lselectreferecevolum   = Label(self.Frame_submission, text="Select Reference Volume:", anchor=E)
        self.lselectreferecevolum.place(relx=0.00,rely=0.45,relheight=0.15, relwidth=0.75)
        self.reference   = Entry(self.Frame_submission,textvariable=self.str_referene)
        self.reference.place(relx=0.01,rely=0.6,relheight=0.14, relwidth=0.75)
  
        self.breference  = Button(self.Frame_submission,text = "Browse", command=self.ReferenceCallback)
        self.breference.place(relx=0.77,rely=0.59,relheight=0.17, relwidth=0.22)
        

        self.breference  = Button(self.Frame_submission,text = "Save and Submit", command=self.SubmitCallback)
        self.breference.place(relx=0.01,rely=0.83,relheight=0.17, relwidth=0.99)


        self.lstep3 = Label(self.master, text="3:Select AnalysisModel",  anchor=W)
        self.lstep3.place(relx=0.01, rely=0.52, relheight=0.05, relwidth=0.15)
        self.lstep4 = Label(self.master, text="4:Run",  anchor=W)
        self.lstep4.place(relx=0.48, rely=0.52, relheight=0.05, relwidth=0.15)


        self.set_tooltip_messages()

    def update_entries(self):
        self.itisme = True
        index = map(int, self.listbox.curselection())
        if not index:
            self.disable_entries()
            return
        self.enable_entries()
        index = index[0]

        self.str_in.set(self.list_entries.sin[index])
        self.str_out.set(self.list_entries.out[index])
        self.str_prefix.set(self.list_entries.prefix[index])
        self.str_physio.set(self.list_entries.physio[index])
        self.str_customreg.set(self.list_entries.customreg[index])
        self.str_structrual.set(self.list_entries.structrual[index])
        self.drop.delete(0,END)    
        self.drop.insert(0,self.list_entries.drop[index])    

        self.str_taskinfo.set(self.list_entries.taskinfo[index])
        self.str_ctype.set(self.list_entries.xtype[index])
        self.ctrmsec.delete(0,END)            
        self.ctrmsec.insert(0,self.list_entries.tr_msec[index])
        self.str_cunit.set(self.list_entries.unit[index])
        self.str_cseed.set(self.list_entries.seed[index])
        
        self.cconditions['menu'].delete(0, 'end')
        for i in self.list_entries.conditions[index].keys():
            self.cconditions['menu'].add_command(label=i, command=lambda i=i:self.str_cname.set(i))                    

        name = self.str_cname.get()
        if not self.str_cname.get() in self.list_entries.conditions[index].keys():
            if self.list_entries.conditions[index].keys():
                name = self.list_entries.conditions[index].keys()[0]
            else:
                name = ''
        self.itisme = True
        self.str_cname.set(name)
        

        
           
        self.textpad.delete("1.0",END)
        self.textpad.insert(INSERT,self.list_entries.tasktext[index])

    def disable_entries(self):

         self.ein.config(state='disabled')
         self.eout.config(state='disabled')
         self.eprefix.config(state='disabled')
         self.ephysio.config(state='disabled')
         self.customreg.config(state='disabled')
         self.structural.config(state='disabled')
         self.taskinfo.config(state='disabled')
         self.drop.config(state='disabled')

         self.bin.config(state='disabled')
         self.bout.config(state='disabled')
         self.bcustomreg.config(state='disabled')
         self.bphysio.config(state='disabled')
         self.bstructrual.config(state='disabled')
         self.btaskinfo.config(state='disabled')
         self.btaskinfo.config(state='disabled')

         self.cunits.config(state='disabled')
         self.cconditions.config(state='disabled')
         self.baddconditions.config(state='disabled')
         self.bdeleteconditions.config(state='disabled')
         self.cseed.config(state='disabled')
         self.bseed.config(state='disabled')
         self.ctrmsec.config(state='disabled')
         self.cblocklength.config(state='disabled')
         self.consets.config(state='disabled')
         self.str_cname.set("")
         self.cconditions['menu'].delete(0, 'end')
         self.cconditions.unbind("<Button-3>")
        
    def enable_entries(self):

         self.ein.config(state='normal')
         self.eout.config(state='normal')
         self.eprefix.config(state='normal')
         self.ephysio.config(state='normal')
         self.customreg.config(state='normal')
         self.structural.config(state='normal')
         self.taskinfo.config(state='normal')
         self.drop.config(state='normal')

         self.bin.config(state='normal')
         self.bout.config(state='normal')
         self.bcustomreg.config(state='normal')
         self.bphysio.config(state='normal')
         self.bstructrual.config(state='normal')
         self.btaskinfo.config(state='normal')
         self.btaskinfo.config(state='normal')

         self.cunits.config(state='normal')
         self.cconditions.config(state='normal')
         self.baddconditions.config(state='normal')
         self.bdeleteconditions.config(state='normal')
         self.cseed.config(state='normal')
         self.bseed.config(state='normal')
         self.ctrmsec.config(state='normal')
         self.cblocklength.config(state='normal')
         self.consets.config(state='normal')
         self.cconditions.bind("<Button-3>", self.popup_optionmenu_contex)
         
    def listbox_callback(self,evt):
 
    # Note here that Tkinter passes an event object to onselect()
        index = map(int, self.listbox.curselection())
        if not index:
            self.disable_entries()
    # disable all the maps
            return
        index = index[0]
        self.enable_entries()
        
        if self.previously_selected_item>-1:

            self.list_entries.edit(self.previously_selected_item, \
                                   sin=self.str_in.get(),out=self.str_out.get(),prefix=self.str_prefix.get(), \
                                   physio=self.str_physio.get(),customreg=self.str_customreg.get(), \
                                   structrual=self.str_structrual.get(),taskinfo=self.str_taskinfo.get(), \
                                   tr_msec=self.str_ctrmsec.get(),unit=self.str_cunit.get(),seed=self.cseed.get(),  \
                                   tasktext=self.textpad.get("1.0",END))

       
        self.previously_selected_item = index
        if index<len(self.list_entries.sin):
            self.update_entries()

    def move_cursor(self,index):
        prev = self.listbox.curselection()
        if prev:
            self.listbox.selection_clear(self.listbox.curselection())
        self.listbox.selection_set(index)
        self.listbox.activate(index)
        self.listbox.yview_moveto(1)

     # move cursor to the end of list
    '''
    def all_inputs_callback(self,*args):
        index = map(int, self.listbox.curselection())
        index = index[0]
        self.list_entries.edit(index, \
                               sin=self.str_in.get(),out=self.str_out.get(),prefix=self.str_prefix.get(), \
                               physio=self.str_physio.get(),customreg=self.str_customreg.get(), \
                               structrual=self.str_structrual.get(),taskinfo=self.str_taskinfo.get(), \
                               tr_msec=self.str_ctrmsec.get(),unit=self.str_cunit.get(),seed=self.cseed.get(),drop=self.str_drop.get(),  \
                               tasktext=self.textpad.get("1.0",END))
    '''
    def validate_trmsec_entry(self,p):
        if p.isdigit() or not p:
            self.str_ctrmsec.set(p)
            return True
        return False
    def validate_onsets_entry(self,pin,vv):
        if self.itisme:
            self.itisme=False
            #return True
        if not pin: return True
        #p = removeempty(pin.split())
        p = pin.replace(" ",",")
        pl = p.split(",")
        for q in pl:
            if not isfloat(q) and q: return False
        name = self.str_cname.get()
        index = map(int, self.listbox.curselection())
        if not index:
            return True
        number = index[0]
        if name:
            self.list_entries.conditions[number][name][0] = pin
        return True
    def  validate_blocklength_entry(self,pin,vv):
        if self.itisme:
            self.itisme=False
            #return True
        if not pin: return True
        p = pin.replace(" ",",")
        #p = removeempty(p.split(","))
        pl = p.split(',')
        for q in pl:
            if not isfloat(q) and q: return False
        name = self.str_cname.get()
        
        index = map(int, self.listbox.curselection())
        if not index:
            return False
        number = index[0]
        if name:
              self.list_entries.conditions[number][name][1] = pin
        return True
    def  validate_drop_entry(self,pin):
        if not pin: return True
        p = pin.replace(" ",",")
        #p = removeempty(p.split(","))
        pl = p.split(',')
        if len(pl)>2:return False
        for q in pl:
            if not q.isdigit() and q: return False

        index = map(int, self.listbox.curselection())
        if not index:
            return False
        number = index[0]
        self.list_entries.drop[number]=p

        return True
    def  validate_smooth_entry(self,pin):
        if not pin: return True
        p = pin.replace(" ",",")
        #p = removeempty(p.split(","))
        pl = p.split(',')
        for q in pl:
            if not isfloat(q) and q: return False
        self.str_csmooth.set(p)
        return True
    def  validate_detrend_entry(self,pin):
        if not pin: return True
        p = pin.replace(" ",",")
        #p = removeempty(p.split(","))
        pl = p.split(',')
        for q in pl:
            if not q.isdigit() and q: return False
        self.str_cdetrend.set(p)
        return True
    def  validate_censor_entry(self,pin):
        if not pin: return True
        p = pin.replace(" ",",")
        #p = removeempty(p.split(","))
        pl = p.split(',')
        for q in pl:
            if not q.isdigit() and q: return False
            if q and int(q)>3:return False
        self.str_ccensor.set(p)
        return True
    def  validate_voxelsize_entry(self,pin):
        if not pin: return True
        p = pin.replace(" ",",")
        #p = removeempty(p.split(","))
        pl = p.split(',')
        if len(pl)>3: return False
        for q in pl:
            if not isfloat(q) and q: return False
        self.str_voxelsize.set(p)
        return True

    
    def taskinfo_callback(self):
        index = map(int, self.listbox.curselection())
        number = index[0]
        self.list_entries.edit(number,taskinfo=self.str_taskinfo.get())
        self.load_task_info(number)
        self.str_cname.set(self.list_entries.conditions[number].keys()[0])
                               
    def popup(self,event):
        self.menu.post(event.x_root, event.y_root)
    def popup_optionmenu_contex(self,event):
        self.optionmenu_contex.post(event.x_root, event.y_root)
        
    def Add_subject(self):
        self.itisme = True
        x={}
        if self.list_entries.sin:
            cond = self.list_entries.conditions[0]
            for v in cond.keys():
                x[v]=['','']
            self.list_entries.append(drop="0,0",conditions=x.copy())
        else:
            self.list_entries.append(drop="0,0",conditions={})
    
        if self.listbox.curselection():
            self.listbox.selection_clear(self.listbox.curselection())

        self.listbox.insert(END, "%03d-subject/run" % len(self.list_entries.sin))
        self.move_cursor(END)

        self.update_entries()
        
    def Make_copy_subject(self):
        index = self.listbox.curselection()
        if not index:
            return
        index = int(index[0])
        self.list_entries.make_copy(index)
        self.listbox.insert(END, self.listbox.get(index))
        self.move_cursor(END)

        self.update_entries()
 
        
    def Delete_subject(self):
        index = self.listbox.curselection()
        if not index:
            return

        self.previously_selected_item = -1
        self.listbox.delete(index)
        index = int(index[0])

        self.list_entries.delete(index)
        
        if index==len(self.list_entries.sin):
            index = len(self.list_entries.sin)-1

        self.move_cursor(index)
        
        self.previously_selected_item = -1
        
        self.update_entries()


    

# Allow the class to run stand-alone.
    def BrowseNiftiFile(self):
        y = askopenfilename(filetypes=[("functional Nifti file","*.*")])
        index = self.listbox.curselection()
        
        index = int(index[0])
        if y:
            self.str_in.set(y)
            if self.str_prefix.get()=="":
                p,n = os.path.split(self.str_in.get())
                self.str_prefix.set("p"+n)
                self.list_entries.prefix[index]="p"+n
            self.list_entries.sin[index] = y
    def BrowseDir(self):
        y = askdirectory()
        index = self.listbox.curselection()
        index = int(index[0])
        if y:
            self.str_out.set(y)

        self.list_entries.out[index] = y
        
    def BrowsePhysio(self):
        y = askopenfilename(filetypes=[("Physio data","*.1D.*")])
        index = self.listbox.curselection()
        index = int(index[0])
        if y:
            self.str_physio.set(y)
        self.list_entries.physio[index] = y
    def BrowseCustomReg(self):
        y = askopenfilename(filetypes=[("Noise ROI","*.nii")])
        index = self.listbox.curselection()
        index = int(index[0])
        if y:
            self.str_customreg.set(y)
        self.list_entries.customreg[index]=y
    def BrowseStructural(self):
        y = askopenfilename(filetypes=[("Structral nifti file","*.nii")])
        index = self.listbox.curselection()
        index = int(index[0])
        if y:
            self.str_structrual.set(y)    
        self.list_entries.structrual[index]=y

    def BrowseSeed(self):
        y = askopenfilename(filetypes=[("Seed nifti file","*.nii")])
        if y:
            self.str_cseed.set(y)
        
    def BrowseTaskInfo(self):
        y = askopenfilename(filetypes=[("Task info","*.txt")])
        index = self.listbox.curselection()
        index = int(index[0])
        if y:
            self.str_taskinfo.set(y)    
            self.list_entries.taskinfo[index]=y
            self.taskinfo_callback()
          

    def BrowsePipelineList(self):
        y = askopenfilename(filetypes=[("Pipeline List","*.txt")])
        if y:
           self.read_pipeline_file(y) 
    
    def SavePipelineList(self):
        y = asksaveasfilename(filetypes=[("Pipeline List","*.txt")])
        if y:
           self.save_pipeline_file(y)

    def LoadDefaultPipelineList(self):
        l = (self.cmotcor, self.ccensor, self.cretroicor, self.ctimecor, self.csmooth, self.cdetrend, self.ctask,self.cmotreg, self.cgspc1, self.cphyplus, self.clowpass, self.ccustomreg)
        steps = ["MOTCOR","CENSOR","RETROICOR","TIMECOR","SMOOTH","DETREND","TASK","MOTREG","GSPC1","PHYPLUS","LOWPASS","CUSTOMREG"]
        d = ["ON/OFF","0,1","OFF","ON/OFF","6","0,1,2,3,4,5","ON/OFF","ON/OFF","ON/OFF","ON/OFF","ON/OFF","OFF"]
        for i in range(0,len(l)):
            s=l[i]
            if i!=4 and i!=5 and i!=1:
                s["text"]=d[i]
            else:
                s.delete(0,END)
                s.insert(0,d[i])
                
    def BrowseInputFile(self):
        y = askopenfilename(filetypes=[("Input file","*.txt")])
        if y:
            self.list_entries.clear()    
            self.listbox.delete(0,END)
            self.cconditions['menu'].delete(0, 'end')
            self.str_cname.set('')
            
            self.read_input_file(y) 
            self.previously_selected_item = -1    
            self.listbox.select_set(0)
    def AddInputFile(self):
        y = askopenfilename(filetypes=[("Input file","*.txt")])
        if y:
           self.read_input_file(y) 
    def SaveInputFile(self):
        y = asksaveasfilename(filetypes=[("Input file","*.txt")])
        if y:
           self.save_input_file(y) 
              
    def NewInputs(self):
        self.list_entries.clear()
        self.cconditions['menu'].delete(0, 'end')
        self.previously_selected_item = -1
        self.listbox.delete(0, END)
        self.disable_entries()
 

    def MotcorState(self,number):
        l = (self.cmotcor, self.ccensor, self.cretroicor, self.ctimecor, self.str_csmooth, self.str_cdetrend, self.ctask,self.cmotreg, self.cgspc1, self.cphyplus, self.clowpass, self.ccustomreg)
        s = l[number]
            
        if s["text"]=="ON":
            s["text"]="OFF"
            return
        if s["text"]=="OFF":
            s["text"]="ON/OFF"
            return
        if s["text"]=="ON/OFF":
            s["text"]="ON"
            return     
    def ReferenceCallback(self):
        y = askopenfilename(filetypes=[("Anatomical Reference nifti file","*.nii")]) 
        if y:
            self.str_referene.set(y)
            
    def modelmenu_callback(self,*args):

        name = self.str_modelmenu.get()

        self.parametermenu['menu'].delete(0, 'end')
        self.metricmenu['menu'].delete(0, 'end')

        plist = self.parameter_dict[name]
        for p in plist:
            self.parametermenu['menu'].add_command(label=p, command=lambda p=p: self.str_parametermenu.set(p))
        if not plist:
            self.str_parametermenu.set("")
        else:
            self.str_parametermenu.set(plist[0])
        
        met=("dPR","P","R")
        for mm in met:
            self.metricmenu['menu'].add_command(label=mm,command=lambda mm=mm: self.str_metricmenu.set(mm))                    
        self.str_metricmenu.set("dPR")

        if name=="GLM" or name=="erGLM":
            self.metricmenu['menu'].delete(0, 'end')
            met=("R")
            for mm in met:
                self.metricmenu['menu'].add_command(label=mm,command=lambda mm=mm: self.str_metricmenu.set(mm))                    
            self.str_metricmenu.set("R")

        if name=="None":
            self.metricmenu['menu'].delete(0, 'end')
            self.str_metricmenu.set("")

        #if name!="SCONN":
         #   print            self.cseed
          #  self.cseed["state"]="disabled"
           # self.bseed["state"]="disabled"
        #else:
         #   self.cseed["state"]="normal"       
          #  self.bseed["state"]="normal"       
        
    def parametermenu_callback(self,*args):

        name = self.str_parametermenu.get()
        if name=="drf":
            self.parameters.config(textvariable=self.str_parameter_drf)
        if name=="Nblock":
            self.parameters.config(textvariable=self.str_parameter_Nblock)
        if name=="WIND":
            self.parameters.config(textvariable=self.str_parameter_WIND)
        if name=="spm":
            self.parameters.config(textvariable=self.str_parameter_spm)
        if name=="decision_model":
            self.parameters.config(textvariable=self.str_parameter_decision_model)
        if name=="subspace":
            self.parameters.config(textvariable=self.str_parameter_subspace)     
        if name=="":
            self.parameters.config(textvariable=self.str_parameter_null)     

    def design_type_call(self,*args):
        name = self.str_ctype.get()
        for i in range(0,len(self.list_entries.sin)):
                             self.list_entries.edit(i,xtype = name)
        
        if hasattr(self,'modelmenu'):
            self.modelmenu['menu'].delete(0, 'end')
        if not name:
            return
        plist = self.design_dict[name]
        for p in plist:
            self.modelmenu['menu'].add_command(label=p, command=lambda p=p: self.str_modelmenu.set(p))
        if not plist:
            self.str_modelmenu.set("")
        else:
            self.str_modelmenu.set(plist[0])
        self.modelmenu_callback()
    
    def condition_callback(self, *args):
        name = self.str_cname.get()
        index = map(int, self.listbox.curselection())
        if not index:
            return
        number = index[0]
        
        if name:
            self.itisme=True
            self.consets.delete(0,END)
            self.itisme=True
            self.cblocklength.delete(0,END)
            self.itisme=True
            self.list_entries.conditions[number][name][0] = removeempty(self.list_entries.conditions[number][name][0])
            self.consets.insert(0,self.list_entries.conditions[number][name][0])
            self.itisme=True
            self.list_entries.conditions[number][name][1] = removeempty(self.list_entries.conditions[number][name][1])
            self.cblocklength.insert(0,self.list_entries.conditions[number][name][1])
            self.cblocklength.config(state='normal')
            self.consets.config(state='normal')
        else:
            self.itisme=True
            self.consets.delete(0,END)
            self.itisme=True
            self.cblocklength.delete(0,END)
            self.cblocklength.config(state='disabled')
            self.consets.config(state='disabled')

      

    def read_pipeline_file(self,pipeline):
        steps = ["MOTCOR","CENSOR","RETROICOR","TIMECOR","SMOOTH","DETREND","TASK","MOTREG","GSPC1","PHYPLUS","LOWPASS","CUSTOMREG"]
        l = (self.cmotcor, self.str_ccensor, self.cretroicor, self.ctimecor, self.str_csmooth, self.str_cdetrend, self.ctask,self.cmotreg, self.cgspc1, self.cphyplus, self.clowpass, self.ccustomreg)
        stron = ['OFF','ON']
        pl = range(0,12)
        with open(pipeline) as f:
            lines = f.readlines()
            for line in lines:
                ctemp=line.split("=")
                temp=ctemp[0].rstrip()
                k     = [i for i in range(0,12) if (steps[i]==temp)]
                if not k:   
                    print "Unknown preprocessing step in line: "+line.rstrip()
                    print "Check file: "+pipeline
                else:
                    s = l[k[0]]
                    pl.remove(k[0])
                    temp = ctemp[1].rstrip()
                    number = temp[1:-1].split(',')
                    if k[0]!=4 and k[0]!=5 and k[0]!=1:
                        if len(number)==2:
                            s["text"]="ON/OFF"
                        else:
                            if int(number[0])==0:
                                s["text"]="OFF"
                            if int(number[0])==1:
                                s["text"]="ON"
                            if int(number[0])>1:
                                print "ERROR"
                    else:
                        if k[0]==1:
                            self.ccensor.delete(0,END)
                            self.ccensor.insert(0,temp[1:-1])
                        if k[0]==4:
                            self.csmooth.delete(0,END)
                            self.csmooth.insert(0,temp[1:-1])
                        if k[0]==5:
                            self.cdetrend.delete(0,END)
                            self.cdetrend.insert(0,temp[1:-1])    
                        s.set(temp[1:-1])   
        for p in pl:
            s = l[p]
            if p!=4 and p!=5 and p!=1:
                s["text"] = "OFF"
            else:
                s.set("0")
            
    def save_pipeline_file(self,pipeline):
        steps = ["MOTCOR","CENSOR","RETROICOR","TIMECOR","SMOOTH","DETREND","TASK","MOTREG","GSPC1","PHYPLUS","LOWPASS","CUSTOMREG"]
        l = (self.cmotcor, self.str_ccensor, self.cretroicor, self.ctimecor, self.str_csmooth, self.str_cdetrend, self.ctask,self.cmotreg, self.cgspc1, self.cphyplus, self.clowpass, self.ccustomreg)
        F2 = open(pipeline,"w")
        for i in range(0,len(steps)):
            s = l[i]
            if i!=4 and i!=5 and i!=1:
                if s["text"]=="ON":
                    tx = "[1]"
                if s["text"]=="OFF":
                    tx = "[0]"
                if s["text"]=="ON/OFF":
                    tx = "[0,1]"
            else:
                tx=s.get();
            if tx[0]!="[":
                tx = "["+tx
            if tx[-1]!="]":
                tx = tx+"]"
        
            F2.write(steps[i]+"="+tx+"\n")
        F2.close()

     

    def read_input_file(self,inputfile):
 
        
        with open(inputfile) as f:
            lines = f.readlines()
            
            for line in lines:
                current_subject = len(self.list_entries.sin)
                self.list_entries.append(conditions={})
                line=line.rstrip()
                uline = line.split(' ')
                for temp in uline:
                    if ("STRUCT=" in temp.upper()):
                        self.list_entries.structrual[current_subject]=temp[7:]
                    if ("OUT=" in temp.upper()):
                        self.list_entries.out[current_subject]=os.path.split(temp[4:])[0]
                        self.list_entries.prefix[current_subject]=os.path.split(temp[4:])[1]
                    if ("IN=" in  temp.upper()):
                        self.list_entries.sin[current_subject]=temp[3:]
                        xtemp = os.path.split(temp[3:])[0]
                        xtemp = os.path.split(xtemp)
                        self.listbox.insert(END,("%03d" % (current_subject+1))+"-"+xtemp[1])
                    if ("PHYSIO=" in temp.upper()):
                        self.list_entries.physio[current_subject]=temp[7:]
                    if ("TASK=" in temp.upper()):
                        self.list_entries.taskinfo[current_subject]=temp[5:]
                        self.load_task_info(current_subject)
                    if ("CUSTOMREG=" in temp.upper()):
                        self.list_entries.customreg[current_subject]=temp[10:]
                    if ("DROP=" in temp.upper()):
                        self.list_entries.drop[current_subject]=temp[6:-1]
                    

        list_cond = []
        for i in range(0,len(self.list_entries.sin)):
            for j in self.list_entries.conditions[i].keys():
                list_cond.append(j)
        list_cond =list(set(list_cond))
        for i in range(0,len(self.list_entries.sin)):
            for lc in list_cond:
                if not lc in self.list_entries.conditions[i].keys():
                    self.list_entries.conditions[i][lc] = ['','']
                
            
                    
                
        self.listbox_callback([])
        
    def save_input_file(self,inputfile):
        F2 = open(inputfile,"w")
        num_items = len(self.list_entries.sin)
        for i in range(0,num_items):
            if not self.list_entries.customreg[i]:
                cs=""
            else:
                cs = " CUSTOMREG="+self.list_entries.customreg[i]
            s = "IN="+self.list_entries.sin[i]+" OUT="+os.path.join(self.list_entries.out[i],self.list_entries.prefix[i])+" STRUCT="+self.list_entries.structrual[i]+" PHYSIO="+self.list_entries.physio[i]+" TASK="+self.list_entries.taskinfo[i]+" DROP=["+self.list_entries.drop[i]+"]"+cs+"\n"
            F2.write(s)
            self.save_task_info(i)
        F2.close()


    def load_task_info(self,number):

        inputfile = self.list_entries.taskinfo[number]

        
        if not inputfile:
            return
        tx = ""
        c_flag=False
        if not os.path.isfile(inputfile):
            return

        (filename_temp,ext) = os.path.splitext(inputfile)
        if ext==".mat":
            d = loadmat(inputfile)
            '''
            tr_msec = d["split_info"]["TR_MSEC"][0][0][0][0].tolist()
            self.list_entries.tr_msec[number] = tr_msec
            self.list_entries.xtype[number] = d["split_info"]["type"][0][0].tolist()[0]
            if "seed" in d["split_info"].keys():
                self.list_entries.xtype[number] = d["split_info"]["seed"][0][0].tolist()[0]
            
            numconditions = len(d["split_info"]["cond"][0][0]["name"][0])
            for cond_no in range(0,numconditions):
                d["split_info"]["Type"][0][0][0][0].tolist()
            
            '''
            tkMessageBox.showerror("Error","The .mat Task Info is not supported!")
            return        

        with open(inputfile) as f:
            lines = f.readlines()
            for line in lines:
                line = line.rstrip()
                if not line and c_flag:
                    tx=tx+"\n"
                uline = line.split(' ')
                for temp in uline:
                    if ("UNIT=" in temp.upper()):
                        self.list_entries.unit[number]=temp[6:-1]
                        continue
                    if ("TR_MSEC=" in temp.upper()):
                        self.list_entries.tr_msec[number]=temp[9:-1]
                        continue
                    if ("TYPE=" in  temp.upper()):
                        self.list_entries.xtype[number]=temp[6:-1].lower()
                        continue
                    if ("SEED=" in  temp.upper()):
                        self.list_entries.seed[number]=temp[6:]
                        continue
                    if ("NAME=" in  temp.upper()):
                        self.list_entries.conditions[number][temp[6:-1]] = ['','']
                        current_name = temp[6:-1]
                    if ("ONSETS=" in  temp.upper()):
                        self.list_entries.conditions[number][current_name][0] = temp[8:-1]
                    if ("DURATION=" in  temp.upper()):
                        self.list_entries.conditions[number][current_name][1] = temp[10:-1]
                    if not temp:
                        continue
                    c_flag = True
                    tx = tx+temp+"\n"
                    
        self.list_entries.tasktext[number] = tx

    def save_task_info(self,number):
        outputfile = self.list_entries.taskinfo[number]
        if not outputfile:
            return
        c_flag=False
        F2 = open(outputfile,"w")
        if self.list_entries.unit[number]:
            F2.write("UNIT=[%s]\n" % self.list_entries.unit[number])
        if self.list_entries.tr_msec[number]:
            F2.write("TR_MSEC=[%s]\n" % self.list_entries.tr_msec[number])
        if self.list_entries.xtype[number]:
            F2.write("TYPE=[%s]\n" % self.list_entries.xtype[number])
        if self.list_entries.seed[number]:
            F2.write("SEED=[%s]\n" % self.list_entries.seed[number])
        
        F2.write("\n")
        for name in self.list_entries.conditions[number].keys():
            F2.write("NAME=[%s]\n" % name)
            F2.write("ONSETS=[%s]\n" % self.list_entries.conditions[number][name][0].replace(" ",","))
            F2.write("DURATION=[%s]\n" % self.list_entries.conditions[number][name][1].replace(" ",","))             
            F2.write("\n")
                       
        F2.close()

    def read_session_file(self,inputfile):


        d = loadmat(inputfile,struct_as_record=True)
        num_runs = d["session_info"]["num_runs"].tolist()[0][0].tolist()[0][0]
        num_conds = d["session_info"]["num_conditions"].tolist()[0][0].tolist()[0][0]
        cond_list = []
        for cond_no in range(0,num_conds):
            cond_list.append(d["session_info"]["condition"].tolist()[0][0].tolist()[0][cond_no].tolist()[0])
        for run_no in range(0,num_runs):
            data_name= d["session_info"]["run"][0][0]["data_files"][0][run_no].tolist()[0][0].tolist()[0]
            data_path = d["session_info"]["run"][0][0]["data_path"][0][run_no].tolist()[0]
            inp = os.path.join(data_path,data_name)
            pref = "preprocessed_"+data_name
            dr = "0,0"
            un = "TR"
            conds = {}
            for cond_no in range(0,num_conds):
                blk_onsets    = convert_list(d["session_info"]["run"][0][0]["blk_onsets"][0][run_no][0][cond_no].tolist(),1)
                blk_length    = convert_list(d["session_info"]["run"][0][0]["blk_length"][0][run_no][0][cond_no].tolist(),0)
                conds[cond_list[cond_no]] = [blk_onsets,blk_length]
                
            self.list_entries.append(sin=inp,prefix=pref,drop=dr,unit=un,conditions=conds)
            self.listbox.insert(END, "%03d-%s~%s/run%d" % (len(self.list_entries.sin),data_name[0:1],data_name[-8:-4],run_no+1))

        list_cond = []
        for i in range(0,len(self.list_entries.sin)):
            for j in self.list_entries.conditions[i].keys():
                list_cond.append(j)
        list_cond =list(set(list_cond))
        for i in range(0,len(self.list_entries.sin)):
            for lc in list_cond:
                if not lc in self.list_entries.conditions[i].keys():
                    self.list_entries.conditions[i][lc] = ['','']
                
    def LoadPLSSessionFiles(self):
        y = askopenfilenames(filetypes=[("Select PLS sessionfile","*.mat")])
        if y:
            for i in range(0, len(y)):
                self.read_session_file(y[i])

        self.move_cursor(END)     
        self.update_entries()
                       
    def Menus(self):
        root = self.parent
        menubar = Menu(root)
        filemenu = Menu(menubar, tearoff=0)
        filemenu.add_command(label="Clear Input List",command=self.NewInputs) #command=self.donothing)
        filemenu.add_command(label="Load Input Files From Text File",command=self.BrowseInputFile) # command=self.donothing)
        filemenu.add_command(label="Save Input Files As Text File",command=self.SaveInputFile)#, command=self.donothing)
            

        menubar.add_cascade(label="1: Input Files", menu=filemenu)

        pipelinemenu = Menu(menubar, tearoff=0)
        pipelinemenu.add_command(label="Load Default Pipeline list",command=self.LoadDefaultPipelineList) #command=self.donothing)
        pipelinemenu.add_command(label="Load Pipeline list From Text File", command=self.BrowsePipelineList) # command=self.donothing)
        pipelinemenu.add_command(label="Save Pipeline list As Text File", command=self.SavePipelineList)#, command=self.donothing)
  
        menubar.add_cascade(label="2: Pipeline list", menu=pipelinemenu)

        analysismenu = Menu(menubar, tearoff=0)
        menubar.add_cascade(label="3: Select Analysis Model", menu=analysismenu)

        
        runmenu = Menu(menubar, tearoff=0)
        runmenu.add_command(label="Advanced settings", command=self.AdvanceSettings)
        runmenu.add_command(label="Submit", command=self.SubmitCallback)
        menubar.add_cascade(label="4: Run", menu=runmenu)

        exitmenu = Menu(menubar, tearoff=0)
        exitmenu.add_command(label="Exit",command=root.destroy)
        menubar.add_cascade(label="Exit", menu=exitmenu)


        helpmenu = Menu(menubar, tearoff=0)
        helpmenu.add_command(label="About...",command=self.AboutMe)

        menubar.add_cascade(label="Help",menu=helpmenu)

        root.config(menu=menubar)

    def donothin(self):
        debug = 1

    def DeleteCondition(self):
        name = self.str_cname.get()
        if not name:
            return
        t = Toplevel(self)
        t.geometry('{0}x{1}'.format(300, 150))
        t.wm_title("Delete")
        l = Label(t, text="Condition: %s" % self.str_cname.get())
        l.place(relx=0.053,rely=0.135, relheight=0.2, relwidth=0.52)
        l = Button(t, text="Delete from all entries",  command=self.DeleteFromAllEntries)
        l.place(relx=0.053,rely=0.335, relheight=0.2, relwidth=0.9)
        l = Button(t, text="Cancel",command=t.destroy)
        l.place(relx=0.053,rely=0.635, relheight=0.2, relwidth=0.9)
        t.grab_set()
        self.current_window = t
        
    def DeleteFromAllEntries(self):
        name = self.str_cname.get()
        for i in range(0,len(self.list_entries.conditions)):
            if name in self.list_entries.conditions[i].keys():
                del self.list_entries.conditions[i][name]

        index = map(int, self.listbox.curselection())
        index = index[0]       
        self.cconditions['menu'].delete(0, 'end')
        for i in self.list_entries.conditions[index].keys():
            self.cconditions['menu'].add_command(label=i, command=lambda i=i:self.str_cname.set(i))
        cnd = self.list_entries.conditions[index].keys()
        if cnd:
            self.str_cname.set(cnd[0])
        else:
            self.str_cname.set('')
        self.current_window.destroy()

    def RenameCondition(self):
        name = self.str_cname.get()
        if not name:
            return
        t = Toplevel(self)
        t.geometry('{0}x{1}'.format(300, 150))
        t.grab_set()
        t.wm_title("Rename")
        l = Label(t, text="New Condition Name:")
        l.place(relx=0.053,rely=0.135, relheight=0.2, relwidth=0.52)
        l = Entry(t, textvariable=self.str_addcondition)
        l.place(relx=0.053,rely=0.335, relheight=0.2, relwidth=0.9)
        l = Button(t, text="Rename in all entries",command=self.RenameConditionInAllEntries)
        l.place(relx=0.053,rely=0.635, relheight=0.2, relwidth=0.9)
        self.current_window = t    
    
    def RenameConditionInAllEntries(self):
        newname = self.str_addcondition.get()
        name = self.str_cname.get()
        for i in range(0,len(self.list_entries.conditions)):
            if name in self.list_entries.conditions[i].keys():
                temp = self.list_entries.conditions[i][name]
                del self.list_entries.conditions[i][name]
                self.list_entries.conditions[i][newname] = temp

        index = map(int, self.listbox.curselection())
        index = index[0]       
        self.cconditions['menu'].delete(0, 'end')
        for i in self.list_entries.conditions[index].keys():
            self.cconditions['menu'].add_command(label=i, command=lambda i=i:self.str_cname.set(i))
        cnd = newname
        if cnd:
            self.str_cname.set(cnd)
        else:
            self.str_cname.set('')
        self.current_window.destroy()
        self.str_addcondition.set("")


    def AddCondition(self):
        t = Toplevel(self)
        t.geometry('{0}x{1}'.format(300, 150))
        t.grab_set()
        t.wm_title("Add")
        l = Label(t, text="Condition Name:")
        l.place(relx=0.053,rely=0.135, relheight=0.2, relwidth=0.52)
        l = Entry(t, textvariable=self.str_addcondition)
        l.place(relx=0.053,rely=0.335, relheight=0.2, relwidth=0.9)
        l = Button(t, text="Add to all entries",command=self.AddtoAllEntries)
        l.place(relx=0.053,rely=0.635, relheight=0.2, relwidth=0.9)
        self.current_window = t
    def AddtoAllEntries(self):

        self.current_window.destroy()
        if not self.str_addcondition.get() in self.list_entries.conditions[0].keys():
            for i in range(0,len(self.list_entries.sin)):
                self.list_entries.conditions[i][self.str_addcondition.get()]=['','']
                
        self.cconditions['menu'].delete(0, 'end')
        for i in self.list_entries.conditions[0].keys():
            self.cconditions['menu'].add_command(label=i, command=lambda i=i:self.str_cname.set(i))
        self.str_cname.set(self.str_addcondition.get())
        self.str_addcondition.set("")
        self.current_window.destroy()

    def AdvanceSettings(self):
        t = Toplevel(self)
        t.geometry('{0}x{1}'.format(300, 150))
        t.wm_title("Advance Settings")
        Label(t, text="NUMCORES", anchor=E).grid(row=0,column=1)
        Label(t, text="QUEUE",anchor=E).grid(row=1,column=1)
        Label(t, text="ENVIRONMENT",anchor=E).grid(row=2,column=1)
        Label(t, text="AUTODETECT",anchor=E).grid(row=3,column=1)
        Label(t, text="N_resample", anchor=E).grid(row=4,column=1)
        Label(t, text="SUBMISSION NODE",anchor=E).grid(row=5,column=1)
        
        Entry(t,textvariable=self.str_numcores).grid(row=0, column=2)
        Entry(t,textvariable=self.str_queue).grid(row=1, column=2)
        OptionMenu(t,self.str_environment,'matlab','octave').grid(row=2, column=2)
        Checkbutton(t,textvariable='', variable=self.int_autodetect).grid(row=3, column=2)
        Entry(t,textvariable=self.str_nresample).grid(row=4, column=2)        
        Entry(t,textvariable=self.str_submissionnode).grid(row=5, column=2)
        t.grab_set()

    def SubmitCallback(self):
        if len(self.list_entries.sin)==0:
            tkMessageBox.showerror("Error","Specify the dataset first! There is nothing to submit...")
            return
        y = asksaveasfilename(filetypes=[("Input File","*.txt")])
        if not y:
            tkMessageBox.showerror("Error","Input File have to be saved before submission! \n Submission Cancelled!")
            return
        (dirname,filename)= os.path.split(y)
        (fname,ext) = os.path.splitext(filename)
        dirname = os.path.dirname(y) + "/%s_taskinfo/" % fname
        inputfile_name = os.path.dirname(y) +'/'+ fname + '.txt'
        
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        for i in range(0,len(self.list_entries.sin)):
        #    if not self.list_entries.taskinfo[i]:
            self.list_entries.taskinfo[i] = dirname+"%03d.txt" %i
            self.save_task_info(i)
        print inputfile_name
        self.save_input_file(inputfile_name)
        pipelinelist_name = os.path.dirname(y) + "/%s_PipelineList.txt" % fname

        self.save_pipeline_file(pipelinelist_name)

        submissionnode= self.str_submissionnode.get()
        if submissionnode:
            submissionnode = "ssh "+ submissionnode + " \""
        
        codepathname  = os.path.dirname(sys.argv[0])
        codefull_path = os.path.abspath(codepathname)

        ptemp = self.str_partmenu.get()
        if "all parts" in ptemp: part = 0
        if "2" in ptemp: part=2
        if "1" in ptemp: part=1
        if "3" in ptemp: part=3

        inputfile = inputfile_name
        pipeline  = pipelinelist_name
        metric    = self.str_metricmenu.get()
        contrast  = self.str_contrast.get()
        numcores  = self.str_numcores.get()
        queue     = self.str_queue.get()
        keepmean  = int(self.str_checkbox_keepmean.get())
        voxelsize = self.str_voxelsize.get()
        env       = self.str_environment.get()
        analysis_model = self.str_modelmenu.get()
        ref = self.str_referene.get()
        if ref:
            ref = "-r "+ref
        if voxelsize:
            voxelsize = "-v '"+voxelsize+"'"
        if contrast:
            contrast = "--contrast "+contrast
        if analysis_model=="None":
            contrast = "--contrast nocontrast"
        convolve = int(self.str_checkbox_convolve.get())
        
        codefull_path2 = os.path.dirname(codefull_path)
        cmd = submissionnode + codefull_path+"/../Run_Pipelines.py -p {0} -i {1} -c {2} -a {3} {4} --convolve {5} -m {6} {7} -n {8} -q {9} -k {10} {11} -e {12}".format(part,inputfile,pipeline,analysis_model,ref,convolve,metric,contrast,numcores,queue,keepmean,voxelsize,env)
        
        if int(self.str_checkbox_deoblique.get())==1:
            cmd = cmd + " --DEOBLIQUE"
        if self.str_tpatternmenu.get()!="None":
            cmd = cmd + " --TPATTERN" + str_tpatternmenu.get()
            
        parameter_list = self.parameter_dict[analysis_model]
        for name in parameter_list:
            val = self.parameters_rf[name].get()
            cmd = cmd + " --"+name+" "+val
        if self.int_autodetect.get():
            cmd = cmd + " --autodetect"
            cmd = cmd + " --N_resample "+self.str_nresample.get()

        if submissionnode:
            cmd = cmd + "\""
        print cmd
        os.system(cmd)
        
    def set_tooltip_messages(self):

        
        CreateToolTip(self.lin,"Name and location of unprocessed fMRI data you wish to optimize")
        CreateToolTip(self.lout,"Location for final processed & optimized outputs")
        CreateToolTip(self.lprefix,"Name for final processed & optimized outputs")
        CreateToolTip(self.lphysio,"PHYSIO={(path)/(physiological data prefix)} \n Name and location of Cardiac and Respiratory data, converted into .1D design files.\n must be formatted as (path)/(physiological data prefix).puls.1D and \n (path)/(physiological data prefix).resp.1D, respectively")
        CreateToolTip(self.lcustomreg,"CUSTOMREG={(path)/(nifti file)}\nName and location of binary mask, where 1=brain regions containing noise confounds (e.g. white matter). \n The mean time-series of this ROI is computed and regressed out of all voxels in the brain by the CUSTOMREG pipeline step.\n This volume must have the same dimensions as the user's STRUCT file. \n If you want to perform CUSTOMREG, you will need to create the binary ROI mask yourself, \n e.g. using one of the neuroimaging packages that allows you to manually edit images \n (e.g. fslview, mricron, afni). For advanced users only!")
        CreateToolTip(self.lstructrual,"STRUCT={(path)/(T1 anatomical nifti file)} \n Name and location of subject's 3D structural brain image, which is used to \n data to a common MRI template space, if desired. \n * STRUCT field can be omitted if you do not warp subject data to common template")
        CreateToolTip(self.ltaskinfo,"{(path)/(task information)}\n Name and location of formatted text file, describing the experimental paradigm for the selected dataset in \n Input Files list box (left), including task onsets and duration. \n You can create this file yourself (see the manual) or \n you can create/edit the task file or define/remove different conditions and their onsets, durations below.")
        CreateToolTip(self.ldrop,"DROP={# scans to drop at start},{# scans to drop at end}\nDiscards non-equilibrium and instructional scan volumes at the ends of the run\nDROP=0,0 does not discard any volumes")

        CreateToolTip(self.lunit,"Unit in which task onset/duration is measured \n TR, sec or msec(#scan volumes, seconds or milliseconds)")
        CreateToolTip(self.ltr_msec,"Time between scan volumes (repetition time) \n Integer value, in milliseconds")
        CreateToolTip(self.lconditions,"Name of this task condition (its onsets and duration are listed below) \n Use any name you like, but avoid special characters e.g. ""+-[]"". \n The name baseline is reserved for the scans in which the brain is in the ""REST"" condition (See manual for details)")
        CreateToolTip(self.lonset,"List of onset times, in time-units specified by UNIT field \n Must be a comma-separated list of numbers, e.g. ONSETS=1,10,30")
        CreateToolTip(self.lduration,"List of condition durations corresponding to the above onsets, in time-units specified by UNIT field \n Must be a comma-separated list of numbers, e.g. Duration=9,10,8");
        CreateToolTip(self.lseed,"Name of ROI brain mask, for seed-based connectivity analysis.\n*only required for SCONN analysis model")
        CreateToolTip(self.ltype,"Type of task paradigm being processed/analyzed \n block, event or nocontrast \n all the dataset added in the list below must have a similar type")



        CreateToolTip(self.lmotcor,"Rigid-body motion correction")
        CreateToolTip(self.lcensor,"Censoring outlier spikes in data \n Comma-separated list of 0,1,2,3 options \n (0=OFF, 1=BASIC, 2=AGGRESSIVE PCA, 3=AGGRESSIVE ICA)")
        CreateToolTip(self.lretroicor,"Physiological noise correction #1")
        CreateToolTip(self.ltimecor,"Slice-timing correction")
        CreateToolTip(self.lsmooth,"Gaussian spatial smoothing \n Comma-separated list of FWHM smoothing kernels (mm), e.g. 0,6,8")
        CreateToolTip(self.ldetrend,"Temporal detrending \n Comma-separated list of Legendre polynomial detrending order, e.g. 1,2,3")
        CreateToolTip(self.lmotreg,"Motion parameter regression")
        CreateToolTip(self.lgspc1,"Global signal corrected with PC#1")
        CreateToolTip(self.lphyplus,"Physiological noise correction #2")
        CreateToolTip(self.llowpass,"Regression of noise ROI")
        CreateToolTip(self.ltask,"Regression of noise ROI")
        CreateToolTip(self.lcustomreg2,"Type of task paradigm being processed/analyzed \n block, event or nocontrast")


        CreateToolTip(self.lanalysismodel,"Available analysis models in the PRONTO framework \n GNB (Gaussian Naïve Bayes) predictive general linear model of difference between 2 conditions \n LDA (Linear Discriminant Analysis)predictive multivariate model of differences between 2 conditions \n erGNB (Event-Related GNB)predictive general linear model, treating each time-point in HRF as a different condition \n erCVA (Event-Related CVA)predictive multivariate model, treating each time-point in an HRF window as a different condition \n SCONN (Seed-based Connectivity)measures pairwise correlations of all voxels, relative to average signal in a seed Region of Interest (ROI) \n gPCA (Generalized PCA)Principal Component Analysis decomposition of fMRI dataset, identifying most consistent subspace \n\n GLM and erGLM The PRONTO code also has standard General Linear Model (GLM) and event-related GLM (erGLM) models available, \n but we do not recommend them: they can only be used to measure R (Reproducibility), whereas P (Prediction) accuracy cannot be computed. \n If you want to analyze your data assuming independent voxels, as in GLM, \n we recommend using GNB or erGNB, which are the predictive versions of these models \n\n None is used when you do not want to optimize preprocessing pipeline, the code only generates the preprocessed datasets with preprocessing steps chosen in 2: Pipeline list")
        CreateToolTip(self.lmodelparameter,"Optional arguments that are specific to individual analysis models (See the manual) \n The options depend on the Analysis Model")
        CreateToolTip(self.lmetric,"Metric used to choose optimal pipelines standard options \n includes R= reproducibility P= prediction dPR=combined prediction & reproducibility (this is the default) \n The options depend on the Analysis Model")
        CreateToolTip(self.lcontrast,"Task contrast being analyzed and optimized for. \n Necessary when more than two conditions are defined in the task files. \n Syntax: ""CON1-CON2,CON2-CON3"", where CON1, CON2, CON3 are condition names. \n Condition names are listed in Condition(s) Option Menu")
        CreateToolTip(self.lvoxelsize,"Determines the output voxel size of nifti file. Syntax: 3.0 3.0 5.0 gives 3x3x5 mm voxels as final results. \n Default is to keep output voxels the same size as input.")
        CreateToolTip(self.keepmean_checkbox,"Determines whether the output nifti files have temporal means re-added to each voxel after processing \n (default removes the mean scan)")
        CreateToolTip(self.deoblique_checkbox,"Corrects for raw data that are at an oblique angle relative to the cardinal scanning axes, \n using AFNI's 3dWarp program. May improve the quality of registration")
        CreateToolTip(self.convolv_checkbox,"determines whether user-provided design matrix should be convolved with the canonical HRF modeled by two gamma functions \n (AFNI's SPM1 function), only for TASK pipeline step and GLM Analysis Model")
        CreateToolTip(self.ltpattern,"Defines axial slices acquisition order for slice-timing correction, \n if this information is not available in the header. Options include: \n altplus (alternating, + direction - standard interleaved protocol) \n altminus (alternating, negative direction) \n seqplus (sequential, positive direction) \n seqminus (sequential, negative direction)")


        CreateToolTip(self.lselectpipelinepart,"Only do specific parts of full optimization pipeline. \n 1= run all pipelines and produce metrics \n 2= select optimal pipelines based on metrics (must have already run part-1) \n 3= do spatial normalization of optimized results (must have already run parts-1,2) \n Default= do all three steps")
        CreateToolTip(self.lselectreferecevolum,"Path and name of reference volume used in (optional) spatial normalization step")

        CreateToolTip(self.lstep1,"STEP 1: You should specify the unproprocessed data: \n You can add/remove/copy entries in the list box below (Right Click)\n Or you can load the input file from the above menu 1:Input Files \n For each entries you should set/edit the input files and output paths (See Left: Input, Output, PHYSIO, ...) \n")       
        CreateToolTip(self.lstep2,"STEP 2: Here You should specify what preprocessing steps you want to run on your datasets \n You can load/save preprocessing step using the above menu 2:Pipeline list ")
        
        CreateToolTip(self.lstep3,"STEP 3: In this step you have to specify the Analysis Model, and Contrast as two important parameters \n Do not change the other paramters if you are not sure what they are.")
        CreateToolTip(self.lstep4,"STEP 4: This is last step of PRONTO preprocessing! \n You can set some advanced settings in the above menu 4:Run")
       
      
    def AboutMe(self):
        tkMessageBox.showinfo("About Pronto", "This GUI runs PRONTO (PReprocessing OptimizatioN TOolkit),which does fast optimization of preprocessing pipelines for BOLD fMRI (Blood Oxygenation Level Dependent functional MRI).This software package identifies the set of preprocessing steps ('pipeline') specific to each dataset, which optimizes quality metrics of Prediction and Reproducibility for post-processing analysis results (Strother et al., 2002). This procedure has been shown to significantly improve signal detection, reliability of brain activations, and sensitivity to brain-behaviour correlations (Churchill et al., 2012a, 2012b). The pipeline software can also be used for simple automated batch-processing of fMRI datasets, if no appropriate analysis model is currently available to do optimization (e.g. some resting-state connectivity studies). \n\n\n PRONTO GUI designed by Babak Afshin-Pour (bafshinpour@research.baycrest.org), 2015 \n Rotman Research Institue, Baycrest")
        

if __name__ == "__main__": 
    
    root = Tk()
    root.geometry('{0}x{1}'.format(1250, 400))
    root.title("PRONTO (Preprocessing optimization toolkit)")
    Files(master=root)
  #  main = Files(root)
    root.mainloop()
    
    #MainFrame().mainloop()
    
