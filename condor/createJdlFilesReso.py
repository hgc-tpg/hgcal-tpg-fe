import itertools
import os
import sys
import subprocess
import time

#IMPORT MODULES FROM OTHER DIR

iloop= "86"

#samplelist_Ideal = ["SingleEle_Ideal_PU0","SinglePi_Ideal_PU0"]
#samplelist_Ideal = ["SingleEle_Ideal_PU0"]
#samplelist_Ideal = ["SinglePi_Ideal_PU0"]
#ntuple_path_ideal = ["flatpt_10K"]
samplelist_Ideal = ["SinglePi_Ideal_PU0_v17_Thresh0"]

ntuple_path_ideal = ["ntuples"]

samplelist_PU0 = ["SinglePi_realistic_PU0", "SingleEle_realistic_PU0", "VBFHToInvisible_realistic_PU0", "MinBias_realistic_PU0"]

#samplelist_PU0_Emyr = ["singlePion_PU0_Ideal"]
#samplelist_PU0_Emyr = ["doublePhoton_PU0", "singleElectron_PU0_Ideal"]
#samplelist_PU0_Emyr = ["singlePion_PU0", "singlePion_PU0", "doublePhoton_PU0", "vbfHInv_0PU"] #"doubleElectron_PU200", "singlePion_PU200", "singlePion_PU0_Ideal", "singleElectron_PU0_Ideal", 
#samplelist_PU0_Emyr = ["doubleElectron_PU200", "singlePion_PU200", "singlePion_PU0_Ideal", "singleElectron_PU0_Ideal", "vbfHInv_0PU", "vbfHInv_200PU"]
samplelist_PU0_Emyr = ["vbfHInv_0PU_ID"]
#samplelist_PU0_Emyr = ["singlePion_PU0_Realistic"]
#samplelist_PU0_Emyr = ["vbfHInv_0PU", "vbfHInv_0PU_ID", "singlePion_PU0", "singlePion_PU0_Ideal", "singlePion_PU0_Realistic"] #"vbfHInv_0PU", "vbfHInv_0PU_ID", 
#"vbfHInv_0PU", "vbfHInv_200PU" "doubleElectron_PU200", "singlePion_PU200", "singlePion_PU0_Ideal", "singleElectron_PU0_Ideal"
#samplelist_PU0_Emyr = ["vbfHInv_0PU", "doublePhoton_PU0", "singlePion_PU0", "singleElectron_PU0_Ideal", "singlePion_PU0_Ideal", "singlePion_PU0_Realistic"]
#samplelist_PU0_Emyr = ["doublePhoton_PU0", "singleElectron_PU0_Ideal", "vbfHInv_0PU", "singlePion_PU0_Ideal", "singlePion_PU0_Realistic"]

samplelist_PU200 = ["SinglePi_realistic_PU200", "SingleEle_realistic_PU200", "VBFHToInvisible_realistic_PU200", "MinBias_realistic_PU140"]
#samplelist_PU200_Emyr = ["doubleElectron_PU200", "singlePion_PU200", "vbfHInv_200PU"]
#samplelist_PU200_Emyr = ["doubleElectron_PU200"]
#samplelist_PU200_Emyr = ["singlePion_PU200", "vbfHInv_200PU"]
samplelist_PU200_Emyr = ["vbfHInv_200PU"]

ntuple_path = ["ntuples"]

#triangle_side_list = ["0.016", "0.03", "0.045", "0.06", "0.075", "0.090", "0.105"]
triangle_side_list = ["0.016", "0.03", "0.045"]
#triangle_side_list = ["0.0113", "0.016", "0.0226"]
#triangle_side_list = ["0.03", "0.045"]
#triangle_side_list = ["0.0113"]
#triangle_side_list = ["0.016"]

#----------------------------------------
#Create run and log directory
#----------------------------------------
jdlDir = 'tmpLog_s2reso_iter%s'%iloop
if not os.path.exists("%s/logs"%jdlDir):
    os.makedirs("%s/logs"%jdlDir)
condorLogDir = "logs"
os.system("cp runStage2Reso.sh %s/"%jdlDir)

# lxsetup
common_command = \
'executable = runStage2Reso.sh\n\
+MaxRuntime = 10799\n\
request_memory = 4000\n\
output = %s/log_$(cluster)_$(process).stdout\n\
error  = %s/log_$(cluster)_$(process).stderr\n\
log    = %s/log_$(cluster)_$(process).condor\n\n'%(condorLogDir, condorLogDir, condorLogDir)

# lxplus setup
# common_command = \
# 'Universe   = vanilla\n\
# should_transfer_files = YES\n\
# when_to_transfer_output = ON_EXIT\n\
# Transfer_Input_Files = local.tar.gz, runStage2Reso.sh\n\
# x509userproxy = $ENV(X509_USER_PROXY)\n\
# use_x509userproxy = true\n\
# +BenchmarkJob = True\n\
# +MaxRuntime = 10799\n\
# notification = never\n\
# max_transfer_input_mb = 4096\n\
# max_transfer_output_mb = 4096\n\
# request_disk = 4000000\n\
# Output = %s/log_$(cluster)_$(process).stdout\n\
# Error  = %s/log_$(cluster)_$(process).stderr\n\
# Log    = %s/log_$(cluster)_$(process).condor\n\n'%(condorLogDir, condorLogDir, condorLogDir)

#----------------------------------------
#Create tar file
#----------------------------------------
tarFile = "local.tar.gz"
if os.path.exists(tarFile):
    os.system("rm %s"%tarFile)
os.system("tar --exclude-from='../.tarexclude' --exclude-backups -zcvf %s ../../hgcal-tpg-fe "%tarFile)
os.system("cp local.tar.gz %s/"%jdlDir)

year=2026
ofindex=0
infile='/vols/cms/idas/PiPlusMinus/PU0/Pion_test2/SingleEle/pt100GeV_10K/ntuple_Thresh.root'
nevents=1000
sidelength=0.016
ofextn='SinglePi_Ideal'
sampletype='Result_iter%s/SinglePi_Ideal_PU0'%iloop
# iloop=$7
# clusproc=$8

# ----------------------------------------
# Create jdl files
# ----------------------------------------
subFile = open('%s/condorSubmit.sh'%jdlDir,'w')

# jdlName = 'submitJobs_%s.jdl'%(year)
# jdlFile = open('%s/%s'%(jdlDir,jdlName),'w')
# jdlFile.write(common_command)
# #nqueue = subprocess.Popen(['grep -i queue local/%s/%s/Cat1_Inc/Mass%s/condor_run.sub | awk \'{print $2}\''%(year,channel,mass)],shell=True,stdout=subprocess.PIPE).communicate()[0].decode("utf-8")
# run_command =  'Arguments  = %s %s %s %s %s %s %s $(process) \nQueue 1\n\n' %(infile,ofindex,nevents,sidelength,ofextn,sampletype,iloop)
# jdlFile.write(run_command)        
# jdlFile.close()
# subFile.write("condor_submit %s\n"%jdlName)

nevents=0

# for sample in samplelist_Ideal:
#     for ntuplepath in ntuple_path_ideal:
#         for sidelength in triangle_side_list:
#             dirpath = '%s/stage2_emulator_tests/Indra/%s/%s/'%(os.environ["HOME"],sample,ntuplepath)
#             findarg = '%s -name \"ntuple_step3_*.root\"'%(dirpath)
#             p = subprocess.Popen([("find %s "%findarg)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)
#             #print ("findarg: %s"%(findarg))
#             filelist, errors = p.communicate()
#             filelist = filelist.split('\n')
#             filelist.remove("")
#             nFiles = len(filelist)
#             findex  = 0
#             fsidelen = float(sidelength)
#             fsidelen_index = int(fsidelen*1000)
#             jdlName = 'submitJobs_%s_%s_%s.jdl'%(sample,ntuplepath,fsidelen_index)
#             jdlFile = open('%s/%s'%(jdlDir,jdlName),'w')
#             jdlFile.write(common_command)
#             for fname in filelist:
#                 print ("fname: %s, sidelength: %s, index: %s"%(fname,sidelength,findex))
#                 ofextn = '%s_%s'%(ntuplepath,fsidelen_index)
#                 run_command =  'Arguments  = %s %s %s %s %s %s %s $(process) \nQueue 1\n\n' %(fname,findex,nevents,sidelength,ofextn,sample,iloop)
#                 jdlFile.write(run_command)        
#                 findex = findex + 1
#             jdlFile.close()
#             subFile.write("condor_submit %s\n"%jdlName)
            
# for sample in samplelist_Ideal:
#     for ntuplepath in ntuple_path_ideal:
#         for sidelength in triangle_side_list:
#             dirpath = '%s/stage2_emulator_tests/%s/%s/'%(os.environ["HOME"],sample,ntuplepath)
#             findarg = '%s -name \"ntuple_Thresh_cmstwk.root\"'%(dirpath)
#             p = subprocess.Popen([("find %s "%findarg)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)
#             filelist, errors = p.communicate()
#             filelist = filelist.split('\n')
#             filelist.remove("")
#             nFiles = len(filelist)
#             findex  = 0
#             for fname in filelist:
#                 print ("fname: %s, sidelength: %s, index: %s"%(fname,sidelength,findex))
#                 fsidelen = float(sidelength)
#                 fsidelen_index = int(fsidelen*1000)
#                 ofextn = '%s_%s'%(ntuplepath,fsidelen_index)
#                 jdlName = 'submitJobs_%s_%s_%s_%s.jdl'%(sample,ntuplepath,fsidelen_index,findex)
#                 jdlFile = open('%s/%s'%(jdlDir,jdlName),'w')
#                 jdlFile.write(common_command)
#                 run_command =  'Arguments  = %s %s %s %s %s %s %s $(process) \nQueue 1\n\n' %(fname,findex,nevents,sidelength,ofextn,sample,iloop)
#                 jdlFile.write(run_command)        
#                 jdlFile.close()
#                 subFile.write("condor_submit %s\n"%jdlName)
#                 findex = findex + 1
                
# for sample in samplelist_PU0:
#     for ntuplepath in ntuple_path:
#         for sidelength in triangle_side_list:
#             dirpath = '%s/stage2_emulator_tests/%s/%s/'%(os.environ["HOME"],sample,ntuplepath)
#             findarg = '%s -name \"*.root\"'%(dirpath)
#             p = subprocess.Popen([("find %s "%findarg)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)
#             filelist, errors = p.communicate()
#             filelist = filelist.split('\n')
#             filelist.remove("")
#             nFiles = len(filelist)
#             findex  = 0
#             for fname in filelist:
#                 print ("fname: %s, sidelength: %s, index: %s"%(fname,sidelength,findex))
#                 fsidelen = float(sidelength)
#                 fsidelen_index = int(fsidelen*1000)
#                 ofextn = '%s_%s'%(ntuplepath,fsidelen_index)
#                 jdlName = 'submitJobs_%s_%s_%s_%s.jdl'%(sample,ntuplepath,fsidelen_index,findex)
#                 jdlFile = open('%s/%s'%(jdlDir,jdlName),'w')
#                 jdlFile.write(common_command)
#                 run_command =  'Arguments  = %s %s %s %s %s %s %s $(process) \nQueue 1\n\n' %(fname,findex,nevents,sidelength,ofextn,sample,iloop)
#                 jdlFile.write(run_command)        
#                 jdlFile.close()
#                 subFile.write("condor_submit %s\n"%jdlName)
#                 findex = findex + 1

# for sample in samplelist_PU200:
#     for ntuplepath in ntuple_path:
#         for sidelength in triangle_side_list:
#             dirpath = '%s/stage2_emulator_tests/%s/%s/'%(os.environ["HOME"],sample,ntuplepath)
#             findarg = '%s -name \"*.root\"'%(dirpath)
#             p = subprocess.Popen([("find %s "%findarg)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)
#             filelist, errors = p.communicate()
#             filelist = filelist.split('\n')
#             filelist.remove("")
#             nFiles = len(filelist)
#             findex  = 0
#             for fname in filelist:
#                 print ("fname: %s, sidelength: %s, index: %s"%(fname,sidelength,findex))
#                 fsidelen = float(sidelength)
#                 fsidelen_index = int(fsidelen*1000)
#                 ofextn = '%s_%s'%(ntuplepath,fsidelen_index)
#                 jdlName = 'submitJobs_%s_%s_%s_%s.jdl'%(sample,ntuplepath,fsidelen_index,findex)
#                 jdlFile = open('%s/%s'%(jdlDir,jdlName),'w')
#                 jdlFile.write(common_command)
#                 run_command =  'Arguments  = %s %s %s %s %s %s %s $(process) \nQueue 1\n\n' %(fname,findex,nevents,sidelength,ofextn,sample,iloop)
#                 jdlFile.write(run_command)        
#                 jdlFile.close()
#                 subFile.write("condor_submit %s\n"%jdlName)
#                 findex = findex + 1

########## Ntuples from Emyr
# for sample in samplelist_PU0_Emyr:
#     for sidelength in triangle_side_list:
#         dirpath = '%s/stage2_emulator_tests/Emyr/%s/'%(os.environ["HOME"],sample)
#         findarg = '%s -name \"*.root\"'%(dirpath)
#         p = subprocess.Popen([("find %s "%findarg)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)
#         filelist, errors = p.communicate()
#         filelist = filelist.split('\n')
#         filelist.remove("")
#         nFiles = len(filelist)
#         findex  = 0
#         for fname in filelist:
#             print ("fname: %s, sidelength: %s, index: %s"%(fname,sidelength,findex))
#             fsidelen = float(sidelength)
#             fsidelen_index = int(fsidelen*1000)
#             ofextn = 'ntuples_%s'%(fsidelen_index)
#             jdlName = 'submitJobs_%s_ntuples_%s_%s.jdl'%(sample,fsidelen_index,findex)
#             jdlFile = open('%s/%s'%(jdlDir,jdlName),'w')
#             jdlFile.write(common_command)
#             run_command =  'Arguments  = %s %s %s %s %s %s %s $(process) \nQueue 1\n\n' %(fname,findex,nevents,sidelength,ofextn,sample,iloop)
#             jdlFile.write(run_command)        
#             jdlFile.close()
#             subFile.write("condor_submit %s\n"%jdlName)
#             findex = findex + 1
                
# for sample in samplelist_PU200_Emyr:
#     for sidelength in triangle_side_list:
#         dirpath = '%s/stage2_emulator_tests/Emyr/%s/'%(os.environ["HOME"],sample)
#         findarg = '%s -name \"*.root\"'%(dirpath)
#         p = subprocess.Popen([("find %s "%findarg)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)
#         filelist, errors = p.communicate()
#         filelist = filelist.split('\n')
#         filelist.remove("")
#         nFiles = len(filelist)
#         findex  = 0
#         for fname in filelist:
#             print ("fname: %s, sidelength: %s, index: %s"%(fname,sidelength,findex))
#             fsidelen = float(sidelength)
#             fsidelen_index = int(fsidelen*1000)
#             ofextn = 'ntuples_%s'%(fsidelen_index)
#             jdlName = 'submitJobs_%s_ntuples_%s_%s.jdl'%(sample,fsidelen_index,findex)
#             jdlFile = open('%s/%s'%(jdlDir,jdlName),'w')
#             jdlFile.write(common_command)
#             run_command =  'Arguments  = %s %s %s %s %s %s %s $(process) \nQueue 1\n\n' %(fname,findex,nevents,sidelength,ofextn,sample,iloop)
#             jdlFile.write(run_command)        
#             jdlFile.close()
#             subFile.write("condor_submit %s\n"%jdlName)
#             findex = findex + 1

for sample in samplelist_PU200_Emyr:
    for sidelength in triangle_side_list:
        dirpath = '%s/stage2_emulator_tests/Emyr/%s/'%(os.environ["HOME"],sample)
        findarg = '%s -name \"*.root\"'%(dirpath)
        p = subprocess.Popen([("find %s "%findarg)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)
        filelist, errors = p.communicate()
        filelist = filelist.split('\n')
        filelist.remove("")
        nFiles = len(filelist)
        findex  = 0
        fsidelen = float(sidelength)
        fsidelen_index = int(fsidelen*1000)
        jdlName = 'submitJobs_%s_ntuples_%s.jdl'%(sample,fsidelen_index)
        jdlFile = open('%s/%s'%(jdlDir,jdlName),'w')
        jdlFile.write(common_command)
        for fname in filelist:
            #if findex < 200:
            print ("fname: %s, sidelength: %s, index: %s"%(fname,sidelength,findex))
            ofextn = 'ntuples_%s'%(fsidelen_index)
            run_command =  'Arguments  = %s %s %s %s %s %s %s $(process) \nQueue 1\n\n' %(fname,findex,nevents,sidelength,ofextn,sample,iloop)
            jdlFile.write(run_command)        
            findex = findex + 1
        jdlFile.close()
        subFile.write("condor_submit %s\n"%jdlName)
        
subFile.close()

print ("Now run from directory : %s/"%jdlDir)

