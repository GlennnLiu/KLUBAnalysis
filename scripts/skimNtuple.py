#!/usr/bin/env python

import os,sys
import optparse
import fileinput
import commands
import time
import glob
import subprocess
from os.path import basename
import ROOT


def isGoodFile (fileName) :
    ff = ROOT.TFile (fname)
    if ff.IsZombie() : return False
    if ff.TestBit(ROOT.TFile.kRecovered) : return False
    return True
    

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

def parseInputFileList (fileName) :
    filelist = []
    with open (fileName) as fIn:
        for line in fIn:
            line = (line.split("#")[0]).strip()
            if line:
                filelist.append(line)
    return filelist

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


if __name__ == "__main__":

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option ('-i', '--input'     , dest='input'     , help='input folder'                          , default='none')
    parser.add_option ('-x', '--xs'        , dest='xs'        , help='sample xs'                             , default='1.')
    parser.add_option ('-f', '--force'     , dest='force'     , help='replace existing reduced ntuples'      , default=False)
    parser.add_option ('-o', '--output'    , dest='output'    , help='output folder'                         , default='none')
    # parser.add_option ('-q', '--queue'     , dest='queue'     , help='batch queue'                           , default='cms')
    parser.add_option ('-q', '--queue'     , dest='queue'     , help='batch queue'                           , default='short')
    parser.add_option ('-r', '--resub'     , dest='resub'     , help='resubmit failed jobs'                  , default='none')
    parser.add_option ('-v', '--verb'      , dest='verb'      , help='verbose'                               , default=False)
    parser.add_option ('-s', '--sleep'     , dest='sleep'     , help='sleep in submission'                   , default=False)
    parser.add_option ('-d', '--isdata'    , dest='isdata'    , help='data flag'                             , default=False)
    parser.add_option ('-T', '--tag'       , dest='tag'       , help='folder tag name'                       , default='')
    parser.add_option ('-H', '--hadd'      , dest='hadd'      , help='hadd the resulting ntuples'            , default='none')
    parser.add_option ('-c', '--config'    , dest='config'    , help='skim config file'                      , default='none')
    parser.add_option ('-n', '--njobs'     , dest='njobs'     , help='number of skim jobs'                   , default=100, type = int)
    parser.add_option ('-k', '--kinfit'    , dest='dokinfit'  , help='run HH kin fitter'                     , default=True)
    parser.add_option ('-m', '--mt2'       , dest='domt2'     , help='run stransverse mass calculation'      , default=True)
    parser.add_option ('-y', '--xsscale'   , dest='xsscale'   , help='scale to apply on XS for stitching'    , default='1.0')
    parser.add_option ('-z', '--htcut'     , dest='htcut'     , help='HT cut for stitching on inclusive'     , default='-999.0')
    parser.add_option ('-e', '--njets'     , dest='njets'     , help='njets required for stitching on inclusive'     , default='-999')
    parser.add_option ('-t', '--toprew'    , dest='toprew'    , help='is TT bar sample to compute reweight?' , default=False)
    parser.add_option ('-b', '--topstitch' , dest='topstitch' , help='type of TT gen level decay pruning for stitch'        , default='0')
    parser.add_option ('-g', '--genjets'   , dest='genjets'   , help='loop on genjets to determine the number of b hadrons' , default=False)
    parser.add_option ('-w', '--weight'    , dest='weightHH'  , help='histo map for hh reweight'             , default='0')
    parser.add_option ('-a', '--ishhsignal', dest='ishhsignal', help='isHHsignal'                            , default=False)
    parser.add_option ('--kl',               dest='klreweight', help='kl for dynamic reweight'              , default='-999.0')
    parser.add_option ('--kt',               dest='ktreweight', help='kt for dynamic reweight'              , default='-999.0')
    parser.add_option ('--c2',               dest='c2reweight', help='c2 for dynamic reweight'              , default='-999.0')
    parser.add_option ('--cg',               dest='cgreweight', help='cg for dynamic reweight'              , default='-999.0')
    parser.add_option ('--c2g',              dest='c2greweight', help='c2g for dynamic reweight'            , default='-999.0')
    parser.add_option ('--susy',             dest='susyModel' , help='name of susy model to select'         , default='NOTSUSY')
    
    (opt, args) = parser.parse_args()

    currFolder = os.getcwd ()

    # verify the result of the process
    # ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

    if (opt.hadd != 'none') :

        scriptFile = open (opt.output + '/hadder.sh', 'w')
        scriptFile.write ('#!/bin/bash\n')
        scriptFile.write ('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
        scriptFile.write ('cd /data_CMS/cms/govoni/CMSSW_7_4_5/src\n')
        scriptFile.write ('export SCRAM_ARCH=slc6_amd64_gcc472\n')
        scriptFile.write ('eval `scram r -sh`\n')
        scriptFile.write ('cd %s\n'%currFolder)
        scriptFile.write ('source scripts/setup.sh\n')
        scriptFile.write ('mkdir ' + opt.output + '/singleFiles\n')
        scriptFile.write ('mv ' + opt.output + '/* ' + opt.output + '/singleFiles\n')
        scriptFile.write ('hadd ' + opt.output + '/total.root ' + opt.output + '/singleFiles/*.root\n')
        scriptFile.write ('touch ' + opt.output + '/done\n')
        scriptFile.write ('echo "Hadding finished" \n')
        scriptFile.close ()
        os.system ('chmod u+rwx ' + opt.output + '/hadder.sh')
        command = ('/opt/exp_soft/cms/t3/t3submit -q cms \'' +  opt.output + '/hadder.sh\'')
        os.system (command)
        sys.exit (0)

    # verify the result of the process
    # ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

    if (opt.resub != 'none') :
        if (opt.input == 'none') :
            print 'input folder to be checked missing\n'
            print '(this is the folder that contains the jobs to be submitted)'
            sys.exit (1)

        if opt.input[-1] == '/' : opt.input = opt.input[:-1]
        tagname = opt.tag + "/" if opt.tag else ''
        opt.input = tagname + 'SKIM_' + basename (opt.input)
        jobs = [word.replace ('_', '.').split ('.')[1] for word in os.listdir (opt.input) if 'skim' in word]
        missing = []
        
        # check the existence of the done file
        for num in jobs :
            if not os.path.exists (opt.input + '/done_' + num) :
                if opt.verb : print num, ' : missing done file'
                missing.append (num)

        # check the log file
        for num in jobs :
            # get the log file name
            filename = opt.input + '/skimJob_' + num + '.sh'
#            print os.path.exists (filename) 
            with open (filename, 'r') as myfile :
                data = [word for word in myfile.readlines () if 'log' in word]
            rootfile = data[0].split ()[2]
            if not os.path.exists (rootfile) :
                if opt.verb : print num, 'missing root file', rootfile
                missing.append (num)
                continue
            if not isGoodFile (rootfile) :
                if opt.verb : print num, 'root file corrupted', rootfile
                missing.append (num)
                continue
            logfile = data[0].split ()[-1]
            if not os.path.exists (logfile) :
                if opt.verb : print num, 'missing log file'
                missing.append (num)
                continue
            with open (logfile, 'r') as logfile :
                problems = [word for word in logfile.readlines () if 'Error' in word and 'TCling' not in word]
                if len (problems) != 0 :
                    if opt.verb : print num, 'found error ', problems[0]
                    missing.append (num)
        print 'the following jobs did not end successfully:'
        print missing   
        for num in missing :
            command = '`cat ' + opt.input + '/submit.sh | grep skimJob_' + num + '.sh | tr "\'" " "`'
            if opt.verb : print command
        if (opt.resub == 'run') :
            for num in missing :
                command = '`cat ' + opt.input + '/submit.sh | grep skimJob_' + num + '.sh | tr "\'" " "`'
                time.sleep (int (num) % 5)
                os.system (command)
        sys.exit (0)

    # submit the jobs
    # ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

#    skimmer = './bin/skimNtuple.exe'
    # skimmer = 'skimNtupleInclusive_Luca.exe'
    skimmer = 'skimNtuple.exe'
    # skimmer = 'getSelectionEfficiencyNew.exe'

    if opt.config == 'none' :
        print 'config file missing, exiting'
        sys.exit (1)

    if opt.input[-1] == '/' : opt.input = opt.input[:-1]
    if opt.output == 'none' : opt.output = opt.input + '_SKIM'
   
    if not os.path.exists (opt.input) :
        print 'input folder', opt.input, 'not existing, exiting'
        sys.exit (1)
    if not opt.force and os.path.exists (opt.output) :
        print 'output folder', opt.output, 'existing, exiting'
        sys.exit (1)
    elif os.path.exists (opt.output) :
        os.system ('rm -rf ' + opt.output + '/*')
    os.system ('mkdir ' + opt.output)
    os.system ('cp ' + opt.config + " " + opt.output)
    
    #inputfiles = glob.glob (opt.input + '/*.root')    
    inputfiles = parseInputFileList (opt.input)
    if opt.njobs > len (inputfiles) : opt.njobs = len (inputfiles)
    nfiles = (len (inputfiles) + len (inputfiles) % opt.njobs) / opt.njobs
    inputlists = [inputfiles[x:x+nfiles] for x in xrange (0, len (inputfiles), nfiles)]

    tagname = "/" + opt.tag if opt.tag else ''
    jobsDir = currFolder + tagname + '/SKIM_' + basename (opt.input)
    jobsDir = jobsDir.rstrip (".txt")
    if os.path.exists (jobsDir) : os.system ('rm -f ' + jobsDir + '/*')
    else                        : os.system ('mkdir ' + jobsDir)

    # proc = subprocess.Popen ('voms-proxy-info', stdout=subprocess.PIPE)
    # tmp = [word for word in proc.stdout.read ().split ('\n') if 'timeleft' in word]
    # if len (tmp) == 0 or int (tmp[0].split (':')[1]) < 10 : # hours
    #     os.system ('source /opt/exp_soft/cms/t3/t3setup')

    n = int (0)
    commandFile = open (jobsDir + '/submit.sh', 'w')
    for listname in inputlists : 
        #create a wrapper for standalone cmssw job
        listFileName = "filelist_%i.txt" % n
        thisinputlistFile = open(jobsDir + "/" + listFileName, 'w')
        for line in listname:
            thisinputlistFile.write(line+"\n")
        thisinputlistFile.close()
        scriptFile = open ('%s/skimJob_%d.sh'% (jobsDir,n), 'w')
        scriptFile.write ('#!/bin/bash\n')
        scriptFile.write ('export X509_USER_PROXY=~/.t3/proxy.cert\n')
        scriptFile.write ('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
        scriptFile.write ('cd /data_CMS/cms/govoni/CMSSW_7_4_5/src\n')
        scriptFile.write ('export SCRAM_ARCH=slc6_amd64_gcc472\n')
        scriptFile.write ('eval `scram r -sh`\n')
        scriptFile.write ('cd %s\n'%currFolder)
        scriptFile.write ('source scripts/setup.sh\n')
        command = skimmer + ' ' + jobsDir+"/"+listFileName + ' ' + opt.output + '/' + "output_"+str(n)+".root" + ' ' + opt.xs
        if opt.isdata :  command += ' 1 '
        else          :  command += ' 0 '    
        command += ' ' + opt.config + ' '
        if opt.dokinfit=="True" : command += " 1 "
        else                    : command += " 0 "
        command += " " + opt.xsscale
        command += " " + opt.htcut
        if opt.toprew=="True" : command += " 1 "
        else                  : command += " 0 "   
        if opt.genjets=="True": command += " 1 "
        else                  : command += " 0 "   
        command += (" " + opt.weightHH)
        command += " " + opt.topstitch
        if opt.domt2          : command += " 1 " ## inspiegabilmente questo e' un bool
        else                  : command += " 0 "
        if opt.ishhsignal     : command += " 1 "
        else                  : command += " 0 "
        command += (" " + opt.njets)
        command += (" " + opt.klreweight + " " + opt.ktreweight + " " + opt.c2reweight + " " + opt.cgreweight + " " + opt.c2greweight)
        command += (" " + opt.susyModel)

        command += ' >& ' + opt.output + '/' + "output_" + str(n) + '.log\n'
        scriptFile.write (command)
        scriptFile.write ('touch ' + jobsDir + '/done_%d\n'%n)
        scriptFile.write ('echo "All done for job %d" \n'%n)
        scriptFile.close ()
        os.system ('chmod u+rwx %s/skimJob_%d.sh'% (jobsDir,n))

        # command = ('/opt/exp_soft/cms/t3/t3submit -q ' + opt.queue + ' \'' + jobsDir + '/skimJob_' + str (n) + '.sh\'')
        command = '/opt/exp_soft/cms/t3/t3submit -'+opt.queue + ' ' + jobsDir + '/skimJob_' + str (n) + '.sh'
        if opt.sleep : time.sleep (0.1)
        os.system (command)
        commandFile.write (command + '\n')
        n = n + 1
    commandFile.close ()
    


