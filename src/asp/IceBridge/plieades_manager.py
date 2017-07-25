#!/usr/bin/env python
# __BEGIN_LICENSE__
#  Copyright (c) 2009-2013, United States Government as represented by the
#  Administrator of the National Aeronautics and Space Administration. All
#  rights reserved.
#
#  The NGT platform is licensed under the Apache License, Version 2.0 (the
#  "License"); you may not use this file except in compliance with the
#  License. You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
# __END_LICENSE__

# Top level program that 

import os, sys, optparse, datetime, time, subprocess, logging, multiprocessing
import re, shutil, time, getpass

import os.path as P

# The path to the ASP python files and tools
basepath      = os.path.dirname(os.path.realpath(__file__)) # won't change, unlike syspath
pythonpath    = os.path.abspath(basepath + '/../Python')     # for dev ASP
libexecpath   = os.path.abspath(basepath + '/../libexec')    # for packaged ASP
icebridgepath = os.path.abspath(basepath + '/../IceBridge')  # IceBridge tools

# Prepend to Python path
sys.path.insert(0, basepath)
sys.path.insert(0, pythonpath)
sys.path.insert(0, libexecpath)
sys.path.insert(0, icebridgepath)

import icebridge_common, pbs_functions

#asp_system_utils.verify_python_version_is_supported()

# Prepend to system PATH
os.environ["PATH"] = basepath       + os.pathsep + os.environ["PATH"]
os.environ["PATH"] = pythonpath     + os.pathsep + os.environ["PATH"]
os.environ["PATH"] = libexecpath    + os.pathsep + os.environ["PATH"]
os.environ["PATH"] = icebridgepath  + os.pathsep + os.environ["PATH"]




# Constants used in this file
PBS_QUEUE = 'devel'
#NODE_TYPE = 'san' # Sandy bridge = 16 core,  32 GB mem, SBU 1.82
NODE_TYPE = 'ivy' # Ivy bridge   = 20 core,  64 GB mem, SBU 2.52
#NODE_TYPE = 'has' # Haswell      = 24 core, 128 GB mem, SBU 3.34
#NODE_TYPE = '???' # Broadwell    = 28 core, 128 GB mem, SBU 4.04
NUM_ORTHO_NODES     = 1
NUM_ORTHO_PROCESSES = 10
NUM_ORTHO_THREADS   = 2
NUM_BATCH_NODES     = 1
NUM_BATCH_THREADS   = 8 # MGM is limited to 8 threads
NUM_BATCH_PROCESSES = 3
GROUP_ID = 1234 # TODO
MAX_PROCESS_HOURS = 40

# Wait this long between checking for job completion
SLEEP_TIME = 30

ALL_RUN_LIST       = '/u/smcmich1/icebridge/full_run_list.txt'
SKIP_RUN_LIST      = '/u/smcmich1/icebridge/run_skip_list.txt'
COMPLETED_RUN_LIST = '/u/smcmich1/icebridge/completed_run_list.txt'

LOG_FOLDER              = '/u/smcmich1/icebridge/manager_logs'
UNPACK_FOLDER           = '/u/smcmich1/icebridge/data'
CALIBRATION_FILE_FOLDER = '/u/smcmich1/icebridge/calib_files'
REFERENCE_DEM_FOLDER    = '/u/smcmich1/icebridge/reference_dems'
REMOTE_INPUT_FOLDER     = 'lou:/u/oalexan1/projects/data/icebridge'
REMOTE_OUTPUT_FOLDER    = 'lou:/u/oalexan1/projects/data/icebridge'


BUNDLE_LENGTH = 10

# TODO: Move this!
def getOutputFolder(run):
    '''Get the output folder for a run'''
    return os.path.join(UNPACK_FOLDER, run[0] +'_'+ run[1])

def getUser():
    '''Return the current user name.'''
    return getpass.getuser()

def submitJob(jobName, scriptPath, args, logPrefix):
    '''Trivial wrapper to fill in the parts of the job command that never change'''
    return pbs_functions.submitJob(jobName, PBS_QUEUE, MAX_PROCESS_HOURS, GROUP_ID, NODE_TYPE, scriptPath, args, logPrefix)

def sendEmail(address, subject, body):
    '''Send a simple email from the command line'''
    os.system('mail -s '+subject+' '+address+' <<< '+body)

def readRunList(path):
    '''Reads a list of runs in this format: GR 20110411.
       Output is a list of (site, yyyymmdd) pairs.'''

    runList = []
    with open(path, 'r') as f:
        for line in f:
            if line[0] == '#': # Skip comments
                continue
            parts = line.split() # Site, yyyymmdd
            runList.append(parts)
    return runList

def addToRunList(path, run):
    '''Append a run to a list file'''
    with open(path, 'a') as f:
        f.write(path[0]+' '+path[1]+'\n')
    

def getRunsToProcess(allRuns, skipRuns, doneRuns):
    '''Go through the run lists to get the list of runs that
       should be processed.'''
       
    runList = []
    for run in allRuns:
        if (run in skipRuns) or (run in doneRuns):
            continue
        runList.append(run)
    return runList

def isRunReadyForProcessing(run):
    '''Return true if the run is unpacked locally'''
    
    # Just check if all the subfolders are there.
    # - If anything is missing it should be handled in a future processing step
    runFolder = getOutputFolder(run)
    subFolders = ['jpeg', 'lidar', 'ortho']
    for f in subFolders:
        path = os.path.join(runFolder, f)
        if not os.path.exists(path):
            return False
    return True

def retrieveRunData(run):
    '''Retrieve the data for the specified run from Lou.'''

    logger = logging.getLogger(__name__)
    
    # Skip retrieval if we already have the data
    if isRunReadyForProcessing(run):
        logger.info('No need to retrieve run ' + str(run) + ', it is already on disk.')
        return

    logger.info('Retrieving data for run ' + str(run))

    fileName = run[0]+'_'+run[1]+'.tar'
    louPath  = os.path.join(REMOTE_INPUT_FOLDER, fileName)

    cmd = 'shiftc --wait --sync --verify --extract-tar ' + louPath + ' ' + UNPACK_FOLDER
    logger.info(cmd)
    status = os.system(cmd)
    if status != 0:
        raise Exception('Failed to copy data for run ' + str(run))


def getFrameRange(run):
    '''Return the min and max frame which exists on disk for a run'''

    # Loop through all the files in the jpeg folder to get the frame range
    outputFolder = getOutputFolder(run)
    jpegFolder   = os.path.join(outputFolder, 'jpeg')
    if not os.path.exists(jpegFolder):
        raise Exception('Cannot get frame range from run, data not available: ' + str(run))

    minFrame = 9999999
    maxFrame = 0
    jpegFiles = os.listdir(jpegFolder)
    for jpegFile in jpegFiles:
        
        inputPath = os.path.join(jpegFolder, jpegFile)
        
        # Skip non-image files
        ext = os.path.splitext(jpegFile)[1]
        if ext != '.JPG':
            continue
        
        # Update frame range
        frame = icebridge_common.getFrameNumberFromFilename2(inputPath)
        if frame < minFrame:
            minFrame = frame
        if frame > maxFrame:
            maxFrame = frame

    return (minFrame, maxFrame)


# Do we need this function?
#def getConversionJobName(run, startFrame):
#    '''Get the PBS name to use for a conversion job'''
#    # PBS job names are limited to 15 characters
#    name = str(startFrame) + run[0] + run[1][2:]
#    if len(name) > pbs_functions.MAX_PBS_NAME_LENGTH:
#        raise Exception('Computed name ' +name+' for run '+str(run)+' and frame '+startFrame+' but it is too long')

def conversionIsFinished(run):
    '''Return true if this run is present and conversion has finished running on it'''
    outputFolder = getOutputFolder(run)
    return os.path.exists(os.path.join(outputFolder, 'conversion_complete'))
    # TODO: Option to do a full check for files here!

def runConversion(run):
    '''Run the conversion tasks for this run on the supercomputer nodes'''
    
    logger = logging.getLogger(__name__)
    
    # Get the frame range for the data.
    (minFrame, maxFrame) = getFrameRange(run)
    
    # Split the conversions across multiple nodes using frame ranges
    # - The first job submitted will be the only one that converts the lidar data
    numFrames        = maxFrame - minFrame + 1
    numFramesPerNode = numFrames / NUM_ORTHO_NODES
    numFramesPerNode += 1 # Make sure we don't cut off a few frames at the end
    
    # TODO: Is using the default output folder ok here?
    scriptPath = 'full_processing_script.py' # TODO
    args       = ('%s %s --site %s --yyyymmdd %d --skip-fetch --stop-after-convert --num-threads %d --num-processes %d' 
                  % (CALIBRATION_FILE_FOLDER, REFERENCE_DEM_FOLDER, run[0], run[1], NUM_ORTHO_THREADS, NUM_ORTHO_PROCESSES))
    
    baseName = run[0] + run[1][2:] # SITE + YYMMDD = 8 chars, leaves seven for frame digits.

    # Get the location to store the logs    
    outputFolder = getOutputFolder(run)
    pbsLogFolder = os.path.join(outputFolder, 'pbsLogs')
    os.mkdir(pbsLogFolder)
    
    currentFrame = 0
    for i in range(0,NUM_ORTHO_NODES):
        jobName  = str(currentFrame) + baseName
        thisArgs = args
        if i != 0: # Only the first job will convert lidar files
            thisArgs += ' --no-lidar-convert'
        logPrefix = os.path.join(pbsLogFolder + '_convert_' + jobName)
        logger.info('Submitting conversion job with args: '+thisArgs)
        pbs_functions.submitJob(jobName, scriptPath, thisArgs, logPrefix)

    # Wait for conversions to finish
    waitForRunCompletion(baseName)

    # Create a file to mark that conversion is finished for this folder
    os.system('touch '+ os.path.join(outputFolder, 'conversion_complete'))

def generateBatchList(run):
    '''Generate a list of all the processing batches required for a run'''

    # No actual processing is being done here so it can run on the PFE
    outputFolder = getOutputFolder(run)
    cmd = ('python full_processing_script.py --yyyymmdd %s --site %s --output-folder %s %s %s   --num-processes 1 --num-threads 1  --bundle-length %d  --skip-fetch  --skip-convert' % (run[1], run[0], outputFolder, CALIBRATION_FILE_FOLDER, REFERENCE_DEM_FOLDER, BUNDLE_LENGTH))
    os.system(cmd)
    
    listPath = os.path.join(outputFolder, 'processed/batch_commands_log.txt')
    return listPath
    

    

def submitBatchJobs(run, batchListPath):
    '''Read all the batch jobs required for a run and distribute them across job submissions'''

    if not os.path.exists(batchListPath):
        raise Exception('Failed to generate batch list file: ' + batchListPath)

    # Number of batches = number of lines in the file
    p = subprocess.Popen('wc -l '+batchListPath , stdout=subprocess.PIPE)
    textOutput, err = p.communicate()
    numBatches        = int(textOutput.split()[0])
    
    numBatchesPerNode = numBatches / NUM_BATCH_NODES
    numBatchesPerNode += 1 # Make sure we don't cut off a few batches at the end

    baseName = run[0] + run[1][2:] # SITE + YYMMDD = 8 chars, leaves seven for frame digits.

    # Call the tool which just executes commands from a file
    scriptPath = 'multi_process_command.py' # TODO

    outputFolder = getOutputFolder(run)
    pbsLogFolder = os.path.join(outputFolder, 'pbsLogs')

    currentBatch = 0
    for i in range(0,NUM_BATCH_NODES):
        jobName  = str(currentBatch) + baseName

        # Specify the range of lines in the file we want this node to execute
        args = ('%s, %d, %d, %d' % (batchListPath, NUM_BATCH_PROCESSES, currentBatch, currentBatch+numBatchesPerNode))

        logPrefix = os.path.join(pbsLogFolder + '_batch_' + jobName)
        logger.info('Submitting batch job with args: '+args)
        pbs_functions.submitJob(jobName, scriptPath, args, logPrefix)

    # Waiting on these jobs happens outside this function


def waitForRunCompletion(text):
    '''Sleep until all of the submitted jobs containing the provided text have completed'''

    user = getUser()
    stillWorking = True
    while stillWorking:
        os.sleep(SLEEP_TIME)
        stillWorking = False

        # Look through the list for jobs with the run's date in the name        
        jobList = pbs_functions.getActiveJobs(user)
        for job in jobList():
            if text in job:
                stillWorking = True
                break

def checkResults(batchListPath):
    '''Return (numOutputs, numProduced) to help validate our results'''
    
    batchOutputName = 'out-align-DEM.tif'
    
    numOutputs  = 0
    numProduced = 0
    with open(batchListPath, 'r') as f:
        for line in f:
            parts        = line.split()
            outputFolder = parts[4] # This needs to be kept up to date with the file format!
            targetPath = os.path.join(outputFolder, batchOutputName)
            numOutputs += 1
            if os.path.exists(targetPath):
                numProduced += 1
            
    return (numOutputs, numProduced)    


def packAndSendCompletedRun(runFolder, tarPath):
    '''Assembles and compresses the deliverable parts of the run'''
    
    # TODO: What do we want to deliver?
    # - The aligned DEM file from each batch folder
    # - Some low-size diagnostic information about the run
    # - Camera calibration files used?
    # - Information about how we created the run (process date, ASP version, etc)

    # Reassemble the files in a new folder or can we use links to do this?
    assemblyFolder = 'TODO'
    
    
    # Tar up the assembled files and send them at the same time using the shiftc command
    # - No need to use a compression algorithm here
    localFolder = 'TODO'
    
    louPath = os.path.join(REMOTE_OUTPUT_FOLDER, fileName)

    cmd = 'shiftc --wait --sync --extract-tar ' + assemblyFolder + ' ' + louPath
    status = os.system(cmd)
    if status != 0:
        raise Exception('Failed to pack/send results for run ' + str(run))


def main(argsIn):


    try:

        usage = '''usage: plieades_manager.py <options> '''
                      
        parser = optparse.OptionParser(usage=usage)

        ## Run selection
        #parser.add_option("--yyyymmdd",  dest="yyyymmdd", default=None,
        #                  help="Specify the year, month, and day in one YYYYMMDD string.")

        # TODO: What format for defining the run?
        parser.add_option("--manual-run",  dest="manualRun", default=None,
                          help="Manually specify a single run to process.") 
                          
        (options, args) = parser.parse_args(argsIn)

        #if len(args) != 2:
        #    print usage
        #    return -1


    except optparse.OptionError, msg:
        raise Usage(msg)

    logLevel = logging.INFO
    logger   = icebridge_common.setUpLogger(LOG_FOLDER, logLevel, 'plieades_manager_log')

    # Get the list of runs to process
    logger.info('Reading run lists...')
    allRuns  = readRunList(ALL_RUN_LIST      )
    skipRuns = readRunList(SKIP_RUN_LIST     )
    doneRuns = readRunList(COMPLETED_RUN_LIST)
    
    runList = getRunsToProcess(allRuns, skipRuns, doneRuns)

    runList = [('GR','20110504')] # DEBUG

   
    # Loop through the incomplete runs
    for run in runList:

        # TODO: Put this in a try/except block so it keeps going

        # Obtain the data for a run
        retrieveRunData(run)
        
        if not conversionIsFinished(run):
            # TODO: There needs to be a logging system
            logger.info('Converting data for run ' + str(run))
            runConversion(run)
         
        raise Exception('DEBUG')
         
        # Run command to generate the list of batch jobs for this run
        logger.info('Fetching batch list for run ' + str(run))
        batchListPath = generateBatchList(run)
        
        raise Exception('DEBUG')
        
        # Divide up batches into jobs and submit them to machines.
        logger.info('Submitting jobs for run ' + str(run))
        submitBatchJobs(run, batchListPath)
        
        # TODO: Don't wait for one run to finish before starting the next one!
        
        # Wait for all the jobs to finish
        logger.info('Waiting for job completion of run ' + str(run))
        waitForRunCompletion(run[0]) # TODO: Specialize this!
        logger.info('All jobs finished for run '+str(run))
        
        # Log the run as completed
        addToRunList(COMPLETED_RUN_LIST, run)

        # Generate a simple report of the results
        (numOutputs, numProduced) = checkResults(batchListPath)
        resultText = 'Created ' + str(numProduced) + ' out of ' + str(numOutputs) + ' output targets.'
        logger.info(resultText)

        # Notify that the run is done processing
        sendEmail('scott.t.mcmichael@nasa.gov', 'IB run complete: '+run[0]+'_'+run[1], resultText)
        
        raise Exception('DEBUG')
        
    # End loop through runs
    logger.info('==== plieades_manager script has finished! ====')


# Run main function if file used from shell
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


