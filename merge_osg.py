#!/usr/bin/env python3
"""
to merge the mc files in /volatile/clas12/osg2
"""

import numpy as np
import sys, os, subprocess


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-j","--jobNum", help="job number", type = str)
    parser.add_argument("-c","--collectionNum", help="select 0 to 19 for 10k jobs", type = int)
    parser.add_argument("-p","--process", help="process used to generate events", type = str)
    parser.add_argument("-pol","--polarity", help="polarity used for the detector simulation", type = str)
    parser.add_argument("-u","--username", help="user who simulated the file", type = str, default = 'sangbaek')

    args = parser.parse_args()

    jobNum = args.jobNum
    collectionNum = args.collectionNum
    process = args.process
    polarity = args.polarity
    user = args.username

    jobStart = collectionNum * 500
    jobEnd   = collectionNum * 500 + 500
    fileInds = np.linspace(jobStart, jobEnd-1, 500, dtype = int)

    osgDir   = "/volatile/clas12/osg2"
    jobDir   = "{}/{}/job_{}/output/simu_".format(osgDir, username, jobNum)
    print("merge {} {} to {} ".format(jobDir, jobStart, jobEnd))

    outputDir = "/volatile/clas12/sangbaek/regular_backup/temp/"
    targetDir = "{}{}".format(outputDir, jobNum)
    mergedFile = "{}/{}/{}/{}_{}.hipo".format(outputDir, process, polarity, jobNum, collectionNum)
    runCmd = ["hipo-utils", "-filter", "-b", "'RUN::*,RAW::epics,RAW::scaler,HEL::flip,HEL::online,REC::*,RECFT::*,MC::*'", "-merge", "-o", mergedFile]

    print("creating sym links inside {}".format(targetDir))

    subprocess.run(['mkdir','-p',targetDir])
    for fileInd in fileInds:
        thisFile   = "{}{}/dst.hipo".format(jobDir, fileInd)
        symlinkFile = "{}/{}_{}.hipo".format(targetDir, jobNum, fileInd)
        try:
            subprocess.check_output(['ls', thisFile])
        except:
            continue
        run = subprocess.run(['ln','-s', thisFile, symlinkFile])  # quieter
        runCmd.append("{}".format(symlinkFile))

    print("merging the files as {}".format(mergedFile))
    print("running {}".format(runCmd))
    

    subprocess.run(runCmd)