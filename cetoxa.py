#!/usr/bin/env python3

#=====================================================================
# Script to dock one ligand to CETOXA using QVina-2
# Authors: Natasha Kamerlin    <natasha.kamerlin@icm.uu.se>
#          David van der Spoel <david.vanderspoel@icm.uu.se>
#=====================================================================

import os, glob, getpass, math, sys, argparse
from random import randint
import pandas as pd

class BoxSize():
    def __init__(self, dimx, dimy, dimz):
        self.dimx = dimx
        self.dimy = dimy
        self.dimz = dimz

def default_box_size(filename):

    coord_x = []
    coord_y = []
    coord_z = []

    with open(filename, "r") as Structure:
        for line in Structure:
            try:
                if line[0:6] == "HETATM" or line[0:4] == "ATOM":
                    coord_x.append(float(line[30:38]))
                    coord_y.append(float(line[38:46]))
                    coord_z.append(float(line[46:54]))
            except:
                print ("Unknown error detected.")

        N = len(coord_x)

        xmin = min(coord_x)-5.0
        ymin = min(coord_y)-5.0
        zmin = min(coord_z)-5.0

        xmax = max(coord_x)+5.0
        ymax = max(coord_y)+5.0
        zmax = max(coord_z)+5.0

        random = []
        for i in range(3):
            random.append(randint(0,1))

        if random[0] > 0:
            xmax += 5.0
        else:
            xmin -= 5.0

        if random[1] > 0:
            ymax += 5.0
        else:
            ymin -= 5.0

        if random[2] > 0:
            zmax += 5.0
        else:
            zmin -= 5.0

        dimx = xmax - xmin
        dimy = ymax - ymin
        dimz = zmax - zmin

        if dimx < 22.5:
            offsetx= (22.5 - dimx)/2
            xmax += offsetx
            xmin -= offsetx
            dimx = xmax - xmin
        if dimy < 22.5:
            offsety = (22.5 - dimy)/2
            ymax += offsety
            ymin -= offsety
            dimy = ymax - ymin
        if dimz < 22.5:
            offsetz = (22.5 - dimz)/2
            zmax += offsetz
            zmin -= offsetz
            dimz = zmax - zmin

        return BoxSize(dimx, dimy, dimz)

def get_pdbqt(dir):
    wd = os.getcwd()
    os.chdir(dir)
    pdbqt = glob.glob("*.pdbqt")
    os.chdir(wd)
    return pdbqt

def get_bindingsite(filename):
    bindingsites=0
    xcm=[]
    ycm=[]
    zcm=[]

    with open(filename, "r") as Sites:
        for line in Sites:
            bindingsites+=1
            xcm.append(line.split()[0])
            ycm.append(line.split()[1])
            zcm.append(line.split()[2])

    return xcm,ycm,zcm,bindingsites

def dock(pathway,target_dir,targets,binding_dir,ligand):
    for ppp in range(len(targets)):
        p=targets[ppp]
        target_pdbqt    = target_dir + "/" + p
        bindingsites_txt = binding_dir + "/" + p[:-6] + "_bindingsites.txt"

        #Get the box dimensions
        box = default_box_size(ligand)
        xcm,ycm,zcm,bindingsites=get_bindingsite(bindingsites_txt)

        for i in range(int(bindingsites)):

            logfn            = "results/" + p[:-6] + "_site" + str(i) + ".log"
            outfn            = "results/" + p[:-6] + "_site" + str(i) + ".models.pdbqt"
            os.system("  %s --log %s --out %s --exhaustiveness 8 --cpu 1 --receptor %s --ligand %s --center_x %f  --center_y %f  --center_z %f --size_x %f --size_y %f --size_z %f\n" % ( pathway,logfn, outfn, target_pdbqt, ligand, float(xcm[i]), float(ycm[i]), float(zcm[i]), box.dimx, box.dimy, box.dimz ) )


def extract(protInfo,targets,binding_dir,ligand):
    output=[]
    for ppp in range(len(targets)):
        p = targets[ppp]
        target = p[:-6]
        protname=protInfo['NAME'][protInfo['PDBID']==target].to_string(index=False)
        protfam=protInfo['FAMILY'][protInfo['PDBID']==target].to_string(index=False)
        bindingsites_txt = binding_dir + "/" + p[:-6] + "_bindingsites.txt"
        bindingsites = len(open(bindingsites_txt).readlines())

        energy=100.0
        index=-1
        arr=[]
        for sites in range(bindingsites):
            model = "results/" + target + "_site" + str(sites) + ".models.pdbqt"
            protsite = open(model, "r")
            for line in protsite:
                if line.find("REMARK VINA RESULT") >= 0:
                    if (float(line.split()[3]) < energy):
                        energy = float(line.split()[3])
                        index=sites
                    break
            protsite.close()
        if (energy>99):
            print("Warning, error in energy")
            exit()
        output.append([protfam,protname,energy])

    #Sort the energies
    sorted_output=sorted(output,key=lambda x:x[2]) 

    header=['Protein Family','Protein Name','Vina Score']
    row_format="{:<25}{:<55}{:<10}\n"

    outputfile=open(ligand[:-6]+'.dockingresults.txt','w')

    #Print the results
    outputfile.write(row_format.format( *header))

    for data in sorted_output:
         outputfile.write(row_format.format(*data))

    outputfile.close()

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--infile", help="Compound file for docking",   type=str, default=None)
    parser.add_argument("-o", "--outfile", help="Logfile for writing",   type=str,
        default="results.txt")
    for qv in [ "qvina02" ]:
        qvina = which(qv)
        if qvina:
            break
    parser.add_argument("-vina", "--vina", help="Pathway to qvina02 executable", default=qvina)

    args = parser.parse_args()
    return args

# Main processing   
if __name__ == '__main__':
    cetoxa_dir  = os.path.dirname(os.path.realpath(__file__))
    args        = parseArguments()
    if not args.infile:
        print("Please give me an input file or run with -h")
        exit(0)
    if not args.vina:
        print("Please add qvina-w or qvina-2 to your search path or use the -vina flag")
    exit(0)
    target_dir  = cetoxa_dir + "/data/Targets/pdbqt"
    targets     = get_pdbqt(target_dir)
    targetList  = cetoxa_dir + "/data/Targets/Target-list.txt"
    ligand      = args.infile
    binding_dir = cetoxa_dir + "/data/Targets/bindingsites"

    protInfo    = pd.read_csv(targetList,sep='|')

    # Dock
    print("Calling QVina-2 to perform docking")
    dock(args.vina, target_dir, targets, binding_dir, ligand)
    print("Docking complete!")

    # Extract top scores
    print("Extracting scores")
    extract(protInfo, targets, binding_dir, ligand)

    print("Completed!")
