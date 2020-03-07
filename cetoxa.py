#!/usr/bin/env python3

#=====================================================================
# Script to dock one ligand to CETOXA using QVina-2
# Authors: Natasha Kamerlin    <natasha.kamerlin@icm.uu.se>
#          David van der Spoel <david.vanderspoel@icm.uu.se>
#=====================================================================

import os, glob, getpass, math, sys, argparse
from random import randint
import pandas as pd
import tempfile, shutil

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

def dock(pathway, tmpdir, target_dir, targets, binding_dir,
         ligand, ncpu, verbose):
    for ppp in range(len(targets)):
        p=targets[ppp]
        target_pdbqt    = target_dir + "/" + p
        pdbcode = p[:-6]
        bindingsites_txt = binding_dir + "/" + pdbcode + "_bindingsites.txt"

        #Get the box dimensions
        box = default_box_size(ligand)
        xcm,ycm,zcm,bindingsites=get_bindingsite(bindingsites_txt)
        for i in range(int(bindingsites)):

            logfn            = tmpdir + "/" + pdbcode + "_site" + str(i) + ".log"
            outfn            = tmpdir + "/" + pdbcode + "_site" + str(i) + ".models.pdbqt"
            cmd = ( "%s --log %s --out %s --exhaustiveness 8 --receptor %s --ligand %s --center_x %f  --center_y %f  --center_z %f --size_x %f --size_y %f --size_z %f --cpu %d" % ( pathway, logfn, outfn, target_pdbqt, ligand, float(xcm[i]), float(ycm[i]), float(zcm[i]), box.dimx, box.dimy, box.dimz, ncpu ) )
            print("%s site %d" % ( pdbcode, i ) )
            if verbose:
                os.system(cmd)
            else:
                os.system("%s >& /dev/null" % cmd)

class DockingResult:
    def __init__(self, family, name, score, site):
        self.family = family.strip()
        self.name   = name.strip()
        self.score  = score
        self.site   = site
    def __lt__(self, other):
        if self.family == other.family:
            return self.score < other.score
        else:
            return self.family < other.family
    def __gt__(self, other):
        if self.family == other.family:
            return self.score > other.score
        else:
            return self.family > other.family
    def __eq__(self, other):
        return self == other
    def __le__(self, other):
        return __eq__(other) or __lt__(other)
    def __ge__(self, other):
        return __eq__(other) or __gt__(other)
    def __ne__(self, other):
        return not __eq__(other)
    def printCsv(self, outputFile):
        outputFile.write("%s,%s,%g,%s\n" % ( self.family, self.name, self.score, self.site ) )
    
def extract(protInfo, tmpdir, targets, binding_dir, ligand, outfile):
    output = []
    for ppp in range(len(targets)):
        p      = targets[ppp]
        target = p[:-6]
        protname=protInfo['Name'][protInfo['PDBID']==target].to_string(index=False)
        protfam=protInfo['Family'][protInfo['PDBID']==target].to_string(index=False)
        bindingsites_txt = binding_dir + "/" + p[:-6] + "_bindingsites.txt"
        bindingsites     = len(open(bindingsites_txt).readlines())

        energy = 1000.0
        index  = None
        for sites in range(bindingsites):
            model = tmpdir + "/" + target + "_site" + str(sites) + ".models.pdbqt"
            if os.path.exists(model):
                with open(model, "r") as protsite:
                    for line in protsite.readlines():
                        if line.strip().find("REMARK VINA RESULT") >= 0:
                            newenergy = float(line.split()[3])
                            if newenergy < energy:
                                energy = newenergy
                                index  = sites
                            break
            else:
                print("No such output file %s" % model)
        if index != None:
            output.append(DockingResult(protfam, protname, energy, index))
        else:
            # Print a warning, but keep going
            print("Warning, no energy for target %s" % target)

    #Sort the energies
    sorted_output = sorted(output)

    with open(outfile, "w") as outputfile:
        # Print the results
        outputfile.write("%s,%s,%s,%s\n" % ( 'Class', 'Target', 'Score', 'Site' ) )
        for data in sorted_output:
            data.printCsv(outputfile)

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
    for qv in [ "qvina02" ]:
        qvina = which(qv)
        if qvina:
            break
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--infile", help="Compound file for docking in pdbqt format",   type=str, default=None)
    parser.add_argument("-o", "--outfile", help="CSV file for writing table of results",   type=str, default="results.csv")
    parser.add_argument("-vina", "--vina", help="Pathway to qvina02 executable", default=qvina)
    parser.add_argument("-v", "--verbose", help="Write to console and leave temporary data on disk", action="store_true")
    parser.add_argument("-rm", "--remove_temp", help="Remove temporary data", action="store_true")
    parser.add_argument("-analyze", "--analyze", help="Do not run the docking, just analyze existing results", action="store_true")
    parser.add_argument("-ncpu", "--ncpu", help="Number of cores to use", type=int, default=1)
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
        print("Please add qvina02 to your search path or use the -vina flag")
    target_dir  = cetoxa_dir + "/data/Targets/pdbqt"
    targets     = get_pdbqt(target_dir)
    targetList  = cetoxa_dir + "/data/Targets/Target-list.csv"
    ligand      = args.infile
    binding_dir = cetoxa_dir + "/data/Targets/bindingsites"

    protInfo    = pd.read_csv(targetList,sep=',')

    # Dock
    print("A Computational Ecotoxicity Assay")
    print("https://doi.org/10.26434/chemrxiv.11944371.v1")
    if args.verbose:
        print("Calling %s to perform docking" % args.vina)
    tmpdir = "temp"
    os.makedirs(tmpdir, exist_ok=True)
    if not args.analyze:
        dock(args.vina, tmpdir, target_dir, targets, binding_dir,
             ligand, args.ncpu, args.verbose)
    if args.verbose:
        print("Docking complete.")
    
    # Extract top scores
    if args.verbose:
        print("Extracting scores from %s to %s." % (tmpdir, args.outfile))
    extract(protInfo, tmpdir, targets, binding_dir, ligand, args.outfile)
    
    print("A summary of results is in %s." % args.outfile)
    if args.remove_temp:
        print("The intermediate results in %s will be removed." % tmpdir)
        shutil.rmtree(tmpdir)
    
