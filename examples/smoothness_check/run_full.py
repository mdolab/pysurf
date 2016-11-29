import subprocess
import numpy as np

def full_run(wingTranslation):

    # print "Extracting curves.\n"
    # subprocess.call(["python script01_curveExtraction.py"], shell=True)
    #
    # print "Fixing trailing edge.\n"
    # subprocess.call(["python script02_fix_trailing_edge.py"], shell=True)

    print "Replacing values in scripts with desired translation.\n"
    python_string_replacement = "sed -i -e 's/wingTranslation =.*/wingTranslation = [[{}, {}, {}]]/g' script03_example_crm_fullInt.py".format(wingTranslation[0], wingTranslation[1], wingTranslation[2])

    cgns_string_replacement = "sed -i -e 's/cgns_utils translate wing_L1_temp.cgns .*/cgns_utils translate wing_L1_temp.cgns {} {} {}/g' script06_mergeMeshes.sh".format(wingTranslation[0], wingTranslation[1], wingTranslation[2])

    subprocess.call([python_string_replacement], shell=True)
    subprocess.call([cgns_string_replacement], shell=True)

    print "Computing intersections and marching collar mesh.\n"
    subprocess.call(["python script03_example_crm_fullInt.py"], shell=True)

    print "Merging collar meshes.\n"
    subprocess.call(["python script04_mergeCollar.py"], shell=True)

    print "Extruding collar mesh using pyHyp.\n"
    subprocess.call(["python script05_runPyhyp.py"], shell=True)

    print "Merging meshes into single .cgns file.\n"
    subprocess.call(["sh script06_mergeMeshes.sh"], shell=True)

    print "Running ADflow using the combined meshes.\n"
    subprocess.call(["cd ADflow; sh run_check.sh"], shell=True)

if __name__ == "__main__":

    nStates = 2

    wingTranslation = np.zeros((nStates,3))
    wingTranslation[:,2] = np.linspace(-.1, 3., nStates)

    for wingT in wingTranslation:
        full_run(wingT)
