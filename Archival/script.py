import pandas as pd
import math

def RAsex_to_RAdec(fRAsex):
    frah = float(fRAsex[0:2])
    fram = float(fRAsex[2:4])
    fras = float(fRAsex[4:])

    return ((1/1 * frah) + (1/60 *fram) + (1/3600 *fras))* (360/24)
    
      
def DEsex_to_DEdec(fDEsex):
    fded = float(fDEsex[0:3])
    fdem = float(fDEsex[3:5])
    fdes = float(fDEsex[5:])    
    fDEdec = (math.fabs(fded)*3600.0+fdem*60.0+fdes)/3600.0
    if fDEsex[0] == '-':
        fDEdec = fDEdec * -1
    return fDEdec
    
df = pd.read_csv('/fred/oz100/sgoode/SMOTEs/Candidate_Dataframes/2015_01_4hr.csv')
f = open('/fred/oz100/sgoode/SMOTEs/Candidate_Dataframes/2015_01_4hr_linetest.dat', 'w+')
for i,r in df.iterrows():
	f.write('{},{},{}\n'.format(RAsex_to_RAdec(r.Filename[3:13]), DEsex_to_DEdec(r.Filename[13:24]), r["Detection MJD"]))
f.close()
